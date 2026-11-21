# -*- coding: utf-8 -*-
"""[試作] 六極ファミリーの K2 に対する RDT・ADTS の解析式の「多項式モデル」と、それを目的関数にした六極の最適化 (R7)。

考え方: 六極の強さ k = (K2L_1, ..., K2L_F) (ファミリーごとの積分強さ) に対して、
  - 3 次の共鳴駆動項 h_jklm (1 次) とクロマティシティ・δ の項 (p ≥ 1) は k に **線形**:  h(k) = h_B + Σ_f k_f H_f
  - 4 次の項 h_jklm (写像の 4 次の生成関数、`rdt2nd.rdt4`) の六極の 2 次の部分 (Magnus の第 2 項) は k の **二次形式**、
    六極 × 偏向磁石の交差項は k に線形、八極・運動学・フリンジ・偏向磁石の 4 次は k に依らない定数
  - 3 次の項の非共鳴の摂動 (`rdt2nd.nf_correction_h4`、正規形の F4') は h3 の二次形式 = k の二次形式
  - ADTS (∂ν/∂J の 3 成分) は F4' の対角部分の実部 (4/π、2/π、4/π 倍) = k の二次形式 + 定数 (adts_offset: `adts(fourth='r3', optics='calc')` の物理的な項との差)
なので、リング全体の 1 回の組み立て (数秒) のあとは、任意の k に対する値と勾配が numpy の小さな行列計算 (ms) で厳密に (解析式 `resdrive.adts`・`rdt2nd.rdt4` と同じ値) 得られる。
モデルの検証: `SxtModel.check(k)` が解析式を直接呼んだ値と比べる (試験 `test_ddba_sextupole.py`)。

最適化の試作 (`SxtObjective`・`optimize`): 変数 = K2L (ファミリー)、等式制約 = クロマティシティ (ξx, ξy)、不等式 = |K2L| の上限。目的関数は
物理的な重み (基準の作用 J_ref・基準の運動量偏差 δ_ref と、共鳴までの距離の分母 |1−e^{iφ}|) だけで決め、追跡に合わせた調整はしない:
  W3_q = |h_q| (2J_ref)^{3/2} / |1 − e^{iφ_q}|   (3 次の非共鳴の正準変換の振幅 (無次元)。5 個の独立な調和 h21000・h30000・h10110・h10020・h10200)
  W4_q = |F4'_q| (2J_ref)^2 / |1 − e^{iφ_q}|     (4 次の非対角の調和、F4' は非共鳴の 3 次の摂動を含む正規形)
  A_c  = π |∂νx/∂Jx| J_ref²、2π |∂νx/∂Jy| J_ref²、π |∂νy/∂Jy| J_ref²  (c = xx, xy, yy。対角の h (2J)² で、ハミルトニアンの大きさ [rad/周]。W と同じ単位)
  D_q  = |h_q^{(p)}| (2J_ref)^{n/2} δ_ref^p / |1 − e^{iφ_q}|  (δ の 1 次の駆動項 h20001・h00201・h10002)
  目的 = Σ W3² + Σ W4² + Σ A_c² + Σ D_q²  (各項の重み a3・a4・aA・aD は既定 1。変えたときは README に書く)
限界: 5 次以上 (J の 2 次の ADTS、5・6 次の駆動項) と δ 微分の ADTS (∂²ν/∂δ² の 2 次の色収差) は目的に入っていない。
"""
import copy
import math
import re

import numpy as np

from . import rdt2nd as Q
from . import resdrive as R

TWO_PI = 2.0 * math.pi
CUBIC5 = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
CHROM3 = ('h20001', 'h00201', 'h10002')
DELTA_XI = ('h11001', 'h00111')


# ---------------------------------------------------------------------------------------------
# ファミリーと K2 の読み書き
# ---------------------------------------------------------------------------------------------
def sextupole_families(lat, names=None):
    """K2≠0 の六極のファミリー (素子名) の一覧 (最初に現れる順)、{名前: [素子の添字]}、K2L の配列 [1/m²] を返す。"""
    fam = {}
    for i, e in enumerate(lat.elements):
        if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0:
            fam.setdefault(e.name, []).append(i)
    order = list(names) if names is not None else list(fam)
    k = np.array([float(lat.elements[fam[n][0]].K2) * float(lat.elements[fam[n][0]].L) for n in order])
    return order, {n: fam[n] for n in order}, k


def with_k2l(lat, names, k):
    """ファミリー names の K2L を k [1/m²] にした Lattice の複製 (六極の K2 = K2L/L)。"""
    out = copy.deepcopy(lat)
    for n, kk in zip(names, k):
        for e in out.elements:
            if type(e).__name__ == 'Sextupole' and e.name == n:
                e.K2 = float(kk) / float(e.L)
    return out


def replace_k2_in_source(src, names, k):
    """SAD のラティス定義の文字列 src の中の SEXT 定義 `NAME =(L =... K2 =...)` の K2 (K2L) を k に置き換えた文字列 (同名の定義がいくつあっても全部)。
    追跡 (TrackSetup.from_source) に同じ K2 を渡すのに使う。置き換えが 0 回の名前があると ValueError。"""
    out = src
    for n, kk in zip(names, k):
        pat = re.compile(r'(\b%s\s*=\s*\(\s*L\s*=\s*[-+0-9.eE]+\s+K2\s*=\s*)([-+0-9.eE]+)' % re.escape(n))
        out, cnt = pat.subn(lambda m: m.group(1) + ('%.12g' % float(kk)), out)
        if cnt == 0:
            raise ValueError('replace_k2_in_source: 定義 %s が見つからない' % n)
    return out


# ---------------------------------------------------------------------------------------------
# 二次形式のフィット (偏極化)。f: k → 実または複素のベクトル。f が厳密に k の二次式であることを使う
# ---------------------------------------------------------------------------------------------
def quad_fit(fun, n):
    """f(k) = c + b·k + kᵀ S k (S は対称) の c (m,)、b (n, m)、S (n, n, m) を、1 + 2n + n(n−1)/2 回の f の呼び出しで厳密に求める。"""
    c = np.asarray(fun(np.zeros(n)))
    m = c.shape[0]
    cdt = c.dtype
    b = np.zeros((n, m), dtype=cdt)
    S = np.zeros((n, n, m), dtype=cdt)
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1.0
        fp, fm = np.asarray(fun(e)), np.asarray(fun(-e))
        b[i] = 0.5 * (fp - fm)
        S[i, i] = 0.5 * (fp + fm) - c
    for i in range(n):
        for j in range(i + 1, n):
            e = np.zeros(n)
            e[i] = e[j] = 1.0
            f = np.asarray(fun(e))
            S[i, j] = S[j, i] = 0.5 * (f - c - b[i] - b[j] - S[i, i] - S[j, j])
    return c, b, S


def quad_eval(cbS, k):
    c, b, S = cbS
    return c + k @ b + np.einsum('i,ijm,j->m', k, S, k)


def quad_jac(cbS, k):
    """∂f_m/∂k_i = b_im + 2 Σ_j S_ijm k_j。形 (m, n)"""
    c, b, S = cbS
    return (b + 2.0 * np.einsum('ijm,j->im', S, k)).T


# ---------------------------------------------------------------------------------------------
# モデル
# ---------------------------------------------------------------------------------------------
class SxtModel:
    """六極ファミリーの K2L に対する RDT・ADTS の多項式モデル。

    属性: names (ファミリー名)、k_ref (組み立てたときの K2L)、qx・qy、
      .h1(k) → 1 次の項 (CUBIC5 + CHROM3 + DELTA_XI、複素 dict。偏向磁石の 3 次の項を含む)
      .h3(k) → 3 次の 1 周の和 (MONOS3 の順 20 成分。rdt4 の `h3`)
      .h4(k) → 4 次の項 (MONOS4 の順 35 成分。rdt4 の total_arr)、.h4nf(k) → 正規形 F4' (h4 + nf_correction_h4)
      .adts(k) → {'xx','xy','yy'} (F4' の対角から。`resdrive.adts().total` と同じ)
      .xi(k) → (ξx, ξy) 六極込みの全クロマティシティ (自然クロマティシティ + 偏向磁石 + 六極)
    """

    def __init__(self, lat, momentum_eV=5e9, names=None, nslice=8, bend_nslice=None, quad_forms=True, fourth='r3'):
        self.lat = lat
        self.momentum_eV = float(momentum_eV)
        self.fourth = fourth
        self.names, self.fam, self.k_ref = sextupole_families(lat, names)
        F = len(self.names)
        self.F = F
        elements, rows = R._tables(lat)
        self.warnings = []
        self.qx, self.qy = float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI
        # ---- 1 次の項 (線形): resonance_terms の per_element をファミリーごとに集める
        terms1 = CUBIC5 + CHROM3 + DELTA_XI
        r1 = R._resonance_terms_cubic(lat, terms1, nslice=nslice, bend_cubic=True, bend_nslice=bend_nslice, momentum_eV=momentum_eV)
        self.terms1 = terms1
        nm_idx = {}
        for ir, nm in enumerate(r1.elem_names):
            nm_idx.setdefault(nm, []).append(ir)
        H = np.zeros((F, len(terms1)), dtype=complex)
        for f, n in enumerate(self.names):
            H[f] = r1.per_element[nm_idx[n]].sum(axis=0) / self.k_ref[f]
        self.H1 = H
        self.h1_const = np.array([r1.total[t] for t in r1.names]) - self.k_ref @ H
        self.r1 = r1
        # 自然クロマティシティ (六極なし) = 全体 − 六極込みの h/π: 全クロマティシティは resonance_terms の h11001/π・h00111/π と自然の和
        self.xi_nat = R.chromaticity_natural(lat)
        # ---- 4 次の項 (二次形式)
        nl = R.collect_nonlinear(elements, rows, 0.0, (0.0, 0.0), self.warnings)
        nl3 = [e for e in nl if 2 in e.strengths]
        nl4 = [e for e in nl if 3 in e.strengths]
        Asx, IGsx = Q.sext_units(nl3, 'analytic', nslice)
        bu = Q.bend_units(elements, rows, bend_nslice or R.BEND_NSLICE_DEFAULT, momentum_eV, self.warnings)
        sext_at = {e.idx: ie for ie, e in enumerate(nl3)}
        bend_at = {i: a for i, a in bu}
        fam_of = {}
        for f, n in enumerate(self.names):
            for i in self.fam[n]:
                fam_of[i] = f
        cA, cIG, tag = [], [], []
        for i in range(len(elements)):
            if i in sext_at:
                ie = sext_at[i]
                f = fam_of.get(i, F + 1)          # 名前に含めない六極は定数 (tag F+1) として扱う
                kf = self.k_ref[f] if f < F else 1.0
                if f < F:
                    cA.append(Asx[:, ie:ie + 1] / kf)
                    cIG.append(IGsx[:, :, ie:ie + 1] / kf ** 2)
                    tag.append(np.full(1, f))
                else:
                    cA.append(Asx[:, ie:ie + 1])
                    cIG.append(IGsx[:, :, ie:ie + 1])
                    tag.append(np.full(1, F))
            if i in bend_at:
                a = bend_at[i]
                cA.append(a)
                cIG.append(0.5 * a[:, None, :] * a[None, :, :])
                tag.append(np.full(a.shape[1], F))
        A = np.concatenate(cA, axis=1)
        IG = np.concatenate(cIG, axis=2)
        tag = np.concatenate(tag)
        self.A_tag, self.tag = A, tag
        # k_tag = (k_1, …, k_F, 1): 偏向磁石 (と名前のない六極) は強さ 1
        Af = [A[:, tag == t] for t in range(F + 1)]
        # 3 次の和 (線形): h3(k) = Σ_t k_t Σ_u A_u
        self.H3 = np.array([a.sum(axis=1) for a in Af])            # (F+1, 20)
        # 2 次の項: P_ts = Σ_{l ∈ t, e ∈ s, l > e} A_l ⊗ A_e + δ_ts Σ_u IG_u
        Cs = []
        for t in range(F + 1):
            Ct = np.zeros_like(A)
            At = np.where(tag[None, :] == t, A, 0.0)
            Ct = np.cumsum(At, axis=1) - At
            Cs.append(Ct)
        Wts = np.zeros((F + 1, F + 1, len(Q.MONOS4)), dtype=complex)
        for t in range(F + 1):
            At = np.where(tag[None, :] == t, A, 0.0)
            for s in range(F + 1):
                P = At @ Cs[s].T
                if t == s:
                    P = P + IG[:, :, tag == t].sum(axis=2)
                h = np.zeros(len(Q.MONOS4), dtype=complex)
                np.add.at(h, Q._BR_T, 0.5 * Q._BR_C * P[Q._BR_I, Q._BR_J])
                Wts[t, s] = h
        self.W4 = Wts
        # k に依らない 4 次の定数: 8 極・運動学・フリンジ・偏向磁石の 4 次 (rdt4 と同じ)
        r4 = Q.rdt4(lat, order2=False, momentum_eV=momentum_eV, nslice=nslice, bend_nslice=bend_nslice or R.BEND_NSLICE_DEFAULT)
        self.h4_const_fixed = r4.total_arr.copy()
        self.warnings += [w for w in r4.warnings if w not in self.warnings]
        self._forms = None
        # ADTS の定数部分: `resdrive.adts` の既定 (fourth='r3'・optics='calc'、物理的な項だけ) と、rdt4 の対角 (四極の台形則・偏向磁石の 4 次の項を追跡の写像から取り出した分を含む)
        # の差は、六極に依らない定数 (運動学・フリンジ・偏向磁石の 4 次の項の扱いの違い)。組み立てたときの K2L で 1 回だけ求めて、二次形式の対角に足す。
        raw = self._adts_of(self.h4nf(self.k_ref))
        ref = R.adts(lat, momentum_eV=momentum_eV, fourth=fourth, optics='calc').total
        self.adts_offset = np.array([float(ref['xx']), float(ref['xy']), float(ref['yy'])]) - raw
        self._ph3 = np.array([TWO_PI * ((q[0] - q[1]) * self.qx + (q[2] - q[3]) * self.qy) for q in Q.MONOS3])
        self._ph4 = np.array([TWO_PI * ((q[0] - q[1]) * self.qx + (q[2] - q[3]) * self.qy) for q in Q.MONOS4])
        if quad_forms:
            self.build_forms()

    # ---- 値 -------------------------------------------------------------------------------
    def _kt(self, k):
        return np.append(np.asarray(k, dtype=float), 1.0)

    def h1(self, k):
        v = self.h1_const + np.asarray(k, dtype=float) @ self.H1
        return {t: complex(v[i]) for i, t in enumerate(self.terms1)}

    def h3(self, k):
        return self._kt(k) @ self.H3

    def h4(self, k):
        kt = self._kt(k)
        return self.h4_const_fixed + np.einsum('t,tsm,s->m', kt, self.W4, kt)

    def h4nf(self, k):
        return self.h4(k) + Q.nf_correction_h4(self.h3(k), self.qx, self.qy)

    @staticmethod
    def _adts_of(h):
        i = Q._IDX4
        return np.array([4.0 * h[i[(2, 2, 0, 0)]].real / math.pi, 2.0 * h[i[(1, 1, 1, 1)]].real / math.pi, 4.0 * h[i[(0, 0, 2, 2)]].real / math.pi])

    def adts(self, k):
        """ADTS {'xx','xy','yy'} [1/(m rad)]。`resdrive.adts(fourth=self.fourth、optics='calc')` の .total と同じ (既定 'r3' = 物理的な項だけ)"""
        a = self._adts_of(self.h4nf(k)) + self.adts_offset
        return {'xx': a[0], 'xy': a[1], 'yy': a[2]}

    def xi(self, k):
        v = self.h1_const + np.asarray(k, dtype=float) @ self.H1
        ix, iy = self.terms1.index('h11001'), self.terms1.index('h00111')
        return (self.xi_nat[0] + v[ix].real / math.pi, self.xi_nat[1] + v[iy].real / math.pi)

    def xi_matrix(self):
        """ξ(k) = ξ0 + M k (M は 2 × F)。"""
        ix, iy = self.terms1.index('h11001'), self.terms1.index('h00111')
        M = np.vstack([self.H1[:, ix].real, self.H1[:, iy].real]) / math.pi
        xi0 = np.array([self.xi_nat[0] + self.h1_const[ix].real / math.pi, self.xi_nat[1] + self.h1_const[iy].real / math.pi])
        return xi0, M

    # ---- 二次形式 ---------------------------------------------------------------------------
    def build_forms(self):
        """h4nf (複素 35 成分) の二次形式を作る (偏極化)。"""
        self._forms = quad_fit(self.h4nf, self.F)
        return self._forms

    def h4nf_fast(self, k):
        return quad_eval(self._forms, np.asarray(k, dtype=float))

    # ---- 検証 -------------------------------------------------------------------------------
    def check(self, k=None, lat=None):
        """モデルと、解析式の直接の呼び出し (resonance_terms・rdt4・adts) の相対差を返す (k = None なら組み立てたときの K2L)。"""
        k = self.k_ref if k is None else np.asarray(k, dtype=float)
        lat2 = with_k2l(self.lat, self.names, k) if lat is None else lat
        out = {}
        r1 = R._resonance_terms_cubic(lat2, self.terms1, bend_cubic=True, momentum_eV=self.momentum_eV)
        m1 = self.h1(k)
        out['h1'] = max(abs(m1[t] - r1.total[t]) / max(abs(r1.total[t]), 1e-30) for t in self.terms1)
        r4 = Q.rdt4(lat2, order2=True, momentum_eV=self.momentum_eV)
        out['h4'] = float(np.max(np.abs(self.h4(k) - r4.total_arr)) / np.max(np.abs(r4.total_arr)))
        out['h3'] = float(np.max(np.abs(self.h3(k) - r4.h3)) / np.max(np.abs(r4.h3)))
        a = R.adts(lat2, momentum_eV=self.momentum_eV, fourth=self.fourth, optics='calc').total
        am = self.adts(k)
        out['adts'] = max(abs(am[c] - float(a[c])) / max(abs(float(a[c])), 1e-30) for c in ('xx', 'xy', 'yy'))
        if self._forms is not None:
            out['forms'] = float(np.max(np.abs(self.h4nf_fast(k) - self.h4nf(k))) / np.max(np.abs(self.h4nf(k))))
        return out


# ---------------------------------------------------------------------------------------------
# 目的関数
# ---------------------------------------------------------------------------------------------
class SxtObjective:
    """物理的な重みの目的関数 (docs 先頭の説明)。値・勾配・項ごとの内訳。"""

    def __init__(self, model, J_ref=1e-6, delta_ref=0.01, weights=None, terms=('w3', 'w4', 'adts', 'chrom')):
        self.m = model
        self.J = float(J_ref)
        self.d = float(delta_ref)
        self.wt = dict(w3=1.0, w4=1.0, adts=1.0, chrom=1.0)
        if weights:
            self.wt.update(weights)
        self.terms = tuple(terms)
        m = model
        F = m.F
        qx, qy = m.qx, m.qy

        def den(j, k_, l, m_):
            ph = TWO_PI * ((j - k_) * qx + (l - m_) * qy)
            return abs(1.0 - np.exp(1j * ph))

        # w3: CUBIC5 の h1 の添字
        self.i3 = [m.terms1.index(t) for t in CUBIC5]
        self.s3 = np.array([(2 * self.J) ** 1.5 / den(*[int(c) for c in t[1:5]]) for t in CUBIC5])
        self.ic = [m.terms1.index(t) for t in CHROM3]
        self.sc = np.array([(2 * self.J) ** (sum(int(c) for c in t[1:5]) / 2.0) * self.d ** int(t[5]) / den(*[int(c) for c in t[1:5]]) for t in CHROM3])
        # w4: 非対角 (共鳴の分母が小さすぎるものと対角 (k=0 の調和) は除く)
        sel, sc4 = [], []
        for i, q in enumerate(Q.MONOS4):
            if q[0] == q[1] and q[2] == q[3]:
                continue
            dn = den(*q)
            if dn < 1e-3:
                continue
            sel.append(i)
            sc4.append((2 * self.J) ** 2 / dn)
        self.i4 = np.array(sel)
        self.s4 = np.array(sc4)
        # ADTS の項: ハミルトニアン (1 周あたり [rad]) の対角部分 h22000 (2J)² = π (∂νx/∂Jx) J²、h11110 (2J)² = 2π (∂νx/∂Jy) J²、h00220 (2J)² = π (∂νy/∂Jy) J²
        self.sA = np.array([math.pi, TWO_PI, math.pi]) * self.J ** 2

    def parts(self, k):
        k = np.asarray(k, dtype=float)
        m = self.m
        v = m.h1_const + k @ m.H1
        out = {}
        out['w3'] = float(np.sum(np.abs(v[self.i3] * self.s3) ** 2))
        out['chrom'] = float(np.sum(np.abs(v[self.ic] * self.sc) ** 2))
        h = m.h4nf_fast(k)
        out['w4'] = float(np.sum(np.abs(h[self.i4] * self.s4) ** 2))
        a = m._adts_of(h) + m.adts_offset
        out['adts'] = float(np.sum((a * self.sA) ** 2))
        return out

    def value(self, k):
        p = self.parts(k)
        return sum(self.wt[t] * p[t] for t in self.terms)

    def grad(self, k):
        k = np.asarray(k, dtype=float)
        m = self.m
        g = np.zeros(m.F)
        v = m.h1_const + k @ m.H1
        if 'w3' in self.terms:
            for idx, s in zip(self.i3, self.s3):
                g += self.wt['w3'] * 2.0 * s * s * np.real(np.conj(v[idx]) * m.H1[:, idx])
        if 'chrom' in self.terms:
            for idx, s in zip(self.ic, self.sc):
                g += self.wt['chrom'] * 2.0 * s * s * np.real(np.conj(v[idx]) * m.H1[:, idx])
        if 'w4' in self.terms or 'adts' in self.terms:
            h = m.h4nf_fast(k)
            Jc = quad_jac(m._forms, k)                       # (35, F) 複素
            if 'w4' in self.terms:
                for i, s in zip(self.i4, self.s4):
                    g += self.wt['w4'] * 2.0 * s * s * np.real(np.conj(h[i]) * Jc[i])
            if 'adts' in self.terms:
                ii = Q._IDX4
                comps = [((2, 2, 0, 0), 4.0 / math.pi), ((1, 1, 1, 1), 2.0 / math.pi), ((0, 0, 2, 2), 4.0 / math.pi)]
                for n, (q, c) in enumerate(comps):
                    a = c * h[ii[q]].real + m.adts_offset[n]
                    da = c * Jc[ii[q]].real
                    g += self.wt['adts'] * 2.0 * self.sA[n] ** 2 * a * da
        return g


def optimize(obj, k0, xi_target, kmax, method='SLSQP', maxiter=300, ftol=1e-14):
    """ξ = xi_target の等式制約と |K2L| ≤ kmax の範囲で obj.value を最小化 (scipy の SLSQP。勾配は解析)。戻り (k, 結果の OptimizeResult)。"""
    from scipy.optimize import minimize
    xi0, M = obj.m.xi_matrix()
    tgt = np.asarray(xi_target, dtype=float)
    cons = [{'type': 'eq', 'fun': lambda k: xi0 + M @ k - tgt, 'jac': lambda k: M}]
    bnds = [(-kmax, kmax)] * obj.m.F if np.isscalar(kmax) else [(-a, a) for a in kmax]
    sc = max(obj.value(k0), 1e-30)
    res = minimize(lambda k: obj.value(k) / sc, k0, jac=lambda k: obj.grad(k) / sc, method=method, bounds=bnds, constraints=cons,
                   options=dict(maxiter=maxiter, ftol=ftol))
    return res.x, res


def project_to_xi(model, k, xi_target):
    """k を、ξ = xi_target を満たす最も近い点 (最小ノルムの補正) に射影する。"""
    xi0, M = model.xi_matrix()
    r = np.asarray(xi_target) - (xi0 + M @ k)
    return k + M.T @ np.linalg.solve(M @ M.T, r)
