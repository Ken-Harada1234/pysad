# -*- coding: utf-8 -*-
"""Prototype polynomial model of RDT and ADTS versus sextupole-family K2, and sextupole optimization (R7) on it.
[試作] 6 極ファミリーの K2 に対する RDT・ADTS の解析式の「多項式モデル」と、それを目的関数にした6 極の最適化 (R7)。

考え方: 6 極の強さ k = (K2L_1, ..., K2L_F) (ファミリーごとの積分強さ) に対して、
  - 3 次の共鳴駆動項 h_jklm (1 次) とクロマティシティ・δ の項 (p ≥ 1) は k に **線形**:  h(k) = h_B + Σ_f k_f H_f
  - 4 次の項 h_jklm (写像の 4 次の生成関数、`rdt2nd.rdt4`) の6 極の 2 次の部分 (Magnus の第 2 項) は k の **二次形式**、
    6 極 × 偏向磁石の交差項は k に線形、八極・運動学・フリンジ・偏向磁石の 4 次は k に依らない定数
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

v2.4.1 の拡張 (項の指定): 項ごとに 'objective' / 'constraint' / 'ignore'、重み、目標値を指定する `SxtSpec` (辞書・文字列・JSON から作れる。SAD の台本のオプションからも
同じ形で渡せる) と、それを使う目的関数 `SxtSpecObjective`、作動点の周りの共鳴の一覧 `near_resonances`・警告 `resonance_warnings`。
既定の指定 `legacy_spec()` は従来の `SxtObjective` と同じ値・勾配 (1e-12 以内) になる (従来のクラスはそのまま残してある)。
"""
import copy
import dataclasses
import json
import math
import re
import warnings
from dataclasses import dataclass, field

import numpy as np

from . import rdt2nd as Q
from . import resdrive as R

TWO_PI = 2.0 * math.pi
CUBIC5 = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
CHROM3 = ('h20001', 'h00201', 'h10002')
DELTA_XI = ('h11001', 'h00111')
ADTS_TERMS = ('adts_xx', 'adts_xy', 'adts_yy')
XI_TERMS = ('xi_x', 'xi_y')


# ---------------------------------------------------------------------------------------------
# ファミリーと K2 の読み書き
# ---------------------------------------------------------------------------------------------
def sextupole_families(lat, names=None):
    """Return the sextupole families with K2 != 0 (names, element indices) and the K2L array [1/m^2].
    K2≠0 の6 極のファミリー (素子名) の一覧 (最初に現れる順)、{名前: [素子の添字]}、K2L の配列 [1/m²] を返す。"""
    fam = {}
    for i, e in enumerate(lat.elements):
        if type(e).__name__ == 'Sextupole' and float(e.K2) != 0.0:
            fam.setdefault(e.name, []).append(i)
    order = list(names) if names is not None else list(fam)
    k = np.array([float(lat.elements[fam[n][0]].K2) * float(lat.elements[fam[n][0]].L) for n in order])
    return order, {n: fam[n] for n in order}, k


def with_k2l(lat, names, k):
    """Return a copy of the Lattice with the K2L of the families `names` set to k [1/m^2].
    ファミリー names の K2L を k [1/m²] にした Lattice の複製 (6 極の K2 = K2L/L)。"""
    out = copy.deepcopy(lat)
    for n, kk in zip(names, k):
        for e in out.elements:
            if type(e).__name__ == 'Sextupole' and e.name == n:
                e.K2 = float(kk) / float(e.L)
    return out


def replace_k2_in_source(src, names, k):
    """Replace K2 (K2L) of SEXT definitions in a SAD lattice source string with k, for every definition of the same name.
    SAD のラティス定義の文字列 src の中の SEXT 定義 `NAME =(L =... K2 =...)` の K2 (K2L) を k に置き換えた文字列 (同名の定義がいくつあっても全部)。
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
    """Exactly fit f(k) = c + b.k + k^T S k (S symmetric), returning c, b and S from 1 + 2n + n(n-1)/2 evaluations of f.
    f(k) = c + b·k + kᵀ S k (S は対称) の c (m,)、b (n, m)、S (n, n, m) を、1 + 2n + n(n−1)/2 回の f の呼び出しで厳密に求める。"""
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
    """Return the Jacobian df_m/dk_i = b_im + 2 sum_j S_ijm k_j, of shape (m, n).
    ∂f_m/∂k_i = b_im + 2 Σ_j S_ijm k_j。形 (m, n)"""
    c, b, S = cbS
    return (b + 2.0 * np.einsum('ijm,j->im', S, k)).T


# ---------------------------------------------------------------------------------------------
# モデル
# ---------------------------------------------------------------------------------------------
class SxtModel:
    """Polynomial model of RDT and ADTS as functions of the sextupole-family K2L.
    6 極ファミリーの K2L に対する RDT・ADTS の多項式モデル。

    属性: names (ファミリー名)、k_ref (組み立てたときの K2L)、qx・qy、
      .h1(k) → 1 次の項 (CUBIC5 + CHROM3 + DELTA_XI、複素 dict。偏向磁石の 3 次の項を含む)
      .h3(k) → 3 次の 1 周の和 (MONOS3 の順 20 成分。rdt4 の `h3`)
      .h4(k) → 4 次の項 (MONOS4 の順 35 成分。rdt4 の total_arr)、.h4nf(k) → 正規形 F4' (h4 + nf_correction_h4)
      .adts(k) → {'xx','xy','yy'} (F4' の対角から。`resdrive.adts().total` と同じ)
      .xi(k) → (ξx, ξy) 6 極込みの全クロマティシティ (自然クロマティシティ + 偏向磁石 + 6 極)
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
        # 台ごとの K2L のランダム誤差 (相対 rms ε、台ごとに独立) が作る 1 次の項の rms: rms(δh_t)² = ε² Σ_f k_f² B[f, t]、B[f, t] = Σ_{i∈f} |a_it|² (a_it = K2L 1 あたりの 1 台の寄与)。
        # 位相の打ち消しがない (台ごとに独立) ので、強さの二乗の重み。errB は (F, len(terms1)) の実数
        self.errB = np.zeros((F, len(terms1)))
        for f, n in enumerate(self.names):
            a = r1.per_element[nm_idx[n]] / self.k_ref[f]
            self.errB[f] = np.sum(np.abs(a) ** 2, axis=0)
        # 自然クロマティシティ (6 極なし) = 全体 − 6 極込みの h/π: 全クロマティシティは resonance_terms の h11001/π・h00111/π と自然の和
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
                f = fam_of.get(i, F + 1)          # 名前に含めない6 極は定数 (tag F+1) として扱う
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
        # k_tag = (k_1, …, k_F, 1): 偏向磁石 (と名前のない6 極) は強さ 1
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
        # の差は、6 極に依らない定数 (運動学・フリンジ・偏向磁石の 4 次の項の扱いの違い)。組み立てたときの K2L で 1 回だけ求めて、二次形式の対角に足す。
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
        """Compute ADTS {xx, xy, yy} [1/(m rad)], equal to the total of `resdrive.adts` (default r3 physical terms only).
        ADTS {'xx','xy','yy'} [1/(m rad)]。`resdrive.adts(fourth=self.fourth、optics='calc')` の .total と同じ (既定 'r3' = 物理的な項だけ)"""
        a = self._adts_of(self.h4nf(k)) + self.adts_offset
        return {'xx': a[0], 'xy': a[1], 'yy': a[2]}

    def xi(self, k):
        v = self.h1_const + np.asarray(k, dtype=float) @ self.H1
        ix, iy = self.terms1.index('h11001'), self.terms1.index('h00111')
        return (self.xi_nat[0] + v[ix].real / math.pi, self.xi_nat[1] + v[iy].real / math.pi)

    def err_rms(self, k, eps_rms):
        """Compute the expected rms of the first-order terms for random K2L errors (independent per magnet, relative rms eps_rms).
        K2L のランダム誤差 (台ごとに独立、相対 rms eps_rms) による 1 次の項の rms の期待値 {項: rms} (CUBIC5・CHROM3・δ の ξ の項も)。"""
        k = np.asarray(k, dtype=float)
        v = eps_rms * np.sqrt((k ** 2) @ self.errB)
        return {t: float(v[i]) for i, t in enumerate(self.terms1)}

    def xi_matrix(self):
        """Return the chromaticity xi(k) = xi0 + M k, with M of size 2 x F.
        ξ(k) = ξ0 + M k (M は 2 × F)。"""
        ix, iy = self.terms1.index('h11001'), self.terms1.index('h00111')
        M = np.vstack([self.H1[:, ix].real, self.H1[:, iy].real]) / math.pi
        xi0 = np.array([self.xi_nat[0] + self.h1_const[ix].real / math.pi, self.xi_nat[1] + self.h1_const[iy].real / math.pi])
        return xi0, M

    # ---- 二次形式 ---------------------------------------------------------------------------
    def build_forms(self):
        """Build the quadratic forms of the 35-component complex h4nf using the polarization identity.
        h4nf (複素 35 成分) の二次形式を作る (偏極化)。"""
        self._forms = quad_fit(self.h4nf, self.F)
        return self._forms

    def h4nf_fast(self, k):
        return quad_eval(self._forms, np.asarray(k, dtype=float))

    # ---- 検証 -------------------------------------------------------------------------------
    def check(self, k=None, lat=None):
        """Return the relative difference between the model and direct calls of the analytic formulas (resonance_terms, rdt4, adts).
        モデルと、解析式の直接の呼び出し (resonance_terms・rdt4・adts) の相対差を返す (k = None なら組み立てたときの K2L)。"""
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
    """Objective function with physical weights (value, gradient and per-term breakdown).
    物理的な重みの目的関数 (docs 先頭の説明)。値・勾配・項ごとの内訳。"""

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
    """Minimize obj.value under the equality xi = xi_target and bounds |K2L| <= kmax (scipy SLSQP, analytic gradient).
    ξ = xi_target の等式制約と |K2L| ≤ kmax の範囲で obj.value を最小化 (scipy の SLSQP。勾配は解析)。戻り (k, 結果の OptimizeResult)。"""
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
    """Project k onto the nearest point satisfying xi = xi_target (minimum-norm correction).
    k を、ξ = xi_target を満たす最も近い点 (最小ノルムの補正) に射影する。"""
    xi0, M = model.xi_matrix()
    r = np.asarray(xi_target) - (xi0 + M @ k)
    return k + M.T @ np.linalg.solve(M @ M.T, r)


# ---------------------------------------------------------------------------------------------
# v2.4.1: 項の指定 (objective / constraint / ignore、重み、目標値)
# ---------------------------------------------------------------------------------------------
MODES = ('objective', 'constraint', 'ignore')
GROUPS = {
    'w3': CUBIC5,                      # 3 次の 5 個の独立な調和
    'chrom': CHROM3,                   # δ の 1 次の駆動項 (h20001・h00201・h10002)
    'adts': ADTS_TERMS,
    'xi': XI_TERMS,
    'err3': tuple('err_' + t for t in CUBIC5),              # 誤差感度: 台ごとのランダム誤差が作る 3 次の項の rms (期待値)
    'errchrom': tuple('err_' + t for t in CHROM3),
}
# 4 次の項: 名前は 'h' + jklm の 4 桁 + '0' (δ の次数。例 h31000)。対角 (h22000・h11110・h00220) は ADTS なので adts_* で指定する。
# グループ 'w4' は非対角の全調和 (共役の対を両方含む。従来の SxtObjective の w4 と同じ)。個別に指定した 4 次の項は 1 個ごとの重み (共役の対のもう一方は別名で別に指定)。


def _name4(q):
    return 'h%d%d%d%d0' % tuple(q)


def _den(q, qx, qy):
    """共鳴までの分母 |1 − e^{iφ}|、φ = 2π((j−k)νx + (l−m)νy)。作動点のチューンから自動で出す。"""
    ph = TWO_PI * ((q[0] - q[1]) * qx + (q[2] - q[3]) * qy)
    return float(abs(1.0 - np.exp(1j * ph)))


def _w4_names(qx, qy, dmin=1e-3):
    out = []
    for q in Q.MONOS4:
        if q[0] == q[1] and q[2] == q[3]:
            continue
        if _den(q, qx, qy) < dmin:
            continue
        out.append(_name4(q))
    return out


@dataclass
class TermSpec:
    """Specification of one term: mode (objective, constraint or ignore), target value, weight and tolerance.
    1 項の指定。mode: 'objective' (目的関数に入れる: weight × |scale (q − target)|²)、'constraint' (制約。tol=None なら q = target の等式、tol があれば |q − target| ≤ tol)、'ignore'。
    target: 目標値 (RDT の項は複素数も可。単位は項の単位: h は RDT の単位、ADTS は [1/(m rad)]、xi は無次元)。None (既定) は RDT・ADTS では 0 (0 へ近づける)、xi では optimize の xi_target。
    scale: 目的関数の規格化 (None = 自動: 共鳴の分母と基準の作用 J_ref・δ_ref で決まる物理的な重み。SxtObjective と同じ)。"""
    mode: str = 'objective'
    weight: float = 1.0
    target: complex = None
    tol: float = None
    scale: float = None

    def __post_init__(self):
        if self.mode not in MODES:
            raise ValueError('mode は %s のどれか: %r' % (MODES, self.mode))
        if self.weight < 0:
            raise ValueError('weight は 0 以上')


def _term_from_any(v):
    if isinstance(v, TermSpec):
        return dataclasses.replace(v)
    if isinstance(v, str):
        return TermSpec(mode=v.lower())
    if isinstance(v, (int, float)):
        return TermSpec('objective', float(v))
    if isinstance(v, (tuple, list)):
        f = ['mode', 'weight', 'target', 'tol', 'scale']
        kw = {}
        for i, x in enumerate(v):
            kw[f[i]] = str(x).lower() if i == 0 else x
        return TermSpec(**kw)
    if isinstance(v, dict):
        kw = {}
        for key, x in v.items():
            kk = str(key).lower()
            if kk not in ('mode', 'weight', 'target', 'tol', 'scale'):
                raise ValueError('指定のキー %r は不明' % key)
            if kk == 'target' and isinstance(x, (list, tuple)) and len(x) == 2:
                x = complex(x[0], x[1])
            kw[kk] = str(x).lower() if kk == 'mode' else x
        return TermSpec(**kw)
    raise TypeError('項の指定の型が不明: %r' % (v,))


@dataclass
class SxtSpec:
    """Specification of all terms as {name: TermSpec}, with individual terms or groups (w3, chrom, w4, adts, xi).
    項の指定の全体。terms: {名前: TermSpec}。名前は個別の項 (h21000・h30000・h10110・h10020・h10200、h20001・h00201・h10002、4 次 h31000 など、adts_xx・adts_xy・adts_yy、xi_x・xi_y) か
    グループ ('w3'・'chrom'・'w4'・'adts'・'xi')。個別の指定はグループより優先。書かれていない項は 'ignore'。

    作り方 (`SxtSpec.parse` が全部受ける。SAD の台本のオプションと同じ形):
      辞書     {'adts': 'objective', 'h10020': {'mode': 'objective', 'weight': 4}, 'xi': ('constraint', 1, None)}
      文字列   'adts=objective, h10020=objective:4:0, xi=constraint'   (名前=mode[:weight[:target[:tol]]]、区切りは , か ;)
      JSON     辞書の JSON 文字列
      規則の列 [('adts', 'objective'), ('xi_x', ('constraint', 1, 0.0))]  (SAD の Rule の列を変換した形)"""
    terms: dict = field(default_factory=dict)

    @staticmethod
    def parse(x):
        if x is None:
            return legacy_spec()
        if isinstance(x, SxtSpec):
            return SxtSpec({n: dataclasses.replace(t) for n, t in x.terms.items()})
        if isinstance(x, str):
            xs = x.strip()
            if xs.startswith('{'):
                return SxtSpec.parse(json.loads(xs))
            d = {}
            for item in re.split(r'[;,]\s*', xs):
                if not item.strip():
                    continue
                name, _, rest = item.partition('=')
                f = [p.strip() for p in rest.split(':')] if rest else ['objective']
                kw = {'mode': f[0].lower()}
                if len(f) > 1 and f[1] != '':
                    kw['weight'] = float(f[1])
                if len(f) > 2 and f[2] != '':
                    kw['target'] = complex(f[2].replace('i', 'j')) if ('j' in f[2] or 'i' in f[2]) else float(f[2])
                if len(f) > 3 and f[3] != '':
                    kw['tol'] = float(f[3])
                d[name.strip()] = TermSpec(**kw)
            return SxtSpec(d)
        if isinstance(x, dict):
            return SxtSpec({str(n): _term_from_any(v) for n, v in x.items()})
        if isinstance(x, (list, tuple)):
            return SxtSpec({str(n): _term_from_any(v) for n, v in x})
        raise TypeError('SxtSpec.parse: 型が不明 %r' % type(x))

    def to_dict(self):
        out = {}
        for n, t in self.terms.items():
            d = dict(mode=t.mode, weight=t.weight, target=t.target)
            if t.tol is not None:
                d['tol'] = t.tol
            if t.scale is not None:
                d['scale'] = t.scale
            if isinstance(t.target, complex):
                d['target'] = [t.target.real, t.target.imag]
            out[n] = d
        return out

    def to_json(self):
        return json.dumps(self.to_dict(), ensure_ascii=False)

    def expand(self, qx, qy):
        """Expand groups into individual terms as {name: TermSpec}, with individual specifications taking precedence.
        グループを個別の項に展開した {名前: TermSpec} (個別の指定が優先。'ignore' は残す)。"""
        out = {}
        for n, t in self.terms.items():
            if n in GROUPS or n == 'w4':
                names = _w4_names(qx, qy) if n == 'w4' else GROUPS[n]
                for nn in names:
                    out[nn] = dataclasses.replace(t)
        for n, t in self.terms.items():
            if n in GROUPS or n == 'w4':
                continue
            out[n] = dataclasses.replace(t)
        for n in out:
            _check_name(n)
        return out


def _check_name(n):
    if n in CUBIC5 or n in CHROM3 or n in ADTS_TERMS or n in XI_TERMS:
        return
    if n.startswith('err_') and (n[4:] in CUBIC5 or n[4:] in CHROM3):
        return
    if re.fullmatch(r'h\d{4}0', n):
        q = tuple(int(c) for c in n[1:5])
        if q in Q._IDX4:
            if q[0] == q[1] and q[2] == q[3]:
                raise ValueError('%s は対角 (ADTS) の項。adts_xx・adts_xy・adts_yy で指定する' % n)
            return
    raise ValueError('項の名前が不明: %s (使える名前: %s、4 次 h31000 など、%s、%s、グループ %s と w4)' % (n, CUBIC5 + CHROM3, ADTS_TERMS, XI_TERMS, tuple(GROUPS)))


def legacy_spec(weights=None, terms=('w3', 'w4', 'adts', 'chrom')):
    """Return the legacy specification equivalent to `SxtObjective(weights=..., terms=...)`.
    従来の `SxtObjective(weights=…、terms=…)` と同じ指定 (グループごとの重みで目的関数、ξ は optimize の xi_target の等式制約)。"""
    wt = dict(w3=1.0, w4=1.0, adts=1.0, chrom=1.0)
    if weights:
        wt.update(weights)
    d = {g: TermSpec('objective', wt[g]) for g in terms}
    d['xi'] = TermSpec('constraint', 1.0, None)           # target None: optimize の xi_target を使う
    return SxtSpec(d)


# ---------------------------------------------------------------------------------------------
# 作動点の周りの共鳴
# ---------------------------------------------------------------------------------------------
def _drive_terms():
    """{(a, b): 共鳴 aνx + bνy = p を駆動する項}: a = j−k、b = l−m を ± 同一視 (正規化: 最初の 0 でない成分が正)。"""
    tab = {}

    def add(q, name, kind):
        a, b = q[0] - q[1], q[2] - q[3]
        if a < 0 or (a == 0 and b < 0):
            a, b = -a, -b
        if a == 0 and b == 0:
            return
        tab.setdefault((a, b), []).append((name, kind))
    for q in Q.MONOS3:
        if (q[2] + q[3]) % 2:                      # y の奇数次の項は midplane 対称で 6 極・八極では駆動されない
            continue
        add(q, _name4(q), '3 次 (6 極の 1 次)')
    for t in CHROM3:
        add([int(c) for c in t[1:5]], t, 'δ の 1 次 (色収差の駆動)')
    for q in Q.MONOS4:
        if (q[2] + q[3]) % 2:
            continue
        add(q, _name4(q), '4 次 (6 極の 2 次・八極)')
    return tab


def near_resonances(qx, qy, max_order=4, radius=0.05):
    """List resonances a*nux + b*nuy = p up to max_order within radius of the working point (qx, qy) in the tune plane, nearest first.
    作動点 (qx, qy) の周りの共鳴 a νx + b νy = p (|a|+|b| = 次数 ≤ max_order) のうち、チューン面の距離が radius 以内のものの一覧 (近い順)。
    各要素: dict(order, a, b, p, dist (チューン面での距離 |aνx+bνy−p|/√(a²+b²))、detune (aνx+bνy−p、符号つき)、terms ([(項の名前, 種類)])、sextupole_driven (駆動項があるか)、
    midplane_forbidden (b が奇数: 水平面に対称な磁石 (6 極・八極・K2L のランダム誤差を含む) では y の次数が偶数だけなので、何次でも駆動されない線。連成・スキュー・垂直の偏向があるときだけ駆動される))。
    「どの駆動項が効くか」は 3 次の h_jklm (6 極の 1 次)、δ の 1 次の項、4 次の h_jklm (6 極の 2 次・八極) のうち (j−k, l−m) = ±(a, b) のもの。
    3 次・4 次の項の名前は h_jklm の 5 桁 (最後が δ の次数)。駆動項は y の次数 (l+m) が偶数のもの (midplane 対称)。駆動項が空なら、その共鳴は 6 極の 1 次・2 次の調和では駆動されない (5 次以上は 6 極の 3 次以降、または八極・十極)。"""
    tab = _drive_terms()
    out = []
    n = int(max_order)
    for a in range(0, n + 1):
        for b in range(-n, n + 1):
            if a == 0 and b <= 0:
                continue
            order = abs(a) + abs(b)
            if order > n:
                continue
            v = a * qx + b * qy
            p = round(v)
            det = v - p
            dist = abs(det) / math.hypot(a, b)
            if dist > radius:
                continue
            terms = tab.get((a, b), [])
            out.append(dict(order=order, a=a, b=b, p=int(p), dist=float(dist), detune=float(det), terms=list(terms),
                            sextupole_driven=bool(terms), midplane_forbidden=bool(b % 2)))
    out.sort(key=lambda d: (d['dist'], d['order']))
    return out


def format_near_resonances(res):
    rows = ['次数  共鳴 a νx + b νy = p     距離(チューン面)  駆動項']
    for d in res:
        tn = ', '.join('%s[%s]' % (nm, kd.split(' ')[0]) for nm, kd in d['terms']) or ('(b が奇数: 水平面対称では駆動されない)' if d['midplane_forbidden'] else '(3・4 次の駆動項なし。5 次以上では駆動されうる)')
        rows.append('%4d  %+d νx %+d νy = %+d   %.4f           %s' % (d['order'], d['a'], d['b'], d['p'], d['dist'], tn))
    return '\n'.join(rows)


def resonance_warnings(qx, qy, max_order=4, thresholds=None):
    """Return warnings for working points near low-order resonances; integer and half-integer resonances must always be checked.
    作動点が低次の共鳴に近いときの警告文の一覧。しきい値 (チューン面の距離): 1〜2 次 0.03、3 次 0.02、4 次 0.01 (既定)。作動点の整数・半整数共鳴は必ず確認する。
    b が奇数の線 (水平面に対称な磁石では駆動されない) は除く。"""
    thr = thresholds or {1: 0.03, 2: 0.03, 3: 0.02, 4: 0.01}
    out = []
    for d in near_resonances(qx, qy, max_order=max_order, radius=max(thr.values())):
        if d['midplane_forbidden'] and not d['terms']:
            continue                  # b が奇数の線は、水平面対称の磁石では何次でも駆動されない (連成・スキューがあれば別)
        if d['dist'] <= thr.get(d['order'], 0.0):
            nm = ', '.join(t[0] for t in d['terms']) or '駆動項なし'
            out.append('作動点 (%.4f, %.4f) が %d 次共鳴 %+dνx%+dνy=%+d に近い (距離 %.4f。駆動項: %s)。駆動項があれば分母 |1−e^{iφ}| が小さく、その項の重みが大きくなる (なければ 6 極の 1 次・2 次では駆動されない共鳴)。'
                       % (qx, qy, d['order'], d['a'], d['b'], d['p'], d['dist'], nm))
    return out


# ---------------------------------------------------------------------------------------------
# 指定つきの目的関数
# ---------------------------------------------------------------------------------------------
class SxtSpecObjective:
    """Objective function following a `SxtSpec`: value, gradient, per-term breakdown (`table`) and constraints (`scipy_constraints`).
    `SxtSpec` に従う目的関数。値 Σ_{objective} w |scale (q − target)|²、勾配、項ごとの内訳 (`table`)、制約 (`scipy_constraints`)。
    scale は既定で SxtObjective と同じ物理的な重み (共鳴の分母は model の作動点のチューンから自動、J_ref・δ_ref)。

    ξ (xi_x・xi_y): 'constraint' なら ξ = target の等式 (target=None なら optimize の xi_target)、'objective' なら w (ξ − target)² (規格化なし。無次元)。
    制約の扱い: 3 次の RDT・δ の項・ξ の等式は k の 1 次 (線形) なので、実部・虚部の行を作って従属な行を落とす (鏡面対称の軸から始めるリングは複素 1 個が実数 1 個)。
    4 次の項・ADTS の制約は非線形 (二次形式) で、等式は (q−target)/参照値、tol つきは 1 − (|q−target|/tol)² ≥ 0 の形。"""

    def __init__(self, model, spec=None, J_ref=1e-6, delta_ref=0.01, warn=True, err_rms=0.03):
        self.m = model
        self.J = float(J_ref)
        self.d = float(delta_ref)
        self.err_rms = float(err_rms)           # 'err_*' 項 (誤差感度) の K2L の相対 rms ε
        self.spec = SxtSpec.parse(spec) if spec is not None else legacy_spec()
        self.warnings = []
        m = model
        qx, qy = m.qx, m.qy
        ex = self.spec.expand(qx, qy)
        self.items = []
        for n, t in ex.items():
            if t.mode == 'ignore':
                continue
            it = dict(name=n, t=t, den=float('nan'))
            base = n[4:] if n.startswith('err_') else n
            if base in CUBIC5 or base in CHROM3:
                q = tuple(int(c) for c in base[1:5])
                it['kind'] = 'sens' if n.startswith('err_') else 'h1'
                it['idx'] = m.terms1.index(base)
                den = _den(q, qx, qy)
                n = base
                if n in CUBIC5:
                    s = (2 * self.J) ** 1.5 / den
                else:
                    s = (2 * self.J) ** (sum(q) / 2.0) * self.d ** int(n[5]) / den
                it['den'] = den
            elif re.fullmatch(r'h\d{4}0', n):
                q = tuple(int(c) for c in n[1:5])
                den = _den(q, qx, qy)
                if den < 1e-3:
                    self.warnings.append('項 %s は共鳴の分母 %.2e < 1e-3 で、扱わない (作動点が共鳴上)' % (n, den))
                    continue
                it['kind'] = 'h4'
                it['idx'] = Q._IDX4[q]
                s = (2 * self.J) ** 2 / den
                it['den'] = den
            elif n in ADTS_TERMS:
                i = ADTS_TERMS.index(n)
                it['kind'] = 'adts'
                it['idx'] = i
                s = [math.pi, TWO_PI, math.pi][i] * self.J ** 2
            else:
                it['kind'] = 'xi'
                it['idx'] = XI_TERMS.index(n)
                s = 1.0
            it['s'] = float(t.scale) if t.scale is not None else float(s)
            self.items.append(it)
            if it['kind'] in ('h1', 'h4') and t.weight > 0:
                dd = abs(((q[0] - q[1]) * qx + (q[2] - q[3]) * qy + 0.5) % 1.0 - 0.5)
                if dd < 0.03:
                    self.warnings.append('項 %s の共鳴 (%d, %d) までの距離 (チューンの小数部) が %.4f と近い。分母 %.3f' % (n, q[0] - q[1], q[2] - q[3], dd, it['den']))
        self.warnings = resonance_warnings(qx, qy) + self.warnings
        if warn:
            for w in self.warnings:
                warnings.warn(w, RuntimeWarning, stacklevel=2)
        self.xi0, self.M = m.xi_matrix()
        self.obj_items = [it for it in self.items if it['t'].mode == 'objective' and it['t'].weight > 0]
        self.con_items = [it for it in self.items if it['t'].mode == 'constraint']
        self._ii = Q._IDX4
        self._acomp = [((2, 2, 0, 0), 4.0 / math.pi), ((1, 1, 1, 1), 2.0 / math.pi), ((0, 0, 2, 2), 4.0 / math.pi)]
        self.linear_rank = None

    # ---- 項の値 q(k) と ∂q/∂k ---------------------------------------------------------------
    def _target(self, it, xi_target=None):
        t = it['t'].target
        if it['kind'] == 'xi' and t is None:
            if xi_target is None:
                raise ValueError('xi の目標値 (spec の target か optimize の xi_target) がない')
            return float(np.asarray(xi_target)[it['idx']])
        return 0.0 if t is None else t

    def _eval(self, k, need_grad):
        """{名前: (q, dq/dk or None)}"""
        m = self.m
        out = {}
        kinds = {it['kind'] for it in self.items}
        v = h = Jc = None
        if 'sens' in kinds:
            e2 = self.err_rms ** 2 * ((k ** 2) @ m.errB)
        if 'h1' in kinds:
            v = m.h1_const + k @ m.H1
        if 'h4' in kinds or 'adts' in kinds:
            h = m.h4nf_fast(k)
            Jc = quad_jac(m._forms, k) if need_grad else None
        for it in self.items:
            i = it['idx']
            if it['kind'] == 'sens':
                q = math.sqrt(max(float(e2[i]), 1e-300))
                out[it['name']] = (q, self.err_rms ** 2 * k * m.errB[:, i] / q if need_grad else None)
            elif it['kind'] == 'h1':
                out[it['name']] = (complex(v[i]), m.H1[:, i] if need_grad else None)
            elif it['kind'] == 'h4':
                out[it['name']] = (complex(h[i]), Jc[i] if need_grad else None)
            elif it['kind'] == 'adts':
                q, c = self._acomp[i]
                out[it['name']] = (float(c * h[self._ii[q]].real + m.adts_offset[i]), c * Jc[self._ii[q]].real if need_grad else None)
            else:
                out[it['name']] = (float(self.xi0[i] + self.M[i] @ k), self.M[i] if need_grad else None)
        return out

    def value(self, k, xi_target=None):
        k = np.asarray(k, dtype=float)
        ev = self._eval(k, False)
        tot = 0.0
        for it in self.obj_items:
            q = ev[it['name']][0]
            r = it['s'] * (q - self._target(it, xi_target))
            tot += it['t'].weight * abs(r) ** 2
        return float(tot)

    def grad(self, k, xi_target=None):
        k = np.asarray(k, dtype=float)
        ev = self._eval(k, True)
        g = np.zeros(self.m.F)
        for it in self.obj_items:
            q, dq = ev[it['name']]
            r = q - self._target(it, xi_target)
            g += it['t'].weight * 2.0 * it['s'] ** 2 * np.real(np.conj(r) * dq)
        return g

    def table(self, k, xi_target=None):
        """Return the per-term breakdown (name, mode, weight, target, value, residual, scaled residual, contribution, den).
        項ごとの内訳: list of dict(name, mode, weight, target, value, residual |q−target|、scaled |scale·residual|、contribution = weight × scaled² (objective のみ)、den)。"""
        k = np.asarray(k, dtype=float)
        ev = self._eval(k, False)
        rows = []
        for it in self.items:
            q = ev[it['name']][0]
            t = self._target(it, xi_target)
            r = abs(q - t)
            sc = it['s'] * r
            rows.append(dict(name=it['name'], mode=it['t'].mode, weight=it['t'].weight, target=t, value=q, residual=float(r), scaled=float(sc),
                             contribution=float(it['t'].weight * sc * sc) if it['t'].mode == 'objective' else 0.0, den=it['den']))
        return rows

    # ---- 制約 --------------------------------------------------------------------------------
    def scipy_constraints(self, k0, xi_target=None):
        """Return the list of constraint dictionaries for scipy.optimize.minimize.
        scipy.optimize.minimize 用の制約の辞書の列。"""
        m = self.m
        k0 = np.asarray(k0, dtype=float)
        cons = []
        rows, rhs = [], []          # 線形の等式 A k = b
        nonlin = []
        for it in self.con_items:
            t = it['t']
            tgt = self._target(it, xi_target)
            if it['kind'] in ('h1', 'xi') and t.tol is None:
                if it['kind'] == 'h1':
                    c0 = m.h1_const[it['idx']]
                    a = m.H1[:, it['idx']]
                    tc = complex(tgt)
                    rows += [a.real, a.imag]
                    rhs += [tc.real - c0.real, tc.imag - c0.imag]
                else:
                    rows.append(self.M[it['idx']])
                    rhs.append(float(tgt) - self.xi0[it['idx']])
            else:
                nonlin.append((it, tgt))
        if rows:
            A = np.array(rows)
            b = np.array(rhs)
            nrm = np.linalg.norm(A, axis=1)
            keep = nrm > 1e-10 * max(nrm.max(), 1e-300)
            A, b, nrm = A[keep], b[keep], nrm[keep]
            A, b = A / nrm[:, None], b / nrm
            from scipy.linalg import qr
            _, Rr, piv = qr(A.T, pivoting=True)         # 従属な行を落とす
            d = np.abs(np.diag(Rr))
            rank = int(np.sum(d > 1e-9 * d[0])) if len(d) else 0
            sel = np.sort(piv[:rank])
            A, b = A[sel], b[sel]
            self.linear_rank = rank
            cons.append({'type': 'eq', 'fun': lambda k, A=A, b=b: A @ k - b, 'jac': lambda k, A=A: A})
        for it, tgt in nonlin:
            t = it['t']
            name = it['name']
            if t.tol is None:
                ref = abs(self._eval(k0, False)[name][0] - tgt)
                if not ref > 0:
                    ref = 1.0

                def f(k, name=name, tgt=tgt, ref=ref):
                    r = (self._eval(k, False)[name][0] - tgt) / ref
                    return np.array([np.real(r), np.imag(r)]) if np.iscomplexobj(r) or isinstance(r, complex) else np.array([float(r)])

                def jf(k, name=name, ref=ref):
                    q, dq = self._eval(k, True)[name]
                    dq = dq / ref
                    return np.vstack([np.real(dq), np.imag(dq)]) if np.iscomplexobj(dq) else np.array([dq])
                cons.append({'type': 'eq', 'fun': f, 'jac': jf})
            else:
                def f(k, name=name, tgt=tgt, tol=t.tol):
                    q = self._eval(k, False)[name][0]
                    return np.array([1.0 - abs(q - tgt) ** 2 / tol ** 2])

                def jf(k, name=name, tgt=tgt, tol=t.tol):
                    q, dq = self._eval(k, True)[name]
                    return np.array([-2.0 * np.real(np.conj(q - tgt) * dq) / tol ** 2])
                cons.append({'type': 'ineq', 'fun': f, 'jac': jf})
        return cons

    def max_constraint_violation(self, k, xi_target=None):
        """Return the maximum constraint violation over terms (|q - target|, or |q - target| - tol with a tolerance), in term units.
        制約の項ごとの |q − target| (tol つきは |q−target|−tol) の最大値 (項の単位)。"""
        ev = self._eval(np.asarray(k, dtype=float), False)
        v = 0.0
        for it in self.con_items:
            r = abs(ev[it['name']][0] - self._target(it, xi_target))
            v = max(v, r - (it['t'].tol or 0.0))
        return float(v)


def optimize_spec(obj, k0, xi_target=None, kmax=40.0, method='SLSQP', maxiter=400, ftol=1e-14, dk=None):
    """Optimize a `SxtSpecObjective` with scipy SLSQP (analytic gradient) under its constraint terms and |K2L| <= kmax.
    `SxtSpecObjective` の最適化 (scipy の SLSQP。勾配は解析)。制約は obj.spec の 'constraint' の項 (ξ は xi_target を既定の目標値に使う)、|K2L| ≤ kmax (スカラーか配列)。
    戻り (k, OptimizeResult)。k0 は制約を満たさなくてよい (SLSQP が寄せる)。"""
    from scipy.optimize import minimize
    k0 = np.asarray(k0, dtype=float)
    cons = obj.scipy_constraints(k0, xi_target)
    km = np.full(obj.m.F, float(kmax)) if np.isscalar(kmax) else np.asarray(kmax, dtype=float)
    lo, hi = -km, km
    if dk is not None:          # 出発点からの変化の上限 (信頼領域。スカラーかファミリーごとの配列、K2L の単位 [1/m²])
        d = np.full(obj.m.F, float(dk)) if np.isscalar(dk) else np.asarray(dk, dtype=float)
        lo, hi = np.maximum(lo, k0 - d), np.minimum(hi, k0 + d)
    bnds = list(zip(lo, hi))
    sc = max(obj.value(k0, xi_target), 1e-300)
    res = minimize(lambda k: obj.value(k, xi_target) / sc, k0, jac=lambda k: obj.grad(k, xi_target) / sc, method=method, bounds=bnds, constraints=cons,
                   options=dict(maxiter=maxiter, ftol=ftol))
    return res.x, res


def design_k2l(model, spec, k0=None, xi_target=None, kmax=40.0, J_ref=1e-6, delta_ref=0.01, maxiter=400, warn=True, err_rms=0.03, dk=None):
    """Compute the sextupole K2L [1/m^2] per family from a term specification `spec` with a single SLSQP optimization.
    項の指定 spec (SxtSpec・辞書・文字列・JSON。`SxtSpec.parse` が受ける形) を入力して、6 極の K2L [1/m²] (ファミリーごと) を計算して返す。
    項と重みの選択は呼ぶ側 (人) が決める (周波数マップ・RDT の一覧・ADTS を見て)。この関数は反復も探索もしない。1 回の最適化 (SLSQP) だけ。

    k0: 出発点 (既定は model.k_ref)。xi_target: ξ の目標 (None なら k0 の ξ)。spec に xi の項がなければ ξ = xi_target の等式制約を足す ('xi=ignore' で外せる)。kmax: |K2L| の上限 (スカラーか配列)。
    dk: 出発点 k0 からの各ファミリーの変化の上限 (信頼領域、K2L の単位 [1/m²]。スカラーか配列。None で制限なし)。err_rms: 'err_*' 項 (誤差感度) の K2L の相対 rms。
    戻り dict(k, success, message, nit, objective (指定の目的関数の値)、table (項ごとの内訳)、xi (最適化後の ξ)、constraint_violation (項の単位の最大値)、
    obj (SxtSpecObjective)、warnings)。"""
    k0 = model.k_ref.copy() if k0 is None else np.asarray(k0, dtype=float)
    if xi_target is None:
        xi_target = np.array(model.xi(k0))
    sp = SxtSpec.parse(spec)
    if not any(n in sp.terms for n in ('xi', 'xi_x', 'xi_y')):
        sp.terms['xi'] = TermSpec('constraint', 1.0, None)       # ξ を書かなければ xi_target の等式制約 (ξ = 0 の補正を壊さない)
    obj = SxtSpecObjective(model, sp, J_ref=J_ref, delta_ref=delta_ref, warn=warn, err_rms=err_rms)
    k, res = optimize_spec(obj, k0, xi_target, kmax, maxiter=maxiter, dk=dk)
    return dict(k=k, success=bool(res.success), message=str(res.message), nit=int(res.nit), objective=obj.value(k, xi_target),
                table=obj.table(k, xi_target), xi=model.xi(k), constraint_violation=obj.max_constraint_violation(k, xi_target), obj=obj,
                warnings=list(obj.warnings), k0=k0, xi_target=np.asarray(xi_target))
