# -*- coding: utf-8 -*-
"""v2.3.0 R5: 4 次の共鳴駆動項 h_jklm (j+k+l+m = 4、p = 0) の解析式 (6 極の 2 次を含む)。K64 にない新機能。導出と式は docs/49 §10。

定義 (リング先頭 s=0 の 1 周の写像): 角度・作用 ζ = √(2Jx) e^{iφx}、v = √(2Jy) e^{iφy} (正規化座標で ζ = x_n − i p_n)。ハミルトニアン H(s) の
相互作用表示 (線形の流れで動く座標) の 1 周の写像は、関数 g に対して g∘Φ = e^{:F:} g (:f: g = {f, g})、F = −h3 − h4 − … で、
  h3 = ∫ H3 ds                                                       (1 次、resonance_terms の h_jklm の 3 次の項)
  h4 = ∫ H4 ds + ½ Σ_{s_l > s_e} {H3(s_l), H3(s_e)} (ds_l ds_e)         (4 次: 4 次のハミルトニアンの 1 次 + 3 次の 2 次 (Magnus の第 2 項、順序つき))
ポアソン括弧は ζ・ζ̄ について {ζ, ζ̄} = 2i: {ζ^a ζ̄^b, ζ^c ζ̄^d} = 2i (a d − b c) ζ^{a+c−1} ζ̄^{b+d−1} (v も同様)。
h4 の係数が単項式 ζ^j ζ̄^k v^l v̄^m の係数 h_jklm (H の係数そのもの。3 次と同じ符号・規格化)。h_kjml = h̄_jklm。
4 次のハミルトニアンの 1 次: 8 極 (OCT、K3≠0 の MULT)、運動学的な項 (px²+py²)²/8、四極・K1 つき偏向磁石の端部の非線形フリンジ (ttfrin4)。
3 次の 2 次 (対の項): 6 極どうし (厚い要素の内部は閉形式)、偏向磁石の 3 次の項どうし、6 極 × 偏向磁石 (交差項)。偏向磁石の 3 次の項は bendnl.py が追跡から取り出す。
ADTS との関係: h22000・h11110・h00220 から ∂νx/∂Jx = 4 h22000/π、∂νx/∂Jy = 2 h11110/π、∂νy/∂Jy = 4 h00220/π (写像の対角部分)。
`adts()` の値は、これに「3 次の項の非共鳴の摂動」(正規形への変換。cot πx_k × |h3|² の項) を足したもの (`adts_from_rdt4`)。
"""
import math
import numpy as np
from . import resdrive as R
from . import bendnl

TWO_PI = 2.0 * math.pi


def _monos(N):
    return [(j, k, l, N - j - k - l) for j in range(N + 1) for k in range(N + 1 - j) for l in range(N + 1 - j - k)]


MONOS3 = _monos(3)          # 20
MONOS4 = _monos(4)          # 35
_IDX4 = {q: i for i, q in enumerate(MONOS4)}
_IDX3 = {q: i for i, q in enumerate(MONOS3)}


def mono_name4(q):
    return 'h%d%d%d%d0' % tuple(q[:4])


def _bracket_table():
    """{ζ^a..} の 3 次どうしの括弧: [(i3, i3', i4, 係数)] (i3, i3' は MONOS3 の添字、i4 は MONOS4 の添字、係数は 2i (a d − b c))"""
    tab = []
    for i, q in enumerate(MONOS3):
        for i2, q2 in enumerate(MONOS3):
            j, k, l, m = q
            j2, k2, l2, m2 = q2
            c = j * k2 - k * j2
            if c:
                t = (j + j2 - 1, k + k2 - 1, l + l2, m + m2)
                tab.append((i, i2, _IDX4[t], 2j * c))
            c = l * m2 - m * l2
            if c:
                t = (j + j2, k + k2, l + l2 - 1, m + m2 - 1)
                tab.append((i, i2, _IDX4[t], 2j * c))
    return tab


_BR = _bracket_table()
_BR_I = np.array([t[0] for t in _BR])
_BR_J = np.array([t[1] for t in _BR])
_BR_T = np.array([t[2] for t in _BR])
_BR_C = np.array([t[3] for t in _BR])


def omega2_from_units(A, IG):
    """順序つきの単位の列 (s の順、各単位 = 薄い kick または厚い要素) の 2 次の項 ½ Σ_{later > earlier} {H_l, H_e} の h4 の 35 成分 (複素)。
    A: (20, N) 各単位の ∫H3 の単項式の係数 (MONOS3 の順)、IG: (20, 20, N) 各単位の内部の ∬_{s1>s2} F_q(s1) F_q'(s2) (薄い kick は ½ A_q A_q')。
    P[q, q'] = Σ_{l>e} A_q(l) A_q'(e) + Σ_u IG[q, q', u]"""
    C = np.cumsum(A, axis=1) - A
    P = A @ C.T + IG.sum(axis=2)
    h = np.zeros(len(MONOS4), dtype=complex)
    np.add.at(h, _BR_T, 0.5 * _BR_C * P[_BR_I, _BR_J])
    return h


def _monos5(ms):
    return [tuple(q) + (0,) for q in ms]


def sext_units(nl3, method='analytic', nslice=8):
    """6 極 (K2≠0 の厚い・薄い要素) の単位: (A (20, Ne), IG (20, 20, Ne))。素子内の対の積分は resdrive.elem_pair_integrals (conj つき) から
    ∬ F_q(s1) F_q'(s2) = c_q c_q' · pi[q, q̄'] (q̄ = (k, j, m, l)) で作る (係数 c は実数)"""
    m5 = _monos5(MONOS3)
    Ne = len(nl3)
    A = np.zeros((20, Ne), dtype=complex)
    IG = np.zeros((20, 20, Ne), dtype=complex)
    conj_idx = [_IDX3[(q[1], q[0], q[3], q[2])] for q in MONOS3]
    for ie, e in enumerate(nl3):
        dens = R._dens_strengths(e)
        integ = R.elem_integrals(e, m5, method, nslice)
        pi = R.elem_pair_integrals(e, m5, method, nslice)
        cf = np.array([R.mono_coef(dens, q) for q in m5])
        A[:, ie] = cf * integ
        IG[:, :, ie] = np.outer(cf, cf) * pi[:, conj_idx]
    return A, IG


def bend_units(elements, rows, nslice=16, momentum=3e9, warnings=None, delta=0.0):
    """偏向磁石の 3 次の項の薄い kick の単位 (入口の端・nslice 個のスライス・出口の端): list[(素子の添字, a (20, nslice+2))] (a は H3 の係数、2^{3/2} を含まない)"""
    warnings = warnings if warnings is not None else []
    out = []
    skipped = {}
    m5 = _monos5(MONOS3)
    for i, el in enumerate(elements):
        if R._tn(el) != 'Bend':
            continue
        L = float(getattr(el, 'L', 0.0) or 0.0)
        ang = float(getattr(el, 'ANGLE', 0.0) or 0.0)
        if L <= 0 or (ang == 0.0 and float(getattr(el, 'K0', 0.0) or 0.0) == 0.0):
            continue
        if float(getattr(el, 'ROTATE', 0.0) or 0.0) or float(getattr(el, 'ROTATION', 0.0) or 0.0):
            skipped.setdefault('ROTATE・ROTATION≠0', []).append(el.name)
            continue
        try:
            a, info = bendnl.bend_atoms_auto(el, rows[i], rows[i + 1], m5, nslice=nslice, momentum=momentum, delta=delta)
        except Exception as ex:
            warnings.append('偏向磁石 %s: 3 次の項の取り出しに失敗 (%s: %s)' % (el.name, type(ex).__name__, ex))
            continue
        if info['w_err'] > 1e-6:
            warnings.append('偏向磁石 %s: スライスの線形写像で伝播した w が Twiss と %.1e 違う' % (el.name, info['w_err']))
        out.append((i, a / 2.0 ** 1.5))
    for why, names in skipped.items():
        warnings.append('偏向磁石の 3 次の項: %s の偏向磁石 %d 台 (%s ...) は入れていない' % (why, len(names), names[0]))
    return out


def assemble_units(elements, nl3, Asx, IGsx, bunits):
    """6 極の単位と偏向磁石の単位を s の順に並べる。戻り: A (20, N)、IG (20, 20, N)、kind (N,) (0: 6 極、1: 偏向磁石)"""
    sext_at = {e.idx: ie for ie, e in enumerate(nl3)}
    bend_at = {i: a for i, a in bunits}
    cols_A, cols_IG, kinds = [], [], []
    for i in range(len(elements)):
        if i in sext_at:
            ie = sext_at[i]
            cols_A.append(Asx[:, ie:ie + 1])
            cols_IG.append(IGsx[:, :, ie:ie + 1])
            kinds.append(np.zeros(1, dtype=int))
        if i in bend_at:
            a = bend_at[i]
            cols_A.append(a)
            cols_IG.append(0.5 * a[:, None, :] * a[None, :, :])
            kinds.append(np.ones(a.shape[1], dtype=int))
    if not cols_A:
        return np.zeros((20, 0), dtype=complex), np.zeros((20, 20, 0), dtype=complex), np.zeros(0, dtype=int)
    return np.concatenate(cols_A, axis=1), np.concatenate(cols_IG, axis=2), np.concatenate(kinds)


def octupole_h4(nl4, method='analytic', nslice=8):
    """8 極 (K3≠0) の 4 次のハミルトニアンの 1 次: h_jklm = Σ c_jklm ∫ w^j w̄^k w_y^l w̄_y^m ds。戻り (35 成分)"""
    h = np.zeros(35, dtype=complex)
    m5 = _monos5(MONOS4)
    for e in nl4:
        dens = R._dens_strengths(e)
        integ = R.elem_integrals(e, m5, method, nslice)
        cf = np.array([R.mono_coef(dens if e.L > 0 else e.strengths, q) for q in m5])
        h += cf * integ
    return h


_KIN_POLY = {(0, 4, 0, 0): 1.0 / 8, (0, 2, 0, 2): 2.0 / 8, (0, 0, 0, 4): 1.0 / 8}
_GL = np.polynomial.legendre.leggauss(8)


def _mat_key(el):
    return (R._tn(el), float(getattr(el, 'L', 0.0) or 0.0), float(getattr(el, 'K1', 0.0) or 0.0), float(getattr(el, 'ANGLE', 0.0) or 0.0),
            tuple(sorted((int(n), float(v)) for n, v in (getattr(el, 'KnL', {}) or {}).items())))


def _quad_ndiv(el, eps=0.2):
    """K64 の四極の本体 (tsolqu) の分割数: ndiv = 1 + int(L·|K1·L|/eps)、eps = 0.2。本体は分割ごとに (厳密な線形写像 + 運動学的な補正のドリフト) の組で、
    補正のドリフト (x += px (1/pz − 1) ds) は両端の半分 + 内部の節点で置く = ∫ p⁴/8 ds の台形則 (節点 ndiv+1 個、誤差 ~ ds²)"""
    L = float(el.L)
    return 1 + int(abs(L * (float(el.K1) * L) / eps))


def kinematic_h4(elements, rows, warnings=None, ngl=8, quad_rule='k64', eps=0.2, skip=()):
    """運動学的な項 H4 = (px² + py²)²/8 (厳密なドリフトの展開の 4 次) の 1 次: h_jklm = ∫ H4 ds (px = Re[w' ζ] など。素子内の w' は線形の転送行列)。
    積分は Gauss-Legendre ngl 点 (厳密)。quad_rule='k64' (既定) では K64 の四極の積分器の台形則 (節点 `_quad_ndiv(el)+1` 個) で足す
    (K64 の追跡と比べるとき。'exact' は厳密な積分で、四極の積分器の誤差 ~ ds² 分だけ追跡と違う)。
    四極・偏向磁石 (端部なし)・MULT (K1) は内部の転送行列、それ以外はドリフトとして扱う。戻り (35 成分)"""
    nodes, wts = np.polynomial.legendre.leggauss(ngl)
    nodes = 0.5 * (nodes + 1.0)
    wts = 0.5 * wts
    cache = {}
    h = np.zeros(35, dtype=complex)
    m5 = _monos5(MONOS4)
    for i, el in enumerate(elements):
        L = float(getattr(el, 'L', 0.0) or 0.0)
        if L <= 0 or i in skip:
            continue
        if quad_rule == 'k64' and R._tn(el) == 'Quadrupole' and float(el.K1) != 0.0:
            nd = _quad_ndiv(el, eps)
            tt = np.arange(nd + 1) / nd
            ww = np.full(nd + 1, 1.0 / nd)
            ww[0] = ww[-1] = 0.5 / nd
            key = ('q', L, float(el.K1), nd)
        else:
            key = _mat_key(el) + (ngl,)
            tt, ww = nodes, wts
        if key not in cache:
            Ms = []
            for t in tt:
                sl = R._slice_for_optics(el, L * float(t)) if t > 0 else None
                if sl is None:
                    s_ = L * float(t)
                    Ms.append(np.array([[1, s_, 0, 0], [0, 1, 0, 0], [0, 0, 1, s_], [0, 0, 0, 1]], dtype=float))
                else:
                    Ms.append(np.asarray(sl.transfer_matrix(0.0), dtype=float)[:4, :4])
            cache[key] = np.array(Ms)
        Ms = cache[key]
        w0 = bendnl.w_from_row(rows[i])
        wn = np.einsum('nab,b->an', Ms, w0)
        amp = bendnl.harmonic_amplitudes(_KIN_POLY, wn, m5)
        h += (amp * (ww * L)[None, :]).sum(axis=1)
    return h


def fringe_edges(elements, rows, warnings=None, skip=()):
    """四極・K1 つき偏向磁石の端部の非線形フリンジ (K64 の既定 ON) の一覧: [(σ, K1, row)]。resdrive.quad_fringe_adts と同じ選択"""
    g = lambda el, k: float(getattr(el, k, 0.0) or 0.0)
    out = []
    for i, el in enumerate(elements):
        if R._tn(el) not in ('Quadrupole', 'Bend') or i in skip:
            continue
        L = float(el.L)
        k1 = float(el.K1)
        if L <= 0 or k1 == 0.0 or g(el, 'DISFRIN') != 0.0:
            continue
        mf = int(round(g(el, 'FRINGE')))
        for sig, r, on in ((1.0, rows[i], mf > -4 and mf != 2), (-1.0, rows[i + 1], mf > -4 and mf != 1)):
            if on:
                out.append((sig, k1, r))
    return out


def fringe_h4(elements, rows, warnings=None, skip=()):
    """端部のフリンジ H_f = σ (K1/12)[(x³ + 3xy²) px − (y³ + 3x²y) py] δ(s − s端) の 1 次。戻り (35 成分)"""
    h = np.zeros(35, dtype=complex)
    for sig, k1, r in fringe_edges(elements, rows, warnings, skip):
        c = sig * k1 / 12.0
        f = {(3, 1, 0, 0): c, (1, 1, 2, 0): 3 * c, (0, 0, 3, 1): -c, (2, 0, 1, 1): -3 * c}
        h += bendnl.harmonic_amplitudes(f, bendnl.w_from_row(r), _monos5(MONOS4))[:, 0]
    return h



def bend_h4_extracted(elements, rows, bunits, momentum=3e9, warnings=None, delta=0.0):
    """偏向磁石の 4 次のハミルトニアンの 1 次の項 (運動学的な項・端 (E1/E2 の楔と硬い端)・K1 つきの本体と端のフリンジ・積分器の誤差を全部含む):
    要素ごとに、その要素 1 個の追跡の写像の 3 次のテイラー展開から取り出した h4 (`rdttrack.element_h4`) から、要素の内部の 3 次どうしの対の項
    (同じ薄い kick の列 = bend_units の原子による Ω2) を引いた残り。入口の Twiss の w で ζ の単項式にする。
    戻り (h (35 成分)、取り出せた偏向磁石の添字の集合)"""
    from . import rdttrack
    warnings = warnings if warnings is not None else []
    h = np.zeros(35, dtype=complex)
    done = set()
    m5 = _monos5(MONOS4)
    for i, a in bunits:
        el = elements[i]
        try:
            h4, info = rdttrack.element_h4(el, momentum, orbit=bendnl.row_orbit(rows[i]), delta=delta)
        except Exception as ex:
            warnings.append('偏向磁石 %s: 4 次の項の取り出しに失敗 (%s: %s)。運動学的な項だけ入れる' % (el.name, type(ex).__name__, ex))
            continue
        he = bendnl.harmonic_amplitudes(h4, bendnl.w_from_row(rows[i]), m5)[:, 0]
        within = omega2_from_units(a, 0.5 * a[:, None, :] * a[None, :, :])
        h += he - within
        done.add(i)
    return h, done


class RDT4Result:
    def __init__(self):
        self.warnings = []

    def table(self, names=None):
        lines = ['4 次の共鳴駆動項 h_jklm0 (写像の 4 次の生成関数の係数。8 極 1 次 + 運動学 + 端部 + 3 次の 2 次 (6 極・偏向磁石)。条件: %s、δ=0、リング先頭が位相の基準)' % self.method]
        lines.append('  作動点 νx=%.6f νy=%.6f' % (self.qx, self.qy))
        for nm in (names or self.names):
            v = self.total[nm]
            lines.append('  %-8s % .6e %+.6e i   |h| = %.6e' % (nm, v.real, v.imag, abs(v)))
        for w in self.warnings:
            lines.append('  (注意) ' + w)
        return '\n'.join(lines)

    __str__ = table

    def adts(self):
        """対角部分から ∂ν/∂J (写像の h だけ。正規形の補正を含まない): {'xx','xy','yy'}"""
        h = self.total_arr
        return {'xx': 4.0 * h[_IDX4[(2, 2, 0, 0)]].real / math.pi, 'xy': 2.0 * h[_IDX4[(1, 1, 1, 1)]].real / math.pi, 'yy': 4.0 * h[_IDX4[(0, 0, 2, 2)]].real / math.pi}


def rdt4(lat, terms=None, order2=True, octupole=True, kinematic=True, fringe=True, bend=True, bend4=True, method='analytic', nslice=8,
         bend_nslice=16, rows=None, momentum_eV=3e9):
    """4 次の共鳴駆動項 h_jklm0 (j+k+l+m=4)。戻り RDT4Result:
      .total[名前] (複素数)、.total_arr (MONOS4 の順 35 成分)、.parts (辞書: 'sext2' 6 極どうし、'bend2' 偏向磁石の 3 次どうし、'cross' 6 極 × 偏向磁石、
      'octupole'、'kinematic'、'fringe' 各 35 成分)、.names、.warnings。
    order2=True: 3 次 (6 極・偏向磁石) の 2 次の項 (Magnus の第 2 項) を入れる。False で 4 次のハミルトニアンの 1 次 (8 極・運動学・フリンジ) だけ。
    bend=False: 偏向磁石の 3 次の項を入れない (6 極どうしの 2 次だけ)。octupole・kinematic・fringe: それぞれの 1 次の項を入れるか。
    terms: 名前 ('h22000'…) のリスト (既定: 35 項全部)"""
    elements, rows = R._tables(lat, rows=rows)
    res = RDT4Result()
    R._check_uncoupled(rows, res.warnings)
    qx, qy = float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI
    nl = R.collect_nonlinear(elements, rows, 0.0, (0.0, 0.0), res.warnings)
    nl3 = [e for e in nl if 2 in e.strengths]
    nl4 = [e for e in nl if 3 in e.strengths]
    z35 = np.zeros(35, dtype=complex)
    parts = {k: z35.copy() for k in ('sext2', 'bend2', 'cross', 'octupole', 'kinematic', 'fringe', 'bend4')}
    if octupole and nl4:
        parts['octupole'] = octupole_h4(nl4, method, nslice)
    res.h3 = None
    bu = bend_units(elements, rows, bend_nslice, momentum_eV, res.warnings) if (bend and (order2 or bend4)) else []
    skip = set()
    if bend4 and bu:
        parts['bend4'], skip = bend_h4_extracted(elements, rows, bu, momentum_eV, res.warnings)
    if kinematic:
        parts['kinematic'] = kinematic_h4(elements, rows, res.warnings, skip=skip)
    if fringe:
        parts['fringe'] = fringe_h4(elements, rows, res.warnings, skip=skip)
    if order2:
        Asx, IGsx = sext_units(nl3, method, nslice)
        A, IG, kind = assemble_units(elements, nl3, Asx, IGsx, bu)
        if A.shape[1]:
            m0, m1 = kind == 0, kind == 1
            h_all = omega2_from_units(A, IG)
            h_s = omega2_from_units(A[:, m0], IG[:, :, m0]) if m0.any() else z35
            h_b = omega2_from_units(A[:, m1], IG[:, :, m1]) if m1.any() else z35
            parts['sext2'], parts['bend2'], parts['cross'] = h_s, h_b, h_all - h_s - h_b
            res.h3 = A.sum(axis=1)                      # 3 次の 1 周の和 (MONOS3 の順)
    res.parts = parts
    res.total_arr = sum(parts.values())
    res.names = [mono_name4(q) for q in MONOS4]
    res.total = {nm: complex(res.total_arr[i]) for i, nm in enumerate(res.names)}
    res.method, res.qx, res.qy = method, qx, qy
    res.order2 = order2
    return res


def nf_correction(h3, qx, qy, with_i=True):
    """3 次の項の 1 周の和 h3 (MONOS3 の順 20 成分、複素) から、写像の対角部分 (4 h22000/π など) に足す正規形の補正 (3 次の非共鳴の摂動、docs/49 §10.4):
      S_k = (1/(4π sin πx_k)) ∬ F_q(s1) F̄_q'(s2) [e^{−iπx_k} (s1>s2) + e^{+iπx_k} (s1<s2)] = −i P_>/(2π) + e^{iπx_k} T/(4π sin πx_k)
    (P_> = Σ_{s1>s2} F_q F̄_q'、T = (ΣF_q)(ΣF̄_q') = 全周の和の積)。−i P_>/(2π) の部分は Magnus の第 2 項 (写像の対角部分) に等しく、
    残りの (cot πx_k + i) T/(4π) が補正。K = −½ Σ_k Σ_(q,q') (k·∇_J)(J の単項式) · (cot πx_k + i) T_qq'/(4π)、∂ν/∂J = ∂²K/∂J²。戻り {'xx','xy','yy'} [1/(m rad)]"""
    byk = {}
    for q in MONOS3:
        byk.setdefault((q[0] - q[1], q[2] - q[3]), []).append(q)
    Kc = {}
    for kv, qs in byk.items():
        xk = kv[0] * qx + kv[1] * qy
        sn = math.sin(math.pi * xk)
        if abs(sn) < 1e-9:
            continue
        kern = (math.cos(math.pi * xk) / sn + (1j if with_i else 0.0)) / (4.0 * math.pi)
        for q in qs:
            for q2 in qs:
                T = 8.0 * h3[_IDX3[q]] * np.conj(h3[_IDX3[q2]])
                px, py = (q[0] + q[1] + q2[0] + q2[1]) / 2.0, (q[2] + q[3] + q2[2] + q2[3]) / 2.0
                for (ax, ay), co in (((px - 1, py), kv[0] * px), ((px, py - 1), kv[1] * py)):
                    if co == 0:
                        continue
                    key = (int(round(ax)), int(round(ay)))
                    Kc[key] = Kc.get(key, 0j) - 0.5 * co * kern * T
    xx, xy, yy = R._deriv_from_K({k: v.real for k, v in Kc.items()})
    return {'xx': xx, 'xy': xy, 'yy': yy}


def fourth_order_adts(elements, rows, bunits, momentum=3e9, delta=0.0, quad_rule='k64', warnings=None):
    """v2.3.0 R6: `resdrive.adts` の 4 次のハミルトニアンの 1 次の部分 (運動学的な項・端部のフリンジ・偏向磁石の 4 次の項) の ADTS。対角 h22000・h11110・h00220 から
    ∂νx/∂Jx = 4 h22000/π、∂νx/∂Jy = 2 h11110/π、∂νy/∂Jy = 4 h00220/π。
      kin    : 運動学的な項 (px²+py²)²/8。四極は K64 の積分器の台形則 (quad_rule='k64'、'exact' で厳密な積分)。偏向磁石以外の素子
      fringe : 四極の端部のフリンジ (偏向磁石以外)
      bend4  : 偏向磁石 (K0≠0 を含む) の 4 次の項の全部 (運動学・端・K1 の端のフリンジ・積分器の誤差。要素 1 個の追跡の写像の 3 次のテイラー展開から。δ のまわりでも取り出す)
    δ≠0: kin は /(1+δ)³ (H = −√((1+δ)² − p²) の展開)、fringe は δ に依らない (追跡の kick は p0 で規格化した場の強さ。検証: ASTRID2 の ∂(ADTS)/∂δ の追跡との比が 1/(1+δ) を掛けると 0.96→0.93)、
    bend4 は δ の追跡から取り出した値そのもの。戻り dict {'kin', 'fringe', 'bend4'} (各 {'xx','xy','yy'})"""
    warnings = warnings if warnings is not None else []
    h4b, skip = bend_h4_extracted(elements, rows, bunits, momentum, warnings, delta=delta)
    hk = kinematic_h4(elements, rows, warnings, quad_rule=quad_rule, skip=skip)
    hf = fringe_h4(elements, rows, warnings, skip=skip)
    if delta != 0.0:
        hk = hk / (1.0 + delta) ** 3          # 端部のフリンジは場の強さそのもの (追跡で δ に依らない) なので、そのまま

    def adts_of(h):
        return {'xx': 4.0 * float(h[_IDX4[(2, 2, 0, 0)]].real) / math.pi, 'xy': 2.0 * float(h[_IDX4[(1, 1, 1, 1)]].real) / math.pi,
                'yy': 4.0 * float(h[_IDX4[(0, 0, 2, 2)]].real) / math.pi}
    return {'kin': adts_of(hk), 'fringe': adts_of(hf), 'bend4': adts_of(h4b)}


def nf_correction_h4(h3, qx, qy, eps=1e-9):
    """v2.3.0 R6 (非対角の調和への拡張): 3 次の項 h3 (1 周の和、MONOS3 の順 20 成分) の非共鳴の摂動が 4 次の項 (MONOS4 の順 35 成分、複素) に足す補正。
    1 周の写像 M = R e^{:F:}、F = F3 + F4 (F3 = −h3 の符号の規約は rdt2nd の docstring と同じ、F4 = `rdt4().total_arr`) で、F3 の各調和 q (位相 φ_q = 2π[(j−k)νx + (l−m)νy]) を
    正準変換 e^{:W:}、W_q = h_q/(1 − e^{+iφ_q}) で消す (共鳴 φ_q ∈ 2πZ の項は消せず、ここでは 0 とする = 共鳴の近くでは使えない)と、BCH から 4 次の項は
        F4' = F4 + ½ {R⁻¹W, W}    (R⁻¹W の各項は W_q × e^{+iφ_q}、括弧は {ζ^a ζ̄^b, ζ^c ζ̄^d} = 2i(ad − bc) ζ^{a+c−1} ζ̄^{b+d−1}、v も同様)
    と変わる。対角 (k=0: h22000・h11110・h00220) の実部の 4/π・2/π・4/π 倍は `nf_correction` と 1e-14 で一致する (試験)。非対角の調和 (h31000・h20110 など) は、
    共鳴線 (k·ν = 整数) の近くで効く 4 次の駆動項の正規形の値 F4'。戻り: 35 成分の補正 (F4' − F4)"""
    ph = np.array([TWO_PI * ((q[0] - q[1]) * qx + (q[2] - q[3]) * qy) for q in MONOS3])
    den = 1.0 - np.exp(1j * ph)
    ok = np.abs(den) > eps
    W = np.where(ok, h3 / np.where(ok, den, 1.0), 0.0)
    RW = W * np.exp(1j * ph)
    out = np.zeros(len(MONOS4), dtype=complex)
    np.add.at(out, _BR_T, 0.5 * _BR_C * RW[_BR_I] * W[_BR_J])
    return out


def rdt4_normal_form(lat, rows=None, **kw):
    """4 次の項 (`rdt4`) に非共鳴の 3 次の摂動 (`nf_correction_h4`) を足した正規形の値 F4' (35 成分の dict と RDT4Result)。h3 は order2=True で得る"""
    kw.setdefault('order2', True)
    r = rdt4(lat, rows=rows, **kw)
    if r.h3 is None:
        raise ValueError('rdt4_normal_form: order2=True が要る (3 次の 1 周の和 h3 が無い)')
    c = nf_correction_h4(r.h3, r.qx, r.qy)
    tot = r.total_arr + c
    return {nm: complex(tot[i]) for i, nm in enumerate(r.names)}, r


def adts_from_rdt4(lat, rows=None, **kw):
    """4 次の写像の h (対角部分) + 3 次の項の正規形の補正 (`nf_correction`) から ∂ν/∂J を作る。`resdrive.adts` の 6 極 2 次 + 偏向磁石の 3 次 + 交差項に対応
    (8 極・運動学・フリンジ・偏向磁石の 4 次の項の 1 次の分も写像の h に入っているので、adts の .total と比べる)。戻り ({'xx','xy','yy'}, RDT4Result)"""
    r = rdt4(lat, rows=rows, **kw)
    a = r.adts()
    c = nf_correction(r.h3, r.qx, r.qy)
    return {k: a[k] + c[k] for k in a}, r
