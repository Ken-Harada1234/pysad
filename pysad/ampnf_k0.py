# -*- coding: utf-8 -*-
"""v2.4.0 G: K0≠0 の偏向磁石 (挿入光源・逆曲げ・MPW) の 1 周の写像のべき級数 (ampnf の TPSA) への取り込み。

【問題】 K0≠0 の偏向磁石の中の基準軌道 x0(s) は 0 でなく (K0 の双極子の軌道)、磁石の出口は原点でない (x ≠ 0, px ≠ 0)。K0=0 の磁石は原点のまわりの
  ハミルトニアンの項 (曲率 3 次・運動学・端の生成関数) を Lie 級数の流れで入れていたが、原点のまわりの展開は軌道が原点から離れた磁石の内部では使えない
  (F1 の試み: 磁石の出口が原点でなく、写像が定数項を持つので、原点まわりの bendnl.extract_map は成り立たなかった)。
【方法】 K64 の偏向磁石の追跡 (track6n.tbendcore1: K1=0 の ANGLE≠0 の本体、tsteer: ANGLE=0) は閉じた式 (扇形の厳密な写像 + 端の楔 E1・E2 + 硬い端の Maxwell の項) である。
  その式を、そのまま **べき級数 (TPSA) の演算で評価する** (ここの `TP` は定数部分を持つ級数で、sqrt・1/x・asin は定数部分 a0 のまわりの厳密な級数)。
  入力の級数 S = (x, px, y, py) は閉軌道を定数項に持つ (磁石の入口の軌道)。出力の級数は、その軌道のまわりで厳密に (次数 D まで) 展開した写像で、定数項が磁石の出口の軌道。
  つまり磁石の内部の軌道 x0(s) のまわりの展開は、閉じた式を軌道のまわりに展開することそのもの (スライスに分けない・Lie 級数の打ち切りも、数値の積分も無い。
  近似は次数 D の切り捨てだけ)。**この写像は曲率・運動学 (全次数)・端の楔・硬い端の項をすべて含み**、K0=0 の模型の bend_cubic・kin の項の選択 (オプション) は K0≠0 の磁石には効かない。
  運動量は、ampnf の正準運動量 (p0 で規格化) と追跡の内部の px/(1+δ) の変換を入口・出口で行う。δ は 5 番目の変数で、δ の全次数の依存が入る (ξ・J δ の項)。
  ROTATE した K0≠0 の磁石は ampnf の `_expand_rolls` が R⁻¹ ∘ (傾けない磁石) ∘ R に展開する (全次数の θ。垂直の軌道・結合が入る)。
【検証】 `check_against_tracker` (追跡の写像と全点で比較)。試験は pysad/tests/test_v24_G_k0_ampnf.py、記録は docs/history/v24_G.md、docs/25 §14.5。
【G の続き】 K1≠0 の BEND は tbendi の積分器そのものを TPSA で評価 (_tbendi)、F1・FB1・FB2 (tbendcore1・k_steer のフリンジの項)、DX・DY・DROTATE・CHI2 (tbshift・tbrot の TPSA 版。ANGLE=0 は ampnf が shift・rot に展開)、
  ANGLE + K0 = 0 は phib → 0 の極限。bend_cubic・kin・fringe が全部 False なら入口の軌道のまわりの 1 次 (linearize)。
【入れていない効果】 ANGLE=0 で K1≠0 かつ K0≠0 (K64 の tmulti)、K1≠0 かつ K0≠0 の F1・FB1・FB2、ANGLE=0 の CHI2、放射。
  端の hard-edge 構造の δ 依存 (K64 の式にある 1/(1+δ)) も、閉じた式の通りに入る (K0=0 の模型の δ に依らない端の生成関数より完全)。
  K64 の式の独立な基準: example/14_nonlinear_dynamics/normal_form/k0_independent_flow.py (実験室の直交座標の運動方程式の数値積分)。
"""
import math
import numpy as np


class TP:
    """定数部分を持つべき級数 (ampnf.get_tab の 5 変数・次数 D)。演算は + − × ÷ と sqrt・asin"""
    __slots__ = ('tab', 'a')

    def __init__(self, tab, a):
        self.tab, self.a = tab, a

    @staticmethod
    def const(tab, c):
        a = np.zeros(tab.N, dtype=complex)
        a[0] = c
        return TP(tab, a)

    def _w(self, o):
        return o if isinstance(o, TP) else TP.const(self.tab, o)

    def __add__(self, o):
        if isinstance(o, TP):
            return TP(self.tab, self.a + o.a)
        a = self.a.copy()
        a[0] += o
        return TP(self.tab, a)

    __radd__ = __add__

    def __neg__(self):
        return TP(self.tab, -self.a)

    def __sub__(self, o):
        return self + (-o if not isinstance(o, TP) else -o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, TP):
            return TP(self.tab, self.tab.mul(self.a, o.a))
        return TP(self.tab, self.a * o)

    __rmul__ = __mul__

    def inv(self):
        return _series(self, _ser_inv)

    def __truediv__(self, o):
        if isinstance(o, TP):
            return self * o.inv()
        return TP(self.tab, self.a / o)

    def __rtruediv__(self, o):
        return self.inv() * o

    def sqrt(self):
        return _series(self, _ser_sqrt)

    def asin(self):
        return _series(self, _ser_asin)

    @property
    def c0(self):
        return self.a[0]


def _ser_inv(a0, K):
    return np.array([(-1.0) ** k / a0 ** (k + 1) for k in range(K + 1)], dtype=complex)


def _binom(al, k):
    r = 1.0
    for j in range(k):
        r *= (al - j) / (j + 1.0)
    return r


def _ser_sqrt(a0, K):
    return np.array([_binom(0.5, k) * a0 ** (0.5 - k) for k in range(K + 1)], dtype=complex)


def _pmul(p, q, K):
    out = np.zeros(K + 1, dtype=complex)
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if i + j <= K:
                out[i + j] += pi * qj
    return out


def _ser_asin(a0, K):
    """asin(a0 + t) の t の級数。asin' = (1 − x²)^(−1/2) = q(t)^(−1/2)、q = (1−a0²) − 2 a0 t − t²"""
    q0 = 1.0 - a0 * a0
    r = np.zeros(K + 1, dtype=complex)
    if K >= 1:
        r[1] = -2.0 * a0 / q0
    if K >= 2:
        r[2] = -1.0 / q0
    g = np.zeros(K + 1, dtype=complex)
    rk = np.zeros(K + 1, dtype=complex)
    rk[0] = 1.0
    for k in range(K + 1):
        g += _binom(-0.5, k) * rk
        rk = _pmul(rk, r, K)
    g *= q0 ** -0.5
    c = np.zeros(K + 1, dtype=complex)
    c[0] = np.arcsin(a0)
    for k in range(K):
        c[k + 1] = g[k] / (k + 1)
    return c


def _series(T, fser):
    """f(T) = Σ c_k (T − T0)^k (T0 = T の定数部分)。fser(a0, K) が c_k"""
    tab = T.tab
    K = tab.D
    c = fser(T.a[0], K)
    eps = T.a.copy()
    eps[0] = 0.0
    out = np.zeros(tab.N, dtype=complex)
    out[0] = c[0]
    pw = np.zeros(tab.N, dtype=complex)
    pw[0] = 1.0
    for k in range(1, K + 1):
        pw = tab.mul(pw, eps)
        if not np.any(pw):
            break
        out += c[k] * pw
    return TP(tab, out)


# ---------------------------------------------------------------------------------------------
# 素子の写像 (K64 の追跡と同じ式。内部の運動量は px/(1+δ))
# ---------------------------------------------------------------------------------------------
def _tp_sqrt1m(tab, *sq):
    """sqrt(1 − Σ s²) の TP (sq は TP のリスト)"""
    acc = None
    for s in sq:
        t = s * s
        acc = t if acc is None else acc + t
    return (1.0 - acc).sqrt()


def _g(el, k):
    return float(getattr(el, k, 0.0) or 0.0)


def bend_k0_params(el):
    """K0≠0 の BEND (L>0) の写像のパラメータ (track6n.Prog._compile と同じ導出)。戻り dict。未対応の条件は NotImplementedError"""
    from . import lattice as _lat
    L, phi0, k0 = _g(el, 'L'), _g(el, 'ANGLE'), _g(el, 'K0')
    if L <= 0.0:
        raise NotImplementedError('ampnf_k0: L=0 の BEND は薄い kick (素子の線形の写像) で扱う: %s' % el.name)
    fb1, fb2, mfring, fringe = _lat._bend_fringe_params(el)
    if mfring <= 0:
        if mfring != -1:
            fb1 = 0.0
        if mfring != -2:
            fb2 = 0.0
    return dict(L=L, phi0=phi0, k0=k0, k1=_g(el, 'K1'), psi1=_g(el, 'E1'), psi2=_g(el, 'E2'), fb1=fb1, fb2=fb2, mfring=mfring, fringe=bool(fringe),
                dx=_g(el, 'DX'), dy=_g(el, 'DY'), theta=_lat._bend_theta(el), dtheta=_g(el, 'DROTATE'), dchi2=_g(el, 'CHI2'), eps0=_g(el, 'EPS'), name=el.name)


def bend_k0_map(el, v2=None, linear=False):
    """K0≠0 の BEND (L>0) の写像 f(S, tab, sd) → S'。S = [x, px, y, py] (TPSA の係数配列、px は正準 (p0 で規格化))。
    ANGLE≠0・K1=0 は K64 の tbendcore1 (扇形の厳密な写像 + 楔 E1・E2 + 硬い端の Maxwell の項 + F1 のフリンジ)、ANGLE≠0・K1≠0 は tbendi (K64 の数値積分器そのものを TPSA で評価。
    スライス数 ndiv = 1 + |ANGLE|/1e-3)、ANGLE=0 は k_steer (直線の磁石の曲がる軌道 + 楔 + Maxwell + F1)。DX・DY・ROTATE・DROTATE・CHI2 は tbshift・tbrot の TPSA 版
    (ANGLE=0 の磁石は ampnf._expand_rolls が ('shift')・('rot') に展開する)。ANGLE + K0 = 0 (場が 0) は K64 の式が phib = 0 で割るので、phib → 0 の極限 (phib = 1e-9 ANGLE)。
    linear=True: 入口の軌道のまわりの 1 次 (定数項 + ヤコビ行列。δ の依存は全次数) だけ (bend_cubic・kin・fringe を全部切った模型。`linearize`)"""
    from . import compat
    if v2 is None:
        v2 = compat.v2()
    P = bend_k0_params(el)
    phi0, k0 = P['phi0'], P['k0']
    if phi0 != 0.0:
        if P['k1'] != 0.0:
            core = _bendi_core(P, v2)
        else:
            core = _bend0_core(P)
        f = _wrap_shift(core, P)
    else:
        if k0 == 0.0:
            raise ValueError('ampnf_k0: ANGLE = 0 かつ K0 = 0 の BEND はドリフト: %s' % el.name)
        # ANGLE=0 の BEND の CHI2 は K64 の tsteer・tmulti が使わない (追跡は無視する)。ampnf も追跡に合わせて無視する (ampnf.build_ops が警告を出す)
        if P['k1'] != 0.0:
            if P['fb1'] or P['fb2']:
                raise NotImplementedError('ampnf: ANGLE=0 で K1≠0 かつ K0≠0 の BEND の F1・FB1・FB2 は未対応 (K64 の tmulti は gs が SIGSEGV で基準がない): %s' % el.name)
            core = _straight_k1_core(P)

            def f(S, tab, sd):
                return core(S, tab, sd)
            if linear:
                f = linearize(f)
            return f
        const = dict(al=P['L'], phib=-k0, tanp1=math.tan(P['psi1']), tanp2=math.tan(P['psi2']) if v2 else 0.0, fb1=P['fb1'], fb2=P['fb2'], fringe=P['fringe'],
                     psi1=P['psi1'], psi2=P['psi2'], ex1=bool(v2 and P['psi1'] != 0.0), ex2=bool(v2 and P['psi2'] != 0.0))
        # tanp2 は fb2≠0 のとき v2 でなくても入る (k_steer)
        if P['fb2'] != 0.0:
            const['tanp2'] = math.tan(P['psi2'])

        def f(S, tab, sd):
            return _tsteer(S, tab, sd, **const)
    if linear:
        f = linearize(f)
    return f


def linearize(f, s_loc=1e-3):
    """写像 f (S → S') を、入口の軌道 c0 (S の定数項) のまわりの 1 次 (S' = F0(δ) + J(δ)·(S − c0)。F0・J は δ の全次数) にする。
    bend_cubic・kin・fringe を全部切った模型 (磁石自身の非線形を落とす。軌道・線形の光学は残す)。F0・J は局所の級数 (c0 + s_loc u) から取る"""
    def g(S, tab, sd):
        c = [float(np.real(a[0])) for a in S]
        loc = []
        for j in range(4):
            a = tab.var(j, s_loc)
            a = a.copy()
            a[0] += c[j]
            loc.append(a)
        out = f(loc, tab, sd)
        z4 = tab.mons[:, :4].sum(axis=1)
        res = []
        for j in range(4):
            oj = np.asarray(out[j])
            F0 = np.where(z4 == 0, oj, 0.0)
            acc = F0.astype(complex)
            for k in range(4):
                Jk = np.zeros(tab.N, dtype=complex)
                for m in np.nonzero((z4 == 1) & (tab.mons[:, k] == 1))[0]:
                    e = list(tab.mons[m])
                    e[k] = 0
                    Jk[tab.idx[tuple(e)]] = oj[m] / s_loc
                if np.any(Jk):
                    dev = np.asarray(S[k]).copy()
                    dev[0] -= c[k]
                    acc = acc + tab.mul(Jk, dev)
            res.append(acc)
        return res
    return g


def _bend0_core(P):
    phi0, k0 = P['phi0'], P['k0']
    L = P['L']
    phib = phi0 + k0
    if phib == 0.0:
        phib = 1e-9 * phi0
    psi1, psi2 = P['psi1'], P['psi2']
    w = phi0 - psi1 - psi2
    al = L
    fb1, fb2 = P['fb1'], P['fb2']
    if fb1 != 0.0 or fb2 != 0.0:
        al = al - ((phi0 * fb1) ** 2 + (phi0 * fb2) ** 2) / al / 48.0 * math.sin(0.5 * w) / math.sin(0.5 * phi0)
    const = dict(al=al, phi0=phi0, ak0=phib - phi0, phib=phib, psi1=psi1, psi2=psi2, cosw=math.cos(w), sinw=math.sin(w), sinwp1=math.sin(phi0 - psi2),
                 fb1=fb1 if (P['mfring'] > 0 or P['mfring'] == -1) else 0.0, fb2=fb2 if (P['mfring'] > 0 or P['mfring'] == -2) else 0.0, fringe=P['fringe'])
    const['sqwh'] = const['sinw'] ** 2 / (1.0 + const['cosw']) if const['cosw'] >= 0.0 else 1.0 - const['cosw']

    def f(S, tab, sd):
        return _tbendcore(S, tab, sd, **const)
    return f


def _tbendcore(S, tab, sd, al, phi0, ak0, phib, psi1, psi2, cosw, sinw, sinwp1, sqwh, fb1=0.0, fb2=0.0, fringe=True):
    """track6n.tbendcore1 の TPSA 版 (F1 のフリンジ fb1・fb2 つき。z は追跡しない)。S の px・py は正準運動量"""
    dxfr1 = fb1 * fb1 / (al / phib) / 24.0 if fb1 != 0.0 else 0.0
    dyfr1 = fb1 / (al / phib) ** 2 / 6.0 if fb1 != 0.0 else 0.0
    dyfra1 = 4.0 * dyfr1 / (fb1 * fb1) if (fb1 != 0.0 and fringe) else 0.0
    dxfr2 = fb2 * fb2 / (al / phib) / 24.0 if fb2 != 0.0 else 0.0
    dyfr2 = fb2 / (al / phib) ** 2 / 6.0 if fb2 != 0.0 else 0.0
    dyfra2 = 4.0 * dyfr2 / (fb2 * fb2) if (fb2 != 0.0 and fringe) else 0.0
    cosp1, sinp1, cosp2, sinp2 = math.cos(psi1), math.sin(psi1), math.cos(psi2), math.sin(psi2)
    dcosp = cosp1 - cosp2
    sinp2wp1 = sinp2 + sinwp1
    sinp2p1 = sinp2 + sinp1
    cosp1p2 = cosp1 * cosp2
    rhob = al / phib
    rho0 = al / phi0
    drhob = -rho0 * ak0 / phib
    x, pxc, y, pyc = (TP(tab, a) for a in S)
    dp = TP(tab, tab.var(4, sd))
    pp = dp + 1.0
    pinv = pp.inv()
    rhoe = pp * rhob
    pxi = pxc * pinv
    pyi = pyc * pinv
    pzi = _tp_sqrt1m(tab, pxi, pyi)
    dpzi = pzi - 1.0
    sp1 = sinp1 / pzi
    x1 = x / (1.0 - pxi * sp1 + (cosp1 - 1.0))
    dz1 = x1 * sp1
    y1 = y + pyi * dz1
    px1 = pxi * cosp1 + pzi * sinp1
    if fb1 != 0.0:
        x1 = x1 + dp * pinv * dxfr1
        py1 = pyi + (dyfr1 - y1 * y1 * dyfra1) * y1 * pinv * pinv
    else:
        py1 = pyi
    pv1sqi = (1.0 - px1 * px1).inv()
    fa = y1 / rhoe * pv1sqi.sqrt()
    q = y1 * (1.0 / rhob)
    f = (1.0 - q * q * (1.0 / 6.0)) * fa
    ff = (f + fa) * y1 * pv1sqi * 0.25
    x2 = x1 + ff
    py2 = py1 - px1 * f
    ph2 = (1.0 - py2 * py2).sqrt()
    dph2 = ph2 - 1.0
    dpx2 = pxi * cosp1 + (dpzi - dph2) * sinp1
    pz2 = _tp_sqrt1m(tab, px1, py2)
    drho1 = rhob * dp + drhob
    drho = drho1 + rhoe * dph2
    t2 = (px1 + ph2 * sinp1) / (pz2 + ph2 * cosp1)
    dpx3a = (x2 * sinw - drho * sinp2wp1) / rhoe
    dpx3 = dpx3a - dpx2 * (cosw - t2 * sinw)
    px3 = ph2 * sinp2 + dpx3
    sqwht2 = t2 * sinw + sqwh
    px1px3 = ph2 * sinp2p1 + dpx3a + dpx2 * sqwht2
    pz3 = _tp_sqrt1m(tab, px3, py2)
    dpz3 = pz3 - 1.0
    dpz32 = px1px3 * (px1 - px3) / (pz2 + pz3)
    pz3ph2 = pz3 + ph2 * cosp2
    pz2px1 = pz2 * cosp1 + px1 * sinp1
    t3 = (px3 + ph2 * sinp2) / pz3ph2
    t2t3 = (ph2 * sinp2 + px1px3) / pz3ph2 + (ph2 * sinp1 + px1 * (dpz32 - ph2 * dcosp) / pz3ph2) / (pz2 + ph2 * cosp1)
    t4 = (t3 * sinp2 + cosp2) * pz2px1
    x3 = (x2 * (1.0 - rho0 * t3 * sinw / (rhoe * ph2) + (cosw - 1.0))
          + (rho0 * ((t2t3 * cosw + (1.0 - t2 * t3) * sinw)) * dpx2
             - drho * (-(rho0 * t3 / rhoe) * sinp2wp1 - dpz32 - pz2 * sqwh - px1 * sinw)) / ph2)
    arg = (dpx2 * ((t2 * (pz3 * cosp2 + px3 * sinp2) - px1 * cosp2) * sinp1
                   - (t3 * pz2px1 - px3 * cosp1) * sinp2
                   + dpz32 * cosp1p2 + t4 * sqwht2) + dpx3a * t4) / (ph2 * ph2)
    da = arg.asin()
    phi0a = da + phi0
    y3 = y1 + py2 * rhoe * phi0a
    pv2sqi = (1.0 - px3 * px3).inv()
    fa = y3 / rhoe * pv2sqi.sqrt()
    q = y3 * (1.0 / rhob)
    f = (1.0 - q * q * (1.0 / 6.0)) * fa
    ff = (f + fa) * y3 * pv2sqi * 0.25
    py4 = py2 - px3 * f
    x4 = x3 - ff
    if fb2 != 0.0:
        x4 = x4 - dp * pinv * dxfr2
        py4 = py4 + (dyfr2 - y3 * y3 * dyfra2) * y3 * pinv * pinv
    pz4 = _tp_sqrt1m(tab, px3, py4)
    dpz4 = pz4 - 1.0
    pxo = -cosp2 * dpx3 + (dpz4 - dpz3 - dpx3 * t3) * sinp2
    dz4 = x4 * sinp2 / (pz4 * cosp2 + px3 * sinp2)
    xo = x4 * cosp2 + pxo * dz4
    yo = y3 + py4 * dz4
    return [xo.a, (pxo * pp).a, yo.a, (py4 * pp).a]


def _tsteer(S, tab, sd, al, phib, tanp1, tanp2, fb1=0.0, fb2=0.0, fringe=True, psi1=0.0, psi2=0.0, ex1=False, ex2=False):
    """track6n.k_steer の TPSA 版 (ndiv = 1、放射なし)。直線の座標系の磁石 (ANGLE=0、K0≠0)。phib = −K0 (k_steer の符号)。
    ex1・ex2: v2 (B221) の端の面は CALC (tsteee) と同じ厳密な硬い端の鎖 (drift・Maxwell・body。F1・FB があるときは K64 の薄い蹴り)"""
    rhob = al / phib
    x, pxc, y, pyc = (TP(tab, a) for a in S)
    dp = TP(tab, tab.var(4, sd))
    pp = dp + 1.0
    pinv = pp.inv()
    rhoe = pp * rhob
    irhoe = rhoe.inv()
    if fb1 != 0.0:
        dxfr1 = -fb1 ** 2 / rhob / 24.0
        dyfr1 = fb1 / rhob ** 2 / 6.0
        dyfra1 = 4.0 * dyfr1 / fb1 ** 2 if fringe else 0.0
    if ex1:
        # B221: 鎖は最初の断片の tblfre (F1 の dxfr1・dyfr1 の項) の後に通す (CALC の tsteee と同じ順序)
        if fb1 != 0.0:
            x = x + dp * pinv * dxfr1
            pyc = pyc + (pinv * (dyfr1 - y * y * dyfra1)) * y
        x, pxc, y, pyc = _edge_exact(x, pxc, y, pyc, pp, -rhob, psi1, True)
        tanp1 = 0.0
        pxi = pxc * pinv
        py0 = pyc * pinv
    else:
        pxi = pxc * pinv - x * irhoe * tanp1
        py0 = pyc * pinv + y * irhoe * tanp1
        if fb1 != 0.0:
            x = x + dp * pinv * dxfr1
            py0 = py0 + (pinv * pinv * (dyfr1 - y * y * dyfra1)) * y
    if ex1:
        pyi = py0
        xi = x
    else:
        pv1 = _tp_sqrt1m(tab, pxi)
        fa = y * irhoe / pv1
        q = y * (1.0 / rhob)
        f = (1.0 - q * q * (1.0 / 6.0)) * fa
        ff = (0.5 - q * q * (1.0 / 24.0)) * y * fa / (pv1 * pv1)
        pyi = py0 + pxi * f
        xi = x - ff
    pz1 = _tp_sqrt1m(tab, pxi, pyi)
    dpz1 = pz1 - 1.0
    dpx = irhoe * al
    pxf = pxi + dpx
    pz2 = _tp_sqrt1m(tab, pxf, pyi)
    dpz2 = pz2 - 1.0
    # u = asin の引数、w。K64 は pxi pxf < 0 (または d = 0) の枝 (phsq) と d ≠ 0 の枝 (d で割る) を数値の安定のために使い分ける。両者は同じ関数で、
    # d の枝は pxi + pxf = 0 (3 極の中央の磁石など) で 0/0 になり級数にできないので、常に phsq の枝を使う
    phsq = pxi * pxi + pz1 * pz1
    u = (dpx + pxf * dpz1 - pxi * dpz2) / phsq
    w = (pxf * dpz1 - pxi * dpz2 + dpx * pyi * pyi) / phsq
    dl = rhoe * (u.asin() - u + w)
    xo = xi + (pxf + pxi) / (pz2 + pz1) * al
    yo = y + pyi * (dl + al)
    if ex2:
        pxo, pyo = pxf, pyi
    else:
        pv2 = _tp_sqrt1m(tab, pxf)
        fa = yo * irhoe / pv2
        q = yo * (1.0 / rhob)
        f = (1.0 - q * q * (1.0 / 6.0)) * fa
        ff = (0.5 - q * q * (1.0 / 24.0)) * yo * fa / (pv2 * pv2)
        pyo = pyi - pxf * f
        xo = xo + ff
        pxo = pxf - xo * irhoe * tanp2
        pyo = pyo + yo * irhoe * tanp2
    pxoc, pyoc = pxo * pp, pyo * pp
    if ex2:
        xo, pxoc, yo, pyoc = _edge_exact(xo, pxoc, yo, pyoc, pp, -rhob, psi2, False)        # 出口の鎖は tblfre (dxfr2) の前
    if fb2 != 0.0:
        dxfr2 = -fb2 ** 2 / rhob / 24.0
        dyfr2 = fb2 / rhob ** 2 / 6.0
        dyfra2 = 4.0 * dyfr2 / fb2 ** 2 if fringe else 0.0
        xo = xo - dp * pinv * dxfr2
        pyoc = pyoc + (pinv * (dyfr2 - yo * yo * dyfra2)) * yo
    return [xo.a, pxoc.a, yo.a, pyoc.a]


def numeric_map(el, delta=0.0, momentum=3e9):
    """素子 el 1 台の pysad の追跡の写像 (n×4 → n×4、正準運動量。δ 固定)。模型の数値の追跡 (ampnf.model_tracker) と検証に使う"""
    from .track6n import track_particles_fast

    def f(Z):
        c = np.zeros((len(Z), 6))
        c[:, :4] = Z
        c[:, 5] = delta
        o = np.array(track_particles_fast([el], c, 1, momentum, 0.0, False, False, False))
        return o[:, :4]
    return f


def check_against_tracker(el, amp=1e-3, orbit=(0.0, 0.0, 0.0, 0.0), delta=0.0, D=5, sd=1e-2, s=1e-3, npts=24, seed=1):
    """K0≠0 の BEND 1 台のべき級数の写像を、追跡の写像と、入口の軌道 orbit のまわりの無作為な点 (振幅 amp) で比べる。戻り (最大の差, 軌道の出口の値)。
    級数は変数 u (x = orbit + s·u、δ = sd·δ̃) の次数 D までなので、差は (amp/長さの目盛り)^(D+1) の程度"""
    from . import ampnf as A
    tab = A.get_tab(5, D)
    S0 = [tab.const(orbit[j]) + s * tab.var(j) for j in range(4)]
    fn = bend_k0_map_for(el)
    S1 = fn(S0, tab, sd)
    rng = np.random.default_rng(seed)
    u = rng.uniform(-1, 1, (npts, 5)) * np.array([amp / s] * 4 + [delta / sd + 1e-3 / sd if delta == 0.0 else 0.0])
    u[:, 4] = delta / sd
    z = np.array(orbit) + s * u[:, :4]
    ref = numeric_map(el, delta)
    # δ は固定の値で比べる (u[:,4] = delta/sd)
    mono = np.prod(u[:, None, :] ** tab.mons[None, :, :], axis=2)
    got = np.column_stack([(mono * S1[j][None, :]).sum(axis=1).real for j in range(4)])
    want = ref(z)
    orb_out = ref(np.array([orbit], dtype=float))[0]
    return float(np.abs(got - want).max()), orb_out


def bend_k0_map_for(el):
    """ROTATE なしの素子の写像 (TPSA の係数配列のリスト S → S')。ROTATE があれば呼び側 (ampnf._expand_rolls) が先に展開する"""
    return bend_k0_map(el)


# ---------------------------------------------------------------------------------------------
# 据え付け (DX・DY・ROTATE・DROTATE・CHI2): track6n.tbshift1・tbrot1 の TPSA 版 (ANGLE≠0 の BEND)
# ---------------------------------------------------------------------------------------------
def _rot2(x, y, c, s):
    return x * c - y * s, x * s + y * c


def _tbshift(S, tab, sd, dx, dy, phi0, cost, sint, ent):
    """tbshift1 の TPSA 版 (追跡の内部の運動量 px/(1+δ) で書かれているので、正準運動量との変換を挟む)"""
    if dx == 0.0 and dy == 0.0 and sint == 0.0 and cost == 1.0:
        return S
    pp = TP(tab, tab.var(4, sd)) + 1.0
    pinv = pp.inv()
    x, y = TP(tab, S[0]), TP(tab, S[2])
    px, py = TP(tab, S[1]) * pinv, TP(tab, S[3]) * pinv
    if dx != 0.0 or dy != 0.0:
        phih = phi0 * 0.5
        sph = math.sin(phih)
        dcph = -2.0 * math.sin(phih * 0.5) ** 2
        st1 = sint if ent else -sint
        dxa = dx * cost - dy * st1
        dx1 = dx + dcph * cost * dxa
        dy1 = dy - dcph * st1 * dxa
        ds = dxa * sph
        if ent:
            al = _tp_sqrt1m(tab, px, py).inv() * ds
            x1 = x + px * al - dx1
            y = y + py * al - dy1
            x, y = _rot2(x1, y, cost, sint)
            px, py = _rot2(px, py, cost, sint)
        else:
            x1, y = _rot2(x, y, cost, sint)
            px, py = _rot2(px, py, cost, sint)
            al = _tp_sqrt1m(tab, px, py).inv() * ds
            x = x1 + px * al - dx1
            y = y + py * al - dy1
    elif sint != 0.0 or cost != 1.0:
        x, y = _rot2(x, y, cost, sint)
        px, py = _rot2(px, py, cost, sint)
    return [x.a, (px * pp).a, y.a, (py * pp).a]


def _xsincos(x):
    return math.sin(x), math.cos(x)


def _tbrot(S, tab, sd, alg, phig, dtheta, dchi2, ent):
    c0, s0 = math.cos(phig), math.sin(phig)
    s2, c2 = _xsincos(dchi2)
    st, ct = _xsincos(dtheta)
    dc2, dct = c2 - 1.0, ct - 1.0
    dr11 = dct * c0 * c0 + dc2 * s0 * s0 + c0 * s0 * s2 * st
    dr22 = c2 * dct + dc2
    dr33 = dc2 * c0 * c0 + dct * s0 * s0 - c0 * s0 * s2 * st
    r11, r22, r33 = dr11 + 1.0, dr22 + 1.0, dr33 + 1.0
    if ent:
        r12 = s0 * s2 - c0 * c2 * st
        r13 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r21 = -ct * s0 * s2 + c0 * st
        r23 = c0 * ct * s2 + s0 * st
        r31 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r32 = -c0 * s2 - c2 * s0 * st
    else:
        r21 = s0 * s2 - c0 * c2 * st
        r31 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r12 = -ct * s0 * s2 + c0 * st
        r32 = c0 * ct * s2 + s0 * st
        r13 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r23 = -c0 * s2 - c2 * s0 * st
    dl = alg if phig == 0.0 else alg * s0 / phig
    xl0 = dl * s0
    zi = -dl * c0
    dxz, dyz, dzz = r13 * zi, r23 * zi, dr33 * zi
    x, px, y, py = (TP(tab, a) for a in S)
    pzi = _tp_sqrt1m(tab, px, py)
    xi = x + xl0
    xf = xi * dr11 + y * r12 + x + dxz
    yf = xi * r21 + y * dr22 + y + dyz
    zf = xi * r31 + y * r32 + dzz
    pxf = px * r11 + py * r12 + pzi * r13
    pyf = px * r21 + py * r22 + pzi * r23
    pzf = px * r31 + py * r32 + pzi * r33
    ipz = pzf.inv()
    return [(xf - pxf * ipz * zf).a, pxf.a, (yf - pyf * ipz * zf).a, pyf.a]


def _wrap_shift(core, P):
    """tbshift(入口) → [tbrot] → 本体 → [tbrot] → tbshift(出口) (k_bend0・k_bendi の外側)。pxc↔内部の運動量の変換は本体の中"""
    phi0, L = P['phi0'], P['L']
    cost, sint = math.cos(P['theta']), math.sin(P['theta'])
    dx, dy, dth, dchi = P['dx'], P['dy'], P['dtheta'], P['dchi2']
    alg, phig = 0.5 * L, 0.5 * phi0
    al_core = core.__dict__.get('al', L)

    def f(S, tab, sd):
        S = _tbshift(S, tab, sd, dx, dy, phi0, cost, sint, True)
        if dth != 0.0 or dchi != 0.0:
            S = _tbrot_can(S, tab, sd, alg, phig, dth, dchi, True)
        S = core(S, tab, sd)
        if dth != 0.0 or dchi != 0.0:
            S = _tbrot_can(S, tab, sd, alg - L, phig - phi0, dth, dchi, False)
        return _tbshift(S, tab, sd, -dx, -dy, -phi0, cost, -sint, False)
    return f


def _tbrot_can(S, tab, sd, alg, phig, dth, dchi, ent):
    """tbrot は追跡の内部運動量 (px/(1+δ)) で書かれているので、正準運動量との変換を挟む"""
    dp = TP(tab, tab.var(4, sd))
    pp = dp + 1.0
    pinv = pp.inv()
    Si = [S[0], (TP(tab, S[1]) * pinv).a, S[2], (TP(tab, S[3]) * pinv).a]
    So = _tbrot(Si, tab, sd, alg, phig, dth, dchi, ent)
    return [So[0], (TP(tab, So[1]) * pp).a, So[2], (TP(tab, So[3]) * pp).a]


# ---------------------------------------------------------------------------------------------
# K1≠0 かつ K0≠0 (ANGLE≠0): track6n.k_bendi (tbendi の数値積分器) を TPSA で評価する
# ---------------------------------------------------------------------------------------------
def _ent_sums(z, a, nmax=9):
    """z (TP) のべき: dc = Σ_{n≥1} z^n a^{2n}/(2n)! (cosh − 1)、dck = dc/z、sk = Σ_{n≥0} z^n a^{2n+1}/(2n+1)! (sinh/√z)、
    zs = Σ_{n≥0} z^{n+1} a^{2n+1}/(2n+1)! (√z sinh = aks。z の符号に依らない整関数)"""
    tab = z.tab
    one = TP.const(tab, 1.0)
    zn = one
    dc = TP.const(tab, 0.0)
    dck = TP.const(tab, 0.0)
    sk = one * a
    zs = z * a
    for n in range(1, nmax + 1):
        zprev = zn
        zn = zn * z
        dc = dc + zn * (a ** (2 * n) / math.factorial(2 * n))
        dck = dck + zprev * (a ** (2 * n) / math.factorial(2 * n))
        sk = sk + zn * (a ** (2 * n + 1) / math.factorial(2 * n + 1))
        zs = zs + zn * z * (a ** (2 * n + 1) / math.factorial(2 * n + 1))
    return dc, dck, sk, zs


def _ttfrin4(x, px, y, py, ak, al, pinv):
    aki = pinv * (ak / al / 4.0)
    a = aki * x * x
    b = aki * y * y
    ab = a - b
    t = ab * ab * (1.0 / 6.0)
    dx1 = x * (a * (1.0 / 3.0) + b + t)
    dy1 = -(y * (a + b * (1.0 / 3.0) - t))
    xn = x + dx1
    yn = y + dy1
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) * (1.0 / 6.0)
    yy = 1.0 - d + ab * (a - 5.0 * b) * (1.0 / 6.0)
    h = aki * x * y * (1.0 - ab * (1.0 / 3.0)) * 2.0
    f = xx * yy + h * h
    fi = f.inv()
    return xn, (px * yy + py * h) * fi, yn, (py * xx - px * h) * fi


def _edge_exact(x, px, y, py, pr, rhob, psi, entrance):
    """track6n.edge_exact1 (tbfrie の tbedge: drift・maxwell・body) の TPSA 版。px・py は正準運動量 (p0 で規格化)。pr = 1 + δ (TP)"""
    cosp, sinp = math.cos(psi), math.sin(psi)
    sinsq = 2.0 * math.sin(0.5 * psi) ** 2

    def drift(x, px, y, py):
        pzi = (pr * pr - px * px - py * py).sqrt()
        f = 1.0 * cosp - px / pzi * sinp
        x = x / f
        sx = x * sinp
        return x, px * cosp + pzi * sinp, y + py / pzi * sx, py

    def maxwell(x, px, y, py, rb):
        pvi2 = (pr - px) * (pr + px)
        pvi = pvi2.sqrt()
        q = y * (1.0 / rb)
        yr = q * q
        af = (1.0 - yr * (1.0 / 6.0)) * y * (1.0 / rb) / pvi
        bf = (0.5 - yr * (1.0 / 24.0)) * y * y * (1.0 / rb) / (pvi2 * pvi)
        return x + bf * pr, px, y, py - af * px

    def body(x, px, y, py):
        s = px * px + py * py
        sqrs = (pr * pr - s).sqrt()
        dpzi = -s / (pr + sqrs)
        pzi = pr + dpzi
        phsq = pzi * pzi + px * px
        px0 = px * cosp - pzi * sinp
        dpz0 = px * sinp + dpzi * cosp - pr * sinsq
        xi = x
        pxf = px0 + xi * (sinp / rhob)
        s1 = py * py + pxf * pxf
        sqrs1 = (pr * pr - s1).sqrt()
        dpzf = -s1 / (pr + sqrs1)
        ddpz = (s - s1) / (sqrs + sqrs1) + (dpzi + pr) * sinsq - px * sinp
        xn = ddpz * rhob + xi * cosp
        a = pr * xi * (sinp / rhob) + dpz0 * pxf - dpzf * px0
        theta = (a / phsq).asin()
        return xn, pxf, y - theta * py * rhob, py
    if entrance:
        x, px, y, py = drift(x, px, y, py)
        x, px, y, py = maxwell(x, px, y, py, rhob)
        return body(x, px, y, py)
    x, px, y, py = body(x, px, y, py)
    x, px, y, py = maxwell(x, px, y, py, -rhob)
    return drift(x, px, y, py)


def _bendi_core(P, v2):
    from . import compat
    phi0, L, k0 = P['phi0'], P['L'], P['k0']
    phib = phi0 + k0
    if phib == 0.0:
        phib = 1e-9 * phi0
    ak = P['k1'] * L
    psi1, psi2 = P['psi1'], P['psi2']
    fb1, fb2 = P['fb1'], P['fb2']
    al = L
    if fb1 != 0.0 or fb2 != 0.0:           # k_bendi に渡る al は F1 の分だけ短い (Prog._compile)
        w_ = phi0 - psi1 - psi2
        al = L - ((phi0 * fb1) ** 2 + (phi0 * fb2) ** 2) / L / 48.0 * math.sin(0.5 * w_) / math.sin(0.5 * phi0)
    eps = 1e-3 if P['eps0'] == 0.0 else 1e-3 * P['eps0']
    ndiv = min(1000, 1 + int(abs(phi0 / eps)))
    fringe, mfring = P['fringe'], P['mfring']
    face = compat.bend_k1_edge() == 'face'
    klin = compat.bend_k1_lin()
    v2edge = bool(fringe and v2)

    def f(S, tab, sd):
        return _tbendi(S, tab, sd, al, phib, phi0, ak, psi1, psi2, fringe, mfring, v2edge, face, klin, ndiv, fb1, fb2)
    return f


def _tbendi(S, tab, sd, al, phib, phi0, ak, psi1, psi2, fringe, mfring, v2edge, face, klin, ndiv, fb1=0.0, fb2=0.0):
    cosp1, sinp1, cosp2, sinp2 = math.cos(psi1), math.sin(psi1), math.cos(psi2), math.sin(psi2)
    tanp1, tanp2 = sinp1 / cosp1, sinp2 / cosp2
    rhob = al / phib
    rho0 = al / phi0
    af = 1.0 if fringe else 0.0
    x, pxc, y, pyc = (TP(tab, a) for a in S)
    dp = TP(tab, tab.var(4, sd))
    pp = dp + 1.0
    pinv = pp.inv()
    px, py = pxc * pinv, pyc * pinv
    rhoe = pp * rhob
    irhoe = rhoe.inv()
    if fringe and mfring > -4 and mfring != 2:
        x, px, y, py = _ttfrin4(x, px, y, py, ak, al, pinv)
    if fb1 != 0.0:          # F1・FB1 の入口のフリンジ (k_bendi。dyfra1 は DISFRIN に依らない)
        dxfr1 = fb1 * fb1 / rhob / 24.0
        dyfr1 = fb1 / (rhob * rhob) / 6.0
        dyfra1 = 4.0 * dyfr1 / (fb1 * fb1)
        x = x + dp * pinv * dxfr1
        py = py + (pinv * pinv * (dyfr1 - y * y * dyfra1)) * y
    if v2edge and sinp1 != 0.0:
        xe, pxe, ye, pye = _edge_exact(x, px * pp, y, py * pp, pp, rhob, math.atan2(sinp1, cosp1), True)
        m_xi, m_pxi, m_yi, m_pyi = xe, pxe * pinv, ye, pye * pinv
        if face:
            kt = pinv * (ak / al * tanp1)
            m_pxi = m_pxi + kt * m_xi * m_xi
            m_pyi = m_pyi - kt * m_xi * m_yi
    else:
        f = y * irhoe
        fpx = px * af
        m_pyi = py - (fpx + tanp1) * f
        ff = y * f * (0.5 * af)
        m_xi = x + ff
        m_pxi = px + m_xi * irhoe * tanp1
        m_yi = y
    alx = al / ndiv
    akxn = ak * alx / al
    phixn = phi0 * alx / al
    rhosq = rhoe * rho0
    drhop = (rhoe - rho0) / rhosq
    akk = akxn / alx
    akxsq = -(rhosq.inv()) - pinv * akk
    akysq = pinv * akk
    dcx, dcxkx, sxkx, _ = _ent_sums(akxsq, alx)
    dcy, _d, syky, aksy = _ent_sums(akysq, alx)
    irho0 = 1.0 / rho0

    def corr(m_xi, m_pxi, m_yi, m_pyi, ak1c, alc1, phic):
        s = m_pxi * m_pxi + m_pyi * m_pyi
        pzi = (1.0 - s).sqrt()
        dpzi = pzi - 1.0
        m_xin = (pzi * m_xi - dpzi * m_pxi * alc1) / (pzi - m_pxi * phic)
        m_pxin = m_pxi + dpzi * phic
        m_yin = m_yi + (m_xin * phic - dpzi * alc1) * m_pyi / pzi
        xr = m_xin * irho0
        if klin:
            hk = pinv * (ak1c * irho0)
            m_pyin = m_pyi + hk * m_xin * m_yin
            m_pxin = m_pxin - hk * (m_xin * m_xin - m_yin * m_yin * 0.5)
        else:
            m_pyin = m_pyi
            m_pxin = m_pxin - m_xin * xr * (0.5 - xr * (2.0 - xr) * (1.0 / 12.0)) * pinv * ak1c
        return m_xin, m_pxin, m_yin, m_pyin
    akx0 = 0.0
    alx0 = 0.0
    phix0 = 0.0
    for m in range(1, ndiv + 1):
        m_xi, m_pxi, m_yi, m_pyi = corr(m_xi, m_pxi, m_yi, m_pyi, (akxn + akx0) * 0.5, (alx + alx0) * 0.5, (phixn + phix0) * 0.5)
        akx0, alx0, phix0 = akxn, alx, phixn
        xiksq = drhop + akxsq * m_xi
        dxf = dcxkx * xiksq + sxkx * m_pxi
        dpxf = sxkx * xiksq + dcx * m_pxi
        dyf = dcy * m_yi + syky * m_pyi
        dpyf = aksy * m_yi + dcy * m_pyi
        m_xi, m_pxi, m_yi, m_pyi = m_xi + dxf, m_pxi + dpxf, m_yi + dyf, m_pyi + dpyf
    m_xi, m_pxi, m_yi, m_pyi = corr(m_xi, m_pxi, m_yi, m_pyi, akx0 * 0.5, alx0 * 0.5, phix0 * 0.5)
    if v2edge and sinp2 != 0.0:
        if face:
            kt = pinv * (ak / al * tanp2)
            m_pxi = m_pxi + kt * m_xi * m_xi
            m_pyi = m_pyi - kt * m_xi * m_yi
        xo, pxo, yo, pyo = _edge_exact(m_xi, m_pxi * pp, m_yi, m_pyi * pp, pp, rhob, math.atan2(sinp2, cosp2), False)
        pxo, pyo = pxo * pinv, pyo * pinv
    else:
        pxo = m_pxi + m_xi * irhoe * tanp2
        f = m_yi * irhoe
        fpx = pxo * af
        pyo = m_pyi - (tanp2 - fpx) * f
        ff = m_yi * f * (0.5 * af)
        xo = m_xi - ff
        yo = m_yi
    if fb2 != 0.0:
        dxfr2 = fb2 * fb2 / rhob / 24.0
        dyfr2 = fb2 / (rhob * rhob) / 6.0
        dyfra2 = 4.0 * dyfr2 / (fb2 * fb2)
        xo = xo - dp * pinv * dxfr2
        pyo = pyo + (pinv * pinv * (dyfr2 - yo * yo * dyfra2)) * yo
    if fringe and mfring > -4 and mfring != 1:
        xo, pxo, yo, pyo = _ttfrin4(xo, pxo, yo, pyo, -ak, al, pinv)
    return [xo.a, (pxo * pp).a, yo.a, (pyo * pp).a]


# ---------------------------------------------------------------------------------------------
# ANGLE = 0 で K1≠0 かつ K0≠0 (直線の座標系の、双極子 + 四極の磁石)
#   K64 の追跡 (tbend → tmulti) は gs が SIGSEGV する組み合わせで、K64 の基準がない。pysad の tmulti は積分器の精度が 5e-7 (L = 0.5 m・K0 = 0.03 で、実験室の運動方程式の積分との差)。
#   ここは厳密な Hamiltonian H = −√((1+δ)² − px² − py²) + b0 x + k1 (x² − y²)/2 (b0 = K0/L、k1 = K1) を、対称な 4 次の Yoshida の分割 (ドリフト = 厳密、kick = 厳密) の TPSA で積分する
#   (シンプレクティック。スライス数 nstep で精度を決める。既定 24 で 1e-12 程度)。端のフリンジは四極の Lee–Whiting (DISFRIN = 0 のとき)。
# ---------------------------------------------------------------------------------------------
_W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
_W0 = -(2.0 ** (1.0 / 3.0)) * _W1


def _maxwell_edge(x, px, y, py, pr, rb):
    """硬い端の Maxwell の項 (track6n._edge_maxwell1。入口 rb = 1/b0、出口 −1/b0)。px・py は正準運動量"""
    pvi2 = (pr - px) * (pr + px)
    pvi = pvi2.sqrt()
    q = y * (1.0 / rb)
    yr = q * q
    af = (1.0 - yr * (1.0 / 6.0)) * y * (1.0 / rb) / pvi
    bf = (0.5 - yr * (1.0 / 24.0)) * y * y * (1.0 / rb) / (pvi2 * pvi)
    return x + bf * pr, px, y, py - af * px


def _straight_k1_core(P, nstep=24):
    L, k0, k1 = P['L'], P['k0'], P['k1']
    b0 = k0 / L
    ak = k1 * L
    fringe = P['fringe']
    maxwell = fringe          # 双極子の硬い端の Maxwell の項 (k_steer と同じ式。K64 の tmulti は DISFRIN = 0 のときだけ入れる)
    c = [0.5 * _W1, 0.5 * (_W0 + _W1), 0.5 * (_W0 + _W1), 0.5 * _W1]
    d = [_W1, _W0, _W1]

    def f(S, tab, sd):
        x, pxc, y, pyc = (TP(tab, a) for a in S)
        dp = TP(tab, tab.var(4, sd))
        pr = dp + 1.0
        if fringe:
            pinv = pr.inv()
            x, px_, y, py_ = _ttfrin4(x, pxc * pinv, y, pyc * pinv, ak, L, pinv)
            pxc, pyc = px_ * pr, py_ * pr
        if maxwell:
            x, pxc, y, pyc = _maxwell_edge(x, pxc, y, pyc, pr, L / k0)
        px, py = pxc, pyc
        ds = L / nstep
        for _ in range(nstep):
            for i in range(4):
                pz = (pr * pr - px * px - py * py).sqrt()
                ipz = pz.inv()
                x = x + px * ipz * (c[i] * ds)
                y = y + py * ipz * (c[i] * ds)
                if i < 3:
                    w = d[i] * ds
                    px = px - (x * k1 + b0) * w
                    py = py + y * k1 * w
        if maxwell:
            x, px, y, py = _maxwell_edge(x, px, y, py, pr, -L / k0)
        if fringe:
            pinv = pr.inv()
            x, px_, y, py_ = _ttfrin4(x, px * pinv, y, py * pinv, -ak, L, pinv)
            px, py = px_ * pr, py_ * pr
        return [x.a, px.a, y.a, py.a]
    return f
