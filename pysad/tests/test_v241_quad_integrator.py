# -*- coding: utf-8 -*-
"""v2.4.1: 四極の本体の積分器 (K64 の tsolqu・tsolque の bz=0 の枝、docs/21 B226) の共通化と、v2 の 4 次の合成。

K64 の四極の本体は、運動学的な非線形ドリフト D (運動量は変えない) と厳密な線形の本体 B を D(h/2) B(h) D(h/2) の Strang 分割 (2 次) で積分する。
刻みの数 1 + int(L·|K1|L/ε) は軌道のずれの大きさに依らないので、軌道が 3 mm・2 mrad ずれた四極 (L = 0.3 m、K1L = −0.6) のヤコビアンが厳密な写像と
CALC で 1.8e-7、追跡で 7e-7 違う。v2 は S(w1 h) S(w0 h) S(w1 h) (Yoshida の 4 次) にする。k64 は従来のまま (ビット単位で同じ)。

独立な基準 (pysad のコードを使わない): 運動学的な運動方程式 x' = px/pz、px' = −k x、y' = py/pz、py' = +k y、z' = −(1/pz − 1) (pz = √(1 − px² − py²)、
運動量は 1+δ で規格化、k = K1/L/(1+δ)) と、その変分方程式 (ヤコビアン) を solve_ivp (DOP853、rtol 1e-13) で積分する。
限界: δ の行・列と z の列の微分、K1 の x・y の結合 (SK1)、端のフリンジ (tqfrie・tqlfre)、放射のスライス (tradk・tradke) は、ここでは見ない
(刻みの構造は共通なので、放射つきの経路は k64 のビット同一と v2 の python 版・numba 版の一致で確認する)。
"""
import math
import os

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, lattice as lat, quadbody as qb
from . import _p3_indep as X

pytestmark = pytest.mark.filterwarnings('ignore')


def _mode(m):
    old = compat._mode
    compat.set_mode(m)
    return old


# ---------------------------------------------------------------- 独立な基準 (運動学的な運動方程式の数値積分)
def ref_map(c, L, k):
    """(x, px, y, py) の厳密な写像と z。k = K1 の勾配 (K1L/L、x に収束で正)。δ = 0"""
    def rhs(s, v):
        x, px, y, py, z = v
        pz = math.sqrt(1.0 - px * px - py * py)
        return [px / pz, -k * x, py / pz, k * y, -(1.0 / pz - 1.0)]
    sol = solve_ivp(rhs, (0.0, L), list(c[:4]) + [0.0], method='DOP853', rtol=1e-13, atol=1e-16)
    return sol.y[:, -1]


def ref_jac(c, L, k):
    """4×4 のヤコビアン (変分方程式を同時に積分)"""
    def rhs(s, v):
        x, px, y, py = v[:4]
        M = v[4:].reshape(4, 4)
        pz = math.sqrt(1.0 - px * px - py * py)
        A = np.zeros((4, 4))
        A[0, 1] = 1 / pz + px * px / pz ** 3
        A[0, 3] = px * py / pz ** 3
        A[2, 3] = 1 / pz + py * py / pz ** 3
        A[2, 1] = px * py / pz ** 3
        A[1, 0] = -k
        A[3, 2] = k
        return np.concatenate([[px / pz, -k * x, py / pz, k * y], (A @ M).ravel()])
    v0 = np.concatenate([c[:4], np.eye(4).ravel()])
    sol = solve_ivp(rhs, (0.0, L), v0, method='DOP853', rtol=1e-13, atol=1e-16)
    return sol.y[4:, -1].reshape(4, 4)


_ROT_IN = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]], float)     # K1 < 0: (x, px, y, py) → (y, py, −x, −px)


def ref_quad_map(c, L, KL):
    """QUAD (L, K1L = KL) の (x, px, y, py)。KL < 0 は −90° 回転して |K1| で通し、戻す (tquade)"""
    k = abs(KL) / L
    if KL < 0:
        out = ref_map(_ROT_IN @ np.asarray(c[:4]), L, k)
        return _ROT_IN.T @ out[:4]
    return ref_map(c, L, k)[:4]


def ref_quad_jac(c, L, KL):
    k = abs(KL) / L
    if KL < 0:
        return _ROT_IN.T @ ref_jac(_ROT_IN @ np.asarray(c[:4]), L, k) @ _ROT_IN
    return ref_jac(c, L, k)


def _quad_el(L, KL):
    eld = 'QUAD SA=(L=%g K1=%g DISFRIN=1);' % (L, KL)
    f1, _ = X.run('MOMENTUM = 3GEV;\n%s\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n' % eld)
    return eld, [e for e in f1.lat.elements if e.name == 'SA'][0]


# ---------------------------------------------------------------- 刻みの構造 (quadbody.group)
@pytest.mark.parametrize('sch', [0, 1])
@pytest.mark.parametrize('nd', [1, 2, 5])
def test_group_structure(sch, nd):
    """D の長さの和 = 本体の長さの和 = L。放射点は刻みの最後の B のあと (最後の刻みを除く)、最後の D のあとにも 1 つ (CALC の規約)。
    隣り合う刻みの D は 1 つにまとまっている (D の数 = m·nd + 1)"""
    L = 0.7
    aln1 = L / nd
    m = qb.QB_M[sch]
    assert qb.ngroup(nd, sch) == m * nd + 1
    ds, bs, rads = [], [], 0
    for g in range(qb.ngroup(nd, sch)):
        dl, bl, jb, rad, fin = qb.group(g, aln1, nd, sch)
        ds.append(dl)
        bs.append(bl)
        rads += bool(rad)
        assert fin == (g == m * nd)
    assert abs(sum(ds) - L) < 1e-15 and abs(sum(bs) - L) < 1e-15
    assert rads == nd      # 刻みの最後 (nd − 1 個) + 最後の D
    if sch == 0:
        # K64: D = [0.5, 1, ..., 1, 0.5] × aln1、B = aln1 (ビット単位で従来と同じ)
        assert ds == [aln1 * 0.5] + [aln1] * (nd - 1) + [aln1 * 0.5] and bs == [aln1] * nd + [0.0]
    else:
        # 4 次の合成: B の重みは (w1, w0, w1)、w0 < 0。合計は 1
        w = qb.QB_BW[1]
        assert abs(sum(w) - 1.0) < 1e-15 and w[1] < 0 < w[0]
        assert abs(bs[1] - aln1 * w[1]) < 1e-15


def test_scheme_follows_compat_and_override():
    old = _mode('k64')
    try:
        assert qb.scheme() == 0
        compat.set_mode('v2')
        assert qb.scheme() == 1
        prev = qb.set_scheme(0)
        assert qb.scheme() == 0
        qb.set_scheme(prev)
        from pysad import track6n
        assert track6n.make_rs(1)[4] == 1
        compat.set_mode('k64')
        assert track6n.make_rs(1)[4] == 0
    finally:
        compat.set_mode(old)


# ---------------------------------------------------------------- 追跡: 1 要素の写像
@pytest.mark.parametrize('L,KL,c,tol_v2', [(0.3, -0.6, [0.01, 0.005, 0.004, -0.004, 0, 0], 1e-8),
                                           (0.5, 1.5, [0.01, 0.005, 0.004, -0.004, 0, 0], 1e-9),
                                           (1.0, 1.5, [0.01, 0.005, 0.004, -0.004, 0, 0], 2e-10)])
@pytest.mark.parametrize('fast', [False, True])
def test_tracking_quad_map_vs_exact(L, KL, c, tol_v2, fast):
    """追跡 (python 版・numba 版) の四極 1 台の写像と、運動方程式の数値積分の差: k64 は K64 のまま (1e-9 以上)、v2 は k64 の 1/4 以下で tol_v2 以内"""
    eld = 'QUAD SA=(L=%g K1=%g DISFRIN=1);' % (L, KL)
    r = ref_quad_map(c, L, KL)
    errs = {}
    for mode in ('k64', 'v2'):
        old = _mode(mode)
        try:
            o = X.one_element_track(eld, c, False, fast)
            errs[mode] = np.abs(o[:4] - r).max()
        finally:
            compat.set_mode(old)
    assert errs['k64'] > 1e-9, errs                  # K64 の誤りは再現のまま
    assert errs['v2'] < tol_v2 and errs['v2'] < errs['k64'] / 4, errs


@pytest.mark.parametrize('eld', ['QUAD SA=(L=0.3 K1=-0.6 DISFRIN=1);', 'QUAD SA=(L=1.0 K1=1.5);', 'MULT SA=(L=0.6 K1=-0.6 K0=0.003 K2=-30 FRINGE=3 F1=0.03);',
                                 'MULT SA=(L=0.3 K1=0.7 K2=20);'])
@pytest.mark.parametrize('rad', [False, True])
def test_tracking_python_and_numba_agree_v2(eld, rad):
    """v2 の追跡の python 版 (track6) と numba 版 (track6n・track6n_ext) が一致する (放射なし・あり)"""
    old = _mode('v2')
    try:
        for c in ([0.01, 0.004, -0.008, 0.003, 0.0, 0.005], [-0.002, 0.0008, 0.004, -0.002, 0.0, -0.01]):
            a = X.one_element_track(eld, c, rad, False)
            b = X.one_element_track(eld, c, rad, True)
            assert np.abs(a - b).max() < 1e-12, (a, b)
    finally:
        compat.set_mode(old)


def test_tracking_fourth_order_convergence():
    """刻みの数を 2 倍にしたときの誤りの比: k64 は 2 次 (≈ 4)、v2 は 4 次 (≈ 16)。EPS (tsolqu の ε の倍率) で刻みの数を変える。
    四極 L = 1 m、K1L = 1.5、x・y が 10 mm、傾きが 10 mrad (1 刻みの角の 2 乗 k·aln² = 1.5/nd²)"""
    from pysad import track6 as t6
    c = np.array([0.01, 0.01, -0.008, 0.008, 0.0])
    L, KL = 1.0, 1.5
    r = ref_map(c, L, KL / L)
    for sch, lo, hi in ((0, 3.0, 4.5), (1, 12.0, 20.0)):
        errs = []
        for eps0 in (8.0, 4.0, 2.0, 1.0):               # ndiv = 1, 2, 4, 8
            o = np.array(t6._qb_nonrad(*c, 0.0, 0.0, L, KL, eps0, 0.0, 0.0, False, sch))
            errs.append(np.abs(o - r).max())
        ratios = [errs[i] / errs[i + 1] for i in range(1, 3)]       # ndiv = 2→4、4→8 (1 刻みの角が大きい 1→2 は漸近に入っていない)
        assert all(lo < q < hi for q in ratios), (sch, errs, ratios)


# ---------------------------------------------------------------- CALC: ヤコビアン
@pytest.mark.parametrize('elkern', [True, False])
@pytest.mark.parametrize('L,KL,c', [(0.3, -0.6, [0.003, 0.002, 0.0, 0.0, 0.0, 0.0]), (0.3, 0.6, [0.003, 0.002, 0.0, 0.0, 0.0, 0.0]),
                                    (0.5, 1.5, [0.01, 0.005, 0.004, -0.004, 0.0, 0.0])])
def test_calc_jacobian_vs_variational(elkern, L, KL, c):
    """CALC (numba 版 elkern と numpy 版) の四極 1 台の解析ヤコビアンと、変分方程式の数値積分の差: k64 は 2e-8〜4e-7 (K64 のまま)、
    v2 は 1e-9 以下 (4 次の合成 + CALC の ε を半分)"""
    eld, el = _quad_el(L, KL)
    Mref = ref_quad_jac(np.array(c), L, KL)
    old_ek = lat._USE_ELKERN
    lat._USE_ELKERN = elkern
    errs = {}
    try:
        for mode in ('k64', 'v2'):
            old = _mode(mode)
            try:
                M, _ = lat._el_cod_step(el, np.array(c, float), 0.0)
                errs[mode] = np.abs(M[:4, :4] - Mref).max()
            finally:
                compat.set_mode(old)
    finally:
        lat._USE_ELKERN = old_ek
    assert 1e-8 < errs['k64'] < 1e-6, errs
    assert errs['v2'] < 1e-9 and errs['v2'] < errs['k64'] / 20, errs


def test_calc_numba_and_numpy_agree_v2():
    """v2 の CALC の numba 版 (elkern) と numpy 版が、四極 1 台の COD と 4×4 の行列で 1e-12 以内に一致する (z の行は 4D の CALC では使わない)"""
    eld, el = _quad_el(0.5, 1.5)
    c = np.array([0.01, 0.005, 0.004, -0.004, 0.0, 0.003])
    old = _mode('v2')
    old_ek = lat._USE_ELKERN
    try:
        res = []
        for ek in (True, False):
            lat._USE_ELKERN = ek
            res.append(lat._el_cod_step(el, c, 0.003))
        assert np.abs(res[0][0][:4, :4] - res[1][0][:4, :4]).max() < 1e-12 and np.abs(res[0][1][:4] - res[1][1][:4]).max() < 1e-12
    finally:
        lat._USE_ELKERN = old_ek
        compat.set_mode(old)


def test_calc_jacobian_is_symplectic_v2():
    """4 次の合成 (w0 < 0 の逆向きの本体を含む) でも CALC の 4×4 ヤコビアンはシンプレクティック (|M^T J M − J| < 1e-12)"""
    eld, el = _quad_el(0.3, -0.6)
    old = _mode('v2')
    try:
        M, _ = lat._el_cod_step(el, np.array([0.003, 0.002, 0.001, -0.001, 0.0, 0.0]), 0.0)
    finally:
        compat.set_mode(old)
    J = np.zeros((4, 4))
    J[0, 1] = J[2, 3] = 1.0
    J -= J.T
    M4 = M[:4, :4]
    assert np.abs(M4.T @ J @ M4 - J).max() < 1e-12


def test_k64_calc_and_tracking_jacobians_keep_k64_error():
    """k64: 軌道が 3 mm・2 mrad ずれた四極 (L = 0.3、K1L = −0.6) の追跡と CALC のヤコビアンの厳密な写像との差は、K64 のまま
    (追跡 7.1e-7、CALC 1.8e-7)。v2 は CALC が 5e-10、追跡が 1/50 以下 (8e-9) (docs/21 B226)"""
    L, KL = 0.3, -0.6
    eld, el = _quad_el(L, KL)
    c0 = np.array([0.003, 0.002, 0, 0, 0, 0.0])
    Mref = ref_quad_jac(c0, L, KL)
    res = {}
    for mode in ('k64', 'v2'):
        old = _mode(mode)
        try:
            h = 1e-6
            Jt = np.zeros((4, 4))
            for j in range(4):
                cp, cm = c0.copy(), c0.copy()
                cp[j] += h
                cm[j] -= h
                Jt[:, j] = (X.one_element_track(eld, cp, False, True)[:4] - X.one_element_track(eld, cm, False, True)[:4]) / (2 * h)
            Mc, _ = lat._el_cod_step(el, c0, 0.0)
            res[mode] = (np.abs(Jt - Mref).max(), np.abs(Mc[:4, :4] - Mref).max())
        finally:
            compat.set_mode(old)
    assert abs(res['k64'][0] / 7.1e-7 - 1) < 0.1 and abs(res['k64'][1] / 1.83e-7 - 1) < 0.1, res
    assert res['v2'][1] < 1e-9 and res['v2'][0] < res['k64'][0] / 50, res
