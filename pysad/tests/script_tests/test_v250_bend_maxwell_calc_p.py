# -*- coding: utf-8 -*-
"""v2.5.0: 偏向磁石の端の Maxwell フリンジ (tbfrie.f tbedgemaxwell) の COD の更新の因子 p = 1+δ (docs/30 B324)。

K64 の CALC の COD の更新は Δx = bf·p、Δz = −bf·px (bf = b y²/(ρ pv³))。追跡 (tbendcore) の式・(5.36)・CALC 自身の解析行列 trans は Δx = bf·p²、Δz = −bf·px·p。
因子 p が 1 つ足りないので、δ ≠ 0 で y ≠ 0 のとき CALC の COD の更新が追跡の写像と食い違う (ずれは δ y²/(2ρ) 程度)。
v2 は COD の更新を p² 型に直した (lattice._tbedgemaxwell_cod、elkern._tbedgemaxwell、追跡側の track6n._edge_maxwell1・ampnf_k0)。k64 は K64 のまま (gs と同じ)。

独立な基準: (5.36) の閉じた式、COD の更新の数値ヤコビアンと CALC の解析行列 (v2 で一致・シンプレクティック)、
gs (SAD 1.4.4.2k64、WSL) の保存した出力 (scripts/v250_bend_maxwell_calc_p__ring_*.gs.txt)。
gs の台本: 4 台の BEND (ANGLE ±0.2、L = 1 m、ρ = ±5 m、FRINGE = 1、DISFRIN = 0) + QUAD 2 台 (DY = ±3 mm で垂直の閉軌道) の 1 周 (3 GeV/c、NORFSW)。
CALC の閉軌道 c (Twiss DX..DDP) から TrackParticles を 1 周 → c を追跡の写像の不動点 c + dc に直すニュートン 1 回 (dc = (1−M)⁻¹ (T(c) − c)、M は追跡から求めた 4×4)
を出力する (`DCOD`)。dc は CALC の閉軌道と追跡の閉軌道の差。δ = 0 では 1e-12 以下、δ ≠ 0 で k64 は 1e-8〜7e-7。
限界: 4 次元 (x, px, y, py) の閉軌道。Δz は 4 次元の CALC に効かないので解析行列と数値微分で見る。v2 の dc の残り (δ = 0.2 で x 3e-8) は端の Maxwell 以外の δ の項で、ここでは原因を特定していない。
"""
import io
import math
import os
import re

import numpy as np
import pytest

from pysad import compat, lattice as lat, elkern
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')


def _mode(m):
    compat.set_mode(m)


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


# ---------------------------------------------------------------- (i) 閉じた式 (5.36)
RHO = 5.0


def closed_form(cod, rho, p_power):
    """(5.36): Δx = b y² p^pp/(ρ pv³)、Δpy = −a y px/(ρ pv)、Δz = −b y² px p^(pp−1)/(ρ pv³)。p_power = 2 が正しい式 (追跡・ヘルプ)、1 が K64 の CALC"""
    x, px, y, py, z, dl = cod
    p = 1.0 + dl
    pv2 = p * p - px * px
    pv = math.sqrt(pv2)
    a = 1.0 - (y / rho) ** 2 / 6.0
    b = 0.5 - (y / rho) ** 2 / 24.0
    bf = b * y * y / rho / pv ** 3
    return bf * p ** p_power, -a * y / rho * px / pv, -bf * px * p ** (p_power - 1)


CODS = [np.array(v) for v in ([0.0, 0.02, 4e-3, 1e-3, 0.0, 0.05], [1e-3, 0.03, -8e-3, 2e-3, 0.0, 0.2], [0.0, -0.01, 6e-3, 0.0, 0.0, -0.1])]


@pytest.mark.parametrize('rho', [RHO, -RHO])
@pytest.mark.parametrize('cod', CODS, ids=lambda c: 'dl=%g' % c[5])
def test_cod_update_matches_closed_form(rho, cod):
    """v2 は (5.36) の p² 型 (Δx・Δz)、k64 は p¹・p⁰ 型。Δpy は両方で同じ (正しい)"""
    for m, pw in (('v2', 2), ('k64', 1)):
        _mode(m)
        o = lat._tbedgemaxwell_cod(cod, rho)
        dx, dpy, dz = closed_form(cod, rho, pw)
        assert abs((o[0] - cod[0]) - dx) < 1e-15, (m, o[0] - cod[0], dx)
        assert abs((o[3] - cod[3]) - dpy) < 1e-15, m
        assert abs((o[4] - cod[4]) - dz) < 1e-15, m
    # k64 と v2 の差は Δx・Δz とも比 = p
    _mode('k64')
    k = lat._tbedgemaxwell_cod(cod, rho)
    _mode('v2')
    v = lat._tbedgemaxwell_cod(cod, rho)
    p = 1.0 + cod[5]
    assert (v[0] - cod[0]) / (k[0] - cod[0]) == pytest.approx(p, rel=1e-12)
    assert (v[4] - cod[4]) / (k[4] - cod[4]) == pytest.approx(p, rel=1e-12)


def test_delta_zero_same_in_both_modes_bitwise():
    """δ = 0 (p = 1) では v2 と k64 が同じ。他の項 (Δpy・クリップ) も変えていない"""
    for rho in (RHO, -RHO):
        for c in ([0.0, 0.02, 4e-3, 1e-3, 0.0, 0.0], [1e-3, -0.03, -8e-3, 2e-3, 1e-6, 0.0]):
            c = np.array(c)
            _mode('k64')
            k = lat._tbedgemaxwell_cod(c, rho)
            _mode('v2')
            v = lat._tbedgemaxwell_cod(c, rho)
            assert np.array_equal(k, v)


# ---------------------------------------------------------------- (ii) 解析行列との整合・シンプレクティック性
SG = np.zeros((6, 6))
for _i in range(3):
    SG[2 * _i, 2 * _i + 1] = 1.0
    SG[2 * _i + 1, 2 * _i] = -1.0


def _jac(f, c, h=1e-6):
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (f(cp) - f(cm)) / (2 * h)
    return M


@pytest.mark.parametrize('dl', [0.0, 0.1, -0.1])
def test_v2_update_is_consistent_with_analytic_matrix_and_symplectic(dl):
    """v2: COD の更新の数値ヤコビアン = CALC の解析行列 trans (1e-9 以下)、シンプレクティック (|MᵀSM − S| 1e-9 以下)。
    k64: 解析行列と 5e-7 以上食い違い、シンプレクティックでない (δ = 0 でも Δz の δ 列で 1e-6)"""
    c = np.array([0.001, 0.02, 0.01, 0.01, 0.0, dl])
    Ma = lat._tbedgemaxwell_mat(c, 10.0)
    sy = lambda M: np.abs(M.T @ SG @ M - SG).max()
    _mode('v2')
    Mv = _jac(lambda cc: lat._tbedgemaxwell_cod(cc, 10.0), c)
    assert np.abs(Mv - Ma).max() < 1e-9
    assert sy(Mv) < 1e-9
    _mode('k64')
    Mk = _jac(lambda cc: lat._tbedgemaxwell_cod(cc, 10.0), c)
    assert np.abs(Mk - Ma).max() > 5e-7
    assert sy(Mk) > 5e-7


# ---------------------------------------------------------------- (iii) numba 版 (CALC の k_bend) と python 版の一致
@pytest.mark.parametrize('p2', [False, True])
def test_numba_kernel_equals_python(p2):
    """elkern._tbedgemaxwell(c, rhob, p2) の COD と行列が lattice の (k64: p2 = False、v2: p2 = True) と一致"""
    for rho in (RHO, -RHO):
        for c0 in CODS:
            _mode('v2' if p2 else 'k64')
            ref_c = lat._tbedgemaxwell_cod(c0, rho)
            ref_M = lat._tbedgemaxwell_mat(c0, rho)
            c = c0.copy()
            M = elkern._tbedgemaxwell(c, rho, p2)
            assert np.abs(c - ref_c).max() < 1e-15
            assert np.abs(M - ref_M).max() < 1e-13


def _one_bend(psi):
    src = ('MOMENTUM = 3GEV;\nBEND B1 = (L=1.0 ANGLE=0.2 E1=%g E2=%g FRINGE=1 DISFRIN=0);\nMARK ST = ();\nLINE ONE = (ST B1);\nFFS;\nUSE ONE;\n' % (psi, psi))
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return [e for e in ffs.lat.elements if e.name == 'B1'][0]


@pytest.mark.parametrize('psi', [0.0, 0.05])
def test_calc_bend_step_numba_equals_python_chain(psi):
    """BEND の CALC (numba の k_bend) の出口の COD (x, px, y, py) が、python の鎖 (lattice._bend_cod_update = tbedge → 本体 → tbedge) と一致 (k64・v2 とも。δ ≠ 0、y ≠ 0)"""
    for m in ('k64', 'v2'):
        _mode(m)
        el = _one_bend(psi)
        c = np.array([1e-3, 0.01, 6e-3, 2e-3, 0.0, 0.2])
        fast = lat._bend_step_fast(el, c.copy())
        assert fast is not None
        py = lat._bend_cod_update(c.copy(), el.L, el.ANGLE, el.E1, el.E2, el.ANGLE, *lat._bend_fringe_params(el))
        assert np.abs(np.asarray(fast[0])[:4] - py[:4]).max() < 1e-12, (m, psi)


# ---------------------------------------------------------------- (iv) gs の保存した出力との比較
def _read(name):
    return open(os.path.join(SCRIPTS, name), encoding='utf-8', errors='replace').read()


def _vec(t, key):
    m = re.search(r'^%s \{([^}]*)\}' % key, t, re.M)
    return np.array([float(v) for v in m.group(1).replace('\n', ' ').split(',')])


RINGS = {'b1': 'dl=0.05 E=0', 'b2': 'dl=0.2 E=0', 'b3': 'dl=0', 'c1': 'dl=0.2 E=0.05', 'c2': 'dl=-0.1 E=0.05', 'd1': 'dl=0.1 y=6mm'}


def _pysad(name, mode):
    _mode(mode)
    buf = io.StringIO()
    run_sad_string(_read('v250_bend_maxwell_calc_p__ring_%s.sad' % name), out=buf)
    t = buf.getvalue()
    return _vec(t, 'COD'), _vec(t, 'DCOD')


@pytest.mark.parametrize('name', list(RINGS), ids=list(RINGS.values()))
def test_k64_matches_gs_saved_output(name):
    """k64: CALC の閉軌道 COD と、CALC の閉軌道と追跡の閉軌道の差 DCOD が gs の保存した出力と一致 (COD 1e-12、DCOD 1e-13 以下)"""
    g = _read('v250_bend_maxwell_calc_p__ring_%s.gs.txt' % name)
    cod, dcod = _pysad(name, 'k64')
    assert np.abs(cod - _vec(g, 'COD')).max() < 1e-12
    assert np.abs(dcod - _vec(g, 'DCOD')).max() < 1e-13


def test_gs_calc_vs_tracking_cod_disagree_only_for_nonzero_delta():
    """gs (保存した出力): CALC の閉軌道と追跡の閉軌道の差 (x, px) は δ = 0 で 1e-12 以下、δ ≠ 0 で px が 8e-8〜7e-7"""
    d = {n: _vec(_read('v250_bend_maxwell_calc_p__ring_%s.gs.txt' % n), 'DCOD') for n in RINGS}
    assert np.abs(d['b3'][:2]).max() < 1e-12
    for n in ('b1', 'b2', 'c1', 'c2', 'd1'):
        assert abs(d[n][1]) > 7e-8, n
    assert abs(d['d1'][1]) > 6e-7 and abs(d['b2'][1]) > 3e-7


@pytest.mark.parametrize('name', ['b1', 'b2', 'c1', 'c2', 'd1'])
def test_v2_calc_cod_agrees_with_tracking_better(name):
    """v2: CALC の閉軌道と追跡の閉軌道の差の px が、gs (k64 と同じ) の 12% 以下に減る (δ = 0.05 の x は 1e-9 以下)"""
    g = _vec(_read('v250_bend_maxwell_calc_p__ring_%s.gs.txt' % name), 'DCOD')
    _, dv = _pysad(name, 'v2')
    assert abs(dv[1]) < 0.12 * abs(g[1]), (name, dv[1], g[1])
    if name == 'b1':
        assert abs(dv[0]) < 1e-9


def test_v2_delta_zero_calc_tracking_consistent():
    """v2・δ = 0: CALC の閉軌道と追跡の閉軌道の差 (x, px) が 1e-12 以下 (端の更新を直しても δ = 0 は変わらない)"""
    _, dv = _pysad('b3', 'v2')
    assert np.abs(dv[:2]).max() < 1e-12
