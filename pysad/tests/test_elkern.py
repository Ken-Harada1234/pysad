# -*- coding: utf-8 -*-
"""elkern (要素ごとの numba カーネル, docs/history/work_perf2.md) が従来の Python の piece 経路と同じ結果を出すことの回帰テスト。

(1) 要素 1 個 (ドリフト・SEXT (ずれ・ロール)・QUAD (ずれ・ロール, K1<0)・BEND (K1=0/K1≠0, DX/DY/DROTATE/K0)・ステアリング) を
    乱数の cod から 6D (放射・D/B 付き, radcod 有/無) で通し、cod/T/D/B/u0 を比べる
(2) 4D の _el_cod_step (CALC) を同じ要素で比べる
(3) ESRF (誤差入り) の 1 周 tturne と Emittance を Python 経路と比べる
ステアリング (tsteee) は Python 版が数値微分なので行列の許容を緩める。"""
import os

import numpy as np
import pytest

from pysad import lattice as _lat
from pysad import emit6d as E
from pysad.elements import Drift, Quadrupole, Bend, Sextupole

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
PH = E._Phys(3e9, 603.12)


def _els():
    out = []
    out.append(('drift', Drift(name='D', L=0.7)))
    s = Sextupole(name='S', L=0.2, K2=12.0); s.DX = 3e-5; s.DY = -4e-5; s.ROTATION = 2e-4
    out.append(('sext_err', s))
    s2 = Sextupole(name='S2', L=0.0, K2=3.0); s2.DX = 2e-5
    out.append(('sext_thin', s2))
    q = Quadrupole(name='Q', L=0.4, K1=1.2); q.DX = 4e-5; q.DY = 3e-5; q.ROTATION = -1e-4
    out.append(('quad_err', q))
    q2 = Quadrupole(name='Q2', L=0.3, K1=-0.9); q2.DX = -2e-5
    out.append(('quad_neg', q2))
    q3 = Quadrupole(name='Q3', L=0.3, K1=0.8)
    out.append(('quad_plain', q3))
    b = Bend(name='B', L=1.0, ANGLE=0.05, E1=0.01, E2=0.015); b.K0 = 1e-4; b.DROTATE = 1e-4
    out.append(('bend_k0_dth', b))
    b2 = Bend(name='B2', L=0.8, ANGLE=0.04, E1=0.0, E2=0.0, K1=-0.3); b2.DX = 5e-5; b2.DY = -3e-5
    out.append(('bend_k1_dxdy', b2))
    b3 = Bend(name='B3', L=0.6, ANGLE=0.03, E1=0.005, E2=0.005)
    out.append(('bend_plain', b3))
    b4 = Bend(name='HD', L=0.1, ANGLE=0.0, E1=0.0, E2=0.0); b4.K0 = 2e-4
    out.append(('steer', b4))
    b5 = Bend(name='VD', L=0.0, ANGLE=0.0); b5.K0 = -1e-4
    out.append(('steer_thin', b5))
    return out


def _cods():
    rng = np.random.default_rng(5)
    c = rng.normal(size=(3, 6)) * np.array([3e-4, 2e-4, 3e-4, 2e-4, 1e-3, 5e-4])
    c[0] = [1e-4, 2e-4, -1e-4, 1e-4, 2e-3, 1e-3]
    return c


def _run6(el, cod, fast, radcod):
    E._USE_ELKERN = fast
    try:
        st = E.tturne([el], cod, PH, rfsw=False, radcod=radcod, want_db=True, want_jac=True)
    finally:
        E._USE_ELKERN = True
    return st


@pytest.mark.parametrize('radcod', [False, True])
@pytest.mark.parametrize('name,el', _els())
def test_element_6d_matches_python(name, el, radcod):
    tol = 1e-6 if name.startswith('steer') else 1e-10
    for cod in _cods():
        a = _run6(el, cod, False, radcod)
        b = _run6(el, cod, True, radcod)
        np.testing.assert_allclose(b.cod, a.cod, rtol=tol, atol=1e-16, err_msg=name + ' cod')
        np.testing.assert_allclose(b.T, a.T, rtol=tol, atol=tol * 1e-3, err_msg=name + ' T')
        np.testing.assert_allclose(b.D, a.D, rtol=tol, atol=1e-18, err_msg=name + ' D')
        np.testing.assert_allclose(b.B, a.B, rtol=tol, atol=1e-18, err_msg=name + ' B')
        assert abs(b.u0 - a.u0) <= 1e-12 * max(1.0, abs(a.u0)), name


@pytest.mark.parametrize('name,el', _els())
def test_element_4d_matches_python(name, el):
    tol = 1e-10
    for cod in _cods():
        _lat._USE_ELKERN = False
        try:
            Ma, ca = _lat._el_cod_step(el, cod.copy(), 0.0)
        finally:
            _lat._USE_ELKERN = True
        Mb, cb = _lat._el_cod_step(el, cod.copy(), 0.0)
        # z は K1≠0 BEND で 4D の Python 版が仮の運動量 (1 GeV/c) の dv 項を入れていたので 1e-6 (CALC は z を使わない)
        ix = [0, 1, 2, 3, 5]
        np.testing.assert_allclose(cb[ix], ca[ix], rtol=tol, atol=1e-16, err_msg=name + ' cod')
        np.testing.assert_allclose(cb[4], ca[4], rtol=1e-6, atol=1e-9, err_msg=name + ' z')
        # 4D の CALC が使う x,px,y,py の行 (z 行は使わない。K1≠0 BEND の z は 4D では dv 項の扱いが違う)
        np.testing.assert_allclose(Mb[:4], Ma[:4], rtol=tol, atol=1e-14, err_msg=name + ' M')


@pytest.fixture(scope='module')
def esrf_err():
    from pysad.ffscmd import FFS
    from pysad import parser
    pytest.importorskip('numba')
    # 誤差入りのラティス (SAD スクリプト) は環境変数 PYSAD_ERR_LATTICE で指定する (リポジトリには含めていない)
    path = os.environ.get('PYSAD_ERR_LATTICE')
    if not path or not os.path.exists(path):
        pytest.skip('環境変数 PYSAD_ERR_LATTICE (誤差入りラティスの .sad) がない')
    src = open(path, encoding='utf-8').read()
    from pysad.ffscmd import run_sad_string
    import io
    out = io.StringIO()
    ffs, _ = run_sad_string(src + '''
FFS;
USE RING;
SeedRandom[125188771];
   ae = 5E-5;   fe = 5E-4;   te = 1E-4;
    FFS["DELX   ae Q*;"];    FFS["DELY   ae Q*;"];
    FFS["DELX   ae S*;"];    FFS["DELY   ae S*;"];
    FFS["DK     fe B*;"];    FFS["DK     fe Q*;"];    FFS["DK     fe S*;"];
    FFS["DTHETA te B*;"];    FFS["DTHETA te Q*;"];    FFS["DTHETA te S*;"];
''', out=out)
    return ffs


def test_esrf_turn_and_emittance_match_python(esrf_err):
    ffs = esrf_err
    els = ffs.lat.elements
    cod0 = np.array([1e-4, 0.0, -5e-5, 0.0, 0.0, 1e-3])
    E._USE_ELKERN = False
    try:
        a = E.tturne(els, cod0, PH, rfsw=True, radcod=True)
    finally:
        E._USE_ELKERN = True
    b = E.tturne(els, cod0, PH, rfsw=True, radcod=True)
    np.testing.assert_allclose(b.cod, a.cod, rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(b.T, a.T, rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(b.B, a.B, rtol=1e-7, atol=1e-20)
    assert abs(b.u0 - a.u0) <= 1e-10 * abs(a.u0)
    # Emittance (tcod の反復込み)
    E._USE_ELKERN = False
    _lat._USE_ELKERN = False
    try:
        ea = E.emittance6d(ffs.lat, 3e9, rfsw=True, radcod=True)
    finally:
        E._USE_ELKERN = True
        _lat._USE_ELKERN = True
    eb = E.emittance6d(ffs.lat, 3e9, rfsw=True, radcod=True)
    for k in ('emittance_x', 'emittance_y', 'momentum_spread', 'energy_loss_eV'):
        assert abs(eb[k] - ea[k]) <= 1e-6 * abs(ea[k]), (k, ea[k], eb[k])


# ---------------------------------------------------------------------------
# 4D 1 周カーネル (turnkern.turn4 / lattice._one_turn_cod_matrix) と要素ごとの Python 版の一致
# ---------------------------------------------------------------------------

def _turn(els, cod0, delta, kern):
    _lat._USE_TURNKERN = kern
    _lat._EL_CACHE.clear()
    try:
        M, c, em = _lat._one_turn_cod_matrix(els, delta, cod0)
    finally:
        _lat._USE_TURNKERN = True
    return M, c, np.array([m[0] for m in em]), np.array([m[1] for m in em])


def _line_with_fallbacks():
    """要素 11 種 + MARK/MONI + MULT (Python に戻る) + K1=0 の QUAD (解析ヤコビアン不可 → Python)"""
    from pysad.elements import Marker, Monitor, Multipole
    els = [el for _, el in _els()]
    m = Multipole(name='M', L=0.1); m.KnL = {2: 0.3, 3: 1.0}
    q0 = Quadrupole(name='Q0', L=0.2, K1=0.0)
    return els[:3] + [Marker(name='MK'), m] + els[3:7] + [Monitor(name='MO'), q0] + els[7:]


@pytest.mark.parametrize('delta', [0.0, 1e-3])
def test_turn_kernel_matches_python_line(delta):
    pytest.importorskip('numba')
    els = _line_with_fallbacks()
    cods = list(_cods()) + [np.zeros(6)]            # 0 は設計軌道の近道 (T0) を通す
    for cod in cods:
        cod = cod.copy()
        cod[5] = delta
        Ma, ca, ma, cda = _turn(els, cod, delta, False)
        Mb, cb, mb, cdb = _turn(els, cod, delta, True)
        np.testing.assert_array_equal(mb, ma)         # 要素行列と入口 cod は同一 (同じカーネル・同じ式)
        np.testing.assert_array_equal(cdb, cda)
        np.testing.assert_array_equal(cb, ca)
        np.testing.assert_allclose(Mb, Ma, rtol=1e-12, atol=1e-14)   # 一周の積だけ積の順序 (numpy の BLAS と numba) で 1e-13
    # 要素を書き換えたら pack が詰め直される
    els[5].K1 *= 1.01
    Ma, ca, ma, cda = _turn(els, cods[0], delta, False)
    Mb, cb, mb, cdb = _turn(els, cods[0], delta, True)
    np.testing.assert_array_equal(mb, ma)
    np.testing.assert_array_equal(cb, ca)


def test_turn_kernel_matches_python_esrf(esrf_err):
    els = esrf_err.lat.elements
    for cod0, delta in [(np.zeros(6), 0.0), (np.array([1e-4, 0.0, -5e-5, 0.0, 0.0, 0.0]), 0.0),
                        (np.array([1e-4, 2e-5, -5e-5, 1e-5, 0.0, 1e-3]), 1e-3)]:
        Ma, ca, ma, cda = _turn(els, cod0, delta, False)
        Mb, cb, mb, cdb = _turn(els, cod0, delta, True)
        np.testing.assert_array_equal(mb, ma)
        np.testing.assert_array_equal(cdb, cda)
        np.testing.assert_allclose(Mb, Ma, rtol=1e-11, atol=1e-13)
    # qcod (Newton 反復込み) の閉軌道
    _lat._USE_TURNKERN = False
    _lat._EL_CACHE.clear()
    try:
        ca, Ma, _ = _lat._qcod_periodic(els, 0.0)
    finally:
        _lat._USE_TURNKERN = True
    cb, Mb, _ = _lat._qcod_periodic(els, 0.0)
    np.testing.assert_allclose(cb, ca, rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(Mb, Ma, rtol=1e-9, atol=1e-12)
