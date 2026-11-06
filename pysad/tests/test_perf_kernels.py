# -*- coding: utf-8 -*-
"""性能改善 (11/05, docs/history/work_perf2.md) の回帰テスト:
(1) K1≠0 BEND の tbendecorr/tbendebody の numba 版 (perfkern_bend) が numba 化前の Python 版の参照値 (data_bend_kernels.npz) と一致する
(2) tchge (BEND の据え付け誤差) の解析ヤコビアン (_tbrote_jac/_tbgrote_mat/_tchge_jac) が写像の数値微分と一致する
    ((5,6) は SAD の式のままなので除く。入口で th≠0 のときの「回転が入らない」癖は _tbrote_jac と同じ行列になることを見る)
(3) _State.apply の in-place 版が numpy の積と一致する
"""
import os
import itertools

import numpy as np
import pytest

from pysad import lattice as _lat
from pysad import emit6d as E
from pysad import perfkern_bend as K

HERE = os.path.dirname(os.path.abspath(__file__))
PH = E._Phys(3e9, 603.12)


def test_bend_kernels_match_reference():
    d = np.load(os.path.join(HERE, 'data_bend_kernels.npz'))
    cods, par = d['cods'], d['par']
    assert abs(float(d['h0']) - PH.h0) < 1e-6 * PH.h0
    for i in range(len(cods)):
        al, phi, ak1, rho0, rhob = par[i]
        st = E._State(cods[i], PH)
        E._tbendecorr(st, ak1 * 0.5, al * 0.5, rho0)
        np.testing.assert_allclose(st.T, d['tc'][i], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(st.cod, d['cc'][i], rtol=1e-12, atol=1e-15)
        st = E._State(cods[i], PH)
        E._tbendebody(st, al, phi, ak1, al, rhob, rho0, False)
        np.testing.assert_allclose(st.T, d['tb'][i], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(st.cod, d['cb'][i], rtol=1e-12, atol=1e-15)


def _rand_cods(n=3, seed=0):
    rng = np.random.default_rng(seed)
    c = rng.normal(size=(n, 6)) * np.array([2e-3, 1e-3, 2e-3, 1e-3, 1e-2, 2e-3])
    c[0] = [1e-3, -2e-3, 2e-3, 1e-3, 0.01, 3e-3]
    return c


def _num(f, c, eps=1e-6):
    return _lat._num_jacobian(f, c, eps)


def _cmp_rows(Ja, Jn, tol=2e-7):
    """x,px,y,py,z の行 (z の (5,6) を除く) を比べる"""
    m = np.ones((6, 6), dtype=bool)
    m[4, 5] = False
    m[5, :] = False
    d = np.abs(Ja - Jn) / np.maximum(np.abs(Jn), 1.0)
    assert d[m].max() < tol, (d[m].max(), np.unravel_index(np.argmax(d * m), d.shape))


@pytest.mark.parametrize('ent', [True, False])
@pytest.mark.parametrize('phig,alg', [(0.0, 0.3), (0.05, 0.4)])
@pytest.mark.parametrize('dx,dz,dtheta', [(1e-3, 0.0, 0.0), (0.0, 0.0, 2e-4), (1e-3, 5e-4, 3e-4), (0.0, 2e-3, 0.0)])
def test_tbrote_jac_vs_numeric(ent, phig, alg, dx, dz, dtheta):
    for c in _rand_cods():
        f = lambda cc: _lat._tbrote_cod(cc, alg, phig, dx, dz, dtheta, 0.0, ent)
        cod_out, J = _lat._tbrote_jac(c, alg, phig, dx, dz, dtheta, 0.0, ent)
        np.testing.assert_allclose(cod_out, f(c), rtol=1e-13, atol=1e-18)
        _cmp_rows(J, _num(f, c))


@pytest.mark.parametrize('ent', [True, False])
@pytest.mark.parametrize('th', [0.0, 0.2])
@pytest.mark.parametrize('phig,alg', [(0.0, 0.0), (0.05, 0.4)])
@pytest.mark.parametrize('dx0,dy,dz,dtheta', [(1e-3, -2e-3, 0.0, 0.0), (0.0, 0.0, 0.0, 2e-4), (1e-3, 2e-3, 5e-4, 3e-4), (0.0, 1e-3, 0.0, 0.0)])
def test_tchge_jac(ent, th, phig, alg, dx0, dy, dz, dtheta):
    for c in _rand_cods():
        f = lambda cc: _lat._tchge_cod(cc, dx0, dy, dz, th, dtheta, alg, phig, ent)
        J = _lat._tchge_jac(c, dx0, dy, dz, th, dtheta, alg, phig, ent)
        thr = th if phig != 0.0 else th + dtheta
        ct, st = np.cos(thr), np.sin(thr)
        dxr = ct * dx0 - st * dy
        brote = (dxr != 0.0 or dz != 0.0 or dtheta != 0.0)
        if ent and thr != 0.0 and brote:
            # SAD の癖: 回して y をずらした後の cod での tbrote の行列だけ (回転は入らない)
            c2 = c.copy()
            c2[0] = ct * c[0] - st * c[2]; c2[2] = st * c[0] + ct * c[2]
            c2[1] = ct * c[1] - st * c[3]; c2[3] = st * c[1] + ct * c[3]
            c2[2] -= st * dx0 + ct * dy
            fb = lambda cc: _lat._tbrote_cod(cc, alg, phig, dxr, dz, dtheta, 0.0, True)
            _cmp_rows(J, _num(fb, c2))
        else:
            _cmp_rows(J, _num(f, c))


def test_apply_inplace_matches_numpy():
    rng = np.random.default_rng(3)
    for _ in range(5):
        J = rng.normal(size=(6, 6))
        T = rng.normal(size=(6, 6)); D = rng.normal(size=(6, 6)); B = rng.normal(size=(6, 6)); B = B + B.T
        T2, D2, B2 = T.copy(), D.copy(), B.copy()
        K.apply_tdb(J, T2, D2, B2, True)
        np.testing.assert_allclose(T2, J @ T, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(D2, J @ D, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(B2, J @ B @ J.T, rtol=1e-13, atol=1e-14)
        T3, D3, B3 = T.copy(), D.copy(), B.copy()
        K.apply_tdb(J, T3, D3, B3, False)
        np.testing.assert_allclose(T3, J @ T, rtol=1e-13, atol=1e-14)
        assert np.array_equal(D3, D) and np.array_equal(B3, B)
