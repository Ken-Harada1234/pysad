# -*- coding: utf-8 -*-
"""v2.3.0 U2: K1≠0 の BEND の K1 の定義のオプション (compat.set_bend_k1_def('log'|'lin')、B160) と B161 (CALC の dk の修正)。
基準は独立な数値積分 (_u1_indep.py: 実験室系のローレンツ力)。n = 0.4 の扇形磁石 16 台のリング (3 GeV 相当、_u1_indep と同じ)。"""
import math
import os
import sys

import numpy as np
import pytest

from pysad import compat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _u1_indep as U  # noqa: E402

N, LB, NIDX = 16, 1.0, 0.4
ANG = 2 * math.pi / N
RHO = LB / ANG
K1 = -NIDX / RHO**2


@pytest.fixture(autouse=True)
def _reset():
    old = (compat._mode, compat._bend_k1_def, compat._bend_dk_fix)
    compat.set_mode('v2')
    compat.set_bend_k1_def(None)
    compat.set_bend_dk_fix(None)
    yield
    compat._mode, compat._bend_k1_def, compat._bend_dk_fix = old


def _src(el='BEND', extra=''):
    return ('MOMENTUM = 3 GEV;\nMARK M0 = ();\n%s B = (L = %g ANGLE = %.12f K1 = %.12f*%g %s);\nLINE RING = (M0 %d*B);\n'
            % (el, LB, ANG, K1, LB, extra, N))


def _setup(el='BEND', extra='', rad=False):
    from pysad.trksetup import TrackSetup
    return TrackSetup.from_source(_src(el, extra), 'RING', rad=rad, rfsw=False, alost=1.0)


def _frac(v):
    return min(v % 1, 1 - v % 1)


def test_default_and_option_api(monkeypatch):
    assert compat.bend_k1_def() == 'log'
    compat.set_bend_k1_def('lin')
    assert compat.bend_k1_def() == 'lin' and compat.bend_k1_lin()
    compat.set_mode('k64')
    assert compat.bend_k1_def() == 'log' and not compat.bend_k1_lin() and compat.bend_ecorr_flags() == 0   # k64 では無視
    compat.set_mode('v2')
    compat.set_bend_k1_def(None)
    monkeypatch.setenv('PYSAD_BEND_K1_DEF', 'lin')
    assert compat.bend_k1_def() == 'lin'
    with pytest.raises(ValueError):
        compat.set_bend_k1_def('pow')
    assert compat.bend_ecorr_flags() == 7   # bit0,1 (v2.3.0) + bit2 (B250, v2.4.0)


def test_calc_dnu_ddelta_lin_and_log():
    """CALC (解析行列 + 閉軌道) の dν/dδ: 'log' = K64 (−0.17213、+0.21082)、'lin' = B' 一様の解析式・elegant (−0.60246、+0.73786)・独立な積分。"""
    from pysad import resdrive
    _, ffs = _setup()
    r = resdrive.chromaticity_numeric(ffs.lat)['d1']
    assert abs(r[0] + 0.17213) < 2e-5 and abs(r[1] - 0.21082) < 2e-5
    compat.set_bend_k1_def('lin')
    r = resdrive.chromaticity_numeric(ffs.lat)['d1']
    assert abs(r[0] + 0.60246) < 2e-5 and abs(r[1] - 0.73786) < 2e-5
    cx, cy, _, _ = U.dnu_ddelta('lin', RHO, K1, N)
    assert abs(r[0] - cx) < 1e-4 and abs(r[1] - cy) < 1e-4


@pytest.mark.parametrize('model', ['log', 'lin'])
def test_tracking_tunes_vs_integral_5_deltas(model):
    """追跡の 1 周の写像 (閉軌道の周りの線形化) と CALC のチューンが、独立な積分と 5 つの δ で 1e-7 以内 (DISFRIN=1)。"""
    compat.set_bend_k1_def(model)
    st, ffs = _setup(extra='DISFRIN=1')
    for delta in (-0.004, -0.002, 0.0, 0.002, 0.004):
        p = 1 + delta
        xc = U.closed_orbit_radius(U.make_field(model, RHO, K1), RHO, p) - RHO
        s0 = np.array([xc, 0, 0, 0, 0, delta], float)

        def f(c):
            cc = np.array(c, float)
            cc[1] *= p
            cc[3] *= p
            r = st.track(1, cc.reshape(1, 6), 1, record=True)['rec'][0, 0]
            return np.array([r[0], r[1] / p, r[2], r[3] / p])
        m = np.zeros((4, 4))
        for j in range(4):
            d = np.zeros(6)
            d[j] = 1e-7
            m[:, j] = (f(s0 + d) - f(s0 - d)) / 2e-7
        tx = math.acos(0.5 * (m[0, 0] + m[1, 1])) / 2 / math.pi
        ty = math.acos(0.5 * (m[2, 2] + m[3, 3])) / 2 / math.pi
        ix, iy, _ = U.tunes(model, RHO, K1, N, delta)
        assert abs(tx - _frac(ix)) < 1e-7 and abs(ty - _frac(iy)) < 1e-7
        t = ffs.lat.twiss_table(delta=delta)[-1]
        assert abs(_frac(t.mux / 2 / math.pi) - _frac(ix)) < 1e-7 and abs(_frac(t.muy / 2 / math.pi) - _frac(iy)) < 1e-7


@pytest.mark.parametrize('delta', [0.0, 0.004, -0.004])
def test_tracking_map_lin_vs_integral(delta):
    """'lin' の 16 台の写像が B' 一様の独立な積分と一致 (DISFRIN=1)。x 1 mm・y 0.15 mm で 5e-10 以内、x 5 mm・y 7.5 mm で 1e-7 以内
    (kick は 3 次のポテンシャルまで。B' 一様の場の 4 次の項 K1 h² y³/6 は入れていない)。'log' の追跡とは 1e-6 以上違う。"""
    compat.set_bend_k1_def('lin')
    st = _setup(extra='DISFRIN=1')[0]
    p = 1 + delta
    xc = U.closed_orbit_radius(U.make_field('lin', RHO, K1), RHO, p) - RHO
    f = U.make_field('lin', RHO, K1)

    def trk(stt, s0):
        c = np.array([s0[0], s0[1] * p, s0[2], s0[3] * p, 0.0, delta])
        r = stt.track(1, c.reshape(1, 6), 1, record=True)['rec'][0, 0]
        return np.array([r[0], r[1] / p, r[2], r[3] / p])

    def integ(s0):
        s = np.array(s0, float)
        for _ in range(N):
            s = U.sector_map(f, RHO, ANG, p, s)
        return s
    s_small = np.array([xc + 1e-3, 0.5e-3, 1.5e-4, -0.3e-3])
    s_big = np.array([xc + 5e-3, 2.5e-3, 7.5e-3, -1.5e-3])
    assert np.max(np.abs(trk(st, s_small) - integ(s_small))) < 5e-10
    assert np.max(np.abs(trk(st, s_big) - integ(s_big))) < 1e-7
    compat.set_bend_k1_def('log')
    st2 = _setup(extra='DISFRIN=1')[0]
    assert np.max(np.abs(trk(st2, s_big) - integ(s_big))) > 1e-6


def test_k64_ignores_option_bitwise():
    """k64 モードでは 'lin' を設定しても追跡・CALC は K64 のまま (ビット単位)。"""
    compat.set_mode('k64')
    c = np.array([[3e-3, 1e-3, 2e-3, -5e-4, 0.0, 3e-3]])
    out = []
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st, ffs = _setup()
        r = st.track(3, c.copy(), 1, record=True)['rec'].copy()
        t = ffs.lat.twiss_table(delta=0.003)[-1]
        out.append((r, t.mux, t.muy))
    assert np.array_equal(out[0][0], out[1][0]) and out[0][1] == out[1][1] and out[0][2] == out[1][2]


def test_delta0_design_unchanged():
    """δ = 0 の設計軌道では 'lin' でも光学・行列が変わらない (E1・E2≠0 の長方形、DISFRIN=0 も)。δ≠0・振幅ありでは変わる。"""
    from pysad.trksetup import TrackSetup
    src = ('MOMENTUM = 3 GEV;\nMARK M0 = ();\nQUAD Q = (L = 0.3 K1 = 0.5*0.3);\nDRIFT D = (L = 0.8);\n'
           'BEND B = (L = 1 ANGLE = %.12f K1 = -0.07 E1 = 0.5 E2 = 0.5);\nLINE C = (Q D B D);\nLINE RING = (M0 %d*C);\n' % (ANG, N))
    res = {}
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=1.0)
        t = ffs.lat.twiss_table(delta=0.0)
        trk = st.track(2, np.zeros((1, 6)), 1, record=True)['rec'][0, 0].copy()
        off = st.track(2, np.array([[2e-3, 0, 1e-3, 0, 0, 3e-3]]), 1, record=True)['rec'][0, 0].copy()
        res[d] = ((t[-1].mux, t[-1].muy, max(r.betx for r in t), max(r.etx for r in t)), trk, off)
    assert res['log'][0] == res['lin'][0]
    assert np.array_equal(res['log'][1], res['lin'][1])
    assert np.max(np.abs(res['log'][2] - res['lin'][2])) > 1e-8


def test_emittance_delta0_unchanged():
    """放射・RF つきの Emittance[] が 'lin' でもほぼ変わらない (相対差 5e-3 以内。EMITX 1e-4、SIGE 2e-3 程度)。放射の損失で閉軌道が δ ≠ 0 (のこぎり状、1e-3 程度) になるため、'lin' の 6 極成分のずれが 1e-4 程度の差を作る (厳密に同じではない)。"""
    from pysad.trksetup import TrackSetup
    src = ('MOMENTUM = 3 GEV;\nMARK M0 = ();\nQUAD QF = (L = 0.3 K1 = 1.2*0.3);\nQUAD QD = (L = 0.3 K1 = -1.2*0.3);\nDRIFT D = (L = 0.6);\n'
           'BEND B = (L = 1 ANGLE = %.12f K1 = -0.06);\nCAVI C = (L = 0 VOLT = 2e6 FREQ = 5e8);\nLINE CELL = (QF D B D QD D B D);\n'
           'LINE RING = (%d*CELL C);\n' % (math.pi / 16, 16))
    out = {}
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=1.0)
        g = ffs.I.globals
        out[d] = [float(g[k]) for k in ('EMITX', 'EMITY', 'SIGE')]
    for i in (0, 2):           # EMITY は平面リングで 0 (1e-42) なので比べない
        assert abs(out['log'][i] / out['lin'][i] - 1) < 5e-3


def test_bendecorr_jacobian_lin_and_dk_fix():
    """解析行列 t (perfkern_bend.bendecorr) が出口 cod の数値微分と一致 ('lin' も)。
    B161: dk を直す (flags bit1) と x/ρ0 = 0.1 でも 1e-8 以内、K64 の dk (flags 0) は 1e-5 以上ずれる。"""
    from pysad.perfkern_bend import bendecorr

    def jac(cod, ak1, al, rho0, flags):
        J = np.zeros((6, 6))
        for j in range(6):
            d = np.zeros(6)
            d[j] = 1e-7
            J[:, j] = (bendecorr(cod + d, ak1, al, rho0, 0.9999, flags)[1]
                       - bendecorr(cod - d, ak1, al, rho0, 0.9999, flags)[1]) / 2e-7
        return J
    rho0, al = 2.5465, 0.01
    ak1 = -0.0617 * al
    cod = np.array([3e-3, 1e-3, 2e-3, -5e-4, 0.0, 3e-3])
    for fl in (1, 3):
        t = bendecorr(cod, ak1, al, rho0, 0.9999, fl)[0]
        assert np.max(np.abs(t - jac(cod, ak1, al, rho0, fl))) < 1e-9
    cod2 = np.array([0.1 * rho0, 1e-3, 5e-3, -5e-4, 0.0, 3e-3])
    ak2 = -0.6 * al
    j = jac(cod2, ak2, al, rho0, 2)
    assert np.max(np.abs(bendecorr(cod2, ak2, al, rho0, 0.9999, 2)[0] - j)) < 1e-8
    assert np.max(np.abs(bendecorr(cod2, ak2, al, rho0, 0.9999, 0)[0] - j)) > 1e-5


def test_b161_ring_calc_vs_integral():
    """B161: K1 の強い扇形磁石のリングで、v2 (dk を直した CALC) のチューンの δ 依存が独立な積分 ('log') と δ = 2% で 2e-7 以内。
    K64 の dk (set_bend_dk_fix(False)・k64) は 1e-4 ずれる。"""
    _, ffs = _setup()
    d = 0.02
    ix, iy, _ = U.tunes('log', RHO, K1, N, d)
    assert abs(ffs.lat.twiss_table(delta=d)[-1].mux / 2 / math.pi - ix) < 2e-7
    compat.set_bend_dk_fix(False)
    assert abs(ffs.lat.twiss_table(delta=d)[-1].mux / 2 / math.pi - ix) > 5e-5
    compat.set_bend_dk_fix(None)
    compat.set_mode('k64')
    assert abs(ffs.lat.twiss_table(delta=d)[-1].mux / 2 / math.pi - ix) > 5e-5


def test_mult_angle_k1_follows_option():
    """MULT (ANGLE≠0) の K1 も同じ本体: 'lin' でも BEND と同じ dν/dδ (1e-9)。"""
    from pysad import resdrive
    out = {}
    for d in ('log', 'lin'):
        compat.set_bend_k1_def(d)
        for el in ('BEND', 'MULT'):
            _, ffs = _setup(el)
            out[(d, el)] = np.array(resdrive.chromaticity_numeric(ffs.lat)['d1'])
        assert np.allclose(out[(d, 'BEND')], out[(d, 'MULT')], atol=1e-9)
    assert abs(out[('lin', 'MULT')][0] + 0.60246) < 2e-5
    assert abs(out[('log', 'MULT')][0] + 0.17213) < 2e-5


def test_python_tracker_equals_numba_lin():
    """Python 版の追跡 (track6._tbendi) と numba (track6n.k_bendi) が 'lin' でも 1e-12 以内で同じ。"""
    from pysad import track6 as t6
    compat.set_bend_k1_def('lin')
    st, _ = _setup(extra='DISFRIN=1')
    c = np.array([[3e-3, 1e-3, 2e-3, -5e-4, 0.0, 3e-3]])
    rn = np.array(st.track(1, c.copy(), 1, record=True)['rec'][0, 0])
    tr = np.array(t6.track_particles(st.els, c.copy(), 1, st.p0, 0.0, False, False, False))
    assert np.max(np.abs(rn[:4] - tr.reshape(-1)[:4])) < 1e-12
