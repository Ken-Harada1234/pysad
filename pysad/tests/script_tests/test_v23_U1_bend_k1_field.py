# -*- coding: utf-8 -*-
"""v2.3.0 U1: K1≠0 の BEND (複合機能) の K1 の定義 (B160 の切り分け)。

独立な数値積分 (_u1_indep.py: 実験室系のローレンツ力、Bρ0 = 1) と、gs の保存値 (example/15_lma_touschek/touschek_elegant/gs_bend_k1.out.txt の 5 つの δ のチューンの端数) の比較。
磁場の模型は K1 = g'(ρ0)/Bρ0 が共通で、g''(ρ0) だけが違う: 'lin' (B の勾配が x に一様、elegant K2 = 0)、'log' (r ∂B/∂r が一様)、'pow' (磁場指数が一様)。
結論: K64 の BEND の K1 の本体 (硬い端・フリンジなし) は 'log' の磁場の厳密な運動と 1e-10 で一致する。'lin' とは一致しない (1e-5〜1e-4 違う)。
pysad の挙動は変えない (試験は元の挙動の確認だけ)。
"""
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
K1 = -NIDX / RHO**2          # per length [1/m^2]
H = 1.0 / RHO

# gs (K64) の追跡のチューンの端数 |f| (x, y)。example/15_lma_touschek/touschek_elegant/gs_bend_k1.out.txt (400 周、1e-6 m の振動の NAFF)。ν = 1 − |f|
GS_DELTAS = [-0.004, -0.002, 0.0, 0.002, 0.004]
GS_FX = [0.2247161719, 0.2250594092, 0.2254033310, 0.2257479407, 0.2260932435]
GS_FY = [0.3683869964, 0.3679659188, 0.3675444692, 0.3671226428, 0.3667004387]


@pytest.fixture(autouse=True)
def _mode():
    old = compat._mode
    compat.set_mode('k64')
    yield
    compat.set_mode(old)


def test_log_field_matches_gs_tunes():
    """'log' の磁場の厳密な運動のチューンが gs のチューンと 5 つの δ で 1e-7 以内 (どの δ でも)。'lin' は合わない。"""
    for d, fx, fy in zip(GS_DELTAS, GS_FX, GS_FY):
        nx, ny, _ = U.tunes('log', RHO, K1, N, d)
        assert abs(nx - (1 - fx)) < 1e-7 and abs(ny - (1 - fy)) < 1e-7
    nx, ny, _ = U.tunes('lin', RHO, K1, N, 0.004)
    assert abs(nx - (1 - GS_FX[-1])) > 1e-3 and abs(ny - (1 - GS_FY[-1])) > 1e-3


def test_dnu_ddelta_models():
    """dν/dδ: 'log' = gs の (−0.17213, +0.21082)、'lin' = 解析式の (−0.60246, +0.73786)、'pow' (磁場指数が一様) = 0。"""
    cx, cy, _, _ = U.dnu_ddelta('log', RHO, K1, N)
    assert abs(cx + 0.17213) < 2e-5 and abs(cy - 0.21082) < 2e-5
    cx, cy, _, _ = U.dnu_ddelta('lin', RHO, K1, N)
    assert abs(cx + 0.60246) < 1e-4 and abs(cy - 0.73786) < 1e-4
    cx, cy, _, _ = U.dnu_ddelta('pow', RHO, K1, N)
    assert abs(cx) < 1e-5 and abs(cy) < 1e-5


def test_cubic_terms_second_difference():
    """1 台の薄い扇形磁石 (中心角 2 mrad) の px の x²・y² の係数 (2 階差分)。
    'lin' は K1 h の (x³/3 − x y²/2) の項 (px の x² の係数 −hK1 L、y² の係数 +hK1 L/2 を足す)、'log' は (x³/6) だけ (x² は −hK1 L/2、y² は 0)。"""
    th = 2e-3
    L = RHO * th
    a = 2e-3
    c = {}
    for m in ('lin', 'log'):
        f = U.make_field(m, RHO, K1)

        def px(x, y):
            return U.sector_map(f, RHO, th, 1.0, [x, 0, y, 0])[1]
        c[m] = ((px(a, 0) + px(-a, 0) - 2 * px(0, 0)) / (2 * a * a), (px(0, a) + px(0, -a) - 2 * px(0, 0)) / (2 * a * a))
    hk = H * K1 * L
    # 運動学 (曲率) の寄与は両方に共通。差を見る
    assert abs((c['lin'][0] - c['log'][0]) - (-hk / 2)) < 2e-9
    assert abs((c['lin'][1] - c['log'][1]) - (+hk / 2)) < 2e-9
    assert abs(c['log'][1]) < 1e-9                    # 'log': x y² の項なし (K64 の取り出した生成関数と同じ)
    assert abs(c['log'][0] - (-hk / 2)) < 2e-9        # 'log': x³ の項は K1 h/6 → px の x² の係数 −hK1 L/2
    assert abs(c['lin'][0] - (-hk)) < 2e-9            # 'lin': −hK1 L


def _ring_int(model, s0, p):
    f = U.make_field(model, RHO, K1)
    s = np.array(s0, float)
    for _ in range(N):
        s = U.sector_map(f, RHO, ANG, p, s)
    return s


@pytest.mark.parametrize('delta', [0.0, 0.004, -0.004])
def test_pysad_body_equals_log_field_nonlinear(delta):
    """pysad (K64 と同じ BEND の本体。DISFRIN=1 でフリンジを切る) の 16 台の追跡の写像が 'log' の磁場の厳密な積分と 1e-9 で一致し、
    'lin' とは 1e-6 以上違う (振幅 5 mm、δ = 0、±0.4%)。追跡の px・py は P0 規格 (方向余弦 × (1 + δ))。"""
    from pysad.trksetup import TrackSetup
    src = ('MOMENTUM = 3 GEV;\nMARK M0 = ();\nBEND B = (L = %g ANGLE = %.12f K1 = %.12f*%g DISFRIN=1);\nLINE RING = (M0 %d*B);\n'
           % (LB, ANG, K1, LB, N))
    st = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=1.0)[0]
    p = 1.0 + delta
    xc = U.closed_orbit_radius(U.make_field('log', RHO, K1), RHO, p) - RHO
    amp = 5e-3
    s0 = np.array([xc + amp, 0.5 * amp, 1.5 * amp, -0.3 * amp])
    c = np.array([s0[0], s0[1] * p, s0[2], s0[3] * p, 0.0, delta])
    r = st.track(1, c.reshape(1, 6), 1, record=True)['rec'][0, 0]
    ps = np.array([r[0], r[1] / p, r[2], r[3] / p])
    assert np.max(np.abs(_ring_int('log', s0, p) - ps)) < 1e-9
    assert np.max(np.abs(_ring_int('lin', s0, p) - ps)) > 1e-6


def test_mult_angle_k1_same_as_bend_k1():
    """SADHelp の multipole_with_nonzero_ANGLE (A_s の展開) で定義される MULT (ANGLE≠0) の K1 は、BEND の K1 の本体と同じ dν/dδ (1e-9)。"""
    from pysad import resdrive
    from pysad.trksetup import TrackSetup
    out = {}
    for el in ('BEND', 'MULT'):
        src = ('MOMENTUM = 3 GEV;\nMARK M0 = ();\n%s B = (L = %g ANGLE = %.12f K1 = %.12f*%g);\nLINE RING = (M0 %d*B);\n'
               % (el, LB, ANG, K1, LB, N))
        _, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=1.0)
        out[el] = np.array(resdrive.chromaticity_numeric(ffs.lat)['d1'])
    assert np.allclose(out['BEND'], out['MULT'], atol=1e-9)
    assert abs(out['BEND'][0] + 0.17213) < 2e-5 and abs(out['BEND'][1] - 0.21082) < 2e-5
