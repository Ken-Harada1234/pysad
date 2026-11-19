# -*- coding: utf-8 -*-
"""v2next 担当 P3 項目 1: TrackParticles の F1 (ANGLE=0 の BEND の `tsteer`、MULT の `tmulti` の F1・F2・FB1・FB2・FRINGE) と、MULT の端の角 AE1・AE2 の追跡。

(a) K64 の動作 (k64 モード): gs (K64 の 10-19 ビルド、2026-11-19) の追跡 (NORAD・RAD、3 周) と 10 桁表示の全桁で一致 (12 セルの FODO、3 GeV、周長 45.6 m。
    ステアリング 6 通り、MULT 9 通り。大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m) の先頭に補正磁石 2 台を入れたもの 3 通り)。
    F1 の追跡は K64 の動作の追加 (gate なし)。高速追跡 (numba の `k_steer`) と Python (track6) は 1e-13 で一致 (RNG の最後の状態も)。
(b) 独立な基準 (v2・k64 とも同じ値): 軸上の偏向角 = K0/(1+δ) (F1 に依らない)、放射の積分 ∫B²dl (磁場が直線で立ち上がる端: 1 − (fb1+fb2)/(6 al))、
    端の磁場の鉛直収束 (x 方向に一様な Maxwell 場 (By = B − y²B''/2、Bs = yB') の中の厳密な運動、正準方程式の数値積分: Δpy/y = (fb1+fb2)/(6ρ²))、
    x・px は fb1 = fb2 のとき厳密な運動と 1e-10。追跡の一周写像の固有値 (閉軌道のまわり) = CALC のチューン (別の経路: tsteee・tmulte)。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat
from . import _p3_gsdata as G
from . import _p3_indep as X
from . import test_v2_E_steer as T


def _vals(out, tag):
    return np.array(X.tag_values(out, tag))


@pytest.mark.parametrize('case', list(X.TRACK_CASES))
def test_tracking_matches_gs_k64(k64_mode, case):
    """F1・FB1・FB2・FRINGE (−2..3)・AE1・AE2 つきのステアリング・MULT の追跡 (NORAD・RAD) が gs と 10 桁の全桁で一致"""
    kind, a, b = X.TRACK_CASES[case]
    out = X.run(X.ring_script(kind, a, b, X.TRACK_TAIL))[1]
    n_ref, r_ref = G.TRACK[case]
    assert np.allclose(_vals(out, 'N3'), n_ref, rtol=1e-8, atol=1e-17), case
    assert np.allclose(_vals(out, 'R3'), r_ref, rtol=1e-8, atol=1e-17), case


def _cmp_emit(d, ref):
    for k in ('calc', 'tunes'):
        assert np.allclose(d[k], ref[k], rtol=0, atol=3e-9), (k, d[k], ref[k])
    for k in ('damp', 'emit'):
        assert np.allclose(d[k], ref[k], rtol=1e-7, atol=0), (k, d[k], ref[k])
    for k in ('u0', 'alpha', 'eqpos'):
        assert abs(d[k] / ref[k] - 1.0) < 1e-7, (k, d[k], ref[k])


@pytest.mark.parametrize('key', ['ra1', 'ra2', 'ra3'])
def test_large_ring_a_tracking_and_emittance_match_gs_k64(k64_mode, key):
    """大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m) の先頭に補正磁石 2 台 (K0 = ±5e-4、L = 0.2 m) を入れる。ra1: BEND に F1・FB1・FB2・FRINGE=3・AE1・AE2、
    ra2: MULT (K1 = ±0.05) に F1・F2・FB1・FB2・FRINGE=3・AE1・AE2、ra3: BEND の AE だけ (FB2 なし)。追跡 (NORAD・RAD) と CALC・Emittance[] が gs と一致"""
    out = X.run(X.ring_a_script(key, X.RING_A_TAIL))[1]
    n_ref, r_ref = G.TRACK[key]
    assert np.allclose(_vals(out, 'N3'), n_ref, rtol=1e-8, atol=1e-17)
    assert np.allclose(_vals(out, 'R3'), r_ref, rtol=1e-8, atol=1e-17)
    _cmp_emit(T._values(out), G.EMIT[key])


# ---------------------------------------------------------------- 高速追跡 (numba) と Python
def _fast_vs_python(mode, a, b, rad, fluc):
    from pysad import track6 as t6, track6n as tn, emit6d as e6
    compat.set_mode(mode)
    ffs, _ = X.run(X.ring_script('steer', a, b, 'FFS;\nUSE RING; CELL; CALC;\n'))
    els = ffs.lat.elements
    p0 = ffs.momentum()
    ph = e6._Phys(p0, sum(el.L for el in els))
    prog = tn.Prog(els, ph)
    assert any(k == tn.K_STEER for k in prog.kinds)          # フォールバックしていない (numba の k_steer を使う)
    rng = np.random.default_rng(5)
    c = np.zeros((20, 6))
    c[:, 0] = rng.normal(0, 3e-4, 20)
    c[:, 1] = rng.normal(0, 1e-4, 20)
    c[:, 2] = rng.normal(0, 2e-4, 20)
    c[:, 3] = rng.normal(0, 1e-4, 20)
    c[:, 4] = rng.normal(0, 1e-3, 20)
    c[:, 5] = rng.normal(0, 5e-4, 20)
    ra = t6.track_particles(els, c, 5, p0, trf0=1e-4, rad=rad, fluc=fluc, seed=777)
    sa = t6.RNG.seed
    rb = tn.track_particles_fast(els, c, 5, p0, trf0=1e-4, rad=rad, fluc=fluc, seed=777)
    sb = t6.RNG.seed
    return np.abs(np.asarray(ra) - np.asarray(rb)).max(), sa == sb


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('rad,fluc', [(False, False), (True, False), (True, True)])
@pytest.mark.parametrize('a,b', [
    ('K0=0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1', 'K0=-0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1'),
    ('K0=0.003 F1=0.04 FRINGE=-1 AE2=0.3', 'K0=-0.003 F1=0.05 FRINGE=-2 DISFRIN=1 AE1=-0.2')])
def test_steer_f1_fast_equals_python(mode, rad, fluc, a, b):
    """numba の k_steer (F1・FB1・FB2・FRINGE・AE1・AE2 つき) と Python の _tsteer_ens が 1e-13 で一致 (FLUC は乱数の最後の状態も)"""
    old = compat._mode
    try:
        d, same = _fast_vs_python(mode, a, b, rad, fluc)
    finally:
        compat.set_mode(old)
    assert d < 1e-13 and same, (d, same)


# ---------------------------------------------------------------- 独立な基準
_STEER_DEFS = [
    '', 'F1=0.04 FRINGE=3', 'F1=0.06 FB1=0.02 FB2=-0.03 FRINGE=3', 'F1=0.05 FRINGE=1', 'F1=0.05 FRINGE=2', 'F1=0.05 FRINGE=-1', 'F1=0.05 FRINGE=-2',
]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('extra', _STEER_DEFS)
@pytest.mark.parametrize('delta', [0.0, 0.01, -0.02])
def test_steer_on_axis_deflection_is_k0_for_any_f1(mode, extra, delta):
    """軸上の粒子の運動量の変化 (p0 で規格化した px) の合計は F1 (断片の長さと分け方) に依らず −K0 (∫B dl は磁場の立ち上がりの形にも δ にも依らない)。FRINGE ≤ 0 (F1 が効かない) も同じ"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        c = np.zeros(6)
        c[5] = delta
        r = X.one_element_track('BEND SA = (L=0.2 ANGLE=0 K0=0.003 %s);' % extra, c)
    finally:
        compat.set_mode(old)
    assert abs(r[1] + 0.003) < 1e-14 and abs(r[3]) < 1e-18, r


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('fb1,fb2', [(0.04, 0.04), (0.1, 0.06), (0.05, 0.0001)])
def test_steer_radiated_energy_equals_ramp_integral(mode, fb1, fb2):
    """放射の積分 ∫B² dl: 磁場が長さ fb で直線に 0 → B0 と立ち上がると ∫ = B0² (al − (fb1+fb2)/2 + (fb1+fb2)/3) = B0² al (1 − (fb1+fb2)/(6 al))。
    F1 なしの放射量との比が 1e-5 で一致 (放射の刻み alr = f1r の 2 つの断片の B² の積分 (rbl² + rbh²) f1r = 2f1r/3 の確認)。高速追跡と Python の両方"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        al = 0.2

        def u(extra, fast):
            return -X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=0.02 %s);' % (al, extra), np.zeros(6), rad=True, fast=fast)[5]
        u0 = u('', True)
        for fast in (True, False):
            ratio = u('F1=%r FB2=%r FRINGE=3' % (fb1, fb2 - fb1), fast) / u0
            assert abs(ratio - (1.0 - (fb1 + fb2) / (6.0 * al))) < 2e-5, (fast, ratio)
    finally:
        compat.set_mode(old)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('fb1,fb2', [(0.04, 0.04), (0.06, 0.02), (0.1, 0.0001)])
def test_steer_f1_vertical_focusing_matches_maxwell_ramp(mode, fb1, fb2):
    """F1 のあるステアリングの鉛直の kick Δpy = ((fb1+fb2)/(6ρ²)) y (F1 なしとの差)。基準は x 方向に一様な Maxwell 場の端の厳密な運動
    (X.maxwell_steer: ベクトルポテンシャルの正準方程式の数値積分)。厳密な運動の py/y は −B0² al + B0² (fb1+fb2)/6 (解析式 ∫F B' ds) と 1e-3 (数値積分の誤差) で一致"""
    al, K0, y0 = 0.2, 0.003, 1e-5
    B0 = K0 / al
    ex = X.maxwell_steer(fb1, fb2, al, K0, y0)[3] / y0
    assert abs(ex / B0**2 - (-al + (fb1 + fb2) / 6.0)) < 1e-3 * (fb1 + fb2) / 6.0 + 2e-6
    old = compat._mode
    compat.set_mode(mode)
    try:
        c = np.zeros(6)
        c[2] = y0
        k = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r F1=%r FB1=0 FB2=%r FRINGE=3);' % (al, K0, fb1, fb2 - fb1), c)
        k0 = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r);' % (al, K0), c)
    finally:
        compat.set_mode(old)
    d = (k[3] - k0[3]) / y0
    assert abs(d / (B0**2 * (fb1 + fb2) / 6.0) - 1.0) < 1e-4, (d, B0**2 * (fb1 + fb2) / 6.0)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('delta', [0.0, 0.01, -0.02])
def test_steer_f1_horizontal_matches_maxwell_ramp_for_equal_fb(mode, delta):
    """fb1 = fb2 のとき、F1 のあるステアリングの出口の x・px が x 方向に一様な Maxwell 場の厳密な運動 (出口の面 s = al に戻す) と 1e-10 で一致 (δ ≠ 0 の色収差 dxfr も含む)。
    fb1 ≠ fb2 では x が (K0/al)(fb1² − fb2²)/24 だけ違う (K64 の端の模型の仕様。記録に書いた)"""
    al, K0, fb = 0.2, 0.003, 0.04
    x, px, _y, _py = X.maxwell_steer(fb, fb, al, K0, 0.0, delta)
    pz = math.sqrt((1 + delta)**2 - px * px)
    xnom = x - (fb / 2) * px / pz
    old = compat._mode
    compat.set_mode(mode)
    try:
        c = np.zeros(6)
        c[5] = delta
        k = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r F1=%r FRINGE=3);' % (al, K0, fb), c)
    finally:
        compat.set_mode(old)
    assert abs(k[0] - xnom) < 1e-10 and abs(k[1] - px) < 1e-12, (k[:2], xnom, px)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('fb1,fb2', [(0.04, 0.04), (0.06, 0.02)])
def test_mult_fb_vertical_focusing_matches_maxwell_ramp(mode, fb1, fb2):
    """MULT (K0 のみ、FRINGE=3) の FB1・FB2 (tblfri) の鉛直の kick の差も (fb1+fb2)/(6ρ²) y (ρ = L/K0)。同じ厳密な運動を基準にする"""
    al, K0, y0 = 0.2, 0.003, 1e-5
    B0 = K0 / al
    old = compat._mode
    compat.set_mode(mode)
    try:
        c = np.zeros(6)
        c[2] = y0
        k = X.one_element_track('MULT SA = (L=%r K0=%r FRINGE=3 FB1=%r FB2=%r);' % (al, K0, fb1, fb2), c)
        k0 = X.one_element_track('MULT SA = (L=%r K0=%r FRINGE=3);' % (al, K0), c)
    finally:
        compat.set_mode(old)
    d = (k[3] - k0[3]) / y0
    assert abs(d / (B0**2 * (fb1 + fb2) / 6.0) - 1.0) < 1e-4, (d, B0**2 * (fb1 + fb2) / 6.0)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_steer_fringe_zero_ignores_f1(mode):
    """FRINGE=0 (既定) の F1 は tpara が fb1・fb2 を 0 にするので追跡に効かない (gs と同じ: b_e)。MULT の F1・FB は FRINGE=1・2・3 だけ"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        c = np.zeros(6)
        c[0], c[2], c[5] = 1e-3, 1e-3, 0.01
        a = X.one_element_track('BEND SA = (L=0.2 ANGLE=0 K0=0.003 F1=0.05 FB1=0.02);', c, rad=True)
        b = X.one_element_track('BEND SA = (L=0.2 ANGLE=0 K0=0.003);', c, rad=True)
        m1 = X.one_element_track('MULT SA = (L=0.2 K1=0.7 K0=0.002 F1=0.05 FB1=0.02);', c, rad=True)
        m2 = X.one_element_track('MULT SA = (L=0.2 K1=0.7 K0=0.002);', c, rad=True)
    finally:
        compat.set_mode(old)
    assert np.array_equal(a, b) and np.array_equal(m1, m2)


# ---------------------------------------------------------------- 追跡と CALC (別の経路) の一致
@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_steer_f1_tracking_tune_equals_calc(mode):
    """F1・FB1・FB2・FRINGE=3・AE1・AE2 つきのステアリング (FB2 ≠ 0 なので K64 の追跡も AE2 を入れる) の追跡の一周写像の固有値 (閉軌道のまわり)
    が CALC (tsteee) のチューンと 1e-6 で一致"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        a = 'K0=0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01'
        b = 'K0=-0.003 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.02 AE2=-0.01'
        src = X.ring_script('steer', a, b, 'FFS;\nUSE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nExit;\n')
        nu, calc = X.tracking_and_calc_tunes(src)
    finally:
        compat.set_mode(old)
    assert np.allclose(nu, calc, atol=2e-6), (nu, calc)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_mult_f1_tracking_tune_equals_calc(mode):
    """MULT (K1・K0・K2、F1・F2・FB1・FB2・FRINGE=3) の追跡の一周写像の固有値が CALC (tmulte) のチューンと 1e-6 で一致 (AE なし。AE は test_v2_P3_ae.py)"""
    old = compat._mode
    compat.set_mode(mode)
    try:
        a = 'F1=0.05 F2=0.02 FRINGE=3 K0=0.0004 FB1=0.03 FB2=0.02 K2=5'
        b = 'F1=0.04 F2=0.03 FRINGE=3 K0=-0.0004 FB1=0.03 FB2=0.02 K2=-5'
        src = X.ring_script('mult', a, b, 'FFS;\nUSE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nExit;\n')
        nu, calc = X.tracking_and_calc_tunes(src)
    finally:
        compat.set_mode(old)
    assert np.allclose(nu, calc, atol=2e-6), (nu, calc)
