# -*- coding: utf-8 -*-
"""v2next 担当 P3 項目 2・3: MULT の端の角 AE1・AE2 と B108 の残り (K0 つき MULT・K1≠0 の ANGLE=0 BEND の端の角、追跡 (tsteer・tmulti) の端の放射)。

K64 (gs、10-19 ビルド、2026-11-19) の動作 (k64 モード):
- Emittance[]・CALC は MULT の AE を `tmulke` (apsi) で入れる (pysad は以前これを読まず黙って無視していた: 追加。gs と全桁一致: test_v2_P3_f1.py の
  EMIT_CASES と大きなリング A)。出口の端の蹴りは最後の半スライスの前に置く。追跡 (tmulti の tmulk) は半スライスの後。
- 追跡の `tsteer` (ANGLE=0 の BEND) は出口の端の角 AE2 を fb2 ≠ 0 のときだけ入れる (CALC の `tsteee` は常に入れる)。
- 追跡の `tbend` は K1≠0・K0≠0 の ANGLE=0 BEND を tmulti に渡すとき端の角 psi = 0 を渡す (Emittance[]・CALC の tbende → tmulte は AE を渡す)。
  この BEND は gs が USE の中で SIGSEGV するので K64 の値は確認できない (K64 のソースの読みのみ)。
v2 の動作 (根拠は別の経路との一致と、厳密な幾何):
- 追跡の AE2 は常に入れる。MULT の Emittance[]・CALC の出口の端の蹴りは最後の半スライスの後 (追跡と同じ)。K1≠0・K0≠0 の BEND の追跡も AE を使う。
  追跡の一周写像の固有値 (閉軌道のまわり) = CALC のチューン (1e-7 台)。k64 は 3e-5〜1.4e-4 違う。
- B108: 端の角のある面では場の中の行路が x·tan(e) 短い。MULT の Emittance[]・追跡 (tmulti)・BEND ANGLE=0 の追跡 (tsteer) にも入れる (放射の行路)。
  硬い端の一様磁場の極の厳密な幾何 (円弧と直線) + ロビンソンの定理 (test_v2_E_steer._jx_exact) と、AE による ΔJx が一致する (比 1.03〜1.04)。K64 は約 2 倍。
"""
import io
import math

import numpy as np
import pytest

from pysad import compat
from . import _p3_gsdata as G
from . import _p3_indep as X
from . import test_v2_E_steer as T


def _calc_tracking_tunes(mode, src):
    old = compat._mode
    compat.set_mode(mode)
    try:
        return X.tracking_and_calc_tunes(src)
    finally:
        compat.set_mode(old)


_CALC_PRINT = 'FFS;\nUSE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nExit;\n'


@pytest.mark.parametrize('key', list(X.EMIT_CASES))
def test_mult_edge_angle_emittance_matches_gs_k64(k64_mode, key):
    """MULT の AE1・AE2 (K0 つき、F1・K2 つき) の CALC・Emittance[] が gs と全桁一致 (出口の端の蹴りを最後の半スライスの前に置く順序も)"""
    a, b = X.EMIT_CASES[key]
    d = T._values(X.run(X.mult_emit_script(a, b))[1])
    ref = G.EMIT[key]
    for k in ('calc', 'tunes'):
        assert np.allclose(d[k], ref[k], rtol=0, atol=3e-9), (k, d[k], ref[k])
    for k in ('damp', 'emit'):
        assert np.allclose(d[k], ref[k], rtol=1e-7, atol=0), (k, d[k], ref[k])
    for k in ('u0', 'alpha', 'eqpos'):
        assert abs(d[k] / ref[k] - 1.0) < 1e-7, (k, d[k], ref[k])


def test_mult_edge_angle_changes_the_result_v2_and_k64():
    """AE を読むようにした効果: AE のある MULT は AE なしと CALC のチューンが 1e-2 のオーダーで違う (以前の pysad は黙って無視していた)"""
    a, b = X.EMIT_CASES['e_a']
    for mode in ('k64', 'v2'):
        compat.set_mode(mode)
        try:
            d = T._values(X.run(X.mult_emit_script(a, b))[1])
            d0 = T._values(X.run(X.mult_emit_script('K0=0.01', 'K0=-0.01'))[1])
        finally:
            compat.set_mode('k64')
        assert abs(d['calc'][0] - d0['calc'][0]) > 1e-2


# ---------------------------------------------------------------- 追跡と CALC の一致 (別の経路)
_STEER_TUNE = ('K0=0.003 %s', 'K0=-0.003 %s')


@pytest.mark.parametrize('ae', ['AE2=0.3', 'AE1=0.3 AE2=-0.2'])
def test_steer_exit_edge_angle_tracking_equals_calc_v2(ae):
    """v2: FB2 = 0 でも追跡の出口の端の角 (AE2) が入り、追跡の一周写像の固有値 (閉軌道のまわり) が CALC (tsteee) のチューンと 6e-6 で一致
    (AE の 1 次の効果は 1.2e-4。残る差は K0 の 2 次の項: K64 の tsteer と tsteee の端の蹴りの 2 次の違いで K0 に比例して 2 乗で増える
    (K0 = 0.0015・0.003・0.006 で 4e-6・1.6e-5・6e-5)。K64 の仕様で直さない)"""
    src = X.ring_script('steer', 'K0=0.0015 ' + ae, 'K0=-0.0015 ' + ae, _CALC_PRINT)
    nu, calc = _calc_tracking_tunes('v2', src)
    assert np.allclose(nu, calc, atol=6e-6), (nu, calc)


def test_steer_exit_edge_angle_ignored_by_k64_tracking_without_fb2():
    """K64 (k64): 追跡の tsteer は AE2 を FB2 ≠ 0 のときだけ入れるので、AE2 だけのステアリングの追跡のチューンは AE なしと同じで CALC と 2e-4 違う
    (gs でも同じ: 追跡の一周写像から求めた x のチューンは AE なし 0.134095、AE2 だけ 0.134095、FB2=1e-6 を付けると 0.133883 で CALC (0.133471 + 配置のずれ) に近づく)。
    FB2 = 1e-6 (fb2 ≠ 0) にすると追跡と CALC が合う"""
    base = X.ring_script('steer', 'K0=0.003', 'K0=-0.003', _CALC_PRINT)
    ae = X.ring_script('steer', 'K0=0.003 AE2=0.3', 'K0=-0.003 AE2=0.3', _CALC_PRINT)
    aef = X.ring_script('steer', 'K0=0.003 AE2=0.3 FB2=1e-6 FRINGE=3', 'K0=-0.003 AE2=0.3 FB2=1e-6 FRINGE=3', _CALC_PRINT)
    base2 = X.ring_script('steer', 'K0=0.0015 AE2=0.3 FB2=1e-6 FRINGE=3', 'K0=-0.0015 AE2=0.3 FB2=1e-6 FRINGE=3', _CALC_PRINT)
    nu0, _c0 = _calc_tracking_tunes('k64', base)
    nu1, c1 = _calc_tracking_tunes('k64', ae)
    nu2, c2 = _calc_tracking_tunes('k64', aef)
    assert np.allclose(nu1, nu0, atol=1e-9)                   # AE2 が追跡に入っていない
    assert max(abs(np.array(nu1) - np.array(c1))) > 1e-4
    assert max(abs(np.array(nu2) - np.array(c2))) < 2e-5      # fb2 ≠ 0 なら入る (差は 1 次の効果 2e-4 の 1/10 以下。残りは K0 の 2 次)
    nu3, c3 = _calc_tracking_tunes('k64', base2)
    assert np.allclose(nu3, c3, atol=6e-6)


@pytest.mark.parametrize('case', ['AE1=0.3 AE2=0.3', 'AE2=0.3', 'AE1=0.3'])
def test_mult_edge_angle_tracking_equals_calc_v2_and_k64_exit_order_differs(case):
    """MULT (K0 = ±0.01、K1 = 0.7・−0.6、L = 0.3 m): v2 は追跡の一周写像の固有値と CALC のチューンが 1e-7 台で一致。
    k64 の CALC は出口の端の蹴りを最後の半スライスの前に置く (追跡の tmulti は後) ので AE2 のとき 3e-5〜1.2e-4 違う (入口 AE1 だけなら一致)"""
    a = case + ' K0=0.01'
    b = case + ' K0=-0.01'
    src = X.ring_script('mult', a, b, _CALC_PRINT)
    nu, calc = _calc_tracking_tunes('v2', src)
    assert np.allclose(nu, calc, atol=2e-7), (nu, calc)
    nu, calc = _calc_tracking_tunes('k64', src)
    dif = max(abs(np.array(nu) - np.array(calc)))
    if 'AE2' in case:
        assert dif > 2e-5, dif
    else:
        assert dif < 2e-7, dif


def test_bend_k1_k0_edge_angle_tracking_uses_ae_in_v2_only():
    """ANGLE=0・K0≠0・K1≠0 の BEND (K64 では tmulti / tmulte): CALC (tbende → tmulte) は AE を使う。K64 の追跡 (tbend → tmulti) は psi = 0 を渡すので
    k64 の追跡は AE を無視して CALC と 1e-4 以上違う。v2 の追跡は AE を使い CALC と 2e-7 で一致。(gs はこの BEND で SIGSEGV するので K64 の値は未確認)"""
    a = 'K0=0.003 K1=0.05 AE1=0.3 AE2=0.3'
    b = 'K0=-0.003 K1=-0.05 AE1=0.3 AE2=0.3'
    src = X.ring_script('steer', a, b, _CALC_PRINT)
    nu, calc = _calc_tracking_tunes('v2', src)
    assert np.allclose(nu, calc, atol=2e-6), (nu, calc)           # 残りは K0 の 2 次 (tsteer・tmulti と tmulte の端の蹴りの 2 次の違い)
    nu, calc = _calc_tracking_tunes('k64', src)
    assert max(abs(np.array(nu) - np.array(calc))) > 5e-4


# ---------------------------------------------------------------- 大きなリング A (3 GeV、603 m)
@pytest.mark.parametrize('key', ['ra1', 'ra2', 'ra3'])
def test_large_ring_a_tracking_tune_equals_calc_v2(key):
    """大きなリング A に補正磁石 2 台 (F1・FB・AE つきの BEND (ra1)・MULT (ra2)・AE だけ (ra3)) を入れる: v2 の追跡の一周写像のチューンが CALC と 3e-7 で一致。
    ra3 (FB2 = 0、AE2 ≠ 0) では k64 の追跡が CALC と 1e-4 違う (K64 の tsteer は fb2 = 0 の AE2 を入れない)"""
    src = X.ring_a_script(key, _CALC_PRINT)
    nu, calc = _calc_tracking_tunes('v2', src)
    assert np.allclose(nu, calc, atol=3e-7), (nu, calc)
    if key == 'ra3':
        nu, calc = _calc_tracking_tunes('k64', src)
        assert max(abs(np.array(nu) - np.array(calc))) > 5e-5


# ---------------------------------------------------------------- B108: 厳密な幾何 (円弧と直線) + ロビンソンの定理
def _pole_ring_mult(K0, ae, el, k1=0.0):
    s = T._pole_ring(K0, ae, 'EMITDIV=20;')
    s = s.replace('BEND  P1 = (L=0.2 ANGLE=0 K0=%s' % K0, '%s  P1 = (L=0.2 K1=%s K0=%s' % (el, k1, K0))
    return s.replace('BEND  P2 = (L=0.2 ANGLE=0 K0=%s' % (-2 * K0), '%s  P2 = (L=0.2 K1=%s K0=%s' % (el, k1, -2 * K0))


def _emit_jx_and_tracking_j(src):
    from pysad import ffscmd, track6n
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(src, out=buf)
    o = [l for l in buf.getvalue().split(chr(10)) if l.strip()]
    d = [float(x) for x in o[-2].strip('{}').split(',')]
    k = float(o[-1]) / 6e9
    els = ffs.lat.elements
    p = ffs.momentum()

    def f(v):
        return np.asarray(track6n.track_particles_fast(els, np.array([v]), 1, p, 0.0, True, True, False))[0][:6]
    h = np.array([1e-6] * 4 + [1e-5, 1e-6])

    def jac(v):
        J = np.zeros((6, 6))
        for i in range(6):
            dd = np.zeros(6)
            dd[i] = h[i]
            J[:, i] = (f(v + dd) - f(v - dd)) / (2 * h[i])
        return J
    v = np.zeros(6)
    for _ in range(8):
        v = v - np.linalg.solve(jac(v) - np.eye(6), f(v) - v)
    ev = np.linalg.eigvals(jac(v))
    jt = sorted(set(round(-math.log(abs(z)) / k, 9) for z in ev))
    return -d[0] / k, jt[0]


@pytest.mark.parametrize('el', ['MULT', 'BEND'])
@pytest.mark.parametrize('ae', [(0.3, 0.3), (-0.2, 0.4)])
def test_edge_angle_term_v2_matches_exact_geometry_emit_and_tracking(el, ae):
    """B108。極 3 個の組 (K0・−2K0・K0。24 セルの FODO、3 GeV、K0 = 0.005 rad/極) の AE による ΔJx (AE なしとの差) が、Emittance[] と追跡 (RF・放射つきの一周写像の
    固有値) の両方で厳密な幾何 + ロビンソンの定理と一致する (v2: 比 1.00〜1.06)。MULT (K0 つき、tmulte・tmulti) も ANGLE=0 の BEND (tsteee・tsteer) も。
    k64: Emittance[] は約 2 倍 (比 1.8〜2.4)。(K64 の追跡は BEND で出口の AE2 を入れないので比は不定: 比べない)"""
    K0 = 0.005
    ex = T._jx_exact(K0, *ae) - T._jx_exact(K0, 0.0, 0.0)
    assert abs(ex) > 3e-3

    def src(a):
        return _pole_ring_mult(K0, a, el) if el == 'MULT' else T._pole_ring(K0, a, 'EMITDIV=20;')
    out = {}
    for mode in ('v2', 'k64'):
        compat.set_mode(mode)
        try:
            e1, t1 = _emit_jx_and_tracking_j(src(ae))
            e0, t0 = _emit_jx_and_tracking_j(src((0.0, 0.0)))
        finally:
            compat.set_mode('k64')
        out[mode] = (e1 - e0, t1 - t0)
    assert 1.0 < out['v2'][0] / ex < 1.06, (ex, out)
    assert 1.0 < out['v2'][1] / ex < 1.06, (ex, out)
    assert 1.8 < out['k64'][0] / ex < 2.4, (ex, out)
    if el == 'MULT':
        assert 1.8 < out['k64'][1] / ex < 2.4, (ex, out)


def test_edge_angle_wedge_cancels_to_first_order_for_mult_v2(v2_mode):
    """AE1 = −AE2 (くさび): 行路の 1 次の変化が打ち消し、v2 の ΔJx は小さく (K0 = 0.005 で 1e-3 以下)、厳密な幾何と 1e-3 で一致。MULT の極も同じ"""
    K0 = 0.005
    ex = T._jx_exact(K0, 0.3, -0.3) - T._jx_exact(K0, 0.0, 0.0)
    e1, _t = _emit_jx_and_tracking_j(_pole_ring_mult(K0, (0.3, -0.3), 'MULT'))
    e0, _t = _emit_jx_and_tracking_j(_pole_ring_mult(K0, (0.0, 0.0), 'MULT'))
    assert abs(ex) < 2e-3 and abs((e1 - e0) - ex) < 1e-3, (ex, e1 - e0)


@pytest.mark.parametrize('ae', [(0.3, 0.3), (-0.2, 0.4)])
def test_edge_angle_term_v2_matches_exact_geometry_for_combined_function_poles(ae):
    """K1≠0 の MULT の極 (K0・−2K0・K0 の組に K1 = 0.06 (積分値) を足した勾配つき。B(x) = (K0 + K1 x)/L) の AE による ΔJx が、厳密な幾何
    (X.pole_k1: 運動の微分方程式の数値積分、端の面は実験室系で傾く) + ロビンソンの定理と、Emittance[] と追跡の両方で一致 (v2: 比 1.00〜1.06)。k64 は約 2 倍"""
    K0, K1 = 0.005, 0.06
    ex = X.jx_exact_k1(K0, K1, *ae) - X.jx_exact_k1(K0, K1, 0.0, 0.0)
    assert abs(ex) > 3e-3
    out = {}
    for mode in ('v2', 'k64'):
        compat.set_mode(mode)
        try:
            e1, t1 = _emit_jx_and_tracking_j(_pole_ring_mult(K0, ae, 'MULT', K1))
            e0, t0 = _emit_jx_and_tracking_j(_pole_ring_mult(K0, (0.0, 0.0), 'MULT', K1))
        finally:
            compat.set_mode('k64')
        out[mode] = (e1 - e0, t1 - t0)
    assert 1.0 < out['v2'][0] / ex < 1.06 and 1.0 < out['v2'][1] / ex < 1.06, (ex, out)
    assert 1.8 < out['k64'][0] / ex < 2.4, (ex, out)


def test_edge_angle_rotated_mult_emit_equals_tracking_v2(v2_mode):
    """K1 < 0 の MULT (K64 は多極全体を主軸へ回すので端の面の傾きの軸も回る (実験室系の y)): Emittance[] の行路の補正 (回した座標系の tiny・touty) と追跡の補正が
    同じ ΔJx を与える (差 3e-4 以下)。実験室系の x に傾く厳密な幾何 (K1 > 0) とは別の幾何なので、厳密な幾何との比較は K1 > 0 の試験"""
    K0, K1, ae = 0.005, -0.01, (0.3, 0.3)
    e1, t1 = _emit_jx_and_tracking_j(_pole_ring_mult(K0, ae, 'MULT', K1))
    e0, t0 = _emit_jx_and_tracking_j(_pole_ring_mult(K0, (0.0, 0.0), 'MULT', K1))
    assert abs((e1 - e0) - (t1 - t0)) < 3e-4, (e1 - e0, t1 - t0)
