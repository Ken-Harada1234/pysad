# -*- coding: utf-8 -*-
"""v2.4.1 A: ANGLE=0 で K1≠0 かつ K0≠0 の BEND (F1・FB1・FB2 を含む。tbend の tmulti、tbende の tmulte の経路)。

(1) gs (K64 の 10-19 ビルド) の BEND はこの組み合わせを正しく扱わない (docs/30 B290): 追跡は tmulti に rtaper = 0 (掛け算の係数) を渡すので K0・K1 とも
    無視して**ただのドリフト**になり (nmmax の判定も逆)、CALC・Emittance[] は tmulte の引数が 2 個足りず SIGSEGV で落ちる。pysad はこれを再現せず
    (どちらのモードも)、tbend が tmulti に渡すはずの値 (akm(0) = K0、akm(1) = K1·L) で追跡・CALC・Emittance[] を計算する = **同じ K0・K1・L の MULT と同じ**。
    したがって gs との照合は gs の MULT (tmulti・tmulte の正しい呼び出し) で行う (FB1 = F1 + FB1。BEND の F1 は MULT の FB に行く)。
    gs の値は 12 セルの FODO (3 GeV、周長 45.6 m。ステアリング磁石 2 種を 6 セルずつ) の追跡 (NORAD・RAD、3 周)・CALC・Emittance[] (10 桁表示の全桁)。
(2) 高速追跡 (numba の K_MULT) と Python (track6) が 1e-13 で一致 (以前は常に Python に落ちていた)。F1・FB・AE つきは Python のまま。
(3) 独立な基準: 一様な K0 と K1 (B = (K0 + K1 x)/L) の中の運動方程式 (正準方程式、DOP853 rtol 1e-13) と本体が 1e-8 以内、v2 の追跡の一周写像の固有値 =
    CALC のチューン (AE・F1・FB・DX・ROTATE つき。k64 は追跡が AE を無視する K64 の仕様 (B142) と ROTATE の CALC の癖で違う)。
"""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat
from . import _p3_indep as X
from . import test_v2_E_steer as T

# gs の MULT の値 (台本は _MK で作る。MULT K0=±、K1=0.05・−0.04 (12 セル FODO))。tag_values の最初の 18 個
GS = {
    'plain': {
        'N3': [-0.0158609541, -0.0120227366, -0.0213235933, -0.0044457145, -0.0031774306, -0.0051673678, -0.00085095678, 0.0, -0.00027892509, 0.0011755767, -0.0, -0.00055490738, 0.0, 0.0, 0.0, 0.0, 0.001, -0.01],
        'R3': [-0.0160448632, -0.0122813292, -0.0212606486, -0.004449556, -0.0032180543, -0.0050845118, -0.00080876981, 0.0, -0.00031314183, 0.0011578591, -0.0, -0.00053101061, -0.0046284845, -0.0080894089, 0.0313800963, -0.0008674211, 7.75332139e-05, -0.0102886089],
        'emit': {'calc': [3.1126377227, 2.2778149058], 'tunes': [0.1193897296, 0.2829585002, -0.0155385108], 'u0': 779806.611971, 'alpha': 0.02454764, 'eqpos': 0.0329573943, 'emit': [2.50529979e-07, 1.78521317e-05], 'damp': [-0.00012523989, -0.00013017727, -0.00026437997]},
    },
    'fringe': {
        'N3': [-0.0158609546, -0.0120227383, -0.0213235746, -0.0044457147, -0.0031774313, -0.0051673614, -0.00085142055, 0.0, -0.00027838109, 0.0011757594, -0.0, -0.00055527447, 0.0, 0.0, 0.0, 0.0, 0.001, -0.01],
        'R3': [-0.0160448627, -0.0122813299, -0.0212606293, -0.0044495559, -0.0032180545, -0.0050845053, -0.00080924008, 0.0, -0.00031259873, 0.0011580509, -0.0, -0.00053138146, -0.004628361, -0.0080892859, 0.0313802195, -0.0008674191, 7.75351891e-05, -0.010288607],
        'emit': {'calc': [3.1126377227, 2.2777994438], 'tunes': [0.1193895881, 0.2829428703, -0.0155385114], 'u0': 779806.665581, 'alpha': 0.0245476421, 'eqpos': 0.0329573967, 'emit': [2.50530729e-07, 1.78521235e-05], 'damp': [-0.00012523975, -0.00013017728, -0.00026438014]},
    },
}
# (BEND の引数 a・b、同じ磁石の MULT の引数 a・b)
BEND_MULT = {
    'plain': ('K1=0.05 K0=0.002', 'K1=-0.04 K0=-0.003', 'K1=0.05 K0=0.002', 'K1=-0.04 K0=-0.003'),
    'fringe': ('K1=0.05 K0=0.002 F1=0.03 FB1=0.02 FB2=0.01 FRINGE=3', 'K1=-0.04 K0=-0.003 F1=0.03 FB1=0.02 FB2=0.01 FRINGE=3',
               'K1=0.05 K0=0.002 FB1=0.05 FB2=0.04 FRINGE=3', 'K1=-0.04 K0=-0.003 FB1=0.05 FB2=0.04 FRINGE=3'),
}
_TAIL = X.TRACK_TAIL.replace('USE RING; CELL; CALC; NORFSW; NORAD;', 'USE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nNORFSW; NORAD;')
_TAIL = _TAIL.replace('Exit;', '''e = Emittance[];
Print["Tunes ", Tunes/.e];
Print["U0 ", EnergyLossU0/.e, " alpha ", MomentumCompaction/.e, " EqPos ", EquilibriumPosition/.e];
Print["Emittances ", Emittances/.e];
Print["Damp ", DampingRate/.e];
Exit;''')


def _mult_ring(a, b, tail):
    """ステアリング位置 (SA・SB) を MULT (ANGLE なし) にした同じリング"""
    return (X._STEER.replace('BEND SA = (L=0.2 ANGLE=0 %s)', 'MULT SA = (L=0.2 %s)').replace('BEND SB = (L=0.2 ANGLE=0 %s)', 'MULT SB = (L=0.2 %s)')) % (a, b) + tail


def _vals(out, tag):
    return np.array(X.tag_values(out, tag))


@pytest.mark.parametrize('key', list(GS))
def test_bend_k1_k0_matches_gs_mult_k64(k64_mode, key):
    """k64: BEND (ANGLE=0、K1≠0、K0≠0) の追跡 (NORAD・RAD)・CALC・Emittance[] が、gs の MULT (同じ K0・K1・L) と 10 桁表示の全桁で一致。gs の BEND は B290 で使えない"""
    a, b, *_ = BEND_MULT[key]
    out = X.run(X.ring_script('steer', a, b, _TAIL))[1]
    assert np.allclose(_vals(out, 'N3'), GS[key]['N3'], rtol=1e-8, atol=1e-17)
    assert np.allclose(_vals(out, 'R3'), GS[key]['R3'], rtol=1e-8, atol=1e-17)
    T._check(T._values(out), GS[key]['emit'])


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('key', list(GS))
def test_bend_k1_k0_equals_mult_in_both_modes(mode, key):
    """BEND (ANGLE=0、K1≠0、K0≠0) は同じ K0・K1・L の MULT と、追跡・CALC・Emittance[] が全桁同じ (AE・F1 なし。F1 は MULT の FB = F1 + FB)"""
    a, b, ma, mb = BEND_MULT[key]
    old = compat._mode
    compat.set_mode(mode)
    try:
        ob = X.run(X.ring_script('steer', a, b, _TAIL))[1]
        om = X.run(_mult_ring(ma, mb, _TAIL))[1]
    finally:
        compat.set_mode(old)
    assert np.allclose(_vals(ob, 'N3'), _vals(om, 'N3'), rtol=0, atol=1e-15)
    assert np.allclose(_vals(ob, 'R3'), _vals(om, 'R3'), rtol=0, atol=1e-15)
    vb, vm = T._values(ob), T._values(om)
    for k in vb:
        assert np.allclose(vb[k], vm[k], rtol=1e-12, atol=1e-14), k


# ---------------------------------------------------------------- 高速追跡と Python
def _fast_vs_python(mode, a, b, rad, fluc):
    from pysad import track6 as t6, track6n as tn, emit6d as e6
    compat.set_mode(mode)
    ffs, _ = X.run(X.ring_script('steer', a, b, 'FFS;\nUSE RING; CELL; CALC;\n'))
    els = ffs.lat.elements
    p0 = ffs.momentum()
    ph = e6._Phys(p0, sum(el.L for el in els))
    prog = tn.Prog(els, ph)                                   # F1・FB・AE が無ければ NotImplementedError にならない
    assert sum(1 for k in prog.kinds if k == 10) == 12            # K_MULT (numba)
    rng = np.random.default_rng(5)
    c = np.zeros((20, 6))
    for i, s in enumerate((3e-4, 1e-4, 2e-4, 1e-4, 1e-3, 5e-4)):
        c[:, i] = rng.normal(0, s, 20)
    ra = t6.track_particles(els, c, 5, p0, trf0=1e-4, rad=rad, fluc=fluc, seed=777)
    sa = t6.RNG.seed
    rb = tn.track_particles_fast(els, c, 5, p0, trf0=1e-4, rad=rad, fluc=fluc, seed=777)
    sb = t6.RNG.seed
    return np.abs(np.asarray(ra) - np.asarray(rb)).max(), sa == sb


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('rad,fluc', [(False, False), (True, False), (True, True)])
@pytest.mark.parametrize('a,b', [('K1=0.05 K0=0.003', 'K1=-0.04 K0=-0.003'),
                                 ('K1=0.05 K0=0.003 ROTATE=0.1 DX=1e-4', 'K1=-0.04 K0=-0.003 DROTATE=0.2 DY=1e-4 DISFRIN=1')])
def test_bend_k1_k0_fast_equals_python(mode, rad, fluc, a, b):
    """numba の K_MULT と Python の _tmulti_ens が 1e-13 で一致 (FLUC は乱数の最後の状態も)"""
    old = compat._mode
    try:
        d, same = _fast_vs_python(mode, a, b, rad, fluc)
    finally:
        compat.set_mode(old)
    assert d < 1e-13 and same, (d, same)


def test_bend_k1_k0_with_f1_or_ae_falls_back_to_python():
    """F1・FB・AE つきは numba の K_MULT が扱わず (NotImplementedError)、track_particles_fast が Python に落ちて同じ値を返す"""
    from pysad import track6 as t6, track6n as tn, emit6d as e6
    old = compat._mode
    compat.set_mode('v2')
    try:
        ffs, _ = X.run(X.ring_script('steer', 'K1=0.05 K0=0.003 F1=0.03 FRINGE=3 AE1=0.1', 'K1=-0.04 K0=-0.003 AE2=0.1', 'FFS;\nUSE RING; CELL; CALC;\n'))
        els = ffs.lat.elements
        p0 = ffs.momentum()
        with pytest.raises(NotImplementedError):
            tn.Prog(els, e6._Phys(p0, sum(el.L for el in els)))
        c = np.zeros((3, 6))
        c[:, 0] = [1e-4, -2e-4, 3e-4]
        ra = t6.track_particles(els, c, 3, p0, trf0=1e-4, rad=False, fluc=False, seed=1)
        rb = tn.track_particles_fast(els, c, 3, p0, trf0=1e-4, rad=False, fluc=False, seed=1)
    finally:
        compat.set_mode(old)
    assert np.abs(np.asarray(ra) - np.asarray(rb)).max() < 1e-15


# ---------------------------------------------------------------- 独立な基準
def _exact(c, L, K0, K1, d):
    """一様な B = (K0 + K1 x)/L の中の運動 (px は p0 で規格化。H = −√((1+δ)² − px² − py²) + (K0 x + K1 (x² − y²)/2)/L)"""
    def rhs(s, v):
        x, px, y, py = v
        pz = math.sqrt((1 + d)**2 - px * px - py * py)
        return [px / pz, -(K0 + K1 * x) / L, py / pz, K1 * y / L]
    return solve_ivp(rhs, (0, L), c[:4], method='DOP853', rtol=1e-13, atol=1e-18).y[:, -1]


@pytest.mark.parametrize('mode,tol', [('k64', 6e-9), ('v2', 1e-9)])
@pytest.mark.parametrize('L,K0,K1', [(0.3, 0.006, -0.6), (0.3, 0.006, 0.9), (0.2, 0.003, 0.05)])
@pytest.mark.parametrize('delta', [0.0, 0.01])
def test_bend_k1_k0_body_matches_equation_of_motion(mode, tol, L, K0, K1, delta):
    """本体 (DISFRIN=1: フリンジなし) の出口の x・px・y・py が運動方程式の数値積分と一致 (k64 の Strang 分割は 4e-9 まで、v2 の 4 次の合成は 6e-10 まで)。
    粒子は x = 3 mm、px = 2 mrad、y = 2 mm、py = 1 mrad"""
    c = np.array([3e-3, 2e-3, 2e-3, 1e-3, 0.0, delta])
    ex = _exact(c, L, K0, K1, delta)
    old = compat._mode
    compat.set_mode(mode)
    try:
        r = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r K1=%r DISFRIN=1);' % (L, K0, K1), c, fast=True)
    finally:
        compat.set_mode(old)
    assert np.abs(r[:4] - ex).max() < tol, (r[:4] - ex)


_TUNE_TAIL = 'FFS;\nUSE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\n'


@pytest.mark.parametrize('a,b', [
    ('K1=0.05 K0=0.003', 'K1=-0.04 K0=-0.003'),
    ('K1=0.05 K0=0.006 AE1=0.2 AE2=-0.1', 'K1=-0.04 K0=-0.006 AE1=0.2 AE2=-0.1'),
    ('K1=0.05 K0=0.006 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1', 'K1=-0.04 K0=-0.006 F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1'),
    ('K1=0.05 K0=0.006 DX=2e-4 ROTATE=0.1', 'K1=-0.04 K0=-0.006 DROTATE=0.2 DY=1e-4')])
def test_bend_k1_k0_tracking_tunes_equal_calc_v2(v2_mode, a, b):
    """v2: 追跡 (RF・放射なし) の閉軌道のまわりの 4 次元の一周写像の固有値のチューン = CALC のチューン (1e-7 以内。AE・F1・FB・DX・ROTATE つき)。
    追跡と CALC が別の経路 (tmulti 相当と tmulte 相当) なので、閉軌道・光学の独立な確認になる"""
    nu, calc = X.tracking_and_calc_tunes(X.ring_script('steer', a, b, _TUNE_TAIL))
    assert np.allclose(nu[-2:], calc, atol=2e-7), (nu, calc)


def test_bend_k1_k0_ae_is_used_by_v2_tracking_only(v2_mode):
    """v2 は追跡も AE を使う (B142)。AE を付けるとチューンが変わる (K64 は追跡が AE を無視するので CALC と 7e-3 違った)"""
    nu0, _c = X.tracking_and_calc_tunes(X.ring_script('steer', 'K1=0.05 K0=0.006', 'K1=-0.04 K0=-0.006', _TUNE_TAIL))
    nu1, _c = X.tracking_and_calc_tunes(X.ring_script('steer', 'K1=0.05 K0=0.006 AE1=0.2 AE2=-0.1', 'K1=-0.04 K0=-0.006 AE1=0.2 AE2=-0.1', _TUNE_TAIL))
    assert abs(nu1[-1] - nu0[-1]) > 1e-3
