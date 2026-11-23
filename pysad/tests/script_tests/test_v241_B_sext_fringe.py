# -*- coding: utf-8 -*-
"""v2.4.1 B: 6 極 (L≠0 の SEXT・OCT・DECA・DODECA) の端の非線形フリンジ ttfrin を CALC・Emittance[] の COD の写像とヤコビアンに入れる (v2、docs/30 B291)。

K64 の追跡 (tthin) は L≠0 の 6 極以上に非線形のフリンジ (Lee-Whiting 型、DISFRIN=0 が既定) を入れるが、CALC・Emittance[] の tthine は入れない。
フリンジの写像は座標の 3 次なので、軌道が 6 極の中心を通るとき (x = y = 0) はヤコビアンも COD も変わらず、理想のリングの結果は変わらない。
軌道が中心からずれているとき (ステアリング・ずれた 6 極・誤差) は、追跡の一周写像から取ったチューンが CALC・Emittance[] と合わなかった (下の FODO・6 極の
K2 = ±6、K0 = 0.012 のステアリングで 1.9e-3)。v2 は追跡と同じ式で入れる (elkern._fringe_step: 実数の式を複素ステップ微分して解析ヤコビアン)。
k64 は K64 のまま。独立な基準: 追跡 (track6) の一周写像 (閉軌道のまわりの数値微分) の固有値、1 要素の追跡のヤコビアンの DISFRIN の差。
限界: 本体 (4 次の α・β・γ の厚い 6 極) と CALC の tthine (2 回のキック) の積分の次数の違いで、軌道が数 mm〜1 cm ずれた L = 0.15 m の 6 極 (K2 L = 0.9) に
2〜3e-6 (ヤコビアン 1e-5) の差が残る (フリンジの追加とは別。短い 6 極ほど小さい)。
"""
import re

import numpy as np
import pytest

from pysad import compat
from pysad import elkern as ek
from pysad import lattice as lat
from pysad import track6 as t6
from . import _p3_indep as X


def _ring(dis, k0, k2=6.0, kn='K2', sext='SEXT'):
    return f'''MOMENTUM = 3GEV;
DRIFT D1 = (L=0.8) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
{sext} SF = (L=0.15 {kn}={k2} DISFRIN={dis}) SD = (L=0.15 {kn}={-k2} DISFRIN={dis});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND KK = (L=0.1 ANGLE=0 K0={k0});
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF SF D2 B D2 QD SD D1);
LINE RING = (ST KK 12*CELL1 RFC);
FFS;
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
'''


def _gap(mode, dis, k0):
    old = compat._mode
    compat.set_mode(mode)
    try:
        nu, calc = X.tracking_and_calc_tunes(_ring(dis, k0))
    finally:
        compat.set_mode(old)
    return float(np.max(np.abs(np.array(nu[-2:]) - np.array(calc))))


def test_ideal_orbit_is_unchanged_by_the_fringe():
    """軌道が 6 極の中心を通る (K0 = 0) ときは、フリンジの有無で CALC のチューンは変わらない (v2・k64 とも。追跡との差 7e-10)"""
    for mode in ('k64', 'v2'):
        assert _gap(mode, 0, 0.0) < 2e-9
        assert _gap(mode, 1, 0.0) < 2e-9


def test_k64_calc_misses_the_fringe_and_v2_includes_it():
    """K0 = 0.012 のステアリングで軌道が 6 極の中心からずれる。k64 の CALC は追跡と 1.9e-3 違い (K64 の tthine にフリンジが無い)、
    v2 は DISFRIN=1 (フリンジなし) のときと同じ小さな差 (4e-5。積分の次数の違い) まで合う。k64 の DISFRIN=1 の差も同じ 8e-5 台"""
    g_k64 = _gap('k64', 0, 0.012)
    g_v2 = _gap('v2', 0, 0.012)
    g_v2_off = _gap('v2', 1, 0.012)
    assert g_k64 > 1.5e-3
    assert g_v2 < 5e-5
    assert abs(g_v2 - g_v2_off) < 2e-7            # フリンジを入れた v2 の差 = フリンジのない系の差
    assert _gap('k64', 1, 0.012) < 1e-4


def test_emittance_tunes_include_the_fringe_v2_and_not_k64():
    """NORFSW NORAD の Emittance[] のチューン (6 次元の一周写像) = CALC のチューン (v2: 5e-8 以内、フリンジあり)。k64 の Emittance[] は DISFRIN に依らない"""
    tail = ('NORFSW; NORAD;\ne = Emittance[];\nPrint["Tunes ", Tunes/.e];\n')

    def emit(mode, dis):
        old = compat._mode
        compat.set_mode(mode)
        try:
            _f, out = X.run(_ring(dis, 0.012) + tail)
        finally:
            compat.set_mode(old)
        tu = [float(v) for v in re.search(r'\nTunes \{([^}]*)\}', out).group(1).split(',')]
        calc = [float(v) for v in re.search(r'CALC tunes\s+([-\d.eE+]+)\s+([-\d.eE+]+)', out).groups()]
        al = lambda x: min(x % 1.0, 1.0 - x % 1.0)
        return sorted([al(tu[0]), al(tu[1])]), sorted([al(calc[0]), al(calc[1])])
    e, c = emit('v2', 0)
    assert np.allclose(e, c, atol=5e-8), (e, c)
    e0, _c0 = emit('k64', 0)
    e1, _c1 = emit('k64', 1)
    assert np.allclose(e0, e1, atol=1e-12)


def test_fringe_map_equals_tracking_fringe_and_jacobian_is_exact():
    """elkern._fringe_step の写像が追跡 (track6._ttfrin_n) と 1e-18、解析ヤコビアンが数値微分と 1e-11 (n = nord/2 = 3 (6 極)・4・5。δ ≠ 0 を含む。
    CALC の正準運動量 px は追跡の px/(1+δ) に直して同じ式を使う)"""
    for n in (3, 4, 5):
        c0 = np.array([3e-3, 1e-3, -2e-3, 0.5e-3, 1e-4, 0.01])
        ak, al = 60.0 * n, 0.2
        pr = 1.0 + c0[5]
        p = t6._Particle.__new__(t6._Particle)
        p.c = np.array([c0[0], c0[1] / pr, c0[2], c0[3] / pr, c0[4], c0[5]])
        t6._ttfrin_n(p, 2 * n, ak, al)
        ref = np.array([p.c[0], p.c[1] * pr, p.c[2], p.c[3] * pr, p.c[4], p.c[5]])
        c = c0.copy()
        J = ek._fringe_step(c, ak, al, n, True)
        assert np.abs(c - ref).max() < 1e-17
        Jn = np.zeros((6, 6))
        for j in range(6):
            h = 1e-6
            a = c0.copy(); a[j] += h
            b = c0.copy(); b[j] -= h
            ek._fringe_step(a, ak, al, n, False)
            ek._fringe_step(b, ak, al, n, False)
            Jn[:, j] = (a - b) / (2 * h)
        assert np.abs(J - Jn).max() < 1e-11
        c2 = c.copy()
        J2 = ek._fringe_step(c2, -ak, al, n, True)         # 出口 (−ak) は入口とほぼ逆 (同じ座標での対ではないので厳密には逆でない)
        assert np.abs(J2 @ J - np.eye(6)).max() < 1e-4


def _one(sx, kn, K, dis):
    src = 'MOMENTUM = 3GEV;\n%s SA = (L=0.15 %s=%r DISFRIN=%d);\nMARK ST=();\nLINE ONE=(ST SA);\nFFS;\nUSE ONE;\n' % (sx, kn, K, dis)
    ffs, _ = X.run(src)
    return [e for e in ffs.lat.elements if e.name == 'SA'][0], ffs


def _track_jac(ffs, c0):
    def f(v):
        c = np.zeros((1, 6))
        c[0, :4] = v
        return np.asarray(t6.track_particles(ffs.lat.elements, c, 1, ffs.momentum(), 0.0, False, False, False))[0][:4]
    h = 1e-6
    J = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4)
        e[j] = h
        J[:, j] = (f(c0[:4] + e) - f(c0[:4] - e)) / (2 * h)
    return J, f(c0[:4])


@pytest.mark.parametrize('sx,kn,K', [('SEXT', 'K2', 12.0), ('OCT', 'K3', 2000.0)])
def test_fringe_part_of_the_element_jacobian_matches_tracking_v2(v2_mode, sx, kn, K):
    """独立な基準: 1 要素の追跡 (track6、数値微分) のヤコビアン・出口の座標の DISFRIN=0 と DISFRIN=1 の差 (= フリンジの寄与) が、CALC の要素の行列・COD の
    写像の差と一致する (軌道 x = 12 mm、y = −8 mm の 6 極 (K2 L = 1.8) で寄与は 2.4e-5。一致は 4.4e-4 の比 (本体の積分の次数の違いとの交差項)。COD は 1.3e-4 の比)"""
    c0 = np.array([1.2e-2, 1e-3, -8e-3, 0.5e-3, 0.0, 0.0])
    r = {}
    for dis in (0, 1):
        el, ffs = _one(sx, kn, K, dis)
        M, cout = lat._el_cod_step(el, c0, 0.0)
        Jt, ct = _track_jac(ffs, c0)
        r[dis] = (M[:4, :4], cout[:4], Jt, ct)
    dM = r[0][0] - r[1][0]
    dJ = r[0][2] - r[1][2]
    assert np.abs(dJ).max() > 1e-5
    assert np.abs(dM - dJ).max() < 1e-3 * np.abs(dJ).max()
    dc = r[0][1] - r[1][1]
    dt = r[0][3] - r[1][3]
    assert np.abs(dt).max() > 1e-8
    assert np.abs(dc - dt).max() < 5e-4 * np.abs(dt).max()


def test_k64_element_matrix_has_no_fringe():
    """k64 は K64 のまま: DISFRIN の 0・1 で CALC の要素の行列・COD が完全に同じ"""
    c0 = np.array([1.2e-2, 1e-3, -8e-3, 0.5e-3, 0.0, 0.0])
    old = compat._mode
    compat.set_mode('k64')
    try:
        el0, _f = _one('SEXT', 'K2', 12.0, 0)
        el1, _f = _one('SEXT', 'K2', 12.0, 1)
        m0, c_0 = lat._el_cod_step(el0, c0, 0.0)
        m1, c_1 = lat._el_cod_step(el1, c0, 0.0)
    finally:
        compat.set_mode(old)
    assert np.array_equal(m0, m1) and np.array_equal(c_0, c_1)


@pytest.mark.parametrize('dx,th', [(0.0, 0.0), (2e-3, 0.0), (2e-3, 0.3)])
def test_numba_and_python_paths_agree_v2(v2_mode, dx, th):
    """numba の k_thin (フリンジつき) と Python の _thin_cod_matrix が一致 (COD 1e-14、行列 1e-9。DX・ROTATE つきも)"""
    c0 = np.array([1.2e-2, 1e-3, -8e-3, 0.5e-3, 0.0, 0.0])
    L, ak = 0.15, 1.8
    r = lat._thin_step_fast(L, ak, 2, dx, 1e-3, th, c0, True)
    assert r is not None
    cf, Tf = r
    cp, Mp = lat._thin_cod_matrix(L, ak, 2, dx, 1e-3, th, c0, True)
    assert np.abs(cf - cp).max() < 1e-14
    assert np.abs(Tf - Mp).max() < 1e-9
    c1, M1 = lat._thin_cod_matrix(L, ak, 2, dx, 1e-3, th, c0, False)
    assert np.abs(cf - c1).max() > 1e-9              # フリンジの寄与がある


# ---------------------------------------------------------------- MULT の n ≥ 2 (tmulte の ttfrins)
def _mult_ring(dis, k0, extra=''):
    return f'''MOMENTUM = 3GEV;
DRIFT D1 = (L=0.8) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT SF = (L=0.15 K2=6 DISFRIN={dis} {extra}) SD = (L=0.15 K2=-6 DISFRIN={dis} {extra});
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND KK = (L=0.1 ANGLE=0 K0={k0});
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF SF D2 B D2 QD SD D1);
LINE RING = (ST KK 12*CELL1 RFC);
FFS;
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
'''


def _mult_gap(mode, dis, extra=''):
    old = compat._mode
    compat.set_mode(mode)
    try:
        nu, calc = X.tracking_and_calc_tunes(_mult_ring(dis, 0.012, extra))
    finally:
        compat.set_mode(old)
    return float(np.max(np.abs(np.array(nu[-2:]) - np.array(calc))))


@pytest.mark.parametrize('extra', ['', 'K3=20 ROTATE=0.1'])
def test_mult_sextupole_fringe_v2_matches_tracking_and_k64_does_not(extra):
    """MULT の K2 (K3 の ROTATE つきも) の非線形フリンジも v2 は CALC に入れる (追跡の tmulti と同じ ttfrins、複素の係数は主軸に回す)。ステアリング K0 = 0.012 の軌道のずれで、
    k64 の CALC は追跡と 1.8e-3〜6.8e-3 違い、v2 は DISFRIN=1 (フリンジなし) のときと同じ小さな差 (2〜7e-6。本体の積分の次数の違い)"""
    assert _mult_gap('k64', 0, extra) > 1e-3
    g0 = _mult_gap('v2', 0, extra)
    g1 = _mult_gap('v2', 1, extra)
    assert g0 < 1e-5 and abs(g0 - g1) < 1.5e-6, (g0, g1)


# ---------------------------------------------------------------- gs (K64) の確認 (k64)
_GS_RING = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=0.8) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.15 K2=3 DISFRIN=%d) SD = (L=0.15 K2=-3 DISFRIN=%d);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND KK = (L=0.1 ANGLE=0 K0=0.012);
CAVI RFC = (L=0 VOLT=3e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF SF D2 B D2 QD SD D1);
LINE RING = (ST KK 12*CELL1 RFC);
FFS;
$FORM="14.10";
bm={{1e-3,0,1e-3},{1e-4,0,0},{1e-3,0,-1e-3},{0,0,0},{0,0,0},{0,1e-3,-1e-2},{1,1,1}};
USE RING; CELL; CALC;
Print["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];
NORFSW; NORAD;
r=TrackParticles[{1,bm},1,1,3];
Print["N3 ", r[[2]]];
Exit;
'''
_GS_N3 = {
    0: [-0.0121682776, -0.0114377957, 0.0047217414, -0.0011968506, 0.000500405707, 0.0019876427, -0.0017689638, 0.0, 0.0018622433, 0.000608751822, -0.0, -0.00077331532, 0.0, 0.0, 0.0, 0.0, 0.001, -0.01],
    1: [-0.0123502911, -0.0117322488, 0.005040253, -0.0014810396, 0.000116190295, 0.0021520485, -0.0017589789, 0.0, 0.0018534117, 0.000594610772, -0.0, -0.00076633038, 0.0, 0.0, 0.0, 0.0, 0.001, -0.01],
}


@pytest.mark.parametrize('dis', [0, 1])
def test_gs_tracking_has_sextupole_fringe_and_calc_does_not_k64(k64_mode, dis):
    """gs (K64) の確認: 同じ 6 極のリング (K2 = ±3、K0 = 0.012 のステアリング) の追跡 (NORAD、3 周) は DISFRIN=0 と 1 で違い (フリンジが追跡に入る)、CALC のチューン
    (2.7025824135・1.8710392339) は DISFRIN に依らない (CALC・tthine にフリンジが無い)。k64 の pysad は追跡を 10 桁表示の全桁、CALC を 10 桁で再現する"""
    _f, out = X.run(_GS_RING % (dis, dis))
    assert np.allclose(X.tag_values(out, 'N3'), _GS_N3[dis], rtol=1e-8, atol=1e-15)
    m = re.search(r'CALC tunes\s+([-\d.eE+]+)\s+([-\d.eE+]+)', out)
    assert abs(float(m.group(1)) - 2.7025824135) < 2e-10 and abs(float(m.group(2)) - 1.8710392339) < 2e-10
    assert abs(_GS_N3[0][0] - _GS_N3[1][0]) > 1e-4
