# -*- coding: utf-8 -*-
"""v2 の修正 (担当 A): FFS まわり — B11 (4D の CALC の縦方向の Twiss 列)、B15 (UnstableLevel)、B16 (CalculateOptics の stabz・tracez)、
B21 (RF なしの閉軌道の δ)、B36 (L≠0 の TCAVI の Emittance[])、B63 (SOL の長さ)。

条件: 3 GeV、24 セルの FODO リング (周長 91.2 m)、RF 500 MHz 4 MV。
"""
import io
import math

import numpy as np

from pysad import ffscmd, emit6d, track6n

RING = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE RING = (ORIGIN RF 24*CELLQ);
FFS;
USE RING; CELL; RFSW;
'''


def _run(src):
    buf = io.StringIO()
    ffs, _I = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _vals(out, tag):
    line = [l for l in out.splitlines() if l.startswith(tag + ' ')][0]
    return [float(x) for x in line[len(tag) + 1:].replace('{', ' ').replace('}', ' ').replace(',', ' ').split()]


B11 = RING + '''CALC6D; CALC;
Print["A ",{Twiss["AZ",5],Twiss["BZ",5],Twiss["ZX",5],Twiss["ZPX",5]}];
o3=CalculateOptics[1,LINE["LENGTH"],Twiss["*",1],1,3];
Print["O3 ",o3[[2]]," ",o3[[3]]];
CALC4D; CALC;
Print["B ",{Twiss["AZ",5],Twiss["BZ",5],Twiss["ZX",5],Twiss["ZPX",5]}];
o2=CalculateOptics[1,LINE["LENGTH"],Twiss["*",1],1,2];
Print["O2 ",o2[[2]]," ",o2[[3]]];
'''


def test_calc4d_resets_longitudinal_columns_v2(v2_mode):
    """B11: 6D の CALC のあとの 4D の CALC は、縦方向の列を既定値 (AZ=0、BZ=1、ZX=ZPX=0) にする"""
    _f, out = _run(B11)
    a = _vals(out, 'A')
    assert abs(a[1] - 1.0) > 1.0 and a[2] != 0.0           # 6D では BZ・ZX が計算されている
    assert _vals(out, 'B') == [0.0, 1.0, 0.0, 0.0]


def test_calc4d_keeps_stale_columns_k64(k64_mode):
    """K64: 直前の 6D の先頭行の値が全行に残る"""
    _f, out = _run(B11)
    b = _vals(out, 'B')
    assert abs(b[1] - 1.0) > 1.0 and b[2] != 0.0


def test_calculate_optics_idim3_returns_z_trace_v2(v2_mode):
    """B16: idim=3 は縦方向の対角和 2cos(μz) と安定性を返す。idim=2 は stabz=1、tracez=0"""
    _f, out = _run(B11 + 'e=Emittance[];\nPrint["T ",Tunes/.e];\n')
    o3 = _vals(out, 'O3')
    nuz = _vals(out, 'T')[2]
    assert o3[2] == 1.0
    assert abs(o3[5] - 2.0 * math.cos(2 * math.pi * nuz)) < 2e-3      # CALC6D (放射なしの正規形) と Emittance[] のチューンの差の範囲
    assert _vals(out, 'O2')[2] == 1.0 and _vals(out, 'O2')[5] == 0.0


def test_calculate_optics_idim3_k64(k64_mode):
    _f, out = _run(B11)
    o3 = _vals(out, 'O3')
    assert o3[2] == 1.0 and o3[5] == 0.0                    # pysad の k64: 未初期化の代わりに 1 と 0


B15 = RING.replace('USE RING; CELL; RFSW;', 'USE RING; CELL;') + '''QF K1 2.5;
FIT $$$ NY 4.07;
FREE QD;
GO;
Print["U ",{UnstableLevel,MatchingResidual,NetResidual}];
'''


def test_unstable_level_v2(v2_mode):
    """B15: 不安定のまま終わったマッチングで UnstableLevel = nstab (MatchingResidual = 10·nstab + r)"""
    _f, out = _run(B15)
    u = _vals(out, 'U')
    assert u[0] == 2.0
    assert abs(u[1] - (10.0 * u[0] + u[2])) < 1e-9


def test_unstable_level_k64_is_zero(k64_mode):
    _f, out = _run(B15)
    u = _vals(out, 'U')
    assert u[0] == 0.0 and u[1] == 20.0


def test_no_rf_closed_orbit_closes_z():
    """B21: RFSW OFF・RADCOD OFF の閉軌道は「z が閉じる δ」(行路 = 設計周長)。v2 でも変えていない。
    K0 の誤差のある BEND で δ = -1.4e-4 になり、別のコード (追跡) で 1 周して z のずれが 1e-9 m 以下"""
    src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QE = (L=0.3 K1=0.7 DX=1E-3);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5) BE = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5 K0=2E-4);
MARK  ORIGIN = ();
LINE CELL = (QF D2 B D2 QD D1);
LINE CELLE = (QE D2 BE D2 QD D1);
LINE RING = (ORIGIN 20*CELL CELLE 3*CELL);
FFS;
USE RING; CELL; CALC;
'''
    from pysad import compat
    for mode in ('k64', 'v2'):
        old = compat._mode
        compat.set_mode(mode)
        try:
            ffs, _out = _run(src)
            e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False)
            cod = np.array(e['orbit_entrance'])
            assert abs(cod[5] + 1.428e-4) < 2e-6
            r = np.asarray(track6n.track_particles_fast(ffs.lat.elements, np.array([cod]), 1, ffs.momentum(), 0.0,
                                                        True, False, False))[0][:6]
            assert abs(r[4] - cod[4]) < 1e-9
            assert np.abs(r[:4] - cod[:4]).max() < 5e-9
        finally:
            compat.set_mode(old)


def _emit_tc(tc):
    src = RING.replace('MARK  ORIGIN = ();', tc + '\nMARK  ORIGIN = ();').replace('(ORIGIN RF 24*CELLQ)', '(ORIGIN RF TC 24*CELLQ)') + \
        'CALC; RAD;\ne=Emittance[];\nPrint["E ",{EnergyLossU0,Emittances,Tunes,MomentumCompaction,BunchLength}/.e];\n'
    _f, out = _run(src)
    return _vals(out, 'E')


def test_tcavi_with_length_emittance_is_finite_and_has_drift_limit(v2_mode):
    """B36: L≠0 の TCAVI を含むリングの Emittance[] は有限 (gs は NaN)。K0=0 では同じ長さのドリフトと一致し、K0 を入れても連続に変わる"""
    d = _emit_tc('DRIFT TC = (L=0.4);')
    t0 = _emit_tc('TCAVI TC = (L=0.4 K0=0 FREQ=4.5e8 PHI=0.1);')
    t1 = _emit_tc('TCAVI TC = (L=0.4 K0=1e-5 FREQ=4.5e8 PHI=0.1);')
    assert all(math.isfinite(x) for x in t1)
    assert np.allclose(d, t0, rtol=1e-12, atol=1e-40)
    assert abs(t1[0] / d[0] - 1) < 1e-6 and abs(t1[1] / d[1] - 1) < 1e-3 and abs(t1[4] - d[4]) < 1e-5


def test_sol_has_no_length():
    """B63: SOL に長さのキーワードは無く (SADHelp)、L を書いても軌道・光学・S 座標・周長に入らない。
    gs: 下のリングの LENG は 46.6 (SOL の L=0.3 を入れない)、MomentumCompaction .029609026292907528"""
    src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
SOL SO = (L=0.3 BZ=0.0 BOUND=1 GEO=1);
SOL SO2 = (L=0 BZ=0 BOUND=1);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1 SO D1 SO2 6*CELL1 RFC);
FFS;
USE RING; CELL; CALC;
Print["LENG ",{LINE["LENG","$$$"],LINE["S","$$$"]}];
RFSW; RAD;
e=Emittance[];
Print["E ",{MomentumCompaction,EnergyLossU0,BunchLength}/.e];
'''
    from pysad import compat
    for mode in ('k64', 'v2'):
        old = compat._mode
        compat.set_mode(mode)
        try:
            _f, out = _run(src)
            assert np.allclose(_vals(out, 'LENG'), [46.6, 46.6], rtol=1e-13)
            if mode == 'k64':
                assert abs(_vals(out, 'E')[0] / .029609026292907528 - 1) < 1e-9
            else:
                assert abs(_vals(out, 'E')[0] / .029609026292907528 - 1) < 1e-3      # v2 は α を 4 sin²(μ/2) で求める (ALPHA)
        finally:
            compat.set_mode(old)
