# -*- coding: utf-8 -*-
"""v2 の修正 B65: RADLIGHT の軌道 (TrackParticles の {t, x, y, z})。K64 は横方向の座標が入らず、BEND で向きが変わらず、BEND の中を記録しない。

根拠:
  - 振幅 0 の粒子の軌道は設計軌道: 記録した点が幾何 LINE["GX"/"GY"/"GZ"] (ROTATE つき・負の角の BEND を含む) に一致。時刻は行路長/(βc)。
  - 横にずれた粒子は、幾何の横の軸の向きにずれる (始点の距離が √(x²+y²))。
  - 偏向磁石の放射場のピーク: 接線上の遠方の観測点で E = e β²/(4π ε0 ρ R (1-β)²) (Liénard-Wiechert の加速度場、n ∥ β)。
    RadiationField (Feynman-Heaviside の式、gs では SIGSEGV) に v2 の軌道を入れた値がこれに一致する。
"""
import io
import math

SRC1 = r'''MOMENTUM = 10.2069e6;
DRIFT D1 = (L=0.5);
BEND B = (L=0.5 ANGLE=0.5 ROTATE=0.3) B2=(L=0.3 ANGLE=-0.2);
QUAD Q=(L=0.2 K1=0.1);
MARK ST = ();
LINE BL = (ST D1 B D1 Q B2 D1);
FFS;
USE BL;
INS; CALC;
RADLIGHT; NORAD;
$FORM="20.14";
r = TrackParticles[{1, {{0,1e-3},{0,0},{0,2e-3},{0,0},{0,0},{0,0},{1,1}}}, 8, 1, 1];
t=r[[2,1]];
Print["G ",LINE["GX",8]," ",LINE["GY",8]," ",LINE["GZ",8]," ",LINE["GX",4]," ",LINE["GY",4]," ",LINE["GZ",4]];
Print["E ",t[[2,-1]]," ",t[[3,-1]]," ",t[[4,-1]]," ",t[[1,-1]]," ",Length[t[[1]]]];
i=41;
Print["M ",t[[2,i+1]]," ",t[[3,i+1]]," ",t[[4,i+1]]];
t2=r[[2,2]];
Print["P ",t2[[2,1]]," ",t2[[3,1]]," ",t2[[4,1]]];
Exit;
'''


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return {l.split()[0]: [float(v) for v in l.split()[1:]] for l in buf.getvalue().splitlines()
            if l.split() and l.split()[0] in ('G', 'E', 'M', 'P', 'EM')}


def test_b65_trajectory_follows_geometry(v2_mode):
    o = _run(SRC1.replace('RADLIGHT;', 'RADLIGHT; RADLIGHTDIV=4;'))
    g = o['G']
    e = o['E']
    for a, b in zip(e[:3], g[:3]):
        assert abs(a - b) < 1e-11
    beta = math.sqrt(1.0 - 1.0 / (1.0 + (10.2069e6 / 510998.95)**2))
    assert abs(e[3] / (2.5 / (beta * 299792458.0)) - 1.0) < 1e-9
    # 点の数: 始点 1 + D1 1 + B 40 (0.5 rad × γ 20 × 4) + D1 1 + Q 1 + B2 16 + D1 1
    assert int(e[4]) == 61
    # B の出口 (始点 + D1 + 40 点 = 42 番目) は要素 4 の入口の幾何
    for a, b in zip(o['M'], g[3:]):
        assert abs(a - b) < 1e-11
    # 横にずれた粒子: 始点は軸に垂直な面内で √(x²+y²) だけ離れる
    p = o['P']
    assert abs(p[0]) < 1e-15 and abs(math.hypot(p[1], p[2]) - math.hypot(1e-3, 2e-3)) < 1e-15


SRC2 = r'''MOMENTUM = 10.2069e6;
DRIFT D1 = (L=0.5);
BEND B = (L=0.5 ANGLE=0.5);
MARK ST = ();
LINE BL = (ST D1 B D1);
FFS;
USE BL;
INS; CALC;
RADLIGHT; NORAD;
$FORM="20.14";
r = TrackParticles[{1, {{0},{0},{0},{0},{0},{0},{1}}}, 5, 1, 1];
t=r[[2,1]];
R=100;
sg=Sign[LINE["GY",4]];
obs={0.5+Sin[0.25]+R Cos[0.25], sg (1-Cos[0.25])+sg R Sin[0.25], 0};
f=RadiationField[t,obs];
Print["EM ",Max[Sqrt[f[[2]]^2+f[[3]]^2]]," ",Length[f[[2]]]];
Exit;
'''


def test_b65_bend_radiation_peak_field(v2_mode):
    o = _run(SRC2)
    g = math.hypot(1.0, 10.2069e6 / 510998.95)
    b = math.sqrt(1.0 - 1.0 / g**2)
    exp = 1.602176634e-19 / (4 * math.pi * 8.8541878128e-12) * b * b / (1.0 * 100.0 * (1.0 - b)**2)
    assert abs(o['EM'][0] / exp - 1.0) < 5e-3              # RADLIGHTDIV の既定 16 (1/γ あたり 16 点)
    assert int(o['EM'][1]) == 163


def test_b65_k64_mode_has_no_bend_points(k64_mode):
    o = _run(SRC2)
    assert int(o['EM'][1]) == 3                              # K64: 始点と 2 つのドリフトの後だけ
