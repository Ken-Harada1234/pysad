"""TCAVI (横偏向空洞、K64 の ttcav.f / ttcave.f)。参照は SAD K64 (gs) の出力。
tct*: TrackParticles (RFSW ON、全周を 3 周)。tct7、tct8 は RAD FLUC (20 粒子: 乱数列まで一致)。
tce4、tce5: Emittance[] (L=0 の TCAVI。tce5 はずれと回転付きで、COD・チューン・運動量収縮因子だけ比べる)。
L≠0 の Emittance[] は gs が NaN や不自然な値を出す (docs/30 B36) ので比べない。
参照の作り方: local/wake/gentc.py (docs/19_other_elements.md)。"""
import io
import re

import pytest

TCT1_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0.4 K0=1e-4 FREQ=4.5e8 PHI=0.1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
TCT1_REF = [[1.0, 0.0021211997, 0.0010910073, 0.00037940659036, 0.0008683068486], [2.0, 0.0021272585, 0.00034874837749, 0.00012040389165, 0.00023714447459], [3.0, -0.0001211541409, -4.672806647e-05, 4.5216036841e-06, 2.5480275072e-05], [4.0, 0.00030851656895, 9.9515773964e-05, 1.8143983459e-05, -4.947268386e-05], [5.0, 0.0207742688, -0.0021465827, -0.0023375958, 0.0021116427], [6.0, -0.0023425695, 0.00093632915055, 0.00027632316958, 0.00070955311993]]

TCT2_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0.4 K0=1e-4 FREQ=4.5e8 PHI=0.1 DX=2e-4 DY=-1e-4 ROTATE=0.05);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
TCT2_REF = [[1.0, 0.0019942535, 0.0010903574, 0.00037876275152, 0.00086767375397], [2.0, 0.0020774396, 0.00034849017479, 0.00012014933849, 0.00023689355124], [3.0, -0.0001663229958, -4.711463741e-05, 4.2088666856e-06, 2.5220698136e-05], [4.0, 0.00036634550581, 9.9904517254e-05, 1.8511494983e-05, -4.921916173e-05], [5.0, 0.020951836, -0.0021457077, -0.0023367198, 0.0021125499], [6.0, -0.0025728836, 0.00093517481572, 0.00027517239444, 0.00070839714339]]

TCT3_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0.4 K0=1e-4 HARM=64.3 PHI=0.2);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
TCT3_REF = [[1.0, -0.0017545943, 0.0010626453, 0.00035101946947, 0.00085579687425], [2.0, 0.0055202659, 0.00036202496459, 0.00013418568466, 0.00025563759626], [3.0, -0.0001212482865, -4.674591096e-05, 4.5152554094e-06, 2.5490863692e-05], [4.0, 0.00030857756471, 9.9528255739e-05, 1.81476292e-05, -4.948002022e-05], [5.0, 0.0242366752, -0.0021339169, -0.0023271251, 0.0021368213], [6.0, -0.0022254974, 0.00093402146646, 0.00027603168565, 0.00070811380456]]

TCT5_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0 K0=1e-4 FREQ=4.5e8 PHI=0.1 ROTATE=0.5);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
TCT5_REF = [[1.0, -0.000536892153, 0.00087931614086, 0.00032144794148, 0.00087327958094], [2.0, 0.00036538993704, 0.00041399666402, 0.00015696447531, 0.0003996099298], [3.0, -0.0041122212, -0.0002056266731, -7.232621015e-05, 7.4837322854e-05], [4.0, 0.0028996884, 0.00011037471192, 4.8796832086e-05, -3.513873317e-05], [5.0, 0.0200399163, -0.0021421583, -0.0023297971, 0.0021496597], [6.0, -0.0023454024, 0.00093585228384, 0.00027557491605, 0.0007103516939]]

TCT7_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0.4 K0=1e-4 FREQ=4.5e8 PHI=0.1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; RAD; FLUC;
n = 20;
x0  = Table[1e-5*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-5*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TCT7_REF = [[1.0, -0.0031427713, 0.00096176258685, -0.0008185374308, 0.0010434499], [2.0, -0.0011439784, 0.00028428557856, -0.0002602524321, 0.0003254691063], [3.0, -2.622810787e-05, -1.556107493e-05, 1.6969386629e-05, -1.831642224e-05], [4.0, -2.380285466e-06, 2.1173510757e-05, -1.84855832e-05, 1.8187260375e-05], [5.0, 0.0307401028, -0.0019594733, 0.0026230573, -0.0006480910386], [6.0, -0.0054759229, 0.0006982963805, -0.000845698706, 0.00071942812576]]

TCT8_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0.4 K0=3e-4 FREQ=4.5e8 PHI=0.1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
TRPT; INS; RAD; FLUC;
n = 20;
x0  = Table[1e-5*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-5*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.003*Sin[0.13*i]+1e-5*i-0.001, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
RFSW;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 52, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TCT8_REF = [[1.0, -0.0023836886, 0.0010139718, -0.0007848370658, 0.0010899043], [2.0, -0.0007751023117, 0.00031043595782, -0.0002443547165, 0.00034897302673], [3.0, -2.624171687e-05, -1.556581482e-05, 1.6969026934e-05, -1.83192537e-05], [4.0, -2.374676975e-06, 2.1175480977e-05, -1.848543004e-05, 1.8187992883e-05], [5.0, 0.0318793997, -0.0018898655, 0.0026784055, -0.0005866397982], [6.0, -0.0054903572, 0.00070266217971, -0.0008493354591, 0.00072392409916]]

TCE4_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0 K0=1e-4 FREQ=4.5e8 PHI=0.1);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
$FORM="18.14"; PageWidth = 1999;
CELL; CALC; RFSW; NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
TCE4_REF = [[1.167228799024e-07, -1.3902640181e-32, 9.993001301677e-06], [-0.22453092316574, 0.22100388888845, -0.0097094844334], [0.01887451928217, 293123.340749077, 0.01486882031539], [2.393100431013e-05, 9.489386682268e-06, -1.1506194385e-20, 5.50215731928e-21, -9.29015889364e-07, 1.879866717058e-05]]

TCE5_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
TCAVI TC = (L=0 K0=1e-4 FREQ=4.5e8 PHI=0.1 DX=2e-4 DY=-1e-4 ROTATE=0.05);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 TC CELL1 CELL1 TC 3*CELL1 RFC);
FFS;
USE RING;
$FORM="18.14"; PageWidth = 1999;
CELL; CALC; RFSW; NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
TCE5_REF = [[], [0.22196024916761, -0.00973499943373, -0.22196024916761], [9.86511050633349, 293123.2657098092, 0.01486881645914], [2.373334551664e-05, 9.391526100761e-06, -1.20837561599e-06, 6.35833066783e-07, 5.83515478921e-05, 1.88627621052e-05]]


NAMES = ['tct1', 'tct2', 'tct3', 'tct5', 'tct7', 'tct8', 'tce4', 'tce5']


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _nums(s):
    return [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:E[-+]?\d+)?', s)]


@pytest.mark.parametrize('name', [n for n in NAMES if n.startswith('tct')])
def test_tcav_tracking_matches_gs(name):
    ref = globals()[name.upper() + '_REF']
    rows = [[float(v) for v in l.split()] for l in _out(globals()[name.upper() + '_SRC']).splitlines()
            if re.match(r'^ +[0-9]+\.[0-9]+ ', l)]
    assert len(rows) == len(ref) == 6
    for r, q in zip(rows, ref):
        for a, b in zip(r, q):
            assert abs(a - b) <= 2e-9 * max(abs(b), 1e-3) + 1e-12, (name, r, q)


@pytest.mark.parametrize('name', ['tce4', 'tce5'])
def test_tcav_emittance_matches_gs(name):
    ref = globals()[name.upper() + '_REF']
    got = [_nums(l[4:]) for l in _out(globals()[name.upper() + '_SRC']).splitlines()
           if l[:4] in ('EMT ', 'TUN ', 'ALP ', 'COD ')]
    assert len(got) == len(ref) == 4
    for k, (r, q) in enumerate(zip(got, ref)):
        if name == 'tce5' and k == 0:
            continue                  # gs の EMT は NaN
        if k == 1:                    # チューンは符号・順序の違いがあるので絶対値で
            r, q = sorted(abs(v) for v in r), sorted(abs(v) for v in q)
        if name == 'tce5' and k == 2:  # この並びは縦が不安定で、運動量収縮因子の符号だけ流儀の差がある
            r, q = [abs(r[0])] + r[1:], [abs(q[0])] + q[1:]
        for a, b in zip(r, q):
            assert abs(a - b) <= 2e-8 * max(abs(b), 1e-3) + 1e-9, (name, k, r, q)
