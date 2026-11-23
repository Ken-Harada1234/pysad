"""WakeFunction (txwake.f) の追跡。参照値は SAD K64 の gs (WSL) で同じ台本を走らせた出力 (10 桁)。
各台本は粒子 200〜400 個を TrackParticles で通し、成分ごとの総和と 3 粒子の値を表示する。
gs の実測で分かった K64 の癖: 追跡では WakeFunction の表の最後の 1 行が使われない (twxiwp の lwl0 = 行数 - 1)。"""
import io

import pytest

from pysad.ffscmd import run_sad_string


W3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
W3_REF = [[0.3010093897, 0.0010546522, 0.0010599768, 0.0014788142], [-0.1735844718, -0.0006259058275, -0.0006930857471, -0.0014403482], [-0.1883379886, -0.0032387128, -0.0016085325, 0.0011916239], [-0.1105499065, -0.0019962087, -0.0009386596886, 0.00070700790872], [0.0584583581, -0.0187934151, -0.0179981392, 0.0294898416], [-0.0038440346, 0.00095214552167, 0.00030301340003, -9.313428675e-05], [400.0, 1.0, 1.0, 1.0]]

W4_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 3;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+0.5*Mod[i,3]+1e-4*i, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
W4_REF = [[0.301009548, 0.0010546549, 0.0010599753, 0.0014788256], [-0.173584509, -0.0006259073492, -0.0006930844868, -0.0014403544], [-0.1883372163, -0.0032387046, -0.0016085261, 0.0011916232], [-0.1105495177, -0.0019962044, -0.0009386566519, 0.00070700778976], [208.058458363, 0.5012065847, 1.002001861, 0.5494898408], [-0.0038441499, 0.00095214399825, 0.00030301373267, -9.313406316e-05], [400.0, 1.0, 1.0, 1.0]]

W5_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
LOSSAMPL = 0.0025; LOSSDZ = 100;
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
W5_REF = [[0.3010093897, 0.0010546522, 0.0010599768, 0.0014788142], [-0.1735844718, -0.0006259058275, -0.0006930857471, -0.0014403482], [-0.1883379886, -0.0032387128, -0.0016085325, 0.0011916239], [-0.1105499065, -0.0019962087, -0.0009386596886, 0.00070700790872], [0.0584583581, -0.0187934151, -0.0179981392, 0.0294898416], [-0.0038440346, 0.00095214552167, 0.00030301340003, -9.313428675e-05], [293.0, 0.0, 1.0, 1.0]]

C3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
CAVI RF = (VOLT=0 HARM=100 PHI=0.2 L=0.4);
MARK ST = ();
LINE BL = (ST QF DR RF DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "RF"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "RF"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
C3_REF = [[-0.0015223686, 3.9456420645e-05, -2.014326793e-05, 0.00039572105192], [-0.0017476653, -9.267674937e-05, 5.725409155e-05, 0.00063689586325], [0.0088282277, 0.0026935059, 0.00078326872127, -0.0012861579], [-0.0022680207, -0.0006915471325, -0.0002457308629, 0.00031948495619], [0.0408319596, -0.0086042681, -0.0078238836, 0.0176251065], [-0.0025336477, 0.00095234971465, 0.00030321764289, 0.00067342171247]]


C4N_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
CAVI RF = (L=0.6 VOLT=3e7 HARM=100 PHI=0.3);
MARK ST = ();
LINE BL = (ST QF DR RF DR QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "RF"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "RF"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
C4N_REF = [[-0.0019153466, 1.377880065e-05, -4.702804062e-06, 0.00049584444517], [-0.0018691845, -0.0001044320116, 6.4322969485e-05, 0.00069059125454], [0.0091373329, 0.0028401338, 0.00081817583455, -0.0013406832], [-0.0023939799, -0.0007719861468, -0.0002649275223, 0.00034394429164], [0.0408235258, -0.0086043375, -0.0078238898, 0.0176250752], [0.6475981853, -0.0115375007, -0.0121031195, 0.0187046497]]

C4_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
CAVI RF = (L=0.6 VOLT=3e7 HARM=100 PHI=0.3);
MARK ST = ();
LINE BL = (ST QF DR RF DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "RF"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "RF"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
C4_REF = [[-0.0019151256, 1.3780652137e-05, -4.70044305e-06, 0.00049584452131], [-0.0018689696, -0.0001044302358, 6.4325237423e-05, 0.00069059128618], [0.0091373694, 0.0028401321, 0.00081817127905, -0.0013406833], [-0.0023939653, -0.0007719863618, -0.0002649281105, 0.00034394430881], [0.040823532, -0.0086043372, -0.0078238896, 0.017625075], [0.6475930854, -0.0115374845, -0.0121031077, 0.0187045638]]

C6_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
CAVI RF = (VOLT=3e7 HARM=100 PHI=0.3);
MARK ST = ();
LINE BL = (ST QF DR RF DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "RF"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "RF"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
C6_REF = [[-0.0012371215, 8.8223490382e-05, -4.936736716e-05, 0.0001840134457], [-0.0016048197, -6.995351257e-05, 4.3635209122e-05, 0.00053879846798], [0.0079643084, 0.0024565432, 0.00072877394555, -0.0011619945], [-0.0017993047, -0.0005631016029, -0.0002162531816, 0.00025217612226], [0.0408483543, -0.0086041811, -0.0078238793, 0.0176251879], [0.7707364411, -0.0110691546, -0.0121128358, 0.0122577872]]

S1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
SOL SA = (L=0 BZ=0.5 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR SA DR2 QX DR2 SB DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QX"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QX"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 12, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
S1_REF = [[0.0012123278, 0.0012611047, 0.00079751747932, 0.0010946808], [0.00054652364465, 0.00053429903322, 0.0003962336981, 0.0005728055456], [0.0126696472, 0.0039099742, 0.0010574483, -0.0019367024], [-0.0029067829, -0.0008752114379, -0.0002864050398, 0.0004068641942], [0.0405271827, -0.0086066187, -0.0078241806, 0.017624435], [-0.0025336477, 0.00095234971495, 0.00030321764227, 0.0006734217124]]

S3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
SOL SA = (L=0 BZ=0.5 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR SA DR2 QX DR2 SB DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "SB"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "SB"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 12, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
S3_REF = [[0.0012122326, 0.0012611018, 0.00079751344692, 0.0010946808], [0.00054650426097, 0.00053429751917, 0.00039623161649, 0.00057280553033], [0.0126696799, 0.0039099739, 0.0010574475, -0.0019367024], [-0.0029067615, -0.0008752119985, -0.0002864065832, 0.00040686416683], [0.0405271827, -0.0086066187, -0.0078241806, 0.0176244351], [-0.0025336477, 0.00095234971681, 0.00030321764061, 0.00067342171208]]


def _rows(src):
    """台本の表示 (Print の行) を数値の行にして返す"""
    out = io.StringIO()
    run_sad_string(src, out=out)
    rows = []
    for line in out.getvalue().splitlines():
        t = line.split()
        if len(t) >= 5:
            try:
                rows.append([float(v) for v in t])
            except ValueError:
                pass
    return rows


@pytest.mark.parametrize('src,ref', [(W3_SRC, W3_REF), (W4_SRC, W4_REF), (W5_SRC, W5_REF), (C3_SRC, C3_REF)],
                         ids=['two_wake_elements_quad_offset', 'nbunch3', 'with_losses', 'cavity_zero_volt'])
def test_wake_matches_sad(src, ref):
    rows = _rows(src)
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        assert r[0] == int(r[0])
        for a, b in zip(r[1:], q):
            assert a == pytest.approx(b, rel=2e-8, abs=1e-12)


def test_wake_has_an_effect():
    src = W3_SRC.replace(' LWAKE; TWAKE;', '')
    on = _rows(W3_SRC)
    off = _rows(src)
    assert abs(on[5][1] - off[5][1]) > 1e-6          # dp の総和が変わる


def test_wake_last_table_row_is_unused():
    # 最後の 1 行だけを変えても追跡の結果は変わらない (gs で確認した K64 の癖)
    a = W3_SRC
    b = W3_SRC.replace('{0.2, -1e9}};', '{0.2, -7e9}};', 1)
    assert a != b
    ra, rb = _rows(a), _rows(b)
    for x, y in zip(ra, rb):
        assert x[1:] == y[1:]



@pytest.mark.parametrize('src,ref', [(C4N_SRC, C4N_REF), (C4_SRC, C4_REF), (C6_SRC, C6_REF), (S1_SRC, S1_REF), (S3_SRC, S3_REF)],
                         ids=['accel_cavity_no_wake', 'accel_cavity_thick_wake', 'accel_cavity_thin_wake', 'sol_region_inner_wake', 'sol_region_end_wake'])
def test_wake_trpt_acceleration_and_sol_matches_sad(src, ref):
    # TRPT の加速 (vnominal≠0) のある CAVI と、SOL 領域の内側・終端の要素のウェイク (参照は gs)
    rows = _rows(src)
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r[1:], q):
            assert a == pytest.approx(b, rel=2e-8, abs=1e-12)


R1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4 ROTATE=0.3) QD = (L=0.3 K1=-0.5 DY=-1e-4 ROTATE=-0.2);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
R1_REF = [[0.4126550023, 0.0091994766, 0.0032311081, -8.929681209e-05], [0.0248064827, 0.00039127845364, 1.4584347331e-05, 0.00035819276566], [0.3531458151, 0.0094434545, 0.0033172404, 0.00020923862286], [0.3688684605, 0.00931313, 0.0034251031, -3.053179818e-05], [0.0723405185, -0.019379521, -0.0181220364, 0.0299030235], [-0.0038440337, 0.00095214562619, 0.00030301343399, -9.313424894e-05], [400.0, 1.0, 1.0, 1.0]]

T2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.003*((i-Mod[i,5])/5)-0.1, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
T2_REF = [[0.3009466864, 0.0010544467, 0.0010598206, 0.0014786331], [-0.1735524567, -0.0006257965762, -0.0006930027504, -0.0014402504], [-0.1883076961, -0.0032386311, -0.0016084672, 0.0011917009], [-0.1105346616, -0.0019961678, -0.0009386271816, 0.00070704693007], [7.5841942275, -0.1001897428, -0.055174284, 0.1396235788], [-0.0038022738, 0.00095230144826, 0.00030310773818, -9.298532657e-05], [400.0, 1.0, 1.0, 1.0]]



@pytest.mark.parametrize('src,ref', [(R1_SRC, R1_REF), (T2_SRC, T2_REF)],
                         ids=['rotate_on_wake_elements', 'equal_z_ordering'])
def test_wake_rotate_and_equal_z_match_sad(src, ref):
    # ROTATE≠0 の要素 (K64 は ROTATE の値を dy に入れる癖)、同じ z の粒子の並び (spkick.f の spsort) (参照は gs)
    rows = _rows(src)
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r[1:], q):
            assert a == pytest.approx(b, rel=2e-8, abs=1e-12)

MR1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT MR = (L=0.5 K1=0.2 K2=3 VOLT=3e6 HARM=100 PHI=0.3 DISFRIN=1);
MARK ST = ();
LINE BL = (ST QF DR MR DR QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
MR1_REF = [[-0.0018612993, -0.0001932896097, 0.00013222425369, 0.0013273154], [-0.0009779328926, -0.0001601950269, 0.0001033762944, 0.00077355394935], [0.0089219549, 0.0027236155, 0.00076620012884, -0.0013037295], [-0.0019569588, -0.000598940843, -0.000213177464, 0.00027622569679], [0.0407194827, -0.0086049298, -0.0078239512, 0.0176245954], [0.0640869972, -0.0002890051787, -0.000935934981, 0.0023868766]]

MR2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT MR = (L=0 K1=0.2 K2=3 VOLT=3e6 HARM=100 PHI=0.3);
MARK ST = ();
LINE BL = (ST QF DR MR DR QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
MR2_REF = [[-0.0016344494, -0.0001060635513, 7.6480646775e-05, 0.00099108648607], [-0.001005448, -0.0001369601306, 8.7878091518e-05, 0.00068058128215], [0.0080957688, 0.0024688751, 0.00070283454283, -0.0011792269], [-0.0013236828, -0.0004062279537, -0.0001648872051, 0.0001829554826], [0.0407673614, -0.0086045922, -0.0078239239, 0.0176247982], [0.0744353095, -0.0002441715607, -0.0009325688992, 0.0018265037]]

MR3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT MR = (L=0.5 K1=0.2 K2=3 K3=40 DX=1e-4 DY=-2e-4 ROTATE=0.1 VOLT=3e6 DVOLT=2e5 FREQ=5.2e8 PHI=0.3 DPHI=0.1);
MARK ST = ();
LINE BL = (ST QF DR MR DR QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
MR3_REF = [[0.0173352494, 0.00016242173408, 0.00030526663859, 0.0012957547], [0.0087867749, 2.1136108856e-05, 0.00019168796755, 0.00075633842293], [0.0212751338, 0.0027833818, 0.00081571913985, -0.0012572188], [0.00024678577893, -0.0005882351516, -0.0002042978955, 0.0002854585419], [0.0407200851, -0.0086049395, -0.007823983, 0.0176246297], [-0.0434694856, 0.00058202173717, -5.332842923e-05, 0.00080688475472]]

MW1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = ();
MULT WK2 = (L=0.5 K1=0.2 K2=3 VOLT=3e6 HARM=100 PHI=0.3);
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
MW1_REF = [[0.1648091466, 0.00024721423792, 0.0002009412094, 0.0010135122], [-0.1176217729, -0.0003430693563, -0.0004015999008, -0.0011899103], [-0.2526728274, -0.0011251153, -0.0012142055, 5.8702959155e-05], [-0.1416880319, -0.0008726210559, -0.0007175318279, 0.00011104640306], [0.0639658738, -0.0187729752, -0.0179811773, 0.0295208927], [0.0544629606, -0.0002539944104, -0.000884903714, 0.002169646], [400.0, 1.0, 1.0, 1.0]]

MW3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = ();
MULT WK2 = (L=0 K1=0.2 VOLT=3e6 HARM=100 PHI=0.3);
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
MW3_REF = [[0.1971108497, 0.0003197148672, 0.00031201235661, 0.0012507468], [-0.1139928967, -0.0003221929106, -0.0003813030158, -0.0011564243], [-0.2306251957, -0.000324317584, -0.0008932602947, -0.0002708551529], [-0.1319448908, -0.0005207491789, -0.000576398652, -3.251759368e-05], [0.0676794473, -0.0187611913, -0.0179667515, 0.0295488026], [0.0588121323, -0.000269992832, -0.0009041918662, 0.0021877403], [400.0, 1.0, 1.0, 1.0]]

MW4_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = ();
MULT WK2 = (L=0.5 K1=0.2 K2=3 DX=1e-3 ROTATE=0.2 VOLT=3e6 HARM=100 PHI=0.3 DISFRIN=0);
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
MW4_REF = [[0.3178499409, 0.00094894154422, 0.00067639921742, 0.001310618], [-0.2075877282, -0.0007561644951, -0.0006814838805, -0.0013640245], [-0.2945604109, -0.0014030976, -0.0012356577, 0.00090438401266], [-0.1626292512, -0.0010116886, -0.0007281509715, 0.0005353444986], [0.0506041728, -0.0188340831, -0.0180224597, 0.02949471], [0.0544592747, -0.0002539975392, -0.0008849080305, 0.0021696406], [400.0, 1.0, 1.0, 1.0]]

MW2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = ();
MULT WK2 = (L=0.5 K1=0.2 K2=3);
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QF.1"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "WK2"] = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
WakeFunction[Transverse, "WK2"] = {{-0.1, 1e12}, {0.0, 3e12}, {0.1, 0}, {0.3, 0}};
n = 400;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.02, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 17, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 400]]], {k, 7}];
Exit;
'''
MW2_REF = [[0.1486309372, 0.00058352325161, 0.00053279111286, 0.00037423554244], [-0.1131302909, -0.0004362075244, -0.0004935594189, -0.0010129978], [-0.2526180059, -0.0011269544, -0.0012138569, 5.7015066081e-05], [-0.1416565813, -0.0008740001221, -0.0007177932375, 0.00010995431739], [0.0640772352, -0.0187748873, -0.0179830377, 0.0295244198], [-0.0038440346, 0.00095214552166, 0.00030301339994, -9.313428674e-05], [400.0, 1.0, 1.0, 1.0]]

C10_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
CAVI RF = (L=0.6 VOLT=3e7 HARM=100 PHI=0.3 DPHI=0.2);
MARK ST = ();
LINE BL = (ST QF DR RF DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "RF"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "RF"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
C10_REF = [[-0.0018996153, 1.3779002862e-05, -4.700448705e-06, 0.00049526879897], [-0.0018629754, -0.0001044561101, 6.4325184687e-05, 0.00069090301222], [0.0090743252, 0.0028396486, 0.00081817171583, -0.0013394255], [-0.0023389729, -0.0007715582871, -0.0002649284583, 0.00034289983893], [0.0408233262, -0.0086043368, -0.0078238896, 0.0176250761], [0.3762506756, -0.0111266318, -0.0121043686, 0.0208928672]]



@pytest.mark.parametrize('src,ref', [(MR1_SRC, MR1_REF), (MR2_SRC, MR2_REF), (MR3_SRC, MR3_REF), (MW1_SRC, MW1_REF), (MW3_SRC, MW3_REF), (MW4_SRC, MW4_REF), (MW2_SRC, MW2_REF), (C10_SRC, C10_REF)],
                         ids=['mult_rf_trpt', 'mult_rf_thin', 'mult_rf_all_keys', 'mult_rf_wake', 'mult_thin_rf_wake', 'mult_rf_wake_offset_rotate_fringe', 'mult_wake_without_rf', 'cavity_dphi_wake'])
def test_mult_rf_and_wake_match_sad(src, ref):
    # RF のある MULT (tmulti.f tmultiacc: TRPT の加速、薄い MULT、DVOLT/DPHI/FREQ/DX/ROTATE/フリンジ)、MULT の wake、CAVI の DPHI (参照は gs)
    rows = _rows(src)
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r[1:], q):
            assert a == pytest.approx(b, rel=2e-8, abs=1e-12)

SS1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
SOL SA = (L=0 BZ=0.5 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR SA DR2 QX DR2 SB DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "SA"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "SA"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 12, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
SS1_REF = [[0.0012122101, 0.0012611074, 0.0007975203238, 0.0010946809], [0.00054644907343, 0.00053430015331, 0.0003962348652, 0.00057280555226], [0.0126696649, 0.003909974, 0.001057448, -0.0019367024], [-0.0029067893, -0.0008752110841, -0.0002864041446, 0.00040686420813], [0.0405271826, -0.0086066187, -0.0078241806, 0.017624435], [-0.0025336477, 0.00095234971462, 0.00030321764298, 0.00067342171248]]

SS2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
SOL SA = (L=0 BZ=0.5 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR SA DR2 QX DR2 SB DR QD DR);
FFS;
USE BL;
TRPT; INS; LWAKE; TWAKE; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QX"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Transverse, "QX"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
WakeFunction[Longitudinal, "SA"] = {{-0.1, 1e9}, {-0.02, 2e9}, {0, 3e8}, {0.02, -1e9}, {0.1, -1e9}, {0.2, -1e9}};
WakeFunction[Transverse, "SA"] = {{0, 2e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-4*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-4*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 12, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
SS2_REF = [[0.0012129138, 0.0012611108, 0.00079752380417, 0.001094681], [0.00054677751042, 0.00053430173779, 0.00039623643791, 0.00057280556103], [0.0126691863, 0.0039099766, 0.0010574413, -0.0019367025], [-0.0029067504, -0.000875211514, -0.0002864045483, 0.00040686430653], [0.0405271813, -0.0086066187, -0.0078241806, 0.017624435], [-0.0025409993, 0.00095237212904, 0.00030323567641, 0.00067327509125]]



@pytest.mark.parametrize('src,ref', [(SS1_SRC, SS1_REF), (SS2_SRC, SS2_REF)],
                         ids=['sol_region_start_wake', 'sol_region_start_and_inner_wake'])
def test_wake_on_sol_region_start_matches_patched_sad(src, ref):
    # 元の SAD K64 は領域の始点の SOL にウェイクがあると SIGSEGV で落ちる (tturn.f が tsol の後で、領域内で進んだ nwak を使う)。
    # 参照は修正した SAD (docs/history/sad_patch/wake_sigsegv_fixes.patch) の出力
    rows = _rows(src)
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r[1:], q):
            assert a == pytest.approx(b, rel=2e-8, abs=1e-12)
