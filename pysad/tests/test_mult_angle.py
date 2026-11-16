"""MULT (ANGLE≠0) (K64 の tmulta / tmultae)。参照は SAD K64 (gs) の出力。
ma1..ma11: TrackParticles (20 粒子、3 周)。mae*: Emittance[] (12 セル + RFC、MULT が曲げ磁石)。
参照の作り方: local/wake/genma.py genmae.py で .sad を作り gs で実行 (docs/26_mult_angle.md)。"""
import io
import re

import pytest

TRK_SRC = {}
TRK_REF = {}

TRK_SRC[1] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[1] = [[1.0, -0.0001030122734, 0.00016923166889, -7.129506147e-05, 0.00030505059211], [2.0, 9.2532124902e-06, 4.0601296375e-05, 7.6164737355e-06, 0.00011888889643], [3.0, -0.0016571804, 0.00041095988796, -8.616979181e-05, 3.2352525785e-05], [4.0, 0.0006946227098, -0.0001685142293, 3.3211065589e-05, -1.072403066e-05], [5.0, -0.0337774802, -0.0086183342, -0.0013917042, -0.0028509984], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[2] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K2=2 K3=30 K4=500);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[2] = [[1.0, -0.0002095549282, 0.00025305565255, -0.0001621140473, 0.00033822597407], [2.0, -4.391840455e-05, 7.335259117e-05, -2.864113052e-05, 0.0001384741268], [3.0, -0.0002342306502, 7.8392269111e-05, -2.971676964e-05, 2.16772288e-05], [4.0, 0.00033193031236, -8.722941965e-05, 2.1737332435e-05, -1.081365911e-05], [5.0, -0.033778151, -0.0086206815, -0.0013898159, -0.0028508718], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[3] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K0=0.001 K1=-0.1 K2=1 SK2=0.5 SK3=7);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[3] = [[1.0, -0.1692863207, -0.0080901414, -0.0087221827, -0.0080914144], [2.0, -0.0546019209, -0.0026090954, -0.0028007949, -0.002574597], [3.0, 0.00087725790117, -7.456594497e-05, 2.5260595941e-05, 3.1159012773e-05], [4.0, -1.344993012e-05, -1.002672265e-05, 9.7142168055e-06, -1.048595396e-05], [5.0, -0.0253084622, -0.008198426, -0.0009653710832, -0.0024262225], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[4] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 E1=0.5 E2=0.5 K1=0.1 K2=0.5);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[4] = [[1.0, -0.0001058552004, 0.00017103153856, -7.333790787e-05, 0.00030590052919], [2.0, 7.9237293308e-06, 4.1268967925e-05, 6.8474721907e-06, 0.00011938543354], [3.0, -0.0016100472, 0.00040002733932, -8.438196786e-05, 3.2075872062e-05], [4.0, 0.0006855221379, -0.0001664895081, 3.2935847615e-05, -1.073824679e-05], [5.0, -0.0337774738, -0.0086183937, -0.0013916532, -0.0028509994], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[5] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 E1=0.5 E2=0.5 K2=1 FB1=0.2 FB2=0.2 FRINGE=3);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[5] = [[1.0, -0.0002128563818, 0.00025552601811, -0.0001646909768, 0.00033908163578], [2.0, -4.527889421e-05, 7.4371572602e-05, -2.970347155e-05, 0.00013897035545], [3.0, -0.0002103848212, 7.2512011798e-05, -2.874951499e-05, 2.1428652216e-05], [4.0, 0.00032350778985, -8.529964072e-05, 2.1469237952e-05, -1.080875159e-05], [5.0, -0.0337780887, -0.0086207667, -0.0013897505, -0.0028508986], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[6] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5 DX=1e-3 DY=-5e-4 ROTATE=0.01 DROTATE=0.002);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[6] = [[1.0, 0.0077839293, 0.00057295401846, 0.00032210569478, 0.00070254159953], [2.0, 0.0016344408, 0.00012527133638, 8.8457800356e-05, 0.00020138521581], [3.0, 0.0418722835, 0.002581278, 0.0020915378, 0.0022015423], [4.0, -0.0096308112, -0.0006820706387, -0.0004842286184, -0.0005239224734], [5.0, -0.0344815314, -0.0086537568, -0.0014268304, -0.0028862811], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[7] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5 DISFRIN=1);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[7] = [[1.0, -0.0001029770818, 0.00016923563788, -7.129486308e-05, 0.00030505062321], [2.0, 9.2592490153e-06, 4.0602000302e-05, 7.6165206526e-06, 0.00011888890569], [3.0, -0.0016572262, 0.0004109507908, -8.617075304e-05, 3.2352263637e-05], [4.0, 0.00069462874947, -0.000168511826, 3.3211435251e-05, -1.072405212e-05], [5.0, -0.033777478, -0.0086183339, -0.0013917042, -0.0028509984], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[8] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[8] = [[1.0, -0.0001047540457, 0.00016914597388, -7.138405744e-05, 0.00030496162292], [2.0, 8.5703256796e-06, 4.0567601188e-05, 7.5817061291e-06, 0.00011885401091], [3.0, -0.0016572062, 0.00041095620023, -8.616987868e-05, 3.2352369834e-05], [4.0, 0.00069462909436, -0.0001685129688, 3.3211219972e-05, -1.072400334e-05], [5.0, -0.0337774252, -0.0086183314, -0.0013917014, -0.0028509958], [6.0, -0.0002779852581, 0.00095181352814, -0.000564503083, 0.000996022188]]

TRK_SRC[9] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K2=1 E1=0.3 E2=0.3);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[9] = [[1.0, -0.0002143463019, 0.00025430989587, -0.0001638621511, 0.00033861615184], [2.0, -4.57601365e-05, 7.3895421294e-05, -2.936458946e-05, 0.00013873485818], [3.0, -0.0002188484073, 7.4523383454e-05, -2.910899804e-05, 2.1515774855e-05], [4.0, 0.00032649080751, -8.597355991e-05, 2.156692054e-05, -1.081130874e-05], [5.0, -0.0337780151, -0.0086207054, -0.0013897803, -0.0028508663], [6.0, -0.0002786132698, 0.00095156257706, -0.0005652542225, 0.00099610534587]]

TRK_SRC[10] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.3 K2=20 K3=300 E1=0.4 E2=0.4);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[10] = [[1.0, -0.001228551, 0.0014510954, -0.0009349333098, 0.001608917], [2.0, -0.0004404283749, 0.00037160237893, -0.0002335900952, 0.00045636107899], [3.0, 0.00039196801412, -4.312007529e-05, -1.707256988e-06, 1.6306073499e-05], [4.0, 6.0003489208e-05, -3.59158889e-05, 1.2715935274e-05, -1.04748154e-05], [5.0, -0.0335363987, -0.0092441731, -0.0009896726993, -0.0034468084], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

TRK_SRC[11] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5 DROTATE=0.002 CHI2=0.01);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 8, 1, 3];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Exit;
'''
TRK_REF[11] = [[1.0, -3.507395333e-05, 0.00017479229512, -6.788665784e-05, 0.00030893259935], [2.0, 4.2733345252e-05, 4.3171017275e-05, 9.2861083519e-06, 0.00012076991418], [3.0, 0.0289537892, 0.0019417663, 0.0014440265, 0.0015626167], [4.0, -0.005606435, -0.000482894496, -0.000282212076, -0.0003249484667], [5.0, -0.0337861372, -0.0086190925, -0.001392137, -0.0028515056], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

EMT_SRC = {}
EMT_REF = {}

EMT_SRC['mae2'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 K2=2 E1=0.5 E2=0.5);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae2'] = {'EMT': [2.862888364842e-07, -1.0364854729e-36, 2.458179962058e-05], 'TUN': [-0.14318100905425, -0.36458354332279, -0.01137326354254], 'ALP': [0.02975741785275, 781632.856818632, 0.06215632196636]}

EMT_SRC['mae4'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K2=3 K3=100);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
RAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae4'] = {'EMT': [2.711338205361e-07, -3.8603638976e-33, 2.688253240139e-05], 'TUN': [-0.25135852508567, -0.19458417164485, -0.01173114780542], 'ALP': [0.03165964575467, 781632.8568186371, 0.06215632196636]}

EMT_SRC['mae5e'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 K2=2 DX=1e-4 DY=-1e-4 ROTATE=0.003);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae5e'] = {'EMT': [8.387254725965e-06, 2.7978257168e-11, 2.230817303781e-05], 'TUN': [0.00931275064945, 0.34880538941986, -0.0107275945063]}

EMT_SRC['mae5f'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K0=0.001 K1=0.1 K2=2 ROTATE=0.003);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae5f'] = {'EMT': [7.444472066229e-07, 5.18778497144e-12, 2.206511980223e-05], 'TUN': [-0.000392197613648, 0.33939550498311, -0.0107989856248]}

EMT_SRC['mae5g'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K0=0.001 K1=0.1 K2=2 DX=1e-4 DY=-1e-4);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae5g'] = {'EMT': [7.582227705978e-07, 1.38455631995e-12, 2.201236790978e-05], 'TUN': [-0.000913034971094, 0.33885562733179, -0.01080261982536]}

EMT_SRC['mae5b'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 K2=2 DX=1e-4 DY=-1e-4);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae5b'] = {'EMT': [7.777446522218e-06, 1.32519029231e-12, 2.227947338116e-05], 'TUN': [0.00926282463643, 0.3488723306832, -0.01073052820523]}

EMT_SRC['mae5c'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 K2=2 ROTATE=0.003);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae5c'] = {'EMT': [3.112898750671e-05, 5.0300387168e-12, 2.301533099789e-05], 'TUN': [0.00984652673696, 0.34941805662368, -0.01067249510875]}

EMT_SRC['mae7'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MULT B = (L=1.2 ANGLE=0.2617993877991494 K1=0.1 K2=2 SK2=1 SK3=40 K4=800 E1=0.4 E2=0.6);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW;
NORAD; NOFLUC;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["COD ", OrbitAtExit/.em];
Exit;
'''
EMT_REF['mae7'] = {'EMT': [2.851681307433e-07, 9.31174991393e-38, 2.445251647823e-05], 'TUN': [-0.13233554887425, -0.34906236468872, -0.01133630715066]}

CALC_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
BEND B1 = (L=1 ANGLE=0.1 E1=0.5 E2=0.5);
MULT M1 = (L=1 ANGLE=0.05 K1=0.1 K2=0.5);
MARK ST = ();
LINE BL = (ST QF DR M1 DR2 QD DR B1 DR2 QF DR M1 DR2 QD DR B1 DR);
FFS;
USE BL;
INS; NORAD; NOFLUC;
CALC;
Print["TUN ", Twiss["NX", "$$$"]/(2 Pi), " ", Twiss["NY", "$$$"]/(2 Pi)];
Print["BX  ", Twiss["BX", "$$$"], " ", Twiss["BY", "$$$"], " ", Twiss["EX", "$$$"]];
Print["ALP ", Twiss["AX","$$$"]];
Exit;
'''
CALC_REF = {'TUN': [0.7918745354374561, 0.26496578708000484], 'BX': [99.22268132795514, 2.0822579971314257, 0.7659249741321728], 'ALP': [-32.04020396747074]}


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _nums(s):
    return [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:E[-+]?\d+)?', s)]


@pytest.mark.parametrize('i', range(1, 12))
def test_mult_angle_tracking_matches_gs(i):
    ref = TRK_REF[i]
    rows = [[float(v) for v in l.split()] for l in _out(TRK_SRC[i]).splitlines() if re.match(r'^ +[0-9]+\.[0-9]+ ', l)]
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r, q):
            assert abs(a - b) <= 5e-9 * max(abs(b), 1e-3) + 1e-12, (i, r, q)


@pytest.mark.parametrize('name', sorted(EMT_REF))
def test_mult_angle_emittance_matches_gs(name):
    ref = EMT_REF[name]
    got = {l.split()[0]: _nums(l.split(None, 1)[1]) for l in _out(EMT_SRC[name]).splitlines() if l.split()[:1] in (['EMT'], ['TUN'], ['ALP'])}
    for key in ref:
        assert len(got[key]) == len(ref[key])
        for a, b in zip(got[key], ref[key]):
            assert abs(a - b) <= 2e-8 * abs(b) + 2e-10, (name, key, got[key], ref[key])


def test_mult_angle_calc_matches_gs():
    got = {l.split()[0]: _nums(l.split(None, 1)[1]) for l in _out(CALC_SRC).splitlines() if l.split()[:1] in (['TUN'], ['BX'], ['ALP'])}
    for key in ('TUN', 'ALP'):
        for a, b in zip(got[key], CALC_REF[key]):
            assert abs(a - b) <= 1e-10 * max(abs(b), 1.0), (key, a, b)
    assert abs(got['BX'][0] - CALC_REF['BX'][0]) <= 1e-9 * CALC_REF['BX'][0]


def test_mult_angle_unsupported_cases_raise():
    with pytest.raises(NotImplementedError):
        _out(TRK_SRC[1].replace('ANGLE=0.05', 'ANGLE=0.05 BZ=0.1'))
    with pytest.raises(NotImplementedError):
        _out(TRK_SRC[1].replace('ANGLE=0.05', 'ANGLE=0.05 SK0=0.001'))
