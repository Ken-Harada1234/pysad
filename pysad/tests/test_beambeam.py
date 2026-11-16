"""BEAMBEAM (強弱ビームのビーム・ビーム相互作用)。参照は SAD K64 (gs) の TrackParticles の出力 (全周を 5 周、200 粒子 (bs10 は 300) の 1・77・200 番目と成分ごとの和)。
BSTRL (bs1、bs10) は SeedRandom[7] と 330 粒子以下で乱数列まで gs と一致する (330 粒子を超えると gs の乱数列が変わるので、分布だけが一致する: docs/22_beambeam.md)。
b6・b9 は K64 の傾き回転の癖 (cost=sin) で座標がほぼ 0 に潰れる場合 (相対 1e-6 で比べる)。
lm* は LuminosityMonitor (シンボルの値と、利用者が定義した関数の呼び出し)。
参照の作り方: local/wake/genbb.py で b*.sad を作り gs で実行 (docs/22_beambeam.md)。"""
import io

import pytest

B1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B1_REF = [[-0.0593647929, -0.0162782798, -0.0135741312, 0.0140166342], [-0.0278419297, -0.0076419695, -0.0064008285, 0.0065502218], [-0.2126114516, -0.0583820774, -0.0496124923, 0.0491420083], [0.0115163184, 0.0031472695, 0.002671569, -0.0026240885], [-0.0278496838, -0.0090822584, -0.0081685477, 0.0172858657], [-0.0025272805, 0.00095234195851, 0.00030321105632, 0.00067351081554]]

B2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.01 BY=0.5 EMIX=1e-10 EMIY=1e-8 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B2_REF = [[34.6705494015, -1.0049315627, 3.0, 0.0227878733], [-6.144976365, -0.1904189234, 0.00011645421971, -0.0167704143], [8.4235472714, 0.8995574523, -3.0, -0.4713617394], [1.9067573512, 0.3915931281, 0.00051795120333, 0.1110478037], [50881780884.0, 1.4697704066, 71.0950573222, -0.0303627536], [-0.0025266206, 0.00095234659164, 0.00030321704997, 0.00067351348765]]

B3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.3 BY=0.3 EMIX=1e-9 EMIY=1e-9 SIGZ=0.005 DP=1e-3 SLICE=4 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B3_REF = [[-0.0309631723, -0.0104807173, 1.6975793152e-05, 0.0073276873], [-0.0145270012, -0.0048143458, -0.0001016663513, 0.0032824732], [-0.0999538201, -0.0344223158, -0.0018156017, 0.0213061092], [0.0052720857, 0.0017943107, 0.00014142308114, -0.001081155], [0.0296977639, -0.0087718379, -0.0078245694, 0.0175587948], [-0.0025214954, 0.00095235546971, 0.00030321336601, 0.00067355751789]]

B4_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10 XANGLE=0.01);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B4_REF = [[0.1557855929, -0.0076440594, -0.0055579865, 0.00047102854858], [0.0736797912, -0.003650009, -0.002682882, 0.00028974080109], [0.5791310138, -0.028861306, -0.0217863092, 0.0024258876], [-0.0315060721, 0.0015426996, 0.0011738634, -0.0001609747966], [0.0116183311, -0.0087196566, -0.007889852, 0.017624708], [-0.0025381168, 0.0009535368038, 0.00030451334231, 0.00067170740781]]

B5_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=3 NP=1e10 DX=2e-5 DY=1e-6 DPX=1e-6 DPY=2e-7 DZ=1e-3);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B5_REF = [[0.4543335715, -0.0158498751, -0.0066680239, 0.0140507456], [0.2086943564, -0.0074642334, -0.0031840001, 0.0065507629], [1.5277920323, -0.057505589, -0.0252607114, 0.0488106201], [-0.0814151016, 0.0031039626, 0.0013532742, -0.0026013851], [-0.0240762374, -0.0090673489, -0.0079127296, 0.0172900089], [-0.0025273878, 0.00095234182307, 0.0003032050602, 0.00067351583096]]

B6_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10 R3=0.3 R4=0.1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B6_REF = [[4.067937029e-43, 3.320938159e-44, -4.13630801e-44, -9.54350164e-44], [2.138478212e-43, 1.748315576e-44, -2.17634154e-44, -5.03967188e-44], [2.328113255e-43, 1.897292716e-44, -2.36181401e-44, -5.42573344e-44], [7.169084726e-44, 5.871665049e-45, -7.30347165e-45, -1.70216228e-44], [0.0409343903, -0.0086036571, -0.0078238444, 0.0176255839], [-0.0025286772, 0.00095235751588, 0.00030322788416, 0.00067346600861]]

B7_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 AX=0.2 AY=-0.1 EX=0.1 EPX=0.01 EY=0.02 EPY=0.005 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 AZ=2e-6 SLICE=5 NP=1e10 XANGLE=0.005);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B7_REF = [[-0.0003142522955, 4.4111394005e-05, 4.7212869566e-05, -8.404607101e-05], [-8.691110318e-05, 6.8574352328e-06, 9.7425722516e-06, -1.100562698e-05], [7.0155191127e-05, -3.224933537e-07, -5.211293377e-06, -3.89340278e-06], [6.5883692375e-06, 1.9258015067e-06, 6.4936534213e-07, -4.888960896e-06], [0.0409340929, -0.0086036584, -0.0078238457, 0.0176255789], [-0.0025297625, 0.00095293197227, 0.00030361798693, 0.00067224304932]]

B8_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=1 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B8_REF = [[-0.0605516918, -0.0162099774, -0.0149134108, 0.0151035952], [-0.0284822906, -0.007609972, -0.0070289612, 0.0070609439], [-0.2165080306, -0.0581355407, -0.0544648703, 0.0530827545], [0.01174301, 0.003133288, 0.0029434938, -0.0028443201], [-0.0318292541, -0.0090782053, -0.0082396278, 0.0172290845], [-0.0025271607, 0.00095234179568, 0.00030321364527, 0.00067351384491]]

B9_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10 R3=0.3 R4=0.1 XANGLE=0.01 DX=1e-5 DPX=1e-6);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B9_REF = [[0.0005631880143, 6.5986512242e-05, 6.0368165034e-05, -0.0001220864959], [-0.0002755481939, -9.835038897e-07, -1.025839663e-06, -2.148515506e-06], [0.0003421260802, 2.2326634662e-05, 2.0510763696e-05, -3.907301788e-05], [0.00013848990691, 9.0548921214e-06, 8.3111618787e-06, -1.584276117e-05], [0.0409332241, -0.0086036647, -0.0078238508, 0.0176255614], [-0.0025144851, 0.00095407734959, 0.00030449794838, 0.00067009423459]]

B10_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B10_REF = [[-0.0193785915, -0.0007448212705, -0.0022718982, 0.0035381829], [-0.0089365513, -0.0002588604163, -0.0011369536, 0.0015458746], [-0.0599348534, -0.0011279062, -0.0103541646, 0.0063687641], [0.0031946006, 4.2974530623e-05, 0.00056954068195, -0.0003073664159], [0.0343583401, -0.0086041638, -0.0078387922, 0.0176163457], [-0.0025281095, 0.00095233124239, 0.00030320549173, 0.00067351525339]]

B11_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.3 BY=0.3 EMIX=1e-9 EMIY=1e-9 SIGZ=0.005 DP=1e-3 SLICE=4 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-4*Sin[0.37*(i-1)], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Sin[0.23*(i-1)], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B11_REF = [[-0.0499961795, 0.0010142565, 0.0055289153, 0.0020062853], [-0.0224390393, 0.00063716493169, 0.0026626332, 0.00077000576282], [-0.1655768186, 0.0035906039, 0.0214155743, 0.0024264785], [0.0085370321, -0.000241255352, -0.0011523432, -8.822132868e-05], [0.0301025963, -0.0086057759, -0.0078876606, 0.0176229974], [-0.0025216096, 0.00095233951208, 0.00030320415186, 0.00067355535159]]

B12_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.01 BY=0.5 EMIX=1e-10 EMIY=1e-8 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 200;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[1e-4*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
B12_REF = [[19.169608846, 3.0, -0.1500968267, -2.89871122], [3.0617970321, 1.4198381543e-08, -0.5733782645, 0.1565172164], [1.9114065879, 3.0, 3.0, 3.0], [7.552215858, -9.430606232e-07, -0.1540801776, -0.0683795728], [20000994034.0, 4702.9457791029, 1.2346216882, -1.7213935202], [-0.0025281061, 0.00095232983861, 0.00030320554136, 0.00067351444613]]

BS1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11 BSTRL=1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
SeedRandom[7]; n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
BS1_REF = [[0.2959506715, 0.0940457721, -0.0220071013, 0.0320347856], [0.1305227646, 0.0416398614, -0.010681424, 0.0140058213], [1.0930239359, 0.3493474928, -0.0870986487, 0.0863071645], [-0.1309265383, -0.0409781569, 0.00498243, -0.004551162], [-2.6383404736, -0.0302175817, -0.0088913187, 0.016479832], [-0.0023689266, 0.00095338575677, 0.00030317224938, 0.00067507468652]]

BS10_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=2e12 BSTRL=1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
SeedRandom[7]; n = 300;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 300]]], {k, 6}];
Exit;
'''
BS10_REF = [[0.5682415534, 0.1070195067, 0.0135634676, -0.0636277397], [0.2828189209, 0.0371478732, 0.0039973355, -0.036108932], [0.9905491493, 0.4372652902, 0.4282022715, 0.3326821802], [-0.390193171, -0.3965759922, -0.2533304653, -0.1764339713], [-72.7959707891, -0.301975589, -0.2183650563, -0.172596598], [0.00044753270893, 0.00096261102619, 0.00030632933856, 0.000286788324]]

LM2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
cnt = 0; tot = 0; LuminosityMonitor[x_] := (cnt = cnt + 1; tot = tot + x); n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 5, 1];
$FORM = "16.10";
res = r[[2]];
Print["LUM ", LuminosityMonitor, " ", cnt, " ", tot];
Exit;
'''
LM2_REF = [441603719.01268, 1.0, 441603719.01268]

LM3_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=1e10);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
cnt = 0; tot = 0; LuminosityMonitor[x_] := (cnt = cnt + 1; tot = tot + x); n = 200;
x0  = Table[1e-5*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-7*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 3];
$FORM = "16.10";
res = r[[2]];
Print["LUM ", LuminosityMonitor, " ", cnt, " ", tot];
Exit;
'''
LM3_REF = [9465679.7972342, 3.0, 1999645363.9777]

BSR1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11 BSTRL=1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC;
SeedRandom[7]; n = 200;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 200]]], {k, 6}];
Exit;
'''
BSR1_REF = [[0.29581015, 0.0940027632, -0.0220064989, 0.0320343596], [0.1305231667, 0.0416370195, -0.0106812585, 0.0140057452], [1.0929929334, 0.3493268995, -0.0870978881, 0.0863065596], [-0.13031717, -0.0407959886, 0.0049801738, -0.0045489991], [-2.6311506015, -0.0301688015, -0.0088911758, 0.0164799662], [-0.004304739, 0.00091651248184, 0.0003023490364, 0.0006718365524]]

BSR10_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5) QX = (L=0.3 K1=0.4 ROTATE=0.2);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=2e12 BSTRL=1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QX DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC;
SeedRandom[7]; n = 300;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 5];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 77]], " ", res[[k, 300]]], {k, 6}];
Exit;
'''
BSR10_REF = [[0.5608098794, 0.1083793762, 0.0165519761, -0.0591518809], [0.2783959604, 0.0378144501, 0.0054087312, -0.0336358928], [0.9669956579, 0.4377312115, 0.4304791259, 0.3403594942], [-0.3909532711, -0.3992085899, -0.2576321302, -0.1803382027], [-73.0509463029, -0.3039673184, -0.2203231293, -0.1727275758], [-0.0364753375, 0.00087862125997, 6.2022946686e-05, 0.00016603466889]]

RF0R2_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC;
SeedRandom[7];
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 3, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Print["RND ", Random[], " ", Random[]]; Exit;
'''
RF0R2_REF = [0.5861273285, 0.4596477298]

RF0R_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC;
SeedRandom[7];
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 6}];
Print["RND ", Random[], " ", Random[]]; Exit;
'''
RF0R_REF = [0.6510187115, 0.0227274296]

PH4_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
BEAMBEAM BB = (BX=0.5 BY=0.01 EMIX=1e-8 EMIY=1e-10 SIGZ=0.005 DP=1e-3 SLICE=5 NP=3e11 BSTRL=1);
MARK ST = ();
LINE BL = (ST QF DR BB DR2 QD DR);
FFS;
USE BL;
TRPT; INS; RAD; FLUC; PHOTONS;
SeedRandom[7];
n = 20;
x0  = Table[1e-4*Sin[0.37*i], {i, n}];
px0 = Table[1e-6*Cos[0.11*i], {i, n}];
y0  = Table[3e-6*Cos[0.23*i], {i, n}];
py0 = Table[1e-7*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 7, 1];
$FORM = "16.10";
Print["NPH ", Length[PhotonList]];
Do[Print["PA ", PhotonList[[k, 1]], " ", PhotonList[[k, 2]], " ", PhotonList[[k, 3]], " ", PhotonList[[k, 4]], " ", PhotonList[[k, 5]], " ", PhotonList[[k, 6]]], {k, Length[PhotonList]}];
Do[Print["PB ", PhotonList[[k, 7]], " ", PhotonList[[k, 8]], " ", PhotonList[[k, 9]], " ", PhotonList[[k, 10]], " ", PhotonList[[k, 11]], " ", PhotonList[[k, 12]]], {k, Length[PhotonList]}];
Exit;
'''
PH4_REF = [[152.6092405739, 1.2969462079, -3.153182883e-05, 1.1497769242e-06, 0.9999991813, 0.001244315], [972.0408713847, 1.3035533676, 3.655127411e-06, 2.7360210421e-06, 0.9999997679, -0.0004837313926], [15762.015132928, 1.3003177906, 1.530924221e-05, -1.96852792e-06, 0.9999994472, -0.0009664024202], [53817.967154491, 1.3018302454, 2.1657362632e-05, -3.422220412e-06, 0.9999978205, -0.0018430343], [1459.3134251685, 1.299284794, 3.1358088559e-05, 1.701216592e-06, 0.9999947745, -0.0012079058], [1423.2880821511, 1.2990650749, 2.1287808542e-05, 1.4877166745e-06, 0.9999948937, -0.0008491842532], [92.2089874929, 1.2959927091, -1.174004115e-05, 2.7433420528e-06, 0.9999985203, 0.00024699664132], [0.00029844849173, 0.0, 0.0, 1.0, 5.0, 4.0], [-0.0004797727685, 0.0, 0.0, 1.0, 9.0, 4.0], [0.00041438493586, 0.0, 0.0, 1.0, 10.0, 4.0], [0.00098088864607, 0.0, 0.0, 1.0, 11.0, 4.0], [-0.0029986707, 0.0, 0.0, 1.0, 13.0, 4.0], [-0.0030808347, 0.0, 0.0, 1.0, 15.0, 4.0], [-0.0017024918, 0.0, 0.0, 1.0, 18.0, 4.0]]

CASES = ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'bs1', 'bs10', 'bsr1', 'bsr10']


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _run(src):
    rows = []
    for l in _out(src).splitlines():
        t = l.split()
        if len(t) == 5 and t[0].replace('.', '').isdigit():
            rows.append([float(v) for v in t[1:]])
    return rows


@pytest.mark.parametrize('name', CASES)
def test_beambeam_matches_gs(name):
    src = globals()[name.upper() + '_SRC']
    ref = globals()[name.upper() + '_REF']
    rows = _run(src)
    assert len(rows) == len(ref) == 6
    scale = max(abs(v) for r in ref for v in r)
    for r, q in zip(rows, ref):
        for a, b in zip(r, q):
            tol = 1e-6 * scale if name in ('b6', 'b9') else (3e-8 if name.startswith('bs') else 2e-9) * max(abs(b), 1e-3 * scale)
            assert abs(a - b) <= tol, (name, r, q)


@pytest.mark.parametrize('name', ['lm2', 'lm3'])
def test_luminosity_monitor_matches_gs(name):
    ref = globals()[name.upper() + '_REF']
    line = [l for l in _out(globals()[name.upper() + '_SRC']).splitlines() if l.startswith('LUM')][0]
    got = [float(v) for v in line.split()[1:]]
    assert len(got) == len(ref)
    for a, b in zip(got, ref):
        assert abs(a - b) <= 1e-9 * abs(b)


def test_beambeam_has_an_effect():
    off = _run(B1_SRC.replace('SLICE=5', 'SLICE=0'))
    on = _run(B1_SRC)
    assert abs(on[0][1] - off[0][1]) > 1e-8


def test_beambeam_bstrl_too_large_stops():
    # K64 は rnbs>0.5 で stop する。pysad は SAD 風のエラーメッセージを出して結果を作らない
    assert 'Nbs is too large' in _out(B1_SRC.replace('NP=1e10', 'NP=1e13 BSTRL=1'))


@pytest.mark.parametrize('name', ['rf0r2', 'rf0r'])
def test_random_state_after_rad_fluc_matches_gs(name):
    """TRPT の放射の分割 (tffs.f ndivrad は TRPT で 10 倍) が gs と同じなので、追跡のあとの乱数の状態が一致する"""
    line = [l for l in _out(globals()[name.upper() + '_SRC']).splitlines() if l.startswith('RND')][0]
    got = [float(v) for v in line.split()[1:]]
    for a, b in zip(got, globals()[name.upper() + '_REF']):
        assert abs(a - b) < 2e-10


def test_photon_list_matches_gs():
    """PHOTONS: BEAMBEAM の光子 (PhotonList)。他の要素の光子 (tphrec) は含まない。この台本は他の光子が出ない"""
    out = _out(PH4_SRC).splitlines()
    got = [[float(v) for v in l.split()[1:]] for l in out if l.startswith('PA ') or l.startswith('PB ')]
    ref = PH4_REF
    assert len(got) == len(ref)
    for r, q in zip(got, ref):
        for a, b in zip(r, q):
            assert abs(a - b) <= 5e-9 * max(abs(b), 1e-12) + 1e-14


def test_photons_require_rad_fluc():
    assert 'Track::photons' in _out(PH4_SRC.replace('RAD; FLUC; PHOTONS;', 'NORAD; NOFLUC; PHOTONS;'))
