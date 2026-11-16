"""UND (アンジュレータ、K64 の undulator.f)。参照は SAD K64 (gs) の出力。
und1: TrackParticles (20 粒子、3 周)。undo: INS の Twiss と LINE["G.."] (UND は恒等)。unde: Emittance[] (UND は恒等)。
参照の作り方: local/wake/und1.sad などを gs で実行 (docs/25_und.md)。"""
import io
import re

import pytest

UND1_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
UND U1 = (L=1 KX=0.1 QY=0.2 POLE=10 SLICE=4 PHI=30);
MARK ST = ();
LINE BL = (ST QF DR U1 DR2 QD DR);
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
UND1_REF = [[1.0, -0.0001432008549, -1.843720977e-05, -2.943880252e-06, 5.2543413725e-05], [2.0, -2.735941901e-05, -3.058103703e-07, 1.4306230522e-05, 6.0850545755e-05], [3.0, -0.0002342346384, 7.7939650991e-05, -2.980626548e-05, 2.1667617463e-05], [4.0, 0.00033180887469, -8.714027285e-05, 2.1748954178e-05, -1.081789859e-05], [5.0, -0.0337865057, -0.0086036814, -0.0014049668, -0.0028450065], [6.0, -0.0002676093611, 0.00095233356989, -0.0005639850576, 0.00099654209702]]

UNDO_SRC = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
UND CO = (L=1 KX=0.1 QY=0.2 POLE=10 SLICE=4 PHI=30);
MARK ST = ();
LINE BL = (ST QF DR B DR CO DR2 QD DR B DR);
FFS;
USE BL;
INS; CALC;
$FORM = "16.10";
Do[Print[k, " ", Twiss["BX",k], " ", Twiss["BY",k], " ", Twiss["NX",k], " ", Twiss["R1",k]," ", Twiss["R3",k], " ", Twiss["EX",k], " ", Twiss["DX",k], " ", Twiss["DY",k]], {k, 1, 12}];
Do[Print[k, " ", LINE["GX",k], " ", LINE["GY",k], " ", LINE["GZ",k], " ", LINE["GCHI1",k], " ", LINE["GCHI2",k], " ", LINE["GCHI3",k]], {k, 1, 12}];
Exit;
'''
UNDO_REF = [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [3.0, 0.915272053, 1.2865938391, 0.3091566337, 0.0, 0.0, 0.0], [4.0, 1.5541799398, 4.8836921619, 1.3035791111, 0.0, 0.0, 0.0], [5.0, 4.5030303087, 11.6196981644, 1.6905389177, 0.0, 0.0, 0.0499583472], [6.0, 9.718551186, 21.4811414066, 1.842284062, 0.0, 0.0, 0.1498167432], [7.0, 9.718551186, 21.4811414066, 1.842284062, 0.0, 0.0, 0.1498167432], [8.0, 13.1834919993, 27.5867030027, 1.8864711465, 0.0, 0.0, 0.1997459413], [9.0, 17.8522886158, 27.3055787851, 1.9065201022, 9.412163796e-18, 6.118641707e-17, 0.2456271814], [10.0, 51.9414209648, 14.5667198989, 1.9393654774, 7.059858086e-17, 6.118641707e-17, 0.4554576318], [11.0, 103.1184917463, 5.8886371158, 1.953007043, 1.297956874e-16, 6.018119069e-17, 0.7128488472], [12.0, 170.6325191936, 1.2886845764, 1.9605458901, 1.89976878e-16, 6.018119069e-17, 1.0175004672], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [3.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0], [4.0, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0], [5.0, 2.2983341665, 0.0499583472, 0.0, 0.1, 0.0, 0.0], [6.0, 3.2933383317, 0.1497917639, 0.0, 0.1, 0.0, 0.0], [7.0, 4.288342497, 0.2496251805, 0.0, 0.1, 0.0, 0.0], [8.0, 4.7858445797, 0.2995418888, 0.0, 0.1, 0.0, 0.0], [9.0, 5.0843458292, 0.3294919138, 0.0, 0.1, 0.0, 0.0], [10.0, 6.0793499945, 0.4293253305, 0.0, 0.1, 0.0, 0.0], [11.0, 7.067709136, 0.5787012048, 0.0, 0.2, 0.0, 0.0], [12.0, 8.0477757138, 0.7773705356, 0.0, 0.2, 0.0, 0.0]]

UNDE_SRC = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.19634954084936207 E1=0.0981747704 E2=0.0981747704);
CAVI RFC = (VOLT=1.5e6 HARM=64);
UND CA = (L=0.5 KX=0.1 POLE=4 SLICE=2);
DRIFT CB = (L=0.2);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1 CA CELL1 CELL1 CB 3*CELL1 RFC);
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
UNDE_REF = [[1.173205048517e-07, -2.5511736867e-31, 1.003472792376e-05], [-0.20628736801641, 0.23870613953742, -0.00952897059415], [0.01817924942795, 293112.3213069903, 0.015210615141], [0.0, 0.0, 0.0, 0.0, -3.469446952e-18, -2.7755576018e-20]]


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _nums(s):
    return [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:E[-+]?\d+)?', s)]


@pytest.mark.parametrize('name', ['und1', 'undo'])
def test_und_matches_gs(name):
    ref = globals()[name.upper() + '_REF']
    rows = [[float(v) for v in l.split()] for l in _out(globals()[name.upper() + '_SRC']).splitlines()
            if re.match(r'^ +[0-9]+\.[0-9]+ ', l)]
    assert len(rows) == len(ref)
    for r, q in zip(rows, ref):
        for a, b in zip(r, q):
            assert abs(a - b) <= 5e-9 * max(abs(b), 1e-3) + 1e-12, (name, r, q)


def test_und_emittance_matches_gs():
    got = [_nums(l[4:]) for l in _out(UNDE_SRC).splitlines() if l[:4] in ('EMT ', 'TUN ', 'ALP ', 'COD ')]
    assert len(got) == len(UNDE_REF) == 4
    for k, (r, q) in enumerate(zip(got, UNDE_REF)):
        for a, b in zip(r, q):
            assert abs(a - b) <= 5e-9 * max(abs(b), 1e-3) + 1e-9, (k, r, q)


def test_und_field_keywords_are_undefined_in_k64():
    """K64 は FBX FBY BX BY を Undefined keyword にする (kwFBX が未定義の変数)。pysad も受けない"""
    with pytest.raises(NotImplementedError):
        _out(UND1_SRC.replace('UND U1 = (L=1 ', 'UND U1 = (L=1 FBX=0.5 '))
