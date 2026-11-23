"""MAP 要素 (ExternalMap の OPTICS / EMIT / TRACK)。参照は SAD K64 (gs) の出力。
mp1: OPTICS と EMIT (CALC のチューンと軌道)。mp2c: EMIT の {cod1, trans} (Emittance[])、mp2d: {cod1, trans, dtrans, dbeam}。
mp3: TRACK (座標の書き換え、フラグ 0 で粒子を落とす、落ちた粒子の復活)。
参照の作り方: local/wake/mp*.sad を gs で実行 (docs/19_other_elements.md)。"""
import io
import re

import pytest

SRC = {}
REF = {}

SRC['mp1'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MAP P1 = (L=0.5);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 P1);
FFS;
$FORM="18.14";
USE RING;
CELL;
ExternalMap["OPTICS",LINE["POSITION","P1"],cod_]:=(
  {cod+{0,0.0003,0,0,0,0}, {{1,0.1,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}}});
ExternalMap["EMIT",LINE["POSITION","P1"],cod_]:=(
  {cod+{0,0.0003,0,0,0,0}, {{1,0.1,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0},{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}}});
CALC;
Print["TUN ", Twiss["NX","$$$"]/(2 Pi), " ", Twiss["NY","$$$"]/(2 Pi)];
Print["COD ", Twiss["DX","$$$"], " ", Twiss["DPX","$$$"]];
Print["BX ", Twiss["BX","$$$"], " ", Twiss["BX", 3]];
Exit;
'''
REF['mp1'] = [('TUN', [2.7156833010379, 1.86689215267466]), ('COD', [-0.000570723210672, -9.56051932975e-05])]

SRC['mp2c'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MAP P1 = (L=0.5);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1 P1 6*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
mm = {{0.99875026039496628,0.099875026039496634,0.049979169270678331,0.0049979169270678338,0,0},{0,0.99875026039496628,0,0.049979169270678331,0,0},{-0.049979169270678331,0.0024989584635339169,0.99875026039496628,-0.049937513019748317,0,0},{0,-0.049979169270678331,0,0.99875026039496628,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};
ExternalMap["EMIT",LINE["POSITION","P1"],cod_]:=({mm.cod+{0,0.0003,0,0,0,0},mm});
RFSW;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["DMP ", DampingRate/.em];
Print["COD ", OrbitAtExit/.em];
Print["SPR ", {MomentumSpread/.em, BunchLength/.em}];
Exit;
'''
REF['mp2c'] = [('EMT', [2.635809660299e-07, 1.750998030298e-08, 2.795843174682e-05]), ('TUN', [-0.28534830304845, -0.13540600394458, -0.01186657763182]), ('ALP', [0.03239140715536, 781409.766002063, 0.06213656285128]), ('DMP', [-0.000124703406952, -0.000130484101894, -0.000265668702221]), ('COD', [0.0009193706906851, 0.0004101715223099, -0.000128747032329, 8.217532302146e-05, -1.4818620071e-13, -0.000152743209303]), ('SPR', [0.00118843466889, 0.02353836099367])]

SRC['mp2d'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MAP P1 = (L=0.5);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1 P1 6*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
mm = {{0.99875026039496628,0.099875026039496634,0.049979169270678331,0.0049979169270678338,0,0},{0,0.99875026039496628,0,0.049979169270678331,0,0},{-0.049979169270678331,0.0024989584635339169,0.99875026039496628,-0.049937513019748317,0,0},{0,-0.049979169270678331,0,0.99875026039496628,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}};
ExternalMap["EMIT",LINE["POSITION","P1"],cod_]:=({mm.cod+{0,0.0003,0,0,0,0},mm,{{0,0,0,0,0,0},{0,-2e-5,0,0,0,0},{0,0,0,0,0,0},{0,0,0,-3e-5,0,0},{0,0,0,0,0,0},{0,0,0,0,0,-1e-5}},{{0,0,0,0,0,0},{0,1e-14,0,2e-15,0,3e-14},{0,0,0,0,0,0},{0,0,0,2e-14,0,1e-15},{0,0,0,0,0,0},{0,0,0,0,0,5e-13}}});
RFSW;
em = Emittance[];
Print["EMT ", Emittances/.em];
Print["TUN ", Tunes/.em];
Print["ALP ", {MomentumCompaction/.em, EnergyLossU0/.em, EquilibriumPosition/.em}];
Print["DMP ", DampingRate/.em];
Print["COD ", OrbitAtExit/.em];
Print["SPR ", {MomentumSpread/.em, BunchLength/.em}];
Exit;
'''
REF['mp2d'] = [('EMT', [2.441501499527e-07, 1.581422082962e-08, 2.745132807098e-05]), ('TUN', [-0.28534830304845, -0.13540600394458, -0.01186657763182]), ('ALP', [0.03239140715536, 781409.766002063, 0.06213656285128]), ('DMP', [-0.000134751577342, -0.000145428760045, -0.00027066483761]), ('COD', [0.0009193706906851, 0.0004101715223099, -0.000128747032329, 8.217532302146e-05, -1.4818620071e-13, -0.000152743209303]), ('SPR', [0.00117760836642, 0.02332386900523])]

SRC['mp3'] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.6) QD = (L=0.3 K1=-0.5);
MAP P1 = (L=0.2);
MARK ST = ();
LINE BL = (ST QF DR P1 DR2 QD DR P1 DR);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
n = 20;
x0  = Table[1e-3*Sin[0.37*i], {i, n}];
px0 = Table[1e-5*Cos[0.11*i], {i, n}];
y0  = Table[1e-3*Cos[0.23*i], {i, n}];
py0 = Table[1e-5*Sin[0.71*i], {i, n}];
z0  = Table[0.01*Sin[0.13*i]+1e-4*i-0.01, {i, n}];
dp0 = Table[1e-3*Cos[0.31*i], {i, n}];
fl0 = Table[1, {i, n}];
ExternalMap["TRACK", LINE["POSITION","P1"], nt_, x_] := (Print["MAPA ", nt, " ", Length[x[[1]]], " ", Plus@@x[[1]], " ", Plus@@x[[7]]];
  {x[[1]]*0.9+1e-4, x[[2]]+2e-5, x[[3]], x[[4]]*1.1, x[[5]], x[[6]], x[[7]]});
ExternalMap["TRACK", LINE["POSITION","P1.2"], nt_, x_] := (Print["MAPB ", nt, " ", Length[x[[1]]], " ", Plus@@x[[1]], " ", Plus@@x[[7]]];
  {x[[1]], x[[2]], x[[3]], x[[4]], x[[5]], x[[6]], If[Abs[#]<1.2e-3,1,0]&/@x[[1]]});
r = TrackParticles[{1, {x0, px0, y0, py0, z0, dp0, fl0}}, 10, 1, 4];
$FORM = "16.10";
res = r[[2]];
Do[Print[k, " ", Plus@@res[[k]], " ", res[[k, 1]], " ", res[[k, 7]], " ", res[[k, 20]]], {k, 7}];
Exit;
'''
REF['mp3'] = [('MAPA', [1.0, 20.0, 0.0007244572433262667, 20.0]), ('MAPB', [1.0, 20.0, 0.002750226764543056, 20.0]), ('MAPA', [2.0, 20.0, 0.0016109970242288055, 20.0]), ('MAPB', [2.0, 20.0, 0.0031896297262501176, 20.0]), ('MAPA', [3.0, 20.0, -0.0016161574953568244, 16.0]), ('MAPB', [3.0, 20.0, 0.0032850078720438203, 16.0]), ('MAPA', [4.0, 20.0, 0.00405920137659319, 20.0]), ('MAPB', [4.0, 20.0, 0.007081826474673882, 20.0]), ('1.0000000000', [0.0088177573, 0.00047840973948, 0.00073524242822, 0.0013525168]), ('2.0000000000', [0.0017371436, 9.1177348133e-05, 0.00015043571838, 0.00029335671473]), ('3.0000000000', [0.0013428423, -0.0011457198, 4.5490737252e-05, 0.00013333826213]), ('4.0000000000', [6.2701984991e-06, -2.682971047e-05, 1.2515344479e-05, -1.021616653e-05]), ('5.0000000000', [-0.0383242336, -0.0086051241, -0.0014052117, -0.0028459803]), ('6.0000000000', [0.0025284101, 0.00095233356989, -0.0005639850576, 0.00099654209702]), ('7.0000000000', [20.0, 1.0, 1.0, 1.0])]


_NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
_KEYS = ('TUN', 'COD', 'EMT', 'DMP', 'SPR', 'ALP', 'MAPA', 'MAPB')


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _rows(txt):
    out = []
    for l in txt.splitlines():
        t = l.split()
        if not t:
            continue
        if t[0] in _KEYS or re.match(r'^[1-7]\.0000000000$', t[0]):
            out.append((t[0], [float(v) for v in re.findall(_NUM, l.split(None, 1)[1] if len(t) > 1 else '')]))
    return out


@pytest.mark.parametrize('name', ['mp1', 'mp2c', 'mp2d', 'mp3'])
def test_map_matches_gs(name):
    got = _rows(_run(SRC[name]))
    ref = REF[name]
    assert [g[0] for g in got] == [r[0] for r in ref]
    for (k, v), (k2, w) in zip(got, ref):
        assert len(v) == len(w), (name, k)
        for a, b in zip(v, w):
            assert abs(a - b) <= 1e-8 * abs(b) + 2e-12, (name, k, v, w)


def test_map_default_is_identity_with_zero_length():
    """ExternalMap が定義されていなければ恒等で長さも 0 (K64: tgmap の 9000、MAP の L は使われない)"""
    base = SRC['mp3'].split('ExternalMap')[0]
    a = _run(base + 'CELL; CALC; Print["S ", Twiss["S", "$$$"]]; Exit;')
    assert abs(float(re.search(r'S\s+(' + _NUM + ')', a).group(1)) - 4.1) < 1e-12


def test_lossmap_records_map_element_number():
    """v2.4.0: MAP (ExternalMap["TRACK"]) が落とした粒子の LOSSMAP の要素番号 (temap.f: kptbl(j,4) = l、kptbl(j,5) = nt)。
    gs: 粒子 1・18 が 4 周目に MAP (要素 8) で落ちる。以前は要素番号が 0 だった"""
    from pysad.ffscmd import run_sad_string
    s = SRC['mp3'].replace('TRPT; INS; NORAD; NOFLUC;', 'TRPT; INS; NORAD; NOFLUC; LOSSMAP;').replace('x0  = Table[1e-3*Sin', 'x0  = Table[3e-3*Sin')
    s = s[:s.index('$FORM = "16.10";')] + 'Print["F ", r[[2,7]]]; Print["L ", r[[2,8]]]; Print["T ", r[[2,9]]];\nExit;\n'
    b = io.StringIO()
    run_sad_string(s, out=b)
    got = {l[0]: l for l in b.getvalue().splitlines() if l[:2] in ('F ', 'L ', 'T ')}
    assert got['F'].replace(' ', '') == 'F{0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1}'
    assert got['L'].replace(' ', '') == 'L{8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0}'
    assert got['T'].replace(' ', '') == 'T{4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,4,5,5}'
