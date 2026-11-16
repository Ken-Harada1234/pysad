"""DELCOR (tcorr.f / tvcorr.f): リングに沿った相関のあるランダムなずれ。参照は SAD K64 (gs) の出力 (乱数も gs と同じ列)。
dc1: CELL。相関長 2 m、相関長 0 (ランダムウォーク)、ADD とワイルドカード。dc2: INS、MARK・SEXT・BEND を含み、`*` で全要素。
MARK には DX DY ではなく DETR と DPX が入る (K64 の latt+15, +17)。参照の作り方: local/wake/dc*.sad を gs で実行 (docs/29_delcor.md)。"""
import io
import re

import pytest

SRC = {}
REF = {}

SRC['dc1'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
$FORM="18.14";
USE RING;
CELL;
DELCOR 2 1e-4 2e-4 QF QD;
CALC;
Print["TUN ", Twiss["NX","$$$"]/(2 Pi), " ", Twiss["NY","$$$"]/(2 Pi)];
Print["COD ", Twiss["DX","$$$"], " ", Twiss["DY","$$$"], " ", Twiss["DX",5], " ", Twiss["DY",12]];
Print["EL ", LINE["DX","QF.3"], " ", LINE["DY","QF.3"], " ", LINE["DX","QD.7"], " ", LINE["DY","QD.7"]];
DELCOR 0 1e-4 2e-4 B;
CALC;
Print["TU2 ", Twiss["NX","$$$"]/(2 Pi), " ", Twiss["NY","$$$"]/(2 Pi)];
Print["CO2 ", Twiss["DX","$$$"], " ", Twiss["DY","$$$"]];
Print["EL2 ", LINE["DX","B.3"], " ", LINE["DY","B.5"]];
DELCOR ADD 1 1e-4 1e-4 Q*;
CALC;
Print["TU3 ", Twiss["NX","$$$"]/(2 Pi), " ", Twiss["NY","$$$"]/(2 Pi)];
Print["EL3 ", LINE["DX","QF.3"], " ", LINE["DX","QD.7"]];
Exit;
'''
REF['dc1'] = [('TUN', [2.70675708579474, 1.86687143519596]), ('COD', [1.29782659648e-05, 1.09420667009e-05, 3.598559038786e-06, -7.77406304187e-06]), ('TU2', [2.70675360335206, 1.86686813817295]), ('CO2', [4.24423510175e-06, 6.224452081767e-05]), ('TU3', [2.70677007546656, 1.86688680624028])]

SRC['dc2'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.2 K2=3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 SF D2 B D2 QD D1);
LINE RING = (ST 4*CELL1);
FFS;
$FORM="18.14";
USE RING;
INS;
DELCOR 3 1e-3 2e-3 ST QF SF;
CALC;
Print["TUN ", Twiss["NX","$$$"]/(2 Pi), " ", Twiss["NY","$$$"]/(2 Pi)];
Print["COD ", Twiss["DX","$$$"], " ", Twiss["DY","$$$"], " ", Twiss["DX",5], " ", Twiss["DY",12]];
Print["CO2 ", Twiss["DX","ST"], " ", Twiss["DY","ST"]];
DELCOR 1 1e-3 1e-3 B* *;
CALC;
Print["CO3 ", Twiss["DX","$$$"], " ", Twiss["DY","$$$"], " ", Twiss["DX",9], " ", Twiss["DY",14]];
Exit;
'''
REF['dc2'] = [('TUN', [1.18776392271283, 0.6588668538414]), ('COD', [5.981605880572e-05, 0.0004755258174011, 5.10240120526e-05, -0.000479408701361]), ('CO2', [0.0, 0.0]), ('CO3', [0.00190917567511, 0.00123610873858, 0.00495113311102, -0.000695276118903])]


_NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
_KEYS = ('TUN', 'COD', 'TU2', 'CO2', 'CO3', 'TU3')


def _rows(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return [(l.split()[0], [float(v) for v in re.findall(_NUM, l.split(None, 1)[1])])
            for l in buf.getvalue().splitlines() if l.split() and l.split()[0] in _KEYS]


@pytest.mark.parametrize('name', ['dc1', 'dc2'])
def test_delcor_matches_gs(name):
    got = _rows(SRC[name])
    ref = REF[name]
    assert [g[0] for g in got] == [r[0] for r in ref]
    for (k, v), (k2, w) in zip(got, ref):
        for a, b in zip(v, w):
            assert abs(a - b) <= 1e-9 * abs(b) + 1e-12, (name, k, v, w)


def test_delcor_cavity_writes_harm_and_phi_like_gs():
    """CAVI への DELCOR は K64 の添字のずれで DX・DY ではなく HARM・PHI に書く (tcorr.f の le=latt+3)。1 回目の値は gs と一致 (docs/29、B62)"""
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(r"""MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64 DX=1e-5 DY=2e-5);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1 RFC 6*CELL1);
FFS;
PageWidth=1999;
USE RING;
DELCOR 3 1e-4 1e-4 RFC;
Print["B ",Element[{"DX","DY","VOLT","HARM","PHI"},"RFC"]];
Exit;
""", out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('B ')][0]
    nums = [float(x) for x in re.findall(r'[-+]?\d*\.?\d+(?:e[-+]?\d+)?', line[2:])]
    want = [1e-5, 2e-5, 1500000.0, -.00015145432812411818, .00011480572629164664]
    for g, w in zip(nums, want):
        assert abs(g - w) < 1e-9 * abs(w)
