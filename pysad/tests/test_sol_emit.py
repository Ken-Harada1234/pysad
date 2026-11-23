"""SOL 領域の Emittance[] (放射は bz 付きの tradke、写像は追跡のカーネルの数値微分)。参照は SAD K64 (gs) の出力。
e1: 12 セルのリングの 1 セルに SOL 領域 (SA GEO=1 BZ=2、回転四極、領域内の SOL SM)。QFE の DX DY で軌道が領域を斜めに通る。
e2: 強い場 (BZ=5)・長い四極 (K1=0.3, L=0.5)・大きい軌道 (放射の分割数が増える)。e1z: BZ=0 の領域 (普通の要素になる)。
参照の作り方: local/sol/e*.sad を gs で実行 (docs/17_sol.md)。SOL には GEO=1 が要る (K64 の tsgeo。なければ S 座標が壊れ周長が狂う)。"""
import io
import re

import pytest

SRC = {}
REF = {}

SRC['e1'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.25);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QFE = (L=0.3 K1=0.7 DX=3e-4 DY=-2e-4) QS = (L=0.2 K1=0.05 ROTATE=0.3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
SOL SA = (L=0 BZ=2 BOUND=1 GEO=1) SM = (L=0 BZ=0.8) SB = (L=0 BZ=0 BOUND=1);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELLS = (QFE D2 B D2 QD D3 SA D3 QS D3 SB D3);
LINE RING = (ST 11*CELL1 CELLS RFC);
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
REF['e1'] = {'EMT': [2.706114500404e-07, 2.245979871823e-08, 2.778701689678e-05], 'TUN': [-0.26504588346481, -0.12445663744341, -0.01173531080122], 'ALP': [0.03168015702303, 781502.932093495, 0.06241737927533], 'COD': [0.0004636825843481, 0.0002959927995334, -5.88837617083e-05, 0.0001895161502015, -3.1606229567e-14, -9.15330628223e-05]}

SRC['e2'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.2);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QFE = (L=0.3 K1=0.7 DX=2e-3 DY=-1.5e-3) QS = (L=0.5 K1=0.3 ROTATE=0.3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
SOL SA = (L=0 BZ=5 BOUND=1 GEO=1) SM = (L=0 BZ=2.5) SB = (L=0 BZ=0 BOUND=1);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELLS = (QFE D2 B D2 QD D3 SA D3 QS D3 SM D3 QD D3 SB D3);
LINE RING = (ST 11*CELL1 CELLS RFC);
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
REF['e2'] = {'EMT': [4.811746951828e-07, 7.984466593388e-08, 2.986125136093e-05], 'TUN': [-0.35517144849815, 0.0, -0.01200735787346], 'ALP': [0.0331586643553, 781038.9608041726, 0.06346564914798], 'COD': [0.00270094066643, 0.00219949118538, 0.00116781476859, 0.0006072307014011, -5.8156185184e-11, -0.000954677884633]}

SRC['e1z'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.25);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6) QFE = (L=0.3 K1=0.7 DX=3e-4 DY=-2e-4) QS = (L=0.2 K1=0.05 ROTATE=0.3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
SOL SA = (L=0 BZ=0 BOUND=1 GEO=1) SM = (L=0 BZ=0) SB = (L=0 BZ=0 BOUND=1);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELLS = (QFE D2 B D2 QD D3 SA D3 QS D3 SB D3);
LINE RING = (ST 11*CELL1 CELLS RFC);
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
REF['e1z'] = {'EMT': [2.726962061573e-07, 7.02684400415e-09, 2.775676479358e-05], 'TUN': [-0.26559059046515, -0.12826554367745, -0.01175894403078], 'ALP': [0.03180815083801, 781520.5234392975, 0.06241894416425], 'COD': [0.0005119456473012, 0.0003065937446587, -0.000148573548866, 0.0002735228152398, -2.3185592149e-14, -8.24297775399e-05]}


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


_NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'


def _vals(txt):
    out = {}
    for l in txt.splitlines():
        t = l.split()
        if t and t[0] in ('EMT', 'TUN', 'ALP', 'COD'):
            out[t[0]] = [float(v) for v in re.findall(_NUM, l.split(None, 1)[1])]
    return out


# (相対許容, 絶対許容): ε_x ε_y は 1e-7、ε_z は 1e-4 (z の行は数値微分)、チューンは絶対 1e-6、軌道は 1e-6
TOL = {'EMT': [(1e-7, 0), (1e-6, 0), (1e-4, 0)], 'TUN': [(0, 1e-6)] * 3, 'ALP': [(1e-4, 0), (1e-6, 0), (1e-5, 0)],
       'COD': [(1e-6, 1e-9)] * 6}


@pytest.mark.parametrize('name', ['e1', 'e2', 'e1z'])
def test_sol_emittance_matches_gs(name):
    got = _vals(_out(SRC[name]))
    ref = REF[name]
    for key in ('EMT', 'TUN', 'ALP', 'COD'):
        assert len(got[key]) == len(ref[key]), (name, key)
        for a, b, (rt, at) in zip(got[key], ref[key], TOL[key]):
            assert abs(a - b) <= rt * abs(b) + at + 1e-300, (name, key, got[key], ref[key])
