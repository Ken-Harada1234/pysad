"""SOL 領域の Twiss (CALC、INS と CELL)。参照は SAD K64 (gs) の出力。
s1: 単純な領域 (チューン、β)。s2: 領域内の各要素の Twiss (β、α、分散、結合 R2 R3) と回転四極・領域内の SOL。s3: INS・RAD ON。
参照の作り方: local/sol/s*.sad を gs で実行 (docs/27_sol_optics.md)。"""
import io
import re

import pytest

SRC = {}
REF = {}

SRC['s1'] = r'''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.3) QD = (L=0.3 K1=-0.25);
BEND B = (L=1 ANGLE=0.2);
SOL SA = (L=0 BZ=0.3 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR B DR2 QD DR SA DR2 QF DR2 QD DR2 SB DR B DR);
FFS;
USE BL;
CELL; NORAD; NOFLUC;
CALC;
Print["TUN ", Twiss["NX", "$$$"]/(2 Pi), " ", Twiss["NY", "$$$"]/(2 Pi)];
Print["BX  ", Twiss["BX", "$$$"], " ", Twiss["BY", "$$$"]];
Print["EX  ", Twiss["EX", "$$$"], " ", Twiss["EY", "$$$"]];
Exit;
'''
REF['s1'] = [('TUN', '', [0.2927588438819682, 0.10396199217183737]), ('BX', '', [10.241248080000528, 9.006699810727321]), ('EX', '', [1.5123914044930997, 0.018709928673425448])]

SRC['s2'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.7) D4 = (L=0.4) D5 = (L=0.6);
QUAD QF = (L=0.3 K1=0.5) QD = (L=0.3 K1=-0.45) QS = (L=0.2 K1=0.1 ROTATE=0.4) QT = (L=0.25 K1=-0.2);
BEND B = (L=1.2 ANGLE=0.15 E1=0.075 E2=0.075);
SOL SA = (L=0 BZ=1.2 BOUND=1) SM = (L=0 BZ=0.5) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF D1 B D2 QD D1 SA D3 QS D4 SM D5 QT D2 SB D1 B D2);
FFS;
USE BL;
CELL; NORAD; NOFLUC;
CALC;
Print["TUN ", Twiss["NX", "$$$"]/(2 Pi), " ", Twiss["NY", "$$$"]/(2 Pi)];
Print["TW D3 ", Twiss["BX","D3"], " ", Twiss["BY","D3"], " ", Twiss["AX","D3"], " ", Twiss["EX","D3"], " ", Twiss["R2","D3"], " ", Twiss["R3","D3"]];
Print["TW QS ", Twiss["BX","QS"], " ", Twiss["BY","QS"], " ", Twiss["AX","QS"], " ", Twiss["EX","QS"], " ", Twiss["R2","QS"], " ", Twiss["R3","QS"]];
Print["TW D4 ", Twiss["BX","D4"], " ", Twiss["BY","D4"], " ", Twiss["AX","D4"], " ", Twiss["EX","D4"], " ", Twiss["R2","D4"], " ", Twiss["R3","D4"]];
Print["TW SM ", Twiss["BX","SM"], " ", Twiss["BY","SM"], " ", Twiss["AX","SM"], " ", Twiss["EX","SM"], " ", Twiss["R2","SM"], " ", Twiss["R3","SM"]];
Print["TW D5 ", Twiss["BX","D5"], " ", Twiss["BY","D5"], " ", Twiss["AX","D5"], " ", Twiss["EX","D5"], " ", Twiss["R2","D5"], " ", Twiss["R3","D5"]];
Print["TW QT ", Twiss["BX","QT"], " ", Twiss["BY","QT"], " ", Twiss["AX","QT"], " ", Twiss["EX","QT"], " ", Twiss["R2","QT"], " ", Twiss["R3","QT"]];
Print["TW D2 ", Twiss["BX","D2"], " ", Twiss["BY","D2"], " ", Twiss["AX","D2"], " ", Twiss["EX","D2"], " ", Twiss["R2","D2"], " ", Twiss["R3","D2"]];
Print["TW SB ", Twiss["BX","SB"], " ", Twiss["BY","SB"], " ", Twiss["AX","SB"], " ", Twiss["EX","SB"], " ", Twiss["R2","SB"], " ", Twiss["R3","SB"]];
Exit;
'''
REF['s2'] = [('TUN', '', [0.3801951603117019, 0.5624994182984239]), ('TW', 'D3', [1.556052649617704, 220.3843103323368, 0.07628496495573156, 0.4512657284475775, -0.7440131994583328, 0.002934760660367761]), ('TW', 'QS', [1.7575597783197598, 187.3653354162725, -0.40442804359035783, 0.6173599904515145, -1.1700675122839832, 0.0014885367434375604]), ('TW', 'D4', [1.7516116654944058, 165.47935064629834, -0.3409200970146661, 0.6313727539197576, -1.237986732057121, -0.02212113351868325]), ('TW', 'SM', [2.139740491483131, 156.87533035164927, -0.6162963197204796, 0.7011874380755544, -1.4694825140810996, -0.02308130320389755]), ('TW', 'D5', [2.139740491483131, 156.87533035164927, -0.6162963197204796, 0.7011874380755543, -1.4694825140810996, -0.02308130320389755]), ('TW', 'QT', [3.1408793243296675, 141.75567207970948, -1.020520746104661, 0.8057725104333203, -1.7966242547407019, -0.023393475670954778]), ('TW', 'D2', [6.005111109691054, 184.17047319412725, 3.2697997150187756, 0.08865287436454039, 0.1781255156126979, 0.2788418492629534]), ('TW', 'SB', [6.282249203464716, 94.49712253155775, -2.6265782317734065, 1.0490254152078613, -2.1194008299813896, -0.23451190070501948])]

SRC['s3'] = r'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5) D3 = (L=0.7) D4 = (L=0.4) D5 = (L=0.6);
QUAD QF = (L=0.3 K1=0.5) QD = (L=0.3 K1=-0.45) QS = (L=0.2 K1=0.1 ROTATE=0.4) QT = (L=0.25 K1=-0.2);
BEND B = (L=1.2 ANGLE=0.15 E1=0.075 E2=0.075);
SOL SA = (L=0 BZ=1.2 BOUND=1) SM = (L=0 BZ=0.5) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF D1 B D2 QD D1 SA D3 QS D4 SM D5 QT D2 SB D1 B D2);
FFS;
USE BL;
INS; RAD; NOFLUC;
DP=0.002;

CALC;
Print["TUN ", Twiss["NX", "$$$"]/(2 Pi), " ", Twiss["NY", "$$$"]/(2 Pi)];
Print["TW D3 ", Twiss["BX","D3"], " ", Twiss["BY","D3"], " ", Twiss["AX","D3"], " ", Twiss["EX","D3"], " ", Twiss["R2","D3"], " ", Twiss["R3","D3"]];
Print["TW QS ", Twiss["BX","QS"], " ", Twiss["BY","QS"], " ", Twiss["AX","QS"], " ", Twiss["EX","QS"], " ", Twiss["R2","QS"], " ", Twiss["R3","QS"]];
Print["TW D4 ", Twiss["BX","D4"], " ", Twiss["BY","D4"], " ", Twiss["AX","D4"], " ", Twiss["EX","D4"], " ", Twiss["R2","D4"], " ", Twiss["R3","D4"]];
Print["TW SM ", Twiss["BX","SM"], " ", Twiss["BY","SM"], " ", Twiss["AX","SM"], " ", Twiss["EX","SM"], " ", Twiss["R2","SM"], " ", Twiss["R3","SM"]];
Print["TW D5 ", Twiss["BX","D5"], " ", Twiss["BY","D5"], " ", Twiss["AX","D5"], " ", Twiss["EX","D5"], " ", Twiss["R2","D5"], " ", Twiss["R3","D5"]];
Print["TW QT ", Twiss["BX","QT"], " ", Twiss["BY","QT"], " ", Twiss["AX","QT"], " ", Twiss["EX","QT"], " ", Twiss["R2","QT"], " ", Twiss["R3","QT"]];
Print["TW D2 ", Twiss["BX","D2"], " ", Twiss["BY","D2"], " ", Twiss["AX","D2"], " ", Twiss["EX","D2"], " ", Twiss["R2","D2"], " ", Twiss["R3","D2"]];
Print["TW SB ", Twiss["BX","SB"], " ", Twiss["BY","SB"], " ", Twiss["AX","SB"], " ", Twiss["EX","SB"], " ", Twiss["R2","SB"], " ", Twiss["R3","SB"]];
Exit;
'''
REF['s3'] = [('TUN', '', [0.29477636002692437, 0.5074643094495893]), ('TW', 'D3', [31.406098060594395, 11.740953714674696, -13.069815167584059, 0.4572646047808035, -7.3258141645236e-17, 5.507562798683353e-17]), ('TW', 'QS', [52.306241468769244, 8.760108312352344, -16.769998661159242, 0.6216833066064101, -1.734723475976807e-16, 5.603969979034451e-17]), ('TW', 'D4', [58.40381866204118, 8.120948707818721, -13.553389460085693, 0.6638699158700763, 0.0008704124069437279, -0.06482156976309505]), ('TW', 'SM', [69.8117959168413, 7.148913007785871, -14.825750163228465, 0.7389746238692476, 0.016473640883945886, -0.06486827605999296]), ('TW', 'D5', [69.8117959168413, 7.148913007785871, -14.825750163228465, 0.7389746238692475, 0.016473640883945886, -0.06486827605999296]), ('TW', 'QT', [88.8740001460526, 5.867094024365463, -16.773273418904875, 0.8517792843127536, 0.0788850697860075, -0.06491268757823443]), ('TW', 'D2', [5.375859730753708, 11.98542008967992, -2.150043613974747, 0.08983137651166206, 0.0, 0.0]), ('TW', 'SB', [143.74810371541773, 3.4057692131281745, -44.758224556212646, 1.1035056481397376, 0.23678729773363905, -0.09967235133273888])]


def _out(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return buf.getvalue()


def _rows(txt, tags):
    return [(l.split()[0], l.split()[1] if l.split()[0] == 'TW' else '', [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', l.split(None, 2 if l.split()[0] == 'TW' else 1)[-1])])
            for l in txt.splitlines() if l.split()[:1] and l.split()[0] in tags]


@pytest.mark.parametrize('name', ['s1', 's2', 's3'])
def test_sol_optics_matches_gs(name):
    src, ref = SRC[name], REF[name]
    got = _rows(_out(src), ('TUN', 'BX', 'EX', 'TW'))
    assert len(got) == len(ref)
    for (t, n, v), (t2, n2, w) in zip(got, ref):
        assert (t, n) == (t2, n2)
        assert len(v) == len(w)
        for a, b in zip(v, w):
            assert abs(a - b) <= 2e-9 * max(abs(b), 1.0) + 1e-9, (name, t, n, v, w)


def test_sol_calc6d_is_refused():
    """CALC6D・RADCOD つきの CALC は未対応 (Emittance[]・EMIT は solemit。tests/test_sol_emit.py)"""
    src = SRC['s1'].split('CELL;')[0] + 'CELL; RFSW; NORAD; NOFLUC; CALC6D; CALC; Exit;'
    out = _out(src)
    assert '未対応' in out or 'Error' in out or '???' in out
