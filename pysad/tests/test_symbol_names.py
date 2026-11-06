"""裸のシンボルを要素名として使う (`qlist = {Q1A, Q2A}; LINE["K1", qlist[i]]`: tfgetstrns はシンボルの名前を使う) と、
要素の型名と同じ名前の LINE (`TEST = (-N2 2*NCELL)`) の読み込み (docs/history/work_trpt_acc.md 7 節、pf2 disp_q_K1_r)。"""
import io

from pysad.ffscmd import run_sad_string

SRC = '''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5);
QUAD QF = (L=0.3 K1=1.2) QD = (L=0.3 K1=-1.1);
SEXT SF = (L=0.1 K2=3) SD = (L=0.1 K2=-4);
MARK MK = ();
LINE CELL = (MK QF D SF D QD D SD D)
     N2 = (CELL CELL)
     TEST = (-N2 2*CELL)
     RING = (8*CELL);
FFS;
USE RING;
CELL;
CALC;
qlist = {QF, QD};
slist = {SF, SD};
k1 = Table[LINE["K1", qlist[i]], {i, 1, Length[qlist]}];
k2 = Table[LINE["K2", slist[i]] * 2, {i, 1, Length[slist]}];
bx = Twiss["BX", QF];
bxs = Twiss["BX", "QF"];
el = Element["K1", QD];
'''


def test_bare_symbol_as_element_name():
    out = io.StringIO()
    ffs, I = run_sad_string(SRC, out=out)
    assert [float(x) for x in I.get('k1')] == [1.2, -1.1]
    assert [float(x) for x in I.get('k2')] == [6.0, -8.0]
    assert I.get('bx') == I.get('bxs') and float(I.get('bx')) > 0     # シンボルは文字列と同じ扱い
    assert float(I.get('el')) == -1.1


def test_line_named_like_element_type():
    out = io.StringIO()
    ffs, I = run_sad_string(SRC + 'USE TEST;\nINS;\nCALC;\n', out=out)
    assert ffs.name == 'TEST' and ffs.nlat == 4 * 9 + 1
    assert 'Warning' not in out.getvalue()
