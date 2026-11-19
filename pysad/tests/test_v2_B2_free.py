# -*- coding: utf-8 -*-
"""成分ごとの FREE (K64 の tffsfreefix の comp 分岐、`FREE QF.1 K1;`・`FREE B.3 K1;`)。参照は gs (SAD K64、2026-10-19 ビルド、2026-11-18) の出力。

K64 の規則 (gs の実測で確かめた):
  - 既定のキー (QUAD の K1、BEND の ANGLE など) の成分の変数は、同じ名前の要素全体の値 (成分は先頭の成分に結合している)。
    先頭の成分 (QF.1) は要素全体の変数と同じで、結合している成分を Info で知らせる。先頭でない成分 (QF.3) の変数は、結合しているので値を持たず、フィットに効かない。
  - 既定でないキー (BEND の K1 など) の変数は、その成分だけの値。NOKEEPEXP のとき CALC・GO は展開しない (tfvcomp)。
使った条件: 3 GeV、周長 18.4 m の小さな FODO リング (QF・QD 各 4 台、BEND 8 台)。実機の大きなリングでは確かめていない。
"""
import io
import re

import pytest

_NUM = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

_LAT = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ST = ();
BEND B = (L=0.5 ANGLE=0.1);
LINE CELL1 = (QF D2 B D1 QD D2 B D1);
LINE RING = (ST 4*CELL1);
'''


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_LAT + 'FFS;\nPageWidth=1999;\nUSE RING; CELL; CALC;\n' + body + '\nExit;\n', out=buf)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _line(lines, tag):
    for l in lines:
        if l.startswith(tag):
            return l[len(tag):]
    raise AssertionError('no line %r in %r' % (tag, lines[-15:]))


def _nums(s):
    return [float(t) for t in re.findall(_NUM, s)]


def test_component_variable_of_nondefault_key_fits_like_gs():
    # gs: 3 回で収束 (残差 2.1729E-03、4.8421E-08、3.2593E-17)。B.3 だけの K1 で NX = 1.30 にし、NY は 4.915571256742764
    r = _run('FIT $$$ NX 1.30; FREE B.3 K1; GO;\nPrint["N ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]];')
    nx, ny = _nums(_line(r, 'N '))
    assert abs(nx - 8.168140905042474) < 1e-9
    assert abs(ny - 4.915571256742764) < 1e-9
    its = [l for l in r if '(NEWTON)' in l]
    assert [l.split()[2] for l in its] == ['2.1729E-03', '4.8421E-08']
    assert any('Matched. ( 3.2593E-17)' in l for l in r)
    # 要素全体の変数 (FREE B K1) とは別の解 (全部の BEND の K1 を動かすので NY が違う)
    r2 = _run('FIT $$$ NX 1.30; FREE B K1; GO;\nPrint["N ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]];')
    assert abs(_nums(_line(r2, 'N '))[1] - 4.5271083671751) > 0.01 or abs(_nums(_line(r2, 'N '))[1] - 4.915571256742764) > 0.01


def test_component_of_default_key_head_and_coupled_components():
    r = _run('FREE QF.1 K1;\nFREE QF.3 K1;\nFREE QD K1;')
    infos = [l for l in r if 'Info-Component' in l]
    # gs: QF.1 (先頭) は QF.2..4 を結合先として知らせ、QF.3 (先頭でない) は QF.1 に結合していると知らせる
    assert infos == [' ???-FFS-Info-Component QF.2 is coupled to QF.1 .', ' ???-FFS-Info-Component QF.3 is coupled to QF.1 .',
                     ' ???-FFS-Info-Component QF.4 is coupled to QF.1 .', ' ???-FFS-Info-Component QF.3 is coupled to QF.1 .']


def test_non_head_component_variable_has_no_effect():
    # gs: FREE QF.3 K1 だけでは、GO は何も変えない (残差 0.4415 のまま、K1 は 0.7)。pysad はこの変数を作らないので ?No variable. を出す
    r = _run('FIT $$$ NX 1.30; FREE QF.3 K1; GO;\nPrint["N ",Twiss["NX","$$$"]," ",LINE["K1","QF.1"]," ",LINE["K1","QF.3"]];')
    nx, k1a, k1c = _nums(_line(r, 'N '))
    assert abs(nx - 7.503685419504841) < 1e-12 and k1a == 0.7 and k1c == 0.7


def test_mixed_variables_fit_and_table_order_like_gs():
    r = _run('FIT $$$ NX 1.30 NY 0.95;\nFREE QF.1 K1; FREE B.3 K1; FREE B.6 K1; FREE QF.3 K1; FREE B.3 K1; FREE B K1;\nVAR;\nGO;\n'
             'Print["N ",Twiss["NX","$$$"]," ",Twiss["NY","$$$"]," ",LINE["K1","QF.1"]," ",LINE["K1","QF.3"]," ",LINE["K1","B.3"]," ",LINE["K1","B.6"]," ",LINE["K1","B.4"]];')
    # gs: 要素全体の B K1 は、成分 B.6 が使われているのでエラー (最後に FREE した成分の名前)
    assert ' ???-FFS-Error-A component has been already used as variable: B.6' in r
    rows = [l.split()[0] for l in r if l.startswith(' ') and re.match(r' (QF|B)\S*\s+K1\s', l)][:4]
    assert rows == ['QF.3', 'QF.1', 'B.6', 'B.3']          # 同じ要素の変数は、後から FREE したものが先 (gs: QF.3、QF.1、B.6、B.3 の順。P4 で何もしない変数 QF.3 を足した)
    v = _nums(_line(r, 'N '))
    # gs: 8.168140923841754 5.9690399039285 .768018850721068 .768018850721068 -.08865336845062985 -.18088701204518765 0 (反復は 5・6・7 回目)
    assert abs(v[0] - 8.168140923841754) < 1e-8 and abs(v[1] - 5.9690399039285) < 1e-8
    assert abs(v[2] - 0.768018850721068) < 1e-8 and abs(v[3] - v[2]) < 1e-14
    assert abs(v[4] + 0.08865336845062985) < 1e-8 and abs(v[5] + 0.18088701204518765) < 1e-8 and v[6] == 0.0
    assert [l.split()[2] for l in r if '(NEWTON)' in l] == ['6.7769E-02', '1.9660E-03', '9.1155E-06']


def test_fix_component_and_family_like_gs():
    r = _run('FREE B.2 K1; FREE B.5 K1; FREE QF K1; FREE QF.2 ROTATE; FREE QF ROTATE;\nVAR;\nFIX B.2 K1;\nVAR;\nFIX B K1;\nVAR;\nFIX B*;\nVAR;\nFIX *;\nVAR;')
    assert ' ???-FFS-Error-A component has been already used as variable: QF.2' in r
    blocks = []
    cur = None
    for l in r:
        if l.startswith('!Variable'):
            cur = []
            blocks.append(cur)
        elif cur is not None and re.match(r' \S+\s+\S+\s+-?[.\d]', l):
            cur.append(l.split()[0])
    # gs: 同じ要素の変数は後から FREE したものが先。FIX B.2 は B.2 だけ、FIX B (要素全体の語) は成分の変数に触れない、FIX B* は成分も消し、FIX * は全部
    assert blocks[0] == ['QF.2', 'QF', 'B.5', 'B.2']
    assert 'B.2' not in blocks[1] and 'B.5' in blocks[1]
    assert 'B.5' in blocks[2]
    assert 'B.5' not in blocks[3] and 'QF' in blocks[3]
    assert blocks[4] == []


def test_nokeepexp_with_component_variable_of_nondefault_key():
    # gs: NOKEEPEXP のとき、既定でないキーの成分の変数 (B.3 K1) があると、GO は展開しない (NOCALEXP と Info)。QF.1 K1 (既定のキー) では展開する
    r = _run('NOKEEPEXP;\nFIT $$$ NX 1.30;\nFREE QF K1;\nGO;\nFIX *;\nFREE B.3 K1;\nGO;\nFIX *;\nFREE QF.1 K1;\nGO;')
    matched = [l for l in r if l.startswith(' Matched.')][1:]            # 先頭は USE のあとの CALC
    assert matched[0].rstrip().endswith('CALEXP') and not matched[0].rstrip().endswith('NOCALEXP')
    assert matched[1].rstrip().endswith('NOCALEXP')
    assert matched[2].rstrip().endswith('CALEXP') and not matched[2].rstrip().endswith('NOCALEXP')
    assert ' ???-FFS-Info-Element values are not expanded.' in r
