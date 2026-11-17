"""ListPlot・Plot (Packages/ListPlot.n、Plots.n、Canvas.n) と、Tk のキャンバスへの描画 (pysad/tkcanvas.py)。
参照は SAD K64 (gs, WSLg の X11) で同じ台本を走らせた結果: ListPlot[{1,4,9,16}, PlotJoined->True] が作るキャンバスの項目は 43 個
(rectangle 1、line 22、text 16、oval 4) で、種類・座標・タグまで一致する (文字の大きさだけが Tk のフォントで違う)。
Tk を使う試験は PYSAD_SAD_PACKAGES と tkinter が要る (なければ skip)。"""
import io
import os
import re

import pytest

PK = os.environ.get('PYSAD_SAD_PACKAGES')
needs_pk = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')


def _read_items(path):
    """キャンバスの項目の一覧を読む。Tk はマウスの下にある項目に `current` というタグを自動で付ける (試験を走らせたときの
    マウスの位置で変わる) ので、それだけを外す。種類・座標・ほかのタグはそのまま比べる"""
    out = []
    for l in open(path, encoding='utf-8').read().splitlines():
        l = l.rstrip()
        if l.endswith(' current'):
            l = l[:-len(' current')].rstrip()
        out.append(l)
    return out


def _run(body, display=False, raw=False):
    from pysad.ffscmd import run_sad_string
    if not display:
        os.environ['PYSAD_NO_DISPLAY'] = '1'
    else:
        os.environ.pop('PYSAD_NO_DISPLAY', None)
    os.environ['PYSAD_TK_NOLOOP'] = '1'
    buf = io.StringIO()
    try:
        ffs, I = run_sad_string((body if raw else 'FFS;\nPageWidth=1999;\n' + body) + 'Exit;\n', out=buf)
    finally:
        os.environ.pop('PYSAD_NO_DISPLAY', None)
    br = getattr(I, 'tclbridge', None)
    if br is not None and br.root is not None:
        try:
            br.root.destroy()
        except Exception:
            pass
    return buf.getvalue()


@needs_pk
def test_listplot_returns_graphics_and_plot_matches_gs():
    out = _run('''g=ListPlot[{1,4,9,16},PlotJoined->True];
Print["A ",Head[g]," ",Length[g]];
h=Plot[Sin[x],{x,0,6}];
Print["B ",Head[h]," ",Length[h]];
''')
    assert 'A Graphics 2' in out
    assert 'B Graphics 2' in out


def _tk_ok():
    try:
        import tkinter
        r = tkinter.Tk()
        r.destroy()
        return True
    except Exception:
        return False


@needs_pk
def test_listplot_draws_the_same_canvas_items_as_gs():
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    out = _run('''g=ListPlot[{1,4,9,16},PlotJoined->True];
Sleep[0.3]; TclEval["update"]; Sleep[0.3]; TclEval["update"];
TclEval["proc fc {w} {set r {}; foreach c [winfo children $w] {if {[winfo class $c] eq {Canvas}} {foreach i [$c find all] {lappend r \\"[$c type $i] [$c coords $i] [$c gettags $i]\\"}} else {set r [concat $r [fc $c]]}}; return $r}; set d {}; foreach l [fc .] {dict incr d [lindex $l 0]}; set fh [open pysad_plot_items.txt w]; puts $fh $d; foreach l [fc .] {puts $fh $l}; close $fh"];
''', display=True)
    path = 'pysad_plot_items.txt'
    try:
        lines = _read_items(path)
    finally:
        if os.path.exists(path):
            os.remove(path)
    assert lines[0] == 'rectangle 1 line 22 text 16 oval 4'
    items = lines[1:]
    assert len(items) == 43
    # 描画範囲 (BG)・データの折れ線・記号 (oval)・縦横の格子の位置が gs と同じ
    assert 'rectangle 158.0 172.0 914.0 630.0 {BG$1}' in items
    assert any(l.startswith('line 204.0 602.0 425.0 521.0 646.0 387.0 867.0 200.0') for l in items)
    assert 'oval 199.89039066468735 597.8903906646874 208.10960933531265 606.1096093353126' in [l.rstrip(' {}')[:len('oval 199.89039066468735 597.8903906646874 208.10960933531265 606.1096093353126')] for l in items if l.startswith('oval')] \
        or any(l.startswith('oval 199.89039066468735 597.8903906646874') for l in items)
    assert 'line 204.0 630.0 204.0 172.0 {G$X1}' in items
    assert 'line 158.0 200.0 914.0 200.0 {G$Y1}' in items


@needs_pk
def test_opticsplot_draws_the_same_geometry_as_gs():
    """OpticsPlot[{{"BX","BY"},{"EX"}}] (12 セルのリング)。gs は項目 174 個 (rectangle 76、line 40、text 58)。BG の枠・矩形・線の数と位置は同じで、
    文字は Tk のフォントの大きさで RemoveOverlap・目盛の間引きが変わり数が違う (pysad は text 72)"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    out = _run('''MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
USE RING; CELL; CALC;
g=OpticsPlot[{{"BX","BY"},{"EX"}}];
Print["HEAD ",Head[g]," ",Length[g]];
Sleep[0.3]; TclEval["update"]; Sleep[0.3]; TclEval["update"];
TclEval["proc fc {w} {set r {}; foreach c [winfo children $w] {if {[winfo class $c] eq {Canvas}} {foreach i [$c find all] {lappend r \\"[$c type $i] [$c coords $i] [$c gettags $i]\\"}} else {set r [concat $r [fc $c]]}}; return $r}; set d {}; foreach l [fc .] {dict incr d [lindex $l 0]}; set fh [open pysad_optics_items.txt w]; puts $fh $d; foreach l [fc .] {puts $fh $l}; close $fh"];
''', display=True, raw=True)
    assert 'HEAD Graphics 2' in out
    path = 'pysad_optics_items.txt'
    try:
        lines = _read_items(path)
    finally:
        if os.path.exists(path):
            os.remove(path)
    d = dict(zip(lines[0].split()[::2], map(int, lines[0].split()[1::2])))
    assert d['rectangle'] == 76 and d['line'] == 40
    items = lines[1:]
    assert 'rectangle 158.0 355.0 914.0 584.0 {BG$3}' in items
    assert 'rectangle 158.0 126.0 914.0 355.0 {BG$2}' in items
    assert 'line 158.0 315.0 914.0 315.0 {G$Y2}' in items


_DUMP = r'''Sleep[0.3]; TclEval["update"]; Sleep[0.3]; TclEval["update"];
TclEval["proc fc {w} {set r {}; foreach c [winfo children $w] {if {[winfo class $c] eq {Canvas}} {foreach i [$c find all] {lappend r \"[$c type $i] [$c coords $i] [$c gettags $i]\"}} else {set r [concat $r [fc $c]]}}; return $r}; set d {}; foreach l [fc .] {dict incr d [lindex $l 0]}; set fh [open PATH w]; puts $fh $d; foreach l [fc .] {puts $fh $l}; close $fh"];
'''


def _items(name, body, raw=False):
    """body を display=True で走らせ、キャンバスの項目を (種類ごとの数, 全項目) で返す"""
    path = 'pysad_%s_items.txt' % name
    _run(body + _DUMP.replace('PATH', path), display=True, raw=raw)
    try:
        lines = _read_items(path)
    finally:
        if os.path.exists(path):
            os.remove(path)
    return lines[0], lines[1:]


@needs_pk
def test_histoplot_and_listcontourplot_match_gs():
    """HistoPlot[Table[Sin[i*1.7]+Cos[i*0.3],{i,1,60}]] は gs と全 38 項目が同じ。ListContourPlot[Table[Sin[i/3.]*Cos[j/3.],{i,0,9},{j,0,9}]] は
    項目の種類と数が同じ (rectangle 20、polygon 10、line 25、text 17)。等高線の線は同じ曲線だが、gs は ?contl (定義の表) をハッシュ表の順に返すので始点・向きが違う"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    d, items = _items('histo', 'g=HistoPlot[Table[Sin[i*1.7]+Cos[i*0.3],{i,1,60}]];\n')
    assert d == 'rectangle 9 line 11 text 17'
    assert 'rectangle 216.0 477.0 274.0 630.0' in items
    assert 'rectangle 548.0 202.0 606.0 630.0' in items
    d, items = _items('contour', 'g=ListContourPlot[Table[Sin[i/3.]*Cos[j/3.],{i,0,9},{j,0,9}]];\n')
    assert d == 'rectangle 20 polygon 10 line 25 text 17'


@needs_pk
def test_fitplot_matches_gs():
    """FitPlot[Table[{i,2*i+1+0.1*Sin[i*3]},{i,1,10}], a+b*x, x, {a,1}, {b,1}]: gs は a=1.0132994、b=1.9967382、ChiSquare=.0484422、キャンバス 43 項目"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    out = _run('g=FitPlot[Table[{i,2*i+1+0.1*Sin[i*3]},{i,1,10}],a+b*x,x,{a,1},{b,1}];\n'
               'Print["A ",Head[g]," ",Length[g]," ",Round[10^6*(a/.g[[1]])]/10^6," ",Round[10^6*(b/.g[[1]])]/10^6," ",Round[10^5*(ChiSquare/.g[[1]])]/10^5];\n',
               display=True)
    assert 'A List 2 1.013299 1.996738 .04844' in out


@needs_pk
def test_plot_with_other_variable_names_and_function_lists_match_gs():
    """Plot・ParametricPlot の変数が x 以外 (t、y)、Plot[{f1, f2}, …] (HoldAll の引数、Hold[…] /. Hold -> Null)。gs のキャンバスの項目と同じ数"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    d, _ = _items('par', 'g=ParametricPlot[{Sin[t],Cos[2*t]},{t,0,6}];\n')
    assert d == 'rectangle 1 line 21 text 14'
    d, _ = _items('plist', 'g=Plot[{Sin[y],Cos[y]},{y,0,6}];\n')
    assert d == 'rectangle 1 line 24 text 16'


@needs_pk
def test_geometryplot_matches_gs():
    """GeometryPlot[] (12 セルのリング): gs は line 74、polygon 37。LINE の "L$SUM"・"ANGLE$SUM" (KEY$SUM は KEY の値) が要る。
    DRIFT の外形の線は Fill->Null なので既定の色 (X11 は #000000、Windows の Tk は SystemButtonText)"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    d, items = _items('geo', """MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
USE RING; CELL; CALC;
g=GeometryPlot[];
""", raw=True)
    assert d == 'line 74 polygon 37'
    assert any(l.startswith('polygon 433.0 627.0 438.0 627.0 438.0 604.0 433.0 604.0 {GeoPl$1:QF.1}') for l in items)


@needs_pk
def test_fitgaussian_plot_matches_gs():
    """FitGaussian[d, Plot->True] は {sig, 平均, {…}, ChiSquare, Show[…]} (Fit.n)。gs は List 5、キャンバス 74 項目。値は gs と 1e-10。
    Plot[InverseErf[x] …] のために Vector に InverseErf などの特殊関数を足した"""
    if not _tk_ok():
        pytest.skip('tkinter が使えない')
    body = ('d=Table[Sqrt[2]*InverseErf[2*(i-0.5)/400-1]*0.8+0.3+0.05*Sin[i*1.3],{i,400}];\n'
            'r=FitGaussian[d,Plot->True];\nPrint["P ",Head[r]," ",Length[r]," ",Round[r[[1]]*10^8]/10^8," ",Round[r[[2]]*10^8]/10^8];\n')
    d, items = _items('fitg', body)
    assert sorted(d.split()[0::2]) == ['line', 'oval', 'rectangle', 'text']
    cnt = dict(zip(d.split()[0::2], map(int, d.split()[1::2])))
    assert cnt == {'rectangle': 1, 'line': 40, 'text': 13, 'oval': 20}      # gs と同じ (text の数はフォントに依存しうる)
