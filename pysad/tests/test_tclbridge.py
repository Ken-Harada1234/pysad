"""Tcl/Tk の接点 (pysad/tclbridge.py) と、SAD の Tkinter.n / KBMainFrame.n を pysad のインタプリタで動かした結果の試験。
SAD のパッケージの場所を環境変数 PYSAD_SAD_PACKAGES で指定し、tkinter が使える画面があるときだけ走る (無ければ skip)。
ウィジェットの構成は SAD (WSL の gs、Tk 8.6 X11) で同じスクリプトを実行して取った構成と比べる。"""
import io
import os

import pytest

PK = os.environ.get('PYSAD_SAD_PACKAGES')
pytestmark = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')


@pytest.fixture()
def interp():
    os.environ['PYSAD_TK_NOLOOP'] = '1'
    import time
    err = None
    for _ in range(5):                                       # Windows の Tcl は、まれに起動直後の初期化に失敗するのでやり直す
        try:
            import tkinter
            r = tkinter.Tk()
            r.destroy()
            err = None
            break
        except Exception as e:                               # 画面が無い
            err = e
            time.sleep(0.3)
    if err is not None:
        pytest.skip('tkinter が使えない: %s' % err)
    from pysad.sadscript import Interp
    I = Interp()
    I.out = io.StringIO()
    yield I
    br = getattr(I, 'tclbridge', None)
    if br is not None and br.root is not None:
        try:
            br.root.destroy()
        except Exception:
            pass


def test_tk_version_and_font_system(interp):
    interp.run('v = Tk$FullVersion; fs = $FontSystem;', 'script')
    assert str(interp.get('v')).startswith(('8.6', '9.0'))     # Python 3.12 以前は Tk 8.6、3.14 は同梱の Tk 9.0
    assert interp.get('fs') == 'X11'


def test_simple_widgets_get_their_options(interp):
    interp.run('''
w = Window[];
b = Button[w, Text->"Hi"];
c = CheckButton[w, Text->"Chk"];
l = TextLabel[w, Text->"Lab"];
TclEval["update"];
info = TclEval["set r {}; foreach c [winfo children .] {lappend r [winfo class $c] [$c cget -text]}; return $r", String];
''', 'script')
    assert interp.get('info') == 'Button Hi Checkbutton Chk Label Lab'


def test_kbmainframe_menu_bar_structure(interp):
    """KBMainFrame は '.' の下にメニューバー (Frame: Menubutton 6 個 (KEK、File、Edit、Window、Task、Help) と Menu、仕切りの Frame、時計の Label) を作る (SAD の構成と同じ)"""
    interp.run('''
mf = KBMainFrame["TkTest", f, Title->"Test"];
TclEval["update"];
nmb = TclEval["set n 0; foreach c [winfo children [lindex [winfo children .] 0]] {if {[winfo class $c] eq {Menubutton}} {incr n}}; return $n", String];
''', 'script')
    assert str(interp.get('nmb')) == '6'


def test_kbframe_component_frame_with_label_widget(interp):
    """KBFComponentFrame (KBMainFrame.n) が LabelFrame のラベルにチェックボタンを埋め込む構成を作る。
    SAD (WSL の gs、Tk 8.6 X11) の同じ台本の構成 (Checkbutton 4 個 = Enable と A, B, C、Labelframe 2 個、Menubutton 6 個) と比べる"""
    interp.run('''
w = KBMainFrame["LabelFrameDemo", f, Title->"LabelFrame Demo"];
cb[0] = CheckButton[Text->"Enable", Variable:>cbv[0],
  Command:>(Scan[(cb[#][State] = If[cbv[0], "normal", "disabled"])&, Range[3]])];
cf = KBFComponentFrame[f, Add->{
  KBFHook[Command:>(cbv[0] = True)],
  KBFGroup[WidgetOption->{LabelWidget->cb[0]}],
  KBFCheckButton[Text->"A", Variable:>cbv[1], WidgetVariable:>cb[1]],
  KBFCheckButton[Text->"B", Variable:>cbv[2], WidgetVariable:>cb[2]],
  KBFCheckButton[Text->"C", Variable:>cbv[3], WidgetVariable:>cb[3]]}];
TclEval["update"];
cnt = TclEval["proc cnt {w c} {set n [expr {[winfo class $w] eq $c}]; foreach k [winfo children $w] {incr n [cnt $k $c]}; return $n}; list [cnt . Checkbutton] [cnt . Labelframe] [cnt . Menubutton]", String];
''', 'script')
    assert str(interp.get('cnt')) == '4 2 6'
    assert '???' not in interp.out.getvalue()


def test_widget_values_and_methods_like_gs(interp):
    """gs: w=Window[]・e=Entry[w] の値は Widget[…]、e[Insert]・lb[Insert] が働く (Tkinter.n の `a$:=Widget[a$,w[master$]]`)。
    ネストした規則の呼び出しで With[{a$=…}] が呼び出し元の同名の引数を見てしまい、ユーザーの w に遅延した値が付かなかった"""
    interp.run('''
w = Window[];
ev = "abc";
e = Entry[w, TextVariable:>ev, Width->20, Side->"top"];
lb = ListBox[w, Height->3, Side->"top"];
sv = 5;
s = Scale[w, From->0, To->10, Orient->"horizontal", Variable:>sv, Side->"top"];
hw = Head[w]; he = Head[e];
e[Insert] = {"end", "XY"};
lb[Insert] = {"end", "item1"};
lb[Insert] = {"end", "item2"};
TclEval["update"];
en = TclEval["return [[lindex [ff . Entry] 0] get]", String];
''', 'script')
    interp.run('''TclEval["proc ff {w cls} {set r {}; foreach c [winfo children $w] {if {[winfo class $c] eq $cls} {lappend r $c}; set r [concat $r [ff $c $cls]]}; return $r}"];
en = TclEval["return [[lindex [ff . Entry] 0] get]", String];
ls = TclEval["return [[lindex [ff . Listbox] 0] size]", String];''', 'script')
    assert interp.get('hw').name == 'Widget' and interp.get('he').name == 'Widget'
    assert interp.get('en') == 'abcXY'
    assert str(interp.get('ls')) in ('2', '2.0')


def test_canvas_widget_assignment_runs_the_package_upvalue(interp):
    """gs: `Canvas$Widget = c1` で Canvas.n の `(Canvas$Widget=Widget[c_,_Canvas])^:=…` が働き、描く先が c1 になる
    (既定の窓を作らない)。UpValue の照合が、遅延した定義を持つ記号 (The$CanvasObject := …) を評価してしまい既定の Canvas を作っていた"""
    interp.run('''
w = KBMainFrame["plot", plotframe, Title->"Trajectory Plot"];
c1 = Canvas[plotframe, Width->800, Height->500];
$DisplayFunction = CanvasDrawer;
Canvas$Widget = c1;
x = CanvasGraphics@CanvasObject[c1];
''', 'script')
    out = interp.out.getvalue()
    assert 'tclerror' not in out and 'DefaultCanvas' not in out
    assert interp.get('x') is not None and str(interp.get('x')) != 'Null'


def test_tkcanvaslastitem_keeps_counting_after_delete_all(interp):
    """K64 の TkCanvasLastItem は Tk の nextId-1。`delete all` のあとも id は続くので、Canvas.n の `StartID=CanvasLastID+1`
    (BelowFrame の `lower id StartID`) が古い 1 にならない (以前は ListContourPlot で tagOrId "1" doesn't match any items)"""
    interp.run('''
w = Window[];
c = Canvas[w, Width->100, Height->100];
p = TkPathName[c];
TclEval[p//" create line 0 0 1 1"];
TclEval[p//" create line 0 0 1 1"];
TclEval[p//" create line 0 0 1 1"];
l1 = TkCanvasLastItem[p];
TclEval[p//" delete all"];
l2 = TkCanvasLastItem[p];
TclEval[p//" create line 0 0 1 1"];
l3 = TkCanvasLastItem[p];
''', 'script')
    assert (interp.get('l1'), interp.get('l2'), interp.get('l3')) == (3, 3, 4)
