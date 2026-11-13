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
    assert str(interp.get('v')).startswith('8.6')
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
