# -*- coding: utf-8 -*-
"""v2.4.1: `$DisplayFunction` = MatplotlibDrawer (pysad/mpl_display.py)。Graphics[…] (Plot・ListPlot・OpticsPlot ほか) を matplotlib に描く。

Agg (窓を出さない・Tk 不要)。Graphics の Line・Point・Rectangle・Polygon・Circle・Text の個数と座標が、matplotlib 側のデータ
(Line2D・Patch・Text) と一致することを数値で確かめる。切替 (環境変数 PYSAD_DISPLAY)・保存 (PYSAD_DISPLAY_FILE、2 回目以降は番号)・
matplotlib が無いとき・未対応の要素の警告も試験する。SAD の Packages (ListPlot.n・Plots.n・OpticsPlot.n) が要る試験は
PYSAD_SAD_PACKAGES が無ければ skip。
"""
import io
import math
import os

import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from pysad import draw_mpl as dm                      # noqa: E402
from pysad import mpl_display as md                   # noqa: E402
from pysad.ffscmd import run_sad_string                # noqa: E402
from pysad.sadscript import Expr, Rule, Sym            # noqa: E402

PK = os.environ.get('PYSAD_SAD_PACKAGES')
needs_pk = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')

FODO = '''MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.2 K2=0);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
'''
PNG = b'\x89PNG\r\n\x1a\n'


def _run(body, env=None, exit_=True):
    """台本を実行して (ffs, I, 出力) を返す。env は一時的な環境変数 (None で削除)"""
    env = dict(env or {})
    env.setdefault('PYSAD_NO_DISPLAY', '1')
    old = {k: os.environ.get(k) for k in env}
    for k, v in env.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    try:
        buf = io.StringIO()
        ffs, I = run_sad_string(body + ('Exit;\n' if exit_ else ''), out=buf)
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    return ffs, I, buf.getvalue()


def _mpl(tmp_path, body, name='a.png', extra=None):
    env = {'PYSAD_DISPLAY': 'mpl', 'PYSAD_DISPLAY_FILE': str(tmp_path / name)}
    env.update(extra or {})
    return _run(body if 'FFS;' in body else 'FFS;\n' + body, env)


def _opt(g, name):
    for r in g.args[1]:
        if isinstance(r, Rule) and getattr(r.lhs, 'name', None) == name:
            return r.rhs
    return None


def _count(e, head):
    n = 0
    if isinstance(e, Expr):
        n += 1 if str(getattr(e.head, 'name', e.head)) == head else 0
        for a in e.args:
            n += _count(a, head)
    elif isinstance(e, list):
        for a in e:
            n += _count(a, head)
    return n


def _frames():
    return md.last_figure()._pysad_display['frames']


def _xy(line):
    return list(line.get_xdata()), list(line.get_ydata())


# --------------------------------------------------------------------------------------------- 色・文字・幾何 (SAD 不要)
def test_sad_color():
    import matplotlib.colors as mc
    assert md.sad_color('dodger blue') == pytest.approx(mc.to_rgb('dodgerblue'))
    assert md.sad_color('forest green') == pytest.approx(mc.to_rgb('forestgreen'))
    t = mc.to_rgb('tomato')
    assert md.sad_color('tomato3') == pytest.approx(tuple(v * 0.804 for v in t))      # X11 の末尾の数字 3 = 80.4%
    assert md.sad_color('yellow4') == pytest.approx((0.545, 0.545, 0.0))
    assert md.sad_color('gray40') == pytest.approx((0.4, 0.4, 0.4))
    assert md.sad_color('#ff8000') == pytest.approx((1.0, 128 / 255, 0.0))
    assert md.sad_color('#f80') == pytest.approx((1.0, 136 / 255, 0.0))
    assert md.sad_color('SystemButtonText') == pytest.approx((0, 0, 0))
    assert md.sad_color('') is None
    assert md.sad_color(Sym('Null')) is None
    assert md.sad_color('nosuchcolorname') is None


def test_sad_text_markup():
    assert md.sad_text('plain') == 'plain'
    assert md.sad_text('`fh`n`dx`n') == 'η$_{x}$'                      # η_x
    assert md.sad_text('`f\xd6`ob`n`dx`O`n') == '√β$_{x}$'            # √β_x
    assert md.sad_text('`f\xd6`ob`n`dx`O`n, `f\xd6`ob`n`dy`O`n (`f\xd6`n`om`O)') == '√β$_{x}$, √β$_{y}$ (√m)'
    assert md.sad_text('10`u3') == '10$^{3}$'
    assert md.sad_text('`fs`n`dxx`n') == 'σ$_{xx}$'
    bs = chr(92)
    # LaTeX の文字 (FitPlot の題): 先頭の記号と centering を取り、二重のバックスラッシュを改行に。$…$ はそのまま
    assert md.sad_text(bs + '@' + bs + 'centering Fit: a,' + bs * 2 + ' $x^2$') == 'Fit: a,' + chr(10) + ' $x^2$'
    # mathtext に読めない $ は打ち消す (描画で落ちない)
    assert md.sad_text('cost $5') == 'cost ' + bs + '$5'


def test_make_window_matches_canvas_n():
    # MakeWindow[Automatic, Automatic]: 窓の幅 9.6 (原点 2)、縦横比 GoldenRatio、窓の中心は生の窓の中心
    x1, x2, y1, y2 = md.make_window([[0, 1], [0, 1]], None, None)
    assert (x2 - x1) / (y2 - y1) == pytest.approx(md.GOLDEN)
    assert (x1, x2) == pytest.approx((2.0, 11.6))
    assert (y1 + y2) / 2 == pytest.approx(1.8 + md.RAW[1] / 2)
    # 縦横比 1: 高さに合わせて幅が縮み、中心は変わらない
    a1, a2, b1, b2 = md.make_window([[0, 1], [0, 1]], None, 1.0)
    assert a2 - a1 == pytest.approx(b2 - b1)
    assert (a1 + a2) / 2 == pytest.approx(6.8)
    # rect (PlotRegion の窓に対する割合) の枠: OpticsPlot の 2 枠目 {{0,0.1},{1,0.6}} は高さ 0.5 倍
    r = md.make_window([[0, 1], [0, 1]], [[0.0, 0.1], [1.0, 0.6]], None)
    assert r[0] == pytest.approx(2.0) and r[1] == pytest.approx(11.6)
    assert r[2] == pytest.approx(1.8 + 0.1 * md.RAW[1]) and r[3] == pytest.approx(1.8 + 0.6 * md.RAW[1])


def test_frame_range_matches_canvas_n():
    # 両方 Automatic: データ範囲の外に 7% の余白
    lo, hi = md.frame_range('Linear', (None, None), (1.0, 5.0))
    assert (lo, hi) == pytest.approx((1.0 - 0.07 * 4.0, 5.0 + 0.07 * 4.0))
    assert md.frame_range('Linear', (0.0, 10.0), (1.0, 5.0)) == (0.0, 10.0)
    lo, hi = md.frame_range('Linear', (None, 10.0), (1.0, 5.0))          # 下だけ Automatic: min(b1, r2) と r2 から
    assert (lo, hi) == pytest.approx((1.0 - 0.07 * 9.0, 10.0))
    lo, hi = md.frame_range('Log', (None, None), (1.0, 100.0))           # Log は log10 で 7%
    assert (lo, hi) == pytest.approx((0.0 - 0.14, 2.0 + 0.14))
    assert md.frame_range('Linear', (None, None), (0.0, 0.0)) == (-0.001, 0.001)
    lo, hi = md.frame_range('Linear', (3.0, 3.0), (0.0, 1.0))
    assert lo < 3.0 < hi


def test_display_mode_from_env(monkeypatch):
    for v, want in (('mpl', 'mpl'), ('MatplotLib', 'mpl'), ('tk', 'tk'), ('none', 'none'), ('', None), ('xyz', None)):
        monkeypatch.setenv('PYSAD_DISPLAY', v)
        assert md.display_mode_from_env() == want
    monkeypatch.delenv('PYSAD_DISPLAY')
    assert md.display_mode_from_env() is None


# --------------------------------------------------------------------------------------------- 切替 ($DisplayFunction)
@needs_pk
def test_display_function_switch():
    _, I, _ = _run('FFS;\n', {'PYSAD_DISPLAY': 'mpl'})
    assert I.globals['$DisplayFunction'] == Sym('MatplotlibDrawer')
    _, I, _ = _run('FFS;\n', {'PYSAD_DISPLAY': 'none'})
    assert I.globals['$DisplayFunction'] == Sym('Identity')
    _, I, _ = _run('FFS;\n', {'PYSAD_DISPLAY': None})                  # 従来どおり (PYSAD_NO_DISPLAY=1 は Identity)
    assert I.globals['$DisplayFunction'] == Sym('Identity')
    assert 'MatplotlibDrawer' in I.builtins                           # 組み込みの関数はいつも登録される


@needs_pk
def test_identity_makes_no_file(tmp_path):
    _, I, txt = _run('FFS;\ng=Plot[Sin[x],{x,0,6}];\nPrint["head=",Head[g]];\n',
                     {'PYSAD_DISPLAY_FILE': str(tmp_path / 'z.png')})
    assert 'head=Graphics' in txt
    assert not list(tmp_path.iterdir())


# --------------------------------------------------------------------------------------------- Plot・ListPlot
@needs_pk
def test_plot_sin_matches_graphics(tmp_path):
    ffs, I, txt = _mpl(tmp_path, 'g=Plot[Sin[x],{x,0,6}];\nPrint["head=",Head[g]];\n')
    g = I.globals['g']
    assert 'head=Graphics' in txt
    # ListPlot.n 自身が Show を 2 回呼ぶので、同じ内容は 1 つにまとめる: 保存は 1 枚
    assert sorted(p.name for p in tmp_path.iterdir()) == ['a.png']
    assert (tmp_path / 'a.png').read_bytes()[:8] == PNG
    (fr,) = [f for f in _frames() if f.ax is not None]
    pts = dm.graphics_lines(g)
    assert len(pts) == 1 and len(fr.rec['lines']) == 1
    xs, ys = _xy(fr.rec['lines'][0])
    assert len(xs) == len(pts[0]) == 505
    assert max(abs(a - p[0]) for a, p in zip(xs, pts[0])) < 1e-12
    assert max(abs(a - p[1]) for a, p in zip(ys, pts[0])) < 1e-12
    dr = _opt(g, 'DataRange')
    x0, x1 = dr[0]
    assert fr.ax.get_xlim() == pytest.approx((x0 - 0.07 * (x1 - x0), x1 + 0.07 * (x1 - x0)))      # PlotRange Automatic → 7% の余白
    assert fr.ax.get_ylim() == pytest.approx(tuple(_opt(g, 'PlotRange')))                          # y は Plot が決めた PlotRange
    assert fr.rec['lines'][0].get_color() == pytest.approx(md.sad_color('dodger blue'))
    assert fr.rec['lines'][0].get_linewidth() == pytest.approx(2.0 * md.LW_PT)                     # Thickness->2
    # 枠: 既定は FrameTicks {True,True,Ticks,Ticks}、格子 (Automatic) は主目盛りの線
    ax = fr.ax
    assert ax.xaxis.get_major_ticks()[0].label1.get_visible() and ax.yaxis.get_major_ticks()[0].label1.get_visible()
    assert ax.xaxis.get_major_ticks()[0].tick2line.get_visible()          # 上の辺の目盛りの線 (Ticks)
    assert not ax.xaxis.get_major_ticks()[0].label2.get_visible()
    assert any(gl.get_visible() for gl in ax.get_xgridlines())


@needs_pk
def test_second_plot_gets_numbered_file(tmp_path):
    _mpl(tmp_path, 'a1=Plot[Sin[x],{x,0,6}];\na2=Plot[Cos[x],{x,0,6}];\na3=Plot[Cos[x],{x,0,6}];\n')
    # 1 つ目、2 つ目は別の内容で 2 枚。3 つ目は 2 つ目と同じ内容が続くので出さない
    assert sorted(p.name for p in tmp_path.iterdir()) == ['a.png', 'a_2.png']


@needs_pk
def test_two_curves_and_options(tmp_path):
    ffs, I, txt = _mpl(tmp_path, 'g=Plot[{Sin[y],Cos[y]},{y,0,6},PlotLabel->"two curves",FrameLabel->{"y axis","value"}];\n')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 2
    pts = dm.graphics_lines(I.globals['g'])
    assert len(pts) == 2
    for ln, p in zip(fr.rec['lines'], pts):
        xs, ys = _xy(ln)
        assert max(abs(a - q[0]) for a, q in zip(xs, p)) < 1e-12 and max(abs(a - q[1]) for a, q in zip(ys, p)) < 1e-12
    assert fr.ax.get_xlabel() == 'y axis' and fr.ax.get_ylabel() == 'value'
    assert [t.get_text() for t in md.last_figure().texts] == ['two curves']
    assert fr.rec['lines'][0].get_color() != fr.rec['lines'][1].get_color()                       # PlotColor の表の別の色
    assert fr.rec['lines'][1].get_linestyle() == '--'                                              # 2 本目は破線 (Dashing)


@needs_pk
def test_listplot_joined_and_points(tmp_path):
    ffs, I, txt = _mpl(tmp_path, 'g=ListPlot[{1,4,9,16},PlotJoined->True];\n')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 1 and len(fr.rec['points']) == 1                  # ListPlot の Plot->True (既定) で、折れ線に点も描く
    assert _xy(fr.rec['lines'][0]) == ([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0])
    assert _xy(fr.rec['points'][0]) == ([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0])
    assert fr.ax.get_xlim() == pytest.approx((1 - 0.07 * 3, 4 + 0.07 * 3))
    assert fr.ax.get_ylim() == pytest.approx((1 - 0.07 * 15, 16 + 0.07 * 15))
    _mpl(tmp_path, 'g=ListPlot[{1,4,9,16},PlotJoined->True,Plot->False];\n', 'a0.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 1 and not fr.rec['points']
    # 点だけ (PlotJoined 既定 False): Point[…] → 記号つきの点列。線は無い
    ffs, I, txt = _mpl(tmp_path, 'g=ListPlot[{{1,2},{2,5},{3,3}},PointSymbol->"SQ"];\n', 'b.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert not fr.rec['lines'] and len(fr.rec['points']) == 1
    ln = fr.rec['points'][0]
    assert _xy(ln) == ([1.0, 2.0, 3.0], [2.0, 5.0, 3.0])
    assert ln.get_marker() == 's' and ln.get_linestyle() == 'None'
    # 点つきの折れ線 Plot->True
    ffs, I, txt = _mpl(tmp_path, 'g=ListPlot[{{1,2},{2,5}},PlotJoined->True,Plot->True,PointSymbol->"UT",Thickness->3];\n', 'c.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 1 and len(fr.rec['points']) == 1
    assert fr.rec['points'][0].get_marker() == '^' and fr.rec['lines'][0].get_linewidth() == pytest.approx(3 * md.LW_PT)


@needs_pk
def test_listplot_options_range_scale_grid_frame(tmp_path):
    # PlotRange 指定・格子なし・Log
    _mpl(tmp_path, 'g=ListPlot[{{1,2},{2,5},{3,3}},PlotJoined->True,PlotRange->{{0,4},{0,7}},GridLines->None];\n')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert fr.ax.get_xlim() == (0.0, 4.0) and fr.ax.get_ylim() == (0.0, 7.0)
    assert not any(gl.get_visible() for gl in fr.ax.get_xgridlines() + fr.ax.get_ygridlines())
    _mpl(tmp_path, 'g=ListPlot[{{1,1},{10,100},{100,10000}},Scale->{Log,Log},PlotJoined->True];\n', 'b.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert fr.ax.get_xscale() == 'log' and fr.ax.get_yscale() == 'log'
    lo, hi = fr.ax.get_xlim()
    assert lo == pytest.approx(10 ** (0 - 0.14)) and hi == pytest.approx(10 ** (2 + 0.14))       # log10 の 7% の余白
    assert _xy(fr.rec['lines'][0])[1] == [1.0, 100.0, 10000.0]
    # Frame->False: 枠・目盛り・格子なし
    _mpl(tmp_path, 'g=ListPlot[{1,4,9},PlotJoined->True,Frame->False];\n', 'c.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert not any(sp.get_visible() for sp in fr.ax.spines.values())
    assert not fr.ax.xaxis.get_major_ticks()[0].label1.get_visible()
    assert not fr.ax.xaxis.get_major_ticks()[0].tick1line.get_visible()
    # FrameTicks {Ticks, True, Ticks, Ticks}: 下の文字なし、左の文字あり、(OpticsPlot の上の枠と同じ)
    _mpl(tmp_path, 'g=ListPlot[{1,4,9},PlotJoined->True,FrameTicks->{Ticks,True,Ticks,Ticks}];\n', 'd.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert not fr.ax.xaxis.get_major_ticks()[0].label1.get_visible()
    assert fr.ax.yaxis.get_major_ticks()[0].label1.get_visible()
    assert fr.ax.yaxis.get_major_ticks()[0].tick2line.get_visible()
    # FrameLabel の 4 辺と PlotLabel
    _mpl(tmp_path, 'g=ListPlot[{1,4,9},PlotJoined->True,FrameLabel->{"bottom","left","top","right"},PlotLabel->"title"];\n', 'e.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert fr.ax.get_xlabel() == 'bottom' and fr.ax.get_ylabel() == 'left'
    ann = [t.get_text() for t in fr.ax.texts]
    assert 'top' in ann and 'right' in ann
    assert [t.get_text() for t in md.last_figure().texts] == ['title']


@needs_pk
def test_listplot_errorbars_legend_and_step(tmp_path):
    _mpl(tmp_path, 'g=ListPlot[{{1,2,0.3,0.5},{2,5,0.2,0.4}},PlotRange->{{0,3},{0,7}},Legend->"data"];\n')
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['errorbars']) == 1
    segs = fr.rec['errorbars'][0].lines[2]                          # (横棒, 縦棒) の線分
    xerr = segs[0].get_segments()
    assert [(round(s[0][0], 12), round(s[1][0], 12)) for s in xerr] == [(0.7, 1.3), (1.8, 2.2)]
    yerr = segs[1].get_segments()
    assert [(round(s[0][1], 12), round(s[1][1], 12)) for s in yerr] == [(1.5, 2.5), (4.6, 5.4)]
    lg = fr.rec['legends'][0]
    assert [t.get_text() for t in lg.get_texts()] == ['data']
    # Step: StepRatio 0.5 → 横・縦の順の階段 (Canvas$StepData)
    _mpl(tmp_path, 'g=ListPlot[{{0,1},{2,3},{4,2}},PlotJoined->Step,StepRatio->0.5];\n', 'b.png')
    (fr,) = [f for f in _frames() if f.ax is not None]
    xs, ys = _xy(fr.rec['lines'][0])
    assert xs == [0.0, 1.0, 1.0, 3.0, 3.0, 4.0] and ys == [1.0, 1.0, 3.0, 3.0, 2.0, 2.0]


@needs_pk
def test_histoplot_bars(tmp_path):
    ffs, I, txt = _mpl(tmp_path, 'g=HistoPlot[Table[Sin[i*1.7]+Cos[i*0.3],{i,1,60}]];\n')
    g = I.globals['g']
    n = _count(g, 'Rectangle')
    frs = [f for f in _frames() if f.ax is not None]
    assert n > 5 and sum(len(f.rec['rects']) for f in frs) == n
    r = frs[0].rec['rects'][0]
    assert r.get_facecolor()[:3] == pytest.approx(md.sad_color('royal blue'))


# --------------------------------------------------------------------------------------------- OpticsPlot・GeometryPlot
@needs_pk
def test_opticsplot_frames_lines_boxes_names(tmp_path):
    ffs, I, txt = _mpl(tmp_path, FODO + 'USE RING; CELL; CALC;\ng=OpticsPlot[{{"BX","BY"},{"EX"}},PlotLabel->"FODO"];\n')
    g = I.globals['g']
    assert (tmp_path / 'a.png').read_bytes()[:8] == PNG and sorted(p.name for p in tmp_path.iterdir()) == ['a.png']
    frames = _frames()
    assert len(frames) == 4                                            # 外側の空の枠 + BX・BY + EX + ラティス
    assert frames[0].ax is None
    f_bxby, f_ex, f_lat = frames[1], frames[2], frames[3]
    # 点列: Graphics の Line と全点一致 (BX・BY は √β、EX は mm に直したあと)
    gl = dm.graphics_lines(g)
    assert len(gl) == 3
    lines = f_bxby.rec['lines'] + f_ex.rec['lines']
    assert len(lines) == 3
    for ln, p in zip(lines, gl):
        xs, ys = _xy(ln)
        assert len(xs) == len(p)
        assert max(abs(a - q[0]) for a, q in zip(xs, p)) < 1e-12
        assert max(abs(a - q[1]) for a, q in zip(ys, p)) < 1e-12 * max(1.0, max(abs(q[1]) for q in p))
    # 枠の積み方: 上から BX・BY、EX、ラティス。x の範囲・窓の左右は同じ
    pos = [f.ax.get_position() for f in (f_bxby, f_ex, f_lat)]
    assert pos[0].y0 == pytest.approx(pos[1].y1) and pos[1].y0 > pos[2].y1
    assert all(p.x0 == pytest.approx(pos[0].x0) and p.x1 == pytest.approx(pos[0].x1) for p in pos)
    xl = f_bxby.ax.get_xlim()
    circ = float(ffs.line_value('S', ffs.nlat, 0.0)) - float(ffs.line_value('S', 1, 0.0))
    assert xl == pytest.approx((-0.001, circ + 0.001))                                              # PlotRange は {-0.001, 周長 + 0.001}
    assert f_ex.ax.get_xlim() == pytest.approx(xl) and f_lat.ax.get_xlim() == pytest.approx(xl)
    # 凡例 (2 本 + 1 本) と軸の文字
    assert [t.get_text() for t in f_bxby.rec['legends'][0].get_texts()] == ['√β$_{x}$', '√β$_{y}$']
    assert [t.get_text() for t in f_ex.rec['legends'][0].get_texts()] == ['η$_{x}$']
    assert f_ex.ax.get_ylabel() == 'η$_{x}$ (mm)'
    assert f_bxby.ax.get_ylabel() == '√β$_{x}$, √β$_{y}$ (√m)'
    assert [t.get_text() for t in md.last_figure().texts] == ['FODO']
    # 下の枠 (EX) だけ x の文字、上の枠は目盛りの線だけ
    assert not f_bxby.ax.xaxis.get_major_ticks()[0].label1.get_visible()
    assert f_ex.ax.xaxis.get_major_ticks()[0].label1.get_visible()
    # 配色: BX・BY は別の色 (PlotColorDefault の表の順)
    assert f_bxby.rec['lines'][0].get_color() == pytest.approx(md.sad_color('royal blue'))
    assert f_bxby.rec['lines'][1].get_color() == pytest.approx(md.sad_color('tomato3'))
    # ラティスの箱: Graphics の Rectangle (中に Graphics を持たないもの) と個数・座標が一致
    boxes = dm.graphics_rectangles(g)
    assert len(f_lat.rec['rects']) == len(boxes) and len(boxes) > 50
    for r, (x0, y0, x1, y1) in zip(f_lat.rec['rects'], boxes):
        assert r.get_x() == pytest.approx(min(x0, x1)) and r.get_y() == pytest.approx(min(y0, y1))
        assert r.get_width() == pytest.approx(abs(x1 - x0)) and r.get_height() == pytest.approx(abs(y1 - y0))
    colors = {tuple(round(c, 3) for c in r.get_facecolor()[:3]) for r in f_lat.rec['rects'] if r.get_facecolor()[3] > 0}
    assert tuple(round(c, 3) for c in md.sad_color('limegreen')) in colors and tuple(round(c, 3) for c in md.sad_color('gold')) in colors
    # 名前: Graphics の Text の部分集合 (RemoveOverlap で重なるものを消す)。残った名前は互いに重ならない
    names = dm.graphics_texts(g)
    kept = [(t.get_text(), t.get_position()[0]) for t in f_lat.rec['texts']]
    assert len(kept) + len(f_lat.rec['texts_removed']) == len(names)
    assert 0 < len(kept) < len(names)
    it = iter(names)
    assert all(any(k[0] == n[0] and k[1] == pytest.approx(n[1]) for n in it) for k in kept)        # 順序つきの部分列
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    rend = FigureCanvasAgg(md.last_figure()).get_renderer()
    bbs = [t.get_window_extent(rend) for t in f_lat.rec['texts']]
    for i in range(len(bbs)):
        for j in range(i + 1, len(bbs)):
            assert not bbs[i].overlaps(bbs[j]), (kept[i], kept[j])
    assert all(t.get_rotation() == 270 and t.get_ha() == 'center' and t.get_va() == 'top' for t in f_lat.rec['texts'])


@needs_pk
def test_opticsplot_region_and_unit_label(tmp_path):
    ffs, I, txt = _mpl(tmp_path, FODO + 'USE RING; CELL; CALC;\ng=OpticsPlot[{{"BX"}},Region->{1,20}];\n')
    frames = [f for f in _frames() if f.ax is not None]
    assert len(frames) == 2                                            # BX の枠 + ラティス
    g = I.globals['g']
    p = dm.graphics_lines(g)[0]
    xs, ys = _xy(frames[0].rec['lines'][0])
    assert xs == [q[0] for q in p] and ys == [q[1] for q in p]
    ann = [t.get_text() for t in frames[0].ax.texts]
    assert 'm' in ann                                                   # FrameUnitLabel (s の単位) は軸の右下


@needs_pk
def test_geometryplot_lines_and_polygons(tmp_path):
    ffs, I, txt = _mpl(tmp_path, FODO + 'USE RING; CELL; CALC;\nq=GeometryPlot[];\n')
    q = I.globals['q']
    frs = [f for f in _frames() if f.ax is not None]
    assert len(frs) == 1
    nl, npoly = _count(q, 'Line'), _count(q, 'Polygon')
    assert nl >= 1 and npoly > 10                                                   # 12 セルのリング: 基準軌道の線 1 本 + 磁石の多角形 73
    assert len(frs[0].rec['lines']) == nl and len(frs[0].rec['polys']) == npoly
    assert not any(sp.get_visible() for sp in frs[0].ax.spines.values())      # Frame->False
    gl = dm.graphics_lines(q)
    for ln, p in zip(frs[0].rec['lines'], gl):
        xs, ys = _xy(ln)
        assert xs == [r[0] for r in p] and ys == [r[1] for r in p]


# --------------------------------------------------------------------------------------------- 手で組んだ Graphics の各要素
SHOW = 'Show[Graphics[{%s},PlotRange->{{0,10},{0,5}}%s],DisplayFunction->MatplotlibDrawer];\n'


@needs_pk
def test_primitives_line_polygon_circle_rectangle_text(tmp_path):
    prims = ('Line[{{0,0},{5,5},{10,0}},PlotColor->"red",Thickness->3,Dashing->{0.2,0.1},PlotJoined->True,Plot->False],'
             'Polygon[{{1,1},{3,1},{2,3}},FillColor->"light blue",PlotColor->"navy"],'
             'Circle[{7,3},0.5,FillColor->"yellow4",PlotColor->"black"],'
             'Rectangle[{4,3},{6,4.5},FillColor->"pink",PlotColor->"gray40"],'
             'Text[{"hello",{Scaled[2],Scaled[4]}},TextSize->1.5,Anchor->"nw",PlotColor->"blue",TextRotate->30],'
             'Text[{"east",{Scaled[9],Scaled[4]}},Anchor->"e"],'
             'Text[{"win",{WindowScaled[2.0],WindowScaled[1.8]}},Anchor->"sw"],'
             'Point[{{1,4},{2,4.5}},PointSymbol->"DT",PointColor->"green"]')
    ffs, I, txt = _mpl(tmp_path, SHOW % (prims, ''))
    assert not [ln for ln in txt.splitlines() if 'MatplotlibDrawer:' in ln and '保存しました' not in ln]          # 警告なし
    (fr,) = [f for f in _frames() if f.ax is not None]
    ax = fr.ax
    assert ax.get_xlim() == (0.0, 10.0) and ax.get_ylim() == (0.0, 5.0)
    (ln,) = fr.rec['lines']
    assert _xy(ln) == ([0.0, 5.0, 10.0], [0.0, 5.0, 0.0])
    assert ln.get_color() == pytest.approx(md.sad_color('red')) and ln.get_linewidth() == pytest.approx(3 * md.LW_PT)
    on, off = ln._dash_pattern[1]
    assert on == pytest.approx(0.2 * md.PT_PER_INCH) and off == pytest.approx(0.1 * md.PT_PER_INCH)
    (poly,) = fr.rec['polys']
    assert [tuple(v) for v in poly.get_xy()[:3]] == [(1.0, 1.0), (3.0, 1.0), (2.0, 3.0)]
    assert poly.get_facecolor()[:3] == pytest.approx(md.sad_color('light blue'))
    assert poly.get_edgecolor()[:3] == pytest.approx(md.sad_color('navy'))
    (ov,) = fr.rec['ovals']
    assert ov.center == (7.0, 3.0) and ov.width == pytest.approx(1.0) and ov.height == pytest.approx(1.0)
    assert ov.get_facecolor()[:3] == pytest.approx(md.sad_color('yellow4'))
    (rc,) = fr.rec['rects']
    assert (rc.get_x(), rc.get_y(), rc.get_width(), rc.get_height()) == (4.0, 3.0, 2.0, 1.5)
    assert rc.get_edgecolor()[:3] == pytest.approx(md.sad_color('gray40'))
    t_hello, t_east, t_win = fr.rec['texts']
    assert t_hello.get_text() == 'hello' and t_hello.get_position() == (2.0, 4.0)
    assert (t_hello.get_ha(), t_hello.get_va()) == ('left', 'top')             # Tk の anchor nw: 文字の北西の角が点 (文字は点の右下)
    assert t_hello.get_rotation() == 30 and t_hello.get_color() == pytest.approx(md.sad_color('blue'))
    assert t_hello.get_fontsize() == pytest.approx(1.5 * md.FONT_TEXT)
    assert (t_east.get_ha(), t_east.get_va()) == ('right', 'center')
    assert (t_win.get_ha(), t_win.get_va()) == ('left', 'bottom')
    # WindowScaled[2.0], [1.8] は生の窓の左下 (原点 {2, 1.8}) = 窓の左下の隅
    px, py = ax.transAxes.inverted().transform(t_win.get_transform().transform(t_win.get_position()))
    assert px == pytest.approx(0.0) and py == pytest.approx(0.0)
    (pts,) = fr.rec['points']
    assert _xy(pts) == ([1.0, 2.0], [4.0, 4.5]) and pts.get_marker() == 'v'


@needs_pk
def test_anchor_table_matches_tk():
    # Tk の anchor は「文字のどの点を座標に置くか」。n = 文字の上辺の中央が点 (文字は点の下)
    assert md._ANCHOR['n'] == ('center', 'top') and md._ANCHOR['s'] == ('center', 'bottom')
    assert md._ANCHOR['w'] == ('left', 'center') and md._ANCHOR['e'] == ('right', 'center')
    assert md._ANCHOR['c'] == ('center', 'center') and md._ANCHOR['sw'] == ('left', 'bottom')


@needs_pk
def test_clip_and_background_and_frame_color(tmp_path):
    prims = 'Line[{{0,0},{12,6}},PlotColor->"black"],Line[{{0,5},{12,5}},Clip->False]'
    _mpl(tmp_path, SHOW % (prims, ',Background->"#eeeeff",FrameColor->"red"'))
    (fr,) = [f for f in _frames() if f.ax is not None]
    l1, l2 = fr.rec['lines']
    assert l1.get_clip_on() and not l2.get_clip_on()                                # Clip->False は枠で切らない
    assert fr.ax.patch.get_facecolor()[:3] == pytest.approx((238 / 255, 238 / 255, 1.0))
    assert fr.ax.spines['left'].get_edgecolor()[:3] == pytest.approx((1.0, 0.0, 0.0))


@needs_pk
def test_unsupported_elements_warn_and_skip(tmp_path):
    prims = 'Line[{{0,0},{5,5}}],Foo[1,2],Hook[Command->Identity],Line[{{0,1},{5,1}},Arrow->1],"stray string"'
    ffs, I, txt = _mpl(tmp_path, SHOW % (prims, ''))
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 2                                                # Foo・Hook・文字列は飛ばす (Arrow つきの線は描く)
    assert 'Foo[…] は未対応' in txt and 'Arrow' in txt
    assert 'Hook' not in txt                                                        # Hook (対話の束縛) は黙って飛ばす
    w = md.last_figure()._pysad_display['warnings']
    assert any('Foo' in m for m in w) and any('Arrow' in m for m in w) and any('str' in m for m in w)
    assert txt.count('Foo[…]') == 1                                                 # 同じ種類の警告は 1 回


@needs_pk
def test_unsupported_scale_text_position_and_tags(tmp_path):
    prims = 'Line[{{0,0},{5,5}}],Text[{"x",{LeftTags["a",1],Scaled[1]}}]'
    ffs, I, txt = _mpl(tmp_path, SHOW % (prims, ',Scale->{Date,Linear}'))
    assert 'Scale->Date は未対応' in txt and 'LeftTags' in txt
    (fr,) = [f for f in _frames() if f.ax is not None]
    assert len(fr.rec['lines']) == 1 and not fr.rec['texts']


# --------------------------------------------------------------------------------------------- 保存・出力
@needs_pk
@pytest.mark.parametrize('ext,sig', [('png', PNG), ('pdf', b'%PDF'), ('svg', b'<?xml')])
def test_save_formats(tmp_path, ext, sig):
    _mpl(tmp_path, 'g=ListPlot[{1,4,9},PlotJoined->True];\n', 'f.' + ext)
    assert (tmp_path / ('f.' + ext)).read_bytes()[:len(sig)] == sig


@needs_pk
def test_file_numbering_with_braces_and_sad_option(tmp_path):
    _mpl(tmp_path, 'a=ListPlot[{1,4,9},PlotJoined->True];\nb=ListPlot[{2,4,8},PlotJoined->True];\n', 'p{n}.png')
    assert sorted(p.name for p in tmp_path.iterdir()) == ['p1.png', 'p2.png']
    out = tmp_path / 'sub' / 'sad_option.png'
    _run('FFS;\nMatplotlibDrawer$Option={File->"%s",DPI->50};\ng=ListPlot[{1,4,9},PlotJoined->True];\n' % out.as_posix(),
         {'PYSAD_DISPLAY': 'mpl'})
    assert out.read_bytes()[:8] == PNG                                              # SAD 側の File->… (フォルダは作る)


@needs_pk
def test_bad_extension_is_reported(tmp_path):
    ffs, I, txt = _mpl(tmp_path, 'g=ListPlot[{1,4,9},PlotJoined->True];\nPrint["after"];\n', 'f.xyz')
    assert '保存できません' in txt and 'after' in txt
    assert not list(tmp_path.iterdir())


@needs_pk
def test_without_file_agg_shows_nothing(tmp_path):
    ffs, I, txt = _run('FFS;\ng=ListPlot[{1,4,9},PlotJoined->True];\nPrint["done"];\n',
                       {'PYSAD_DISPLAY': 'mpl', 'PYSAD_DISPLAY_FILE': None})
    assert 'done' in txt and '保存しました' not in txt
    assert md.last_figure() is not None and not list(tmp_path.iterdir())


@needs_pk
def test_python_api_render(tmp_path):
    ffs, I, txt = _run('FFS;\ng=ListPlot[{1,4,9,16},PlotJoined->True,DisplayFunction->Identity];\n', exit_=False)
    g = I.globals['g']
    fig = md.render(I, g)
    (fr,) = [f for f in fig._pysad_display['frames'] if f.ax is not None]
    assert _xy(fr.rec['lines'][0]) == ([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0])
    assert I.globals['$DisplayFunction'] == Sym('Identity')                          # もとの $DisplayFunction は戻す
    fig2 = md.render(I, g, save=str(tmp_path / 'r.png'))                              # 同じ式をもう一度: 指紋で捨てない
    assert fig2 is not None and (tmp_path / 'r.png').read_bytes()[:8] == PNG
    with pytest.raises(Exception):
        md.render(I, 123)


@needs_pk
def test_matplotlib_missing_warns_only(monkeypatch, tmp_path):
    monkeypatch.setattr(dm, '_import_mpl', lambda: None)
    ffs, I, txt = _mpl(tmp_path, 'a=ListPlot[{1,4,9},PlotJoined->True];\nb=Plot[Sin[x],{x,0,6}];\nPrint["still running ",Head[b]];\n')
    assert txt.count('matplotlib がありません') == 1                                   # 警告は 1 回だけ。評価は続く
    assert 'still running Graphics' in txt and not list(tmp_path.iterdir())
    with pytest.raises(ImportError):
        md.render(I, I.globals['a'])


@needs_pk
def test_draw_command_still_works_with_display_mpl(tmp_path):
    # DRAW は draw_mpl のまま (PYSAD_DISPLAY=mpl でも保存先は PYSAD_DRAW_FILE)
    ffs, I, txt = _run(FODO + 'USE RING; CELL; CALC;\nDRAW BX BY;\n',
                       {'PYSAD_DISPLAY': 'mpl', 'PYSAD_DRAW_FILE': str(tmp_path / 'd.png')})
    assert (tmp_path / 'd.png').read_bytes()[:8] == PNG
    assert 'MatplotlibDrawer' not in txt


@needs_pk
def test_global_display_function_assignment(tmp_path):
    # 環境変数なしでも、台本で $DisplayFunction=MatplotlibDrawer; と書けば使える
    ffs, I, txt = _run('FFS;\n$DisplayFunction=MatplotlibDrawer;\ng=ListPlot[{1,4,9},PlotJoined->True];\n',
                       {'PYSAD_DISPLAY_FILE': str(tmp_path / 's.png')})
    assert (tmp_path / 's.png').read_bytes()[:8] == PNG
    assert math.isfinite(_xy(_frames()[0].rec['lines'][0])[1][0])
