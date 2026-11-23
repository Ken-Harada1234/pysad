# -*- coding: utf-8 -*-
"""v2.4.1: FFS の DRAW (matplotlib)。pysad/draw_mpl.py。

- 構文解析 (CANVASDRAW.n の SeparateAmpsnd・r1 r2・& ・LAT・要素名) と、
- matplotlib に渡す点列・ラティスの箱・名前が、`OpticsPlot[…, DisplayFunction->Identity]` が返す Graphics (Canvas.n の経路。gs と同じ式) の
  Line・Rectangle・Text と一致すること、
- 保存 (png・pdf・svg、環境変数 PYSAD_DRAW_FILE、Draw$Option の File)、matplotlib が無いときの動作。
matplotlib は Agg (窓を出さない)。Graphics の照合は PYSAD_SAD_PACKAGES が要る (なければ skip)。
"""
import io
import os

import pytest

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')

from pysad import draw_mpl as dm                      # noqa: E402
from pysad.ffscmd import run_sad_string                # noqa: E402

PK = os.environ.get('PYSAD_SAD_PACKAGES')
needs_pk = pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')

FODO = '''MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.2 K2=0);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 6*CELL1);
FFS;
'''
ASTRID2 = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')


def _run(body, display_off=True, env=None, exit_=True):
    old = {}
    env = dict(env or {})
    if display_off:
        env.setdefault('PYSAD_NO_DISPLAY', '1')
    for k, v in env.items():
        old[k] = os.environ.get(k)
        os.environ[k] = v
    try:
        buf = io.StringIO()
        ffs, I = run_sad_string(body + 'Exit;\n', out=buf)
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
    return ffs, I, buf.getvalue()


@pytest.fixture(scope='module')
def fodo():
    ffs, I, _ = _run(FODO + 'USE RING; CELL; CALC;\n')
    return ffs, I


# --------------------------------------------------------------------------------------------- 構文解析
def test_separate_ampsnd_like_canvasdraw():
    assert dm.separate_ampsnd('BX') == ('BX', '')
    assert dm.separate_ampsnd('BX&BY') == ('BX', '&BY')
    assert dm.separate_ampsnd('&BY') == ('&', 'BY')
    assert dm.separate_ampsnd('BX;') == ('BX', ';')
    assert dm.separate_ampsnd(';') == (';', '')
    assert dm.separate_ampsnd('') == ('', '')


def test_is_label():
    for n in ('BX', 'by', 'EX', 'NX', 'PEX', 'PEPY', 'AZ', 'DETR', 'ZX', 'RBX', 'DBX', 'BMAGX', 'BZS', 'SIGX', 'SIZEXX'):
        assert dm.is_label(n), n
    for n in ('LAT', 'LATTICE', 'QF', 'S', 'GAMMA', 'B*', '&', ''):
        assert not dm.is_label(n), n


@pytest.mark.parametrize('stmt,panels,r1,r2,lat,names,left', [
    ('DRAW BX BY & EX LAT', [['BX', 'BY'], ['EX']], 1, 'N', True, None, ''),
    ('DRAW BX&BY&EX LATTICE', [['BX'], ['BY'], ['EX']], 1, 'N', True, None, ''),
    ('DRAW 3 10 BX BY', [['BX', 'BY']], 3, 10, False, None, ''),
    ('DRAW QF.2 QD.3 BX & EX', [['BX'], ['EX']], None, None, False, None, ''),
    ('DRAW 5 BX', [['BX']], 5, 'N', False, None, ''),
    ('DRAW BX BY & EX B*', [['BX', 'BY'], ['EX']], 1, 'N', True, 'B*', ''),
    ('DRAW BX CALC', [['BX']], 1, 'N', False, None, 'CALC'),
    ('DRAW bx by', [['BX', 'BY']], 1, 'N', False, None, ''),
    ('DRAW 1.5 BX', [['BX']], 1, 'N', False, None, ''),        # 実数は整数でも要素名でもない: r1=1 の既定で、その語は関数として読む → 1.5 は関数でない
])
def test_parse_draw(fodo, stmt, panels, r1, r2, lat, names, left):
    ffs, _ = fodo
    words = stmt.split()
    spec = dm.parse_draw(ffs, words, 1)
    n = ffs.nlat
    if stmt.startswith('DRAW 1.5'):
        # 1.5 は関数でないので関数の並びは空になり、LAT でも要素名でもないので語は FFS に返る
        assert spec.panels == [[]] and spec.leftover == '1.5'
        return
    assert spec.panels == panels
    if r1 is not None:
        assert spec.r1 == r1
    else:
        assert spec.r1 == ffs.positions('QF.2')[0] and spec.r2 == ffs.positions('QD.3')[0]
    if r2 is not None:
        assert spec.r2 == (n if r2 == 'N' else r2)
    assert spec.lattice is lat and spec.names == names and spec.leftover == left
    assert spec.consumed == len(words) - 1 - (1 if left else 0)


def test_parse_draw_bad_r1_falls_back_to_default(fodo):
    ffs, _ = fodo
    spec = dm.parse_draw(ffs, ['DRAW', 'XYZ', 'BX'], 1)
    assert spec.r1 == 1 and spec.panels == [[]] and spec.leftover == 'XYZ'      # XYZ は要素でも関数でもない


def test_unit_and_label_tables():
    assert dm.twiss_unit('BX') == 'SqrtMeter' and dm.twiss_unit('EX') == 'Meter' and dm.twiss_unit('NX') == 'Pi2'
    assert dm.twiss_unit('EPX') == 'Ratio' and dm.twiss_unit('R3') == 'InvMeter' and dm.twiss_unit('AX') == '1'
    assert dm._optics_unit('Meter', 0.02) == 'Milimeter' and dm._optics_unit('Meter', 2.0) == 'Meter'
    assert dm.twiss_label('BX') == r'$\sqrt{\beta_{x}}$'


# --------------------------------------------------------------------------------------------- Graphics との照合
GRAPHICS_CASES = [
    ([['BX', 'BY'], ['EX']], (None, None), True, '*'),
    ([['BX', 'BY'], ['EX'], ['NX', 'NY']], (3, 20), True, 'B*'),
    ([['BX'], ['EY', 'EX']], ('QF.2', 'QD.4'), False, '*'),
    ([['BX', 'BY']], (30, 4), True, 'Q*'),                  # 周回をまたぐ領域 (pr2 <= pr1)
    ([['AX', 'AY']], (5, 12), False, '*'),
]


def _opt_text(region, lat, names):
    r1, r2 = region
    f = lambda r: '"%s"' % r if isinstance(r, str) else str(r)
    s = ''
    if r1 is not None:
        s += ',Region->{%s,%s}' % (f(r1), f(r2))
    s += ',Lattice->%s' % ('True' if lat else 'False')
    if lat:
        s += ',Names->"%s"' % names
    return s


def _check_against_graphics(prefix, panels, region, lat, names):
    """pysad の OpticsPlot (Canvas.n の経路、Identity) が返す Graphics の点列・箱・名前と、matplotlib に渡すデータの一致"""
    lst = '{' + ','.join('{' + ','.join('"%s"' % f for f in p) + '}' for p in panels) + '}'
    ffs, I, _ = _run(prefix + 'g=OpticsPlot[%s%s,DisplayFunction->Identity];\n' % (lst, _opt_text(region, lat, names)))
    g = I.globals['g']
    lines = dm.graphics_lines(g)
    data = dm.build_plot_data(ffs, panels, region, lat, names)
    ser = [s for p in data.panels for s in p.series]
    assert len(lines) == len(ser) == sum(len(p) for p in panels)
    for a, (xs, ys) in zip(lines, ser):
        assert len(a) == len(xs) > 10
        for (gx, gy), x, y in zip(a, xs, ys):
            assert gx == pytest.approx(x, abs=1e-12, rel=1e-12)
            assert gy == pytest.approx(y, abs=1e-12, rel=1e-12)
    if lat:
        rects = dm.graphics_rectangles(g)
        assert len(rects) == len(data.lattice.rects) > 3
        for a, b in zip(rects, data.lattice.rects):
            assert a == pytest.approx(b[:4], abs=1e-12)
        texts = dm.graphics_texts(g)
        assert [t[0] for t in texts] == [n for _, n in data.lattice.names]
        for (_, gx), (x, _) in zip(texts, data.lattice.names):
            assert gx == pytest.approx(x, abs=1e-12)
    return data


@needs_pk
@pytest.mark.parametrize('panels,region,lat,names', GRAPHICS_CASES)
def test_points_match_opticsplot_graphics_fodo(panels, region, lat, names):
    _check_against_graphics(FODO + 'USE RING; CELL; CALC;\n', panels, region, lat, names)


@needs_pk
@pytest.mark.skipif(not os.path.exists(ASTRID2), reason='ASTRID2 のラティスが無い')
def test_points_match_opticsplot_graphics_astrid2():
    pre = open(ASTRID2, encoding='utf-8').read() + '\nFFS;\nUSE RING; CALC;\n'
    d = _check_against_graphics(pre, [['BX', 'BY'], ['EX']], (None, None), True, '*')
    assert len(d.lattice.rects) > 50
    _check_against_graphics(pre, [['BX', 'BY'], ['EX', 'PEX']], (10, 150), True, 'Q*')    # 領域・名前のパターンつき


@needs_pk
def test_unit_scaling_matches_graphics(fodo):
    """EX は mm に直す (Optics$Unit)。Graphics の y と一致 (上の照合) し、軸のラベルに単位が付く"""
    ffs, I = fodo
    d = dm.build_plot_data(ffs, [['EX'], ['BX']], (None, None), False)
    assert d.panels[0].unit == 'Milimeter' and d.panels[0].unit_label == 'mm'
    assert d.panels[0].ylabel.endswith('(mm)')
    assert d.panels[1].ylim[0] == 0.0 and d.panels[1].ylim[1] > max(d.panels[1].series[0][1])     # 正の関数は 0 から


# --------------------------------------------------------------------------------------------- 描画・保存
def test_draw_returns_figure_and_saves(fodo, tmp_path):
    ffs, _ = fodo
    out = tmp_path / 'x.png'
    fig = dm.draw(ffs, 'BX BY & EX LAT', save=str(out))
    assert out.exists() and out.stat().st_size > 5000
    axes = fig._pysad_draw['axes']
    assert len(axes) == 3                                  # 2 枠 + ラティス
    assert axes[0].get_xlim() == axes[2].get_xlim()        # s 軸を共有
    assert len(axes[0].get_lines()) == 2 and len(axes[1].get_lines()) == 1
    import matplotlib.pyplot as plt
    plt.close(fig)
    for ext in ('pdf', 'svg'):
        p = tmp_path / ('y.' + ext)
        plt.close(dm.draw(ffs, 'BX', save=str(p)))
        assert p.exists() and p.stat().st_size > 1000


def test_draw_region_and_title(fodo):
    ffs, _ = fodo
    import matplotlib.pyplot as plt
    fig = dm.draw(ffs, 'BX', region=(3, 10), title='T1', lattice=True, names='B*')
    assert fig._suptitle.get_text() == 'T1'
    x0, x1 = fig.axes[0].get_xlim()
    s3 = float(ffs.line_value('S', 3))
    s10 = float(ffs.line_value('S', 10))
    assert x0 == pytest.approx(s3 - 0.001) and x1 == pytest.approx(s10 + 0.001)
    plt.close(fig)


def test_draw_info_text(fodo):
    ffs, _ = fodo
    t = dm.element_info(ffs, 1.0)
    assert t.startswith('Entrance of ') and 'beta_x,y' in t and 'nu_x,y' in t
    fig = dm.draw(ffs, 'BX', tunes=True)
    assert any('nu' in x.get_text() for x in fig.texts)
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_draw_without_calc_computes_optics(tmp_path):
    """CALC していなくても (pysad の Twiss[] と同じく) 先に光学を計算する。gs の OpticsPlot は計算していない光学の値 (0) を描く"""
    ffs, I, _ = _run(FODO + 'USE RING; CELL;\n')
    d = dm.build_plot_data(ffs, [['BX']], (None, None), False)
    assert max(d.panels[0].series[0][1]) > 1.0


def test_draw_command_saves_file_and_returns_rest(tmp_path):
    out = tmp_path / 'cmd.png'
    ffs, I, txt = _run(FODO + 'USE RING; CELL; CALC;\nDRAW BX BY & EX LAT;\nDRAW 3 12 BX;\n',
                       env={'PYSAD_DRAW_FILE': str(out)})
    assert out.exists()
    assert (tmp_path / 'cmd_2.png').exists()               # 2 回目は番号つき (上書きしない)
    assert '保存しました' in txt and '未対応' not in txt


def test_draw_command_draw_option_file(tmp_path):
    out = (tmp_path / 'opt.svg').as_posix()
    ffs, I, txt = _run(FODO + 'Draw$Option={File->"%s",Region->{2,9},InfoLabel->False};\nUSE RING; CELL; CALC;\nDRAW BX LAT;\n' % out)
    assert os.path.exists(out)


def test_draw_command_identity_does_nothing(tmp_path):
    """K64: $DisplayFunction が CanvasDrawer でなければ DRAW は何もしない (保存先も無いとき)"""
    ffs, I, txt = _run(FODO + 'USE RING; CELL; CALC;\nDRAW BX BY;\n', env={'PYSAD_NO_DISPLAY': '1'})
    assert '未対応' not in txt and 'DRAW' not in txt


def test_draw_command_passes_leftover_to_ffs(tmp_path):
    """DRAW が関数でない語に当たったら、その語以降は FFS の次のコマンド (CALC)"""
    out = tmp_path / 'lo.png'
    ffs, I, txt = _run(FODO + 'USE RING; CELL;\nDRAW BX CALC;\n', env={'PYSAD_DRAW_FILE': str(out)})
    assert out.exists() and 'Matched' in txt              # CALC が実行された


def test_draw_without_matplotlib_warns(monkeypatch, tmp_path):
    monkeypatch.setattr(dm, '_import_mpl', lambda: None)
    ffs, I, txt = _run(FODO + 'USE RING; CELL; CALC;\nDRAW BX;\n', env={'PYSAD_DRAW_FILE': str(tmp_path / 'n.png')})
    assert 'matplotlib がありません' in txt and not (tmp_path / 'n.png').exists()


def test_draw_no_function_message(tmp_path):
    ffs, I, txt = _run(FODO + 'USE RING; CELL; CALC;\nDRAW LAT;\n', env={'PYSAD_DRAW_FILE': str(tmp_path / 'e.png')})
    assert '描く関数がありません' in txt


def test_draw_different_units_dropped(fodo):
    ffs, _ = fodo
    msgs = []
    d = dm.build_plot_data(ffs, [['BX', 'EX']], (None, None), False, warn=msgs.append)
    assert d.panels[0].funcs == ['BX'] and d.panels[0].dropped == ['EX'] and msgs


# --------------------------------------------------------------------------------------------- gs (K64) との照合
# gs は表示環境が無いと DRAW を描かないので、OpticsPlot[…, DisplayFunction->Identity] が返す Graphics の Line・Rectangle・Text を
# 書き出した参照 (scripts/v25_draw_gs_*.gs.txt。台本は同名の .sad、gs は 2026-11 の WSL ビルド) と、matplotlib に渡すデータを比べる。
_SCR = os.path.join(os.path.dirname(__file__), 'scripts')


def _read_gs(name):
    import re
    lines, rects, texts = [], [], []
    num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
    for ln in open(os.path.join(_SCR, name), encoding='utf-8').read().splitlines():
        if ln.startswith('LINE '):
            v = [float(x) for x in re.findall(num, ln[ln.index('{'):])]
            lines.append((v[0::2], v[1::2]))
        elif ln.startswith('RECT '):
            v = [float(x) for x in re.findall(num, ln.split(None, 2)[2])]
            rects.append(tuple(v))
        elif ln.startswith('TEXT '):
            p = ln.split()
            texts.append((p[2], float(p[3])))
    return lines, rects, texts


def _compare_gs(prefix, ref, panels, region, names, rtol=1e-10):
    ffs, I, _ = _run(prefix)
    data = dm.build_plot_data(ffs, panels, region, True, names)
    lines, rects, texts = _read_gs(ref)
    ser = [s for p in data.panels for s in p.series]
    assert len(lines) == len(ser)
    for (gx, gy), (xs, ys) in zip(lines, ser):
        assert len(gx) == len(xs) > 10
        for a, b in zip(gx, xs):
            assert a == pytest.approx(b, rel=1e-12, abs=1e-12)
        for a, b in zip(gy, ys):
            assert a == pytest.approx(b, rel=rtol, abs=1e-9 * max(abs(v) for v in gy))     # gs の出力は 15 桁
    assert len(rects) == len(data.lattice.rects)
    for a, b in zip(rects, data.lattice.rects):
        assert a == pytest.approx(b[:4], abs=1e-10)
    assert [t[0] for t in texts] == [n for _, n in data.lattice.names]
    for (_, gx), (x, _) in zip(texts, data.lattice.names):
        assert gx == pytest.approx(x, rel=1e-12)


def test_matches_gs_graphics_fodo():
    _compare_gs(FODO + 'USE RING; CELL; CALC;\n', 'v25_draw_gs_fodo.gs.txt', [['BX', 'BY'], ['EX']], (3, 20), 'B*')


@pytest.mark.skipif(not os.path.exists(ASTRID2), reason='ASTRID2 のラティスが無い')
def test_matches_gs_graphics_astrid2():
    pre = open(ASTRID2, encoding='utf-8').read() + '\nFFS;\nUSE RING; CALC;\n'
    _compare_gs(pre, 'v25_draw_gs_astrid2.gs.txt', [['BX', 'BY'], ['EX']], (None, None), '*', rtol=1e-7)
