# -*- coding: utf-8 -*-
"""SADScript の `Graphics[…]` (Plot・ListPlot・OpticsPlot・GeometryPlot ほか) を matplotlib に描く `$DisplayFunction` (v2.4.1)。

SAD の `$DisplayFunction` は、`Show1$` (ListPlot.n) が `df[opt, list, init, w, rect, last, ini]` の形で呼ぶ関数で、K64 では `CanvasDrawer`
(Canvas.n、Tk のキャンバス) が実体。pysad はこの関数の名前 `MatplotlibDrawer` を組み込み関数として足し、同じ引数を受けて
Graphics の要素 (Line・Point・Rectangle・Polygon・Circle・Text。`Rectangle[…, Graphics[…]]` の入れ子の枠を含む) を matplotlib の
Figure に描く。評価器 (sadscript.py) は変えない。Canvas.n の `CanvasDrawer` と同じく Show1$ の流れで枠ごとに呼ばれるので、
PlotRange・DataRange・Scale・GridLines・FrameTicks などは Show1$ が決めた値 (opt) を使う。

切替 (`pkgload.py`):
    環境変数 PYSAD_DISPLAY=mpl (または matplotlib)   $DisplayFunction を MatplotlibDrawer にする
    SADScript で `$DisplayFunction=MatplotlibDrawer;` または `Plot[…, DisplayFunction->MatplotlibDrawer]`
    PYSAD_DISPLAY=tk (CanvasDrawer)・none (Identity) も受ける。未設定は従来どおり (Windows か DISPLAY があれば CanvasDrawer、
    PYSAD_NO_DISPLAY=1 で Identity)

出力:
    環境変数 PYSAD_DISPLAY_FILE (拡張子 png・pdf・svg)。同じ名前の 2 回目以降は `名前_2.png` のように番号を付ける
    (`{n}` があれば番号に置き換える。DRAW の PYSAD_DRAW_FILE と同じ流儀、番号の数え方も draw_mpl と共通)。
    SAD 側は大域変数 `MatplotlibDrawer$Option={File->"x.png", Show->True, DPI->120}` (File・Show・DPI)。
    保存先がなければ、対話の backend のときだけ plt.show (Agg など非対話では何も出さない)。

Python から:
    from pysad import mpl_display
    fig = mpl_display.render(I, g, save="x.png")       # I: インタプリタ、g: Graphics[…] の式 (DisplayFunction->Identity で得たもの)
    mpl_display.last_figure()                           # 最後に描いた Figure

matplotlib は遅延 import。無ければ警告を出すだけで何もしない (評価は続く)。未対応の要素は警告して飛ばす (docs/42 §7c)。
"""
from __future__ import annotations

import math
import os
import re
from typing import Dict, List, Optional, Tuple

from .sadscript import Expr, Null, Rule, Sym, S, is_num

__all__ = ['install', 'render', 'last_figure', 'sad_color', 'sad_text', 'display_mode_from_env', 'make_window',
           'frame_range', 'MatplotlibDrawer']

# ---------------------------------------------------------------------------------------------------------------------
# 1. Canvas.n の幾何 (キャンバス 13 x 10 インチ。窓の原点 {2, 1.8}、幅 9.6、縦横比 GoldenRatio)
# ---------------------------------------------------------------------------------------------------------------------
GOLDEN = (1.0 + math.sqrt(5.0)) / 2.0
CANVAS = (13.0, 10.0)                    # $DefaultCanvasSize
ORIGIN = (2.0, 1.8)                      # $DefaultWindowOrigin
RAW = (9.6, 9.6 / GOLDEN)                # $DefaultWindowWidth, 高さ = 幅/GoldenRatio
MARGIN = ((0.02, 0.98), (0.02, 0.98))
FIG_SCALE = 0.65                         # キャンバスの 1 インチ → 図の 0.65 インチ
PT_PER_INCH = 72.0 * FIG_SCALE           # 破線の長さ (インチ) → 点
LW_PT = 0.8                              # Tk の線の太さ 1 px → 点
FONT_TEXT, FONT_TICK, FONT_LABEL, FONT_TITLE, FONT_LEGEND = 8.5, 8.0, 9.0, 11.0, 8.0
EPS = 1e-300


def _clamp(v: float, lo: float, hi: float) -> float:
    return min(max(v, lo), hi)


def make_window(plot_region, rect, ar) -> Tuple[float, float, float, float]:
    """Canvas.n の MakeWindow。plot_region {{px1,px2},{py1,py2}}・rect (Automatic か {{x1,y1},{x2,y2}}、PlotRegion の窓の割合)・
    ar (縦横比 = 幅/高さ。実数でなければ Automatic)。返り値は窓 (x1, x2, y1, y2) (キャンバスのインチ、y は下から)"""
    (px1, px2), (py1, py2) = plot_region
    wx1, wx2 = ORIGIN[0], ORIGIN[0] + RAW[0]
    wy1, wy2 = ORIGIN[1], ORIGIN[1] + RAW[1]
    x1 = _clamp(px1 * (wx2 - wx1) + wx1, CANVAS[0] * MARGIN[0][0], CANVAS[0] * MARGIN[0][1])
    x2 = _clamp(px2 * (wx2 - wx1) + wx1, CANVAS[0] * MARGIN[0][0], CANVAS[0] * MARGIN[0][1])
    y1 = _clamp(py1 * (wy2 - wy1) + wy1, CANVAS[1] * MARGIN[1][0], CANVAS[1] * MARGIN[1][1])
    y2 = _clamp(py2 * (wy2 - wy1) + wy1, CANVAS[1] * MARGIN[1][0], CANVAS[1] * MARGIN[1][1])
    real_ar = is_num(ar)
    if isinstance(rect, list) and len(rect) == 2:
        (rx1, ry1), (rx2, ry2) = rect
        wx1, wx2 = rx1 * (x2 - x1) + x1, rx2 * (x2 - x1) + x1
        wy1, wy2 = ry1 * (y2 - y1) + y1, ry2 * (y2 - y1) + y1
        dar = (wx2 - wx1) / (wy2 - wy1) if wy2 != wy1 else GOLDEN
        xc, yc = (wx1 + wx2) / 2, (wy1 + wy2) / 2
        if real_ar and ar > 0:
            if ar <= dar:
                wx1, wx2 = xc - (wy2 - wy1) * ar / 2, xc + (wy2 - wy1) * ar / 2
            else:
                wy1, wy2 = yc - (wx2 - wx1) / ar / 2, yc + (wx2 - wx1) / ar / 2
        return wx1, wx2, wy1, wy2
    wx1, wx2, wy1, wy2 = x1, x2, y1, y2
    ar0 = float(ar) if real_ar and ar > 0 else GOLDEN
    xc, yc = (wx1 + wx2) / 2, (wy1 + wy2) / 2
    if ar0 <= (wx2 - wx1) / (wy2 - wy1):
        wx1, wx2 = xc - (wy2 - wy1) * ar0 / 2, xc + (wy2 - wy1) * ar0 / 2
    else:
        wy1, wy2 = yc - (wx2 - wx1) / ar0 / 2, yc + (wx2 - wx1) / ar0 / 2
    return wx1, wx2, wy1, wy2


def _fr1(scale: str, r, x1: float, x2: float) -> float:
    """FrameRange1 (r が Automatic のときは x1 から 7% の余白、実数ならその値)。Log は log10"""
    if scale == 'Log':
        if r is None:
            a, b = math.log10(max(x1, EPS)), math.log10(max(x2, EPS))
            return a + 0.07 * (a - b)
        return math.log10(max(EPS, r))
    return x1 + 0.07 * (x1 - x2) if r is None else float(r)


def frame_range(scale: str, r, b) -> Tuple[float, float]:
    """FrameRange + FrameRange2。r = (r1, r2) (None は Automatic)、b = DataRange の (b1, b2)。返り値は軸の範囲 (Log は log10)"""
    r1, r2 = r
    b1, b2 = b
    if r1 is None and r2 is None:
        lo, hi = _fr1(scale, None, b1, b2), _fr1(scale, None, b2, b1)
    elif r1 is None:
        lo, hi = _fr1(scale, None, min(b1, r2), r2), _fr1(scale, r2, r2, b1)
    elif r2 is None:
        lo, hi = _fr1(scale, r1, r1, b2), _fr1(scale, None, max(b2, r1), r1)
    else:
        lo, hi = _fr1(scale, r1, r1, r2), _fr1(scale, r2, r2, r1)
    if lo == 0 and hi == 0:                                   # FrameRange2
        return -0.001, 0.001
    if lo == hi:
        return min(lo * 0.999, lo * 1.001), max(lo * 0.999, lo * 1.001)
    if hi - lo < 1e-100:
        m = (lo + hi) / 2
        return m * 0.999, m * 1.001
    return lo, hi


# ---------------------------------------------------------------------------------------------------------------------
# 2. 色・文字
# ---------------------------------------------------------------------------------------------------------------------
_X11_SHADE = {1: 1.0, 2: 0.932, 3: 0.804, 4: 0.545}      # X11 の色名の末尾の数字 (1〜4) の明るさ
_SYSTEM_COLORS = {'systembuttontext': 'black', 'systemwindowtext': 'black', 'systemwindow': 'white',
                  'systembuttonface': 'lightgray', 'systemhighlight': 'royalblue', 'systemhighlighttext': 'white'}
PALETTE = ['royal blue', 'tomato3', 'forest green', 'dark goldenrod', 'dodger blue', 'orange3', 'blue violet',
           'lime green', 'hot pink', 'yellow4']          # Plot$DefaultColor (ListPlot.n)


def sad_color(name) -> Optional[Tuple[float, float, float]]:
    """Tk (X11) の色名・#rrggbb を matplotlib の RGB にする。空文字・文字列でないもの・未知の名前は None。
    `dodger blue` のような空白つき、`tomato3` のような末尾の数字 (1〜4 = 明るさ 100・93.2・80.4・54.5%)、`gray40`、
    `#rgb`・`#rrggbb`・`#rrrgggbbb` を受ける"""
    if not isinstance(name, str):
        return None
    s = name.strip()
    if not s:
        return None
    try:
        import matplotlib.colors as mc
    except Exception:
        return None
    if s.startswith('#'):
        h = s[1:]
        if len(h) in (3, 6, 9, 12) and re.fullmatch(r'[0-9a-fA-F]+', h):
            n = len(h) // 3
            return tuple(int(h[k * n:(k + 1) * n], 16) / (16 ** n - 1) for k in range(3))
        return None
    low = re.sub(r'\s+', '', s).lower()
    if low in _SYSTEM_COLORS:
        low = _SYSTEM_COLORS[low]
    m = re.fullmatch(r'gr[ae]y(\d{1,3})', low)
    if m:
        v = min(int(m.group(1)), 100) / 100.0
        return (v, v, v)
    try:
        return tuple(mc.to_rgb(low))
    except ValueError:
        pass
    m = re.fullmatch(r'([a-z]+?)([1-4])', low)
    if m:
        try:
            r, g, b = mc.to_rgb(m.group(1))
        except ValueError:
            return None
        f = _X11_SHADE[int(m.group(2))]
        return (r * f, g * f, b * f)
    return None


# SAD の文字の書式 (` で始まる 2 文字の制御): `f = Symbol フォント (ギリシャ文字)、`n = 通常に戻る、`d = 下付き、`u = 上付き、
# `o…`O = 重ね書き (ここでは並べて書く)。そのほかの制御は読み飛ばす
_SYMBOL = dict(zip('abcdefghijklmnopqrstuvwxyz', 'αβχδεφγηιϕκλμνοπθρστυϖωξψζ'))
_SYMBOL.update(zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'ΑΒΧΔΕΦΓΗΙϑΚΛΜΝΟΠΘΡΣΤΥςΩΞΨΖ'))
_SYMBOL.update({'\xd6': '√', '\xb4': '×', '\xb0': '°', '\xb1': '±', '\xa5': '∞', '\xb6': '∂', '\xd1': '∇', '\xb7': '·'})
_MATH_ESC = {'_': r'\_', '^': r'\^{}', '{': r'\{', '}': r'\}', '\\': r'\backslash ', '#': r'\#', '%': r'\%', '&': r'\&', '$': r'\$'}


def _mathtext_ok(txt: str) -> bool:
    try:
        from matplotlib.mathtext import MathTextParser
        pr = MathTextParser('agg')
        for ln in txt.split(chr(10)):                                # 複数行は 1 行ずつ (matplotlib も行ごとに mathtext にする)
            if '$' in ln:
                pr.parse(ln)
        return True
    except Exception:                                               # noqa: BLE001
        return False


def sad_text(s) -> str:
    r"""SAD の文字の書式を matplotlib の mathtext にする (`fh`n`dx`n → η$_{x}$)。書式がなければそのまま。
    バックスラッシュ + @ で始まる LaTeX の文字 (FitPlot の題など) は、先頭のそれと centering の指定を取り、二重のバックスラッシュを改行にして
    `$…$` を mathtext に任せる。
    mathtext に読めない文字は `$` を打ち消した平文にする (描画で落ちない)"""
    if not isinstance(s, str):
        s = '' if s is None else str(s)
    out = _sad_markup(s)
    if '$' in out and not _mathtext_ok(out):
        out = out.replace(BS + '$', '$').replace('$', BS + '$')
    return out


BS = chr(92)


def _sad_markup(s: str) -> str:
    if s.startswith(BS + '@'):
        return s[2:].replace(BS + 'centering', '').replace(BS * 2, chr(10)).strip()
    if '`' not in s:
        return s
    pieces: List[Tuple[str, str]] = []
    mode, sym, i = 'n', False, 0
    while i < len(s):
        c = s[i]
        if c == '`' and i + 1 < len(s):
            k = s[i + 1]
            i += 2
            if k == 'n':
                mode, sym = 'n', False
            elif k == 'f':
                sym = True
            elif k == 'd':
                mode = 'sub'
            elif k == 'u':
                mode = 'sup'
            continue
        pieces.append((mode, _SYMBOL.get(c, c) if sym else c))
        i += 1
    out: List[str] = []
    k = 0
    while k < len(pieces):
        md = pieces[k][0]
        j = k
        while j < len(pieces) and pieces[j][0] == md:
            j += 1
        txt = ''.join(p[1] for p in pieces[k:j])
        if md == 'n':
            out.append(txt.replace('$', r'\$'))
        else:
            out.append('$' + ('_' if md == 'sub' else '^') + '{' + ''.join(_MATH_ESC.get(ch, ch) for ch in txt) + '}$')
        k = j
    return ''.join(out)


# ---------------------------------------------------------------------------------------------------------------------
# 3. 値の取り出し
# ---------------------------------------------------------------------------------------------------------------------
def _name(v) -> Optional[str]:
    return v.name if isinstance(v, Sym) else None


def _truth(v, default=False) -> bool:
    """SAD の真偽 (True は 1.0、False は 0.0。記号 True・False・None も受ける)"""
    if is_num(v):
        return v != 0
    n = _name(v)
    if n in ('True', 'Yes', 'All'):
        return True
    if n in ('False', 'None', 'No', 'Null'):
        return False
    return default


def _head(e) -> str:
    h = getattr(e, 'head', None)
    return _name(h) or str(h)


def _is_auto(v) -> bool:
    return _name(v) == 'Automatic'


def _opts(items, ev=None) -> Dict[str, object]:
    """Rule のリスト → {名前: 値}。先に書いたものが優先 (Override)。遅延規則 (:>) は ev があれば評価する"""
    out: Dict[str, object] = {}
    for r in items or []:
        if isinstance(r, Rule):
            k = _name(r.lhs) or (r.lhs if isinstance(r.lhs, str) else str(r.lhs))
            if k in out:
                continue
            v = r.rhs
            if r.delayed and ev is not None:
                try:
                    v = ev(v)
                except Exception:                                   # noqa: BLE001
                    pass
            out[k] = v
        elif isinstance(r, list):
            out.update({k: v for k, v in _opts(r, ev).items() if k not in out})
    return out


def _flat(items):
    for x in items:
        if isinstance(x, list):
            yield from _flat(x)
        else:
            yield x


def _first(v, default=None):
    """リストなら先頭、スカラーならそのもの"""
    if isinstance(v, list):
        return v[0] if v else default
    return v


# ---------------------------------------------------------------------------------------------------------------------
# 4. 描画
# ---------------------------------------------------------------------------------------------------------------------
_MARKERS = {'1O': 'o', 'CI': 'o', 'CL': 'o', '6O': '<', 'LT': '<', '7O': '>', 'RT': '>', '8O': 'v', 'DT': 'v',
            '9O': '^', 'UT': '^', 'BX': 's', 'SQ': 's', 'RH': 'D', 'DM': 'D', 'DI': 'D', 'PL': 'P', 'TI': 'X', 'ML': 'X',
            'MU': 'X'}
_ANCHOR = {'n': ('center', 'top'), 's': ('center', 'bottom'), 'e': ('right', 'center'), 'w': ('left', 'center'),
           'c': ('center', 'center'), 'ne': ('right', 'top'), 'nw': ('left', 'top'), 'se': ('right', 'bottom'),
           'sw': ('left', 'bottom'), 'left': ('left', 'center'), 'right': ('right', 'center'), 'center': ('center', 'center'),
           'top': ('center', 'top')}


class _Unsupported(Exception):
    pass


class _Frame:
    """1 つの枠 (Canvas の Canvas$Frame)。ax は範囲が正しければ作る (OpticsPlot の外側の枠のように範囲が空なら None で、入れ子の枠だけ持つ)"""

    def __init__(self, idx: int, op: Dict[str, object]):
        self.idx = idx
        self.op = op
        self.ax = None
        self.xlim: Tuple[float, float] = (0.0, 1.0)
        self.ylim: Tuple[float, float] = (0.0, 1.0)
        self.scales = ('Linear', 'Linear')
        self.window: Optional[Tuple[float, float, float, float]] = None
        self.legend: List[Tuple[str, dict]] = []
        self.texts: List[Tuple[object, str]] = []
        self.auto = 0                                              # 色が評価されていない線に順に割り当てる色の番号
        self.rmo = ''
        self.zc = 0
        self.rec: Dict[str, list] = {'lines': [], 'points': [], 'rects': [], 'polys': [], 'ovals': [], 'texts': [],
                                     'texts_removed': [], 'errorbars': []}
        self.font_scale = 1.0
        self.legend_scale = 1.0
        self.legend_bg = None

    def z(self) -> float:
        self.zc += 1
        return 3.0 + 0.0001 * self.zc


class _Session:
    """1 回の Show (ini=True) = 1 つの Figure"""

    def __init__(self, drv: '_Driver', op: Dict[str, object]):
        import matplotlib.pyplot as plt
        self.drv = drv
        self.fig = plt.figure(figsize=(CANVAS[0] * FIG_SCALE, CANVAS[1] * FIG_SCALE))
        self.fig.patch.set_facecolor('white')
        pr = op.get('PlotRegion')
        self.plot_region = pr if (isinstance(pr, list) and len(pr) == 2 and all(isinstance(a, list) and len(a) == 2 for a in pr)) \
            else [[0.0, 1.0], [0.0, 1.0]]
        self.plot_region = [[float(a[0]), float(a[1])] for a in self.plot_region]
        self.frames: List[_Frame] = []
        self.cur: Optional[_Frame] = None
        self.title: Optional[str] = None
        self.title_scale = 1.0
        self.warned = set()
        self.warnings: List[str] = []
        self.parts: List[int] = []                 # 呼び出しごとの (opt, list) の指紋。直前の Show と同じ内容か調べる


class _Driver:
    """インタプリタ I に 1 つ。MatplotlibDrawer の本体"""

    def __init__(self, I):
        self.I = I
        self.sess: Optional[_Session] = None
        self.override: Dict[str, object] = {}          # render() が一時的に与える File・Show・DPI
        self.last_fig = None
        self.last_fp = None                            # 直前に出力した Show の指紋 (ListPlot.n 自身が Show を 2 回呼ぶので、続く同じ内容は 1 つにする)
        self.warned_nompl = False

    # ---------------------------------------------------------------- 入口
    def write(self, msg: str):
        try:
            self.I.out.write(msg + '\n')
        except Exception:                                           # noqa: BLE001
            pass

    def warn(self, key: str, msg: str):
        s = self.sess
        text = f'[pysad] MatplotlibDrawer: {msg}'
        if s is None:
            self.write(text)
            return
        if key in s.warned:
            return
        s.warned.add(key)
        s.warnings.append(msg)
        self.write(text)

    def call(self, args):
        from . import draw_mpl
        if draw_mpl._import_mpl() is None:
            if not self.warned_nompl:
                self.warned_nompl = True
                self.write('[pysad] MatplotlibDrawer: matplotlib がありません。何も描きません (pip install matplotlib)')
            return Null
        try:
            return self._call(args)
        except Exception as e:                                      # noqa: BLE001
            self.write(f'[pysad] MatplotlibDrawer: 描けません: {type(e).__name__}: {e}')
            self._abort()
            return Null

    def _abort(self):
        s, self.sess = self.sess, None
        if s is not None:
            try:
                import matplotlib.pyplot as plt
                plt.close(s.fig)
            except Exception:                                       # noqa: BLE001
                pass

    def _call(self, args):
        if len(args) != 7:
            raise ValueError('MatplotlibDrawer[opt, list, init, w, rect, last, ini] の引数は 7 個')
        opt, lst, init, w, rect, last, ini = args
        w, ini, last = _truth(w), _truth(ini), _truth(last)
        op = _opts(opt if isinstance(opt, list) else [], self.I.eval)
        fp = hash((repr(opt), repr(lst), repr(rect)))
        if w:
            if ini:
                self._abort()                                       # 前の Show が最後まで届かなかった場合
                if not _truth(init, True) and _name(init) != 'Frame':
                    pass
                self.sess = _Session(self, op)
            elif self.sess is None:
                self.sess = _Session(self, op)
            frame = self._new_frame(op, rect)
            self.sess.cur = frame
        else:
            if self.sess is None:
                self.sess = _Session(self, op)
            if self.sess.cur is None:
                self.sess.cur = self._new_frame(op, rect)
            frame = self.sess.cur
        self.sess.parts.append(fp)
        if w:
            pl = op.get('PlotLabel')
            if isinstance(pl, Expr) and _head(pl) == 'Text' and pl.args and isinstance(pl.args[0], list) and pl.args[0]:
                pl = pl.args[0][0]
            if isinstance(pl, list) and pl and isinstance(pl[0], str):
                pl = pl[0]
            if isinstance(pl, str) and pl != '':
                self.sess.title = pl
                tfs = op.get('TitleFontScale')
                self.sess.title_scale = float(tfs) if is_num(tfs) else 1.0
        for p in _flat(lst if isinstance(lst, list) else [lst]):
            self._primitive(frame, p)
        if ini and last:
            self._finish()
        return 1.0

    # ---------------------------------------------------------------- 枠
    def _new_frame(self, op: Dict[str, object], rect) -> _Frame:
        s = self.sess
        fr = _Frame(len(s.frames) + 1, op)
        s.frames.append(fr)
        fs = op.get('FontScale')
        fr.font_scale = float(fs) if is_num(fs) else 1.0
        lfs = op.get('LegendFontScale')
        fr.legend_scale = float(lfs) if is_num(lfs) else 1.0
        fr.legend_bg = op.get('LegendBG')
        rmo = op.get('RemoveOverlap')
        fr.rmo = rmo if isinstance(rmo, str) else ''
        # Scale・PlotRange・DataRange (Canvas.n の OverrideOptions と SetupFrame)
        sca = op.get('Scale', [Sym('Linear'), Sym('Linear')])
        if isinstance(sca, Sym):
            sca = [Sym('Linear'), sca]
        scn = [(_name(a) or 'Linear') for a in (sca if isinstance(sca, list) and len(sca) == 2 else [Sym('Linear')] * 2)]
        for k, a in enumerate(scn):
            if a not in ('Linear', 'Log'):
                self.warn('scale-' + a, f'Scale->{a} は未対応 (Linear で描く)')
                scn[k] = 'Linear'
        fr.scales = (scn[0], scn[1])
        plr = op.get('PlotRange', Sym('Automatic'))
        auto2 = [[None, None], [None, None]]

        def _pair(v):
            if isinstance(v, list) and len(v) == 2:
                return [None if (_is_auto(a) or not is_num(a)) else float(a) for a in v]
            return [None, None]
        if _is_auto(plr) or plr is Default:
            plr2 = auto2
        elif isinstance(plr, list) and len(plr) == 2 and all(isinstance(a, list) for a in plr):
            plr2 = [_pair(plr[0]), _pair(plr[1])]
        elif isinstance(plr, list) and len(plr) == 2 and isinstance(plr[0], list):
            plr2 = [_pair(plr[0]), [None, None]]
        elif isinstance(plr, list) and len(plr) == 2 and _is_auto(plr[0]) and isinstance(plr[1], list):
            plr2 = [[None, None], _pair(plr[1])]
        elif isinstance(plr, list) and len(plr) == 2 and isinstance(plr[1], list):
            plr2 = [_pair(plr[0]) if isinstance(plr[0], list) else [None, None], _pair(plr[1])]
        elif isinstance(plr, list) and len(plr) == 2:
            plr2 = [[None, None], _pair(plr)]
        else:
            plr2 = auto2
        dr = op.get('DataRange', [[1e100, -1e100], [1e100, -1e100]])
        try:
            drv = [[float(dr[0][0]), float(dr[0][1])], [float(dr[1][0]), float(dr[1][1])]]
        except Exception:                                           # noqa: BLE001
            drv = [[1e100, -1e100], [1e100, -1e100]]
        rng = [frame_range(fr.scales[a], tuple(plr2[a]), tuple(drv[a])) for a in (0, 1)]
        valid = rng[0][0] < rng[0][1] and rng[1][0] < rng[1][1] and all(math.isfinite(v) for r in rng for v in r)
        if not valid:
            return fr
        lims = []
        for a in (0, 1):
            lo, hi = rng[a]
            if fr.scales[a] == 'Log':
                lo, hi = 10.0 ** lo, 10.0 ** hi
            lims.append((lo, hi))
        fr.xlim, fr.ylim = lims[0], lims[1]
        ar = op.get('AspectRatio')
        fr.window = make_window(s.plot_region, rect if isinstance(rect, list) else None, ar)
        self._make_axes(fr, op)
        return fr

    def _make_axes(self, fr: _Frame, op: Dict[str, object]):
        from matplotlib.ticker import NullFormatter
        s = self.sess
        x1, x2, y1, y2 = fr.window
        ax = s.fig.add_axes([x1 / CANVAS[0], y1 / CANVAS[1], (x2 - x1) / CANVAS[0], (y2 - y1) / CANVAS[1]])
        fr.ax = ax
        if fr.scales[0] == 'Log':
            ax.set_xscale('log')
        if fr.scales[1] == 'Log':
            ax.set_yscale('log')
        ax.set_xlim(*fr.xlim)
        ax.set_ylim(*fr.ylim)
        ax.set_axisbelow(True)
        ax.minorticks_on()
        bg = op.get('Background')
        c = sad_color(bg) if isinstance(bg, str) else None
        ax.patch.set_facecolor(c if c is not None else 'none')
        ax.patch.set_alpha(1.0 if c is not None else 0.0)
        frame_on = _truth(op.get('Frame', 1.0), True)
        fcol = sad_color(op.get('FrameColor')) or (0.0, 0.0, 0.0)
        tcol = sad_color(op.get('TickColor')) or fcol
        tfs = op.get('TickFontScale')
        tfs = float(tfs) if is_num(tfs) else 1.0
        ffs = op.get('FrameFontScale')
        ffs = [float(a) for a in ffs] if isinstance(ffs, list) and ffs and all(is_num(a) for a in ffs) \
            else ([float(ffs)] * 4 if is_num(ffs) else [1.0] * 4)
        ffs = (ffs + [1.0] * 4)[:4]
        ts = op.get('TickSize')
        tsz = float(ts) if is_num(ts) else 1.0
        fs = fr.font_scale
        if not frame_on:
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.tick_params(which='both', bottom=False, top=False, left=False, right=False,
                           labelbottom=False, labelleft=False, labeltop=False, labelright=False)
            axes_opt = op.get('Axes')
            if _truth(axes_opt, False) or (isinstance(axes_opt, list) and any(_truth(a) for a in axes_opt)):
                self.warn('axes-ext', 'Axes は Canvas.n にない拡張: Frame->False のとき原点を通る軸の線だけ描く')
                ax_x, ax_y = (axes_opt if isinstance(axes_opt, list) and len(axes_opt) == 2 else [axes_opt, axes_opt])
                if _truth(ax_x) and fr.ylim[0] < 0 < fr.ylim[1]:
                    ax.axhline(0.0, color='#606060', linewidth=0.8, zorder=1.5)
                if _truth(ax_y) and fr.xlim[0] < 0 < fr.xlim[1]:
                    ax.axvline(0.0, color='#606060', linewidth=0.8, zorder=1.5)
            return
        for sp in ax.spines.values():
            sp.set_color(fcol)
            sp.set_linewidth(0.8)
            sp.set_zorder(5)
        ft = _expand_frame_ticks(op.get('FrameTicks', Sym('Automatic')))
        kinds = [_tick_kind(v) for v in ft]
        for k, v in enumerate(ft):
            if isinstance(v, list) and len(v) >= 2 and _name(v[0]) == 'Ticks':
                self.warn('ticks-func', 'FrameTicks の {Ticks, 関数・リスト} (目盛りの文字を関数で与える形) は未対応 (目盛りの線だけ、文字は自動)')
        bottom, left, top, right = kinds
        ax.tick_params(axis='x', which='both', direction='in', color=tcol, labelsize=FONT_TICK * tfs * fs,
                       bottom=True, top=top != 'none', labelbottom=bottom == 'labels', labeltop=top == 'labels')
        ax.tick_params(axis='y', which='both', direction='in', color=tcol, labelsize=FONT_TICK * tfs * fs,
                       left=left != 'none', right=right != 'none', labelleft=left == 'labels', labelright=right == 'labels')
        ax.tick_params(which='major', length=4.5 * tsz, width=0.8)
        ax.tick_params(which='minor', length=2.2 * tsz, width=0.6)
        # 格子 (DrawGrid: Automatic・Major = 主目盛りの線、Minor、All・Both = 両方 + 軸、Axis = 原点の線、Axis[x…] = その位置の線)
        gl = op.get('GridLines', [Sym('Automatic'), Sym('Automatic')])
        if isinstance(gl, Sym) or is_num(gl):
            gl = [gl, gl]
        if not (isinstance(gl, list) and len(gl) == 2):
            gl = [Sym('None'), Sym('None')]
        for axis, g in ((0, gl[0]), (1, gl[1])):
            self._grid(fr, axis, g)
        # 枠のラベル {下, 左, 上, 右}
        fl = op.get('FrameLabel')
        if isinstance(fl, list):
            fl = (list(fl) + [''] * 4)[:4]
            pts = fs * FONT_LABEL
            if isinstance(fl[0], str) and fl[0]:
                ax.set_xlabel(sad_text(fl[0]), fontsize=pts * ffs[0])
            if isinstance(fl[1], str) and fl[1]:
                ax.set_ylabel(sad_text(fl[1]), fontsize=pts * ffs[1])
            if isinstance(fl[2], str) and fl[2]:
                ax.annotate(sad_text(fl[2]), (0.5, 1.0), xycoords='axes fraction', xytext=(0, 4), textcoords='offset points',
                            ha='center', va='bottom', fontsize=pts * ffs[2])
            if isinstance(fl[3], str) and fl[3]:
                ax.annotate(sad_text(fl[3]), (1.0, 0.5), xycoords='axes fraction',
                            xytext=(34 if right == 'labels' else 6, 0), textcoords='offset points', rotation=270,
                            ha='left', va='center', fontsize=pts * ffs[3])
        ful = op.get('FrameUnitLabel')
        if isinstance(ful, list) and ful and isinstance(ful[0], str) and ful[0] and bottom == 'labels':
            ax.annotate(sad_text(ful[0]), (1.0, 0.0), xycoords='axes fraction', xytext=(0, -15), textcoords='offset points',
                        ha='right', va='top', fontsize=FONT_TICK * tfs * fs)
        if isinstance(ful, list) and len(ful) > 1 and isinstance(ful[1], str) and ful[1]:
            self.warn('unitlabel-y', 'FrameUnitLabel の左 (y) 側は未対応')

    def _grid(self, fr: _Frame, axis: int, g):
        ax = fr.ax
        items = list(_flat(g if isinstance(g, list) else [g]))
        names = {_name(i) for i in items}
        if is_num(g) and not isinstance(g, bool):
            names = {'Automatic'} if g != 0 else {'None'}
        xy = 'x' if axis == 0 else 'y'
        major = bool(names & {'Automatic', 'All', 'Both', 'Major', 'True'})
        minor = bool(names & {'All', 'Both', 'Minor'})
        axis_line = bool(names & {'Automatic', 'All', 'Axis'})
        if major:
            ax.grid(axis=xy, which='major', color='#a0a0a0', linestyle=(0, (1.0, 2.5)), linewidth=0.6, zorder=0.5)
        if minor:
            ax.grid(axis=xy, which='minor', color='#606060', linestyle=(0, (1.0, 3.0)), linewidth=0.4, zorder=0.5)
        lim = fr.xlim if axis == 0 else fr.ylim

        def line(v):
            if lim[0] < v < lim[1]:
                (ax.axvline if axis == 0 else ax.axhline)(v, color='#606060', linewidth=0.7, zorder=0.6)
        if axis_line:
            line(0.0)
        for it in items:
            if isinstance(it, Expr) and _head(it) == 'Axis':
                for v in _flat(it.args):
                    if is_num(v):
                        line(float(v))
            elif is_num(it) and not isinstance(it, bool) and not is_num(g):
                # 位置のリスト {1,2,5}: Canvas.n の DrawGrid は何も描かないが、Mathematica の流儀どおりその位置に線を引く (拡張)
                line(float(it))

    # ---------------------------------------------------------------- 要素
    def _primitive(self, fr: _Frame, p):
        if not isinstance(p, Expr):
            self.warn('prim-' + type(p).__name__, f'Graphics の要素 {p!r} (型 {type(p).__name__}) は未対応で飛ばした')
            return
        h = _head(p)
        fn = getattr(self, '_p_' + h, None)
        if h == 'Hook':
            return                                    # 対話の束縛 (OpticsPlot の InfoLabel など)。静的な図では要らない
        if fn is None:
            self.warn('prim-' + h, f'Graphics の要素 {h}[…] は未対応で飛ばした')
            return
        try:
            fn(fr, p)
        except _Unsupported as e:
            self.warn('unsup-' + str(e), f'{h}: {e} は未対応で飛ばした')

    def _axes(self, fr: _Frame):
        """要素を描く Axes。枠が空 (範囲が空) なら、窓いっぱいの見えない Axes を作る"""
        if fr.ax is None:
            s = self.sess
            if fr.window is None:
                fr.window = make_window(s.plot_region, None, None)
            x1, x2, y1, y2 = fr.window
            ax = s.fig.add_axes([x1 / CANVAS[0], y1 / CANVAS[1], (x2 - x1) / CANVAS[0], (y2 - y1) / CANVAS[1]])
            ax.set_xlim(*fr.xlim)
            ax.set_ylim(*fr.ylim)
            ax.axis('off')
            fr.ax = ax
        return fr.ax

    def _col(self, fr: _Frame, c, default='black', idx=None):
        """色。文字列は X11 の名前、評価されていない式 (PlotColorDefault[…]) やリストは、既定の色の表から枠の中の順番で選ぶ"""
        if isinstance(c, str):
            if c == '':
                return None
            r = sad_color(c)
            if r is None:
                self.warn('color-' + c, f'色 {c!r} を解釈できない (黒にした)')
                return (0.0, 0.0, 0.0)
            return r
        if isinstance(c, list) and c and all(isinstance(a, str) for a in c):
            i = fr.auto if idx is None else idx
            return sad_color(c[i % len(c)])
        if isinstance(c, Expr) or isinstance(c, list):
            i = fr.auto if idx is None else idx
            return sad_color(PALETTE[i % len(PALETTE)])
        r = sad_color(default) if isinstance(default, str) else default
        return r

    @staticmethod
    def _unresolved(c) -> bool:
        return isinstance(c, Expr) or (isinstance(c, list) and not all(isinstance(a, str) for a in c))

    @staticmethod
    def _points(body) -> Optional[List[List[float]]]:
        if not isinstance(body, list) or not body:
            return None
        if not all(isinstance(r, list) and len(r) >= 2 for r in body):
            return None
        return body

    @staticmethod
    def _dash(dash, lw: float):
        """Dashing → matplotlib の linestyle。1・{1} は実線、実数 d は {d, d/3}、リストは (線, 空き, …) のインチ"""
        if is_num(dash):
            dash = [1.0] if dash == 1 else [float(dash), float(dash) / 3.0]
        if isinstance(dash, list) and dash and isinstance(dash[0], list):
            dash = dash[0]
        if not isinstance(dash, list) or not dash or dash == [1] or dash == [1.0] or not all(is_num(a) for a in dash):
            return '-'
        seq = [max(float(a), 1e-3) * PT_PER_INCH / max(lw, 0.1) for a in dash]
        if len(seq) % 2:
            seq = seq * 2
        return (0, tuple(seq))

    def _lw(self, thick) -> float:
        t = float(thick) if is_num(thick) else 1.0
        return max(t, 0.2) * LW_PT

    def _add(self, fr: _Frame, art, kind: str, tag: str = ''):
        art.set_zorder(fr.z())
        fr.rec[kind].append(art)

    def _p_Line(self, fr: _Frame, p: Expr):
        body = p.args[0] if p.args else []
        o = _opts(p.args[1:], self.I.eval)
        pts = self._points(body)
        if pts is None:
            return
        ax = self._axes(fr)
        pj = o.get('PlotJoined', 1.0)
        pjn = _name(pj)
        plot_pts = _truth(o.get('Plot', 0.0))
        pc_raw = o.get('PlotColor', 'Black')
        ptc_raw = o.get('PointColor', 'forest green')
        idx = fr.auto
        if self._unresolved(pc_raw) or self._unresolved(ptc_raw):
            fr.auto += 1
        pc = self._col(fr, pc_raw, 'black', idx)
        lw = self._lw(o.get('Thickness', 1.0))
        ls = self._dash(o.get('Dashing', [1.0]), lw)
        clip = _truth(o.get('Clip', 1.0), True)
        xs = [float(r[0]) for r in pts]
        ys = [float(r[1]) for r in pts]
        if pjn == 'Step':
            sr = o.get('StepRatio', 1.0)
            sr = float(sr) if is_num(sr) else 1.0
            sx, sy = _step(xs, ys, sr)
        else:
            sx, sy = xs, ys
        if pj is True or pjn == 'Step' or (is_num(pj) and pj != 0) or pjn == 'True':
            (ln,) = ax.plot(sx, sy, color=pc, linewidth=lw, linestyle=ls, solid_capstyle='butt', clip_on=clip)
            ln.set_gid('line')
            self._add(fr, ln, 'lines')
        elif not (pj is False or (is_num(pj) and pj == 0) or pjn == 'False'):
            raise _Unsupported(f'PlotJoined->{pj!r}')
        arrow = o.get('Arrow', Null)
        if arrow is not Null and _name(arrow) not in ('Null', 'None', 'False') and not (is_num(arrow) and arrow == 0):
            self.warn('arrow', 'Line の Arrow・ArrowShape は未対応 (線だけ描く)')
        if plot_pts:
            self._draw_points(fr, ax, xs, ys, o, idx)
        self._errorbars(fr, ax, pts, o, pc)
        self._legend(fr, o, pc, ls, lw, plot_pts, pj, idx)

    def _p_Point(self, fr: _Frame, p: Expr):
        body = p.args[0] if p.args else []
        o = _opts(p.args[1:], self.I.eval)
        pts = self._points(body)
        if pts is None:
            return
        ax = self._axes(fr)
        pc_raw = o.get('PlotColor', 'black')
        ptc_raw = o.get('PointColor', 'forest green')
        idx = fr.auto
        if self._unresolved(ptc_raw) or self._unresolved(pc_raw):
            fr.auto += 1
        xs = [float(r[0]) for r in pts]
        ys = [float(r[1]) for r in pts]
        self._draw_points(fr, ax, xs, ys, o, idx)
        pc = self._col(fr, pc_raw, 'black', idx)
        self._errorbars(fr, ax, pts, o, pc)
        self._legend(fr, o, pc, '-', 1.0, True, 0.0, idx)

    def _draw_points(self, fr, ax, xs, ys, o, idx):
        sym = o.get('PointSymbol', '1O')
        if isinstance(sym, list):
            sym = sym[0] if sym else '1O'
        bar = _name(sym) == 'Bar' or sym in ('BA', 'BR')
        size = _first(o.get('PointSize', 1.0), 1.0)
        size = float(size) if is_num(size) else 1.0
        ptc = self._col(fr, o.get('PointColor', 'forest green'), 'forest green', idx)
        pbc_raw = o.get('PointBorderColor', Sym('Automatic'))
        pbc = self._col(fr, pbc_raw, 'black', idx) if (isinstance(pbc_raw, str) and pbc_raw) else ptc
        if bar:
            real = [(x, y) for x, y in zip(xs, ys) if abs(x) < 1e250 and abs(y) < 1e250]
            if real:                               # 位置が 1e300 の点は ColumnPlot の凡例の見本だけ (描かない)
                for x, y in real:
                    v = ax.vlines([x], [0.0], [y], colors=[ptc], linewidth=max(size, 1.0) * LW_PT * 2)
                    v.set_zorder(fr.z())
                    fr.rec['lines'].append(v)
                self.warn('bar', 'PointSymbol->Bar は基線 0 からの垂線で近似 (Canvas.n の基線 Offset は未対応)')
            return
        mk = _MARKERS.get(str(sym))
        if mk is None:
            self.warn('symbol-' + str(sym), f'PointSymbol {sym!r} は未対応 (○ にした)')
            mk = 'o'
        ms = 6.0 * size * math.sqrt(512.0 / max(512.0, len(xs)))
        (ln,) = ax.plot(xs, ys, linestyle='none', marker=mk, markersize=ms, markerfacecolor=ptc, markeredgecolor=pbc,
                        markeredgewidth=0.6, clip_on=True)
        ln.set_gid('points')
        self._add(fr, ln, 'points')

    def _errorbars(self, fr, ax, pts, o, pc):
        ncol = min(len(r) for r in pts)
        if ncol < 3:
            return
        if ncol == 3:
            self.warn('errbar3', '3 列のデータの誤差棒は未対応 ({x,y,dx,dy} の 4 列だけ)')
            return
        ebc = o.get('ErrorBarColor', Sym('Default'))
        col = self._col(fr, ebc, 'black') if isinstance(ebc, str) else pc
        ebl = _first(o.get('ErrorBarTickSize', 1.0), 1.0)
        ebt = _first(o.get('ErrorBarThickness', 1.0), 1.0)
        eb = ax.errorbar([float(r[0]) for r in pts], [float(r[1]) for r in pts], xerr=[abs(float(r[2])) for r in pts],
                         yerr=[abs(float(r[3])) for r in pts], fmt='none', ecolor=col, elinewidth=self._lw(ebt),
                         capsize=2.0 * (float(ebl) if is_num(ebl) else 1.0))
        for a in ([eb.lines[0]] if eb.lines[0] is not None else []) + list(eb.lines[1]) + list(eb.lines[2]):
            a.set_zorder(fr.z())
        fr.rec['errorbars'].append(eb)

    def _legend(self, fr, o, pc, ls, lw, plot_pts, pj, idx):
        leg = o.get('Legend', 0.0)
        if leg is False or (is_num(leg) and leg == 0) or _name(leg) in ('False', 'None', 'Null') or leg is Null:
            return
        label = leg
        extra = {}
        if isinstance(leg, list) and leg:
            label = leg[0]
            extra = _opts(leg[1:], self.I.eval)
        if not isinstance(label, str):
            return
        joined = not (pj is False or (is_num(pj) and pj == 0) or _name(pj) == 'False')
        ptc = self._col(fr, o.get('PointColor', 'forest green'), 'forest green', idx)
        sym = o.get('PointSymbol', '1O')
        sym = sym[0] if isinstance(sym, list) and sym else sym
        bar = _name(sym) == 'Bar' or sym in ('BA', 'BR')
        entry = dict(label=label, color=pc, ls=ls, lw=lw, joined=joined, points=bool(plot_pts), ptc=ptc,
                     marker=_MARKERS.get(str(sym), 'o'), extra=extra, bar=bar)
        if fr.legend and fr.legend[-1][0] == label:                # AppendLegend: 同じ名前が続けば 1 つにまとめる
            prev = fr.legend[-1][1]
            prev['joined'] = prev['joined'] or joined
            prev['points'] = prev['points'] or bool(plot_pts)
            return
        fr.legend.append((label, entry))

    def _p_Rectangle(self, fr: _Frame, p: Expr):
        a = p.args
        if len(a) >= 2 and any(isinstance(x, Expr) and _head(x) == 'Graphics' or
                               (isinstance(x, list) and any(isinstance(y, Expr) and _head(y) == 'Graphics' for y in x))
                               for x in a[2:]):
            self._nested(fr, p)
            return
        if len(a) < 2 or not (isinstance(a[0], list) and isinstance(a[1], list) and len(a[0]) >= 2 and len(a[1]) >= 2):
            raise _Unsupported('Rectangle の引数')
        o = _opts(a[2:], self.I.eval)
        from matplotlib.patches import Rectangle
        ax = self._axes(fr)
        x0, y0 = float(a[0][0]), float(a[0][1])
        x1, y1 = float(a[1][0]), float(a[1][1])
        fc = o.get('FillColor', Null)
        fcol = sad_color(fc) if isinstance(fc, str) else None
        pcn = o.get('PlotColor', '')
        ecol = sad_color(pcn) if isinstance(pcn, str) and pcn else None
        if isinstance(o.get('MeshStyle', Null), str):
            self.warn('stipple', 'Rectangle の MeshStyle (stipple の模様) は未対応 (塗りつぶしだけ)')
        r = Rectangle((min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0), facecolor=fcol if fcol is not None else 'none',
                      edgecolor=ecol if ecol is not None else 'none', linewidth=0.6 if ecol is not None else 0.0, clip_on=True)
        r.set_gid('rect')
        ax.add_patch(r)
        self._add(fr, r, 'rects')

    def _nested(self, fr: _Frame, p: Expr):
        """Rectangle[min, max, Graphics[…] (かそのリスト), 選択肢…]: Canvas.n の Primitive[Rectangle[min,max,g__,opt]] と同じく
        Show1$[False, {min,max}, g, opt] を呼び、枠 (rect は PlotRegion の窓の割合) をその中に描く"""
        a = p.args
        s = self.sess
        rules = [x for x in a[2:] if isinstance(x, Rule)]
        g = [x for x in a[2:] if not isinstance(x, Rule)]
        pos = _opts(rules).get('Position')
        rect = [list(a[0]), list(a[1])]
        if isinstance(pos, Expr) and _head(pos) in ('Right', 'Left', 'Above', 'Below') and fr.window is not None:
            sep = float(pos.args[0]) if pos.args and is_num(pos.args[0]) else 0.0
            rect = self._rect_position(fr, rect, _head(pos), sep)
        elif pos is not None and _name(pos) not in ('Default', 'Null'):
            self.warn('rect-position', f'Rectangle の Position->{pos!r} は未対応 (Default にした)')
        saved = s.cur
        try:
            self.I.eval(Expr(S('Show1$'), [0.0, rect] + g + rules))
        finally:
            s.cur = saved

    def _rect_position(self, fr: _Frame, rect, how: str, sep: float):
        """RectPosition (Canvas.n): いまの窓 (PlotRegion の窓に対する割合 w) の右・左・上・下に、rect と同じ大きさの箱を置く (カラーバーなど)"""
        (px1, px2), (py1, py2) = self.sess.plot_region
        x1p = _clamp(px1 * RAW[0] + ORIGIN[0], CANVAS[0] * MARGIN[0][0], CANVAS[0] * MARGIN[0][1])
        x2p = _clamp(px2 * RAW[0] + ORIGIN[0], CANVAS[0] * MARGIN[0][0], CANVAS[0] * MARGIN[0][1])
        y1p = _clamp(py1 * RAW[1] + ORIGIN[1], CANVAS[1] * MARGIN[1][0], CANVAS[1] * MARGIN[1][1])
        y2p = _clamp(py2 * RAW[1] + ORIGIN[1], CANVAS[1] * MARGIN[1][0], CANVAS[1] * MARGIN[1][1])
        wx1, wx2, wy1, wy2 = fr.window
        w = [[(wx1 - x1p) / (x2p - x1p), (wx2 - x1p) / (x2p - x1p)], [(wy1 - y1p) / (y2p - y1p), (wy2 - y1p) / (y2p - y1p)]]
        mn, mx = rect
        if how == 'Right':
            return [[w[0][1] + sep, w[1][0]], [w[0][1] + sep + mx[0] - mn[0], w[1][1]]]
        if how == 'Left':
            return [[mn[0] - mx[0] - sep + w[0][0], w[1][0]], [w[0][0] - sep, w[1][1]]]
        if how == 'Above':
            return [[w[0][0], w[1][1] + sep], [w[0][1], w[1][1] + sep + mx[1] - mn[1]]]
        return [[w[0][0], w[1][0] - sep + mn[1] - mx[1]], [w[0][1], w[1][0] - sep]]

    def _p_Polygon(self, fr: _Frame, p: Expr):
        body = p.args[0] if p.args else []
        o = _opts(p.args[1:], self.I.eval)
        pts = self._points(body)
        if pts is None:
            return
        from matplotlib.patches import Polygon
        ax = self._axes(fr)
        fc = o.get('FillColor', Null)
        fcol = sad_color(fc) if isinstance(fc, str) else None
        pc_raw = o.get('PlotColor', 'black')
        idx = fr.auto
        if self._unresolved(pc_raw):
            fr.auto += 1
        pc = self._col(fr, pc_raw, 'black', idx)
        lw = self._lw(o.get('Thickness', 1.0))
        pj = o.get('PlotJoined', 1.0)
        border = pj is True or (is_num(pj) and pj != 0) or _name(pj) == 'True'
        poly = Polygon([(float(r[0]), float(r[1])) for r in pts], closed=True,
                       facecolor=fcol if fcol is not None else 'none', edgecolor=pc if border else 'none',
                       linewidth=lw if border else 0.0, linestyle=self._dash(o.get('Dashing', [1.0]), lw), clip_on=True)
        poly.set_gid('polygon')
        ax.add_patch(poly)
        self._add(fr, poly, 'polys')
        if _truth(o.get('Plot', 0.0)) and _name(o.get('Plot')) != 'False':
            self._draw_points(fr, ax, [float(r[0]) for r in pts], [float(r[1]) for r in pts], o, idx)

    def _p_Circle(self, fr: _Frame, p: Expr):
        a = p.args
        if len(a) < 2 or not isinstance(a[0], list) or len(a[0]) < 2:
            raise _Unsupported('Circle の引数')
        from matplotlib.patches import Ellipse
        o = _opts(a[2:], self.I.eval)
        ax = self._axes(fr)
        cx, cy = float(a[0][0]), float(a[0][1])
        r = a[1]
        rx, ry = (float(r[0]), float(r[1])) if isinstance(r, list) and len(r) >= 2 else (float(r), float(r))
        fc = o.get('FillColor', Null)
        pc = o.get('PlotColor', Null)
        th = o.get('Thickness', Null)
        e = Ellipse((cx, cy), 2 * rx, 2 * ry, facecolor=sad_color(fc) if isinstance(fc, str) and fc else 'none',
                    edgecolor=(sad_color(pc) if isinstance(pc, str) and pc else (0, 0, 0)),
                    linewidth=self._lw(th) if is_num(th) else LW_PT, clip_on=True)
        e.set_gid('oval')
        ax.add_patch(e)
        self._add(fr, e, 'ovals')

    def _coord(self, fr: _Frame, v, axis: int):
        """Text の座標の 1 成分 → (値, 種類)。Scaled = データ座標、WindowScaled = 窓の座標、RectScaled・CanvasScaled・数 = キャンバスのインチ・画素"""
        if isinstance(v, Expr):
            h = _head(v)
            if h in ('Scaled', 'WindowScaled', 'RectScaled', 'CanvasScaled') and v.args and is_num(v.args[0]):
                return float(v.args[0]), {'Scaled': 'data', 'WindowScaled': 'win', 'RectScaled': 'rect', 'CanvasScaled': 'canvas'}[h]
            raise _Unsupported(h)
        if is_num(v):
            return float(v), 'inch'
        raise _Unsupported(type(v).__name__)

    def _trans(self, fr: _Frame, ax, v, axis: int):
        """(値, 種類) → (その軸の変換, 変換後の値)"""
        val, kind = v
        fig = ax.figure
        if kind == 'data':
            return ax.transData, val
        if kind == 'win':
            return ax.transAxes, (val - ORIGIN[axis]) / RAW[axis]
        if kind == 'rect':
            return fig.transFigure, (val * RAW[axis] + ORIGIN[axis]) / CANVAS[axis]
        if kind == 'inch':
            return fig.transFigure, val / CANVAS[axis]
        # canvas (画素。DPI 72 とみなす。y は上から)
        inch = val / 72.0
        return fig.transFigure, (inch / CANVAS[0]) if axis == 0 else (CANVAS[1] - inch) / CANVAS[1]

    def _p_Text(self, fr: _Frame, p: Expr):
        a = p.args
        if not a or not isinstance(a[0], list) or len(a[0]) < 2:
            raise _Unsupported('Text の引数')
        body = a[0]
        txt, pos = body[0], body[1]
        if txt == '' or (isinstance(txt, str) and txt == ''):
            return
        if not isinstance(pos, list) or len(pos) < 2:
            raise _Unsupported('Text の位置')
        o = _opts(a[1:], self.I.eval)
        ax = self._axes(fr)
        cx, cy = self._coord(fr, pos[0], 0), self._coord(fr, pos[1], 1)
        if _truth(o.get('Clip', 0.0)) and cx[1] == 'data' and cy[1] == 'data':
            if not (fr.xlim[0] <= cx[0] <= fr.xlim[1] and fr.ylim[0] <= cy[0] <= fr.ylim[1]):
                return
        from matplotlib.transforms import blended_transform_factory
        tx, vx = self._trans(fr, ax, cx, 0)
        ty, vy = self._trans(fr, ax, cy, 1)
        tr = blended_transform_factory(tx, ty)
        size = o.get('TextSize', 1.0)
        size = float(size) if is_num(size) else 1.0
        rot = o.get('Rotation', None)
        if not is_num(rot):
            rot = o.get('TextRotate', 0.0)
        rot = float(rot) % 360.0 if is_num(rot) else 0.0
        anch = o.get('Anchor')
        if not isinstance(anch, str):
            al = o.get('TextAlign', 'w')
            anch = str(al).lower() if isinstance(al, str) else 'w'
        ha, va = _ANCHOR.get(anch.lower(), ('left', 'center'))
        just = o.get('Justify', o.get('TextJustify', Null))
        ma = just if isinstance(just, str) and just in ('left', 'right', 'center') else ha
        col = self._col(fr, o.get('PlotColor', 'black'), 'black', 0)
        tag = o.get('Tags')
        tag = tag if isinstance(tag, str) else ''
        if _truth(o.get('UseBG', 0.0)):
            self.warn('usebg', 'Text の UseBG->True (文字の背景) は未対応 (背景なし)')
        if o.get('TextFont') not in (None, Null) and isinstance(o.get('TextFont'), (str, list)) and \
                o.get('TextFont') != 'helvetica neue':
            self.warn('textfont', 'Text の TextFont (フォントの名前・書体) は未対応 (既定のフォント)')
        t = ax.text(vx, vy, sad_text(txt if isinstance(txt, str) else str(txt)), transform=tr, ha=ha, va=va,
                    rotation=rot, rotation_mode='default', fontsize=FONT_TEXT * size * fr.font_scale,
                    color=col or 'black', multialignment=ma, clip_on=False)
        t.set_gid('text')
        t.set_zorder(fr.z())
        fr.rec['texts'].append(t)
        fr.texts.append((t, tag))

    # ---------------------------------------------------------------- 仕上げ
    def _finish(self):
        s, self.sess = self.sess, None
        if s is None:
            return
        import matplotlib.pyplot as plt
        fig = s.fig
        for fr in s.frames:
            self._draw_legend(fr)
        self._remove_overlaps(s)
        self._remove_tick_overlaps(s)
        self._draw_title(s)
        fig._pysad_display = {'session': s, 'frames': s.frames, 'title': s.title, 'warnings': s.warnings}
        self.last_fig = fig
        global _LAST
        _LAST = fig
        fp = tuple(s.parts)
        if fp == self.last_fp and not self.override:
            plt.close(fig)                                          # 直前と同じ内容 (ListPlot.n の 2 回目の Show): 出力しない
            return
        self.last_fp = fp
        self._output(fig)

    def _draw_legend(self, fr: _Frame):
        if not fr.legend or fr.ax is None:
            return
        from matplotlib.lines import Line2D
        handles, labels = [], []
        for label, e in fr.legend:
            mk = e['marker'] if e['points'] else ''
            ls = e['ls'] if e['joined'] else 'none'
            if e.get('bar'):
                handles.append(Line2D([], [], color=e['ptc'], linestyle='-', linewidth=7.0))
                labels.append(sad_text(label))
                continue
            handles.append(Line2D([], [], color=e['color'], linestyle=ls, linewidth=e['lw'], marker=mk,
                                  markersize=5.0, markerfacecolor=e['ptc'] if e['points'] else None,
                                  markeredgecolor=e['ptc'] if e['points'] else None))
            labels.append(sad_text(label))
        bg = fr.legend_bg
        fc = sad_color(bg) if isinstance(bg, str) else None
        lg = fr.ax.legend(handles, labels, loc='best', fontsize=FONT_LEGEND * fr.legend_scale, framealpha=1.0 if fc else 0.85,
                          facecolor=fc if fc is not None else 'white', edgecolor='#666666', handlelength=2.6)
        lg.set_zorder(9)
        fr.rec.setdefault('legends', []).append(lg)

    def _remove_overlaps(self, s: _Session):
        """Remove$Overlap[RemoveOverlap の Tags]: その Tags の文字を書いた順に見て、残した文字の外接矩形の和に重なるものを消す"""
        todo = [fr for fr in s.frames if fr.rmo and any(t == fr.rmo for _, t in fr.texts)]
        if not todo:
            return
        try:
            rend = s.fig.canvas.get_renderer()
        except Exception:                                           # noqa: BLE001
            return
        for fr in todo:
            items = [t for t, tag in fr.texts if tag == fr.rmo]
            keep = []
            union = None
            for t in items:
                bb = t.get_window_extent(rend)
                box = (bb.x0, bb.y0, bb.x1, bb.y1)
                if union is not None and not (box[0] > union[2] or union[0] > box[2] or box[1] > union[3] or union[1] > box[3]):
                    t.remove()
                    fr.rec['texts'].remove(t)
                    fr.rec['texts_removed'].append(t)
                    continue
                keep.append(t)
                union = box if union is None else (min(union[0], box[0]), min(union[1], box[1]),
                                                  max(union[2], box[2]), max(union[3], box[3]))
            fr.texts = [(t, tag) for t, tag in fr.texts if t in keep or tag != fr.rmo]

    def _remove_tick_overlaps(self, s: _Session):
        """DeleteOverlapLabel: 枠をまたいで重なる目盛りの文字 (上の枠の下端と下の枠の上端など) は、あとの方を消す"""
        axs = [fr.ax for fr in s.frames if fr.ax is not None and fr.ax.axison]
        if len(axs) < 2:
            return
        try:
            s.fig.canvas.draw()
            rend = s.fig.canvas.get_renderer()
        except Exception:                                           # noqa: BLE001
            return
        kept: List[tuple] = []
        for ax in axs:
            for axis in (ax.xaxis, ax.yaxis):
                for tk in axis.get_major_ticks():
                    for lab in (tk.label1, tk.label2):
                        if not lab.get_visible() or not lab.get_text() or not lab.get_figure():
                            continue
                        bb = lab.get_window_extent(rend)
                        box = (bb.x0, bb.y0, bb.x1, bb.y1)
                        if any(not (box[0] > k[2] or k[0] > box[2] or box[1] > k[3] or k[1] > box[3]) for k in kept):
                            lab.set_visible(False)
                            continue
                        kept.append(box)

    def _draw_title(self, s: _Session):
        if not s.title:
            return
        tops = [fr.ax.get_position().y1 for fr in s.frames if fr.ax is not None]
        lefts = [(fr.ax.get_position().x0 + fr.ax.get_position().x1) / 2 for fr in s.frames if fr.ax is not None]
        y = (max(tops) + 0.012) if tops else 0.95
        x = lefts[0] if lefts else 0.5
        s.fig.text(x, y, sad_text(s.title), ha='center', va='bottom', fontsize=FONT_TITLE * s.title_scale)

    def _output(self, fig):
        from . import draw_mpl
        import matplotlib.pyplot as plt
        opts = self._global_options()
        file_ = self.override.get('File', opts.get('FILE'))
        file_ = file_ if isinstance(file_, str) and file_ else os.environ.get('PYSAD_DISPLAY_FILE', '')
        show = self.override.get('Show', opts.get('SHOW'))
        dpi = self.override.get('DPI', opts.get('DPI'))
        dpi = float(dpi) if is_num(dpi) else float(os.environ.get('PYSAD_DISPLAY_DPI', '100'))
        mpl = draw_mpl._import_mpl()
        interactive = draw_mpl._interactive_backend(mpl)
        if file_:
            path = draw_mpl._next_path(file_)
            try:
                ext = os.path.splitext(path)[1].lower().lstrip('.')
                if ext not in ('png', 'pdf', 'svg', 'jpg', 'jpeg', 'eps', 'ps'):
                    raise ValueError(f'保存の拡張子は png・pdf・svg ({path})')
                os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
                fig.savefig(path, dpi=dpi, bbox_inches='tight', pad_inches=0.15)
                self.write(f'[pysad] MatplotlibDrawer: {path} に保存しました')
            except Exception as e:                                  # noqa: BLE001
                self.write(f'[pysad] MatplotlibDrawer: 保存できません: {e}')
            if _truth(show, False) and interactive:
                draw_mpl._show(fig, mpl)
            else:
                plt.close(fig)
        elif (show is None or _truth(show, True)) and interactive:
            draw_mpl._show(fig, mpl)
        else:
            plt.close(fig)                                          # Agg など非対話: 何も出さない (last_figure() で取り出せる)

    def _global_options(self) -> Dict[str, str]:
        lst = self.I.globals.get('MatplotlibDrawer$Option')
        out: Dict[str, object] = {}
        for r in lst if isinstance(lst, list) else []:
            if isinstance(r, Rule):
                k = str(_name(r.lhs) or r.lhs).upper()
                out.setdefault(k, r.rhs)
        return out


Default = Sym('Default')
_LAST = None


def _step(xs, ys, sr):
    """Canvas$StepData: 階段 (StepRatio sr: 0 = 先に縦、1 = 先に横)"""
    n = len(xs)
    if n < 2:
        return xs, ys
    if sr < 0:
        sr = 0.0
    if sr > 1:
        sr = 1.0
    j = [xs[i] + sr * (xs[i + 1] - xs[i]) for i in range(n - 1)]
    ax_ = [xs[0]] + j
    bx_ = j + [xs[-1]]
    ox, oy = [], []
    for i in range(n):
        ox += [ax_[i], bx_[i]]
        oy += [ys[i], ys[i]]
    return ox, oy


def _expand_frame_ticks(ft):
    """ExpandFrameTicks (Canvas.n)"""
    if _is_auto(ft):
        return [True, True, Sym('Ticks'), Sym('Ticks')]
    if not isinstance(ft, list):
        return [ft] * 4
    f = list(ft)
    n = len(f)
    if n >= 1 and _is_auto(f[0]):
        f[0] = True
    if n >= 2 and _is_auto(f[1]):
        f[1] = True
    if n >= 3 and _is_auto(f[2]):
        f[2] = Sym('Ticks')
    if n >= 4 and _is_auto(f[3]):
        f[3] = Sym('Ticks')
    if n == 1:
        f += [True, Sym('Ticks'), Sym('Ticks')]
    elif n == 2:
        f += [Sym('Ticks'), Sym('Ticks')]
    elif n == 3:
        f += [Sym('Ticks')]
    return f[:4]


def _tick_kind(v) -> str:
    """FrameTicks の 1 辺: True・Both = 目盛りと文字、Ticks (と {Ticks, 関数}) = 目盛りの線だけ、False = なし"""
    if isinstance(v, list) and v and _name(v[0]) == 'Ticks':
        return 'ticks'
    n = _name(v)
    if n == 'Ticks':
        return 'ticks'
    if n == 'Both' or v is True:
        return 'labels'
    if is_num(v):
        return 'labels' if v != 0 else 'none'
    if n == 'True':
        return 'labels'
    return 'none'


# ---------------------------------------------------------------------------------------------------------------------
# 5. インタプリタへの登録・環境変数・Python API
# ---------------------------------------------------------------------------------------------------------------------
def display_mode_from_env() -> Optional[str]:
    """環境変数 PYSAD_DISPLAY → 'mpl'・'tk'・'none'。未設定・知らない値は None (従来の決め方)"""
    v = os.environ.get('PYSAD_DISPLAY', '').strip().lower()
    if v in ('mpl', 'matplotlib'):
        return 'mpl'
    if v in ('tk', 'canvas', 'canvasdrawer'):
        return 'tk'
    if v in ('none', 'off', 'identity'):
        return 'none'
    return None


def install(I):
    """組み込み関数 MatplotlibDrawer を I に登録する (matplotlib はまだ読まない)。何度呼んでもよい"""
    drv = getattr(I, '_mpl_display', None)
    if drv is None:
        drv = _Driver(I)
        I._mpl_display = drv
    I.builtins['MatplotlibDrawer'] = drv.call
    return drv


def MatplotlibDrawer(I, *args):
    """Python から直接 (opt, list, init, w, rect, last, ini) を渡す"""
    return install(I).call(list(args))


def last_figure():
    """最後に描いた Figure (保存して閉じたあとでも、線・箱・文字は `fig._pysad_display['frames']` の rec から取れる)"""
    return _LAST


def render(I, g, save: Optional[str] = None, show: Optional[bool] = None, dpi: Optional[float] = None):
    """Graphics[…] の式 g (Plot・ListPlot・OpticsPlot などを DisplayFunction->Identity で作ったもの) を matplotlib に描いて
    Figure を返す (matplotlib が無ければ ImportError)。save があればその名前で保存 (番号つき)。show=True で対話表示"""
    from . import draw_mpl
    if draw_mpl._import_mpl() is None:
        raise ImportError('matplotlib がありません (pip install matplotlib)')
    drv = install(I)
    drv.override = {k: v for k, v in (('File', save), ('Show', show), ('DPI', dpi)) if v is not None}
    if save is None and show is None:
        drv.override['Show'] = False            # render() は保存も表示もしないとき、図を閉じて返すだけ
    old = I.globals.get('$DisplayFunction')
    I.globals['$DisplayFunction'] = Sym('MatplotlibDrawer')
    drv.last_fig = None
    drv.last_fp = None
    try:
        I.eval(Expr(S('Show'), [g]))
    finally:
        drv.override = {}
        if old is None:
            I.globals.pop('$DisplayFunction', None)
        else:
            I.globals['$DisplayFunction'] = old
    return drv.last_fig
