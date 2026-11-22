# -*- coding: utf-8 -*-
"""FFS の `DRAW` コマンドを matplotlib で実装する (v2.4.1)。

SAD K64 の `DRAW` は `src/tffsa.f` が `System`CANVASDRAW[]` を評価し、実体は `Packages/CANVASDRAW.n` (K. Oide) が
`OpticsPlot[list, Region->{r1,r2}, Lattice->…, Names->…, GridLines->{None,Automatic}, InfoLabel->True, PlotLabel->TITLE]` を呼ぶ。
`$DisplayFunction` が `CanvasDrawer` (Tk) でなければ何もしない。pysad は Tk のキャンバスの代わりに matplotlib に描く。

構成 (matplotlib は遅延 import。無ければ警告して何もしない。データの作成まではこのモジュールだけで動く)
  1. parse_draw          CANVASDRAW.n の構文解析の移植 (SeparateAmpsnd、r1・r2、`&`、LAT・要素名)
  2. build_plot_data     OpticsPlot.n の OpticsList / OpticsPlot1$ / LatticePlot の移植 (点列・単位・ラティスの箱と名前)
  3. render              PlotData を matplotlib の Figure にする
  4. draw / run_draw     Python API と FFS コマンドの入口 (保存・表示)

OpticsPlot.n が `Twiss[f, i+j/ndiv]` で作る値と同じ呼び出し (FFS.twiss_value) で点列を作るので、`OpticsPlot[…, DisplayFunction->Identity]`
が返す Graphics の Line の点列と一致する (試験 `test_v25_draw_mpl.py`)。

使い方 (Python):
    from pysad.ffscmd import run_sad_string
    from pysad import draw_mpl
    ffs, I = run_sad_string(script)              # USE・CALC まで
    fig = draw_mpl.draw(ffs, "BX BY & EX LAT", region=(1, 40), save="optics.png")

FFS から: `DRAW [r1 [r2]] 関数 [関数 …] [& 関数 …] … [LAT|LATTICE|要素名];`
保存先は環境変数 `PYSAD_DRAW_FILE` (拡張子 png・pdf・svg)、または SAD 側の `Draw$Option={File->"x.png"}`。
"""
from __future__ import annotations

import atexit
import contextlib
import math
import os
import re
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = ['parse_draw', 'DrawSpec', 'build_plot_data', 'PlotData', 'render', 'draw', 'run_draw', 'is_label',
           'element_info', 'twiss_unit', 'graphics_lines', 'graphics_rectangles', 'graphics_texts']

# ---------------------------------------------------------------------------------------------------------------------
# 1. 関数名・ラベル・単位 (OpticsPlot.n の Twiss$Functions / Twiss$Label / Twiss$Unit)
# ---------------------------------------------------------------------------------------------------------------------
_TWISS_X = ['AX', 'BX', 'NX', 'EX', 'EPX', 'PEX', 'PEPX', 'DX', 'DPX', 'GMX', 'ZX', 'ZPX', 'PZX', 'PZPX']
_TWISS_Y = [n[:-1] + 'Y' if n.endswith('X') else n for n in _TWISS_X]
_TWISS_Y = ['AY', 'BY', 'NY', 'EY', 'EPY', 'PEY', 'PEPY', 'DY', 'DPY', 'GMY', 'ZY', 'ZPY', 'PZY', 'PZPY']
_TWISS_OTHER = ['AZ', 'BZ', 'NZ', 'GMZ', 'DZ', 'DDP', 'R1', 'R2', 'R3', 'R4', 'DETR']
# OpticsPlot.n の Twiss$Functions と同じ集合 (順は違う)
TWISS_FUNCTIONS = frozenset(_TWISS_X + _TWISS_Y + _TWISS_OTHER)
# Twiss$Functions に R・D を付けたもの (参照光学・参照との差。点は要素ごとの 1 点だけ)
RTWISS_FUNCTIONS = frozenset('R' + n for n in TWISS_FUNCTIONS)
DTWISS_FUNCTIONS = frozenset('D' + n for n in TWISS_FUNCTIONS)
XTWISS_FUNCTIONS = RTWISS_FUNCTIONS | DTWISS_FUNCTIONS | {'BMAGX', 'BMAGY', 'BMAGZ'}
LINE_FUNCTIONS = frozenset(['GAMMA', 'GAMMABETA', 'S', 'BZS', 'BZT', 'GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3',
                            'OGX', 'OGY', 'OGZ', 'OCHI1', 'OCHI2', 'OCHI3'])
# Twiss$Label が空でない名前 (DRAW が関数として受け取る語)。LINE$Functions のうち GAMMA・GAMMABETA・S にラベルは無い
_LABELLED_LINE = LINE_FUNCTIONS - {'GAMMA', 'GAMMABETA', 'S'}
_SIZE_RX = re.compile(r'^(SIG|SIZE)([XYZPDE]{1,4})$')


def is_label(name: str) -> bool:
    """`Twiss$Label[name] <=> ""` (CANVASDRAW.n が関数名かどうかを決める)。SIG・SIZE で始まるサイズの名前も含む"""
    n = name.upper()
    return n in TWISS_FUNCTIONS or n in XTWISS_FUNCTIONS or n in _LABELLED_LINE or bool(_SIZE_RX.match(n))


# 単位 (Twiss$Unit)。値は OpticsPlot.n の単位の名前
_U_METER = {'EX', 'EY', 'R2', 'DX', 'DY', 'S', 'PEX', 'PEY', 'DZ', 'GX', 'GY', 'GZ', 'OGX', 'OGY', 'OGZ'}
_U_RATIO = {'DPX', 'DPY', 'EPX', 'EPY', 'PEPX', 'PEPY', 'DDP', 'DBX', 'DBY', 'DBZ', 'DGMX', 'DGMY', 'DGMZ'}
_U_INVM = {'R3', 'ZPX', 'ZPY', 'PZPX', 'PZPY', 'GMX', 'GMY', 'GMZ', 'BZS'}
_U_SQRT = {'BX', 'BY', 'BZ', 'RBX', 'RBY', 'RBZ'}
_U_RAD = {'GCHI1', 'GCHI2', 'GCHI3', 'OCHI1', 'OCHI2', 'OCHI3'}
_U_PI2 = {'NX', 'NY', 'NZ'}


def _with_rd(s):
    out = set(s)
    for n in s:
        out.add('R' + n)
        out.add('D' + n)
    return out


_U_METER = _with_rd(_U_METER - {'GX', 'GY', 'GZ', 'OGX', 'OGY', 'OGZ', 'S'}) | {'GX', 'GY', 'GZ', 'OGX', 'OGY', 'OGZ', 'S'}
_U_RATIO = _with_rd(_U_RATIO - {'DDP', 'DBX', 'DBY', 'DBZ', 'DGMX', 'DGMY', 'DGMZ'}) | {'DDP', 'DBX', 'DBY', 'DBZ', 'DGMX', 'DGMY', 'DGMZ'}
_U_INVM = _with_rd(_U_INVM - {'GMX', 'GMY', 'GMZ', 'BZS'}) | {'GMX', 'GMY', 'GMZ', 'RGMX', 'RGMY', 'RGMZ', 'BZS'}
_U_PI2 = _with_rd(_U_PI2)


def _size_unit(name: str) -> str:
    """MakeSize$Label$Unit の単位 (SIG・SIZE の後ろの文字数で決まる)"""
    t = _SIZE_RX.match(name).group(2)
    if len(t) == 1:
        return 'Meter'
    if len(t) == 2:
        if t[0] == 'D':
            return 'Ratio'
        if t[0] == 'P':
            return 'Radian'
        return 'SquareMeter'
    if len(t) == 3:
        return 'Meter'
    return '1'


def twiss_unit(name: str) -> str:
    """Twiss$Unit[name] (OpticsPlot.n)。単位なしは '1'"""
    n = name.upper()
    if _SIZE_RX.match(n):
        return _size_unit(n)
    if n in _U_SQRT:
        return 'SqrtMeter'
    if n in _U_METER:
        return 'Meter'
    if n in _U_RATIO:
        return 'Ratio'
    if n in _U_INVM:
        return 'InvMeter'
    if n in _U_RAD:
        return 'Radian'
    if n == 'BZT':
        return 'Tesla'
    if n in _U_PI2:
        return 'Pi2'
    return '1'


_POSITIVE = {'BX', 'BY', 'BZ', 'GMX', 'GMY', 'GMZ', 'SIGX', 'SIGY', 'SIGZ', 'SIGXX', 'SIGYY', 'SIGZZ',
             'SIZEX', 'SIZEY', 'SIZEZ', 'SIZEXX', 'SIZEYY', 'SIZEZZ'}


def _optics_unit(unit: str, am: float) -> str:
    """Optics$Unit[unit, am]: 振幅 am に合わせて単位を選ぶ"""
    if unit == 'Meter':
        if am < 1e-15:
            return 'Attometer'
        if am < 1e-12:
            return 'Femtometer'
        if am < 1e-9:
            return 'Picometer'
        if am < 1e-6:
            return 'Nanometer'
        if am < 1e-3:
            return 'Micron'
        if am < 1:
            return 'Milimeter'
        if am > 1000:
            return 'Kilometer'
        return 'Meter'
    if unit == 'SquareMeter':
        if am < 1e-9:
            return 'SquareMicrometer'
        if am < 1e-3:
            return 'SquareMilimeter'
        return 'SquareMeter'
    if unit == 'InvMeter':
        return 'InvKilometer' if am < 1e-3 else 'InvMeter'
    if unit == 'Radian':
        if am < 1e-3:
            return 'Microradian'
        if am < 1:
            return 'Miliradian'
        return 'Radian'
    if unit == 'Ratio':
        if am < 1e-7:
            return 'PPB'
        if am < 1e-4:
            return 'PPM'
        if am < 1:
            return 'Percent'
        return 'Ratio'
    return unit


_UNIT_FACTOR = {'Kilometer': 1e3, 'Milimeter': 1e-3, 'Micron': 1e-6, 'Nanometer': 1e-9, 'Picometer': 1e-12,
                'Femtometer': 1e-15, 'Attometer': 1e-18, 'InvKilometer': 1e-3, 'SquareMilimeter': 1e-6,
                'SquareMicrometer': 1e-12, 'Pi2': 2 * math.pi, 'Miliradian': 1e-3, 'Microradian': 1e-6,
                'PPB': 1e-9, 'PPM': 1e-6, 'Percent': 1e-2}

# 単位の文字 (mathtext)。Twiss$Unit$Label の $HaveLaTeX 版に対応
_UNIT_LABEL = {'Meter': 'm', 'Kilometer': 'km', 'Milimeter': 'mm', 'Miliradian': 'mrad', 'Radian': 'rad',
               'Nanometer': 'nm', 'Picometer': 'pm', 'Femtometer': 'fm', 'Attometer': 'am', 'Ratio': '',
               'PPM': 'ppm', 'PPB': 'ppb', 'Percent': '%', 'Micron': r'$\mu$m', 'Microradian': r'$\mu$rad',
               'SqrtMeter': r'$\sqrt{\mathrm{m}}$', 'SquareMeter': r'm$^2$', 'SquareMilimeter': r'mm$^2$',
               'SquareMicrometer': r'$\mu$m$^2$', 'InvMeter': r'm$^{-1}$', 'InvKilometer': r'km$^{-1}$',
               'Pi2': r'$2\pi$', 'Tesla': 'T'}

# 関数のラベル (Twiss$Label の $HaveLaTeX 版を mathtext に)。{r} は参照光学の添字 R、{s} は `_{R}` (R が無ければ空)
_LABEL = {
    'AX': r'\alpha_{{x{r}}}', 'AY': r'\alpha_{{y{r}}}', 'AZ': r'\alpha_{{z{r}}}',
    'GMX': r'\gamma_{{x{r}}}', 'GMY': r'\gamma_{{y{r}}}', 'GMZ': r'\gamma_{{z{r}}}',
    'BX': r'\sqrt{{\beta_{{x{r}}}}}', 'BY': r'\sqrt{{\beta_{{y{r}}}}}', 'BZ': r'\sqrt{{\beta_{{z{r}}}}}',
    'NX': r'\nu_{{x{r}}}', 'NY': r'\nu_{{y{r}}}', 'NZ': r'\nu_{{z{r}}}',
    'EX': r'\eta_{{x{r}}}', 'EY': r'\eta_{{y{r}}}',
    'EPX': r'\eta_{{px{r}}}', 'EPY': r'\eta_{{py{r}}}',
    'PEX': r'\eta_{{x\mathrm{{P}}{r}}}', 'PEY': r'\eta_{{y\mathrm{{P}}{r}}}',
    'PEPX': r'\eta_{{px\mathrm{{P}}{r}}}', 'PEPY': r'\eta_{{py\mathrm{{P}}{r}}}',
    'DX': r'\Delta x{s}', 'DY': r'\Delta y{s}', 'DZ': r'\Delta z{s}',
    'DPX': r'\Delta p_{{x{r}}}', 'DPY': r'\Delta p_{{y{r}}}',
    'DDP': r'\Delta p{s}/p_0',
    'R1': r'R_{{1{r}}}', 'R2': r'R_{{2{r}}}', 'R3': r'R_{{3{r}}}', 'R4': r'R_{{4{r}}}',
    'DETR': r'\det R{s}',
    'ZX': r'\zeta_{{x{r}}}', 'ZY': r'\zeta_{{y{r}}}',
    'ZPX': r'\zeta_{{px{r}}}', 'ZPY': r'\zeta_{{py{r}}}',
    'PZX': r'\zeta_{{x\mathrm{{P}}{r}}}', 'PZY': r'\zeta_{{y\mathrm{{P}}{r}}}',
    'PZPX': r'\zeta_{{px\mathrm{{P}}{r}}}', 'PZPY': r'\zeta_{{py\mathrm{{P}}{r}}}',
}
_LABEL_LINE = {'BMAGX': r'B_{\mathrm{MAG}x}', 'BMAGY': r'B_{\mathrm{MAG}y}', 'BMAGZ': r'B_{\mathrm{MAG}z}',
               'BZS': r'B_z/B\rho', 'BZT': r'B_z', 'GX': 'G_x', 'GY': 'G_y', 'GZ': 'G_z',
               'OGX': r'G_{x,\mathrm{orb}}', 'OGY': r'G_{y,\mathrm{orb}}', 'OGZ': r'G_{z,\mathrm{orb}}',
               'GCHI1': r'\chi_1', 'GCHI2': r'\chi_2', 'GCHI3': r'\chi_3',
               'OCHI1': r'\chi_{1,\mathrm{orb}}', 'OCHI2': r'\chi_{2,\mathrm{orb}}', 'OCHI3': r'\chi_{3,\mathrm{orb}}'}


def twiss_label(name: str) -> str:
    """Twiss$Label[name] (mathtext。`$…$` で囲んだ文字列)。未定義の名前は名前そのもの"""
    n = name.upper()
    if n in _LABEL_LINE:
        return '$' + _LABEL_LINE[n] + '$'
    m = _SIZE_RX.match(n)
    if m:
        t = m.group(2).lower().replace('dp', r'\delta')
        return r'$\sigma_{' + t + '}$' if m.group(1) == 'SIG' else r'$\sigma_{' + t + '}$'
    if n in _LABEL:
        return '$' + _LABEL[n].format(r='', s='') + '$'
    if n[:1] == 'R' and n[1:] in _LABEL:
        return '$' + _LABEL[n[1:]].format(r=r'\mathrm{R}', s=r'_{\mathrm{R}}') + '$'
    if n[:1] == 'D' and n[1:] in _LABEL:
        b = n[1:]
        if b in ('BX', 'BY', 'BZ'):
            return r'$(\delta\beta/\beta)_' + b[1].lower() + '$'
        return r'$\delta ' + _LABEL[b].format(r='', s='') + '$'
    return n


# ---------------------------------------------------------------------------------------------------------------------
# 2. DRAW の構文解析 (CANVASDRAW.n の移植)
# ---------------------------------------------------------------------------------------------------------------------
def separate_ampsnd(w: str) -> Tuple[str, str]:
    """CANVASDRAW.n の SeparateAmpsnd: 先頭から最初の & か ; の手前までと、残り (& か ; から) に分ける。
    & か ; が先頭なら {その 1 文字, 残り}、無ければ {w, ""}"""
    ps = [p for p in (w.find('&'), w.find(';')) if p >= 0]
    if not ps:
        return w, ''
    p = min(ps)
    if p == 0:
        return w[0], w[1:]
    return w[:p], w[p:]


class DrawSpec:
    """`DRAW` の解析結果。r1・r2 は位置 (整数)。panels は図ごとの関数名のリスト、lattice・names はラティス図と名前の表示。
    consumed はコマンドの語 (DRAW の次から数えて) のうち使った数、leftover は DRAW が読み過ぎて FFS に返す語 (無ければ '')"""

    def __init__(self):
        self.r1 = 1
        self.r2 = None
        self.panels: List[List[str]] = []
        self.lattice = False
        self.names: Optional[str] = None
        self.consumed = 0
        self.leftover = ''

    def __repr__(self):
        return (f'DrawSpec(r1={self.r1}, r2={self.r2}, panels={self.panels}, lattice={self.lattice}, names={self.names!r}, '
                f'consumed={self.consumed}, leftover={self.leftover!r})')


class _Reader:
    """`Read[-1, Word]`: 文の語を 1 つずつ返す。文の終わり (SAD では ; か行末) は ''"""

    def __init__(self, words, pos):
        self.words = words
        self.pos = pos

    def read(self) -> str:
        if self.pos < len(self.words):
            w = self.words[self.pos]
            self.pos += 1
            return w
        return ''


def _ffs_of(lat):
    """FFS の状態 (twiss_value を持つもの) を返す。(ffs, I) の組、`.ffs` を持つものも受ける"""
    if hasattr(lat, 'twiss_value'):
        return lat
    if isinstance(lat, tuple) and lat and hasattr(lat[0], 'twiss_value'):
        return lat[0]
    f = getattr(lat, 'ffs', None)
    if f is not None and hasattr(f, 'twiss_value'):
        return f
    raise TypeError('draw_mpl: FFS の状態 (run_sad_string が返す ffs) を渡してください')


def _int_or_position(ffs, w: str):
    """`IntegerQ[Check[r=ToExpression[w],]] || (r=LINE["POSITION",ToUpperCase[w]]) <=> {}`。(成否, 位置) を返す。
    整数の語はそのまま (負は末尾から、0 は 1 に。OpticsPlot が LINE["POSITION",pr]/.{0->1} とする)。名前は最初に一致した位置"""
    if re.fullmatch(r'[+-]?\d+', w):
        return True, int(w)
    try:
        v = ffs.I.run(w, 'script')
        if isinstance(v, int) and not isinstance(v, bool):
            return True, int(v)
    except Exception:
        pass
    try:
        pos = ffs.positions(w.upper())
    except Exception:
        pos = []
    if pos:
        return True, int(pos[0])
    return False, None


def parse_draw(ffs, words: Sequence[str], start: int = 0) -> DrawSpec:
    """CANVASDRAW.n の `Draw[]` の構文解析。words[start:] は `DRAW` の次の語から (先頭が DRAW でも読み飛ばす)。

    [r1 [r2]] 関数… [& 関数…]… [LAT|LATTICE|要素名]
      - r1: 整数か要素名 (LINE["POSITION"])。読めなければ r1=1 でその語は関数として読む。r2 も同様 (失敗は LINE["LENGTH"])
      - 関数: Twiss$Label が空でない語。& で図を分ける。& や ; は語に付いていても分ける (SeparateAmpsnd)
      - 関数でない語が来たら関数の並びは終わり。LAT・LATTICE ならラティス図、要素名 (パターン可) ならラティス図 + その名前
      - それ以外の語は DRAW が読まずに FFS に返す (leftover)"""
    ffs = _ffs_of(ffs)
    spec = DrawSpec()
    nlat = int(ffs.nlat)
    spec.r2 = nlat
    rd = _Reader(list(words), start)
    w = rd.read()
    if w.upper() == 'DRAW':
        w = rd.read()
    if w != '':
        ok, v = _int_or_position(ffs, w)
        if ok:
            spec.r1 = v
            w = rd.read()
            if w != '':
                ok2, v2 = _int_or_position(ffs, w)
                if ok2:
                    spec.r2 = v2
                    w = rd.read()
                else:
                    spec.r2 = nlat
        else:
            spec.r1 = 1
    w = w.upper()
    w1 = '&'
    panels: List[List[str]] = []
    while w1 == '&':
        fig: List[str] = []
        while True:
            w1, w = separate_ampsnd(w)
            if w1 == '':
                w = rd.read().upper()
                if w == '':
                    w1 = ';'
                    break
                continue
            if w1 in ('&', ';'):
                break
            w1u = w1.upper()
            if not is_label(w1u):
                w = w1
                break
            fig.append(w1u)
        panels.append(fig)
    spec.panels = panels
    if w1 != ';' and w != '':
        if w in ('LAT', 'LATTICE'):
            spec.lattice = True
            w1 = ';'
        else:
            try:
                hit = ffs.positions(w.upper())
            except Exception:
                hit = []
            if hit:
                spec.lattice = True
                spec.names = w
                w1 = ';'
    spec.leftover = '' if w1 == ';' else w
    spec.consumed = rd.pos - start - (1 if spec.leftover else 0)
    return spec


# ---------------------------------------------------------------------------------------------------------------------
# 3. データ (OpticsPlot.n の OpticsList・OpticsPlot1$・LatticePlot の移植)
# ---------------------------------------------------------------------------------------------------------------------
_DLMIN = 1e-6
_RESOL = 2048
_PLOT_DIVISION = 128                       # PkgConfig$DefaultList["OpticsPlot"]["PlotDivision"]
_ELEMENT_HEIGHT = {1: 0.08, 2: 0.4, 4: 0.8, 6: 0.6, 8: 0.6, 10: 0.6, 12: 0.6, 21: 0.4, 31: 0.3, 32: 0.3}


def _element_height(ffs, typ: int) -> float:
    if typ == 22:
        return 0.5 if ffs.flags.get('TRPT', True) else 0.8       # If[?TRPT,0.5,0.8,0.8]
    return _ELEMENT_HEIGHT.get(typ, 0.08)


def _element_color(typ: int, v: float):
    if typ == 4:
        return 'limegreen'
    if typ == 2:
        return 'tomato' if v == 0 else 'gold'
    if typ == 1:
        return None
    return 'skyblue'


class PanelData:
    """1 つの図 (枠) の中身。series は関数ごとの (x のリスト, y のリスト) (x は s×ss、y は単位 unit で割った値)"""

    def __init__(self):
        self.funcs: List[str] = []
        self.series: List[Tuple[List[float], List[float]]] = []
        self.labels: List[str] = []
        self.unit = '1'
        self.unit_label = ''
        self.ylabel = ''
        self.ylim: Optional[Tuple[float, float]] = None
        self.dropped: List[str] = []              # 単位が違って描かない関数 (OpticsPlot::units)


class LatticeData:
    """ラティス図。rects: (x0, y0, x1, y1, 塗りの色 (無しは None), 要素名)、names: (表示位置 x (s×ss), 名前)"""

    def __init__(self):
        self.rects: List[Tuple[float, float, float, float, Optional[str], str]] = []
        self.names: List[Tuple[float, str]] = []


class PlotData:
    def __init__(self):
        self.ss = 1.0
        self.slabel = 'm'
        self.pr: Tuple[int, int] = (1, 1)
        self.xlim: Tuple[float, float] = (0.0, 1.0)
        self.panels: List[PanelData] = []
        self.lattice: Optional[LatticeData] = None
        self.title = ''


@contextlib.contextmanager
def _cached_rows(ffs):
    """row_frac (要素の中の Twiss。入口の閉軌道から伝搬する計算) を (k, fr) ごとに保存して、関数ごとの再計算を避ける"""
    orig = ffs.row_frac
    cache: Dict[tuple, object] = {}
    had = 'row_frac' in getattr(ffs, '__dict__', {})

    def rf(k, fr, rows=None):
        if rows is not None:
            return orig(k, fr, rows)
        key = (k, fr)
        r = cache.get(key)
        if r is None:
            r = cache[key] = orig(k, fr)
        return r
    ffs.row_frac = rf
    try:
        yield
    finally:
        if not had:
            try:
                del ffs.row_frac
            except AttributeError:
                pass
        else:
            ffs.row_frac = orig


def _line(ffs, key, k, fr=0.0):
    return ffs.line_value(key, k, fr)


def _function_value(ffs, f: str, kind: str, v: float):
    """OpticsList の fun[i + j/ndiv]: Twiss[f, v] (v は実数の位置)。v の整数部と小数部は Twiss[] の位置の解釈と同じ"""
    k = int(math.floor(v))
    fr = v - k
    if k > ffs.nlat:
        k = ((k - 1) % ffs.nlat) + 1
    if kind == 'line':
        return float(ffs.line_value(f, k, fr))
    x = float(ffs.twiss_value(f, k, fr))
    if kind == 'sqrt':
        return math.sqrt(x) if x >= 0 else float('nan')
    return x


def _function_kind(f: str):
    """(種類, 要素の中を分割するか)。OpticsList の Which。種類: sqrt (√Twiss)・twiss・line"""
    if f in ('BX', 'BY', 'BZ'):
        return 'sqrt', True
    if f in TWISS_FUNCTIONS:
        return 'twiss', True
    if f in ('RBX', 'RBY', 'RBZ'):
        return 'sqrt', False
    if f in XTWISS_FUNCTIONS:
        return 'twiss', False
    if f in LINE_FUNCTIONS or _SIZE_RX.match(f):
        return 'line' if f in LINE_FUNCTIONS else 'twiss', True
    raise ValueError(f'OpticsPlot::wronglist: {f}')


def optics_list(ffs, f: str, i1: int, i2: int, dl: float, off: float):
    """OpticsList[f, i1, i2, dl, off, …]: 位置 i1..i2 の各要素を ndiv 個に割り、(S + L j/ndiv + off, f(i + j/ndiv)) を返す"""
    kind, div = _function_kind(f)
    xs: List[float] = []
    ys: List[float] = []
    for i in range(i1, i2 + 1):
        if not div:
            nd = 1
        elif str(_line(ffs, 'TYPENAME', i)) == 'MARK':
            nd = 0
        else:
            nd = int(math.floor(float(_line(ffs, 'L', i)) / dl)) + 1
        if nd <= 0:
            continue
        s_i = float(_line(ffs, 'S', i))
        l_i = float(_line(ffs, 'L', i))
        for j in range(nd):
            xs.append(s_i + l_i * j / nd + off)
            ys.append(_function_value(ffs, f, kind, i + j / nd))
    return xs, ys


def resolve_region(ffs, r1, r2) -> Tuple[int, int]:
    """OpticsPlot の `pr=LINE["POSITION",pr]/.{0->1,{}->LINE["LENGTH"]}`。整数・要素名・None を位置の整数にする"""
    nlat = int(ffs.nlat)
    out = []
    for r, dflt in ((r1, 1), (r2, nlat)):
        if r is None:
            out.append(dflt)
        elif isinstance(r, str):
            pos = ffs.positions(r.upper())
            out.append(int(pos[0]) if pos else nlat)
        else:
            k = int(r)
            if k == 0:
                k = 1
            elif k < 0:
                k = nlat + k + 1
            while k > nlat:
                k -= nlat
            out.append(max(1, k))
    return out[0], out[1]


def _optics_series(ffs, f: str, pr1: int, pr2: int):
    """OpticsList1$: 領域 pr1..pr2 (pr2<=pr1 は周回をまたぐ) の点列。返り値 (xs, ys, x0, leng)。x は s (m)。x0 は領域の始点の s"""
    nlat = int(ffs.nlat)
    s = lambda k: float(_line(ffs, 'S', k))
    ndiv = _PLOT_DIVISION
    if pr2 <= pr1:
        ndiv = min(_RESOL, max(ndiv, nlat + pr2 - pr1))
        leng = s(nlat) - s(1) + s(pr2) - s(pr1)
        dl = max(_DLMIN, leng / ndiv)
        x0 = s(pr1) - s(nlat) + s(1)
        xa, ya = optics_list(ffs, f, pr1, nlat, dl, s(1) - s(nlat))
        xb, yb = optics_list(ffs, f, 1, pr2, dl, 0.0)
        return xa + xb, ya + yb, x0, leng
    ndiv = min(_RESOL, max(ndiv, pr2 - pr1))
    leng = s(pr2) - s(pr1)
    dl = max(_DLMIN, leng / ndiv)
    xs, ys = optics_list(ffs, f, pr1, pr2, dl, 0.0)
    return xs, ys, s(pr1), leng


def _unit_label_text(unit: str) -> str:
    return _UNIT_LABEL.get(unit, '')


def _build_panel(ffs, funcs: Sequence[str], pr: Tuple[int, int], ss: float, warn=None, ylim=None) -> Tuple[PanelData, float, float]:
    """OpticsPlot1$: 1 枠の関数のリストから点列・単位・軸の範囲・ラベルを作る。返り値 (PanelData, x0, leng)"""
    pd = PanelData()
    if not funcs:
        return pd, 0.0, 0.0
    unit0 = twiss_unit(funcs[0])
    data = []
    x0 = leng = 0.0
    for f in funcs:
        xs, ys, x0_, leng_ = _optics_series(ffs, f, pr[0], pr[1])
        if twiss_unit(f) != unit0:                                   # Cases[l,{{_,unit,_,_,_},_,_},1]: 先頭と同じ単位だけ描く
            pd.dropped.append(f)
            if warn:
                warn(f'OpticsPlot::units: Different units in one frame: {twiss_unit(f)} {unit0} ({f} は描かない)')
            continue
        data.append((f, xs, ys))
        x0, leng = x0_, leng_
    allv = [y for _, _, ys in data for y in ys if y == y]
    mn, mx = (min(allv), max(allv)) if allv else (0.0, 0.0)
    am = max(abs(mn), abs(mx))
    unit = _optics_unit(unit0, am)
    ru = _UNIT_FACTOR.get(unit, 1.0)
    for f, xs, ys in data:
        pd.funcs.append(f)
        pd.series.append(([x * ss for x in xs], [(y / ru) if ru != 1.0 else y for y in ys]))
        pd.labels.append(twiss_label(f))
    pd.unit = unit
    pd.unit_label = _unit_label_text(unit)
    if ylim is not None:
        pd.ylim = (ylim[0] / ru, ylim[1] / ru)
    elif all(f in _POSITIVE for f in pd.funcs):
        pd.ylim = (0.0, am / ru * 1.1)
    lab = ', '.join(pd.labels)
    pd.ylabel = f'{lab} ({pd.unit_label})' if pd.unit_label else lab
    return pd, x0, leng


def _wild_match(pat: str, name: str) -> bool:
    """StringMatchQ[name, pat] (* と ? のワイルドカード、大文字小文字は区別しない)"""
    rx = '^' + ''.join('.*' if c == '*' else ('.' if c in '?%' else re.escape(c)) for c in pat) + '$'
    return re.match(rx, name, re.I) is not None


def lattice_data(ffs, pr: Tuple[int, int], ss: float, names: Optional[str] = '*') -> LatticeData:
    """LatticePlot (LatticePlot1$ と LatticeNames)。領域 pr (周回をまたぐなら 2 つに分ける)。
    names: 名前を出す要素のパターン (StringMatchQ。既定 "*" = 全部)。None・False なら名前なし"""
    ld = LatticeData()
    nlat = int(ffs.nlat)
    s = lambda k: float(_line(ffs, 'S', k))
    pr1, pr2 = pr
    if pr2 <= pr1:
        ltotal = s(nlat) - s(1)
        segs = [((pr1, nlat), -ltotal), ((1, pr2), 0.0)]
        s0_ref = s(pr1) - ltotal
    else:
        segs = [((pr1, pr2), 0.0)]
        s0_ref = s(pr1)
    leng = (s(pr2) - s0_ref)
    for (a, b), off in segs:
        _lattice_rects(ffs, ld, a, b, off, ss)
    if isinstance(names, str) and names != '':
        for (a, b), off in segs:
            _lattice_names(ffs, ld, a, b, off, leng / 100.0, names, ss)
    return ld


def _lattice_rects(ffs, ld: LatticeData, pr1: int, pr2: int, off: float, ss: float):
    """LatticePlot1$: 要素の箱。薄い要素は (s の幅が 1/400 以下で、前に描いた箱より低ければ) 描かない"""
    s1 = float(_line(ffs, 'S', pr1))
    ds = (float(_line(ffs, 'S', pr2 - 1)) - s1) / 400.0 if pr2 - 1 >= pr1 else 0.0
    h1 = None
    for i in range(pr1, pr2):
        typ = int(float(_line(ffs, 'TYPE', i)))
        h = _element_height(ffs, typ)
        si = float(_line(ffs, 'S', i))
        si1 = float(_line(ffs, 'S', i + 1))
        if typ == 4 or typ == 22 or si1 > s1 + ds or (h1 is not None and h > h1):
            el = ffs.els[i - 1] if i - 1 < len(ffs.els) else None
            v = 0.0
            if el is not None:
                dk = ffs.default_key(el)
                if dk:
                    try:
                        v = float(_line(ffs, dk, i))
                    except Exception:
                        v = 0.0
            h1, s1 = h, si1
            yoff = (math.copysign(1.0, v) if v != 0 else 0.0) * 0.1 if h > 0.1 * 1.1 else 0.0
            ld.rects.append(((si + off) * ss, yoff - h / 2, (si1 + off) * ss, yoff + h / 2, _element_color(typ, v),
                             str(_line(ffs, 'NAME', i))))


def _lattice_names(ffs, ld: LatticeData, pr1: int, pr2: int, off: float, w: float, pattern: str, ss: float):
    """LatticeNames: 要素名 (LINE["NAME"] が pattern に合うもの)。直前に出した名前と同じ要素 (LINE["ELEMENT"]) は出さない"""
    s0 = float(_line(ffs, 'S', pr1)) + off
    e = ''
    for i in range(pr1, pr2):
        nm = str(_line(ffs, 'NAME', i))
        el = str(_line(ffs, 'ELEMENT', i))
        smid = (float(_line(ffs, 'S', i)) + float(_line(ffs, 'S', i + 1))) / 2 + off
        if _wild_match(pattern, nm) and e != el and smid >= s0:
            s0 = max(smid - w / 2, s0)
            e = el
            ld.names.append((s0 * ss, nm))
        else:
            e = ''


def build_plot_data(lat, panels: Sequence[Sequence[str]], region=(None, None), lattice: bool = False,
                    names: Optional[str] = '*', title: str = '', warn=None, ylim=None) -> PlotData:
    """OpticsPlot[panels, Region->region, Lattice->lattice, Names->names] が描く点・箱・名前を作る (描画はしない)。
    CALC していなければ先に光学を計算する (pysad の Twiss[] と同じ。gs の OpticsPlot は計算していない光学の 0 を描く)"""
    ffs = _ffs_of(lat)
    if getattr(ffs, 'lat', None) is None:
        raise RuntimeError('DRAW: ビームラインがありません (USE が先)')
    ffs.ensure_rows()
    pr = resolve_region(ffs, region[0], region[1])
    pd = PlotData()
    pd.pr = pr
    length = abs(float(_line(ffs, 'S', pr[1])) - float(_line(ffs, 'S', pr[0])))
    pd.ss, pd.slabel = (0.001, 'km') if length >= 10000 else (1.0, 'm')
    pd.title = title or ''
    with _cached_rows(ffs):
        x0 = None
        for funcs in panels:
            p, x0_, leng_ = _build_panel(ffs, list(funcs), pr, pd.ss, warn=warn,
                                         ylim=(ylim if ylim is not None and not isinstance(ylim, dict) else (ylim or {}).get(len(pd.panels))))
            pd.panels.append(p)
            if funcs:
                x0, leng = x0_, leng_
        if x0 is None:                                             # 関数が 1 つも無いときも領域は決める
            s = lambda k: float(_line(ffs, 'S', k))
            if pr[1] <= pr[0]:
                x0, leng = s(pr[0]) - s(ffs.nlat) + s(1), s(ffs.nlat) - s(1) + s(pr[1]) - s(pr[0])
            else:
                x0, leng = s(pr[0]), s(pr[1]) - s(pr[0])
        pd.xlim = ((x0 - 0.001) * pd.ss, (x0 + leng + 0.001) * pd.ss)
        if lattice:
            pd.lattice = lattice_data(ffs, pr, pd.ss, names)
    return pd


# ---------------------------------------------------------------------------------------------------------------------
# 3b. InfoLabel (OP$InfoLabel: 位置をクリックしたときの要素と光学関数の表)
# ---------------------------------------------------------------------------------------------------------------------
def _f10(x) -> str:
    """$FORM="10.5" 相当の数の表示 (有効数字 5 桁)"""
    try:
        return f'{float(x):.5g}'
    except Exception:
        return str(x)


def element_info(lat, s: float) -> str:
    """OP$InfoLabel のポップアップの文字: s (m) の位置を含む要素の入口の名前・位置・種類の値・光学関数 (βx,y・αx,y・νx,y・ηx,y・η′・R)。
    OpticsPlot.n は R3,4 の行の 2 つ目に R2 を出している (R4 の誤記と思われる)。ここでは R4 を出す"""
    ffs = _ffs_of(lat)
    nlat = int(ffs.nlat)
    circ = float(_line(ffs, 'S', nlat)) - float(_line(ffs, 'S', 1))
    if s < float(_line(ffs, 'S', 1)):
        s = s + circ
    idx = nlat
    for k in range(1, nlat + 1):
        if float(_line(ffs, 'S', k)) >= s:
            idx = k
            break
    p = max(1, idx - 1)
    tt = str(_line(ffs, 'TYPENAME', p))
    te = ''
    names = {'QUAD': ('K1', 'K1 = {} m^-1'), 'SEXT': ('K2', 'K2 = {} m^-2'), 'OCT': ('K3', 'K3 = {} m^-3'),
             'DECA': ('K4', 'K4 = {} m^-4'), 'DODECA': ('K5', 'K5 = {} m^-5'), 'CAVI': ('VOLT', 'Vacc = {} V'),
             'TCAVI': ('K0', 'K0 = {} rad')}
    if tt == 'BEND':
        if float(_line(ffs, 'ANGLE', p)) == 0:
            tt = 'STEER'
            te = f'K0 = {_f10(_line(ffs, "K0", p))} rad'
        else:
            te = f'angle = {_f10(_line(ffs, "ANGLE", p))} rad'
    elif tt in names:
        key, fmt = names[tt]
        te = fmt.format(_f10(_line(ffs, key, p)))
    tw = lambda k: ffs.twiss_value(k, p)
    lines = [f'Entrance of {_line(ffs, "NAME", p)}:   s = {_f10(float(_line(ffs, "S", p)))} m']
    if te:
        lines.append(f'{tt:<6}: {te}')
    lines += [f'beta_x,y = {_f10(tw("BX"))}, {_f10(tw("BY"))} m',
              f'alpha_x,y = {_f10(tw("AX"))}, {_f10(tw("AY"))}',
              f'nu_x,y = {_f10(tw("NX") / 2 / math.pi)}, {_f10(tw("NY") / 2 / math.pi)}',
              f'eta_x,y = {_f10(tw("EX"))}, {_f10(tw("EY"))} m',
              f"eta'_x,y = {_f10(tw('EPX'))}, {_f10(tw('EPY'))}",
              f'R1,2 = {_f10(tw("R1"))}, {_f10(tw("R2"))} m',
              f'R3,4 = {_f10(tw("R3"))}, {_f10(tw("R4"))} m^-1']
    return '\n'.join(lines)


def tune_text(lat) -> str:
    """図に添える情報: 周回のチューン νx, νy (Twiss["NX","$$$"]/2π)・周長"""
    ffs = _ffs_of(lat)
    n = int(ffs.nlat)
    nx = ffs.twiss_value('NX', n) / 2 / math.pi
    ny = ffs.twiss_value('NY', n) / 2 / math.pi
    return f'$\\nu_x$ = {nx:.4f}   $\\nu_y$ = {ny:.4f}   C = {float(_line(ffs, "S", n)) - float(_line(ffs, "S", 1)):.4g} m'


# ---------------------------------------------------------------------------------------------------------------------
# 4. matplotlib での描画
# ---------------------------------------------------------------------------------------------------------------------
_COLORS = ['royalblue', 'tomato', 'forestgreen', 'darkgoldenrod', 'dodgerblue', 'darkorange', 'blueviolet',
           'limegreen', 'hotpink', 'olive']                       # Plot$DefaultColor (tomato3・orange3・yellow4 は近い色に)
_DASHES = [(0, ()), (0, (6, 1.5)), (0, (3, 1)), (0, (1.5, 1)), (0, (1, 0.8)), (0, (6, 0.8, 0.8, 0.8)), (0, (3, 0.8, 0.8, 0.8))]
_NONINTERACTIVE = {'agg', 'pdf', 'ps', 'svg', 'cairo', 'template', 'pgf'}


def _import_mpl():
    """matplotlib を遅延 import。無ければ None"""
    try:
        import matplotlib
        import matplotlib.pyplot as plt          # noqa: F401
        return matplotlib
    except Exception:
        return None


def _interactive_backend(matplotlib) -> bool:
    b = str(matplotlib.get_backend()).lower()
    return b not in _NONINTERACTIVE and 'inline' not in b


def render(data: PlotData, lat=None, info: bool = True, grid=(False, True), figsize=None, tunes: bool = False,
           names_overlap: bool = True):
    """PlotData を matplotlib の Figure にする (pyplot の Figure。表示も保存もしない)。
    枠は縦に積み、s の軸を共有する。ラティス図は一番下 (枠の高さの 0.23 倍。OpticsPlot.n の LatticeRegion {-0.15, 0.08})。
    info=True なら、図をクリックするとその位置の要素と光学関数を出す (OP$InfoLabel)。lat はそのときと tunes に使う"""
    mpl = _import_mpl()
    if mpl is None:
        raise ImportError('matplotlib がありません (pip install matplotlib)')
    import matplotlib.pyplot as plt
    nf = max(1, len(data.panels))
    has_lat = data.lattice is not None
    ratios = [1.0 / nf] * nf + ([0.23] if has_lat else [])
    nrow = len(ratios)
    if figsize is None:
        figsize = (9.0, max(3.2, 2.1 * nf + (1.5 if has_lat else 0.6)))
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrow, 1, height_ratios=ratios, hspace=0.22)
    axes = []
    for r in range(nrow):
        axes.append(fig.add_subplot(gs[r, 0], sharex=axes[0] if axes else None))
    dax = axes[:nf]
    for ax, p in zip(dax, data.panels):
        for n, (xs, ys) in enumerate(p.series):
            ax.plot(xs, ys, color=_COLORS[n % len(_COLORS)], linestyle=_DASHES[n % len(_DASHES)][1] and (0, _DASHES[n % len(_DASHES)][1]) or '-',
                    linewidth=1.2, label=p.labels[n])
        ax.set_ylabel(p.ylabel)
        if p.ylim is not None:
            ax.set_ylim(*p.ylim)
        if len(p.series) > 1:
            ax.legend(loc='upper right', fontsize=8, framealpha=0.8, ncol=len(p.series))
        for axis_, on_ in (('x', grid[0]), ('y', grid[1])):
            if on_:
                ax.grid(axis=axis_, alpha=0.35)
        ax.tick_params(labelbottom=False)
    ax_last = dax[-1]
    lat_ax = axes[-1] if has_lat else None
    if has_lat:
        _draw_lattice(lat_ax, data, names_overlap)
        lat_ax.tick_params(labelbottom=False, bottom=False)
    ax_last.tick_params(labelbottom=True)
    ax_last.set_xlabel(f's ({data.slabel})')
    axes[0].set_xlim(*data.xlim)
    if data.title:
        fig.suptitle(data.title, fontsize=11)
    if tunes and lat is not None:
        try:
            fig.text(0.99, 0.005, tune_text(lat), ha='right', va='bottom', fontsize=8)
        except Exception:
            pass
    nm = max((len(n) for _, n in data.lattice.names), default=0) if has_lat else 0
    bottom = 0.10 + min(0.22, 0.010 * nm) if has_lat else 0.12
    fig.subplots_adjust(left=0.10, right=0.97, top=0.93 if data.title else 0.97, bottom=bottom)
    fig._pysad_draw = {'data': data, 'axes': axes}
    if info and lat is not None:
        _connect_info(fig, axes, lat)
    return fig


def _draw_lattice(ax, data: PlotData, names_overlap: bool):
    from matplotlib.patches import Rectangle
    ld = data.lattice
    for x0, y0, x1, y1, color, name in ld.rects:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=color if color else 'none', edgecolor='black',
                               linewidth=0.5, gid=name))
    ax.set_ylim(-0.6, 0.6)
    ax.set_yticks([])
    for sp in ('left', 'right', 'top'):
        ax.spines[sp].set_visible(False)
    shown = ld.names
    if names_overlap:
        # RemoveOverlap: 文字 (縦書き、高さ = フォントの大きさ) が重なる名前は消す。s の軸の長さと図の幅から重なりを決める
        fs = 6.0
        width_pt = ax.figure.get_figwidth() * 0.87 * 72.0
        span = max(data.xlim[1] - data.xlim[0], 1e-12)
        keep = []
        last = -1e30
        for x, n in shown:
            if (x - last) / span * width_pt >= fs * 1.15:
                keep.append((x, n))
                last = x
        shown = keep
    ld_names = []
    for x, n in shown:
        ld_names.append(ax.text(x, -0.65, n, rotation=270, ha='center', va='top', fontsize=6.0, clip_on=False))
    ax._pysad_names = ld_names
    ax.set_xlim(*data.xlim)


def _connect_info(fig, axes, lat):
    """InfoLabel->True: 図をクリックするとその s の要素・光学関数を注釈に出す (もう一度クリックすると位置を変える。右クリックで消す)"""
    state = {'ann': None}

    def on_click(ev):
        if ev.inaxes is None or ev.xdata is None:
            return
        ax = ev.inaxes
        if state['ann'] is not None:
            try:
                state['ann'].remove()
            except Exception:
                pass
            state['ann'] = None
        if getattr(ev, 'button', 1) == 3:
            fig.canvas.draw_idle()
            return
        data = fig._pysad_draw['data']
        try:
            txt = element_info(lat, ev.xdata / data.ss)
        except Exception as e:                                      # noqa: BLE001
            txt = f'(情報なし: {e})'
        state['ann'] = ax.annotate(txt, (ev.xdata, ev.ydata), xytext=(12, 12), textcoords='offset points', fontsize=7,
                                   family='monospace', bbox=dict(boxstyle='round', fc='lightyellow', ec='gray', alpha=0.95),
                                   zorder=10)
        fig.canvas.draw_idle()
    fig._pysad_info_cid = fig.canvas.mpl_connect('button_press_event', on_click)


# ---------------------------------------------------------------------------------------------------------------------
# 5. 入口: draw (Python API)・run_draw (FFS の DRAW)
# ---------------------------------------------------------------------------------------------------------------------
_saved_paths: Dict[str, int] = {}
_atexit_registered = [False]


def _next_path(path: str) -> str:
    """同じ保存先に 2 回目以降の DRAW は `名前_2.拡張子` のように番号を付ける (上書きしない)。{n} があれば番号に置き換える"""
    if '{n}' in path:
        n = _saved_paths.get(path, 0) + 1
        _saved_paths[path] = n
        return path.replace('{n}', str(n))
    n = _saved_paths.get(path, 0) + 1
    _saved_paths[path] = n
    if n == 1:
        return path
    root, ext = os.path.splitext(path)
    return f'{root}_{n}{ext}'


def _save(fig, path: str, dpi: Optional[float] = None):
    ext = os.path.splitext(path)[1].lower().lstrip('.')
    if ext not in ('png', 'pdf', 'svg', 'jpg', 'jpeg', 'eps', 'ps'):
        raise ValueError(f'DRAW: 保存の拡張子は png・pdf・svg ({path})')
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    fig.savefig(path, dpi=dpi or float(os.environ.get('PYSAD_DRAW_DPI', '90')))


def _show(fig, mpl):
    import matplotlib.pyplot as plt
    plt.show(block=False)
    try:
        plt.pause(0.001)
    except Exception:
        pass
    if not _atexit_registered[0] and os.environ.get('PYSAD_DRAW_BLOCK', '1') != '0':
        _atexit_registered[0] = True
        atexit.register(lambda: plt.get_fignums() and plt.show())


def draw(lat, spec='', region=None, save: Optional[str] = None, show: Optional[bool] = None, lattice: Optional[bool] = None,
         names=None, title: Optional[str] = None, grid=(False, True), info: bool = True, tunes: bool = False,
         figsize=None, dpi: Optional[float] = None, warn=None):
    """`DRAW` を Python から呼ぶ。Figure を返す (matplotlib が無ければ ImportError)。

    lat:     FFS の状態 (`run_sad_string` が返す ffs)。USE・CALC のあとに呼ぶ (CALC していなければ計算する)
    spec:    DRAW の引数と同じ構文の文字列 `"[r1 [r2]] BX BY & EX [LAT|名前]"` (リストでも可。先頭の DRAW は付けても付けなくてもよい)
    region:  (r1, r2) 位置の整数か要素名 (None は spec の r1 r2、無ければ 1 と LINE["LENGTH"])
    lattice: None は spec の LAT・名前による。True・False で上書き
    names:   ラティス図の名前のパターン (None は spec か、LAT なら "*" = 全部。False で名前なし)
    save:    ファイル名 (png・pdf・svg)。None は環境変数 PYSAD_DRAW_FILE を使わず保存しない (DRAW コマンドは環境変数を見る)
    show:    None は保存しなければ対話の backend のときだけ plt.show。True・False で上書き
    title:   PlotLabel (None は SAD の TITLE。ここでは何も出さない)
    tunes:   True でチューンと周長の文字を図の右下に出す (SAD の OpticsPlot にはない追加)"""
    ffs = _ffs_of(lat)
    words = spec.split() if isinstance(spec, str) else list(spec)
    ps = parse_draw(ffs, words, 0)
    r1, r2 = (ps.r1, ps.r2) if region is None else region
    use_lat = ps.lattice if lattice is None else bool(lattice)
    nm = names if names is not None else (ps.names if ps.names else '*')
    if nm is False:
        nm = None
    if not ps.panels or all(not p for p in ps.panels):
        raise ValueError('DRAW: 描く関数がありません (BX BY EX … を書く)')
    mpl = _import_mpl()
    if mpl is None:
        raise ImportError('matplotlib がありません (pip install matplotlib)')
    pdat = build_plot_data(ffs, ps.panels, (r1, r2), use_lat, nm, title or '', warn=warn)
    fig = render(pdat, ffs, info=info, grid=grid, figsize=figsize, tunes=tunes)
    if save:
        _save(fig, save, dpi)
    if show is None:
        show = (not save) and _interactive_backend(mpl)
    if show:
        _show(fig, mpl)
    return fig


def _draw_options(ffs) -> Dict[str, object]:
    """SAD の `Draw$Option` (Rule のリスト) を {名前 (大文字): 値} にする。先に書いた Rule が優先 (Override の規則)"""
    try:
        lst = ffs.I.globals.get('Draw$Option')
    except Exception:
        return {}
    out: Dict[str, object] = {}
    if not isinstance(lst, list):
        return out
    for r in lst:
        if hasattr(r, 'lhs') and hasattr(r, 'rhs'):
            k = getattr(r.lhs, 'name', r.lhs)
            k = str(k).upper()
            if k not in out:
                out[k] = r.rhs
    return out


def _sym_name(v) -> str:
    return str(getattr(v, 'name', v))


def _display_mode(ffs) -> str:
    """'off' (K64 と同じく何もしない)・'mpl'。`$DisplayFunction` が Identity なら off、CanvasDrawer・MatplotlibDrawer なら mpl。
    環境変数 PYSAD_DISPLAY=mpl は $DisplayFunction に依らず mpl"""
    if os.environ.get('PYSAD_DISPLAY', '').lower() in ('mpl', 'matplotlib'):
        return 'mpl'
    try:
        df = ffs.I.globals.get('$DisplayFunction')
    except Exception:
        df = None
    n = _sym_name(df) if df is not None else 'Identity'
    return 'off' if n == 'Identity' else 'mpl'


def _title_of(ffs) -> str:
    try:
        t = ffs.I.globals.get('TITLE')
    except Exception:
        return ''
    return t if isinstance(t, str) else ''


def run_draw(ffs, words: List[str], i: int) -> int:
    """FFS の `DRAW` の本体 (ffscmd.py から呼ぶ)。words[i] が DRAW。使った語の数を返す (DRAW が読まない語は残して次のコマンドに回す)。

    CANVASDRAW.n と同じ: 領域・関数・LAT・要素名を読み、Draw$Option (SAD の Rule のリスト) を最優先にして OpticsPlot 相当を描く。
    出力: 保存先 (環境変数 PYSAD_DRAW_FILE か Draw$Option の `File->"x.png"`) があれば保存 (表示しない)。無ければ対話の backend なら plt.show。
    `$DisplayFunction` が Identity (K64 では DRAW が何もしない) で保存先も無ければ何もしない"""
    ps = parse_draw(ffs, words, i + 1)
    used = 1 + ps.consumed
    if ps.leftover and ps.leftover.upper() != words[i + used].upper():
        words[i + used] = ps.leftover                       # SeparateAmpsnd で割った語の前半を次のコマンドに返す
    opts = _draw_options(ffs)
    file_ = opts.get('FILE')
    file_ = str(file_) if isinstance(file_, str) else os.environ.get('PYSAD_DRAW_FILE', '')
    if not file_ and _display_mode(ffs) == 'off':
        return used
    if not ps.panels or all(not p for p in ps.panels):
        ffs.write('[pysad] DRAW: 描く関数がありません (BX BY EX … を書く)\n')
        return used
    if _import_mpl() is None:
        ffs.write('[pysad] DRAW: matplotlib がありません。何も描きません (pip install matplotlib)\n')
        return used
    if getattr(ffs, 'lat', None) is None:
        ffs.write(' ???-DRAW-Missing beamline in USE.\n')
        return used
    # Draw$Option (最優先)・DRAW の引数・既定の順
    region = opts.get('REGION')
    if isinstance(region, list) and len(region) == 2:
        r1, r2 = region
    else:
        r1, r2 = ps.r1, ps.r2
    lat_opt = opts.get('LATTICE')
    lattice = ps.lattice if lat_opt is None else (_sym_name(lat_opt) == 'True')
    names = opts.get('NAMES', ps.names if ps.names else '*')
    if not isinstance(names, str):
        names = None
    grid = (False, True)
    gl = opts.get('GRIDLINES')
    if isinstance(gl, list) and len(gl) == 2:
        grid = tuple(_sym_name(g) != 'None' for g in gl)
    info = _sym_name(opts.get('INFOLABEL', 'True')) != 'False'
    title = opts.get('PLOTLABEL')
    title = title if isinstance(title, str) else _title_of(ffs)
    try:
        data = build_plot_data(ffs, ps.panels, (r1, r2), lattice, names, title, warn=lambda m: ffs.write(f'[pysad] {m}\n'))
        fig = render(data, ffs, info=info, grid=grid, tunes=_sym_name(opts.get('TUNES', 'False')) == 'True')
    except Exception as e:                                  # noqa: BLE001
        ffs.write(f'[pysad] DRAW: 描けません: {e}\n')
        return used
    import matplotlib.pyplot as plt
    mpl = _import_mpl()
    if file_:
        path = _next_path(file_)
        try:
            _save(fig, path)
            ffs.write(f'[pysad] DRAW: {path} に保存しました\n')
        except Exception as e:                              # noqa: BLE001
            ffs.write(f'[pysad] DRAW: 保存できません: {e}\n')
        if _sym_name(opts.get('SHOW', 'False')) == 'True' and _interactive_backend(mpl):
            _show(fig, mpl)
        else:
            plt.close(fig)
    elif _interactive_backend(mpl):
        _show(fig, mpl)
    else:
        plt.close(fig)                                      # Agg など: 何も出さない
    return used


# ---------------------------------------------------------------------------------------------------------------------
# 6. 照合の道具: OpticsPlot が返す Graphics (pysad の Canvas.n の経路) から点列・箱・名前を取り出す
# ---------------------------------------------------------------------------------------------------------------------
def _walk(e):
    from .sadscript import Expr
    if isinstance(e, Expr):
        yield e
        for a in e.args:
            yield from _walk(a)
    elif isinstance(e, (list, tuple)):
        for a in e:
            yield from _walk(a)


def _headname(e) -> str:
    return str(getattr(e.head, 'name', e.head))


def graphics_lines(g) -> List[List[Tuple[float, float]]]:
    """Graphics の Line[{{x,y},…}] の点列を、現れる順に全部返す (OpticsPlot では 枠 1 → 枠 2 … の関数の順)"""
    out = []
    for e in _walk(g):
        if _headname(e) == 'Line' and e.args and isinstance(e.args[0], list):
            pts = e.args[0]
            if pts and isinstance(pts[0], list) and len(pts[0]) == 2:
                out.append([(float(p[0]), float(p[1])) for p in pts])
    return out


def graphics_rectangles(g) -> List[Tuple[float, float, float, float]]:
    """ラティス図の Rectangle[{x0,y0},{x1,y1}, 選択肢…] (中に Graphics を持たない箱)。(x0, y0, x1, y1) を返す"""
    out = []
    for e in _walk(g):
        if _headname(e) == 'Rectangle' and len(e.args) >= 2 and all(isinstance(a, list) and len(a) == 2 for a in e.args[:2]):
            if any(_headname(x) == 'Graphics' for a in e.args[2:] for x in _walk(a)):
                continue
            (a0, b0), (a1, b1) = e.args[0], e.args[1]
            out.append((float(a0), float(b0), float(a1), float(b1)))
    return out


def graphics_texts(g) -> List[Tuple[str, float]]:
    """ラティス図の名前 Text[{name,{Scaled[x],Scaled[y]}}, …]。(名前, x) を返す"""
    out = []
    for e in _walk(g):
        if _headname(e) == 'Text' and e.args and isinstance(e.args[0], list) and len(e.args[0]) == 2:
            nm, pos = e.args[0]
            if isinstance(nm, str) and isinstance(pos, list) and pos and hasattr(pos[0], 'args') and pos[0].args:
                out.append((nm, float(pos[0].args[0])))
    return out
