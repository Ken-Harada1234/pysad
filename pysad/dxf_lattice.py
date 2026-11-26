"""Write a plan view of the lattice as a text DXF (ASCII, AutoCAD R12 equivalent) file.
dxf_lattice.py — ラティスの平面図 (上から見た図) をテキストの DXF (ASCII、AutoCAD R12 相当) で書き出す (v2.4.1)

K. Oide の `drawmagvw3d.n` (SAD の 3GSR/esrf。CAD ソフト Vectorworks の VectorScript を書き出す `DrawMagnets3D`) の
幾何の考え方 (LINE["GEO"]・LINE["TYPE"]・LINE["ANGLE"]・LINE["E1"/"E2"]・LINE["DIR"] から偏向磁石・四極などの外形を作る)
を参考にした自前版。既定は平面図 (2D)。`mode='3d'` で磁石の立体の外形 (3DFACE の直方体・楔形) も書ける (下の「立体」)。
同じ計算を SADScript で書いた `example/dxf_lattice/DrawMagnetsDXF.n` (gs でも pysad でも走る) の出力と一致する (試験 `test_dxf_lattice.py`)。

使い方 (USE 済みの FFS を渡す):
    from pysad.dxf_lattice import load_ffs, write_dxf
    ffs = load_ffs('example/00_lattice/FODO_ring.sad', 'RING')
    info = write_dxf(ffs, 'ring.dxf', region=None, scale=1.0)

外形 (すべて局所座標 (u, v) = (入口の進行方向, 左) で作ってから、入口の GEO (gx, gy, chi1) で実験室の (X, Y) に写す。
beam の x は -v なので、ANGLE>0 の偏向は左に曲がる):
  * 偏向磁石 (BEND。ANGLE≠0 の MULT も): 連続する同じ名前の要素を 1 台の磁石にまとめ、入口・出口の GEO から弦の長さ lc0 と偏向角を求め、
    drawmagvw3d.n の DrawBend と同じ式の 4 角形 (辺は弦に平行、端面は E1・E2 の傾き、軌道は辺の中央から ±hs ずれる) を作る。
    磁石の幅は Width、PoleWidth、Type (Htype: 左右対称、Ctype1・Ctype2: PoleWidth/2 が軌道の外側で残りが内側の C 型)。
    PhysicalLength があれば物理長を使い、弦の中央に置く。
  * QUAD・SEXT・OCT・CAVI など: 長さ × 幅 (Width) の矩形。SOL は長さ 0 の境界なので、BZ≠0 の SOL から次の SOL までを矩形 (レイヤー SOL) にする
  * DRIFT・その他: 基準軌道の線だけ。MARK: 軌道に直角な短い線 (MarkSize)
  * 基準軌道: ANGLE≠0 の要素は GEO の分数位置 (LINE["GEO", i+f]) で ArcStep (rad) 以下の折れ線に近似、直線部は端点だけ
レイヤーと色番号: BEND 1 (赤)・QUAD 5 (青)・SEXT 3 (緑)・OCT 6 (マゼンタ)・MULT 30 (橙)・SOL 4 (シアン)・CAVI 2 (黄)・
  MARK 8 (灰)・TRAJECTORY 7 (白/黒)・TEXT 7。
単位は m (INSUNITS=6)。scale=1000 なら mm (INSUNITS=4)、100 なら cm (5)。座標は 6 桁 (scale 後) に丸める。
磁石の寸法の取り方 (優先順は 1 つの量ごと): `spec` の要素名 > `spec` の型名 > `MagnetSpec[MagnetType[要素名]]` (drawmagvw3d.n や draw_cad.sad が使う定義。
USE した台本が定義していれば自動で読む。`magnet_spec=False` で無視) > 既定値。
要素名の文字の重なり: `overlap='remove'` で、文字高と位置から見積もった外接矩形が重なる名前を、偏向磁石 > 四極 > 6 極 > 8 極 > その他 > MARK の優先で間引く
(既定 'keep' は従来の出力のまま)。
立体 (mode='3d'。既定の 'plan' の出力は 1 行も変わらない):
  各磁石を、2D の外形の 4 角形 (偏向磁石は入口・出口の面が E1・E2 に傾いた楔形) を z = -Height/2 … +Height/2 に押し出した 6 面体 (3DFACE 6 枚)
  にする。底面の頂点の (X, Y) は 2D の外形と同じ値。Height は `spec` / `MagnetSpec` の Height (なければ BEND 0.406 m (drawmagvw3d.n の既定)、
  他は幅 Width)。基準軌道の面 (z = 0) が磁石の中心。レイヤー・色は 2D と同じ規則。`parts='split'` で偏向磁石 (ANGLE≠0 の MULT も) と四極を
  ヨークとポールに分ける (ポールはレイヤー POLE、色 9)。コイル・端板・ベローズは未対応。
限界は docs/22_lattice_dxf.md。
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

# SAD の要素型コード (LINE["TYPE"]) → (レイヤー, 既定の幅 m)。drawmagvw3d.n の既定は BEND 0.63・QUAD 0.45。他は本実装の既定
_MAGNET_TYPES = {2: ('BEND', 0.63), 4: ('QUAD', 0.45), 6: ('SEXT', 0.30), 8: ('OCT', 0.30), 10: ('MULT', 0.30),
                 12: ('MULT', 0.30), 22: ('MULT', 0.30), 20: ('SOL', 0.40), 31: ('CAVI', 0.40), 32: ('CAVI', 0.40)}
_TYPE_NAME = {2: 'BEND', 4: 'QUAD', 6: 'SEXT', 8: 'OCT', 10: 'DECA', 12: 'DODECA', 22: 'MULT', 20: 'SOL', 31: 'CAVI', 32: 'TCAV'}
_BEND_POLEWIDTH = 0.13           # drawmagvw3d.n の既定
_CHAR_W = 0.6                    # 要素名の外接矩形の見積もり: 幅 = 文字数 × 文字高 × 0.6、高さ = 文字高 (フォントに依らない定数。gs と pysad で同じ結果にする)
_TEXT_PRIO_MARK = 6              # 重なり消しの優先順位: 1 偏向磁石 (ANGLE≠0 の MULT も)、2 四極、3 6 極、4 8 極、5 その他の磁石、6 MARK (数が小さいほど残す)
_TEXT_PRIO = {2: 1, 4: 2, 6: 3, 8: 4}
LAYERS = [('0', 7), ('BEND', 1), ('QUAD', 5), ('SEXT', 3), ('OCT', 6), ('MULT', 30), ('SOL', 4), ('CAVI', 2),
          ('MARK', 8), ('TRAJECTORY', 7), ('TEXT', 7)]
_COMMENT = 'plan view of the lattice (DrawMagnetsDXF); units: see $INSUNITS'
_COMMENT_3D = '3D magnet solids of the lattice (DrawMagnetsDXF Mode 3D); units: see $INSUNITS'
LAYERS_SPLIT = LAYERS + [('POLE', 9)]       # 3D の parts='split' のときだけ。ポールのレイヤー
_BEND_HEIGHT = 0.406             # 3D の偏向磁石の既定の高さ (drawmagvw3d.n の既定)。他の磁石の既定の高さは幅 Width
_BEND_GAP = 0.04                 # 3D split: 偏向磁石の磁極の間隔の既定 (drawmagvw3d.n の Gap の既定)


def _fmt(x: float, scale: float) -> str:
    v = x * scale
    if abs(v) < 5e-7:
        v = 0.0
    return '%.6f' % v


class _Slot:
    """overlap='remove' で、要素名の TEXT を後で入れる場所 (k は候補の番号)"""
    __slots__ = ('k',)

    def __init__(self, k):
        self.k = k


class _Writer:
    def __init__(self, scale: float):
        self.scale = scale
        self.out: List[str] = []
        self.n = {'POLYLINE': 0, 'TEXT': 0, 'LINE': 0}

    def num(self, x: float) -> str:
        return _fmt(x, self.scale)

    def polyline(self, layer: str, pts: Sequence[Tuple[float, float]], closed: bool):
        o = self.out
        o += ['  0', 'POLYLINE', '  8', layer, ' 66', '1', ' 10', '0.000000', ' 20', '0.000000', ' 30', '0.000000',
              ' 70', '1' if closed else '0']
        for (x, y) in pts:
            o += ['  0', 'VERTEX', '  8', layer, ' 10', self.num(x), ' 20', self.num(y), ' 30', '0.000000']
        o += ['  0', 'SEQEND', '  8', layer]
        self.n['POLYLINE'] += 1

    def line(self, layer: str, p0, p1):
        o = self.out
        o += ['  0', 'LINE', '  8', layer, ' 10', self.num(p0[0]), ' 20', self.num(p0[1]), ' 30', '0.000000',
              ' 11', self.num(p1[0]), ' 21', self.num(p1[1]), ' 31', '0.000000']
        self.n['LINE'] += 1

    def face(self, layer: str, vs):
        """3DFACE 1 枚。vs = 4 頂点の (x, y, z) の文字列 (num 済み)"""
        o = self.out
        o += ['  0', '3DFACE', '  8', layer]
        for k, (x, y, z) in enumerate(vs):
            o += [' 1%d' % k, x, ' 2%d' % k, y, ' 3%d' % k, z]
        self.n['3DFACE'] = self.n.get('3DFACE', 0) + 1

    def solid(self, layer: str, q, z0: float, z1: float):
        """4 角形 q ((X, Y)、上から見て反時計回り) を z0 〜 z1 に押し出した 6 面体 (3DFACE 6 枚、外向きの法線)。頂点は 8 個を一度だけ数にする"""
        top = [(self.num(x), self.num(y), self.num(z1)) for (x, y) in q]
        bot = [(self.num(x), self.num(y), self.num(z0)) for (x, y) in q]
        self.face(layer, (top[0], top[1], top[2], top[3]))
        self.face(layer, (bot[3], bot[2], bot[1], bot[0]))
        for i in range(4):
            j = (i + 1) % 4
            self.face(layer, (bot[i], bot[j], top[j], top[i]))

    def text(self, layer: str, p, height: float, rot_deg: float, s: str):
        self.out += self.text_lines(layer, p, height, rot_deg, s)
        self.n['TEXT'] += 1

    def text_lines(self, layer: str, p, height: float, rot_deg: float, s: str):
        x, y = self.num(p[0]), self.num(p[1])
        return ['  0', 'TEXT', '  8', layer, ' 10', x, ' 20', y, ' 30', '0.000000', ' 40', self.num(height), '  1', s,
              ' 50', '%.4f' % (0.0 if abs(rot_deg) < 5e-5 else rot_deg), ' 72', '1', ' 11', x, ' 21', y, ' 31', '0.000000', ' 73', '2']


def _insunits(scale: float) -> str:
    return {1.0: '6', 1000.0: '4', 100.0: '5'}.get(float(scale), '0')


def _header(scale: float, layers=None, comment=None) -> List[str]:
    layers = LAYERS if layers is None else layers
    h = ['999', _COMMENT if comment is None else comment, '  0', 'SECTION', '  2', 'HEADER', '  9', '$ACADVER', '  1', 'AC1009',
         '  9', '$INSUNITS', ' 70', _insunits(scale), '  0', 'ENDSEC',
         '  0', 'SECTION', '  2', 'TABLES',
         '  0', 'TABLE', '  2', 'LTYPE', ' 70', '1',
         '  0', 'LTYPE', '  2', 'CONTINUOUS', ' 70', '0', '  3', 'Solid line', ' 72', '65', ' 73', '0', ' 40', '0.0',
         '  0', 'ENDTAB',
         '  0', 'TABLE', '  2', 'LAYER', ' 70', str(len(layers))]
    for name, col in layers:
        h += ['  0', 'LAYER', '  2', name, ' 70', '0', ' 62', str(col), '  6', 'CONTINUOUS']
    h += ['  0', 'ENDTAB',
          '  0', 'TABLE', '  2', 'STYLE', ' 70', '1',
          '  0', 'STYLE', '  2', 'STANDARD', ' 70', '0', ' 40', '0.0', ' 41', '1.0', ' 50', '0.0', ' 71', '0', ' 42', '0.2',
          '  3', 'txt', '  4', '',
          '  0', 'ENDTAB',
          '  0', 'ENDSEC',
          '  0', 'SECTION', '  2', 'ENTITIES']
    return h


_FOOTER = ['  0', 'ENDSEC', '  0', 'EOF']


# ---------------------------------------------------------------- ラティスの読み込み
def load_ffs(src: str, use: Optional[str] = None, cwd: Optional[str] = None):
    """Read a SAD script (file name or text), apply USE and return the FFS object (no CALC).
    SAD の台本 (ファイル名または本文) を読み、USE して FFS オブジェクトを返す (CALC はしない。幾何だけなら不要)。
    src が存在するファイル名ならその内容、そうでなければ本文とみなす。`FFS;` の行が無ければ補う。use は LINE 名。"""
    import io
    import re
    from .ffscmd import run_sad_string
    is_file = len(src) < 400 and os.path.isfile(src)
    text = open(src, encoding='utf-8').read() if is_file else src
    if cwd is None and is_file:
        cwd = os.path.dirname(os.path.abspath(src))
    add = '\n'
    if re.search(r'(?mi)^\s*FFS\s*;', text) is None:
        add += 'FFS;\n'
    if use:
        add += 'USE %s;\n' % use
    ffs, _ = run_sad_string(text + add, out=io.StringIO(), cwd=cwd)
    return ffs


# ---------------------------------------------------------------- 幾何
def _wrap_pi(a: float) -> float:
    while a > math.pi:
        a -= 2.0 * math.pi
    while a <= -math.pi:
        a += 2.0 * math.pi
    return a


def _region_index(ffs, x, nlat: int, first: bool) -> int:
    if x is None:
        return 1 if first else nlat - 1
    if isinstance(x, str):
        v = ffs.fn_line(['POSITION', x])
        if isinstance(v, list):
            if not v:
                raise ValueError('region の位置が見つからない: %r' % x)
            v = v[0] if first else v[-1]
        return int(v)
    return int(x)


class _Plan:
    """LINE[...] を引く部品 (SAD 台本と同じキーを使う)。幾何は ffs.line_value('GEO', k, fr) = LINE["GEO", k+fr]"""

    def __init__(self, ffs):
        self.f = ffs
        self.nlat = int(ffs.nlat)

    def lv(self, key, k, fr=0.0):
        return self.f.line_value(key, k, fr)

    def geo(self, k, fr=0.0):
        g = self.lv('GEO', k, fr)
        return g[0][0], g[0][1], g[1][0]          # gx, gy, chi1


def _frame(g):
    """入口の幾何 (gx, gy, chi1) から (u, v) → (X, Y) の写像"""
    gx, gy, c1 = g
    cu, su = math.cos(c1), math.sin(c1)

    def f(u, v):
        return (gx + u * cu - v * su, gy + u * su + v * cu)
    return f, c1


def _bend_polygon(spec, angle, lc0, e1, e2, polewidth, width, typ):
    """drawmagvw3d.n DrawBend の式 (局所の (u, v))。angle = 偏向角 (符号つき)、lc0 = 入口から出口までの弦、e1・e2 = 端面の角 (rad)。
    戻り (4 角形の頂点, 物理長の弦 lc, 弦方向のずらし dshift)。辺は弦に平行、wl = 左の半幅、wr = 右の半幅。
    出口の r3・r4 は辺が弦に平行になるよう入口と同じ wl・wr (drawmagvw3d.n は wl・wr の入れ替わりがある。非対称の磁石で違う)。"""
    phys = spec.get('PhysicalLength')
    if phys is None:
        la = lc0 if angle == 0 else angle / 2.0 / math.sin(angle / 2.0) * lc0
    else:
        la = float(phys)
    lc = la if angle == 0 else 2.0 * math.sin(angle / 2.0) / angle * la
    hs = 0.0 if abs(angle) < 1e-10 else la / angle * math.sin(angle / 4.0) ** 2
    if typ == 'Ctype1':
        if angle >= 0:
            wr = polewidth / 2.0
            wl = width - wr
        else:
            wl = polewidth / 2.0
            wr = width - wl
    elif typ == 'Ctype2':
        if angle < 0:
            wr = polewidth / 2.0
            wl = width - wr
        else:
            wl = polewidth / 2.0
            wr = width - wl
    else:
        wr = wl = width / 2.0
    ah = angle / 2.0
    r1 = (wl - hs) / math.cos(e1 - ah)
    r2 = (wr + hs) / math.cos(e1 - ah)
    r3 = (wr + hs) / math.cos(e2 - ah)
    r4 = (wl - hs) / math.cos(e2 - ah)
    p1 = (-r1 * math.sin(e1), r1 * math.cos(e1))
    p2 = (r2 * math.sin(e1), -r2 * math.cos(e1))
    pe = (math.cos(ah) * lc, math.sin(ah) * lc)
    p3 = (pe[0] + r3 * math.sin(angle - e2), pe[1] - r3 * math.cos(angle - e2))
    p4 = (pe[0] - r4 * math.sin(angle - e2), pe[1] + r4 * math.cos(angle - e2))
    # 物理長の弦を、軌道の弦の中央に置く (弦方向に (lc0-lc)/2 ずらす)
    ds = (lc0 - lc) / 2.0
    c, s = math.cos(ah), math.sin(ah)
    pts = [(p[0] + ds * c, p[1] + ds * s) for p in (p1, p2, p3, p4)]
    return pts, lc, ds


def _bend_strip3d(angle, lc0, phys, e1, e2, vr, vl):
    """3D 用。偏向磁石の楔形の、横位置 vr … vl (入口の軌道からの左右の距離。左が正。vr < vl) の帯の 4 角形 (局所の (u, v))。
    頂点の順は `_bend_polygon` と同じ (左入口・右入口・右出口・左出口)。vr = -wr、vl = wl なら `_bend_polygon` の出力と同じ。
    弦に平行な辺の横位置は v - hs、入口・出口の面は E1・E2 の傾きの直線 (入口の面は軌道の入口の点を通る)"""
    la = (lc0 if angle == 0 else angle / 2.0 / math.sin(angle / 2.0) * lc0) if phys is None else float(phys)
    lc = la if angle == 0 else 2.0 * math.sin(angle / 2.0) / angle * la
    hs = 0.0 if abs(angle) < 1e-10 else la / angle * math.sin(angle / 4.0) ** 2
    ah = angle / 2.0
    pe = (math.cos(ah) * lc, math.sin(ah) * lc)

    def ent(v):
        r = (v - hs) / math.cos(e1 - ah)
        return (-r * math.sin(e1), r * math.cos(e1))

    def ext(v):
        r = (v - hs) / math.cos(e2 - ah)
        return (pe[0] - r * math.sin(angle - e2), pe[1] + r * math.cos(angle - e2))
    ds = (lc0 - lc) / 2.0
    c, s = math.cos(ah), math.sin(ah)
    return [(p[0] + ds * c, p[1] + ds * s) for p in (ent(vl), ent(vr), ext(vr), ext(vl))]


def _box_strip3d(ds, lc, vr, vl):
    """3D 用。直線の磁石 (QUAD など) の、横位置 vr … vl の帯の 4 角形 (局所の (u, v))。頂点の順は `_bend_strip3d` と同じ"""
    return [(ds, vl), (ds, vr), (ds + lc, vr), (ds + lc, vl)]


def _bend_blocks(wl, wr, polew, height, gap, pole_h, leg_w):
    """3D の parts='split': 偏向磁石 (H 型・C 型) を、ヨークの上下の板・側脚と上下のポールの直方体に分ける。
    戻り [(部品 'yoke'|'pole', vr, vl, z0, z1)]。z は軌道の面が 0。ポールは軌道の左右 PoleWidth/2 の幅、間隔 gap の上下に厚み pole_h。
    ヨークの板はその外側の厚み (height - gap)/2 - pole_h。側脚は半幅が PoleWidth/2 より大きい側 (H 型は両側、C 型は広い側) に幅 leg_w で、
    上下の板をつなぐ。寸法が入らないときは ValueError"""
    zi = gap / 2.0 + pole_h
    if not (gap > 0 and pole_h > 0 and zi < height / 2.0 - 1e-12 and polew > 0 and leg_w > 0):
        raise ValueError('parts="split" の偏向磁石: Gap=%g・PoleHeight=%g・PoleWidth=%g・LegWidth=%g は Height=%g に入らない '
                         '(Gap/2 + PoleHeight < Height/2 が必要)' % (gap, pole_h, polew, leg_w, height))
    if wl < polew / 2.0 - 1e-9 or wr < polew / 2.0 - 1e-9:
        raise ValueError('parts="split" の偏向磁石: PoleWidth=%g が Width に入らない (左右の半幅 %g・%g)' % (polew, wl, wr))
    blocks = [('yoke', -wr, wl, zi, height / 2.0), ('yoke', -wr, wl, -height / 2.0, -zi)]
    for side, hw in ((1, wl), (-1, wr)):
        if hw > polew / 2.0 + 1e-9:
            if hw - leg_w < polew / 2.0 - 1e-9:
                raise ValueError('parts="split" の偏向磁石: LegWidth=%g が大きすぎる (ポールの幅 %g に当たる)' % (leg_w, polew))
            blocks.append(('yoke', hw - leg_w, hw, -zi, zi) if side > 0 else ('yoke', -hw, -hw + leg_w, -zi, zi))
    blocks += [('pole', -polew / 2.0, polew / 2.0, gap / 2.0, zi), ('pole', -polew / 2.0, polew / 2.0, -zi, -gap / 2.0)]
    return blocks


def _quad_blocks(width, height, aperture, leg_w, polew):
    """3D の parts='split': 四極を、ヨークの 4 本の枠 (上・下・左・右) と 4 つのポール (上・下・左・右の向きに内側へ) の直方体に分ける
    (模式図。実機の四極は 45 度の向きの 4 極。ボア = aperture の直径)。戻り [(部品, vr, vl, z0, z1)] (左 = v が正)"""
    m = min(width, height)
    if not (aperture > 0 and leg_w > 0 and polew > 0 and aperture / 2.0 < m / 2.0 - leg_w - 1e-12 and polew <= aperture + 1e-12):
        raise ValueError('parts="split" の四極: Aperture=%g・LegWidth=%g・PoleWidth=%g は Width=%g・Height=%g に入らない '
                         '(Aperture/2 < min(Width, Height)/2 - LegWidth、PoleWidth <= Aperture が必要)' % (aperture, leg_w, polew, width, height))
    w2, h2, a2, zi, vi = width / 2.0, height / 2.0, aperture / 2.0, height / 2.0 - leg_w, width / 2.0 - leg_w
    return [('yoke', -w2, w2, zi, h2), ('yoke', -w2, w2, -h2, -zi), ('yoke', vi, w2, -zi, zi), ('yoke', -w2, -vi, -zi, zi),
            ('pole', -polew / 2.0, polew / 2.0, a2, zi), ('pole', -polew / 2.0, polew / 2.0, -zi, -a2),
            ('pole', a2, vi, -polew / 2.0, polew / 2.0), ('pole', -vi, -a2, -polew / 2.0, polew / 2.0)]


def _text_rect(x, y, h, rot_deg, s):
    """要素名の外接矩形 (中心 (x, y)、TEXT の中央寄せ 72=1・73=2 に合わせる): [x, y, 半幅, 半高, cos, sin]。幅は文字数 × 文字高 × _CHAR_W"""
    r = rot_deg * math.pi / 180.0
    return [x, y, h * _CHAR_W * len(s) / 2.0, h / 2.0, math.cos(r), math.sin(r)]


def _rects_overlap(a, b) -> bool:
    """向きのついた 2 つの矩形が重なるか (分離軸の判定。接するだけは重ならない)。矩形は _text_rect の形"""
    dx, dy = b[0] - a[0], b[1] - a[1]
    for (nx, ny) in ((a[4], a[5]), (-a[5], a[4]), (b[4], b[5]), (-b[5], b[4])):
        ra = a[2] * abs(nx * a[4] + ny * a[5]) + a[3] * abs(-nx * a[5] + ny * a[4])
        rb = b[2] * abs(nx * b[4] + ny * b[5]) + b[3] * abs(-nx * b[5] + ny * b[4])
        if abs(dx * nx + dy * ny) >= ra + rb:
            return False
    return True


def _remove_overlaps(cands) -> List[bool]:
    """要素名の候補 [(x, y, h, rot_deg, s, prio)] (m 単位、描く順) のうち残すものの真偽のリスト。
    優先順位 prio の小さいものから順に (同じ優先なら描く順に)、すでに残したものと重なるなら捨てる。
    近傍の判定は一辺 2 × (最大の外接円の半径) の格子で行う (全部の組を調べない)"""
    n = len(cands)
    rects = [_text_rect(c[0], c[1], c[2], c[3], c[4]) for c in cands]
    keep = [False] * n
    if n == 0:
        return keep
    rmax = max(math.hypot(r[2], r[3]) for r in rects)
    cs = 2.0 * rmax if rmax > 0 else 1.0
    grid: Dict[Tuple[int, int], List[int]] = {}
    for lev in range(1, _TEXT_PRIO_MARK + 1):
        for k in range(n):
            if cands[k][5] != lev:
                continue
            ix, iy = int(math.floor(rects[k][0] / cs)), int(math.floor(rects[k][1] / cs))
            ok = True
            for jx in (ix - 1, ix, ix + 1):
                for jy in (iy - 1, iy, iy + 1):
                    for j in grid.get((jx, jy), ()):
                        if _rects_overlap(rects[k], rects[j]):
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    break
            if ok:
                keep[k] = True
                grid.setdefault((ix, iy), []).append(k)
    return keep


def _interp_magnet_spec(ffs, ename: str, cache: Dict) -> Dict:
    """USE した台本の `MagnetSpec[MagnetType[要素名]]` を pysad のインタプリタで評価して {キー: 値} にする (drawmagvw3d.n と同じ定義)。
    `MagnetType[]` が未定義なら要素名そのもの (文字列、次いで Symbol[要素名]) をキーにする。リスト以外 (未定義) は {}。
    値は数ならそのまま、Symbol (Ctype1 など) は文字列、None は None。規則 (Rule) でない項は無視"""
    if ename in cache:
        return cache[ename]
    res: Dict = {}
    I = getattr(ffs, 'I', None)
    if I is not None and '"' not in ename and '\\' not in ename:
        import io
        old, old_stopped = I.out, getattr(I, 'stopped', False)
        I.out = io.StringIO()
        I.stopped = False                # 台本が Exit・STOP で終わっていても (stopped のままだと run が何もしない) 定義は読める
        try:
            I.run('DXFmsr = Module[{k, r}, k = MagnetType["%s"]; If[Head[k] === MagnetType, k = "%s"]; r = MagnetSpec[k]; '
                  'If[Head[r] === List, r, MagnetSpec[Symbol["%s"]]]];' % (ename, ename, ename), 'ffs')
            v = None
            for k in list(I.globals):
                if k.lower() == 'dxfmsr':
                    v = I.globals.pop(k)
            if isinstance(v, list):
                for r in v:
                    if hasattr(r, 'lhs') and hasattr(r, 'rhs'):
                        val = r.rhs
                        if not isinstance(val, (int, float)):
                            val = None if str(val) == 'None' else str(val)
                        res[str(r.lhs)] = val
        except Exception:      # 評価できない定義は無いものとする
            res = {}
        finally:
            I.out = old
            I.stopped = old_stopped
    cache[ename] = res
    return res


def build_dxf(ffs, region=None, scale: float = 1.0, text_height: float = 0.15, names: bool = True,
              mark_size: float = 0.1, mark_text: bool = False, arc_step: float = 0.05, text_side: float = -1,
              spec: Optional[Dict] = None, overlap: str = 'keep', magnet_spec=None, magnet_type=None,
              mode: str = 'plan', parts: str = 'body') -> Tuple[List[str], Dict]:
    """Return the DXF text lines and a dict of information such as the number of elements; see `write_dxf`.
    DXF の本文 (行のリスト) と、要素数などの情報の dict を返す。引数は write_dxf を見よ"""
    spec = spec or {}
    mode = str(mode).lower()
    if mode in ('2d', 'plan'):
        mode = 'plan'
    if mode not in ('plan', '3d'):
        raise ValueError("mode は 'plan' (2D、'2d') か '3d': %r" % (mode,))
    parts = str(parts).lower()
    if parts not in ('body', 'split'):
        raise ValueError("parts は 'body' か 'split': %r" % (parts,))
    solid3d = (mode == '3d')
    split3d = solid3d and parts == 'split'
    if overlap not in ('keep', 'remove'):
        raise ValueError("overlap は 'keep' か 'remove': %r" % (overlap,))
    ms_cache: Dict = {}

    def mspec_of(ename):
        """要素名 → MagnetSpec の規則の dict。magnet_spec: None=台本の MagnetSpec[] を自動で読む、False=読まない、dict={キー: {規則}}"""
        if magnet_spec is False:
            return {}
        if isinstance(magnet_spec, dict):
            key = magnet_type(ename) if callable(magnet_type) else ename
            d = magnet_spec.get(key)
            if d is None:
                d = magnet_spec.get(ename)
            return d if isinstance(d, dict) else {}
        return _interp_magnet_spec(ffs, ename, ms_cache)
    P = _Plan(ffs)
    nlat = P.nlat
    if region is None:
        region = (None, None)
    li = _region_index(ffs, region[0], nlat, True)
    lf = min(_region_index(ffs, region[1], nlat, False), nlat - 1)
    if li < 1 or li > lf:
        raise ValueError('region が空: %r' % (region,))
    W = _Writer(scale)
    info = {'li': li, 'lf': lf, 'nmagnets': 0, 'nmarks': 0, 'ntraj': 0, 'magnet_spec_used': {}}
    if solid3d:
        info.update({'mode': '3d', 'parts': parts, 'nsolids': 0, 'solids': []})
    cands: List[Tuple] = []            # overlap='remove' の要素名の候補 (x, y, h, rot, 名前, 優先順位)。DXF の中の場所は W.out の _Slot

    # ---- 基準軌道 (LINE["GEO", i + f] の折れ線。直線部は曲がり・COORD の前後だけ頂点を置く)
    pts: List[Tuple[float, float]] = []

    def add(g):
        p = (g[0], g[1])
        if not pts or pts[-1] != p:
            pts.append(p)
    for i in range(li, lf + 1):
        ang_i = P.lv('ANGLE', i)
        ang_p = P.lv('ANGLE', i - 1) if i > li else 0.0
        typ_i = int(P.lv('TYPE', i))
        typ_p = int(P.lv('TYPE', i - 1)) if i > li else 0
        if i == li or ang_i != 0 or ang_p != 0 or typ_i == 35 or typ_p == 35:
            add(P.geo(i))
        if ang_i != 0:
            n = max(1, int(math.ceil(abs(ang_i) / arc_step - 1e-9)))
            for j in range(1, n):
                add(P.geo(i, j / n))
    add(P.geo(lf + 1))
    W.polyline('TRAJECTORY', pts, False)
    info['ntraj'] = len(pts)

    # ---- 磁石・MARK
    def lookup(name_el, typ_code, key, default):
        # 優先順: spec の要素名 > spec の型名 > MagnetSpec[MagnetType[要素名]] > 既定値 (量ごと)
        for k in (name_el, _TYPE_NAME.get(typ_code, '?')):
            d = spec.get(k)
            if isinstance(d, dict) and key in d:
                return d[key]
        ms = mspec_of(name_el)
        if ms.get(key) is not None:
            info['magnet_spec_used'].setdefault(name_el, {})[key] = ms[key]
            return ms[key]
        return default

    def put_text(p, h, rot, s_, prio):
        if overlap == 'remove':
            W.out.append(_Slot(len(cands)))
            cands.append((p[0], p[1], h, rot, s_, prio))
        else:
            W.text('TEXT', p, h, rot, s_)

    def draw_group(is_, ia_, it):
        layer, wdef = _MAGNET_TYPES[it]
        ename = P.lv('ELEMENT', is_)
        g0 = P.geo(is_)
        g1 = P.geo(ia_ + 1)
        lc0 = math.hypot(g1[0] - g0[0], g1[1] - g0[1])
        if lc0 < 1e-12:
            return
        fr, c1 = _frame(g0)
        width = float(lookup(ename, it, 'Width', wdef))
        lookup(ename, it, 'Height', None)            # 平面図では描かない。MagnetSpec の値を info['magnet_spec_used'] に残すだけ
        angle = 0.0
        bendlike = (it == 2) or (it == 22 and P.lv('ANGLE', is_) != 0)
        if bendlike:
            angle = _wrap_pi(g1[2] - g0[2])
        if bendlike:
            polew = float(lookup(ename, it, 'PoleWidth', _BEND_POLEWIDTH))
            typ = lookup(ename, it, 'Type', 'Htype')
            if it == 2:
                ang_s = P.lv('ANGLE', is_)
                e1 = (P.lv('E1', is_) if P.lv('DIR', is_) == 1 else P.lv('E2', is_)) * ang_s
                ang_e = P.lv('ANGLE', ia_)
                e2 = (P.lv('E2', ia_) if P.lv('DIR', ia_) == 1 else P.lv('E1', ia_)) * ang_e
            else:
                e1 = e2 = 0.0
            sp = {'PhysicalLength': lookup(ename, it, 'PhysicalLength', None)}
            poly, lc, ds = _bend_polygon(sp, angle, lc0, e1, e2, polew, width, typ)
            strip = lambda vr, vl: _bend_strip3d(angle, lc0, sp['PhysicalLength'], e1, e2, vr, vl)
        else:
            phys = lookup(ename, it, 'PhysicalLength', None)
            lc = lc0 if phys is None else float(phys)
            ds = (lc0 - lc) / 2.0
            poly = [(ds, width / 2.0), (ds, -width / 2.0), (ds + lc, -width / 2.0), (ds + lc, width / 2.0)]
            strip = lambda vr, vl: _box_strip3d(ds, lc, vr, vl)
        if solid3d:
            hgt = float(lookup(ename, it, 'Height', _BEND_HEIGHT if bendlike else width))
            if hgt <= 0:
                raise ValueError('Height は正: %s の Height = %g' % (ename, hgt))
            blocks = [('body', -width / 2.0, width / 2.0, -hgt / 2.0, hgt / 2.0)]
            if split3d and bendlike:
                if (typ == 'Ctype1' and angle >= 0) or (typ == 'Ctype2' and angle < 0):
                    wr_ = polew / 2.0
                    wl_ = width - wr_
                elif typ in ('Ctype1', 'Ctype2'):
                    wl_ = polew / 2.0
                    wr_ = width - wl_
                else:
                    wl_ = wr_ = width / 2.0
                gap = float(lookup(ename, it, 'Gap', _BEND_GAP))
                pole_h = float(lookup(ename, it, 'PoleHeight', (hgt - gap) / 4.0))
                leg_w = float(lookup(ename, it, 'LegWidth', (width - polew) / 4.0))
                blocks = _bend_blocks(wl_, wr_, polew, hgt, gap, pole_h, leg_w)
            elif split3d and it == 4:
                apert = float(lookup(ename, it, 'Aperture', 0.4 * min(width, hgt)))
                leg_w = float(lookup(ename, it, 'LegWidth', (min(width, hgt) - apert) / 4.0))
                pole_w = float(lookup(ename, it, 'PoleWidth', apert))
                blocks = _quad_blocks(width, hgt, apert, leg_w, pole_w)
            for part, vr, vl, z0, z1 in blocks:
                q = poly if part == 'body' else strip(vr, vl)        # body は 2D の外形そのもの (頂点は 2D と同じ値)
                lay = 'POLE' if part == 'pole' else layer
                W.solid(lay, [fr(u, v) for (u, v) in q], z0, z1)
                info['nsolids'] += 1
                info['solids'].append((ename, part, lay, z1 - z0))
        else:
            xy = [fr(u, v) for (u, v) in poly]
            W.polyline(layer, xy, True)
        info['nmagnets'] += 1
        if names:
            cu = sum(p[0] for p in poly) / 4.0
            cv = sum(p[1] for p in poly) / 4.0
            vt = cv + text_side * (width / 2.0 + text_height) if text_side else cv
            # 向き: 弦の方向 (入口の向き + angle/2)。読める向き (-90, 90] にする
            rot = math.degrees(c1 + angle / 2.0)
            rot = (rot + 180.0) % 360.0 - 180.0
            if rot > 90.0:
                rot -= 180.0
            elif rot <= -90.0:
                rot += 180.0
            put_text(fr(cu, vt), text_height, rot, ename, 1 if bendlike else _TEXT_PRIO.get(it, 5))

    it = 0
    gname = ''
    start = False
    is_ = 0
    sol0 = 0
    for i in range(li, lf + 2):
        last = i > lf
        t = 0 if last else int(P.lv('TYPE', i))
        nm = '' if last else P.lv('ELEMENT', i)
        if start:
            same = (not last) and t == it and nm == gname
            nxt = min(i + 1, nlat)
            same_next = (not last) and i + 1 <= lf and int(P.lv('TYPE', nxt)) == it and P.lv('ELEMENT', nxt) == gname
            stop = (not same) and (last or P.lv('L', i) != 0 or not same_next)
            if stop:
                draw_group(is_, i - 1, it)
                start = False
                it = 0
        if (not last) and t == 20:
            # SOL は長さ 0 の境界 (BOUND=1)。BZ≠0 の SOL から次の SOL までを、ソレノイドの領域の矩形 (レイヤー SOL) にする
            if sol0 > 0:
                draw_group(sol0, i - 1, 20)
                sol0 = 0
            if P.lv('BZ', i) != 0:
                sol0 = i
        if last and sol0 > 0:
            draw_group(sol0, lf, 20)
        if (not start) and (not last) and t in _MAGNET_TYPES and t != 20:
            it, gname, is_, start = t, nm, i, True
        if (not last) and t == 41:
            g = P.geo(i)
            fr, _ = _frame(g)
            h = mark_size / 2.0
            a, b = fr(0.0, h), fr(0.0, -h)
            W.line('MARK', a, b)
            info['nmarks'] += 1
            if mark_text:
                put_text(fr(0.0, h + text_height), text_height, 0.0, nm, _TEXT_PRIO_MARK)

    if overlap == 'remove':
        keep = _remove_overlaps(cands)
        out2: List[str] = []
        for o in W.out:
            if isinstance(o, _Slot):
                c = cands[o.k]
                if keep[o.k]:
                    out2 += W.text_lines('TEXT', (c[0], c[1]), c[2], c[3], c[4])
                    W.n['TEXT'] += 1
            else:
                out2.append(o)
        W.out = out2
        info['ntext_candidates'] = len(cands)
        info['ntext_removed'] = len(cands) - sum(keep)
        info['text_removed'] = [c[4] for c, kp in zip(cands, keep) if not kp]
    if solid3d:
        body = _header(scale, LAYERS_SPLIT if split3d else LAYERS, _COMMENT_3D) + W.out + _FOOTER
    else:
        body = _header(scale) + W.out + _FOOTER
    info.update({'nentities': sum(W.n.values()), 'counts': dict(W.n), 'scale': scale})
    return body, info


def write_dxf(src, path: str, region=None, scale: float = 1.0, text_height: float = 0.15, names: bool = True,
              mark_size: float = 0.1, mark_text: bool = False, arc_step: float = 0.05, text_side: float = -1,
              spec: Optional[Dict] = None, overlap: str = 'keep', magnet_spec=None, magnet_type=None,
              mode: str = 'plan', parts: str = 'body') -> Dict:
    """Write a plan view of the lattice to a DXF (ASCII, R12 equivalent) file and return a dict of entity counts.
    ラティスの平面図を DXF (ASCII、R12 相当) に書く。戻り値は {'nmagnets','nmarks','ntraj','nentities',...}

    src        USE 済みの `FFS` オブジェクト (`load_ffs(...)` または `run_sad_string(...)[0]`)
    region     (始点, 終点) の要素 (位置の番号または `NAME.n` の名前)。終点の要素も含む。None なら全周 (1 〜 LINE["LENGTH"]-1)
    scale      m からの倍率 (1: m、INSUNITS=6。1000: mm、4。100: cm、5)
    text_height  要素名の文字の高さ (m 単位。scale は掛かる)。names=False で要素名を出さない
    mark_size  MARK の線の長さ (m)。mark_text=True で MARK の名前も出す
    arc_step   BEND の軌道の折れ線の最大の角度 (rad)
    text_side  要素名を置く側 (-1: 進行方向の右 = 周回する CCW のリングの外側、+1: 左、0: 磁石の上)
    spec       磁石の寸法 {名前 または 型名 ("BEND","QUAD","SEXT","OCT","SOL","CAVI","MULT"): {"Width":..,
               "PoleWidth":.., "Type": "Htype"|"Ctype1"|"Ctype2", "PhysicalLength":..}} (m)。要素名の指定が型名より優先。
               無いときの既定は BEND の Width 0.63・PoleWidth 0.13・Htype、QUAD 0.45、SEXT・OCT・MULT 0.30、SOL・CAVI 0.40
    overlap    'keep' (既定: 従来どおり全部出す) か 'remove' (要素名の外接矩形が重なるものを、偏向磁石 > 四極 > 6 極 > 8 極 > その他 > MARK の
               優先で間引く。矩形は文字高と位置から見積もる: 幅 = 文字数 × 文字高 × 0.6、高さ = 文字高、向きは TEXT の回転。
               戻り値に ntext_candidates・ntext_removed・text_removed)
    magnet_spec  磁石の寸法の定義 `MagnetSpec[MagnetType[要素名]]` (drawmagvw3d.n・draw_cad.sad と同じ。規則のリスト
               {Width->, PoleWidth->, Type->, PhysicalLength->, Height->(平面図では使わない)})。None (既定): src の台本が定義していれば自動で読む
               (pysad のインタプリタで評価。定義がなければ何も変わらない)。False: 読まない。dict: {キー: {規則}} を直接与える
               (キーは magnet_type(要素名) の戻り、なければ要素名)。magnet_type は要素名 → キーの関数 (dict のときだけ使う)。
               優先順は量ごとに spec の要素名 > spec の型名 > MagnetSpec > 既定値。使った定義は戻り値の magnet_spec_used
    mode       'plan' ('2d'。既定: 平面図。出力は mode を足す前と 1 行も変わらない) か '3d' (磁石の立体: 2D の外形を z = -Height/2 … +Height/2 に
               押し出した 6 面体を 3DFACE 6 枚で書く。偏向磁石は入口・出口の面が E1・E2 に傾いた楔形。2D の閉じた POLYLINE の代わり。
               基準軌道・MARK・要素名はそのまま (z = 0)。戻り値に nsolids・solids)。Height は spec・MagnetSpec の Height (なければ BEND 0.406 m、
               他は Width)。コイル・端板・ベローズは描かない
    parts      mode='3d' のときだけ。'body' (既定: 磁石 1 台 = 6 面体 1 個) か 'split' (偏向磁石と四極をヨークとポールに分ける。ポールはレイヤー POLE
               (色 9)。偏向磁石: 上下のヨークの板・側脚・上下のポール。四極: 枠 4 本・ポール 4 個 (模式図)。spec・MagnetSpec の Gap・PoleHeight・LegWidth
               (偏向磁石)、Aperture・LegWidth・PoleWidth (四極) で寸法を与える。無ければ既定の割合 (docs/22 §5)。入らない寸法は ValueError)
    """
    body, info = build_dxf(src, region, scale, text_height, names, mark_size, mark_text, arc_step, text_side, spec,
                           overlap, magnet_spec, magnet_type, mode, parts)
    with open(path, 'w', encoding='ascii', newline='\n') as fh:
        fh.write('\n'.join(body) + '\n')
    info['path'] = path
    return info


# ---------------------------------------------------------------- 最小の読み戻し (検証用。ezdxf は使わない)
def read_dxf(path: str) -> Dict:
    """Minimal parser that reads back a DXF written by `write_dxf`, returning insunits, layers and entities.
    write_dxf の DXF を読み戻す最小のパーサ。戻り {'insunits','layers': {名前: 色}, 'entities': [dict ...]}
    entity: {'type','layer','pts' (POLYLINE・LINE の頂点 (x, y)。3DFACE は 4 頂点の (x, y, z)), 'closed', 'text', 'height', 'rot', 'at'}"""
    with open(path, encoding='ascii') as fh:
        raw = [ln.rstrip('\r\n') for ln in fh]
    if len(raw) % 2:
        raise ValueError('DXF の行数が奇数')
    pairs = [(int(raw[k].strip()), raw[k + 1].strip()) for k in range(0, len(raw), 2)]
    res = {'insunits': None, 'layers': {}, 'entities': [], 'acadver': None, 'comment': None}
    sec = None
    cur = None
    table = None
    k = 0
    n = len(pairs)
    var = None
    layer_cur = None
    while k < n:
        c, v = pairs[k]
        if c == 999:
            res['comment'] = v
        elif c == 0 and v == 'SECTION':
            sec = pairs[k + 1][1]
            k += 1
        elif c == 0 and v == 'ENDSEC':
            sec = None
        elif c == 0 and v == 'EOF':
            res['eof'] = True
            break
        elif sec == 'HEADER':
            if c == 9:
                var = v
            elif var == '$INSUNITS' and c == 70:
                res['insunits'] = int(v)
            elif var == '$ACADVER' and c == 1:
                res['acadver'] = v
        elif sec == 'TABLES':
            if c == 0 and v == 'TABLE':
                table = pairs[k + 1][1]
                k += 1
            elif c == 0 and v == 'LAYER':
                layer_cur = {}
            elif table == 'LAYER' and layer_cur is not None:
                if c == 2:
                    layer_cur['name'] = v
                elif c == 62:
                    layer_cur['color'] = int(v)
                    res['layers'][layer_cur['name']] = int(v)
        elif sec == 'ENTITIES':
            if c == 0:
                if v in ('POLYLINE', 'LINE', 'TEXT', '3DFACE'):
                    cur = {'type': v, 'layer': None, 'pts': [], 'closed': False}
                    res['entities'].append(cur)
                elif v == 'VERTEX':
                    cur = res['entities'][-1]
                    cur['_vx'] = True
                elif v == 'SEQEND':
                    cur = None
            elif cur is not None:
                e = cur
                if c == 8:
                    e['layer'] = v
                elif e['type'] == 'POLYLINE':
                    if c == 70 and not e.get('_vx'):
                        e['closed'] = bool(int(v) & 1)
                    elif e.get('_vx'):
                        if c == 10:
                            e['pts'].append([float(v), 0.0])
                        elif c == 20:
                            e['pts'][-1][1] = float(v)
                elif e['type'] == '3DFACE':
                    if 10 <= c <= 13 or 20 <= c <= 23 or 30 <= c <= 33:
                        i, ax = c % 10, c // 10 - 1
                        while len(e['pts']) <= i:
                            e['pts'].append([0.0, 0.0, 0.0])
                        e['pts'][i][ax] = float(v)
                elif e['type'] == 'LINE':
                    if c == 10:
                        e['pts'] = [[float(v), 0.0], [0.0, 0.0]]
                    elif c == 20:
                        e['pts'][0][1] = float(v)
                    elif c == 11:
                        e['pts'][1][0] = float(v)
                    elif c == 21:
                        e['pts'][1][1] = float(v)
                elif e['type'] == 'TEXT':
                    if c == 10:
                        e['at'] = [float(v), 0.0]
                    elif c == 20:
                        e['at'][1] = float(v)
                    elif c == 40:
                        e['height'] = float(v)
                    elif c == 1:
                        e['text'] = v
                    elif c == 50:
                        e['rot'] = float(v)
        k += 1
    for e in res['entities']:
        e.pop('_vx', None)
        e['pts'] = [tuple(p) for p in e['pts']]
    return res


def check_solids(entities, faces_per_solid: int = 6) -> List[Dict]:
    """Check that groups of faces_per_solid 3DFACE entities form closed solids (for verification).
    3DFACE を、ファイルの順に `faces_per_solid` 枚 (6 面体 = 6) ずつ 1 個の立体とみなして、閉じた立体か確かめる (検証用)。
    `entities` は read_dxf の戻りの 'entities' (3DFACE 以外は無視)。立体ごとに dict を返す:
      'layer'・'nfaces'・'nverts' (相異なる頂点の数)・'closed' (すべての有向辺が 1 回ずつ現れ、逆向きの辺も 1 回ずつ現れる = 各辺が 2 面で共有され、
      向きが外向きにそろう)・'volume' (符号つきの体積。外向きの法線なら正)・'bbox'。頂点は小数 6 桁に丸めた値で同一とする"""
    faces = [e for e in entities if e['type'] == '3DFACE']
    res = []
    for k0 in range(0, len(faces), faces_per_solid):
        chunk = faces[k0:k0 + faces_per_solid]
        key = lambda p: (round(p[0], 6) + 0.0, round(p[1], 6) + 0.0, round(p[2], 6) + 0.0)
        edges: Dict = {}
        verts = set()
        vol = 0.0
        for f in chunk:
            pts = [key(p) for p in f['pts']]
            verts.update(pts)
            for a, b in zip(pts, pts[1:] + pts[:1]):
                edges[(a, b)] = edges.get((a, b), 0) + 1
            for t in ((pts[0], pts[1], pts[2]), (pts[0], pts[2], pts[3])):      # 4 角形 = 三角形 2 枚。原点からの四面体の符号つき体積の和
                (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = t
                vol += (x1 * (y2 * z3 - y3 * z2) - y1 * (x2 * z3 - x3 * z2) + z1 * (x2 * y3 - x3 * y2)) / 6.0
        closed = all(n == 1 and edges.get((b, a), 0) == 1 for (a, b), n in edges.items())
        xs = [v[0] for v in verts]
        ys = [v[1] for v in verts]
        zs = [v[2] for v in verts]
        res.append({'layer': chunk[0]['layer'], 'nfaces': len(chunk), 'nverts': len(verts), 'closed': closed and len(chunk) == faces_per_solid,
                    'volume': vol, 'bbox': (min(xs), max(xs), min(ys), max(ys), min(zs), max(zs))})
    return res
