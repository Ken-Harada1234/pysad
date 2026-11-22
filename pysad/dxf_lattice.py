"""
dxf_lattice.py — ラティスの平面図 (上から見た図) をテキストの DXF (ASCII、AutoCAD R12 相当) で書き出す (v2.4.1)

K. Oide の `drawmagvw3d.n` (SAD の 3GSR/esrf。CAD ソフト Vectorworks の VectorScript を書き出す `DrawMagnets3D`) の
幾何の考え方 (LINE["GEO"]・LINE["TYPE"]・LINE["ANGLE"]・LINE["E1"/"E2"]・LINE["DIR"] から偏向磁石・四極などの外形を作る)
を参考にした自前版。出力を DXF にし、立体 (3D) は作らない (平面図だけ)。
同じ計算を SADScript で書いた `example/DrawMagnetsDXF.sad` (gs でも pysad でも走る) の出力と一致する (試験 `test_dxf_lattice.py`)。

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
限界は docs/51_lattice_dxf.md。
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
LAYERS = [('0', 7), ('BEND', 1), ('QUAD', 5), ('SEXT', 3), ('OCT', 6), ('MULT', 30), ('SOL', 4), ('CAVI', 2),
          ('MARK', 8), ('TRAJECTORY', 7), ('TEXT', 7)]
_COMMENT = 'plan view of the lattice (DrawMagnetsDXF); units: see $INSUNITS'


def _fmt(x: float, scale: float) -> str:
    v = x * scale
    if abs(v) < 5e-7:
        v = 0.0
    return '%.6f' % v


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

    def text(self, layer: str, p, height: float, rot_deg: float, s: str):
        o = self.out
        x, y = self.num(p[0]), self.num(p[1])
        o += ['  0', 'TEXT', '  8', layer, ' 10', x, ' 20', y, ' 30', '0.000000', ' 40', self.num(height), '  1', s,
              ' 50', '%.4f' % (0.0 if abs(rot_deg) < 5e-5 else rot_deg), ' 72', '1', ' 11', x, ' 21', y, ' 31', '0.000000', ' 73', '2']
        self.n['TEXT'] += 1


def _insunits(scale: float) -> str:
    return {1.0: '6', 1000.0: '4', 100.0: '5'}.get(float(scale), '0')


def _header(scale: float) -> List[str]:
    h = ['999', _COMMENT, '  0', 'SECTION', '  2', 'HEADER', '  9', '$ACADVER', '  1', 'AC1009',
         '  9', '$INSUNITS', ' 70', _insunits(scale), '  0', 'ENDSEC',
         '  0', 'SECTION', '  2', 'TABLES',
         '  0', 'TABLE', '  2', 'LTYPE', ' 70', '1',
         '  0', 'LTYPE', '  2', 'CONTINUOUS', ' 70', '0', '  3', 'Solid line', ' 72', '65', ' 73', '0', ' 40', '0.0',
         '  0', 'ENDTAB',
         '  0', 'TABLE', '  2', 'LAYER', ' 70', str(len(LAYERS))]
    for name, col in LAYERS:
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
    """SAD の台本 (ファイル名または本文) を読み、USE して FFS オブジェクトを返す (CALC はしない。幾何だけなら不要)。
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


def build_dxf(ffs, region=None, scale: float = 1.0, text_height: float = 0.15, names: bool = True,
              mark_size: float = 0.1, mark_text: bool = False, arc_step: float = 0.05, text_side: float = -1,
              spec: Optional[Dict] = None) -> Tuple[List[str], Dict]:
    """DXF の本文 (行のリスト) と、要素数などの情報の dict を返す。引数は write_dxf を見よ"""
    spec = spec or {}
    P = _Plan(ffs)
    nlat = P.nlat
    if region is None:
        region = (None, None)
    li = _region_index(ffs, region[0], nlat, True)
    lf = min(_region_index(ffs, region[1], nlat, False), nlat - 1)
    if li < 1 or li > lf:
        raise ValueError('region が空: %r' % (region,))
    W = _Writer(scale)
    info = {'li': li, 'lf': lf, 'nmagnets': 0, 'nmarks': 0, 'ntraj': 0}

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
        for k in (name_el, _TYPE_NAME.get(typ_code, '?')):
            d = spec.get(k)
            if isinstance(d, dict) and key in d:
                return d[key]
        return default

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
        else:
            phys = lookup(ename, it, 'PhysicalLength', None)
            lc = lc0 if phys is None else float(phys)
            ds = (lc0 - lc) / 2.0
            poly = [(ds, width / 2.0), (ds, -width / 2.0), (ds + lc, -width / 2.0), (ds + lc, width / 2.0)]
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
            W.text('TEXT', fr(cu, vt), text_height, rot, ename)

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
                W.text('TEXT', fr(0.0, h + text_height), text_height, 0.0, nm)

    body = _header(scale) + W.out + _FOOTER
    info.update({'nentities': sum(W.n.values()), 'counts': dict(W.n), 'scale': scale})
    return body, info


def write_dxf(src, path: str, region=None, scale: float = 1.0, text_height: float = 0.15, names: bool = True,
              mark_size: float = 0.1, mark_text: bool = False, arc_step: float = 0.05, text_side: float = -1,
              spec: Optional[Dict] = None) -> Dict:
    """ラティスの平面図を DXF (ASCII、R12 相当) に書く。戻り値は {'nmagnets','nmarks','ntraj','nentities',...}

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
    """
    body, info = build_dxf(src, region, scale, text_height, names, mark_size, mark_text, arc_step, text_side, spec)
    with open(path, 'w', encoding='ascii', newline='\n') as fh:
        fh.write('\n'.join(body) + '\n')
    info['path'] = path
    return info


# ---------------------------------------------------------------- 最小の読み戻し (検証用。ezdxf は使わない)
def read_dxf(path: str) -> Dict:
    """write_dxf の DXF を読み戻す最小のパーサ。戻り {'insunits','layers': {名前: 色}, 'entities': [dict ...]}
    entity: {'type','layer','pts' (POLYLINE・LINE の頂点 (x, y)), 'closed', 'text', 'height', 'rot', 'at'}"""
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
                if v in ('POLYLINE', 'LINE', 'TEXT'):
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
