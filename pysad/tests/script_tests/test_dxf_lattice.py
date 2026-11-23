"""ラティスの平面図の DXF (pysad/dxf_lattice.py、example/dxf_lattice/DrawMagnetsDXF.n)。
(1) 読み戻し (自前の最小パーサ) のエンティティ数・レイヤー・色・単位、(2) 偏向磁石の外形の幾何 (入口・出口の GEO が端面の上にある、
辺が弦に平行、幅)、(3) 閉じたリングで始点と終点が合う、(4) Region・scale・Spec・PhysicalLength、
(5) SADScript の `DrawMagnetsDXF[]` (pysad で走らせる) と Python 関数の出力が一致。gs との一致は
example/dxf_lattice/run_dxf_examples.py (gs を WSL で走らせる) で確かめる。"""
import io
import math
import os

import pytest

from pysad.dxf_lattice import build_dxf, load_ffs, read_dxf, write_dxf

HERE = os.path.dirname(os.path.abspath(__file__))
EXDIR = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'example', 'dxf_lattice'))

LAT = r"""MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
OCT   OC = (L=0.1 K3=0);
SOL   SA = (L=0 BZ=0.1 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0.4 VOLT=1.5E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1 OC DS SA D2 SB D2 RF D2);
LINE RING = (ORIGIN 24*CELL);
"""

# 矩形の偏向磁石 (E1 = E2 = 0.5: 端面が弦に直角)・逆向きの偏向・直線だけ
LAT_RECT = r"""MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0);
BEND BR = (L=1.0 ANGLE=0.2 E1=0.5 E2=0.5) BS = (L=1.0 ANGLE=-0.1 E1=0 E2=0) B0 = (L=0.5 ANGLE=0);
MARK M1 = ();
LINE T = (M1 D1 BR D1 BS D1 B0 D1 M1);
"""


def _ffs(text=LAT, use='RING'):
    return load_ffs(text, use)


def _lab(poly_pts):
    return poly_pts


def test_roundtrip_counts_layers_units(tmp_path):
    f = str(tmp_path / 'ring.dxf')
    info = write_dxf(_ffs(), f)
    d = read_dxf(f)
    assert d['acadver'] == 'AC1009' and d['insunits'] == 6 and d.get('eof')
    assert d['layers'] == {'0': 7, 'BEND': 1, 'QUAD': 5, 'SEXT': 3, 'OCT': 6, 'MULT': 30, 'SOL': 4, 'CAVI': 2,
                           'MARK': 8, 'TRAJECTORY': 7, 'TEXT': 7}
    cnt = {}
    for e in d['entities']:
        key = (e['type'], e['layer'])
        cnt[key] = cnt.get(key, 0) + 1
    assert cnt[('POLYLINE', 'TRAJECTORY')] == 1
    assert cnt[('POLYLINE', 'BEND')] == 24
    assert cnt[('POLYLINE', 'QUAD')] == 48
    assert cnt[('POLYLINE', 'SEXT')] == 48
    assert cnt[('POLYLINE', 'OCT')] == 24
    assert cnt[('POLYLINE', 'SOL')] == 24
    assert cnt[('POLYLINE', 'CAVI')] == 24
    assert cnt[('LINE', 'MARK')] == 1
    assert cnt[('TEXT', 'TEXT')] == 24 * 8
    assert info['nentities'] == len(d['entities'])
    for e in d['entities']:
        if e['type'] == 'POLYLINE' and e['layer'] != 'TRAJECTORY':
            assert e['closed'] and len(e['pts']) == 4
    names = sorted({e['text'] for e in d['entities'] if e['type'] == 'TEXT'})
    assert names == ['B', 'OC', 'QD', 'QF', 'RF', 'SA', 'SD', 'SF']


def test_ring_closes(tmp_path):
    """(d) 閉じたリングで基準軌道の始点と終点が合う (BEND の円弧と GEO の整合)。軌道の折れ線の長さは周長に近い"""
    f = str(tmp_path / 'ring.dxf')
    write_dxf(_ffs(), f)
    tr = [e for e in read_dxf(f)['entities'] if e['layer'] == 'TRAJECTORY'][0]
    (x0, y0), (x1, y1) = tr['pts'][0], tr['pts'][-1]
    assert math.hypot(x1 - x0, y1 - y0) < 1e-5
    ffs = _ffs()
    g = ffs.line_value('GEO', ffs.nlat)
    assert abs(g[1][0]) < 1e-9 and math.hypot(g[0][0], g[0][1]) < 1e-9       # SAD の GEO 自体が閉じる (ANGLE の和が 2π)
    # 折れ線の長さは弧の弦の和なので周長より僅かに短い (ArcStep 0.05 rad: 弦/弧 = 1 - a²/24 ≈ 1e-4)
    L = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(tr['pts'][:-1], tr['pts'][1:]))
    Ls = ffs.line_value('S', ffs.nlat)
    assert 0.9995 * Ls < L <= Ls + 1e-6


def _poly_of(d, layer, k=0):
    return [e for e in d['entities'] if e['type'] == 'POLYLINE' and e['layer'] == layer][k]['pts']


def _dist_to_seg(p, a, b):
    ax, ay = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * ax + (p[1] - a[1]) * ay) / (ax * ax + ay * ay)
    t = max(0.0, min(1.0, t))
    return math.hypot(p[0] - (a[0] + t * ax), p[1] - (a[1] + t * ay))


def test_bend_outline_geometry(tmp_path):
    """偏向磁石の外形 (SAD の幾何から独立に確かめる): 入口・出口の GEO の点が端面の上、端面の向きが E1・E2、
    左右の辺が弦に平行で間隔 = Width、軌道は左右の辺の中央から ±hs"""
    ffs = _ffs(LAT_RECT, 'T')
    f = str(tmp_path / 't.dxf')
    write_dxf(ffs, f, spec={'BEND': {'Width': 0.6}}, names=False)
    d = read_dxf(f)
    # BR: 矩形 (ANGLE 0.2、E=0.5)
    i0 = [i for i in range(1, ffs.nlat) if ffs.line_value('ELEMENT', i) == 'BR'][0]
    gin, gout = ffs.line_value('GEO', i0), ffs.line_value('GEO', i0 + 1)
    pin, pout = tuple(gin[0][:2]), tuple(gout[0][:2])
    P = _poly_of(d, 'BEND', 0)
    assert _dist_to_seg(pin, P[0], P[1]) < 1e-6
    assert _dist_to_seg(pout, P[2], P[3]) < 1e-6
    ang = 0.2
    chord = (pout[0] - pin[0], pout[1] - pin[1])
    cl = math.hypot(*chord)
    ux, uy = chord[0] / cl, chord[1] / cl
    # 端面は弦に直角 (矩形): P0-P1 と弦の内積 = 0
    for a, b in ((P[0], P[1]), (P[2], P[3])):
        assert abs((b[0] - a[0]) * ux + (b[1] - a[1]) * uy) < 2e-6
    # 左右の辺は弦に平行: P1-P2 (右の辺)・P3-P0 (左の辺)
    for a, b in ((P[1], P[2]), (P[3], P[0])):
        cr = (b[0] - a[0]) * uy - (b[1] - a[1]) * ux
        assert abs(cr) < 2e-6                                   # 出力は 1e-6 に丸めてある
    # 幅 = 0.6、長さ = 弦
    assert abs(math.hypot(P[1][0] - P[0][0], P[1][1] - P[0][1]) - 0.6) < 2e-6
    assert abs(math.hypot(P[2][0] - P[1][0], P[2][1] - P[1][1]) - cl) < 2e-6
    # 軌道は箱の中心線から、入口・出口で左に hs、中央で右に hs (hs = la/angle · sin²(angle/4)、la = 弧の長さ = L = 1.0)
    hs = 1.0 / ang * math.sin(ang / 4) ** 2
    dl = math.hypot(P[0][0] - pin[0], P[0][1] - pin[1])          # 入口の左の辺までの距離 = 幅/2 - hs
    dr = math.hypot(P[1][0] - pin[0], P[1][1] - pin[1])          # 右の辺まで = 幅/2 + hs
    assert abs(dl - (0.3 - hs)) < 2e-5 and abs(dr - (0.3 + hs)) < 2e-5
    # 逆向きの偏向 BS (ANGLE<0): 弧は弦の左、箱は左にずれる
    P2 = _poly_of(d, 'BEND', 1)
    i1 = [i for i in range(1, ffs.nlat) if ffs.line_value('ELEMENT', i) == 'BS'][0]
    a, b = ffs.line_value('GEO', i1), ffs.line_value('GEO', i1 + 1)
    pin2, pout2 = a[0][:2], b[0][:2]
    ch = (pout2[0] - pin2[0], pout2[1] - pin2[1])
    c2 = math.hypot(*ch)
    ux, uy = ch[0] / c2, ch[1] / c2
    cen2 = [sum(p[0] for p in P2) / 4, sum(p[1] for p in P2) / 4]
    mid2 = [(pin2[0] + pout2[0]) / 2, (pin2[1] + pout2[1]) / 2]
    off2 = (cen2[0] - mid2[0]) * (-uy) + (cen2[1] - mid2[1]) * ux
    assert off2 > 0
    # B0 (ANGLE=0): 矩形 (長さ 0.5、幅 0.6) で軌道が中央
    P3 = _poly_of(d, 'BEND', 2)
    assert abs(math.hypot(P3[2][0] - P3[1][0], P3[2][1] - P3[1][1]) - 0.5) < 2e-6
    assert abs(math.hypot(P3[1][0] - P3[0][0], P3[1][1] - P3[0][1]) - 0.6) < 2e-6


def test_sector_bend_edges_follow_E1_E2(tmp_path):
    """E1・E2 (ANGLE に対する比) を変えると端面が回る。E1=E2=0 (セクター) の端面は入口・出口の向きに直角"""
    t = r"""MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0);
BEND BA = (L=1.0 ANGLE=0.3 E1=0 E2=0) BB = (L=1.0 ANGLE=0.3 E1=0.8 E2=0.2);
LINE T = (D1 BA D1 BB D1);
"""
    ffs = _ffs(t, 'T')
    f = str(tmp_path / 't.dxf')
    write_dxf(ffs, f, names=False)
    d = read_dxf(f)
    for k, name, e1r, e2r in ((0, 'BA', 0.0, 0.0), (1, 'BB', 0.8, 0.2)):
        i = [j for j in range(1, ffs.nlat) if ffs.line_value('ELEMENT', j) == name][0]
        gi, go = ffs.line_value('GEO', i), ffs.line_value('GEO', i + 1)
        P = _poly_of(d, 'BEND', k)
        # 端面の向き (左 - 右) と入口の進行方向 chi1 の角度: 直角から e1 = E1·ANGLE だけ傾く (SAD の定義)
        # 入口の面は P[0] (左) と P[1] (右)
        v = (P[0][0] - P[1][0], P[0][1] - P[1][1])
        ang_in = (math.atan2(v[1], v[0]) - (gi[1][0] + math.pi / 2) + math.pi) % (2 * math.pi) - math.pi
        assert abs(ang_in - e1r * 0.3) < 1e-6, (name, ang_in)
        v = (P[3][0] - P[2][0], P[3][1] - P[2][1])         # 出口の面: 右 P[2] → 左 P[3]
        ang_out = (math.atan2(v[1], v[0]) - (go[1][0] + math.pi / 2) + math.pi) % (2 * math.pi) - math.pi
        assert abs(ang_out + e2r * 0.3) < 1e-6, (name, ang_out)


def test_scale_region_and_spec(tmp_path):
    ffs = _ffs()
    f1, f2 = str(tmp_path / 'm.dxf'), str(tmp_path / 'mm.dxf')
    write_dxf(ffs, f1, region=('QF.1', 'QF.2'))
    write_dxf(ffs, f2, region=('QF.1', 'QF.2'), scale=1000)
    d1, d2 = read_dxf(f1), read_dxf(f2)
    assert d1['insunits'] == 6 and d2['insunits'] == 4
    # Region: QF.1 から QF.2 まで (終点を含む) = 1 セル分の磁石 + 次のセルの SF・(QF) だけ
    nm = [e['text'] for e in d1['entities'] if e['type'] == 'TEXT']
    assert nm.count('QF') == 2 and nm.count('B') == 1 and nm.count('SF') == 1 and nm.count('SD') == 1
    for e1, e2 in zip(d1['entities'], d2['entities']):
        assert e1['type'] == e2['type'] and len(e1['pts']) == len(e2['pts'])
        for p, q in zip(e1['pts'], e2['pts']):
            assert abs(q[0] - 1000 * p[0]) < 2e-3 and abs(q[1] - 1000 * p[1]) < 2e-3
    # 数字で指定しても同じ
    f3 = str(tmp_path / 'n.dxf')
    i1 = int(ffs.fn_line(['POSITION', 'QF.1']))
    i2 = int(ffs.fn_line(['POSITION', 'QF.2']))
    write_dxf(ffs, f3, region=(i1, i2))
    assert open(f3).read() == open(f1).read()
    # Spec: 名前 > 型名 > 既定
    f4 = str(tmp_path / 's.dxf')
    write_dxf(ffs, f4, spec={'QUAD': {'Width': 0.2}, 'QF': {'Width': 0.5, 'PhysicalLength': 0.2}}, names=False)
    d4 = read_dxf(f4)
    qs = [e for e in d4['entities'] if e['layer'] == 'QUAD']
    qf, qd = qs[0], qs[1]
    wid = lambda e: math.hypot(e['pts'][1][0] - e['pts'][0][0], e['pts'][1][1] - e['pts'][0][1])
    ln = lambda e: math.hypot(e['pts'][2][0] - e['pts'][1][0], e['pts'][2][1] - e['pts'][1][1])
    assert abs(wid(qf) - 0.5) < 1e-5 and abs(ln(qf) - 0.2) < 1e-5     # 物理長 0.2 (L=0.3 の中央)
    assert abs(wid(qd) - 0.2) < 1e-5 and abs(ln(qd) - 0.3) < 1e-5
    # 物理長が中央: 矩形の中心は L=0.3 の入口から 0.15 のところ
    gq = ffs.line_value('GEO', int(ffs.fn_line(['POSITION', 'QF.1'])))[0]
    cx = sum(p[0] for p in qf['pts']) / 4
    cy = sum(p[1] for p in qf['pts']) / 4
    assert abs(math.hypot(cx - gq[0], cy - gq[1]) - 0.15) < 1e-5


def test_ctype_asymmetry(tmp_path):
    """Ctype1: ANGLE>0 は PoleWidth/2 が右 (外側)、残りが左 (内側)。Ctype2 は逆。Htype は左右対称"""
    ffs = _ffs(LAT_RECT, 'T')
    res = {}
    for typ in ('Htype', 'Ctype1', 'Ctype2'):
        f = str(tmp_path / (typ + '.dxf'))
        write_dxf(ffs, f, spec={'BR': {'Width': 0.5, 'PoleWidth': 0.1, 'Type': typ}}, names=False)
        P = _poly_of(read_dxf(f), 'BEND', 0)
        # 局所の左右: 入口の面 P[0] (左)・P[1] (右) の軌道の点からの距離
        i = [j for j in range(1, ffs.nlat) if ffs.line_value('ELEMENT', j) == 'BR'][0]
        g = ffs.line_value('GEO', i)
        c = (g[0][0], g[0][1])
        res[typ] = (math.hypot(P[0][0] - c[0], P[0][1] - c[1]), math.hypot(P[1][0] - c[0], P[1][1] - c[1]))
    hs = 1.0 / 0.2 * math.sin(0.2 / 4) ** 2
    exp = {'Htype': (0.25 - hs, 0.25 + hs), 'Ctype1': (0.45 - hs, 0.05 + hs), 'Ctype2': (0.05 - hs, 0.45 + hs)}
    for typ in exp:
        for a_, b_ in zip(res[typ], exp[typ]):
            assert abs(a_ - b_) < 2e-5, (typ, res[typ], exp[typ])


def test_script_matches_python(tmp_path):
    """SADScript の DrawMagnetsDXF[] (pysad で走らせる) と Python 関数の出力が一致する (3 通りのオプション)"""
    from pysad.ffscmd import run_sad_string
    spec_s = 'Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1}, "QF"->{Width->.3, PhysicalLength->.26}, "QUAD"->{Width->.36}}'
    spec_p = {'B': {'Width': .5, 'PoleWidth': .15, 'Type': 'Ctype1'}, 'QF': {'Width': .3, 'PhysicalLength': .26}, 'QUAD': {'Width': .36}}
    cases = [('a.dxf', '', {}),
             ('b.dxf', ', ' + spec_s + ', MarkText->True, TextSide->1', {'spec': spec_p, 'mark_text': True, 'text_side': 1}),
             ('c.dxf', ', Region->{"QF.2","B.5"}, Scale->100, Names->False, ArcStep->.1',
              {'region': ('QF.2', 'B.5'), 'scale': 100, 'names': False, 'arc_step': .1})]
    prog = LAT + 'FFS;\nUSE RING;\nGet["%s"];\n' % os.path.join(EXDIR, 'DrawMagnetsDXF.n').replace('\\', '/')
    for fn, opts, kw in cases:
        out = str(tmp_path / fn).replace('\\', '/')
        buf = io.StringIO()
        ffs, _ = run_sad_string(prog + 'DrawMagnetsDXF["%s"%s];\nExit;\n' % (out, opts), out=buf)
        assert '???' not in buf.getvalue(), buf.getvalue()[-500:]
        py = str(tmp_path / ('py_' + fn))
        write_dxf(ffs, py, **kw)
        a = open(out, newline='').read().replace('\r\n', '\n')
        b = open(py, newline='').read().replace('\r\n', '\n')
        assert a == b, fn


def test_example_files_and_closure():
    """example/dxf_lattice の DXF (追跡対象のもの) を読み戻し、レイヤー・エンティティ数・基準軌道が整合している"""
    p = os.path.join(EXDIR, 'DXF_FODO_ring.dxf')
    if not os.path.exists(p):
        pytest.skip('example/dxf_lattice/DXF_FODO_ring.dxf なし')
    d = read_dxf(p)
    assert d['insunits'] == 6
    ncl = sum(1 for e in d['entities'] if e['type'] == 'POLYLINE' and e['closed'])
    assert ncl == 120
    tr = [e for e in d['entities'] if e['layer'] == 'TRAJECTORY'][0]
    assert math.hypot(tr['pts'][0][0] - tr['pts'][-1][0], tr['pts'][0][1] - tr['pts'][-1][1]) < 1e-5


def test_empty_region_error():
    with pytest.raises(ValueError):
        build_dxf(_ffs(), region=(10, 5))


def test_ezdxf_readback(tmp_path):
    """ezdxf が環境にあれば、ezdxf でも読めて、エンティティ数とレイヤーが自前のパーサと同じ (無ければ skip)"""
    ezdxf = pytest.importorskip('ezdxf')
    f = str(tmp_path / 'ring.dxf')
    write_dxf(_ffs(), f)
    doc = ezdxf.readfile(f)
    msp = doc.modelspace()
    mine = read_dxf(f)
    assert len(list(msp)) == len(mine['entities'])
    assert {l.dxf.name for l in doc.layers} >= set(mine['layers'])
    for e, m in zip(msp, mine['entities']):
        assert e.dxftype() == m['type'] and e.dxf.layer == m['layer']


# ---------------------------------------------------------------- v2.4.1: MagnetSpec[] の読み込みと要素名の重なり消し
MS_DEFS = r"""MagnetType[n_] := Switch[n, "QF" | "QD", "QQ", "B", "BX", _, n];
MagnetSpec["QQ"] = {Width -> .31, Height -> .2, PhysicalLength -> .26};
MagnetSpec["BX"] = {Width -> .55, PoleWidth -> .14, Type -> Ctype1, Height -> .4};
MagnetSpec["SF"] = {Width -> .22, PhysicalLength -> None};
"""


def _ms_ffs(defs=MS_DEFS, text=LAT, use='RING'):
    return load_ffs(text + 'FFS;\n' + defs, use)


def _rect_polys(d):
    """TEXT エンティティの外接矩形の頂点 (テスト用の独立な実装: 文字数 × 高さ × 0.6 の矩形を回して頂点にする)"""
    out = []
    for e in d['entities']:
        if e['type'] != 'TEXT':
            continue
        h = e['height']
        hw, hh = h * 0.6 * len(e['text']) / 2.0, h / 2.0
        a = math.radians(e.get('rot', 0.0))
        c, s = math.cos(a), math.sin(a)
        pts = [(e['at'][0] + c * u - s * v, e['at'][1] + s * u + c * v) for u, v in ((-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh))]
        out.append((e['text'], tuple(e['at']), pts))
    return out


def _polys_overlap(p, q, eps=1e-9):
    """2 つの凸多角形が (面積を持って) 重なるか: 辺の法線を分離軸とする"""
    for poly in (p, q):
        for k in range(len(poly)):
            x0, y0 = poly[k]
            x1, y1 = poly[(k + 1) % len(poly)]
            nx, ny = y1 - y0, -(x1 - x0)
            pa = [nx * x + ny * y for x, y in p]
            pb = [nx * x + ny * y for x, y in q]
            if max(pa) <= min(pb) + eps * math.hypot(nx, ny) or max(pb) <= min(pa) + eps * math.hypot(nx, ny):
                return False
    return True


def _wid(e):
    return math.hypot(e['pts'][1][0] - e['pts'][0][0], e['pts'][1][1] - e['pts'][0][1])


def _len(e):
    return math.hypot(e['pts'][2][0] - e['pts'][1][0], e['pts'][2][1] - e['pts'][1][1])


def test_magnet_spec_read_and_priority(tmp_path):
    """MagnetSpec[MagnetType[要素名]] を読む。優先順は量ごとに Spec の要素名 > Spec の型名 > MagnetSpec > 既定値。
    無い要素・量は従来の既定。magnet_spec=False で無視"""
    ffs = _ms_ffs()

    def draw(**kw):
        f = str(tmp_path / 'q.dxf')
        info = write_dxf(ffs, f, names=False, **kw)
        d = read_dxf(f)
        lay = lambda n: [e for e in d['entities'] if e['layer'] == n]
        return lay('QUAD'), lay('BEND'), lay('SEXT'), info
    q, b, sx, info = draw()
    # QF・QD は MagnetType → "QQ": 幅 .31、物理長 .26 (L=.3 の中央)
    assert abs(_wid(q[0]) - .31) < 2e-6 and abs(_len(q[0]) - .26) < 2e-6
    assert abs(_wid(q[1]) - .31) < 2e-6 and abs(_len(q[1]) - .26) < 2e-6
    # B → "BX": 幅 .55、Ctype1 (PoleWidth .14: ANGLE>0 は右の半幅 .07、左 .48)
    P = b[0]['pts']
    i0 = [i for i in range(1, ffs.nlat) if ffs.line_value('ELEMENT', i) == 'B'][0]
    c = ffs.line_value('GEO', i0)[0][:2]
    dl, dr = math.hypot(P[0][0] - c[0], P[0][1] - c[1]), math.hypot(P[1][0] - c[0], P[1][1] - c[1])
    an = 0.2617993877991494
    hs = 1.2 / an * math.sin(an / 4) ** 2
    cs = math.cos(0.1308996938995747 * an - an / 2)          # 端面の傾き (E1 は ANGLE に対する比) による r = (半幅 ∓ hs) / cos(e1 - ANGLE/2)
    assert abs(dl - (.48 - hs) / cs) < 3e-5 and abs(dr - (.07 + hs) / cs) < 3e-5
    # 使った定義の記録 (Height は描かないが読んでいる)
    assert info['magnet_spec_used']['QF']['Width'] == pytest.approx(.31) and info['magnet_spec_used']['B']['Height'] == pytest.approx(.4)
    # SF は PhysicalLength->None (= 指定なし) と Width .22。MagnetSpec に無い SD は既定 (SEXT 0.30)
    assert abs(_wid(sx[0]) - .22) < 2e-6 and abs(_len(sx[0]) - .2) < 2e-6
    assert abs(_wid(sx[1]) - .30) < 2e-6
    # Spec の要素名 > Spec の型名 > MagnetSpec: QF だけ要素名で .5、他の QUAD は型名で .2 (幅だけ。物理長は MagnetSpec のまま)
    q, _b, _s, info = draw(spec={'QUAD': {'Width': .2}, 'QF': {'Width': .5}})
    assert abs(_wid(q[0]) - .5) < 2e-6 and abs(_wid(q[1]) - .2) < 2e-6
    assert abs(_len(q[0]) - .26) < 2e-6 and abs(_len(q[1]) - .26) < 2e-6
    # magnet_spec=False: 従来の既定 (QUAD 0.45、長さ = L)
    q, _b, _s, info = draw(magnet_spec=False)
    assert abs(_wid(q[0]) - .45) < 2e-6 and abs(_len(q[0]) - .3) < 2e-6 and info['magnet_spec_used'] == {}


def test_magnet_spec_dict_and_default_unchanged(tmp_path):
    """dict で与えても (magnet_type つき) 台本から読んだのと同じ DXF。MagnetType なしならキーは要素名 (文字列でも Symbol でも)。
    MagnetSpec を定義しない台本は magnet_spec=False・overlap='keep' と同じ出力 (従来のまま)"""
    ffs = _ms_ffs()
    f1, f2 = str(tmp_path / 'a.dxf'), str(tmp_path / 'b.dxf')
    write_dxf(ffs, f1)
    ms = {'QQ': {'Width': .31, 'Height': .2, 'PhysicalLength': .26}, 'BX': {'Width': .55, 'PoleWidth': .14, 'Type': 'Ctype1', 'Height': .4},
          'SF': {'Width': .22, 'PhysicalLength': None}}
    write_dxf(ffs, f2, magnet_spec=ms, magnet_type=lambda n: 'QQ' if n in ('QF', 'QD') else 'BX' if n == 'B' else n)
    assert open(f1).read() == open(f2).read()
    # MagnetType なし: dict は要素名がキー。台本は MagnetSpec["QF"] (文字列) でも MagnetSpec[QD] (Symbol) でもよい
    f4, f5 = str(tmp_path / 'd.dxf'), str(tmp_path / 'e.dxf')
    write_dxf(ffs, f4, magnet_spec={'QF': {'Width': .31, 'PhysicalLength': .26}, 'QD': {'Width': .31, 'PhysicalLength': .26}})
    ffs_n = load_ffs(LAT + 'FFS;\nMagnetSpec["QF"] = {Width -> .31, PhysicalLength -> .26};\nMagnetSpec[QD] = {Width -> .31, PhysicalLength -> .26};\n', 'RING')
    write_dxf(ffs_n, f5)
    assert open(f4).read() == open(f5).read()
    q = [e for e in read_dxf(f5)['entities'] if e['layer'] == 'QUAD']
    assert abs(_wid(q[0]) - .31) < 2e-6 and abs(_wid(q[1]) - .31) < 2e-6
    # 定義のない台本 = magnet_spec=False = 従来の出力
    plain = _ffs()
    g1, g2 = str(tmp_path / 'p1.dxf'), str(tmp_path / 'p2.dxf')
    write_dxf(plain, g1)
    write_dxf(plain, g2, magnet_spec=False, overlap='keep')
    assert open(g1).read() == open(g2).read()


def test_overlap_remove_valid_and_maximal(tmp_path):
    """overlap='remove': 残った名前どうしは重ならず (独立な多角形の判定)、消した名前はそれぞれ、優先が同じか高い残った名前と重なる。
    残るのは 'keep' (既定) の名前の部分列で、磁石の図形は変わらない"""
    ffs = _ffs()
    fk, fr = str(tmp_path / 'k.dxf'), str(tmp_path / 'r.dxf')
    ik = write_dxf(ffs, fk, text_height=.4)
    ir = write_dxf(ffs, fr, text_height=.4, overlap='remove')
    dk, dr = read_dxf(fk), read_dxf(fr)
    assert ik['counts']['TEXT'] == 192 and 'ntext_removed' not in ik
    assert ir['ntext_candidates'] == 192 and ir['ntext_removed'] == 192 - ir['counts']['TEXT'] > 0
    assert len(ir['text_removed']) == ir['ntext_removed']
    shape = lambda d: [(e['type'], e['layer'], e['pts']) for e in d['entities'] if e['type'] != 'TEXT']
    assert shape(dk) == shape(dr)
    rk, rr = _rect_polys(dk), _rect_polys(dr)
    prio = {'B': 1, 'QF': 2, 'QD': 2, 'SF': 3, 'SD': 3, 'OC': 4, 'SA': 5, 'RF': 5}
    kept_idx, j = [], 0
    for i, (nm, at, _p) in enumerate(rk):                 # 残ったものは 'keep' の部分列 (同じ名前・位置)
        if j < len(rr) and rr[j][0] == nm and rr[j][1] == at:
            kept_idx.append(i)
            j += 1
    assert j == len(rr)
    for a in range(len(rr)):
        for b in range(a + 1, len(rr)):
            assert not _polys_overlap(rr[a][2], rr[b][2]), (rr[a][0], rr[b][0])
    kept = set(kept_idx)
    for i, (nm, at, p) in enumerate(rk):
        if i in kept:
            continue
        assert any(prio[rk[k][0]] <= prio[nm] and _polys_overlap(p, rk[k][2]) for k in kept_idx), (i, nm)


def test_overlap_priority_order_and_ties(tmp_path):
    """優先順は偏向磁石 > 四極 > 6 極 > 8 極 > その他 > MARK。描く順が先の四極より後の偏向磁石の名前を残す。同じ優先なら描く順が先を残す"""
    t = r"""MOMENTUM = 3GEV;
DRIFT DS = (L=0.02) D1 = (L=0.5);
QUAD QA = (L=0.1 K1=0.3) QB = (L=0.1 K1=-0.3);
SEXT SX = (L=0.1 K2=0);
BEND BE = (L=0.1 ANGLE=0.01 E1=0 E2=0);
MARK MM = ();
LINE T = (QA DS BE DS SX DS QB);
LINE T2 = (QA DS QB);
LINE T3 = (QB DS QA);
LINE T4 = (D1 MM QA D1);
"""
    kw = dict(text_height=.3, overlap='remove', text_side=0)       # text_side=0: 名前は磁石の中心 (長さ 0.1 が 0.12 ピッチ、名前の幅 0.36 以上 → 重なる)
    info = write_dxf(_ffs(t, 'T'), str(tmp_path / 'o.dxf'), **kw)
    assert 'BE' not in info['text_removed'] and sorted(info['text_removed']) == ['QA', 'QB', 'SX']
    assert info['counts']['TEXT'] == 1
    assert write_dxf(_ffs(t, 'T2'), str(tmp_path / 'o2.dxf'), **kw)['text_removed'] == ['QB']      # 同じ優先: 先を残す
    assert write_dxf(_ffs(t, 'T3'), str(tmp_path / 'o3.dxf'), **kw)['text_removed'] == ['QA']
    # SEXT より QUAD、QUAD より BEND、BEND を後ろに置いても同じ (QA SX BE)
    t5 = t.replace('LINE T = (QA DS BE DS SX DS QB);', 'LINE T = (SX DS QA DS BE);')
    i5 = write_dxf(_ffs(t5, 'T'), str(tmp_path / 'o5.dxf'), **kw)
    assert i5['text_removed'] == ['SX', 'QA']
    # MARK の名前 (mark_text=True) は最後: 軌道に直角に h + 文字高 だけ離れるので、磁石の名前と離れていれば残り、重なれば消える
    i6 = write_dxf(_ffs(t, 'T4'), str(tmp_path / 'o6.dxf'), mark_text=True, mark_size=0.02, **kw)
    assert i6['text_removed'] == []
    i7 = write_dxf(_ffs(t, 'T4'), str(tmp_path / 'o7.dxf'), mark_text=True, mark_size=0.02, text_height=.3, overlap='remove', text_side=1)
    assert i7['text_removed'] == ['MM']                        # QA の名前は左 (Width/2 + 0.3 = 0.525) で、MM の 0.31 と重なる
    with pytest.raises(ValueError):
        build_dxf(_ffs(t, 'T4'), overlap='thin')


def test_rects_overlap_function():
    """向きのついた矩形の判定は、独立な多角形の判定とランダムな配置で一致する。接するだけは重ならない"""
    import random
    from pysad.dxf_lattice import _rects_overlap, _text_rect
    rng = random.Random(7)
    poly = lambda r: [(r[0] + r[4] * u - r[5] * v, r[1] + r[5] * u + r[4] * v) for u, v in ((-r[2], -r[3]), (r[2], -r[3]), (r[2], r[3]), (-r[2], r[3]))]
    for _ in range(400):
        a = _text_rect(rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(.1, .6), rng.uniform(-90, 90), 'ABCQ'[:rng.randint(1, 4)])
        b = _text_rect(rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(.1, .6), rng.uniform(-90, 90), 'XYZW'[:rng.randint(1, 4)])
        assert _rects_overlap(a, b) == _polys_overlap(poly(a), poly(b))
    assert not _rects_overlap(_text_rect(0, 0, .1, 0, 'AB'), _text_rect(.12, 0, .1, 0, 'AB'))
    assert _rects_overlap(_text_rect(0, 0, .1, 0, 'AB'), _text_rect(.1199, 0, .1, 0, 'AB'))


def test_script_magnet_spec_and_overlap_match_python(tmp_path):
    """SADScript の DrawMagnetsDXF[] (MagnetSpec・MagnetType の定義つき。Overlap->"remove"、UseMagnetSpec->False など) と Python 関数の出力が一致する"""
    from pysad.ffscmd import run_sad_string
    cases = [('a.dxf', ', TextHeight->.4, Overlap->"remove"', {'text_height': .4, 'overlap': 'remove'}),
             ('b.dxf', ', UseMagnetSpec->False, TextHeight->.4, Overlap->"remove", MarkText->True',
              {'magnet_spec': False, 'text_height': .4, 'overlap': 'remove', 'mark_text': True}),
             ('c.dxf', ', Spec->{"QF"->{Width->.5}}, TextSide->0, TextHeight->.25, Overlap->"remove", Region->{"QF.2","B.5"}',
              {'spec': {'QF': {'Width': .5}}, 'text_side': 0, 'text_height': .25, 'overlap': 'remove', 'region': ('QF.2', 'B.5')}),
             ('d.dxf', ', Overlap->"keep"', {'overlap': 'keep'})]
    prog = LAT + 'FFS;\nUSE RING;\n' + MS_DEFS + 'Get["%s"];\n' % os.path.join(EXDIR, 'DrawMagnetsDXF.n').replace('\\', '/')
    n_removed = 0
    for fn, opts, kw in cases:
        out = str(tmp_path / fn).replace('\\', '/')
        buf = io.StringIO()
        ffs, _ = run_sad_string(prog + 'DrawMagnetsDXF["%s"%s];\nExit;\n' % (out, opts), out=buf)
        assert '???' not in buf.getvalue(), buf.getvalue()[-500:]
        py = str(tmp_path / ('py_' + fn))
        info = write_dxf(ffs, py, **kw)
        n_removed += info.get('ntext_removed', 0)
        a = open(out, newline='').read().replace('\r\n', '\n')
        b = open(py, newline='').read().replace('\r\n', '\n')
        assert a == b, fn
    assert n_removed > 0
    # Overlap の値が不正なら画面にメッセージを出して "keep" 扱い
    buf = io.StringIO()
    out = str(tmp_path / 'e.dxf').replace('\\', '/')
    run_sad_string(prog + 'DrawMagnetsDXF["%s", Overlap->"thin"];\nExit;\n' % out, out=buf)
    assert 'Overlap must be' in buf.getvalue()


def test_example_v2_files():
    """example/dxf_lattice の v2 の DXF (KEK-LS の標準セル) を Python 関数で作り直して全行一致する。既存の DXF_KEKLS_cell.dxf も既定のまま一致する"""
    lat = os.path.abspath(os.path.join(EXDIR, '..', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad'))
    need = ['DXF_KEKLS_cell.dxf', 'DXF_KEKLS_cell_v2_remove.dxf', 'DXF_KEKLS_cell_v2_spec.dxf', 'DXF_KEKLS_cell_v2_spec_remove.dxf',
            'DXF_KEKLS_cell_v2_spec_h30.dxf', 'DXF_KEKLS_cell_v2_spec_remove_h30.dxf']
    if not os.path.exists(lat) or not all(os.path.exists(os.path.join(EXDIR, n)) for n in need):
        pytest.skip('KEK-LS のラティスまたは v2 の DXF なし')
    import tempfile
    text = open(lat, encoding='utf-8').read() + '\nFFS;\nGet["MagnetSpec_KEKLS.n"];\n'
    ffs = load_ffs(text, 'ARCA', cwd=EXDIR)
    rd = lambda p: open(p, newline='').read().replace('\r\n', '\n')
    with tempfile.TemporaryDirectory() as td:
        for fn, kw in (('DXF_KEKLS_cell.dxf', {'magnet_spec': False}),                      # 既存 (v2.4.1 の最初の版) と同じ
                       ('DXF_KEKLS_cell_v2_remove.dxf', {'magnet_spec': False, 'overlap': 'remove'}),
                       ('DXF_KEKLS_cell_v2_spec.dxf', {}),
                       ('DXF_KEKLS_cell_v2_spec_remove.dxf', {'overlap': 'remove'}),
                       ('DXF_KEKLS_cell_v2_spec_h30.dxf', {'text_height': .3}),
                       ('DXF_KEKLS_cell_v2_spec_remove_h30.dxf', {'text_height': .3, 'overlap': 'remove'})):
            p = os.path.join(td, fn)
            write_dxf(ffs, p, **dict({'text_height': .12}, **kw))
            assert rd(p) == rd(os.path.join(EXDIR, fn)), fn
    d0, d1 = read_dxf(os.path.join(EXDIR, 'DXF_KEKLS_cell.dxf')), read_dxf(os.path.join(EXDIR, 'DXF_KEKLS_cell_v2_remove.dxf'))
    nt = lambda d: sum(e['type'] == 'TEXT' for e in d['entities'])
    assert nt(d1) == 48 < nt(d0) == 52
    assert nt(read_dxf(os.path.join(EXDIR, 'DXF_KEKLS_cell_v2_spec_remove_h30.dxf'))) == 24


# ---------------------------------------------------------------- v2.4.1: 3D (磁石の立体外形。mode='3d'、Mode->"3D"、parts='split')
LAT_STR = r"""MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0);
QUAD  QF = (L=0.3 K1=0.7);
SEXT  SF = (L=0.2 K2=0);
BEND  B0 = (L=0.5 ANGLE=0);
LINE T = (D1 B0 D1 QF D1 SF D1);
"""
SPEC3 = {'BEND': {'Width': .5, 'PoleWidth': .12, 'Type': 'Htype', 'Height': .3}, 'QUAD': {'Width': .36, 'Height': .3}}


def _solids(path):
    from pysad.dxf_lattice import check_solids
    d = read_dxf(path)
    return d, check_solids(d['entities']), [e for e in d['entities'] if e['type'] == '3DFACE']


def _area_signed(P):
    return sum(P[i][0] * P[(i + 1) % 4][1] - P[(i + 1) % 4][0] * P[i][1] for i in range(4)) / 2.0


def _area(P):
    return abs(_area_signed(P))


def test_3d_plan_default_unchanged(tmp_path):
    """mode を与えない出力・mode='2d'・'plan'・parts='split' (平面図では無視) は全行同じで、3DFACE も POLE レイヤーもない"""
    ffs = _ffs()
    fs = [str(tmp_path / n) for n in ('a.dxf', 'b.dxf', 'c.dxf', 'd.dxf')]
    write_dxf(ffs, fs[0])
    write_dxf(ffs, fs[1], mode='2d')
    write_dxf(ffs, fs[2], mode='plan')
    write_dxf(ffs, fs[3], mode='2D', parts='split')
    txt = [open(f).read() for f in fs]
    assert txt[0] == txt[1] == txt[2] == txt[3]
    assert '3DFACE' not in txt[0] and 'POLE' not in txt[0] and 'plan view of the lattice' in txt[0]


def test_3d_body_closed_solids_match_plan(tmp_path):
    """mode='3d' (body): 磁石 1 台 = 3DFACE 6 枚の閉じた立体 (各有向辺が 1 回・逆向きも 1 回)、体積 = 2D の面積 × Height > 0 (法線が外向き)、
    上面・下面の (X, Y) は 2D の外形と同じ値、z = ±Height/2。基準軌道・要素名・MARK は平面図と同じ。レイヤー・色は 2D と同じ"""
    ffs = _ffs()
    f2, f3 = str(tmp_path / 'p.dxf'), str(tmp_path / 's.dxf')
    write_dxf(ffs, f2)
    info = write_dxf(ffs, f3, mode='3d')
    d2 = read_dxf(f2)
    d3, sol, faces = _solids(f3)
    assert d3['acadver'] == 'AC1009' and d3['eof'] and d3['layers'] == d2['layers'] and 'POLE' not in d3['layers']
    plan = [e for e in d2['entities'] if e['type'] == 'POLYLINE' and e['layer'] != 'TRAJECTORY']
    assert info['nsolids'] == len(sol) == len(plan) == 192 and len(faces) == 1152 and info['counts']['3DFACE'] == 1152
    assert not [e for e in d3['entities'] if e['type'] == 'POLYLINE' and e['closed']]           # 閉じた POLYLINE は 3DFACE に置き換わる
    for k, (s, pl) in enumerate(zip(sol, plan)):
        assert s['closed'] and s['nfaces'] == 6 and s['nverts'] == 8 and s['layer'] == pl['layer']
        top, bot = faces[6 * k]['pts'], faces[6 * k + 1]['pts']
        assert [(p[0], p[1]) for p in top] == [tuple(p) for p in pl['pts']]
        assert [(p[0], p[1]) for p in bot] == [tuple(p) for p in pl['pts'][::-1]]
        h = 0.406 if pl['layer'] == 'BEND' else {'QUAD': .45, 'SEXT': .3, 'OCT': .3, 'CAVI': .4, 'SOL': .4}[pl['layer']]
        assert s['bbox'][4] == pytest.approx(-h / 2) and s['bbox'][5] == pytest.approx(h / 2)
        assert s['volume'] == pytest.approx(_area(pl['pts']) * h, rel=1e-9)
    rest = lambda d: [(e['type'], e['layer'], e['pts'], e.get('text'), e.get('at')) for e in d['entities']
                      if not (e['type'] == '3DFACE' or (e['type'] == 'POLYLINE' and e['layer'] != 'TRAJECTORY'))]
    assert rest(d2) == rest(d3)


def test_3d_wedge_bend_faces_vertical(tmp_path):
    """偏向磁石は入口・出口の面が E1・E2 に傾いた楔形。側面・端面は鉛直 (上面と下面の (X, Y) が同じ)、4 つの側面は上下の辺が同じ (X, Y)、
    入口の面 (上面の第 1・2 頂点) の向きは E1·ANGLE (2D の試験と同じ式)"""
    t = r"""MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0);
BEND BA = (L=1.0 ANGLE=0.3 E1=0 E2=0) BB = (L=1.0 ANGLE=0.3 E1=0.8 E2=0.2) BC = (L=1.0 ANGLE=-0.2 E1=0.5 E2=0.5);
LINE T = (D1 BA D1 BB D1 BC D1);
"""
    ffs = _ffs(t, 'T')
    f = str(tmp_path / 'w.dxf')
    write_dxf(ffs, f, mode='3d', names=False, spec={'BEND': {'Width': .4, 'Height': .25}})
    _d, sol, faces = _solids(f)
    assert len(sol) == 3 and all(s['closed'] and s['volume'] > 0 for s in sol)
    for k in range(3):
        fc = faces[6 * k:6 * k + 6]
        top, bot = fc[0]['pts'], fc[1]['pts']
        assert [(p[0], p[1]) for p in top] == [(p[0], p[1]) for p in bot[::-1]]
        for side in fc[2:]:
            p = side['pts']
            assert p[0][2] == p[1][2] and p[2][2] == p[3][2] and p[2][2] - p[0][2] == pytest.approx(.25)
            assert (p[0][0], p[0][1]) == (p[3][0], p[3][1]) and (p[1][0], p[1][1]) == (p[2][0], p[2][1])
    for k, (name, e1r, ang_s) in enumerate((('BA', 0.0, 0.3), ('BB', 0.8, 0.3), ('BC', 0.5, -0.2))):
        i = [j for j in range(1, ffs.nlat) if ffs.line_value('ELEMENT', j) == name][0]
        gi = ffs.line_value('GEO', i)
        top = faces[6 * k]['pts']
        v = (top[0][0] - top[1][0], top[0][1] - top[1][1])
        ang_in = (math.atan2(v[1], v[0]) - (gi[1][0] + math.pi / 2) + math.pi) % (2 * math.pi) - math.pi
        assert abs(ang_in - e1r * ang_s) < 1e-6, name


def test_3d_default_and_magnet_spec_height(tmp_path):
    """Height: Spec の要素名 > Spec の型名 > MagnetSpec > 既定 (BEND 0.406、他は Width)。MagnetSpec の Height を読む"""
    ffs = _ms_ffs()

    def zext(**kw):
        f = str(tmp_path / 'h.dxf')
        write_dxf(ffs, f, mode='3d', names=False, **kw)
        _d, sol, _fc = _solids(f)
        return [round(s['bbox'][5] - s['bbox'][4], 6) for s in sol], [s['layer'] for s in sol]

    def by(h, lay, name):
        return [hh for hh, ll in zip(h, lay) if ll == name]
    h, lay = zext()
    assert set(by(h, lay, 'QUAD')) == {0.2}                 # MagnetSpec["QQ"] の Height
    assert set(by(h, lay, 'BEND')) == {0.4}                 # MagnetSpec["BX"]
    assert set(by(h, lay, 'SEXT')) == {0.22, 0.3}           # SF: Width .22 (Height なし → Width)、SD: 既定 0.30
    assert set(by(h, lay, 'OCT')) == {0.3} and set(by(h, lay, 'CAVI')) == {0.4}
    h, lay = zext(spec={'QUAD': {'Height': .1}, 'QF': {'Height': .15}, 'B': {'Height': .35}})
    q = by(h, lay, 'QUAD')
    assert q[0] == 0.15 and q[1] == 0.1 and set(q) == {0.15, 0.1}
    assert set(by(h, lay, 'BEND')) == {0.35}
    h, lay = zext(magnet_spec=False)
    assert set(by(h, lay, 'BEND')) == {0.406} and set(by(h, lay, 'QUAD')) == {0.45}
    with pytest.raises(ValueError):
        write_dxf(ffs, str(tmp_path / 'e.dxf'), mode='3d', spec={'QUAD': {'Height': -.1}})


def test_3d_strip_equals_polygon():
    """3D の帯の 4 角形 `_bend_strip3d` は、全幅 (−wr … wl) なら 2D の外形 `_bend_polygon` と同じ (ランダムな角度・端面・物理長・Ctype)"""
    import random
    from pysad.dxf_lattice import _bend_polygon, _bend_strip3d
    rng = random.Random(3)
    for _ in range(200):
        ang = rng.uniform(-0.5, 0.5)
        e1, e2 = rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)
        lc0, w, pw = rng.uniform(.3, 2.0), rng.uniform(.3, .8), rng.uniform(.05, .25)
        typ = rng.choice(['Htype', 'Ctype1', 'Ctype2'])
        phys = rng.choice([None, lc0 * rng.uniform(.5, 1.0)])
        poly, lc, ds = _bend_polygon({'PhysicalLength': phys}, ang, lc0, e1, e2, pw, w, typ)
        if typ == 'Ctype1':
            wr, wl = (pw / 2, w - pw / 2) if ang >= 0 else (w - pw / 2, pw / 2)
        elif typ == 'Ctype2':
            wr, wl = (pw / 2, w - pw / 2) if ang < 0 else (w - pw / 2, pw / 2)
        else:
            wr = wl = w / 2
        q = _bend_strip3d(ang, lc0, phys, e1, e2, -wr, wl)
        for a, b in zip(poly, q):
            assert abs(a[0] - b[0]) < 1e-12 and abs(a[1] - b[1]) < 1e-12


def _overlap_vol(a, b):
    d = 1.0
    for k in (0, 2, 4):
        d *= max(0.0, min(a[k + 1], b[k + 1]) - max(a[k], b[k]))
    return d


def test_3d_split_volumes_and_nonoverlap(tmp_path):
    """parts='split' (直線の偏向磁石 B0 と四極 QF。軸に平行なので直方体の体積を解析的に出せる): 各部品は閉じた立体で体積が正、
    ヨーク・ポールの体積は解析式と一致、互いに重ならない (AABB の重なり体積 0)、全体 (body) の bbox に入る。
    ポールはレイヤー POLE (色 9)、ヨークは型のレイヤー。SEXT は分けない (body 1 個)"""
    ffs = _ffs(LAT_STR, 'T')
    f = str(tmp_path / 'sp.dxf')
    info = write_dxf(ffs, f, mode='3d', parts='split', names=False, spec=SPEC3)
    d, sol, _fc = _solids(f)
    assert d['layers']['POLE'] == 9 and all(s['closed'] and s['volume'] > 0 for s in sol)
    parts = [(e, p, l) for (e, p, l, _h) in info['solids']]
    assert [x[1] for x in parts if x[0] == 'B0'] == ['yoke', 'yoke', 'yoke', 'yoke', 'pole', 'pole']       # 板 2 + 側脚 2 (H 型) + ポール 2
    assert [x[1] for x in parts if x[0] == 'QF'] == ['yoke'] * 4 + ['pole'] * 4
    assert [x[1] for x in parts if x[0] == 'SF'] == ['body']
    assert [x[2] for x in parts if x[1] == 'pole'] == ['POLE'] * 6 and [x[2] for x in parts if x[0] == 'B0' and x[1] == 'yoke'] == ['BEND'] * 4
    assert [x[2] for x in parts if x[0] == 'QF' and x[1] == 'yoke'] == ['QUAD'] * 4
    # 偏向磁石 B0: L=0.5、Width .5、Height .3、PoleWidth .12、Gap .04 (既定)、PoleHeight = (.3-.04)/4 = .065、LegWidth = (.5-.12)/4 = .095
    L, W, H, pw, g = .5, .5, .3, .12, .04
    ph, lw = (H - g) / 4, (W - pw) / 4
    zi = g / 2 + ph
    b = [s for (e, _p, _l), s in zip(parts, sol) if e == 'B0']
    yoke = L * W * (H / 2 - zi) * 2 + 2 * L * lw * 2 * zi
    pole = 2 * L * pw * ph
    assert sum(s['volume'] for s in b[:4]) == pytest.approx(yoke, rel=1e-6) and sum(s['volume'] for s in b[4:]) == pytest.approx(pole, rel=1e-6)
    # 四極 QF: Width .36、Height .3 → m = .3、Aperture = .12、LegWidth = (.3-.12)/4 = .045、PoleWidth = Aperture
    Wq, Hq, Lq = .36, .3, .3
    ap = .4 * min(Wq, Hq)
    lq = (min(Wq, Hq) - ap) / 4
    zq, vq = Hq / 2 - lq, Wq / 2 - lq
    q = [s for (e, _p, _l), s in zip(parts, sol) if e == 'QF']
    yoke_q = Lq * (2 * Wq * lq + 2 * lq * 2 * zq)
    pole_q = Lq * (2 * ap * (zq - ap / 2) + 2 * ap * (vq - ap / 2))
    assert sum(s['volume'] for s in q[:4]) == pytest.approx(yoke_q, rel=1e-6) and sum(s['volume'] for s in q[4:]) == pytest.approx(pole_q, rel=1e-6)
    for grp in (b, q):
        for i in range(len(grp)):
            for j in range(i + 1, len(grp)):
                assert _overlap_vol(grp[i]['bbox'], grp[j]['bbox']) < 1e-12, (i, j)
    fb = str(tmp_path / 'bd.dxf')
    write_dxf(ffs, fb, mode='3d', names=False, spec=SPEC3)
    _d2, body, _f2 = _solids(fb)
    bb = {'B0': body[0]['bbox'], 'QF': body[1]['bbox']}
    for (e, _p, _l), s in zip(parts, sol):
        if e in bb:
            assert all(s['bbox'][2 * k] >= bb[e][2 * k] - 1e-6 and s['bbox'][2 * k + 1] <= bb[e][2 * k + 1] + 1e-6 for k in range(3)), e


def test_3d_split_ctype_bend_and_errors(tmp_path):
    """C 型の偏向磁石 (Ctype1・ANGLE>0) の側脚は広い側 (左) の 1 本だけ。曲がった偏向磁石の split の部品も閉じて体積が正で、2D の外形の 4 角形の中に入る。
    入らない寸法 (Gap が大きい・PoleWidth が Width より大きい・LegWidth が大きい・PoleHeight が大きい・Aperture が大きい) と不正な mode・parts は ValueError"""
    ffs = _ffs(LAT_RECT, 'T')
    f = str(tmp_path / 'c.dxf')
    sp = {'BEND': {'Width': .5, 'PoleWidth': .12, 'Type': 'Ctype1', 'Height': .3}}
    info = write_dxf(ffs, f, mode='3d', parts='split', names=False, spec=sp)
    _d, sol, faces = _solids(f)
    assert all(s['closed'] and s['volume'] > 0 for s in sol)
    ps = [p for (e, p, _l, _h) in info['solids'] if e == 'BR']
    assert ps == ['yoke', 'yoke', 'yoke', 'pole', 'pole']                     # 板 2 + 側脚 1 + ポール 2
    plan_f = str(tmp_path / 'p.dxf')
    write_dxf(ffs, plan_f, spec=sp, names=False)
    P = [e for e in read_dxf(plan_f)['entities'] if e['layer'] == 'BEND'][0]['pts']
    pp = P if _area_signed(P) > 0 else P[::-1]

    def inside(pt):
        return all((pp[(i + 1) % 4][0] - pp[i][0]) * (pt[1] - pp[i][1]) - (pp[(i + 1) % 4][1] - pp[i][1]) * (pt[0] - pp[i][0]) >= -2e-6 for i in range(4))
    for fc in faces[:6 * len(ps)]:
        for p in fc['pts']:
            assert inside(p), p
    for bad in ({'Gap': .4}, {'PoleWidth': .6}, {'LegWidth': .5}, {'PoleHeight': .2}):
        with pytest.raises(ValueError):
            write_dxf(ffs, str(tmp_path / 'x.dxf'), mode='3d', parts='split', spec={'BEND': dict(sp['BEND'], **bad)})
    with pytest.raises(ValueError):
        write_dxf(_ffs(LAT_STR, 'T'), str(tmp_path / 'x.dxf'), mode='3d', parts='split', spec={'QUAD': {'Width': .3, 'Aperture': .3}})
    with pytest.raises(ValueError):
        write_dxf(ffs, str(tmp_path / 'x.dxf'), mode='solid')
    with pytest.raises(ValueError):
        write_dxf(ffs, str(tmp_path / 'x.dxf'), mode='3d', parts='thin')


def test_check_solids_detects_open_solid(tmp_path):
    """検査関数自身の確認: 面を 1 枚欠くと closed が偽、向きを 1 枚反転しても偽 (閉じて向きのそろった立体だけ真)"""
    from pysad.dxf_lattice import check_solids
    f = str(tmp_path / 's.dxf')
    write_dxf(_ffs(LAT_STR, 'T'), f, mode='3d', names=False)
    ents = [e for e in read_dxf(f)['entities'] if e['type'] == '3DFACE'][:6]
    assert check_solids(ents)[0]['closed']
    assert not check_solids(ents[:5], faces_per_solid=5)[0]['closed']
    flipped = [dict(e) for e in ents]
    flipped[2] = dict(flipped[2], pts=flipped[2]['pts'][::-1])
    assert not check_solids(flipped)[0]['closed']


def test_script_3d_matches_python(tmp_path):
    """SADScript の DrawMagnetsDXF[..., Mode->"3D"] (pysad で走らせる) と Python 関数の出力が一致する (body・split・MagnetSpec の Height・scale・重なり消し)。
    不正な Mode・Parts は画面にメッセージを出す。Mode なし = Mode->"2D" = 従来の平面図"""
    from pysad.ffscmd import run_sad_string
    spec_s = 'Spec->{"B"->{Width->.5, PoleWidth->.15, Type->Ctype1, Height->.3}, "QUAD"->{Width->.36, Height->.3}}'
    spec_p = {'B': {'Width': .5, 'PoleWidth': .15, 'Type': 'Ctype1', 'Height': .3}, 'QUAD': {'Width': .36, 'Height': .3}}
    cases = [('a.dxf', ', Mode->"3D"', {'mode': '3d'}),
             ('b.dxf', ', Mode->"3D", Parts->"split", ' + spec_s + ', TextHeight->.2',
              {'mode': '3d', 'parts': 'split', 'spec': spec_p, 'text_height': .2}),
             ('c.dxf', ', Mode->"3D", Parts->"split", Region->{"QF.2","B.5"}, Scale->1000, Names->False',
              {'mode': '3d', 'parts': 'split', 'region': ('QF.2', 'B.5'), 'scale': 1000, 'names': False}),
             ('d.dxf', ', Mode->"3D", Parts->"split", Overlap->"remove", TextHeight->.4', {'mode': '3d', 'parts': 'split', 'overlap': 'remove', 'text_height': .4}),
             ('e.dxf', ', Mode->"2D"', {})]
    prog = LAT + 'FFS;\nUSE RING;\n' + MS_DEFS + 'Get["%s"];\n' % os.path.join(EXDIR, 'DrawMagnetsDXF.n').replace('\\', '/')
    for fn, opts, kw in cases:
        out = str(tmp_path / fn).replace('\\', '/')
        buf = io.StringIO()
        ffs, _ = run_sad_string(prog + 'DrawMagnetsDXF["%s"%s];\nExit;\n' % (out, opts), out=buf)
        assert '???' not in buf.getvalue(), buf.getvalue()[-500:]
        py = str(tmp_path / ('py_' + fn))
        write_dxf(ffs, py, **kw)
        a = open(out, newline='').read().replace('\r\n', '\n')
        b = open(py, newline='').read().replace('\r\n', '\n')
        assert a == b, fn
    assert '3DFACE' in open(str(tmp_path / 'a.dxf')).read() and 'POLE' in open(str(tmp_path / 'b.dxf')).read() and 'POLE' not in open(str(tmp_path / 'a.dxf')).read()
    out0 = str(tmp_path / 'f.dxf').replace('\\', '/')
    run_sad_string(prog + 'DrawMagnetsDXF["%s"];\nExit;\n' % out0, out=io.StringIO())
    assert open(out0, newline='').read() == open(str(tmp_path / 'e.dxf'), newline='').read()
    buf = io.StringIO()
    run_sad_string(prog + 'DrawMagnetsDXF["%s", Mode->"4D", Parts->"thin"];\nExit;\n' % out0, out=buf)
    assert 'Mode must be' in buf.getvalue() and 'Parts must be' in buf.getvalue()


def test_example_3d_files():
    """example/dxf_lattice の 3D の DXF (FODO・KEK-LS の標準セル) を Python 関数で作り直して全行一致する。全立体が閉じて体積が正。
    body の頂点は同じ引数の平面図の外形と一致する"""
    need = ['DXF_FODO_cell_3d.dxf', 'DXF_FODO_cell_3d_split.dxf', 'DXF_FODO_ring_3d.dxf', 'DXF_KEKLS_cell_3d.dxf', 'DXF_KEKLS_cell_3d_split.dxf']
    lat = os.path.abspath(os.path.join(EXDIR, '..', 'kek_ls', 'lattice', 'kek_ls_v03_68.sad'))
    fodo = os.path.abspath(os.path.join(EXDIR, '..', '00_lattice', 'FODO_ring.sad'))
    if not os.path.exists(lat) or not os.path.exists(fodo) or not all(os.path.exists(os.path.join(EXDIR, n)) for n in need):
        pytest.skip('ラティスまたは 3D の DXF なし')
    import tempfile
    spec = {'B': {'Width': .5, 'PoleWidth': .15, 'Type': 'Ctype1', 'Height': .3}, 'QUAD': {'Width': .36, 'Height': .36}, 'SEXT': {'Width': .24, 'Height': .24}}
    ffs_f = load_ffs(fodo, 'RING')
    ffs_k = load_ffs(open(lat, encoding='utf-8').read() + '\nFFS;\nGet["MagnetSpec_KEKLS.n"];\n', 'ARCA', cwd=EXDIR)
    rd = lambda p: open(p, newline='').read().replace('\r\n', '\n')
    cases = (('DXF_FODO_cell_3d.dxf', ffs_f, {'region': ('SF.1', 'SD.1'), 'spec': spec, 'mode': '3d'}),
             ('DXF_FODO_cell_3d_split.dxf', ffs_f, {'region': ('SF.1', 'SD.1'), 'spec': spec, 'mode': '3d', 'parts': 'split'}),
             ('DXF_FODO_ring_3d.dxf', ffs_f, {'mode': '3d'}),
             ('DXF_KEKLS_cell_3d.dxf', ffs_k, {'text_height': .12, 'mode': '3d'}),
             ('DXF_KEKLS_cell_3d_split.dxf', ffs_k, {'text_height': .12, 'mode': '3d', 'parts': 'split'}))
    with tempfile.TemporaryDirectory() as td:
        for fn, ffs, kw in cases:
            p = os.path.join(td, fn)
            write_dxf(ffs, p, **kw)
            assert rd(p) == rd(os.path.join(EXDIR, fn)), fn
            _d, sol, faces = _solids(os.path.join(EXDIR, fn))
            assert sol and all(s['closed'] and s['volume'] > 0 for s in sol), fn
            if kw.get('parts') != 'split':
                pp = os.path.join(td, 'plan_' + fn)
                write_dxf(ffs, pp, **{k: v for k, v in kw.items() if k not in ('mode', 'parts')})
                plan = [e for e in read_dxf(pp)['entities'] if e['type'] == 'POLYLINE' and e['closed']]
                assert len(plan) == len(sol)
                for k, pl in enumerate(plan):
                    assert [(q[0], q[1]) for q in faces[6 * k]['pts']] == [tuple(q) for q in pl['pts']]


def test_plot_dxf3d_runs(tmp_path):
    """mplot3d の図 (example/dxf_lattice/plot_dxf3d.py) が PNG を作る (matplotlib がなければ skip)"""
    pytest.importorskip('matplotlib')
    import sys
    sys.path.insert(0, EXDIR)
    try:
        import plot_dxf3d
    finally:
        sys.path.remove(EXDIR)
    f, png = str(tmp_path / 's.dxf'), str(tmp_path / 's.png')
    write_dxf(_ffs(), f, mode='3d', parts='split', region=('QF.1', 'QF.3'), names=False)
    plot_dxf3d.plot3d(f, png, 't', panels=((25, -60, None), (30, -30, (-1, 8, -1, 4))), size=(6, 3), dpi=40)
    assert os.path.getsize(png) > 1000
