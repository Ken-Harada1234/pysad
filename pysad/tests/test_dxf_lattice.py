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
EXDIR = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'dxf_lattice'))

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
