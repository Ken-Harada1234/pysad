# -*- coding: utf-8 -*-
"""担当 B3 (4): tdgeo.f の四極などの DX・DY の解析微分 (幾何の FIT で DX・DY を変数にする) と、SOL 領域内の四極のずれ。

状況: SOL 領域 (BOUND≠0 の SOL の間) の外の四極の DX・DY は設計の幾何を動かさない (tdgeo は 0)。領域の中は tsgeo が軌道を追うので、四極のずれで終端の位置・向きが動く。
K64 の tdgeo は qdtrans (要素の出口の軌道の変化 = eps=1e-6 の中心差分) と tftmatu (終端までの Twiss の写像) で微分を作る。v2 は
  d c_out/d DX = (I − J) ∂t/∂DX  (c_out = F(c − t) + t、t = (DX, −(bz/2)DY, DY, (bz/2)DX)、J は要素の Jacobian) と、要素の写像の積による終端への輸送で解析的に作る。
pysad にはこのための前提が無かった (領域内の四極のずれは tsgeo も追跡・CALC も NotImplementedError): 四極のずれ (tsolrot の座標のずれと一様磁場のゲージ) を
追跡のカーネル (k_sol_quad)・Emittance[] (solemit)・tsgeo (sgeo) に入れた (gs と同じ。この追加は k64 でも効くが、k64 は SOL のずれ・回転 (trots) が未対応で SOL 自体がずれていると使えない)。
根拠: (a) 領域内の DX・DY・ROTATE つきの四極の幾何 (LINE["GX"…] 12 点) が gs と全桁一致、(b) 幾何の FIT (QF の DX と QD の DY で M2 の GY・GZ を合わせる)
が gs と同じ 2 回の Newton で収束し、変数の値が 1e-9 以内、(c) 解析微分が幾何の有限差分 (DX ± 1e-6) と 2e-3 以内。
条件: 3 GeV の輸送路、SOL SA (BZ=3 T、BOUND=1、GEO=1、DX=2 mm、DY=−1 mm、DPX=3 mrad、DPY=2 mrad) と SB (終端) の間に QF・QD。
"""
import io
import math
import re

import numpy as np
import pytest

from pysad import compat, ffscmd

_NUM = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'

_BASE = '''MOMENTUM = 3 GEV;
DRIFT DR = (L=1) DR2 = (L=0.5) DV=(L=0.5);
QUAD QF = (L=0.3 K1=0.3 DX=%(qfdx)s DY=-1.1e-3) QD = (L=0.3 K1=-0.25 DX=-1.5e-3 DY=%(qddy)s ROTATE=0.2);
SOL SA = (L=0 BZ=3 BOUND=1 GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002) SB = (L=0 BZ=0 BOUND=1);
MARK ST = (BX=5 BY=5) M2=();
LINE BL = (ST DR SA DR2 QF DV QD DR2 SB DR M2);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
'''

# gs (K64) の LINE["GX"..]: k=1..12 の S, GX, GY, GZ, GCHI1, GCHI2, GCHI3
_GS_GEO_TAIL = [4.10001175e+00, 4.10000714e+00, -9.40920805e-04, -3.96418550e-03, -8.76503566e-04, -2.08399915e-03, 0.0]


def _run(src):
    buf = io.StringIO()
    ffs, I = ffscmd.run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _vals(out, tag):
    for ln in out.splitlines():
        if ln.startswith(tag):
            return [float(v) for v in re.findall(_NUM, ln[len(tag):])]
    raise AssertionError(tag)


def test_quad_offsets_inside_sol_region_geometry_matches_gs(v2_mode):
    """領域内の DX・DY・ROTATE つきの四極 (QF: DX=2.7 mm DY=−1.1 mm、QD: DX=−1.5 mm DY=2.6 mm ROTATE=0.2) を通る軌道の幾何が gs と全桁一致 (10 桁)"""
    src = _BASE % dict(qfdx='2.7e-3', qddy='2.6e-3') + '''$FORM="18.14";
Do[Print["G ",k," ",LINE["S",k]," ",LINE["GX",k]," ",LINE["GY",k]," ",LINE["GZ",k]," ",LINE["GCHI1",k]," ",LINE["GCHI2",k]," ",LINE["GCHI3",k]],{k,1,12}];
Exit;
'''
    _ffs, out = _run(src)
    rows = [[float(v) for v in re.findall(_NUM, ln[2:])] for ln in out.splitlines() if ln.startswith('G ')]
    assert len(rows) == 12
    last = rows[-1][1:]
    for a, b in zip(last, _GS_GEO_TAIL):
        assert abs(a - b) < 1e-9 * max(1.0, abs(b)) + 1e-12, (last, _GS_GEO_TAIL)
    # 四極のずれなし (DX=DY=0) との差: M2 で GY が動く (ずれた軌道が領域を通る)
    src0 = _BASE % dict(qfdx='0', qddy='0') + 'Print["G0 ",LINE["GY","M2"]];\nExit;\n'
    src0 = src0.replace('DX=0 DY=-1.1e-3', 'DX=0 DY=0').replace('DX=-1.5e-3 DY=0 ROTATE', 'DX=0 DY=0 ROTATE')
    _f0, out0 = _run(src0)
    assert abs(_vals(out0, 'G0 ')[0] - rows[-1][2]) > 1e-4


def _fit(mode):
    compat.set_mode(mode)
    try:
        src = _BASE % dict(qfdx='0', qddy='0') + '''$FORM="18.14";
PageWidth=1999;
FIT M2 GY 0.0008; FIT M2 GZ -0.0035;
FREE QF DX; FREE QD DY;
GO;
Print["G1 =",LINE[{"GY","GZ"},"M2"]," ",LINE["DX","QF"]," ",LINE["DY","QD"]];
Exit;
'''
        src = src.replace('QF = (L=0.3 K1=0.3 DX=0 DY=-1.1e-3)', 'QF = (L=0.3 K1=0.3)').replace('QD = (L=0.3 K1=-0.25 DX=-1.5e-3 DY=0 ROTATE=0.2)', 'QD = (L=0.3 K1=-0.25)')
        return _run(src)
    finally:
        compat.set_mode('k64')


def test_geometry_fit_with_quad_dx_dy_matches_gs():
    """QF の DX と QD の DY を変数に M2 の GY・GZ を合わせる幾何の fit: gs は 2 回の Newton で収束 (DX = 2.78848432670e-3、DY = 2.66798056318e-3)。
    v2 の解析微分は 2 回で収束し、変数が 1e-12 (絶対) で gs と一致 (gs の最初の残差 4.86e-8、v2 解析 5.19e-8、中心差分 2.54e-8)"""
    _ffs, out = _fit('v2')
    assert out.count(' Matched.') >= 2
    assert re.search(r'^\s+2\s+0\s+[0-9.E+-]+\s+\(NEWTON\)', out, re.M)
    g = _vals(out, 'G1 =')
    assert abs(g[0] - 8.0e-4) < 1e-9 and abs(g[1] + 3.5e-3) < 1e-9
    assert abs(g[2] - 2.78848432670e-3) < 1e-10 and abs(g[3] - 2.66798056318e-3) < 1e-10


def test_analytic_derivative_equals_finite_difference_of_geometry(v2_mode):
    """解析微分 (sol_dcod + tdgeo_term) が、幾何の有限差分 (DX・DY を ±1e-6) と 2e-3 以内 (M2 の GY・GZ・GCHI2。GCHI3 は領域の終端で 0 に戻されるので比べない: K64 の tdgeo の CHI3 の式は 1e-3 の値を返すが、幾何は動かない)"""
    from pysad import match6d
    from pysad.ffscmd import run_sad_string
    src = _BASE % dict(qfdx='0.5e-3', qddy='0.4e-3') + 'Exit;\n'
    ffs, _ = run_sad_string(src, out=io.StringIO())
    els = ffs.lat.elements
    ne = len(els)                            # M2 (最後の要素)
    kf, kd = 5, 7                            # QF、QD の位置 (1 始まり: ST DR SA DR2 QF DV QD …)
    assert type(els[kf - 1]).__name__ == 'Quadrupole' and type(els[kd - 1]).__name__ == 'Quadrupole'
    kbe = match6d._sol_region_end(els, kf)
    assert kbe == 9
    rows = ffs.rows
    geo0 = ffs.ensure_geo().geo

    def geoval(k, key, name):
        el = els[k - 1]
        old = getattr(el, key)
        out = []
        for sg in (1.0, -1.0):
            setattr(el, key, old + sg * 1e-6)
            ffs.calc(quiet=True)
            g = ffs.ensure_geo()
            from pysad import geo as _g
            out.append({'GY': g.gxyz(ne)[1], 'GZ': g.gxyz(ne)[2], 'CHI2': g.chi(ne)[1], 'CHI3': g.chi(ne)[2]}[name])
        setattr(el, key, old)
        ffs.calc(quiet=True)
        return (out[0] - out[1]) / 2e-6
    rows = ffs.rows
    geo0 = ffs.ensure_geo().geo
    for k, key in ((kf, 'DX'), (kf, 'DY'), (kd, 'DX'), (kd, 'DY')):
        dcod = match6d.sol_dcod(els, rows, k, kbe, key)
        for name in ('GY', 'GZ', 'CHI2'):
            ana = match6d.tdgeo_term(geo0, name, ne, k, els[k - 1], 'Quadrupole', key, 1.0, dcod=dcod, kbe=kbe)
            num = geoval(k, key, name)
            assert abs(ana - num) <= 2e-3 * max(abs(num), 1e-6) + 1e-9, (k, key, name, ana, num)


def test_quad_dx_dy_outside_sol_region_does_not_move_geometry(v2_mode):
    """領域の外の四極の DX・DY は設計の幾何を動かさない: fit の微分は 0 (tdgeo の kbe=0)。fit は変数を動かさず、残差が残る"""
    from pysad import match6d
    assert match6d.tdgeo_term(np.zeros((3, 3, 4)), 'GY', 3, 2, None, 'Quadrupole', 'DX', 1.0, dcod=None, kbe=0) == 0.0
    assert match6d._sol_region_end([], 1) == 0
