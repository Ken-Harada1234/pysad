"""幾何の fit (GX/GY/GZ/CHI1-3)・DISP G/OG・ビームサイズ Twiss["SIG.."] のテスト (docs/history/work_geo_sig.md)。
参照値は SAD K64 (<work-perf>/geo_tiny1.sad.out, geo_tiny2.sad.out) から。"""
import io
import math
import os

import numpy as np
import pytest

from pysad.ffscmd import run_sad_string
from pysad.sadscript import SADError

TINY = '''
MOMENTUM = 1GEV;
DRIFT D1 = (L=1) D2 = (L=2) D3 = (L=1.5);
BEND B1 = (L=1 ANGLE=0.3 E1=0.15 E2=0.15) B2 = (L=1 ANGLE=-0.2 E1=-0.1 E2=-0.1) B3 = (L=0.5 ANGLE=0.1 ROTATE=90 DEG);
QUAD QF = (L=0.3 K1=1) QD = (L=0.3 K1=-1);
MARK M1 = (BX=5 BY=5) M2 = ();
LINE L1 = (M1 D1 QF D2 B1 D3 QD D1 B3 D2 B2 D1 M2);
FFS;
USE L1;
INS;
CALC;
$FORM = "20.15";
'''


def _run(src):
    out = io.StringIO()
    ffs, I = run_sad_string(src, out=out)
    return ffs, I, out.getvalue()


def test_disp_geo_matches_sad():
    ffs, I, text = _run(TINY + '''
disp geo;
disp ogeo;
''')
    lines = text.splitlines()
    # SAD の disp geo (<work-perf>/geo_tiny1.sad.out)
    assert ' Element              Gx         Gy           Gz            s        Length     Value        Chi1        Chi2        Chi3      #' in lines
    assert ' D3                    4.285067     .148878     .0         4.300000  1.500000  1.5000000   17.188734     .0          .0          6' in lines
    assert ' B2                    9.338010    1.711937    -.224646    9.600000  1.000000  -.2000000   17.188734   -5.729578 -1.7615E-17    11' in lines
    assert ' $$$                  11.302029    2.007195    -.421659   11.600000   .0      0             5.673576   -5.614997    1.141950    14' in lines
    assert ' Element             OGx        OGy          OGz            s        Length     Value       OChi1       OChi2       OChi3      #' in lines
    # LINE["GX".."GCHI3", $$$]
    v = [ffs.line_value(k, ffs.nlat) for k in ('GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3')]
    ref = [11.30202922541477, 2.007194890032480, -.421658592291833, .099022575150505, -.098000185923169, .019930782663914]
    np.testing.assert_allclose(v, ref, rtol=1e-12, atol=1e-15)


def test_geometry_fit_matches_sad():
    ffs, I, text = _run(TINY + '''
fit M2 GX 11.5;
fit M2 GZ -0.5;
free D2 L;
free D3 L;
go;
r1 = {LINE["GX","$$$"], LINE["GY","$$$"], LINE["GZ","$$$"], LINE["L","D2"], LINE["L","D3"]};
fit;
fit $$$ CHI1 6;
fit $$$ CHI3 1.2;
free B3 ANGLE;
free B2 ANGLE;
go;
r2 = {LINE["GCHI1","$$$"], LINE["GCHI2","$$$"], LINE["GCHI3","$$$"], LINE["ANGLE","B3"], LINE["ANGLE","B2"]};
''')
    r1 = [float(x) for x in I.get('r1')]
    r2 = [float(x) for x in I.get('r2')]
    # SAD (<work-perf>/geo_tiny2.sad.out)。微分は数値微分なので最後の桁は違う
    np.testing.assert_allclose(r1, [11.5, 1.825691685480795, -0.5, 2.784721292123146, .105017001818998], rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(r2, [.104719755119660, -.106098377968199, .020943951023932, .108138187840613, -.194168101296658],
                               rtol=1e-9, atol=1e-12)
    # 条件表に GX/GZ/CHI1/CHI3 の行が出る
    assert 'M2        f GX   11.5      1  11.500000' in text
    assert '$$$       f CHI1 6         1   6.000000' in text
    # 反復の表示は SAD と同じ列 (3 行目まで)
    assert '       3        0    29.47      (NEWTON)   6.2500E-02      4' in text
    assert '      12        0   4.5146E-08  (NEWTON)    1.000          4' in text


RING = '''
MOMENTUM = 1GEV;
DRIFT D = (L=0.5);
BEND B = (L=1 ANGLE=0.39269908169872414 E1=0.5 E2=0.5);
QUAD QF = (L=0.3 K1=1.2) QD = (L=0.3 K1=-1.1);
MARK MK = ();
LINE CELL = (MK QF D B D QD D B D);
LINE RING = (8*CELL);
FFS;
USE RING;
CELL;
CALC;
'''


def test_beam_size_requires_emit_and_uses_twiss():
    ffs, I, text = _run(RING)
    with pytest.raises(SADError):
        ffs.twiss_value('SIGX', ffs.nlat)
    ffs, I, text = _run(RING + '''
emit;
sx = Twiss["SIGX","$$$"]; spx = Twiss["SIGPX","$$$"]; sy = Twiss["SIGY","$$$"]; spy = Twiss["SIGPY","$$$"];
lx = LINE["SIGX","$$$"]; sxy = Twiss["SIGXY","$$$"]; sdp = Twiss["SIGDP","$$$"]; m = Twiss["SIG","$$$"];
sx2 = Twiss["SIGX","QF.2"]; sxf = Twiss["SIGX","QF.2+0.5"];
''')
    g = I.globals
    emx, emy, sige = g['EMITX'], g['EMITY'], g['SIGE']
    assert emx > 0 and sige > 0, text[-2000:]
    r = ffs.row(ffs.nlat)
    # 無結合のリング: σx² = βx εx + ηx² σE²、σy² = βy max(εy, 0.01(εx+εy)) (MINCOUP の下限)
    emyy = max(emy, 0.01 * (emx + emy))
    assert abs(float(I.get('sx')) - math.sqrt(r.betx * emx + r.etx**2 * sige**2)) < 1e-15
    assert abs(float(I.get('spx')) - math.sqrt((1 + r.alfx**2) / r.betx * emx + r.etpx**2 * sige**2)) < 1e-15
    assert abs(float(I.get('sy')) - math.sqrt(r.bety * emyy)) < 1e-15
    assert abs(float(I.get('spy')) - math.sqrt((1 + r.alfy**2) / r.bety * emyy)) < 1e-15
    assert float(I.get('lx')) == float(I.get('sx'))                     # LINE["SIGX"] は同じ値 (CALC6D でない)
    assert abs(float(I.get('sxy'))) < 1e-30                             # Σxy = ηx ηy σE² = 0
    assert abs(float(I.get('sdp')) - sige) < 1e-18
    m = I.get('m')
    assert len(m) == 6 and len(m[0]) == 6 and abs(float(m[0][0]) - float(I.get('sx'))**2) < 1e-25
    r2 = ffs.row(ffs.locate('QF.2') if hasattr(ffs, 'locate') else 10)
    assert float(I.get('sx2')) > 0 and float(I.get('sxf')) > 0 and float(I.get('sx2')) != float(I.get('sxf'))


def test_beamsize_command_needs_emit():
    ffs, I, text = _run(RING + 'beamsize;\n')
    assert 'EMIT' in text
    ffs, I, text = _run(RING + 'emit;\ncell calc beamsize;\n')
    assert 'BEAMSIZE' not in text or 'EMIT を先に' not in text
