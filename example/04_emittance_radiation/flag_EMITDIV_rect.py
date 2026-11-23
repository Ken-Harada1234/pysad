# -*- coding: utf-8 -*-
"""v2 の Emittance[] と RADINT の差 (矩形磁石、SP8 (大きなリング B) で 1.4e-3) の原因 = EMITDIV の既定 (1) の離散化の誤差 (k64・v2・扇形・矩形とも同じ)。
SP8 (大きなリング B、8 GeV、1436 m。SPring-8 型、デザインレポートより再現。ファイルは追跡外なのでパスは書かない) の数値 (2026-11-19): εx (RFSW、RADCOD) を EMITDIV = 1・2・3・4・8・16 で
  k64 の矩形 (E1 = E2 = 0.5): 6.955790e-9 6.947071e-9 6.946637e-9 6.946573e-9 6.946554e-9 6.946556e-9   (RADINT の 6.955384e-9 に既定で 1e-4: 離散化 +1.3e-3 と端の項が 2 倍の誤り −1.3e-3 が相殺)
  v2 の矩形:              6.965116e-9 6.955988e-9 6.955497e-9 6.955414e-9 6.955378e-9 6.955375e-9   (収束値は RADINT と 1.3e-6。既定は +1.4e-3)
  k64 の扇形 (E = 0):     6.877552e-9 6.868536e-9 6.868086e-9 6.868019e-9 6.867998e-9 6.868000e-9   (既定の誤差 +1.4e-3 は矩形と同じ)
ここでは、小さな 24 セルの FODO (3 GeV、矩形) の収束を示す。実行: py -3.12 example/04_emittance_radiation/flag_EMITDIV_rect.py (出力は flag_EMITDIV_rect.out.txt)"""
import io, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'pysad', 'tests'))
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
from pysad import ffscmd, compat, emit6d
import _a2_indep as A

SRC = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=%g E2=%g);
CAVI RF = (L=0 VOLT=4e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 24*CELL1 RF);
FFS;
USE RING; CELL; CALC;
'''
for e in (0.5, 0.0):
    compat.set_mode('k64')
    ffs, _ = ffscmd.run_sad_string(SRC % (e, e), out=io.StringIO())
    I, C, _nu = A.integrals_rich(ffs.lat.elements, n=80)
    d = A.derived(I, C, ffs.momentum())
    print('E1 = E2 = %g·ANGLE: 平面リングの放射積分 (RADINT と同じ式、端の項は 1 倍) Jx = %.5f' % (e, d['Jx']))
    for mode in ('k64', 'v2'):
        compat.set_mode(mode)
        row = []
        for div in (1, 2, 4, 16):
            emit6d.EMITDIV[:] = [div] * 3
            em = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
            row.append(em['emittance_x'])
        emit6d.EMITDIV[:] = [4, 4, 4]
        e0 = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=False)
        jx = -2 * e0['damping_rate'][0] * ffs.momentum() / e0['energy_loss_eV']
        print('   %-4s εx (RFSW) EMITDIV=1,2,4,16: %s  (1 の 16 に対する差 %+.1e)   RFSW OFF の Jx (EMITDIV=4) = %.5f' % (
            mode, ' '.join('%.6e' % r for r in row), row[0] / row[3] - 1, jx))
    emit6d.EMITDIV[:] = [1, 1, 1]
