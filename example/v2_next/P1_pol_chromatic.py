# -*- coding: utf-8 -*-
"""v2next P1 項目 1: 4.2 GeV・12 MV の 12 セル FODO (B3_pol_sideband.sad と同じ) の Polarization の脱偏極率に、ベータトロンのクロマティックな変調を入れる前後。
追跡 (FLUC つきスピン追跡、A2・B3、480 粒子) は (10.4 +- 1.6)e-8 /周。使い方: py -3.12 example/v2_next/P1_pol_chromatic.py [RF電圧 MV (既定 12)]  (約 1 分)。結果は P1_pol_chromatic.out.txt。"""
import io
import sys

from pysad import compat, emit6d, ffscmd, spin

compat.set_mode('v2')
V = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0
SRC = f'''MOMENTUM = 4.2GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD QFE = (L=0.3 K1=0.7 DY=3e-3) QDE = (L=0.3 K1=-0.6 DY=-3e-3);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT={V}e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELL2 = (QFE D2 B D2 QD D1);
LINE CELL3 = (QF D2 B D2 QDE D1);
LINE RING = (ST 2*CELL1 CELL2 4*CELL1 CELL3 4*CELL1 RFC);
FFS;
USE RING; CELL; CALC; RFSW; RAD; NOFLUC; POL;
Exit;
'''
ffs, _ = ffscmd.run_sad_string(SRC, out=io.StringIO())
for chrom in (False, True):
    spin.CHROMATIC = chrom
    sp = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)['spin']
    print('RF %g MV  クロマティックな変調 %-5s  脱偏極率 %.4e /周 (線形 %.4e)  P = %.5f  変調指数 (スピン %.2f, ベータトロン %.2f)' %
          (V, chrom, sp['depol_rate'], sp['depol_rate_linear'], sp['equpol'][0], sp['mod_spin'], sp['mod_chrom']))
print('参考 (12 MV): FLUC つきスピン追跡 (10.4 +- 1.6)e-8、線形の軌道の MC 5.0e-8、クロマティックな軌道の MC 1.24e-7 (3000 粒子 x 30000 周)')
