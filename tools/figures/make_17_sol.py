# -*- coding: utf-8 -*-
"""docs/17: SOL 領域 (BZ=0.3 T、1.7 GeV) を含むリングの光学 (DRAW と同じ描画)。docs/17 の最小例と同じラティス。
実行: py -3.12 tools/figures/make_17_sol.py"""
from _common import *
from pysad import draw_mpl
SRC = '''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1) DR2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.3) QD = (L=0.3 K1=-0.25);
BEND B = (L=1 ANGLE=0.2);
SOL SA = (L=0 BZ=0.3 BOUND=1) SB = (L=0 BZ=0 BOUND=1);
MARK ST = ();
LINE BL = (ST QF DR B DR2 QD DR SA DR2 QF DR2 QD DR2 SB DR B DR);
FFS;
USE BL;
CELL; NORAD; NOFLUC;
CALC;
'''
ffs, I = run(SRC)
print(val(I, 'Twiss["NX","$$$"]/(2 Pi)'), val(I, 'Twiss["NY","$$$"]/(2 Pi)'))
draw_mpl.draw(ffs, "BX BY & EX EY LAT", save=str(OUT / "17_sol_ring_optics.png"), show=False,
              title="SOL 領域 (SA〜SB、BZ=0.3 T) を含むリングの光学", figsize=(6.5, 4.8), dpi=115)
