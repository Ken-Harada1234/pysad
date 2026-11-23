# -*- coding: utf-8 -*-
"""DDBA の図を作る。リポジトリの直下から: py -3.12 example/ddba/figures_make_ddba.py
  fig_DDBA_28_cell_optics.png   DDBA_28 (5 GeV、周長 749.48 m) のリングの先頭 36.05 m (長直線部の中心 MKLS.1 から MKMC.3 まで) の β・η とラティス (周期解)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from figlib import *

D = "example/ddba/"
ffs = run_file(D + "FFS_cell_calc_emit_DDBA_28.sad")      # Exit; を除いた台本 (USE RING; cell calc emit;)
k1, k2 = 2, 175          # MKLS.1 (s = 0、長直線部の中心) から MKMC.3 (s = 36.05 m) まで。長いセルの半分 + 標準の短いセル 1 個
for nmf in ("BX", "BY", "EX"):
    s, v = series(ffs, nmf, 0.05, (k1, k2))
    if nmf != "EX":
        v = [x * x for x in v]
    print(nmf, min(v), max(v), "s_end", s[-1])
optics(ffs, "BX BY & EX LAT", D + "fig_DDBA_28_cell_optics.png", "DDBA_28 (5 GeV): 長直線部の中心から 36 m", region=(k1, k2), figsize=(7.4, 4.6))
