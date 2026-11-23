# -*- coding: utf-8 -*-
"""ASTRID2 のラティスの図を作る。リポジトリの直下から: py -3.12 example/astrid/lattice/figures_make_astrid_lattice.py
  fig_ASTRID2_cell_optics.png   標準セル NCELL (7.6174 m) の β・η とラティス (ASTRID2_v01_01.sad、周期解)
  fig_ASTRID2_ring_beta.png     リング全周 (45.7044 m) の βy。MPW なし (RING) と MPW あり (RINGMPW) の比較 (ASTRID2_v01_02_id.sad)
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "tools"))
from figlib import *

D = "example/astrid/lattice/"
L1 = D + "ASTRID2_v01_01.sad"; L2 = D + "ASTRID2_v01_02_id.sad"

# 1) セルの光学 (リングを CELL モードで解いた周期解のうち最初のセル)
ffs = run_file(L1, tail="\nFFS;\nUSE RING; CELL; CALC;\n")
i1, i2 = draw_mpl.resolve_region(ffs, None, None)
names = [str(ffs.line_value("NAME", k, 0.0)) for k in range(i1, i2 + 1)]
ke = [k for k in range(i1, i2 + 1) if str(ffs.line_value("TYPENAME", k, 0.0)) == "MARK"]
print("marks", [(k, names[k - i1]) for k in ke][:8])
c_end = float(ffs.line_value("S", ke[1], 0.0)) if len(ke) > 1 else 7.6174
optics(ffs, "BX BY & EX LAT", D + "fig_ASTRID2_cell_optics.png", "ASTRID2 (0.58 GeV) の標準セルの光学", region=(ke[0], ke[1]), figsize=(6.5, 4.4))

# 2) MPW の有無でのリング全周の βy
fa = run_file(L2, tail="\nFFS;\nUSE RING; CELL; CALC;\n")
fb = run_file(L2, tail="\nFFS;\nUSE RINGMPW; CELL; CALC;\n")
sa, ba = series(fa, "BY", 0.05); sb, bb = series(fb, "BY", 0.05)
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(sa, [v * v for v in ba], color=GRAY, ls="--", label="MPW なし (RING)")
ax.plot(sb, [v * v for v in bb], color=BLUE, label="MPW あり (RINGMPW)")
ax.set_xlabel("s [m]"); ax.set_ylabel("βy [m]"); ax.set_ylim(bottom=0)
ax.set_title("ASTRID2: MPW 入りリングの βy"); ax.legend(fontsize=8)
save(fig, D + "fig_ASTRID2_ring_beta.png")
