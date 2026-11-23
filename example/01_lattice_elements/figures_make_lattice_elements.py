# -*- coding: utf-8 -*-
"""01_lattice_elements の図を作る。リポジトリの直下から: py -3.12 example/01_lattice_elements/figures_make_lattice_elements.py
  fig_elem_QUAD_beta.png   elem_QUAD.sad の 1 本の四極を通る β (drift + Q + drift)
  fig_FODO_cell_optics.png FODO_ring の 1 セルの β・η とラティス (周期解)
  fig_FODO_ring_layout.png FODO_ring の設計軌道 (GX-GY) と偏向磁石の位置 (FFS_ORG の幾何)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from figlib import *

D = "example/01_lattice_elements/"
RING_TAIL = "\nFFS;\nUSE RING; CELL; CALC;\n"
# 1) elem_QUAD: START(β=10) - Q - FIN
ffs = run_file(D + "elem_QUAD.sad", tail="\nFFS;\nUSE T; INS; CALC;\n")
optics(ffs, "BX BY LAT", D + "fig_elem_QUAD_beta.png", "elem_QUAD: 四極 1 本 (L=0.3 m) を通る β", figsize=(6.5, 4.0), tunes=False)

# 2) FODO 1 セル
ffs = run_file("example/00_lattice/FODO_ring.sad", tail=RING_TAIL)
optics(ffs, "BX BY & EX LAT", D + "fig_FODO_cell_optics.png", "FODO_ring 1 セル (3 GeV): β と分散", region=(3, 13), figsize=(6.5, 4.4))

# 3) リング配置 (FFS_ORG と同じ FODO_ring)
ffs = run_file("example/00_lattice/FODO_ring.sad", tail=RING_TAIL)
i1, i2 = draw_mpl.resolve_region(ffs, None, None)
gx, gz, bend = [], [], []
for k in range(i1, i2 + 2):
    gx.append(float(ffs.line_value("GX", k, 0.0))); gz.append(float(ffs.line_value("GY", k, 0.0)))
    bend.append(str(ffs.line_value("TYPENAME", k, 0.0)) == "BEND" if k <= i2 else False)
fig, ax = plt.subplots(figsize=(5.2, 5.0))
ax.plot(gx, gz, "-", color=GRAY, lw=1.0, label="設計軌道")
bx = [(gz[i], gx[i]) for i in range(len(bend)) if bend[i]]
ax.plot([p[1] for p in bx], [p[0] for p in bx], ".", color=RED, ms=5, label="偏向磁石の入口")
ax.plot(gx[0], gz[0], "s", color=BLUE, label="s = 0")
ax.set_aspect("equal"); ax.set_xlabel("GX [m]"); ax.set_ylabel("GY [m]")
ax.set_title("FODO_ring (24 セル, 周長 105.6 m) の設計軌道"); ax.legend(loc="center", fontsize=8)
save(fig, D + "fig_FODO_ring_layout.png")
