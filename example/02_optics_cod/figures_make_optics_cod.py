# -*- coding: utf-8 -*-
"""02_optics_cod の図を作る。リポジトリの直下から: py -3.12 example/02_optics_cod/figures_make_optics_cod.py
  fig_flag_CELL_INS_beta.png   CELL (周期解) と INS (先頭の BX=10) の βx の 1 周の比較 (flag_CELL_INS.sad)
  fig_FFS_BEAMSIZE_sigma.png   EMIT 後のビームサイズ σx・σy の 1 セル (FFS_BEAMSIZE.sad の RING)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from figlib import *

D = "example/02_optics_cod/"
HEAD = "FFS;\nUSE RING; CELL; CALC;\n"
base = "example/00_lattice/FODO_ring.sad"

# 1) CELL と INS
f1 = run_file(base, tail="\n" + HEAD)
s1, b1 = series(f1, "BX", 0.1); b1 = [v * v for v in b1]   # series は √β を返す
f2 = run_file(base, tail="\n" + HEAD + "INS;\nBXI 10; AXI 0; BYI 10; AYI 0; EXI 0; EPXI 0;\nCALC;\n")
s2, b2 = series(f2, "BX", 0.1); b2 = [v * v for v in b2]
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(s1, b1, color=BLUE, lw=1.2, label="CELL: 周期解 (先頭 βx = 5.03 m)")
ax.plot(s2, b2, color=ORANGE, lw=1.2, ls="--", label="INS: 先頭 BXI=10, AXI=0")
ax.set_xlabel("s [m]"); ax.set_ylabel("βx [m]"); ax.set_xlim(0, 105.6); ax.set_ylim(bottom=0)
ax.set_title("FODO_ring: CELL (周期解) と INS (初期値) の βx"); ax.legend(fontsize=8, loc="upper left")
save(fig, D + "fig_flag_CELL_INS_beta.png")

# 2) ビームサイズ (EMIT 後)
f3 = run_file(base, tail="\nFFS;\nUSE RING; CELL; CALC;\nRF VOLT 6.0E6;\nCALC;\nRFSW; RAD; NOFLUC;\nEMIT;\nMINCOUP = 0.01;\nBEAMSIZE;\n")
i1, i2 = 3, 13
sx = series(f3, "SIGX", 0.02, (i1, i2)); sy = series(f3, "SIGY", 0.02, (i1, i2))
fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(sx[0], [v * 1e3 for v in sx[1]], color=BLUE, label="σx")
ax.plot(sy[0], [v * 1e3 for v in sy[1]], color=ORANGE, ls="--", label="σy")
ax.set_xlabel("1 セルの始点からの距離 s [m]"); ax.set_ylabel("ビームサイズ σ [mm]")
ax.set_title("FODO_ring 1 セルのビームサイズ (6 MV, RAD ON, EMIT 後)"); ax.legend(fontsize=8)
save(fig, D + "fig_FFS_BEAMSIZE_sigma.png")
