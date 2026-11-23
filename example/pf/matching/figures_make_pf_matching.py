# -*- coding: utf-8 -*-
"""PF のマッチングの図を作る。リポジトリの直下から: py -3.12 example/pf/matching/figures_make_pf_matching.py
  fig_FIT_injbump_PF.png  FFS_FIT_injbump_PF.sad の入射バンプ (GO のあとの DX [mm]、キッカー KK1〜KK4 と入射点 INJ の位置)"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "tools"))
from figlib import *

D = "example/pf/matching/"
ffs = run_file(D + "FFS_FIT_injbump_PF.sad")
n = draw_mpl.resolve_region(ffs, None, None)[1]
names = [str(ffs.line_value("NAME", k, 0.0)) for k in range(1, n + 1)]
s, dx = series(ffs, "DX", 0.1)
S = lambda nm: float(ffs.line_value("S", names.index(nm) + 1, 0.0))
print("line length", s[-1], "n", n)
fig, ax = plt.subplots(figsize=(7.0, 3.8))
ax.plot(s, [v * 1e3 for v in dx], color=BLUE, label="DX (バンプ軌道)")
for nm, c in (("KK1", RED), ("KK2", RED), ("KK3", RED), ("KK4", RED)):
    ax.axvline(S(nm), color=c, ls=":", lw=1.0)
    ax.text(S(nm), 16.5, nm, color=c, fontsize=8, ha="center")
ax.plot([S("INJ")], [15], "o", color=ORANGE, label="INJ (目標 DX = 15 mm, DPX = 1 mrad)")
ax.set_xlabel("入射部分の s [m]"); ax.set_ylabel("DX [mm]"); ax.set_ylim(top=19)
ax.set_title("PF 入射バンプ (キッカー 4 台、2.5 GeV)"); ax.legend(fontsize=8, loc="lower right")
save(fig, D + "fig_FIT_injbump_PF.png")
