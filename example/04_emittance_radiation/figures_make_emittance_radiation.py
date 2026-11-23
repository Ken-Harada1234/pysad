# -*- coding: utf-8 -*-
"""04_emittance_radiation の図を作る。リポジトリの直下から: py -3.12 example/04_emittance_radiation/figures_make_emittance_radiation.py
  fig_FFS_RADINT_H.png   RADINT.sad のリング (RINGQ) の 1 セルの分散関数 H(s) と I5 の被積分関数 H/ρ^3 (偏向磁石の中だけ)"""
import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from figlib import *

D = "example/04_emittance_radiation/"
ffs = run_file(D + "FFS_RADINT.sad")
i1, i2 = 3, 12                      # RINGQ の最初のセル (QF D2 B D2 QD DS SD D1 SF DS)
dl = 0.01
s, rb = series(ffs, "BX", dl, (i1, i2)); _, ax_ = series(ffs, "AX", dl, (i1, i2))
_, eta = series(ffs, "EX", dl, (i1, i2)); _, etap = series(ffs, "EPX", dl, (i1, i2))
beta = [x * x for x in rb]
H = [(e * e + (b * ep + a * e) ** 2) / b for b, a, e, ep in zip(beta, ax_, eta, etap)]
# 偏向磁石 (B) の範囲
rho = 1.2 / 0.2617993877991494
sb = float(ffs.line_value("S", i1 + 2, 0.0)); eb = sb + 1.2
inb = [sb <= x <= eb for x in s]
print("H in bend: min %.4g max %.4g m" % (min(h for h, f in zip(H, inb) if f), max(h for h, f in zip(H, inb) if f)))
# I5 (リング全体 = 24 セル、各セルの偏向磁石 1 台) を同じ点列から台形則で
I5 = 0.0
for k in range(len(s) - 1):
    if inb[k] and inb[k + 1]:
        I5 += 0.5 * (H[k] + H[k + 1]) / rho ** 3 * (s[k + 1] - s[k])
print("I5 (24 cells, trapezoid) = %.6g  (RADINT I5U = 0.01736053)" % (24 * I5))
fig, ax = plt.subplots(figsize=(6.5, 3.9))
ax.plot(s, [h * 1e3 for h in H], color=BLUE, label="H(s) (全域)")
ax.axvspan(sb, eb, color=ORANGE, alpha=0.2, label="偏向磁石 B (ρ = 4.58 m)")
ax.set_xlabel("1 セルの始点からの距離 s [m]"); ax.set_ylabel("H = (η² + (βη' + αη)²) / β  [mm]")
ax.set_title("FODO_ring (RINGQ) 1 セルの H 関数"); ax.legend(fontsize=8); ax.set_ylim(bottom=0)
save(fig, D + "fig_FFS_RADINT_H.png")
