# -*- coding: utf-8 -*-
"""ASTRID2 の放射積分の図を作る。リポジトリの直下から: py -3.12 example/astrid/tune_radiation/figures_make_astrid_radiation.py
  fig_ASTRID2_H.png  標準セル (MKS.1 から MKE.1) の H(s) = (η² + (βη' + αη)²)/β と、偏向磁石の範囲 (I5 の被積分関数 H/|ρ|³ は磁石の中だけ)。
                     リング全体 (標準セル 6 個) の I5 も同じ点列から求めて calc_RadInt_ASTRID2.pysad.txt の値と比べる"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "tools"))
from figlib import *

D = "example/astrid/tune_radiation/"
ffs = run_file("example/astrid/lattice/ASTRID2_v01_01.sad", tail="\nFFS;\nUSE RING; CELL; CALC;\n")
k1, k2 = 1, 33                                   # MKS.1 .. MKE.1
dl = 0.005
s, rb = series(ffs, "BX", dl, (k1, k2)); _, a = series(ffs, "AX", dl, (k1, k2))
_, e = series(ffs, "EX", dl, (k1, k2)); _, ep = series(ffs, "EPX", dl, (k1, k2))
H = [(x * x + (b * b * y + c * x) ** 2) / (b * b) for b, c, x, y in zip(rb, a, e, ep)]
bends = []
for k in range(k1, k2 + 1):
    if str(ffs.line_value("TYPENAME", k, 0.0)) == "BEND" and abs(float(ffs.line_value("ANGLE", k, 0.0))) > 1e-9:
        bends.append((float(ffs.line_value("S", k, 0.0)), float(ffs.line_value("L", k, 0.0)), float(ffs.line_value("ANGLE", k, 0.0))))
print("bends (s, L, angle):", bends)
I5 = 0.0
for k in range(len(s) - 1):
    for (sb, L, th) in bends:
        if sb <= s[k] and s[k + 1] <= sb + L + 1e-9:
            I5 += 0.5 * (H[k] + H[k + 1]) * (s[k + 1] - s[k]) / (L / th) ** 3
print("I5 (6 cells) = %.6g   (calc_RadInt_ASTRID2.pysad.txt: 0.1367070724)" % (6 * I5))
fig, ax = plt.subplots(figsize=(6.5, 3.9))
ax.plot(s, [h * 1e3 for h in H], color=BLUE, label="H(s)")
for i, (sb, L, th) in enumerate(bends):
    ax.axvspan(sb, sb + L, color=ORANGE, alpha=0.2, label="偏向磁石 BH (ρ = %.3f m)" % (L / th) if i == 0 else None)
ax.set_xlabel("セルの始点 (MKS) からの距離 s [m]"); ax.set_ylabel("H = (η² + (βη' + αη)²) / β  [mm]")
ax.set_title("ASTRID2 (0.58 GeV) 標準セルの H 関数"); ax.legend(fontsize=8, loc="upper left"); ax.set_ylim(bottom=0)
save(fig, D + "fig_ASTRID2_H.png")
