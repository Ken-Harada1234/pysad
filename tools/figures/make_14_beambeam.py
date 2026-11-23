# -*- coding: utf-8 -*-
"""docs/14: BEAMBEAM (弱-強、薄いレンズ) を FODO リングの先頭に 1 つ置いたときの、振幅によるチューンシフト (追跡)。
実行: py -3.12 tools/figures/make_14_beambeam.py"""
from _common import *
import numpy as np
SRC = (ROOT / "example/00_lattice/bb_ring.sad").read_text(encoding="utf-8") if (ROOT / "example/00_lattice/bb_ring.sad").exists() else None
ffs, I = run(SRC + "\nFFS;\nNORFSW; NORAD; NOFLUC;\nUSE RINGB; CELL; CALC;\n")
bx, by = val(I, 'Twiss["BX","ORIGIN"]'), val(I, 'Twiss["BY","ORIGIN"]')
sig = 1e-4         # BB の強ビームの σ = sqrt(EMIX*BX) = 1e-4 m
amps = np.array([0.05, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5])
def tune(plane, a):
    # 3 周の座標から cos(mu) を出す (docs/example と同じ方法)
    ffs_do(I, "a = %r; b = {{a},{0},{0},{0},{0},{0},{1}}; xs = {};"
              "Do[ xs = Append[xs, b[[1,1]]]; b = TrackParticles[{1, b}, 1, n, n][[2]], {n, 1, 3}];" % float(a))
    xs = get(I, "xs"); return xs
res = []
for a in amps:
    ffs_do(I, "a = %r; b = {{a},{0},{0},{0},{0},{0},{1}}; xs = {};"
              "Do[ xs = Append[xs, b[[1,1]]]; b = TrackParticles[{1, b}, 1, n, n][[2]], {n, 1, 3}];" % float(a * sig))
    x = get(I, "xs"); res.append(1 - np.arccos(np.clip((x[2] + x[0]) / (2 * x[1]), -1, 1)) / (2 * np.pi))
res = np.array(res)
print(bx, by, res)
fig, ax = plt.subplots(figsize=(5.2, 3.5))
ax.plot(amps, res, "o-", color=BLUE)
ax.axhline(0.696329, color=GRAY, ls="--", lw=1, label="BB なしの νx の小数部 (0.6963)")
ax.set_xlabel("水平振幅 [σ] (強ビームの σ = 0.1 mm)"); ax.set_ylabel("νx の小数部 (追跡)")
ax.set_title("ビーム・ビームによる振幅依存チューンシフト", fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); save(fig, "14_beambeam_tune_vs_amplitude.png")
