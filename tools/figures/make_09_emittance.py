# -*- coding: utf-8 -*-
"""docs/09: 放射減衰と量子励起のつり合い。RF 4 MV の FODO リングで、広がりゼロのバンチを RAD+FLUC で追跡し、
σδ・σz が Emittance 関数の平衡値へ近づく様子。実行: py -3.12 tools/figures/make_09_emittance.py"""
import time
from _common import *
import _common
import numpy as np
_common.FODO = _common.FODO.replace("VOLT=1.5E6", "VOLT=4.0E6")
ffs, I = start("RFSW; RADCOD; RAD; FLUC;\ne = Emittance[];\nsige = MomentumSpread/.e; sigz = BunchLength/.e; emx = (Emittances/.e)[[1]];\n")
sige, sigz, emx = val(I, "sige"), val(I, "sigz"), val(I, "emx")
print(sige, sigz, emx)
n, nturn, step = 300, 6000, 250
ffs_do(I, "NPARA = 1; n = %d; c0 = Table[0, {i, n}]; c = {c0, c0, c0, c0, c0, c0, Table[1, {i, n}]};" % n)
sd, sz, sx2, turns = [], [], [], []
t0 = time.time()
for k in range(nturn // step):
    ffs_do(I, "c = TrackParticles[{1, c}, Length[LINE[\"NAME\"]], %d, 1][[2]]; sx = Plus@@(c[[1]]^2)/n; sdd = Sqrt[Plus@@(c[[6]]^2)/n]; szz = Sqrt[Plus@@(c[[5]]^2)/n];" % step if False else
              "r = TrackParticles[{1, c}, 243, 1, %d]; c = r[[2]];" % step)
    c = get(I, "c")
    sd.append(c[5].std()); sz.append(c[4].std()); turns.append((k + 1) * step)
    if k == 0: print("time per block", time.time() - t0)
turns = np.array(turns)
fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.4))
ax[0].plot(turns, np.array(sd) * 1e3, "o-", ms=3, color=BLUE, label="追跡 (300 粒子)")
ax[0].axhline(sige * 1e3, color=ORANGE, ls="--", label="Emittance 関数の平衡値")
ax[0].set_ylabel("σδ [1e-3]"); ax[0].set_xlabel("周回数"); ax[0].legend(fontsize=7)
ax[1].plot(turns, np.array(sz) * 1e3, "o-", ms=3, color=BLUE)
ax[1].axhline(sigz * 1e3, color=ORANGE, ls="--")
ax[1].set_ylabel("σz [mm]"); ax[1].set_xlabel("周回数")
fig.suptitle("RAD+FLUC の追跡: バンチの広がりが平衡値へ", fontsize=9)
fig.tight_layout()
save(fig, "09_emittance_equilibrium_growth.png")
