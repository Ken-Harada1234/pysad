# -*- coding: utf-8 -*-
"""docs/16: FODO リング (3 GeV、RF 4 MV) のスピン追跡 (POL; RAD; NOFLUC)。水平面内に向けたスピンの累積の歳差角と、
SpinTune の直線からのずれ (シンクロトロン振動 δ(n) による)。実行: py -3.12 tools/figures/make_16_spin_pol.py"""
from _common import *
import _common
import numpy as np
_common.FODO = _common.FODO.replace("VOLT=1.5E6", "VOLT=4.0E6")
ffs, I = start("RFSW; RAD; NOFLUC; POL;\nem = Emittance[]; st = SpinTune/.em; nst = NominalSpinTune/.em; u0 = EnergyLossU0/.em;\n")
st, nst, u0 = val(I, "st"), val(I, "nst"), val(I, "u0")
z0 = float(np.arcsin(u0 / 4e6) / (2 * np.pi * 5e8 / 2.99792458e8))     # 放射損失とつり合う同期位相の z
N = 200
ffs_do(I, "b = {{0},{0},{0},{0},{%r},{0},{0},{0},{1},{1}}; dl = {}; ang = {};"
          "Do[ r = TrackParticles[{1, b}, 243, n, n]; b = r[[2]]; dl = Append[dl, b[[6,1]]]; ang = Append[ang, b[[7,1]]], {n, 1, %d}];" % (-z0, N))
dl, ang = get(I, "dl"), get(I, "ang")
n = np.arange(1, N + 1)
turns = -np.unwrap(ang) / (2 * np.pi)          # 累積の回転数 (向きは SpinTune の符号に合わせる)
print(st, nst, np.polyfit(n, -turns, 1)[0])
fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.5))
ax[0].plot(n, -turns - st * n, color=BLUE)
ax[0].set_xlabel("周回数"); ax[0].set_ylabel("累積の回転数 − SpinTune × 周回数")
ax[0].set_title("SpinTune の直線からのずれ", fontsize=9)
ax[1].plot(n, dl * 1e3, color=ORANGE)
ax[1].set_xlabel("周回数"); ax[1].set_ylabel("δ [1e-3]"); ax[1].set_title("エネルギーのずれ (シンクロトロン振動)", fontsize=9)
fig.tight_layout(); save(fig, "16_spin_pol_precession.png")
