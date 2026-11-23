# -*- coding: utf-8 -*-
"""docs/13: 縦ウェイクがバンチ内の粒子のエネルギーに与える変化 (TrackParticles、LWAKE の ON/OFF の差)。
docs/13 §2.1 のラインと同じ。実行: py -3.12 tools/figures/make_13_wake.py"""
from _common import *
import numpy as np
LAT = '''MOMENTUM = 1.7GEV;
DRIFT DR = (L=1);
QUAD QF = (L=0.3 K1=0.6 DX=2e-4) QD = (L=0.3 K1=-0.5 DY=-1e-4);
BEND B = (L=1 ANGLE=0.1 E1=0.05 E2=0.05);
MARK ST = () WK2 = ();
LINE BL = (ST QF DR B DR QD DR WK2 DR B DR QF DR QD DR QF);
FFS;
USE BL;
TRPT; INS; NORAD; NOFLUC;
PBUNCH = 1E12; NBUNCH = 1;
WakeFunction[Longitudinal, "QF.1"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};
WakeFunction[Longitudinal, "WK2"]  = {{-0.2, 1e9}, {0, 3e9}, {0.2, -1e9}, {0.3, 0}};
n = 201;
zz = Table[-0.1 + 0.2 (i-1)/(n-1), {i, n}];
zero = Table[0, {i, n}];
'''
ffs, I = run(LAT + "NOLWAKE;\nr0 = TrackParticles[{1, {zero, zero, zero, zero, zz, zero, Table[1, {i, n}]}}, 17, 1];\n"
                   "LWAKE;\nr1 = TrackParticles[{1, {zero, zero, zero, zero, zz, zero, Table[1, {i, n}]}}, 17, 1];\n"
                   "d0 = r0[[2,6]]; d1 = r1[[2,6]]; zf = r1[[2,5]];\n")
z = get(I, "zz"); d0, d1 = get(I, "d0"), get(I, "d1")
print(np.abs(d0).max(), np.abs(d1).max())
fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.5))
zt = np.array([-0.2, 0, 0.2, 0.3]); wt = np.array([1, 3, -1, 0]) * 1e9
zq = np.array([-0.1, -0.02, 0, 0.02, 0.1, 0.2]); wq = np.array([2, 1, 0.5, -1, -2, -1]) * 1e9
ax[0].plot(zq, wq / 1e9, "o-", color=BLUE, label="QF.1")
ax[0].plot(zt, wt / 1e9, "s-", color=ORANGE, label="WK2")
ax[0].set_xlabel("z [m]"); ax[0].set_ylabel("縦ウェイク関数 W [GV/C]"); ax[0].legend(fontsize=8)
ax[0].set_title("入力のウェイク関数 (表)", fontsize=9)
ax[1].plot(z * 1e3, (d1 - d0) * 1e6, color=BLUE)
ax[1].set_xlabel("初期の z [mm]"); ax[1].set_ylabel("δ の変化 (LWAKE ON − OFF) [1e-6]")
ax[1].set_title("追跡後の δ のウェイクによる変化", fontsize=9)
fig.tight_layout(); save(fig, "13_wake_longitudinal_kick.png")
