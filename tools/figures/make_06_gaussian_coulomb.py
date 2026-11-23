# -*- coding: utf-8 -*-
"""docs/06: GaussianCoulomb[x, y, sigx, sigy] (ガウス分布の電荷が作る 2 次元の力) の形。丸いビームと平たいビーム。
SAD の台本から GaussianCoulomb を呼んで描く。実行: py -3.12 tools/figures/make_06_gaussian_coulomb.py"""
from _common import *
import numpy as np
ffs, I = run("FFS;\nx = 0;\n")
def force(x, y, sx, sy):
    ffs_do(I, "ff = GaussianCoulomb[%r, %r, %r, %r];" % (float(x), float(y), float(sx), float(sy)))
    v = I.globals["ff"]
    return complex(v)
r = np.linspace(0.02, 10, 160)
round_ = np.array([-force(a, 0, 1, 1).real for a in r])
fx = np.array([-force(a * 3, 0, 3, 1).real for a in r])        # 平たいビーム σx=3σy: x 軸上、a は x/σx
fy = np.array([-force(0, a, 3, 1).imag for a in r])            # y 軸上、a は y/σy
fig, ax = plt.subplots(figsize=(5.8, 3.7))
ax.plot(r, round_, color=BLUE, label="丸いビーム (σx = σy = 1)")
ax.plot(r, 2 * (1 - np.exp(-r ** 2 / 2)) / r, "k:", lw=1, label="解析式 2(1 − exp(−r²/2σ²))/r")
ax.plot(r, fx, color=ORANGE, label="平たいビーム (σx = 3σy): x 軸上 (横軸 x/σx)")
ax.plot(r, fy, color=GREEN, ls="--", label="平たいビーム (σx = 3σy): y 軸上 (横軸 y/σy)")
ax.set_xlabel("位置 [σ]"); ax.set_ylabel("力の大きさ |f| (σy = 1 を与えた値)")
ax.set_title("GaussianCoulomb: ガウス分布の電荷が作る力", fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); save(fig, "06_sadscript_library_gaussian_coulomb_force.png")
print(round_[:3], fx.max(), fy.max())
