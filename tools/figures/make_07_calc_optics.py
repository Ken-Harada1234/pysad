# -*- coding: utf-8 -*-
"""docs/07: FODO リングの 1 セルの β・η とラティス (DRAW と同じ描画)、チューンの運動量ずれ DP0 依存。
実行: py -3.12 tools/figures/make_07_calc_optics.py"""
from _common import *
from pysad import draw_mpl
import numpy as np

ffs, I = start()
fig0 = draw_mpl.draw(ffs, "BX BY & EX LAT", region=(3, 13), show=False,
                     title="FODO リング 1 セル (3 GeV)", figsize=(6.5, 4.6))
save(fig0, "07_calc_optics_fodo_cell.png")

# チューンの DP0 依存 (6 極なし / SF=+5, SD=-8)
dps = np.linspace(-0.01, 0.01, 11)
def scan():
    r = []
    for d in dps:
        ffs_do(I, "DP0 = %r; CALC;" % float(d))
        r.append((val(I, 'Twiss["NX","$$$"]/2/Pi'), val(I, 'Twiss["NY","$$$"]/2/Pi')))
    ffs_do(I, "DP0 = 0; CALC;")
    return np.array(r)
a = scan()
ffs_do(I, "SF K2 5; SD K2 -8; CALC;")
b = scan()
fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.6))
for k, (t, nm) in enumerate((("x", "水平"), ("y", "垂直"))):
    ax[k].plot(dps * 100, a[:, k], "o-", color=BLUE, label="6 極なし (SF=SD=0)")
    ax[k].plot(dps * 100, b[:, k], "s-", color=ORANGE, label="SF K2=5, SD K2=-8")
    ax[k].set_xlabel("運動量のずれ DP0 [%]"); ax[k].set_ylabel("ν%s" % t); ax[k].set_title(nm + "チューン")
ax[0].legend(fontsize=7)
fig.tight_layout()
save(fig, "07_calc_optics_tune_vs_dp.png")
