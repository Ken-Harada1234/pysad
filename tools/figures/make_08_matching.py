# -*- coding: utf-8 -*-
"""docs/08: マッチング (FIT NX 6.30 NY 3.30; FREE QF QD; GO) の前後の β 関数 (FODO 1 セル)。
実行: py -3.12 tools/figures/make_08_matching.py"""
from _common import *
ffs, I = start()
s0, bx0 = series(ffs, "BX", 3, 13); _, by0 = series(ffs, "BY", 3, 13)
nx0, ny0 = val(I, 'Twiss["NX","$$$"]/2/Pi'), val(I, 'Twiss["NY","$$$"]/2/Pi')
k0 = (val(I, 'LINE["K1","QF"]'), val(I, 'LINE["K1","QD"]'))
ffs_do(I, "FIT NX 6.30 NY 3.30; FREE QF QD; GO;")
s1, bx1 = series(ffs, "BX", 3, 13); _, by1 = series(ffs, "BY", 3, 13)
nx1, ny1 = val(I, 'Twiss["NX","$$$"]/2/Pi'), val(I, 'Twiss["NY","$$$"]/2/Pi')
k1 = (val(I, 'LINE["K1","QF"]'), val(I, 'LINE["K1","QD"]'))
print(nx0, ny0, k0, nx1, ny1, k1)
fig, ax = plt.subplots(1, 2, figsize=(6.8, 3.6), sharey=True)
for a, (b0, b1, t) in zip(ax, ((bx0, bx1, "x"), (by0, by1, "y"))):
    a.plot(s0, b0, color=BLUE, label="マッチング前" if t == "x" else None)
    a.plot(s1, b1, color=ORANGE, ls="--", label="マッチング後" if t == "x" else None)
    a.set_xlabel("s [m] (セル 1)"); a.set_title("β%s" % t)
ax[0].set_ylabel("β [m]")
ax[0].legend(fontsize=8, loc="upper left")
ax[1].text(0.97, 0.04, "ν=(%.3f, %.3f) → (%.2f, %.2f)\nQF K1 %.3f→%.3f, QD K1 %.3f→%.3f" % (nx0, ny0, nx1, ny1, k0[0], k1[0], k0[1], k1[1]),
           transform=ax[1].transAxes, fontsize=6.5, va="bottom", ha="right")
fig.tight_layout()
save(fig, "08_matching_fodo_beta_before_after.png")
