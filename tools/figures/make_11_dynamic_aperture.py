# -*- coding: utf-8 -*-
"""docs/11: PF 2.5 GeV リングの x-y ダイナミックアパーチャ (24 方向、1000 周、運動量 0 と ±2%)。
データは example/pf/tracking/fn_DynamicApertureSurvey_da_xy_PF.pysad.txt (pysad の出力)。新しい計算はしない。
実行: py -3.12 tools/figures/make_11_dynamic_aperture.py"""
import re
from _common import *
import numpy as np
txt = (ROOT / "example/pf/tracking/fn_DynamicApertureSurvey_da_xy_PF.pysad.txt").read_text(encoding="utf-8", errors="replace").splitlines()
blocks = {}
cur = None
for ln in txt:
    m = re.match(r"dp/p =\s*(\S+) \(in sigma_E\)", ln)
    if m:
        cur = float(m.group(1)); blocks[cur] = []; continue
    if cur is not None:
        p = ln.split()
        if len(p) == 2:
            try: blocks[cur].append((float(p[0]), float(p[1])))
            except ValueError: pass
fig, ax = plt.subplots(figsize=(5.6, 4.2))
for (d, pts), c, mk in zip(sorted(blocks.items(), key=lambda t: abs(t[0])), (BLUE, ORANGE, GREEN), ("o-", "s", "^")):
    a = np.array(pts); a = np.vstack([a, a[:1]]) if abs(d) < 1 else a
    ax.plot(a[:, 0], a[:, 1], mk, color=c, ms=3, lw=1.2, label="δ = %+.1f σE (%+.1f %%)" % (d, d * 0.01 / 27.4423 * 1.0) if False else "運動量のずれ %+.1f σE" % d)
ax.set_xlabel("x [σx]"); ax.set_ylabel("y [σy]"); ax.set_title("PF 2.5 GeV: x-y ダイナミックアパーチャ (1000 周)")
ax.set_aspect("equal"); ax.legend(fontsize=7, loc="upper right")
fig.tight_layout(); save(fig, "11_dynamic_aperture_pf_xy.png")
print({k: len(v) for k, v in blocks.items()})
