# -*- coding: utf-8 -*-
"""KEK-LS (v03_68) の COD 応答の図を作る。リポジトリの直下から: py -3.12 example/kek_ls/errors/figures_make_kekls_errors.py  (約 20 秒)
  fig_KEKLS_cod_response.png  水平ステアリング HD.1 に K0 = 1 urad を与えたときの閉軌道 (SAD の COD) と、
                              応答行列の式 -K0 * 0.5 sqrt(βm βc) cos(πν − |Δμ|) / sin(πν) のモニター位置の値 (KEKLS_cod_response.sad)"""
import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "tools"))
from figlib import *

D = "example/kek_ls/errors/"
src = (ROOT / D / "KEKLS_cod_response.sad").read_text(encoding="utf-8")
head = src[:src.index("Responsex = Table")]                    # ラティス + USE RING; cell calc emit; USE RING NOEXPAND
ffs = run_text(head + '\nLINE["K0","HD.1"] = 1E-6;\nFFS["cell calc noexpand"];\n')
n = draw_mpl.resolve_region(ffs, None, None)[1]
names = [str(ffs.line_value("NAME", k, 0.0)) for k in range(1, n + 1)]
mon = [k for k, nm in enumerate(names, 1) if nm.startswith("PM.")]
hd = [k for k, nm in enumerate(names, 1) if nm.startswith("HD.")]
print("monitors", len(mon), "steerings HD", len(hd))
kc = [k for k, nm in enumerate(names, 1) if nm == "HD.1"][0]
S = lambda k: float(ffs.line_value("S", k, 0.0))
tv = lambda key, k: float(ffs.twiss_value(key, k, 0.0))
nux = tv("NX", n) / 2 / math.pi
sd, dx = series(ffs, "DX", 0.5)
cod = [tv("DX", k) * 1e6 for k in mon]
# 応答行列の式 (-K0 × Response)。β・位相は K0 を与えたあとの Twiss (結合・K0 の影響は小さい)
resp = []
for k in mon:
    r = 0.5 * math.sqrt(tv("BX", k) * tv("BX", kc)) * math.cos(math.pi * nux - abs(tv("NX", k) - tv("NX", kc))) / math.sin(math.pi * nux)
    resp.append(-1e-6 * r * 1e6)
print("max|COD| um", max(abs(v) for v in cod), "max|COD - formula| um", max(abs(a - b) for a, b in zip(cod, resp)))
fig, ax = plt.subplots(figsize=(7.0, 3.8))
ax.plot(sd, [v * 1e6 for v in dx], color=BLUE, lw=1.0, label="SAD の COD (DX)")
ax.plot([S(k) for k in mon], resp, ".", color=ORANGE, ms=4, label="式 −K0·Response (モニター PM)")
ax.axvline(S(kc), color=GRAY, ls=":", lw=1.0); ax.text(S(kc) + 3, ax.get_ylim()[1] * 0.85, "HD.1 (K0 = 1 μrad)", fontsize=8)
ax.set_xlabel("s [m]"); ax.set_ylabel("閉軌道 x [μm]"); ax.set_xlim(0, 570.72)
ax.set_title("KEK-LS: ステアリング 1 個の COD 応答"); ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
save(fig, D + "fig_KEKLS_cod_response.png")
