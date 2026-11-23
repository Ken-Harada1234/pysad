# -*- coding: utf-8 -*-
"""03_matching の図を作る。リポジトリの直下から: py -3.12 example/03_matching/figures_make_matching.py
  fig_FFS_FIT_before_after.png        FFS_FIT.sad: 位相進み 90 度のマッチング前後の β (1 セル)
"""
import sys, pathlib, re
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from figlib import *

D = "example/03_matching/"


def before_after(sad, marker, title, out, region=None, labels=("マッチング前", "マッチング後"), xl="s [m]"):
    """sad の marker (最初の FREE の行) より前までを「前」、全体を「後」として √β ではなく β を比べる。"""
    src = (ROOT / sad).read_text(encoding="utf-8")
    src = re.sub(r"(?im)^\s*(Exit|Abort)\s*;.*$", "", src)
    k = src.index(marker)
    fb = run_text(src[:k]); fa = run_text(src)
    fig, axs = plt.subplots(2, 1, figsize=(6.5, 5.0), sharex=True)
    info = []
    for ax, name, c in zip(axs, ("BX", "BY"), (BLUE, ORANGE)):
        for f, lab, ls, col in ((fb, labels[0], "--", GRAY), (fa, labels[1], "-", c)):
            s, v = series(f, name, 0.02, region or (None, None))
            ax.plot(s, [x * x for x in v], ls, color=col, label=lab)
            info.append((name, lab, min(x * x for x in v), max(x * x for x in v)))
        ax.set_ylabel("β%s [m]" % name[1].lower()); ax.legend(fontsize=8, loc="upper right"); ax.set_ylim(bottom=0)
    axs[1].set_xlabel(xl); axs[0].set_title(title)
    save(fig, out)
    for i in info:
        print(i)


before_after("example/03_matching/FFS_FIT.sad", "FREE QFH", "FFS_FIT: 位相進み 90° へのマッチング (1 セル、3 GeV)", D + "fig_FFS_FIT_before_after.png")
