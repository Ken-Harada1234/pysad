# -*- coding: utf-8 -*-
"""例の図のスクリプト (各フォルダの figures_make_*.py) が使う共通部分。

- SAD の台本 (.sad) を pysad で走らせて FFS の状態を得る (run_file・run_text)。
- `pysad.draw_mpl.draw` で光学とラティス配置を描く (optics)。点列は `series` で取り出せる。
- 保存 (save) は PNG を 1 枚 300 KB 以下に収める (dpi を下げて再保存する)。
実行は、リポジトリの直下から `py -3.12 example/<フォルダ>/figures_make_<名前>.py`。環境変数 PYSAD_SAD_PACKAGES が無ければ
D:/work/SAD-master/SAD-master/Packages を使う。
"""
import io
import os
import pathlib
import re
import sys

os.environ.setdefault("PYSAD_SAD_PACKAGES", "D:/work/SAD-master/SAD-master/Packages")
ROOT = pathlib.Path(__file__).resolve().parents[2]
EX = ROOT / "example"
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

# 色覚に配慮した配色 (Okabe-Ito)
BLUE, ORANGE, GREEN, RED, SKY, PURPLE, GRAY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#555555"
for _f in ("Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP"):
    if any(_f == x.name for x in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = _f
        break
plt.rcParams.update({"axes.unicode_minus": False, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 100})

from pysad.ffscmd import run_sad_string  # noqa: E402
from pysad import draw_mpl  # noqa: E402


def run_text(src):
    """台本の文字列を実行して FFS の状態を返す (画面表示は捨てる)。"""
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def run_file(path, tail="", head=""):
    """台本ファイル (リポジトリの直下からの相対パスか絶対パス) を実行する。head・tail は前後に足す文字列。"""
    p = pathlib.Path(path)
    if not p.is_absolute():
        p = ROOT / p
    src = p.read_text(encoding="utf-8", errors="replace")
    src = re.sub(r"(?im)^\s*(Exit|Abort)\s*;.*$", "", src)      # Exit; があると USE の状態が消えるので除く
    return run_text(head + src + tail)


def series(ffs, name, dl=0.05, region=(None, None)):
    """OpticsPlot と同じ点列 (s [m], 値)。name は BX・EX・SIGX など DRAW の関数名。
    BX・BY は OpticsPlot の単位どおり √β [√m] で返る (β が欲しいときは 2 乗する)。"""
    i1, i2 = draw_mpl.resolve_region(ffs, *region)
    return draw_mpl.optics_list(ffs, name, i1, i2, dl, 0.0)


def optics(ffs, spec, out, title="", region=None, figsize=(7.4, 4.4), tunes=True):
    """DRAW 相当 (`spec` は DRAW の引数。例 "BX BY & EX LAT") を PNG に保存する。"""
    fig = draw_mpl.draw(ffs, spec, region=region, title=title, show=False, info=False, figsize=figsize, tunes=tunes)
    save(fig, out)


def save(fig, out, dpi=110):
    """PNG に保存する。300 KB を超えたら dpi を下げる。out は相対パスなら呼び出し側のフォルダ基準でなく ROOT 基準。"""
    p = pathlib.Path(out)
    if not p.is_absolute():
        p = ROOT / p
    p.parent.mkdir(parents=True, exist_ok=True)
    while True:
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        if p.stat().st_size <= 300 * 1024 or dpi <= 60:
            break
        dpi -= 10
    plt.close(fig)
    print(p.relative_to(ROOT), p.stat().st_size // 1024, "KB", "dpi", dpi)
