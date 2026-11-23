# -*- coding: utf-8 -*-
"""図のスクリプトの共通部分 (環境・色・保存)。リポジトリの直下から `py -3.12 tools/figures/make_<名前>.py` で実行する。"""
import os, sys, pathlib
os.environ.setdefault("PYSAD_SAD_PACKAGES", "D:/work/SAD-master/SAD-master/Packages")
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
BLUE, ORANGE, GREEN, RED, GRAY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#666666"
for f in ("Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP"):
    from matplotlib import font_manager as fm
    if any(f == x.name for x in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({"axes.unicode_minus": False, "figure.dpi": 100, "axes.grid": True, "grid.alpha": 0.3})

import figsave   # 折れ線系は SVG、そうでなければ PNG (tools/figures/figsave.py)

def save(fig, name, vector=True):
    """name は .png で与える。SVG にできる図は .svg を書き、PNG は置かない (figsave.save_figure)。"""
    return figsave.save_figure(fig, OUT / name, dpi=115, vector=vector)

import io

def run(script):
    """SAD 台本の文字列を実行して (ffs, I) を返す (画面出力は捨てる)。FFS 部は `FFS;` の次の行から書く。"""
    from pysad.ffscmd import run_sad_string
    return run_sad_string(script, out=io.StringIO())

def ffs_do(I, text):
    """FFS の文を追加で実行する (出力は捨てる)。"""
    old = I.out
    I.out = io.StringIO()
    try:
        I.run(text, 'ffs')
    finally:
        I.out = old

def val(I, expr):
    """SAD の式の値 (数) を返す。"""
    ffs_do(I, "pysad$tmp = %s;" % expr)
    return float(I.globals["PYSAD$TMP"] if "PYSAD$TMP" in I.globals else I.globals["pysad$tmp"])

def start(body=""):
    return run(FODO + "\nFFS;\nUSE RING;\nCELL;\nCALC;\n" + body)

FODO = (ROOT / "example" / "00_lattice" / "FODO_ring.sad").read_text(encoding="utf-8")

def series(ffs, f, r1, r2):
    """位置 r1..r2 の間の関数 f (BX・EX など、素の値) の (s [m], 値) (DRAW と同じ点列)。s は r1 の位置から測る。"""
    import numpy as np
    from pysad import draw_mpl
    xs, ys, x0, leng = draw_mpl._optics_series(ffs, f, r1, r2)
    ys = np.array(ys)
    if f in ('BX', 'BY'):
        ys = ys ** 2          # DRAW は √β を描くので 2 乗して β [m] に戻す
    return np.array(xs) - x0, ys

def get(I, name):
    """SAD の変数を numpy 配列にして返す。"""
    import numpy as np
    return np.array(I.globals[name], dtype=float)
