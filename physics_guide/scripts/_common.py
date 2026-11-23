# -*- coding: utf-8 -*-
"""物理ガイドの確認スクリプトの共通部。

リポジトリのルートから `py -3.12 physics_guide/scripts/quad_01_h2_exact.py` のように動かす。
結果は標準出力と `physics_guide/outputs/<スクリプト名>.out.txt` (UTF-8) に出る。
"""
import io
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUT_DIR = os.path.join(ROOT, 'physics_guide', 'outputs')
FIG_DIR = os.path.join(ROOT, 'physics_guide', 'figures')

try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass


class _Tee:
    def __init__(self, path):
        self.f = open(path, 'w', encoding='utf-8', newline='\n')
        self.o = sys.stdout

    def write(self, s):
        self.o.write(s)
        self.f.write(s)

    def flush(self):
        self.o.flush()
        self.f.flush()


def start(script_file):
    """標準出力を outputs/<名前>.out.txt にも書く。戻り値: 名前"""
    name = os.path.splitext(os.path.basename(script_file))[0]
    os.makedirs(OUT_DIR, exist_ok=True)
    sys.stdout = _Tee(os.path.join(OUT_DIR, name + '.out.txt'))
    return name


def font():
    """matplotlib の日本語フォント設定 (Windows の既定フォントから選ぶ)"""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager as fm
    names = {f.name for f in fm.fontManager.ttflist}
    for cand in ('Yu Gothic', 'Meiryo', 'MS Gothic', 'IPAexGothic', 'Noto Sans CJK JP'):
        if cand in names:
            matplotlib.rcParams['font.family'] = cand
            break
    matplotlib.rcParams['axes.unicode_minus'] = False
    return matplotlib


def set_mode(mode):
    """'k64' (K64 と同じ動作) か 'v2' (既定)。四極の積分の方式 (quadbody.scheme) もこれに従う"""
    from pysad import compat
    old = compat._mode
    compat.set_mode(mode)
    return old


def make_quad(L, K1, extra='', momentum='3GEV'):
    """QUAD 1 台だけのラティスを作り、(ffs, 要素) を返す。K1 は SAD の積分強度 K1 (= B'L/Bρ [1/m])。
    要素の属性 el.K1 は長さあたり K1/L [1/m^2] になる"""
    from pysad.ffscmd import run_sad_string
    src = ('MOMENTUM = %s;\nQUAD SA=(L=%.15g K1=%.15g %s);\nMARK ST = ();\nLINE ONE = (ST SA);\nFFS;\nUSE ONE;\n'
           % (momentum, L, K1, extra))
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs, [e for e in ffs.lat.elements if e.name == 'SA'][0]


def track_one(ffs, c):
    """(canonical な x, px, y, py, z, δ) を、ラティス全体 (ST SA) で追跡 (RFSW ON で z を保つ。RFSW OFF の追跡は 1 周ごとに z = 0 に戻す。放射なし)。戻り値: 6 成分"""
    import numpy as np
    from pysad import track6
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6.track_particles(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]
