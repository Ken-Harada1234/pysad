# -*- coding: utf-8 -*-
"""primer_*.py の共通部分: FODO リング (example/00_lattice/FODO_ring.sad、6 極 SF・SD の K2 = 4, -6 を入れる) を読む。
追跡は pysad の k64 モード (gs と同じ)、RF 空洞は止め、放射なし。"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
os.environ.setdefault('MPLBACKEND', 'Agg')
from pysad import compat
from pysad.parser import load_sad_string
compat.set_mode(os.environ.get('PYSAD_COMPAT', 'k64'))

P_FODO = 3e9      # 運動量 [eV/c]


def fodo(k2=(4.0, -6.0)):
    """FODO リングの Lattice と Twiss の表 (rows)、リング先頭の ((αx, βx), (αy, βy))、線形のチューン (νx, νy)、周長 [m]"""
    src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
    src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g)' % tuple(k2))
    lat = load_sad_string(src).build_lattice('RING')
    rows = lat.twiss_table()
    r0 = rows[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    nu = (float(rows[-1].mux) / (2 * math.pi), float(rows[-1].muy) / (2 * math.pi))
    return lat, rows, tw, nu, float(rows[-1].s)


def setup_plot():
    """matplotlib の日本語フォント (あるものを使う)"""
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    for nm in ('Yu Gothic', 'Meiryo', 'MS Gothic', 'IPAexGothic', 'Noto Sans CJK JP'):
        if nm in have:
            plt.rcParams['font.family'] = nm
            break
    plt.rcParams['axes.unicode_minus'] = False
    return plt
