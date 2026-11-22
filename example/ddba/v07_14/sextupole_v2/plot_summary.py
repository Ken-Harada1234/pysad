# -*- coding: utf-8 -*-
"""最終の比較の図: 候補ごと・誤差の大きさごとの DA 面積 (シードごとの点と平均)。pysad (二分あり、_c の一式) と gs (格子 9 光線 x 96 点、二分なし) の 2 枚組。
  py -3.12 plot_summary.py results/fig_da_compare.png [名前 ...]"""
import glob
import json
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
out = sys.argv[1]
names = sys.argv[2:] or ['final', 'r1a', 'r6a', 'r6b']
py = defaultdict(list)
for f in glob.glob(os.path.join(V.OUT, 'eval_*_c.json')):
    d = json.load(open(f, encoding='utf-8'))
    py[(d['name'], round(100 * d['rms']))].append(d['da']['area_mm2'])
gs = defaultdict(list)
for f in glob.glob(os.path.join(V.OUT, 'gs_*.json')):
    d = json.load(open(f, encoding='utf-8'))
    gs[(d['name'], round(100 * d['rms']))].append(d['gs']['area_mm2'])
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
cols = {'final': '#444444', 'r1a': '#d95f02', 'r6a': '#1b9e77', 'r6b': '#7570b3', 'r9c': '#e7298a'}
for ax, data, title in ((axs[0], py, 'pysad 追跡 (二分あり、512 周)'), (axs[1], gs, 'gs TrackParticles (格子 9×96 点、512 周)')):
    for j, lv in enumerate((0, 3, 10)):
        for i, n in enumerate(names):
            v = data.get((n, lv), [])
            if not v:
                continue
            x = j + (i - (len(names) - 1) / 2) * 0.18
            ax.scatter([x] * len(v), v, s=14, color=cols.get(n, 'k'), alpha=0.6)
            ax.plot([x - 0.07, x + 0.07], [np.mean(v)] * 2, color=cols.get(n, 'k'), lw=2.5, label=n if j == 0 else None)
            ax.plot([x - 0.04, x + 0.04], [np.min(v)] * 2, color=cols.get(n, 'k'), lw=1, ls=':')
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['誤差なし', 'K2L 誤差 3% (5 シード)', 'K2L 誤差 10% (5 シード)'])
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
axs[0].set_ylabel('x-y の DA 面積 [mm²] (横線 = 平均、点線 = 最悪)')
axs[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out, dpi=110)
