# -*- coding: utf-8 -*-
"""保存した周波数マップ (results/<tag>_fmap.npz) のチューンの足跡 (νx, νy) を D で色分けし、近い共鳴線 (次数 ≤ order) を重ねる。
  py -3.12 plot_fmap.py OUT.png TAG [TAG ...] [--order 6]   (TAG は eval_final_e10_s1 など。_fmap.npz は自動で付ける)
黒い丸 = 失われた点 (チューンが取れた分)、色 = D (log10 の拡散、青 = 安定、赤 = カオス)。"""
import argparse
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.family'] = ['Yu Gothic', 'MS Gothic', 'sans-serif']
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

ap = argparse.ArgumentParser()
ap.add_argument('out')
ap.add_argument('tags', nargs='+')
ap.add_argument('--order', type=int, default=6)
ap.add_argument('--qint', default='47,16')
args = ap.parse_args()
qi = [int(v) for v in args.qint.split(',')]
n = len(args.tags)
nc = min(n, 3)
nr = (n + nc - 1) // nc
fig, axs = plt.subplots(nr, nc, figsize=(5.2 * nc, 4.8 * nr), squeeze=False)
for ax, tag in zip(axs.ravel(), args.tags):
    r = V.Saved(os.path.join(V.OUT, tag + '_fmap.npz'))
    fx, fy = r.nux1, r.nuy1
    fin = np.isfinite(fx) & np.isfinite(fy)
    lost = r.lost_turn > 0
    ok = fin & ~lost & np.isfinite(r.D)
    sc = ax.scatter(fx[ok], fy[ok], c=r.D[ok], cmap='jet', vmin=-7, vmax=-2, s=14)
    ax.scatter(fx[fin & lost], fy[fin & lost], c='k', s=10, marker='x')
    x0, x1 = np.nanmin(fx[fin]) - 0.01, np.nanmax(fx[fin]) + 0.01
    y0, y1 = np.nanmin(fy[fin]) - 0.01, np.nanmax(fy[fin]) + 0.01
    for (m, nn) in V.lines_of_order(args.order):
        for p in range(-80, 400):
            # 整数部を含む共鳴 m(47+fx) + nn(16+fy) = p
            if nn == 0:
                xv = (p - nn * qi[1]) / m - qi[0]
                if x0 < xv < x1:
                    ax.axvline(xv, color='gray', lw=1.2 if m <= 2 else 0.5)
            elif m == 0:
                yv = p / nn - qi[1]
                if y0 < yv < y1:
                    ax.axhline(yv, color='gray', lw=1.2 if abs(nn) <= 2 else 0.5)
            else:
                xs = np.array([x0, x1])
                ys = (p - m * (qi[0] + xs)) / nn - qi[1]
                if min(ys) < y1 and max(ys) > y0 and (p - m * (qi[0] + x0)) / nn - qi[1] > y0 - 1 and abs(m) + abs(nn) <= args.order:
                    ax.plot(xs, ys, color='gray', lw=1.0 if abs(m) + abs(nn) <= 4 else 0.4)
                    if y0 < ys[0] < y1:
                        ax.text(xs[0], ys[0], '(%d,%d)' % (m, nn), fontsize=6, color='k')
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_title(tag, fontsize=9)
    ax.set_xlabel('νx 端数 (整数部 %d)' % qi[0])
    ax.set_ylabel('νy 端数 (整数部 %d)' % qi[1])
    fig.colorbar(sc, ax=ax, label='D')
fig.tight_layout()
fig.savefig(args.out, dpi=110)
