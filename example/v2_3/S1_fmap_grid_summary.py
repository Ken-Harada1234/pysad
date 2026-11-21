# -*- coding: utf-8 -*-
"""v2.3.0 S1 (f): 格子の細かさへの依存。S1_fmap_variants.py grid astrid2 81 61 の表 (S1_fmap_grid_astrid2_81x61.table.txt、4941 点、1024 周) から、
1 つおき (step 2 = 41 × 31、Q1 と同じ格子)、3 つおき (step 4 = 21 × 16)、7 つおき (step 8 = 11 × 8) の部分格子を取り出し、D の分布の統計を比べる。
(部分格子の点は 81 × 61 の点そのもの = 同じ追跡・同じ D。差は標本の密度だけ。)  使い方: py -3.12 example/v2_3/S1_fmap_grid_summary.py"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from pysad.freqmap import FMapResult, line_distance

r = FMapResult.load_table(os.path.join(HERE, 'S1_fmap_grid_astrid2_81x61.table.txt'))
nx, ny = r.grid[0], r.grid[1]
idx_all = np.arange(nx * ny).reshape(nx, ny)
print('ASTRID2 位置 2、x = 0.5〜18 mm、y = 0.5〜14 mm、1024 周、81 × 61 の格子から部分格子を取り出す')
print('%-12s %6s %8s %8s %10s %10s %12s %12s %s' % ('格子', '点数', 'D 中央値', 'D 95%', 'D>-5 割合', 'D>-4 割合', 'D>-5 の点数×(4941/点数)', 'D>-4 の点数×(4941/点数)', '近い/遠い D>-5 の割合'))
for step in (8, 4, 2, 1):
    idx = idx_all[::step, ::step].ravel()
    D = r.D[idx]
    ok = np.isfinite(D)
    d = line_distance(r.nux1[idx][ok], r.nuy1[idx][ok], 4)
    near = d < 0.003
    sc = 4941.0 / len(idx)
    print('%-12s %6d %8.2f %8.2f %10.3f %10.3f %12.0f %12.0f   %.3f / %.3f' % ('%d x %d' % (len(range(0, nx, step)), len(range(0, ny, step))), len(idx), np.median(D[ok]), np.percentile(D[ok], 95),
          np.mean(D[ok] > -5), np.mean(D[ok] > -4), (D[ok] > -5).sum() * sc, (D[ok] > -4).sum() * sc,
          np.mean(D[ok][near] > -5) if near.any() else np.nan, np.mean(D[ok][~near] > -5)))
# 1 つおきの格子の D の、隣の点からの補間誤差: 81 × 61 の奇数番の点の D を、両隣 (x 方向) の偶数番の点の平均で予測したときの差
D2 = r.D.reshape(nx, ny)
pred = 0.5 * (D2[0:-2:2, :] + D2[2::2, :])
act = D2[1:-1:2, :]
ok = np.isfinite(pred) & np.isfinite(act)
print('x 方向の隣の点 (間隔 0.43 mm) の D の平均で、間の点の D を予想した差 |ΔD|: 中央値 %.2f、90%% %.2f、|ΔD| > 1 の割合 %.3f' % (
    np.median(np.abs(pred - act)[ok]), np.percentile(np.abs(pred - act)[ok], 90), np.mean(np.abs(pred - act)[ok] > 1)))
