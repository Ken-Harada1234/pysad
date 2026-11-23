# -*- coding: utf-8 -*-
"""格子の細かさへの依存の再確認。fmap_<ring>.table.txt の 41 × 31 の表 (ASTRID2・KEK-LS) から 1 つおき・2 つおき・3 つおきの部分格子を取り出し、D の分布 (百分位・D > −5・D > −4・失われた点) を比べる
(同じ追跡の部分集合なので、格子の細かさで D の値は変わらない。分布の統計がどこまで粗い格子で再現するかを見る)。
使い方: py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_check_grid.py     出力: 標準出力"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
from pysad.freqmap import FMapResult

for ring in ('astrid2', 'kekls'):
    r = FMapResult.load_table(os.path.join(HERE, 'fmap_%s.table.txt' % ring))
    nx, ny = r.grid[0], r.grid[1]
    D = r.D.reshape(nx, ny)
    print('%s (Q1 の表 %d × %d、1024 周)' % (ring, nx, ny))
    for st in (1, 2, 3, 4, 5):
        d = D[::st, ::st].ravel()
        nl = int(np.isnan(d).sum())
        d = d[np.isfinite(d)]
        print('  %d つおき: %3d × %-3d = %4d 点、失われた %3d、D の百分位 5/25/50/75/95 = %s、D > -5: %.3f、D > -4: %.3f' % (
            st, D[::st, ::st].shape[0], D[::st, ::st].shape[1], D[::st, ::st].size, nl,
            ' '.join('%6.2f' % v for v in np.percentile(d, [5, 25, 50, 75, 95])), np.mean(d > -5), np.mean(d > -4)))
