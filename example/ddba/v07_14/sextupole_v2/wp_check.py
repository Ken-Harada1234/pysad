# -*- coding: utf-8 -*-
"""作動点の整数・半整数共鳴の確認: 周波数マップ (_c の一式) の生き残った点の端数のチューンの範囲と、整数 (νx, νy = 整数)・半整数 (2νx, 2νy = 整数) の線までの最短距離 (チューン面の距離) を、候補・誤差の大きさごとに出す。
線形のチューン (47.865, 16.655) は誤差・6 極の K2L に依らない (6 極は 1 次ではチューンを変えない)。
  py -3.12 wp_check.py [--suffix _c]"""
import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

ap = argparse.ArgumentParser()
ap.add_argument('--suffix', default='_c')
args = ap.parse_args()
res = defaultdict(list)
for f in sorted(glob.glob(os.path.join(V.OUT, 'eval_*%s_fmap.npz' % args.suffix))):
    tag = os.path.basename(f)[len('eval_'):-len('_fmap.npz')]
    parts = tag.split('_')
    name = '_'.join(parts[:-3]) if len(parts) > 3 else parts[0]
    # tag: <name>_e<rms>_s<seed>_c
    name = tag[:tag.index('_e')]
    e = tag[tag.index('_e') + 2:].split('_')[0]
    r = V.Saved(f)
    fx, fy = r.nux1, r.nuy1
    ok = np.isfinite(fx) & np.isfinite(fy)
    fx, fy = fx[ok], fy[ok]
    d_int_x = np.minimum(fx, 1 - fx)
    d_int_y = np.minimum(fy, 1 - fy)
    d_half_x = np.abs(fx - 0.5)
    d_half_y = np.abs(fy - 0.5)
    res[(name, e)].append((fx.min(), fx.max(), fy.min(), fy.max(), d_int_x.min(), d_half_x.min(), d_int_y.min(), d_half_y.min()))
lines = ['| 候補 | 誤差 rms % | 点の数 (シード) | νx 端数の範囲 | νy 端数の範囲 | 整数 νx までの最短距離 | 半整数 2νx までの最短距離 (νx = 0.5 からの距離) | 整数 νy まで | 半整数 2νy まで (νy = 0.5 からの距離) |', '|---|---|---|---|---|---|---|---|---|']
for (name, e), v in sorted(res.items(), key=lambda kv: (float(kv[0][1]), kv[0][0])):
    a = np.array(v)
    lines.append('| %s | %s | %d | %.3f〜%.3f | %.3f〜%.3f | %.3f | %.3f | %.3f | %.3f |' % (name, e, len(v), a[:, 0].min(), a[:, 1].max(), a[:, 2].min(), a[:, 3].max(), a[:, 4].min(), a[:, 5].min(), a[:, 6].min(), a[:, 7].min()))
sys.stdout.buffer.write(('\n'.join(lines)).encode('utf-8') + bytes([10]))
