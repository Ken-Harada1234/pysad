# -*- coding: utf-8 -*-
"""results/gs_*.json (gs_check.py) を候補 × 誤差の大きさごとにまとめる: gs の DA 面積 (格子 9 光線 x 96 点)・x・y の広さの平均と最悪、シードごと。
pysad と同じ点を追跡した gs_*.json は生き残りの一致も出す。
  py -3.12 summarize_gs.py [--md results/summary_gs.md]"""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

ap = argparse.ArgumentParser()
ap.add_argument('--md', default=None)
args = ap.parse_args()
rows = defaultdict(list)
for f in sorted(glob.glob(os.path.join(V.OUT, 'gs_*.json'))):
    d = json.load(open(f, encoding='utf-8'))
    rows[(d['rms'], d['name'])].append(d)
lines = ['| 候補 | 誤差 rms | シード数 | gs の DA 面積 mm² 平均 | 最悪 | 各シード | x の広さ mm 平均 (最悪) | y の広さ mm 平均 (最悪) | pysad と同じ点の一致 |', '|---|---|---|---|---|---|---|---|---|']
order = {'final': 0, 'r1a': 1, 'r6a': 2, 'r6b': 3}
for (rms, name), ds in sorted(rows.items(), key=lambda kv: (kv[0][0], order.get(kv[0][1], 9), kv[0][1])):
    ds = sorted(ds, key=lambda d: d['seed'])
    a = np.array([d['gs']['area_mm2'] for d in ds])
    x = np.array([d['gs']['x_mm'] for d in ds])
    y = np.array([d['gs']['y_mm'] for d in ds])
    ag = [d for d in ds if 'agree' in d]
    agree = ', '.join('%.4f (%d/%d 点が違う)' % (d['agree'], d['n_disagree'], d['n']) for d in ag) or '-'
    lines.append('| %s | %g%% | %d | %.1f | %.1f | %s | %.2f (%.2f) | %.2f (%.2f) | %s |' % (name, 100 * rms, len(ds), a.mean(), a.min(), ' '.join('%.1f' % v for v in a), x.mean(), x.min(), y.mean(), y.min(), agree))
txt = '\n'.join(lines)
sys.stdout.buffer.write(txt.encode('utf-8') + bytes([10]))
if args.md:
    open(args.md, 'w', encoding='utf-8').write(txt + '\n')
