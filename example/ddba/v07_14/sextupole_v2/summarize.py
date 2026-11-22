# -*- coding: utf-8 -*-
"""results/eval_*.json を候補 × 誤差の大きさごとにまとめる: シードごとの DA (面積・x・y)、平均・最悪 (最小)、x-δ の有効な δ の範囲、周波数マップの「悪い点」の割合。
  py -3.12 summarize.py [--md results/summary.md]
"""
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
ap.add_argument('--names', default=None)
ap.add_argument('--suffix', default='', help='_c など。空なら suffix なしのファイルだけ')
args = ap.parse_args()
rows = defaultdict(list)
for f in sorted(glob.glob(os.path.join(V.OUT, 'eval_*.json'))):
    base = os.path.basename(f)[:-5]
    if (args.suffix and not base.endswith(args.suffix)) or (not args.suffix and base.endswith('_c')):
        continue
    d = json.load(open(f, encoding='utf-8'))
    rows[(d['name'], d['rms'])].append(d)
want = args.names.split(',') if args.names else None
lines = ['| 候補 | 誤差 rms | シード数 | DA 面積 mm² 平均 | 最悪 | 各シード | x の広さ mm 平均 (最悪) | y の広さ mm 平均 (最悪) | x-δ の広さ mm (δ 平均) 平均 (最悪) | 有効 δ の幅 % 平均 (最悪) | 周波数マップ悪い点の割合 平均 (最悪) |',
         '|---|---|---|---|---|---|---|---|---|---|---|']


def key(kv):
    (n, r), _ = kv
    return (n.startswith('r') and n[1:2].isdigit(), n, r)


for (name, rms), ds in sorted(rows.items(), key=lambda kv: (kv[0][1], kv[0][0])):
    if want and name not in want:
        continue
    ds = sorted(ds, key=lambda d: d['seed'])
    a = np.array([d['da']['area_mm2'] for d in ds if 'da' in d])
    x = np.array([d['da']['x_da_mm'] for d in ds if 'da' in d])
    y = np.array([d['da']['y_da_mm'] for d in ds if 'da' in d])
    xw = np.array([d['xdp']['mean_width_mm'] for d in ds if 'xdp' in d])
    dp = np.array([(d['xdp']['valid_dp_range'][1] - d['xdp']['valid_dp_range'][0]) * 100 for d in ds if 'xdp' in d])
    fb = np.array([d['fmap']['frac_bad'] for d in ds if 'fmap' in d])

    def ms(v, f='%.1f', worst=min):
        return '-' if len(v) == 0 else (f + ' (' + f + ')') % (np.mean(v), worst(v))
    lines.append('| %s | %g%% | %d | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        name, 100 * rms, len(ds), '%.1f' % a.mean() if len(a) else '-', '%.1f' % a.min() if len(a) else '-', ' '.join('%.1f' % v for v in a),
        ms(x, '%.2f'), ms(y, '%.2f'), ms(xw, '%.1f'), ms(dp, '%.1f'), ms(fb, '%.2f', max)))
txt = '\n'.join(lines)
sys.stdout.buffer.write(txt.encode('utf-8') + bytes([10]))
if args.md:
    open(args.md, 'w', encoding='utf-8').write(txt + '\n')
