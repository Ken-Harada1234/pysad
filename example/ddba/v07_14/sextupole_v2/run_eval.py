# -*- coding: utf-8 -*-
"""候補の 6 極 (results/cands.json の K2L) を、K2L のランダム誤差 (rms 相対、複数のシード) つきで追跡して評価する。1 候補・1 誤差・1 シードごとに 1 ファイル
results/eval_<候補>_e<rms%>_s<シード>.json (既にあれば飛ばす = 上書きしない)。シード 0 は誤差なし。
  py -3.12 run_eval.py NAME [NAME ...] --rms 0.03 --seeds 1,2,3 [--parts da,xdp,fmap] [--nturn 512]
簡略化した設定 (sextupole/sextupole_tracking_check.py と同じ判定): RF・放射なし、δ = 0、リング先頭から px = py = 0 で出発、LOSSAMPL = ±50 mm。
  da: x-y の 9 本の光線 (規格化した振幅 u ≤ 4e-3 [√m] まで、二分 2 回)、512 周。DA 面積・x 方向・y 方向の広さ。xdp: δ = −3〜+3% で x の広さ (有効な δ の範囲 = 広さ > 0 の δ の幅)。
  fmap: x ≤ 4 mm、y ≤ 3 mm の 16 × 10 点、1024 周 (DA の内側を細かく見る。外側の失われる点はこの設定では少ない)。
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

ap = argparse.ArgumentParser()
ap.add_argument('names', nargs='+')
ap.add_argument('--rms', type=float, default=0.03)
ap.add_argument('--seeds', default='1,2,3')
ap.add_argument('--parts', default='da,xdp,fmap')
ap.add_argument('--nturn', type=int, default=512)
ap.add_argument('--cands', default=os.path.join(V.OUT, 'cands.json'))
ap.add_argument('--quick', action='store_true')
ap.add_argument('--suffix', default='', help='出力のファイル名の末尾 (最終の比較の一式は _c。既存の結果と混ぜない)')
args = ap.parse_args()
V.lower_priority()
import sxt_track as T

src, lat, m, ks, xt = V.load()
cands = json.load(open(args.cands, encoding='utf-8'))
for name in args.names:
    k = np.array(cands[name]['k'])
    for seed in [int(s) for s in args.seeds.split(',')]:
        rms = 0.0 if seed == 0 else args.rms
        tag = 'eval_%s_e%g_s%d%s' % (name, 100 * rms, seed, args.suffix)
        path = os.path.join(V.OUT, tag + '.json')
        if os.path.exists(path):
            print('skip', tag)
            continue
        t0 = time.time()
        eps = None if seed == 0 else V.sext_error_vector(lat, rms, seed)
        st = V.make_setup(src, m, k, eps)
        out = dict(name=name, rms=rms, seed=seed, nturn=args.nturn, k=k.tolist())
        for p in args.parts.split(','):
            t1 = time.time()
            if p == 'da':
                out['da'] = T.da_rays(st, 1, args.nturn, args.quick)
            elif p == 'xdp':
                out['xdp'] = T.da_xdp(st, 1, args.nturn, args.quick)
                w = np.array(out['xdp']['width_mm'])
                dps = np.array(out['xdp']['dp'])
                out['xdp']['valid_dp_range'] = [float(dps[w > 0].min()), float(dps[w > 0].max())] if (w > 0).any() else [0.0, 0.0]
            elif p == 'fmap':
                r = V.fmap_run(st, nx=16, ny=10, xmax=4e-3, ymax=3e-3)
                out['fmap'] = V.fmap_summary(r)
                V.save_fmap(os.path.join(V.OUT, tag + '_fmap.npz'), r)
                out['fmap']['lines'] = V.line_scores(r)[:12]
            out.setdefault('time', {})[p] = time.time() - t1
        json.dump(out, open(path, 'w', encoding='utf-8'), ensure_ascii=False)
        d = out.get('da')
        print('%s: %.0f s DA面積 %s x %s y %s' % (tag, time.time() - t0, '%.1f' % d['area_mm2'] if d else '-', '%.2f' % d['x_da_mm'] if d else '-', '%.2f' % d['y_da_mm'] if d else '-'), flush=True)
