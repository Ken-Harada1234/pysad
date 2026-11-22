# -*- coding: utf-8 -*-
"""候補の局所運動量アクセプタンス (LMA) の簡略化した確認 (誤差なし): sextupole/sxt_track.lma と同じ設定 (RF 6.5 MV・放射あり、LOSSAMPL = ±20 mm、シンクロトロン振動 2 周期、
リング全周に約 190 m おき 4 点)。x-y の DA の改善が、運動量方向 (LMA) を犠牲にしていないかの確認用。
  py -3.12 run_lma.py NAME [NAME ...]    → results/lma_<名前>.json (あれば飛ばす)"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

V.lower_priority()
import sxt_track as T

src, lat, m, ks, xt = V.load()
cands = json.load(open(os.path.join(V.OUT, 'cands.json'), encoding='utf-8'))
for name in sys.argv[1:]:
    path = os.path.join(V.OUT, 'lma_%s.json' % name)
    if os.path.exists(path):
        print('skip', name)
        continue
    k = np.array(cands[name]['k'])
    t0 = time.time()
    r = T.lma(V.S.replace_k2_in_source(src, m.names, k), False)
    r['name'] = name
    json.dump(r, open(path, 'w', encoding='utf-8'), ensure_ascii=False)
    print('%s: %.0f s | LMA の 4 点の min(δ+, |δ-|) の平均 %.2f%%、δ+ 平均 %.2f%%、|δ-| 平均 %.2f%% (RF バケット高さ %s)' % (name, time.time() - t0, 100 * r['mean_min'], 100 * r['mean_plus'], 100 * r['mean_minus'], r.get('bucket')), flush=True)
