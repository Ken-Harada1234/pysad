# -*- coding: utf-8 -*-
"""∂(ADTS)/∂δ の多項式モデル (R7): 追跡の写像の数値の正規形 (nform.adts_at_delta、δ = ±3e-3 の中心差分) を、6 極の K2L の二次式として偏極化でフィットする。
使い方: py -3.12 sxt_delta_forms.py [--check]    出力: results/dadts_forms.npz (c (3,)、b (F,3)、S (F,F,3))、results/dadts_forms.out.txt
δ のまわりの ADTS は、閉軌道のずれによる6 極のフィードダウン (K1 = K2 x0) を通して K2 の 3 次以上にも依存するので、厳密な二次式ではない。
--check: 乱数の K2L (最終の周り 10%・30%) と最適化の組で、モデルと直接の値の相対差 (|差|/|値|、3 成分のノルム) を出す。
1 点 2.6 s × (1 + 2×20 + 190) 点 = 10 分。
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
C.lower_priority()
from pysad import nform
from pysad import sxtopt as S

ap = argparse.ArgumentParser()
ap.add_argument('--check', action='store_true')
ap.add_argument('--delta0', type=float, default=3e-3)
args = ap.parse_args()
src, lat = C.load_lattice('RING')
CAND = json.load(open(os.path.join(C.OUT, 'candidates.json'), encoding='utf-8'))
names = CAND['family_names']
kf = np.array(CAND['candidates']['final']['k'])


def dadts(k):
    lat2 = S.with_k2l(lat, names, k)
    d = nform.adts_delta_derivative(lat2, momentum_eV=C.P0, delta0=args.delta0, method='nform', amp3=1e-5)
    return np.array([d['xx'], d['xy'], d['yy']], dtype=float)


path = os.path.join(C.OUT, 'dadts_forms.npz')
if args.check and os.path.exists(path):
    z = np.load(path)
    cbS = (z['c'], z['b'], z['S'])
else:
    t0 = time.time()
    cbS = S.quad_fit(dadts, len(names))
    np.savez(path, c=cbS[0], b=cbS[1], S=cbS[2])
    print('フィット %.0f s' % (time.time() - t0))
if args.check:
    out = []
    for nm in ('final', 'opt_full_1', 'opt_full_k125', 'retune1', 'rev9', 'v07_12', 'rnd10_0', 'rnd10_3', 'rnd30_5', 'opt_adts'):
        k = np.array(CAND['candidates'][nm]['k'])
        a, b = dadts(k), S.quad_eval(cbS, k)
        out.append('%-14s 直接 %s  モデル %s  相対差 %.3f' % (nm, np.round(a, -3).tolist(), np.round(b, -3).tolist(), np.linalg.norm(a - b) / np.linalg.norm(a)))
    print('\n'.join(out))
    open(os.path.join(C.OUT, 'dadts_forms.out.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')
