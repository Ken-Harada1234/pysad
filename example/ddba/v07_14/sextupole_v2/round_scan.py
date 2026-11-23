# -*- coding: utf-8 -*-
"""感度の走査 (人手の試行錯誤の道具): 基準の候補 BASE の主な項 (3 次の RDT 5 項、δ の 3 項、ADTS 3 成分、h20200、h20020) を 1 つずつ 0.3 倍へ動かし (目的)、他の項は BASE の値の ±tol の制約、
変化の上限 dk、ξ は BASE の値の等式、で sxtopt.design_k2l を呼んで候補 TAG_<項> を results/cands.json に足す (同名があれば止まる)。どの項の変更が追跡の DA に効くかを見るために使う (走査の後の追跡は run_eval.py)。
  py -3.12 round_scan.py BASE TAG [--tol 0.10] [--dk 4] [--factor 0.3] [--terms h21000,h00201,...]
反復も探索も自動化していない: 項の選択・重み・目標値は人が決める。この走査は 1 回の入力 (13 通りの指定) を作って関数を 13 回呼ぶだけ。"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxv2_common as V

S = V.S
ap = argparse.ArgumentParser()
ap.add_argument('base')
ap.add_argument('tag')
ap.add_argument('--tol', type=float, default=0.10)
ap.add_argument('--dk', type=float, default=4.0)
ap.add_argument('--factor', type=float, default=0.3)
ap.add_argument('--terms', default=None)
args = ap.parse_args()
src, lat, m, ks, xt = V.load()
cands = json.load(open(os.path.join(V.OUT, 'cands.json'), encoding='utf-8'))
k0 = np.array(cands[args.base]['k'])
ev = {r['name']: r for r in S.SxtSpecObjective(m, 'w3=objective, chrom=objective, adts=objective, h20200=objective, h20020=objective', warn=False).table(k0, xt)}
base = {n: r['value'] for n, r in ev.items()}
names = args.terms.split(',') if args.terms else list(base)
for t in names:
    spec = {}
    for n, v in base.items():
        if n == t:
            spec[n] = ('objective', 1.0, args.factor * v)
        else:
            spec[n] = ('constraint', 1.0, v, args.tol * abs(v) + (1e-3 if n.startswith('h') and abs(v) < 1 else 0.0))
    r = S.design_k2l(m, spec, k0=k0, xi_target=xt, kmax=V.KMAX, dk=args.dk, warn=False)
    nm = '%s_%s' % (args.tag, t)
    V.add_cand(nm, r['k'], 'ラウンドの走査: %s から出発して %s を %.1f 倍へ (目的)、他の主な項は ±%d%% の制約、変化の上限 dk=%g、ξ は %s の値' % (args.base, t, args.factor, round(100 * args.tol), args.dk, args.base),
               'scan %s x%.2f tol=%.2f dk=%g' % (t, args.factor, args.tol, args.dk))
    print('%s ok=%s viol=%.2g maxdk=%.2f | %s' % (nm, r['success'], r['constraint_violation'], np.max(np.abs(r['k'] - k0)), V.show_k(m, '', r['k'])), flush=True)
