# -*- coding: utf-8 -*-
"""v2.4.0 A6: 6 極の検討 (R7) の解析量に、ampnf の J² の項・δ 依存 (交差項) を足すと、追跡の判定との順位相関が良くなるかを調べる。
候補 (results/candidates.json の 32 組) ごとに、ampnf (6 極・8 極だけの模型、D = 5、厚い 6 極 4 分割。偏向磁石・端部・運動学的な項は入れない: 6 極の組に依る部分だけ。
偏向磁石との交差項は入っていない (J¹ で 3% の程度)) の
  J² の係数 (c20, c11, c02 の νx・νy)、∂ADTS/∂δ (c_{a,b,1} × 1!)、∂²ADTS/∂δ² (2! c_{a,b,2})、ξ2 (2 c_002)
を出し、基準の振幅 J_ref (Jx = Jy) での J² によるチューンのずれ |Δν_J2(J_ref)| などを作る。追跡の判定 (results/track_<名前>.json、R7 の結果) との Spearman の順位相関を、
R7 の解析量 (目的関数・ADTS の項・|∂ADTS/∂δ| (nform)) と並べる。
使い方: PYTHONUTF8=1 py -3.12 sextupole_j2_check.py   (results/j2_table.json、results/j2_correlation.out.txt に保存)
"""
import glob
import json
import os
import re
import sys
import time

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sxt_common as C
from pysad import ampnf as A
from pysad import sxtopt as S

OUTF = os.path.join(C.OUT, 'j2_table.json')
JREFS = (3e-7, 1e-6)


def compute():
    src, lat = C.load_lattice('RING')
    cand = json.load(open(os.path.join(C.OUT, 'candidates.json'), encoding='utf-8'))
    names = cand['family_names']
    rows = {}
    t0 = time.time()
    for nm, c in cand['candidates'].items():
        lat2 = S.with_k2l(lat, names, np.array(c['k']))
        r = A.amplitude_tune(lat2, D=5, nslice=4)
        d = dict(J1=[r.adts['xx'], r.adts['xy'], r.adts['yy']], J2=[r.adts2[k] for k in ('x20', 'x11', 'x02', 'y20', 'y11', 'y02')],
                 d1=[r.adts_delta(0, 1, 0, 1), r.adts_delta(0, 0, 1, 1), r.adts_delta(1, 0, 1, 1)],
                 d2=[r.adts_delta(0, 1, 0, 2), r.adts_delta(0, 0, 1, 2), r.adts_delta(1, 0, 1, 2)],
                 xi2=[r.chrom(0, 2), r.chrom(1, 2)], min_den=r.min_den)
        for j in JREFS:
            d['dnuJ2_%g' % j] = [float(abs(r.coef(0, 2, 0) * j ** 2 + r.coef(0, 1, 1) * j ** 2 + r.coef(0, 0, 2) * j ** 2)),
                                 float(abs(r.coef(1, 2, 0) * j ** 2 + r.coef(1, 1, 1) * j ** 2 + r.coef(1, 0, 2) * j ** 2))]
        rows[nm] = d
        print('%-14s %.0f s  J1 %s  |J2|max %.2e' % (nm, time.time() - t0, ['%.0f' % v for v in d['J1']], max(abs(v) for v in d['J2'])), flush=True)
    json.dump(rows, open(OUTF, 'w'), indent=1)
    return rows


def correlate(rows):
    rdt = json.load(open(os.path.join(C.OUT, 'rdt_table.json'), encoding='utf-8'))
    tr = {}
    for f in glob.glob(os.path.join(C.OUT, 'track_*.json')):
        if f.endswith('_quick.json') or re.search(r'_(err\d+|n\d+)(_|\.json)', os.path.basename(f)):
            continue
        d = json.load(open(f, encoding='utf-8'))
        tr[d['name']] = d
    TRACK = {'DA面積': (lambda d: d['da']['area_mm2'], 1), 'x_DA': (lambda d: d['da']['x_da_mm'], 1), 'y_DA': (lambda d: d['da']['y_da_mm'], 1),
             'xδ幅': (lambda d: d['xdp']['mean_width_mm'], 1), 'FMカオス割合': (lambda d: d['fmap']['frac_chaotic'], -1), 'LMA': (lambda d: d['lma']['mean_min'], 1)}
    n2 = lambda v: float(np.sqrt(np.sum(np.square(v))))
    ANA = {
        '目的関数(全項) [R7]': lambda n: rdt[n]['obj'],
        '|ADTS| [R7、J¹]': lambda n: rdt[n]['adts_norm'],
        '|∂ADTS/∂δ| [R7、nform]': lambda n: rdt[n]['dadts_norm'],
        '|ξ2| (CALC) [R7]': lambda n: float(np.hypot(*rdt[n]['xi2_calc'])),
        '|J²| (ampnf、6 極のみ、‖c‖)': lambda n: n2(rows[n]['J2']),
        '|Δν_J2| (J_ref 3e-7 と 1e-6 は同じ順位)': lambda n: n2(rows[n]['dnuJ2_3e-07']),
        '|∂ADTS/∂δ| (ampnf、6 極のみ)': lambda n: n2(rows[n]['d1']),
        '|∂²ADTS/∂δ²| (ampnf、6 極のみ)': lambda n: n2(rows[n]['d2']),
    }
    GROUPS = {'all': lambda n: True, 'near': lambda n: rdt[n]['obj_rel'] < 30 and n not in ('opt_rdt1', 'opt_adts', 'opt_rdt4'),
              'optimized': lambda n: rdt[n]['group'] == 'optimized'}
    names_all = [n for n in tr if n in rdt and n in rows]
    out = []

    def P(s=''):
        print(s)
        out.append(s)
    P('追跡できた組 %d 個' % len(names_all))
    for g, sel in GROUPS.items():
        ns = [n for n in names_all if sel(n)]
        P('\n[群 %s: %d 組] Spearman ρ (解析量は小さいほど良い、追跡は「良いほど大きい」に符号をそろえる。ρ < 0 が「解析量が小さい組ほど追跡が良い」)' % (g, len(ns)))
        P('| 解析量 | ' + ' | '.join(TRACK) + ' |')
        P('|' + '---|' * (len(TRACK) + 1))
        for an, fa in ANA.items():
            cells = []
            for tn, (ft, sg) in TRACK.items():
                xs, ys = [], []
                for n in ns:
                    try:
                        y = sg * ft(tr[n])
                        x = fa(n)
                    except (KeyError, TypeError):
                        continue
                    if y is None or not np.isfinite(y) or not np.isfinite(x):
                        continue
                    xs.append(x)
                    ys.append(y)
                if len(xs) < 5 or len(set(ys)) < 2:
                    cells.append('-')
                    continue
                rho, p = spearmanr(xs, ys)
                cells.append('%+.2f (%.0e)' % (rho, p))
            P('| %s | %s |' % (an, ' | '.join(cells)))
    open(os.path.join(C.OUT, 'j2_correlation.out.txt'), 'w', encoding='utf-8').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    rows = json.load(open(OUTF)) if os.path.exists(OUTF) and '--again' not in sys.argv else compute()
    correlate(rows)
