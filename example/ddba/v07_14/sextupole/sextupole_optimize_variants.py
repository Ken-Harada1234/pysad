# -*- coding: utf-8 -*-
"""[試作] 最適化の変種 (R7): 基準の作用 J_ref を変えた全項の最適化、δ 依存 (∂ADTS/∂δ) の項を足した最適化。results/candidates.json に追記する。
使い方: py -3.12 sextupole_optimize_variants.py [--jref 3e-7,3e-6,1e-5]
**注意**: これらは「物理的な重みの選び方を変えて比べる」実験で、どれが追跡の DA に良いかは追跡で見ている。最良のものを選ぶと、追跡を使って重みを選んだことになる (README で区別して書く)。
J_ref の物理的な意味: 目的関数の各項を、その作用での 1 周あたりのハミルトニアンの大きさ [rad] にそろえている (docs: sxtopt.py の説明)。
  1e-6 m rad は x = 2.3 mm (βx = 5.4 m) 程度、3e-6 は 4 mm、1e-5 は 7 mm。
(∂ADTS/∂δ の項を足す変種は試したが、δ のまわりの ADTS が K2 の二次式にならない (sxt_delta_forms.py --check: 二次形式のフィットが直接の値と 1e2 倍違う) ので入れていない。)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()
from pysad import resdrive as R
from pysad import sxtopt as S

ap = argparse.ArgumentParser()
ap.add_argument('--jref', default='3e-7,3e-6,1e-5')
ap.add_argument('--nstart', type=int, default=12)
args = ap.parse_args()
path = os.path.join(C.OUT, 'candidates.json')
CJ = json.load(open(path, encoding='utf-8'))
src, lat = C.load_lattice('RING')
m = S.SxtModel(lat, C.P0)
assert m.names == CJ['family_names']
kf = np.array(CJ['candidates']['final']['k'])
xi0, M = m.xi_matrix()
xi_an = np.array(CJ['xi_an_final'])
xi_calc_final = np.array(CJ['xi_calc_final'])
kmax = CJ['kmax']
best, pair = 0, (0, 1)
for i in range(m.F):
    for j in range(i + 1, m.F):
        d = abs(np.linalg.det(M[:, [i, j]] * np.abs(kf[[i, j]])))
        if d > best:
            best, pair = d, (i, j)


def trim_xi(k, nit=4):
    k = k.copy()
    Mp = M[:, list(pair)]
    for _ in range(nit):
        xi = np.array(R.chromaticity_numeric(S.with_k2l(lat, m.names, k), order=3, dmax=2e-3)['d1'])
        r = xi_calc_final - xi
        if np.max(np.abs(r)) < 1e-4:
            break
        k[list(pair)] += np.linalg.solve(Mp, r)
    return k


def run(obj, tag, desc):
    sols = []
    for sd in range(args.nstart):
        rng = np.random.default_rng(500 + sd)
        k0 = kf.copy() if sd == 0 else S.project_to_xi(m, np.clip(kf * (1 + 0.3 * rng.standard_normal(m.F)), -kmax, kmax), xi_an)
        try:
            k, r = S.optimize(obj, k0, xi_an, kmax, maxiter=600)
        except Exception:
            continue
        sols.append((obj.value(k), k))
    sols.sort(key=lambda t: t[0])
    v, k = sols[0]
    k = trim_xi(k)
    a = m.adts(k)
    print('%s: 最適化の目的 %.3e (最終の組 %.3e)、ADTS = %s' % (tag, v, obj.value(kf), np.round([float(a[c]) for c in ('xx', 'xy', 'yy')], 0).tolist()))
    CJ['candidates'][tag] = dict(k=[float(x) for x in k], desc=desc, group='optimized', trim=None)


for jr in [float(x) for x in args.jref.split(',') if x]:
    run(S.SxtObjective(m, J_ref=jr), 'opt_J%g' % jr, '全項の最適化 (J_ref = %g m rad、kmax = 最終の最大)。重みの選び方の変種' % jr)
json.dump(CJ, open(path, 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
