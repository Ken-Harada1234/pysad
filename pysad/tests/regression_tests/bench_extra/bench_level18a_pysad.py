#!/usr/bin/env python3
"""
bench_level18a_pysad.py — Level 18: CALC6D (pysad.calc6d = qcell6d/temit/tfetwiss の直訳) の照合
対応する SAD スクリプト: bench_level18a.sad (RFSW RAD), 18b (NORFSW RAD), 18c (RFSW NORAD)
  USE RING; CELL; CALC6D; CALC; Twiss[28 成分, 位置] を $$$, 1, QFA1H.3, SFA.5, MKC.2, RFC, 100, 200.5 で比較
使い方: python bench_level18a_pysad.py [18a|18b|18c]
"""
import os
import re
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.calc6d import calc6d, twiss_frac

tag = sys.argv[1] if len(sys.argv) > 1 else '18a'
flags = {'18a': dict(rfsw=True), '18b': dict(rfsw=False), '18c': dict(rfsw=True)}[tag]
KEYS = ["AX", "BX", "NX", "AY", "BY", "NY", "EX", "EPX", "EY", "EPY", "R1", "R2", "R3", "R4", "DETR",
        "DX", "DPX", "DY", "DPY", "DZ", "DDP", "AZ", "BZ", "NZ", "ZX", "ZPX", "ZY", "ZPY"]
FIELDS = ('alfx', 'betx', 'mux', 'alfy', 'bety', 'muy', 'etx', 'etpx', 'ety', 'etpy',
          'r1', 'r2', 'r3', 'r4', 'detr', 'x_cod', 'px_cod', 'y_cod', 'py_cod', 'z_cod',
          'delta_cod', 'alfz', 'betz', 'muz', 'zx', 'zpx', 'zy', 'zpy')

s = open(os.path.join(HERE, f'bench_level{tag}.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadvec(label):
    m = re.search(r'^' + label + r' \{(.*?)\}', s, re.M)
    return np.array([float(v) for v in m.group(1).split(',')])


def sadlist(label):
    m = re.search(r'^' + label + r' \{(.*?)\}', s, re.M)
    return [t.strip().strip('"') for t in m.group(1).split(',')]


sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
els = lat.elements
nlat = int(float(re.search(r'^NLAT\s+([\d.]+)', s, re.M).group(1)))
print('NLAT SAD', nlat, 'pysad', len(els) + 1)
names = sadlist('NAMES')
print('NAMES SAD', names, 'pysad', [els[99].name, els[199].name, els[200].name])


def pos(name, n):
    """name の n 番目 (1 始まり) の SAD 位置 l"""
    k = 0
    for i, e in enumerate(els):
        if e.name == name:
            k += 1
            if k == n:
                return i + 1
    raise KeyError(name)


t0 = time.time()
res = calc6d(lat, 3e9, out=sys.stdout, **flags)
print(f'calc6d {time.time() - t0:.1f} s  stable={res.stable}  tunes={res.tunes}  trf0={res.trf0}')
rows = res.rows
cases = [('TW_END', rows[-1]), ('TW_ORG', rows[0]),
         ('TW_QFA1H3', rows[pos('QFA1H', 3) - 1]), ('TW_SFA5', rows[pos('SFA', 5) - 1]),
         ('TW_MKC2', rows[pos('MKC', 2) - 1]), ('TW_RFC', rows[pos('RFC', 1) - 1]),
         ('TW_P100', rows[99]), ('TW_P200H', twiss_frac(res, 200, 0.5))]
worst = 0.0
for label, row in cases:
    ref = sadvec(label)
    val = np.array([getattr(row, f) for f in FIELDS])
    print(f'--- {label} ({row.name})')
    for k, a, b in zip(KEYS, val, ref):
        d = a - b
        rel = abs(d) / max(abs(b), 1e-300)
        flag = ''
        if abs(d) > 1e-9 and rel > 1e-9:
            flag = '  <<<'
            worst = max(worst, min(abs(d), rel))
        print(f'  {k:5s} pysad {a: .14e}  SAD {b: .14e}  diff {d: .2e}{flag}')
print('WORST (min of abs/rel diff over flagged):', worst)
if res.emit is not None:
    e = res.emit
    print('EMIT pysad', e['emittance_x'], e['emittance_y'], e['emittance_z'], e['momentum_spread'], e['bunch_length'])
