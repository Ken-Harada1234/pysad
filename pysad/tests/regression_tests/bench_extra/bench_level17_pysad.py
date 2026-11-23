#!/usr/bin/env python3
"""
bench_level17_pysad.py — Level 17: INTRA (IBS + Touschek 寿命表) の照合 (pysad.touschek = temit/tintrb の直訳)
対応する SAD スクリプト: bench_level17a.sad (RFSW RAD FLUC INTRA), bench_level17b.sad (+RADCOD)
  参照値: bench_level17a.result / bench_level17b.result
使い方: python bench_level17_pysad.py [a|b]
"""
import ast
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.emit6d import emittance6d
from pysad.touschek import touschek_table, touschek_lifetime

which = sys.argv[1] if len(sys.argv) > 1 else 'a'
radcod = (which == 'b')
res = os.path.join(HERE, f'bench_level17{which}.result')
s = open(res, encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadval(tag):
    m = re.search(r'^' + tag + r'\s+(.*)$', s, re.M)
    txt = m.group(1).strip()
    txt = txt.replace('{', '[').replace('}', ']')
    txt = re.sub(r'(?<![\d.])\.(\d)', r'0.\1', txt)     # " .004" → "0.004"
    txt = re.sub(r'-\.(\d)', r'-0.\1', txt)
    return ast.literal_eval(txt)


def sadline(tag):
    m = re.search(r'^' + tag + r'\s+(.*)$', s, re.M)
    return m.group(1)


m = re.search(r'^PBUNCH\s+([-+.\dE]+)\s+MINCOUP\s+([-+.\dE]+)', s, re.M)
PBUNCH = float(m.group(1))
MINCOUP = float(m.group(2))
print(f'PBUNCH={PBUNCH!r} MINCOUP={MINCOUP!r} radcod={radcod}')

sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
t0 = time.time()
em = emittance6d(lat, 3.0e9, rfsw=True, radcod=radcod)
print(f'emittance6d: {time.time() - t0:.1f} s  emx={em["emittance_x"]:.6e}')
t0 = time.time()
out = touschek_table(lat, 3.0e9, em, PBUNCH, mincoup=MINCOUP, rfsw=True, radcod=radcod)
print(f'touschek_table: {time.time() - t0:.1f} s  it={out["iterations"]} converged={out["converged"]}')
tab = out['table']


def rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = np.maximum(np.abs(b), 1e-300)
    return float(np.max(np.abs(a - b) / den))


# EMITX, EMITY, EMITZ, SIGE, SIGZ, EMITXC, EMITYC, EMITZC
ref = sadval('EMIT')
mine = [out['emittance_x'], out['emittance_y'], out['emittance_z'], out['momentum_spread'],
        out['bunch_length'], out['EMITXC'], out['EMITYC'], out['EMITZC']]
names = ['EMITX', 'EMITY', 'EMITZ', 'SIGE', 'SIGZ', 'EMITXC', 'EMITYC', 'EMITZC']
print('--- グローバル ---')
for n, a, b in zip(names, mine, ref):
    print(f'  {n:7s} pysad={a: .12e} SAD={b: .12e} rel={abs(a - b) / abs(b) if b else abs(a):.2e}')

print('--- TouschekTable ---')
print('  id:', tab[0], 'SAD', sadval('TT1'))
tt2 = sadval('TT2')
print(f'  [[2,1]] dp     max rel diff = {rel(tab[1][0], tt2[0]):.2e}')
print(f'  [[2,2]] 1/tau  max rel diff = {rel(tab[1][1], tt2[1]):.2e}')
print(f'  [[3,1]] ampl   max rel diff = {rel(tab[2][0], sadval("TT31")):.2e}')
print(f'  [[3,2]] mx     max rel diff = {rel(tab[2][1], sadval("TT32")):.2e}')
print(f'  [[3,3]] my     max rel diff = {rel(tab[2][2], sadval("TT33")):.2e}')
t4 = np.array(tab[3])
tt4s = sadval('TT4S')
print(f'  [[4]] 要素和 (120)   max rel diff = {rel(t4.sum(axis=0), tt4s):.2e}')
tt4e = np.array(sadval('TT4E'))
e4 = t4.sum(axis=1)
nz = tt4e != 0
print(f'  [[4]] dp 和 (nlat={len(tt4e)}, pysad {len(e4)})  max rel diff (非零 {nz.sum()} 要素) = '
      f'{rel(e4[nz], tt4e[nz]):.2e}  零要素の pysad 最大 = {np.max(np.abs(e4[~nz])) if (~nz).any() else 0:.2e}')
for tag, col in (('TT4C1', 0), ('TT4C2', 1), ('TT4C100', 99), ('TT4CL', -1)):
    print(f'  [[4,{col + 1 if col >= 0 else -1}]]  max rel diff = {rel(t4[col], sadval(tag)):.2e}')

print('--- TouschekLifetime ---')
INF = None
sige = out['momentum_spread']
cases = {
    'TL1': (INF, INF, 3), 'TL2': (10, INF, 3), 'TL3': (INF, 10, 3), 'TL4': (INF, INF, 0.03 / sige),
    'TL5': (3, INF, 8), 'TL6': (0.5, INF, 3), 'TL7': (INF, 0.5, 2), 'TL8': (INF, INF, 0.5),
    'TL9': (1000, INF, 200), 'TL10': (INF, 1000, 200), 'TL11': (INF, INF, 100), 'TL12': (INF, INF, 50),
    'TL13': (7.3, INF, 12.5), 'TL14': (INF, 7.3, 12.5), 'TL15': (INF, INF, 0),
}
for tag, (nx, ny, nzv) in cases.items():
    ref = float(sadline(tag).strip())
    v = touschek_lifetime(tab, nx, ny, nzv, sige)
    d = abs(v - ref) / abs(ref) if ref != 0 else abs(v)
    print(f'  {tag:5s} {str((nx, ny, nzv)):28s} pysad={v: .12e} SAD={ref: .12e} rel={d:.2e}')
