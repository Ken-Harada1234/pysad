#!/usr/bin/env python3
"""
bench_level13b_pysad.py — Level 13B: 実ラティス (ESRF v03_12) のマッチング (pysad.match6 = tffsmatch の直訳)
対応する SAD スクリプト: bench_level13b.sad (参照値: bench_level13b.result)
  fit NX 47.60, NY 17.60, MKC BX 1.2 AX 0 BY 7.5 AY 0; free Q*; go
  (SAD: 2 反復 NEWTON で Matched, 残差 3.3724e-10)
各要素 (Q*) の K1 (積分値) と MKC/$$$ の Twiss を照合する。
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
from pysad.match6 import Matcher

NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
s = open(os.path.join(HERE, 'bench_level13b.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadlist(tag):
    i = s.index('\n' + tag + ' ')
    j = s.index('}', i)
    body = s[i + len(tag) + 2:j].replace('{', '')
    return [t.strip().strip('"') for t in body.split(',') if t.strip()]


def sadvec(tag, n):
    m = re.search(r'^' + tag + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
qels = [e for e in lat.elements if e.name.startswith('Q')]
qnames = sadlist('QNAMES')
assert [q.name for q in qels] == [n.split('.')[0] for n in qnames], 'Q* の要素列が違う'
k10 = [float(v) for v in sadlist('K10')]
print('K10 最大差:', max(abs(q.K1 * q.L - r) for q, r in zip(qels, k10)))


def tw():
    tab = lat.twiss_table()
    mkc = [r for r in tab if r.name.startswith('MKC')][0]
    return [tab[-1].mux / 2 / np.pi, tab[-1].muy / 2 / np.pi, mkc.betx, mkc.alfx, mkc.bety, mkc.alfy]


print('TW0 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(tw(), sadvec('TW0', 6))))
print('=== 13B: fit NX NY (ring) and BX/AX/BY/AY at MKC, free Q* ===')
m = Matcher(lat)
m.fit_clear()
m.fit('NX', 47.60)
m.fit('NY', 17.60)
m.fit('BX', 1.2, loc='MKC')
m.fit('AX', 0)
m.fit('BY', 7.5)
m.fit('AY', 0)
m.free('Q*')
print('変数:', [v['family'] for v in m.variables])
t0 = time.time()
m.go()
print(f'go {time.time() - t0:.1f} s, 反復 {m.iterations}')
k1b = [float(v) for v in sadlist('K1B')]
print('K1B 最大差:', max(abs(q.K1 * q.L - r) for q, r in zip(qels, k10 and k1b)))
print('TWB 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(tw(), sadvec('TWB', 6))))
print(f'RESB pysad {m.matching_residual:.5e} SAD {sadvec("RESB", 3)[0]:.5e}')
