#!/usr/bin/env python3
"""
bench_level13a_pysad.py — Level 13A: マッチング (FIT/FREE/GO, pysad.match6 = tffsmatch の直訳) の小さなリング版
対応する SAD スクリプト: bench_level13a.sad (参照値: bench_level13a.result)
  13A-1: fit NX 7.55, NY 7.18; free QF QD; go        (SAD: 2 反復で Matched, 残差 5.08e-16)
  13A-2: + fit ORG BX 17.5, BY 1.0; free QF QD QF2    (SAD: 収束せず Residual 2.5977e-3)
"""
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.match6 import Matcher

NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
s = open(os.path.join(HERE, 'bench_level13a.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadvec(tag, n):
    m = re.search(r'^' + tag + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


sc = parser.load_sad_file(os.path.join(HERE, 'bench_level13a.sad'))
lat = sc.build_lattice('RING')
tab = lat.twiss_table()
org = tab[0]
py = [tab[-1].mux / 2 / np.pi, tab[-1].muy / 2 / np.pi, org.betx, org.bety]
print('TUNE0 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(py, sadvec('TUNE0', 4))))

print('=== 13A-1: fit NX NY, free QF QD ===')
m = Matcher(lat)
m.fit_clear()
m.fit('NX', 7.55)
m.fit('NY', 7.18)
m.free('QF QD'.replace(' ', ','))
m.go()
k1 = [m.values().get('QF K1'), m.values().get('QD K1'), lat._instances('QF2')[0].K1 * 0.3]
print('K1A pysad', k1)
print('K1A 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(k1, sadvec('K1A', 3))))
tab = lat.twiss_table()
py = [tab[-1].mux / 2 / np.pi, tab[-1].muy / 2 / np.pi]
print('TUNEA 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(py, sadvec('TUNEA', 2))))
print(f'RESA pysad {m.matching_residual:.4e} SAD {sadvec("RESA", 3)[0]:.4e}  反復 {m.iterations}')

print('=== 13A-2: fit NX NY and BX/BY at ORG, free QF QD QF2 ===')
m.fix('QF,QD')
m.fit_clear()
m.fit('NX', 7.55)
m.fit('NY', 7.18)
m.fit('BX', 17.5, loc='ORG')
m.fit('BY', 1.0, loc='ORG')
m.free('QF,QD,QF2')
m.go()
v = m.values()
k1 = [v['QF K1'], v['QD K1'], v['QF2 K1']]
print('K1B pysad', k1)
print('K1B 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(k1, sadvec('K1B', 3))))
tab = lat.twiss_table()
org = tab[0]
py = [tab[-1].mux / 2 / np.pi, tab[-1].muy / 2 / np.pi, org.betx, org.bety]
print('TWB 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(py, sadvec('TWB', 4))))
print(f'RESB pysad {m.matching_residual:.4e} SAD {sadvec("RESB", 3)[0]:.4e}  反復 {m.iterations}')
