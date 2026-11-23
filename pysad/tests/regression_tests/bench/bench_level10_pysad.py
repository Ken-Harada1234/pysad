#!/usr/bin/env python3
"""
bench_level10_pysad.py — Level 10 (放射の揺らぎ RAD+FLUC の粒子追跡) の pysad 側
対応する SAD スクリプト: bench_level10.sad (参照値: bench_level10.result)

Level 6 のリング。emit で trf0 を決め、SeedRandom[17] 相当 (seed=17) で 5 粒子を 1, 10, 100 周追跡。
SAD の乱数プラグイン "SAD" (合同法) を pysad/sadrandom.py に直訳し、tradkfn/tdusrn/drndsr で
SAD と同じ順に乱数を引く。追跡後の乱数の内部状態 (seed) も SAD の SeedRandom[] の値と比べる。
"""
import os
import re
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad.elements import make_element
from pysad.lattice import Lattice
from pysad.emit6d import emittance6d
from pysad import track6
from pysad.track6 import track_particles

A = 0.157079632679490
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
PROG = int(os.environ.get("PROG", "0"))


def cell(with_cav):
    els = [make_element('MARK', 'ORG')]
    if with_cav:
        ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
        ca.HARM = 234.0
        els.append(ca)
    els += [make_element('QUAD', 'QF', L=0.3, K1=1.8),
            make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SF', L=0.2, K2=0.6 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2),
            make_element('QUAD', 'QD', L=0.3, K1=-1.7),
            make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SD', L=0.2, K2=-1.2 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


els = cell(True)
for _ in range(19):
    els += cell(False)

s = open(os.path.join(HERE, 'bench_level10.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


def sadseed(tag):
    m = re.search(r'^SEED ' + tag + r' \{"SAD",\s*' + NUM, s, re.M)
    return int(round(float(m.group(1))))


r = emittance6d(Lattice(els), 2.5e9, rfsw=True)
trf0 = r['trf0']
print(f'trf0 = {trf0:.10e}', flush=True)

x0 = np.array([[1e-3, 0, 0, 0, 0, 0],
               [3e-3, 0, 0, 0, 0, 0],
               [0, 0, 1e-3, 0, 0, 0],
               [1e-3, 1e-4, 1e-3, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)

worst = 0.0
for nt in (1, 10, 100):
    t0 = time.time()
    out = track_particles(els, x0, nt, 2.5e9, trf0=trf0, rfsw=True, rad=True, fluc=True, seed=17,
                          progress=PROG)
    seed_py = track6.RNG.seed % (1 << 32)
    seed_sad = sadseed(f'T{nt}')
    print(f'  F T{nt}  ({time.time() - t0:.1f} s)  乱数の状態: pysad {seed_py}  SAD {seed_sad}  '
          + ('一致' if seed_py == seed_sad else '不一致'), flush=True)
    for i in range(5):
        ref = sadrow(f'F T{nt}', i + 1)
        d = out[i] - ref
        worst = max(worst, np.max(np.abs(d)))
        print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
              + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
print(f'最大差 = {worst:.2e}')
