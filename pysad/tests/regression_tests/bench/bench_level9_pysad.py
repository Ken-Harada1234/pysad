#!/usr/bin/env python3
"""
bench_level9_pysad.py — Level 9 (ソレノイド領域の粒子追跡) の pysad 側
対応する SAD スクリプト: bench_level9.sad (参照値: bench_level9.result)

Level 6 のリングの CELL0 で B の後の D1 QD を SOL (BZ=2 T, BOUND=1) で挟む。
emit6d は SOL 未対応なので trf0 には SAD の表示値 (Equilibrium position dz) を使う。
5 粒子を 1, 100 周追跡 (RAD OFF / RAD ON NOFLUC)。
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
from pysad.track6 import track_particles

A = 0.157079632679490
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
PROG = int(os.environ.get("PROG", "0"))


def sol(name, bz, bound, geo=0.0):
    e = make_element('SOL', name, L=0.0)
    e.BZ = bz
    e.BOUND = bound
    e.GEO = geo
    return e


def cell(first):
    els = [make_element('MARK', 'ORG')]
    if first:
        ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
        ca.HARM = 234.0
        els.append(ca)
    els += [make_element('QUAD', 'QF', L=0.3, K1=1.8),
            make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SF', L=0.2, K2=0.6 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2)]
    if first:
        els.append(sol('SB', 2.0, 1.0, 1.0))
    els += [make_element('DRIFT', 'D1', L=1.2),
            make_element('QUAD', 'QD', L=0.3, K1=-1.7)]
    if first:
        els.append(sol('SE', 2.0, 1.0))
    els += [make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SD', L=0.2, K2=-1.2 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


els = cell(True)
for _ in range(19):
    els += cell(False)

s = open(os.path.join(HERE, 'bench_level9.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


trf0 = float(re.search(r'Equilibrium position dz =\s*' + NUM, s).group(1)) * 1e-3
print(f'trf0 = {trf0:.10e} [SAD の表示値 (emit6d は SOL 未対応)]', flush=True)

x0 = np.array([[1e-3, 0, 0, 0, 0, 0],
               [3e-3, 0, 0, 0, 0, 0],
               [0, 0, 1e-3, 0, 0, 0],
               [1e-3, 1e-4, 1e-3, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)

worst = 0.0
for rt, rad in (('N', False), ('R', True)):
    for nt in (1, 100):
        t0 = time.time()
        out = track_particles(els, x0, nt, 2.5e9, trf0=trf0, rfsw=True, rad=rad, progress=PROG)
        print(f'  S {rt} T{nt}  ({time.time() - t0:.1f} s)', flush=True)
        for i in range(5):
            ref = sadrow(f'S {rt} T{nt}', i + 1)
            d = out[i] - ref
            worst = max(worst, np.max(np.abs(d)))
            print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
                  + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
print(f'最大差 = {worst:.2e}')
