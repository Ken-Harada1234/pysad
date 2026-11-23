#!/usr/bin/env python3
"""
bench_level11_pysad.py — Level 11 (MULT の K0/K1/K2 とずれ・回転, ANGLE=0) の pysad 側
対応する SAD スクリプト: bench_level11.sad (参照値: bench_level11.result)

Level 6 のリングで QF/QD → MULT K1、SF/SD → MULT K2 (MSF は K0=2E-4 付き)、D1 の先頭に MK0 (L=0.1 K0=1E-4)。
MQD に DX=1E-4、MSD に DY=1E-4 と ROTATE=1E-3。
emit (emittance6d) で trf0 を決めてから 5 粒子を 1, 100 周追跡 (RAD OFF / RAD ON NOFLUC)。
emit6d が未対応のときは SAD の表示値 (dz) を使う。
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
from pysad.track6 import track_particles

A = 0.157079632679490
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
PROG = int(os.environ.get("PROG", "0"))


def mult(name, L, knl, **kw):
    e = make_element('MULT', name, L=L, KnL=knl)
    for k, v in kw.items():
        setattr(e, k, v)
    return e


def cell(first):
    els = [make_element('MARK', 'ORG')]
    if first:
        ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
        ca.HARM = 234.0
        els.append(ca)
    els += [mult('MQF', 0.3, {1: 1.8 * 0.3}),
            make_element('DRIFT', 'DQS', L=0.1),
            mult('MSF', 0.2, {0: 2e-4, 2: 0.6}),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            mult('MK0', 0.1, {0: 1e-4}),
            make_element('DRIFT', 'D1S', L=1.1),
            mult('MQD', 0.3, {1: -1.7 * 0.3}, DX=1e-4),
            make_element('DRIFT', 'DQS', L=0.1),
            mult('MSD', 0.2, {2: -1.2}, DY=1e-4, ROTATION=1e-3),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


els = cell(True)
for _ in range(19):
    els += cell(False)

s = open(os.path.join(HERE, 'bench_level11.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


dz_sad = float(re.search(r'Equilibrium position dz =\s*' + NUM, s).group(1)) * 1e-3
try:
    r = emittance6d(Lattice(els), 2.5e9, rfsw=True)
    trf0 = r['trf0']
    src = 'pysad emit'
except NotImplementedError as e:
    trf0 = dz_sad
    src = f'SAD の表示値 (pysad emit 未対応: {e})'
print(f'trf0 = {trf0:.10e} [{src}]  (SAD dz = {dz_sad:.10e}, 差 {trf0 - dz_sad: .2e})', flush=True)

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
        print(f'  K {rt} T{nt}  ({time.time() - t0:.1f} s)', flush=True)
        for i in range(5):
            ref = sadrow(f'K {rt} T{nt}', i + 1)
            d = out[i] - ref
            worst = max(worst, np.max(np.abs(d)))
            print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
                  + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
print(f'最大差 = {worst:.2e}')
