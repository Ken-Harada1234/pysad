#!/usr/bin/env python3
"""
bench_level6_pysad.py — Level 6 (粒子追跡 TrackParticles) の pysad 側
対応する SAD スクリプト: bench_level6.sad (参照値: bench_level6.result)

ラティスは Level 4 と同じ (CELL0 に空洞)。emit で trf0 を決めてから 5 粒子を 1/10/100 周追跡。
(1) RAD OFF (N1 行)  (2) RAD ON, FLUC OFF (R1 行)
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
lat = Lattice(els)

s = open(os.path.join(HERE, 'bench_level6.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)          # SAD の折り返し行をつなぐ


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


r = emittance6d(lat, 2.5e9, rfsw=True)
trf0 = r['trf0']
print(f'trf0 = {trf0:.12e}', flush=True)

x0 = np.array([[1e-3, 0, 0, 0, 0, 0],
               [3e-3, 0, 0, 0, 0, 0],
               [0, 0, 1e-3, 0, 0, 0],
               [1e-3, 1e-4, 1e-3, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)

PROG = int(os.environ.get("PROG", "10"))
names = ['x', 'px', 'y', 'py', 'z', 'dp']
for tag, rad in (('N1', False), ('R1', True)):
    print(f'=== {tag}: RAD {"ON (FLUC OFF)" if rad else "OFF"} ===')
    for nt in (1, 10, 100):
        t0 = time.time()
        out = track_particles(els, x0, nt, 2.5e9, trf0=trf0, rfsw=True, rad=rad, progress=PROG)
        print(f'  T{nt}  ({time.time() - t0:.1f} s)')
        for i in range(5):
            ref = sadrow(f'{tag} T{nt}', i + 1)
            d = out[i] - ref
            print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
                  + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
