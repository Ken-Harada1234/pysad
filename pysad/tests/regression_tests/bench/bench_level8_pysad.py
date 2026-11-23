#!/usr/bin/env python3
"""
bench_level8_pysad.py — Level 8 (追跡: 四極フリンジ、長さのある空洞、MULT) の pysad 側
対応する SAD スクリプト: bench_level8.sad (参照値: bench_level8.result)

8-1 (F 行): 四極に F1=0.05, FRINGE=3
8-2 (C 行): 空洞 L=0.5 (CELLC0 の先頭)。QF K1=1.76*0.3 (半整数共鳴を避けるため)
8-3 (M 行): 六極 SF/SD を MULT (K2 のみ) に
どれも emit (emittance6d) で trf0 を決めてから 5 粒子を 1, 100 周追跡 (RAD OFF / RAD ON NOFLUC)。
環境変数 CASES (例 "F,M") で実行するケースを選べる。PROG で進捗表示。
"""
import os
import re
import sys
import time
import traceback
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
CASES = os.environ.get("CASES", "F,C,M").split(',')


def cell(first, kind):
    els = [make_element('MARK', 'ORG')]
    if first:
        if kind == 'C':
            ca = make_element('CAVI', 'CAL', L=0.5, VOLT=1.5e6, PHI=0.0)
        else:
            ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
        ca.HARM = 234.0
        els.append(ca)

    def quad(name, k1):
        q = make_element('QUAD', name, L=0.3, K1=k1)
        if kind == 'F':
            q.F1 = 0.05
            q.FRINGE = 3
        return q

    def sext(name, k2):
        if kind == 'M':
            return make_element('MULT', name, L=0.2, KnL={2: k2})
        return make_element('SEXT', name, L=0.2, K2=k2 / 0.2)
    kqf = 1.76 if kind == 'C' else 1.8
    els += [quad('QF', kqf),
            make_element('DRIFT', 'DQS', L=0.1),
            sext('SF', 0.6),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2),
            quad('QD', -1.7),
            make_element('DRIFT', 'DQS', L=0.1),
            sext('SD', -1.2),
            make_element('DRIFT', 'DS', L=0.9),
            make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


def ring(kind):
    els = cell(True, kind)
    for _ in range(19):
        els += cell(False, kind)
    return els


s = open(os.path.join(HERE, 'bench_level8.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)          # SAD の折り返し行をつなぐ


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


dz_sad = [float(v) for v in re.findall(r'Equilibrium position dz =\s*' + NUM, s)]
IDX = {'F': 0, 'C': 1, 'M': 2}

x0 = np.array([[1e-3, 0, 0, 0, 0, 0],
               [3e-3, 0, 0, 0, 0, 0],
               [0, 0, 1e-3, 0, 0, 0],
               [1e-3, 1e-4, 1e-3, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)

worst = {}
for kind in CASES:
    els = ring(kind)
    try:
        r = emittance6d(Lattice(els), 2.5e9, rfsw=True)
        trf0 = r['trf0']
        src = 'pysad emit'
    except NotImplementedError as e:
        trf0 = dz_sad[IDX[kind]] * 1e-3
        src = f'SAD の表示値 (pysad emit 未対応: {e})'
    print(f'=== {kind}: trf0 = {trf0:.10e} [{src}]  (SAD dz = {dz_sad[IDX[kind]] * 1e-3:.10e})', flush=True)
    w = 0.0
    try:
        for rt, rad in (('N', False), ('R', True)):
            for nt in (1, 100):
                t0 = time.time()
                out = track_particles(els, x0, nt, 2.5e9, trf0=trf0, rfsw=True, rad=rad, progress=PROG)
                print(f'  {kind} {rt} T{nt}  ({time.time() - t0:.1f} s)', flush=True)
                for i in range(5):
                    ref = sadrow(f'{kind} {rt} T{nt}', i + 1)
                    d = out[i] - ref
                    w = max(w, np.max(np.abs(d)))
                    print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
                          + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
        worst[kind] = f'{w:.2e}'
    except NotImplementedError as e:
        print(f'  未対応: {e}')
        worst[kind] = '未対応'
print('最大差:', worst)
