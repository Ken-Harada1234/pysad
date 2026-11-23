#!/usr/bin/env python3
"""
bench_level3b_pysad.py — Level 3B (六極の中を通る閉軌道) の pysad 側
対応する SAD スクリプト: bench_level3b.sad (参照値: bench_level3b.result)

誤差は SAD が乱数で入れた値 (E 行) を読み、set_error で同じ 1 台だけに与える。
比較: 閉軌道 (^^^, QD.1, SF.1, ORG.11 の入口), COD ありの光学, 色収差 (DP0 スキャン + 3 次当てはめ)。
"""
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad.elements import make_element
from pysad.lattice import Lattice

A = 0.157079632679490
L_S = 0.2
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'


def make_ring():
    els = []
    for _ in range(20):
        els += [make_element('MARK', 'ORG'),
                make_element('QUAD', 'QF', L=0.3, K1=1.8),
                make_element('DRIFT', 'DQS', L=0.1),
                make_element('SEXT', 'SF', L=L_S, K2=0.6 / L_S),
                make_element('DRIFT', 'DS', L=0.9),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2),
                make_element('QUAD', 'QD', L=0.3, K1=-1.8),
                make_element('DRIFT', 'DQS', L=0.1),
                make_element('SEXT', 'SD', L=L_S, K2=-1.2 / L_S),
                make_element('DRIFT', 'DS', L=0.9),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2)]
    return Lattice(els)


s = open(os.path.join(HERE, 'bench_level3b.result'), encoding='utf-8', errors='replace').read()


def one(pat):
    return float(re.search(pat + r'\s+' + NUM, s, re.M).group(1))


def vec(pat, n):
    m = re.search(pat + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


dx_qd1 = one(r'^\s*E DX QD\.1')
dy_sf1 = one(r'^\s*E DY SF\.1')
print(f'誤差 (SAD の値): QD.1 DX = {dx_qd1:.10e}, SF.1 DY = {dy_sf1:.10e}')

lat = make_ring()
lat.twiss()
lat.set_error('QD.1', 'DX', dx_qd1)
lat.set_error('SF.1', 'DY', dy_sf1)
tw = lat.twiss_table()

# 要素の入口の行 = table[要素の通し番号] (table[0] = リング入口)
NC = 13
idx = {'^^^': 0, 'QD.1': 7, 'SF.1': 3, 'ORG.11': 10 * NC}

print('\n--- 閉軌道 (x, px, y, py) ---')
for name, i in idx.items():
    r = tw[i]
    py_ = [r.x_cod, r.px_cod, r.y_cod, r.py_cod]
    sad = vec(r'^\s*O ' + re.escape(name), 4)
    print(f'{name:7s} pysad ' + ' '.join(f'{v: .9e}' for v in py_))
    print(f'{"":7s} SAD   ' + ' '.join(f'{v: .9e}' for v in sad))

print('\n--- COD ありの光学 (入口) ---')
r0 = tw[0]
for key, val in (('BX', r0.betx), ('BY', r0.bety), ('AX', r0.alfx), ('AY', r0.alfy),
                 ('EX', r0.etx), ('EY', r0.ety),
                 ('QX', tw[-1].mux / 2 / np.pi), ('QY', tw[-1].muy / 2 / np.pi)):
    ref = one(r'^\s*T ' + key)
    print(f'{key}: pysad {val: .10f}  SAD {ref: .10f}  差 {val - ref: .2e}')
rr = [r0.r1, r0.r2, r0.r3, r0.r4]
print('R1-R4 pysad ' + ' '.join(f'{v: .6e}' for v in rr))
print('R1-R4 SAD   ' + ' '.join(f'{v: .6e}' for v in vec(r'^\s*T R', 4)))

print('\n--- COD ありの色収差 ---')
dps = np.round(np.arange(-0.01, 0.01 + 1e-12, 0.001), 6)
q = []
for dp in dps:
    t = lat.twiss(delta=float(dp))
    q.append((t[-1]['nx'] / 2 / np.pi, t[-1]['ny'] / 2 / np.pi))
q = np.array(q)
for ax, col in (('X', 0), ('Y', 1)):
    f = np.polynomial.polynomial.polyfit(dps, q[:, col], 3)
    sad = vec(r'^\s*FIT' + ax + r' \{', 4)
    print(f'{ax}: c0..c3 pysad ' + ', '.join(f'{v: .10f}' for v in f))
    print(f'{ax}: c0..c3 SAD   ' + ', '.join(f'{v: .10f}' for v in sad))
