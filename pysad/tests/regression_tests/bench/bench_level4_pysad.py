#!/usr/bin/env python3
"""
bench_level4_pysad.py — Level 4 (RF 空洞と 6D) の pysad 側
対応する SAD スクリプト: bench_level4.sad (参照値: bench_level4.result)

比較するのは (1) RFSW ON, RADCOD OFF の Emittance[] (M1 行)。
(2) RADCOD ON (M2 行): 6D 閉軌道は Newton で解く (tcod の直訳ではない)。z と trf0 は z+trf0 で比較。
空洞: L=0, VOLT=1.5 MV, HARM=234, PHI=0 (SAD 既定)。pysad Cavity の PHI 既定値は π/2 なので明示する。
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
from pysad.emit6d import emittance6d

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

s = open(os.path.join(HERE, 'bench_level4.result'), encoding='utf-8', errors='replace').read()


def one(pat):
    return float(re.search(pat + r'\s+' + NUM, s, re.M).group(1))


def vec(pat, n):
    m = re.search(pat + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


r = emittance6d(lat, 2.5e9, rfsw=True)
emits = vec(r'^\s*M1 EMITS', 3)
tunes = vec(r'^\s*M1 TUNES', 3)
damp = vec(r'^\s*M1 DAMP', 3)
rows = [('εx', r['emittance_x'], emits[0]), ('εz', r['emittance_z'], emits[2]),
        ('σE', r['momentum_spread'], one(r'^\s*M1 SIGE')),
        ('σz', r['bunch_length'], one(r'^\s*M1 SIGZ')),
        ('Qx', r['tunes'][0], tunes[0]), ('Qy', r['tunes'][1], tunes[1]), ('Qz', r['tunes'][2], tunes[2]),
        ('damp x', r['damping_rate'][0], damp[0]), ('damp y', r['damping_rate'][1], damp[1]),
        ('damp z', r['damping_rate'][2], damp[2]),
        ('U0', r['energy_loss_eV'], one(r'^\s*M1 U0')),
        ('Vc', r['rf_voltage_eff'], one(r'^\s*M1 VC')),
        ('trf0', r['trf0'], one(r'^\s*M1 TRF0')),
        ('alpha', r['momentum_compaction'], one(r'^\s*M1 ALPHA')),
        ('bucket', r['bucket_height'], one(r'^\s*M1 BH'))]
print('--- (1) RFSW ON, RADCOD OFF ---')
for name, v, ref in rows:
    rel = (v - ref) / abs(ref) if ref != 0 else float('nan')
    print(f'{name:7s}: pysad {v: .9e}  SAD {ref: .9e}  相対差 {rel: .2e}')

r2 = emittance6d(lat, 2.5e9, rfsw=True, radcod=True)
emits = vec(r'^\s*M2 EMITS', 3)
tunes = vec(r'^\s*M2 TUNES', 3)
damp = vec(r'^\s*M2 DAMP', 3)
orb = vec(r'^\s*M2 ORB', 6)
trf0_sad = one(r'^\s*M2 TRF0')
co = r2['closed_orbit']
beta = 1.0  # z/β の β (2.5 GeV 電子なのでほぼ 1。比較は z + trf0 で行う)
rows = [('εx', r2['emittance_x'], emits[0]), ('εz', r2['emittance_z'], emits[2]),
        ('σE', r2['momentum_spread'], one(r'^\s*M2 SIGE')),
        ('σz', r2['bunch_length'], one(r'^\s*M2 SIGZ')),
        ('Qx', r2['tunes'][0], tunes[0]), ('Qy', r2['tunes'][1], tunes[1]), ('Qz', r2['tunes'][2], tunes[2]),
        ('damp x', r2['damping_rate'][0], damp[0]), ('damp z', r2['damping_rate'][2], damp[2]),
        ('U0', r2['energy_loss_eV'], one(r'^\s*M2 U0')),
        ('cod x', co[0], orb[0]), ('cod px', co[1], orb[1]), ('cod δ', co[5], orb[5]),
        ('z+trf0', co[4] + r2['trf0'], orb[4] + trf0_sad)]
print('\n--- (2) RFSW ON, RADCOD ON ---')
for name, v, ref in rows:
    rel = (v - ref) / abs(ref) if ref != 0 else float('nan')
    print(f'{name:7s}: pysad {v: .9e}  SAD {ref: .9e}  相対差 {rel: .2e}')
print(f'(参考) 閉軌道の z と trf0 の分け方は SAD と異なる: pysad z={co[4]:.6e}, trf0={r2["trf0"]:.9e} / '
      f'SAD z={orb[4]:.6e}, trf0={trf0_sad:.9e}')
