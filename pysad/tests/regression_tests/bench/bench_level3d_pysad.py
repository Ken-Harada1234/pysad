#!/usr/bin/env python3
"""
bench_level3d_pysad.py — Level 3D (強さの誤差 DK, DDK) の pysad 側
対応する SAD スクリプト: bench_level3d.sad (参照値: bench_level3d.result)

SAD が CALC 後に持つ値 (E 行) をそのまま 1 台ずつ与える:
  K1 (積分値) → pysad Quadrupole.K1 (長さあたり) = K1/L
  B.3 の K0 (BEND の追加の曲げ角, phib = ANGLE + K0) → pysad Bend.K0
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
LQ = 0.3
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'


def make_ring():
    els = []
    for _ in range(20):
        els += [make_element('MARK', 'ORG'),
                make_element('QUAD', 'QF', L=LQ, K1=1.8),
                make_element('DRIFT', 'DQS', L=0.1),
                make_element('SEXT', 'SF', L=0.2, K2=0.6 / 0.2),
                make_element('DRIFT', 'DS', L=0.9),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2),
                make_element('QUAD', 'QD', L=LQ, K1=-1.7),
                make_element('DRIFT', 'DQS', L=0.1),
                make_element('SEXT', 'SD', L=0.2, K2=-1.2 / 0.2),
                make_element('DRIFT', 'DS', L=0.9),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2)]
    return Lattice(els)


s = open(os.path.join(HERE, 'bench_level3d.result'), encoding='utf-8', errors='replace').read()


def one(pat):
    return float(re.search(pat + r'\s+' + NUM, s, re.M).group(1))


def vec(pat, n):
    m = re.search(pat + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


lat = make_ring()
lat.twiss()
lat.set_error('QF.2', 'K1', one(r'^\s*E K1 QF\.2') / LQ)
for m in re.finditer(r'^\s*E K1 QD\.\s*' + NUM + r'\s+' + NUM, s, re.M):
    i = int(round(float(m.group(1))))
    lat.set_error(f'QD.{i}', 'K1', float(m.group(2)) / LQ)
k0 = one(r'^\s*E K0 B\.3 ')
lat.set_error('B.3', 'K0', k0)
print(f'誤差: QF.2 K1={one(r"^\s*E K1 QF\.2"):.10f}, QD.1..20 の K1, B.3 K0={k0:.10e}')
tw = lat.twiss_table()

NC = 13
print('\n--- 閉軌道 (x, px, y, py) ---')
for name, i in (('^^^', 0), ('ORG.11', 10 * NC)):
    r = tw[i]
    print(f'{name:7s} pysad ' + ' '.join(f'{v: .9e}' for v in (r.x_cod, r.px_cod, r.y_cod, r.py_cod)))
    print(f'{"":7s} SAD   ' + ' '.join(f'{v: .9e}' for v in vec(r'^\s*O ' + re.escape(name), 4)))

print('\n--- 光学 ---')
r0 = tw[0]
r11 = tw[10 * NC]
for key, val in (('BX', r0.betx), ('BY', r0.bety), ('AX', r0.alfx), ('AY', r0.alfy),
                 ('EX', r0.etx), ('EPX', r0.etpx), ('BX11', r11.betx), ('BY11', r11.bety),
                 ('QX', tw[-1].mux / 2 / np.pi), ('QY', tw[-1].muy / 2 / np.pi)):
    ref = one(r'^\s*T ' + key)
    print(f'{key:4s}: pysad {val: .10f}  SAD {ref: .10f}  差 {val - ref: .2e}')

print('\n--- エミッタンス (Emittance[] = pysad.emit6d) ---')
from pysad.emit6d import emittance6d
e6 = emittance6d(lat, 2.5e9)
for key, val in (('EMITX', e6['emittance_x']), ('SIGE', e6['momentum_spread']),
                 ('U0', e6['energy_loss_eV'])):
    ref = one(r'^\s*M ' + key)
    print(f'{key:5s}: pysad {val: .9e}  SAD {ref: .9e}  相対差 {(val - ref) / abs(ref): .2e}')
print(f'EMITY: pysad {e6["emittance_y"]: .3e}  SAD {one(r"^\s*M EMITY"): .3e}  (結合なしなので 0 付近)')
