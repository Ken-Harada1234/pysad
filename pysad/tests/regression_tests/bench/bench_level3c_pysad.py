#!/usr/bin/env python3
"""
bench_level3c_pysad.py — Level 3C (x-y 結合) の pysad 側
対応する SAD スクリプト: bench_level3c.sad (参照値: bench_level3c.result)

誤差は SAD が乱数で入れた値 (E 行) を読み、set_error で同じ 1 台だけに与える:
  QF.3 の ROTATE (DTHETA, 四極ロール) → pysad では ROTATION
  B.5  の DROTATE (DTHETA, BEND ロール誤差) → pysad では DROTATE
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
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'


def make_ring():
    els = []
    for _ in range(20):
        els += [make_element('MARK', 'ORG'),
                make_element('QUAD', 'QF', L=0.3, K1=1.8),
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
    return Lattice(els)


s = open(os.path.join(HERE, 'bench_level3c.result'), encoding='utf-8', errors='replace').read()


def one(pat):
    return float(re.search(pat + r'\s+' + NUM, s, re.M).group(1))


def vec(pat, n):
    m = re.search(pat + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


rot_qf3 = one(r'^\s*E ROTATE  QF\.3')
drot_b5 = one(r'^\s*E DROTATE B\.5 ')
print(f'誤差 (SAD の値): QF.3 ROTATE = {rot_qf3:.10e}, B.5 DROTATE = {drot_b5:.10e}')

lat = make_ring()
lat.twiss()
lat.set_error('QF.3', 'ROTATION', rot_qf3)
lat.set_error('B.5', 'DROTATE', drot_b5)
tw = lat.twiss_table()

NC = 13
print('\n--- 閉軌道 (x, px, y, py) ---')
for name, i in (('^^^', 0), ('ORG.11', 10 * NC)):
    r = tw[i]
    print(f'{name:7s} pysad ' + ' '.join(f'{v: .9e}' for v in (r.x_cod, r.px_cod, r.y_cod, r.py_cod)))
    print(f'{"":7s} SAD   ' + ' '.join(f'{v: .9e}' for v in vec(r'^\s*O ' + re.escape(name), 4)))

print('\n--- 結合ありの光学 (入口) ---')
r0 = tw[0]
for key, val in (('BX', r0.betx), ('BY', r0.bety), ('AX', r0.alfx), ('AY', r0.alfy),
                 ('EX', r0.etx), ('EPX', r0.etpx), ('EY', r0.ety), ('EPY', r0.etpy),
                 ('QX', tw[-1].mux / 2 / np.pi), ('QY', tw[-1].muy / 2 / np.pi)):
    ref = one(r'^\s*T ' + key)
    print(f'{key:3s}: pysad {val: .10f}  SAD {ref: .10f}  差 {val - ref: .2e}')
print('R1-R4 pysad ' + ' '.join(f'{v: .8e}' for v in (r0.r1, r0.r2, r0.r3, r0.r4)))
print('R1-R4 SAD   ' + ' '.join(f'{v: .8e}' for v in vec(r'^\s*T R ', 4)))

print('\n--- エミッタンス (Emittance[] = pysad.emit6d) ---')
from pysad.emit6d import emittance6d
e6 = emittance6d(lat, 2.5e9)
for key, val in (('EMITX', e6['emittance_x']), ('EMITY', e6['emittance_y']),
                 ('SIGE', e6['momentum_spread']), ('U0', e6['energy_loss_eV'])):
    ref = one(r'^\s*M ' + key)
    print(f'{key:5s}: pysad {val: .9e}  SAD {ref: .9e}  相対差 {(val - ref) / abs(ref): .2e}')
