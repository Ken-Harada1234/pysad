#!/usr/bin/env python3
"""
bench_level7_pysad.py — Level 7 (誤差入り要素と K1≠0 BEND の粒子追跡) の pysad 側
対応する SAD スクリプト: bench_level7.sad (参照値: bench_level7.result)

7-1 (K 行): BEND を K1=-0.05 (SAD の積分値) の BK に替えたリング。誤差なし。
7-2 (E 行): Level 6 のリングに SAD が乱数で入れた誤差 (ERR 行の値) を set_error で同じ 1 台に入れる。
どちらも emit (emittance6d) で trf0 を決めてから 5 粒子を 1, 100 周追跡 (RAD OFF / RAD ON NOFLUC)。
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


def cell(with_cav, bk1=0.0):
    els = [make_element('MARK', 'ORG')]
    if with_cav:
        ca = make_element('CAVI', 'CA', L=0.0, VOLT=1.5e6, PHI=0.0)
        ca.HARM = 234.0
        els.append(ca)
    bname = 'BK' if bk1 else 'B'

    def b():
        return make_element('BEND', bname, L=0.8, ANGLE=A, E1=A / 2, E2=A / 2, K1=bk1 / 0.8)
    els += [make_element('QUAD', 'QF', L=0.3, K1=1.8),
            make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SF', L=0.2, K2=0.6 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            b(),
            make_element('DRIFT', 'D1', L=1.2),
            make_element('QUAD', 'QD', L=0.3, K1=-1.7),
            make_element('DRIFT', 'DQS', L=0.1),
            make_element('SEXT', 'SD', L=0.2, K2=-1.2 / 0.2),
            make_element('DRIFT', 'DS', L=0.9),
            b(),
            make_element('DRIFT', 'D1', L=1.2)]
    return els


def ring(bk1=0.0):
    els = cell(True, bk1)
    for _ in range(19):
        els += cell(False, bk1)
    return els


s = open(os.path.join(HERE, 'bench_level7.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)          # SAD の折り返し行をつなぐ


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' +
                  r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


def one(pat):
    return float(re.search(pat + r'\s*' + NUM, s, re.M).group(1))


dz_sad = [float(v) for v in re.findall(r'Equilibrium position dz =\s*' + NUM, s)]

x0 = np.array([[1e-3, 0, 0, 0, 0, 0],
               [3e-3, 0, 0, 0, 0, 0],
               [0, 0, 1e-3, 0, 0, 0],
               [1e-3, 1e-4, 1e-3, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)

# 7-1: K1≠0 BEND
els_k = ring(-0.05)
# 7-2: 誤差入り
els_e = ring(0.0)
lat_e = Lattice(els_e)
errs = [('QD.1', 'DX', one(r'^ERR QD\.1 DX ')), ('SF.1', 'DY', one(r'^ERR SF\.1 DY ')),
        ('QF.3', 'ROTATION', one(r'^ERR QF\.3 ROTATE ')), ('B.5', 'DROTATE', one(r'^ERR B\.5 DROTATE ')),
        ('B.2', 'DX', one(r'^ERR B\.2 DX '))]
for inst, key, v in errs:
    lat_e.set_error(inst, key, v)
print('誤差 (SAD の値): ' + ', '.join(f'{i} {k} = {v:.10e}' for i, k, v in errs))

worst = 0.0
for prefix, els, idz in (('K', els_k, 0), ('E', els_e, 1)):
    lat = Lattice(els) if prefix == 'K' else lat_e
    r = emittance6d(lat, 2.5e9, rfsw=True)
    trf0 = r['trf0']
    print(f'=== {prefix}: trf0 = {trf0:.10e}  (SAD dz = {dz_sad[idz] * 1e-3:.10e}, 差 {trf0 - dz_sad[idz] * 1e-3: .2e})',
          flush=True)
    for rt, rad in (('N', False), ('R', True)):
        for nt in (1, 100):
            t0 = time.time()
            out = track_particles(els, x0, nt, 2.5e9, trf0=trf0, rfsw=True, rad=rad, progress=PROG)
            print(f'  {prefix} {rt} T{nt}  ({time.time() - t0:.1f} s)', flush=True)
            for i in range(5):
                ref = sadrow(f'{prefix} {rt} T{nt}', i + 1)
                d = out[i] - ref
                worst = max(worst, np.max(np.abs(d)))
                print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .2e}' for v in d)
                      + f'   |SAD| max {np.max(np.abs(ref)):.2e}')
print(f'最大差 = {worst:.2e}')
