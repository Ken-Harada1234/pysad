#!/usr/bin/env python3
"""
bench_level3a_pysad.py — Level 3A (六極と色収差) の pysad 側
対応する SAD スクリプト: bench_level3a.sad (参照値: bench_level3a.result)

SAD と同じ手順:
  DP0 = -0.01..0.01 (0.001 刻み) で CELL 計算 → チューン → 3 次多項式で当てはめ
  (SAD 側は自作 MyPolyFit。ここでは numpy の最小二乗で同じことをする)
SEXT: SAD の K2 は積分値。pysad Sextupole.K2 は長さあたりなので K2_SAD/L を与える。
補正 (FIT) の比較はまだしない: SAD が求めた補正後の K2 を与えて、補正後の色収差を比べる。
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


def make_ring(k2f_sad, k2d_sad):
    els = []
    for _ in range(20):
        els += [make_element('MARK', 'ORG'),
                make_element('QUAD', 'QF', L=0.3, K1=1.8),
                make_element('SEXT', 'SF', L=L_S, K2=k2f_sad / L_S),
                make_element('DRIFT', 'DS', L=1.0),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2),
                make_element('QUAD', 'QD', L=0.3, K1=-1.8),
                make_element('SEXT', 'SD', L=L_S, K2=k2d_sad / L_S),
                make_element('DRIFT', 'DS', L=1.0),
                make_element('BEND', 'B', L=0.8, ANGLE=A, E1=A / 2, E2=A / 2),
                make_element('DRIFT', 'D1', L=1.2)]
    return Lattice(els)


def scan(lat):
    dps = np.round(np.arange(-0.01, 0.01 + 1e-12, 0.001), 6)
    q = []
    for dp in dps:
        tw = lat.twiss(delta=float(dp))
        q.append((tw[-1]['nx'] / 2 / np.pi, tw[-1]['ny'] / 2 / np.pi))
    lat.twiss(delta=0.0)
    q = np.array(q)
    fx = np.polynomial.polynomial.polyfit(dps, q[:, 0], 3)
    fy = np.polynomial.polynomial.polyfit(dps, q[:, 1], 3)
    return dps, q, fx, fy


def read_sad(path):
    s = open(path, encoding='utf-8', errors='replace').read()
    num = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
    out = {}
    for tag in ('B', 'A'):
        out[tag + 'Q'] = {float(m[0]): (float(m[1]), float(m[2])) for m in
                          re.findall(r'^\s*' + tag + r' Q \{\s*' + num + r',\s*' + num + r',\s*' + num + r'\}', s, re.M)}
        for ax in ('X', 'Y'):
            m = re.search(r'^\s*' + tag + ' FIT' + ax + r' \{\{\s*' + num + r',\s*' + num + r',\s*' + num + r',\s*' + num, s, re.M)
            out[tag + 'FIT' + ax] = [float(v) for v in m.groups()] if m else None
    for fam in ('SF', 'SD'):
        m = re.search(r'^\s*C K2 ' + fam + r'\s+' + num, s, re.M)
        out['K2' + fam] = float(m.group(1)) if m else None
    return out


sad = read_sad(os.path.join(HERE, 'bench_level3a.result'))

for tag, (kf, kd) in (('B', (0.6, -1.2)), ('A', (sad['K2SF'], sad['K2SD']))):
    title = '補正前' if tag == 'B' else '補正後 (SAD の補正結果 K2 を使用)'
    print('=' * 70)
    print(f'{title}: SF K2={kf}, SD K2={kd} (積分値)')
    print('=' * 70)
    lat = make_ring(kf, kd)
    dps, q, fx, fy = scan(lat)
    dq = []
    for dp, (qx, qy) in zip(dps, q):
        sq = sad[tag + 'Q'].get(round(float(dp), 6))
        if sq:
            dq.append((qx - sq[0], qy - sq[1]))
    dq = np.array(dq)
    print(f'  チューン 21 点の SAD との差 (最大): Qx {np.abs(dq[:, 0]).max():.2e}, Qy {np.abs(dq[:, 1]).max():.2e}')
    for ax, f in (('X', fx), ('Y', fy)):
        s = sad[tag + 'FIT' + ax]
        print(f'  {ax}: c0..c3 pysad = ' + ', '.join(f'{v: .10f}' for v in f))
        if s:
            print(f'     c0..c3 SAD   = ' + ', '.join(f'{v: .10f}' for v in s))
