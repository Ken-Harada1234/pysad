#!/usr/bin/env python3
"""
bench_level5_pysad.py — Level 5 (輸送ライン INS とフリンジ) の pysad 側
対応する SAD スクリプト: bench_level5.sad (参照値: bench_level5.result)

ケース B (F1=0: ハードエッジのフリンジのみ) と C (DISFRIN=1) を比較する。
ケース A (F1=0.1, FRINGE=3: tqlfre, tblfre/tbendef1) も比較する。
入口の Twiss と軌道は先頭 MARK ORG の値 (SAD スクリプトと同じ)。
"""
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad.elements import make_element
from pysad.lattice import ins_twiss

NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
ANG = 0.1


def line(disfrin, f1=0.0, fringe=0):
    def q(name, k1):
        e = make_element('QUAD', name, L=0.5, K1=k1 / 0.5)
        e.DISFRIN = disfrin
        e.F1 = f1
        e.FRINGE = fringe
        return e

    def b():
        e = make_element('BEND', 'B', L=1.0, ANGLE=ANG, E1=0.5 * ANG, E2=0.5 * ANG)
        e.DISFRIN = disfrin
        e.F1 = f1
        e.FRINGE = fringe
        return e
    d1 = lambda: make_element('DRIFT', 'D1', L=1.0)
    d2 = lambda: make_element('DRIFT', 'D2', L=0.5)
    return [make_element('MARK', 'ORG'), d1(), q('QF', 0.4), d2(), b(), d2(), q('QD', -0.4),
            d1(), q('QF', 0.4), d2(), b(), d2(), q('QD', -0.4), d1(), make_element('MARK', 'EE')]


tw0 = dict(bx=10.0, ax=-1.0, by=5.0, ay=0.5, ex=0.1, epx=0.01)
cod0 = [1e-3, 1e-4, 0.5e-3, -1e-4]

s = open(os.path.join(HERE, 'bench_level5.result'), encoding='utf-8', errors='replace').read()


def vec(pat, n):
    m = re.search(pat + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


for tag, dis, f1, fr in (('LA', 0.0, 0.1, 3), ('LB', 0.0, 0.0, 0), ('LC', 1.0, 0.0, 0)):
    els = line(dis, f1, fr)
    out = ins_twiss(els, tw0, cod0)
    ee = out[len(els) - 1]      # EE の入口 = 最後の D1 の出口
    e6 = out[5]                 # 6 番目の要素 (D2) の入口
    print(f'=== {tag} ({"F1=0.1, FRINGE=3" if f1 else ("F1=0, フリンジ (tqfrie/tbedge) あり" if dis == 0 else "DISFRIN=1")}) ===')
    for lab, t, key in (('T EE', ee, r'T EE'), ('T 6 ', e6, r'T 6 ')):
        py = [t['bx'], t['ax'], t['by'], t['ay'], t['ex'], t['epx']]
        sad = vec(r'^\s*' + tag + ' ' + key, 6)
        print(f'  {lab} BX,AX,BY,AY,EX,EPX 差: ' + ' '.join(f'{a - b: .2e}' for a, b in zip(py, sad)))
    for lab, t, key in (('O EE', ee, r'O EE'), ('O 6 ', e6, r'O 6 ')):
        sad = vec(r'^\s*' + tag + ' ' + key, 4)
        c = t['cod']
        print(f'  {lab} x,px,y,py      差: ' + ' '.join(f'{a - b: .2e}' for a, b in zip(c[:4], sad)))
    sad = vec(r'^\s*' + tag + r' NX NY', 2)
    print(f'  NX,NY 差: {ee["nx"] - sad[0]: .2e} {ee["ny"] - sad[1]: .2e}  (SAD {sad[0]:.10f} {sad[1]:.10f})')
