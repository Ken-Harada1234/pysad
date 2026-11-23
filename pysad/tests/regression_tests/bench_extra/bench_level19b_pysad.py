#!/usr/bin/env python3
"""
bench_level19b_pysad.py — Level 19B: 幾何の分岐テスト (ROTATE 付き BEND, 垂直偏向, フリンジ長補正 (FRINGE/F1/FB1/FB2),
MARK OFFSET, 部分要素 (i+0.37), ORG) を pysad.geo で照合
対応する SAD スクリプト: bench_level19b.sad (参照値 bench_level19b.result)
使い方: PYTHONIOENCODING=utf-8 python -u bench_level19b_pysad.py
"""
import ast
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.geo import compute_geo

s = open(os.path.join(HERE, 'bench_level19b.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def _num(txt):
    txt = txt.replace('{', '[').replace('}', ']')
    txt = re.sub(r'(?<![\d.])\.(\d)', r'0.\1', txt)
    txt = re.sub(r'-\.(\d)', r'-0.\1', txt)
    return ast.literal_eval(txt)


ln = [(float(m.group(1)), _num(m.group(2))) for m in re.finditer(r'^LN\s+(\S+)\s+(\{.*\})\s*$', s, re.M)]
org = [(float(m.group(1)), _num(m.group(2))) for m in re.finditer(r'^ORG\s+(\S+)\s+(\{.*\})\s*$', s, re.M)]
pos_sad = _num(re.search(r'^POS\s+(.*)$', s, re.M).group(1))
names_sad = _num(re.search(r'^NAMES\s+(.*)$', s, re.M).group(1))

sc = parser.load_sad_file(os.path.join(HERE, 'bench_level19b.sad'))
els = sc.build_lattice('L1').elements
print('elements:', [e.name for e in els])
assert len(els) + 1 == len(names_sad)
g = compute_geo(els)
nlat = g.nlat
keys = ['GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3']

worst = 0.0
print('--- LINE[{GX..GCHI3}, i] / i+0.37 (ORG なし)')
for v, ref in ln:
    k = int(math.floor(v))
    fr = v - k
    mine = [g.line_value(key, k, fr) for key in keys]
    d = max(abs(a - b) for a, b in zip(mine, ref))
    worst = max(worst, d)
    flag = '' if d < 1e-10 else '   <-- NG'
    print(f'{v:7.2f} {names_sad[k - 1]:>6s} ' + ' '.join(f'{x:17.12f}' for x in mine) + f'  {d:.1e}{flag}')
    if d >= 1e-10:
        print(' ' * 15 + ' '.join(f'{x:17.12f}' for x in ref) + '  <- SAD')
dpos = max(abs(g.pos[k] - pos_sad[k]) for k in range(nlat))
print(f'POS max|diff| = {dpos:.1e}')
worst = max(worst, dpos)

# ---- ORG B3 1.5 -2.5 0.5 10 20 30 (角度は度)。
# SAD の getvad は "1.5 -2.5" を式 1.5-2.5 = -1.0 と読むので実際は gx0=-1.0, gy0=0.5, gz0=10, chi0=(20,30,[欠落→初期値 0]) 度
# となり "?ORG-Missing parameters" が出る (bench_level19c で確認)。ここではその解釈で照合する。
iorg = [e.name for e in els].index('B3') + 1
print(f'--- ORG (iorg={iorg}; gx0=-1.0 gy0=0.5 gz0=10, chi0=(20,30,0) deg)')
g2 = compute_geo(els, gx0=-1.0, gy0=0.5, gz0=10.0, chi0=(20 * math.pi / 180, 30 * math.pi / 180, 0.0), iorg=iorg)
d = 0.0
for v, ref in org:
    k = int(v)
    mine = [g2.line_value(key, k) for key in keys]
    dk = max(abs(a - b) for a, b in zip(mine, ref))
    d = max(d, dk)
    print(f'{v:7.2f} {names_sad[k - 1]:>6s} ' + ' '.join(f'{x:17.12f}' for x in mine) + f'  {dk:.1e}')
print(f'ORG max|diff| = {d:.1e}')
best = (d, None)
worst = max(worst, best[0])
print(f'\n全体の最大差 = {worst:.2e} (目標 1e-10) -> {"OK" if worst < 1e-10 else "NG"}')
