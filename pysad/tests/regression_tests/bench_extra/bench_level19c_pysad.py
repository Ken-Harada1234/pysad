#!/usr/bin/env python3
"""
bench_level19c_pysad.py — Level 19C: ORG コマンド (原点・向きの付け替え: tfgeo の末尾の変換) の照合
対応する SAD スクリプト: bench_level19c.sad (参照値 bench_level19c.result)
SAD の getvad は "1.5 -2.5 0.5" を "1.5-2.5" "0.5" と式として読むので、各 ORG 行の実際の引数は
  A: ORG B3 1.5 -2.5 0.5 0 0 0   → gx0=-1.0 gy0=0.5 gz0=0  chi0=(0,0,0)       (chi3 欠落 → 初期値 0)
  B: ORG B3 0 0 0 10 0 0         → chi0=(10,0,0) 度
  C: ORG B3 0 0 0 0 20 0         → chi0=(0,20,0)
  D: ORG B3 0 0 0 0 0 30         → chi0=(0,0,30)
  E: ORG B3 1.5 -2.5 0.5 10 20 30 → gx0=-1.0 gy0=0.5 gz0=10 chi0=(20,30,0)  (chi3 欠落 → 0: 直前の D の 30 は残らない)
  F: ORG D2.1 ...                → 同上、iorg=D2.1
(tffsa.f の chi0(i)=getvad(exist,word,chi0(i)) は既定値が chi0(i) だが、実測では欠落分は 0 になった)
使い方: PYTHONIOENCODING=utf-8 python -u bench_level19c_pysad.py
"""
import ast
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.geo import compute_geo

s = open(os.path.join(HERE, 'bench_level19c.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def _num(txt):
    txt = txt.replace('{', '[').replace('}', ']')
    txt = re.sub(r'(?<![\d.])\.(\d)', r'0.\1', txt)
    txt = re.sub(r'-\.(\d)', r'-0.\1', txt)
    return ast.literal_eval(txt)


sc = parser.load_sad_file(os.path.join(HERE, 'bench_level19c.sad'))
els = sc.build_lattice('L1').elements
names = [e.name for e in els]
keys = ['GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3']
D = math.pi / 180
cases = {
    'A': dict(iorg=names.index('B3') + 1, gx0=-1.0, gy0=0.5, gz0=0.0, chi0=(0, 0, 0)),
    'B': dict(iorg=names.index('B3') + 1, chi0=(10 * D, 0, 0)),
    'C': dict(iorg=names.index('B3') + 1, chi0=(0, 20 * D, 0)),
    'D': dict(iorg=names.index('B3') + 1, chi0=(0, 0, 30 * D)),
    'E': dict(iorg=names.index('B3') + 1, gx0=-1.0, gy0=0.5, gz0=10.0, chi0=(20 * D, 30 * D, 0.0)),
    'F': dict(iorg=names.index('D2') + 1, gx0=-1.0, gy0=0.5, gz0=10.0, chi0=(20 * D, 30 * D, 0.0)),
}
worst = 0.0
for tag, kw in cases.items():
    rows = [(int(float(m.group(1))), _num(m.group(2))) for m in re.finditer(r'^ORG' + tag + r'\s+(\S+)\s+(\{.*\})\s*$', s, re.M)]
    g = compute_geo(els, **kw)
    d = 0.0
    for k, ref in rows:
        mine = [g.line_value(key, k) for key in keys]
        d = max(d, max(abs(a - b) for a, b in zip(mine, ref)))
    worst = max(worst, d)
    print(f'ORG {tag}: {kw}  max|diff| = {d:.1e}  ({len(rows)} 点)')
print(f'\n全体の最大差 = {worst:.2e} (目標 1e-10) -> {"OK" if worst < 1e-10 else "NG"}')
