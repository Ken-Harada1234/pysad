#!/usr/bin/env python3
"""
bench_level19a_pysad.py — Level 19A: 要素幾何 (GX/GY/GZ/CHI1..3) の照合 (pysad.geo = tfgeo/tfgeo1/geolib の直訳)
対応する SAD スクリプト: bench_level19a.sad (pysad/tests/kek-ls_v03_12.lat + FFS; USE RING; CELL; CALC;
  LINE[{"GX","GY","GZ","GCHI1","GCHI2","GCHI3"}, loc] を 21 箇所、LINE["GEO",...]、LINE["S",...])
  参照値: bench_level19a.result
使い方: PYTHONIOENCODING=utf-8 python -u bench_level19a_pysad.py
目標: 1e-10 (絶対差)
"""
import ast
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.geo import compute_geo

res = os.path.join(HERE, 'bench_level19a.result')
s = open(res, encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def _num(txt):
    txt = txt.replace('{', '[').replace('}', ']')
    txt = re.sub(r'(?<![\d.])\.(\d)', r'0.\1', txt)
    txt = re.sub(r'-\.(\d)', r'-0.\1', txt)
    return ast.literal_eval(txt)


# ---- SAD の出力を読む
ln_rows = []          # (loc 文字列, [6 値])
for m in re.finditer(r'^LN\s+(\S+)\s+(\{.*\})\s*$', s, re.M):
    ln_rows.append((m.group(1), _num(m.group(2))))
geo2 = _num(re.search(r'^GEO2\s+(.*)$', s, re.M).group(1))
geo3 = _num(re.search(r'^GEO3\s+(.*)$', s, re.M).group(1))
pos_sad = _num(re.search(r'^POS\s+(.*)$', s, re.M).group(1))
nlat_sad = int(float(re.search(r'^NLAT\s+(.*)$', s, re.M).group(1)))

# ---- pysad: ラティスを読み、幾何を計算
sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
els = lat.elements
g = compute_geo(els)
nlat = g.nlat
print(f'nlat pysad={nlat} SAD={nlat_sad}')
assert nlat == nlat_sad


def resolve(loc: str):
    """SAD の位置指定 → (k, fr)。'$$$' は nlat、'^^^' は 1、'NAME.n' は n 台目、数値は整数部+小数部"""
    if loc == '$$$':
        return nlat, 0.0
    if loc == '^^^':
        return 1, 0.0
    try:
        v = float(loc)
        k = int(np.floor(v))
        if k > nlat:
            k = k - nlat                       # SAD の LINE[] は nlat を超える番号を周回させる (4000 -> 436 = 4000 - nlat)
        k = max(1, k)
        return k, v - np.floor(v) if k < nlat else 0.0
    except ValueError:
        pass
    m = re.match(r'^([^.]+)(?:\.(\d+))?$', loc)
    name, ordn = m.group(1), int(m.group(2) or 1)
    cnt = 0
    for i, el in enumerate(els):
        if el.name == name:
            cnt += 1
            if cnt == ordn:
                return i + 1, 0.0
    raise KeyError(loc)


keys = ['GX', 'GY', 'GZ', 'GCHI1', 'GCHI2', 'GCHI3']
worst = 0.0
print(f'{"loc":>18s} {"k":>5s}  ' + '  '.join(f'{k:>18s}' for k in keys) + '   max|diff|')
for loc, vals in ln_rows:
    k, fr = resolve(loc)
    mine = [g.line_value(key, k, fr) for key in keys]
    d = max(abs(a - b) for a, b in zip(mine, vals))
    worst = max(worst, d)
    print(f'{loc:>18s} {k:5d}  ' + '  '.join(f'{v:18.12f}' for v in mine) + f'   {d:.2e}')
    if d > 1e-9:
        print(' ' * 26 + '  '.join(f'{v:18.12f}' for v in vals) + '   <- SAD')

# LINE["GEO", ...] (位置 3 成分 + CHI 3 成分)
for tag, ref, loc in (('GEO2', geo2, 'BPI6.3'), ('GEO3', geo3, '$$$')):
    k, fr = resolve(loc)
    mine = g.line_value('GEO', k, fr)
    d = max(abs(a - b) for a, b in zip(mine[0] + mine[1], ref[0] + ref[1]))
    worst = max(worst, d)
    print(f'{tag} {loc}: pysad {mine}  diff {d:.2e}')

# LINE["S", locs] (入口の経路長 pos)
locs = [loc for loc, _ in ln_rows]
dpos = 0.0
for loc, ps in zip(locs, pos_sad):
    k, fr = resolve(loc)
    dpos = max(dpos, abs(g.pos_frac(k, fr) - ps))
print(f'POS (LINE["S"]) max|diff| = {dpos:.2e}')
worst = max(worst, dpos)

print(f'\n全体の最大差 = {worst:.2e}  (目標 1e-10) -> {"OK" if worst < 1e-10 else "NG"}')
