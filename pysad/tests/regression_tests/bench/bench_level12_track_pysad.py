#!/usr/bin/env python3
"""
bench_level12_track_pysad.py — Level 12 (実ラティス ESRF v03_12) の粒子追跡
対応する SAD スクリプト: bench_level12_track.sad (参照値: bench_level12_track.result)

emit (RFSW) で trf0 を決め、5 粒子を 1 / 20 周: NORAD (N 行), RAD NOFLUC (R 行), RAD FLUC seed 17 (F 行)。
F では追跡後の乱数の状態 (seed) も比べる。
"""
import os
import re
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))   # リポジトリのルート
sys.path.insert(0, ROOT)
from pysad import parser
from pysad.emit6d import emittance6d
from pysad import track6
from pysad.track6 import track_particles

NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
PROG = int(os.environ.get("PROG", "0"))
s = open(os.path.join(HERE, 'bench_level12_track.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)


def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' + r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])


def sadseed(tag):
    return int(round(float(re.search(r'^SEED ' + tag + r' \{"SAD",\s*' + NUM, s, re.M).group(1))))


sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
els = lat.elements
em = emittance6d(lat, 3.0e9, rfsw=True)
trf0 = em['trf0']
dz_sad = float(re.search(r'Equilibrium position dz =\s*' + NUM, s).group(1)) * 1e-3
print(f'trf0 = {trf0:.10e} (SAD dz {dz_sad:.10e}, 差 {trf0 - dz_sad: .1e})', flush=True)

x0 = np.array([[1e-4, 0, 0, 0, 0, 0],
               [5e-4, 0, 0, 0, 0, 0],
               [0, 0, 1e-4, 0, 0, 0],
               [1e-4, 1e-5, 1e-4, 0, 0, 1e-3],
               [0, 0, 0, 0, 1e-3, 1e-3]], dtype=float)
worst = {}
for tag, rad, fluc in (('N', False, False), ('R', True, False), ('F', True, True)):
    w = 0.0
    for nt in (1, 20):
        t0 = time.time()
        out = track_particles(els, x0, nt, 3.0e9, trf0=trf0, rfsw=True, rad=rad, fluc=fluc,
                              seed=17 if fluc else None, progress=PROG)
        extra = ''
        if fluc:
            sp, ss = track6.RNG.seed % (1 << 32), sadseed(f'T{nt}')
            extra = f'  乱数の状態 pysad {sp} SAD {ss} ' + ('一致' if sp == ss else '不一致')
        print(f'  {tag} T{nt}  ({time.time() - t0:.1f} s){extra}', flush=True)
        for i in range(5):
            d = out[i] - sadrow(f'{tag} T{nt}', i + 1)
            w = max(w, np.max(np.abs(d)))
            print(f'    P{i + 1} 差: ' + ' '.join(f'{v: .1e}' for v in d))
    worst[tag] = f'{w:.1e}'
print('最大差:', worst)
