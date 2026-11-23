#!/usr/bin/env python3
"""bench_level12_track_fast.py — track6n (numba) を track6 と SAD (bench_level12_track.result) に照合"""
import os, re, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..')); sys.path.insert(0, ROOT)   # ROOT = リポジトリのルート
from pysad import parser, track6
from pysad.emit6d import emittance6d
from pysad.track6n import track_particles_fast
NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
s = open(os.path.join(HERE, 'bench_level12_track.result'), encoding='utf-8', errors='replace').read()
s = re.sub(r',\s*\n\s*', ', ', s)
def sadrow(tag, i):
    m = re.search(r'^' + tag + r' P\s+' + str(i) + r'\.0+\s*\{' + r',\s*'.join([r'\s*' + NUM] * 7), s, re.M)
    return np.array([float(v) for v in m.groups()[:6]])
def sadseed(tag):
    return int(round(float(re.search(r'^SEED ' + tag + r' \{"SAD",\s*' + NUM, s, re.M).group(1))))
sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat')); lat = sc.build_lattice('RING'); els = lat.elements
em = emittance6d(lat, 3.0e9, rfsw=True); trf0 = em['trf0']
x0 = np.array([[1e-4,0,0,0,0,0],[5e-4,0,0,0,0,0],[0,0,1e-4,0,0,0],[1e-4,1e-5,1e-4,0,0,1e-3],[0,0,0,0,1e-3,1e-3]], dtype=float)
for tag, rad, fluc in (('N', False, False), ('R', True, False), ('F', True, True)):
    for nt in (1, 20):
        t0 = time.time()
        out = track_particles_fast(els, x0, nt, 3.0e9, trf0=trf0, rfsw=True, rad=rad, fluc=fluc, seed=17 if fluc else None)
        dt = time.time() - t0
        extra = ''
        if fluc:
            sp, ss = track6.RNG.seed % (1 << 32), sadseed(f'T{nt}')
            extra = f'  乱数 pysad {sp} SAD {ss} ' + ('一致' if sp == ss else '不一致')
        w = max(np.max(np.abs(out[i] - sadrow(f'{tag} T{nt}', i + 1))) for i in range(5))
        print(f'  {tag} T{nt:2d}  numba {dt:6.2f} s  SAD との最大差 {w:.1e}{extra}', flush=True)
# track6 (純 Python) との同一性
out_f = track_particles_fast(els, x0, 3, 3.0e9, trf0=trf0, rfsw=True, rad=True, fluc=True, seed=17)
out_s = track6.track_particles(els, x0, 3, 3.0e9, trf0=trf0, rfsw=True, rad=True, fluc=True, seed=17)
print('track6 vs track6n (RAD FLUC 3 turns) 最大差:', np.max(np.abs(out_f - out_s)))
t0 = time.time(); x1 = np.repeat(x0[:1], 500, axis=0) * np.linspace(0.1, 2.0, 500)[:, None]
track_particles_fast(els, x1, 10, 3.0e9, trf0=trf0, rfsw=True, rad=True, fluc=True); dt = time.time() - t0
print(f'500 粒子 × 10 周 (RAD FLUC): {dt:.1f} s → {dt / 500 / 10 / len(els) * 1e9:.0f} ns/粒子・要素')
