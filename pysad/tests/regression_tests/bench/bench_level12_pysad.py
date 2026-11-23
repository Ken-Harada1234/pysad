#!/usr/bin/env python3
"""
bench_level12_pysad.py — Level 12: 実ラティス (ESRF 型 3GLS, sad/3GSR/esrf/20cells/v03_12, 3 GeV, 3564 要素) の pysad 側
対応する SAD スクリプト: bench_level12.sad (= pysad/tests/kek-ls_v03_12.lat + FFS 部; 参照値: bench_level12.result)

(1) ラティス読み込み (pysad.parser): 要素名・型・L・K1・ANGLE・K2・K3・E1・ROTATE を SAD の LINE[] と照合
(2) CELL CALC: 末尾 ($$$) と MKC の Twiss、チューン
(3) emit (RFSW): εx, σE, σz, α, U0, trf0, チューン, 減衰
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

NUM = r'([-+]?\d*\.?\d+(?:E[-+]?\d+)?)'
s = open(os.path.join(HERE, 'bench_level12.result'), encoding='utf-8', errors='replace').read()


def sadlist(tag):
    i = s.index('\n' + tag + ' ')
    j = s.index('}', i)
    body = s[i + len(tag) + 2:j].replace('{', '').replace('\n', ' ')
    return [t.strip().strip('"') for t in body.split(',') if t.strip()]


def sadvec(tag, n):
    m = re.search(r'^' + tag + r'\s*\{' + r',\s*'.join([r'\s*' + NUM] * n), s, re.M)
    return [float(v) for v in m.groups()]


def one(pat):
    return float(re.search(pat + r'\s*=?\s*' + NUM, s, re.M).group(1))


t0 = time.time()
sc = parser.load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'))
lat = sc.build_lattice('RING')
els = lat.elements
nlat_sad = int(float(re.search(r'^NLAT\s+' + NUM, s, re.M).group(1)))
print(f'ラティス読み込み {time.time() - t0:.1f} s: {len(els)} 要素, 周長 {sum(e.L for e in els):.6f} m '
      f'(SAD: {nlat_sad - 1} 要素, {one(r"^LENG"):.6f} m)')

names = sadlist('NAMES')
cols = {k: [float(v) for v in sadlist(k)] for k in ('LS', 'K1S', 'ANGLES', 'K2S', 'K3S', 'E1S', 'ROTS')}
g = lambda e, k: float(getattr(e, k, 0.0) or 0.0)
mism = 0
mx = {}
for i, e in enumerate(els):
    if e.name != names[i].split('.')[0]:
        mism += 1
    tn = type(e).__name__
    L = e.L
    d = {'L': abs(L - cols['LS'][i]), 'ANGLE': abs(g(e, 'ANGLE') - cols['ANGLES'][i]),
         'ROTATE': abs(g(e, 'ROTATION') + (g(e, 'ROTATE') if tn == 'Bend' else 0.0) - cols['ROTS'][i]) if tn != 'Drift' else 0.0,   # BEND の ROT は Bend.ROTATE (ANGLE=0 の VD も Bend のまま)
         'K1': abs(g(e, 'K1') * (L or 1) - cols['K1S'][i]) if tn in ('Quadrupole', 'Bend') else 0.0,
         'K2': abs(g(e, 'K2') * (L or 1) - cols['K2S'][i]) if tn == 'Sextupole' else 0.0,
         'K3': abs(g(e, 'K3') * (L or 1) - cols['K3S'][i]) if tn == 'Octupole' else 0.0,
         'E1': abs(g(e, 'E1') - cols['E1S'][i] * cols['ANGLES'][i]) if tn == 'Bend' else 0.0}
    for k, v in d.items():
        if v > mx.get(k, 0.0):
            mx[k] = v
print('要素名の不一致:', mism, ' パラメータの最大差 (SAD は 10 桁表示):', ' '.join(f'{k} {v:.1e}' for k, v in mx.items()))

t0 = time.time()
tab = lat.twiss_table()
print(f'CELL CALC (twiss_table) {time.time() - t0:.1f} s')
end = tab[-1]
ref = sadvec('TW_END', 8)
py = [end.betx, end.alfx, end.bety, end.alfy, end.etx, end.etpx, end.mux / 2 / np.pi, end.muy / 2 / np.pi]
print('  $$$ BX,AX,BY,AY,EX,EPX,NX,NY 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(py, ref)))
mkc = [r for r in tab if r.name.startswith('MKC')][0]
ref = sadvec('TW_MKC', 6)
py = [mkc.betx, mkc.alfx, mkc.bety, mkc.alfy, mkc.etx, mkc.etpx]
print('  MKC BX,AX,BY,AY,EX,EPX 差:', ' '.join(f'{a - b: .1e}' for a, b in zip(py, ref)))

t0 = time.time()
em = emittance6d(lat, 3.0e9, rfsw=True)
print(f'emit (RFSW) {time.time() - t0:.1f} s')
rows = [('εx', em['emittance_x'], one(r'^Emittance X')),
        ('σE', em['momentum_spread'], one(r'Energy spread')),
        ('σz', em['bunch_length'], one(r'^Bunch Length') * 1e-3),
        ('α', em['momentum_compaction'], one(r'Momentum compact\. alpha')),
        ('U0', em['energy_loss_eV'], one(r'^Energy loss per turn U0') * 1e6),
        ('trf0', em['trf0'], one(r'^Equilibrium position dz') * 1e-3)]
tun = sadvec('  Real tune:', 3) if re.search(r'^  Real tune:\s*\{', s, re.M) else None
m = re.search(r'Real tune:\s*' + NUM + r'\s+' + NUM + r'\s+' + NUM, s)
tun = [float(v) for v in m.groups()]
m = re.search(r'Damping per one revolution:\s*X :\s*' + NUM + r'\s+Y :\s*' + NUM + r'\s+Z :\s*' + NUM, s)
dmp = [float(v) for v in m.groups()]
rows += [('Qx', em['tunes'][0], tun[0]), ('Qy', em['tunes'][1], tun[1]), ('Qz', em['tunes'][2], tun[2]),
         ('damp x', em['damping_rate'][0], dmp[0]), ('damp y', em['damping_rate'][1], dmp[1]),
         ('damp z', em['damping_rate'][2], dmp[2])]
for name, v, r in rows:
    print(f'  {name:7s} pysad {v: .7e}  SAD {r: .7e}  相対差 {(v - r) / r: .1e}')
