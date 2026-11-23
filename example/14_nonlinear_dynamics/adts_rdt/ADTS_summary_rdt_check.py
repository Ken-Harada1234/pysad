# -*- coding: utf-8 -*-
"""ADTS_summary.md の RDT の節 (文献の式との照合)。小さな FODO (example/00_lattice/FODO_ring.sad、SF・SD の K2 = 4, -6 [積分値]、3 GeV、105.6 m) の
6 極を長さ 0 (薄レンズ) にして、Bengtsson 型の薄い 6 極の式 (h21000 = -(1/16) ΣK2L βx^(3/2) e^(iμx) など) と pysad の resonance_terms の比を出す。
比が -1 なら pysad の h は H の係数そのもの (文献は符号を反転した定義)。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/adts_rdt/ADTS_summary_rdt_check.py"""
import os, sys
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, ROOT)
from pysad import resdrive as R, compat
from pysad.parser import load_sad_string
compat.set_mode('v2')
src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0 K2=4) SD = (L=0 K2=-6)')
lat = load_sad_string(src).build_lattice('RING')
rows = lat.twiss_table()
T = ('h21000', 'h30000', 'h10110', 'h10020', 'h10200')
r = R.resonance_terms(lat, T, bend_cubic=False)      # 文献型の式は 6 極だけ (resonance_terms の既定は偏向磁石の 3 次の項つき。ここでは bend_cubic=False)
h = {t: 0j for t in T}
n = 0
for i, e in enumerate(lat.elements):
    if type(e).__name__ != 'Sextupole':
        continue
    n += 1
    k2l, rw = e.K2, rows[i]
    bx, by, mx, my = rw.betx, rw.bety, rw.mux, rw.muy
    h['h21000'] += -k2l * bx ** 1.5 * np.exp(1j * mx) / 16
    h['h30000'] += -k2l * bx ** 1.5 * np.exp(3j * mx) / 48
    h['h10110'] += k2l * bx ** .5 * by * np.exp(1j * mx) / 8
    h['h10020'] += k2l * bx ** .5 * by * np.exp(1j * (mx - 2 * my)) / 16
    h['h10200'] += k2l * bx ** .5 * by * np.exp(1j * (mx + 2 * my)) / 16
print('薄い 6 極 %d 台 (FODO、K2L = 4, -6)。h = 全周の和。比 = pysad / 文献型の式' % n)
print('%-8s %-28s %-28s %s' % ('項', 'pysad resonance_terms', '文献型の式', '比 (pysad / 文献)'))
for t in T:
    a, b = r.total[t], h[t]
    print('%-8s %+.9f %+.9fi   %+.9f %+.9fi   %+.8f %+.1ei' % (t, a.real, a.imag, b.real, b.imag, (a / b).real, (a / b).imag))
