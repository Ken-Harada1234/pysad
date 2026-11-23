# -*- coding: utf-8 -*-
"""v2 の Polarization の緩和率 (Sokolov-Ternov の厳密式) を、独立な DK の数値積分 (pysad/tests/script_tests/_p1_dk.py) と比べる。
リング: 0.5 GeV の 12 セル FODO + (a) 垂直の偏向の組 (±0.1 rad、ROTATE=90°)、(b) ソレノイド (BZ = 1 T、0.8 m)、(c) 四極の縦ずれ (垂直の分散・軌道)、
0.58 GeV の ASTRID2 + MPW の MPW を ROTATE = 0.3 rad / 90° にしたもの。使い方: py -3.12 example/16_spin_polarization/pol_st_dk_check.py (約 6 分)。結果は pol_st_dk_check.out.txt。"""
import io
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'pysad', 'tests', 'script_tests'))
from pysad import compat, ffscmd
import test_v2_P1_st as T
import _p1_dk as DK

compat.set_mode('v2')
print('A2 の Σ r (b·n)² の τ / 厳密な τ は、垂直の偏向・n の傾き・水平磁場の MPW で変わる。v2 = 厳密式。')
for name, line, kw in [('垂直の偏向 ±0.1 rad', 'LINE RING = (ST 6*CELLQ VDOG 6*CELLQ RFC);', dict(vb=0.1, qd=-0.55)),
                       ('ソレノイド 1 T', 'LINE RING = (ST SOLR 12*CELLQ RFC);', dict(bz=1.0, qd=-0.55)),
                       ('四極の縦ずれ (垂直の分散)', 'LINE RING = (ST 2*CELLQ CELLE1 4*CELLQ CELLE2 4*CELLQ RFC);', {})]:
    sp, r = T._both(T._fodo(line, **kw))
    print('%-24s τ v2 %.6g DK %.6g (比 %.6f) | A2 の式の τ / 厳密 = %.4f | P v2 %.5f DK %.5f (比 %.4f) | n の傾き 2/9<(n·v)²> = %.2g' %
          (name, sp['taup'], r['tau_st'], sp['taup'] / r['tau_st'], r['tau_b2'] / r['tau_st'], sp['equpol'][0], r['P_dk'], sp['equpol'][0] / r['P_dk'], r['two_ninth']))
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
for chi in (0.0, 0.3, math.pi / 2):
    src = open(path, encoding='utf-8').read() + '\nFFS;\nUSE RINGMPW;\nCELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    for e in ffs.lat.elements:
        if e.name.startswith('BMPW'):
            e.ROTATE = chi
    ffs.calc(quiet=True)
    sp, r = T._both(ffs)
    print('ASTRID2 + MPW ROTATE=%.4f τ v2 %.6g DK %.6g (比 %.6f) | A2 の式の τ / 厳密 = %.4f | P v2 %.5f DK %.5f (比 %.4f)' %
          (chi, sp['taup'], r['tau_st'], sp['taup'] / r['tau_st'], r['tau_b2'] / r['tau_st'], sp['equpol'][0], r['P_dk'], sp['equpol'][0] / r['P_dk']))
