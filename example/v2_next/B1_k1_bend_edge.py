# -*- coding: utf-8 -*-
"""v2next B1 項目 2: K1≠0 の BEND の CALC (tbende) と追跡 (tbendi) の写像の差。どちらが正しいかを、独立な「厳密な写像」(ローレンツ力の数値積分) との
比較で出す。
条件: 3 GeV、L = 1.2 m、ANGLE = 15 度、K1 = -0.05 (積分値)、端の角 E1 = E2 = e·ANGLE (e = 0.5: 矩形、0.25)。写像は (x, px)、軌道は (x0, px0, δ)。
結論: CALC は厳密な写像と 1e-11 で一致。K64 の追跡 (tbendi) は薄い端の写像で 1e-9〜1e-8 ずれる (誤りは追跡の側)。v2 の追跡は CALC と同じ厳密な端 (tbedge) を使い 1e-11。
使い方: py -3.12 example/v2_next/B1_k1_bend_edge.py   (v2 と k64 の両方の追跡の誤差を出す)"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pysad import compat
from pysad.tests.test_v2_B1_optics import _exact_map, _k1_bend, _calc_and_track, _TH, _ORBITS

print('軌道 (x0, px0, δ)        |CALC − 厳密|   |追跡 v2 − 厳密|   |追跡 k64 − 厳密|   (x, px の最大の差)')
for e in (0.5, 0.25):
    print('--- e = %s (E1 = E2 = %.4f rad) ---' % (e, e * _TH))
    rows = {}
    for mode in ('v2', 'k64'):
        compat.set_mode(mode)
        el, ffs = _k1_bend(-0.05, e)
        for c3 in _ORBITS:
            c = [c3[0], c3[1], 0, 0, 0, c3[2]]
            cc, _, t = _calc_and_track(el, ffs, c)
            ex = _exact_map(c, el.K1, e * _TH, e * _TH)
            rows.setdefault(c3, {})[mode] = (np.abs(cc[:2] - ex).max(), np.abs(t[:2] - ex).max())
    for c3, r in rows.items():
        print('%-24s %.1e        %.1e          %.1e' % (c3, r['v2'][0], r['v2'][1], r['k64'][1]))
compat.set_mode(None)
