# -*- coding: utf-8 -*-
"""v2next E 項目 2 (B108): ANGLE=0 の BEND (補正磁石・挿入光源の極) の端の角 AE1・AE2 と放射減衰の分配 Jx。
K64 は端の蹴りを偏向角に入れるが場の中の行路を変えないので、端の項 (-eta tan(e)/rho^2) が約 2 倍になる (ANGLE≠0 の BEND の B75 と同じ構造)。v2 は行路を x tan(e) だけ短くする。
独立な値 = 硬い端の一様磁場の極の厳密な幾何 (円弧と直線) とロビンソンの定理 Je = 2 + dln(∮B² dl)/dδ、Jx = 3 − Je (pysad の光学・放射のコードを使わない)。
リング: 24 セルの FODO (3 GeV、QF D2 B D2 QD DS [P1 DW P2 DW P1] DS D1)。極 3 個の組 (K0・−2K0・K0。軌道は極の外で 0 に戻る)。
表は「AE なし」との差 dJx (AE の効果)。使い方: py -3.12 example/v2_next/E_b108_poles.py (約 2 分)。結果は E_b108_poles.out.txt。"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'pysad', 'tests'))
from pysad import compat
import test_v2_E_steer as T

print('K0 [rad/極]  (AE1, AE2) [rad]   dJx 厳密な幾何    dJx v2     (v2/厳密)    dJx k64    (k64/厳密)')
for K0 in (0.0025, 0.005):
    base = {}
    for mode in ('v2', 'k64'):
        compat.set_mode(mode)
        base[mode] = T._emit_jx(K0, (0.0, 0.0))
    ex0 = T._jx_exact(K0, 0.0, 0.0)
    for ae in ((0.3, 0.3), (-0.2, 0.4), (0.3, -0.3)):
        ex = T._jx_exact(K0, *ae) - ex0
        d = {}
        for mode in ('v2', 'k64'):
            compat.set_mode(mode)
            d[mode] = T._emit_jx(K0, ae) - base[mode]
        print('%6.4f   %s   %+.5f   %+.5f   (%.3f)   %+.5f   (%.3f)' % (K0, ae, ex, d['v2'], d['v2'] / ex, d['k64'], d['k64'] / ex))
