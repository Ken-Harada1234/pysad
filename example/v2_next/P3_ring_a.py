# -*- coding: utf-8 -*-
"""v2next P3: 大きなリング A (pysad/tests/kek-ls_v03_12.lat、3 GeV、周長 603 m、3563 要素) の先頭に補正磁石 2 台 (K0 = +-5e-4 rad、L = 0.2 m) を入れて確かめる。
  ra1: BEND ANGLE=0 に F1=0.04 FB1=0.03 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1
  ra2: MULT (K1 = +-0.05) に F1=0.03 F2=0.01 FB1=0.02 FB2=0.02 FRINGE=3 AE1=0.2 AE2=-0.1
  ra3: BEND ANGLE=0 に AE1・AE2 だけ (FB2 なし。K64 の追跡は AE2 を入れない)
(1) k64 モード: 追跡 (NORAD・RAD、3 粒子 3 周) と CALC・Emittance[] の gs (K64 の 10-19 ビルド) との差 (gs の出力は P3_ring_a_ra*.gs.txt)。
(2) 追跡の一周写像 (4 次元、閉軌道のまわりの数値微分) の固有値から求めたチューンの小数部と CALC のチューンの差 (k64 と v2)。
使い方: py -3.12 example/v2_next/P3_ring_a.py (数十秒)。結果は P3_ring_a.out.txt。
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import numpy as np
from pysad import compat
from pysad.tests import _p3_gsdata as G
from pysad.tests import _p3_indep as X
from pysad.tests import test_v2_E_steer as T

CALC_PRINT = 'FFS;\nUSE RING; CELL; CALC;\nPrint["CALC tunes ", Twiss["NX","$$$"]/(2Pi), " ", Twiss["NY","$$$"]/(2Pi)];\nExit;\n'
print('(1) k64 モード: gs との差 (最大絶対差。追跡は 10 桁表示の全桁一致 = 0)')
compat.set_mode('k64')
for k in ('ra1', 'ra2', 'ra3'):
    out = X.run(X.ring_a_script(k, X.RING_A_TAIL))[1]
    n = np.array(X.tag_values(out, 'N3'))
    r = np.array(X.tag_values(out, 'R3'))
    d = T._values(out)
    ref = G.EMIT[k]
    print('  %s: NORAD %.1e  RAD %.1e | CALC tunes %.1e Emittance[] tunes %.1e U0 相対 %.1e 減衰率 相対 %.1e' %
          (k, np.abs(n - np.array(G.TRACK[k][0])).max(), np.abs(r - np.array(G.TRACK[k][1])).max(),
           np.abs(np.array(d['calc']) - ref['calc']).max(), np.abs(np.array(d['tunes']) - ref['tunes']).max(),
           abs(d['u0'] / ref['u0'] - 1), np.abs(np.array(d['damp']) / np.array(ref['damp']) - 1).max()))
print('(2) 追跡の一周写像のチューン (小数部、昇順) と CALC のチューンの差 (追跡 - CALC、x・y の順に小さい方から)')
for k in ('ra1', 'ra2', 'ra3'):
    for mode in ('k64', 'v2'):
        compat.set_mode(mode)
        nu, calc = X.tracking_and_calc_tunes(X.ring_a_script(k, CALC_PRINT))
        print('  %s %s: 追跡 %s CALC %s 差 %s' % (k, mode, ['%.8f' % v for v in nu], ['%.8f' % v for v in calc], ['%+.1e' % v for v in (np.array(nu) - np.array(calc))]))
