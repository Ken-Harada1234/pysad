# -*- coding: utf-8 -*-
"""(B108 の残り) 端の角 AE1・AE2 のある極 (K0 を持つ ANGLE=0 の要素) の ΔJx を、厳密な幾何 (円弧・微分方程式) + ロビンソンの定理と比べる。
24 セルの FODO (3 GeV、周長 100.8 m、QF D2 B D2 QD DS [P1 DW P2 DW P1] DS D1。極 3 個の組は K0・-2K0・K0 で、軌道は極の外で 0 に戻る)。AE なしとの差 dJx。
独立な値: 硬い端の磁場の極の厳密な運動 (面は実験室系で z = x tan a。K1 つきは微分方程式 B(x) = (K0 + K1 x)/L を数値積分) と Je = 2 + dln(sum ∫B^2 dl)/d(delta)、Jx = 3 - Je。
比べるもの: pysad の Emittance[] (EMITDIV=20) と、追跡 (RF・放射つきの一周写像の固有値から減衰率) を、ANGLE=0 の BEND (tsteee・tsteer) と MULT (tmulte・tmulti) の両方で。
使い方: py -3.12 example/01_lattice_elements/elem_AE_exact_geometry.py (数分)。結果は elem_AE_exact_geometry.out.txt。
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pysad import compat
from pysad.tests.script_tests import _p3_indep as X
from pysad.tests.script_tests import test_v2_E_steer as T
from pysad.tests.script_tests import test_v2_P3_ae as A

K0 = 0.005
print('極 = MULT。K0 = %g rad/極。比 = pysad の dJx / 厳密な dJx (v2 は 1.03〜1.04、k64 は約 2 倍)' % K0)
for K1 in (0.0, 0.06):
    e00 = X.jx_exact_k1(K0, K1, 0.0, 0.0)
    for ae in ((0.3, 0.3), (-0.2, 0.4), (0.3, -0.3)):
        ex = X.jx_exact_k1(K0, K1, *ae) - e00
        res = {}
        for mode in ('v2', 'k64'):
            compat.set_mode(mode)
            a, ta = A._emit_jx_and_tracking_j(A._pole_ring_mult(K0, ae, 'MULT', K1))
            b, tb = A._emit_jx_and_tracking_j(A._pole_ring_mult(K0, (0.0, 0.0), 'MULT', K1))
            res[mode] = (a - b, ta - tb)
        compat.set_mode('k64')
        print('K1 %.2f AE %s: 厳密 dJx %+.5f | v2 Emittance[] %+.5f (%.3f) 追跡 %+.5f (%.3f) | k64 Emittance[] %+.5f (%.3f) 追跡 %+.5f (%.3f)' %
              ((K1, ae, ex, res['v2'][0], res['v2'][0] / ex, res['v2'][1], res['v2'][1] / ex, res['k64'][0], res['k64'][0] / ex, res['k64'][1], res['k64'][1] / ex)), flush=True)
print('極 = ANGLE=0 の BEND (K1 = 0)。K64 の追跡は出口の AE2 を FB2 != 0 のときだけ入れるので k64 の追跡の比は不定')
e00 = T._jx_exact(K0, 0.0, 0.0)
for ae in ((0.3, 0.3), (-0.2, 0.4), (0.3, -0.3)):
    ex = T._jx_exact(K0, *ae) - e00
    res = {}
    for mode in ('v2', 'k64'):
        compat.set_mode(mode)
        a, ta = A._emit_jx_and_tracking_j(T._pole_ring(K0, ae, 'EMITDIV=20;'))
        b, tb = A._emit_jx_and_tracking_j(T._pole_ring(K0, (0.0, 0.0), 'EMITDIV=20;'))
        res[mode] = (a - b, ta - tb)
    compat.set_mode('k64')
    print('AE %s: 厳密 dJx %+.5f | v2 Emittance[] %+.5f (%.3f) 追跡 %+.5f (%.3f) | k64 Emittance[] %+.5f (%.3f) 追跡 %+.5f (%.3f)' %
          ((ae, ex, res['v2'][0], res['v2'][0] / ex, res['v2'][1], res['v2'][1] / ex, res['k64'][0], res['k64'][0] / ex, res['k64'][1], res['k64'][1] / ex)), flush=True)
print('残る 2.7% の切り分け (v2 の Emittance[]、MULT でなく BEND の極、K0 = 0.005): AE の大きさ・EMITDIV・K0 への依存')
compat.set_mode('v2')
j00 = A._emit_jx_and_tracking_j(T._pole_ring(K0, (0.0, 0.0), 'EMITDIV=20;'))[0]
for ae in ((0.05, 0.05), (0.1, 0.1), (0.3, 0.3), (0.5, 0.5)):
    ex = T._jx_exact(K0, *ae) - e00
    d = A._emit_jx_and_tracking_j(T._pole_ring(K0, ae, 'EMITDIV=20;'))[0] - j00
    print('  AE %s: 比 %.4f' % (ae, d / ex))
for div in (1, 5, 80, 200):
    ex = T._jx_exact(K0, 0.3, 0.3) - e00
    a = A._emit_jx_and_tracking_j(T._pole_ring(K0, (0.3, 0.3), 'EMITDIV=%d;' % div))[0]
    b = A._emit_jx_and_tracking_j(T._pole_ring(K0, (0.0, 0.0), 'EMITDIV=%d;' % div))[0]
    print('  EMITDIV=%d: 比 %.4f' % (div, (a - b) / ex))
for k0 in (0.0025, 0.03):
    ex0 = T._jx_exact(k0, 0.0, 0.0)
    ex = T._jx_exact(k0, 0.1, 0.1) - ex0
    a = A._emit_jx_and_tracking_j(T._pole_ring(k0, (0.1, 0.1), 'EMITDIV=20;'))[0]
    b = A._emit_jx_and_tracking_j(T._pole_ring(k0, (0.0, 0.0), 'EMITDIV=20;'))[0]
    print('  K0=%g AE=(0.1, 0.1): 比 %.4f' % (k0, (a - b) / ex))
