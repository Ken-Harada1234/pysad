# -*- coding: utf-8 -*-
"""補正磁石 (ANGLE=0 の BEND、K64 の tsteer) の F1 (端の磁場が直線で立ち上がる長さ) の追跡を、独立な基準 (x 方向に一様な Maxwell 場の端の厳密な運動) と比べる。
磁場 By = B(s) - y^2 B''/2、Bs = y B'(s) (div B = 0、rot B = 0) のベクトルポテンシャル Ax = F(s) - y^2 B'/2 (F' = B) の正準方程式を solve_ivp (DOP853、rtol 1e-13) で積分する
(pysad/tests/script_tests/_p3_indep.maxwell_steer)。B(s) は入口 [-fb1/2, fb1/2]・出口 [al-fb2/2, al+fb2/2] で 0 -> B0 -> 0 に直線で変わる (fb = F1 + FB1 など。B0 = K0/al)。
比べるもの (3 GeV、L = 0.2 m、K0 = 0.003 rad の補正磁石、FRINGE = 3):
 (1) 鉛直の kick の差 dpy/y (F1 あり - なし) = (fb1 + fb2)/(6 rho^2) (rho = al/K0)。厳密な運動の py/y は -B0^2 al + B0^2 (fb1+fb2)/6 (解析式 y ∫F B' ds)。
 (2) 出口の x・px (fb1 = fb2 のとき、δ ≠ 0 の色収差も含めて 1e-10 で一致)。fb1 ≠ fb2 では x が (K0/al)(fb1^2 - fb2^2)/24 だけ違う (K64 の端の模型の仕様)。
 (3) 放射の積分 ∫B^2 dl の比 = 1 - (fb1+fb2)/(6 al)。
使い方: py -3.12 example/01_lattice_elements/elem_BEND_steer_F1_maxwell_ramp.py (数秒)。結果は elem_BEND_steer_F1_maxwell_ramp.out.txt。v2・k64 とも同じ値 (F1 の追跡は K64 の動作の追加)。
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
import numpy as np
from pysad import compat
from pysad.tests.script_tests import _p3_indep as X

compat.set_mode('v2')
al, K0, y0 = 0.2, 0.003, 1e-5
B0 = K0 / al
print('(1) 鉛直の kick の差 dpy/y [1/m]: (F1 あり - なし)。K64 (pysad 追跡)、解析式 (fb1+fb2)/(6 rho^2)、厳密な運動 (py/y の fb 依存部分)')
for fb1, fb2 in ((0.04, 0.04), (0.06, 0.02), (0.1, 0.0001)):
    ex = X.maxwell_steer(fb1, fb2, al, K0, y0)[3] / y0
    c = np.zeros(6)
    c[2] = y0
    k = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r F1=%r FB1=0 FB2=%r FRINGE=3);' % (al, K0, fb1, fb2 - fb1), c)
    k0 = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r);' % (al, K0), c)
    print('  fb1 %.4f fb2 %.4f: K64 %.6e  解析式 %.6e  厳密な運動 py/y = %.6e (= -B0^2 al + (fb1+fb2)/6 B0^2 = %.6e)' %
          (fb1, fb2, (k[3] - k0[3]) / y0, (fb1 + fb2) / 6 * B0**2, ex, B0**2 * (-al + (fb1 + fb2) / 6)))
print('(2) 出口の x (出口の面 s = al、K64 の追跡 と 厳密な運動、δ 依存)')
for fb1, fb2, d in ((0.04, 0.04, 0.0), (0.04, 0.04, 0.01), (0.04, 0.04, -0.02), (0.1, 0.06, 0.02), (0.3, 0.0001, 0.01)):
    x, px, _y, _py = X.maxwell_steer(fb1, fb2, al, K0, 0.0, d)
    xnom = x - (fb2 / 2) * px / math.sqrt((1 + d)**2 - px * px)
    c = np.zeros(6)
    c[5] = d
    k = X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=%r F1=%r FB1=0 FB2=%r FRINGE=3);' % (al, K0, fb1, fb2 - fb1), c)
    print('  fb1 %.4f fb2 %.4f delta %+.2f: x K64 %.10e  厳密 %.10e  差 %+.2e  (K0/al (fb1^2-fb2^2)/24 = %.2e)  px K64 %.10f 厳密 %.10f' %
          (fb1, fb2, d, k[0], xnom, k[0] - xnom, B0 * (fb1**2 - fb2**2) / 24, k[1], px))
print('(3) 放射の積分: エネルギー損失の比 (F1 あり / なし)。高速追跡と Python の追跡')
for fb1, fb2 in ((0.04, 0.04), (0.1, 0.06), (0.05, 0.0001)):
    r = []
    for fast in (True, False):
        u0 = -X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=0.02);' % al, np.zeros(6), rad=True, fast=fast)[5]
        u = -X.one_element_track('BEND SA = (L=%r ANGLE=0 K0=0.02 F1=%r FB2=%r FRINGE=3);' % (al, fb1, fb2 - fb1), np.zeros(6), rad=True, fast=fast)[5]
        r.append(u / u0)
    print('  fb1 %.4f fb2 %.4f: 高速 %.7f Python %.7f  解析式 1 - (fb1+fb2)/(6 al) = %.7f' % (fb1, fb2, r[0], r[1], 1 - (fb1 + fb2) / (6 * al)))
