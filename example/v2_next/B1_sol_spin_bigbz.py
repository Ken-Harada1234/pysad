# -*- coding: utf-8 -*-
"""v2next B1 項目 6: SOL 領域内の POL (スピン)・加速 + RAD・bz ≥ 3 1/m。
(a) 領域内の DRIFT (0.4 m) のスピン (POL + RAD の追跡、軸を外れた粒子) と Thomas-BMT の数値積分の差。bz = 3・4.5・9 1/m。v2 は DRIFT を bz·al/ndiv ≤ 0.02 に刻む。
(b) POL ON の Emittance[] の SpinTune: 平面リングの SOL 領域 (DRIFT・MULT・CAVI) と解析式 cos(μ/2) = cos(π ν0) cos(φ/2)。
(c) 領域内の RF つき MULT の加速 + RAD: 放射損失 (RAD ON − OFF) と (2/3) r_e mc² (bz γ0² θ0)² L / E2 の比 (加速 0〜3.1%)。
使い方: py -3.12 example/v2_next/B1_sol_spin_bigbz.py"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from pysad import compat
from pysad.tests import test_v2_B1_sol as T

print('(a) 領域内の DRIFT のスピン (0.4 m、δ = 1e-3・−5e-3 の粒子 3 つ × 3 つのスピン): Thomas-BMT との最大の差')
for mode in ('v2', 'k64'):
    compat.set_mode(mode)
    for p, bzT in ((0.1, 1.5), (0.1, 3.0), (0.3, 3.0)):
        ffs, _ = T._run(T._BASE % (p, 'DRIFT E1 = (L=0.4);', bzT))
        bz = bzT * T.CL / (p * 1e9)
        fld = lambda x_, y_, s_: np.array([0.0, 0.0, bz])
        worst = 0.0
        for s0 in T._S0:
            b, s = T._track(ffs, T._P, 1, len(ffs.els) + 1, [s0] * len(T._P))
            for c, sv in zip(T._P, s):
                ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
                u = T._tbmt(ck, s0, 0.4, fld, p * 1e9 / T.ME)
                worst = max(worst, float(np.abs(sv - u[6:]).max()))
        print('  %-3s p = %.1f GeV  bz = %.2f 1/m  最大の差 %.1e' % (mode, p, bz, worst))

compat.set_mode('v2')
print('(b) POL の Emittance[] の SpinTune: ソレノイドの領域 (L = 1 m) の解析式との差 (bz = 4.5・9 1/m は 0.1 GeV)')
import math
for name, region, extra in (('drift', 'DS DS', ''), ('MULT', 'DS MM DS', 'MULT MM=(L=0.5 K1=0.0 DISFRIN=1);'), ('CAVI VOLT=0', 'DS CC DS', 'CAVI CC=(L=0.5 VOLT=0 HARM=0 FREQ=5e8);')):
    for p, bzT in ((0.5, 1.0), (0.1, 1.5), (0.1, 3.0)):
        L = 0.5 if name == 'drift' else 1.0
        nu0 = T._spin_tune(region, extra, 0.0, p)
        nu = T._spin_tune(region, extra, bzT, p)
        bz = bzT * T.CL / (p * 1e9)
        half = math.cos(math.pi * nu0) * math.cos((1 + T.GS) * bz * L / 2)
        ana = math.acos(half) / math.pi
        print('  %-12s p = %.1f GeV  bz = %.2f 1/m  SpinTune %.12f  解析 %.12f  差 %.1e' % (name, p, bz, abs(nu), min(ana, 1 - ana) if abs(nu) < 0.5 else ana,
                                                                                     min(abs(abs(nu) - ana), abs(abs(nu) - (1 - ana)))))

print('(c) 領域内の RF つき MULT (3 GeV、θ0 = 1e-2、L = 0.5 m) の放射損失 (RAD ON − OFF) / 式')
for bz in (3.0, 9.0):
    for volt in (0.0, 1e6, 2e7, 1e8):
        P = 3.0e9
        el = T._elem('MULT E1 = (L=0.5 VOLT=%g FREQ=1.3E9 PHI=-1.2 DISFRIN=1);' % volt, '3GEV')

        def run(rad):
            ph = T._e._Phys(P, 100.0)
            col = T._Col()
            T.x._compile_sol_v2(col, el, 'Multipole', bz, ph, 2)
            k, pr = col.items[0]
            S = T.tn.make_state(np.array([[0, 1e-2, 0, 0, 0, 0]], float), ph)
            ph2 = T.x.step_kind(k, pr, S, 1, ph, 0.0, True, rad, False, T.tn.make_rs(1), tm=10.0)
            return T.tn.state_out(S, 1)[0], (ph2 or ph)
        on, ph2 = run(True)
        off, _ = run(False)
        g0 = math.sqrt(1 + (P / T.ME) ** 2)
        exp = -(2.0 / 3.0) * T.RE * T.ME * (bz * g0**2 * 1e-2)**2 * 0.5 / (ph2.p0 * T.ME)
        print('  bz = %.0f 1/m  VOLT = %.0e V (p2/p0 = %.4f)  比 %.5f' % (bz, volt, ph2.p0 / (P / T.ME), (on[5] - off[5]) / exp))
compat.set_mode(None)
