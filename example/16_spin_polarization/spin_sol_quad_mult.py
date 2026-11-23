# SOL 領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD) の Thomas-BMT の積分との差 (最大の成分)。
# 使い方: py -3.12 example/16_spin_polarization/spin_sol_quad_mult.py [v2|k64]。条件: 0.4 m の要素、DISFRIN=1、軸上と軸を外れた粒子 3 つ、3 つのスピンの向き。
# v2 はどの bz でも 1e-6 以下。k64: QUAD は K64 の bsi の式 (1 スライス分足りない) で 6e-2〜0.19、MULT は領域内が未対応 (NotImplementedError)。
import math
import sys

import numpy as np

import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
from pysad import compat
compat.set_mode(sys.argv[1] if len(sys.argv) > 1 else 'v2')
from pysad.tests.script_tests.test_v2_B1_sol import CL, ME, _BASE, _P, _S0, _run, _tbmt, _track

CASES = {
    'quad': ('QUAD E1 = (L=0.4 K1=0.3 DISFRIN=1);', lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.75 * x, bz]))),
    'mult': ('MULT E1 = (L=0.4 K0=0.02 K1=0.3 DISFRIN=1);', lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.05 + 0.75 * x, bz]))),
}
print('mode', compat._mode)
for case, (el, fld) in CASES.items():
    for p, bzT in [(1.0, 1.0), (0.3, 3.0), (0.1, 1.5), (0.1, 3.0)]:
        bz = bzT * CL / (p * 1e9)
        try:
            ffs, _ = _run(_BASE % (p, el, bzT))
            worst = 0.0
            for s0 in _S0:
                b, s = _track(ffs, _P, 1, len(ffs.els) + 1, [s0] * len(_P))
                for c, sv in zip(_P, s):
                    ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
                    u = _tbmt(ck, s0, 0.4, fld(bz), p * 1e9 / ME)
                    worst = max(worst, float(np.abs(sv - u[6:]).max()))
            print('%-5s p = %.1f GeV  bz = %.2f 1/m  spin error vs Thomas-BMT %.2e' % (case, p, bz, worst))
        except NotImplementedError as e:
            print('%-5s p = %.1f GeV  bz = %.2f 1/m  NotImplementedError (k64: region MULT unsupported)' % (case, p, bz))
