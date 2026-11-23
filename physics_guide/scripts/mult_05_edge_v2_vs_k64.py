# -*- coding: utf-8 -*-
r"""MULT の端面の角 (AE1) の kick: K64 の tmulk (08_mult.md の式 (8.9)) と、pysad v2 の厳密な面の幾何 (_mult_edge_exact + 縦の蹴り) の比較 (§4.4)。

短い MULT (L = 0.02 m、積分強度 K1 = 0.3) の入口にだけ角 AE1 = 0.05 を付け、(x, y) = (a, a) の粒子の通過後の (px, py) を、
k64 (tmulk) と v2 で比べる。差が edge の kick の定義の違い。式 (8.9) の K1: Δpx = tK1(x²-y²)/L、Δpy = -2tK1xy/L (t = tan AE1)。
v2 は Δpx = tK1x²/L (面までの余分な経路) で、y² の項が無い。Δpy は合計が一致する (面の幾何の分 + Maxwell の分)。
"""
import io
import numpy as np
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6, compat


def run(mode, extra, c, L=0.02):
    compat.set_mode(mode)
    src = 'MOMENTUM = 3GEV;\nMULT MP=(L=%g %s);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n' % (L, extra)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6.track_particles(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


L = 0.02
tt = np.tan(0.05)
print('== L = %g m、AE1 = 0.05 (t = tan = %.6f)。初期値 (x, y) = (a, a)、px = py = 0。差 = v2 - k64 ==' % (L, tt))
for extra, lab, K in (('K0=0.01 AE1=0.05', 'K0 = 0.01', 0.01), ('K1=0.3 AE1=0.05', 'K1 = 0.3', 0.3)):
    for amp in (1e-3, 2e-3):
        c = [amp, 0, amp, 0, 0, 0]
        a = run('v2', extra, c)
        b = run('k64', extra, c)
        print('  %-9s a = %g: px (v2, k64, 差) = %.4e %.4e %.3e | py = %.4e %.4e %.3e' % (lab, amp, a[1], b[1], a[1] - b[1], a[3], b[3], a[3] - b[3]))
print('  比較: tK1/L × a² = %.4e (K1 の px の差は +tK1 y²/L に一致する。a² に比例)。K0 は差が 1e-11 以下 (K0 では v2 = K64 の kick)。' % (tt * 0.3 / L * 1e-6))
compat.set_mode('v2')
