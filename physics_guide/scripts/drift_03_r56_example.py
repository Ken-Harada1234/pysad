# -*- coding: utf-8 -*-
"""最短の例: ドリフトだけの運動量依存の遅れ (02_drift.md §9)。R56 に当たる量 ∂z/∂δ = L h0/h1^3 ≈ L/γ0^2。

(a) 3 GeV の電子 (γ0 = 5871) と 10 MeV の電子 (γ0 = 20.6) で、δ ごとに、長さ 1 m のドリフトを k = 1〜4 回通った後の z (pysad の追跡) と、
    k·(-dv L) (01 章の式 (1.9)) の比較。z は回数に比例して遅れる/進む。
(b) 同じことを 1 周 (4 つのドリフトのラティス) の行列で: CALC の行列の T(5,6) の積。
(c) z = -v (t - t0) との対応: 速い粒子 (δ > 0) は z > 0 (先頭)。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
from pysad import lattice as lat

MC2 = 0.51099895e6
Jm = None


def dv_of(dp, pev):
    p0 = pev / MC2
    h0 = math.sqrt(1 + p0**2)
    h1 = math.sqrt(1 + (p0 * (1 + dp))**2)
    return -dp * (2 + dp) / (h1 * (h1 + (1 + dp) * h0)), p0, h0, h1


for label, mom in (('3GEV', 3.0e9), ('10MEV', 10.0e6)):
    p0 = mom / MC2
    g0 = math.sqrt(1 + p0**2)
    print('== (a) %s の電子: p0 = %.3f、γ0 = %.3f。長さ 1 m のドリフトを k 回通った後の z (δ ごと。x = px = 0) ==' % (label, p0, g0))
    print('   δ        dv (1.9)          k=1: z (追跡)    k=2: z (追跡)    k=3: z (追跡)    k=4: z (追跡)    k=4 の z と -4 dv L の差   L δ/γ0^2 (近似)')
    ffs4, _ = C.make_drift(1.0, '', label, n=4)
    ffs1, _ = C.make_drift(1.0, '', label, n=1)
    for dp in (1e-3, 1e-2, -1e-2):
        dv, _, _, _ = dv_of(dp, mom)
        zs = []
        for k in (1, 2, 3):
            ffsk, _ = C.make_drift(1.0, '', label, n=k)
            zs.append(C.track_one(ffsk, np.array([0, 0, 0, 0, 0, dp]))[4])
        z4 = C.track_one(ffs4, np.array([0, 0, 0, 0, 0, dp]))[4]
        zs.append(z4)
        print('  %+.0e  %+.6e  %+.6e  %+.6e  %+.6e  %+.6e   %.1e            %+.3e' % (dp, dv, zs[0], zs[1], zs[2], zs[3], abs(z4 + 4 * dv), dp / g0**2))
    print()

print('== (b) 4 つのドリフトの行列 (CALC 側。3 GeV、COD = 設計軌道) ==')
mom = 3.0e9
p0 = mom / MC2
h0 = math.sqrt(1 + p0**2)
T = np.eye(6)
for _ in range(4):
    T = lat._tdrife1_mat(np.zeros(6), 0.0, h0 / h0**3, 1.0) @ T
print('  T(1,2) = %.6f (= 4 L)、T(5,6) = %.6e (= 4 L/γ0^2 = %.6e)' % (T[0, 1], T[4, 5], 4 / h0**2))
print('  δ = 1e-3 の粒子の z の変化 (行列) = T(5,6) δ = %.4e。追跡 (a) の k = 4 の z = %.4e (差は δ の 2 次の項)' %
      (T[4, 5] * 1e-3, C.track_one(C.make_drift(1.0, '', '3GEV', n=4)[0], np.array([0, 0, 0, 0, 0, 1e-3]))[4]))

print('\n== (c) 符号 ==')
print('  δ > 0 (速い) では dv < 0、z\' = -dv > 0: z が増える = 先頭 (z = -v (t - t0) で、早く着いた粒子が z > 0)')
print('  3 GeV: δ = +1e-3 の粒子は 4 m で z = %+.3e m' % (4 * -dv_of(1e-3, 3.0e9)[0]))
print('  10 MeV: 同じ δ で 4 m で z = %+.3e m (γ0 が小さいほど R56 が大きい)' % (4 * -dv_of(1e-3, 10.0e6)[0]))
