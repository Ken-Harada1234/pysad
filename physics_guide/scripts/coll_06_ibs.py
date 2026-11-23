# -*- coding: utf-8 -*-
r"""IBS (INTRA) の導出の確認 (17_collective.md §12)。

K64 の tintrb の拡散行列の核 e_i = f_j + f_k - 2 f_i (f_i = 2 λ_i ∫_0^1 fintrbu(u; λ_i/λ_j, λ_i/λ_k) du) を、
Rutherford 散乱の小角の極限から導いた  d<p_i p_j>/dt = 4π r_e² c n lnΛ <(q² δ_ij - 3 q_i q_j)/|q|³>  (q = p1 - p2 ~ N(0, 2P)) とモンテカルロで比べる。
τ_ii ≡ (4π)^{3/2} √det(P) <(q² - 3 q_i²)/|q|³> (P の主軸で) に対し、K64 の e_i が τ_ii/(4π) に等しいはず。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
from pysad import touschek as T

rng = np.random.default_rng(1)
nsamp = 6_000_000
print('λ = P の固有値 (p0 単位の運動量の分散)。e_i: K64 の式 (bint = 台形則 + 外挿)、τ_ii/(4π): モンテカルロ (%.0e 個)' % nsamp)
print(' λ1, λ2, λ3                 i   e_i (K64)     τ_ii/(4π) (MC)   比')
for lam in ((1.0, 0.3, 0.1), (1.0, 1.0, 0.2), (4.0, 1.0, 0.05), (1.0, 1e-2, 1e-3)):
    lam = np.array(lam)
    f = []
    for i in range(3):
        j, k = [m for m in range(3) if m != i]
        a1, a2 = lam[i] / lam[j], lam[i] / lam[k]
        f.append(2.0 * T.bint(lambda u: T.fintrbu(u, a1, a2), 0.0, 1.0, 1e-3, 3e-17) * lam[i])
    e = [f[1] + f[2] - 2 * f[0], f[2] + f[0] - 2 * f[1], f[0] + f[1] - 2 * f[2]]
    q = rng.normal(size=(nsamp, 3)) * np.sqrt(2 * lam)
    q2 = (q**2).sum(axis=1)
    qa = np.sqrt(q2)
    pref = (4 * math.pi)**1.5 * math.sqrt(lam.prod())
    print('  Σ e_i (K64) = %.2e (トレースレス: 静止系の運動量空間のエネルギーは保存)' % sum(e))
    for i in range(3):
        v = (q2 - 3 * q[:, i]**2) / qa**3
        tau = pref * v.mean() / (4 * math.pi)
        err = pref * v.std() / math.sqrt(nsamp) / (4 * math.pi)
        print(' %-24s %d  %+.5e  %+.5e ± %.1e  %.4f' % (str(tuple(float(x) for x in lam)), i + 1, e[i], tau, err, e[i] / tau))
