# -*- coding: utf-8 -*-
r"""02 章 (ドリフト) の途中式を sympy で確かめる (本文に書き足した展開の確認)。

(a) §7.1: 運動学的ドリフト D の流れ (ΔH の偏微分) が
    Δx = q (1/ζ - 1) L、Δz = -(2+ζ)(1-ζ)^2/(2ζ) L (ζ = sqrt(1 - q^2 - r^2)) になること。
(b) §7.1: 小さい傾きで 1/ζ - 1 = (q^2 + r^2)/2 + (3/8)(q^2+r^2)^2 + ... (Δx ≃ q(q^2+r^2) L/2)。
(c) §7.1: 近軸 + D = 厳密 (z の和)。
(d) §9: M^T J M = J の (j,6) 成分から (2.15) を作る。T_66 = 1、T_6j = 0 を使うこと。
    ドリフトの 6×6 行列 (2.13) が (2.15) を満たすこと、および T_66 = 1, T_6j = 0 でない一般の
    シンプレクティック行列では (2.15) が成り立たないこと (簡単な例)。
"""
import sympy as sp
import _common as C

C.start(__file__)

q, r, p = sp.symbols('q r p', positive=True)
px, py = sp.symbols('px py', positive=True)
pz = sp.sqrt(p**2 - px**2 - py**2)

# (a) ------------------------------------------------------------------
print('(a) D の流れ')
dH = p - pz - (px**2 + py**2) / (2 * p)
dx = sp.diff(dH, px)
dz = sp.diff(dH, p)
zeta = pz / p
sub = {px: q * p, py: r * p}
zz = sp.sqrt(1 - q**2 - r**2)
print('  dΔH/dpx = q (1/ζ - 1): ', sp.simplify(dx.subs(sub) - q * (1 / zz - 1)) == 0)
print('  dΔH/dp  = -(2+ζ)(1-ζ)^2/(2ζ): ',
      sp.simplify(dz.subs(sub) + (2 + zz) * (1 - zz)**2 / (2 * zz)) == 0)
Z = sp.symbols('Z', positive=True)
print('  1 - 1/Z + (1 - Z^2)/2 の通分: ', sp.simplify(1 - 1 / Z + (1 - Z**2) / 2 + (2 + Z) * (1 - Z)**2 / (2 * Z)) == 0)
print('  Z^3 - 3Z + 2 = (Z-1)^2 (Z+2): ', sp.expand((Z - 1)**2 * (Z + 2) - (Z**3 - 3 * Z + 2)) == 0)

# (b) ------------------------------------------------------------------
print('(b) 小さい傾きの展開')
u = sp.symbols('u')   # u = q^2 + r^2
print('  1/sqrt(1-u) - 1 =', sp.series(1 / sp.sqrt(1 - u) - 1, u, 0, 3))
print('  -(2+ζ)(1-ζ)^2/(2ζ) =', sp.series(-(2 + sp.sqrt(1 - u)) * (1 - sp.sqrt(1 - u))**2 / (2 * sp.sqrt(1 - u)), u, 0, 4))

# (c) ------------------------------------------------------------------
print('(c) 近軸 + D = 厳密')
paraz = -(q**2 + r**2) / 2
Dz = -(2 + zz) * (1 - zz)**2 / (2 * zz)
print('  -(q^2+r^2)/2 + Dz = -(1-ζ)/ζ: ', sp.simplify(paraz + Dz + (1 - zz) / zz) == 0)

# (d) ------------------------------------------------------------------
print('(d) (2.15) の導出と成立条件')
M = sp.Matrix(6, 6, lambda i, j: sp.Symbol('M%d%d' % (i + 1, j + 1)))
J = sp.zeros(6, 6)
for k in range(3):
    J[2 * k, 2 * k + 1] = 1
    J[2 * k + 1, 2 * k] = -1
cond = (M.T * J * M)
j = sp.symbols('j')
# (j,6) 成分 (j = 1..4): sum_ab M_aj J_ab M_b6 = 0
for jj in range(1, 5):
    expr = cond[jj - 1, 5]
    print('  (%d,6) 成分 =' % jj, expr)
# ドリフトの行列
L, pr_ = sp.symbols('L pr', positive=True)
a = px**2 + py**2
pzi = sp.sqrt(pr_**2 - a)
ale = L / pzi
alz = ale / pzi**2
dvdp = sp.Symbol('dvdp')
T = sp.eye(6)
T[0, 1] = ale + px**2 * alz; T[0, 3] = px * py * alz; T[0, 5] = -px * pr_ * alz
T[2, 1] = T[0, 3]; T[2, 3] = ale + py**2 * alz; T[2, 5] = -py * pr_ * alz
T[4, 1] = T[0, 5]; T[4, 3] = T[2, 5]; T[4, 5] = dvdp * L + a * alz
print('  ドリフトの T は M^T J M = J:', sp.simplify(T.T * J * T - J) == sp.zeros(6, 6))
ok = True
for jj in range(1, 5):
    rhs = -T[0, jj - 1] * T[1, 5] + T[1, jj - 1] * T[0, 5] - T[2, jj - 1] * T[3, 5] + T[3, jj - 1] * T[2, 5]
    ok &= sp.simplify(T[4, jj - 1] - rhs) == 0
print('  (2.15) が T_5j (j=1..4) を再現:', ok)
# 反例: 一般のシンプレクティック行列 (T_66 != 1, T_6j != 0) では (2.15) が成り立たない (数値)
import numpy as np
from scipy.linalg import expm
rng = np.random.default_rng(1)
Jn = np.array(J.tolist(), dtype=float)
S = rng.normal(size=(6, 6)); S = 0.2 * (S + S.T)
Mn = expm(Jn @ S)
print('  ランダムなシンプレクティック行列: |M^T J M - J| = %.1e、T_66 = %.3f' % (abs(Mn.T @ Jn @ Mn - Jn).max(), Mn[5, 5]))
for jj in range(4):
    rhs = -Mn[0, jj] * Mn[1, 5] + Mn[1, jj] * Mn[0, 5] - Mn[2, jj] * Mn[3, 5] + Mn[3, jj] * Mn[2, 5]
    # 一般の式: M5j M66 - M6j M56 = rhs  (T_66 = 1, T_6j = 0 のとき M5j = rhs)
    gen = Mn[4, jj] * Mn[5, 5] - Mn[5, jj] * Mn[4, 5]
    print('  j=%d: (2.15) の右辺 %.4f、T_5j %.4f (一致せず)、一般形 T_5j T_66 - T_6j T_56 = %.4f (右辺と一致)' % (jj + 1, rhs, Mn[4, jj], gen))
