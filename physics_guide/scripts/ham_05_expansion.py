# -*- coding: utf-8 -*-
"""P0・δ による規格化と 2 次の展開、偏向磁石・四極の運動方程式 (01a_hamiltonian_derivation.md §9〜§10) を sympy と数値積分で確かめる。

H~/P0 = -a_s - (1 + x/ρ) sqrt((1 + δ)^2 - p_x^2 - p_y^2)、a_s = e A_s/P0。
(a) 横方向の量 (x, y, p_x, p_y) を ε 倍して ε の級数に展開 (δ は厳密に残す)。
(b) 磁場とベクトルポテンシャル: 偏向磁石 A_s = -(x + x^2/(2ρ)) B0 → B_y = B0 (B_x = 0)、四極 A_s = -B'(x^2 - y^2)/2 → B_x = B'y、B_y = B'x。
(c) 偏向磁石: H2 = x^2/(2ρ^2) - δ x/ρ + (p_x^2 + p_y^2)/(2(1+δ))、p_x' = -x/ρ^2 + δ/ρ、x' = p_x/(1+δ)。
(d) 四極: x'' = -K x、y'' = +K y。
(e) 厳密な H~ の運動方程式 (DOP853) と 2 次の展開の差が、振幅の 3 乗で小さくなること。
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)

x, y, px, py, eps = sp.symbols('x y p_x p_y epsilon', real=True)
d = sp.Symbol('delta', positive=True)
rho, K, B0, Bp, Pe = sp.symbols('rho K B_0 B_p P_0_over_e', positive=True)   # Pe = P0/e
h = 1 + x / rho

print('== (a) 展開 (x, y, p_x, p_y を ε 倍) ==')
sq = sp.sqrt((1 + d)**2 - eps**2 * (px**2 + py**2))
Hkin = -(1 + eps * x / rho) * sq                         # -(1 + x/ρ) sqrt(...)
ser = sp.series(Hkin, eps, 0, 4).removeO()
for k in range(4):
    print('  ε^%d の係数: %s' % (k, sp.simplify(ser.coeff(eps, k))))
print('  → ε^0 は定数 -(1+δ) (運動に効かない)、ε^1 は -(x/ρ)(1+δ)、ε^2 は (p_x^2 + p_y^2)/(2(1+δ))、ε^3 は x(p_x^2 + p_y^2)/(2ρ(1+δ)) (3 次。落とす)')

print('\n== (b) 磁場とベクトルポテンシャル (e A_s/P0 = a_s、B/(P0/e) = 規格化した磁場) ==')
As_dip = -(x + x**2 / (2 * rho)) * B0
Bx = sp.diff(As_dip, y) / h
By = -sp.diff(As_dip, x) / h
print('  偏向磁石 A_s = -(x + x^2/(2ρ)) B0: B_x = %s、B_y = %s' % (sp.simplify(Bx), sp.simplify(By)))
As_q = -Bp / 2 * (x**2 - y**2)
print('  四極 A_s = -B′(x^2 - y^2)/2 (ρ→∞): B_x = %s、B_y = %s' % (sp.diff(As_q, y), -sp.diff(As_q, x)))

print('\n== (c) 偏向磁石 ==')
# e B0/P0 = 1/ρ  ⇒ -e A_s/P0 = (x + x^2/(2ρ))/ρ
aS = -(x + x**2 / (2 * rho)) / rho                          # a_s = e A_s/P0 (e B0/P0 = 1/ρ を使った)
Hfull = -aS - h * sq.subs(eps, 1)
H2 = ((x**2 / (2 * rho**2)) - d * x / rho + (px**2 + py**2) / (2 * (1 + d)))
Hser = sp.series(Hfull.subs({x: eps * x, y: eps * y, px: eps * px, py: eps * py}), eps, 0, 3).removeO()

print('  H2 との差 = %s' % sp.simplify(sp.expand((Hser - Hser.coeff(eps, 0)).subs(eps, 1) - H2)))
print('  x′ = ∂H2/∂p_x = %s' % sp.diff(H2, px))
print('  p_x′ = -∂H2/∂x = %s' % sp.simplify(-sp.diff(H2, x)))
print('  p_y′ = -∂H2/∂y = %s' % sp.simplify(-sp.diff(H2, y)))
print('  δ ≪ 1 で x′ ≈ p_x′ = -x/ρ^2 + δ/ρ。符号: x が曲率の中心から外向き、δ > 0 は偏向が弱いので外へ向かう。')

print('\n== (d) 四極 (K = e B′/P0) ==')
aSq = -K / 2 * (x**2 - y**2)                                 # a_s = e A_s/P0
Hq = -aSq + (px**2 + py**2) / (2 * (1 + d))
print('  H2q = %s' % Hq)
print('  p_x′ = %s、p_y′ = %s、x′ (δ=0) = %s、y′ (δ=0) = %s' % (-sp.diff(Hq, x), -sp.diff(Hq, y), sp.simplify((-sp.diff(Hq, x)).subs(d, 0)), sp.simplify((-sp.diff(Hq, y)).subs(d, 0))))

print('\n== (e) 厳密な H~ の運動方程式と 2 次の展開の差 (偏向磁石、ρ = 10 m、δ = 0) ==')
rho_v = 10.0
dv_ = 0.0
Hex_expr = (-aS - h * sp.sqrt((1 + d)**2 - px**2 - py**2)).subs({rho: rho_v, d: dv_})
f_ex = [sp.lambdify((x, y, px, py), q, 'numpy') for q in (sp.diff(Hex_expr, px), sp.diff(Hex_expr, py), -sp.diff(Hex_expr, x), -sp.diff(Hex_expr, y))]
H2n = H2.subs({rho: rho_v, d: dv_})
f_2 = [sp.lambdify((x, y, px, py), q, 'numpy') for q in (sp.diff(H2n, px), sp.diff(H2n, py), -sp.diff(H2n, x), -sp.diff(H2n, y))]


def run(fs, y0, S):
    sol = solve_ivp(lambda s_, v: [f(*v) for f in fs], (0, S), y0, method='DOP853', rtol=1e-13, atol=1e-16)
    return sol.y[:, -1]


S = 3.0
print('  3 m 進んだ後の (x, p_x) の差 (厳密 - 2 次の展開)。初期値は (x, y, p_x, p_y) = (a, 0, a/2, 0) で振幅 a を 2 倍ずつ (差は a の 2 乗で増えるはず):')
prev = None
for a in (1e-3, 2e-3, 4e-3, 8e-3):
    y0 = [a, 0.0, 0.5 * a, 0.0]
    ye = run(f_ex, y0, S)
    y2 = run(f_2, y0, S)
    dd = max(abs(ye[0] - y2[0]), abs(ye[2] - y2[2]))
    print('   a = %.0e: |Δx|, |Δp_x| の最大 = %.3e%s' % (a, dd, '' if prev is None else '  (前の値の %.2f 倍)' % (dd / prev)))
    prev = dd
print('  → 倍率 4 (= 2^2): 落とした項は H の 3 次 (x p^2/(2ρ(1+δ)) など) で、その運動方程式への寄与 (力) は振幅の 2 次。δ x/ρ などの 2 次の項は H2 に厳密に残してある')
