# -*- coding: utf-8 -*-
"""ルジャンドル変換: L = -m0 c^2/γ - e(φ - A·v) から H = c sqrt((p - eA)^2 + m0^2 c^2) + eφ を作る
(01a_hamiltonian_derivation.md §5)。sympy で確かめる。

(a) p = γ m0 v + eA から v を p で解く: v = (p - eA)/(γ m0)、γ = sqrt(1 + |p - eA|^2/(m0 c)^2)。
(b) H = p·v - L がその v で c sqrt(...) + eφ になる。
(c) ∂H/∂p = v、-∂H/∂r = ∂L/∂r (v 固定の偏微分と、p 固定の偏微分の関係)。
"""
import sympy as sp
import _common as C

C.start(__file__)

m0, c, e = sp.symbols('m_0 c e', positive=True)
phi = sp.Symbol('phi', real=True)
x, y, z = sp.symbols('x y z', real=True)
R = (x, y, z)
p = sp.symbols('p_x p_y p_z', real=True)
A = [sp.Function('A_' + n)(x, y, z) for n in 'xyz']
pi = [p[i] - e * A[i] for i in range(3)]               # 機械運動量 γ m0 v
pi2 = sum(q**2 for q in pi)
S = sp.sqrt(pi2 + m0**2 * c**2)                          # = γ m0 c
gam = S / (m0 * c)
v = [pi[i] / (gam * m0) for i in range(3)]               # (a) で解いた v

print('== (a) p = γ m0 v + eA を v について解く ==')
print('  γ = sqrt(1 + |p - eA|^2/(m0 c)^2) と置いて v = (p - eA)/(γ m0) とすると、')
v2 = sum(q**2 for q in v)
print('  1 - v^2/c^2 - 1/γ^2 = %s  (つまり γ の定義 1/sqrt(1 - v^2/c^2) と矛盾しない)' % sp.simplify(1 - v2 / c**2 - 1 / gam**2))
print('  γ m0 v + eA - p = %s' % [sp.simplify(gam * m0 * v[i] + e * A[i] - p[i]) for i in range(3)])

print('\n== (b) H = p·v - L ==')
L = -m0 * c**2 / gam - e * (phi - sum(A[i] * v[i] for i in range(3)))
H_leg = sum(p[i] * v[i] for i in range(3)) - L
H = c * S + e * phi
print('  H_ルジャンドル - [c sqrt((p-eA)^2 + m0^2 c^2) + eφ] = %s' % sp.simplify(H_leg - H))
print('  途中: p·v - L = γ m0 v·v + m0 c^2/γ + eφ (eA·v は相殺)。')
print('        つまり γ m0 v^2 + m0 c^2/γ = γ m0 c^2 (差 %s)' % sp.simplify(gam * m0 * v2 + m0 * c**2 / gam - gam * m0 * c**2))
print('        γ m0 c^2 = c sqrt((p - eA)^2 + m0^2 c^2) (差 %s)' % sp.simplify(gam * m0 * c**2 - c * S))

print('\n== (c) ハミルトンの方程式との整合 ==')
print('  ∂H/∂p_i - v_i = %s' % [sp.simplify(sp.diff(H, p[i]) - v[i]) for i in range(3)])
vs = sp.symbols('v_x v_y v_z', real=True)
Lv = -m0 * c**2 * sp.sqrt(1 - sum(q**2 for q in vs) / c**2) - e * (phi - sum(A[i] * vs[i] for i in range(3)))
for i in range(3):
    dL = sp.diff(Lv, R[i]).subs(dict(zip(vs, v)), simultaneous=True)
    dH = -sp.diff(H, R[i])
    print('  成分 %s: ∂L/∂r_i (v 固定、後で v(p) を代入) - [-∂H/∂r_i (p 固定)] = %s' % ('xyz'[i], sp.simplify(dL - dH)))
print('  → ṗ = ∂L/∂r = -∂H/∂r。ラグランジュの方程式とハミルトンの方程式は同じ運動')
