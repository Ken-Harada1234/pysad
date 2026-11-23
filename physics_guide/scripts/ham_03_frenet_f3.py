# -*- coding: utf-8 -*-
"""Frenet-Serret 座標への正準変換 (母関数 F3(p, Q) = -p·r(Q)、01a_hamiltonian_derivation.md §7) を sympy で確かめる。

水平面内の円弧 (曲率半径 ρ、ねじれ 0): r(X, Y, s) = r0(s) + X n(s) + Y b(s)。
(a) ∂r/∂X = n、∂r/∂Y = b、∂r/∂s = (1 + X/ρ) t。
(b) P = p·∂r/∂Q = (p·n, p·b, (1 + X/ρ) p·t)。
(c) (q, p) → (Q, P) がシンプレクティック (ポアソン括弧)。
(d) |p|^2 = P_X^2 + P_Y^2 + P_s^2/(1 + X/ρ)^2。
(e) ベクトルポテンシャル A も同じ規則で A_s = (1 + X/ρ) A·t。磁場の式 B_X = ∂_Y A_s/(1+X/ρ)、B_Y = -∂_X A_s/(1+X/ρ) を、
    デカルト座標の rot A を直接計算した値と比べる (具体的な多項式の A_t で、いくつかの点)。
"""
import sympy as sp
import _common as C

C.start(__file__)

X, Y, s = sp.symbols('X Y s', real=True)
rho = sp.Symbol('rho', positive=True)
th = s / rho
t_ = sp.Matrix([-sp.sin(th), 0, sp.cos(th)])
n_ = sp.Matrix([sp.cos(th), 0, sp.sin(th)])
b_ = t_.cross(n_)
r0 = sp.Matrix([rho * sp.cos(th), 0, rho * sp.sin(th)])
r = r0 + X * n_ + Y * b_

print('== (a) 基底ベクトルと ∂r/∂Q ==')
print('  t = dr0/ds - %s' % sp.simplify(r0.diff(s) - t_).T)
print('  n = -ρ d^2 r0/ds^2 の差: %s' % sp.simplify(-rho * r0.diff(s, 2) - n_).T)
print('  b = t × n = %s' % sp.simplify(b_).T)
print('  dt/ds + n/ρ = %s、dn/ds - t/ρ = %s、db/ds = %s' % (sp.simplify(t_.diff(s) + n_ / rho).T, sp.simplify(n_.diff(s) - t_ / rho).T, sp.simplify(b_.diff(s)).T))
print('  ∂r/∂X - n = %s、∂r/∂Y - b = %s' % (sp.simplify(r.diff(X) - n_).T, sp.simplify(r.diff(Y) - b_).T))
print('  ∂r/∂s - (1 + X/ρ) t = %s' % sp.simplify(r.diff(s) - (1 + X / rho) * t_).T)

print('\n== (b) P = p·∂r/∂Q (F3 = -p·r(Q) から P = -∂F3/∂Q) ==')
px, py, pz = sp.symbols('p_x p_y p_z', real=True)     # デカルト成分
p = sp.Matrix([px, py, pz])
Q = (X, Y, s)
P = [sp.simplify(p.dot(r.diff(q))) for q in Q]
pn, pb, pt = sp.simplify(p.dot(n_)), sp.simplify(p.dot(b_)), sp.simplify(p.dot(t_))
print('  P_X - p·n = %s、P_Y - p·b = %s、P_s - (1 + X/ρ) p·t = %s' % (sp.simplify(P[0] - pn), sp.simplify(P[1] - pb), sp.simplify(P[2] - (1 + X / rho) * pt)))

print('\n== (c) シンプレクティック性: 古い (q, p) = (r0 + ..., p) のデカルト座標の変数で、ポアソン括弧 ==')
# 古い座標: q = r(X, Y, s) = (qx, qy, qz) をデカルト成分で書き直す。新しい座標 Q は q の関数 (逆写像)。
# ここでは逆に、q = r(Q) と p = (J^-T) P の関係で、{Q_i, P_j}_(q,p) = δ_ij を Jacobian から示す。
J = r.jacobian(sp.Matrix(Q))                       # dq/dQ
Jinv = sp.simplify(J.inv())                        # dQ/dq
# P = J^T p → p = J^-T P。{Q_i, P_j}_(q,p) = Σ_k ∂Q_i/∂q_k ∂P_j/∂p_k = Σ_k (Jinv)_{ik} (J^T)_{jk}
PB_QP = sp.simplify(Jinv * J)                       # (Jinv)_{ik} J_{kj}
print('  {Q_i, P_j} の行列 (単位行列なら OK): %s' % PB_QP.tolist())
# {P_i, P_j}: P_i = Σ_k J_{ki}(q) p_k。{P_i,P_j} = Σ_k (∂P_i/∂q_k ∂P_j/∂p_k - ∂P_i/∂p_k ∂P_j/∂q_k) をデカルトの q について計算する。
# q を (X,Y,s) の関数として扱うため、∂/∂q_k = Σ_m (Jinv)_{mk} ∂/∂Q_m
def dq(f, k):
    return sum(Jinv[m, k] * sp.diff(f, Q[m]) for m in range(3))
def dpk(f, k):
    return sp.diff(f, p[k])
PBPP = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        PBPP[i, j] = sp.simplify(sum(dq(P[i], k) * dpk(P[j], k) - dpk(P[i], k) * dq(P[j], k) for k in range(3)))
print('  {P_i, P_j} = %s (零行列なら OK)' % PBPP.tolist())
PBQQ = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        PBQQ[i, j] = sp.simplify(sum(dq(Q[i], k) * dpk(Q[j], k) - dpk(Q[i], k) * dq(Q[j], k) for k in range(3)))
print('  {Q_i, Q_j} = %s (Q は q だけの関数なので 0)' % PBQQ.tolist())

print('\n== (d) 運動量の大きさ ==')
P2 = P[0]**2 + P[1]**2 + P[2]**2 / (1 + X / rho)**2
print('  |p|^2 - [P_X^2 + P_Y^2 + P_s^2/(1 + X/ρ)^2] = %s' % sp.simplify(p.dot(p) - P2))

print('\n== (e) ベクトルポテンシャルと磁場 ==')
print('  A のカルテシアン成分 A_cart を、Frenet 成分 (A_X, A_Y, A_t) から A = A_X n + A_Y b + A_t t で作ると、')
print('  A_s ≡ A·∂r/∂s = (1 + X/ρ) A_t (P と同じ規則)。')
print('  A_X = A_Y = 0 の場合に、B = rot A をデカルト座標で直接計算し (X, Y, s は x, y, z から作る)、')
print('  B_X = B·n = ∂_Y A_s/(1+X/ρ)、B_Y = B·b = -∂_X A_s/(1+X/ρ)、B·t = 0 と比べる:')
xx, yy, zz = sp.symbols('x y z', real=True)
Xc = sp.sqrt(xx**2 + zz**2) - rho
sc = rho * sp.atan2(zz, xx)
Yc = yy
At_poly = lambda Xv, Yv: sp.Rational(3, 7) * Xv + sp.Rational(5, 3) * Xv**2 - sp.Rational(2, 5) * Yv**2 + Xv * Yv**2 + sp.Rational(1, 4) * Xv**3
Atc = At_poly(Xc, Yc)
th_c = sc / rho
tvec = sp.Matrix([-sp.sin(th_c), 0, sp.cos(th_c)])
nvec = sp.Matrix([sp.cos(th_c), 0, sp.sin(th_c)])
bvec = sp.Matrix([0, 1, 0])
Acart = Atc * tvec
Bcart = sp.Matrix([sp.diff(Acart[2], yy) - sp.diff(Acart[1], zz), sp.diff(Acart[0], zz) - sp.diff(Acart[2], xx), sp.diff(Acart[1], xx) - sp.diff(Acart[0], yy)])
As = (1 + X / rho) * At_poly(X, Y)
BX_f = sp.diff(As, Y) / (1 + X / rho)
BY_f = -sp.diff(As, X) / (1 + X / rho)
rho_v = sp.Rational(10)
worst = 0
for (Xv, Yv, sv) in ((sp.Rational(1, 100), sp.Rational(-2, 100), sp.Rational(3)), (sp.Rational(-3, 100), sp.Rational(1, 50), sp.Rational(7, 2)), (sp.Rational(1, 50), sp.Rational(1, 30), sp.Rational(1, 2))):
    # (X, Y, s) → デカルト点
    pt = r.subs({rho: rho_v, X: Xv, Y: Yv, s: sv})
    subs_c = {xx: pt[0], yy: pt[1], zz: pt[2], rho: rho_v}
    Bn = (Bcart.T * nvec)[0].subs(subs_c).evalf(30)
    Bb = (Bcart.T * bvec)[0].subs(subs_c).evalf(30)
    Bt = (Bcart.T * tvec)[0].subs(subs_c).evalf(30)
    f_n = BX_f.subs({rho: rho_v, X: Xv, Y: Yv}).evalf(30)
    f_b = BY_f.subs({rho: rho_v, X: Xv, Y: Yv}).evalf(30)
    worst = max(worst, abs(Bn - f_n), abs(Bb - f_b), abs(Bt))
    print('  (X, Y, s) = (%s, %s, %s): B·n = %.12f (式 %.12f)、B·b = %.12f (式 %.12f)、B·t = %.2e' % (Xv, Yv, sv, Bn, f_n, Bb, f_b, Bt))
print('  最大差 (B·t は 0 のはず) = %.2e' % worst)
