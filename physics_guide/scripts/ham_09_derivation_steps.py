# -*- coding: utf-8 -*-
r"""01a 章の途中式の確認 (本文に書き足した展開)。

(a) (1a.13): L = T - U をラグランジュの方程式に入れた結果が、
    [d/dt dT/dq' - dT/dq] - [-dU/dq + d/dt dU/dq'] になること (T, U は一般の関数)。
(b) (1a.38): 平面曲線 (角度 θ(s)、曲率 1/ρ = θ') で、t' = -n/ρ、n' = t/ρ、b' = 0。
(c) (1a.21) の x 成分の組み立て: [a×(∇×b)]_x + [b×(∇×a)]_x + (a·∇)b_x + (b·∇)a_x = ∂_x(a·b)。
(d) (1a.60a)(1a.60b)(1a.65): 曲線座標の回転 (スケール因子 1, 1, h) から B、および偏向・四極の A_s の積分。
"""
import sympy as sp
import _common as C

C.start(__file__)

# (a) ------------------------------------------------------------------
print('(a) (1a.13)')
t = sp.symbols('t')
q = sp.Function('q')(t)
qs, qds = sp.symbols('qs qds')
T = sp.Function('T')
U = sp.Function('U')
Tq = T(qs, qds); Uq = U(qs, qds)
Lq = Tq - Uq
sub = {qs: q, qds: q.diff(t)}
EL_L = sp.diff(sp.diff(Lq, qds).subs(sub), t) - sp.diff(Lq, qs).subs(sub)
EL_T = sp.diff(sp.diff(Tq, qds).subs(sub), t) - sp.diff(Tq, qs).subs(sub)
rhs = -sp.diff(Uq, qs).subs(sub) + sp.diff(sp.diff(Uq, qds).subs(sub), t)
print('  EL(L) = [d/dt dT/dq\' - dT/dq] - [-dU/dq + d/dt dU/dq\']:', sp.simplify(sp.expand(EL_L - (EL_T - rhs))) == 0)

# (b) ------------------------------------------------------------------
print('(b) Frenet-Serret (平面曲線)')
s = sp.symbols('s')
th = sp.Function('theta')(s)
rho = 1 / th.diff(s)
r0p = sp.Matrix([sp.cos(th), 0, sp.sin(th)])          # t = r0'
tvec = r0p
n_vec = sp.simplify(-rho * tvec.diff(s))                 # n = -ρ r0''
b_vec = tvec.cross(n_vec)
print('  n =', list(n_vec), '、b =', list(sp.simplify(b_vec)))
print('  |n| = 1:', sp.simplify(n_vec.dot(n_vec)) == 1)
print('  t\' = -n/ρ:', sp.simplify(tvec.diff(s) + n_vec / rho) == sp.zeros(3, 1))
print('  n\' = t/ρ:', sp.simplify(n_vec.diff(s) - tvec / rho) == sp.zeros(3, 1))
print('  b\' = 0:', sp.simplify(b_vec.diff(s)) == sp.zeros(3, 1))
print('  n\'·t = -n·t\' = 1/ρ、n\'·b = 0:', sp.simplify(n_vec.diff(s).dot(tvec) - 1 / rho) == 0, sp.simplify(n_vec.diff(s).dot(b_vec)) == 0)

# (c) ------------------------------------------------------------------
print('(c) (1a.21) の x 成分')
x, y, z = sp.symbols('x y z')
a = [sp.Function('a%d' % i)(x, y, z) for i in range(3)]
b = [sp.Function('b%d' % i)(x, y, z) for i in range(3)]
X = [x, y, z]
d = lambda f, i: sp.diff(f, X[i])
curl = lambda v: [d(v[2], 1) - d(v[1], 2), d(v[0], 2) - d(v[2], 0), d(v[1], 0) - d(v[0], 1)]
cross = lambda u, w: [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
adotgrad = lambda u, f: sum(u[i] * d(f, i) for i in range(3))
ca = cross(a, curl(b)); cb = cross(b, curl(a))
tot = ca[0] + cb[0] + adotgrad(a, b[0]) + adotgrad(b, a[0])
grad_ab_x = d(sum(a[i] * b[i] for i in range(3)), 0)
print('  x 成分の和 = ∂_x(a·b):', sp.simplify(sp.expand(tot - grad_ab_x)) == 0)

# (d) ------------------------------------------------------------------
print('(d) 曲線座標の回転と A_s の積分')
X_, Y_, rho_, B0, Bp = sp.symbols('X Y rho B0 Bp', positive=True)
h = 1 + X_ / rho_
As = sp.Function('As')(X_, Y_)
# 直交曲線座標 (スケール因子 hX=1, hY=1, hs=h)。A の物理成分は A_X=A_Y=0、A_t = As/h
# (∇×A)_X = (1/(hY hs)) [∂(hs A_t)/∂Y - ∂(hY A_Y)/∂s]、(∇×A)_Y = (1/(hs hX)) [∂(hX A_X)/∂s - ∂(hs A_t)/∂X]
BX = sp.diff(h * (As / h), Y_) / h
BY = -sp.diff(h * (As / h), X_) / h
print('  B_X = (1/h) ∂A_s/∂Y:', sp.simplify(BX - sp.diff(As, Y_) / h) == 0, '、B_Y = -(1/h) ∂A_s/∂X:', sp.simplify(BY + sp.diff(As, X_) / h) == 0)
# 偏向: A_s = -(X + X^2/(2ρ)) B0
As_d = -(X_ + X_**2 / (2 * rho_)) * B0
print('  偏向: B_X =', sp.simplify(sp.diff(As_d, Y_) / h), '、B_Y =', sp.simplify(-sp.diff(As_d, X_) / h))
# 四極 (ρ→∞): A_s = -(B'/2)(x^2 - y^2)
x_, y_ = sp.symbols('x_ y_')
As_q = -Bp / 2 * (x_**2 - y_**2)
print('  四極: B_X = ∂_y A_s =', sp.diff(As_q, y_), '、B_Y = -∂_x A_s =', -sp.diff(As_q, x_))
