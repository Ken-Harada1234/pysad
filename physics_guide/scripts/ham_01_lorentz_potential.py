# -*- coding: utf-8 -*-
"""ローレンツ力が一般化ポテンシャル U = e(φ - A·v) で表せること、ラグランジュの運動方程式と一般化運動量 (01a_hamiltonian_derivation.md §3〜§4) を sympy で確かめる。

(a) ベクトル恒等式を成分で: ∇(a·b) の 4 項。v が r と独立なら 2 項が消える。v×(∇×A) = ∇(v·A) - (v·∇)A。
(b) 全微分 dA/dt = ∂A/∂t + (v·∇)A (A は r(t) と t の関数)。
(c) 一般化力 Q = -∂U/∂r + d/dt(∂U/∂v) が e(E + v×B) と一致する (E = -∇φ - ∂A/∂t, B = ∇×A)。
(d) p = ∂L/∂v = γ m0 v + eA。ラグランジュの運動方程式が d(γ m0 v)/dt = e(E + v×B) になる。
"""
import sympy as sp
import _common as C

C.start(__file__)

x, y, z, t = sp.symbols('x y z t', real=True)
vx, vy, vz = sp.symbols('v_x v_y v_z', real=True)
ax, ay, az = sp.symbols('a_x a_y a_z', real=True)       # 加速度 dv/dt
e, m0, c = sp.symbols('e m_0 c', positive=True)
R = (x, y, z)
V = (vx, vy, vz)
Acc = (ax, ay, az)

phi = sp.Function('phi')(x, y, z, t)
A = [sp.Function('A_' + n)(x, y, z, t) for n in 'xyz']


def grad(f):
    return [sp.diff(f, q) for q in R]


def curl(F):
    return [sp.diff(F[2], y) - sp.diff(F[1], z), sp.diff(F[0], z) - sp.diff(F[2], x), sp.diff(F[1], x) - sp.diff(F[0], y)]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(i * j for i, j in zip(a, b))


def zero(expr):
    return sp.simplify(sp.expand(expr)) == 0


print('== (a) ベクトル恒等式 ∇(a·b) = (a·∇)b + (b·∇)a + a×(∇×b) + b×(∇×a) を成分で ==')
# a = u(r) (r に依存する一般のベクトル場)、b = A
u = [sp.Function('u_' + n)(x, y, z) for n in 'xyz']
Ar = [sp.Function('B_' + n)(x, y, z) for n in 'xyz']      # r だけの関数 (恒等式の確認用)
lhs = grad(dot(u, Ar))
adv = lambda a, b: [sum(a[j] * sp.diff(b[i], R[j]) for j in range(3)) for i in range(3)]
t1 = adv(u, Ar)
t2 = adv(Ar, u)
t3 = cross(u, curl(Ar))
t4 = cross(Ar, curl(u))
ok = all(zero(lhs[i] - t1[i] - t2[i] - t3[i] - t4[i]) for i in range(3))
print('  4 項の恒等式 (3 成分): %s' % ('成立' if ok else '不成立'))
print('  x 成分の中身: (a·∇)b_x + (b·∇)a_x + [a×(∇×b)]_x + [b×(∇×a)]_x。第 2・4 項は a の r 微分を含む')

print('\n  v が r と独立 (∂v/∂x_j = 0): 第 2・4 項が消えて ∇(v·A) = (v·∇)A + v×(∇×A) になる')
Ax = A
lhs = grad(dot(V, Ax))                      # V は記号 (r に依存しない)
rhs = [adv(V, Ax)[i] + cross(V, curl(Ax))[i] for i in range(3)]
t2v = adv(Ax, V)
t4v = cross(Ax, curl(V))
print('  (A·∇)v の 3 成分: %s、A×(∇×v) の 3 成分: %s' % (t2v, t4v))
print('  ∇(v·A) - [(v·∇)A + v×(∇×A)] = %s' % [sp.simplify(lhs[i] - rhs[i]) for i in range(3)])
print('  x 成分の展開: [v×(∇×A)]_x = v_y(∂_x A_y - ∂_y A_x) + v_z(∂_x A_z - ∂_z A_x)')
vxB = cross(V, curl(Ax))
chk = vy * (sp.diff(A[1], x) - sp.diff(A[0], y)) + vz * (sp.diff(A[2], x) - sp.diff(A[0], z))
print('    手で書いた式との差: %s' % sp.simplify(vxB[0] - chk))

print('\n== (b) 全微分 dA/dt = ∂A/∂t + (v·∇)A (A は r(t) を通じて t に依存) ==')
tt = sp.Symbol('tau', real=True)
rt = [sp.Function('r_' + n)(tt) for n in 'xyz']
a1, a2, a3 = sp.symbols('a1 a2 a3')
# 成分ごとに別の形の関数を使い、A(r(τ), τ) を τ で微分した結果が ∂_t A + (v·∇)A と同じことを確かめる
samples = [x**2 * y + a1 * z * t + x * t**2, sp.sin(x * y) + a2 * z**3 * t, sp.exp(z) * x + a3 * y * t**3]
for i, Af in enumerate(samples):
    lhs = sp.diff(Af.subs({x: rt[0], y: rt[1], z: rt[2], t: tt}, simultaneous=True), tt)
    vv = [sp.diff(q, tt) for q in rt]
    rhs = (sp.diff(Af, t) + sum(sp.diff(Af, R[j]) * vv[j] for j in range(3))).subs({x: rt[0], y: rt[1], z: rt[2], t: tt}, simultaneous=True)
    print('  例 A_%s = %s : 全微分 - (∂_t A + v·∇A) = %s' % ('xyz'[i], Af, sp.simplify(lhs - rhs)))
print('  → dA/dt = ∂_t A + v_x ∂_x A + v_y ∂_y A + v_z ∂_z A (成分ごとの鎖則)')

print('\n== (c) 一般化力 Q = -∂U/∂r + d/dt(∂U/∂v)、U = e(φ - A·v) ==')
U = e * (phi - dot(A, V))


def Dt(f):
    """全時間微分 (f は r, v, t の関数。∂/∂t + v·∇_r + a·∇_v)"""
    return sp.diff(f, t) + sum(V[j] * sp.diff(f, R[j]) for j in range(3)) + sum(Acc[j] * sp.diff(f, V[j]) for j in range(3))


E = [-g - sp.diff(A[i], t) for i, g in enumerate(grad(phi))]
B = curl(A)
F = [e * (E[i] + cross(V, B)[i]) for i in range(3)]
Q = [-sp.diff(U, R[i]) + Dt(sp.diff(U, V[i])) for i in range(3)]
print('  ∂U/∂v_i = -e A_i : %s' % [sp.simplify(sp.diff(U, V[i]) + e * A[i]) for i in range(3)])
for i in range(3):
    print('  成分 %s: Q - e(E + v×B) = %s' % ('xyz'[i], sp.simplify(sp.expand(Q[i] - F[i]))))
print('  第 2 項だけ (d/dt ∂U/∂v) を落とすと、差は x 成分で %s' % sp.simplify(sp.expand(-sp.diff(U, x) - F[0])))
print('  → 落とした差は e dA/dt = e(∂_t A + (v·∇)A) そのもの。追加項 d/dt(∂U/∂v) = -e dA/dt が必要')

print('\n== (d) 相対論的ラグランジアン L = -m0 c^2 sqrt(1 - v^2/c^2) - e(φ - A·v) ==')
v2 = dot(V, V)
L = -m0 * c**2 * sp.sqrt(1 - v2 / c**2) - e * (phi - dot(A, V))
gam = 1 / sp.sqrt(1 - v2 / c**2)
p = [sp.diff(L, V[i]) for i in range(3)]
print('  p_i = ∂L/∂v_i - (γ m0 v_i + e A_i) = %s' % [sp.simplify(p[i] - gam * m0 * V[i] - e * A[i]) for i in range(3)])
print('  → 正準運動量 p = γ m0 v + e A (機械運動量 γ m0 v に e A が加わる)')
EL = [Dt(p[i]) - sp.diff(L, R[i]) for i in range(3)]       # = 0 が運動方程式
mech = [Dt(gam * m0 * V[i]) for i in range(3)]               # d(γ m0 v)/dt
for i in range(3):
    print('  成分 %s: [d/dt(∂L/∂v) - ∂L/∂r] - [d(γ m0 v)/dt - e(E + v×B)] = %s' % ('xyz'[i], sp.simplify(sp.expand(EL[i] - (mech[i] - F[i])))))
print('  → EL 方程式は d(γ m0 v)/dt = e(E + v×B) (ローレンツ力の運動方程式) と同じ')
