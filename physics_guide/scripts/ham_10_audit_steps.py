# -*- coding: utf-8 -*-
r"""付録 C の途中式の確認 (本文に書き足した展開)。

(a) §6.2: BCH の 3 次。u3、u1u2、u2u1、u1^3 の語ごとの係数の表と、その合計 (8 語)。
(b) §7.1: :f: の括弧に対するライプニッツ則、:f:^n[a,b] の二項展開、exp(:f:)[a,b] = [exp(:f:)a, exp(:f:)b]
    (多項式の f, a, b で 6 次まで)、[Z_x, Z_p] = 1 (f の非線形のべきの次数ごと)。
(c) §7 (C.15): (-:H:)^n z_k = Σ_m [(JA)^n]_{km} z_m (2 自由度、一般の対称行列 A)。
(d) §8.2: ヘルプの式 (QUAD・SEXT・BEND の V) の偏微分。
"""
import itertools
import sympy as sp
import _common as C

C.start(__file__)
R = sp.Rational

# (a) ------------------------------------------------------------------
print('(a) BCH の 3 次の語ごとの係数: u3 - (1/2)(u1 u2 + u2 u1) + (1/3) u1^3')
X, Y = sp.symbols('X Y', commutative=False)


def coeff(expr, word):
    expr = sp.expand(expr)
    tot = 0
    for term in sp.Add.make_args(expr):
        c, nc = term.args_cnc()
        w = ''
        for f in nc:
            b, e = f.as_base_exp()
            w += str(b) * int(e)
        if w == word:
            tot += sp.Mul(*c)
    return tot


u1 = X + Y
u2 = X * X / 2 + X * Y + Y * Y / 2
u3 = X**3 / 6 + X * X * Y / 2 + X * Y * Y / 2 + Y**3 / 6
u1u2 = sp.expand(u1 * u2)
u2u1 = sp.expand(u2 * u1)
u13 = sp.expand(u1 * u1 * u1)
print('  語   u3     u1u2   u2u1   u1^3   合計')
for w in map(''.join, itertools.product('XY', repeat=3)):
    a, b, c, d = (coeff(e, w) for e in (u3, u1u2, u2u1, u13))
    print('  %s  %-6s %-6s %-6s %-6s %s' % (w, a, b, c, d, a - R(1, 2) * (b + c) + R(1, 3) * d))
cm = lambda a, b: sp.expand(a * b - b * a)
e1 = cm(X, cm(X, Y))
e2 = cm(Y, cm(Y, X))
print('  [X,[X,Y]] の語:', e1, '/ [Y,[Y,X]] の語:', e2)

# (b) ------------------------------------------------------------------
print('\n(b) :f: のライプニッツ則と exp(:f:) が括弧を保つこと (1 自由度、多項式の f, a, b)')
x, p, eps = sp.symbols('x p eps')


def br(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x))


f = x**3 + x * p**2 + p**3 / 3
a = x**2 * p + p**3
b = x * p + x**4
ok1 = sp.expand(br(f, br(a, b)) - br(br(f, a), b) - br(a, br(f, b))) == 0
print('  :f:[a,b] = [:f:a, b] + [a, :f:b]:', ok1)


def op(n, g):
    for _ in range(n):
        g = br(f, g)
    return g


allok = True
for n in range(0, 7):
    lhs = op(n, br(a, b))
    rhs = sum(sp.binomial(n, k) * br(op(k, a), op(n - k, b)) for k in range(n + 1))
    allok = allok and sp.expand(lhs - rhs) == 0
print('  n = 0..6 で :f:^n[a,b] = Σ_k C(n,k) [:f:^k a, :f:^(n-k) b]:', allok)
# exp(eps :f:) を eps^6 まで
def expop(g, N=6):
    tot, cur = g, g
    for n in range(1, N + 1):
        cur = br(f, cur)
        tot = tot + eps**n * cur / sp.factorial(n)
    return sp.expand(tot)


# f -> eps f の流れとして、:f: に eps をかける (eps は 1 次の次数のカウント)
def op_eps(g, N=6):
    tot, cur = g, g
    for n in range(1, N + 1):
        cur = br(eps * f, cur)
        tot = tot + cur / sp.factorial(n)
    return sp.expand(tot)


lhs = op_eps(br(a, b))
rhs = br(op_eps(a), op_eps(b))
d = sp.expand(lhs - rhs)
low = sum(d.coeff(eps, k) * eps**k for k in range(0, 7))
print('  exp(eps:f:)[a,b] - [exp(eps:f:)a, exp(eps:f:)b] の eps^0..eps^6 の項 (展開の打ち切り eps^6):', sp.simplify(low) == 0)
Zx, Zp = op_eps(x), op_eps(p)
dd = sp.expand(br(Zx, Zp) - 1)
print('  [Z_x, Z_p] - 1 の eps^1..eps^6 の係数 (すべて 0 のはず):', [sp.simplify(dd.coeff(eps, k)) == 0 for k in range(1, 7)])

# (c) ------------------------------------------------------------------
print('\n(c) (C.15) の帰納法: (-:H:)^n z_k = Σ_m [(JA)^n]_{km} z_m (2 自由度。z = (x, px, y, py)、H = z^T A z / 2)')
xs = sp.symbols('x px y py')
A = sp.Matrix([[2, 1, 0, 3], [1, 5, -1, 0], [0, -1, 7, 2], [3, 0, 2, 4]])
H = sp.expand((sp.Matrix([xs]) * A * sp.Matrix(xs))[0] / 2)
J = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])


def PB(f, g):
    return sp.expand(sum(sp.diff(f, xs[2 * i]) * sp.diff(g, xs[2 * i + 1]) - sp.diff(f, xs[2 * i + 1]) * sp.diff(g, xs[2 * i]) for i in range(2)))


JA = J * A
ok = True
for k in range(4):
    cur = xs[k]
    for n in range(1, 6):
        cur = sp.expand(-PB(H, cur))
        pred = sum((JA**n)[k, m] * xs[m] for m in range(4))
        ok = ok and sp.expand(cur - pred) == 0
print('  k = 0..3、n = 1..5 で (-:H:)^n z_k = Σ_m [(JA)^n]_km z_m:', ok)
print('  z_k\' = [z_k, H] = Σ_l J_kl ∂H/∂z_l (n = 1):', all(sp.expand(PB(xs[k], H) - sum(J[k, l] * sp.diff(H, xs[l]) for l in range(4))) == 0 for k in range(4)))

# (d) ------------------------------------------------------------------
print('\n(d) §8.2 の確認: ヘルプの式の偏微分')
X_, Y_, K1, K2, n_, F1, rb, P = sp.symbols('x y K1 K2 n F1 rho_b p', positive=True)
px, py = sp.symbols('px py')
H2n = K1 * (X_**2 - Y_**2) / (2 * n_)
print('  (ii) QUAD H_2n の K1(x^2-y^2)/(2n): Δpx = -∂H/∂x =', sp.diff(-H2n, X_), '、Δpy = -∂H/∂y =', sp.diff(-H2n, Y_))
Hs = K2 / 6 * sp.re(sp.expand((X_ - sp.I * Y_)**3))
print('  (iii) SEXT: Re(x-iy)^3 =', sp.expand(sp.re(sp.expand((X_ - sp.I * Y_)**3))), '、Δpx =', sp.factor(sp.diff(-Hs, X_)), '、Δpy =', sp.factor(sp.diff(-Hs, Y_)))
V = -(F1**2 / (24 * rb) * px + F1 / (12 * rb**2) * Y_**2 - Y_**4 / (6 * F1 * rb**2)) / P
print('  (iv) BEND: Δx = ∂V/∂px =', sp.simplify(sp.diff(V, px)), '、Δpy = -∂V/∂y =', sp.factor(sp.simplify(-sp.diff(V, Y_))))
print('       y(F1^2 - 4y^2)/(6 F1 p rho_b^2) との差:', sp.simplify(-sp.diff(V, Y_) - Y_ * (F1**2 - 4 * Y_**2) / (6 * F1 * P * rb**2)))
pp, pxx, pyy, Lg = sp.symbols('p px py L', positive=True)
Hd = -sp.sqrt(pp**2 - pxx**2 - pyy**2) * Lg
print('  (i) DRIFT: Δx = ∂H/∂px =', sp.simplify(sp.diff(Hd, pxx)), '(= L px/pz)')
