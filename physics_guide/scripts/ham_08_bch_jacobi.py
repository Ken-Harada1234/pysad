# -*- coding: utf-8 -*-
r"""付録 C の途中式の確認 (本文に書き足した展開)。

(a) ヤコビ恒等式 [f,[g,h]] + [g,[h,f]] + [h,[f,g]] = 0 (1 自由度の一般の関数 f, g, h で sympy)。
    これから [:f:, :g:] = :[f,g]: (リー演算子の交換子)。
(b) BCH の 3 次まで: 非可換記号 X, Y で log(e^X e^Y) の語ごとの係数を数え、
    (1/2)[X,Y] と (1/12)([X,[X,Y]] + [Y,[Y,X]]) に一致すること。
(c) (C.7) の帰納法: g^{(n)} = (-:H:)^n g をドリフトで確かめる (H = p^2/2、g = x^2)。
(d) (C.9) の偶奇の規則: :f:^{2m} x = (-k s^2)^m x、:f:^{2m+1} x = (-k s^2)^m s p など。
"""
import itertools
import sympy as sp
import _common as C

C.start(__file__)

x, p = sp.symbols('x p')


def br(f, g):
    return sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x)


# (a) ------------------------------------------------------------------
print('(a) ヤコビ恒等式と [:f:,:g:] = :[f,g]:')
f, g, h = (sp.Function(n)(x, p) for n in 'fgh')
jac = br(f, br(g, h)) + br(g, br(h, f)) + br(h, br(f, g))
print('  Jacobi (一般の f, g, h):', sp.simplify(sp.expand(jac)) == 0)
lhs = br(f, br(g, h)) - br(g, br(f, h))     # [:f:,:g:] h
rhs = br(br(f, g), h)                        # :[f,g]: h
print('  [:f:,:g:] h = :[f,g]: h:', sp.simplify(sp.expand(lhs - rhs)) == 0)

# (b) ------------------------------------------------------------------
print('(b) BCH (3 次)')
X, Y = sp.symbols('X Y', commutative=False)


def trunc(expr, n):
    """非可換多項式を n 次以下に切る"""
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        c, nc = term.args_cnc()
        deg = sum(abs(e) if e.is_number else 0 for b, e in (a.as_base_exp() for a in nc))
        if deg <= n:
            out += term
    return sp.expand(out)


def exp_series(Z, n):
    out, term = 1, 1
    for k in range(1, n + 1):
        term = sp.expand(term * Z / k)
        out += term
    return out


E = trunc(exp_series(X, 3) * exp_series(Y, 3), 3)
u = E - 1
logE = trunc(u - u * u / 2 + u * u * u / 3, 3)
comm = lambda a, b: sp.expand(a * b - b * a)
bch = sp.expand(X + Y + comm(X, Y) / 2 + (comm(X, comm(X, Y)) + comm(Y, comm(Y, X))) / 12)
print('  log(e^X e^Y) の 3 次までの各語の係数:')
for term in sp.Add.make_args(logE):
    print('    ', term)
print('  log(e^X e^Y) - [X + Y + [X,Y]/2 + ([X,[X,Y]] + [Y,[Y,X]])/12] = ', sp.expand(logE - bch))
print('  [X,[X,Y]] =', comm(X, comm(X, Y)), '、[Y,[Y,X]] =', comm(Y, comm(Y, X)))

# (c) ------------------------------------------------------------------
print('(c) (C.7) の帰納法 (H = p^2/2、g = x^2): d^n g/ds^n = (-:H:)^n g')
H = p**2 / 2
gg = x**2
cur = gg
s = sp.symbols('s')
xs, ps = x + s * p, p        # 厳密解 x(s) = x + s p
gs = gg.subs(x, xs)
for n in range(0, 4):
    dn = sp.diff(gs, s, n).subs(s, 0)
    print('  n=%d: d^n g/ds^n|0 = %s、(-:H:)^n g = %s' % (n, sp.expand(dn), sp.expand(cur)))
    cur = sp.expand(-br(H, cur))

# (d) ------------------------------------------------------------------
print('(d) (C.9) の偶奇の規則 (f = -(s/2)(k x^2 + p^2))')
k = sp.symbols('k')
fq = -s / 2 * (k * x**2 + p**2)
cx, cp = x, p
ok = True
for n in range(0, 7):
    m = n // 2
    if n % 2 == 0:
        ex, ep = (-k * s**2)**m * x, (-k * s**2)**m * p
    else:
        ex, ep = (-k * s**2)**m * s * p, (-k * s**2)**m * (-k * s) * x
    ok &= sp.simplify(cx - ex) == 0 and sp.simplify(cp - ep) == 0
    cx, cp = br(fq, cx), br(fq, cp)
print('  n = 0..6 で規則と一致:', ok)
