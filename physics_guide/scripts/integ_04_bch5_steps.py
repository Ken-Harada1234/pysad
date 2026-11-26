# -*- coding: utf-8 -*-
r"""03 章 §5・§8.3 の 5 次の BCH を、手で追える形に展開した確認 (自由 Lie 代数)。

(a) 次数 5 の基底 c1..c6 と、ヤコビ恒等式による簡約の関係 (J1)〜(J5)。
(b) 補題 S: log(e^{X/2} e^{Y} e^{X/2}) の 5 次の項 C5 (6 個の係数)。
    (4,1) と (1,4) は 1 つの文字について線形な厳密式 (x/2)/sinh(x/2)、(x/2)coth(x/2) から。
    (3,2) は Y の 2 次の摂動、(2,3) は X の 2 次の摂動 (微分方程式 Z' = F(ad_Z) h、F(u) = u/(e^u - 1))。
    各項を 1 つずつ書き出し、基底に直した値を示す。直接の log (自由代数) と照合する。
(c) Yoshida の S4 = S2(w1) S2(w0) S2(w1) の 5 次 (3.7): 補題 S の C3 の交差項だけが残る。
(d) 力勾配型 (3.9) を 3 段 (内側・中間・外側) で組み、5 次の 6 個の係数を t, c3, c5 の式で出す。(3.13) と照合。
(e) W5 についての主張: [B,[B,[A,[A,B]]]] = b4 - b5、多自由度の W と W5 の関数。
(f) 多極の補正キックを厳密なドリフトの分割と比べるときの係数 (質量 → p_z)。
"""
import itertools
import sympy as sp
import _common as C

C.start(__file__)
R = sp.Rational
MAXD = 5


# ---- 自由結合代数 (語 → 係数) ----
def padd(a, b, sa=1, sb=1):
    r = {}
    for w, c in a.items():
        r[w] = r.get(w, 0) + sa * c
    for w, c in b.items():
        r[w] = r.get(w, 0) + sb * c
    return {w: sp.expand(c) for w, c in r.items() if sp.expand(c) != 0}


def pmul(a, b):
    r = {}
    for w1, c1 in a.items():
        for w2, c2 in b.items():
            if len(w1) + len(w2) <= MAXD:
                r[w1 + w2] = r.get(w1 + w2, 0) + c1 * c2
    return {w: sp.expand(c) for w, c in r.items() if sp.expand(c) != 0}


def sc(a, s):
    return {w: sp.expand(c * s) for w, c in a.items() if sp.expand(c * s) != 0}


def comm(a, b):
    return padd(pmul(a, b), pmul(b, a), 1, -1)


def lsum(*terms):
    r = {}
    for t in terms:
        r = padd(r, t)
    return r


ONE = {'': sp.Integer(1)}


def pexp(a):
    r = dict(ONE)
    term = dict(ONE)
    for n in range(1, MAXD + 1):
        term = sc(pmul(term, a), R(1, n))
        r = padd(r, term)
    return r


def plog(s):
    X = padd(s, ONE, 1, -1)
    r = {}
    term = dict(ONE)
    for n in range(1, MAXD + 1):
        term = pmul(term, X)
        r = padd(r, sc(term, R((-1)**(n + 1), n)))
    return r


def deg(a, n):
    return {w: c for w, c in a.items() if len(w) == n}


def bideg(a, nx, ny, lx='X', ly='Y'):
    return {w: c for w, c in a.items() if w.count(lx) == nx and w.count(ly) == ny}


def ad(a, b):
    return comm(a, b)


def adn(a, n, b):
    for _ in range(n):
        b = comm(a, b)
    return b


def is_zero(a):
    return all(sp.simplify(c) == 0 for c in a.values())


def letters(lx, ly):
    return {lx: sp.Integer(1)}, {ly: sp.Integer(1)}


# ---- 次数 5 の基底 (Xe, Ye は Lie 元) ----
def parts(Xe, Ye):
    U = comm(Xe, Ye)
    V = comm(Xe, U)
    W = comm(Ye, U)
    c = [comm(Xe, comm(Xe, V)),   # c1
         comm(Ye, comm(Xe, V)),   # c2
         comm(U, V),              # c3
         comm(Xe, comm(Ye, W)),   # c4
         comm(U, W),              # c5
         comm(Ye, comm(Ye, W))]   # c6
    return U, V, W, c


def coords(L, basis, syms=None):
    """L = Σ k_i basis[i] を満たす k_i (L は基底の 1 次結合のはず)"""
    ks = sp.symbols('k1:%d' % (len(basis) + 1))
    expr = {}
    for k, b in zip(ks, basis):
        expr = padd(expr, sc(b, k))
    diff = padd(L, expr, 1, -1)
    eqs = list(diff.values())
    sol = sp.solve(eqs, ks, dict=True)
    if not sol:
        return None
    s = sol[0]
    return [sp.simplify(s.get(k, k)) for k in ks]


def fmt(v):
    return str(sp.nsimplify(v)) if v.is_Number else str(sp.factor(v))


X, Y = letters('X', 'Y')

# =====================================================================
print('== (a) 次数 5 の基底とヤコビ恒等式 ==')
U, V, W, cb = parts(X, Y)
print('  U=[X,Y]、V=[X,U]、W=[Y,U]。次数 4: [X,V]、[Y,V]、[X,W]、[Y,W]')
J1 = padd(comm(X, W), comm(Y, V), 1, -1)
print('  (J1) [X,W] = [Y,V]:', is_zero(J1))
J2 = padd(comm(X, comm(Y, V)), padd(cb[2], cb[1]), 1, -1)
print('  (J2) [X,[Y,V]] = c3 + c2:', is_zero(J2))
J3 = padd(comm(Y, comm(Y, V)), padd(cb[3], cb[4], 1, -1), 1, -1)
print('  (J3) [Y,[Y,V]] = c4 - c5:', is_zero(J3))
J4 = padd(comm(X, comm(X, W)), padd(cb[2], cb[1]), 1, -1)
print('  (J4) [X,[X,W]] = c3 + c2:', is_zero(J4))
J5 = padd(comm(Y, comm(X, W)), padd(cb[3], cb[4], 1, -1), 1, -1)
print('  (J5) [Y,[X,W]] = c4 - c5:', is_zero(J5))
# 次数 5 の Lie 元の次元 = 6 (基底が語の空間で 1 次独立)
ks = sp.symbols('q1:7')
expr = {}
for k, b in zip(ks, cb):
    expr = padd(expr, sc(b, k))
print('  c1..c6 は 1 次独立 (Σ q_i c_i = 0 の解は自明のみ):', sp.solve(list(expr.values()), ks, dict=True) == [{k: 0 for k in ks}])
# 次数 5 の右ネストした全ブラケット (2^5 通り) が c1..c6 の 1 次結合
bad = 0
for ws in itertools.product('XY', repeat=5):
    el = {ws[-1]: sp.Integer(1)}
    for ch in reversed(ws[:-1]):
        el = comm({ch: sp.Integer(1)}, el)
    if el and coords(el, cb) is None:
        bad += 1
print('  右ネストの 32 個の括弧がすべて c1..c6 の 1 次結合になる (外れ %d 個):' % bad, bad == 0)

# =====================================================================
print('\n== (b) 補題 S: log(e^{X/2} e^{Y} e^{X/2}) の 5 次 ==')
S = pmul(pmul(pexp(sc(X, R(1, 2))), pexp(Y)), pexp(sc(X, R(1, 2))))
LS = plog(S)
print('  2 次 = 0:', deg(LS, 2) == {}, '、4 次 = 0:', deg(LS, 4) == {})
C3d = deg(LS, 3)
chk = padd(C3d, padd(sc(V, R(-1, 24)), sc(W, R(-1, 12))), 1, -1)
print('  3 次 = -V/24 - W/12 (V=[X,[X,Y]]、W=[Y,[X,Y]]):', chk == {})
C5d = deg(LS, 5)
k_direct = coords(C5d, cb)
print('  直接計算 (自由代数) の 5 次の係数 (c1..c6):', [fmt(k) for k in k_direct])

# (b-1) (4,1): Y について線形
u = sp.symbols('u')
print('  (4,1): (u/2)/sinh(u/2) の級数:', sp.series((u / 2) / sp.sinh(u / 2), u, 0, 6))
print('  (1,4): (u/2)coth(u/2) の級数:', sp.series((u / 2) / sp.tanh(u / 2), u, 0, 6))
kk = [R(7, 5760), None, None, None, None, R(1, 720)]
# [Y,[Y,[Y,[Y,X]]]] = -c6
print('  ad_Y^4 X = -c6:', is_zero(padd(adn(Y, 4, X), cb[5])))
print('  ad_X^4 Y = c1:', is_zero(padd(adn(X, 4, Y), cb[0], 1, -1)))


# (b-2) (3,2): Y の 2 次。Z' = F(ad_Z) g、g = e^{ad_X/2} Y
F = {0: R(1), 1: R(-1, 2), 2: R(1, 12), 3: R(0), 4: R(-1, 720)}
Z1_Y = {0: R(1), 2: R(-1, 24), 4: R(7, 5760)}           # Z1 = Σ z_m ad_X^m Y
g_Y = {k: R(1, 2**k * sp.factorial(k)) for k in range(0, 5)}   # g = Σ g_k ad_X^k Y
print('\n  -- (3,2): S(ε) = e^{X/2} e^{εY} e^{X/2}、Z = X + εZ1 + ε²Z2。2 Z2 = Σ_n f_n Σ_{i+j=n-1} ad_X^i ad_{Z1} ad_X^j g --')
print('     Z1 = Σ z_m ad_X^m Y (z0=1, z2=-1/24, z4=7/5760)、g = Σ g_k ad_X^k Y (g_k = 1/(2^k k!))')
terms32 = []
for n in (1, 2, 4):
    for i in range(0, n):
        j = n - 1 - i
        for m, zm in Z1_Y.items():
            for k, gk in g_Y.items():
                if i + j + m + k != 3:
                    continue
                z1 = sc(adn(X, m, Y), zm)
                gg = sc(adn(X, k, Y), gk)
                el = adn(X, j, gg)
                el = comm(z1, el)
                el = adn(X, i, el)
                el = sc(el, F[n])
                if not el:
                    continue
                cc = coords(el, cb)
                terms32.append((n, i, j, m, k, el, cc))
tot32 = {}
for (n, i, j, m, k, el, cc) in terms32:
    print('     f_%d=%s  ad_X^%d [z%d Z1part, ad_X^%d (g%d)]: 係数込みの基底での値 (c1..c6) = %s' % (
        n, F[n], i, m, j, k, [fmt(v) for v in cc]))
    tot32 = padd(tot32, el)
half32 = sc(tot32, R(1, 2))
cc32 = coords(half32, cb)
print('     Z2 の (3,2) 部分 = (1/2)Σ: c1..c6 の係数 =', [fmt(v) for v in cc32])
print('     直接計算の (3,2) 部分      =', [fmt(v) for v in coords(bideg(C5d, 3, 2), cb)])

# (b-3) (2,3): X の 2 次。Z(0) = Y、Z' = F(ad_Z) h、h = X/2 + e^{ε ad_X/2} e^{ad_Y} X/2
Z1_X = {0: R(1), 2: R(1, 12), 4: R(-1, 720)}           # Z1 = Σ z_m ad_Y^m X ((u/2)coth(u/2))
h0_X = {0: R(1)}
for k in range(1, 5):
    h0_X[k] = R(1, 2 * sp.factorial(k))                  # h0 = (X + e^{ad_Y} X)/2
print('\n  -- (2,3): S(ε) = e^{εX/2} e^{Y} e^{εX/2}、Z = Y + εZ1 + ε²Z2。dS S^-1 = h = X/2 + e^{ε ad_X/2} e^{ad_Y} X/2 --')
print('     Z1 = Σ z_m ad_Y^m X (z0=1, z2=1/12, z4=-1/720)、h0 = Σ η_k ad_Y^k X (η0=1, η_k=1/(2 k!))')
print('     2 Z2 = Σ_n f_n Σ_{i+j=n-1} ad_Y^i ad_{Z1} ad_Y^j h0 + Σ_n f_n ad_Y^n h1、h1 = (1/4) ad_X e^{ad_Y} X')
terms23 = []
tot23 = {}
for n in (1, 2, 4):
    for i in range(0, n):
        j = n - 1 - i
        for m, zm in Z1_X.items():
            for k, hk in h0_X.items():
                if i + j + m + k != 3:
                    continue
                z1 = sc(adn(Y, m, X), zm)
                hh = sc(adn(Y, k, X), hk)
                el = adn(Y, j, hh)
                el = comm(z1, el)
                el = adn(Y, i, el)
                el = sc(el, F[n])
                if not el:
                    continue
                cc = coords(el, cb)
                print('     [I]  f_%d=%s  ad_Y^%d [z%d, ad_Y^%d (η%d)]: c1..c6 = %s' % (n, F[n], i, m, j, k, [fmt(v) for v in cc]))
                tot23 = padd(tot23, el)
# 第 2 部: F(ad_Y) h1、h1 = (1/4) ad_X Σ_k ad_Y^k X / k!
for n in (0, 1, 2, 4):
    for k in range(0, 5):
        if n + k != 3:
            continue
        h1 = sc(comm(X, adn(Y, k, X)), R(1, 4) * R(1, sp.factorial(k)))
        el = sc(adn(Y, n, h1), F[n])
        if not el:
            continue
        cc = coords(el, cb)
        print('     [II] f_%d=%s  ad_Y^%d [X, ad_Y^%d X]/(4 %d!): c1..c6 = %s' % (n, F[n], n, k, k, [fmt(v) for v in cc]))
        tot23 = padd(tot23, el)
cc23 = coords(sc(tot23, R(1, 2)), cb)
print('     Z2 の (2,3) 部分 = (1/2)Σ: c1..c6 の係数 =', [fmt(v) for v in cc23])
print('     直接計算の (2,3) 部分      =', [fmt(v) for v in coords(bideg(C5d, 2, 3), cb)])

K5 = [R(7, 5760), cc32[1], cc32[2], cc23[3], cc23[4], R(1, 720)]
print('  補題 S の 5 次: C5 = Σ κ_i c_i、κ =', [fmt(v) for v in K5])
print('  直接計算との一致:', [sp.simplify(a - b) for a, b in zip(K5, k_direct)] == [0] * 6)
print('  c1 = [X,[X,V]]、c2 = [Y,[X,V]]、c3 = [U,V]、c4 = [X,[Y,W]]、c5 = [U,W]、c6 = [Y,[Y,W]] (U=[X,Y]、V=[X,U]、W=[Y,U])')


def sym(Xe, Ye):
    """補題 S: Xe + Ye + C3 + C5 (5 次まで。Xe, Ye は次数 1 以上の Lie 元)"""
    U_, V_, W_, c_ = parts(Xe, Ye)
    out = lsum(Xe, Ye, sc(V_, R(-1, 24)), sc(W_, R(-1, 12)))
    for k, b in zip(K5, c_):
        out = padd(out, sc(b, k))
    return out


# 補題 S の検算 (X, Y を一般の Lie 元にしても成立: 例 Xe = 2A + 3 [A,B]、Ye = B)
A, B = letters('A', 'B')
Xe = lsum(sc(A, 2), sc(comm(A, B), 3))
Ye = B
Sg = pmul(pmul(pexp(sc(Xe, R(1, 2))), pexp(Ye)), pmul(pexp(sc(Xe, R(1, 2))), {'': sp.Integer(1)}))
print('  (確認) Xe = 2A + 3[A,B]、Ye = B でも補題 S が成り立つ (5 次まで):', is_zero(padd(plog(Sg), sym(Xe, Ye), 1, -1)))

# =====================================================================
print('\n== (c) Yoshida (3.7) ==')
w1s = 1 / (2 - sp.cbrt(2))
w0s = 1 - 2 * w1s
Z1 = lsum(A, B)
Uab, Vab, Wab, cab = parts(A, B)
C3ab = lsum(sc(Vab, R(-1, 24)), sc(Wab, R(-1, 12)))
C5ab = {}
for k, b in zip(K5, cab):
    C5ab = padd(C5ab, sc(b, k))
print('  C3(A,B) = -V/24 - W/12、C5(A,B) = Σ κ_i c_i (補題 S で X=A, Y=B)')
a1, b1, a3, b3 = 2 * w1s, w0s, 2 * w1s**3, w0s**3
print('  a1 = 2w1、b1 = w0、a3 = 2w1^3、b3 = w0^3 とおく。X\' = 2Ω1 = a1 Z1 + a3 C3 + ...、Y\' = Ω0 = b1 Z1 + b3 C3 + ...')
print('  a1 b3 - a3 b1 = %s = -2 w1^3: %s' % (sp.nsimplify(sp.simplify(a1 * b3 - a3 * b1)), sp.simplify(a1 * b3 - a3 * b1 + 2 * w1s**3) == 0))
coef_c3 = sp.simplify((a1 * b3 - a3 * b1) * (-a1 / 24 - b1 / 12))
print('  [Z1,[Z1,C3]] の係数 (a1 b3 - a3 b1)(-a1/24 - b1/12) = %s = w1^3 (w1 + w0)/6: %s' % (
    sp.N(coef_c3, 10), sp.simplify(coef_c3 - w1s**3 * (w1s + w0s) / 6) == 0))
c5coef = sp.simplify(2 * w1s**5 + w0s**5)
print('  2 w1^5 + w0^5 = %s (= %.10f)' % (sp.nsimplify(c5coef), float(c5coef)))
# 検証: 自由代数で S4 の 5 次を直接
S4 = pmul(pmul(pexp(sc(Z1, 1)), {'': sp.Integer(1)}), {'': sp.Integer(1)})


def S2w(w):
    return pmul(pmul(pexp(sc(A, w / 2)), pexp(sc(B, w))), pexp(sc(A, w / 2)))


w1n = sp.Float(1, 40) / (2 - sp.Float(2, 40)**R(1, 3))
w0n = 1 - 2 * w1n
S4 = pmul(pmul(S2w(w1n), S2w(w0n)), S2w(w1n))
L4 = plog(S4)
Z1ZC3 = comm(Z1, comm(Z1, C3ab))
pred = padd(sc(C5ab, 2 * w1n**5 + w0n**5), sc(Z1ZC3, w1n**3 * (w1n + w0n) / 6))
dd = padd(deg(L4, 5), pred, 1, -1)
print('  S4 の 5 次 (自由代数の直接計算) - [(2w1^5+w0^5) C5 + w1^3(w1+w0)/6 [Z1,[Z1,C3]]] の最大の係数: %.1e' % max([abs(float(sp.N(v, 30))) for v in dd.values()] + [0.0]))

# 関数表現 (1 自由度、A = -p^2/2、B = -V)
x, p = sp.symbols('x p', real=True)
Vf = sp.Function('V')(x)


def PB(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x))


FUN = {'A': -p**2 / 2, 'B': -Vf}


def fun_of(builder):
    """builder(br, a, b) で Lie 元を関数表現に直す (括弧 [f,g] → PB(f,g))"""
    return sp.expand(builder(PB, FUN['A'], FUN['B']))


def basis_fun(Xf, Yf):
    Uf = PB(Xf, Yf)
    Vf_ = PB(Xf, Uf)
    Wf_ = PB(Yf, Uf)
    return [Uf, Vf_, Wf_,
            [PB(Xf, PB(Xf, Vf_)), PB(Yf, PB(Xf, Vf_)), PB(Uf, Vf_), PB(Xf, PB(Yf, Wf_)), PB(Uf, Wf_), PB(Yf, PB(Yf, Wf_))]]


Uf, Vff, Wff, cf = basis_fun(FUN['A'], FUN['B'])
d1, d2, d3, d4 = [sp.diff(Vf, x, n) for n in (1, 2, 3, 4)]
print('  関数表現 (A=-p^2/2、B=-V):')
print('    U = [A,B] ↦', sp.factor(Uf), '、V ↦', sp.factor(Vff), '、W ↦', sp.factor(Wff))
for i, ci in enumerate(cf):
    print('    c%d ↦ %s' % (i + 1, sp.collect(sp.expand(ci), p)))


def pcoeffs(f):
    P = sp.Poly(sp.expand(f), p)
    return {m[0]: sp.simplify(c) for m, c in zip(P.monoms(), P.coeffs())}


C5f = sum(k * c for k, c in zip(K5, cf))
C3f = -Vff / 24 - Wff / 12
Z1f = FUN['A'] + FUN['B']
ZZC3f = PB(Z1f, PB(Z1f, C3f))
print('  Strang の C3 の関数 = -V/24 - W/12 ↦', sp.collect(sp.expand(C3f), p))
print('  Strang の C5 の関数 (p のべきごと):', {k: v for k, v in pcoeffs(C5f).items()})
print('  [Z1,[Z1,C3]] の関数 (p のべきごと):', {k: v for k, v in pcoeffs(ZZC3f).items()})
H5y = sp.expand(c5coef * C5f + w1s**3 * (w1s + w0s) / 6 * ZZC3f)
pc = pcoeffs(H5y)
print('  s5 = 2 w1^5 + w0^5 = %.10f、r = w1^3 (w1 + w0)/6 = %.10f' % (float(c5coef), float(w1s**3 * (w1s + w0s) / 6)))
print('  H5(Yoshida) = s5 * (Strang の C5 の関数) + r * ([Z1,[Z1,C3]] の関数):')
# 係数ごとに取り出して数値表示
def num_coeffs(expr):
    out = {}
    e = sp.expand(expr)
    for pw in (0, 2, 4):
        cp = e.coeff(p, pw)
        out[pw] = sp.collect(sp.expand(cp), [d1 * d1 * d2, d1 * d3, d2 * d2, d4])
    return out
for pw, cp in num_coeffs(H5y).items():
    terms = sp.Poly(cp, d1, d2, d3, d4)
    print('    p^%d:' % pw, ', '.join('%s * %s' % (sp.N(cf_, 30).round(6), sp.factor(sp.prod([g**e for g, e in zip((d1, d2, d3, d4), mon)]))) for mon, cf_ in zip(terms.monoms(), terms.coeffs())))
print('  厳密な値 (w1 = 1/(2 - 2^(1/3))): p^4 の V\'\'\'\' の係数 =', sp.nsimplify(sp.simplify(pcoeffs(H5y)[4] / d4)), '=', sp.N(pcoeffs(H5y)[4] / d4, 8))

# =====================================================================
print('\n== (d) 力勾配型 (3.9) の 5 次を 3 段で組む ==')
t, c3s, c5s = sp.symbols('t c3 c5')
aa = R(1, 2) - t
Uab, Vab, Wab, cab = parts(A, B)                 # 以降の基底は X→A、Y→B
Wb = comm(B, comm(A, B))                         # = [B,U]
U_ = comm(A, B)
W5b = comm(Wb, U_)
Cel = padd(sc(Wb, c3s), sc(W5b, c5s))
print('  W = [B,[A,B]] = [B,U]、W5 = [W,U]、C = c3 W + c5 W5')
# 内側: e^{tA} e^{C} e^{tA} = e^{Y1}、Y1 = sym(2tA, C)
Y1 = sym(sc(A, 2 * t), Cel)
Y1_dir = plog(pmul(pmul(pexp(sc(A, t)), pexp(Cel)), pexp(sc(A, t))))
print('  内側 Y1 = 2tA + C - (t^2/6)[A,[A,C]] + O(7): 補題 S と直接計算が一致:', is_zero(padd(Y1_dir, Y1, 1, -1)))
print('    Y1 の 5 次 = c5 W5 - (t^2 c3/6)[A,[A,W]] (他に 5 次なし):',
      is_zero(padd(deg(Y1, 5), padd(sc(W5b, c5s), sc(comm(A, comm(A, Wb)), -t**2 * c3s / 6)), 1, -1)))
# 中間: e^{B/2} e^{Y1} e^{B/2} = e^{Y2}、Y2 = sym(B, Y1)
Y2 = sym(B, Y1)
Y2_dir = plog(pmul(pmul(pexp(sc(B, R(1, 2))), pexp(Y1)), pexp(sc(B, R(1, 2)))))
print('  中間 Y2: 補題 S と直接計算が一致:', is_zero(padd(Y2_dir, Y2, 1, -1)))
# 外側: e^{aA} e^{Y2} e^{aA}
Z = sym(sc(A, 2 * aa), Y2)
Zdir = plog(pmul(pmul(pexp(sc(A, aa)), pexp(Y2)), pexp(sc(A, aa))))
print('  外側 Z = log S_FG: 補題 S と直接計算が一致:', is_zero(padd(Zdir, Z, 1, -1)))
# Y2, Z の 3 次、5 次を基底で
def coords_ab(L, n):
    if n == 3:
        return coords(L, [comm(A, comm(A, B)), comm(B, comm(B, A))])
    return coords(L, cab)
print('  3 次: Y2 の係数 ([A,[A,B]], [B,[B,A]]) =', [fmt(v) for v in coords_ab(deg(Y2, 3), 3)])
print('        Z の係数 =', [fmt(sp.simplify(v)) for v in coords_ab(deg(Z, 3), 3)], '(本文 (3.11): (12t^2-1)/24、1/12 - t/4 - c3)')
# 5 次の寄与の内訳 (外側)
print('  5 次 (外側 Z) の内訳: Z5 = [Y2 の 5 次] + [C3(2aA, Y2) の 5 次の部分] + [C5(2aA, B + 2tA)]')
Y2_3 = deg(Y2, 3)
Y2_1 = lsum(B, sc(A, 2 * t))
Xo = sc(A, 2 * aa)
Uo, Vo, Wo, co = parts(Xo, Y2_1)
c5_lead = {}
for k, b in zip(K5, co):
    c5_lead = padd(c5_lead, sc(b, k))
UoV, VoV, WoV, _ = parts(Xo, Y2_3)
# C3(X, Y2) の次数 5: Y2 の 3 次部分を 1 つ含む項
C3_cross = padd(sc(comm(Xo, comm(Xo, Y2_3)), R(-1, 24)),
                sc(padd(comm(Y2_3, comm(Y2_1, Xo)), comm(Y2_1, comm(Y2_3, Xo))), R(1, 12)))
parts_list = [('Y2 の 5 次そのもの', deg(Y2, 5)), ('C3 の交差項 (Y2 の 3 次を含む)', C3_cross), ('C5(2aA, B + 2tA)', c5_lead)]
tot = {}
for name, el in parts_list:
    print('    %s: c1..c6 = %s' % (name, [fmt(sp.simplify(v)) for v in coords(el, cab)]))
    tot = padd(tot, el)
tv0 = sp.sqrt(3) / 6
c3v0 = R(1, 12) - tv0 / 4
print('    (t = 1/sqrt(12)、c3 = 1/12 - t/4 を代入した値)')
for name, el in parts_list:
    print('    %s: b1..b6 = %s' % (name, [sp.expand(sp.simplify(v.subs({t: tv0, c3s: c3v0}))) for v in coords(el, cab)]))
print('    合計が Z の 5 次と一致:', is_zero(padd(tot, deg(Z, 5), 1, -1)))
Z5c = [sp.simplify(v) for v in coords(deg(Z, 5), cab)]
print('  Z の 5 次: b1..b6 の係数 (b1=[A,[A,V]]、b2=[B,[A,V]]、b3=[U,V]、b4=[A,[B,W]]、b5=[U,W]、b6=[B,[B,W]]、ここで U=[A,B]、V=[A,U]、W=[B,U])')
print('   ※ 補題 S の基底の名前 (X=A, Y=B): c1=[A,[A,V]]、c2=[B,[A,V]]、c3=[U,V]、c4=[A,[B,W]]、c5=[U,W]、c6=[B,[B,W]]')
for i, v in enumerate(Z5c):
    print('    b%d: %s' % (i + 1, sp.factor(v)))
tv = sp.sqrt(3) / 6
c3v = R(1, 12) - tv / 4
Z5v = [sp.simplify(v.subs({t: tv, c3s: c3v})) for v in Z5c]
print('  t = 1/sqrt(12)、c3 = 1/12 - t/4 を代入:')
for i, v in enumerate(Z5v):
    print('    b%d: %s (= %.9g)' % (i + 1, sp.nsimplify(v), float(v.subs(c5s, 0)) if True else 0))
# 関数表現
bf = basis_fun(FUN['A'], FUN['B'])[3]
H5f = sp.expand(sum(v * f for v, f in zip(Z5v, bf)))
print('  関数表現 (b_i ↦):')
for i, ci in enumerate(bf):
    print('    b%d ↦ %s' % (i + 1, sp.collect(sp.expand(ci), p)))
print('  H5(FG) (p のべきごと):')
for pw in (0, 2, 4):
    cp = sp.expand(H5f).coeff(p, pw)
    print('    p^%d: %s' % (pw, sp.simplify(cp)))
ref0 = (-720 * c5s + 9 - 5 * sp.sqrt(3)) / 360 * d1**2 * d2
ref2 = (sp.sqrt(3) / 288 - R(1, 120)) * d1 * d3 + (R(1, 240) - sp.sqrt(3) / 288) * d2**2
ref4 = d4 / 4320
print('  (3.13) との差: p^0 %s、p^2 %s、p^4 %s' % (
    sp.simplify(sp.expand(H5f).coeff(p, 0) - ref0), sp.simplify(sp.expand(H5f).coeff(p, 2) - ref2), sp.simplify(sp.expand(H5f).coeff(p, 4) - ref4)))

# =====================================================================
print('\n== (e) W5 についての主張 ==')
Ub = comm(A, B)
Vb = comm(A, Ub)
Wbb = comm(B, Ub)
lhs = comm(B, comm(B, Vb))
rhs = padd(comm(A, comm(B, Wbb)), comm(Ub, Wbb), 1, -1)
print('  [B,[B,[A,[A,B]]]] = [A,[B,W]] - [U,W] (J3):', is_zero(padd(lhs, rhs, 1, -1)))
fU = PB(FUN['A'], FUN['B'])
fW = PB(FUN['B'], fU)
print('  関数表現: B = -V、W ↦ %s はどちらも x だけの関数で [B,W] ↦ %s' % (sp.factor(fW), PB(FUN['B'], fW)))
fW5 = PB(fW, fU)
fl = PB(FUN['B'], PB(FUN['B'], PB(FUN['A'], PB(FUN['A'], FUN['B']))))
print('  [B,[B,[A,[A,B]]]] ↦ %s、W5 ↦ %s、比 = %s' % (sp.factor(fl), sp.factor(fW5), sp.simplify(fl / fW5)))
# 多自由度
x1, x2, p1, p2 = sp.symbols('x1 x2 p1 p2', real=True)
Vm = sp.Function('V')(x1, x2)
xs, ps = (x1, x2), (p1, p2)


def PB2(f, g):
    return sp.expand(sum(sp.diff(f, xs[j]) * sp.diff(g, ps[j]) - sp.diff(f, ps[j]) * sp.diff(g, xs[j]) for j in range(2)))


T2 = (p1**2 + p2**2) / 2
fAB = PB2(T2, Vm)
fW2 = -PB2(Vm, fAB)
fW52 = PB2(fW2, fAB)
Vj = [sp.diff(Vm, xj) for xj in xs]
Tm = sp.symbols('m', positive=True)
print('  2 自由度での b1..b6 の関数 (T = (p1^2+p2^2)/2、A = -:T:、B = -:V:):')
fAm, fBm = -T2, -Vm
fU2 = PB2(fAm, fBm)
fV2 = PB2(fAm, fU2)
fW2b = PB2(fBm, fU2)
fb = [PB2(fAm, PB2(fAm, fV2)), PB2(fBm, PB2(fAm, fV2)), PB2(fU2, fV2), PB2(fAm, PB2(fBm, fW2b)), PB2(fU2, fW2b), PB2(fBm, PB2(fBm, fW2b))]
for i, fi in enumerate(fb):
    P_ = sp.Poly(sp.expand(fi), p1, p2) if fi != 0 else None
    degs = sorted({sum(m_) for m_ in P_.monoms()}) if P_ is not None else []
    print('    b%d ↦ 運動量の次数 %s%s' % (i + 1, degs if degs else 'なし (0)', '' if fi != 0 else ''))
print('    b5 ↦ 2 Σ V_j V_k V_jk との差:', sp.simplify(fb[4] - 2 * sum(Vj[j] * Vj[k] * sp.diff(Vm, xs[j], xs[k]) for j in range(2) for k in range(2))))
print('  2 自由度: [T,V] =', sp.simplify(fAB), '(= -Σ p_j V_j)')
print('           W ↦', sp.simplify(fW2), '(= Σ V_j^2)')
Gx = sum(Vj[j] * Vj[k] * sp.diff(Vm, xs[j], xs[k]) for j in range(2) for k in range(2))
print('           W5 ↦ (-2 Σ V_j V_k V_jk との差):', sp.simplify(fW52 + 2 * Gx))

# =====================================================================
print('\n== (f) §8.4: A を厳密なドリフトにしたときの W (1 自由度) ==')
pz = sp.symbols('pz', positive=True)
pp = sp.symbols('P', positive=True)
Tex = -sp.sqrt(pp**2 - p**2)   # 厳密なドリフト H_a = -p_z (p_x = p、全運動量 P。P は定数として扱う)
fABe = sp.simplify(PB(Tex, Vf))
fWe = sp.simplify(-PB(Vf, fABe))
print('  H_a = -sqrt(P^2 - p^2): [H_a, V] =', fABe, '、W = -[V, [H_a,V]] =', fWe)
print('  展開: W = V\'^2 / p_z、p_z = sqrt(P^2 - p^2)。1/p_z = (1/P)(1 + p^2/(2P^2) + ...):', sp.series(1 / sp.sqrt(pp**2 - p**2), p, 0, 4))
