# -*- coding: utf-8 -*-
r"""07_thin_multipole.md の途中式の補完 (sympy)。

(A) 補正キックの母関数 f_C (式 (3.15)) から Δq_x, Δq_y, Δz を一般の k = 1..5 で作り、(7.12) の係数 f1..f5 の式と比べる。
    K64 の tthin の 3 つの分岐 (kord == 2、kord == 1、一般) を書き写して同じ式になることも確かめる。
    途中の恒等式 (7.10) (7.11)、∂_x rcx1 = (k+1) Re cx、∂_x r^k なども k = 1..5 で確かめる。
(B) リー代数の表現: A = -:T:、B = -:V: とし、[A,B]、W = [B,[A,B]]、W5 = [W,[A,B]] が |∇V|^2/m、-2 V_j V_k V_jk/m^2 になること。
(C) pysad の追跡 (DECA、DODECA を含む k = 2..5) と、(7.12) の一般 k の式の比較。
(D) 端部の非線形フリンジ: (7.16) のラプラス方程式、(7.17) の導出の中間式、K64 の cp1 の式と J^{-T} p の一致、6x6 のシンプレクティック性。
(E) フィードダウンの一般 k の式 (7.14) の展開、CALC のキックのヤコビアンと 4x4 行列の成分。
(F) 6 極の追跡 (小さい傾き) の厳密な積分との差の符号 (成分ごと)。
"""
import io
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

x, y, p, kap, L = sp.symbols('x y p kappa L', positive=True)
beta, gamma = sp.symbols('beta gamma', positive=True)
c3, c5 = beta / 2, gamma / 2
zb = x - sp.I * y


def Vk(k):
    return sp.expand(kap / sp.factorial(k + 1) * sp.re(sp.expand(zb**(k + 1))))


print('== (A) f_C から作った Δq_x, Δq_y, Δz と (7.12) の比較 (一般の k) ==')
print('   f_C = c3 L^3 |∇V|^2/p - 2 c5 L^5 V_j V_k V_jk / p^2、c3 = β/2、c5 = γ/2、Δq = (1/p) ∂f_C/∂x、Δz = -∂f_C/∂p')
for k in range(1, 6):
    V = Vk(k)
    Vx, Vy = sp.diff(V, x), sp.diff(V, y)
    Vxx, Vxy, Vyy = sp.diff(V, x, 2), sp.diff(V, x, y), sp.diff(V, y, 2)
    gradsq = sp.expand(Vx**2 + Vy**2)
    G = sp.expand(Vx * Vx * Vxx + 2 * Vx * Vy * Vxy + Vy * Vy * Vyy)      # V_j V_k V_jk
    r = x**2 + y**2
    cx = sp.expand(zb**k)
    Rcx, Icx = sp.re(cx), sp.im(cx)
    rcx1 = sp.expand(x * Rcx + y * Icx)
    a0 = kap / sp.factorial(k)
    d710 = sp.simplify(sp.expand(gradsq - a0**2 * r**k))
    d711 = sp.simplify(sp.expand(G - k * a0**3 * r**(k - 1) * rcx1))
    d_rcx = sp.simplify(sp.expand(sp.diff(rcx1, x) - (k + 1) * Rcx))
    d_rcy = sp.simplify(sp.expand(sp.diff(rcx1, y) - (k + 1) * Icx))
    fC = c3 * L**3 * gradsq / p - 2 * c5 * L**5 * G / p**2
    dqx = sp.diff(fC, x) / p
    dqy = sp.diff(fC, y) / p
    dz = -sp.diff(fC, p)
    aki = kap * L / (sp.factorial(k) * p)
    f1 = beta * k * L
    f2 = 2 * gamma * k * (k - 1) * L**2
    f3 = gamma * k * (k + 1) * L**2
    f4 = beta * L / 2
    f5 = 2 * gamma * k * L**2
    rk1 = r**(k - 2) if k >= 2 else 0
    rk = r**(k - 1)
    t712x = aki**2 * (f1 * rk * x - aki * (f2 * rk1 * x * rcx1 + f3 * rk * Rcx))
    t712y = aki**2 * (f1 * rk * y - aki * (f2 * rk1 * y * rcx1 + f3 * rk * Icx))
    t712z = aki**2 * rk * (f4 * r - f5 * aki * rcx1)
    e = [sp.simplify(sp.expand(dqx - t712x)), sp.simplify(sp.expand(dqy - t712y)), sp.simplify(sp.expand(dz - t712z))]
    print('  k = %d: (7.10) 差 %s、(7.11) 差 %s、∂_x rcx1 - (k+1)Re cx = %s、∂_y rcx1 - (k+1)Im cx = %s' % (k, d710, d711, d_rcx, d_rcy))
    print('         (7.12) との差 Δq_x: %s、Δq_y: %s、Δz: %s' % tuple(e))
    al = L
    if k == 2:
        F = (beta * 2 * al, 4 * gamma * al**2, 6 * gamma * al**2, sp.Rational(1, 2) * beta * al, 4 * gamma * al**2)
        kx = aki**2 * (F[0] * r * x - aki * (F[1] * x * rcx1 + F[2] * r * Rcx))
        ky = aki**2 * (F[0] * r * y - aki * (F[1] * y * rcx1 + F[2] * r * Icx))
        kz = aki**2 * r * (F[3] * r - F[4] * aki * rcx1)
        br = 'kord == 2 の分岐'
    elif k == 1:
        F = (beta * al, 0, 2 * gamma * al**2, beta * al / 2, 2 * gamma * al**2)
        kx = aki**2 * (F[0] - aki * F[2]) * x
        ky = aki**2 * (F[0] + aki * F[2]) * y
        kz = aki**2 * (F[3] * r - F[4] * aki * rcx1)
        br = 'kord == 1 の分岐 (aki はソースでこの分岐の前に設定されない。ここでは aki = kL/(k! p) と仮定)'
    else:
        F = (beta * k * al, 2 * gamma * k * (k - 1) * al**2, gamma * k * (k + 1) * al**2, beta * al / 2, 2 * gamma * k * al**2)
        rk1_ = r**(k - 2)
        rk_ = r * rk1_
        kx = aki**2 * (F[0] * rk_ * x - aki * (F[1] * rk1_ * x * rcx1 + F[2] * rk_ * Rcx))
        ky = aki**2 * (F[0] * rk_ * y - aki * (F[1] * rk1_ * y * rcx1 + F[2] * rk_ * Icx))
        kz = aki**2 * rk_ * (F[3] * r - F[4] * aki * rcx1)
        br = '一般の分岐 (kord >= 3)'
    e2 = [sp.simplify(sp.expand(dqx - kx)), sp.simplify(sp.expand(dqy - ky)), sp.simplify(sp.expand(dz - kz))]
    print('         K64 の %s との差: %s、%s、%s' % (br, *e2))

print('\n   f_C の偏微分を 1 つずつ (k = 3 の例):')
k = 3
V = Vk(k)
gradsq = sp.expand(sp.diff(V, x)**2 + sp.diff(V, y)**2)
a0 = kap / sp.factorial(k)
print('   |∇V|^2 = %s' % sp.factor(gradsq))
print('   ∂_x |∇V|^2 - a0^2 2k r^(k-1) x = %s' % sp.simplify(sp.diff(gradsq, x) - a0**2 * 2 * k * (x**2 + y**2)**(k - 1) * x))

print('\n== (B) リー代数の表現 (C 付録の規約 [f,g] = Σ(∂q f ∂p g - ∂p f ∂q g)、:f:g = [f,g]) ==')
X, Y, PX, PY, m = sp.symbols('X Y PX PY m', real=True)
Vg = sp.Function('V')(X, Y)


def PB(f, g):
    return sum(sp.diff(f, q) * sp.diff(g, pp) - sp.diff(f, pp) * sp.diff(g, q) for q, pp in ((X, PX), (Y, PY)))


T = (PX**2 + PY**2) / (2 * m)
fA, fB = -T, -Vg                    # A = :fA:、B = :fB: (exp(:f:) の f = -∫H)
fAB = sp.simplify(PB(fA, fB))       # [:f:, :g:] = :[f, g]: なので、交換子は括弧の関数で表す
fW = sp.simplify(PB(fB, fAB))
fW5 = sp.simplify(PB(fW, fAB))
print('  [A,B] <-> [fA, fB] = %s' % fAB)
print('  W = [B,[A,B]] <-> %s' % sp.simplify(fW))
print('  W5 = [W,[A,B]] <-> %s' % sp.simplify(fW5))
Vx_, Vy_ = sp.diff(Vg, X), sp.diff(Vg, Y)
VVV = Vx_ * Vx_ * sp.diff(Vg, X, 2) + 2 * Vx_ * Vy_ * sp.diff(Vg, X, Y) + Vy_ * Vy_ * sp.diff(Vg, Y, 2)
print('  W - |∇V|^2/m = %s、W5 - (-2 V_j V_k V_jk/m^2) = %s' % (sp.simplify(fW - (Vx_**2 + Vy_**2) / m), sp.simplify(fW5 + 2 * VVV / m**2)))
print('  V_j V_k V_jk - (1/2) ∇V·∇|∇V|^2 = %s' % sp.simplify(VVV - sp.Rational(1, 2) * (Vx_ * sp.diff(Vx_**2 + Vy_**2, X) + Vy_ * sp.diff(Vx_**2 + Vy_**2, Y))))

print('\n== (C) pysad の追跡 (SEXT、OCT、DECA、DODECA) と (7.12) の一般 k の式 ==')
from pysad.ffscmd import run_sad_string
ALPHA = 1.0 / math.sqrt(12.0)
ALPHA1 = 0.5 - ALPHA
BETA = 1.0 / 6.0 - ALPHA / 2.0
GAMMA = 1.0 / 40.0 - 1.0 / (24.0 * math.sqrt(3.0))


def drift(c, h):
    c = c.copy()
    pz = math.sqrt(1.0 - c[1]**2 - c[3]**2)
    a = c[1]**2 + c[3]**2
    c[0] += c[1] / pz * h
    c[2] += c[3] / pz * h
    c[4] -= a / (1.0 + pz) / pz * h
    return c


def kick(c, akt, k):
    c = c.copy()
    g = akt / math.factorial(k) * complex(c[0], -c[2])**k
    c[1] -= g.real
    c[3] -= g.imag
    return c


def tthin_gen(c, Lh, ak, k, corr=True):
    """(7.12) の一般 k の式 (k >= 2)。delta = 0、p = 1。corr=False で補正キックなし"""
    c = drift(c, Lh * ALPHA1)
    c = kick(c, ak * 0.5, k)
    c = drift(c, Lh * ALPHA)
    xx, yy = c[0], c[2]
    r = xx * xx + yy * yy
    cx = complex(xx, -yy)**k
    rcx1 = xx * cx.real + yy * cx.imag
    aki = ak / math.factorial(k)
    f1, f2, f3 = BETA * k * Lh, 2 * GAMMA * k * (k - 1) * Lh**2, GAMMA * k * (k + 1) * Lh**2
    f4, f5 = 0.5 * BETA * Lh, 2 * GAMMA * k * Lh**2
    rk1 = r**(k - 2)
    rk = r * rk1
    c = c.copy()
    if corr:
        c[1] += aki**2 * (f1 * rk * xx - aki * (f2 * rk1 * xx * rcx1 + f3 * rk * cx.real))
        c[3] += aki**2 * (f1 * rk * yy - aki * (f2 * rk1 * yy * rcx1 + f3 * rk * cx.imag))
        c[4] += aki**2 * rk * (f4 * r - f5 * aki * rcx1)
    c = drift(c, Lh * ALPHA)
    c = kick(c, ak * 0.5, k)
    return drift(c, Lh * ALPHA1)


NAMES = {2: ('SEXT', 'K2'), 3: ('OCT', 'K3'), 4: ('DECA', 'K4'), 5: ('DODECA', 'K5')}
KAPPA = {2: 400.0, 3: 4.0e4, 4: 4.0e6, 5: 4.0e8}
c0 = np.array([0.020, 0.0010, -0.010, -0.0005, 0.0, 0.0])
for k in (2, 3, 4, 5):
    name, key = NAMES[k]
    worst = 0.0
    for Lh in (0.1, 0.4):
        akt = KAPPA[k] * Lh
        src = ('MOMENTUM = 3GEV;\n%s MP=(L=%.15g %s=%.15g DISFRIN=1);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n' % (name, Lh, key, akt))
        ffs, _ = run_sad_string(src, out=io.StringIO())
        tr = C.track_one(ffs, c0)
        mine = tthin_gen(c0, Lh, akt, k)
        nocorr = tthin_gen(c0, Lh, akt, k, corr=False)
        worst = max(worst, np.abs(tr - mine).max())
    print('  %s (k = %d、kappa = %g、L = 0.4 まで): pysad の追跡と (7.12) の式の最大差 %.2e (x, px, y, py, z)。補正キックを抜くと差は %.1e' % (
        name, k, KAPPA[k], worst, np.abs(tr - nocorr).max()))

print('\n== (D) 端部の非線形フリンジ ==')
w, wb = sp.symbols('w wb')
for m_ in (2, 3, 4, 5, 6):
    G = w * wb * (w**m_ - wb**m_) / (2 * sp.I)               # |w|^2 Im w^m
    lap = sp.expand(4 * sp.diff(sp.diff(G, w), wb))
    chk1 = sp.simplify(lap - 4 * (m_ + 1) * (w**m_ - wb**m_) / (2 * sp.I))
    d2 = sp.expand(2 * sp.diff(G, wb))
    chk2 = sp.simplify(d2 - (-sp.I) * w * (w**m_ - (m_ + 1) * wb**m_))
    print('  m = %d: ∇⊥^2 (|w|^2 Im w^m) - 4(m+1) Im w^m = %s、2∂_w̄ (|w|^2 Im w^m) - (-i) w [w^m - (m+1) w̄^m] = %s' % (m_, chk1, chk2))
s = sp.symbols('s', real=True)
kf = sp.Function('kappa')(s)
for m_ in (2, 3, 4):
    xx_, yy_ = sp.symbols('xx yy', real=True)
    wv = xx_ + sp.I * yy_
    imw = sp.im(sp.expand(wv**m_))
    Phi = kf * imw / sp.factorial(m_) - sp.diff(kf, s, 2) * (xx_**2 + yy_**2) * imw / (4 * sp.factorial(m_ + 1))
    lap3 = sp.expand(sp.diff(Phi, xx_, 2) + sp.diff(Phi, yy_, 2) + sp.diff(Phi, s, 2))
    print('  m = %d: ∇^2 Φ の残り (κ の 4 階微分の項だけ残る) = %s' % (m_, sp.factor(lap3)))
for m_ in (2, 3, 4):
    lhs = -w * (w**m_ - (m_ + 1) * wb**m_) / (4 * sp.factorial(m_ + 1))
    rhs = -w * (w**m_ / (m_ + 1) - wb**m_) / (4 * sp.factorial(m_))
    print('  m = %d: (7.17) の 2 つの形の差 = %s' % (m_, sp.simplify(lhs - rhs)))
xs, ys = sp.symbols('xs ys', real=True)
wq = xs + sp.I * ys
dwq = -(1 / sp.Integer(8)) * wq * (wq**2 / 3 - sp.conjugate(wq)**2)
print('  m = 2: Re Δw (kappa0/p = 1) = %s' % sp.factor(sp.re(sp.expand(dwq))))
print('         Im Δw = %s' % sp.factor(sp.im(sp.expand(dwq))))

a_, xv, yv, pxv, pyv = sp.symbols('aki xv yv pxv pyv', real=True)
for m_ in (3, 4, 5):
    wv = xv + sp.I * yv
    dw = -a_ * wv * (wv**m_ / (m_ + 1) - sp.conjugate(wv)**m_)
    xn = xv + sp.re(sp.expand(dw))
    yn = yv + sp.im(sp.expand(dw))
    J = sp.Matrix([[sp.diff(xn, xv), sp.diff(xn, yv)], [sp.diff(yn, xv), sp.diff(yn, yv)]])
    detJ = sp.simplify(J.det())
    Pn = sp.simplify(J.T.inv() * sp.Matrix([pxv, pyv]))
    cp = pxv - sp.I * pyv
    aa = 2 * a_ * sp.im(sp.expand(wv**m_))
    d = 1 + aa**2 - m_**2 * a_**2 * (xv**2 + yv**2)**m_
    cz1 = wv**(m_ - 1)
    cp1 = ((1 + sp.I * aa) * cp - m_ * a_ * cz1 * sp.conjugate(wv * cp)) / d
    cpn = Pn[0] - sp.I * Pn[1]
    print('  m = %d: det J - d = %s、(J^{-T} p) の px - i py と K64 の cp1 の差 = %s' % (
        m_, sp.simplify(sp.expand(detJ - d)), sp.simplify(sp.expand(cpn - cp1))))
# z と 6x6 のシンプレクティック性
m_ = 3
X6 = sp.symbols('x y Px Py z dl', real=True)
xv, yv, Pxv, Pyv, zv, dl = X6
k0 = sp.Rational(7, 3)
aki6 = k0 / (1 + dl)
wv = xv + sp.I * yv
dw = -aki6 * wv * (wv**m_ / (m_ + 1) - sp.conjugate(wv)**m_)
xn = xv + sp.re(sp.expand(dw))
yn = yv + sp.im(sp.expand(dw))
J = sp.Matrix([[sp.diff(xn, xv), sp.diff(xn, yv)], [sp.diff(yn, xv), sp.diff(yn, yv)]])
Pn = J.T.inv() * sp.Matrix([Pxv, Pyv])
for sign, label in ((-1, "z' = z - dr.P'/(1+dl)  [(7.18)]"), (+1, "z' = z + dr.P'/(1+dl)")):
    zn = zv + sign * (sp.re(sp.expand(dw)) * Pn[0] + sp.im(sp.expand(dw)) * Pn[1]) / (1 + dl)
    newv = sp.Matrix([xn, Pn[0], yn, Pn[1], zn, dl])
    old = sp.Matrix([xv, Pxv, yv, Pyv, zv, dl])
    Jac = newv.jacobian(old)
    pt = {xv: sp.Rational(1, 100), yv: -sp.Rational(3, 200), Pxv: sp.Rational(1, 1000), Pyv: -sp.Rational(1, 2000), zv: 0, dl: sp.Rational(1, 200)}
    Jn = np.array(Jac.subs(pt).evalf(30).tolist(), dtype=float)
    Om = np.zeros((6, 6))
    for i in range(3):
        Om[2 * i, 2 * i + 1] = 1.0
        Om[2 * i + 1, 2 * i] = -1.0
    print('  m = 3 (係数 7/3、delta = 0.005): %s のとき max|M^T Om M - Om| = %.2e' % (label, np.abs(Jn.T @ Om @ Jn - Om).max()))

print('\n== (E) フィードダウンと CALC のキックのヤコビアン ==')
dxs, dys, K = sp.symbols('d_x d_y K', real=True)
zbar = x - sp.I * y
dbar = dxs - sp.I * dys
for k in (2, 3, 4):
    exact_ = sp.expand(-K / sp.factorial(k) * (zbar - dbar)**k)
    ser = 0
    for j in range(k + 1):
        Kj = K * (-dbar)**j / sp.factorial(j)
        ser += -Kj / sp.factorial(k - j) * zbar**(k - j)
    print('  k = %d: (7.14) の二項展開の差 = %s' % (k, sp.simplify(sp.expand(exact_ - ser))))
    for j in range(0, k + 1):
        Kj = sp.expand(K * (-dbar)**j / sp.factorial(j))
        print('      j = %d: 次数 %d の多極、強さ K_%d = (%s) + i(%s)' % (j, k - j, k - j, sp.expand(sp.re(Kj)), sp.expand(sp.im(Kj))))
for k in (2, 3):
    Pc = -K / sp.factorial(k) * zbar**k
    Dpx, Dpy = sp.re(sp.expand(Pc)), sp.im(sp.expand(Pc))
    gz = sp.expand(k * zbar**(k - 1))
    gr, gi = sp.re(gz), sp.im(gz)
    ak_ = K / sp.factorial(k)
    chk = [sp.simplify(sp.diff(Dpx, x) + ak_ * gr), sp.simplify(sp.diff(Dpx, y) + ak_ * gi), sp.simplify(sp.diff(Dpy, x) + ak_ * gi), sp.simplify(sp.diff(Dpy, y) - ak_ * gr)]
    print('  k = %d: (7.13) の 4 成分の差 = %s' % (k, chk))
l1, l2, h1, h2 = sp.symbols('l1 l2 h1 h2')
Dm = lambda l: sp.Matrix([[1, l], [0, 1]])
Km = lambda h: sp.Matrix([[1, 0], [-h, 1]])
Mx = sp.expand(Dm(l1) * Km(h2) * Dm(l2) * Km(h1) * Dm(l1))
print('  水平面 (x, px) の 2x2 (y = 0、px = 0 の COD): M = D(l1) K(h2) D(l2) K(h1) D(l1)、K(h) = [[1,0],[-h,1]]、h = (ak/k!/2) k x^(k-1) (キックの位置の x):')
for (i_, j_) in ((0, 0), (0, 1), (1, 0), (1, 1)):
    print('     M[%d,%d] = %s' % (i_, j_, sp.collect(Mx[i_, j_], [h1, h2])))
Ls_ = sp.symbols('L_s')
Mc = sp.expand(Mx.subs({h2: h1}).subs({l1: Ls_ / 6, l2: 2 * Ls_ / 3}))
print('  h1 = h2 = h、l1 = L/6、l2 = 2L/3 のとき:')
for (i_, j_) in ((0, 0), (0, 1), (1, 0), (1, 1)):
    print('     M[%d,%d] = %s' % (i_, j_, sp.expand(Mc[i_, j_].subs(h1, sp.Symbol('h')))))
from pysad import elkern
for k in (2, 3):
    xc, Lh = 0.004, 0.1
    akt = KAPPA[k] * Lh
    c = np.array([xc, 0.0, 0.0, 0.0, 0.0, 0.0])
    ok, Jc = elkern._jac_thin(c.copy(), Lh, akt, k, 0.0)
    ah = akt / math.factorial(k) * 0.5
    hh1 = ah * k * xc**(k - 1)
    px1 = -ah * xc**k
    x2 = xc + Lh / 1.5 * px1                      # 1 つめのキックのあと、ドリフト 2L/3 を進んだ位置
    hh2 = ah * k * x2**(k - 1)
    sub = {l1: Lh / 6, l2: Lh / 1.5}
    Mn = np.array(Mx.subs(sub).subs({h1: hh1, h2: hh2}).tolist(), dtype=float)
    Mn_same = np.array(Mx.subs(sub).subs({h1: hh1, h2: hh1}).tolist(), dtype=float)
    Mv = np.array(Mx.subs(sub).subs({h1: -hh1, h2: -hh2}).tolist(), dtype=float)
    print('  k = %d、x = 4 mm、L = 0.1: elkern._jac_thin の (x, px) ブロックとの差: 2 つめのキックの位置 x2 を使う %.2e、(h2 = h1 とした場合 %.2e)、(y, py) ブロック (h -> -h) %.2e' % (
        k, np.abs(Jc[:2, :2] - Mn).max(), np.abs(Jc[:2, :2] - Mn_same).max(), np.abs(Jc[2:4, 2:4] - Mv).max()))

print('\n== (F) 6 極 (傾きが小さい初期値) の tthin の誤差の符号 (成分ごと、厳密な積分は DOP853) ==')
from scipy.integrate import solve_ivp


def exact(c, Lh, kappa, k):
    def rhs(s_, v):
        xx, px, yy, py, zz = v
        pz = math.sqrt(1.0 - px * px - py * py)
        g = kappa / math.factorial(k) * complex(xx, -yy)**k
        return [px / pz, -g.real, py / pz, -g.imag, 1.0 - 1.0 / pz]
    sol = solve_ivp(rhs, (0, Lh), list(c[:5]), method='DOP853', rtol=1e-13, atol=1e-18)
    return np.array(list(sol.y[:, -1]) + [0.0])


k = 2
print('  L [m]   x, px, y, py, z の (tthin - 厳密) [符号つき]')
for Lh in (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125):
    e = tthin_gen(c0, Lh, KAPPA[k] * Lh, k) - exact(c0, Lh, KAPPA[k], k)
    print('  %6.4f  ' % Lh + '  '.join('%+.2e' % v for v in e[:5]))


print('\n== (G) 力勾配型の c3、c5 の導出 (1 自由度、H = p^2/2 + V(x)、初期値 x = 0、p = 0。V_n = V の n 階微分) ==')
print("   分割: D(aL) K(1/2) D(tL) C D(tL) K(1/2) D(aL)、a = 1/2 - t、t = sqrt(3)/6 = 1/sqrt(12)、C: dp = df_C/dx、f_C = c3 h^3 V'^2 - 2 c5 h^5 V'^2 V''")
print('   厳密な流れとの差 (分割 - 厳密) を h のべきで展開 (s3 = sqrt(3))')
def _g07():
    h, x, p, c3, c5 = sp.symbols('h x p c3 c5')
    N = 5
    V = sp.symbols('V0:%d' % (N + 2))   # V[n] = n 階微分 at x
    s3 = sp.Symbol("s3"); t = s3/6
    a = sp.Rational(1, 2) - t
    X = sp.Symbol('X')

    # 1 自由度: H = p^2/2 + V(x)。関数は (x, p) と V_n の多項式。 d/dx V_n = V_{n+1}
    def dx(f):
        r = sp.diff(f, x)
        for n in range(N + 1):
            r += sp.diff(f, V[n]) * V[n + 1]
        return r
    def Dop(f):  # ハミルトンの流れの微分 = p ∂_x f - V1 ∂_p f
        return p * dx(f) - V[1] * sp.diff(f, p)
    def red(e):
        e = sp.expand(e)
        return sp.expand(sp.rem(sp.Poly(e, s3), sp.Poly(s3**2 - 3, s3)).as_expr()) if e.has(s3) else e
    def trunc(e, n=5):
        e = sp.expand(e)
        return sum(red(e.coeff(h, i)) * h**i for i in range(n + 1))
    def flow_exact(f, n=5):
        s, term = 0, f
        for i in range(n + 1):
            s += term * h**i / sp.factorial(i)
            term = Dop(term)
        return trunc(s, n)
    xe, pe = [flow_exact(f).subs({p: 0, x: 0}) for f in (x, p)]

    def Vder(n, d):   # V^{(n)}(x + d)
        return sum(V[n + j] * d**j / sp.factorial(j) for j in range(N + 1 - n))

    def apply(state, kind, arg):
        xs, ps = state
        if kind == 'D':
            return trunc(xs + arg * h * ps), ps
        if kind == 'K':        # p -= arg * h * V'(x)
            return xs, trunc(ps - arg * h * Vder(1, xs))
        if kind == 'C':
            d = xs
            v1, v2, v3 = Vder(1, d), Vder(2, d), Vder(3, d)
            # f_C = c3 h^3 V'^2 - 2 c5 h^5 V'^2 V''、Δp = ∂f/∂x
            dp = c3 * h**3 * 2 * v1 * v2 - 2 * c5 * h**5 * (2 * v1 * v2 * v2 + v1**2 * v3)
            return xs, trunc(ps + dp)
    seq = [('D', a), ('K', sp.Rational(1, 2)), ('D', t), ('C', 0), ('D', t), ('K', sp.Rational(1, 2)), ('D', a)]
    st = (sp.Integer(0), sp.Integer(0))
    for k_, ar in seq:
        st = apply(st, k_, ar)
    dxx = sp.expand(st[0] - xe)
    dpp = sp.expand(st[1] - pe)
    for n in range(0, 6):
        cx_, cp_ = (dxx.coeff(h, n)), (dpp.coeff(h, n))
        print('h^%d: x の差 =' % n, (cx_.subs(p, 0)), '|  p の差 =', sp.collect(sp.expand((cp_)), [p]) if n <= 3 else None)
    # p = 0 の 3 次: c3 を求める
    e3 = (dpp.coeff(h, 3).subs(p, 0))
    print('h^3 (p=0) 差:', e3)
    c3v = sp.solve(e3.subs(s3, sp.sqrt(3)), c3)
    print('c3 =', c3v, ' = 1/12 - t/4 =', sp.Rational(1,12) - sp.sqrt(3)/24)
    e5 = sp.expand(dpp.coeff(h, 5).subs({p: 0}).subs(s3, sp.sqrt(3)).subs(c3, c3v[0]))
    print('h^5 (p=0) の p の差:', sp.collect(e5, [V[1], V[2], V[3]]))
    sol = sp.solve(sp.Poly(e5, V[1], V[2], V[3]).coeffs(), c5, dict=True)
    print('c5 の解 (全係数を同時に 0):', sol, ' (9-5√3)/720 =', ((9 - 5 * sp.sqrt(3)) / 720))
    e5x = sp.simplify(dxx.coeff(h, 5).subs({p: 0}).subs(s3, sp.sqrt(3)).subs(c3, c3v[0]))
    print('h^5 (p=0) の x の差 :', e5x, ' h^4:', sp.simplify(dxx.coeff(h, 4).subs({p: 0}).subs(s3, sp.sqrt(3)).subs(c3, c3v[0])))

_g07()
