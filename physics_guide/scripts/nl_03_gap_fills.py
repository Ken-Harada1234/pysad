# -*- coding: utf-8 -*-
r"""18 章の補足の途中式の確認 (§3.4 のヤコビ恒等式と 3 つの指数の積、§3.8 の起点の独立性、§3.10 の共鳴近傍のハミルトニアン、§3.12 の Hann 窓の裾)。

(a) [:a:, :b:] g = :[a,b]: g (ヤコビ恒等式) を一般の 3 次式 a, b と g について sympy で確認。
(b) e^a e^b e^c = exp(a + b + c + ½([a,b]+[a,c]+[b,c])) (ζ への作用を 3 次まで) を一般の 3 次式 a, b, c で確認。
(c) 1 次元の K_2200 (18.35) が、要素の並びを巡回的にずらした (起点を動かした) ときに変わらないこと。ランダムな要素、整数から離れたチューン。
(d) 3 ν = N + 3δ の近くの 1 周の写像: 共鳴の調和だけを残した F' から作った写像 e^{:F':} R の周期 3 の不動点が、(18.38) の K の不動点と一致すること (数値)。
(e) Hann 窓の核の括弧 [D(Δ) - ½D(Δ-1/N) - ½D(Δ+1/N)] の位相・符号の因子と、裾 -1/(π N^2 Δ^3)。
"""
import itertools
import math
import random
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import _common as C

C.start(__file__)
z, zb = sp.symbols('z zb')


def pb(f, g):
    """(18.9): [f, g] = 2i (f_z g_zb - f_zb g_z)"""
    return sp.expand(2 * sp.I * (sp.diff(f, z) * sp.diff(g, zb) - sp.diff(f, zb) * sp.diff(g, z)))


def rand_poly(deg, tag):
    return sum(sp.Symbol('%s%d%d' % (tag, i, deg - i)) * z**i * zb**(deg - i) for i in range(deg + 1))


print('== (a) ヤコビ: [a,[b,g]] - [b,[a,g]] - [[a,b],g] ==')
a3, b3, g3 = rand_poly(3, 'a'), rand_poly(3, 'b'), rand_poly(3, 'g')
print('差 =', sp.simplify(pb(a3, pb(b3, g3)) - pb(b3, pb(a3, g3)) - pb(pb(a3, b3), g3)))


def trunc(expr, maxdeg):
    p = sp.Poly(sp.expand(expr), z, zb)
    return sum(c * z**m[0] * zb**m[1] for m, c in p.terms() if m[0] + m[1] <= maxdeg)


def expL(f, g, maxdeg, nmax=3):
    """e^{:f:} g = g + [f,g] + ½[f,[f,g]] + ... (maxdeg 次まで)"""
    out, term = g, g
    for n in range(1, nmax + 1):
        term = pb(f, term) / n
        out += term
    return trunc(out, maxdeg)


print('\n== (b) e^a e^b e^c = exp(a+b+c+½([a,b]+[a,c]+[b,c])) (ζ への作用、3 次まで) ==')
cc3 = rand_poly(3, 'c')
lhs = expL(a3, expL(b3, expL(cc3, z, 3), 3), 3)
dgen = a3 + b3 + cc3 + sp.Rational(1, 2) * (pb(a3, b3) + pb(a3, cc3) + pb(b3, cc3))
rhs = expL(dgen, z, 3)
print('差 (ζ の像、3 次まで) =', sp.simplify(sp.expand(lhs - rhs)))

print('\n== (c) 起点の独立性: K_2200 の非 h4 部分 -18 ImP3 - 6 ImP1 - 9|h30|^2 cot3πν - 3|h21|^2 cotπν ==')
random.seed(3)
nu = 0.2347
mu = 2 * math.pi * nu


def kfun(elems):
    """elems = [(位相 μ_e, A30_e, A21_e)] の順 (最初に通る要素が先頭)。c30 = A30 e^{3iμ_e}、c21 = A21 e^{iμ_e}"""
    c30 = [A * np.exp(3j * m) for m, A, B in elems]
    c21 = [B * np.exp(1j * m) for m, A, B in elems]
    P3 = sum(c30[l] * np.conj(c30[e]) for l in range(len(elems)) for e in range(l))
    P1 = sum(c21[l] * np.conj(c21[e]) for l in range(len(elems)) for e in range(l))
    h30, h21 = sum(c30), sum(c21)
    return (-18 * P3.imag - 6 * P1.imag - 9 * abs(h30)**2 / math.tan(3 * mu / 2) - 3 * abs(h21)**2 / math.tan(mu / 2))


n = 5
phs = sorted(random.uniform(0, mu) for _ in range(n))     # 起点から測った位相 (0 < μ_e < μ)
Ams = [random.uniform(0.3, 1.2) for _ in range(n)]
elems = [(phs[i], Ams[i], 0.37 * Ams[i] * (1 + 0.2 * i)) for i in range(n)]
print('起点 0:', kfun(elems))
for k in range(1, n):
    # 起点を k 個目の要素の後ろに動かす: 位相を μ_k だけ引き、先頭の k 個は μ を足して最後に回す
    sh = elems[k - 1][0]
    new = [(m - sh, A, B) for (m, A, B) in elems[k:]] + [(m - sh + mu, A, B) for (m, A, B) in elems[:k]]
    print('起点を要素 %d の後ろへ:' % k, kfun(new))

print('\n  先頭の 1 要素だけを最後に回した ΔP3 = X e^{iθ} - conj(X)、Δ|h|^2 = 2Re[X(e^{iθ}-1)] (X = c1 conj(S\'))')
c30 = [A * np.exp(3j * m) for m, A, B in elems]
c1, Sp = c30[0], sum(c30[1:])
X = c1 * np.conj(Sp)
th = 3 * mu
P3_old = sum(c30[l] * np.conj(c30[e]) for l in range(n) for e in range(l))
c30n = c30[1:] + [c1 * np.exp(1j * th)]
P3_new = sum(c30n[l] * np.conj(c30n[e]) for l in range(n) for e in range(l))
print('  ΔP3 (直接) =', P3_new - P3_old, ' X e^{iθ} - conj(X) =', X * np.exp(1j * th) - np.conj(X))
print('  Im ΔP3 =', (P3_new - P3_old).imag, ' Im[X(e^{iθ}+1)] =', (X * (np.exp(1j * th) + 1)).imag)
print('  Δ|h|^2 (直接) =', abs(sum(c30n))**2 - abs(sum(c30))**2, ' 2Re[X(e^{iθ}-1)] =', 2 * (X * (np.exp(1j * th) - 1)).real)

print('\n== (d) 共鳴近傍: 1 周の写像 e^{:F\':} R の周期 3 の不動点 vs (18.38) の K の不動点 ==')
# ν = 1/3 + δ (N = 1)。F' = -(h ζ^3 + c.c.) - c J^2 (K = 2πδ J + c J^2 + h ζ^3 + c.c.)
h, chi, cc, delta = 0.02, 0.4, -0.5, 0.001
mu_ = 2 * math.pi * (1 / 3 + delta)


def Kham(J, ph):
    return 2 * math.pi * delta * J + cc * J**2 + 2 * h * (2 * J)**1.5 * math.cos(3 * ph + chi)


def flow(J, ph):
    """H = c J^2 + 2h(2J)^{3/2} cos(3φ+χ) の時間 1 の流れ ((φ, J) 正準対)"""
    def rhs(t, y):
        Jv, pv = y
        dJ = -(-6 * h * (2 * Jv)**1.5 * math.sin(3 * pv + chi))
        dp = 2 * cc * Jv + 2 * h * 3 * math.sqrt(2) * math.sqrt(Jv) * math.cos(3 * pv + chi)
        return [dJ, dp]
    s = solve_ivp(rhs, [0, 1], [J, ph], rtol=1e-12, atol=1e-14)
    return s.y[0, -1], s.y[1, -1]


def turn(J, ph):
    J1, p1 = flow(J, ph)
    return J1, p1 + mu_       # R: φ -> φ + μ


def three(v):
    J, ph = v
    for _ in range(3):
        J, ph = turn(J, ph)
    return [J - v[0], (ph - v[1]) / (2 * math.pi) - round((ph - v[1]) / (2 * math.pi))]


# K の不動点: ∂K/∂φ = 0 (3φ+χ = π)、∂K/∂J = 0
def Kfix(J):
    return 2 * math.pi * delta + 2 * cc * J - 6 * math.sqrt(2) * h * math.sqrt(J) * 1.0
from scipy.optimize import brentq
roots = []
Js = np.linspace(1e-6, 0.05, 20001)
vals = [Kfix(j) for j in Js]
for i in range(len(Js) - 1):
    if vals[i] * vals[i + 1] < 0:
        roots.append(brentq(Kfix, Js[i], Js[i + 1]))
print('K の不動点 (3φ+χ = π) の J:', roots, ' (c を無視した (18.39): %.6e)' % (math.pi**2 * delta**2 / (18 * h**2)))
for J0 in roots:
    ph0 = (math.pi - chi) / 3
    sol = fsolve(three, [J0 * 1.02, ph0 + 0.01], xtol=1e-13)
    print('  写像 M^3 の不動点 J = %.6e (K の J_fp = %.6e、相対差 %.2e)、φ = %.6f (K: %.6f)' % (sol[0], J0, sol[0] / J0 - 1, sol[1] % (2 * math.pi / 3), ph0 % (2 * math.pi / 3)))

print('\n== (e) Hann 窓の括弧 ==')
for N in (256, 1024):
    for k in (2.5, 5.5, 20.5, 80.5):
        Dl = k / N
        Dk = lambda d: np.exp(1j * math.pi * d * (N - 1)) * math.sin(math.pi * N * d) / (N * math.sin(math.pi * d))
        S = Dk(Dl) - 0.5 * Dk(Dl - 1 / N) - 0.5 * Dk(Dl + 1 / N)
        br = 1 / math.sin(math.pi * Dl) - 0.5 * np.exp(1j * math.pi / N) / math.sin(math.pi * (Dl - 1 / N)) - 0.5 * np.exp(-1j * math.pi / N) / math.sin(math.pi * (Dl + 1 / N))
        S2 = np.exp(1j * math.pi * Dl * (N - 1)) * math.sin(math.pi * N * Dl) / N * br
        print('N=%d NΔ=%5.1f: |S|=%.4e、括弧の式から %.4e、予想 1/(π(NΔ)^3)=%.4e、括弧 %.4e vs -1/(πN²Δ³) = %.4e' % (
            N, k, abs(S), abs(S2), 1 / (math.pi * k**3), br.real, -1 / (math.pi * N**2 * Dl**3)))
