# -*- coding: utf-8 -*-
r"""Touschek 効果の導出の確認 (17_collective.md §10)。

(a) 非相対論的 Møller 断面積 (スピン平均) の整理: 1/s^4 + 1/c^4 - 1/(s^2 c^2) = 16/sin^4(chi) - 12/sin^2(chi) (s = sin(chi/2)、c = cos(chi/2)) を sympy で。
(b) 散乱後の縦運動量が閾値を超える立体角の積分 W(kappa) = 2 pi (kappa^2 - 1 - ln kappa) (片側) を、
    (i) 極軸を縦方向にとって sympy で厳密に、(ii) 極軸を入射の方向にとった (chi, phi) の 2 重積分の数値で。
(c) 相対運動量 q の 2 次元ガウス分布での平均 <|q| sigma_loss> を、モンテカルロと、J(c) の 1 次元積分 (Bessel I0) で比べる。
(d) 係数: 1 つの断面 (スライス) の寿命の式 1/tau = r_e^2 c N J / (gamma V* sqrt(l1 l2)) を、pysad の piwinski_slice_rate と比べる。
(e) 横 1 方向だけ運動量の広がりがある極限 (Bruck の平坦ビーム) の関数 H(eps) と、K64 の touckf (当てはめ) の形を比べる。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
import sympy as sp
from scipy import integrate, special

# ------------------------------------------------------------------ (a)
chi = sp.symbols('chi', positive=True)
s, c = sp.sin(chi / 2), sp.cos(chi / 2)
lhs = 1 / s**4 + 1 / c**4 - 1 / (s**2 * c**2)
rhs = 16 / sp.sin(chi)**4 - 12 / sp.sin(chi)**2
print('== (a) Møller (スピン平均): 1/s^4 + 1/c^4 - 1/(s^2 c^2) - (16/sin^4 - 12/sin^2) の数値 (chi = 0.3, 1.1, 2.5) ==')
print([float((lhs - rhs).subs(chi, v)) for v in (0.3, 1.1, 2.5)])
S2 = sp.symbols('S2', positive=True)   # S2 = sin^2(chi) とおくと s^2 c^2 = S2/4、s^4 + c^4 = 1 - S2/2
print('  s^2 c^2 = S2/4、s^4 + c^4 = 1 - 2 s^2 c^2 = 1 - S2/2 から:', sp.simplify((1 - S2 / 2) / (S2 / 4)**2 - 1 / (S2 / 4) - (16 / S2**2 - 12 / S2)), '(0 なら一致)')

# ------------------------------------------------------------------ (b)
u, psi, kap = sp.symbols('u psi kappa', positive=True)
sth = sp.sqrt(1 - u**2)                      # sin(theta)、u = |cos(theta)| = |n_z|
# 方位角 psi の平均: sin^2 chi = 1 - sin^2(theta) cos^2(psi)
a_ = sp.symbols('a_', positive=True)         # a_ = sin^2(theta) < 1
print('\n== (b) 立体角の積分 ==')
# 手計算の結果 2 pi/sqrt(1-a) と (1 + u^2) pi/u^3 を数値で確認 (a = 1 - u^2)
for uu in (0.2, 0.5, 0.9):
    aa = 1 - uu**2
    n1 = integrate.quad(lambda p: 1 / (1 - aa * math.cos(p)**2), 0, 2 * math.pi)[0]
    n2 = integrate.quad(lambda p: 1 / (1 - aa * math.cos(p)**2)**2, 0, 2 * math.pi)[0]
    print('  u = %.1f: ∫dψ/(..) = %.6f (2π/u = %.6f)、∫dψ/(..)^2 = %.6f (π(1+u²)/u³ = %.6f)' % (uu, n1, 2 * math.pi / uu, n2, math.pi * (1 + uu**2) / uu**3))
integrand = 4 * sp.pi * (1 + u**2) / u**3 - 3 * 2 * sp.pi / u
print('  4∫dψ/(..)^2 - 3∫dψ/(..) = ', sp.simplify(integrand))
Wk = sp.integrate(sp.expand(sp.simplify(integrand)), (u, 1 / kap, 1))
print('  W(κ) (片側: n_z > 1/κ) = ∫_{1/κ}^1 du (...) =', sp.simplify(Wk), '  [2π(κ²-1-ln κ) との差:', sp.simplify(Wk - 2 * sp.pi * (kap**2 - 1 - sp.log(kap))), ']')


def Wnum(k):
    """(chi, phi) の 2 重積分: n = (cos chi, sin chi cos phi, sin chi sin phi)、n_z = sin chi sin phi > 1/kappa の側だけ"""
    f = lambda ph, ch: math.sin(ch) * (4 / math.sin(ch)**4 - 3 / math.sin(ch)**2) * (1.0 if math.sin(ch) * math.sin(ph) > 1 / k else 0.0)
    # 被積分関数に段差があるので phi の範囲を陽に書く: sin(phi) > 1/(k sin chi)
    def inner(ch):
        sc = math.sin(ch)
        r = 1 / (k * sc)
        if r >= 1:
            return 0.0
        p1 = math.asin(r)
        return (math.pi - 2 * p1) * math.sin(ch) * (4 / sc**4 - 3 / sc**2)
    chmin = math.asin(1 / k)
    return integrate.quad(inner, chmin, math.pi - chmin, epsabs=0, epsrel=1e-11, limit=400)[0]


for k in (1.5, 3.0, 10.0):
    print('  κ = %5.1f: (χ,φ) の 2 重積分 = %.8f、2π(κ²-1-ln κ) = %.8f' % (k, Wnum(k), 2 * math.pi * (k**2 - 1 - math.log(k))))

# ------------------------------------------------------------------ (c)
def W(k):
    return np.where(k > 1, 2 * np.pi * (k**2 - 1 - np.log(np.maximum(k, 1e-300))), 0.0)


def ive0(z):
    return special.ive(0, z)


l1, l2, cth = 4.0e-6, 1.0e-6, 1.2e-3                 # P の固有値 (mc 単位の運動量の分散)、閾値 c = 2 beta a
a1, a2 = 1 / (4 * l1), 1 / (4 * l2)
print('\n== (c) <|q| σ_loss> (片側): モンテカルロ と J(c) の 1 次元積分 ==')
print('l1 = %.1e、l2 = %.1e [(mc)^2]、c = %.2e (閾値)' % (l1, l2, cth))
rng = np.random.default_rng(3)
nsamp = 4_000_000
q = rng.normal(size=(nsamp, 2)) * np.sqrt([2 * l1, 2 * l2])      # 相対運動量 q = p1 - p2 ~ N(0, 2P)
qm = np.hypot(q[:, 0], q[:, 1])
# σ_loss = 4 r_e^2 W(κ)/k^4 (片側)、|q| σ = 4 r_e^2 W(κ)/q^3: r_e^2 = 1 の単位で
mc_val = np.mean(4 * W(qm / cth) / qm**3)
mc_err = np.std(4 * W(qm / cth) / qm**3) / math.sqrt(nsamp)


def J(cc):
    f = lambda k: W(np.array(k / cc)) / k**2 * math.exp(-k**2 * a1) * ive0(k**2 * (a2 - a1) / 2)
    return integrate.quad(f, cc, math.sqrt(60 / a1), epsabs=0, epsrel=1e-10, limit=500)[0]


Jv = J(cth)
quad_val = 2.0 / math.sqrt(l1 * l2) * Jv
print('モンテカルロ (%.0e 個): %.6e ± %.1e' % (nsamp, mc_val, mc_err))
print('1 次元積分 (2/√(l1 l2)) J(c): %.6e' % quad_val)

# ------------------------------------------------------------------ (d)
print('\n== (d) 1/τ の係数: 自作の式と pysad.touschek.piwinski_slice_rate (n = 1 スライス) ==')
from pysad import touschek as T
E = 3e9
mc2 = 0.51099895e6
g = E / mc2
p0 = math.sqrt(g**2 - 1)
beta = p0 / g
sx, sy, sz = 100e-6, 10e-6, 5e-3                       # 実験室系の大きさ
xx = np.diag([sx**2, sy**2, (g * sz)**2])              # 静止系: z* = γ z
pl = np.diag([(1e-4)**2, (3e-5)**2, 1e-12])           # 静止系の運動量の共分散 / p0^2 (横 x', y'、縦は無視)
Np, ap, al_, circ = 1e10, 0.02, 1.0, 1.0
sl = dict(p0=p0, h0=g, xx=xx, pl=pl, al=al_)
rate_py = T.piwinski_slice_rate(sl, ap, ap, circ, Np)
re = 2.8179403262e-15
cc = 2.99792458e8
P = pl[:2, :2] * p0**2
lam1, lam2 = P[0, 0], P[1, 1]
if lam2 > lam1:
    lam1, lam2 = lam2, lam1
a1, a2 = 1 / (4 * lam1), 1 / (4 * lam2)
Vs = (4 * math.pi)**1.5 * math.sqrt(np.linalg.det(xx))
Jc = 2 * J(2 * beta * ap)                              # 正負の受容 (同じ) の 2 つ
rate_me = re**2 * cc * Np / (g * Vs * math.sqrt(lam1 * lam2)) * Jc
print('自作: %.6e [1/s]、pysad: %.6e [1/s]、比 %.6f' % (rate_me, rate_py, rate_me / rate_py))

# ------------------------------------------------------------------ (e)
print('\n== (e) 横 1 方向 (l2 -> 0) の極限: H(ε) = ∫_ε^∞ dt t^-2 (t/ε - 1 - (1/2) ln(t/ε)) e^-t と touckf(ε) ==')
def H(eps):
    f = lambda t: t**-2 * (t / eps - 1 - 0.5 * math.log(t / eps)) * math.exp(-t)
    return integrate.quad(f, eps, eps + 60, epsabs=0, epsrel=1e-10, limit=500)[0]


print('  ε        H(ε)        touckf(ε)   H/touckf')
for eps in (1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0, 10.0):
    hv, tk = H(eps), float(T.touschek_f(eps))
    print('  %-8g %.5e %.5e %.5f' % (eps, hv, tk, hv / tk))
