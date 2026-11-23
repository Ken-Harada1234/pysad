# -*- coding: utf-8 -*-
"""放射の章 (13_radiation.md) の補完。

(a) 放射スペクトルの 3 つの積分 (13.1) を、メリン変換 (ガンマ関数) で閉じた形にする。
(b) dg = -de/(1 + 2 de) の導出: 曲率半径が一定のとき dp/ds ∝ -p⁴ の解 (1 + 3 de)^(-1/3) - 1 との 2 次までの一致。
(c) dv の分母: tradk1 の (h2 + p2) と (1.9) の (h2 + pr h0) の差 (低エネルギーで効く)。
(d) tradke の微分 (減衰の行列): 連鎖律で dg の微分を導き、K64 のコードの各項と比べる (sympy)。B138 の 4 点のうち導出できるもの。
(e) 分割数の式 (tbende.f:462-465): 1 + N_γ になること。
"""
import math
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)

# ---------------------------------------------------------------------------------------------------- (a)
print('== (a) 放射スペクトルのモーメント (13.1) のメリン変換 ==')
mp.mp.dps = 20


def Fnum(mu):
    """∫0^∞ ξ^mu F(ξ) dξ、F(ξ) = ξ ∫_ξ^∞ K_{5/3}(t) dt (scipy の数値積分)"""
    from scipy.integrate import quad
    from scipy.special import kv
    G = lambda x: quad(lambda t: kv(5.0 / 3.0, t), x, np.inf, epsabs=1e-14, epsrel=1e-12)[0]
    return quad(lambda x: x**(mu + 1) * G(x), 0, 60, epsabs=1e-13, epsrel=1e-11, limit=200)[0]


def Fexact(mu):
    """部分積分: ∫ ξ^{s-1} G dξ = (1/s) ∫ t^s K dt、s = mu + 2。∫0^∞ t^s K_ν dt = 2^{s-1} Γ((1+s+ν)/2) Γ((1+s-ν)/2)"""
    s = mu + 2
    nu = mp.mpf(5) / 3
    return (1 / mp.mpf(s)) * 2**(s - 1) * mp.gamma((1 + s + nu) / 2) * mp.gamma((1 + s - nu) / 2)


for mu, lab, closed in ((0, '∫F dξ', 8 * mp.pi / (9 * mp.sqrt(3))), (-1, '∫F/ξ dξ', 5 * mp.pi / 3), (1, '∫ξ F dξ', 55 * mp.pi / 81)):
    print('  %-9s: 数値積分 = %s、ガンマ関数の式 = %s、閉じた形 = %s' % (lab, '%.12f' % Fnum(mu), mp.nstr(Fexact(mu), 12), mp.nstr(closed, 12)))
print('  (ガンマ関数の式: 2^{mu+1}/(mu+2) Γ(mu/2 + 7/3) Γ(mu/2 + 2/3)。Γ(x)Γ(1-x) = π/sin(πx) で 8π/(9√3)、5π/3、55π/81 になる。)')

# ---------------------------------------------------------------------------------------------------- (b)
print('\n== (b) dg = -de/(1 + 2 de) の導出 ==')
de = sp.symbols('de', positive=True)
exact = (1 + 3 * de) ** sp.Rational(-1, 3) - 1
k64 = -de / (1 + 2 * de)
print('  厳密 (曲率半径一定、dp/ds ∝ -p⁴): (1 + 3 de)^(-1/3) - 1 = %s' % sp.series(exact, de, 0, 4))
print('  K64: -de/(1 + 2 de)                                 = %s' % sp.series(k64, de, 0, 4))


def solve_u(de_v):
    f = lambda t, u: [-de_v * u[0]**4]
    return solve_ivp(f, [0, 1], [1.0], rtol=1e-13, atol=1e-16).y[0, -1] - 1


for dv_ in (5e-5, 1e-2, 5e-2, 0.2):
    print('  de = %g: ODE (du/dτ = -de u⁴) の dg = %.10e、(1+3de)^(-1/3)-1 = %.10e、K64 = %.10e (K64 - ODE = %.1e)' % (
        dv_, solve_u(dv_), float(exact.subs(de, dv_)), float(k64.subs(de, dv_)), float(k64.subs(de, dv_)) - solve_u(dv_)))

# ---------------------------------------------------------------------------------------------------- (c)
print('\n== (c) dv の分母: tradk1 の h2 + p2 と、(1.9) の h2 + pr h0 ==')


def dv_exact(p0, pr):
    h0 = math.sqrt(1 + p0**2)
    h = math.sqrt(1 + (p0 * pr)**2)
    return 1 - (p0 * pr / h) / (p0 / h0)          # 1 - β/β0


def dv_19(p0, pr):
    h0 = math.sqrt(1 + p0**2)
    h = math.sqrt(1 + (p0 * pr)**2)
    return -(pr - 1) * (1 + pr) / h / (h + pr * h0)


def dv_tradk(p0, pr):
    h = math.sqrt(1 + (p0 * pr)**2)
    return -(pr - 1) * (1 + pr) / h / (h + p0 * pr)


for p0, pr, lab in ((5871.0, 1 - 1e-4, '3 GeV (p0 = 5871 mc)、δ = -1e-4'), (50.0, 0.99, 'p0 = 50 mc (25 MeV)、δ = -0.01'), (1.0, 1.1, 'p0 = 1 mc (0.51 MeV/c)、δ = +0.1')):
    e, a, b = dv_exact(p0, pr), dv_19(p0, pr), dv_tradk(p0, pr)
    print('  %-34s: 1 - β/β0 = %+.6e、(1.9) = %+.6e、tradk1 の式 = %+.6e (tradk1 - 厳密 = %+.2e、相対 %.1e)' % (lab, e, a, b, b - e, abs((b - e) / e)))
print('  (2 つの式の差: h + pr h0 と h + p0 pr の差 pr (h0 - p0) = pr/(h0 + p0)。相対 O(1/γ²)。)')

# ---------------------------------------------------------------------------------------------------- (d)
print('\n== (d) tradke の微分 (連鎖律) ==')
theta, h1, al1, pr, p0s, anrad, cuc, cave = sp.symbols('theta h1 al1 pr p0 anrad cuc cave', positive=True)
anp = anrad * h1 * theta
uct = cuc * h1**3 / p0s / al1
uc = uct * theta
dg = -cave * anp * uc
# h1 = sqrt(1 + (p0 pr)²)
h1_of_pr = sp.sqrt(1 + (p0s * pr)**2)
dh1_dpr = sp.diff(h1_of_pr, pr)
dh1r = sp.simplify(dh1_dpr / h1_of_pr)
print('  dh1r = (1/h1) dh1/dpr = %s   (コード: p p0/h1² = p0² pr/h1²)' % dh1r)
# 連鎖律: d(dg) = ∂dg/∂θ dθ + ∂dg/∂al1 dal + ∂dg/∂pr dpr
ddg_theta = sp.diff(dg, theta)
ddg_al = sp.diff(dg, al1)
# コードの式
danp_code = anrad * h1 * 1                  # danp = anrad*h1*dtheta (dθ の係数)
duc_code = uct                              # duc = uct*dtheta (dθ の係数)
ddg_theta_code = -cave * (danp_code * uc + anp * duc_code)
print('  ∂dg/∂θ - コードの -cave (danp uc + anp duc)|_θ = %s' % sp.simplify(ddg_theta - ddg_theta_code))
duc_al_code = -uc / al1                     # duc = ... - uc*dal/al1
ddg_al_code = -cave * (anp * duc_al_code)
print('  ∂dg/∂al1 - コードの -cave anp (-uc/al1)       = %s' % sp.simplify(ddg_al - ddg_al_code))
# pr 微分 (h1 経由) : danp(6) += anp dh1r、duc(6) += 3 uc dh1r
dg_h = dg.subs(h1, h1_of_pr)
ddg_pr = sp.diff(dg_h, pr)
ddg_pr_code = (-cave * (anp * dh1r * uc + anp * 3 * uc * dh1r)).subs(h1, h1_of_pr)
print('  ∂dg/∂pr (h1 を通して) - コードの -cave (anp dh1r uc + anp 3 uc dh1r) = %s' % sp.simplify(ddg_pr - ddg_pr_code))
print('  → ddg = -cave (danp uc + anp duc) は、dg = -cave anp uc (anp = anrad h1 θ、uc = cuc h1³ θ/(p0 al1)) の連鎖律どおり。')

# xpa の微分 (B138 (4))
a_ = sp.symbols('a')
xp = sp.Matrix([sp.Function('X')(a_), sp.Function('Y')(a_), sp.Function('Z')(a_)])
prf = sp.Function('P')(a_)
xpa_f = sp.sqrt(xp.dot(xp)) / prf**2
d_xpa = sp.diff(xpa_f, a_)
dxp = xp.diff(a_)
claim = (xp.dot(dxp)) / (xpa_f * prf**4) - 2 * xpa_f * sp.diff(prf, a_) / prf
print('  xpa = |xp|/pr² の微分: d xpa - [(xp·dxp)/(xpa pr⁴) - 2 xpa dpr/pr] = %s' % sp.simplify(d_xpa - claim))
print('  K64 の dxpa = (xp·dxp)/xpa/pr² (第 1 項が pr² で、pr⁴ ではない。B138 (4))。δ = 0 (pr = 1) では同じ。')

# dddpx = ½ [(d px - d pxr0) dg + dpx ddg] (B138 (2))
print('  ddpx = ½ dpx dg の微分: d ddpx = ½ [(d px - d pxr0) dg + dpx ddg]。K64 は第 2 項が ddpx*ddg (= ½ dpx dg ddg) で、dpx ddg ではない (B138 (2))。')

# ---------------------------------------------------------------------------------------------------- (e)
print('\n== (e) 分割数 (tbende.f:462-465) ==')
epsrad = 1e-6
# 例: 3 GeV、L = 1.2 m、ANGLE = 2π/16 の偏向磁石
anrad_v = 5 / (2 * math.sqrt(3)) / 137.035999
h0_v = 3.0e9 / 0.51099895e6
phib = 2 * math.pi / 16
Ng = phib * h0_v * anrad_v
print('  N_γ = φ h0 anrad = %.3f。1 + int(N_γ/epsrad/1e6 · emidiv·emidib) = %d (EMITDIV = EMITDIVB = 1、epsrad = 1e-6)' % (Ng, 1 + int(abs(Ng) / epsrad / 1e6)))
loss = 4.8850e-5
print('  1 台の相対損失 (追跡) = %.4e。損失 / epsrad = %.1f (nrad = 1 + int(L ρ_b crad/epsr1 (h0 b)²) が 1 刻みの相対損失を epsrad 以下にする数だとすると、これに近い)' % (loss, loss / epsrad))
