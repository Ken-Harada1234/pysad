# -*- coding: utf-8 -*-
r"""01 章 (変数と約束) の途中式を sympy で確かめる。

(a) (1.2) -> (1.3): L = L_t d(ct)/ds の変形。
(b) (1.4) の偏微分、H_t = x'p_x + y'p_y + t'p_t - L の整理、(1.5) の根号の恒等式。
(c) (1.8): dH/dp、E の微分。
(d) (1.9): v/v0 と dv の 2 通りの式、dv の δ 微分 = -h0/h1^3、小さい δ の係数。
(e) 5 節: ρ = ∞, a = 0 で x' = dH/dp_x = q_x / sqrt(1 - q_x^2 - q_y^2)。
"""
import sympy as sp
import _common as C

C.start(__file__)

# (a) ------------------------------------------------------------------
print('(a) (1.2) -> (1.3)')
ctp, xp, yp, u = sp.symbols('ctp xp yp u', positive=True)   # ctp = d(ct)/ds, u = 1 + x/rho
xd, yd, sd = xp / ctp, yp / ctp, 1 / ctp                      # d/d(ct) = (1/ctp) d/ds
Lt_root = sp.sqrt(1 - xd**2 - yd**2 - u**2 * sd**2)
Pi = sp.sqrt(ctp**2 - xp**2 - yp**2 - u**2)
print('  ctp * sqrt(1 - xd^2 - yd^2 - u^2 sd^2) = sqrt(ctp^2 - xp^2 - yp^2 - u^2) : ',
      sp.simplify((ctp * Lt_root)**2 - Pi**2) == 0)
print('  a_x xd * ctp = a_x xp (自明): ', sp.simplify(sp.Symbol('ax') * xd * ctp - sp.Symbol('ax') * xp) == 0)

# (b) ------------------------------------------------------------------
print('(b) (1.4), H_t, (1.5)')
tp, mh, ax, ay, as_, c = sp.symbols('tp mh ax ay a_s c', positive=True)  # tp = t' = dt/ds
Pi = sp.sqrt((c * tp)**2 - xp**2 - yp**2 - u**2)
L = -mh * Pi + ax * xp + ay * yp + u * as_
px = sp.diff(L, xp); py = sp.diff(L, yp); pt = sp.diff(L, tp)
print('  dPi/dx\' = -x\'/Pi: ', sp.simplify(sp.diff(Pi, xp) + xp / Pi) == 0)
print('  px - ax = mh x\'/Pi: ', sp.simplify(px - ax - mh * xp / Pi) == 0)
print('  pt = -mh c^2 t\'/Pi: ', sp.simplify(pt + mh * c**2 * tp / Pi) == 0)
Ht = xp * px + yp * py + tp * pt - L
print('  x\'(px-ax)+y\'(py-ay)+t\' pt = -mh (Pi^2 + u^2)/Pi: ',
      sp.simplify(xp * (px - ax) + yp * (py - ay) + tp * pt + mh * (Pi**2 + u**2) / Pi) == 0)
print('  H_t = -mh u^2/Pi - u a_s: ', sp.simplify(Ht - (-mh * u**2 / Pi - u * as_)) == 0)
rad = (pt / c)**2 - (px - ax)**2 - (py - ay)**2 - mh**2
print('  (pt/c)^2-(px-ax)^2-(py-ay)^2-mh^2 = mh^2 u^2/Pi^2: ', sp.simplify(rad - mh**2 * u**2 / Pi**2) == 0)

# (c) ------------------------------------------------------------------
print('(c) (1.8)')
p, pxx, pyy, m, c_, p0m, vv0 = sp.symbols('p pxx pyy m c p0 vv0', positive=True)
E = sp.sqrt(m**2 * c_**4 + p0m**2 * p**2 * c_**2)
print('  dE/dp = p0^2 p c^2 / E: ', sp.simplify(sp.diff(E, p) - p0m**2 * p * c_**2 / E) == 0)
root = sp.sqrt(p**2 - pxx**2 - pyy**2)
print('  d(-root)/dp = -p/root: ', sp.simplify(sp.diff(-root, p) + p / root) == 0)

# (d) ------------------------------------------------------------------
print('(d) (1.9)')
d, P0 = sp.symbols('delta P0', positive=True)        # P0 = beta0 gamma0 (mc 単位)
h0 = sp.sqrt(1 + P0**2); h1 = sp.sqrt(1 + P0**2 * (1 + d)**2)
dv = 1 - (1 + d) * h0 / h1
dv2 = -d * (2 + d) / (h1 * (h1 + (1 + d) * h0))
print('  h1^2 - (1+d)^2 h0^2 = -d(2+d): ', sp.simplify(h1**2 - (1 + d)**2 * h0**2 + d * (2 + d)) == 0)
print('  dv == 閉じた式: ', sp.simplify(dv - dv2) == 0)
print('  d(dh1/dd): dh1/dd = P0^2 (1+d)/h1: ', sp.simplify(sp.diff(h1, d) - P0**2 * (1 + d) / h1) == 0)
print('  d(dv)/d(delta) + h0/h1^3 = 0: ', sp.simplify(sp.diff(dv, d) + h0 / h1**3) == 0)
print('  d=0 で d(dv)/d(delta) = -1/gamma0^2: ', sp.simplify(sp.diff(dv, d).subs(d, 0) + 1 / h0**2) == 0)
print('  v/v0 = (1+d) h0/h1 : beta = p1/h1, beta0 = P0/h0, p1 = P0 (1+d) -> ',
      sp.simplify(((P0 * (1 + d) / h1) / (P0 / h0)) - (1 + d) * h0 / h1) == 0)
g0 = 3.0e9 / 0.51099895e6
print('  3 GeV 電子: gamma0 = %.2f、1/gamma0^2 = %.4e' % (g0, 1 / g0**2))

# (e) ------------------------------------------------------------------
print('(e) §5: H = -sqrt(p^2 - px^2 - py^2), x\' = dH/dpx, q = px/p')
qx, qy = sp.symbols('qx qy')
H = -sp.sqrt(p**2 - pxx**2 - pyy**2)
xprime = sp.diff(H, pxx).subs({pxx: qx * p, pyy: qy * p})
print('  x\' = qx/sqrt(1-qx^2-qy^2): ', sp.simplify(xprime - qx / sp.sqrt(1 - qx**2 - qy**2)) == 0)
