# -*- coding: utf-8 -*-
r"""01 章の途中式の確認 (本文に書き足した展開)。

(a) §3: (1.2) の導出 (速さの式と規格化)。デカルト座標の -mc^2/γ から出る根号が (1.3) の根号 (s を独立変数にしたもの) になること。
(b) §3: p_t = -E/P0 (Π = ct'/γ)、p^2 = p_t^2/c^2 - m^2c^2/P0^2 = (P/P0)^2。
(c) §5: 四極の p_x' = -K x ((1.1) で ρ = ∞、a_s = -K(x^2 - y^2)/2)。
(d) §2: E/(P0 v0) = h0 h1 / p0^2 (P0 = p0 m c)。
"""
import sympy as sp
import _common as C

C.start(__file__)

# (a) ------------------------------------------------------------------
print('(a) 速さと根号')
xd, yd, sd, hx, c, m, P0 = sp.symbols('xd yd sd hx c m P0', positive=True)   # xd = dx/d(ct)、hx = 1 + x/ρ
v2c2 = xd**2 + yd**2 + hx**2 * sd**2
tp = sp.symbols('tp', positive=True)   # ct' = d(ct)/ds
xp, yp = sp.symbols('xp yp', real=True)
# x' = xd * ct'、y' = yd * ct'、1 = sd * ct' (ds/ds = 1)
sub = {xd: xp / tp, yd: yp / tp, sd: 1 / tp}
root_t = sp.sqrt(1 - v2c2)               # (1.2) の根号
root_s = sp.sqrt(tp**2 - xp**2 - yp**2 - hx**2) / tp   # (1.3) の根号を ct' で割ったもの
print('  (1.2) の根号を (x\', y\', ct\') で書いたもの - (1.3) の根号/ct\' =', sp.simplify(root_t.subs(sub) - root_s))
print('  L = L_t ct\' → (mc/P0) sqrt(...) ct\' = (mc/P0) sqrt((ct\')^2 - x\'^2 - y\'^2 - hx^2): 一致 (上の式に ct\' をかけただけ)')

# (b) ------------------------------------------------------------------
print('\n(b) p_t = -E/P0')
a_x, a_y, a_s = sp.symbols('a_x a_y a_s', real=True)
Pi = sp.sqrt(tp**2 - xp**2 - yp**2 - hx**2)
mh = m * c / P0
L = -mh * Pi + a_x * xp + a_y * yp + hx * a_s
pt = sp.diff(L, tp)
gamma = 1 / sp.sqrt(1 - (xp**2 + yp**2 + hx**2) / tp**2)
E = gamma * m * c**2
print('  p_t = ∂L/∂(ct\')... ここでは t\' = d t/ds ではなく ct\' で微分: p_{ct} =', sp.simplify(pt))
# p_t (t' で微分) = c * ∂L/∂(ct')
pt_t = c * pt
print('  p_t = c ∂L/∂(ct\') と -E/P0 の差:', sp.simplify(pt_t + E / P0))
print('  Π = ct\'/γ との差:', sp.simplify(Pi - tp / gamma))
Pt = -E / P0
p2 = Pt**2 / c**2 - mh**2
Pmom2 = (E**2 - m**2 * c**4) / (c**2)   # P^2
print('  p_t^2/c^2 - (mc/P0)^2 - P^2/P0^2 =', sp.simplify(p2 - Pmom2 / P0**2))

# (c) ------------------------------------------------------------------
print('\n(c) 四極の p_x\' = -K x')
x, y, px, py, p, K = sp.symbols('x y px py p K', real=True)
a_s_q = -K * (x**2 - y**2) / 2
H = -(sp.sqrt(p**2 - px**2 - py**2) + a_s_q)    # ρ = ∞、a_x = a_y = 0
print('  p_x\' = -∂H/∂x =', sp.diff(-H, x), '、p_y\' = -∂H/∂y =', sp.diff(-H, y))

# (d) ------------------------------------------------------------------
print('\n(d) E/(P0 v0) = h0 h1/p0^2')
p0, dl = sp.symbols('p0 delta', positive=True)
h0 = sp.sqrt(1 + p0**2)
h1 = sp.sqrt(1 + p0**2 * (1 + dl)**2)
beta0 = p0 / h0
lhs = (h1 * m * c**2) / (p0 * m * c * beta0 * c)
print('  h1 mc^2/(p0 mc · β0 c) - h0 h1/p0^2 =', sp.simplify(lhs - h0 * h1 / p0**2))
