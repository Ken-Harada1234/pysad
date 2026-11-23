# -*- coding: utf-8 -*-
r"""04_quad.md の第 2 回の監査で足した途中式 (Q-01〜Q-12) を sympy・mpmath で確かめる。

(A) §7.4: 出口の符号 (g_out(u) = g_in(-u))、直線の立ち上がりと tanh の I_n、K(s) = K0 (s/ℓ)^2 の I_n、出口の z。
(B) §8 (4.36): x_f, q_xf, y_f, q_yf, z の w 微分と、Fortran の trans1(1:5,6) の形。
(C) §8 (4.37): 正準運動量での δ 微分の合成。
(D) §8 の z の行: D の段の (5,2)(5,4)(5,6) と (4.38) の一致。
(E) §4.1: solution-H2 の B_z = 0 への縮約。
(F) §4.2: K < 0 の回転。
(G) §7.2: p_z -> p の相対誤差と力の式の符号。
(H) §6.2: 交換子の次数。
"""
import mpmath as mp
import sympy as sp
import warnings
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
warnings.filterwarnings('ignore')
mp.mp.dps = 40
rng = np.random.default_rng(7)


def rnd(n):
    return [mp.mpf(float(v)) for v in rng.uniform(0.3, 1.3, n)]


# ---------------------------------------------------------------- (A)
print('== (A) §7.4 の I_n (硬い端との差のモーメント) ==')
s, l, s0, u, d = sp.symbols('s ell s0 u d', positive=True)
# 直線の立ち上がり g = s/ell (0..ell)
I0 = sp.integrate(s / l, (s, 0, l)) - (l - s0)
print('  直線: I0 =', sp.simplify(I0), '→ s0 =', sp.solve(I0, s0))
s0v = l / 2
I1 = sp.integrate((s - s0) * s / l, (s, 0, l)) - sp.integrate(s - s0, (s, s0, l))
I2 = sp.integrate((s - s0)**2 * s / l, (s, 0, l)) - sp.integrate((s - s0)**2, (s, s0, l))
print('  直線 I1 (一般 s0) =', sp.expand(I1), '、s0 = ell/2 で', sp.simplify(I1.subs(s0, s0v)))
print('  直線 I2 (s0 = ell/2) =', sp.simplify(I2.subs(s0, s0v)))
print('  直線 alpha = -24 I1 =', sp.simplify(-24 * I1.subs(s0, s0v)), '→ F1 = ell、F2 = I2 =', sp.simplify(I2.subs(s0, s0v)))
# tanh: g = (1 + tanh(u/d))/2、s0 = 0
for dd in (mp.mpf('0.5'), mp.mpf('2')):
    gm = lambda x: (1 + mp.tanh(x / dd)) / 2
    dg = lambda x: gm(x) - (1 if x > 0 else 0)
    n0 = mp.quad(dg, [-mp.inf, 0]) + mp.quad(dg, [0, mp.inf])
    n1 = mp.quad(lambda x: x * dg(x), [-mp.inf, 0]) + mp.quad(lambda x: x * dg(x), [0, mp.inf])
    n2 = mp.quad(lambda x: x**2 * dg(x), [-mp.inf, 0]) + mp.quad(lambda x: x**2 * dg(x), [0, mp.inf])
    print('  tanh d = %s: I0 = %s、I1 = %s (-pi^2 d^2/24 = %s)、I2 = %s' % (dd, mp.nstr(n0, 5), mp.nstr(n1, 15), mp.nstr(-mp.pi**2 * dd**2 / 24, 15), mp.nstr(n2, 5)))
print('  ∫0^∞ t/(e^t+1) dt = %s、pi^2/12 = %s、Σ(-1)^(k-1)/k^2 = %s' % (mp.nstr(mp.quad(lambda t: t / (mp.e**t + 1), [0, mp.inf]), 15), mp.nstr(mp.pi**2 / 12, 15), mp.nstr(mp.nsum(lambda k: (-1)**(k - 1) / k**2, [1, mp.inf]), 15)))
# K(s) = K0 (s/ell)^2
s0q = 2 * l / 3
J0 = sp.integrate(s**2 / l**2, (s, 0, l)) - (l - s0)
print('  K0 (s/ell)^2: I0 = 0 → s0 =', sp.solve(J0, s0))
J1 = sp.integrate((s - s0) * s**2 / l**2, (s, 0, l)) - sp.integrate(s - s0, (s, s0, l))
J2 = sp.integrate((s - s0)**2 * s**2 / l**2, (s, 0, l)) - sp.integrate((s - s0)**2, (s, s0, l))
print('  I1 (s0 = 2 ell/3) =', sp.simplify(J1.subs(s0, s0q)), '、I2 =', sp.simplify(J2.subs(s0, s0q)))
print('  I2 の途中: ∫(s-s0)^2 s^2/ell^2 ds =', sp.simplify(sp.integrate((s - s0)**2 * s**2 / l**2, (s, 0, l)).subs(s0, s0q)), '、∫_{s0}^{ell}(s-s0)^2 ds =', sp.simplify(sp.integrate((s - s0)**2, (s, s0, l)).subs(s0, s0q)))
print('  I1 の途中: ∫(s-s0) s^2/ell^2 ds =', sp.simplify(sp.integrate((s - s0) * s**2 / l**2, (s, 0, l)).subs(s0, s0q)), '、∫_{s0}^{ell}(s-s0) ds =', sp.simplify(sp.integrate(s - s0, (s, s0, l)).subs(s0, s0q)))
print('  α = -24 I1 / ell^2 の平方根 =', sp.sqrt(-24 * J1.subs(s0, s0q) / l**2), '= %.5f' % float(sp.sqrt(sp.Rational(2, 3))))

# 出口の符号: 立ち上がり K0 (s/ell)^2 を逆向きにした立ち下がり g_out(s) = ((ell - s)/ell)^2 (0..ell)、硬い端は s < s0o で K0
ell = 0.1
K0 = 1e-4
p = 1.2
s0o = ell / 3                      # 入口の s0 = 2ell/3 を s → ell - s で写した位置
gin = lambda x: 0.0 if x < 0 else (1.0 if x > ell else (x / ell)**2)
gout = lambda x: 0.0 if x > ell else (1.0 if x < 0 else ((ell - x) / ell)**2)
from scipy.integrate import quad
Iin = [quad(lambda x: (x - 2 * ell / 3)**n * (gin(x) - (1.0 if x > 2 * ell / 3 else 0.0)), 0, ell, points=[2 * ell / 3], epsabs=1e-16, epsrel=1e-13)[0] for n in range(3)]
Iout = [quad(lambda x: (x - s0o)**n * (gout(x) - (1.0 if x < s0o else 0.0)), 0, ell, points=[s0o], epsabs=1e-16, epsrel=1e-13)[0] for n in range(3)]
print('  K0 (s/ell)^2 の入口 I_n =', ['%.6e' % v for v in Iin], '、鏡像の出口 I_n =', ['%.6e' % v for v in Iout], '→ (-1)^n 倍の確認: I0, -I1, I2 =', ['%.6e' % v for v in (Iin[0], -Iin[1], Iin[2])])
# 出口の F = [[1-a, b],[0, 1+a]] を数値積分と比較 (K の 1 次なので K0 は小さく取る)
S1 = 0.2
def prof_out(x):
    return K0 * gout(x)
def rhs(x, v):
    M = v.reshape(2, 2)
    k = prof_out(x) / p
    return np.array([[M[1, 0], M[1, 1]], [-k * M[0, 0], -k * M[0, 1]]]).ravel()
Mx = solve_ivp(rhs, (-S1, ell + S1), np.eye(2).ravel(), method='DOP853', rtol=1e-13, atol=1e-17, max_step=0.01).y[:, -1].reshape(2, 2)
# 硬い端: s < s0o で K0。本体 (長さ A = S1 + s0o) の後に F、ドリフト (ell - s0o + S1)
w = np.sqrt(K0 / p)
A_, B_ = S1 + s0o, S1 + ell - s0o
body = lambda L_: np.array([[np.cos(w * L_), np.sin(w * L_) / w], [-w * np.sin(w * L_), np.cos(w * L_)]])
dr = lambda L_: np.array([[1, L_], [0, 1]])
a_o = K0 * Iout[1] / p
b_o = K0 * Iout[2] / p
Fo2 = np.array([[np.exp(a_o), b_o], [0, np.exp(-a_o)]])      # F = [[1 + a_o, b],[0, 1 - a_o]] の指数形
for nm, F in (('F = [[e^{a_o}, b],[0, e^{-a_o}]] (a_o = K0 I1^out/p = -K0 I1^in/p = -a)', Fo2),):
    m = dr(B_) @ F @ body(A_)
    m0 = dr(B_) @ body(A_)
    print('  出口の数値積分との差: フリンジなし %.2e、%s: %.2e' % (np.abs(Mx - m0).max(), nm, np.abs(Mx - m).max()))
print('  (入口の a_in = K0 I1_in/p = %.4e、出口の a_o = K0 I1_out/p = %.4e、b: %.4e / %.4e)' % (K0 * Iin[1] / p, a_o, K0 * Iin[2] / p, b_o))

# 出口の z: F2 = e^{-a} x1 p2 + (1/2) b' e^{-a} p2^2 (b' = b/p。a ∝ 1/p、b' ∝ 1/p^2)
x1, p2, pp, a1, b1, Q2 = sp.symbols('x1 p2 pp a1 b1 Q2')
aa, bb = a1 / pp, b1 / pp**2
F2o = sp.exp(-aa) * x1 * p2 + sp.Rational(1, 2) * bb * sp.exp(-aa) * p2**2
dz_q = sp.simplify(sp.diff(F2o, pp).subs(p2, pp * Q2))
# 出口: p_x1 = ∂F2/∂x1 = e^{-a} p_x2 より q_x1 = e^{-a} q_x2。本文の式: z2 - z1 = (a x1 - b (1 - a/2) q_x2) q_x1 (b = p b' = b1/pp)
zz = (aa * x1 - (b1 / pp) * (1 - aa / 2) * Q2) * (sp.exp(-aa) * Q2)
print('  出口の z (x 平面): ∂F2/∂p と (a x1 - b (1 - a/2) q_x2) q_x1 の差 =', sp.simplify(dz_q - zz))
# 出口の x2 = ∂F2/∂p2 と p_x1
p_x1 = sp.diff(F2o, x1)
x2_ = sp.diff(F2o, p2)
print('  出口の x2 - (e^{-a} x1 + b q_x1) (q_x1 = p_x1/p、b = b1/pp) =', sp.simplify(x2_ - (sp.exp(-aa) * x1 + (b1 / pp) * p_x1 / pp)))

# ---------------------------------------------------------------- (B)
print('\n== (B) §8 (4.36): w 微分 ==')
x, qx, y, qy, w_, l_ = sp.symbols('x qx y qy w ell_')
ph = w_ * l_
xf = x * sp.cos(ph) + qx * sp.sin(ph) / w_
qxf = -x * w_ * sp.sin(ph) + qx * sp.cos(ph)
yf = y * sp.cosh(ph) + qy * sp.sinh(ph) / w_
qyf = y * w_ * sp.sinh(ph) + qy * sp.cosh(ph)
dc1 = sp.cos(ph) - 1
dch2 = sp.cosh(ph) - 1
xs1 = ph - sp.sin(ph)
xsh2 = ph - sp.sinh(ph)
s1, sh2, c1, ch2 = sp.sin(ph), sp.sinh(ph), sp.cos(ph), sp.cosh(ph)
vals = dict(zip((x, qx, y, qy, w_, l_), (0.7, 0.4, -0.5, 0.3, 1.3, 0.6)))
def num(e):
    return float(sp.N(e.subs(vals), 30))
f_xf = -l_ * x * s1 + qx * (ph * dc1 + xs1) / w_**2
f_qxf = -x * (ph * c1 + s1) - l_ * qx * s1
f_yf = l_ * y * sh2 + qy * (ph * dch2 + xsh2) / w_**2
f_qyf = y * (ph * ch2 + sh2) + l_ * qy * sh2
for nm, e, f in (('x_f', xf, f_xf), ('q_xf', qxf, f_qxf), ('y_f', yf, f_yf), ('q_yf', qyf, f_qyf)):
    print('  ∂%s/∂w: sympy の微分と手で書いた式の差 %.1e' % (nm, abs(num(sp.diff(e, w_)) - num(f))))
print('  φ cosφ - sinφ = φ dc1 + xs1:', abs(num(ph * c1 - s1) - num(ph * dc1 + xs1)))
print('  φ coshφ - sinhφ = φ dch2 + xsh2 (xsh2 = φ - sinhφ):', abs(num(ph * ch2 - sh2) - num(ph * dch2 + xsh2)))
# z: (4.13) と (4.14)
dv = sp.symbols('dv')
Cx = qx**2 + w_**2 * x**2
Cy = qy**2 - w_**2 * y**2
z13 = -dv * l_ - sp.Rational(1, 4) * ((Cx + Cy) * l_ + (xf * qxf - x * qx) + (yf * qyf - y * qy))
u1 = xf - x
v1 = yf - y
z14 = -sp.Rational(1, 4) * (w_ * (x**2 * xs1 - y**2 * xsh2) + (qx**2 + qy**2) * l_ + u1 * qxf + x * qx * dc1 + v1 * qyf + y * qy * dch2) - dv * l_
vals[dv] = 0.01
print('  (4.13) と (4.14) の差 %.1e' % abs(num(z13) - num(z14)))
# ∂z/∂w: (4.14) の各項の w 微分 (本文の表)
dz_w = sp.diff(z14, w_)
terms = (x**2 * (xs1 - ph * dc1) - y**2 * (xsh2 - ph * dch2) + l_ * (-x * qx * s1 + y * qy * sh2)
         + sp.diff(xf, w_) * qxf + u1 * sp.diff(qxf, w_) + sp.diff(yf, w_) * qyf + v1 * sp.diff(qyf, w_))
print('  ∂/∂w[w x^2 xs1] = x^2 (xs1 - φ dc1):', abs(num(sp.diff(w_ * x**2 * xs1, w_)) - num(x**2 * (xs1 - ph * dc1))))
print('  ∂/∂w[-w y^2 xsh2] = -y^2 (xsh2 - φ dch2):', abs(num(sp.diff(-w_ * y**2 * xsh2, w_)) - num(-y**2 * (xsh2 - ph * dch2))))
print('  ∂/∂w[x qx dc1] = -ℓ x qx s1、∂/∂w[y qy dch2] = +ℓ y qy sh2:', abs(num(sp.diff(x * qx * dc1, w_)) - num(-l_ * x * qx * s1)), abs(num(sp.diff(y * qy * dch2, w_)) - num(l_ * y * qy * sh2)))
print('  ∂z/∂w = -(1/4){...} (Fortran trans1(5,6) の括弧) の差 %.1e' % abs(num(dz_w) - num(-sp.Rational(1, 4) * terms)))
# (4.13) からの導出 (別経路)
dz13 = sp.diff(z13, w_)
alt = -sp.Rational(1, 4) * (2 * w_ * (x**2 - y**2) * l_ + sp.diff(xf, w_) * qxf + xf * sp.diff(qxf, w_) + sp.diff(yf, w_) * qyf + yf * sp.diff(qyf, w_))
print('  (4.13) の w 微分 (別経路) と (4.14) の w 微分の差 %.1e' % abs(num(dz13) - num(dz_w)), '、別経路の式との差 %.1e' % abs(num(dz13) - num(alt)))
print('  q_xf = q_x + u2 (Fortran の u2 + pxi):', abs(num(qxf - (qx + (-x * w_ * s1 + qx * dc1)))))

# ---------------------------------------------------------------- (C)
print('\n== (C) (4.37): 正準運動量での δ 微分 ==')
dl, K, px_, py_ = sp.symbols('delta K px py')
pr = 1 + dl
wd = sp.sqrt(K / pr)
phd = wd * l_
qx_c = px_ / pr
qy_c = py_ / pr
xf_c = x * sp.cos(phd) + qx_c * sp.sin(phd) / wd
pxf_c = pr * (-x * wd * sp.sin(phd) + qx_c * sp.cos(phd))
yf_c = y * sp.cosh(phd) + qy_c * sp.sinh(phd) / wd
vv = {x: 0.7e-3, px_: 0.4e-3, y: -0.5e-3, py_: 0.3e-3, K: 2.0, l_: 0.3, dl: 0.01}
nv = lambda e: sp.N(e.subs(vv), 30)
wv = nv(wd)
phv = nv(phd)
qxv = nv(qx_c)
qyv = nv(qy_c)
prv = nv(pr)
w_p = -wv / (2 * prv)             # ∂w/∂δ
# 手で書いた式 ((4.36)(4.37))
sx, cx_ = sp.sin(phv), sp.cos(phv)
xv, yv = vv[x], vv[y]
lv = vv[l_]
dxf_dq = sx / wv
dxf_dw = -lv * xv * sx + qxv * (phv * cx_ - sx) / wv**2
hand_x = dxf_dw * w_p - qxv / prv * dxf_dq
print('  ∂x_f/∂δ|_{p_x}: 直接微分と (4.37) の差 %.1e' % abs(float(nv(sp.diff(xf_c, dl)) - hand_x)))
dqxf_dw = -xv * (phv * cx_ + sx) - lv * qxv * sx
dqxf_dq = cx_
hand_px = qxv * 0 + (-xv * wv * sx) + prv * dqxf_dw * w_p    # q_xf - q_x cosφ = -x w sinφ
print('  ∂p_xf/∂δ|_{p_x}: 直接微分と (4.37)+積の微分の差 %.1e' % abs(float(nv(sp.diff(pxf_c, dl)) - hand_px)))
q_xf_v = nv(-x * wd * sp.sin(phd) + qx_c * sp.cos(phd))
print('  q_xf = %.6e、q_x cosφ + (-x w sinφ) = %.6e' % (float(q_xf_v), float(qxv * cx_ - xv * wv * sx)))
sh_, ch_ = sp.sinh(phv), sp.cosh(phv)
dyf_dw = lv * yv * sh_ + qyv * (phv * ch_ - sh_) / wv**2
hand_y = dyf_dw * w_p - qyv / prv * sh_ / wv
print('  ∂y_f/∂δ|_{p_y}: 直接微分と (4.37) の差 %.1e' % abs(float(nv(sp.diff(yf_c, dl)) - hand_y)))
# 一般の合成: F(x, q_x, q_y, δ) の任意関数で (4.37) を確認
Fg = sp.Function('F')
Fe = x**2 * sp.cos(2 * qx_c) * sp.exp(qy_c * dl) * (1 + dl)
Fq = lambda: None
qx_s, qy_s = sp.symbols('qx_s qy_s')
Fe_q = x**2 * sp.cos(2 * qx_s) * sp.exp(qy_s * dl) * (1 + dl)
lhs = sp.diff(Fe_q.subs({qx_s: px_ / pr, qy_s: py_ / pr}), dl)
rhs = (sp.diff(Fe_q, dl) - qx_s / pr * sp.diff(Fe_q, qx_s) - qy_s / pr * sp.diff(Fe_q, qy_s)).subs({qx_s: px_ / pr, qy_s: py_ / pr})
print('  一般の F での (4.37) の差 %.1e' % abs(float((lhs - rhs).subs(vv))))

# ---------------------------------------------------------------- (D)
print('\n== (D) D の段の z の行 ==')
Px, Py, P, ell3 = sp.symbols('Px Py P ell3', positive=True)
pz = sp.sqrt(P**2 - Px**2 - Py**2)
r = ell3 * (1 / pz - 1 / P)
Zd = ell3 * (1 - P / pz + (Px**2 + Py**2) / (2 * P**2))        # Δz = ℓ z' (4.15)
vv2 = {Px: 0.3, Py: 0.2, P: 1.1, ell3: 0.5}
nv2 = lambda e: float(sp.N(e.subs(vv2), 30))
print('  ∂Δz/∂p_x = p_x ∂r/∂p = ∂x/∂δ (T16):', abs(nv2(sp.diff(Zd, Px)) - nv2(Px * sp.diff(r, P))))
print('  ∂Δz/∂p_y = p_y ∂r/∂p:', abs(nv2(sp.diff(Zd, Py)) - nv2(Py * sp.diff(r, P))))
print('  ∂r/∂p = ℓ(1/p^2 - p/p_z^3):', abs(nv2(sp.diff(r, P)) - nv2(ell3 * (1 / P**2 - P / pz**3))))
print('  ∂Δz/∂p = ℓ a (1/p_z^3 - 1/p^3) (a = p_x^2 + p_y^2):', abs(nv2(sp.diff(Zd, P)) - nv2(ell3 * (Px**2 + Py**2) * (1 / pz**3 - 1 / P**3))))
aa_ = Px**2 + Py**2
gg = (pz - P) * (1 + pz / P + pz**2 / P**2)
ff = ell3 / pz**3
print('  Fortran の g f p_x = T16、-a g f/p = T56 (D):', abs(nv2(gg * ff * Px - Px * sp.diff(r, P))), abs(nv2(-aa_ * gg * ff / P - sp.diff(Zd, P))))
# (4.38) を D の段で: T = [[1,r+..,...]]
Tm = sp.zeros(6, 6)
# x = x0 + r Px; 微分行列 (5 行は求める対象なので 0)
Tm[0, 0] = 1
Tm[0, 1] = sp.diff(r * Px, Px)
Tm[0, 3] = sp.diff(r * Px, Py)
Tm[0, 5] = sp.diff(r * Px, P)
Tm[1, 1] = 1
Tm[2, 2] = 1
Tm[2, 1] = sp.diff(r * Py, Px)
Tm[2, 3] = sp.diff(r * Py, Py)
Tm[2, 5] = sp.diff(r * Py, P)
Tm[3, 3] = 1
Tm[5, 5] = 1
# 最終段の運動量は p_x, p_y 不変、δ も不変なので T26 = T46 = 0
for j in range(4):
    t5j = -Tm[0, j] * Tm[1, 5] + Tm[1, j] * Tm[0, 5] - Tm[2, j] * Tm[3, 5] + Tm[3, j] * Tm[2, 5]
    ref = [sp.Integer(0), sp.diff(Zd, Px), sp.Integer(0), sp.diff(Zd, Py)][j]
    print('  (4.38) の T5%d = %.6e、直接微分 %.6e' % (j + 1, nv2(t5j), nv2(ref)))

# ---------------------------------------------------------------- (E)
print('\n== (E) solution-H2 の B_z = 0 への縮約 ==')
Kk, pp_, x0, qx0, y0, qy0 = sp.symbols('K p x0 qx0 y0 qy0', positive=True)
w1 = w2 = sp.sqrt(Kk / pp_)
bp = 0
wp = w1 + w2
wdd = bp / wp
V = w1**2 + w2**2
a_h = (w2 * wp * x0 - bp * qy0) / V
bw_h = (wp * qx0 - bp * w2 * y0) / V
c_h = (w1 * wdd * x0 + qy0) / V
dw_h = (-wdd * qx0 + w1 * y0) / V
print('  V =', sp.simplify(V), '、U = 1/V =', sp.simplify(1 / V), '、w_+ =', wp)
print('  a =', sp.simplify(a_h), '、bw =', sp.simplify(bw_h), '→ b = w1 bw =', sp.simplify(w1 * bw_h))
print('  c =', sp.simplify(c_h), '、dw =', sp.simplify(dw_h), '→ d = w2 dw =', sp.simplify(w2 * dw_h))
phi_ = sp.symbols('phi')
dx_h = a_h * (sp.cos(phi_) - 1) + bw_h * sp.sin(phi_)
dy_h = wp * (c_h * sp.sinh(phi_) + dw_h * (sp.cosh(phi_) - 1))
print('  Δx = a dc1 + bw s1 =', sp.simplify(dx_h), '、Δy = w_+ (c sh2 + dw dch2) =', sp.simplify(dy_h))
print('  c - p q_y0/(2K) = c - p_y0/(2K):', sp.simplify(c_h - (qy0 * pp_) / (2 * Kk)), '(p_y0 = p q_y0 を入れる)')

# ---------------------------------------------------------------- (F)
print('\n== (F) K < 0: 90 度の回転 ==')
X, PX, Y, PY = sp.symbols('X PX Y PY')
wq, ll = sp.symbols('wq ll', positive=True)
xx0, qxx0, yy0, qyy0 = sp.symbols('xx0 qxx0 yy0 qyy0')
# 回転後 (X, QX, Y, QY) = (y, qy, -x, -qx)、|K| の集束: X は cos/sin、Y は cosh/sinh
Xf = yy0 * sp.cos(wq * ll) + qyy0 * sp.sin(wq * ll) / wq
QXf = -yy0 * wq * sp.sin(wq * ll) + qyy0 * sp.cos(wq * ll)
Yin, QYin = -xx0, -qxx0
Yf = Yin * sp.cosh(wq * ll) + QYin * sp.sinh(wq * ll) / wq
QYf = Yin * wq * sp.sinh(wq * ll) + QYin * sp.cosh(wq * ll)
xf_ = -Yf
qxf_ = -QYf
yf_ = Xf
qyf_ = QXf
print('  x_f =', sp.simplify(xf_), '、q_xf =', sp.simplify(qxf_))
print('  y_f =', sp.simplify(yf_), '、q_yf =', sp.simplify(qyf_))

# ---------------------------------------------------------------- (G)
print('\n== (G) §7.2: p_z -> p の近似と力の式 ==')
qv = sp.symbols('q', positive=True)
print('  1/sqrt(1 - q^2) - 1 =', sp.series(1 / sp.sqrt(1 - qv**2) - 1, qv, 0, 5))
for qq in (1e-3, 1e-2, 3e-2):
    print('  q = %.0e mrad: 相対誤差 q^2/2 = %.2e (厳密 1/ζ - 1 = %.4e)' % (qq * 1e3, qq**2 / 2, 1 / np.sqrt(1 - qq**2) - 1))
bx, by, bs, p_x, p_y, p_z = sp.symbols('bx by bs p_x p_y p_z')
cr = sp.Matrix([p_x, p_y, p_z]).cross(sp.Matrix([bx, by, bs]))
print('  (p × b)_x =', cr[0], '、(p × b)_y =', cr[1])
Ksym, xs_, ys_ = sp.symbols('K x y')
print('  b_y = Kx, b_s = 0: dp_x/ds =', sp.simplify((cr[0] / p_z).subs({by: Ksym * xs_, bs: 0, bx: Ksym * ys_})), '、b_x = Ky: dp_y/ds =', sp.simplify((cr[1] / p_z).subs({by: Ksym * xs_, bs: 0, bx: Ksym * ys_})))

# ---------------------------------------------------------------- (H)
print('\n== (H) §6.2: 交換子の次数 ==')
Lc, C3 = sp.symbols('L C3', commutative=False)
h, k_i, k_j = sp.symbols('h k_i k_j')
Xi = k_i * h * Lc + k_i**3 * h**3 * C3
Xj = k_j * h * Lc + k_j**3 * h**3 * C3
cm = sp.expand(Xi * Xj - Xj * Xi)
print('  [X_i, X_j] の h^2 の係数 =', sp.factor(cm.coeff(h, 2)), '(L と L の交換子は 0)')
claim = k_i * k_j * (k_j**2 - k_i**2) * (Lc * C3 - C3 * Lc)
print('  h^4 の係数 - k_i k_j (k_j^2 - k_i^2)[L, C3] =', sp.expand(cm.coeff(h, 4) - claim), '、h^6 以上は X_i の h^5 以上の項から (C3 どうしは可換)')

# ---------------------------------------------------------------- (C2)
print('\n== (C2) δ の列の閉じた形 (正準運動量、p_x を固定) ==')
pyf_c = pr * (y * wd * sp.sinh(phd) + qy_c * sp.cosh(phd))
phi_v = phv
w_v = wv
p_v = prv
closed_x = (phi_v * xv * sx - qxv * (phi_v * cx_ + sx) / w_v) / (2 * p_v)
closed_y = -(phi_v * yv * sh_ + qyv * (phi_v * ch_ + sh_) / w_v) / (2 * p_v)
closed_px = w_v / 2 * (xv * (phi_v * cx_ - sx) + lv * qxv * sx)
closed_py = w_v / 2 * (yv * (sh_ - phi_v * ch_) - lv * qyv * sh_)
print('  ∂x_f/∂δ = [φ x sinφ - q_x (φ cosφ + sinφ)/w]/(2p): 差 %.1e' % abs(float(nv(sp.diff(xf_c, dl)) - closed_x)))
print('  ∂y_f/∂δ = -[φ y sinhφ + q_y (φ coshφ + sinhφ)/w]/(2p): 差 %.1e' % abs(float(nv(sp.diff(yf_c, dl)) - closed_y)))
print('  ∂p_xf/∂δ = (w/2)[x (φ cosφ - sinφ) + ℓ q_x sinφ]: 差 %.1e' % abs(float(nv(sp.diff(pxf_c, dl)) - closed_px)))
print('  ∂p_yf/∂δ = (w/2)[y (sinhφ - φ coshφ) - ℓ q_y sinhφ]: 差 %.1e' % abs(float(nv(sp.diff(pyf_c, dl)) - closed_py)))
# 運動量 p_x の行の他の成分
print('  ∂p_xf/∂x = p ∂q_xf/∂x = -p w sinφ:', abs(float(nv(sp.diff(pxf_c, x)) - (-p_v * w_v * sx))))
print('  ∂p_xf/∂p_x = cosφ:', abs(float(nv(sp.diff(pxf_c, px_)) - cx_)))
