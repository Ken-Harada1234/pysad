# -*- coding: utf-8 -*-
r"""04_quad.md の途中式を sympy・数値で確かめる (本文の「展開」の各行に対応)。

(A) Strang 分割の BCH の係数 (4.19) と、Yoshida の合成 (4.20)〜(4.22)。
(B) 運動学的ドリフト (4.15)〜(4.17) の展開。
(C) H2 の解 (4.9) と z の更新 (4.12)〜(4.14) の恒等式、(4.36)。
(D) 非線形フリンジ: 力の積分 (4.26)、生成関数 (4.27) の符号 (C_poisson_lie.md の約束)、母関数 (4.28) から点変換 (4.29)・(4.30)・z。
(E) 線形フリンジ: 摂動の積分 (4.32) (I0 != 0 の一般形)、端の位置をずらしたときの I_n の変換、F2 型母関数から z (4.33)。
(F) 転送行列 (4.35)〜(4.38)。
"""
import numpy as np
import sympy as sp
from scipy.linalg import expm, logm
from scipy.integrate import quad
import _common as C

C.start(__file__)

# ------------------------------------------------------------------ (A)
print('== (A) BCH: S2(h) = exp(hD/2) exp(hB) exp(hD/2) の log の h^3 の係数 ==')
rng = np.random.default_rng(1)
n = 4
Bm, Dm = rng.normal(size=(n, n)), rng.normal(size=(n, n))
com = lambda X, Y: X @ Y - Y @ X
c_dd = com(Dm, com(Dm, Bm))          # [D,[D,B]]
c_bb = com(Bm, com(Bm, Dm))          # [B,[B,D]]
res = []
for h in (0.04, 0.02, 0.01):
    S2 = expm(h * Dm / 2) @ expm(h * Bm) @ expm(h * Dm / 2)
    res.append((logm(S2) - h * (Bm + Dm)) / h**3)
L3 = (4 * res[2] - res[1]) / 3.0     # Richardson (h^2 の誤差を消す)
Amat = np.stack([c_dd.ravel(), c_bb.ravel()], 1)
sol, *_ = np.linalg.lstsq(Amat, L3.ravel(), rcond=None)
print('  h^3 の係数 (最小二乗): [D,[D,B]] の係数 = %.6f (-1/24 = %.6f)、[B,[B,D]] の係数 = %.6f (1/12 = %.6f)' % (sol[0], -1 / 24, sol[1], 1 / 12))
print('  最小二乗の残差 %.2e' % np.abs(Amat @ sol - L3.ravel()).max())

print('\n  Yoshida: S4 = S2(w1 h) S2(w0 h) S2(w1 h) の log の h^3 の係数が消えること')
w1 = 1 / (2 - 2**(1 / 3.0))
w0 = 1 - 2 * w1
print('  w1 = %.9f、w0 = %.9f、2 w1 + w0 - 1 = %.1e、2 w1^3 + w0^3 = %.1e、w1 + w0 = %.6f ((w1+w0)/2 = %.6f、w1/2 = %.6f)' % (
    w1, w0, 2 * w1 + w0 - 1, 2 * w1**3 + w0**3, w1 + w0, (w1 + w0) / 2, w1 / 2))
S2f = lambda h: expm(h * Dm / 2) @ expm(h * Bm) @ expm(h * Dm / 2)
lst = []
for h in (0.04, 0.02, 0.01):
    S4 = S2f(w1 * h) @ S2f(w0 * h) @ S2f(w1 * h)
    lst.append(np.abs(logm(S4) - h * (Bm + Dm)).max())
print('  max|log S4 - h(B+D)| (h = 0.04, 0.02, 0.01): %.3e, %.3e, %.3e  比 %.2f, %.2f (h^5 なら 32)' % (lst[0], lst[1], lst[2], lst[0] / lst[1], lst[1] / lst[2]))
h = 0.03
S4h, S4m = S2f(w1 * h) @ S2f(w0 * h) @ S2f(w1 * h), S2f(-w1 * h) @ S2f(-w0 * h) @ S2f(-w1 * h)
print('  S4(h) S4(-h) - 1 = %.2e (対称な積なので S4(-h) = S4(h)^-1、log S4 は h の奇数次だけ)' % np.abs(S4h @ S4m - np.eye(n)).max())

# ------------------------------------------------------------------ (B)
print('\n== (B) ΔH のドリフト (4.15)〜(4.17) ==')
ell = sp.symbols('ell', positive=True)
p, px, py = sp.symbols('p p_x p_y', positive=True)
pz = sp.sqrt(p**2 - px**2 - py**2)
dH = p - pz - (px**2 + py**2) / (2 * p)
print('  ∂ΔH/∂px - px(1/pz - 1/p) = %s、∂ΔH/∂py - py(1/pz - 1/p) = %s、∂ΔH/∂p - (1 - p/pz + (px^2 + py^2)/(2p^2)) = %s' % (
    sp.simplify(sp.diff(dH, px) - px * (1 / pz - 1 / p)), sp.simplify(sp.diff(dH, py) - py * (1 / pz - 1 / p)),
    sp.simplify(sp.diff(dH, p) - (1 - p / pz + (px**2 + py**2) / (2 * p**2)))))
zeta = sp.symbols('zeta', positive=True)
zexpr = 1 - 1 / zeta + (1 - zeta**2) / 2          # (4.15) の z' を ζ = pz/p、q^2 = 1 - ζ^2 で
print('  z\' = 1 - 1/ζ + (1 - ζ^2)/2 の通分: 分子 = %s、(ζ-1)^2(ζ+2) との差 %s' % (
    sp.expand(sp.simplify(zexpr * 2 * zeta)), sp.simplify(sp.expand(-(zeta - 1)**2 * (zeta + 2)) - sp.expand(sp.simplify(zexpr * 2 * zeta)))))
print('  -(2+ζ)(1-ζ)^2/(2ζ) との差 = %s' % sp.simplify(zexpr + (2 + zeta) * (1 - zeta)**2 / (2 * zeta)))
a = sp.symbols('a', positive=True)
zt = sp.sqrt(1 - a**2)
print('  1/ζ - 1 = %s' % sp.series(1 / zt - 1, a, 0, 6))
print('  Δx/ℓ = a (1/ζ - 1) = %s' % sp.series(a * (1 / zt - 1), a, 0, 6))
print('  Δz/ℓ = -(2+ζ)(1-ζ)^2/(2ζ) = %s' % sp.series(-(2 + zt) * (1 - zt)**2 / (2 * zt), a, 0, 7))

# ------------------------------------------------------------------ (C)
print('\n== (C) H2 の解と z の更新 (4.9)(4.12)〜(4.14)(4.36) ==')
s, w, x0, q0, y0, r0 = sp.symbols('s w x_0 q_0 y_0 r_0', positive=True)
xs = x0 * sp.cos(w * s) + q0 * sp.sin(w * s) / w
qs = sp.diff(xs, s)
print('  水平: x\'\' + w^2 x = %s、q = x\' = %s' % (sp.simplify(sp.diff(xs, s, 2) + w**2 * xs), sp.simplify(qs)))
ys = y0 * sp.cosh(w * s) + r0 * sp.sinh(w * s) / w
rs = sp.diff(ys, s)
print('  垂直: y\'\' - w^2 y = %s、q_y = y\' = %s' % (sp.simplify(sp.diff(ys, s, 2) - w**2 * ys), sp.simplify(rs)))
Cx = qs**2 + w**2 * xs**2
Cy = rs**2 - w**2 * ys**2
print('  d(q^2 + w^2 x^2)/ds = %s、d(q_y^2 - w^2 y^2)/ds = %s' % (sp.simplify(sp.diff(Cx, s)), sp.simplify(sp.diff(Cy, s))))
print('  d(x q)/ds - (q^2 - w^2 x^2) = %s、d(y q_y)/ds - (q_y^2 + w^2 y^2) = %s' % (
    sp.simplify(sp.diff(xs * qs, s) - (qs**2 - w**2 * xs**2)), sp.simplify(sp.diff(ys * rs, s) - (rs**2 + w**2 * ys**2))))
# (4.13): ∫(q_x^2 + q_y^2) ds = (1/2)[(C_x + C_y) ℓ + [x q] + [y q_y]]  (数値で検査)
vals = {w: 1.7, x0: 0.3, q0: -0.2, y0: 0.5, r0: 0.4, ell: 0.8}
lhs = sp.Integral(qs**2 + rs**2, (s, 0, ell)).subs(vals).evalf(20)
rhs = (sp.Rational(1, 2) * ((Cx + Cy).subs(s, 0) * ell + (xs * qs).subs(s, ell) - (xs * qs).subs(s, 0)
                           + (ys * rs).subs(s, ell) - (ys * rs).subs(s, 0))).subs(vals).evalf(20)
print('  ∫(q_x^2 + q_y^2) ds (数値積分) = %s、(4.13) の右辺 = %s、差 %s' % (lhs, rhs, lhs - rhs))
# (4.14) の括弧 = (4.13) の括弧 (水平・垂直を別々に)
ph = sp.symbols('varphi', positive=True)
xf, qf = xs.subs(s, ell), qs.subs(s, ell)
u1 = xf - x0
dc1, s1, xs1 = sp.cos(w * ell) - 1, sp.sin(w * ell), w * ell - sp.sin(w * ell)
lhsx = Cx.subs(s, 0) * ell + xf * qf - x0 * q0
rhsx = w * x0**2 * xs1 + q0**2 * ell + u1 * qf + x0 * q0 * dc1
print('  水平: C_x ℓ + [x q] - (w x^2 xs1 + q^2 ℓ + u1 q\' + x q dc1) = %s' % sp.simplify((lhsx - rhsx).rewrite(sp.exp)))
yf, rf = ys.subs(s, ell), rs.subs(s, ell)
v1 = yf - y0
dch, shh, xsh = sp.cosh(w * ell) - 1, sp.sinh(w * ell), w * ell - sp.sinh(w * ell)
lhsy = Cy.subs(s, 0) * ell + yf * rf - y0 * r0
rhsy = -w * y0**2 * xsh + r0**2 * ell + v1 * rf + y0 * r0 * dch
print('  垂直: C_y ℓ + [y q_y] - (-w y^2 xsh2 + q_y^2 ℓ + v1 q_y\' + y q_y dch2) = %s' % sp.simplify((lhsy - rhsy).rewrite(sp.exp)))
# 中間の恒等式: x_f q_f - x q - u1 q_f = x u2、 u2 = -x w sin + q dc
u2 = qf - q0
print('  x u2 - x(-x w s1 + q dc1) = %s、  y v2 - y(y w sh2 + q_y dch2) = %s' % (
    sp.simplify((x0 * u2 - x0 * (-x0 * w * s1 + q0 * dc1)).rewrite(sp.exp)), sp.simplify((y0 * (rf - r0) - y0 * (y0 * w * shh + r0 * dch)).rewrite(sp.exp))))
# (4.36)
K, pp = sp.symbols('K p', positive=True)
wd = sp.sqrt(K / pp)
xf_p = x0 * sp.cos(wd * ell) + q0 * sp.sin(wd * ell) / wd
wp = sp.diff(wd, pp)
mine = (-ell * x0 * sp.sin(wd * ell) + q0 * (wd * ell * sp.cos(wd * ell) - sp.sin(wd * ell)) / wd**2) * wp
print('  ∂x_f/∂p - (4.36) = %s、∂w/∂p + w/(2p) = %s' % (sp.simplify(sp.diff(xf_p, pp) - mine), sp.simplify(wp + wd / (2 * pp))))
print('  φ cos φ - sin φ = φ dc1 + xs1 の確認: %s' % sp.simplify(ph * sp.cos(ph) - sp.sin(ph) - (ph * (sp.cos(ph) - 1) + (ph - sp.sin(ph)))))

# ------------------------------------------------------------------ (D)
print('\n== (D) 非線形フリンジ ==')
x, y, K, p = sp.symbols('x y K p', real=True)
px, py = sp.symbols('p_x p_y', real=True)
Sx = sp.symbols('s', real=True)
Kf = sp.Function('K')(Sx)
V = Kf * x * y - sp.diff(Kf, Sx, 2) * x * y * (x**2 + y**2) / 12
lap = sp.expand(sp.diff(V, x, 2) + sp.diff(V, y, 2) + sp.diff(V, Sx, 2))
print('  ∇^2 V = %s (残りは K\'\'\'\' の項で x y の 5 次。捨てる)' % lap)
Kpp = sp.diff(Kf, Sx, 2)
bx_ = sp.diff(V, x)
by_ = sp.diff(V, y)
bs_ = sp.diff(V, Sx)
print('  b_x = ∂V/∂x = %s' % sp.expand(bx_))
print('  b_y = ∂V/∂y = %s' % sp.expand(by_))
print('  b_s = ∂V/∂s = %s (K\'\'\' の項は 5 次で捨てる → K\' x y)' % sp.expand(bs_))
# 力の積分
xp, yp = sp.symbols('xp yp')                    # x' = px/p、y' = py/p
g = (x**3 + 3 * x * y**2)
gd = sp.diff(g, x) * xp + sp.diff(g, y) * yp    # 軌道に沿った微分 g' (自由な運動)
dpx_Kpp = sp.expand(-K / 12 * gd)               # ∫(+K''/12) g ds (K'' = K δ') = -(K/12) g'(0)
dpx_bs = K * x * y * yp                         # ∫ (p_y/p_z) b_s ds = K x y p_y/p   (K' = K δ)
want = K * (-sp.Rational(1, 4) * (x**2 + y**2) * xp + sp.Rational(1, 2) * x * y * yp)
print('  Δp_x の K\'\' の項 = %s、b_s の項 = %s' % (sp.factor(dpx_Kpp), dpx_bs))
print('  合計 - (4.26) = %s' % sp.simplify(dpx_Kpp + dpx_bs - want))
hh2 = (3 * x**2 * y + y**3)
hd = sp.diff(hh2, x) * xp + sp.diff(hh2, y) * yp
dpy_Kpp = sp.expand(+K / 12 * hd)              # b_x の -K''/12 h: ∫(-K''/12) h ds = +(K/12) h'
dpy_bs = -K * x * y * xp                       # -(p_x/p_z) b_s
print('  Δp_y の K\'\' の項 = %s、b_s の項 = %s、合計 = %s' % (sp.factor(dpy_Kpp), dpy_bs, sp.factor(sp.expand(dpy_Kpp + dpy_bs))))
# 座標の階段: dpx/ds ⊃ (K/12) δ'(s) g → px ⊃ (K/12) δ(s) g → Δx = (K/12p) g
# 生成関数 f_in の符号 (C_poisson_lie.md: [f,g] = Σ(∂q f ∂p g - ∂p f ∂q g)、exp(:f:): x → x - ∂f/∂px、px → px + ∂f/∂x)
fin = -(K / (12 * p)) * (x**3 * px + 3 * x * y**2 * px - 3 * x**2 * y * py - y**3 * py)
def pb(f_, g_):
    return sum(sp.diff(f_, q_) * sp.diff(g_, p_) - sp.diff(f_, p_) * sp.diff(g_, q_) for q_, p_ in ((x, px), (y, py)))
print('  [f, x] + ∂f/∂px = %s、[f, px] - ∂f/∂x = %s (:f: x = -∂f/∂px、:f: px = +∂f/∂x)' % (sp.simplify(pb(fin, x) + sp.diff(fin, px)), sp.simplify(pb(fin, px) - sp.diff(fin, x))))
sq = {px: p * xp, py: p * yp}
print('  Δx = -∂f/∂px - (K/12p)(x^3 + 3 x y^2) = %s' % sp.simplify(-sp.diff(fin, px) - K / (12 * p) * g))
print('  Δy = -∂f/∂py + (K/12p)(3 x^2 y + y^3) = %s' % sp.simplify(-sp.diff(fin, py) + K / (12 * p) * hh2))
print('  Δp_x = +∂f/∂x - (4.26) = %s' % sp.simplify(sp.diff(fin, x).subs(sq) - (dpx_Kpp + dpx_bs)))
print('  Δp_y = +∂f/∂y - (K\'\' + b_s の合計) = %s' % sp.simplify(sp.diff(fin, y).subs(sq) - (dpy_Kpp + dpy_bs)))

# 母関数 (4.28)
x1, y1, px2, py2 = sp.symbols('x_1 y_1 p_x2 p_y2', real=True)
c = K / (12 * p)
dX = c * (x1**3 + 3 * x1 * y1**2)
dY = -c * (3 * x1**2 * y1 + y1**3)
aa, bb = K * x1**2 / (4 * p), K * y1**2 / (4 * p)
print('\n  Δx1 - x1(a/3 + b) = %s、Δy1 + y1(a + b/3) = %s' % (sp.simplify(dX - x1 * (aa / 3 + bb)), sp.simplify(dY + y1 * (aa + bb / 3))))
H0 = px2 * dX + py2 * dY
G = H0 + sp.Rational(1, 2) * (sp.diff(H0, x1) * sp.diff(H0, px2) + sp.diff(H0, y1) * sp.diff(H0, py2))
print('  ∂Δx1/∂x1 = %s = a + b: %s、∂Δx1/∂y1 = %s = 2 b x1 y1... → 6c x1 y1 = %s' % (
    sp.factor(sp.diff(dX, x1)), sp.simplify(sp.diff(dX, x1) - (aa + bb)), sp.factor(sp.diff(dX, y1)), sp.factor(6 * c * x1 * y1)))
V2x = sp.Rational(1, 2) * (sp.diff(dX, x1) * dX + sp.diff(dX, y1) * dY)
V2y = sp.Rational(1, 2) * (sp.diff(dY, x1) * dX + sp.diff(dY, y1) * dY)
tt = (aa - bb)**2 / 6
print('  V2x - t x1 = %s、V2y - t y1 = %s (t = (a - b)^2/6)' % (sp.simplify(V2x - tt * x1), sp.simplify(V2y - tt * y1)))
print('  (括弧の中の整理) V2x/x1 = (1/2)[(a+b)(a/3+b) - 2b(a+b/3)] = %s' % sp.expand(sp.Rational(1, 2) * ((aa + bb) * (aa / 3 + bb) - 2 * bb * (aa + bb / 3))))
Phi_x = x1 + dX + V2x
Phi_y = y1 + dY + V2y
Jm = sp.Matrix([[sp.diff(Phi_x, x1), sp.diff(Phi_x, y1)], [sp.diff(Phi_y, x1), sp.diff(Phi_y, y1)]])
d_ = aa + bb
ab_ = aa - bb
xx = 1 + d_ + ab_ * (5 * aa - bb) / 6
yy = 1 - d_ + ab_ * (aa - 5 * bb) / 6
hh = K / (2 * p) * x1 * y1 * (1 - ab_ / 3)
print('  x_x - (1 + d + ab(5a - b)/6) = %s、y_y - (1 - d + ab(a - 5b)/6) = %s' % (sp.simplify(Jm[0, 0] - xx), sp.simplify(Jm[1, 1] - yy)))
print('  ∂x2/∂y1 - h = %s、∂y2/∂x1 + h = %s' % (sp.simplify(Jm[0, 1] - hh), sp.simplify(Jm[1, 0] + hh)))
print('  x2 = x1 + ∂G/∂px2 - x1(1 + a/3 + b + t) = %s、y2 - y1(1 - a - b/3 + t) = %s' % (
    sp.simplify(x1 + sp.diff(G, px2) - x1 * (1 + aa / 3 + bb + tt)), sp.simplify(y1 + sp.diff(G, py2) - y1 * (1 - aa - bb / 3 + tt))))
Fp = sp.Matrix([sp.diff(G, x1) + px2, sp.diff(G, y1) + py2])
print('  ∂F/∂(x1, y1) - J^T (px2, py2) = %s (G は p2 に線形)' % sp.simplify(Fp - Jm.T * sp.Matrix([px2, py2])))
dGdp = sp.diff(G, p)
zform = -(px2 * (dX + 2 * tt * x1) + py2 * (dY + 2 * tt * y1)) / p
print('  ∂G/∂p - (-[px2(Δx1 + 2 t x1) + py2(Δy1 + 2 t y1)]/p) = %s、(Δx1 + 2 t x1 = (x2 - x1) + t x1)' % sp.simplify(dGdp - zform))
print('  逆行列: J^-T - (1/f)[[y_y, h], [-h, x_x]] (f = x_x y_y + h^2) = %s' % sp.simplify(
    Jm.T.inv() - sp.Matrix([[Jm[1, 1], -Jm[1, 0]], [-Jm[0, 1], Jm[0, 0]]]) / (Jm[0, 0] * Jm[1, 1] - Jm[0, 1] * Jm[1, 0])))

# ------------------------------------------------------------------ (E)
print('\n== (E) 線形フリンジ ==')
u, A, B, k0, p = sp.symbols('u A B k_0 p')
I0s, I1s, I2s = sp.symbols('I_0 I_1 I_2')
D = lambda l_: sp.Matrix([[1, l_], [0, 1]])
N = sp.Matrix([[0, 0], [1, 0]])
integrand = (-(D(B - u) * N * D(u + A)) * k0 / p).applyfunc(sp.expand)


def mom(e):
    e = sp.Poly(sp.expand(e), u)
    return sum(co * [I0s, I1s, I2s][m_[0]] for m_, co in e.terms())


print('  被積分関数 -D(B-u) N D(u+A) k0/p =')
sp.pprint(integrand)
dM = integrand.applyfunc(mom)
print('  δM (一般の I0) = ∫ ... δg du (u^n → I_n):')
sp.pprint(dM.applyfunc(sp.factor))
a_s, b_s = sp.symbols('a b')
dM0 = dM.subs(I0s, 0)
Fm = sp.Matrix([[1 + a_s, b_s], [0, 1 - a_s]])
ref = D(B) * (Fm - sp.eye(2)) * D(A)
print('  D(B)(F - 1)D(A) =')
sp.pprint(ref)
solab = sp.solve([dM0[0, 0] - ref[0, 0], dM0[0, 1] - ref[0, 1]], [a_s, b_s])
print('  I0 = 0: a = %s、b = %s。残りの成分の差: (1,1) %s、(1,0) %s' % (
    solab[a_s], sp.simplify(solab[b_s]), sp.simplify(dM0[1, 1] - ref[1, 1].subs(solab)), sp.simplify(dM0[1, 0] - ref[1, 0])))
print('  I0 の項: δM - δM(I0=0) + (k0 I0/p) D(B) N D(A) = %s' % sp.simplify(dM - dM0 + k0 * I0s / p * D(B) * N * D(A)))
Dl, uu = sp.symbols('Delta uu')
J0 = I0s + Dl
J1 = I1s - Dl * I0s + sp.integrate((uu - Dl), (uu, 0, Dl))
J2 = I2s - 2 * Dl * I1s + Dl**2 * I0s + sp.integrate((uu - Dl)**2, (uu, 0, Dl))
print('  端を Δ ずらしたときの I_n: I0\' = %s、I1\' = %s、I2\' = %s' % (J0, sp.expand(J1), sp.expand(J2)))
shf = {Dl: -I0s}
print('  I0\' = 0 になる Δ = -I0 のとき: I1\' = %s、I2\' = %s' % (sp.expand(J1.subs(shf)), sp.expand(J2.subs(shf))))
# 数値例
ell0 = 0.1
prof = lambda s_: 0.0 if s_ < 0 else (1.0 if s_ > ell0 else (s_ / ell0)**2)


def moments(s0):
    return [quad(lambda s_: (s_ - s0)**m * (prof(s_) - (1.0 if s_ > s0 else 0.0)), 0, ell0, points=[s0], epsabs=1e-16, epsrel=1e-13)[0] for m in range(3)]


s00 = 2 * ell0 / 3
mom0 = moments(s00)
print('  数値例 K(s) = K0 (s/l)^2、l = 0.1: I0 = 0 の位置 s00 = 2l/3 = %.6f。そこでの F1|F1| = -24 I1 = %.6e、F2 = I2 = %.6e' % (s00, -24 * mom0[1], mom0[2]))
for dss in (0.0, 0.01, -0.01):
    I0, I1, I2 = moments(s00 + dss)
    print('   s0 = s00 %+.2f: I0 = %+.4e、I1 = %.6e、I2 = %.6e。 -24(I1 + I0^2/2) = %.6e、I2 + 2 I0 I1 + 2 I0^3/3 = %.6e | ヘルプの形 24(I0^2/2 - I1) = %.6e、I2 - I0^3/3 = %.6e' % (
        dss, I0, I1, I2, -24 * (I1 + I0**2 / 2), I2 + 2 * I0 * I1 + 2 * I0**3 / 3, 24 * (I0**2 / 2 - I1), I2 - I0**3 / 3))
# F2 型母関数から線形フリンジの z (4.33)
x1, p2, pp = sp.symbols('x_1 p_2 p', real=True)
a0, b0 = sp.symbols('a_0 b_0', real=True)
aa = a0 / pp
bq = b0 / pp                              # (4.33) の b (q の単位)
bc = b0 / pp**2                            # 正準運動量での係数 b'
F2 = sp.exp(aa) * x1 * p2 + sp.Rational(1, 2) * bc * sp.exp(aa) * p2**2
x2 = sp.diff(F2, p2)
p1 = sp.diff(F2, x1)
q2 = p2 / pp
q1 = sp.exp(aa) * q2
print('  x2 = ∂F/∂p2 - (e^a x1 + b q1) = %s、p1 - e^a p2 = %s' % (sp.simplify(x2 - (sp.exp(aa) * x1 + bq * q1)), sp.simplify(p1 - sp.exp(aa) * p2)))
print('  z2 - z1 = ∂F/∂p - [-(a x1 + b(1 + a/2) q2) q1] = %s' % sp.simplify(sp.diff(F2, pp) + (aa * x1 + bq * (1 + aa / 2) * q2) * q1))

# ------------------------------------------------------------------ (F)
print('\n== (F) 転送行列の導出 ==')
ell, pxs, pys, p = sp.symbols('ell p_x p_y p', positive=True)
pzz = sp.sqrt(p**2 - pxs**2 - pys**2)
r = ell * (1 / pzz - 1 / p)
xxd = pxs * r
print('  ∂x/∂px - (r + ℓ px^2/pz^3) = %s、∂x/∂py - ℓ px py/pz^3 = %s、∂r/∂p - ℓ(-p/pz^3 + 1/p^2) = %s' % (
    sp.simplify(sp.diff(xxd, pxs) - (r + ell * pxs**2 / pzz**3)), sp.simplify(sp.diff(xxd, pys) - ell * pxs * pys / pzz**3),
    sp.simplify(sp.diff(r, p) - ell * (-p / pzz**3 + 1 / p**2))))
Msym = sp.Matrix(6, 6, lambda i, j: sp.Symbol('M%d%d' % (i + 1, j + 1)))
J6 = sp.zeros(6, 6)
for i in range(3):
    J6[2 * i, 2 * i + 1] = 1
    J6[2 * i + 1, 2 * i] = -1
cond = (Msym.T * J6 * Msym - J6)
subs6 = {sp.Symbol('M66'): 1, **{sp.Symbol('M6%d' % k): 0 for k in range(1, 6)}}
Msy = lambda a__, b__: sp.Symbol('M%d%d' % (a__, b__))
print('  (M^T J M - J)_{1,6} (j = 1 の例、M6k = 0、M66 = 1) = %s' % sp.expand(cond[0, 5].subs(subs6)))
for j in range(1, 5):
    e = cond[j - 1, 5].subs(subs6)
    sol5 = sp.solve(e, Msy(5, j))[0]
    mine = -Msy(1, j) * Msy(2, 6) + Msy(2, j) * Msy(1, 6) - Msy(3, j) * Msy(4, 6) + Msy(4, j) * Msy(3, 6)
    print('  j = %d: M5%d = %s、(4.38) との差 %s' % (j, j, sol5, sp.simplify(sol5 - mine)))
# 連鎖律 (4.37)
pxx, pyy, dd_ = sp.symbols('P_x P_y d')
Fexpr = sp.Function('F')(pxx / (1 + dd_), pyy / (1 + dd_), dd_)
print('  連鎖律: ∂F/∂P_x = %s' % sp.diff(Fexpr, pxx))
print('  連鎖律: ∂F/∂δ|_Px = %s' % sp.diff(Fexpr, dd_))
