# -*- coding: utf-8 -*-
"""空洞の章 (10_cavity.md) の補完: シンプレクティック性の議論に出る式の確認。

(a) 薄い空洞のエネルギー利得 (t, E) -> (t, E + f(t)) のヤコビ行列と、(t, E) -> (z, δ) の座標変換のヤコビ行列 (sympy、§4)。
(b) 横のキックの符号 s = +1 (正しい) と s = -1 (逆) について、正準座標 (x, p_x, t, P = -E) のヤコビ行列から
    D = J^T Ω J - Ω を解析的に取り出す (sympy、§5.1)。
(c) (z, δ) で測った同じ量を、薄い空洞の写像 (10.1)(10.2)(10.5)(10.6)(10.8) の数値ヤコビアンと、
    (b) を座標変換した予想 D_SAD = (1/p0) A^{-T} D A^{-1} で比べる (numpy、§5.1)。
(d) tcavfrin の q の再スケール欠落による欠陥: J^T S J - S = g e_x^T - e_x g^T (g = ∇η、η = -(p_x + X) dpr) の確認 (numpy、§7.2)。
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

MC2 = 0.51099895e6
c_light = 299792458.0

# ---------------------------------------------------------------------------------------------------- (a)
print('== (a) (t, E) -> (t, E + f(t)) と、(t, E) -> (z, δ) の座標変換 (sympy) ==')
t, E, p0 = sp.symbols('t E p0', positive=True)
f = sp.Function('f')
pE = sp.sqrt(E**2 - 1)                      # p(E) (mc 単位)
bE = pE / E                                 # β(E)
G = sp.Matrix([-t * bE, pE / p0 - 1])       # (z, δ) = (-t β(E), p(E)/p0 - 1)
M = G.jacobian([t, E])                      # ∂(z, δ)/∂(t, E)
print('  ∂(z, δ)/∂(t, E) = %s' % sp.simplify(M))
print('  det [∂(z, δ)/∂(t, E)] = %s    (定数 -1/p0 になるはず)' % sp.simplify(M.det()))
Jte = sp.Matrix([[1, 0], [sp.diff(f(t), t), 1]])
print('  (t, E) -> (t, E + f(t)) のヤコビ行列 = %s、det = %s' % (Jte.tolist(), Jte.det()))
# (z, δ) 空間の写像のヤコビ行列 J_zδ = M(F(q)) J_tE M(q)^{-1}
MF = M.subs(E, E + f(t))
Jzd = MF * Jte * M.inv()
dJ = sp.simplify(Jzd.det())
print('  J_zδ = M(F) J_tE M^{-1} の det = %s' % dJ)
Sm = sp.Matrix([[0, 1], [-1, 0]])
print('  2 x 2 では J^T S J = det(J) S。したがって J_zδ^T S J_zδ - S = (det - 1) S = 0')
# 数値で成分を確認 (f(t) = 3e-4 sin(7 t)、t = 0.03、E = 5870 の電子 3 GeV 付近ではなく、E = 10)
vals = {p0: 4.0, t: 0.03, E: 10.0}
fv = lambda tt: 3e-4 * sp.sin(7 * tt)
Jn = sp.N(Jzd.subs(sp.Function('f')(t), fv(t)).doit().subs(vals), 15)
print('  数値 (f = 3e-4 sin 7t、t = 0.03、E = 10、p0 = 4): J_zδ^T S J_zδ - S の最大成分 = %.2e' % float(max(abs(x) for x in (Jn.T * Sm * Jn - Sm))))

# ---------------------------------------------------------------------------------------------------- (b)
print('\n== (b) 横のキックの符号 s (= ±1) と欠陥 D = J^T Ω J - Ω (sympy。正準座標 (x, p_x, t, P = -E)) ==')
x, tt, s = sp.symbols('x t s', real=True)
K = sp.Function('K')(x, tt)
# x' = x、p_x' = p_x - s ∂_x K、t' = t、P' = P - ∂_t K  (P = -E なので E' = E + ∂_t K)
u = sp.symbols('x px t P')
Kx, Kt = sp.diff(K, x), sp.diff(K, tt)
Fm = sp.Matrix([u[0], u[1] - s * Kx.subs({x: u[0], tt: u[2]}), u[2], u[3] - Kt.subs({x: u[0], tt: u[2]})])
Jc = Fm.jacobian(u)
Om = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
Dm = sp.simplify(Jc.T * Om * Jc - Om)
print('  J = %s' % Jc.applyfunc(lambda e: e).tolist())
print('  s = +1: D = %s' % Dm.subs(s, 1).tolist())
print('  s = -1: D = %s' % Dm.subs(s, -1).tolist())
print('  一般の s: D の 0 でない成分 = %s' % {(i, j): sp.factor(Dm[i, j]) for i in range(4) for j in range(4) if Dm[i, j] != 0})

# ---------------------------------------------------------------------------------------------------- (c)
print('\n== (c) 薄い空洞 (10.1)(10.2)(10.5)(10.6)(10.8) の数値ヤコビアン: s = +1 と s = -1 ==')
S6 = np.zeros((6, 6))
for i in (0, 2, 4):
    S6[i, i + 1] = 1
    S6[i + 1, i] = -1


def h_of(p):
    return math.sqrt(1.0 + p * p)


def thin(c, par, sgn):
    """c = (x, px, y, py, z, δ)。px、py は正準運動量を p0 で規格化したもの。offset = offset1 = 0 (RAD ON / TRPT の形)"""
    p0_, w_, dphis, v_, v1, v20, v02, v11 = par
    X, PX, Y, PY, Z, G_ = c
    p1 = p0_ * (1.0 + G_)
    h1 = h_of(p1)
    tt_ = -Z * h1 / p1                                          # (10.1)
    th = w_ * tt_ - dphis
    V = v_ + v1 * X + .5 * v20 * X * X + .5 * v02 * Y * Y + v11 * X * Y
    Vx = v1 + v20 * X + v11 * Y
    Vy = v02 * Y + v11 * X
    dh = -V * math.sin(th)                                      # (10.2)
    h2 = h1 + dh
    p2 = math.sqrt(h2 * h2 - 1.0)                               # (10.5)
    Z2 = -tt_ * p2 / h2                                         # (10.6)
    PX2 = PX + sgn * (-Vx * math.cos(th) / (w_ * p0_))          # (10.8)
    PY2 = PY + sgn * (-Vy * math.cos(th) / (w_ * p0_))
    return np.array([X, PX2, Y, PY2, Z2, p2 / p0_ - 1.0])


def jac(fn, c, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h * max(1.0, abs(c[j])) if j not in (0, 2) else 1e-6
        J[:, j] = (fn(c + d) - fn(c - d)) / (2 * d[j])
    return J


f_rf = 500e6
w_ = 2 * math.pi * f_rf / c_light
for lab, pe in (('3 GeV/c', 3.0e9), ('2 MeV/c', 2.0e6)):
    p0_ = pe / MC2
    v_ = 2.0e6 / MC2
    # V1 = 2e6 V/m、V20 = 2e8 V/m²、V02 = -1e8 V/m²、V11 = 1e8 V/m² を mc² で割った値 (v 単位、長さは m)
    par = (p0_, w_, -0.3, v_, 2e6 / MC2, 2e8 / MC2, -1e8 / MC2, 1e8 / MC2)
    c0 = np.array([2e-3, 1e-4, 1e-3, -5e-5, 0.05, 1e-3])
    for sg in (+1, -1):
        J = jac(lambda cc: thin(cc, par, sg), c0)
        Dm_ = J.T @ S6 @ J - S6
        print('  %-8s s = %+d: |J^T S J - S| = %.2e (最大成分)' % (lab, sg, np.abs(Dm_).max()))
    # 解析の予想: D_SAD = (1/p0) A^{-T} D_can A^{-1}  (s = -1)。A = ∂W/∂u、u = (x, p_x, y, p_y, t, P = -E)、W = (x, p_x/p0, y, p_y/p0, z, δ)
    X, PX, Y, PY, Z, G_ = c0
    p1 = p0_ * (1 + G_)
    h1 = h_of(p1)
    t0 = -Z * h1 / p1
    th = w_ * t0 - par[2]
    V = v_ + par[4] * X + .5 * par[5] * X * X + .5 * par[6] * Y * Y + par[7] * X * Y
    Vx = par[4] + par[5] * X + par[7] * Y
    Vy = par[6] * Y + par[7] * X
    wi = 1.0 / w_
    # K = V(x, y) cos(θ)/ω。K_ab (a、b = x、y、t)
    Vxx, Vyy, Vxy = par[5], par[6], par[7]
    cs, sn = math.cos(th), math.sin(th)
    Kxx, Kyy, Kxy = Vxx * cs * wi, Vyy * cs * wi, Vxy * cs * wi
    Kxt, Kyt = -Vx * sn, -Vy * sn
    Ktt = -w_ * V * cs
    # u = (x, px, y, py, t, P): x' = x、px' = px - s K_x、y' = y、py' = py - s K_y、t' = t、P' = P - K_t
    def Jcan(sg):
        Jm = np.eye(6)
        Jm[1, 0] = -sg * Kxx
        Jm[1, 2] = -sg * Kxy
        Jm[1, 4] = -sg * Kxt
        Jm[3, 0] = -sg * Kxy
        Jm[3, 2] = -sg * Kyy
        Jm[3, 4] = -sg * Kyt
        Jm[5, 0] = -Kxt
        Jm[5, 2] = -Kyt
        Jm[5, 4] = -Ktt
        return Jm
    Om6 = np.zeros((6, 6))
    for i in (0, 2, 4):
        Om6[i, i + 1] = 1
        Om6[i + 1, i] = -1
    # A = ∂W/∂u: x、y はそのまま、p_x → px_s = p_x/p0、(t, P) → (z, δ)。z = -t β(E)、δ = p(E)/p0 - 1、E = h1、P = -E
    E_ = h1
    bt = p1 / h1
    dz_dt = -bt
    dz_dE = -t0 * (1.0 / (h1 * h1 * h1)) / bt * 1.0      # dβ/dE = 1/(E³ β) (β = p/E、p = √(E² - 1))
    dd_dE = 1.0 / (bt * p0_)
    A = np.zeros((6, 6))
    A[0, 0] = 1.0
    A[1, 1] = 1.0 / p0_
    A[2, 2] = 1.0
    A[3, 3] = 1.0 / p0_
    A[4, 4] = dz_dt
    A[4, 5] = -dz_dE                                     # P = -E なので ∂/∂P = -∂/∂E
    A[5, 5] = -dd_dE
    Ainv = np.linalg.inv(A)
    for sg in (+1, -1):
        Dcan = Jcan(sg).T @ Om6 @ Jcan(sg) - Om6
        Dpred = (1.0 / p0_) * Ainv.T @ Dcan @ Ainv
        Jn = jac(lambda cc: thin(cc, par, sg), c0)
        Dnum = Jn.T @ S6 @ Jn - S6
        print('  %-8s s = %+d: 予想 (1/p0) A^-T D A^-1 の最大成分 %.3e、数値の最大成分 %.3e、差の最大 %.1e' % (
            lab, sg, np.abs(Dpred).max(), np.abs(Dnum).max(), np.abs(Dpred - Dnum).max()))
    print('  %-8s ∂_x∂_t K = -V_x sin(ωt - dphis) = %.3e、2 |K_xt| = %.3e (s = -1 の D_can の x-t 成分の大きさ)' % (lab, Kxt, 2 * abs(Kxt)))


# ---------------------------------------------------------------------------------------------------- (d)
print('\n== (d) tcavfrin (入口の端だけ) の q の再スケール欠落: J^T S J - S = g e_x^T - e_x g^T + (y 方向の同じ項)、g = ∇η ==')


def fringe(c, par, rescale):
    p0_, h0_, w_, dphis, offset, al, v_ = par
    X, PX, Y, PY, Z, G_ = c
    p1r = 1.0 + G_
    p1 = p0_ * p1r
    h1 = h_of(p1)
    qx, qy = PX / p1r, PY / p1r
    vf_ = v_ / al / p0_ * .5
    s0_ = math.sin(dphis) - offset
    tt_ = -Z * h1 / p1
    ph = w_ * tt_
    dpt = vf_ * (math.sin(ph - dphis) + s0_) / p1r
    qx2 = qx + X * dpt
    qy2 = qy + Y * dpt
    dh = -w_ * p0_ * vf_ * (X * X + Y * Y) * .5 * math.cos(ph - dphis)
    h2 = h1 + dh
    a = max(dh * (h1 + h2) / p1**2, -1.0)
    dpr = a / (1.0 + math.sqrt(1.0 + a))
    dp2r = G_ + p1r * dpr
    p2r = 1.0 + dp2r
    Z2 = Z * p2r / p1r * h1 / h2
    if rescale:
        qx2 = qx * p1r / p2r + X * dpt * p1r / p2r
        qy2 = qy * p1r / p2r + Y * dpt * p1r / p2r
    return np.array([X, qx2 * p2r, Y, qy2 * p2r, Z2, dp2r])


def eps_fn(c, par, idx=1):
    return fringe(c, par, True)[idx] - fringe(c, par, False)[idx]


V, al, phi = 1.0e5, 0.1, 0.5
vv = V / MC2
print('  CAVI 入口の端だけ (L = %.2f m、VOLT = %.0e V、500 MHz、PHI = %.1f、dphis = -PHI)。z = 0.02 m、δ = 1e-3' % (al, V, phi))
print('  η ≡ (q の再スケールあり) - (なし) の p_x 出力 = -(p_x + X_k) dpr (X_k = x v_f [sin + s0])。確認: η の式との差、')
print('  D の成分 D[px,x] = ∂η/∂p_x = -dpr、D[z,x]、D[δ,x] と予想 g e_x^T - e_x g^T との差')
for lab, pe in (('3 GeV/c', 3.0e9), ('2 MeV/c', 2.0e6), ('0.6 MeV/c', 0.6e6)):
    p0_ = pe / MC2
    par = (p0_, h_of(p0_), w_, -phi, math.sin(-phi), al, vv)
    for pxv in (0.01, 0.0):
        for xv in (1e-3, 1e-2):
            c0 = np.array([xv, pxv, 0.0, 0.0, 0.02, 1e-3])
            Jk = jac(lambda cc: fringe(cc, par, False), c0)
            Dk = Jk.T @ S6 @ Jk - S6
            g = np.zeros(6)
            gy = np.zeros(6)
            for j in range(6):
                d = np.zeros(6)
                d[j] = 1e-7 * max(1.0, abs(c0[j])) if j not in (0, 2) else 1e-6
                g[j] = (eps_fn(c0 + d, par, 1) - eps_fn(c0 - d, par, 1)) / (2 * d[j])
                gy[j] = (eps_fn(c0 + d, par, 3) - eps_fn(c0 - d, par, 3)) / (2 * d[j])
            ex = np.zeros(6)
            ex[0] = 1.0
            ey = np.zeros(6)
            ey[2] = 1.0
            # p_y の出力にも同じ欠落 η_y = -(p_y + Y_k) dpr がある。D = g e_x^T - e_x g^T + g_y e_y^T - e_y g_y^T
            Dpred = np.outer(g, ex) - np.outer(ex, g) + np.outer(gy, ey) - np.outer(ey, gy)
            # η の式
            p1r = 1 + c0[5]
            p1 = p0_ * p1r
            h1 = h_of(p1)
            vf_ = vv / al / p0_ * .5
            th = w_ * (-c0[4] * h1 / p1) + phi
            dh = -w_ * p0_ * vf_ * xv**2 * .5 * math.cos(th)
            a = dh * (h1 + h1 + dh) / p1**2
            dpr = a / (1 + math.sqrt(1 + a))
            Xk = xv * vf_ * math.sin(th)                      # s0 = sin(dphis) - offset = 0
            eps_formula = -(c0[1] + Xk) * dpr
            print('   %-10s p_x = %.2f x = %4.0f mm: |D| = %.2e、|D - 予想| = %.1e、D[px,x] = %.3e (-dpr = %.3e)、g_z = %.2e、g_δ = %.2e、η = %.3e (式 %.3e)' % (
                lab, pxv, xv * 1e3, np.abs(Dk).max(), np.abs(Dk - Dpred).max(), Dk[1, 0], -dpr, g[4], g[5], eps_fn(c0, par), eps_formula))
print('  (中心差分の精度は 1e-10 程度。D[px,x] が -dpr に合い、p_x = 0 でも同じ大きさで残る。)')
