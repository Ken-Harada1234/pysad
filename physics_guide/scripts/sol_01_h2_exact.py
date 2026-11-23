# -*- coding: utf-8 -*-
"""ソレノイド + 四極の H2 の厳密解 (09_solenoid.md §3〜§5)。

(a) sympy: 運動方程式、特性方程式 λ⁴ + (b/p)² λ² - (K/p)² = 0 の根 (w1, w2) の関係式。
(b) 数値: 行列指数 exp(A ℓ) (4×4、正準変数) と、ヘルプの式 (w1, w2, w+, U, a, b, c, d)、K64 の tsolqu の式 (pysad の _TZ) の一致。
(c) z の更新: 厳密な ODE の -∫ (q_x² + q_y²)/2 ds と、K64 の tztafs を使った dz の一致。
(d) ΔH の解 (ヘルプの solution-dH): 回転角 w ℓ、z の式 (3/2 - p/pz - pz²/(2p²))ℓ が、tsolqu のコード (r、(3 + dpz) ap/(2(2 + dpz)) r) と一致 (sympy)。数値の流れ (ODE) との比較。
"""
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

print('== (a) 運動方程式と特性方程式 (sympy) ==')
b_, p_, K_, lam = sp.symbols('b p K lambda', positive=True)
bp, kp = b_ / p_, K_ / p_
# x'' = bp y' - kp x、y'' = -bp x' + kp y。 x = X e^{λ s}、y = Y e^{λ s}
M = sp.Matrix([[lam**2 + kp, -bp * lam], [bp * lam, lam**2 - kp]])
poly = sp.expand(M.det())
print('  det[[λ²+K/p, -(b/p)λ], [(b/p)λ, λ²-K/p]] = %s' % poly)
Vv = sp.sqrt(bp**4 + 4 * kp**2)
w1 = sp.sqrt((bp**2 + Vv) / 2)
w2 = kp / w1
print('  w1 = sqrt(((b/p)² + V)/2)、w2 = (K/p)/w1、V = sqrt((b/p)⁴ + 4 (K/p)²)')
print('  根: λ² = -w1² → det = %s、λ² = +w2² → det = %s' % (sp.simplify(poly.subs(lam, sp.I * w1)), sp.simplify(poly.subs(lam, w2))))
print('  w1² + w2² - V = %s、w1² - w2² - (b/p)² = %s、w1 w2 - K/p = %s' % (
    sp.simplify(w1**2 + w2**2 - Vv), sp.simplify(w1**2 - w2**2 - bp**2), sp.simplify(w1 * w2 - kp)))

print('\n== (b) 行列指数・ヘルプの式・K64 (_TZ) の比較 ==')
Bz, K, ell, dl = 1.3, 0.9, 0.4, 0.01      # b_z = Bz/Bρ [1/m]、K = K1/L [1/m²]、長さ ℓ [m]、δ
p = 1 + dl
x0, px0, y0, py0 = 0.003, 0.002, -0.002, 0.0015     # 正準 (p0 で規格化)


def amat(Bz, K, p):
    """d/ds (x, px, y, py)。P = p_c - a、a = (-Bz y/2, Bz x/2)"""
    A = np.zeros((4, 4))
    A[0, 1] = 1 / p
    A[0, 2] = Bz / (2 * p)
    A[2, 3] = 1 / p
    A[2, 0] = -Bz / (2 * p)
    # px' = (Bz/(2p)) (py - Bz x/2) - K x、py' = -(Bz/(2p)) (px + Bz y/2) + K y
    A[1, 3] = Bz / (2 * p)
    A[1, 0] = -Bz**2 / (4 * p) - K
    A[3, 1] = -Bz / (2 * p)
    A[3, 2] = -Bz**2 / (4 * p) + K
    return A


s0 = np.array([x0, px0, y0, py0])
ref = expm(amat(Bz, K, p) * ell) @ s0


def help_form(Bz, K, p, ell, s0):
    x0, px0, y0, py0 = s0
    bp = Bz / p
    V = math.sqrt(bp**4 + 4 * (K / p)**2)
    w1 = math.sqrt((bp**2 + V) / 2)
    w2 = K / (p * w1)
    wp = w1 + w2
    U = 1 / V
    pxm, pym = px0 + Bz / 2 * y0, py0 - Bz / 2 * x0
    a = U * (w2 * wp * x0 - Bz / p**2 * pym)
    b = w1 * U / p * (wp * pxm - Bz * w2 * y0)
    c = U / p * (w1 * Bz / wp * x0 + pym)
    d = U * w2 * (-Bz / (p**2 * wp) * pxm + w1 * y0)
    f1, f2 = w1 * ell, w2 * ell
    u1 = a * w1 * (math.cos(f1) - 1) + b * math.sin(f1)
    u2 = -a * w1 * math.sin(f1) + b * (math.cos(f1) - 1)
    v1 = c * w2 * (math.cosh(f2) - 1) + d * math.sinh(f2)
    v2 = c * w2 * math.sinh(f2) + d * (math.cosh(f2) - 1)
    dx = u1 / w1 + v1 * Bz / (p * w2)
    dy = u2 * Bz / (p * w1 * wp) + wp * v2 / w2
    return np.array([x0 + dx, px0 + p * u2 + Bz * (v2 - dy / 2), y0 + dy, py0 + p * wp * v1 + Bz * (-u1 / wp + dx / 2)])


hf = help_form(Bz, K, p, ell, s0)


def k64_step(Bz, K, p, ell, s0):
    """tsolqu.f の H2 の部分 (行 117〜142、運動学的なドリフト r = 0)。px, py は q = p_c/p"""
    tz = t6._TZ(p - 1, ell, K, Bz)
    x, qx, y, qy = s0[0], s0[1] / p, s0[2], s0[3] / p
    bzp = tz.bzp
    pxi, pyi = qx + bzp * y * .5, qy - bzp * x * .5            # ラーモア系 (運動学的運動量/p)
    a = (tz.w2 * tz.ws * x - bzp * pyi) * tz.wss
    bw = (tz.ws * pxi - bzp * tz.w2 * y) * tz.wss
    c = (tz.w1 * tz.wd * x + pyi) * tz.wss
    dw = (-tz.wd * pxi + tz.w1 * y) * tz.wss
    u1w = a * tz.dc1 + bw * tz.s1
    u2w = -a * tz.s1 + bw * tz.dc1
    v1w = c * tz.dch2 + dw * tz.sh2
    v2w = c * tz.sh2 + dw * tz.dch2
    xn = x + u1w + v1w * bzp
    pxn = pxi + u2w * tz.w1 + v2w * tz.w2 * bzp
    yn = y + tz.wd * u2w + tz.ws * v2w
    pyn = pyi - tz.wd * u1w * tz.w1 + tz.ws * v1w * tz.w2
    awu = a / tz.ws * tz.w1
    dwu = dw * tz.w2
    dz1 = tz.tztafs(awu, pxi, pyi, tz.aw1, tz.ws, tz.w12, tz.wss, tz.g1)
    dz2 = tz.tztafs(-dwu, -pyi, pxi, tz.aw2, -tz.w12, tz.ws, -tz.wss, tz.g2)
    dz = bzp * (-((awu * dwu * tz.dxs**2) / tz.akkp) + tz.ca1 * pxi * pyi * tz.wss) + dz1 + dz2
    # q → 正準: px = (pxn - bzp y/2) p
    return np.array([xn, (pxn - bzp * yn * .5) * p, yn, (pyn + bzp * xn * .5) * p]), dz, tz


k64, dzk, tz = k64_step(Bz, K, p, ell, s0)
print('  Bz = %.1f /m、K = %.1f /m²、ℓ = %.1f m、δ = %.2f (w1 = %.6f、w2 = %.6f、w+ = %.6f /m)' % (Bz, K, ell, dl, tz.w1, tz.w2, tz.ws))
print('  行列指数 = %s' % np.array2string(ref, precision=12))
print('  |ヘルプの式 - 行列指数| = %.2e   |K64 (tsolqu の式) - 行列指数| = %.2e' % (np.abs(hf - ref).max(), np.abs(k64 - ref).max()))
print('  対応 (コードの変数 ↔ ヘルプ): a = a、bw = b/w1、c = c、dw = d/w2、u1w = u1/w1、u2w = u2/w1、v1w = v1/w2、v2w = v2/w2、wss = U、wd = (Bz/p)/w+')
for (Bz_, K_, dl_) in ((2.5, 1.2, -0.02), (0.3, 4.0, 0.02)):
    p2 = 1 + dl_
    r2 = expm(amat(Bz_, K_, p2) * ell) @ s0
    h2 = help_form(Bz_, K_, p2, ell, s0)
    k2, _, _ = k64_step(Bz_, K_, p2, ell, s0)
    print('  (Bz, K, δ) = (%.1f, %.1f, %+.2f): |ヘルプ - 行列指数| = %.1e、|K64 - 行列指数| = %.1e' % (Bz_, K_, dl_, np.abs(h2 - r2).max(), np.abs(k2 - r2).max()))

print('\n== (c) z の更新 ==')


def zode(Bz, K, p, ell, s0):
    A = amat(Bz, K, p)

    def rhs(s, v):
        st = v[:4]
        dd = A @ st
        x, px, y, py = st
        Px, Py = px + Bz * y / 2, py - Bz * x / 2
        return list(dd) + [-(Px**2 + Py**2) / (2 * p**2)]
    return solve_ivp(rhs, (0, ell), list(s0) + [0.0], rtol=1e-13, atol=1e-18, method='DOP853').y[4, -1]


zo = zode(Bz, K, p, ell, s0)
print('  Δz (ODE: -∫ (q_x² + q_y²)/2 ds、q = P/p) = %.10e' % zo)
print('  Δz (K64: dz = bzp(…) + dz1 + dz2)       = %.10e   (差 %.1e)' % (dzk, dzk - zo))
for (Bz_, K_, dl_) in ((2.5, 1.2, -0.02), (0.3, 4.0, 0.02)):
    p2 = 1 + dl_
    _, dzk2, _ = k64_step(Bz_, K_, p2, ell, s0)
    zo2 = zode(Bz_, K_, p2, ell, s0)
    print('  (Bz, K, δ) = (%.1f, %.1f, %+.2f): ODE %.8e、K64 %.8e、差 %.1e' % (Bz_, K_, dl_, zo2, dzk2, dzk2 - zo2))

print('\n== (d) ΔH の解 ==')
d = sp.symbols('d', real=True)                 # dpz = pz/p - 1 (tsolqu の dpz)
ell_s = sp.symbols('ell', positive=True)
pzr = 1 + d                                    # pz/p
z_help = (sp.Rational(3, 2) - 1 / pzr - pzr**2 / 2) * ell_s
ap = 1 - pzr**2                                # q_x² + q_y² = 1 - (pz/p)²
r = -d / (1 + d) * ell_s
z_code = -(3 + d) * ap / 2 / (2 + d) * r
print('  z: ヘルプ (3/2 - p/pz - pz²/(2p²))ℓ と tsolqu の -(3 + dpz) ap r/(2(2 + dpz)) の差 = %s' % sp.simplify(z_help - z_code))
bsym, pr_, pz_ = sp.symbols('b p p_z', positive=True)
w_help = bsym * (pr_ - pz_) / (pr_ * pz_)
w_code = bsym / pr_ * (-d / (1 + d))
print('  回転角: ヘルプ w = b(p - pz)/(p pz) と、コード r·bzp ((-dpz/(1+dpz)) b/p) の差 = %s (pz = p(1 + dpz) とおく)' % sp.simplify(w_help.subs(pz_, pr_ * (1 + d)) - w_code))

# 数値: ΔH = p - pz - P²/(2p) (P = p_c - a) の流れ
Bz2, ell2 = 2.0, 0.5
s1_ = np.array([0.004, 0.02, -0.003, 0.015])
pp = 1.02


def dh_rhs(s, v):
    """dx/ds = ∂H/∂px、dp/ds = -∂H/∂x、dz/ds = ∂H/∂p (固定 p_c、x)。
    H = p - pz - (Px² + Py²)/(2p)。∂H/∂Px = Px/pz - Px/p。∂H/∂p = 1 - p/pz + P²/(2p²)"""
    x, px, y, py, z = v
    Px, Py = px + Bz2 * y / 2, py - Bz2 * x / 2
    pz = math.sqrt(pp * pp - Px * Px - Py * Py)
    fx = Px / pz - Px / pp
    fy = Py / pz - Py / pp
    return [fx, fy * Bz2 / 2, fy, -fx * Bz2 / 2, 1 - pp / pz + (Px**2 + Py**2) / (2 * pp**2)]


v = solve_ivp(dh_rhs, (0, ell2), list(s1_) + [0.0], rtol=1e-13, atol=1e-18, method='DOP853').y[:, -1]
Px0, Py0 = s1_[1] + Bz2 * s1_[2] / 2, s1_[3] - Bz2 * s1_[0] / 2
pz0 = math.sqrt(pp * pp - Px0**2 - Py0**2)
w_ = Bz2 * (pp - pz0) / (pp * pz0)
sh = {}
for name, fac in (('ヘルプの式のまま (1/B_z なし)', 1.0), ('1/B_z をつけた式', 1.0 / Bz2)):
    dxh = fac * (Px0 * math.sin(w_ * ell2) - Py0 * (math.cos(w_ * ell2) - 1))
    dyh = fac * (Px0 * (math.cos(w_ * ell2) - 1) + Py0 * math.sin(w_ * ell2))
    xh, yh = s1_[0] + dxh, s1_[2] + dyh
    pxh, pyh = s1_[1] + Bz2 / 2 * dyh, s1_[3] - Bz2 / 2 * dxh
    zh = (1.5 - pp / pz0 - pz0**2 / (2 * pp**2)) * ell2
    sh[name] = np.array([xh, pxh, yh, pyh, zh])
print('  Bz = %.1f、ℓ = %.1f、p = %.2f: ODE (ΔH の流れ) = %s' % (Bz2, ell2, pp, np.array2string(v, precision=12)))
for name, arr in sh.items():
    print('  %s: 差の最大 = %.2e (x, px, y, py, z)' % (name, np.abs(v - arr).max()))
print('  → Δx、Δy は 1/B_z が要る (速度 = P (1/pz - 1/p) = P w/B_z。回転角 w = B_z (p - pz)/(p pz))。tsolqu は (a24 pxi - a14 py)/bzp と 1/bzp を持つ。z は一致。')

print('\n== (e) 固有モードと a, bw, c, dw の意味 ==')
# 固有ベクトル: (λ² + K/p) X = (b/p) λ Y から Y/X = (λ² + kp)/(bp λ)
bp_s, kp_s = sp.symbols('b_p k_p', positive=True)
w1s, w2s = sp.symbols('w1 w2', positive=True)
sub = {kp_s: w1s * w2s}
bp2 = w1s**2 - w2s**2
for name, lamv in (('+i w1', sp.I * w1s), ('-i w1', -sp.I * w1s), ('+w2', w2s), ('-w2', -w2s)):
    ratio = sp.simplify(((lamv**2 + w1s * w2s) / (sp.sqrt(bp2) * lamv)))
    wd_s = sp.sqrt(bp2) / (w1s + w2s)
    print('  λ = %s: Y/X = %s  (wd = (b/p)/w+ = %s。±i wd または ±1/wd と比べる: %s、%s)' % (
        name, ratio, sp.simplify(wd_s), sp.simplify(ratio / wd_s), sp.simplify(ratio * wd_s)))
tz = t6._TZ(dl, ell, K, Bz)
bzp = tz.bzp
x, qx, y, qy = x0, px0 / p, y0, py0 / p
pxi, pyi = qx + bzp * y * .5, qy - bzp * x * .5
# 4 つの初期条件で A, B, C, D を解く: x = A cos + B sin + (C cosh + D sinh)、y = wd(-A sin + B cos) + (C sinh + D cosh)/wd
Mx = np.array([[1, 0, 1, 0], [0, tz.wd, 0, 1 / tz.wd], [0, tz.w1, 0, tz.w2], [-tz.wd * tz.w1, 0, tz.w2 / tz.wd, 0]])
ABCD = np.linalg.solve(Mx, np.array([x, y, pxi, pyi]))
a = (tz.w2 * tz.ws * x - bzp * pyi) * tz.wss
bw = (tz.ws * pxi - bzp * tz.w2 * y) * tz.wss
c = (tz.w1 * tz.wd * x + pyi) * tz.wss
dw = (-tz.wd * pxi + tz.w1 * y) * tz.wss
print('  初期条件から解いた (A, B, C, D) = %s' % np.array2string(ABCD, precision=10))
print('  K64 の (a, bw, bzp c, bzp dw)    = %s   (差 %.1e)' % (np.array2string(np.array([a, bw, bzp * c, bzp * dw]), precision=10), np.abs(ABCD - np.array([a, bw, bzp * c, bzp * dw])).max()))
