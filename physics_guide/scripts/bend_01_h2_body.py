# -*- coding: utf-8 -*-
"""偏向磁石の本体の H2 と、その厳密解 (05_bend.md §3〜§5、式 (5.1)〜(5.20))。

(a) A_s から磁場を出す (一様な B0)、a_s、ANGLE ≠ K0 (ρ0 ≠ ρb) の H の 2 次の展開 (sympy)。
(b) 偏微分を 1 項ずつ、運動方程式。Fortran の akxsq・drhop・akysq との一致。
(c) 水平の厳密解が運動方程式を満たすこと、分散 (R16・R26)・R56 の項 (sympy)。
(d) K64 の tbendibody の式 (x, q, y, z の更新) と、H2 の運動方程式の DOP853 の積分との比較。K64 の xsxkx の符号 (akxsq < 0 の枝)。
(e) pysad の track6._tbendiinit・_tbendibody (k64・v2) との比較。
(f) H2 の運動方程式と、元の H (根号つき) の運動方程式の差 (展開の次数): 振幅を 2 倍にすると差が何倍になるか。
"""
import math
import types
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)

x, y, px, py = sp.symbols('x y p_x p_y', real=True)
p = sp.Symbol('p', positive=True)
h, hb, K, s, d = sp.symbols('h h_b K s delta', real=True)    # h = 1/ρ0、hb = 1/ρb (場の曲率)

print('== (a) A_s と磁場、a_s、H の展開 (ρb ≠ ρ0 を許す) ==')
B0 = sp.Symbol('B0', real=True)
As = -(x + h * x**2 / 2) * B0
Bx = sp.diff(As, y) / (1 + h * x)
By = -sp.diff(As, x) / (1 + h * x)
print('  A_s = -(x + h x^2/2) B0 から B_x = %s、B_y = %s  (曲線座標の回転の式 B_x = ∂_y A_s/(1+hx)、B_y = -∂_x A_s/(1+hx))' % (sp.simplify(Bx), sp.simplify(By)))
print('  a_s = e A_s/P0 = -(x + h x^2/2) hb  (hb = eB0/P0 = 1/ρb)')
eps = sp.Symbol('epsilon', positive=True)
Q = px**2 + py**2
Hfull = (x + h * x**2 / 2) * hb - (1 + h * x) * sp.sqrt(p**2 - Q)
# x, px, py を ε 倍して展開
Hs = Hfull.subs({x: eps * x, px: eps * px, py: eps * py}, simultaneous=True)
ser = sp.series(Hs, eps, 0, 4).removeO()
print('  H = (x + h x^2/2) hb - (1 + h x) sqrt(p^2 - px^2 - py^2) を (x, px, py) の次数で展開:')
cf = {}
for n in range(4):
    cf[n] = sp.expand(ser.coeff(eps, n))
    print('    %d 次: %s' % (n, sp.simplify(cf[n])))
H2 = cf[0] + cf[1] + cf[2]
H2_expect = -p + (hb - h * p) * x + h * hb * x**2 / 2 + Q / (2 * p)
print('  0〜2 次の和が -p + (hb - h p) x + h hb x^2/2 + (px^2 + py^2)/(2p) に等しい: 差 = %s' % sp.simplify(H2 - H2_expect))
print('  hb = h (ρb = ρ0) のとき: (hb - h p) x = -h δ x、h hb x^2/2 = x^2/(2ρ^2) で、01a の式 (1a.62) と同じ。')

print('\n== (b) 偏微分を 1 項ずつ、運動方程式 (定数項 -p は落とす) ==')
Kq = sp.Symbol('K', real=True)
terms = [('(hb - h p) x', (hb - h * p) * x), ('h hb x^2/2', h * hb * x**2 / 2), ('K (x^2 - y^2)/2 (複合機能の K1)', Kq * (x**2 - y**2) / 2),
         ('px^2/(2p)', px**2 / (2 * p)), ('py^2/(2p)', py**2 / (2 * p))]
print('  %-34s %-12s %-18s %-10s %s' % ('項', '∂/∂px', '∂/∂x', '∂/∂py', '∂/∂y'))
for name, t in terms:
    print('  %-34s %-12s %-18s %-10s %s' % (name, sp.diff(t, px), sp.diff(t, x), sp.diff(t, py), sp.diff(t, y)))
Ht = sum(t for _, t in terms)
print('  合計: x\' = %s、px\' = %s、y\' = %s、py\' = %s' % (sp.diff(Ht, px), -sp.diff(Ht, x), sp.diff(Ht, py), -sp.diff(Ht, y)))
# Fortran: q = px/p、rhoe = rhob p、rhosq = rho0 rhoe、drhop = (rhoe - rho0)/rhosq、akxsq = -1/rhosq - akk/p、akysq = akk/p
rho0, rhob = 1 / h, 1 / hb
rhoe = rhob * p
rhosq = rho0 * rhoe
drhop = (rhoe - rho0) / rhosq
akxsq = -1 / rhosq - Kq / p
akysq = Kq / p
qx = sp.Symbol('q', real=True)
qdot = sp.simplify((-sp.diff(Ht, x)).subs(px, p * qx) / p)    # q' = px'/p
print('  q = px/p の方程式: q\' = %s' % qdot)
print('  Fortran の drhop + akxsq x と q\' の差 (x, q の係数と定数): %s' % sp.simplify(qdot - (drhop + akxsq * x)))
print('  Fortran の akysq y と q_y\' の差: %s' % sp.simplify((sp.diff(Ht, y) * -1).subs(py, p * sp.Symbol('qy')) / p - akysq * y))

print('\n== (c) 水平の厳密解、分散、R56 ==')
om = sp.Symbol('omega', positive=True)       # ω^2 = -akxsq (K=0 では h hb/p)
xe, x0, q0 = sp.symbols('x_e x_0 q_0', real=True)
xs = xe + (x0 - xe) * sp.cos(om * s) + q0 / om * sp.sin(om * s)
print('  x(s) = x_e + (x0 - x_e) cos ωs + (q0/ω) sin ωs は x\'\' = -ω^2 (x - x_e) を満たす: 残差 %s、初期条件 x(0) - x0 = %s、x\'(0) - q0 = %s' % (
    sp.simplify(sp.diff(xs, s, 2) + om**2 * (xs - xe)), sp.simplify(xs.subs(s, 0) - x0), sp.simplify(sp.diff(xs, s).subs(s, 0) - q0)))
print('  平衡点 x_e = drhop/(-akxsq) = %s  (= ρe - ρ0。K = 0 のとき)' % sp.simplify((drhop / (-akxsq)).subs(Kq, 0)))
# 分散: hb = h (ρb = ρ0 = ρ)、p = 1 + δ、一次の δ
rr = sp.Symbol('rho', positive=True)
L = sp.Symbol('L', positive=True)
dl = sp.Symbol('delta', real=True)
pp = 1 + dl
rhoe_ = rr * pp
xeq = rhoe_ - rr
omg = 1 / sp.sqrt(rr * rhoe_)
xf = xeq + (0 - xeq) * sp.cos(omg * L)      # x0 = q0 = 0
qf = (0 - xeq) * (-omg) * sp.sin(omg * L)    # q = x'
R16 = sp.simplify(sp.diff(xf, dl).subs(dl, 0))
R26 = sp.simplify(sp.diff(qf * pp, dl).subs(dl, 0))     # 正準 px = p q
print('  x0 = q0 = 0、1 次の δ: R16 = ∂x/∂δ = %s (= ρ(1 - cos(L/ρ)))、R26 = ∂px/∂δ = %s (= sin(L/ρ))' % (R16, R26))
# z の更新: z' = -h x - (q^2)/2 - dv。x0 = q0 = 0 の軌道で、δ の 1 次
zint = sp.integrate(sp.series((-(xeq + (-xeq) * sp.cos(omg * s)) / rr - qf.subs(L, s)**2 / 2), dl, 0, 2).removeO(), (s, 0, L))
R56_h = sp.simplify(sp.diff(zint, dl).subs(dl, 0))
print('  z\' = -h x - q^2/2 - dv の積分の δ の 1 次: ∂z/∂δ = %s - dv の δ 微分 × L  (= -(L - ρ sin(L/ρ)) + L/γ0^2)' % R56_h)
print('  展開: R16 = %s、R56(幾何) = %s (L/ρ の 3 次まで)' % (sp.series(R16, L, 0, 5), sp.series(R56_h, L, 0, 6)))

print('\n== (d) K64 の tbendibody の式と、H2 の運動方程式の DOP853 の積分 ==')


def fortran_body(rho0, rhob, K1L, al, dp, c0, xsxkx_k64):
    """tbend.f tbendiinit・tbendibody の transcription。c0 = (x, q, y, qy)。戻り値: (x, q, y, qy, dz)。xsxkx_k64: True なら akxsq < 0 の枝の符号が K64 のまま"""
    p_ = 1 + dp
    rhoe_ = rhob * p_
    rhosq_ = rho0 * rhoe_
    drhop_ = (rhoe_ - rho0) / rhosq_
    akk = K1L / al
    akxsq_ = -1 / rhosq_ - akk / p_
    akysq_ = akk / p_
    if akxsq_ < 0:
        akx = math.sqrt(-akxsq_)
        ph = akx * al
        dcx = -2 * math.sin(.5 * ph)**2
        spx = math.sin(ph)
        xspx = ph - spx
        dcxkx = dcx / akxsq_
        sxkx = spx / akx
        xsxkx = (-xspx / akx / akxsq_) if xsxkx_k64 else (xspx / akx / akxsq_)
    elif akxsq_ > 0:
        akx = math.sqrt(akxsq_)
        ph = akx * al
        dcx = 2 * math.sinh(.5 * ph)**2
        spx = math.sinh(ph)
        xspx = ph - spx
        dcxkx = dcx / akxsq_
        sxkx = spx / akx
        xsxkx = xspx / akx / akxsq_
    else:
        dcx, dcxkx, sxkx, xsxkx = 0.0, .5 * al**2, al, -al**3 / 6
    if akysq_ > 0:
        aky = math.sqrt(akysq_)
        dcy = 2 * math.sinh(.5 * aky * al)**2
        syky = math.sinh(aky * al) / aky
        aksy = aky * math.sinh(aky * al)
    elif akysq_ < 0:
        aky = math.sqrt(-akysq_)
        dcy = -2 * math.sin(.5 * aky * al)**2
        syky = math.sin(aky * al) / aky
        aksy = -aky * math.sin(aky * al)
    else:
        dcy, syky, aksy = 0.0, al, 0.0
    xi, pxi, yi, pyi = c0
    xiksq = drhop_ + akxsq_ * xi
    dxf = dcxkx * xiksq + sxkx * pxi
    dpxf = sxkx * xiksq + dcx * pxi
    dyf = dcy * yi + syky * pyi
    dpyf = aksy * yi + dcy * pyi
    hi = .5 * (pxi**2 + pyi**2 - akxsq_ * xi**2 - akysq_ * yi**2) - drhop_ * xi
    dz = -.25 * (dxf * pxi + xi * dpxf + dxf * dpxf + dyf * pyi + yi * dpyf + dyf * dpyf
                 - (5 / rho0 - 1 / rhoe_) * (drhop_ * xsxkx - sxkx * xi - dcxkx * pxi)) - hi * al * .5
    return xi + dxf, pxi + dpxf, yi + dyf, pyi + dpyf, dz


def ode_h2(rho0, rhob, K1L, al, dp, c0):
    p_ = 1 + dp
    rhoe_ = rhob * p_
    rhosq_ = rho0 * rhoe_
    drhop_ = (rhoe_ - rho0) / rhosq_
    akk = K1L / al
    akxsq_ = -1 / rhosq_ - akk / p_
    akysq_ = akk / p_

    def rhs(s_, v):
        return [v[1], drhop_ + akxsq_ * v[0], v[3], akysq_ * v[2], -v[0] / rho0 - .5 * (v[1]**2 + v[3]**2)]
    sol = solve_ivp(rhs, (0, al), list(c0) + [0.0], rtol=1e-13, atol=1e-17, method='DOP853')
    return sol.y[:, -1]


c0 = (0.002, 0.001, 0.001, -0.0005)
print('  (rho0 = rhob = 10 m、L = 1.3 m、δ = 0.01、振幅 2 mm・1 mrad。ODE は x\' = q、q\' = drhop + akxsq x、y\' = qy、qy\' = akysq y、z\' = -x/ρ0 - (q^2 + qy^2)/2。dv の項は別)')
print('  %-34s %-12s %-12s %-12s %-12s %s' % ('場合', '|Δx|', '|Δq|', '|Δy|', '|Δqy|', '|Δz| (K64 の式 / 符号を直した式)'))
for K1L, name in ((0.0, 'K1 = 0 (akxsq < 0)'), (0.065, 'K1L = 0.065 (集束、akxsq < 0)'), (-0.065, 'K1L = -0.065 (発散、akxsq > 0)')):
    ref = ode_h2(10.0, 10.0, K1L, 1.3, 0.01, c0)
    out = []
    for k64 in (True, False):
        r = fortran_body(10.0, 10.0, K1L, 1.3, 0.01, c0, k64)
        out.append(r)
    dd_ = [abs(out[1][i] - ref[i]) for i in range(4)]
    print('  %-34s %-12.1e %-12.1e %-12.1e %-12.1e %.2e / %.2e' % (name, dd_[0], dd_[1], dd_[2], dd_[3], abs(out[0][4] - ref[4]), abs(out[1][4] - ref[4])))
print('  → x, q, y, qy は 1e-17 台で一致 (H2 の厳密解そのもの)。z は akxsq > 0 の枝は K64 の式のままで一致するが、akxsq < 0 の枝 (K1 ≥ 0) は K64 の式だと δ L^3/(3ρ^2) 程度ずれ、`xsxkx` の符号を反転すると 1e-18 で一致 (docs/30 の B250)。')

print('\n== (e) pysad の track6._tbendiinit・_tbendibody (K1 ≠ 0 の本体 1 刻み) ==')
from pysad import track6, compat
for mode in ('k64', 'v2'):
    old = C.set_mode(mode)
    for K1L in (0.065, -0.065):
        m = types.SimpleNamespace(dp=0.01, rho0=10.0, rhoe=10.0 * 1.01, xi=c0[0], pxi=c0[1], yi=c0[2], pyi=c0[3], zi=0.0)
        track6._tbendiinit(m, K1L, 1.3)
        track6._tbendibody(m, 1.3)
        ref = ode_h2(10.0, 10.0, K1L, 1.3, 0.01, c0)
        print('  %s K1L=%+.3f: |Δx| %.1e |Δq| %.1e |Δy| %.1e |Δqy| %.1e |Δz| %.1e' % (mode, K1L, abs(m.xi - ref[0]), abs(m.pxi - ref[1]), abs(m.yi - ref[2]), abs(m.pyi - ref[3]), abs(m.zi - ref[4])))
    C.set_mode('v2')

print('\n== (f) H2 の運動方程式と元の H の運動方程式の差 (展開の次数) ==')
# 元の H (K1 = 0、ρb = ρ0 = 10 m、δ = 0)。ODE: x' = (1+hx) px/Pz、px' = -hb(1+hx) + h Pz、y' = (1+hx) py/Pz、py' = 0、Pz = sqrt(p^2 - px^2 - py^2)
rho = 10.0
hh = 1 / rho


def exact_h(c, al, dp=0.0):
    pp_ = 1 + dp

    def rhs(s_, v):
        xx, pxx, yy, pyy = v
        Pz = math.sqrt(pp_**2 - pxx**2 - pyy**2)
        return [(1 + hh * xx) * pxx / Pz, -hh * (1 + hh * xx) + hh * Pz, (1 + hh * xx) * pyy / Pz, 0.0]
    return solve_ivp(rhs, (0, al), c, rtol=1e-13, atol=1e-17, method='DOP853').y[:, -1]


def lin_h2(c, al, dp=0.0):
    pp_ = 1 + dp
    r = ode_h2(rho, rho, 0.0, al, dp, (c[0], c[1] / pp_, c[2], c[3] / pp_))
    return np.array([r[0], r[1] * pp_, r[2], r[3] * pp_])


print('  厳密な H (根号つき) の写像と H2 の写像の差 (L = 3 m、δ = 0。(x, px, y, py) = a (1, 0.5, 1, 0.5)):')
prev = None
for a in (1e-3, 2e-3, 4e-3, 8e-3):
    c = np.array([a, 0.5 * a, a, 0.5 * a])
    de = np.abs(exact_h(c, 3.0) - lin_h2(c, 3.0)).max()
    print('    a = %.0f mm: 最大の差 %.2e%s' % (a * 1e3, de, '' if prev is None else '、前の振幅との比 %.2f' % (de / prev)))
    prev = de
print('  → 比 4 = 振幅の 2 乗。H2 が落とした項 (H の 3 次の項 x(px^2+py^2)/(2ρ) など) は、運動方程式では振幅の 2 次の力になる。')

print('\n== (g) 設計軌道の行列の分散・R56 の項 (扇形、L = 1 m、ANGLE = 0.1 rad、電子 3 GeV) ==')
from pysad import lattice as lat
ffs, el = C.make_bend(1.0, 0.1, 'E1=0 E2=0')
Mc, _ = lat._el_cod_step(el, np.zeros(6), 0.0)
th = 0.1
rho_ = 1.0 / th
p0c = 3e9 / 0.51099895e6
gam0 = math.sqrt(1 + p0c**2)
print('  式: R16 = ρ(1 - cos θ) = %.9f、R26 = sin θ = %.9f、R56 = -(L - ρ sin θ) + L/γ0^2 = %.9e (γ0 = %.2f)' % (rho_ * (1 - math.cos(th)), math.sin(th), -(1.0 - rho_ * math.sin(th)) + 1.0 / gam0**2, gam0))
hh_ = 1e-6
dz_dd = (C.track_one(ffs, [0, 0, 0, 0, 0, hh_])[4] - C.track_one(ffs, [0, 0, 0, 0, 0, -hh_])[4]) / (2 * hh_)
print('  CALC の行列 (lattice._el_cod_step): M[x,δ] = %.9f、M[px,δ] = %.9f。R56 は追跡の中心差分: ∂z/∂δ = %.9e' % (Mc[0, 5], Mc[1, 5], dz_dd))
print('  R11 = R22 = cos θ = %.9f (CALC %.9f %.9f)、R12 = ρ sin θ = %.9f (CALC %.9f)、R21 = -sin θ/ρ = %.9f (CALC %.9f)、R33 = 1 (CALC %.3f)、R34 = L = 1 (CALC %.9f)' % (
    math.cos(th), Mc[0, 0], Mc[1, 1], rho_ * math.sin(th), Mc[0, 1], -math.sin(th) / rho_, Mc[1, 0], Mc[2, 2], Mc[2, 3]))
