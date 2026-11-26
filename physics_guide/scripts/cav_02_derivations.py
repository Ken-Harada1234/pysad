# -*- coding: utf-8 -*-
"""空洞の章 (10_cavity.md) の補完。

(a) 端の項 H_f から ΔE = ∂_t H_f、Δp_x = -∂_x H_f を取り、K64 の tcavfrin の式 (dh、dpt) と比べる (sympy、§7.2)。
(b) tcavfrin を numpy で組み直し (入口の端だけ)、運動学的な傾き q の再スケール q *= p1r/p2r の有無でシンプレクティック性
    |J^T S J - S| を比べる。K64 の tcavfrin は再スケールが無い (本体 tcav の px=(px*p1r+dpx)/p2r との対比)。
(c) pysad の CAVI (フリンジあり) の追跡の同じ量 (端付きの空洞全体)。k64 (K64 の動作) と v2 (B330 の 1/β と B332 の q の再スケールを直した版) の両方。
(d) 端の収束 Δpx の β 依存 (§7.1 (f) の見直し): 端だけを取り出して、K64 の式と比べる。k64 と v2 の両方。
(e) 重み ws: wl = ln(1+r)、dpepe の導出 (§6.2)。
"""
import io
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6
from pysad import compat

MC2 = 0.51099895e6
c_light = 299792458.0

# ---------------------------------------------------------------------------------------------------- (a)
print('== (a) 端の項 H_f から tcavfrin の式を作る (sympy) ==')
x, y, t, w, L, v, dphi, s0 = sp.symbols('x y t omega L v dphi s0', real=True)
p0 = sp.symbols('p0', positive=True)
# H_f の s の積分 (delta 関数の係数): K = -(v/L)[sin(ωt - dphi) + s0](x² + y²)/4   (s0 = sin dphi - offset)
K = -(v / L) * (sp.sin(w * t - dphi) + s0) * (x**2 + y**2) / 4
dE = sp.diff(K, t)                  # ΔE = ∂K/∂t (§5.1。mc² 単位、c = 1)
dpx = -sp.diff(K, x)                # Δp_x = -∂K/∂x (mc 単位)
vf = v / L / p0 / 2
dh_code = -w * p0 * vf * (x**2 + y**2) / 2 * sp.cos(w * t - dphi)
dpt_num = vf * (sp.sin(w * t - dphi) + s0)       # コード: dpt = vf*(sin(ph - dphis) + s0)/p1r、Δ(px/p0) = x*dpt*p1r
print('  ΔE = ∂_t K = %s' % sp.simplify(dE))
print('  コードの dh = -w p0 vf (x² + y²)/2 cos(ωt - dphis) との差 = %s' % sp.simplify(dE - dh_code))
print('  Δp_x/p0 = -∂_x K / p0 = %s' % sp.simplify(dpx / p0))
print('  コードの x·vf (sin + s0) との差 = %s' % sp.simplify(dpx / p0 - x * dpt_num))

# ---------------------------------------------------------------------------------------------------- (b)
print('\n== (b) tcavfrin (入口の端だけ) の q の再スケールの有無とシンプレクティック性 ==')


def h_of(p):
    return math.sqrt(1.0 + p * p)


def fringe(c, par, rescale, v_sign=+1):
    """c = (x, px, y, py, z, δ) (px・py は正準運動量を p0 で規格化したもの)。K64 tcavfrin の式。
    rescale = True なら q ← q p1r/p2r (本体 tcav の px=(px*p1r+dpx)/p2r と同じ扱い) を足す"""
    p0_, h0_, w_, dphis, offset, al, v_ = par
    v_ = v_ * v_sign
    X, PX, Y, PY, Z, G = c
    p1r = 1.0 + G
    p1 = p0_ * p1r
    h1 = h_of(p1)
    qx, qy = PX / p1r, PY / p1r
    vf_ = v_ / al / p0_ * .5
    s0_ = math.sin(dphis) - offset
    tt = -Z * h1 / p1
    ph = w_ * tt
    dpt = vf_ * (math.sin(ph - dphis) + s0_) / p1r
    qx2 = qx + X * dpt
    qy2 = qy + Y * dpt
    dh = -w_ * p0_ * vf_ * (X * X + Y * Y) * .5 * math.cos(ph - dphis)
    h2 = h1 + dh
    a = max(dh * (h1 + h2) / p1**2, -1.0)
    dpr = a / (1.0 + math.sqrt(1.0 + a))
    dp2r = G + p1r * dpr
    p2r = 1.0 + dp2r
    Z2 = Z * p2r / p1r * h1 / h2
    if rescale:
        qx2 = qx * p1r / p2r + X * dpt * p1r / p2r
        qy2 = qy * p1r / p2r + Y * dpt * p1r / p2r
    return np.array([X, qx2 * p2r, Y, qy2 * p2r, Z2, dp2r])


S = np.zeros((6, 6))
for i in (0, 2, 4):
    S[i, i + 1] = 1
    S[i + 1, i] = -1


def jac(f, c, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h * max(1.0, abs(c[j])) if j not in (0, 2) else h * 1e-3
        J[:, j] = (f(c + d) - f(c - d)) / (2 * d[j])
    return J


V, f_rf, phi, al = 1.0e5, 500.0e6, 0.5, 0.1
w_ = 2 * math.pi * f_rf / c_light
v_ = V / MC2
print('  CAVI (L = %.2f m、VOLT = %.0e V、FREQ = 500 MHz、PHI = %.1f)、PHI から dphis = -PHI (trf0 = 0)、offset = sin(dphis)' % (al, V, phi))
print('  |J^T S J - S| (中心差分)  [入口の端だけ。(x, px) = (1 mm または 10 mm, 0.01)、z = 0.02 m]')
rows = {}
for lab, pe in (('3 GeV/c', 3.0e9), ('2 MeV/c', 2.0e6), ('0.6 MeV/c', 0.6e6)):
    p0_ = pe / MC2
    par = (p0_, h_of(p0_), w_, -phi, math.sin(-phi), al, v_)
    for xv in (1e-3, 1e-2):
        c0 = np.array([xv, 0.01, 0.0, 0.0, 0.02, 1e-3])
        e0 = np.abs((lambda J: J.T @ S @ J - S)(jac(lambda cc: fringe(cc, par, False), c0))).max()
        e1 = np.abs((lambda J: J.T @ S @ J - S)(jac(lambda cc: fringe(cc, par, True), c0))).max()
        print('   %-10s x = %4.0f mm: q の再スケールなし (K64 の式) %.1e、あり %.1e' % (lab, xv * 1e3, e0, e1))
print('  (中心差分の精度が 1e-10 〜 1e-11 なので、それ以下は 0 と読む。)')

# ---------------------------------------------------------------------------------------------------- (c)
print('\n== (c) pysad の CAVI 全体 (入口・出口の端つき、L = 0.1 m) のヤコビアン。k64 (K64 の動作) と v2 ==')


def mk(el, momentum):
    src = 'MOMENTUM = %s;\n%s\nMARK ST = ();\nLINE ONE = (ST CV);\nFFS;\nUSE ONE;\n' % (momentum, el)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def trk(ffs, c):
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6.track_particles(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


for mode in ('k64', 'v2'):
    compat.set_mode(mode)
    print('  -- 互換モード %s' % mode)
    for lab, mom in (('3 GeV/c', '3GEV'), ('2 MeV/c', '2MEV'), ('0.6 MeV/c', '600KEV')):
        for fr in ('', 'DISFRIN=1'):
            ffs = mk('CAVI CV=(L=%g VOLT=%g FREQ=%g PHI=%g %s);' % (al, V, f_rf, phi, fr), mom)
            for xv in (1e-3, 1e-2):
                c0 = np.array([xv, 0.01, 0.0, 0.0, 0.02, 1e-3])
                J = jac(lambda cc: trk(ffs, cc), c0)
                print('   %-10s %-10s x = %4.0f mm: |J^T S J - S| = %.1e' % (lab, 'フリンジあり' if not fr else 'DISFRIN=1', xv * 1e3, np.abs(J.T @ S @ J - S).max()))
compat.set_mode(None)

# ---------------------------------------------------------------------------------------------------- (d)
print('\n== (d) 端の収束の β 依存 (§7.1 (f) の見直し): pysad の入口の端だけ (FRINGE = 1 と DISFRIN = 1 の差) ==')
print('  本体のキックが 0 になる位相 (cos(ωt - dphis) = 0、すなわち ωt - dphis = π/2) に粒子を置く。')
print('  そのとき K64 の端の式は Δp_x = x v/(2 L p0) (sin(ωt - dphis) = 1、s0 = 0)。δ = 0、x = 1 mm、p_x = 0。')
for mode in ('k64', 'v2'):
    compat.set_mode(mode)
    print('  -- 互換モード %s' % mode)
    for lab, mom, pe in (('3 GeV/c', '3GEV', 3.0e9), ('2 MeV/c', '2MEV', 2.0e6), ('0.6 MeV/c', '600KEV', 0.6e6)):
        p0_ = pe / MC2
        beta = p0_ / h_of(p0_)
        tt = (math.pi / 2 - phi) / w_              # ωt - dphis = π/2、dphis = -PHI
        z0 = -beta * tt
        c0 = np.array([1e-3, 0.0, 0.0, 0.0, z0, 0.0])
        o1 = trk(mk('CAVI CV=(L=%g VOLT=%g FREQ=%g PHI=%g FRINGE=1);' % (al, V, f_rf, phi), mom), c0)
        o0 = trk(mk('CAVI CV=(L=%g VOLT=%g FREQ=%g PHI=%g DISFRIN=1);' % (al, V, f_rf, phi), mom), c0)
        dpx_end = o1[1] - o0[1]
        k64 = c0[0] * v_ / al / p0_ * .5
        print('   %-10s β = %.4f: Δp_x(端) = %.5e、K64 の式 = %.5e、比 = %.6f、(1/β をかけた式との比 = %.4f)' % (lab, beta, dpx_end, k64, dpx_end / k64, dpx_end / (k64 / beta)))
compat.set_mode(None)
print('  (本体のキックが 0 の位相では、k64 の比は 1 (エネルギーによらない。1/β は入っていない)、v2 は 1/β (粒子の β) の値に一致する (B330)。')
print('   §7.1 の (f) の k64 の 0.9993・0.9734 は、端のキックで動いた x に、本体の x 依存のキック (V20・V02 の補正) が反応する分が混ざったもの。)')

# ---------------------------------------------------------------------------------------------------- (e)
print('\n== (e) 重みの wl と dpepe (§6.2) ==')
r = sp.symbols('r', positive=True)
wl = sp.asinh(r * (2 + r) / (2 * (1 + r)))
for rv in (1e-3, 0.05, 0.5, 3.0):
    print('   r = %g: asinh(r(2+r)/(2(1+r))) = %.15f、ln(1 + r) = %.15f' % (rv, float(wl.subs(r, rv)), math.log(1 + rv)))
h0s, ves, pes, p0s = sp.symbols('h0 v_n pe p0', positive=True)
hes = h0s + ves
# pe² - p0² = he² - h0² = v_n (he + h0)  → (pe - p0)/pe = v_n (he + h0)/((pe + p0) pe) = dpepe
lhs = (pes - p0s) / pes
rhs = ves * (hes + h0s) / ((pes + p0s) * pes)
num = sp.simplify(((pes**2 - p0s**2) - ves * (hes + h0s)))
print('   pe² - p0² = v_n (he + h0) を使うと (pe - p0)/pe = v_n (he + h0)/((pe + p0) pe): 条件 pe² - p0² - v_n(he + h0) = %s (he² - 1 = pe²、h0² - 1 = p0² のとき 0)' % sp.simplify(num.subs(pes, sp.sqrt(hes**2 - 1)).subs(p0s, sp.sqrt(h0s**2 - 1)).subs(ves, hes - h0s)))
