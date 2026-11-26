# -*- coding: utf-8 -*-
"""09 章 (09_solenoid.md) の途中式の確認 (本文に書いた中間式を 1 つずつ確かめる)。

(A) §4.2: 4 元連立 (x0, y0, q_xm, q_ym) の Cramer の解 (行列式 -V/b_p)、A, B, C, D、a, bw, c, dw。
(B) §3: K=0 の検算 (ζ = x + iy、回転系)。
(C) §5.3: (9.8a) を (9.6) の微分に入れた x', y' の式 (9.8b)、10 個の単項式の係数を、手で導いた閉じた形・数値積分・K64 の tztafs と比較。
(D) §5.2: 交差項の 4 つの積分の微分 (V 倍が被積分関数、s=0 の定数項)。
(E) §8.2: R2'(-χ2) R1'(-χ1) の積の成分、出口の成分式、運動量の回転。
"""
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)

# ------------------------------------------------------------------------------------------------ (A)
print('== (A) §4.2 の 4 元連立の解 ==')
w1, w2 = sp.symbols('w1 w2', positive=True)
wp, wm = w1 + w2, w1 - w2
V = w1**2 + w2**2
bp = sp.sqrt(wp * wm)                         # b_p^2 = w+ w-
wd = bp / wp                                    # w_d = b_p / w+
x0, y0, qx, qy = sp.symbols('x0 y0 qxm qym')
A_, B_, C_, D_ = sp.symbols('A B C D')
eqs = [sp.Eq(x0, A_ + C_), sp.Eq(y0, wd * B_ + D_ / wd), sp.Eq(qx, w1 * B_ + w2 * D_), sp.Eq(qy, -wd * w1 * A_ + w2 / wd * C_)]
sol = sp.solve(eqs, [A_, B_, C_, D_], dict=True)[0]
a_f = (w2 * wp * x0 - bp * qy) / V
bw_f = (wp * qx - bp * w2 * y0) / V
c_f = (w1 * wd * x0 + qy) / V
dw_f = (-wd * qx + w1 * y0) / V
for nm, got, ref in (('A - a', sol[A_], a_f), ('B - bw', sol[B_], bw_f), ('C - bp c', sol[C_], bp * c_f), ('D - bp dw', sol[D_], bp * dw_f)):
    print('  %-10s = %s' % (nm, sp.simplify(got - ref)))
det2 = sp.simplify(wd * w2 - w1 / wd)
print('  y0, q_xm の 2×2 の行列式 w_d w2 - w1/w_d = %s、-V/b_p との差 = %s' % (det2, sp.simplify(det2 + V / bp)))
print('  w_d^2 w1 + w2 - V/w+ = %s' % sp.simplify(wd**2 * w1 + w2 - V / wp))
print('  V - w2 w+ - w1 w- = %s' % sp.simplify(V - w2 * wp - w1 * wm))
print('  w1 w- x0 + bp q_ym = bp (w1 w_d x0 + q_ym): 差 = %s' % sp.simplify(w1 * wm * x0 + bp * qy - bp * (w1 * wd * x0 + qy)))

# ------------------------------------------------------------------------------------------------ (B)
print('\n== (B) K=0 の検算 ==')
s = sp.symbols('s', real=True)
bpp = sp.symbols('b_p', positive=True)
al, be = sp.symbols('alpha beta')
zeta = al + be * sp.exp(-sp.I * bpp * s)
print('  ζ = α + β e^{-i b_p s}: ζ\'\' + i b_p ζ\' = %s' % sp.simplify(sp.diff(zeta, s, 2) + sp.I * bpp * sp.diff(zeta, s)))
eta = sp.exp(sp.I * bpp * s / 2) * zeta
print('  η = e^{i b_p s/2} ζ: η\'\' + (b_p/2)^2 η = %s (ζ\'\' = -i b_p ζ\' を使って)' % sp.simplify(sp.diff(eta, s, 2) + (bpp / 2)**2 * eta))
xx, yy = sp.Function('x')(s), sp.Function('y')(s)
zz = xx + sp.I * yy
eq_zeta = sp.expand(sp.diff(zz, s, 2) - (bpp * sp.diff(yy, s) - sp.I * bpp * sp.diff(xx, s)) * 0)
print('  (9.4) (K=0): x\'\' = b_p y\', y\'\' = -b_p x\'  →  ζ\'\' = b_p y\' - i b_p x\' = -i b_p ζ\': 恒等式 (b_p y\' - i b_p x\') + i b_p (x\' + i y\') = %s' % sp.simplify(bpp * sp.diff(yy, s) - sp.I * bpp * sp.diff(xx, s) + sp.I * bpp * (sp.diff(xx, s) + sp.I * sp.diff(yy, s))))

# ------------------------------------------------------------------------------------------------ (C)
print('\n== (C) x\', y\' の (9.8b) と 10 個の単項式 ==')
mp.mp.dps = 30
W1, W2, L = mp.mpf('0.9'), mp.mpf('0.35'), mp.mpf('0.4')
WP, WM = W1 + W2, W1 - W2
VV = W1**2 + W2**2
BP = mp.sqrt(WP * WM)
WD = BP / WP
KP = W1 * W2


def basis(sv):
    S, Cc, Sh, Ch = mp.sin(W1 * sv), mp.cos(W1 * sv), mp.sinh(W2 * sv), mp.cosh(W2 * sv)
    return S, Cc, Sh, Ch


def xyprime(sv, awu, dwu, pxi, pyi):
    S, Cc, Sh, Ch = basis(sv)
    xp = awu * (-WP * S + WM * Sh) + dwu * BP * (Ch - Cc) + pxi * Cc + WD * pyi * Sh
    yp = awu * BP * (Ch - Cc) + dwu * (WM * S + WP * Sh) - WD * pxi * S + pyi * Ch
    return xp, yp


def xyprime_ABCD(sv, awu, dwu, pxi, pyi):
    A = WP / W1 * awu
    B = (pxi - BP * dwu) / W1
    Cq = WM / W2 * awu + WD / W2 * pyi
    D = BP / W2 * dwu
    S, Cc, Sh, Ch = basis(sv)
    xp = W1 * (-A * S + B * Cc) + W2 * (Cq * Sh + D * Ch)
    yp = WD * W1 * (-A * Cc - B * S) + W2 / WD * (Cq * Ch + D * Sh)
    return xp, yp


import random
random.seed(2)
v = [mp.mpf(random.uniform(-1, 1)) for _ in range(4)]
d1 = max(abs(a - b) for a, b in zip(xyprime(mp.mpf('0.27'), *v), xyprime_ABCD(mp.mpf('0.27'), *v)))
print('  (9.8b) の x\', y\' と、(9.8a) を (9.6) の微分に入れたものの差 (s = 0.27) = %s' % mp.nstr(d1, 3))
print('  s = 0: x\'(0) - pxi = %s、y\'(0) - pyi = %s' % tuple(mp.nstr(a - b, 3) for a, b in zip(xyprime(mp.mpf(0), *v), (v[2], v[3]))))


def Qexact(awu, dwu, pxi, pyi):
    return -mp.quad(lambda sv: sum(c**2 for c in xyprime(sv, awu, dwu, pxi, pyi)), [0, L]) / 2


s1, c1, sh, ch = mp.sin(W1 * L), mp.cos(W1 * L), mp.sinh(W2 * L), mp.cosh(W2 * L)
p1, p2 = W1 * L, W2 * L
dc1, dch2, xs1, xsh2 = c1 - 1, ch - 1, p1 - s1, p2 - sh
akkp = KP
aw1 = L + (sh * ch - p2) * W1 / W2 / WP
aw2 = L + (s1 * c1 - p1) * W2 / W1 / WP
cr2 = c1 * W2 / W1
cr3 = ch * W1 / W2
csw1 = akkp * L * (-dc1 + dch2)
csws = (ch * W1 * p1 + c1 * W2 * p2) / VV
ca1 = -ch * dc1 - dch2 - s1 * sh * WM / WP
dxs = W2 * xs1 - W1 * xsh2
dcw1 = dc1 * p2
dcw2 = dch2 * p1
g1 = (s1 * W2)**2 / akkp
g2 = -(sh * W1)**2 / akkp
bzp = BP


def tztafs(zwu, pxt, pyt, zw1, wst, w12t, swss, g1t):
    w12xsh = w12t * xsh2
    wsxs = wst * xs1
    cwsd = c1 * w12xsh - ch * wsxs
    cswa = (csw1 + cwsd) * swss
    wsm = 4 * swss * akkp
    ws1 = (-1 + wsm) * wst
    ws2 = w12t * (1 + wsm)
    cr1 = bzp**2 * swss
    f1 = cswa + csws - zw1
    f2 = ca1 * cr1 + g1t
    f3 = 2 * cr1 * cwsd + cr3 * w12xsh + dcw1 * ws1 + dcw2 * ws2 - cr2 * wsxs
    return zwu * (zwu * f3 / 2 + f2 * pxt + bzp * f1 * pyt) - zw1 * pyt**2 / 2


def Qk64(awu, dwu, pxi, pyi):
    dz1 = tztafs(awu, pxi, pyi, aw1, WP, WM, 1 / VV, g1)
    dz2 = tztafs(-dwu, -pyi, pxi, aw2, -WM, WP, -1 / VV, g2)
    return bzp * (-awu * dwu * dxs**2 / akkp + ca1 * pxi * pyi / VV) + dz1 + dz2


E = WM * s1 * sh + WP * (c1 * ch - 1)
f3_1 = -(2 * W1 * W2 + W1**2 - W2**2) * L + W2 * WP / W1 * s1 * c1 - W1 * WM / W2 * sh * ch - 2 * WP * WM / VV * (WM * c1 * sh - WP * s1 * ch)
f3_2 = (2 * W1 * W2 - W1**2 + W2**2) * L - W2 * WM / W1 * s1 * c1 - W1 * WP / W2 * sh * ch + 2 * WP * WM / VV * (WP * c1 * sh + WM * s1 * ch)
hand = {
    'pxi^2': -aw2 / 2,
    'pyi^2': -aw1 / 2,
    'awu^2': f3_1 / 2,
    'dwu^2': f3_2 / 2,
    'awu pxi': W2 * s1**2 / W1 + BP**2 * ca1 / VV,
    'dwu pyi': -(W1 * sh**2 / W2 + BP**2 * ca1 / VV),
    'awu pyi': BP * (-aw1 - (WM * c1 * sh - WP * s1 * ch) / VV),
    'dwu pxi': BP * (aw2 - (WP * c1 * sh + WM * s1 * ch) / VV),
    'awu dwu': -BP / KP * (W1 * sh - W2 * s1)**2,
    'pxi pyi': BP * ca1 / VV,
}
unit = {'awu': 0, 'dwu': 1, 'pxi': 2, 'pyi': 3}


def vec(*idx):
    t = [mp.mpf(0)] * 4
    for i in idx:
        t[i] += 1
    return t


print('  単項式の係数: 厳密 (数値積分) / 手で導いた閉じた形 / K64 (tztafs) と、その差')
worst = mp.mpf(0)
for name, hv in hand.items():
    parts = name.replace('^2', '').split()
    if '^2' in name:
        i = unit[parts[0]]
        ex = Qexact(*vec(i))
        k6 = Qk64(*vec(i))
    else:
        i, j = unit[parts[0]], unit[parts[1]]
        ex = Qexact(*vec(i, j)) - Qexact(*vec(i)) - Qexact(*vec(j))
        k6 = Qk64(*vec(i, j)) - Qk64(*vec(i)) - Qk64(*vec(j))
    worst = max(worst, abs(ex - hv), abs(hv - k6))
    print('  %-8s: 厳密 %s、手 %s、K64 %s (手 - 厳密 = %s、手 - K64 = %s)' % (name, mp.nstr(ex, 14), mp.nstr(hv, 14), mp.nstr(k6, 14), mp.nstr(hv - ex, 3), mp.nstr(hv - k6, 3)))
print('  最大の差 = %s' % mp.nstr(worst, 3))
print('  補助の関係: E = w- s1 sh + w+ (c1 ch - 1) = -w+ ca1: 差 = %s、b_p^2 = w+ w-: 差 = %s、cswa + csws の和 (call 1) と -(w- c1 sh - w+ s1 ch)/V の差 = %s' % (
    mp.nstr(E + WP * ca1, 3), mp.nstr(BP**2 - WP * WM, 3),
    mp.nstr(((csw1 + c1 * WM * xsh2 - ch * WP * xs1) / VV + csws) + (WM * c1 * sh - WP * s1 * ch) / VV, 3)))
print('  w_d b_p = w-: 差 = %s、w-/w_d = b_p: 差 = %s、b_p/w_d = w+: 差 = %s' % (mp.nstr(WD * BP - WM, 3), mp.nstr(WM / WD - BP, 3), mp.nstr(BP / WD - WP, 3)))
print('  dxs = w1 sh - w2 s1 (xs1 = φ1 - s1、xsh2 = φ2 - sh、w1 w2 ℓ が打ち消す): 差 = %s' % mp.nstr(dxs - (W1 * sh - W2 * s1), 3))

# ------------------------------------------------------------------------------------------------ (D)
print('\n== (D) §5.2 の 4 つの交差積分 ==')
S_, Cc_, Sh_, Ch_ = sp.sin(w1 * s), sp.cos(w1 * s), sp.sinh(w2 * s), sp.cosh(w2 * s)
tab = {
    'cos cosh': ((w1 * S_ * Ch_ + w2 * Cc_ * Sh_) / V, Cc_ * Ch_),
    'cos sinh': ((w1 * S_ * Sh_ + w2 * Cc_ * Ch_) / V, Cc_ * Sh_),
    'sin cosh': ((-w1 * Cc_ * Ch_ + w2 * S_ * Sh_) / V, S_ * Ch_),
    'sin sinh': ((-w1 * Cc_ * Sh_ + w2 * S_ * Ch_) / V, S_ * Sh_),
}
for nm, (F, f) in tab.items():
    print('  %-9s: dF/ds - f = %s、F(0) = %s' % (nm, sp.simplify(sp.diff(F, s) - f), sp.simplify(F.subs(s, 0))))

# ------------------------------------------------------------------------------------------------ (E)
print('\n== (E) §8.2 の行列の積と出口の成分 ==')
c1_, s1_, c2_, s2_ = sp.symbols('c1 s1 c2 s2')
R1p = sp.Matrix([[c1_, 0, -s1_], [0, 1, 0], [s1_, 0, c1_]])
R2p = sp.Matrix([[1, 0, 0], [0, c2_, -s2_], [0, s2_, c2_]])
R1m = R1p.subs(s1_, -s1_)
R2m = R2p.subs(s2_, -s2_)
P = sp.expand(R2m * R1m)
print('  R2\'(-χ2) R1\'(-χ1) = %s' % P.tolist())
print('  期待 [[c1,0,s1],[-s1 s2,c2,s2 c1],[-s1 c2,-s2,c1 c2]] との差 = %s' % sp.simplify(P - sp.Matrix([[c1_, 0, s1_], [-s1_ * s2_, c2_, s2_ * c1_], [-s1_ * c2_, -s2_, c1_ * c2_]])).tolist())
qx_, qy_, qz_ = sp.symbols('qx qy qz')
qz1 = s1_ * qx_ + c1_ * qz_
qvec = R2p * R1p * sp.Matrix([qx_, qy_, qz_])
print('  運動量の回転 (入口): q_z1 = s1 qx + c1 qz、q_x\' = c1 qx - s1 qz、q_y\' = c2 qy - s2 q_z1: 差 = %s' % sp.simplify(
    sp.Matrix([c1_ * qx_ - s1_ * qz_, c2_ * qy_ - s2_ * qz1, s2_ * qy_ + c2_ * qz1]) - qvec).tolist())
x_, y_, s0, dz_ = sp.symbols('x y s0s dz')
vin = sp.Matrix([x_, y_, s0 - dz_])
r1 = R1p * vin
r2 = R2p * r1
ds1 = sp.simplify(r1[2] - s0)
ds2 = sp.simplify(r2[2] - s0)
code_ds1 = s1_ * x_ + (c1_ - 1) * s0 - c1_ * dz_
code_ds2 = s2_ * y_ + c2_ * code_ds1 + (c2_ - 1) * s0
print('  入口: ds1 - コード = %s、ds2 - コード = %s、x1 - コード = %s、y\' - コード = %s' % (
    sp.simplify(ds1 - code_ds1), sp.simplify(ds2 - code_ds2), sp.simplify(r1[0] - (c1_ * x_ - s1_ * (s0 - dz_))),
    sp.simplify(r2[1] - (c2_ * y_ - s2_ * (code_ds1 + s0)))))
# 出口: (R2' R1')^T を (x, y, -s0s) に施す
vout = sp.Matrix([x_, y_, -s0])
u1 = R2p.T * vout
u2 = R1p.T * u1
code_y1 = c2_ * y_ - s2_ * s0
code_ds1x = -s2_ * y_ - (c2_ - 1) * s0
code_x1 = c1_ * x_ + s1_ * (code_ds1x - s0)
code_ds2x = -s1_ * x_ + c1_ * code_ds1x - (c1_ - 1) * s0
print('  出口: y1 - コード = %s、(z1 + s0s) - ds1 = %s、x1 - コード = %s、(z2 + s0s) - ds2 = %s' % (
    sp.simplify(u1[1] - code_y1), sp.simplify(u1[2] + s0 - code_ds1x), sp.simplify(u2[0] - code_x1), sp.simplify(u2[2] + s0 - code_ds2x)))
pz0 = sp.Symbol('pz0')
qq = R1p.T * (R2p.T * sp.Matrix([qx_, qy_, pz0]))
pz1o = -s2_ * qy_ + c2_ * pz0
print('  出口の運動量: pz1 = -s2 qy + c2 qz、q_y\' = c2 qy + s2 qz、q_x\' = c1 qx + s1 pz1、pz2 = -s1 qx + c1 pz1: 差 = %s' % sp.simplify(
    sp.Matrix([c1_ * qx_ + s1_ * pz1o, c2_ * qy_ + s2_ * pz0, -s1_ * qx_ + c1_ * pz1o]) - qq).tolist())
