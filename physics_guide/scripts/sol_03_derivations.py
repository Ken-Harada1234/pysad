# -*- coding: utf-8 -*-
r"""ソレノイドの章 (09_solenoid.md) の展開の補完。

(a) H2 を P で展開した形 (§3)。
(b) 初期条件から (A, B, C, D) を解く (§4.2) と、K64 の a, bw, c, dw との対応 (sympy、恒等式)。
(c) z の更新 (9.8) の積分を 4 つのモードで閉じた形にし (§5)、K64 の tztafs (係数 aw1, aw2, g1, g2, ca1, cr2, cr3, csw1, csws, dcw1, dcw2, dxs)
    が、その単項式ごとの係数と一致することを sympy の恒等式で示す。
(d) tsolqu の ΔH の運動量の回転 (pyi の更新順序) の写像の行列と行列式 (付録 B)。
(e) pysad (v2 既定) の SOL + QUAD + SOL の追跡の 6×6 ヤコビアンのシンプレクティック性の欠陥 |J^T S J - S| (付録 B)。
"""
import io
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

# ---------------------------------------------------------------------------------------------------------
print('== (a) H2 を P = p_c - a で展開 (§3) ==')
x, y, px, py, p, b, K = sp.symbols('x y p_x p_y p b K', real=True)
Px, Py = px + b * y / 2, py - b * x / 2
quad = sp.expand((Px**2 + Py**2) / (2 * p) + K / 2 * (x**2 - y**2))
pphi = x * py - y * px
claim = (px**2 + py**2) / (2 * p) - b / (2 * p) * pphi + b**2 / (8 * p) * (x**2 + y**2) + K / 2 * (x**2 - y**2)
print('  (P_x² + P_y²)/(2p) + K(x² - y²)/2 - [ (p_x² + p_y²)/(2p) - (b/2p) p_phi + (b²/8p)(x² + y²) + K(x² - y²)/2 ] = %s' % sp.simplify(quad - claim))
print('  展開: P_x² + P_y² = p_x² + p_y² + b(y p_x - x p_y) + (b²/4)(x² + y²)  → %s' % sp.expand(Px**2 + Py**2 - (px**2 + py**2) + b * pphi - b**2 / 4 * (x**2 + y**2)))

# ---------------------------------------------------------------------------------------------------------
print('\n== (b) 初期条件から (A, B, C, D) を解く (§4.2)==')
w1, w2, wd, ell = sp.symbols('w1 w2 w_d ell', positive=True)
x0, y0, qx, qy = sp.symbols('x0 y0 qxm qym', real=True)   # qxm, qym: 運動学的な傾き P/p (= 正準 q + b_p·(∓y, x)/2)
A, B, Cc, D = sp.symbols('A B C D')
ws_, V = w1 + w2, w1**2 + w2**2
bp = wd * ws_                                # b_p = w_d w_+  (w_d² = (w1 - w2)/(w1 + w2) ⇔ b_p² = w1² - w2²)
sol = sp.solve([A + Cc - x0, wd * B + D / wd - y0, w1 * B + w2 * D - qx, -wd * w1 * A + w2 * Cc / wd - qy], [A, B, Cc, D], dict=True)[0]
wd2 = (w1 - w2) / (w1 + w2)


def red(e):
    """w_d² を (w1 - w2)/(w1 + w2) に置き換えて簡約"""
    e = sp.together(sp.simplify(e))
    n, d = sp.fraction(e)
    n = sp.rem(sp.expand(n), wd**2 - wd2 * 1, wd) if False else n
    return sp.simplify(e.subs(wd, sp.sqrt(wd2)))


a_k = (w2 * ws_ * x0 - bp * qy) / V
bw_k = (ws_ * qx - bp * w2 * y0) / V
c_k = (w1 * wd * x0 + qy) / V
dw_k = (-wd * qx + w1 * y0) / V
for name, mine, k64 in (('A = a', sol[A], a_k), ('B = bw', sol[B], bw_k), ('C = b_p c', sol[Cc], bp * c_k), ('D = b_p dw', sol[D], bp * dw_k)):
    print('  %s: 解 - K64 の式 = %s' % (name, red(mine - k64)))
print('  解 A = %s' % sp.factor(sol[A]))
print('  (w_d² = (w1 - w2)/(w1 + w2) を代入して簡約した差)')

# ---------------------------------------------------------------------------------------------------------
print('\n== (c) z の積分の閉じた形と tztafs の係数 (§5) ==')
s = sp.symbols('s', positive=True)
Ap, Bp, Cp, Dp = sp.symbols('Ap Bp Cp Dp')
xp = w1 * (-Ap * sp.sin(w1 * s) + Bp * sp.cos(w1 * s)) + w2 * (Cp * sp.sinh(w2 * s) + Dp * sp.cosh(w2 * s))
yp = wd * w1 * (-Ap * sp.cos(w1 * s) - Bp * sp.sin(w1 * s)) + (w2 / wd) * (Cp * sp.cosh(w2 * s) + Dp * sp.sinh(w2 * s))
integrand = sp.expand(xp**2 + yp**2)
terms = sp.Poly(integrand, Ap, Bp, Cp, Dp)
# 各単項式の係数 (s の関数) を積分する
S1, C1, SH, CH = sp.symbols('s1 c1 sh2 ch2')   # sin φ1, cos φ1, sinh φ2, cosh φ2 (φ = w ℓ)


def to_sym(e):
    e = sp.expand(sp.expand_trig(e))
    e = e.subs({sp.sin(w1 * ell): S1, sp.cos(w1 * ell): C1, sp.sinh(w2 * ell): SH, sp.cosh(w2 * ell): CH})
    return e


# 閉じた形の原始関数 (微分して確かめる)
prim = {}
for name, f in (('cc', sp.cos(w1 * s)**2), ('ss', sp.sin(w1 * s)**2), ('sc', sp.sin(w1 * s) * sp.cos(w1 * s)),
                ('chch', sp.cosh(w2 * s)**2), ('shsh', sp.sinh(w2 * s)**2), ('shch', sp.sinh(w2 * s) * sp.cosh(w2 * s)),
                ('c_ch', sp.cos(w1 * s) * sp.cosh(w2 * s)), ('c_sh', sp.cos(w1 * s) * sp.sinh(w2 * s)),
                ('s_ch', sp.sin(w1 * s) * sp.cosh(w2 * s)), ('s_sh', sp.sin(w1 * s) * sp.sinh(w2 * s))):
    prim[name] = (f, sp.integrate(f, s))
# 直接 sympy に積分させ、微分して元に戻ることを確認 (原始関数の形は sympy 任せ)
tot = 0
for mon, coef in terms.terms():
    mono = Ap**mon[0] * Bp**mon[1] * Cp**mon[2] * Dp**mon[3]
    F = sp.integrate(sp.expand(coef), (s, 0, ell))
    tot += mono * F
tot = sp.expand(-tot / 2)
print('  Δz = -(1/2)∫(x'' ² + y'' ²) ds を (A, B, C, D) の 2 次形式として積分した (項数 %d)' % len(sp.Poly(tot, Ap, Bp, Cp, Dp).terms()))
# 数値に落として ODE の直接積分と比べる
vals = {w1: 0.9, w2: 0.35, ell: 0.4, Ap: 0.003, Bp: 0.0007, Cp: -0.0012, Dp: 0.0004}
vals[wd] = math.sqrt((0.9 - 0.35) / (0.9 + 0.35))
from scipy.integrate import quad as _quad
f_num = sp.lambdify(s, (xp**2 + yp**2).subs(vals), 'numpy')
print('  数値積分 -(1/2)∫(x'' ² + y'' ²) = %.14e、閉じた形 = %.14e' % (-0.5 * _quad(f_num, 0, 0.4, epsabs=1e-16, epsrel=1e-14)[0], float(tot.subs(vals))))


W1E = w2 * (1 + wd**2) / (1 - wd**2)          # w_d² = (w1 - w2)/(w1 + w2) <=> w1 = w2 (1 + w_d²)/(1 - w_d²)


def reduce_zero(expr):
    """w1 を (w_d, w2) で書き、sin²+cos²=1、cosh²-sinh²=1 で簡約した分子 (0 なら恒等式)"""
    e = sp.together(expr.subs(w1, W1E))
    n = sp.expand(sp.numer(e))
    n = sp.expand(sp.rem(n, S1**2 + C1**2 - 1, S1))
    n = sp.expand(sp.rem(n, CH**2 - SH**2 - 1, CH))
    return sp.simplify(n)

# K64 の tztafs と dz の式 (Fortran をそのまま sympy に写す)
wr1, wr2 = w1 / w2, w2 / w1
phi1, phi2 = w1 * ell, w2 * ell
dc1, dch2 = C1 - 1, CH - 1
xs1, xsh2 = phi1 - S1, phi2 - SH              # xsincos: xs = φ - sin φ、xsincosh: xsh = φ - sinh φ
wss = 1 / V
akkp = w1 * w2
bzp = bp
cxs1, cxs2 = dc1 * S1 - xs1, dch2 * SH - xsh2
aw1 = ell + cxs2 * wr1 / ws_
aw2 = ell + cxs1 * wr2 / ws_
g1 = (S1 * w2)**2 / akkp
g2 = -(SH * w1)**2 / akkp
cr2, cr3 = C1 * wr2, CH * wr1
csw1 = akkp * ell * (-dc1 + dch2)
csws = (CH * w1 * phi1 + C1 * w2 * phi2) * wss
ca1 = -(CH * dc1) - dch2 - (S1 * SH * (w1 - w2)) / ws_
dxs = w2 * xs1 - w1 * xsh2
dcw1, dcw2 = dc1 * phi2, dch2 * phi1


def tztafs(zwu, pxt, pyt, zw1, wst, w12t, swss, g1t):
    w12xsh = w12t * xsh2
    wsxs = wst * xs1
    cwsd = C1 * w12xsh - CH * wsxs
    cswa = (csw1 + cwsd) * swss
    wsm = 4 * swss * akkp
    ws1 = (-1 + wsm) * wst
    ws2 = w12t * (1 + wsm)
    cr1 = bzp**2 * swss
    f1 = cswa + csws - zw1
    f2 = ca1 * cr1 + g1t
    f3 = 2 * cr1 * cwsd + cr3 * w12xsh + dcw1 * ws1 + dcw2 * ws2 - cr2 * wsxs
    return zwu * (zwu * f3 / 2 + f2 * pxt + bzp * f1 * pyt) - (zw1 * pyt**2) / 2


# (x0, y0, pxi, pyi) → (a, bw, c, dw)
pxi, pyi = qx, qy
a_ = (w2 * ws_ * x0 - bp * pyi) * wss
bw_ = (ws_ * pxi - bp * w2 * y0) * wss
c_ = (w1 * wd * x0 + pyi) * wss
dw_ = (-wd * pxi + w1 * y0) * wss
awu = a_ / ws_ * w1
dwu = dw_ * w2
dz_code = (bzp * (-((awu * dwu * dxs**2) / akkp) + ca1 * pxi * pyi * wss)
           + tztafs(awu, pxi, pyi, aw1, ws_, (w1 - w2), wss, g1)
           + tztafs(-dwu, -pyi, pxi, aw2, -(w1 - w2), ws_, -wss, g2))
# 厳密な Δz を (x0, y0, pxi, pyi) で
sub_abcd = {Ap: a_, Bp: bw_, Cp: bp * c_, Dp: bp * dw_}
dz_exact = to_sym(tot.subs(sub_abcd))
print('  tztafs (K64 の dz 全体) - 厳密な Δz の分子 (s1² + c1² = 1、ch2² - sh2² = 1、w_d² = (w1 - w2)/(w1 + w2) で簡約) = %s' % reduce_zero(dz_code - dz_exact))

# 単項式ごとの対応 (awu, dwu, pxi, pyi の 2 次形式)
Pa, Pd, Pxx, Pyy = sp.symbols('awu dwu pxi pyi')
# x0, y0 を (awu, dwu, pxi, pyi) で解く
xs_sol = sp.solve([sp.Eq(Pa, awu), sp.Eq(Pd, dwu)], [x0, y0], dict=True)[0]
sub2 = {x0: xs_sol[x0], y0: xs_sol[y0], qx: Pxx, qy: Pyy}
sub2 = {k: v.subs({qx: Pxx, qy: Pyy}) for k, v in sub2.items()}


def in_modes(expr):
    e = expr.subs({x0: xs_sol[x0].subs({qx: Pxx, qy: Pyy}), y0: xs_sol[y0].subs({qx: Pxx, qy: Pyy})}).subs({qx: Pxx, qy: Pyy})
    return sp.expand(e)


ex_modes = in_modes(dz_exact)
code_modes = sp.expand(
    bzp * (-((Pa * Pd * dxs**2) / akkp) + ca1 * Pxx * Pyy * wss)
    + tztafs(Pa, Pxx, Pyy, aw1, ws_, (w1 - w2), wss, g1)
    + tztafs(-Pd, -Pyy, Pxx, aw2, -(w1 - w2), ws_, -wss, g2))
mons = [('awu²', Pa**2), ('dwu²', Pd**2), ('pxi²', Pxx**2), ('pyi²', Pyy**2), ('awu·dwu', Pa * Pd), ('awu·pxi', Pa * Pxx), ('awu·pyi', Pa * Pyy), ('dwu·pxi', Pd * Pxx), ('dwu·pyi', Pd * Pyy), ('pxi·pyi', Pxx * Pyy)]
Pex = sp.Poly(ex_modes, Pa, Pd, Pxx, Pyy)
Pco = sp.Poly(code_modes, Pa, Pd, Pxx, Pyy)
print('  単項式ごとの係数の差 (厳密な積分 - K64):')
for nm, m in mons:
    mon = sp.Poly(m, Pa, Pd, Pxx, Pyy).monoms()[0]
    ce, cc = Pex.coeff_monomial(mon), Pco.coeff_monomial(mon)
    print('    %-8s: 厳密 - K64 の分子の簡約 = %s' % (nm, reduce_zero(ce - cc)))
    # 数値でも
    nv = {w1: 0.9, w2: 0.35, ell: 0.4, wd: math.sqrt((0.9 - 0.35) / (0.9 + 0.35)),
          S1: math.sin(0.9 * 0.4), C1: math.cos(0.9 * 0.4), SH: math.sinh(0.35 * 0.4), CH: math.cosh(0.35 * 0.4)}
    print('             数値 (w1 = 0.9, w2 = 0.35, ℓ = 0.4): 厳密 = %.12e、K64 = %.12e' % (float(ce.subs(nv)), float(cc.subs(nv))))

# 係数の読み: pyi²、pxi² は (積分の定数部分)
print('  読み: pyi² の係数 (K64) = %s' % sp.simplify(Pco.coeff_monomial((0, 0, 0, 2))))
print('        pxi² の係数 (K64) = %s' % sp.simplify(Pco.coeff_monomial((0, 0, 2, 0))))

# ---------------------------------------------------------------------------------------------------------

# (c2) 変数 (awu, dwu, pxi, pyi) で見た (A, B, C, D) と、pxi²・pyi² の項 (単一モードの積分)
print('\n== (c2) (A, B, C, D) を (awu, dwu, pxi, pyi) で書く (§5) ==')
bp_ = wd * ws_
claims = {'A': (ws_ / w1) * Pa, 'B': (Pxx - bp_ * Pd) / w1, 'C': ((w1 - w2) / w2) * Pa + wd / w2 * Pyy, 'D': bp_ * Pd / w2}
real = {'A': a_, 'B': bw_, 'C': bp_ * c_, 'D': bp_ * dw_}
for k in 'ABCD':
    e = real[k].subs({x0: xs_sol[x0].subs({qx: Pxx, qy: Pyy}), y0: xs_sol[y0].subs({qx: Pxx, qy: Pyy})}).subs({qx: Pxx, qy: Pyy})
    print('  %s = %s   [差 = %s]' % (k, claims[k], reduce_zero(sp.together(e - claims[k]))))
# pxi² の項: B = pxi/w1 だけのモード、pyi² の項: C = w_d pyi/w2 だけのモード
Bm = Pxx / w1
Ix = -sp.Rational(1, 2) * sp.integrate(sp.expand((w1 * Bm * sp.cos(w1 * s))**2 + (wd * w1 * Bm * sp.sin(w1 * s))**2), (s, 0, ell))
Cm = wd * Pyy / w2
Iy = -sp.Rational(1, 2) * sp.integrate(sp.expand((w2 * Cm * sp.sinh(w2 * s))**2 + (w2 / wd * Cm * sp.cosh(w2 * s))**2), (s, 0, ell))
Ix_ = to_sym(Ix)
Iy_ = to_sym(Iy)
print('  pxi² の係数 (単一モード B): %s' % sp.simplify(sp.expand(Ix_ / Pxx**2).subs(wd, sp.sqrt(wd2))))
print('  pyi² の係数 (単一モード C): %s' % sp.simplify(sp.expand(Iy_ / Pyy**2).subs(wd, sp.sqrt(wd2))))
print('  K64: -aw1/2 = %s、-aw2/2 = %s' % (sp.simplify(-aw1 / 2), sp.simplify(-aw2 / 2)))
print('  差 (pyi²: -aw1/2 - 単一モード C) = %s、(pxi²: -aw2/2 - 単一モード B) = %s' % (
    reduce_zero(-aw1 / 2 - Iy_ / Pyy**2), reduce_zero(-aw2 / 2 - Ix_ / Pxx**2)))

print('\n== (c3) 積分の閉じた形 (§5.3): 微分して被積分関数に戻るか ==')
Vv = w1**2 + w2**2
cl = {
    'cos^2': (sp.cos(w1 * s)**2, ell / 2 + S1 * C1 / (2 * w1)),
    'sin^2': (sp.sin(w1 * s)**2, ell / 2 - S1 * C1 / (2 * w1)),
    'sin cos': (sp.sin(w1 * s) * sp.cos(w1 * s), S1**2 / (2 * w1)),
    'cosh^2': (sp.cosh(w2 * s)**2, ell / 2 + SH * CH / (2 * w2)),
    'sinh^2': (sp.sinh(w2 * s)**2, -ell / 2 + SH * CH / (2 * w2)),
    'sinh cosh': (sp.sinh(w2 * s) * sp.cosh(w2 * s), SH**2 / (2 * w2)),
    'cos cosh': (sp.cos(w1 * s) * sp.cosh(w2 * s), (w1 * S1 * CH + w2 * C1 * SH) / Vv),
    'cos sinh': (sp.cos(w1 * s) * sp.sinh(w2 * s), (w1 * S1 * SH + w2 * C1 * CH - w2) / Vv),
    'sin cosh': (sp.sin(w1 * s) * sp.cosh(w2 * s), (-w1 * C1 * CH + w2 * S1 * SH + w1) / Vv),
    'sin sinh': (sp.sin(w1 * s) * sp.sinh(w2 * s), (-w1 * C1 * SH + w2 * S1 * CH) / Vv),
}
for nm_, (f_, F_) in cl.items():
    Fs = F_.subs({S1: sp.sin(w1 * ell), C1: sp.cos(w1 * ell), SH: sp.sinh(w2 * ell), CH: sp.cosh(w2 * ell)})
    num_ = sp.integrate(f_, (s, 0, ell))
    chk = [float((num_ - Fs).subs({w1: a_v, w2: b_v, ell: l_v})) for (a_v, b_v, l_v) in ((0.9, 0.35, 0.4), (1.7, 0.2, 1.3))]
    print('  ∫0^ℓ %-9s ds = %s   [sympy の積分との差 %.1e, %.1e]' % (nm_, F_, chk[0], chk[1]))

print('\n== (d) tsolqu の ΔH の運動量の回転 (pyi の更新順序) ==')
a22, a24 = sp.symbols('a22 a24')
pxn_ = a22 * px + a24 * py
pyn_ = -a24 * pxn_ + a22 * py                  # K64: 更新後の pxi を使う
Mmap = sp.Matrix([[sp.diff(pxn_, px), sp.diff(pxn_, py)], [sp.diff(pyn_, px), sp.diff(pyn_, py)]])
print('  K64 の写像 (pxi, py) → (pxi'', pyi'') の行列 = %s' % Mmap.applyfunc(sp.expand))
print('  行列式 = %s → a22² + a24² = 1 のとき det = a22² = cos²(wℓ) (シンプレクティックなら 1)' % sp.expand(Mmap.det()))
pyo_ = -a24 * px + a22 * py                    # 正しい回転
Mok = sp.Matrix([[a22, a24], [-a24, a22]])
print('  正しい回転: 行列式 = %s' % sp.expand(Mok.det()))

# ---------------------------------------------------------------------------------------------------------
print('\n== (e) SOL + QUAD + SOL の追跡の 6×6 ヤコビアン |J^T S J - S| (pysad v2 既定) ==')
from pysad.ffscmd import run_sad_string
from pysad import compat

print('  compat mode = %s' % getattr(compat, '_mode', '?'))
brho = 3e9 / 299792458.0
Kq, Lq = 0.36, 0.4


def build(bz_inv_m):
    src = """MOMENTUM = 3GEV;
SOL S1=(BOUND=1 BZ=%.15g GEO=1 DISFRIN=1);
QUAD Q=(L=%.15g K1=%.15g DISFRIN=1);
SOL S2=(BOUND=1 BZ=0 DISFRIN=1);
MARK ST=();
LINE ONE=(ST S1 Q S2);
FFS;
USE ONE;
""" % (bz_inv_m * brho, Lq, Kq * Lq)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


Sj = np.zeros((6, 6))
for i in range(3):
    Sj[2 * i, 2 * i + 1] = 1
    Sj[2 * i + 1, 2 * i] = -1


def jac(ffs, c0, h=1e-6):
    J = np.zeros((6, 6))
    for k in range(6):
        cp, cm = list(c0), list(c0)
        cp[k] += h
        cm[k] -= h
        J[:, k] = (np.array(C.track_one(ffs, cp)) - np.array(C.track_one(ffs, cm))) / (2 * h)
    return J


rows = []
for bz_ in (1.3, 5.0, 10.0):
    ffs = build(bz_)
    for (pxv, pyv) in ((0.003, 0.0), (0.05, 0.0), (0.1, 0.0), (0.1, 0.05)):
        c0 = [0.001, pxv, -0.0005, pyv, 0.0, 0.0]
        J = jac(ffs, c0)
        err = np.abs(J.T @ Sj @ J - Sj).max()
        print('  b_z = %5.1f /m、(p_x, p_y) = (%.3f, %.2f): |J^T S J - S| = %.2e、det J - 1 = %+.2e' % (bz_, pxv, pyv, err, np.linalg.det(J) - 1))
