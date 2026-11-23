# -*- coding: utf-8 -*-
r"""BEAMBEAM の交差角の写像 (16_undulator_beambeam.md §8、(16.12)) の導出の確認。

構成: (1) 4 元運動量 (E = pn、横の運動量 px, py、縦 pz) に、座標軸の回転 (角 φ) とその回転後の x 方向の Lorentz boost (速度 β = sin φ) を施し、
運動量を新しい基準のエネルギー E2ref = cos φ で割る。(2) 位置: 事象 (x, y, s = 0, ct = -z) に同じ変換を施し、新しい運動量の向きに直進させて s2 = 0 の面に戻す。
(3) 結果が bbstrhl.f の入口の式 (16.12) と一致することを、記号計算 (sympy) と 40 桁の数値で確かめる。(4) 出口の式が入口の逆写像であること。
(5) 入口の写像の 6×6 ヤコビアン J が J^T S J = S/cos φ (共形シンプレクティック。因子 1/cos φ は運動量の単位の変更 p0 -> p0 cos φ) であること。
(6) 小振幅の極限。
"""
import random
import sympy as sp
import mpmath as mp
import _common as C

C.start(__file__)

s, c = sp.symbols('s c', positive=True)       # sin φ, cos φ (s^2 + c^2 = 1 は最後に使う)
t = s / c
x, y, z, px, py, pz, pn = sp.symbols('x y z p_x p_y p_z p_n')
print('== (1) 運動量: 回転 + boost + 規格化 (pn、px、py、pz は独立な記号。pn^2 = px^2+py^2+pz^2 は使わない) ==')
px1 = c * px + s * pz                  # 座標軸を φ だけ回した系の x 運動量
ps1 = c * pz - s * px                  # 同、縦運動量
E2 = (pn - s * px1) / c                # boost (γ = 1/c、γβ = tan φ): E' = γ (E - β px1)
px2 = (px1 - s * pn) / c               # px' = γ (px1 - β E)
h1 = pn - pz


def red(e):
    return sp.simplify(sp.expand(e).subs(c**2, 1 - s**2))


print('E2 - ps1 - h1/c =', red(E2 - ps1 - h1 / c), ' (E2 - ps1 = h1/cos φ)')
dprime = E2 / c - 1                   # 新しい δ = E2/E2ref - 1、E2ref = cos φ
fort_delta = (pn - 1) - t * px + t**2 * h1
print('δ_new (boost) - δ_new (Fortran) =', red(dprime - fort_delta))
px_new = px2 / c
print('px_new (boost) - (px - tan φ h1)/cos φ =', red(px_new - (px - t * h1) / c))
print('py_new = py/cos φ (boost は y 方向の運動量を変えない: py2 = py)。')

print('\n== (2) 位置: 事象 (x, y, s=0, ct=-z) -> 回転 -> boost -> 直進で s2 = 0 の面へ ==')
ct = -z
x1 = c * x                               # 回転: x1 = c x + s S (S = 0)
x2 = (x1 - s * ct) / c                   # boost: γ (x1 - β ct)
ct2 = (ct - s * x1) / c                  # boost: γ (ct - β x1)
ps2 = ps1
# s2 = s1 = -s x を 0 にするため、新しい運動量の向き (px2, py, ps2; エネルギー E2) に Δs = s x だけ進む
ds_ = s * x
xn = x2 + ds_ * px2 / ps2
yn = y + ds_ * py / ps2
ctn = ct2 + ds_ * E2 / ps2
zn = -ctn
psn = ps2 / c                            # 新しい規格化での縦運動量。Fortran の 1/pz
x_F = t * z + (1 + s * px_new / psn) * x
y_F = y + s * (py / c) / psn * x
z_F = z / c - t / c * h1 * (1 / psn) * x     # Fortran: z/cosx - tanx/cosx * h1 * xx * pz (h1 は変換前の h1)
h1_new = dprime + 1 - psn
print("h1_new = (1+δ') - ps_new - h1/cos^2 φ =", red(h1_new - h1 / c**2))
print('x_new (導出) - x_F =', red(xn - x_F))
print('y_new (導出) - y_F =', red(yn - y_F))
print('z_new (導出) - z_F =', red(zn - z_F))

print('\n== (3) 数値 (40 桁): Fortran の式を逐語的に写した関数と、導出した写像 (回転 + boost + 直進) を、質量殻上の点で比べる ==')
mp.mp.dps = 40


def fortran_entry(v, phi):
    xx_, px_, y_, py_, z_, g_ = v
    sinx, cosx, tanx = mp.sin(phi), mp.cos(phi), mp.tan(phi)
    pxy2 = px_**2 + py_**2
    pn_ = 1 + g_
    pz_ = mp.sqrt(pn_**2 - pxy2)
    h1_ = pxy2 / (pn_ + pz_)
    g_ = g_ - tanx * px_ + tanx * tanx * h1_
    px_ = (px_ - tanx * h1_) / cosx
    py_ = py_ / cosx
    pxy2 = px_**2 + py_**2
    pz_ = 1 / mp.sqrt((1 + g_)**2 - pxy2)
    xo = xx_
    x_n = tanx * z_ + (1 + sinx * px_ * pz_) * xo
    y_n = y_ + sinx * py_ * xo * pz_
    z_n = z_ / cosx - tanx / cosx * h1_ * xo * pz_
    return [x_n, px_, y_n, py_, z_n, g_]


def fortran_exit(v, phi):
    xx_, px_, y_, py_, z_, g_ = v
    sinx, cosx, tanx = mp.sin(phi), mp.cos(phi), mp.tan(phi)
    pxy2 = px_**2 + py_**2
    pn_ = 1 + g_
    pz_ = mp.sqrt(pn_**2 - pxy2)
    h1_ = pxy2 / (pn_ + pz_)
    pz_ = 1 / pz_
    x_n = (xx_ - sinx * z_) / (1 + (sinx * (px_ + sinx * h1_)) * pz_)
    y_n = y_ - sinx * py_ * pz_ * x_n
    z_n = z_ * cosx + sinx * cosx * h1_ * pz_ * x_n
    h1_ = h1_ * cosx**2
    xo = px_
    px_n = xo * cosx + tanx * h1_
    py_n = py_ * cosx
    g_n = g_ + sinx * xo
    return [x_n, px_n, y_n, py_n, z_n, g_n]


def derived_entry(v, phi):
    sv, cv = mp.sin(phi), mp.cos(phi)
    xx_, px_, y_, py_, z_, g_ = v
    pn_ = 1 + g_
    pz_ = mp.sqrt(pn_**2 - px_**2 - py_**2)
    p1 = cv * px_ + sv * pz_
    ps = cv * pz_ - sv * px_
    E = (pn_ - sv * p1) / cv
    p2 = (p1 - sv * pn_) / cv
    ctv = -z_
    X1 = cv * xx_
    X2 = (X1 - sv * ctv) / cv
    CT2 = (ctv - sv * X1) / cv
    dsv = sv * xx_
    xn_ = X2 + dsv * p2 / ps
    yn_ = y_ + dsv * py_ / ps
    ctn_ = CT2 + dsv * E / ps
    return [xn_, p2 / cv, yn_, py_ / cv, -ctn_, E / cv - 1]


random.seed(5)
phi = mp.mpf('0.0137')
maxd = 0
maxr = 0
for k in range(5):
    v = [mp.mpf(random.uniform(-1, 1)) * 1e-3, mp.mpf(random.uniform(-1, 1)) * 3e-3, mp.mpf(random.uniform(-1, 1)) * 1e-3,
         mp.mpf(random.uniform(-1, 1)) * 3e-3, mp.mpf(random.uniform(-1, 1)) * 1e-2, mp.mpf(random.uniform(-1, 1)) * 1e-2]
    a = fortran_entry(v, phi)
    b = derived_entry(v, phi)
    r = fortran_exit(a, phi)
    maxd = max(maxd, max(abs(a[i] - b[i]) for i in range(6)))
    maxr = max(maxr, max(abs(r[i] - v[i]) for i in range(6)))
print('入口: |Fortran の式 - 導出した写像| の最大 (5 点) = %s' % mp.nstr(maxd, 3))
print('出口 ∘ 入口 = 恒等: |出口(入口(v)) - v| の最大 = %s' % mp.nstr(maxr, 3))

print('\n== (4) 入口の写像のヤコビアン: J^T S J = S/cos φ (共形シンプレクティック) ==')
X = sp.symbols('X PX Y PY Z G')
phi_s = sp.symbols('phi')
sinx, cosx, tanx = sp.sin(phi_s), sp.cos(phi_s), sp.tan(phi_s)
pxy2 = X[1]**2 + X[3]**2
pn_e = 1 + X[5]
pz_e = sp.sqrt(pn_e**2 - pxy2)
h1_e = pxy2 / (pn_e + pz_e)
g_n = X[5] - tanx * X[1] + tanx**2 * h1_e
px_n = (X[1] - tanx * h1_e) / cosx
py_n = X[3] / cosx
pz_n = 1 / sp.sqrt((1 + g_n)**2 - px_n**2 - py_n**2)
x_n = tanx * X[4] + (1 + sinx * px_n * pz_n) * X[0]
y_n = X[2] + sinx * py_n * X[0] * pz_n
z_n = X[4] / cosx - tanx / cosx * h1_e * X[0] * pz_n
Fm = sp.Matrix([x_n, px_n, y_n, py_n, z_n, g_n])
Jm = Fm.jacobian(sp.Matrix(X))
pt = {X[0]: sp.Rational(3, 10000), X[1]: sp.Rational(2, 1000), X[2]: -sp.Rational(1, 10000), X[3]: sp.Rational(-15, 10000),
      X[4]: sp.Rational(7, 1000), X[5]: sp.Rational(5, 1000), phi_s: sp.Rational(137, 10000)}
Jn = Jm.subs(pt).evalf(40)
Sm = sp.zeros(6, 6)
for k in range(3):
    Sm[2 * k, 2 * k + 1] = 1
    Sm[2 * k + 1, 2 * k] = -1
cphi = sp.cos(pt[phi_s]).evalf(40)
dev = Jn.T * Sm * Jn - Sm / cphi
print('max |J^T S J - S/cos φ| = %s (φ = 0.0137、1/cos φ - 1 = %s)' % (sp.N(max(abs(e) for e in dev), 3), sp.N(1 / cphi - 1, 3)))
print('(参考) max |J^T S J - S| = %s' % sp.N(max(abs(e) for e in (Jn.T * Sm * Jn - Sm)), 3))

print('\n== (5) 小振幅の極限 (原点まわりの 1 次) ==')
J0 = Jm.subs({X[0]: 0, X[1]: 0, X[2]: 0, X[3]: 0, X[4]: 0, X[5]: 0})
print('J(0) (行・列の順は x, px, y, py, z, δ) =')
print(sp.simplify(J0))
v = [mp.mpf(0), mp.mpf('1e-3'), mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(0)]
xx1 = fortran_entry(v, phi)
print('px = 1e-3、他は 0: px_new = %s、px/cos φ = %s、差 = %s、-tan φ px^2/(2 cos φ) = %s' % (
    mp.nstr(xx1[1], 15), mp.nstr(mp.mpf('1e-3') / mp.cos(phi), 15), mp.nstr(xx1[1] - mp.mpf('1e-3') / mp.cos(phi), 6),
    mp.nstr(-mp.tan(phi) * mp.mpf('1e-6') / 2 / mp.cos(phi), 6)))
v = [mp.mpf(0)] * 5 + [mp.mpf('1e-3')]
print('δ = 1e-3、他は 0: px_new = %s (δ の 1 次の項はない。h1 = 0)' % mp.nstr(fortran_entry(v, phi)[1], 6))
print('px = 1e-3、δ = 0: δ_new = %s、-tan φ px = %s' % (mp.nstr(xx1[5], 10), mp.nstr(-mp.tan(phi) * mp.mpf('1e-3'), 10)))
