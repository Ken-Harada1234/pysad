# -*- coding: utf-8 -*-
"""取り付け誤差と COORD の章 (11_misalignment_and_coord.md) の補完。

(a) COORD の DIR != 0 の K64 の行列 = R1 R2m R3 (ヘルプは R3 R2m R1)。逆回転 (R3 R2 R1)^-1 との関係 (sympy)。
(b) BEND の tbrot の回転行列 (Fortran の rr) = 軸 t のまわりの回転 R_t(Θ) と軸 n のまわりの回転 R_n(-χ2) の積 (sympy、ロドリゲスの式)。
    入口の行列と出口 (= 転置) の関係、(11.10) の位置の式が「支点まわりの回転」であること。
(c) tbgrote の 6x6 行列 (Fortran の trans の各成分) と、写像の数値微分 (中心差分) の比較。
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

# ------------------------------------------------------------------------------------------------------ (a)
print('== (a) COORD の DIR != 0: K64 の行列の積の順序 (sympy) ==')
c1, s1, c2, s2, c3, s3 = sp.symbols('c1 s1 c2 s2 c3 s3')
R3 = sp.Matrix([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
R2p = sp.Matrix([[1, 0, 0], [0, c2, s2], [0, -s2, c2]])
R2m = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
R1 = sp.Matrix([[c1, 0, s1], [0, 1, 0], [-s1, 0, c1]])
M1_code = sp.Matrix([[c1*c3 + s1*s2*s3, s1*s2*c3 - c1*s3, s1*c2],
                     [c2*s3, c2*c3, -s2],
                     [-s1*c3 + c1*s2*s3, c1*s2*c3 + s1*s3, c1*c2]])
z33 = sp.zeros(3, 3)
ok = lambda M: M.applyfunc(sp.expand) == z33
print('  K64 の DIR != 0 の行列 - R1 R2m R3  = 0 か: %s' % ok(M1_code - R1 * R2m * R3))
print('  K64 の DIR != 0 の行列 - R3 R2m R1  = 0 か (ヘルプの式): %s' % ok(M1_code - R3 * R2m * R1))
print('  (R2m は R2 の s2 の符号を反転したもの)')
M0 = R3 * R2p * R1
print('  DIR = 0 の行列 M0 = R3 R2 R1 の逆 (転置) = R1^T R2^T R3^T = R1(-χ1) R2(-χ2) R3(-χ3)。符号を替えた記号 s_i → -s_i で書く:')
sub_neg = {s1: -s1, s2: -s2, s3: -s3}
Minv = (M0.T).applyfunc(sp.expand)
print('    M0^T - R1 R2m R3 (s1 → -s1, s3 → -s3 を入れる前) = 0 か: %s' % ok(Minv - R1.subs(s1, -s1) * R2m * R3.subs(s3, -s3)))
print('  → K64 の DIR != 0 は、真の逆 M0^T = R1(-χ1) R2m R3(-χ3) の、χ1 と χ3 の符号を反転して入れたもの')
print('     (R1(-χ1) → R1(χ1)、R3(-χ3) → R3(χ3))。逆変換にならない (ヘルプの "NOT the inverse") ことの代数的な理由。')

# ------------------------------------------------------------------------------------------------------ (b)
print('\n== (b) tbrot の回転行列 = R_t(Θ) R_n(-χ2) (ロドリゲスの式) ==')
phi, th, ch = sp.symbols('phi Theta chi2')
c0, s0 = sp.cos(phi), sp.sin(phi)
t_ax = sp.Matrix([-s0, 0, c0])
n_ax = sp.Matrix([c0, 0, s0])


def rodrigues(a, ang):
    K = sp.Matrix([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return sp.eye(3) + sp.sin(ang) * K + (1 - sp.cos(ang)) * K * K


# Fortran の rr (ent)
st, ct, s2_, c2_ = sp.sin(th), sp.cos(th), sp.sin(ch), sp.cos(ch)
dct, dc2 = ct - 1, c2_ - 1
rr = sp.zeros(3, 3)
rr[0, 0] = 1 + dct * c0**2 + dc2 * s0**2 + c0 * s0 * s2_ * st
rr[1, 1] = 1 + c2_ * dct + dc2
rr[2, 2] = 1 + dc2 * c0**2 + dct * s0**2 - c0 * s0 * s2_ * st
rr[0, 1] = s0 * s2_ - c0 * c2_ * st
rr[0, 2] = c0 * ((dct - dc2) * s0 - c0 * s2_ * st)
rr[1, 0] = -ct * s0 * s2_ + c0 * st
rr[1, 2] = c0 * ct * s2_ + s0 * st
rr[2, 0] = s0 * (c0 * (dct - dc2) + s0 * s2_ * st)
rr[2, 1] = -c0 * s2_ - c2_ * s0 * st
cands = {
    'R_t(+Θ) R_n(-χ2)': rodrigues(t_ax, th) * rodrigues(n_ax, -ch),
    'R_n(-χ2) R_t(+Θ)': rodrigues(n_ax, -ch) * rodrigues(t_ax, th),
    'R_t(+Θ) R_n(+χ2)': rodrigues(t_ax, th) * rodrigues(n_ax, ch),
    'R_t(-Θ) R_n(-χ2)': rodrigues(t_ax, -th) * rodrigues(n_ax, -ch),
}
vals = {phi: 0.37, th: 0.21, ch: -0.13}
for nm, M in cands.items():
    d = np.abs(np.array((M - rr).subs(vals).evalf(), dtype=float)).max()
    print('  Fortran の rr (入口) と %s の差の最大 (φ = 0.37、Θ = 0.21、χ2 = -0.13): %.1e' % (nm, d))
M = cands['R_t(+Θ) R_n(-χ2)']
dd = (M - rr).applyfunc(lambda e: sp.simplify(sp.expand_trig(e)))
print('  R_t(+Θ) R_n(-χ2) - rr (入口) を sympy で簡約: %s' % ('すべて 0' if dd == z33 else dd))
print('  出口の rr は入口の転置 (Fortran の ent = .false. の代入: rr(2,1) と rr(1,2) などを入れ替え) か (数値):')
rr_ext = rr.copy()
for (i, j) in ((0, 1), (0, 2), (1, 2)):
    rr_ext[i, j], rr_ext[j, i] = rr[j, i], rr[i, j]
print('     転置との差 = %.1e' % np.abs(np.array((rr_ext - rr.T).subs(vals).evalf(), dtype=float)).max())
# Θ のみのとき R11 = 1 + (cosΘ - 1) c0²
rr_t = rr.subs(ch, 0)
print('  Θ のみ (χ2 = 0): R11 = %s、 R_t(Θ) の x 成分の 2 乗 = %s (R11 = 1 - (1 - cosΘ)(1 - t_x²)... ではなく 1 + (cosΘ - 1) c0² で、t_x² = s0²、1 - t_x² = c0²)' % (
    sp.simplify(rr_t[0, 0]), sp.simplify(t_ax[0]**2)))

# (11.10): 位置は r' = r_p + R (r - r_p)
print('\n  (11.10) の位置の式: x_f = x + dr11 (x + xl0) + r12 y + r13 zi = x + [(R - 1)(r - r_p)]_x (r_p = dl (-s0, 0, c0)、zi = -dl c0、xl0 = dl s0)')
dl_, x_, y_ = sp.symbols('dl x y')
r_vec = sp.Matrix([x_, y_, 0])
r_p = dl_ * t_ax
rp_new = (r_p + M * (r_vec - r_p))
xl0 = dl_ * s0
zi = -dl_ * c0
xf = x_ + (rr[0, 0] - 1) * (x_ + xl0) + rr[0, 1] * y_ + rr[0, 2] * zi
yf = y_ + rr[1, 0] * (x_ + xl0) + (rr[1, 1] - 1) * y_ + rr[1, 2] * zi
zf = rr[2, 0] * (x_ + xl0) + rr[2, 1] * y_ + (rr[2, 2] - 1) * zi
vv = {phi: 0.37, th: 0.21, ch: -0.13, dl_: 0.6, x_: 0.003, y_: -0.002}
print('     (xf, yf, zf) - (r_p + R(r - r_p)) の差の最大 = %.1e' % max(abs(float((a - b).subs(vv))) for a, b in zip((xf, yf, zf), rp_new)))
print('     (zf は面 z = 0 から測った位置の z 成分。r_p の z 成分は dl c0 で、zf は元の面 z = 0 に対して。)')

# ------------------------------------------------------------------------------------------------------ (c)
print('\n== (c) tbgrote の 6×6 行列と数値微分 ==')


def pxy2dpz(px, py):
    return math.sqrt(1.0 - px * px - py * py) - 1.0


def tbgrote_map(cod, dr, dx, dz):
    xi = cod[0] + dx
    yi = cod[2]
    pxi, pyi = cod[1], cod[3]
    pr = 1.0 + cod[5]
    dpzi = pxy2dpz(pxi / pr, pyi / pr) * pr
    pzi = pr + dpzi
    zi = dz
    xf = dr[0, 0] * xi + dr[0, 1] * yi + dr[0, 2] * zi + cod[0]
    yf = dr[1, 0] * xi + dr[1, 1] * yi + dr[1, 2] * zi + yi
    zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * zi
    r11, r22, r33 = 1 + dr[0, 0], 1 + dr[1, 1], 1 + dr[2, 2]
    pxf = r11 * pxi + dr[0, 1] * pyi + dr[0, 2] * pzi
    pyf = dr[1, 0] * pxi + r22 * pyi + dr[1, 2] * pzi
    dpzf = dr[2, 0] * pxi + dr[2, 1] * pyi + r33 * dpzi + dr[2, 2] * pr
    pzf = pr + dpzf
    out = np.array(cod, float)
    out[1] = pxf
    out[3] = pyf
    out[0] = xf - pxf * zf / pzf
    out[2] = yf - pyf * zf / pzf
    out[4] = cod[4] + pr * zf / pzf
    return out


def tbgrote_trans(cod, dr, dx, dz):
    """Fortran の trans (6x6) の成分を、そのまま写す"""
    xi = cod[0] + dx
    yi = cod[2]
    pxi, pyi = cod[1], cod[3]
    pr = 1.0 + cod[5]
    dpzi = pxy2dpz(pxi / pr, pyi / pr) * pr
    pzi = pr + dpzi
    zi = dz
    zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * zi
    r11, r22, r33 = 1 + dr[0, 0], 1 + dr[1, 1], 1 + dr[2, 2]
    pxf = r11 * pxi + dr[0, 1] * pyi + dr[0, 2] * pzi
    pyf = dr[1, 0] * pxi + r22 * pyi + dr[1, 2] * pzi
    dpzf = dr[2, 0] * pxi + dr[2, 1] * pyi + r33 * dpzi + dr[2, 2] * pr
    pzf = pr + dpzf
    pxfpzf, pyfpzf, zfpzf, prpzf = pxf / pzf, pyf / pzf, zf / pzf, pr / pzf
    dpzidpx, dpzidpy, dpzidp = -pxi / pzi, -pyi / pzi, pr / pzi
    T = np.zeros((6, 6))
    T[1, 1] = r11 + dr[0, 2] * dpzidpx
    T[1, 3] = dr[0, 1] + dr[0, 2] * dpzidpy
    T[1, 5] = dr[0, 2] * dpzidp
    T[3, 1] = dr[1, 0] + dr[1, 2] * dpzidpx
    T[3, 3] = r22 + dr[1, 2] * dpzidpy
    T[3, 5] = dr[1, 2] * dpzidp
    T[0, 0] = r11 - dr[2, 0] * pxfpzf
    T[0, 1] = -zfpzf * (T[1, 1] - pxfpzf * (dr[2, 0] + r33 * dpzidpx))
    T[0, 2] = dr[0, 1] - dr[2, 1] * pxfpzf
    T[0, 3] = -zfpzf * (T[1, 3] - pxfpzf * (dr[2, 1] + r33 * dpzidpy))
    T[0, 5] = -zfpzf * (dr[0, 2] - r33 * pxfpzf) * dpzidp
    T[2, 0] = dr[1, 0] - dr[2, 0] * pyfpzf
    T[2, 1] = -zfpzf * (T[3, 1] - pyfpzf * (dr[2, 0] + r33 * dpzidpx))
    T[2, 2] = r22 - dr[2, 1] * pyfpzf
    T[2, 3] = -zfpzf * (T[3, 3] - pyfpzf * (dr[2, 1] + r33 * dpzidpy))
    T[2, 5] = -zfpzf * (dr[1, 2] - r33 * pyfpzf) * dpzidp
    T[4, 0] = dr[2, 0] * prpzf
    T[4, 1] = -zfpzf * prpzf * (dr[2, 0] + r33 * dpzidpx)
    T[4, 2] = dr[2, 1] * prpzf
    T[4, 3] = -zfpzf * prpzf * (dr[2, 1] + r33 * dpzidpy)
    T[4, 4] = 1.0
    T[4, 5] = r33 * (dpzf * prpzf + dpzi / pzi) - dr[2, 2]
    T[5, 5] = 1.0
    return T


def make_dr(phig, dtheta, dchi2):
    c0_, s0_ = math.cos(phig), math.sin(phig)
    s2v, c2v = math.sin(dchi2), math.cos(dchi2)
    stv, ctv = math.sin(dtheta), math.cos(dtheta)
    dctv, dc2v = ctv - 1, c2v - 1
    d = np.zeros((3, 3))
    d[0, 0] = dctv * c0_**2 + dc2v * s0_**2 + c0_ * s0_ * s2v * stv
    d[1, 1] = c2v * dctv + dc2v
    d[2, 2] = dc2v * c0_**2 + dctv * s0_**2 - c0_ * s0_ * s2v * stv
    d[0, 1] = s0_ * s2v - c0_ * c2v * stv
    d[0, 2] = c0_ * ((dctv - dc2v) * s0_ - c0_ * s2v * stv)
    d[1, 0] = -ctv * s0_ * s2v + c0_ * stv
    d[1, 2] = c0_ * ctv * s2v + s0_ * stv
    d[2, 0] = s0_ * (c0_ * (dctv - dc2v) + s0_ * s2v * stv)
    d[2, 1] = -c0_ * s2v - c2v * s0_ * stv
    return d


phig, dth, dch = 0.0654, 0.02, -0.015
alg = 0.6
dl = alg * math.sin(phig) / phig
dr = make_dr(phig, dth, dch)
dxp, dzp = dl * math.sin(phig), -dl * math.cos(phig)
S = np.zeros((6, 6))
for i in (0, 2, 4):
    S[i, i + 1] = 1
    S[i + 1, i] = -1
rng = np.random.default_rng(3)
worst = 0.0
worst_s = 0.0
t56 = []
for k in range(5):
    cod = np.array([rng.normal(0, 2e-3), rng.normal(0, 3e-3), rng.normal(0, 2e-3), rng.normal(0, 3e-3), rng.normal(0, 1e-2), rng.normal(0, 2e-3)])
    T = tbgrote_trans(cod, dr, dxp, dzp)
    Jn = np.zeros((6, 6))
    for j in range(6):
        h = 1e-6
        d = np.zeros(6)
        d[j] = h
        Jn[:, j] = (tbgrote_map(cod + d, dr, dxp, dzp) - tbgrote_map(cod - d, dr, dxp, dzp)) / (2 * h)
    mask = np.ones((6, 6), bool)
    mask[4, 5] = False                      # (5,6) は別に (B109)
    worst = max(worst, np.abs((T - Jn)[mask]).max())
    # (5,6) の 3 通り: Fortran の現在の式、コメントアウトされた式、数値微分
    xi = cod[0] + dxp
    yi = cod[2]
    pr = 1 + cod[5]
    pxi, pyi = cod[1], cod[3]
    dpzi = pxy2dpz(pxi / pr, pyi / pr) * pr
    pzi = pr + dpzi
    zf = dr[2, 0] * xi + dr[2, 1] * yi + dr[2, 2] * dzp
    r33 = 1 + dr[2, 2]
    dpzf = dr[2, 0] * pxi + dr[2, 1] * pyi + r33 * dpzi + dr[2, 2] * pr
    pzf = pr + dpzf
    old = zf / pzf * (1 - pr / pzf * r33 * pr / pzi)
    t56.append((T[4, 5], old, Jn[4, 5], pxi**2 / pr**2 + pyi**2 / pr**2))
    worst_s = max(worst_s, np.abs(T.T @ S @ T - S).max())
print('  5 つの乱数点 (x, y ~ 2 mm、p_x, p_y ~ 3 mrad、z ~ 1 cm、δ ~ 2e-3、φ = %.4f、Θ = %.3f、χ2 = %.3f)で、' % (phig, dth, dch))
print('  Fortran の trans の成分 ((5,6) 以外の 35 成分) と、写像の中心差分 (刻み 1e-6) の最大の差 = %.1e' % worst)
print('  (5,6) 成分 [∂z/∂δ]: Fortran の現在の式 r33(dpzf pr/pzf + dpzi/pzi) - dr33、コメントアウトされた式 zf/pzf (1 - pr r33 dpzidp/pzf)、数値微分、q² = (px² + py²)/pr²')
for a_, b_, c_, d_ in t56:
    print('     %.3e   %.3e   %.3e   (q² = %.2e)' % (a_, b_, c_, d_))
print('  → 現在の式は数値微分と一致せず (dr = 0 でも dpzi (pr + 1)/pzi が残り、大きさは q² の程度。B109)、コメントアウトされた式は数値微分と一致する。')

# ------------------------------------------------------------------------------------------------------ (d)
print('\n== (d) (11.3)(11.4) の 1 次の展開 (§3.5): 原点での 1 次の係数 (sympy) ==')
X, PX, Y, PY, Z, DL = sp.symbols('x p_x y p_y z delta')
K1, K2, K3, DX, DY, DZ = sp.symbols('chi1 chi2 chi3 DX DY DZ')
pp = 1 + DL
pzz = sp.sqrt(pp**2 - PX**2 - PY**2)
uu = sp.Matrix([PX, PY, pzz]) / pp
rr0 = sp.Matrix([X, Y, 0]) + Z * uu - sp.Matrix([DX, DY, DZ])
Rm = (sp.Matrix([[sp.cos(K3), -sp.sin(K3), 0], [sp.sin(K3), sp.cos(K3), 0], [0, 0, 1]])
      * sp.Matrix([[1, 0, 0], [0, sp.cos(K2), sp.sin(K2)], [0, -sp.sin(K2), sp.cos(K2)]])
      * sp.Matrix([[sp.cos(K1), 0, sp.sin(K1)], [0, 1, 0], [-sp.sin(K1), 0, sp.cos(K1)]]))
r1_ = Rm * rr0
P1_ = Rm * sp.Matrix([PX, PY, pzz])
z1_ = pp * r1_[2] / P1_[2]
x1_ = r1_[0] - (P1_[0] / pp) * z1_
y1_ = r1_[1] - (P1_[1] / pp) * z1_
out = [x1_, P1_[0], y1_, P1_[1], z1_, sp.sqrt(P1_[0]**2 + P1_[1]**2 + P1_[2]**2) - 1]
names = ["x'", "p_x'", "y'", "p_y'", "z'", "delta'"]
pz_s = sp.symbols('p_z', positive=True)           # p_z = sqrt(p² - p_x² - p_y²)
pars = [K1, K2, K3, DX, DY, DZ]
zero_pars = {v: 0 for v in pars}
print('  誤差のパラメータでの 1 次の係数 (∂出力/∂パラメータ、パラメータ = 0。x, p_x, y, p_y, z, δ は一般の値。p_z = sqrt(p² - p_x² - p_y²)、p = 1 + δ):')
for nm, f in zip(names, out):
    terms = []
    for v in pars:
        c = sp.simplify(sp.diff(f, v).subs(zero_pars))
        c = sp.simplify(c.subs(sp.sqrt(-PX**2 - PY**2 + (DL + 1)**2), pz_s))
        if c != 0:
            terms.append('(%s)·%s' % (c, v))
    print("   %-7s = (誤差なしの値) + %s" % (nm, ' + '.join(terms) if terms else '0'))
print('  (原点 x = p_x = ... = 0 の線形項は x_new = x - DX、p_x_new = p_x + chi1 など。)')
