# -*- coding: utf-8 -*-
r"""UND の途中式の確認 (16_undulator_beambeam.md §3)。

(a) B = ∇×A の各成分 (F0 の項、G0 の項) を、条件 Ky² = Kz² + Kx²、Qx² = Kz² + Qy² を使わずに求め、∇·B ≡ 0 と、∇×B = 0 の成分が条件のときだけ 0 になることを示す。
(b) 周期の積分と、端の A。
(c) H の展開 (傍軸)、正準方程式、z の更新。
(d) スライスの写像: 陰的な式の 2×2 の線形系 (16.5) の解と、不動点反復の解が一致すること。写像の 6×6 ヤコビアンが任意の ds でシンプレクティックであること。
(e) 磁場なしの UND の z の更新の欠落: R56 = L/γ² との比較。
"""
import math
import sympy as sp
import _common as C

C.start(__file__)
x, y, s = sp.symbols('x y s', real=True)
F0, G0, Kz, Kx, Ky, Qx, Qy, ph = sp.symbols('F0 G0 Kz Kx Ky Qx Qy phi', positive=True)


def curl(V):
    return sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], s), sp.diff(V[0], s) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])


# ---------------------------------------------------------------- (a)
print('== (a) 磁場の成分 (条件を使わない) ==')
A_F = sp.Matrix([F0 / Kz * sp.cos(Kx * x) * sp.cosh(Ky * y) * sp.sin(Kz * s),
                 Kx * F0 / (Ky * Kz) * sp.sinh(Ky * y) * sp.sin(Kx * x) * sp.sin(Kz * s), 0])
# B = ∇×A (A_s = 0): Bx = -∂s Ay、By = ∂s Ax、Bz = ∂x Ay - ∂y Ax
B_F = sp.Matrix([-sp.diff(A_F[1], s), sp.diff(A_F[0], s), sp.diff(A_F[1], x) - sp.diff(A_F[0], y)])
for nm, e in zip('xyz', B_F):
    print('F0 の項 B_%s =' % nm, sp.simplify(e))
print('∇·B (F0 の項、条件なし) =', sp.simplify(sum(sp.diff(B_F[i], v) for i, v in enumerate((x, y, s)))))
cB = curl(B_F)
for nm, e in zip('xyz', cB):
    print('(∇×B)_%s (F0 の項、条件なし) =' % nm, sp.simplify(e))
cond_F = Ky**2 - Kx**2 - Kz**2
print('(∇×B)_x に Ky² = Kz² + Kx² を代入 =', sp.simplify(cB[0].subs(Ky, sp.sqrt(Kz**2 + Kx**2))))

A_G = sp.Matrix([-Qy * G0 / (Qx * Kz) * sp.sinh(Qx * x) * sp.sin(Qy * y) * sp.sin(Kz * s + ph),
                 -G0 / Kz * sp.cos(Qy * y) * sp.cosh(Qx * x) * sp.sin(Kz * s + ph), 0])
B_G = sp.Matrix([-sp.diff(A_G[1], s), sp.diff(A_G[0], s), sp.diff(A_G[1], x) - sp.diff(A_G[0], y)])
for nm, e in zip('xyz', B_G):
    print('G0 の項 B_%s =' % nm, sp.simplify(e))
print('∇·B (G0 の項、条件なし) =', sp.simplify(sum(sp.diff(B_G[i], v) for i, v in enumerate((x, y, s)))))
cG = curl(B_G)
for nm, e in zip('xyz', cG):
    print('(∇×B)_%s (G0 の項、条件なし) =' % nm, sp.simplify(e))
print('(∇×B)_x に Qx² = Kz² + Qy² を代入 =', sp.simplify(cG[0].subs(Qx, sp.sqrt(Kz**2 + Qy**2))))
print('(∇×B)_y に Qx² = Kz² + Qy² を代入 =', sp.simplify(cG[1].subs(Qx, sp.sqrt(Kz**2 + Qy**2))))
print('(∇×B)_z に Qx² = Kz² + Qy² を代入 =', sp.simplify(cG[2].subs(Qx, sp.sqrt(Kz**2 + Qy**2))))

# ---------------------------------------------------------------- (b)
print('\n== (b) 周期の積分と端の A ==')
lam = sp.symbols('lambda_u', positive=True)
print('∫_0^{2π N/Kz} cos(Kz s) ds (N 周期) =', sp.simplify(sp.integrate(sp.cos(Kz * s), (s, 0, 2 * sp.pi / Kz))))
print('軸上の A_x(s=0) = %s、A_y(s=0) = %s (φ ≠ 0 のとき -G0 sin φ/Kz)' % (A_F[0].subs({x: 0, y: 0, s: 0}), sp.simplify(A_G[1].subs({x: 0, y: 0, s: 0}))))

# ---------------------------------------------------------------- (c)
print('\n== (c) H の展開 ==')
u, pnv = sp.symbols('u p_n', positive=True)       # u = π² (運動学的な横運動量の 2 乗)
Hs = -sp.sqrt(pnv**2 - u)
print('H = -√(pn² - π²) の π² による展開:', sp.series(Hs, u, 0, 3))
px, py, Z, dl = sp.symbols('p_x p_y z delta')
Ax = sp.Function('A_x')(x, y, s)
Ay = sp.Function('A_y')(x, y, s)
H2 = ((px - Ax)**2 + (py - Ay)**2) / (2 * (1 + dl))
print("x' = ∂H/∂px =", sp.simplify(sp.diff(H2, px)))
print("px' = -∂H/∂x =", sp.simplify(-sp.diff(H2, x)))
print("z' = ∂H/∂δ =", sp.simplify(sp.diff(H2, dl)))

# ---------------------------------------------------------------- (d)
print('\n== (d) スライスの写像 ==')
mp = sp
xv, yv, pxv, pyv, zv, dv_ = sp.symbols('X Y PX PY Z DL')
Fn, Gn, KZn, KXn, KYn, QXn, QYn, PHn, dsn, sn = [sp.Rational(v) for v in ('3/10', '1/5', '31/10', '4', '0', '0', '5', '3/10', '1/50', '7/10')]
KYn = sp.sqrt(KZn**2 + KXn**2)
QXn = sp.sqrt(KZn**2 + QYn**2)
ff = (Fn / KZn * sp.cos(KXn * xv) * sp.cosh(KYn * yv) * sp.sin(KZn * sn)
      - QYn * Gn / QXn / KZn * sp.sinh(QXn * xv) * sp.sin(QYn * yv) * sp.sin(KZn * sn + PHn))
gg = (-Gn / KZn * sp.cos(QYn * yv) * sp.cosh(QXn * xv) * sp.sin(KZn * sn + PHn)
      + KXn * Fn / KYn / KZn * sp.sinh(KYn * yv) * sp.sin(KXn * xv) * sp.sin(KZn * sn))
pnn = 1 + dv_
dpds = dsn / pnn
fxh, fyh, gxh, gyh = [dpds * sp.diff(e, v) for e, v in ((ff, xv), (ff, yv), (gg, xv), (gg, yv))]
det1a = (1 - fxh) * (1 - gyh) - fyh * gxh
rx = pxv - (fxh * ff + gxh * gg)
ry = pyv - (fyh * ff + gyh * gg)
Pxn = ((1 - gyh) * rx + gxh * ry) / det1a          # 逆行列 (1/det) [[1-ĝy, ĝx],[f̂y, 1-f̂x]]
Pyn = (fyh * rx + (1 - fxh) * ry) / det1a
# 陰的な式 P = p - ds ∂q H(q, P) が満たされているか
pi_x = Pxn - ff
pi_y = Pyn - gg
res_x = Pxn - (pxv + dpds * (pi_x * sp.diff(ff, xv) + pi_y * sp.diff(gg, xv)))
res_y = Pyn - (pyv + dpds * (pi_x * sp.diff(ff, yv) + pi_y * sp.diff(gg, yv)))
pt = {xv: sp.Rational(3, 1000), yv: -sp.Rational(2, 1000), pxv: sp.Rational(1, 1000), pyv: sp.Rational(-5, 10000), dv_: sp.Rational(2, 1000), zv: 0}
print('陰的な式の残差 (P - p + ds ∂qH): x = %s、y = %s' % (sp.N(res_x.subs(pt), 5), sp.N(res_y.subs(pt), 5)))
xn = xv + dpds * pi_x
yn = yv + dpds * pi_y
zn = zv - dsn * (pi_x**2 + pi_y**2) / (2 * pnn**2)
Fmap = sp.Matrix([xn, Pxn, yn, Pyn, zn, dv_])
Jm = Fmap.jacobian(sp.Matrix([xv, pxv, yv, pyv, zv, dv_]))
Jn = Jm.subs(pt).evalf(40)
Sm = sp.zeros(6, 6)
for k in range(3):
    Sm[2 * k, 2 * k + 1] = 1
    Sm[2 * k + 1, 2 * k] = -1
print('ds = 0.02 (大きな刻み)、1 スライスの写像 (x, px, y, py, z, δ): max|J^T S J - S| = %s、det J = %s' % (
    sp.N(max(abs(e) for e in (Jn.T * Sm * Jn - Sm)), 3), sp.N(Jn.det(), 20)))

# ---------------------------------------------------------------- (e)
print('\n== (e) -dv L と R56 = L/γ² ==')
me = 0.51099895e6
for pgev in (3e9, 0.5e9):
    p0 = pgev / me
    h0 = math.sqrt(1 + p0**2)
    d = 1e-3
    pn = 1 + d
    h1 = math.sqrt(1 + (pn * p0)**2)
    dvv = -d * (2 + d) / h1 / (h1 + pn * h0)
    print('E = %.1f GeV: dv = -δ(2+δ)/[h1(h1+pn h0)] = %.4e、-δ/γ² = %.4e、比 %.5f' % (pgev / 1e9, dvv, -d / h0**2, dvv / (-d / h0**2)))
