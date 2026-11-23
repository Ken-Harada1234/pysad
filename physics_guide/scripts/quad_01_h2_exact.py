# -*- coding: utf-8 -*-
"""四極の H2 の厳密解 (04_quad.md 式 (4.7)〜(4.14)) を、独立な数値計算 (行列指数・数値積分) で確かめ、pysad の行列・本体の写像と比べる。

(a) H2 の正準方程式 x' = px/p, px' = -K x の行列指数 (scipy.linalg.expm) と、閉じた式 (cos/cosh) の一致。K -> 0 の極限。
(b) pysad の Quadrupole.transfer_matrix (CALC の設計軌道上の行列) との一致。
(c) z の更新 (式 (4.14)): -1/2 ∫(x'^2+y'^2) ds の数値積分と、Fortran の式 (tsolqu の -0.25(...)) の一致。
(d) quadbody の本体 B (track6._qb_body) が (a)(c) の式と一致。
sympy は使わない (環境に無い)。導出の途中は 04_quad.md に手で書き、ここでは結果を独立な数値で確かめる。
"""
import math
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
import _common as C

C.start(__file__)
from pysad import track6 as t6


def Mx_closed(K, L, p):
    """集束 (K>0): (x, px) の行列 (px は p0 で規格化した正準運動量)"""
    w = math.sqrt(K / p)
    return np.array([[math.cos(w * L), math.sin(w * L) / (p * w)], [-p * w * math.sin(w * L), math.cos(w * L)]])


def My_closed(K, L, p):
    w = math.sqrt(K / p)
    return np.array([[math.cosh(w * L), math.sinh(w * L) / (p * w)], [p * w * math.sinh(w * L), math.cosh(w * L)]])


print('== (a) H2 の正準方程式の行列指数 exp(A L) と、閉じた式の比較 ==')
for p, K, L in ((1.0, 2.0, 0.3), (1.01, 5.0, 0.7), (0.98, 0.1, 2.0)):
    Ax = np.array([[0, 1 / p], [-K, 0]])         # (x, px): x' = px/p, px' = -K x
    Ay = np.array([[0, 1 / p], [K, 0]])          # (y, py): y' = py/p, py' = +K y
    ex = np.abs(expm(Ax * L) - Mx_closed(K, L, p)).max()
    ey = np.abs(expm(Ay * L) - My_closed(K, L, p)).max()
    print('  p = %.2f, K = %.1f, L = %.1f: |expm - 閉じた式| x: %.1e、y: %.1e、det = %.15f、%.15f'
          % (p, K, L, ex, ey, np.linalg.det(Mx_closed(K, L, p)), np.linalg.det(My_closed(K, L, p))))
print('  K -> 0 の極限 (K = 1e-12): x の行列 = [[%.12f, %.12f], [%.3e, %.12f]] (厳密なドリフトは [[1, L/p], [0, 1]]、L/p = %.12f)'
      % (*Mx_closed(1e-12, 0.3, 1.0)[0], Mx_closed(1e-12, 0.3, 1.0)[1, 0], Mx_closed(1e-12, 0.3, 1.0)[1, 1], 0.3))
print('  K < 0 の x は y の式で K -> |K| (cos と cosh の入れ替え)。K = 0 の前後で連続:  x の (0,0) 成分は K = +1e-9 で %.15f、K = -1e-9 で %.15f'
      % (Mx_closed(1e-9, 0.3, 1.0)[0, 0], My_closed(1e-9, 0.3, 1.0)[0, 0]))

print('\n== (b) pysad の Quadrupole.transfer_matrix との比較 (3 GeV、L = 0.3 m、K1 = ±0.6 (積分強度)、δ = 0 と 0.01) ==')
for K1 in (0.6, -0.6):
    for dl in (0.0, 0.01):
        ffs, el = C.make_quad(0.3, K1)
        M = el.transfer_matrix(dl)
        kk, pp = el.K1, 1.0 + dl                    # el.K1 は長さあたり [1/m^2]
        if kk > 0:
            Mx, My = Mx_closed(kk, 0.3, pp), My_closed(kk, 0.3, pp)
        else:
            Mx, My = My_closed(-kk, 0.3, pp), Mx_closed(-kk, 0.3, pp)
        err = max(abs(M[:2, :2] - Mx).max(), abs(M[2:4, 2:4] - My).max())
        print('  K1 = %+.1f, δ = %.2f: 最大差 %.2e' % (K1, dl, err))

print('\n== (c) z の更新: -1/2 ∫(x\'^2+y\'^2) ds の数値積分と、Fortran の式 -0.25(...) の比較 ==')


def z_fortran(w1, L, x0, q0, y0, r0):
    """tsolqu.f の z の式から dv の項を除いたもの (q = px/p, r = py/p は規格化した運動量)"""
    phi = w1 * L
    c1, s1, ch2, sh2 = math.cos(phi), math.sin(phi), math.cosh(phi), math.sinh(phi)
    dc1, xs1, dch2, xsh2 = c1 - 1, phi - s1, ch2 - 1, phi - sh2     # SAD の xsincos: (cos - 1, x - sin x)
    u1 = x0 * dc1 + q0 * s1 / w1
    u2 = -x0 * w1 * s1 + q0 * dc1
    v1 = y0 * dch2 + r0 * sh2 / w1
    v2 = y0 * w1 * sh2 + r0 * dch2
    ap = q0**2 + r0**2
    return -0.25 * (w1 * (x0**2 * xs1 - y0**2 * xsh2) + ap * L + u1 * (q0 + u2) + x0 * q0 * dc1 + v1 * (r0 + v2) + y0 * r0 * dch2)


def z_integral(w1, L, x0, q0, y0, r0):
    qx = lambda s: -x0 * w1 * math.sin(w1 * s) + q0 * math.cos(w1 * s)
    qy = lambda s: y0 * w1 * math.sinh(w1 * s) + r0 * math.cosh(w1 * s)
    return -0.5 * quad(lambda s: qx(s)**2 + qy(s)**2, 0, L, epsabs=1e-16, epsrel=1e-14)[0]


for args in ((1.3, 0.37, 0.011, -0.004, -0.007, 0.003), (2.0, 0.15, 0.02, 0.01, 0.015, -0.01)):
    a, b = z_integral(*args), z_fortran(*args)
    print('  数値積分 %.15e、Fortran 型の式 %.15e、差 %.2e' % (a, b, abs(a - b)))

print('\n== (d) quadbody の本体 B (track6._qb_body) と (a)(c) の式の比較 (乱数 5 点、δ = 0.003、dv = 1e-5、長さ 0.15 m) ==')
rng = np.random.default_rng(1)
worst = 0.0
for _ in range(5):
    x, px, y, py = rng.normal(0, 0.01, 4)
    dp, akk, ln, dv = 0.003, 2.0, 0.15, 1e-5
    out = t6._qb_body(x, px, y, py, 0.0, dv, t6._tzsetparam0(dp, ln, akk), ln, 0.0, 0.0, False)
    wv = math.sqrt(akk / (1 + dp))
    mx = Mx_closed(akk, ln, 1 + dp)
    my = My_closed(akk, ln, 1 + dp)
    pr = 1 + dp
    X = mx @ np.array([x, px * pr])          # 規格化した px (= px/p) を正準運動量に直して行列をかけ、戻す
    Y = my @ np.array([y, py * pr])
    ref = [X[0], X[1] / pr, Y[0], Y[1] / pr, z_integral(wv, ln, x, px, y, py) - dv * ln]
    err = max(abs(np.array(out) - np.array(ref)))
    worst = max(worst, err)
print('  最大差 %.2e' % worst)

print('\n== (e) 例: 3 GeV、L = 0.3 m、K1 = 0.6 (積分強度 [1/m]) の四極 (δ = 0、設計軌道) ==')
ffs, el = C.make_quad(0.3, 0.6)
M = el.transfer_matrix(0.0)
kk = el.K1
print('  長さあたりの強さ K = K1/L = %.3f [1/m^2]、w = sqrt(K) = %.6f [1/m]、位相 φ = w L = %.6f rad (= %.2f 度)、K1 L = φ^2 = %.3f' % (kk, math.sqrt(kk), math.sqrt(kk) * 0.3, math.degrees(math.sqrt(kk) * 0.3), 0.18))
print('  x の行列 [[%.6f, %.6f], [%.6f, %.6f]] (焦点距離 -1/M21 = %.4f m、薄肉近似 1/K1 = %.4f m)' % (M[0, 0], M[0, 1], M[1, 0], M[1, 1], -1 / M[1, 0], 1 / 0.6))
print('  y の行列 [[%.6f, %.6f], [%.6f, %.6f]]' % (M[2, 2], M[2, 3], M[3, 2], M[3, 3]))
print('  ∂z/∂δ = L/γ0^2 = %.4e (01 章の式 (1.9) から。CALC の T(5,6) の dvdp·L。γ0 = %.1f)' % (0.3 / (3e9 / 0.51099895e6)**2 * 1.0, 3e9 / 0.51099895e6))
