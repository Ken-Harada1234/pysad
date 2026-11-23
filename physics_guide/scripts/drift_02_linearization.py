# -*- coding: utf-8 -*-
"""閉軌道 (COD) のまわりの 6×6 (02_drift.md §7〜§8、式 (2.13)〜(2.20)) を確かめる。

(a) sympy: 厳密な写像 (x += px L/pz、…、z += -(px^2+py^2)L/((p+pz)pz) - dv(δ) L) を正準変数で偏微分し、tdrife1 の式 (toplvl.f:1368) と比べる。
(b) pysad の行列 (lattice._tdrife1_mat、elkern._jac_drift) と、pysad の追跡 (DRIFT 1 台) の数値微分の比較。
(c) シンプレクティック条件 M^T J M - J と、z の行をシンプレクティック条件から組む式。
(d) 設計軌道 (COD = 0) での値 (R12 = L/p、R56 = L h0/h1^3 ≈ L/γ0^2) と、COD が傾いたときの R16 = -px p L/pz^3 (経路長の分散)。
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad import lattice as lat
from pysad import elkern

MC2 = 0.51099895e6
p0 = 3.0e9 / MC2
h0 = math.sqrt(1 + p0**2)
L = 1.0
Jm = np.zeros((6, 6))
for i in range(3):
    Jm[2 * i, 2 * i + 1] = 1.0
    Jm[2 * i + 1, 2 * i] = -1.0

print('== (a) sympy: 厳密な写像の偏微分と tdrife1 の式 ==')
X, PX, Y, PY, Z, D = sp.symbols('x p_x y p_y z delta', real=True)
Lsym = sp.symbols('L', positive=True)
P0s, H0s = sp.symbols('p0 h0', positive=True)
pr = 1 + D
pz = sp.sqrt(pr**2 - PX**2 - PY**2)
a = PX**2 + PY**2
h1 = sp.sqrt(1 + (P0s * pr)**2)
dv = 1 - pr * H0s / h1
mapf = sp.Matrix([X + PX * Lsym / pz, PX, Y + PY * Lsym / pz, PY, Z - a / ((pr + pz) * pz) * Lsym - dv * Lsym, D])
Jsym = mapf.jacobian([X, PX, Y, PY, Z, D])
# tdrife1 の式 (toplvl.f:1368-1390)。ale = L/pz、alz = ale/pz^2、dvdp = h0/h1^3
ale = Lsym / pz
alz = ale / pz**2
dvdp = H0s / h1**3
T = sp.eye(6)
T[0, 1] = ale + PX**2 * alz
T[0, 3] = PX * PY * alz
T[0, 5] = -PX * pr * alz
T[2, 1] = T[0, 3]
T[2, 3] = ale + PY**2 * alz
T[2, 5] = -PY * pr * alz
T[4, 1] = T[0, 5]
T[4, 3] = T[2, 5]
T[4, 5] = dvdp * Lsym + a * alz
vals = {X: 0.001, PX: 0.013, Y: -0.002, PY: -0.021, Z: 0.0, D: 0.004, Lsym: 1.3, P0s: p0, H0s: h0}
diff = (Jsym - T).subs(vals)
print('  ヤコビアンと tdrife1 の式の差 (数値代入、最大) = %.1e' % max(abs(float(v)) for v in diff))
nonzero = [(i + 1, j + 1) for i in range(6) for j in range(6) if sp.simplify(Jsym[i, j] - T[i, j]) != 0]
print('  記号の差が 0 でない成分 (simplify): %s' % (nonzero if nonzero else 'なし'))
print('  z の行の 2 つの成分を、手で導いた式と比べる (記号の差):')
print('    ∂z/∂px - (-p px L/pz^3) = %s' % sp.simplify(Jsym[4, 1] + pr * PX * Lsym / pz**3))
print('    ∂z/∂δ - (L (px^2+py^2)/pz^3 + L h0/h1^3) = %s' % sp.simplify(Jsym[4, 5] - (Lsym * a / pz**3 + Lsym * H0s / h1**3)))

print('\n== (b) pysad の行列と、追跡の数値微分 (3 GeV、L = 1 m) ==')
ffs, el = C.make_drift(L)
worst1 = worst2 = 0.0
for c in (np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
          np.array([0.001, 0.01, -0.002, 0.02, 0.0, 0.003]),
          np.array([0.0, 0.1, 0.0, -0.05, 0.0, -0.02]),
          np.array([0.002, 0.3, 0.0, 0.2, 0.001, 0.01])):
    pr_ = 1 + c[5]
    dvv = elkern.__dict__.get('_dv', None)
    # 追跡の数値微分 (中心差分)
    Jt = np.zeros((6, 6))
    h = 1e-6
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        Jt[:, j] = (C.track_one(ffs, cp) - C.track_one(ffs, cm)) / (2 * h)
    # pysad の CALC の行列 (dvdp は -∂dv/∂δ = h0/h1^3)
    h1n = math.sqrt(1 + (p0 * pr_)**2)
    Mt = lat._tdrife1_mat(c, 0.0, h0 / h1n**3, L)
    ok, Me = elkern._jac_drift(c, L, h0 / h1n**3)
    e1 = np.abs(Mt - Jt).max()
    e2 = np.abs(Me - Jt).max()
    worst1, worst2 = max(worst1, e1), max(worst2, e2)
    print('  cod = %s: lattice._tdrife1_mat - 追跡の数値微分 %.1e、elkern._jac_drift - 追跡の数値微分 %.1e (数値微分の精度 ~1e-9)' %
          (np.array2string(c, precision=3), e1, e2))

print('\n== (c) シンプレクティック性と z の行 ==')
for c in (np.array([0.001, 0.01, -0.002, 0.02, 0.0, 0.003]), np.array([0.002, 0.3, 0.0, 0.2, 0.001, 0.01])):
    pr_ = 1 + c[5]
    h1n = math.sqrt(1 + (p0 * pr_)**2)
    M = lat._tdrife1_mat(c, 0.0, h0 / h1n**3, L)
    print('  cod = %s: |M^T J M - J| = %.1e' % (np.array2string(c, precision=3), np.abs(M.T @ Jm @ M - Jm).max()))
    # z の行をシンプレクティック条件で組む: T(5,j) = -T(1,j) T(2,6) + T(2,j) T(1,6) - T(3,j) T(4,6) + T(4,j) T(3,6)  (j = 1..4)
    zrow = np.zeros(4)
    for j in range(4):
        zrow[j] = -M[0, j] * M[1, 5] + M[1, j] * M[0, 5] - M[2, j] * M[3, 5] + M[3, j] * M[2, 5]
    print('    z の行 T(5,1..4) の式 %s、シンプレクティック条件で組んだ値 %s、差 %.1e' % (np.array2string(M[4, :4], precision=6), np.array2string(zrow, precision=6), np.abs(M[4, :4] - zrow).max()))
    # (5,5) = 1、(6,6) = 1、δ の行は恒等
    print('    T(5,5) = %.1f、T(6,6) = %.1f、δ の行 %s' % (M[4, 4], M[5, 5], M[5, :]))

print('\n== (d) 設計軌道と、COD が傾いたときの成分 ==')
Mz = lat._tdrife1_mat(np.zeros(6), 0.0, h0 / h0**3, L)
print('  COD = 0: R12 = %.6f (= L)、R34 = %.6f、R16 = %.1e、R56 = %.6e (= L h0/h1^3 = L/γ0^2 = %.6e)' % (Mz[0, 1], Mz[2, 3], Mz[0, 5], Mz[4, 5], L / h0**2))
for dp in (0.01, -0.01):
    pr_ = 1 + dp
    h1n = math.sqrt(1 + (p0 * pr_)**2)
    M = lat._tdrife1_mat(np.array([0, 0, 0, 0, 0, dp]), 0.0, h0 / h1n**3, L)
    print('  COD = (0, 0, 0, 0, 0, δ = %+.2f): R12 = %.6f (= L/p = %.6f)、R56 = %.6e (L/γ0^2 = %.6e との比 %.5f)' % (dp, M[0, 1], L / pr_, M[4, 5], L / h0**2, M[4, 5] / (L / h0**2)))
for px_ in (1e-3, 1e-2):
    c = np.array([0, px_, 0, 0, 0, 0.0])
    M = lat._tdrife1_mat(c, 0.0, h0 / h0**3, L)
    print('  COD = (0, px = %.0e, 0, 0, 0, 0): R12 = %.8f、R16 = %.4e (= -px p L/pz^3 ≈ -px L = %.4e)、R52 = %.4e (= R16)、R56 = %.6e (L/γ0^2 + px^2 L = %.6e)' %
          (px_, M[0, 1], M[0, 5], -px_ * L, M[4, 1], M[4, 5], L / h0**2 + px_**2 * L))
