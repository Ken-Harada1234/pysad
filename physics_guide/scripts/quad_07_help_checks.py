# -*- coding: utf-8 -*-
"""ヘルプ (SADHelp.html) の QUAD・solution-H2 の式の確認 (04_quad.md 付録 A)。
 (a) 線形フリンジの b = K1 F2/L: a には 1/p1 があるが b には無い。q = px/p の単位の線形運動 x'' = -(K(s)/p) x の
     数値積分と比べ、b も 1/p が要ることを示す。
 (b) solution-H2 の z の式 (複素正規座標 u, pu, v, pv を使う): Bz = 0 で、z' = -(q²/2) - dv の積分と比べる。
 (c) F1 = sgn(a) sqrt(a) (a < 0 のとき根号が虚数)。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.linalg import expm
import _common as C

C.start(__file__)

print('== (a) 線形フリンジの a, b の p 依存 (K(s) = K0 (s/ℓ)^2 の立ち上がり、ℓ = 0.1 m、K0 = 1e-3、q = px/p の単位) ==')
ell = 0.1
s0 = 2 * ell / 3
prof = lambda s: 0.0 if s < 0 else (1.0 if s > ell else (s / ell)**2)
I = lambda n: quad(lambda s: (s - s0)**n * (prof(s) - (1.0 if s > s0 else 0.0)), 0, ell, points=[s0], epsabs=1e-16, epsrel=1e-13)[0]
I1, I2 = I(1), I(2)
F1 = math.sqrt(-24 * I1)
F2 = I2
K0, S1 = 1e-3, 0.2
print('   F1 = %.6f、F2 = %.6e' % (F1, F2))
for p in (1.0, 1.2, 0.8):
    Ke = K0 / p

    def rhs(s, v):
        k = Ke * prof(s)
        M = v.reshape(2, 2)
        return np.array([[M[1, 0], M[1, 1]], [-k * M[0, 0], -k * M[0, 1]]]).ravel()
    Mx = solve_ivp(rhs, (-S1, ell + S1), np.eye(2).ravel(), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1].reshape(2, 2)
    w = math.sqrt(Ke)
    drift = lambda L_: np.array([[1, L_], [0, 1]])
    body = lambda L_: np.array([[math.cos(w * L_), math.sin(w * L_) / w], [-w * math.sin(w * L_), math.cos(w * L_)]])
    fr = lambda a, b: np.array([[math.exp(a), b], [0, math.exp(-a)]])
    Lin, Lout = S1 + s0, S1 + ell - s0
    a_h = -K0 * F1 * abs(F1) / (24 * p)             # ヘルプの a (1/p1 あり)
    e = {}
    for name, b_ in (('b = K0 F2/p (導出・K64 の _tqlfr)', K0 * F2 / p), ('b = K0 F2 (ヘルプの式のまま)', K0 * F2)):
        m = body(Lout) @ fr(a_h, b_) @ drift(Lin)
        e[name] = np.abs(Mx - m).max()
    m0 = body(Lout) @ drift(Lin)
    print('   p = %.1f: フリンジなし %.2e、' % (p, np.abs(Mx - m0).max()) + '、'.join('%s: %.2e' % (k, v) for k, v in e.items()))
print('   → p ≠ 1 では、b にも 1/p が要る (ヘルプの b = K1 F2/L は 1/p1 を略している)。F2 = 0 の端 (直線・tanh) では効かない。')

print('\n== (b) solution-H2 の z (Bz = 0): ヘルプの複素正規座標の式 ==')
K, p, ell2 = 0.9, 1.01, 0.4
dvv = -0.0123   # 任意の dv (Δv = -dv)
w = math.sqrt(K / p)
w1 = w2 = w
wp = w1 + w2
Uu = 1.0 / (w1**2 + w2**2)
s = math.sqrt(Uu)
r = math.sqrt(w1**2 - w2**2)            # = 0
def Rhelp(fix):
    """ヘルプの 4×4 の複素正規座標への行列 (Bz = 0)。fix = True は、4 行 3 列の -i p/4 (w1+w2)² の i を除いたもの"""
    c43 = -p / 4 * wp**2 if fix else -1j * p / 4 * wp**2
    return s * np.array([
        [wp / 2, -1j / p, -r / 2, (1 / p) * math.sqrt((w1 - w2) / wp)],
        [-1j * p / 4 * wp**2, wp / 2, -p * (w1 - w2) / 4 * r, -(w1 - w2) / 2 * math.sqrt((w1 - w2) / wp)],
        [-r / 2, (1 / p) * math.sqrt((w1 - w2) / wp), wp / 2, 1 / p],
        [-p * (w1 - w2) / 4 * r, -1j * (w1 - w2) / 2 * math.sqrt((w1 - w2) / wp), c43, wp / 2]])


Jm = np.zeros((4, 4)); Jm[0, 1] = Jm[2, 3] = 1; Jm -= Jm.T
for fix in (False, True):
    Rt = Rhelp(fix)
    print('   行列 (%s): |R^T J R - J| = %.2e (正準変換なら 0)' % ('原文のまま' if not fix else '4 行 3 列の i を除く', np.abs(Rt.T @ Jm @ Rt - Jm).max()))
R = Rhelp(True)
c0 = np.array([1.2e-3, 0.7e-3, -0.9e-3, 0.4e-3])         # x, px, y, py
u, pu, v, pv = R @ c0
dz_help1 = (-1j * u * pu * (-w1 / (2 * p)) - v * pv * (-w2 / (2 * p)) + (-dvv)) * ell2
dz_help2 = (Uu / p) * (1j * w1**3 * u * pu + w2**3 * v * pv + (-dvv)) * ell2
# 厳密: z' = -(q²/2) - dv、q = (px, py)/p。x' = px/p、px' = -K x、y' = py/p、py' = K y
A = np.array([[0, 1 / p, 0, 0], [-K, 0, 0, 0], [0, 0, 0, 1 / p], [0, 0, K, 0]])
f = lambda sx: expm(A * sx) @ c0
zq = quad(lambda sx: -0.5 * ((f(sx)[1] / p)**2 + (f(sx)[3] / p)**2), 0, ell2, epsabs=1e-18, epsrel=1e-12)[0] - dvv * ell2
print('   厳密 Δz = %.10e (うち 2 次の部分 -∫q²/2 ds = %.6e、-dv ℓ = %.6e)' % (zq, zq + dvv * ell2, -dvv * ell2))
q1 = (-1j * u * pu * (-w1 / (2 * p)) - v * pv * (-w2 / (2 * p))) * ell2
q2 = (Uu / p) * (1j * w1**3 * u * pu + w2**3 * v * pv) * ell2
print('   ヘルプの 2 次の部分: 1 行目 %.6e、2 行目 (Δv を括弧の外に出した形) %.6e' % (q1.real, q2.real))
print('   ヘルプの 1 行目 (-i u pu ∂w1/∂p - v pv ∂w2/∂p + Δv) ℓ (Δv = -dv): %.10e + %.1e i' % (dz_help1.real, dz_help1.imag))
print('   ヘルプの 2 行目 (U/p)(i w1³ u0 pu0 + w2³ v0 pv0 + Δv) ℓ: %.10e + %.1e i' % (dz_help2.real, dz_help2.imag))
print('   → 差 (1 行目) %.1e、(2 行目) %.1e' % (abs(dz_help1 - zq), abs(dz_help2 - zq)))
print('   結論: Δv は原文に定義が無く、-dv = v/v0 - 1 と読んで dv の部分は一致する。2 次の部分 (u pu, v pv の項) は、厳密な -∫q²/2 ds と一致しない')
print('         (符号も大きさも違う: ヘルプの式は H_2u を (u, pu, v, pv) を固定して p で微分した形で、正規座標への変換自体の p 依存の寄与と、')
print('          振動する成分が入っていない可能性。未特定)。2 行目 (U/p)(… + Δv) は Δv まで U/p がかかるので、括弧の位置も誤りとみられる (2 行目全体の差が大きい)。')
print('         行列は、4 行 3 列の i を除いた形 (正準になる形) で使った。')
print('   Bz ≠ 0 の行列は、原文のままでも、1 項の符号・i を変えるだけの修正でも正準でなかった (原因は未特定)。solution-H2 の x, px, y, py の式そのものは 09 章 sol_01 で行列指数と一致。')

print('\n== (c) F1 = sgn(a) sqrt(a)、a = 24 (I0²/2 - I1) ==')
for a_ in (4.0e-3, -4.0e-3):
    try:
        val = math.copysign(math.sqrt(a_), a_)
    except ValueError:
        val = 'sqrt(負) は実数でない (sgn(a) sqrt|a| と読む)'
    print('   a = %+.1e → %s' % (a_, val))
