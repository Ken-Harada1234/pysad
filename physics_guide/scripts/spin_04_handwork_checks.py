# -*- coding: utf-8 -*-
r"""スピン (15_spin.md) の、紙の上で追う途中式の検算 (sympy・数値)。本文の手導出の各行を確かめる。

(a) 場の変換 F' = Λ F Λ^T の成分 (E ≠ 0 も含めて): B'_x = B_x、B'_y = γ(B_y + β E_z)、B'_z = γ(B_z - β E_y)  (§3.1)
(b) Thomas 歳差の手導出 (§3.2): P = Λ(v+δv) Λ(-v) を δv の 1 次まで成分で書く。
    空間 3×3 部分 = I + [w]× (対称部分は 0)、w = (γ-1)/β (0, δv_z, -δv_y)、時間成分の部分は対称 (純ブースト)
(c) 向きの確認: 一様な円運動 (トルクなし) の Fermi-Walker 輸送で、Λ(v(t)) S_lab(t) が -(γ-1)ω z^ で回ること
(d) 規格化した軸の内積 n^(x)·n^(x+ξ) = 1 - |b|²/2 の展開 (§6.2)
(e) 2 次の相関和 (§6.5): n0·A2 A1 n0 = -Re(ω_c2 conj(ω_c1))、(e2 + i e3)^T R0^{-1} = e^{iμ} (e2 + i e3)^T
(f) 中点則の誤差 (§4.2): Magnus 展開の 2 項で、1 刻みの回転ベクトルの誤差 = Δ³/24 ω'' + Δ³/12 (ω' × ω_m)
"""
import math
import numpy as np
import sympy as sp
import scipy.linalg as sl
import scipy.integrate as si
import _common as C

C.start(__file__)

# ---------------------------------------------------------------- (a)
print('== (a) F\' = Λ F Λ^T (x 軸方向のブースト、E ≠ 0 も入れる) ==')
be, ga = sp.symbols('beta gamma', positive=True)
Ex, Ey, Ez, Bx, By, Bz = sp.symbols('Ex Ey Ez Bx By Bz', real=True)
# 座標 (t, x, y, z)、F^{0i} = -E_i、F^{ij} = -ε_{ijk} B_k
F = sp.Matrix([[0, -Ex, -Ey, -Ez],
               [Ex, 0, -Bz, By],
               [Ey, Bz, 0, -Bx],
               [Ez, -By, Bx, 0]])
L = sp.Matrix([[ga, -ga * be, 0, 0], [-ga * be, ga, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
Fp = (L * F * L.T).subs(ga, 1 / sp.sqrt(1 - be**2))
Fp = Fp.applyfunc(sp.simplify)
Bp = [-Fp[2, 3], -Fp[3, 1], -Fp[1, 2]]
g_ = 1 / sp.sqrt(1 - be**2)
expect = [Bx, g_ * (By + be * Ez), g_ * (Bz - be * Ey)]
print("B' - 期待 [Bx, γ(By + βEz), γ(Bz - βEy)] =", [sp.simplify(b - e) for b, e in zip(Bp, expect)])
print("  (β × E)_y = -β Ez、(β × E)_z = β Ey なので、γ(B - β×E)_⊥ と同じ。E = 0 では B'_⊥ = γ B_⊥、B'_∥ = B_∥")

# ---------------------------------------------------------------- (b)
print('\n== (b) Thomas 歳差: P = Λ(v + δv) Λ(-v) の 1 次 (v = (β,0,0)、δv = (0, ε_y, ε_z)) ==')
ey, ez, t = sp.symbols('eps_y eps_z t', real=True)


def boost(u):
    """実験室系の座標 -> 速度 u で動く系の座標 (t, x, y, z)、c = 1"""
    u = sp.Matrix(u)
    u2 = u.dot(u)
    gm = 1 / sp.sqrt(1 - u2)
    M = sp.zeros(4, 4)
    M[0, 0] = gm
    for i in range(3):
        M[0, i + 1] = -gm * u[i]
        M[i + 1, 0] = -gm * u[i]
        for j in range(3):
            M[i + 1, j + 1] = (1 if i == j else 0) + (gm - 1) * u[i] * u[j] / u2
    return M


v = [be, 0, 0]
P = boost([be, t * ey, t * ez]) * boost([-be, 0, 0])      # t = 微小量のカウンタ (δv = t ε)
P1 = P.applyfunc(lambda e: sp.simplify(sp.diff(e, t).subs(t, 0)))
P0 = P.applyfunc(lambda e: sp.simplify(e.subs(t, 0)))
print('P(t=0) = 単位行列か:', sp.simplify(P0 - sp.eye(4)) == sp.zeros(4, 4))
gam = 1 / sp.sqrt(1 - be**2)
print('1 次の係数 dP/dt (そのまま):')
for i_ in range(4):
    print('  ', [str(sp.simplify(P1[i_, j_])) for j_ in range(4)])
S3 = P1[1:, 1:]
print('空間 3×3 の対称部分 (S3 + S3^T)/2 =', ((S3 + S3.T) / 2).applyfunc(sp.simplify).tolist())
print('手計算の値: P_yx = ε_y (1-γ)/β、P_xy = ε_y (γ-1)/β (z も同様):')
print('  P_yx - ε_y(1-γ)/β =', sp.simplify(P1[2, 1] - ey * (1 - gam) / be),
      '、P_xy - ε_y(γ-1)/β =', sp.simplify(P1[1, 2] - ey * (gam - 1) / be))
print('  P_zx - ε_z(1-γ)/β =', sp.simplify(P1[3, 1] - ez * (1 - gam) / be),
      '、P_xz - ε_z(γ-1)/β =', sp.simplify(P1[1, 3] - ez * (gam - 1) / be))
# 回転ベクトル w: ([w]×)_{yx} = w_z、([w]×)_{zx} = -w_y
wz = S3[1, 0]
wy = -S3[2, 0]
wx = S3[2, 1]       # ([w]×)_{zy} = w_x
print('w = (%s, %s, %s)' % (sp.simplify(wx), sp.simplify(wy), sp.simplify(wz)))
print('  期待 (γ-1)/β (0, ε_z, -ε_y) との差:',
      [sp.simplify(wx), sp.simplify(wy - (gam - 1) / be * ez), sp.simplify(wz + (gam - 1) / be * ey)])
Kb = P1[0, 1:]
Kc = P1[1:, 0]
print('時間-空間の部分: 行 0 =', [sp.simplify(x) for x in Kb], '、列 0 =', [sp.simplify(x) for x in Kc],
      '(対称 = 純ブースト、ε_y・ε_z の成分は -γε)')
# P = R Q (Q: 隣り合う静止系の間の純ブースト) の確認: P Q^{-1} の空間部分が回転、時間成分が 0
Q = sp.eye(4) + t * sp.Matrix(4, 4, lambda i, j: (P1[i, j] if (i == 0 or j == 0) else 0))
R = (P * Q.inv()).applyfunc(lambda e: sp.simplify(sp.diff(e, t).subs(t, 0)))
print('d(P Q^{-1})/dt (1 次) の時間行・列は 0、空間部分は反対称か:',
      [sp.simplify(x) for x in R[0, :]], [sp.simplify(x) for x in R[:, 0]],
      sp.simplify(R[1:, 1:] + R[1:, 1:].T) == sp.zeros(3, 3))
# w = γ²/(γ+1) δv × v
vv = sp.Matrix(v)
dv = sp.Matrix([0, ey, ez])
wform = (gam**2 / (gam + 1)) * dv.cross(vv)
print('w - γ²/(γ+1) δv × v =', [sp.simplify(a - b) for a, b in zip([wx, wy, wz], wform)])

# ---------------------------------------------------------------- (c)
print('\n== (c) Thomas 歳差の向き: トルクのない円運動の Fermi-Walker 輸送を数値で ==')
beta, om = 0.6, 1.0
gamma = 1 / math.sqrt(1 - beta**2)


def u4(th):
    return gamma * np.array([1, -beta * math.sin(th), beta * math.cos(th), 0])


def a4(th):          # 4 元加速度 (反変成分)。a_3 = ω×v = -βω (cosθ, sinθ)、a^μ = (0, γ² a_3)
    return np.array([0, -gamma**2 * beta * om * math.cos(th), -gamma**2 * beta * om * math.sin(th), 0])


eta = np.diag([1, -1, -1, -1])


def rhs(t_, S):
    th = om * t_
    u, a = u4(th), a4(th)
    dS_dtau = -(S @ eta @ a) * u          # D S/dτ = -(a·S) u (S·u = 0 を保つ Fermi-Walker 輸送、署名 + - - -)
    return dS_dtau / gamma


def boost_np(vel):
    b2 = vel @ vel
    g = 1 / math.sqrt(1 - b2)
    M = np.zeros((4, 4))
    M[0, 0] = g
    M[0, 1:] = -g * vel
    M[1:, 0] = -g * vel
    M[1:, 1:] = np.eye(3) + (g - 1) * np.outer(vel, vel) / b2
    return M


v0 = beta * np.array([0, 1, 0.0])
S0 = np.linalg.inv(boost_np(v0)) @ np.array([0, 1, 0, 0.0])          # 静止系で s = x^ の向き
T = 2 * math.pi / om * 0.25
sol = si.solve_ivp(rhs, [0, T], S0, rtol=1e-12, atol=1e-13)
SB = boost_np(beta * np.array([-math.sin(om * T), math.cos(om * T), 0])) @ sol.y[:, -1]
ang = math.atan2(SB[2], SB[1]) - 0.0
print('t = T/4 (軌道を 1/4 周) での、ブースト軸の成分 S_B = Λ(v(t)) S_lab(t) の空間部分:', SB[1:], '(時間成分 %.2e)' % SB[0])
print('  z 軸まわりの回転角 = %.6f rad、予想 -(γ-1) ω T = %.6f rad' % (ang, -(gamma - 1) * om * T))
print('  → 回転は軌道の向きと逆 (Thomas 歳差)。ω_T = γ²/(γ+1) a × v の向きと一致 (a × v の z 成分 = -β²ω)')

# ---------------------------------------------------------------- (d)
print('\n== (d) 規格化した軸の内積 (§6.2) ==')
a1, a2, b1, b2, h = sp.symbols('a1 a2 b1 b2 h', real=True)
# n0 = ŷ のとき a, b は xz 平面 (成分 1, 2)。|a| ~ h、|b| ~ 1 次までの展開のため a -> h a
av = sp.Matrix([h * a1, h * a2])
bv = sp.Matrix([b1, b2])
num = 1 + av.dot(av) + av.dot(bv)
den = sp.sqrt(1 + av.dot(av)) * sp.sqrt(1 + (av + bv).dot(av + bv))
ser = sp.series(num / den, h, 0, 2).removeO()
print('h の 1 次までの展開:', sp.simplify(ser))
print('  h^0 の項 - (1 - |b|²/2 の厳密な形 1/√(1+b²)) =', sp.simplify(ser.subs(h, 0) - 1 / sp.sqrt(1 + bv.dot(bv))))
print('  1/√(1+b²) の b の 2 次まで =', sp.series(1 / sp.sqrt(1 + b1**2), b1, 0, 4))

# ---------------------------------------------------------------- (e)
print('\n== (e) 2 次の相関和の骨格 (§6.5) ==')
w1 = sp.Matrix(sp.symbols('w1x w1y w1z', real=True))
w2 = sp.Matrix(sp.symbols('w2x w2y w2z', real=True))


def cross_mat(w):
    return sp.Matrix([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


n0 = sp.Matrix([0, 1, 0])
lhs = (n0.T * cross_mat(w2) * cross_mat(w1) * n0)[0]
# n0 = ŷ のとき、垂直面は (z, x)。e2 = ẑ、e3 = n0 × e2 = x̂
e2 = sp.Matrix([0, 0, 1])
e3 = n0.cross(e2)
wc1 = e2.dot(w1) + sp.I * e3.dot(w1)
wc2 = e2.dot(w2) + sp.I * e3.dot(w2)
print('n0·([w2]× [w1]× n0) + Re(ω_c2 conj ω_c1) =', sp.simplify(lhs + sp.re(sp.expand(wc2 * sp.conjugate(wc1)))))
print('  (ω_c = e2·ω + i e3·ω、e3 = n0 × e2。ω_∥ = n0·ω は寄与しない)')
mu = sp.symbols('mu', real=True)
# K64 の約束 R0 = exp(-μ [n0]×) : R0 e2 = cosμ e2 - sinμ e3、R0 e3 = sinμ e2 + cosμ e3
R0 = sp.simplify((-mu * cross_mat(n0)).exp())
ec = (e2 + sp.I * e3).T
print('(e2 + i e3)^T R0^{-1} - e^{iμ} (e2 + i e3)^T =',
      [sp.simplify(sp.expand_complex(sp.expand(x.rewrite(sp.cos)))) for x in (ec * R0.inv() - sp.exp(sp.I * mu) * ec)])
print('  → R0^{-(n+1)} D x の複素成分 = e^{i(n+1)μ} d x_n (d = (e2 + i e3)^T D)')
print('S̃_{n+1} = exp(-[R0^{-(n+1)} D x_n]×) S̃_n の導出は、B exp(-[v]×) B^{-1} = exp(-[Bv]×) (B = R0^{n+1}) による')

# ---------------------------------------------------------------- (f)
print('\n== (f) 中点則: Magnus 展開の 2 項で 1 刻みの回転ベクトルの誤差 ==')
rng = np.random.default_rng(3)
wm = rng.normal(size=3)
wp = rng.normal(size=3)
wpp = rng.normal(size=3)


def omega(s):
    return wm + wp * s + 0.5 * wpp * s**2


def hat(w):
    return np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])


def exact_prop(D):
    f = lambda s, y: (hat(omega(s - D / 2)) @ y.reshape(3, 3)).ravel()
    sol = si.solve_ivp(f, [0, D], np.eye(3).ravel(), rtol=1e-13, atol=1e-15)
    return sol.y[:, -1].reshape(3, 3)


def rotvec(Rm):
    L_ = sl.logm(Rm).real
    return np.array([L_[2, 1], L_[0, 2], L_[1, 0]])


print("  D        |Ω_exact - D ω_m|     予想 D³(ω''/24 + ω'×ω_m/12) との差      比 (D を半分にしたときの誤差の比)")
prev = None
for D in (0.4, 0.2, 0.1, 0.05):
    Om = rotvec(exact_prop(D))
    err = Om - D * wm
    pred = D**3 * (wpp / 24 + np.cross(wp, wm) / 12)
    n_err = np.linalg.norm(err)
    print('  %.3f  %.4e   %.3e' % (D, n_err, np.linalg.norm(err - pred)), '  %s' % ('%.2f' % (prev / n_err) if prev else ''))
    prev = n_err
print('  → 1 刻みの誤差は D³ (比 8)。N = 1/D 刻みを足し合わせると全体は D² (比 4)。')
