# -*- coding: utf-8 -*-
r"""スピン (15_spin.md) の途中式の確認 (sympy・数値)。

(a) 磁場の Lorentz 変換 B' = γ[B - γ/(γ+1) β (β·B)] (E = 0)  (15.1)
(b) Thomas 歳差: Λ(v+δv) Λ(v)^{-1} の空間の回転部分 = γ²/(γ+1) (δv × v) 型  (15.3 の前提)
(c) Rodrigues の式: [g]×^3 = -|g|² [g]×、exp(-[g]×) S = cos|g| S + sin|g|/|g| S×g + (1-cos|g|)/|g|² (S·g) g  (15.10)
(d) (15.11) の恒等式: (sz - s sx')/c = c sz - s sx
(e) スピン反転率から P_ST と 1/τ_ST を出す算術  (15.17)
(f) γ ∂/∂γ = γ²/(p0 p) ∂/∂δ  (7.3)
(g) D の積算 (15.12): スライスごとの D <- R_i (D + G_i M_{0->i}) が、有限差分の微分と合う。誤差はスライスを細かくすると 1 次で減る
(h) 逐次の N の式 (15.13') が、1 周の式 (15.13) に畳み込まれる
(i) ガウス分布のスピンチューン広がりを入れた和の積分 ∫0^∞ exp(-c t - a² t²) dt = (√π/2a) exp(c²/4a²) erfc(c/2a)、
    および exp(cb²)·erfc(cb) の倍精度での桁あふれ (K64 の spintrm の B55 の NaN の機構)
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())

# ---------------------------------------------------------------- (a)
print('== (a) 磁場の Lorentz 変換 (E = 0) ==')
b1, b2, b3 = sp.symbols('b1 b2 b3', real=True)
B1, B2, B3 = sp.symbols('B1 B2 B3', real=True)
bet = sp.Matrix([b1, b2, b3])
b2s = bet.dot(bet)
gam = 1 / sp.sqrt(1 - b2s)


def boost(beta, gm):
    """実験室系 -> 速度 beta で動く系 (4×4)。座標 (t, x, y, z)、c = 1。"""
    L = sp.zeros(4, 4)
    L[0, 0] = gm
    for i in range(3):
        L[0, i + 1] = -gm * beta[i]
        L[i + 1, 0] = -gm * beta[i]
        for j in range(3):
            L[i + 1, j + 1] = (1 if i == j else 0) + (gm - 1) * beta[i] * beta[j] / b2s
    return L


Bv = sp.Matrix([B1, B2, B3])
F = sp.zeros(4, 4)    # F^{ij} = -ε_{ijk} B_k (E = 0)
F[1, 2], F[2, 1] = -B3, B3
F[2, 3], F[3, 2] = -B1, B1
F[3, 1], F[1, 3] = -B2, B2
L = boost(bet, gam)
Fp = L * F * L.T
Bp = sp.Matrix([-Fp[2, 3], -Fp[3, 1], -Fp[1, 2]])
Bexp = gam * (Bv - gam / (gam + 1) * bet * bet.dot(Bv))
res = []
for k in range(3):
    d = sp.simplify(sp.together(Bp[k] - Bexp[k]))
    res.append(d)
print('B\' (F\' = Λ F Λ^T から) - γ[B - γ/(γ+1) β(β·B)] の各成分:', res)
# 縦・横
Bpar = sp.Matrix([0, 0, B3])
print('β = (0,0,β) のとき B\'∥ = B∥、B\'⊥ = γ B⊥ は (15.1) の式から:')
bb = sp.symbols('beta', positive=True)
g_ = 1 / sp.sqrt(1 - bb**2)
Bexp2 = g_ * (Bv - g_ / (g_ + 1) * sp.Matrix([0, 0, bb]) * (bb * B3))
print('  ', [sp.simplify(x) for x in Bexp2], '(期待: γB1, γB2, B3)')

# ---------------------------------------------------------------- (b)
print('\n== (b) Thomas 歳差: 速度が v から v+δv に変わるときの Λ(v+δv)Λ(v)^{-1} の回転部分 ==')
be = sp.symbols('beta', positive=True)
dx_, dy_, dz_ = sp.symbols('dx dy dz', real=True)
eps = sp.symbols('eps', real=True)
gm = 1 / sp.sqrt(1 - be**2)
v0 = sp.Matrix([be, 0, 0])


def boost2(v):
    s2 = v.dot(v)
    g = 1 / sp.sqrt(1 - s2)
    M = sp.zeros(4, 4)
    M[0, 0] = g
    for i in range(3):
        M[0, i + 1] = -g * v[i]
        M[i + 1, 0] = -g * v[i]
        for j in range(3):
            M[i + 1, j + 1] = (1 if i == j else 0) + (g - 1) * v[i] * v[j] / s2
    return M


dv = sp.Matrix([dx_, dy_, dz_])
L0 = boost2(v0)
L1 = boost2(v0 + eps * dv)
E1 = sp.simplify((sp.diff(L1, eps).subs(eps, 0)) * L0.inv())     # 1 次の生成子 (Lie 代数の元)
spatial = E1[1:, 1:]
anti = sp.simplify((spatial - spatial.T) / 2)
# 回転ベクトル w: [w]× = R の無限小部分 -> anti_{ij} = -ε_{ijk} w_k、w = (anti[2,1], anti[0,2], anti[1,0])
w = sp.Matrix([anti[2, 1], anti[0, 2], anti[1, 0]])
w = sp.simplify(w)
cross_dv_v = dv.cross(v0)
ref = sp.simplify(gm**2 / (gm + 1) * cross_dv_v)
print('Λ(v+δv)Λ(v)^{-1} の空間の回転ベクトル w =', list(w))
print('γ²/(γ+1) δv×v =', list(ref))
print('差 (w - γ²/(γ+1) δv×v):', [sp.simplify(x) for x in (w - ref)], '、和 (w + ...):', [sp.simplify(x) for x in (w + ref)])

# ---------------------------------------------------------------- (c)
print('\n== (c) Rodrigues の式 ==')
g1, g2, g3 = sp.symbols('g1 g2 g3', real=True)
gv = sp.Matrix([g1, g2, g3])
K = sp.Matrix([[0, -g3, g2], [g3, 0, -g1], [-g2, g1, 0]])     # [g]× v = g × v
gg = gv.dot(gv)
print('K³ + |g|² K = ', sp.simplify(K**3 + gg * K).tolist())
print('K² - (g gᵀ - |g|² I) = ', sp.simplify(K**2 - (gv * gv.T - gg * sp.eye(3))).tolist())
th = sp.sqrt(gg)
Rod = sp.eye(3) - sp.sin(th) / th * K + (1 - sp.cos(th)) / gg * K**2
S1, S2, S3 = sp.symbols('S1 S2 S3', real=True)
Sv = sp.Matrix([S1, S2, S3])
form = sp.cos(th) * Sv + sp.sin(th) / th * Sv.cross(gv) + (1 - sp.cos(th)) / gg * Sv.dot(gv) * gv
print('Rod·S - (15.10 の右辺) = ', [sp.simplify(x) for x in (Rod * Sv - form)])
import scipy.linalg as sl
rng = np.random.default_rng(3)
gn = rng.normal(size=3)
Kn = np.array([[0, -gn[2], gn[1]], [gn[2], 0, -gn[0]], [-gn[1], gn[0], 0]])
Rn = np.array(Rod.subs({g1: gn[0], g2: gn[1], g3: gn[2]}).evalf(), dtype=float)
print('Rod の閉じた式と scipy.linalg.expm(-[g]×) の差 = %.2e' % np.abs(Rn - sl.expm(-Kn)).max())
ca = sp.symbols('c s', real=True)
cc_, ss_ = ca
sx, sz, sxn = sp.symbols('sx sz sxn')
print('(d) (15.11): sx\' = c sx + s sz のとき (sz - s sx\')/c = c sz - s sx (c² + s² = 1): ',
      sp.simplify(((sz - ss_ * (cc_ * sx + ss_ * sz)) / cc_ - (cc_ * sz - ss_ * sx)).subs(ss_, sp.sqrt(1 - cc_**2))))

# ---------------------------------------------------------------- (e)
print('\n== (e) スピン反転率から P_ST と 1/τ_ST ==')
w0 = sp.symbols('w0', positive=True)        # = r_e ƛ_C c γ⁵/|ρ|³
wup = sp.Rational(5, 16) * sp.sqrt(3) * w0 * (1 + 8 * sp.sqrt(3) / 15)
wdn = sp.Rational(5, 16) * sp.sqrt(3) * w0 * (1 - 8 * sp.sqrt(3) / 15)
print('w↑↓ + w↓↑ =', sp.simplify(wup + wdn), '(= 5√3/8 w0)')
print('(w↑↓ - w↓↑)/(w↑↓ + w↓↑) =', sp.simplify((wup - wdn) / (wup + wdn)), '=', sp.nsimplify(sp.simplify((wup - wdn) / (wup + wdn))), ' (= 8/(5√3) =', sp.simplify(8 / (5 * sp.sqrt(3))), ')')

# ---------------------------------------------------------------- (f)
print('\n== (f) γ ∂/∂γ と ∂/∂δ ==')
p0, d = sp.symbols('p0 delta', positive=True)
p = p0 * (1 + d)
gam_p = sp.sqrt(1 + p**2)
dgdd = sp.diff(gam_p, d)
print('dγ/dδ - p0 p/γ =', sp.simplify(dgdd - p0 * p / gam_p))
print('γ/(dγ/dδ) - γ²/(p0 p) =', sp.simplify(gam_p / dgdd - gam_p**2 / (p0 * p)), '。δ = 0 で γ²/(p0²) = 1/β0² =', sp.simplify((gam_p**2 / (p0 * p)).subs(d, 0) - (1 + p0**2) / p0**2))

# ---------------------------------------------------------------- (g)
print('\n== (g) D の積算 (15.12): D <- R_i (D + G_i M_{0→i}) と、1 周の回転の有限差分 ==')


def rotm(v):
    th = np.linalg.norm(v)
    if th < 1e-300:
        return np.eye(3)
    k = v / th
    Kk = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) - math.sin(th) * Kk + (1 - math.cos(th)) * Kk @ Kk      # exp(-[v]×)


def logrot(Rm):
    """exp(-[v]×) = Rm の v (|v| < π)"""
    ang = math.acos(max(-1, min(1, (np.trace(Rm) - 1) / 2)))
    a = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]])
    return -ang * a / (2 * math.sin(ang))


rg = np.random.default_rng(11)
Mtot_angle = 1.3                  # 1 周のスピンの回転 (全スライスの |g| の合計の目安)
for Nsl in (4, 8, 16, 32, 64):
    rg = np.random.default_rng(11)
    Ti = []
    gi = []
    Gi = []
    # 全体の「1 周」を Nsl 個のスライスに分ける。各スライスの g_i = (全体)/Nsl、座標の依存 G_i = (全体)/Nsl
    # (スライスを細かくしたとき同じ全体の写像になるように、全体の量を固定して 1/Nsl に割る)
    g_total = np.array([0.4, 1.2, 0.3])
    G_total = np.random.default_rng(5).normal(size=(3, 6)) * 0.5
    Tx = np.eye(6)
    ang = 0.05
    Rs = np.eye(6)
    for k in range(3):
        c_, s_ = math.cos(ang * (k + 1)), math.sin(ang * (k + 1))
        Rs[2 * k:2 * k + 2, 2 * k:2 * k + 2] = [[c_, s_], [-s_, c_]]
    # 各スライスの座標の写像: Rs^(1/Nsl) (回転角 ang/Nsl)
    Ts = np.eye(6)
    for k in range(3):
        c_, s_ = math.cos(ang * (k + 1) / Nsl), math.sin(ang * (k + 1) / Nsl)
        Ts[2 * k:2 * k + 2, 2 * k:2 * k + 2] = [[c_, s_], [-s_, c_]]
    gslice = g_total / Nsl
    Gslice = G_total / Nsl

    def total_rot(x0):
        x = x0.copy()
        Rm = np.eye(3)
        for i in range(Nsl):
            gv_ = gslice + Gslice @ x
            Rm = rotm(gv_) @ Rm
            x = Ts @ x
        return Rm

    # 公式 D
    D = np.zeros((3, 6))
    Rm0 = np.eye(3)
    Mi = np.eye(6)
    for i in range(Nsl):
        Ri = rotm(gslice)
        D = Ri @ (D + Gslice @ Mi)
        Rm0 = Ri @ Rm0
        Mi = Ts @ Mi
    # 有限差分: R(x0) R0^{-1} = exp(-[D x0]×) の v
    err = 0.0
    scale = 1e-6
    for j in range(6):
        x0 = np.zeros(6)
        x0[j] = scale
        v = logrot(total_rot(x0) @ total_rot(np.zeros(6)).T)
        err = max(err, np.abs(v / scale - D[:, j]).max())
    print('スライス数 %3d: |D (15.12 の積算) - 有限差分| の最大 = %.3e   (|D| の最大 = %.3f)' % (Nsl, err, np.abs(D).max()))

# ---------------------------------------------------------------- (h)
print('\n== (h) 逐次の N の式の畳み込み: N_{i+1} T_i = R_i N_i + [n_{i+1}]× D_i を 1 周でつなぐと N M = R0 N + [n0]× D ==')
rg = np.random.default_rng(21)
nel = 6
Ts_, Rs_, Ds_ = [], [], []
for i in range(nel):
    A = rg.normal(size=(6, 6)) * 0.2 + np.eye(6)
    Ts_.append(A)
    Rs_.append(rotm(rg.normal(size=3) * 0.7))
    Ds_.append(rg.normal(size=(3, 6)) * 0.3)
# 1 周: T_tot = T_n ... T_1、R_tot = R_n ... R_1、D_tot = Σ R_n..R_{i+1} D_i T_{i-1}...T_1 (D_i は出口の座標系、R_i の後。ここでは「R_i exp(-[D_i x])」ではなく 要素 i を S -> exp(-[D_i x_in]) R_i S とする)
Ttot = np.eye(6)
Rtot = np.eye(3)
Dtot = np.zeros((3, 6))
for i in range(nel):
    Dtot = Dtot + 0  # 位置合わせ
    Dtot = Dtot * 1.0
    # S -> exp(-[D_i x_in]) R_i S: 次の要素の左に残るので、既存の D は R_i を左からかける (R_i exp(-[D x]) R_i^{-1} = exp(-[R_i D x]))
    Dtot = Rs_[i] @ Dtot + Ds_[i] @ Ttot
    Ttot = Ts_[i] @ Ttot
    Rtot = Rs_[i] @ Rtot
# n0: Rtot の軸
w_, V_ = np.linalg.eig(Rtot)
n0 = np.real(V_[:, np.argmin(abs(w_ - 1))])
n0 /= np.linalg.norm(n0)
cm = lambda v: np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
# 1 周の式 (15.13): N M = R0 N + [n0]× D (ここで 1 要素の式は S -> exp(-[D x]) R S、n(x_out) = R n(x_in) 型)
Kmat = np.eye(18) - np.kron(np.linalg.inv(Ttot).T, Rtot)
rhs = cm(n0) @ Dtot @ np.linalg.inv(Ttot)
Nvec = np.linalg.lstsq(Kmat, rhs.reshape(-1, order='F'), rcond=None)[0]
N0 = Nvec.reshape((3, 6), order='F')
print('(15.13) の残差 |N M - R0 N - [n0]× D| = %.2e (n0·N = %.1e。注: ここで M はシンプレクティックとは限らない任意の行列で、Rtot は固有値 1 を持つ)' % (np.abs(N0 @ Ttot - Rtot @ N0 - cm(n0) @ Dtot).max(), np.abs(n0 @ N0).max()))
# 逐次の式で N_i をつないで、1 周後に N_0 に戻るか
Ni = N0.copy()
ni = n0.copy()
Tcum = np.eye(6)
Dcum = None
for i in range(nel):
    ni1 = Rs_[i] @ ni
    # 要素 i: S -> exp(-[D_i x]) R_i S。軸: n_{i+1}(T x) = exp(-[D_i x]) R_i n_i(x) ≈ n_{i+1} + R_i N x + n_{i+1} × (D_i x)
    Ni1 = (Rs_[i] @ Ni + cm(ni1) @ Ds_[i]) @ np.linalg.inv(Ts_[i])
    ni, Ni = ni1, Ni1
print('逐次の式 N_{i+1} = (R_i N_i + [n_{i+1}]× D_i) T_i^{-1} で 1 周つないだ後の |N_end - N_0| = %.2e、|n_end - n_0| = %.1e' % (np.abs(Ni - N0).max(), np.abs(ni - n0).max()))

# ---------------------------------------------------------------- (i)
print('\n== (i) ガウス分布のスピンチューン広がりを入れた和 (K64 の spintrm の形) ==')
import mpmath as mp
mp.mp.dps = 30
for c_, a_ in ((0.3, 0.5), (1.0, 0.1), (2.0 + 1.5j, 0.4)):
    lhs = mp.quad(lambda t: mp.e**(-c_ * t - a_**2 * t**2), [0, mp.inf])
    rhs_ = mp.sqrt(mp.pi) / (2 * a_) * mp.e**((c_ / (2 * a_))**2) * mp.erfc(c_ / (2 * a_))
    print('c = %s, a = %s: ∫ = %s、(√π/2a) exp(c²/4a²) erfc(c/2a) = %s' % (c_, a_, mp.nstr(lhs, 12), mp.nstr(rhs_, 12)))
print('c, a の記号で: d/dc の微分方程式を満たすことの確認 (sympy):')
cs, as_, ts = sp.symbols('c a t', positive=True)
Iexpr = sp.sqrt(sp.pi) / (2 * as_) * sp.exp(cs**2 / (4 * as_**2)) * sp.erfc(cs / (2 * as_))
# I(c) = ∫ e^{-ct - a²t²} dt の満たす式: dI/dc = -∫ t e^{-ct-a²t²} = -(1 - c I)/(2a²)   (部分積分: ∫ (2a²t + c) e^{...} dt = 1)
lhs_ode = sp.diff(Iexpr, cs)
rhs_ode = -(1 - cs * Iexpr) / (2 * as_**2)
print('  dI/dc + (1 - c I)/(2a²) = ', sp.simplify(lhs_ode - rhs_ode), '、I(0) =', sp.simplify(Iexpr.subs(cs, 0)), '(期待 √π/2a)')
print('\nexp(x²)·erfc(x) の倍精度での計算 (erfc の漸近式 exp(-x²)(1 - 1/(2x²) + ...)/(x√π) と exp(x²) を別々に掛ける):')
for x in (5.0, 18.0, 20.0, 26.0, 26.7, 27.0, 30.0):
    e1 = math.exp(-x * x) if x * x < 745 else 0.0
    try:
        e2 = math.exp(x * x)
    except OverflowError:
        e2 = float('inf')
    ser = 1 - 1 / (2 * x * x) + 3 / (4 * x**4)
    prod = e2 * (e1 * ser / (x * math.sqrt(math.pi)))
    print('  x = %5.1f: exp(x²) = %.3e、漸近式の erfc = %.3e、積 = %s (正確な値 1/(x√π)(1 - 1/2x² + ...) = %.6e)' % (
        x, e2, e1 * ser / (x * math.sqrt(math.pi)), 'NaN' if math.isnan(prod) else '%.6e' % prod, ser / (x * math.sqrt(math.pi))))
print('exp の引数が 709.78 (= ln 1.8e308) を超えると exp が桁あふれ: x > %.2f' % math.sqrt(709.78))

# ---------------------------------------------------------------- (j)
print('\n== (j) K64 の sprot: pzm = √(1 - pxm² - pym²) による b = sinθ/pzm と、弦 2 sin(θ/2) の相対差 (px0 = 0 の平面の偏向) ==')
for th_ in (0.05, 0.1, 0.2617993877991494, 0.5):
    pxm_ = math.sin(th_) / 2
    pzm_ = math.sqrt(1 - pxm_**2)
    bK = math.sin(th_) / pzm_
    chord = 2 * math.sin(th_ / 2)
    print('θ = %.4f: sinθ/pzm = %.8f、2sin(θ/2) = %.8f、相対差 = %.3e、θ⁴/32 = %.3e、(弦 - θ)/θ = %.3e (= -θ²/24 = %.3e)' % (
        th_, bK, chord, bK / chord - 1, th_**4 / 32, chord / th_ - 1, -th_**2 / 24))

# ---------------------------------------------------------------- (k)
print('\n== (k) N の正規座標での形: z = η^T N (η = e2 + i e3) は z M = e^{-iμs} z + i d (d = η^T D)。複素固有モードごとに z_k = i e^{iμs} d_k/(e^{-a_k ± iμ_k + iμs} - 1) ==')
rg = np.random.default_rng(5)
n0v = np.array([0.0, 1.0, 0.0])
mus = 2 * math.pi * 0.403
# n0 = ŷ のまわりの R0 = exp(-μ[n0]×)、e2 = x̂, e3 = n0×e2 = ŷ×x̂ = -ẑ
R0 = rotm(mus * n0v)
e2 = np.array([1.0, 0.0, 0.0])
e3 = np.cross(n0v, e2)
Mm = np.zeros((6, 6))
mu_k = 2 * math.pi * np.array([0.27, 0.19, 0.03])
a_k = np.array([1e-3, 2e-3, 5e-3])
for k in range(3):
    c, s_ = math.cos(mu_k[k]), math.sin(mu_k[k])
    Mm[2 * k:2 * k + 2, 2 * k:2 * k + 2] = math.exp(-a_k[k]) * np.array([[c, s_], [-s_, c]])
Dm = rg.normal(size=(3, 6)) * 0.3
Mi = np.linalg.inv(Mm)
Kk = np.eye(18) - np.kron(Mi.T, R0)
Nn = np.linalg.lstsq(Kk, (cm(n0v) @ Dm @ Mi).reshape(-1, order='F'), rcond=None)[0].reshape((3, 6), order='F')
eta = e2 + 1j * e3
z_direct = eta @ Nn
d_ = eta @ Dm
lam, V = np.linalg.eig(Mm)
Vi = np.linalg.inv(V)
z_mode = 1j * d_ @ V @ np.diag(1 / (lam - np.exp(-1j * mus))) @ Vi
print('z = η^T N (連立方程式の解) と モードごとの公式の差 = %.2e' % np.abs(z_direct - z_mode).max())
print('固有値 λ_k = e^{-a_k ± iμ_k}、分母 |e^{-a_k ± iμ_k + iμs} - 1| = ', np.array2string(np.abs(lam * np.exp(1j * mus) - 1), precision=4))

# ---------------------------------------------------------------- (l)
print('\n== (l) 2 次の相関和と SLIM: 定常の軌道の揺らぎ x_n (<x_{n+m} x_n> = Λ^m Σ) による脱偏極率 ==')
r_, ph_ = sp.symbols('r phi', positive=True)
q_ = r_ * sp.exp(sp.I * ph_)
expr = sp.re(sp.simplify((1 + q_) / (1 - q_))) if False else None
qq = sp.symbols('qr qi', real=True)
q = qq[0] + sp.I * qq[1]
lhs = sp.re(sp.simplify(sp.expand_complex((1 + q) / (1 - q))))
rhs = (1 - (qq[0]**2 + qq[1]**2)) / ((1 - qq[0])**2 + qq[1]**2)
print('Re[(1+q)/(1-q)] - (1 - |q|²)/|1 - q|² =', sp.simplify(lhs - rhs))
print('  (1/2 + Σ_{m>=1} q^m の実部 = (1/2) Re[(1+q)/(1-q)] = (1/2)(1-|q|²)/|1-q|²。q = e^{-a + i(μs ± μk)}、<|x|²> = ε とすると、')
print('   率 = ε |d|² (1/2)(1 - e^{-2a})/|1-q|² = (1/2) B̄ |d|²/|1-q|² (B̄ = ε(1 - e^{-2a}) = 1 周の励起)。(15.14)-(15.15) の SLIM と一致。')
print('   ε を残したまま (1 - e^{-2a}) を落とすと率は 1/(1 - e^{-2a}) ≈ 1/(2a) 倍に過大になる。例: a = 5e-4 なら %.0f 倍' % (1 / (1 - math.exp(-1e-3))))

# ---------------------------------------------------------------- (m)
print('\n== (m) K64 の spintrm の a ≠ 0 の枝 (cm = -mat·2π√π exp(cb² - cφ) erfc(cb)/a、cb = -(e^{-cφ} - 1)/a) の a → 0 の極限 ==')
mp.mp.dps = 40
cphi = mp.mpc(-0.002, 0.3)
cm0 = -1 / (mp.e**cphi - 1)                  # a = 0 の枝: cm = cm0 = -mat/(e^{cφ} - 1) (mat = 1)
for a_ in (1e-2, 1e-3, 1e-4, 1e-5):
    a_ = mp.mpf(a_)
    cb = -(mp.e**(-cphi) - 1) / a_
    cm = -2 * mp.pi * mp.sqrt(mp.pi) * mp.e**(cb**2 - cphi) * mp.erfc(cb) / a_
    print('a = %.0e: |cb| = %.1f、cm/cm0 = %s (2π = %.7f)' % (float(a_), float(abs(cb)), mp.nstr(cm / cm0, 8), 2 * math.pi))
print('→ a → 0 で cm/cm0 → 2π (a = 0 の枝は cm = cm0)。a ≠ 0 の枝の 2π は、a → 0 の極限で a = 0 の枝と合わない因子になっている。')
print('連続の近似 ∫0^∞ exp(-c t - a_g² t²) dt (a_g = a/2、ガウス e^{-a²m²/4} = e^{-(ind σ)² m²/2}) = (√π/a) exp(cb²) erfc(cb) (cb = c/a)')
c_ = mp.mpf('0.05')
a_g = mp.mpf('0.2')
Iint = mp.quad(lambda t: mp.e**(-c_ * t - (a_g / 2)**2 * t**2), [0, mp.inf])
cb_ = c_ / a_g
print('c = 0.05, a = 0.2: ∫ = %s、(√π/a) exp(cb²) erfc(cb) = %s' % (mp.nstr(Iint, 10), mp.nstr(mp.sqrt(mp.pi) / a_g * mp.e**(cb_**2) * mp.erfc(cb_), 10)))
