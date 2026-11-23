# -*- coding: utf-8 -*-
"""四極の端部 (フリンジ) の式 (04_quad.md §7、式 (4.24)〜(4.34)) を確かめる。

(a) 非線形フリンジの母関数 G (ヘルプの QUAD の transformation) から x2, y2, px, py, z を数値で作り (複素ステップ微分)、
    K64 の式 (track6._ttfrin4 = ttfrin 4 次、tqfrie の COD 部分) と一致すること。1 次は H_f = K/(12p)(x^3 px + 3 x y^2 px - 3 x^2 y py - y^3 py) の流れ。
(b) Maxwell 方程式を満たす磁場 (K(s) の 3 次まで) が tanh で立ち上がる端を、ローレンツ力で厳密に積分した結果と、
    「ドリフト → 線形・非線形フリンジ → 一様な四極」(硬い端 + フリンジ) の結果を比べる。
(c) 線形フリンジ F1, F2 (式 (4.31)〜(4.33)): 一般の立ち上がりの K(s) の線形運動を数値積分し、I_n から a, b を作った硬い端の写像と K の 1 次まで比べる。
(d) 両フリンジの 6x6 ヤコビアン (正準運動量) のシンプレクティック性。
"""
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import solve_ivp, quad
import _common as C

C.start(__file__)
from types import SimpleNamespace
from pysad import track6 as t6

J6 = np.zeros((6, 6))
for i in range(3):
    J6[2 * i, 2 * i + 1] = 1.0
    J6[2 * i + 1, 2 * i] = -1.0


# ---------------------------------------------------------------- (a) 母関数
def G(x1, y1, px2, py2, p, K):
    """ヘルプの母関数 G = H0 + (∂H0/∂x1 ∂H0/∂px2 + ∂H0/∂y1 ∂H0/∂py2)/2、H0 = px2 Δx1 + py2 Δy1。
    Δx1 = x1 (a/3 + b) = K (x1^3 + 3 x1 y1^2)/(12 p)、Δy1 = -y1 (a + b/3) = -K (3 x1^2 y1 + y1^3)/(12 p) (a = K x1^2/(4p)、b = K y1^2/(4p)、
    K は長さあたりの K1/L)。偏微分は手で書いた式 (複素ステップ微分の入れ子を避けるため)"""
    c = K / (12 * p)
    dx, dy = c * (x1**3 + 3 * x1 * y1**2), -c * (3 * x1**2 * y1 + y1**3)
    dx_x, dx_y = c * 3 * (x1**2 + y1**2), c * 6 * x1 * y1                   # ∂Δx1/∂x1、∂Δx1/∂y1
    dy_x, dy_y = -c * 6 * x1 * y1, -c * 3 * (x1**2 + y1**2)                 # ∂Δy1/∂x1、∂Δy1/∂y1
    H0_x = px2 * dx_x + py2 * dy_x                                          # ∂H0/∂x1
    H0_y = px2 * dx_y + py2 * dy_y                                          # ∂H0/∂y1
    return px2 * dx + py2 * dy + 0.5 * (H0_x * dx + H0_y * dy)              # ∂H0/∂px2 = Δx1、∂H0/∂py2 = Δy1


def dd(f, k, v, h=1e-30):
    """複素ステップ微分 ∂f/∂v[k] (f は解析的)"""
    v = [complex(a) for a in v]
    v[k] += 1j * h
    return f(*v).imag / h


def map_from_G(c, K):
    """母関数 F = x1 px2 + y1 py2 + G から (x2, px2, y2, py2, z2): p_old = ∂F/∂x1 (px2 について線形なので解く)、x_new = ∂F/∂px2、z2 = z1 + ∂G/∂p。
    c = (x, px, y, py, z, δ) は canonical (px は p0 で規格化した正準運動量)"""
    x1, px1, y1, py1, z1, dl = c
    p = 1.0 + dl
    Gf = lambda X, Y, PX, PY, P: G(X, Y, PX, PY, P, K)

    def oldp(px2, py2):
        v = [x1, y1, px2, py2, p]
        return np.array([px2 + dd(Gf, 0, v), py2 + dd(Gf, 1, v)])         # 型 2 の母関数 F = x1 px2 + y1 py2 + G: p_old = ∂F/∂x1
    f0 = oldp(0.0, 0.0)
    A = np.array([oldp(1.0, 0.0) - f0, oldp(0.0, 1.0) - f0]).T
    px2, py2 = np.linalg.solve(A, np.array([px1, py1]) - f0)
    v = [x1, y1, px2, py2, p]
    return np.array([x1 + dd(Gf, 2, v), px2, y1 + dd(Gf, 3, v), py2, z1 + dd(Gf, 4, v), dl])


def k64_fringe(c, K, sign=+1.0):
    """track6._ttfrin4 (ttfrin、nord = 4)。canonical → 規格化した運動量に直して適用し、戻す。ak = sign K、al = 1 (ak/al が長さあたりの K)"""
    pr = 1.0 + c[5]
    q = SimpleNamespace(c=np.array([c[0], c[1] / pr, c[2], c[3] / pr, c[4], c[5]], float))
    t6._ttfrin4(q, sign * K, 1.0)
    o = q.c
    return np.array([o[0], o[1] * pr, o[2], o[3] * pr, o[4], o[5]])


print('== (a) 母関数 G から作った写像と、K64 の式 (ttfrin 4 次) の比較 (K = 3 [1/m^2]) ==')
rng = np.random.default_rng(3)
K0 = 3.0
worst = 0.0
for _ in range(4):
    r = rng.normal(0, 0.01, 4)
    c = np.array([r[0], r[1], r[2], r[3], 0.0, rng.normal(0, 0.005)])
    worst = max(worst, np.abs(map_from_G(c, K0) - k64_fringe(c, K0)).max())
print('  乱数 4 点 (x, y ~ 10 mm、δ ~ 0.5%%) での最大差 %.2e (x, px, y, py, z)' % worst)
c = np.array([0.01, 0.004, -0.008, 0.003, 0.0, 0.0])
dk = k64_fringe(c, K0) - c
fx = K0 / 12 * (c[0]**3 + 3 * c[0] * c[2]**2)
fy = -K0 / 12 * (3 * c[0]**2 * c[2] + c[2]**3)
print('  x の変化: K64 %.6e、H_f の 1 次の流れ (K/12)(x^3 + 3 x y^2) = %.6e (差 %.1e は t x = (a-b)^2 x/6 の項)' % (dk[0], fx, dk[0] - fx))
print('  y の変化: K64 %.6e、H_f の 1 次の流れ -(K/12)(3 x^2 y + y^3) = %.6e (差 %.1e)' % (dk[2], fy, dk[2] - fy))


# ---------------------------------------------------------------- (b) Maxwell 場の端
def maxwell_rhs(s, v, K0, d, sgn):
    """運動量 p = 1 の粒子 (x, y, px, py, z)。K(s) = K0 (1 + sgn tanh(s/d))/2、sgn = +1 (入口) / -1 (出口)。
    場 (Bρ で規格化): bx = K y - K''(3x^2 y + y^3)/12、by = K x - K''(x^3 + 3 x y^2)/12、bs = K' x y (div B = 0、rot B = 0 を 3 次まで満たす)。
    運動方程式: d(px, py, pz)/ds = (p × b)/pz"""
    x, y, px, py, z = v
    t = math.tanh(s / d)
    sech2 = 1.0 - t * t
    K = K0 * (1 + sgn * t) / 2
    K1 = sgn * K0 * sech2 / (2 * d)
    K2 = -sgn * K0 * t * sech2 / d**2
    bx = K * y - K2 * (3 * x * x * y + y**3) / 12
    by = K * x - K2 * (x**3 + 3 * x * y * y) / 12
    bs = K1 * x * y
    pz = math.sqrt(max(1.0 - px * px - py * py, 1e-12))
    return [px / pz, py / pz, (py * bs - pz * by) / pz, (pz * bx - px * bs) / pz, 1.0 - 1.0 / pz]


def rhs_lin(s, v, K0, d, sgn):
    """同じ K(s) で、線形な場 bx = K y、by = K x だけ (K''、K' の項なし。Maxwell 方程式は満たさない) の運動"""
    x, y, px, py, z = v
    K = K0 * (1 + sgn * math.tanh(s / d)) / 2
    pz = math.sqrt(max(1.0 - px * px - py * py, 1e-12))
    return [px / pz, py / pz, -K * x, K * y, 1.0 - 1.0 / pz]


def kin_body(c, S, K):
    """一様な四極 (長さ S、勾配 K) の運動学的な厳密な写像 (DOP853)。c = (x, y, px, py, z)"""
    def rhs(s, v):
        x, y, px, py, z = v
        pz = math.sqrt(1.0 - px * px - py * py)
        return [px / pz, py / pz, -K * x, K * y, 1.0 - 1.0 / pz]
    return solve_ivp(rhs, (0, S), list(c), method='DOP853', rtol=1e-13, atol=1e-16).y[:, -1]


def kin_drift(c, S):
    x, y, px, py, z = c
    pz = math.sqrt(1.0 - px * px - py * py)
    return np.array([x + px / pz * S, y + py / pz * S, px, py, z + (1.0 - 1.0 / pz) * S])


def lin_fringe(c, a, b):
    """track6._tqlfr (p = 1)。c = (x, y, px, py, z)。a, b は f1, f2 (式 (4.33) の a, b)"""
    q = SimpleNamespace(c=np.array([c[0], c[2], c[1], c[3], c[4], 0.0]))
    t6._tqlfr(q, a, b)
    return np.array([q.c[0], q.c[2], q.c[1], q.c[3], q.c[4]])


def nonlin_fringe(c, K, sign):
    o = k64_fringe(np.array([c[0], c[2], c[1], c[3], c[4], 0.0]), K, sign)
    return np.array([o[0], o[2], o[1], o[3], o[4]])


print('\n== (b) Maxwell 場の端 (tanh の立ち上がり、幅 d) を厳密に積分した結果との比較 (K = 3 [1/m^2]、入口、S = 0.06 m) ==')
print('   線形な場 (K(s) x, K(s) y のみ) の ODE は、線形フリンジの確認になる。Maxwell の場 (K\'\', K\' の項つき) との差が非線形フリンジの効果')
print('   Δ(ODE) = ODE(Maxwell) - ODE(線形な場)、Δ(モデル) = (非線形フリンジあり) - (なし)。r は x, y の大きさ')
S = 0.06
print('   r [mm]  d [mm] | 線形: フリンジなし  F1 あり | 非線形: Δ(ODE)   Δ(モデル)   差/Δ(モデル)')
rowsb = {}
for scale in (0.2, 0.1):
    for d in (0.0005, 0.001, 0.002, 0.005):
        c0 = np.array([0.010, -0.008, 0.005, -0.003, 0.0]) * scale
        F1 = math.pi * d                                      # tanh の端: I0 = I2 = 0、I1 = -π² d²/24 → F1 = π d、F2 = 0
        f1 = -K0 * F1 * F1 / 24
        kw = dict(method='DOP853', rtol=1e-13, atol=1e-17)
        ex_m = solve_ivp(maxwell_rhs, (-S, S), c0, args=(K0, d, +1.0), **kw).y[:, -1]
        ex_l = solve_ivp(rhs_lin, (-S, S), c0, args=(K0, d, +1.0), **kw).y[:, -1]
        cc = kin_drift(c0, S)
        m_none = kin_body(cc, S, K0)
        m_lin = kin_body(lin_fringe(cc, f1, 0.0), S, K0)
        m_all = kin_body(nonlin_fringe(lin_fringe(cc, f1, 0.0), K0, +1.0), S, K0)
        dex, dmod = ex_m - ex_l, m_all - m_lin
        rel = np.abs(dex - dmod).max() / np.abs(dmod).max()
        rowsb[(scale, d)] = rel
        print('   %5.1f   %5.1f  |   %.2e      %.2e   |  %.3e   %.3e     %.4f' % (
            scale * 10, d * 1e3, np.abs(ex_l - m_none).max(), np.abs(ex_l - m_lin).max(), np.abs(dex).max(), np.abs(dmod).max(), rel))
print('   (r を半分にしても「差/Δ(モデル)」は変わらない = 3 次の効果は合っている。残りは d に比例し、K に比例 (約 2 K d [m] の程度)。')
print('    原因は未確認。硬い端の極限 d -> 0 で消え、振幅に依らないので、K の 2 次の補正 (フリンジの式は K の 1 次まで) と推定)')

print('\n   出口 (K が tanh で下がる端): 一様な四極(S) → 非線形フリンジ (出口) → 線形フリンジ (出口) → ドリフト(S)')
for d in (0.001, 0.002):
    c0 = np.array([0.010, -0.008, 0.005, -0.003, 0.0]) * 0.1
    kw = dict(method='DOP853', rtol=1e-13, atol=1e-17)
    ex_m = solve_ivp(maxwell_rhs, (-S, S), c0, args=(K0, d, -1.0), **kw).y[:, -1]
    ex_l = solve_ivp(rhs_lin, (-S, S), c0, args=(K0, d, -1.0), **kw).y[:, -1]
    f1 = -K0 * (math.pi * d)**2 / 24
    cc = kin_body(c0, S, K0)
    m_lin = kin_drift(lin_fringe(cc, -f1, 0.0), S)            # 出口: tqlfri(-f1out、f2out) (a = -f1out/p)
    m_all = kin_drift(lin_fringe(nonlin_fringe(cc, K0, -1.0), -f1, 0.0), S)     # 出口: ttfrin(-ak) のあと線形フリンジ
    dex, dmod = ex_m - ex_l, m_all - m_lin
    print('   d = %.0f mm: 線形 (F1 あり) の残り %.2e、非線形 Δ(ODE) %.3e、Δ(モデル) %.3e、差/Δ(モデル) %.4f' % (
        d * 1e3, np.abs(ex_l - m_lin).max(), np.abs(dex).max(), np.abs(dmod).max(), np.abs(dex - dmod).max() / np.abs(dmod).max()))

# ---------------------------------------------------------------- (c) 線形フリンジ
print('\n== (c) 線形フリンジ: 立ち上がり K(s) = K0 (s/ℓ)^2 (0 <= s <= ℓ、ℓ = 0.1 m)。K0 を小さくして K の 1 次を見る ==')
ell = 0.1
s0 = 2 * ell / 3                                           # I0 = 0 になる実効的な端の位置


def prof(s):
    return 0.0 if s < 0 else (1.0 if s > ell else (s / ell)**2)


def Iint(n):
    return quad(lambda s: (s - s0)**n * (prof(s) - (1.0 if s > s0 else 0.0)), 0, ell, points=[s0], epsabs=1e-16, epsrel=1e-13)[0]


I0, I1, I2 = Iint(0), Iint(1), Iint(2)
aa = 24 * (I0**2 / 2 - I1)
F1 = math.copysign(math.sqrt(abs(aa)), aa)
F2 = I2 - I0**3 / 3
print('  s0 = %.6f、I0 = %.2e、I1 = %.6e (-ℓ^2/36 = %.6e)、I2 = %.6e (ℓ^3/405 = %.6e)' % (s0, I0, I1, -ell**2 / 36, I2, ell**3 / 405))
print('  F1 = %.6f、F2 = %.6e' % (F1, F2))
S1 = 0.2
print('   K0       フリンジなしの差   a, b を入れた差   (x-px の 2x2 行列の最大差)')
prev = None
for K in (1e-3, 2e-3, 4e-3):
    def lin_rhs(s, v):
        k = K * prof(s)
        M = v.reshape(2, 2)
        return np.array([[M[1, 0], M[1, 1]], [-k * M[0, 0], -k * M[0, 1]]]).ravel()
    Mx = solve_ivp(lin_rhs, (-S1, ell + S1), np.eye(2).ravel(), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1].reshape(2, 2)
    a_map, b_map = -K * F1 * abs(F1) / 24, K * F2
    drift = lambda L_: np.array([[1, L_], [0, 1]])
    w = math.sqrt(K)
    body = lambda L_: np.array([[math.cos(w * L_), math.sin(w * L_) / w], [-w * math.sin(w * L_), math.cos(w * L_)]])
    fr = lambda a, b: np.array([[math.exp(a), b], [0, math.exp(-a)]])
    Lin, Lout = S1 + s0, S1 + ell - s0
    m0 = body(Lout) @ drift(Lin)
    m1 = body(Lout) @ fr(a_map, b_map) @ drift(Lin)
    e0, e1 = np.abs(Mx - m0).max(), np.abs(Mx - m1).max()
    print('   %.0e      %.3e         %.3e   (a = %.3e、b = %.3e)%s' % (K, e0, e1, a_map, b_map, '' if prev is None else '  比 %.2f' % (e1 / prev)))
    prev = e1
print('   (a, b を入れた差は K0 の 2 乗に比例して小さい (K0 が 2 倍で比 4 付近): フリンジの式は K の 1 次まで)')

# ---------------------------------------------------------------- (d) シンプレクティック性
print('\n== (d) フリンジ写像の 6x6 ヤコビアンのシンプレクティック性 (正準運動量、中心差分) ==')
c = np.array([0.01, 0.004, -0.008, 0.003, 0.0, 0.004])


def jac(f, c, h=1e-6):
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (f(cp) - f(cm)) / (2 * h)
    return M


def lin_f(c_, a1=-3.0 * 0.05**2 / 24, b1=3.0 * 0.01):
    pr = 1.0 + c_[5]
    q = SimpleNamespace(c=np.array([c_[0], c_[1] / pr, c_[2], c_[3] / pr, c_[4], c_[5]], float))
    t6._tqlfr(q, a1, b1)
    o = q.c
    return np.array([o[0], o[1] * pr, o[2], o[3] * pr, o[4], o[5]])


Jn = jac(lambda v: k64_fringe(v, K0), c)
Jl = jac(lin_f, c)
print('  非線形フリンジ: |M^T J M - J| = %.1e、線形フリンジ: %.1e' % (np.abs(Jn.T @ J6 @ Jn - J6).max(), np.abs(Jl.T @ J6 @ Jl - J6).max()))
