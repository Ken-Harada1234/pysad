# -*- coding: utf-8 -*-
"""スピン (15_spin.md): Thomas-BMT と 1 ステップの回転、定数、SLIM の N、リングの偏極時間・スピンチューン。

(a) sympy: Thomas 歳差 + 静止系の歳差 = BMT (E = 0)。B の縦・横の係数 (1+a)/γ、(1+aγ)/γ。
(b) 一様磁場のらせん軌道 (縦成分あり) で、BMT の厳密な積分 (DOP853) と、K64 の 1 ステップ回転
    Ω_step = (1+aγ) θ_vec + (1+a) ψ p̂_mid の積 (刻み数 N) を比べる。N を 2 倍にしたときの誤差の比 (次数)。
(c) 定数の恒等式: cave·anrad·cuc² = r_e ƛ_C、(sflc·anrad·cuc の積) → Sokolov-Ternov の係数 cfpd、DK の 11/18。
(d) 玩具の 1 周写像で、N M = R0 N + [n0]× D (SLIM の式) が R(x) n(x) = n(Mx) を 1 次で満たす。
(e) 平面リング (3 GeV、12 セル、偏向の合計 π): pysad の Emittance[] の SpinTune、PolarizationTime、PolarizationVector、
    Polarization (v2) を、aγθ/2π・Sokolov-Ternov の式・8/(5√3) と比べる。
"""
import io
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())

mc2 = 0.51099895e6
re = 2.8179403262e-15
alpha = 1 / 137.035999084
lamC = re / alpha
cc = 299792458.0
a_g = 0.001159652181280002

# ---------------------------------------------------------------- (a)
print('== (a) Thomas-BMT: 静止系の歳差 + Thomas 歳差 = BMT (E = 0, sympy) ==')
import sympy as sp
a, gam, q, m = sp.symbols('a gamma q m', positive=True)
Bx, By, Bz, b1, b2, b3 = sp.symbols('Bx By Bz b1 b2 b3', real=True)
B = sp.Matrix([Bx, By, Bz])
bet = sp.Matrix([b1, b2, b3])
b2s = bet.dot(bet)
# 静止系の磁場 B' = γ [B - γ/(γ+1) β (β·B)] (E = 0)。時間は実験室系 (dt' = dt/γ)
Omega_rest = -((1 + a) * q / m) * (B - gam / (gam + 1) * bet * bet.dot(B))
# 加速度 a = (q/(γ m)) v × B (c = 1)。Thomas: ω_T = γ²/(γ+1) (a × v)
acc = (q / (gam * m)) * bet.cross(B)
omega_T = gam**2 / (gam + 1) * acc.cross(bet)
Omega_sum = Omega_rest + omega_T
Omega_bmt = -(q / m) * ((a + 1 / gam) * B - a * gam / (gam + 1) * bet * bet.dot(B))
diff = (Omega_sum - Omega_bmt).subs(gam, 1 / sp.sqrt(1 - b2s))
print('Ω_rest + ω_T - Ω_BMT (β² = b1²+b2²+b3² として γ を消去) の各成分:', [sp.simplify(d) for d in diff])
# B ∥ β̂ と B ⊥ β̂ の係数
bb = sp.symbols('beta', positive=True)
Bpar, Bperp = sp.symbols('Bpar Bperp', real=True)
bv = sp.Matrix([0, 0, bb])
for name, Bv in (('B∥', sp.Matrix([0, 0, Bpar])), ('B⊥', sp.Matrix([Bperp, 0, 0]))):
    Om = -(q / m) * ((a + 1 / gam) * Bv - a * gam / (gam + 1) * bv * bv.dot(Bv) / 1)
    Om = Om.subs(bb, sp.sqrt(1 - 1 / gam**2))
    Om = sp.simplify(Om / (-(q / m)))
    print('%s に対する Ω/(-q/m) = %s' % (name, list(Om)))
print('→ 縦: (1+a)/γ、横: (1+aγ)/γ (期待: (1+a)/γ, (1+aγ)/γ)')
print('   差 (縦): %s、差 (横): %s' % (
    sp.simplify(sp.simplify(((a + 1 / gam) * Bpar - a * gam / (gam + 1) * (1 - 1 / gam**2) * Bpar)) - (1 + a) / gam * Bpar),
    sp.simplify((a + 1 / gam) * Bperp - (1 + a * gam) / gam * Bperp)))

# ---------------------------------------------------------------- (b)
print('\n== (b) らせん軌道 (一様磁場、B と運動量の角 α) での BMT の厳密な積分と K64 の 1 ステップ回転の積 ==')
from scipy.integrate import solve_ivp


def rot_matrix(axis_vec):
    """回転ベクトル (右手系、ベクトル = 軸 × 角) の回転行列 (Rodrigues)"""
    th = np.linalg.norm(axis_vec)
    if th == 0:
        return np.eye(3)
    k = axis_vec / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(th) * K + (1 - math.cos(th)) * K @ K


gamma = 3000e6 / mc2 + 1.0
gamma = math.sqrt(1 + (3000e6 / mc2) ** 2)
ag = a_g * gamma
Bhat = np.array([0.0, 1.0, 0.0])
alpha_t = 0.3                     # p̂ と B̂ の角 [rad]
p0hat = np.array([math.sin(alpha_t), 0.0, 0.0])
# p̂(0) を B̂ とα の角にする: B̂ = ŷ とし、p̂(0) = (sin α, cos α, 0)ではなく B̂ から測る
p0hat = np.array([math.sin(alpha_t), math.cos(alpha_t), 0.0])
kappa = 1.0                        # 運動量の回転率 q B/p [rad/m]
Ltot = 2.0                         # 経路長 [m]。p̂ は B̂ のまわりに κ·L = 2 rad 回る
S0 = np.array([0.3, 0.2, math.sqrt(1 - 0.09 - 0.04)])


def phat(s):
    return rot_matrix(kappa * s * Bhat) @ p0hat


def rhs(s, S):
    p = phat(s)
    cosa = Bhat @ p
    Om = kappa * ((1 + ag) * (Bhat - cosa * p) + (1 + a_g) * cosa * p)   # Ω/v [rad/m]
    return np.cross(Om, S)


sol = solve_ivp(rhs, [0, Ltot], S0, method='DOP853', rtol=1e-13, atol=1e-14)
S_ref = sol.y[:, -1]
print('γ = %.2f、aγ = %.5f、α = %.2f rad、κL = %.2f rad。厳密な積分 (DOP853): S = %s' % (gamma, ag, alpha_t, kappa * Ltot, S_ref))
prev = None
for N in (10, 20, 40, 80, 160):
    ds = Ltot / N
    S = S0.copy()
    for i in range(N):
        pin = phat(i * ds)
        pout = phat((i + 1) * ds)
        pm = pin + pout
        pm = pm / np.linalg.norm(pm)
        cr = np.cross(pin, pout)
        sinth = np.linalg.norm(cr)
        theta = math.asin(sinth)
        th_vec = theta * cr / sinth                       # 運動量の向きの回転ベクトル (横)
        psi = kappa * (Bhat @ pm) * ds                    # ∫B∥ ds / (Bρ)
        Om_step = (1 + ag) * th_vec + (1 + a_g) * psi * pm
        S = rot_matrix(Om_step) @ S
    err = np.linalg.norm(S - S_ref)
    print('N = %4d: |S - S_ref| = %.3e%s' % (N, err, '' if prev is None else '   (前の N の誤差との比 %.2f)' % (prev / err)))
    prev = err

# ---------------------------------------------------------------- (b2) 平面の偏向と縦磁場の単純な場合
print('\n-- 単純な場合: 平面の偏向 (B ⊥ p): S と p のなす角の変化 = aγθ。縦磁場 (B ∥ p): S の回転 = (1+a)ψ --')
for th in (0.2617993877991494, 1.0):
    # 厳密な積分 (α = π/2): p̂ は B̂ のまわりに θ 回る
    pp = np.array([0.0, 0.0, 1.0])
    p0hat_s = np.array([0.0, 0.0, 1.0])
    Bh = np.array([0.0, 1.0, 0.0])
    S_in = np.array([1.0, 0.0, 0.0])     # 水平 (p に垂直)

    def rhs2(s, S):
        p = rot_matrix(s * Bh) @ p0hat_s
        cosa = Bh @ p
        Om = (1 + ag) * (Bh - cosa * p) + (1 + a_g) * cosa * p
        return np.cross(Om, S)
    sl = solve_ivp(rhs2, [0, th], S_in, method='DOP853', rtol=1e-13, atol=1e-14)
    Sf = sl.y[:, -1]
    pf = rot_matrix(th * Bh) @ p0hat_s
    # p̂ の回転と同じ向きの回転角を、p̂ の回転後の座標軸 (x̂' = 回した x̂) から測る
    xh = rot_matrix(th * Bh) @ S_in
    zh = pf
    ang = math.atan2(-(Sf @ zh), Sf @ xh)      # ŷ のまわりの右手系の角 (z → x が正)
    d = (ang - ag * th + math.pi) % (2 * math.pi) - math.pi
    print('θ = %.4f: S の p̂ に対する角 = %.6f rad (2π の整数倍を除く)、aγθ = %.6f (差 mod 2π = %.1e)' % (th, ang, ag * th, d))
psi = 0.5
pz = np.array([0.0, 0.0, 1.0])
sl = solve_ivp(lambda s, S: np.cross((1 + a_g) * kappa * pz, S), [0, psi / kappa], np.array([1.0, 0, 0]), method='DOP853', rtol=1e-13, atol=1e-14)
print('縦磁場 ψ = %.2f: S の回転角 = %.8f、(1+a)ψ = %.8f' % (psi, math.atan2(sl.y[1, -1], sl.y[0, -1]), (1 + a_g) * psi))

print('\n-- K64 の刻みの角: |ppx,ppy,ppz|/pzm = 2 sin(θ/2) (弦)。偏向角 θ を n 刻みにしたときの歳差の角 (1+aγ)·n·2 sin(θ/2n) と (1+aγ)θ の差 --')
for th_, n_ in ((0.2617993877991494, 1), (0.2617993877991494, 3), (0.2617993877991494, 10)):
    chord = n_ * 2 * math.sin(th_ / (2 * n_))
    print('θ = %.4f、%2d 刻み: 角の差 = %.3e rad (相対 %.3e、θ²/(24 n²) = %.3e)' % (
        th_, n_, (1 + ag) * (chord - th_), (chord - th_) / th_, -(th_ / n_) ** 2 / 24))

# ---------------------------------------------------------------- (c)
print('\n== (c) 定数の恒等式 (β = 1 の極限。anrad = 5α/(2√3)、cuc = 1.5 ƛ_C、cave = 8/(15√3)、cuu = 11/27) ==')
anrad = 5 * alpha / (2 * math.sqrt(3))
cuc = 1.5 * lamC
cave = 8 / (15 * math.sqrt(3))
cuu = 11 / 27
cfpd = 5 * math.sqrt(3) / 8
sflc = 0.75 * (re / alpha) ** 2
print('cave·anrad·cuc² = %.6e m²、r_e ƛ_C = %.6e m² (比 %.12f)' % (cave * anrad * cuc**2, re * lamC, cave * anrad * cuc**2 / (re * lamC)))
print('sflc·anrad = %.6e、cfpd·r_e ƛ_C = %.6e (比 %.12f)  [gintd の係数 sflc·anp·(...)² = Sokolov-Ternov の 1/τ の係数]' % (
    sflc * anrad, cfpd * re * lamC, sflc * anrad / (cfpd * re * lamC)))
print('(1/2)·cuu·anrad·cuc² / (sflc·anrad) = %.9f、DK の 11/18 = %.9f   [SLIM の δ の寄与 (1/2) B_δδ と DK の (11/18) の一致]' % (
    0.5 * cuu * anrad * cuc**2 / (sflc * anrad), 11 / 18))
print('pst = 8/(5√3) = %.9f、cfpd = 5√3/8 = %.9f、積 = %.9f' % (8 / (5 * math.sqrt(3)), cfpd, 8 / (5 * math.sqrt(3)) * cfpd))

# ---------------------------------------------------------------- (d)
print('\n== (d) 玩具の 1 周写像で N M = R0 N + [n0]× D を確かめる ==')
rng = np.random.default_rng(7)


def cross_mat(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


n0 = np.array([0.0, 1.0, 0.0])
mu_s = 2 * math.pi * 0.403
R0 = rot_matrix(mu_s * n0)
M = np.zeros((6, 6))
for k, nu in enumerate((0.27, 0.19, 0.03)):
    c, s_ = math.cos(2 * math.pi * nu), math.sin(2 * math.pi * nu)
    M[2 * k:2 * k + 2, 2 * k:2 * k + 2] = 0.9 * np.array([[c, s_], [-s_, c]])
Dm = rng.normal(size=(3, 6)) * 0.3          # 1 周の回転 R(x) = exp(-[D x]×) R0 (K64 の約束: Rot(g) = exp(-[g]×))
# N: N M = R0 N + [n0]× D → vec 方程式 (R0 ⊗ ... )。n0 ⊥ の成分 (2 次元) だけ解く
rhs_mat = cross_mat(n0) @ Dm @ np.linalg.inv(M)
# N - R0 N M^{-1} = [n0]× D M^{-1}: vec(N) (列優先) = (I - M^{-T} ⊗ R0)^{-1} vec(rhs)
Minv = np.linalg.inv(M)
Kmat = np.eye(18) - np.kron(Minv.T, R0)
Nvec = np.linalg.lstsq(Kmat, rhs_mat.reshape(-1, order='F'), rcond=None)[0]
N = Nvec.reshape((3, 6), order='F')
print('n0·N の大きさ (0 のはず): %.2e' % np.abs(n0 @ N).max())
for eps in (1e-2, 1e-3, 1e-4):
    x0 = eps * rng.normal(size=6)
    n_x0 = n0 + N @ x0
    n_x0 /= np.linalg.norm(n_x0)
    # R(x0) は 1 次まで exp(-[D x0]×) R0
    lhs = rot_matrix(-(Dm @ x0)) @ R0 @ n_x0
    n_x1 = n0 + N @ (M @ x0)
    n_x1 /= np.linalg.norm(n_x1)
    print('|x0| = %.0e: |R(x0) n(x0) - n(M x0)| = %.3e   (|N x0| = %.2e)' % (eps, np.linalg.norm(lhs - n_x1), np.linalg.norm(N @ x0)))
# 回転の向きの注意: rot_matrix(-v) = exp(-[v]×) = K64 の Rot(g)。

# ---------------------------------------------------------------- (e)
print('\n== (e) 平面リング (3 GeV、12 セル、偏向 15 度 × 12 = π) の Emittance[] ==')
from pysad.ffscmd import run_sad_string
src = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW; RAD; NOFLUC;
POL;
em = Emittance[];
Print["OUT ", {PolarizationTime, SpinTune, NominalSpinTune, Polarization, PolarizationVector}/.em];
Exit;
'''
buf = io.StringIO()
run_sad_string(src, out=buf)
for line in buf.getvalue().splitlines():
    if line.startswith('OUT'):
        print(line)
        txt = line[3:].strip()
        import re as _re
        nums = _re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|NaN|nan', txt)
        tau, nus, nus0, Pol = float(nums[0]), float(nums[1]), float(nums[2]), nums[3]
thetaB = 12 * 0.2617993877991494
rho = 1.2 / 0.2617993877991494
C0 = 12 * (0.3 + 0.5 + 1.2 + 0.5 + 0.3 + 1.0)
g0 = math.sqrt(1 + (3000e6 / mc2) ** 2)
print('周長 = %.1f m、偏向の合計 = %.5f rad (= π)、ρ = %.5f m、γ = %.3f、aγ = %.6f' % (C0, thetaB, rho, g0, a_g * g0))
print('aγ·θ/2π = %.6f → 端数 %.6f。Emittance[]: SpinTune = %.8f (NominalSpinTune = h0·g = %.8f)' % (
    a_g * g0 * thetaB / (2 * math.pi), (a_g * g0 * thetaB / (2 * math.pi)) % 1, nus, nus0))
rho3 = (12 * 1.2 / C0) / rho**3
rate = cfpd * re * lamC * cc * g0**5 * rho3
print('Sokolov-Ternov: 1/τ = (5√3/8) r_e ƛ_C c γ⁵ <ρ⁻³> = %.6e /s → τ = %.4f s。Emittance[]: PolarizationTime = %.4f s (比 %.6f)' % (
    rate, 1 / rate, tau, tau * rate))
print('Polarization = %s (v2)、8/(5√3) = %.6f' % (Pol, 8 / (5 * math.sqrt(3))))
