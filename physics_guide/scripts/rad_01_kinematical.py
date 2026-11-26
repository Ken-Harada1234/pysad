# -*- coding: utf-8 -*-
"""放射 (13_radiation.md): kinematical method の定数と、1 要素の平均・分散を、古典の式と比べる。

(a) 定数: cuc = 1.5 r_e/rcratio = 1.5 λ_C (換算コンプトン波長)、anrad = 5/(2√3) α (rcratio = α)。
(b) シンクロトロン放射の関数 F(ξ) = ξ ∫_ξ^∞ K_{5/3} の積分から、平均光子エネルギー cave = 8/(15√3)、2 次モーメント cuu = 11/27
    (u_c で規格化) と、光子数・パワーの比。
(c) K64 の乱数 drndsr (有理式の近似) の 1 次・2 次モーメントを、乱数 2e6 個の標本で。
(d) 偏向磁石 1 台 (3 GeV、ANGLE = 2π/16、L = 1.2 m) の追跡 (RAD ON): 平均の損失 = (2/3) r_e mc² γ⁴ θ/ρ、FLUC ON の標準偏差 = √(N_γ cuu) u_c。
(e) 軸からずれた四極 (曲率 = K1 x): 同じ式を、四極内の軌道 x(s) = x0 cos(√k s) で積分した値と、追跡の損失の比較 (どんな要素でも同じ式で扱える)。
(f) 1 スライスの放射による横運動量・位置の補正 -(1/2) dpx dg、x += (-(1/2) dpx dg) al/3 を、連続的な放射の積分から導く (数値積分)。
(g) リング全体: U0 = C_γ E⁴ I2/(2π) と、16 台の偏向磁石の追跡の損失の和。pysad の放射積分 I2。
"""
import io
import math
import numpy as np
from scipy import integrate, special
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6

print('pysad:', C.pysad_src())
mc2 = 0.51099895e6
re = 2.8179403262e-15
alpha = 1 / 137.035999084
lam = re / alpha                              # 換算コンプトン波長 ħ/(mc)
print('== (a) 定数 ==')
print('rcratio = r_e/λ_C = α = %.9f (1/%.4f)、λ_C = %.5e m、cuc = 1.5 λ_C = %.5e m、anrad = 5/(2√3) α = %.6e' % (
    re / lam, 1 / (re / lam), lam, 1.5 * lam, 5 / (2 * math.sqrt(3)) * alpha))

print('\n== (b) 放射スペクトルのモーメント ==')


def F(x):
    return x * integrate.quad(lambda t: special.kv(5.0 / 3.0, t), x, np.inf)[0]


m0 = integrate.quad(lambda x: F(x) / x, 0, np.inf, limit=200)[0]      # 光子数に比例
m1 = integrate.quad(lambda x: F(x), 0, np.inf, limit=200)[0]         # パワーに比例
m2 = integrate.quad(lambda x: x * F(x), 0, np.inf, limit=200)[0]
print('∫F/ξ dξ = %.8f (5π/3 = %.8f)、∫F dξ = %.8f (8π/(9√3) = %.8f)、∫ξF dξ = %.8f (55π/81 = %.8f)' % (
    m0, 5 * math.pi / 3, m1, 8 * math.pi / (9 * math.sqrt(3)), m2, 55 * math.pi / 81))
print('<u>/u_c = ∫F/∫(F/ξ) = %.8f (cave = 8/(15√3) = %.8f)、<u²>/u_c² = ∫ξF/∫(F/ξ) = %.8f (cuu = 11/27 = %.8f)' % (
    m1 / m0, 8 / (15 * math.sqrt(3)), m2 / m0, 11 / 27))


def drndsr(p):
    p = np.asarray(p)
    out = np.empty_like(p)
    a = p < 0.658
    p2 = p[a] ** 2
    out[a] = (((0.41841680 * p2 + 0.14175015) * p2 + 0.30528148) * p2 + 0.53520052) * p2 * p[a]
    q = 1 - p
    b = (~a) & (q > 0.0297)
    qb = q[b]
    out[b] = ((-0.32809490 * qb + 0.20653610) * qb + 0.011920005) / (((qb + 0.88765651) * qb + 0.19270658) * qb + 0.0033141692)
    c = (~a) & (~b)
    t = -np.log(q[c])
    out[c] = (674.99117 * t + 148.32871) / ((t - 225.45670) * t - 692.23986) + t
    return out


print('\n== (c) drndsr (有理式の近似) のモーメント (乱数 2e6 個) ==')
rng = np.random.default_rng(5)
y = drndsr(rng.random(2_000_000))
print('<y> = %.6f (8/(15√3) x 1.0000136 = %.6f)、<y²> = %.6f (11/27 x 1.0000153 = %.6f)。標本の誤差 ~ %.1e' % (
    y.mean(), 8 / (15 * math.sqrt(3)) * 1.0000136, (y * y).mean(), 11 / 27 * 1.0000153, y.std() / math.sqrt(len(y))))

print('\n== (d) 偏向磁石 1 台の放射 (3 GeV) ==')
th = 2 * math.pi / 16
src = 'MOMENTUM = 3GEV;\nBEND B=(L=1.2 ANGLE=%.15g);\nMARK ST=();\nLINE ONE=(ST B);\nFFS;\nUSE ONE;\n' % th
ffs, _ = run_sad_string(src, out=io.StringIO())
els = ffs.lat.elements
a0 = np.zeros((20000, 6))
r1 = np.asarray(track6.track_particles(els, a0, 1, ffs.momentum(), 0.0, False, True, False))
r2 = np.asarray(track6.track_particles(els, a0.copy(), 1, ffs.momentum(), 0.0, False, True, True, seed=7))
g = 3e9 / mc2
rho = 1.2 / th
U = 2 / 3 * re * mc2 * g**4 * th / rho
uc = 1.5 * lam * g**3 / rho * mc2 / 3e9
N = 5 / (2 * math.sqrt(3)) * alpha * g * th
print('θ = %.5f rad、ρ = %.4f m、γ = %.1f、N_γ = anrad γ θ = %.3f 個、u_c/(p0 c) = %.4e (臨界エネルギー %.2f keV)' % (th, rho, g, N, uc, uc * 3e9 / 1e3))
print('平均の損失 δ: 追跡 (FLUC OFF) %.6e、古典 (2/3) r_e mc² γ⁴ θ/ρ / (p0 c) = %.6e、<u> N = cave N u_c = %.6e' % (
    r1[:, 5].mean(), U / 3e9, N * 8 / (15 * math.sqrt(3)) * uc))
print('FLUC ON: 平均 %.6e、標準偏差 %.4e (理論 √(N cuu) u_c = %.4e。標本誤差 ~ %.1e)' % (
    r2[:, 5].mean(), r2[:, 5].std(), math.sqrt(N * 11 / 27) * uc, math.sqrt(N * 11 / 27) * uc / math.sqrt(2 * 20000)))

print('\n== (e) 軸からずれた四極 (曲率 = K1 x) ==')
Lq, Kq = 0.3, 0.35
kq = Kq / Lq
srcq = 'MOMENTUM = 3GEV;\nQUAD Q=(L=%g K1=%g);\nMARK ST=();\nLINE ONE=(ST Q);\nFFS;\nUSE ONE;\n' % (Lq, Kq)
ffq, _ = run_sad_string(srcq, out=io.StringIO())
for x0 in (1e-3, 2e-3, 4e-3):
    a = np.zeros((1, 6))
    a[0, 0] = x0
    r = np.asarray(track6.track_particles(ffq.lat.elements, a, 1, ffq.momentum(), 0.0, False, True, False))
    sk = math.sqrt(kq)
    I2 = kq**2 * x0**2 * (Lq / 2 + math.sin(2 * sk * Lq) / (4 * sk))        # ∫ k² x0² cos²(√k s) ds
    pred = 2 / 3 * re * mc2 * g**4 * I2 / 3e9
    print('x0 = %.0e: 追跡の損失 -δ = %.6e、古典 (2/3) r_e mc² γ⁴ ∫(1/ρ)² ds / (p0 c) = %.6e、比 %.5f' % (x0, -r[0, 5], pred, -r[0, 5] / pred))

print('\n== (f) 1 スライスの放射による横運動量・位置の補正 ==')
# 一定の力 F (単位長さあたり)、一定の放射 (g(s) = dg s/L)。px(s) = q(s) p(s) の ODE: d px/ds = F + (dg/L) px/(1+g) ... を q = px/(1+g) で積分
dg_, Lc, dpx_ = -1e-3, 1.0, 2e-3
from scipy.integrate import solve_ivp


def rhs(s, y):
    q, x = y
    gs = dg_ * s / Lc
    return [dpx_ / Lc / (1 + gs), q]          # 力 F = dpx/Lc は運動量に一定の衝撃、q = px/(1+g) の変化は F/(1+g)


sol = solve_ivp(rhs, [0, Lc], [0.0, 0.0], rtol=1e-12, atol=1e-16)
q1, x1 = sol.y[0, -1], sol.y[1, -1]
print('数値積分: q(L) - dpx = %.4e、公式 -(1/2) dpx dg = %.4e' % (q1 - dpx_, -0.5 * dpx_ * dg_))
print('          x(L) - dpx L/2 = %.4e、公式 (-(1/2) dpx dg) L/3 = %.4e (一様な力のときの x = dpx L/2 への加算分)' % (x1 - dpx_ * Lc / 2, -0.5 * dpx_ * dg_ * Lc / 3))

print('\n== (h) 初期の px ≠ 0 でも δ の分散は N cuu u_c² (tradke の dee = anp uc² cuu) ==')
px0 = 3e-3
a1 = np.zeros((20000, 6))
a1[:, 1] = px0
rn = np.asarray(track6.track_particles(els, a1.copy(), 1, ffs.momentum(), 0.0, False, True, False))
rf = np.asarray(track6.track_particles(els, a1.copy(), 1, ffs.momentum(), 0.0, False, True, True, seed=11))
dpx = rf[:, 1] - rf[:, 1].mean()
dd = rf[:, 5] - rf[:, 5].mean()
var_d = (dd * dd).mean()
cov = (dpx * dd).mean()
print('px0 = %.0e: Var(δ) = %.4e (N cuu u_c² = %.4e)' % (px0, var_d, N * 11 / 27 * uc**2))

print('\n== (g) リング全体: U0 = C_γ E⁴ I2/(2π) ==')
Cg = 4 * math.pi / 3 * re / (mc2 / 1e9) ** 3                    # [m/GeV³]
I2 = 16 * th**2 / 1.2
U0 = Cg * 3.0**4 * I2 / (2 * math.pi) * 1e9                  # eV
print('C_γ = 4π r_e/(3 (mc²)³) = %.5e m/GeV³、I2 = 16 θ²/L = %.5f 1/m、U0 = C_γ E⁴ I2/(2π) = %.4f MeV' % (Cg, I2, U0 / 1e6))
print('16 台の偏向磁石の追跡の損失の和 = 16 x %.6e x 3 GeV = %.4f MeV' % (-r1[:, 5].mean(), -16 * r1[:, 5].mean() * 3e9 / 1e6))
