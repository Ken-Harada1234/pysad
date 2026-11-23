# -*- coding: utf-8 -*-
r"""空間電荷 SPAC・PSPAC の導出の確認 (17_collective.md §9、§10)。

(a) 一様な円柱: 実験室系の電場・磁場から x'' = K x/a²、K = 2 r_e λ/(β²γ³) (sympy)。ガウス分布の小振幅の極限との整合。
(b) 丸いガウス分布の電位: 壁 (半径 b) を基準にした電位が (λ'/4π ε0)[Ein(b²/2σ²) - Ein(r²/2σ²)] で、K64 の U (U(0) = 0) との差が縦の力の基準 (sympy と数値)。
(c) 3 次元ガウス分布の電位 φ = (Q/4πε0) (2/√π) ∫_0^∞ dτ Π (1 + 2σ_i²τ²)^{-1/2} exp(-x_i²τ²/(1 + 2σ_i²τ²)) が Poisson 方程式を満たすこと (数値)、
    σ_ζ → ∞ で 2 次元ガウス分布 × 線密度に近づくこと (WSPAC・PSPAC の近似の条件)。
(d) 接地した矩形の箱の中の 2 次元 Poisson 方程式の正弦級数解 (PSPAC の基礎) を、自由空間の解と比べる (像の電荷の効果)。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
import sympy as sp
from scipy import integrate, special

# ------------------------------------------------------------------ (a)
print('== (a) 一様な円柱 ==')
re_, lam, a, r, beta, gam, eps0, e, m, c = sp.symbols('r_e lambda a r beta gamma epsilon_0 e m c', positive=True)
lamp = lam / gam                                   # 静止系の線密度
Er_rest = e * lamp * r / (2 * sp.pi * eps0 * a**2)   # 一様な円柱の内部の電場 (Gauss の法則。λ' は粒子の線密度 [1/m]、電荷 e λ')
F = e * (gam * Er_rest - beta * c * (beta * gam * Er_rest / c))   # 実験室系: E = γE'、B = βγE'/c
dpp = F * (1 / (beta * c)) / (gam * m * beta * c)  # Δp/p0 per length (Δt = L/(βc)、p0 = γ m β c)
r_e_expr = e**2 / (4 * sp.pi * eps0 * m * c**2)
K = sp.simplify(dpp.subs(eps0, e**2 / (4 * sp.pi * re_ * m * c**2)) * a**2 / r)
print('x\'\' = K x/a² (Δp/p0 の長さあたり): K =', K)
gsub = 1 / sp.sqrt(1 - beta**2)
print('  K - 2 r_e λ/(β²γ³) (γ → 1/√(1-β²) を入れて):', sp.simplify((K - 2 * re_ * lam / (beta**2 * gam**3)).subs(gam, gsub)))

# ------------------------------------------------------------------ (b)
print('\n== (b) 丸いガウス分布の電位と壁の基準 ==')
rr, sg, b = sp.symbols('rr sigma b', positive=True)
x_ = rr**2 / (2 * sg**2)
Er = (1 - sp.exp(-x_)) / rr                         # (λ'/2πε0) = 1 の単位
Ein = lambda xx: sp.Integral((1 - sp.exp(-sp.Symbol('tt'))) / sp.Symbol('tt'), (sp.Symbol('tt'), 0, xx))
# φ(r) = ∫_r^b Er dr' = (1/2)[Ein(b²/2σ²) - Ein(r²/2σ²)]   (λ'/2πε0 の単位。Ein(x) = ∫_0^x (1 - e^{-t})/t dt、r' → t = r'²/2σ² で ∫ (1 - e^{-t})/(2t) dt)
for rv, bv in ((0.5, 20.0), (1.5, 20.0), (3.0, 20.0)):
    num = integrate.quad(lambda t: (1 - math.exp(-t**2 / 2)) / t, rv, bv, epsabs=0, epsrel=1e-12)[0]
    ein = lambda xx: integrate.quad(lambda t: (1 - math.exp(-t)) / t, 0, xx, epsabs=0, epsrel=1e-13)[0]
    print('  r/σ = %.1f、b/σ = %.0f: ∫_r^b E_r dr = %.10f、(1/2)[Ein(b²/2σ²) - Ein(r²/2σ²)] = %.10f' % (rv, bv, num, 0.5 * (ein(bv**2 / 2) - ein(rv**2 / 2))))
bv = 20.0
print('  Ein(b²/2σ²) と ln(b²/2σ²) + γ_E の差 (b/σ = 20):', '%.3e' % (integrate.quad(lambda t: (1 - math.exp(-t)) / t, 0, bv**2 / 2, epsabs=0, epsrel=1e-13)[0] - (math.log(bv**2 / 2) + float(sp.EulerGamma))))
print('  → φ(0) - φ(壁) = (λ\'/4πε0)[ln(b²/2σ²) + γ_E]: 縦の電場の主要項 ln(b/σ)。K64 の U (U(0) = 0) はこの定数を落とす。')

# ------------------------------------------------------------------ (c)
print('\n== (c) 3 次元ガウス分布の電位 ==')
# 単位: Q/(4π ε0) = 1。φ = (2/√π) ∫_0^∞ dτ Π_i (1 + 2σ_i² τ²)^{-1/2} exp(-x_i² τ²/(1 + 2σ_i² τ²))
def phi3(x, y, z, sx, sy, sz):
    f = lambda t: (1 + 2 * sx**2 * t**2)**-0.5 * (1 + 2 * sy**2 * t**2)**-0.5 * (1 + 2 * sz**2 * t**2)**-0.5 * math.exp(
        -x**2 * t**2 / (1 + 2 * sx**2 * t**2) - y**2 * t**2 / (1 + 2 * sy**2 * t**2) - z**2 * t**2 / (1 + 2 * sz**2 * t**2))
    v = 0.0
    lo = 0.0
    for hi in (1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4, 1e5, np.inf):
        v += integrate.quad(f, lo, hi, epsabs=0, epsrel=1e-12, limit=200)[0]
        lo = hi
    return 2 / math.sqrt(math.pi) * v


# 1/|r| のガウス表示の確認: 点電荷 (σ → 0) で 1/r
print('  σ → 0 の極限: φ(1, 2, 2) × 3 =', '%.8f' % (phi3(1.0, 2.0, 2.0, 1e-6, 1e-6, 1e-6) * 3), '(1/r = 1/3 から 1 になる)')
# Poisson 方程式: ∇²φ = -4π ρ (Q/(4π ε0) = 1、Q = 1)
sx, sy, sz = 1.0, 0.6, 1.7
x0, y0, z0 = 0.5, 0.3, 0.2
h = 0.05
lap = (phi3(x0 + h, y0, z0, sx, sy, sz) + phi3(x0 - h, y0, z0, sx, sy, sz) + phi3(x0, y0 + h, z0, sx, sy, sz) + phi3(x0, y0 - h, z0, sx, sy, sz)
       + phi3(x0, y0, z0 + h, sx, sy, sz) + phi3(x0, y0, z0 - h, sx, sy, sz) - 6 * phi3(x0, y0, z0, sx, sy, sz)) / h**2
rho = math.exp(-x0**2 / (2 * sx**2) - y0**2 / (2 * sy**2) - z0**2 / (2 * sz**2)) / ((2 * math.pi)**1.5 * sx * sy * sz)
print('  ∇²φ = %.6f、-4π ρ = %.6f (差分の刻み h = %.2f σ)' % (lap, -4 * math.pi * rho, h))
# 長いバンチの極限: x 方向の電場 E_x(x, y, 0) を 3 次元の式の数値微分と、(線密度 × 2 次元の場) で比べる
print('  長いバンチの極限 (σx = 1、σy = 0.6、E_x(0.8, 0.4, 0)): 3 次元 / [λ(0) × 2 次元ガウス] の比')
for szz in (3.0, 10.0, 30.0, 100.0):
    hh = 1e-4
    Ex3 = -(phi3(0.8 + hh, 0.4, 0.0, 1.0, 0.6, szz) - phi3(0.8 - hh, 0.4, 0.0, 1.0, 0.6, szz)) / (2 * hh)
    # 2 次元の場 (単位 Q/(4πε0) = 1、線密度 λ(0) = 1/(√(2π) σζ)): E_x = ∂U/∂x (U は (16.6) の U。単位線電荷の場 2/r (1 - e^-x))
    def Ex2(x, y):
        f = lambda t: (2 * x / (2 * 1.0**2 + t)) * math.exp(-x**2 / (2 + t) - y**2 / (2 * 0.36 + t)) / math.sqrt((2 + t) * (2 * 0.36 + t))
        v = 0.0
        lo = 0.0
        for hi in (1e-3, 1e-2, 1e-1, 1, 10, 100, np.inf):
            v += integrate.quad(f, lo, hi, epsabs=0, epsrel=1e-12, limit=200)[0]
            lo = hi
        return v
    lam0 = 1 / (math.sqrt(2 * math.pi) * szz)
    print('    σζ = %5.0f σx: E_x(3D) = %.6e、λ(0) E_x(2D) = %.6e、比 %.5f' % (szz, Ex3, lam0 * Ex2(0.8, 0.4), Ex3 / (lam0 * Ex2(0.8, 0.4))))

# ------------------------------------------------------------------ (d)
print('\n== (d) 接地した正方形の箱 (一辺 a、中心に丸いガウス分布 σ) の正弦級数と、自由空間の解 ==')
# ε0 = 1、線電荷 1。∇²Φ = -ρ、Φ = 0 (壁)。Φ = Σ Φ_mn sin(mπx/a) sin(nπy/a)、ρ_mn = g_m g_n、g_m = (2/a) sin(mπ/2) exp(-m²π²σ²/(2a²))
def box_field(x, y, aa, sig, M=401):
    m = np.arange(1, M + 1, 2)                      # 奇数だけ (中心の対称性)
    gm = (2 / aa) * np.sin(m * np.pi / 2) * np.exp(-(m * np.pi * sig / aa)**2 / 2)
    kx = m * np.pi / aa
    Mm, Nn = np.meshgrid(m, m, indexing='ij')
    Phimn = (gm[:, None] * gm[None, :]) / ((np.pi / aa)**2 * (Mm**2 + Nn**2))
    sx_, cx_ = np.sin(kx * x), np.cos(kx * x)
    sy_, cy_ = np.sin(kx * y), np.cos(kx * y)
    Ex = np.einsum('m,n,mn->', kx * cx_, sy_, Phimn)
    Ey = np.einsum('m,n,mn->', sx_, kx * cy_, Phimn)
    return Ex, Ey                                    # E = -∇Φ なので E_x = +Σ k cos sin Φ... (符号: -∂/∂x [sin] = -k cos、Φ の符号を含めて下で整える)


sig = 1.0
for aa in (20.0, 40.0, 80.0):
    out = []
    for rv in (1.0, 2.0, 4.0):
        xx, yy = aa / 2 + rv, aa / 2
        Ex, Ey = box_field(xx, yy, aa, sig)
        Efree = (1 - math.exp(-rv**2 / (2 * sig**2))) / (2 * math.pi * rv)
        out.append((-Ex, Efree))                     # 中心からの距離 rv の x 方向。-∂Φ/∂x = -Σ k cos(kx) sin(ky) Φ_mn (cos は中心で符号が交代する) で符号は級数から
    print('  a/σ = %3.0f: ' % (aa / sig) + '、'.join('r = %.0fσ: 箱 %.6e / 自由空間 %.6e (比 %.5f)' % (rv, o[0], o[1], o[0] / o[1]) for rv, o in zip((1, 2, 4), out)))
print('  (符号は級数の実計算の結果をそのまま出した。比が 1 に近づけば規格化と係数 1/(π²(m²+n²)/a²) が正しい。差は壁の像の電荷で、4 回対称の正方形の箱では (r/a)^4 に比例して効く: 上の差は a を 2 倍にすると約 16 分の 1)')

# ------------------------------------------------------------------ (e)
print('\n== (e) 楕円のガウス分布の小振幅の線形レンズ (17.8) の係数 ==')
t_ = sp.symbols('t', positive=True)
A_, B_ = sp.symbols('A B', positive=True)
Iab = sp.integrate((A_ + t_)**sp.Rational(-3, 2) * (B_ + t_)**sp.Rational(-1, 2), (t_, 0, sp.oo))
print('∫_0^∞ dt (A+t)^{-3/2} (B+t)^{-1/2} =', sp.simplify(Iab), ' (2/(√A(√A+√B)) との差:', sp.simplify(Iab - 2 / (sp.sqrt(A_) * (sp.sqrt(A_) + sp.sqrt(B_)))), ')')
sxs, sys_ = sp.symbols('sigma_x sigma_y', positive=True)
print('A = 2σx²、B = 2σy² を入れて、∂U/∂x = 2x × 積分 = 2x/(σx(σx+σy)):', sp.simplify(2 * Iab.subs({A_: 2 * sxs**2, B_: 2 * sys_**2}) - 2 / (sxs * (sxs + sys_))))
# Laslett: 薄いレンズ (x'' = +k x) のあとのチューンシフト
mu0, bt, al_, kk = sp.symbols('mu0 beta alpha k', real=True)
M0 = sp.Matrix([[sp.cos(mu0) + al_ * sp.sin(mu0), bt * sp.sin(mu0)], [-(1 + al_**2) / bt * sp.sin(mu0), sp.cos(mu0) - al_ * sp.sin(mu0)]])
L = sp.Matrix([[1, 0], [kk, 1]])
tr = sp.simplify((M0 * L).trace() / 2)
print('cos μ\' = tr(M0 L)/2 =', tr, ' → cos μ0 に足されるのは', sp.simplify(tr - sp.cos(mu0)), '(= β k sin μ0 / 2)')

# ------------------------------------------------------------------ (f)
print()
print("== (f) 円筒の壁のグリーン関数: 加法定理 K_0(k|r - r'|) = Σ ε_m K_m(k r_>) I_m(k r_<) cos m(θ - θ') の数値確認 ==")
kq, rr1, rr2, th_ = 3.0, 0.9, 0.4, 0.7
d_ = math.sqrt(rr1**2 + rr2**2 - 2 * rr1 * rr2 * math.cos(th_))
ssum = special.k0(kq * rr1) * special.iv(0, kq * rr2) + sum(2 * special.kv(m, kq * rr1) * special.iv(m, kq * rr2) * math.cos(m * th_) for m in range(1, 40))
print('  K_0(k d) = %.12f、級数 = %.12f' % (special.k0(kq * d_), ssum))
print("  壁 (r = R) で消える像の項 -Σ ε_m K_m(kR)/I_m(kR) I_m(kr) I_m(kr') cos m(θ-θ') を足した結果 (r = R): 自由空間の K_0 - 像 = 0 になるか")
RR = 1.0
rr3 = 0.4
for thv in (0.0, 1.0, 2.5):
    dd_ = math.sqrt(RR**2 + rr3**2 - 2 * RR * rr3 * math.cos(thv))
    free = special.k0(kq * dd_)
    img = special.k0(kq * RR) * special.iv(0, kq * rr3) + sum(2 * special.kv(m, kq * RR) / special.iv(m, kq * RR) * special.iv(m, kq * RR) * special.iv(m, kq * rr3) * math.cos(m * thv) for m in range(1, 40))
    print("  θ - θ' = %.1f: 自由空間 K_0 - 像 = %.3e" % (thv, free - img))
