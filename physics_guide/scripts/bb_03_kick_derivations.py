# -*- coding: utf-8 -*-
r"""BEAMBEAM の途中式の確認 (16_undulator_beambeam.md §6〜§9)。

(a) ポテンシャル U (16.6) が 2 次元の Poisson 方程式 ∇²U = 4π ρ を満たすこと (被積分関数の t 微分の恒等式)。
(b) 丸いビーム: U = Ein(r²/2σ²)、∂U_phys/∂σ² = e^h/σ² (sympy)。
(c) (16.8) の導出の各段: 変数 ρ = √(b/a) への変換、w 関数への書き換え。数値積分との比較。
(d) D = 時間 1 の流れ exp(F) の級数が 1 項と 2 項で止まること (sympy)。
(e) 条件付きガウス分布 (§7.2): 精度行列のブロックの逆行列 = Schur 補行列。傾き角の式。
(f) スライスの重心 (16.11) と境界。
(g) ビームストラールング: 係数の恒等式、光子スペクトルのモーメント 8/(15√3)、11/27 (mpmath)。
(h) 軌道のずれの符号 (x -= cod1、px += cod2): 強ビームの位置が弱粒子の座標系で X - P s* になっていること (pysad の BEAMBEAM で確認)。
"""
import math
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.integrate import quad
from scipy.special import wofz
import _common as C

C.start(__file__)

print('== (a) ∇²U = 4πρ ==')
xs, ys, t, sxx, syy = sp.symbols('x y t sigma_x sigma_y', positive=True)
a = 2 * sxx**2 + t
b = 2 * syy**2 + t
Eg = sp.exp(-xs**2 / a - ys**2 / b) / sp.sqrt(a * b)
g = -Eg                                          # U の被積分関数 = (1 - e^...) / √(ab) の x, y に依る部分
lap = sp.diff(g, xs, 2) + sp.diff(g, ys, 2)
print('∇²(被積分関数) + 4 d/dt[e^(..)/√(ab)] =', sp.simplify(lap + 4 * sp.diff(Eg, t)))
rho = sp.exp(-xs**2 / (2 * sxx**2) - ys**2 / (2 * syy**2)) / (2 * sp.pi * sxx * syy)
t0 = Eg.subs(t, 0)
print('4 e^(..)/√(ab)|_{t=0} - 4π ρ =', sp.simplify(4 * t0 - 4 * sp.pi * rho),
      ' (∫0^∞ d/dt[..] dt = 0 - [..]_{t=0} なので ∇²U = 4π ρ)')

print('\n== (b) 丸いビーム ==')
r2, sg2 = sp.symbols('r2 s2', positive=True)
w = r2 / (2 * sg2)
Ein = lambda w_: sp.integrate((1 - sp.exp(-sp.Symbol('q', positive=True))) / sp.Symbol('q', positive=True), (sp.Symbol('q', positive=True), 0, w_))
# U_round(r) = ∫0^∞ (1 - e^{-r²/(2σ²+t)}) /(2σ²+t) dt,  q = r²/(2σ²+t): dt = -r² dq/q²、1/(2σ²+t) = q/r²  ->  ∫0^{w} (1-e^{-q})/q dq
q = sp.Symbol('q', positive=True)
dUdw = (1 - sp.exp(-w)) / w
Uphys_d = dUdw * sp.diff(w, sg2) + sp.diff(sp.log(4 * sg2), sg2)
print('d/dσ² [Ein(r²/2σ²) + ln(4σ²)] - e^{-r²/2σ²}/σ² =', sp.simplify(Uphys_d - sp.exp(-w) / sg2))
mp.mp.dps = 20
for r_, s_ in ((0.3, 0.1), (2.0, 1.0)):
    f = lambda tt: (1 - mp.e**(-r_**2 / (2 * s_**2 + tt))) / (2 * s_**2 + tt)
    val = mp.quad(f, [0, 1, 100, mp.inf])
    ein = mp.quad(lambda qq: (1 - mp.e**(-qq)) / qq, [0, r_**2 / (2 * s_**2)])
    print('r = %.1f σ = %.1f: ∫dt = %s、Ein(r²/2σ²) = %s、遠方の形 ln(r²/2σ²) + γ_E = %s' % (
        r_, s_, mp.nstr(val, 12), mp.nstr(ein, 12), mp.nstr(mp.log(r_**2 / (2 * s_**2)) + mp.euler, 12)))

print('\n-- (b2) 楕円のビームの遠方: U + 2 ln(σx+σy) - ln r² -> ln 2 + γ_E (σ に依らない) --')
mp.mp.dps = 20
for sx, sy in ((1, 1), (2, 1), (3, 0.5), (1, 4)):
    row = []
    for (X_, Y_) in ((300.0, 0.0), (0.0, 300.0), (210.0, 210.0)):
        f = lambda tt: (1 - mp.e**(-X_**2 / (2 * sx**2 + tt) - Y_**2 / (2 * sy**2 + tt))) / mp.sqrt((2 * sx**2 + tt) * (2 * sy**2 + tt))
        v = mp.quad(f, [0, 10, 1e3, 1e5, 1e7, mp.inf])
        row.append(mp.nstr(v + 2 * mp.log(sx + sy) - mp.log(X_**2 + Y_**2), 8))
    print('σ = (%g, %g): (300,0)、(0,300)、(210,210) で' % (sx, sy), row, '(ln2 + γ_E = %s)' % mp.nstr(mp.log(2) + mp.euler, 8))

print('\n== (c) (16.8) の導出: t の積分 -> ρ = √(b/a) の積分 -> w 関数 ==')
sx_, sy_ = 0.2, 0.05
x0, y0 = 0.15, 0.08
S = math.sqrt(2 * (sx_**2 - sy_**2))


def dU_t(x_, y_):
    fx = lambda tt: 2 * x_ * math.exp(-x_**2 / (2 * sx_**2 + tt) - y_**2 / (2 * sy_**2 + tt)) / ((2 * sx_**2 + tt)**1.5 * (2 * sy_**2 + tt)**0.5)
    fy = lambda tt: 2 * y_ * math.exp(-x_**2 / (2 * sx_**2 + tt) - y_**2 / (2 * sy_**2 + tt)) / ((2 * sx_**2 + tt)**0.5 * (2 * sy_**2 + tt)**1.5)
    pieces = (0, 1e-4, 1e-3, 1e-2, 1e-1, 1, np.inf)
    vx = sum(quad(fx, pieces[i], pieces[i + 1], epsabs=0, epsrel=1e-12)[0] for i in range(len(pieces) - 1))
    vy = sum(quad(fy, pieces[i], pieces[i + 1], epsabs=0, epsrel=1e-12)[0] for i in range(len(pieces) - 1))
    return vx, vy


rho0 = sy_ / sx_
Eexp = lambda rr: math.exp(-(1 - rr**2) * (x0**2 + y0**2 / rr**2) / S**2)
vx_r = 4 * x0 / S**2 * quad(Eexp, rho0, 1, epsabs=0, epsrel=1e-13)[0]
vy_r = 4 * y0 / S**2 * quad(lambda rr: Eexp(rr) / rr**2, rho0, 1, epsabs=0, epsrel=1e-13)[0]
vx_t, vy_t = dU_t(x0, y0)
print('∂U/∂x: t の積分 = %.12e、ρ の積分 (4x/S²)∫e^E dρ = %.12e' % (vx_t, vx_r))
print('∂U/∂y: t の積分 = %.12e、ρ の積分 (4y/S²)∫e^E ρ^-2 dρ = %.12e' % (vy_t, vy_r))
# 指数の恒等式: -(1-ρ²)(x² + y²/ρ²)/S² = -[(x - i y)² + (y/ρ + i ρ x)²]/S²
rr = 0.7
lhs = -(1 - rr**2) * (x0**2 + y0**2 / rr**2) / S**2
rhs = -(((x0 - 1j * y0)**2) + (y0 / rr + 1j * rr * x0)**2) / S**2
print('指数の恒等式: 左辺 = %.12f、右辺 = %.12f %+.1ej' % (lhs, rhs.real, rhs.imag))
z1 = (x0 + 1j * y0) / S
z2 = (x0 * sy_ / sx_ + 1j * y0 * sx_ / sy_) / S
h = -x0**2 / (2 * sx_**2) - y0**2 / (2 * sy_**2)
cf = 2 * math.sqrt(math.pi) / S * (wofz(z1) - math.exp(h) * wofz(z2))
print('(16.8): (2√π/S)[w(ζ1) - e^h w(ζ2)] = %.12e %+.12ej  (∂U/∂y + i ∂U/∂x = %.12e %+.12ej)' % (cf.real, cf.imag, vy_t, vx_t))

print('\n== (d) D の流れ: F = s*(z)(px²+py²)/2、s* = (z - zs)/2、{g,F} = Σ(∂g/∂q ∂F/∂p - ∂g/∂p ∂F/∂q) ==')
X, Px, Y, Py, Z, Dl, zs = sp.symbols('x p_x y p_y z delta z_s')
Fd = (Z - zs) / 2 * (Px**2 + Py**2) / 2


def br(f, g_):
    return sum(sp.diff(f, q_) * sp.diff(g_, p_) - sp.diff(f, p_) * sp.diff(g_, q_) for q_, p_ in ((X, Px), (Y, Py), (Z, Dl)))


for name, f in (('x', X), ('y', Y), ('px', Px), ('py', Py), ('z', Z), ('δ', Dl)):
    L1 = br(f, Fd)
    L2 = br(L1, Fd)
    print('%-2s: {g,F} = %s、{{g,F},F} = %s' % (name, sp.simplify(L1), sp.simplify(L2)))

print('\n== (e) 条件付きガウス分布: 精度行列のブロックの逆 = Schur 補行列 ==')
rng = np.random.default_rng(2)
Ar = rng.normal(size=(6, 6))
Sig = Ar @ Ar.T + 0.5 * np.eye(6)
P = np.linalg.inv(Sig)
iz = 4                                    # z の添字 (x, px, y, py, z, δ)
idx = [0, 1, 2, 3, 5]
Saa = Sig[np.ix_(idx, idx)]
Saz = Sig[idx, iz]
cov_schur = Saa - np.outer(Saz, Saz) / Sig[iz, iz]
cov_prec = np.linalg.inv(P[np.ix_(idx, idx)])
mean_sig = Saz / Sig[iz, iz]
mean_prec = -np.linalg.solve(P[np.ix_(idx, idx)], P[idx, iz])
print('共分散: |(P_aa)^-1 - (Σ_aa - Σ_az Σ_zz^-1 Σ_za)| = %.2e' % np.abs(cov_prec - cov_schur).max())
print('平均 (z = 1 あたり): |-(P_aa)^-1 P_az - Σ_az/Σ_zz| = %.2e' % np.abs(mean_prec - mean_sig).max())
# 傾き: 回転で Σ13 を消す角
S2 = np.array([[3.0, 0.4], [0.4, 1.0]])
tt = 0.5 * math.atan2(2 * 0.4, 3.0 - 1.0)
Rm = np.array([[math.cos(tt), math.sin(tt)], [-math.sin(tt), math.cos(tt)]])
print('回転角 t = ½ atan2(2Σ13, Σ11-Σ33) = %.6f で R Σ R^T の非対角 = %.2e' % (tt, (Rm @ S2 @ Rm.T)[0, 1]))
print('Σ13 = 0、Σ11 < Σ33 のとき t = ½ atan2(0, 負) = %.6f (= π/2 = %.6f)' % (0.5 * math.atan2(0.0, -1.0), math.pi / 2))

print('\n== (f) スライスの重心 (16.11) ==')
from scipy.stats import norm
n = 5
edges = [-np.inf] + [norm.ppf(i / n) for i in range(1, n)] + [np.inf]
cen_formula = [n / math.sqrt(2 * math.pi) * (math.exp(-edges[i]**2 / 2) - math.exp(-edges[i + 1]**2 / 2)) for i in range(n)]
cen_int = [n * quad(lambda zz: zz * norm.pdf(zz), edges[i], edges[i + 1])[0] for i in range(n)]
print('n = 5: 公式 =', np.array2string(np.array(cen_formula), precision=6), '\n        数値積分 n∫z φ(z)dz =', np.array2string(np.array(cen_int), precision=6))

print('\n== (g) ビームストラールング ==')
alpha = 1 / 137.035999084
a_nrad_c = 5 * alpha / (2 * math.sqrt(3))
print('a_nrad c_ave c_uc / (α λ̄_C) = (5/(2√3)) (8/(15√3)) (3/2) = %.12f (= 2/3 = %.12f)' % (5 / (2 * math.sqrt(3)) * 8 / (15 * math.sqrt(3)) * 1.5, 2 / 3))
mp.mp.dps = 25
Mk = lambda nu: 2**(nu - 1) * mp.gamma((nu + mp.mpf(8) / 3) / 2) * mp.gamma((nu - mp.mpf(2) / 3) / 2)     # ∫0^∞ x^nu K_{5/3}(x) dx
for nu in (1, 2):
    print('∫x^%d K_{5/3} dx: 公式 %s、数値積分 %s' % (nu, mp.nstr(Mk(nu), 12), mp.nstr(mp.quad(lambda xx: xx**nu * mp.besselk(mp.mpf(5) / 3, xx), [0, 1, 10, 80]), 12)))
# S(ξ) = ∫_ξ^∞ K_{5/3} dx。部分積分: ∫0^∞ ξ^μ S(ξ) dξ = (1/(μ+1)) ∫0^∞ x^{μ+1} K_{5/3}(x) dx
mom = lambda m_: Mk(m_ + 1) / (m_ + 1)
print('<ξ> = ∫ξ S/∫S = %s (8/(15√3) = %s)、<ξ²> = %s (11/27 = %s)' % (
    mp.nstr(mom(1) / mom(0), 12), mp.nstr(mp.mpf(8) / (15 * mp.sqrt(3)), 12), mp.nstr(mom(2) / mom(0), 12), mp.nstr(mp.mpf(11) / 27, 12)))

print('\n== (h) 軌道のずれの符号 ==')
from pysad import beambeam as BB
mom = 1.7e9
base = dict(BX=1.0, BY=1.0, EMIX=1e-8, EMIY=1e-8, AX=0.0, AY=0.0, DP=1e-3, SIGZ=0.01, SLICE=1, NP=1e10)


def kick(prm, row):
    bb = BB.setup(prm, mom)
    S_ = np.zeros((1, 16))
    pn_ = 1 + row[5]
    S_[0, :6] = [row[0], row[1] / pn_, row[2], row[3] / pn_, row[4], row[5]]
    BB.apply(S_, 1, bb)
    return S_[0, 1] * (1 + S_[0, 5]), S_[0, 3] * (1 + S_[0, 5])


zp = 0.004                                # SLICE = 1 の重心は z_s = 0 なので s* = z/2 = 0.002
P = 1e-3
sz = zp / 2
x_in = 1.5e-4
kA = kick(dict(base, DPX=P), (x_in, 0, 0, 0, zp, 0))                 # 強ビームの軌道 (X, P) = (0, P)
kB = kick(dict(base), (x_in + P * sz, 0, 0, 0, zp, 0))              # 強ビームが軸上で、弱粒子の x が P s* だけずれている (X - P s* = -P s* を相対位置に直した)
kC = kick(dict(base), (x_in - P * sz, 0, 0, 0, zp, 0))              # 逆符号
print('A: DPX = %.0e (強ビームの傾き)、x = %.1e: Δpx = %.8e' % (P, x_in, kA[0]))
print('B: DPX = 0、x = x + P s* = %.4e:       Δpx = %.8e (A と一致するはず)' % (x_in + P * sz, kB[0]))
print('C: DPX = 0、x = x - P s* = %.4e:       Δpx = %.8e (符号を逆にすると一致しない)' % (x_in - P * sz, kC[0]))
