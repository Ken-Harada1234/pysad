# -*- coding: utf-8 -*-
r"""CSR (自由空間・定常) と Haissinski 方程式の導出の確認 (17_collective.md §11、§12)。

CSR:
 (a) Schwinger のスペクトル F(y) = y ∫_y^∞ K_{5/3}: 全パワーの規格化 ∫F = 8π/(9√3)、y << 1 の係数 F ≈ 4π/(√3 Γ(1/3)) (y/2)^{1/3}。
 (b) 単一電子の全パワー P = e² c γ⁴/(6π ε0 R²) と ∫ P_ω dω の一致。
 (c) インピーダンスの係数と、時間領域の係数 1/(2π ε0 3^{1/3} R^{2/3}) の一致 (Γ(1/3)Γ(2/3) = 2π/√3)。
 (d) ガウス分布のバンチの平均の損失を、周波数領域と時間領域の 2 通りで数値積分して比べる。
Haissinski:
 (e) Fokker-Planck の定常解 ψ = exp(-δ²/2σ² - U/(rσ²)) を sympy で。
 (f) pysad の csrhaissin (K64 の Newton 法) の解と、実空間の不動点反復の解 (方程式 ρ = K exp(-z²/2σ² - aw_H V)、V' = W*ρ - 平均) を比べる。
 (g) その解が、pysad の csrtrack (K64 の追跡) の平衡と合うか。aw_H = inten/(r65 σz²) (K64 の tfcsrhaissin の関係) を使う。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())
import mpmath as mp
import sympy as sp
from scipy import integrate

# ------------------------------------------------------------------ (a)
mp.mp.dps = 15
F = lambda y: y * mp.quad(lambda t: mp.besselk(mp.mpf(5) / 3, t), [y, mp.inf])
print('== (a) F(y) = y ∫_y^∞ K_{5/3}(t) dt ==')
tot = mp.quad(lambda t: t**2 * mp.besselk(mp.mpf(5) / 3, t) / 2, [0, 1, 10, mp.inf])   # ∫_0^∞ F dy = ∫_0^∞ dt K_{5/3}(t) ∫_0^t y dy = ∫ t² K_{5/3}/2
print('∫_0^∞ F dy = %s、8π/(9√3) = %s' % (mp.nstr(tot, 12), mp.nstr(8 * mp.pi / (9 * mp.sqrt(3)), 12)))
for y in (1e-2, 1e-3, 1e-4):
    asym = 4 * mp.pi / (mp.sqrt(3) * mp.gamma(mp.mpf(1) / 3)) * (mp.mpf(y) / 2)**(mp.mpf(1) / 3)
    print('y = %g: F = %s、4π/(√3Γ(1/3)) (y/2)^{1/3} = %s、比 %s' % (y, mp.nstr(F(y), 8), mp.nstr(asym, 8), mp.nstr(F(y) / asym, 8)))

# ------------------------------------------------------------------ (b)
print('小さい y の展開: K_{5/3}(t) ≈ (1/2)Γ(5/3) 2^{5/3} t^{-5/3} → F ≈ (3/4)Γ(5/3) 2^{5/3} y^{1/3}: 係数 %s、4π/(√3Γ(1/3)) 2^{-1/3} = %s' % (
    mp.nstr(mp.mpf(3) / 4 * mp.gamma(mp.mpf(5) / 3) * 2**(mp.mpf(5) / 3), 12), mp.nstr(4 * mp.pi / (mp.sqrt(3) * mp.gamma(mp.mpf(1) / 3)) * 2**(-mp.mpf(1) / 3), 12)))
print('\n== (b) 単一電子の全パワー: P_ω = (√3 e² γ/(8π² ε0 R)) F(ω/ω_c)、ω_c = 3γ³c/(2R) ==')
ratio = mp.sqrt(3) / (8 * mp.pi**2) * mp.mpf(3) / 2 * (8 * mp.pi / (9 * mp.sqrt(3))) * 6 * mp.pi
print('∫P_ω dω / (e²cγ⁴/(6π ε0 R²)) =', mp.nstr(ratio, 12), '(1 なら一致)')

# ------------------------------------------------------------------ (c)
print('\n== (c) 係数の恒等式 ==')
G13, G23 = mp.gamma(mp.mpf(1) / 3), mp.gamma(mp.mpf(2) / 3)
print('Γ(1/3)Γ(2/3) = %s、2π/√3 = %s' % (mp.nstr(G13 * G23, 12), mp.nstr(2 * mp.pi / mp.sqrt(3), 12)))
print('1/(√3 Γ(1/3)Γ(2/3)) = %s、1/(2π) = %s' % (mp.nstr(1 / (mp.sqrt(3) * G13 * G23), 12), mp.nstr(1 / (2 * mp.pi), 12)))
s_, k_, eps_ = sp.symbols('s k epsilon', positive=True)
mu = sp.Rational(2, 3)
val = sp.integrate(s_**(mu - 1) * sp.exp(-(eps_ - sp.I * k_) * s_), (s_, 0, sp.oo), conds='none')
print('∫_0^∞ s^{μ-1} e^{-(ε - ik)s} ds (μ = 2/3) =', sp.simplify(val), ' (= Γ(μ)(ε - ik)^{-μ})')
print('(-i)^{1/3} = e^{-iπ/6}: cos(π/6) = %.6f (= √3/2 = %.6f)、sin(π/6) = %.6f' % (float(sp.cos(sp.pi / 6)), math.sqrt(3) / 2, float(sp.sin(sp.pi / 6))))

# ------------------------------------------------------------------ (d)
print('\n== (d) ガウス分布のバンチの平均の損失: 周波数領域と時間領域 (単位 N e²/(4π ε0), R = 1, σ = 1) ==')
sig = 1.0
coef = 4 * math.pi / (2 * float(G13) * 3**(1 / 3))              # Re ŵ(k) = coef k^{1/3} ((e²/4πε0) の単位)
fd = -integrate.quad(lambda k: math.exp(-k**2 * sig**2) * k**(1 / 3), 0, np.inf)[0] * coef / math.pi
lam = lambda z: math.exp(-z * z / (2 * sig**2)) / (math.sqrt(2 * math.pi) * sig)
dlam = lambda z: -z / sig**2 * lam(z)


def dEds(z):
    # dE/ds(z) = +(2/3^{1/3}) ∫_0^∞ λ'(z+s) s^{-1/3} ds  ((e²N/4πε0) の単位)。s = t³: ∫ ... = 3 ∫_0^∞ t λ'(z + t³) dt
    return 2 / 3**(1 / 3) * 3 * integrate.quad(lambda t: t * dlam(z + t**3), 0, 6, epsabs=1e-13, epsrel=1e-12, limit=400)[0]


td = integrate.quad(lambda z: lam(z) * dEds(z), -8, 8, epsabs=1e-12, epsrel=1e-10, limit=400)[0]
closed = -float(G23) / (float(G13) * 3**(1 / 3))
print('周波数領域: %.8f、時間領域: %.8f、閉じた式 -Γ(2/3)/(3^{1/3}Γ(1/3)) = %.8f' % (fd, td, closed))
print('頭 (z = +2σ) と尾 (z = -2σ) の dE/ds (同じ単位、N e²/(4πε0) を除く): %.5f、%.5f' % (dEds(2.0), dEds(-2.0)))

# ------------------------------------------------------------------ (e)
print('\n== (e) Fokker-Planck の定常解 ==')
z, dl = sp.symbols('z delta', real=True)
r, d, sg = sp.symbols('r d sigma', positive=True)
U = sp.Function('U')
psi = sp.exp(-dl**2 / (2 * sg**2) - U(z) / (r * sg**2))
F_ = -sp.diff(U(z), z)
L = -r * dl * sp.diff(psi, z) - F_ * sp.diff(psi, dl) + 2 * d * sp.diff(dl * psi, dl) + 2 * d * sg**2 * sp.diff(psi, dl, 2)
print('L ψ / ψ =', sp.simplify(L / psi))

# ------------------------------------------------------------------ (f)(g)
from pysad import csr
nr = 512
dz = 0.15e-3
ell = 1.0e-3                                          # 指数関数ウェイクの減衰長
j = np.arange(nr)
wk = np.where(j >= 1, np.exp(-j * dz / ell), 0.0)     # 配列の向き: 源が z の小さい側 (csrwake の畳み込み v[n] = Σ w[n-m] ρ[m])
wk[0] = 0.5
wk[nr // 2:] = 0.0
zl = csr.trftr(wk.copy(), True)
zc = (nr // 4) * dz
zgrid = dz * np.arange(1, nr // 2 + 1) - zc


def solve_mine(awH, sigz0, nit=40000, relax=0.15):
    u0 = -(zgrid / sigz0)**2 / 2
    rho = np.exp(u0)
    rho /= rho.sum() * dz
    nrh = nr // 2
    for _ in range(nit):
        full = np.zeros(nr)
        full[:nrh] = rho
        v = np.real(np.fft.ifft(np.fft.fft(wk) * np.fft.fft(full)))
        dv = v - v.mean()                               # k = 0 を除く
        V = np.cumsum(dv) * dz
        un = u0 - awH * V[:nrh]
        new = np.exp(un - un.max())
        new /= new.sum() * dz
        if np.abs(new - rho).max() < 1e-13 * new.max():
            rho = new
            break
        rho = (1 - relax) * rho + relax * new
    return rho


def solve_root(awH, sigz0, guess):
    """同じ方程式を scipy.optimize.root (Newton 型) で解く。未知数 u = ln(ρ) (ρ は密度)。"""
    from scipy.optimize import root
    u0 = -(zgrid / sigz0)**2 / 2
    nrh = nr // 2

    def resid(u):
        lz = np.log(np.exp(u - u.max()).sum() * dz) + u.max()
        rho = np.exp(u - lz)
        full = np.zeros(nr)
        full[:nrh] = rho
        v = np.real(np.fft.ifft(np.fft.fft(wk) * np.fft.fft(full)))
        V = np.cumsum(v - v.mean()) * dz
        t = u0 - awH * V[:nrh]
        lzt = np.log(np.exp(t - t.max()).sum() * dz) + t.max()
        return (u - lz) - (t - lzt)
    sol = root(resid, np.log(np.maximum(guess, 1e-300)), method='hybr', tol=1e-11)
    rho = np.exp(sol.x)
    return rho / (rho.sum() * dz), sol.success


def moments(rho):
    m = np.sum(rho * zgrid) * dz
    return m, math.sqrt(np.sum(rho * zgrid**2) * dz - m * m)


mus = 0.04                                            # 1 周のシンクロトロン位相 [rad]
sde = 1.0e-3
sigz0 = 3.0e-3
r56 = 2 * math.sin(mus / 2) * sigz0 / sde             # r56 σδ/(2 sin(μ/2)) = σz0
r65 = -4 * math.sin(mus / 2)**2 / r56
damp = 4.0e-3                                         # 1 周の振幅減衰率
print('\n== (f)(g) Haissinski 方程式と K64 の追跡 ==')
print('μ_s = %.3f rad、r56 = %.4f、r65 = %.5f、σz0 = r56 σδ/(2 sin(μ/2)) = %.3e m、σδ = %.1e、damp = %.1e、格子 dz = %.2e m、ウェイク e^{-s/l} (l = %.1e m)' % (
    mus, r56, r65, sigz0, sde, damp, dz, ell))
npart = 20000
rng = np.random.default_rng(7)
print('\n  inten   aw_H = inten/(r65 σz0²)   csrhaissin: σz [m]   Newton 解: σz [m]   最大差/最大 ρ   追跡 σz [m] (最後の 1500 周の平均)   σz/σz0')
for inten in (0.0, -5.0e-8, -1.5e-7, -3.0e-7):
    awH = inten / (r65 * sigz0**2)
    rh = None if inten == 0 else csr.csrhaissin(zl, awH, dz, sigz0)
    rm = solve_mine(awH, sigz0)
    rm, okroot = solve_root(awH, sigz0, rh[:nr // 2] / (rh[:nr // 2].sum() * dz) if rh is not None else rm)
    if rh is None:
        s_h, dd = float('nan'), float('nan')
        if inten == 0:
            s_h, dd = sigz0, 0.0
    else:
        s_h = moments(rh[:nr // 2])[1]
        dd = np.abs(rh[:nr // 2] - rm).max() / rm.max()
    sm = moments(rm)[1]
    cm = moments(rm)[0]
    ch = float('nan') if rh is None else moments(rh[:nr // 2])[0]
    # 追跡 (K64 の csrtrack の直訳)。初期分布はウェイクなしのガウス分布
    zp = rng.normal(0, sigz0, npart)
    dpp = rng.normal(0, sde, npart)
    nt = int(12 / damp)
    gauss = lambda n: rng.normal(size=n)
    dat = csr.csrtrack(zp, dpp, zl, dz, 1, inten / npart, r56, r65, 0.0, 0.0, damp, sde, nt, gauss)
    st = dat[0, -1500:].mean()
    print('  %+.1e  %+.4e        %.4e          %.4e         %.2e       %.4e                 %.4f   <z>: csrhaissin %.3e、反復 %.3e (root %s)' % (inten, awH, s_h, sm, dd, st, st / sigz0, ch, cm, okroot))

# ------------------------------------------------------------------ (h)
print()
print("== (h) 追跡の揺らぎ: 定常解が安定なとき追跡の σz は揺らがず、不安定になると揺らぐ。σz と σδ の最後の区間の平均と、σz の時間変動の相対値 ==")
for awH_t, mus2, damp2 in ((0.4167, 0.04, 4.0e-3), (1.2501, 0.04, 4.0e-3), (2.5002, 0.04, 4.0e-3), (2.5002, 0.02, 2.0e-3), (3.75, 0.04, 4.0e-3)):
    r56b = 2 * math.sin(mus2 / 2) * sigz0 / sde
    r65b = -4 * math.sin(mus2 / 2)**2 / r56b
    inten2 = awH_t * r65b * sigz0**2
    rng2 = np.random.default_rng(21)
    zp = rng2.normal(0, sigz0, npart)
    dpp = rng2.normal(0, sde, npart)
    nt = int(12 / damp2)
    dat = csr.csrtrack(zp, dpp, zl, dz, 1, inten2 / npart, r56b, r65b, 0.0, 0.0, damp2, sde, nt, lambda n: rng2.normal(size=n))
    tail = int(1500 * 4e-3 / damp2)
    zt, dt = dat[0, -tail:], dat[1, -tail:]
    rh_t = csr.csrhaissin(zl, awH_t, dz, sigz0)
    sh = moments(rh_t[:nr // 2])[1] if rh_t is not None else float('nan')
    print('  aw_H = %.4f、μ_s = %.2f、d = %.1e: 追跡 σz = %.4e m (時間変動 %.2f%%)、σδ = %.4e (ウェイクなし %.1e)、csrhaissin の σz = %.4e m (比 %.4f)' % (
        awH_t, mus2, damp2, zt.mean(), 100 * zt.std() / zt.mean(), dt.mean(), sde, sh, zt.mean() / sh))
