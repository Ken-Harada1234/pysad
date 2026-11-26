# -*- coding: utf-8 -*-
"""放射の章 (13_radiation.md) の途中式の確認 (本文の各式を 1 行ずつ確かめる)。

(a) dv の式の差 (tradk1 の h2 + p2 と (1.9) の h2 + pr h0) を、桁落ちしない mpmath (40 桁) で。見積もり 1/(4γ²) との比較。
(b) Liénard のパワー (相対論的 Larmor から) と、スペクトルの規格化 (パワーの総和 = P)、光子の発生率 P/<u> との一致。
(c) 複合ポアソン過程の平均と分散 (母関数・全分散の法則・モンテカルロ)。
(d) drndsr の確率密度の規格化定数 3/(5π)。
(e) O(1/γ²) の補正: dp/ds ∝ -p⁴ (1 + 1/p²) の常微分方程式の数値解と、u⁴ のみの解の差。
(f) tradke の微分の途中式: dpz、xpa の微分 (有限差分との比較)、dθ = dxpa/√(1 - xpa²)。
(g) (13.7) の途中式: 軌道の s 平均、ソレノイドの項、dbeam(8) の Fortran の式と正しい式との差 (sympy)、
    Δδ と p が独立なときの 2 次モーメントのモンテカルロ。
(h) 薄い偏向磁石 (N_γ ≈ 0.06) の FLUC で ⟨Δp_x Δδ⟩ / Var(δ) など (13.7) の係数の確認 (pysad の追跡)。
"""
import io
import math
import numpy as np
import sympy as sp
import mpmath as mp
from scipy import integrate, special
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6

mc2 = 0.51099895e6
re = 2.8179403262e-15
alpha = 1 / 137.035999084

# ---------------------------------------------------------------------------------------------------- (a)
print('== (a) dv の式の差 (桁落ちしない計算、mpmath 40 桁) ==')
mp.mp.dps = 40


def dv_set(p0, pr):
    p0, pr = mp.mpf(p0), mp.mpf(pr)
    h0 = mp.sqrt(1 + p0**2)
    h = mp.sqrt(1 + (p0 * pr)**2)
    ex = 1 - (p0 * pr / h) / (p0 / h0)                                  # 厳密 1 - β/β0
    d19 = -(pr - 1) * (1 + pr) / h / (h + pr * h0)                      # (1.9)
    dtr = -(pr - 1) * (1 + pr) / h / (h + p0 * pr)                      # tradk1 の式
    return ex, d19, dtr, h0


for p0, dl, lab in ((3e9 / mc2, -1e-4, '3 GeV、δ = -1e-4'), (50.0, -0.01, 'p0 = 50 mc、δ = -0.01'), (1.0, 0.1, 'p0 = 1 mc、δ = +0.1')):
    ex, d19, dtr, h0 = dv_set(p0, 1 + dl)
    est = 1 / (4 * mp.mpf(p0)**2)
    print('  %-22s: (1.9) - 厳密 = %s (相対 %s)、tradk1 - (1.9) の相対 = %s、見積もり 1/(4 p0²) = %s' % (
        lab, mp.nstr(d19 - ex, 3), mp.nstr(abs((d19 - ex) / ex), 3), mp.nstr((dtr - d19) / d19, 5), mp.nstr(est, 5)))
print('  (解析: (tradk1 - (1.9))/(1.9) = (h + pr h0)/(h + p0 pr) - 1 = pr (h0 - p0)/(h + p0 pr)。h0 - p0 = 1/(h0 + p0) ≈ 1/(2 p0)、h + p0 pr ≈ 2 p0 で 1/(4 p0²)。)')
ex, d19, dtr, h0 = dv_set(3e9 / mc2, 1 - 1e-4)
print('  通常の倍精度での 1 - β/β0 (3 GeV): %.6e (厳密 %s)。倍精度の桁落ちで相対 1e-5 程度ずれる。' % (
    1 - (3e9 / mc2 * (1 - 1e-4) / math.sqrt(1 + (3e9 / mc2 * (1 - 1e-4))**2)) / (3e9 / mc2 / math.sqrt(1 + (3e9 / mc2)**2)), ex))

# ---------------------------------------------------------------------------------------------------- (b)
print('\n== (b) Liénard のパワーとスペクトルの規格化 ==')
g_, b_, a_, rho_, c_, k_ = sp.symbols('gamma beta a rho c k', positive=True)
v2a2 = sp.symbols('v2a2')
# 相対論的 Larmor: P = (2/3) k γ⁶ [a² - |v×a|²/c²]  (k = e²/(4π ε0 c³))。a ⊥ v のとき |v×a| = v a
P_larmor = sp.Rational(2, 3) * k_ * g_**6 * (a_**2 - (b_ * c_ * a_)**2 / c_**2)
P_simpl = sp.simplify(P_larmor.subs(b_, sp.sqrt(1 - 1 / g_**2)))
print('  P = (2/3) k γ⁶ [a² - (v a)²/c²] (a ⊥ v、β² = 1 - 1/γ²) = %s' % P_simpl)
P_circ = sp.simplify(P_simpl.subs(a_, b_ * b_ * c_**2 / rho_).subs(b_, sp.sqrt(1 - 1 / g_**2)))
print('  a = β² c²/ρ を入れる: P = %s  (k = e²/(4π ε0 c³) = r_e mc²/c³ を使うと (2/3) r_e mc² c γ⁴β⁴/ρ²)' % P_circ)
mp.mp.dps = 20


def Fsc(x):
    return x * integrate.quad(lambda t: special.kv(5.0 / 3.0, t), x, np.inf, epsabs=1e-14, epsrel=1e-12)[0]


pref = 9 * math.sqrt(3) / (8 * math.pi)
m1 = integrate.quad(Fsc, 0, 60, epsabs=1e-13, epsrel=1e-11, limit=200)[0]
print('  ∫ F dξ (scipy の直接積分) = %.10f、8π/(9√3) = %.10f、(9√3/8π) ∫F dξ = %.10f' % (m1, 8 * math.pi / (9 * math.sqrt(3)), pref * m1))
m0 = integrate.quad(lambda x: Fsc(x) / x, 0, 60, epsabs=1e-13, epsrel=1e-11, limit=200)[0]
print('  ∫ F/ξ dξ = ∫ G dξ = %.10f (5π/3 = %.10f)、(9√3/8π) ∫F/ξ dξ = %.8f (15√3/8 = %.8f)。' % (m0, 5 * math.pi / 3, pref * m0, 15 * math.sqrt(3) / 8))
print('  光子の発生率 Ṅ = (P/u_c)(15√3/8) = P/<u> (<u> = 8 u_c/(15√3)): 15√3/8 × 8/(15√3) = %s' % mp.nstr(mp.mpf(1), 5))
rho0, gam = 3.0558, 3e9 / mc2
Pw = 2 / 3 * re * mc2 * gam**4 / rho0**2     # eV·(c)/m、c をかけて eV/s にする
print('  3 GeV、ρ = 3.0558 m: u_c = (3/2) λ_C mc² γ³/ρ = %.4f keV' % (1.5 * re / alpha * mc2 * gam**3 / rho0 / 1e3))
print('  漸近形: ξ ≪ 1 で F ≈ 2.1495 ξ^(1/3)、ξ ≫ 1 で F ≈ √(πξ/2) e^(-ξ) (数値: ξ = 1e-3 で F/ξ^(1/3) = %.4f、ξ = 20 で F/[√(πξ/2) e^(-ξ)] = %.4f)' % (Fsc(1e-3) / 1e-3**(1 / 3), Fsc(20) / (math.sqrt(math.pi * 20 / 2) * math.exp(-20))))

# ---------------------------------------------------------------------------------------------------- (c)
print('\n== (c) 複合ポアソン過程の分散 ==')
t, Nn, m1s, m2s = sp.symbols('t N m1 m2')
M = 1 + m1s * t + m2s * t**2 / 2                          # 1 個の光子のエネルギー u の積率母関数 E[e^{tu}] の 2 次まで
K = Nn * (M - 1)                                           # 複合ポアソンの累積母関数 log E[e^{t Σu}] = N (E[e^{tu}] - 1)
print('  累積母関数 K(t) = N (M(t) - 1): K\'(0) = %s、K\'\'(0) = %s' % (sp.diff(K, t).subs(t, 0), sp.diff(K, t, 2).subs(t, 0)))
Eu, Vu, EN, VN = sp.symbols('Eu Vu EN VN')
print('  全分散の法則: Var = E[n] Var(u) + Var(n) E[u]² で E[n] = Var(n) = N のとき = %s' % sp.simplify((Nn * (m2s - m1s**2) + Nn * m1s**2)))
# モンテカルロ: u/u_c の分布 (K5/3 の確率密度) を数表の逆関数で
# G(y) = ∫_y^∞ K_{5/3} dt を数表にする (y = w³ の等間隔の w で。pp dy = pp 3w² dw は w → 0 で有限)
w = np.linspace(1e-6, 40**(1 / 3), 4001)
yw = w**3
Gw = np.array([integrate.quad(lambda t: special.kv(5 / 3, t), y, np.inf, epsabs=1e-13, epsrel=1e-10)[0] for y in yw])
pp = 3 / (5 * math.pi) * Gw                                          # pp(y) = (3/(5π)) G(y)
ppw = pp * 3 * w**2                                                  # pp dy/dw
cdf = np.concatenate(([0], np.cumsum((ppw[1:] + ppw[:-1]) / 2 * np.diff(w))))
print('  数表: ∫ pp dy (y ≤ 40) = %.6f (規格化 1)。数表の平均 = %.5f (8/(15√3) = %.5f)、2 次モーメント = %.5f (11/27 = %.5f)' % (
    cdf[-1], np.trapezoid(ppw * yw, w), 8 / (15 * math.sqrt(3)), np.trapezoid(ppw * yw**2, w), 11 / 27))
rng = np.random.default_rng(3)
Ng = 24.283
ntr = 200000
nph = rng.poisson(Ng, ntr)
tot = nph.sum()
us = np.interp(rng.random(tot), cdf / cdf[-1], yw)
idx = np.repeat(np.arange(ntr), nph)
sums = np.bincount(idx, weights=us, minlength=ntr)
print('  N = %.3f、試行 %d: 平均 %.4f (N <u> = %.4f)、分散 %.4f (N <u²> = %.4f、N (<u²> - <u>²) = %.4f)' % (
    Ng, ntr, sums.mean(), Ng * 8 / (15 * math.sqrt(3)), sums.var(), Ng * 11 / 27, Ng * (11 / 27 - (8 / (15 * math.sqrt(3)))**2)))
print('  → 分散は N <u²> (n を固定した N <u²> - N <u>² ではない)。標本誤差 ~ %.2f%%' % (100 * math.sqrt(2 / ntr) * 1.5))

# ---------------------------------------------------------------------------------------------------- (d)
print('\n== (d) drndsr の確率密度 pp(y) = c ∫_y^∞ K_{5/3} の規格化定数 ==')
mp.mp.dps = 20
tK = mp.quad(lambda x: x * mp.besselk(mp.mpf(5) / 3, x), [0, 1, 10, 60])
print('  ∫_0^∞ t K_{5/3}(t) dt = %s (Γ(11/6)Γ(1/6) = %s、5π/3 = %s) → c = 3/(5π) = %s' % (
    mp.nstr(tK, 12), mp.nstr(mp.gamma(mp.mpf(11) / 6) * mp.gamma(mp.mpf(1) / 6), 12), mp.nstr(5 * mp.pi / 3, 12), mp.nstr(3 / (5 * mp.pi), 12)))
print('  ∫_0^∞ dy ∫_y^∞ K dt = ∫_0^∞ dt K(t) ∫_0^t dy = ∫ t K dt (積分の順序を入れ替え)。')
print('  c ∫_0^∞ y pp 型の平均: ∫ y G dy / ∫ G dy = (1/2)∫ t² K / ∫ t K = %s (8/(15√3) = %s)' % (
    mp.nstr(mp.quad(lambda x: x**2 * mp.besselk(mp.mpf(5) / 3, x), [0, 1, 10, 60]) / 2 / tK, 10), mp.nstr(8 / (15 * mp.sqrt(3)), 10)))

# ---------------------------------------------------------------------------------------------------- (e)
print('\n== (e) γ ∝ p からの補正 (O(1/γ²)) ==')
# mc 単位の運動量 p、h = √(1 + p²)、β = p/h。曲率半径 ρ 一定で P = (2/3) r_e mc² c γ⁴β⁴/ρ²、dp/ds = -P/(v² (mc)) より
# dp/ds = -(2/3)(r_e/ρ²) γ⁴β² = -(2/3)(r_e/ρ²) h² p² = -(2/3)(r_e/ρ²) p⁴ (1 + 1/p²)
ps = sp.symbols('p', positive=True)
hh = sp.sqrt(1 + ps**2)
bb = ps / hh
print('  γ⁴β²/p⁴ - (1 + 1/p²) = %s (h² = 1 + p²、β = p/h)' % sp.simplify(hh**4 * bb**2 / ps**4 - (1 + 1 / ps**2)))
print('  展開: (1 + 1/p²) = 1 + 1/p² (厳密。補正は 1/p² の 1 次で終わる)')
from scipy.integrate import solve_ivp


def dg_ode(de_, p1, exact):
    def f(tau, u):
        uu = u[0]
        fac = ((1 + 1 / (uu * p1)**2) / (1 + 1 / p1**2)) if exact else 1.0
        return [-de_ * uu**4 * fac]
    return solve_ivp(f, [0, 1], [1.0], rtol=1e-13, atol=1e-17).y[0, -1] - 1


for p1, de_ in ((5870.85, 5e-5), (50.0, 5e-3), (10.0, 5e-3), (3.0, 5e-3)):
    a1 = dg_ode(de_, p1, False)
    a2 = dg_ode(de_, p1, True)
    print('  p1 = %8.2f、de = %.0e: dg (u⁴ のみ) = %.10e、dg (補正 (1 + 1/p²) つき) = %.10e、差 = %.2e (相対 %.1e)、de·|1/u²-1|/p1² の見積もり %.1e' % (
        p1, de_, a1, a2, a2 - a1, abs((a2 - a1) / a1), abs(1 / (1 + a1)**2 - 1) / p1**2))

# ---------------------------------------------------------------------------------------------------- (f)
print('\n== (f) tradke の微分の途中式 ==')
pxs, pys, prs = sp.symbols('px py pr', positive=True)
pz_ = sp.sqrt(prs**2 - pxs**2 - pys**2)
dpz_code = (prs - pxs * 0 - pys * 0)  # 使わない
# dpz = (pr dpr - px dpx - py dpy)/pz  (ddpz の式)
dpx_, dpy_, dpr_ = sp.symbols('dpx dpy dpr')
lhs = sp.diff(pz_, pxs) * dpx_ + sp.diff(pz_, pys) * dpy_ + sp.diff(pz_, prs) * dpr_
rhs_ = (prs * dpr_ - pxs * dpx_ - pys * dpy_) / pz_
print('  pz = √(pr² - px² - py²): dpz - (pr dpr - px dpx - py dpy)/pz = %s' % sp.simplify(lhs - rhs_))
# 外積の微分と xpa を数値微分と比べる
rngf = np.random.default_rng(5)


def vecs(pxi, pyi, px, py, pr, sph, cph):
    pz00 = math.sqrt(pr**2 - pxi**2 - pyi**2)
    pxr0 = cph * pxi + sph * pz00
    pz0 = -sph * pxi + cph * pz00
    pz = math.sqrt(pr**2 - px**2 - py**2)
    xp = np.array([py * pz0 - pz * pyi, pz * pxr0 - px * pz0, px * pyi - py * pxr0])
    return xp, xp.copy(), pz00, pxr0, pz0, pz


def xpa_of(pxi, pyi, px, py, pr, sph, cph):
    xp = vecs(pxi, pyi, px, py, pr, sph, cph)[0]
    return np.linalg.norm(xp) / pr**2


pxi, pyi, px, py, pr = 0.012, -0.004, 0.017, 0.006, 1 - 3e-3
phir0 = 0.2
sph, cph = math.sin(phir0), math.cos(phir0)
xp, _, pz00, pxr0, pz0, pz = vecs(pxi, pyi, px, py, pr, sph, cph)
xpa = np.linalg.norm(xp) / pr**2
# 解析的な微分 (本文): 入力は (pxi, pyi, px, py, pr)。tr2 に相当する px, py の微分は単位 (入口と出口を独立変数にする)
names = ['pxi', 'pyi', 'px', 'py', 'pr']
res = []
for k_ in range(5):
    dpxi_ = 1.0 if k_ == 0 else 0.0
    dpyi_ = 1.0 if k_ == 1 else 0.0
    dpx2 = 1.0 if k_ == 2 else 0.0
    dpy2 = 1.0 if k_ == 3 else 0.0
    dpr2 = 1.0 if k_ == 4 else 0.0
    dpz = (pr * dpr2 - px * dpx2 - py * dpy2) / pz
    dpz00 = (-pxi * dpxi_ - pyi * dpyi_ + pr * dpr2) / pz00
    dpxr0 = cph * dpxi_ + sph * dpz00
    dpz0 = (-pxr0 * dpxr0 - pyi * dpyi_ + pr * dpr2) / pz0
    # コードの式 (dxpx = tr2(4)*pz0 + py*dpz0 - ddpz*pyi - pz*dpyi など)
    dxpx = dpy2 * pz0 + py * dpz0 - dpz * pyi - pz * dpyi_
    dxpy = dpz * pxr0 + pz * dpxr0 - dpx2 * pz0 - px * dpz0
    dxpz = dpx2 * pyi + px * dpyi_ - dpy2 * pxr0 - py * dpxr0
    dxpa_ana = (xp[0] * dxpx + xp[1] * dxpy + xp[2] * dxpz) / (xpa * pr**4) - (2 * xpa / pr) * dpr2
    h = 1e-6
    args = [pxi, pyi, px, py, pr]
    ap = list(args); am = list(args)
    ap[k_] += h; am[k_] -= h
    num = (xpa_of(*ap, sph, cph) - xpa_of(*am, sph, cph)) / (2 * h)
    # K64 の式 (第 1 項が /xpa/pr**2)
    dxpa_k64 = (xp[0] * dxpx + xp[1] * dxpy + xp[2] * dxpz) / xpa / pr**2 - (2 * xpa / pr) * dpr2
    res.append((names[k_], dxpa_ana, num, dxpa_k64))
for nm, ana, num, k64 in res:
    print('  d xpa/d %-3s: 本文の式 %+.8e、中心差分 %+.8e、(差 %+.1e)、K64 の式 (/xpa/pr²) %+.8e' % (nm, ana, num, ana - num, k64))
print('  (xpa = %.6e、pr = %.4f。pr ≠ 1 で K64 の式は xp 方向の微分が pr² 倍ずれる: B138 (4))' % (xpa, pr))
print('  dθ = dxpa/√(1 - xpa²) (θ = asin xpa): 差 = %s' % sp.simplify(sp.diff(sp.asin(sp.Symbol('x')), sp.Symbol('x')) - 1 / sp.sqrt(1 - sp.Symbol('x')**2)))

# ---------------------------------------------------------------------------------------------------- (g)
print('\n== (g) (13.7) の途中式 ==')
tt, pxi_, px_, pyi_, py_ = sp.symbols('t pxi px pyi py')
lin_x = pxi_ + (px_ - pxi_) * tt
lin_y = pyi_ + (py_ - pyi_) * tt
mx2 = sp.integrate(lin_x**2, (tt, 0, 1))
mxy = sp.integrate(lin_x * lin_y, (tt, 0, 1))
pxm, pym = pxi_ + px_, pyi_ + py_
print('  ∫_0^1 px(t)² dt - (pxm² + pxi² + px²)/6 = %s' % sp.simplify(mx2 - (pxm**2 + pxi_**2 + px_**2) / 6))
print('  ∫_0^1 px(t) py(t) dt - (pxm pym + pxi pyi + px py)/6 = %s' % sp.simplify(mxy - (pxm * pym + pxi_ * pyi_ + px_ * py_) / 6))
print('  ∫_0^1 px(t) dt - pxm/2 = %s' % sp.simplify(sp.integrate(lin_x, (tt, 0, 1)) - pxm / 2))
# ビーム行列の 21 成分の並び (x, px, y, py, z, δ) の下三角を行ごと。beam(i) の意味
order = ['x', 'px', 'y', 'py', 'z', 'dl']
idx = {}
n_ = 1
for i, a in enumerate(order):
    for j in range(i + 1):
        idx[(order[j], a)] = n_
        idx[(a, order[j])] = n_
        n_ += 1
print('  beam(n) の並び (行ごとの下三角): ' + '、'.join('beam(%d)=<%s %s>' % (idx[(order[j], order[i])], order[j], order[i]) for i in range(6) for j in range(i + 1) if (i, j) in [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (3, 3), (5, 1), (5, 3), (5, 5)]))
b = sp.Symbol('b')
beam = {n: sp.Symbol('beam%d' % n) for n in range(1, 22)}
B = lambda u, v: beam[idx[(u, v)]]
corr8 = B('px', 'py') - b * B('x', 'px') + b * B('y', 'py') - b**2 * B('x', 'y')                 # <(px + b y)(py - b x)> のゆらぎ部分
fort8 = beam[8] - b * (beam[2] - beam[10] + b * beam[4])
print('  正しい式 <(px + b y)(py - b x)> = %s' % sp.expand(corr8))
print('  Fortran dbeam(8) の括弧 = %s' % sp.expand(fort8))
print('  Fortran - 正しい式 = %s  (b ≠ 0 のときだけ残る: beam(10) = <py²> が beam(9) = <y py> の代わり)' % sp.expand(fort8 - corr8))
corr3 = B('px', 'px') + 2 * b * B('px', 'y') + b**2 * B('y', 'y')
fort3 = beam[3] + b * (2 * beam[5] + b * beam[6])
corr10 = B('py', 'py') - 2 * b * B('x', 'py') + b**2 * B('x', 'x')
fort10 = beam[10] + b * (-2 * beam[7] + b * beam[1])
print('  dbeam(3): Fortran - <(px + b y)²> = %s、dbeam(10): Fortran - <(py - b x)²> = %s' % (sp.simplify(fort3 - corr3), sp.simplify(fort10 - corr10)))
# 独立性の仮定のモンテカルロ
rg = np.random.default_rng(11)
n = 400000
pxv = rg.normal(0.004, 0.003, n)
pyv = 0.4 * pxv + rg.normal(-0.002, 0.002, n)
dd = rg.standard_normal(n)**2 * 1e-4 - 1e-4                  # p に依らない Δδ (平均 0 の偏った分布)
print('  Δδ が p に依らない: ⟨Δpx Δpy⟩ = %.4e、⟨px py⟩⟨Δδ²⟩ = %.4e、⟨Δpx Δδ⟩ = %.4e、⟨px⟩⟨Δδ²⟩ = %.4e' % (
    np.mean(pxv * pyv * dd * dd), np.mean(pxv * pyv) * np.mean(dd * dd), np.mean(pxv * dd * dd), np.mean(pxv) * np.mean(dd * dd)))
dd2 = dd * (1 + 300 * pxv)                                        # Δδ が px に依る (仮定を破る)
print('  Δδ が px に依る (仮定を破る例): ⟨Δpx²⟩ = %.4e、⟨px²⟩⟨Δδ²⟩ = %.4e (相対 %.1f%% 違う)' % (
    np.mean(pxv**2 * dd2**2), np.mean(pxv**2) * np.mean(dd2**2), 100 * (np.mean(pxv**2 * dd2**2) / (np.mean(pxv**2) * np.mean(dd2**2)) - 1)))

# ---------------------------------------------------------------------------------------------------- (h)
print('\n== (h) 薄い偏向磁石の FLUC: ⟨Δpx Δδ⟩ = ⟨px⟩ Var(δ)、⟨Δpx²⟩ = ⟨px²⟩ Var(δ)、⟨Δpx Δpy⟩ = ⟨px py⟩ Var(δ) ==')
thin = 1.0e-3
Lth = 0.0006
src = 'MOMENTUM = 3GEV;\nBEND B=(L=%.15g ANGLE=%.15g);\nMARK ST=();\nLINE ONE=(ST B);\nFFS;\nUSE ONE;\n' % (Lth, thin)
ffs, _ = run_sad_string(src, out=io.StringIO())
els = ffs.lat.elements
gam = 3e9 / mc2
Nth = 5 / (2 * math.sqrt(3)) * alpha * gam * thin
uc_th = 1.5 * (re / alpha) * gam**3 / (Lth / thin) * mc2 / 3e9
print('  θ = %.1e rad、L = %.1e m、N_γ = %.4f、u_c/(p0 c) = %.3e、理論 Var(δ) = N (11/27) u_c² = %.4e' % (thin, Lth, Nth, uc_th, Nth * 11 / 27 * uc_th**2))
npart = 400000
for px0, py0 in ((0.0, 0.0), (0.2, 0.1), (-0.15, 0.2)):
    a1 = np.zeros((npart, 6))
    a1[:, 1] = px0
    a1[:, 3] = py0
    rf = np.asarray(track6.track_particles(els, a1.copy(), 1, ffs.momentum(), 0.0, False, True, True, seed=21))
    dxp = rf[:, 1] - rf[:, 1].mean()
    dyp = rf[:, 3] - rf[:, 3].mean()
    ddl = rf[:, 5] - rf[:, 5].mean()
    var = (ddl**2).mean()
    print('  px0 = %+.2f、py0 = %+.2f: Var(δ) = %.4e、⟨Δpx Δδ⟩/Var = %+.4f (px0 = %+.2f)、⟨Δpx²⟩/Var = %.5f (px0² = %.5f)、⟨Δpx Δpy⟩/Var = %+.5f (px0 py0 = %+.5f)' % (
        px0, py0, var, (dxp * ddl).mean() / var, px0, (dxp**2).mean() / var, px0**2, (dxp * dyp).mean() / var, px0 * py0))
print('  (偏向磁石の薄い 1 台で、軌道の px は入口 px0 から出口まで θ 程度 (1e-3) しか変わらないので ⟨px⟩ ≈ px0 とみなせる。標本誤差は Var(δ) の分布の裾で 2〜3%。)')
