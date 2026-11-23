# -*- coding: utf-8 -*-
"""多極の端部の非線形フリンジ (07_thin_multipole.md §9、式 (7.20)〜(7.27)) を確かめる。

2m 極 (m = k + 1: 四極 m = 2、6 極 m = 3、8 極 m = 4、…) のスカラーポテンシャル Φ = κ(s) Im w^m/m! - κ''(s)|w|^2 Im w^m/(4 (m+1)!) (w = x + i y) の場で、
κ = κ0 Θ(s) (硬い端) の端を通る粒子の位置のずれは、K の 1 次で Δw = -aki w (w^m/(m+1) - w̄^m)、aki = κ0/(4 m! p) (式 (7.22))。K64 の ttfrin の ca と同じ。
(a) sympy: Φ は調和 (∇²Φ + ∂s²Φ = κ'' ... の相殺)、複素形の微分で Δw の式を出す。
(b) pysad の ttfrin (track6._ttfrin_n、m = 2〜6) の位置のずれと Δw の式の比較 (κ0 を 1/2 にすると差が 1/4: K の 2 次)。
(c) 運動量は点変換の逆転置 (正準変換)、z は -Δr·p。シンプレクティック性。
(d) Maxwell の場 (tanh の端) をローレンツ力で積分した結果との比較 (6 極・8 極)。
"""
import math
import warnings
warnings.filterwarnings("ignore")
from types import SimpleNamespace
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

print('== (a) sympy: 場と端の位置のずれ ==')
x, y = sp.symbols('x y', real=True)
kap, kap2 = sp.symbols('kappa kappa2', real=True)       # kappa2 = κ''(s)
w = x + sp.I * y
wb = x - sp.I * y
for m in (2, 3, 4):
    Phi_xy = kap / sp.factorial(m) * sp.im(sp.expand(w**m)) - kap2 * sp.expand(w * wb) * sp.im(sp.expand(w**m)) / (4 * sp.factorial(m + 1))
    lap = sp.simplify(sp.diff(Phi_xy, x, 2) + sp.diff(Phi_xy, y, 2))
    print('  m = %d: ∇⊥²Φ = %s、∂s²Φ = κ\'\' Im w^m/m! を足すと 0 (∇⊥²Φ + κ\'\' Im w^m/m! = %s)' % (
        m, sp.factor(lap), sp.simplify(lap + kap2 * sp.im(sp.expand(w**m)) / sp.factorial(m))))
    # 力の κ'' の部分: Fx = -∂Φ/∂y、Fy = +∂Φ/∂x (pz ≈ p)。κ'' = κ0 δ'(s) → 位置のずれ Δ(x + i y) = (κ0/p)(-∂y + i ∂x) の係数
    G = -kap2 * sp.expand(w * wb) * sp.im(sp.expand(w**m)) / (4 * sp.factorial(m + 1))
    dx_ = -sp.diff(G, y)
    dy_ = sp.diff(G, x)
    dw = sp.expand(dx_ + sp.I * dy_)
    aki = 1 / (4 * sp.factorial(m))
    claim = sp.expand(-aki * kap2 * w * (w**m / (m + 1) - wb**m))
    print('       Δ(x + i y) = (Fx + i Fy の κ\'\' の係数) と -aki κ\'\' w (w^m/(m+1) - w̄^m) (aki = 1/(4 m!)) の差: %s' % sp.simplify(sp.expand(dw - claim)))
print('  (κ(s) = κ0 Θ(s) で κ\'\' = κ0 δ\'(s)。力が δ\'(s) なので、運動量は δ(s) 型にスパイクし、位置は端で階段状にずれる: Δw = (1/p) × (κ\'\' の係数) × κ0。')
print('   運動量の変化 (K の 1 次) は ∫δ\'(s) g ds = -g\'(0) から、軌道に沿った微分で決まる (四極の章 (4.26) と同じ手順)。')

print('\n== (b) pysad の ttfrin (入口) の位置のずれと Δw の式 ==')


def ttfrin_pysad(c, nord, ak, al):
    q = SimpleNamespace(c=np.array([c[0], c[1], c[2], c[3], c[4], c[5]], float))
    t6._ttfrin_n(q, nord, ak, al)
    return q.c.copy()


rng = np.random.default_rng(11)
cpt = np.array([0.004, 0.0003, -0.003, 0.0002, 0.0, 0.0])
for m in (2, 3, 4, 5, 6):
    out = []
    k_base = 1e-2 * 4 * math.factorial(m) / 5e-3**m          # aki r^m = 0.01 になる強さ (r = 5 mm)
    for kappa0 in (k_base, k_base / 2):
        o = ttfrin_pysad(cpt, 2 * m, kappa0, 1.0)
        dw_code = complex(o[0] - cpt[0], o[2] - cpt[2])
        wv = complex(cpt[0], cpt[2])
        aki = kappa0 / (4 * math.factorial(m))
        dw_form = -aki * wv * (wv**m / (m + 1) - wv.conjugate()**m)
        out.append((kappa0, dw_code, dw_form))
    e1 = abs(out[0][1] - out[0][2])
    e2 = abs(out[1][1] - out[1][2])
    print('  m = %d (nord = %d): κ0 = %g: |Δw| = %.4e、式 %.4e、差 %.2e。κ0 を 1/2 にすると差 %.2e (比 %.1f。K の 2 次なら 4)' % (
        m, 2 * m, out[0][0], abs(out[0][1]), abs(out[0][2]), e1, e2, e1 / e2 if e2 > 0 else float('nan')))

print('\n== (c) 点変換 + 逆転置の運動量 + z。シンプレクティック性 (m = 3、6 極、κ0 = 500 /m^3、δ = 0.005) ==')
m = 3
kappa0 = 500.0
c1 = np.array([0.004, 0.0003, -0.003, 0.0002, 0.0, 0.005])


def pt(xy):
    """位置の点変換 (運動量によらない)"""
    o = ttfrin_pysad(np.array([xy[0], 0.0, xy[1], 0.0, 0.0, c1[5]]), 2 * m, kappa0, 1.0)
    return np.array([o[0], o[2]])


h = 1e-7
Jp = np.zeros((2, 2))
for j in range(2):
    a_, b_ = np.array([c1[0], c1[2]]), np.array([c1[0], c1[2]])
    a_[j] += h
    b_[j] -= h
    Jp[:, j] = (pt(a_) - pt(b_)) / (2 * h)
o = ttfrin_pysad(c1, 2 * m, kappa0, 1.0)
pnew = np.array([o[1], o[3]])
pold = np.array([c1[1], c1[3]])
pred = np.linalg.solve(Jp.T, pold)
print('  運動量: pysad の ttfrin の p_new と、点変換のヤコビアンの逆転置 J^{-T} p_old の差 %.2e' % np.abs(pnew - pred).max())
dr = np.array([o[0] - c1[0], o[2] - c1[2]])
print('  z: Δz = %.6e、-Δr·p_new = %.6e (差 %.2e。q 変数なので p はすでに 1/(1+δ) で規格化)' % (o[4] - c1[4], -dr @ pnew, abs(o[4] - c1[4] + dr @ pnew)))

Jm = np.zeros((6, 6))
for i in range(3):
    Jm[2 * i, 2 * i + 1] = 1.0
    Jm[2 * i + 1, 2 * i] = -1.0


def canon(c):
    """canonical → q 変数 → ttfrin → canonical"""
    pr = 1.0 + c[5]
    q = np.array([c[0], c[1] / pr, c[2], c[3] / pr, c[4], c[5]])
    o = ttfrin_pysad(q, 2 * m, kappa0, 1.0)
    return np.array([o[0], o[1] * pr, o[2], o[3] * pr, o[4], o[5]])


Jf = np.zeros((6, 6))
for j in range(6):
    cp, cm = c1.copy(), c1.copy()
    cp[j] += 1e-6
    cm[j] -= 1e-6
    Jf[:, j] = (canon(cp) - canon(cm)) / 2e-6
print('  6x6 ヤコビアン (正準運動量、δ ≠ 0 を含む) のシンプレクティック性 |M^T J M - J| = %.1e (中心差分 1e-6 の精度)' % np.abs(Jf.T @ Jm @ Jf - Jm).max())

print('\n== (d) Maxwell の場の端 (tanh) をローレンツ力で積分した結果との比較 (p = 1、入口、S = 0.06 m) ==')
xs, ys = sp.symbols('xs ys', real=True)


def make_fields(m):
    """Φ = κ Φ0 - κ'' ψ (Φ0 = Im w^m/m!、ψ = |w|^2 Im w^m/(4 (m+1)!))。b = (∂xΦ, ∂yΦ, ∂sΦ)。b_s は κ' Φ0 だけ (κ の 3 階微分の項は 3 次の展開の外)。
    戻り値: 場の関数 (bx, by, bs = f(x, y, κ, κ', κ'')) と、線形な場だけ (bx = κ ∂xΦ0、by = κ ∂yΦ0、bs = 0) の関数"""
    ww = xs + sp.I * ys
    Phi0 = sp.im(sp.expand(ww**m)) / sp.factorial(m)
    psi = sp.expand(ww * sp.conjugate(ww)) * sp.im(sp.expand(ww**m)) / (4 * sp.factorial(m + 1))
    f = sp.lambdify((xs, ys), [sp.diff(Phi0, xs), sp.diff(Phi0, ys), sp.diff(psi, xs), sp.diff(psi, ys), Phi0], 'math')

    def full(x_, y_, k0, k1, k2):
        p0x, p0y, psx, psy, p0 = f(x_, y_)
        return k0 * p0x - k2 * psx, k0 * p0y - k2 * psy, k1 * p0

    def lin(x_, y_, k0, k1, k2):
        p0x, p0y, psx, psy, p0 = f(x_, y_)
        return k0 * p0x, k0 * p0y, 0.0
    return full, lin


def rhs_factory(bf, kappa0, d, sgn):
    def rhs(s, v):
        x_, px_, y_, py_, z_ = v
        t = math.tanh(s / d)
        sech2 = 1 - t * t
        k0 = kappa0 * (1 + sgn * t) / 2
        k1 = sgn * kappa0 * sech2 / (2 * d)
        k2 = -sgn * kappa0 * t * sech2 / d**2
        bx, by, bs = bf(x_, y_, k0, k1, k2)
        pz = math.sqrt(max(1.0 - px_**2 - py_**2, 1e-12))
        return [px_ / pz, (py_ * bs - pz * by) / pz, py_ / pz, (pz * bx - px_ * bs) / pz, 1.0 - 1.0 / pz]
    return rhs


def kin_body_mult(c, S, kappa0, m):
    k = m - 1

    def rhs(s, v):
        x_, px_, y_, py_, z_ = v
        g = kappa0 / math.factorial(k) * complex(x_, -y_)**k
        pz = math.sqrt(1.0 - px_**2 - py_**2)
        return [px_ / pz, -g.real, py_ / pz, -g.imag, 1.0 - 1.0 / pz]
    return solve_ivp(rhs, (0, S), list(c), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1]


def kin_drift(c, S):
    x_, px_, y_, py_, z_ = c
    pz = math.sqrt(1.0 - px_**2 - py_**2)
    return np.array([x_ + px_ / pz * S, px_, y_ + py_ / pz * S, py_, z_ + (1.0 - 1.0 / pz) * S])


S = 0.06
print('   Δ(ODE) = ODE(Maxwell の場) - ODE(κ(s) Im w^m/m! だけの場)、Δ(モデル) = (ドリフト → ttfrin → 一様な多極) - (ドリフト → 一様な多極)')
for m, kappa0 in ((3, 400.0), (4, 4.0e4)):
    bf_m, bf_l = make_fields(m)
    print('  m = %d (%s、κ0 = %g /m^%d):' % (m, '6 極' if m == 3 else '8 極', kappa0, m))
    print('    r [mm]  d [mm] | Δ(ODE)    Δ(モデル)   差/Δ(モデル)')
    for scale in (1.0, 0.5):
        for d in (0.0005, 0.001, 0.002):
            c0 = np.array([0.004, 0.0003, -0.003, 0.0002, 0.0]) * scale
            kw = dict(method='DOP853', rtol=1e-13, atol=1e-17)
            ex_m = solve_ivp(rhs_factory(bf_m, kappa0, d, +1.0), (-S, S), c0, **kw).y[:, -1]
            ex_l = solve_ivp(rhs_factory(bf_l, kappa0, d, +1.0), (-S, S), c0, **kw).y[:, -1]
            cc = kin_drift(c0, S)
            mod_n = kin_body_mult(cc, S, kappa0, m)
            of = ttfrin_pysad(np.array([cc[0], cc[1], cc[2], cc[3], cc[4], 0.0]), 2 * m, kappa0, 1.0)
            mod_f = kin_body_mult(of[:5], S, kappa0, m)
            dex, dmod = ex_m - ex_l, mod_f - mod_n
            print('    %5.1f   %5.1f  | %.3e  %.3e   %.4f' % (4.0 * scale, d * 1e3, np.abs(dex).max(), np.abs(dmod).max(), np.abs(dex - dmod).max() / np.abs(dmod).max()))
print('   (r = 振幅 x。r を半分にしても「差/Δ(モデル)」が変わらなければ、この次数の効果として合っている。残りは K の 2 次以上と端の幅 d の補正)')
