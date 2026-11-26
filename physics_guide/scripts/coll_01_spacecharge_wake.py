# -*- coding: utf-8 -*-
"""集団効果 (17_collective.md): 空間電荷 (WSPAC) のキックの係数と、ウェイクのキックの規格化。

(a) sympy: 静止系の電場 E' を実験室系に直し、E と B の力の和が E'/γ になること。キックの係数 1/(β²γ³)。
(b) pysad の WSPAC (spch.WSpace.kick) のキックを、2 次元ポテンシャルの数値積分 (16 章) から作った式と比べる: 横のキック (丸い・楕円)、
    δ のキック (ポテンシャル U(x,y) の z 微分)、小振幅の線形の傾き 2 f/(σx(σx+σy))。
(c) ウェイク: ガウス分布のバンチ + 指数関数のウェイク W(s) = W0 exp(-s/l) (s ≥ 0) の平均のエネルギー損失が、解析式
    k = (1/2) exp(σ²/l²) erfc(σ/l) (自己項の 1/2 を含む) と合うか。pysad の txwake の直訳 (wake.WakeKick) で追跡する。
    横ウェイクは、バンチ全体を x0 ずらしたときの平均のキック。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())

# ---------------------------------------------------------------- (a)
import sympy as sp
beta, gam_s, Ep, e_, m_, c_ = sp.symbols('beta gamma Ep e m c', positive=True)
gm = 1 / sp.sqrt(1 - beta**2)
E_lab = gm * Ep                       # 実験室系の横の電場 (E⊥ = γ E'⊥)
B_lab = beta * gm * Ep / c_           # B⊥ = (β/c) γ E'⊥
force = e_ * (E_lab - beta * c_ * B_lab)
print('== (a) 実験室系の力 F = e(E - v B) = e E\'/γ ==')
print('F γ / (e E\') - 1 =', sp.simplify(force * gm / (e_ * Ep) - 1))
# Δp/p0: Δp = F Δs/(βc)、p0 = γ m β c。E' = 2 r_e... は r_e = e²/(4πε0 m c²) とすれば、Δp/p0 = (F/e) e Δs/(β c γ m β c)
dpp = (force / e_) * 1 / (beta * c_) / (gm * m_ * beta * c_)
print("Δp/p0 (= F/e/(βc)/(γ m βc)) × γ² β² m c² / E' - 1 =", sp.simplify(dpp * gm**2 * beta**2 * m_ * c_**2 / Ep - 1), "(e E' Δs/(γ² β² m c²) の形。E' ∝ λ' = λ/γ を入れて全体 1/(β²γ³))")
print('→ λ\' = λ/γ、E\' = 2 r_e (mc²/e) λ\' (1 - e^{-r²/2σ²})/r から Δp/p0 = 2 r_e λ Δs (1 - e^{-r²/2σ²}) / (r β² γ³)')

# ---------------------------------------------------------------- (b)
from pysad import spch
from scipy.integrate import quad
re = 2.8179403262e-15
mc2 = 0.51099895e6


def U_pot(x, y, sx2, sy2):
    f = lambda t: (1 - math.exp(-x * x / (2 * sx2 + t) - y * y / (2 * sy2 + t))) / math.sqrt((2 * sx2 + t) * (2 * sy2 + t))
    v = 0.0
    a = 0.0
    for b in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, np.inf):
        v += quad(f, a, b, epsabs=0, epsrel=1e-12, limit=200)[0]
        a = b
    return v


print('\n== (b) WSPAC のキック (pysad v2: spch.WSpace) と 2 次元ポテンシャルの数値積分 ==')
Ek = 20e6                              # 運動エネルギー [eV]
g0 = 1 + Ek / mc2
p0 = math.sqrt(g0**2 - 1)
h0 = g0
b0 = p0 / h0
N = 1e10
sigz = 5e-3
al = 0.3
print('20 MeV の電子: γ = %.4f、β = %.6f、1/(β²γ³) = %.6e、N = %.0e、σz = %.0e m、長さ %.1f m' % (g0, b0, 1 / (b0**2 * g0**3), N, sigz, al))
for label, (sxm, sym_) in (('丸い (σ = 1 mm)', (1e-3, 1e-3)), ('楕円 (σx = 1.5 mm、σy = 0.5 mm)', (1.5e-3, 0.5e-3))):
    Sig = np.zeros((6, 6))
    Sig[0, 0], Sig[1, 1] = sxm**2, 1e-10
    Sig[2, 2], Sig[3, 3] = sym_**2, 1e-10
    Sig[4, 4], Sig[5, 5] = sigz**2, 1e-6
    ws = spch.WSpace([Sig], np.zeros((1, 6)), re * N)
    pts = [(1.2e-3, 0.4e-3, 3e-3), (-0.5e-3, 1.1e-3, -6e-3), (2.5e-3, -0.8e-3, 0.0)]
    S = np.zeros((len(pts), 16))
    for i, (xx, yy, zz) in enumerate(pts):
        S[i, 0], S[i, 2], S[i, 4] = xx, yy, zz
    S0 = S.copy()
    ws.kick(S, len(pts), 0, al, p0, h0)
    print('\n%s' % label)
    for i, (xx, yy, zz) in enumerate(pts):
        # 独立な式: Δp/p0 = A(z) ∇U、A(z) = N r_e L exp(-z²/2σz²)/(β²γ³ √(2π) σz)。丸い・楕円とも U は 2 次元ポテンシャルの数値積分 (16 章 (16.6))
        A = N * re * al * math.exp(-0.5 * zz**2 / sigz**2) / (b0**2 * g0**3 * math.sqrt(2 * math.pi) * sigz)
        h = 1e-5 * max(abs(xx), abs(yy))
        sx2, sy2 = sxm**2, sym_**2
        dUx = (U_pot(xx + h, yy, sx2, sy2) - U_pot(xx - h, yy, sx2, sy2)) / (2 * h)
        dUy = (U_pot(xx, yy + h, sx2, sy2) - U_pot(xx, yy - h, sx2, sy2)) / (2 * h)
        # Δδ = (z/σz²) A (U_ref - U)。K64 の基準 (rwall = 0): -U_ref = ln((1 + r)/2)、r = σ小/σ大 (twspu の log(0.5 + 0.5 r)。丸いビームでは 0、軸上で
        # 0 でない。gs で確認)。S の px は p/p0 の内部変数ではなく
        # px/(1+δ) (追跡の内部変数) だが、δ0 = 0 の入力なので、Δpx は新しい δ で割る前の値と δ の差 (~1e-12) 以内で等しい
        u_ref = -math.log(0.5 + 0.5 * min(sxm, sym_) / max(sxm, sym_))
        dd_ex = (zz / sigz**2) * A * (u_ref - U_pot(xx, yy, sx2, sy2))
        print('  (x, y, z) = (%+.1e, %+.1e, %+.1e): Δpx = %+.6e (式 %+.6e、差 %.1e)、Δpy = %+.6e (式 %+.6e、差 %.1e)、Δδ = %+.5e (式 %+.5e、差 %.1e)' % (
            xx, yy, zz, S[i, 1], A * dUx, S[i, 1] - A * dUx, S[i, 3], A * dUy, S[i, 3] - A * dUy, S[i, 5], dd_ex, S[i, 5] - dd_ex))
    kx, ky = ws.linear_k(0, al, p0)
    fl = 2 * re * N * al / (b0**2 * g0**3 * math.sqrt(2 * math.pi) / 1 * sigz)
    print('  線形の傾き: linear_k = (%.6e, %.6e)、式 2 f/(σx(σx+σy)) = %.6e、2 f/(σy(σx+σy)) = %.6e' % (
        kx, ky, 2 * re * N * al / (b0**2 * g0**3 * math.sqrt(2 * math.pi) * sigz) / (sxm * (sxm + sym_)),
        2 * re * N * al / (b0**2 * g0**3 * math.sqrt(2 * math.pi) * sigz) / (sym_ * (sxm + sym_))))

# ---------------------------------------------------------------- (c)
print('\n== (c) ウェイク: ガウス分布 + 指数関数 W(s) = W0 exp(-s/l) ==')
from pysad import wake as W
ELEMCH = 1.602176634e-19
p0w = math.sqrt((1.7e9 / mc2) ** 2 + 0.0)    # 1.7 GeV (運動量)
g0w = math.sqrt(1 + p0w**2)
npart = 20000
rng = np.random.default_rng(5)
zz = rng.normal(0.0, sigz, npart)
W0 = 1e12                                    # V/C
Nb = 1e10
print('1.7 GeV: p0 = βγ = %.2f、粒子 %d、σz = %.1f mm、W0 = %.0e V/C、N = %.0e' % (p0w, npart, sigz * 1e3, W0, Nb))
for ell in (2e-3, 5e-3, 20e-3):
    # 表: z (= 源 - 試験粒子、先行する源が正) の昇順。z < 0 は 0 (非因果)。表の最後の行は追跡で使われない
    zt = np.concatenate([[-0.05, -1e-5], np.linspace(0.0, 0.2, 4001), [0.21]])
    wt = np.where(zt >= 0.0, W0 * np.exp(-zt / ell), 0.0)
    tab = np.stack([zt, wt], axis=1)
    S = np.zeros((npart, 16))
    S[:, 4] = zz
    fw = W.kick_factor(-1, Nb, 1, npart)
    wk = W.WakeKick(p0w, g0w, fw, 1)
    wk.apply(S, npart, tab, None, 0.0, 0.0, True)       # 要素の前の半分
    wk.apply(S, npart, tab, None, 0.0, 0.0, False)      # 要素の後の半分
    mean_d = S[:, 5].mean()
    k = 0.5 * math.exp(sigz**2 / ell**2) * math.erfc(sigz / ell)
    exp_d = -(ELEMCH * Nb * W0 * k / mc2) / p0w
    print('l = %4.0f mm: <Δδ> = %+.6e、解析式 -(e N W0 k/mc²)/p0 = %+.6e (k = %.5f)、比 %.4f' % (ell * 1e3, mean_d, exp_d, k, mean_d / exp_d))

print('\n-- 横ウェイク: バンチ全体を x0 = 1 mm ずらす (W_T(s) = W0t exp(-s/l)、V/C/m) --')
W0t = 1e13
ell = 5e-3
zt = np.concatenate([[-0.05, -1e-5], np.linspace(0.0, 0.2, 4001), [0.21]])
wt = np.where(zt >= 0.0, W0t * np.exp(-zt / ell), 0.0)
tab = np.stack([zt, wt], axis=1)
S = np.zeros((npart, 16))
S[:, 4] = zz
S[:, 0] = 1e-3
fw = W.kick_factor(-1, Nb, 1, npart)
wk = W.WakeKick(p0w, g0w, fw, 1)
wk.apply(S, npart, None, tab, 0.0, 0.0, True)
wk.apply(S, npart, None, tab, 0.0, 0.0, False)
k = 0.5 * math.exp(sigz**2 / ell**2) * math.erfc(sigz / ell)
exp_p = (ELEMCH * Nb * W0t * 1e-3 * k / mc2) / p0w
print('<Δpx> = %+.6e、解析式 +(e N W0t x0 k/mc²)/p0 = %+.6e (比 %.4f)。Δδ = %.1e (縦ウェイクなし)' % (S[:, 1].mean(), exp_p, S[:, 1].mean() / exp_p, S[:, 5].mean()))

print('\n-- 1 つの粒子の自己項: 点電荷 1 個 (n = 1) の損失 = (1/2) e N... は k → 1/2 (σ → 0)。l = 5 mm で σ → 0 の k =', 0.5 * math.exp(0.0) * math.erfc(0.0))
