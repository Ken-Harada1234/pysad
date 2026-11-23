# -*- coding: utf-8 -*-
r"""1 周の写像の 2 次の摂動 (18_nonlinear_bridge.md §3.7〜3.8) の式を、薄い 6 極・8 極の列の写像の追跡と比べる。

模型 (pysad を使わない自己完結の確認): 正規化座標 (x^, p^, y^, q^) で、薄い 6 極・8 極の kick と線形の回転を交互に並べた 1 周の写像。
  回転: ζ → e^{iμ} ζ (x^' = x^ cos μ + p^ sin μ、p^' = -x^ sin μ + p^ cos μ)。
  6 極 (K2L): Δp_x = -(K2L/2)(x² - y²)、Δp_y = K2L x y。x = √β_x x^、p^ = √β p (Δp^ = √β Δp)。
  8 極 (K3L): Δp_x = -(K3L/6)(x³ - 3xy²)、Δp_y = (K3L/6)(3x²y - y³)。
式 (本文の導出どおり):
  K_diag = [h4 の対角 = Magnus の対の項 (1/2) Σ_{l>e} {H3(l), H3(e)} + 8 極の H4] - (1/2){g, Rg} の対角、
  g_q = F_q/(e^{i m·μ} - 1) (F = -h3、m = (j-k, l-m))、(Rg)_q = e^{i m·μ} g_q。
  ∂νx/∂Jx = 4 K_{22 00}/π、∂νx/∂Jy = 2 K_{1111}/π、∂νy/∂Jy = 4 K_{0022}/π。
追跡: 1024〜4096 周、±位相の対、NAFF (Hann 窓)、J の格子の最小二乗フィット (2 次)。
"""
import itertools
import numpy as np
from scipy.optimize import minimize_scalar
import _common as C

C.start(__file__)

TWO_PI = 2 * np.pi


# ---------- 多項式 (ζ, ζbar, v, vbar の単項式の辞書) ----------
def pmul(A, B):
    R = {}
    for ka, ca in A.items():
        for kb, cb in B.items():
            k = tuple(i + j for i, j in zip(ka, kb))
            R[k] = R.get(k, 0) + ca * cb
    return R


def padd(A, B, s=1.0):
    R = dict(A)
    for k, c in B.items():
        R[k] = R.get(k, 0) + s * c
    return R


def pscale(A, s):
    return {k: s * c for k, c in A.items()}


def pbracket(A, B):
    """{A, B} = 2i Σ (∂A/∂ζ ∂B/∂ζbar - ∂A/∂ζbar ∂B/∂ζ + (v も同様))"""
    R = {}
    for (a, b, c, d), ca in A.items():
        for (a2, b2, c2, d2), cb in B.items():
            fx = a * b2 - b * a2
            fy = c * d2 - d * c2
            if fx:
                k = (a + a2 - 1, b + b2 - 1, c + c2, d + d2)
                R[k] = R.get(k, 0) + 2j * fx * ca * cb
            if fy:
                k = (a + a2, b + b2, c + c2 - 1, d + d2 - 1)
                R[k] = R.get(k, 0) + 2j * fy * ca * cb
    return R


def kick_polys(bx, by, mux, muy, K2L, K3L):
    """位置の位相進み (mux, muy) [rad]、β での H3 (6 極) と H4 (8 極) の ζ 単項式。x = (w ζ + wbar ζbar)/2、w = √β e^{iμ}"""
    wx, wy = np.sqrt(bx) * np.exp(1j * mux), np.sqrt(by) * np.exp(1j * muy)
    X = {(1, 0, 0, 0): wx / 2, (0, 1, 0, 0): np.conj(wx) / 2}
    Y = {(0, 0, 1, 0): wy / 2, (0, 0, 0, 1): np.conj(wy) / 2}
    X2, Y2 = pmul(X, X), pmul(Y, Y)
    X3, X4, Y4 = pmul(X2, X), pmul(X2, X2), pmul(Y2, Y2)
    H3 = padd(pscale(X3, K2L / 6), pscale(pmul(X, Y2), -K2L / 2))
    H4 = padd(padd(pscale(X4, K3L / 24), pscale(pmul(X2, Y2), -K3L / 4)), pscale(Y4, K3L / 24))
    return H3, H4


def adts_formula(kicks, nux, nuy, verbose=False):
    """kicks: [(bx, by, mux, muy, K2L, K3L)] (位相進みは 1 周の起点からの値 [rad])。戻り値: (xx, xy, yy) = (∂νx/∂Jx, ∂νx/∂Jy, ∂νy/∂Jy)"""
    H3s, H4s = [], []
    for kk in kicks:
        a, b = kick_polys(*kk)
        H3s.append(a)
        H4s.append(b)
    # h3 = Σ H3、h4 = (1/2) Σ_{l>e} {H3(l), H3(e)} + Σ H4
    h3 = {}
    for H in H3s:
        h3 = padd(h3, H)
    h4 = {}
    for H in H4s:
        h4 = padd(h4, H)
    for l in range(len(H3s)):
        for e in range(l):
            h4 = padd(h4, pscale(pbracket(H3s[l], H3s[e]), 0.5))
    # 正規形: F = -h3、g_q = F_q/(e^{i m.μ} - 1)、Rg
    g, Rg = {}, {}
    for (a, b, c, d), cf in h3.items():
        m = (a - b) * TWO_PI * nux + (c - d) * TWO_PI * nuy
        ph = np.exp(1j * m)
        g[(a, b, c, d)] = -cf / (ph - 1)
        Rg[(a, b, c, d)] = -cf / (ph - 1) * ph
    gRg = pbracket(g, Rg)
    K = padd(h4, gRg, -0.5)
    xx = K.get((2, 2, 0, 0), 0)
    xy = K.get((1, 1, 1, 1), 0)
    yy = K.get((0, 0, 2, 2), 0)
    if verbose:
        print('    K の対角: ζ²ζbar² %s、ζζbar v vbar %s、v²vbar² %s' % (xx, xy, yy))
    return np.array([4 * xx / np.pi, 2 * xy / np.pi, 4 * yy / np.pi]).real, max(abs(xx.imag), abs(xy.imag), abs(yy.imag))


# ---------- 追跡 ----------
def one_turn(state, kicks, nux, nuy):
    x, p, y, q = state
    mx0 = my0 = 0.0
    for (bx, by, mux, muy, K2L, K3L) in kicks:
        dmx, dmy = mux - mx0, muy - my0
        mx0, my0 = mux, muy
        cx, sx, cy, sy = np.cos(dmx), np.sin(dmx), np.cos(dmy), np.sin(dmy)
        x, p = x * cx + p * sx, -x * sx + p * cx
        y, q = y * cy + q * sy, -y * sy + q * cy
        X, Y = np.sqrt(bx) * x, np.sqrt(by) * y
        dpx = -K2L / 2 * (X**2 - Y**2) - K3L / 6 * (X**3 - 3 * X * Y**2)
        dpy = K2L * X * Y + K3L / 6 * (3 * X**2 * Y - Y**3)
        p = p + np.sqrt(bx) * dpx
        q = q + np.sqrt(by) * dpy
    dmx, dmy = TWO_PI * nux - mx0, TWO_PI * nuy - my0
    cx, sx, cy, sy = np.cos(dmx), np.sin(dmx), np.cos(dmy), np.sin(dmy)
    x, p = x * cx + p * sx, -x * sx + p * cx
    y, q = y * cy + q * sy, -y * sy + q * cy
    return x, p, y, q


def naff(sig, nu0=None):
    """Hann 窓つきの内積の最大の周波数 (符号つき)。sig = x^ - i p^ (位相が +2πν n で進む)"""
    N = len(sig)
    n = np.arange(N)
    w = 1 + np.cos(np.pi * (n - N / 2) / (N / 2))
    f = lambda nu: -abs(np.sum(w * sig * np.exp(-1j * TWO_PI * nu * n)))
    sp_ = np.fft.fft(sig * w)
    k = int(np.argmax(abs(sp_)))
    c = k / N if k < N / 2 else (k - N) / N
    r = minimize_scalar(f, bracket=(c - 1.0 / N, c, c + 1.0 / N), tol=1e-14)
    return r.x


def track_tunes(kicks, nux, nuy, Jx, Jy, nturn=2048):
    """±位相の対で追跡し、(νx, νy) の平均を返す"""
    res = []
    for sgn in (+1, -1):
        x = np.array([sgn * np.sqrt(2 * Jx)])
        p = np.array([0.0])
        y = np.array([sgn * np.sqrt(2 * Jy)])
        q = np.array([0.0])
        st = (x, p, y, q)
        X, Pp, Yy, Qq = [], [], [], []
        for _ in range(nturn):
            st = one_turn(st, kicks, nux, nuy)
            X.append(st[0][0]); Pp.append(st[1][0]); Yy.append(st[2][0]); Qq.append(st[3][0])
        X, Pp, Yy, Qq = map(np.array, (X, Pp, Yy, Qq))
        a = naff(X - 1j * Pp) if Jx > 0 else nux
        b = naff(Yy - 1j * Qq) if Jy > 0 else nuy
        res.append((a, b))
    return np.mean(res, axis=0)


def fit_adts(kicks, nux, nuy, J0, nturn=2048):
    """J の格子 {0..3}·J0 (Jx, Jy) で ν(Jx, Jy) を 2 次の多項式フィットし、(∂νx/∂Jx, ∂νx/∂Jy, ∂νy/∂Jy) を返す"""
    pts = [(i, j) for i in range(4) for j in range(4) if (i, j) != (0, 0)]
    rows, nx_, ny_, rx_, ry_ = [], [], [], [], []
    for i, j in pts:
        Jx, Jy = i * J0, j * J0
        a, b = track_tunes(kicks, nux, nuy, Jx, Jy, nturn)
        r = [Jx, Jy, Jx**2, Jx * Jy, Jy**2]
        if i > 0:
            rx_.append(r); nx_.append(a - nux)
        if j > 0:
            ry_.append(r); ny_.append(b - nuy)
    cx = np.linalg.lstsq(np.array(rx_), np.array(nx_), rcond=None)[0]
    cy = np.linalg.lstsq(np.array(ry_), np.array(ny_), rcond=None)[0]
    return np.array([cx[0], cx[1], cy[1]]), cy[0]


nux, nuy = 0.2347, 0.3141
phase = [(0.05, 0.03), (0.40, 0.22), (0.71, 0.60)]
betas = [(6.0, 2.0), (2.0, 7.0), (9.0, 3.0)]
K2 = [0.8, -1.1, 0.9]
print('チューン (νx, νy) = (%.4f, %.4f)。3νx = %.4f、νx ± 2νy = %.4f、%.4f (整数から離れている)' % (nux, nuy, 3 * nux, nux + 2 * nuy, nux - 2 * nuy))


def make(lam, octupole=0.0):
    ks = []
    for (px_, py_), (bx, by), k2 in zip(phase, betas, K2):
        ks.append((bx, by, TWO_PI * px_, TWO_PI * py_, lam * k2, 0.0))
    if octupole:
        bx, by = betas[1]
        ks[1] = (bx, by, TWO_PI * phase[1][0], TWO_PI * phase[1][1], lam * K2[1], octupole)
    return ks


print('\n== (A) 1 個の薄い 6 極: ∂ν/∂J = -a²/(64π)[cot 3πν + 3 cot πν] (a = K2L β^{3/2}) ==')
bx0, K2L0 = 5.0, 0.6
a0 = K2L0 * bx0**1.5
cotf = lambda t: np.cos(t) / np.sin(t)
for nu1 in (0.2347, 0.41):
    pred = -a0**2 / (64 * np.pi) * (cotf(3 * np.pi * nu1) + 3 * cotf(np.pi * nu1))
    kk = [(bx0, 3.0, 0.0, 0.0, K2L0, 0.0)]
    form, im = adts_formula(kk, nu1, 0.3141)
    trk, _ = fit_adts(kk, nu1, 0.3141, 1e-5)
    print('  νx = %.4f: 閉じた式 xx = %.4f、一般の式 xx = %.4f (虚部 %.1e)、追跡 xx = %.4f (相対差 %.2e)' % (nu1, pred, form[0], im, trk[0], trk[0] / pred - 1))

print('\n== (B) 3 個の薄い 6 極 (位相・β が違う): 式 (Magnus の対の項 + 正規形の補正) と追跡 ==')
for lam in (1.0, 0.5, 0.25):
    kk = make(lam)
    form, im = adts_formula(kk, nux, nuy)
    trk, yyc = fit_adts(kk, nux, nuy, 1e-5 / lam**2 if lam < 1 else 1e-5)
    print('  強さ ×%.2f: 式 (xx, xy, yy) = (%.4f, %.4f, %.4f) (虚部 %.1e)' % (lam, *form, im))
    print('            追跡          = (%.4f, %.4f, %.4f)、相対差 = (%.2e, %.2e, %.2e)' % (*trk, *(trk / form - 1)))
print('  式は強さの 2 乗に比例 (×0.5 → ×0.25)。追跡との相対差は強さによらず同じ (J の格子 J0 ∝ 1/λ² の J の高次・フィットの誤差。強さの高次ではない)。')

print('\n== (B2) xx の明示の式 (1 次元): K_2200 = -18 Im P3 - 6 Im P1 - 9|h30000|² cot 3πν - 3|h21000|² cot πν、∂νx/∂Jx = 4 K_2200/π ==')
kk = make(1.0)
c30 = np.array([k[4] * k[0]**1.5 / 48 * np.exp(3j * k[2]) for k in kk])
c21 = np.array([k[4] * k[0]**1.5 / 16 * np.exp(1j * k[2]) for k in kk])
P3 = sum(c30[l] * np.conj(c30[e]) for l in range(len(kk)) for e in range(l))
P1 = sum(c21[l] * np.conj(c21[e]) for l in range(len(kk)) for e in range(l))
K2200 = -18 * P3.imag - 6 * P1.imag - 9 * abs(c30.sum())**2 / np.tan(3 * np.pi * nux) - 3 * abs(c21.sum())**2 / np.tan(np.pi * nux)
print('  Im P3 = %.6e、Im P1 = %.6e、|h30000|² = %.6e、|h21000|² = %.6e' % (P3.imag, P1.imag, abs(c30.sum())**2, abs(c21.sum())**2))
print('  明示の式 xx = %.6f、一般の式 (多項式の括弧と正規形) xx = %.6f' % (4 * K2200 / np.pi, adts_formula(kk, nux, nuy)[0][0]))

print('\n== (C) 6 極に 8 極 (1 個、K3L = 40) を足す: 8 極は 1 次 (h4 に足すだけ) ==')
kk = make(1.0, 40.0)
form, im = adts_formula(kk, nux, nuy)
trk, _ = fit_adts(kk, nux, nuy, 1e-5)
print('  式 = (%.4f, %.4f, %.4f)、追跡 = (%.4f, %.4f, %.4f)、相対差 = (%.2e, %.2e, %.2e)' % (*form, *trk, *(trk / form - 1)))
b8 = betas[1]
oct_only = np.array([40 * b8[0]**2 / (16 * np.pi), -40 * b8[0] * b8[1] / (8 * np.pi), 40 * b8[1]**2 / (16 * np.pi)])
kk0 = make(1.0)
f0, _ = adts_formula(kk0, nux, nuy)
print('  8 極だけの式 (K3L β²/(16π)、-K3L βxβy/(8π)、K3L βy²/(16π)) = (%.4f, %.4f, %.4f)、(C) の式 - (B) の式 = (%.4f, %.4f, %.4f)' % (*oct_only, *(form - f0)))

print('\n== (D) 対の項の順序: 6 極の並びを入れ替えると (位相を保って順序だけ入れ替え) 式が変わらないこと (位相・β を保つ) ==')
kk = make(1.0)
perm = [kk[2], kk[0], kk[1]]
# 起点を移す: 位相の基準を変える (起点を 1 番目の後ろに) → 各 kick の位相から μ1 を引き、1 周の位相を足す
mu1x, mu1y = kk[0][2], kk[0][3]
shifted = []
for (bx, by, mux, muy, k2, k3) in kk[1:]:
    shifted.append((bx, by, mux - mu1x, muy - mu1y, k2, k3))
bx, by, mux, muy, k2, k3 = kk[0]
shifted.append((bx, by, TWO_PI * nux, TWO_PI * nuy, k2, k3))     # 起点を最初の 6 極の直後に移すと、それは 1 周の最後 (位相 2πν)
fa, _ = adts_formula(kk, nux, nuy)
fb, _ = adts_formula(shifted, nux, nuy)
print('  起点を最初の 6 極の位置に移した式 = (%.6f, %.6f, %.6f)、元の起点 = (%.6f, %.6f, %.6f)' % (*fb, *fa))


print('\n== (E) 1 個の薄い 6 極と 3νx 共鳴の近く: 固定点 (周期 3 の軌道) の振幅 ==')
from scipy.optimize import fsolve
bx0, K2L0 = 5.0, 0.5
a0 = K2L0 * bx0**1.5
h30, h21 = a0 / 48, a0 / 16
kk = [(bx0, 3.0, 0.0, 0.0, K2L0, 0.0)]
for dlt in (0.004, 0.008):
    nu1 = 1 / 3 + dlt
    c = -12 * h21**2 / np.tan(np.pi * nu1)          # K ∋ c J² (非共鳴の 21 調和の ADTS)
    # 2π δ + 2 c J - 6√2 |h| √J = 0 (cos ψ = -1) を √J について解く
    r = np.roots([2 * c, -6 * np.sqrt(2) * h30, 2 * np.pi * dlt])
    r = r[np.isreal(r)].real
    r = r[r > 0]
    Jth = r**2
    Jsimple = np.pi**2 * dlt**2 / (18 * h30**2)
    found = []
    for phi0 in (np.pi / 3, np.pi, 5 * np.pi / 3):
        for Jg in Jth:
            z0 = np.array([np.sqrt(2 * Jg) * np.cos(phi0), -np.sqrt(2 * Jg) * np.sin(phi0)])
            def f(z):
                st = (np.array([z[0]]), np.array([z[1]]), np.array([0.0]), np.array([0.0]))
                for _ in range(3):
                    st = one_turn(st, kk, nu1, 0.3141)
                return [st[0][0] - z[0], st[1][0] - z[1]]
            sol, info, ier, msg = fsolve(f, z0, full_output=True, xtol=1e-13)
            if ier == 1:
                found.append(0.5 * (sol[0]**2 + sol[1]**2))
    fs = sorted(set(round(j, 12) for j in found))
    print('  δ = ν - 1/3 = %.3f: 理論の J_fp = %.4e (J² の項 c = %.3e を含む)、c を無視した π²δ²/(18|h|²) = %.4e' % (dlt, Jth[0], c, Jsimple))
    print('     追跡 (M³ の不動点) の J = %s (周期 3 の軌道の 3 点で J は少し違う)、その平均 %.4e → 理論との相対差 %.2e' % (['%.4e' % j for j in fs], np.mean(fs), np.mean(fs) / Jth[0] - 1))

print('\n== (F) NAFF の核: Hann 窓 = 3 つの矩形窓の核の和、裾の減り方、2 本の線があるときのずれ ==')
N = 256
n = np.arange(N)


def Dk(d):
    """矩形窓の核 (1/N) Σ e^{2πi d n}"""
    return np.sum(np.exp(2j * np.pi * d * n)) / N


def Sh(d):
    """Hann 窓 w_n = 1 - cos(2πn/N) の核 (1/N) Σ w_n e^{2πi d n}"""
    return np.sum((1 - np.cos(2 * np.pi * n / N)) * np.exp(2j * np.pi * d * n)) / N


dd = 3.37 / N
print('  Hann 窓の核 - [D(Δ) - D(Δ - 1/N)/2 - D(Δ + 1/N)/2] = %.2e (Δ = %.3f/N)' % (abs(Sh(dd) - (Dk(dd) - 0.5 * Dk(dd - 1 / N) - 0.5 * Dk(dd + 1 / N))), 3.37))
for m in (2.5, 5.5, 20.5, 80.5):
    print('  Δ = %.1f/N: 矩形の核 |D| = %.3e、Hann の核 |S|/|S(0)| = %.3e (S(0) = %.3f)' % (m, abs(Dk(m / N)), abs(Sh(m / N)) / abs(Sh(0.0)), abs(Sh(0.0))))
nu0, nu1, eps = 0.3036, 0.3281, 1e-3
res = {}
for N in (256, 1024, 4096):
    nn = np.arange(N)
    sig = np.exp(2j * np.pi * nu0 * nn) + eps * np.exp(2j * np.pi * nu1 * nn)
    w = 1 + np.cos(np.pi * (nn - N / 2) / (N / 2))
    f_h = lambda nu: -abs(np.sum(w * sig * np.exp(-2j * np.pi * nu * nn)))
    f_r = lambda nu: -abs(np.sum(sig * np.exp(-2j * np.pi * nu * nn)))
    kpk = int(np.argmax(abs(np.fft.fft(sig * w))))
    c = kpk / N
    eh = minimize_scalar(f_h, bracket=(c - 1.0 / N, c, c + 1.0 / N), tol=1e-14).x - nu0
    er = minimize_scalar(f_r, bracket=(c - 1.0 / N, c, c + 1.0 / N), tol=1e-14).x - nu0
    print('  N = %4d: 2 本目 (振幅 1e-3、距離 %.4f) があるときの ν0 の誤差: 矩形窓 %.2e、Hann 窓 %.2e' % (N, nu1 - nu0, er, eh))
