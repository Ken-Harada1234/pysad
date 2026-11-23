# -*- coding: utf-8 -*-
"""偏向磁石の硬い端の Maxwell フリンジ: 因子 p の独立な再確認 (06_edge_fringe.md §7、§12 の B-1 の疑い)。

(a) 追跡 (tbendcore の f・fa・ff)、ヘルプの式、(5.36) の p² 型は同じ式 (ランダムな入力で 1e-16 台)。
(b) 端の磁場 (直線座標、Maxwell を y⁴ まで満たす) をローレンツ力で積分 (DOP853)。硬い端との差の d → 0 の外挿を、
    p² 型 (追跡・ヘルプ)、p¹ 型 (CALC の COD の更新)と比べる。δ = 0, 0.3, -0.2、q = px/p = 0.1, 0.3 (p_y = 0)。
(c) CALC の解析行列 (tbedgemaxwell の trans) が、COD の更新 (p¹・p⁰) の微分でなく、p² 型の更新の微分であること。シンプレクティック性。
"""
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import lattice as lat

rho = 50.0
y0 = 0.0003


def prof(s, d):
    """g = (1 + tanh(s/d))/2 とその 1〜3 階微分"""
    t = math.tanh(s / d)
    return 0.5 * (1 + t), 0.5 * (1 - t * t) / d, -t * (1 - t * t) / d**2, -(1 - t * t) * (1 - 3 * t * t) / d**3


def g4(s, d):
    """g の 4 階微分 (g = (1 + tanh(s/d))/2)。sympy で作った式を数値化"""
    return _g4(s / d) / d**4


_u = sp.symbols('u')
_g4 = sp.lambdify(_u, sp.diff((1 + sp.tanh(_u)) / 2, _u, 4), 'math')


def run(d, c0, p, smooth, s0=-0.05, s1=0.05):
    """直線座標の運動量 P = (Px, Py, Pz)、|P| = p のローレンツ力。b = (0, g - y² g''/2 + y⁴ g''''/24, y g' - y³ g'''/6)/ρ。
    戻り値: (x, Px, y, Py, 経路長 ∫ p/Pz ds)"""
    def rhs(s, v):
        xx, Px, yy, Py, L = v
        if smooth:
            gg, g1, g2, g3 = prof(s, d)
            g4v = g4(s, d)
        else:
            gg, g1, g2, g3, g4v = (1.0 if s >= 0 else 0.0), 0.0, 0.0, 0.0, 0.0
        byy = (gg - yy**2 * g2 / 2 + yy**4 * g4v / 24) / rho
        bss = (yy * g1 - yy**3 * g3 / 6) / rho
        Pz = math.sqrt(p * p - Px * Px - Py * Py)
        return [Px / Pz, (Py * bss - Pz * byy) / Pz, Py / Pz, -Px * bss / Pz, p / Pz]
    kw = dict(rtol=1e-12, atol=1e-18, method='DOP853')
    if smooth:
        return solve_ivp(rhs, (s0, s1), c0, max_step=d / 5, **kw).y[:, -1]
    r1 = solve_ivp(rhs, (s0, 0.0), c0, **kw).y[:, -1]
    return solve_ivp(rhs, (0.0, s1), r1, **kw).y[:, -1]


# ---- 閉じた式 (入口)。座標は canonical (p0 で規格化)、p = 1 + δ
def k64_track(x, px, y, py, dl):
    """tbendcore の f・fa・ff (q = p/p の単位) → canonical の Δx, Δpy, Δz"""
    p = 1 + dl
    q = px / p
    pv1sqi = 1 / (1 - q**2)
    fa = y / (rho * p) * math.sqrt(pv1sqi)
    f = (1 - (y / rho)**2 / 6) * fa
    ff = .25 * (f + fa) * y * pv1sqi
    return ff, -q * f * p, -q * ff * p / p * p / p     # Δz = -px1*ff (q 単位。canonical の z と同じ量)


def help_form(x, px, y, py, dl):
    p = 1 + dl
    pv = math.sqrt(p * p - px * px)
    a = 1 - y**2 / (6 * rho**2)
    b = .5 - y**2 / (24 * rho**2)
    return b * y**2 / rho * p**2 / pv**3, -a * y / rho * px / pv, -b * y**2 / rho * px * p / pv**3


def calc_form(x, px, y, py, dl):
    cod = np.array([x, px, y, py, 0.0, dl])
    o = lat._tbedgemaxwell_cod(cod, rho)
    return o[0] - x, o[3] - py, o[4]


print('== (a) 追跡の式・ヘルプの式の同値 ==')
rng = np.random.default_rng(3)
mx = 0.0
for _ in range(2000):
    dl = rng.uniform(-0.3, 0.3)
    p = 1 + dl
    px = rng.uniform(-0.3, 0.3) * p
    y = rng.uniform(-0.01, 0.01)
    k = k64_track(0, px, y, 0, dl)
    h = help_form(0, px, y, 0, dl)
    mx = max(mx, abs(k[0] - h[0]), abs(k[1] - h[1]), abs(k[2] - h[2]))
print('  ランダム 2000 点 (δ ∈ ±0.3、px/p ∈ ±0.3、y ∈ ±10 mm、ρ = 50 m) での最大差 = %.2e (追跡の式 = ヘルプの式 = p² 型)' % mx)

print('\n== (b) ローレンツ力の積分 (ρ = %.0f m、y = %.1f mm、端の形 tanh、d = 1, 0.5 mm → d → 0 に外挿、窓は ±8d) ==' % (rho, y0 * 1e3))
print('  (場の y の展開は y ≪ d でだけ収束するので y = 0.3 mm。窓を ±8d に短くするのは、p_y ≠ 0 で y が動くと展開が使えなくなるため。p_y = 0 だけを調べる。)')
print('  表: 比 = (ODE の外挿値)/(式)。Δx・Δpy・Δz の順。p² 型 = 追跡・ヘルプ、p 型 = CALC の COD の更新 (Δx は p¹、Δz は p⁰)')
for q in (0.1, 0.3):
    for dl in (0.0, 0.3, -0.2):
        p = 1 + dl
        c0 = [0.0, q * p, y0, 0.0, 0.0]
        res = {}
        for d in (0.001, 0.0005):
            a = run(d, c0, p, True, s0=-8 * d, s1=8 * d)
            b = run(d, c0, p, False, s0=-8 * d, s1=8 * d)
            res[d] = np.array([a[0] - b[0], a[3] - b[3], -(a[4] - b[4])])      # Δx, Δpy, Δz (z = -経路長の差)
        ext = (4 * res[0.0005] - res[0.001]) / 3          # d² の項を消す
        H = np.array(help_form(0, q * p, y0, 0.0, dl))
        Cc = np.array(calc_form(0, q * p, y0, 0.0, dl))
        print('  q = %.1f、δ = %+.1f: ODE = (%.4e, %.4e, %.4e)' % (q, dl, *ext))
        print('        比 p² 型 = (%.4f, %.4f, %.4f)   p 型 (CALC) = (%.4f, %.4f, %.4f)' % (*(ext / H), *(ext / Cc)))
print('  (Δx は p² 型と 1.000。CALC の p 型は δ ≠ 0 で 1 + δ 倍ずれる。Δz は ODE の経路長の差で、直線座標の z の定義 (z = -経路長の差) による値。)')

print('\n== (c) CALC の解析行列と COD の更新 ==')
Sg = np.zeros((6, 6))
for i in range(3):
    Sg[2 * i, 2 * i + 1] = 1
    Sg[2 * i + 1, 2 * i] = -1


def jac(f, c, h=1e-6):
    M = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        M[:, j] = (f(cp) - f(cm)) / (2 * h)
    return M


def cod_p2(cc):
    """p² 型 (Δx = bf p²、Δz = -bf px p) の COD の更新 (追跡の式と同じ)"""
    pr = 1 + cc[5]
    pxi = cc[1]
    pv2 = pr * pr - pxi * pxi
    pv = math.sqrt(pv2)
    yi = cc[2]
    yr = (yi / rho)**2
    af = (1 - yr / 6) * yi / rho / pv
    bf = (.5 - yr / 24) * yi**2 / rho / pv2 / pv
    o = cc.copy()
    o[0] += bf * pr * pr
    o[3] -= af * pxi
    o[4] -= bf * pxi * pr
    return o


for dl in (0.0, 0.1, -0.1):
    c = np.array([0.001, 0.02, 0.01, 0.01, 0.0, dl])
    Ma = lat._tbedgemaxwell_mat(c, rho)
    Mc = jac(lambda cc: lat._tbedgemaxwell_cod(cc, rho), c)
    M2 = jac(cod_p2, c)
    sy = lambda M: np.abs(M.T @ Sg @ M - Sg).max()
    print('  δ = %+.1f: |解析行列 - COD 更新(p 型)の微分| = %.2e   |解析行列 - p² 型 (Δz = -bf px p) の微分| = %.2e   シンプレクティック性 |MᵀSM - S|: 解析行列 %.1e、COD 更新 %.1e、p² 型 %.1e' % (
        dl, np.abs(Ma - Mc).max(), np.abs(Ma - M2).max(), sy(Ma), sy(Mc), sy(M2)))
    diffs = np.abs(Ma - Mc)
    idx = np.unravel_index(np.argmax(diffs), diffs.shape)
    print('        解析行列と COD 更新の微分が最も違う成分 = (%d, %d)、行列 %.6e、更新の微分 %.6e' % (idx[0] + 1, idx[1] + 1, Ma[idx], Mc[idx]))
print('  → CALC の解析行列は p² 型 (Δx) と、Δz = -bf px p (p¹) の微分に一致し、COD の更新 (p¹・p⁰) と合わない。')

print('\n== (d) CALC の更新が因子 p を落とすことによる 1 つの端面の Δx のずれの大きさ ==')
for yy_, dl_, rho_ in ((0.010, 0.05, 10.0), (0.002, 0.05, 50.0), (0.010, 0.01, 10.0)):
    dx_true = yy_**2 / (2 * rho_)          # 桁の見積り: Δx ≈ y²/(2ρ_b) (p_x = 0 では (y²/2ρ_b)/p。p ≈ 1 とした)
    print('  y = %.0f mm、δ = %.2f、ρ_b = %.0f m: δ y²/(2ρ_b) = %.2e m  (y²/(2ρ_b) = %.2e m)' % (yy_ * 1e3, dl_, rho_, dl_ * dx_true, dx_true))
print('  (入口と出口で f が逆符号なので、Δx のずれも逆符号。1 つの端面のずれは上の値。入口と出口の合計は、対称な磁石ではほぼ相殺する (未確認))')
