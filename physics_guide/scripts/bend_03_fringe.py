# -*- coding: utf-8 -*-
"""偏向磁石の端部 (Maxwell フリンジと線形フリンジ F1) (05_bend.md §7〜§8、式 (5.41)〜(5.60))。

(a) 端の磁場の展開 (∇·b = 0、∇×b = 0)。
(b) 硬い端の Maxwell フリンジの閉じた式 (Δx、Δpy) を、端の幅 d の tanh の磁場の中のローレンツ力の積分 (直線座標、DOP853) と比べる。
    δ = 0 と δ = 0.3 で、Δx の p^2 の因子 (追跡・ヘルプ) と p の因子 (CALC の tbedgemaxwell) のどちらかを見分ける。
(c) CALC の tbedgemaxwell と追跡の式のシンプレクティック性 (6 変数、中心差分) と、δ ≠ 0 の差。
(d) 線形フリンジ: F1 = 6 ∫ (g - g^2) ds (ヘルプの定義) と y・y^3 の項、x の δ 依存の項 (ローレンツ力の積分)。
(e) tbendal の rbh・rbl、DISFRIN が効く場所。
"""
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad
import _common as C

C.start(__file__)

print('== (a) 端の磁場 (b = eB/P0、直線座標 (x, y, s)、一様な B0 = 1/ρ の端) ==')
y_, s_ = sp.symbols('y s', real=True)
g = sp.Function('g')
V = sum((-1)**k * y_**(2 * k + 1) * sp.diff(g(s_), s_, 2 * k) / sp.factorial(2 * k + 1) for k in range(4))
by = sp.diff(V, y_)
bs = sp.diff(V, s_)
print('  スカラーポテンシャル V = y g - y^3 g\'\'/6 + y^5 g\'\'\'\'/120 - ... (∇²V = 0 を満たす。b_x = ∂_x V = 0 の直線の端)')
print('  b_y = ∂_y V = %s' % sp.series(sp.expand(by), y_, 0, 5).removeO())
print('  b_s = ∂_s V = %s' % sp.series(sp.expand(bs), y_, 0, 4).removeO())
div = sp.expand(sp.diff(by, y_) + sp.diff(bs, s_))
print('  ∇·b = ∂_y b_y + ∂_s b_s = %s  (y の 5 次 (V の打ち切り) だけが残る)' % sp.simplify(div))
print('  ∇×b: (∇×b)_x = ∂_y b_s - ∂_s b_y = %s (同じく V の打ち切りだけ)、b_x = 0 から他の成分は 0' % sp.simplify(sp.expand(sp.diff(bs, y_) - sp.diff(by, s_))))

# ------------------------------------------------------------------
rho = 50.0
B0 = 1 / rho


def prof(s, d):
    t = math.tanh(s / d)
    return 0.5 * (1 + t), 0.5 * (1 - t * t) / d, -t * (1 - t * t) / d**2, -(1 - t * t) * (1 - 3 * t * t) / d**3


def run(d, c0, p, smooth=True, s0=None, s1=0.05):
    """直線座標の機械運動量 P = (Px, Py, Pz)、|P| = p のローレンツ力。b_x = 0、b_y = B0 (g - y^2 g''/2)、b_s = B0 (y g' - y^3 g'''/6)。
    smooth = False は s = 0 に硬い端 (b_y = B0 Θ(s)、b_s = 0)。c0 = (x, Px, y, Py)"""
    if s0 is None:
        s0 = -0.05

    def rhs(s, v):
        xx, Px, yy, Py = v
        if smooth:
            gg, g1, g2, g3 = prof(s, d)
        else:
            gg, g1, g2, g3 = (1.0 if s >= 0 else 0.0), 0.0, 0.0, 0.0
        byy = B0 * (gg - yy**2 * g2 / 2)
        bss = B0 * (yy * g1 - yy**3 * g3 / 6)
        Pz = math.sqrt(p * p - Px * Px - Py * Py)
        return [Px / Pz, (Py * bss - Pz * byy) / Pz, Py / Pz, -Px * bss / Pz]
    if smooth:
        return solve_ivp(rhs, (s0, s1), c0, rtol=1e-12, atol=1e-16, method='DOP853', max_step=d / 5).y[:, -1]
    r1 = solve_ivp(rhs, (s0, 0.0), c0, rtol=1e-12, atol=1e-16, method='DOP853').y[:, -1]
    return solve_ivp(rhs, (0.0, s1), r1, rtol=1e-12, atol=1e-16, method='DOP853').y[:, -1]


print('\n== (b) 硬い端の Maxwell フリンジ: 閉じた式とローレンツ力の積分 ==')
print('  (直線座標、1/ρ = %.3f /m (ρ = %.0f m)、入射 x = 0、y = 2 mm、q = px/p = 0.1、py = 0。端の形 g = (1 + tanh(s/d))/2、d = 2・1 mm)' % (B0, rho))
print('  Δ = (磁場が滑らかな端) - (硬い端 s = 0) の、s = +5 cm での x と py の差。Δx は d^2 に比例する項 (線形フリンジの F1^2/(24ρ) 型) と、d に依らない項 A の和なので、d = 2・1 mm の 2 点から A を外挿する。')
y0 = 0.002
for dl in (0.0, 0.3):
    p = 1 + dl
    q = 0.1
    c0 = [0.0, q * p, y0, 0.0]
    pv = math.sqrt(p * p - (q * p)**2)
    res = {}
    for d in (0.002, 0.001):
        a = run(d, c0, p, True)
        b = run(d, c0, p, False)
        res[d] = (a[0] - b[0], a[3] - b[3])
    A = (4 * res[0.001][0] - res[0.002][0]) / 3          # Δx = A + c d^2 → A = (4 Δx(d/2) - Δx(d))/3
    dpy = res[0.001][1]
    fx2 = y0**2 / (2 * rho) * p**2 / pv**3
    fx1 = y0**2 / (2 * rho) * p / pv**3
    fpy = -(q * p) / pv * y0 / rho
    print('  δ = %.1f: Δx の d → 0 の外挿 A = %.4e   閉じた式 y^2/(2ρ) p^2/pv^3 = %.4e (比 %.4f)、p の因子だけにした式 = %.4e (比 %.4f)' % (dl, A, fx2, A / fx2, fx1, A / fx1))
    print('         Δpy (d = 1 mm) = %.4e   閉じた式 -(px/pv) y/ρ = %.4e (比 %.5f)   (pv = sqrt(p^2 - px^2))' % (dpy, fpy, dpy / fpy))

print('\n== (c) CALC の tbedgemaxwell と追跡の式 (tbendcore の f・ff) のシンプレクティック性 ==')
from pysad import lattice as lat
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


rhob = 10.0


def trk(cc, ent=True):
    """追跡の式 (tbend.f tbendcore の f, fa, ff) を正準変数で。q = px/p"""
    x, px, yy, py, z, dp = cc
    pr = 1 + dp
    q, qy = px / pr, py / pr
    pv1sqi = 1 / max(1e-4, 1 - q**2)
    fa = yy / (rhob * pr) * math.sqrt(pv1sqi)
    f = (1 - (yy / rhob)**2 / 6) * fa
    ff = .25 * (f + fa) * yy * pv1sqi
    sg = 1 if ent else -1
    return np.array([x + sg * ff, px, yy, (qy - sg * q * f) * pr, z - sg * q * ff, dp])


print('  入口 (ρb = %.0f m): 6 変数のヤコビアンの |J^T S J - S| (中心差分 1e-6。10^-10 程度までが数値微分の精度)' % rhob)
print('  %-34s %-18s %-18s' % ('(x, px, y, py, z, δ)', 'CALC tbedgemaxwell', '追跡の式'))
for c in ([0.001, 0.02, 0.02, 0.01, 0.0, 0.0], [0.001, 0.02, 0.02, 0.01, 0.0, 0.1], [0.001, 0.02, 0.02, 0.01, 0.0, -0.1]):
    c = np.array(c)
    Mc = jac(lambda cc: lat._tbedgemaxwell_cod(cc, rhob), c)
    Mt = jac(lambda cc: trk(cc, True), c)
    print('  %-34s %-18.2e %-18.2e' % (str([float(v) for v in c]), np.abs(Mc.T @ Sg @ Mc - Sg).max(), np.abs(Mt.T @ Sg @ Mt - Sg).max()))
c = np.array([0.001, 0.02, 0.02, 0.01, 0.0, 0.1])
dc = lat._tbedgemaxwell_cod(c, rhob) - c
dt = trk(c, True) - c
print('  δ = 0.1 の 1 回の更新: CALC (Δx, Δz) = (%.4e, %.4e)、追跡 (Δx, Δz) = (%.4e, %.4e)。比 = %.4f、1/(1 + δ) = %.4f。Δpy は CALC %.4e、追跡 %.4e (同じ)' % (
    dc[0], dc[4], dt[0], dt[4], dc[0] / dt[0], 1 / 1.1, dc[3], dt[3]))

print('\n== (d) 線形フリンジ F1 ==')
print('  ヘルプの定義 F1 = 6 ∫ (g - g^2) ds (g = B_y(s)/B0)。端の形ごとの F1、∫g\'^2 ds、1/F1:')
d_ = 0.002
Ftanh = 6 * quad(lambda s: prof(s, d_)[0] * (1 - prof(s, d_)[0]), -0.1, 0.1, points=[0], epsabs=1e-14)[0]
G2tanh = quad(lambda s: prof(s, d_)[1]**2, -0.1, 0.1, points=[0], epsabs=1e-14)[0]
ell = 0.006
Framp = 6 * ell / 6
G2ramp = 1 / ell
sig = 0.002
gerf = lambda s: 0.5 * (1 + math.erf(s / (math.sqrt(2) * sig)))
Ferf = 6 * quad(lambda s: gerf(s) * (1 - gerf(s)), -0.05, 0.05, points=[0], epsabs=1e-14)[0]
G2erf = quad(lambda s: (math.exp(-s * s / (2 * sig * sig)) / (math.sqrt(2 * math.pi) * sig))**2, -0.05, 0.05, epsabs=1e-14)[0]
print('  %-28s F1 = %.5f m   ∫g\'^2 ds = %.3f /m   1/F1 = %.3f /m   比 (∫g\'^2 F1) = %.4f' % ('直線の立ち上がり (長さ 6 mm)', Framp, G2ramp, 1 / Framp, G2ramp * Framp))
print('  %-28s F1 = %.5f m (= 3d)   ∫g\'^2 ds = %.3f /m   1/F1 = %.3f /m   比 = %.4f' % ('tanh (d = 2 mm)', Ftanh, G2tanh, 1 / Ftanh, G2tanh * Ftanh))
print('  %-28s F1 = %.5f m   ∫g\'^2 ds = %.3f /m   1/F1 = %.3f /m   比 = %.4f' % ('誤差関数 (σ = 2 mm)', Ferf, G2erf, 1 / Ferf, G2erf * Ferf))
print('  → ヘルプの Δy_fra = (2/3)/(F1 ρ^2) は、(2/3)(1/ρ^2) ∫g\'^2 ds を F1 で表した形。直線・tanh では厳密に一致、誤差関数では %.1f%% 違う。' % (100 * (G2erf * Ferf - 1)))
print('  ローレンツ力の積分: 入射 px = 0、x = 0、s = ±5 cm の間。tanh (d = 2 mm、F1 = 3d = 6 mm)、ρ = %.0f m。Δpy を y と y^3 に分ける (y = 0.5・1 mm の 2 点)。' % rho)
for dl in (0.0, 0.3):
    p = 1 + dl
    vals = {}
    for yv in (0.0005, 0.001):
        a = run(d_, [0.0, 0.0, yv, 0.0], p, True)
        vals[yv] = a[3]
    # Δpy = a1 y - a3 y^3
    A = np.array([[0.0005, -0.0005**3], [0.001, -0.001**3]])
    a1, a3 = np.linalg.solve(A, np.array([vals[0.0005], vals[0.001]]))
    f1 = Ftanh
    print('  δ = %.1f: a1 = %.5e  (F1/(6ρ^2) = %.5e、比 %.4f。p を掛けた値 a1 p = %.5e)、a3 = %.5e  ((2/3)/(F1ρ^2) = %.5e、比 %.4f)' % (
        dl, a1, f1 / (6 * rho**2), a1 / (f1 / (6 * rho**2)), a1 * p, a3, (2 / 3) / (f1 * rho**2), a3 / ((2 / 3) / (f1 * rho**2))))
print('  (Δpy は機械運動量の p0 規格化。式: canonical Δpy = y (F1/(6ρ^2) - y^2 (2/3)/(F1ρ^2))/p。比 1 に近く、δ = 0.3 では a1 p が δ = 0 の a1 と同じ)')
print('  x のずれ (入射 px = 0、y = 0)、硬い端 (s = 0) との差 ΔX(p)。SAD の式: 入口で x に dxfr δ/p (dxfr = F1^2/(24 ρb)) を足す。')
dxfr = Ftanh**2 / (24 * rho)


def run_ramp(ell, c0, p, smooth=True, s1=0.05):
    """直線の立ち上がり g = (s + ell/2)/ell (|s| < ell/2) の端。y = 0 の水平運動だけ (b_y = B0 g)"""
    def mk(gf):
        def rhs(s, v):
            xx, Px = v
            Pz = math.sqrt(p * p - Px * Px)
            return [Px / Pz, -B0 * gf(s)]
        return rhs
    v = list(c0)
    s0 = -0.05
    if smooth:
        for (a, b, gf) in ((s0, -ell / 2, lambda s: 0.0), (-ell / 2, ell / 2, lambda s: (s + ell / 2) / ell), (ell / 2, s1, lambda s: 1.0)):
            v = solve_ivp(mk(gf), (a, b), v, rtol=1e-13, atol=1e-17, method='DOP853').y[:, -1]
    else:
        for (a, b, gf) in ((s0, 0.0, lambda s: 0.0), (0.0, s1, lambda s: 1.0)):
            v = solve_ivp(mk(gf), (a, b), v, rtol=1e-13, atol=1e-17, method='DOP853').y[:, -1]
    return v


Sx, Sr = {}, {}
for dl in (0.0, 0.3, -0.2):
    p = 1 + dl
    Sx[dl] = run(d_, [0.0, 0.0, 0.0, 0.0], p, True)[0] - run(d_, [0.0, 0.0, 0.0, 0.0], p, False)[0]
    Sr[dl] = run_ramp(ell, [0.0, 0.0], p, True)[0] - run_ramp(ell, [0.0, 0.0], p, False)[0]
dxr = Framp**2 / (24 * rho)
print('  %-8s %-26s %-26s' % ('δ', 'tanh (d = 2 mm、F1 = 3d)', '直線の立ち上がり (6 mm、F1 = 6 mm)'))
for dl in (0.0, 0.3, -0.2):
    p = 1 + dl
    print('  %+.1f     ΔX p/dxfr = %+.4f           ΔX p/dxfr = %+.4f' % (dl, Sx[dl] * p / dxfr, Sr[dl] * p / dxr))
print('  ΔX(δ) - ΔX(0) と dxfr δ/p の比 (δ = 0.3・-0.2):  tanh %.4f・%.4f、直線 %.4f・%.4f   (π^2/9 = %.4f)' % (
    (Sx[0.3] - Sx[0.0]) / (dxfr * 0.3 / 1.3), (Sx[-0.2] - Sx[0.0]) / (dxfr * (-0.2) / 0.8), (Sr[0.3] - Sr[0.0]) / (dxr * 0.3 / 1.3), (Sr[-0.2] - Sr[0.0]) / (dxr * (-0.2) / 0.8), math.pi**2 / 9))
print('  → ΔX = -F1^2/(24 ρ p)。直線の立ち上がりでは厳密 (比 1)。tanh では F1 を 6∫g(1-g) で定義すると π^2/9 = 1.097 倍 (x のずれは、1 次のモーメント I1 (04 章の定義 π d) で決まり、y の収束は ∫g(1-g) で決まる別の量)。')
print("    δ = 0 の分 -dxfr は L' の補正と ρb の定義で設計軌道に吸収され、K64 の式は δ 依存分 dxfr δ/p = dxfr (1 - 1/p) だけを入れる。")

print('\n== (e) tbendal の rbh・rbl、DISFRIN の効き方 ==')
rbh = .5 + 1 / math.sqrt(12)
rbl = 1 / 6 / rbh
print('  rbh = 1/2 + 1/√12 = %.10f、rbl = 1/(6 rbh) = %.10f、rbh + rbl = %.12f、rbh rbl = %.12f (= 1/6)。区間 [0,1] の 2 点 Gauss-Legendre の節に一致 (1/2 ± 1/(2√3) = %.10f, %.10f)' % (
    rbh, rbl, rbh + rbl, rbh * rbl, .5 - 1 / (2 * math.sqrt(3)), .5 + 1 / (2 * math.sqrt(3))))
import types
for extra, name in (('E1=0.5 E2=0.5', 'F1 = 0 (FRINGE = 0)'), ('E1=0.5 E2=0.5 F1=0.2 FRINGE=1', 'F1 = 0.2、FRINGE = 1')):
    ffs0, _ = C.make_bend(1.0, 0.1, extra)
    ffs1, _ = C.make_bend(1.0, 0.1, extra + ' DISFRIN=1')
    c = np.array([0.001, 0.0005, 0.01, 0.0002, 0.0, 0.01])
    print('  %-24s DISFRIN = 0 と 1 の追跡の差 (y = 10 mm): %s' % (name, np.array2string(C.track_one(ffs0, c) - C.track_one(ffs1, c), precision=2)))
print('  → DISFRIN は F1 の y^3 の項 (dyfra) だけを止める。硬い端の Maxwell フリンジ (f, ff) は DISFRIN でも止まらない (K1 = 0 の BEND。tbendcore に判定が無い)。')
