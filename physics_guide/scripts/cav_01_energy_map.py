# -*- coding: utf-8 -*-
"""空洞 CAVI・TCAVI の写像 (10_cavity.md): 追跡 (K64 の tcav・ttcav の pysad 版) を、導出した式で組み直した写像と比べる。電子。
 (a) 薄い CAVI (L = 0) の δ・z・px の更新: ΔE = -eV sin(ωt - dphis) の形、時間 t = -z h1/p1、h2 = h1 + Δh、p2 = √(h2² - 1)、z2 = -t p2/h2、q2 = q1 p1/p2。
 (b) シンプレクティック性 (6×6 のヤコビアン、中心差分): 薄い CAVI (V1・V20 つき)、L ≠ 0 でフリンジ (H_f) つき。
 (c) 非相対論の補正 (VOLT+DVOLT)(2π FREQ/c)²/(γβ)²/4: 低エネルギー (2 MeV/c の電子) で、ΔE の x 依存が I0 の展開 1 + (ω x/(βγ c))²/4 に合う。
 (d) 長い空洞の分割の重み ws の和が 1、e^{wl} = 1 + r。
 (e) TCAVI: 写像の式、シンプレクティック性、ハミルトニアン K x sin(φ - ω z/β) の流れとの比較 (z の更新の違い)。
"""
import io
import math
import numpy as np
from scipy.special import iv
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6, track6n_ext

MC2 = 0.51099895e6
c_light = 299792458.0


def mk(el, momentum='3GEV'):
    src = 'MOMENTUM = %s;\n%s\nMARK ST = ();\nLINE ONE = (ST CV);\nFFS;\nUSE ONE;\n' % (momentum, el)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def trk(ffs, c, trf0=0.0, ext=False):
    a = np.zeros((1, 6))
    a[0, :] = c
    fn = track6n_ext.track_particles_ext if ext else track6.track_particles
    return np.asarray(fn(ffs.lat.elements, a, 1, ffs.momentum(), trf0, True, False, False))[0]


def jac(f, c, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h * max(1.0, abs(c[j]))
        J[:, j] = (f(c + d) - f(c - d)) / (2 * d[j])
    return J


S = np.zeros((6, 6))
for i in (0, 2, 4):
    S[i, i + 1] = 1
    S[i + 1, i] = -1
# 座標順 (x, px, y, py, z, δ) の正準対は (0,1)(2,3)(4,5)

print('== (a) 薄い CAVI の写像: 3 GeV の電子、VOLT = 2 MV、FREQ = 500 MHz、PHI = 0.3 rad、DISFRIN = 1 ==')
p_eV = 3.0e9
p0 = p_eV / MC2
h0 = math.sqrt(1 + p0**2)
V, f_rf, phi = 2.0e6, 500.0e6, 0.3
w = 2 * math.pi * f_rf / c_light
v = V / MC2
ffs = mk('CAVI CV=(L=0 VOLT=%g FREQ=%g PHI=%g DISFRIN=1);' % (V, f_rf, phi))


def my_cav(c, trf0, rad=False):
    x, px, y, py, z, dl = c
    p1r = 1 + dl
    p1 = p0 * p1r
    h1 = math.sqrt(1 + p1**2)
    t = -z * h1 / p1                       # 位置 z から時間 (c t [m])。z = -v (t - t0)、v = p1/h1 (c = 1)
    phic = phi
    phis = w * trf0
    dphis = phis - phic
    offset = math.sin(dphis) if rad else math.sin(dphis) - math.sin(phis)
    dh = -v * (math.sin(w * t - dphis) + math.sin(dphis) - offset)
    h2 = h1 + dh
    p2 = math.sqrt(h2**2 - 1)
    return np.array([x, px, y, py, -t * p2 / h2, p2 / p0 - 1])        # 正準運動量 px, py は変わらない (追跡の内部変数 q = p/(1+δ) は p1/p2 倍になる)


worst = 0.0
for trf0 in (0.0, 0.4 / w):
    for c in (np.array([1e-3, 2e-4, -1e-3, 1e-4, -0.1, 0.002]), np.array([0, 0, 0, 0, 0.05, -0.004]), np.array([2e-3, -1e-4, 5e-4, 3e-4, 0.3, 0.0])):
        a = trk(ffs, c, trf0)
        b = my_cav(c, trf0)
        d = np.abs(a - b).max()
        worst = max(worst, d)
        print('   trf0 = %.3e: z = %+.3f m、δ: %+.5f → pysad %+.9f、式 %+.9f、最大差 %.1e' % (trf0, c[4], c[5], a[5], b[5], d))
print('   最大差 %.1e (x, px, y, py, z, δ)。ΔE = -eV[sin(ωt - dphis) + sin(phis)] (RAD OFF、offset = sin(dphis) - sin(phis))' % worst)

print('\n== (b) シンプレクティック性 (中心差分のヤコビアン) ==')
c = np.array([1e-3, 2e-4, -1e-3, 1e-4, -0.1, 0.002])
for label, el, ext in (
        ('薄い CAVI (V1 = V20 = 0)', 'CAVI CV=(L=0 VOLT=%g FREQ=%g PHI=%g DISFRIN=1);' % (V, f_rf, phi), False),
        ('薄い CAVI (V1, V20, V11, V02 あり)', 'CAVI CV=(L=0 VOLT=%g FREQ=%g PHI=%g DISFRIN=1 V1=2E7 V20=-3E8 V11=1E8 V02=2E8);' % (V, f_rf, phi), True),
        ('L = 0.5、フリンジ (H_f) あり', 'CAVI CV=(L=0.5 VOLT=%g FREQ=%g PHI=%g);' % (V, f_rf, phi), False),
        ('L = 0.5、フリンジなし (DISFRIN = 1)', 'CAVI CV=(L=0.5 VOLT=%g FREQ=%g PHI=%g DISFRIN=1);' % (V, f_rf, phi), False)):
    ffs_b = mk(el)
    J = jac(lambda cc: trk(ffs_b, cc, 0.0, ext), c)
    print('   %-34s |J^T S J - S| = %.2e、det J = %.12f' % (label, np.abs(J.T @ S @ J - S).max(), np.linalg.det(J)))

print('\n== (c) 非相対論の補正: 2 MeV/c の電子、VOLT = 1 MV、FREQ = 500 MHz、z = 0 (t = 0)、PHI = π/2 の ΔE の x 依存 ==')
pe_ = 2.0e6 / MC2
ffs_c = mk('CAVI CV=(L=0 VOLT=1E6 FREQ=500E6 PHI=%.15g DISFRIN=1);' % (math.pi / 2), momentum='2MEV')
h_c = math.sqrt(1 + pe_**2)


def dE(x):
    out = trk(ffs_c, np.array([x, 0, 0, 0, 0.0, 0.0]), 0.0)
    p2 = pe_ * (1 + out[5])
    return math.sqrt(1 + p2**2) - h_c                       # mc² 単位


d0 = dE(0.0)
for x in (0.01, 0.03, 0.05):
    ratio = dE(x) / d0
    k = w * x / pe_
    print('   x = %.2f m: ΔE(x)/ΔE(0) = %.6f、1 + (ωx/(βγc))²/4 = %.6f、I0(ωx/(βγc)) = %.6f' % (x, ratio, 1 + k * k / 4, iv(0, k)))
print('   (βγ = %.4f。ΔE の x 依存は、ヘルプの「V20・V02 に (VOLT+DVOLT)(2π FREQ/c)²/(γβ)²/4 が自動で足される」に合う。)' % pe_)

print('\n== (d) 長い空洞の分割の重み ws (vnominal ≠ 0 のとき。tcav.f) ==')
for r, ndiv in ((1e-3, 10), (0.05, 8), (0.5, 12), (3.0, 20)):
    wl = math.asinh(r * (2 + r) / 2 / (1 + r))
    r1 = wl / ndiv / 2
    ws1 = 2 * math.exp(r1) * math.sinh(r1)
    we = 1 + ws1
    ws = [ws1 / r]
    for i in range(1, ndiv):
        ws.append(ws[-1] * we)
    print('   r = %-6g ndiv = %2d: Σ ws = %.15f、exp(wl) = %.12f、1 + r = %.12f、ws(n)/ws(1) = %.4f' % (r, ndiv, sum(ws), math.exp(wl), 1 + r, ws[-1] / ws[0]))

print('\n== (e) TCAVI: K0 = 1e-4 rad (kick の角)、FREQ = 500 MHz、PHI = 0.4、3 GeV ==')
ffs_t = mk('TCAVI CV=(L=0 K0=1E-4 FREQ=500E6 PHI=0.4);')
K = 1.0e-4


def my_tcav(c):
    x, px, y, py, z, dl = c
    p1r = 1 + dl
    p1 = p0 * p1r
    h1 = math.sqrt(1 + p1**2)
    t = -z * h1 / p1
    ph = w * t + 0.4
    dh = K * p0 * w * x * math.cos(ph)
    px2 = px - K * math.sin(ph)           # 正準の px (p0 規格化)。内部の q では px_q - K sin/p1r、その後 p1r/p2r 倍
    h2 = h1 + dh
    p2 = math.sqrt(h2**2 - 1)
    p2r = p2 / p0
    return np.array([x, px2, y, py, -t * p2 / h2, p2r - 1])


worst = 0.0
cs = (np.array([1e-3, 2e-4, -1e-3, 1e-4, -0.1, 0.002]), np.array([3e-3, 0, 2e-3, 0, 0.2, -0.003]))
for c in cs:
    a = trk(ffs_t, c)
    b = my_tcav(c)
    worst = max(worst, np.abs(a - b).max())
print('   pysad の追跡と式 (px -= K sin(ωt + φ)/p1r、Δh = K p0 ω x cos(ωt + φ)、時間を保って z を更新) の最大差 %.1e' % worst)
J = jac(lambda cc: trk(ffs_t, cc), cs[0])
print('   シンプレクティック性 |J^T S J - S| = %.2e、det J = %.12f' % (np.abs(J.T @ S @ J - S).max(), np.linalg.det(J)))


def flow(c):
    """H = K p0 x sin(φ - ω z/β(δ)) の時間 1 の流れ (px は p0 規格化の正準運動量、δ は z と共役。dβ/dδ = p0/h³)"""
    def rhs(s, u):
        x, px, y, py, z, dl = u
        p1 = p0 * (1 + dl)
        h1 = math.sqrt(1 + p1**2)
        b = p1 / h1
        db = p0 / h1**3
        arg = 0.4 - w * z / b
        return [0.0, -K * math.sin(arg), 0.0, 0.0, K * x * math.cos(arg) * (w * z / b**2) * db, K * x * math.cos(arg) * w / b]
    return solve_ivp(rhs, (0, 1), list(c), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1]


a = trk(ffs_t, cs[0])
b = flow(cs[0])
print('   追跡 − H の流れ (x, px, y, py, z, δ) の最大 %.1e (z の変化の大きさ %.1e、δ の変化 %.1e)' % (np.abs(a - b).max(), abs(a[4] - cs[0][4]), abs(a[5] - cs[0][5])))
print('   (3 GeV では z の変化が 1e-15 以下で、H の流れと「時間を保つ z の更新」の違いは見えない。δ・px の一致だけが確認。z の更新の違いは未確認。)')

print('\n== (f) 入口の端 (tcavfrin) の収束: Δpx = x · v/(2 L p0) · (sin(ωt - dphis) + sin(phis)) / 1 (正準 px、p0 規格化) ==')
print('   Maxwell から (10 章 §7): Δp_r = -(e r/2)(E_0/(β c)) cos(ωt + φ_E) で、β = 1 の近似。K64 の式は 1/β を持たない。FRINGE = 1 (入口だけ)、L = 0.1 m、DISFRIN = 0 と 1 の差。')
for mom, name in ((3.0e9, '3GEV'), (2.0e6, '2MEV'), (0.6e6, '600KEV')):
    pp = mom / MC2
    hh = math.sqrt(1 + pp * pp)
    beta = pp / hh
    Lc = 0.1
    ph_ = 0.3
    e1 = mk('CAVI CV=(L=%g VOLT=1E5 FREQ=500E6 PHI=%g FRINGE=1);' % (Lc, ph_), momentum=name)
    e0 = mk('CAVI CV=(L=%g VOLT=1E5 FREQ=500E6 PHI=%g DISFRIN=1);' % (Lc, ph_), momentum=name)
    cc = np.array([1e-3, 0.0, 0.0, 0.0, -0.02, 0.0])          # t > 0
    d = trk(e1, cc)[1] - trk(e0, cc)[1]
    # 入口の端の収束だけを取り出すため、同じ L の DISFRIN=1 との差 (本体・出口の端は同じ。FRINGE = 1 は入口だけ)
    tt = -cc[4] * hh / pp
    ww = 2 * math.pi * 500e6 / c_light
    vv = 1e5 / MC2
    formula = cc[0] * vv / (2 * Lc * pp) * (math.sin(ww * tt - (0 - ph_)) + 0.0)
    print('   p = %-7s βγ = %.3f、β = %.4f: Δpx(端) = %.4e、K64 の式 x v/(2 L p0)(sin(ωt - dphis)) = %.4e (比 %.4f)、1/β をかけた式 %.4e' % (
        name, pp, beta, d, formula, d / formula, formula / beta))
