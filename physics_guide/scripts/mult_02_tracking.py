# -*- coding: utf-8 -*-
"""MULT の追跡 (08_mult.md): pysad の MULT を、独立な厳密な運動方程式 (DOP853) と比べる。電子 3 GeV、δ = 0。
 (a) ANGLE = 0、K2 (6 極) のみ: 刻み数 ndiv の式と、誤差が ndiv^-2 (2 次) であること (EPS を変えて ndiv を変える)。
 (b) ANGLE = 0、K1 + K2: tsolqu (厳密な線形本体) + 6 極キックの Strang 分割の誤差。
 (c) 反対称 (SK_n) = ROTATE = 90°/(n+1) (ヘルプの SK_n の説明)、DX・DY、ROTATE の向き。
 (d) ANGLE ≠ 0、K2: A_s の展開 (mult_01) から作った厳密な ODE (曲線座標) と pysad の tmulta。
"""
import io
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string

c0 = np.array([0.004, 0.0008, -0.003, -0.0004, 0.0, 0.0])


def make(src_el, L, extra):
    src = ('MOMENTUM = 3GEV;\n%s MP=(L=%.15g %s);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n' % (src_el, L, extra))
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def pm(ffs):
    from pysad import lattice as lat
    el = [e for e in ffs.lat.elements if e.name == 'MP'][0]
    return lat._mult_params(el)


def exact_straight(c, L, k1, kn):
    """H = -√(p²-px²-py²) + K1 (x²-y²)/(2L) + Re Σ_n K_n z^(n+1)/((n+1)! L)  (δ = 0)。K は積分強度、kn = {n: K_n}"""
    def U_grad(x, y):
        z = complex(x, y)
        g = 0j
        for n, K in kn.items():
            g += K / L * z**n / math.factorial(n)         # ∂_x U = Re(g)、∂_y U = -Im(g)
        gx = g.real + k1 / L * x
        gy = -g.imag - k1 / L * y
        return gx, gy

    def rhs(s, v):
        x, px, y, py = v
        pz = math.sqrt(1 - px * px - py * py)
        gx, gy = U_grad(x, y)
        return [px / pz, -gx, py / pz, -gy]
    return solve_ivp(rhs, (0, L), c[:4].tolist(), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1]


print('== (a) ANGLE = 0、K2 のみ。L = 1 m、K2 = 40 /m² (積分強度)。ndiv = 1 + int(√(ampmaxm^(n-1)/(6 (n-1)! eps) |K_n| L))、ampmaxm = 0.05、eps = 0.005 EPS ==')
L, K2 = 1.0, 40.0
ref = exact_straight(c0, L, 0.0, {2: K2})
prev = None
for eps0 in (4.0, 2.0, 1.0, 0.5, 0.25):
    ffs = make('MULT', L, 'K2=%g DISFRIN=1 EPS=%g' % (K2, eps0))
    p = pm(ffs)
    tr = C.track_one(ffs, c0)
    err = np.abs(tr[:4] - ref).max()
    formula = 1 + int(math.sqrt(0.05**1 / 6 / math.factorial(1) / (0.005 * eps0) * K2 * L))
    print('   EPS = %-4g ndiv(pysad) = %3d、式 = %3d   |追跡 - 厳密| = %.3e%s' % (eps0, p['ndiv'], formula, err,
          '' if prev is None else '  (前の誤差との比 %.2f、ndiv の比の 2 乗 %.2f)' % (prev[0] / err, (p['ndiv'] / prev[1])**2)))
    prev = (err, p['ndiv'])
print('   → 誤差は ndiv の 2 乗に反比例 (2 次の Strang 分割)。')

print('\n== (b) ANGLE = 0、K1 + K2 (四極 + 6 極)。L = 1 m、K1 = 0.8、K2 = 40 ==')
ref = exact_straight(c0, L, 0.8, {2: K2})
for eps0 in (2.0, 1.0, 0.5):
    ffs = make('MULT', L, 'K1=0.8 K2=%g DISFRIN=1 EPS=%g' % (K2, eps0))
    p = pm(ffs)
    tr = C.track_one(ffs, c0)
    print('   EPS = %-4g ndiv = %3d  |追跡 - 厳密| = %.3e' % (eps0, p['ndiv'], np.abs(tr[:4] - ref).max()))

print('\n== (c) SK_n = ROTATE 90°/(n+1)、DX・DY、ROTATE の向き (L = 0.3、K_n 積分強度 = 3 × 10^n) ==')
Lc = 0.3
for n in (1, 2, 3):
    K = 3.0 * 10**n
    a = C.track_one(make('MULT', Lc, 'SK%d=%g DISFRIN=1' % (n, K)), c0)
    b = C.track_one(make('MULT', Lc, 'K%d=%g ROTATE=%.15g DEG DISFRIN=1' % (n, K, 90.0 / (n + 1))), c0)
    print('   n = %d: |SK%d - (K%d、ROTATE = %g°)| = %.2e (x, px, y, py, z の最大差)' % (n, n, n, 90.0 / (n + 1), np.abs(a - b).max()))
K = 30.0
dx, dy = 0.002, -0.001
t_d = C.track_one(make('MULT', Lc, 'K2=%g DX=%g DY=%g DISFRIN=1' % (K, dx, dy)), c0)
ci = c0.copy(); ci[0] -= dx; ci[2] -= dy
t_0 = C.track_one(make('MULT', Lc, 'K2=%g DISFRIN=1' % K), ci); t_0[0] += dx; t_0[2] += dy
print('   DX・DY: |MULT(DX, DY) - (x - DX、y - DY で通して DX、DY を足す)| = %.2e' % np.abs(t_d - t_0).max())
th = 0.3
t_r = C.track_one(make('MULT', Lc, 'K2=%g ROTATE=%g DISFRIN=1' % (K, th)), c0)
for sgn in (+1, -1):
    ct, st = math.cos(sgn * th), math.sin(sgn * th)
    ci = c0.copy()
    ci[0], ci[2] = ct * c0[0] - st * c0[2], st * c0[0] + ct * c0[2]
    ci[1], ci[3] = ct * c0[1] - st * c0[3], st * c0[1] + ct * c0[3]
    t0 = C.track_one(make('MULT', Lc, 'K2=%g DISFRIN=1' % K), ci)
    out = t0.copy()
    out[0], out[2] = ct * t0[0] + st * t0[2], -st * t0[0] + ct * t0[2]
    out[1], out[3] = ct * t0[1] + st * t0[3], -st * t0[1] + ct * t0[3]
    print('   ROTATE = %g: 入口で座標を %+g だけ回し、出口で戻す: 差 %.2e' % (th, sgn * th, np.abs(t_r - out).max()))
print('   (差が 1e-15 のほうが tsolrot の向き。入口で (x, y) → (cosθ x - sinθ y, sinθ x + cosθ y) の回転。)')

print('\n== (d) ANGLE ≠ 0、K2 のみ (K0 = K1 = 0)。L = 1 m、ANGLE = 0.2 (ρ = 5 m)、E1 = E2 = 0 (扇形) ==')
print('   K64 の本体 (tbend) は、硬い端の Maxwell フリンジ (y²/ρ の項、05 章 §6) を入れる。H の ODE にはそれが無いので、')
print('   K2 の効果だけ (K2 ≠ 0 の結果 - K2 = 0 の結果) を、追跡と ODE で比べる。')
Lb, ang = 1.0, 0.2
rho = Lb / ang
h = 1.0 / rho
from pysad import track6n_ext
x, y = sp.symbols('x y', real=True)
z = x + sp.I * y


def gk(n, k):
    g = sp.Integer(1)
    for kk in range(1, k + 1):
        g *= sp.Rational((2 * kk - 3) * (2 * kk + 1), 8 * kk * (kk + n + 1))
    return g


KMAX = 6


def ode(K2b, nn=2):
    As = sp.re(sp.expand(sum(gk(nn, k) * K2b / sp.factorial(nn + 1) / Lb * (rho + x)**(sp.Rational(1, 2) - k) * z**(nn + k + 1) / sp.sqrt(rho)
                              for k in range(KMAX + 1)))) if K2b != 0 else 0 * x
    Ax = sp.lambdify((x, y), sp.diff(As, x), 'math')
    Ay = sp.lambdify((x, y), sp.diff(As, y), 'math')

    def rhs(s, v):
        xx, pxx, yy, pyy = v
        pzz = math.sqrt(1 - pxx**2 - pyy**2)
        f = 1 + h * xx
        # H = -(1+hx) pz + h (x + h x²/2) + A_s(x, y)  (δ = 0。曲率に合わせた一様な双極子場。05 章 bend_02 と同じ)
        return [f * pxx / pzz, -(-h * pzz + h * (1 + h * xx) + Ax(xx, yy)), f * pyy / pzz, -Ay(xx, yy)]
    return solve_ivp(rhs, (0, Lb), c0[:4].tolist(), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1]


def trk(K2b, eps0):
    ffs = make('MULT', Lb, 'ANGLE=%g K2=%g E1=0 E2=0 DISFRIN=1 EPS=%g' % (ang, K2b, eps0))
    a = np.zeros((1, 6)); a[0, :] = c0
    return np.asarray(track6n_ext.track_particles_ext(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0][:4]


r0 = ode(0.0)
t0 = trk(0.0, 1.0)
print('   K2 = 0 (扇形の偏向だけ): |追跡 - ODE| = %.2e (Maxwell フリンジの分)' % np.abs(t0 - r0).max())
for K2b in (40.0, 80.0):
    d_ode = ode(K2b) - r0
    print('   K2 = %g: K2 の効果の大きさ (ODE) = %.3e' % (K2b, np.abs(d_ode).max()))
    for eps0 in (4.0, 2.0, 1.0, 0.5, 0.25):
        t1 = trk(K2b, eps0)
        print('     EPS = %-5g |追跡の K2 の効果 - ODE の K2 の効果| = %.3e' % (eps0, np.abs((t1 - t0) - d_ode).max()))
