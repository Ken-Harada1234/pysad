# -*- coding: utf-8 -*-
"""アンジュレータ UND (16_undulator_beambeam.md): A_s の模型、スライスの写像。

(a) sympy: A = (ff, gg, 0) から B = ∇×A を作り、∇×B = 0 (真空の磁場) と軸上の磁場 (By = F0 cos kz s、Bx = G0 cos(kz s + φ))。
    条件 Ky² = Kz² + Kx²、Qx² = Kz² + Qy²。
(b) スライスの写像 (undulator.f の陰的な 1 次のシンプレクティック Euler、v2 の dpds を 1 回だけ掛ける形) と、同じ H の
    DOP853 の積分の比較。スライス数を 2 倍にしたときの誤差の比 (1 次)。
(c) 軸上の軌道の最大角 F0/Kz = K/γ と、追跡の最大の px。
(d) UND と DRIFT の z の更新の違い (UND は −dv L の項を持たない)。
"""
import io
import math
import numpy as np
import _common as C

C.start(__file__)
print('pysad:', C.pysad_src())

# ---------------------------------------------------------------- (a)
import sympy as sp
x, y, s = sp.symbols('x y s', real=True)
F0, G0, Kz, Kx, Qy, ph = sp.symbols('F0 G0 Kz Kx Qy phi', positive=True)
Ky = sp.sqrt(Kz**2 + Kx**2)
Qx = sp.sqrt(Kz**2 + Qy**2)
ff = (F0 / Kz * sp.cos(Kx * x) * sp.cosh(Ky * y) * sp.sin(Kz * s)
      - Qy * G0 / Qx / Kz * sp.sinh(Qx * x) * sp.sin(Qy * y) * sp.sin(Kz * s + ph))
gg = (-G0 / Kz * sp.cos(Qy * y) * sp.cosh(Qx * x) * sp.sin(Kz * s + ph)
      + Kx * F0 / Ky / Kz * sp.sinh(Ky * y) * sp.sin(Kx * x) * sp.sin(Kz * s))
A = sp.Matrix([ff, gg, 0])
X = (x, y, s)


def curl(V):
    return sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], s), sp.diff(V[0], s) - sp.diff(V[2], x), sp.diff(V[1], x) - sp.diff(V[0], y)])


Bv = curl(A)
print('== (a) B = ∇×A (A_s = 0) ==')
print('∇×B の各成分 (sympy で 0 になるか。数値代入 10 点の最大値):')
cB = curl(Bv)
rng = np.random.default_rng(3)
vals = []
for k in range(10):
    sub = {x: rng.normal() * 0.01, y: rng.normal() * 0.01, s: rng.normal(), F0: 0.3, G0: 0.2, Kz: 62.8, Kx: 5.0, Qy: 7.0, ph: 0.7}
    vals.append(max(abs(float(c.subs(sub).evalf(30))) for c in cB))
print('  max |∇×B| = %.2e (B の大きさ ~ 0.3)' % max(vals))
print('  ∇·B (恒等的に 0): %s' % sp.simplify(sum(sp.diff(Bv[i], X[i]) for i in range(3))))
Bax = [sp.simplify(Bv[i].subs({x: 0, y: 0})) for i in range(3)]
print('軸上 (x = y = 0) の B = (Bx, By, Bz) =', Bax)

# ---------------------------------------------------------------- (b)
from pysad import und as U

pgev = 3e9
p0 = pgev / 0.51099895e6      # p0 = βγ (単位 m c)
h0 = math.sqrt(1 + p0**2)
Bt = 1.0                       # [T]
Brho = pgev / 299792458.0     # [T m]
F0n = Bt / Brho
G0n = 0.4 * Bt / Brho
L = 1.0
npole = 10
Kzn = 2 * math.pi * npole / L
Kxn = 4.0
Qyn = 5.0
dphn = 0.3
Kyn = math.sqrt(Kzn**2 + Kxn**2)
Qxn = math.sqrt(Kzn**2 + Qyn**2)
print('\n== (b) スライスの写像と DOP853 (3 GeV、By = %.1f T (F0 = %.5f /m)、Bx = %.2f T、L = %.0f m、%d 周期、Kx = %.0f、Qy = %.0f、φ = %.1f) ==' % (
    Bt, F0n, 0.4 * Bt, L, npole, Kxn, Qyn, dphn))
from scipy.integrate import solve_ivp


def Aff(xv, yv, sv):
    sinkz = math.sin(Kzn * sv)
    sinkzp = math.sin(Kzn * sv + dphn)
    f = (F0n / Kzn * math.cos(Kxn * xv) * math.cosh(Kyn * yv) * sinkz
         - Qyn * G0n / Qxn / Kzn * math.sinh(Qxn * xv) * math.sin(Qyn * yv) * sinkzp)
    g = (-G0n / Kzn * math.cos(Qyn * yv) * math.cosh(Qxn * xv) * sinkzp
         + Kxn * F0n / Kyn / Kzn * math.sinh(Kyn * yv) * math.sin(Kxn * xv) * sinkz)
    return f, g


def dA(xv, yv, sv, h=1e-7):
    fx = (Aff(xv + h, yv, sv)[0] - Aff(xv - h, yv, sv)[0]) / (2 * h)
    fy = (Aff(xv, yv + h, sv)[0] - Aff(xv, yv - h, sv)[0]) / (2 * h)
    gx = (Aff(xv + h, yv, sv)[1] - Aff(xv - h, yv, sv)[1]) / (2 * h)
    gy = (Aff(xv, yv + h, sv)[1] - Aff(xv, yv - h, sv)[1]) / (2 * h)
    return fx, fy, gx, gy


delta = 2e-3
pn = 1 + delta


def rhs(sv, w):
    xv, yv, px, py, zv = w
    f, g = Aff(xv, yv, sv)
    fx, fy, gx, gy = dA(xv, yv, sv)
    ux, uy = px - f, py - g               # 運動学的な運動量
    # H = (ux² + uy²)/(2 pn)。dx/ds = ux/pn、dpx/ds = -∂H/∂x = (ux fx + uy gx)/pn
    return [ux / pn, uy / pn, (ux * fx + uy * gx) / pn, (ux * fy + uy * gy) / pn, -(ux**2 + uy**2) / (2 * pn**2)]


# 入口: 外の運動量 (px_out/pn) から。入口の s = 0 で ff = gg = 0 (sin(0) の項と sin(φ) の項: φ ≠ 0 なので gg ≠ 0)
x0, px_out, y0, py_out = 3e-5, 1e-4, -2e-5, 5e-5
a0 = Aff(x0, y0, 0.0)
w0 = [x0, y0, px_out * pn, py_out * pn, 0.0]
sol = solve_ivp(rhs, [0, L], w0, method='DOP853', rtol=1e-13, atol=1e-16)
xe, ye, pxe, pye, ze = sol.y[:, -1]
ae = Aff(xe, ye, L)
print('入口の A = (%.3e, %.3e)、出口の A = (%.3e, %.3e) (φ ≠ 0 なので 0 でない: 端でポテンシャルの飛び = 蹴り)' % (a0[0], a0[1], ae[0], ae[1]))
ref = np.array([xe, pxe / pn, ye, pye / pn, ze])
print('DOP853 の出口 (x, px_canon/pn, y, py_canon/pn, z) =', ref)
prev = None
for ns in (500, 1000, 2000, 4000, 8000):
    ds = L / ns
    S = np.zeros((1, 16))
    S[0, :6] = [x0, px_out, y0, py_out, 0.0, delta]
    U.und_seg1(0, S, 0, ns, ns, F0n, G0n, dphn, Kzn, Kxn, Kyn, Qxn, Qyn, ds, False, False, np.array([1], dtype=np.int64), p0, 1.0, 0.0)
    # und_seg1 の出口は正準運動量 / pn (j1 = nslice のとき A を引かない)
    out = np.array([S[0, 0], S[0, 1], S[0, 2], S[0, 3], S[0, 4]])
    err = np.abs(out - ref)
    print('slice = %5d: |Δx| = %.2e、|Δpx| = %.2e、|Δy| = %.2e、|Δpy| = %.2e、|Δz| = %.2e%s' % (
        ns, err[0], err[1], err[2], err[3], err[4], '' if prev is None else '   (px の誤差の前との比 %.2f)' % (prev / err[1])))
    prev = err[1]

# ---------------------------------------------------------------- (c)
print('\n== (c) 平面のアンジュレータ (G0 = 0、Kx = 0) の軸上の軌道の最大角 ==')
F1 = 0.3
Kz1 = 2 * math.pi * 10 / 1.0
S = np.zeros((1, 16))
S[0, :6] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
mx = 0.0
ds = 1.0 / 4000
for j in range(1, 4001, 20):
    S[0, :6] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    U.und_seg1(0, S, 0, j, 4000, F1, 0.0, 0.0, Kz1, 0.0, Kz1, Kz1, 0.0, ds, False, False, np.array([1], dtype=np.int64), p0, 1.0, 0.0)
    # 途中のとき (j1 < nslice) は px = 正準 − A を引いた力学的な運動量/pn
    mx = max(mx, abs(S[0, 1]))
print('F0 = %.2f /m、Kz = %.3f /m: F0/Kz = %.5e (= K/γ)、追跡の最大の px (力学的) = %.5e' % (F1, Kz1, F1 / Kz1, mx))
print('  K = (F0/Kz) γ = %.4f、周期 %.2f m、B = F0 Bρ = %.3f T' % (F1 / Kz1 * h0, 1.0 / 10, F1 * Brho))

# ---------------------------------------------------------------- (d)
print('\n== (d) UND (磁場なし) と DRIFT の z の更新 ==')
from pysad.ffscmd import run_sad_string
from pysad import compat
res = {}
for mode, name, el in [(m_, n_, e_) for m_ in ('k64', 'v2') for n_, e_ in (('UND', 'UND U=(L=1 KX=0.0 QY=0.0 POLE=10 SLICE=4);'), ('DRIFT', 'DRIFT U=(L=1);'))]:
    compat.set_mode(mode)
    src = ('MOMENTUM = 3GEV;\n' + el + '\nMARK ST=();\nLINE BL=(ST U);\nFFS;\nUSE BL;\nTRPT; INS; NORAD; NOFLUC;\n$FORM="18.14";\n'
           'x0={1e-4};px0={2e-4};y0={0};py0={0};z0={0};dp0={1e-3};fl0={1};\nr=TrackParticles[{1,{x0,px0,y0,py0,z0,dp0,fl0}},3,1];\nPrint["R ",r[[2]]];\nExit;\n')
    b = io.StringIO()
    run_sad_string(src, out=b)
    for l in b.getvalue().splitlines():
        if l.startswith('R'):
            import re as _re
            nums = [float(t) for t in _re.findall(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?', l)]
            res[(mode, name)] = nums
            print('[%s]' % mode, name, 'x, px, y, py, z =', nums[:5])
compat.set_mode(None)
dvL = 1e-3 * (2 + 1e-3) / math.sqrt(1 + (p0 * (1 + 1e-3))**2) / (math.sqrt(1 + (p0 * (1 + 1e-3))**2) + (1 + 1e-3) * h0)
for mode in ('k64', 'v2'):
    dz = res[(mode, 'UND')][4] - res[(mode, 'DRIFT')][4]
    print('[%s] z(UND) - z(DRIFT) = %.3e m。−dv L (dv = -δ(2+δ)/h1/(h1+pn h0)、L = 1) の大きさ = %.3e m' % (mode, dz, dvL))
