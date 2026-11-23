# -*- coding: utf-8 -*-
"""CAVI の端の収束の 1/β (10_cavity.md §7.1): 軸の近くの Maxwell の場の中のローレンツ力の数値積分 (DOP853) と、追跡 (k64 と v2) の比較。電子。
 場 (c = 1、mc² = 1、長さ m): qE_z = −(v/L) g(z) sin(ωt − dphis) − (r²/4)(∂_z² − ∂_t²)(同)、qE_x = −(x/2)∂_z(qE_z)、qB_y = (x/2)∂_t(qE_z)
 (∇·E = 0 と ∇×B = c⁻²∂_t E の軸の近くの展開)。入口: g = ½(1 + tanh(z/d)) (d = 1e-5 m)、E0 = 1e5 V/m、L = 2e-4 m、500 MHz、PHI = 0.3・1.2。
 (a) 入口の px/x、(b) 端のエネルギー変化の x² の項、(c) δ ≠ 0 の粒子 (粒子の β)、(d) gs の保存した出力 (pysad/tests/script_tests/scripts/v250_cavi_edge_beta__*.gs.txt) との比。
 k64 の追跡は gs と同じ。v2 は 1/β (粒子の β = p/h) を入れた版。
"""
import io
import math
import os
import re
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6, compat

MC2 = 0.51099895e6
W = 2 * math.pi * 500e6 / 299792458.0
L, E0 = 2e-4, 1e5
V = E0 * L
VN = V / MC2


def mk(el, mom):
    ffs, _ = run_sad_string('MOMENTUM = %s;\n%s\nMARK ST = ();\nLINE ONE = (ST CV);\nFFS;\nUSE ONE;\n' % (mom, el), out=io.StringIO())
    return ffs


def trk(ffs, c):
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6.track_particles(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


def lorentz(p0, delta, x0, phi, with_r2=False, d=1e-5):
    dphis = -phi
    pz0 = p0 * (1 + delta)
    beta = pz0 / math.sqrt(1 + pz0 * pz0)
    A = -(VN / L)

    def rhs(t, y):
        x, z, px, pz, E = y
        gam = math.sqrt(1 + px * px + pz * pz)
        s = math.tanh(z / d)
        g, gp, gpp = 0.5 * (1 + s), 0.5 * (1 - s * s) / d, -s * (1 - s * s) / d**2
        sn, cs = math.sin(W * t - dphis), math.cos(W * t - dphis)
        Qz = A * g * sn - ((x * x / 4) * A * (gpp * sn + g * W * W * sn) if with_r2 else 0.0)
        Qx, By = -(x / 2) * A * gp * sn, (x / 2) * A * g * W * cs
        vx, vz = px / gam, pz / gam
        return [vx, vz, Qx - vz * By, Qz, Qz * vz + Qx * vx]
    zs = -40 * d
    ev = lambda t, y: y[1] - L
    ev.terminal = True
    r = solve_ivp(rhs, (zs / beta, zs / beta + 10), [x0, zs, 0.0, pz0, 0.0], method='DOP853', rtol=1e-13, atol=1e-19, events=ev)
    y = r.y[:, -1]
    return y[2] / p0, y[4]


print('== (a) 入口の px/x (FRINGE=1)。数値積分 / k64 / v2 と、比 (追跡/数値) と β ==')
for mom, name in ((2.0e6, '2MEV'), (0.6e6, '600KEV')):
    p0 = mom / MC2
    beta = p0 / math.sqrt(1 + p0 * p0)
    for phi in (0.3, 1.2):
        ffs = mk('CAVI CV=(L=%g VOLT=%g FREQ=500E6 PHI=%g FRINGE=1);' % (L, V, phi), name)
        x0 = 1e-5
        pn, _ = lorentz(p0, 0.0, x0, phi)
        r = {}
        for m in ('k64', 'v2'):
            compat.set_mode(m)
            r[m] = trk(ffs, np.array([x0, 0, 0, 0, 0, 0.0]))[1] / x0
        print('  %-6s β=%.4f PHI=%.1f: 数値 %.5e  k64 %.5e (比 %.4f)  v2 %.5e (比 %.4f)' % (name, beta, phi, pn / x0, r['k64'], r['k64'] / (pn / x0), r['v2'], r['v2'] / (pn / x0)))

print('\n== (b) 端のエネルギー変化の x² の項 (x = 3e-4、E_z の r² の項まで入れた場)。ΔE(x) − ΔE(0) ==')
for mom, name in ((2.0e6, '2MEV'), (0.6e6, '600KEV')):
    p0 = mom / MC2
    ffs = mk('CAVI CV=(L=%g VOLT=%g FREQ=500E6 PHI=0.3 FRINGE=1);' % (L, V), name)
    _, Ea = lorentz(p0, 0.0, 3e-4, 0.3, True)
    _, Eb = lorentz(p0, 0.0, 0.0, 0.3, True)
    out = {}
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        a, b = trk(ffs, np.array([3e-4, 0, 0, 0, 0, 0.0])), trk(ffs, np.zeros(6))
        out[m] = math.sqrt(1 + (p0 * (1 + a[5]))**2) - math.sqrt(1 + (p0 * (1 + b[5]))**2)
    print('  %-6s 数値 %.4e  k64 %.4e (比 %.4f)  v2 %.4e (比 %.4f)' % (name, Ea - Eb, out['k64'], out['k64'] / (Ea - Eb), out['v2'], out['v2'] / (Ea - Eb)))

print('\n== (c) δ ≠ 0 の粒子 (0.6 MeV/c、PHI = 0.3)。比 = 追跡/数値、粒子の β = p0(1+δ)/h ==')
p0 = 0.6e6 / MC2
ffs = mk('CAVI CV=(L=%g VOLT=%g FREQ=500E6 PHI=0.3 FRINGE=1);' % (L, V), '600KEV')
for dl in (-0.05, 0.0, 0.05, 0.15):
    pn, _ = lorentz(p0, dl, 1e-5, 0.3)
    bp = p0 * (1 + dl) / math.sqrt(1 + (p0 * (1 + dl))**2)
    r = {}
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        r[m] = trk(ffs, np.array([1e-5, 0, 0, 0, 0, dl]))[1] / 1e-5
    print('  δ = %+.2f: β = %.4f  k64 比 %.4f  v2 比 %.4f' % (dl, bp, r['k64'] / (pn / 1e-5), r['v2'] / (pn / 1e-5)))

print('\n== (d) gs の保存した出力 (L = 1e-4、VOLT = 10 V、PHI = 0.3、x = 1e-5) の (FRINGE=1 − DISFRIN=1) の px と、k64 の追跡 ==')
for tag, name in (('3gev', '3GEV'), ('2mev', '2MEV'), ('600kev', '600KEV')):
    txt = open(os.path.join(C.ROOT, 'pysad', 'tests', 'script_tests', 'scripts', 'v250_cavi_edge_beta__%s.gs.txt' % tag), encoding='utf-8').read()
    g = {m.group(1): float(m.group(2)) for m in re.finditer(r'(ONE|TWO) PX \{\s*(\S+),', txt)}
    compat.set_mode('k64')
    k = []
    for fr in ('FRINGE=1', 'DISFRIN=1'):
        f = mk('CAVI CV=(L=0.0001 VOLT=10 FREQ=500E6 PHI=0.3 %s);' % fr, name)
        k.append(trk(f, np.array([1e-5, 0, 0, 0, 0, 0.0]))[1])
    pp = {'3GEV': 3e9, '2MEV': 2e6, '600KEV': 0.6e6}[name] / MC2
    print('  %-6s β=%.4f: gs %.12e  k64 %.12e (差 %.1e)  gs / (x v sin(PHI)/(2 L p0)) = %.6f' % (
        name, pp / math.sqrt(1 + pp * pp), g['ONE'] - g['TWO'], k[0] - k[1], (g['ONE'] - g['TWO']) - (k[0] - k[1]),
        (g['ONE'] - g['TWO']) / (1e-5 * (10.0 / MC2) / (2e-4 * pp) * math.sin(0.3) * 1.0)))
