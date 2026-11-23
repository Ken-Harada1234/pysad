# -*- coding: utf-8 -*-
"""偏向磁石の中で鉛直方向の傾き py ≠ 0 の粒子のスピン回転を、厳密な Thomas-BMT の積分 (実験室系の一様磁場、出口は基準軌道に垂直な面) と、
Emittance[] (tturne) の srot・追跡 (track6、刻み 0.01 m) で比べる (K64 の誤りの候補)。
扇形 (E1 = E2 = 0)・K1 = 0・L = 1.2 m・15°・0.5 GeV。使い方: py -3.12 example/16_spin_polarization/spin_bend_vertical_angle.py (数十秒)"""
import io
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
from pysad import compat, emit6d, ffscmd, track6

compat.set_mode('v2')
GS = 0.001159652181280002
ME = 0.51099895e6
E, L, TH = 0.5, 1.2, 0.2617993877991494
src = f'MOMENTUM = {E} GEV;\nDRIFT D1=(L=1);\nBEND E1 = (L={L} ANGLE={TH} E1=0 E2=0);\nMARK ST = ();\nLINE RING = (ST E1 D1);\nFFS;\nUSE RING;\nExit;\n'
ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
els = [e for e in ffs.lat.elements if e.name == 'E1']
p = ffs.momentum()
b = TH / L
rho = L / TH
rex = np.array([-rho * (1 - math.cos(TH)), 0.0, rho * math.sin(TH)])
zex = np.array([-math.sin(TH), 0.0, math.cos(TH)])
xex = np.array([math.cos(TH), 0.0, math.sin(TH)])
yex = np.array([0.0, 1.0, 0.0])


def exact(c, s0, perp_only=False):
    pt = 1 + c[5]
    pv = np.array([c[1], c[3], math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)])
    gam = math.sqrt(1 + (p / ME * pt) ** 2)
    B = np.array([0.0, b, 0.0])

    def rhs(l, u):
        pp, S = u[3:6], u[6:]
        ph = pp / pt
        Bl = np.dot(B, ph) * ph
        Om = (1 + GS * gam) * (B - Bl) / pt if perp_only else ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt
        return np.concatenate([ph, np.cross(pp, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: (u[:3] - rex) @ zex
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2 * L + 1), np.concatenate([[c[0], c[2], 0.0], pv, s0]), events=ev, rtol=1e-12, atol=1e-15, method='DOP853')
    S = sol.y_events[0][0][6:]
    return np.array([S @ xex, S @ yex, S @ zex])


def R_of(f, c):
    return np.array([f(c, np.eye(3)[j]) for j in range(3)]).T


def trk(c, s0, ds=0.01):
    import copy
    cc = np.array(c, float)
    s = s0
    n = int(math.ceil(L / ds))
    for i in range(n):
        sl = copy.copy(els[0])
        sl.L = L / n
        sl.ANGLE = TH / n
        out, s1 = track6.track_particles([sl], np.array([cc]), 1, p, 0.0, False, True, False, spins=[s], radpol=False, circ=100.0)
        cc = out[0]
        s = s1[0]
    return s


c0 = np.array([0, 0, 0, 5e-3, 0, 0.0])
Rx = R_of(exact, c0)
Rp = R_of(lambda c, s: exact(c, s, True), c0)
Rt = R_of(trk, c0)
emit6d.POL_MODE = True
st = emit6d.tturne(els, c0.copy(), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
Re = st.srot[:, :3]
print('py = 5e-3 (他は 0)。スピン回転行列の R[0,1] (z のまわりの回転 −φ_z): 厳密 %.4e | 厳密 (B∥ なし) %.4e | 追跡 %.4e | Emittance[] %.4e' % (Rx[0, 1], Rp[0, 1], Rt[0, 1], Re[0, 1]))
print('  厳密 − (B∥ なし) = (1+a)θ py = %.4e (B∥ = B·v̂ の歳差。K64 は入れない)' % (Rx[0, 1] - Rp[0, 1]))
print('  追跡・Emittance[] − 厳密 = %.4e、(2+a) θ py = %.4e (B∥ なしの分 (1+a)θ py と、さらに θ py)' % (Rt[0, 1] - Rx[0, 1], -(2 + GS) * TH * 5e-3))
for lab, c in (('軸上', np.zeros(6)), ('x = 3 mm', np.array([3e-3, 0, 0, 0, 0, 0.0])), ('δ = 5e-3', np.array([0, 0, 0, 0, 0, 5e-3])), ('y = 5 mm', np.array([0, 0, 5e-3, 0, 0, 0.0]))):
    st = emit6d.tturne(els, c.copy(), emit6d._Phys(p, 100.0), rfsw=False, radcod=False)
    print('%-9s |R_Emittance − R_厳密| = %.2e   |R_追跡(0.01 m) − R_厳密| = %.2e' % (lab, np.abs(st.srot[:, :3] - R_of(exact, c)).max(), np.abs(R_of(trk, c) - R_of(exact, c)).max()))
