# -*- coding: utf-8 -*-
"""B165 の独立な基準 (pysad のコードを使わない): Maxwell 方程式を満たす硬い端 (薄いフリンジ) の扇形 BEND の中の Thomas-BMT の数値積分。
場は スカラーポテンシャル Φ = −b y G(r) (G = f(入口面からの距離) f(−出口面からの距離)、f = (1 + tanh(d/w))/2) から B = b (G ŷ + y ∇G) (y の 1 次)。
∫B_s ds (端の縦磁場) と本体の B∥ が両方入り、全体の ∫B·dr = 0。w → 0 (4e-3・1e-3・5e-4 m) で結果は 1e-6 以内で動かない。"""
import math

import numpy as np
from scipy.integrate import solve_ivp

GS = 0.001159652181280002; ME = 0.51099895e6
def make(E, L, TH, w):
    p = None
    return None
def maxwell_R(Ep, L, TH, PY, w=2e-3, zs=1.0e-2, dp=0.0, perp_only=False):
    """Ep: momentum eV/c.  Returns 3x3 spin map in (x,y,z) of entrance / exit frames"""
    b = TH / L; rho = L / TH
    C = np.array([-rho, 0.0, 0.0])
    zin = np.array([0, 0, 1.0]); zex = np.array([-math.sin(TH), 0, math.cos(TH)])
    xex = np.array([math.cos(TH), 0, math.sin(TH)]); yex = np.array([0, 1.0, 0])
    f = lambda d: 0.5 * (1 + math.tanh(d / w))
    fp = lambda d: 0.5 / w * (1 - math.tanh(d / w) ** 2)
    pt = 1 + dp
    gam = math.sqrt(1 + (Ep / ME * pt) ** 2)
    def field(r):
        din = r[2]; dout = (r - C) @ zex
        G = f(din) * f(-dout)
        gG = fp(din) * f(-dout) * zin - f(din) * fp(-dout) * zex
        return b * (G * np.array([0, 1.0, 0]) + r[1] * gG)
    pz = math.sqrt(pt * pt - PY * PY)
    def rhs(l, u):
        r, pp, S = u[:3], u[3:6], u[6:]
        B = field(r); ph = pp / pt
        Bl = (B @ ph) * ph
        Om = ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt if not perp_only else (1 + GS * gam) * (B - Bl) / pt
        return np.concatenate([ph, np.cross(pp, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: (u[:3] - C) @ zex - zs
    ev.terminal = True
    def one(s0):
        u0 = np.concatenate([[0, -PY * zs / pz, -zs], [0, PY, pz], s0])
        sol = solve_ivp(rhs, (0, 4 * L + 1), u0, events=ev, rtol=1e-12, atol=1e-15, method='DOP853', max_step=w / 4)
        S = sol.y_events[0][0][6:]
        return np.array([S @ xex, S @ yex, S @ zex])
    return np.array([one(np.eye(3)[j]) for j in range(3)]).T
