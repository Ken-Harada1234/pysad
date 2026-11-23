# -*- coding: utf-8 -*-
"""v2.3.0 U1 の試験の補助: 一様な磁場指数の扇形 BEND (複合機能) の独立な数値積分 (pysad・SAD のコードを使わない)。

実験室系 (円筒座標 r, y, θ) でローレンツ力を直接積分する。Bρ0 = 1 の単位、運動量 p = 1 + δ。
 dr/dθ = r ur/uθ, dy/dθ = r uy/uθ, dur/dθ = -r By/p + uθ, duy/dθ = r Br/p, uθ = sqrt(1 - ur² - uy²)
 (ur・uy は方向余弦、θ は扇形磁石の中心角。SAD の x = r - ρ0、px = ur、y = y、py = uy。)
磁場は中央面の By(r, y=0) = g(r) で決まる (軸対称・静磁場・∇×B = 0)。
 By = g - y²/2 (g'' + g'/r) + ...、Br = y g' - y³/6 (g'' + g'/r)' + ... (ラプラス方程式 φ_rr + φ_r/r + φ_yy = 0 の解)。
 磁場の模型 (K1 = g'(ρ0)/Bρ0 は共通、g''(ρ0) だけが違う):
  'lin'  : g = h + K1 (r - ρ0)             B の勾配が x に一様 (elegant CSBEND K2 = 0。g'' = 0)
  'log'  : g = h + (K1/h) ln(r/ρ0)         r ∂B/∂r が一様 = 曲線座標の勾配 (1+hx)∂b/∂x が一様 (g'' = -h K1)
  'pow'  : g = h (ρ0/r)^n                  磁場指数 n が一様 (g'' = -(n+1) h K1)
 端は硬い端・フリンジなし (扇形磁石、E1 = E2 = 0)。
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def make_field(model, rho0, k1, nidx=None):
    h = 1.0 / rho0

    if model == 'lin':
        def f(r, y):
            g = h + k1 * (r - rho0)
            return g - 0.5 * y * y * k1 / r, y * k1 + y**3 * k1 / (6 * r * r)
    elif model == 'log':
        def f(r, y):
            return h + (k1 / h) * math.log(r / rho0), y * k1 * rho0 / r
    elif model == 'pow':
        n = -k1 * rho0 * rho0 if nidx is None else nidx

        def f(r, y):
            g = h * (rho0 / r)**n
            gp = -n * g / r
            gpp = n * (n + 1) * g / r**2
            lg = gpp + gp / r
            return g - 0.5 * y * y * lg, y * gp   # y³ の項は省略 (線形・y の 2 次までの比較用)
    else:
        raise ValueError(model)
    return f


def sector_map(field, rho0, theta0, p, state, rtol=1e-13, atol=1e-15):
    """扇形磁石 1 台 (中心角 theta0) の写像。state = (x, px, y, py) (x = r - rho0)。"""
    x, px, y, py = state

    def rhs3(th, s):
        r, ur, yy, uy = s
        ut = math.sqrt(1.0 - ur * ur - uy * uy)
        by, br = field(r, yy)
        return [r * ur / ut, -r * by / p + ut, r * uy / ut, r * br / p]
    sol = solve_ivp(rhs3, (0.0, theta0), [rho0 + x, px, y, py], method='DOP853', rtol=rtol, atol=atol)
    r, ur, yy, uy = sol.y[:, -1]
    return np.array([r - rho0, ur, yy, uy])


def closed_orbit_radius(field, rho0, p):
    return brentq(lambda r: r * field(r, 0.0)[0] - p, 0.5 * rho0, 1.5 * rho0, xtol=1e-15, rtol=1e-15)


def tunes(model, rho0, k1, nmag, delta, nidx=None, eps=1e-7):
    """16 台 (nmag 台) の扇形磁石のリングの (νx, νy) (オフモメンタムの円軌道の周りの線形写像)"""
    fld = make_field(model, rho0, k1, nidx)
    theta0 = 2 * math.pi / nmag
    p = 1.0 + delta
    rc = closed_orbit_radius(fld, rho0, p)
    xc = rc - rho0
    s0 = np.array([xc, 0.0, 0.0, 0.0])
    m = np.zeros((4, 4))
    for j in range(4):
        d = np.zeros(4)
        d[j] = eps
        fp = sector_map(fld, rho0, theta0, p, s0 + d)
        fm = sector_map(fld, rho0, theta0, p, s0 - d)
        m[:, j] = (fp - fm) / (2 * eps)
    out = []
    for a, b in ((0, 1), (2, 3)):
        tr = 0.5 * (m[a, a] + m[b, b])
        out.append(nmag * math.acos(tr) / (2 * math.pi))
    return out[0], out[1], m


def dnu_ddelta(model, rho0, k1, nmag, dmax=0.004, npts=5, nidx=None):
    ds = np.linspace(-dmax, dmax, npts)
    t = np.array([tunes(model, rho0, k1, nmag, d, nidx)[:2] for d in ds])
    cx = np.polyfit(ds, t[:, 0], 1)[0]
    cy = np.polyfit(ds, t[:, 1], 1)[0]
    return cx, cy, ds, t
