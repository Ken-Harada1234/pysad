# -*- coding: utf-8 -*-
"""v2.4.0 D2 の独立な基準 (SAD・pysad の光学・追跡のコードを使わない)。

硬い端の直線の極 (ANGLE = 0、K0 [rad]、長さ L) の厳密な写像。磁場は一様 (B = K0/L)、端の面は実験室の直線座標で
入口 z = X tan(e1)、出口 z = L − X tan(e2) (X = その面を横切る点の x)。基準軌道に垂直な面 z = 0 の (x, px, δ) から z = L の面の (x, px, δ, z) を返す。
運動の式は z を独立変数にして dpx/dz = −F (F = K0/L)。pz = √((1+δ)² − px² − py²)。円弧は解析解、面との交点は Brent 法。
z は経路の遅れ ∫(p/pz − 1) dz (K64 の z とは符号が逆)。放射つきは `pole_radiative` (RK4。放射の反作用 dδ/dz = −c_r (1+δ)² F² (p/pz))。
"""
import math

import numpy as np
from scipy.optimize import brentq

RE = 2.8179403262e-15
MC2 = 0.51099895e6


def pole_exact(x0, px0, d0, L, K0, e1, e2, py0=0.0, y0=0.0):
    """放射なし。戻り: (x, px, y, py, zlag, δ)"""
    F = K0 / L
    p = 1.0 + d0
    t1, t2 = math.tan(e1), math.tan(e2)
    q2 = p * p - py0 * py0
    q = math.sqrt(q2)
    pz0 = math.sqrt(q2 - px0 * px0)
    zf = x0 * t1 / (1 - px0 / pz0 * t1)           # 入口の面までの直進
    x = x0 + px0 / pz0 * zf
    y = y0 + py0 / pz0 * zf
    zl = zf * (p / pz0 - 1.0)

    def state(z):
        dz = z - zf
        px = px0 - F * dz
        pz = math.sqrt(q2 - px * px)
        xx = x + (pz - pz0) / F
        iy = -(math.asin(px / q) - math.asin(px0 / q)) / F
        ip = -(p / F) * (math.asin(px / q) - math.asin(px0 / q)) - dz
        return xx, px, pz, y + py0 * iy, zl + ip

    z = brentq(lambda zz: zz - (L - state(zz)[0] * t2), zf - 1.0, L + 1.0, xtol=1e-16, rtol=1e-15)
    xx, px, pz, yy, zlag = state(z)
    dz = L - z
    return np.array([xx + px / pz * dz, px, yy + py0 / pz * dz, py0, zlag + dz * (p / pz - 1.0), d0])


def pole_radiative(c, L, K0, e1, e2, E_eV, nz=400):
    """放射つき (x, px, δ, zlag) の 4 次元写像 (RK4 + 面との交点の反復)。c = (x, px, δ, zlag)。電子 (γ = E/mc²)、
    放射の反作用 dδ/dz = −(2/3) r_e γ³ (1+δ)² F² (p/pz) (古典放射の平均、光子の離散性なし)。戻り (x, px, δ, zlag)"""
    g = E_eV / MC2
    cr = (2 / 3) * RE * g ** 3
    F = K0 / L
    t1, t2 = math.tan(e1), math.tan(e2)

    def f(y):
        x, px, d, z = y
        pz = math.sqrt((1 + d) ** 2 - px * px)
        dl = (1 + d) / pz
        dd = -cr * (1 + d) ** 2 * F * F * dl
        return np.array([px / pz, -F + px * dd / (1 + d), dd, dl - 1.0])

    x, px, d, z = c
    pz = math.sqrt((1 + d) ** 2 - px * px)
    zf = x * t1 / (1 - (px / pz) * t1)
    y = np.array([x + (px / pz) * zf, px, d, z + zf * ((1 + d) / pz - 1)])
    zend = L
    for _ in range(8):
        h = (zend - zf) / nz
        yy = y.copy()
        for _i in range(nz):
            k1 = f(yy)
            k2 = f(yy + h / 2 * k1)
            k3 = f(yy + h / 2 * k2)
            k4 = f(yy + h * k3)
            yy = yy + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        znew = L - yy[0] * t2
        if abs(znew - zend) < 1e-15:
            break
        zend = znew
    x, px, d, z = yy
    pz = math.sqrt((1 + d) ** 2 - px * px)
    dz = L - zend
    return np.array([x + (px / pz) * dz, px, d, z + dz * ((1 + d) / pz - 1.0)])


def numeric_jacobian(fn, c, idx, h=1e-6):
    n = len(idx)
    J = np.zeros((n, n))
    for j, i in enumerate(idx):
        cp = list(c)
        cm = list(c)
        cp[i] += h
        cm[i] -= h
        J[:, j] = (np.asarray(fn(cp)) - np.asarray(fn(cm))) / (2 * h)
    return J
