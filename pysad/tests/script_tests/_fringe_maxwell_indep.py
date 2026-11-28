# -*- coding: utf-8 -*-
"""BEND の端のソフトフリンジの y に 1 次の py のキックの運動量 (1+δ) 依存の独立な基準 (pysad のコードを使わない)。docs/31 §9.1。

場 (x に依らない、Laplace を満たす): By = B0[b - y^2/2 b'' + ...]、Bz = B0[y b' - y^3/6 b'''' + ...]、b(z) = erf の立ち上がり (幅 σ)。
単位: p0 = 1、e B0 / p0 = g = 1 (曲率半径 1)。運動量 p = 1+δ (kinetic、p0 で規格化)。独立変数は z。
   dx/dz = px/pz, dy/dz = py/pz, dpx/dz = py Bz/pz - By, dpy/dz = -px Bz/pz,  pz = sqrt(p^2 - px^2 - py^2)
ゲージ: A_y = 0 (B_x = 0 の場) なので py は正準運動量に等しい。フリンジの後の一様場 (By = B0, Bz = 0) では py は保存される。
測るもの: Δpy = -κ y の κ(δ) を y = ±1e-4 の差分で求め、κ(δ)/κ(0) = (1+δ)^-k の k を取り出す。K = ∫ b(1-b) dz とすると κ(0) = -K g^2。
結果 (σ = 0.04, 0.02, 0.01, 0.005): k = 1.003〜1.004, 1.001, 1.0002, 1.0000〜1.0001 (δ = -0.10〜+0.10 のどこでも)、
κ(0)/K = -1.0017, -1.0004, -1.0001, -0.9999 (符号は Δpy = +K g^2 y)。
使い方: py -3.12 _fringe_maxwell_indep.py (数秒)
"""
from math import erf, sqrt, pi, exp, log

from scipy.integrate import solve_ivp, quad


def run(sigma, delta, y0, zspan=9.0):
    """(x, y, px, py) のフリンジの後の値。初期値は x = px = py = 0、y = y0"""
    s = sigma
    p = 1.0 + delta

    def b(z):
        return 0.5 * (1.0 + erf(z / (sqrt(2.0) * s)))

    def phi(z):
        return exp(-z * z / (2 * s * s)) / (s * sqrt(2 * pi))

    def f(z, u):
        x, y, px, py = u
        ph = phi(z)
        d1 = ph                                  # b'
        d2 = -z / s**2 * ph                      # b''
        d3 = (z * z / s**4 - 1.0 / s**2) * ph    # b'''
        By = b(z) - 0.5 * y * y * d2
        Bz = y * d1 - y**3 / 6.0 * d3
        pz = sqrt(p * p - px * px - py * py)
        return [px / pz, py / pz, py * Bz / pz - By, -px * Bz / pz]
    sol = solve_ivp(f, (-zspan * s, zspan * s), [0.0, y0, 0.0, 0.0], method='DOP853', rtol=1e-13, atol=1e-16)
    return sol.y[:, -1]


def kappa(sigma, delta, y0=1e-4):
    """Δpy = -κ y の κ (y の奇関数の 1 次の係数)"""
    up = run(sigma, delta, y0)
    dn = run(sigma, delta, -y0)
    return -(up[3] - dn[3]) / (2 * y0)


def fringe_integral(sigma):
    """K = ∫ b (1 - b) dz (erf の立ち上がり)"""
    return sigma * quad(lambda z: 0.5 * (1 + erf(z / sqrt(2))) * (1 - 0.5 * (1 + erf(z / sqrt(2)))), -12, 12)[0]


if __name__ == '__main__':
    deltas = [-0.10, -0.05, -0.02, 0.02, 0.05, 0.10]
    for sigma in (0.04, 0.02, 0.01, 0.005):
        K = fringe_integral(sigma)
        k0 = kappa(sigma, 0.0)
        row = ['sigma=%.3f kappa0/K=%.6f' % (sigma, k0 / K)]
        for d in deltas:
            row.append('delta=%+.2f k=%.4f' % (d, -log(kappa(sigma, d) / k0) / log(1 + d)))
        print(' | '.join(row), flush=True)
