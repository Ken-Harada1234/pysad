# -*- coding: utf-8 -*-
"""ソレノイドの端の非線形フリンジ (06_edge_fringe.md §8): H = -(Bz/(8 p² Bρ)) p_φ p_r。

(a) K64 の tsfrin (pysad の _tsfrin) は、f = -(b/8) p_φ p_r (p = 1) の厳密な流れ (ハミルトン方程式の DOP853) と一致するか。
(b) Maxwell を満たす端の磁場 (b_z = b g - r² b g''/4、b_r = -r b g'/2 + r³ b g'''/16、g = (1 + tanh(s/d))/2) の中のローレンツ力の積分と、
    モデル「自由な運動 → 硬い端 (正準運動量が連続) → 一様な Bz の中」の差を、tsfrin の写像と比べる (d → 0 に外挿)。
"""
import math
import warnings
warnings.filterwarnings("ignore")
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

print('== (a) tsfrin と f = -(b/8) p_φ p_r の流れ ==')
b = 1.7


def tsf(c0, dbz):
    q = SimpleNamespace(c=np.array(list(c0) + [0.0], float))
    t6._tsfrin(q, dbz)
    return q.c[:5].copy()


def flow(c0, dbz):
    def rhs(tau, v):
        x, px, y, py = v[:4]
        pphi = x * py - y * px
        pr = x * px + y * py
        f_px = -(dbz / 8) * (-y * pr + pphi * x)
        f_py = -(dbz / 8) * (x * pr + pphi * y)
        f_x = -(dbz / 8) * (py * pr + pphi * px)
        f_y = -(dbz / 8) * (-px * pr + pphi * py)
        return [f_px, -f_x, f_py, -f_y, -(dbz / 8) * 0.0]
    r = solve_ivp(rhs, (0, 1), [c0[0], c0[1], c0[2], c0[3], 0.0], rtol=1e-13, atol=1e-18, method='DOP853')
    v = r.y[:, -1]
    # z: dz/dτ = ∂f/∂p (f = -(b/(8p²)) p_φ p_r → ∂f/∂p = +(b/(4p³)) p_φ p_r)。p_φ p_r は τ で変わるので積分
    return v


def flow_z(c0, dbz):
    def rhs(tau, v):
        x, px, y, py, z = v
        pphi = x * py - y * px
        pr = x * px + y * py
        return [-(dbz / 8) * (-y * pr + pphi * x), (dbz / 8) * (py * pr + pphi * px),
                -(dbz / 8) * (x * pr + pphi * y), (dbz / 8) * (-px * pr + pphi * py), (dbz / 4) * pphi * pr]
    return solve_ivp(rhs, (0, 1), [c0[0], c0[1], c0[2], c0[3], 0.0], rtol=1e-13, atol=1e-18, method='DOP853').y[:, -1]


c0 = [0.003, 0.02, -0.002, 0.015, 0.0]
a_, b_ = tsf(c0, b), flow_z(c0, b)
print('  入射 (x, px, y, py) = %s、b = %.1f /m、δ = 0' % (c0[:4], b))
print('  tsfrin = %s' % np.array2string(a_, precision=10))
print('  流れ   = %s' % np.array2string(b_, precision=10))
print('  差の最大 = %.2e (x, px, y, py, z。厳密な流れと一致)' % np.abs(a_ - b_).max())

print('\n== (b) Maxwell の磁場の中のローレンツ力の積分 ==')
bz0 = 1.0
x0, y0 = 3e-4, -2e-4
px0, py0 = 3e-3, 2e-3     # δ = 0、p = 1 (規格化)
Pz = lambda Px, Py: math.sqrt(1 - Px * Px - Py * Py)


def lorentz(d, s0, s1, state, uniform=False, bb=bz0):
    def rhs(s, v):
        x, y, Px, Py = v
        r2 = x * x + y * y
        if uniform:
            gg, g1, g2, g3 = 1.0, 0.0, 0.0, 0.0
        else:
            t = math.tanh(s / d)
            gg, g1, g2, g3 = 0.5 * (1 + t), 0.5 * (1 - t * t) / d, -t * (1 - t * t) / d**2, -(1 - t * t) * (1 - 3 * t * t) / d**3
        bzv = bb * (gg - r2 * g2 / 4)
        br_over_r = bb * (-g1 / 2 + r2 * g3 / 16)       # b_r / r
        bxv, byv = br_over_r * x, br_over_r * y
        pz = Pz(Px, Py)
        # dP/ds = P × b / pz
        dPx = (Py * bzv - pz * byv) / pz
        dPy = (pz * bxv - Px * bzv) / pz
        return [Px / pz, Py / pz, dPx, dPy]
    return solve_ivp(rhs, (s0, s1), state, rtol=1e-13, atol=1e-20, method='DOP853', max_step=d / 5 if not uniform else 1.0).y[:, -1]


def canon(x, y, Px, Py, bb):
    """正準運動量 p = P + a、a = (-b y/2, b x/2)"""
    return Px + (-bb * y / 2), Py + bb * x / 2


def run(d, with_f, bz0=1.0):
    L = 8 * d
    global_bz = bz0
    if with_f is None:    # ローレンツ力の積分
        r = lorentz(d, -L, L, [x0, y0, px0, py0], bb=bz0)
        bend = bz0 * 0.5 * (1 + math.tanh(L / d))
        pxc, pyc = canon(r[0], r[1], r[2], r[3], bend)
        return np.array([r[0], pxc, r[1], pyc])
    # モデル: 自由な運動 -L → 0 (b = 0 なので直線)
    pz = Pz(px0, py0)
    x1, y1 = x0 + px0 / pz * L, y0 + py0 / pz * L
    c1 = [x1, px0, y1, py0, 0.0]
    if with_f:
        c1 = list(tsf(c1, bz0)[:4]) + [0.0]
    # 一様な Bz の中を 0 → L (正準運動量が連続: P = p - a)
    Px, Py = c1[1] - (-bz0 * c1[2] / 2), c1[3] - bz0 * c1[0] / 2
    r = lorentz(d, 0.0, L, [c1[0], c1[2], Px, Py], uniform=True, bb=bz0)
    pxc, pyc = canon(r[0], r[1], r[2], r[3], bz0)
    return np.array([r[0], pxc, r[1], pyc])


def odd(d, with_f):
    """b → -b で符号の変わる部分 (O(b))。b の偶数次 (線形な端の集束 F1 型など) は消える"""
    return (run(d, with_f, +1.0) - run(d, with_f, -1.0)) / 2


res = {}
for d in (1e-3, 5e-4, 2.5e-4):
    ode = odd(d, None)
    m0 = odd(d, False)
    m1 = odd(d, True)
    res[d] = (ode - m0, m1 - m0)
    print('  d = %.2f mm: b の奇数部分 (ODE - 硬い端のみ) = %s' % (d * 1e3, ' '.join('%.4e' % v for v in (ode - m0))))
    print('               b の奇数部分 (tsfrin の効果)      = %s' % ' '.join('%.4e' % v for v in (m1 - m0)))
r1 = [(4 * res[a][0] - res[b][0]) / 3 for a, b in ((5e-4, 1e-3), (2.5e-4, 5e-4))]
ext_f = res[2.5e-4][1]
print('  d → 0 (d² の項を消す、d = 0.5/1 と 0.25/0.5): ODE = %s' % ' '.join('%.4e' % v for v in r1[0]))
print('                                              ODE = %s' % ' '.join('%.4e' % v for v in r1[1]))
print('                                              tsfrin = %s' % ' '.join('%.4e' % v for v in ext_f))
print('  比 (x, px, y, py) = %s' % ' '.join('%.4f' % v for v in r1[1] / ext_f))
print('  (振幅 0.3 mm、角 3 mrad、|b| = 1 /m (Bz/Bρ)、y ≪ d の展開が使える範囲。窓は ±8d)')
