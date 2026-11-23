# -*- coding: utf-8 -*-
r"""04_quad.md §7.2 の検証: 非線形フリンジと Maxwell の場の端の ODE の差 (相対 0.1〜3%) の中身を調べる。

quad_05_fringe.py (b) の同じ設定 (K = 3 /m^2、入口、S = 0.06 m、tanh の端、幅 d) で、
(1) ODE の許容誤差 rtol を 1e-12, 1e-13, 1e-14 と変えて、差 (ODE - モデル) が動くか (= 積分の雑音か) を見る。
(2) 差を振幅 s (粒子の初期座標の倍率) の線形成分 α s と 3 次成分 β s^3 に分ける。R(s) = α s + β s^3 とすると、
    R(2s) - 2 R(s) = 6 β s^3 (線形成分が消える)。線形成分は、線形フリンジ (F1 あり) の残り (K の 2 次の線形な項) と比べる。
"""
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import solve_ivp
from types import SimpleNamespace
import _common as C

C.start(__file__)
from pysad import track6 as t6


def maxwell_rhs(s, v, K0, d, sgn):
    x, y, px, py, z = v
    t = math.tanh(s / d)
    sech2 = 1.0 - t * t
    K = K0 * (1 + sgn * t) / 2
    K1 = sgn * K0 * sech2 / (2 * d)
    K2 = -sgn * K0 * t * sech2 / d**2
    bx = K * y - K2 * (3 * x * x * y + y**3) / 12
    by = K * x - K2 * (x**3 + 3 * x * y * y) / 12
    bs = K1 * x * y
    pz = math.sqrt(max(1.0 - px * px - py * py, 1e-12))
    return [px / pz, py / pz, (py * bs - pz * by) / pz, (pz * bx - px * bs) / pz, 1.0 - 1.0 / pz]


def rhs_lin(s, v, K0, d, sgn):
    x, y, px, py, z = v
    K = K0 * (1 + sgn * math.tanh(s / d)) / 2
    pz = math.sqrt(max(1.0 - px * px - py * py, 1e-12))
    return [px / pz, py / pz, -K * x, K * y, 1.0 - 1.0 / pz]


def kin_body(c, S, K, rtol):
    def rhs(s, v):
        x, y, px, py, z = v
        pz = math.sqrt(1.0 - px * px - py * py)
        return [px / pz, py / pz, -K * x, K * y, 1.0 - 1.0 / pz]
    return solve_ivp(rhs, (0, S), list(c), method='DOP853', rtol=rtol, atol=1e-17).y[:, -1]


def kin_drift(c, S):
    x, y, px, py, z = c
    pz = math.sqrt(1.0 - px * px - py * py)
    return np.array([x + px / pz * S, y + py / pz * S, px, py, z + (1.0 - 1.0 / pz) * S])


def k64_fringe(c, K, sign):
    q = SimpleNamespace(c=np.array([c[0], c[2], c[1], c[3], c[4], 0.0]))
    t6._ttfrin4(q, sign * K, 1.0)
    o = q.c
    return np.array([o[0], o[2], o[1], o[3], o[4]])


def lin_fringe(c, a, b):
    q = SimpleNamespace(c=np.array([c[0], c[2], c[1], c[3], c[4], 0.0]))
    t6._tqlfr(q, a, b)
    return np.array([q.c[0], q.c[2], q.c[1], q.c[3], q.c[4]])


K0, S = 3.0, 0.06


def resid(scale, d, rtol):
    """(ODE の差 dex, モデルの差 dmod, 線形フリンジの残り)"""
    c0 = np.array([0.010, -0.008, 0.005, -0.003, 0.0]) * scale
    f1 = -K0 * (math.pi * d)**2 / 24
    kw = dict(method='DOP853', rtol=rtol, atol=1e-17)
    ex_m = solve_ivp(maxwell_rhs, (-S, S), c0, args=(K0, d, +1.0), **kw).y[:, -1]
    ex_l = solve_ivp(rhs_lin, (-S, S), c0, args=(K0, d, +1.0), **kw).y[:, -1]
    cc = kin_drift(c0, S)
    cl = lin_fringe(cc, f1, 0.0)
    m_lin = kin_body(cl, S, K0, rtol)
    m_all = kin_body(k64_fringe(cl, K0, +1.0), S, K0, rtol)
    return ex_m - ex_l, m_all - m_lin, ex_l - m_lin


print('== (1) ODE の許容誤差 rtol を変えたときの 差 R = max|Δ(ODE) - Δ(モデル)| [m] (r = 2 mm = 倍率 0.2) ==')
print('   d [mm] |  rtol=1e-12    rtol=1e-13    rtol=1e-14  | 3 つの最大のばらつき  R の大きさ  ばらつき/R')
for d in (0.0005, 0.001, 0.002, 0.005):
    Rs = []
    for rt in (1e-12, 1e-13, 1e-14):
        dex, dmod, _ = resid(0.2, d, rt)
        Rs.append(dex - dmod)
    nz = [np.abs(r).max() for r in Rs]
    spread = max(np.abs(Rs[i] - Rs[j]).max() for i in range(3) for j in range(i))
    print('   %5.1f  |  %.3e   %.3e   %.3e |   %.2e          %.3e   %.1e' % (d * 1e3, nz[0], nz[1], nz[2], spread, nz[2], spread / nz[2]))
print('   (ばらつきは R の 1% 以下: R は積分の雑音ではない)')
print('\n== (2) 差 R(s) = α s + β s^3 + γ s^5 の分離 (rtol = 1e-14。成分ごとに s = 0.05, 0.1, 0.2 の 3 点で解き、s = 0.4 (r = 4 mm) で検算) ==')
print('   α: 線形成分、β: 3 次成分、γ: 5 次成分。表は r = 2 mm (s = 0.2) での各成分の最大値 [m]')
print('   d [mm] | Δ(モデル) | α s (線形) | 線形フリンジの残り | β s^3 (3 次) | γ s^5 (5 次) | s=0.4 での検算誤差/R(0.4)')
scs = (0.05, 0.1, 0.2)
for d in (0.0005, 0.001, 0.002, 0.005):
    Rv, Dm, Lr = {}, {}, {}
    for sc in scs + (0.4,):
        dex, dmod, lr = resid(sc, d, 1e-14)
        Rv[sc], Dm[sc], Lr[sc] = dex - dmod, dmod, lr
    A = np.array([[s_, s_**3, s_**5] for s_ in scs])
    Y = np.array([Rv[s_] for s_ in scs])
    coef = np.linalg.solve(A, Y)
    pred04 = np.array([0.4, 0.4**3, 0.4**5]) @ coef
    chk = np.abs(pred04 - Rv[0.4]).max() / np.abs(Rv[0.4]).max()
    t = np.array([0.2, 0.2**3, 0.2**5])[:, None] * coef
    print('   %5.1f  |  %.3e |  %.3e |  %.3e       |  %.3e  |  %.3e  |  %.1e' % (
        d * 1e3, np.abs(Dm[0.2]).max(), np.abs(t[0]).max(), np.abs(Lr[0.2]).max(), np.abs(t[1]).max(), np.abs(t[2]).max(), chk))
print('   読み方 (d <= 1 mm では 5 次成分が 3 次成分と同じ桁。d = 0.5 mm では 5 次が上回る): 線形成分 α は 0 (線形フリンジの残りは Δ(ODE) と Δ(モデル) の両方に入って打ち消す)。')
print('   3 次成分 β s^3 は d に比例し (下の β s^3 / d)、5 次成分 γ s^5 は d に反比例する (下の γ s^5 d)。')
print('   d <= 1 mm では 5 次成分が 3 次成分と同じ桁以上で、r を半分にすると 3 次成分は 1/8、5 次成分は 1/32 になるので、')
print('   「差/Δ(モデル)」は振幅によって変わる。d >= 2 mm では 3 次成分が主で、比は振幅によらない。')
print('   d [mm] | β s^3 / d [m/m] | (β s^3 /Δ(モデル)) / (2 K d) | γ s^5 d [m^2]')
for d in (0.0005, 0.001, 0.002, 0.005):
    Rv, Dm = {}, {}
    for sc in scs:
        dex, dmod, lr = resid(sc, d, 1e-14)
        Rv[sc], Dm[sc] = dex - dmod, dmod
    A = np.array([[s_, s_**3, s_**5] for s_ in scs])
    coef = np.linalg.solve(A, np.array([Rv[s_] for s_ in scs]))
    b3 = np.abs(0.2**3 * coef[1]).max()
    g5 = np.abs(0.2**5 * coef[2]).max()
    print('   %5.1f  |   %.3e      |   %.3f                       |   %.3e' % (d * 1e3, b3 / d, b3 / np.abs(Dm[0.2]).max() / (2 * K0 * d), g5 * d))
