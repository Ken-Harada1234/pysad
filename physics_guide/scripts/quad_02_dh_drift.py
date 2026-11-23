# -*- coding: utf-8 -*-
"""運動学的な残り ΔH の解 (04_quad.md 式 (4.15)〜(4.17)) を確かめる。

ΔH(px, py, p) = p - sqrt(p^2 - px^2 - py^2) - (px^2 + py^2)/(2p)   は座標を含まないので px, py, p は一定で、
x' = ∂ΔH/∂px、y' = ∂ΔH/∂py、z' = ∂ΔH/∂p (= ∂ΔH/∂δ) が定数 (「非線形ドリフト」D)。
(a) 偏微分を複素ステップ微分で数値的に作り、手で導いた式 (4.16)(4.17) と比べる。
(b) D を DOP853 で正準方程式から積分した結果と、pysad の D (track6._qb_drift、K64 の tsolqu の「x+=px r」の行) を比べる。
(c) ΔH の大きさ (傾きの 4 乗のオーダー) を表にする。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6


def dH(px, py, p):
    return p - np.sqrt(p * p - px * px - py * py) - (px * px + py * py) / (2 * p)


def cs_grad(px, py, p, h=1e-30):
    """複素ステップ微分 (桁落ちなし): (∂ΔH/∂px, ∂ΔH/∂py, ∂ΔH/∂p)"""
    f = lambda a, b, c: dH(a, b, c)
    return (f(px + 1j * h, py + 0j, p + 0j).imag / h, f(px + 0j, py + 1j * h, p + 0j).imag / h, f(px + 0j, py + 0j, p + 1j * h).imag / h)


print('== (a) ΔH の偏微分: 複素ステップ微分 と 手で導いた式 ==')
print('   x\' = px (1/pz - 1/p)、y\' = py (1/pz - 1/p)、z\' = 1 - p/pz + (px^2+py^2)/(2 p^2)  (pz = sqrt(p^2-px^2-py^2))')
for px, py, p in ((1e-3, -2e-3, 1.0), (0.01, 0.02, 1.003), (0.05, -0.03, 0.97)):
    pz = math.sqrt(p * p - px * px - py * py)
    hand = (px * (1 / pz - 1 / p), py * (1 / pz - 1 / p), 1 - p / pz + (px * px + py * py) / (2 * p * p))
    num = cs_grad(px, py, p)
    print('  (px, py, p) = (%.3f, %.3f, %.3f): 最大差 %.2e (相対 %.1e)' % (px, py, p, max(abs(a - b) for a, b in zip(hand, num)),
          max(abs(a - b) / abs(b) for a, b in zip(hand, num))))

print('\n== (b) D の解: DOP853 (正準方程式) と pysad の D の比較 (長さ 0.3 m。座標は規格化した運動量 px/p) ==')
print('   ここの D は ΔH の流れだけ (H2 の項は含まない)。pysad の D は z に dv の項を含まない (dv は本体 B 側)')
worst = 0.0
for px_n, py_n, dp in ((1e-3, -2e-3, 0.0), (0.01, 0.02, 0.003), (0.05, -0.03, -0.03)):
    p = 1.0 + dp
    px, py = px_n * p, py_n * p                       # 正準運動量
    Ln = 0.3

    def rhs(s, v):
        g = cs_grad(px, py, p)
        return [g[0], g[1], g[2]]
    sol = solve_ivp(rhs, (0, Ln), [0.02, -0.01, 0.0], method='DOP853', rtol=1e-13, atol=1e-16)
    ref = sol.y[:, -1]
    x, y, z = t6._qb_drift(0.02, px_n, -0.01, py_n, 0.0, Ln)
    err = max(abs(x - ref[0]), abs(y - ref[1]), abs(z - ref[2]))
    worst = max(worst, err)
    print('  px/p = %.3f, py/p = %.3f, δ = %+.3f: x 変化 %.4e、z 変化 %.4e、pysad との差 %.1e' % (px_n, py_n, dp, ref[0] - 0.02, ref[2], err))
print('  最大差 %.2e' % worst)

print('\n== (c) ΔH の効果の大きさ (長さ 1 m、p = 1、py = 0 の粒子; 傾き a = px/p) ==')
print('   傾き a [mrad]   x の変化 a(1/pz-1) [m]   (傾きの 3 乗/2)    z の変化 [m]      (-3/8 a^4)')
for a in (1e-3, 5e-3, 1e-2, 3e-2, 1e-1):
    x, y, z = t6._qb_drift(0.0, a, 0.0, 0.0, 0.0, 1.0)
    print('   %8.1f        %14.4e        %14.4e   %14.4e   %14.4e' % (a * 1e3, x, a**3 / 2, z, -3 * a**4 / 8))
