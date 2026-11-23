# -*- coding: utf-8 -*-
"""四極の非線形フリンジの母関数 (06_edge_fringe.md §6): リー変換 exp(:f:) の 2 次までの打ち切りと、z の更新。

(a) sympy: f = p·V(x) (V は x の関数のベクトル) は点変換の生成子。母関数 G = H0 + (∂H0/∂x ∂H0/∂px + ∂H0/∂y ∂H0/∂py)/2 から
    x2 = ∂F/∂p2 = x + V + (V·∇)V/2 が出る。四極の (4.29) の x2 = x(1 + a/3 + b + t) と一致する。
(b) 数値: f = (4.27) の厳密な流れ (ハミルトン方程式、τ = 0 → 1、DOP853) と K64 の式 (pysad の ttfrin 4 次) の差は K³ に比例 (K を 2 倍で 8 倍)。
(c) z の更新: z2 = z1 + ∂G/∂p。K64 は -(p2·(dx1 + t x))/p (dx1 に t を含むので、t の項が 2 回)。
    「-p2·(x2 - x1)/p」(t の項が 1 回) とは違う。
"""
import math
import warnings
warnings.filterwarnings("ignore")
from types import SimpleNamespace
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

x, y, px, py, p, K = sp.symbols('x y p_x p_y p K', real=True)
c = K / (12 * p)
Vx = c * (x**3 + 3 * x * y**2)
Vy = -c * (3 * x**2 * y + y**3)
H0 = px * Vx + py * Vy
G = H0 + sp.Rational(1, 2) * (sp.diff(H0, x) * sp.diff(H0, px) + sp.diff(H0, y) * sp.diff(H0, py))

print('== (a) 母関数 G から出る点変換 ==')
x2 = sp.expand(x + sp.diff(G, px))
y2 = sp.expand(y + sp.diff(G, py))
VV_x = Vx * sp.diff(Vx, x) + Vy * sp.diff(Vx, y)
VV_y = Vx * sp.diff(Vy, x) + Vy * sp.diff(Vy, y)
print('  x2 - (x + V_x + (V·∇)V_x/2) = %s、y2 - (y + V_y + (V·∇)V_y/2) = %s' % (
    sp.simplify(x2 - (x + Vx + VV_x / 2)), sp.simplify(y2 - (y + Vy + VV_y / 2))))
a = K * x**2 / (4 * p)
b = K * y**2 / (4 * p)
t = (a - b)**2 / 6
print('  (4.29) との差: x2 - x(1 + a/3 + b + t) = %s、y2 - y(1 - a - b/3 + t) = %s' % (
    sp.simplify(x2 - x * (1 + a / 3 + b + t)), sp.simplify(y2 - y * (1 - a - b / 3 + t))))
print('  V = (Δx1, Δy1) = ∂f/∂(px, py) (f = H0 = (4.27))。G の第 2 項は (V·∇)V/2 を作るためのもので、1 次 (K の 1 次の Maxwell) には効かない。')

print('\n== (b) 厳密な流れ exp(:f:) と K64 の式 (ttfrin) ==')
f_ = H0
fun = sp.lambdify((x, px, y, py, p, K), [sp.diff(f_, px), -sp.diff(f_, x), sp.diff(f_, py), -sp.diff(f_, y), sp.diff(f_, p)], 'math')


def flow(c0, Kv):
    def rhs(tau, v):
        return fun(v[0], v[1], v[2], v[3], 1.0, Kv)
    r = solve_ivp(rhs, (0, 1), [c0[0], c0[1], c0[2], c0[3], c0[4]], rtol=1e-13, atol=1e-18, method='DOP853')
    return r.y[:, -1]


def ttf(c0, Kv):
    q = SimpleNamespace(c=np.array([c0[0], c0[1], c0[2], c0[3], c0[4], 0.0], float))
    t6._ttfrin_n(q, 4, Kv, 1.0)
    return q.c[:5].copy()


c0 = np.array([0.004, 0.003, -0.003, 0.002, 0.0])
prev = None
print('  入射 (x, px, y, py, z) = %s、δ = 0' % c0.tolist())
for Kv in (4e3, 8e3, 1.6e4):
    d = np.abs(flow(c0, Kv) - ttf(c0, Kv))
    # 1 次の項の大きさ
    first = np.abs(ttf(c0, Kv) - c0)
    print('  K = %.1e: |K64 - exact flow| = %s   (|変化| = %s)' % (Kv, ' '.join('%.2e' % v for v in d), ' '.join('%.2e' % v for v in first)))
    if prev is not None:
        print('           差の比 (前の K の差で割る) = %s' % ' '.join('%.2f' % (v / w if w > 0 else float('nan')) for v, w in zip(d, prev)))
    prev = d
print('  → 差は K³ に比例 (比 8)。K64 の式は exp(:f:) の 2 次までと一致し、3 次以上は違う (1 次の Maxwell の導出の範囲では意味のない次数)。')

print('\n== (c) z の更新 ==')
Gf = sp.lambdify((x, y, px, py, p, K), G, 'math')
dGdp = sp.lambdify((x, y, px, py, p, K), sp.diff(G, p), 'math')
Kv = 1.6e4
for dl in (0.0, 0.01):
    pv = 1 + dl
    # 新しい運動量 (K64 の式で)
    q = SimpleNamespace(c=np.array([0.004, 0.003, -0.003, 0.002, 0.0, dl], float))
    t6._ttfrin_n(q, 4, Kv, 1.0)
    X2, PX2, Y2, PY2, Z2 = q.c[:5]
    # q 変数 (px/p) なので、正準運動量に直す
    PX2c, PY2c = PX2 * pv, PY2 * pv
    # ∂G/∂p: G は p2 (新しい運動量) の関数。K64 の x1 = 入射の x
    zG = dGdp(0.004, -0.003, PX2c, PY2c, pv, Kv)
    t_ = ((Kv * 0.004**2 - Kv * 0.003**2) / (4 * pv))**2 / 6
    dx1 = X2 - 0.004
    dy1 = Y2 - (-0.003)
    z_one = -(PX2c * dx1 + PY2c * dy1) / pv
    z_two = -(PX2c * (dx1 + t_ * 0.004) + PY2c * (dy1 + t_ * (-0.003))) / pv
    print('  δ = %.2f: K64 の Δz = %.6e   ∂G/∂p = %.6e   -p2·(x2-x1)/p = %.6e   -p2·(dx1 + t x)/p = %.6e' % (dl, Z2, zG, z_one, z_two))
print('  → K64 の Δz は ∂G/∂p と一致する。「-p2·(x2 - x1)/p」は t の項が 1 回で、K の 2 次 (x²/(4p) の 2 乗) だけ違う。')
