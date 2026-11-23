# -*- coding: utf-8 -*-
"""薄い多極のハミルトニアンとキック (07_thin_multipole.md §3〜§5、式 (7.1)〜(7.12)) を確かめる。

(a) 多極のポテンシャル V_k = K/(k+1)! Re (x - i y)^(k+1) (k = 1: 四極、2: 6 極、…): 調和関数 (∇²V = 0)、k = 1…5 で中面の場 V_x(y=0) = K x^k/k!、
    K が b_y の x に関する k 階微分 (積分強度の定義)、複素形 V_x + i V_y = (K/k!)(x - i y)^k、すなわち Δ(px + i py) = -(K/k!)(x - i y)^k。
(b) exp(:f:) (f = -V、座標だけの関数): ポアソン括弧の級数が 2 項で止まる (k = 2: 6 極、k = 3: 8 極)。シンプレクティック性 ([x, px'] = 1、[px', py'] = 0)。
(c) pysad の薄い多極 (L = 0) の追跡と、閉じた式 (正準運動量で px -= K/k! Re cx) の比較 (δ ≠ 0 でも)。
(d) フィードダウン (DX、DY): 閉じた式を展開して、中心がずれた 6 極のキック = 四極 + 二極 (+ スキュー) になること、pysad の追跡との比較。
"""
import io
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import track6

x, y, px, py, K = sp.symbols('x y p_x p_y K', real=True)
zb = x - sp.I * y                                           # z̄ = x - i y


def Vk(k):
    return sp.expand(K / sp.factorial(k + 1) * sp.re(sp.expand(zb**(k + 1))))


def PB(f, g):
    """[f, g] = Σ (∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q)、q = (x, y)、p = (px, py)"""
    return sp.expand(sum(sp.diff(f, q) * sp.diff(g, pp) - sp.diff(f, pp) * sp.diff(g, q) for q, pp in ((x, px), (y, py))))


print('== (a) 多極のポテンシャル V_k = K/(k+1)! Re (x - i y)^(k+1) ==')
names = {1: '四極', 2: '6 極', 3: '8 極', 4: '10 極 (DECA)', 5: '12 極 (DODECA)'}
for k in range(1, 6):
    V = Vk(k)
    lap = sp.simplify(sp.diff(V, x, 2) + sp.diff(V, y, 2))
    by0 = sp.simplify(sp.diff(V, x).subs(y, 0))                           # b_y(y=0) = ∂V/∂x
    bx0 = sp.simplify(sp.diff(V, y).subs(y, 0))                           # b_x(y=0) = ∂V/∂y (中面で 0)
    kth = sp.simplify(sp.diff(by0, x, k))
    cplx = sp.simplify(sp.expand(sp.diff(V, x) + sp.I * sp.diff(V, y) - K / sp.factorial(k) * zb**k))
    print('  k = %d (%s): ∇²V = %s、b_y(y=0) = %s、b_x(y=0) = %s、∂^k b_y/∂x^k = %s、(V_x + i V_y) - (K/k!) z̄^k = %s' % (k, names[k], lap, by0, bx0, kth, cplx))
print('  → V_x + i V_y = (K/k!)(x - i y)^k (z̄ = x - i y)。キック Δpx = -V_x、Δpy = -V_y から Δ(px + i py) = -(K/k!) z̄^k、すなわち '
      'Δpx = -(K/k!) Re z̄^k、Δpy = -(K/k!) Im z̄^k (K64 の px -= aki Re cx、py -= aki Im cx、cx = (x - i y)^k、aki = K/k!)')
for k in (2, 3):
    V = Vk(k)
    ck = sp.expand(zb**k)
    print('  k = %d: -V_x - (-(K/k!) Re cx) = %s、-V_y - (-(K/k!) Im cx) = %s' % (k, sp.simplify(-sp.diff(V, x) + K / sp.factorial(k) * sp.re(ck)), sp.simplify(-sp.diff(V, y) + K / sp.factorial(k) * sp.im(ck))))

print('\n== (b) exp(:f:) の級数: f = -V (座標だけの関数) ==')
for k in (2, 3):
    f = -Vk(k)
    t1x, t1y = PB(f, px), PB(f, py)
    t2x, t2y = PB(f, t1x), PB(f, t1y)
    print('  k = %d: [f, x] = %s、[f, y] = %s、[f, px] = %s' % (k, PB(f, x), PB(f, y), sp.factor(t1x)))
    print('         [f, py] = %s、[f,[f, px]] = %s、[f,[f, py]] = %s → 級数は 2 項で止まる' % (sp.factor(t1y), t2x, t2y))
    pxn, pyn = px + t1x, py + t1y
    print('         exp(:f:) px = px + %s、exp(:f:) py = py + %s' % (sp.factor(t1x), sp.factor(t1y)))
    print('         シンプレクティック: [x, px\'] = %s、[y, py\'] = %s、[x, py\'] = %s、[y, px\'] = %s、[px\', py\'] = %s' % (
        PB(x, pxn), PB(y, pyn), PB(x, pyn), PB(y, pxn), sp.simplify(PB(pxn, pyn))))

print('\n== (c) pysad の薄い多極 (L = 0) の追跡と閉じた式 ==')


def thin_ffs(name, key, ak, L=0.0, extra='', momentum='3GEV'):
    src = ('MOMENTUM = %s;\n%s MP=(L=%.15g %s=%.15g %s);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n' % (momentum, name, L, key, ak, extra))
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def closed_kick(c, K_, k, dx=0.0, dy=0.0):
    """Δ(px, py) = -K/k! (Re, Im) (x - dx - i (y - dy))^k (canonical)。x, y は磁石の中心から測る"""
    z = complex(c[0] - dx, -(c[2] - dy))**k
    out = np.array(c, float)
    out[1] -= K_ / math.factorial(k) * z.real
    out[3] -= K_ / math.factorial(k) * z.imag
    return out


rng = np.random.default_rng(7)
for name, key, k in (('SEXT', 'K2', 2), ('OCT', 'K3', 3), ('DECA', 'K4', 4), ('DODECA', 'K5', 5)):
    K_ = 40.0 * 10**(k - 2)                    # 積分強度 K (薄いので ak = K)。大きさは 6 極 40、8 極 400、… (1/m^k)
    ffs = thin_ffs(name, key, K_)
    worst = 0.0
    for dl in (0.0, 0.01):
        for _ in range(3):
            c = np.concatenate([rng.normal(0, 3e-3, 4), [0.0, dl]])
            c[4] = 0.0
            tr = C.track_one(ffs, c)
            ref = closed_kick(c, K_, k)
            worst = max(worst, np.abs(tr - ref).max())
    print('  %s (k = %d、K = %g、L = 0): pysad の追跡と閉じた式 px -= K/k! Re (x - i y)^k の最大差 %.2e (δ = 0 と 0.01、乱数 6 点)' % (name, k, K_, worst))

print('\n== (d) フィードダウン: 中心がずれた 6 極 ==')
Ks, dxs, dys = sp.symbols('K d_x d_y', real=True)
Vs = Ks / 6 * sp.re(sp.expand((x - dxs - sp.I * (y - dys))**3))
dpx = sp.expand(-sp.diff(Vs, x))
dpy = sp.expand(-sp.diff(Vs, y))
print('  Δpx = %s' % sp.collect(dpx, (x, y)))
print('  Δpy = %s' % sp.collect(dpy, (x, y)))
quad_n = sp.expand(dpx.coeff(x, 1).subs(y, 0))
print('  展開: 6 極項 -(K/2)(x^2 - y^2)、四極項: Δpx ⊃ K d_x x - K d_y y、Δpy ⊃ -K d_y x - K d_x y、二極項: Δpx = -(K/2)(d_x^2 - d_y^2)、Δpy = K d_x d_y')
print('  = 四極 (K1 = -K d_x: Δpx = -K1 x、Δpy = +K1 y) + スキュー四極 (強さ K d_y) + 二極 (キック)。ずれ d_x が大きいほど 6 極が四極として働く (フィードダウン)')
chk = sp.simplify(dpx - (-Ks / 2 * (x**2 - y**2) + Ks * (dxs * x - dys * y) - Ks / 2 * (dxs**2 - dys**2)))
chk2 = sp.simplify(dpy - (Ks * x * y - Ks * (dys * x + dxs * y) + Ks * dxs * dys))
print('  展開式との差: Δpx %s、Δpy %s' % (chk, chk2))
K2v, dxv, dyv = 40.0, 1.5e-3, -0.8e-3
ffs_d = thin_ffs('SEXT', 'K2', K2v, 0.0, 'DX=%.15g DY=%.15g' % (dxv, dyv))
worst = 0.0
for _ in range(5):
    c = np.concatenate([rng.normal(0, 3e-3, 4), [0.0, 0.002]])
    c[4] = 0.0
    tr = C.track_one(ffs_d, c)
    ref = closed_kick(c, K2v, 2, dxv, dyv)
    worst = max(worst, np.abs(tr - ref).max())
print('  pysad の SEXT (K2 = %g、DX = %.1f mm、DY = %.1f mm) の追跡と、(x - DX, y - DY) を使う閉じた式の最大差 %.2e' % (K2v, dxv * 1e3, dyv * 1e3, worst))
tr0 = C.track_one(ffs_d, np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
print('  軌道が x = y = 0 (磁石から (-DX, -DY) ずれ) のときの二極のキック (pysad): Δpx = %.4e、Δpy = %.4e (式: -(K/2)(dx^2 - dy^2) = %.4e、K dx dy = %.4e)' % (
    tr0[1], tr0[3], -K2v / 2 * (dxv**2 - dyv**2), K2v * dxv * dyv))
