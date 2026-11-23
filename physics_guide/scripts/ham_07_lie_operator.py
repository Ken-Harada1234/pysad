# -*- coding: utf-8 -*-
"""ポアソン括弧とリー演算子 exp(:f:) の初歩 (C_poisson_lie.md) を sympy で確かめる。

規約: [f, g] = Σ_i (∂f/∂q_i ∂g/∂p_i - ∂f/∂p_i ∂g/∂q_i)、:f: g = [f, g]、exp(:f:) = Σ :f:^n / n!。
(a) 正準方程式: q' = [q, H]、p' = [p, H]。任意の g(q, p) について dg/ds = [g, H]。
(b) 例 1: ドリフト f = -s p^2/2。
(c) 例 2: f = -s (k x^2 + p^2)/2 (集束四極の線形マップ)。
(d) 例 3: 6 極のキック f = -k2 x^3/3 (級数が 2 項で止まる)。
(e) 合成の順序: exp(:f:) exp(:g:) の座標への作用は「f の変換 → g の変換」の順。
(f) [:f:, :g:] = :[f, g]: (ヤコビ恒等式) と BCH 公式の 3 次まで。
(g) 2 次の f ⇔ 線形写像 exp(sF)。シンプレクティック性 (行列と、(d) の写像のポアソン括弧)。
"""
import sympy as sp
import _common as C

C.start(__file__)

x, p, s, k, k2, eps, w = sp.symbols('x p s k k2 epsilon w', positive=True)


def PB(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x))


def colon(f, g, n=1):
    for _ in range(n):
        g = PB(f, g)
    return g


def lie_exp(f, g, nmax=14):
    """exp(:f:) g = Σ :f:^n g / n!。:f:^n g が 0 になったところで止める。戻り: (和、項の数)"""
    tot = 0
    term = g
    n = 0
    while True:
        tot += term / sp.factorial(n)
        n += 1
        term = PB(f, term)
        if term == 0 or n > nmax:
            break
    return sp.expand(tot), n


print('== (a) ポアソン括弧と正準方程式 ==')
print('  [x, p] = %s、[x, x] = %s、[p, p] = %s' % (PB(x, p), PB(x, x), PB(p, p)))
V = sp.Function('V')(x)
Hm = p**2 / 2 + V
print('  H = p^2/2 + V(x): [x, H] = %s (= ∂H/∂p = x\')、[p, H] = %s (= -∂H/∂x = p\')' % (PB(x, Hm), PB(p, Hm)))
g = sp.Function('g')(x, p)
dg = sp.diff(g, x) * PB(x, Hm) + sp.diff(g, p) * PB(p, Hm)      # 鎖則: dg/ds = g_x x' + g_p p'
print('  dg/ds (鎖則) - [g, H] = %s' % sp.simplify(dg - PB(g, Hm)))
print('  → dg/ds = [g, H] = -:H: g。解は g(s) = exp(-s :H:) g(0) (H が s によらないとき)')

print('\n== (b) 例 1: ドリフト f = -s p^2/2 ==')
f = -s * p**2 / 2
print('  :f: x = %s、:f:^2 x = %s' % (colon(f, x), colon(f, x, 2)))
print('  :f: p = %s' % colon(f, p))
print('  exp(:f:) x = %s (項数 %d)、exp(:f:) p = %s' % (lie_exp(f, x)[0], lie_exp(f, x)[1], lie_exp(f, p)[0]))

print('\n== (c) 例 2: f = -s (k x^2 + p^2)/2、w = sqrt(k) ==')
f = -s * (k * x**2 + p**2) / 2
print('  :f: x = %s、:f: p = %s' % (colon(f, x), colon(f, p)))
print('  :f:^2 x = %s、:f:^3 x = %s、:f:^4 x = %s' % (colon(f, x, 2), colon(f, x, 3), colon(f, x, 4)))
print('  :f:^2 p = %s、:f:^3 p = %s' % (colon(f, p, 2), colon(f, p, 3)))
N = 12
ex = lie_exp(f, x, N)[0]
ep = lie_exp(f, p, N)[0]
tx = sp.series(x * sp.cos(sp.sqrt(k) * s) + p * sp.sin(sp.sqrt(k) * s) / sp.sqrt(k), s, 0, N).removeO()
tp = sp.series(-sp.sqrt(k) * x * sp.sin(sp.sqrt(k) * s) + p * sp.cos(sp.sqrt(k) * s), s, 0, N).removeO()
print('  級数 (s の %d 次未満) - [x cos(√k s) + p sin(√k s)/√k] の級数 = %s' % (N, sp.simplify(sp.expand(sp.series(ex, s, 0, N).removeO() - tx))))
print('  級数 - [-√k x sin(√k s) + p cos(√k s)] の級数 = %s' % sp.simplify(sp.expand(sp.series(ep, s, 0, N).removeO() - tp)))
print('  最初の項: exp(:f:) x = x + s p - (k s^2/2) x - (k s^3/6) p + ... = %s' % sp.series(ex, s, 0, 4).removeO())
print('  → x → x cos(√k s) + p sin(√k s)/√k、p → -√k x sin(√k s) + p cos(√k s) (集束四極の線形写像)')
print('  k < 0 なら cos → cosh、sin → sinh (発散)。k -> 0 の極限はドリフト (b)')

print('\n== (d) 例 3: 6 極のキック f = -k2 x^3/3 (k2 は係数。SAD の K2 との関係は未確認) ==')
f = -k2 * x**3 / 3
print('  :f: x = %s、:f: p = %s、:f:^2 p = %s' % (colon(f, x), colon(f, p), colon(f, p, 2)))
rx, nx = lie_exp(f, x)
rp, npp = lie_exp(f, p)
print('  exp(:f:) x = %s (項数 %d)、exp(:f:) p = %s (項数 %d)  → 級数が止まる (多項式)' % (rx, nx, rp, npp))
print('  確認: ハミルトニアン H = k2 x^3/3 の 1 ステップ (長さ 1 の薄いキック) は Δp = -∂H/∂x = -k2 x^2')
newx, newp = rx, rp
print('  新しい (x, p) のポアソン括弧 [x_new, p_new] = %s (1 なら正準変換 = シンプレクティック)' % sp.simplify(PB(newx, newp)))

print('\n== (e) 合成の順序 ==')
f = -s * p**2 / 2          # ドリフト
g = -k * x**2 / 2          # 薄いキック p -> p - k x
fx, _ = lie_exp(f, x)
fp, _ = lie_exp(f, p)
gx, _ = lie_exp(g, x)
gp_, _ = lie_exp(g, p)
print('  exp(:f:): x -> %s, p -> %s (f = ドリフト)、exp(:g:): x -> %s, p -> %s (g = キック)' % (fx, fp, gx, gp_))
# 演算子の積 exp(:f:) exp(:g:) を座標関数 x, p に作用 (右の演算子が先に関数に作用)
Ux, _ = lie_exp(f, lie_exp(g, x)[0])
Up, _ = lie_exp(f, lie_exp(g, p)[0])
print('  exp(:f:) exp(:g:) x = %s、exp(:f:) exp(:g:) p = %s' % (Ux, Up))
# 粒子の追跡: 先にドリフト f、後でキック g
tx_ = fx
tp_ = fp
trk_x = gx.subs({x: tx_, p: tp_}, simultaneous=True)
trk_p = gp_.subs({x: tx_, p: tp_}, simultaneous=True)
print('  粒子を f の変換 → g の変換の順に動かした結果: x = %s、p = %s' % (sp.expand(trk_x), sp.expand(trk_p)))
print('  差: x = %s、p = %s → exp(:f:)exp(:g:) の座標への作用は「f が先、g が後」' % (sp.simplify(Ux - trk_x), sp.simplify(Up - trk_p)))
Vx, _ = lie_exp(g, lie_exp(f, x)[0])
Vp, _ = lie_exp(g, lie_exp(f, p)[0])
print('  逆の積 exp(:g:) exp(:f:): x = %s、p = %s (上と違う。非可換)' % (sp.expand(Vx), sp.expand(Vp)))

print('\n== (f) [:f:, :g:] = :[f, g]: と BCH ==')
f1 = -p**2 / 2 + x * p**2
g1 = -x**3 / 3 + x**2 * p
test = [x, p, x**2 * p, x**3 + p**2]
bad = 0
for h in test:
    lhs = colon(f1, colon(g1, h)) - colon(g1, colon(f1, h))
    rhs = colon(PB(f1, g1), h)
    d = sp.simplify(lhs - rhs)
    bad += (d != 0)
print('  4 つの多項式 h で [:f:, :g:] h - :[f, g]: h の不一致の数 = %d (0 ならヤコビ恒等式が成り立つ)' % bad)

fe, ge = -eps * p**2 / 2, -eps * x**2 / 2           # f = -p^2/2 (ドリフト)、g = -x^2/2 (キック k=1)


def series_eps(expr, n):
    expr = sp.expand(expr)
    return sum(expr.coeff(eps, i) * eps**i for i in range(n + 1))


def lie_exp_eps(G, h, n):
    tot, term = 0, h
    for m in range(n + 1):
        tot += term / sp.factorial(m)
        term = series_eps(PB(G, term), n)
    return series_eps(tot, n)


Gbch2 = fe + ge + PB(fe, ge) / 2
Gbch3 = Gbch2 + (PB(fe, PB(fe, ge)) + PB(ge, PB(ge, fe))) / 12
for hname, h in (('x', x), ('p', p)):
    lhs = series_eps(lie_exp_eps(fe, lie_exp_eps(ge, h, 3), 3), 3)        # exp(:f:) exp(:g:) h
    r2 = lie_exp_eps(Gbch2, h, 3)
    r3 = lie_exp_eps(Gbch3, h, 3)
    d2 = sp.expand(lhs - r2)
    d3 = sp.expand(lhs - r3)
    print('  h = %s: exp(:f:)exp(:g:) h - exp(:f+g+[f,g]/2:) h = %s' % (hname, d2))
    print('        3 次の項 ([f,[f,g]] + [g,[g,f]])/12 を足すと差 = %s' % d3)
print('  → exp(:f:) exp(:g:) = exp(:f + g + [f, g]/2 + ([f,[f,g]] + [g,[g,f]])/12 + ...:) (BCH)。f と g が可換 ([f, g] = 0) でない限り exp(:f+g:) とは違う')

print('\n== (g) 2 次の f と線形写像 ==')
z = sp.Matrix([x, p])
f = -s * (k * x**2 + p**2) / 2
F = sp.Matrix([[sp.diff(colon(f, zi), zj) for zj in (x, p)] for zi in (x, p)])
print('  :f: (x, p) = F (x, p) の行列 F = %s' % F.tolist())
M = sp.simplify(F.exp().subs(sp.sqrt(-k * s**2), sp.I * s * sp.sqrt(k)))
Mt = sp.Matrix([[sp.cos(sp.sqrt(k) * s), sp.sin(sp.sqrt(k) * s) / sp.sqrt(k)], [-sp.sqrt(k) * sp.sin(sp.sqrt(k) * s), sp.cos(sp.sqrt(k) * s)]])
print('  行列指数 exp(F) - [[cos, sin/√k], [-√k sin, cos]] = %s' % sp.simplify((M - Mt).applyfunc(sp.simplify)))
J = sp.Matrix([[0, 1], [-1, 0]])
print('  det M = %s、M^T J M - J = %s (0 行列ならシンプレクティック)' % (sp.simplify(Mt.det()), sp.simplify(Mt.T * J * Mt - J).tolist()))
print('  2 次の f (H = (k x^2 + p^2)/2 の -s 倍) で :f: が (x, p) の 1 次式になる → exp(:f:) は行列 exp(F) になる')
