# -*- coding: utf-8 -*-
"""端部の一般論 (06_edge_fringe.md §3〜§5): ベクトルポテンシャルと、端のリー生成関数 f = ∫H₁ ds (sympy)。

(a) 偏向磁石の直線の端: A = (A_x, 0, A_s)、A_x = -(y²/2) g' + (y⁴/24) g'''、A_s = -x g (全体に 1/ρ)。∇×A が (5.31) の場に一致。
(b) 四極: A_x = -K'(x³ + 3xy²)/12、A_y = K'(3x²y + y³)/12、A_s = -K(x²-y²)/2。∇×A が (4.25) の場に一致。
(c) 2m 極 (m = 2, 3, 4): A_s = -κ Re w^m/m!、a_w = a_x + i a_y = κ' w (w^m/(m+1) - w̄^m)/(4 m!)。∇×A = ∇Φ (7.16)。
(d) 硬い端 (κ' = κ0 δ(s)) の生成関数 f = -(1/p) Re(p_w conj(a_w)) の積分 → (4.27)、(7.17)。偏向磁石は f = (y²/2ρ) p_x/p_z。
(e) 偏向磁石の y⁴ の項 (g''' = δ'') の寄与 (p_y に依存する項) と、p_y² の補正。
(f) ソレノイド: A_φ = r b/2 - r³ b''/16。b'' = b0 δ'(s) の項から f = -(b0/(8 p p_z)) p_φ p_r (ヘルプの H の paraxial 形)。
"""
import sympy as sp
import _common as C

C.start(__file__)

x, y, s = sp.symbols('x y s', real=True)
px, py, p = sp.symbols('p_x p_y p', real=True)


def curl(Ax, Ay, As):
    return (sp.diff(As, y) - sp.diff(Ay, s), sp.diff(Ax, s) - sp.diff(As, x), sp.diff(Ay, x) - sp.diff(Ax, y))


print('== (a) 偏向磁石の直線の端 (b = e B/P0 を 1/ρ で割ったもの。g(s) = B_y(s)/B_0) ==')
g = sp.Function('g', real=True)(s)
Ax = -y**2 / 2 * sp.diff(g, s) + y**4 / 24 * sp.diff(g, s, 3)
Ay = 0
As = -x * g
bx, by, bs = curl(Ax, Ay, As)
by_ref = g - y**2 / 2 * sp.diff(g, s, 2) + y**4 / 24 * sp.diff(g, s, 4)
bs_ref = y * sp.diff(g, s) - y**3 / 6 * sp.diff(g, s, 3)
print('  A_x = %s、A_s = %s' % (Ax, As))
print('  b_x = %s' % sp.simplify(bx))
print('  b_y - (g - y² g\'\'/2 + y⁴ g\'\'\'\'/24) = %s' % sp.simplify(by - by_ref))
print('  b_s - (y g\' - y³ g\'\'\'/6)           = %s' % sp.simplify(bs - bs_ref))
print('  → ∇×A = b (5.31) (y⁴・y³ まで)。A_s = -x g は本体の偏向 (一様な b_y = g)、A_x は端 (g\' ≠ 0) だけに現れる。')

print('\n== (b) 四極: K(s) ==')
K = sp.Function('K', real=True)(s)
Ax = -sp.diff(K, s) * (x**3 + 3 * x * y**2) / 12
Ay = sp.diff(K, s) * (3 * x**2 * y + y**3) / 12
As = -K * (x**2 - y**2) / 2
bx, by, bs = curl(Ax, Ay, As)
bx_ref = K * y - sp.diff(K, s, 2) / 12 * (3 * x**2 * y + y**3)
by_ref = K * x - sp.diff(K, s, 2) / 12 * (x**3 + 3 * x * y**2)
bs_ref = sp.diff(K, s) * x * y
print('  b_x - (4.25) = %s、b_y - (4.25) = %s、b_s - K\' x y = %s' % (sp.simplify(bx - bx_ref), sp.simplify(by - by_ref), sp.simplify(bs - bs_ref)))

print("\n== (c) 2m 極 (m = 2, 3, 4): Φ = κ Im w^m/m! - κ'' |w|² Im w^m/(4 (m+1)!) ==")
kk = sp.symbols('k0:6', real=True)            # κ, κ', κ'', κ''', ...


def ds(e):
    """s の全微分 (κ の導関数の記号 k0, k1, ... を通して)"""
    return sum(sp.diff(e, kk[i]) * kk[i + 1] for i in range(5))


w = x + sp.I * y
wb = x - sp.I * y
for m in (2, 3, 4):
    wm = sp.expand(w**m)
    Phi = kk[0] * sp.im(wm) / sp.factorial(m) - kk[2] * sp.expand(w * wb) * sp.im(wm) / (4 * sp.factorial(m + 1))
    a_w = kk[1] * sp.expand(w * (w**m / (m + 1) - wb**m)) / (4 * sp.factorial(m))
    Ax_ = sp.re(sp.expand(a_w))
    Ay_ = sp.im(sp.expand(a_w))
    As_ = -kk[0] * sp.re(wm) / sp.factorial(m)
    bx_ = sp.diff(As_, y) - ds(Ay_)
    by_ = ds(Ax_) - sp.diff(As_, x)
    bs_ = sp.diff(Ay_, x) - sp.diff(Ax_, y)
    res = [sp.expand(bx_ - sp.diff(Phi, x)), sp.expand(by_ - sp.diff(Phi, y)), sp.expand(bs_ - ds(Phi))]
    print('  m = %d: ∇×A - ∇Φ = %s' % (m, [sp.factor(r) for r in res]))
print("  (残りは κ''' (k3) の項だけ。b_s = ∂Φ/∂s の κ''' の項で、A の次の次数の補正が要る。κ' (k1)・κ'' (k2) の次数では 0)")

print('\n== (d) 硬い端の生成関数 f = ∫H₁ ds、H₁ = -(p_x a_x + p_y a_y)/p (p_z ≈ p)。κ\' = κ0 δ(s) ==')
k0 = sp.symbols('kappa0', real=True)
for m in (2, 3, 4):
    a_w = sp.expand(w * (w**m / (m + 1) - wb**m)) / (4 * sp.factorial(m))          # κ' の係数
    ax_, ay_ = sp.re(a_w), sp.im(a_w)
    f = sp.expand(-k0 * (px * ax_ + py * ay_) / p)
    dx_, dy_ = sp.diff(f, px), sp.diff(f, py)
    dw_formula = -k0 / p * a_w
    diff_ = sp.simplify(sp.expand(dx_ + sp.I * dy_ - dw_formula))
    print('  m = %d: Δw = ∂f/∂p_x + i ∂f/∂p_y と -(κ0/p) a_w の差 = %s' % (m, diff_))
    if m == 2:
        f27 = k0 / (12 * p) * (x**3 * px + 3 * x * y**2 * px - 3 * x**2 * y * py - y**3 * py)
        print('         (4.27) の f_in との差 = %s' % sp.simplify(f - f27))
rho = sp.symbols('rho', positive=True)
pz = sp.sqrt(p**2 - px**2 - py**2)
f_bend = sp.Rational(1, 2) * y**2 / rho * px / pz
print('  偏向磁石: a_x = -y² g\'/(2ρ) → f = -(p_x a_x/p_z) の ∫ds = (y²/2ρ) p_x/p_z。Δx = ∂f/∂p_x = %s' % sp.simplify(sp.diff(f_bend, px)))
print('     Δp_y = -∂f/∂y = %s、Δz = ∂f/∂p = %s' % (sp.simplify(-sp.diff(f_bend, y)), sp.simplify(sp.diff(f_bend, p))))

print('\n== (e) 偏向磁石の y⁴ の項: A_x の y⁴ g\'\'\'/(24ρ) (g\'\'\' = δ\'\') の寄与 ==')
sx = sp.symbols('sx', real=True)
y_s = y + py / pz * sx                      # 自由な運動での y(s) (ε = 0 次)
F = -px * y_s**4 / (24 * rho * pz)
f4 = sp.simplify(sp.diff(F, sx, 2).subs(sx, 0))        # ∫δ''(s) F(s) ds = F''(0)
print('  f₄ = F\'\'(0) = %s' % f4)
f2 = f_bend
pv = sp.sqrt(p**2 - px**2)
ser = sp.series((f2 + f4).subs(py, sp.Symbol('t') * py), sp.Symbol('t'), 0, 3).removeO().subs(sp.Symbol('t'), 1)
ref = y**2 / (2 * rho) * px / pv * (1 - py**2 / (2 * pv**2))
print('  f₂ + f₄ の p_y² までの展開 - (y²/2ρ)(p_x/p_v)(1 - p_y²/(2 p_v²)) = %s   (p_v = √(p² - p_x²))' % sp.simplify(ser - ref))
print('  → K64 の f (5.34)・(5.36) は p_y を p_z に入れない (p_v = √(p² - p_x²))。この 4 次 (y² p_y² p_x/ρ) の項は入れていない。')

print('\n== (f) ソレノイドの端: A_φ = r b/2 - r³ b\'\'/16 (b(s) = 軸上の B_z/Bρ) ==')
b = sp.Function('b', real=True)(s)
r2 = x**2 + y**2
Ax = -y * (b / 2 - r2 * sp.diff(b, s, 2) / 16)
Ay = x * (b / 2 - r2 * sp.diff(b, s, 2) / 16)
bx, by, bs = curl(Ax, Ay, 0)
bz_ref = b - r2 * sp.diff(b, s, 2) / 4
bx_ref = x * (-sp.diff(b, s) / 2 + r2 * sp.diff(b, s, 3) / 16)
by_ref = y * (-sp.diff(b, s) / 2 + r2 * sp.diff(b, s, 3) / 16)
print('  b_x - b_r x/r = %s、b_y - b_r y/r = %s、b_s - (b - r² b\'\'/4) = %s' % (
    sp.simplify(bx - bx_ref), sp.simplify(by - by_ref), sp.simplify(sp.expand(bs - bz_ref))))
print('  (b_r = -r b\'/2 + r³ b\'\'\'/16、b_z = b - r² b\'\'/4。b\'\'\'\' 以上は打ち切り。)')
pphi = x * py - y * px
# 端の項: H₁ = +p_φ r² b''/(16 p_z)。 b'' = b0 δ'(s): ∫δ' F ds = -F'(0)
b0 = sp.symbols('b0', real=True)
xs_ = x + px / pz * sx
ys_ = y + py / pz * sx
F = pphi * (xs_**2 + ys_**2) / (16 * pz)
f_sol = sp.simplify(-b0 * sp.diff(F, sx).subs(sx, 0))
pr_ = x * px + y * py
print('  f = -b0 F\'(0) = %s' % f_sol)
print('  f - (-(b0/(8 p_z²)) p_φ p_r) = %s    (p_z² → p² (paraxial) がヘルプの H = -Bz/(8 p² Bρ) p_φ p_r)' % sp.simplify(f_sol + b0 / (8 * pz**2) * pphi * pr_))
print('  (p_φ = x p_y - y p_x は自由な運動で保存するので、微分は r² だけに掛かる。)')
