# -*- coding: utf-8 -*-
r"""端部 (06_edge_fringe.md) の途中式の確認 (sympy)。

約束 (C_poisson_lie.md): exp(:f:) の流れは Δx = -∂f/∂p_x、Δp_x = +∂f/∂x、Δz = -∂f/∂p。f = -∫H_ham ds。

(a) 偏向磁石の端の磁場 (6.1): b_x = 0 として ∇·b = 0、∇×b = 0 を y のべき級数で解く (漸化式)。
(b) 2m 極の端の磁場 (6.2): κ'' の係数 c = -1/(4 (m+1)!) と、∇²_⊥(|w|² Im w^m) = 4(m+1) Im w^m。
(c) ソレノイドの端の磁場 (6.3): ∇·b = 0、∇×b = 0 が b''' (b'''') まで成り立つ。
(d) CAVI の端 (§3.4): ガウスの法則から E_r = -(r/2) E g'。スカラーポテンシャルのある H の展開 ∂H/∂(eφ) = 1/β_z。
(e) H の a_⊥ による展開 (6.6) の 1 次・2 次の係数。
(f) 1 次の摂動の生成関数 (6.7): H0 = p²/2、H1 = ε x² h(s) の 1 自由度のおもちゃで、厳密な 1 次の応答が
    「自由な運動 ∘ exp(:f:)」(f = -∫H1(自由な運動) ds) と一致すること。
(g) 偏向磁石の端のローレンツ力の積分 (6.8)(§8.3 (i) の別導出): ΔP_x、ΔP_y の力の積分から Δx、Δp_y。
(h) K64 の Maxwell フリンジの係数 a = 1 - y²/(6ρ²)、b = 1/2 - y²/(24ρ²) は sin(y/ρ)/(y/ρ)、(1 - cos(y/ρ))/(y/ρ)² の展開 (導出ではなく一致の確認)。
(i) ソレノイドの端の流れ: 極座標 (u = ln r, φ) で f = k p_φ p_u は p_φ、p_u を保存し、流れが閉じた形になる。
"""
import sympy as sp
import _common as C

C.start(__file__)

x, y, s = sp.symbols('x y s', real=True)
px, py, p = sp.symbols('p_x p_y p', real=True)

# ---------------------------------------------------------------- (a)
print('== (a) 偏向磁石の端の磁場 (6.1): b_x = 0、b_y(y, s)、b_s(y, s) を y のべき級数で ==')
g = sp.Function('g')(s)
N = 6
cy = [sp.Function('cy%d' % n)(s) for n in range(N + 1)]
cs = [sp.Function('cs%d' % n)(s) for n in range(N + 1)]
by = sum(cy[n] * y**n for n in range(N + 1))
bs = sum(cs[n] * y**n for n in range(N + 1))
# 条件: 中央面 y = 0 で b_y = g(s)、b_s = 0 (b_s は y の奇関数)。b_y は偶関数。
# ∇×b = 0 (x 成分): ∂_y b_s - ∂_s b_y = 0、∇·b = 0: ∂_y b_y + ∂_s b_s = 0
eq1 = sp.expand(sp.diff(bs, y) - sp.diff(by, s))
eq2 = sp.expand(sp.diff(by, y) + sp.diff(bs, s))
sol = {cy[0]: g}
# y^n の係数: (n+1) cs_{n+1} = cy_n'、(n+1) cy_{n+1} = -cs_n'  (n = 0, 1, 2, ...)
cs_sol = {}
cy_sol = {0: g}
cs_sol[0] = sp.Integer(0)
for n in range(0, N):
    cs_sol[n + 1] = sp.diff(cy_sol[n], s) / (n + 1)
    cy_sol[n + 1] = -sp.diff(cs_sol[n], s) / (n + 1)
print('  漸化式 (y^n の係数): (n+1) cs_{n+1} = cy_n\'、(n+1) cy_{n+1} = -cs_n\'  (cs_0 = 0、cy_0 = g)')
for n in range(0, 5):
    print('    y^%d: cy_%d = %s、 cs_%d = %s' % (n, n, cy_sol[n], n, cs_sol[n]))
ref_by = g - y**2 / 2 * sp.diff(g, s, 2) + y**4 / 24 * sp.diff(g, s, 4)
ref_bs = y * sp.diff(g, s) - y**3 / 6 * sp.diff(g, s, 3)
by_s = sum(cy_sol[n] * y**n for n in range(5))
bs_s = sum(cs_sol[n] * y**n for n in range(5))
print('  b_y - (6.1) = %s、b_s - (6.1) = %s' % (sp.simplify(by_s - ref_by), sp.simplify(bs_s - ref_bs)))

# ---------------------------------------------------------------- (b)
print('\n== (b) 2m 極: Φ = κ Im w^m/m! + c κ\'\' |w|² Im w^m、∇²Φ = 0 から c ==')
kk = sp.symbols('k0:5', real=True)       # κ, κ', κ'', κ''', κ''''


def ds(e):
    return sum(sp.diff(e, kk[i]) * kk[i + 1] for i in range(4))


w = x + sp.I * y
wb = x - sp.I * y
cc = sp.symbols('c')
for m in (1, 2, 3, 4):
    wm = sp.expand(w**m)
    Imwm = sp.im(wm)
    lap_perp = lambda F: sp.diff(F, x, 2) + sp.diff(F, y, 2)
    ident = sp.simplify(lap_perp(sp.expand(w * wb) * Imwm) - 4 * (m + 1) * Imwm)
    Phi = kk[0] * Imwm / sp.factorial(m) + cc * kk[2] * sp.expand(w * wb) * Imwm
    # ∇²Φ = ∂_s²Φ + ∇⊥²Φ の κ'' の項
    lap = ds(ds(Phi)) + lap_perp(Phi)
    coef = sp.simplify(sp.expand(lap).coeff(kk[2]))
    csol = sp.solve(coef, cc)
    print('  m = %d: ∇⊥²(|w|² Im w^m) - 4(m+1) Im w^m = %s;  κ\'\' の項を消す c = %s  (-1/(4 (m+1)!) = %s)' % (
        m, ident, csol, -sp.Rational(1, 4) / sp.factorial(m + 1)))
print('  (∇⊥² = 4 ∂_w ∂_w̄。|w|² Im w^m = (w^{m+1} w̄ - w w̄^{m+1})/(2i)、∂_w ∂_w̄ = ((m+1) w^m - (m+1) w̄^m)/(2i)。)')

# ---------------------------------------------------------------- (c)
print('\n== (c) ソレノイド: b_z = b - r² b\'\'/4、b_r = -r b\'/2 + r³ b\'\'\'/16 (r = √(x²+y²)) ==')
r = sp.symbols('r', positive=True)
b = sp.Function('b')(s)
bz = b - r**2 * sp.diff(b, s, 2) / 4 + r**4 * sp.diff(b, s, 4) / 64
br = -r * sp.diff(b, s) / 2 + r**3 * sp.diff(b, s, 3) / 16 - r**5 * sp.diff(b, s, 5) / 384
div = sp.simplify(sp.diff(r * br, r) / r + sp.diff(bz, s))
curl = sp.simplify(sp.diff(bz, r) - sp.diff(br, s))
print('  ∇·b = (1/r)∂_r(r b_r) + ∂_s b_z = %s、(∇×b)_φ = ∂_r b_z - ∂_s b_r = %s  (b⁽⁶⁾ 以上を無視)' % (div, curl))
A_phi = r * b / 2 - r**3 * sp.diff(b, s, 2) / 16
print('  ベクトルポテンシャル A_φ = r b/2 - r³ b\'\'/16: (1/r)∂_r(r A_φ) = %s、-∂_s A_φ = %s' % (
    sp.simplify(sp.diff(r * A_phi, r) / r), sp.simplify(-sp.diff(A_phi, s))))

# ---------------------------------------------------------------- (d)
print('\n== (d) CAVI の端 ==')
E = sp.Function('E')(s)
Ez = sp.Function('Ez')(s)
Er = -r / 2 * sp.diff(Ez, s)
print('  ∇·E = (1/r)∂_r(r E_r) + ∂_s E_z (E_r = -(r/2) E_z\') = %s' % sp.simplify(sp.diff(r * Er, r) / r + sp.diff(Ez, s)))
phi, e, Etot, m_, c_ = sp.symbols('phi e E m c', positive=True)
pz2 = (Etot - e * phi)**2 - m_**2 - px**2 - py**2        # c = 1
Hphi = -sp.sqrt(pz2)
dH = sp.simplify(sp.diff(Hphi, phi).subs(phi, 0))
print('  H = -√((E - eφ)² - m² - p⊥²) (c = 1): ∂H/∂φ|_(φ=0) = %s = e E/p_z = e/β_z' % dH)
print('   → 運動量のキック ΔP = -∫∂H/∂x ds = -(e/β_z) ∫∂_xφ ds: 速度で割る因子 1/β_z が H から自動的に出る (H_f の係数に 1/β は書かれていない)。')

# ---------------------------------------------------------------- (e)
print('\n== (e) H = -√(p² - (p_x - a_x)² - (p_y - a_y)²) の a_⊥ による展開 ==')
ax, ay = sp.symbols('a_x a_y', real=True)
eps = sp.symbols('epsilon')
H = -sp.sqrt(p**2 - (px - eps * ax)**2 - (py - eps * ay)**2)
pz = sp.sqrt(p**2 - px**2 - py**2)
H1 = sp.simplify(sp.diff(H, eps).subs(eps, 0))
print('  ∂H/∂ε|_0 = %s    (H₁ = -(p_x a_x + p_y a_y)/p_z との差 = %s)' % (H1, sp.simplify(H1 + (px * ax + py * ay) / pz)))
H2 = sp.simplify(sp.diff(H, eps, 2).subs(eps, 0) / 2)
H2_ref = (ax**2 + ay**2) / (2 * pz) + (px * ax + py * ay)**2 / (2 * pz**3)
print('  (1/2) ∂²H/∂ε²|_0 = %s    (+a²/(2p_z) + (p·a)²/(2p_z³) との差 = %s)' % (H2, sp.simplify(H2 - H2_ref)))

# ---------------------------------------------------------------- (f)
print('\n== (f) 1 次の摂動: H0 = p²/2、H1 = ε x² h(s) (h(s) = s + s²、区間 0〜S) ==')
S = sp.symbols('S', positive=True)
q0, p0 = sp.symbols('q0 p0', real=True)
h = s + s**2
t = sp.symbols('t', real=True)
x0 = q0 + p0 * t
# 厳密な 1 次の応答: x' = p、p' = -2 ε x h  → x = x0 + ε x1、x1'' = -2 x0 h、x1(0) = x1'(0) = 0
x1S = sp.integrate((S - t) * (-2 * x0 * h.subs(s, t)), (t, 0, S))
x1pS = sp.integrate(-2 * x0 * h.subs(s, t), (t, 0, S))
# exp(:f:) の f = -∫ H1(自由な運動) ds = -ε ∫ (q0 + p0 t)² h dt。出力の写像 = 自由な運動 (S) ∘ exp(:f:)
f = -sp.integrate((q0 + p0 * t)**2 * h.subs(s, t), (t, 0, S))     # ε の係数
dq = -sp.diff(f, p0)
dp = sp.diff(f, q0)
x_out = (q0 + dq) + S * (p0 + dp)
p_out = p0 + dp
print('  x_out - x_exact(1 次) = %s、p_out - p_exact(1 次) = %s' % (
    sp.simplify(sp.expand(x_out - (q0 + p0 * S) - x1S)), sp.simplify(sp.expand(p_out - p0 - x1pS))))
print('  (ε の 1 次。出力 = 自由な運動 ∘ exp(:f:)、f = -∫H1(自由な運動) ds。階段の端が s = 0 にあれば自由な運動は恒等で f = -H1 の積分そのもの。)')

# ---------------------------------------------------------------- (g)
print('\n== (g) 偏向磁石の端のローレンツ力の積分 (p_y = 0 から始める。b_y = ε(g - y² g\'\'/2)、b_s = ε y g\'、ε = 1/ρ_b) ==')
rho = sp.symbols('rho', positive=True)
Px, Py, Pz = sp.symbols('P_x P_y P_z', real=True)
gp = sp.symbols('gp', real=True)         # g'(s)
# dP/ds = P × b / P_z (b = (0, b_y, b_s)):  (P × b)_x = P_y b_s - P_z b_y、(P × b)_y = P_z b_x - P_x b_s = -P_x b_s
by_f = (1 / rho) * (-y**2 / 2) * sp.Symbol('gpp')      # 端の項だけ (本体 (1/ρ) g は本体の偏向)
bs_f = (1 / rho) * y * gp
dPx = (Py * bs_f - Pz * by_f) / Pz
dPy = (-Px * bs_f) / Pz
print('  端の項: dP_x/ds = %s (P_y = 0)、dP_y/ds = %s' % (sp.simplify(dPx.subs(Py, 0)), sp.simplify(dPy)))
print('   ∫ g\' ds = 1 (0 → 1 の立ち上がり)、∫ g\'\' ds = 0。P_x の端の項の積分: ∫ΔP_x(s) ds, ΔP_x(s) = +(y²/(2ρ)) g\'(s) (g\'\' の積分から)')
print('   → ΔP_y = -(P_x/P_z)(y/ρ)、∫ΔP_x ds = y²/(2ρ)。dx/ds = P_x/P_z の摂動: Δx = ∂(P_x/P_z)/∂P_x ∫ΔP_x ds')
Pxs = sp.symbols('P_xs', real=True)
dxdPx = sp.simplify(sp.diff(Pxs / sp.sqrt(p**2 - Pxs**2 - py**2), Pxs))
print('     ∂(P_x/P_z)/∂P_x = %s = (p² - p_y²)/P_z³' % dxdPx)
print('     Δx = (y²/(2ρ)) (p² - p_y²)/p_z³ = -∂f/∂p_x (f = -(y²/2ρ) p_x/p_z) と一致 (edge_01 (d))')

# ---------------------------------------------------------------- (h)
print('\n== (h) K64 の Maxwell フリンジの係数 a、b と三角関数の展開 ==')
th = sp.symbols('theta')
print('  sin θ/θ = %s' % sp.series(sp.sin(th) / th, th, 0, 6))
print('  (1 - cos θ)/θ² = %s' % sp.series((1 - sp.cos(th)) / th**2, th, 0, 6))
print('  K64: a = 1 - (y/ρ)²/6、b = 1/2 - (y/ρ)²/24 (θ = y/ρ)。θ² までの係数が一致 (θ⁴ 以上は K64 に無い)。')
print('  すなわち Δp_y = -(p_x/p_v) sin(y/ρ)、Δx = ρ(1 - cos(y/ρ)) p²/p_v³ の θ² までの展開。物理的な導出は未 (本文 §8.2)。')

# ---------------------------------------------------------------- (i)
print('\n== (i) ソレノイドの端: f = k p_φ p_u (u = ln r、p_u = x p_x + y p_y、k = b0/(8p²)) の流れ ==')
r_, ph, pu, pph, kk_, pp = sp.symbols('r phi p_u p_phi k p', real=True)
# (x, y, p_x, p_y) → (u, φ, p_u, p_φ) は正準変換: p_u = x p_x + y p_y、p_φ = x p_y - y p_x
xx, yy, pxx, pyy = sp.symbols('xx yy pxx pyy', real=True)
u = sp.log(sp.sqrt(xx**2 + yy**2))
phi_ = sp.atan2(yy, xx)
pu_ = xx * pxx + yy * pyy
pph_ = xx * pyy - yy * pxx
pb = lambda A, B: (sp.diff(A, xx) * sp.diff(B, pxx) + sp.diff(A, yy) * sp.diff(B, pyy)
                   - sp.diff(A, pxx) * sp.diff(B, xx) - sp.diff(A, pyy) * sp.diff(B, yy))
print('  [u, p_u] = %s、[φ, p_φ] = %s、[u, p_φ] = %s、[φ, p_u] = %s、[u, φ] = %s、[p_u, p_φ] = %s' % tuple(
    sp.simplify(pb(A, B)) for A, B in ((u, pu_), (phi_, pph_), (u, pph_), (phi_, pu_), (u, phi_), (pu_, pph_))))
print('  f = k p_φ p_u は (u, φ) を含まないので p_u、p_φ は保存。Δu = -∂f/∂p_u = -k p_φ、Δφ = -∂f/∂p_φ = -k p_u、Δz = -∂f/∂p = +2k p_φ p_u/p (k ∝ 1/p²)。')
print('  K64: bq = -dbz/(8p) (q 単位: p_φ/p を使うので k p_φ = bq p_φ,q): r1 = r e^{bq p_φ}、φ1 = φ + bq p_r、z1 = z + bp p_r p_φ (bp = -2 bq)。')
