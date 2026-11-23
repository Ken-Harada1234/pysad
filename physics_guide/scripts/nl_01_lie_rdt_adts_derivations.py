# -*- coding: utf-8 -*-
r"""非線形への橋渡し (18_nonlinear_bridge.md) の式の導出を sympy で確かめる。

規約: ポアソン括弧 {f, g} = Σ (∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q) (C_poisson_lie.md の (C.1))、:f: g = {f, g}、
要素の写像は exp(-s:H:) (z∘Φ = exp(:F:) z、F = -s H)。複素変数 ζ = x^ - i p^ (ζ^+)、ζbar = x^ + i p^ (ζ^-)。
(a) {ζ, ζbar} = 2i と {ζ^a ζbar^b, ζ^c ζbar^d} = 2i (ad - bc) ζ^(a+c-1) ζbar^(b+d-1)。(φ, J) が正準対で {F, J} = ∂F/∂φ。
(b) 6 極・8 極・運動学・端部フリンジのハミルトニアンと kick。
(c) 正規化座標: x = (w ζ + wbar ζbar)/2、p_x = (w' ζ + wbar' ζbar)/2。
(d) 薄い 6 極の h_jklmp (分散の項を含む) と 1 次の ΔJ = -∂h/∂φ。
(e) 1 次の ADTS: 8 極、運動学的な項、端部フリンジ。4 次の h の対角から ∂ν/∂J = 4h/π など。
(f) 写像の合成と Magnus 展開 (2 つの 3 次の項の合成を、明示の写像で確かめる)、BCH (3 つの指数の積)。
(g) Dragt-Finn: {f,{f,z}} = (u.∇)u。2 次のテイラー係数から f_3 を復元。
(h) 正規形: F4' = F4 + (1/2){g, Rg}、1 次元の対角の係数 (cot の式)。
(i) 薄い四極のチューンシフトと、6 極のクロマティシティ ξ = h11001/π。
(j) 共鳴の近く (1 本の共鳴): 固定点と島の幅 (振り子)。RF バケットの分離線。
"""
import itertools
import sympy as sp
import _common as C

C.start(__file__)

I = sp.I
zt, zb = sp.symbols('zt zb')            # ζ, ζbar (独立な記号として扱う)
vt, vb = sp.symbols('vt vb')            # y の複素変数
xh, ph, yh, qh = sp.symbols('xh ph yh qh', real=True)   # x^, p^ (正規化座標)


def pb_zeta(f, g):
    """複素変数での括弧 {f, g}: {ζ, ζbar} = 2i、{v, vbar} = 2i。"""
    return sp.expand(2 * I * (sp.diff(f, zt) * sp.diff(g, zb) - sp.diff(f, zb) * sp.diff(g, zt)
                              + sp.diff(f, vt) * sp.diff(g, vb) - sp.diff(f, vb) * sp.diff(g, vt)))


import sys
print('== (a) 複素変数のポアソン括弧 ==')
x_, p_ = xh, ph
zeta = x_ - I * p_
zbar = x_ + I * p_
print('  {ζ, ζbar} (x^, p^ で {x^, p^} = 1) = %s' % sp.simplify(sp.diff(zeta, x_) * sp.diff(zbar, p_) - sp.diff(zeta, p_) * sp.diff(zbar, x_)))


def pb_xp(f, g):
    return sp.expand(sp.diff(f, xh) * sp.diff(g, ph) - sp.diff(f, ph) * sp.diff(g, xh))


bad = 0
for a, b, c, d in itertools.product(range(0, 4), repeat=4):
    lhs = pb_xp(zeta**a * zbar**b, zeta**c * zbar**d)
    rhs = 2 * I * (a * d - b * c) * zeta**(a + c - 1) * zbar**(b + d - 1) if (a + c >= 1 and b + d >= 1) else 0
    if sp.expand(lhs - rhs) != 0:
        bad += 1
print('  {ζ^a ζbar^b, ζ^c ζbar^d} = 2i (ad - bc) ζ^(a+c-1) ζbar^(b+d-1): a, b, c, d = 0..3 の 256 通りで不一致 %d 個' % bad)
# (φ, J) が正準対: ζ = sqrt(2J) e^{iφ}
J, phi = sp.symbols('J phi', positive=True)
zJ = sp.sqrt(2 * J) * sp.exp(I * phi)
zbJ = sp.sqrt(2 * J) * sp.exp(-I * phi)
pb_phiJ = lambda f, g: sp.simplify(sp.diff(f, phi) * sp.diff(g, J) - sp.diff(f, J) * sp.diff(g, phi))
print('  {ζ, ζbar}_(φ,J) = %s (括弧を {f, g} = f_φ g_J - f_J g_φ としたとき。(φ, J) を (q, p) とみなす)' % pb_phiJ(zJ, zbJ))
# x^ = sqrt(2J) cos φ, p^ = -sqrt(2J) sin φ から {x^, p^}_(φ,J) = 1 を確かめる
print('  {x^, p^}_(φ,J) = %s' % pb_phiJ(sp.sqrt(2 * J) * sp.cos(phi), -sp.sqrt(2 * J) * sp.sin(phi)))
F = sp.Function('F')(phi, J)
print('  {F, J}_(φ,J) = %s = ∂F/∂φ、{F, φ}_(φ,J) = %s = -∂F/∂J' % (pb_phiJ(F, J), pb_phiJ(F, phi)))

print('\n== (b) 要素のハミルトニアン ==')
x, y, K2, K3, px, py = sp.symbols('x y K2 K3 px py', real=True)
H3 = sp.expand(K2 / 6 * sp.re(sp.expand((x + I * y)**3)))
H3 = sp.expand(K2 / 6 * (x**3 - 3 * x * y**2))
chk = sp.expand(sp.re(sp.expand((x + I * y)**3, complex=True)) - (x**3 - 3 * x * y**2))
print('  Re (x + iy)^3 - (x^3 - 3 x y^2) = %s' % chk)
print('  6 極: H3 = (K2/6)(x^3 - 3xy^2) → p_x\' = -∂H/∂x = %s、p_y\' = -∂H/∂y = %s' % (sp.factor(-sp.diff(H3, x)), sp.factor(-sp.diff(H3, y))))
H4 = sp.expand(K3 / 24 * (x**4 - 6 * x**2 * y**2 + y**4))
chk4 = sp.expand(sp.re(sp.expand((x + I * y)**4, complex=True)) - (x**4 - 6 * x**2 * y**2 + y**4))
print('  Re (x + iy)^4 - (x^4 - 6 x^2 y^2 + y^4) = %s。8 極: p_x\' = %s' % (chk4, sp.factor(-sp.diff(H4, x))))
pp = sp.symbols('pp', positive=True)
print('  運動学: -sqrt(1 - p^2) の展開 = %s (p^2 = p_x^2 + p_y^2)' % sp.series(-sp.sqrt(1 - pp**2), pp, 0, 7).removeO())
Kq = sp.symbols('K1', real=True)

print('\n== (c) 正規化座標 ==')
bx, ax, mu = sp.symbols('beta alpha mu', real=True)       # β, α, μ (位相進み)
# x^ = x/√β、p^ = (αx + β p_x)/√β → x = √β x^、p_x = (p^ - α x^)/√β
xs = sp.sqrt(bx) * xh
pxs = (ph - ax * xh) / sp.sqrt(bx)
# 位相進み μ だけ回転した ζ: ζ_local = e^{iμ} ζ (x^ = Re ζ_local)
w = sp.sqrt(bx) * sp.exp(I * mu)
wb = sp.sqrt(bx) * sp.exp(-I * mu)
wp = (-ax + I) / sp.sqrt(bx) * sp.exp(I * mu)       # w' = dw/ds (β' = -2α、μ' = 1/β)
wpb = (-ax - I) / sp.sqrt(bx) * sp.exp(-I * mu)
# 検算: dw/ds
s_ = sp.symbols('s_', real=True)
bf, muf = sp.Function('b')(s_), sp.Function('m')(s_)
wf = sp.sqrt(bf) * sp.exp(I * muf)
dw = sp.diff(wf, s_).subs({sp.Derivative(bf, s_): -2 * ax, sp.Derivative(muf, s_): 1 / bx}).subs({bf: bx, muf: mu})
print('  dw/ds - w\' = %s (β\' = -2α、μ\' = 1/β を使った)' % sp.simplify(dw - wp))
# 位置 s の局所の正規化座標 (x^_loc, p^_loc) = 起点の (x^, p^) を μ だけ回転したもの (ζ_loc = e^{iμ} ζ)
xl = xh * sp.cos(mu) + ph * sp.sin(mu)
pl = -xh * sp.sin(mu) + ph * sp.cos(mu)
x_from_w = (w * (xh - I * ph) + wb * (xh + I * ph)) / 2
px_from_w = (wp * (xh - I * ph) + wpb * (xh + I * ph)) / 2
print('  (w ζ + wbar ζbar)/2 - √β x^_loc = %s' % sp.simplify(sp.expand(x_from_w - sp.sqrt(bx) * xl, complex=True)))
print("  (w' ζ + wbar' ζbar)/2 - (p^_loc - α x^_loc)/√β = %s" % sp.simplify(sp.expand(px_from_w - (pl - ax * xl) / sp.sqrt(bx), complex=True)))
print('  |w\'|^2 = %s (= γ = (1 + α^2)/β)' % sp.simplify(wp * wpb))

print('\n== (d) 薄い 6 極の h_jklmp ==')
wx, wxb, wy, wyb, ex, exb = sp.symbols('wx wxb wy wyb ex exb')
dl = sp.symbols('delta')
bxs, bys, etas, K2L = sp.symbols('beta_x beta_y eta K2L', positive=True)
X = (wx * zt + wxb * zb) / 2
Y = (wy * vt + wyb * vb) / 2
Hthin = sp.expand(K2L / 6 * ((X + etas * dl)**3 - 3 * (X + etas * dl) * Y**2))
P = sp.Poly(Hthin, zt, zb, vt, vb, dl)
subsw = {wx: sp.sqrt(bxs) * sp.exp(I * sp.Symbol('mux')), wxb: sp.sqrt(bxs) * sp.exp(-I * sp.Symbol('mux')),
         wy: sp.sqrt(bys) * sp.exp(I * sp.Symbol('muy')), wyb: sp.sqrt(bys) * sp.exp(-I * sp.Symbol('muy'))}
for key in [(3, 0, 0, 0, 0), (2, 1, 0, 0, 0), (1, 0, 1, 1, 0), (1, 0, 0, 2, 0), (1, 0, 2, 0, 0), (1, 1, 0, 0, 1), (0, 0, 1, 1, 1), (2, 0, 0, 0, 1), (0, 0, 2, 0, 1), (1, 0, 0, 0, 2)]:
    c = P.coeff_monomial(key)
    print('  h%d%d%d%d%d = %s' % (*key, sp.simplify(c.subs(subsw))))
print('  (参考: docs の式 h30000 = K2L β_x^{3/2} e^{3iμx}/48、h21000 = K2L β_x^{3/2} e^{iμx}/16、h10110 = -K2L √β_x β_y e^{iμx}/8、h10020 = -K2L √β_x β_y e^{i(μx-2μy)}/16、h11001 = K2L η β_x/4、h00111 = -K2L η β_y/4)')
# 1 次の ΔJ = -∂h/∂φ (1 つの薄い 6 極、y なし、δ = 0) と、実際の kick の J の変化
a_ = sp.symbols('a_', positive=True)    # a = K2L β^{3/2}
Hn = a_ / 6 * xh**3                     # 正規化座標での H3 (x = √β x^)
dJ_kick = sp.expand(ph * (-a_ / 2 * xh**2))   # 実際の kick Δp^ = -(a/2) x^2 → ΔJ = p^ Δp^ + (Δp^)^2/2 の 1 次
xs_, ps_ = sp.sqrt(2 * J) * sp.cos(phi), -sp.sqrt(2 * J) * sp.sin(phi)
lhs = sp.simplify(dJ_kick.subs({xh: xs_, ph: ps_}))
rhs = sp.simplify(-sp.diff(Hn.subs(xh, xs_), phi))
print('  1 次の ΔJ: kick の p^ Δp^ = %s、-∂h/∂φ = %s、差 %s' % (sp.simplify(lhs), sp.simplify(rhs), sp.simplify(lhs - rhs)))

print('\n== (e) 1 次の ADTS: ∂ν/∂J = 4h22000/π、2 h11110/π、4 h00220/π ==')
K3L = sp.symbols('K3L', positive=True)
H4thin = K3L / 24 * (X**4 - 6 * X**2 * Y**2 + Y**4)
P4 = sp.Poly(sp.expand(H4thin), zt, zb, vt, vb)
b22 = sp.simplify(P4.coeff_monomial((2, 2, 0, 0)).subs(subsw))
b11 = sp.simplify(P4.coeff_monomial((1, 1, 1, 1)).subs(subsw))
b00 = sp.simplify(P4.coeff_monomial((0, 0, 2, 2)).subs(subsw))
print('  8 極: h22000 = %s、h11110 = %s、h00220 = %s' % (b22, b11, b00))
print('  ∂νx/∂Jx = 4h22000/π = %s、∂νx/∂Jy = 2h11110/π = %s、∂νy/∂Jy = 4h00220/π = %s' % (sp.simplify(4 * b22 / sp.pi), sp.simplify(2 * b11 / sp.pi), sp.simplify(4 * b00 / sp.pi)))
# 運動学的な項 (px^2 + py^2)^2 / 8、p = Re(w' ζ)
gx, gy = sp.symbols('gamma_x gamma_y', positive=True)
mx, my = sp.symbols('mx my', real=True)
Pkx = (sp.sqrt(gx) * sp.exp(I * mx) * zt + sp.sqrt(gx) * sp.exp(-I * mx) * zb) / 2   # |w'|^2 = γ
Pky = (sp.sqrt(gy) * sp.exp(I * my) * vt + sp.sqrt(gy) * sp.exp(-I * my) * vb) / 2
Hk = sp.expand((Pkx**2 + Pky**2)**2 / 8)
Pk = sp.Poly(Hk, zt, zb, vt, vb)
k22, k11, k00 = [sp.simplify(Pk.coeff_monomial(k)) for k in ((2, 2, 0, 0), (1, 1, 1, 1), (0, 0, 2, 2))]
print('  運動学 (px^2 + py^2)^2/8: h22000 = %s、h11110 = %s、h00220 = %s' % (k22, k11, k00))
print('   → ∂νx/∂Jx = %s ∮ds、∂νx/∂Jy = %s ∮ds、∂νy/∂Jy = %s ∮ds (docs: 3γx²/(16π)、γxγy/(8π)、3γy²/(16π))' % (sp.simplify(4 * k22 / sp.pi), sp.simplify(2 * k11 / sp.pi), sp.simplify(4 * k00 / sp.pi)))
# 端部フリンジ
sg, K1, axx, ayy, bxx, byy = sp.symbols('sigma K1 ax ay bx by', real=True)
bxp, byp = sp.symbols('bxp byp', positive=True)
Xf = (sp.sqrt(bxp) * sp.exp(I * mx) * zt + sp.sqrt(bxp) * sp.exp(-I * mx) * zb) / 2
Yf = (sp.sqrt(byp) * sp.exp(I * my) * vt + sp.sqrt(byp) * sp.exp(-I * my) * vb) / 2
# p_x = (w' ζ + c.c.)/2、w' = (-α + i)/√β e^{iμ}
Pxf = ((-axx + I) / sp.sqrt(bxp) * sp.exp(I * mx) * zt + (-axx - I) / sp.sqrt(bxp) * sp.exp(-I * mx) * zb) / 2
Pyf = ((-ayy + I) / sp.sqrt(byp) * sp.exp(I * my) * vt + (-ayy - I) / sp.sqrt(byp) * sp.exp(-I * my) * vb) / 2
Hf = sp.expand(sg * K1 / 12 * ((Xf**3 + 3 * Xf * Yf**2) * Pxf - (Yf**3 + 3 * Xf**2 * Yf) * Pyf))
Pf = sp.Poly(Hf, zt, zb, vt, vb)
f22, f11, f00 = [sp.simplify(Pf.coeff_monomial(k)) for k in ((2, 2, 0, 0), (1, 1, 1, 1), (0, 0, 2, 2))]
print('  端のフリンジ H_f: ∂νx/∂Jx = %s、∂νx/∂Jy = %s、∂νy/∂Jy = %s' % (sp.simplify(4 * f22 / sp.pi), sp.simplify(2 * f11 / sp.pi), sp.simplify(4 * f00 / sp.pi)))
print('   (docs: -σK1 αx βx/(8π)、σK1 (αy βx - αx βy)/(8π)、σK1 αy βy/(8π)。ここで β は √β^2 = β の形)')

print('\n== (f) 写像の合成と Magnus 展開 ==')
q, pp_ = sp.symbols('q pp_')
c1, c2 = sp.symbols('c1 c2')


def pbqp(f, g):
    return sp.expand(sp.diff(f, q) * sp.diff(g, pp_) - sp.diff(f, pp_) * sp.diff(g, q))


def trunc(expr, deg):
    """q, pp_ の次数 deg 以下の部分"""
    P_ = sp.Poly(sp.expand(expr), q, pp_)
    return sum(c * q**m[0] * pp_**m[1] for m, c in P_.terms() if m[0] + m[1] <= deg)


def lie_apply(f, g, nmax=6, deg=3):
    """exp(:f:) g。f は 3 次以上の同次式の和なので、括弧のたびに次数が上がる。deg 次より上は捨てる"""
    tot, term = g, g
    for n in range(1, nmax + 1):
        term = trunc(pbqp(f, term), deg) if deg else pbqp(f, term)
        if term == 0:
            break
        tot += term / sp.factorial(n)
    return sp.expand(tot)


H1 = c1 * q**3
H2 = c2 * pp_**3
# 点の写像: Φ1 (早い、H1 = c1 q^3): z∘Φ1 = exp(-:H1:) z → (q, p - 3 c1 q^2)。Φ2: (q + 3 c2 p^2, p)
print('  z∘Φ1 = exp(-:H1:) z: q → %s、p → %s' % (lie_apply(-H1, q), lie_apply(-H1, pp_)))
q1, p1 = lie_apply(-H1, q), lie_apply(-H1, pp_)
q2 = q1 + 3 * c2 * p1**2
p2 = p1
Fm = -H1 - H2 + sp.Rational(1, 2) * pbqp(H1, H2)          # H1 が早い: F = -(H1 + H2) + (1/2){H1, H2}
tq, tp = trunc(q2, 3), trunc(p2, 3)
lq, lp = trunc(lie_apply(Fm, q, 4), 3), trunc(lie_apply(Fm, pp_, 4), 3)
print('  合成の写像 (Φ2∘Φ1、3 次まで) と exp(:F:) z、F = -(H1 + H2) + (1/2){H1, H2} (H1 が早い) の差: q %s、p %s' % (sp.expand(tq - lq), sp.expand(tp - lp)))
Fm_wrong = -H1 - H2 - sp.Rational(1, 2) * pbqp(H1, H2)
lq_w, lp_w = trunc(lie_apply(Fm_wrong, q, 4), 3), trunc(lie_apply(Fm_wrong, pp_, 4), 3)
print('  対照: 括弧の符号を逆にした F の差: q %s、p %s (0 でない = 順序・符号の約束が上のとおりに決まる)' % (sp.expand(tq - lq_w), sp.expand(tp - lp_w)))
# [:a:, :b:] = :{a, b}: (ヤコビ恒等式)
aa = c1 * q**3 + c2 * q * pp_**2
bb = q**2 * pp_ + 2 * pp_**3
tf = q**2 * pp_ + q * pp_**3
comm = pbqp(aa, pbqp(bb, tf)) - pbqp(bb, pbqp(aa, tf))
print('  [:a:, :b:] g - :{a, b}: g = %s (ヤコビ恒等式)' % sp.expand(comm - pbqp(pbqp(aa, bb), tf)))
# 3 つの指数の積の BCH (3 次の項どうしの括弧まで): exp(:a:) exp(:b:) exp(:c:) = exp(:a + b + c + (1/2)({a,b} + {a,c} + {b,c}):)
ca = sp.symbols('ca0:4'); cb = sp.symbols('cb0:4'); cc = sp.symbols('cc0:4')
cub = lambda cs: cs[0] * q**3 + cs[1] * q**2 * pp_ + cs[2] * q * pp_**2 + cs[3] * pp_**3
a3, b3, c3 = cub(ca), cub(cb), cub(cc)
lhs_q = lie_apply(a3, lie_apply(b3, lie_apply(c3, q)))          # exp(:a:) exp(:b:) exp(:c:) q (c が最初に作用)
lhs_p = lie_apply(a3, lie_apply(b3, lie_apply(c3, pp_)))
dd = a3 + b3 + c3 + sp.Rational(1, 2) * (pbqp(a3, b3) + pbqp(a3, c3) + pbqp(b3, c3))
rhs_q, rhs_p = lie_apply(dd, q), lie_apply(dd, pp_)
print('  BCH (a, b, c は一般の 3 次式): 3 次までの差 q %s、p %s' % (sp.expand(trunc(lhs_q - rhs_q, 3)), sp.expand(trunc(lhs_p - rhs_p, 3))))

print('\n== (g) Dragt-Finn: {f,{f,z}} = (u.∇)u、2 次のテイラー係数から f_3 ==')
cs = sp.symbols('d0:4')
f3 = cub(cs)
u = (pbqp(f3, q), pbqp(f3, pp_))        # u_i = {f, z_i}
uq, up = u
ff_q = pbqp(f3, pbqp(f3, q))
ff_p = pbqp(f3, pbqp(f3, pp_))
udu_q = uq * sp.diff(uq, q) + up * sp.diff(uq, pp_)
udu_p = uq * sp.diff(up, q) + up * sp.diff(up, pp_)
print('  {f,{f,q}} - (u.∇)u_q = %s、{f,{f,p}} - (u.∇)u_p = %s' % (sp.expand(ff_q - udu_q), sp.expand(ff_p - udu_p)))
print('  u = ({f,q}, {f,p}) = (-∂f/∂p, ∂f/∂q): 差 %s、%s' % (sp.expand(uq + sp.diff(f3, pp_)), sp.expand(up - sp.diff(f3, q))))
print('  f3 の復元 (オイラー): (1/3)(q ∂f/∂q + p ∂f/∂p) - f3 = %s' % sp.expand((q * sp.diff(f3, q) + pp_ * sp.diff(f3, pp_)) / 3 - f3))
# 6 極の kick Δp = -(K2L/2) q^2 → f3 = -K2L q^3/6 = -H3
K2s = sp.symbols('K2s')
print('  6 極の kick Δp = -(K2L/2) q^2 から: ∂f/∂q = %s → f = %s (= -H3)' % (-K2s / 2 * q**2, sp.integrate(-K2s / 2 * q**2, q)))

print('\n== (h) 正規形: F4\' = F4 + (1/2){g, Rg} と 1 次元の対角の係数 ==')
mu_ = sp.symbols('mu_', real=True)
fj = {(3, 0): sp.Symbol('F30'), (0, 3): sp.Symbol('F03'), (2, 1): sp.Symbol('F21'), (1, 2): sp.Symbol('F12')}


def pbz(f, g):
    return sp.expand(2 * I * (sp.diff(f, zt) * sp.diff(g, zb) - sp.diff(f, zb) * sp.diff(g, zt)))


F3 = sum(c * zt**j * zb**k for (j, k), c in fj.items())
g3 = sum(c / (sp.exp(I * (j - k) * mu_) - 1) * zt**j * zb**k for (j, k), c in fj.items())      # g_jk = F_jk/(e^{i(j-k)μ} - 1)
Rg3 = sum(c / (sp.exp(I * (j - k) * mu_) - 1) * sp.exp(I * (j - k) * mu_) * zt**j * zb**k for (j, k), c in fj.items())
print('  F3 + g - Rg = %s (非共鳴の 3 次の項が消える)' % sp.simplify(F3 + g3 - Rg3))
print('  代数: (1/2){g + Rg, F3} - (1/2){g, Rg} - (1/2){g, Rg} (F3 = Rg - g を使って) = %s' % sp.simplify(sp.Rational(1, 2) * pbz(g3 + Rg3, Rg3 - g3) - sp.Rational(1, 2) * pbz(g3, Rg3) - sp.Rational(1, 2) * pbz(g3, Rg3)))
diag = sp.expand(sp.Rational(1, 2) * pbz(g3, Rg3))
cd = sp.Poly(diag, zt, zb).coeff_monomial((2, 2))
cd = cd.subs({fj[(0, 3)]: sp.conjugate(fj[(3, 0)]), fj[(1, 2)]: sp.conjugate(fj[(2, 1)])})
cot = lambda t: sp.cos(t) / sp.sin(t)
expect = 9 * sp.Abs(fj[(3, 0)])**2 * cot(3 * mu_ / 2) + 3 * sp.Abs(fj[(2, 1)])**2 * cot(mu_ / 2)
num = sp.simplify((cd - expect).rewrite(sp.exp))
vals = {fj[(3, 0)]: 0.7 + 0.2 * I, fj[(2, 1)]: -0.3 + 0.5 * I, mu_: 1.3}
print('  (1/2){g, Rg} の ζ²ζbar² の係数 = 9|F30|² cot(3μ/2) + 3|F21|² cot(μ/2): 数値代入での差 %.2e' % abs(complex(sp.N((cd - expect).subs(vals)))))
print('  → K = -F4\' (対角) = -[9|F30|² cot(3πν) + 3|F21|² cot(πν)] (2J)²、∂ν/∂J = (1/2π)∂K/∂J = -(4/π)[9|h30000|² cot 3πν + 3|h21000|² cot πν]')
a3s, nu_ = sp.symbols('a3s nu_', positive=True)
dnu = -(4 / sp.pi) * (9 * (a3s / 48)**2 * cot(3 * sp.pi * nu_) + 3 * (a3s / 16)**2 * cot(sp.pi * nu_))
print('  1 個の薄い 6 極 (a = K2L β^{3/2}、h30000 = a/48、h21000 = a/16): 上の式と -a²/(64π)[cot 3πν + 3 cot πν] の差 = %s' % sp.simplify(dnu - (-(a3s**2 / (64 * sp.pi)) * (cot(3 * sp.pi * nu_) + 3 * cot(sp.pi * nu_)))))

print('\n== (i) 薄い四極のチューンシフトとクロマティシティ ==')
mm, kb, dm = sp.symbols('mm kb dm')
# 1 周の行列 R(μ) と薄い四極 [[1,0],[-k,1]] (正規化座標で k β): cos μ' = cos μ - (kβ/2) sin μ
Rm = sp.Matrix([[sp.cos(mm), sp.sin(mm)], [-sp.sin(mm), sp.cos(mm)]])
Kq_ = sp.Matrix([[1, 0], [-kb, 1]])       # kb = k β (正規化座標での kick: Δp^ = -kβ x^)
tr = (Rm * Kq_).trace() / 2
print('  cos μ\' = tr(R K)/2 = %s' % sp.simplify(tr))
mus = sp.symbols('mus')
solmu = sp.series(sp.acos(tr.subs(mm, mus)) - mus, kb, 0, 2).removeO()
print('  μ\' - μ の kβ の 1 次: %s → Δν = kβ/(4π)' % sp.simplify(solmu))
print('  6 極 (分散 η): 軌道が ηδ ずれて kick -(K2L/2)(x + ηδ)^2 の x の 1 次 = -K2L η δ x → k = K2L η δ、Δν = K2L η β δ/(4π) = (h11001/π) δ (h11001 = K2L ηβ/4)')

print('\n== (j) 1 本の共鳴: 固定点・島の幅、RF バケット ==')
Jj, d_, hh, ps = sp.symbols('J delta_ h_abs psi', positive=True)
Kres = 2 * sp.pi * d_ * Jj + 2 * hh * (2 * Jj)**sp.Rational(3, 2) * sp.cos(ps)     # ψ = 3φ + χ - 2πN n
dKdJ = sp.diff(Kres, Jj)
Jfp = sp.solve(sp.Eq(dKdJ.subs(sp.cos(ps), -1), 0), Jj)
print('  K = 2π δ J + 2|h| (2J)^{3/2} cos ψ (ψ = 3φ + 位相)。固定点 (cos ψ = -1、δ > 0 側) ∂K/∂J = 0: J = %s' % Jfp)
print('  (単独の 3νx 共鳴では J_fp = π²δ²/(18|h|²)。δ = ν - N/3)')
# 一般の 1 本の共鳴: K(I, θ) = π G I² + V cos θ (V > 0、G > 0)。θ = 0 が不安定、θ = π が安定な固定点
Gg, Vv, Ii, th = sp.symbols('G V I theta', positive=True)
Kp = sp.pi * Gg * Ii**2 + Vv * sp.cos(th)
Kux = Kp.subs({Ii: 0, th: 0})                 # 不安定固定点の K = V
Imax = sp.solve(sp.Eq(Kp.subs(th, sp.pi), Kux), Ii)
print('  K = πG I² + V cos θ。不安定固定点 (I, θ) = (0, 0) の K = %s。分離線 K = V は θ = π で I = %s' % (Kux, Imax))
print('   → 島の半幅 (チューン平面) W = G I_max = %s (docs/26 §3.5 の W = sqrt(2 |G0| V/π) と同じ)' % sp.simplify(Gg * Imax[0]))
# 1 本の 3νx 共鳴 (δ = ν - N/3): K = 2πδ J + c J² + 2|h|(2J)^{3/2} cos ψ の固定点 (cos ψ = -1): 2πδ + 2cJ - 6√2|h|√J = 0
cc_, sJ = sp.symbols('c_ sJ', real=True)
fp = sp.solve(sp.Eq(2 * sp.pi * d_ + 2 * cc_ * sJ**2 - 6 * sp.sqrt(2) * hh * sJ, 0), sJ)
print('  J² の項 c を含む固定点 (cos ψ = -1): √J = %s' % fp)
# RF バケット
psi, phs, aa_, kap = sp.symbols('psi phi_s a_ kappa_', real=True)
Vpot = aa_ * (-sp.cos(phs + psi) - psi * sp.sin(phs))
dVu = sp.simplify(Vpot.subs(psi, sp.pi - 2 * phs) - Vpot.subs(psi, 0))
print('  RF バケット: V(ψ) = a[-cos(φs + ψ) - ψ sin φs]。dV/dψ = %s' % sp.simplify(sp.diff(Vpot, psi)))
print('  V(π - 2φs) - V(0) = %s' % sp.simplify(dVu))
print('  (1/2) κ δ_rf² = V(ψu) - V(0) = a [2 cos φs - (π - 2φs) sin φs]: 差 %s' % sp.simplify(dVu - aa_ * (2 * sp.cos(phs) - (sp.pi - 2 * phs) * sp.sin(phs))))
