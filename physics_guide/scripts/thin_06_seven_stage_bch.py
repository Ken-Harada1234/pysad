# -*- coding: utf-8 -*-
r"""07_thin_multipole.md の §5・§6.2・§7・§9.1 の途中式の補完 (sympy、軽い常微分方程式)。

(A) 1 自由度 H = p^2/2 + V(x) の 7 段 (ドリフト aℓ、キック、ドリフト tℓ、補正、ドリフト tℓ、キック、ドリフト aℓ) の
    x, p を h のべきで段ごとに書き、厳密解 (x'' = -V'(x) の Picard 反復) との差から h^0〜h^5 を出す (T-01)。
(B) 7 因子の積の BCH を三つ組の BCH に分けて足し合わせる (T-02)。自由代数 (非可換記号) で検算。
    リー演算子の積の順序を逆にしても同じ (T-03) ことも、順序を逆にした積で確かめる。
(C) CALC の x 面 2x2 行列 M_x を 2 段に分けた積と、l1 = L/6、l2 = 2L/3 を入れた成分 (T-04)。
(D) 非線形フリンジの b_s の項の大きさ: ローレンツ力の積分で b_s を入れた場合・落とした場合の差 (T-05)。
"""
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

# ---------------------------------------------------------------- (A)
print('== (A) 1 自由度の 7 段と厳密解の差 (T-01) ==')
print('   H = p^2/2 + V(x)。初期値 x = 0、p = 0。V_n = V の n 階微分 (x = 0 で)。a = 1/2 - t。補正 C: dp = ∂f_C/∂x、f_C = c3 h^3 V\'^2 - 2 c5 h^5 V\'^2 V\'\'')
h, t, c3, c5 = sp.symbols('h t c3 c5')
N = 5
V = sp.symbols('V0:%d' % (N + 4))
a_ = sp.Rational(1, 2) - t


def trunc(e, n=N):
    e = sp.expand(e)
    return sum(sp.expand(e.coeff(h, i)) * h**i for i in range(n + 1))


def Vder(n, d):       # V^{(n)}(x = d) = Σ V_{n+j} d^j / j!
    return sum(V[n + j] * d**j / sp.factorial(j) for j in range(N + 2 - n))


def stage(state, kind, arg):
    xs, ps = state
    if kind == 'D':
        return trunc(xs + arg * h * ps), ps
    if kind == 'K':
        return xs, trunc(ps - arg * h * Vder(1, xs))
    v1, v2, v3 = Vder(1, xs), Vder(2, xs), Vder(3, xs)
    dp = c3 * h**3 * 2 * v1 * v2 - 2 * c5 * h**5 * (2 * v1 * v2 * v2 + v1**2 * v3)
    return xs, trunc(ps + dp)


seq = [('D', a_, 'ドリフト aℓ'), ('K', sp.Rational(1, 2), 'キック 1/2'), ('D', t, 'ドリフト tℓ'),
       ('C', 0, '補正 C'), ('D', t, 'ドリフト tℓ'), ('K', sp.Rational(1, 2), 'キック 1/2'), ('D', a_, 'ドリフト aℓ')]
st = (sp.Integer(0), sp.Integer(0))
for i, (kd, ar, name) in enumerate(seq, 1):
    st = stage(st, kd, ar)
    xs, ps = st
    print('  段 %d (%s):' % (i, name))
    print('     x = %s' % sp.collect(sp.expand(sp.series(xs, h, 0, 6).removeO()), h))
    print('     p = %s' % sp.collect(sp.expand(sp.series(ps, h, 0, 6).removeO()), h))
xs7, ps7 = st
# 厳密解: Picard 反復 x(s) = -∫_0^s (s - u) V'(x(u)) du。s のべきで 6 次まで
s, u = sp.symbols('s u')
xe = sp.Integer(0)
print("  厳密解 (Picard 反復。x_{n+1}(s) = -∫_0^s (s - u) V'(x_n(u)) du、s^6 まで):")


def trunc_s(e, n=6):
    e = sp.expand(e)
    return sum(sp.expand(e.coeff(s, i)) * s**i for i in range(n + 1))


for it in range(1, 5):
    xu = xe.subs(s, u)
    integrand = sum(V[1 + j] * xu**j / sp.factorial(j) for j in range(0, 5))
    integrand = sp.expand(integrand)
    integrand = sum(integrand.coeff(u, i) * u**i for i in range(0, 5))
    xe = trunc_s(-sp.integrate(sp.expand((s - u) * integrand), (u, 0, s)))
    print('    x_%d(s) = %s' % (it, sp.collect(xe, s)))
pe = sp.diff(xe, s)
xe_h, pe_h = [sp.expand(q.subs(s, h)) for q in (xe, pe)]
print('  厳密解 p(h) = %s' % sp.collect(trunc(pe_h), h))
dx = sp.expand(xs7 - trunc(xe_h))
dp = sp.expand(ps7 - trunc(pe_h))
print('  (7 段の結果) - (厳密解)、t と c3 と c5 は記号のまま:')
for n in range(0, 6):
    print('    h^%d: x の差 = %s   |   p の差 = %s' % (n, sp.factor(dx.coeff(h, n)), sp.factor(dp.coeff(h, n))))
print('  h^3 の p の差 = V1 V2 (2 c3 - 1/6 + t/2) か:', sp.simplify(dp.coeff(h, 3) - V[1] * V[2] * (2 * c3 - sp.Rational(1, 6) + t / 2)) == 0)
# t = 1/sqrt(12)
t_val = 1 / sp.sqrt(12)
c3_val = sp.Rational(1, 12) - t_val / 4
print('  t = 1/√12 = √3/6、c3 = 1/12 - t/4 = 1/12 - √3/24 を入れる:')
sub3 = {t: t_val, c3: c3_val}
for n in range(0, 5):
    print('    h^%d: x の差 = %s   |   p の差 = %s' % (n, sp.simplify(dx.coeff(h, n).subs(sub3)), sp.simplify(dp.coeff(h, n).subs(sub3))))
e5 = sp.expand(sp.simplify(dp.coeff(h, 5).subs(sub3)))
print('    h^5: p の差 =', sp.collect(sp.factor(e5), [V[1]]))
cand = (sp.Rational(1, 40) - sp.sqrt(3) / 72 - 2 * c5) * (V[1]**2 * V[3] + 2 * V[1] * V[2]**2)
print('    h^5 の p の差 = (1/40 - √3/72 - 2 c5)(V1² V3 + 2 V1 V2²) か:', sp.simplify(e5 - cand) == 0)
c5_val = sp.Rational(1, 80) - sp.sqrt(3) / 144
print('    c5 = 1/80 - √3/144 = (9 - 5√3)/720 を入れると h^5 の p の差:', sp.simplify(e5.subs(c5, c5_val)),
      ' x の差 (h^4, h^5):', sp.simplify(dx.coeff(h, 4).subs(sub3)), sp.simplify(dx.coeff(h, 5).subs(sub3).subs(c5, c5_val)))
print('    数値: c3 = %.10f (β/2、β = 1/6 - 1/√48 = %.10f)、c5 = %.10e (γ/2、γ = 1/40 - 1/(24√3) = %.10e)' % (
    float(c3_val), float(sp.Rational(1, 6) - 1 / sp.sqrt(48)), float(c5_val), float(sp.Rational(1, 40) - 1 / (24 * sp.sqrt(3)))))
q = sp.Symbol('q'); Vq = sp.Function('Vq')(q)
print('    d/dq (V\'^2 V\'\') = %s' % sp.expand(sp.diff(sp.diff(Vq, q)**2 * sp.diff(Vq, q, 2), q)))

# ---------------------------------------------------------------- (B)
print('\n== (B) 7 因子の BCH を三つ組の BCH に分けて足し合わせる (T-02、T-03) ==')
A, B = sp.symbols('A B', commutative=False)
hh = sp.Symbol('hh')


def tr(e, n=3):
    e = sp.expand(e)
    return sum(e.coeff(hh, i) * hh**i for i in range(n + 1))


def mexp(X):
    return tr(1 + X + X * X / 2 + X * X * X / 6)


def mlog(P):
    Y = tr(P - 1)
    return tr(Y - Y * Y / 2 + Y * Y * Y / 3)


def cm(X, Y):
    return sp.expand(X * Y - Y * X)


def nest(X, Y, Z):   # [X,[Y,Z]]
    return cm(X, cm(Y, Z))


def same(X, Y):
    return sp.expand(X - Y) == 0


# 三つ組の公式 (X1 = 端の因子、Z = 中央): log(e^{X1} e^{Z} e^{X1}) = 2X1 + Z + (1/6)[Z,[Z,X1]] - (1/6)[X1,[X1,Z]]
X1, Z1 = sp.symbols('X1 Z1', commutative=False)
lhs = mlog(mexp(hh * X1) * mexp(hh * Z1) * mexp(hh * X1))
rhs = 2 * hh * X1 + hh * Z1 + hh**3 * (nest(Z1, Z1, X1) / 6 - nest(X1, X1, Z1) / 6)
print('  (i) 三つ組: log(e^{X}e^{Z}e^{X}) = 2X + Z + (1/6)[Z,[Z,X]] - (1/6)[X,[X,Z]] + O(h^5) (h^2 の項は 0):', same(tr(lhs - rhs), 0))
print('      h^2 の係数 =', sp.expand(lhs.coeff(hh, 2)))
# 段 1: Q = e^{tA} e^{C} e^{tA}、C = c3 h^3 [B,[A,B]]
Cc = c3 * hh**3 * nest(B, A, B)
Q = mlog(mexp(t * hh * A) * mexp(Cc) * mexp(t * hh * A))
print('  (ii) Q = log(e^{tA}e^{C}e^{tA}) = 2t hA + C + O(h^5):', same(tr(Q - (2 * t * hh * A + Cc)), 0))
# 段 2: W = log(e^{B/2} e^{Q} e^{B/2})
W = mlog(mexp(hh * B / 2) * mexp(Q) * mexp(hh * B / 2))
Wc = hh * B + 2 * t * hh * A + Cc + hh**3 * (t**2 / 3 * nest(A, A, B) - t / 12 * nest(B, B, A))
print('  (iii) W = B + 2tA + C + (t²/3)[A,[A,B]] - (t/12)[B,[B,A]] (h をつけて):', same(tr(W - Wc), 0))
# 段 3: log(e^{aA} e^{W} e^{aA})
P = mlog(mexp(a_ * hh * A) * mexp(W) * mexp(a_ * hh * A))
cA = (12 * t**2 - 1) / 24
cB = sp.Rational(1, 12) - t / 4 - c3
Pc = hh * (A + B) + hh**3 * (cA * nest(A, A, B) + cB * nest(B, B, A))
print('  (iv) 全体: log(P) = h(A + B) + h^3 [ (12t²-1)/24 [A,[A,B]] + (1/12 - t/4 - c3) [B,[B,A]] ] + O(h^5):', same(tr(P - Pc), 0))
print('       h^2 の係数 =', sp.expand(P.coeff(hh, 2)))
# 手計算の途中: 外側の段で (1/6)[W,[W,aA]] - (1/6)[aA,[aA,W]] の W を B + 2tA で評価
Wl = B + 2 * t * A
mid = a_ * (nest(Wl, Wl, A) / 6) - a_**2 * nest(A, A, Wl) / 6
# 基底 [A,[A,B]]、[B,[B,A]] に分解
exp_mid = sp.expand(mid)
cA_tot = sp.simplify(t**2 / 3 - a_ * t / 3 - a_**2 / 6)
cB_tot = sp.simplify(-t / 12 + a_ / 6 - c3)
print('  (v) 手計算の和: [A,[A,B]]: t²/3 - a t/3 - a²/6 =', sp.factor(cA_tot), '、[B,[B,A]]: -t/12 + a/6 - c3 =', sp.factor(cB_tot),
      '(a = 1/2 - t)。(iv) と一致:', sp.simplify(cA_tot - cA) == 0 and sp.simplify(cB_tot - cB) == 0)
# 外側の公式を W の主要部 B + 2tA で評価した 3 次の項が、手計算 (t²/3 の寄与を除く) と一致
manual = hh**3 * (-t / 12 * nest(B, B, A) + t**2 / 3 * nest(A, A, B) + mid)
print('  (vi) 手計算の 3 次の全体 (W の 3 次の項 + 外側の公式) と (iv) の差 (C を除く):',
      same(sp.expand(tr(manual - hh**3 * (cA * nest(A, A, B) + (cB + c3) * nest(B, B, A)))), 0))
# t の値
print('  (vii) t = 1/√12 で第 1 係数 =', sp.simplify(cA.subs(t, 1 / sp.sqrt(12))), '、c3 = 1/12 - t/4 で第 2 係数 =', sp.simplify(cB.subs({t: 1 / sp.sqrt(12), c3: sp.Rational(1, 12) - 1 / (4 * sp.sqrt(12))})))
print('  (viii) CALC 型: t = 1/3、C = 0 (c3 = 0): 第 1 係数 =', cA.subs(t, sp.Rational(1, 3)), '、第 2 係数 =', cB.subs({t: sp.Rational(1, 3), c3: 0}),
      '。Strang (a = 1/2、中央なし、t = 0 で C = 0 の 3 因子 e^{A/2} e^{B} e^{A/2}) は別の式:')
S2 = mlog(mexp(hh * A / 2) * mexp(hh * B) * mexp(hh * A / 2))
print('        Strang の 3 次 = -(1/24)[A,[A,B]] + (1/12)[B,[B,A]] か:', same(S2.coeff(hh, 3), -nest(A, A, B) / 24 + nest(B, B, A) / 12))
# 符号反転の具体確認: ψ: X1 X2 → X2 X1 (反準同型)。ψ([A,[A,B]]) = [[B,A],A]... を直接作る
# T-03: 積の順序の規約を逆にすることは、非可換多項式の反準同型 ψ(X1 X2) = X2 X1 に当たる。
# 因子の並びは回文なので ψ(積) = 積。ψ は [X,Y] → -[X,Y] (2 重の交換子では 2 回反転して不変)。
def psi(e):
    e = sp.expand(e)
    out = 0
    for term in sp.Add.make_args(e):
        c, nc = term.args_cnc()
        out += sp.Mul(*c) * sp.Mul(*reversed(nc))
    return sp.expand(out)
prod7 = tr(mexp(a_ * hh * A) * mexp(hh * B / 2) * mexp(t * hh * A) * mexp(Cc) * mexp(t * hh * A) * mexp(hh * B / 2) * mexp(a_ * hh * A))
print('  (ix) 7 因子の積 (3 次まで) は ψ で不変 (因子の並びが回文):', same(psi(prod7), prod7))
print('       反準同型 ψ(X1 X2) = X2 X1 で ψ([A,[A,B]]) = [A,[A,B]]:', same(psi(nest(A, A, B)), nest(A, A, B)),
      '、ψ([B,[B,A]]) = [B,[B,A]]:', same(psi(nest(B, B, A)), nest(B, B, A)),
      '、ψ([A,B]) = -[A,B]:', same(psi(cm(A, B)), -cm(A, B)))

# ---------------------------------------------------------------- (C)
print('\n== (C) CALC の x 面 2x2 行列 M_x を 2 段に分けた積 (T-04) ==')
l1, l2, h1, h2, Ls, hs = sp.symbols('l1 l2 h1 h2 L h')


def Dm(l):
    return sp.Matrix([[1, l], [0, 1]])


def Km(hv):
    return sp.Matrix([[1, 0], [-hv, 1]])


N1 = sp.expand(Km(h1) * Dm(l1))
N2 = sp.expand(Dm(l2) * N1)
N3 = sp.expand(Km(h2) * N2)
M = sp.expand(Dm(l1) * N3)
for nm, mat in (('N1 = K(h1) D(l1)', N1), ('N2 = D(l2) N1', N2), ('N3 = K(h2) N2', N3), ('M = D(l1) N3', M)):
    print('  %s =' % nm, mat.tolist())
Mtext = sp.Matrix([
    [1 - h1 * (l1 + l2) - h2 * l1 + h1 * h2 * l1 * l2, 2 * l1 + l2 - h1 * l1 * (l1 + l2) - h2 * l1 * (l1 + l2) + h1 * h2 * l1**2 * l2],
    [-h1 - h2 + h1 * h2 * l2, 1 - h1 * l1 - h2 * (l1 + l2) + h1 * h2 * l1 * l2]])
print('  M の 4 成分が本文の式と一致 (差の最大):', sp.simplify(M - Mtext) == sp.zeros(2, 2))
Mh = sp.simplify(M.subs({h1: hs, h2: hs, l1: Ls / 6, l2: 2 * Ls / 3}))
print('  h1 = h2 = h、l1 = L/6、l2 = 2L/3 を入れた各成分:')
for (i, j) in ((0, 0), (0, 1), (1, 0), (1, 1)):
    print('    M_%d%d = %s' % (i + 1, j + 1, sp.expand(Mh[i, j])))
print('    行列式 =', sp.simplify(Mh.det()))
Mt = sp.Matrix([[1 - Ls * hs + Ls**2 * hs**2 / 9, Ls - 5 * Ls**2 * hs / 18 + Ls**3 * hs**2 / 54],
                [-2 * hs + 2 * Ls * hs**2 / 3, 1 - Ls * hs + Ls**2 * hs**2 / 9]])
print('  本文の M_x と一致:', sp.simplify(Mh - Mt) == sp.zeros(2, 2))
print('  途中の成分 (l1 = L/6、l2 = 2L/3):  l1 + l2 = %s、2 l1 + l2 = %s、l1² l2 = %s、l1 l2 = %s、l1 (l1 + l2) = %s' % (
    sp.Rational(5, 6) * Ls, Ls, Ls**3 / 54, Ls**2 / 9, sp.Rational(5, 36) * Ls**2))

# ---------------------------------------------------------------- (D)
print('\n== (D) 非線形フリンジで b_s の項を入れた場合と落とした場合 (T-05) ==')
from scipy.integrate import solve_ivp
xs_, ys_ = sp.symbols('xs ys', real=True)


def make_fields(m):
    ww = xs_ + sp.I * ys_
    Phi0 = sp.im(sp.expand(ww**m)) / sp.factorial(m)
    psi_ = sp.expand(ww * sp.conjugate(ww)) * sp.im(sp.expand(ww**m)) / (4 * sp.factorial(m + 1))
    f = sp.lambdify((xs_, ys_), [sp.diff(Phi0, xs_), sp.diff(Phi0, ys_), sp.diff(psi_, xs_), sp.diff(psi_, ys_), Phi0], 'math')

    def fld(x_, y_, k0, k1, k2, use_bs=True, use_k2=True):
        p0x, p0y, psx, psy, p0 = f(x_, y_)
        kk2 = k2 if use_k2 else 0.0
        return k0 * p0x - kk2 * psx, k0 * p0y - kk2 * psy, (k1 * p0 if use_bs else 0.0)
    return fld


def rhs_factory(fld, kappa0, d, **kw):
    def rhs(s_, v):
        x_, px_, y_, py_, z_ = v
        th = math.tanh(s_ / d)
        sech2 = 1 - th * th
        k0 = kappa0 * (1 + th) / 2
        k1 = kappa0 * sech2 / (2 * d)
        k2 = -kappa0 * th * sech2 / d**2
        bx, by, bs = fld(x_, y_, k0, k1, k2, **kw)
        pz = math.sqrt(max(1.0 - px_**2 - py_**2, 1e-12))
        return [px_ / pz, (py_ * bs - pz * by) / pz, py_ / pz, (pz * bx - px_ * bs) / pz, 1.0 - 1.0 / pz]
    return rhs


S = 0.06
print('   入口 (tanh の端、-S < s < S、S = 0.06 m、p = 1)。「線形」= κ(s) Im w^m/m! だけの場。')
print('   全て: κ\'\' の力と b_s の両方。b_s なし: b_s = 0 (κ\'\' の力は残す)。κ\'\' なし: κ\'\' の力 = 0 (b_s は残す)。')
print('   Δ_全て、Δ_{b_s なし}、Δ_{κ\'\' なし} = (その場の積分) - (線形の積分)。成分 (x, px, y, py) の最大値を比べる。')
print('   m  κ0      r[mm] d[mm] | |Δ_全て|    |Δ_bsなし|   |Δ_κ\'\'なし|  | 位置: Δ_bsなし/Δ_全て (x,y)  | 運動量: Δ_bsなし/Δ_全て (px,py)')
for m, kappa0 in ((3, 400.0), (4, 4.0e4)):
    fld = make_fields(m)
    for scale in (1.0, 0.5):
        for d in (0.0005, 0.002):
            c0 = np.array([0.004, 0.0003, -0.003, 0.0002, 0.0]) * scale
            kw = dict(method='DOP853', rtol=1e-13, atol=1e-17)

            def run(**k):
                return solve_ivp(rhs_factory(fld, kappa0, d, **k), (-S, S), c0, **kw).y[:, -1]
            full = run()
            nobs = run(use_bs=False)
            nok2 = run(use_k2=False)
            # 線形な場
            f0 = sp.lambdify((xs_, ys_), [sp.diff(sp.im(sp.expand((xs_ + sp.I * ys_)**m)) / sp.factorial(m), xs_),
                                          sp.diff(sp.im(sp.expand((xs_ + sp.I * ys_)**m)) / sp.factorial(m), ys_)], 'math')

            def flin(x_, y_, k0, k1, k2, **kk):
                g = f0(x_, y_)
                return k0 * g[0], k0 * g[1], 0.0
            lin = solve_ivp(rhs_factory(flin, kappa0, d), (-S, S), c0, **kw).y[:, -1]
            D_full, D_nobs, D_nok2 = full - lin, nobs - lin, nok2 - lin
            ix, ip = [0, 2], [1, 3]
            print('   %d  %-7g %4.1f  %4.1f  | %.3e  %.3e  %.3e | %.3f  | %.3f' % (
                m, kappa0, 4.0 * scale, d * 1e3, np.abs(D_full[[0, 1, 2, 3]]).max(), np.abs(D_nobs[[0, 1, 2, 3]]).max(), np.abs(D_nok2[[0, 1, 2, 3]]).max(),
                np.abs(D_nobs[ix]).max() / np.abs(D_full[ix]).max(), np.abs(D_nobs[ip]).max() / np.abs(D_full[ip]).max()))
print('   局所の比 (b_s の力 / κ\'\' の力) の目安: 4 θ d / r (θ = p_⊥/p ~ 3e-4、d = 0.5〜2 mm、r = 4 mm) =',
      ', '.join('%.1e' % (4 * 3e-4 * dd / 4e-3) for dd in (0.5e-3, 2e-3)), ' (無次元。d、r は m で)')

print('\n   端のすぐ後 (s = 10 d) での位置のずれ Δ(x + iy) (複素数。単位 m)。(7.17) = -κ0/(4 m! p) w (w^m/(m+1) - w̄^m)、w は端 (s = 0) での位置 (p = 1)')
print('   m  r[mm] d[mm] | |Δ_全て|     |Δ_bsなし|    |(7.17)|    | (Δ_全て - (7.17))/(7.17) | (Δ_bsなし - (7.17))/(7.17) | 運動量のずれの b_s の寄与 |Δp_全て - Δp_bsなし|/|Δp_全て|')
for m, kappa0 in ((3, 400.0), (4, 4.0e4)):
    fld = make_fields(m)
    f0 = sp.lambdify((xs_, ys_), [sp.diff(sp.im(sp.expand((xs_ + sp.I * ys_)**m)) / sp.factorial(m), xs_),
                                  sp.diff(sp.im(sp.expand((xs_ + sp.I * ys_)**m)) / sp.factorial(m), ys_)], 'math')

    def flin(x_, y_, k0, k1, k2, **kk):
        g = f0(x_, y_)
        return k0 * g[0], k0 * g[1], 0.0
    for scale in (1.0, 0.5):
        for d in (0.0005, 0.002):
            c0 = np.array([0.004, 0.0003, -0.003, 0.0002, 0.0]) * scale
            kw = dict(method='DOP853', rtol=1e-13, atol=1e-17)
            send = 10 * d

            def run(fl, **k):
                return solve_ivp(rhs_factory(fl, kappa0, d, **k), (-S, send), c0, **kw).y[:, -1]
            full, nobs, lin = run(fld), run(fld, use_bs=False), run(flin)
            # 端 (s = 0) での位置 (κ0 の効果の前: 線形場の積分で s = 0 まで)
            e0 = solve_ivp(rhs_factory(flin, kappa0, d), (-S, 0.0), c0, **kw).y[:, -1]
            w = complex(e0[0], e0[2])
            dw = -kappa0 / (4 * math.factorial(m)) * w * (w**m / (m + 1) - w.conjugate()**m)
            Df = complex(full[0] - lin[0], full[2] - lin[2])
            Dn = complex(nobs[0] - lin[0], nobs[2] - lin[2])
            pf = complex(full[1] - lin[1], full[3] - lin[3])
            pn = complex(nobs[1] - lin[1], nobs[3] - lin[3])
            print('   %d  %4.1f  %4.1f  | %.3e  %.3e  %.3e | %+.3f | %+.3f | %.3f' % (
                m, 4.0 * scale, d * 1e3, abs(Df), abs(Dn), abs(dw), abs(Df - dw) / abs(dw), abs(Dn - dw) / abs(dw), abs(pf - pn) / abs(pf)))
print('   (b_s なしの積分でも位置のずれは (7.17) と一致し、b_s は位置のずれに関与しない。b_s の有無で変わるのは運動量 (K の 1 次) で、運動量は §9.3 の正準変換から決める)')
