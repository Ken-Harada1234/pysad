# -*- coding: utf-8 -*-
"""14_emittance.md の補足の代数 (sympy): tsymp の収束、btr の対角の置き換え、縦 2x2 の κ と νs、有限の νs の効果で使う 1 次の摂動の式。

(a) tsymp: t' = t (3/2 - t^s t/2) → t'^s t' = 1 - 3E²/4 + E³/4 (E = t^s t - 1)。E^s = E の確認 (2 次元の具体例で)。
(b) btr の対角: Λ' = P(μ) + rd、rd = -a P(μ) + S (S は対称でトレース 0) の 2x2 ブロックで、
    「和」の変数 (Σ_XX + Σ_PP)/√2 の係数 1 - (Λ'⊗Λ' の対角) と、tecalc の式の比較。
(c) 縦の 2x2: [[1, R56],[κ, 1 + κ R56]] の固有値と νs: κ R56 = -4 sin²(μs/2)。
(d) 回転で 1 周の変換を打ち消す座標で (Σ)_(sum) だけを残す理由: 2x2 の Σ → ΛΣΛ^T の和・差・交差の 3x3 行列。
"""
import sympy as sp
import _common as C

C.start(__file__)

print('== (a) tsymp ==')
E = sp.symbols('E')
lhs = sp.expand((1 - E / 2)**2 * (1 + E))
print('(1 - E/2)² (1 + E) =', lhs, ' → 1 - 3E²/4 + E³/4 と一致:', sp.simplify(lhs - (1 - sp.Rational(3, 4) * E**2 + E**3 / 4)) == 0)
# E^s = E: 4x4 の一般の行列 t で t^s = -J t^T J。(t^s t)^s = t^s t
n = 4
t = sp.Matrix(n, n, lambda i, j: sp.Symbol('t%d%d' % (i, j)))
J = sp.zeros(n, n)
J[0, 1], J[1, 0], J[2, 3], J[3, 2] = 1, -1, 1, -1
ts = lambda a: -J * a.T * J
Ee = ts(t) * t
print('(t^s t)^s - t^s t = 0 :', sp.simplify(ts(Ee) - Ee) == sp.zeros(n, n))
print('(AB)^s = B^s A^s (A,B は 4x4 の記号行列) の 1 成分の確認:')
B = sp.Matrix(n, n, lambda i, j: sp.Symbol('b%d%d' % (i, j)))
print('  ', sp.simplify(ts(t * B) - ts(B) * ts(t)) == sp.zeros(n, n))

print('\n== (b) btr の対角ブロック ==')
a, mu, s1, s2 = sp.symbols('a mu s1 s2', real=True)
c, s = sp.cos(mu), sp.sin(mu)
P = sp.Matrix([[c, s], [-s, c]])                 # tecalc: rdb(i,i+6) = cos, rdb(i,i+7) = sin、rdb(i+1,i+6) = -sin、rdb(i+1,i+7) = cos
Sm = sp.Matrix([[s1, s2], [s2, -s1]])
rd = -a * P + Sm
expr = (-(rd[0, 0]**2 + rd[0, 1]**2 + rd[1, 0]**2 + rd[1, 1]**2) / 2
        - c * (rd[0, 0] + rd[1, 1]) - s * (rd[0, 1] - rd[1, 0]))
print('tecalc の btr(k1,k1) (rd = -aP + S) =', sp.simplify(sp.expand(expr)))
expr0 = sp.simplify(sp.expand(expr.subs({s1: 0, s2: 0})))
print('S = 0 のとき =', expr0, ' (= 2a - a² = 1 - (1-a)²)')
# 2x2 の Σ → Λ Σ Λ^T、Λ = P + rd = (1-a)P + S の 3x3 行列 (基底: 和 u = (Σ11 + Σ22)/√2、差 v = (Σ22 - Σ11)/√2、交差 w = √2 Σ12)
Lam = P + rd
sx, sy, sz = sp.symbols('sx sy sz')
Sig = sp.Matrix([[sx, sz], [sz, sy]])
Sig2 = Lam * Sig * Lam.T
def basis(S_):
    return sp.Matrix([(S_[0, 0] + S_[1, 1]) / sp.sqrt(2), (S_[1, 1] - S_[0, 0]) / sp.sqrt(2), sp.sqrt(2) * S_[0, 1]])
vec_in = basis(Sig)
vec_out = basis(Sig2)
u, v, w = sp.symbols('u v w')
sol = sp.solve([vec_in[0] - u, vec_in[1] - v, vec_in[2] - w], [sx, sy, sz], dict=True)[0]
Mrot = sp.Matrix([[sp.simplify(sp.diff(vec_out[i].subs(sol), q)) for q in (u, v, w)] for i in range(3)])
print('Σ → ΛΣΛ^T の (和,差,交差) の 3x3 行列で、S = 0 のとき:')
sp.pprint(sp.simplify(Mrot.subs({s1: 0, s2: 0})))
print('和-和の成分 (1 - 係数) =', sp.simplify(1 - Mrot[0, 0].subs({s1: 0, s2: 0})), ' = 1 - (1-a)² =', sp.expand(1 - (1 - a)**2))
print('和-和の成分、S あり: 1 - M[0,0] =', sp.simplify(sp.expand(1 - Mrot[0, 0])))
print('  tecalc の式との差 (S あり) =', sp.simplify(sp.expand(1 - Mrot[0, 0] - expr)))

print('\n== (c) 縦の 2x2 ==')
R56, kap, lamb = sp.symbols('R56 kappa lamb')
M2 = sp.Matrix([[1, R56], [kap, 1 + kap * R56]])
print('det =', sp.simplify(M2.det()), '、trace =', M2.trace(), ' → 2 cos μ = 2 + κ R56')
mus = sp.symbols('mu_s')
print('κ R56 = 2cos μ - 2 = -4 sin²(μ/2):', sp.simplify(2 * sp.cos(mus) - 2 + 4 * sp.sin(mus / 2)**2) == 0)
print('小さい μ: -4 sin²(μ/2) の級数 =', sp.series(-4 * sp.sin(mus / 2)**2, mus, 0, 6))

print('\n== (d) 1 周の変換行列が rotation の 2x2 ブロックだけの Λ = e^{-a}P の場合: Σ の方程式の固有値 ==')
Lam1 = (1 - a) * P
vec2 = basis(Lam1 * Sig * Lam1.T)
sol1 = sp.solve([basis(Sig)[0] - u, basis(Sig)[1] - v, basis(Sig)[2] - w], [sx, sy, sz], dict=True)[0]
M1 = sp.Matrix([[sp.simplify(sp.diff(vec2[i].subs(sol1), q)) for q in (u, v, w)] for i in range(3)])
sp.pprint(M1)
print('固有値:', [sp.simplify(x) for x in M1.eigenvals().keys()])
