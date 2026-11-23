# -*- coding: utf-8 -*-
"""分割の一般論 (03_symplectic_integrators.md §3〜§6) の BCH 計算を、自由代数 (非可換な多項式) で行う。

2 つの生成元 A (ドリフトの流れ)、B (キックの流れ)。積 S = exp(..)exp(..)... の log を、語の長さ 5 まで、非可換な級数で求める。
Lie 多項式は Dynkin の写像 (長さ n の語 a1…an → (1/n)[a1,[a2,[…,an]]]) で、具体的な関数 (A = -p^2/2、B = -V(x)。1 自由度) のポアソン括弧に直す。
(a) Strang S2 = e^{A/2} e^{B} e^{A/2}: 3 次の項 (BCH の係数)。
(b) Yoshida (3 つの S2): 3 次が消える条件と、5 次の残り。
(c) 力勾配型: e^{aA} e^{B/2} e^{tA} e^{C} e^{tA} e^{B/2} e^{aA} (a = 1/2 - t、C = c3 W + c5 W5、W = [B,[A,B]]、W5 = [W,[A,B]])。
    3 次が消える条件から t、c3、5 次の x だけの項 (運動量を含まない項) が消える条件から c5 を決め、K64 の a、b、β、γ と比べる。
(d) 2 自由度で、K64 の tthin の補正キック (kord = 2、3) と、導いた c3、c5 から作る補正の関数形が一致するか。
"""
import itertools
import math
import sympy as sp
import _common as C

C.start(__file__)

MAXD = 5
x, p, y, q = sp.symbols('x p y q', real=True)
V1 = sp.Function('V')(x)


# ---- 自由結合代数 (語 → 係数) ----
def padd(a, b, sa=1, sb=1):
    r = {}
    for w, c in a.items():
        r[w] = r.get(w, 0) + sa * c
    for w, c in b.items():
        r[w] = r.get(w, 0) + sb * c
    return {w: sp.expand(c) for w, c in r.items() if sp.expand(c) != 0}


def pmul(a, b):
    r = {}
    for w1, c1 in a.items():
        for w2, c2 in b.items():
            if len(w1) + len(w2) <= MAXD:
                r[w1 + w2] = r.get(w1 + w2, 0) + c1 * c2
    return {w: sp.expand(c) for w, c in r.items() if sp.expand(c) != 0}


def pscale(a, s):
    return {w: sp.expand(c * s) for w, c in a.items() if sp.expand(c * s) != 0}


ONE = {'': sp.Integer(1)}


def comm(a, b):
    return padd(pmul(a, b), pmul(b, a), 1, -1)


def pexp(a):
    """exp(a) (a は定数項なし)。語の長さ MAXD まで"""
    r = dict(ONE)
    term = dict(ONE)
    for n in range(1, MAXD + 1):
        term = pscale(pmul(term, a), sp.Rational(1, n))
        r = padd(r, term)
    return r


def plog(s):
    """log(S) (S = 1 + X)"""
    X = padd(s, ONE, 1, -1)
    r = {}
    term = dict(ONE)
    for n in range(1, MAXD + 1):
        term = pmul(term, X)
        r = padd(r, pscale(term, sp.Rational((-1)**(n + 1), n)))
    return r


def degree_part(a, n):
    return {w: c for w, c in a.items() if len(w) == n}


Aw = {'A': sp.Integer(1)}
Bw = {'B': sp.Integer(1)}


def cA(c):
    return pexp(pscale(Aw, c))


def cB(c):
    return pexp(pscale(Bw, c))


# ---- 具体的な表現 (1 自由度。A = -p^2/2、B = -V(x)) ----
def PB(f, g):
    return sp.expand(sp.diff(f, x) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, x))


FUN = {'A': -p**2 / 2, 'B': -V1}


def theta(word):
    """右にネストしたポアソン括弧 [a1,[a2,[…,an]]]"""
    f = FUN[word[-1]]
    for ch in reversed(word[:-1]):
        f = PB(FUN[ch], f)
    return f


_theta_cache = {}


def concrete(lie_poly, n):
    """長さ n の Lie 多項式 (語の係数) → 関数 (Dynkin: (1/n) Σ c_w θ(w))"""
    tot = 0
    for w, c in degree_part(lie_poly, n).items():
        if w not in _theta_cache:
            _theta_cache[w] = theta(w)
        tot += c * _theta_cache[w]
    return sp.expand(tot / n)


def show(f, label):
    print('  %s = %s' % (label, sp.simplify(f) if f != 0 else 0))


def pcoeffs(f):
    """p のべきで整理した係数のリスト (p の次数, 係数)"""
    P = sp.Poly(sp.expand(f), p)
    return {m[0]: c for m, c in zip(P.monoms(), P.coeffs())}


# ---------------------------------------------------------------------------
print('== (a) Strang: S2 = e^{A/2} e^{B} e^{A/2} の log ==')
S2 = pmul(pmul(cA(sp.Rational(1, 2)), cB(1)), cA(sp.Rational(1, 2)))
L2 = plog(S2)
print('  1 次 = A + B か: %s' % (padd(degree_part(L2, 1), {'A': 1, 'B': 1}, 1, -1) == {}))
print('  2 次の項 (対称な積なので 0 のはず): %s' % (degree_part(L2, 2) == {}))
# 3 次: Lie 要素 [A,[A,B]]、[B,[B,A]] の係数
e1 = comm(Aw, comm(Aw, Bw))
e2 = comm(Bw, comm(Bw, Aw))
c1, c2 = sp.symbols('c1 c2')
sol = sp.solve([padd(degree_part(L2, 3), padd(pscale(e1, c1), pscale(e2, c2)), 1, -1).get(w, 0) for w in ('AAB', 'ABA', 'BAA', 'ABB', 'BAB', 'BBA')], [c1, c2], dict=True)
print('  3 次の項 = c1 [A,[A,B]] + c2 [B,[B,A]]: %s' % sol)
print('  → log S2 = (A + B) + (-1/24)[A,[A,B]] + (1/12)[B,[B,A]] + O(5 次)  (ε を付ければ ε^3)')

# ---------------------------------------------------------------------------
print('\n== (b) Yoshida: S4 = S2(w1) S2(w0) S2(w1)  (S2(w) は A, B を w 倍) ==')
w1 = sp.Float(1, 40) / (2 - sp.Float(2, 40)**(sp.Rational(1, 3)))
w0 = 1 - 2 * w1


def S2w(w):
    return pmul(pmul(cA(sp.Rational(1, 2) * w), cB(w)), cA(sp.Rational(1, 2) * w))


def maxabs(poly):
    return max([abs(float(sp.N(c, 30))) for c in poly.values()] + [0.0])


S4 = pmul(pmul(S2w(w1), S2w(w0)), S2w(w1))
L4 = plog(S4)
print('  w1 = 1/(2 - 2^(1/3)) = %.9f、w0 = 1 - 2 w1 = %.9f、2 w1 + w0 = %.1e、2 w1^3 + w0^3 = %.2e' % (float(w1), float(w0), float(2 * w1 + w0), float(2 * w1**3 + w0**3)))
print('  1 次 - (A + B) の最大の係数 %.1e、2 次 %.1e、3 次 %.1e (消える)、4 次 %.1e (対称な積なので 0)' % (
    maxabs(padd(degree_part(L4, 1), {'A': 1, 'B': 1}, 1, -1)), maxabs(degree_part(L4, 2)), maxabs(degree_part(L4, 3)), maxabs(degree_part(L4, 4))))
d5 = degree_part(L4, 5)
print('  5 次の項は残る (語の数 %d、最大の係数 %.2e。Yoshida は 4 次)。具体的な表現 (A = -p^2/2、B = -V) での 5 次の関数 (p のべきごとの係数):' % (len(d5), maxabs(d5)))
H5 = concrete(L4, 5)
pc = pcoeffs(H5)
for k in sorted(pc):
    print('    p^%d: %s' % (k, sp.N(sp.expand(pc[k]), 6)))

# ---------------------------------------------------------------------------
print('\n== (c) 力勾配型: e^{aA} e^{B/2} e^{tA} e^{C} e^{tA} e^{B/2} e^{aA} (a = 1/2 - t) ==')
t, c3, c5 = sp.symbols('t c3 c5')
W = comm(Bw, comm(Aw, Bw))                   # [B,[A,B]]
W5 = comm(W, comm(Aw, Bw))                   # [[B,[A,B]],[A,B]]
Cpoly = padd(pscale(W, c3), pscale(W5, c5))
eC = pexp(Cpoly)
a_ = sp.Rational(1, 2) - t
Sfg = cA(a_)
for fac in (cB(sp.Rational(1, 2)), cA(t), eC, cA(t), cB(sp.Rational(1, 2)), cA(a_)):
    Sfg = pmul(Sfg, fac)
Lfg = plog(Sfg)
print('  1 次 = A + B か: %s' % (padd(degree_part(Lfg, 1), {'A': 1, 'B': 1}, 1, -1) == {}))
d3 = degree_part(Lfg, 3)
cond3 = [sp.expand(padd(d3, padd(pscale(e1, c1), pscale(e2, c2)), 1, -1).get(w, 0)) for w in ('AAB', 'ABB')]
solc = sp.solve(cond3, [c1, c2], dict=True)[0]
print('  3 次の係数: c1 ([A,[A,B]]) = %s、c2 ([B,[B,A]]) = %s' % (sp.factor(solc[c1]), sp.factor(solc[c2])))
sol_t = sp.solve([solc[c1], solc[c2]], [t, c3], dict=True)
print('  3 次が消える解 (t, c3): %s' % sol_t)
pick = [s for s in sol_t if s[t].is_real and s[t] > 0 and s[t] < sp.Rational(1, 2)]
pick = pick[0] if pick else sol_t[0]
t_val, c3_val = sp.nsimplify(pick[t]), sp.nsimplify(pick[c3])
alpha = 1 / sp.sqrt(12)
beta = sp.Rational(1, 6) - alpha / 2
print('  t = %s = %.12f (K64: alb/al = alpha = 1/sqrt(12) = %.12f)' % (sp.simplify(t_val), float(t_val), float(alpha)))
print('  a = 1/2 - t = %.12f (K64: ala/al = alpha1 = 1/2 - alpha = %.12f)' % (float(sp.Rational(1, 2) - t_val), float(sp.Rational(1, 2) - alpha)))
print('  c3 = %s = %.12f、beta/2 (beta = 1/6 - alpha/2 = %.12f) = %.12f' % (sp.simplify(c3_val), float(c3_val), float(beta), float(beta / 2)))
# 5 次: 具体的な表現で p の 0 次 (x だけの項)
Lsub = {w: c.subs({t: t_val, c3: c3_val}) for w, c in Lfg.items()}
H5 = concrete({w: c for w, c in Lsub.items() if len(w) == 5}, 5)
pc = pcoeffs(H5)
print('  5 次の関数 (t, c3 を代入) の p のべきごとの係数 (c5 を含む):')
for k in sorted(pc):
    print('    p^%d: %s' % (k, sp.simplify(pc[k])))
sol5 = sp.solve(pc[0], c5)
gamma = sp.Rational(1, 40) - 1 / (24 * sp.sqrt(3))
print('  x だけの項 (p^0) を 0 にする c5 = %s = %.12f' % (sp.simplify(sol5[0]), float(sol5[0])))
print('  K64 の gamma = 1/40 - 1/(24 sqrt 3) = %.12f、c5/gamma = %.9f、c5 の候補: gamma/2 = %.12f' % (float(gamma), float(sol5[0] / gamma), float(gamma / 2)))
c5_val = sp.nsimplify(sol5[0])
H5r = sp.expand(H5.subs(c5, c5_val))
pcr = pcoeffs(H5r)
print('  c5 を入れた 5 次の関数 (残る誤差。p のべきごと):')
for k in sorted(pcr):
    print('    p^%d: %s' % (k, sp.simplify(pcr[k])))

# ---------------------------------------------------------------------------
print('\n== (d) 2 自由度、全運動量 p (質量に当たる) つき: K64 の tthin の補正キックと z の更新との比較 (kord = 2: 6 極、3: 8 極) ==')
print('   T = (px^2 + py^2)/(2p) として、導いた補正 f_C = c3 |∇V|^2/p - 2 c5 G/p^2 (G = V_j V_k V_jk) から: Δp = ∂f_C/∂x、Δq = Δp/p (q = px/p)、Δz = -∂f_C/∂p。')
X2, Y2 = sp.symbols('X Y', real=True)
K, pp = sp.symbols('K pp', positive=True)
for kord in (2, 3):
    Vm = sp.expand(K / sp.factorial(kord + 1) * sp.re(sp.expand((X2 - sp.I * Y2)**(kord + 1))))
    gradV = [sp.diff(Vm, X2), sp.diff(Vm, Y2)]
    Wf = sp.expand(gradV[0]**2 + gradV[1]**2)
    Gf = sp.expand(sum(gradV[j] * gradV[k2] * sp.diff(Vm, [X2, Y2][j], [X2, Y2][k2]) for j in range(2) for k2 in range(2)))
    fC = c3_val * Wf / pp - 2 * c5_val * Gf / pp**2
    dqx_mine = sp.diff(fC, X2) / pp
    dqy_mine = sp.diff(fC, Y2) / pp
    dz_mine = -sp.diff(fC, pp)
    al = 1
    aki = K / sp.factorial(kord) / pp                    # K64: aki = akf/pr*2 = ak/kord!/pr
    r = X2**2 + Y2**2
    cx = sp.expand((X2 - sp.I * Y2)**kord)
    rcx1 = X2 * sp.re(cx) + Y2 * sp.im(cx)
    if kord == 2:
        f1, f2, f3 = beta * 2 * al, 4 * gamma * al**2, 6 * gamma * al**2
        f4, f5 = sp.Rational(1, 2) * beta * al, f2
        rk, rk1 = r, 1
    else:
        f1 = beta * kord * al
        f2 = 2 * gamma * kord * (kord - 1) * al**2
        f3 = gamma * kord * (kord + 1) * al**2
        f4 = sp.Rational(1, 2) * beta * al
        f5 = 2 * gamma * kord * al**2
        rk1 = r**(kord - 2)
        rk = r * rk1
    dqx_k = aki**2 * (f1 * rk * X2 - aki * (f2 * rk1 * X2 * rcx1 + f3 * rk * sp.re(cx)))
    dqy_k = aki**2 * (f1 * rk * Y2 - aki * (f2 * rk1 * Y2 * rcx1 + f3 * rk * sp.im(cx)))
    dz_k = aki**2 * rk * (f4 * r - f5 * aki * rcx1)
    res = [sp.simplify(sp.expand(a_ - b_)) for a_, b_ in ((dqx_mine, dqx_k), (dqy_mine, dqy_k), (dz_mine, dz_k))]
    print('  kord = %d: 導いた補正と K64 の tthin の差: Δqx %s、Δqy %s、Δz %s' % (kord, res[0], res[1], res[2]))

print('\n== (e) CALC の tthine の分割 (ドリフト L/6、キック 1/2、ドリフト 2L/3、キック 1/2、ドリフト L/6。補正なし) の 3 次の係数 ==')
# 上の (c) の式で t = 1/3 (中央の 2/3 を 2 つに分ける)、c3 = 0
c1_th = sp.simplify(solc[c1].subs({t: sp.Rational(1, 3)}))
c2_th = sp.simplify(solc[c2].subs({t: sp.Rational(1, 3), c3: 0}))
print('  t = 1/3、c3 = 0: c1 ([A,[A,B]]) = %s、c2 ([B,[B,A]]) = %s' % (c1_th, c2_th))
print('  (Strang は c1 = -1/24、c2 = +1/12。tthine は [B,[B,A]] の項だけ消した 2 次の分割)')
