# -*- coding: utf-8 -*-
"""ANGLE ≠ 0 の MULT の A_s の展開 (08_mult.md §6): ヘルプの g_kn と K64 (`gkninit`) の g_kn を、Maxwell 方程式から確かめる。

曲線座標 (x, y, s)、h = 1/ρ、スケール因子 (1, 1, 1 + h x)。s によらない場で A = A_s ŝ だけとすると、
  B_x = ∂_y A、B_y = -(1/(1+hx)) ∂_x((1+hx) A)  (∇·B = 0 は恒等的に満たされる)
  (∇×B)_s = ∂_x B_y - ∂_y B_x = 0  →  ∂_x[∂_x A + A/(ρ+x)] + ∂_y² A = 0     (M)
H の中の項は -(1+x/ρ) a_s = 「A_s(x, y)」 (ヘルプ)。A_s := (ρ+x) a_s/ρ とおく。(M) は (1+hx) a_s の方程式にすると
  ∂_x[ (1/(1+hx)) ∂_x((1+hx) a) ] + ∂_y² a = 0 、 a = -A_s/(1+hx)  (符号は線形なので無関係)。
ここでは F = (1+hx) a (= ∝ A_s) について解く:  (M') ∂_x[(1/(1+hx)) ∂_x F] + ∂_y²[F/(1+hx)] = 0
"""
import sympy as sp
import _common as C

C.start(__file__)
x, y, rho = sp.symbols('x y rho', positive=True)
z = x + sp.I * y
h = 1 / rho


def residual(Fexpr):
    """(M'): ∂_x[(1/(1+hx)) ∂_x F] + ∂_y²[F/(1+hx)]"""
    g = 1 + h * x
    return sp.diff(sp.diff(Fexpr, x) / g, x) + sp.diff(Fexpr / g, y, 2)


def gk_help(n, k):
    """ヘルプ (multipole_with_nonzero_ANGLE): g_kn = -(2k-1)!!(2k-3)!!(n+1)!/(8^k (n+k+1)! k!)"""
    return -sp.factorial2(2 * k - 1) * sp.factorial2(2 * k - 3) * sp.factorial(n + 1) / (8**k * sp.factorial(n + k + 1) * sp.factorial(k))


def gk_k64(n, k):
    """K64 gkninit: gkn(n,0) = 1、gkn(n,k) = gkn(n,k-1) (2k-3)(2k+1)/(8 k (k+n+1))"""
    g = sp.Integer(1)
    for kk in range(1, k + 1):
        g *= sp.Rational((2 * kk - 3) * (2 * kk + 1), 8 * kk * (kk + n + 1))
    return g


print('== (a) g_kn の比較 (ヘルプの式と K64 の漸化式) ==')
print('   (-1)!! = 1、(-3)!! = -1 (= (-1)!!/(-1)) を使う')
sp.factorial2  # 負の引数は手で扱う


def dfact(m):
    return {-3: -1, -1: 1}.get(m, None) if m < 0 else sp.factorial2(m)


def gk_help2(n, k):
    return -dfact(2 * k - 1) * dfact(2 * k - 3) * sp.factorial(n + 1) / (8**k * sp.factorial(n + k + 1) * sp.factorial(k))


for n in (0, 1, 2):
    print('   n = %d: k = 0..3  ヘルプ: %s   K64: %s' % (n, [gk_help2(n, k) for k in range(4)], [gk_k64(n, k) for k in range(4)]))

print('\n== (b) (M\') の残差: F = Σ_k g_k (ρ+x)^(1/2-k) z^(m+k) / √ρ の、k = 0..K の和 (n = 1: 4 極) ==')
print('   z の冪が n + k + 1 (K64 の kick の冪 n+k に対応) か、ヘルプの n + k か、を残差の最低次で見る。展開は x, y の冪 (ρ = 1 とする)。')
K = 3
for n in (0, 1, 2):
    for name, gfun in (('K64 の g_kn', gk_k64), ('ヘルプの g_kn', gk_help2)):
        for shift in (1, 0):       # 冪 = n + k + shift
            F = sum(gfun(n, k) * (1 + x)**(sp.Rational(1, 2) - k) * z**(n + k + shift) for k in range(K + 1))
            res = sp.simplify(residual(F).subs(rho, 1))
            t = sp.symbols('t')
            ser = sp.series(sp.expand(res.subs({x: t * x, y: t * y})), t, 0, n + K + shift).removeO()
            ser = sp.simplify(sp.expand(ser))
            # 残差の最低次
            lowest = None
            for p in range(0, n + K + shift + 1):
                c = sp.simplify(sp.expand(ser).coeff(t, p))
                if c != 0:
                    lowest = p
                    break
            print('   n = %d、%s、冪 n+k+%d: 残差の最低次 = %s (打ち切り k ≤ %d の次の項 z^(n+k+%d+1) が埋めるので、n+K+%d-1 = %d 次まで 0 なら正しい)'
                  % (n, name, shift, lowest, K, shift, shift, n + K + shift - 1))

print('\n== (c) K64 の漸化式の閉じた形: g_kn = -(2k+1)!!(2k-3)!!(n+1)!/(8^k (n+k+1)! k!) ((2k-1)!! でなく (2k+1)!!) ==')
def gk_closed(n, k):
    return -dfact(2 * k + 1) * dfact(2 * k - 3) * sp.factorial(n + 1) / (8**k * sp.factorial(n + k + 1) * sp.factorial(k))
ok = all(gk_closed(n, k) == gk_k64(n, k) for n in range(0, 8) for k in range(0, 8))
print('   n = 0..7、k = 0..7 で K64 の漸化式と一致: %s' % ok)
print('   → ヘルプの (2k-1)!! は (2k+1)!! の誤り、z の冪 n+k は n+k+1 の誤り (残差が 0 になるのは K64 の g_kn と冪 n+k+1 の組だけ)。')
print('   直線の極限 ρ → ∞: k = 0 の項だけが残り、A_s → Re[(K_n + i SK_n) z^(n+1)]/(n+1)! (07 章の多極のポテンシャル)。')
