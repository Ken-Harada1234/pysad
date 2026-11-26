# -*- coding: utf-8 -*-
r"""02 章の途中式の確認 (本文に書き足した展開)。

(a) §3 (2.2): E/(P0 v0) = h0 h1/p0^2 (P0 = p0 m c、β0 = p0/h0)。
(b) §6: pxy2dpz = sqrt(1-x) - 1 の級数 (K64 は x の 3 次まで)、打ち切りの誤差と閾値 8e-6、有理形 -x/(1 + sqrt(1-x))。
(c) §7.1: 小さい傾きの展開 (u = q^2 + r^2): 1/ζ の二項係数、1/ζ - 1、D の z' = 1 - 1/ζ + u/2。
"""
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)
R = sp.Rational

print('(a) (2.2)')
p0, m, c = sp.symbols('p0 m c', positive=True)
dl = sp.symbols('delta', positive=True)
h0 = sp.sqrt(1 + p0**2)
h1 = sp.sqrt(1 + p0**2 * (1 + dl)**2)
beta0 = p0 / h0
P0 = p0 * m * c
v0 = beta0 * c
E = h1 * m * c**2
print('  P0 v0 = p0 mc · (p0/h0) c = (p0^2/h0) m c^2:', sp.simplify(P0 * v0 - p0**2 / h0 * m * c**2) == 0)
print('  E/(P0 v0) - h0 h1/p0^2 =', sp.simplify(E / (P0 * v0) - h0 * h1 / p0**2))

print('\n(b) sqrt(1-x) - 1 の展開と pxy2dpz')
x = sp.symbols('x')
ser = sp.series(sp.sqrt(1 - x) - 1, x, 0, 6).removeO()
print('  sqrt(1-x) - 1 =', ser, '+ O(x^6)')
print('  二項係数 C(1/2,k)(-1)^k (k = 1..5):', [sp.binomial(R(1, 2), k) * (-1)**k for k in range(1, 6)])
print('  K64 の級数 -x(1/2 + x/8 + x^2/16) との差 = 次の項:', sp.expand(ser - (-x * (R(1, 2) + x * (R(1, 8) + x * R(1, 16))))))
mp.mp.dps = 40
for xv in (mp.mpf('8e-6'), mp.mpf('1e-6'), mp.mpf('1e-3')):
    exact = mp.sqrt(1 - xv) - 1
    ser3 = -xv * (mp.mpf(1) / 2 + xv * (mp.mpf(1) / 8 + xv * mp.mpf(1) / 16))
    rat = -xv / (1 + mp.sqrt(1 - xv))
    print('  x = %s: 級数 3 項の絶対誤差 %.2e (相対 %.2e、5x^3/64 = %.2e)、有理形の絶対誤差 %.1e' % (
        mp.nstr(xv, 3), float(abs(ser3 - exact)), float(abs((ser3 - exact) / exact)), float(5 * xv**3 / 64), float(abs(rat - exact))))
print('  倍精度の丸め 2^-53 = %.2e。x = 8e-6 の相対誤差は丸めより小さい。' % 2.0**-53)
print('  直接式 sqrt(1-x) - 1 を倍精度で計算した相対誤差の目安 (x = 8e-6): 丸め/x = %.1e' % (2.0**-53 / 8e-6))
print('  有理形: sqrt(1-x) - 1 = ((1-x) - 1)/(sqrt(1-x) + 1) の差 =', sp.simplify(sp.sqrt(1 - x) - 1 - (-x / (1 + sp.sqrt(1 - x)))))

print('\n(c) 小さい傾きの展開 (u = q^2 + r^2、ζ = sqrt(1-u))')
u = sp.symbols('u')
print('  1/ζ = (1-u)^(-1/2) の二項係数 C(-1/2,k)(-1)^k:', [sp.binomial(-R(1, 2), k) * (-1)**k for k in range(0, 5)])
print('  1/ζ =', sp.series(1 / sp.sqrt(1 - u), u, 0, 5))
z = sp.sqrt(1 - u)
zprime = 1 - 1 / z + (1 - z**2) / 2
print('  (1 - ζ^2)/2 = u/2:', sp.simplify((1 - z**2) / 2 - u / 2) == 0)
print('  z\' = 1 - 1/ζ + u/2 =', sp.series(zprime, u, 0, 5))
print('  -(2+ζ)(1-ζ)^2/(2ζ) との差:', sp.simplify(zprime + (2 + z) * (1 - z)**2 / (2 * z)))
