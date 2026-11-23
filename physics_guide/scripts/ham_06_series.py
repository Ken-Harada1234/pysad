# -*- coding: utf-8 -*-
"""展開の途中式を sympy の series で 1 段ずつ確かめる (01a_hamiltonian_derivation.md §9.2〜§9.5、04_quad.md §3.2)。

(a) 基本の展開: sqrt(1 + u)、1/(1 + δ)。
(b) H = -a_s - (1+x/ρ) sqrt((1+δ)^2 - q^2) (q^2 = p_x^2 + p_y^2) の展開を、因子ごとに。次数 (x, p_x, p_y を 1 次と数える) ごとの項。
(c) 偏向磁石: H の項を 1 つずつ、∂H/∂p_x、∂H/∂x を 1 項ずつ。
(d) 四極: 同じく 1 項ずつ。
(e) x'' = p_x'/(1+δ) の δ の展開。
(f) ΔH = p - sqrt(p^2 - q^2) - q^2/(2p) の q^2 の展開 (04_quad.md 式 (4.6))。
"""
import sympy as sp
import _common as C

C.start(__file__)

u = sp.Symbol('u')
d, eps = sp.symbols('delta epsilon', positive=True)
x, y, px, py, rho, K, p = sp.symbols('x y p_x p_y rho K p', positive=True)
q2 = px**2 + py**2

print('== (a) 基本の展開 ==')
print('  sqrt(1 + u) = %s' % sp.series(sp.sqrt(1 + u), u, 0, 4))
print('  sqrt(1 - u) = %s' % sp.series(sp.sqrt(1 - u), u, 0, 4))
print('  1/(1 + δ)   = %s' % sp.series(1 / (1 + d), d, 0, 4))
print('  → u = (p_x^2 + p_y^2)/(1+δ)^2 は p の 2 次なので、u の 1 次までが p の 2 次、u^2 は p の 4 次 (落とす)')

print('\n== (b) H の展開 (因子ごと) ==')
print('  (i) sqrt((1+δ)^2 - q^2) = (1+δ) sqrt(1 - u)、u = q^2/(1+δ)^2:')
sq = sp.sqrt((1 + d)**2 - eps**2 * q2)
ser1 = sp.series(sq, eps, 0, 6).removeO()
for n in (0, 2, 4):
    print('      ε^%d の項: %s' % (n, sp.simplify(ser1.coeff(eps, n))))
print('      (ε は p_x, p_y を ε 倍するための目印。ε^2 までが 2 次、ε^4 = p の 4 次は落とす)')
print('  (ii) (1 + x/ρ) の因子: 1 + ε x/ρ (これ以上の展開はない。x/ρ ≪ 1 の仮定は、掛け算の項を次数で分類するために使う)')
Hk = -(1 + eps * x / rho) * sq
ser2 = sp.series(Hk, eps, 0, 4).removeO()
print('  (iii) 掛け合わせた -(1 + x/ρ) sqrt(...) の、次数 (ε の冪) ごとの項:')
for n in range(4):
    print('      ε^%d: %s' % (n, sp.simplify(ser2.coeff(eps, n))))
print('      ε^0 は定数 (運動に効かない)、ε^1 は x の 1 次、ε^2 は p の 2 次、ε^3 は x と p^2 の積 = 3 次 (落とす)。落とした項の大きさは、2 次の項の x/ρ 倍')

print('\n== (c) 偏向磁石の H2 = x^2/(2ρ^2) - δ x/ρ + q^2/(2(1+δ)) を 1 項ずつ ==')
dd = sp.Symbol('delta', positive=True)
terms = [('x^2/(2ρ^2)', x**2 / (2 * rho**2)), ('-δ x/ρ', -dd * x / rho), ('p_x^2/(2(1+δ))', px**2 / (2 * (1 + dd))), ('p_y^2/(2(1+δ))', py**2 / (2 * (1 + dd)))]
print('  項                  ∂/∂p_x            ∂/∂x              ∂/∂p_y            ∂/∂y')
for name, t in terms:
    print('  %-18s %-17s %-17s %-17s %s' % (name, sp.diff(t, px), sp.diff(t, x), sp.diff(t, py), sp.diff(t, y)))
Htot = sum(t for _, t in terms)
print('  合計: ∂H/∂p_x = %s、∂H/∂x = %s、∂H/∂p_y = %s、∂H/∂y = %s' % (sp.diff(Htot, px), sp.diff(Htot, x), sp.diff(Htot, py), sp.diff(Htot, y)))
print('  → x\' = ∂H/∂p_x = p_x/(1+δ)、p_x\' = -∂H/∂x = -x/ρ^2 + δ/ρ、y\' = p_y/(1+δ)、p_y\' = 0')

print('\n== (d) 四極の H = K (x^2 - y^2)/2 + q^2/(2(1+δ)) を 1 項ずつ ==')
termsq = [('K x^2/2', K * x**2 / 2), ('-K y^2/2', -K * y**2 / 2), ('p_x^2/(2(1+δ))', px**2 / (2 * (1 + dd))), ('p_y^2/(2(1+δ))', py**2 / (2 * (1 + dd)))]
print('  項                  ∂/∂p_x            ∂/∂x              ∂/∂p_y            ∂/∂y')
for name, t in termsq:
    print('  %-18s %-17s %-17s %-17s %s' % (name, sp.diff(t, px), sp.diff(t, x), sp.diff(t, py), sp.diff(t, y)))
print('  → x\' = p_x/(1+δ)、p_x\' = -K x、y\' = p_y/(1+δ)、p_y\' = +K y')

print('\n== (e) x\'\' = p_x\'/(1+δ) の δ の展開 (偏向磁石) ==')
expr = (-x / rho**2 + dd / rho) / (1 + dd)
print('  p_x\'/(1+δ) = (-x/ρ^2 + δ/ρ)/(1+δ) の δ の級数 (2 次まで): %s' % sp.expand(sp.series(expr, dd, 0, 3).removeO()))
print('  δ の 1 次までで x\'\' ≈ -x/ρ^2 + δ/ρ。落とした項: δ x/ρ^2、δ^2/ρ など (x, δ について 2 次)')

print('\n== (f) ΔH = p - sqrt(p^2 - q^2) - q^2/(2p) の q^2 の展開 (04_quad.md 式 (4.6)) ==')
Q = sp.Symbol('Q', positive=True)       # Q = q^2 = p_x^2 + p_y^2
dH = p - sp.sqrt(p**2 - Q) - Q / (2 * p)
print('  sqrt(p^2 - Q) = %s' % sp.series(sp.sqrt(p**2 - Q), Q, 0, 4))
print('  ΔH = %s' % sp.series(dH, Q, 0, 4))
print('  → 最初の項は Q^2/(8 p^3) = (p_x^2+p_y^2)^2/(8 p^3) (横方向の運動量の 4 乗)')
