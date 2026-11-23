# -*- coding: utf-8 -*-
"""ヘルプ (SADHelp.html、K64 版同梱) の exp(:H:) の記法の確認 (02_drift.md 付録 A、C_poisson_lie.md §8.2)。

ヘルプは、DRIFT・QUAD・SEXT 等の transformation を exp(:H:) の積で書くが、ポアソン括弧の定義も、
:H: の符号も書いていない。ここでは、ヘルプの式そのものを sympy に入れ、
 (1) exp(:H:) を「H (長さをかけた量) を時間 1 のハミルトニアンとする流れ」 (x → x + ∂H/∂p、p → p - ∂H/∂x) と読むと、
     ヘルプの各式 (DRIFT の √ の項、QUAD の K1 の項、SEXT の H_2、BEND の線形フリンジの V) が物理的に正しい変換になること、
     逆向き (標準の括弧 [f,g] = ∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q で exp(:H:) g = g + [H,g]) だと全部逆になること、
 (2) DRIFT の H の E/v0 の符号だけは、この読みでも他の式 (QUAD の H_2n、dynamics の H、BEND の body の z) と合わないこと、
を示す。
"""
import math
import sympy as sp
import _common as C

C.start(__file__)

x, y, px, py, z, p = sp.symbols('x y p_x p_y z p', real=True)
K1, K2, L, A, B = sp.symbols('K1 K2 L A B', real=True)
n = sp.Symbol('n', positive=True)
P0, H0 = sp.symbols('p0 h0', positive=True)
h1 = sp.sqrt(1 + (P0 * p)**2)
E_over_p0v0 = H0 * h1 / P0**2         # E/(p0 v0) (E = h1 mc^2、v0 = (p0/h0) c)  [02 章 (2.2)]
vv0 = p * H0 / h1                     # v/v0 = ∂[E/(p0v0)]/∂p

pz = sp.sqrt(p**2 - px**2 - py**2)
q = (x, px, y, py, z, p)


def flow_forward(Hh):
    """H を時間 1 のハミルトニアンとする無限小の変化 (Δq = ∂H/∂p、Δp = -∂H/∂q)。z と共役なのは p (z' = ∂H/∂p)"""
    return {'x': sp.diff(Hh, px), 'px': -sp.diff(Hh, x), 'y': sp.diff(Hh, py), 'py': -sp.diff(Hh, y),
            'z': sp.diff(Hh, p)}


def bracket_std(f, g):
    """標準の括弧 [f,g] = Σ (∂f/∂q ∂g/∂p - ∂f/∂p ∂g/∂q)。(x,px)、(y,py)、(z,p) が共役"""
    r = 0
    for qq, pp in ((x, px), (y, py), (z, p)):
        r += sp.diff(f, qq) * sp.diff(g, pp) - sp.diff(f, pp) * sp.diff(g, qq)
    return r


vals = {px: 0.013, py: -0.021, p: 1.004, x: 1e-3, y: -2e-3, P0: 3.0e9 / 0.51099895e6}
vals[H0] = math.sqrt(1 + float(vals[P0])**2)
vals[L] = 1.0
vals[K1] = 0.7; vals[K2] = 0.9; vals[n] = 1.0

print('== (1) ヘルプの各式を、exp(:H:) = 「H を時間 1 のハミルトニアンとする流れ」と読む ==')
print('  標準の括弧での :H: g = [H, g] は、流れの逆向き (Δx = [H,x] = -∂H/∂px) になる。確認 (H = px^2/2): [H,x] = %s、∂H/∂px = %s'
      % (bracket_std(px**2 / 2, x), sp.diff(px**2 / 2, px)))

print('\n -- DRIFT: H = (-√(p²-px²-py²) + 1 - E/v0) L  (ヘルプ原文。平文・TeX とも)')
Hdrift_help = (-pz + 1 - E_over_p0v0) * L
f = flow_forward(Hdrift_help)
print('   Δx = ∂H/∂px = %s  → x は px/pz × L だけ進む (正しい)' % sp.simplify(f['x']))
print('   Δpx = -∂H/∂x = %s' % f['px'])
zz = sp.simplify(f['z'])
print('   Δz = ∂H/∂p = %s' % zz)
Hdrift_plus = (-pz + E_over_p0v0) * L
zz_plus = sp.simplify(flow_forward(Hdrift_plus)['z'])
print('   E/v0 の符号を + に直した H = (-pz + E/(p0v0)) L では Δz = %s' % zz_plus)
claim_body = (-p / pz + vv0) * L      # BEND の body の z の式で ANGLE → 0 の極限: z2 = z1 - s p/√(p²-py²) + (v1/v0) L'、s → L px/p... は py=0 での直線の経路長
# 直線のドリフトの z は、経路長 L p/pz に対し  z' = -p/pz + v/v0  (01 章 (1.9) の dv = 1 - v/v0 と z' = -dv - (p/pz - 1))
print('   物理的な Δz = (-p/pz + v/v0) L (= -(p/pz - 1) L - dv L):')
print('     + の符号の H との差 (数値) = %.2e' % abs(float((zz_plus - claim_body).subs(vals))))
print('     ヘルプの H (- の符号) との差 (数値) = %.2e  (= 2 (v/v0) L = %.6f)' % (
    abs(float((zz - claim_body).subs(vals))), float((2 * vv0 * L).subs(vals))))

print('\n -- QUAD: H_2n = ((-p + (px²+py²)/(2p) + E/v0) L + K1 (x²-y²)/2)/n  (ヘルプ原文) ----')

H2n = ((-p + (px**2 + py**2) / (2 * p) + E_over_p0v0) * L + K1 * (x**2 - y**2) / 2) / n
f = flow_forward(H2n)
print('   Δpx = -∂H/∂x = %s (K1 > 0 で x 方向に収束: 力は -K1 x)' % f['px'])
print('   Δpy = %s (垂直には発散)' % f['py'])
print('   Δz = ∂H/∂p = %s' % sp.simplify(f['z']))
chk = sp.simplify(f['z'] * n / L - (-(px**2 + py**2) / (2 * p**2) - 1 + vv0))
print('   z\' = -1 - (px²+py²)/(2p²) + v/v0 (= -(q²/2) - dv) との差: %s  (E/v0 の符号は + で、物理的に正しい)' % chk)
print('   (DRIFT の H の符号は、この H_2n と逆。同じ表記 E/v0 なので、DRIFT の方が誤植とみられる。)')

print('\n -- SEXT: H_2 = K2/3! Re(x - i y)^3 (薄いキック) ----')
H2 = K2 / 6 * sp.re(sp.expand((x - sp.I * y)**3))
f = flow_forward(sp.simplify(H2))
print('   Δpx = -∂H/∂x = %s' % sp.simplify(f['px']))
print('   Δpy = -∂H/∂y = %s' % sp.simplify(f['py']))
print('   → Δpx = -(K2/2)(x²-y²)、Δpy = K2 x y: 6 極の標準のキック。K2 > 0 は「水平収束」の符号 (ヘルプの K2 の説明) と整合')
print('   標準の括弧で exp(:H:) g = g + [H,g] と読むと符号が逆: [H, px] = %s' % sp.simplify(bracket_std(H2, px)))

print('\n -- BEND の線形フリンジ: V = -(1/p)(F1²/(24ρb) px - F1/(12ρb²) y² ... ) (平文版) ----')
F1, rb = sp.symbols('F1 rho_b', positive=True)
V_ascii = -F1**2 / rb * px / p / 24 - F1 / rb**2 * y**2 / p / 12 + 1 / rb**2 / F1 * y**4 / p / 6
f = flow_forward(V_ascii)
dxfr = F1**2 / (24 * rb)
dyfr = F1 / (6 * rb**2)
dyfra = sp.Rational(2, 3) / (F1 * rb**2)
print('   Δx = %s  (ヘルプの x2 = x1 + dxfr (p-p0)/p の p 依存部分 -dxfr/p と一致)' % sp.simplify(f['x']))
print('   Δpy = %s' % sp.simplify(f['py']))
print('   ヘルプの py2 = py1 + y (dyfr - dyfra y²)/p1² の y(dyfr - dyfra y²)/p との差 (q の単位 py/p なら 1/p²): %s' %
      sp.simplify(f['py'] - (dyfr - dyfra * y**2) * y / p))
print('   → 平文版の V は、ヘルプ自身の線形フリンジの式 (x2、py2) と、「流れ」の読みで一致する。')

print('\n== (2) 結論 ==')
print('  ヘルプの exp(:H:) は、H (長さをかけた量) を時間 1 のハミルトニアンとする流れ。')
print('  C_poisson_lie.md の記法 (標準の括弧、転送は exp(-s:H_ham:)) では、ヘルプの :f: は -:f: に当たる。')
print('  ヘルプに括弧の定義は無く、この対応は式どうしの整合から確定した (DRIFT の E/v0 の符号だけは例外で誤植とみられる)。')
