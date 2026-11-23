# -*- coding: utf-8 -*-
r"""16 章の補足の途中式の確認 (§3.3 の 1 次精度、§6.1 の遠方の U、§6.3 の K と z、§7.2 の Schur 補行列、§8.2 の ct = -z)。

(a) シンプレクティック Euler の 1 刻みの局所誤差 (一般の H(q,p) で、厳密な流れとの差が ds^2/2 (H_pq H_p + H_pp H_q))。
(b) 楕円ビームの U の遠方: ∫0^T dt/√((t+A)(t+B)) の閉じた形、U + 2ln(σx+σy) - ln r^2 -> ln2 + γ_E の誤差が 1/r^2 で減ること。
(c) K の生成関数 W(x,y,z) と z のポアソン括弧が 0、D の s* が K で不変。
(d) ブロック逆行列: Λ_aa^{-1} = Σ_aa - Σ_az Σ_zz^{-1} Σ_za、Λ_aa^{-1}Λ_az = -Σ_az Σ_zz^{-1} (sympy)。
(e) ct = -z: 01 章 (1.7) の z = -v (t - t0) で v = c、t0 = 0 (面 S = 0 に設計粒子が来る時刻) とおく。符号を逆にすると x' の tan φc z の符号が逆になる。
"""
import math
import sympy as sp
import mpmath as mp
import _common as C

C.start(__file__)

print('== (a) シンプレクティック Euler の局所誤差 ==')
q, p, ds = sp.symbols('q p ds')
Hf = sp.Function('H')
H = Hf(q, p)
Hq, Hp = sp.diff(H, q), sp.diff(H, p)
# 厳密な流れの ds^2 までの級数: q̇ = H_p, ṗ = -H_q。q̈ = H_pq q̇ + H_pp ṗ
qdd = sp.diff(Hp, q) * Hp + sp.diff(Hp, p) * (-Hq)
pdd = -(sp.diff(Hq, q) * Hp + sp.diff(Hq, p) * (-Hq))
q_ex = q + ds * Hp + ds**2 / 2 * qdd
p_ex = p - ds * Hq + ds**2 / 2 * pdd
# Euler: P = p - ds H_q(q,P)、q' = q + ds H_p(q,P)。P = p + ds P1 + ds^2 P2 と展開して解く
P1, P2, dP = sp.symbols('P1 P2 dP')
Pser = p + ds * P1 + ds**2 * P2
# H_q(q, P) を P = p + dP で dP の 1 次まで展開 (dP = O(ds))
Hq_P = Hq + sp.diff(Hq, p) * (Pser - p)
eq = sp.expand(Pser - (p - ds * Hq_P))
sol1 = sp.solve(eq.coeff(ds, 1), P1)[0]
sol2 = sp.solve(eq.coeff(ds, 2).subs(P1, sol1), P2)[0]
P_eu = p + ds * sol1 + ds**2 * sol2
Hp_P = Hp + sp.diff(Hp, p) * (P_eu - p)
q_eu = q + ds * Hp_P
dq = sp.simplify(sp.expand(q_ex - q_eu).coeff(ds, 2))
dp_ = sp.simplify(sp.expand(p_ex - P_eu).coeff(ds, 2))
print('P_eu = p - ds H_q + ds^2 (H_qp H_q) + ...: ds^2 の係数 sol2 =', sp.simplify(sol2))
print('q: (厳密 - Euler) の ds^2 の係数 =', dq)
print('   期待値 1/2 (H_pq H_p + H_pp H_q) との差 =',
      sp.simplify(dq - sp.Rational(1, 2) * (sp.diff(Hp, q) * Hp + sp.diff(Hp, p) * Hq)))
print('p: (厳密 - Euler) の ds^2 の係数 =', dp_)
print('ds^1 の係数は両者とも H_p, -H_q で一致:',
      sp.simplify(sp.expand(q_ex - q_eu).coeff(ds, 1)), sp.simplify(sp.expand(p_ex - P_eu).coeff(ds, 1)))

print('\n== (b) 楕円ビームの U の遠方 ==')
A, B, T = sp.symbols('A B T', positive=True)
tt = sp.symbols('t', positive=True)
F = 2 * sp.log((sp.sqrt(tt + A) + sp.sqrt(tt + B)) / (sp.sqrt(A) + sp.sqrt(B)))
print('d/dT で被積分関数に戻るか:', sp.simplify(sp.diff(F, tt) - 1 / sp.sqrt((tt + A) * (tt + B))),
      ' F(0) =', F.subs(tt, 0))
print('T->∞: F - [ln(4T) - 2 ln(√A+√B)] ->', sp.limit(F - (sp.log(4 * tt) - 2 * sp.log(sp.sqrt(A) + sp.sqrt(B))), tt, sp.oo))
mp.mp.dps = 30
def U_num(x, y, sx, sy):
    a0, b0 = 2 * sx**2, 2 * sy**2
    f = lambda t: (1 - mp.e**(-x**2 / (a0 + t) - y**2 / (b0 + t))) / mp.sqrt((a0 + t) * (b0 + t))
    return mp.quad(f, [0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e8, mp.inf])
c0 = mp.log(2) + mp.euler
print('ln2 + γ_E =', mp.nstr(c0, 10))
for (sx, sy) in [(1, 1), (3, 0.5), (1, 4)]:
    row = []
    for r in (10, 30, 100, 300):
        x, y = r * 0.6, r * 0.8
        d = U_num(x, y, sx, sy) + 2 * mp.log(sx + sy) - mp.log(r**2) - c0
        row.append('%.2e (r^2 d = %.2f)' % (float(d), float(d * r**2)))
    print('σ = (%g, %g): 誤差 U + 2ln(σx+σy) - ln r^2 - (ln2+γ_E) at r = 10, 30, 100, 300:' % (sx, sy), row)

print('\n== (c) K と z のポアソン括弧 ==')
x, y, z, px, py, dl = sp.symbols('x y z p_x p_y delta')
W = sp.Function('W')(x, y, z)
def pb(f, g):   # [f,g] = Σ (∂q f ∂p g - ∂p f ∂q g)、(x,px)(y,py)(z,δ)
    return sum(sp.diff(f, qq) * sp.diff(g, pp) - sp.diff(f, pp) * sp.diff(g, qq)
               for qq, pp in ((x, px), (y, py), (z, dl)))
print('[W(x,y,z), z] =', sp.simplify(pb(W, z)), ' (K は z を変えない)')
print('[z, F] with F = s*(z) (px²+py²)/2:', sp.simplify(pb(z, (z - sp.Symbol('zs')) / 2 * (px**2 + py**2) / 2)), ' (D も z を変えない)')

print('\n== (d) ブロック逆行列 ==')
s = sp.symbols('s11 s12 s13 s22 s23 s33')
Sg = sp.Matrix([[s[0], s[1], s[2]], [s[1], s[3], s[4]], [s[2], s[4], s[5]]])
Lam = Sg.inv()
Saa, Saz, Szz = Sg[:2, :2], Sg[:2, 2:], Sg[2, 2]
Schur = Saa - Saz * Saz.T / Szz
Laa, Laz = Lam[:2, :2], Lam[:2, 2:]
print('Λ_aa^{-1} - (Σ_aa - Σ_az Σ_zz^{-1} Σ_za) =', sp.simplify(Laa.inv() - Schur).tolist())
print('Λ_aa^{-1} Λ_az + Σ_az Σ_zz^{-1} =', sp.simplify(Laa.inv() * Laz + Saz / Szz).tolist())
print('Λ_zz - 1/Σ_zz は Σ_za Schur^{-1} Σ_az / Σ_zz^2 に等しい:',
      sp.simplify(Lam[2, 2] - 1 / Szz - (Saz.T * Schur.inv() * Saz)[0] / Szz**2))

print('\n== (e) ct = -z ==')
zz, xx, c_, s_ = sp.symbols('z x c s')
t0 = 0
ct_event = -(zz)      # z = -c (t - t0), t0 = 0
for sign, label in ((1, 'ct = -z (導出)'), (-1, 'ct = +z (逆)')):
    ct = sign * (-zz)
    x1 = c_ * xx
    x2 = (x1 - s_ / c_ * 0) / 1   # placeholder
    gam, beta = 1 / c_, s_
    x2 = gam * (x1 - beta * ct)
    print(label, ': boost 後の x_2 =', sp.simplify(x2.subs(s_, sp.sqrt(1 - c_**2))), '  (= x + tanφc z なら z の係数 s/c)')
