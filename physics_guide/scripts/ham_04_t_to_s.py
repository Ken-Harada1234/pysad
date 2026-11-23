# -*- coding: utf-8 -*-
"""独立変数 t → s の入れ替え (拡張位相空間、01a_hamiltonian_derivation.md §8) を sympy で確かめる。

座標 (x, y, s)、時間 t。Frenet-Serret 座標のビームのラグランジアン (φ = 0):
  L = -m0 c^2 sqrt(1 - (ẋ^2 + ẏ^2 + h^2 ṡ^2)/c^2) + e (A_x ẋ + A_y ẏ + A_s ṡ)、h = 1 + x/ρ。
(a) 新しいラグランジアン L~ = L t' (独立変数 s、座標 (x, y, t)、q' = dq/ds、t' = dt/ds、s' = 1 なので q̇ = q'/t'、ṡ = 1/t')。
    P_t = ∂L~/∂t' = -H、P_x = ∂L~/∂x' = p_x、P_y = p_y。
(b) H~ = P_x x' + P_y y' + P_t t' - L~ = -p_s。
(c) H~ を (x, y, P_x, P_y, s, P_t) の関数として作り (P_s を H = E で解き、E = -P_t)、そのハミルトンの方程式が
    元の方程式 (dx/dt など) を dṡ で割ったもの (dx/ds = ẋ/ṡ、dP_x/ds = Ṗ_x/ṡ、dt/ds = 1/ṡ) になる。
"""
import sympy as sp
import _common as C

C.start(__file__)

m0, c, e = sp.symbols('m_0 c e', positive=True)
rho = sp.Symbol('rho', positive=True)
x, y, s = sp.symbols('x y s', real=True)
h = 1 + x / rho
Ax, Ay, As = [sp.Function(n)(x, y) for n in ('A_x', 'A_y', 'A_s')]

xd, yd, sd = sp.symbols('xd yd sd', real=True)          # ẋ, ẏ, ṡ (時間微分)
L = -m0 * c**2 * sp.sqrt(1 - (xd**2 + yd**2 + h**2 * sd**2) / c**2) + e * (Ax * xd + Ay * yd + As * sd)
px, py, ps = [sp.diff(L, q) for q in (xd, yd, sd)]
Ham = px * xd + py * yd + ps * sd - L

print("== (a) L~ = L t' と共役運動量 ==")
Tp, xp, yp = sp.symbols('Tp xp yp', real=True)          # t'、x'、y' (s で微分)
Lt = Tp * L.subs({xd: xp / Tp, yd: yp / Tp, sd: 1 / Tp}, simultaneous=True)
Pt = sp.diff(Lt, Tp)
Px = sp.diff(Lt, xp)
Py = sp.diff(Lt, yp)
back = {xp: xd * Tp, yp: yd * Tp}                         # x' = t' ẋ
sub_back = lambda f: f.subs(back, simultaneous=True).subs(Tp, 1 / sd)     # さらに t' ṡ = 1
print('  P_t + H = %s' % sp.simplify(sub_back(Pt) + Ham))
print('  P_x - p_x = %s、P_y - p_y = %s' % (sp.simplify(sub_back(Px) - px), sp.simplify(sub_back(Py) - py)))
print('  → t の共役運動量は -H。x, y の共役運動量は (独立変数を替えただけなので) 元の p_x, p_y のまま')

print("\n== (b) 新しいハミルトニアン H~ = Σ P Q' - L~ ==")
Ht = Px * xp + Py * yp + Pt * Tp - Lt
print('  H~ + p_s = %s' % sp.simplify(sub_back(Ht) + ps))
print("  → H~ = -p_s (s' = 1 かつ t' ṡ = 1 を使った)")

print('\n== (c) H~ を P_s について解いた形と、ハミルトンの方程式 ==')
PX, PY, Ps = sp.symbols('P_X P_Y P_s', real=True)
E = sp.Symbol('E', positive=True)
Pt_s = sp.Symbol('P_t', real=True)
kin = (PX - e * Ax)**2 + (PY - e * Ay)**2
Hfull = c * sp.sqrt(kin + (Ps - e * As)**2 / h**2 + m0**2 * c**2)             # 元の H (φ = 0、Frenet 座標)
Ps_sol = e * As + h * sp.sqrt(E**2 / c**2 - m0**2 * c**2 - kin)                # H = E を P_s について解いたもの
Ht_new = -Ps_sol.subs(E, -Pt_s)                                                  # H~ (P_t = -E)
chk = sp.simplify(sp.expand(Hfull.subs(Ps, Ps_sol) - E))
print('  元の H に P_s の解を代入すると E に戻る: H(P_s = 解) - E = %s' % chk)

ok = True
dHdPs = sp.diff(Hfull, Ps)
pairs = [('dx/ds = ∂H~/∂P_X', sp.diff(Ht_new, PX), sp.diff(Hfull, PX) / dHdPs),
         ('dy/ds = ∂H~/∂P_Y', sp.diff(Ht_new, PY), sp.diff(Hfull, PY) / dHdPs),
         ('dP_X/ds = -∂H~/∂x', -sp.diff(Ht_new, x), -sp.diff(Hfull, x) / dHdPs),
         ('dP_Y/ds = -∂H~/∂y', -sp.diff(Ht_new, y), -sp.diff(Hfull, y) / dHdPs),
         ('dt/ds = ∂H~/∂P_t', sp.diff(Ht_new, Pt_s), 1 / dHdPs)]
for name, new, old in pairs:
    # 元の側は P_s と E を含むので、E = -P_t を通して P_s の解を入れて比べる
    old_s = old.subs(Ps, Ps_sol).subs(E, -Pt_s)
    new_s = new
    d = sp.simplify(sp.expand(new_s - old_s))
    if d != 0:
        # 数値でも確かめる
        vals = {PX: sp.Rational(3, 1000), PY: sp.Rational(-2, 1000), Pt_s: -sp.Rational(5871, 1) * m0 * c, x: sp.Rational(1, 100), y: sp.Rational(-1, 200),
                m0: 1, c: 1, e: 1, rho: 10}
        fa = {Ax: sp.Rational(1, 1000) * x * y, Ay: sp.Rational(2, 1000) * y, As: -sp.Rational(1, 500) * x}
        n1 = new_s.subs(fa).doit().subs(vals)
        n2 = old_s.subs(fa).doit().subs(vals)
        d = sp.N(n1 - n2, 20)
    print('  %-20s の ∂ 側 - 元の方程式/ṡ = %s' % (name, d))
    ok = ok and (d == 0 or abs(d) < 1e-15)
print('  → H~ のハミルトンの方程式は、元の方程式の時間微分を ṡ = ∂H/∂P_s で割ったもの (陰関数の微分)')
print('  P_t = -E は H~ が t を含まない (静的な場) ので保存する' + (' (全部一致)' if ok else ' (不一致あり)'))
