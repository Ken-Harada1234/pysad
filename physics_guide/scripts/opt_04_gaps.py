# -*- coding: utf-8 -*-
"""12 章 (12_linear_optics.md) の導出の途中式の確認 (本文の補足の式 (12.11) 以降と、§3.2・§5・§6・§7・§8 の途中式)。

(A) BAB^{-1} = P(μ) の行列の積を成分ごとに (§6.1)。
(B) A2 = T A1 T^{-1} の (1,1)、(2,2)、(1,2) 成分と α2 の導出 (§6.2)。
(C) 2×2 の恒等式 r J r^T = det r J、(F1)〜(F5) (§5.2・§5.3)。
(D) 厚い四極の色収差: ∫β ds の閉じた式、端の平均との差の展開、例のリングの数値 (§7)。
(E) 閉軌道の Newton 法の収束係数 C (r_{k+1} ≈ C r_k²) (§3.2)。
(F) (12.10) の H のシンプレクティック性: 手で示すブロックの中間式と、6 つの全ブロックの確認 (§8)。
(G) qmdiag.f の特別な場合の r1..r4 の式が X_f N^{-1}、Y_f (...) であること (§5.3)。
"""
import io
import math
import os
import re
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string

print('pysad:', C.pysad_src())

# ------------------------------------------------------------------ (A)
print('\n== (A) B A B^{-1} = P(μ): 行列の積 ==')
al, be, mu = sp.symbols('alpha beta mu', positive=True)
ga = (1 + al**2) / be
A = sp.Matrix([[sp.cos(mu) + al * sp.sin(mu), be * sp.sin(mu)], [-ga * sp.sin(mu), sp.cos(mu) - al * sp.sin(mu)]])
B = sp.Matrix([[1 / sp.sqrt(be), 0], [al / sp.sqrt(be), sp.sqrt(be)]])
Binv = sp.Matrix([[sp.sqrt(be), 0], [-al / sp.sqrt(be), 1 / sp.sqrt(be)]])
print('B B^{-1} - 1 =', sp.simplify(B * Binv - sp.eye(2)).tolist())
BA = (B * A).applyfunc(sp.simplify)
BA_hand = sp.Matrix([[(sp.cos(mu) + al * sp.sin(mu)) / sp.sqrt(be), sp.sqrt(be) * sp.sin(mu)],
                     [(al * sp.cos(mu) - sp.sin(mu)) / sp.sqrt(be), sp.sqrt(be) * sp.cos(mu)]])
print('B A - (本文の成分) =', sp.simplify(BA - BA_hand).tolist())
P = sp.Matrix([[sp.cos(mu), sp.sin(mu)], [-sp.sin(mu), sp.cos(mu)]])
print('(B A) B^{-1} - P(μ) =', sp.simplify(BA_hand * Binv - P).tolist())

# ------------------------------------------------------------------ (B)
print('\n== (B) A2 = T A1 T^{-1} の成分と α2 ==')
t11, t12, t21, t22, a11, a12, a21, a22 = sp.symbols('t11 t12 t21 t22 a11 a12 a21 a22')
T = sp.Matrix([[t11, t12], [t21, t22]])
A1 = sp.Matrix([[a11, a12], [a21, a22]])
Tinv = sp.Matrix([[t22, -t12], [-t21, t11]])      # det T = 1 のときの T^{-1}
A2 = (T * A1 * Tinv).applyfunc(sp.expand)
hand11 = t22 * (t11 * a11 + t12 * a21) - t21 * (t11 * a12 + t12 * a22)
hand22 = -t12 * (t21 * a11 + t22 * a21) + t11 * (t21 * a12 + t22 * a22)
hand12 = -t22 * 0 + (-t12) * (t11 * a11 + t12 * a21) + t11 * (t11 * a12 + t12 * a22)
print('(A2)_11 - 手計算 =', sp.expand(A2[0, 0] - hand11), '、(A2)_22 - 手計算 =', sp.expand(A2[1, 1] - hand22),
      '、(A2)_12 - 手計算 =', sp.expand(A2[0, 1] - hand12))
print('(A2)_12 - 本文の式 t11² a12 - t11 t12 (a11 - a22) - t12² a21 =',
      sp.expand(A2[0, 1] - (t11**2 * a12 - t11 * t12 * (a11 - a22) - t12**2 * a21)))
diff = sp.expand(hand11 - hand22)
print('(A2)_11 - (A2)_22 を a_ij で整理 =', sp.collect(diff, [a11, a12, a21, a22]))
s_ = sp.symbols('s')                                   # sin μ
b1, a1, g1 = sp.symbols('beta1 alpha1 gamma1')
diff_book = (t11 * t22 + t12 * t21) * (a11 - a22) - 2 * t11 * t21 * a12 + 2 * t12 * t22 * a21
print('  本文の式 (t11 t22 + t12 t21)(a11 - a22) - 2 t11 t21 a12 + 2 t12 t22 a21 との差 =', sp.expand(diff - diff_book))
sub = {a12: b1 * s_, a21: -g1 * s_}
d2 = sp.expand(diff_book.subs(a11, a22 + 2 * a1 * s_).subs(sub))
print('a11 - a22 = 2 α1 s、a12 = β1 s、a21 = -γ1 s を入れる: (A2)_11 - (A2)_22 =', sp.collect(d2, s_))
alpha2 = sp.expand(d2 / (2 * s_))
alpha2_book = -t11 * t21 * b1 + (t11 * t22 + t12 * t21) * a1 - t12 * t22 * g1
print('α2 = 差 / (2 s) =', alpha2, '、(12.9) の式との差:', sp.expand(alpha2 - alpha2_book))
beta2 = sp.expand(A2[0, 1].subs(a11, a22 + 2 * a1 * s_).subs(sub) / s_)
print('β2 = (A2)_12 / s =', beta2, '、(12.9) の式との差:', sp.expand(beta2 - (t11**2 * b1 - 2 * t11 * t12 * a1 + t12**2 * g1)))
# (A2)_11 + (A2)_22 = tr A1 (det T = 1 を使う)
tr_diff = sp.expand(A2[0, 0] + A2[1, 1] - a11 - a22)
print('(A2)_11 + (A2)_22 - (a11 + a22) = (t11 t22 - t12 t21 - 1)(a11 + a22) の形: ',
      sp.expand(tr_diff - (t11 * t22 - t12 * t21 - 1) * (a11 + a22)))

# ------------------------------------------------------------------ (C)
print('\n== (C) 2×2 の恒等式 ==')
a, b, c, d = sp.symbols('a b c d')
Ar = sp.Matrix([[a, b], [c, d]])
J = sp.Matrix([[0, 1], [-1, 0]])
print('A J =', (Ar * J).tolist(), '、A J A^T =', sp.expand(Ar * J * Ar.T).tolist())
print('(A J A^T)_12 = a(a c... ) の成分:', sp.expand((Ar * J * Ar.T)[0, 1]), '、det A J の (1,2) =', Ar.det())
print('A J A^T - det A J =', sp.expand(Ar * J * Ar.T - Ar.det() * J).tolist())


def bar(X):
    return J * X.T * J


e, f, g, h = sp.symbols('e f g h')
Br = sp.Matrix([[e, f], [g, h]])
print('bar A =', bar(Ar).tolist(), '、-adj A =', (-Ar.adjugate()).tolist())
print('(F1) A barA + det A I =', sp.expand(Ar * bar(Ar) + Ar.det() * sp.eye(2)).tolist())
print('(F2) bar(AB) + bar(B) bar(A) =', sp.expand(bar(Ar * Br) + bar(Br) * bar(Ar)).tolist())
print('  途中: J (AB)^T J - (J B^T J)(J^{-1} A^T J) =',
      sp.expand(J * (Ar * Br).T * J - (J * Br.T * J) * (J.inv() * Ar.T * J)).tolist())
print('  J^{-1} + J =', (J.inv() + J).tolist(), '(J^{-1} = -J = J^T)')
print('(F3) bar(bar A) - A =', sp.expand(bar(bar(Ar)) - Ar).tolist())
print('(F4) tr bar A + tr A =', sp.expand(bar(Ar).trace() + Ar.trace()))
print('(F5) A - bar A - tr A I =', sp.expand(Ar - bar(Ar) - Ar.trace() * sp.eye(2)).tolist())

# ------------------------------------------------------------------ (D)
print('\n== (D) 厚い四極の色収差 ==')
K, L, b0, a0 = sp.symbols('K L beta0 alpha0', positive=True)
g0 = (1 + a0**2) / b0
s = sp.symbols('s', positive=True)
w = sp.sqrt(K)
Cc, Ss = sp.cos(w * s), sp.sin(w * s) / w
beta_s = Cc**2 * b0 - 2 * Cc * Ss * a0 + Ss**2 * g0
al_s = -sp.diff(beta_s, s) / 2
gam_s = (1 + al_s**2) / beta_s
print("β(s) = C² β0 - 2 C S α0 + S² γ0 (C = cos √K s、S = sin(√K s)/√K)、α = -β'/2、γ = (1+α²)/β")
print("  β'' - (2γ - 2Kβ) =", sp.simplify(sp.diff(beta_s, s, 2) - (2 * gam_s - 2 * K * beta_s)))
print("  α' - (Kβ - γ) =", sp.simplify(sp.diff(al_s, s) - (K * beta_s - gam_s)))
print('  β(0) = β0、β\'(0) = -2α0:', beta_s.subs(s, 0), sp.simplify(sp.diff(beta_s, s).subs(s, 0)))
I_closed = (b0 * (L / 2 + sp.sin(2 * w * L) / (4 * w)) - a0 * sp.sin(w * L)**2 / K
            + g0 / K * (L / 2 - sp.sin(2 * w * L) / (4 * w)))
# 閉じた式の微分が β(L) に等しいことで確認 (d/dL ∫0^L β ds = β(L))
print("  d/dL (閉じた式) - β(L) =", sp.simplify(sp.diff(I_closed, L) - beta_s.subs(s, L)), '、L = 0 で 0:', sp.limit(I_closed, L, 0))
beta_end = beta_s.subs(s, L)
bar_beta = I_closed / L
avg_end = (b0 + beta_end) / 2
ser_b = sp.series(bar_beta, L, 0, 4).removeO()
ser_e = sp.series(avg_end, L, 0, 4).removeO()
print('β̄ = ∫β ds / L の L 展開:', sp.collect(sp.expand(ser_b), L))
print('(β0 + β(L))/2 の L 展開:', sp.collect(sp.expand(ser_e), L))
dser = sp.expand(sp.series(bar_beta - avg_end, L, 0, 4).removeO())
print('β̄ - (β0 + β(L))/2 =', sp.collect(dser, L))
print('  L² の係数 + (γ0 - Kβ0)/6 =', sp.simplify(dser.coeff(L, 2) + (g0 - K * b0) / 6))
print('  L³ の係数 =', sp.simplify(dser.coeff(L, 3)), "(= -β'''/24、β''' = 8Kα0 (β' = -2α、γ' = 2Kα、β''' = 2γ' - 2Kβ'))")


def beta_mean(b0_, a0_, Kk, L_, n=4000):
    g0_ = (1 + a0_**2) / b0_
    ss = np.linspace(0, L_, n + 1)
    if Kk > 0:
        wv = math.sqrt(Kk); CC = np.cos(wv * ss); SS = np.sin(wv * ss) / wv
    else:
        wv = math.sqrt(-Kk); CC = np.cosh(wv * ss); SS = np.sinh(wv * ss) / wv
    bb = CC**2 * b0_ - 2 * CC * SS * a0_ + SS**2 * g0_
    return np.trapezoid(bb, ss) / L_


def beta_mean_closed(b0_, a0_, Kk, L_):
    """∫β ds / L の閉じた式 (K > 0)。K < 0 は φ → iφ の双曲線関数。"""
    g0_ = (1 + a0_**2) / b0_
    if Kk > 0:
        wv = math.sqrt(Kk); ph = wv * L_
        return (b0_ * (0.5 + math.sin(2 * ph) / (4 * ph)) - a0_ * L_ * (math.sin(ph) / ph)**2
                + g0_ * L_**2 / ph**2 * (0.5 - math.sin(2 * ph) / (4 * ph)))
    wv = math.sqrt(-Kk); ph = wv * L_
    return (b0_ * (0.5 + math.sinh(2 * ph) / (4 * ph)) - a0_ * L_ * (math.sinh(ph) / ph)**2
            + g0_ * L_**2 / ph**2 * (math.sinh(2 * ph) / (4 * ph) - 0.5))


srcq = '''MOMENTUM = 3GEV;
DRIFT D1=(L=0.8) D2=(L=0.3) B=(L=1.2);
QUAD QF=(L=0.3 K1=0.35) QD=(L=0.3 K1=-0.33);
MARK ST=();
LINE CELL=(QF D1 B D2 QD D1 B D2);
LINE RING=(ST CELL CELL CELL CELL CELL CELL CELL CELL);
FFS;
USE RING;
'''
ffq, _ = run_sad_string(srcq, out=io.StringIO())
elq = ffq.lat.elements
twq = ffq.lat.twiss()
print('\n例のリング (四極だけ、偏向なし、opt_01 (e2) と同じ)。K = pysad の el.K1 (= SAD の K1 / L) [1/m²]、y 面は K の符号を反転')
print('  面 種類   K       φ²=KL²   β端    α端      β端平均   β̄(厚い,数値積分)  β̄(閉じた式)  比 β̄/β端平均   1+φ²/6')
tot = {'x': [0.0, 0.0], 'y': [0.0, 0.0]}
seen = set()
for k, el in enumerate(elq):
    if type(el).__name__ != 'Quadrupole':
        continue
    Kk, Ll = el.K1, el.L
    for pl, sg in (('x', 1.0), ('y', -1.0)):
        b0_, a0_, b1_ = twq[k]['b' + pl], twq[k]['a' + pl], twq[k + 1]['b' + pl]
        Kp = sg * Kk
        bm = beta_mean(b0_, a0_, Kp, Ll)
        bc = beta_mean_closed(b0_, a0_, Kp, Ll)
        tot[pl][0] += Kp * Ll * (b0_ + b1_) / 2
        tot[pl][1] += Kp * Ll * bm
        key = (pl, round(Kk, 6))
        if key not in seen:
            seen.add(key)
            print('  %s  %s  %7.4f  %7.4f  %7.4f  %7.4f  %9.5f  %9.5f  %9.5f  %9.5f  %9.5f' % (
                pl, 'QF' if Kk > 0 else 'QD', Kp, Kp * Ll**2, b0_, a0_, (b0_ + b1_) / 2, bm, bc, bm / ((b0_ + b1_) / 2),
                1 + Kp * Ll**2 / 6 * (1 - ((1 + a0_**2) / b0_) / (Kp * b0_))))
for pl in 'xy':
    print('  %s 面: -(1/4π) Σ K L β: 端の平均 = %.4f、厚い四極 (∫β ds / L) = %.4f、比 (端平均 / 厚い) = %.4f' % (
        pl, -tot[pl][0] / (4 * math.pi), -tot[pl][1] / (4 * math.pi), tot[pl][0] / tot[pl][1]))
print('  (追跡の Jacobian: dν/dδ = (-1.2144, -1.1799) (opt_01 (e2))。厚い四極の ∫ と一致)')

# ------------------------------------------------------------------ (E)
print('\n== (E) 閉軌道の Newton 法の収束係数 ==')
here = os.path.dirname(os.path.abspath(__file__))
txt = open(os.path.join(here, '..', 'outputs', 'opt_01_linear_optics.out.txt'), encoding='utf-8').read()
m = re.search(r'反復ごとの残差 max\|f\(c\) - c\|:\s*(.*)', txt)
rs = [float(x) for x in m.group(1).split()]
print('opt_01 (a) の残差 r_k =', ' '.join('%.2g' % x for x in rs[:5]))
for k in range(3):
    print('  k = %d -> %d: r_{k+1}/r_k² = %.0f、r_k = %.2g の 2 乗 = %.2g、× 250 = %.2g (実測 r_{k+1} = %.2g)' % (
        k + 1, k + 2, rs[k + 1] / rs[k]**2, rs[k], rs[k]**2, 250 * rs[k]**2, rs[k + 1]))
print('  k = 4 -> 5: r_4² × 250 = %.1g だが実測 r_5 = %.1g (丸め誤差の床。2 乗則の予想より上)' % (250 * rs[3]**2, rs[4]))

# ------------------------------------------------------------------ (F)
print('\n== (F) (12.10) の H のシンプレクティック性 ==')
x11, x12, x21, x22, y11, y12, y21, y22, aa_ = sp.symbols('x11 x12 x21 x22 y11 y12 y21 y22 a')
Hx = sp.Matrix([[x11, x12], [x21, x22]])
Hy = sp.Matrix([[y11, y12], [y21, y22]])
dx, dy = Hx.det(), Hy.det()
u = 1 + aa_
I2 = sp.eye(2)
H = sp.BlockMatrix([[(1 - dx / u) * I2, Hx * J * Hy.T * J / u, -Hx],
                    [Hy * J * Hx.T * J / u, (1 - dy / u) * I2, -Hy],
                    [-J * Hx.T * J, -J * Hy.T * J, aa_ * I2]]).as_explicit()
print('X^T J X - det X J (X = Hx) =', sp.expand(Hx.T * J * Hx - dx * J).tolist(), '、det(-J Hx^T J) - det Hx =', sp.expand((J * Hx.T * J).det() - dx))
print('(3,3): (-Hx)^T J (-Hx) + (-Hy)^T J (-Hy) + a² J - J = (dx + dy + a² - 1) J ->',
      sp.expand(Hx.T * J * Hx + Hy.T * J * Hy + aa_**2 * J - J - (dx + dy + aa_**2 - 1) * J).tolist())
X2 = Hy * J * Hx.T * J
print('(1,1) の第 2 項: (Hy J Hx^T J)^T J (Hy J Hx^T J) - dx dy J =', sp.expand(X2.T * J * X2 - dx * dy * J).tolist())
numer = sp.expand((u - dx)**2 + dx * dy + dx * u**2 - u**2)
print('(1,1): u² × (係数 - 1) = (u - dx)² + dx dy + dx u² - u² - dx (dx + dy + a² - 1) =', sp.expand(numer - dx * (dx + dy + aa_**2 - 1)))
print('  (展開: u = 1 + a。-2u dx + dx² + dx dy + dx u² = dx (dx + dy + u² - 2u) = dx (dx + dy + a² - 1))')
print('(1,3): -(1/u) (Hy J Hx^T J)^T J Hy - J Hx 項: J Hx の係数の和 (-u + dx + dy + a u) / u -> 本文の分子 - (dx + dy + a² - 1) =',
      sp.expand(-u + dx + dy + aa_ * u - (dx + dy + aa_**2 - 1)))
print('(1,3) の各項を明示: H11^T J H13 + H21^T J H23 + H31^T J H33 - J Hx (-(u-dx)/u + dy/u + a) =',
      sp.simplify(((1 - dx / u) * I2).T * J * (-Hx) + (Hy * J * Hx.T * J / u).T * J * (-Hy) + (-J * Hx.T * J).T * J * (aa_ * I2)
                  - J * Hx * (-(u - dx) / u + dy / u + aa_)).tolist())
J6 = sp.diag(J, J, J)
M = (H.T * J6 * H - J6)
# det 条件を y11 について解いて代入 (dx + dy = 1 - a² => y11 = (1 - a² - dx + y12 y21) / y22)
subsd = {y11: (1 - aa_**2 - dx + y12 * y21) / y22}
names = {(0, 0): '(1,1)', (0, 1): '(1,2)', (0, 2): '(1,3)', (1, 1): '(2,2)', (1, 2): '(2,3)', (2, 2): '(3,3)'}
for (bi, bj), nm in names.items():
    blk = M[2 * bi:2 * bi + 2, 2 * bj:2 * bj + 2].subs(subsd).applyfunc(sp.simplify)
    print('  記号計算 (y11 を a² + det Hx + det Hy = 1 から消去) ブロック %s の H^T J6 H - J6 =' % nm, blk.tolist())
rng = np.random.default_rng(5)
worst = 0
for trial in range(30):
    hx = [sp.Rational(int(rng.integers(-5, 6)), 7) for _ in range(4)]
    h12, h21 = sp.Rational(int(rng.integers(-5, 6)), 7), sp.Rational(int(rng.integers(-5, 6)), 7)
    h22 = sp.Rational(int(rng.integers(1, 6)), 7)
    dxv = hx[0] * hx[3] - hx[1] * hx[2]
    h11 = (1 - sp.Rational(9, 25) - dxv + h12 * h21) / h22
    vals = {x11: hx[0], x12: hx[1], x21: hx[2], x22: hx[3], y11: h11, y12: h12, y21: h21, y22: h22, aa_: sp.Rational(3, 5)}
    Mv = M.subs(vals)
    worst = max(worst, max(abs(v) for v in Mv))
print('有理数 30 組 (a = 3/5、a² + det Hx + det Hy = 1 を厳密に満たす) で H^T J6 H - J6 の全 36 成分の最大 |値| =', worst, '(厳密な有理数の計算)')

# (F 続き) (1,2) ブロックの中間式
barX = lambda X: J * X.T * J
Hbx, Hby = barX(Hx), barX(Hy)
px_, py_ = 1 - dx / u, 1 - dy / u
t1 = ((px_ * I2).T * J * (Hx * Hby / u))
t2 = ((Hy * Hbx / u).T * J * (py_ * I2))
t3 = ((-Hbx).T * J * (-Hby))
print('(1,2) の 3 項: 第 1 項 - (p_x/u) J Hx barHy =', sp.simplify(t1 - px_ / u * J * Hx * Hby).tolist(),
      '、第 2 項 - (p_y/u) J Hx barHy =', sp.simplify(t2 - py_ / u * J * Hx * Hby).tolist(),
      '、第 3 項 + J Hx barHy =', sp.simplify(t3 + J * Hx * Hby).tolist())
print('  p_x + p_y - u (a² + dx + dy = 1 のもと) =', sp.simplify((px_ + py_ - u).subs(y11, (1 - aa_**2 - dx + y12 * y21) / y22)))
print('  bar X の転置 = J X J、det(bar X) = det X:', sp.expand(Hbx.T - J * Hx * J).tolist(), sp.expand(Hbx.det() - dx))
print('  H_x J_2 H_y^T J_2 = H_x bar H_y:', sp.expand(Hx * J * Hy.T * J - Hx * Hby).tolist())
print('  (1,3) の第 2 項: H21^T J H23 - (dy/u) J Hx =', sp.simplify((Hy * Hbx / u).T * J * (-Hy) - dy / u * J * Hx).tolist())

# ------------------------------------------------------------------ (G)
print('\n== (G) qmdiag.f の det N ≠ 0 の分岐: r = X_f N^{-1} ==')
M13, M14, M23, M24, A_, B_, C_ = sp.symbols('M13 M14 M23 M24 aaa bbb ccc')
detn_ = M13 * M24 - M14 * M23
xx11, xx12, xx21, xx22 = A_, B_ + C_, -B_ + C_, -A_
r_fortran = sp.Matrix([[(xx11 * M24 - xx12 * M23) / detn_, (-xx11 * M14 + xx12 * M13) / detn_],
                       [(xx21 * M24 - xx22 * M23) / detn_, (-xx21 * M14 + xx22 * M13) / detn_]])
Xf = sp.Matrix([[xx11, xx12], [xx21, xx22]])
Nm = sp.Matrix([[M13, M14], [M23, M24]])
print('Fortran の r1..r4 - X_f N^{-1} =', sp.simplify(r_fortran - Xf * Nm.inv()).tolist())
yy11, yy12, yy21, yy22 = A_ + B_, C_, C_, A_ - B_
r_f2 = sp.Matrix([[yy11 * (-M23) + yy12 * (-M24), yy11 * M13 + yy12 * M14], [yy21 * (-M23) + yy22 * (-M24), yy21 * M13 + yy22 * M14]])
Yf = sp.Matrix([[yy11, yy12], [yy21, yy22]])
print('detn = 0 の分岐: Fortran の r1..r4 - Y_f [[-M23, M13], [-M24, M14]] =', sp.simplify(r_f2 - Yf * sp.Matrix([[-M23, M13], [-M24, M14]])).tolist())
