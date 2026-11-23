# -*- coding: utf-8 -*-
"""14_emittance.md の途中式の補足 (sympy・numpy。pysad は使わない)。

(a) §4.2 (14.6a): スライスの減衰行列 D_i の跡 = -4u (p_x, p_y の縮み -u ずつ、δδ が -2u)。
(b) §4.4 (14.8b): rd = -aP + S を Fortran の btr(k,k) の式に代入する手計算の各項、Λ'^TΛ' の分解による「和-和」の係数。
(c) §4.4 (14.8a): Σ = [[A,b],[b,C]] への PΣP^T の成分 (和・差・交差)。
(d) §4.2 (14.6): ℓ D r と K64 の cd の関係 (共役)、λ = e^{-iμ}。
(e) §3.1 (14.2b): 級数の項ずらし Σ = Δ + MΣM^T の数値確認。
(f) §5: RF の 1 周の写像 (φ = ωz/c の規約)、ハミルトニアン、固定点、δ_max²、数値 (3 GeV、3 MV、h = 64)、
    写像の追跡でセパラトリクスの内外を確かめる。縦の 2x2 の P(μ) の符号と R56 の符号。
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)

print('== (a) スライスの減衰行列 D_i の跡 ==')
# 1 スライス: 入口の変数 (x, px, y, py, δ)。q = p_x/p は光子放出で変わらない (光子は速度の向き)。
# p = 1 + δ。スライスの損失 dg = -u(1+δ)^2 (u は δ=0 での損失率)。
# 根拠: dg ∝ h1^4 θ²/(p0 al1) で h1 ≈ p0 p (γ >> 1)、θ = Δp_⊥/p (磁場の衝撃 Δp_⊥ は p に依らない) なので dg ∝ p^4 p^-2 = p²。
x, px, y, py, dl, u = sp.symbols('x p_x y p_y delta u')
p = 1 + dl
dg = -u * p**2
qx, qy = px / p, py / p
# 出口: δ' = δ + dg、p_x' = q_x (1 + δ')
dl2 = dl + dg
px2 = qx * (1 + dl2)
py2 = qy * (1 + dl2)
# 放射による増分 (放射なしの写像は恒等として、増分のヤコビアン D_i = d(出口 - 入口)/d(入口))
inc = sp.Matrix([0, px2 - px, 0, py2 - py, dl2 - dl])
vars5 = [x, px, y, py, dl]
Dm = inc.jacobian(vars5).subs(dl, 0)
print('D_i (δ = 0 で評価、変数 x, px, y, py, δ の順):')
sp.pprint(sp.simplify(Dm))
print('D_pxpx =', sp.simplify(Dm[1, 1]), '、D_pypy =', sp.simplify(Dm[3, 3]), '、D_δδ =', sp.simplify(Dm[4, 4]))
print('Tr D_i =', sp.simplify(Dm.trace()), ' (z は D_zz = O(u/γ²) で省く。-4u)')
# δ ≠ 0 の一般式
Dg = inc.jacobian(vars5)
print('δ ≠ 0 の跡 =', sp.simplify(Dg.trace()), ' (δ = 0 で', sp.simplify(Dg.trace().subs(dl, 0)), ')')
# h1 を厳密に (h1 = sqrt(1 + p0² p²)) として D_δδ = d dg/dδ: d ln(h1^4/p²)/dδ = 4 p0² p²/h1² - 2 → γ >> 1 で 2
p0 = sp.symbols('p0', positive=True)
h1 = sp.sqrt(1 + p0**2 * p**2)
dgx = -u * h1**4 / p**2 / (1 + p0**2)**2   # δ=0 で -u
Dd = sp.simplify(sp.diff(dgx, dl).subs(dl, 0) / (-u))
print('h1 を厳密にしたときの D_δδ/(-u) =', Dd, '(p0 = 5870 (3 GeV): ', float(Dd.subs(p0, 5870.0)), ')')

print('h1 を厳密にしたときの D_δδ = -u (2 - 4/(1+p0²)) と γ >> 1 の -2u との差は 4/γ² (相対)')
# 力の補正項 ddpx = -½ dp_x dg、x += ddpx L/3 (13 章 §4.3)。スライス内の磁場の衝撃 dp_x = -k L x (x だけに依る)
k, L = sp.symbols('k L', positive=True)
dpx_kick = -k * L * x
ddpx = -sp.Rational(1, 2) * dpx_kick * dg
inc2 = sp.Matrix([ddpx * L / 3, px2 - px + ddpx, 0, py2 - py, dl2 - dl])
D2 = inc2.jacobian(vars5).subs(dl, 0)
print('力の補正を含めた D_i (x, px, y, py, δ):')
sp.pprint(sp.simplify(D2))
print('Tr D_i =', sp.simplify(D2.trace()), ' → -4u (1 + kL²/24): 補正は (スライス長/焦点距離)² 程度で、薄いスライスで 0')

print('\n== (b) (14.8b) の手計算 ==')
a, mu, s1, s2 = sp.symbols('a mu s1 s2', real=True)
c, s = sp.cos(mu), sp.sin(mu)
P = sp.Matrix([[c, s], [-s, c]])
S = sp.Matrix([[s1, s2], [s2, -s1]])
rd = -a * P + S
for (i, j) in ((0, 0), (0, 1), (1, 0), (1, 1)):
    print('rd[%d,%d] =' % (i, j), sp.expand(rd[i, j]))
sq = sp.expand(rd[0, 0]**2 + rd[0, 1]**2 + rd[1, 0]**2 + rd[1, 1]**2)
print('4 成分の 2 乗和 =', sp.simplify(sq), ' (= 2a² + 2(s1² + s2²))')
t1 = sp.simplify(-sq / 2)
t2 = sp.simplify(-c * (rd[0, 0] + rd[1, 1]))
t3 = sp.simplify(-s * (rd[0, 1] - rd[1, 0]))
print('第 1 項 -½Σrd² =', t1)
print('rd_ii + rd_{i+1,i+1} =', sp.simplify(rd[0, 0] + rd[1, 1]), '、第 2 項 -cosμ(...) =', t2)
print('rd_{i,i+1} - rd_{i+1,i} =', sp.simplify(rd[0, 1] - rd[1, 0]), '、第 3 項 -sinμ(...) =', t3)
print('合計 =', sp.simplify(t1 + t2 + t3))
# Λ' = P + rd = (1-a)P + S の Λ'^T Λ'
Lam = P + rd
LL = sp.simplify(Lam.T * Lam)
print('Λ\'^T Λ\' =')
sp.pprint(LL)
print('  (1,1)+(2,2) の半分 =', sp.simplify((LL[0, 0] + LL[1, 1]) / 2), ' (= (1-a)² + s1² + s2²)')
print('  トレースなしの部分の (1,1) - (2,2) 半分 =', sp.simplify((LL[0, 0] - LL[1, 1]) / 2))
print('  → tr(Λ\'ΣΛ\'^T) = tr(Λ\'^TΛ\' Σ) = [(1-a)² + s1² + s2²] trΣ + (トレースなしの対称行列)·Σ')
print('  Λ\'^TΛ\' - [(1-a)²+s1²+s2²] 1 のトレース =', sp.simplify(LL.trace() - 2 * ((1 - a)**2 + s1**2 + s2**2)))

print('\n== (c) (14.8a): PΣP^T の成分 ==')
A, b, Cc = sp.symbols('A b C')
Sig = sp.Matrix([[A, b], [b, Cc]])
Q = P * Sig * P.T
q11, q22, q12 = [sp.expand(sp.simplify(Q[i, j])) for (i, j) in ((0, 0), (1, 1), (0, 1))]
print('(PΣP^T)_11 =', q11)
print('(PΣP^T)_22 =', q22)
print('(PΣP^T)_12 =', q12)
print('和 22+11 =', sp.simplify(q11 + q22))
print('差 22-11 =', sp.simplify(q22 - q11), ' vs cos2μ(C-A) - 2 sin2μ b:', sp.simplify(q22 - q11 - (sp.cos(2 * mu) * (Cc - A) - 2 * sp.sin(2 * mu) * b)) == 0)
print('交差 12 =', sp.simplify(q12), ' vs sin2μ(C-A)/2 + cos2μ b:', sp.simplify(q12 - (sp.sin(2 * mu) * (Cc - A) / 2 + sp.cos(2 * mu) * b)) == 0)
u_, v_, w_ = sp.symbols('u v w')
print('v\' = cos2μ v - sin2μ w、w\' = sin2μ v + cos2μ w (v = (C-A)/√2、w = √2 b):')
vp = (q22 - q11) / sp.sqrt(2)
wp = sp.sqrt(2) * q12
print('  v\' の差 =', sp.simplify(vp - (sp.cos(2 * mu) * (Cc - A) / sp.sqrt(2) - sp.sin(2 * mu) * sp.sqrt(2) * b)))
print('  w\' の差 =', sp.simplify(wp - (sp.sin(2 * mu) * (Cc - A) / sp.sqrt(2) + sp.cos(2 * mu) * sp.sqrt(2) * b)))

print('\n== (d) ℓ D r と cd の関係 ==')
rng = np.random.default_rng(7)
muv = 0.83
Pn = np.array([[math.cos(muv), math.sin(muv)], [-math.sin(muv), math.cos(muv)]])
Dn = 1e-3 * rng.standard_normal((2, 2))
RX = np.array([1.0, 0.0]); RP = np.array([0.0, 1.0])
r = (RX - 1j * RP) / 2
ell = np.array([1.0, 1j])
lam = np.exp(-1j * muv)
print('T r = λ r (λ = e^{-iμ}):', np.allclose(Pn @ r, lam * r), '、ℓ r =', ell @ r)
dl_over_lam = (ell @ Dn @ r) / lam
num = 0.5 * ((Dn[0, 0] + Dn[1, 1]) + 1j * (Dn[0, 1] - Dn[1, 0]))
cc = np.exp(1j * muv)
cd = num / cc
print('δλ/λ = ℓDr/λ =', dl_over_lam, '、K64 の cd =', cd, '、共役:', np.allclose(cd, np.conj(dl_over_lam)))
print('a = -Re cd =', -cd.real, '、-Re(δλ/λ) =', -dl_over_lam.real, '、Im cd = -Im(δλ/λ):', np.isclose(cd.imag, -dl_over_lam.imag))
# 固有値の厳密な変化との比較
w_ex = np.linalg.eigvals(Pn + Dn)
k = int(np.argmin(abs(w_ex - lam)))
print('厳密な固有値の変化 δλ/λ =', (w_ex[k] - lam) / lam, '(1 次の ℓDr/λ との差は 2 次)')

print('\n== (e) (14.2b) の項ずらし ==')
M = np.array([[0.6, 0.5, 0, 0], [-0.5, 0.6, 0, 0], [0, 0, 0.7, 0.2], [0, 0, -0.2, 0.8]]) * 0.97
M[0, 2] = 0.05; M[3, 1] = -0.03
Dl = np.diag([1.0, 2.0, 0.5, 1.5]); Dl[0, 1] = Dl[1, 0] = 0.3
Ssum = np.zeros((4, 4)); term = Dl.copy()
for n in range(4000):
    Ssum += term
    term = M @ term @ M.T
print('Σ_{n>=0} M^n Δ M^nT - [Δ + M Σ M^T] の最大 =', abs(Ssum - (Dl + M @ Ssum @ M.T)).max())
vec = np.linalg.solve(np.eye(16) - np.kron(M, M), Dl.reshape(-1))
print('(1 - M⊗M) vecΣ = vecΔ の解との差 =', abs(vec.reshape(4, 4) - Ssum).max())

print('\n== (f) RF の力学: φ = ωz/c の規約 ==')
# (1) 規約: z = -v(t - t0) (01 章 (1.7))、φ ≡ ω_RF z/c = 2π h z / C = -ω_RF t。10 章 (10.4): Δh = -v sin(ωt - φs) = v sin(φs + φ)
# (2) 経路: z → z - α_p C δ (R56eff = -α_p C) なので φ → φ - 2π h α_p δ
V, E, h, Ccirc = 3e6, 3e9, 64, 41.6
U0 = 2.3449e6
ap_I1 = 16.5806 / Ccirc
ap_c6 = 0.39780
sphi = U0 / V
phis = math.asin(sphi)
print('sinφs = U0/V = %.5f、φs = %.5f rad、cosφs = %.5f' % (sphi, phis, math.cos(phis)))
kap = V / E * math.cos(phis) * 2 * math.pi * h / Ccirc
print('κ = d(Δδ)/dz = (V/E) cosφs 2πh/C = %.5e /m (14 章の dvcacc/pgev = 6.0293e-3)' % kap)
# ハミルトニアン H = -π h α_p δ² + (V/E)[cos(φs+φ) + φ sinφs]
phs, ph_, dd, apv, Vv, Ev, hh = sp.symbols('phi_s phi delta alpha_p V E h', positive=True)
H = -sp.pi * hh * apv * dd**2 + Vv / Ev * (sp.cos(phs + ph_) + ph_ * sp.sin(phs))
print('dφ/dn = ∂H/∂δ =', sp.diff(H, dd), '、dδ/dn = -∂H/∂φ =', sp.simplify(-sp.diff(H, ph_)))
print('∂²H/∂φ² (φ=0) =', sp.simplify(sp.diff(H, ph_, 2).subs(ph_, 0)), ' (< 0 で、-πhα_pδ² (< 0) と同じ向きに凸: 安定固定点 (鞍点でなく最大)')
phu = sp.pi - 2 * phs
print('不安定固定点: sin(φs+φu) - sinφs =', sp.simplify(sp.sin(phs + phu) - sp.sin(phs)))
Hu = sp.simplify(H.subs({dd: 0, ph_: phu}))
H0 = sp.simplify(H.subs({dd: 0, ph_: 0}))
print('H(φu, 0) - H(0, 0) =', sp.simplify(Hu - H0))
# セパラトリクス上で φ = 0: -πhα_p δmax² + H0 = Hu
dmax2 = sp.solve(sp.Eq(H.subs(ph_, 0), Hu), dd**2) if False else None
dm2 = sp.symbols('dm2')
sol = sp.solve(sp.Eq(-sp.pi * hh * apv * dm2 + Vv / Ev * sp.cos(phs), Hu), dm2)[0]
print('δ_max² =', sp.simplify(sol))
print('  = V/(π h α_p E) [2cosφs - (π-2φs) sinφs]:', sp.simplify(sol - Vv / (sp.pi * hh * apv * Ev) * (2 * sp.cos(phs) - (sp.pi - 2 * phs) * sp.sin(phs))) == 0)
def dm(ap):
    return math.sqrt(V / (math.pi * h * ap * E) * (2 * math.cos(phis) - (math.pi - 2 * phis) * sphi))
print('δ_max (α_p = I1/C = %.5f) = %.5e (14 章 §5: 1.55828e-3)、δ_max (calc6d の α_p = %.5f) = %.5e (1.55977e-3)' % (ap_I1, dm(ap_I1), ap_c6, dm(ap_c6)))
# 写像を追跡して、セパラトリクスの内外を確かめる (1 周の写像: δ += (V/E)(sin(φs+φ) - sinφs)、φ -= 2π h α_p δ)
def survive(d0, nturn=60000):
    ph, de = 0.0, d0
    for n in range(nturn):
        de += V / E * (math.sin(phis + ph) - sphi)
        ph -= 2 * math.pi * h * ap_I1 * de
        if abs(ph) > 12:
            return False, n
    return True, nturn
dmax = dm(ap_I1)
for f in (0.8, 0.95, 0.99, 1.02, 1.1):
    ok, n = survive(f * dmax)
    print('  δ0 = %.2f δ_max (φ = 0 から): %s' % (f, ('有界 (%d 周)' % n) if ok else ('逃げる (%d 周で |φ| > 12)' % n)))
# 縦の 2x2 の回転の向き: P(μ) の sinμ と R56 の符号
print('\n縦の 2x2 [[1, R56],[κ, 1+κR56]] を R P(μ) R^{-1} (R = [[a, 0],[c, 1/a]]) と書くと M12 = a² sinμ:')
aa, cc2, mm = sp.symbols('a c mu', real=True)
R2 = sp.Matrix([[aa, 0], [cc2, 1 / aa]])
Pm = sp.Matrix([[sp.cos(mm), sp.sin(mm)], [-sp.sin(mm), sp.cos(mm)]])
M2 = sp.simplify(R2 * Pm * R2.inv())
print('  M12 =', sp.simplify(M2[0, 1]), ' → R56 の符号 = sinμ の符号 (α_p > 0 なら sinμ < 0、Im cd(6) = μ < 0)')
