# -*- coding: utf-8 -*-
r"""08_mult.md の途中式の確認 (sympy / numpy)。

(a) Horner 形 (tmulk・tmulti の cx) = C(z)、tmulk の cx = -(tan psi/L) z C(z)
(b) 回転 (8.4)
(c) 端面の角の式 (8.9): 解析性 (コーシー・リーマン)・ポテンシャル・K0、K1 の展開、「長さ x tan(psi) × 場」だけの式との違い
(d) ずれ (8.10): tsolrot の式 (K64 tsol.f、bz = 0) が R2 R1 の剛体変換に一致すること、入口と出口が逆変換であること
(e) A_s の漸化式 (8.14) の導出の係数、閉じた形、キック (8.16)
(f) Strang 分割の誤差 (調和振動子)。ndiv の式 (8.8) の係数 6 の位置づけ
"""
import math
import numpy as np
import sympy as sp
import _common as C

C.start(__file__)
x, y = sp.symbols('x y', real=True)
z = x + sp.I * y
Lq, t = sp.symbols('L t', positive=True)

print('== (a) Horner 形 ==')
N = 4
a = sp.symbols('a0:%d' % (N + 1))           # 複素数 K_n + i SK_n
zz = sp.symbols('z')
cx = 0
for k in range(N, 0, -1):
    cx = (cx + a[k]) * zz / k                # tmulk: cx=(cx+akn(kord))*cx0*aninv(kord)
Cz = sum(a[n] * zz**n / sp.factorial(n) for n in range(N + 1))
print('   ループ後の cx + a0 - C(z) =', sp.simplify(cx + a[0] - Cz))
fe = -t / Lq
cx_tmulk = fe * (cx + a[0]) * zz            # tmulk: cx=fe*(cx+akn(0))*cx0
print('   tmulk の cx - (-(tan psi/L) z C(z)) =', sp.simplify(cx_tmulk + t / Lq * zz * Cz))

print('\n== (b) 回転 (8.4): z -> z\' = e^{i theta} z で、係数が a e^{-i(n+1) theta} になること ==')
th = sp.symbols('theta', real=True)
for n in (1, 2, 3):
    ar, ai = sp.symbols('ar ai', real=True)
    A = ar + sp.I * ai
    lhs = sp.re(sp.expand(A * z**(n + 1)))
    zp = sp.exp(sp.I * th) * z
    rhs = sp.re(sp.expand(A * sp.exp(-sp.I * (n + 1) * th) * zp**(n + 1)))
    print('   n = %d: 差 =' % n, sp.simplify(lhs - rhs))

print('\n== (c) 端面の角の式 (8.9) ==')
F = t / Lq * zz * Cz                          # F = Δpx - iΔpy
Fxy = sp.expand(F.subs(zz, z))
u = sp.re(Fxy)
v = -sp.im(Fxy)                               # Δpx = Re F、Δpy = -Im F
cxm = sp.expand(cx_tmulk.subs(zz, z))
dpx_tmulk = -sp.re(cxm)
dpy_tmulk = sp.im(cxm)
print('   Δpx (tmulk) - Re F =', sp.simplify(dpx_tmulk - u), '   Δpy (tmulk) - (-Im F) =', sp.simplify(dpy_tmulk - v))
print('   コーシー・リーマン: ∂x Re F - ∂y Im F =', sp.simplify(sp.diff(u, x) - sp.diff(sp.im(Fxy), y)),
      '、∂y Re F + ∂x Im F =', sp.simplify(sp.diff(u, y) + sp.diff(sp.im(Fxy), x)))
print('   rot(Δp) = ∂x Δpy - ∂y Δpx =', sp.simplify(sp.diff(v, x) - sp.diff(u, y)), ' (0: 勾配。運動量だけの kick なのでシンプレクティック)')
print('   div(Δp) = ∂x Δpx + ∂y Δpy =', sp.simplify(sp.diff(u, x) + sp.diff(v, y)), ' (0: 調和)')
Phi = -t / Lq * sp.re(sp.expand(sum(a[n] * z**(n + 2) / ((n + 2) * sp.factorial(n)) for n in range(N + 1))))
print('   Φ = -(t/L) Re Σ a_n z^(n+2)/((n+2) n!) として ∂xΦ + Δpx =', sp.simplify(sp.diff(Phi, x) + u), '、∂yΦ + Δpy =', sp.simplify(sp.diff(Phi, y) + v), ' (Δp = -∇Φ)')
print('   ラプラシアン Φ =', sp.simplify(sp.diff(Phi, x, 2) + sp.diff(Phi, y, 2)))
k0, k1 = sp.symbols('K0 K1', real=True)
only0 = {a[0]: k0, a[1]: 0, a[2]: 0, a[3]: 0, a[4]: 0}
only1 = {a[0]: 0, a[1]: k1, a[2]: 0, a[3]: 0, a[4]: 0}
print('   K0 のみ: Δpx, Δpy =', sp.simplify(u.subs(only0)), ',', sp.simplify(v.subs(only0)))
print('   K1 のみ: Δpx, Δpy =', sp.simplify(u.subs(only1)), ',', sp.simplify(v.subs(only1)))
Cxy = sp.expand(Cz.subs(zz, z))
nx = sp.simplify((-x * t) * (-sp.re(Cxy)) / Lq)     # 「余分な長さ (-x tanψ) × 力 (-Re C, +Im C)/L」だけの式
ny = sp.simplify((-x * t) * sp.im(Cxy) / Lq)
print('   「余分な長さ × 力」だけの式 (K1 のみ): Δpx, Δpy =', sp.simplify(nx.subs(only1)), ',', sp.simplify(ny.subs(only1)))
rot_n = sp.diff(ny, x) - sp.diff(nx, y)
print('   その rot: K1 のみ =', sp.simplify(rot_n.subs(only1)), '、K0 のみ =', sp.simplify(rot_n.subs(only0)))
allk = {a[1]: k1, a[2]: 0, a[3]: 0, a[4]: 0, a[0]: k0}
print('   y = 0 (中央面): (8.9) と「余分な長さ × 力」の差 Δpx:', sp.simplify((u - nx).subs(y, 0).subs(allk)), '、Δpy:', sp.simplify((v - ny).subs(y, 0).subs(allk)))

print('\n== (d) ずれ (8.10) と tsolrot (bz = 0) ==')
dx_, dy_, dz_, alg = sp.symbols('DX DY DZ alg', real=True)
ch1, ch2 = sp.symbols('chi1 chi2', real=True)
c1, s1, c2, s2 = sp.cos(ch1), sp.sin(ch1), sp.cos(ch2), sp.sin(ch2)
R1 = sp.Matrix([[c1, 0, -s1], [0, 1, 0], [s1, 0, c1]])
R2 = sp.Matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
X, Y = sp.symbols('X Y', real=True)
s0 = -alg
vv0 = sp.Matrix([X - dx_, Y - dy_, s0 - dz_])
vv2 = R2 * R1 * vv0
xc = X - dx_
yc = Y - dy_
ds1 = s1 * xc + (c1 - 1) * s0 - c1 * dz_
xn = c1 * xc - s1 * (s0 - dz_)
ds2 = s2 * yc + c2 * ds1 + (c2 - 1) * s0
yn = c2 * yc - s2 * (ds1 + s0)
print('   x: コード - (R2 R1 v0)_x =', sp.simplify(xn - vv2[0]), '、y:', sp.simplify(yn - vv2[1]),
      '、第 3 成分: ds2 - ((R2 R1 v0)_s - s0) =', sp.simplify(ds2 - (vv2[2] - s0)))
px_, py_ = sp.symbols('px py', real=True)
pz0 = sp.sqrt(1 - px_**2 - py_**2)
pv = R2 * R1 * sp.Matrix([px_, py_, pz0])
pz1 = s1 * px_ + c1 * pz0
print('   px: コード (c1 px - s1 pz0) - (R2 R1 p)_x =', sp.simplify(c1 * px_ - s1 * pz0 - pv[0]),
      '、py: (c2 py - s2 pz1) - (R2 R1 p)_y =', sp.simplify(c2 * py_ - s2 * pz1 - pv[1]))


def ent(vec, alg, dx, dy, dz, chi1, chi2):
    X_, px, Y_, py, zc = vec
    s0_ = -alg
    C1, S1, C2, S2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    X_ -= dx
    Y_ -= dy
    p0 = math.sqrt(1 - px * px - py * py)
    q1 = S1 * px + C1 * p0
    px = C1 * px - S1 * p0
    py = C2 * py - S2 * q1
    d1 = S1 * X_ + (C1 - 1) * s0_ - C1 * dz
    X_ = C1 * X_ - S1 * (s0_ - dz)
    d2 = S2 * Y_ + C2 * d1 + (C2 - 1) * s0_
    Y_ = C2 * Y_ - S2 * (d1 + s0_)
    q2 = math.sqrt(1 - px * px - py * py)
    aa = d2 / q2
    return [X_ - aa * px, px, Y_ - aa * py, py, zc + aa]


def ext(vec, alg, dx, dy, dz, chi1, chi2):          # tsol.f の出口 (引数 alg に alg - al を渡す)
    X_, px, Y_, py, zc = vec
    s0_ = alg
    C1, S1, C2, S2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    p0 = math.sqrt(1 - px * px - py * py)
    q1 = -S2 * py + C2 * p0
    pyi = C2 * py + S2 * p0
    q2 = -S1 * px + C1 * q1
    pxi = C1 * px + S1 * q1
    y1 = C2 * Y_ - S2 * s0_
    d1 = -S2 * Y_ - (C2 - 1) * s0_
    d2 = -S1 * X_ + C1 * d1 - (C1 - 1) * s0_ + dz
    x1 = C1 * X_ + S1 * (d1 - s0_)
    aa = d2 / q2
    return [x1 - aa * pxi + dx, pxi, y1 - aa * pyi + dy, pyi, zc + aa]


print('   入口 → (要素内で長さ L の直進) → 出口 が、実験室での長さ L の直進に一致するか (空の要素。alg = L/2、出口は alg - L を渡す):')
Lel = 2.0
for (dx, dy, dz, c1v, c2v) in [(0, 0, 0, 0.01, 0), (0, 0, 0, 0, 0.01), (0, 0, 0, 0.01, 0.02), (1e-3, 2e-3, 3e-3, 0.01, 0.02), (1e-3, 2e-3, 3e-3, -0.03, 0.05)]:
    v0_ = [1e-3, 2e-3, -1e-3, 1e-3, 0.0]
    w_ = ent(v0_, Lel / 2, dx, dy, dz, c1v, c2v)
    pzz = math.sqrt(1 - w_[1]**2 - w_[3]**2)
    w_ = [w_[0] + Lel * w_[1] / pzz, w_[1], w_[2] + Lel * w_[3] / pzz, w_[3], w_[4]]
    w_ = ext(w_, Lel / 2 - Lel, dx, dy, dz, c1v, c2v)
    pz0_ = math.sqrt(1 - v0_[1]**2 - v0_[3]**2)
    ref = [v0_[0] + Lel * v0_[1] / pz0_, v0_[1], v0_[2] + Lel * v0_[3] / pz0_, v0_[3]]
    print('     (DX, DY, DZ, chi1, chi2) = %s: 最大差 (x, px, y, py) = %.1e' % ((dx, dy, dz, c1v, c2v), max(abs(w_[i] - ref[i]) for i in range(4))))
print('   (差 1e-17 以下: 入口と出口は、中心 (alg = L/2) を支点とする同じ剛体変換の順・逆である。)')

print('\n== (e) A_s の漸化式 (8.14)・閉じた形・キック (8.16) ==')
a_, m_ = sp.symbols('a m')
xs, ys, rho = sp.symbols('x y rho', positive=True)
Gxy = (rho + xs)**a_ * (xs + sp.I * ys)**m_
lap = sp.diff(Gxy, xs, 2) + sp.diff(Gxy, ys, 2)
rhs = sp.diff(Gxy, xs) / (rho + xs)
res = sp.simplify(sp.expand_power_base(sp.simplify((lap - rhs) / ((rho + xs)**(a_ - 2) * (xs + sp.I * ys)**(m_ - 1))), force=True))
expect = a_ * (a_ - 2) * (xs + sp.I * ys) + m_ * (2 * a_ - 1) * (rho + xs)
print('   [∇² - ∂x/u](u^a z^m) / (u^(a-2) z^(m-1)) - [a(a-2) z + m(2a-1) u] =', sp.simplify(sp.expand(res - expect)))
kk, nn = sp.symbols('k n', integer=True, nonnegative=True)
ak_ = sp.Rational(1, 2) - kk
rec = sp.simplify(-(ak_ * (ak_ - 2)) / ((nn + 1 + kk + 1) * (2 * (ak_ - 1) - 1)))
print('   g_{k+1}/g_k = -a_k(a_k-2)/((m+k+1)(2a_{k+1}-1)) (m = n+1) =', sp.factor(rec))
print('   a_k(a_k-2) =', sp.factor(sp.expand(ak_ * (ak_ - 2))), '、2a_{k+1}-1 =', sp.simplify(2 * (ak_ - 1) - 1))
ok = True
for n in range(0, 8):
    for k in range(0, 8):
        prod = sp.Integer(1)
        for j in range(1, k + 1):
            prod *= sp.Rational((2 * j - 3) * (2 * j + 1), 8 * j * (j + n + 1))

        def df(mm):
            return {-3: -1, -1: 1}[mm] if mm < 0 else sp.factorial2(mm)
        closed = (-df(2 * k + 1) * df(2 * k - 3) * sp.factorial(n + 1) / (8**k * sp.factorial(k) * sp.factorial(n + k + 1))) if k >= 1 else sp.Integer(1)
        ok = ok and (prod == closed)
print('   n, k ≤ 7 で、積 Π_{j=1..k}(2j-3)(2j+1)/(8j(j+n+1)) と閉じた形 -(2k+1)!!(2k-3)!!(n+1)!/(8^k k!(n+k+1)!) が一致:', ok)
dxF = sp.simplify(sp.diff(Gxy, xs) - (a_ * Gxy / (rho + xs) + m_ * (xs + sp.I * ys)**(m_ - 1) * (rho + xs)**a_))
dyF = sp.simplify(sp.diff(Gxy, ys) - sp.I * m_ * (xs + sp.I * ys)**(m_ - 1) * (rho + xs)**a_)
print('   ∂x(u^a z^m) = a u^(a-1) z^m + m u^a z^(m-1)、の差 =', dxF, '、∂y(u^a z^m) = i m u^a z^(m-1)、の差 =', dyF)

print('\n== (f) Strang 分割 (調和振動子) ==')
w, h = sp.symbols('omega h', positive=True)
D = lambda s_: sp.Matrix([[1, s_], [0, 1]])
Kk = lambda s_: sp.Matrix([[1, 0], [-w**2 * s_, 1]])
M = D(h / 2) * Kk(h) * D(h / 2)
Mex = sp.Matrix([[sp.cos(w * h), sp.sin(w * h) / w], [-w * sp.sin(w * h), sp.cos(w * h)]])
print('   M = D(h/2) K(h) D(h/2) =', M.applyfunc(sp.expand).tolist())
print('   M11 - 厳密 =', sp.series(M[0, 0] - Mex[0, 0], h, 0, 6).removeO(), '、M12 - 厳密 =', sp.series(M[0, 1] - Mex[0, 1], h, 0, 5).removeO())
cs = (M[0, 0] + M[1, 1]) / 2
print('   cos θ = tr M / 2 =', sp.expand(cs), ' → θ の展開:', sp.series(sp.acos(cs), h, 0, 5).removeO())
print('   n 個の刻みの位相誤差 = n ω h (ωh)²/24 = (ωL)³/(24 ndiv²)。(ωL)² = Q = a^(n-1)|K_n| L/(n-1)! とおくと Q^{3/2}/(24 ndiv²)。')
print('   (8.8) は ε = Q/(6 ndiv²)。次数 (ndiv^-2) と Q への依存の仕方 (∝ Q の冪) は合うが、係数 1/6 と Q の冪 (1 対 3/2) はこの見積りから出ない。係数 6 は未導出。')

print('\n== (f2) 対称な 3 つの積の誤差項 (BCH): exp(hX/2) exp(hY) exp(hX/2) = exp(h(X+Y) + h^3 [ -[X,[X,Y]]/24 + [Y,[Y,X]]/12 ] + ...) ==')
from scipy.linalg import expm, logm
rng = np.random.default_rng(3)
Xm = rng.normal(size=(4, 4))
Ym = rng.normal(size=(4, 4))
comm = lambda A, B: A @ B - B @ A
for hh in (0.02, 0.01):
    S = expm(hh * Xm / 2) @ expm(hh * Ym) @ expm(hh * Xm / 2)
    lhs = logm(S) - hh * (Xm + Ym)
    rhs = hh**3 * (-comm(Xm, comm(Xm, Ym)) / 24 + comm(Ym, comm(Ym, Xm)) / 12)
    print('   h = %g: |log S - h(X+Y) - h^3(...)| / |h^3(...)| = %.3e' % (hh, np.abs(lhs - rhs).max() / np.abs(rhs).max()))
print('   (h を半分にすると残りの比が 1/4 になる = 次の項が h^5。X = 本体、Y = キックとして、1 刻みの誤差は h^3 に比例。)')

print('\n== (f3) ndiv の数値と、調和振動子の見積りとの比較 (K2 = 40 /m²、L = 1 m、a_max = 0.05、ε = 0.005 EPS) ==')
Kn, Ll, amax, n2 = 40.0, 1.0, 0.05, 2
Q = amax**(n2 - 1) * Kn * Ll / math.factorial(n2 - 1)
print('   Q = a_max^(n-1) |K_n| L/(n-1)! = %g' % Q)
for eps0 in (4, 2, 1, 0.5, 0.25):
    nd = int(math.sqrt(Q / 6 / (0.005 * eps0))) + 1
    print('   EPS = %-5g ndiv = 1 + floor(sqrt(Q/(6 ε))) = %d' % (eps0, nd))
aa = 0.004
q_act = Kn * aa * Ll
print('   実際の振幅 a = 4 mm: (ωL)² = K_n a L = %g。位相誤差の見積り a (ωL)³/(24 ndiv²):' % q_act)
for nd, meas in ((5, 6.4e-7), (17, 5.5e-8)):
    print('     ndiv = %d: 見積り %.2e、測定 (mult_02 (a) の厳密解との差) %.2e、比 %.2f' % (nd, aa * q_act**1.5 / (24 * nd**2), meas, meas / (aa * q_act**1.5 / (24 * nd**2))))

print('\n== (e2) g_kn の表 (漸化式 (8.14)) ==')
for n in (0, 1, 2):
    row = []
    g_ = sp.Integer(1)
    row.append(g_)
    for k in range(1, 4):
        g_ = g_ * sp.Rational((2 * k - 3) * (2 * k + 1), 8 * k * (k + n + 1))
        row.append(g_)
    print('   n = %d: g_0n .. g_3n =' % n, row)
