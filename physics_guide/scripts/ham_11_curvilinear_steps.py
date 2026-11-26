# -*- coding: utf-8 -*-
r"""01a 章 §9 の途中式の確認 (本文に書き足した展開)。

(a) §9.2: sqrt(1-u) の二項係数 (1/2 choose k)(-u)^k。
(b) §9.3: 曲線座標 (X, Y, s) の微小な面での周回積分から磁場の成分を作る式。
    (i) 周回積分の展開 (記号)。(ii) 円弧の基準軌道で、実際の線積分 (4 辺) を数値で作り、デカルト座標の rot A と比べる。
(c) §9.4: (1a.57) から、設計軌道の条件 (1a.61)。ローレンツ力の向き (t × b = -n)。
(d) §9.4: 分散の方程式 η'' + η/ρ² = 1/ρ の解 (初期条件の解と、周期的な特解 η = ρ、一般解)。
(e) §9.6: 落とした 3 次の項 H3 の運動方程式への寄与の見積もり。
"""
import mpmath as mp
import sympy as sp
import _common as C

C.start(__file__)
R = sp.Rational

# (a) ------------------------------------------------------------------
print('(a) sqrt(1-u) の二項展開: Σ_k C(1/2,k) (-u)^k')
u = sp.symbols('u')
for k in range(0, 5):
    b = sp.binomial(R(1, 2), k)
    print('  k=%d: C(1/2,%d) = %s、C(1/2,%d)(-u)^%d = %s' % (k, k, b, k, k, sp.simplify(b * (-u)**k)))
print('  sympy series:', sp.series(sp.sqrt(1 - u), u, 0, 5))
x5 = sp.symbols('x5', positive=True)
print('  x^4 の項 5/128: 捨てる x < 8e-6 で |5 x^4/128| = %.2e' % float(R(5, 128) * (8e-6)**4))

# (b) ------------------------------------------------------------------
print('\n(b-i) 周回積分の展開 (共変成分 A_Y, A_s、微小な長方形 dY × ds、X 一定)')
Y, s, dY, ds = sp.symbols('Y s dY ds')
AY = sp.Function('A_Y')(Y, s)
As = sp.Function('A_s')(Y, s)
# (Y,s) -> (Y+dY,s) -> (Y+dY,s+ds) -> (Y,s+ds) -> (Y,s)
# 辺ごとに中点則 (2 次精度) で書く
side1 = AY.subs({Y: Y + dY / 2}) * dY                       # (Y,s) -> (Y+dY,s): A_Y dY
side2 = As.subs({Y: Y + dY, s: s + ds / 2}, simultaneous=True) * ds   # (Y+dY,s) -> (Y+dY,s+ds): A_s ds
side3 = -AY.subs({Y: Y + dY / 2, s: s + ds}, simultaneous=True) * dY   # (Y+dY,s+ds) -> (Y,s+ds): -A_Y dY
side4 = -As.subs({s: s + ds / 2}) * ds                      # (Y,s+ds) -> (Y,s): -A_s ds
tot = side1 + side2 + side3 + side4
lead = sp.series(sp.series(tot, dY, 0, 2).removeO(), ds, 0, 2).removeO()
print('  周回積分 = (∂_Y A_s - ∂_s A_Y) dY ds + (高次):', sp.simplify(sp.expand(lead) - (sp.diff(As, Y) - sp.diff(AY, s)) * dY * ds) == 0)

print('\n(b-ii) 円弧 r(X,Y,s) = (ρ+X)(cos θ, 0, sin θ) + Y (0,1,0)、θ = s/ρ: 小さい長方形の周回積分 / 面積 と n·rot A')
mp.mp.dps = 30
rho = mp.mpf(10)


def rvec(X, Yv, sv):
    th = sv / rho
    return mp.matrix([(rho + X) * mp.cos(th), Yv, (rho + X) * mp.sin(th)])


def Acart(P):
    x, y, z = P[0], P[1], P[2]
    return mp.matrix([y * z + x**2, x * z - y**3 + z * x, x * y + y * z**2 + x * z])


def rotA(P, h=mp.mpf('1e-10')):
    def d(i, j):
        Pp = P.copy()
        Pm = P.copy()
        Pp[j] += h
        Pm[j] -= h
        return (Acart(Pp)[i] - Acart(Pm)[i]) / (2 * h)
    return mp.matrix([d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)])


def line_int(path):
    """path(τ) (τ=0..1) に沿う ∫A·dr"""
    h = mp.mpf('1e-12')

    def integrand(tau):
        dr = (path(tau + h) - path(tau - h)) / (2 * h)
        A = Acart(path(tau))
        return A[0] * dr[0] + A[1] * dr[1] + A[2] * dr[2]
    return mp.quad(integrand, [0, 1])


X0, Y0, s0 = mp.mpf('0.002'), mp.mpf('0.001'), mp.mpf('1.3')
th0 = s0 / rho
nvec = mp.matrix([mp.cos(th0), 0, mp.sin(th0)])
bvec = mp.matrix([0, 1, 0])
tvec = mp.matrix([-mp.sin(th0), 0, mp.cos(th0)])
rc = rotA(rvec(X0, Y0, s0))
Bn_ref = (nvec.T * rc)[0]
Bb_ref = (bvec.T * rc)[0]
Bt_ref = (tvec.T * rc)[0]


def loop(dq1, dq2, which):
    """which=0: (Y,s) 平面 (n 方向の rot)、1: (s,X) 平面 (b 方向)、2: (X,Y) 平面 (t 方向)。戻り値: 周回積分 / (面積)"""
    hs = 1 + X0 / rho
    if which == 0:
        P = lambda a, b: rvec(X0, Y0 + a, s0 + b)
        area = dq1 * hs * dq2
    elif which == 1:
        P = lambda a, b: rvec(X0 + b, Y0, s0 + a)
        area = dq1 * hs * dq2
    else:
        P = lambda a, b: rvec(X0 + a, Y0 + b, s0)
        area = dq1 * dq2
    corners = [(0, 0), (dq1, 0), (dq1, dq2), (0, dq2), (0, 0)]
    tot = mp.mpf(0)
    for (a0, b0), (a1, b1) in zip(corners[:-1], corners[1:]):
        path = lambda tau, a0=a0, b0=b0, a1=a1, b1=b1: P(a0 + (a1 - a0) * tau, b0 + (b1 - b0) * tau)
        tot += line_int(path)
    return tot / area


for dq in (mp.mpf('1e-3'), mp.mpf('1e-4')):
    vals = [loop(dq, dq, w) for w in (0, 1, 2)]
    print('  辺 %.0e: B_n = %.8f (rot A の n 成分 %.8f)、B_b = %.8f (%.8f)、B_t = %.8f (%.8f)' % (
        float(dq), float(vals[0]), float(Bn_ref), float(vals[1]), float(Bb_ref), float(vals[2]), float(Bt_ref)))
print('  (面は始点 (Y0, s0, X0) の角から辺 dq の長方形。誤差は dq に比例して小さくなる)')
print('  (n, b, t) は右手系: n × b = t:', [float(v) for v in (mp.matrix([nvec[1] * bvec[2] - nvec[2] * bvec[1], nvec[2] * bvec[0] - nvec[0] * bvec[2], nvec[0] * bvec[1] - nvec[1] * bvec[0]]) - tvec)])

# (c) ------------------------------------------------------------------
print('\n(c) 設計軌道の条件 (1a.61)')
x, px, py, dl, rho_s, B0, P0e = sp.symbols('x px py delta rho B0 e_over_P0')
a_s = -(x + x**2 / (2 * rho_s)) * B0 * P0e   # eB0/P0 を B0*P0e と書く
H = -a_s - (1 + x / rho_s) * sp.sqrt((1 + dl)**2 - px**2 - py**2)
dHdx = sp.diff(H, x).subs({x: 0, px: 0, py: 0, dl: 0})
print('  ∂H/∂x|設計 = %s (= eB0/P0 - 1/ρ)' % sp.simplify(dHdx))
print('  p_x\' = -∂H/∂x|設計 = %s → 0 とおいて eB0/P0 = 1/ρ' % sp.simplify(-dHdx))
tt = sp.Matrix([1, 0, 0])
bb = sp.Matrix([0, 1, 0])
print('  t = (1,0,0)、b = (0,1,0)、n = -b×t…: t × b = %s、n = b × t = %s (n = -(t × b))' % (list(tt.cross(bb)), list(bb.cross(tt))))

# (d) ------------------------------------------------------------------
print('\n(d) 分散の方程式 η\'\' + η/ρ² = 1/ρ')
sv, rho2, a, b = sp.symbols('s rho a b', positive=True)
eta1 = rho2 * (1 - sp.cos(sv / rho2))
eta2 = rho2
for name, eta in (('η = ρ(1 - cos(s/ρ))', eta1), ('η = ρ (定数)', eta2), ('一般解 ρ + a cos + b sin', rho2 + a * sp.cos(sv / rho2) + b * sp.sin(sv / rho2))):
    print('  %s: η\'\' + η/ρ² - 1/ρ = %s' % (name, sp.simplify(sp.diff(eta, sv, 2) + eta / rho2**2 - 1 / rho2)))
print('  η1(0) = %s、η1\'(0) = %s。η1 - ρ = -ρ cos(s/ρ) (同次方程式の解)' % (eta1.subs(sv, 0), sp.diff(eta1, sv).subs(sv, 0)))

# (e) ------------------------------------------------------------------
print('\n(e) §9.6: 落とした H3 = x (px^2 + py^2)/(2ρ(1+δ)) の寄与 (ρ=10 m、δ=0、x = a、px = a/2、py = 0)')
for a_ in (1e-3, 2e-3, 4e-3, 8e-3):
    xx, ppx = a_, a_ / 2
    dxp = xx * ppx / 10          # ∂H3/∂px = x px/(ρ(1+δ))
    dpxp = -ppx**2 / (2 * 10)    # -∂H3/∂x
    print('  a = %.0f mm: ∂H3/∂px = %.2e (x\' への寄与/m)、-∂H3/∂x = %.2e (px\' への寄与/m)、3 m で %.2e、%.2e' % (a_ * 1e3, dxp, dpxp, 3 * dxp, 3 * dpxp))
