# -*- coding: utf-8 -*-
r"""端部 (06_edge_fringe.md) の第 2 回監査で足した途中式の確認 (sympy / numpy / scipy)。

約束 (C_poisson_lib.md): exp(:f:) の流れは Δx = -∂f/∂p_x、Δp_x = +∂f/∂x、Δz = -∂f/∂p。f = -∫H_ham ds。

(a) E-01: 相互作用表示。H = H0 + εH1(s) で、1 次の厳密な応答 = 「自由な運動 M0(S) ∘ exp(ε:f:)」(f = -∫H1(M0(s)z0)ds)。
         2 自由度 (H0 = (px²+py²)/2、H1 = x² y h(s))。および M0 が正準変換のとき [M0 a, M0 b] = M0 [a, b]。
(b) E-02: (V·∇)V の y 成分 = 3c² y (x²-y²)²、t y = ½(V·∇)V_y。
(c) E-03: K64 の四極の写像 (x + V + ½DV と運動量の逆転置) は、ある f' = -p·U の厳密な exp(:f':)。U = V - (1/6) D²V + O(K⁴)。
(d) E-04: CAVI。スカラーポテンシャルの端の項、H_f、ΔP_x、ΔE = ∂H_f/∂t、z の更新、E_r の力の積分 (ODE)。
(e) E-05: MULT。複素強度 κ_c = ak1 + i ak2 の多極 Φ = Im(κ_c w^m)/m!、ttfrins の回転の向き、tblfri (垂直 = 水平を +90 度回したもの)。
(f) E-06: ソレノイドの端の流れから (x, y, p_x, p_y) を戻す式。
(g) E-08: edge_04 (b) の 2 つの外挿の差の定義と値。
"""
import math
import warnings
warnings.filterwarnings("ignore")
from types import SimpleNamespace
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad import track6 as t6

# ---------------------------------------------------------------- (a)
print('== (a) E-01: 相互作用表示 (一般の導出の確認) ==')
x, y, px, py, s, S, eps = sp.symbols('x y p_x p_y s S epsilon', real=True)
q = [x, y]
pp = [px, py]


def pb(f, g):
    return sum(sp.diff(f, q[i]) * sp.diff(g, pp[i]) - sp.diff(f, pp[i]) * sp.diff(g, q[i]) for i in range(2))


# (a1) M0 が正準変換 (自由な運動 = 時間 s の流れ) のとき、関数の合成は括弧を保つ: (M0 a)(z) = a(M0 z)
M0 = {x: x + px * s, y: y + py * s, px: px, py: py}


def comp(f):
    return sp.expand(f.subs(M0, simultaneous=True))


a_f = x**2 * py + y**3
b_f = y * px**2 + x * py
lhs = pb(comp(a_f), comp(b_f))
rhs = comp(pb(a_f, b_f))
print('  [M0 a, M0 b] - M0 [a, b] = %s (a = x² p_y + y³、b = y p_x² + x p_y)' % sp.simplify(lhs - rhs))

# (a2) 1 次の厳密な応答 と M0 ∘ exp(ε:f:) の 1 次
h = 1 + s + s**2
H1 = x**2 * y * h
z0 = [x, px, y, py]
H1s = comp(H1)                          # H1(M0(s) z0)
f = -sp.integrate(H1s, (s, 0, S))
# 厳密な 1 次: δx' = δp_x、δp_x' = -∂H1/∂x(自由軌道)、δy' = δp_y、δp_y' = -∂H1/∂y
dHdx = comp(sp.diff(H1, x))
dHdy = comp(sp.diff(H1, y))
t_ = sp.symbols('t_', real=True)
dpx_t = sp.integrate((-dHdx).subs(s, t_), (t_, 0, s))
dpy_t = sp.integrate((-dHdy).subs(s, t_), (t_, 0, s))
dx_S = sp.integrate((S - s) * (-dHdx), (s, 0, S))
dy_S = sp.integrate((S - s) * (-dHdy), (s, 0, S))
dpx_S = dpx_t.subs(s, S)
dpy_S = dpy_t.subs(s, S)
# M0(S) ∘ exp(ε:f:) の ε の 1 次: z_out = M0(S)(z0 + ε [f, z0])
fx = -sp.diff(f, px)       # Δx = -∂f/∂p_x
fpx = sp.diff(f, x)        # Δp_x = +∂f/∂x
fy = -sp.diff(f, py)
fpy = sp.diff(f, y)
out_x = fx + S * fpx       # x(S) の 1 次: Δx + S Δp_x
out_y = fy + S * fpy
print('  Δx(S) の厳密な 1 次 - (Δx_f + S Δp_f)   = %s' % sp.simplify(sp.expand(dx_S - out_x)))
print('  Δy(S) の厳密な 1 次 - (Δy_f + S Δp_yf)  = %s' % sp.simplify(sp.expand(dy_S - out_y)))
print('  Δp_x(S) の厳密な 1 次 - Δp_xf           = %s' % sp.simplify(sp.expand(dpx_S - fpx)))
print('  Δp_y(S) の厳密な 1 次 - Δp_yf           = %s' % sp.simplify(sp.expand(dpy_S - fpy)))
# 順序を逆 (exp(:f:) が後) にした f~ = -∫H1(M0(s - S)... は別の f (終点で評価): 違うことを示す
f_end = -sp.integrate(H1.subs({x: x + px * (s - S), y: y + py * (s - S)}, simultaneous=True), (s, 0, S))
print('  (参考) 終点に引き戻す f_end - f = %s (≠ 0。順序で f が変わる)' % sp.simplify(sp.expand(f_end - f)))

# ---------------------------------------------------------------- (b)
print('\n== (b) E-02: (V·∇)V の y 成分 ==')
c = sp.symbols('c', real=True)
Vx = c * x * (x**2 + 3 * y**2)
Vy = -c * y * (3 * x**2 + y**2)
DVx = Vx * sp.diff(Vx, x) + Vy * sp.diff(Vx, y)
DVy = Vx * sp.diff(Vy, x) + Vy * sp.diff(Vy, y)
print('  (V·∇)V_x = %s' % sp.factor(DVx))
print('  (V·∇)V_y = %s' % sp.factor(DVy))
print('  差 (V·∇)V_y - 3c² y (x²-y²)² = %s' % sp.simplify(DVy - 3 * c**2 * y * (x**2 - y**2)**2))
a_, b_ = c * x**2 * 3, c * y**2 * 3          # a = Kx²/(4p) = 3c x²、b = 3c y² (c = K/(12p))
t_q = (a_ - b_)**2 / 6
print('  t y - ½(V·∇)V_y = %s (t = (a-b)²/6)' % sp.simplify(t_q * y - DVy / 2))
print('  (4.29): y(1 - a - b/3 + t) の y 成分のずれ = V_y + ½(V·∇)V_y: 差 = %s' % sp.simplify(y * (-a_ - b_ / 3 + t_q) - (Vy + DVy / 2)))

# ---------------------------------------------------------------- (c)
print('\n== (c) E-03: K64 の写像 = exp(:f\':)、f\' = -p·U、U = V - (1/6)D²V + O(K⁴) ==')
e = sp.symbols('e', real=True)          # K の次数を数える (V → eV)
Vv = [e * Vx, e * Vy]


def Dop(F, U):
    return [sum(U[j] * sp.diff(Fi, [x, y][j]) for j in range(2)) for Fi in F]


def flow_series(U, order):
    """x + U + ½ D_U U + (1/6) D_U² U + ... (D_U = U·∇。ベクトル場 U の時間 1 の流れの位置)。e の order 次まで"""
    term = [x, y]
    tot = [x, y]
    for n in range(1, order + 1):
        term = [sp.expand(t / n) for t in Dop(term, U)]
        tot = [tot[i] + term[i] for i in range(2)]
    return [sp.series(sp.expand(ti), e, 0, order + 1).removeO() for ti in tot]


W = [x + Vv[0] + Dop(Vv, Vv)[0] / 2, y + Vv[1] + Dop(Vv, Vv)[1] / 2]      # K64 (= x + V + ½(V·∇)V)
D1 = Dop(Vv, Vv)
D2 = Dop(D1, Vv)
U = [Vv[i] - D2[i] / 6 for i in range(2)]
fl = flow_series(U, 4)
res3 = [sp.expand(fl[i] - W[i]) for i in range(2)]
print('  U = V - (1/6)(V·∇)²V の時間 1 の流れ - K64 の位置の写像 を e の冪で: ')
for i, nm in enumerate('xy'):
    pol = sp.Poly(res3[i], e)
    print('   %s 成分: e^1,e^2,e^3 の係数 = %s、%s、%s (0 のはず)、e^4 の係数は 0 でない: %s' % (
        nm, sp.simplify(pol.coeff_monomial(e)), sp.simplify(pol.coeff_monomial(e**2)), sp.simplify(pol.coeff_monomial(e**3)),
        sp.simplify(pol.coeff_monomial(e**4)) != 0))
print('  (V·∇)²V_x = %s' % sp.factor(D2[0] / e**3))
print('  (V·∇)²V_y = %s' % sp.factor(D2[1] / e**3))
# 運動量の持ち上げ: ハミルトニアン H = p·U の流れで p(1) = J^{-T} p0 (dJ/dτ = (∂U)J から)
Jm = sp.Matrix(2, 2, lambda i, j: sp.diff(W[i], [x, y][j]))
print('  運動量の持ち上げ: d(J^{-T})/dτ = -(∂U)^T J^{-T} は、d J/dτ = (∂U) J から (一般の恒等式。ここでは確認しない)')

# ---------------------------------------------------------------- (d)
print('\n== (d) E-04: CAVI の端 ==')
r, sg, t, E0, Edot, g = sp.symbols('r s t E E_dot g', real=True)
G_ = sp.Function('G')(sg)
gs = sp.diff(G_, sg)
Et = sp.Function('E')(t)
# スカラーポテンシャル: φ = -E(t) G(s) + (r²/4) E g'(s) + O(r⁴)、g = G'
phi = -Et * G_ + r**2 / 4 * Et * sp.diff(gs, sg)
Er = -sp.diff(phi, r)
Ez = -sp.diff(phi, sg)
print('  φ = -E G + (r²/4) E g\'、E_r = -∂φ/∂r = %s' % sp.simplify(Er))
print('  E_z = -∂φ/∂s = %s (= E g - (r²/4) E g\'\')' % sp.simplify(Ez))
div = sp.simplify(sp.diff(r * Er, r) / r + sp.diff(Ez, sg))
print('  ∇·E = %s (= -(r²/4) E g\'\'\'。r^0 の項は 0。r² の項は φ の r⁴ の項で打ち消す)' % div)
# H = -sqrt((E_k - eφ)² - m² - p⊥²) (c = 1)。∂H/∂φ|_0 = e E/p_z = e/β
Ek, m_, pz, ee, ph_ = sp.symbols('E_k m p_z e phi_', positive=True)
Hh = -sp.sqrt((Ek - ee * ph_)**2 - m_**2 - pz**2 + pz**2)   # p_z² = E² - m² - p⊥² は φ=0 のとき
dHdphi = sp.diff(Hh, ph_).subs(ph_, 0)
print('  ∂H/∂φ (φ = 0) = %s (= e E_k/p_z = e/β_z)' % sp.simplify(dHdphi))
# 端 g' = δ(s - s0) の項: H_f = (e E(t) r²/(4 β)) δ(s - s0)。ここで e E(t) = -(eV/L)(sin(ωt - dφ) + sin dφ - offset)
om, V, L, dphi, off, beta = sp.symbols('omega V L dphi offset beta', positive=True)
xx, yy, tt = sp.symbols('x y t', real=True)
eE = -(V / L) * (sp.sin(om * tt - dphi) + sp.sin(dphi) - off)       # e = 1
Hf = eE * (xx**2 + yy**2) / (4 * beta)                               # δ の係数
print('  H_f = %s (δ の係数)。β = 1 でヘルプの (6.13)' % sp.simplify(Hf))
print('  Δp_x = -∂H_f/∂x = %s' % sp.simplify(-sp.diff(Hf, xx)))
print('  ΔE   = +∂H_f/∂t = %s' % sp.simplify(sp.diff(Hf, tt)))
# 力の積分: ΔP_x = ∫ E_r dt (e = 1) = -(x/2) E / β (g: 0 → 1)
print('  力の積分: ΔP_x = ∫ E_r ds/(β) = -(x/2)(E/β)(g(+) - g(-)) = -(x/2)E/β → -∂H_f/∂x と同じ')
# z の更新: t を保ち E を dE だけ変える。z = -β t (c=1)。β = p/h
h1, h2, p1, p2, z1 = sp.symbols('h1 h2 p1 p2 z1', positive=True)
tdef = -z1 * h1 / p1                  # K64 の定義 t = -z h/p (c = 1)
z2 = sp.simplify(-tdef * p2 / h2)      # t 一定で p/h だけ変わる: z = -t p/h
print('  z の更新: t を保つと z_2 = -t p_2/h_2 = %s (K64 の z → z p2/p1 · h1/h2。z_1 = %s の比 %s)' % (z2, z1, sp.simplify(z2 / z1)))

# ODE: 静的な E_z = E g、E_r = -(r/2) E g' + (r³/16) E g''' の中を走る粒子 (m = c = e = 1)。p_x の跳び / (-E x/(2β)) を β = 0.3, 0.6, 0.95 で
print('  ODE (静的な場。m = c = e = 1、E = 1e-3、d = 1e-3、x = 1e-4):')
Ec = 1e-3
dd = 1e-3
x_in = 1e-4


def gfun(sv):
    return 0.5 * (1 + math.tanh(sv / dd))


def g1(sv):
    return 0.5 / dd / math.cosh(sv / dd)**2


def g2(sv):
    u = sv / dd
    return -math.tanh(u) / dd**2 / math.cosh(u)**2


def g3(sv):
    u = sv / dd
    ch2 = 1 / math.cosh(u)**2
    return (-ch2 + 2 * math.tanh(u)**2 * ch2) / dd**3 * 1.0 * 1.0  # d/du(-tanh sech²) = -sech⁴ + 2 tanh² sech²


def rhs(tau, v):
    X, Y, Z, Px, Py, Pz = v
    en = math.sqrt(1 + Px * Px + Py * Py + Pz * Pz)
    r2 = X * X + Y * Y
    Er_over_r = -0.5 * Ec * g1(Z) + r2 / 16 * Ec * g3(Z)
    Ezv = Ec * gfun(Z) - r2 / 4 * Ec * g2(Z)
    return [Px / en, Py / en, Pz / en, X * Er_over_r, Y * Er_over_r, Ezv]


for bet in (0.3, 0.6, 0.95):
    gam = 1 / math.sqrt(1 - bet**2)
    pz0 = gam * bet
    s0 = -12 * dd
    # 速度 β の粒子を s0 から s1 まで。(τ は時間)
    T = 24 * dd / bet
    sol = solve_ivp(rhs, (0, T), [x_in, 0.0, s0, 0.0, 0.0, pz0], rtol=1e-12, atol=1e-18, method='DOP853')
    dpx = sol.y[3, -1]
    # 本体の収束 (一様な E_z の中の E_r = 0 なので p_x は端でだけ変わる)。kick = dpx
    print('   β = %.2f: Δp_x / (-E x/(2β)) = %.5f' % (bet, dpx / (-Ec * x_in / (2 * bet))))

# ---------------------------------------------------------------- (e)
print('\n== (e) E-05: MULT の複素強度・ttfrins の回転・tblfri ==')
w_ = sp.symbols('w', complex=True)
xr, yr = sp.symbols('x y', real=True)
a1, a2, kp = sp.symbols('ak1 ak2 kappa_prime', real=True)
for m in (2, 3, 4):
    wc = xr + sp.I * yr
    wb = xr - sp.I * yr
    kc = a1 + sp.I * a2
    kcb = a1 - sp.I * a2
    a_w = kp / (4 * sp.factorial(m)) * wc * (kc * wc**m / (m + 1) - kcb * wb**m)     # κ_c' で a_w (実線形の拡張)
    ax = sp.re(sp.expand(a_w))
    ay = sp.im(sp.expand(a_w))
    bs = sp.simplify(sp.diff(ay, xr) - sp.diff(ax, yr))
    target = kp * sp.im(sp.expand(kc * wc**m)) / sp.factorial(m)
    print('  m = %d: b_s = ∂_x a_y - ∂_y a_x - κ\' Im(κ_c w^m)/m! = %s' % (m, sp.simplify(sp.expand(bs - target))))
    # 横成分の κ'' 項 (Φ_2 = -κ_c'' ... 実部・虚部ごとに線形): 実線形なので ak1、ak2 それぞれで (6.2) と同じ構造
# ttfrins: 回転 w' = e^{iθ} w (x' = cosθ x - sinθ y、y' = sinθ x + cosθ y)、θ = arg(ak)·2/nord、
# 回転した系で正規成分の強さ |ak| を与え、戻す: Δw = e^{-iθ} Δw'(w')。Δw'(w') = -(|ak|/(4 m! p)) w'(w'^m/(m+1) - conj(w')^m)
m = 3
rng = np.random.default_rng(3)
maxd = 0.0
for trial in range(5):
    ak = complex(rng.normal(), rng.normal()) * 3.0
    al, dl = 1.7, rng.uniform(-0.02, 0.02)
    x0, y0, px0, py0 = rng.normal(size=4) * 2e-3
    qp = SimpleNamespace(c=np.array([x0, px0, y0, py0, 0.0, dl], float))
    t6._ttfrins(qp, 2 * m, ak, al)
    wcn = complex(x0, y0)
    pr = 1 + dl
    kc_ = ak / al
    dw = -(1 / (4 * math.factorial(m) * pr)) * wcn * (kc_ * wcn**m / (m + 1) - kc_.conjugate() * wcn.conjugate()**m)
    maxd = max(maxd, abs(complex(qp.c[0], qp.c[2]) - (wcn + dw)))
print('  pysad の _ttfrins (m = 3、複素 ak) の位置 - [w + Δw、Δw = -(1/(4 m! p)) w (κ_c w^m/(m+1) - conj(κ_c) conj(w)^m)] の最大差 = %.2e' % maxd)
print('   (κ_c = ak/al = ak1 + i ak2 は Φ = Im(κ_c w^m)/m! の複素強度)。差が 1e-18 台なら回転の向き (+θ で入れ -θ で戻す) と複素強度の規約が一致')
# 逆向き (-θ で入れる) だと一致しないことの確認
maxd2 = 0.0
for trial in range(5):
    ak = complex(rng.normal(), rng.normal()) * 3.0
    al = 1.7
    x0, y0 = rng.normal(size=2) * 2e-3
    wcn = complex(x0, y0)
    kc_ = ak / al
    th = math.atan2(ak.imag, ak.real) * 2.0 / (2 * m)
    wr = wcn * np.exp(-1j * th)       # 逆向き
    dwr = -(abs(kc_) / (4 * math.factorial(m))) * wr * (wr**m / (m + 1) - wr.conjugate()**m)
    wrev = (wr + dwr) * np.exp(1j * th)
    wfwd = wcn + (-(1 / (4 * math.factorial(m))) * wcn * (kc_ * wcn**m / (m + 1) - kc_.conjugate() * wcn.conjugate()**m))
    maxd2 = max(maxd2, abs(wrev - wfwd))
print('  (参考) 回転を逆向きに取ったときの同じ量との差 = %.2e (0 でない: 向きが決まる)' % maxd2)

# tblfri: 垂直の曲がり (ck0.imag) = 水平の式を +90 度回したもの
def tbl(c0, ck0, fb1, al):
    qp = SimpleNamespace(c=np.array(c0, float))
    t6._tblfri(qp, al, ck0, fb1)
    return qp.c.copy()


c0 = [1.3e-3, 2.1e-3, -0.8e-3, 1.1e-3, 0.0, 0.03]
al, fb1, k0 = 1.0, 0.004, 0.02
ch = tbl(c0, complex(k0, 0.0), fb1, al)                       # 水平の曲がり
# 回転 (x', px', y', py') = (-y, -py, x, px) で入れて水平の式 → 戻す (x, px, y, py) = (y', py', -x', -px')
cr = [-c0[2], -c0[3], c0[0], c0[1], c0[4], c0[5]]
chr_ = tbl(cr, complex(k0, 0.0), fb1, al)
back = [chr_[2], chr_[3], -chr_[0], -chr_[1], chr_[4], chr_[5]]
cv = tbl(c0, complex(0.0, k0), fb1, al)                       # 垂直の曲がり
print('  tblfri: 垂直 (ck0 = i k0) - [水平を +90° 回して入れ戻したもの] の最大差 = %.2e' % np.abs(cv - np.array(back)).max())
# (5.42) の f_lin を回した f から Δ を出す: f = p_y dxfr (1 - 1/p) + (1/p)(dyfr x²/2 - dyfra x⁴/4)
xs_, pxs_, ys_, pys_, ps_ = sp.symbols('x p_x y p_y p', real=True)
dxfr, dyfr, dyfra = sp.symbols('dxfr dyfr dyfra', real=True)
fH = -pxs_ * dxfr * (1 - 1 / ps_) + (dyfr * ys_**2 / 2 - dyfra * ys_**4 / 4) / ps_
# 回転: x' = -y、y' = x、p_x' = -p_y、p_y' = p_x
fV = fH.subs({xs_: -ys_, ys_: xs_, pxs_: -pys_}, simultaneous=True)
print('  f_V = f_H(x\'=-y, y\'=x, p_x\'=-p_y) = %s' % sp.simplify(fV))
print('   Δy = -∂f_V/∂p_y = %s、Δp_x = ∂f_V/∂x = %s' % (sp.simplify(-sp.diff(fV, pys_)), sp.simplify(sp.diff(fV, xs_))))
print('   Δx = -∂f_V/∂p_x = %s (0)、Δz = -∂f_V/∂p = %s' % (sp.simplify(-sp.diff(fV, pxs_)), sp.simplify(-sp.diff(fV, ps_))))

# ---------------------------------------------------------------- (f)
print('\n== (f) E-06: ソレノイドの端の流れから (x, y, p_x, p_y) を戻す ==')
X, Y, PX, PY, k = sp.symbols('X Y P_x P_y k', real=True)
R2 = X**2 + Y**2
PU = X * PX + Y * PY
PPHI = X * PY - Y * PX
print('  p_x - (x p_u - y p_φ)/r² = %s' % sp.simplify(PX - (X * PU - Y * PPHI) / R2))
print('  p_y - (y p_u + x p_φ)/r² = %s' % sp.simplify(PY - (Y * PU + X * PPHI) / R2))
# 数値: tsfrin と、この戻し方
b0 = 1.7
c0 = [0.003, 0.02, -0.002, 0.015, 0.0, 0.01]
qp = SimpleNamespace(c=np.array(c0, float))
t6._tsfrin(qp, b0)
pr_ = 1 + c0[5]
kk = b0 / (8 * pr_)        # q 単位での k (= b0/(8p²) を p で割ったもの = b0/(8p)): Δu = -k p_φ(q)
x_, px_, y_, py_ = c0[0], c0[1], c0[2], c0[3]
pu = x_ * px_ + y_ * py_
pphi = x_ * py_ - y_ * px_
r0 = math.hypot(x_, y_)
ph0 = math.atan2(y_, x_)
r1 = r0 * math.exp(-kk * pphi)
ph1 = ph0 - kk * pu
x1, y1 = r1 * math.cos(ph1), r1 * math.sin(ph1)
px1 = (x1 * pu - y1 * pphi) / r1**2
py1 = (y1 * pu + x1 * pphi) / r1**2
z1 = c0[4] + 2 * kk * pphi * pu
mine = np.array([x1, px1, y1, py1, z1])
print('  この戻し方 (r1 = r e^{-k p_φ}、φ1 = φ - k p_u、p_x = (x p_u - y p_φ)/r² …) - tsfrin = %.2e (最大)' % np.abs(mine - qp.c[:5]).max())

# ---------------------------------------------------------------- (g)
print('\n== (g) E-08: edge_04 (b) の 2 つの外挿の差 ==')
ode1 = np.array([-5.7554e-11, 5.9074e-10, 9.2178e-12, 1.1405e-10])
ode2 = np.array([-5.7605e-11, 5.6647e-10, 1.0096e-11, 8.5775e-11])
tsf_ = np.array([-5.7624e-11, 5.8153e-10, 9.4860e-12, 1.0275e-10])
print('  成分 (x, px, y, py)')
print('  tsfrin に対する 外挿 1 (d = 0.5/1 mm) の比   = %s' % np.array2string(ode1 / tsf_, precision=4))
print('  tsfrin に対する 外挿 2 (d = 0.25/0.5 mm) の比 = %s' % np.array2string(ode2 / tsf_, precision=4))
print('  |外挿 1 - 外挿 2| / |tsfrin|                  = %s' % np.array2string(np.abs(ode1 - ode2) / np.abs(tsf_), precision=4))

# ---------------------------------------------------------------- (h)
print('\n== (h) E-07 の補足: 偏向の端の場 b_y + i b_s = g(s + i y) (6.1 の級数の和) ==')
s_, y_ = sp.symbols('s y', real=True)
gf = sp.Function('g')
# 実部・虚部を y の冪で取り出す (g は実関数なので Taylor 展開の係数は g^(n)(s))
re_part = sp.expand(sum(sp.Derivative(gf(s_), (s_, n)) * (sp.I * y_)**n / sp.factorial(n) for n in (0, 2, 4)))
im_part = sp.expand(sum(sp.Derivative(gf(s_), (s_, n)) * (sp.I * y_)**n / sp.factorial(n) for n in (1, 3, 5)) / sp.I)
print('  Re g(s+iy) = %s' % re_part)
print('  Im g(s+iy) = %s' % im_part)
print('  (6.1): b_y = g - y² g\'\'/2 + y⁴ g\'\'\'\'/24、b_s = y g\' - y³ g\'\'\'/6 + y⁵ g^(5)/120 と一致')
print('  tanh の端 g = (1 + tanh(ζ/d))/2 の極: tanh ζ/d の分母 cosh(ζ/d) = 0 → ζ = ±iπ d/2: |y| < π d/2 でだけ場は有限')
dd_ = sp.symbols('d', positive=True)
print('  cosh(ζ/d) at ζ = i π d/2: %s' % sp.simplify(sp.cosh(sp.I * sp.pi / 2)))
