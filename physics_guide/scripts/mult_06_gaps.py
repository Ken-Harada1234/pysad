# -*- coding: utf-8 -*-
r"""08_mult.md の第 2 回監査で足した途中式の確認 (sympy / numpy / scipy)。

(a) M-01: 刻み数 (8.8) の ε の読みと係数 1/6、Q の冪、実測が見積りの約 1.5 倍になる理由。
    Strang 分割の実効ハミルトニアン H~ = H + h²[(1/12)|∇V|² - (1/24) p_i p_j ∂_i∂_j V] (BCH、sympy)。H~ の厳密な流れが分割の追跡の誤差を再現する。
    調和振動子: 位相 (1/24) と運動量の振幅 (1/12) を合わせて、運動量の跳びの相対誤差 (ωh)²/6 = Q/(6 ndiv²)。
(b) M-02: dofr(n) の量 D_n = (K_n/L) a^(n+1)/(n+1) と、(7.17) の端の位置のずれ (Δw/a) との比 4 n! (n+2)。
(c) M-03: 端面の角のキック。「長さ × 力」だけの式の rot = -(t/L) Im C(z)、(8.9) との差 (i t y/L) C(z)。
(d) M-04: 出口の写像 (CHI3 → R1^-1 R2^-1 → 面 → +DX、+DY)。tsol.f の出口と同じであること。
(e) M-05: MULT が tsolrot に -CHI2 を渡す理由: tbrot(+CHI2) (ANGLE ≠ 0 の経路) = tsolrot(-CHI2) (ANGLE = 0 の経路)。
(f) M-06: ANGLE ≠ 0 の MULT の L' の二重補正。実験室の平面の ODE (線形の立ち上がりの端) と、MULT・BEND の追跡の比較。
"""
import io
import math
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)

# ---------------------------------------------------------------- (a)
print('== (a) M-01: ndiv の式 (8.8) の ε の読み ==')
x, y, px, py, h, K, L = sp.symbols('x y p_x p_y h K L', real=True)
V = K / (6 * L) * (x**3 - 3 * x * y**2)          # K2 (n = 2) の力 = -∇V、∇V = (K/(2L)(x² - y²), -K x y / L)
T = (px**2 + py**2) / 2
qs, ps = [x, y], [px, py]


def pb(f, g):
    return sum(sp.diff(f, qs[i]) * sp.diff(g, ps[i]) - sp.diff(f, ps[i]) * sp.diff(g, qs[i]) for i in range(2))


fX, fY = -T, -V                                   # f = -H (C 付録の約束)
ferr = h**3 * (-sp.Rational(1, 24) * pb(fX, pb(fX, fY)) + sp.Rational(1, 12) * pb(fY, pb(fY, fX)))
dH = sp.simplify(-ferr / h)                       # 1 刻み h あたり f_err = -h δH
claim = h**2 * ((sp.diff(V, x)**2 + sp.diff(V, y)**2) / 12
                - (px**2 * sp.diff(V, x, 2) + 2 * px * py * sp.diff(V, x, y) + py**2 * sp.diff(V, y, 2)) / 24)
print('  δH (BCH: -1/24 [X,[X,Y]] + 1/12 [Y,[Y,X]] から) - h²[|∇V|²/12 - p_i p_j V_ij/24] = %s' % sp.simplify(dH - claim))
w2, hh = sp.symbols('omega2 hh', positive=True)
Vh = w2 * x**2 / 2
dHh = (sp.diff(Vh, x)**2 / 12 - px**2 * sp.diff(Vh, x, 2) / 24) * hh**2
Hh = sp.expand(px**2 / 2 + Vh + dHh)
A = sp.simplify(Hh.coeff(px, 2) * 2)
B = sp.simplify(Hh.coeff(x, 2) * 2)
print('  調和振動子: H~ = A p²/2 + B x²/2、A = %s、B = %s' % (A, B))
print('   Ω² = A B = ω²(1 + ω²h²/12 + ...)、Ω/ω = 1 + ω²h²/24: %s' % sp.series(sp.sqrt(A * B) / sp.sqrt(w2), hh, 0, 3))
print('   運動量の振幅: p = -x0 Ω sin(Ω t)/A、Ω/A の ω に対する比 = %s' % sp.series(sp.sqrt(A * B) / A / sp.sqrt(w2), hh, 0, 3))
print('   → ωL ≪ 1 の運動量の跳び (≈ -x0 ω² L) の相対誤差 = Ω/A の変化 ω²h²/8 (= 1/24 + 1/12) と、sin(Ωt) ≈ Ωt の位相の変化 ω²h²/24 の和 = (ωh)²/6')

# 追跡との比較
def force(x_, y_, K_, L_):
    z = complex(x_, y_)
    g = K_ / L_ * z**2 / 2
    return g.real, -g.imag


def exact(c, L_, K_):
    def rhs(s, v):
        xx, pxx, yy, pyy = v
        pz = math.sqrt(1 - pxx * pxx - pyy * pyy)
        gx, gy = force(xx, yy, K_, L_)
        return [pxx / pz, -gx, pyy / pz, -gy]
    return solve_ivp(rhs, (0, L_), list(c), method='DOP853', rtol=1e-13, atol=1e-17).y[:, -1]


def drift(c, hv):
    xx, pxx, yy, pyy = c
    pz = math.sqrt(1 - pxx * pxx - pyy * pyy)
    return np.array([xx + pxx / pz * hv, pxx, yy + pyy / pz * hv, pyy])


def kick(c, K_, L_, wv):
    xx, pxx, yy, pyy = c
    gx, gy = force(xx, yy, K_, L_)
    return np.array([xx, pxx - gx * wv * L_, yy, pyy - gy * wv * L_])


def strang(c, L_, K_, n):
    hv = L_ / n
    c = drift(np.array(c, float), hv / 2)
    for i in range(n):
        c = kick(c, K_, L_, 1.0 / n)
        c = drift(c, hv if i < n - 1 else hv / 2)
    return c


dHn = claim
Hf = sp.lambdify((x, y, px, py, h, K, L), [sp.diff(T + V + dHn, px), -sp.diff(T + V + dHn, x), sp.diff(T + V + dHn, py), -sp.diff(T + V + dHn, y)], 'numpy')


def flow_mod(c, L_, K_, n):
    hv = L_ / n
    f = lambda s, v: Hf(v[0], v[2], v[1], v[3], hv, K_, L_)
    return solve_ivp(f, (0, L_), list(c), method='DOP853', rtol=1e-13, atol=1e-18).y[:, -1]


Lv, Kv = 1.0, 40.0
c0 = [0.004, 0.0008, -0.003, -0.0004]
print('  K2 = 40 /m²、L = 1 m、初期値 (4 mm, 0.8 mrad, -3 mm, -0.4 mrad)。(x, px, y, py) の誤差:')
for n in (5, 17):
    ex = exact(c0, Lv, Kv)
    st = strang(c0, Lv, Kv, n) - ex
    md = flow_mod(c0, Lv, Kv, n) - ex
    print('   ndiv = %2d: Strang の追跡 - 厳密 = %s' % (n, np.array2string(st, precision=3)))
    print('            H~ の厳密な流れ - 厳密 = %s (差の最大 %.1e)' % (np.array2string(md, precision=3), np.abs(st - md).max()))
# 運動量の誤差の予想 Δp = -(h² L/6) (∇∇V)(∇V) を経路の中点で評価
for n in (5, 17):
    hv = Lv / n
    xm, ym = c0[0] + c0[1] * Lv / 2, c0[2] + c0[3] * Lv / 2
    gx, gy = Kv / (2 * Lv) * (xm**2 - ym**2), -Kv / Lv * xm * ym
    Hxx, Hxy, Hyy = Kv * xm / Lv, -Kv * ym / Lv, -Kv * xm / Lv
    dpx = -(hv**2 * Lv / 6) * (Hxx * gx + Hxy * gy)
    dpy = -(hv**2 * Lv / 6) * (Hxy * gx + Hyy * gy)
    st = strang(c0, Lv, Kv, n) - exact(c0, Lv, Kv)
    print('   ndiv = %2d: Δp の予想 (-(h²L/6)(∇∇V)(∇V)、経路の中点) = (%.3e, %.3e)、測定 (px, py) = (%.3e, %.3e)' % (n, dpx, dpy, st[1], st[3]))
# 調和の見積りとの比
Q = 0.16
a_ = 0.004
for n in (5, 17):
    est_phase = a_ * Q**1.5 / (24 * n**2)
    est_kick = a_ * Q**2 / (6 * n**2)
    st = np.abs(strang(c0, Lv, Kv, n) - exact(c0, Lv, Kv)).max()
    print('   ndiv = %2d: 位相の見積り a (ωL)³/(24 n²) = %.3e (測定/見積り = %.2f)、運動量の跳びの見積り a Q²/(6 n²) = %.3e (測定/見積り = %.2f)、'
          '見積りの比 4 ωL = %.2f' % (n, est_phase, st / est_phase, est_kick, st / est_kick, 4 * math.sqrt(Q)))
print('  (8.8) の ndiv の下で、運動量の跳びの相対誤差 (ωh)²/6 = Q/(6 ndiv²) が ε になる: 例 Q = 2、EPS = 1 (ε = 0.005): ndiv = %d、Q/(6 ndiv²) = %.5f' % (
    1 + int(math.sqrt(2 / (6 * 0.005))), 2 / (6 * (1 + int(math.sqrt(2 / (6 * 0.005))))**2)))
print('  (参考) (8.8) を「全位相誤差 Q^{3/2}/(24 ndiv²) ≤ ε」と読むと ndiv² = Q^{3/2}/(24 ε)。(8.8) の ndiv² との比 = 4/√Q = %.2f (Q = 2)' % (4 / math.sqrt(2)))

# ---------------------------------------------------------------- (b)
print('\n== (b) M-02: dofr(n) の判定量 ==')
nn = sp.symbols('n', positive=True, integer=True)
Dn = sp.Symbol('K_L_a')                                       # (K_n/L) a^(n+1)
dw = Dn / (4 * sp.factorial(nn + 1) * (nn + 2))              # |Δw|/a の第 1 項: κ a^(m+1)/(4 m! (m+1)) / a、m = n + 1
crit = Dn / (nn + 1)
print('  D_n = (K_n/L) a^(n+1)/(n+1) と Δw/a (第 1 項) の比 = %s' % sp.simplify(crit / dw))
for n_ in (2, 3, 4, 5, 6):
    print('   n = %d: 4 n!(n+2) = %d、判定 D_n > ε は Δw/a > ε/%d (ε = 0.005 で %.1e)' % (n_, 4 * math.factorial(n_) * (n_ + 2), 4 * math.factorial(n_) * (n_ + 2),
                                                                                     0.005 / (4 * math.factorial(n_) * (n_ + 2))))

# ---------------------------------------------------------------- (c)
print('\n== (c) M-03: 「長さ × 力」だけの式の rot ==')
xr, yr, t_, L_s = sp.symbols('x y t L', real=True)
zc = xr + sp.I * yr
Cf = sp.Function('C')
for nm, Cz in (('C = K0 + K1 z', sp.Symbol('K0', real=True) + sp.Symbol('K1', real=True) * zc),
               ('C = K0 + K1 z + K2 z²/2', sp.Symbol('K0', real=True) + sp.Symbol('K1', real=True) * zc + sp.Symbol('K2', real=True) * zc**2 / 2)):
    dpx = t_ * xr / L_s * sp.re(sp.expand(Cz))
    dpy = -t_ * xr / L_s * sp.im(sp.expand(Cz))
    rot = sp.simplify(sp.diff(dpy, xr) - sp.diff(dpx, yr))
    claim_rot = sp.simplify(-t_ / L_s * sp.im(sp.expand(Cz)))
    div = sp.simplify(sp.diff(dpx, xr) + sp.diff(dpy, yr))
    print('  %s: rot = %s、-(t/L) Im C(z) との差 = %s、div = %s' % (nm, rot, sp.simplify(rot - claim_rot), div))
    full = t_ / L_s * zc * Cz                 # Δp_x - i Δp_y (8.9)
    ans = t_ * xr / L_s * Cz
    diffc = sp.simplify(sp.expand(full - ans - sp.I * t_ * yr / L_s * Cz))
    print('     (8.9) - 「長さ×力」 - (i t y/L) C = %s' % diffc)
K1s = sp.Symbol('K1', real=True)
print('  K1 だけ: (i t y/L) K1 z の実部・虚部 = %s、%s (Δp_x に -t K1 y²/L、Δp_y に -t K1 x y/L を足す)' % (
    sp.simplify(sp.re(sp.expand(sp.I * t_ * yr / L_s * K1s * zc))), sp.simplify(-sp.im(sp.expand(sp.I * t_ * yr / L_s * K1s * zc)))))

# ---------------------------------------------------------------- (d)
print('\n== (d) M-04: 出口の写像 ==')


def Rm(chi1, chi2):
    c1, s1, c2, s2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    R1 = np.array([[c1, 0, -s1], [0, 1, 0], [s1, 0, c1]])
    R2 = np.array([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
    return R1, R2


def ent(vec, alg, dx, dy, dz, chi1, chi2):
    X_, px_, Y_, py_, zc_ = vec
    s0_ = -alg
    C1, S1, C2, S2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    X_ -= dx
    Y_ -= dy
    p0 = math.sqrt(1 - px_ * px_ - py_ * py_)
    q1 = S1 * px_ + C1 * p0
    px_ = C1 * px_ - S1 * p0
    py_ = C2 * py_ - S2 * q1
    d1 = S1 * X_ + (C1 - 1) * s0_ - C1 * dz
    X_ = C1 * X_ - S1 * (s0_ - dz)
    d2 = S2 * Y_ + C2 * d1 + (C2 - 1) * s0_
    Y_ = C2 * Y_ - S2 * (d1 + s0_)
    q2 = math.sqrt(1 - px_ * px_ - py_ * py_)
    aa = d2 / q2
    return [X_ - aa * px_, px_, Y_ - aa * py_, py_, zc_ + aa]


def ext(vec, alg, dx, dy, dz, chi1, chi2):          # K64 の tsolrot の出口の枝を書き写したもの (引数 alg に alg - al を渡す)
    X_, px_, Y_, py_, zc_ = vec
    s0_ = alg
    C1, S1, C2, S2 = math.cos(chi1), math.sin(chi1), math.cos(chi2), math.sin(chi2)
    p0 = math.sqrt(1 - px_ * px_ - py_ * py_)
    q1 = -S2 * py_ + C2 * p0
    pyi = C2 * py_ + S2 * p0
    q2 = -S1 * px_ + C1 * q1
    pxi = C1 * px_ + S1 * q1
    y1 = C2 * Y_ - S2 * s0_
    d1 = -S2 * Y_ - (C2 - 1) * s0_
    d2 = -S1 * X_ + C1 * d1 - (C1 - 1) * s0_ + dz
    x1 = C1 * X_ + S1 * (d1 - s0_)
    aa = d2 / q2
    return [x1 - aa * pxi + dx, pxi, y1 - aa * pyi + dy, pyi, zc_ + aa]


def ext_matrix(vec, Lel, dx, dy, dz, chi1, chi2):
    """出口を (8.10) の逆変換として 3 行で書いたもの。要素の中心を原点とし、出口の面の第 3 成分は +L/2。"""
    X_, px_, Y_, py_, zc_ = vec
    R1, R2 = Rm(chi1, chi2)
    Rinv = R1.T @ R2.T                                # (R2 R1)^-1 = R1^-1 R2^-1 (CHI3 は 0)
    pz_ = math.sqrt(1 - px_ * px_ - py_ * py_)
    v = Rinv @ np.array([X_, Y_, Lel / 2])            # 行 1: 位置・運動量を逆回転
    p = Rinv @ np.array([px_, py_, pz_])
    a = (v[2] + dz - Lel / 2) / p[2]                  # 行 2: 実験室の面 s = L/2 まで直進 (後ろへ a)、z は +a
    return [v[0] - a * p[0] / 1.0 + dx, p[0], v[1] - a * p[1] / 1.0 + dy, p[1], zc_ + a]   # 行 3: +DX、+DY


rng = np.random.default_rng(5)
maxd = 0.0
for trial in range(6):
    dx_, dy_, dz_ = rng.normal(size=3) * 2e-3
    c1v, c2v = rng.normal(size=2) * 0.03
    vec = list(rng.normal(size=4) * 2e-3) + [0.0]
    vv = [vec[0], vec[1], vec[2], vec[3], vec[4]]
    a1 = ext(vv, -1.0 / 2, dx_, dy_, dz_, c1v, c2v)          # Lel = 1: alg - al = -1/2
    a2 = ext_matrix(vv, 1.0, dx_, dy_, dz_, c1v, c2v)
    maxd = max(maxd, max(abs(a1[i] - a2[i]) for i in range(5)))
print('  tsol.f の出口 (書き写し) と、3 行 (R1^-1 R2^-1、面まで a、+DX +DY) の最大差 = %.2e (x, px, y, py, z。6 組のランダムな DX..CHI2)' % maxd)
print('  z の符号: 点を a だけ後ろへ戻す (x → x - a p_x/p_z) と z → z + a。z = -v(t - t0) なので時刻を早める (遅い点から面に戻す) と z が増える。入口も出口も同じ式')

# ---------------------------------------------------------------- (e)
print('\n== (e) M-05: tsolrot(-CHI2) = tbrot(+CHI2) (phig = 0、dtheta = 0、入口) ==')


def tbrot_ent(vec, alg, dchi2):
    """tbend.f の tbrot (ent = .true.、phig = 0、dtheta = 0): c0 = 1、s0 = 0、st = 0、ct = 1"""
    X_, px_, Y_, py_, zc_ = vec
    s2, c2 = math.sin(dchi2), math.cos(dchi2)
    dr11, dr22, dr33 = 0.0, c2 - 1.0 + 0.0, c2 - 1.0
    # dr11 = dct c0² + dc2 s0² + ... = 0、dr22 = c2 dct + dc2 = c2 - 1 (dct = 0)、dr33 = dc2 c0² + ... = c2 - 1
    r11, r22, r33 = 1.0, c2, c2
    r12, r13, r21, r23, r31, r32 = 0.0, 0.0, 0.0, s2, 0.0, -s2
    pzi = math.sqrt(1 - px_ * px_ - py_ * py_)
    dl = alg
    zi = -dl
    dyz = r23 * zi
    dzz = dr33 * zi
    xi = X_
    yi = Y_
    xf = dr11 * xi + r12 * yi + 0.0 + X_
    yf = r21 * xi + dr22 * yi + dyz + yi
    zf = r31 * xi + r32 * yi + dzz
    pxf = r11 * px_ + r12 * py_ + r13 * pzi
    pyf = r21 * px_ + r22 * py_ + r23 * pzi
    pzf = r31 * px_ + r32 * py_ + r33 * pzi
    return [xf - pxf / pzf * zf, pxf, yf - pyf / pzf * zf, pyf, zc_ + zf / pzf]


maxd = 0.0
for trial in range(6):
    chi2 = rng.normal() * 0.03
    vec = list(rng.normal(size=4) * 2e-3) + [0.0]
    a1 = tbrot_ent(vec, 0.5, chi2)
    a2 = ent(vec, 0.5, 0.0, 0.0, 0.0, 0.0, -chi2)
    maxd = max(maxd, max(abs(a1[i] - a2[i]) for i in range(5)))
print('  tbrot(+CHI2) と tsolrot(-CHI2) の最大差 = %.2e (6 組)' % maxd)
maxd2 = 0.0
for trial in range(6):
    chi2 = rng.normal() * 0.03
    vec = list(rng.normal(size=4) * 2e-3) + [0.0]
    a1 = tbrot_ent(vec, 0.5, chi2)
    a2 = ent(vec, 0.5, 0.0, 0.0, 0.0, 0.0, +chi2)
    maxd2 = max(maxd2, max(abs(a1[i] - a2[i]) for i in range(5)))
print('  (参考) tbrot(+CHI2) と tsolrot(+CHI2) の最大差 = %.2e (0 でない: -CHI2 が要る)' % maxd2)

# ---------------------------------------------------------------- (f)
print('\n== (f) M-06: ANGLE ≠ 0 の MULT の L\' の二重補正 (実験室の平面 ODE との比較) ==')
from pysad.ffscmd import run_sad_string
from pysad import track6n_ext
from pysad import compat
compat.set_mode('k64')      # この比較は K64 の動作 (B331 の二重補正) を測る。v2 の結果は末尾

L0, phi, ell = 1.0, 0.2, 0.2
c1_ = ((phi * ell)**2 * 2) / L0 / 48.0
Lp = L0 - c1_
rho_p = Lp / phi
rho0 = L0 / phi
Ed = 2 * rho0 * math.sin(phi / 2) * np.array([-math.sin(phi / 2), math.cos(phi / 2)])
d_out = np.array([-math.sin(phi), math.cos(phi)])
x_out = np.array([math.cos(phi), math.sin(phi)])


def g_in(s):
    return min(1.0, max(0.0, s / ell + 0.5))


def bfield(X, Z):
    s_in = Z
    s_out = (X - Ed[0]) * d_out[0] + (Z - Ed[1]) * d_out[1]
    return g_in(s_in) * g_in(-s_out) / rho_p


def ode_exit(x0, px0, delta):
    p = 1 + delta
    pz = math.sqrt(p * p - px0 * px0)
    z_start = -ell
    X0 = x0 - px0 / pz * (0 - z_start)

    def rhs(tau, v):
        X, Z, PX, PZ = v
        b = bfield(X, Z)
        return [PX / p, PZ / p, -b * PZ / p, b * PX / p]

    def ev(tau, v):
        return (v[0] - Ed[0]) * d_out[0] + (v[1] - Ed[1]) * d_out[1] - ell
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 5.0), [X0, z_start, px0, pz], method='DOP853', rtol=1e-13, atol=1e-16, events=ev, max_step=0.01)
    X, Z, PX, PZ = sol.y[:, -1]
    s_out = (X - Ed[0]) * d_out[0] + (Z - Ed[1]) * d_out[1]
    P = np.array([PX, PZ])
    back = s_out / (P @ d_out) * p
    pos = np.array([X, Z]) - back * P / p
    return (pos - Ed) @ x_out, P @ x_out


src = ('MOMENTUM = 3GEV;\n%s MP=(L=%.15g %s);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n')
ex = 'ANGLE=%g E1=0 E2=0 FRINGE=3 FB1=%g FB2=%g DISFRIN=1' % (phi, ell, ell)


def track(kind, c):
    ffs, _ = run_sad_string(src % (kind, L0, ex), out=io.StringIO())
    a = np.zeros((1, 6))
    a[0, :] = c
    return np.asarray(track6n_ext.track_particles_ext(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


print('  条件: 扇形 (E1 = E2 = 0)、ANGLE = %g、L = %g m、FB1 = FB2 = %g m (線形の立ち上がり、FRINGE = 3)、平面 (y = 0)。' % (phi, L0, ell))
print('  磁場の円の半径 ρ\' = L\'/ANGLE (L\' = L - c1、c1 = %.4e)、設計の弦 E = 2 (L/ANGLE) sin(ANGLE/2)' % c1_)
r0 = ode_exit(0.0, 0.0, 0.0)
print('  設計粒子 (δ = 0) の ODE の出口: x = %.2e、px = %.2e (硬い端の幾何に L\' が合っている)' % r0)
res = {}
for delta in (0.01, -0.01):
    r = ode_exit(0.0, 0.0, delta)
    m = track('MULT', [0, 0, 0, 0, 0, delta])
    b = track('BEND', [0, 0, 0, 0, 0, delta])
    res[delta] = (r, m, b)
    print('  δ = %+.2f: ODE (x, px) = (%.6e, %.6e)' % (delta, r[0], r[1]))
    print('            MULT = (%.6e, %.6e)、BEND = (%.6e, %.6e)' % (m[0], m[1], b[0], b[1]))
    print('            MULT - ODE = (%.2e, %.2e)、BEND - ODE = (%.2e, %.2e)、MULT - BEND = (%.2e, %.2e)' % (
        m[0] - r[0], m[1] - r[1], b[0] - r[0], b[1] - r[1], m[0] - b[0], m[1] - b[1]))
disp = lambda i: (res[0.01][i][0] - res[-0.01][i][0]) / 2 / 0.01     # x の δ の奇数部分 / δ = 分散 R16
c2_ = ((phi * ell)**2 * 2) / Lp / 48.0
pred = -c2_ / phi * (1 - math.cos(phi))                              # 円弧が c2 だけ短い (角は同じ): ρ が c2/φ だけ小さい → R16 = ρ (1 - cos φ) が変わる量
print('  分散 (δ の奇数部分 / δ) の誤差: MULT - ODE = %.3e、BEND - ODE = %.3e、MULT - BEND = %.3e' % (
    disp(1) - disp(0), disp(2) - disp(0), disp(1) - disp(2)))
print('  2 回目の補正 c2 = %.4e を円弧の長さから引くことによる R16 の変化の予想 = -c2 (1 - cos φ)/φ = %.3e' % (c2_, pred))
compat.set_mode('v2')
mv = {d: (track('MULT', [0, 0, 0, 0, 0, d]), track('BEND', [0, 0, 0, 0, 0, d])) for d in (0.01, -0.01)}
dv = lambda i: (mv[0.01][i][0] - mv[-0.01][i][0]) / 2 / 0.01
print('  [v2] 分散の誤差: MULT - ODE = %.3e、BEND - ODE = %.3e、MULT - BEND = %.3e (v2 は 2 回目の補正を省く。B331)' % (
    dv(0) - disp(0), dv(1) - disp(0), dv(0) - dv(1)))
compat.set_mode(None)
