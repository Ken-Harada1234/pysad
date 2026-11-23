# -*- coding: utf-8 -*-
"""長さのある多極 (L ≠ 0) の 4 次の分割 tthin (07_thin_multipole.md §6〜§8) を、厳密な運動方程式の数値積分と比べる。

H = -sqrt(p^2 - px^2 - py^2) + κ Re (x - i y)^(k+1)/(k+1)!  (δ = 0、p = 1。z' = 1 - 1/pz)。κ = K/L は長さあたりの強さ。
(a) pysad の追跡 (tthin。SEXT の DISFRIN=1 でフリンジなし) と、K64 の式を書き直した関数 (この script の tthin_k64) の一致。
(b) 誤差の次数: L を 2 倍にしたときの誤差 (厳密な積分との差) の比。tthin (4 次: 2^5 = 32)、tthine 型 (2 次: 2^3 = 8)。γ の項の有無。
(c) ヤコビアン: 追跡の数値微分、CALC の tthine の行列 (elkern._jac_thin)、厳密な積分の数値微分。
(d) 6 極 (k = 2) と 8 極 (k = 3)。
"""
import io
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
from pysad.ffscmd import run_sad_string
from pysad import elkern

ALPHA = 1.0 / math.sqrt(12.0)
ALPHA1 = 0.5 - ALPHA
BETA = 1.0 / 6.0 - ALPHA / 2.0
GAMMA = 1.0 / 40.0 - 1.0 / (24.0 * math.sqrt(3.0))


def V_grad(x, y, kappa, k):
    """∇V、V = κ/(k+1)! Re (x - i y)^(k+1)。V_x + i V_y = κ/k! (x - i y)^k"""
    g = kappa / math.factorial(k) * complex(x, -y)**k
    return g.real, g.imag


def exact(c, L, kappa, k):
    """運動学的な厳密な流れ (δ = 0)"""
    def rhs(s, v):
        x, px, y, py, z = v
        pz = math.sqrt(1.0 - px * px - py * py)
        gx, gy = V_grad(x, y, kappa, k)
        return [px / pz, -gx, py / pz, -gy, 1.0 - 1.0 / pz]
    sol = solve_ivp(rhs, (0, L), list(c[:5]), method='DOP853', rtol=1e-13, atol=1e-17)
    return np.array(list(sol.y[:, -1]) + [0.0])


def drift(c, h):
    c = c.copy()
    pz = math.sqrt(1.0 - c[1]**2 - c[3]**2)
    a = c[1]**2 + c[3]**2
    c[0] += c[1] / pz * h
    c[2] += c[3] / pz * h
    c[4] -= a / (1.0 + pz) / pz * h
    return c


def kick(c, ak_total, k, L=None):
    """全体の積分強度 ak_total の分の薄いキック (px -= ak/k! Re (x - i y)^k)"""
    c = c.copy()
    g = ak_total / math.factorial(k) * complex(c[0], -c[2])**k
    c[1] -= g.real
    c[3] -= g.imag
    return c


def tthin_k64(c, L, ak, k, gamma_scale=1.0, beta_scale=1.0):
    """K64 の tthin (tquad.f:118-、フリンジなし、δ = 0、p = 1) を書き直したもの。一般の kord。γ・β の倍率でその項を切れる"""
    beta, gamma = BETA * beta_scale, GAMMA * gamma_scale
    ala, alb = L * ALPHA1, L * ALPHA
    akf = ak / math.factorial(k) * 0.5
    c = drift(c, ala)
    c = kick(c, ak * 0.5, k)
    # 中央: ドリフト alb → 補正キック → ドリフト alb
    c = drift(c, alb)
    x, y = c[0], c[2]
    r = x * x + y * y
    cx = complex(x, -y)**k
    rcx1 = x * cx.real + y * cx.imag
    aki = akf * 2.0
    if k == 2:
        f1, f2, f3, f4, f5 = beta * 2 * L, 4 * gamma * L**2, 6 * gamma * L**2, 0.5 * beta * L, 4 * gamma * L**2
        rk1, rk = 1.0, r
    else:
        f1 = beta * k * L
        f2 = 2 * gamma * k * (k - 1) * L**2
        f3 = gamma * k * (k + 1) * L**2
        f4 = 0.5 * beta * L
        f5 = 2 * gamma * k * L**2
        rk1 = r**(k - 2)
        rk = r * rk1
    c = c.copy()
    c[1] += aki**2 * (f1 * rk * x - aki * (f2 * rk1 * x * rcx1 + f3 * rk * cx.real))
    c[3] += aki**2 * (f1 * rk * y - aki * (f2 * rk1 * y * rcx1 + f3 * rk * cx.imag))
    c[4] += aki**2 * rk * (f4 * r - f5 * aki * rcx1)
    c = drift(c, alb)
    c = kick(c, ak * 0.5, k)
    return drift(c, ala)


def tthine_type(c, L, ak, k):
    """CALC の tthine と同じ 2 次の分割 (L/6、キック 1/2、2L/3、キック 1/2、L/6)"""
    c = drift(c, L / 6.0)
    c = kick(c, ak * 0.5, k)
    c = drift(c, L / 1.5)
    c = kick(c, ak * 0.5, k)
    return drift(c, L / 6.0)


def strang(c, L, ak, k):
    """比較用: ドリフト L/2、キック (全体)、ドリフト L/2 (2 次)"""
    c = drift(c, L / 2.0)
    c = kick(c, ak, k)
    return drift(c, L / 2.0)


def thin_ffs(name, key, ak, L):
    src = ('MOMENTUM = 3GEV;\n%s MP=(L=%.15g %s=%.15g DISFRIN=1);\nMARK ST = ();\nLINE ONE = (ST MP);\nFFS;\nUSE ONE;\n' % (name, L, key, ak))
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


NAMES = {2: ('SEXT', 'K2'), 3: ('OCT', 'K3')}
c0 = np.array([0.020, 0.0010, -0.010, -0.0005, 0.0, 0.0])
KAPPA = {2: 400.0, 3: 4.0e4}            # 長さあたりの強さ [1/m^(k+1)]。6 極 400 /m^3、8 極 4e4 /m^4

print('== (a) pysad の追跡 (tthin) と、K64 の式を書き直した関数 (tthin_k64) の比較 ==')
for k in (2, 3):
    name, key = NAMES[k]
    worst = 0.0
    for L in (0.1, 0.4):
        ak = KAPPA[k] * L
        ffs = thin_ffs(name, key, ak, L)
        tr = C.track_one(ffs, c0)
        mine = tthin_k64(c0, L, ak, k)
        worst = max(worst, np.abs(tr - mine).max())
    print('  %s (k = %d、κ = %g): 最大差 %.2e (x, px, y, py, z)' % (name, k, KAPPA[k], worst))

def run_b(c_init, Ls, label):
    for k in (2, 3):
        name, key = NAMES[k]
        print('  %s (k = %d、κ = %g /m^%d)、%s:' % (name, k, KAPPA[k], k + 1, label))
        rows = {}
        for L in Ls:
            ak = KAPPA[k] * L
            ref = exact(c_init, L, KAPPA[k], k)
            ffs = thin_ffs(name, key, ak, L)
            e = {
                'pysad の追跡 (tthin、4 次)': np.abs(C.track_one(ffs, c_init) - ref).max(),
                'tthin から γ の項を除く': np.abs(tthin_k64(c_init, L, ak, k, gamma_scale=0.0) - ref).max(),
                'tthin から β, γ の項を除く (補正なし)': np.abs(tthin_k64(c_init, L, ak, k, gamma_scale=0.0, beta_scale=0.0) - ref).max(),
                'tthine 型 (CALC、2 次)': np.abs(tthine_type(c_init, L, ak, k) - ref).max(),
                'Strang (L/2、キック、L/2)': np.abs(strang(c_init, L, ak, k) - ref).max(),
            }
            for key_, v in e.items():
                rows.setdefault(key_, []).append(v)
        for key_, vals in rows.items():
            ratios = ['%.1f' % (vals[i + 1] / vals[i]) for i in range(len(vals) - 1)]
            print('    %-34s %s  (比 %s)' % (key_, ', '.join('%.2e' % v for v in vals), ', '.join(ratios)))


print('\n== (b) 厳密な積分 (DOP853) との差と、L を 2 倍にしたときの比 (κ 一定) ==')
Ls = (0.0125, 0.025, 0.05, 0.1, 0.2)
print('  L = %s [m]' % ', '.join(str(v) for v in Ls))
run_b(c0, Ls, '初期値 x = 20 mm、y = -10 mm、px = 1 mrad、py = -0.5 mrad (傾きが小さい)')
c_steep = np.array([0.020, 0.030, -0.010, -0.020, 0.0, 0.0])
run_b(c_steep, Ls, '初期値 x = 20 mm、y = -10 mm、px = 30 mrad、py = -20 mrad (傾きが大きい)')
print('  (4 次の方法は 1 刻みの誤差が L^5 なので比は 32、2 次は L^3 で 8。「補正なし」は β, γ の補正キックを落とした分割で、ドリフト・キックの長さの割り振り (α) だけが')
print('   [A,[A,B]] の項 (運動量の 2 乗に比例) を消し、残る [B,[B,A]] の項 (力の勾配 |∇V|^2、傾きによらない) で 2 次。tthine 型は逆に [B,[B,A]] の項を消し、')
print('   [A,[A,B]] の項 (傾きの 2 乗に比例) が残る: 傾きが小さいと 5 次の項が見え (比 32)、傾きが大きいと 3 次の項が見える (比 8))\n')

print('\n== (c) ヤコビアン (L = 0.1、κ = %g、COD は c0): 追跡の数値微分と、CALC の tthine の行列と、厳密な積分の数値微分 ==' % KAPPA[2])
k = 2
L = 0.1
ak = KAPPA[k] * L
name, key = NAMES[k]
ffs = thin_ffs(name, key, ak, L)


def jac(f, c, h=1e-6):
    Jm = np.zeros((6, 6))
    for j in range(6):
        cp, cm = c.copy(), c.copy()
        cp[j] += h
        cm[j] -= h
        Jm[:, j] = (f(cp) - f(cm)) / (2 * h)
    return Jm


BLK = (slice(0, 5), slice(0, 4))        # x, px, y, py, z の行、x, px, y, py の列 (厳密な積分は δ を含まないので δ の列は比べない)
J_ex = jac(lambda c: exact(c, L, KAPPA[k], k), c0)
J_tr = jac(lambda c: C.track_one(ffs, c), c0)
ok, J_calc = elkern._jac_thin(c0.copy(), L, ak, k, 0.0)
print('  厳密な積分のヤコビアン (x, px, y, py, z の行・x, px, y, py の列) との最大差: pysad の追跡 (tthin) %.2e、CALC の tthine の行列 %.2e' % (np.abs((J_tr - J_ex)[BLK]).max(), np.abs((J_calc - J_ex)[BLK]).max()))
print('  追跡のヤコビアンと CALC の行列の最大差 (同じ範囲) %.2e' % np.abs((J_tr - J_calc)[BLK]).max())
Jm = np.zeros((6, 6))
for i in range(3):
    Jm[2 * i, 2 * i + 1] = 1.0
    Jm[2 * i + 1, 2 * i] = -1.0
print('  シンプレクティック性 |M^T J M - J|: 追跡 %.1e (中心差分の精度)、CALC の行列 %.1e' % (np.abs(J_tr.T @ Jm @ J_tr - Jm).max(), np.abs(J_calc.T @ Jm @ J_calc - Jm).max()))
print('  設計軌道 (COD = 0) では、追跡と CALC の行列の差 (同じ範囲): %.1e (6 極のキックが 0 で、どちらもドリフト)' % np.abs((jac(lambda c: C.track_one(ffs, c), np.zeros(6)) - elkern._jac_thin(np.zeros(6), L, ak, k, 0.0)[1])[BLK]).max())
