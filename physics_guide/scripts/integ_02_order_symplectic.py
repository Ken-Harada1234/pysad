# -*- coding: utf-8 -*-
"""分割の一般論 (03_symplectic_integrators.md §6〜§8) の数値確認: 誤差の次数、シンプレクティック性、長時間の保存量。

H = p^2/2 + V(x)、V = x^2/2 + x^4/4 (1 自由度。ドリフト D = 運動量だけの流れと、キック K = 座標だけの流れの分割)。
(a) 誤差の次数: Strang (2 次)、Yoshida (4 次)、力勾配型 (4 次。tthin の構成)、RK4 (4 次。シンプレクティックでない)。刻みの数 n を 2 倍にしたときの誤差の比。
(b) シンプレクティック性: 1 刻みのヤコビアン (中心差分) の det と |M^T J M - J|。
(c) 長時間の H の誤差: シンプレクティックな方法は有界 (振動)、RK4 は系統的にずれる。
(d) 計算量 (力の評価の回数) あたりの誤差。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)

ALPHA = 1.0 / math.sqrt(12.0)
BETA = 1.0 / 6.0 - ALPHA / 2.0
GAMMA = 1.0 / 40.0 - 1.0 / (24.0 * math.sqrt(3.0))
C3 = BETA / 2.0          # 03 章の (c3、c5)。integ_01 で導出
C5 = GAMMA / 2.0
W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
W0 = 1.0 - 2.0 * W1


def V(x):
    return x * x / 2 + x**4 / 4


def V1(x):
    return x + x**3


def V2(x):
    return 1 + 3 * x * x


def V3(x):
    return 6 * x


def Hn(x, p):
    return p * p / 2 + V(x)


def drift(x, p, h):
    return x + h * p, p


def kick(x, p, h):
    return x, p - h * V1(x)


def strang(x, p, h):
    x, p = drift(x, p, h / 2)
    x, p = kick(x, p, h)
    return drift(x, p, h / 2)


def yoshida(x, p, h):
    for w in (W1, W0, W1):
        x, p = strang(x, p, w * h)
    return x, p


def force_gradient(x, p, h):
    """D(a h) K(h/2) D(α h) C D(α h) K(h/2) D(a h)、a = 1/2 - α。C: Δp = c3 h^3 ∂x(V'^2) - 2 c5 h^5 ∂x(V'^2 V'')"""
    a = 0.5 - ALPHA
    x, p = drift(x, p, a * h)
    x, p = kick(x, p, h / 2)
    x, p = drift(x, p, ALPHA * h)
    dW = 2 * V1(x) * V2(x)                                  # ∂x (V'^2)
    dG = 2 * V1(x) * V2(x)**2 + V1(x)**2 * V3(x)            # ∂x (V'^2 V'')
    p = p + C3 * h**3 * dW - 2 * C5 * h**5 * dG
    x, p = drift(x, p, ALPHA * h)
    x, p = kick(x, p, h / 2)
    return drift(x, p, a * h)


def force_gradient_c3only(x, p, h):
    """c5 を落としたもの (γ の項なし)"""
    a = 0.5 - ALPHA
    x, p = drift(x, p, a * h)
    x, p = kick(x, p, h / 2)
    x, p = drift(x, p, ALPHA * h)
    p = p + C3 * h**3 * 2 * V1(x) * V2(x)
    x, p = drift(x, p, ALPHA * h)
    x, p = kick(x, p, h / 2)
    return drift(x, p, a * h)


def rk4(x, p, h):
    f = lambda xx, pp: (pp, -V1(xx))
    k1 = f(x, p)
    k2 = f(x + h / 2 * k1[0], p + h / 2 * k1[1])
    k3 = f(x + h / 2 * k2[0], p + h / 2 * k2[1])
    k4 = f(x + h * k3[0], p + h * k3[1])
    return x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), p + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])


METHODS = (('Strang (2 次)', strang, 1, 0), ('Yoshida (4 次)', yoshida, 3, 0), ('力勾配型 (4 次、c3 と c5)', force_gradient, 2, 1),
           ('力勾配型 (c3 だけ)', force_gradient_c3only, 2, 1), ('RK4 (4 次、非シンプレクティック)', rk4, 4, 0))


def run(m, x0, p0, T, n):
    h = T / n
    x, p = x0, p0
    for _ in range(n):
        x, p = m(x, p, h)
    return x, p


x0, p0, T = 0.5, 0.3, 2.0
sol = solve_ivp(lambda s, v: [v[1], -V1(v[0])], (0, T), [x0, p0], method='DOP853', rtol=1e-13, atol=1e-15)
ref = sol.y[:, -1]

print('== (a) 誤差の次数 (H = p^2/2 + x^2/2 + x^4/4、x = 0.5、p = 0.3 から s = 2) ==')
print('   参照 (DOP853): x = %.12f、p = %.12f' % tuple(ref))
for name, m, nf, ng in METHODS:
    print('  %s (1 刻みの力の評価 %d 回、力の勾配の評価 %d 回):' % (name, nf, ng))
    prev = None
    line = []
    for n in (2, 4, 8, 16, 32, 64):
        xx, pp = run(m, x0, p0, T, n)
        err = max(abs(xx - ref[0]), abs(pp - ref[1]))
        line.append('n=%d: %.2e%s' % (n, err, '' if prev is None else ' (比 %.1f)' % (prev / err)))
        prev = err
    print('    ' + '、'.join(line))

print('\n== (b) シンプレクティック性 (1 刻み h = 0.3 のヤコビアン。中心差分 1e-6) ==')
for name, m, nf, ng in METHODS:
    xx, pp = 0.7, -0.2
    eps = 1e-6
    Jm = np.zeros((2, 2))
    for j in range(2):
        a = [xx, pp]
        b = [xx, pp]
        a[j] += eps
        b[j] -= eps
        ra, rb = m(a[0], a[1], 0.3), m(b[0], b[1], 0.3)
        Jm[:, j] = [(ra[0] - rb[0]) / (2 * eps), (ra[1] - rb[1]) / (2 * eps)]
    print('  %s: det J - 1 = %+.2e' % (name, np.linalg.det(Jm) - 1))
print('  (2 自由度以上では det = 1 だけでは足りず、M^T J M = J。D・K・力勾配のキックはどれも H の流れ (座標だけ、運動量だけの関数の流れ) なので厳密にシンプレクティック。')
print('   RK4 の det J - 1 は 1 刻みの誤差のオーダー ~h^5 で 0 でない)')

print('\n== (c) 長時間の H の誤差 (h = 0.25、5 万刻み = s = 12500。H0 = %.6f) ==' % Hn(0.5, 0.3))
H0 = Hn(0.5, 0.3)
for name, m, nf, ng in METHODS:
    x, p = 0.5, 0.3
    h = 0.25
    dmax_first = 0.0
    dmax_all = 0.0
    n = 50000
    for i in range(1, n + 1):
        x, p = m(x, p, h)
        d = abs(Hn(x, p) - H0)
        if i <= n // 10:
            dmax_first = max(dmax_first, d)
        dmax_all = max(dmax_all, d)
    print('  %s: 最初の 5000 刻みの max|ΔH| = %.2e、全体の max|ΔH| = %.2e (比 %.1f。有界なら 1 に近い)、最後の ΔH = %+.2e' % (name, dmax_first, dmax_all, dmax_all / dmax_first, Hn(x, p) - H0))

print('\n== (d) 計算量あたりの誤差 (s = 2 を、力の評価の総数が約 64 回になる刻みの数で) ==')
for name, m, nf, ng in METHODS:
    n = max(1, int(round(64 / (nf + ng))))
    xx, pp = run(m, x0, p0, T, n)
    err = max(abs(xx - ref[0]), abs(pp - ref[1]))
    print('  %s: n = %d (評価 %d 回) 誤差 %.2e' % (name, n, n * (nf + ng), err))
