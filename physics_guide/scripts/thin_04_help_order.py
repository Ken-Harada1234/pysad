# -*- coding: utf-8 -*-
"""ヘルプ (SADHelp.html) の SEXT・OCT・DECA・DODECA の transformation の式 (並びと補正項 V) の確認 (07 章・03 章 付録 A)。

ヘルプの式 (TeX・平文とも同じ):
  exp(:Fin:) exp(:aL:) exp(:H/2:) exp(:bL:) exp(:V:) exp(:aL:) exp(:H/2:) exp(:bL:) exp(:Fout:)、a = 1/2 - 1/√12、b = 1/2 - a、
  V = Σ_{j,k} [ -β/2 (H,k)^2 + γ H,j H,k H,jk ]、β = 1/6 - 1/√48、γ = 1/40 - 1/(24√3)。
exp(:H:) は「H を時間 1 のハミルトニアンとする流れ」 (drift_04_help_convention.py)。
1 自由度 (x, px)、H_ham = px^2/2 + v(x)、v = κ x^3/6 (6 極)、p = 1。厳密な流れ (DOP853) との差の L 依存で次数を見る。
 (A) 前後対称の並び a, H/2, b, V, b, H/2, a (K64 の追跡、07 章・03 章)
 (B) ヘルプの並びを、そのまま読んだもの a, H/2, b, V, a, H/2, b (非対称)
補正 V は、キックの強さを L v(x) (H は積分強度) として V = -(β/2) L^3 (v')^2 + γ L^5 (v')^2 v'' (= -(β/2) L (H')^2 + γ L^2 H'^2 H'' / ... の積分強度 H = L v での書き換え。L・1/p の因子つき)。
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
t = 1 / math.sqrt(12.0)
a = 0.5 - t
beta = 1 / 6 - 1 / math.sqrt(48.0)
gamma = 1 / 40 - 1 / (24 * math.sqrt(3.0))
kap = 40.0
v1 = lambda x: kap * x**2 / 2      # v'
v2 = lambda x: kap * x              # v''


def exact(x, px, L):
    s = solve_ivp(lambda s, w: [w[1], -v1(w[0])], (0, L), [x, px], method='DOP853', rtol=1e-13, atol=1e-16)
    return s.y[:, -1]


def drift(w, l):
    return [w[0] + w[1] * l, w[1]]


def kick(w, f, L):
    return [w[0], w[1] - f * L * v1(w[0])]


def corr(w, L, sg=1.0, use_g=True):
    # W = -(β/2) L^3 v'^2 + γ L^5 v'^2 v''  (流れの Hamiltonian の項)、Δpx = -∂W/∂x
    x = w[0]
    dW = -(beta / 2) * L**3 * 2 * v1(x) * v2(x) + (gamma * L**5 * (2 * v1(x) * v2(x) * v2(x) + v1(x)**2 * kap) if use_g else 0.0)
    return [x, w[1] - sg * dW]


def seq(order, w, L, sg=1.0, use_g=True):
    for kind, arg in order:
        if kind == 'D':
            w = drift(w, arg * L)
        elif kind == 'K':
            w = kick(w, arg, L)
        else:
            w = corr(w, L, sg, use_g)
    return np.array(w)


A_ = [('D', a), ('K', .5), ('D', t), ('C', 0), ('D', t), ('K', .5), ('D', a)]
B_ = [('D', a), ('K', .5), ('D', t), ('C', 0), ('D', a), ('K', .5), ('D', t)]
x0, p0 = 0.3, 0.1
print('x0 = %.2f, px0 = %.2f, 6 極 κ = %.0f (v = κ x^3/6)。誤差 = |分割 - 厳密| (x, px の最大)' % (x0, p0, kap))
Ls = [0.4, 0.2, 0.1, 0.05]
for name, order, sg, ug in (('(A) 対称 a,H/2,b,V,b,H/2,a   (V の符号 正しい)', A_, 1.0, True),
                            ('(A) 対称、V の符号を逆', A_, -1.0, True),
                            ('(A) 対称、γ の項なし', A_, 1.0, False),
                            ('(B) ヘルプの並び a,H/2,b,V,a,H/2,b', B_, 1.0, True)):
    errs = []
    for L in Ls:
        e = np.max(np.abs(seq(order, [x0, p0], L, sg, ug) - exact(x0, p0, L)))
        errs.append(e)
    sl = [math.log(errs[i] / errs[i + 1]) / math.log(2) for i in range(len(Ls) - 1)]
    print(name)
    print('   L = ' + ', '.join('%.2f' % L for L in Ls) + ' の誤差: ' + ', '.join('%.2e' % e for e in errs)
          + '   次数 (対数の傾き) = ' + ', '.join('%.2f' % s for s in sl))
print('\n[β, γ] ヘルプの β = 1/6 - 1/√48 = %.12f (07・03 章 c3 = β/2 = %.12f)、γ = %.12e (c5 = γ/2)' % (beta, beta / 2, gamma))
print('3 次の条件 c3 = 1/12 - t/4 = %.12f' % (1 / 12 - t / 4))
