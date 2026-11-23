# -*- coding: utf-8 -*-
r"""BEND の章 (05_bend.md) §6.4 (5.43) の L' の補正を、一般の ψ1 ≠ ψ2・斜め入射・出口まで広げて確かめる (numpy・scipy・mpmath・sympy。pysad は使わない)。

(a) 弦の幾何: 入口の面の点 -a1 e_x1、出口の面の点 E - a2 e_f2 を通る円の弦の変化 (射影) と横ずれ
    (a1 = a2 = s0 で ΔC = -2 s0 sin(ω/2) cos((ψ1-ψ2)/2)、横ずれ = 2 s0 sin(ω/2) sin((ψ1-ψ2)/2)、
    ヘルプの式は |ΔC| = 2 s0 sin(ω/2)。差 (ψ1-ψ2)^2/8 の数値)
(b) 面の座標系での端の通過: 面の法線から ψ の角で入る粒子の x のずれ。(5.41) の斜め入射への拡張 ΔX = I1 p^2/(rho p_z^3)
    (1 次の摂動を求積し、1/rho の 1 次で外挿する)。出口の ΔX の符号 (入口と反対で、円が通る点は E - s0 e_f2)
(c) 実験室系で、入口の面と出口の面の端の場 b_y = g(s1) g(-s2)/rho_b を水平面の ODE (DOP853) で通し、
    出口の基準面で p_xr = 0 (出口の向きが設計どおり) になる rho_b から L' を求めて、射影の式・ヘルプの式 (5.43) と比べる。
    ヘルプの L' を使ったときの出口の x_out・p_xr の残りも出す
"""
import math
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import solve_ivp
import _common as C

C.start(__file__)
mp.mp.dps = 30

print('== (a) 弦の幾何 (sympy・数値) ==')
phi, p1, p2, s0 = sp.symbols('phi psi1 psi2 s0', real=True)
ex1 = sp.Matrix([sp.cos(p1), sp.sin(p1)])
ef2 = sp.Matrix([sp.cos(phi - p2), sp.sin(phi - p2)])
Eh = sp.Matrix([-sp.sin(phi / 2), sp.cos(phi / 2)])
Ep = sp.Matrix([sp.cos(phi / 2), sp.sin(phi / 2)])
D = s0 * (ex1 - ef2)          # 円が通る 2 点の差の変化 (E は除く): (E - s0 ef2) - (-s0 ex1) = E + s0 (ex1 - ef2)
om = phi - p1 - p2
dC = sp.simplify(sp.expand_trig(D.dot(Eh)))
dP = sp.simplify(sp.expand_trig(D.dot(Ep)))
print('  ΔC (E_hat 方向) - (-2 s0 sin(ω/2) cos((ψ1-ψ2)/2)) =', sp.simplify(sp.expand_trig(dC + 2 * s0 * sp.sin(om / 2) * sp.cos((p1 - p2) / 2))))
print('  横ずれ (E_perp 方向) - (2 s0 sin(ω/2) sin((ψ1-ψ2)/2)) =', sp.simplify(sp.expand_trig(dP - 2 * s0 * sp.sin(om / 2) * sp.sin((p1 - p2) / 2))))
print('  |D|^2 - (2 s0 sin(ω/2))^2 =', sp.simplify(sp.expand_trig(D.dot(D) - (2 * s0 * sp.sin(om / 2))**2)))
print('  d1 = (0,1) の接線と E_hat の角、e_x1·E_hat = -sin(φ/2-ψ1)、e_f2·E_hat = sin(φ/2-ψ2) の確認:',
      sp.simplify(sp.expand_trig(ex1.dot(Eh) + sp.sin(phi / 2 - p1))), sp.simplify(sp.expand_trig(ef2.dot(Eh) - sp.sin(phi / 2 - p2))))
print('  ヘルプの式 (|D| = 2 s0 sin(ω/2)) に対する射影の比 cos((ψ1-ψ2)/2) と 1 - (ψ1-ψ2)^2/8:')
for dpsi in (0.04, 0.1, 0.2, 0.4):
    print('    ψ1-ψ2 = %.2f: cos(Δψ/2) = %.8f、1 - Δψ^2/8 = %.8f、差 %.2e' % (dpsi, math.cos(dpsi / 2), 1 - dpsi**2 / 8, math.cos(dpsi / 2) - (1 - dpsi**2 / 8)))

print('\n== (b) 面の座標系での端の通過: 斜め入射の ΔX (tanh の端) ==')
mp.mp.dps = 30
d = mp.mpf('0.01')        # tanh の幅 (F1 = 3 d)


def Gsoft(s):
    """G = ∫_{-∞}^s g ds'、g = (1 + tanh(s/d))/2"""
    return s / 2 + d / 2 * mp.log(mp.cosh(s / d)) + d / 2 * mp.log(2)


def Ghard(s):
    return s if s > 0 else mp.mpf(0)


def F(px, pp):
    return px / mp.sqrt(pp**2 - px**2)


I1 = -mp.pi**2 * d**2 / 24          # I1 = ∫ s (g - Θ) ds (tanh)


def shift(sinpsi_p, pp, eps, mode):
    """入口 (mode = 'in'): px(s) = px0 - eps G(s)。出口 (mode = 'out'): g_out(s) = g(-s)。
    戻り値: 軟らかい端と硬い端の x の差 (ramp の外で ∫ (F(px) - F(px_hard)) ds)"""
    px0 = sinpsi_p
    lim = 14 * d
    if mode == 'in':
        f = lambda s: F(px0 - eps * Gsoft(s), pp) - F(px0 - eps * Ghard(s), pp)
        return mp.quad(f, [-lim, 0, lim])
    # 出口: 入る前の px0 (本体の中) から。G_out(s) = ∫_{-lim}^s g(-s') ds'、硬い端 ∫_{-lim}^s Θ(-s') ds'
    Gt = lambda s: Gsoft(lim) - Gsoft(-s)           # ∫_{-lim}^s g(-s') ds' = ∫_{-s}^{lim} g du
    Gh = lambda s: (-lim - 0) * 0 + (min(s, 0) + lim)
    f = lambda s: F(px0 - eps * Gt(s), pp) - F(px0 - eps * Gh(s), pp)
    return mp.quad(f, [-lim, 0, lim])


print('  ΔX = (ε d F/dp_x ∫D) の 1 次: ΔX/(ε I1 p^2/p_z^3) (入口)、-ΔX/(ε I1 p^2/p_z^3) (出口)。ε = 1e-6、d = %s、p = 1.0 と 1.03' % mp.nstr(d, 3))
for pp in (mp.mpf(1), mp.mpf('1.03')):
    for psi in (0, 0.1, 0.3, 0.6):
        px0 = pp * mp.sin(psi)
        pz = mp.sqrt(pp**2 - px0**2)
        eps = mp.mpf('1e-6')
        pred = I1 * pp**2 / pz**3
        xin = shift(px0, pp, eps, 'in') / eps
        xout = shift(px0, pp, eps, 'out') / eps
        print('    p = %s ψ = %.1f: 入口 比 = %s、出口 -比 = %s    (1/cos^3 ψ = %.6f、ΔX/ε = %s)' %
              (mp.nstr(pp, 3), psi, mp.nstr(xin / pred, 8), mp.nstr(-xout / pred, 8), 1 / math.cos(psi)**3, mp.nstr(xin, 6)))

print('\n== (c) 実験室系の ODE: 入口・出口の端の場で L\' を数値的に求める ==')
RHO0 = 10.0
PHI = 0.1
L0 = RHO0 * PHI


def gfun(s, dd):
    return 0.5 * (1.0 + math.tanh(s / dd))


def run(rhob, psi1, psi2, dd):
    """p = 1 の電子。実験室系 (X, Z)。入口の基準面 Z = 0 (面は原点を通り法線 (-sin ψ1, cos ψ1))、出口は E を通り法線 (-sin(φ-ψ2), cos(φ-ψ2))。
    戻り値: 出口の基準面での (x_out, p_xr)"""
    E = np.array([-RHO0 * (1 - math.cos(PHI)), RHO0 * math.sin(PHI)])
    n1 = np.array([-math.sin(psi1), math.cos(psi1)])
    n2 = np.array([-math.sin(PHI - psi2), math.cos(PHI - psi2)])
    e2 = np.array([math.cos(PHI - psi2), math.sin(PHI - psi2)])
    nref = np.array([-math.sin(PHI), math.cos(PHI)])
    eref = np.array([math.cos(PHI), math.sin(PHI)])
    z0 = -22.0 * dd / math.cos(psi1)
    y0 = np.array([0.0, z0, 0.0, 1.0])

    def rhs(t, y):
        X, Z, PX, PZ = y
        s1 = X * n1[0] + Z * n1[1]
        s2 = (X - E[0]) * n2[0] + (Z - E[1]) * n2[1]
        b = gfun(s1, dd) * gfun(-s2, dd) / rhob
        return [PX, PZ, -PZ * b, PX * b]

    def ev(t, y):
        return ((y[0] - E[0]) * n2[0] + (y[1] - E[1]) * n2[1]) - 22.0 * dd
    ev.terminal = True
    sol = solve_ivp(rhs, [0, 5.0], y0, method='DOP853', rtol=1e-13, atol=1e-15, events=ev)
    X, Z, PX, PZ = sol.y[:, -1]
    r = np.array([X, Z])
    P = np.array([PX, PZ])
    t = ((E - r) @ nref) / (P @ nref)
    rr = r + t * P
    return (rr - E) @ eref, P @ eref


def root_pxr(psi1, psi2, dd):
    f = lambda rb: run(rb, psi1, psi2, dd)[1]
    a, b = RHO0 * (1 - 1e-4), RHO0 * (1 - 2e-4)
    fa, fb = f(a), f(b)
    for _ in range(60):
        c = b - fb * (b - a) / (fb - fa)
        a, fa = b, fb
        b, fb = c, f(c)
        if abs(b - a) < 1e-13:
            break
    return b


def models(psi1, psi2, sI):
    """sI = -I1/rho (tanh では pi^2 d^2/(24 rho))。戻り値 (help, M1 (斜め 1/cos^3)、M1 (面に垂直な式)、M1 の残り x_out) すべて sI 単位の L'-L と x_out"""
    om_ = PHI - psi1 - psi2
    k = PHI / (2 * math.sin(PHI / 2))
    out = {}
    out['help'] = -PHI * math.sin(om_ / 2) / math.sin(PHI / 2)
    for name, (a1, a2) in (('M1', (1 / math.cos(psi1)**3, 1 / math.cos(psi2)**3)), ('M1n', (1.0, 1.0))):
        dCp = -a1 * math.sin(om_) / math.cos(PHI / 2 - psi2)
        out[name] = dCp * k
        ex1v = np.array([math.cos(psi1), math.sin(psi1)])
        ef2v = np.array([math.cos(PHI - psi2), math.sin(PHI - psi2)])
        Ehv = np.array([-math.sin(PHI / 2), math.cos(PHI / 2)])
        eref = np.array([math.cos(PHI), math.sin(PHI)])
        out[name + '_x'] = (-a1 * ex1v + dCp * Ehv + a2 * ef2v) @ eref
    return out


print("  M1 : 円弧の角 = φ (p_xr = 0) の条件 ΔC' = -a1 sin ω/cos(φ/2-ψ2)、a_i = s_I/cos^3ψ_i (斜め入射)。M1n: a_i = s_I (面に垂直な式)")
print("  すべて sI = -I1/rho = π^2 d^2/(24 rho) の単位 (tanh の端)。F1 = 0.06, 0.03, 0.015 の線形外挿")
cases = [(0.0, 0.0), (0.02, 0.02), (0.3, 0.3), (0.3, 0.0), (0.3, -0.1), (-0.05, 0.25)]
for (psi1, psi2) in cases:
    vals = []
    for F1 in (0.06, 0.03, 0.015):
        dd = F1 / 3.0
        sI = math.pi**2 * dd**2 / (24 * RHO0)
        rb = root_pxr(psi1, psi2, dd)
        xr_, pr_ = run(rb, psi1, psi2, dd)
        dl = (rb * PHI - L0) / sI
        F1e = math.pi * dd       # F1^2/24 を -I1 = π^2 d^2/24 に置き換えた実効の F1 (tanh の端の形の因子 π^2/9 を除く)
        rbh = (L0 - (PHI * F1e)**2 / (24 * L0) * math.sin((PHI - psi1 - psi2) / 2) / math.sin(PHI / 2)) / PHI
        xh, ph = run(rbh, psi1, psi2, dd)
        vals.append((dl, xr_ / sI, xh / sI, ph / sI))
    ex = [2 * vals[2][i] - vals[1][i] for i in range(4)]
    print('      (参考) F1 = 0.06, 0.03, 0.015 の (L-prime - L)/sI の生の値: ' + ', '.join('%.5f' % v[0] for v in vals) + '  (F1 に比例してずれる = F1^2 の次の項)')
    m = models(psi1, psi2, 1.0)
    print("  ψ1 = %5.2f、ψ2 = %5.2f (ω = %5.2f):" % (psi1, psi2, PHI - psi1 - psi2))
    print("      (L'-L)/sI: ODE %.5f | M1 %.5f | M1n %.5f | ヘルプ (5.43) %.5f (F1^2/24 を -I1 = π^2 d^2/24 に置き換えて、sI 単位)" % (ex[0], m['M1'], m['M1n'], m['help']))
    print("      p_xr = 0 の ρ_b での x_out/sI: ODE %.5f | M1 の予測 %.5f" % (ex[1], m['M1_x']))
    print("      ヘルプの L' を使った ρ_b での ODE: x_out/sI = %.5f、p_xr/sI = %.5f" % (ex[2], ex[3]))
