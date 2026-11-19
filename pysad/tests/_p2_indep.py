# -*- coding: utf-8 -*-
"""担当 P2 の独立な計算 (SAD・pysad の光学を使わない): ローレンツ力の数値積分 (DOP853、rtol 1e-13) による硬い端の BEND の厳密な写像。
test_v2_B1_optics._exact_map (K64 の模型: 四極は半径方向の面の間) を、四極も端面の間とする模型 (quad='face') を選べるように拡張したもの。"""
import math
import numpy as np
from scipy.integrate import solve_ivp

_L = 1.2
_TH = 0.2617993877991494


def exact_map(c, k, e1, e2, phi=_TH, L=_L, cq=-0.5, quad='radial'):
    """硬い端の厳密な写像 (ローレンツ力の数値積分、DOP853 rtol 1e-13)。磁石の形: 双極磁場 B_y = h (入口・出口の端面の間、端の角 e1, e2 [rad])、
    勾配の磁場 B_y = k x + cq k h x² (x = r − ρ0。中面。cq = −1/2 は K64 の tbendicorr の蹴り −K1 x² h/2 から)、B_r = k y (半径方向、y の 1 次)
    が、半径方向の入口・出口の面の間。b = 軸上の曲率中心側を +r 向きにとる (ANGLE>0 は SAD の座標と同じ x を外向きにとる)。
    戻り (x, px) = 半径方向の出口の面での外向きの座標と運動量。端面の外はドリフト。
    quad='radial': 四極は半径方向の入口・出口の面の間 (K64 の模型)。quad='face': 四極も端面の間 (端面から始まる硬い端の磁石)"""
    rho0 = L / phi
    h = 1 / rho0
    x, px, y, py, z, dl = c
    p = 1 + dl
    pz = math.sqrt(p * p - px * px - py * py)
    Pin = np.array([rho0, 0.0, 0.0])
    Pout = np.array([rho0 * math.cos(phi), 0.0, rho0 * math.sin(phi)])
    tin = np.array([0, 0, 1.0])
    rin = np.array([1.0, 0, 0])
    tout = np.array([-math.sin(phi), 0, math.cos(phi)])
    rout = np.array([math.cos(phi), 0, math.sin(phi)])
    n1 = math.cos(e1) * tin - math.sin(e1) * rin
    n2 = math.cos(e2) * tout + math.sin(e2) * rout
    s1 = lambda w: (w[0] - Pin[0]) * n1[0] + (w[2] - Pin[2]) * n1[2]
    s2 = lambda w: (w[0] - Pout[0]) * n2[0] + (w[2] - Pout[2]) * n2[2]
    r1 = lambda w: w[2]
    r2 = lambda w: -w[0] * math.sin(phi) + w[2] * math.cos(phi)
    P = np.array([rho0 + x, y, 0.0])
    pv = np.array([px, py, pz])
    u = pv / np.linalg.norm(pv)
    w = np.concatenate([P - u * 0.4, pv])
    planes = [s1, r1, s2, r2]
    inside_dip = lambda q: s1(q) > 0 and s2(q) < 0
    inside_quad = (lambda q: r1(q) > 0 and r2(q) < 0) if quad == 'radial' else inside_dip      # 'face': 四極の磁場も端面の間
    for _ in range(12):
        uu = w[3:] / np.linalg.norm(w[3:])
        q = w.copy()
        q[:3] = q[:3] + uu * 1e-9
        dip, quad = inside_dip(q), inside_quad(q)

        def f(l, ww):
            X, Y, Z, pX, pY, pZ = ww
            pp = math.sqrt(pX * pX + pY * pY + pZ * pZ)
            ux, uy, uz = pX / pp, pY / pp, pZ / pp
            rr = math.hypot(X, Z)
            xx = rr - rho0
            By = (h if dip else 0.0) + (k * xx + cq * k * h * xx * xx if quad else 0.0)
            Br = k * Y if quad else 0.0
            Bx, Bz = Br * X / rr, Br * Z / rr
            return [ux, uy, uz, uy * Bz - uz * By, uz * Bx - ux * Bz, ux * By - uy * Bx]
        evs = []
        for pl in planes:
            sg = np.sign(pl(q))

            def mk(pl=pl, sg=sg):
                def e(l, ww):
                    return pl(ww) * sg
                e.terminal = True
                e.direction = -1
                return e
            evs.append(mk())
        sol = solve_ivp(f, [0, 3.0], w, method='DOP853', rtol=1e-13, atol=1e-15, events=evs)
        hit = [j for j in range(4) if len(sol.t_events[j])]
        if not hit:
            w = sol.y[:, -1]
            break
        j = min(hit, key=lambda j: sol.t_events[j][0])
        w = sol.y_events[j][0]
        if s2(w) > 0 and r2(w) > 0 and j in (2, 3):
            break
    u = w[3:] / np.linalg.norm(w[3:])
    lstr = -np.dot(w[:3], tout) / np.dot(u, tout)
    Pf = w[:3] + u * lstr
    rr = math.hypot(Pf[0], Pf[2])
    return np.array([rr - rho0, np.dot(w[3:], rout)])




# ---------------------------------------------------------------------------------------------
# 放射反作用つきのローレンツ力の数値積分 (一様磁場の扇形の偏向磁石)。tradke・tradkn の ddpx の符号を決める独立な基準 (v2next P2 項目 3)。
# ---------------------------------------------------------------------------------------------
def bend_with_radiation(x0, px0, g0, theta, L, k):
    """水平面の一様磁場の扇形の偏向磁石 (曲率 h0 = theta/L、基準運動量 1) を、連続な放射の反作用つきで通る。無次元:
      dϑ/ds = -h0/pr,  dx/ds = sinϑ,  dz/ds = cosϑ,  dpr/ds = -k h0² pr²     (損失は運動量の方向、∝ γ² B²。高エネルギーの極限)
    曲線座標 (設計軌道の出口の面) の (x, px, δ) と経路長を返す。k = 0 は放射なし (放射のない写像と 1e-12 で一致する)。"""
    h0 = theta / L
    rho = 1.0 / h0
    pr0 = 1.0 + g0
    th0 = np.arcsin(px0 / pr0)
    Pe = np.array([-rho * (1 - np.cos(theta)), rho * np.sin(theta)])
    te = np.array([-np.sin(theta), np.cos(theta)])
    ne = np.array([np.cos(theta), np.sin(theta)])

    def rhs(s, y):
        x, z, th, pr = y
        return [np.sin(th), np.cos(th), -h0 / pr, -k * h0 ** 2 * pr ** 2]

    def ev(s, y):
        return (np.array([y[0], y[1]]) - Pe) @ te
    ev.terminal = True
    ev.direction = 1
    sol = solve_ivp(rhs, [0, 3 * L], [x0, 0.0, th0, pr0], events=ev, rtol=1e-13, atol=1e-16, method='DOP853')
    x, z, th, pr = sol.y_events[0][0]
    xc = (np.array([x, z]) - Pe) @ ne
    pxc = pr * (np.sin(th) * ne[0] + np.cos(th) * ne[1])
    return xc, pxc, pr - 1.0, sol.t_events[0][0]


def radiation_reaction_k_for_dg(x0, px0, g0, theta, L, dg_target):
    """bend_with_radiation の放射の強さ k を、出口までの δ の変化が dg_target になるように求める (1 次の量なので割線 4 回で十分)"""
    k = 1e-3
    for _ in range(5):
        _, _, dge, _ = bend_with_radiation(x0, px0, g0, theta, L, k)
        k *= dg_target / (dge - g0)
    return k
