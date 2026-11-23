# -*- coding: utf-8 -*-
"""担当 C3 の独立な計算 (pysad の光学・放射・追跡のコードを使わない)。
  Path・lw_field: 円弧と直線をつないだ軌道に対する Liénard-Wiechert の場 (遅延時間を brentq で解く)。RadiationField の検証。
  planar_quad_offset_jx: 平面リングの放射積分 (3x3 行列) に、軌道が四極の中心からずれた分 (h = k x_o、I4 の 2 k² x_o η) を足した Jx。"""
import math

import numpy as np
from scipy.optimize import brentq

C = 299792458.0
E0 = 1.602176634e-19
EPS0 = 8.8541878128e-12


class Path:
    """segs = [('arc', 長さ, 曲率 κ[1/m] (正: +x へ曲がる)), ('line', 長さ)] の列。電子は z 方向 (φ=0) から (x, z) = (0, 0) を出る。
    速度 = βc (sin φ, 0, cos φ)、加速度 = β²c² κ (cos φ, 0, −sin φ)。時刻 t' = s/(βc)"""

    def __init__(self, segs, gamma):
        self.gamma = gamma
        self.beta = math.sqrt(1.0 - 1.0 / gamma ** 2)
        self.segs = []
        s0 = 0.0
        x, z, phi = 0.0, 0.0, 0.0
        for sg in segs:
            kind, ln = sg[0], sg[1]
            kap = sg[2] if kind == 'arc' else 0.0
            self.segs.append((s0, ln, kap, x, z, phi))
            if kap != 0.0:
                phin = phi + kap * ln
                x += (math.cos(phi) - math.cos(phin)) / kap
                z += (math.sin(phin) - math.sin(phi)) / kap
                phi = phin
            else:
                x += ln * math.sin(phi)
                z += ln * math.cos(phi)
            s0 += ln
        self.length = s0

    def state(self, s):
        """s (路程) での位置 (3,)、向き φ、曲率 κ"""
        s = min(max(s, 0.0), self.length)
        for (s0, ln, kap, x0, z0, phi0) in self.segs:
            if s <= s0 + ln or (s0, ln, kap, x0, z0, phi0) == self.segs[-1]:
                d = s - s0
                if kap != 0.0:
                    phi = phi0 + kap * d
                    x = x0 + (math.cos(phi0) - math.cos(phi)) / kap
                    z = z0 + (math.sin(phi) - math.sin(phi0)) / kap
                else:
                    phi = phi0
                    x = x0 + d * math.sin(phi0)
                    z = z0 + d * math.cos(phi0)
                return np.array([x, 0.0, z]), phi, kap
        raise RuntimeError

    def sample(self, m):
        """時刻の等間隔の m 点 (区間の境目は 1 点ずつ。重なり点なし) → (m, 4) = [t, x, y, z]"""
        ss = np.linspace(0.0, self.length, m)
        return ss, np.array([[s / (self.beta * C), *self.state(s)[0]] for s in ss])


def lw_field(path, tau, obs, q=-E0, r0=None):
    """観測時刻 τ (= t_obs − r0/c、r0 = |obs|) での電場 (Liénard-Wiechert、遅延時間を解く)。戻り: (len(tau), 3)"""
    obs = np.array(obs, float)
    if r0 is None:
        r0 = float(np.linalg.norm(obs))
    b = path.beta
    T = path.length / (b * C)
    out = []
    for ta in tau:
        def f(tp):
            R2 = float(np.sum((obs - path.state(b * C * tp)[0]) ** 2))
            Rm = math.sqrt(R2)
            return tp + (R2 - r0 * r0) / (Rm + r0) / C - ta       # R − r0 を桁落ちなしに
        tp = brentq(f, 0.0, T, xtol=1e-30, rtol=1e-15, maxiter=300)
        pos, phi, kap = path.state(b * C * tp)
        Rv = obs - pos
        Rm = float(np.linalg.norm(Rv))
        n = Rv / Rm
        bv = b * np.array([math.sin(phi), 0.0, math.cos(phi)])
        bd = b * b * C * kap * np.array([math.cos(phi), 0.0, -math.sin(phi)])
        kp = 1.0 - n @ bv
        E = (n - bv) / (path.gamma ** 2 * kp ** 3 * Rm ** 2) + np.cross(n, np.cross(n - bv, bd)) / (C * kp ** 3 * Rm)
        out.append(q / (4.0 * math.pi * EPS0) * E)
    return np.array(out)


def planar_quad_offset_jx(dxf, dxd, n_cell=24, nsl=60):
    """平面 FODO (QF D2 B D2 QD DS SD D1 SF DS。QF L=0.3 K1=0.7 (積分値)、QD L=0.3 K1=-0.6、B L=1.2 ANGLE=15 度 E1=E2=ANGLE/2·0.2618/… (試験の台本と同じ)) の
    Jx = 1 − I4/I2。QF・QD の据え付けのずれ dxf・dxd [m] (四極の中心が +x へずれる) による閉軌道のずれが四極の中を通ることで I2・I4 に入る分:
    四極の中の軌道の中心からの位置 x_o、局所の曲率 h = k x_o、I4 の被積分関数 η (h³ + 2 k h) = η (k³x_o³ + 2k²x_o)。
    閉軌道 (δ = 0。RF の周波数は軌道長に合わせたとして、経路長の変化による δ のずれは入れない) は 3x3 のアフィン写像から求める。"""
    kf, kd = 0.7 / 0.3, -0.6 / 0.3
    ang = 0.2617993877991494
    e1 = 0.1308996938995747 * ang
    cell = [('quad', 0.3, kf, dxf), ('drift', 0.5), ('bend', 1.2, ang / 1.2, e1, e1), ('drift', 0.5), ('quad', 0.3, kd, dxd),
            ('drift', 0.1), ('drift', 0.2), ('drift', 1.0), ('drift', 0.2), ('drift', 0.1)]
    P = cell * n_cell

    def drift(l):
        return np.array([[1, l, 0], [0, 1, 0], [0, 0, 1.0]])

    def quad(k, l):
        w = math.sqrt(abs(k))
        if k > 0:
            return np.array([[math.cos(w * l), math.sin(w * l) / w, 0], [-w * math.sin(w * l), math.cos(w * l), 0], [0, 0, 1.0]])
        return np.array([[math.cosh(w * l), math.sinh(w * l) / w, 0], [w * math.sinh(w * l), math.cosh(w * l), 0], [0, 0, 1.0]])

    def bend(l, h):
        kx = h * h
        w = math.sqrt(kx)
        c, s = math.cos(w * l), math.sin(w * l)
        return np.array([[c, s / w, h * (1 - c) / kx], [-w * s, c, h * s / w], [0, 0, 1.0]])

    def edge(h, e):
        return np.array([[1, 0, 0], [math.tan(e) * h, 1, 0], [0, 0, 1.0]])

    def mat(p):
        if p[0] == 'drift':
            return drift(p[1]), np.zeros(3)
        if p[0] == 'quad':
            m = quad(p[2], p[1])
            d = np.array([p[3], 0, 0])
            return m, d - m @ d
        _, L, h, a, b = p
        return edge(h, b) @ bend(L, h) @ edge(h, a), np.zeros(3)
    M = np.eye(3)
    c = np.zeros(3)
    for p in P:
        m, cc = mat(p)
        M = m @ M
        c = m @ c + cc
    eta = np.linalg.solve(np.eye(2) - M[:2, :2], M[:2, 2])
    orb = np.linalg.solve(np.eye(2) - M[:2, :2], c[:2])
    X = np.array([eta[0], eta[1], 1.0])
    O = np.array([orb[0], orb[1], 0.0])
    I2 = I4 = 0.0
    for p in P:
        if p[0] == 'quad':
            L, k, d = p[1], p[2], p[3]
            ds = L / nsl
            m = quad(k, ds)
            for _ in range(nsl):
                O1 = m @ (O - np.array([d, 0, 0])) + np.array([d, 0, 0])
                xo = 0.5 * (O[0] + O1[0]) - d
                X1 = m @ X
                h = k * xo
                I2 += h * h * ds
                I4 += 0.5 * (X[0] + X1[0]) * (h ** 3 + 2 * k * h) * ds
                X, O = X1, O1
        elif p[0] == 'bend':
            _, L, h, a, b = p
            I4 += -X[0] * math.tan(a) * h * h
            X = edge(h, a) @ X
            O = edge(h, a) @ O
            n = int(L * 80)
            ds = L / n
            mb = bend(ds, h)
            for _ in range(n):
                X1 = mb @ X
                O = mb @ O
                I2 += h * h * ds
                I4 += 0.5 * (X[0] + X1[0]) * h * (h * h) * ds
                X = X1
            I4 += -X[0] * math.tan(b) * h * h
            X = edge(h, b) @ X
            O = edge(h, b) @ O
        else:
            X = drift(p[1]) @ X
            O = drift(p[1]) @ O
    return 1.0 - I4 / I2


# ---------------------------------------------------------------------------------------------
# 結合リングの平衡エミッタンス (Chao の 6x6 の方法): 放射の拡散行列と減衰行列から Lyapunov 方程式
# ---------------------------------------------------------------------------------------------
RE = 2.8179403262e-15
LC = 3.8615926796e-13          # ħ/(mc)
MC2 = 0.51099895e6


def coupled_equilibrium(els, E_eV, ndiv=40):
    """平面の偏向磁石 (K1・端の角つき)・四極 (傾きつき)・ドリフト・RF のリングの、結合のある平衡エミッタンス。
    SAD・pysad の光学を使わない。要素の 6x6 行列 (x,x',y,y',z,δ) を自前で作り、放射を 1 周の写像の摂動として扱う:
      減衰: 偏向磁石の刻みごとに Δδ = −ε(X)、Δ(x',y') = −ε0(x',y') (ε = K(1+δ)²[(h+kx)²+(ky)²](1+hx)ds)。M_rad = M0(1 − Σ T⁻¹ G T)
      拡散: 各刻みの光子の放出 δ → δ − u が、その刻みから 1 周の終わりまでの写像 S_i で運ばれる。D = Σ ⟨u²⟩ S_i e_δ e_δᵀ S_iᵀ、
            ⟨u²⟩ = (55/(24√3)) r_e λ̄_C γ⁵ |h|³ ds
      平衡: Σ = M_rad Σ M_radᵀ + D (Kronecker 積で厳密に解く)。固有エミッタンス = Σ·S の固有値の虚部 (S = シンプレクティック行列)。
    戻り: (固有エミッタンス 3 つ (チューンの順に x 様・y 様・z 様)、σδ、U0 [eV])。RF 電圧が足りなければ ValueError。
    RF は 'rf' 要素 (VOLT [V]、HARM) として全周 C に対し FREQ = HARM c/C。"""
    from _a2_indep import _parts, _elem, _edge
    gam = E_eV / MC2
    Kc = (2.0 / 3.0) * RE * gam ** 3                     # 放射の energy loss 係数 (A2 と同じ)
    Ku = 55.0 / (24.0 * math.sqrt(3.0)) * RE * LC * gam ** 5
    P = _parts(els)
    S4 = np.zeros((4, 4))
    S4[0, 1] = S4[2, 3] = 1
    S4 = S4 - S4.T

    def addz(M, hL):
        A, d = M[:4, :4], M[:4, 5]
        for k in range(4):
            M[4, k] = -float(d @ S4 @ A[:, k])
        M[4, 5] = -hL * float(d[0]) / 3.0
        return M
    circ = sum((getattr(x, 'L', 0.0) or 0.0) for x in els)
    steps = []
    I2 = 0.0
    for p in P:
        if p[0] == 'rf':
            steps.append(('rf', p[1], 0.0))
        elif p[0] == 'bend':
            _, L, h, k, e1, e2 = p
            G = np.zeros((6, 6))
            G[5, 0] = -Kc * h * h * math.tan(e1)
            steps.append((_edge(h, e1, 6), G, 0.0))
            ds = L / ndiv
            Mb = addz(_elem('bend', ds, (h, k), 0.0, 6), h * ds)
            eps0 = Kc * h * h * ds
            G = np.zeros((6, 6))
            G[5, 0] = eps0 * (2 * k / h + h)
            G[5, 5] = 2 * eps0
            G[1, 1] = eps0
            G[3, 3] = eps0
            u2 = Ku * abs(h) ** 3 * ds
            for _ in range(ndiv):
                steps.append((Mb, G, u2))
                I2 += h * h * ds
            G = np.zeros((6, 6))
            G[5, 0] = -Kc * h * h * math.tan(e2)
            steps.append((_edge(h, e2, 6), G, 0.0))
        elif p[0] == 'quad':
            steps.append((addz(_elem('quad', p[1], p[2], p[3], 6), 0.0), np.zeros((6, 6)), 0.0))
        else:
            steps.append((addz(_elem('drift', p[1], 0, 0.0, 6), 0.0), np.zeros((6, 6)), 0.0))
    U0 = Kc * I2 * E_eV
    # 1 周の写像 M0 と、各刻みまでの累積 T_i
    M = np.eye(6)
    Ts = []
    for s in steps:
        if isinstance(s[0], str):
            cav = s[1]
            phis = math.asin(U0 / cav.VOLT)
            freq = float(getattr(cav, 'FREQ', 0.0) or 0.0)
            if freq == 0.0:
                freq = float(getattr(cav, 'HARM', 0.0)) * C / circ
            R = np.eye(6)
            R[5, 4] = -cav.VOLT * 2 * math.pi * freq * math.cos(phis) / (E_eV * C)
            M = R @ M
        else:
            M = s[0] @ M
            if np.any(s[1]) or s[2] != 0.0:
                Ts.append((M.copy(), s[1], s[2]))
    Gs = np.zeros((6, 6))
    for T, G, u2 in Ts:
        Gs += np.linalg.inv(T) @ G @ T
    Mrad = M @ (np.eye(6) - Gs)
    D = np.zeros((6, 6))
    e6 = np.zeros(6)
    e6[5] = 1.0
    for T, G, u2 in Ts:
        if u2 == 0.0:
            continue
        v = (M @ np.linalg.inv(T)) @ e6                   # 刻み i の直後から 1 周の終わりまで
        D += u2 * np.outer(v, v)
    # 平衡の分散行列 Σ = Mrad Σ Mradᵀ + D
    A = np.eye(36) - np.kron(Mrad, Mrad)
    Sig = np.linalg.solve(A, D.reshape(-1)).reshape(6, 6)
    Sig = 0.5 * (Sig + Sig.T)
    S6 = np.zeros((6, 6))
    for i in (0, 2, 4):
        S6[i, i + 1] = 1.0
        S6[i + 1, i] = -1.0
    w, _ = np.linalg.eig(Sig @ S6)
    eps = sorted(abs(w.imag))[::2]
    # チューンの順 (x 様・y 様・z 様) に並べる: M0 の固有値の位相から
    ev, _evec = np.linalg.eig(M)
    nus = []
    for z in ev:
        if z.imag > 0:
            nus.append(math.atan2(z.imag, z.real) / 2 / math.pi)
    return eps, math.sqrt(Sig[5, 5]), U0, nus, Sig
