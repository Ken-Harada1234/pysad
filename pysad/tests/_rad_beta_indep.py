# -*- coding: utf-8 -*-
"""β < 1 の放射 (B167) の独立な計算 (pysad の光学・放射・追跡のコードを使わない)。試験 test_v24_D_radbeta.py から使う。
  schott_pn_direct: 第 n 高調波の放射パワー P_n を、円軌道の粒子の遅延 (リエナール・ウィーヘルト) の放射場の時間領域の数値積分から直接求める
                    (Schott の閉じた式 radbeta.schott_s を使わない。観測方向の極角 ψ にわたって積分し、周期の位相積分でフーリエ係数を取る)。
  Ring3B: (x, px, δ) の放射つき 1 周写像。運動方程式を RK4 で積分。放射は粒子ごとの局所の γ_l・β_l の厳密な古典式
          dE = (2/3) r_c β_l³ γ_l⁴ κ² dl [mc²] (κ = 軌道の曲率)、運動量の損失 dp = dE/β_l (dE = v dp)。量子励起は 2 次モーメント
          N<u²> = α λ̄² κ³ r2(γ_l) (55/(24√3)) γ_l⁷ dl を厳密なスペクトル (r2 は引数) から取り、δ の分散 N<u²>/(β_l p0)² とする。
          エネルギー損失は 1 周の終わりに一定の G (空洞の代わり。δ に依らない) で補う: δ の減衰は損失の δ 依存性だけ。共分散は Lyapunov の反復で解く。"""
import math

import numpy as np

ALPHA = 1.0 / 137.035999084
CFAC = 55.0 / (24.0 * math.sqrt(3.0))


def schott_pn_direct(n, beta, nt=4000, npsi=400):
    """第 n 高調波の放射パワー P_n を、時間領域の遅延放射場から (e²/(4πε0 c) = 1、ρ = 1、c = 1 の単位で P_n = ω0² n S_n/β = β n S_n)。
    円軌道 (半径 ρ = 1、ω0 = β、c = 1)。放射場の方向成分 E ∝ n̂×((n̂−β)×β̇)/((1−n̂·β)³ R) の観測時間 t_obs = t' − n̂·r(t') に関するフーリエ係数
    A_n(ψ) = (1/T) ∫ E(t') (1 − n̂·β) e^{i n ω0 t_obs(t')} dt'。単位立体角あたりの第 n 高調波のパワー ∝ |A_n|²・(2 を掛けて正の周波数に)。
    全パワー規格化: 同じ計算の n 総和が Larmor の式 (2/3) β⁴γ⁴ になるよう、定数 K をこの関数の外で決めず、
    ここでは E の振幅を e²/(4π ε0 c) = 1 の単位 (dP/dΩ = |n̂×((n̂−β)×β̇)|²/(4π (1−n̂·β)⁶)・(t_obs の密度 (1−n̂·β)))で作る"""
    # 1 周 (T = 2π/β) の放射。軌道 r = (cos φ, sin φ, 0)、φ = β t'、β_vec = β(−sin φ, cos φ, 0)、β̇ = −β² (cos φ, sin φ, 0)
    phis = (np.arange(nt) + 0.5) * (2.0 * math.pi / nt)
    tp = phis / beta
    rr = np.stack([np.cos(phis), np.sin(phis), np.zeros(nt)], axis=1)
    bv = beta * np.stack([-np.sin(phis), np.cos(phis), np.zeros(nt)], axis=1)
    bd = -beta * beta * rr
    psis, wpsi = np.polynomial.legendre.leggauss(npsi)
    psis = psis * (math.pi / 2.0)
    wpsi = wpsi * (math.pi / 2.0)
    tot = 0.0
    T = 2.0 * math.pi / beta
    dt = T / nt
    for psi, wp in zip(psis, wpsi):
        nh = np.array([math.cos(psi), 0.0, math.sin(psi)])
        nb = nh[None, :] - bv
        v = np.cross(nh[None, :], np.cross(nb, bd))
        kap = 1.0 - bv @ nh
        tobs = tp - rr @ nh
        # E(t') (1 − n̂·β) dt'/T の n 次フーリエ係数。位相 exp(i n ω0 t_obs)、ω0 = β
        ph = np.exp(1j * n * beta * tobs)
        a = ((v * (1.0 / kap**2)[:, None]) * ph[:, None]).sum(axis=0) * dt / T
        # dP/dΩ の第 n 高調波 (正の周波数側 2|A_n|²)、立体角 dΩ = cos ψ dψ dφ_obs (φ_obs の積分は 2π: 軸対称)
        tot += wp * math.cos(psi) * 2.0 * (np.abs(a)**2).sum()
    return tot * 2.0 * math.pi / (4.0 * math.pi)      # dP/dΩ の 1/(4π) (e²/(4πε0 c) = 1 の単位)。結果は β n S_n (ω0 = β、ρ = 1、c = 1 で P_n = ω0² n S_n/β)


class Ring4B:
    """(x, px, δ, z) の放射つき 1 周写像 + 薄い空洞 (先頭)。z = β0 c Δt (時間遅れ)。空洞は Δδ = V sin(φs + k z)/(p0 mc²) (V [eV]、K64 の cavity と同じく Δp c = ΔE)。
    elems: dict(type in drift|quad|bend, L, h (曲率、bend)、k1)。p0 = βγ (mc 単位)、rc = 古典半径 [m]、lam = ħ/(mc) [m]、r2fun(γ) = 厳密な 2 次モーメント比"""

    def __init__(self, p0, rc, lam, elems, r2fun, V_over_mc2, harm, ds=0.02):
        self.p0 = p0
        self.g0 = math.sqrt(1.0 + p0 * p0)
        self.b0 = p0 / self.g0
        self.rc = rc
        self.lam = lam
        self.r2fun = r2fun
        self.Vn = V_over_mc2
        self.S = []
        for e in elems:
            L = e['L']
            if e['type'] == 'drift':
                self.S.append(e)
                continue
            n = max(1, int(math.ceil(L / ds)))
            for _ in range(n):
                hh = e.get('h', 0.0) if e['type'] == 'bend' else 0.0
                self.S.append(dict(type=e['type'], L=L / n, h=hh, hf=hh, k1=e.get('k1', 0.0)))
        self.C = sum(e['L'] for e in elems)
        self.k = 2.0 * math.pi * harm / self.C
        self.phis = 0.0
        self.rf = True

    def _loc(self, d):
        p = self.p0 * (1.0 + d)
        g = math.sqrt(1.0 + p * p)
        return p, g, p / g

    def _f(self, y, s):
        x, px, d, z = y
        h = s['h']
        F = s['hf'] + s['k1'] * x
        pz = math.sqrt((1.0 + d) ** 2 - px * px)
        dl = (1.0 + h * x) * (1.0 + d) / pz
        kap = F / (1.0 + d)
        p, g, b = self._loc(d)
        dE = (2.0 / 3.0) * self.rc * b ** 3 * g ** 4 * kap * kap * dl       # [mc²]
        dd = -dE / (b * self.p0)
        return np.array([(1.0 + h * x) * px / pz, h * pz - (1.0 + h * x) * F + px * dd / (1.0 + d), dd,
                         self.b0 / b * dl - 1.0]), dl, F

    def slice(self, y, s):
        L = s['L']
        if s['type'] == 'drift':
            x, px, d, z = y
            pz = math.sqrt((1.0 + d) ** 2 - px * px)
            dl = (1.0 + d) / pz
            p, g, b = self._loc(d)
            return np.array([x + L * px / pz, px, d, z + L * (self.b0 / b * dl - 1.0)]), 0.0
        if s['type'] == 'quad':
            s = dict(s, h=0.0, hf=0.0)
        f = lambda v: self._f(v, s)[0]
        k1 = f(y)
        k2 = f(y + L / 2 * k1)
        k3 = f(y + L / 2 * k2)
        k4 = f(y + L * k3)
        yn = y + L / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ym = 0.5 * (y + yn)
        _, dl, F = self._f(ym, s)
        kap = F / (1.0 + ym[2])
        p, g, b = self._loc(ym[2])
        m2 = self.rc * self.lam * abs(kap) ** 3 * self.r2fun(g) * CFAC * g ** 7 * dl * L
        return yn, m2 / (b * self.p0) ** 2

    def turn(self, y, rf=True):
        if rf and self.Vn == 0.0:
            y = self._noscale_turn(y)
            y = y.copy()
            y[2] += self.G
            return y
        if rf:
            y = y.copy()
            y[2] = self._cavity(y[2], math.sin(self.phis + self.k * y[3]))
        for s in self.S:
            y, _ = self.slice(y, s)
        return y

    def _cavity(self, d, sinphi):
        """薄い空洞: 粒子のエネルギーが Vn sin φ [mc²] 増える (空洞の仕事)。運動量は E dE = p dp から (Δp = ΔE/(β c)。B232 の独立な基準)。
        energy=False (K64 の旧式) は Δp c = ΔE"""
        if not getattr(self, 'energy_kick', True):
            return d + self.Vn / self.p0 * sinphi
        p = self.p0 * (1.0 + d)
        g = math.sqrt(1.0 + p * p) + self.Vn * sinphi
        return math.sqrt(g * g - 1.0) / self.p0 - 1.0

    def _noscale_turn(self, y):
        for s in self.S:
            y, _ = self.slice(y, s)
        return y

    def solve_rfless(self):
        """空洞なし (Vn = 0): 1 周の損失 G (δ に依らない定数) を終わりに足して閉軌道を解く。縦は δ の減衰 (損失の δ 依存) だけ。3x3 (x, px, δ)"""
        self.G = 0.0
        y = np.zeros(4)
        y1 = self._noscale_turn(y)
        self.G = -(y1[2] - y[2])
        self.u0_delta = self.G
        for _ in range(30):
            f = (self.turn(y) - y)[:3]
            if np.max(np.abs(f)) < 1e-18:
                break
            J = self.jac(self.turn, y)[:3, :3]
            y[:3] -= np.linalg.solve(J - np.eye(3), f)
        self.orbit = y
        self.M = self.jac(self.turn, y)
        return np.linalg.eigvals(self.M[:3, :3])

    def equilibrium_rfless(self):
        y = self.orbit.copy()
        Qt = np.zeros((3, 3))
        for s in self.S:
            J = self.jac(lambda v: self.slice(v, s)[0], y)[:3, :3]
            yn, var = self.slice(y, s)
            w = np.array([0.0, y[1] / (1.0 + y[2]), 1.0])
            Qt = J @ Qt @ J.T + var * np.outer(w, w)
            y = yn
        M = self.M[:3, :3]
        w, V = np.linalg.eig(M)
        Vi = np.linalg.inv(V)
        Qc = Vi @ Qt @ Vi.conj().T
        return np.real(V @ (Qc / (1.0 - np.outer(w, w.conj()))) @ V.conj().T)

    def jac(self, fn, y, h=1e-6):
        n = len(y)
        J = np.zeros((n, n))
        for i in range(n):
            d = np.zeros(n)
            d[i] = h
            J[:, i] = (fn(y + d) - fn(y - d)) / (2 * h)
        return J

    def solve(self):
        y = np.zeros(4)
        for _ in range(20):
            f = self.turn(y, rf=False) - y
            J = np.zeros((2, 2))
            for i in range(2):
                d = np.zeros(4)
                d[i] = 1e-7
                J[:, i] = ((self.turn(y + d, False) - self.turn(y - d, False)) / 2e-7)[:2]
            if np.max(np.abs(f[:2])) < 1e-18:
                break
            y[:2] -= np.linalg.solve(J - np.eye(2), f[:2])
        y1 = self.turn(y, rf=False)
        self.u0_delta = -(y1[2] - y[2])
        eg = lambda d: math.sqrt(1.0 + (self.p0 * (1.0 + d)) ** 2)
        self.u0_energy = eg(y[2]) - eg(y1[2])                      # 1 周のエネルギー損失 [mc²]
        a0 = math.asin(self.u0_energy / self.Vn) if getattr(self, 'energy_kick', True) else math.asin(self.u0_delta * self.p0 / self.Vn)
        y0 = y.copy()
        for cand in (a0, math.pi - a0):                    # 安定な位相 (遷移エネルギーの上下で違う): 縦の固有値が複素で |λ| < 1 の方
            self.phis = cand
            y = y0.copy()
            for _ in range(40):
                f = self.turn(y) - y
                if np.max(np.abs(f)) < 1e-17:
                    break
                y = y - np.linalg.solve(self.jac(self.turn, y) - np.eye(4), f)
            ev = np.linalg.eigvals(self.jac(self.turn, y))
            if np.max(np.abs(ev)) < 1.0:
                break
        self.orbit = y
        self.M = self.jac(self.turn, y)
        return np.linalg.eigvals(self.M)

    def equilibrium(self):
        """1 周の共分散 Σ = M Σ Mᵀ + Q を、固有分解による Lyapunov (M の固有値の絶対値 < 1) で解く (閉じた解)"""
        y = self.orbit.copy()
        y[2] = self._cavity(y[2], math.sin(self.phis + self.k * y[3]))
        Qt = np.zeros((4, 4))
        for s in self.S:
            J = self.jac(lambda v: self.slice(v, s)[0], y)
            yn, var = self.slice(y, s)
            w = np.array([0.0, y[1] / (1.0 + y[2]), 1.0, 0.0])
            Qt = J @ Qt @ J.T + var * np.outer(w, w)
            y = yn
        # Σ = Σ_k M^k Q (M^T)^k: M = V Λ V⁻¹ → Σ = V [ (V⁻¹ Q V⁻ᴴ)_ij / (1 − λi λj*) ] Vᴴ
        w, V = np.linalg.eig(self.M)
        Vi = np.linalg.inv(V)
        Qc = Vi @ Qt @ Vi.conj().T
        Sc = Qc / (1.0 - np.outer(w, w.conj()))
        self._modal = (w, V, Sc)
        return np.real(V @ Sc @ V.conj().T)

    def mode_emittance(self, tune_lo, tune_hi):
        """モード (固有位相の絶対値が tune_lo..tune_hi の共役対) だけの共分散の (x, px) ブロックの sqrt(det) = そのモードの水平エミッタンス。
        δ・z の分散との相関を引く条件つき共分散 (病的に条件が悪い逆行列) を使わない"""
        w, V, Sc = self._modal
        ph = np.abs(np.angle(w)) / (2 * math.pi)
        idx = [i for i in range(len(w)) if tune_lo < ph[i] < tune_hi]
        Sx = np.zeros((4, 4), dtype=complex)
        for i in idx:
            for j in idx:
                Sx += Sc[i, j] * np.outer(V[:, i], V[:, j].conj())
        Sx = np.real(Sx)[:2, :2]
        return math.sqrt(np.linalg.det(Sx))


def beam_sizes(S):
    """σxx、σδδ、σzz、σxδ"""
    return S[0, 0], S[2, 2], S[3, 3], S[0, 2]
