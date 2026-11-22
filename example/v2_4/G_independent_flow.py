# -*- coding: utf-8 -*-
"""v2.4.0 G: K0≠0 の偏向磁石の K64 の追跡の式 (tbendcore1・k_steer) を、K64 の式を使わない独立な基準と比べる。
基準 = 実験室の直交座標 (X, Y, Z) での相対論的な荷電粒子の運動方程式 (dr/ds = p/|p|、dp/ds = (p/|p|) × b、p は p0 で規格化、b = qB/p0 [1/m]) の数値積分 (scipy DOP853、rtol 1e-13)。
  磁場: 一様な双極子 b = (ANGLE + K0)/L の領域を 2 枚の面 (入口の面 E1・出口の面 E2、tanh の端の幅 w) で切った場 G(X, Z) = Θ(s1) Θ(−s2)。
  縦方向 Y の依存はスカラーポテンシャル φ = −b [Y G − Y³ ∇²G/6] (∇²: X-Z 面) から Maxwell の式を満たす場 (Bx・By・Bz とも Y³ まで)。
  w → 0 で K64 の硬い端 (楔 E1・E2 + Maxwell の項) の極限になるはずで、w = 4e-3・2e-3・1e-3 の差が w に比例して 0 に向かうことを見る。
  参照の軌道 (幾何) は ANGLE の弧 (半径 ρ0 = L/ANGLE)、場の曲がりは b。入口・出口の座標は参照軌道に垂直な面 (入口 Z = 0、出口は O2 を通る面)。
  K1 = 0 の磁石だけ (K1 の場の独立な基準は未作成)。使い方: py -3.12 example/v2_4/G_independent_flow.py"""
import math
import os
import sys
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)


def _th(u, k):
    """Θ(u) = (1 + tanh(u))/2 の k 階微分"""
    t = np.tanh(u)
    if k == 0:
        return 0.5 * (1 + t)
    s2 = 1 - t * t
    if k == 1:
        return 0.5 * s2
    if k == 2:
        return -s2 * t
    if k == 3:
        return -s2 * (1 - 3 * t * t)
    raise ValueError(k)


class Field:
    """面 1・2 (法線 n1・n2 は X-Z 面内の単位ベクトル、s1 = n1·(r − r1)・s2 = n2·(r − r2)) で切った一様な双極子 b の Maxwell 場 (Y³ まで)"""

    def __init__(self, b, n1, r1, n2, r2, w):
        self.b, self.n1, self.r1, self.n2, self.r2, self.w = b, np.array(n1), np.array(r1), np.array(n2), np.array(r2), w
        self.c12 = float(np.dot(self.n1, self.n2))

    def terms(self, r):
        w = self.w
        u1 = float(np.dot(self.n1, [r[0] - self.r1[0], r[2] - self.r1[1]])) / w
        u2 = float(np.dot(self.n2, [r[0] - self.r2[0], r[2] - self.r2[1]])) / w
        a = [_th(u1, k) / w ** k for k in range(4)]
        c = [_th(-u2, k) * (-1) ** k / w ** k for k in range(4)]       # Θ(−s2) の s2 微分
        return a, c

    def B(self, r):
        a, c = self.terms(r)
        n1, n2, c12 = self.n1, self.n2, self.c12
        G = a[0] * c[0]
        Gx = a[1] * n1[0] * c[0] + a[0] * c[1] * n2[0]
        Gz = a[1] * n1[1] * c[0] + a[0] * c[1] * n2[1]
        # ∇²(a_p c_q) = a_{p+2} c_q + 2 c12 a_{p+1} c_{q+1} + a_p c_{q+2}
        lap = a[2] * c[0] + 2 * c12 * a[1] * c[1] + a[0] * c[2]
        # ∂x ∇²G = Σ: ∂x(a_p c_q) = a_{p+1} n1x c_q + a_p c_{q+1} n2x
        lapx = (a[3] * c[0] + 2 * c12 * a[2] * c[1] + a[1] * c[2]) * n1[0] + (a[2] * c[1] + 2 * c12 * a[1] * c[2] + a[0] * c[3]) * n2[0]
        lapz = (a[3] * c[0] + 2 * c12 * a[2] * c[1] + a[1] * c[2]) * n1[1] + (a[2] * c[1] + 2 * c12 * a[1] * c[2] + a[0] * c[3]) * n2[1]
        Y = r[1]
        b = self.b
        return np.array([b * (Y * Gx - Y ** 3 / 6.0 * lapx), b * (G - Y * Y / 2.0 * lap), b * (Y * Gz - Y ** 3 / 6.0 * lapz)])


def geometry(el, sig=(1.0, 1.0)):
    """入口の面・出口の面と出口の座標系 (O2, t2, x̂2)。ANGLE の弧 (X < 0 側へ曲がる)。E1・E2 の符号は sig で選ぶ (sig = ±1)"""
    L, phi = float(el.L), float(el.ANGLE)
    psi1, psi2 = float(getattr(el, 'E1', 0.0) or 0.0), float(getattr(el, 'E2', 0.0) or 0.0)
    if phi != 0.0:
        rho0 = L / phi
        O2 = np.array([-rho0 * (1 - math.cos(phi)), 0.0, rho0 * math.sin(phi)])
        t2 = np.array([-math.sin(phi), 0.0, math.cos(phi)])
        x2 = np.array([math.cos(phi), 0.0, math.sin(phi)])
    else:
        O2 = np.array([0.0, 0.0, L])
        t2 = np.array([0.0, 0.0, 1.0])
        x2 = np.array([1.0, 0.0, 0.0])
    # 面の法線 (X-Z 面): 入口は進行方向 ẑ を sig1 ψ1 回した向き (内側が s1 > 0)、出口は t2 を sig2 ψ2 回した向き (内側が s2 < 0)
    n1 = np.array([sig[0] * math.sin(psi1), math.cos(psi1)])
    n2 = np.array([sig[1] * math.sin(psi2) * x2[0] + math.cos(psi2) * t2[0], sig[1] * math.sin(psi2) * x2[2] + math.cos(psi2) * t2[2]])
    return O2, t2, x2, n1, n2


def lorentz_map(el, Z, delta=0.0, w=1e-3, sig=(1.0, 1.0), dm=None):
    """el の K64 と同じ座標 (x, px, y, py) の入力 Z (n×4、正準運動量) を、運動方程式の数値積分で出口の座標にする (n×4)"""
    L, phi, k0 = float(el.L), float(el.ANGLE), float(el.K0)
    b = (phi + k0) / L
    O2, t2, x2, n1, n2 = geometry(el, sig)
    f = Field(b, n1, [0.0, 0.0], n2, [O2[0], O2[2]], w)
    dm = 30 * w if dm is None else dm
    out = []
    for z in np.atleast_2d(Z):
        x, px, y, py = z
        pr = 1.0 + delta
        pz = math.sqrt(pr * pr - px * px - py * py)
        p = np.array([px, py, pz])                     # (X, Y, Z) 成分 (入口の座標系 = 実験室)
        r = np.array([x, y, 0.0]) - p / pz * dm        # 面の前の Z = −dm まで自由に戻る

        def rhs(s, v):
            pp = v[3:]
            tt = pp / np.linalg.norm(pp)
            return np.concatenate([tt, np.cross(tt, f.B(v[:3]))])

        def ev(s, v):
            return float(np.dot(v[:3] - O2, t2)) - dm
        ev.terminal = True
        ev.direction = 1
        sol = solve_ivp(rhs, [0.0, 4 * L + 3 * dm], np.concatenate([r, p]), method='DOP853', rtol=1e-13, atol=1e-16, events=ev)
        v = sol.y[:, -1]
        pos, pm = v[:3], v[3:]
        tt = pm / np.linalg.norm(pm)
        ds = -float(np.dot(pos - O2, t2)) / float(np.dot(tt, t2))
        pos = pos + tt * ds
        d = pos - O2
        out.append([float(np.dot(d, x2)), float(np.dot(pm, x2)), d[1], pm[1]])
    return np.array(out)


def compare(el, pts, delta=0.0, ws=(4e-3, 2e-3, 1e-3), sig=(1.0, 1.0)):
    from pysad import ampnf_k0 as K
    ref = K.numeric_map(el, delta)(np.array(pts))
    res = []
    for w in ws:
        got = lorentz_map(el, pts, delta, w, sig)
        res.append(float(np.abs(got - ref).max()))
    return res, ref


def P(*a):
    print(*a, flush=True)


def main():
    from pysad.elements import Bend
    from pysad import compat
    mode = sys.argv[1] if len(sys.argv) > 1 else 'v2'
    compat.set_mode(mode)
    P('互換モード:', compat.mode())
    cases = {'steer': dict(L=0.3, K0=0.06), 'steer_edge': dict(L=0.3, K0=0.06, E1=0.1, E2=-0.05),
             'bend': dict(L=0.5, ANGLE=0.1, K0=0.03), 'bend_edge': dict(L=0.5, ANGLE=0.1, K0=-0.02, E1=0.04, E2=0.03),
             'reverse': dict(L=0.5, ANGLE=-0.1, K0=0.05, E1=0.04, E2=0.03)}
    rng = np.random.default_rng(5)
    for nm, kw in cases.items():
        el = Bend(name='B', **kw)
        for delta in (0.0, 2e-3):
            # 水平面 (y = 0): 面の符号 sig を、v2 (E2 を常に入れる) の線形の項 (小さな振幅) で決める。k64 の steer は出口の楔を入れないので sig の決定には使わない
            compat.set_mode('v2')
            best = None
            for sig in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                e, _ = compare(el, [[1e-3, 2e-3, 0.0, 0.0]], delta, ws=(2e-5,), sig=sig)
                if best is None or e[0] < best[0]:
                    best = (e[0], sig)
            sig = best[1]
            compat.set_mode(mode)
            pts_h = np.column_stack([rng.uniform(-3e-3, 3e-3, 6), rng.uniform(-3e-3, 3e-3, 6), np.zeros(6), np.zeros(6)])
            eh, _ = compare(el, pts_h, delta, ws=(2e-5, 1e-5), sig=sig)
            pts_v = np.column_stack([rng.uniform(-1e-3, 1e-3, 4), rng.uniform(-1e-3, 1e-3, 4), rng.uniform(-1.5e-4, 1.5e-4, 4), rng.uniform(-1e-3, 1e-3, 4)])
            ev, ref = compare(el, pts_v, delta, ws=(4e-3, 2e-3, 1e-3), sig=sig)
            scale = float(np.abs(ref).max())
            P('%-11s δ=%-6g 面の符号 %s: 水平面 (y=0, w=2e-5・1e-5) の差 %.1e・%.1e、(x, px, y, py の 4 次元、w = 4e-3・2e-3・1e-3) の差 %.1e・%.1e・%.1e (出力の最大 %.1e)'
              % (nm, delta, sig, eh[0], eh[1], ev[0], ev[1], ev[2], scale))


if __name__ == '__main__':
    main()


# ---------------------------------------------------------------------------------------------
# K1≠0 の磁石の本体 (水平面 y = 0、端の面なし = E1 = E2 = 0・DISFRIN = 1 でフリンジなし) の独立な基準
#   曲がった磁石 (ANGLE≠0): 幾何の曲率中心のまわりの軸対称な場 B_Y(r) = b0 + K1 ρ0 ln(r/ρ0) ('log' = (1+hx) ∂b/∂x が一様。K64 の暗黙の定義)
#     または b0 + K1 (r − ρ0) ('lin' = B' 一様)。r = |位置 − 中心|。b0 = (ANGLE + K0)/L。入口 (Z = 0) から出口の面 (O2 を通り t2 に垂直) まで場が一様にある。
#   直線の磁石 (ANGLE = 0): B_Y = b0 + K1 X (b0 = K0/L)。
# ---------------------------------------------------------------------------------------------
def lorentz_body_planar(el, Z, delta=0.0, kdef='log'):
    L, phi, k0, k1 = float(el.L), float(el.ANGLE), float(el.K0), float(getattr(el, 'K1', 0.0) or 0.0)
    if phi != 0.0:
        rho0 = L / phi
        b0 = (phi + k0) / L
        O2 = np.array([-rho0 * (1 - math.cos(phi)), 0.0, rho0 * math.sin(phi)])
        t2 = np.array([-math.sin(phi), 0.0, math.cos(phi)])
        x2 = np.array([math.cos(phi), 0.0, math.sin(phi)])
        cen = np.array([-rho0, 0.0])

        def by(X, Zc):
            r = math.hypot(X - cen[0], Zc - cen[1])
            return b0 + (k1 * rho0 * math.log(r / abs(rho0)) if kdef == 'log' else k1 * (r - abs(rho0)) * (1 if rho0 > 0 else -1))
    else:
        b0 = k0 / L
        O2 = np.array([0.0, 0.0, L])
        t2 = np.array([0.0, 0.0, 1.0])
        x2 = np.array([1.0, 0.0, 0.0])

        def by(X, Zc):
            return b0 + k1 * X
    out = []
    for z in np.atleast_2d(Z):
        x, px, y, py = z
        pr = 1.0 + delta
        pz = math.sqrt(pr * pr - px * px - py * py)
        v0 = np.array([x, y, 0.0, px, py, pz])

        def rhs(s, v):
            pp = v[3:]
            tt = pp / np.linalg.norm(pp)
            B = np.array([0.0, by(v[0], v[2]), 0.0])
            return np.concatenate([tt, np.cross(tt, B)])

        def ev(s, v):
            return float(np.dot(v[:3] - O2, t2))
        ev.terminal = True
        ev.direction = 1
        sol = solve_ivp(rhs, [0.0, 4 * L], v0, method='DOP853', rtol=1e-13, atol=1e-16, events=ev)
        v = sol.y[:, -1]
        pos, pm = v[:3], v[3:]
        d = pos - O2
        out.append([float(np.dot(d, x2)), float(np.dot(pm, x2)), d[1], pm[1]])
    return np.array(out)


def compare_body(el, pts, delta=0.0, kdef='log'):
    from pysad import ampnf_k0 as K
    ref = K.numeric_map(el, delta)(np.array(pts))
    return float(np.abs(lorentz_body_planar(el, pts, delta, kdef) - ref).max()), ref


# ---------------------------------------------------------------------------------------------
# ANGLE = 0 で K0≠0 かつ K1≠0 の直線の磁石 (端 E = 0。面は Z = 0 と Z = L、端の幅 w)。Maxwell を満たす場 (Y・X の 3 次まで):
#   スカラーポテンシャル φ = −b0 (Y G − Y³ G''/6) − k1 [X Y G − (X³ Y + X Y³) G''/12] (G = G(Z))。B = −∇φ。
# ---------------------------------------------------------------------------------------------
def lorentz_straight_k1(el, Z, delta=0.0, w=1e-3):
    L, k0, k1 = float(el.L), float(el.K0), float(el.K1)
    b0 = k0 / L
    dm = 30 * w

    def G(zc):
        a = [_th((zc) / w, k) / w ** k for k in range(5)] if False else None
        return None

    def fld(r):
        X, Y, Zc = r
        u1, u2 = Zc / w, -(Zc - L) / w
        a = [_th(u1, k) / w ** k for k in range(4)]
        c = [_th(u2, k) * (-1) ** k / w ** k for k in range(4)]
        # G = a0 c0 の Z 微分 (n 階): Σ C(n, j) a_j c_{n-j}
        g = [sum(math.comb(n, j) * a[j] * c[n - j] for j in range(n + 1)) for n in range(4)]
        Bx = b0 * 0.0 + k1 * (Y * g[0] - (3 * X * X * Y + Y ** 3) * g[2] / 12.0)
        By = b0 * (g[0] - Y * Y * g[2] / 2.0) + k1 * (X * g[0] - (X ** 3 + 3 * X * Y * Y) * g[2] / 12.0)
        Bz = b0 * (Y * g[1] - Y ** 3 * g[3] / 6.0) + k1 * (X * Y * g[1] - (X ** 3 * Y + X * Y ** 3) * g[3] / 12.0)
        return np.array([Bx, By, Bz])
    out = []
    for z in np.atleast_2d(Z):
        x, px, y, py = z
        pr = 1.0 + delta
        pz = math.sqrt(pr * pr - px * px - py * py)
        p = np.array([px, py, pz])
        r = np.array([x, y, 0.0]) - p / pz * dm

        def rhs(s, v):
            pp = v[3:]
            tt = pp / np.linalg.norm(pp)
            return np.concatenate([tt, np.cross(tt, fld(v[:3]))])

        def ev(s, v):
            return v[2] - (L + dm)
        ev.terminal = True
        ev.direction = 1
        sol = solve_ivp(rhs, [0.0, 4 * L + 3 * dm], np.concatenate([r, p]), method='DOP853', rtol=1e-13, atol=1e-16, events=ev)
        v = sol.y[:, -1]
        pos, pm = v[:3], v[3:]
        tt = pm / np.linalg.norm(pm)
        pos = pos + tt * ((L - pos[2]) / tt[2])
        out.append([pos[0], pm[0], pos[1], pm[1]])
    return np.array(out)


# ---------------------------------------------------------------------------------------------
# K1≠0 の曲がった磁石 (ANGLE≠0) の全 4 次元・端つき (楔 E1・E2、Maxwell の端、四極のフリンジ) の独立な基準
#   曲率中心 C = (−ρ0, 0) のまわりの軸対称な中心面の場 b(r) ('log': b0 + K1 ρ0 ln(r/ρ0)、'lin': b0 + K1 (r − ρ0)) に、入口・出口の面 (法線 n1・n2、tanh の端の幅 w) の
#   関数 G(X, Z) = Θ(s1) Θ(−s2) を掛け、スカラーポテンシャル φ = Y φ0 + Y³ φ1 (φ0 = −b G、φ1 = Δ⊥(b G)/6、Δ⊥ = ∂X² + ∂Z²。∇²φ = 0 を Y³ まで満たす) から B = −∇φ。
#   ∇φ は複素ステップ (厳密) で求める。ρ0 = L/ANGLE、b0 = (ANGLE + K0)/L。
# ---------------------------------------------------------------------------------------------
def lorentz_curved_k1(el, Z, delta=0.0, w=1e-3, sig=(-1.0, 1.0), kdef='log', dm=None, lw=True):
    L, phi_a, k0, k1 = float(el.L), float(el.ANGLE), float(el.K0), float(getattr(el, 'K1', 0.0) or 0.0)
    rho0 = L / phi_a
    b0 = (phi_a + k0) / L
    O2, t2, x2, n1, n2 = geometry(el, sig)
    cX, cZ = -rho0, 0.0
    r1 = (0.0, 0.0)
    r2 = (O2[0], O2[2])
    c12 = float(np.dot(n1, n2))
    dm = 30 * w if dm is None else dm

    def theta(u, k):
        t = np.tanh(u)
        if k == 0:
            return 0.5 * (1 + t)
        s2 = 1 - t * t
        return [None, 0.5 * s2, -s2 * t, -s2 * (1 - 3 * t * t)][k]

    def phi0(X, Zc):
        """中心面の場の関数 −B_Y(X, 0, Z) = φ0。lw=True: 多極 (調和) 展開の端 (Lee–Whiting の前提) = 直線の四極の B_Y = k1 (X G − X³ G''/12) と同じ x³ ΔG の項を足す
        (x = r − ρ0 は参照の円からの動径方向の距離)。lw=False: 中心面の場 = b(r) G (階段) だけ"""
        dX, dZ = X - cX, Zc - cZ
        r = np.sqrt(dX * dX + dZ * dZ)
        u1 = (n1[0] * (X - r1[0]) + n1[1] * (Zc - r1[1])) / w
        u2 = (n2[0] * (X - r2[0]) + n2[1] * (Zc - r2[1])) / w
        a = [theta(u1, k) / w ** k for k in range(3)]
        c = [theta(-u2, k) * (-1) ** k / w ** k for k in range(3)]
        G = a[0] * c[0]
        lapG = a[2] * c[0] + 2 * c12 * a[1] * c[1] + a[0] * c[2]
        if kdef == 'log':
            b = b0 + k1 * rho0 * np.log(r / abs(rho0))
        else:
            b = b0 + k1 * (r - abs(rho0)) * (1.0 if rho0 > 0 else -1.0)
        val = -b * G
        if lw:
            xr = r - abs(rho0)
            val = val + (k1 / 12.0) * xr ** 3 * lapG
        return val

    def potential(X, Y, Zc):
        h = 2e-5
        p0 = phi0(X, Zc)
        lap = (phi0(X + h, Zc) + phi0(X - h, Zc) + phi0(X, Zc + h) + phi0(X, Zc - h) - 4 * p0) / (h * h)
        return Y * p0 - Y ** 3 * lap / 6.0

    def fld(r):
        h = 1e-30
        out = np.empty(3)
        for k in range(3):
            q = [complex(v) for v in r]
            q[k] = complex(r[k], h)
            out[k] = -potential(*q).imag / h
        return out
    res = []
    for z in np.atleast_2d(Z):
        x, px, y, py = z
        pr = 1.0 + delta
        pz = math.sqrt(pr * pr - px * px - py * py)
        p = np.array([px, py, pz])
        r = np.array([x, y, 0.0]) - p / pz * dm

        def rhs(s, v):
            pp = v[3:]
            tt = pp / np.linalg.norm(pp)
            return np.concatenate([tt, np.cross(tt, fld(v[:3]))])

        def ev(s, v):
            return float(np.dot(v[:3] - O2, t2)) - dm
        ev.terminal = True
        ev.direction = 1
        sol = solve_ivp(rhs, [0.0, 4 * L + 3 * dm], np.concatenate([r, p]), method='DOP853', rtol=1e-12, atol=1e-15, events=ev)
        v = sol.y[:, -1]
        pos, pm = v[:3], v[3:]
        tt = pm / np.linalg.norm(pm)
        pos = pos + tt * (-float(np.dot(pos - O2, t2)) / float(np.dot(tt, t2)))
        d = pos - O2
        res.append([float(np.dot(d, x2)), float(np.dot(pm, x2)), d[1], pm[1]])
    return np.array(res)
