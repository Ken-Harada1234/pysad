# -*- coding: utf-8 -*-
"""v2.4.1: ampnf の限界を広げる拡張 (5・6 次の共鳴の駆動項の表・1 本の共鳴の島)。`ampnf.py` の正規形の枠組みの延長。解説は docs/49 §15。

【J² と 5・6 次の駆動項】 J² (∂²ν/∂J²) は `ampnf.amplitude_tune` の D ≥ 5 で既に入る (`r.adts2`、`r.d2nu(plane)`、`r.nu(..., jorder=)`)。J³ は D ≥ 7、J^k は D ≥ 2k+1。
  5・6 次の共鳴 (nν = 整数、n = 5・6) の駆動項は、1 周の写像の次数 n−1 の単項式 ζ̄^(n−1) の係数 (D ≥ n−1)。Birkhoff の正規形で、その共鳴の単項式を消さずに残すと
  (`resonance_normal_form`・`drive_table`)、残った係数 r が、その共鳴の駆動項 (低い次数の非共鳴の項の変換との合成を含む)。これは写像の単項式の係数で、
  Hamiltonian の h_jklm と同じ量ではない。強さのべきでは、5 次が 6 極 K2 の 1〜5 次・8 極との積、6 次がその 6 次までの積。
【島 (1 本の共鳴の包絡)】 共鳴 k·ν = p の 1 本だけを残した正規形から、振り子の Hamiltonian H = H0(J) + A cos(k·θ + ψ) (H0' = 2π ν(J)) の島の幅を出す (`island_theory`)。
  1 周の写像の ζ' = λ(J) ζ + r ζ^α (α = 共鳴の単項式) の作用の増分 ΔJ_i = k_i A sinψ (k = 共鳴のベクトル) から A、M = kᵀ (2π ∂ν/∂J) k から
  島の半幅 ΔJ = 2|k| √(A/|M|) (k の方向の作用)。固定点の作用 J_r は ν(J_r) の級数 (J²、…) で k·ν = p を解く。
  入れた効果: 正規形の次数 D までの全ての非共鳴の項 (ν(J) の高次を含む)、共鳴の駆動項は最低次 (α の次数の単項式 1 個、(ζζ̄)^m 倍の補正は無視)。
  入れていない効果: 他の共鳴との重なり (カオス層)、駆動項の J 依存、島の非対称、δ ≠ 0、2 本目の共鳴。検証は結合なし 1 自由度の水平面 (純粋な x の共鳴) のみ
  (垂直面・和差の共鳴は式だけで追跡と比べていない)。**リングの対称性で禁止された共鳴 (超周期 P の lat で nν が P の倍数でない) にも、駆動項を線形のチューンで評価するので
  有限の幅を返す。対称性のある lat では `drive_table` の |r| を対称性を破った lat と比べて先に確かめる。**
  `pendulum_ok` (島の半幅 / J_r < 0.3) が False のときは幅は桁の目安 (振り子の近似の外)。
【検証用】 `track_island`: 追跡 (pysad k64、RF・放射なし) の周期 n の固定点 O・X と、X から出る不安定多様体 (分離線) の正規形の作用の範囲 (1 自由度)。
"""
import math

import numpy as np

from . import ampnf as A


def resonance_keep(plane, alpha):
    """plane (0 = x、1 = y) の方程式の単項式 alpha = (a, b, c, d) (ζx^a ζ̄x^b ζy^c ζ̄y^d) とその複素共役の方程式の単項式の集合 (normal_form の keep)"""
    a, b, c, d = alpha
    return {(2 * plane, (a, b, c, d, 0)), (2 * plane + 1, (b, a, d, c, 0))}


def pure_alpha(n, plane=0):
    """純粋な n 次の共鳴 (n ν_plane = 整数) の単項式: ζ̄^(n−1) (方程式 ζ' = λζ + r ζ̄^(n−1))"""
    return (0, n - 1, 0, 0) if plane == 0 else (0, 0, 0, n - 1)


def resonance_normal_form(res, plane, alpha):
    """res (`ampnf.one_turn_series`) の正規形を、共鳴の単項式 alpha だけ残して作る。戻り nf (`ampnf.normal_form` と同じ。nf['keep_index'] = (i, m))"""
    nf = A.normal_form(res, keep=resonance_keep(plane, alpha))
    tab = res['tab']
    nf['keep_index'] = (2 * plane, tab.idx[tuple(alpha) + (0,)])
    return nf


def drive_table(res, orders=(3, 4, 5, 6), planes=(0, 1)):
    """純粋な n 次の共鳴 (n ν_plane = 整数) の駆動項の係数の表 (3〜6 次。5・6 次は D ≥ 5・6 が要る: n 次の共鳴の単項式は次数 n−1、H の n 次)。
    全部の共鳴の単項式を同時に残した 1 回の正規形で作る (残した項どうしの交差は、次数 n1 + n2 − 1 以上の変換に入るので最低次の係数 r には効かない)。
    戻り dict {(plane, n): dict(r (複素)、abs (|r|、ζ̃ の規格化、s = res['s'])、dist (共鳴の分母 |λ^α − λ| = |e^{i(n−1)... } − e^{iμ}|)、nu_dist (n ν − 最寄りの整数))}。
    |r| は線形のチューンでの 1 周の写像の係数 (対称性で禁止された共鳴は 0 に近い)"""
    D = res['D']
    keep = set()
    for pl in planes:
        for n in orders:
            if n - 1 <= D:
                keep |= resonance_keep(pl, pure_alpha(n, pl))
    nf = A.normal_form(res, keep=keep)
    tab = res['tab']
    out = {}
    for pl in planes:
        mu = res['mux'] if pl == 0 else res['muy']
        for n in orders:
            if n - 1 > D:
                continue
            al = pure_alpha(n, pl)
            m = tab.idx[tuple(al) + (0,)]
            rr = nf['Gn'][2 * pl][m]
            nu = mu / (2 * math.pi)
            out[(pl, n)] = dict(r=rr, abs=float(abs(rr)), dist=float(abs(nf['lam'][2 * pl] ** (-(n - 1)) - nf['lam'][2 * pl])), nu_dist=float(n * nu - round(n * nu)))
    return out


def _nu_all(r_nu, jx, jy, jorder=None):
    return np.array([r_nu.nu(0, jx, jy, 0.0, jorder), r_nu.nu(1, jx, jy, 0.0, jorder)])


def _dnu(r_nu, jx, jy, jorder=None):
    """∂ν_i/∂J_j (2×2) の級数の解析微分"""
    out = np.zeros((2, 2))
    for pl in range(2):
        ser = r_nu.nux if pl == 0 else r_nu.nuy
        for (a, b, p), v in ser.items():
            if p or (jorder is not None and a + b > jorder):
                continue
            if a:
                out[pl, 0] += v * a * jx ** (a - 1) * jy ** b
            if b:
                out[pl, 1] += v * b * jx ** a * jy ** (b - 1)
    return out


def island_theory(r_nu, nf, plane, alpha, jother=0.0, solve_for=None, jmax=1e-3, jorder=None):
    """1 本の共鳴の島 (包絡の理論)。r_nu = AmpTuneResult (ν の級数と res)、nf = resonance_normal_form の正規形、plane・alpha = 共鳴の単項式。
    jother: 共鳴の面でない方 (または固定する方) の作用。solve_for: 固定点の作用を解く面 (既定: k が 0 でない面。k の成分が両方 0 でなければ指定する)。
    戻り dict: k (共鳴のベクトル)、p (整数)、Jr = (Jx, Jy) 固定点の作用、r (駆動項の係数)、A (Hamiltonian の共鳴項の振幅 [rad])、M、
    half_J = k の方向の島の半幅 (ΔJx, ΔJy)、Jlo・Jhi = (固定点の作用 ∓ 半幅) の共鳴の面の値、nu_lock = 島の中で固定される ν の幅"""
    res = r_nu.res
    tab = res['tab']
    s = res['s']
    a, b, c, d = alpha
    kx, ky = (a - b - 1, c - d) if plane == 0 else (a - b, c - d - 1)
    i, m = 2 * plane, tab.idx[tuple(alpha) + (0,)]
    r = nf['Gn'][i][m]
    nu0 = np.array(r_nu.nu0)
    p = int(round(kx * nu0[0] + ky * nu0[1]))
    k = np.array([kx, ky], dtype=float)
    if solve_for is None:
        solve_for = 0 if kx != 0 else 1
        if kx != 0 and ky != 0:
            solve_for = plane
    fixed = jother

    def J(j):
        return (j, fixed) if solve_for == 0 else (fixed, j)

    def f(j):
        jx, jy = J(j)
        return k @ _nu_all(r_nu, jx, jy, jorder) - p

    grid = np.concatenate([[0.0], np.logspace(-10, math.log10(jmax), 600)])
    vals = np.array([f(j) for j in grid])
    jr = None
    for q in range(len(grid) - 1):
        if vals[q] == 0.0:
            jr = grid[q]
            break
        if vals[q] * vals[q + 1] < 0:
            lo, hi = grid[q], grid[q + 1]
            flo = vals[q]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                fm = f(mid)
                if fm * flo > 0:
                    lo, flo = mid, fm
                else:
                    hi = mid
            jr = 0.5 * (lo + hi)
            break
    out = dict(k=k, p=p, r=r, alpha=tuple(alpha), plane=plane, ok=False, Jr=None, A=None, half_J=None)
    if jr is None:
        out['reason'] = 'ν(J) の級数が J ≤ %.1e で k·ν = %d に達しない (共鳴は作る振幅の範囲にない)' % (jmax, p)
        out['detune'] = float(k @ nu0 - p)
        return out
    jx, jy = J(jr)
    kick = s * s * abs(r) * (2 * jx / s ** 2) ** ((a + b + (1 if plane == 0 else 0)) / 2.0) * (2 * jy / s ** 2) ** ((c + d + (1 if plane == 1 else 0)) / 2.0)
    # 作用の増分 ΔJ_plane = k_plane · A sinψ の振幅 = kick (|ζ'|² の増分の 1/2 × s²: 追跡の 1 周の写像 ζ' = λζ + r ζ^α)
    kp = abs(kx) if plane == 0 else abs(ky)
    Aamp = kick / kp
    H2 = 2 * math.pi * _dnu(r_nu, jx, jy, jorder)
    H2 = 0.5 * (H2 + H2.T)
    M = float(k @ H2 @ k)
    sh = 2.0 * math.sqrt(Aamp / abs(M))
    half = np.abs(k) * sh
    out.update(ok=True, Jr=(jx, jy), A=Aamp, M=M, half_J=(float(half[0]), float(half[1])), kick=kick,
               nu_lock=float(abs(k @ (_nu_all(r_nu, jx + k[0] * sh, jy + k[1] * sh, jorder) - _nu_all(r_nu, jx - k[0] * sh, jy - k[1] * sh, jorder)))),
               asym_dnu=float(np.linalg.norm(_dnu(r_nu, jx, jy, jorder).diagonal())))
    pl = solve_for
    j0 = (jx, jy)[pl]
    out['Jlo'], out['Jhi'] = j0 - half[pl], j0 + half[pl]
    out['rel_width'] = float(half[pl] / j0) if j0 > 0 else float('inf')
    # 振り子の近似が成り立つ目安: 島の半幅が固定点の作用の 30% 以下 (これを超えると、J 依存の駆動項・ν の高次・他の共鳴の影響が効き、幅は桁の目安)
    out['pendulum_ok'] = bool(out['rel_width'] < 0.3)
    return out


# ---------------------------------------------------------------------------------------------
# 検証用: 実際の追跡の島 (周期 n の固定点 O・X と、X から出る不安定多様体 (分離線) の作用の範囲)。1 自由度 (純粋な x の共鳴: y = py = 0、δ = 0)
# ---------------------------------------------------------------------------------------------
def _plane_map(lat, momentum_eV, plane=0):
    from .track6n import track_particles_fast
    ie = (0, 1) if plane == 0 else (2, 3)

    def f(Z):
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
        c = np.zeros((len(Z), 6))
        c[:, ie[0]] = Z[:, 0]
        c[:, ie[1]] = Z[:, 1]
        o = np.array(track_particles_fast(lat.elements, c, 1, momentum_eV, 0.0, False, False, False))
        return o[:, list(ie)]
    return f


def _fn(f, z, n):
    z = np.atleast_2d(z)
    for _ in range(n):
        z = f(z)
    return z


def _jac(f, z, n, h=1e-7):
    z = np.asarray(z, dtype=float)
    J = np.zeros((2, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = h
        J[:, j] = (_fn(f, z + e, n)[0] - _fn(f, z - e, n)[0]) / (2 * h)
    return J


def action_nf(r_nu, Z, plane=0):
    """面内の座標 Z (m×2) の正規形の作用 (その面だけ、他の面は 0)。`ampnf.nf_action` の T の逆を m 点まとめて"""
    res, nf = r_nu.res, r_nu.nf
    tab, s, sd = res['tab'], res['s'], res['sd']
    Ainv = res['Ainv']
    Z = np.atleast_2d(Z)
    coords = np.zeros((len(Z), 4))
    coords[:, 2 * plane] = Z[:, 0]
    coords[:, 2 * plane + 1] = Z[:, 1]
    zn = (Ainv @ coords.T).T
    xh, ph = zn[:, 2 * plane], zn[:, 2 * plane + 1]
    z = np.zeros((len(Z), 4), dtype=complex)
    z[:, 2 * plane] = (xh - 1j * ph) / s
    z[:, 2 * plane + 1] = (xh + 1j * ph) / s
    w = z.copy()
    mons = tab.mons
    P = np.array([nf['phi'][i] for i in range(4)])
    nz = [np.nonzero(P[i])[0] for i in range(4)]
    for _ in range(80):
        full = np.concatenate([w, np.zeros((len(Z), 1))], axis=1)
        mono = np.ones((len(Z), tab.N), dtype=complex)
        for v in range(5):
            mono *= full[:, v:v + 1] ** mons[None, :, v]
        zz = np.stack([mono @ P[i] for i in range(4)], axis=1)
        w = w - (zz - z)
    return abs(w[:, 2 * plane]) ** 2 * s * s / 2


def track_island(lat, r_nu, n, plane=0, momentum_eV=3e9, jr=None, nturn_sep=40, nth=8, eps=3e-6, nturn_max=3000):
    """追跡で島を測る。周期 n の固定点 (n 周の写像の不動点) を、理論の J_r のまわりの位相の格子から Newton 法で探し、n 周の写像のヤコビ行列の跡で O (|tr| < 2)・X (> 2) に分け、
    X から不安定多様体の方向に eps ずらして nturn_sep 周期 (n × nturn_sep 周) 追跡した点の正規形の作用の範囲 (分離線の J の最大・最小) を返す。
    戻り dict: O・X (各 (座標, 作用, 跡))、sep・Jmin・Jmax (分離線の作用の範囲)。nturn_sep は周期数の下限 (X の不安定な固有値から必要な周期数を見積もり、nturn_max まで増やす)"""
    f = _plane_map(lat, momentum_eV, plane)
    res = r_nu.res
    pts = []
    for q in range(nth * n):
        ph = 2 * math.pi * (q + 0.5) / (nth * n)
        c = A.initial_from_nf(r_nu, jr if plane == 0 else 0.0, jr if plane == 1 else 0.0, ph if plane == 0 else 0.0, ph if plane == 1 else 0.0, 0.0)
        z = np.array([c[2 * plane], c[2 * plane + 1]])
        ok = False
        try:
            for it in range(60):
                g = _fn(f, z, n)[0] - z
                if np.abs(g).max() < 1e-13:
                    ok = True
                    break
                Jm = _jac(f, z, n) - np.eye(2)
                dz = np.linalg.solve(Jm, g)
                if np.abs(dz).max() > 0.1:
                    break
                z = z - dz
        except (np.linalg.LinAlgError, FloatingPointError):
            ok = False
        if not ok:
            continue
        # 最小周期 n
        per = n
        for kk in range(1, n):
            if n % kk == 0 and np.abs(_fn(f, z, kk)[0] - z).max() < 1e-9:
                per = kk
                break
        if per != n:
            continue
        if any(np.abs(z - zz).max() < 1e-8 for zz, _ in pts):
            continue
        # 同じ軌道 (周期 n の他の点) も別の固定点として数える (n 周の写像の不動点 n 個は 1 つの軌道)
        Jm = _jac(f, z, n)
        pts.append((z, float(np.trace(Jm))))
    O = [(z, t) for z, t in pts if abs(t) < 2]
    X = [(z, t) for z, t in pts if abs(t) > 2]
    out = dict(n_O=len(O), n_X=len(X), O=[], X=[])
    for lst, key in ((O, 'O'), (X, 'X')):
        if lst:
            Js = action_nf(r_nu, np.array([z for z, _ in lst]), plane)
            out[key] = [(z, float(j), t) for (z, t), j in zip(lst, Js)]
    if X:
        z0, _, _ = out['X'][0]
        Jm = _jac(f, z0, n)
        w, v = np.linalg.eig(Jm)
        vu = np.real(v[:, int(np.argmax(np.abs(w)))])
        vu = vu / np.linalg.norm(vu)
        # X の不安定な固有値 λu (n 周あたり) が 1 に近い (弱い共鳴の細い島) と、eps から島の大きさまで離れるのに時間がかかる: 周期数を ln(1e-3/eps)/ln λu の 8 倍以上にする
        lam_u = float(np.max(np.abs(w)))
        if lam_u > 1.0 + 1e-12:
            nturn_sep = int(min(max(nturn_sep, math.ceil(8.0 * math.log(1e-3 / eps) / math.log(lam_u))), nturn_max))
        mx = []
        for sgn in (1.0, -1.0):
            z = z0 + sgn * eps * vu
            zs = [z.copy()]
            c = z.reshape(1, 2)
            zmax = 3.0 * float(np.abs(z0).max())
            escaped_flag = False
            try:
                for _ in range(n * nturn_sep):
                    c = f(c)
                    if not np.all(np.isfinite(c)) or np.abs(c).max() > zmax:
                        escaped_flag = True
                        break
                    zs.append(c[0].copy())
            except FloatingPointError:
                escaped_flag = True
            zs = np.array(zs)
            if len(zs) > 4000:       # 点が多いときは等間隔に間引く (T の逆の反復は点数 × 単項式数のメモリを使う)。分離線は連続な曲線なので端の作用は 1e-3 以内で拾える
                zs = zs[np.linspace(0, len(zs) - 1, 4000).astype(int)]
            Js = np.concatenate([action_nf(r_nu, zs[i:i + 500], plane) for i in range(0, len(zs), 500)])
            Js = np.where(np.isfinite(Js) & (Js < 20.0 * np.median(Js[np.isfinite(Js)])), Js, np.nan)       # T の逆の反復が発散した点 (島の外へ出た点) は除く
            mx.append((float(np.nanmin(Js)), float(np.nanmax(Js)), float(np.abs(zs).max()), bool(escaped_flag)))
        out['sep'] = mx          # 各枝: (J の最小, J の最大, |座標| の最大, 脱出して途中で止めたか)
        out['escaped'] = any(m[3] for m in mx)
        out['Jmin'] = min(m[0] for m in mx)
        out['Jmax'] = max(m[1] for m in mx)
    return out
