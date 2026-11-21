# -*- coding: utf-8 -*-
"""v2.3.0 R4: 追跡から共鳴駆動項 (RDT) h_jklmp を取り出す (K64 にない新機能。`resdrive.resonance_terms` の検証用)。

方法 (a) 1 周の写像の 2 次のテイラー展開 (有限振幅の追跡の差分)
  リング 1 周の追跡の写像 M(z) = R z + T(z, z) + ... (z = (x, px, y, py)、RF・放射なし) の 2 次の部分 T は、各要素の kick の 2 次の部分の和そのもの
  (2 つの kick の積は 3 次以上): T = Σ_i R_{i→end} τ_i(R_{0→i} z)。したがって 1 周の kick の生成関数 f_eff (3 次、∇f = −S R⁻¹T) は
  f_eff(z) = Σ_i f_i(R_{0→i} z) で、角度・作用の単項式 z^j z̄^k v^l v̄^m (z = √(2Jx) e^{iφx}、リング先頭の w) の係数が全周の和の h_jklm に等しい (1 次の摂動で厳密)。
  正規化座標 (x/√β, (αx + βpx)/√β) で取り出す (ζ = x_n − i p_n、w = (1, i, 1, i))。
  δ の項 (h11001 など、p=1) は、閉軌道 z_co(δ) でのヤコビ行列の δ 微分から: R0⁻¹ ∂J(z_co; δ)/∂δ = S Q、g2 = ½ zᵀ Q z (全体 = 自然なクロマティックな項
  + 6 極 + 偏向磁石。6 極の分は強さ ×1 と ×0 の差)。δ² の項 h10002 は 2 次の閉軌道の係数 (固定点の式 z* = (1−R)⁻¹ R S∇g δ²、同様に全体)。
方法 (b) 共鳴線の近くの周期点 (`periodic_points`): νx を 3νx = 整数の近くに動かし、3 周写像の周期点の振幅 J を |h30000| と ADTS から予測した値と比べる。
"""
import math
import numpy as np
from . import bendnl as B

S4 = B._S4
CUBIC_MONOS = [(j, k, l, m, 0) for j in range(4) for k in range(4) for l in range(4) for m in range(4) if j + k + l + m == 3]
W_NORM = np.array([1.0, 1.0j, 1.0, 1.0j])          # 正規化座標 (x_n, p_n, y_n, q_n) の w: x_n = Re ζ、p_n = −Im ζ (ζ = √(2J) e^{iφ})


def mono_name(q):
    return 'h%d%d%d%d%d' % tuple(q)


def normalizing_matrix(row):
    """z = A zn (zn = (x/√β, (αx + βpx)/√β, y/√β_y, (α_y y + β_y py)/√β_y))。row: Twiss の行 (betx, alfx, bety, alfy)"""
    bx, ax, by, ay = float(row.betx), float(row.alfx), float(row.bety), float(row.alfy)
    A = np.array([[math.sqrt(bx), 0, 0, 0], [-ax / math.sqrt(bx), 1 / math.sqrt(bx), 0, 0],
                  [0, 0, math.sqrt(by), 0], [0, 0, -ay / math.sqrt(by), 1 / math.sqrt(by)]])
    return A, np.linalg.inv(A)


def one_turn(elements, momentum, delta=0.0):
    """リング 1 周の追跡 f(Z (n×4)) → (n×4)。RF・放射なし、δ は一定"""
    from .track6n import track_particles_fast

    def f(Z):
        Z = np.atleast_2d(np.asarray(Z, dtype=float))
        c = np.zeros((len(Z), 6))
        c[:, :4] = Z
        c[:, 5] = delta
        return np.array(track_particles_fast(elements, c, 1, momentum, 0.0, False, False, False))[:, :4]
    return f


def jacobian(f, z, h=1e-6):
    """z でのヤコビ行列 (中心差分の Richardson)"""
    z = np.asarray(z, dtype=float)

    def cen(hh):
        return np.array([(f(z + hh * e)[0] - f(z - hh * e)[0]) / (2 * hh) for e in np.eye(4)]).T
    return (4.0 * cen(h) - cen(2 * h)) / 3.0


def closed_orbit(elements, momentum, delta=0.0, tol=1e-14, maxit=40):
    """4 次元の閉軌道 (Newton 法)。戻り (z_co, f)"""
    f = one_turn(elements, momentum, delta)
    z = np.zeros(4)
    for _ in range(maxit):
        r = f(z)[0] - z
        if np.abs(r).max() < tol:
            break
        z = z - np.linalg.solve(jacobian(f, z) - np.eye(4), r)
    return z, f


class RingMap:
    """正規化座標の 1 周の写像 (R: 線形、T: 2 次、orbit: 展開の中心 (物理座標)、A: z = A zn)"""
    pass


def ring_map(elements, momentum, row0, delta=0.0, amp=1e-5, hlin=1e-6, orbit=None, tracker=None):
    """リング先頭の Twiss の行 row0 の正規化座標で、閉軌道 (または orbit) のまわりの 1 周の写像の 2 次のテイラー展開 (bendnl.extract_map)。
    amp: 正規化座標 (√m) での振幅 (Richardson の振幅 amp と 2·amp。誤差 ~ amp⁴、大きなリング A (603 m) で amp = 1e-4: 8e-4、3e-5: 6e-6、1e-5: 4e-7)。J = amp²/2 の程度の振幅の追跡。
    tracker: 1 周の追跡 f(Z (n×4)) → (n×4) (物理座標、δ=0) を差し替える (gs との照合用。orbit も与える)"""
    A, Ai = normalizing_matrix(row0)
    f = tracker if tracker is not None else one_turn(elements, momentum, delta)
    if orbit is None:
        orbit, _ = closed_orbit(elements, momentum, delta)
    orbit = np.asarray(orbit, dtype=float)

    def fn(Zn):
        Z = (A @ np.atleast_2d(Zn).T).T + orbit
        return (Ai @ (f(Z) - orbit).T).T
    m = B.extract_map(fn, amp=amp, hlin=hlin)
    rm = RingMap()
    rm.R, rm.T, rm.orbit, rm.A, rm.delta, rm.amp = m.R, m.T, orbit, A, delta, amp
    return rm


def cubic_generator(rm):
    """1 周の kick の生成関数 f_eff (正規化座標の 3 次式 {(ex, ep, ey, eq): 係数}): ∇f = −S R⁻¹ T、f = (1/3) z·∇f"""
    Tk = np.einsum('ij,jab->iab', np.linalg.inv(rm.R), rm.T)
    return B.cubic_from_kick(Tk)


def rdt_from_map(rm, monos=CUBIC_MONOS):
    """1 周の写像から取り出した h_jklm (p=0 の 3 次の項、全部の和 = 6 極 + 偏向磁石 + その他の 2 次の kick)。戻り {名前: 複素数}"""
    amp = B.harmonic_amplitudes(cubic_generator(rm), W_NORM, monos)[:, 0]
    return {mono_name(q): complex(a) for q, a in zip(monos, amp)}


def quad_to_h(Q, monos):
    """2 次の生成関数 g2 = ½ zᵀ Q z (正規化座標) の h_(j,k,l,m,1) (monos の各要素は p=1)。戻り {名前: 複素数}"""
    f2 = {}
    for a in range(4):
        for b in range(4):
            e = [0] * 4
            e[a] += 1
            e[b] += 1
            f2[tuple(e)] = f2.get(tuple(e), 0.0) + 0.5 * Q[a, b]
    amp = B.harmonic_amplitudes(f2, W_NORM, [tuple(q[:4]) + (0,) for q in monos])[:, 0]
    return {mono_name(q): complex(a) for q, a in zip(monos, amp)}


DELTA1_MONOS = [(1, 1, 0, 0, 1), (0, 0, 1, 1, 1), (2, 0, 0, 0, 1), (0, 2, 0, 0, 1), (0, 0, 2, 0, 1), (0, 0, 0, 2, 1), (1, 0, 1, 0, 1), (1, 0, 0, 1, 1)]


def delta_quadratic_total(elements, momentum, row0, d=1e-4):
    """δ の 1 次の 2 次形式 (閉軌道でのヤコビ行列の δ 微分): R0⁻¹ ∂J(z_co(δ); δ)/∂δ = S Q (正規化座標、中心差分の Richardson、δ = d と 2d)、g2 = ½ zᵀ Q z。
    **全体** = 四極・偏向磁石の自然なクロマティックな項 + 6 極 + 偏向磁石の曲率・端 (分散による軌道のずれ)。6 極の分は強さを ×1 と ×0 にした結果の差
    (T は kick の和で、この 1 次の量は 6 極の強さに厳密に比例)。戻り (Q (4×4)、dict h (名前 → 複素数))"""
    A, Ai = normalizing_matrix(row0)
    f0 = one_turn(elements, momentum, 0.0)
    R0 = Ai @ jacobian(f0, np.zeros(4)) @ A

    def dj(dd):
        zp, fp = closed_orbit(elements, momentum, dd)
        zm, fm = closed_orbit(elements, momentum, -dd)
        return Ai @ (jacobian(fp, zp) - jacobian(fm, zm)) @ A / (2 * dd)
    X = (4.0 * dj(d) - dj(2 * d)) / 3.0
    Q = -S4 @ np.linalg.solve(R0, X)
    return Q, quad_to_h(Q, DELTA1_MONOS)


def delta2_orbit_h10002(elements, momentum, row0, d=1e-4):
    """2 次の閉軌道の係数 c2 = ½ ∂²z_co/∂δ² (正規化座標) から、固定点の式 z* = (1−R)⁻¹ R S ∇g δ² の線形項 g = ∇g·z の h10002 (複素数) を返す。
    **全体 (自然な 2 次の分散 + 6 極 + 偏向磁石)**。6 極の分は 6 極の強さを ×1 と ×0 にした結果の差 (呼び出し側)。戻り (h10002, c2)"""
    A, Ai = normalizing_matrix(row0)
    f0 = one_turn(elements, momentum, 0.0)
    R = Ai @ jacobian(f0, np.zeros(4)) @ A
    zs = []
    for dd in (d, 2 * d):
        zp, _ = closed_orbit(elements, momentum, dd)
        zm, _ = closed_orbit(elements, momentum, -dd)
        zs.append(Ai @ (0.5 * (zp + zm)) / (dd * dd))
    c2 = (4 * zs[0] - zs[1]) / 3.0
    gradg = -S4 @ np.linalg.solve(R, (np.eye(4) - R) @ c2)       # S∇g = R⁻¹(1−R) c2
    amp = B.harmonic_amplitudes({(1, 0, 0, 0): gradg[0], (0, 1, 0, 0): gradg[1], (0, 0, 1, 0): gradg[2], (0, 0, 0, 1): gradg[3]},
                                W_NORM, [(1, 0, 0, 0)])[0, 0]
    return complex(amp), c2


def tracked_rdt_table(elements, momentum, row0, with_delta=True, amp=1e-5, d=1e-4):
    """追跡から取り出した h の表 (全体 = 6 極 + 偏向磁石 + 自然な項): {'cubic': {名前: 複素数 (p=0)}、'delta1': {名前: 複素数 (p=1、自然なクロマティックな項を含む)}、
    'h10002': 複素数 (p=2、自然な 2 次の分散の分を含む)}。6 極だけの分は強さ ×1 と ×0 の差"""
    rm = ring_map(elements, momentum, row0, amp=amp)
    out = {'cubic': rdt_from_map(rm), 'R': rm.R}
    if with_delta:
        Q, h1 = delta_quadratic_total(elements, momentum, row0, d=d)
        out['delta1'] = h1
        out['h10002'], _ = delta2_orbit_h10002(elements, momentum, row0, d=d)
    return out


# ---------------------------------------------------------------------------------------------
# 方法 (b): 共鳴線の近くの周期点
# ---------------------------------------------------------------------------------------------
def periodic_points(elements, momentum, row0, nper, jguess, nphase=24, tol=1e-12):
    """y = 0 の x の運動で、nper 周の写像 T^nper の固定点 (周期点) を Newton 法で探す (位相 nphase 通りの初期値)。
    戻り: list[(J, φ)] (J = アクション (x = √(2Jβ) cos φ)、重複は除く)"""
    from scipy.optimize import fsolve
    f = one_turn(elements, momentum, 0.0)
    bx, ax = float(row0.betx), float(row0.alfx)

    def T(v):
        z = np.array([[v[0], v[1], 0.0, 0.0]])
        for _ in range(nper):
            z = f(z)
        return np.array([z[0, 0], z[0, 1]]) - v
    pts = []
    for ph in np.arange(nphase) * 2 * math.pi / nphase:
        x0 = np.array([math.sqrt(2 * jguess * bx) * math.cos(ph), -math.sqrt(2 * jguess / bx) * (math.sin(ph) + ax * math.cos(ph))])
        try:
            sol, info, ier, msg = fsolve(T, x0, full_output=True, xtol=1e-12)
        except Exception:
            continue
        if np.abs(T(sol)).max() > tol:
            continue
        x, px = sol
        J = (x * x + (ax * x + bx * px) ** 2) / (2 * bx)
        if J < 1e-3 * jguess:                 # 原点 (自明な固定点)
            continue
        phi = math.atan2(-(ax * x + bx * px), x)
        if not any(abs(J - p[0]) < 1e-6 * J and abs(math.remainder(phi - p[1], 2 * math.pi)) < 1e-3 for p in pts):
            pts.append((J, phi))
    return pts


def predict_periodic_J(eps, a, hmag, nint=3):
    """3 周写像の周期点の振幅 J の予測: 1 周の生成関数 G = 2π ν J + π a J² + 2|h|(2J)^{3/2} cos(3φ+θ) (ν = 共鳴の整数/3 + ε) の固定点の条件
    2π(ε + aJ) = ∓ 6|h| √(2J) (a = ∂ν/∂J)。戻り: 昇順の J のリスト (最大 2 つ)"""
    from scipy.optimize import brentq
    out = []
    for sig in (+1, -1):
        g = lambda u: 2 * math.pi * (eps + a * u * u / 2) + 6 * abs(hmag) * u * sig
        try:
            if g(1e-9) * g(1.0) < 0:
                u = brentq(g, 1e-9, 1.0)
                out.append(u * u / 2)
        except ValueError:
            pass
    return sorted(out)


# ---------------------------------------------------------------------------------------------
# v2.3.0 R5: 1 周の写像の 3 次のテイラー展開 → 4 次の生成関数 f4 (Dragt–Finn。h4 の数値の取り出し)
# ---------------------------------------------------------------------------------------------
_MON4VAR3 = [e for e in [(a, b, c, d) for a in range(4) for b in range(4) for c in range(4) for d in range(4)] if sum(e) == 3]


def _evalmon(Z, exps):
    return np.prod(Z[:, None, :] ** np.array(exps)[None, :, :], axis=2)


def ring_map3(elements, momentum, row0, delta=0.0, amp=1e-5, amp3=2e-4, hlin=1e-6, orbit=None, tracker=None, ndir=48, seed=1):
    """リング先頭の Twiss の行 row0 の正規化座標で、閉軌道のまわりの 1 周の写像の 3 次までのテイラー展開 (R、T (2 次)、T3 (3 次))。
    R・T は `ring_map` (amp)。T3[i] は 3 次の同次式の係数 {(ea, eb, ec, ed): 係数} (i = 0..3 の出力成分、正規化座標)。取り出し: 方向 v (ndir 個の乱数の単位ベクトル) で
    奇関数 o(a) = ½(f(a v) − f(−a v)) (= a R v + a³ T3(v) + a⁵ T5(v) …) の差 u1 = o(2a) − 2 o(a) = 6a³T3 + 30a⁵T5、u2 = o(3a) − 3 o(a) = 24a³T3 + 240a⁵T5 から
    T3(v) = (8 u1 − u2)/(24 a³) (誤差 ~ a⁴ T7、線形項は差で厳密に消える。丸め誤差 ~ 1e-16/a²)。a = amp3 (正規化座標 √m)。T3 は 3 次の同次式に最小二乗フィット"""
    rm = ring_map(elements, momentum, row0, delta=delta, amp=amp, hlin=hlin, orbit=orbit, tracker=tracker)
    A, Ai = rm.A, np.linalg.inv(rm.A)
    f = tracker if tracker is not None else one_turn(elements, momentum, delta)
    orbit = rm.orbit

    def fn(Zn):
        Z = (A @ np.atleast_2d(Zn).T).T + orbit
        return (Ai @ (f(Z) - orbit).T).T
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(ndir, 4))
    V /= np.linalg.norm(V, axis=1)[:, None]

    def odd(a):
        return 0.5 * (fn(a * V) - fn(-a * V))
    o1, o2, o3 = odd(amp3), odd(2 * amp3), odd(3 * amp3)
    u1, u2 = o2 - 2 * o1, o3 - 3 * o1
    T3v = (8 * u1 - u2) / (24 * amp3 ** 3)             # (ndir, 4)
    Mx = _evalmon(V, _MON4VAR3)                         # (ndir, 20)
    coef, res, rk, sv = np.linalg.lstsq(Mx, T3v, rcond=None)
    fit_res = float(np.abs(Mx @ coef - T3v).max() / max(np.abs(T3v).max(), 1e-300))
    rm.T3 = [{e: float(coef[m, i]) for m, e in enumerate(_MON4VAR3)} for i in range(4)]
    rm.amp3, rm.fit_res3, rm.cond3 = amp3, fit_res, float(sv[0] / sv[-1])
    return rm


def _padd(p, q, s=1.0):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0.0) + s * v
    return r


def _pmulp(p, q):
    r = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            r[k] = r.get(k, 0.0) + c1 * c2
    return r


def _pdiff(p, j):
    r = {}
    for k, c in p.items():
        if k[j]:
            kk = list(k)
            kk[j] -= 1
            r[tuple(kk)] = r.get(tuple(kk), 0.0) + c * k[j]
    return r


def _quad_poly(Tk_i):
    """τ_i(z) = zᵀ Tk[i] z → {(e..): 係数}"""
    p = {}
    for a in range(4):
        for b in range(4):
            e = [0] * 4
            e[a] += 1
            e[b] += 1
            p[tuple(e)] = p.get(tuple(e), 0.0) + Tk_i[a, b]
    return p


def f4_from_map3(rm):
    """3 次のテイラー展開から 4 次の生成関数 h4 (正規化座標の 4 次式 {(ex, ep, ey, eq): 係数}) を作る。kick の 3 次の部分 τ3 = R⁻¹ T3 から、2 次の kick u = R⁻¹T (= S∇h3) の流れの
    2 次 ½(u·∇)u を引き (Dragt–Finn: 写像 = e^{−:h3 + h4:}、3 次の部分 = S∇h4 + ½(u·∇)u)、∇h4 = −S v (v = S∇h4)、h4 = (1/4) z·∇h4。
    戻り (h4, info): info['symp'] = ∇h4 の可積分条件 (∂_k g_j = ∂_j g_k) の相対残差 (写像が 4 次までシンプレクティックなら ~0)"""
    Ri = np.linalg.inv(rm.R)
    Tk = np.einsum('ij,jab->iab', Ri, rm.T)
    u = [_quad_poly(Tk[i]) for i in range(4)]
    t3 = [dict() for _ in range(4)]
    for j in range(4):
        for e, c in rm.T3[j].items():
            for i in range(4):
                t3[i][e] = t3[i].get(e, 0.0) + Ri[i, j] * c
    v = []
    for i in range(4):
        uu = {}
        for j in range(4):
            uu = _padd(uu, _pmulp(u[j], _pdiff(u[i], j)))
        v.append(_padd(t3[i], uu, -0.5))
    g = [dict() for _ in range(4)]
    for j in range(4):
        for i in range(4):
            if S4[j, i] != 0.0:
                g[j] = _padd(g[j], v[i], -S4[j, i])
    h4 = {}
    for j in range(4):
        for e, c in g[j].items():
            ee = list(e)
            ee[j] += 1
            h4[tuple(ee)] = h4.get(tuple(ee), 0.0) + 0.25 * c
    num = den = 0.0
    for j in range(4):
        for k in range(j + 1, 4):
            d = _padd(_pdiff(g[j], k), _pdiff(g[k], j), -1.0)
            num = max(num, max([abs(c) for c in d.values()] + [0.0]))
            den = max(den, max([abs(c) for c in _pdiff(g[j], k).values()] + [0.0]))
    return h4, {'symp': num / max(den, 1e-300), 'v': v}


MONOS4V = [(j, k, l, 4 - j - k - l) for j in range(5) for k in range(5 - j) for l in range(5 - j - k)]


def rdt4_from_map(rm, monos=MONOS4V):
    """追跡の 1 周の写像の 3 次のテイラー展開 (ring_map3) から取り出した h4_jklm (p=0。全体 = 8 極 + 運動学 + フリンジ + 3 次の 2 次 + 四極・偏向磁石の本体の 4 次の項)。
    戻り ({名前: 複素数}, info)"""
    h4, info = f4_from_map3(rm)
    amp = B.harmonic_amplitudes(h4, W_NORM, [tuple(q) + (0,) for q in monos])[:, 0]
    return {'h%d%d%d%d0' % tuple(q): complex(a) for q, a in zip(monos, amp)}, info


_ELEM_H4_CACHE = {}


def element_h4(el, momentum=3e9, amp3=2e-4, orbit=None, delta=0.0):
    """1 要素の写像 (入口の物理座標 (x, px, y, py)、閉軌道 0) の 4 次の生成関数 h4 (物理座標の 4 次式 {(ex, epx, ey, epy): 係数})。
    N = R_el⁻¹ M_el (要素自身の線形写像を除いた相互作用表示の写像) の 3 次のテイラー展開から `f4_from_map3` で取り出す。
    要素の内部の 3 次どうしの 2 次の項 (Magnus の第 2 項) と 4 次のハミルトニアンの 1 次 (運動学・端のフリンジなど) の両方を含む。
    入口の Twiss に依らない (物理座標) ので、素子の種類ごとに 1 回 (キャッシュ)。戻り (h4, info)
    v2.3.0 R6: orbit (入口の閉軌道 (x, px, y, py)。K0≠0 の偏向磁石 = 挿入光源は内部の軌道が傾いているので必須)・delta を渡すと、その軌道・δ のまわりの
    展開 (出口の軌道は要素の写像で決まる。写像の出力から出口の軌道を引く)"""
    from types import SimpleNamespace
    orb = np.zeros(4) if orbit is None else np.asarray(orbit, dtype=float)
    key = (type(el).__name__, tuple((k, repr(v)) for k, v in sorted(vars(el).items()) if k != 'name'), momentum, amp3, B._mode_key(), float(delta), tuple(np.round(orb, 13)))
    if key in _ELEM_H4_CACHE:
        return _ELEM_H4_CACHE[key]
    row = SimpleNamespace(betx=1.0, alfx=0.0, mux=0.0, bety=1.0, alfy=0.0, muy=0.0)
    f_el = one_turn([el], momentum, delta)
    z_out = f_el(orb[None, :])[0]
    rm = ring_map3([el], momentum, row, amp3=amp3, orbit=orb, delta=delta, tracker=(lambda Z: f_el(Z) - z_out + orb))
    h4, info = f4_from_map3(rm)
    info = {'symp': info['symp'], 'fit': rm.fit_res3}
    _ELEM_H4_CACHE[key] = (h4, info)
    return h4, info


# ---------------------------------------------------------------------------------------------
# v2.3.0 R6: 追跡の光学 (δ のまわり): CALC ではなく、追跡の写像のヤコビ行列から閉軌道・Twiss を作る
# ---------------------------------------------------------------------------------------------
def tracked_rows(elements, momentum, delta=0.0, base_rows=None, h=1e-5):
    """δ の閉軌道を追跡 (Newton 法) で求め、各素子の追跡の写像 (`bendnl._tracker_map`) のヤコビ行列 (閉軌道の各点で、中心差分の Richardson) を掛けて
    各素子の入口の Twiss (β、α、位相進み μ) と閉軌道を作る (結合なしの 2×2 ブロック)。戻り: 行のリスト (長さ N+1、`resdrive` が使う属性を持つ SimpleNamespace)。
    CALC の δ の光学 (`lat.twiss_table(delta=)`) との違い: 追跡の閉軌道は 6 極の 2 次の kick (2 次の分散) を含み、ヤコビ行列は追跡そのもの。
    ASTRID2 (δ = 0.002) で CALC と追跡のチューンは 5e-5〜8e-5 違う (CALC は 2 次の分散・6 極の 2 次のクロマティシティの分を含まない)。
    base_rows: 位相進みの整数部・s・分散 (etx 等) を借りる Twiss の行 (CALC の δ=0 または δ)。整数部は base_rows の mux・muy に最も近いものを選ぶ"""
    from types import SimpleNamespace
    z0, f_ring = closed_orbit(elements, momentum, delta)
    n = len(elements)
    zs = [z0]
    Js = []
    for el in elements:
        f = B._tracker_map([el], momentum, delta)
        z = zs[-1]
        Js.append(B.jac_at(f, z, h=1e-4))
        zs.append(f(z[None, :])[0])
    M = np.eye(4)
    for J in Js:
        M = J @ M

    def block(M, a):
        return M[a:a + 2, a:a + 2]
    tw = []
    for a in (0, 2):
        B2 = block(M, a)
        cm = 0.5 * (B2[0, 0] + B2[1, 1])
        sm = math.copysign(math.sqrt(max(1.0 - cm * cm, 1e-300)), B2[0, 1])
        tw.append((B2[0, 1] / sm, (B2[0, 0] - B2[1, 1]) / (2.0 * sm), math.atan2(sm, cm)))
    rows = []
    state = [(tw[0][0], tw[0][1], 0.0), (tw[1][0], tw[1][1], 0.0)]
    for i in range(n + 1):
        r = SimpleNamespace()
        r.betx, r.alfx, r.mux = state[0]
        r.bety, r.alfy, r.muy = state[1]
        r.x_cod, r.px_cod, r.y_cod, r.py_cod = [float(v) for v in zs[i]]
        r.r1 = r.r2 = r.r3 = r.r4 = 0.0
        r.s = 0.0
        r.etx = r.etpx = r.ety = r.etpy = 0.0
        if base_rows is not None:
            b = base_rows[i]
            r.s = float(b.s)
            r.etx, r.etpx = float(b.etx), float(b.etpx)
            r.ety, r.etpy = float(getattr(b, 'ety', 0.0) or 0.0), float(getattr(b, 'etpy', 0.0) or 0.0)
        rows.append(r)
        if i == n:
            break
        J = Js[i]
        new = []
        for k, a in enumerate((0, 2)):
            b0, a0, mu0 = state[k]
            m11, m12, m21, m22 = J[a, a], J[a, a + 1], J[a + 1, a], J[a + 1, a + 1]
            g0 = (1.0 + a0 * a0) / b0
            b1 = m11 * m11 * b0 - 2 * m11 * m12 * a0 + m12 * m12 * g0
            a1 = -m11 * m21 * b0 + (m11 * m22 + m12 * m21) * a0 - m12 * m22 * g0
            dmu = math.atan2(m12, m11 * b0 - m12 * a0)
            if dmu < -1e-9:
                dmu += 2 * math.pi
            new.append((b1, a1, mu0 + dmu))
        state = new
    # 位相進みの整数部: 1 周の合計が tw の μ (mod 2π) と一致するよう、周ごとの増分の合計を整数部で補正する (増分は連続なので合計を base_rows に合わせる)
    for k, name in ((0, 'mux'), (1, 'muy')):
        tot = getattr(rows[-1], name)
        frac = tw[k][2]
        # 合計は frac + 2π m。m は、増分の総和から決まる
        m = round((tot - frac) / (2 * math.pi))
        sc = (frac + 2 * math.pi * m) / tot if tot != 0 else 1.0
        if abs(sc - 1.0) > 1e-6:
            raise RuntimeError('tracked_rows: 位相進みの合計 %.6f と 1 周の写像のチューン %.6f が合わない' % (tot / (2 * math.pi), (frac + 2 * math.pi * m) / (2 * math.pi)))
    return rows
