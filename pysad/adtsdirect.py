# -*- coding: utf-8 -*-
"""v2.3.0 R2: 実際の強さの追跡 + NAFF から ADTS・チューンの δ 依存を直接求める (多項式フィット)。K64 にない機能 (解析値 resdrive.adts との比較用)。

Q2 の比較 (6 極の強さを ×(+1,-1,0) にした偶数次の差) は交差項を消す間接的な方法だった。ここは**強さを変えず**に、
アクション (Jx, Jy) の格子 (各 ±位相の対) を追跡し、各粒子の NAFF のチューンを J の多項式でフィットして ∂ν/∂J を取る。

J の定義 (解析側・追跡側とも): x = √(2Jx βx) cos φx、px = −√(2Jx/βx)(sin φx + αx cos φx) (y も同様)。単位 m·rad。
ADTS = ∂νx/∂Jx、∂νx/∂Jy、∂νy/∂Jx、∂νy/∂Jy [1/(m rad)] (νx = ν0 + Σ c_ab Jx^a Jy^b の c_10 など)。

  grid_actions(jmax, nstep)          Jx, Jy ∈ {0, 1, ..., nstep}·jmax/nstep の格子 (x だけ・y だけ・両方を含む)
  track_record(...)                  追跡 (pysad の k64、または tracker 関数)。座標列 (ターン × 粒子 × 6)
  record_tunes(rec, twiss0, ...)     各粒子の NAFF のチューン (符号つき端数)。窓・周数・初期位相・探索幅
  fit_tune_surface(jx, jy, nu, ...)  ν(Jx, Jy) の多項式 (全次数 ≤ order) の最小二乗。係数 c_ab とフィットの残差
  adts_from_record(...)              上をまとめ、ADTS 4 個と ν0 を返す
  fit_delta_curve(delta, nu, order)  ν(δ) の多項式フィットで ∂ν/∂δ、∂²ν/∂δ²

±位相の対: (φx, φy) と (φx+π, φy+π) の 2 粒子のチューンを平均する (座標の符号反転 = 振幅の奇数次 (J^{3/2} など) の項が消える)。
"""
import math
import numpy as np

TWO_PI = 2.0 * math.pi


def fr(a):
    """端数 (−0.5〜0.5)"""
    return (np.asarray(a) + 0.5) % 1.0 - 0.5


def grid_actions(jmax, nstep=5):
    """Jx, Jy ∈ {0..nstep}·jmax/nstep の (nstep+1)² 点 (Jx, Jy の配列)。(0, 0) は含めない"""
    a = np.arange(nstep + 1) * jmax / nstep
    jx, jy = np.meshgrid(a, a, indexing='ij')
    jx, jy = jx.ravel(), jy.ravel()
    keep = (jx > 0) | (jy > 0)
    return jx[keep], jy[keep]


def initial_coords(twiss0, jx, jy, phix=0.0, phiy=0.0, delta=0.0, orbit=None):
    """アクション・位相から (x, px, y, py, z, δ)。twiss0 = ((αx, βx), (αy, βy))。orbit = (x, px, y, py) の閉軌道 (δ≠0 の追跡)"""
    (ax, bx), (ay, by) = twiss0
    jx = np.atleast_1d(np.asarray(jx, dtype=float))
    jy = np.atleast_1d(np.asarray(jy, dtype=float))
    ph = np.atleast_1d(np.asarray(phix, dtype=float)) * np.ones_like(jx)
    pv = np.atleast_1d(np.asarray(phiy, dtype=float)) * np.ones_like(jx)
    c = np.zeros((len(jx), 6))
    c[:, 0] = np.sqrt(2 * jx * bx) * np.cos(ph)
    c[:, 1] = -np.sqrt(2 * jx / bx) * (np.sin(ph) + ax * np.cos(ph))
    c[:, 2] = np.sqrt(2 * jy * by) * np.cos(pv)
    c[:, 3] = -np.sqrt(2 * jy / by) * (np.sin(pv) + ay * np.cos(pv))
    c[:, 5] = delta
    if orbit is not None:
        c[:, :4] += np.asarray(orbit, dtype=float)
    return c


def track_record(elements, momentum_eV, coords, nturn=1024, tracker=None):
    """coords (粒子数 × 6) を nturn 周追跡 (RF 空洞は止める rfsw=False、放射なし) し、各周の始点の座標 (nturn, n, 6) を返す。
    tracker(c) は 1 周の写像 (省略すると pysad の track6n = k64 の追跡)"""
    from .track6n import track_particles_fast
    tr = tracker or (lambda cc: track_particles_fast(elements, cc, 1, momentum_eV, 0.0, False, False, False))
    c = np.array(coords, dtype=float)
    rec = np.zeros((nturn, len(c), 6))
    for t in range(nturn):
        rec[t] = c
        c = tr(c)
    return rec


def record_tunes(rec, twiss0, window=1, guess=None, halfwidth=0.02, nturn=None, skip=0, orbit=None):
    """座標列 rec (ターン × 粒子 × 6) の各粒子の NAFF のチューン (νx, νy) (符号つき端数、粒子数 × 2)。
    窓の次数 window (Laskar の Hann^p)、周数 nturn (先頭から)、skip (先頭の捨てる周数 = 初期位相の位置をずらす)。
    guess = (fx0, fy0) (線形のチューンの端数) を与えると fx0 ± halfwidth の範囲のピークを探す (強制振動のピークを避ける)。
    orbit は閉軌道の座標 (差し引く。省略すると平均を引く)"""
    from . import naff
    (ax, bx), (ay, by) = twiss0
    n = rec.shape[1]
    t1 = rec.shape[0] if nturn is None else skip + nturn
    out = np.zeros((n, 2))
    for i in range(n):
        for p, (al, be) in enumerate(((ax, bx), (ay, by))):
            x = rec[skip:t1, i, 2 * p]
            px = rec[skip:t1, i, 2 * p + 1]
            if orbit is not None:
                x = x - orbit[2 * p]
                px = px - orbit[2 * p + 1]
            xh, ph = naff.twiss_normalize(x, px, al, be)
            z = (xh - xh.mean()) - 1j * (ph - ph.mean())
            if guess is None:
                f, _a = naff.naff(z, 1, window=window)
            else:
                g = guess[p]
                f, _a = naff.naff(z, 1, window=window, fmin=g - halfwidth, fmax=g + halfwidth)
            out[i, p] = f[0]
    return out


def _terms(order, need_x, need_y):
    """全次数 1..order の単項式 (a, b): Jx^a Jy^b"""
    return [(a, k - a) for k in range(1, order + 1) for a in range(k, -1, -1)]


def fit_tune_surface(jx, jy, nu, order=3, jscale=None, weights=None):
    """ν(Jx, Jy) = ν0 + Σ_{1≤a+b≤order} c_ab Jx^a Jy^b の最小二乗。戻り dict: nu0、c (キー (a,b) → c_ab、単位 1/(m rad)^(a+b))、rms (フィットの残差の rms)、
    maxres、n (点数)、terms。jscale (J の規格化、既定 max J) で条件数を保つ"""
    jx = np.asarray(jx, dtype=float)
    jy = np.asarray(jy, dtype=float)
    nu = np.asarray(nu, dtype=float)
    s = jscale if jscale is not None else max(jx.max(), jy.max())
    terms = _terms(order, True, True)
    A = np.column_stack([np.ones_like(jx)] + [(jx / s) ** a * (jy / s) ** b for a, b in terms])
    w = np.ones_like(jx) if weights is None else np.asarray(weights, dtype=float)
    sol, *_ = np.linalg.lstsq(A * w[:, None], nu * w, rcond=None)
    res = nu - A @ sol
    c = {t: sol[1 + i] / s ** (t[0] + t[1]) for i, t in enumerate(terms)}
    return {'nu0': float(sol[0]), 'c': c, 'rms': float(np.sqrt(np.mean(res ** 2))), 'maxres': float(np.abs(res).max()), 'n': len(nu),
            'terms': terms, 'res': res}


def pair_average(nu, npair_axis=0):
    """±位相の対 (前半 = (φ, φ)、後半 = (φ+π, φ+π)) のチューンの平均 (粒子数が偶数、前半と後半が対)"""
    n = len(nu) // 2
    return 0.5 * (nu[:n] + nu[n:2 * n])


def adts_from_tunes(jx, jy, nu, order=3, nu_lin=None):
    """nu (粒子数 × 2、符号つき端数。jx, jy はその粒子のアクション) から ADTS を求める。
    νx のフィットは Jx>0 の点、νy は Jy>0 の点 (その面の信号がある粒子) だけを使う。nu_lin (線形のチューン端数 (fx, fy)) があれば、
    フィットの前に ν から fr(ν − nu_lin) + nu_lin にして整数・折り返しを揃える。戻り dict (ν0、xx = ∂νx/∂Jx、xy = ∂νx/∂Jy、yx = ∂νy/∂Jx、yy = ∂νy/∂Jy、rms、フィット)"""
    jx = np.asarray(jx, dtype=float)
    jy = np.asarray(jy, dtype=float)
    nu = np.array(nu, dtype=float)
    if nu_lin is not None:
        nu[:, 0] = fr(nu[:, 0] - nu_lin[0])
        nu[:, 1] = fr(nu[:, 1] - nu_lin[1])
    mx, my = jx > 0, jy > 0
    s = max(jx.max(), jy.max())
    fx = fit_tune_surface(jx[mx], jy[mx], nu[mx, 0], order, s)
    fy = fit_tune_surface(jx[my], jy[my], nu[my, 1], order, s)
    off = (0.0, 0.0) if nu_lin is None else tuple(nu_lin)
    return {'nu0': (fx['nu0'] + off[0], fy['nu0'] + off[1]), 'xx': fx['c'][(1, 0)], 'xy': fx['c'][(0, 1)], 'yx': fy['c'][(1, 0)], 'yy': fy['c'][(0, 1)],
            'rms': (fx['rms'], fy['rms']), 'fitx': fx, 'fity': fy, 'order': order, 'jmax': float(s)}


def fit_delta_curve(delta, nu, order=3):
    """ν(δ) = Σ_{n≤order} c_n δ^n の最小二乗。戻り dict: nu0、d1 = ∂ν/∂δ、d2 = ∂²ν/∂δ²、d3 (n! 倍)、rms"""
    delta = np.asarray(delta, dtype=float)
    nu = np.asarray(nu, dtype=float)
    s = np.abs(delta).max()
    co = np.polynomial.polynomial.polyfit(delta / s, nu, order)
    res = nu - np.polynomial.polynomial.polyval(delta / s, co)
    out = {'nu0': float(co[0]), 'rms': float(np.sqrt(np.mean(res ** 2)))}
    for n in range(1, order + 1):
        out['d%d' % n] = float(co[n] / s ** n * math.factorial(n))
    return out


# ---------------------------------------------------------------------------------------------
# パラキシャルなハミルトニアンだけの追跡 (resdrive の解析式と同じ模型): 線形の転送行列 + 薄い 6 極・8 極の kick
# ---------------------------------------------------------------------------------------------
def paraxial_ops(lat, nslice=40, bend_cubic=False, bend_edge=False):
    """lat の各素子を、線形の 4×4 転送行列 (δ=0、素子の transfer_matrix の x, px, y, py の部分) と、薄い 6 極・8 極の kick
    に分けた操作列 [(M4, kind, k)]。SEXT・OCT・薄い MULT (K2L・K3L) が kick。厚い SEXT・OCT は nslice 個の薄レンズ (中点則: ドリフト h/2, kick, h/2 …)。
    kick は pysad の追跡と同じ向き: 6 極 Δpx = −(K2L/2)(x²−y²)、Δpy = +K2L·x·y。8 極 Δpx = −(K3L/6)(x³−3xy²)、Δpy = −(K3L/6)(y³−3x²y)。
    厚い MULT (K1≠0 で K2・K3 つき)、スキュー・回転・DX・DY つきの多極は未対応 (例外)。他の素子 (BEND・QUAD・DRIFT・...) は線形の行列だけ"""
    from .elements import Drift
    ops = []
    M = np.eye(4)

    def lin(el):
        return np.asarray(el.transfer_matrix(0.0), dtype=float)[:4, :4]

    for el in lat.elements:
        tn = type(el).__name__
        L = float(getattr(el, 'L', 0.0) or 0.0)
        k2l = k3l = 0.0
        if tn == 'Bend' and float(getattr(el, 'K0', 0.0) or 0.0) != 0.0:
            raise NotImplementedError('paraxial_ops: K0 のある BEND (挿入光源): %s' % el.name)
        if bend_cubic and tn == 'Bend' and L > 0 and float(el.ANGLE) != 0.0:
            # 曲率の 3 次の項 H3b = (h/2) x (px² + py²) の流れを、偏向磁石を nslice 分割して Strang 分割で入れる (端部 E1=E2=0 のセクター型だけ)
            import dataclasses
            if el.E1 or el.E2:
                raise NotImplementedError('paraxial_ops(bend_cubic): 端 E1・E2≠0 の偏向磁石は未対応: %s' % el.name)
            h = float(el.ANGLE) / L
            ds = L / nslice
            if bend_edge:
                ops.append((M, 8, (h, 0.0)))      # 入口の端 (E1=0 の硬い端): H = +(h/2) px y²。ここまでの行列 M を先に掛けてから
                M = np.eye(4)
            kw0_ = {'K0': float(el.K0) * ds / 2 / L} if float(getattr(el, 'K0', 0.0) or 0.0) != 0.0 else {}      # v2 (B237): K0 も角度と同じ比で縮める
            Mh = np.asarray(dataclasses.replace(el, L=ds / 2, ANGLE=float(el.ANGLE) * ds / 2 / L, **kw0_).transfer_matrix(0.0), dtype=float)[:4, :4]
            M = Mh @ M
            for i in range(nslice):
                ops.append((M, 7, (h * ds, 0.0)))
                M = Mh @ Mh if i < nslice - 1 else Mh
            if bend_edge:
                ops.append((M, 8, (-h, 0.0)))     # 出口の端: H = −(h/2) px y²
                M = np.eye(4)
            continue
        if tn == 'Sextupole':
            k2l = float(el.K2) * (L if L > 0 else 1.0)
        elif tn == 'Octupole':
            k3l = float(el.K3) * (L if L > 0 else 1.0)
        elif tn == 'Multipole':
            knl = getattr(el, 'KnL', {}) or {}
            ksl = getattr(el, 'KsL', {}) or {}
            if any(v for v in ksl.values() if v) or float(getattr(el, 'ROTATION', 0.0) or 0.0):
                raise NotImplementedError('paraxial_ops: スキュー・回転つきの MULT: %s' % el.name)
            k2l, k3l = float(knl.get(2, 0.0)), float(knl.get(3, 0.0))
            if (k2l or k3l) and L > 0 and knl.get(1, 0.0):
                raise NotImplementedError('paraxial_ops: K1≠0 の厚い MULT: %s' % el.name)
        if not (k2l or k3l):
            M = lin(el) @ M
            continue
        if L > 0:
            h = L / nslice
            D = lin(Drift(name='_d', L=h / 2))
            D2 = D @ D
            M = D @ M
            for i in range(nslice):
                ops.append((M, 6, (k2l / nslice, k3l / nslice)))
                M = D2 if i < nslice - 1 else D
        else:
            M = lin(el) @ M
            ops.append((M, 6, (k2l, k3l)))
            M = np.eye(4)
    ops.append((M, 0, (0.0, 0.0)))
    return ops


def paraxial_tracker(lat, nslice=40, bend_cubic=False, bend_edge=False):
    """paraxial_ops から 1 周の写像 tracker(coords (n × 6)) → coords を作る (z, δ は変えない)"""
    ops = paraxial_ops(lat, nslice, bend_cubic, bend_edge)

    def tr(c):
        c = np.array(c, dtype=float)
        v = c[:, :4].copy()
        for M, kind, (k2l, k3l) in ops:
            v = v @ M.T
            if kind == 7:                       # 曲率の 3 次の項の流れ (RK4、1 区間 hds = h·ds)
                hds = k2l
                def f(u):
                    return np.column_stack([u[:, 0] * u[:, 1], -0.5 * (u[:, 1] ** 2 + u[:, 3] ** 2), u[:, 0] * u[:, 3], np.zeros(len(u))])
                a1 = f(v); a2 = f(v + 0.5 * hds * a1); a3 = f(v + 0.5 * hds * a2); a4 = f(v + hds * a3)
                v = v + hds / 6.0 * (a1 + 2 * a2 + 2 * a3 + a4)
            elif kind == 8:                     # 硬い端 (sector、E=0) の 3 次の項: H = (c/2) px y² の流れ (厳密): x += c y²/2、py −= c px y
                c_ = k2l
                v[:, 0] += 0.5 * c_ * v[:, 2] ** 2
                v[:, 3] -= c_ * v[:, 1] * v[:, 2]
            elif kind:
                x, y = v[:, 0], v[:, 2]
                if k2l:
                    v[:, 1] -= 0.5 * k2l * (x * x - y * y)
                    v[:, 3] += k2l * x * y
                if k3l:
                    v[:, 1] -= k3l / 6.0 * (x ** 3 - 3 * x * y * y)
                    v[:, 3] -= k3l / 6.0 * (y ** 3 - 3 * x * x * y)
        c[:, :4] = v
        return c
    return tr


# ---------------------------------------------------------------------------------------------
# まとめ: ラティスから直接 (追跡 + NAFF + フィット) の ADTS
# ---------------------------------------------------------------------------------------------
def lattice_start(lat, delta=0.0, rows=None):
    """リング先頭の (twiss0 = ((αx, βx), (αy, βy)), orbit = (x, px, y, py) の閉軌道, 線形のチューンの端数 guess, 線形のチューン (νx, νy))。
    delta≠0 ならそのエネルギーの CALC (閉軌道・光学)"""
    rows = rows if rows is not None else lat.twiss_table(delta=delta)
    r0 = rows[0]
    tw = ((float(r0.alfx), float(r0.betx)), (float(r0.alfy), float(r0.bety)))
    orb = (float(r0.x_cod), float(r0.px_cod), float(r0.y_cod), float(r0.py_cod))
    nu = (float(rows[-1].mux) / TWO_PI, float(rows[-1].muy) / TWO_PI)
    return tw, orb, (float(fr(nu[0])), float(fr(nu[1]))), nu


def direct_adts(lat, momentum_eV, jmax, delta=0.0, nturn=1024, order=2, nstep=5, window=1, phase=(0.0, 0.0), tracker=None, halfwidth=0.02, rows=None):
    """ラティス lat の δ における ADTS を追跡 + NAFF + 多項式フィットで直接求める (J の格子 (nstep+1)² 点 × ±位相の対、RF 止め、放射なし)。
    戻り dict: adts_from_tunes の結果 + 'nu_lin'。xy は ∂νx/∂Jy と ∂νy/∂Jx の平均"""
    tw, orb, g, nu = lattice_start(lat, delta, rows)
    jx, jy = grid_actions(jmax, nstep)
    JX, JY = np.r_[jx, jx], np.r_[jy, jy]
    px = np.r_[np.full(len(jx), phase[0]), np.full(len(jx), phase[0] + math.pi)]
    py = np.r_[np.full(len(jx), phase[1]), np.full(len(jx), phase[1] + math.pi)]
    c = initial_coords(tw, JX, JY, px, py, delta=delta, orbit=orb)
    rec = track_record(lat.elements, momentum_eV, c, nturn, tracker=tracker)
    nuh = record_tunes(rec, tw, window=window, guess=g, halfwidth=halfwidth, orbit=None)
    nuh = pair_average(nuh)
    r = adts_from_tunes(jx, jy, nuh, order=order, nu_lin=g)
    r['xy_mean'] = 0.5 * (r['xy'] + r['yx'])
    r['nu_lin'] = nu
    return r


def tune_vs_delta(lat, momentum_eV, deltas, j=1e-12, nturn=1024, window=1, halfwidth=0.05):
    """微小振幅 (J=j、x・y とも) の NAFF のチューンの δ 依存。粒子は δ ごとの閉軌道から出す。戻り (n, 2) のチューン (符号つき端数) と CALC の線形のチューン (n, 2)"""
    from .track6n import track_particles_fast
    out, ref = [], []
    for d in deltas:
        rows = lat.twiss_table(delta=float(d))
        tw, orb, g, nu = lattice_start(lat, float(d), rows)
        c = initial_coords(tw, [j], [j], [0.3], [0.9], delta=float(d), orbit=orb)
        rec = track_record(lat.elements, momentum_eV, c, nturn)
        t = record_tunes(rec, tw, guess=g, halfwidth=halfwidth)
        out.append(t[0])
        ref.append(nu)
    return np.array(out), np.array(ref)
