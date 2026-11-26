# -*- coding: utf-8 -*-
"""Quantitative link between frequency maps, resonance driving terms (RDT) and ADTS.
周波数マップと共鳴駆動項 (RDT)・ADTS の定量的な接続 (v2.3.0 S1、K64 にない新機能。互換モードによらず同じ)。

目的: 周波数マップ (freqmap.py) の D の大きい領域が、RDT (resdrive.resonance_terms) の大きさと ADTS (resdrive.adts) から予想される共鳴の
位置と太さに入るか、を数で調べる。共鳴 m νx + n νy = p のまわりの「その振幅でその共鳴に入る」初期条件 (x, y) の領域を計算する。

式 (1 周ごとの写像を作用・角度で見る。H = 2π[Q·J + ½ J·A·J] + V(J, φ)、A = ADTS の行列 (∂νi/∂Jj、1/(m·rad))):
  1 周の生成関数の 6 極の 1 次の項 V = Σ_q 2 Re[ h_q (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) e^{i(a1 φx + a2 φy)} ] (調和 a = (j−k, l−m)。
  h_q は resonance_terms の h_jklm0)。調和 a = g (m0, n0) (g = gcd) の項をまとめて F(J) = Σ_q h_q (2Jx)^.. (2Jy)^.. とすると V = 2|F| cos(a·φ + arg F)。
  共鳴の方向 (m0, n0) の作用 I と、それと直交する保存量 K への変換 (Jx = m0 I + u K、Jy = n0 I + v K、m0 v − n0 u = 1) で、
  H = 2π∫Δ0(I) dI + V cos(g ψ + arg F)、ψ = m0 φx + n0 φy、Δ0(J) = m0 νx(J) + n0 νy(J) − p/g、νi(J) = νi0 + Σj Aij Jj。
  Δ0 は I について線形で、傾き G0 = (m0, n0)·A·(m0, n0)ᵀ。振り子の式: 共鳴の位置 I_r = I − Δ0/G0 (J_r = J + (I_r − I)(m0, n0)、成分が負なら共鳴は存在しない)、
  捕らえられる条件 π Δ0²/|G0| < |V(J_r)| (1 − sgn(G0) cos θ)、θ = a·φ0 + arg F(J_r) (φ0 = 初期条件の角度)。位相を問わない外側の包絡 (cos θ = −sgn(G0) の最良の位相) は
  |Δ0| < W = √(2 |G0| V(J) / π) (チューン平面の幅。V は 2|F|)。
  φ0: 位置 loc で px = py = 0 の初期条件は tan θ = −α、θ = φ + μ(loc) (x < 0 なら π を足す)。

使い方:
  pr = ResonancePredictor.from_lattice(lat, loc_twiss, rows)   (lat = Lattice、rows = lat.twiss_table(delta))
  m = pr.predict(x, y)       → {'any': bool 配列, 'envelope': bool 配列, 'by_harmonic': {(a1, a2): ...}}
  coverage_stats(D, mask, grid, thr=-5.0, nshift=...)  共鳴の予想の領域に D > thr の点が入る割合と、領域を格子の上でずらした対照との比
"""
import math

import numpy as np

from . import resdrive as _R

# 6 極 (3 次) の単項式: (j, k, l, m) の和 3
_MONOS3 = [(j, k, l, m, 0) for j in range(4) for k in range(4) for l in range(4) for m in range(4) if j + k + l + m == 3]


def _canon(a1, a2):
    """調和 (a1, a2) と (−a1, −a2) は共役。a1 > 0、または a1 = 0 で a2 > 0 の方を代表にする"""
    return a1 > 0 or (a1 == 0 and a2 > 0)


def first_order_terms(lat, rows=None, delta=0.0, bend_cubic=None, amin=1e-6):
    """Collect the first-order RDTs of sextupoles (and third order of bends) per harmonic, as {(a1, a2): [(j, k, l, m, h)]}.
    6 極 (と偏向磁石の 3 次) の 1 次の RDT を調和ごとにまとめる。戻り {(a1, a2): [(j, k, l, m, h)]} (代表の調和だけ。|h| < amin × (最大の |h|) の項は除く)"""
    res = _R.resonance_terms(lat, terms=_MONOS3, rows=rows, delta=delta, bend_cubic=bend_cubic)
    out = {}
    hmax = max(abs(v) for v in res.total.values())
    for nm, q in zip(res.names, res.terms):
        j, k, l, m, _p = q
        a = (j - k, l - m)
        if not _canon(*a):
            continue
        h = res.total[nm]
        if abs(h) < amin * hmax:
            continue
        out.setdefault(a, []).append((j, k, l, m, h))
    return out, res


def second_order_terms(lat, rows=None, amin=1e-3, order2=True, **kw):
    """Collect the fourth-order RDTs h_jklm0 per harmonic, excluding the diagonal ADTS harmonic.
    v2.3.0 R5: 4 次の RDT h_jklm0 (j+k+l+m = 4。6 極の 2 次 + 8 極 + 運動学 + フリンジ + 偏向磁石の 4 次の項、`rdt2nd.rdt4`) を調和 (j−k, l−m) ごとにまとめる。
    戻り {(a1, a2): [(j, k, l, m, h)]} (代表の調和。対角 (0, 0) は ADTS なので除く。|h| < amin × (最大の |h|) は除く)。
    注意: 写像の 4 次の生成関数の係数で、正規形の補正 (3 次の項の非共鳴の摂動、`rdt2nd.nf_correction` は対角だけ) は非対角の調和には入れていない"""
    from . import rdt2nd
    r = rdt2nd.rdt4(lat, rows=rows, order2=order2, **kw)
    out = {}
    hmax = max(abs(v) for v in r.total.values())
    for q in rdt2nd.MONOS4:
        a = (q[0] - q[1], q[2] - q[3])
        if a == (0, 0) or not _canon(*a):
            continue
        h = r.total['h%d%d%d%d0' % q]
        if abs(h) < amin * hmax:
            continue
        out.setdefault(a, []).append((q[0], q[1], q[2], q[3], h))
    return out, r


def _egcd(a, b):
    if b == 0:
        return (1, 0) if a >= 0 else (-1, 0)
    q, r = divmod(a, b)
    x, y = _egcd(b, r)
    return y, x - q * y


class ResonancePredictor:
    """Prediction region for each resonance, holding the harmonics, ADTS matrix A, linear tunes nu0 and Twiss at the location.
    共鳴ごとの予想の領域。属性: harms {(a1, a2): [(j, k, l, m, h)]}、A (2 × 2)、nu0 (νx0, νy0) 整数部を含む線形のチューン、
    tw (αx, βx, αy, βy, μx, μy: 位置 loc の Twiss)"""

    def __init__(self, harms, A, nu0, tw):
        self.harms = harms
        self.A = np.array(A, dtype=float)
        self.nu0 = tuple(nu0)
        self.tw = tuple(tw)

    @classmethod
    def from_lattice(cls, lat, loc, rows=None, delta=0.0, bend_cubic=None, order_adts=True, order2=False, amin2=1e-3, model=None):
        """Build the predictor from a Lattice at a location, optionally adding fourth-order RDTs or returning the S3 model.
        lat: Lattice。loc: 位置 (1 始まり、要素の入口)。rows: lat.twiss_table(delta=delta) (省略可)。RDT は resonance_terms、A は adts().total。
        order2=True (v2.3.0 R5): 4 次の RDT (second_order_terms) も調和ごとに足す (1 次だけでは予想できなかった νx−νy・2νx±2νy・4νx・νx±3νy などの島の領域)
        model='s3' (v2.3.0 S3): `fmappred.ResonancePredictorS3` を返す (4 次の h に非共鳴の 3 次の補正、連続のスコア、信頼度のフラグ。docs/25 §12)"""
        if model == 's3':
            from .fmappred import ResonancePredictorS3
            return ResonancePredictorS3.from_lattice(lat, loc, rows, delta=delta, bend_cubic=bend_cubic)
        if rows is None:
            rows = lat.twiss_table(delta=delta)
        harms, _res = first_order_terms(lat, rows, delta, bend_cubic)
        if order2:          # v2.3.0 R5: 4 次の RDT (6 極の 2 次など) を同じ調和にまとめる (δ = 0 だけ)
            h2, _r2 = second_order_terms(lat, rows, amin=amin2)
            for a, t in h2.items():
                harms[a] = harms.get(a, []) + t
        a = _R.adts(lat, rows=rows, delta=delta, bend_cubic=bend_cubic)
        tot = a.total
        A = [[tot['xx'], tot['xy']], [tot['xy'], tot['yy']]]
        r = rows[loc - 1]
        tw = (float(r.alfx), float(r.betx), float(r.alfy), float(r.bety), float(r.mux), float(r.muy))
        nu0 = (float(rows[-1].mux) / (2 * math.pi), float(rows[-1].muy) / (2 * math.pi))
        return cls(harms, A, nu0, tw)

    # ------------------------------------------------------------------
    def actions(self, x, y):
        """Return the actions (Jx, Jy) of initial conditions (x, y) with px = py = 0 and the angles at the lattice start.
        px = py = 0 の初期条件 (x, y) [m] の作用 (Jx, Jy) と、基準 (リングの先頭) での角度 (φx, φy)"""
        ax, bx, ay, by, mux, muy = self.tw
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        Jx = (1 + ax * ax) * x * x / (2 * bx)
        Jy = (1 + ay * ay) * y * y / (2 * by)
        thx = np.arctan(-ax) + np.where(x < 0, math.pi, 0.0)
        thy = np.arctan(-ay) + np.where(y < 0, math.pi, 0.0)
        return Jx, Jy, thx - mux, thy - muy

    def tunes_lin(self, Jx, Jy):
        A = self.A
        return (self.nu0[0] + A[0, 0] * Jx + A[0, 1] * Jy, self.nu0[1] + A[1, 0] * Jx + A[1, 1] * Jy)

    def _F(self, terms, Jx, Jy):
        F = 0j
        for (j, k, l, m, h) in terms:
            F = F + h * (2 * Jx) ** ((j + k) / 2.0) * (2 * Jy) ** ((l + m) / 2.0)
        return F

    def predict(self, x, y, nu=None):
        """Return, for initial conditions (x, y), whether each harmonic traps them (with phase) or contains them in its envelope.
        初期条件 (x, y) の配列 (同じ形) に対し、各調和で「捕らえられる (位相つき)」「包絡 (位相を問わない)」の真偽の配列と、Δ0・W・p を返す。
        nu: (νx, νy) の配列 (整数部を含む) を渡すと、ν(J) を線形の式の代わりにそれ (実測のチューン) にする (G0 は A から)。
        戻り {'any': 位相つきの和集合、'envelope': 包絡の和集合、'by_harmonic': {a: dict(trapped, envelope, delta0, W, p, G0)}}"""
        x = np.asarray(x, dtype=float)
        shape = x.shape
        Jx, Jy, phx, phy = self.actions(x, y)
        if nu is None:
            nux, nuy = self.tunes_lin(Jx, Jy)
        else:
            nux, nuy = nu
        anyt = np.zeros(shape, dtype=bool)
        anye = np.zeros(shape, dtype=bool)
        by = {}
        A = self.A
        for a, terms in self.harms.items():
            g = math.gcd(abs(a[0]), abs(a[1]))
            m0, n0 = a[0] // g, a[1] // g
            G0 = m0 * m0 * A[0, 0] + 2 * m0 * n0 * A[0, 1] + n0 * n0 * A[1, 1]
            s = (m0 * nux + n0 * nuy) * g
            p = np.round(s)
            d0 = (s - p) / g
            if abs(G0) < 1e-12:
                trapped = np.zeros(shape, dtype=bool)
                env = np.zeros(shape, dtype=bool)
                W = np.full(shape, np.nan)
            else:
                dI = -d0 / G0
                Jxr, Jyr = Jx + dI * m0, Jy + dI * n0
                ok = (Jxr > 0) & (Jyr > 0)
                F = self._F(terms, Jx, Jy)          # V は初期条件の J で評価する (島が J に比べて小さいときの振り子の近似。共鳴の位置 J_r で評価すると、
                                                    # J から遠い共鳴 (|J_r − J| ≳ J) の V を過大に見積もる)
                V = 2.0 * np.abs(F)
                th = a[0] * phx + a[1] * phy + np.angle(F)
                sg = 1.0 if G0 > 0 else -1.0
                en = math.pi * d0 * d0 / abs(G0)
                trapped = ok & (en < V * (1.0 - sg * np.cos(th)))
                W = np.sqrt(2.0 * abs(G0) * V / math.pi)
                env = ok & (np.abs(d0) < W)
            by[a] = dict(trapped=trapped, envelope=env, delta0=d0, W=W, p=p, G0=G0, g=g)
            anyt |= trapped
            anye |= env
        return dict(any=anyt, envelope=anye, by_harmonic=by)


def coverage_stats(D, mask, grid, thr=-5.0, nshift=400, seed=1):
    """Compute coverage statistics of a resonance prediction mask against D: area, hit count, coverage and enrichment.
    D (点のリスト、格子の順 = meshgrid(ij) の ravel、(nx, ny))、mask (同じ形の真偽の配列、共鳴の予想の領域)。
    戻り dict: area (領域の割合)、n_hi (D > thr の点数)、cov (D > thr の点が領域に入る割合)、enrich (cov/area)、p_in・p_out (領域の中・外の D > thr の割合)、
    shift: 領域を格子の上で巡回ずらし (nshift 回、領域と同じ形・面積) して得た cov の平均・標準偏差・最大、pval (ずらしの cov ≥ cov の割合)。nan の点は除く"""
    nx, ny = grid[0], grid[1]
    D = np.asarray(D, dtype=float).reshape(nx, ny)
    M = np.asarray(mask, dtype=bool).reshape(nx, ny)
    ok = np.isfinite(D)
    hi = ok & (D > thr)
    n_ok = int(ok.sum())
    n_hi = int(hi.sum())
    area = float((M & ok).sum()) / max(n_ok, 1)
    cov = float((M & hi).sum()) / max(n_hi, 1)
    p_in = float((M & hi).sum()) / max(int((M & ok).sum()), 1)
    p_out = float((~M & hi).sum()) / max(int((~M & ok).sum()), 1)
    rng = np.random.default_rng(seed)
    covs = []
    for _ in range(nshift):
        sx, sy = int(rng.integers(0, nx)), int(rng.integers(0, ny))
        Ms = np.roll(np.roll(M, sx, axis=0), sy, axis=1)
        covs.append(float((Ms & hi).sum()) / max(n_hi, 1))
    covs = np.array(covs)
    return dict(n_ok=n_ok, n_hi=n_hi, area=area, cov=cov, enrich=(cov / area if area > 0 else float('nan')), p_in=p_in, p_out=p_out,
                shift_mean=float(covs.mean()), shift_std=float(covs.std()), shift_max=float(covs.max()),
                pval=float((covs >= cov - 1e-12).mean()) if n_hi else float('nan'))


# ---------------------------------------------------------------------------------------------------------------------
# 共鳴線までの距離・図・誤差の与え方 (v2.3.0 S1)
# ---------------------------------------------------------------------------------------------------------------------
FIRST_ORDER_MN = [(3, 0), (1, 0), (1, 2), (1, -2)]          # 6 極 (normal) の 1 次の RDT が駆動する調和 (代表)


def lines_distance(nux, nuy, mn_list):
    """Return the shortest tune-plane distance from each point (nux, nuy) to the resonance lines m nux + n nuy = p in mn_list.
    各点 (νx, νy) (端数でも整数部つきでもよい) から、mn_list = [(m, n)] の共鳴線 m νx + n νy = p (どの整数 p でも) までの最短距離 (チューン平面の距離)"""
    nux = np.asarray(nux, dtype=float)
    nuy = np.asarray(nuy, dtype=float)
    best = np.full(nux.shape, np.inf)
    for (m, n) in mn_list:
        v = m * nux + n * nuy
        best = np.minimum(best, np.abs(v - np.round(v)) / math.hypot(m, n))
    return best


def other_mn(order=4, exclude=FIRST_ORDER_MN):
    """Return the (m, n) up to a given order that the first-order RDT does not drive (difference resonance, fourth order, etc.).
    order 次以下で、1 次の RDT が駆動しない (m, n) (差の共鳴 νx − νy、4 次など)。(m, n) と (−m, −n) は 1 つ"""
    from .freqmap import _mn_list
    ex = set(exclude) | set((-a, -b) for a, b in exclude)
    return [mn for mn in _mn_list(order) if mn not in ex and (-mn[0], -mn[1]) not in ex]


def apply_alignment_errors(setup, seed=1, quad_dxy=0.0, quad_rot=0.0, sext_dxy=0.0, bend_dxy=0.0, bend_rot=0.0, sext_rot=0.0, oct_dxy=0.0, streams='single'):
    """Return a copy of a TrackSetup with random position offsets (DX, DY) and rolls applied independently to each element.
    TrackSetup の複製に、乱数 (正規分布、標準偏差 = 引数) の位置のずれ (DX・DY [m])・ロール (四極・6 極の ROTATION・偏向磁石の DROTATE [rad]) の誤差を、
    要素ごとに別々に (周期性を使わずに) 与える。戻り (新しい TrackSetup, 辞書: 要素の数と誤差の rms)。元の setup は変えない (要素を複製して書き換える)。
    閉軌道は frequency_map が (RFSW OFF では) 4 次元の Newton 法で求め直す。Twiss (setup.rows) は誤差なしの設計値のまま (frequency_map の規格化は標本から求める)。
    streams: 'single' (S1 の既定、乱数の列が 1 本で、要素の順に四極のずれ・ロール、6 極のずれ、偏向磁石のずれ・ロールを引く。どれかの引数を 0 にすると列がずれる)、
    'type' (v2.3.0 S2: 誤差の種類 (quad_dxy・quad_rot・sext_dxy・sext_rot・oct_dxy・bend_dxy・bend_rot) ごとに別の乱数の列 (seed と種類で決まる) を使う。
    ある種類だけを入れた結果が、全部入れた結果の中のその種類の誤差と同じになり (切り分けに使う)、強さを k 倍すると同じ乱数の k 倍になる (強さの走査に使う)。
    sext_rot・oct_dxy は 'type' のときだけ (単数の列では S1 の結果を変えないため)。八極 (Octupole) は oct_dxy のとき 'DX'・'DY' を与える)"""
    import copy
    st = setup.copy()
    st.els = list(setup.els)
    cnt = {'quad': 0, 'sext': 0, 'bend': 0}
    if streams == 'type':
        cnt['oct'] = 0
        kinds = ('quad_dxy', 'quad_rot', 'sext_dxy', 'sext_rot', 'oct_dxy', 'bend_dxy', 'bend_rot')
        rg = {k: np.random.default_rng([int(seed), j]) for j, k in enumerate(kinds)}
        sg = dict(quad_dxy=quad_dxy, quad_rot=quad_rot, sext_dxy=sext_dxy, sext_rot=sext_rot, oct_dxy=oct_dxy, bend_dxy=bend_dxy, bend_rot=bend_rot)

        def nrm(k, n=1):
            # 標準偏差 0 でも乱数を引く (列の位置を変えない)
            v = rg[k].normal(0.0, 1.0, n) * sg[k]
            return [float(t) for t in v]
        for i, e in enumerate(st.els):
            tn = type(e).__name__
            if tn == 'Quadrupole':
                dxy, rot = nrm('quad_dxy', 2), nrm('quad_rot')[0]
                if not (quad_dxy or quad_rot):
                    continue
                e2 = copy.copy(e)
                e2.DX, e2.DY = dxy
                e2.ROTATION = float(getattr(e, 'ROTATION', 0.0)) + rot
                cnt['quad'] += 1
            elif tn == 'Sextupole':
                dxy, rot = nrm('sext_dxy', 2), nrm('sext_rot')[0]
                if not (sext_dxy or sext_rot):
                    continue
                e2 = copy.copy(e)
                e2.DX, e2.DY = dxy
                e2.ROTATION = float(getattr(e, 'ROTATION', 0.0)) + rot
                cnt['sext'] += 1
            elif tn == 'Octupole':
                dxy = nrm('oct_dxy', 2)
                if not oct_dxy:
                    continue
                e2 = copy.copy(e)
                e2.DX, e2.DY = dxy
                cnt['oct'] += 1
            elif tn == 'Bend' and getattr(e, 'L', 0.0) > 0:
                dxy, rot = nrm('bend_dxy', 2), nrm('bend_rot')[0]
                if not (bend_dxy or bend_rot):
                    continue
                e2 = copy.copy(e)
                e2.DX, e2.DY = dxy
                e2.DROTATE = rot
                cnt['bend'] += 1
            else:
                continue
            st.els[i] = e2
        return st, cnt
    rng = np.random.default_rng(seed)
    for i, e in enumerate(st.els):
        tn = type(e).__name__
        if tn == 'Quadrupole' and (quad_dxy or quad_rot):
            e2 = copy.copy(e)
            e2.DX, e2.DY = float(rng.normal(0, quad_dxy)), float(rng.normal(0, quad_dxy))
            e2.ROTATION = float(getattr(e, 'ROTATION', 0.0)) + float(rng.normal(0, quad_rot))
            cnt['quad'] += 1
        elif tn == 'Sextupole' and sext_dxy:
            e2 = copy.copy(e)
            e2.DX, e2.DY = float(rng.normal(0, sext_dxy)), float(rng.normal(0, sext_dxy))
            cnt['sext'] += 1
        elif tn == 'Bend' and getattr(e, 'L', 0.0) > 0 and (bend_dxy or bend_rot):
            e2 = copy.copy(e)
            e2.DX, e2.DY = float(rng.normal(0, bend_dxy)), float(rng.normal(0, bend_dxy))
            e2.DROTATE = float(rng.normal(0, bend_rot))
            cnt['bend'] += 1
        else:
            continue
        st.els[i] = e2
    return st, cnt


def plot_with_prediction(res, masks, path, vmin=-8.0, vmax=-2.0, title=None, axes=None):
    """Plot the D map of a frequency map (x-y, x-delta or y-delta) with contours of the predicted resonance regions overlaid.
    周波数マップ (res: FMapResult、x–y・x–δ・y–δ のどれか) の D の図に、共鳴の予想の領域の輪郭を重ねる。
    masks: [(ラベル, 真偽の配列 (点の順), 色, 線種)]。axes: 軸の名前 2 つ ('x','y')・('x','dp')・('y','dp') (None なら格子の長さが 1 でない 2 つ)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .freqmap import _plane_names
    names = axes or _plane_names(res.grid)
    cols = {'x': (res.x * 1e3, 'x [mm]'), 'y': (res.y * 1e3, 'y [mm]'), 'dp': (res.dp * 100, 'dp [%]')}
    shape = tuple(g for g in res.grid if g > 1) or (res.grid[0], res.grid[1])
    if len(shape) != 2:
        shape = (res.grid[0], res.grid[1])
    X = cols[names[0]][0].reshape(shape)
    Y = cols[names[1]][0].reshape(shape)
    D = res.D.reshape(shape)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    pc = ax.pcolormesh(X, Y, np.ma.masked_invalid(D), cmap='turbo', vmin=vmin, vmax=vmax, shading='nearest')
    lost = (res.lost_turn > 0).reshape(shape)
    if lost.any():
        ax.scatter(X[lost], Y[lost], c='k', s=6, marker='x', linewidths=0.6)
    handles = []
    for (lab, m, col, ls) in masks:
        M = np.asarray(m, dtype=float).reshape(shape)
        if M.max() > 0 and M.min() < 1:
            ax.contour(X, Y, M, levels=[0.5], colors=[col], linestyles=ls, linewidths=1.3)
        handles.append(plt.Line2D([], [], color=col, ls=ls, label=lab))
    if handles:
        ax.legend(handles=handles, fontsize=8, loc='upper right')
    ax.set_xlabel(cols[names[0]][1])
    ax.set_ylabel(cols[names[1]][1])
    ax.set_title(title or 'D (colour) and predicted resonance regions (contours)', fontsize=9)
    fig.colorbar(pc, ax=ax, label='D = log10 |dnu|')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
