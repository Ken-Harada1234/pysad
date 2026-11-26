# -*- coding: utf-8 -*-
"""Helpers for trusting frequency maps: error statistics, feed-down estimates and warnings.
周波数マップを信頼して使うための補助 (v2.3.0 S2、K64 にない新機能。互換モードによらず同じ)。

  誤差 (位置のずれ・ロール) の統計と、フィードダウンの見積もり
    closed_orbit_all(setup, dp)            全要素の入口の閉軌道 (RFSW OFF の 4 次元)
    linear_tunes(setup, dp)                閉軌道まわりの 1 周の行列 (4 × 4) の固有値から、結合を含むチューン (端数)
    feeddown_estimate(setup, dp)           6 極のフィードダウン (閉軌道・ずれ) と四極のロールから、チューンのずれ・結合 |c−| を見積もる
    errors_summary(setup)                  閉軌道の rms・最大値、チューン、結合
  周波数マップの結果の統計と、周回数に依らない指標
    map_stats(res, thr)                    D > thr の割合、失われた点の割合、D の百分位
    rate_index(res, nturn)                 拡散率 log10(|Δν| / (nturn/2)^κ) など (κ は scaling_exponent で求める)
  使えない条件・信頼できない条件の警告
    check_conditions(setup, nturn, res)    警告の文字列のリスト (周回数、チューンと整数・半整数・低次の共鳴の近さ、周期性、RF・放射・FLUC、…)
"""
import math

import numpy as np

from . import trksetup as _ts
from .track6n import state_out


# ---------------------------------------------------------------------------------------------------------------------
# 閉軌道・線形のチューン
# ---------------------------------------------------------------------------------------------------------------------
def closed_orbit_all(setup, dp=0.0, locs=None):
    """Compute the 4D closed orbit with RFSW OFF at delta = dp by Newton's method and transport it to the entrance of each element.
    RFSW OFF の δ = dp の閉軌道 (4 次元) を、位置 1 で Newton 法で求め、全要素の入口 (locs、1 始まり。既定は 1〜nlat) へ運ぶ。戻り (len(locs), 6)"""
    c1 = setup.refine_orbit(np.zeros(6), dp=dp, loc=1)
    if locs is None:
        locs = list(range(1, setup.nlat + 1))
    out = np.zeros((len(locs), 6))
    tr = setup.make_run(c1[None, :], alost=1e30, zlost=1e30)
    cur = 1
    for k, l in enumerate(locs):
        if l > cur:
            with np.errstate(all='ignore'):
                tr.segment(cur, l)
            cur = l
        out[k] = state_out(tr.S, tr.n)[0] if l > 1 else c1
    return out


def _jacobian4(setup, c, loc=1, h=1e-8):
    pts = [c.copy()]
    for k in range(4):
        p = c.copy()
        p[k] += h
        pts.append(p)
    out = setup._map_batch(np.array(pts), loc)
    return np.array([(out[1 + j] - out[0]) / h for j in range(4)]).T[:4, :4]


def linear_tunes(setup, dp=0.0, loc=1, design=None):
    """Compute the fractional tunes (including coupling) from the eigenvalues of the one-turn 4x4 matrix around the closed orbit.
    閉軌道まわりの 1 周の行列 (4 × 4、数値微分) の固有値から (νa, νb) (端数、design = (νx, νy) の端数に近い順に x 的・y 的)。
    戻り dict: nu (x 的、y 的の端数)、trace_x・trace_y (固有値の絶対値からのずれ 最大 = 線形の安定性)。design 省略は setup の Twiss の端数"""
    c = setup.refine_orbit(np.zeros(6), dp=dp, loc=loc if loc == 1 else 1)
    M = _jacobian4(setup, c, 1)
    w = np.linalg.eigvals(M)
    pos = [z for z in w if z.imag > 0]
    ang = sorted(math.atan2(z.imag, z.real) / (2 * math.pi) for z in pos)       # 0 < a < 0.5
    if len(ang) != 2:
        return dict(nu=(float('nan'), float('nan')), stable=False, orbit=c)
    dx, dy = (v % 1.0 for v in (design or setup.tunes_full()))
    out = []
    for d in (dx, dy):
        cand = []
        for a in ang:
            cand += [a, 1.0 - a]
        out.append(min(cand, key=lambda v: abs(v - d)))
    # x 的と y 的が同じ固有値を選ばないように (結合が強いとき)
    # v2.4.0: RAD ON (放射減衰) や RFSW ON (z との結合) の追跡の 1 周の行列は固有値の絶対値が 1 から少し下がる (1 周の減衰 ~1e-3)。それは不安定ではない。
    # 不安定 = 絶対値が 1 を超える (成長)、または減衰が 1 周で 5% を超える (RAD ON でも起こらない) / RAD・RFSW OFF では 1e-6 を超えるずれ
    aw = np.abs(w)
    low = 0.05 if (getattr(setup, 'rad', False) or getattr(setup, 'rfsw', False)) else 1e-6
    stable = bool(np.max(aw) - 1.0 < 1e-6 and np.min(aw) > 1.0 - low)
    return dict(nu=tuple(out), stable=stable, orbit=c)


# ---------------------------------------------------------------------------------------------------------------------
# フィードダウンの見積もり
# ---------------------------------------------------------------------------------------------------------------------
def feeddown_estimate(setup, dp=0.0, orbit=None):
    """Estimate tune shifts and coupling strength |c-| from sextupole/octupole feed-down and quadrupole roll (thin-lens formulas).
    6 極 (と 8 極) のフィードダウンによる実効の四極・スキュー四極と、四極のロールによるスキュー四極から、線形のチューンのずれ (Δνx, Δνy) と結合の強さ |c−|
    (差の共鳴) を見積もる。式 (薄いレンズ。6 極の中心の β は入口・出口の平均):
      実効の四極の積分強さ ΔK1L = K2L·x_rel (x_rel = 6 極の中心の閉軌道 − 6 極のずれ DX)、スキュー四極 ΔK1sL = K2L·y_rel、四極のロール θ → 2·K1L·θ。
      Δνx = Σ βx ΔK1L/(4π)、Δνy = −Σ βy ΔK1L/(4π)、c− = (1/2π) Σ ΔK1sL √(βxβy) exp(i[μx − μy − (νx − νy − p)·2πs/C])
    orbit: closed_orbit_all の結果 (省略すると計算する)。戻り dict: dnux、dnuy、cminus (複素数の和)、cminus_fd (6 極のフィードダウン分)、cminus_roll (四極のロール分)、
    rms_orbit_sext (6 極の中心の閉軌道の rms (x, y))"""
    els = setup.els
    rows = setup.rows
    if orbit is None:
        orbit = closed_orbit_all(setup, dp)
    nux, nuy = setup.tunes_full()
    p = round(nux - nuy)
    C = setup.circ
    dnux = dnuy = 0.0
    cfd = 0j
    croll = 0j
    xr, yr = [], []
    for i, e in enumerate(els):
        tn = type(e).__name__
        L = float(getattr(e, 'L', 0.0) or 0.0)
        if tn == 'Sextupole' and L > 0:
            k2l = float(e.K2) * L
            xm = 0.5 * (orbit[i, 0] + orbit[i + 1, 0]) - float(getattr(e, 'DX', 0.0) or 0.0)
            ym = 0.5 * (orbit[i, 2] + orbit[i + 1, 2]) - float(getattr(e, 'DY', 0.0) or 0.0)
            xr.append(xm)
            yr.append(ym)
            kn, ks = k2l * xm, k2l * ym
        elif tn == 'Octupole' and L > 0:
            k3l = float(e.K3) * L
            xm = 0.5 * (orbit[i, 0] + orbit[i + 1, 0]) - float(getattr(e, 'DX', 0.0) or 0.0)
            ym = 0.5 * (orbit[i, 2] + orbit[i + 1, 2]) - float(getattr(e, 'DY', 0.0) or 0.0)
            kn, ks = 0.5 * k3l * (xm * xm - ym * ym), k3l * xm * ym         # 八極のフィードダウン: 実効の四極 K3L(x²−y²)/2、スキュー四極 K3L·x·y
        elif tn == 'Quadrupole' and L > 0:
            k1l = float(e.K1) * L
            th = float(getattr(e, 'ROTATION', 0.0) or 0.0)
            kn, ks = 0.0, 2.0 * k1l * th
        else:
            continue
        bx = 0.5 * (rows[i][1] + rows[i + 1][1])
        by = 0.5 * (rows[i][3] + rows[i + 1][3])
        smid = 0.5 * (rows[i][0] + rows[i + 1][0])
        mux = 0.5 * (rows[i][7] + rows[i + 1][7])
        muy = 0.5 * (rows[i][8] + rows[i + 1][8])
        dnux += bx * kn / (4 * math.pi)
        dnuy -= by * kn / (4 * math.pi)
        z = ks * math.sqrt(bx * by) * np.exp(1j * (mux - muy - (nux - nuy - p) * 2 * math.pi * smid / C)) / (2 * math.pi)
        if tn in ('Sextupole', 'Octupole'):
            cfd += z
        else:
            croll += z
    return dict(dnux=dnux, dnuy=dnuy, cminus=cfd + croll, cminus_fd=cfd, cminus_roll=croll, p=p,
                rms_orbit_sext=(float(np.sqrt(np.mean(np.square(xr)))) if xr else 0.0, float(np.sqrt(np.mean(np.square(yr)))) if yr else 0.0))


def center_sextupoles(setup, dp=0.0, iters=6, tol=1e-9, kinds=('Sextupole',)):
    """Shift sextupole offsets DX and DY onto the closed orbit at their centers to remove feed-down, returning a new TrackSetup.
    フィードダウンを切る切り分け用: 6 極 (kinds に Octupole も入れれば 8 極も) のずれ DX・DY を、その中心の閉軌道 (入口・出口の平均) に合わせる (x_rel = y_rel = 0)。
    6 極のずれで閉軌道が少し変わるので、不動点の反復 (iters 回、ずれの変化が tol 以下で止める)。戻り 新しい TrackSetup (元は変えない)。
    6 極を通る粒子の軌道は閉軌道に載るので、閉軌道のずれ・四極のロールなどの他の誤差は残したまま、6 極のフィードダウン (実効の四極・スキュー四極) だけが消える"""
    import copy
    st = setup.copy()
    st.els = list(setup.els)
    idx = [i for i, e in enumerate(st.els) if type(e).__name__ in kinds and float(getattr(e, 'L', 0.0) or 0.0) > 0]
    for _ in range(iters):
        orb = closed_orbit_all(st, dp)
        change = 0.0
        for i in idx:
            e = copy.copy(st.els[i])
            nx = 0.5 * (orb[i, 0] + orb[i + 1, 0])
            ny = 0.5 * (orb[i, 2] + orb[i + 1, 2])
            change = max(change, abs(nx - e.DX), abs(ny - e.DY))
            e.DX, e.DY = float(nx), float(ny)
            st.els[i] = e
        if change < tol:
            break
    return st


def errors_summary(setup, dp=0.0, loc=None):
    """Summarize the linear properties of an error-laden TrackSetup: closed-orbit rms and maximum, tune shift, coupling, feed-down.
    誤差つきの TrackSetup の線形の性質: 閉軌道 (全要素の入口) の rms・最大 (x, y)、線形のチューンのずれ、結合 (チューンの最小の距離の見積もり)、
    フィードダウンの見積もり。loc: 閉軌道を報告する位置 (既定 2)"""
    nux, nuy = setup.tunes_full()
    try:
        orb = closed_orbit_all(setup, dp)
        lt = linear_tunes(setup, dp)
    except Exception:           # 閉軌道が求まらない (Newton 法が発散、粒子が失われる: 線形の運動が不安定で閉軌道がない)
        nan = float('nan')
        return dict(orbit_rms=(nan, nan), orbit_max=(nan, nan), orbit_loc=(nan, nan), dnu_lin=(nan, nan), stable=False, closed_orbit=False, dnu_pred=(nan, nan),
                    cminus_pred=nan, cminus_fd=nan, cminus_roll=nan, split_pred=nan, split_num=nan, nu_lin=(nan, nan), rms_orbit_sext=(nan, nan))
    fd = feeddown_estimate(setup, dp, orbit=orb)
    n = setup.nlat - 1
    rms = (float(np.sqrt(np.mean(orb[:n, 0] ** 2))), float(np.sqrt(np.mean(orb[:n, 2] ** 2))))
    mx = (float(np.max(np.abs(orb[:n, 0]))), float(np.max(np.abs(orb[:n, 2]))))
    dnu = (lt['nu'][0] - (nux % 1.0), lt['nu'][1] - (nuy % 1.0))
    # 差の共鳴の近さ: 結合のない仮の値 (設計 + フィードダウン分のずれ) の距離から、結合込みの距離
    d0 = (nux + fd['dnux']) - (nuy + fd['dnuy']) - fd['p']
    split_pred = math.sqrt(d0 * d0 + abs(fd['cminus']) ** 2)
    split_num = abs((lt['nu'][0] - lt['nu'][1]) - (fd['p'] - (math.floor(nux) - math.floor(nuy)))) if lt['stable'] else float('nan')
    return dict(orbit_rms=rms, orbit_max=mx, orbit_loc=tuple(float(v) for v in orb[(loc or 2) - 1, [0, 2]]), dnu_lin=dnu, stable=lt['stable'], closed_orbit=True,
                dnu_pred=(fd['dnux'], fd['dnuy']), cminus_pred=abs(fd['cminus']), cminus_fd=abs(fd['cminus_fd']), cminus_roll=abs(fd['cminus_roll']),
                split_pred=split_pred, split_num=split_num, nu_lin=lt['nu'], rms_orbit_sext=fd['rms_orbit_sext'])


# ---------------------------------------------------------------------------------------------------------------------
# 結果の統計
# ---------------------------------------------------------------------------------------------------------------------
def map_stats(res, thr=-5.0):
    """Compute statistics of an FMapResult: counts, lost fraction, fraction with D above a threshold, and D percentiles.
    FMapResult の統計: n、nlost、lost_frac、frac_hi (D > thr の割合、残った点のうち)、frac_hi_all (失われた点を D > thr とみなした割合)、D の中央値・最小・最大・百分位"""
    ok = np.isfinite(res.D)
    n = len(res)
    nl = int((res.lost_turn > 0).sum())
    D = res.D[ok]
    out = dict(n=n, nlost=nl, lost_frac=nl / max(n, 1), n_ok=int(ok.sum()))
    if len(D):
        out.update(frac_hi=float((D > thr).mean()), frac_hi_all=float(((D > thr).sum() + nl) / max(n, 1)), Dmed=float(np.median(D)),
                   Dmin=float(D.min()), Dmax=float(D.max()), Dq=[float(v) for v in np.percentile(D, [5, 25, 50, 75, 95])])
    else:
        out.update(frac_hi=float('nan'), frac_hi_all=1.0 if nl else float('nan'), Dmed=float('nan'), Dmin=float('nan'), Dmax=float('nan'), Dq=[float('nan')] * 5)
    return out


# ---------------------------------------------------------------------------------------------------------------------
# 周回数に依らない指標
# ---------------------------------------------------------------------------------------------------------------------
D_THRESHOLD_NTURN = 1024            # D > −5 の閾値の基準の周回数
SCALING_KAPPA = 1.0                 # 規則的な点の D の周回数の依存: 周回数 2 倍で D が下がる量 (S2 の測定: ASTRID2 で中央値 −0.96〜−1.04、256〜4096 周)


def threshold_for_nturn(nturn, thr1024=-5.0, kappa=None):
    """Convert the D threshold between regular and chaotic points to the given number of turns, using thr1024 - kappa log2(nturn/1024).
    D の閾値 (規則的/カオス的の境) の周回数への換算: 規則的な点の D は周回数が 2 倍になると約 κ 下がる (NAFF の誤差が周回数のべきで下がるため)
    ので、閾値も thr1024 − κ log2(nturn / 1024)。カオス的な点の D は周回数にほとんど依らない"""
    k = SCALING_KAPPA if kappa is None else kappa
    return thr1024 - k * math.log2(float(nturn) / D_THRESHOLD_NTURN)


# ---------------------------------------------------------------------------------------------------------------------
# 警告
# ---------------------------------------------------------------------------------------------------------------------
def _has_errors(els):
    """要素にずれ・ロールの誤差 (DX・DY・DROTATE、四極・6 極の ROTATION) があるか"""
    for e in els:
        tn = type(e).__name__
        if tn not in ('Bend', 'Quadrupole', 'Sextupole', 'Octupole'):
            continue
        if float(getattr(e, 'L', 0.0) or 0.0) <= 0.0:
            continue
        for k in ('DX', 'DY', 'DROTATE') + (('ROTATION',) if tn != 'Bend' else ()):
            if float(getattr(e, k, 0.0) or 0.0) != 0.0:
                return True
    return False


def check_setup(setup, nturn=1024, window=1, loc=None, thr_tune=0.03):
    """Return warnings about conditions known before tracking (turns, window, RFSW, RAD, FLUC, tunes near resonances, errors, orbit).
    追跡の前に分かる条件の警告 (setup と周回数・窓から)。戻り [(コード, 文)]。コード:
    nturn_short (周回数 < 512)、window_none (窓なし)、rfsw (RFSW ON)、rad (RAD ON)、fluc (FLUC ON)、tune_integer (作動点が整数・半整数に近い)、
    tune_lowres (作動点が 3 次・4 次の共鳴に近い)、errors (誤差つき: 周期性なし・位置依存)、orbit (閉軌道が収束していない)、coupling (結合が強い)、orbit_large (閉軌道が大きい)"""
    from .freqmap import resonance_check
    w = []
    if nturn < 512:
        w.append(('nturn_short', '周回数 %d が少ない: NAFF の周波数の精度が悪く D の下限が上がり、D > −5 の閾値の目安 (1024 周) と比べられない。512 周以上を使う' % nturn))
    if window == 0:
        w.append(('window_none', 'NAFF の窓なし (window=0): D の下限が −5.7 まで上がり (Hann の窓では −8 以下)、D > −5 の割合が 5 倍になる。使わない (window=1 以上)'))
    if setup.rfsw:
        w.append(('rfsw', 'RFSW ON: 空洞が働き z が 1 周ごとに 0 に戻らず、δ が (自然クロマティシティ) × z の変調を受ける。D はシンクロトロン振動 (νs 周期、周回数が 1/νs 程度以上) の影響を含み、'
                  '振幅だけの拡散の指標ではなくなる。gs の TrackParticles と同じ設定 (FFS に入ると RFSW は ON)。x–y の周波数マップは通常 NORFSW (setup.copy(rfsw=False))'))
    if setup.rad:
        w.append(('rad', 'RAD ON: 放射減衰・励起で振幅が周ごとに変わり、前半・後半のチューンの差は減衰によるチューンの動きで、D は振幅に依る拡散ではなく減衰の指標になる'
                  '(ASTRID2: D の中央値 −7.9 → −3.5、D > −5 が 0.03 → 0.99)。RF なし (RFSW OFF) だと放射損失で δ が周ごとに下がり (1024 周で約 −1%)、さらに意味がない。通常 NORAD'))
    if setup.fluc:
        w.append(('fluc', 'FLUC ON: 量子励起の乱数で振幅が確率的に動き、D は乱数の種に依る。同じ条件でも種で変わる (統計が必要)。通常 NOFLUC'))
    try:
        nux, nuy = setup.tunes_full()
        rc = resonance_check(nux, nuy, 4, maxn=6)
        for (m, n, p, d, o, nm) in rc:
            if o <= 2 and nm.startswith(('integer', 'half-integer')) and d < thr_tune:
                w.append(('tune_integer', '作動点 (%.4f, %.4f) が %s の共鳴 %d νx %+d νy = %d に近い (距離 %.4f < %.2f): 線形の運動が不安定に近く、D が大きくなる・NAFF の窓が足りない可能性' % (
                    nux, nuy, nm, m, n, p, d, thr_tune)))
            elif o == 2 and d < 0.02:
                w.append(('tune_lowres', '作動点 (%.4f, %.4f) が 2 次の%s共鳴 %d νx %+d νy = %d に近い (距離 %.4f < 0.02): 結合・6 極の 2 次の項で D が大きい領域になる可能性' % (
                    nux, nuy, '和・差の' if m * n else '', m, n, p, d)))
            elif o == 3 and d < 0.01:
                w.append(('tune_lowres', '作動点 (%.4f, %.4f) が 3 次の共鳴 %d νx %+d νy = %d に近い (距離 %.4f): 小振幅でも D が大きい領域になる可能性' % (nux, nuy, m, n, p, d)))
    except Exception:
        pass
    r = setup.info.get('orbit_resid')
    if r is not None and r == r and r > 1e-9:
        w.append(('orbit', '閉軌道が収束していない (残差 %.2g)' % r))
    if _has_errors(setup.els):
        w.append(('errors', '誤差つき (ずれ・ロールあり): 周期性がなく、結果は位置 %s の初期条件に固有 (他のセル・位置では違う)。閉軌道はずれ、D > −5 の割合は誤差なしより大きい (ASTRID2 で 0.03 → 0.5) のが普通。'
                  '誤差の種 (乱数) ごとにばらつくので、複数の種の統計で判断する' % (loc if loc is not None else '')))
        try:
            es = errors_summary(setup, loc=loc)
            if not es['stable']:
                w.append(('unstable', '誤差つきの線形の運動が不安定、または閉軌道が求まらない (閉軌道まわりの 1 周の行列の固有値が円の上にない)。周波数マップは意味を持たない'))
                return w
            d0 = (setup.tunes_full()[0] + es['dnu_pred'][0]) - (setup.tunes_full()[1] + es['dnu_pred'][1])
            d0 = abs(d0 - round(d0))
            if es['cminus_pred'] > 0.5 * max(d0, 1e-3):
                w.append(('coupling', '結合が強い: |c−| の見積もり %.4f が差の共鳴 νx − νy までの距離 %.4f の半分より大きい。νx・νy の対応が崩れ、D の下限が上がる' % (es['cminus_pred'], d0)))
            mo = max(es['orbit_max'])
            if mo > 0.2 * setup.alost:
                w.append(('orbit_large', '閉軌道のずれが大きい (最大 %.2f mm、LOSSAMPL %.1f mm の 20%% 超): 初期条件が閉軌道からの相対で、アパーチャが非対称になる' % (mo * 1e3, setup.alost * 1e3)))
        except Exception:
            pass
    return w


def check_result(res, nturn=None, thr_lost=0.2, thr_naff=0.0, thr_edge=0.02):
    """Return warnings about a frequency-map result (many lost points, all lost, NAFF failure, tunes near 0 or 0.5).
    周波数マップの結果の警告。コード: lost_many (失われた点が多い)、lost_all (全部失われた)、naff_fail (NAFF が収束しない・非有限)、tune_edge (チューンが 0・0.5 に近い点が多い)"""
    w = []
    n = len(res)
    nl = int((res.lost_turn > 0).sum())
    if nl == n:
        w.append(('lost_all', '全部の点 (%d) が失われた: D は求まらない。振幅・周回数・LOSSAMPL を確かめる' % n))
        return w
    if nl > thr_lost * n:
        w.append(('lost_many', '失われた点が多い (%d / %d = %.0f%%): D の統計 (D > −5 の割合など) は残った点だけの値で、実際はもっと悪い (失われた点を D 大とみなした割合は %.2f)。'
                  '作動点が共鳴に近い・アパーチャ・振幅の範囲を確かめる' % (nl, n, 100.0 * nl / n, (int((res.D > -5).sum()) + nl) / n)))
    nf = int(((res.lost_turn == 0) & ~np.isfinite(res.D)).sum())
    if nf > thr_naff * n:
        w.append(('naff_fail', 'NAFF が収束しない・非有限の点が %d (残った点 %d のうち): D は nan' % (nf, n - nl)))
    ok = np.isfinite(res.nux1) & np.isfinite(res.nuy1)
    if ok.any():
        def edge(v):
            return np.minimum(np.abs(v - np.round(v)), np.abs(v - np.round(2 * v) / 2.0)) < 0.01
        e = (edge(res.nux1[ok]) | edge(res.nuy1[ok])).mean()
        if e > thr_edge:
            w.append(('tune_edge', 'チューンの端数が 0・0.5 の 0.01 以内の点が %.0f%% (整数・半整数の共鳴の近く)。NAFF の分解能 (1/周回数) と同程度で D の信頼度が低い' % (100 * e)))
    return w
