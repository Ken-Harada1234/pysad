# -*- coding: utf-8 -*-
"""Local momentum acceptance (LMA) by tracking with radiation and RF, found by coarse scan and bisection.
局所運動量アクセプタンス (LMA: local momentum acceptance) (v2.3.0 Q1)。

位置 s (要素の入口) ごとに、そこで運動量偏差 δ を正・負に変えた粒子 (横方向は閉軌道、z は閉軌道) を放射・RF つきで追跡し、
nturn 周のあいだ失われない δ の限界 δ+(s)・δ−(s) を、粗い走査と二分法で求める。損失の判定は SAD の TrackParticles と同じ
(tapert: |x|,|y| < LOSSAMPL、|z| ≤ LOSSDZ、APERT 要素。trksetup.TrackSetup.track)。

  momentum_acceptance(setup, locs, ...)   → LMAResult
  LMAResult.table()/save()/plot()/touschek_input()/touschek_lifetime(...)

アルゴリズム (位置ごと、正負それぞれ):
  1. δ = 0 と、δ = dpmax·k/ncoarse (k = 1..ncoarse) を同時に追跡する (1 つの追跡の束)。外向きに最初に失われた点を挟む区間 (lo, hi) を作る
     (lo: 残った最大の δ、hi: 失われた最小の δ)。全部残れば δ = dpmax (走査の上限で打ち切り、'saturated')。
  2. 区間 (lo, hi) を nsub+1 等分して nsub 点を追跡し、最初に失われた点と、その手前の残った点で区間を狭める。hi − lo ≤ tol まで繰り返す
     (nsub = 1 なら普通の二分法)。結果は lo (確かに残った δ)。真の限界は (lo, hi) のどこか (幅 ≤ tol)。
  外向きの最初の損失を限界とする (走査の点のあいだの狭い損失の窓・その外側の残る島は見ない: 限界 → 保守的)。
RF バケット: bucket=True (既定。RFSW ON のとき) は、|z| がバケットの長さ (trksetup.bucket_zlimit) を超えた粒子も失われたとする
(バケットの外へ出た粒子は位相がすべる)。損失の理由は 'aper' (横の |x|,|y|、APERT) か 'rf' (|z|) で記録する。
周回数 (nturn の既定): シンクロトロン振動 nsynch (既定 2) 周期分 = ceil(nsynch/|νs|)。nsynch < 1 は許さない (nturn を直接指定したときを除く)。
並列: 位置ごとに独立なので、子プロセス (workers 個、既定 2、優先度を下げる) に分ける。
"""
import math
import os
import time

import numpy as np

from . import trksetup as _ts

_WSET = {}


def _init_worker(setup):
    import sys
    from . import compat
    compat.set_mode(setup.mode)
    try:
        compat.set_bend_k1_edge(setup.bend_k1_edge)
    except Exception:
        pass
    _ts.lower_priority()
    _WSET['s'] = setup


def _work(args):
    return _lma_location(_WSET['s'], *args)


def _first_lost(lost_turn):
    """外向きの順に並んだ δ の損失 (lost_turn > 0 が失われた) の最初の添字。なければ None"""
    w = np.nonzero(lost_turn > 0)[0]
    return int(w[0]) if len(w) else None


def _bisect(setup, loc, cod, nturn, nsub, tol, zb, x0, y0, signs, alost, st):
    """区間 (lo, hi) を hi − lo ≤ tol まで狭める (nsub 点ずつ)。追跡した粒子の数を返す"""
    npart = 0
    while True:
        pts = []
        for sg in signs:
            s = st[sg]
            if s['hi'] is None or s['hi'] - s['lo'] <= tol:
                continue
            for j in range(1, nsub + 1):
                pts.append((sg, s['lo'] + (s['hi'] - s['lo']) * j / (nsub + 1)))
        if not pts:
            break
        res = _track_points(setup, loc, cod, pts, nturn, zb, x0, y0, alost)
        npart += len(pts)
        for sg in signs:
            s = st[sg]
            if s['hi'] is None or s['hi'] - s['lo'] <= tol:
                continue
            mine = [p for p in pts if p[0] == sg]
            lt = np.array([res[p][0] for p in mine])
            k = _first_lost(lt)
            if k is None:
                s['lo'] = mine[-1][1]
            else:
                s['hi'] = mine[k][1]
                s['reason'] = res[mine[k]][1]
                s['lost_turn'] = res[mine[k]][0]
                if k > 0:
                    s['lo'] = mine[k - 1][1]
    return npart


def _lma_location(setup, loc, cod, nturn, dpmax, ncoarse, nsub, tol, zb, x0, y0, signs=(1, -1), alost=None, verify=0, turncheck=0):
    """1 つの位置の δ+・δ−。戻り dict (dp[±]、lo・hi、理由、追跡した粒子の数)。
    verify (点の数): 二分で見つけた限界 lo (全部残れば dpmax) の手前 [0, lo] を verify 点の等間隔で再追跡し、粗い走査の点のあいだの狭い損失の窓を探す。
    窓があれば最初に失われた点を新しい上の端として二分をやり直し (最大 5 回)、見つけた窓の点の数と最初の窓の δ を記録する。
    turncheck (倍率): 限界 lo より tol だけ手前の粒子を nturn × 倍率 周追跡し、失われれば 'tdep' (限界が周回数に依存し、増やすと tol より大きく下がる) とする"""
    t0 = time.time()
    st = {}
    npart = 0
    for sg in signs:
        st[sg] = dict(lo=0.0, hi=None, reason='', sat=False)
    # --- 粗い走査 (δ = 0 も入れる)
    grid = [0.0] + [dpmax * k / ncoarse for k in range(1, ncoarse + 1)]
    pts = []
    for sg in signs:
        for d in grid:
            if sg == signs[0] or d > 0.0:
                pts.append((sg, d))
    res = _track_points(setup, loc, cod, pts, nturn, zb, x0, y0, alost)
    npart += len(pts)
    zero_lost = False
    for sg in signs:
        d_list = [p for p in pts if p[0] == sg]
        rr = [res[(p[0], p[1])] for p in d_list]
        if sg != signs[0]:
            rr = [res[(signs[0], 0.0)]] + rr
            d_list = [(sg, 0.0)] + d_list
        lt = np.array([r[0] for r in rr])
        k = _first_lost(lt)
        s = st[sg]
        if k is None:
            s['lo'], s['hi'], s['sat'] = dpmax, None, True
        elif k == 0:
            s['lo'], s['hi'], s['reason'] = 0.0, 0.0, rr[0][1]
            s['zero'] = True
            zero_lost = True
        else:
            s['lo'], s['hi'], s['reason'] = d_list[k - 1][1], d_list[k][1], rr[k][1]
            s['lost_turn'] = rr[k][0]
    # --- 二分 (nsub 点ずつ)
    npart += _bisect(setup, loc, cod, nturn, nsub, tol, zb, x0, y0, signs, alost, st)
    # --- 窓の確認 (verify)
    for sg in signs:
        st[sg]['nwin'] = 0
        st[sg]['win'] = float('nan')
    if verify and verify > 0:
        for _it in range(5):
            pts = []
            for sg in signs:
                s = st[sg]
                top = dpmax if s['hi'] is None else s['lo']
                if top <= 0.0 or s.get('zero'):
                    continue
                pts.extend((sg, float(d)) for d in np.linspace(0.0, top, int(verify) + 2)[1:-1])
            if not pts:
                break
            res = _track_points(setup, loc, cod, pts, nturn, zb, x0, y0, alost)
            npart += len(pts)
            found = False
            for sg in signs:
                mine = [p for p in pts if p[0] == sg]
                if not mine:
                    continue
                lt = np.array([res[p][0] for p in mine])
                k = _first_lost(lt)
                if k is None:
                    continue
                s = st[sg]
                s['nwin'] += int((lt > 0).sum())
                if not (s['win'] == s['win']):
                    s['win'] = mine[k][1]
                s['lo'] = mine[k - 1][1] if k > 0 else 0.0
                s['hi'] = mine[k][1]
                s['sat'] = False
                s['reason'] = res[mine[k]][1]
                s['lost_turn'] = res[mine[k]][0]
                found = True
            if not found:
                break
            npart += _bisect(setup, loc, cod, nturn, nsub, tol, zb, x0, y0, signs, alost, st)
    # --- 周回数の確認 (turncheck)
    if turncheck and turncheck > 0:
        nt2 = int(math.ceil(nturn * float(turncheck)))
        pts = [(sg, max(st[sg]['lo'] - tol, 0.0)) for sg in signs]      # 限界より tol だけ手前 (tol 以内の変化は周回数の依存とみなさない)
        res = _track_points(setup, loc, cod, pts, nt2, zb, x0, y0, alost)
        npart += len(pts)
        for sg in signs:
            st[sg]['tdep'] = bool(res[(sg, max(st[sg]['lo'] - tol, 0.0))][0] > 0)
    out = dict(loc=loc, npart=npart, seconds=time.time() - t0)
    for sg, nm in ((1, 'plus'), (-1, 'minus')):
        if sg in st:
            s = st[sg]
            out['dp_' + nm] = s['lo']
            out['hi_' + nm] = s['hi'] if s['hi'] is not None else float('nan')
            out['reason_' + nm] = 'sat' if s['sat'] else s['reason']
            out['nwin_' + nm] = s['nwin']
            out['win_' + nm] = s['win']
            out['tdep_' + nm] = s.get('tdep', False)
            out['lostturn_' + nm] = s.get('lost_turn', 0)
    return out


def _track_points(setup, loc, cod, pts, nturn, zb, x0, y0, alost):
    c = np.tile(cod, (len(pts), 1))
    for i, (sg, d) in enumerate(pts):
        c[i, 5] += sg * d
        c[i, 0] += x0
        c[i, 2] += y0
    r = setup.track(loc, c, nturn, zlost=zb, alost=alost)
    out = {}
    for i, p in enumerate(pts):
        lt = int(r['lost_turn'][i])
        why = ''
        if lt > 0:
            f = r['final'][i]
            why = 'rf' if (zb is not None and abs(f[4] - cod[4]) > zb * 0.999) else 'aper'
        out[p] = (lt, why)
    return out


class LMAResult:
    """Result of `momentum_acceptance`: arrays of location, s, name, delta+ and delta-, loss reasons and the conditions used.
    momentum_acceptance の結果。配列 (位置の順): loc、s、name、dp_plus (δ+ ≥ 0)、dp_minus (δ− ≤ 0 の値を負で)、
    hi_plus・hi_minus (最初に失われた δ の絶対値。nan は上限まで残った)、reason_plus・reason_minus ('aper'・'rf'・'sat')、
    cond (使った条件の辞書)"""

    def __init__(self, loc, s, name, dpp, dpm, hip, him, rp, rm, cond, setup=None, extra=None, warnings=None):
        self.loc = np.asarray(loc, dtype=int)
        self.s = np.asarray(s, dtype=float)
        self.name = list(name)
        self.dp_plus = np.asarray(dpp, dtype=float)
        self.dp_minus = -np.asarray(dpm, dtype=float)       # 負で持つ
        self.hi_plus = np.asarray(hip, dtype=float)
        self.hi_minus = -np.asarray(him, dtype=float)
        self.reason_plus = list(rp)
        self.reason_minus = list(rm)
        self.cond = dict(cond)
        self.setup = setup
        ex = dict(extra or {})
        n = len(self.loc)
        z = lambda v: np.asarray(ex.get(v, np.zeros(n)))
        self.nwin_plus, self.nwin_minus = z('nwin_plus').astype(int), z('nwin_minus').astype(int)         # verify で見つけた窓の点の数
        self.win_plus, self.win_minus = np.asarray(ex.get('win_plus', np.full(n, np.nan)), dtype=float), np.asarray(ex.get('win_minus', np.full(n, np.nan)), dtype=float)
        self.tdep_plus, self.tdep_minus = z('tdep_plus').astype(bool), z('tdep_minus').astype(bool)       # turncheck: 周回数を増やすと限界が下がる
        self.lostturn_plus, self.lostturn_minus = z('lostturn_plus').astype(int), z('lostturn_minus').astype(int)   # 限界のすぐ上の粒子が失われた周
        self.warnings = list(warnings or [])

    def flags(self, i):
        """Return the caution flags for location i (W, T, S and 0) that mark possibly unreliable acceptance limits.
        位置 i の注意の旗: W+/W− (verify で窓を見つけた)、T+/T− (周回数を 2 倍にすると、限界より tol 手前の粒子が失われる = 限界が tol より大きく下がる)、S+/S− (走査の上限まで残った)、0 (δ = 0 で失われた)"""
        f = []
        for sg, nm, rs, dp in (('+', 'plus', self.reason_plus, self.dp_plus), ('-', 'minus', self.reason_minus, self.dp_minus)):
            if getattr(self, 'nwin_' + nm)[i] > 0:
                f.append('W' + sg)
            if getattr(self, 'tdep_' + nm)[i]:
                f.append('T' + sg)
            if rs[i] == 'sat':
                f.append('S' + sg)
            if dp[i] == 0.0:
                f.append('0' + sg)
        return ','.join(f) if f else '-'

    def __len__(self):
        return len(self.loc)

    def table(self):
        """Return a text table of loc, s, name, delta+, delta- and loss reasons, with the conditions in the header.
        列 (loc, s, name, dp_plus, dp_minus, reason_plus, reason_minus) の文字列表 (ヘッダーに条件)"""
        lines = ['# LMA: local momentum acceptance (dp_plus >= 0, dp_minus <= 0). bracket width <= tol, value = largest surviving dp']
        for k, v in self.cond.items():
            lines.append('# %s = %s' % (k, v))
        for w in self.warnings:
            lines.append('# WARNING: ' + w)
        lines.append('# loc  s[m]  name  dp_plus  dp_minus  hi_plus  hi_minus  reason_plus  reason_minus  flags')
        for i in range(len(self)):
            lines.append('%5d %12.6f %-10s %12.6e %12.6e %12.6e %12.6e %-5s %-5s %s' % (
                self.loc[i], self.s[i], self.name[i], self.dp_plus[i], self.dp_minus[i], self.hi_plus[i], self.hi_minus[i],
                self.reason_plus[i], self.reason_minus[i], self.flags(i)))
        return '\n'.join(lines) + '\n'

    def save(self, path):
        with open(path, 'w', encoding='utf-8') as f:
            f.write(self.table())

    def touschek_input(self):
        """Return the acceptance as (s, delta+, delta-) in the form accepted by `touschek_lifetime_local`.
        既存の Touschek の機能 (touschek_lifetime_local) に渡す形: (s, δ+, δ−) (δ− は負)"""
        return self.s.copy(), self.dp_plus.copy(), self.dp_minus.copy()

    def touschek_lifetime(self, table, ssel=None, period=None, s0=0.0, sig_by_cell=None, **kw):
        """Call `touschek_lifetime_local` with the interpolation of this result, periodic if a period is given.
        touschek_lifetime_local(table, 全要素の s, ...) を、この結果の補間 (period があれば周期的) で呼ぶ"""
        return touschek_lifetime_local(table, self, setup=self.setup, period=period, s0=s0, **kw)

    def plot(self, path, setup=None, title=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        _jfont(plt)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(self.s, self.dp_plus * 100, '.-', ms=3, lw=0.8, label='dp+ (positive)')
        ax.plot(self.s, self.dp_minus * 100, '.-', ms=3, lw=0.8, label='dp- (negative)')
        b = self.cond.get('bucket_height')
        if b:
            ax.axhline(b * 100, color='gray', ls='--', lw=0.8, label='RF bucket height')
            ax.axhline(-b * 100, color='gray', ls='--', lw=0.8)
        ax.set_xlabel('s [m]')
        ax.set_ylabel('local momentum acceptance [%]')
        ax.set_title(title or 'LMA')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(path, dpi=130)
        plt.close(fig)


def _jfont(plt):
    try:
        from matplotlib import font_manager as fm
        have = {f.name for f in fm.fontManager.ttflist}
        for nm in ('Yu Gothic', 'Meiryo', 'MS Gothic', 'IPAexGothic', 'Noto Sans CJK JP'):
            if nm in have:
                plt.rcParams['font.family'] = nm
                break
        plt.rcParams['axes.unicode_minus'] = False
    except Exception:
        pass


def momentum_acceptance(setup, locs=None, nturn=None, nsynch=2.0, dpmax=None, ncoarse=12, nsub=3, tol=1e-4, bucket=True,
                        bucket_margin=1.0, x0=0.0, y0=0.0, workers=2, progress=0, alost=None, signs=(1, -1), check_every=None,
                        verify=0, turncheck=2.0, quiet=False):
    """Compute the local momentum acceptance at the given locations of a `TrackSetup` by scan and bisection.
    局所運動量アクセプタンスを求める。
    setup: trksetup.TrackSetup。locs: 位置 (要素の入口、1 始まり) のリスト (None なら要素の境界すべて: select_locations(setup))。
    nturn: 周回数 (None なら nsynch 周期分: ceil(nsynch/|νs|))。RF なし (RFSW OFF) のときは nturn を指定する (既定 512)。
    dpmax: 走査の上限 (None なら RF ありは バケット高さの 1.25 倍、なしは 0.05)。ncoarse: 粗い走査の点数。nsub: 二分の 1 回の追跡点数 (1 は二分法)。
    tol: 区間の幅の終了条件 (δ の絶対値)。bucket: RF バケットから出た粒子 (|z| > バケットの長さ) を失われたとする。
    x0・y0: 初期の x・y のずれ (閉軌道から。既定 0)。workers: 並列の子プロセスの数 (既定 2。1 なら親プロセスで順に)。
    progress: 0 以外なら、その数の位置ごとに経過を出す。alost: |x|,|y| の限界 (None なら setup.alost = LOSSAMPL)。
    check_every: 損失の判定の間隔 (要素の数)。**既定 (None) は 1 = 要素ごとの判定** (setup.check_every が 0 でないときはその値)。0 は SAD の tapert のまま (最初の要素・
    15 要素ごと・区間の終わり、gs の TrackParticles と同じ判定)。アパーチャが限界のとき、SAD の判定は判定点のあいだのずれの最大を見逃し、δ の限界を大きく見積もる
    (ASTRID2 の例で最大 1.5 倍) ので、既定は安全側の 1。0 を選ぶと警告が出る (アパーチャが限界の位置があるとき)。
    verify: 0 でなければ、限界の手前 [0, δ] をその数の点で再追跡し、粗い走査のあいだの狭い損失の窓を探す (窓があれば限界を下げる。LMAResult.nwin_*)。
    turncheck: 0 でなければ、限界の粒子をその倍率の周回数 (既定 2 倍) 追跡し、失われたら限界が周回数に依存する (tdep_*、警告)。quiet: True なら警告を出さない (cond['warnings'] には残る)."""
    explicit_ce_arg = check_every
    if check_every is None:
        check_every = setup.check_every if setup.check_every else 1
    setup = setup.copy(check_every=check_every)
    if locs is None:
        locs = _ts.select_locations(setup)
    locs = [int(l) for l in locs]
    nt_syn = setup.nturn_for_synchrotron(nsynch)
    nturn_given = nturn is not None
    explicit_ce = None if explicit_ce_arg is None else int(explicit_ce_arg)
    if nturn is None:
        nturn = nt_syn if nt_syn is not None else 512
    bh = setup.bucket.get('height', 0.0) if setup.rfsw else 0.0
    if dpmax is None:
        dpmax = 1.25 * bh if bh > 0.0 else 0.05
    zb = None
    if bucket and setup.rfsw:
        zb = setup.bucket_zlimit(bucket_margin)
    if zb is not None:
        zb = min(zb, setup.zlost)
    cods = setup.closed_orbit(locs)
    args = [(l, cods[i], int(nturn), float(dpmax), int(ncoarse), int(nsub), float(tol), zb, float(x0), float(y0), tuple(signs), alost,
             int(verify or 0), float(turncheck or 0))
            for i, l in enumerate(locs)]
    t0 = time.time()
    results = []
    workers = max(1, min(int(workers), max(1, (os.cpu_count() or 2) - 1), len(args)))
    if workers == 1:
        for i, a in enumerate(args):
            results.append(_lma_location(setup, *a))
            if progress and (i + 1) % progress == 0:
                print('    [LMA] %d/%d (%.0f s)' % (i + 1, len(args), time.time() - t0), flush=True)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(setup,)) as ex:
            for i, r in enumerate(ex.map(_work, args, chunksize=1)):
                results.append(r)
                if progress and (i + 1) % progress == 0:
                    print('    [LMA] %d/%d (%.0f s)' % (i + 1, len(args), time.time() - t0), flush=True)
    g = lambda k, d=float('nan'): [r.get(k, d) for r in results]
    cond = setup.conditions()
    cond.update(nturn=int(nturn), nsynch_default=nsynch, dpmax=float(dpmax), ncoarse=ncoarse, nsub=nsub, tol=tol, bucket=bool(zb is not None),
                bucket_zlimit_m=zb, x0=x0, y0=y0, npositions=len(locs), workers=workers, seconds=round(time.time() - t0, 1),
                tracked_particles=int(sum(r['npart'] for r in results)), alost_used=setup.alost if alost is None else alost,
                verify=int(verify or 0), turncheck=float(turncheck or 0))
    extra = {k: g(k, 0) for k in ('nwin_plus', 'nwin_minus', 'tdep_plus', 'tdep_minus', 'lostturn_plus', 'lostturn_minus')}
    extra.update(win_plus=g('win_plus'), win_minus=g('win_minus'))
    res = LMAResult(locs, [setup.s_of(l) for l in locs], [setup.names[l - 1] for l in locs],
                    g('dp_plus', 0.0), g('dp_minus', 0.0), g('hi_plus'), g('hi_minus'),
                    g('reason_plus', ''), g('reason_minus', ''), cond, setup=setup, extra=extra)
    res.warnings = lma_warnings(setup, res, check_every_explicit=explicit_ce, alost_arg=alost, nturn_given=nturn_given)
    res.cond['warnings'] = list(res.warnings)
    if not quiet:
        import warnings as _w
        for m in res.warnings:
            _w.warn('LMA: ' + m, stacklevel=2)
    return res


def lma_warnings(setup, res, check_every_explicit=None, alost_arg=None, nturn_given=False):
    """Build warning strings for unreliable or unusable conditions from the result and conditions.
    結果と条件から、信頼できない・使えない条件の警告 (文字列のリスト) を作る (momentum_acceptance が呼ぶ)"""
    w = []
    c = res.cond
    n = len(res)
    al = setup.alost if alost_arg is None else alost_arg
    if al >= 0.1:
        w.append('LOSSAMPL = %.3g m: no physical aperture is set; the limit is the RF bucket (or the %.3g m cut), not an aperture limit' % (al, al))
    if check_every_explicit == 0 and n and any(r == 'aper' for r in list(res.reason_plus) + list(res.reason_minus)):
        w.append('check_every = 0 (SAD tapert, every 15 elements): at aperture-limited positions the limit is overestimated (ASTRID2: up to 1.5 times); '
                 'use check_every = 1 (the default)')
    if setup.fluc:
        w.append('FLUC is ON: the bisection is stochastic; the result is one sample (repeat with several seeds, use the mean and spread)')
    if not setup.rfsw:
        w.append('RFSW is OFF: no RF bucket, nturn is not tied to the synchrotron period; the limit depends on nturn only through betatron/chromatic motion')
    else:
        nsyn = setup.nturn_for_synchrotron(1.0, nmin=1)
        if nsyn and c.get('nturn', 0) < nsyn:
            w.append('nturn = %d is shorter than one synchrotron period (%d turns)' % (c['nturn'], nsyn))
    nsat = int(sum(1 for r in list(res.reason_plus) + list(res.reason_minus) if r == 'sat'))
    if nsat:
        w.append('%d limits reached the scan top dpmax = %.4g without a loss (the acceptance is larger than dpmax)' % (nsat, c.get('dpmax', float('nan'))))
    nz = int(sum(1 for v in list(res.dp_plus) + list(-res.dp_minus) if v == 0.0))
    if nz:
        w.append('%d limits are 0: the particle on the closed orbit itself is lost (aperture too small or unstable orbit)' % nz)
    nt = int(res.tdep_plus.sum() + res.tdep_minus.sum())
    if nt:
        w.append('%d limits are turn dependent (a particle tol below the limit is lost when nturn is multiplied by %.3g): '
                 'the limit falls by more than tol with more turns (near the RF bucket edge, nturn below ~2 synchrotron periods)' % (nt, c.get('turncheck', 0)))
    nw = int(res.nwin_plus.sum() + res.nwin_minus.sum())
    if nw:
        w.append('verify found %d lost points below the bisection limit (narrow loss windows); the limit was lowered to the first window' % nw)
    r = c.get('orbit_residual')
    if r is not None and r == r and r > 1e-9:
        w.append('closed orbit residual %.2e (> 1e-9): the closed orbit is not converged' % r)
    try:
        from .trksetup import find_periodicity
        pass
    except Exception:
        pass
    return w


def periodicity_report(setup, cod_tol=None, warn=True):
    """Check whether one-cell results may be extended to the full ring, adding a closed-orbit comparison.
    1 セルの結果を全周に延ばしてよいかの確認 (v2.3.0 L2)。find_periodicity の結果に、閉軌道 (x, px, y, py) が各セルの最初の位置で同じかを足す。
    閉軌道の許容 cod_tol (None なら 1e-7 + 2 U0/E0 max(|ηx|, |η'x|, |ηy|, |η'y|)): RF が 1 か所なので、閉軌道の δ がセルごとに U0/E0 ずつ変わり、分散の分だけ
    x が変わる (正常な周期的なリングでもこの大きさの差が出る) ことを許す。
    戻り dict: nper、starts、cell_length、s_starts、ok (セルの長さが同じ)、periodic (nper > 1 かつ ok かつ閉軌道が同じ)、cod_spread (セルの最初の位置の閉軌道の
    x, px, y, py の最大の差)、cod_tol (使った許容)、delta_sawtooth (U0/E0 = 周期性が崩れる RF の位置による δ の差の見積もり: 1 セルの結果と別のセルの差の目安)。
    誤差・補正・空洞が 2 か所以上などで周期的でないとき、warn=True なら警告を出す (1 セルの結果は全周に使えない: 全周の位置で計算する)"""
    from .trksetup import find_periodicity
    per = find_periodicity(setup.els, setup.rows)
    rep = dict(per)
    rep['cod_spread'] = float('nan')
    b = setup.bucket
    rep['delta_sawtooth'] = (b['U0'] / setup.p0) if b and b.get('U0') else 0.0
    msgs = []
    if cod_tol is None:
        eta = max(max(abs(r[5]), abs(r[6])) for r in setup.rows)
        cod_tol = 1e-7 + 2.0 * rep['delta_sawtooth'] * eta
    if per['nper'] <= 1:
        msgs.append('the element sequence is not periodic (errors, an insertion device, a different cell, ...): the one-cell LMA cannot be extended to the ring; compute the whole ring')
    elif not per.get('ok', True):
        msgs.append('the cells have different lengths: not periodic')
    else:
        try:
            co = setup.closed_orbit(per['starts'])[:, :4]
            rep['cod_spread'] = float(np.max(np.abs(co - co[0])))
            if rep['cod_spread'] > cod_tol:
                msgs.append('the closed orbit at the cell starts differs by %.2e (> %.0e): the ring is not periodic in the orbit (misalignments?); '
                            'the one-cell LMA cannot be extended' % (rep['cod_spread'], cod_tol))
        except Exception as e:      # 閉軌道が求まらない
            msgs.append('closed orbit at the cell starts failed (%s)' % e)
    rep['cod_tol'] = cod_tol
    rep['periodic'] = not msgs
    rep['messages'] = msgs
    if warn:
        import warnings as _w
        for m in msgs:
            _w.warn('LMA periodicity: ' + m, stacklevel=2)
    return rep


def cell_locations(setup, ds=0.0, step=1, warn=True):
    """Return the element-boundary locations of one cell if periodic, otherwise warn and return the full-ring locations, with a report.
    1 セル (周期的なとき) の要素の境界の位置。周期的でなければ警告して全周の位置を返す。戻り (locs, report)。report['periodic'] が真なら、結果を
    `period = report['cell_length']`、`s0 = report['s_starts'][0]` で全周に延ばせる (touschek_lifetime_local の period・s0)。
    1 セルの結果と別のセルの同じ位置の差の目安は report['delta_sawtooth'] (= U0/E0。RF の位置が 1 か所なので、閉軌道の δ がセルごとに U0/E0 ずつ変わる)"""
    from .trksetup import select_locations
    rep = periodicity_report(setup, warn=warn)
    if rep['periodic']:
        locs = select_locations(setup, rep['starts'][0], rep['starts'][1] - 1, step, ds)
    else:
        locs = select_locations(setup, 1, None, step, ds)
    return locs, rep


# ---------------------------------------------------------------------------------------------------------------------
# 既存の Touschek (touschek.py: TouschekTable) に渡す
# ---------------------------------------------------------------------------------------------------------------------

def acceptance_at(s_query, s, dpp, dpm, period=None, s0=0.0):
    """Linearly interpolate (delta+, |delta-|) at position s_query from the LMA points, periodic if a period is given.
    位置 s_query での (δ+, |δ−|) を、LMA の点 (s, δ+, δ−) から線形補間する。period (セルの長さ) を渡すと周期的 (s0 から period の範囲を
    1 セルとして、s の外は折り返す)。period なしは範囲の外で端の値"""
    s = np.asarray(s, dtype=float)
    q = np.asarray(s_query, dtype=float)
    o = np.argsort(s)
    s, dpp, dpm = s[o], np.asarray(dpp)[o], np.abs(np.asarray(dpm))[o]
    if period:
        q = s0 + np.mod(q - s0, period)
        ss = np.concatenate([s - period, s, s + period])
        return np.interp(q, ss, np.tile(dpp, 3)), np.interp(q, ss, np.tile(dpm, 3))
    return np.interp(q, s, dpp), np.interp(q, s, dpm)


def touschek_lifetime_local(table, lma, setup=None, period=None, s0=0.0, s_el=None, dp_floor=None):
    """Compute the Touschek lifetime using position-dependent momentum acceptances delta+(s) and delta-(s).
    位置ごとの運動量アクセプタンス (δ+(s), δ−(s)) を使った Touschek 寿命。
    table: touschek.touschek_table()['table'] (または ffs.touschek['table']): table[3] は要素ごと (nlat × 120) の損失率、table[1][0] は
    一様な運動量アクセプタンス dp (0.004 … 0.242、0.002 刻み) の表。SAD の TouschekTable は、各要素 l の損失率 R_l(a) を
    一様な受容 a (正負対称) の関数として持つ。ここでは要素 l に正負の受容 δ+(l)・|δ−(l)| がある場合を、
        1/τ = Σ_l ½ [R_l(δ+_l) + R_l(|δ−_l|)]
    として求める (対称な受容 a での損失は、散乱された 2 粒子のうち +Δ を受け取る粒子と −Δ を受け取る粒子の両方の損失を含むので、片側ずつ半分)。
    R_l(a) は dp の表の上で log-log 補間 (表の外は a^−2 で外挿。a が表の最小 0.004 より小さいと過小評価でなく、外挿の値)。
    δ±(l): 要素 l の中央の s での lma の線形補間 (period があれば周期的)。s_el: 要素の s (None なら setup.rows から要素の中央)。
    戻り dict: rate (1/s)、lifetime (s)、rate_per_element、dp_plus・dp_minus (要素ごと)、extrapolated (表の外を使った要素の数)。
    対称で一様な a のとき (δ+ = |δ−| = a が表の点) は table[1][1] の値と一致する"""
    dp_tab = np.asarray(table[1][0], dtype=float)
    R = np.asarray(table[3], dtype=float)          # (nlat, 120)
    nl = R.shape[0]
    if s_el is None:
        rows = setup.rows if setup is not None else lma.setup.rows
        s_all = np.array([r[0] for r in rows])
        s_el = np.empty(nl)
        for i in range(nl):
            s_el[i] = 0.5 * (s_all[i] + s_all[i + 1]) if i + 1 < len(s_all) else s_all[i]
    ap, am = acceptance_at(s_el, lma.s, lma.dp_plus, lma.dp_minus, period=period, s0=s0)
    if dp_floor:
        ap, am = np.maximum(ap, dp_floor), np.maximum(am, dp_floor)
    ldp = np.log(dp_tab)
    nextra = 0
    tot = np.zeros(nl)
    for sgn, a in ((0, ap), (1, am)):
        la = np.log(np.maximum(a, 1e-12))
        for l in range(nl):
            r = R[l]
            if not np.any(r > 0):
                continue
            lr = np.log(np.maximum(r, 1e-300))
            if la[l] > ldp[-1]:
                v = lr[-1] - 2.0 * (la[l] - ldp[-1])
                nextra += 1
            elif la[l] < ldp[0]:
                v = lr[0] - 2.0 * (la[l] - ldp[0])
                nextra += 1
            else:
                v = np.interp(la[l], ldp, lr)
            tot[l] += 0.5 * math.exp(v)
    rate = float(tot.sum())
    return dict(rate=rate, lifetime=(1.0 / rate if rate > 0 else float('inf')), rate_per_element=tot, dp_plus=ap, dp_minus=am,
                extrapolated=nextra, s_el=s_el)


def touschek_lifetime_local_piwinski(tt, lma, setup=None, period=None, s0=0.0, s_el=None, dp_floor=None, npan=40, bunch_delta=None, bunch_floor=0.25, bunch_qs=None):
    """Compute the local Touschek lifetime like `touschek_lifetime_local`, but with the Piwinski-type loss rate instead of the K64 table.
    touschek_lifetime_local と同じ受け渡し (位置ごとの δ+(s)・δ−(s)) で、損失率に K64 の表ではなく Piwinski 型の式 (touschek.touschek_piwinski、v2.3.0 T2) を使う。
    tt: touschek_table(..., record_rest=True) の戻り (FFS の INTRA ON の Emittance[] の結果 ffs.touschek)。
    要素ごとの受容は touschek_lifetime_local と同じ (要素の中央の s での線形補間、period があれば周期的)。戻り dict: rate・lifetime・dp_plus・dp_minus・s_el。
    対称な一様 a (δ+ = |δ−| = a) では touschek_piwinski(tt, a=a) と一致する。
    bunch_delta=True (v2.4.0 A3): バンチ内の δ の分布 (対の重心の δ_c ~ N(0, σδ²/2)) を入れた損失率 (touschek_bunch.touschek_piwinski_bunch。損失率は増え、寿命は短くなる)。
    **v2.4.0 F2 (B200): 既定 (None) は v2 = True (δ_old = 0 の基準は近似で、バンチ内の δ の分布を入れる方が正しい)、k64 = False。bunch_qs も同じ (相対運動量の縦成分 q_s の効果。数十秒)。False で切る**"""
    from .touschek import touschek_piwinski
    from . import compat as _compat
    if bunch_delta is None:
        bunch_delta = _compat.v2()
    if bunch_qs is None:
        bunch_qs = _compat.v2()
    nl = int(tt['rest_info']['nlat'])
    if s_el is None:
        rows = setup.rows if setup is not None else lma.setup.rows
        s_all = np.array([r[0] for r in rows])
        s_el = np.empty(nl)
        for i in range(nl):
            s_el[i] = 0.5 * (s_all[i] + s_all[i + 1]) if i + 1 < len(s_all) else s_all[i]
    ap, am = acceptance_at(s_el, lma.s, lma.dp_plus, lma.dp_minus, period=period, s0=s0)
    if dp_floor:
        ap, am = np.maximum(ap, dp_floor), np.maximum(am, dp_floor)
    if bunch_delta:
        # v2.4.0 A3: バンチ内の δ の分布 (対の重心の δ_c ~ N(0, σδ²/2)) を入れる。v2.4.0 F2: v2 の既定は入れる (touschek_bunch.py の docstring)
        from .touschek_bunch import touschek_piwinski_bunch
        r = touschek_piwinski_bunch(tt, ap=ap, am=am, npan=npan, floor=bunch_floor, qs=bunch_qs)
    else:
        r = touschek_piwinski(tt, ap=ap, am=am, npan=npan)
    r.update(dp_plus=ap, dp_minus=am, s_el=s_el)
    return r


# ---------------------------------------------------------------------------------------------------------------------
# 線形の見積もり (検証用)
# ---------------------------------------------------------------------------------------------------------------------

def linear_momentum_limit(setup, locs, aper=None):
    """Estimate the delta limit from the horizontal aperture in linear optics, for verification.
    線形光学での、水平の開口 ±aper による δ の限界の見積もり (検証用)。位置 s で運動量が δ ずれた粒子は、閉軌道 (x=x_co, px=px_co) から出発すると、
    δ の閉軌道 (η δ, η' δ) のまわりに振幅 √(H(s) β(s')) δ のベータトロン振動をする (H = γη² + 2αηη' + βη'² は分散の不変量)。
    位置 s' の |x| の最大は δ (|η(s')| + √(H(s) β(s')))。これが aper を超える最小の |δ|:
        δ_lin(s) = aper / max_{s'} (|η(s')| + √(H(s) β(s')))
    s' は Twiss の全位置。垂直方向の分散・結合・非線形 (6 極のクロマティック・振幅依存) は入れない。
    戻り: δ_lin (位置ごと)。aper=None なら setup.alost。損失の判定は要素の境界の tapert (最初の要素・15 要素ごと・周末) で行われるので、
    判定点のあいだの最大は見えない: 追跡の限界は、この見積もりより少し大きく出ることがある"""
    a = setup.alost if aper is None else float(aper)
    rows = np.array(setup.rows)
    s_, bx, ax, by, ay, ex, epx = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4], rows[:, 5], rows[:, 6]
    out = []
    for l in locs:
        b, al, e, ep = bx[l - 1], ax[l - 1], ex[l - 1], epx[l - 1]
        H = (1.0 + al * al) / b * e * e + 2.0 * al * e * ep + b * ep * ep
        amp = np.abs(ex) + np.sqrt(max(H, 0.0) * bx)
        out.append(a / amp.max())
    return np.array(out)
