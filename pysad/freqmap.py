# -*- coding: utf-8 -*-
"""周波数マップ (frequency map、FMA) (v2.3.0 Q1)。

1 つの位置 loc (要素の入口) から、初期条件の格子 (x, y) (または x–δ、y–δ) を追跡し、NAFF (pysad/naff.py) で、追跡の前半・後半のチューン
(νx, νy) を求め、チューンの拡散の指標
    D = log10 sqrt( (νx2 − νx1)² + (νy2 − νy1)² )
を出す (νi1: 前半 nturn/2 周、νi2: 後半 nturn/2 周)。途中で粒子が失われた点 (TrackParticles と同じ損失判定) は区別して持つ (D = nan、lost_turn > 0)。

  frequency_map(setup, loc, x=…, y=…, dp=…, nturn=1024)  → FMapResult
  FMapResult.table()/save()/plot(path)                   表 (ファイル・文字列)、図 (チューン平面・初期条件平面)
  resonance_lines(order, window)                          共鳴線 m νx + n νy = p のリスト
  resonance_check(nux, nuy, order)                        作動点に近い共鳴線 (整数・半整数を含む) の表

初期条件: x, y, dp は数または 1 次元の配列。配列が 2 つなら、その 2 つの平面の格子 (3 つなら 3 次元の格子)。位置 loc の閉軌道からの
ずれ (px = py = 0、z = 0)。δ ≠ 0 の点は、RFSW OFF なら δ の閉軌道 (4 次元、Newton 法) の上に、RFSW ON なら線形の分散 (η, η') の上に置く
(x, y は δ の軌道からのずれ)。チューンの単位は 1/周回、端数 (0〜1)。符号は SAD の位相の向き (naff.tune_xp)。
"""
import math
import os
import time
import warnings

import numpy as np

from . import naff as _naff
from . import trksetup as _ts


class FMapWarning(UserWarning):
    """周波数マップを信頼して使えない・信頼度が低い条件の警告 (v2.3.0 S2)。warnings.simplefilter('ignore', FMapWarning) で消せる。結果の FMapResult.warnings にも入る"""

_WSET = {}


def _init_worker(setup):
    from . import compat
    compat.set_mode(setup.mode)
    try:
        compat.set_bend_k1_edge(setup.bend_k1_edge)
    except Exception:
        pass
    _ts.lower_priority()
    _WSET['s'] = setup


def _cs_params(x, p):
    """標本の 2 次のモーメントから Courant–Snyder のパラメータ (α, β) (平均を引いた (x, p) の楕円)。振幅 0 なら None"""
    x = x - x.mean()
    p = p - p.mean()
    sxx, spp, sxp = np.mean(x * x), np.mean(p * p), np.mean(x * p)
    eps2 = sxx * spp - sxp * sxp
    if not (eps2 > 0.0):
        return None
    eps = math.sqrt(eps2)
    return -sxp / eps, sxx / eps


def _tunes_of(series, tw, window, norm='data'):
    """series: (nturn, 4) の (x, px, y, py)。前半・後半のチューン [νx1, νy1, νx2, νy2] (端数 0〜1) と振幅 (前半の x・y、後半の x・y の |a|)。
    norm: 'data' (既定) は、(x, px) の楕円の規格化 (α, β) を、その半分の追跡の標本の 2 次のモーメントから求める (δ ≠ 0 の軌道の光学関数は Twiss の
    δ = 0 の値と違い、規格化がずれると NAFF の周波数に 1e-6 程度の系統誤差が入るため)。'twiss' は引数 tw = ((αx, βx), (αy, βy)) を使う"""
    n = series.shape[0]
    h = n // 2
    (ax, bx), (ay, by) = tw
    out = []
    amp = []
    for (a, b) in ((0, h), (h, 2 * h)):
        for (ix, al, be) in ((0, ax, bx), (2, ay, by)):
            if norm == 'data':
                cs = _cs_params(series[a:b, ix], series[a:b, ix + 1])
                if cs is not None:
                    al, be = cs
            f, am = _naff.tune_xp(series[a:b, ix], series[a:b, ix + 1], al, be, window=window)
            out.append(f % 1.0)
            amp.append(am)
    return out, amp


def _fmap_chunk(args):
    setup, loc, coords, nturn, window, tw, norm = args[:7]
    scaling = args[7] if len(args) > 7 else False
    r = setup.track(loc, coords, nturn, record=True)
    n = len(coords)
    res = np.full((n, 9), np.nan)
    h = nturn // 2
    for i in range(n):
        lt = r['lost_turn'][i]
        if lt > 0:
            continue
        sr = r['rec'][:, i, :4]
        if not np.all(np.isfinite(sr)):
            continue
        try:
            f, a = _tunes_of(sr, tw, window, norm)
        except Exception:
            continue
        res[i, :8] = f + a
        if scaling and nturn >= 256:
            try:
                fh, _ah = _tunes_of(sr[:h], tw, window, norm)           # 最初の nturn/2 周 (前半・後半 nturn/4 周) の拡散 (D の周回数への依存、スケーリング)
                res[i, 8] = math.hypot(((fh[2] - fh[0] + 0.5) % 1.0) - 0.5, ((fh[3] - fh[1] + 0.5) % 1.0) - 0.5)
            except Exception:
                pass
    return r['lost_turn'], res, r['final']


def _chunk_job(args):
    return _fmap_chunk((_WSET['s'],) + tuple(args))


class FMapResult:
    """frequency_map の結果。配列 (点の数 n): x, y, dp (閉軌道からの初期のずれ)、nux1, nuy1, nux2, nuy2 (端数のチューン)、D (log10 の拡散)、
    lost_turn (0: 最後まで残った、k>0: k 周目に失われた)、amp (前半の x・y の振幅の規格化した値)。grid: (nx, ny, ndp)、cond: 使った条件"""

    def __init__(self, x, y, dp, nu, lost, grid, cond, setup=None, loc=None):
        self.x, self.y, self.dp = x, y, dp
        self.nux1, self.nuy1, self.nux2, self.nuy2 = nu[:, 0], nu[:, 1], nu[:, 2], nu[:, 3]
        self.amp = nu[:, 4:8]
        self.lost_turn = lost
        with np.errstate(all='ignore'):
            d = np.hypot(self.nux2 - self.nux1, self.nuy2 - self.nuy1)
            d = np.hypot(((self.nux2 - self.nux1 + 0.5) % 1.0) - 0.5, ((self.nuy2 - self.nuy1 + 0.5) % 1.0) - 0.5)
            self.D = np.where(d > 0, np.log10(d), -16.0)
        self.D[lost > 0] = np.nan
        self.D[~np.isfinite(self.nux1) | ~np.isfinite(self.nux2)] = np.nan
        # v2.3.0 S2: 周回数を半分にした解析 (最初の nturn/2 周) の D_half と、周回数への依存の傾き slope = (D − D_half)/log10(2)。
        # 規則的な点 (NAFF の誤差だけ) は周回数 2 倍で D が約 1.1 下がる (slope ≈ −3.7)、カオス的な点は ≈ 0 (D が周回数にほとんど依らない)
        self.D_half = np.full(len(self.x), np.nan)
        if nu.shape[1] > 8:
            dh = nu[:, 8]
            with np.errstate(all='ignore'):
                self.D_half = np.where(dh > 0, np.log10(dh), np.where(dh == 0, -16.0, np.nan))
            self.D_half[~np.isfinite(self.D)] = np.nan
        self.warnings = []
        self.grid = grid
        self.cond = dict(cond)
        self.setup = setup
        self.loc = loc

    @classmethod
    def load_table(cls, path):
        """table() / save() が書いた表から作る (条件の行 `# key = value` は cond に文字列で入る。setup は None)"""
        cond = {}
        rows = []
        for l in open(path, encoding='utf-8'):
            if l.startswith('# ') and ' = ' in l:
                k, v = l[2:].strip().split(' = ', 1)
                cond[k] = v
            elif not l.startswith('#') and l.strip():
                rows.append([float(t) for t in l.split()])
        a = np.array(rows)
        lost = a[:, 8].astype(np.int64)
        nu = np.column_stack([a[:, 3], a[:, 4], a[:, 5], a[:, 6], np.full((len(a), 4), np.nan)])
        g = cond.get('grid_nx_ny_ndp', '(1, 1, 1)')
        grid = tuple(int(t) for t in g.strip('()').split(','))
        res = cls(a[:, 0], a[:, 1], a[:, 2], nu, lost, grid, cond, setup=None, loc=cond.get('loc'))
        if a.shape[1] > 9:
            res.D_half = a[:, 9].copy()
        res.D = a[:, 7].copy()              # 表の D (チューンの桁数が 10 桁なので、D は表の値を使う。失われた点は nan)
        return res

    @property
    def slope(self):
        """D の周回数への依存の傾き (D − D_half)/log10(2) (周回数 1 桁あたりの D の変化)。D_half がない点は nan"""
        with np.errstate(all='ignore'):
            return (self.D - self.D_half) / math.log10(2.0)

    @property
    def nturn(self):
        try:
            return int(self.cond.get('nturn', 1024))
        except Exception:
            return 1024

    @property
    def D_norm(self):
        """1024 周に換算した D: D + κ log2(nturn/1024) (κ = fmapcheck.SCALING_KAPPA = 1.0。規則的な点の D は周回数が 2 倍になると約 1 下がる。S2 の測定 (512〜4096 周))。
        周回数の違う結果を比べるとき、D > −5 などの閾値は、周回数に依らず D_norm に対して使う (D > −5 − κ log2(nturn/1024) と同じ)"""
        from .fmapcheck import SCALING_KAPPA
        return self.D + SCALING_KAPPA * math.log2(self.nturn / 1024.0)

    def chaotic(self, thr=-5.0):
        """周回数に依らないカオス的な点の判定: D_norm > thr (1024 周に換算した D が thr より大きい)。点ごとの傾き slope は単独では判定に使えない
        (点ごとの推定の揺れが大きく、周回数を変えたときの判定の一致率が低い。D の閾値を周回数で換算する方が一致する、S2)"""
        with np.errstate(all='ignore'):
            return np.isfinite(self.D) & (self.D_norm > thr)

    @property
    def survived(self):
        return self.lost_turn == 0

    def __len__(self):
        return len(self.x)

    def table(self):
        lines = ['# frequency map: D = log10 sqrt(dnux^2 + dnuy^2), first half vs second half of the turns (NAFF); lost: lost_turn > 0 (D = nan)']
        for k, v in self.cond.items():
            lines.append('# %s = %s' % (k, v))
        for w in self.warnings:
            lines.append('# warning = %s' % w[1].replace(chr(10), ' '))
        lines.append('# x[m]  y[m]  dp  nux1  nuy1  nux2  nuy2  D  lost_turn  D_half')
        for i in range(len(self)):
            lines.append('%13.6e %13.6e %13.6e %.10f %.10f %.10f %.10f %9.5f %6d %9.5f' % (
                self.x[i], self.y[i], self.dp[i], self.nux1[i], self.nuy1[i], self.nux2[i], self.nuy2[i], self.D[i], self.lost_turn[i], self.D_half[i]))
        return '\n'.join(lines) + '\n'

    def save(self, path):
        with open(path, 'w', encoding='utf-8') as f:
            f.write(self.table())

    def as_columns(self):
        return dict(x=self.x, y=self.y, dp=self.dp, nux1=self.nux1, nuy1=self.nuy1, nux2=self.nux2, nuy2=self.nuy2, D=self.D,
                    lost_turn=self.lost_turn)

    def plot(self, path, vmin=-8.0, vmax=-2.0, order=4, lines=True, title=None, axes=None, tune_window=None, scale=1e3):
        """2 枚の図 (1 つのファイル): 左 = チューン平面 (νx1, νy1) を D で色づけ (共鳴線 order 次以下を重ねる)、右 = 初期条件の平面を D で
        色づけ (失われた点は黒い x)。axes: 初期条件平面の 2 軸 ('x','y')・('x','dp')・('y','dp') (None なら配列の長さが 1 でない 2 つ)"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from .lma import _jfont
        _jfont(plt)
        ax_names = axes or _plane_names(self.grid)
        cols = {'x': (self.x * scale, 'x [mm]'), 'y': (self.y * scale, 'y [mm]'), 'dp': (self.dp * 100, 'dp [%]')}
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
        ok = np.isfinite(self.D)
        sc = a1.scatter(self.nux1[ok], self.nuy1[ok], c=self.D[ok], s=7, cmap='turbo', vmin=vmin, vmax=vmax)
        if tune_window is None:
            if ok.any():
                # 外れた点 (2%〜98% の外) は窓に入れない (図の共鳴線が見えなくなるので)
                px, py = self.nux1[ok], self.nuy1[ok]
                x0, x1 = np.percentile(px, [2, 98])
                y0, y1 = np.percentile(py, [2, 98])
                mx, my = 0.08 * max(x1 - x0, 1e-3), 0.08 * max(y1 - y0, 1e-3)
                tune_window = ((x0 - mx, x1 + mx), (y0 - my, y1 + my))
            else:
                tune_window = ((0, 1), (0, 1))
        if lines:
            for (m, n, p) in resonance_lines(order, tune_window):
                col = {1: 'k', 2: 'k', 3: 'tab:purple', 4: 'tab:green'}.get(abs(m) + abs(n), 'gray')
                lw = {1: 1.2, 2: 1.0, 3: 0.8, 4: 0.6}.get(abs(m) + abs(n), 0.4)
                _draw_line(a1, m, n, p, tune_window, col, lw)
        a1.set_xlim(*tune_window[0])
        a1.set_ylim(*tune_window[1])
        a1.set_xlabel('nu_x (fractional, first half)')
        a1.set_ylabel('nu_y (fractional, first half)')
        a1.set_title('tune plane (resonance lines up to order %d)' % order)
        fig.colorbar(sc, ax=a1, label='D = log10 |dnu|')
        u, ul = cols[ax_names[0]]
        v, vl = cols[ax_names[1]]
        a2.scatter(u[ok], v[ok], c=self.D[ok], s=14, marker='s', cmap='turbo', vmin=vmin, vmax=vmax)
        lost = self.lost_turn > 0
        if lost.any():
            a2.scatter(u[lost], v[lost], c='k', s=14, marker='x', linewidths=0.8, label='lost (%d)' % lost.sum())
            a2.legend(fontsize=8, loc='upper right')
        a2.set_xlabel(ul)
        a2.set_ylabel(vl)
        a2.set_title(title or 'initial conditions at loc %s' % self.loc)
        fig.colorbar(sc, ax=a2, label='D = log10 |dnu|')
        fig.tight_layout()
        fig.savefig(path, dpi=130)
        plt.close(fig)


def _plane_names(grid):
    names = [n for n, g in zip(('x', 'y', 'dp'), grid) if g > 1]
    if len(names) == 2:
        return tuple(names)
    if len(names) == 0:
        return ('x', 'y')
    if len(names) == 1:
        return (names[0], 'y' if names[0] != 'y' else 'x')
    return ('x', 'y')


def _arr(v):
    a = np.atleast_1d(np.asarray(v, dtype=float)).ravel()
    return a


def frequency_map(setup, loc, x=0.0, y=0.0, dp=0.0, nturn=1024, window=1, workers=1, chunk=400, progress=0, orbit='auto', norm='data', scaling=True, warn=True):
    """周波数マップ。setup: trksetup.TrackSetup (通常 rfsw=False・rad=False の複製: setup.copy(rad=False, rfsw=False))。loc: 位置 (1 始まり)。
    x, y, dp: 数または配列 (閉軌道からの初期のずれ [m, m, δ])。nturn: 周回数の合計 (前半・後半に nturn//2 周ずつ。偶数にする)。
    window: NAFF の窓の次数 (1 = Hann)。workers: 子プロセスの数 (1 なら親プロセスで)。chunk: 一度に追跡する粒子の数。
    norm: 'data' (既定。(x, px) の楕円の規格化を標本から求める) または 'twiss' (位置 loc の Twiss の δ = 0 の値)。
    orbit: 'auto' (RFSW OFF は δ の閉軌道を Newton 法で、ON は線形の分散)、'linear' (常に線形の分散)、'none' (閉軌道に δ を足すだけ)。
    scaling: True (v2.3.0 S2 の既定) なら最初の nturn/2 周の D_half も求め、結果の slope・chaotic() (周回数に依らないカオスの判定) に使う (NAFF は追跡の 1% の時間)。
    warn: True (既定) なら、信頼できない・信頼度の低い条件 (周回数が少ない、窓なし、RFSW・RAD・FLUC ON、作動点が整数・半整数に近い、誤差つき・結合が強い、
    失われた点が多い、NAFF の失敗など。fmapcheck.check_setup・check_result) を FMapWarning で知らせ、結果の warnings (リスト [(コード, 文)]) にも入れる"""
    xs, ys, dps = _arr(x), _arr(y), _arr(dp)
    nturn = int(nturn)
    if nturn < 64:
        raise ValueError('nturn が少ない (64 以上)')
    nturn -= nturn % 2
    wlist = []
    if warn:
        from . import fmapcheck as _fc
        wlist = _fc.check_setup(setup, nturn, window, loc)
    # 並列でない (親プロセスで走る) ときも、設定の K1 の定義・端面の扱い (compat) に合わせる (並列の子プロセスは _init_worker・TrackSetup.__setstate__ が合わせる)
    _saved = _compat_apply(setup) if int(workers) <= 1 else None
    try:
        return _frequency_map_body(setup, loc, xs, ys, dps, nturn, window, workers, chunk, progress, orbit, norm, scaling, warn, wlist)
    finally:
        _compat_restore(_saved)


def _compat_apply(setup):
    from . import compat
    saved = {}
    for nm, val in (('bend_k1_edge', getattr(setup, 'bend_k1_edge', None)), ('bend_k1_def', getattr(setup, 'bend_k1_def', None))):
        get, st = getattr(compat, nm, None), getattr(compat, 'set_' + nm, None)
        if get is None or st is None or val is None:
            continue
        try:
            saved[nm] = get()
            if saved[nm] != val and compat.mode() == 'v2':
                st(val)
        except Exception:
            saved.pop(nm, None)
    return saved


def _compat_restore(saved):
    from . import compat
    for nm, val in (saved or {}).items():
        try:
            getattr(compat, 'set_' + nm)(val)
        except Exception:
            pass


def _frequency_map_body(setup, loc, xs, ys, dps, nturn, window, workers, chunk, progress, orbit, norm, scaling, warn, wlist):
    X, Y, DP = np.meshgrid(xs, ys, dps, indexing='ij')
    X, Y, DP = X.ravel(), Y.ravel(), DP.ravel()
    n = len(X)
    base = setup.closed_orbit([loc])[0] if setup.rfsw else np.zeros(6)
    tw = setup.twiss_at(loc)
    orb = {}
    coords = np.zeros((n, 6))
    for i in range(n):
        d = float(DP[i])
        if d not in orb:
            if setup.rfsw or orbit == 'linear':
                c = base.copy()
                c[5] += d
                c[0] += tw['etx'] * d
                c[1] += tw['etpx'] * d
            elif orbit == 'none':
                c = base.copy()
                c[5] += d
            else:
                try:
                    c = setup.closed_orbit_dp(loc, d)
                except Exception as e:
                    raise RuntimeError('δ = %g の閉軌道が求まらない (線形の運動が不安定、誤差・δ が大きすぎる可能性): %r' % (d, e))
            orb[d] = c
        coords[i] = orb[d]
        coords[i, 0] += X[i]
        coords[i, 2] += Y[i]
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))
    jobs = [(loc, coords[a:a + chunk], nturn, window, twp, norm, bool(scaling)) for a in range(0, n, chunk)]
    t0 = time.time()
    outs = []
    workers = max(1, min(int(workers), max(1, (os.cpu_count() or 2) - 1), len(jobs)))
    if workers == 1:
        for k, j in enumerate(jobs):
            outs.append(_fmap_chunk((setup,) + j))
            if progress:
                print('    [FMAP] %d/%d (%.0f s)' % (k + 1, len(jobs), time.time() - t0), flush=True)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=(setup,)) as ex:
            for k, o in enumerate(ex.map(_chunk_job, jobs, chunksize=1)):
                outs.append(o)
                if progress:
                    print('    [FMAP] %d/%d (%.0f s)' % (k + 1, len(jobs), time.time() - t0), flush=True)
    lost = np.concatenate([o[0] for o in outs])
    nu = np.concatenate([o[1] for o in outs])
    cond = setup.conditions()
    cond.update(loc=loc, s_m=setup.s_of(loc), name=setup.names[loc - 1], nturn=nturn, window=window, grid_nx_ny_ndp=(len(xs), len(ys), len(dps)),
                x_range=(float(xs.min()), float(xs.max())), y_range=(float(ys.min()), float(ys.max())), dp_range=(float(dps.min()), float(dps.max())),
                npoints=n, nlost=int((lost > 0).sum()), workers=workers, seconds=round(time.time() - t0, 1),
                betx=tw['betx'], bety=tw['bety'], orbit=orbit, norm=norm, tunes_full=tuple(round(v, 6) for v in setup.tunes_full()))
    try:
        cond['periodicity_cells'] = _ts.find_periodicity(setup.els)['nper']       # セルの数 (1 = 周期的でない。誤差つきは 1)
    except Exception:
        pass
    res = FMapResult(X, Y, DP, nu, lost, (len(xs), len(ys), len(dps)), cond, setup=setup, loc=loc)
    if warn:
        from . import fmapcheck as _fc
        wlist = list(wlist) + _fc.check_result(res, nturn)
        res.warnings = wlist
        res.cond['warnings'] = ','.join(w[0] for w in wlist) or 'none'
        for (_code, msg) in wlist:
            warnings.warn('FrequencyMap: ' + msg, FMapWarning, stacklevel=4)
    return res


# ---------------------------------------------------------------------------------------------------------------------
# 共鳴線
# ---------------------------------------------------------------------------------------------------------------------

def _mn_list(order):
    """order 次以下 (|m|+|n| ≤ order) の (m, n) (m ≥ 0、m = 0 なら n > 0。(m, n) と (−m, −n) は同じ線)"""
    out = []
    for o in range(1, order + 1):
        for m in range(0, o + 1):
            for n in {o - m, -(o - m)}:
                if m == 0 and n <= 0:
                    continue
                if (m, n) not in out:
                    out.append((m, n))
    return out


def resonance_lines(order, window=((0.0, 1.0), (0.0, 1.0))):
    """m νx + n νy = p の線で、(νx, νy) の窓 window = ((x0, x1), (y0, y1)) と交わるもの。戻り [(m, n, p)]。order 次以下 (|m|+|n| ≤ order)"""
    (x0, x1), (y0, y1) = window
    out = []
    for (m, n) in _mn_list(order):
        corners = [m * x + n * y for x in (x0, x1) for y in (y0, y1)]
        for p in range(int(math.ceil(min(corners) - 1e-12)), int(math.floor(max(corners) + 1e-12)) + 1):
            out.append((m, n, p))
    return out


def _draw_line(ax, m, n, p, window, col, lw):
    (x0, x1), (y0, y1) = window
    if n == 0:
        xv = p / m
        if x0 <= xv <= x1:
            ax.plot([xv, xv], [y0, y1], color=col, lw=lw, alpha=0.8)
    elif m == 0:
        yv = p / n
        if y0 <= yv <= y1:
            ax.plot([x0, x1], [yv, yv], color=col, lw=lw, alpha=0.8)
    else:
        xs = np.array([x0, x1])
        ys = (p - m * xs) / n
        ax.plot(xs, ys, color=col, lw=lw, alpha=0.8)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)


def resonance_check(nux, nuy, order=4, maxn=12):
    """作動点 (νx, νy) (整数部を含む) に近い共鳴線 m νx + n νy = p (|m|+|n| ≤ order) の表。距離 = |m νx + n νy − p| / sqrt(m²+n²) (チューン平面の
    距離)。距離の小さい順に maxn 個 (整数 (νx = p, νy = p)・半整数 (2νx = p, 2νy = p) の 4 本は、順位に関わらず必ず入れる)。整数・半整数の名前は 'integer'・'half-integer'。
    戻り [(m, n, p, 距離, 次数, 名前)]"""
    rows = []
    for (m, n) in _mn_list(order):
        v = m * nux + n * nuy
        p = round(v)
        rows.append((m, n, int(p), abs(v - p) / math.hypot(m, n), abs(m) + abs(n)))
    rows.sort(key=lambda r: r[3])
    must = [r for r in rows if (r[0], r[1]) in ((1, 0), (0, 1), (2, 0), (0, 2))]
    sel = rows[:maxn] + [r for r in must if r not in rows[:maxn]]
    out = []
    for (m, n, p, d, o) in sel:
        nm = ('integer' if o == 1 else 'half-integer' if (o == 2 and (m == 0 or n == 0)) else 'sum' if (m * n > 0) else 'difference'
              if (m and n) else 'normal') + ' order %d' % o
        out.append((m, n, p, d, o, nm))
    return out


def line_distance(nux, nuy, order=4):
    """各点 (νx, νy) (配列、端数) から、order 次以下の共鳴線 (どの p でも) までの最短距離 (チューン平面の距離)"""
    nux = np.asarray(nux, dtype=float)
    nuy = np.asarray(nuy, dtype=float)
    best = np.full(nux.shape, np.inf)
    for (m, n) in _mn_list(order):
        v = m * nux + n * nuy
        d = np.abs(v - np.round(v)) / math.hypot(m, n)
        best = np.minimum(best, d)
    return best
