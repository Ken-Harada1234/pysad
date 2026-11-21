# -*- coding: utf-8 -*-
"""v2.3.0 S1: 周波数マップの追加の検証 (b)〜(f)。モード:
  k64 <ring>            (c) 互換モード k64 (PYSAD_COMPAT=k64 と同じ) で同じ周波数マップを走らせ、v2 (Q1 の表) と点ごとに比べる
  errors <ring> <level> (d) 四極・6 極・偏向磁石の位置のずれ・ロールの乱数の誤差 (要素ごとに別々、周期性を使わない全周の追跡) つきの周波数マップ
  mpw                   (e) 挿入光源 (多極ウィグラー MPW) つきの ASTRID2 (ASTRID2_v01_02_id.sad の RINGMPW)
  delta <ring> <plane>  (b) x–δ・y–δ の周波数マップ (RFSW・RAD OFF、δ の閉軌道の上)。δ ごとの線形のチューン・RDT・ADTS から共鳴の予想と比べる
  naff <ring>           (f) NAFF の窓 (Hann の次数)・周回数への依存 (追跡を 1 回して、周回数・窓を変えて NAFF をかけ直す)
  grid <ring> <nx> <ny> (f) 格子の細かさへの依存 (細かい格子の周波数マップ)
  lambda <ring>         (a の原因の追及) 6 極の強さを λ 倍した周波数マップ (D > -5 の点の数と、それが載る共鳴線)
使い方 (リポジトリ直下から): py -3.12 example/v2_3/S1_fmap_variants.py <モード> <引数> [--workers N]    出力: 標準出力、S1_fmap_<モード>_*.table.txt・.png
共通の条件: v2 モード (k64 の節を除く)、RFSW OFF・RAD OFF、LOSSAMPL = 0.1 m、px = py = z = 0 から出発、NAFF の窓 Hann (次数 1)、規格化 = 標本の楕円。
  ASTRID2 (0.58 GeV、周長 45.7 m): 位置 2、x = 0.5〜18 mm、y = 0.5〜14 mm。kek-ls_v03_12.lat (3 GeV、603 m): 位置 3、x = 0.2〜6 mm、y = 0.2〜8 mm。
  格子は full = 41 × 31 (Q1 と同じ)、coarse = 21 × 16 (full の 1 つおき)。周回数 1024 (前半・後半 512 周)。
"""
import copy
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

RINGS = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), line='RING', loc=2, x=(0.0005, 0.018), y=(0.0005, 0.014),
                    q1=os.path.join(HERE, 'Q1_fmap_astrid2.table.txt'), label='ASTRID2 (0.58 GeV、45.7 m)'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), line='RING', loc=3, x=(0.0002, 0.006), y=(0.0002, 0.008),
                  q1=os.path.join(HERE, 'Q1_fmap_kekls.table.txt'), label='kek-ls_v03_12.lat (3 GeV、603 m)'),
    'mpw': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), line='RINGMPW', loc=2, x=(0.0005, 0.018), y=(0.0005, 0.014),
                q1=None, label='ASTRID2 + MPW (ASTRID2_v01_02_id.sad の RINGMPW、0.58 GeV)'),
}
NTURN = 1024


def grid_of(ring, kind):
    c = RINGS[ring]
    nx, ny = (41, 31) if kind == 'full' else (21, 16)
    return np.linspace(*c['x'], nx), np.linspace(*c['y'], ny)


def make_setup(ring):
    from pysad.trksetup import TrackSetup
    c = RINGS[ring]
    src = open(c['path'], encoding='utf-8').read()
    st, _ffs = TrackSetup.from_source(src, c['line'], rad=False, rfsw=False, alost=0.1)
    return st, src


def describe(r, tag=''):
    from pysad.freqmap import line_distance
    ok = np.isfinite(r.D)
    D = r.D[ok]
    q = np.percentile(D, [5, 25, 50, 75, 95]) if ok.any() else [np.nan] * 5
    print('%s点 %d、失われた点 %d、D の範囲 %.2f 〜 %.2f、D の百分位 (5, 25, 50, 75, 95%%) = %s、D > -5: %d 点 (%.3f)、D > -4: %d 点 (%.3f)' % (
        tag, len(r), (r.lost_turn > 0).sum(), D.min(), D.max(), ', '.join('%.2f' % v for v in q), (D > -5).sum(), (D > -5).mean(), (D > -4).sum(), (D > -4).mean()))
    d = line_distance(r.nux1[ok], r.nuy1[ok], 4)
    near = d < 0.003
    print('%s  次数 4 以下の共鳴線から 0.003 以内の点 %d、それ以外 %d。D > -5 の割合 近い / 遠い = %.3f / %.3f' % (
        tag, near.sum(), (~near).sum(), np.mean(D[near] > -5) if near.any() else np.nan, np.mean(D[~near] > -5) if (~near).any() else np.nan))


def top_lines(r, thr=-4.0, n=5, order=5):
    from collections import Counter
    from pysad.freqmap import resonance_lines
    cnt = Counter()
    ok = np.isfinite(r.D)
    for i in np.nonzero(ok & (r.D > thr))[0]:
        best = None
        for (m, k, p) in resonance_lines(order, ((r.nux1[i] - 0.05, r.nux1[i] + 0.05), (r.nuy1[i] - 0.05, r.nuy1[i] + 0.05))):
            dd = abs(m * r.nux1[i] + k * r.nuy1[i] - p) / np.hypot(m, k)
            if best is None or dd < best[0]:
                best = (dd, m, k, p)
        if best:
            cnt[(best[1], best[2])] += 1
    print('  D > %g の点 (%d 点) の最も近い共鳴線 (次数 %d 以下、端数の平面) の上位: %s' % (
        thr, int((ok & (r.D > thr)).sum()), order, ', '.join('%d νx %+d νy: %d 点' % (k[0], k[1], v) for k, v in cnt.most_common(n))))


def load_ref(ring):
    from pysad.freqmap import FMapResult
    return FMapResult.load_table(RINGS[ring]['q1'])


def sub(r, step=2):
    """full の格子の結果から 1 つおきの点 (coarse の格子) を取り出す"""
    from pysad.freqmap import FMapResult
    nx, ny = r.grid[0], r.grid[1]
    idx = np.arange(nx * ny).reshape(nx, ny)[::step, ::step].ravel()
    nu = np.column_stack([r.nux1[idx], r.nuy1[idx], r.nux2[idx], r.nuy2[idx], np.full((len(idx), 4), np.nan)])
    o = FMapResult(r.x[idx], r.y[idx], r.dp[idx], nu, r.lost_turn[idx], ((nx + 1) // step, (ny + 1) // step, 1), r.cond, loc=r.loc)
    o.D = r.D[idx].copy()
    return o


def compare(a, b, la, lb):
    both = np.isfinite(a.D) & np.isfinite(b.D)
    dD = np.abs(a.D[both] - b.D[both])
    dnux = np.abs(((a.nux1[both] - b.nux1[both] + 0.5) % 1) - 0.5)
    dnuy = np.abs(((a.nuy1[both] - b.nuy1[both] + 0.5) % 1) - 0.5)
    lost_diff = int(((a.lost_turn > 0) != (b.lost_turn > 0)).sum())
    print('%s と %s の比較: 両方で残った点 %d、失われた点の食い違い %d (%s: %d、%s: %d)、|ΔD| の最大 %.3g・中央値 %.3g・|ΔD| > 0.1 の点 %d、|Δνx1| の最大 %.3g、|Δνy1| の最大 %.3g' % (
        la, lb, both.sum(), lost_diff, la, (a.lost_turn > 0).sum(), lb, (b.lost_turn > 0).sum(), dD.max() if len(dD) else np.nan, np.median(dD) if len(dD) else np.nan,
        (dD > 0.1).sum(), dnux.max() if len(dnux) else np.nan, dnuy.max() if len(dnuy) else np.nan))


def run_map(st, loc, xs, ys, workers, dp=0.0, nturn=NTURN):
    from pysad.freqmap import frequency_map
    t0 = time.time()
    r = frequency_map(st, loc, x=xs, y=ys, dp=dp, nturn=nturn, workers=workers)
    print('追跡と NAFF: %.0f s (%d 点、%d 並列)' % (time.time() - t0, len(r), r.cond['workers']))
    return r


# ---------------------------------------------------------------------------------------------------------------------
def mode_k64(ring, workers):
    from pysad import compat
    compat.set_mode('k64')
    c = RINGS[ring]
    print('(c) 互換モード k64 の周波数マップ: %s、位置 %d' % (c['label'], c['loc']))
    kind = 'full' if ring == 'astrid2' else 'coarse'
    st, _ = make_setup(ring)
    print('compat mode =', compat.mode(), '、作動点 νx, νy = %.6f, %.6f' % st.tunes_full())
    xs, ys = grid_of(ring, kind)
    r = run_map(st, c['loc'], xs, ys, workers)
    r.save(os.path.join(HERE, 'S1_fmap_k64_%s.table.txt' % ring))
    r.plot(os.path.join(HERE, 'S1_fmap_k64_%s.png' % ring), order=4)
    describe(r, 'k64: ')
    ref = load_ref(ring)
    if kind == 'coarse':
        ref = sub(ref)
    describe(ref, 'v2 (Q1 の表): ')
    compare(r, ref, 'k64', 'v2')
    top_lines(r)
    top_lines(ref)


def mode_errors(ring, level, workers):
    from pysad import compat, fmapres as FR
    compat.set_mode('v2')
    c = RINGS[ring]
    lev = float(level)
    sig = dict(quad_dxy=50e-6 * lev, quad_rot=0.2e-3 * lev, sext_dxy=50e-6 * lev, bend_dxy=50e-6 * lev, bend_rot=0.1e-3 * lev)
    print('(d) 誤差つきの周波数マップ: %s、位置 %d、誤差のレベル ×%g (標準偏差 [m, rad]: %s)、乱数の種 1、要素ごとに別々の誤差 (周期性なし)' % (
        c['label'], c['loc'], lev, sig))
    st0, _ = make_setup(ring)
    st, cnt = FR.apply_alignment_errors(st0, seed=1, **sig)
    print('誤差を与えた要素の数:', cnt)
    c0 = st0.closed_orbit_dp(c['loc'], 0.0)
    c1 = st.closed_orbit_dp(c['loc'], 0.0)
    print('位置 %d の閉軌道 (x, px, y, py) [mm, mrad]: 誤差なし %s、誤差つき %s、Newton 法の残差 %.1e' % (
        c['loc'], np.round(c0[:4] * 1e3, 5), np.round(c1[:4] * 1e3, 4), st.info.get('orbit_resid', float('nan'))))
    kind = 'full' if (ring == 'astrid2' and lev <= 1.0) else 'coarse'
    xs, ys = grid_of(ring, kind)
    r = run_map(st, c['loc'], xs, ys, workers)
    tag = 'L%g' % lev
    r.save(os.path.join(HERE, 'S1_fmap_errors_%s_%s.table.txt' % (ring, tag)))
    r.plot(os.path.join(HERE, 'S1_fmap_errors_%s_%s.png' % (ring, tag)), order=4)
    describe(r, '誤差つき: ')
    ref = load_ref(ring)
    if kind == 'coarse':
        ref = sub(ref)
    describe(ref, '誤差なし (Q1 の表): ')
    compare(r, ref, '誤差つき', '誤差なし')
    i0 = 0
    print('最小振幅の点 (x = %.2f mm、y = %.2f mm) のチューン: 誤差なし (%.6f, %.6f)、誤差つき (%.6f, %.6f)、D: %.2f → %.2f' % (
        r.x[i0] * 1e3, r.y[i0] * 1e3, ref.nux1[i0], ref.nuy1[i0], r.nux1[i0], r.nuy1[i0], ref.D[i0], r.D[i0]))
    top_lines(r)


def mode_mpw(workers):
    from pysad import compat, fmapres as FR
    from pysad.freqmap import resonance_check
    compat.set_mode('v2')
    c = RINGS['mpw']
    st, src = make_setup('mpw')
    print('(e) 挿入光源 (MPW) つきの周波数マップ: %s、位置 %d、nlat = %d' % (c['label'], c['loc'], st.nlat))
    nux, nuy = st.tunes_full()
    print('作動点 (CALC): νx = %.6f、νy = %.6f (MPW なし 5.179990・2.140178)、周長 %.4f m' % (nux, nuy, st.circ))
    print('作動点に近い共鳴線 (次数 4 以下):')
    for (m, n, p, d, o, nm) in resonance_check(nux, nuy, 4, maxn=8):
        print('    %2d νx %+2d νy = %4d  距離 %.4f  %s' % (m, n, p, d, nm))
    xs, ys = grid_of('mpw', 'full')
    r = run_map(st, c['loc'], xs, ys, workers)
    r.save(os.path.join(HERE, 'S1_fmap_mpw.table.txt'))
    r.plot(os.path.join(HERE, 'S1_fmap_mpw.png'), order=4)
    describe(r, 'MPW つき: ')
    ref = load_ref('astrid2')
    describe(ref, 'MPW なし (Q1 の表): ')
    compare(r, ref, 'MPW つき', 'MPW なし')
    top_lines(r)
    # 失われた点の分布
    lost = r.lost_turn > 0
    if lost.any():
        print('失われた点 %d: x の範囲 %.1f〜%.1f mm、y の範囲 %.1f〜%.1f mm、落ちた周の最小・中央値 %d・%d' % (
            lost.sum(), r.x[lost].min() * 1e3, r.x[lost].max() * 1e3, r.y[lost].min() * 1e3, r.y[lost].max() * 1e3, r.lost_turn[lost].min(), np.median(r.lost_turn[lost])))
    # RDT・ADTS の予想 (偏向磁石 K0 ≠ 0 の MPW の内部は resdrive の対象外)
    try:
        from pysad.parser import load_sad_string
        lat = load_sad_string(src).build_lattice(c['line'])
        rows = lat.twiss_table()
        pr = FR.ResonancePredictor.from_lattice(lat, c['loc'], rows)
        ip = (np.floor(pr.nu0[0]), np.floor(pr.nu0[1]))
        nu_meas = (r.nux1 + ip[0], r.nuy1 + ip[1])
        print('MPW つきの ADTS (全部入り) xx, xy, yy = %.1f, %.1f, %.1f、1 次の RDT の調和: %s' % (
            pr.A[0, 0], pr.A[0, 1], pr.A[1, 1], {a: round(max(abs(q[4]) for q in t), 3) for a, t in pr.harms.items()}))
        P = pr.predict(r.x, r.y, nu=nu_meas)
        for thr in (-5.0, -4.0):
            s = FR.coverage_stats(r.D, P['any'], r.grid, thr=thr)
            print('  1 次の RDT の共鳴の予想 (位相つき、追跡のチューン) / D > %g: area %.3f、cov %.3f、対照 %.3f ± %.3f、p = %.3f' % (
                thr, s['area'], s['cov'], s['shift_mean'], s['shift_std'], s['pval']))
    except Exception as e:      # noqa
        print('MPW つきの RDT・ADTS の予想は計算できなかった:', repr(e))


def mode_delta(ring, plane, workers):
    from pysad import compat, fmapres as FR
    from pysad.freqmap import line_distance, resonance_check
    from pysad.parser import load_sad_string
    compat.set_mode('v2')
    c = RINGS[ring]
    st, src = make_setup(ring)
    nx = 41
    ndp = 31
    dpmax = 0.02
    dps = np.linspace(-dpmax, dpmax, ndp)
    if plane == 'xdp':
        xs, ys = np.linspace(*c['x'], nx), [2e-3 if ring == 'astrid2' else 1e-3]
    else:
        xs, ys = [2e-3 if ring == 'astrid2' else 1e-3], np.linspace(*c['y'], nx)
    print('(b) δ を含む周波数マップ (%s): %s、位置 %d、RFSW・RAD OFF、δ = ±%g%% (%d 点、δ の閉軌道 (4 次元の Newton 法) の上)、固定した %s = %.1f mm、%d 周' % (
        plane, c['label'], c['loc'], dpmax * 100, ndp, 'y' if plane == 'xdp' else 'x', (ys if plane == 'xdp' else xs)[0] * 1e3, NTURN))
    tabp = os.path.join(HERE, 'S1_fmap_delta_%s_%s.table.txt' % (ring, plane))
    if '--reuse' in sys.argv and os.path.exists(tabp):
        from pysad.freqmap import FMapResult
        r = FMapResult.load_table(tabp)
        print('保存済みの表 %s を読む' % os.path.basename(tabp))
    else:
        r = run_map(st, c['loc'], xs, ys, workers, dp=dps)
        r.save(tabp)
    describe(r, 'δ マップ: ')
    top_lines(r)
    # δ ごとの線形のチューン・ADTS・RDT (lat.twiss_table(delta=δ))
    lat = load_sad_string(src).build_lattice(c['line'])
    ip = None
    shape = (nx, ndp)
    X = r.x.reshape(shape)
    Y = r.y.reshape(shape)
    DP = r.dp.reshape(shape)
    nuxm = r.nux1.reshape(shape)
    nuym = r.nuy1.reshape(shape)
    Dm = r.D.reshape(shape)
    pred_x = np.full(shape, np.nan)
    pred_y = np.full(shape, np.nan)
    lin_x = np.full(ndp, np.nan)
    lin_y = np.full(ndp, np.nan)
    masks = {k: np.zeros(shape, dtype=bool) for k in ('trapped', 'envelope')}
    t0 = time.time()
    for j, d in enumerate(dps):
        rows = lat.twiss_table(delta=float(d))
        try:
            pr = FR.ResonancePredictor.from_lattice(lat, c['loc'], rows, delta=float(d))
        except Exception as e:      # noqa
            print('  δ = %+.4f の予想が計算できなかった: %r' % (d, e))
            continue
        if ip is None:
            ip = (math.floor(pr.nu0[0]), math.floor(pr.nu0[1]))
        lin_x[j], lin_y[j] = pr.nu0
        xcol, ycol = X[:, j], Y[:, j]
        # x, y は δ の閉軌道からのずれ。位置 loc の Twiss (δ) の光学で作用を求める
        nu_pred = pr.tunes_lin(*pr.actions(xcol, ycol)[:2])
        pred_x[:, j], pred_y[:, j] = nu_pred[0] - ip[0], nu_pred[1] - ip[1]
        P = pr.predict(xcol, ycol, nu=(nuxm[:, j] + ip[0], nuym[:, j] + ip[1]))
        masks['trapped'][:, j] = P['any']
        masks['envelope'][:, j] = P['envelope']
    print('δ ごとの RDT・ADTS・線形のチューン: %.0f s' % (time.time() - t0))
    # (1) 線形のクロマティックなチューン ν0(δ) と、追跡で測ったチューン (最小振幅の列) の比較
    imin = 0 if plane == 'xdp' else 0
    print('線形のチューン (lat.twiss_table(delta)) と、追跡 (NAFF) の最小振幅の点のチューンの比較 (整数部を除く) [δ, ν0x(δ), νx(追跡), ν0y(δ), νy(追跡)]:')
    errs = []
    for j in range(0, ndp, 5):
        print('   δ = %+6.3f%%: νx %.5f / %.5f、νy %.5f / %.5f' % (dps[j] * 100, lin_x[j] - ip[0], nuxm[imin, j], lin_y[j] - ip[1], nuym[imin, j]))
    ex = (lin_x - ip[0]) - nuxm[imin, :]
    ey = (lin_y - ip[1]) - nuym[imin, :]
    print('   差 (線形 − 追跡) の最大 |Δνx| = %.4f、|Δνy| = %.4f (最小振幅の点も固定した %s の振幅を持つので、その ADTS 分を含む)' % (
        np.nanmax(np.abs(ex)), np.nanmax(np.abs(ey)), 'y' if plane == 'xdp' else 'x'))
    ex2 = pred_x[imin, :] - nuxm[imin, :]
    ey2 = pred_y[imin, :] - nuym[imin, :]
    print('   線形のチューン + ADTS (その点の J) との差の最大 |Δνx| = %.4f、|Δνy| = %.4f (δ ごとの ADTS、6 極 2 次 + 8 極 + 運動学・フリンジ。偏向磁石の 3 次は δ ≠ 0 では入らない)' % (np.nanmax(np.abs(ex2)), np.nanmax(np.abs(ey2))))
    # (2) 予想したチューンと追跡のチューンの差、共鳴線の近く
    ok = np.isfinite(Dm)
    errp = np.hypot(pred_x - nuxm, pred_y - nuym)
    print('線形のチューン + ADTS (δ ごと) の予想と追跡のチューンの差: 中央値 %.4f、0.003 以内の点の割合 %.2f' % (np.nanmedian(errp[ok]), np.mean(errp[ok] < 0.003)))
    dist_pred = line_distance(pred_x, pred_y, 4)
    near_pred = (dist_pred < 0.003) & ok
    near_meas = (line_distance(nuxm, nuym, 4) < 0.003) & ok
    for lab, M in (('予想した作動点 (線形 + ADTS) が次数 4 以下の共鳴線の 0.003 以内', near_pred), ('追跡のチューンが次数 4 以下の共鳴線の 0.003 以内', near_meas)):
        for thr in (-5.0, -4.0):
            s = FR.coverage_stats(r.D, M.ravel(), (nx, ndp), thr=thr)
            print('  [%s / D > %g] area %.3f、cov %.3f、cov/area %.2f、中 %.3f / 外 %.3f、対照 (ずらし) %.3f ± %.3f、p = %.3f' % (
                lab, thr, s['area'], s['cov'], s['enrich'], s['p_in'], s['p_out'], s['shift_mean'], s['shift_std'], s['pval']))
    for key in ('trapped', 'envelope'):
        for thr in (-5.0, -4.0):
            s = FR.coverage_stats(r.D, masks[key].ravel(), (nx, ndp), thr=thr)
            print('  [1 次の RDT の共鳴の予想 (δ ごとの RDT・ADTS、追跡のチューン) %s / D > %g] area %.3f、cov %.3f、対照 %.3f ± %.3f、p = %.3f' % (
                key, thr, s['area'], s['cov'], s['shift_mean'], s['shift_std'], s['pval']))
    # δ = 0 の列 (Q1 の x–y の格子の y = 固定値の行と同じはず) との整合: x–δ の δ = 0 の列
    j0 = int(np.argmin(np.abs(dps)))
    print('δ = %+.4f の列の D の中央値 %.2f、最大 %.2f; |δ| = %.1f%% の端の列の中央値 %.2f・%.2f' % (
        dps[j0], np.nanmedian(Dm[:, j0]), np.nanmax(Dm[:, j0]), dpmax * 100, np.nanmedian(Dm[:, 0]), np.nanmedian(Dm[:, -1])))
    # 図: D のマップ + 予想した共鳴線の 0.003 以内 (等高線)
    ax_names = ('x', 'dp') if plane == 'xdp' else ('y', 'dp')
    FR.plot_with_prediction(r, [('predicted tune within 0.003 of a line (order<=4)', near_pred, 'white', '-'), ('first-order RDT: trapped (measured nu)', masks['trapped'], 'magenta', '--')],
                            os.path.join(HERE, 'S1_fmap_delta_%s_%s.png' % (ring, plane)), axes=ax_names,
                            title='%s %s: D and resonance prediction (white: predicted nu(delta, J))' % (ring, plane))


def mode_naff(ring, workers):
    from pysad import compat
    from pysad.freqmap import _tunes_of
    compat.set_mode('v2')
    c = RINGS[ring]
    st, _ = make_setup(ring)
    xs, ys = grid_of(ring, 'coarse')
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    X, Y = X.ravel(), Y.ravel()
    coords = np.zeros((len(X), 6))
    coords[:, 0], coords[:, 2] = X, Y
    NMAX = 4096
    print('(f) NAFF の窓・周回数への依存: %s、位置 %d、coarse の格子 %d × %d = %d 点、追跡は %d 周を 1 回 (記録して、周回数・窓を変えて NAFF をかけ直す)' % (
        c['label'], c['loc'], len(xs), len(ys), len(X), NMAX))
    tw = st.twiss_at(c['loc'])
    twp = ((tw['alfx'], tw['betx']), (tw['alfy'], tw['bety']))
    t0 = time.time()
    res = st.track(c['loc'], coords, NMAX, record=True)
    print('追跡 %.0f s、失われた点 %d' % (time.time() - t0, (res['lost_turn'] > 0).sum()))
    rec = res['rec']
    alive = res['lost_turn'] == 0

    def Dof(N, win, norm='data'):
        D = np.full(len(X), np.nan)
        nu1 = np.full((len(X), 2), np.nan)
        for i in range(len(X)):
            if not alive[i]:
                continue
            f, _a = _tunes_of(rec[:N, i, :4], twp, win, norm)
            d = np.hypot(((f[2] - f[0] + .5) % 1) - .5, ((f[3] - f[1] + .5) % 1) - .5)
            D[i] = math.log10(d) if d > 0 else -16.0
            nu1[i] = (f[0], f[1])
        return D, nu1
    base, nu_base = Dof(1024, 1)
    np.savetxt(os.path.join(HERE, 'S1_fmap_naff_%s.base.txt' % ring), np.column_stack([X, Y, base]), header='x y D (nturn 1024, Hann p=1)')
    print('基準: 1024 周、Hann (p = 1)。D > -5: %d 点 (%.3f)、D > -4: %d 点' % ((base > -5).sum(), np.nanmean(base > -5), (base > -4).sum()))
    print('D の分布の比較 (基準との差は同じ点どうし。前半・後半の周回数 = 周回数/2):')
    print('  条件                      D の百分位 (5, 25, 50, 75, 95)          D>-5 の割合  D>-4 の割合  基準との D の差の中央値  |差|>0.5 の点  基準と D>-5 の判定が違う点   相関  順位相関  上位 10% の一致')
    rows = []
    ok = np.isfinite(base)
    for lab, N, win in (('512 周 Hann p=1', 512, 1), ('1024 周 Hann p=1 (基準)', 1024, 1), ('2048 周 Hann p=1', 2048, 1), ('4096 周 Hann p=1', 4096, 1),
                        ('1024 周 窓なし p=0', 1024, 0), ('1024 周 Hann p=2', 1024, 2), ('1024 周 Hann p=3', 1024, 3), ('2048 周 Hann p=2', 2048, 2)):
        D, _nu = Dof(N, win)
        okk = ok & np.isfinite(D)
        q = np.percentile(D[okk], [5, 25, 50, 75, 95])
        diff = D[okk] - base[okk]
        flip = int(((D[okk] > -5) != (base[okk] > -5)).sum())
        corr = float(np.corrcoef(D[okk], base[okk])[0, 1])
        from scipy.stats import spearmanr
        rho = float(spearmanr(D[okk], base[okk])[0])
        k10 = max(1, int(0.1 * okk.sum()))
        top = set(np.argsort(-np.where(okk, D, -99))[:k10])
        topb = set(np.argsort(-np.where(okk, base, -99))[:k10])
        print('  %-26s %s  %.3f  %.3f  %+.2f  %d  %d  %.3f  %.3f  %.2f' % (lab, ', '.join('%6.2f' % v for v in q), np.mean(D[okk] > -5), np.mean(D[okk] > -4), np.median(diff),
                                                           int((np.abs(diff) > 0.5).sum()), flip, corr, rho, len(top & topb) / k10))
        rows.append((lab, D))
    # D の下限 (線形の小振幅) は周回数でどう変わるか
    print('最小振幅の 5 点の D (x = %.2f mm、y = %.2f〜%.2f mm): ' % (xs[0] * 1e3, ys[0] * 1e3, ys[4] * 1e3))
    for lab, D in rows:
        print('   %-26s %s' % (lab, ', '.join('%.2f' % D[j] for j in range(5))))


def mode_grid(ring, nx, ny, workers):
    from pysad import compat
    compat.set_mode('v2')
    c = RINGS[ring]
    st, _ = make_setup(ring)
    nx, ny = int(nx), int(ny)
    print('(f) 格子の細かさ: %s、位置 %d、%d × %d = %d 点、%d 周' % (c['label'], c['loc'], nx, ny, nx * ny, NTURN))
    r = run_map(st, c['loc'], np.linspace(*c['x'], nx), np.linspace(*c['y'], ny), workers)
    r.save(os.path.join(HERE, 'S1_fmap_grid_%s_%dx%d.table.txt' % (ring, nx, ny)))
    r.plot(os.path.join(HERE, 'S1_fmap_grid_%s_%dx%d.png' % (ring, nx, ny)), order=4)
    describe(r, '%d x %d: ' % (nx, ny))


def mode_lambda(ring, workers):
    from pysad import compat
    compat.set_mode('v2')
    c = RINGS[ring]
    st0, _ = make_setup(ring)
    xs, ys = grid_of(ring, 'coarse')
    print('(a の原因の追及) 6 極の強さを λ 倍した周波数マップ: %s、位置 %d、coarse の格子 %d × %d、%d 周 (6 極の K2 を λ 倍。作動点は λ で動く)' % (c['label'], c['loc'], len(xs), len(ys), NTURN))
    for lam in (0.0, 0.5, 1.0, 1.5):
        st = st0.copy()
        st.els = list(st0.els)
        for i, e in enumerate(st.els):
            if type(e).__name__ == 'Sextupole':
                e2 = copy.copy(e)
                e2.K2 = e.K2 * lam
                st.els[i] = e2
        r = run_map(st, c['loc'], xs, ys, workers)
        r.save(os.path.join(HERE, 'S1_fmap_lambda_%s_%g.table.txt' % (ring, lam)))
        print('λ = %g:' % lam)
        describe(r, '  ')
        top_lines(r, thr=-5.0)
        ok = np.isfinite(r.D)
        print('  D の最大 %.2f、D > -3: %d 点' % (np.nanmax(r.D), (ok & (r.D > -3)).sum()))


if __name__ == '__main__':
    a = [t for t in sys.argv[1:] if not t.startswith('--')]
    w = 1
    if '--workers' in sys.argv:
        w = int(sys.argv[sys.argv.index('--workers') + 1])
        a = [t for t in a if t != str(w)]
    mode = a[0]
    if mode == 'k64':
        mode_k64(a[1], w)
    elif mode == 'errors':
        mode_errors(a[1], a[2], w)
    elif mode == 'mpw':
        mode_mpw(w)
    elif mode == 'delta':
        mode_delta(a[1], a[2], w)
    elif mode == 'naff':
        mode_naff(a[1], w)
    elif mode == 'grid':
        mode_grid(a[1], a[2], a[3], w)
    elif mode == 'lambda':
        mode_lambda(a[1], w)
    else:
        raise SystemExit(__doc__)
