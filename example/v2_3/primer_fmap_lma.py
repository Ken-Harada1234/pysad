# -*- coding: utf-8 -*-
"""入門 §6(c)(d)・§8(4): ASTRID2 (0.58 GeV、45.704 m、example/astrid/lattice/ASTRID2_v01_01.sad) の周波数マップと局所運動量アクセプタンス (LMA) を、
入門用に小さく (粗い格子・少ない位置で) 走らせる。条件と手順は docs/48 と同じで、点の数だけ少ない (Q1_fmap_ring.py・Q1_lma_ring.py は 41×31 と 25 位置)。
  scan : 位置 2 (長直線部の入口、RF が限界) と位置 17 (分散の大きい位置、アパーチャが限界) で δ = 0 … +2% を 0.1% 刻みに追跡し、落ちた周を図にする。
         LMA の二分法がこの図の「最初に落ちる δ」を探していることを見る。放射・RF つき、LOSSAMPL = ±10 mm、バケット判定、1312 周 (= νs の 2 周期分)。
  lma  : 最初のセルの位置 (1 m 以上の間隔) の δ±(s) (二分法 tol = 1e-4、check_every = 1、2 並列)
  fmap : 位置 2、x = 0.5 … 18 mm (21 点) × y = 0.5 … 14 mm (16 点)、1024 周、RFSW・RAD OFF、δ = 0、NAFF Hann p=1。表と図
使い方 (リポジトリ直下から。PYTHONUTF8=1 が日本語の出力のリダイレクトに要る):
  py -3.12 example/v2_3/primer_fmap_lma.py scan|lma|fmap|all      出力: primer_lma_scan.png、primer_lma_astrid2.png、primer_fmap_astrid2.png (+ .table.txt)
"""
import os
import sys
import time
import numpy as np
import primer_common as C
plt = C.setup_plot()

SRC = os.path.join(C.ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')
NTURN_FMAP = 1024


def part_scan():
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    src = open(SRC, encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.010)
    st = st.copy(check_every=1)
    zb = st.bucket_zlimit()
    nt = st.nturn_for_synchrotron(2.0)
    print('ASTRID2 0.58 GeV、RFSW ON、RAD ON、FLUC なし、LOSSAMPL = ±10 mm、check_every = 1、バケット判定 |z| < %.4f m、周回数 %d (= νs の 2 周期分、νs = %.5f)、δ_rf = %.4f%%' % (
        zb, nt, 2.0 / nt, 100 * st.bucket['height']))
    dps = np.arange(0.0, 0.0201, 0.001)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for a, loc in zip(axs, (2, 17)):
        cod = st.closed_orbit([loc])[0]
        c = np.tile(cod, (len(dps), 1))
        c[:, 5] += dps
        r = st.track(loc, c, nt, zlost=zb)
        lt = np.where(r['lost_turn'] == 0, nt + 1, r['lost_turn'])
        first = next((d for d, l in zip(dps, r['lost_turn']) if l > 0), None)
        if first is not None:
            k = [i for i, l in enumerate(r['lost_turn']) if l > 0][0]
            fx = r['final'][k]
            print('  最初に落ちた粒子 (δ = %.1f%%): 落ちた周 %d、落ちた時点の x = %.2f mm、z = %.3f m (バケットの限界 |z| = %.3f m、アパーチャ ±10 mm)' % (
                100 * dps[k], r['lost_turn'][k], 1e3 * fx[0], fx[4], zb))
        print('位置 %d (s = %.3f m, %s): 最初に落ちる δ = %s。δ ごとの落ちた周 (0 = 残った): %s' % (
            loc, st.s_of(loc), st.names[loc - 1], 'なし (2% まで残る)' if first is None else '%.1f%%' % (100 * first),
            ' '.join('%d' % l for l in r['lost_turn'])))
        a.plot(dps * 100, lt, 'o-', ms=4)
        a.axhline(nt, color='gray', ls=':')
        a.set_xlabel('δ [%]'); a.set_ylabel('落ちた周 (残った粒子は %d + 1 で表示)' % nt)
        a.set_title('位置 %d (%s): δ の走査' % (loc, '長直線部の入口' if loc == 2 else '分散の大きい位置'), fontsize=10)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(C.HERE, 'primer_lma_scan.png'), dpi=110)


def part_lma(workers=2):
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity, select_locations
    from pysad.lma import momentum_acceptance
    compat.set_mode('v2')
    src = open(SRC, encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=0.010)
    st = st.copy(check_every=1)
    per = find_periodicity(st.els, st.rows)
    print('周期: セルの数 %d、セルの長さ %.4f m、全部同じ: %s' % (per['nper'], per['cell_length'], per['ok']))
    locs = select_locations(st, per['starts'][0], per['starts'][1] - 1, 1, 1.0)
    print('位置の数 %d: %s' % (len(locs), locs))
    t0 = time.time()
    r = momentum_acceptance(st, locs, workers=workers)
    print('LMA: %.0f s、追跡した粒子 %d、周回数 %d、バケット δ_rf = %.4f%%' % (time.time() - t0, r.cond['tracked_particles'], r.cond['nturn'], 100 * st.bucket['height']))
    print('%5s %9s %-8s %9s %9s  %s' % ('loc', 's[m]', 'name', 'dp+[%]', 'dp-[%]', '損失の理由 (+, −)'))
    for i in range(len(r.loc)):
        print('%5d %9.3f %-8s %9.4f %9.4f  %s, %s' % (r.loc[i], r.s[i], r.name[i], 100 * r.dp_plus[i], 100 * r.dp_minus[i], r.reason_plus[i], r.reason_minus[i]))
    r.save(os.path.join(C.HERE, 'primer_lma_astrid2.table.txt'))
    r.plot(os.path.join(C.HERE, 'primer_lma_astrid2.png'))


def part_fmap():
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.freqmap import frequency_map, resonance_check
    compat.set_mode('v2')
    src = open(SRC, encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=False, rfsw=False, alost=0.1)
    nux, nuy = st.tunes_full()
    print('ASTRID2 作動点 νx = %.6f、νy = %.6f' % (nux, nuy))
    for (m, n, p, d, o, nm) in resonance_check(nux, nuy, 3):
        print('  %s: %d νx + %d νy = %d (次数 %d)、距離 %.4f' % (nm, m, n, p, o, d))
    t0 = time.time()
    r = frequency_map(st, 2, x=np.linspace(5e-4, 18e-3, 21), y=np.linspace(5e-4, 14e-3, 16), nturn=NTURN_FMAP, workers=1)
    print('周波数マップ: %d 点 × %d 周、%.0f s' % (len(r.D), NTURN_FMAP, time.time() - t0))
    ok = np.isfinite(r.D)
    print('失われた点 %d、D の範囲 %.2f 〜 %.2f、D > −4 の点 %d、D < −6 の点 %d (規則的)' % ((~ok).sum(), np.nanmin(r.D), np.nanmax(r.D), (r.D > -4).sum(), (r.D < -6).sum()))
    r.save(os.path.join(C.HERE, 'primer_fmap_astrid2.table.txt'))
    r.plot(os.path.join(C.HERE, 'primer_fmap_astrid2.png'), order=4)


if __name__ == '__main__':
    w = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if w in ('scan', 'all'):
        part_scan()
    if w in ('lma', 'all'):
        part_lma()
    if w in ('fmap', 'all'):
        part_fmap()
