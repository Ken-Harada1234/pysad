# -*- coding: utf-8 -*-
"""v2.3.0 S1 (a): 周波数マップの D の大きい領域を、共鳴駆動項 (RDT) と ADTS から予想した共鳴の領域と定量的に比べる。
使い方 (リポジトリ直下から): py -3.12 example/v2_3/S1_fmap_rdt.py astrid2|kekls|fodo [--recompute] [--workers N]
  astrid2・kekls は Q1 の周波数マップの表 (Q1_fmap_<ring>.table.txt、1024 周) を読む (--recompute で追跡し直して S1_ の名前で保存)。
  fodo は新しく追跡する (24 セルの FODO、QF の K1 = 0.697 で νx = 20/3 の少し下、SF・SD の K2 = 4, −6、3 GeV、位置 1 から、1024 周)。
予想 (pysad/fmapres.py): 6 極 (と偏向磁石の 3 次) の 1 次の RDT h_jklm (resonance_terms) と ADTS (adts().total) の行列 A から、各調和 (3νx、νx、νx ± 2νy) の
  共鳴 m νx + n νy = p について、初期条件 (x, y) (px = py = 0) の作用 J から、その振幅で共鳴の島に捕らえられるか (振り子の式、位相つき)、位相を問わない包絡に入るか、を計算し、
  D のマップに重ねる。D > −5 (または −4) の点が予想の領域に入る割合を、領域を格子の上でずらした (同じ形・同じ面積の) 対照と比べる。
  ν(J) は (1) 線形の ADTS ν0 + A·J と (2) 追跡で測ったチューン (NAFF の前半のチューン) の 2 通り。
出力: S1_fmap_rdt_<ring>.out.txt (標準出力)、.png、(fodo は .table.txt)。
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), line='RING', loc=2, x=(0.0005, 0.018), y=(0.0005, 0.014),
                    table=os.path.join(HERE, 'Q1_fmap_astrid2.table.txt'),
                    label='ASTRID2 (0.58 GeV、周長 45.7 m)、位置 2 (長直線部の入口、βx=3.89 m、βy=1.62 m)'),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), line='RING', loc=3, x=(0.0002, 0.006), y=(0.0002, 0.008),
                  table=os.path.join(HERE, 'Q1_fmap_kekls.table.txt'),
                  label='kek-ls_v03_12.lat (3 GeV、周長 603 m)、位置 3 (長直線部の入口、βx=4.64 m、βy=6.23 m)'),
    'fodo': dict(path=os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), line='RING', loc=1, x=(0.0005, 0.014), y=(0.0005, 0.008), table=None,
                 label='FODO (24 セル、3 GeV、周長 105.6 m、QF の K1 = 0.697、SF・SD の K2 = 4, −6)、位置 1 (リングの先頭、βx=5.03 m、αx=−2.17)'),
}
NX, NY, NTURN = 41, 31, 1024


def source(name):
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    if name == 'fodo':
        src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)').replace('QF = (L=0.3 K1=0.7)', 'QF = (L=0.3 K1=0.697)')
    return src


def main(name, recompute=False, workers=1):
    from pysad import compat
    from pysad import fmapres as FR
    from pysad.freqmap import FMapResult, frequency_map, resonance_check
    from pysad.parser import load_sad_string
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    cs = CASES[name]
    src = source(name)
    lat = load_sad_string(src).build_lattice(cs['line'])
    rows = lat.twiss_table()
    loc = cs['loc']
    print(cs['label'])
    nu0 = (rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi)
    print('作動点 (CALC): nu_x = %.6f, nu_y = %.6f' % nu0)
    tab = cs['table'] or os.path.join(HERE, 'S1_fmap_rdt_%s.table.txt' % name)
    if os.path.exists(tab) and not recompute:
        r = FMapResult.load_table(tab)
        print('周波数マップ: 保存済みの表 %s を読む (格子 %s、nturn %s、失われた点 %s)' % (os.path.basename(tab), r.cond.get('grid_nx_ny_ndp'), r.cond.get('nturn'), r.cond.get('nlost')))
    else:
        st = TrackSetup.from_source(src, cs['line'], rad=False, rfsw=False, alost=0.1)[0]
        t0 = time.time()
        r = frequency_map(st, loc, x=np.linspace(*cs['x'], NX), y=np.linspace(*cs['y'], NY), nturn=NTURN, workers=workers, progress=0)
        print('周波数マップ: 追跡 %.0f s、条件 %s' % (time.time() - t0, r.cond))
        r.save(os.path.join(HERE, 'S1_fmap_rdt_%s.table.txt' % name))
    nx, ny, _ = r.grid
    ok = np.isfinite(r.D)
    print('点 %d、失われた点 %d、D の範囲 %.2f 〜 %.2f、D > -5: %d 点、D > -4: %d 点' % (len(r), (r.lost_turn > 0).sum(), np.nanmin(r.D), np.nanmax(r.D),
                                                                      (ok & (r.D > -5)).sum(), (ok & (r.D > -4)).sum()))
    pr = FR.ResonancePredictor.from_lattice(lat, loc, rows)
    A = pr.A
    print('ADTS (解析、全部入り) [1/(m rad)]: xx = %.1f, xy = %.1f, yy = %.1f' % (A[0, 0], A[0, 1], A[1, 1]))
    print('1 次の RDT |h| (調和 (m, n) = (j−k, l−m) ごと、単項式の |h|、h は resonance_terms の値) と、線形の作動点からの共鳴線 m νx + n νy = p までの距離 (チューン平面):')
    for a, terms in sorted(pr.harms.items()):
        v = a[0] * nu0[0] + a[1] * nu0[1]
        d = abs(v - round(v)) / math.hypot(*a)
        print('   調和 %-8s 共鳴 %d νx %+d νy = %d: 距離 %.4f、%s' % (str(a), a[0], a[1], round(v), d,
                                                              ', '.join('h%d%d%d%d0 = %.3f' % (q[0], q[1], q[2], q[3], abs(q[4])) for q in terms)))
    # 追跡のチューンと線形の ADTS の差
    ip = np.floor(nu0[0]), np.floor(nu0[1])
    nu_meas = (r.nux1 + ip[0], r.nuy1 + ip[1])
    Jx, Jy, _px, _py = pr.actions(r.x, r.y)
    nl = pr.tunes_lin(Jx, Jy)
    ex, ey = nl[0] - nu_meas[0], nl[1] - nu_meas[1]
    good = ok & (np.abs(ex) < 0.003) & (np.abs(ey) < 0.003)
    srt = np.argsort(Jx + Jy)
    third = [i for i in srt if ok[i]][:int(ok.sum() / 3)]
    print('線形の ADTS ν0 + A·J と、追跡 (NAFF) のチューンの差: 全点の中央値 |Δνx| = %.4f、|Δνy| = %.4f、(0.003 以内の点) %d / %d、'
          '小振幅の 1/3 (Jx + Jy が小さい順) の中央値 %.4f・%.4f、小振幅の 1/3 で 0.003 以内 %.2f' % (
              np.median(np.abs(ex[ok])), np.median(np.abs(ey[ok])), good.sum(), ok.sum(), np.median(np.abs(ex[third])), np.median(np.abs(ey[third])),
              np.mean((np.abs(ex[third]) < 0.003) & (np.abs(ey[third]) < 0.003))))
    Pl = pr.predict(r.x, r.y)
    Pm = pr.predict(r.x, r.y, nu=nu_meas)
    print('共鳴の予想の領域に D > thr の点が入る割合 (cov)、領域の面積の割合 (area)、cov/area、領域の中・外の D > thr の割合、'
          '領域を格子の上で巡回ずらしした対照 (400 回) の cov の平均 ± 標準偏差・最大、p 値 (ずらしの cov ≥ cov の割合):')
    for vname, P in (('線形の ADTS の ν(J)', Pl), ('追跡のチューン ν', Pm)):
        for key, kn in (('any', '位相つきで捕らえられる'), ('envelope', '包絡 (位相を問わない)')):
            for thr in (-5.0, -4.0):
                s = FR.coverage_stats(r.D, P[key], r.grid, thr=thr)
                print('  [%s / %s / D > %g] 点 %d (D>thr %d): area %.3f、cov %.3f、cov/area %.2f、中 %.3f / 外 %.3f、対照 %.3f ± %.3f (最大 %.3f)、p = %.3f' % (
                    vname, kn, thr, s['n_ok'], s['n_hi'], s['area'], s['cov'], s['enrich'], s['p_in'], s['p_out'], s['shift_mean'], s['shift_std'], s['shift_max'], s['pval']))
    print('調和ごとの予想 (追跡のチューンを使う。位相つき / 包絡の点の数、G0 = (m0, n0)·A·(m0, n0) [1/(m rad)]、幅 W の中央値 [チューン]):')
    for a, v in sorted(Pm['by_harmonic'].items()):
        s = FR.coverage_stats(r.D, v['trapped'], r.grid, thr=-5.0, nshift=100)
        print('   調和 %-8s: %4d / %4d 点、G0 = %.0f、W の中央値 %.3g、この調和だけの D > -5 の cov %.3f (area %.3f)' % (
            str(a), v['trapped'].sum(), v['envelope'].sum(), v['G0'], np.nanmedian(v['W']), s['cov'], s['area']))
    # D の大きい点の最も近い共鳴線が、1 次の RDT が駆動する調和か、そうでないか
    oth = FR.other_mn(4)
    nuf = (r.nux1, r.nuy1)
    d1 = FR.lines_distance(nuf[0], nuf[1], FR.FIRST_ORDER_MN)
    d2 = FR.lines_distance(nuf[0], nuf[1], oth)
    for thr in (-5.0, -4.0):
        sel = ok & (r.D > thr)
        n = int(sel.sum())
        near1 = sel & (d1 < 0.003) & (d1 <= d2)
        near2 = sel & (d2 < 0.003) & (d2 < d1)
        none = sel & (d1 >= 0.003) & (d2 >= 0.003)
        print('D > %g の %d 点の、最も近い共鳴線 (チューン平面の距離 < 0.003): 1 次の RDT の調和 (3νx・νx・νx±2νy) %d (%.2f)、それ以外の 4 次以下 (νx−νy、4 次など) %d (%.2f)、どちらにも近くない %d (%.2f)' % (
            thr, n, near1.sum(), near1.sum() / max(n, 1), near2.sum(), near2.sum() / max(n, 1), none.sum(), none.sum() / max(n, 1)))
    a1 = (d1 < 0.003)[ok].mean()
    a2 = (d2 < 0.003)[ok].mean()
    print('参考: 全点のうち 1 次の調和の線から 0.003 以内 %.3f、それ以外の 4 次以下の線から 0.003 以内 %.3f' % (a1, a2))
    print('y の列ごとの、D > -5 の点の x の範囲と、予想 (追跡のチューン、位相つき・包絡) の x の範囲 [mm] (D > -5 の点がある列だけ):')
    X2, Y2 = r.x.reshape(nx, ny), r.y.reshape(nx, ny)
    D2, T2, E2 = r.D.reshape(nx, ny), Pm['any'].reshape(nx, ny), Pm['envelope'].reshape(nx, ny)
    for j in range(ny):
        hi = np.isfinite(D2[:, j]) & (D2[:, j] > -5)
        if not hi.any():
            continue
        def rng(m):
            return '%.2f〜%.2f' % (X2[m, j].min() * 1e3, X2[m, j].max() * 1e3) if m.any() else '(なし)'
        print('   y = %5.2f mm: D>-5 %2d 点 x = %s、予想 位相つき %s、包絡 %s' % (Y2[0, j] * 1e3, hi.sum(), rng(hi), rng(T2[:, j]), rng(E2[:, j])))
    # 3νx の共鳴の中心 (Δ0 = 0) の x と、D の尾根 (列ごとの D の最大の x) の比較
    if (3, 0) in pr.harms and Pm['by_harmonic'][(3, 0)]['trapped'].any():
        d0l = Pl['by_harmonic'][(3, 0)]['delta0'].reshape(nx, ny)
        d0m = Pm['by_harmonic'][(3, 0)]['delta0'].reshape(nx, ny)
        diffs_l, diffs_m = [], []
        print('3νx の共鳴の中心 x(Δ0 = 0) (線形の ADTS / 追跡のチューン) と D の尾根 (列の D の最大の x) [mm] (D > -5 の点がある列):')
        for j in range(ny):
            col = D2[:, j]
            if not (np.isfinite(col) & (col > -5)).any():
                continue
            xs_ = X2[:, j] * 1e3

            def zero(d0):
                idx = np.nonzero(np.sign(d0[:-1, j]) * np.sign(d0[1:, j]) < 0)[0]
                if len(idx) == 0:
                    return float('nan')
                i = idx[0]
                return xs_[i] + (xs_[i + 1] - xs_[i]) * d0[i, j] / (d0[i, j] - d0[i + 1, j])
            zl, zm = zero(d0l), zero(d0m)
            ridge = xs_[np.nanargmax(np.where(np.isfinite(col), col, -99))]
            diffs_l.append(zl - ridge)
            diffs_m.append(zm - ridge)
            if j % 3 == 0:
                print('   y = %5.2f mm: 線形 ADTS %.2f、追跡のチューン %.2f、尾根 %.2f' % (Y2[0, j] * 1e3, zl, zm, ridge))
        dl, dm = np.array(diffs_l), np.array(diffs_m)
        print('   中心 − 尾根の平均 ± 標準偏差 [mm]: 線形 ADTS %.2f ± %.2f、追跡のチューン %.2f ± %.2f (%d 列、x の格子の間隔 %.2f mm)' % (
            np.nanmean(dl), np.nanstd(dl), np.nanmean(dm), np.nanstd(dm), len(dl), (X2[1, 0] - X2[0, 0]) * 1e3))
    FR.plot_with_prediction(r, [('trapped (linear ADTS nu(J))', Pl['any'], 'white', '-'), ('trapped (measured nu)', Pm['any'], 'magenta', '--'),
                                ('envelope (measured nu)', Pm['envelope'], 'k', ':')],
                            os.path.join(HERE, 'S1_fmap_rdt_%s.png' % name), title='%s: D (colour) and predicted first-order resonance regions (contours)' % name)


if __name__ == '__main__':
    nm = sys.argv[1] if len(sys.argv) > 1 else 'astrid2'
    w = 1
    if '--workers' in sys.argv:
        w = int(sys.argv[sys.argv.index('--workers') + 1])
    main(nm, '--recompute' in sys.argv, w)
