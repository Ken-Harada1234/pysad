# -*- coding: utf-8 -*-
"""周波数マップの D の大きい領域を、1 次の RDT だけ (fmap_rdt.py。出力の表記 S1) と 2 次の RDT (4 次の h_jklm、6 極の 2 次など) も入れた共鳴の領域の予想で比べる。
使い方 (リポジトリ直下から): py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_rdt4.py astrid2|kekls|fodo
周波数マップは保存済みの表を読む (astrid2・kekls: fmap_<ring>.table.txt、fodo: fmap_rdt_fodo.table.txt。1024 周、v2、RFSW・RAD OFF、LOSSAMPL 0.1 m)。
予想 (fmapres.ResonancePredictor): V = 2|F|、F = Σ h (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) (調和 a = (j−k, l−m) ごと)。1 次だけ (S1) は 3 次の h だけ、2 次も (R5) は 4 次の h (map の f4) も同じ調和に足す。
  振り子の式 (位相つき)・包絡の定義、ν(J) = 線形の ADTS (または追跡のチューン)、対照 = 領域を格子の上で巡回ずらし (400 回)、は fmap_rdt.py と同じ。
出力: fmap_rdt4_<ring>.out.txt、fmap_rdt4_<ring>.png"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import fmap_rdt as S1

OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(s)


def main(name):
    from pysad import compat
    from pysad import fmapres as FR
    from pysad.freqmap import FMapResult
    from pysad.parser import load_sad_string
    compat.set_mode('v2')
    cs = S1.CASES[name]
    lat = load_sad_string(S1.source(name)).build_lattice(cs['line'])
    rows = lat.twiss_table()
    loc = cs['loc']
    tab = cs['table'] or os.path.join(HERE, 'fmap_rdt_%s.table.txt' % name)
    r = FMapResult.load_table(tab)
    nx, ny, _ = r.grid
    ok = np.isfinite(r.D)
    nu0 = (rows[-1].mux / 2 / math.pi, rows[-1].muy / 2 / math.pi)
    P('R5 周波数マップの共鳴線の予想: %s' % cs['label'])
    P('作動点 (CALC) νx=%.6f νy=%.6f。周波数マップ: %s (格子 %s、nturn %s)。点 %d、D > -5: %d 点、D > -4: %d 点' % (
        nu0[0], nu0[1], os.path.basename(tab), r.cond.get('grid_nx_ny_ndp'), r.cond.get('nturn'), len(r), (ok & (r.D > -5)).sum(), (ok & (r.D > -4)).sum()))
    pr1 = FR.ResonancePredictor.from_lattice(lat, loc, rows)
    pr2 = FR.ResonancePredictor.from_lattice(lat, loc, rows, order2=True)
    P('1 次の調和 %s、2 次を足した調和 %d 個' % (sorted(pr1.harms), len(pr2.harms)))
    P('調和ごとの 4 次の RDT の最大 |h| と、線形の作動点からの共鳴線までの距離 (チューン平面):')
    for a in sorted(pr2.harms):
        t4 = [t for t in pr2.harms[a] if t[0] + t[1] + t[2] + t[3] == 4]
        if not t4:
            continue
        v = a[0] * nu0[0] + a[1] * nu0[1]
        P('   調和 %-8s 共鳴 %d νx %+d νy = %d: 距離 %.4f、最大 |h4| = %.3g (%s)' % (str(a), a[0], a[1], round(v), abs(v - round(v)) / math.hypot(*a), max(abs(t[4]) for t in t4), 'h%d%d%d%d0' % max(t4, key=lambda t: abs(t[4]))[:4]))
    ip = np.floor(nu0[0]), np.floor(nu0[1])
    nu_meas = (r.nux1 + ip[0], r.nuy1 + ip[1])
    P('')
    P('D > thr の点が予想の領域に入る割合 (cov)、領域の面積 (area)、cov/area、対照 (領域を格子の上でずらした 400 回) の cov の平均 ± 標準偏差 (最大)、p 値')
    rows_out = []
    for vname, nu in (('線形の ADTS の ν(J)', None), ('追跡のチューン ν', nu_meas)):
        for lab, pr in (('S1 (1 次のみ)', pr1), ('R5 (2 次も)', pr2)):
            Pm = pr.predict(r.x, r.y, nu=nu)
            for key, kn in (('any', '位相つき'), ('envelope', '包絡')):
                for thr in (-5.0, -4.0):
                    s = FR.coverage_stats(r.D, Pm[key], r.grid, thr=thr)
                    P('  [%-14s %-12s %-8s D>%g] D>thr %d 点: area %.3f、cov %.3f、cov/area %.2f、対照 %.3f ± %.3f (最大 %.3f)、p = %.3f' % (
                        vname, lab, kn, thr, s['n_hi'], s['area'], s['cov'], s['enrich'], s['shift_mean'], s['shift_std'], s['shift_max'], s['pval']))
                    rows_out.append((vname, lab, key, thr, s))
    # 調和ごと (R5、追跡のチューン、包絡)
    P('')
    P('調和ごとの寄与 (R5、追跡のチューン、包絡、D > -5): 調和 / 領域の点数 / その調和だけの cov (area)、対照の平均')
    Pm = pr2.predict(r.x, r.y, nu=nu_meas)
    for a, v in sorted(Pm['by_harmonic'].items()):
        if not v['envelope'].any():
            continue
        s = FR.coverage_stats(r.D, v['envelope'], r.grid, thr=-5.0, nshift=200)
        P('   %-8s %4d 点、cov %.3f (area %.3f)、対照 %.3f ± %.3f、G0 = %.0f、W 中央値 %.3g' % (str(a), v['envelope'].sum(), s['cov'], s['area'], s['shift_mean'], s['shift_std'], v['G0'], np.nanmedian(v['W'])))
    # 幅の係数: 領域 |Δ0| < c·W (W = 振り子の式の幅) の c を変えたときの、2 次の調和 (R5、追跡のチューン) の領域の統計。c < 1 で面積が小さくなるときの enrich を見る
    P('')
    P('幅の係数 c (領域 |Δ0| < c W、R5 の全調和の和、追跡のチューン、D > -5):')
    for c in (1.0, 0.5, 0.25, 0.125, 0.0625):
        m = np.zeros(r.x.shape, dtype=bool)
        for a, v in Pm['by_harmonic'].items():
            with np.errstate(invalid='ignore'):
                m |= v['envelope'] & (np.abs(v['delta0']) < c * v['W'])
        s = FR.coverage_stats(r.D, m, r.grid, thr=-5.0)
        P('   c = %-6g area %.3f、cov %.3f、cov/area %.2f、対照 %.3f ± %.3f、p = %.3f' % (c, s['area'], s['cov'], s['enrich'], s['shift_mean'], s['shift_std'], s['pval']))
    # 最も近い共鳴線 (4 次以下) の内訳
    oth = FR.other_mn(4)
    d1 = FR.lines_distance(r.nux1, r.nuy1, FR.FIRST_ORDER_MN)
    d2 = FR.lines_distance(r.nux1, r.nuy1, oth)
    for thr in (-5.0, -4.0):
        sel = ok & (r.D > thr)
        n = int(sel.sum())
        P('D > %g の %d 点: 最も近い共鳴線 (距離 < 0.003) は 1 次の RDT の調和 %d、それ以外の 4 次以下 %d、どちらにも近くない %d' % (
            thr, n, (sel & (d1 < 0.003) & (d1 <= d2)).sum(), (sel & (d2 < 0.003) & (d2 < d1)).sum(), (sel & (d1 >= 0.003) & (d2 >= 0.003)).sum()))
    # 4 次以下の各共鳴線の「その線に最も近い D > -5 の点」の数と、R5 の予想の領域に入る数
    names = {}
    from pysad.freqmap import _mn_list
    for mn in _mn_list(4):
        v = mn[0] * r.nux1 + mn[1] * r.nuy1
        dist = np.abs(v - np.round(v)) / math.hypot(*mn)
        sel = ok & (r.D > -5) & (dist < 0.003)
        if sel.sum():
            names[mn] = int(sel.sum())
    P('D > -5 の点のうち、各共鳴線 (m, n) (次数 4 以下) の 0.003 以内に入る点の数 (1 点が複数の線に入りうる): %s' % ', '.join('%s: %d' % (k, v) for k, v in sorted(names.items(), key=lambda kv: -kv[1])[:10]))
    # より高次の共鳴: 測ったチューンが次数 N 以下のどれかの線の 0.003 以内にある点の割合 (D > -5 の点 / 全点)
    P('測ったチューンが次数 N 以下の共鳴線 (どの (m, n)) の 0.003 以内にある点の割合: D > -5 の点 / 全点 (比 = 線への近さの D への効き方)')
    from pysad.freqmap import _mn_list
    for N in (3, 4, 5, 6, 8):
        d = FR.lines_distance(r.nux1, r.nuy1, _mn_list(N))
        near = ok & (d < 0.003)
        hi = ok & (r.D > -5)
        P('   N = %d: %.3f / %.3f (比 %.2f)' % (N, (near & hi).sum() / max(hi.sum(), 1), near.sum() / max(ok.sum(), 1), ((near & hi).sum() / max(hi.sum(), 1)) / max(near.sum() / max(ok.sum(), 1), 1e-9)))
    FR.plot_with_prediction(r, [('S1: order 1 (measured nu, envelope)', pr1.predict(r.x, r.y, nu=nu_meas)['envelope'], 'white', '-'),
                                ('R5: + order 2 (measured nu, envelope)', pr2.predict(r.x, r.y, nu=nu_meas)['envelope'], 'magenta', '--')],
                            os.path.join(HERE, 'fmap_rdt4_%s.png' % name), title='%s: D and predicted resonance regions (order 1 vs 1+2)' % name)
    with open(os.path.join(HERE, 'fmap_rdt4_%s.out.txt' % name), 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2')
