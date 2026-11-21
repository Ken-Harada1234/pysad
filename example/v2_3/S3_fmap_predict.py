# -*- coding: utf-8 -*-
"""v2.3.0 S3: 周波数マップの共鳴線の予想の作り直し (pysad/fmappred.py) の判別力の評価と、改善しない原因の切り分け。
使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/v2_3/S3_fmap_predict.py astrid2|kekls|fodo|mpw
周波数マップは保存済みの表 (astrid2・kekls: Q1_fmap_*.table.txt、fodo: S1_fmap_rdt_fodo.table.txt、mpw: S1_fmap_mpw.table.txt。1024 周、v2、RFSW・RAD OFF、LOSSAMPL 0.1 m、41×31) を読む。
予想: S1 (3 次の h)、R5 (+ 4 次の h の写像の係数)、S3 (+ 非共鳴の 3 次の摂動の補正 F4' = rdt4_normal_form)。ν(J) は解析の線形の ADTS (adts(fourth='r3', optics='calc'))。
  当てはめの項・追跡に合わせた補正は入れない。追跡のチューンを使う行 (診断) は ν(J) の誤差の効果を測るもので、予想ではない。
出力: S3_fmap_predict_<ring>.out.txt、.png"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import S1_fmap_rdt as S1

S1.CASES['mpw'] = dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), line='RINGMPW', loc=2,
                       table=os.path.join(HERE, 'S1_fmap_mpw.table.txt'), label='ASTRID2 + MPW (ASTRID2_v01_02_id.sad の RINGMPW、0.58 GeV、45.7 m)')
OUT = []


def P(s=''):
    print(s, flush=True)
    OUT.append(s)


def fmt(e):
    c, s = e['cov'], e['score']
    return 'area %.3f cov %.3f 対照 %.3f±%.3f p=%.3f | AUC %.3f ρ %.3f 上位10%%の D>-5 率 %.2f (全体 %.2f) D 中央値 %.2f (全体 %.2f) | 振幅だけの AUC %.3f、偏相関 %.3f' % (
        c['area'], c['cov'], c['shift_mean'], c['shift_std'], c['pval'], s['auc'], s['rho'], s['top_frac_pos'], s['base_frac_pos'], s['top_Dmed'], s['base_Dmed'], s['amp_auc'], s['rho_partial'])


def _red(a):
    g = math.gcd(abs(a[0]), abs(a[1]))
    m0 = (a[0] // g, a[1] // g)
    if m0[0] < 0 or (m0[0] == 0 and m0[1] < 0):
        m0 = (-m0[0], -m0[1])
    return m0


def s1_score(mm):
    """S1・R5 の予想 (連続のスコアが無い) から、包絡の中の log10(W/|Δ0|) の最大を作る (比較用)"""
    sc = np.full(np.shape(mm['envelope']), -1.0)
    for a, v in mm['by_harmonic'].items():
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.log10(v['W'] / np.maximum(np.abs(v['delta0']), 1e-300))
        t = np.where(np.isfinite(t), t, -1.0)
        sc = np.maximum(sc, np.where(v['envelope'], t, -1.0))
    return sc


def main(name):
    from pysad import compat
    from pysad import fmapres as FR
    from pysad import fmappred as FP
    from pysad.freqmap import FMapResult, _mn_list
    from pysad.parser import load_sad_string
    compat.set_mode('v2')
    cs = S1.CASES[name]
    src = open(cs['path'], encoding='utf-8').read() if name == 'mpw' else S1.source(name)
    lat = load_sad_string(src).build_lattice(cs['line'])
    rows = lat.twiss_table()
    loc = cs['loc']
    tab = cs['table'] or os.path.join(HERE, 'S1_fmap_rdt_%s.table.txt' % name)
    r = FMapResult.load_table(tab)
    ok = np.isfinite(r.D)
    P('S3 周波数マップの共鳴線の予想: %s' % cs['label'])
    P('周波数マップ: %s (格子 %s、nturn %s)。点 %d、失われた点 %d、D > -5: %d 点、D > -4: %d 点' % (
        os.path.basename(tab), r.cond.get('grid_nx_ny_ndp'), r.cond.get('nturn'), len(r), (r.lost_turn > 0).sum(), (ok & (r.D > -5)).sum(), (ok & (r.D > -4)).sum()))
    s3 = FP.ResonancePredictorS3.from_lattice(lat, loc, rows)
    s1 = FR.ResonancePredictor.from_lattice(lat, loc, rows)
    r5 = FR.ResonancePredictor.from_lattice(lat, loc, rows, order2=True)
    ip = (np.floor(s3.nu0[0]), np.floor(s3.nu0[1]))
    nm = (r.nux1 + ip[0], r.nuy1 + ip[1])
    P('作動点 νx=%.5f νy=%.5f、ADTS xx=%.1f xy=%.1f yy=%.1f [1/(m rad)]、1 次の RDT の共鳴線までの最短距離 %.4f' % (
        s3.nu0[0], s3.nu0[1], s3.A[0, 0], s3.A[0, 1], s3.A[1, 1], s3.rowsinfo['d_first_min']))
    P('駆動される調和 (3 次 + 4 次): %s' % sorted(s3.harms))
    P('')
    P('(1) 判別力 (D > -5 を正例。包絡 |Δ0| < W の領域、スコア = log10(W/|Δ0|) の調和にわたる最大)')
    Jx, Jy, _a, _b = s3.actions(r.x, r.y)
    amp = np.sqrt(Jx / Jx.max() + Jy / Jy.max())
    res = {}
    for vname, nu in (('解析 ν(J) (予想)', None), ('追跡のチューン (診断)', nm)):
        for lab, pr, kw in (('S1  (3 次)', s1, {}), ('R5  (+4 次の写像の係数)', r5, {}), ('S3  (+正規形の補正)', s3, dict(nf=True)), ('S3  補正なし', s3, dict(nf=False))):
            if pr is s3:
                m = pr.predict(r.x, r.y, nu=nu, **kw)
            else:
                mm = pr.predict(r.x, r.y, nu=nu)
                m = dict(envelope=mm['envelope'], any=mm['any'], score=s1_score(mm), Jx=Jx, Jy=Jy)
            e = FP.evaluate_prediction(r, m, amp=amp, nshift=300)
            res[(vname, lab)] = (m, e)
            P('  [%-18s %-22s] %s' % (vname, lab, fmt(e)))
    P('  (ずらしの対照 = 領域を格子の上で巡回ずらし 300 回。AUC 0.5 = 判別なし。振幅だけのスコアの AUC は、D が振幅とともに増えることによる見かけの判別力の目安。偏相関 = スコアと D の順位相関から振幅の順位の効果を除いたもの)')
    base = FP.line_distance_score(nm[0], nm[1], 4)
    ee = FP.score_stats(r.D, base, amp=amp)
    P('  比較 (予想の強さを使わない): 測ったチューンの次数 4 以下の共鳴線への近さ −log10(距離): AUC %.3f、ρ %.3f、偏相関 %.3f' % (ee['auc'], ee['rho'], ee['rho_partial']))
    m3 = res[('解析 ν(J) (予想)', 'S3  (+正規形の補正)')][0]
    P('')
    P('(2) 信頼度のフラグ (S3、解析 ν(J)): 点の数 / D > -5 の割合 / D の中央値 / スコアの AUC (クラスの中)')
    rel = m3['reliability']
    for k, nmk in ((2, '信頼できる'), (1, '注意'), (0, '信頼できない')):
        sel = ok & (rel == k)
        if sel.sum() == 0:
            P('   %-8s 0 点' % nmk)
            continue
        pos = r.D[sel] > -5
        P('   %-8s %4d 点、D > -5 の割合 %.3f、D の中央値 %.2f、AUC %.3f (正例 %d)' % (nmk, sel.sum(), pos.mean(), np.median(r.D[sel]), FP.auc_score(m3['score'][sel], pos), pos.sum()))
    for k in range(0, 4):
        sel = ok & (m3['n_env'] == k)
        if sel.sum():
            P('   包絡に入る調和の数 %d: %4d 点、D > -5 の割合 %.3f' % (k, sel.sum(), (r.D[sel] > -5).mean()))
    for w in s3.warnings_for(m3):
        P('   警告: ' + w)
    P('')
    P('(3) 測ったチューンから距離 < tol にある共鳴線の最小次数の分布 (D > -5 の点 / 全点 (比))')
    for tol in (0.003, 0.006):
        o = FP.nearest_order_stats(r.D, r.nux1, r.nuy1, tol=tol)
        P('   tol=%.3f: ' % tol + '、'.join('%s: %.2f/%.2f (%.2f)' % (('%s次' % a if a != 'none' else '近くない'), h, b, c) for a, h, b, c in zip(o['orders'], o['frac_hi'], o['frac_all'], o['enrich'])))
    drv = set(_red(a) for a in s3.harms)
    hi = ok & (r.D > -5)
    keys = ['駆動項あり (次数 ≤ 4)', '駆動項なし (次数 ≤ 4、対称性で消える)', '次数 5', '次数 6', '近くない (次数 ≤ 6)']
    cnt = dict.fromkeys(keys, 0)
    allc = dict.fromkeys(keys, 0)
    lines = _mn_list(6)
    for i in np.nonzero(ok)[0]:
        best = None
        for (mm_, nn_) in lines:
            v = mm_ * r.nux1[i] + nn_ * r.nuy1[i]
            d = abs(v - round(v)) / math.hypot(mm_, nn_)
            if d < 0.003:
                o_ = abs(mm_) + abs(nn_)
                if best is None or o_ < best[0]:
                    best = (o_, (mm_, nn_))
        if best is None:
            k = keys[4]
        elif best[0] <= 4:
            k = keys[0] if _red(best[1]) in drv else keys[1]
        else:
            k = '次数 %d' % best[0]
        allc[k] += 1
        if hi[i]:
            cnt[k] += 1
    P('   最も低い次数の線 (tol = 0.003) の種類 (駆動項あり = その線の調和が 3・4 次の h に出る):')
    for k in keys:
        P('   %-36s D > -5: %3d (%.2f)、全点 %4d (%.2f)' % (k, cnt[k], cnt[k] / max(hi.sum(), 1), allc[k], allc[k] / max(ok.sum(), 1)))
    P('')
    P('(4) 調和ごとの寄与 (S3、解析 ν(J)、包絡、D > -5): 調和 / 次数 / 領域の点数 / cov (area) / 対照 / G0 / W (領域の中の中央値)')
    for a, v in sorted(m3['by_harmonic'].items()):
        if not v['envelope'].any():
            continue
        s = FR.coverage_stats(r.D, v['envelope'], r.grid, thr=-5.0, nshift=200)
        P('   %-8s %d次 %4d 点、cov %.3f (area %.3f)、対照 %.3f±%.3f、G0 = %.0f、W = %.3g' % (str(a), abs(a[0]) + abs(a[1]), v['envelope'].sum(), s['cov'], s['area'], s['shift_mean'], s['shift_std'], np.nanmedian(v['G0']), np.nanmedian(v['W'][v['envelope']])))
    P('')
    dn = np.hypot(nm[0] - m3['nu_model'][0], nm[1] - m3['nu_model'][1])
    P('(5) 解析の ν(J) と追跡のチューンの差 |Δν|: 全点の中央値 %.4f、90%% 点 %.4f。0.003 以内の点 %.2f。小振幅の 1/3 の中央値 %.4f' % (
        np.nanmedian(dn[ok]), np.nanpercentile(dn[ok], 90), np.mean(dn[ok] < 0.003), np.nanmedian(dn[ok & (amp < np.quantile(amp, 1 / 3))])))
    P('')
    e3 = res[('解析 ν(J) (予想)', 'S3  (+正規形の補正)')][1]['score']
    P('(5b) 振幅の 3 層 (小・中・大) ごとの S3 のスコアの AUC (層の中の判別力。振幅の効果を除く): %s、順位相関 ρ: %s' % (
        ', '.join('%.3f' % v for v in e3['auc_layers']), ', '.join('%.3f' % v for v in e3['rho_layers'])))
    P('')
    P('(6) 共鳴への捕捉 (チューンの固着) の検査 = D を使わない島の予想の検証: 調和ごとに、予想の「捕らえられる」点 (位相つき) と、その外側の輪 (W ≤ |Δ0| < 3W) で、')
    P('    測ったチューンの Δ0 (追跡) が解析の Δ0 より十分小さく線に寄っている (|Δ0meas| < 0.003 かつ < 0.25 |Δ0model|) 点の割合')
    for a, v in sorted(m3['by_harmonic'].items()):
        g_ = v['g']
        m0, n0 = a[0] // g_, a[1] // g_
        s_ = (m0 * nm[0] + n0 * nm[1]) * g_
        d0m = (s_ - np.round(s_)) / g_
        d0 = v['delta0']
        inn = v['trapped'] & ok
        ring = ok & ~v['trapped'] & (np.abs(d0) >= v['W']) & (np.abs(d0) < 3 * v['W'])
        lock = (np.abs(d0m) < 0.003) & (np.abs(d0m) < 0.25 * np.abs(d0))
        if inn.sum() < 5:
            continue
        P('   %-8s %d次: 捕らえられる %4d 点で固着 %.2f、外側の輪 %4d 点で %.2f' % (str(a), abs(a[0]) + abs(a[1]), inn.sum(), lock[inn].mean(), ring.sum(), (lock[ring].mean() if ring.sum() else float('nan'))))
    FR.plot_with_prediction(r, [('S1: 3rd order h (linear nu)', res[('解析 ν(J) (予想)', 'S1  (3 次)')][0]['envelope'], 'white', '-'),
                                ('S3: + 4th order, normal form (linear nu)', m3['envelope'], 'magenta', '--')],
                            os.path.join(HERE, 'S3_fmap_predict_%s.png' % name), title='%s: D and predicted resonance regions (S1 vs S3)' % name)
    with open(os.path.join(HERE, 'S3_fmap_predict_%s.out.txt' % name), 'w', encoding='utf-8') as f:
        f.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2')
