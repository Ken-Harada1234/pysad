# -*- coding: utf-8 -*-
"""局所運動量アクセプタンス (LMA) の例 (pysad/lma.py)。1 セル (1 周期) の位置を計算し、全周の別のセルの数点で 1 セルの結果と合うことを確かめ、
線形の見積もり・RF バケットの高さと比べ、Touschek 寿命 (既存の TouschekTable) に渡す。
使い方 (リポジトリ直下から): py -3.12 example/15_lma_touschek/lma/lma_ring.py astrid2   または   kekls
出力: example/15_lma_touschek/lma/lma_<名前>.table.txt (表)、.png (図)、標準出力 (.out.txt に保存した)。
条件 (どちらも FLUC なし、RAD ON、RFSW ON、v2 モード): アパーチャ ±APER (SAD の LOSSAMPL。台本の例の値で、実機の値ではない)、損失の判定は check 要素ごと (ASTRID2 は 1、kek-ls は 5)、
放射・RF つき、周回数 = シンクロトロン振動 2 周期分、バケット判定 ON、粗い走査 12 点 + 3 点ずつの二分 (幅 1e-4 まで)。
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)

CASES = {
    'astrid2': dict(path=os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), ring='RING', aper=0.010, ds=0.0,
                    label='ASTRID2 (0.58 GeV、周長 45.7 m、6 セル、RF 50 kV)', pbunch=1e9, check=1, extra_cells=(2, 3, 5)),
    'kekls': dict(path=os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), ring='RING', aper=0.010, ds=0.5,
                  label='kek-ls_v03_12.lat (3 GeV、周長 603 m、20 セル、RF 1.4 MV)', pbunch=1e9, check=5, extra_cells=(3, 9, 15)),
}


def main(name, workers=2):
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity, select_locations
    from pysad.lma import momentum_acceptance, linear_momentum_limit, touschek_lifetime_local
    compat.set_mode('v2')
    cs = CASES[name]
    src = open(cs['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, cs['ring'], rad=True, rfsw=True, alost=cs['aper'])
    st0 = st                              # SAD の tapert のままの判定 (check_every=0)
    st = st.copy(check_every=cs['check'])
    print(cs['label'])
    print('条件:', st.conditions())
    per = find_periodicity(st.els, st.rows)
    print('周期: セルの数 %d、セルの長さ %.4f m、セルの長さが全部同じ: %s、最初のセルの位置 %d 〜 %d' % (
        per['nper'], per['cell_length'], per['ok'], per['starts'][0], per['starts'][1] - 1))
    ls, le = per['starts'][0], per['starts'][1] - 1
    locs = select_locations(st, ls, le, 1, cs['ds'])
    print('1 セルの位置の数:', len(locs))
    t0 = time.time()
    r = momentum_acceptance(st, locs, workers=workers, progress=10)
    print('1 セル: %.0f s、追跡した粒子 %d' % (time.time() - t0, r.cond['tracked_particles']))
    r.save(os.path.join(HERE, 'lma_%s.table.txt' % name))
    # SAD の tapert のままの判定 (最初の要素・15 要素ごと・区間の終わり) との違い: 同じ位置の数点
    pick0 = [locs[len(locs) // 5], locs[len(locs) // 2], locs[(4 * len(locs)) // 5]]
    r0 = momentum_acceptance(st0, pick0, workers=workers)
    print()
    print('--- 損失の判定の間隔: check_every=%d (上) と 0 (SAD の tapert のまま) ---' % cs['check'])
    for i, l in enumerate(pick0):
        j = list(r.loc).index(l)
        print('loc %4d s=%8.3f  dp+ %.5f (check_every=%d) / %.5f (0)   dp- %.5f / %.5f' % (l, r.s[j], r.dp_plus[j], cs['check'], r0.dp_plus[i], r.dp_minus[j], r0.dp_minus[i]))
    # 全周の別のセルの数点
    L = per['cell_length']
    pick = [locs[len(locs) // 5], locs[len(locs) // 2], locs[(4 * len(locs)) // 5]]
    print('\n--- 周期性の確認: 別のセルの同じ位置 (s + k·セル長) ---')
    print('%4s %10s %4s %12s %12s %12s %12s' % ('cell', 's[m]', 'loc', 'dp+', 'dp-', 'd(dp+)', 'd(dp-)'))
    worst = 0.0
    for k in cs['extra_cells']:
        locs_k = []
        for l in pick:
            lk = per['starts'][k - 1] + (l - per['starts'][0])
            assert st.names[lk - 1] == st.names[l - 1] and abs(st.s_of(lk) - st.s_of(l) - (per['s_starts'][k - 1] - per['s_starts'][0])) < 1e-6
            locs_k.append(lk)
        rk = momentum_acceptance(st, locs_k, workers=workers)
        for i, l in enumerate(pick):
            j = list(r.loc).index(l)
            dpp, dpm = rk.dp_plus[i] - r.dp_plus[j], rk.dp_minus[i] - r.dp_minus[j]
            worst = max(worst, abs(dpp), abs(dpm))
            print('%4d %10.4f %4d %12.6f %12.6f %12.2e %12.2e' % (k, rk.s[i], locs_k[i], rk.dp_plus[i], rk.dp_minus[i], dpp, dpm))
    print('最大の差 |Δδ| = %.2e (二分の幅 tol = %.1e)' % (worst, r.cond['tol']))
    # 線形の見積もり、バケット
    lin = linear_momentum_limit(st, r.loc, cs['aper'])
    bh = st.bucket['height']
    print('\nRF バケットの高さ δ_rf = %.5f。δ+ の最大 %.5f (δ+/δ_rf の最大 %.3f)、δ− の最小 %.5f' % (bh, r.dp_plus.max(), r.dp_plus.max() / bh, r.dp_minus.min()))
    print('損失の理由 (δ+): rf %d 点、aper %d 点、sat %d 点' % tuple(sum(1 for x in r.reason_plus if x == k) for k in ('rf', 'aper', 'sat')))
    print('線形の見積もり δ_lin (水平の開口 ±%.3f m) と δ+: δ+ ≤ min(δ_rf, δ_lin) × (1+ε) の ε の最大 = %.3f' % (
        cs['aper'], np.max(r.dp_plus / np.minimum(bh, lin)) - 1))
    # 図
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pysad.lma import _jfont
    _jfont(plt)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(r.s, r.dp_plus * 100, 'o-', ms=3, lw=0.8, label='dp+ (tracking)')
    ax.plot(r.s, r.dp_minus * 100, 'o-', ms=3, lw=0.8, label='dp- (tracking)')
    ax.plot(r.s, np.minimum(lin, bh) * 100, 'k--', lw=0.8, label='min(RF bucket, linear estimate)')
    ax.plot(r.s, -np.minimum(lin, bh) * 100, 'k--', lw=0.8)
    ax.set_xlabel('s [m] (one cell)')
    ax.set_ylabel('local momentum acceptance [%]')
    ax.set_title('%s, aperture +-%.0f mm, nturn=%d' % (name, cs['aper'] * 1e3, r.cond['nturn']))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'lma_%s.png' % name), dpi=130)
    plt.close(fig)
    # Touschek
    ffs.flags['INTRA'] = True
    ffs.I.globals['PBUNCH'] = cs['pbunch']
    ffs.flags['RFSW'] = True
    ffs.flags['RAD'] = True
    ffs.run_emit(quiet=True)
    tab = ffs.touschek['table']
    res = touschek_lifetime_local(tab, r, setup=st, period=per['cell_length'], s0=per['s_starts'][0])
    from pysad.touschek import touschek_lifetime
    sige = float(ffs.I.globals['SIGE'])
    t_rf = touschek_lifetime(tab, None, None, bh / sige, sige)
    t_min = touschek_lifetime(tab, None, None, min(r.dp_plus.min(), -r.dp_minus.max()) / sige, sige)
    print('\nTouschek (PBUNCH=%.1e 粒子/バンチ、INTRA ON、SIGE=%.4e、表の dp は 0.004〜0.242):' % (cs['pbunch'], sige))
    print('  局所の δ±(s) (LMA、周期的に延ばす): 寿命 %.4f s (表の外を使った要素 %d)' % (res['lifetime'], res['extrapolated']))
    print('  一様な受容 = RF バケットの高さ %.4f (TouschekLifetime[INF,INF,nz]): %.4f s' % (bh, t_rf))
    print('  一様な受容 = LMA の最小値 %.4f: %.4f s' % (min(r.dp_plus.min(), -r.dp_minus.max()), t_min))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', workers=int(sys.argv[2]) if len(sys.argv) > 2 else 2)
