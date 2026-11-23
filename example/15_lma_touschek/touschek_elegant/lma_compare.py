# -*- coding: utf-8 -*-
"""pysad の LMA (lma.py) と elegant の momentum_aperture の局所運動量アクセプタンスを同じ位置で比べる。

使い方 (リポジトリ直下から。PYTHONUTF8=1 を付ける):
  PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_compare.py astrid2 [--pysad]     ASTRID2: pysad の表は lma_astrid2.table.txt (check_every = 1)。--pysad で pysad_astrid2.lma.txt を作り直す
  PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_compare.py kekls --pysad [--workers 3]   kek-ls: pysad を elegant と同じ位置 (1 セル + 別のセルの 3 点) で走らせてから比べる
先に elegant_lma.py で elegant_<名前>.lma.txt を作る。
出力: lma_compare_<名前>.table.txt (位置ごとの比較)、lma_compare_<名前>.png (図)、標準出力 (.out.txt に保存した)。
pysad の条件: アパーチャ ±10 mm、check_every = 1 (要素ごとの判定)、周回数 = シンクロトロン振動 2 周期、RF バケットの判定あり、粗い走査 12 点 + 3 点ずつの二分 (tol 1e-4)。
elegant の条件は elegant_lma.py の冒頭。
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

RF_FRAC = 0.97            # |δ| がバケットの高さのこの割合以上なら「RF が限界」(それより小さければアパーチャが限界)。elegant の損失の周 (lostOnPass) は
                          # 限界ぎりぎりのアパーチャ損失でもシンクロトロン振動の半周期ほど遅れるので、判定に使わない


def read_pysad(path):
    rows = []
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        rows.append(dict(loc=int(t[0]), s=float(t[1]), name=t[2], dpp=float(t[3]), dpm=float(t[4]), hip=float(t[5]), him=float(t[6]),
                         rp=t[7], rm=t[8]))
    return rows


def read_elegant(path):
    rows = {}
    for ln in open(path, encoding='utf-8'):
        if ln.startswith('#') or not ln.strip():
            continue
        t = ln.split()
        rows[int(t[0])] = dict(loc=int(t[0]), s=float(t[1]), dpp=float(t[2]), dpm=float(t[3]), pp=int(t[4]), pm=int(t[5]),
                               slp=float(t[6]), slm=float(t[7]))
    return rows


def run_pysad(name, locs, workers, aper=0.010, check=1):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    from pysad.lma import momentum_acceptance
    import sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    st = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True, alost=aper)[0].copy(check_every=check)
    t0 = time.time()
    r = momentum_acceptance(st, locs, workers=workers, progress=10)
    print('pysad LMA: %d 位置、%.0f s、追跡した粒子 %d' % (len(locs), time.time() - t0, r.cond['tracked_particles']))
    r.save(os.path.join(HERE, 'pysad_%s.lma.txt' % name))
    return st


def main(name, do_pysad=False, workers=3):
    import sad2ele as c2e
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity
    compat.set_mode('v2')
    st0 = TrackSetup.from_source(open(c2e.CASES[name]['path'], encoding='utf-8').read(), 'RING', rad=True, rfsw=True, alost=0.010)[0]
    per = find_periodicity(st0.els, st0.rows)
    import elegant_lma as el
    if do_pysad and not os.path.exists(os.path.join(HERE, 'elegant_%s.lma.txt' % name)):
        run_pysad(name, el.lma_locs(name, st0, per), workers)
        return
    ele = read_elegant(os.path.join(HERE, 'elegant_%s.lma.txt' % name))
    ls, le = per['starts'][0], per['starts'][1] - 1
    cell = [l for l in sorted(ele) if ls <= l <= le]
    if name == 'astrid2':
        ppath = os.path.join(HERE, '..', 'lma', 'lma_astrid2.table.txt')
        if do_pysad:
            run_pysad(name, cell, workers)
            ppath = os.path.join(HERE, 'pysad_astrid2.lma.txt')
    else:
        if do_pysad:
            run_pysad(name, sorted(ele), workers)
        ppath = os.path.join(HERE, 'pysad_%s.lma.txt' % name)
    py = read_pysad(ppath)
    pyd = {r['loc']: r for r in py}
    bh = st0.bucket['height']
    nturn = st0.nturn_for_synchrotron(2.0)
    out = []
    P = print
    P('%s: LMA の比較 (pysad = %s、elegant = elegant_%s.lma.txt)' % (name, os.path.basename(ppath), name))
    P('RF バケットの高さ δ_rf = %.5f、周回数 %d、アパーチャ ±10 mm' % (bh, nturn))
    lines = ['# loc  s[m]  name  pysad_dp+  elegant_dp+  rel_diff+[%]  pysad_dp-  elegant_dp-  rel_diff-[%]  limit+(pysad/elegant)  limit-(pysad/elegant)  elegant_lostpass+ -']
    cmp = []
    for l in cell:
        if l not in pyd:
            continue
        a, e = pyd[l], ele[l]
        kp_e = 'R' if e['dpp'] >= RF_FRAC * bh else 'A'
        km_e = 'R' if -e['dpm'] >= RF_FRAC * bh else 'A'
        a['kp'] = 'R' if a['dpp'] >= RF_FRAC * bh else 'A'
        a['km'] = 'R' if -a['dpm'] >= RF_FRAC * bh else 'A'
        rp = 100 * (e['dpp'] / a['dpp'] - 1)
        rm = 100 * (e['dpm'] / a['dpm'] - 1)
        cmp.append((l, a, e, rp, rm, kp_e, km_e))
        lines.append('%5d %9.4f %-9s %10.6f %10.6f %7.2f  %10.6f %10.6f %7.2f  %s/%s  %s/%s  %d %d' % (
            l, a['s'], a['name'], a['dpp'], e['dpp'], rp, a['dpm'], e['dpm'], rm, a['rp'], kp_e, a['rm'], km_e, e['pp'], e['pm']))
    open(os.path.join(HERE, 'lma_compare_%s.table.txt' % name), 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')
    P('\n位置の数 (1 セル): %d' % len(cmp))
    d = np.array([[c[2]['dpp'] - c[1]['dpp'], c[2]['dpm'] - c[1]['dpm']] for c in cmp])
    rel = np.array([[c[3], c[4]] for c in cmp])
    P('最大の差 |Δδ+| = %.2e、|Δδ−| = %.2e (elegant の最終刻み %.2e、pysad の tol 1e-4)' % (np.abs(d[:, 0]).max(), np.abs(d[:, 1]).max(),
                                                                          {'astrid2': 0.0025 / 32, 'kekls': 0.005 / 64}[name]))
    P('相対差 (elegant/pysad − 1) の平均 [%%]: δ+ %.2f、δ− %.2f、最大 |.|: %.2f、%.2f' % (rel[:, 0].mean(), rel[:, 1].mean(), np.abs(rel[:, 0]).max(), np.abs(rel[:, 1]).max()))
    # 限界の種類ごと
    def grp(sign):
        g = {}
        for c in cmp:
            a, e = c[1], c[2]
            kp = (a['kp'] if sign > 0 else a['km'])
            ke = (c[5] if sign > 0 else c[6])
            g.setdefault((kp, ke), []).append(c[3] if sign > 0 else c[4])
        return g
    for sign, lab in ((1, 'δ+'), (-1, 'δ−')):
        for (kp, ke), v in sorted(grp(sign).items()):
            P('  %s: 限界の種類: pysad %s / elegant %s (A: アパーチャ、R: RF バケット。|δ| ≥ %.2f δ_rf): %2d 点、相対差の平均 %+.2f %%、範囲 %+.2f 〜 %+.2f %%' % (
                lab, kp, ke, RF_FRAC, len(v), np.mean(v), min(v), max(v)))
    dmax = np.array([[c[1]['dpp'], c[2]['dpp'], c[1]['dpm'], c[2]['dpm']] for c in cmp])
    P('δ+ の最小 (pysad / elegant): %.5f / %.5f、最大 %.5f / %.5f。δ− の最小の絶対値 %.5f / %.5f、最大 %.5f / %.5f' % (
        dmax[:, 0].min(), dmax[:, 1].min(), dmax[:, 0].max(), dmax[:, 1].max(),
        np.abs(dmax[:, 2]).min(), np.abs(dmax[:, 3]).min(), np.abs(dmax[:, 2]).max(), np.abs(dmax[:, 3]).max()))
    # elegant の周期性: 別のセルの同じ位置との差
    off = []
    for k in range(2, per['nper'] + 1):
        for l in cell:
            lk = per['starts'][k - 1] + (l - per['starts'][0])
            if lk in ele:
                off.append((k, l, ele[lk]['dpp'] - ele[l]['dpp'], ele[lk]['dpm'] - ele[l]['dpm']))
    if off:
        o = np.array([[x[2], x[3]] for x in off])
        P('\nelegant の周期性 (別のセルの同じ位置 %d 点と最初のセルの差): 最大 |Δδ+| = %.2e、|Δδ−| = %.2e、二乗平均 %.2e・%.2e' % (
            len(off), np.abs(o[:, 0]).max(), np.abs(o[:, 1]).max(), np.sqrt((o[:, 0] ** 2).mean()), np.sqrt((o[:, 1] ** 2).mean())))
    # 図
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    s = [c[1]['s'] for c in cmp]
    ax[0].plot(s, [c[1]['dpp'] * 100 for c in cmp], 'o-', ms=3, lw=0.8, label='pysad dp+')
    ax[0].plot(s, [c[2]['dpp'] * 100 for c in cmp], 's--', ms=3, lw=0.8, label='elegant dp+')
    ax[0].plot(s, [-c[1]['dpm'] * 100 for c in cmp], 'o-', ms=3, lw=0.8, label='pysad |dp-|')
    ax[0].plot(s, [-c[2]['dpm'] * 100 for c in cmp], 's--', ms=3, lw=0.8, label='elegant |dp-|')
    ax[0].axhline(bh * 100, color='gray', ls=':', lw=0.8, label='RF bucket height')
    ax[0].set_ylabel('local momentum acceptance [%]')
    ax[0].legend(fontsize=7, ncol=3)
    ax[0].grid(alpha=0.3)
    ax[1].plot(s, [c[3] for c in cmp], 'o-', ms=3, lw=0.8, label='dp+')
    ax[1].plot(s, [c[4] for c in cmp], 's-', ms=3, lw=0.8, label='dp-')
    ax[1].set_ylabel('elegant / pysad - 1 [%]')
    ax[1].set_xlabel('s [m] (first cell)')
    ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)
    ax[0].set_title('%s: local momentum acceptance, aperture +-10 mm, nturn=%d' % (name, nturn))
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'lma_compare_%s.png' % name), dpi=120)
    plt.close(fig)


if __name__ == '__main__':
    w = int(sys.argv[sys.argv.index('--workers') + 1]) if '--workers' in sys.argv else 3
    main(sys.argv[1], do_pysad='--pysad' in sys.argv, workers=w)
