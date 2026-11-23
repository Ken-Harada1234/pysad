# -*- coding: utf-8 -*-
"""1 セルの結果を全周に延ばす影響。RF 空洞は 1 か所なので、全セルの同じ位置の LMA は同じにならない (周期性が崩れる)。
最初のセルの 3 位置を、全部のセル (ASTRID2 6、kek-ls 20) の同じ位置で計算し、空洞からの距離との関係・最初のセルとの差を U0/E0 (閉軌道の δ のセルごとの変化) と比べる。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_cell.py astrid2|kekls [workers] [tol]  (tol を小さくすると U0/E0 の大きさの差が見える)"""
import sys
import time
import numpy as np
from lma_common import *


def main(name, workers=2, tol=1e-4):
    from pysad.lma import momentum_acceptance, periodicity_report
    from pysad.trksetup import select_locations
    st, ffs = load(name)
    cs = CASES[name]
    rep = periodicity_report(st)
    print(cs['label'], 'セル %d、セルの長さ %.4f m、閉軌道のセル間の差 %.2e、U0/E0 = %.3e、nturn %d' % (
        rep['nper'], rep['cell_length'], rep['cod_spread'], rep['delta_sawtooth'], st.nturn_for_synchrotron(2.0)))
    cav = [i for i, e in enumerate(st.els) if type(e).__name__ == 'Cavity'][0] + 1
    s_cav = st.s_of(cav)
    print('空洞: 要素 %d、s = %.4f m' % (cav, s_cav))
    first = select_locations(st, rep['starts'][0], rep['starts'][1] - 1, 1, cs['ds'])
    if name == 'astrid2':
        picks = [2, 17, first[len(first) // 2]]
    else:
        ele = elegant_kekls_locs()
        cand = [l for l in first if l in ele]
        picks = [cand[1], cand[len(cand) // 2], cand[-2]]
    picks = list(dict.fromkeys(picks))
    locs, tag = [], []
    for k in range(rep['nper']):
        for l in picks:
            lk = rep['starts'][k] + (l - rep['starts'][0])
            assert st.names[lk - 1] == st.names[l - 1]
            locs.append(lk)
            tag.append((k + 1, l))
    t0 = time.time()
    r = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0, tol=tol)
    print('位置 %d、%.0f s、粒子 %d' % (len(locs), time.time() - t0, r.cond['tracked_particles']))
    L = st.circ
    res = {}
    for i, (k, l) in enumerate(tag):
        res[(k, l)] = (r.s[i], r.dp_plus[i], r.dp_minus[i], r.reason_plus[i])
    for l in picks:
        print('\n位置 loc %d (最初のセル、%s): セル k、s、空洞の下流の距離 d = (s − s_cav) mod C、δ+、δ−、最初のセルとの差' % (l, st.names[l - 1]))
        b = res[(1, l)]
        dd_p, dd_m, dist = [], [], []
        for k in range(1, rep['nper'] + 1):
            s, dp, dm, rs = res[(k, l)]
            d = (s - s_cav) % L
            dist.append(d); dd_p.append(dp - b[1]); dd_m.append(dm - b[2])
            print('  k=%2d s=%8.3f d=%8.3f  δ+ %.5f δ− %.5f  Δδ+ %+.1e Δδ− %+.1e  %s' % (k, s, d, dp, dm, dp - b[1], dm - b[2], rs))
        ap = np.array([res[(k, l)][1] for k in range(1, rep['nper'] + 1)])
        am = np.array([res[(k, l)][2] for k in range(1, rep['nper'] + 1)])
        print('  全セルの δ+ の幅 (最大−最小) %.2e、δ− %.2e。最初のセルとの差の最大: δ+ %.2e、δ− %.2e (tol %.0e、U0/E0 %.2e)' % (
            ap.max() - ap.min(), am.max() - am.min(), np.max(np.abs(dd_p)), np.max(np.abs(dd_m)), r.cond['tol'], rep['delta_sawtooth']))
    r.save(os.path.join(HERE, 'lma_cell_%s%s.lma.txt' % (name, '' if tol == 1e-4 else '_tol%g' % tol)))
    allp = np.array([max(abs(res[(k, l)][1] - res[(1, l)][1]), abs(res[(k, l)][2] - res[(1, l)][2])) for (k, l) in res])
    print('\n全体の最大の差 %.2e (tol %.0e、U0/E0 %.2e)' % (allp.max(), r.cond['tol'], rep['delta_sawtooth']))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 2, float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4)
