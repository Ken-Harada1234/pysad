# -*- coding: utf-8 -*-
"""v2.3.0 L2 項目 1: check_every (損失の判定の間隔) の影響と計算時間。elegant (要素ごとの MAXAMP = check_every 1 と同じ) と比べる。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/L2_check_every.py astrid2|kekls [workers]
check_every = 0 (SAD の tapert: 最初の要素・15 要素ごと・区間の終わり)、1 (要素ごと、既定)、5 の δ± と時間を、最初のセルの位置で比べる。"""
import sys
import time
import numpy as np
from L2_common import *


def main(name, workers=2):
    from pysad.lma import momentum_acceptance
    from pysad.trksetup import find_periodicity, select_locations
    st, ffs = load(name)
    cs = CASES[name]
    per = find_periodicity(st.els, st.rows)
    locs = select_locations(st, per['starts'][0], per['starts'][1] - 1, 1, cs['ds'])
    if name == 'kekls':
        ele = elegant_kekls_locs()
        locs = [l for l in locs if l in ele]          # T1 の elegant と同じ位置 (最初のセルの 11 点)
        eref = ele
    else:
        eref = elegant_astrid2_table()
    print(cs['label'], '位置', len(locs), 'LOSSAMPL ±%.0f mm' % (cs['aper'] * 1e3))
    R = {}
    for ce in (0, 1, 5):
        t0 = time.time()
        R[ce] = momentum_acceptance(st, locs, check_every=ce, workers=workers, quiet=True, turncheck=0)
        print('check_every = %d: %.0f s、粒子 %d、nturn %d' % (ce, time.time() - t0, R[ce].cond['tracked_particles'], R[ce].cond['nturn']))
    print('\n loc   s[m]   | δ+: ce=0   ce=1   ce=5  elegant | δ−: ce=0   ce=1  ce=5  elegant | reason(ce=0,1)+')
    for i, l in enumerate(locs):
        e = eref.get(l)
        print('%5d %8.3f | %.5f %.5f %.5f %s | %.5f %.5f %.5f %s | %s %s' % (
            l, R[1].s[i], R[0].dp_plus[i], R[1].dp_plus[i], R[5].dp_plus[i], ('%.5f' % e[1]) if e else '   -   ',
            R[0].dp_minus[i], R[1].dp_minus[i], R[5].dp_minus[i], ('%.5f' % e[2]) if e else '   -   ', R[0].reason_plus[i], R[1].reason_plus[i]))
    a, b = R[0], R[1]
    print('\n(ce=0) − (ce=1): δ+ の最大 %.4f (最小 %.4f)、δ− の最大 %.4f。比 (ce0/ce1) の最大: δ+ %.3f、δ− %.3f' % (
        np.max(a.dp_plus - b.dp_plus), np.min(a.dp_plus - b.dp_plus), np.max(-(a.dp_minus - b.dp_minus)),
        np.max(a.dp_plus / b.dp_plus), np.max(a.dp_minus / b.dp_minus)))
    print('(ce=5) − (ce=1): δ+ の最大 %.2e、δ− の最大 %.2e' % (np.max(np.abs(R[5].dp_plus - b.dp_plus)), np.max(np.abs(R[5].dp_minus - b.dp_minus))))
    if eref:
        m = [i for i, l in enumerate(locs) if l in eref]
        for ce in (0, 1):
            dp = [R[ce].dp_plus[i] - eref[locs[i]][1] for i in m]
            dm = [R[ce].dp_minus[i] - eref[locs[i]][2] for i in m]
            print('pysad(ce=%d) − elegant: δ+ の最大 %.2e (範囲 %.2e〜%.2e)、δ− の最大 %.2e' % (ce, np.max(np.abs(dp)), min(dp), max(dp), np.max(np.abs(dm))))
    for ce in (0, 1, 5):
        R[ce].save(os.path.join(HERE, 'L2_check_every_%s_ce%d.lma.txt' % (name, ce)))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 2)
