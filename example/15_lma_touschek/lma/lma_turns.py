# -*- coding: utf-8 -*-
"""LMA の周回数への依存 (kek-ls)。elegant_lma.py の elegant (1500 周、14 位置) と同じ位置で、周回数 491・1000・1500・2000 の δ± を求める。
check_every = 5 (kek-ls は全位置がバケットの限界で、check_every 0・1・5 の結果が同じ: lma_check_every.py。1 が既定だが 2.2 倍遅い)。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_turns.py [workers]
  (a) bucket 判定あり (既定): |z| がバケットを出たら損失。周回数で δ± がどう変わるか。1 次元の写像・6 次元 (アパーチャ無限大) と同じ周回数で比べる。
  (b) bucket 判定なし (bucket=False、dpmax = 5%): elegant の momentum_aperture と同じ判定 (アパーチャに当たるまで)。同じ 1500 周で elegant と比べる。"""
import sys
import time
import numpy as np
from lma_common import *


def main(workers=2, turns=(491, 1000, 1500, 2000)):
    from pysad.lma import momentum_acceptance
    from pysad import lmabucket as lb
    st, ffs = load('kekls')
    ele = elegant_kekls_locs()
    first = [l for l in sorted(ele) if l < st.nlat // 10]       # 最初のセルの位置 (9 番目のセルの 3 点を除く)
    locs = first
    bh, nus = st.bucket['height'], abs(st.bucket['nus'])
    print('kek-ls: δ_rf = %.5f、νs = %.5f (周期 %.0f 周)、elegant 1500 周、位置 %d (最初のセル)' % (bh, nus, 1 / nus, len(locs)))
    par = lb.longitudinal_params(st)
    R = {}
    for N in turns:
        t0 = time.time()
        R[N] = momentum_acceptance(st, locs, nturn=N, workers=workers, quiet=True, turncheck=0, check_every=5)
        m1 = lb.bucket_limits_1d(st, N, par=par)
        print('(a) nturn = %d (%.2f 周期): %.0f s。1 次元の写像 δ+ = %.5f、|δ−| = %.5f' % (N, N * nus, time.time() - t0, m1['dp_plus'], m1['dp_minus']))
        R[N].save(os.path.join(HERE, 'lma_turns_kekls_n%d.lma.txt' % N))
    print('\n loc  s[m] | δ+ (491, 1000, 1500, 2000) | δ− | reason+(491, 1500)')
    for i, l in enumerate(locs):
        print('%5d %7.2f | %s | %s | %s %s' % (l, R[turns[0]].s[i], ' '.join('%.5f' % R[N].dp_plus[i] for N in turns), ' '.join('%.5f' % R[N].dp_minus[i] for N in turns),
                                              R[491].reason_plus[i], R[1500].reason_plus[i]))
    for N in turns:
        print('nturn %d: δ+ 平均 %.5f (最小 %.5f 最大 %.5f)、|δ−| 平均 %.5f (最小 %.5f 最大 %.5f)、理由 δ+ %s' % (
            N, R[N].dp_plus.mean(), R[N].dp_plus.min(), R[N].dp_plus.max(), -R[N].dp_minus.mean(), -R[N].dp_minus.max(), -R[N].dp_minus.min(),
            {k: R[N].reason_plus.count(k) for k in set(R[N].reason_plus)}))
    # 6 次元 (アパーチャ無限大) のバケットの限界
    print('\n6 次元 (アパーチャ無限大)・1 位置 (loc %d) のバケットの限界の周回数依存:' % locs[0])
    for N in turns:
        rf = lb.bucket_limits_6d(st, [locs[0]], nturn=N)
        print('  nturn %d: δ+ = %.5f、δ− = %.5f' % (N, rf.dp_plus[0], rf.dp_minus[0]))
    # (b) elegant と同じ判定
    N = 1500
    t0 = time.time()
    rb = momentum_acceptance(st, locs, nturn=N, workers=workers, quiet=True, bucket=False, dpmax=0.05, turncheck=0, check_every=5)
    print('\n(b) bucket 判定なし、1500 周、dpmax = 5%%: %.0f s' % (time.time() - t0))
    rb.save(os.path.join(HERE, 'lma_turns_kekls_n1500_nobucket.lma.txt'))
    dp, dm, dpb, dmb = [], [], [], []
    for i, l in enumerate(locs):
        e = ele[l]
        dp.append(R[1500].dp_plus[i] - e[1]); dm.append(R[1500].dp_minus[i] - e[2])
        dpb.append(rb.dp_plus[i] - e[1]); dmb.append(rb.dp_minus[i] - e[2])
        print('%5d  elegant %.5f %.5f | bucket 判定あり %.5f %.5f | なし %.5f %.5f (理由 %s %s)' % (
            l, e[1], e[2], R[1500].dp_plus[i], R[1500].dp_minus[i], rb.dp_plus[i], rb.dp_minus[i], rb.reason_plus[i], rb.reason_minus[i]))
    print('1500 周: pysad(bucket 判定あり) − elegant: δ+ %.2e〜%.2e、δ− %.2e〜%.2e' % (min(dp), max(dp), min(dm), max(dm)))
    print('1500 周: pysad(bucket 判定なし) − elegant: δ+ %.2e〜%.2e、δ− %.2e〜%.2e' % (min(dpb), max(dpb), min(dmb), max(dmb)))


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
