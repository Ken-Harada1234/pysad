# -*- coding: utf-8 -*-
"""シンクロトロン振動 1 周期未満・数周期での収束。周回数をシンクロトロン周期の 0.25〜3 倍にして、バケットで決まる δ+ と、
1 次元の写像の限界 (バケットのみ)、LMA (aper 判定) を比べる。ASTRID2 (位置 2: δ+ はアパーチャ、δ− はバケット付近)、kek-ls (位置 3: バケット)。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/lma/lma_turns_short.py [workers]"""
import sys
import numpy as np
from lma_common import *


def main(workers=2):
    from pysad.lma import momentum_acceptance
    from pysad import lmabucket as lb
    for name, loc in (('astrid2', 2), ('kekls', 3)):
        st, ffs = load(name)
        nus = abs(st.bucket['nus'])
        per = 1.0 / nus
        par = lb.longitudinal_params(st)
        print('\n%s: δ_rf = %.5f、νs = %.5f (周期 %.0f 周)、位置 loc %d' % (name, st.bucket['height'], nus, per, loc))
        print('  周期数  nturn | LMA δ+ δ− (理由, 旗) | 6D バケット δ+ δ− | 1D 写像 δ+ δ−')
        for f in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0):
            N = int(round(f * per))
            r = momentum_acceptance(st, [loc], nturn=N, workers=workers, quiet=True, check_every=5 if name == 'kekls' else 1, turncheck=2.0, dpmax=2.0 * st.bucket['height'])
            b = lb.bucket_limits_6d(st, [loc], nturn=N, dpmax=2.0 * st.bucket['height'], check_every=5 if name == 'kekls' else 0)
            m = lb.bucket_limits_1d(st, N, par=par, dp_top=2.0 * st.bucket['height'])
            f1 = lambda v: '   none ' if v is None else '%.5f' % v
            print('  %4.2f  %5d | %.5f %.5f (%s %s, %s) | %.5f %.5f | %s %s' % (
                f, N, r.dp_plus[0], r.dp_minus[0], r.reason_plus[0], r.reason_minus[0], r.flags(0), b.dp_plus[0], b.dp_minus[0], f1(m['dp_plus']), f1(m['dp_minus'])), flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 2)
