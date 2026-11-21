# -*- coding: utf-8 -*-
"""v2.3.0 L2 項目 3: 二分法が最初の損失しか見ないことの影響。`verify=n` (限界の手前 [0, δ] を n 点で再追跡) で、粗い走査の点のあいだの狭い損失の窓を数える。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/L2_verify.py astrid2|kekls [n] [workers]
n の既定: astrid2 100 (刻み 1.6e-4)、kekls 150 (刻み 2e-4)。窓があれば限界が下がり (nwin > 0、旗 W)、なければ限界は verify なしと同じ。"""
import sys
import time
import numpy as np
from L2_common import *


def main(name, n, workers=2):
    from pysad.lma import momentum_acceptance
    from pysad.trksetup import find_periodicity, select_locations
    st, ffs = load(name)
    cs = CASES[name]
    per = find_periodicity(st.els, st.rows)
    locs = select_locations(st, per['starts'][0], per['starts'][1] - 1, 1, cs['ds'])
    if name == 'kekls':
        locs = [l for l in locs if l in elegant_kekls_locs()]
    print(cs['label'], '位置', len(locs), 'verify = %d' % n)
    t0 = time.time()
    r0 = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0, verify=0)
    t_a = time.time() - t0
    t0 = time.time()
    r1 = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0, verify=n)
    t_b = time.time() - t0
    print('verify なし: %.0f s、粒子 %d。verify = %d: %.0f s、粒子 %d' % (t_a, r0.cond['tracked_particles'], n, t_b, r1.cond['tracked_particles']))
    nw = int(r1.nwin_plus.sum() + r1.nwin_minus.sum())
    print('窓の点 (verify で失われた手前の点) の合計 %d。窓のある限界: δ+ %d、δ− %d (全 %d 位置 × 2)' % (nw, (r1.nwin_plus > 0).sum(), (r1.nwin_minus > 0).sum(), len(locs)))
    d = np.maximum(np.abs(r1.dp_plus - r0.dp_plus), np.abs(r1.dp_minus - r0.dp_minus))
    print('verify あり − なし: 最大の差 %.2e (tol %.0e)' % (d.max(), r0.cond['tol']))
    for i in np.nonzero((r1.nwin_plus > 0) | (r1.nwin_minus > 0))[0]:
        print('  窓: loc %d s=%.3f δ+ %.5f→%.5f (最初の窓 %.5f、点 %d)、δ− %.5f→%.5f (最初の窓 %.5f、点 %d)' % (
            r1.loc[i], r1.s[i], r0.dp_plus[i], r1.dp_plus[i], r1.win_plus[i], r1.nwin_plus[i], r0.dp_minus[i], r1.dp_minus[i], r1.win_minus[i], r1.nwin_minus[i]))
    r1.save(os.path.join(HERE, 'L2_verify_%s.lma.txt' % name))


if __name__ == '__main__':
    nm = sys.argv[1] if len(sys.argv) > 1 else 'astrid2'
    main(nm, int(sys.argv[2]) if len(sys.argv) > 2 else (100 if nm == 'astrid2' else 150), int(sys.argv[3]) if len(sys.argv) > 3 else 2)
