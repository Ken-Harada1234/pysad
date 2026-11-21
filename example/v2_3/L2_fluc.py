# -*- coding: utf-8 -*-
"""v2.3.0 L2 項目 6: FLUC (量子励起) つきの LMA。複数の乱数の種 (既定 6) で pysad の δ± の平均・標準偏差を求め、FLUC なしと、elegant の ISR つき
(L2_elegant_isr.py の出力 L2_elegant_isr_<名前>.lma.txt があれば) と比べる。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/L2_fluc.py astrid2|kekls [同時に走らせる種の数] [nseed] [tol] [nturn (0 = 既定)] [位置 例 2,17]"""
import json
import sys
import time
import numpy as np
from L2_common import *

PICK = {'astrid2': [2, 17, 40, 100], 'kekls': None}


def _one(a):
    """1 つの種 (並列の子プロセス): (seed, dp_plus, dp_minus, reasons, 秒)"""
    name, seed, locs, tol, nturn = a
    from pysad.lma import momentum_acceptance
    st, ffs = load(name)
    sf = st.copy(fluc=True)
    sf.seed = seed
    t0 = time.time()
    r = momentum_acceptance(sf, locs, workers=1, quiet=True, turncheck=0, tol=tol, nturn=nturn)
    return seed, r.dp_plus.copy(), r.dp_minus.copy(), list(r.reason_plus), time.time() - t0


def main(name, workers=2, nseed=6, tol=1e-4, nturn=None, pick=None):
    from pysad.lma import momentum_acceptance
    st, ffs = load(name)
    cs = CASES[name]
    if pick:
        locs = pick
    elif name == 'astrid2':
        locs = PICK[name]
    else:
        ele = elegant_kekls_locs()
        c = [l for l in sorted(ele) if l < st.nlat // 10]
        locs = [c[1], c[len(c) // 2]]
    print(cs['label'], '位置', locs, '種', nseed, 'σδ = %s' % ffs.emit.get('momentum_spread') if hasattr(ffs, 'emit') and ffs.emit else '')
    t0 = time.time()
    r0 = momentum_acceptance(st, locs, workers=workers, quiet=True, turncheck=0, tol=tol, nturn=nturn)
    print('FLUC なし: %.0f s' % (time.time() - t0))
    S = []
    # 種ごとに別のプロセス (workers 個まで同時、位置は子の中で順に)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max(1, workers)) as ex:
        for seed, dpp, dpm, rs, sec in ex.map(_one, [(name, sd, locs, tol, nturn) for sd in range(1, nseed + 1)]):
            S.append((dpp, dpm, rs))
            print('seed %d: %.0f s  δ+ %s  δ− %s' % (seed, sec, ' '.join('%.5f' % v for v in dpp), ' '.join('%.5f' % v for v in dpm)), flush=True)
    P = np.array([s[0] for s in S])
    M = np.array([s[1] for s in S])
    eleg = {}
    p = os.path.join(HERE, 'L2_elegant_isr_%s.lma.txt' % name)
    if os.path.exists(p):
        for ln in open(p, encoding='utf-8'):
            if ln.startswith('#') or not ln.strip():
                continue
            t = ln.split()
            eleg.setdefault(int(t[0]), []).append((float(t[3]), float(t[4])))
    print('\n loc  s | FLUC なし δ+ δ− | pysad FLUC 平均 δ+ (σ) δ− (σ) | elegant ISR 平均 δ+ (σ) δ− (σ) (種の数) | 差の平均 pysad(FLUC) − なし')
    out = {}
    for i, l in enumerate(locs):
        e = np.array(eleg.get(l, []))
        es = ('%.5f (%.1e) %.5f (%.1e) (%d)' % (e[:, 0].mean(), e[:, 0].std(ddof=1), e[:, 1].mean(), e[:, 1].std(ddof=1), len(e))) if len(e) > 1 else '-'
        print('%5d %7.2f | %.5f %.5f | %.5f (%.1e) %.5f (%.1e) | %s | %+.2e %+.2e' % (
            l, r0.s[i], r0.dp_plus[i], r0.dp_minus[i], P[:, i].mean(), P[:, i].std(ddof=1), M[:, i].mean(), M[:, i].std(ddof=1), es,
            P[:, i].mean() - r0.dp_plus[i], M[:, i].mean() - r0.dp_minus[i]))
        out[l] = dict(nofluc=[r0.dp_plus[i], r0.dp_minus[i]], fluc_p=P[:, i].tolist(), fluc_m=M[:, i].tolist())
    json.dump(out, open(os.path.join(HERE, 'L2_fluc_%s%s.json' % (name, '' if (tol == 1e-4 and nturn is None) else '_tol%g_n%s' % (tol, nturn))), 'w'))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 2, int(sys.argv[3]) if len(sys.argv) > 3 else 6,
         float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4, int(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5] != '0' else None,
         [int(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else None)
