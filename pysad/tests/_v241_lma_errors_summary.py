# -*- coding: utf-8 -*-
"""v2.4.1: 誤差つき ASTRID2 の LMA (gs・pysad の二分法の区間) のまとめ。scripts/v241_lma_err_s<seed>_t<nturn>.{gs,pysad,pysadk64}.txt を読む。
区間 (lo, hi) が同じなら 'same'、1 つ隣 (幅 0.02/512 = 3.9e-5) なら 'adj'、それ以上は 'diff'。"""
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SP = os.path.join(HERE, 'scripts')
BIN = 0.02 / 512


def parse_b(path):
    if not os.path.exists(path):
        return None
    out = {}
    for l in open(path, encoding='utf-8', errors='replace'):
        if l.startswith('B '):
            v = l.split()
            out[(int(v[1]), int(v[2]))] = (float(v[3]), float(v[4]))
    return out or None


def compare(a, b):
    k = sorted(set(a) & set(b))
    d = np.array([round((a[x][0] - b[x][0]) / BIN) for x in k])
    return k, d


def main(seeds=(1, 2, 3, 4), nturn=2624):
    tot = {}
    allg, allp = [], []
    for s in seeds:
        base = os.path.join(SP, 'v241_lma_err_s%d_t%d' % (s, nturn))
        g, p, k = parse_b(base + '.gs.txt'), parse_b(base + '.pysad.txt'), parse_b(base + '.pysadk64.txt')
        if g is None or p is None:
            continue
        for nm, o in (('v2', p), ('k64', k)):
            if o is None:
                continue
            keys, d = compare(g, o)
            t = tot.setdefault(nm, [0, 0, 0])
            t[0] += int((d == 0).sum()); t[1] += int((np.abs(d) == 1).sum()); t[2] += int((np.abs(d) > 1).sum())
            print('seed %d gs vs pysad %s: same %d adj %d diff %d (of %d)' % (s, nm, (d == 0).sum(), (np.abs(d) == 1).sum(), (np.abs(d) > 1).sum(), len(d)))
        allg.append(g); allp.append(p)
    print('合計', tot)
    for sg, name in ((1, 'δ+'), (-1, 'δ-')):
        vg = np.array([abs(0.5 * (g[k][0] + g[k][1])) for g in allg for k in g if k[1] == sg])
        vp = np.array([abs(0.5 * (p[k][0] + p[k][1])) for p in allp for k in p if k[1] == sg])
        print('%s: gs 平均 %.5f 標準偏差 %.5f | pysad(v2) 平均 %.5f 標準偏差 %.5f (n=%d)' % (name, vg.mean(), vg.std(ddof=1), vp.mean(), vp.std(ddof=1), len(vg)))


if __name__ == '__main__':
    main(nturn=int(sys.argv[1]) if len(sys.argv) > 1 else 2624)
