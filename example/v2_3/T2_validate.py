# -*- coding: utf-8 -*-
"""v2.3.0 T2: touschek_piwinski (1 次元の式) を、独立な 2 つの数値積分 (T2_check3d: 3 次元の一般の積分、T2_mc: 実験室系の対のモンテカルロ、厳密な相対論) と、スライスごとに比べる。
使い方: PYTHONUTF8=1 py -3.12 example/v2_3/T2_validate.py astrid2|kekls [pbunch:mincoup ...]   (T2_validate_<名前>.out.txt に保存した)
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)


def main():
    name = sys.argv[1]
    sets = [tuple(float(v) for v in s.split(':')) for s in sys.argv[2:]] or [(1e9, 0.01), (1e9, 1.0)]
    amps = [0.004, 0.016] if name == 'astrid2' else [0.004, 0.020]
    from T2_scan import gather
    from pysad.touschek import piwinski_slice_rate
    from T2_check3d import slice_rate3d
    from T2_mc import mc_rate
    for pb, mc in sets:
        tt, tt2, ref = gather(name, pb, mc)
        info = tt2['rest_info']
        rest = tt2['rest']
        n = len(rest)
        picks = [n // 7, n // 3, (2 * n) // 3]
        print('\n%s PBUNCH=%g MINCOUP=%g (εy/εx = %.3f, σδ = %.3e): スライス %s' % (name, pb, mc, ref['EMITYC'] / ref['EMITXC'], ref['SIGE'], picks), flush=True)
        print('  slice  l   a      1D 式 [1/s]   3D (±)            3D/1D     MC 厳密運動学 (±)    MC/1D')
        for j in picks:
            s = rest[j]
            for a in amps:
                r1 = piwinski_slice_rate(s, a, a, info['circ'], pb)
                r3, e3 = slice_rate3d(s['al'], None, s['pl'], s['p0'], s['h0'], s['xx'], info['circ'], pb, a, a, nu=4096)
                rm, em = mc_rate(s['lab'], s['al'], info['circ'], pb, s['p0'], a, a, npairs=8192)
                print('  %5d %4d %.4f %12.5e  %12.5e (%.0e)  %.4f   %12.5e (%.0e)   %.4f' % (j, s['l'], a, r1, r3, e3, r3 / r1, rm, em, rm / r1), flush=True)


if __name__ == '__main__':
    main()
