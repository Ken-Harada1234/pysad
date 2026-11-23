# -*- coding: utf-8 -*-
"""局所の運動量アクセプタンス (LMA) の Touschek 寿命を、K64 の表と Piwinski 型の式 (touschek_lifetime_local_piwinski) で出す。
LMA の表は lma_ring.py・lma_compare.py の保存済みの表 (ASTRID2: lma_astrid2.table.txt、KEK-LS: pysad_kekls.lma.txt。1 セルを周期的に延ばす)。IBS 込み PBUNCH = 1e9、MINCOUP = 0.01。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/touschek_formula_lma_lifetime.py astrid2|kekls   (touschek_formula_lma_lifetime.out.txt に保存した)
"""
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)


def main(name):
    from pysad import compat
    from pysad.trksetup import TrackSetup, find_periodicity
    from pysad.lma import touschek_lifetime_local, touschek_lifetime_local_piwinski
    from pysad.touschek import touschek_table
    import sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES[name]['path'], encoding='utf-8').read()
    st, ffs = TrackSetup.from_source(src, 'RING', rad=True, rfsw=True)
    g = ffs.I.globals
    ffs.flags['INTRA'] = True
    g['PBUNCH'] = 1e9
    ffs.flags['RFSW'] = True
    ffs.flags['RAD'] = True
    ffs.run_emit(quiet=True)
    tt = ffs.touschek
    per = find_periodicity(st.els, st.rows)
    tp = os.path.join(HERE, *(('..', 'lma', 'lma_astrid2.table.txt') if name == 'astrid2' else ('pysad_kekls.lma.txt',)))
    rows = [ln.split() for ln in open(tp, encoding='utf-8') if ln.strip() and not ln.startswith('#')]
    cs = [r for r in rows if float(r[1]) < per['cell_length'] + 1e-9]
    lm = SimpleNamespace(s=np.array([float(r[1]) for r in cs]), dp_plus=np.array([float(r[3]) for r in cs]), dp_minus=np.array([float(r[4]) for r in cs]), setup=st)
    k = touschek_lifetime_local(tt['table'], lm, setup=st, period=per['cell_length'], s0=per['s_starts'][0])
    p = touschek_lifetime_local_piwinski(tt, lm, setup=st, period=per['cell_length'], s0=per['s_starts'][0])
    print('%s: 局所の LMA (%d 位置、1 セルを周期的に) の寿命: K64 の表 %.1f s、Piwinski 型 %.1f s (Piwinski/K64 − 1 = %+.1f%%)' % (
        name, len(cs), k['lifetime'], p['lifetime'], (p['lifetime'] / k['lifetime'] - 1) * 100))


if __name__ == '__main__':
    main(sys.argv[1])
