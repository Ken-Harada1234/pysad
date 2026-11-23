# -*- coding: utf-8 -*-
"""elegant の momentum_aperture (ISR つき = 量子励起あり、SYNCH_RAD = 1・ISR = 1) を、位置ごと・種ごとに serial の elegant で走らせる (WSL、最大 4 本同時)。
使い方: PYTHONUTF8=1 py -3.12 example/15_lma_touschek/touschek_elegant/elegant_isr.py astrid2|kekls [nseed] [np] [周回数]
位置は lma_fluc.py と同じ。lte は sad2ele の変換 (K2 = −h K1、EDGE_ORDER = 2) の ISR = 0 を ISR = 1 にしたもの。種は run_setup の random_number_seed。
出力: elegant_isr_<名前>.lma.txt (loc、s、種、δ+、δ−、損失の周)。作業は WSL の ~/l2/ 。終了後に elegant の残りを ps で確かめる (timeout つき)。"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lma'))
sys.path.insert(0, HERE)
from lma_common import *
HERE = os.path.dirname(os.path.abspath(__file__))  # lma_common の HERE (lma/) で上書きしない
import sad2ele as c2e
import elegant_lma as te


def main(name, nseed=6, nproc=4, turns=None):
    from pysad.trksetup import TrackSetup
    from pysad import compat
    compat.set_mode('v2')
    cs = te.CASES[name]
    st, ffs = load(name)
    nel = len(st.els)
    if name == 'astrid2':
        locs = [2, 17, 40, 100]
    else:
        ele = elegant_kekls_locs()
        c = [l for l in sorted(ele) if l < st.nlat // 10]
        locs = [c[1], c[len(c) // 2]]
    ap_sel = set(((l - 2) % nel) + 1 for l in locs)
    path, _, _ = c2e.convert(name, HERE, st.conditions()['U0_eV'], aper=cs['aper'], ap_sel=ap_sel, tag='_L2isr', am_unique=True)
    txt = open(path, encoding='utf-8').read().replace('ISR=0', 'ISR=1')
    open(path, 'w', encoding='utf-8', newline='\n').write(txt)
    lte = os.path.basename(path)
    nturn = cs.get('turns') or st.nturn_for_synchrotron(2.0)
    if name == 'kekls':
        nturn = 491
    if turns:
        nturn = int(turns)
    jlist = sorted(ap_sel)
    occ2loc = {i + 1: ((j % nel) + 1) for i, j in enumerate(jlist)}
    scratch = os.path.join(ROOT, 'local', 'scratch_L2', 'isr_' + name + (('_n%d' % turns) if turns else ''))
    os.makedirs(scratch, exist_ok=True)
    for f in os.listdir(scratch):
        os.remove(os.path.join(scratch, f))
    tags = []
    for k in range(len(jlist)):
        for seed in range(1, nseed + 1):
            ptag = 'L2isr_%s%s_p%03d_s%02d' % (name, ('_n%d' % turns) if turns else '', k, seed)
            t = te.ele_text(name, lte, cs['pmev'], nturn, cs['step'], cs['splits'], cs['div'], cs['limit'], ptag, pattern='AM%03d' % (k + 1))
            t = t.replace('output_mode = 1', 'output_mode = 0, soft_failure = 1').replace('rootname = %s,' % ptag, 'rootname = %s,\n  random_number_seed = %d,' % (ptag, 1000 * seed + 7))
            open(os.path.join(scratch, ptag + '.ele'), 'w', newline='\n').write(t)
            tags.append((k, seed, ptag))
    wd = '~/l2/%s%s' % (name, ('_n%d' % turns) if turns else '')
    cmd = ('source ~/elegant/env.sh; mkdir -p %s; cd %s; rm -f L2isr_*; cp "%s" .; cp "%s"/*.ele .; '
           "ls L2isr_*.ele | xargs -P %d -I{} bash -c 'timeout -k 20 1800 elegant {} > {}.log 2>&1'; echo done" % (
               wd, wd, te.wsl_path(os.path.join(HERE, lte)), te.wsl_path(scratch), nproc))
    t0 = time.time()
    te.wsl(cmd, timeout=14400)
    print('elegant (serial、同時 %d 本、%d 本): %.0f s' % (nproc, len(tags), time.time() - t0))
    print('終了後の elegant・mpirun の残り:', te.leftovers() or 'なし')
    lines = ['# elegant momentum_aperture with ISR (SYNCH_RAD=1, ISR=1), %s, nturn=%d, MAXAMP +-%.1f mm; loc s seed dp+ dp- lostpass+ lostpass-' % (name, nturn, cs['aper'] * 1e3)]
    for k, seed, ptag in tags:
        cmd = ('source ~/elegant/env.sh; cd %s; sdds2stream -col=ElementName,s,deltaPositive,deltaNegative,lostOnPassPositive,lostOnPassNegative %s.mmap' % (wd, ptag))
        p = te.wsl(cmd)
        for ln in p.stdout.splitlines():
            t = ln.split()
            if len(t) >= 6 and t[0].startswith('AM'):
                loc = occ2loc[int(t[0][2:])]
                lines.append('%5d %10.4f %3d %10.6f %10.6f %6d %6d' % (loc, float(t[1]), seed, float(t[2]), float(t[3]), int(t[4]), int(t[5])))
    open(os.path.join(HERE, 'elegant_isr_%s%s.lma.txt' % (name, ('_n%d' % turns) if turns else '')), 'w', encoding='utf-8', newline='\n').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'astrid2', int(sys.argv[2]) if len(sys.argv) > 2 else 6, int(sys.argv[3]) if len(sys.argv) > 3 else 4,
         int(sys.argv[4]) if len(sys.argv) > 4 else None)
