# -*- coding: utf-8 -*-
"""v2.4.1: 誤差つき ASTRID2 の LMA を gs と pysad で比べる (同じ誤差の値を台本に書く)。
gs の台本: TrackParticles を 1 粒子ずつ呼ぶ二分法 (9 回、区間 [0, 0.02]、RFSW・RAD ON、FLUC なし、LOSSAMPL・LOSSDZ)。同じ台本を pysad の台本評価器にも渡す。
pysad の API (lma.momentum_acceptance、check_every=0 = tapert と同じ判定) の区間とも比べる。
使い方: py -3.12 pysad/tests/script_tests/_v241_lma_errors_run.py SEED NTURN [--nogs] [--locs 2,17,...] [--id 番号]
出力: pysad/tests/script_tests/scripts/v241_lma_err_s<SEED>_t<NTURN>.{sad,pysad.txt,gs.txt}、標準出力 (表)
"""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _v241_lma_tools as T

APER = 0.01
LOCS = [2, 17, 40, 75, 110, 150]


def make_script(seed, nturn, locs, level=1.0, rad=True):
    from pysad import compat
    from pysad.trksetup import TrackSetup
    compat.set_mode('v2')
    src = T.read(T.ASTRID2)
    text, info = T.error_source(src, seed, level)
    st, ffs = TrackSetup.from_source(text, 'RINGERR', rad=True, rfsw=True, alost=APER)
    zb = st.bucket_zlimit()
    L = ['FFS;', 'USE RINGERR; CELL; CALC; RFSW; RAD; RADCOD; NOFLUC;', 'Emittance[];', 'LOSSAMPL = %.17g; LOSSDZ = %.17g;' % (APER, zb), '$FORM = "20.14";', 'NPARA = 1;']
    for loc in locs:
        cod = st.closed_orbit([loc])[0]
        for sg in (1, -1):
            L.append('cod = {%s}; lo = 0; hi = 0.02;' % ','.join('%.17g' % v for v in cod))
            L.append('Do[mid = (lo + hi)/2; b1 = {{cod[[1]]}, {cod[[2]]}, {cod[[3]]}, {cod[[4]]}, {cod[[5]]}, {cod[[6]] + (%d)*mid}, {1}};' % sg)
            L.append('  r = TrackParticles[{%d, b1}, %d, 1, %d]; If[r[[2,7,1]] == 1, lo = mid, hi = mid], {k, 1, 9}];' % (loc, loc, nturn))
            L.append('Print["B %d %d ", lo, " ", hi];' % (loc, sg))
    L.append('Exit;')
    return text + '\n' + '\n'.join(L) + '\n', st, info, zb


def parse_b(t):
    return {(int(l.split()[1]), int(l.split()[2])): (float(l.split()[3]), float(l.split()[4])) for l in t.splitlines() if l.startswith('B ')}


def main():
    a = sys.argv[1:]
    seed, nturn = int(a[0]), int(a[1])
    locs = [int(v) for v in a[a.index('--locs') + 1].split(',')] if '--locs' in a else LOCS
    uid = a[a.index('--id') + 1] if '--id' in a else '1'
    tag = 'v241_lma_err_s%d_t%d' % (seed, nturn)
    script, st, info, zb = make_script(seed, nturn, locs)
    sad = os.path.join(T.SP, tag + '.sad')
    open(sad, 'w', encoding='utf-8').write(script)
    print('seed %d、誤差 %s、nturn %d、LOSSDZ %.4f m' % (seed, info, nturn, zb), flush=True)
    py_txt, tpy = T.run_pysad_text(script, os.path.join(T.SP, tag + '.pysad.txt'))
    print('pysad の台本の二分法: %.0f s' % tpy, flush=True)
    from pysad.lma import momentum_acceptance
    import time
    t0 = time.time()
    res = momentum_acceptance(st, locs, nturn=nturn, dpmax=0.02, ncoarse=8, nsub=3, tol=2e-5, alost=APER, check_every=0, workers=1, quiet=True, turncheck=0)
    print('LMA API: %.0f s' % (time.time() - t0), flush=True)
    gtxt = None
    if '--nogs' not in a:
        gsname = os.path.join(T.SP, '_gs_lma_%s_%s.sad' % (uid, tag))
        open(gsname, 'w', encoding='utf-8').write(script)
        gtxt, tgs, left = T.run_gs(gsname, os.path.join(T.SP, tag + '.gs.txt'))
        if os.path.exists(gsname):
            os.remove(gsname)
        print('gs: %.0f s、自分の台本の sad1.exe の残り: %s' % (tgs, left or 'なし'), flush=True)
    bpy, bgs = parse_b(py_txt), (parse_b(gtxt) if gtxt else {})
    print('loc sg | pysad 台本 (lo, hi) | gs 台本 (lo, hi) | LMA API (lo, hi)')
    for i, loc in enumerate(res.loc):
        for sg, lo, hi in ((1, res.dp_plus[i], res.hi_plus[i]), (-1, -res.dp_minus[i], -res.hi_minus[i])):
            p = bpy[(loc, sg)]
            g = bgs.get((loc, sg))
            print('%4d %+d | %.6f %.6f | %s | %.6f %.6f' % (loc, sg, p[0], p[1], ('%.6f %.6f' % g) if g else '-', lo, hi))


if __name__ == '__main__':
    main()
