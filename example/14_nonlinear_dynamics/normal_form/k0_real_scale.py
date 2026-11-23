# -*- coding: utf-8 -*-
"""実機規模のリングで K0≠0 の偏向磁石 (多極ウィグラー MPW) つきの ν(J, δ)・ADTS を ampnf で求め、実際の追跡 + NAFF と比べる。
  astrid  : ASTRID2 (0.58 GeV、周長 45.7 m、222 要素) + 2 T の MPW (K0 = ±0.0640435 rad の 23 極。example/astrid/lattice/ASTRID2_v01_02_id.sad の RINGMPW)。
            tilt = MPW の全極を 20 mrad 傾ける (ROTATE)。
  kekls   : KEK-LS (3 GeV、周長 603 m、3563 要素。pysad/tests/kek-ls_v03_12.lat) の 1 つ目のセルの長直線部に、K0 = ±0.00875 rad (3 GeV で約 1 T、極長 0.0875 m) の
            23 極の MPW (合計 1.93 m。K0 の和 = 0、軌道は局所のバンプ) を挿入したもの。
使い方: PYTHONIOENCODING=utf-8 py -3.12 example/14_nonlinear_dynamics/normal_form/k0_real_scale.py [astrid | astrid_tilt | kekls] [D]"""
import math
import os
import sys
import time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import ampnf as A, resdrive as R, compat
from pysad.parser import load_sad_file, load_sad_string


def P(*a):
    print(*a, flush=True)


def astrid(tilt=0.0):
    f = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
    s = open(f, encoding='utf-8', errors='replace').read()
    if tilt:
        s = s.replace('BMPW1I = ( L = .175/4 K0 =  0.0640435140143602/2           E2 =1)', 'BMPW1I = ( L = .175/4 K0 =  0.0640435140143602/2           E2 =1 ROTATE=%r)' % tilt)
        s = s.replace('BMPW1  = ( L = .175/2 K0 =  0.0640435140143602    E1 =.5    E2 =.5)', 'BMPW1  = ( L = .175/2 K0 =  0.0640435140143602    E1 =.5    E2 =.5 ROTATE=%r)' % tilt)
        s = s.replace('BMPW2  = ( L = .175/2 K0 = -0.0640435140143602    E1 =.5    E2 =.5)', 'BMPW2  = ( L = .175/2 K0 = -0.0640435140143602    E1 =.5    E2 =.5 ROTATE=%r)' % tilt)
        s = s.replace('BMPW1F = ( L = .175/4 K0 =  0.0640435140143602/2  E1 =1)', 'BMPW1F = ( L = .175/4 K0 =  0.0640435140143602/2  E1 =1 ROTATE=%r)' % tilt)
        assert s.count('ROTATE=%r' % tilt) == 4
    return load_sad_string(s).build_lattice('RINGMPW'), 0.58e9, 5e-9


def kekls(k=0.00875, lp=0.0875):
    f = os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat')
    s = open(f, encoding='utf-8', errors='replace').read()
    mpw = (' BEND MPI=(L=%r K0=%r) MP1=(L=%r K0=%r) MP2=(L=%r K0=%r) MPF=(L=%r K0=%r);\n'
           ' LINE MPWL=(MPI MP2 10*(MP1 MP2) MPF);\n' % (lp / 2, k / 2, lp, k, lp, -k, lp / 2, k / 2))
    total = lp / 2 + lp + 20 * lp + lp / 2
    drida = 2.85 - total
    mpw += ' DRIFT DRIDA=(L=%r);\n' % drida
    # 1 つ目の ARCAH の先頭の DRID1 を DRIDA + MPWL に置き換えた ARCAHM を作り、RING の 1 つ目の ARCA を ARCAM にする
    i0 = s.index('LINE' + chr(10) + '  ARCAH    =(MKS')
    i = s.index('ARCAH    =(MKS')
    j = s.index('ARCAH2    =(MKS')
    bodym = s[i:j].replace('ARCAH    =(MKS', 'ARCAHM    =(MKS', 1).replace('DRID1 BPM DRID2', 'DRIDA MPWL BPM DRID2', 1).rstrip()
    nl = chr(10)
    s = s[:i0] + mpw + 'LINE' + nl + '  ' + bodym + ';' + nl + s[i0:]
    s = s.replace(' ARCA      = (ARCAH -ARCAH2)', ' ARCA      = (ARCAH -ARCAH2)' + nl + ' ARCAM     = (ARCAHM -ARCAH2)', 1)
    old_ring = 'ARCA, ARCA, ARCA, ARCA,' + nl + '             ARCA, ARCA, ARCA, ARCA,' + nl + '             ARCA, ARCA, ARCA, ARCA, RFC,'
    assert old_ring in s
    s = s.replace(old_ring, old_ring.replace('ARCA, ARCA, ARCA, ARCA,', 'ARCAM, ARCA, ARCA, ARCA,', 1), 1)
    assert 'ARCAM, ARCA' in s
    return load_sad_string(s).build_lattice('RING'), 3e9, 1e-8


def main():
    kind = sys.argv[1] if len(sys.argv) > 1 else 'astrid'
    D = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    compat.set_mode('k64')
    if kind == 'astrid':
        lat, p0, j0 = astrid()
    elif kind == 'astrid_tilt':
        lat, p0, j0 = astrid(0.02)
    else:
        lat, p0, j0 = kekls()
    mp = [e for e in lat.elements if type(e).__name__ == 'Bend' and float(getattr(e, 'K0', 0) or 0) != 0]
    P('## %s、D = %d: 要素 %d、K0≠0 の BEND %d 台、運動量 %.2f GeV' % (kind, D, len(lat.elements), len(mp), p0 / 1e9))
    opt = dict(nslice=8, nslice_bend=16, bend_cubic=True, kin=True, fringe=True)       # ASTRID2 の BH (0.85 m、0.52 rad) は nslice_bend = 4 だと δ = ±1e-3 の ν が 1e-5 外れる (O(1/n²)。16 で 3e-7)
    if kind == 'kekls':
        opt.update(nslice=4, nslice_bend=4, nslice_quad=2)
    t = time.time()
    r = A.amplitude_tune(lat, D=D, **opt)
    P('  ampnf 全部入り %.0f s、ν0 = (%.5f, %.5f)、最小の分母 %.3g、閉軌道 %.1e、結合 %.1e' % (time.time() - t, r.nu0[0], r.nu0[1], r.min_den, np.abs(r.res['orbit']).max(), r.res['coupling']))
    P('  J¹ = xx %.1f xy %.1f yy %.1f、ξ1 = %.4f %.4f' % (r.adts['xx'], r.adts['xy'], r.adts['yy'], r.chrom(0), r.chrom(1)))
    for (jx, jy, d) in ((j0, j0, 0.0), (3 * j0, j0, 0.0), (j0, j0, 1e-3), (j0, j0, -1e-3)):
        t = time.time()
        tr = A.tracked_tunes(lat, r, jx, jy, delta=d, momentum_eV=p0)
        P('  (Jx, Jy, δ) = (%g, %g, %g): |追跡 − 級数| = %.1e (νx)、%.1e (νy) (追跡 %.0f s)' % (jx, jy, d, abs(tr[0] - r.nu(0, jx, jy, d)), abs(tr[1] - r.nu(1, jx, jy, d)), time.time() - t))
    t = time.time()
    tr = A.tracked_adts(lat, r, j=j0, momentum_eV=p0)
    P('  J¹ (追跡 + NAFF のフィット) xx %.1f xy %.1f yy %.1f、ampnf との相対差 %.1e %.1e %.1e (%.0f s)' % (
        tr['xx'], tr['xy'], tr['yy'], *(abs(r.adts[k] - tr[k]) / abs(tr[k]) for k in ('xx', 'xy', 'yy')), time.time() - t))
    t = time.time()
    a = R.adts(lat)
    P('  adts().total xx %.1f xy %.1f yy %.1f (%.0f s)、追跡との相対差 %.1e %.1e %.1e' % (
        a.total['xx'], a.total['xy'], a.total['yy'], time.time() - t, *(abs(a.total[k] - tr[k]) / abs(tr[k]) for k in ('xx', 'xy', 'yy'))))
    if hasattr(a, 'rotate') and getattr(a, 'rolled', None):
        P('  .rotate = %s' % {k: round(float(v), 2) for k, v in a.rotate.items()})


if __name__ == '__main__':
    main()
