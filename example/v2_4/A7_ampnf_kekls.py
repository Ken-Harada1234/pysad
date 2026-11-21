# -*- coding: utf-8 -*-
"""v2.4.0 A7: ampnf の大きなリング (kek-ls_v03_12.lat、3 GeV、603 m) の全部入り (bend_cubic・kin・fringe) の計測と、ADTS・ξ の解析式との比較。
使い方: PYTHONUTF8=1 py -3.12 example/v2_4/A7_ampnf_kekls.py [D]"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
from pysad import ampnf as A, resdrive as R, compat
from pysad.parser import load_sad_file
compat.set_mode('v2')
D = int(sys.argv[1]) if len(sys.argv) > 1 else 5
lat = load_sad_file(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat')).build_lattice('RING')
print('要素 %d' % len(lat.elements), flush=True)
t = time.time()
r = A.amplitude_tune(lat, D=D, nslice=4, nslice_bend=4, nslice_quad=2, bend_cubic=True, kin=True, fringe=True)
print('ampnf D=%d 全部入り %.0f s、ν0 = (%.5f, %.5f)、min_den %.3g' % (D, time.time() - t, r.nu0[0], r.nu0[1], r.min_den), flush=True)
print('J1 xx,xy,yy', r.adts['xx'], r.adts['xy'], r.adts['yy'])
an = R.adts(lat)
print('adts().total', an.total['xx'], an.total['xy'], an.total['yy'])
print('J2', r.adts2)
print('xi1,2,3 x', [r.chrom(0, p) for p in (1, 2, 3)], 'y', [r.chrom(1, p) for p in (1, 2, 3)])
ref = R.chromaticity_numeric(lat, order=3)
print('CALC', ref['d1'], ref['d2'], ref['d3'])
