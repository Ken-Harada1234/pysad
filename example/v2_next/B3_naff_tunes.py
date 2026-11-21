# -*- coding: utf-8 -*-
"""v2next B3 項目 5: 追跡の座標列からチューンを取り出す NAFF (pysad/naff.py)。小さなリングと大きなリングで CALC のチューンと比べ、
MAIN の TRACK の wtunem (FFT のピークの放物線補間) と精度を比べる。
使い方 (リポジトリ直下から): py -3.12 example/v2_next/B3_naff_tunes.py     (1024 周、小さな振幅 x = 2e-7 m √β、放射・RF なし)
出力: リングごとに CALC のチューン (端数)、NAFF、wtunem、それぞれの差。"""
import io, math, os, sys
import numpy as np
from pysad import compat, ffscmd, naff
from pysad.maintrack import wtunem
from pysad.track6n import track_particles_fast
compat.set_mode('v2')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
TW = '''FFS;
USE %s; CELL; CALC;
Print["TW ",{Twiss["NX","$$$"],Twiss["NY","$$$"],Twiss["AX",1],Twiss["BX",1],Twiss["AY",1],Twiss["BY",1]}];
Exit;
'''
FODO = '''MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT SF = (L=0.2 K2=20) SD = (L=0.2 K2=-30);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 SF D2 B D2 SD D2 QD D1);
LINE RING = (ST 12*CELL1);
'''
RINGS = [('12 セル FODO (3 GeV、45.6 m、6 極つき)', FODO, 'RING'),
         ('大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m、3563 要素)', open(os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), encoding='utf-8').read(), 'RING'),
         ('PF リング (2.5 GeV、187 m)', open(os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), encoding='utf-8').read(), 'RING'),
         ('ASTRID2 (0.58 GeV、45.7 m)', open(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read(), 'RING')]
fr = lambda a: (a + 0.5) % 1.0 - 0.5
N = 1024
for name, src, ring in RINGS:
    buf = io.StringIO()
    ffs, _ = ffscmd.run_sad_string(src + '\n' + TW % ring, out=buf)
    line = [l for l in buf.getvalue().splitlines() if l.startswith('TW ')][0]
    nx, ny, ax, bx, ay, by = [float(v) for v in line[3:].strip('{}').split(',')]
    nx, ny = nx / (2 * math.pi), ny / (2 * math.pi)
    els, p = ffs.lat.elements, ffs.momentum()
    amp = 2e-7
    cc = np.array([[amp * math.sqrt(bx), 0.0, 0.7 * amp * math.sqrt(by), 0.0, 0.0, 0.0]])
    rec = []
    for _ in range(N):
        rec.append(cc[0].copy())
        cc = track_particles_fast(els, cc, 1, p, 0.0, False, False, False)
    rec = np.array(rec)
    (fx, fy, _fz), _a = naff.tunes_from_coords(rec, twiss=((ax, bx), (ay, by), (0.0, 1.0)))
    z = np.zeros((N, 3), dtype=complex)
    xh, ph = naff.twiss_normalize(rec[:, 0], rec[:, 1], ax, bx)
    yh, qh = naff.twiss_normalize(rec[:, 2], rec[:, 3], ay, by)
    z[:, 0] = xh + 1j * ph
    z[:, 1] = yh + 1j * qh
    f1, _fa = wtunem(z)
    print(name)
    print('   CALC   nu_x = %.10f   nu_y = %.10f   (端数)' % (fr(nx), fr(ny)))
    print('   NAFF   差 %.2e   %.2e' % (fr(fx - nx), fr(fy - ny)))
    print('   wtunem 差 %.2e   %.2e   (FFT の放物線補間: MAIN の TRACK、TunesByTracking)' % (fr(f1[0][0] - nx), fr(f1[1][0] - ny)))
