# -*- coding: utf-8 -*-
"""ASTRID2 (0.58 GeV) の RADPOL 追跡で、P=1 からの緩和が ST の値の 1.58 倍になる原因 = 倍精度の丸め。
使い方 (リポジトリ直下から): py -3.12 example/16_spin_polarization/radpol_roundoff.py [周数]   (既定 1000。1 周 約 20 ms)
出力: 初期の偏極 P0 = 1、0.5、0.25 から始めた追跡の 1 周あたりの dP と、ST の値 (PST − P0)·sdamp の比。P0=1 だけ 1.5 倍になる (増分 < ulp)。"""
import io, math, os, sys
import numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
from pysad import compat, emit6d, ffscmd, track6
compat.set_mode('v2')
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, '..', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), encoding='utf-8').read() + \
    "\nFFS;\nUSE RING;\nCELL; CALC;\nRFSW; RAD; NOFLUC; POL;\nExit;\n"
ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
p = ffs.momentum(); els = ffs.lat.elements
e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True, pol=True)
sp = e['spin']; circ = emit6d._emit_circ(els); trf0 = e.get('trf0', 0.0); cod = np.array(e['orbit_entrance'])
nt = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
PST = 8 / (5 * math.sqrt(3))
n0 = np.array([0.0, 1.0, 0.0])
print('ASTRID2 0.58 GeV: tau_p = %.2f s (sdamp = %.4e /turn)' % (sp['taup'], sp['sdamp']))
for P0 in (1.0, 0.5, 0.25):
    _o, s_ = track6.track_particles(els, np.array([cod]), nt, p, trf0, True, True, False, spins=[P0 * n0], radpol=True, circ=circ)
    dP = (np.array(s_[0]) @ n0 - P0) / nt
    exp = (PST - P0) * sp['sdamp']
    print('P0 = %.2f: dP/turn = %.4e, ST = %.4e, ratio = %.3f' % (P0, dP, exp, dP / exp))
