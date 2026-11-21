# -*- coding: utf-8 -*-
"""v2next P2 項目 1: SynchroBetaEmittance と Emittance[] の比 (AzimuthalModes = 3)。実行: py -3.12 example/v2_next/P2_sb_rings.py (出力は P2_sb_rings.out.txt)
リング: 大きなリング A (pysad/tests/kek-ls_v03_12.lat、3 GeV、周長 603 m)、PF リング (example/pf/lattice/v9_16_asym.sad、2.5 GeV)、
ASTRID2 (example/astrid/lattice/ASTRID2_v01_01.sad、0.58 GeV)。K64 (k64 モード) は大きなリング A で 40 倍 (原因: 始点の分散 η = 0.0445 m)、v2 は 1.0084 倍。
v2 の SB は線形の Emittance[] との差が 1% 以内 (A・ASTRID2) か 0.6% (PF)。"""
import io, math, os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
# 環境変数 PYSAD_SAD_PACKAGES (K64 版 SAD の Packages のフォルダ) を指定して実行する
from pysad import ffscmd, compat

RINGS = [('large ring A (kek-ls_v03_12.lat, 3 GeV, 603 m)', os.path.join(ROOT, 'pysad', 'tests', 'kek-ls_v03_12.lat'), 'USE RING;'),
         ('PF ring (v9_16_asym.sad, 2.5 GeV, 187 m)', os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), 'USE RING;'),
         ('ASTRID2 (ASTRID2_v01_01.sad, 0.58 GeV, 45.7 m)', os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'USE RING;')]
TAIL = '''
FFS;
PageWidth = 1999;
%s CELL; RFSW; RAD; CALC;
em = Emittance[]; nus = Abs[(Tunes/.em)[[3]]];
r = SynchroBetaEmittance[{nus, nus, 0.001}, AzimuthalModes->%d];
Print["RES ", (Emittances/.em)[[1]], " ", r[[1,2]], " ", nus];
Exit;
'''
for name, path, use in RINGS:
    src = open(path, encoding='utf-8', errors='replace').read()
    for modes in (3, 12) if 'large' in name else (3,):
        for mode in ('k64', 'v2'):
            compat.set_mode(mode)
            buf = io.StringIO()
            ffscmd.run_sad_string(src + TAIL % (use, modes), out=buf)
            l = [x for x in buf.getvalue().splitlines() if x.startswith('RES ')][0].split()
            em, sb = float(l[1]), float(l[2])
            print('%-52s modes %2d %-4s Emittance[] %.4e  SB %.4e  SB/EM %.4f' % (name, modes, mode, em, sb, sb / em))
