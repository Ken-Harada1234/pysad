# -*- coding: utf-8 -*-
"""v2.3.0 U2: ASTRID2 (30 度の複合機能偏向磁石 K1 = -1.678 /m^2) の色収差 ξ・∂²ν/∂δ² を、pysad の K1 の定義 'log' (K64)・'lin' (B' 一様) と
elegant (CSBEND K2 = -h K1 と K2 = 0、EDGE_ORDER = 1・2) で比べる。B160 の「∂²ν/∂δ² の y の 3 倍の差が解消するか」の確認。
使い方 (リポジトリ直下から): PYTHONUTF8=1 py -3.12 example/v2_3/U2_astrid2_chrom.py   (標準出力を U2_astrid2_chrom.out.txt に保存した)
  elegant は WSL (~/elegant/env.sh、serial、timeout 300 s)。lte は T1_astrid2_all.lte (T1_sad2ele.py の変換) の B_BH の K2・EDGE_ORDER と KSEXT の K2 だけ書き換える。
  ラティス: example/astrid/lattice/ASTRID2_v01_01.sad (0.58 GeV、45.7 m)。「自然」は 6 極 (SF・SD) の強さを 0 にしたもの。
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
ENV = dict(os.environ, MSYS_NO_PATHCONV='1')
SCR = os.path.join(ROOT, 'local', 'scratch_U2')
DMAX, NPTS = 3e-3, 11


def pysad_side():
    from pysad import compat, resdrive
    from pysad.trksetup import TrackSetup
    import T1_sad2ele as c2e
    compat.set_mode('v2')
    src = open(c2e.CASES['astrid2']['path'], encoding='utf-8').read()
    nat = src.replace('K2 =-11.729486011842393', 'K2 = 0').replace('K2 =8.681418961367605', 'K2 = 0')
    # DISFRIN=1: 四極・偏向磁石・6 極の端のフリンジ (K64 の既定は Lee-Whiting 型の四極のフリンジなど) を切る。elegant の KQUAD・CSBEND (HGAP = 0)・KSEXT にはフリンジがない
    nof = lambda t: (t.replace('QF3 = ( L =0.1800', 'QF3 = ( DISFRIN=1 L =0.1800').replace('QF2 = ( L =0.1800', 'QF2 = ( DISFRIN=1 L =0.1800')
                      .replace('BH = (L = 0.8500', 'BH = (DISFRIN=1 L = 0.8500').replace('SF      =(L =.075', 'SF      =(DISFRIN=1 L =.075')
                      .replace('SD      =(L =.15', 'SD      =(DISFRIN=1 L =.15'))
    assert nof(src) != src
    out = {}
    for lab, s in (('full', src), ('natural', nat), ('full_nofringe', nof(src)), ('natural_nofringe', nof(nat))):
        st, ffs = TrackSetup.from_source(s, 'RING', rad=False, rfsw=False, alost=0.01)
        for d in ('log', 'lin'):
            compat.set_bend_k1_def(d)
            r = resdrive.chromaticity_numeric(ffs.lat, dmax=DMAX)
            out[(lab, d)] = (r['d1'], r['d2'])
    compat.set_bend_k1_def(None)
    return out


def elegant_side():
    base = open(os.path.join(HERE, 'T1_astrid2_all.lte'), encoding='utf-8').read()
    os.makedirs(SCR, exist_ok=True)
    wp = os.path.abspath(SCR).replace(os.sep, '/')
    wp = '/mnt/' + wp[0].lower() + wp[2:]
    res = {}
    for lab in ('full', 'natural'):
        for k2lab, eo in (('log', 2), ('lin', 2), ('lin_eo1', 1)):
            s = base
            m = re.search(r'B_BH: CSBEND, (.*)', s)
            line = m.group(0)
            new = re.sub(r'K2=[-0-9.e+]+', 'K2=%s' % ('1.03376922564706' if k2lab == 'log' else '0'), line)
            new = re.sub(r'EDGE_ORDER=\d', 'EDGE_ORDER=%d' % eo, new)
            s = s.replace(line, new)
            if lab == 'natural':
                s = re.sub(r'(S_S[FD]: KSEXT, L=[0-9.]+, )K2=[-0-9.e+]+', r'\1K2=0', s)
            name = 'u2_%s_%s' % (lab, k2lab)
            open(os.path.join(SCR, name + '.lte'), 'w', newline='\n').write(s)
            ele = ('&run_setup lattice=%s.lte, use_beamline=RING, p_central_mev=580, default_order=3, rootname=%s, &end\n&run_control &end\n'
                   '&twiss_output filename=%%s.twi, radiation_integrals=1, higher_order_chromaticity=1, higher_order_chromaticity_points=%d, '
                   'higher_order_chromaticity_range=%g &end\n&bunched_beam n_particles_per_bunch=1 &end\n&track &end\n' % (name, name, NPTS, DMAX))
            open(os.path.join(SCR, name + '.ele'), 'w', newline='\n').write(ele)
            cmd = ('source ~/elegant/env.sh; mkdir -p ~/u2/chrom; cd ~/u2/chrom; cp "%s/%s.lte" "%s/%s.ele" .; timeout 280 elegant %s.ele > %s.log 2>&1; '
                   'sddsprintout %s.twi -par=dnux/dp -par=dnuy/dp -par=dnux/dp2 -par=dnuy/dp2') % (wp, name, wp, name, name, name, name)
            p = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', cmd], capture_output=True, text=True, timeout=320, env=ENV)
            v = re.findall(r'dnu[xy]/dp2?[^=]*=\s*([-0-9.e+]+)', p.stdout)
            res[(lab, k2lab)] = (float(v[0]), float(v[1]), float(v[2]), float(v[3]))
    ps = subprocess.run(['wsl', '-d', 'Ubuntu', '--', 'bash', '-lc', 'ps -eo pid,etime,cmd | grep -E "elegant|sad1" | grep -v grep'], capture_output=True, text=True, env=ENV)
    print('# elegant の残り (ps):', repr(ps.stdout.strip()))
    return res


def main():
    P = pysad_side()
    print('# ASTRID2 色収差 ξ = ∂ν/∂δ、∂²ν/∂δ² (pysad は CALC の δ = ±%g の多項式フィット、elegant は higher_order_chromaticity)' % DMAX)
    print('# 6 極つき (full) と 6 極なし (natural)。pysad: log = K64 の K1 の定義、lin = B\' 一様 (v2 のオプション)。elegant: K2 = -hK1 (log に相当)、K2 = 0 (lin に相当)')
    E = elegant_side()
    for lab in ('full', 'natural', 'full_nofringe', 'natural_nofringe'):
        print('\n[%s]  (νx, νy) の ξ | ∂²ν/∂δ²' % lab)
        for k in ('log', 'lin'):
            d1, d2 = P[(lab, k)]
            print('  pysad %-3s            ξ = (%9.4f, %9.4f)  d2 = (%8.2f, %8.2f)' % (k, d1[0], d1[1], d2[0], d2[1]))
        for k in ('log', 'lin', 'lin_eo1'):
            e = E[(lab.split('_')[0], k)]
            print('  elegant %-8s EO=%d  ξ = (%9.4f, %9.4f)  d2 = (%8.2f, %8.2f)' % (k, 1 if k.endswith('eo1') else 2, e[0], e[1], e[2], e[3]))
        for k in ('log', 'lin'):
            d1, d2 = P[(lab, k)]
            e = E[(lab.split('_')[0], k)]
            print('  比 (pysad %s / elegant %s EO=2): d2 = (%.3f, %.3f)、ξ の差 = (%.4f, %.4f)' % (k, k, d2[0] / e[2], d2[1] / e[3], d1[0] - e[0], d1[1] - e[1]))


if __name__ == '__main__':
    main()
