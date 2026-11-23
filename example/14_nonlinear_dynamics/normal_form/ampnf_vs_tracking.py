# -*- coding: utf-8 -*-
"""ampnf (1 周の写像のべき級数の正規形) の ν(Jx, Jy, δ) を、(1) 同じ模型の数値の追跡、(2) pysad の追跡 (k64、実際の素子) の NAFF と比べる。
使い方: py -3.12 example/14_nonlinear_dynamics/normal_form/ampnf_vs_tracking.py fodo|astrid  (出力は標準出力。保存は ampnf_vs_tracking_*.out.txt)
(2) は、追跡の初期座標の J (Twiss の J) と正規形の J が O(√J) 違うので、追跡の初期座標を正規形の J に直してから (ampnf.nf_action) 級数を評価する。
条件: pysad k64、RF 止め・放射なし・結合なし、δ = 0 と δ = ±1e-3 (閉軌道のまわり)、1024 周、NAFF Hann p=1。"""
import math
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import adtsdirect as D, ampnf as A, compat, resdrive as R
from pysad.parser import load_sad_file, load_sad_string
compat.set_mode('k64')


def load(name):
    if name == 'fodo':
        src = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
        src = src.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
        return load_sad_string(src).build_lattice('RING'), 3e9, [(2e-7, 2e-7), (1e-6, 1e-6), (3e-6, 1e-6)]
    if name == 'astrid':
        lat = load_sad_file(os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad')).build_lattice('RING')
        return lat, 0.58e9, [(2e-9, 2e-9), (1e-8, 1e-8), (3e-8, 1e-8)]
    raise SystemExit('fodo|astrid')


def P(*a):
    print(*a, flush=True)


BEND = True
NSB = 8


def main(name, deltas=(0.0, 1e-3, -1e-3), nturn=1024, Dord=6):
    lat, p0, amps = load(name)
    t0 = A.amplitude_tune(lat, D=Dord, nslice=8, bend_cubic=BEND, nslice_bend=NSB, kin=True, fringe=True)
    P('# %s: ampnf D=%d  nu0 = (%.6f, %.6f)  最小の共鳴の分母 %.3g  警告 %s' % (name, Dord, t0.nu0[0], t0.nu0[1], t0.min_den, t0.warnings))
    an = R.adts(lat)
    P('J1: ampnf xx,xy,yy = %.3f %.3f %.3f | adts().sext+octu = %.3f %.3f %.3f | adts().total = %.3f %.3f %.3f'
      % (t0.adts['xx'], t0.adts['xy'], t0.adts['yy'], an.sext['xx'] + an.octu['xx'], an.sext['xy'] + an.octu['xy'], an.sext['yy'] + an.octu['yy'],
         an.total['xx'], an.total['xy'], an.total['yy']))
    P('J2 (c20,c11,c02 of nu_x | nu_y):', ['%.4g' % t0.adts2[k] for k in ('x20', 'x11', 'x02', 'y20', 'y11', 'y02')])
    P('chrom xi1 (x,y) = %.4f %.4f  (6 極 + 自然 + 分散。偏向磁石の曲率の項は入れていない)' % (t0.chrom(0, 1), t0.chrom(1, 1)))
    ex = {k: 0.0 for k in ('xx', 'xy', 'yy')}
    P('追跡 (pysad k64) との比較: Δν = ν_追跡 − ν_級数 (模型 = 線形 + 6 極・8 極 + 偏向磁石の曲率・端 + 運動学的な項 + 四極・K1 の偏向磁石の端部のフリンジ、D = %d)' % Dord)
    P('δ, (Jx0, Jy0), plane, J_nf, Δν = 追跡 − 級数, J² の寄与, Δν / J² の寄与')
    for d in deltas:
        rows = lat.twiss_table(delta=d)
        tw, orb, g, nu = D.lattice_start(lat, d, rows)
        for (jx, jy) in amps:
            c = D.initial_coords(tw, [jx], [jy], [0.3], [0.9], delta=d, orbit=orb)
            rec = D.track_record(lat.elements, p0, c, nturn)
            nuh = D.record_tunes(rec, tw, guess=g, halfwidth=0.05, orbit=orb)[0]
            jxn, jyn, _, _, err = A.nf_action(t0, c[0, :4] - 0 * np.array(orb), d)
            for pl in (0, 1):
                ser = t0.nu(pl, jxn, jyn, d)
                f0 = D.fr(nuh[pl] - g[pl]) + g[pl]
                ser_f = D.fr(ser - g[pl]) + g[pl]
                add = (ex['xx'] * jxn + ex['xy'] * jyn) if pl == 0 else (ex['xy'] * jxn + ex['yy'] * jyn)
                j2 = sum(t0.coef(pl, a, b, 0) * jxn ** a * jyn ** b for a, b in ((2, 0), (1, 1), (0, 2)))
                P('%+.0e (%.0e,%.0e) %s J=(%.3e,%.3e) dnu=%+.3e J2=%+.3e dnu/J2=%+.2f' % (d, jx, jy, 'xy'[pl], jxn, jyn, f0 - ser_f, j2, (f0 - ser_f) / j2 if j2 else 0))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fodo')
