# -*- coding: utf-8 -*-
"""1 本の共鳴の島 (包絡の理論 `ampnf_ext.island_theory`) を追跡の島 (`ampnf_ext.track_island`) と比べる。
リング: FODO (3 GeV、105.6 m、24 セル、example/00_lattice/FODO_ring.sad) の 24 セルのうち 12 セルだけに 6 極 (K2 = 4・−6、SF・SD) を置いた非周期の配置
(AABABBABAABBBAABABBAABBA。24 回対称のままだと 24 の倍数でない共鳴 nνx = p は禁止され、駆動項は 1e-4 以下になるので島が作れない)。QF の K1 でチューンを調整。
面: 水平 (y = py = 0、δ = 0、結合なし 1 自由度)。共鳴: n νx = p (n = 4・5・6)、作動点は共鳴の下側 δν = p/n − νx。
理論: ampnf (D = 6、bend_cubic・kin・fringe) の正規形から共鳴の単項式を 1 本残し、固定点の作用 J_r と島の半幅を出す。
追跡: pysad (k64、RF・放射なし) の周期 n の固定点 (O・X) と X から出る不安定多様体の J の範囲 (正規形の作用)。
使い方: py -3.12 example/14_nonlinear_dynamics/normal_form/island_vs_tracking.py [n:dnu ...]  (既定 6:0.001 6:0.002 6:0.004 6:0.008 5:0.001 5:0.002 5:0.004 5:0.008 4:0.004)"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import ampnf as A, ampnf_ext as E, compat
from pysad.parser import load_sad_string

compat.set_mode('k64')
BASE = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
PAT = 'AABABBABAABBBAABABBAABBA'
OPT = dict(bend_cubic=True, kin=True, fringe=True, nslice_bend=8)


def mk(kf, kd=-0.6, k2f=4.0, k2d=-6.0, pat=PAT):
    src = BASE.replace('SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=%g) SD = (L=0.2 K2=%g) SF0 = (L=0.2 K2=0) SD0 = (L=0.2 K2=0)' % (k2f, k2d))
    src = src.replace('K1=0.7', 'K1=%r' % kf).replace('K1=-0.6', 'K1=%r' % kd)
    src = src.replace('LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);', 'LINE CELLA = (SF DS QF D2 B D2 QD DS SD D1); LINE CELLB = (SF0 DS QF D2 B D2 QD DS SD0 D1);')
    src = src.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF %s);' % ' '.join('CELL' + c for c in pat))
    return load_sad_string(src).build_lattice('RING')


def nux(kf):
    return A.amplitude_tune(mk(kf), D=3, nslice=2).nu0[0]


def find_kf(target, k0):
    k = k0
    for _ in range(10):
        v = nux(k)
        sl = (nux(k + 1e-4) - v) / 1e-4
        e = ((target - v + 0.5) % 1) - 0.5
        k += e / sl
        if abs(e) < 1e-7:
            break
    return k


KF0 = {4: 0.7038, 5: 0.6137, 6: 0.7115}
TARGET = {4: 3 / 4, 5: 3 / 5, 6: 5 / 6}


def run(n, dnu, nturn_sep=300):
    kf = find_kf(TARGET[n] - dnu, KF0[n])
    lat = mk(kf)
    t0 = time.time()
    r = A.amplitude_tune(lat, D=6, nslice=8, **OPT)
    al = E.pure_alpha(n, 0)
    nf = E.resonance_normal_form(r.res, 0, al)
    th = E.island_theory(r, nf, 0, al)
    print('## n = %d  δν = %.4f  kf = %.6f  ν0 = (%.5f, %.5f)  xx = %.1f  x20 = %.3g  最小の分母 %.3f  (ampnf %.0f s)' % (n, dnu, kf, r.nu0[0], r.nu0[1], r.adts['xx'], r.adts2['x20'], r.min_den, time.time() - t0), flush=True)
    if not th['ok']:
        print('  理論: ', th['reason'])
        return
    print('  理論: J_r = %.4e  駆動項 |r| = %.3e  A = %.3e  半幅 ΔJ = %.3e (%.1f%% of J_r)  [J_lo, J_hi] = [%.4e, %.4e]  振り子の近似 %s'
          % (th['Jr'][0], abs(th['r']), th['A'], th['half_J'][0], 100 * th['rel_width'], th['Jlo'], th['Jhi'], 'OK' if th['pendulum_ok'] else '外れる (幅が J_r の 30% 超)'), flush=True)
    t1 = time.time()
    tr = E.track_island(lat, r, n, 0, 3e9, jr=th['Jr'][0], nturn_sep=nturn_sep)
    Jo = np.mean([j for _, j, _ in tr['O']]) if tr['O'] else float('nan')
    Jx = np.mean([j for _, j, _ in tr['X']]) if tr['X'] else float('nan')
    print('  追跡: O %d 個 (J 平均 %.4e)・X %d 個 (J 平均 %.4e)  X の跡 %.4f' % (tr['n_O'], Jo, tr['n_X'], Jx, tr['X'][0][2] if tr['X'] else float('nan')), flush=True)
    if 'sep' in tr:
        print('  追跡の分離線の J の範囲 [%.4e, %.4e]  幅/2 = %.3e  (理論の半幅の %.3f 倍)  脱出 %s  (追跡 %.0f s)'
              % (tr['Jmin'], tr['Jmax'], 0.5 * (tr['Jmax'] - tr['Jmin']), 0.5 * (tr['Jmax'] - tr['Jmin']) / th['half_J'][0], tr['escaped'], time.time() - t1), flush=True)
        print('  比: J_r の理論/追跡 (O・X の平均) = %.4f、J_lo の理論/追跡 = %.3f、J_hi の理論/追跡 = %.3f'
              % (th['Jr'][0] / (0.5 * (Jo + Jx)), th['Jlo'] / tr['Jmin'], th['Jhi'] / tr['Jmax']), flush=True)


if __name__ == '__main__':
    args = sys.argv[1:] or ["6:0.001", "6:0.002", "6:0.004", "6:0.008", "5:0.001", "5:0.002", "5:0.004", "5:0.008", "4:0.004"]
    for a in args:
        n, d = a.split(':')
        run(int(n), float(d))
