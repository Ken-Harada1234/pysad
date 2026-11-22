# -*- coding: utf-8 -*-
"""v2.4.0 G: K0≠0 の偏向磁石 (挿入光源・MPW・傾けた MPW) の ν(J, δ)・ADTS・共鳴駆動項を ampnf (TPSA) で求め、追跡と比べる。
リング: 小さな FODO (3 GeV、105.6 m、24 セル、SF・SD の K2L = 4, -6、example/00_lattice/FODO_ring.sad) の 1 セルの長直線部 (D1) に、3 極の挿入光源 (K0 = a, −2a, a [rad]。
各極 L = 0.3 m、極のあいだ・前後 0.025 m (3 極 + 4 ドリフトで D1 と同じ 1.0 m) のドリフト) を入れる。K0 の和 = 0 で、軌道は局所のバンプ (挿入光源の外に漏れない)。
使い方: py -3.12 example/v2_4/G_k0_ampnf.py [mpw | mpw_big | mpw_ang | mpw_edge | tilt | pi | pi_neg]  (出力は G_k0_ampnf.out.txt、記録は docs/history/v24_G.md)"""
import math
import os
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import ampnf as A, compat
from pysad.parser import load_sad_string

FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read().replace(
    'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)', 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)')
OPT = dict(bend_cubic=True, kin=True, fringe=True, nslice_bend=4)


def pole_defs(a=0.06, ang=0.0, e=(0.0, 0.0), rot=0.0, pat=(1.0, -2.0, 1.0), L=0.3, extra='', k1=0.0):
    out = []
    for i, f in enumerate(pat):
        x = 'L=%r K0=%r' % (L, a * f)
        if ang:
            x += ' ANGLE=%r E1=%r E2=%r' % (ang * f, e[0], e[1])
        if rot:
            x += ' ROTATE=%r RAD' % rot
        if k1:
            x += ' K1=%r' % k1
        x += ' ' + extra
        out.append('P%d = (%s)' % (i + 1, x))
    return out


def build(a=0.06, ang=0.0, e=(0.0, 0.0), rot=0.0, pat=(1.0, -2.0, 1.0), L=0.3, extra='', k1=0.0):
    """FODO の最後のセルの長直線部に 3 極 (K0 = a·pat)。a = 0 かつ ang = 0 なら挿入光源なし (ドリフトに置き換え)"""
    defs = pole_defs(a, ang, e, rot, pat, L, extra, k1)
    s = FODO.replace('LINE CELL =', 'BEND ' + ' '.join(defs) + ';\nDRIFT D1X = (L=0.025);\n'
                     'LINE CELLX = (SF DS QF D2 B D2 QD DS SD D1X P1 D1X P2 D1X P3 D1X);\nLINE CELL =')
    s = s.replace('LINE RING = (ORIGIN RF 24*CELL);', 'LINE RING = (ORIGIN RF 23*CELL CELLX);')
    return load_sad_string(s).build_lattice('RING')


def P(*a):
    print(*a, flush=True)


KINDS = {'mpw': dict(a=0.06), 'mpw_big': dict(a=0.075), 'mpw_ang': dict(a=0.06, ang=0.01), 'mpw_edge': dict(a=0.06, ang=0.01, e=(0.5, 0.5)),
         'tilt': dict(a=0.06, rot=0.2), 'pi': dict(a=0.06, rot=math.pi), 'pi_neg': dict(a=-0.06)}


def run(kind, D=5):
    """ampnf の級数 (D 次) と、pysad の実際の追跡 + NAFF の ν (4 点) と ADTS (J¹ のフィット)、1 周の写像の 2 次のテイラー展開の h (3 次 20 項・δ 9 項)"""
    from pysad import rdttrack as T
    lat = KINDS[kind] and build(**KINDS[kind])
    r = A.amplitude_tune(lat, D=D, nslice=8, **OPT)
    P('## %s %s' % (kind, KINDS[kind]))
    P('  閉軌道 (リング先頭) |c| = %.2e、ν0 = (%.6f, %.6f)、結合 %.1e、最小の分母 %.3f' % (np.abs(r.res['orbit']).max(), r.nu0[0], r.nu0[1], r.res['coupling'], r.min_den))
    P('  J¹ (ampnf) xx xy yy = %.2f %.2f %.2f、ξ1 = %.3f %.3f' % (r.adts['xx'], r.adts['xy'], r.adts['yy'], r.chrom(0), r.chrom(1)))
    for (jx, jy, d) in ((1e-6, 1e-6, 0.0), (2e-6, 5e-7, 0.0), (1e-6, 1e-6, 1e-3), (1e-6, 1e-6, -1e-3)):
        t = A.tracked_tunes(lat, r, jx, jy, delta=d)
        P('  (Jx, Jy, δ) = (%g, %g, %g): |追跡 − 級数| = %.1e (νx)、%.1e (νy)' % (jx, jy, d, abs(t[0] - r.nu(0, jx, jy, d)), abs(t[1] - r.nu(1, jx, jy, d))))
    tr = A.tracked_adts(lat, r)
    P('  J¹ (追跡 + NAFF のフィット) xx xy yy = %.2f %.2f %.2f、ampnf との相対差 %.1e %.1e %.1e' % (
        tr['xx'], tr['xy'], tr['yy'], *(abs(r.adts[k] - tr[k]) / abs(tr[k]) for k in ('xx', 'xy', 'yy'))))
    res3 = A.one_turn_series(lat, D=3, nslice=16, bend_cubic=True, nslice_bend=8)
    h = A.rdt_series(res3)
    h.pop('_imag')
    rm = T.ring_map(lat.elements, 3e9, res3['Amat'], amp=1e-5)
    ht = T.rdt_from_map(rm)
    tt = T.tracked_rdt_table(lat.elements, 3e9, res3['Amat'])
    ht.update(tt['delta1'])
    ht['h10002'] = tt['h10002']
    rel = {k: abs(h[k] - ht[k]) / max(abs(ht[k]), 1e-3) for k in h if k in ht}
    worst = sorted(rel, key=rel.get, reverse=True)[:3]
    P('  h (ampnf 対 追跡の 2 次のテイラー展開) の相対差の最大: ' + ' '.join('%s %.1e' % (k, rel[k]) for k in worst))


if __name__ == '__main__':
    compat.set_mode('k64')
    for k in (sys.argv[1:] or list(KINDS)):
        run(k)
