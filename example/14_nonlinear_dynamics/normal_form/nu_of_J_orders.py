# -*- coding: utf-8 -*-
"""ν(J) の級数 (ampnf) の J の次数 (J¹ = ADTS のみ・J²・J³) ごとの誤差を、追跡の ν と比べる。
追跡: pysad (k64、RF・放射なし、1024 周、NAFF Hann p=1)。基準: SAD (gs) の TrackParticles (同じ初期座標、NORFSW、数点)。
初期座標は正規形の作用 (Jx, Jy) と位相から作る (`ampnf.initial_from_nf`。Twiss の J から出発すると O(√J) ずれるので、級数の J と対応が付かない)。
点は (Jx, Jy) = (J, J/2) の J を 2 倍刻みで。級数は bend_cubic・kin・fringe つきの模型 (docs/25 §13)。
使い方: py -3.12 example/14_nonlinear_dynamics/normal_form/nu_of_J_orders.py fodo|astrid|kekls [--gs] [--D 7]"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'adts_rdt'))
os.environ.setdefault('PYSAD_COMPAT', 'k64')
from pysad import ampnf as A, compat, naff
from pysad.parser import load_sad_string

compat.set_mode('k64')
import adts_gs as Q

CFG = {   # key: (D, 6 極の分割、偏向磁石の分割、四極の分割、J 列)
    'fodo': (7, 8, 8, 4, [5e-7 * 2 ** k for k in range(6)]),
    'astrid': (7, 8, 8, 4, [1e-8 * 2 ** k for k in range(7)]),
    'kekls': (5, 4, 4, 2, [2e-8 * 2 ** k for k in range(6)]),
}
PH = (0.3, 0.9)


def P(*a):
    print(*a, flush=True)


def nu_from_rec(r, rec4, orbit=None):
    """rec4 (nturn, 4) 物理座標 → (νx, νy) (NAFF。r.nu0 に近い側の端数。`ampnf.tracked_tunes` と同じ)"""
    zr = (r.res['Ainv'] @ (rec4 - r.res.get('orbit', np.zeros(4))).T).T
    out = []
    for p in range(2):
        xh, ph = zr[:, 2 * p], zr[:, 2 * p + 1]
        zz = (xh - xh.mean()) - 1j * (ph - ph.mean())
        f0 = (r.nu0[p] + 0.5) % 1.0 - 0.5
        f, _ = naff.naff(zz, 1, window=1, fmin=f0 - 0.05, fmax=f0 + 0.05)
        out.append(f[0] % 1.0)
    return np.array(out)


def main(key, use_gs=False, Dord=None, nturn=1024):
    D, ns, nb, nq, Js = CFG[key]
    D = Dord or D
    src, line, p0, _ = Q.read_ring(key)
    lat = load_sad_string(src).build_lattice(line)
    t0 = time.time()
    r = A.amplitude_tune(lat, D=D, nslice=ns, bend_cubic=True, nslice_bend=nb, kin=True, fringe=True, nslice_quad=nq)
    P('# %s: D=%d  ν0 = (%.6f, %.6f)  最小の共鳴の分母 %.3f  (ampnf %.0f s)' % (key, D, r.nu0[0], r.nu0[1], r.min_den, time.time() - t0))
    P('# J¹: xx xy yy = %.4g %.4g %.4g  J²: x20 x11 x02 = %.3g %.3g %.3g | y20 y11 y02 = %.3g %.3g %.3g' % (
        r.adts['xx'], r.adts['xy'], r.adts['yy'], r.adts2['x20'], r.adts2['x11'], r.adts2['x02'], r.adts2['y20'], r.adts2['y11'], r.adts2['y02']))
    P('# J³ の係数 (x30 x21 x12 x03): %s' % ['%.3g' % r.coef(0, a, 3 - a) for a in (3, 2, 1, 0)] if D >= 7 else '# (D < 7: J³ なし)')
    P('# 点: (Jx, Jy) = (J, J/2)、位相 %s。Δν = ν(追跡) − ν(級数 jorder)' % (PH,))
    jmaxo = 3 if D >= 7 else 2
    hdr = 'J(m rad)  plane  ν_track        Δν(J¹)      Δν(J²)' + ('      Δν(J³)' if jmaxo == 3 else '') + '   |J² の項|'
    P(hdr)
    recs, pts = {}, []
    nus = []
    for J in Js:
        jx, jy = J, J / 2
        t1 = time.time()
        nu = np.array(A.tracked_tunes(lat, r, jx, jy, phase=PH, nturn=nturn, momentum_eV=p0))
        nus.append(nu)
        pts.append((jx, jy))
        for pl in range(2):
            ser = [r.nu(pl, jx, jy, 0.0, jorder=o) for o in range(1, jmaxo + 1)]
            d = [((nu[pl] - s + 0.5) % 1.0) - 0.5 for s in ser]
            j2 = r.nu(pl, jx, jy, 0.0, jorder=2) - r.nu(pl, jx, jy, 0.0, jorder=1)
            P('%.2e  %s  %.8f  %+.3e  %+.3e' % (J, 'xy'[pl], nu[pl], d[0], d[1]) + ('  %+.3e' % d[2] if jmaxo == 3 else '') + '  %.3e' % abs(j2))
    if use_gs:
        sel = [0, len(Js) // 2, len(Js) - 1]
        x0 = np.array([A.initial_from_nf(r, pts[i][0], pts[i][1], PH[0], PH[1], 0.0) for i in sel])
        x6 = np.zeros((len(sel), 6))
        x6[:, :4] = x0
        sext = {e.name: float(e._raw['K2']) for e in lat.elements if type(e).__name__ == 'Sextupole'}
        octu = {e.name: float(e._raw['K3']) for e in lat.elements if type(e).__name__ == 'Octupole'}
        scr = Q.make_script(src, line, p0, 1.0, 1.0, x6, nturn, sext, octu, by_text=True)
        t2 = time.time()
        out, f = Q.run_gs(scr, timeout=3600)
        rec = Q.parse(out, len(sel), nturn)
        P('# gs (TrackParticles、NORFSW、%d 周、%d 点) %.0f s。gs の残り: %s' % (nturn, len(sel), time.time() - t2, Q.check_gs_left() or 'なし'))
        if np.isnan(rec).any():
            P('gs 出力が不完全', out[-300:])
            return
        P('J        plane  ν_gs            ν_pysad_track    差 (gs − pysad)')
        for q, i in enumerate(sel):
            nug = nu_from_rec(r, rec[:, q, :])
            for pl in range(2):
                P('%.2e  %s  %.10f  %.10f  %+.2e' % (Js[i], 'xy'[pl], nug[pl], nus[i][pl], nug[pl] - nus[i][pl]))


if __name__ == '__main__':
    key = sys.argv[1] if len(sys.argv) > 1 else 'fodo'
    Dord = int(sys.argv[sys.argv.index('--D') + 1]) if '--D' in sys.argv else None
    main(key, use_gs='--gs' in sys.argv, Dord=Dord)
