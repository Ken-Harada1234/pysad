# -*- coding: utf-8 -*-
"""追跡の判定の関数 (sextupole_tracking_check.py と sextupole_settings_table.py が使う)。設定の説明は sextupole_tracking_check.py の docstring。"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sxt_common as C

C.setup_pysad()


def track_survive(st, loc, coords, nturn):
    r = st.track(loc, coords, nturn)
    return r['lost_turn'] == 0


def da_rays(st, loc, nturn, quick):
    row = st.rows[loc - 1]
    bx, by = row[1], row[3]
    th = np.deg2rad(np.linspace(0.0, 90.0, 9 if not quick else 3))
    nc = 24 if not quick else 8
    umax = 4.0e-3
    ugrid = np.linspace(umax / nc, umax, nc)

    def pts(u, t):
        return u * np.cos(t) * math.sqrt(bx), u * np.sin(t) * math.sqrt(by) + 1e-5

    cod = st.refine_orbit(np.zeros(6))              # 誤差なしなら 0
    coords = []
    for t in th:
        for u in ugrid:
            x, y = pts(u, t)
            coords.append([x + cod[0], cod[1], y + cod[2], cod[3], 0.0, 0.0])
    alive = track_survive(st, loc, np.array(coords), nturn).reshape(len(th), nc)
    res = []
    for i, t in enumerate(th):
        lost_idx = np.where(~alive[i])[0]
        j = lost_idx[0] if len(lost_idx) else nc          # 最初に失われた点
        lo = ugrid[j - 1] if j > 0 else 0.0
        hi = ugrid[j] if j < nc else ugrid[-1]
        for _ in range(2 if j < nc else 0):                # 二分 2 回 (各 6 点)
            sub = np.linspace(lo, hi, 8)[1:-1]
            cc = np.array([[pts(u, t)[0] + cod[0], cod[1], pts(u, t)[1] + cod[2], cod[3], 0.0, 0.0] for u in sub])
            a = track_survive(st, loc, cc, nturn)
            bad = np.where(~a)[0]
            k = bad[0] if len(bad) else len(sub)
            new_lo = sub[k - 1] if k > 0 else lo
            new_hi = sub[k] if k < len(sub) else hi
            lo, hi = new_lo, new_hi
        res.append(lo)
    res = np.array(res)
    xs = res * np.cos(th) * math.sqrt(bx)
    ys = res * np.sin(th) * math.sqrt(by)
    poly_x = np.concatenate([[0.0], xs, [0.0]])
    poly_y = np.concatenate([[0.0], ys, [0.0]])
    area_q = 0.5 * abs(np.sum(poly_x[:-1] * poly_y[1:] - poly_x[1:] * poly_y[:-1]))
    return dict(theta_deg=np.rad2deg(th).tolist(), u=res.tolist(), x_mm=(xs * 1e3).tolist(), y_mm=(ys * 1e3).tolist(),
                capped=int(np.sum(res >= ugrid[-1] - 1e-12)), area_mm2=float(4 * area_q * 1e6), x_da_mm=float(xs[0] * 1e3), y_da_mm=float(ys[-1] * 1e3), beta=[bx, by])


def da_xdp(st, loc, nturn, quick):
    dps = np.array([-0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03]) if not quick else np.array([-0.01, 0.0, 0.01])
    xgrid = np.linspace(1e-3, 24e-3, 24 if not quick else 8)
    out = []
    for dp in dps:
        try:
            cod = st.closed_orbit_dp(loc, dp)
        except Exception:
            out.append([0.0, 0.0])
            continue
        w = []
        for sg in (+1, -1):
            cc = np.zeros((len(xgrid), 6))
            cc[:] = cod
            cc[:, 0] += sg * xgrid
            cc[:, 2] += 1e-4
            al = track_survive(st, loc, cc, nturn)
            bad = np.where(~al)[0]
            j = bad[0] if len(bad) else len(xgrid)
            lo = xgrid[j - 1] if j > 0 else 0.0
            hi = xgrid[j] if j < len(xgrid) else xgrid[-1]
            if j < len(xgrid):
                sub = np.linspace(lo, hi, 8)[1:-1]
                cc2 = np.zeros((len(sub), 6))
                cc2[:] = cod
                cc2[:, 0] += sg * sub
                cc2[:, 2] += 1e-4
                a2 = track_survive(st, loc, cc2, nturn)
                b2 = np.where(~a2)[0]
                k = b2[0] if len(b2) else len(sub)
                lo = sub[k - 1] if k > 0 else lo
            w.append(lo)
        out.append(w)
    out = np.array(out)
    width = (out[:, 0] + out[:, 1]) * 1e3
    return dict(dp=dps.tolist(), x_plus_mm=(out[:, 0] * 1e3).tolist(), x_minus_mm=(out[:, 1] * 1e3).tolist(), width_mm=width.tolist(),
                mean_width_mm=float(np.mean(width)), width0_mm=float(width[len(dps) // 2]))


def fmap(st, loc, quick):
    from pysad.freqmap import frequency_map
    x = np.linspace(0.3e-3, 7e-3, 14 if not quick else 5)
    y = np.linspace(0.3e-3, 5e-3, 8 if not quick else 3)
    r = frequency_map(st, loc, x=x, y=y, nturn=1024 if not quick else 256, workers=1, warn=False)
    lost = r.lost_turn > 0
    D = r.D
    ok = np.isfinite(D)
    return dict(nx=len(x), ny=len(y), frac_lost=float(np.mean(lost)), frac_chaotic=float(np.mean(lost | (ok & (D > -3.0)) | ~np.isfinite(D))),
                meanD=float(np.nanmean(D[ok])) if ok.any() else float('nan'), D=[None if not np.isfinite(v) else float(v) for v in D], lost=r.lost_turn.tolist(),
                x=r.x.tolist(), y=r.y.tolist())


def lma(src, quick, ring='RING'):
    from pysad.lma import momentum_acceptance
    from pysad.trksetup import select_locations
    st = C.make_track_setup(src, rad=True, rfsw=True, alost=0.02, ring=ring)
    locs_all = select_locations(st, ds=190.0)
    locs = locs_all[:4] if not quick else locs_all[:2]
    res = momentum_acceptance(st, locs=locs, workers=1, tol=2e-4, nsub=3, ncoarse=12, quiet=True, progress=0, turncheck=0)
    dpp, dpm = np.array(res.dp_plus), -np.array(res.dp_minus)
    return dict(loc=list(map(int, res.loc)), s=res.s.tolist(), dp_plus=dpp.tolist(), dp_minus=dpm.tolist(), reason_plus=res.reason_plus, reason_minus=res.reason_minus,
                mean_min=float(np.mean(np.minimum(dpp, dpm))), mean_plus=float(np.mean(dpp)), mean_minus=float(np.mean(dpm)),
                bucket=res.cond.get('bucket_height'), nturn=res.cond.get('nturn'))


