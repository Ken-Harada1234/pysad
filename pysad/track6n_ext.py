"""
track6n_ext.py — track6n.py の拡張: MULT (tmulti) とソレノイド領域 (tsol) の numba カーネル

track6.py の _tmulti_ens / _tsol_region (SAD tmulti.f / tsol.f / tdrift.f / tsolqu.f の直訳) を、
track6n.py と同じ流儀 (S(n,10) 配列、@njit(error_model='numpy')、fmax/fmin、RNG は int64 配列 rs) で
粒子配列に対するカーネルにしたもの。演算の順序は track6.py と同じ (結果は 1e-13 以下で一致、FLUC の乱数消費も同じ)。

対応:
  MULT  k_mult: K0..K_n (複素 akr)、DX/DY/ROTATE(+DROTATE+akang)、ndiv 分割 (tsolqum (放射なし) / tsolqu+tradk (放射あり))、
        K0 付きは tdrift (tsolconv) の厳密ドリフト、フリンジ ttfrins (dofr)、L=0 の薄い多極。
        ANGLE/F1/FRINGE/VOLT 付きの MULT は未対応 (NotImplementedError、track6 も同じ)。
  SOL   k_sol_frin (tsfrin)、k_sol_drift (tdrift_solenoid、放射は tradk1 相当)、
        k_sol_quad (tquad bz≠0: 回転 → ttfrin4(bz) → tsolqu/tsolqur の bz≠0 分岐 → 戻し)。
        tsolqur (RAD ON, bz≠0) の最後の半ステップは「最後の粒子の bzp を全粒子に使う」SAD の癖を再現。
        SOL の DX/DY/DZ/CHI1-3/F1、BOUND=0 の SOL、領域内の BEND/MULT/CAVI は未対応 (track6 と同じ)。

track6n.py は別担当が編集中なので変更せず、ここでは
  ProgExt (track6n.Prog の派生: MULT とソレノイド領域を compile)、turn_ext (track6n.turn + 新しい kind の分岐)、
  track_particles_ext (track6n.track_particles_fast と同じ仕様、フォールバックなし)
を置く。track6n.py への組み込み案は docs/history/work_track6n_ext.md 参照。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import emit6d as _e
from . import lattice as _lat
from . import track6 as _t6
from . import track6n as _n
from .track6n import (NJ, fmax, fmin, sqrt1, pxy2dpz, asinz, xsincos, xsincosh, drift1, rot1, setr0_1,
                      tradk_core1, tradk_all, tzsetparam0, tsolqu1, tsolqur_all, ttfrin_n_1, ndivrad,
                      k_apert, K_DRIFT, K_NONE, SM_NONE, SM_START, SM_IN, SM_END)
from .elements import Drift, Quadrupole

# 追加する要素の種類 (track6n の K_* = 0..6 と重ならないように 10 から)
K_MULT, K_SOL_FRIN, K_SOL_DRIFT, K_SOL_QUAD = 10, 11, 12, 13

# ソレノイド領域の la (tapert) の扱い (SM_NONE 等) は track6n.py で定義 (Prog.sol_mode)


# ============================================================================ tdrift (tsolconv 付きの双極キック)
@njit(**NJ)
def tsolconv(pz0, plz, pbz, bpr, phi):
    """tdrift.f tsolconv: phi の Newton 反復。戻りの sin 等は最後に評価した phi (更新前) の値 (SAD のまま)"""
    sn, xs, cs, dc = xsincos(phi)
    if plz == 0.0:
        return phi, sn, xs, cs, dc
    for _ in range(15):
        sn, xs, cs, dc = xsincos(phi)
        sv = plz * xs + pz0 * sn - pbz * dc
        u = -plz * dc + pz0 * cs + pbz * sn
        if u != 0.0:
            dphi = (bpr - sv) / u
        else:
            dphi = (bpr - sv) / pz0
        phi0 = phi
        phi = phi + dphi
        if phi0 == phi or abs(dphi) <= 1e-15 * abs(phi):
            break
    return phi, sn, xs, cs, dc


@njit(**NJ)
def tdrift_kick1(i, S, al, ak0x, ak0y):
    """tdrift.f tdrift (bz=0, ak0x/ak0y≠0, 放射なし): 一定の双極キック付きの厳密ドリフト (track6 _tdrift_kick)"""
    if ak0x == 0.0 and ak0y == 0.0:
        drift1(i, S, al)
        return
    b = math.hypot(ak0x, ak0y)
    phix = ak0y / b
    phiy = ak0x / b
    phiz = 0.0
    dphizsq = phix * phix + phiy * phiy
    pr = 1.0 + S[i, 5]
    alb = al * pr / b
    px0 = S[i, 1]; py0 = S[i, 3]
    dpz0 = pxy2dpz(px0, py0)
    pz0 = 1.0 + dpz0
    dpl = px0 * phix + py0 * phiy + dpz0 * phiz
    pl = phiz + dpl
    plx = pl * phix; ply = pl * phiy; plz = pl * phiz
    ptx = px0 - plx
    pty = py0 - ply
    ptz = dpz0 - dpl * phiz + dphizsq
    pbx = pty * phiz - ptz * phiy
    pby = ptz * phix - ptx * phiz
    pbz = ptx * phiy - pty * phix
    bpr = b / pr
    db = bpr - pbz
    ap = math.hypot(pz0, pbz)
    dphi = -math.atan2(pbz, pz0)
    phi = asinz(db / ap) - dphi
    phi, sinphi, xsinphi, cosphi, dcosphi = tsolconv(pz0, plz, pbz, bpr, phi)
    S[i, 0] += (plx * phi + ptx * sinphi - pbx * dcosphi) * alb
    S[i, 2] += (ply * phi + pty * sinphi - pby * dcosphi) * alb
    S[i, 4] += ((dpl * phiz - dphizsq) * xsinphi + dpz0 * sinphi - pbz * dcosphi) * alb - S[i, 6] * al
    S[i, 1] = px0 + ptx * dcosphi + pbx * sinphi
    S[i, 3] = py0 + pty * dcosphi + pby * sinphi


# ============================================================================ tsolqu / tsolqum (MULT 用: K0 付き)
@njit(**NJ)
def tsolqu_m1(i, S, al, ak, eps0, ak0x, ak0y, qum):
    """tsolqu.f tsolqu / tsolqum (bz=0, 放射なし) の MULT 用 (track6 _tsolqu と同じ)。
    ak = |K1|L ≥ 0、ak0x/ak0y は分割あたりの K0 (dx0 = ak0x/ak)。ak = 0 なら tdrift_kick1。
    qum=True は tsolqum (放射なしの tmulti): z の dv の項が −0.25·(ap − dv)·aln0 (SAD のまま)"""
    if ak == 0.0:
        tdrift_kick1(i, S, al, ak0x, ak0y)
        return
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    ndiv = 1 + int(abs(al * ak / eps))
    aln0 = al / ndiv
    dx0 = ak0x / ak
    dy0 = ak0y / ak
    akk = ak / al
    w1, s1, xs1, dc1, sh2, xsh2, dch2 = tzsetparam0(S[i, 5], aln0, akk)
    x = S[i, 0]; px = S[i, 1]; y = S[i, 2]; py = S[i, 3]; z = S[i, 4]
    dv = S[i, 6]
    for n in range(ndiv):
        pxi = px; pyi = py
        ap = pxi * pxi + pyi * pyi
        dpz = sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * aln0
        if n == 0:
            r *= 0.5
        x += pxi * r
        y += pyi * r
        xi = x + dx0; yi = y + dy0
        u1 = xi * dc1 + pxi * s1 / w1
        u2 = -xi * w1 * s1 + pxi * dc1
        v1 = yi * dch2 + pyi * sh2 / w1
        v2 = yi * w1 * sh2 + pyi * dch2
        x += u1
        px = pxi + u2
        y += v1
        py = pyi + v2
        if qum:
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi * xi * xs1 - yi * yi * xsh2) + (ap - dv) * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2))
        else:
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi * xi * xs1 - yi * yi * xsh2) + ap * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2)
                 - dv * aln0)
    ap = px * px + py * py
    dpz = sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0 * 0.5
    x += px * r
    y += py * r
    z -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
    S[i, 0] = x; S[i, 1] = px; S[i, 2] = y; S[i, 3] = py; S[i, 4] = z


# ============================================================================ フリンジ (ttfrins: 複素 ak)
@njit(**NJ)
def ttfrins1(i, S, nord, akre, akim, al):
    """tquad.f ttfrins: 複素 ak (回した座標系の係数) のフリンジ。主軸へ回して ttfrin、戻す"""
    if (akre == 0.0 and akim == 0.0) or al == 0.0:
        return
    if akre == 0.0:
        theta = math.pi / nord
        aka = akim
    elif akim == 0.0:
        theta = 0.0
        aka = akre
    else:
        theta = math.atan2(akim, akre) * 2.0 / nord
        aka = math.hypot(akre, akim)
    if theta != 0.0:
        ct = math.cos(theta)
        st = math.sin(theta)
        rot1(i, S, ct, st)
        ttfrin_n_1(i, S, nord, aka, al)
        rot1(i, S, ct, -st)
    else:
        ttfrin_n_1(i, S, nord, aka, al)


# ============================================================================ 多極キック
@njit(**NJ)
def mult_kick1(i, S, akrm, nmmax, nmin):
    """tmulti の多極キック (Horner)。nmin=2: 長さあり (n≥2 の項、最後に ·cx1/(1+δ))、
    nmin=1: 薄い (L=0) 多極 (n≥1 の項 + akrm[0]、/(1+δ))"""
    cx1 = complex(S[i, 0], S[i, 2])
    cx = akrm[nmmax] * cx1 / float(nmmax)
    for n in range(nmmax - 1, nmin - 1, -1):
        cx = (cx + akrm[n]) * cx1 / float(n)
    if nmin == 2:
        cx = cx * cx1 / (1.0 + S[i, 5])
    else:
        cx = (cx + akrm[0]) / (1.0 + S[i, 5])
    S[i, 1] -= cx.real
    S[i, 3] += cx.imag


# ============================================================================ MULT (tmulti)
@njit(**NJ)
def k_mult(S, n, al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad,
           krad, fluc, rs, p0, h0, anrad, cuc):
    """tmulti.f tmulti (ANGLE=0, BZ=0, F1/RF なし) = track6 _tmulti_ens。
    akr: 複素 K_n·cr1^(n+1) (n=0..nmmax)、dofr: フリンジを掛ける次数、ndiv0: tpara の ndiv、eps_rad: EPS (ndivrad 用)"""
    for i in range(n):                   # tsolrot (ent, bz=0): ずらす → 回す
        S[i, 0] -= dx
        S[i, 2] -= dy
        rot1(i, S, ct, st)
    akr1 = akr[1].real if nmmax >= 1 else 0.0
    akr0 = akr[0]
    ndiv = ndiv0
    nz = al != 0.0
    if nz:
        if krad:
            for i in range(n):
                setr0_1(i, S)
            ndiv = max(ndiv, ndivrad(abs(akr0), akr1, 0.0, eps_rad, h0))
        ndiv = min(1000, ndiv)
        if fringe:
            for i in range(n):
                for m in range(nmmax + 1):
                    if dofr[m]:
                        ttfrins1(i, S, m * 2 + 2, akr[m].real, akr[m].imag, al)
    wi = 1.0 / ndiv
    ak1 = akr1 * wi * 0.5
    al1 = al * wi * 0.5
    ak01 = akr0 * wi * 0.5
    akrm = np.empty(nmmax + 1, dtype=np.complex128)
    for m in range(nmmax + 1):
        akrm[m] = akr[m] * wi
    for m in range(1, ndiv + 1):
        if nz:
            for i in range(n):
                tsolqu_m1(i, S, al1, ak1, 0.0, ak01.real, ak01.imag, not krad)   # 放射なしは tsolqum
            if krad:
                tradk_all(S, n, al1, 0.0, fluc, rs, p0, anrad, cuc)
            wm = wi * 0.5 if m == ndiv else wi
            al1 = al * wm
            ak1 = akr1 * wm
            ak01 = akr0 * wm
            if nmmax >= 2:
                for i in range(n):
                    mult_kick1(i, S, akrm, nmmax, 2)
        elif nmmax >= 1:
            for i in range(n):
                mult_kick1(i, S, akrm, nmmax, 1)
    if nz:
        for i in range(n):
            tsolqu_m1(i, S, al1, ak1, 0.0, ak01.real, ak01.imag, not krad)
        if fringe:
            for i in range(n):
                for m in range(nmmax, -1, -1):
                    if dofr[m]:
                        ttfrins1(i, S, m * 2 + 2, -akr[m].real, -akr[m].imag, al)
        if krad:
            tradk_all(S, n, al1, 0.0, fluc, rs, p0, anrad, cuc)
    for i in range(n):                   # tsolrot (exit)
        rot1(i, S, ct, -st)
        S[i, 0] += dx
        S[i, 2] += dy


# ============================================================================ ソレノイド: tsfrin, tdrift_solenoid
@njit(**NJ)
def tsfrin1(i, S, dbz):
    """tsol.f tsfrin: ソレノイドの端のフリンジ (bz の変化 dbz)"""
    x0 = S[i, 0]; px0 = S[i, 1]; y0 = S[i, 2]; py0 = S[i, 3]
    bq = -0.125 * dbz / (1.0 + S[i, 5])
    bp = -2.0 * bq
    pr = x0 * px0 + y0 * py0
    pphi = x0 * py0 - y0 * px0
    phi0 = bq * pr
    w = math.exp(bq * pphi)
    sn, xs, cs, dc = xsincos(phi0)
    S[i, 0] = (cs * x0 - sn * y0) * w
    S[i, 2] = (cs * y0 + sn * x0) * w
    S[i, 1] = (cs * px0 - sn * py0) / w
    S[i, 3] = (cs * py0 + sn * px0) / w
    S[i, 4] = S[i, 4] + bp * pr * pphi


@njit(**NJ)
def k_sol_frin(S, n, dbz):
    if dbz == 0.0:
        return
    for i in range(n):
        tsfrin1(i, S, dbz)


@njit(**NJ)
def tdrift_solenoid1(i, S, al, bz, enarad, fluc, rs, p0, anrad, cuc):
    """tdrift.f tdrift_solenoid (放射は tradk1: 基準 (pxi, pyi, zi)、cphi0=1)"""
    pr = 1.0 + S[i, 5]
    bzp = bz / pr
    pxi = S[i, 1] + bzp * S[i, 2] * 0.5
    pyi = S[i, 3] - bzp * S[i, 0] * 0.5
    zi = S[i, 4]
    dpzi = pxy2dpz(pxi, pyi)
    pzi = 1.0 + dpzi
    al1 = al / pzi
    phi = al1 * bzp
    a24, a12, a22, a14 = xsincos(phi)
    S[i, 0] += (a24 * pxi - a14 * pyi) / bzp
    S[i, 2] += (a14 * pxi + a24 * pyi) / bzp
    S[i, 4] += dpzi * al1 - S[i, 6] * al
    px1 = a22 * pxi + a24 * pyi
    py1 = -a24 * pxi + a22 * pyi
    if enarad:
        S[i, 1] = px1; S[i, 3] = py1
        S[i, 7] = pxi; S[i, 8] = pyi; S[i, 9] = zi
        tradk_core1(i, S, al, 0.0, fluc, rs, p0, anrad, cuc)
        px1 = S[i, 1]; py1 = S[i, 3]
    pr = 1.0 + S[i, 5]
    bzp = bz / pr
    S[i, 1] = px1 - bzp * S[i, 2] * 0.5
    S[i, 3] = py1 + bzp * S[i, 0] * 0.5


@njit(**NJ)
def k_sol_drift(S, n, al, bz, krad, fluc, rs, p0, anrad, cuc):
    """ソレノイド領域内の DRIFT: bz≠0 は tdrift_solenoid、bz=0 は tdrift_free"""
    if bz == 0.0:
        for i in range(n):
            drift1(i, S, al)
        return
    for i in range(n):
        tdrift_solenoid1(i, S, al, bz, krad, fluc, rs, p0, anrad, cuc)


# ============================================================================ ソレノイド内の四極
@njit(**NJ)
def ttfrin4_bz1(i, S, ak, al, bz):
    """tquad.f ttfrin (nord=4, bz≠0)"""
    if al == 0.0 or ak == 0.0:
        return
    pr = 1.0 + S[i, 5]
    aki = ak / al / 4.0 / pr
    bzph = 0.5 * bz / pr
    xi = S[i, 0]; yi = S[i, 2]
    a = aki * (xi * xi)
    b = aki * (yi * yi)
    ab = a - b
    t = (ab * ab) / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    S[i, 0] = xi + dx1
    S[i, 2] = yi + dy1
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h * h
    if bz == 0.0:
        px1 = (S[i, 1] * yy + S[i, 3] * h) / f
        py1 = (S[i, 3] * xx - S[i, 1] * h) / f
        S[i, 4] = S[i, 4] - px1 * (dx1 + t * xi) - py1 * (dy1 + t * yi)
        S[i, 1] = px1
        S[i, 3] = py1
        return
    pxa = S[i, 1] + bzph * S[i, 2]
    pya = S[i, 3] - bzph * S[i, 0]
    S[i, 1] = (pxa * yy + pya * h) / f - bzph * yi
    S[i, 3] = (pya * xx - pxa * h) / f + bzph * xi
    S[i, 4] = S[i, 4] - S[i, 1] * (dx1 + t * xi) - S[i, 3] * (dy1 + t * yi)


@njit(**NJ)
def tztafs(zwu, pxt, pyt, zw1, wst, w12t, swss, g1t,
           xsh2, xs1, c1, ch2, csw1, akkp, bzp, csws, ca1, cr3, dcw1, dcw2, cr2):
    """tsolque.f tztafs → dz (track6 _TZ.tztafs)"""
    w12xsh = w12t * xsh2
    wsxs = wst * xs1
    cwsd = c1 * w12xsh - ch2 * wsxs
    cswa = (csw1 + cwsd) * swss
    wsm = 4.0 * swss * akkp
    ws1 = (-1.0 + wsm) * wst
    ws2 = w12t * (1.0 + wsm)
    cr1 = (bzp * bzp) * swss
    f1 = cswa + csws - zw1
    f2 = ca1 * cr1 + g1t
    f3 = (2.0 * cr1 * cwsd + cr3 * w12xsh
          + dcw1 * ws1 + dcw2 * ws2 - cr2 * wsxs)
    return zwu * (zwu * f3 / 2.0 + f2 * pxt + bzp * f1 * pyt) - (zw1 * (pyt * pyt)) / 2.0


@njit(**NJ)
def tz_step1(i, S, aln0, akk, bz, first):
    """tsolqu/tsolqur の bz≠0 分岐の 1 スライス (tzsetparams + 本体)。戻り: この粒子の bzp"""
    # --- tzsetparams (track6 _TZ.__init__)
    pr = 1.0 + S[i, 5]
    bzp = bz / pr
    akkp = akk / pr
    w1 = math.sqrt(((bzp * bzp) + math.hypot(bzp * bzp, 2.0 * akkp)) * 0.5)
    w2 = akkp / w1
    wss = 1.0 / ((w1 * w1) + (w2 * w2))
    ws = w1 + w2
    w12 = w1 - w2
    wd = bzp / ws
    phi1 = aln0 * w1
    phi2 = fmin(100.0, fmax(-100.0, aln0 * w2))
    s1, xs1, c1, dc1 = xsincos(phi1)
    sh2, xsh2, ch2, dch2 = xsincosh(phi2)
    t1 = s1 * w2
    g1 = (t1 * t1) / akkp
    t2 = sh2 * w1
    g2 = -(t2 * t2) / akkp
    wr1 = w1 / w2
    wr2 = w2 / w1
    cxs1 = dc1 * s1 - xs1
    cxs2 = dch2 * sh2 - xsh2
    aw1 = aln0 + (cxs2 * wr1) / ws
    aw2 = aln0 + (cxs1 * wr2) / ws
    cr2 = c1 * wr2
    cr3 = ch2 * wr1
    csw1 = akkp * aln0 * (-dc1 + dch2)
    csws = (ch2 * w1 * phi1 + c1 * w2 * phi2) * wss
    ca1 = -(ch2 * dc1) - dch2 - (s1 * sh2 * w12) / ws
    dxs = w2 * xs1 - w1 * xsh2
    dcw1 = dc1 * phi2
    dcw2 = dch2 * phi1
    # --- 本体
    px = S[i, 1]; py = S[i, 3]
    ap = px * px + py * py
    dpz = sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0
    if first:
        r *= 0.5
    a24, a12, a22, a14 = xsincos(r * bzp)
    pxi = S[i, 1]
    S[i, 0] += (a24 * pxi - a14 * S[i, 3]) / bzp
    S[i, 2] += (a14 * pxi + a24 * S[i, 3]) / bzp
    pxi = a22 * pxi + a24 * S[i, 3]
    pyi = -a24 * pxi + a22 * S[i, 3]          # SAD のまま: 更新後の pxi を使う
    S[i, 4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
    xi = S[i, 0]; yi = S[i, 2]
    a = (w2 * ws * xi - bzp * pyi) * wss
    bw = (ws * pxi - bzp * w2 * yi) * wss
    cc = (w1 * wd * xi + pyi) * wss
    dw = (-wd * pxi + w1 * yi) * wss
    u1w = a * dc1 + bw * s1
    u2w = -a * s1 + bw * dc1
    v1w = cc * dch2 + dw * sh2
    v2w = cc * sh2 + dw * dch2
    S[i, 0] = S[i, 0] + u1w + v1w * bzp
    S[i, 1] = pxi + u2w * w1 + v2w * w2 * bzp
    S[i, 2] = S[i, 2] + wd * u2w + ws * v2w
    S[i, 3] = pyi - wd * u1w * w1 + ws * v1w * w2
    awu = a / ws * w1
    dwu = dw * w2
    dz1 = tztafs(awu, pxi, pyi, aw1, ws, w12, wss, g1,
                 xsh2, xs1, c1, ch2, csw1, akkp, bzp, csws, ca1, cr3, dcw1, dcw2, cr2)
    dz2 = tztafs(-dwu, -pyi, pxi, aw2, -w12, ws, -wss, g2,
                 xsh2, xs1, c1, ch2, csw1, akkp, bzp, csws, ca1, cr3, dcw1, dcw2, cr2)
    S[i, 4] = (S[i, 4] + bzp * (-((awu * dwu * (dxs * dxs)) / akkp) + ca1 * pxi * pyi * wss)
               + dz1 + dz2 - aln0 * S[i, 6])
    return bzp


@njit(**NJ)
def tsolqu_bz_final1(i, S, aln0, bzp):
    """tsolqu/tsolqur の bz≠0 分岐の最後の半ステップ (bzp は呼び出し側が決める: 癖の再現はそこで)"""
    px = S[i, 1]; py = S[i, 3]
    ap = px * px + py * py
    dpz = sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0 * 0.5
    a24, a12, a22, a14 = xsincos(r * bzp)
    pxi = S[i, 1]
    S[i, 0] += (a24 * pxi - a14 * S[i, 3]) / bzp
    S[i, 2] += (a14 * pxi + a24 * S[i, 3]) / bzp
    S[i, 1] = a22 * pxi + a24 * S[i, 3] - bzp * S[i, 2] * 0.5
    S[i, 3] = -a24 * pxi + a22 * S[i, 3] + bzp * S[i, 0] * 0.5
    S[i, 4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r


@njit(**NJ)
def tsolqu_bz_all(S, n, al, ak, bz, eps0, krad, fluc, rs, p0, h0, anrad, cuc):
    """tsolqu.f tsolqu / tsolqur の bz≠0 分岐 (ak0x=ak0y=0)。スライス外側・粒子内側、各スライス後に全粒子 tradk。
    最後の半ステップ: tsolqur (krad) では最後に tzsetparams を呼んだ粒子 (最後のスライスの最後の粒子) の bzp を
    全粒子に使う (SAD tsolqu.f:519-531 の癖)。tsolqu (放射なし) は粒子ごとの値"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    aka = math.hypot(ak, bz)
    ndiv = 1 + int(abs(al) * aka / eps)
    if krad:
        ndiv = min(1000, max(ndiv, ndivrad(0.0, ak, bz, eps0, h0)))
    aln0 = al / ndiv
    akk = ak / al
    for i in range(n):
        bzp = bz / (1.0 + S[i, 5])
        S[i, 1] += bzp * S[i, 2] * 0.5
        S[i, 3] -= bzp * S[i, 0] * 0.5
    bzp_last = 0.0
    for m in range(ndiv):
        for i in range(n):
            bzp_last = tz_step1(i, S, aln0, akk, bz, m == 0)
        if krad:
            tradk_all(S, n, aln0, 0.0, fluc, rs, p0, anrad, cuc)
    if krad:
        for i in range(n):
            tsolqu_bz_final1(i, S, aln0, bzp_last)
    else:
        for i in range(n):
            tsolqu_bz_final1(i, S, aln0, bz / (1.0 + S[i, 5]))


@njit(**NJ)
def k_sol_quad(S, n, L, ak, ct, st, fringe, bz, krad, fluc, rs, p0, h0, anrad, cuc):
    """tsol.f の icQUAD → tquad (bz=bzs; DX=DY=F1=FRINGE=0) = track6 _tquad_sol。ak = |K1|L > 0"""
    for i in range(n):
        rot1(i, S, ct, st)               # tsolrot (dx=dy=0): 回転のみ
        if krad:
            setr0_1(i, S)
        if fringe:
            ttfrin4_bz1(i, S, ak, L, bz)
        if krad:
            setr0_1(i, S)
    if bz == 0.0:                        # tsolqu/tsolqur の bz=0 分岐 (EPS は使わない: track6 と同じ)
        if krad:
            tsolqur_all(S, n, L, ak, 0.0, fluc, rs, p0, h0, anrad, cuc)
        else:
            for i in range(n):
                tsolqu1(i, S, L, ak, 0.0, False)
    else:
        tsolqu_bz_all(S, n, L, ak, bz, 0.0, krad, fluc, rs, p0, h0, anrad, cuc)
    for i in range(n):
        if fringe:
            ttfrin4_bz1(i, S, -ak, L, bz)
        rot1(i, S, ct, -st)


# ============================================================================ ラティスの表 (track6n.Prog の拡張)
def compile_mult(el):
    """MULT → (K_MULT, params) / nmmax<0 → (K_DRIFT, (L,))。track6 _tmulti_ens と同じ前処理 (_lat._mult_params)"""
    g = _t6._g
    pm = _lat._mult_params(el)
    al, nmmax, akr, ndiv = pm['al'], pm['nmmax'], pm['akr'], pm['ndiv']
    if g(el, 'ANGLE') or g(el, 'F1') or g(el, 'FRINGE') or g(el, 'VOLT'):
        raise NotImplementedError('track6n_ext: MULT の曲げ・F1・RF は未対応')
    if nmmax < 0:
        return K_DRIFT, (al,)
    eps = pm['eps']
    if al != 0.0:
        dofr = [akr[m] != 0 and (m <= 2 or abs(akr[m]) / al * 0.05**(m + 1) / (m + 1) > eps)
                for m in range(nmmax + 1)]
    else:
        dofr = [False] * (nmmax + 1)     # al=0 ではフリンジを掛けない (track6 では未使用)
    ct, st = math.cos(pm['theta2']), math.sin(pm['theta2'])
    return K_MULT, (al, nmmax, np.array(akr, dtype=np.complex128), np.array(dofr, dtype=np.bool_), ndiv,
                    pm['dx'], pm['dy'], ct, st, pm['fringe'], g(el, 'EPS'),
                    g(el, 'RAD') == 0.0 and al != 0.0)


def compile_sol_quad(el, bzs):
    """ソレノイド内の QUAD → (K_SOL_QUAD, params) (track6 _tquad_sol と同じ制限)"""
    g = _t6._g
    L = el.L
    K1L = g(el, 'K1') * L
    if L == 0.0 or K1L == 0.0:
        raise NotImplementedError('track6n_ext: ソレノイド内の L=0 / K1=0 の四極')
    if any(g(el, k) for k in ('DX', 'DY', 'F1', 'FRINGE')):
        raise NotImplementedError('track6n_ext: ソレノイド内の四極のずれ / F1 / FRINGE')
    ak = abs(K1L)
    th2 = g(el, 'ROTATION') + (0.5 * math.pi if K1L < 0 else 0.0)
    return K_SOL_QUAD, (L, ak, math.cos(th2), math.sin(th2), g(el, 'DISFRIN') == 0.0, bzs,
                        g(el, 'RAD') == 0.0 and L != 0.0)


def compile_sol_region(prog, elements, i, ph):
    """elements[i] (BOUND≠0 の SOL) から終点 (次の BOUND≠0 の SOL) までを compile して prog に足し、終点の添字を返す
    (track6 _tsol_region と同じ制限)。track6n.Prog.__init__ から呼ばれる"""
    g = _t6._g
    sb = elements[i]
    if any(g(sb, k) for k in ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')):
        raise NotImplementedError('track6n_ext: SOL のずれ・回転 (trots)')
    if g(sb, 'F1') != 0.0:
        raise NotImplementedError('track6n_ext: SOL の F1 (tradks)')
    ke = None
    for j in range(i + 1, len(elements)):
        e = elements[j]
        if type(e).__name__ == 'Solenoid' and g(e, 'BOUND') != 0.0:
            ke = j
            break
    if ke is None:
        raise ValueError('track6n_ext: ソレノイドの終端 (SOL BOUND=1) がない')
    bzs = _t6._tfbzs(sb, ph)
    prog._add(K_SOL_FRIN, (bzs if g(sb, 'DISFRIN') == 0.0 else 0.0,), sb.name, True, SM_START)
    for l in range(i + 1, ke + 1):
        el = elements[l]
        name = type(el).__name__
        sm = SM_END if l == ke else SM_IN
        if isinstance(el, Drift):
            prog._add(K_SOL_DRIFT, (el.L, bzs, g(el, 'RAD') == 0.0), el.name, True, sm)
        elif isinstance(el, Quadrupole):
            k, p = compile_sol_quad(el, bzs)
            prog._add(k, p, el.name, True, sm)
        elif name == 'Solenoid':
            bz1 = 0.0 if l == ke else _t6._tfbzs(el, ph)
            if any(g(el, k) for k in ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')) or g(el, 'F1') != 0.0:
                raise NotImplementedError('track6n_ext: SOL のずれ・回転・F1')
            prog._add(K_SOL_FRIN, ((bz1 - bzs) if g(el, 'DISFRIN') == 0.0 else 0.0,), el.name, True, sm)
            bzs = bz1
        elif el.L == 0.0:
            prog._add(K_NONE, (), el.name, True, sm)
        else:
            raise NotImplementedError(f'track6n_ext: ソレノイド内の要素 {name} は未対応 (DRIFT/QUAD のみ)')
    return ke


# track6n に統合済み: ProgExt / turn_ext は track6n.Prog / track6n.turn の別名 (互換のため残す)
ProgExt = _n.Prog


def step_kind(kind, pr, S, n, ph, trf0, rfsw, rad, fluc, rs):
    """1 要素分のカーネル呼び出し (track6n.turn の分岐 + 新しい kind)"""
    p0, h0, anrad, cuc = ph.p0, ph.h0, ph.anrad, ph.cuc
    if kind == K_MULT:
        (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad, radok) = pr
        k_mult(S, n, al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad,
               rad and radok, fluc, rs, p0, h0, anrad, cuc)
    elif kind == K_SOL_FRIN:
        k_sol_frin(S, n, pr[0])
    elif kind == K_SOL_DRIFT:
        (al, bz, radok) = pr
        k_sol_drift(S, n, al, bz, rad and radok, fluc, rs, p0, anrad, cuc)
    elif kind == K_SOL_QUAD:
        (L, ak, ct, st, fringe, bz, radok) = pr
        k_sol_quad(S, n, L, ak, ct, st, fringe, bz, rad and radok, fluc, rs, p0, h0, anrad, cuc)
    elif kind == _n.K_DRIFT:
        _n.k_drift(S, n, pr[0])
    elif kind == _n.K_QUAD:
        (L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0, radok) = pr
        _n.k_quad(S, n, L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0,
                  rad and radok, fluc, rs, p0, h0, anrad, cuc)
    elif kind == _n.K_THIN:
        (L, ak, kord, dx, dy, theta, ct, st, fringe, radok) = pr
        _n.k_thin(S, n, L, ak, kord, dx, dy, theta, ct, st, fringe, rad and radok, fluc, rs, p0, anrad, cuc)
    elif kind == _n.K_BEND:
        (al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
         fb1, fb2, mfring, fringe, eps0, dx, dy, cost, sint, dtheta, dchi2, alg, phig, ak1, theta, radok) = pr
        krad = rad and radok
        if ak1 != 0.0:
            _n.k_bendi(S, n, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak1, dx, dy, theta, dtheta, cost, sint,
                       dchi2, alg, phig, fb1, fb2, mfring, krad, fringe, eps0, True, True,
                       fluc, rs, p0, h0, anrad, cuc)
        else:
            _n.k_bend0(S, n, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
                       fb1, fb2, mfring, fringe, krad, True, True, eps0, dx, dy, cost, sint, dtheta, dchi2, alg, phig,
                       fluc, rs, p0, h0, anrad, cuc)
    elif kind == _n.K_CAV:
        (al, v, w, wi, phic, fringe, mfring) = pr
        phis = w * trf0
        dphis = phis - phic
        if rad:
            offset = math.sin(dphis)
            offset1 = 0.0
        else:
            offset = math.sin(dphis) - math.sin(phis)
            offset1 = math.sin(phis)
        tlim = 1e4 if rfsw else 0.0
        _n.k_cav(S, n, al, v, w, wi, dphis, offset, offset1, tlim, fringe, mfring, p0, h0)
    elif kind == _n.K_STEER0:
        _n.k_steer0(S, n, pr[0], pr[1], pr[2], pr[3], pr[4], pr[5])
    elif kind == _n.K_STEER:
        (al, phib, dx, dy, theta, cost, sint, tanp1, eps0, radok) = pr
        _n.k_steer(S, n, al, phib, dx, dy, theta, cost, sint, tanp1, rad and radok, eps0, fluc, rs, ph.p0, ph.h0, ph.anrad, ph.cuc)


turn_ext = _n.turn


def track_particles_ext(elements, coords, nturn, momentum_eV, trf0=0.0, rfsw=True, rad=False, fluc=False,
                        seed=None, charge=1.0, amass=None, progress=0, fshift=0.0):
    """track6.track_particles / track6n.track_particles_fast と同じ仕様 (MULT/SOL 対応、フォールバックなし)"""
    import time
    circ = sum(el.L for el in elements)
    kw = {} if amass is None else {'amass': amass}
    ph = _e._Phys(momentum_eV, circ, charge=charge, **kw)
    ph.dvfs = fshift if rfsw else 0.0
    prog = ProgExt(elements, ph)
    if seed is not None:
        _t6.RNG.restore(seed)
    rs = np.array([_t6.RNG.seed], dtype=np.int64)
    fl = bool(rad and fluc)
    S = _n.make_state(coords, ph)
    n = S.shape[0]
    t0 = time.time()
    for it in range(nturn):
        turn_ext(prog, S, n, ph, trf0, rfsw, rad, fl, rs)
        if not np.all(np.isfinite(S[:n, :6])):
            raise FloatingPointError(f'track6n_ext: {it + 1} 周目で発散')
        if not rfsw:
            S[:n, 4] = 0.0
        if progress and ((it + 1) % progress == 0 or it + 1 == nturn):
            print(f'    [track6n_ext] turn {it + 1}/{nturn} ({time.time() - t0:.1f} s)', flush=True)
    _t6.RNG.seed = int(rs[0])
    return _n.state_out(S, n)
