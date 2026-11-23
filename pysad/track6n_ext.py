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

import copy
import math

import numpy as np
from numba import njit

from . import compat as _compat
from . import emit6d as _e
from . import lattice as _lat
from . import quadbody as _qb
from . import track6 as _t6
from . import track6n as _n
from .track6n import (NJ, p2h, fmax, fmin, sqrt1, pxy2dpz, asinz, xsincos, xsincosh, drift1, rot1, setr0_1,
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
def tsolqu_m1(i, S, al, ak, eps0, ak0x, ak0y, qum, sch):
    """tsolqu.f tsolqu / tsolqum (bz=0, 放射なし) の MULT 用 (track6 _tsolqu と同じ)。
    ak = |K1|L ≥ 0、ak0x/ak0y は分割あたりの K0 (dx0 = ak0x/ak)。ak = 0 なら tdrift_kick1。
    qum=True は tsolqum (tsolqu と違い z の dv の項が −0.25·(ap − dv)·aln0。SAD のまま)。
    本体の積分は quadbody (sch: 0 = K64、1 = 4 次の合成。rs[4])"""
    if ak == 0.0:
        tdrift_kick1(i, S, al, ak0x, ak0y)
        return
    x, px, y, py, z = _n._qb_nonrad(S[i, 0], S[i, 1], S[i, 2], S[i, 3], S[i, 4], S[i, 5], S[i, 6],
                                    al, ak, eps0, ak0x, ak0y, qum, sch)
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
                # 放射なしは tsolqum (K64)。v2 (B14, rs[1]=1) は tsolqu と同じ z の式 (-dv·aln0) を使う
                tsolqu_m1(i, S, al1, ak1, 0.0, ak01.real, ak01.imag, (not krad) and not (rs.shape[0] > 1 and rs[1] == 1), _qb.sch_of(rs))
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
            tsolqu_m1(i, S, al1, ak1, 0.0, ak01.real, ak01.imag, (not krad) and not (rs.shape[0] > 1 and rs[1] == 1), _qb.sch_of(rs))
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
    if krad and (rs.shape[0] > 1 and rs[1] == 1):
        # v2 (B1): 放射は 1 刻みの前後の向きの変化 (弦) から損失を出すので、らせんの回転角 bz·al が大きいと小さく出る (k_sol_bend と同じ。
        # bz·al/ndiv ≤ 0.02)。K64 の tdrift_solenoid は 1 刻み (k64 モードは今まで通り)
        ndiv = 1 + int(abs(bz * al) / 0.02)
        if ndiv > 1:
            w = al / ndiv
            for m in range(ndiv):
                for i in range(n):
                    tdrift_solenoid1(i, S, w, bz, krad, fluc, rs, p0, anrad, cuc)
            return
    for i in range(n):
        tdrift_solenoid1(i, S, al, bz, krad, fluc, rs, p0, anrad, cuc)


@njit(**NJ)
def tdrift_bz1(i, S, al, bz, ak0x, ak0y, enarad, fluc, rs, p0, anrad, cuc):
    """tdrift.f tdrift (bz と ak0x・ak0y の両方がある一様磁場の中の厳密な運動: らせん)。tsol.f の領域内の BEND が呼ぶ
    (tdrift(al, bzs, phiy, phix, krad))。S の px, py は正準運動量 (ソレノイドのゲージ)、中では力学的運動量にして使う。
    放射は tradk1 (基準は入口の力学的運動量)。出口で正準運動量に戻すときの p は放射の後の値 (K64 は放射の前の値)"""
    if ak0x == 0.0 and ak0y == 0.0:
        if bz == 0.0:
            drift1(i, S, al)
        else:
            tdrift_solenoid1(i, S, al, bz, enarad, fluc, rs, p0, anrad, cuc)
        return
    b = math.sqrt(ak0x * ak0x + ak0y * ak0y + (bz * al)**2)
    phix = ak0y / b
    phiy = ak0x / b
    phiz = bz * al / b
    dphizsq = phix * phix + phiy * phiy
    pr = 1.0 + S[i, 5]
    alb = al * pr / b
    bzp = bz / pr
    px0 = S[i, 1] + bzp * S[i, 2] * 0.5
    py0 = S[i, 3] - bzp * S[i, 0] * 0.5
    z0 = S[i, 4]
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
    px1 = px0 + ptx * dcosphi + pbx * sinphi
    py1 = py0 + pty * dcosphi + pby * sinphi
    if enarad:
        S[i, 1] = px1; S[i, 3] = py1
        S[i, 7] = px0; S[i, 8] = py0; S[i, 9] = z0
        tradk_core1(i, S, al, 0.0, fluc, rs, p0, anrad, cuc)
        px1 = S[i, 1]; py1 = S[i, 3]
    pr = 1.0 + S[i, 5]
    bzp = bz / pr
    S[i, 1] = px1 - bzp * S[i, 2] * 0.5
    S[i, 3] = py1 + bzp * S[i, 0] * 0.5


@njit(**NJ)
def k_sol_bend(S, n, al, bz, ak0x, ak0y, krad, fluc, rs, p0, anrad, cuc):
    """ソレノイド領域内の BEND (直線の座標のまま一様な双極磁場: tsol.f の icBEND → tdrift(al, bzs, phiy, phix))"""
    if not krad:
        for i in range(n):
            tdrift_bz1(i, S, al, bz, ak0x, ak0y, False, fluc, rs, p0, anrad, cuc)
        return
    # 放射: tradk は刻みの前後の向きの変化 (弦) から損失を出すので、磁場の向きが回る (らせん) と 1 刻みでは小さく出る
    # (相対 (bz·al)²/12)。bz·al/ndiv ≤ 0.02 に刻む (K64 の tdrift は 1 刻み)
    ndiv = 1 + int(abs(bz * al) / 0.02)
    w = 1.0 / ndiv
    for m in range(ndiv):
        for i in range(n):
            tdrift_bz1(i, S, al * w, bz, ak0x * w, ak0y * w, True, fluc, rs, p0, anrad, cuc)


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
    if krad and not (rs.shape[0] > 1 and rs[1] == 1):
        for i in range(n):
            tsolqu_bz_final1(i, S, aln0, bzp_last)
    else:
        # 放射なし、または v2 (B08): 粒子ごとの bzp (K64 の tsolqur は最後の粒子の値を全粒子に使う)
        for i in range(n):
            tsolqu_bz_final1(i, S, aln0, bz / (1.0 + S[i, 5]))


@njit(**NJ)
def k_sol_quad(S, n, L, ak, ct, st, fringe, bz, krad, fluc, rs, p0, h0, anrad, cuc, dx, dy):
    """tsol.f の icQUAD → tquad (bz=bzs; F1=FRINGE=0) = track6 _tquad_sol。ak = |K1|L > 0。
    dx, dy (v2next B3): 四極のずれ。tsolrot の (dx, dy) = 座標をずらし (x − dx、y − dy)、ソレノイドの一様磁場のベクトルポテンシャルのゲージが
    変わる分だけ正準運動量もずらす (px + (bz/2)dy、py − (bz/2)dx、bz は運動量で規格化)。出口は逆。回転 (ROTATE) はずらしの後 (出口は前)"""
    for i in range(n):
        if dx != 0.0 or dy != 0.0:
            bzp = bz / (1.0 + S[i, 5])
            S[i, 0] -= dx
            S[i, 2] -= dy
            S[i, 1] += 0.5 * bzp * dy
            S[i, 3] -= 0.5 * bzp * dx
        rot1(i, S, ct, st)               # tsolrot: 回転
        if krad:
            setr0_1(i, S)
        if fringe:
            ttfrin4_bz1(i, S, ak, L, bz)
        if krad:
            setr0_1(i, S)
            if bz != 0.0 and rs.shape[0] > 1 and rs[1] == 1:
                # v2 (C2): 放射の基準の向きは力学的運動量。K64 (tquad → tsolqur) は正準運動量のまま pxr0 に入れるので、
                # 最初のスライスの放射に (bz/2)(y, −x) の見かけの偏向が入る (軸を外れた粒子の損失が増える)
                bzp = bz / (1.0 + S[i, 5])
                S[i, 7] = S[i, 1] + bzp * S[i, 2] * 0.5
                S[i, 8] = S[i, 3] - bzp * S[i, 0] * 0.5
    if bz == 0.0:                        # tsolqu/tsolqur の bz=0 分岐 (EPS は使わない: track6 と同じ)
        if krad:
            tsolqur_all(S, n, L, ak, 0.0, fluc, rs, p0, h0, anrad, cuc)
        else:
            for i in range(n):
                tsolqu1(i, S, L, ak, 0.0, False, _qb.sch_of(rs))
    else:
        tsolqu_bz_all(S, n, L, ak, bz, 0.0, krad, fluc, rs, p0, h0, anrad, cuc)
    for i in range(n):
        if fringe:
            ttfrin4_bz1(i, S, -ak, L, bz)
        rot1(i, S, ct, -st)
        if dx != 0.0 or dy != 0.0:
            bzp = bz / (1.0 + S[i, 5])
            S[i, 0] += dx
            S[i, 2] += dy
            S[i, 1] -= 0.5 * bzp * dy
            S[i, 3] += 0.5 * bzp * dx


# ============================================================================ SOL の端 (v2、C2): ずれ・回転 trots と F1 の放射 tradks
K_SOL_EDGE = 17


@njit(**NJ)
def trots1(i, S, r, dx, dy, dz, ent):
    """tsol.f trots: SOL の端の座標変換。r (3×3) は tpara の R11..R33 (CHI1..3 から)、dx dy dz は SOL の DX DY DZ。
    入口 (ent): 回して、新しい z=0 の面まで直進し、(dx, dy) を足す。出口: (dx, dy) を引いて回し、直進する。
    運動量 (S[:,1], S[:,3] は p で割った値) は回すだけ。SOL の外側 (bz=0) で使うので正準運動量 = 力学的運動量。
    z: 設計の行路長が dz だけ伸びる分は dz·(1 − dv) (z = −v(t − t0) で、基準粒子の時刻が dz/v0 進む)。
    K64 は出口だけ −dv·dz を入れ、入口は dz だけ (入口の DZ≠0 は GEO が出口側の SOL にあるとき)。v2 は両方に入れる"""
    pxi = S[i, 1]
    pyi = S[i, 3]
    pzi = 1.0 + pxy2dpz(pxi, pyi)
    if ent:
        xi = S[i, 0]
        yi = S[i, 2]
    else:
        xi = S[i, 0] - dx
        yi = S[i, 2] - dy
    xf = r[0, 0] * xi + r[0, 1] * yi
    yf = r[1, 0] * xi + r[1, 1] * yi
    zf = r[2, 0] * xi + r[2, 1] * yi
    pxf = r[0, 0] * pxi + r[0, 1] * pyi + r[0, 2] * pzi
    pyf = r[1, 0] * pxi + r[1, 1] * pyi + r[1, 2] * pzi
    pzf = r[2, 0] * pxi + r[2, 1] * pyi + r[2, 2] * pzi
    S[i, 1] = pxf
    S[i, 3] = pyf
    if ent:
        S[i, 0] = xf - pxf / pzf * zf + dx
        S[i, 2] = yf - pyf / pzf * zf + dy
    else:
        S[i, 0] = xf - pxf / pzf * zf
        S[i, 2] = yf - pyf / pzf * zf
    S[i, 4] = S[i, 4] + zf / pzf + dz - S[i, 6] * dz


@njit(**NJ)
def k_sol_edge(S, n, rot, r, dx, dy, dz, dbz, bz0, bz1, f1, krad, fluc, rs, p0, anrad, cuc):
    """SOL の端 (tsol.f): rot = 1 入口 (trots → tsfrin → 放射)、−1 出口 (tsfrin → 放射 → trots)、0 領域内の SOL (回転なし)。
    dbz はフリンジの bz の段差 (DISFRIN なら 0)、bz0・bz1 は端の前後の bz。
    放射 (krad: RAD かつ F1≠0 かつ RAD=0 かつ bz が変わる): 端の前の力学的運動量の向きを基準に、tsfrin の後の向きとの差を
    長さ F1 の間の一様な偏向として tradk に渡す (tradks)"""
    for i in range(n):
        if rot == 1:
            trots1(i, S, r, dx, dy, dz, True)
        if krad:
            bzp = 0.5 * bz0 / (1.0 + S[i, 5])
            S[i, 7] = S[i, 1] + bzp * S[i, 2]
            S[i, 8] = S[i, 3] - bzp * S[i, 0]
            S[i, 9] = S[i, 4]
        if dbz != 0.0:
            tsfrin1(i, S, dbz)
        if krad:
            bzp = 0.5 * bz1 / (1.0 + S[i, 5])
            S[i, 1] += bzp * S[i, 2]
            S[i, 3] -= bzp * S[i, 0]
            tradk_core1(i, S, f1, 0.0, fluc, rs, p0, anrad, cuc)
            bzp = 0.5 * bz1 / (1.0 + S[i, 5])
            S[i, 1] -= bzp * S[i, 2]
            S[i, 3] += bzp * S[i, 0]
        if rot == -1:
            trots1(i, S, r, dx, dy, dz, False)


# ============================================================================ ソレノイド領域内の MULT・CAVI (v2)
K_SOL_BEND, K_SOL_MULT = 15, 16


@njit(**NJ)
def sol_body_all(S, n, al, ak1, bz, ak0x, ak0y, rs):
    """領域内の MULT の本体 (放射なし): K1=0 は一様磁場 (ソレノイド + 双極) の厳密解 tdrift_bz1、
    K1≠0 は ソレノイド + 四極 の厳密解 tsolqu の bz≠0 の枝 (双極は四極の中心のずれとして入れる)"""
    if al == 0.0:
        return
    if ak1 == 0.0:
        for i in range(n):
            tdrift_bz1(i, S, al, bz, ak0x, ak0y, False, False, rs, 0.0, 0.0, 0.0)
    else:
        # 四極 + 双極 = 中心が (-a, -b) にずれた四極: By + i Bx = (K1/L)((x + a) + i (y + b))、a = Re K0/K1、b = Im K0/K1。
        # 座標をずらすと、一様な bz のベクトルポテンシャルのゲージが変わる分だけ正準運動量もずれる (px - bz b/2、py + bz a/2)
        a = ak0x / ak1
        b = ak0y / ak1
        if a != 0.0 or b != 0.0:
            for i in range(n):
                bzp = bz / (1.0 + S[i, 5])
                S[i, 0] += a
                S[i, 2] += b
                S[i, 1] -= bzp * b * 0.5
                S[i, 3] += bzp * a * 0.5
        tsolqu_bz_all(S, n, al, ak1, bz, 0.0, False, False, rs, 0.0, 0.0, 0.0, 0.0)
        if a != 0.0 or b != 0.0:
            for i in range(n):
                bzp = bz / (1.0 + S[i, 5])
                S[i, 0] -= a
                S[i, 2] -= b
                S[i, 1] += bzp * b * 0.5
                S[i, 3] -= bzp * a * 0.5


@njit(**NJ)
def sol_k0_kick_all(S, n, kre, kim):
    for i in range(n):
        pr = 1.0 + S[i, 5]
        S[i, 1] -= kre / pr
        S[i, 3] += kim / pr


@njit(**NJ)
def sol_setr0_kin(S, n, bz):
    """放射の基準点 (力学的運動量の向き)"""
    for i in range(n):
        bzp = bz / (1.0 + S[i, 5])
        S[i, 7] = S[i, 1] + bzp * S[i, 2] * 0.5
        S[i, 8] = S[i, 3] - bzp * S[i, 0] * 0.5
        S[i, 9] = S[i, 4]


@njit(**NJ)
def sol_tradk_all(S, n, al, bz, fluc, rs, p0, anrad, cuc):
    """tradk を力学的運動量で (正準 → 力学 → tradk → 基準点の更新 → 正準。正準へ戻すときの p は放射の後の値)"""
    for i in range(n):
        bzp = bz / (1.0 + S[i, 5])
        S[i, 1] += bzp * S[i, 2] * 0.5
        S[i, 3] -= bzp * S[i, 0] * 0.5
        tradk_core1(i, S, al, 0.0, fluc, rs, p0, anrad, cuc)
        setr0_1(i, S)
        bzp = bz / (1.0 + S[i, 5])
        S[i, 1] -= bzp * S[i, 2] * 0.5
        S[i, 3] += bzp * S[i, 0] * 0.5


def sol_mult_full(S, n, ph, pr, bz, trf0, trpt, rfsw, rad, fluc, rs, fr_in=True, fr_out=True, hook=None, kick=None, wk=None):
    """v2: ソレノイド領域 (bz≠0) の中の MULT と CAVI (tsol.f の icMULT → tmulti1、icCAVI → tmultiacc を bz つきで)。
    K64 の tmulti の骨組み (半分の本体、[多極の蹴り、RF の蹴り、本体] × ndiv、端のフリンジ) はそのままで、本体を
    「ソレノイド + 四極」または「ソレノイド + 双極」の厳密解にする。蹴りは位置を変えないので正準運動量にそのまま足せる。
    放射 (RAD) は本体と蹴りを合わせた向きの変化から、分割ごとに力学的運動量で tradk。
    C2: TRPT の加速 (vnominal = VOLT·sin(PHI)≠0: tmultiacc と同じ重み ws の分割、空洞の後の基準運動量 p2 に規格化し直して
    新しい ph を返す。bz はこの要素の入口の運動量で規格化した値のまま)、wake (kick・wk: 分割ごとに txwake)、
    AUTOPHI (rf[5]: 位相の基準を設計軌道の到着時刻にする。rf[3] の phic に w·DZ が入っている)。
    戻り値: 加速したときは新しい ph、そうでなければ None"""
    (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps0, radok, rf) = pr
    p0, h0 = ph.p0, ph.h0
    use_rf = rf is not None and (trpt or rfsw)
    vnominal = 0.0
    aphi = False
    if use_rf:
        v, w, wi, phic, vnom_tp = rf[:5]
        aphi = len(rf) > 5 and bool(rf[5])
        vnominal = vnom_tp if trpt else 0.0
    krad = rad and radok
    nz = al != 0.0
    mult_rot_in(S, n, dx, dy, ct, st)
    akr1 = akr[1].real if nmmax >= 1 else 0.0
    if akr1 < 0.0:
        raise NotImplementedError('SOL 領域内の MULT: K1 の符号 (theta2 の回転後に負)')
    akr0 = akr[0]
    ndiv = max(ndiv0, 1 + int(abs(bz * al) / 0.1))
    if nz and krad:
        sol_setr0_kin(S, n, bz)
        if hook is not None:
            hook.setr0(S, n, bz)
        ndiv = max(ndiv, 1 + int(abs(bz * al) / 0.02))     # 放射の刻み (k_sol_bend と同じ理由)
        ndiv = max(ndiv, ndivrad(abs(akr0), akr1, bz, eps0, h0 * (10.0 if trpt else 1.0)))
    if use_rf:
        ndiv = max(ndiv, 1 + int(abs(w * al)))
    ndiv = min(1000, ndiv)
    uniform = True
    if vnominal != 0.0:
        # tmultiacc: 加速があるときは、運動量が指数関数的に変わるのに合わせて重み ws を等比にする
        epsr = 1.0 if eps0 == 0.0 else eps0
        he = h0 + vnominal
        r = abs(vnominal * (1.0 / h0 + 1.0 / he))
        ndiv = min(1000, max(ndiv, 1 + int(min(abs(w * al), math.sqrt(r**2 / 3.0 / (epsr * 1e-3))))))
        if abs(r) > 1e-6:
            wl = math.asinh(r * (2.0 + r) / 2.0 / (1.0 + r))
            ndiv = min(1000, max(ndiv, int(math.floor(wl * 10.0 / epsr + 0.5))))
            r1 = wl / ndiv / 2.0
            ws1 = 2.0 * math.exp(r1) * math.sinh(r1)
            we = 1.0 + ws1
            ws = [0.0] * (ndiv + 2)
            ws[1] = ws1 / r
            for i in range(2, ndiv + 1):
                ws[i] = ws[i - 1] * we
            uniform = False
    if uniform:
        ws = [0.0] + [1.0 / ndiv] * ndiv + [0.0]
    ws[ndiv + 1] = 0.0
    if use_rf:
        pe = _n._h2p(h0 + vnominal) if vnominal != 0.0 else p0
        vcorr = v * (w * (0.5 / p0 + 0.5 / pe))**2 / 4.0
        phis = 0.0 if (trpt or aphi) else w * trf0
        dphis = phis - phic
        if rad or trpt or aphi:
            offset = math.sin(dphis)
            offset1 = 0.0
        else:
            offset = math.sin(dphis) - math.sin(phis)
            offset1 = math.sin(phis)
    ak1 = akr1 * ws[1] * 0.5
    al1 = al * ws[1] * 0.5
    ak01 = akr0 * ws[1] * 0.5
    if nz and fringe and fr_in:
        if use_rf:
            _n.cav_frin_all(S, n, al, v, w, dphis, offset, p0, h0)
        mult_fringe_all(S, n, akr, dofr, nmmax, al, 1.0, False)
    fw0 = 2.0 * kick.fw if kick is not None else 0.0          # kick_factor は 1/2 込み
    fw_save = kick.fw if kick is not None else 0.0
    sv = 0.0
    for m in range(1, ndiv + 1):
        wsm = ws[m]
        if kick is not None:
            # wake の重み: K64 (tcav.f・tmulti.f) は fw0·(ws(n) + ws(n−1)·0.5) で、分割が 2 以上だと合計が 1.5 倍になる。v2 (C2) は台形の重み
            kick.fw = fw0 * ws[1] * 0.5 if m == 1 else fw0 * (wsm + ws[m - 1]) * 0.5
            kick.apply(S, n, wk[0], wk[1], 0.0, 0.0, m == 1, p0, h0)
        akrm = akr[:nmmax + 1] * wsm
        if nz:
            sol_body_all(S, n, al1, ak1, bz, ak01.real, ak01.imag, rs)
            if krad:
                if hook is not None:
                    hook.tradk(S, n, al1, bz)              # track6 の経路 (スピン追跡): 放射とスピンの歳差を Python の tradk で
                else:
                    sol_tradk_all(S, n, al1, bz, fluc, rs, p0, ph.anrad, ph.cuc)
            wm = 0.5 * (ws[m] + ws[m + 1])
            al1 = al * wm
            ak1 = akr1 * wm
            ak01 = akr0 * wm
            if nmmax >= 2:
                mult_kick_all(S, n, akrm, nmmax, 2)
        elif nmmax >= 1:
            mult_kick_all(S, n, akrm, nmmax, 1)
        elif akrm[0] != 0.0:
            sol_k0_kick_all(S, n, akrm[0].real, akrm[0].imag)
        if use_rf:
            dzn = 0.0
            if vnominal != 0.0:
                sv += vnominal * wsm
                h2 = h0 + sv
                p2 = _n._h2p(h2)
                dp2 = sv * (h2 + h0) / (p2 + p0) / p0
                pr2 = 1.0 + dp2
                dvn = -dp2 * (1.0 + pr2) / h2 / (h2 + pr2 * h0)
                dzn = -dvn * al * wsm * 0.5
            mult_rf_kick_all(S, n, v, vcorr, w, wi, dphis, offset, offset1, wsm, dzn, p0, h0)
    if kick is not None:
        kick.fw = fw0 * ws[ndiv] * 0.5
        kick.apply(S, n, wk[0], wk[1], 0.0, 0.0, False, p0, h0)
        kick.fw = fw_save
    if nz:
        sol_body_all(S, n, al1, ak1, bz, ak01.real, ak01.imag, rs)
        if fringe and fr_out:
            mult_fringe_all(S, n, akr, dofr, nmmax, al, -1.0, True)
            if use_rf:
                _n.cav_frin_all(S, n, al, -v, w, dphis, offset, p0, h0)
        if krad:
            if hook is not None:
                hook.last = True
                hook.tradk(S, n, al1, bz)
            else:
                sol_tradk_all(S, n, al1, bz, fluc, rs, p0, ph.anrad, ph.cuc)
    mult_rot_out(S, n, dx, dy, ct, st)
    if vnominal != 0.0:
        h2 = h0 + vnominal
        p2 = _n._h2p(h2)
        dp2p2 = vnominal * (h2 + h0) / (p2 + p0) / p2
        # 正準運動量 px = (力学的) − bz·y/2 の bz は入口の p0 で規格化した値。出口では p2 で規格化し直す: S[:,1] は p で割った値
        # (粒子の p は変わらない) なので、そのままで領域の次の要素 (bz·p0/p2) の正準運動量になっている
        mult_exit_momentum(S, n, p0, p2, dp2p2, h2)
        ph2 = copy.copy(ph)
        ph2.set_pgev(p2 * ph.amass)
        return ph2
    return None


# ============================================================================ ラティスの表 (track6n.Prog の拡張)
def compile_mult(el, ph=None, dz=0.0):
    """MULT → (K_MULT, params) / nmmax<0 → (K_DRIFT, (L,))。track6 _tmulti_ens と同じ前処理 (_lat._mult_params)。
    VOLT+DVOLT≠0 の MULT は ph (基準運動量 ph.omega0、charge、amass) を渡すと RF の引数 rf (tmulti1 → tmultiacc) を付ける。
    実行時に (TRPT か RFSW) なら tmultiacc (mult_full)、そうでなければ RF なしの tmulti (k_mult)"""
    g = _t6._g
    pm = _lat._mult_params(el)
    al, nmmax, akr, ndiv = pm['al'], pm['nmmax'], pm['akr'], pm['ndiv']
    if g(el, 'ANGLE'):
        from . import multa as _ma
        return _ma.compile_multa(el, ph)
    if g(el, 'ANGLE') or g(el, 'F1') or g(el, 'FRINGE') or any(g(el, k) for k in ('F2', 'F1K1F', 'F1K1B', 'F2K1F', 'F2K1B', 'FB1', 'FB2', 'AE1', 'AE2')):
        # F1・FRINGE・FB1・FB2・AE1・AE2 つきの MULT は track6 (Python の _tmulti_ens) が扱う (呼び出し側が track6 に落とす。v2next P3)
        raise NotImplementedError('track6n_ext: MULT の曲げ・F1・FRINGE・AE は未対応')
    vc = g(el, 'VOLT') + g(el, 'DVOLT')
    rf = None
    if vc != 0.0:
        if ph is None:
            raise NotImplementedError('track6n_ext: MULT の RF は未対応')
        aphi = g(el, 'AUTOPHI') != 0.0
        if aphi and not _compat.v2():
            raise NotImplementedError('track6n_ext: MULT の AUTOPHI は未対応 (位置の z が要る)')
        harm = g(el, 'HARM')
        w = 2.0 * math.pi * g(el, 'FREQ') / _n.CVELOC if harm == 0.0 else ph.omega0 * harm / _n.CVELOC
        sgn = 1.0 if ph.charge >= 0 else -1.0
        # tffs.f:1393: vnominal = VOLT/amass |charge| sin(-PHI sign(charge)) (TRPT のときだけ。DVOLT は含まない)
        rf = (vc / ph.amass * abs(ph.charge), w, (1.0 / w if w != 0.0 else 0.0), (g(el, 'PHI') + g(el, 'DPHI')) * ph.charge,
              g(el, 'VOLT') / ph.amass * abs(ph.charge) * math.sin(-g(el, 'PHI') * sgn))
        if aphi:
            # v2 (C2): AUTOPHI (tturn.f: ph = DPHI + DZ·w)
            rf = (rf[0], rf[1], rf[2], rf[3] + dz * w * ph.charge, rf[4], True, dz)
    if nmmax < 0 and rf is None:
        return K_DRIFT, (al,)
    if nmmax < 0:
        nmmax = 0
        akr = [0.0]
    eps = pm['eps']
    if al != 0.0:
        dofr = [akr[m] != 0 and (m <= 2 or abs(akr[m]) / al * 0.05**(m + 1) / (m + 1) > eps)
                for m in range(nmmax + 1)]
    else:
        dofr = [False] * (nmmax + 1)     # al=0 ではフリンジを掛けない (track6 では未使用)
    ct, st = math.cos(pm['theta2']), math.sin(pm['theta2'])
    return K_MULT, (al, nmmax, np.array(akr, dtype=np.complex128), np.array(dofr, dtype=np.bool_), ndiv,
                    pm['dx'], pm['dy'], ct, st, pm['fringe'], g(el, 'EPS'),
                    g(el, 'RAD') == 0.0 and al != 0.0, rf)


# ============================================================================ MULT の RF (tmulti.f tmultiacc, ANGLE=0)
@njit(**NJ)
def mult_qu_all(S, n, al1, ak1, eps0, ak0x, ak0y, sch):
    for i in range(n):
        tsolqu_m1(i, S, al1, ak1, eps0, ak0x, ak0y, False, sch)           # tmultiacc は tsolqu (tsolqum ではない)


@njit(**NJ)
def mult_kick_all(S, n, akrm, nmmax, nmin):
    for i in range(n):
        mult_kick1(i, S, akrm, nmmax, nmin)


@njit(**NJ)
def mult_rot_in(S, n, dx, dy, ct, st):
    for i in range(n):
        S[i, 0] -= dx
        S[i, 2] -= dy
        rot1(i, S, ct, st)


@njit(**NJ)
def mult_rot_out(S, n, dx, dy, ct, st):
    for i in range(n):
        rot1(i, S, ct, -st)
        S[i, 0] += dx
        S[i, 2] += dy


@njit(**NJ)
def mult_fringe_all(S, n, akr, dofr, nmmax, al, sgn, rev):
    """ttfrins: 入口は n=0..nmmax の順に +akr、出口は n=nmmax..0 の順に −akr"""
    for i in range(n):
        if rev:
            for m in range(nmmax, -1, -1):
                if dofr[m]:
                    ttfrins1(i, S, m * 2 + 2, sgn * akr[m].real, sgn * akr[m].imag, al)
        else:
            for m in range(nmmax + 1):
                if dofr[m]:
                    ttfrins1(i, S, m * 2 + 2, sgn * akr[m].real, sgn * akr[m].imag, al)


@njit(**NJ)
def mult_rf_kick_all(S, n, v, v20a, w, wi, dphis, offset, offset1, wsn, dzn, p0, h0):
    """tmultiacc の do concurrent (電圧キック 1 分割。v1a=v11a=0、v20a=v02a)"""
    pmin = 1e-10
    oneev = 1.0 + 1e-6
    tlim = 1e4
    for i in range(n):
        dp1r = S[i, 5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        if S[i, 6] > 0.1:
            h1 = p2h(p1)
        else:
            h1 = p1r * h0 / (1.0 - S[i, 6] + S[i, 10])
        S[i, 4] -= dzn
        t = fmin(tlim, fmax(-tlim, -S[i, 4] * h1 / p1))
        ph = 0.5 * w * t
        dh = fmax(oneev - h1, (v + v20a * S[i, 0] * S[i, 0] + v20a * S[i, 2] * S[i, 2])
                  * (-2.0 * math.sin(ph) * math.cos(ph - dphis) + offset) * wsn)
        h2 = h1 + dh
        ah = fmax(dh * (h1 + h2) / (p1 * p1), -1.0 + pmin)
        dpr = sqrt1(ah)
        dp2r = fmax(dp1r + p1r * dpr, pmin - 1.0)
        p2r = 1.0 + dp2r
        S[i, 5] = dp2r
        S[i, 6] = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + S[i, 10]
        S[i, 4] = -t * p2r * p0 / h2 - dzn
        av = -(math.cos(2.0 * ph - dphis) * wi - offset1 * t) * wsn / p0
        dpx = 2.0 * v20a * S[i, 0] * av
        dpy = 2.0 * v20a * S[i, 2] * av
        S[i, 1] = (S[i, 1] * p1r + dpx) / p2r
        S[i, 3] = (S[i, 3] * p1r + dpy) / p2r


@njit(**NJ)
def mult_exit_momentum(S, n, p0, p2, dp2p2, h2):
    pmin = 1e-10
    for i in range(n):
        dp = fmax(S[i, 5], pmin - 1.0)
        dp1 = dp * p0 / p2 - dp2p2
        pr1 = 1.0 + dp1
        S[i, 5] = dp1
        h1 = p2h(p2 * pr1)
        S[i, 6] = -dp1 * (1.0 + pr1) / h1 / (h1 + pr1 * h2) + S[i, 10]


def mult_full(S, n, ph, pr, trf0, trpt, rad, fluc, rs, kick, wk):
    """tmulti.f tmultiacc (ANGLE=0、F1/FRINGE なし)。RF のある MULT (TRPT か RFSW のとき)。
    kick (WakeKick) と wk があれば分割ごとに txwake。新しい ph (空洞の後の基準運動量) を返す (tturn が周の終わりに元へ戻す)"""
    import copy
    (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps0, radok, rf) = pr
    v, w, wi, phic, vnom_tp = rf[:5]
    aphi = len(rf) > 5 and bool(rf[5])          # v2 (C2): AUTOPHI (位相の基準は設計軌道の到着時刻。phic に w·DZ が入っている)
    vnominal = vnom_tp if trpt else 0.0
    p0, h0 = ph.p0, ph.h0
    krad = rad and radok
    nz = al != 0.0
    mult_rot_in(S, n, dx, dy, ct, st)
    akr1 = akr[1].real if nmmax >= 1 else 0.0
    akr0 = akr[0]
    ndiv = ndiv0
    epsr = 1.0 if eps0 == 0.0 else eps0
    ndivmax = 1000
    if nz:
        if krad:
            for i in range(n):
                setr0_1(i, S)
            ndiv = max(ndiv, ndivrad(abs(akr0), akr1, 0.0, eps0, h0 * (10.0 if trpt else 1.0)))
        ndiv = min(ndivmax, ndiv)
    he = h0 + vnominal
    pe = _n._h2p(he)
    vcorr = v * (w * (0.5 / p0 + 0.5 / pe))**2 / 4.0
    r = abs(vnominal * (1.0 / h0 + 1.0 / he))
    ndiv = min(ndivmax, max(ndiv, 1 + int(min(abs(w * al), math.sqrt(r**2 / 3.0 / (epsr * 1e-3))))))
    if abs(r) > 1e-6:
        wl = math.asinh(r * (2.0 + r) / 2.0 / (1.0 + r))
        ndiv = min(ndivmax, max(ndiv, int(math.floor(wl * 10.0 / epsr + 0.5))))
        r1 = wl / ndiv / 2.0
        ws1 = 2.0 * math.exp(r1) * math.sinh(r1)
        we = 1.0 + ws1
        ws = [0.0] * (ndiv + 2)
        ws[1] = ws1 / r
        for i in range(2, ndiv + 1):
            ws[i] = ws[i - 1] * we
    else:
        ws = [0.0] + [1.0 / ndiv] * ndiv + [0.0]
    ws[ndiv + 1] = 0.0
    phis = 0.0 if (trpt or aphi) else w * trf0
    dphis = phis - phic
    if rad or trpt or aphi:
        offset = math.sin(dphis)
        offset1 = 0.0
    else:
        offset = math.sin(dphis) - math.sin(phis)
        offset1 = math.sin(phis)
    ak1 = akr1 * ws[1] * 0.5
    al1 = al * ws[1] * 0.5
    ak01 = akr0 * ws[1] * 0.5
    if nz and fringe:
        _n.cav_frin_all(S, n, al, v, w, dphis, offset, p0, h0)
        mult_fringe_all(S, n, akr, dofr, nmmax, al, 1.0, False)
    fw0 = 2.0 * kick.fw if kick is not None else 0.0          # kick_factor は 1/2 込み
    fw_save = kick.fw if kick is not None else 0.0
    sv = 0.0
    for m in range(1, ndiv + 1):
        wsm = ws[m]
        if kick is not None:
            # wake の重み: K64 (tcav.f・tmulti.f) は fw0·(ws(n) + ws(n−1)·0.5) で、分割が 2 以上だと合計が 1.5 倍になる。v2 (C2) は台形の重み
            kick.fw = fw0 * ws[1] * 0.5 if m == 1 else (fw0 * (wsm + ws[m - 1]) * 0.5 if _compat.v2() else fw0 * (wsm + ws[m - 1] * 0.5))
            kick.apply(S, n, wk[0], wk[1], 0.0, 0.0, m == 1, p0, h0)
        akrm = akr[:nmmax + 1] * wsm
        if nz:
            mult_qu_all(S, n, al1, ak1, eps0, ak01.real, ak01.imag, _qb.sch_of(rs))
            if krad:
                tradk_all(S, n, al1, 0.0, fluc, rs, p0, ph.anrad, ph.cuc)
            wm = 0.5 * (ws[m] + ws[m + 1])
            al1 = al * wm
            ak1 = akr1 * wm
            ak01 = akr0 * wm
            if nmmax >= 2:
                mult_kick_all(S, n, akrm, nmmax, 2)
        elif nmmax >= 1:
            mult_kick_all(S, n, akrm, nmmax, 1)
        dzn = 0.0
        if vnominal != 0.0:
            sv += vnominal * wsm
            h2 = h0 + sv
            p2 = _n._h2p(h2)
            dp2 = sv * (h2 + h0) / (p2 + p0) / p0
            pr2 = 1.0 + dp2
            dvn = -dp2 * (1.0 + pr2) / h2 / (h2 + pr2 * h0)
            dzn = -dvn * al * wsm * 0.5
        mult_rf_kick_all(S, n, v, vcorr, w, wi, dphis, offset, offset1, wsm, dzn, p0, h0)
    if kick is not None:
        kick.fw = fw0 * ws[ndiv] * 0.5
        kick.apply(S, n, wk[0], wk[1], 0.0, 0.0, False, p0, h0)
        kick.fw = fw_save
    if nz:
        mult_qu_all(S, n, al1, ak1, eps0, ak01.real, ak01.imag, _qb.sch_of(rs))
        if fringe:
            mult_fringe_all(S, n, akr, dofr, nmmax, al, -1.0, True)
            _n.cav_frin_all(S, n, al, -v, w, dphis, offset, p0, h0)
        if krad:
            tradk_all(S, n, al1, 0.0, fluc, rs, p0, ph.anrad, ph.cuc)
    mult_rot_out(S, n, dx, dy, ct, st)
    if vnominal != 0.0:
        h2 = h0 + vnominal
        p2 = _n._h2p(h2)
        dp2p2 = vnominal * (h2 + h0) / (p2 + p0) / p2
        mult_exit_momentum(S, n, p0, p2, dp2p2, h2)
        ph2 = copy.copy(ph)
        ph2.set_pgev(p2 * ph.amass)
        return ph2
    return ph


def compile_sol_quad(el, bzs):
    """ソレノイド内の QUAD → (K_SOL_QUAD, params) (track6 _tquad_sol と同じ制限)"""
    g = _t6._g
    L = el.L
    K1L = g(el, 'K1') * L
    if L == 0.0 or K1L == 0.0:
        raise NotImplementedError('track6n_ext: ソレノイド内の L=0 / K1=0 の四極')
    if any(g(el, k) for k in ('F1', 'FRINGE')):
        raise NotImplementedError('track6n_ext: ソレノイド内の四極の F1 / FRINGE')
    ak = abs(K1L)
    th2 = g(el, 'ROTATION') + (0.5 * math.pi if K1L < 0 else 0.0)
    return K_SOL_QUAD, (L, ak, math.cos(th2), math.sin(th2), g(el, 'DISFRIN') == 0.0, bzs,
                        g(el, 'RAD') == 0.0 and L != 0.0, g(el, 'DX'), g(el, 'DY'))


def compile_sol_region(prog, elements, i, ph):
    """elements[i] (BOUND≠0 の SOL) から終点 (次の BOUND≠0 の SOL) までを compile して prog に足し、終点の添字を返す
    (track6 _tsol_region と同じ制限)。track6n.Prog.__init__ から呼ばれる"""
    g = _t6._g
    sb = elements[i]
    v2 = _compat.v2()
    if not v2:
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
    try:
        ph.p0_ref = ph.p0                    # この表の bz を規格化した基準運動量 (step_kind が加速の後に p0_ref/p0 を掛ける)
    except AttributeError:
        pass
    fent = fext = None
    if v2:
        # C2: SOL のずれ・回転 (tsgeo が決める DX DY DZ CHI1..3) と F1 の放射
        from . import sgeo as _sg
        fent, fext = _sg.region_frames(elements, i, ph)

    def edge(el, rot, frame, bz0, bz1):
        dbz = (bz1 - bz0) if g(el, 'DISFRIN') == 0.0 else 0.0
        f1 = g(el, 'F1') if v2 else 0.0
        if frame is None and f1 == 0.0:
            return K_SOL_FRIN, (dbz,)
        r, dx, dy, dz = frame if frame is not None else (np.eye(3), 0.0, 0.0, 0.0)
        return K_SOL_EDGE, (rot if frame is not None else 0, r, dx, dy, dz, dbz, bz0, bz1, f1,
                            g(el, 'RAD') == 0.0 and f1 != 0.0 and bz0 != bz1)

    k_, p_ = edge(sb, 1, fent, 0.0, bzs)
    prog._add(k_, p_, sb.name, True, SM_START)
    for l in range(i + 1, ke + 1):
        el = elements[l]
        name = type(el).__name__
        sm = SM_END if l == ke else SM_IN
        _t6._sol_check_type(el)
        if isinstance(el, Drift) and getattr(el, '_apert', False):
            # tsol.f の icAprt → tapert1 (領域の外と同じ。x, y だけを見るので正準運動量のままでよい)。la は変えない (track6n.turn)
            prog._add(_n.K_APERT, (_n.apert_params(el),), el.name, True, sm)
            continue
        if isinstance(el, Drift):
            prog._add(K_SOL_DRIFT, (el.L, bzs, g(el, 'RAD') == 0.0), el.name, True, sm)
        elif isinstance(el, Quadrupole):
            if v2 and g(el, 'K1') * el.L == 0.0 and not any(g(el, k) for k in ('F1', 'FRINGE')):
                # C2: K1=0 (か L=0) の四極はソレノイドの中のドリフト (tquad の ak=0)
                prog._add(K_SOL_DRIFT, (el.L, bzs, g(el, 'RAD') == 0.0), el.name, True, sm)
            else:
                k, p = compile_sol_quad(el, bzs)
                prog._add(k, p, el.name, True, sm)
        elif name == 'Solenoid':
            bz1 = 0.0 if l == ke else _t6._tfbzs(el, ph)
            if not v2 and (any(g(el, k) for k in ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')) or g(el, 'F1') != 0.0):
                raise NotImplementedError('track6n_ext: SOL のずれ・回転・F1')
            k_, p_ = edge(el, -1, fext if l == ke else None, bzs, bz1)
            prog._add(k_, p_, el.name, True, sm)
            bzs = bz1
        elif name in ('Bend', 'Multipole', 'Cavity') and _compat.v2():
            _compile_sol_v2(prog, el, name, bzs, ph, sm, l)
        elif name == 'BeamBeam':
            raise NotImplementedError('track6n_ext: ソレノイド領域内の BEAMBEAM (K64 の tsol が SIGSEGV で落ちるので比較できない)')
        elif el.L == 0.0:
            prog._add(K_NONE, (), el.name, True, sm)
        else:
            raise NotImplementedError(f'track6n_ext: ソレノイド内の要素 {name} は未対応 (DRIFT/QUAD のみ)')
    return ke


def _compile_sol_v2(prog, el, name, bzs, ph, sm, idx=None):
    """v2: ソレノイド領域内の BEND・MULT・CAVI (docs/30 B69、docs/31)。K64 の tsol.f の呼び出しの意図どおり:
    BEND は直線の座標のまま一様な双極磁場 (tdrift(al, bzs, phiy, phix))、MULT は tmulti1、CAVI は tmultiacc (どちらも bz つき)。
    gs (この K64 のビルド) は領域内の偏向をほぼ無視するので、根拠は一様磁場の中の運動方程式の数値積分との一致 (test_v2_C_sol.py)"""
    g = _t6._g
    if name == 'Bend':
        if g(el, 'K1') != 0.0 or g(el, 'DX') != 0.0 or g(el, 'DY') != 0.0:
            raise NotImplementedError(f'SOL 領域内の BEND {el.name}: K1・DX・DY は未対応 (tsol.f も使わない。四極成分は QUAD か MULT で入れる)')
        theta = g(el, 'ROTATE') + g(el, 'DROTATE')
        phi = float(el.ANGLE) + g(el, 'K0')
        al = float(el.L)
        if al == 0.0:
            # 長さ 0 の双極: 蹴りだけ
            prog._add(K_SOL_MULT, ((0.0, 0, np.array([complex(phi * math.cos(theta), phi * math.sin(theta))]), np.array([False]), 1,
                                    0.0, 0.0, 1.0, 0.0, False, 0.0, False, None), bzs), el.name, True, sm)
            return
        prog._add(K_SOL_BEND, (al, bzs, phi * math.cos(theta), phi * math.sin(theta), g(el, 'RAD') == 0.0), el.name, True, sm)
        return
    if name == 'Multipole':
        if g(el, 'ANGLE') != 0.0:
            raise NotImplementedError(f'SOL 領域内の MULT {el.name}: ANGLE≠0 は未対応')
        k, p = compile_mult(el, ph, dz=(_n.autophi_dz(getattr(prog, 'l0', 0) + idx) if idx is not None else 0.0))
        if k == K_DRIFT:
            prog._add(K_SOL_DRIFT, (p[0], bzs, g(el, 'RAD') == 0.0), el.name, True, sm)
        elif bzs == 0.0:
            prog._add(k, p, el.name, True, sm)               # bz=0 の領域内は普通の MULT
        else:
            prog._add(K_SOL_MULT, (p, bzs), el.name, True, sm)
        return
    # Cavity
    sub = _n.Prog([el], ph, l0=-10**9)                # (AUTOPHI の DZ はここでは足さない: 下で領域の中の位置の DZ を足す)
    k, p = sub.kinds[0], sub.params[0]
    aphi = g(el, 'AUTOPHI') != 0.0 and k == _n.K_CAV
    dz_ = _n.autophi_dz(getattr(prog, 'l0', 0) + idx) if (aphi and idx is not None) else 0.0
    if bzs == 0.0 or k != _n.K_CAV:
        if aphi:
            p = p[:4] + (p[4] + dz_ * p[2] * ph.charge,) + p[5:]
            if hasattr(prog, 'autophi'):
                prog.autophi[len(prog.kinds)] = True
        prog._add(k, p, el.name, True, sm)
        return
    (al, v, w, wi, phic, fringe, mfring, vnom, vq) = p[:9]
    if any(q != 0.0 for q in vq):
        raise NotImplementedError(f'SOL 領域内の CAVI {el.name}: V1・V20・V11・V02 は未対応')
    rot = g(el, 'ROTATION')
    # v2 (E): 領域内の CAVI も放射とスピンの歳差の刻み (hook) を通す。K64 の tsol は CAVI を krad=.false. で tmultiacc に渡すので、
    # ソレノイド磁場の中の空洞のスピンが歳差しない (縦磁場 (1+a)bz·L のぶん、Emittance[] (solemit) と追跡が食い違う)。v2 は物理的に歳差・放射させる
    prm = (al, 0, np.array([0j]), np.array([False]), 1, g(el, 'DX'), g(el, 'DY'), math.cos(rot), math.sin(rot),
           fringe, 0.0, al != 0.0, ((v, w, wi, phic + dz_ * w * ph.charge, vnom, True, dz_) if aphi else (v, w, wi, phic, vnom)))
    prog._add(K_SOL_MULT, (prm, bzs), el.name, True, sm)


# track6n に統合済み: ProgExt / turn_ext は track6n.Prog / track6n.turn の別名 (互換のため残す)
ProgExt = _n.Prog


def step_kind(kind, pr, S, n, ph, trf0, rfsw, rad, fluc, rs, tm=1.0, kick=None, wk=None):
    """1 要素分のカーネル呼び出し (track6n.turn の分岐 + 新しい kind)。
    C2: 領域内の CAVI・MULT が TRPT で加速したときは新しい ph を返す (ほかは None)。SOL 領域の bz は compile のときの基準運動量
    (ph.p0_ref) で規格化してあるので、上流で加速していれば p0_ref/p0 を掛ける (K64 の tfbzs(l) は要素 l の入口の運動量で規格化)"""
    p0, h0, anrad, cuc = ph.p0, ph.h0, ph.anrad, ph.cuc
    hr = h0 * tm                                  # ndivrad だけに渡す (TRPT は 10 倍)
    pref = getattr(ph, 'p0_ref', None)
    if pref is not None and pref != p0 and kind in (K_SOL_FRIN, K_SOL_EDGE, K_SOL_DRIFT, K_SOL_BEND, K_SOL_MULT, K_SOL_QUAD):
        sc = pref / p0
        if kind == K_SOL_FRIN:
            pr = (pr[0] * sc,)
        elif kind == K_SOL_EDGE:
            pr = pr[:5] + (pr[5] * sc, pr[6] * sc, pr[7] * sc) + pr[8:]
        elif kind == K_SOL_DRIFT:
            pr = (pr[0], pr[1] * sc, pr[2])
        elif kind == K_SOL_BEND:
            pr = (pr[0], pr[1] * sc) + pr[2:]
        elif kind == K_SOL_MULT:
            pr = (pr[0], pr[1] * sc)
        else:
            pr = pr[:5] + (pr[5] * sc,) + pr[6:]
    if kind == K_MULT:
        (al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad, radok, rf) = pr
        k_mult(S, n, al, nmmax, akr, dofr, ndiv0, dx, dy, ct, st, fringe, eps_rad,
               rad and radok, fluc, rs, p0, hr, anrad, cuc)
    elif kind == K_SOL_FRIN:
        k_sol_frin(S, n, pr[0])
    elif kind == K_SOL_EDGE:
        (rot, r, dx, dy, dz, dbz, bz0, bz1, f1, radok) = pr
        k_sol_edge(S, n, rot, r, dx, dy, dz, dbz, bz0, bz1, f1, rad and radok, fluc, rs, p0, anrad, cuc)
    elif kind == K_SOL_DRIFT:
        (al, bz, radok) = pr
        k_sol_drift(S, n, al, bz, rad and radok, fluc, rs, p0, anrad, cuc)
    elif kind == K_SOL_BEND:
        (al, bz, ak0x, ak0y, radok) = pr
        k_sol_bend(S, n, al, bz, ak0x, ak0y, rad and radok, fluc, rs, p0, anrad, cuc)
    elif kind == K_SOL_MULT:
        return sol_mult_full(S, n, ph, pr[0], pr[1], trf0, tm != 1.0, rfsw, rad, fluc, rs, kick=kick, wk=wk)
    elif kind == _n.K_APERT:
        pass                                  # 行列・COD の写像 (solmap) では絞りは恒等。追跡は track6n.turn が扱う
    elif kind == K_SOL_QUAD:
        (L, ak, ct, st, fringe, bz, radok, qdx, qdy) = pr
        k_sol_quad(S, n, L, ak, ct, st, fringe, bz, rad and radok, fluc, rs, p0, hr, anrad, cuc, qdx, qdy)
    elif kind == _n.K_UND:
        from . import und as _un
        if pr[1] != 0.0 or pr[2] != 0.0:
            _un.k_und_field(S, n, *pr, p0, h0, rad, fluc, rs, anrad, cuc)
        else:
            _un.k_und(S, n, *pr, p0, h0)
    elif kind == _n.K_DRIFT:
        _n.k_drift(S, n, pr[0])
    elif kind == _n.K_QUAD:
        (L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0, radok) = pr
        _n.k_quad(S, n, L, K1L, dx, dy, ct, st, fringe, mfring, f1in, f2in, f1out, f2out, eps0,
                  rad and radok, fluc, rs, p0, hr, anrad, cuc)
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
                       fluc, rs, p0, hr, anrad, cuc)
        else:
            _n.k_bend0(S, n, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
                       fb1, fb2, mfring, fringe, krad, True, True, eps0, dx, dy, cost, sint, dtheta, dchi2, alg, phig,
                       fluc, rs, p0, hr, anrad, cuc)
    elif kind == _n.K_CAV:
        (al, v, w, wi, phic, fringe, mfring) = pr[:7]
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
        (al, phib, dx, dy, theta, cost, sint, tanp1, eps0, radok, sfb1, sfb2, sfrg, tanp2e, sv2) = pr
        _n.k_steer(S, n, al, phib, dx, dy, theta, cost, sint, tanp1, rad and radok, eps0, fluc, rs, ph.p0, ph.h0 * tm, ph.anrad, ph.cuc,
                   sfb1, sfb2, sfrg, tanp2e, sv2)


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
    rs = _n.make_rs(_t6.RNG.seed)
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
