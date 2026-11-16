"""
SynchroBetaEmittance (K64 の tfemit.f tfemits → temits.f temits・temits1・tefsetup・tesolvd・tevdif・tesumb …) の直訳。

シンクロトロン運動と色収差の下での平衡エミッタンス。運動量のずれ δ ごとに閉軌道・1 周写像・拡散行列を求め (2·ndims+1 点)、
(アクション J, 方位角モード m) の空間で平衡のビーム行列を反復で解く。
配列は Fortran の添字 (1 始まり) に合わせるため、内部の関数は 1 始まりの添字を引数に取り、numpy の配列には 0 始まりで入れる。
pysad に K64 の中間配列を取り出す手段がないので、検証は SynchroBetaEmittance の最終の 6 つの値 (nus, emitx, emity, emitxp, emityp, conv) だけ。
"""
from __future__ import annotations

import math

import numpy as np

PI2 = 2.0 * math.pi
C_LIGHT = 299792458.0


def iaidx(m, n):
    """packed 対称行列の添字 (1 始まり、m >= n)"""
    return m * (m - 1) // 2 + n


# 4×4 の下三角を packed 10 要素にする順番: (1,1)(2,1)(2,2)(3,1)(3,2)(3,3)(4,1)(4,2)(4,3)(4,4)
PK4 = [(k, l) for k in range(1, 5) for l in range(1, k + 1)]


def pack10(M):
    return np.array([M[k - 1, l - 1] for k, l in PK4])


def tetrig(amu, mphi2):
    m = np.arange(mphi2)
    phim = m * amu
    return np.cos(phim), np.sin(phim)


# ---------------------------------------------------------------------------
# tefsetup・補間
# ---------------------------------------------------------------------------

def teintp1(ip, f, tws, tw0, ri, ndims, helpers):
    """teintp1: tws の ip, ip+1 (Fortran の -ndims..ndims の添字) を f で補間した tr1(5,5) と h1(4)"""
    etwiss2ri, tinv6, mfit = helpers
    twf = f * tws[ip + 1 + ndims] + (1.0 - f) * tws[ip + ndims]
    twf = twf.copy()
    twf[14] = twf[10] * twf[13] - twf[11] * twf[12]
    rxi, _normal = etwiss2ri(twf)
    dnx = twf[2] - tw0[2]
    dny = twf[5] - tw0[5]
    dnz = twf[23] - tw0[23]
    cx, sx, cy, sy, cz, sz = math.cos(dnx), math.sin(dnx), math.cos(dny), math.sin(dny), math.cos(dnz), math.sin(dnz)
    rt = np.zeros((6, 6))
    rt[0] = cx * ri[0] + sx * ri[1]
    rt[1] = -sx * ri[0] + cx * ri[1]
    rt[2] = cy * ri[2] + sy * ri[3]
    rt[3] = -sy * ri[2] + cy * ri[3]
    rt[4] = cz * ri[4] + sz * ri[5]
    rt[5] = -sz * ri[4] + cz * ri[5]
    a = tinv6(rxi) @ rt
    tr1 = tmov65(a)
    h1 = twf[15:19].copy()
    return tr1, h1


def tmov65(a):
    b = np.zeros((5, 5))
    b[0:4, 0:4] = a[0:4, 0:4]
    b[0:4, 4] = a[0:4, 5]
    b[4, 4] = a[5, 5]
    return b


def tefsetup(ndp, mphi, mphi2, nz, ndims, beams, trads, tws, dj, dpndim, sigea, tw0, ri, helpers):
    """tefsetup: ha(4,4,mphi,ndp), hb(4,mphi2,ndp), ba(10,10,mphi,ndp), bb(10,mphi2,ndp), bd(10,4,mphi,ndp) (0 始まりで [.., m-1, j-1])"""
    nj = 2
    ba = np.zeros((10, 10, mphi, ndp))
    bb = np.zeros((10, mphi2, ndp))
    bd = np.zeros((10, 4, mphi, ndp))
    ha = np.zeros((4, 4, mphi, ndp))
    hb = np.zeros((4, mphi2, ndp))
    btr = np.zeros((10, 10))
    for j in range(1, ndp + 1):
        for j1 in range(-nj, nj + 1):
            aj = (j - 0.5) * dj + dj / (2 * nj + 1) * j1
            pj = math.sqrt(aj * 2.0)
            for i in range(1, nz + 1):
                phi = (i - 1) * math.pi / (nz - 1)
                cs = math.cos(phi)
                p = pj * cs
                ip = int((p + sigea) / dpndim) - ndims
                f = p / dpndim - ip
                tr1, h1 = teintp1(ip, f, tws, tw0, ri, ndims, helpers)
                # teintm
                f1 = (1.0 - f)
                trd = f1 * trads[ip + ndims] + f * trads[ip + 1 + ndims]
                tr1 = tr1 + trd
                # teintb
                beama = f1 * beams[ip + ndims] + f * beams[ip + 1 + ndims]
                for k in range(1, 5):
                    for k1 in range(1, k + 1):
                        kk = iaidx(k, k1) - 1
                        for l in range(1, 5):
                            for l1 in range(1, l + 1):
                                ll = iaidx(l, l1) - 1
                                if l == l1:
                                    btr[kk, ll] = tr1[k - 1, l - 1] * tr1[k1 - 1, l - 1]
                                else:
                                    btr[kk, ll] = (tr1[k - 1, l - 1] * tr1[k1 - 1, l1 - 1]
                                                   + tr1[k - 1, l1 - 1] * tr1[k1 - 1, l - 1])
                for m in range(1, mphi2 + 1):
                    phim = phi * (m - 1)
                    csm = math.cos(phim) / (nz - 1) / (2 * nj + 1)
                    if i == 1 or i == nz:
                        csm *= 0.5
                    for k in range(1, 5):
                        ha[0:4, k - 1, m - 1, j - 1] += tr1[0:4, k - 1] * csm
                        hb[k - 1, m - 1, j - 1] += h1[k - 1] * csm
                        for l in range(1, k + 1):
                            kk = iaidx(k, l) - 1
                            bb[kk, m - 1, j - 1] += h1[k - 1] * h1[l - 1] * csm
                            bd[kk, :, m - 1, j - 1] += (tr1[k - 1, 0:4] * h1[l - 1] + tr1[l - 1, 0:4] * h1[k - 1]) * csm
                    ba[:, :, m - 1, j - 1] += btr * csm
                    bb[:, m - 1, j - 1] += beama * csm
                for m in range(mphi2 + 1, mphi + 1):
                    phim = phi * (m - 1)
                    csm = math.cos(phim) / (nz - 1) / (2 * nj + 1)
                    if i == 1 or i == nz:
                        csm *= 0.5
                    for k in range(1, 5):
                        ha[:, k - 1, m - 1, j - 1] += tr1[0:4, k - 1] * csm
                        for l in range(1, k + 1):
                            kk = iaidx(k, l) - 1
                            bd[kk, :, m - 1, j - 1] += (tr1[k - 1, 0:4] * h1[l - 1] + tr1[l - 1, 0:4] * h1[k - 1]) * csm
                    ba[:, :, m - 1, j - 1] += btr * csm
    return ha, hb, ba, bb, bd


# ---------------------------------------------------------------------------
# 拡散の演算 (tevdif・tesolvd)
# ---------------------------------------------------------------------------

def _diff2(d, eps, aj, m):
    diff = d * eps * (m - 1) ** 2 * 0.125 / aj
    diff3 = -d * (m - 3) * (1.0 + eps / aj * 0.5 * (m - 1)) * 0.125
    return diff, diff3


def tespl(bx, nd, ndp, dj, e, damp, j0, m):
    """tespl: bx[nd, mphi2, ndp] の (j0, m) での (bdx, bddx)。j0, m は 1 始まり"""
    d = abs(2.0 * damp)
    ajj = (j0 - 0.5) * dj
    f2 = d * 0.5
    f1 = e * f2
    ajj1 = ajj - dj
    ajj2 = ajj + dj
    if j0 == 1:
        bx1 = bx[:, m - 1, 0]
        bx2 = bx[:, m - 1, 1]
    elif j0 == ndp:
        bx1 = bx[:, m - 1, j0 - 2]
        bx2 = np.zeros(nd)
    else:
        bx1 = bx[:, m - 1, j0 - 2]
        bx2 = bx[:, m - 1, j0]
    dy1 = (bx2 - bx1) / dj * 0.5
    bdx = f1 * dy1
    bddx = f2 * (-(e + ajj) * dy1 + e * (ajj1 * bx1 - 2.0 * ajj * bx[:, m - 1, j0 - 1] + ajj2 * bx2) / dj ** 2)
    if j0 == ndp:
        bddx = np.maximum(0.0, bddx)
    return bdx, bddx


def tevdifj(bc, bs, ndp, mphi2, nd, dj, eps, damp):
    for j in range(1, ndp + 1):
        for m in range(1, mphi2 + 1):
            fbc, dfbc = tespl(bc, nd, ndp, dj, eps, damp, j, m)
            fbs, dfbs = tespl(bs, nd, ndp, dj, eps, damp, j, m)
            if m == 1:
                bc[:, m - 1, j - 1] += 2.0 * dfbc
            else:
                bc[:, m - 1, j - 1] += dfbc
                bs[:, m - 1, j - 1] += dfbs
            m1 = m + 2
            if m1 <= mphi2:
                bc[:, m1 - 1, j - 1] += dfbc - m * fbc
                bs[:, m1 - 1, j - 1] += dfbs - m * fbs


def tevdifj1(bc, bs, ndp, mphi2, nd, dj, eps, damp):
    for j in range(ndp, 0, -1):
        for m in range(mphi2, 0, -1):
            fbc, dfbc = tespl(bc, nd, ndp, dj, eps, damp, j, m)
            fbs, dfbs = tespl(bs, nd, ndp, dj, eps, damp, j, m)
            m2 = m - 2
            if m == 3:
                bc[:, m2 - 1, j - 1] += dfbc + 2.0 * fbc
            elif m > 3:
                bc[:, m2 - 1, j - 1] += dfbc + (m - 2) * fbc
                bs[:, m2 - 1, j - 1] += dfbs + (m - 2) * fbs
            bc[:, m - 1, j - 1] += dfbc
            bs[:, m - 1, j - 1] += dfbs


def tesetdm(j, bd, hc, hs, mphi, mphi2):
    """tesetdm: bfb(20*mphi2)"""
    nsb = 10 * mphi2
    bfb = np.zeros(2 * nsb)
    for m in range(1, mphi2 + 1):
        mb = (m - 1) * 10
        ms = mb + nsb
        for i1 in range(1, 5):
            hci1 = hc[i1 - 1, 0, j - 1]
            bfb[mb:mb + 10] += bd[:, i1 - 1, m - 1, j - 1] * hci1
        for m1 in range(2, mphi2 + 1):
            k1 = abs(m - m1) + 1
            k2 = abs(m + m1 - 2) + 1
            for i1 in range(1, 5):
                hci1 = hc[i1 - 1, m1 - 1, j - 1]
                hsi1 = hs[i1 - 1, m1 - 1, j - 1]
                bfb[mb:mb + 10] += (bd[:, i1 - 1, k1 - 1, j - 1] + bd[:, i1 - 1, k2 - 1, j - 1]) * hci1
                bfb[ms:ms + 10] += (bd[:, i1 - 1, k1 - 1, j - 1] - bd[:, i1 - 1, k2 - 1, j - 1]) * hsi1
    return bfb


def _fwd_diag(bc, bs, aj, d, eps, mphi2):
    """j ごとの対角の拡散: m=2 と m>=3 (bc, bs は (nd, mphi2) のビュー)"""
    m = 2
    diff, diff3 = _diff2(d, eps, aj, m)
    bc[:, m - 1] += (diff + diff3) * bc[:, m - 1]
    bs[:, m - 1] += (diff - diff3) * bs[:, m - 1]
    for m in range(3, mphi2 + 1):
        diff = 1.0 + d * eps * (m - 1) ** 2 * 0.125 / aj
        bc[:, m - 1] *= diff
        bs[:, m - 1] *= diff


def tevdif(bb, bc, bs, ba, amuj, bd, hc, hs, ndp, mphi2, mphi, nd, dj, damp0, sige):
    """tevdif: bc, bs (nd, mphi2, ndp) をその場で更新する (1 周ぶんの拡散・移送)"""
    nsb = nd * mphi2
    eps = sige ** 2
    damp = damp0 * 0.5
    d = 2.0 * damp
    tevdifj(bc, bs, ndp, mphi2, nd, dj, eps, damp)
    tevdifj1(bc, bs, ndp, mphi2, nd, dj, eps, damp)
    for j in range(1, ndp + 1):
        aj = (j - 0.5) * dj
        cj, sj = bc[:, :, j - 1], bs[:, :, j - 1]       # ビュー
        _fwd_diag(cj, sj, aj, d, eps, mphi2)
        for m in range(1, mphi2 - 1):
            m2 = m + 2
            diff2 = -d * (m + 1) * (-1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
            cj[:, m2 - 1] += diff2 * cj[:, m - 1]
            sj[:, m2 - 1] += diff2 * sj[:, m - 1]
        for m in range(mphi2, 3, -1):
            m1 = m - 2
            diff1 = -d * (m - 3) * (1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
            cj[:, m1 - 1] += diff1 * cj[:, m - 1]
            sj[:, m1 - 1] += diff1 * sj[:, m - 1]
        _fwd_diag(cj, sj, aj, d, eps, mphi2)
    for j in range(1, ndp + 1):
        bfb = np.zeros(2 * nsb)
        if nd == 10:
            bfb = tesetdm(j, bd, hc, hs, mphi, mphi2)
        bfb[:nsb] += bb[:, :, j - 1].T.reshape(-1)
        cmu, smu = tetrig(amuj[j - 1], mphi2)
        for m in range(1, mphi2 + 1):
            mb = (m - 1) * nd
            ms = mb + nsb
            for i1 in range(1, nd + 1):
                bci1 = bc[i1 - 1, 0, j - 1]
                bsi1 = bs[i1 - 1, 0, j - 1]
                bfb[mb:mb + nd] += ba[:, i1 - 1, m - 1, j - 1] * bci1
                bfb[ms:ms + nd] += ba[:, i1 - 1, m - 1, j - 1] * bsi1
            for m1 in range(2, mphi2 + 1):
                k1 = abs(m - m1) + 1
                k2 = abs(m + m1 - 2) + 1
                for i1 in range(1, nd + 1):
                    bci1 = bc[i1 - 1, m1 - 1, j - 1]
                    bsi1 = bs[i1 - 1, m1 - 1, j - 1]
                    bfb[mb:mb + nd] += (ba[:, i1 - 1, k1 - 1, j - 1] + ba[:, i1 - 1, k2 - 1, j - 1]) * bci1
                    bfb[ms:ms + nd] += (ba[:, i1 - 1, k1 - 1, j - 1] - ba[:, i1 - 1, k2 - 1, j - 1]) * bsi1
        for m in range(1, mphi2 + 1):
            mb = (m - 1) * nd
            ms = mb + nsb
            dcb = bfb[mb:mb + nd].copy()
            dsb = bfb[ms:ms + nd].copy()
            bc[:, m - 1, j - 1] = cmu[m - 1] * dcb + smu[m - 1] * dsb
            bs[:, m - 1, j - 1] = -smu[m - 1] * dcb + cmu[m - 1] * dsb
    for j in range(1, ndp + 1):
        aj = (j - 0.5) * dj
        cj, sj = bc[:, :, j - 1], bs[:, :, j - 1]
        _fwd_diag(cj, sj, aj, d, eps, mphi2)
        for m in range(mphi2, 3, -1):
            m1 = m - 2
            diff1 = -d * (m - 3) * (1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
            cj[:, m1 - 1] += diff1 * cj[:, m - 1]
            sj[:, m1 - 1] += diff1 * sj[:, m - 1]
        for m in range(1, mphi2 - 1):
            m2 = m + 2
            diff2 = -d * (m + 1) * (-1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
            cj[:, m2 - 1] += diff2 * cj[:, m - 1]
            sj[:, m2 - 1] += diff2 * sj[:, m - 1]
        _fwd_diag(cj, sj, aj, d, eps, mphi2)
    tevdifj1(bc, bs, ndp, mphi2, nd, dj, eps, damp)
    tevdifj(bc, bs, ndp, mphi2, nd, dj, eps, damp)


def tesolvd(bb, bc, bs, ba, amuj, bd, hc, hs, ndp, mphi2, mphi, nd, dj, damp0, sige, tsolva):
    """tesolvd: 各 j で (2·nd·mphi2) の連立方程式を解いて bc, bs を求める"""
    nsb = nd * mphi2
    eps = sige ** 2
    damp = damp0 * 0.5
    d = 2.0 * damp
    bff = np.zeros((2 * nsb, 2 * nsb))
    for j in range(1, ndp + 1):
        aj = (j - 0.5) * dj
        cmu, smu = tetrig(amuj[j - 1], mphi2)
        for ma in range(1, mphi2 + 1):
            for ia in range(1, nd + 1):
                cj = np.zeros((nd, mphi2))
                sj = np.zeros((nd, mphi2))
                cj[ia - 1, ma - 1] = 1.0
                sj[ia - 1, ma - 1] = 1.0
                m = ma
                if m == 2:
                    diff, diff3 = _diff2(d, eps, aj, m)
                    cj[ia - 1, m - 1] += (diff + diff3) * cj[ia - 1, m - 1]
                    sj[ia - 1, m - 1] += (diff - diff3) * sj[ia - 1, m - 1]
                else:
                    diff = 1.0 + d * eps * (m - 1) ** 2 * 0.125 / aj
                    cj[ia - 1, m - 1] *= diff
                    sj[ia - 1, m - 1] *= diff
                m2 = ma + 2
                if m2 <= mphi2:
                    diff2 = -d * (m + 1) * (-1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
                    cj[ia - 1, m2 - 1] += diff2 * cj[ia - 1, m - 1]
                    sj[ia - 1, m2 - 1] += diff2 * sj[ia - 1, m - 1]
                    diff1 = -d * (m2 - 3) * (1.0 + eps / aj * 0.5 * (m2 - 1)) * 0.25
                    cj[ia - 1, ma - 1] += diff1 * cj[ia - 1, m2 - 1]
                    sj[ia - 1, ma - 1] += diff1 * sj[ia - 1, m2 - 1]
                if ma > 3:
                    m1 = ma - 2
                    diff1 = -d * (ma - 3) * (1.0 + eps / aj * 0.5 * (ma - 1)) * 0.25
                    cj[ia - 1, m1 - 1] += diff1 * cj[ia - 1, ma - 1]
                    sj[ia - 1, m1 - 1] += diff1 * sj[ia - 1, ma - 1]
                for m in range(ma - 2, ma + 3, 2):
                    if m == 2:
                        diff, diff3 = _diff2(d, eps, aj, m)
                        cj[ia - 1, m - 1] += (diff + diff3) * cj[ia - 1, m - 1]
                        sj[ia - 1, m - 1] += (diff - diff3) * sj[ia - 1, m - 1]
                    elif 3 <= m <= mphi2:
                        diff = 1.0 + d * eps * (m - 1) ** 2 * 0.125 / aj
                        cj[ia - 1, m - 1] *= diff
                        sj[ia - 1, m - 1] *= diff
                bfb = np.zeros(2 * nsb)
                for m in range(1, mphi2 + 1):
                    mb = (m - 1) * nd
                    ms = mb + nsb
                    for m1 in range(ma - 2, ma + 3, 2):
                        if m1 == 1:
                            bci1 = cj[ia - 1, 0]
                            bsi1 = sj[ia - 1, 0]
                            bfb[mb:mb + nd] += ba[:, ia - 1, m - 1, j - 1] * bci1
                            bfb[ms:ms + nd] += ba[:, ia - 1, m - 1, j - 1] * bsi1
                        elif 0 < m1 <= mphi2:
                            k1 = abs(m - m1) + 1
                            k2 = abs(m + m1 - 2) + 1
                            bci1 = cj[ia - 1, m1 - 1]
                            bsi1 = sj[ia - 1, m1 - 1]
                            bfb[mb:mb + nd] += (ba[:, ia - 1, k1 - 1, j - 1] + ba[:, ia - 1, k2 - 1, j - 1]) * bci1
                            bfb[ms:ms + nd] += (ba[:, ia - 1, k1 - 1, j - 1] - ba[:, ia - 1, k2 - 1, j - 1]) * bsi1
                maci = (ma - 1) * nd + ia - 1          # 0 始まり
                masi = maci + nsb
                for m in range(1, mphi2 + 1):
                    mc = (m - 1) * nd
                    ms = mc + nsb
                    dcb = bfb[mc:mc + nd]
                    dsb = bfb[ms:ms + nd]
                    bff[mc:mc + nd, maci] = -cmu[m - 1] * dcb
                    bff[ms:ms + nd, maci] = smu[m - 1] * dcb
                    bff[mc:mc + nd, masi] = -smu[m - 1] * dsb
                    bff[ms:ms + nd, masi] = -cmu[m - 1] * dsb
                cj = np.zeros((nd, mphi2))
                sj = np.zeros((nd, mphi2))
                cj[ia - 1, ma - 1] = 1.0
                sj[ia - 1, ma - 1] = 1.0
                if ma == 2:
                    m = 2
                    diff, diff3 = _diff2(d, eps, aj, m)
                    cj[ia - 1, m - 1] /= (1.0 + (diff + diff3))
                    sj[ia - 1, m - 1] /= (1.0 + (diff - diff3))
                elif ma > 2:
                    m = ma
                    diff = 1.0 + d * eps * (m - 1) ** 2 * 0.125 / aj
                    cj[ia - 1, m - 1] /= diff
                    sj[ia - 1, m - 1] /= diff
                if ma <= mphi2 - 2:
                    m = ma
                    m2 = m + 2
                    diff2 = -d * (m + 1) * (-1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
                    cj[ia - 1, m2 - 1] -= diff2 * cj[ia - 1, m - 1]
                    sj[ia - 1, m2 - 1] -= diff2 * sj[ia - 1, m - 1]
                    if m2 > 3:
                        diff1 = -d * (m2 - 3) * (1.0 + eps / aj * 0.5 * (m2 - 1)) * 0.25
                        cj[ia - 1, ma - 1] -= diff1 * cj[ia - 1, m2 - 1]
                        sj[ia - 1, ma - 1] -= diff1 * sj[ia - 1, m2 - 1]
                if ma > 3:
                    m = ma
                    m1 = m - 2
                    diff1 = -d * (m - 3) * (1.0 + eps / aj * 0.5 * (m - 1)) * 0.25
                    cj[ia - 1, m1 - 1] -= diff1 * cj[ia - 1, m - 1]
                    sj[ia - 1, m1 - 1] -= diff1 * sj[ia - 1, m - 1]
                for m in range(ma - 2, ma + 3, 2):
                    if m == 2:
                        diff, diff3 = _diff2(d, eps, aj, m)
                        cj[ia - 1, m - 1] /= (1.0 + (diff + diff3))
                        sj[ia - 1, m - 1] /= (1.0 + (diff - diff3))
                    elif 3 <= m <= mphi2:
                        diff = 1.0 + d * eps * (m - 1) ** 2 * 0.125 / aj
                        cj[ia - 1, m - 1] /= diff
                        sj[ia - 1, m - 1] /= diff
                bff[0:nsb, maci] += cj.T.reshape(-1)
                bff[nsb:2 * nsb, masi] += sj.T.reshape(-1)
        if nd == 10:
            bfb = tesetdm(j, bd, hc, hs, mphi, mphi2)
        else:
            bfb = np.zeros(2 * nsb)
        bfb[:nsb] += bb[:, :, j - 1].T.reshape(-1)
        for m in range(1, mphi2 + 1):
            mc = (m - 1) * nd
            ms = mc + nsb
            dcb = bfb[mc:mc + nd].copy()
            dsb = bfb[ms:ms + nd].copy()
            bfb[mc:mc + nd] = cmu[m - 1] * dcb + smu[m - 1] * dsb
            bfb[ms:ms + nd] = -smu[m - 1] * dcb + cmu[m - 1] * dsb
        bfx = tsolva(bff, bfb)
        bfx = np.where(np.isfinite(bfx), bfx, 0.0)             # resetnan
        bc[:, :, j - 1] = bfx[:nsb].reshape(mphi2, nd).T
        bs[:, :, j - 1] = bfx[nsb:2 * nsb].reshape(mphi2, nd).T


def tesumb(bc, hc, mphi2, ndp, w, dj, sige, dispp, fj):
    """tesumb: 平衡のビーム行列 beam(21) (packed 6×6)"""
    e = sige ** 2
    beam = np.zeros(21)
    for j in range(1, ndp + 1):
        f = fj[j - 1] * dj / e / w
        beam[0:10] += bc[:, 0, j - 1] * f
    for j in range(1, ndp + 1):
        aj = (j - 0.5) * dj
        f = fj[j - 1] * dj / e / w * math.sqrt(2.0 * aj)
        beam[15:19] += hc[:, 1, j - 1] * f
    for i in range(1, 5):
        for i1 in range(1, i + 1):
            ii = iaidx(i, i1) - 1
            beam[ii] += dispp[i - 1] * dispp[i1 - 1] * e
    beam[20] = e
    return beam


def beam6(beam21):
    """packed (iaidx) 21 要素 → 6×6 対称行列"""
    M = np.zeros((6, 6))
    for m in range(1, 7):
        for n in range(1, m + 1):
            M[m - 1, n - 1] = M[n - 1, m - 1] = beam21[iaidx(m, n) - 1]
    return M


def physdisp(tw):
    """tgfun.f tgetphysdispu: tws の 1 点から物理の分散 (4 成分)"""
    r1, r2, r3, r4 = tw[10], tw[11], tw[12], tw[13]
    detr = r1 * r4 - r2 * r3
    cc = math.sqrt(1.0 - detr)
    ex, epx, ey, epy = tw[6], tw[7], tw[8], tw[9]
    if tw[14] < 1.0:
        return np.array([cc * ex + r4 * ey - r2 * epy,
                         cc * epx - r3 * ey + r1 * epy,
                         cc * ey - r1 * ex - r2 * epx,
                         cc * epy - r3 * ex - r4 * epx])
    return np.array([cc * ey - r1 * ex - r2 * epx,
                     cc * epy - r3 * ex - r4 * epx,
                     cc * ex + r4 * ey - r2 * epy,
                     cc * epx - r3 * ey + r1 * epy])


def _tsolva(a, b, eps=1.0e-20):
    """tsolva = tsvd(.false.) の最小ノルム解 (特異値 < 最大 × eps を捨てる)"""
    U, sv, Vt = np.linalg.svd(a, full_matrices=False)
    thr = sv.max() * eps if sv.size else 0.0
    sinv = np.where(sv > thr, 1.0 / np.where(sv > thr, sv, 1.0), 0.0)
    return Vt.T @ (sinv * (U.T @ b))


def synchro_beta_emittance(lat, momentum, amus0, amus1, amusstep, mphi2, verbose=False):
    """tfemits → temits → temits1。amus0, amus1, amusstep は [rad] (SAD の引数 × 2π)。戻り: [[nus, emix, emiy, emixp, emiyp, conv], …]"""
    from .emit6d import emittance6d, tturne, _tcod_radcod, _Phys, _tinv6, ELMASS
    from .calc6d import tfetwiss, etwiss2ri
    mphi2 = int(max(1, min(32, mphi2)))
    mphiz = mphi2 * 2 - 1
    ndims = mphiz
    ndps = max(15, mphiz)
    nzz = mphiz * 2 + 1
    mphi = mphiz
    ndp = ndps
    nz = nzz
    itmax = 20
    ndiffh = 16
    resconv = 1.0e-4
    minconv = 1.0e-7

    els = lat.elements
    circ = sum(el.L for el in els)
    ph = _Phys(momentum, circ, amass=ELMASS, charge=1.0)
    # temit (RFSW・RADCOD・RADTAPER・TRPT・INTRA は temits が OFF にする)
    out = emittance6d(lat, momentum, rfsw=False, radcod=False)
    r = np.asarray(out['normal_r'], dtype=float)
    ri = np.asarray(out['normal_ri'], dtype=float)
    sige = float(out['momentum_spread'])
    damp = float(out['damping_rate'][2])
    alphap = float(out['momentum_compaction'])
    heff = float(out['heff'])
    omega0 = float(out['omega0'])
    trf0 = float(out['trf0'])
    dz = float(out.get('dz', 0.0) or 0.0)
    codin = np.array(out['orbit_entrance'], dtype=float)
    ext = np.array(out['ext_twiss'], dtype=float)
    dispp = r[:, 5].copy()
    rx = np.eye(6)
    rxi = np.eye(6)
    rx[0:4, 5] = dispp[0:4]
    rxi[0:4, 5] = -dispp[0:4]
    conv = max(resconv * abs(damp), minconv)
    sigea = math.sqrt(12.0) * sige
    dpndim = sigea / ndims
    tw0 = ext.copy()
    tw0[2] = 0.0
    tw0[5] = 0.0
    tw0[23] = 0.0
    helpers = (etwiss2ri, _tinv6, None)

    beams = np.zeros((2 * ndims + 1, 10))
    trads = np.zeros((2 * ndims + 1, 5, 5))
    tws = np.zeros((2 * ndims + 1, 28))
    i = 0
    for _ in range(2 * ndims + 1):
        if i == 0:
            cod = codin.copy()
        elif i > 0:
            disppi = np.zeros(6)
            disppi[0:4] = physdisp(tws[i - 1 + ndims])
            disppi[4] = 0.0
            disppi[5] = 1.0
            cod = tws[i - 1 + ndims][15:21] + disppi * dpndim
        else:
            disppi = np.zeros(6)
            disppi[0:4] = physdisp(tws[i + 1 + ndims])
            disppi[4] = 0.0
            disppi[5] = 1.0
            cod = tws[i + 1 + ndims][15:21] - disppi * dpndim
        cod = np.array(cod, dtype=float)
        cod[4] = 0.0
        cod0 = cod.copy()
        cod_found, fnd, _trf, _st, _stf = _tcod_radcod(els, ph, cod, False, radcod=False)
        if fnd:
            cod = np.array(cod_found, dtype=float)
        if not fnd and (i > 1 or i < -1):
            if i > 1:
                cod[0:4] = 2.0 * tws[i - 1 + ndims][15:19] - tws[i - 2 + ndims][15:19]
            else:
                cod[0:4] = 2.0 * tws[i + 1 + ndims][15:19] - tws[i + 2 + ndims][15:19]
            cod[4] = 0.0
            cod_found, fnd, _trf, _st, _stf = _tcod_radcod(els, ph, cod, False, radcod=False)
            if fnd:
                cod = np.array(cod_found, dtype=float)
            else:
                cod = cod0.copy()
        st = tturne(els, cod, ph, trf0=trf0, rfsw=False, radcod=False)
        codx = np.array(st.cod, dtype=float)
        if fnd:
            codx[0:4] = codx[0:4] - codin[0:4]
        else:
            codx[0:4] = cod0[0:4] - codin[0:4]
        rm = rxi @ (st.D @ rx)
        trads[i + ndims] = tmov65(rm)
        tws[i + ndims] = tfetwiss(ri @ _tinv6(st.T), codx, True)
        Bn = rxi @ st.B @ rxi.T
        beams[i + ndims] = pack10(Bn[0:4, 0:4])
        if i <= 0:
            i = -i + 1
        else:
            i = -i
    ns = 1 if amusstep == 0.0 else int((amus1 - amus0) / amusstep + 1.1)
    dj = sigea ** 2 / ndp * 0.5
    w = 0.0
    fj = np.zeros(ndp)
    for j in range(1, ndp + 1):
        fj[j - 1] = math.exp(-(j - 0.5) * dj / sige ** 2)
        w += fj[j - 1] * dj / sige ** 2
    ha, hb, ba, bb, bd = tefsetup(ndp, mphi, mphi2, nz, ndims, beams, trads, tws, dj, dpndim, sigea, tw0, ri, helpers)
    phi0s = heff * omega0 * (trf0 + dz) / C_LIGHT
    fz = math.cos(phi0s) / (1.0 + math.cos(phi0s)) * 0.25 * (heff * alphap * PI2) ** 2
    dispp0 = physdisp(tws[ndims])
    results = []
    for kk in range(1, ns + 1):
        amus = amus0 + amusstep * (kk - 1)
        amuj = np.array([amus - fz * (j - 0.5) * dj / amus for j in range(1, ndp + 1)])
        hc = np.zeros((4, mphi2, ndp))
        hs = np.zeros((4, mphi2, ndp))
        bc = np.zeros((10, mphi2, ndp))
        bs = np.zeros((10, mphi2, ndp))
        tesolvd(hb, hc, hs, ha, amuj, bd, hc, hs, ndp, mphi2, mphi, 4, dj, damp, sige, _tsolva)
        nd1 = max(1, min(int(abs(0.05 / damp)), ndiffh))
        for _ in range(nd1):
            tevdif(hb, hc, hs, ha, amuj, bd, hc, hs, ndp, mphi2, mphi, 4, dj, damp, sige)
        tesolvd(bb, bc, bs, ba, amuj, bd, hc, hs, ndp, mphi2, mphi, 10, dj, damp, sige, _tsolva)
        beam = tesumb(bc, hc, mphi2, ndp, w, dj, sige, dispp0, fj)
        beamr = ri @ beam6(beam) @ ri.T
        vx = beamr[0, 0] * beamr[1, 1] - beamr[0, 1] ** 2
        vy = beamr[2, 2] * beamr[3, 3] - beamr[2, 3] ** 2
        emix0 = math.copysign(math.sqrt(abs(vx)), vx)
        emiy0 = math.copysign(math.sqrt(abs(vy)), vy)
        res0 = 1.0e100
        res = 1.0
        it = 0
        while res > conv:
            tevdif(hb, hc, hs, ha, amuj, bd, hc, hs, ndp, mphi2, mphi, 4, dj, damp, sige)
            tevdif(bb, bc, bs, ba, amuj, bd, hc, hs, ndp, mphi2, mphi, 10, dj, damp, sige)
            beam = tesumb(bc, hc, mphi2, ndp, w, dj, sige, dispp0, fj)
            beamr = ri @ beam6(beam) @ ri.T
            vx = beamr[0, 0] * beamr[1, 1] - beamr[0, 1] ** 2
            vy = beamr[2, 2] * beamr[3, 3] - beamr[2, 3] ** 2
            emix = math.copysign(math.sqrt(abs(vx)), vx)
            emiy = math.copysign(math.sqrt(abs(vy)), vy)
            res = ((emix - emix0) ** 2 + (emiy - emiy0) ** 2) / (emix + emiy) ** 2
            emix0 = emix
            emiy0 = emiy
            it += 1
            if it > itmax or (it > 1 and (res > resconv or res > res0)):
                res0 = -res
                res = 0.0
            else:
                res0 = res
        B6 = beam6(beam)
        vx = B6[0, 0] * B6[1, 1] - B6[0, 1] ** 2
        vy = B6[2, 2] * B6[3, 3] - B6[2, 3] ** 2
        emixp = math.copysign(math.sqrt(abs(vx)), vx)
        emiyp = math.copysign(math.sqrt(abs(vy)), vy)
        results.append([amus / PI2, emix, emiy, emixp, emiyp, res0])
    return results
