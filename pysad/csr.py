"""
CSRTrack / CSRHaissin (csrtrack.f) の直訳: CSR ウェイクによる縦方向のマクロ粒子追跡と Haissinski 方程式の解。
Packages/CSR.n (Track、InitializeParticles) が呼ぶ。CSRInit / CSRMatrix / CSROy* は別 (docs/36)。

実FFT は K64 の ftr モジュール (tcftr / trftr / tftrr) の配置: a(0:m-1) に (Re f0, Re f(m/2)) と (Re, Im) の組を詰める。
乱数 tgauss_array は SadRandom.gauss(np, GCUT) (ffs の乱数と共通)。
CSRTrack は z・dp のリストをその場で更新する (K64 はスタック上のリスト本体を直接書き換える)。
"""
from __future__ import annotations

import math

import numpy as np

CS = 299792458.0
DPMAX = 0.1


# ---------------------------------------------------------------------------
# tcftr.f: 複素 FFT と実 FFT (ftr モジュール)
# ---------------------------------------------------------------------------

def _tcftr(ca, conj):
    """tcftr: conj=False は exp(+2πi jk/n)、conj=True は exp(−2πi jk/n) (どちらも規格化なし)"""
    n = len(ca)
    return np.fft.fft(ca) if conj else np.fft.ifft(ca) * n


def trftr(a, conj):
    """trftr(a,m,conj): 実数列 a(0:m-1) の実FFT (K64 の詰め方)。a を書き換えて返す"""
    m = len(a)
    ca = a[0::2] + 1j * a[1::2]
    ca = _tcftr(ca, conj)
    a[0::2] = ca.real
    a[1::2] = ca.imag
    w = (-2.0 if conj else 2.0) * math.pi / m
    i = np.arange(1, m // 4 + 1)
    ci = np.cos(w * i)
    si = np.sin(w * i)
    i2 = i * 2
    j2 = m - i2
    ai2 = a[i2].copy()
    ai21 = a[i2 + 1].copy()
    aj2 = a[j2].copy()
    aj21 = a[j2 + 1].copy()
    a[i2] = 0.5 * (ai2 + aj2 + ci * (ai21 + aj21) + si * (ai2 - aj2))
    a[j2] = 0.5 * (ai2 + aj2 - ci * (ai21 + aj21) - si * (ai2 - aj2))
    a[i2 + 1] = 0.5 * (ai21 - aj21 + si * (ai21 + aj21) - ci * (ai2 - aj2))
    a[j2 + 1] = 0.5 * (aj21 - ai21 + si * (ai21 + aj21) - ci * (ai2 - aj2))
    a1 = a[1]
    a[1] = a[0] - a1
    a[0] = a[0] + a1
    return a


def tftrr(a, conj):
    """tftrr(a,m,conj): trftr の逆 (tftrr(trftr(f)) = m f)。a を書き換えて返す"""
    m = len(a)
    w = (2.0 if conj else -2.0) * math.pi / m
    i = np.arange(1, m // 4 + 1)
    ci = np.cos(w * i)
    si = np.sin(w * i)
    i2 = i * 2
    j2 = m - i2
    ai2 = a[i2].copy()
    ai21 = a[i2 + 1].copy()
    aj2 = a[j2].copy()
    aj21 = a[j2 + 1].copy()
    a[i2] = ai2 + aj2 - ci * (ai21 + aj21) + si * (ai2 - aj2)
    a[j2] = ai2 + aj2 + ci * (ai21 + aj21) - si * (ai2 - aj2)
    a[i2 + 1] = ai21 - aj21 + si * (ai21 + aj21) + ci * (ai2 - aj2)
    a[j2 + 1] = aj21 - ai21 + si * (ai21 + aj21) + ci * (ai2 - aj2)
    a1 = a[1]
    a[1] = a[0] - a1
    a[0] = a[0] + a1
    ca = a[0::2] + 1j * a[1::2]
    ca = _tcftr(ca, conj)
    a[0::2] = ca.real
    a[1::2] = ca.imag
    return a


# ---------------------------------------------------------------------------
# spline.f spline1 (dx=1 の自然スプライン。ddy = y''/6)
# ---------------------------------------------------------------------------

def spline1(y):
    """spline1(np,y,ddy,work,0,0): y(1..np) → ddy(1..np) (mode1=mode2=0)。配列は 0 始まり (y[0] = y(1))"""
    n = len(y)
    ddy = np.zeros(n)
    if n <= 2:
        return ddy
    work = np.zeros(n)
    ddy[1] = y[2] - 2.0 * y[1] + y[0]
    work[1] = 4.0
    for i in range(1, n - 2):                      # Fortran i=2..np-2
        work[i + 1] = 4.0 - 1.0 / work[i]
        ddy[i + 1] = y[i + 2] - 2.0 * y[i + 1] + y[i] - ddy[i] / work[i]
    ddy[n - 1] = 0.0
    ddy[n - 2] = ddy[n - 2] / work[n - 2]
    for i in range(n - 3, 0, -1):                  # Fortran i=np-2..2
        ddy[i] = (ddy[i] - ddy[i + 1]) / work[i]
    return ddy


# ---------------------------------------------------------------------------
# csrwake / csrrho
# ---------------------------------------------------------------------------

def csrwake(rho, zl):
    """rho ← IFFT(zl · FFT(rho)) / nr (rho はその場で更新)"""
    nr = len(rho)
    trftr(rho, True)
    i = np.arange(2, nr // 2 + 1)
    i2 = i * 2 - 1                                 # Fortran の 1 始まり → 0 始まりの添字 i2-1
    re = rho[i2 - 1].copy()
    im = rho[i2].copy()
    zr = zl[i2 - 1]
    zi = zl[i2]
    rho[i2 - 1] = zr * re - zi * im
    rho[i2] = zr * im + zi * re
    rho[0] = zl[0] * rho[0]
    rho[1] = zl[1] * rho[1]
    tftrr(rho, False)
    rho /= nr
    return rho


def csrrho(z, dz, nr):
    ic = nr // 4
    zc = dz * ic
    rmin = 1.0
    rmax = nr / 2 - 1e-10
    rho = np.zeros(nr)
    r = (np.asarray(z) + zc) / dz
    ok = (r >= rmin) & (r <= rmax)
    r = r[ok]
    l = np.floor(r).astype(int)
    a = r - l
    np.add.at(rho, l - 1, 1.0 - a)                 # Fortran rho(l) → 0 始まり rho[l-1]
    np.add.at(rho, l, a)
    return rho / dz


# ---------------------------------------------------------------------------
# csrtrack
# ---------------------------------------------------------------------------

def csrtrack(z, dp, zl, dz, ndiv, aw, r56, r65, frf, phirf, damp, sigmae, nt, gauss):
    """z, dp: np.ndarray (その場で更新)。zl: 実FFT 形式のインピーダンス。gauss(np) は乱数 (tgauss_array)。
    戻り: data (2, nt) = {σz, σδ} の各周回 (発散したらそれ以降は 0)"""
    np_ = len(z)
    nr = len(zl)
    awn = aw / ndiv
    ad = 1.0 - 2.0 * damp / ndiv
    bd = math.sqrt(4.0 * damp / ndiv * (1.0 - damp / ndiv)) * sigmae
    zc = dz * nr / 4.0
    nr1 = nr // 2
    r56n = r56 / ndiv
    if frf != 0.0:
        wr = 2.0 * math.pi * frf / CS
        vr = r65 / wr
        dpr = vr * math.sin(phirf)
        zlim = 0.5 * CS / frf
    else:
        wr = vr = dpr = 0.0
        zlim = 1e100
    nc = (ndiv - 1) // 2 + 1
    data = np.zeros((2, nt))
    for i in range(nt):
        for n in range(1, ndiv + 1):
            z += r56n * dp
            rho = csrrho(z, dz, nr)
            csrwake(rho, zl)
            ddrho = spline1(rho[:nr1])
            rand = np.asarray(gauss(np_), dtype=float)
            zk = (z + zc) / dz
            m = zk.astype(np.int64)                # int() は 0 方向への切り捨て
            ins = (m > 0) & (m < nr1)
            mi = np.where(ins, m, 1)
            b = zk - mi
            a = 1.0 - b
            c = -a * b * (a + 1.0)
            d = -a * b * (b + 1.0)
            w = a * rho[mi - 1] + b * rho[mi] + c * ddrho[mi - 1] + d * ddrho[mi]
            kick = np.where(ins, awn * w, 0.0)
            if n == nc:
                if frf == 0.0:
                    dp[:] = (dp + kick + r65 * z) * ad + bd * rand
                else:
                    dp[:] = (dp + kick + vr * np.sin(wr * z + phirf) - dpr) * ad + bd * rand
            else:
                dp[:] = (dp + kick) * ad + bd * rand
        if frf == 0.0:
            sz = z.sum() / np_
            sdp = dp.sum() / np_
            dp -= sdp
            z -= sz
            vz = float((z ** 2).sum())
            vdp = float((dp ** 2).sum())
        else:
            lost = dp < -DPMAX
            z[lost] = -zlim * 10.0
            inb = (~lost) & (z > -zlim) & (z < zlim)
            sz = z[inb].sum() / np_
            sdp = dp[inb].sum() / np_
            inb2 = (z > -zlim) & (z < zlim)
            vz = float(((z[inb2] - sz) ** 2).sum())
            vdp = float(((dp[inb2] - sdp) ** 2).sum())
        data[0, i] = math.sqrt(vz / np_)
        data[1, i] = math.sqrt(vdp / np_)
        if data[1, i] > DPMAX or data[0, i] > zc:
            data[:, i + 1:] = 0.0
            return data
    return data


# ---------------------------------------------------------------------------
# csrhaissin
# ---------------------------------------------------------------------------

def csrhaissin(zl, aw, dz, sigz):
    """Haissinski 方程式を Newton 法で解く。戻り: rho (nr)。収束しなければ None (K64 は irtc=-1)"""
    from .tsvd import tsvd
    nr = len(zl)
    imax = 100
    thre = 1e-5
    pi = math.pi
    fact = 1.0
    sdf0 = 1e100
    nr1 = nr // 2
    zc = (nr // 4) * dz
    a = dz / sigz / math.sqrt(2.0 * pi)
    u = np.zeros(nr1 + 1)
    u0 = np.zeros(nr1 + 1)
    u[0] = math.log(a)
    k = np.arange(1, nr1 + 1)
    u0[1:] = -((dz * k - zc) / sigz) ** 2 / 2.0
    u[1:] = u[0] + u0[1:]
    rho = np.zeros(nr)
    rho[:nr1] = np.exp(u[1:])                       # rho(i) i=1..nr1 → rho[i-1]
    dk = pi / dz / nr
    zi = np.zeros(nr)
    for i in range(3, nr, 2):                       # Fortran i=3,nr-1,2 (1 始まり)
        ak = dk * (i - 1)
        zi[i - 1] = zl[i] / ak                      # zi(i) = zl(i+1)/ak
        zi[i] = -zl[i - 1] / ak                     # zi(i+1) = -zl(i)/ak
    wi = zi / nr / dz
    tftrr(wi, False)
    du = np.zeros(nr1 + 1)
    ulast = u.copy()
    iter_ = 0
    for _ in range(imax):
        v = rho / dz
        csrwake(v, zi)
        f = np.zeros(nr1 + 1)
        g = np.zeros((nr1 + 1, nr1 + 1))
        f[1:] = u[1:] - u[0] - u0[1:] + aw * v[:nr1]
        sdf = float(np.abs(f[1:]).sum())
        # g(k,m) = aw*wi(k-m+1)*rho(m) (m<=k)、aw*wi(nr-(m-k)+1)*rho(m) (m>k)。wi, rho は 1 始まり
        kk = np.arange(1, nr1 + 1)[:, None]
        mm = np.arange(1, nr1 + 1)[None, :]
        idx = np.where(mm <= kk, kk - mm + 1, nr - (mm - kk) + 1)
        g[1:, 1:] = aw * wi[idx - 1] * rho[:nr1][None, :]
        g[kk[:, 0], kk[:, 0]] += 1.0
        g[1:, 0] = -1.0
        g[0, 1:] = rho[:nr1]
        if sdf > sdf0:
            fact *= 0.5
        else:
            if sdf < thre:
                return rho / dz
            elif iter_ > imax:
                return None
            sdf0 = sdf
            fact = min(1.0, fact * 2.0)
            f[0] = 0.0
            g[0, 0] = 0.0
            du = np.asarray(tsvd(g, f, 1e-8))
            ulast = u.copy()
        iter_ += 1
        un = ulast[1:] - du[1:] * fact
        rn = np.exp(un)
        sr = float(rn.sum())
        rho[:nr1] = rn
        rho[nr1:] = 0.0
        du0 = -math.log(sr)
        u[0] = ulast[0] - du[0] * fact + du0
        rho[:nr1] = rho[:nr1] / sr
        u[1:] = un + du0
    return None


# ---------------------------------------------------------------------------
# csroy.f: Oide-Yokoya の CSR (作用角 J の格子で線形化したウェイク行列)
# ---------------------------------------------------------------------------

def csrwkick(nrh, z, zc, dz, v, ddv, aw, r65):
    """csrwkick: ウェイクとRFの 1 回分の dp のキック (z はスカラー)"""
    zk = (z + zc) / dz
    m = int(zk)
    if 0 < m < nrh:
        b = zk - m
        a = 1.0 - b
        c = -a * b * (a + 1.0)
        d = -a * b * (b + 1.0)
        w = a * v[m - 1] + b * v[m] + c * ddv[m - 1] + d * ddv[m]
        return aw * w + r65 * z
    return r65 * z


def csrgetj(nrh, z1, zc, dz, v, ddv, aw, r56, r65, dt):
    """csrgetj: 位相空間の周回軌道 (z1 から出発) の作用 J と角周波数 omega"""
    z = z1
    p = 0.5 * dt * csrwkick(nrh, z, zc, dz, v, ddv, aw, r65)
    p0 = p
    aj = 0.0
    n = 0
    f = 1.0
    while p * p0 > 0.0:
        z = z + p * r56 * dt
        f = csrwkick(nrh, z, zc, dz, v, ddv, aw, r65)
        aj = aj + p * r56 * (p + dt * f * 0.5)
        p = p + dt * f
        n += 1
    t0 = 2.0 * (n * dt + 0.5 * dt - p / f)
    omega = 2.0 * math.pi / t0
    aj = abs(aj * dt / 2.0 / math.pi)
    return aj, omega


def csroysetup(rho, zl, aw, r56, r65, dz, nj, mphi, zspan, dt):
    """戻り: (wi, aj, omega, zjn(nj,mphi), dpjn(nj,mphi), rhoj, irtc)。
    K64 は dpjn を zjn と同じ記憶域に割り当てる (tfcsroysetup: c_loc(rlist(kazjn)) を 2 回) ので、zjn(k,i) の後に dpjn(k,i) を書くと
    zjn が dp の値で上書きされる (docs/21 B59)。返す zjn は dp の値、dpjn の記憶域は書かれない (0)"""
    nr = len(rho)
    nrh = nr // 2
    dk = math.pi / dz / nr
    zi = np.zeros(nr)
    for i in range(3, nr, 2):
        ak = dk * (i - 1)
        zi[i - 1] = zl[i] / ak
        zi[i] = -zl[i - 1] / ak
    wi = zi / nr / dz
    tftrr(wi, False)
    v = np.array(rho, dtype=float)
    csrwake(v, np.asarray(zl, dtype=float))
    ddv = spline1(v[:nrh])
    zc = 0.0
    s = 0.0
    for i in range(1, nrh + 1):
        zc += rho[i - 1] * csrwkick(nrh, i * dz, 0.0, dz, v, ddv, aw, r65)
        s += rho[i - 1]
    zc = zc / s / r65
    rhoth = math.exp(-zspan ** 2 / 2.0) / math.sqrt(2.0 * math.pi) / dz
    i1 = 0
    for i in range(1, nrh + 1):
        if rho[i - 1] >= rhoth:
            i1 = i
            break
    f1 = csrwkick(nrh, dz - zc, zc, dz, v, ddv, aw, r65)
    ajz = np.zeros(nrh)
    for i in range(1, nrh + 1):
        z = dz * i - zc
        f = csrwkick(nrh, z, zc, dz, v, ddv, aw, r65)
        if f * f1 >= 0.0:
            ajz[i - 1], _ = csrgetj(nrh, z, zc, dz, v, ddv, aw, r56, r65, dt)
    aj = np.zeros(nj)
    omega = np.zeros(nj)
    zjn = np.zeros((nj, mphi))        # [j, k] = zjn(k+1, j+1)
    dpjn = np.zeros((nj, mphi))
    from . import compat as _compat
    _v2 = _compat.v2()
    if _v2 and i1 == 0:
        # v2: しきい値 rhoth 以上の点が無い (K64 の csrjasign は i1 が未初期化のまま進む)。失敗として返す
        return wi, np.zeros(nj), np.zeros(nj), zjn, dpjn, np.zeros(nj), -1
    rhoj = np.zeros(nj)
    # csrjasign
    dj = ajz[i1 - 1] / nj if i1 > 0 else 0.0
    k = 0
    m = 1
    m0 = 0
    done = False
    while m0 != m:
        ajn = dj * (k + 0.5)
        m0 = m
        for i in range(1, nrh):
            if ajz[i - 1] > ajn:
                if ajz[i] > 0.0 and ajz[i] <= ajn:
                    zjn[m - 1, 0] = i * dz - zc + (ajz[i - 1] - ajn) / (ajz[i - 1] - ajz[i]) * dz
                    aj[m - 1], omega[m - 1] = csrgetj(nrh, zjn[m - 1, 0], zc, dz, v, ddv, aw, r56, r65, dt)
                    m += 1
                    if m > nj:
                        done = True
                        break
        if done:
            break
        k += 1
    if not done:
        return wi, aj, omega, zjn, np.zeros((nj, mphi)), rhoj, -1
    srho = 0.0
    for i in range(nj):
        dt1 = math.pi / omega[i] / mphi
        z = zjn[i, 0]
        p = 0.5 * dt1 * csrwkick(nrh, z, zc, dz, v, ddv, aw, r65)
        zi_ = (z + zc) / dz
        iz = int(zi_)
        rhoj[i] = rho[iz - 1] + (zi_ - iz) * (rho[iz] - rho[iz - 1])
        srho += rhoj[i]
        for kk in range(mphi):
            if _v2:
                zjn[i, kk] = z             # v2 (B59 の修正): zjn は位置、dpjn は運動量 (別の記憶域)
                dpjn[i, kk] = p
            else:
                zjn[i, kk] = p             # zjn(k,i)=z の後の dpjn(k,i)=p が同じ記憶域を上書きする (B59)
            z = z + p * r56 * dt1
            f = csrwkick(nrh, z, zc, dz, v, ddv, aw, r65)
            p = p + dt1 * f
    srho = 1.0 / 2.0 / math.pi / srho
    rhoj = np.sqrt(rhoj * srho)
    return wi, aj, omega, zjn, dpjn, rhoj, 0


def csroymat(wi0, omega, zjn, rhoj, bw, dz, nj, nphi, nr):
    """csroymat: bm(nj·nphi, nj·nphi)。zjn は (nj, mphi) の配列 (K64 の zjn(mphi,nj))"""
    mphi = nphi * 16
    dphi = math.pi / mphi
    nrh = nr // 2
    wi = np.empty(nr)
    wi[:nrh] = bw * wi0[nrh:nr]
    wi[nrh:] = bw * wi0[:nrh]
    ddwi = spline1(wi)
    n = np.arange(1, nphi + 1)
    m = np.arange(1, mphi + 1)
    x = n[None, :] * (m[:, None] - 0.5) / mphi            # (m, n)
    cosmn = -n[None, :] * np.cos(np.mod(x, 2.0) * math.pi)
    wj = rhoj * omega * dphi                                # (j)
    Z = zjn.reshape(nj * mphi)                              # (j,m) の順
    zd = (Z[:, None] - Z[None, :]) / dz + nrh + 1.0         # [(j,m),(j1,m1)]
    mz = zd.astype(np.int64)
    b = zd - mz
    a = 1.0 - b
    c = -a * b * (a + 1.0)
    d = -a * b * (b + 1.0)
    F = a * wi[mz - 1] + b * wi[mz] + c * ddwi[mz - 1] + d * ddwi[mz]
    F = F.reshape(nj, mphi, nj, mphi)                       # [j,m,j1,m1]
    # bm[(j1,n1),(j,n)] = Σ_{m,m1} wj[j] wj[j1] cosmn[m,n] F[j,m,j1,m1] cosmn[m1,n1]
    t = np.einsum('jmkl,mn->jnkl', F, cosmn)                # [j,n,j1,m1]
    t = np.einsum('jnkl,lp->jnkp', t, cosmn)                # [j,n,j1,n1]
    t = t * wj[:, None, None, None] * wj[None, None, :, None]
    bm = t.transpose(2, 3, 0, 1).reshape(nj * nphi, nj * nphi)   # 行 = (j1,n1)、列 = (j,n)
    bm = np.array(bm)
    for j in range(nj):
        for nn in range(1, nphi + 1):
            jn = nphi * j + nn - 1
            bm[jn, jn] += (nn * omega[j]) ** 2
    return bm


# ---------------------------------------------------------------------------
# csrtest.f: 矩形ダクト断面の格子 (CSRConvert / CSRMatrix の差分行列)
# ---------------------------------------------------------------------------

def csrind(pipe):
    """csrind: ind(0:m+1, 0:n+1) (行 i = x 方向の列、j = y 方向) と点の総数 mm"""
    m = len(pipe)
    n = int(max(pipe)) if m else 0
    ind = np.zeros((m + 2, n + 2), dtype=np.int64)
    k = 0
    for i in range(1, m + 1):
        for j in range(1, int(pipe[i - 1]) + 1):
            k += 1
            ind[i, j] = k
    return ind, k, n


def csrmat(ind, dx, dy, omega, rho, fudge0, fudge1, fudge2, fudge3, tilt, mode, mm, m, n, ic):
    """csrmat: 差分行列 acm(mm,mm)。mode 110/220/330 (抵抗壁) は第 1 列だけ使う"""
    acm = np.zeros((mm, mm))
    if mode == 110:
        pass
    elif mode == 220:
        ac = 0.5 / omega
        cfy = ac * 2.0 / dy ** 3
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                jm = ind[i, j]
                if jm == 0:
                    break
                if ind[i, j + 1] == 0:
                    acm[jm - 1, 0] = cfy
    elif mode == 330:
        ac = 0.5 / omega
        cfx = ac * 2.0 / dx ** 3
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                jm = ind[i, j]
                if jm == 0:
                    break
                if ind[i - 1, j] == 0:
                    acm[jm - 1, 0] = cfx
                if ind[i + 1, j] == 0:
                    acm[jm - 1, 0] += cfx
    elif mode in (22, 33):
        orr = omega * rho
        for i in range(1, m + 1):
            dx1 = dx - (i - ic - 0.5) * tilt / m * dx
            dy1 = dy + (i - ic - 0.5) * tilt / m * dy
            if rho == 0.0:
                rdx = 1.0 / dx1
                rdy = 1.0 / dy1
                ce = 0.0
                cg = 0.0
                ac = 0.5 / omega
            else:
                x = (i - ic - 0.5) * dx
                r = rho + x
                rdx = r / dx1
                rdy = r / dy1
                if mode == 22:
                    ac = 0.5 * omega / (orr ** 2 - 1.0)
                    ce = (3.0 + omega ** 2 * x * (r + rho)) * ac
                    cg = 1.5 * rdx * ac
                else:
                    ac = 0.5 / orr / rho
                    ce = omega ** 2 * x * (r + rho) * ac
                    cg = 0.5 * rdx * ac
            cfx = rdx ** 2 * ac
            cfy = rdy ** 2 * ac
            if mode == 22:
                cfy1, cfya, cfxa1, cfxa2 = cfy, -cfy, cfx - cg, cfx + cg
            else:
                cfy1, cfya, cfxa1, cfxa2 = -cfy, cfy, -cfx + cg, -cfx - cg
            for j in range(1, n + 1):
                jm = ind[i, j]
                if jm == 0:
                    break
                j1 = ind[i, j + 1]
                im1 = ind[i - 1, j]
                i1 = ind[i + 1, j]
                if (im1 == 0 or i1 == 0) and j1 == 0:
                    fc = fudge0
                elif j1 != 0 and im1 != 0 and i1 != 0 and (ind[i + 1, j + 1] == 0 or ind[i - 1, j + 1] == 0):
                    fc = fudge2
                else:
                    fc = 1.0
                acm[jm - 1, jm - 1] = ce - 2.0 * (cfx + cfy) * fc
                if j == 1:
                    acm[jm - 1, jm - 1] += cfy1
                    acm[jm - 1, j1 - 1] = cfy
                else:
                    jm1 = ind[i, j - 1]
                    acm[jm - 1, jm1 - 1] = cfy
                    if j1 != 0:
                        acm[jm - 1, j1 - 1] = cfy
                    else:
                        acm[jm - 1, jm - 1] += cfya * fc
                if im1 != 0:
                    acm[jm - 1, im1 - 1] = cfx - cg
                else:
                    acm[jm - 1, jm - 1] += cfxa1 * fc
                if i1 != 0:
                    acm[jm - 1, i1 - 1] = cfx + cg
                else:
                    acm[jm - 1, jm - 1] += cfxa2 * fc
    elif mode in (23, 32):
        orr = omega * rho
        for i in range(1, m + 1):
            dx1 = dx - (i - ic - 0.5) * tilt / m * dx
            dy1 = dy + (i - ic - 0.5) * tilt / m * dy
            if rho == 0.0:
                rdx = 1.0 / dx1
                rdy = 1.0 / dy1
                ac = 0.5 / omega
                cgy = 0.0
            else:
                x = (i - ic - 0.5) * dx
                r = rho + x
                rdx = r / dx1
                rdy = r / dy1
                if mode == 23:
                    ac = 0.5 * omega / (orr ** 2 - 1.0)
                    cgy = rdy * ac
                else:
                    ac = 0.5 / orr / rho
                    cgy = 0.0
            cfxy = rdx * rdy * ac
            for j in range(1, n + 1):
                jm = ind[i, j]
                if jm == 0:
                    break
                j1 = ind[i, j + 1]
                jm1 = ind[i, j - 1]
                im1 = ind[i - 1, j]
                i1 = ind[i + 1, j]
                if j == 1:
                    acm[jm - 1, jm - 1] = cgy
                if jm1 != 0:
                    acm[jm - 1, jm1 - 1] = -cgy
                if j1 == 0:
                    acm[jm - 1, jm - 1] += cgy
                    if im1 == 0:
                        acm[jm - 1, jm - 1] -= cfxy * fudge1
                    if i1 == 0:
                        acm[jm - 1, jm - 1] += cfxy * fudge1
                else:
                    acm[jm - 1, j1 - 1] = cgy
                    if im1 != 0 and i1 != 0:
                        if ind[i - 1, j + 1] == 0:
                            acm[jm - 1, jm - 1] -= cfxy * fudge3
                        if ind[i + 1, j + 1] == 0:
                            acm[jm - 1, jm - 1] += cfxy * fudge3
    return acm


# ---------------------------------------------------------------------------
# SADScript 組み込み
# ---------------------------------------------------------------------------

def csrphimesh(pipe, dx, dy, ic, nphimax):
    """csrtest.f csrphimesh: ビームパイプ (列ごとの高さ pipe) の境界の点列 (phix, phiy) (1 始まり。0 番は未使用)。戻り: (phix, phiy, nphi, my)。
    nphimax を超えたら途中で打ち切り (nphi >= nphimax)"""
    mx = len(pipe)
    phix = np.zeros(nphimax + 2)
    phiy = np.zeros(nphimax + 2)
    state = {'k': 1, 'full': False}

    def put(x, y):
        state['k'] += 1
        if state['k'] >= nphimax:
            state['full'] = True
            return False
        phix[state['k']] = x
        phiy[state['k']] = y
        return True

    jy = 0
    my = 0
    for i in range(1, mx + 1):
        x = dx * (i - ic - 1)
        l = int(pipe[i - 1])
        my = max(my, l)
        for j in range(jy + 1, l + 1):
            if not put(x, dy * (j - 0.5)):
                return phix, phiy, state['k'], my
        for j in range(jy, l, -1):
            if not put(x, dy * (j - 0.5)):
                return phix, phiy, state['k'], my
        if l > 0:
            if not put((i - ic - 0.5) * dx, l * dy):
                return phix, phiy, state['k'], my
        jy = l
    x = dx * (mx - ic)
    if jy > 0:
        for j in range(jy, 0, -1):
            if not put(x, dy * (j - 0.5)):
                return phix, phiy, state['k'], my
    phix[1] = phix[2]
    phiy[1] = -phiy[2]
    k = state['k']
    put(phix[k], -phiy[k])
    return phix, phiy, state['k'], my


def csr2dpoisson(phix, phiy, dx, dy, pipe, nphi, ic, mx, my, rho0, itern):
    """csrtest.f csr2dpoisson: 境界 (phix, phiy: 1 始まり) と電荷分布 (itern 回の繰り返しで空間電荷 q1 を入れる) から、
    メッシュ上の電場 ex(mx, my)・ey(mx, my)。境界の電荷は tsolva (= tsvd(.false.)) の最小ノルム解"""
    from .tsvd import tsvd
    pipe = np.asarray(pipe, dtype=float)
    ex = np.zeros((mx, my))
    ey = np.zeros((mx, my))
    q1 = np.zeros((mx, my))
    eps = min(1.0e-2, max(1.0 / (nphi - 1) ** 2, 1.0e-4))
    pi = math.pi
    px = np.asarray(phix[1:nphi + 1], dtype=float)
    py = np.asarray(phiy[1:nphi + 1], dtype=float)
    mxp = 0.5 * (px[:-1] + px[1:])
    myp = 0.5 * (py[:-1] + py[1:])
    for it in range(1, itern + 1):
        if it > 1:
            for k in range(1, mx + 1):
                rf = (dx * (k - ic - 0.5) + rho0) * (2.0 * pi) / (dx * dy)
                for l in range(1, int(pipe[k - 1]) + 1):
                    q1[k - 1, l - 1] = -ex[k - 1, l - 1] / rf
        a = np.zeros((nphi - 1, 2 * nphi - 3))
        b = np.zeros(nphi - 1)
        for i in range(2, nphi):
            xf = px[i - 1]
            yf = py[i - 1]
            bs = -math.log(xf ** 2 + yf ** 2)
            a[i - 2, 0:nphi - 1] = -np.log(((xf - mxp) ** 2 + (yf - myp) ** 2) * ((xf - mxp) ** 2 + (yf + myp) ** 2))
            for j in range(2, nphi):
                if j != i:
                    xp = px[j - 1]
                    yp = py[j - 1]
                    a[i - 2, j + nphi - 3] = -math.log(((xf - xp) ** 2 + (yf - yp) ** 2) * ((xf - xp) ** 2 + (yf + yp) ** 2))
            if it > 1:
                for k in range(1, mx + 1):
                    xp = dx * (k - ic - 0.5)
                    for l in range(1, int(pipe[k - 1]) + 1):
                        yp = dy * (l - 0.5)
                        bs = bs - math.log(((xf - xp) ** 2 + (yf - yp) ** 2) * ((xf - xp) ** 2 + (yf + yp) ** 2)) * q1[k - 1, l - 1]
            b[i - 2] = bs
        a[nphi - 2, :] = 1.0
        b[nphi - 2] = 0.5
        rho = tsvd(a, b, eps, False)
        for k in range(1, mx + 1):
            xf = dx * (k - ic - 0.5)
            for l in range(1, int(pipe[k - 1]) + 1):
                yf = dy * (l - 0.5)
                exa = xf / (xf ** 2 + yf ** 2)
                eya = yf / (xf ** 2 + yf ** 2)
                ddx = xf - mxp
                ddy1 = yf - myp
                ddy2 = yf + myp
                rsqi1 = 1.0 / (ddx ** 2 + ddy1 ** 2)
                rsqi2 = 1.0 / (ddx ** 2 + ddy2 ** 2)
                exa = exa - float(np.sum(rho[:nphi - 1] * ddx * (rsqi1 + rsqi2)))
                eya = eya - float(np.sum(rho[:nphi - 1] * (ddy1 * rsqi1 + ddy2 * rsqi2)))
                xs = px[1:nphi - 1]
                ys = py[1:nphi - 1]
                ddx = xf - xs
                ddy1 = yf - ys
                ddy2 = yf + ys
                rsqi1 = 1.0 / (ddx ** 2 + ddy1 ** 2)
                rsqi2 = 1.0 / (ddx ** 2 + ddy2 ** 2)
                rr = rho[nphi - 1:2 * nphi - 3]
                exa = exa - float(np.sum(rr * ddx * (rsqi1 + rsqi2)))
                eya = eya - float(np.sum(rr * (ddy1 * rsqi1 + ddy2 * rsqi2)))
                if it > 1:
                    for i in range(1, mx + 1):
                        x = dx * (i - ic - 0.5)
                        ddx0 = xf - x
                        for j in range(1, int(pipe[i - 1]) + 1):
                            y = dy * (j - 0.5)
                            ddy1_ = yf - y
                            rsq1 = ddx0 ** 2 + ddy1_ ** 2
                            rsqi1_ = 0.0 if rsq1 == 0.0 else 1.0 / rsq1
                            ddy2_ = yf + y
                            rsqi2_ = 1.0 / (ddx0 ** 2 + ddy2_ ** 2)
                            q = q1[i - 1, j - 1]
                            exa = exa + q * ddx0 * (rsqi1_ + rsqi2_)
                            eya = eya + q * (ddy1_ * rsqi1_ + ddy2_ * rsqi2_)
                ex[k - 1, l - 1] = exa
                ey[k - 1, l - 1] = eya
            for l in range(int(pipe[k - 1]) + 1, my + 1):
                ex[k - 1, l - 1] = 0.0
                ey[k - 1, l - 1] = 0.0
    return ex, ey


def _phi1x(x, y, sa):
    """csrtest.f phi1x (スカラーまたは配列)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r = x ** 2 + y ** 2
    with np.errstate(divide='ignore', invalid='ignore'):
        v = np.where(r == 0.0, 0.5 - math.log(sa / math.pi) * 0.5, -np.log(np.where(r == 0.0, 1.0, r)) * 0.5 + sa / 24.0 / np.where(r == 0.0, 1.0, r))
    return v


def _phi1y(x, y):
    """csrtest.f phi1y (-atan(x/y)、y=0 は 0)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        v = np.where(y == 0.0, 0.0, -np.arctan(x / np.where(y == 0.0, 1.0, y)))
    return v


def csrdoublephimesh(phix, phiy, nphi):
    """csrtest.f csrdoublephimesh: 点列の間に中点を足して 2 倍 (配列は 0 始まりの添字 0..2nphi)。戻り: (phix, phiy, nphi_new)"""
    px = np.zeros(nphi * 2 + 1)
    py = np.zeros(nphi * 2 + 1)
    px[:nphi + 1] = phix[:nphi + 1]
    py[:nphi + 1] = phiy[:nphi + 1]
    for i in range(nphi, 0, -1):
        ii = i * 2
        px[ii] = px[i]
        py[ii] = py[i]
        px[ii - 1] = (px[i] + px[i - 1]) * 0.5
        py[ii - 1] = (py[i] + py[i - 1]) * 0.5
    return px, py, nphi * 2 - 1


def _tsolvm_minnorm(a, b, eps):
    """tsolvm(svd=.true.): a x = b の最小ノルム解 (特異値 < 最大 × eps を捨てる。b は複数の右辺)"""
    U, sv, Vt = np.linalg.svd(a, full_matrices=False)
    thr = sv.max() * eps if sv.size else 0.0
    sinv = np.where(sv > thr, 1.0 / np.where(sv > thr, sv, 1.0), 0.0)
    return Vt.T @ (sinv[:, None] * (U.T @ b))


def csrhelmat(pipe, dx, dy, ic, nm):
    """csrtest.f csrhelmat + csr2dpoisson1: らせん断面 (CSRMatrix mode 1000) の (2nm × 2nm) の行列 psi。
    境界の点列 (中点を足した 2 倍の細かさ) から、メッシュの各セル (xc, yc) の電場への応答を境界要素法で求める"""
    mx = len(pipe)
    nphimax = 4000
    phix0, phiy0, nphi, _my = csrphimesh(pipe, dx, dy, ic, nphimax)
    # K64 は phix(0:nphimax) を csrphimesh の phix(*) に渡すので、csrphimesh の k 番目の点は呼び出し側の phix(k-1)
    phix0 = phix0[1:]
    phiy0 = phiy0[1:]
    nphi = nphi - 1
    phix, phiy, nphi = csrdoublephimesh(phix0, phiy0, nphi)         # 0 始まりの添字 0..nphi (Fortran の phix(0:nphi))
    tx = np.zeros(nphi + 1)
    ty = np.zeros(nphi + 1)
    for i in range(1, nphi + 1):
        sx = phix[i + 1] - phix[i - 1] if i + 1 < len(phix) else 0.0
        sy = phiy[i + 1] - phiy[i - 1] if i + 1 < len(phiy) else 0.0
        r = math.hypot(sx, sy)
        tx[i] = sx / r if r != 0.0 else 0.0
        ty[i] = sy / r if r != 0.0 else 0.0
    xc = np.zeros(nm)
    yc = np.zeros(nm)
    xp1 = np.zeros(nm)
    xp2 = np.zeros(nm)
    yp1 = np.zeros(nm)
    ind = 0
    for i in range(1, mx + 1):
        x = (i - ic - 0.5) * dx
        for j in range(1, int(pipe[i - 1]) + 1):
            y = (j - 0.5) * dy
            xc[ind] = x
            yc[ind] = y
            yp1[ind] = pipe[i - 1] * dy
            xp1[ind] = -ic * dx
            xp2[ind] = (mx - ic) * dx
            for k in range(j - 1, 0, -1):
                if int(pipe[k - 1]) < j:
                    xp1[ind] = (k - ic) * dx
                    break
            for k in range(j + 1, mx + 1):
                if int(pipe[k - 1]) < j:
                    xp2[ind] = (k - ic - 1) * dx
                    break
            ind += 1
    cs = dx * dy / 2.0 / math.pi
    sa = cs * 2.0 * math.pi
    ws2 = sa / 16.0
    eps = 1.0e-4
    n1 = nphi - 1
    a1 = np.zeros((n1, nphi))
    b1 = np.zeros((n1, nm * 2))
    px = phix[0:nphi + 1]
    py = phiy[0:nphi + 1]
    for i in range(2, nphi):                     # Fortran i = 2..nphi-1 (行 i)
        xf = px[i]
        yf = py[i]
        b1[i - 1, 0:nm] = (-tx[i] * (_phi1x(xf - xc, yf - yc, sa) + _phi1x(xf - xc, yf + yc, sa))
                           - ty[i] * (-_phi1y(yf - yc, xf - xc) - _phi1y(yf + yc, xf - xc)))
        b1[i - 1, nm:2 * nm] = (-tx[i] * (_phi1y(yf - yc, xf - xc) - _phi1y(yf + yc, xf - xc))
                                - ty[i] * (_phi1x(yf - yc, xf - xc, sa) - _phi1x(yf + yc, xf - xc, sa)))
        for j in range(1, nphi + 1):
            xp = px[j]
            yp = py[j]
            if i != j:
                rsqi1 = 1.0 / ((xf - xp) ** 2 + (yf - yp) ** 2)
                w1 = rsqi1 * (1.0 + ws2 * rsqi1)
            else:
                w1 = 0.0
            rsqi2 = 1.0 / ((xf - xp) ** 2 + (yf + yp) ** 2)
            w2 = rsqi2 * (1.0 + ws2 * rsqi2)
            a1[i - 1, j - 1] = tx[i] * (xf - xp) * (w1 + w2) + ty[i] * ((yf - yp) * w1 + (yf + yp) * w2)
        a1[0, i - 1] = 1.0
    a1[0, 0] = 1.0
    a1[0, nphi - 1] = 1.0
    b1[0, 0:nm] = -yp1 + yc
    b1[0, nm:2 * nm] = yc
    er1 = _tsolvm_minnorm(a1, b1, eps)          # (nphi, 2nm)
    psi = np.zeros((nm * 2, nm * 2))
    for k in range(nm):
        xf = xc[k]
        yf = yc[k]
        exa1 = _phi1x(xf - xc, yf - yc, sa) + _phi1x(xf - xc, yf + yc, sa)
        eya1 = -_phi1y(yf - yc, xf - xc) - _phi1y(yf + yc, xf - xc)
        exa2 = _phi1y(yf - yc, xf - xc) - _phi1y(yf + yc, xf - xc)
        eya2 = _phi1x(yf - yc, xf - xc, sa) - _phi1x(yf + yc, xf - xc, sa)
        for j in range(1, nphi + 1):
            xp = px[j]
            yp = py[j]
            rsqi1 = 1.0 / ((xf - xp) ** 2 + (yf - yp) ** 2)
            rsqi2 = 1.0 / ((xf - xp) ** 2 + (yf + yp) ** 2)
            w1 = rsqi1 * (1.0 + ws2 * rsqi1)
            w2 = rsqi2 * (1.0 + ws2 * rsqi2)
            exa1 = exa1 + (xf - xp) * (w1 + w2) * er1[j - 1, 0:nm]
            eya1 = eya1 + ((yf - yp) * w1 + (yf + yp) * w2) * er1[j - 1, 0:nm]
            exa2 = exa2 + (xf - xp) * (w1 + w2) * er1[j - 1, nm:2 * nm]
            eya2 = eya2 + ((yf - yp) * w1 + (yf + yp) * w2) * er1[j - 1, nm:2 * nm]
        psi[k, 0:nm] = exa1 * cs
        psi[k, nm:2 * nm] = exa2 * cs
        psi[nm + k, 0:nm] = eya1 * cs
        psi[nm + k, nm:2 * nm] = eya2 * cs
    return psi


def install(I):
    from .sadscript import SADError, is_num
    B = I.builtins

    def reallist(x, what):
        if not isinstance(x, list) or not all(is_num(v) for v in x):
            raise SADError('General::wrongtype: Wrong type of argument in ' + what)
        return x

    def gcut():
        g = I.get('GCUT', 1e300) if hasattr(I, 'get') else 1e300
        return float(g) if is_num(g) and g > 0 else 1e300

    def gauss(n):
        if I.rng is None:
            return np.random.normal(size=n)
        v = I.rng.gauss(n, gcut())
        return [v] if n == 1 and not isinstance(v, list) else v

    def csrtrack_f(a):
        if len(a) != 5:
            raise SADError('General::narg: Number of arguments is expected to be 5 in CSRTrack')
        zl_, dpl_ = reallist(a[0], 'CSRTrack'), reallist(a[1], 'CSRTrack')
        if len(dpl_) != len(zl_):
            raise SADError('General::wrongleng: Length of z is expected to be equal to length of dp in CSRTrack')
        zlist = reallist(a[2], 'CSRTrack')
        p = reallist(a[3], 'CSRTrack')
        if len(p) != 9:
            raise SADError('General::wrongtype: Wrong type of argument in CSRTrack')
        npart = len(zl_)
        dz, ndiv = p[0], int(max(1.0, p[1]))
        aw = p[2] / npart
        z = np.array(zl_, dtype=float)
        dp = np.array(dpl_, dtype=float)
        nt = int(a[4])
        data = csrtrack(z, dp, np.array(zlist, dtype=float), dz, ndiv, aw, p[3], p[4], p[5], p[6], p[7], p[8], nt,
                        lambda n: gauss(n))
        zl_[:] = [float(v) for v in z]             # K64 はリスト本体を直接書き換える
        dpl_[:] = [float(v) for v in dp]
        return [float(v) for v in data.T.reshape(-1)]

    def csrhaissin_f(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in CSRHaissin')
        zl = np.array(reallist(a[0], 'CSRHaissin'), dtype=float)
        p = reallist(a[1], 'CSRHaissin')
        if len(p) != 4:
            raise SADError('General::wrongtype: Wrong type of argument in CSRHaissin')
        dz, r65, sigz = p[0], p[2], p[3]
        aw = p[1] / r65 / sigz ** 2
        rho = csrhaissin(zl, aw, dz, sigz)
        if rho is None:
            raise SADError('General::toomany: Too many iterations in CSRHaissin')
        return [float(v) for v in rho]

    def csroysetup_f(a):
        if len(a) != 3:
            raise SADError('General::narg: Number of arguments is expected to be 3 in CSROySetup')
        rho = np.array(reallist(a[0], 'CSROySetup'), dtype=float)
        zl = np.array(reallist(a[1], 'CSROySetup'), dtype=float)
        p = reallist(a[2], 'CSROySetup')
        if len(p) != 7:
            raise SADError('General::wrongtype: Wrong type of argument in CSROySetup')
        dz, nj, nphi, aw, r56, r65, zspan = p[0], int(p[1]), int(p[2]), p[3], p[4], p[5], p[6]
        mphi = nphi * 16
        dt = min(1.0, 0.01 / math.sqrt(abs(r56 * r65)))
        wi, aj, omega, zjn, dpjn, rhoj, irtc = csroysetup(rho, zl, aw, r56, r65, dz, nj, mphi, zspan, dt)
        f = lambda x: [float(v) for v in np.asarray(x).reshape(-1)]
        res = [f(wi), f(aj), f(omega), f(zjn), f(dpjn), f(rhoj)]
        if irtc != 0:
            raise SADError('CSR::irreg: Irregular distribution in CSROySetup')
        return res

    def csroymat_f(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in CSROyMatrix')
        las, p = a
        if not isinstance(las, list) or len(las) != 6:
            raise SADError('General::wrongtype: Wrong type of argument in CSROyMatrix')
        wi = np.array(las[0], dtype=float)
        omega = np.array(las[2], dtype=float)
        zj = np.array(las[3], dtype=float)
        rhoj = np.array(las[5], dtype=float)
        nr, nj = len(wi), len(omega)
        nphi = len(zj) // nj // 16
        p = reallist(p, 'CSROyMatrix')
        if len(p) != 5:
            raise SADError('General::wrongtype: Wrong type of argument in CSROyMatrix')
        dz, anus, sigz, sige = p[0], p[2], p[3], p[4]
        bw = p[1] * 4.0 / math.pi / (2.0 * math.pi * anus * sige * sigz)
        bm = csroymat(wi, omega, zj.reshape(nj, nphi * 16), rhoj, bw, dz, nj, nphi, nr)
        return [[float(v) for v in row] for row in bm]

    B['CSRSetupOY'] = csroysetup_f
    B['CSROYMatrix'] = csroymat_f
    def csrconvert_f(a):
        if len(a) != 2:
            raise SADError('General::narg: Number of arguments is expected to be 2 in CSRConvert')
        kl, pipe = a
        wrong = 'General::wrongtype: Wrong type of argument in CSRConvert (Vector or Matrix, Pipe)'
        if not isinstance(kl, list) or not isinstance(pipe, list) or not all(is_num(v) for v in pipe):
            raise SADError(wrong)
        ind, mm, n1 = csrind(pipe)
        m1 = len(pipe)
        if len(kl) == m1:
            out = [0.0] * mm
            for i in range(1, m1 + 1):
                row = kl[i - 1]
                if not isinstance(row, list) or not all(is_num(v) for v in row):
                    raise SADError(wrong)
                for j in range(1, n1 + 1):
                    jm = ind[i, j]
                    if jm != 0:
                        out[jm - 1] = float(row[j - 1])
            return out
        if len(kl) == mm:
            return [[float(kl[ind[i, j] - 1]) if ind[i, j] != 0 else 0.0 for j in range(1, n1 + 1)]
                    for i in range(1, m1 + 1)]
        raise SADError(wrong)

    def csrmatrix_f(a):
        wrong = ('General::wrongtype: Wrong type of argument in CSRMatrix (Pipe, {dx, dy, fudge0, fudge1, fudge2, fudge3, tilt},'
                 ' omega/c, rho, center, mode)')
        if len(a) != 6:
            raise SADError('General::narg: Number of arguments is expected to be 6 in CSRMatrix')
        pipe, p = a[0], a[1]
        if not isinstance(pipe, list) or not all(is_num(v) for v in pipe) or not isinstance(p, list) or len(p) != 7:
            raise SADError(wrong)
        if not all(is_num(v) for v in a[2:]) or not all(is_num(v) for v in p):
            raise SADError(wrong)
        dx, dy, f0, f1, f2, f3, tilt = p
        omega, rho, ic, mode = a[2], a[3], int(a[4]), int(a[5])
        if mode == 1000:
            nm_ = int(sum(int(v) for v in pipe))
            psi = csrhelmat(pipe, dx, dy, ic, nm_)
            return [[float(v) for v in row] for row in psi]
        ind, mm, n1 = csrind(pipe)
        acm = csrmat(ind, dx, dy, omega, rho, f0, f1, f2, f3, tilt, mode, mm, len(pipe), n1, ic)
        if mode >= 100:
            return [float(v) for v in acm[:, 0]]
        return [[float(v) for v in row] for row in acm]


    def csrinit_f(a):
        wrong = 'General::wrongtype: {dx or {dx, dy}, center, rho, pipe:{y1, y2, ...}} for #1'
        if len(a) != 1:
            raise SADError('General::narg: Number of arguments is expected to be 1 in CSRInit')
        kl = a[0]
        if not isinstance(kl, list) or len(kl) != 4:
            raise SADError(wrong)
        d0 = kl[0]
        if isinstance(d0, list):
            if len(d0) != 2 or not all(is_num(v) for v in d0):
                raise SADError(wrong)
            dx, dy = float(d0[0]), float(d0[1])
        elif is_num(d0):
            dx = 0.0                      # K64 (tfcsrinit) は実数 1 つのとき dx を読まず (dx=0 のまま) dy=dx とするので、常に wrongtype になる
            dy = dx
        else:
            raise SADError(wrong)
        if dx <= 0.0 or dy <= 0.0 or not is_num(kl[1]) or int(kl[1]) <= 0 or not is_num(kl[2]) or not isinstance(kl[3], list) \
                or not all(is_num(v) for v in kl[3]):
            raise SADError(wrong)
        ic = int(kl[1])
        rho = float(kl[2])
        pipe = [float(v) for v in kl[3]]
        mx = len(pipe)
        if mx <= ic:
            raise SADError('General::wrongval: Center must be within beam pipe')
        nphimax = 16384
        phix, phiy, nphi, my = csrphimesh(pipe, dx, dy, ic, nphimax)
        if nphi >= nphimax:
            raise SADError('General::wrongval: Too large number of meshes')
        ex, ey = csr2dpoisson(phix, phiy, dx, dy, pipe, nphi, ic, mx, my, rho, 3)
        return [[[float(v) for v in row] for row in ex], [[float(v) for v in row] for row in ey]]

    B['CSRInit'] = csrinit_f

    B['CSRConvert'] = csrconvert_f
    B['CSRMatrix'] = csrmatrix_f
    B['CSRTrack'] = csrtrack_f
    B['CSRHaissin'] = csrhaissin_f
