# -*- coding: utf-8 -*-
"""Direct k64-mode translation of PSPAC space charge (tpspac in K64 tspch_.c, spch.c, psn.c), including the K64 shape function.
spch_k64.py — PSPAC (K64 の tspch_.c tpspac + spch.c + psn.c) の k64 モードの直訳 (形状関数まで K64 のまま)。

v2 の PSPAC (`spch.Solver2D`・`kick2d`) は K64 と同じ物理モデル (横 2 次元のポアソン方程式 × 縦の線密度) を、別の形状関数
(2 次の B スプライン) と正しい係数で書き直したもの。粗いメッシュ (横の rms が数セル以下、縦が数セル以下) では gs (K64) と値が違う。
k64 モードでは K64 の形状関数をそのまま使い、gs と同じ値を出す。

K64 の手順 (spch.c):
  1. 電荷の配分 (rhocal): 横は CIC (双線形。粒子を囲む 4 ノード)、縦は「ほぼ最近接」(セルの前半 (dzic < 0.5) は重み 1 をすべて下のノードへ、
     後半は 2(1−dzic)² を下のノード、残りを上のノードへ)。横の範囲 (|x| ≤ dx (nx/2 − 1)、|y| ≤ dy (ny/2 − 1)) の外の粒子は電荷に入らない。
     縦は |z| ≤ dz (nz/2 − 1) の外は線密度に入らない。
  2. 横のポアソン方程式 (psn.c の FACR): 5 点差分 φ(i±1,j) + (dx/dy)² φ(i,j±1) − 2(1 + (dx/dy)²) φ = −2π (dx/dy) ρ、境界で φ = 0
     (∇²φ = −2π ρ/(dx dy)、ρ は個数。FACR は離散方程式の厳密解なので、ここでは DST で同じ解を求める)。
  3. 電位の補間 (makebcu、epotkick): 各セル (i, j) の 4 隅の φ・∂φ/∂i・∂φ/∂j・∂²φ/∂i∂j (中心差分、セル単位) から双 3 次補間 (Numerical Recipes の bcucof)。
  4. 線密度 λ(z) (rhoz): ノードの個数を 3 次スプライン (端の傾き 0、ノード数 nz) で補間 (splint)。
  5. キック (epotkick): Kfac = K L /(np² dz)、K = 2 N_b r_0/(β²γ³) (r_0 = r_e m_e/mass、γ は最初の呼び出しの h0)。
        px_canon −= Kfac λ(z) ∂φ/∂x、py_canon −= Kfac λ ∂φ/∂y、δ −= Kfac λ'(z) φ  (np = そのときの粒子数)。
     (px, py) は正準運動量 (1 + δ) を掛けた値で、キックのあと新しい δ で割り戻す。dv = −δ (2 + δ)/h1/(h1 + (1+δ) h0) (dvfs は足さない)。
     z は変えない (v2 は β の比で z を直す)。
  6. PSPACNTURNCALC ターンを超えると電位 (SP[isp]) は凍結して使い回す。λ (rhoz) は大域変数で、最後に計算した SPCH のもののまま。
     Kfac の np も最後に計算したときの粒子数のまま。

K64 の癖 (再現しないもの):
  - gs は PSPAC の状態を C の static 変数に持つので、1 セッションで TrackParticles を繰り返すと 2 回目から働かない (SPC8)。ここでは TrackParticles ごとに新しく始める。
  - PSPACNTURN ターン後にポテンシャルソルバーのメモリを解放するので、それ以降の呼び出しは (C の規則では) 未定義の動作。gs は解放済みのメモリの中身
    (最後に計算した電位 SP[isp] と λ) をそのまま読み、凍結した電位と同じ値になる (2 周のうち PSPACNTURN=1 で gs と全桁一致を確かめた)。ここでも同じ。
  - PSPACNZ=1 (既定) は gs が SIGSEGV (B93)。ここでは PSPACNZ < 3 をエラーにする。
  - 電位のメッシュの端のセル (ii = nx − 1 または jj = ny − 1) の bcucof は配列の外を読む。ここでは外側の φ を 0 とする。
  - 縦のメッシュの kk = nz − 1 のセルの 3 次スプライン係数 rhoz_sp[nz] は未初期化。ここでは 0 とする。
"""
from __future__ import annotations

import math

import numpy as np
from scipy import fft as _fft

from .sadscript import SADError

M_E = 0.51099906e6               # spch.c: m_e [eV]
R_E = 2.81794092e-15             # spch.c: r_e [m]

_WT = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [-3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0],
    [2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1],
    [0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1],
    [-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0],
    [9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2],
    [-6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2],
    [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0],
    [-6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1],
    [4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1]], dtype=float)


def splined(x, y, yp1=0.0, ypn=0.0):
    """Cubic spline coefficients as in `splined` of spch.c (Numerical Recipes spline with end first derivatives yp1, ypn).
    spch.c の splined (Numerical Recipes の spline。端の 1 階微分が yp1・ypn (0 なら傾き 0 の 3 次スプライン)。n = len(x))"""
    n = len(x)
    y2 = np.zeros(n)
    u = np.zeros(n)
    if yp1 > 0.99e33:
        y2[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)
    for i in range(1, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
    if ypn > 0.99e30:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]))
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)
    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
    return y2


class K64Pspac:
    """K64 PSPAC space-charge solver (tpspac_), called from `SpaceCharge._kick2d`; one instance per TrackParticles run.
    K64 の PSPAC (tpspac_)。SpaceCharge._kick2d から呼ぶ。1 回の TrackParticles に 1 つ"""

    def __init__(self, nx, ny, nz, dx, dy, dz, rnp, amass, nturnend=1, nturncalc=0, warn=None):
        nx, ny, nz = max(1, int(nx)), max(1, int(ny)), max(1, int(nz))
        if nx < 4 or ny < 4:
            raise SADError('PSPAC: PSPACNX, PSPACNY は 4 以上')
        if nz < 3:
            raise SADError('PSPAC (k64): PSPACNZ < 3 は K64 (gs) が SIGSEGV で落ちる (B93)。PSPACNZ=32 などにする')
        self.nx, self.ny, self.nz = nx, ny, nz
        self.dx, self.dy, self.dz = float(dx), float(dy), float(dz)
        self.rnp, self.mass = float(rnp), float(amass)
        self.nturnend = max(1, int(nturnend))
        self.nturncalc = int(nturncalc) if int(nturncalc) > 0 else self.nturnend
        self.warn = warn
        self.ixc, self.jyc, self.kzc = nx // 2, ny // 2, nz // 2
        self.xc, self.yc, self.zc = self.ixc * self.dx, self.jyc * self.dy, self.kzc * self.dz
        self.xm, self.ym, self.zm = self.dx * (self.ixc - 1), self.dy * (self.jyc - 1), self.dz * (self.kzc - 1)
        self.z_sp = (np.arange(nz + 1) - self.kzc) * self.dz
        px = np.arange(1, nx)
        py = np.arange(1, ny)
        ex = 2.0 - 2.0 * np.cos(math.pi * px / nx)
        ey = 2.0 - 2.0 * np.cos(math.pi * py / ny)
        self._eig = (ex / self.dx ** 2)[:, None] + (ey / self.dy ** 2)[None, :]
        self.nturn = 1
        self.isp = 0
        self.ispend = 0
        self.kcoeff = None               # SpaceChargeBPSet (最初の呼び出しの h0)
        self.btab = {}                   # SP[isp]: (nx + 1, ny + 1, 16)
        self.rhoz = np.zeros(nz + 1)
        self.np_psn = 1
        self._warned = False
        self.nkick = 0

    # ---- 1 回の呼び出し (tpspac_)
    def kick(self, S, n, al, ph):
        """Apply the space-charge kick of length al to the first n particles of S; returns n.
        S の先頭 n 個にキック (長さ al)。戻り: n。SpaceCharge の 1 回の SPCH 要素 (または区間の終わり) ごとに 1 回呼ぶ"""
        if self.kcoeff is None:                              # SpaceChargeBPSet(rnp, Mass, h0): 最初の呼び出しの h0
            r0 = R_E * M_E / self.mass
            beta = math.sqrt(1.0 - 1.0 / ph.h0 / ph.h0)
            self.kcoeff = 2.0 * self.rnp * r0 / beta / beta / ph.h0 ** 3
        p0, h0 = ph.p0, ph.h0
        x, y, z = S[:n, 0], S[:n, 2], S[:n, 4]
        g = S[:n, 5]
        pr = 1.0 + g
        px = S[:n, 1] * pr
        py = S[:n, 3] * pr
        isp = self.isp
        calc = self.nturn <= self.nturncalc
        if calc:
            self._rhocal(x, y, z)
            self.np_psn = n
            self.btab[isp] = self._bcu(self._poisson(self._rhoxy))
        btab = self.btab.get(isp)
        if btab is None:                                     # 凍結のときに無い isp (SPCH の数が増えた): キックなし
            self._advance(False)
            return n
        y2 = splined(self.z_sp[:self.nz], self.rhoz[:self.nz])
        y2 = np.append(y2, 0.0)
        self._epotkick(x, y, z, g, px, py, btab, y2, al)
        pr = 1.0 + g
        S[:n, 1] = px / pr
        S[:n, 3] = py / pr
        h1 = np.sqrt(1.0 + p0 * p0 * pr * pr)
        S[:n, 6] = -g * (1.0 + pr) / h1 / (h1 + pr * h0)
        self.nkick += 1
        self._advance(False)
        return n

    def end_of_turn(self):
        """Reset isp to 0 and advance the turn, as in the nlat == nlatend branch of tpspac_.
        tpspac_ の「nlat == nlatend」: isp を 0 に戻し、ターンを進める"""
        self.ispend = self.isp
        self.isp = 0
        self.nturn += 1

    def _advance(self, _):
        self.isp += 1

    # ---- 電荷の配分 (rhocal)
    def _rhocal(self, x, y, z):
        nx, ny, nz = self.nx, self.ny, self.nz
        rhoxy = np.zeros((nx + 1) * (ny + 1))
        rhoz = np.zeros(nz + 1)
        sel = ~((np.abs(x) > self.xm) | (np.abs(y) > self.ym))      # C: `if(fabs>..) continue` (NaN は continue しない)
        xs, ys, zs = x[sel], y[sel], z[sel]
        ii = np.trunc((xs + self.xc) / self.dx).astype(np.int64)
        jj = np.trunc((ys + self.yc) / self.dy).astype(np.int64)
        xic = (ii - self.ixc) * self.dx
        yjc = (jj - self.jyc) * self.dy
        rxin = np.abs((xs - xic) / self.dx)
        rxi = 1.0 - rxin
        ryjn = np.abs((ys - yjc) / self.dy)
        ryj = 1.0 - ryjn
        iijj = ii * (ny + 1) + jj
        sz = rhoxy.size
        rhoxy += np.bincount(iijj, weights=rxi * ryj, minlength=sz)[:sz]
        rhoxy += np.bincount(iijj + ny + 1, weights=rxin * ryj, minlength=sz)[:sz]
        rhoxy += np.bincount(iijj + ny + 2, weights=rxin * ryjn, minlength=sz)[:sz]
        rhoxy += np.bincount(iijj + 1, weights=rxi * ryjn, minlength=sz)[:sz]
        self._rhoxy = rhoxy.reshape(nx + 1, ny + 1)
        selz = np.abs(zs) <= self.zm
        zz = zs[selz]
        kz = np.trunc((zz + self.zc) / self.dz).astype(np.int64)
        zic = (kz - self.kzc) * self.dz
        dzic = (zz - zic) / self.dz
        lo = dzic < 0.5
        rzin = np.where(lo, 2.0 * dzic * dzic, 1.0 - 2.0 * (1.0 - dzic) ** 2)
        rzi = np.where(lo, 1.0 - rzin, 2.0 * (1.0 - dzic) ** 2)
        # dzic < 0.5: rhoz[ii] += rzi + rzin (= 1)。それ以外: rhoz[ii] += rzi、rhoz[ii+1] += rzin
        rhoz += np.bincount(kz, weights=np.where(lo, rzi + rzin, rzi), minlength=nz + 1)[:nz + 1]
        up = ~lo
        rhoz += np.bincount(kz[up] + 1, weights=rzin[up], minlength=nz + 1)[:nz + 1]
        self.rhoz = rhoz

    # ---- ポアソン方程式 (psn.c の FACR と同じ離散方程式の解)
    def _poisson(self, rho):
        nx, ny = self.nx, self.ny
        rhs = 2.0 * math.pi * rho[1:nx, 1:ny] / (self.dx * self.dy)       # 5 点差分 ∇²φ = −2π ρ/(dx dy)
        phi = np.zeros((nx + 1, ny + 1))
        phi[1:nx, 1:ny] = _fft.idstn(_fft.dstn(rhs, type=1) / self._eig, type=1)
        return phi

    # ---- 双 3 次補間の係数 (makebcu)
    def _bcu(self, phi):
        nx, ny = self.nx, self.ny
        P = np.zeros((nx + 4, ny + 4))
        P[1:nx + 2, 1:ny + 2] = phi              # P[i + 1, j + 1] = φ(i, j)。外側は 0 (K64 は配列の外を読む)

        def f(di, dj):
            return P[2 + di:nx + 1 + di, 2 + dj:ny + 1 + dj]      # i = 1..nx−1, j = 1..ny−1 のセル → φ(i + di, j + dj)
        F = np.zeros((nx - 1, ny - 1, 16))
        F[..., 0] = f(0, 0)
        F[..., 1] = f(1, 0)
        F[..., 2] = f(1, 1)
        F[..., 3] = f(0, 1)
        F[..., 4] = (f(1, 0) - f(-1, 0)) * 0.5
        F[..., 5] = (f(2, 0) - f(0, 0)) * 0.5
        F[..., 6] = (f(2, 1) - f(0, 1)) * 0.5
        F[..., 7] = (f(1, 1) - f(-1, 1)) * 0.5
        F[..., 8] = (f(0, 1) - f(0, -1)) * 0.5
        F[..., 9] = (f(1, 1) - f(1, -1)) * 0.5
        F[..., 10] = (f(1, 2) - f(1, 0)) * 0.5
        F[..., 11] = (f(0, 2) - f(0, 0)) * 0.5
        F[..., 12] = (f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) * 0.25
        F[..., 13] = (f(2, 1) - f(2, -1) - f(0, 1) + f(0, -1)) * 0.25
        F[..., 14] = (f(2, 2) - f(2, 0) - f(0, 2) + f(0, 0)) * 0.25
        F[..., 15] = (f(1, 2) - f(1, 0) - f(-1, 2) + f(-1, 0)) * 0.25
        B = np.zeros((nx + 1, ny + 1, 16))
        B[1:nx, 1:ny, :] = F @ _WT.T
        return B

    # ---- キック (epotkick)
    def _epotkick(self, x, y, z, g, px, py, btab, y2, al):
        nx, ny, nz = self.nx, self.ny, self.nz
        dx, dy, dz = self.dx, self.dy, self.dz
        kfac = self.kcoeff * al / self.np_psn / dz / self.np_psn
        ok = ~((np.abs(x) > self.xm) | (np.abs(y) > self.ym))
        idx = np.nonzero(ok)[0]
        ii = np.trunc((x[idx] + self.xc) / dx).astype(np.int64)
        jj = np.trunc((y[idx] + self.yc) / dy).astype(np.int64)
        kk = np.trunc((z[idx] + self.zc) / dz).astype(np.int64)
        good = (ii >= 1) & (ii <= nx - 1) & (jj >= 1) & (jj <= ny - 1) & (kk >= 0) & (kk < nz)
        idx, ii, jj, kk = idx[good], ii[good], jj[good], kk[good]
        if idx.size == 0:
            return
        xl = (ii - self.ixc) * dx
        yl = (jj - self.jyc) * dy
        zl = (kk - self.kzc) * dz
        t = (x[idx] - xl) / dx
        u = (y[idx] - yl) / dy
        v = (z[idx] - zl) / dz
        b = btab[ii, jj, :]                                   # (m, 16)

        def c(m):
            return b[:, m]
        a4 = ((c(15) * u + c(14)) * u + c(13)) * u + c(12)
        a3 = ((c(11) * u + c(10)) * u + c(9)) * u + c(8)
        a2 = ((c(7) * u + c(6)) * u + c(5)) * u + c(4)
        a1 = ((c(3) * u + c(2)) * u + c(1)) * u + c(0)
        ph = a1 + (a2 + (a3 + a4 * t) * t) * t
        a4 = (3.0 * c(15) * u + 2.0 * c(14)) * u + c(13)
        a3 = (3.0 * c(11) * u + 2.0 * c(10)) * u + c(9)
        a2 = (3.0 * c(7) * u + 2.0 * c(6)) * u + c(5)
        a1 = (3.0 * c(3) * u + 2.0 * c(2)) * u + c(1)
        phy = a1 + (a2 + (a3 + a4 * t) * t) * t
        a4 = (3.0 * c(15) * t + 2.0 * c(11)) * t + c(7)
        a3 = (3.0 * c(14) * t + 2.0 * c(10)) * t + c(6)
        a2 = (3.0 * c(13) * t + 2.0 * c(9)) * t + c(5)
        a1 = (3.0 * c(12) * t + 2.0 * c(8)) * t + c(4)
        phx = a1 + (a2 + (a3 + a4 * u) * u) * u
        phx = phx / dx
        phy = phy / dy
        rz = self.rhoz
        w = 1.0 - v
        lam = w * rz[kk] + v * rz[kk + 1] + ((w * w * w - w) * y2[kk] + (v * v * v - v) * y2[kk + 1]) * (dz * dz) / 6.0
        dlam = (rz[kk + 1] - rz[kk]) / dz - ((3.0 * w * w - 1.0) * y2[kk] - (3.0 * v * v - 1.0) * y2[kk + 1]) * dz / 6.0
        px[idx] -= kfac * lam * phx
        py[idx] -= kfac * lam * phy
        g[idx] -= kfac * dlam * ph
