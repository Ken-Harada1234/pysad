"""WakeFunction (縦・横のウェイク場) の追跡への適用。K64 の txwake.f (module wakez) の直訳。

SAD の使い方 (SADHelp の WakeFunction):
    WakeFunction[Longitudinal, "要素名"] = {{z1, wl1}, ..., {zn, wln}};   (z: m、w: V/C)
    WakeFunction[Transverse,   "要素名"] = {{z1, wt1}, ..., {zn, wtn}};   (w: V/C/m)
  TRPT かつ LWAKE (縦) / TWAKE (横) が ON のとき、TrackParticles だけで効く (tftrack.f:120)。
  PBUNCH (バンチあたり粒子数)、NBUNCH (バンチ数)、CHARGE を使う。

追跡側の呼び方 (tturn.f:256-269, 562-575): 要素 l にウェイクがあれば、
    要素の前に txwake(fw, init=True)、要素を通し、要素の後に txwake(fw, init=False)    fw = |charge| e pbunch nbunch / amass / np0 / 2
  ただし CAVI は tcav.f の中で分割ごとに呼ばれる (track6n.cav_full / k_cav_in・mid・out)。MULT は VOLT≠0 のとき tmulti の中で呼ばれるが、MULT の RF の追跡自体が無いので未対応。
  SOL の領域の内側は tsol.f が要素ごとに呼ぶ (turn の通常の前後と同じ。dy=DY、theta=ROTATE)。

K64 のソースのままにした癖 (未確認の項目は docs/23_wake.md):
  - tturn.f:263-264 の `dy=cmp%value(krot)` (ROTATE を dy に代入し、実際の dy は捨てる)、rot は未初期化の局所変数 → ここでは theta=0、dy=ROTATE の値。
  - txwdefslice のバンチ分割は、バンチ数 nb より多い間隔があると配列の外に書く (mb(nb+1)=np)。ここでは mb(nb+1)=np が読めるとして扱う。
  - 表の最後の 1 行は追跡では使われない (twxiwp の lwl0 = 行数 - 1。gs の実測で確認)。
  - 粒子の整列 (spsort) は z の昇順 (同じ z の並びは SAD の spsort と違うことがある)。
  - txwake の dv は pa (キック前の 1+g) で h1 を作る (K64 のまま)。
"""
import math

import numpy as np
from numba import njit

from .emit6d import ELEMCH, ELMASS

NJ = dict(cache=True, error_model='numpy', nogil=True)

DZLIM = 1e-5
PMIN = 1e-10
ZMIN = -1e30


@njit(**NJ)
def _spsort(zz):
    """spkick.f spsort の直訳 (クイックソート。同じ z の粒子の並びを SAD と同じにする)。返り値は 0 始まりの粒子の並び"""
    n0 = zz.shape[0]
    itab = np.arange(n0)
    stack_off = np.zeros(2 * 64 + 2 * n0 // 2 + 8, dtype=np.int64)
    stack_n = np.zeros(2 * 64 + 2 * n0 // 2 + 8, dtype=np.int64)
    sp = 0
    stack_off[0] = 0
    stack_n[0] = n0
    sp = 1
    while sp > 0:
        sp -= 1
        off = stack_off[sp]
        n = stack_n[sp]
        if n <= 1:
            continue
        i1 = itab[off]
        i2 = itab[off + n - 1]
        if zz[i1] > zz[i2]:
            i1, i2 = i2, i1
        if n == 2:
            itab[off] = i1
            itab[off + 1] = i2
            continue
        m = (n + 1) // 2
        im = itab[off + m - 1]
        if zz[i1] < zz[im]:
            if zz[im] > zz[i2]:
                im, i2 = i2, im
        elif zz[i1] != zz[im]:
            im, i1 = i1, im
        itab[off] = i1
        itab[off + n - 1] = i2
        itab[off + m - 1] = im
        if n == 3:
            continue
        itab[off + m - 1] = itab[off + 1]
        itab[off + 1] = im
        ip1 = 3
        ip2 = n - 1
        while ip1 <= ip2:
            while ip1 <= ip2 and zz[itab[off + ip1 - 1]] <= zz[im]:
                ip1 += 1
            while ip1 <= ip2 and zz[im] <= zz[itab[off + ip2 - 1]]:
                ip2 -= 1
            if ip2 > ip1:
                t = itab[off + ip1 - 1]
                itab[off + ip1 - 1] = itab[off + ip2 - 1]
                itab[off + ip2 - 1] = t
                ip1 += 1
                ip2 -= 1
        ip1 -= 1
        t = itab[off + ip1 - 1]
        itab[off + ip1 - 1] = im
        itab[off + 1] = t
        if sp + 2 > stack_off.shape[0]:
            return itab                      # 起こらない (スタックは十分大きい)
        stack_off[sp] = off
        stack_n[sp] = ip1 - 1
        sp += 1
        stack_off[sp] = off + ip1
        stack_n[sp] = n - ip1
        sp += 1
    return itab


@njit(**NJ)
def _defslice(np_, zs, nb, izs):
    """txwdefslice。zs[1..np_] は z の昇順 (1 始まり)。izs[1..np_] に切片の番号を返す"""
    za = zs[np_] - zs[1]
    if nb > 1:
        sb = za / (nb - 1) * 0.5
    else:
        sb = za
    mb = np.zeros(nb + 2, dtype=np.int64)          # Fortran は mb(0:nb)。nb+1 は配列外への書き込み (上記の癖) の置き場
    m = 1
    izs[1] = 1
    i = 2
    while i <= np_:
        if zs[i] - zs[i - 1] > sb:
            mb[m] = i - 1
            m += 1
            if m > nb:
                for j in range(i, np_ + 1):
                    izs[j] = m
                mb[m] = np_
                break
        izs[i] = m
        i += 1
    if m == 1:
        mb[1] = np_
    for mm in range(nb, 0, -1):
        if mb[mm] != 0:
            nm = mb[mm] - mb[mm - 1]
            m1 = int(math.sqrt(float(nm)))
            mb[mm] = (nm + m1 - 1) // m1
    m0 = 1
    nm = mb[1] - 1
    ks = 1
    for i in range(2, np_ + 1):
        if izs[i] != m0:
            m0 += 1
            ks += 1
            izs[i] = ks
            nm = mb[m0] - 1
        else:
            izs[i] = ks
            nm -= 1
            if nm <= 0:
                nm = mb[m0]
                ks += 1


@njit(**NJ)
def _kick(S, np_, itab, izs, dzwr, fw, p0, h0, wakel, waket, lwl, lwt):
    """txwake 本体 (並べ替えと切片の定義は済んでいるとする)。S の列: 0 x, 1 px, 2 y, 3 py, 4 z, 5 g, 6 dv"""
    fwp = fw / p0
    ns1 = izs[np_]
    ns = ns1 + 1
    ws = np.zeros(ns + 1)
    xs = np.zeros(ns + 1)
    ys = np.zeros(ns + 1)
    zs = np.zeros(ns + 1)
    wx = np.zeros(ns + 1)
    wy = np.zeros(ns + 1)
    wz = np.zeros(ns + 1)
    for i in range(1, np_ + 1):
        n = izs[i]
        ws[n] += 1.0
        k = itab[i] - 1
        xs[n] += S[k, 0]
        ys[n] += S[k, 2]
        zs[n] += S[k, 4]
    for n in range(1, ns1 + 1):
        if ws[n] != 0.0:
            zs[n] = zs[n] / ws[n]
    zs[0] = S[itab[1] - 1, 4] - dzwr / np_
    zs[ns] = S[itab[np_] - 1, 4] + dzwr / np_
    if lwl > 0:
        zw0l = min(0.0, wakel[0, 0])
    else:
        zw0l = 1e300
    if lwt > 0:
        zw0t = min(0.0, waket[0, 0])
    else:
        zw0t = 1e300
    for n in range(0, ns + 1):
        for m in range(0, ns + 1):
            dz = zs[m] - zs[n]
            if dz >= zw0l:
                for l in range(1, lwl + 1):
                    if wakel[l - 1, 0] >= dz:
                        if n != m or zw0l < 0.0:
                            l1 = min(max(l - 1, 1), lwl - 1)
                            w = (dz - wakel[l1 - 1, 0]) / (wakel[l1, 0] - wakel[l1 - 1, 0])
                            wz[n] += (wakel[l1 - 1, 1] * (1.0 - w) + wakel[l1, 1] * w) * ws[m]
                        else:
                            wz[n] += wakel[0, 1] * ws[n] * 0.5
                        break
            if dz >= zw0t:
                for l in range(1, lwt + 1):
                    if waket[l - 1, 0] >= dz:
                        if n != m or zw0t < 0.0:
                            l1 = min(max(l - 1, 1), lwt - 1)
                            w = (dz - waket[l1 - 1, 0]) / (waket[l1, 0] - waket[l1 - 1, 0])
                            wx[n] += (waket[l1 - 1, 1] * (1.0 - w) + waket[l1, 1] * w) * xs[m]
                            wy[n] += (waket[l1 - 1, 1] * (1.0 - w) + waket[l1, 1] * w) * ys[m]
                        else:
                            wx[n] += waket[0, 1] * xs[n] * 0.5
                            wy[n] += waket[0, 1] * ys[n] * 0.5
                        break
    for i in range(1, np_ + 1):
        n = izs[i]
        k = itab[i] - 1
        zk = S[k, 4]
        dzk = zk - zs[n]
        if dzk < 0.0:
            n -= 1
            dzk = zk - zs[n]
        w = dzk / (zs[n + 1] - zs[n])
        dwx = (1.0 - w) * wx[n] + w * wx[n + 1]
        dwy = (1.0 - w) * wy[n] + w * wy[n + 1]
        dwz = (1.0 - w) * wz[n] + w * wz[n + 1]
        pa = 1.0 + S[k, 5]
        S[k, 5] = max(S[k, 5] - fwp * dwz, PMIN - 1.0)
        pb = 1.0 + S[k, 5]
        S[k, 1] = (pa * S[k, 1] + fwp * dwx) / pb
        S[k, 3] = (pa * S[k, 3] + fwp * dwy) / pb
        h1 = math.sqrt(1.0 + (pa * p0) ** 2)
        S[k, 6] = -S[k, 5] * pa / h1 / (h1 + pa * h0)


def _rot(S, n, cost, sint):
    """TENT.inc / TEXIT.inc の回転 (sint の符号だけ逆)。x, y, px, py を回す。|cos|≥|sin| のときは Fortran と同じ代入の順"""
    x, px, y, py = S[:n, 0], S[:n, 1], S[:n, 2], S[:n, 3]
    if abs(cost) >= abs(sint):
        x[:] = cost * x - sint * y
        y[:] = (sint * x + y) / cost
        px[:] = cost * px - sint * py
        py[:] = (sint * px + py) / cost
    else:
        xi, pxi = x.copy(), px.copy()
        x[:] = xi * cost - y * sint
        y[:] = xi * sint + y * cost
        px[:] = pxi * cost - py * sint
        py[:] = pxi * sint + py * cost


class WakeKick:
    """追跡 1 回 (tturn 1 区間) の間の txwake の状態 (itab、izs、dzwr)。区間の始めに reset() する (tturn.f:200 の dzwr=0)"""

    def __init__(self, p0, h0, fw, nb):
        self.p0, self.h0, self.fw, self.nb = p0, h0, fw, max(int(nb), 1)
        self.itab = None
        self.izs = None
        self.dzwr = 0.0
        self.np = -1

    def reset(self):
        self.dzwr = 0.0
        self.np = -1

    def apply(self, S, n, wl, wt, dx, dy, init, p0=None, h0=None, theta=0.0):
        """wl, wt: (行数, 2) の配列か None (そのウェイクなし)。dx, dy: 要素のずれ (TENT/TEXIT)"""
        # twxiwp (tffs.f:1686): lwl0 = (ブロック長-2)/2 は表の行数より 1 少ない (ブロック長 = 2 行数 + 1)。追跡では最後の行を使わない (gs で確認)
        wl = None if wl is None else wl[:-1]
        wt = None if wt is None else wt[:-1]
        lwl = 0 if wl is None else len(wl)
        lwt = 0 if wt is None else len(wt)
        if (lwl == 0 and lwt == 0) or n <= 0:
            return
        if lwl > 0 and lwl < 2 or lwt > 0 and lwt < 2:
            raise ValueError('WakeFunction: 追跡に使う行 (最後の 1 行を除く) が 2 行以上の表が必要')
        wl_ = np.ascontiguousarray(wl, dtype=float) if lwl else np.zeros((1, 2))
        wt_ = np.ascontiguousarray(wt, dtype=float) if lwt else np.zeros((1, 2))
        cost = sint = 0.0
        if dx != 0.0:
            S[:n, 0] -= dx
        if dy != 0.0:
            S[:n, 2] -= dy
        if theta != 0.0:                          # inc/TENT.inc
            cost, sint = math.cos(theta), math.sin(theta)
            _rot(S, n, cost, sint)
        z = S[:n, 4]
        resort = init or self.itab is None or self.np != n
        if not resort:
            resort = abs(z[self.itab[n] - 1] - z[self.itab[1] - 1] - self.dzwr) >= DZLIM * self.dzwr
        if resort:
            bad = np.isnan(z)
            if bad.any():
                z[bad] = ZMIN
            order = _spsort(np.ascontiguousarray(z))      # spkick.f の spsort (同じ z の粒子の並びもそれに合わせる)
            self.itab = np.zeros(n + 1, dtype=np.int64)
            self.itab[1:] = order + 1
            zs = np.zeros(n + 1)
            zs[1:] = z[order]
            self.izs = np.zeros(n + 1, dtype=np.int64)
            _defslice(n, zs, self.nb, self.izs)
            self.dzwr = float(zs[n] - zs[1])
            self.np = n
        _kick(S, n, self.itab, self.izs, self.dzwr, self.fw, self.p0 if p0 is None else p0,
              self.h0 if h0 is None else h0, wl_, wt_, lwl, lwt)
        if theta != 0.0:                          # inc/TEXIT.inc
            _rot(S, n, cost, -sint)
        if dy != 0.0:
            S[:n, 2] += dy
        if dx != 0.0:
            S[:n, 0] += dx


def kick_factor(charge, pbunch, nbunch, np0, amass=ELMASS):
    """fw = (|charge| e pbunch nbunch / amass) / np0 * .5 (tturn.f:259)"""
    return abs(charge) * ELEMCH * pbunch * nbunch / amass / np0 * 0.5
