"""WAKEOPT: 光学計算 (CALC) での wake (tffswake.f)。

InitialOrbits の各軌道 (と中心の軌道) を INS で伝搬し、WakeFunction を持つ要素の前後で、全軌道の z・x・y から
各軌道の dpx, dpy, ddp へのキックを与える。軌道は FFS の idp (nfam1..nfam、中心 0) で区別する。

K64 の tffswakekick の直訳。縦 (lwake) は dz = z_i − z_j の表引き (z_i − z_j ≥ −sigg の滑らかな境界)、
横 (twake) は dz + 0.3989 sigg の表引き (wbufxy、wbufxyl は表が同じで z の配置が (ほぼ) 同じなら前の値を使い回す)。
キックの前後で同じ表のとき使い回す条件も同じ (dsig = 0.1 sigg)。
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np

from . import lattice as _lat
from .emit6d import ELEMCH

_FACT = 0.398942280401433            # 1/sqrt(2 pi)
_TW_KEYS = ('ax', 'bx', 'nx', 'ay', 'by', 'ny', 'ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4', 'detr')


def _bracket(tab_z: np.ndarray, dz: float) -> int:
    """tab_z[i1] <= dz < tab_z[i1+1] の i1 (tffswakekick の二分探索と同じ区間。表は z の昇順)"""
    i1 = int(np.searchsorted(tab_z, dz, side='right')) - 1
    return min(max(i1, 0), len(tab_z) - 2)


class _Kicker:
    """tffswakekick の状態 (wbuf 系の配列と、前回の表・z) を持つ。1 回の CALC (tffswake) の間だけ生きる"""

    def __init__(self, nslot: int, sigz: float, nbunch: float, charge: float, pbunch: float, lwake: bool, twake: bool):
        self.np = nslot
        self.shift = False                    # nfam1 = 1−nfam (軌道の数が奇数) のとき True (縦のキックの添字のずれ)
        self.sigz = sigz
        self.nbunch = nbunch
        self.charge = charge
        self.pbunch = pbunch
        self.lwake = lwake
        self.twake = twake
        self.key_l = None                     # iwsl
        self.key_t = None                     # iwst
        self.wzl = np.zeros(nslot)
        self.wzt = np.zeros(nslot)
        self.wbuf = np.zeros(nslot)
        self.wbufxy = np.zeros((nslot, nslot))
        self.wbufxyl = np.zeros((nslot, nslot))

    def kick(self, ut: np.ndarray, wl: Optional[np.ndarray], wt: Optional[np.ndarray], dx: float, dy: float, p: float,
             key_l, key_t) -> np.ndarray:
        """ut: (nslot, 6) = (x, px, y, py, z, dp) の各軌道。キックを加えたコピーを返す。p = 2 gammab amass [eV]"""
        n = self.np
        out = ut.copy()
        cp0 = abs(self.charge * ELEMCH * self.pbunch / n / p)
        sigg = max(1e-5, self.sigz / n / self.nbunch)
        dsig = 0.1 * sigg
        z = ut[:, 4]
        if wl is not None and self.lwake:
            if key_l == self.key_l:
                ddz0 = z[0] - self.wzl[0]
                re = bool(np.any(np.abs(z[1:] - self.wzl[1:] - ddz0) > dsig))
            else:
                re = True
            if re:
                self.key_l = key_l
                wbuf = np.zeros(n)
                wz = wl[:, 0]
                wv = wl[:, 1]
                nlw = len(wl)
                for i in range(n):
                    zi = z[i]
                    self.wzl[i] = zi
                    for j in range(n):
                        dz = zi - z[j]
                        if dz >= -sigg:
                            w = min(1.0, (dz + sigg) / sigg * 0.5)
                            if dz <= wz[0]:
                                wbuf[j] += w * wv[0]
                            elif dz < wz[nlw - 1]:
                                i1 = _bracket(wz, dz)
                                a = (dz - wz[i1]) / (wz[i1 + 1] - wz[i1])
                                wbuf[j] += w * ((1.0 - a) * wv[i1] + a * wv[i1 + 1])
                self.wbuf = wbuf
            wb = self.wbuf
            if self.shift:
                # K64 の癖: ut1(mfitddp,:)=ut0(mfitddp,:)-cp0*wbuf は配列全体 (−nfam..nfam) への代入。wbuf は nfam1..nfam の np 個なので、
                # nfam1=1−nfam (軌道の数 n が奇数) のとき要素の対応が 1 つずれる: 軌道 s は 1 つ後ろの軌道 s+1 のキックを受け、
                # 最後の軌道は wbuf の外 (未定義の値、実測では無視できる大きさ) を読む。ここでは 0 とする
                wb = np.append(self.wbuf[1:], 0.0)
            out[:, 5] = ut[:, 5] - cp0 * wb
        if wt is not None and self.twake:
            re = False
            if key_t == self.key_t:
                ddz0 = z[0] - self.wzt[0]
                re = not bool(np.any(np.abs(z[1:] - self.wzt[1:] - ddz0) > dsig))
            wbufx = np.zeros(n)
            wbufy = np.zeros(n)
            wbufxl = np.zeros(n)
            wbufyl = np.zeros(n)
            xs = ut[:, 0] - dx
            ys = ut[:, 2] - dy
            if re:
                for i in range(n):
                    wbufx += xs[i] * self.wbufxy[i, :]
                    wbufy += ys[i] * self.wbufxy[i, :]
                    wbufxl -= xs[i] * self.wbufxyl[i, :]
                    wbufyl -= ys[i] * self.wbufxyl[i, :]
            else:
                self.key_t = key_t
                wz = wt[:, 0]
                wv = wt[:, 1]
                ntw = len(wt)
                dz0 = sigg * _FACT
                for i in range(n):
                    zi = z[i]
                    self.wzt[i] = zi
                    for j in range(n):
                        dz = zi - z[j] + dz0
                        if wz[0] < dz < wz[ntw - 1]:
                            i1 = _bracket(wz, dz)
                            d = (wv[i1 + 1] - wv[i1]) / (wz[i1 + 1] - wz[i1])
                            w = wv[i1] + d * (dz - wz[i1])
                            wbufx[j] += xs[i] * w
                            wbufy[j] += ys[i] * w
                            wbufxl[j] -= xs[i] * d
                            wbufyl[j] -= ys[i] * d
                            self.wbufxy[i, j] = w
                            self.wbufxyl[i, j] = d
                        else:
                            self.wbufxy[i, j] = 0.0
                            self.wbufxyl[i, j] = 0.0
            out[:, 1] = ut[:, 1] + cp0 * wbufx
            out[:, 3] = ut[:, 3] + cp0 * wbufy
            out[:, 5] = out[:, 5] + cp0 * (wbufxl * xs + wbufyl * ys)
        return out


def _state_cod(st):
    return st[1]


def _propagate(els, lb: int, le: int, st, rows_tw: dict):
    """位置 lb から le まで (要素 lb..le−1) を INS で伝搬 (qcell1(ibound, ii, ...))。状態 st=(tw 辞書, cod[6]) を返し、
    rows_tw[位置] に各位置の (tw, cod) を書く (位置 lb は書き換え)"""
    tw, cod = st
    if le == lb:
        rows_tw[lb] = dict(tw, cod=np.array(cod, dtype=float))
        return st
    tws = _lat.ins_twiss(els[lb - 1:le - 1], tw, cod[:5], float(cod[5]))
    for k, t in enumerate(tws):
        rows_tw[lb + k] = t
    t = tws[-1]
    return ({k: t[k] for k in _TW_KEYS if k in t}, np.array(t['cod'], dtype=float))


def wake_orbits(els, slots: Dict[int, tuple], tables: Dict[int, tuple], *, sigz: float, nbunch: float, charge: float,
                pbunch: float, lwake: bool, twake: bool, gammab, amass: float) -> Dict[int, Dict[int, dict]]:
    """tffswake。slots: {idp: (tw, cod[6])} (昇順に並べて使う)、tables: {要素番号: (縦の表, 横の表, dx, dy)}、
    gammab(k): 位置 k の基準運動量 (p/mc)。戻り: {idp: {位置: 状態の辞書 (cod 付き)}} (位置 1..nlat)"""
    nlat = len(els) + 1
    ids = sorted(slots)
    n = len(ids)
    kk = _Kicker(n, sigz, nbunch, charge, pbunch, lwake, twake)
    kk.shift = ids[0] != -ids[-1]
    state = {i: slots[i] for i in ids}
    rows = {i: {} for i in ids}
    keys: Dict[int, tuple] = {}

    def key_of(a):
        return None if a is None else (a.shape, a.tobytes())

    lb = 1
    for le in sorted(k for k in tables if 1 <= k < nlat):
        if le < lb:
            continue
        wl, wt, dx, dy = tables[le]
        # lb..le の伝搬 (wake 要素の入口まで)
        for i in ids:
            state[i] = _propagate(els, lb, le, state[i], rows[i])
        ut = np.array([state[i][1] for i in ids])
        ut1 = kk.kick(ut, wl, wt, dx, dy, 2.0 * gammab(le) * amass, key_of(wl), key_of(wt))
        for a, i in enumerate(ids):
            state[i] = (state[i][0], ut1[a])
            rows[i][le] = dict(state[i][0], cod=ut1[a].copy())
        # 要素 le の伝搬
        for i in ids:
            state[i] = _propagate(els, le, le + 1, state[i], rows[i])
        ut = np.array([state[i][1] for i in ids])
        ut2 = kk.kick(ut, wl, wt, dx, dy, 2.0 * gammab(le + 1) * amass, key_of(wl), key_of(wt))
        for a, i in enumerate(ids):
            state[i] = (state[i][0], ut2[a])
            rows[i][le + 1] = dict(state[i][0], cod=ut2[a].copy())
        lb = le + 1
    for i in ids:
        if lb < nlat:
            state[i] = _propagate(els, lb, nlat, state[i], rows[i])
        elif lb == nlat and nlat not in rows[i]:
            pass
    return rows
