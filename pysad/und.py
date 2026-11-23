"""
UND 要素 (アンジュレータ、K64 の undulator.f: K.Ohmi と E.Forest、KEK Report 92-14)。追跡だけ。

K64 では UND に磁場のキーワードが無い: initb1.f は `kytbl(kwFBX,icUND)=ky_FBX_UND` と書くが kwFBX は未定義の変数 (=0) で、
`FBX` `FBY` `BX` `BY` はどれも "Undefined keyword" になる。そのため F0 = G0 = 0 で、undulator は磁場なしの (正準運動量の
近似の) ドリフトを刻みごとに行う。Twiss・Emittance[] の写像は無い (恒等、長さは幾何にだけ入る)。
直訳は F0, G0 を保ったままにしてあるが、利用者は設定できない。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from .track6n import tradk_core1

NJ = dict(cache=True, error_model='numpy')
_BAD = ('FBX', 'FBY', 'BX', 'BY')
_C = 299792458.0


def compile_und(el, pgev):
    """undinit。戻り値は k_und の引数 (nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds)"""
    p = el.params
    g = lambda k: float(p.get(k, 0.0) or 0.0)
    from . import compat
    if compat.v2():
        # v2 (B37 の修正): FBX・FBY を受ける (undinit: F0 = p_in(2)/pgev*c、G0 = p_in(3)/pgev*c。単位は T、F0・G0 は 1/m)。
        # BX・BY は K64 の initb1.f にも無い名前なので受けない
        for k in ('BX', 'BY'):
            if k in p:
                raise ValueError(f'UND: キーワード {k} は無い (磁場は FBX、FBY [T])')
        f0 = g('FBX') / pgev * _C
        g0 = g('FBY') / pgev * _C
    else:
        for k in _BAD:
            if k in p:
                raise NotImplementedError(f'UND: キーワード {k} は K64 に無い (Undefined keyword)')
        f0 = 0.0 * pgev            # p_in(2)/pgev*c: FBX が設定できないので 0
        g0 = 0.0
    length = float(el.L)
    nsli = g('SLICE')
    npole = g('POLE')
    dph = g('PHI') * math.pi / 180.0
    with np.errstate(all='ignore'):
        kz = 2.0 * math.pi * npole / length if length != 0.0 else float('nan')
        kx = g('KX')
        ky = math.sqrt(kz * kz + kx * kx)
        qy = g('QY')
        qx = math.sqrt(kz * kz + qy * qy)
        ds = length / nsli if nsli != 0.0 else float('inf')
    return (int(nsli + 0.01), f0, g0, dph, kz, kx, ky, qx, qy, ds)


@njit(**NJ)
def k_und(S, n, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, p0, h0, b2=False):
    for i in range(n):
        pn = 1.0 + S[i, 5]
        S[i, 1] = S[i, 1] * pn
        S[i, 3] = S[i, 3] * pn
    for i in range(n):
        x = S[i, 0]
        px = S[i, 1]
        y = S[i, 2]
        py = S[i, 3]
        z = S[i, 4]
        g = S[i, 5]
        for j in range(1, nslice + 1):
            kx1 = Kx * x
            kx2 = Ky * y
            kx3 = Qx * x
            kx4 = Qy * y
            sinkz = math.sin(Kz * j * ds)
            sinkzp = math.sin(Kz * j * ds + dph)
            ff = (F0 / Kz * math.cos(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Qx / Kz * math.sinh(kx3) * math.sin(kx4) * sinkzp)
            gg = (-G0 / Kz * math.cos(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Ky / Kz * math.sinh(kx2) * math.sin(kx1) * sinkz)
            fx = (-Kx * F0 / Kz * math.sin(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Kz * math.cosh(kx3) * math.sin(kx4) * sinkzp)
            fy = (Ky * F0 / Kz * math.cos(kx1) * math.sinh(kx2) * sinkz
                  - Qy * Qy * G0 / Qx / Kz * math.sinh(kx3) * math.cos(kx4) * sinkzp)
            gx = (-G0 * Qx / Kz * math.cos(kx4) * math.sinh(kx3) * sinkzp
                  + Kx * Kx * F0 / Ky / Kz * math.sinh(kx2) * math.cos(kx1) * sinkz)
            gy = (G0 * Qy / Kz * math.sin(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Kz * math.cosh(kx2) * math.sin(kx1) * sinkz)
            dpds = ds / (1.0 + g)
            fx = dpds * fx
            gx = dpds * gx
            fy = dpds * fy
            gy = dpds * gy
            det1a = 1.0 - fx - gy + fx * gy - fy * gx
            fgx = dpds * (ff * fx + gg * gx)
            fgy = dpds * (ff * fy + gg * gy)
            tmp = fx
            fx = (1.0 - gy) / det1a
            gx = gx / det1a
            fy = fy / det1a
            gy = (1.0 - tmp) / det1a
            tmp = px - fgx
            px = fx * tmp + gx * (py - fgy)
            py = fy * tmp + gy * (py - fgy)
            fgx = px - ff
            fgy = py - gg
            x = x + dpds * fgx
            y = y + dpds * fgy
            z = z - 0.5 * (fgx * fgx + fgy * fgy) * (dpds * dpds) / ds
        S[i, 0] = x
        S[i, 1] = px
        S[i, 2] = y
        S[i, 3] = py
        S[i, 4] = z
    for i in range(n):
        pn = 1.0 + S[i, 5]
        S[i, 1] = S[i, 1] / pn
        S[i, 3] = S[i, 3] / pn
        h1 = math.sqrt(1.0 + (pn * p0)**2)
        S[i, 6] = -S[i, 5] * (1.0 + pn) / h1 / (h1 + pn * h0) + S[i, 10]
        if b2:
            S[i, 4] -= S[i, 6] * (nslice * ds)       # v2 (B334): 速度の差による z の更新 −dv·L (K64 の undulator.f には無い。DRIFT は持つ)


@njit(**NJ)
def _und_a(x, y, sj, F0, G0, dph, Kz, Kx, Ky, Qx, Qy):
    """ベクトルポテンシャル (ax, ay) = (ff, gg) (undulator.f と同じ式)"""
    sinkz = math.sin(Kz * sj)
    sinkzp = math.sin(Kz * sj + dph)
    ff = (F0 / Kz * math.cos(Kx * x) * math.cosh(Ky * y) * sinkz
          - Qy * G0 / Qx / Kz * math.sinh(Qx * x) * math.sin(Qy * y) * sinkzp)
    gg = (-G0 / Kz * math.cos(Qy * y) * math.cosh(Qx * x) * sinkzp
          + Kx * F0 / Ky / Kz * math.sinh(Ky * y) * math.sin(Kx * x) * sinkz)
    return ff, gg


@njit(**NJ)
def und_seg1(i, S, j0, j1, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, rad, fluc, rs, p0, anrad, cuc):
    """1 粒子を刻み j0+1..j1 だけ進める。S[i, 1], S[i, 3] は出入りとも力学的運動量の向き (px/p)。
    アンジュレータの両端 (j0=0、j1=nslice) では undulator.f と同じく正準運動量を外の運動量につなぐ (端でのポテンシャルの飛びは蹴りになる)"""
    if True:
        g = S[i, 5]
        pn = 1.0 + g
        x = S[i, 0]
        y = S[i, 2]
        z = S[i, 4]
        vx = S[i, 1]
        vy = S[i, 3]
        a0x, a0y = _und_a(x, y, j0 * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
        if j0 == 0:
            px = vx * pn
            py = vy * pn
            vx = (px - a0x) / pn             # 入口の向き (力学的運動量 / p)
            vy = (py - a0y) / pn
        else:
            px = vx * pn + a0x
            py = vy * pn + a0y
        for j in range(j0 + 1, j1 + 1):
            kx1 = Kx * x
            kx2 = Ky * y
            kx3 = Qx * x
            kx4 = Qy * y
            sinkz = math.sin(Kz * j * ds)
            sinkzp = math.sin(Kz * j * ds + dph)
            ff = (F0 / Kz * math.cos(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Qx / Kz * math.sinh(kx3) * math.sin(kx4) * sinkzp)
            gg = (-G0 / Kz * math.cos(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Ky / Kz * math.sinh(kx2) * math.sin(kx1) * sinkz)
            fx = (-Kx * F0 / Kz * math.sin(kx1) * math.cosh(kx2) * sinkz
                  - Qy * G0 / Kz * math.cosh(kx3) * math.sin(kx4) * sinkzp)
            fy = (Ky * F0 / Kz * math.cos(kx1) * math.sinh(kx2) * sinkz
                  - Qy * Qy * G0 / Qx / Kz * math.sinh(kx3) * math.cos(kx4) * sinkzp)
            gx = (-G0 * Qx / Kz * math.cos(kx4) * math.sinh(kx3) * sinkzp
                  + Kx * Kx * F0 / Ky / Kz * math.sinh(kx2) * math.cos(kx1) * sinkz)
            gy = (G0 * Qy / Kz * math.sin(kx4) * math.cosh(kx3) * sinkzp
                  + Kx * F0 / Kz * math.cosh(kx2) * math.sin(kx1) * sinkz)
            dpds = ds / pn
            fx = dpds * fx
            gx = dpds * gx
            fy = dpds * fy
            gy = dpds * gy
            det1a = 1.0 - fx - gy + fx * gy - fy * gx
            # 陰的な刻み p' = p - ds ∂H/∂q(q, p') の右辺。undulator.f は `fgx=dpds*(ff*fx+gg*gx)` と、すでに dpds を掛けた
            # fx.. にもう一度 dpds を掛けている (収束力が ds に比例して消える誤り)。v2 は 1 回だけ掛ける
            fgx = ff * fx + gg * gx
            fgy = ff * fy + gg * gy
            tmp = fx
            fx = (1.0 - gy) / det1a
            gx = gx / det1a
            fy = fy / det1a
            gy = (1.0 - tmp) / det1a
            tmp = px - fgx
            px = fx * tmp + gx * (py - fgy)
            py = fy * tmp + gy * (py - fgy)
            fgx = px - ff
            fgy = py - gg
            z0 = z
            x = x + dpds * fgx
            y = y + dpds * fgy
            z = z - 0.5 * (fgx * fgx + fgy * fgy) * (dpds * dpds) / ds
            if rad:
                # 刻みの前後の向き (vx, vy) → (fgx/pn, fgy/pn) から放射。刻みの終わりの位置のポテンシャルで正準運動量に戻す
                S[i, 7] = vx
                S[i, 8] = vy
                S[i, 9] = z0
                S[i, 0] = x
                S[i, 1] = fgx / pn
                S[i, 2] = y
                S[i, 3] = fgy / pn
                S[i, 4] = z
                S[i, 5] = g
                tradk_core1(i, S, ds, 0.0, fluc, rs, p0, anrad, cuc)
                x = S[i, 0]
                y = S[i, 2]
                z = S[i, 4]
                g = S[i, 5]
                pn = 1.0 + g
                vx = S[i, 1]
                vy = S[i, 3]
                a1x, a1y = _und_a(x, y, j * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
                px = vx * pn + a1x
                py = vy * pn + a1y
            else:
                vx = fgx / pn
                vy = fgy / pn
        if j1 < nslice:
            a1x, a1y = _und_a(x, y, j1 * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
            px = px - a1x
            py = py - a1y
        S[i, 0] = x
        S[i, 2] = y
        S[i, 4] = z
        S[i, 5] = g
        S[i, 1] = px / pn
        S[i, 3] = py / pn


@njit(**NJ)
def k_und_field(S, n, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, p0, h0, rad, fluc, rs, anrad, cuc, b2=False):
    """v2: 磁場のある UND (F0 か G0 が 0 でない)。写像は undulator.f の刻み (k_und と同じ式)。
    H = [(px-ax)² + (py-ay)²]/(2(1+δ))、ax = ff、ay = gg。By(軸上) = F0·cos(kz s)、Bx(軸上) = G0·cos(kz s + dph) (1/ρ の単位)。
    rad が真なら、刻みごとに運動量の向きの変化から放射 (tradk_core1: SAD の tradk と同じ式、FLUC は tdusrn) を入れる。
    K64 の undulator.f に放射は無い (v2 の追加)。"""
    for i in range(n):
        und_seg1(i, S, 0, nslice, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, rad, fluc, rs, p0, anrad, cuc)
        pn = 1.0 + S[i, 5]
        h1 = math.sqrt(1.0 + (pn * p0)**2)
        S[i, 6] = -S[i, 5] * (1.0 + pn) / h1 / (h1 + pn * h0) + S[i, 10]
        if b2:
            S[i, 4] -= S[i, 6] * (nslice * ds)       # v2 (B334): −dv·L (k_und と同じ)


# ---------------------------------------------------------------------------- 光学・エミッタンス (v2)
def annotate(els, momentum_eV):
    """v2: UND (SLICE>0) に el._solmap = (K_UND, params, momentum) を付ける (CALC の Twiss は solmap.step が追跡の写像を数値微分し、
    Emittance[] は solemit.emit → und.emit が刻みごとに放射 tradke を入れる)。v2 では光学と追跡が同じ写像になる (磁場なしでも長さ L のドリフト)。
    k64 モードでは付けない (K64 の qtwiss・tturne は UND を恒等として扱い、追跡だけがドリフトさせる)"""
    from . import compat
    from . import track6n as _n
    for el in els:
        if type(el).__name__ != 'Undulator':
            continue
        if hasattr(el, '_solmap'):
            del el._solmap
            el.__dict__.pop('transfer_matrix', None)
            el.__dict__.pop('constant_kick', None)
        if not compat.v2():
            continue
        from . import emit6d as _e
        from . import solmap as _sm
        pr = compile_und(el, _e._Phys(momentum_eV, 1.0).pgev)
        if pr[0] <= 0:
            if pr[1] != 0.0 or pr[2] != 0.0:
                raise ValueError(f'UND {el.name}: 磁場を使うには SLICE (刻み数) が要る (POLE の 20 倍以上を勧める)')
            continue                          # SLICE=0: 追跡も何もしない (恒等)
        el._solmap = (_n.K_UND, pr, momentum_eV)
        el.transfer_matrix = (lambda delta=0.0, e=el: _sm.step(e, np.array([0, 0, 0, 0, 0, delta], dtype=float))[0])
        el.constant_kick = (lambda delta=0.0: np.zeros(6))


def emit(st, pr):
    """Emittance[]・EMIT の UND (v2、磁場あり)。1 周期を 64 個以上に分けた区間ごとに、写像 (und_seg1 の数値微分) → tradke。
    tradke (photons.f の放射の式) は区間の前後の運動量の向きの変化から、エネルギー損失・減衰・量子励起を入れる"""
    from . import solemit as _se
    nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds = pr
    ph = st.ph
    rs = np.array([1], dtype=np.int64)
    per = 2.0 * math.pi / Kz / ds                      # 1 周期の刻み数
    step = max(1, int(per / 64.0)) if (F0 != 0.0 or G0 != 0.0) else nslice
    j0 = 0
    st.setr0()
    while j0 < nslice:
        j1 = min(nslice, j0 + step)
        _se._piece(st, lambda S, i, j0=j0, j1=j1: und_seg1(i, S, j0, j1, nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds,
                                                           False, False, rs, ph.p0, ph.anrad, ph.cuc))
        if j0 == 0 and G0 * math.sin(dph) != 0.0:
            # 入口の正準運動量のつなぎ (ポテンシャルの飛び) は放射させない: 基準点を区間の後の向きから作り直すのではなく、
            # 最初の区間は放射を入れず基準点だけ更新する (全体の 1/(64·POLE) 以下)
            st.setr0()
        else:
            st.tradke((j1 - j0) * ds, 0.0)
        j0 = j1
