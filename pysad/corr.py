"""
軌道補正の組み込み関数: ResponseMatrix (K64 の pgrmat.f)、LinearSolveConditioned (pgsolvcond.f、pmicad.f)。
Packages/correction.n・bump.n・corutility.n (CorrectOrbit、Monitor、Steer 等の SADScript) がこれらを呼ぶ。

ResponseMatrix[mon, st, stt, opt, periodic, kind (, tune)]: mon (要素の位置。リストなら 1 要素ずつ)、st・stt (補正磁石の始端・終端の位置のリストまたは実数)、
opt (CalculateOptics の戻り)、periodic (リングは 1)、kind ("X" "Y" "PX" "PY" と "MX" "MY" "MPX" "MPY")、tune ({νx, νy} の目標値。0 は opt のまま)。
未対応: TRPT (加速)。
"""
from __future__ import annotations

import math

import numpy as np

from . import compat as _compat
from .sadscript import SADError, is_num, S, Sym, Expr, Null
from .tsvd import tsvd

# twiss 配列の添字 (tffitcode の mfit*。0 始まり)
AX, BX, NX, AY, BY, NY, EX, EPX, EY, EPY, R1, R2, R3, R4, DETR = range(15)
_PI2 = 2.0 * math.pi


def _sign1(n):
    return 1.0 if n >= 0 else -1.0


def response_row(im, st, stt, rows, rot, periodic, keyword, psix_t=0.0, psiy_t=0.0, bend_idx=None):
    """pgrmat1 の 1 モニタ分。im は 1 始まり。st・stt は 1 始まりの位置の配列。rows[k-1] は要素 k の optics (28 要素)。
    rot(js): 要素 js が BEND なら -ROTATE、そうでなければ None。戻り: 長さ len(st) の配列"""
    kw = keyword.upper()
    normalmode = kw[0] == 'M'
    if kw[0] == 'P':
        angle, xplane = True, kw == 'PX'
    elif kw[0] == 'M':
        xplane = kw in ('MX', 'MPX')
        angle = kw in ('MPX', 'MPY')
    else:
        xplane = kw == 'X'
        angle = False
    n = len(rows)
    first, last = rows[0], rows[n - 1]
    psix0 = last[NX] - first[NX]
    psiy0 = last[NY] - first[NY]
    psix = psix_t * _PI2
    psiy = psiy_t * _PI2
    if psix == 0.0:
        psix = psix0
    if psiy == 0.0:
        psiy = psiy0
    twsi = rows[im - 1]
    axi, bxi, ayi, byi = twsi[AX], twsi[BX], twsi[AY], twsi[BY]
    ri11, ri12, ri21, ri22 = twsi[R1], twsi[R2], twsi[R3], twsi[R4]
    detri = twsi[DETR]
    pinx = 0.5 * psix
    piny = 0.5 * psiy
    if periodic:
        fx = 0.5 / math.sin(pinx)
        fy = 0.5 / math.sin(piny)
    ar = np.zeros(len(st))
    for i in range(len(st)):
        js = int(st[i])
        rr = rot(js)
        if rr is not None:
            cost = math.cos(rr)
            sint = math.sin(rr)
        elif xplane:
            cost, sint = 1.0, 0.0
        else:
            cost, sint = 0.0, 1.0
        tot = 0.0
        for i1 in (1, 2):
            if i1 == 2:
                js = int(stt[i]) + 1
            entrance = (im == js and i1 == 1)
            if js < 1 or js > n:
                continue
            kls = rows[js - 1]
            r11, r12, r21, r22 = kls[R1], kls[R2], kls[R3], kls[R4]
            detr = kls[DETR]
            det = r11 * r22 - r12 * r21
            um = math.sqrt(1.0 - det) if det <= 1.0 else float('nan')
            if detr > 1.0:
                tx = r12 * cost
                tpx = -r11 * cost + um * sint
                ty = r12 * sint
                tpy = um * cost + r22 * sint
            else:
                tx = r12 * sint
                tpx = um * cost - r11 * sint
                ty = r12 * cost
                tpy = r22 * cost + um * sint
            axj, bxj = kls[AX], kls[BX]
            dphix = twsi[NX] - kls[NX]
            ayj, byj = kls[AY], kls[BY]
            dphiy = twsi[NY] - kls[NY]
            if periodic:
                phx = pinx - abs(dphix) * psix / psix0
                phy = piny - abs(dphiy) * psiy / psiy0
                sx, cx, sy, cy = math.sin(phx), math.cos(phx), math.sin(phy), math.cos(phy)
                sg = -1.0 if entrance else _sign1(im - js)
            else:
                sx, cx, sy, cy = math.sin(dphix), math.cos(dphix), math.sin(dphiy), math.cos(dphiy)
                sg = None
                if im < js or entrance:
                    continue

            def xfun():
                if periodic:
                    return fx * (tx * math.sqrt(bxi / bxj) * (sg * sx + axj * cx) + tpx * math.sqrt(bxi * bxj) * cx)
                return tx * math.sqrt(bxi / bxj) * (cx + axj * sx) + tpx * math.sqrt(bxi * bxj) * sx

            def pxfun():
                if periodic:
                    return fx * (-tx / math.sqrt(bxi * bxj) * ((axi - axj) * sg * sx + (1.0 + axi * axj) * cx)
                                 + tpx * math.sqrt(bxj / bxi) * (sg * sx - axi * cx))
                return (-tx / math.sqrt(bxi * bxj) * ((axi - axj) * cx + (1.0 + axi * axj) * sx)
                        + tpx * math.sqrt(bxj / bxi) * (cx - axi * sx))

            def yfun():
                if periodic:
                    return fy * (ty * math.sqrt(byi / byj) * (sg * sy + ayj * cy) + tpy * math.sqrt(byi * byj) * cy)
                return ty * math.sqrt(byi / byj) * (cy + ayj * sy) + tpy * math.sqrt(byi * byj) * sy

            def pyfun():
                if periodic:
                    return fy * (-ty / math.sqrt(byi * byj) * ((ayi - ayj) * sg * sy + (1.0 + ayi * ayj) * cy)
                                 + tpy * math.sqrt(byj / byi) * (sg * sy - ayi * cy))
                return (-ty / math.sqrt(byi * byj) * ((ayi - ayj) * cy + (1.0 + ayi * ayj) * sy)
                        + tpy * math.sqrt(byj / byi) * (cy - ayi * sy))
            if normalmode:
                if angle:
                    v = pxfun() if xplane else pyfun()
                else:
                    v = xfun() if xplane else yfun()
            else:
                x, px, y, py = xfun(), pxfun(), yfun(), pyfun()
                det = ri11 * ri22 - ri12 * ri21
                um = math.sqrt(1.0 - det) if det <= 1.0 else float('nan')
                if detri > 1.0:
                    if angle:
                        v = (-ri21 * x + ri11 * px) + um * py if xplane else um * px + (-ri21 * y - ri22 * py)
                    else:
                        v = (ri22 * x - ri12 * px) + um * y if xplane else um * x + (-ri12 * y - ri22 * py)
                else:
                    if angle:
                        v = um * px - ri21 * y + ri11 * py if xplane else -ri21 * x - ri22 * px + um * py
                    else:
                        v = um * x + ri22 * y - ri12 * py if xplane else -ri11 * x - ri12 * px + um * y
            tot += v
        ar[i] = -0.5 * tot
    return ar


def psvdcond(a, b, c, d, eps):
    """pgsolvcond.f psvdcond: 条件 c.x=d の下での最小二乗。a, b, c は破壊されてよいコピー"""
    epsc = 1e-13
    nc, m = c.shape
    c1 = c.copy()
    u = np.zeros(nc)
    for l in range(nc):
        t = np.array([float(np.dot(c[l], c[i])) / u[i] for i in range(l)])
        if l:
            c[l] = c[l] - t @ c[:l]
        s = float(np.dot(c[l], c[l]))
        if s < m * epsc ** 2:
            c[l] = 0.0
            u[l] = 1.0
        else:
            u[l] = s
    xc = tsvd(c1, np.array(d, dtype=float), 1e-8)
    b = b - a @ xc
    for l in range(a.shape[0]):
        t = np.array([float(np.dot(a[l], c[i])) / u[i] for i in range(nc)])
        a[l] = a[l] - t @ c
    x = tsvd(a, b, eps)
    return x + xc


def pmicad(a, b, nx0, norm, eps, solvf=True):
    """pmicad.f の直訳 (ベストコレクタ法)。a は (m, n) で破壊される。戻り x (n)"""
    prec = 1e-6
    m, n = a.shape
    nx = min(int(nx0), m, n)
    iw = list(range(1, n + 1))
    tm = 0.0
    small = 0.0
    k = 0
    sk = 0.0
    if norm:
        tm = max(float(np.sum(a[:, j] ** 2)) for j in range(n))
    for l in range(1, nx + 1):
        if norm:
            small = tm * prec ** 2
        sm = 0.0
        tm = 0.0
        for j in range(l, n + 1):
            col = a[l - 1:m, j - 1]
            s = float(np.dot(col, b[l - 1:m]))
            sd = float(np.dot(col, col))
            s = s * s
            if sd > tm:
                tm = sd
            if norm:
                s = s / (sd + small)
            if s >= sm:
                sm = s
                k = j
                sk = sd
        iw[l - 1], iw[k - 1] = iw[k - 1], iw[l - 1]
        a[:, [l - 1, k - 1]] = a[:, [k - 1, l - 1]]
        sk = math.copysign(math.sqrt(sk), a[l - 1, l - 1])
        h = 1.0 / (sk * (sk + a[l - 1, l - 1]))
        a[l - 1, l - 1] += sk
        al = a[l - 1:m, l - 1].copy()
        for j in range(l + 1, n + 1):
            s = float(np.dot(al, a[l - 1:m, j - 1])) * h
            a[l - 1:m, j - 1] -= al * s
        s = float(np.dot(al, b[l - 1:m])) * h
        b[l - 1:m] -= al * s
        a[l - 1, l - 1] = -sk
        a[l:m, l - 1] = 0.0
    x = np.ones(n)
    if solvf:
        if _compat.v2():
            # v2: 選んだ nx 本の補正量を、ハウスホルダー変換後の上三角部分 a[:nx,:nx]・b[:nx] だけで解く (旧 `tsolva(a,b,x,nx,nx,mdim,eps)`)。
            # K64 は新しい tsvd (配列の大きさを size() で取る) に替えたとき nx を渡し忘れ、全 n 本の最小二乗解を解いて先頭の nx 個だけ使う
            # (選ばなかった補正量の寄与が落ちるので、軌道がむしろ悪化する)
            x = np.zeros(n)
            x[:nx] = tsvd(a[:nx, :nx].copy(), b[:nx].copy(), eps)
        else:
            x = tsvd(a, b, eps)
    x[nx:] = 0.0
    if not solvf:
        x[:nx] = 1.0
    # msortn(iw, x, b, nx, 1, nx): iw を昇順に並べ替え、x (と b) を同じ順に動かす
    order = sorted(range(nx), key=lambda i: iw[i])
    iwn = [iw[i] for i in order]
    xn = x[:nx][order].copy()
    x[:nx] = xn
    iw[:nx] = iwn
    for i in range(nx, 0, -1):
        if iw[i - 1] == i:
            return x
        x[iw[i - 1] - 1] = x[i - 1]
        x[i - 1] = 0.0
    return x


def linear_solve_conditioned(a, b, eps=1e-8, nx=0, norm=True, c=None, d=None):
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    n, m = a.shape
    if len(b) != n:
        raise SADError('General::equalleng: #1 and #2')
    cond = c is not None
    if cond:
        c = np.array(c, dtype=float).reshape(-1, m)
        d = np.array(d, dtype=float)
        if c.shape[0] != len(d):
            raise SADError('General::equalleng: matrices')
    micado = int(nx) != 0
    if cond:
        if micado:
            aa, bb, cc = a.copy(), b.copy(), c.copy()
            x = pmicad(a, b, nx, norm, eps, False)
            sel = [i for i in range(m) if x[i] != 0.0]
            a1 = aa[:, sel]
            c1 = cc[:, sel]
            xs = psvdcond(a1.copy(), bb.copy(), c1.copy(), d, eps)
            x = np.zeros(m)
            for k, j in enumerate(sel):
                x[j] = xs[k]
            return x
        return psvdcond(a.copy(), b.copy(), c.copy(), d, eps)
    if micado:
        return pmicad(a, b, nx, norm, eps, True)
    return tsvd(a, b, eps)


def install(I, line_fn=None):
    B = I.builtins

    def _lsc(a):
        if len(a) < 2 or len(a) > 7 or len(a) == 6:
            raise SADError('General::narg: Number of arguments is expected to be 2, 3, 4, 5, or 7 in LinearSolveConditioned')
        eps = float(a[2]) if len(a) >= 3 else 1e-8
        nx = float(a[3]) if len(a) >= 4 else 0
        norm = bool(a[4]) if len(a) >= 5 else True
        c = a[5] if len(a) == 7 else None
        d = a[6] if len(a) == 7 else None
        return [float(v) for v in linear_solve_conditioned(a[0], a[1], eps, nx, norm, c, d)]
    B['LinearSolveConditioned'] = _lsc
