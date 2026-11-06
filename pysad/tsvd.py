"""
tsvd.py — SAD の特異値分解ソルバ tsvd (src/tsolvg.f) の直訳 (svd=.false. の経路)

FFS のマッチング (tffsmatch.f の tfsolv → tsolva → tsvd) が使う線形最小二乗の解法。numpy の SVD (特異値の硬い打ち切り) とは
解が違う: SAD は
  (1) a を Householder/Givens で 2 重対角化し、対角要素の最小値 xmin が enorm = anorm*eps を超えれば、
      特異値を落とさずに 2 重対角系を直接解く (quick path)。条件数が 1e6〜1e8 でも打ち切らず、
      n < m (変数が条件より多い) のとき最小ノルム解とは限らない (後ろの m-n-1 成分を 0 にした厳密解)。
  (2) そうでなければ暗黙 QR 法で特異値を求め、特異値 s が anorm=max|s|*eps 以下なら硬い打ち切りでなく
      w = (s*v/anorm²)² の滑らかな減衰で解を作る。
そのためマッチングの Newton 反復の 1 歩目が (特に条件数の大きい INS のフィットで) numpy 版と大きく変わる。
(pf2 の matching_v30_12 の 2 段目で SAD の 1 歩目 K1 を再現するために移植。docs/history/work_item4.md §E)

Fortran の 1 始まりの添字をそのまま使うため、配列は 0 番目を使わない形 (a[1..n, 1..m], b[1..n], v[0..], x[1..m]) で持つ。
"""
from __future__ import annotations

import math

import numpy as np


def _sign(a: float, b: float) -> float:
    """Fortran の sign(a, b): |a| に b の符号 (b=0 は正)"""
    return abs(a) if b >= 0.0 else -abs(a)


def tsvd(a_in, b_in, epslon: float, svd: bool = False, itmax: int = 256, warn=None):
    """a_in: (n, m) 配列、b_in: (n,)。戻り: x (m,)。svd=True の分岐 (V, W^-1 の返却) は未対応"""
    if svd:
        raise NotImplementedError('tsvd(svd=True) は未移植')
    a0 = np.asarray(a_in, dtype=float)
    n, m = a0.shape
    if n == 0 or m == 0:
        return np.zeros(m)
    # 1 始まりに
    a = np.zeros((n + 2, m + 2))
    a[1:n + 1, 1:m + 1] = a0
    b = np.zeros(n + 2)
    b[1:n + 1] = np.asarray(b_in, dtype=float)
    mn = min(n, m)
    v = np.zeros(2 * max(n, m) + 3)
    aam = np.zeros(m + 2)
    bbm = np.zeros(n + 2)
    lsep = np.zeros(n + 3, dtype=int)
    x = np.ones(m + 2)
    v[1:n + 1] = 1.0

    # ---- 2 重対角化
    for i in range(1, mn + 1):
        i1 = i + 1
        if i < n:
            for j in range(i1, n + 1):
                if abs(a[i, i]) > abs(a[j, i]):
                    p = a[j, i] / a[i, i]
                    h1 = v[i] + v[j] * p ** 2
                    q = v[j] * p / h1
                    v[j] = v[i] * v[j] / h1
                    v[i] = h1
                    a[j, i1:m + 1] = a[j, i1:m + 1] - p * a[i, i1:m + 1]
                    a[i, i1:m + 1] = a[i, i1:m + 1] + q * a[j, i1:m + 1]
                    a[j, i] = 0.0
                    b[j] = b[j] - p * b[i]
                    b[i] = b[i] + q * b[j]
                elif a[j, i] != 0.0:
                    p = a[i, i] / a[j, i]
                    h1 = v[j] + v[i] * p ** 2
                    q = v[i] * p / h1
                    v[j] = v[i] * v[j] / h1
                    v[i] = h1
                    aam[i1:m + 1] = a[j, i1:m + 1]
                    a[j, i1:m + 1] = p * aam[i1:m + 1] - a[i, i1:m + 1]
                    a[i, i1:m + 1] = aam[i1:m + 1] - q * a[j, i1:m + 1]
                    a[i, i] = a[j, i]
                    a[j, i] = 0.0
                    bb = b[j]
                    b[j] = p * bb - b[i]
                    b[i] = bb - q * b[j]
        if i1 < m:
            for j in range(i1 + 1, m + 1):
                r1 = x[i1] * a[i, i1]
                r2 = x[j] * a[i, j]
                r = _sign(math.hypot(r1, r2), r1)
                if r != 0.0:
                    c = r1 / r
                    s = r2 / r
                    if abs(c) > abs(s):
                        a[i, j] = 0.0
                        h1 = x[i1] / c
                        h2 = x[j] * c
                        p = s * x[i1] / h2
                        q = s * x[j] / h1
                        x[i1] = h1
                        x[j] = h2
                        a[i1:n + 1, j] = a[i1:n + 1, j] - p * a[i1:n + 1, i1]
                        a[i1:n + 1, i1] = a[i1:n + 1, i1] + q * a[i1:n + 1, j]
                    else:
                        a[i, i1] = a[i, j]
                        a[i, j] = 0.0
                        h1 = x[j] / s
                        h2 = x[i1] * s
                        p = c * x[j] / h2
                        q = c * x[i1] / h1
                        x[i1] = h1
                        x[j] = h2
                        bbm[i1:n + 1] = a[i1:n + 1, j]
                        a[i1:n + 1, j] = p * bbm[i1:n + 1] - a[i1:n + 1, i1]
                        a[i1:n + 1, i1] = bbm[i1:n + 1] - q * a[i1:n + 1, j]
                else:
                    c = 1.0
                    s = 0.0
                a[i, j] = s
                if j < n + 2:
                    a[j - 1, i] = c
                else:
                    if c <= abs(s) and c != 0.0:
                        a[i, j] = _sign(1.0 / c, s)

    anorm = 0.0
    xmin = 1e300
    for i in range(1, mn + 1):
        p = math.sqrt(v[i])
        v[i + mn] = x[i]
        a[i, i] = a[i, i] * p
        x[i] = a[i, i] * x[i]
        b[i] = b[i] * p
        if i < m:
            a[i, i + 1] = a[i, i + 1] * p
            v[i] = a[i, i + 1] * x[i + 1]
        else:
            v[i] = 0.0
        an = abs(x[i]) + abs(v[i])
        anorm = max(anorm, an)
        if an != 0.0:
            xmin = min(xmin, abs(x[i]))

    enorm = anorm * epslon
    # ---- quick path: 特異値を落とさずに 2 重対角系を直接解く
    if xmin > enorm:
        if mn + 2 <= m:
            x[mn + 2:m + 1] = 0.0
        r = math.hypot(x[mn], v[mn])
        if r != 0.0:
            w = 1.0 / r
            d = x[mn] * w
            b[mn] = b[mn] * w
            v[mn] = v[mn] * w
            if mn < m:
                x[mn + 1] = v[mn] * b[mn]
            x[mn] = d * b[mn]
            ma = m
        else:
            ma = mn
            d = 0.0
            v[mn] = 0.0
            b[mn] = 0.0
            x[mn] = 0.0
            if mn < m:
                x[mn + 1] = 0.0
        for i in range(mn - 1, 0, -1):
            if x[i] != 0.0:
                p = v[i] * d
                ra = abs(x[i]) + abs(v[i])
                r = (x[i] / ra) ** 2 + ((v[i] - p) / ra) * ((v[i] + p) / ra)
                w = 1.0 / math.sqrt(r) / ra
                b[i] = (b[i] - p * b[i + 1]) * w
                u = -p * w
                hi_v = min(mn, ma - 1)
                if hi_v >= i + 1:
                    v[i + 1:hi_v + 1] = u * v[i + 1:hi_v + 1]
                hi_x = min(mn + 1, ma)
                if hi_x >= i + 2:
                    x[i + 2:hi_x + 1] = x[i + 2:hi_x + 1] + b[i] * v[i + 1:hi_v + 1]
                v[i] = (v[i] - p * d) * w
                x[i + 1] = x[i + 1] + b[i] * v[i]
                d = x[i] * w
                x[i] = d * b[i]
            else:
                ma = i
                d = 0.0
                v[i] = 0.0
                b[i] = 0.0
                x[i] = 0.0
        for i in range(min(mn, m - 2), 0, -1):
            for j in range(m, i + 1, -1):
                s = a[i, j]
                if j < n + 2:
                    c = a[j - 1, i]
                else:
                    if abs(s) > 1.0:
                        c = abs(1.0 / s)
                        s = _sign(1.0 - c ** 2 / (1.0 + math.sqrt((1.0 - c) * (1.0 + c))), s)
                    else:
                        c = 1.0 - s ** 2 / (1.0 + math.sqrt((1.0 - s) * (1.0 + s)))
                aa = x[i + 1]
                x[i + 1] = c * aa - s * x[j]
                x[j] = s * aa + c * x[j]
        return x[1:m + 1].copy()

    # ---- 通常の経路: V を作り、暗黙 QR 法で特異値を求める
    v[2 * mn + 1:mn + m + 1] = x[mn + 1:m + 1]
    for i in range(min(mn, m - 2), 0, -1):
        i1 = i + 1
        i1mn = i1 + mn
        for j in range(m, i + 1, -1):
            s = a[i, j]
            if j < n + 2:
                c = a[j - 1, i]
                a[j - 1, i] = 0.0
            else:
                if abs(s) > 1.0:
                    c = abs(1.0 / s)
                    s = _sign(1.0 - c ** 2 / (1.0 + math.sqrt((1.0 - c) * (1.0 + c))), s)
                else:
                    c = 1.0 - s ** 2 / (1.0 + math.sqrt((1.0 - s) * (1.0 + s)))
            if abs(c) > abs(s):
                h1 = v[i1mn] * c
                h2 = v[j + mn] / c
                p = s * v[j + mn] / h1
                q = s * v[i1mn] / h2
                v[i1mn] = h1
                v[j + mn] = h2
                a[i, j] = q * a[i, i1]
                a[i1:mn + 1, i1] = a[i1:mn + 1, i1] - p * a[i1:mn + 1, j]
                a[i1:mn + 1, j] = a[i1:mn + 1, j] + q * a[i1:mn + 1, i1]
            else:
                h1 = v[j + mn] / s
                h2 = v[i1mn] * s
                p = c * v[j + mn] / h2
                q = c * v[i1mn] / h1
                v[i1mn] = h1
                v[j + mn] = h2
                a[i, j] = a[i, i1]
                a[i, i1] = a[i, i1] * q
                for k in range(i1, mn + 1):
                    aa = a[k, j]
                    a[k, j] = p * aa + a[k, i1]
                    a[k, i1] = q * a[k, j] - aa

    lsep[0] = 0
    lsep[1] = 1
    isep = 1
    ibegin = 1
    iend = mn
    v[mn + 1:2 * mn + 1] = 1.0
    v[0] = 0.0
    nfail = 4
    done = False
    while not done:                                  # label 1002
        if v[iend] != 0.0:
            f = v[iend]
            v[iend] = 0.0
            for i in range(iend, ibegin - 1, -1):
                if abs(f) + abs(x[i]) != abs(x[i]):
                    p = math.hypot(x[i], f)
                    vv = v[i - 1] / p
                    v[i - 1] = vv * x[i]
                    x[i] = p
                    f = -vv * f
                else:
                    break
        restart = False
        while True:                                  # do3001
            if ibegin >= iend:
                isep -= 1
                if isep <= 0:
                    done = True
                    break
                iend = ibegin - 1
                ibegin = lsep[isep]
                continue
            cycled = False
            for it in range(1, itmax + 1):
                if x[iend] == 0.0:
                    iend -= 1
                    restart = True
                    break
                hit = False
                for i in range(iend - 1, ibegin - 1, -1):
                    an = abs(x[i]) + abs(x[i + 1])
                    if abs(v[i]) <= an * 1e-16:
                        v[i] = 0.0
                        if i == iend - 1:
                            iend -= 1
                        else:
                            isep += 1
                            ibegin = i + 1
                            lsep[isep] = ibegin
                        hit = True
                        break
                if hit:
                    cycled = True
                    break
                # svd=False: 対角がすべて enorm 以上なら、残りを後退代入で片付けて終わる (do1611)
                if all(abs(x[i]) >= enorm for i in range(1, mn + 1)):
                    for i in range(iend, ibegin, -1):
                        p = v[i + mn] * v[i - 1] / v[i - 1 + mn] / x[i]
                        a[i - 1, 1:m + 1] = a[i - 1, 1:m + 1] - p * a[i, 1:m + 1]
                        b[i - 1] = b[i - 1] - p * b[i]
                        v[i - 1] = 0.0
                    iend = ibegin - 1
                    cycled = True
                    break
                # ---- QR 1 ステップ
                w = x[ibegin]
                z = x[iend]
                y = x[iend - 1]
                g = v[iend - 2] if ibegin < iend - 1 else 0.0
                h = v[iend - 1]
                f = ((y - z) * (y + z) + (g - h) * (g + h)) * 0.5
                if w == 0.0:
                    f = 0.0
                else:
                    g = h * y
                    y = f + _sign(math.hypot(f, g), f)
                    if y != 0.0:
                        f = ((w - z) * (w + z) - h ** 2 + g ** 2 / y) / w
                    else:
                        f = ((w - z) * (w + z) - h ** 2) / w
                g = v[ibegin]
                h = g
                finished = True
                for i in range(ibegin, iend):
                    i1 = i + 1
                    z = math.hypot(f, g)
                    v[i - 1] = z
                    if z != 0.0:
                        c = f / z
                        s = g / z
                    else:
                        x[i] = w
                        v[i] = h
                        finished = False
                        break
                    f = w * c + h * s
                    g = -w * s + h * c
                    h = x[i1] * s
                    y = x[i1] * c
                    z = math.hypot(f, h)
                    x[i] = z
                    if z != 0.0:
                        c = f / z
                        s = h / z
                    else:
                        v[i] = g
                        x[i1] = y
                        finished = False
                        break
                    f = c * g + s * y
                    w = -s * g + c * y
                    g = v[i1] * s
                    h = v[i1] * c
                    if abs(c) > abs(s):
                        h1 = v[i + mn] / c
                        h2 = v[i1 + mn] * c
                        r = s * v[i + mn] / h2
                        t = s * v[i1 + mn] / h1
                        v[i + mn] = h1
                        v[i1 + mn] = h2
                        a[i1, 1:m + 1] = a[i1, 1:m + 1] - r * a[i, 1:m + 1]
                        a[i, 1:m + 1] = a[i, 1:m + 1] + t * a[i1, 1:m + 1]
                        b[i1] = b[i1] - r * b[i]
                        b[i] = b[i] + t * b[i1]
                    else:
                        h1 = v[i1 + mn] / s
                        h2 = v[i + mn] * s
                        r = c * v[i1 + mn] / h2
                        t = c * v[i + mn] / h1
                        v[i + mn] = h1
                        v[i1 + mn] = h2
                        aam[1:m + 1] = a[i1, 1:m + 1]
                        a[i1, 1:m + 1] = r * aam[1:m + 1] - a[i, 1:m + 1]
                        a[i, 1:m + 1] = aam[1:m + 1] - t * a[i1, 1:m + 1]
                        bb = b[i1]
                        b[i1] = r * bb - b[i]
                        b[i] = bb - t * b[i1]
                if finished:
                    v[iend - 1] = f
                    x[iend] = w
                v[ibegin - 1] = 0.0
            else:
                # itmax 回で収束しなかった
                if nfail >= 0 and warn is not None:
                    warn(f' TSVD convergence failed: {iend:10d}')
                nfail -= 1
                iend -= 1
                v[iend] = 0.0
                continue
            if restart:
                break
            if cycled:
                continue
        if done:
            break

    anorm = float(np.max(np.abs(x[1:mn + 1]))) * epslon
    for i in range(1, mn + 1):
        s = x[i]
        if abs(s) > anorm:
            w = (v[i + mn] / s) ** 2
        else:
            w = (s * v[i + mn] / anorm ** 2) ** 2
        b[i] = b[i] * w
    return np.array([float(np.sum(a[1:mn + 1, i] * b[1:mn + 1])) for i in range(1, m + 1)])
