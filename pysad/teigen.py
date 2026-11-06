# -*- coding: utf-8 -*-
"""teigen.f (SAD の実行列の固有値・固有ベクトル、K. Oide) の直訳。

teigen(a) = tbal (釣り合いの重み vx のみ求める) → thess (Hessenberg 化, 高速 Givens) → tqr (QR 法)。
固有値の並び・実固有値の対の並び・固有ベクトルの符号と大きさは numpy.linalg.eig とは異なる。
不安定な光学 (実固有値の対 λ, 1/λ) で tnorm の後の並びを SAD と合わせるために必要 (docs/history/work_emit_unstable.md)。

配列は Fortran と同じ 1 始まりで書く (サイズ n+1 の配列の 0 番は使わない)。
"""
import math
import numpy as np


def _tbal(a_in, n):
    """tbal(w=a_in, a=...): 釣り合いの重み vx を求める (行列そのものは変えない)。
    戻り値: (A, W, vx)  A = 入力行列のコピー (Hessenberg 化する行列), W = 単位行列 (固有ベクトルの累積)"""
    A = np.zeros((n + 1, n + 1))
    A[1:, 1:] = a_in
    W = np.zeros((n + 1, n + 1))
    for i in range(1, n + 1):
        W[i, i] = 1.0
    vx = np.ones(n + 1)
    radix = 2.0
    sqrdx = radix ** 2
    rth = 0.95
    last = True
    while last:
        last = False
        for i in range(1, n + 1):
            c = 0.0
            r = 0.0
            for j in range(1, n + 1):
                if j != i:
                    u = vx[j] / vx[i]
                    c = c + abs(A[j, i] * u)
                    r = r + abs(A[i, j] / u)
            if c != 0.0 and r != 0.0:
                g = r / radix
                f = 1.0
                s = c + r
                while c < g:
                    f = f * radix
                    c = c * sqrdx
                g = r * radix
                while c > g:
                    f = f / radix
                    c = c / sqrdx
                if (c + r) / f < rth * s:
                    last = True
                    vx[i] = vx[i] / f
    return A, W, vx


def _thess(a, w, vx, n):
    """thess: 高速 Givens 変換で Hessenberg 形に。a, w, vx は 1 始まり配列で in-place"""
    vx[2:n + 1] = vx[2:n + 1] ** 2
    for i in range(1, n - 1):
        i1 = i + 1
        for j in range(i1 + 1, n + 1):
            if a[j, i] != 0.0:
                if abs(a[i1, i]) > abs(a[j, i]):
                    p = a[j, i] / a[i1, i]
                    h1 = vx[i1] + vx[j] * p ** 2
                    q = vx[j] * p / h1
                    vx[j] = vx[i1] * vx[j] / h1
                    vx[i1] = h1
                    a[j, i] = 0.0
                    a[j, i1:n + 1] = a[j, i1:n + 1] - p * a[i1, i1:n + 1]
                    a[i1, i1:n + 1] = a[i1, i1:n + 1] + q * a[j, i1:n + 1]
                    a[:, i1] = a[:, i1] + p * a[:, j]
                    a[:, j] = a[:, j] - q * a[:, i1]
                    w[1:n + 1, i1] = w[1:n + 1, i1] + p * w[1:n + 1, j]
                    w[1:n + 1, j] = w[1:n + 1, j] - q * w[1:n + 1, i1]
                else:
                    p = a[i1, i] / a[j, i]
                    h1 = vx[j] + vx[i1] * p ** 2
                    q = vx[i1] * p / h1
                    vx[j] = vx[i1] * vx[j] / h1
                    vx[i1] = h1
                    aa = a[j, i1:n + 1].copy()
                    a[j, i1:n + 1] = p * aa - a[i1, i1:n + 1]
                    a[i1, i1:n + 1] = aa - q * a[j, i1:n + 1]
                    a[i1, i] = a[j, i]
                    a[j, i] = 0.0
                    aa = a[:, i1].copy()
                    a[:, i1] = p * aa + a[:, j]
                    a[:, j] = q * a[:, i1] - aa
                    aa = w[1:n + 1, i1].copy()
                    w[1:n + 1, i1] = p * aa + w[1:n + 1, j]
                    w[1:n + 1, j] = q * w[1:n + 1, i1] - aa
    vx[2:n + 1] = np.sqrt(vx[2:n + 1])


def _tqr(a, w, vx, n):
    """tqr: Hessenberg 行列の QR 法 (二重シフト・高速 Givens) と固有ベクトルの逆代入。戻り値 (eig1, eig2)"""
    itmax = 100
    vmax = 1.0e30
    vmin = 1.0 / vmax
    decth = 1.0e-16
    threj = 1.0e-8
    alpha = 1.0e-6
    amth = 1.0e-30
    eig1 = np.zeros(n + 2)
    eig2 = np.zeros(n + 2)
    jordan = False
    ibp = 1
    ibtab = [0] * (n + 3)
    ibtab[1] = 1
    ib = 1
    ie = n
    for i in range(1, ie - 1):
        a[i + 2, i] = 0.0
        a[min(i + 3, ie), i] = 0.0
    anorm = 0.0
    for i in range(1, n + 1):
        for j in range(max(1, i - 1), n + 1):
            anorm = anorm + abs(vx[i] / vx[j] * a[i, j])
    anorm = anorm * 2 / n ** 2
    is2 = 0
    w1 = w2 = w3 = 0.0
    am = am1 = 0.0
    while True:   # do1
        if ie <= ib + 1:
            if ib == 1:
                break
            ie = ib - 1
            ibp = ibp - 1
            ib = ibtab[ibp]
            continue
        ie1 = ie - 1
        it = 0
        itm = min(itmax, 10 * max(ie - ib + 1, 2))
        while True:   # do2 (cycle do1 は break)
            restart = False
            for i in range(ie - 1, ib - 1, -1):
                i1 = i + 1
                s = max(anorm, abs(a[i1, i1]) + abs(a[i, i]))
                if abs(vx[i1] / vx[i] * a[i1, i]) + s == s:
                    a[i1, i] = 0.0
                    if i1 == ie:
                        ie = ie - 1
                        restart = True
                        break
                    ibp = ibp + 1
                    ibtab[ibp] = ib
                    ib = i1
                    restart = True
                    break
            if restart:
                break
            for i in range(ib, ie + 1):
                if abs(vx[i]) < vmin or abs(vx[i]) > vmax:
                    j1 = max(i - 1, ib)
                    a[i, j1:n + 1] = a[i, j1:n + 1] * vx[i]
                    j1 = min(ie, i + 1)
                    a[1:j1 + 1, i] = a[1:j1 + 1, i] / vx[i]
                    w[1:n + 1, i] = w[1:n + 1, i] / vx[i]
                    vx[i] = 1.0
            for is_ in range(ie - 2, ib - 1, -1):
                is1 = is_ + 1
                is2 = is1 + 1
                a1 = a[ie, ie] - a[is_, is_]
                a2 = a[ie1, ie1] - a[is_, is_]
                w1 = ((a1 * a2 - a[ie1, ie] * a[ie, ie1]) / a[is1, is_] + a[is_, is1]) * vx[is_] / vx[is1]
                if w1 == 0.0:
                    w1 = -alpha * a[ie1, ie] * a[ie, ie1] / a[is1, is_] * vx[is_] / vx[is1]
                w2 = a[is1, is1] - a[is_, is_] - a1 - a2
                w3 = a[is2, is1] * vx[is2] / vx[is1]
                if is_ == ib:
                    break
                u = abs(a[is_, is_ - 1] * vx[is_] / vx[is_ - 1]) * (abs(w2) + abs(w3))
                v = abs(w1) * (abs(a[is_ - 1, is_ - 1]) + abs(a[is_, is_]) + abs(a[is1, is1]))
                if u + v == v:
                    a[is1, is_ - 1] = 0.0
                    a[is2, is_ - 1] = 0.0
                    break
            a[is2, is_] = 0.0
            a[min(is2 + 1, ie), is_] = 0.0
            am = math.hypot(w1, w2)
            am1 = math.hypot(am, w3)
            if is_ > ib:
                if am1 != 0.0:
                    a[is_, is_ - 1] = a[is_, is_ - 1] * w1 / am1 * vx[is_]
            for i in range(is_ - 1, ie - 1):
                i1 = i + 1
                i2 = i1 + 1
                i3 = i2 + 1
                if i1 > is_:
                    w1 = vx[i1] * a[i1, i]
                    w2 = vx[i2] * a[i2, i]
                    a[i2, i] = 0.0
                    am = math.hypot(w1, w2)
                    if i3 <= ie:
                        w3 = vx[i3] * a[i3, i]
                        a[i3, i] = 0.0
                        am1 = math.hypot(am, w3)
                    else:
                        w3 = 0.0
                        am1 = 0.0
                if am != 0.0:
                    c = w1 / am
                    s = w2 / am
                    if abs(c) > abs(s):
                        v1 = vx[i1] / c
                        v2 = c * vx[i2]
                        p = s * vx[i1] / v2
                        q = s * vx[i2] / v1
                        vx[i1] = v1
                        if abs(v2) > vmin and abs(v2) < vmax:
                            vx[i2] = v2
                            a[i2, i1:n + 1] = a[i2, i1:n + 1] - p * a[i1, i1:n + 1]
                            a[i1, i1:n + 1] = a[i1, i1:n + 1] + q * a[i2, i1:n + 1]
                            j1 = min(ie, i3)
                            a[1:j1 + 1, i1] = a[1:j1 + 1, i1] + p * a[1:j1 + 1, i2]
                            a[1:j1 + 1, i2] = a[1:j1 + 1, i2] - q * a[1:j1 + 1, i1]
                            w[1:n + 1, i1] = w[1:n + 1, i1] + p * w[1:n + 1, i2]
                            w[1:n + 1, i2] = w[1:n + 1, i2] - q * w[1:n + 1, i1]
                        else:
                            vx[i2] = 1.0
                            aa = a[i2, i1:n + 1] - p * a[i1, i1:n + 1]
                            a[i1, i1:n + 1] = a[i1, i1:n + 1] + q * aa
                            a[i2, i1:n + 1] = v2 * aa
                            j1 = min(ie, i3)
                            a[1:j1 + 1, i1] = a[1:j1 + 1, i1] + p * a[1:j1 + 1, i2]
                            a[1:j1 + 1, i2] = (a[1:j1 + 1, i2] - q * a[1:j1 + 1, i1]) / v2
                            w[1:n + 1, i1] = w[1:n + 1, i1] + p * w[1:n + 1, i2]
                            w[1:n + 1, i2] = (w[1:n + 1, i2] - q * w[1:n + 1, i1]) / v2
                    else:
                        v1 = vx[i2] / s
                        v2 = vx[i1] * s
                        if i1 > is_:
                            a[i1, i] = w2 / vx[i2]
                        p = c * vx[i2] / v2
                        q = c * vx[i1] / v1
                        vx[i1] = v1
                        if abs(v2) > vmin and abs(v2) < vmax:
                            vx[i2] = v2
                            aa = a[i2, i1:n + 1].copy()
                            a[i2, i1:n + 1] = p * aa - a[i1, i1:n + 1]
                            a[i1, i1:n + 1] = aa - q * a[i2, i1:n + 1]
                            j1 = min(ie, i3)
                            aa = a[1:j1 + 1, i1].copy()
                            a[1:j1 + 1, i1] = p * aa + a[1:j1 + 1, i2]
                            a[1:j1 + 1, i2] = q * a[1:j1 + 1, i1] - aa
                            aa = w[1:n + 1, i1].copy()
                            w[1:n + 1, i1] = p * aa + w[1:n + 1, i2]
                            w[1:n + 1, i2] = q * w[1:n + 1, i1] - aa
                        else:
                            vx[i2] = 1.0
                            aa = a[i2, i1:n + 1].copy()
                            aa1 = p * aa - a[i1, i1:n + 1]
                            a[i1, i1:n + 1] = aa - q * aa1
                            a[i2, i1:n + 1] = aa1 * v2
                            j1 = min(ie, i3)
                            aa = a[1:j1 + 1, i1].copy()
                            a[1:j1 + 1, i1] = p * aa + a[1:j1 + 1, i2]
                            a[1:j1 + 1, i2] = (q * a[1:j1 + 1, i1] - aa) / v2
                            aa = w[1:n + 1, i1].copy()
                            w[1:n + 1, i1] = p * aa + w[1:n + 1, i2]
                            w[1:n + 1, i2] = (q * w[1:n + 1, i1] - aa) / v2
                if am1 != 0.0:
                    c = am / am1
                    s = w3 / am1
                    if abs(c) > abs(s):
                        v1 = vx[i1] / c
                        v2 = c * vx[i3]
                        p = s * vx[i1] / v2
                        q = s * vx[i3] / v1
                        vx[i1] = v1
                        vx[i3] = v2
                        a[i3, i1:n + 1] = a[i3, i1:n + 1] - p * a[i1, i1:n + 1]
                        a[i1, i1:n + 1] = a[i1, i1:n + 1] + q * a[i3, i1:n + 1]
                        j1 = min(ie, i3 + 1)
                        a[1:j1 + 1, i1] = a[1:j1 + 1, i1] + p * a[1:j1 + 1, i3]
                        a[1:j1 + 1, i3] = a[1:j1 + 1, i3] - q * a[1:j1 + 1, i1]
                        w[1:n + 1, i1] = w[1:n + 1, i1] + p * w[1:n + 1, i3]
                        w[1:n + 1, i3] = w[1:n + 1, i3] - q * w[1:n + 1, i1]
                    else:
                        v1 = vx[i3] / s
                        v2 = vx[i1] * s
                        if i1 > is_:
                            a[i1, i] = w3 / vx[i3]
                        p = c * vx[i3] / v2
                        q = c * vx[i1] / v1
                        vx[i1] = v1
                        vx[i3] = v2
                        aa = a[i3, i1:n + 1].copy()
                        a[i3, i1:n + 1] = p * aa - a[i1, i1:n + 1]
                        a[i1, i1:n + 1] = aa - q * a[i3, i1:n + 1]
                        j1 = min(ie, i3 + 1)
                        aa = a[1:j1 + 1, i1].copy()
                        a[1:j1 + 1, i1] = p * aa + a[1:j1 + 1, i3]
                        a[1:j1 + 1, i3] = q * a[1:j1 + 1, i1] - aa
                        aa = w[1:n + 1, i1].copy()
                        w[1:n + 1, i1] = p * aa + w[1:n + 1, i3]
                        w[1:n + 1, i3] = q * w[1:n + 1, i1] - aa
            if is_ > ib:
                a[is_, is_ - 1] = a[is_, is_ - 1] / vx[is_]
            it = it + 1
            if it > itm:
                # 収束しなかった: SAD は "TEIGEN convergence failed" を印字して先へ進む
                a[ie, ie - 1] = 0.0
                if math.isnan(a[ie, ie]):
                    a[ie, ie] = 0.0
                    a[ie - 1, ie - 1] = 0.0
                if math.isnan(vx[ie]):
                    vx[ie] = 1.0
                if math.isnan(vx[ie - 1]):
                    vx[ie - 1] = 1.0
                it = 0
                ie = ie - 1
                break
    for i in range(1, n + 1):
        j1 = max(i - 1, 1)
        a[i, j1:n + 1] = a[i, j1:n + 1] * vx[i]
        j1 = min(n, i + 1)
        a[1:j1 + 1, i] = a[1:j1 + 1, i] / vx[i]
        w[1:n + 1, i] = w[1:n + 1, i] / vx[i]
    paired = False
    for i in range(1, n):
        if paired:
            paired = False
            continue
        i1 = i + 1
        if a[i1, i] != 0.0:
            if abs(a[i1, i]) > abs(a[i, i1]):
                aa = a[i, i:n + 1].copy()
                a[i, i:n + 1] = a[i1, i:n + 1]
                a[i1, i:n + 1] = aa
                aa = a[1:i1 + 1, i].copy()
                a[1:i1 + 1, i] = a[1:i1 + 1, i1]
                a[1:i1 + 1, i1] = aa
                aa = w[1:n + 1, i].copy()
                w[1:n + 1, i] = w[1:n + 1, i1]
                w[1:n + 1, i1] = aa
            if a[i1, i] == 0.0:
                eig1[i] = a[i, i]
                eig2[i] = 0.0
                continue
            s = a[i, i] - a[i1, i1]
            dec = s ** 2 + 4.0 * a[i1, i] * a[i, i1]
            s1 = a[i, i] ** 2 + a[i1, i1] ** 2
            if dec >= -decth * s1:
                sqrd = math.sqrt(abs(dec))
                p = 2.0 * a[i1, i] / (s + math.copysign(sqrd, s))
                a[i1, i:n + 1] = a[i1, i:n + 1] - p * a[i, i:n + 1]
                a[1:i1 + 1, i] = a[1:i1 + 1, i] + p * a[1:i1 + 1, i1]
                w[1:n + 1, i] = w[1:n + 1, i] + p * w[1:n + 1, i1]
                eig1[i] = a[i, i]
                eig1[i1] = a[i1, i1]
                eig2[i] = 0.0
                eig2[i1] = 0.0
                a[i1, i] = 0.0
            else:
                p = -s / a[i, i1] * .5
                r = math.copysign(math.sqrt(abs((a[i1, i] - p * (a[i, i] - a[i1, i1] + p * a[i, i1])) / a[i, i1])),
                                  a[i, i1])
                a[i1, i:n + 1] = a[i1, i:n + 1] - p * a[i, i:n + 1]
                a[i, i:n + 1] = a[i, i:n + 1] * r
                a[1:i1 + 1, i] = (a[1:i1 + 1, i] + p * a[1:i1 + 1, i1]) / r
                w[1:n + 1, i] = (w[1:n + 1, i] + p * w[1:n + 1, i1]) / r
                eig1[i] = a[i, i]
                eig1[i1] = a[i1, i1]
                eig2[i] = a[i, i1]
                eig2[i1] = a[i1, i]
            paired = True
        else:
            eig1[i] = a[i, i]
            eig2[i] = 0.0
    if not paired:
        eig1[n] = a[n, n]
        eig2[n] = 0.0
    pk = pl = pm = pn = 0.0
    for i in range(2, n + 1):
        if a[i, i - 1] != 0.0:
            continue
        i1 = i + 1
        if eig2[i] == 0.0:
            for j in range(i - 1, 0, -1):
                if a[j + 1, j] != 0.0:
                    continue
                if eig2[j] == 0.0:
                    da = a[j, j] - a[i, i]
                    sa = abs(a[j, j]) + abs(a[i, i])
                    if abs(da) > threj * sa:
                        p = a[j, i] / da
                        a[j, i:n + 1] = a[j, i:n + 1] + p * a[i, i:n + 1]
                        a[1:j + 1, i] = a[1:j + 1, i] - p * a[1:j + 1, j]
                        w[1:n + 1, i] = w[1:n + 1, i] - p * w[1:n + 1, j]
                    else:
                        jordan = True
                else:
                    j1 = j - 1
                    u = a[i, i] - a[j1, j1]
                    v = a[j1, j]
                    r = u ** 2 + v ** 2
                    pk = (a[j1, i] * u + a[j, i] * v) / r
                    pl = (a[j, i] * u - a[j1, i] * v) / r
                    a[j1, i:n + 1] = a[j1, i:n + 1] - pk * a[i, i:n + 1]
                    a[j, i:n + 1] = a[j, i:n + 1] - pl * a[i, i:n + 1]
                    a[1:j + 1, i] = a[1:j + 1, i] + pk * a[1:j + 1, j1] + pl * a[1:j + 1, j]
                    w[1:n + 1, i] = w[1:n + 1, i] + pk * w[1:n + 1, j1] + pl * w[1:n + 1, j]
        else:
            for j in range(i - 1, 0, -1):
                if a[j + 1, j] != 0.0:
                    continue
                if eig2[j] == 0.0:
                    x = a[i, i] - a[j, j]
                    y = a[i, i1]
                    r = x ** 2 + y ** 2
                    pk = (a[j, i] * x + a[j, i1] * y) / r
                    pl = (a[j, i1] * x - a[j, i] * y) / r
                    a[j, i:n + 1] = a[j, i:n + 1] - pk * a[i, i:n + 1] - pl * a[i1, i:n + 1]
                    a[1:j + 1, i] = a[1:j + 1, i] + pk * a[1:j + 1, j]
                    a[1:j + 1, i1] = a[1:j + 1, i1] + pl * a[1:j + 1, j]
                    w[1:n + 1, i] = w[1:n + 1, i] + pk * w[1:n + 1, j]
                    w[1:n + 1, i1] = w[1:n + 1, i1] + pl * w[1:n + 1, j]
                else:
                    j1 = j - 1
                    cu = complex(a[j1, j1] - a[i, i], -a[i, i1])
                    v = a[j1, j]
                    cr = cu ** 2 + v ** 2
                    if cr != 0j:
                        ca = complex(a[j1, i], a[j1, i1])
                        cc = complex(a[j, i], a[j, i1])
                        ck = (cu * ca - v * cc) / cr
                        cm = (cu * cc + v * ca) / cr
                        pk, pl = ck.real, ck.imag
                        pm, pn = cm.real, cm.imag
                        a[j1, i:n + 1] = a[j1, i:n + 1] + pk * a[i, i:n + 1] + pl * a[i1, i:n + 1]
                        a[j, i:n + 1] = a[j, i:n + 1] + pm * a[i, i:n + 1] + pn * a[i1, i:n + 1]
                        a[1:j + 1, i] = a[1:j + 1, i] - pk * a[1:j + 1, j1] - pm * a[1:j + 1, j]
                        a[1:j + 1, i1] = a[1:j + 1, i1] - pl * a[1:j + 1, j1] - pn * a[1:j + 1, j]
                        w[1:n + 1, i] = w[1:n + 1, i] - pk * w[1:n + 1, j1] - pm * w[1:n + 1, j]
                        w[1:n + 1, i1] = w[1:n + 1, i1] - pl * w[1:n + 1, j1] - pn * w[1:n + 1, j]
    if jordan:
        for i in range(n, 1, -1):
            if eig2[i] == 0.0:
                jm = 0
                for j in range(i - 1, 0, -1):
                    if eig1[j] == eig1[i] and eig2[j] == 0.0:
                        r = a[j, i]
                        if abs(r) > amth:
                            if jm == 0:
                                r1 = math.sqrt(abs(r))
                                r2 = r / r1
                                a[j, j:n + 1] = a[j, j:n + 1] / r1
                                a[1:j + 1, j] = a[1:j + 1, j] * r1
                                w[1:n + 1, j] = w[1:n + 1, j] * r1
                                a[1:i + 1, i] = a[1:i + 1, i] / r2
                                a[i, i:n + 1] = a[i, i:n + 1] * r2
                                w[1:n + 1, i] = w[1:n + 1, i] / r2
                                jm = j
                            else:
                                a[j, jm:n + 1] = a[j, jm:n + 1] - r * a[jm, jm:n + 1]
                                a[1:j + 1, jm] = a[1:j + 1, jm] + r * a[1:j + 1, j]
                                w[1:n + 1, jm] = w[1:n + 1, jm] + r * w[1:n + 1, j]
                        else:
                            a[j, i] = 0.0
    if n % 2 == 0:
        ii = 0
        for i in range(1, n + 1):
            if eig2[i] > 0.0 and ii != 0:
                aa = w[1:n + 1, i - 1].copy()
                w[1:n + 1, i - 1] = w[1:n + 1, i]
                w[1:n + 1, i] = w[1:n + 1, i + 1]
                w[1:n + 1, i + 1] = aa
                ee = eig1[i - 1]
                eig1[i - 1] = eig1[i]
                eig1[i] = eig1[i + 1]
                eig1[i + 1] = ee
                ee = eig2[i - 1]
                eig2[i - 1] = eig2[i]
                eig2[i] = eig2[i + 1]
                eig2[i + 1] = ee
            ii = 1 - ii
    return eig1, eig2


def teigen(a_in):
    """teigen(a,w,eig,n,ndim): 実行列 a_in (n×n) の固有値 (複素配列) と固有ベクトル (列, n×n) を返す。
    複素固有値は SAD と同じ並び (共役の対で虚部が正 → 負)"""
    a_in = np.asarray(a_in, dtype=float)
    n = a_in.shape[0]
    A, W, vx = _tbal(a_in, n)
    # teigen: tbal(a→w, ...) thess(w, a) tqr(w, a): 行列は A、固有ベクトルの累積は W
    _thess(A, W, vx, n)
    e1, e2 = _tqr(A, W, vx, n)
    ceig = np.array([complex(e1[i], e2[i]) for i in range(1, n + 1)])
    return ceig, W[1:, 1:].copy()
