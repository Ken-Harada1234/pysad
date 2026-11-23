"""MAIN の TRACK の FOURIER 解析 (K64 の wjfit.f・wfres.f・wfbin.f・wiord.f の直訳)。

`ON SMEAR; ON FOURIER;` で TRACK した粒子ごとの複素座標 z = (x+i px, y+i py, z+i g) (規格化座標、周ごと) を、
作用変数 J の多項式 + 共鳴項 (mx, my, mz の整数の組ごとの exp(i(mx φx + my φy + mz φz))) でフィットする。
 - `wiord`: 組 (mx, my, mz) の並び (|mx|+|my|+|mz| ≤ MAXORDER)。
 - `wfres`: 各粒子のチューン (wtune) を求め、ADDTERMS 本ぶん (mx, my, mz) の組を、チューンの組み合わせ
   mx νx + my νy + mz νz が整数に近い順 (ADDDENSE で分解能) に追加する (`wfbin`)。
 - `wjfit`: 反復の最小二乗 (収束 1e-10、最大 30 回) で係数 b を求め、2J の rms・最大値と、係数の大きい 16 本を書く。

FFS からは呼べない (MAIN の TRACK だけ)。gs も MAIN の TRACK は SIGSEGV (tracka → tinip) で参照が無いので、
K64 の Fortran 自体 (wjfit・wfres・wfbin・wiord を gfortran で組んだ独立な基準。tfstk に依る else 枝だけを外した)
と比べる (docs/46、試験 test_fourier_fit.py)。未確認: 失った粒子の zbuf の扱い (K64 は失った粒子も 0 を含む全長で当てはめる)、
|mz| の項が lcoeff を超える場合 (ADDTERMS が分解能に比べて大きすぎるとき、K64 は iord の未初期化の値を使う → ここは 0 の組)。
"""
from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np

NBIG = 16
CONV = 1.0e-10
ITMAX = 30
LABEL = ('  |mx|+|my|   ', '|mx|+|my|+|mz|')


def lcoeff_of(maxord: int, nres: int, rfsw: bool) -> int:
    """tracka.f の lcoeff (係数の数 + J の数。J は 2 または 3 本)"""
    if rfsw:
        return ((maxord + 1) * (3 + maxord * (1 + maxord * 2))) // 3 + 2 + nres
    return maxord * (maxord + 1) + 2 + nres


def wiord(maxord: int, rfsw: bool, lcoeff: int | None = None) -> np.ndarray:
    """wiord.f: iord (3, lcoeff) (列 0.. は Fortran の l = 1..)。先頭の ndim 本は J の項 (0, 0, 0)"""
    rows: List[Tuple[int, int, int]] = [(0, 0, 0), (0, 0, 0)]
    if rfsw:
        rows.append((0, 0, 0))
    for m in range(1, maxord + 1, 2):
        for j in range(m, -m, -1):                                  # j = m, ..., 1-m
            rows.append((j, m - abs(j), 0))
        if rfsw:
            for i in range(1, m + 1):
                jm = m - i
                for j in range(-jm, jm + 1):
                    rows.append((j, jm - abs(j), i))
                for j in range(jm - 1, -jm, -1):                    # jm-1, ..., 1-jm
                    rows.append((j, abs(j) - jm, i))
        m1 = m + 1
        if m1 > maxord:
            continue
        if rfsw:
            for i in range(m1, 0, -1):
                jm = m1 - i
                for j in range(1 - jm, jm):
                    rows.append((j, abs(j) - jm, i))
                for j in range(jm, -jm - 1, -1):
                    rows.append((j, jm - abs(j), i))
        for j in range(1 - m1, m1 + 1):
            rows.append((j, m1 - abs(j), 0))
    n = len(rows) if lcoeff is None else max(lcoeff, len(rows))
    iord = np.zeros((3, n), dtype=np.int64)
    iord[:, :len(rows)] = np.array(rows, dtype=np.int64).T
    return iord


def wtune(z: np.ndarray) -> float:
    """wtune.f: FFT (exp(+2πi nk/N)) のパワー最大の 3 点の放物線補間。[-.5, .5] の tune 単位"""
    lsp = z.shape[0]
    X = np.fft.ifft(z) * lsp
    a = X.real ** 2 + X.imag ** 2
    im = int(np.argmax(a))                                          # 最初の最大 (Fortran は a > am の最初)
    am = a[im]
    a1 = a[(im - 1) % lsp]
    a2 = a[(im + 1) % lsp]
    with np.errstate(all='ignore'):
        w = (im + (a2 - a1) / 2.0 / (am * 2.0 - a1 - a2)) / lsp
    if w > 0.5:
        w -= 1.0
    return float(w)


def _wfbin(f1, f2, f3, n1, n2, n3, fmax, ip, iord) -> int:
    """wfbin.f: 組 (n1, n2, n3) の n1 f1 + n2 f2 + n3 f3 が整数に fmax 以内なら iord の ip 番目に入れて ip + 1 を返す"""
    f = n1 * f1 + n2 * f2 + n3 * f3
    df = abs(f - math.trunc(f))
    if df > 0.5:
        df = 1.0 - df
    if df < fmax:
        if ip < iord.shape[1]:
            iord[0, ip] = n1
            iord[1, ip] = n2
            iord[2, ip] = n3
        ip += 1
    return ip


def wfres(z: np.ndarray, minord: int, lcoeff: int, nres: int, iord: np.ndarray, aord: np.ndarray,
          ndim: int, adddense: float):
    """wfres.f: チューン (f1, f2, f3) を返し、nres > 0 なら iord の末尾 nres 本を整数に近い組で埋める (ip は 0 始まり)"""
    f1 = wtune(z[:, 0])
    f2 = wtune(z[:, 1])
    f3 = wtune(z[:, 2]) if ndim >= 3 else 0.0
    if nres <= 0:
        return f1, f2, f3
    lsp = z.shape[0]
    fmax = nres / adddense / lsp
    ip0 = lcoeff - nres                                             # Fortran の ip = lcoeff-nres+1 (1 始まり) → 0 始まり
    ip = [ip0]
    done = False

    def put(n1, n2, n3):
        nonlocal done
        ip[0] = _wfbin(f1, f2, f3, n1, n2, n3, fmax, ip[0], iord)
        if ip[0] >= lcoeff:                                         # Fortran: ip > lcoeff (1 始まり) ⇔ ip0 >= lcoeff
            done = True
        return done

    maxord = 100
    i = minord
    while i <= maxord and not done:
        for k in range(i, -i, -1):                                  # k = i .. 1-i
            if put(k, i - abs(k), 0):
                break
        if done:
            break
        i1 = i + 1
        if ndim >= 3:
            for j in range(1, i + 1):
                jm = i - j
                for k in range(-jm, jm + 1):
                    if put(k, jm - abs(k), j):
                        break
                if done:
                    break
                for k in range(jm - 1, -jm, -1):
                    if put(k, abs(k) - jm, j):
                        break
                if done:
                    break
            if done:
                break
            for j in range(i1, 0, -1):
                jm = i1 - j
                for k in range(1 - jm, jm):
                    if put(k, abs(k) - jm, j):
                        break
                if done:
                    break
                for k in range(jm, -jm - 1, -1):
                    if put(k, jm - abs(k), j):
                        break
                if done:
                    break
            if done:
                break
        for k in range(1 - i1, i1 + 1):
            if put(k, i - abs(k), 0):
                break
        if done:
            break
        i += 2
    # 足りなかった分 (ip < lcoeff) は iord に残った値 (K64 は未初期化) をそのまま aord に写す: ここでは 0 の組
    for l in range(lcoeff - nres, lcoeff):
        aord[:, l] = iord[:, l]
    return f1, f2, f3


def fit_particle(z: np.ndarray, iord: np.ndarray, aord: np.ndarray, lcoeff: int, ndim: int, rfsw: bool):
    """wjfit.f の当てはめ (1 粒子)。z (lsp, 3) 複素。戻り: b (lv)・zave・sj (3)・srms (3)・smax (3)・big (16 本の l の 0 始まり添字)"""
    lsp = z.shape[0]
    lv = (lcoeff - ndim) * 2 + ndim
    b = np.zeros(lv)
    zave = z.mean(axis=0)
    zexp = np.ones((lsp, 3), dtype=complex)
    ajx = np.zeros((lsp, 3))
    for k in range(ndim):
        d = z[:, k] - zave[k]
        aj = np.abs(d)
        nz = aj != 0.0
        zexp[nz, k] = d[nz] / aj[nz]
        ajx[:, k] = aj ** 2
    if not rfsw:
        ajx[:, 2] = 0.0
        zave[2] = z[:, 2].sum() / lsp
    it = 0
    with np.errstate(all='ignore'):
        while True:
            r = 0.0
            s = 0.0
            for l in range(ndim):
                u = ajx[:, l].sum() / lsp
                b[l] += u
                r += u * u
                s += b[l] ** 2
                ajx[:, l] -= u
            aexp = np.ones(lsp, dtype=complex)
            for l in range(ndim, lcoeff):                           # Fortran の l = ndim+1 .. lcoeff
                ll = 2 * (l + 1) - ndim - 1 - 1                     # b(ll) の 0 始まり
                for k in range(ndim):
                    idd = int(iord[k, l] - iord[k, l - 1])
                    if idd == 0:
                        continue
                    if idd == 1:
                        aexp = aexp * zexp[:, k]
                    elif idd == -1:
                        aexp = aexp * np.conj(zexp[:, k])
                    elif idd > 0:
                        aexp = aexp * zexp[:, k] ** idd
                    else:
                        aexp = aexp * np.conj(zexp[:, k]) ** (-idd)
                a1, a2, a3 = aord[0, l], aord[1, l], aord[2, l]
                f = a1 * ajx[:, 0] + a2 * ajx[:, 1] + a3 * ajx[:, 2]
                im = aexp.imag
                re = aexp.real
                p1 = float(np.sum(im ** 2))
                u1 = -float(np.sum(f * im))
                p2 = float(np.sum(re ** 2))
                u2 = -float(np.sum(f * re))
                p3 = float(np.sum(re * im))
                ff = (a1 * a1 + a2 * a2 + a3 * a3) * (p1 * p2 - p3 * p3)
                u = (p2 * u1 - p3 * u2) / ff
                u2 = (p1 * u2 - p3 * u1) / ff
                u1 = u
                b[ll] += u1
                b[ll + 1] += u2
                r += u1 * u1 + u2 * u2
                s += b[ll] ** 2 + b[ll + 1] ** 2
                f = u1 * im + u2 * re
                ajx[:, 0] += a1 * f
                ajx[:, 1] += a2 * f
                ajx[:, 2] += a3 * f
            if math.sqrt(r / s) < CONV:
                break
            it += 1
            if it >= ITMAX:
                break
        bigb = [0.0] * NBIG
        ibigb = [-1] * NBIG
        s1 = b[0]
        s2 = b[1]
        s3 = b[2] if rfsw else 0.0
        r1 = r2 = r3 = 0.0
        for l in range(ndim, lcoeff):
            ll = 2 * (l + 1) - ndim - 1 - 1
            s = math.sqrt(b[ll] ** 2 + b[ll + 1] ** 2)
            s1 += abs(aord[0, l] * s)
            s2 += abs(aord[1, l] * s)
            s3 += abs(aord[2, l] * s)
            r1 += (aord[0, l] * s) ** 2
            r2 += (aord[1, l] * s) ** 2
            r3 += (aord[2, l] * s) ** 2
            s = s * (aord[0, l] ** 2 + aord[1, l] ** 2 + aord[2, l] ** 2)
            for k in range(NBIG):
                if s >= bigb[k]:
                    for kk in range(NBIG - 1, k, -1):
                        bigb[kk] = bigb[kk - 1]
                        ibigb[kk] = ibigb[kk - 1]
                    bigb[k] = s
                    ibigb[k] = l
                    break
        smax = np.array([s1, s2, s3])
        srms = np.sqrt(np.array([r1, r2, r3]) * 0.5)
        sj = np.sqrt((ajx ** 2).sum(axis=0) / lsp)
    return b, zave, sj, srms, smax, ibigb


def _g(x, w, d):
    from .maintrack import fmt_g
    return fmt_g(float(x), w, d)


def fourier_report(zbuf: np.ndarray, lost: np.ndarray, maxord: int, nres: int, rfsw: bool, adddense: float = 6.0):
    """wjfit.f (maxord > 0) の全部。zbuf (np0, lsp, 3) 複素、lost (np0) 真なら失った粒子。
    戻り (行のリスト, 粒子ごとの辞書: tune・b・iord・sj・srms・smax・zave・big)"""
    np0, lsp = zbuf.shape[0], zbuf.shape[1]
    ndim = 3 if rfsw else 2
    lcoeff = lcoeff_of(maxord, nres, rfsw)
    iord0 = wiord(maxord, rfsw, lcoeff)
    lv = (lcoeff - ndim) * 2 + ndim
    lines = ['']
    lines.append(' <Fourier analysis>      ' + LABEL[ndim - 2] + '<=' + f'{maxord:3d}' + '    fitting variables:' + f'{lv - nres * 2:5d}')
    lines.append('                 ' + f'{maxord + 1:5d}' + ' <=' + LABEL[ndim - 2] + '     ' + '    fitting variables:' + f'{nres * 2:5d}')
    recs = []
    for i in range(np0):
        tag = ',lost:' if lost[i] else ':     '
        lines.append('   Particle #' + f'{i + 1:3d}' + tag)
        iord = iord0.copy()
        aord = iord.astype(float)
        f1, f2, f3 = wfres(zbuf[i], maxord + 1, lcoeff, nres, iord, aord, ndim, adddense)
        lines.append('          Tune = ' + ''.join(f'{v:15.7f}' for v in (f1, f2, f3)))
        lines.append(' Last variable = (' + ''.join(f'{int(iord[k, lcoeff - 1]):4d}' for k in range(3)) + ' )')
        b, zave, sj, srms, smax, ibigb = fit_particle(zbuf[i], iord, aord, lcoeff, ndim, rfsw)
        lines.append('      <X>        <PX>         <Y>        <PY>         <Z>         <P>')
        zv = []
        for k in range(3):
            zv += [zave[k].real, zave[k].imag]
        lines.append(''.join(_g(v, 12, 4) for v in zv))
        lines.append('         2J        rms deviation     rms d2J      maximum 2J')
        for k, nm in enumerate('xyz'[:ndim]):
            lines.append(' ' + nm + ''.join(_g(v, 15, 7) for v in (b[k], sj[k], srms[k], smax[k])))
        if ndim == 2:
            lines.append('')                                        # 書式 9101 の ':' の前の '/' が空行を出す
        lines.append('   Big 16 coefficients:')
        lines.append('  mx  my  mz       coefficients       mx  my  mz       coefficients      ')
        row = ''
        for k in range(NBIG):
            l = ibigb[k]
            if l < 0:
                row += f'{0:4d}{0:4d}{0:4d}' + _g(0.0, 12, 4) * 2
            else:
                ll = 2 * (l + 1) - ndim - 1 - 1
                row += f'{int(iord[0, l]):4d}{int(iord[1, l]):4d}{int(iord[2, l]):4d}' + _g(b[ll], 12, 4) + _g(b[ll + 1], 12, 4)
            if k % 2 == 1:
                lines.append(row)
                row = ''
        recs.append({'tune': (f1, f2, f3), 'b': b, 'iord': iord[:, :lcoeff].copy(), 'sj': sj, 'srms': srms,
                     'smax': smax, 'zave': zave, 'big': ibigb})
    return lines, recs
