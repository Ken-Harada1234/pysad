# -*- coding: utf-8 -*-
"""保存量とシンプレクティック性 (04_quad.md §6.1、式 (4.18)(4.19))。

四極の本体の H (z を除く、p = 1、dv = 0) は H = -pz + K (x^2 - y^2)/2 (pz = sqrt(1 - px^2 - py^2))。H は s に陽に依らないので、厳密な流れでは保存する。
(a) B だけ (厳密な線形の流れ) は H2 = (px^2 + py^2)/2 + K (x^2 - y^2)/2 を、D だけは ΔH を、それぞれ厳密に (丸め誤差まで) 保存する。
(b) 分割した写像では H = H2 + ΔH は保存せず、(1 台を通ったときの) 変化が刻みの数 nd とともに Strang は nd^-2、Yoshida は nd^-4 で減る。
(c) 分割した写像の 1 刻みの 4x4 ヤコビアンの行列式とシンプレクティック条件 (D も B も厳密にシンプレクティックなので、合成も丸め誤差まで)。
"""
import math
import numpy as np
import _common as C

C.start(__file__)
from pysad import track6 as t6, quadbody as qb

L, K1L = 1.0, 1.5
K = K1L / L
c0 = (0.01, 0.01, -0.008, 0.008)


def H_full(x, px, y, py):
    return -math.sqrt(1 - px * px - py * py) + K * (x * x - y * y) / 2


def H2(x, px, y, py):
    return (px * px + py * py) / 2 + K * (x * x - y * y) / 2


def dH(px, py):
    return 1 - math.sqrt(1 - px * px - py * py) - (px * px + py * py) / 2


print('== (a) B だけ・D だけの保存量 (長さ 0.3 m の 1 回) ==')
x, px, y, py = c0
par = t6._tzsetparam0(0.0, 0.3, K)
xb, pxb, yb, pyb, _z = t6._qb_body(x, px, y, py, 0.0, 0.0, par, 0.3, 0.0, 0.0, False)
print('  B: H2 の変化 %.2e (H2 = %.6e)' % (H2(xb, pxb, yb, pyb) - H2(x, px, y, py), H2(x, px, y, py)))
xd, yd, _z = t6._qb_drift(x, px, y, py, 0.0, 0.3)
print('  D: px, py は変えない (x, y, z だけ動く) ので ΔH(px, py) = %.6e は厳密に保存' % dH(px, py))
print('  D の前後で H2 (x を使う) は %.3e 変わる: D は H2 を保存しない (H2 と ΔH は交換しない)' % (H2(xd, px, yd, py) - H2(x, px, y, py)))

print('\n== (b) 1 台 (L = %g m、K1L = %g) を通ったときの H の変化 (厳密な流れでは 0) ==' % (L, K1L))
print('   nd    Strang |ΔH|    比       Yoshida |ΔH|   比')


def through(c, nd, sch):
    xx, pxx, yy, pyy = c
    zz = 0.0
    akk = K1L / L
    aln = L / nd
    pa = t6._tzsetparam0(0.0, aln * qb.QB_BW[sch][0], akk)
    pb = t6._tzsetparam0(0.0, aln * qb.QB_BW[sch][1], akk) if sch == 1 else pa
    m = qb.QB_M[sch]
    for g in range(m * nd):
        dl, bl, jb, _r, _f = qb.group(g, aln, nd, sch)
        xx, yy, zz = t6._qb_drift(xx, pxx, yy, pyy, zz, dl)
        xx, pxx, yy, pyy, zz = t6._qb_body(xx, pxx, yy, pyy, zz, 0.0, pb if jb == 1 else pa, bl, 0.0, 0.0, False)
    dl, _b, _j, _r, _f = qb.group(m * nd, aln, nd, sch)
    xx, yy, zz = t6._qb_drift(xx, pxx, yy, pyy, zz, dl)
    return xx, pxx, yy, pyy


prev = [None, None]
for nd in (2, 4, 8, 16, 32):
    row = []
    for sch in (0, 1):
        o = through(c0, nd, sch)
        d = abs(H_full(*o) - H_full(*c0))
        row.append((d, '%.2f' % (prev[sch] / d) if prev[sch] else '   -'))
        prev[sch] = d
    print('  %3d    %.3e   %s    %.3e    %s' % (nd, row[0][0], row[0][1], row[1][0], row[1][1]))

print('\n== (c) 1 台のヤコビアン (4x4、中心差分) のシンプレクティック性: |M^T J M - J| ==')
J4 = np.zeros((4, 4))
J4[0, 1] = J4[2, 3] = 1.0
J4 -= J4.T
for sch, name in ((0, 'Strang (nd = 4)'), (1, 'Yoshida (nd = 4)')):
    M = np.zeros((4, 4))
    h = 1e-6
    for j in range(4):
        cp, cm = list(c0), list(c0)
        cp[j] += h
        cm[j] -= h
        M[:, j] = (np.array(through(cp, 4, sch)) - np.array(through(cm, 4, sch))) / (2 * h)
    print('  %s: |M^T J M - J| = %.1e、det M = %.15f' % (name, np.abs(M.T @ J4 @ M - J4).max(), np.linalg.det(M)))
