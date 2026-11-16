"""
MEASURE / DRAW / TRCOD の統計の部品 (tlum.f)。K64 の tlum / tluma の直訳 (ルミノシティの推定: 二つのビームの粒子数の分布の重なり)。
ttstat.f (MEASURE の統計) から呼ばれる。FFS のコマンドとしての MEASURE / TRCOD (tracke.f、ttinit.f、ttstat.f、trcoda.f) は
gs が動かない (MEASURE は何も出力せず、TRCOD は SIGSEGV で落ちる: docs/19 B42) ので、翻訳していない (参照が取れない)。
bin は real*4 (K64 の `real*4 bin(2,0:nmesh,0:nmesh), alum`)。足し込みは単精度に丸める。
"""
from __future__ import annotations

import math

import numpy as np


def tluma(x, y, xm, ym, sigx, sigy0, nmesh):
    """tluma: x[j, i], y[j, i] (j=0,1: 二つのビーム、i: 粒子)。戻り tluma"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    np_ = x.shape[1]
    s = nmesh / 6.0
    sigy = float(sigy0)
    alum0 = np.float32(0.0)
    f32 = np.float32
    while True:
        bin_ = np.zeros((2, nmesh + 1, nmesh + 1), dtype=np.float32)
        for j in range(2):
            for i in range(np_):
                xa = ((x[j, i] - xm) / sigx + 3.0) * s
                ya = ((y[j, i] - ym) / sigy + 3.0) * s
                if xa > 0.0 and ya > 0.0:
                    ix = int(xa)
                    iy = int(ya)
                    if ix < nmesh and iy < nmesh:
                        fx = xa - ix
                        fy = ya - iy
                        bin_[j, ix, iy] = f32(bin_[j, ix, iy] + (1.0 - fx) * (1.0 - fy))
                        bin_[j, ix + 1, iy] = f32(bin_[j, ix + 1, iy] + fx * (1.0 - fy))
                        bin_[j, ix, iy + 1] = f32(bin_[j, ix, iy + 1] + (1.0 - fx) * fy)
                        bin_[j, ix + 1, iy + 1] = f32(bin_[j, ix + 1, iy + 1] + fx * fy)
        alum = alum0
        for i in range(nmesh + 1):
            for j in range(nmesh + 1):
                alum = f32(alum + f32(bin_[0, j, i]) * f32(bin_[1, j, i]))
        tl = np_**2 / float(alum) * sigx * sigy * (6.0 / nmesh)**2 / 4.0 / math.pi
        if tl < 0.8 * sigx * sigy:
            sigy1 = tl / sigx
            scale = max(0.8, sigy1 / sigy)
            j0 = int((-scale * 3.0 + 3.0) * s)
            j1 = int((scale * 3.0 + 3.0) * s)
            sigy1 = (sigy * (j1 - j0)) / nmesh
            alum = alum0
            for i in range(nmesh + 1):
                for j in range(0, j0):
                    alum = f32(alum + f32(bin_[0, i, j]) * f32(bin_[1, i, j]))
                for j in range(j1 + 1, nmesh + 1):
                    alum = f32(alum + f32(bin_[0, i, j]) * f32(bin_[1, i, j]))
            alum0 = f32(float(alum) * sigy1 / sigy)
            sigy = sigy1
            continue
        return tl


def tlum(x, y, xm, ym, sigx, sigy):
    """tlum: nmesh = max(int(sqrt(np)), 6)。戻り (tlum, nmesh)"""
    np_ = np.asarray(x).shape[1]
    nmesh = max(int(math.sqrt(np_)), 6)
    return tluma(x, y, xm, ym, sigx, sigy, nmesh), nmesh
