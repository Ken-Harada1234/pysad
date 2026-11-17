"""
MEASURE / TRCOD / LOGMEASURE (FFS のコマンド) と、その部品 (tlum.f、ttinit.f、ttstat.f、tbdecoup.f、tracke.f、trackb.f、trcoda.f) の翻訳。

K64 (gs) では MEASURE は何も出力せず、TRCOD は SIGSEGV で落ちる (docs/21 B42) ので、gs の参照は取れない。
ここの実装は「ソースの読みのみ」で、値は解析式 (σ=√(εβ) など) と pysad の独立な経路 (CALC の Twiss、TrackParticles) で確かめた
(pysad/tests/test_v2_C_measure.py、docs/history/v2_changes_C_measure.md)。

- tlum / tluma : tlum.f の直訳 (二つのビームの重なりの実効面積)。bin は real*4 (K64 の `real*4 bin(2,0:nmesh,0:nmesh), alum`)。
- tbdecoup     : tbdecoup.f の直訳 (6×6 のビーム行列 → 射影 / 非結合の emittance・β・α・分散・R1..R4)。
- ttinit       : ttinit.f (先頭 MARK の Twiss・EMITX/EMITY・DP・SIGZ・AZ から初期分布を作る)。正準変数 (x, px/p0, y, py/p0, z, dp/p0) を返す。
                 b41 / center_fix の既定は互換モードで決まる (v2: 修正、k64: K64 のソースどおり)。
- ttstat       : ttstat.f (重心と 6×6 の分散行列)。ttstat_lines は表示 (書式 9001〜9006)。
- cmd_measure / cmd_trcod / cmd_logmeasure : FFS のコマンド (tffsa.f)。**v2 モードだけ** ffscmd から呼ばれる
  (k64 モードは従来どおり「未対応 (無視)」: gs が何も出さない・落ちるので、k64 の動作は変えていない)。
- trcod_fix / trcod_move / trcod_restore : TRCOD FIXCOD / MVCOD / STNCOD (trcoda.f の entry trcodf / trcodm / trcods)。基準軌道の保存。
- cmd_align    : TRCOD ALIGN (nalign.f)。区間の磁石の DX・DY を先頭の磁石の値に揃える。
- cmd_putmeasure / cmd_getmeasure : PUTMEASURE / GETMEASURE (nfgetm.c、tfgetm.f)。測定値を環境変数 SADPOUT / SADPIN のファイルで受け渡す。
- cmd_orbitjitter : ORBITJITTER x y (tfojit.f)。Twiss の DX・DY に要素ごとの正規乱数を足す。
  (以上は v2 だけ。記録は docs/history/v2_changes_C2_measure.md、試験は pysad/tests/test_v2_C2_measure.py)
- track_trpt   : MAIN の TRACK で TRPT ON のとき (tracka.f の trpt の枝: ttinit → Entrance の統計 → 追跡 → Exit の統計)。v2 だけ。
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

from . import compat


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


# =============================================================================================
# tbdecoup.f
# =============================================================================================
def tbdecoup(beam):
    """tbdecoup.f の直訳。beam: 6×6 の分散行列 (x, px, y, py, z, dp)。書き換えない (写しを使う)。
    戻り (param, param1, emix, emiy, emix1, emiy1, sigp)。param / param1 は 1 始まり (添字 0 は使わない、長さ 22):
      param(1)=ax (2)=bx (4)=ay (5)=by (7..10)=ex epx ey epy   … x-y に射影した (結合を無視した) 値
      param1(同上) … x-y を非結合にした値 (u, v)。param1(11..14) = r1..r4
    emix, emiy / emix1, emiy1: 射影 / 非結合の emittance (分散の寄与を引いた 2×2 の行列式の平方根。符号つき)
    sigp = beam(6,6) (dp の分散。標準偏差ではない)"""
    b = np.array(beam, dtype=float)
    param = np.zeros(22)
    param1 = np.zeros(22)
    sigp = b[5, 5]
    if sigp != 0:
        for k in range(4):
            param[7 + k] = b[5, k] / sigp
    for i in range(4):
        for j in range(i + 1):
            b[i, j] = b[i, j] - param[7 + i] * param[7 + j] * sigp
            b[j, i] = b[i, j]

    def _twiss2(bm, o, par):
        d = bm[o, o] * bm[o + 1, o + 1] - bm[o + 1, o] ** 2
        em = math.copysign(math.sqrt(abs(d)), d)
        k = 1 if o == 0 else 4
        if em != 0:
            par[k] = -bm[o + 1, o] / em
            par[k + 1] = bm[o, o] / em
        else:
            par[k] = 0.0
            par[k + 1] = 0.0
        return em

    emix = _twiss2(b, 0, param)
    emiy = _twiss2(b, 2, param)
    p1, p2, p4 = b[1, 1], b[1, 0], b[0, 0]
    q1, q2, q4 = b[3, 3], b[3, 2], b[2, 2]
    m1, m2, m3, m4 = b[3, 1], b[3, 0], b[2, 1], b[2, 0]
    r1 = -m4 * p1 + m3 * p2 + m3 * q2 - m1 * q4
    r2 = m4 * p2 - m3 * p4 - m4 * q2 + m2 * q4
    r3 = -m2 * p1 + m1 * p2 + m3 * q1 - m1 * q2
    r4 = m2 * p2 - m1 * p4 - m4 * q1 + m2 * q2
    detr = r1 * r4 - r2 * r3
    detp = p1 * p4 - p2 ** 2
    detq = q1 * q4 - q2 ** 2
    if detp == detq:
        return param, param1, emix, emiy, 0.0, 0.0, sigp
    f = detr / (detp - detq) ** 2
    g = 4.0 * f + 1.0
    if g < 0.0:
        return param, param1, emix, emiy, 0.0, 0.0, sigp
    a = math.sqrt((g + math.sqrt(g)) / g * 0.5)
    c = (2.0 * a ** 2 - 1.0) / a / (detp - detq)
    r1, r2, r3, r4 = c * r1, c * r2, c * r3, c * r4
    param1[11:15] = (r1, r2, r3, r4)
    t = np.eye(6)
    t[0, 0] = t[1, 1] = t[2, 2] = t[3, 3] = a
    t[0, 2], t[0, 3], t[1, 2], t[1, 3] = -r4, r2, r3, -r1
    t[2, 0], t[2, 1], t[3, 0], t[3, 1] = r1, r2, r3, r4
    b = t @ b
    t[0, 2], t[0, 3], t[1, 2], t[1, 3] = r1, r3, r2, r4
    t[2, 0], t[2, 1], t[3, 0], t[3, 1] = -r4, r3, r2, -r1
    b = b @ t
    emix1 = _twiss2(b, 0, param1)
    emiy1 = _twiss2(b, 2, param1)
    param1[7] = a * param[7] - r4 * param[9] + r2 * param[10]
    param1[8] = a * param[8] + r3 * param[9] - r1 * param[10]
    param1[9] = a * param[9] + r1 * param[7] + r2 * param[8]
    param1[10] = a * param[10] + r3 * param[7] + r4 * param[8]
    return param, param1, emix, emiy, emix1, emiy1, sigp


# =============================================================================================
# ttinit.f
# =============================================================================================
def ttinit(mark, np0, rng, gcut=1.0e300, *, gauss=False, trgauss=True, jitter=True, dp0=0.0, wake=False,
           b41=None, center_fix=None):
    """ttinit.f: 先頭 MARK の値から np0 個の粒子の初期分布を作る。戻り (np0, 6) の配列 (x, px/p0, y, py/p0, z, dp/p0: 正準変数。
    K64 は最後に tconvm で px/p に直すが、pysad の追跡 (TrackParticles) は正準変数を受けるので直さない)。

    mark: dict。AX BX AY BY EX EPX EY EPY R1 R2 R3 R4 DETR DX DPX DY DPY AZ (先頭 MARK の値)、
          EMITX EMITY SIGZ DP (tfsetparam が MARK に写す大域変数。DP は dpmax)、JDX JDPX JDY JDPY (軌道の揺らぎ)、
          ZPX ZY ZPY (K64 の B41 の読み方を再現するときだけ使う)。無いキーは 0 (BX BY は 1)。
    rng : sadrandom.SadRandom (tgauss = gauss(1, gcut)、tran = uniform0)。乱数を引く順は K64 と同じ。
    gauss  : GAUSS (運動量を正規分布にする。OFF なら ±DP の一様分布)。
    trgauss: TRGAUSS (横方向を正規分布。OFF なら単位円内の一様分布 × σ: 実際の rms は σ/2)。
    jitter : JITTER。OFF (QUIET) なら、分布の重心の統計的なずれを引いて重心を設計値に合わせる。
    dp0    : 運動量のずれの中心 (DP0 = 先頭 MARK の DDP)。
    wake   : True なら z の乱数 pk を 0 にする (K64 の TWAKE/LWAKE ON のときの動作: 全粒子 z=0。wake 用の分布はコメントアウトされている)。
    b41    : True なら K64 のとおり、軌道の揺らぎを rlist(kx+26..29) = ZPX ZY ZPY DP から読む (B41)。False なら JDX JDPX JDY JDPY。
             None なら互換モードで決める (k64: True、v2: False)。
    center_fix: JITTER OFF の重心合わせの基準。False なら K64 の式 (x は ex·dp0+dx、px は epx·dp0 − ax/bx·ex·dp0 + dpx。
             px の第 2 項は分布の作り方 (px の平均 = epx·dp0 + dpx) と合わず、結合 R1..R4 も入っていない)。
             True なら分布の作り方と同じ写像で求めた設計の重心。None なら互換モードで決める (k64: False、v2: True)。"""
    if b41 is None:
        b41 = compat.k64()
    if center_fix is None:
        center_fix = compat.v2()

    def g_(k, d=0.0):
        v = mark.get(k, None)
        return d if v is None else float(v)

    axi, ayi = g_('AX'), g_('AY')
    bxi, byi = g_('BX', 1.0), g_('BY', 1.0)
    emx, emy, dpmax = g_('EMITX'), g_('EMITY'), g_('DP')
    exi, epxi, eyi, epyi = g_('EX'), g_('EPX'), g_('EY'), g_('EPY')
    r1i, r2i, r3i, r4i = g_('R1'), g_('R2'), g_('R3'), g_('R4')
    detr = g_('DETR')
    amu = math.sqrt(1.0 - (r1i * r4i - r2i * r3i))
    dxi, dpxi, dyi, dpyi = g_('DX'), g_('DPX'), g_('DY'), g_('DPY')
    sigz, azi = g_('SIGZ'), g_('AZ')
    if b41:
        dxj, dpxj, dyj, dpyj = g_('ZPX'), g_('ZY'), g_('ZPY'), dpmax        # rlist(kx+26..29)
    else:
        dxj, dpxj, dyj, dpyj = g_('JDX'), g_('JDPX'), g_('JDY'), g_('JDPY')  # rlist(kx+34..37)
    tg = lambda: rng.gauss(1, gcut)
    tr = rng.uniform0
    dxi += dxj * tg()
    dpxi += dpxj * tg()
    dyi += dyj * tg()
    dpyi += dpyj * tg()
    sigx = math.sqrt(bxi * emx)
    sigpx = math.sqrt(emx / bxi)
    sigy = math.sqrt(byi * emy)
    sigpy = math.sqrt(emy / byi)
    rz = math.sqrt(1.0 + azi ** 2)

    def couple(xi0, pxi0, yi0, pyi0):
        a_ = amu * xi0 + r4i * yi0 - r2i * pyi0
        b_ = amu * pxi0 - r3i * yi0 + r1i * pyi0
        c_ = amu * yi0 - r1i * xi0 - r2i * pxi0
        d_ = amu * pyi0 - r3i * xi0 - r4i * pxi0
        if detr <= 1.0:
            return a_ + dxi, b_ + dpxi, c_ + dyi, d_ + dpyi
        return c_ + dxi, d_ + dpxi, a_ + dyi, b_ + dpyi

    out = np.zeros((np0, 6))
    for i in range(np0):
        pk = tg()
        if wake:
            pk = 0.0
        if gauss:
            pj = dp0 + dpmax * (tg() - azi * pk) / rz
        else:
            pj = dp0 + dpmax * ((tr() - 0.5) * 2.0 - azi * pk) / rz
        if trgauss:
            xi = sigx * tg()
            xi0 = xi + exi * pj
            pxi0 = sigpx * tg() - axi / bxi * xi + epxi * pj
            yi = sigy * tg()
            yi0 = yi + eyi * pj
            pyi0 = sigpy * tg() - ayi / byi * yi + epyi * pj
        else:
            while True:
                q1 = 2.0 * tr() - 1.0
                q2 = 2.0 * tr() - 1.0
                if q1 ** 2 + q2 ** 2 <= 1.0:
                    break
            xi = sigx * q1
            xi0 = xi + exi * pj
            pxi0 = sigpx * q2 - axi / bxi * xi + epxi * pj
            while True:
                q1 = 2.0 * tr() - 1.0
                q2 = 2.0 * tr() - 1.0
                if q1 ** 2 + q2 ** 2 <= 1.0:
                    break
            yi = sigy * q1
            yi0 = yi + eyi * pj
            pyi0 = sigpy * q2 - ayi / byi * yi + epyi * pj
        out[i, 0], out[i, 1], out[i, 2], out[i, 3] = couple(xi0, pxi0, yi0, pyi0)
        out[i, 4] = sigz * pk
        out[i, 5] = pj
    if not jitter and np0 > 0:
        m = out.sum(axis=0) / np0
        if center_fix:
            tx, tpx, ty, tpy = couple(exi * dp0, epxi * dp0, eyi * dp0, epyi * dp0)
        else:
            tx = exi * dp0 + dxi
            tpx = (epxi * dp0 - axi / bxi * exi * dp0) + dpxi
            ty = eyi * dp0 + dyi
            tpy = (epyi * dp0 - ayi / byi * eyi * dp0) + dpyi
        out[:, 5] -= m[5] - dp0
        out[:, 4] -= m[4]
        out[:, 0] -= m[0] - tx
        out[:, 1] -= m[1] - tpx
        out[:, 2] -= m[2] - ty
        out[:, 3] -= m[3] - tpy
    return out


# =============================================================================================
# ttstat.f
# =============================================================================================
_TAG = ('        x ', '    px/p0 ', '        y ', '    py/p0 ', '        z ', '    dp/p0 ')


def ttstat(coords):
    """ttstat.f の統計: coords (n, 6) = 正準変数 (x, px/p0, y, py/p0, z, dp/p0)。戻り (sa, ss): 重心 (6) と分散行列 (6×6、1/n で割る)。
    K64 は内部変数 (px/p) に (1+dp) を掛けて px/p0 に戻してから積算する。ここは始めから px/p0。"""
    c = np.asarray(coords, dtype=float)
    n = c.shape[0]
    sa = c.sum(axis=0) / n
    ss = (c.T @ c) / n - np.outer(sa, sa)
    return sa, ss


def _e103(v):
    """Fortran の 1P,E10.3 (sngl に丸めてから)"""
    with np.errstate(all='ignore'):
        v = float(np.float32(v))
    if v != v:
        return '       NaN'
    if math.isinf(v):
        return '  Infinity' if v > 0 else ' -Infinity'
    mant, e = f'{v:.3E}'.split('E')
    en = int(e)
    s = mant + (f'E{en:+03d}' if abs(en) < 100 else f'{en:+04d}')
    return s.rjust(10) if len(s) <= 10 else '*' * 10


def ttstat_lines(coords, title, info=None, mat=True, lum=False):
    """ttstat.f の表示 (書式 9001, 9002, 9003, 9006, 9004, 9005)。戻り (行のリスト, sa, ss, es)。
    info: {'RAD','RFSW','GAUSS','DP','DP0','GCUT'} (見出しの 2 行目)。lum=True なら実効面積 (tlum: 奇数番と偶数番の粒子を二つのビームとみなす)
    と垂直の 4・6・8 次のモーメントも出す (K64 の MEASURE・TRACK は lum=.false. で呼ぶので出ない)。es = 実効面積 (lum のときだけ)。"""
    from .maintrack import fmt_g
    from .sadscript import fmt_real
    c = np.asarray(coords, dtype=float)
    n = c.shape[0]
    sa, ss = ttstat(c)
    info = info or {}
    L = []
    es = None
    if title.strip() == '':
        return L, sa, ss, es
    if mat:
        tf = lambda b: 'T' if b else 'F'
        L.append(' Statistics at ' + title.rstrip() + ':' + '   particles =' + f'{n:5d}')
        L.append(' RAD: ' + tf(info.get('RAD')) + ', RFSW: ' + tf(info.get('RFSW')) + ', GAUSS: ' + tf(info.get('GAUSS'))
                 + ', DP =' + fmt_real(float(info.get('DP', 0.0)), '10.6') + ', DP0 =' + fmt_real(float(info.get('DP0', 0.0)), '10.6')
                 + ', GCUT =' + fmt_real(float(info.get('GCUT', 1e300)), '10.6'))
        L.append(' ' * 12 + ''.join(_TAG))
    L.append('   C of M ' + ': ' + ''.join(_e103(v) for v in sa))
    if mat:
        for i in range(6):
            L.append(_TAG[i] + ': ' + ''.join(_e103(ss[i, j]) for j in range(i + 1)))
    sigz = math.sqrt(ss[4, 4]) if ss[4, 4] > 0.0 else 0.0
    ecorr = ss[5, 4] / sigz if ss[4, 4] > 0.0 else 0.0
    p, p1, temix, temiy, temix1, temiy1, sigp = tbdecoup(ss)
    G = lambda v: fmt_g(float(v), 11, 4)
    sq = lambda v: math.sqrt(v) if v > 0.0 else 0.0
    L.append(' x-y projected(coupled) parameters:')
    L.append(' emitx:' + G(temix) + ' bx:' + G(p[2]) + ' ax:' + G(p[1]) + ' ex:' + G(p[7]) + ' epx:' + G(p[8]))
    L.append(' emity:' + G(temiy) + ' by:' + G(p[5]) + ' ay:' + G(p[4]) + ' ey:' + G(p[9]) + ' epy:' + G(p[10]))
    L.append(' x-y decoupled parameters:')
    L.append(' emitu:' + G(temix1) + ' bu:' + G(p1[2]) + ' au:' + G(p1[1]) + ' eu:' + G(p1[7]) + ' epu:' + G(p1[8]))
    L.append(' emitv:' + G(temiy1) + ' bv:' + G(p1[5]) + ' av:' + G(p1[4]) + ' ev:' + G(p1[9]) + ' epv:' + G(p1[10]))
    L.append('    r1:' + G(p1[11]) + ' r2:' + G(p1[12]) + ' r3:' + G(p1[13]) + ' r4:' + G(p1[14])
             + 'detr:' + G(p1[11] * p1[14] - p1[12] * p1[13]))
    L.append('   sigx:' + G(sq(ss[0, 0])) + '  sigy:' + G(sq(ss[2, 2])) + '   tilt:'
             + G(math.atan2(-2.0 * ss[2, 0], ss[0, 0] - ss[2, 2]) * 0.5))
    L.append('  sigpx:' + G(sq(ss[1, 1])) + ' sigpy:' + G(sq(ss[3, 3])))
    L.append(' sigp/p:' + G(sq(sigp)) + '  sigz:' + G(sigz) + ' dp/p/z:' + G(ecorr) + '/sigz')
    if lum and n >= 2 and ss[0, 0] > 0.0 and ss[2, 2] > 0.0:
        dy = (c[:, 2] - sa[2]) ** 2 / ss[2, 2]
        h = n // 2
        x2 = c[:2 * h, 0].reshape(h, 2).T            # x(2, np/2): 奇数番 = ビーム 1、偶数番 = ビーム 2
        y2 = c[:2 * h, 2].reshape(h, 2).T
        alum, nbin = tlum(x2, y2, sa[0], sa[2], math.sqrt(ss[0, 0]), math.sqrt(ss[2, 2]))
        L.append(' ')
        L.append('   Effective beam area:' + fmt_g(alum, 10, 4) + ' mc2' + '     (' + f'{nbin:3d}*{nbin:3d}' + ' bin)')
        L.append('   4-, 6-, 8th Vertical moments:'
                 + ''.join(fmt_g(float(v), 12, 4) for v in ((dy ** 2).sum() / n, (dy ** 3).sum() / n, (dy ** 4).sum() / n)))
        es = alum
    return L, sa, ss, es


# =============================================================================================
# FFS のコマンド (tffsa.f の MEA_SURE / TRC_OD / LOGM_EASURE、trackb.f、tracke.f、trcoda.f)
# =============================================================================================
def mark_values(ffs, main=False):
    """先頭 MARK の値 (ttinit が読む rlist(latt(1)+…))。FFS では tfsetparam のとおり EMITX EMITY SIGZ DP は大域変数から。
    main=True (MAIN の TRACK: tfsetparam を通らない) なら MARK の定義のキーワード EMITX EMITY SIGMAZ DP から。
    先頭が MARK でなければ None (ttinit: 'Missing entrance condition of beam.')"""
    if not ffs.els or type(ffs.els[0]).__name__ != 'Marker':
        return None
    m0 = ffs.els[0]
    raw = ffs.env._elements.get(m0.name, {})
    tw = dict(zip(ffs._MARK28, ffs._mark_tw28(ffs.dp0())))
    mk = {k: tw[k] for k in ('AX', 'BX', 'AY', 'BY', 'EX', 'EPX', 'EY', 'EPY', 'R1', 'R2', 'R3', 'R4', 'DETR',
                             'DX', 'DPX', 'DY', 'DPY', 'AZ', 'ZPX', 'ZY', 'ZPY')}
    for k in ('JDX', 'JDPX', 'JDY', 'JDPY'):
        try:
            mk[k] = float(raw.get(k, 0.0) or 0.0)
        except (TypeError, ValueError):
            mk[k] = 0.0
    g = ffs.I.globals

    def num(k):
        try:
            return float(g.get(k, 0.0) or 0.0)
        except (TypeError, ValueError):
            return 0.0

    if main:
        def rawnum(*keys):
            for k in keys:
                try:
                    if raw.get(k, None) is not None:
                        return float(raw.get(k))
                except (TypeError, ValueError):
                    pass
            return 0.0
        mk['EMITX'], mk['EMITY'], mk['SIGZ'], mk['DP'] = rawnum('EMITX'), rawnum('EMITY'), rawnum('SIGMAZ', 'SIGZ'), rawnum('DP')
        return mk
    mk['EMITX'], mk['EMITY'], mk['SIGZ'] = num('EMITX'), num('EMITY'), num('SIGZ')
    mk['DP'] = ffs.dp()
    return mk


def _np0(ffs):
    try:
        return int(float(ffs.I.globals.get('NP', 0.0) or 0.0))
    except (TypeError, ValueError):
        return 0


def _adjust_l2(ffs, l20):
    """tracker の冒頭: 指定の位置が MARK なら 1 つ先、直前が MARK なら 1 つ手前まで追う (MARK は恒等なので分布は同じ)"""
    nlat = ffs.nlat

    def ismark(k):
        return 1 <= k < nlat and type(ffs.els[k - 1]).__name__ == 'Marker'

    if l20 < nlat:
        if ismark(l20):
            return l20 + 1
        if l20 > 1 and ismark(l20 - 1):
            return l20 - 1
    return l20


class Tracker:
    """tracke.f tracker の状態 (STANDBY / TRACK / CONT / POS / DISP / RESET)。追跡は FFS の TrackParticles (fn_track) と同じ経路。
    粒子は (n, 6) の正準変数と生死のフラグで持つ。統計は生きている粒子だけ (K64 は落ちた粒子を配列の後ろへ回して np を減らす)。"""

    def __init__(self, ffs, jitter=None):
        self.ffs = ffs
        fl = ffs.flags
        rng = ffs.rng
        self._seed = (rng.seed, rng.filled, rng.next)          # NISTACK$FIXSEED@Push[]
        self.mark = mark_values(ffs)
        self.np0 = _np0(ffs)
        self.sa = np.zeros(6)
        self.ss = np.zeros((6, 6))
        # v2: TWAKE/LWAKE でも z の分布を作る (K64 は z=0 にするが、そのための wake 用の分布はコメントアウトされていて使えない)
        self.x0 = ttinit(self.mark, self.np0, rng, ffs.gcut(), gauss=bool(fl['GAUSS']), trgauss=bool(fl['TRGAUSS']),
                         jitter=bool(fl['JITTER']) if jitter is None else jitter, dp0=ffs.dp0(), wake=False)
        self.x = self.x0.copy()
        self.alive0 = np.ones(self.np0)
        self.alive = self.alive0.copy()
        self.l1 = 1
        self.l2s = 1

    def info(self):
        f = self.ffs
        return {'RAD': f.flags['RAD'], 'RFSW': f.flags['RFSW'], 'GAUSS': f.flags['GAUSS'], 'DP': f.dp(), 'DP0': f.dp0(),
                'GCUT': f.gcut()}

    def alive_x(self):
        return self.x[self.alive != 0.0]

    def track(self, l20, cont=True):
        """TRACK (cont=False: 保存した初期分布から) / CONT (前回の位置から) で l2 の入口まで追う"""
        l2 = _adjust_l2(self.ffs, l20)
        if not cont:
            self.x = self.x0.copy()
            self.alive = self.alive0.copy()
        else:
            self.l1 = self.l2s
        self.l2s = l2
        if l2 > self.l1:
            cols = [[float(v) for v in self.x[:, j]] for j in range(6)] + [[float(v) for v in self.alive]]
            res = self.ffs.fn_track([[float(self.l1), cols], float(l2)])
            r = res[1]
            self.x = np.array(r[:6], dtype=float).T
            self.alive = np.array(r[6], dtype=float)

    def stat(self, title, write=None):
        """ttstat (mat=.true., lum=.false.)。戻り (np, sv): sv = (xa, ya, xxa, yya, xya)。np<=1 なら sv は 1e20"""
        xs = self.alive_x()
        n = xs.shape[0]
        if n <= 1:
            return n, [1.0e20] * 5
        lines, self.sa, self.ss, _ = ttstat_lines(xs, title, self.info())
        if write is not None:
            for ln in lines:
                write(ln + chr(10))
        return n, [float(self.sa[0]), float(self.sa[2]), float(self.ss[0, 0]), float(self.ss[2, 2]), float(self.ss[0, 2])]

    def reset(self):
        """RESET と、FIXSEED のときの乱数の種の復元 (NISTACK$FIXSEED@Pop[])"""
        if self.ffs.flags.get('FIXSEED', True):
            rng = self.ffs.rng
            rng.seed, rng.filled, rng.next = self._seed


def _measp_word(ffs, words, i):
    """MEASURE / TRCOD の位置の語 (ielm)。要素ならその位置を measp にして 2 語、そうでなければ measp はそのまま 1 語 (語は次のコマンドとして残す)"""
    n = 1
    if i + 1 < len(words):
        try:
            k, _ = ffs._ielmf(words[i + 1])
        except Exception:
            k = None
        if k is not None:
            ffs._measp = int(k)
            n = 2
    measp = int(getattr(ffs, '_measp', 0) or 0)
    if not (1 <= measp <= ffs.nlat):
        measp = ffs.nlat                                    # tffsa.f:131 flv%measp=nlat
    return measp, n


def _precheck(ffs, cmd):
    if ffs.lat is None:
        ffs.write(f'?{cmd}: No beam line (USE first)' + chr(10))
        return False
    ffs.spread_pending()
    if mark_values(ffs) is None:
        ffs.write(' Missing entrance condition of beam.' + chr(10))          # ttinit.f
        return False
    if _np0(ffs) < 2:
        ffs.write(f'?{cmd}: NP (number of particles) must be 2 or more. Set NP = n; before {cmd}.' + chr(10))
        return False
    return True


def cmd_measure(ffs, words, i):
    """MEA_SURE [位置] (tffsa.f → trackb.f): 先頭 MARK の入口条件 (Twiss、EMITX、EMITY、DP、SIGZ、AZ) の分布の NP 個の粒子を
    位置 (既定は終端。前回の MEASURE/TRCOD の位置を覚える) の入口まで追い、重心・分散行列・射影/非結合の emittance と Twiss を表示する。
    FFS["MEASURE …"] の戻り値は {np, {重心 6}, {分散行列 6×6}}。LOGMEASURE 用に xa ya xxa xya yya を覚える。
    K64 は MEASURE で TRPT を ON にしたままにする。pysad は追跡の間だけ ON にして戻す (CALC などへの副作用を残さない)。"""
    measp, n = _measp_word(ffs, words, i)
    if not _precheck(ffs, 'MEASURE'):
        return n
    name = ffs.elname(measp)
    old_trpt = ffs.flags['TRPT']
    ffs.flags['TRPT'] = True
    tk = None
    try:
        tk = Tracker(ffs)
        tk.track(measp, cont=True)
        np_, sv = tk.stat(name, ffs.write)
    finally:
        ffs.flags['TRPT'] = old_trpt
        if tk is not None:
            tk.reset()
    if np_ <= 1:
        ffs.write(f' MEASURE: {np_} particle(s) survived at {name}; no statistics.' + chr(10))
        ffs._measure_ret = None
    else:
        ffs._measure_ret = [float(np_), [float(v) for v in tk.sa], [[float(v) for v in row] for row in tk.ss]]
    ffs._measure_sv = sv
    return n


def cmd_trcod(ffs, words, i):
    """TRC_OD [位置] (tffsa.f → trcoda.f): 分布 (JITTER は OFF 扱い = 重心を設計値に合わせる) を先頭から位置まで要素ごとに追い、
    各要素の入口での重心 (x, px, y, py) を Twiss の DX DPX DY DPY に、dp との相関 <u·dp>/<dp²> を EX EPX EY EPY に書く
    (dp の分散が 0 なら EX.. は書かない)。最後に位置での統計を表示する。位置より後ろの Twiss は変えない。
    次の CALC で Twiss は計算し直される (K64 と同じ)。
    TRCOD FIXCOD [位置] / MVCOD [位置] / STNCOD / ALIGN 要素 要素 は下の trcod_fix / trcod_move / trcod_restore / cmd_align。"""
    if i + 1 < len(words):
        sub = words[i + 1].upper()
        from .ffscmd import _abbrev
        if sub == 'FIXCOD':
            measp, n = _measp_word(ffs, words, i + 1)
            if _lat_ok(ffs, 'TRCOD FIXCOD'):
                trcod_fix(ffs, measp)
            return 1 + n
        if sub == 'MVCOD':
            k, n = _pos_word(ffs, words, i + 2)
            if _lat_ok(ffs, 'TRCOD MVCOD'):
                trcod_move(ffs, _measp_now(ffs) if k is None else k)
            return 2 + n
        if _abbrev(sub, 'STN_COD'):
            if _lat_ok(ffs, 'TRCOD STNCOD'):
                trcod_restore(ffs)
            return 2
        if sub == 'ALIGN':
            return 1 + cmd_align(ffs, words, i + 1)
    measp, n = _measp_word(ffs, words, i)
    if not _precheck(ffs, 'TRCOD'):
        return n
    name = ffs.elname(measp)
    rows = ffs.ensure_rows()
    tk = Tracker(ffs, jitter=False)
    try:
        def store(k):
            xs = tk.alive_x()
            if xs.shape[0] <= 1:
                return False
            sa, ss = ttstat(xs)
            r = rows[k - 1]
            r.x_cod, r.px_cod, r.y_cod, r.py_cod = (float(v) for v in sa[:4])
            if ss[5, 5] != 0.0:
                r.etx, r.etpx, r.ety, r.etpy = (float(ss[j, 5] / ss[5, 5]) for j in range(4))
            return True

        tk.track(1, cont=False)
        ok = store(1)
        last = 1
        for k in range(2, measp + 1):
            if not ok:
                break
            tk.track(k, cont=True)
            ok = store(k)
            last = k
        if ok:
            np_, sv = tk.stat(name, ffs.write)
        else:
            np_, sv = int((tk.alive != 0.0).sum()), [1.0e20] * 5
            ffs.write(f' TRCOD: {np_} particle(s) survived at {ffs.elname(last)}; the orbit is set only upstream of it.' + chr(10))
    finally:
        tk.reset()
    ffs._measure_sv = sv
    return n


def cmd_logmeasure(ffs, words, i):
    """LOGM_EASURE (tffsa.f 9771): 直前の MEASURE / TRCOD の xa ya xxa xya yya (重心 x, y と分散 <xx> <xy> <yy>) を 1p,5g15.7 で表示"""
    from .maintrack import fmt_g
    sv = getattr(ffs, '_measure_sv', None) or [0.0] * 5
    xa, ya, xxa, yya, xya = sv
    ffs.write(''.join(fmt_g(float(v), 15, 7) for v in (xa, ya, xxa, xya, yya)) + chr(10))
    return 1


def ffs_command(ffs, W, words, i):
    """ffscmd の FFS のコマンド分岐から (v2 モードだけ) 呼ばれる。対象でなければ None"""
    from .ffscmd import _abbrev
    if _abbrev(W, 'MEA_SURE'):
        return cmd_measure(ffs, words, i)
    if _abbrev(W, 'TRC_OD'):
        return cmd_trcod(ffs, words, i)
    if _abbrev(W, 'LOGM_EASURE'):
        return cmd_logmeasure(ffs, words, i)
    if _abbrev(W, 'ORBITJ_ITTER'):
        return cmd_orbitjitter(ffs, words, i)
    if _abbrev(W, 'PUTM_EASURE'):
        return cmd_putmeasure(ffs, words, i)
    if _abbrev(W, 'GET_MEASURE'):
        return cmd_getmeasure(ffs, words, i)
    return None


# =============================================================================================
# TRCOD FIXCOD / MVCOD / STNCOD (trcoda.f の entry trcodf / trcodm / trcods): 基準軌道の保存
# =============================================================================================
# K64 は「基準の軌道 (basic cod)」を twiss(*,ndim,*) に写し、同じものを私用の領域 itw にも取っておく。twiss(*,ndim,*) を読むのは
# 旧い補正パッケージ (mbmp.f の mccor など) と DRAW のモニターの読み (tdrwdt.f の kv=33,34) だけで、K64 ではどちらも FFS から呼ばれない
# (tffsa.f の mcmess はコメントアウト、MON コマンドは無い)。つまり K64 では保存した基準を利用者が読む手段が無い。
# v2 は基準を REFERENCE と同じ参照光学 (K64 の twiss(*,-1,*)、pysad の ffs._ref_rows) の 8 関数に書く:
# Twiss["RDX", k] で基準、Twiss["DDX", k] で 現在 − 基準 が読める (EX EPX EY EPY DX DPX DY DPY。Twiss["DX"] などは絶対値のまま)。
# 私用の写し (itw) は ffs._trcod_sav。REFERENCE が参照光学を上書きした後に TRCOD STNCOD で基準を戻せる。
_COD8 = ('etx', 'etpx', 'ety', 'etpy', 'x_cod', 'px_cod', 'y_cod', 'py_cod')
COD8_KEYS = ('EX', 'EPX', 'EY', 'EPY', 'DX', 'DPX', 'DY', 'DPY')


def _lat_ok(ffs, cmd):
    if ffs.lat is None:
        ffs.write(f'?{cmd}: No beam line (USE first)' + chr(10))
        return False
    return True


def _measp_now(ffs):
    measp = int(getattr(ffs, '_measp', 0) or 0)
    return measp if 1 <= measp <= ffs.nlat else ffs.nlat


def _pos_word(ffs, words, j):
    """words[j] が位置 (ielm) なら (位置, 1)、そうでなければ (None, 0)。measp は変えない"""
    if j < len(words):
        try:
            k, _ = ffs._ielmf(words[j])
        except Exception:
            k = None
        if k is not None:
            return int(k), 1
    return None, 0


def _ref_rows(ffs):
    """参照光学 (REFERENCE の twiss(*,-1,*)) の行。無ければ (または格子の長さが変わっていれば) 0 の行を作る"""
    from .ffscmd import _ZeroRow
    n = len(ffs.ensure_rows())
    ref = getattr(ffs, '_ref_rows', None)
    if ref is None or len(ref) != n:
        ref = [_ZeroRow() for _ in range(n)]
        ffs._ref_rows = ref
    return ref


def trcod_fix(ffs, measp):
    """TRCOD FIXCOD [位置] (trcodf): 位置 1..measp の今の Twiss の EX EPX EY EPY DX DPX DY DPY (CALC か TRCOD が入れた値) を基準にする。
    K64 は添字 7..10 と 15..18 を写す。15..18 は DETR を入れる前の番号で、今の番号では DETR DX DPX DY (DPY が入らない)。
    v2 は意図どおり DX DPX DY DPY (trcodm は mfitdx..mfitdpy の名前で書いてあり、こちらは正しい)。"""
    rows = ffs.ensure_rows()
    nlat = len(rows)
    sav = getattr(ffs, '_trcod_sav', None)
    if sav is not None and sav.shape[0] != nlat:
        ffs.write(f' size of lattice has been changed. {nlat} != {sav.shape[0]}' + chr(10))
        sav = None
    if sav is None:
        sav = np.zeros((nlat, 8))
        ffs._trcod_sav_n = 0
    ffs.write('  trcodf' + f'{measp:12d}' + ' ' + ffs.elname(measp) + chr(10))
    ref = _ref_rows(ffs)
    for k in range(measp):
        for j, a in enumerate(_COD8):
            v = float(getattr(rows[k], a, 0.0))
            setattr(ref[k], a, v)
            sav[k, j] = v
    ffs._trcod_sav = sav
    ffs._trcod_sav_n = max(int(getattr(ffs, '_trcod_sav_n', 0)), measp)


def trcod_move(ffs, pos):
    """TRCOD MVCOD [位置] (trcodm): 位置 pos (省略すると MEASURE / TRCOD の位置) の 1 点だけ、基準の DX DPX DY DPY を今の値にする
    (EX.. は変えない。私用の写しがあればそれも)。"""
    rows = ffs.ensure_rows()
    ref = _ref_rows(ffs)
    sav = getattr(ffs, '_trcod_sav', None)
    if sav is not None and sav.shape[0] != len(rows):
        sav = None
    for j in range(4, 8):
        v = float(getattr(rows[pos - 1], _COD8[j], 0.0))
        setattr(ref[pos - 1], _COD8[j], v)
        if sav is not None:
            sav[pos - 1, j] = v


def trcod_restore(ffs):
    """TRCOD STNCOD (trcods): 私用の写し (FIXCOD・MVCOD が入れた基準) を参照光学の 8 関数へ戻す (REFERENCE で上書きした後など)。
    K64 は戻す範囲を、値の定まらない局所変数 meas0 で決める。v2 は FIXCOD で保存した範囲。FIXCOD がまだなら何もしない (K64 と同じ) が、注意を出す。"""
    sav = getattr(ffs, '_trcod_sav', None)
    rows = ffs.ensure_rows()
    if sav is None or sav.shape[0] != len(rows):
        ffs.write(' TRCOD STNCOD: no reference orbit is saved for this beam line (TRCOD FIXCOD first); nothing is done.' + chr(10))
        return False
    ref = _ref_rows(ffs)
    for k in range(int(getattr(ffs, '_trcod_sav_n', 0))):
        for j, a in enumerate(_COD8):
            setattr(ref[k], a, float(sav[k, j]))
    return True


def trcod_reference(ffs):
    """Python から: TRCOD の基準 (参照光学の EX EPX EY EPY DX DPX DY DPY) を (nlat, 8) の配列で返す。列の順は COD8_KEYS"""
    ref = _ref_rows(ffs)
    return np.array([[float(getattr(r, a, 0.0)) for a in _COD8] for r in ref])


# =============================================================================================
# TRCOD ALIGN 要素1 要素2 (nalign.f)
# =============================================================================================
_ALIGN_K64 = ('Bend', 'Quadrupole', 'Sextupole', 'Octupole')            # nalign.f の id = 2, 4, 6, 8
_ALIGN_V2 = _ALIGN_K64 + ('Decapole', 'Dodecapole', 'Multipole')        # v2 の拡張 (DX・DY を持つ磁石)


def cmd_align(ffs, words, i):
    """TRCOD ALIGN 要素1 要素2 (nalign.f)。i は ALIGN の語の位置。戻りは ALIGN から数えた消費語数。
    要素1 (最初に合う位置 is) の DX・DY を、is+1 から 要素2 (最初に合う位置 ie) までの磁石の DX・DY に書く (区間の磁石を先頭の磁石の
    ずれに揃える。同じ架台に載った磁石の並び)。is > ie なら入れ替える。要素1 が磁石でなければ 0 に揃える (K64 と同じ)。
    対象は BEND QUAD SEXT OCT (K64) と、v2 の拡張で DECA DODECA MULT。分割された磁石 (同じ名前が隙間なく並ぶ) は、
    is をその先頭、ie をその末尾まで広げる。
    K64 との違い: (1) K64 は `if (master(is) .lt. 0) then is=-master(is)` が固定形式では変数 thenis への代入になり、is を先頭へ戻さない。
    (2) K64 の id=20 の枝 (今の番号では SOL。DX DY を「設計値 − ずれ」にする) は旧い型番号の名残と判断して、v2 は SOL・MONI を動かさない。
    (3) 要素が見つからないとき K64 は黙って戻る。v2 は理由を表示する。
    要素2 が見つからなければ K64 と同じく終端までを対象にする (その語は捨てる)。"""
    nlat = ffs.nlat
    if not _lat_ok(ffs, 'TRCOD ALIGN'):
        return len(words) - i
    if i + 1 >= len(words):
        ffs.write('?TRCOD ALIGN: missing elements (usage: TRCOD ALIGN first-element last-element)' + chr(10))
        return 1
    ffs.spread_pending()

    def first(w):
        try:
            p = [k for k in ffs.positions(w) if 1 <= k < nlat]
        except Exception:
            p = []
        return min(p) if p else None

    is_ = first(words[i + 1])
    if is_ is None:
        ffs.write(f'?TRCOD ALIGN: undefined element {words[i + 1]}; nothing is done.' + chr(10))
        return 2
    n = 2
    ie = nlat
    if i + 2 < len(words):
        n = 3
        k = first(words[i + 2])
        if k is None:
            ffs.write(f' TRCOD ALIGN: element {words[i + 2]} is not found; aligned to the end of the beam line.' + chr(10))
        else:
            ie = k
    else:
        ffs.write(' TRCOD ALIGN: the last element is not given; aligned to the end of the beam line.' + chr(10))
    if is_ > ie:
        is_, ie = ie, is_
    master = ffs._master_array()
    if master[is_] < 0:
        is_ = -master[is_]
    if ie < nlat and master[ie] > 0:
        ie = master[ie]
    el0 = ffs.els[is_ - 1]
    if type(el0).__name__ in _ALIGN_V2 and hasattr(el0, 'DX'):
        dx, dy = float(el0.DX or 0.0), float(getattr(el0, 'DY', 0.0) or 0.0)
    else:
        dx = dy = 0.0
    cnt = 0
    for k in range(is_ + 1, min(ie, nlat - 1) + 1):
        el = ffs.els[k - 1]
        if type(el).__name__ in _ALIGN_V2 and hasattr(el, 'DX'):
            el.DX = dx
            el.DY = dy
            cnt += 1
    ffs._align_last = (is_, ie, dx, dy, cnt)
    ffs.invalidate()
    return n


# =============================================================================================
# ORBITJITTER x y (tfojit.f)
# =============================================================================================
def cmd_orbitjitter(ffs, words, i):
    """ORBITJ_ITTER xj yj (tfojit.f): 今の Twiss の DX に xj·tgauss()、DY に yj·tgauss() を、位置 1..nlat の順に (各位置で x、y の順に乱数を
    引いて) 足す。tgauss は GCUT で切った正規乱数 (FFS の乱数列を進める)。DPX・DPY と要素は変えない。次の CALC で Twiss は計算し直される。
    測定の読みの揺らぎ (各位置で独立) を Twiss の軌道に重ねる用途。値が 2 つ無ければ K64 と同じ文を出して行の残りを捨てる。"""
    if not _lat_ok(ffs, 'ORBITJITTER'):
        return len(words) - i
    vals = []
    for j in (1, 2):
        v, n = ffs._getva(words, i + j)
        if n == 0 or isinstance(v, list):
            ffs.write(' ?Missing parameter for ORBITJ_ITTER x y' + chr(10))
            return len(words) - i
        vals.append(float(v))
    xj, yj = vals
    rows = ffs.ensure_rows()
    rng = ffs.rng
    gcut = ffs.gcut()
    for r in rows:
        r.x_cod = float(r.x_cod) + xj * rng.gauss(1, gcut)
        r.y_cod = float(r.y_cod) + yj * rng.gauss(1, gcut)
    return 3


# =============================================================================================
# PUTMEASURE / GETMEASURE (nfgetm.c の nfputm_ / nfgetm_、tfgetm.f)
# =============================================================================================
def _putm_line(v5):
    return ','.join('%15.7g' % float(v) for v in v5)                     # nfgetm.c: "%15.7lg,%15.7lg,%15.7lg,%15.7lg,%15.7lg\n"


def close_measure_files(ffs):
    """PUTMEASURE / GETMEASURE が開いたままにしているファイルを閉じる (Python から。次の PUTMEASURE / GETMEASURE は開き直す:
    GETMEASURE はファイルの先頭から読み直す)"""
    for a in ('_putm', '_getm'):
        st = getattr(ffs, a, None)
        if st is not None:
            try:
                st[1].close()
            except Exception:
                pass
        setattr(ffs, a, None)


def cmd_putmeasure(ffs, words, i):
    """PUTM_EASURE (nfputm_): 直前の MEASURE / TRCOD の xa ya xxa xya yya (重心 x、y と <xx> <xy> <yy>) を、環境変数 SADPOUT の名前の
    ファイル (Unix なら FIFO でもよい) に 1 行 `%15.7g` 5 個のコンマ区切りで追記する。ファイルは最初の PUTMEASURE で追記用に開き、開いたまま
    にする (行ごとに flush)。標準エラーに `Putting data:…` を出す (K64 と同じ)。
    K64 との違い: SADPOUT が無い・開けないとき K64 は落ちる (exit(1)、または NULL の参照)。v2 は理由を表示して続ける。
    SADPOUT の値が変わっていたら開き直す。"""
    sv = getattr(ffs, '_measure_sv', None) or [0.0] * 5
    xa, ya, xxa, yya, xya = sv
    path = os.environ.get('SADPOUT', '')
    if not path:
        ffs.write('?PUTMEASURE: environment variable SADPOUT (name of the file to append the measured values) is not set.' + chr(10))
        return 1
    st = getattr(ffs, '_putm', None)
    if st is None or st[0] != path or st[1].closed:
        if st is not None and not st[1].closed:
            st[1].close()
        sys.stderr.write(f' opening pipe:{path} for write' + chr(10))
        try:
            fh = open(path, 'a', encoding='ascii', newline=chr(10))
        except OSError as e:
            ffs._putm = None
            ffs.write(f'?PUTMEASURE: unable to open {path} for write ({e.strerror})' + chr(10))
            return 1
        sys.stderr.write(f' {path} is opened for write' + chr(10))
        st = ffs._putm = (path, fh)
    line = _putm_line((xa, ya, xxa, xya, yya))
    st[1].write(line + chr(10))
    st[1].flush()
    sys.stderr.write('Putting data:' + line + ' ' + chr(10))
    sys.stderr.flush()
    return 1


def read_measure_record(ffs):
    """nfgetm_: 環境変数 SADPIN のファイルから次の 1 組 (xa ya xxa xya yya: コンマ区切りの 5 個) を読む。戻り 5 個のリスト、無ければ None
    (理由を表示)。ファイルは最初の呼び出しで開き、開いたままにする (呼ぶたびに次の行)。"""
    path = os.environ.get('SADPIN', '')
    if not path:
        ffs.write('?GETMEASURE: environment variable SADPIN (name of the file to read the measured values from) is not set.' + chr(10))
        return None
    st = getattr(ffs, '_getm', None)
    if st is None or st[0] != path or st[1].closed:
        if st is not None and not st[1].closed:
            st[1].close()
        try:
            fh = open(path, 'r', encoding='ascii', errors='replace')
        except OSError as e:
            ffs._getm = None
            ffs.write(f'?GETMEASURE: unable to open {path} for read ({e.strerror})' + chr(10))
            return None
        st = ffs._getm = (path, fh)
    while True:
        line = st[1].readline()
        if line == '':
            ffs.write(f'?GETMEASURE: no more data in {path}' + chr(10))
            return None
        if line.strip():
            break
    try:
        v = [float(t) for t in line.replace('D', 'E').replace('d', 'e').split(',')]
    except ValueError:
        v = []
    if len(v) != 5:
        ffs.write(f'?GETMEASURE: a record of 5 comma-separated numbers is expected in {path}: {line.strip()}' + chr(10))
        return None
    sys.stderr.write('data gotten from :' + _putm_line(v) + ' ' + chr(10))
    sys.stderr.flush()
    return v


def cmd_getmeasure(ffs, words, i):
    """GET_MEASURE [alpha] (tfgetm.f + nfgetm_): SADPIN のファイルから相手の測定値 xa1 ya1 xxa1 xya1 yya1 を 1 組読み、
    MEASURE / TRCOD の位置 measp の 1 点の fit 条件 (FIT measp; DX …; DY …) の目標値を
        DX の目標 = (xa1 − xa)/2·alpha + 今の目標·(1 − alpha)、DY の目標 = (ya1 − ya)/2·alpha + 今の目標·(1 − alpha)
    に書き換える (xa、ya は自分の直前の MEASURE / TRCOD の重心。alpha の既定は 1)。二つのビームをそれぞれ別の SAD で追い、互いの重心の差の半分ずつ
    軌道を寄せる (衝突点の軌道フィードバック) ための道具。条件が無ければ値を読むだけ。
    K64 との違い: (1) K64 は関数の番号 32 と 34 を見る。DX・DY の番号が今と違った頃の番号で、今の番号では PEPY と PZPX に当たり、重心の差を
    書く相手として意味を成さない。v2 は DX と DY (式が x と y の重心の差なので、軌道の関数と判断)。**ソースの読みと式からの判断で、gs での確認は無い。**
    (2) v2 は FFS["GETMEASURE"] の戻り値に読んだ 5 個 {xa1, ya1, xxa1, xya1, yya1} を返す (K64 は Null)。読めなければ目標値は変えない。"""
    v, n = ffs._getva(words, i + 1)
    if n == 0 or isinstance(v, list):
        alphad, n = 1.0, 0
    else:
        alphad = float(v)
    rec = read_measure_record(ffs)
    ffs._getm_last = rec
    if rec is None:
        return 1 + n
    xa1, ya1 = rec[0], rec[1]
    sv = getattr(ffs, '_measure_sv', None) or [0.0] * 5
    xa, ya = float(sv[0]), float(sv[1])
    ffs._measure_ret = [float(x) for x in rec]
    if ffs.lat is None:
        return 1 + n
    measp = _measp_now(ffs)
    m = ffs.matcher_()
    for c in m.conditions:
        if c['kp'] == measp and c['kp1'] == measp and c['mfitp'] != 0 and c['kf'] in ('DX', 'DY'):
            d = (xa1 - xa) if c['kf'] == 'DX' else (ya1 - ya)
            single = c['fmin'] == c['fmax']
            c['fmin'] = d * 0.5 * alphad + float(c['fmin']) * (1.0 - alphad)
            if single:
                c['fmax'] = c['fmin']
    return 1 + n


# =============================================================================================
# MAIN の TRACK … TURNS=-1 (track.f: trpt=.true., rfsw=.true., nturn=1 → tracka.f の trpt の枝)
# =============================================================================================
def track_trpt(ffs, spec, write):
    """MAIN の `TRACK USE=line NPART=n TURNS=-1;` (輸送路の 1 回通過)。maintrack.run_track から v2 モードだけ呼ばれる。
    tracka.f: ttinit (先頭 MARK の Twiss と EMITX EMITY DP SIGMAZ AZ の分布) → `Statistics at Entrance` → tturn(1, nlat) を 1 回
    → `Statistics at Exit`。乱数の種は戻さない (tracka は NISTACK$FIXSEED を使わない)。
    粒子種の語 (PROTON など) は見ない (質量・電荷は大域変数 MASS・CHARGE)。
    戻り {'alive', 'final' (np0, 6: 落ちた粒子は落ちた時点の値), 'entrance': (sa, ss), 'exit': (sa, ss) か None, 'x0'}"""
    import time
    from .maintrack import fmt_f, echo_lists
    from .sadscript import fmt_real
    t_start = time.time()
    np0 = int(spec.get('npart') or 0)
    mk = mark_values(ffs, main=True)
    if mk is None:
        write(' Missing entrance condition of beam.' + chr(10))              # ttinit.f (nturn=0 にして戻る)
        return None
    if spec.get('ptype') not in (None, '', 'ELEC', 'ELECTRON', 'POSI', 'POSITRON'):
        write(f' TRACK TURNS=-1: particle keyword {spec.get("ptype")} is ignored (MASS and CHARGE are used).' + chr(10))
    fl = ffs.flags
    old = (fl['TRPT'], fl['RFSW'])
    fl['TRPT'], fl['RFSW'] = True, True
    try:
        x0 = ttinit(mk, np0, ffs.rng, ffs.gcut(), gauss=bool(fl['GAUSS']), trgauss=bool(fl['TRGAUSS']),
                    jitter=bool(fl['JITTER']), dp0=ffs.dp0(), wake=False)
        info = {'RAD': fl['RAD'], 'RFSW': True, 'GAUSS': fl['GAUSS'], 'DP': mk['DP'], 'DP0': ffs.dp0(), 'GCUT': ffs.gcut()}
        ent = ex = None
        if np0 >= 2:
            lines, sa, ss, _ = ttstat_lines(x0, 'Entrance', info)
            ent = (sa, ss)
            for ln in lines:
                write(ln + chr(10))
        cols = [[float(v) for v in x0[:, j]] for j in range(6)] + [[1.0] * np0]
        r = ffs.fn_track([[1.0, cols], float(ffs.nlat)])[1]
        xf = np.array(r[:6], dtype=float).T
        alive = np.array(r[6], dtype=float) != 0.0
        if int(alive.sum()) >= 2:
            lines, sa, ss, _ = ttstat_lines(xf[alive], 'Exit', info)
            ex = (sa, ss)
            for ln in lines:
                write(ln + chr(10))
        else:
            write(f' TRACK: {int(alive.sum())} particle(s) survived; no statistics at Exit.' + chr(10))
    finally:
        fl['TRPT'], fl['RFSW'] = old
    write(' ' + 'Tracking'.ljust(20) + ' end:  CPU time =' + fmt_f(time.time() - t_start, 10, 3) + ' sec' + chr(10))
    for ln in echo_lists(spec.get('nx'), spec.get('ny'), spec.get('nz'), fmt_real):
        write(ln + chr(10))
    return {'alive': int(alive.sum()), 'final': xf, 'entrance': ent, 'exit': ex, 'x0': x0}
