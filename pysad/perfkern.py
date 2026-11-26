"""Numba performance kernels for lattice.py, giving bitwise the same results as the pure Python versions.
perfkern.py — lattice.py の性能用 numba カーネル (結果は pure Python 版と同じ演算順で、ビット単位で一致させる)

- twiss_chain: lattice._transport_twiss (SAD qmat2twiss の直訳) を全要素に順に適用する。
  Python 版 (dict を要素ごとに作る) は 1 要素 28 µs、これは 3563 要素で 1 ms 以下。
"""
import math

import numpy as np
from numba import njit

NJ = dict(cache=True, error_model='numpy', nogil=True)
_XYTH = 0.375

# 状態ベクトルの並び: ax bx ay by nx ny ex epx ey epy r1 r2 r3 r4 detr
TW_KEYS = ('ax', 'bx', 'ay', 'by', 'nx', 'ny', 'ex', 'epx', 'ey', 'epy', 'r1', 'r2', 'r3', 'r4', 'detr')


@njit(**NJ)
def _step(st, M, clamp, out, v2):
    ax0 = st[0]
    bx0 = st[1]
    ay0 = st[2]
    by0 = st[3]
    nx0 = st[4]
    ny0 = st[5]
    ex0 = st[6]
    epx0 = st[7]
    ey0 = st[8]
    epy0 = st[9]
    r1 = st[10]
    r2 = st[11]
    r3 = st[12]
    r4 = st[13]
    detr0 = st[14]
    normal = detr0 < 1.0

    # t = [M[:4,:4] | M[:4,5]]
    a11 = M[0, 0]; a12 = M[0, 1]; a13 = M[0, 2]; a14 = M[0, 3]; a15 = M[0, 5]
    a21 = M[1, 0]; a22 = M[1, 1]; a23 = M[1, 2]; a24 = M[1, 3]; a25 = M[1, 5]
    a31 = M[2, 0]; a32 = M[2, 1]; a33 = M[2, 2]; a34 = M[2, 3]; a35 = M[2, 5]
    a41 = M[3, 0]; a42 = M[3, 1]; a43 = M[3, 2]; a44 = M[3, 3]; a45 = M[3, 5]
    coup = (a13 != 0.0 or a14 != 0.0 or a23 != 0.0 or a24 != 0.0
            or a31 != 0.0 or a32 != 0.0 or a41 != 0.0 or a42 != 0.0)

    if coup:
        sqa = 1.0 - r1 * r4 + r2 * r3
        sq = math.sqrt(sqa if sqa > 0.0 else 0.0)        # Python の max(0.0, x) と同じ (NaN は 0.0)
        if normal:
            r15 = sq * ex0 + r4 * ey0 - r2 * epy0
            r25 = sq * epx0 - r3 * ey0 + r1 * epy0
            r35 = sq * ey0 - r1 * ex0 - r2 * epx0
            r45 = sq * epy0 - r3 * ex0 - r4 * epx0
            # 行 i の (a1..a5)
            u11 = a11 * sq - a13 * r1 - a14 * r3
            u12 = a12 * sq - a13 * r2 - a14 * r4
            u13 = a11 * r4 - a12 * r3 + a13 * sq
            u14 = -a11 * r2 + a12 * r1 + a14 * sq
            u15 = a11 * r15 + a12 * r25 + a13 * r35 + a14 * r45 + a15
            u21 = a21 * sq - a23 * r1 - a24 * r3
            u22 = a22 * sq - a23 * r2 - a24 * r4
            u23 = a21 * r4 - a22 * r3 + a23 * sq
            u24 = -a21 * r2 + a22 * r1 + a24 * sq
            u25 = a21 * r15 + a22 * r25 + a23 * r35 + a24 * r45 + a25
            u31 = a31 * sq - a33 * r1 - a34 * r3
            u32 = a32 * sq - a33 * r2 - a34 * r4
            u33 = a31 * r4 - a32 * r3 + a33 * sq
            u34 = -a31 * r2 + a32 * r1 + a34 * sq
            u35 = a31 * r15 + a32 * r25 + a33 * r35 + a34 * r45 + a35
            u41 = a41 * sq - a43 * r1 - a44 * r3
            u42 = a42 * sq - a43 * r2 - a44 * r4
            u43 = a41 * r4 - a42 * r3 + a43 * sq
            u44 = -a41 * r2 + a42 * r1 + a44 * sq
            u45 = a41 * r15 + a42 * r25 + a43 * r35 + a44 * r45 + a45
        else:
            r15 = -r1 * ex0 - r2 * epx0 + sq * ey0
            r25 = -r3 * ex0 - r4 * epx0 + sq * epy0
            r35 = sq * ex0 + r4 * ey0 - r2 * epy0
            r45 = sq * epx0 - r3 * ey0 + r1 * epy0
            u11 = -a11 * r1 - a12 * r3 + a13 * sq
            u12 = -a11 * r2 - a12 * r4 + a14 * sq
            u13 = a11 * sq + a13 * r4 - a14 * r3
            u14 = a12 * sq - a13 * r2 + a14 * r1
            u15 = a11 * r15 + a12 * r25 + a13 * r35 + a14 * r45 + a15
            u21 = -a21 * r1 - a22 * r3 + a23 * sq
            u22 = -a21 * r2 - a22 * r4 + a24 * sq
            u23 = a21 * sq + a23 * r4 - a24 * r3
            u24 = a22 * sq - a23 * r2 + a24 * r1
            u25 = a21 * r15 + a22 * r25 + a23 * r35 + a24 * r45 + a25
            u31 = -a31 * r1 - a32 * r3 + a33 * sq
            u32 = -a31 * r2 - a32 * r4 + a34 * sq
            u33 = a31 * sq + a33 * r4 - a34 * r3
            u34 = a32 * sq - a33 * r2 + a34 * r1
            u35 = a31 * r15 + a32 * r25 + a33 * r35 + a34 * r45 + a35
            u41 = -a41 * r1 - a42 * r3 + a43 * sq
            u42 = -a41 * r2 - a42 * r4 + a44 * sq
            u43 = a41 * sq + a43 * r4 - a44 * r3
            u44 = a42 * sq - a43 * r2 + a44 * r1
            u45 = a41 * r15 + a42 * r25 + a43 * r35 + a44 * r45 + a45
        detp = (u11 * u22 - u21 * u12 + u33 * u44 - u43 * u34) * 0.5
        if detp > _XYTH:
            sd = math.sqrt(detp)
            u11 /= sd; u12 /= sd; u21 /= sd; u22 /= sd
            u33 /= sd; u34 /= sd; u43 /= sd; u44 /= sd
            r1n = -u31 * u22 + u32 * u21
            r2n = u31 * u12 - u32 * u11
            r3n = -u41 * u22 + u42 * u21
            r4n = u41 * u12 - u42 * u11
            detr = r1n * r4n - r2n * r3n
            sd = math.sqrt(1.0 - detr)
            ex1 = u15 * sd - r4n * u35 + r2n * u45
            epx1 = u25 * sd + r3n * u35 - r1n * u45
            ey1 = u35 * sd + r1n * u15 + r2n * u25
            epy1 = u45 * sd + r3n * u15 + r4n * u25
            aa = bx0 * u11 - ax0 * u12
            cc = (u12 - ax0 * aa) / bx0
            ax1 = -(u21 * aa + u22 * cc)
            bx1 = u11 * aa + u12 * cc
            dpsix = math.atan2(u12, aa)
            bb = by0 * u33 - ay0 * u34
            dd = (u34 - ay0 * bb) / by0
            ay1 = -(u43 * bb + u44 * dd)
            by1 = u33 * bb + u34 * dd
            dpsiy = math.atan2(u34, bb)
            if v2 and not normal:                      # v2: 入れ替わりから出る跳びは −pi (lattice._exch_phase)
                if aa < 0.0 and abs(u12) < 1e-9 * abs(aa):
                    dpsix = -math.pi
                if bb < 0.0 and abs(u34) < 1e-9 * abs(bb):
                    dpsiy = -math.pi
        else:
            detm = (u31 * u42 - u32 * u41 + u13 * u24 - u14 * u23) * 0.5
            sd = math.sqrt(abs(detm))
            u31 /= sd; u32 /= sd; u41 /= sd; u42 /= sd
            u13 /= sd; u14 /= sd; u23 /= sd; u24 /= sd
            r1n = -u11 * u42 + u12 * u41
            r2n = u11 * u32 - u12 * u31
            r3n = -u21 * u42 + u22 * u41
            r4n = u21 * u32 - u22 * u31
            sd = math.sqrt(1.0 - r1n * r4n + r2n * r3n)
            ex1 = u35 * sd - r4n * u15 + r2n * u25
            epx1 = u45 * sd + r3n * u15 - r1n * u25
            ey1 = u15 * sd + r1n * u35 + r2n * u45
            epy1 = u25 * sd + r3n * u35 + r4n * u45
            detr = 1.0 + _XYTH - r1n * r4n + r2n * r3n
            aa = bx0 * u31 - ax0 * u32
            cc = (u32 - ax0 * aa) / bx0
            ax1 = -(u41 * aa + u42 * cc)
            bx1 = u31 * aa + u32 * cc
            dpsix = math.atan2(u32, aa)
            bb = by0 * u13 - ay0 * u14
            dd = (u14 - ay0 * bb) / by0
            ay1 = -(u23 * bb + u24 * dd)
            by1 = u13 * bb + u14 * dd
            dpsiy = math.atan2(u14, bb)
            if v2 and normal:                          # v2: 入れ替わりに入る跳びは +pi
                if aa < 0.0 and abs(u32) < 1e-9 * abs(aa):
                    dpsix = math.pi
                if bb < 0.0 and abs(u14) < 1e-9 * abs(bb):
                    dpsiy = math.pi
    else:
        if normal:
            t1 = r1 * a22 - r2 * a21
            t2 = -r1 * a12 + r2 * a11
            t3 = r3 * a22 - r4 * a21
            t4 = -r3 * a12 + r4 * a11
            r1n = a33 * t1 + a34 * t3
            r2n = a33 * t2 + a34 * t4
            r3n = a43 * t1 + a44 * t3
            r4n = a43 * t2 + a44 * t4
            detr = r1n * r4n - r2n * r3n
            sd = math.sqrt(1.0 - detr)
            ex1 = a11 * ex0 + a12 * epx0 + sd * a15 - r4n * a35 + r2n * a45
            epx1 = a21 * ex0 + a22 * epx0 + sd * a25 + r3n * a35 - r1n * a45
            ey1 = a33 * ey0 + a34 * epy0 + sd * a35 + r1n * a15 + r2n * a25
            epy1 = a43 * ey0 + a44 * epy0 + sd * a45 + r3n * a15 + r4n * a25
            aa = bx0 * a11 - ax0 * a12
            cc = (a12 - ax0 * aa) / bx0
            ax1 = -(a21 * aa + a22 * cc)
            bx1 = a11 * aa + a12 * cc
            dpsix = math.atan2(a12, aa)
            bb = by0 * a33 - ay0 * a34
            dd = (a34 - ay0 * bb) / by0
            ay1 = -(a43 * bb + a44 * dd)
            by1 = a33 * bb + a34 * dd
            dpsiy = math.atan2(a34, bb)
        else:
            sd = math.sqrt(1.0 - r1 * r4 + r2 * r3)
            t1 = r1 * a44 - r2 * a43
            t2 = -r1 * a34 + r2 * a33
            t3 = r3 * a44 - r4 * a43
            t4 = -r3 * a34 + r4 * a33
            r1n = a11 * t1 + a12 * t3
            r2n = a11 * t2 + a12 * t4
            r3n = a21 * t1 + a22 * t3
            r4n = a21 * t2 + a22 * t4
            detr = 1.0 + _XYTH - r1n * r4n + r2n * r3n
            ex1 = a33 * ex0 + a34 * epx0 + sd * a35 + (-r4n * a15 + r2n * a25)
            epx1 = a43 * ex0 + a44 * epx0 + sd * a45 + (r3n * a15 - r1n * a25)
            ey1 = a11 * ey0 + a12 * epy0 + sd * a15 + (r1n * a35 + r2n * a45)
            epy1 = a21 * ey0 + a22 * epy0 + sd * a25 + (r3n * a35 + r4n * a45)
            aa = bx0 * a33 - ax0 * a34
            cc = (a34 - ax0 * aa) / bx0
            ax1 = -(a43 * aa + a44 * cc)
            bx1 = a33 * aa + a34 * cc
            dpsix = math.atan2(a34, aa)
            bb = by0 * a11 - ay0 * a12
            dd = (a12 - ay0 * bb) / by0
            ay1 = -(a21 * bb + a22 * dd)
            by1 = a11 * bb + a12 * dd
            dpsiy = math.atan2(a12, bb)

    if clamp:
        bx1 = 1e-12 if 1e-12 > bx1 else bx1          # Python の max(bx1, 1e-12) と同じ (NaN は残る)
        by1 = 1e-12 if 1e-12 > by1 else by1
    out[0] = ax1
    out[1] = bx1
    out[2] = ay1
    out[3] = by1
    out[4] = nx0 + dpsix
    out[5] = ny0 + dpsiy
    out[6] = ex1
    out[7] = epx1
    out[8] = ey1
    out[9] = epy1
    out[10] = r1n
    out[11] = r2n
    out[12] = r3n
    out[13] = r4n
    out[14] = detr


@njit(**NJ)
def twiss_chain(tw0, mats, clamp, check, v2=False):
    """Propagate the initial Twiss tw0 through the matrices mats in sequence and return (out, nbad), stopping at the first bad element.
    tw0 (15,) を mats (n,6,6) で順に伝搬。戻り (out (n,15), nbad)。
    check=True (SAD 流の続行): 各要素の入口で βx>0 かつ βy>0 でなければ、または出力が非有限なら、その要素 i で打ち切って nbad=i
    (out[:i] が有効)。check=False: 出力に非有限が出たらその要素 i を nbad として最後まで続ける (呼び出し側が pure Python に戻す)。
    nbad=-1 は異常なし"""
    n = mats.shape[0]
    out = np.empty((n, 15))
    st = tw0.copy()
    nbad = -1
    for i in range(n):
        if check and not (st[1] > 0.0 and st[3] > 0.0):
            return out, i
        _step(st, mats[i], clamp, out[i], v2)
        ok = True
        for j in range(15):
            if not math.isfinite(out[i, j]):
                ok = False
        if not ok:
            if check:
                return out, i
            if nbad < 0:
                nbad = i
        for j in range(15):
            st[j] = out[i, j]
    return out, nbad


# ---------------------------------------------------------------------------
# emit6d._State.tradke の D・B 更新部 (want_db=True のとき)。numpy の小行列演算の連続 (1 回 0.5 ms) を 1 回 20 µs 程度に
# ---------------------------------------------------------------------------

@njit(**NJ)
def _mm(a, b):
    c = np.zeros((6, 6))
    for i in range(6):
        for k in range(6):
            aik = a[i, k]
            if aik != 0.0:
                for j in range(6):
                    c[i, j] += aik * b[k, j]
    return c


@njit(**NJ)
def edge_len(tx, ty, x, y, px, py, dl, sgn):
    """Compute only the path length l of edge_free (no derivatives; used for tracking zr0).
    edge_free の行程 ℓ だけ (微分なし。追跡の zr0 用)"""
    pr = 1.0 + dl
    u = tx * x + ty * y
    pz2 = pr * pr - px * px - py * py
    if pz2 <= 1e-12:
        return u
    pz = math.sqrt(pz2)
    return u * (pr / pz) / (1.0 + sgn * (tx * px + ty * py) / pz)


@njit(**NJ)
def edge_free(tx, ty, x, y, px, py, dl, sgn):
    """Compute the free-flight path l and its coordinate derivatives g before and after a hard-edge face tilted by angle e (v2 B220).
    v2 (B220): 硬い端の面 (面は軸に対して角 e 傾いた平面、tx + ty の方向の tan e) の手前・先の自由な飛行の行程 ℓ とその座標微分 g (x, px, y, py, z, δ)。
    入口 (sgn = -1): 面までの直進 z_f = X tan e = (x + px/pz z_f) tan e → ℓ = z_f p/pz = u/(1 - v/pz)·p/pz。
    出口 (sgn = +1): 面から z = L の面までの直進。x(L) = X (1 + px/pz tan e) → ℓ = u/(1 + v/pz)·p/pz。u = tx x + ty y、v = tx px + ty py。
    tan e → 0 の極限が従来の線形の式 (ℓ = tan e·x)。戻り (ℓ, g)。"""
    g = np.zeros(6)
    pr = 1.0 + dl
    u = tx * x + ty * y
    pz2 = pr * pr - px * px - py * py
    if pz2 <= 1e-12:
        g[0] = tx
        g[2] = ty
        return u, g
    pz = math.sqrt(pz2)
    v = tx * px + ty * py
    f = 1.0 + sgn * v / pz
    a = pr / pz
    ell = u * a / f
    ipz3 = 1.0 / (pz2 * pz)
    dU = (tx, 0.0, ty, 0.0, 0.0, 0.0)
    dV = (0.0, tx, 0.0, ty, 0.0, 0.0)
    dip = (0.0, px * ipz3, 0.0, py * ipz3, 0.0, -pr * ipz3)
    da = (0.0, pr * px * ipz3, 0.0, pr * py * ipz3, 0.0, -(px * px + py * py) * ipz3)
    for k in range(6):
        df = sgn * (dV[k] / pz + v * dip[k])
        g[k] = (dU[k] * a + u * da[k]) / f - u * a * df / (f * f)
    return ell, g


@njit(**NJ)
def tradke_db(T, transr, D, B, pr, px, py, pxi, pyi, pz, pz00, pz0, pxr0, cp, sp, xpx, xpy, xpz, xpa, anp, uc, uct,
              al1, dg, ddpx, ddpy, c1, h1, p, p0, anrad, cave, cuu, bzh=0.0, bzhr0=0.0, tin=0.0, tout=0.0, fixb39=False, rsgn=1.0, exactd=False,
              tiny=0.0, touty=0.0, dla=0.0, dlb=0.0, gin=None, gout=None):
    """Compute radiation-derivative contributions using edge-path and beta < 1 coefficient derivatives; return (changed, D, B).
    gin, gout: v2 (B220) 端の面の手前・先の自由な行程 ℓ (perfkern.edge_free) の座標微分 (始点の座標・終点の座標)。与えると tin・tout の線形の式に代わる。
    dla, dlb: v2 (B167) の β < 1 の係数 a (光子数)・b (光子エネルギー) の d ln/dγ (anp・uct は a・b を含み、anrad は abs(anrad)·a)。
    戻り (changed, D, B)。changed=False (xpa == 0) のとき D・B は変えない (呼び出し側が元のまま使う)"""
    exactd = exactd and bzh == 0.0 and bzhr0 == 0.0      # v2 の厳密な微分は SOL の外だけ (SOL 領域内の tradke は未確認)
    # transi = T @ tinv6(transr)、tinv6(r) = -J r^T J
    tinv = np.empty((6, 6))
    for i in range(6):
        si = 1.0 if i % 2 == 0 else -1.0
        for j in range(6):
            l = j ^ 1
            sl = 1.0 if l % 2 == 0 else -1.0
            tinv[i, j] = -(si * sl * transr[l, i ^ 1])
    if exactd:
        # v2 (B224): 区間の始点の行列の逆行列は、写像がシンプレクティックでなくても厳密な逆行列を使う (硬い端の MULT の端の写像は K1 ≠ 0 で
        # シンプレクティックから 5e-4〜3e-2 ずれ、シンプレクティックの逆行列の式だと微分が数値微分と 2〜3 倍違った)
        tinv = np.linalg.inv(transr)
    tr2 = _mm(T, tinv)
    if xpa == 0.0:
        return False, D, B
    if bzh != 0.0:                       # photons.f: 正準運動量 → 力学的運動量 (ソレノイド内)
        for j in range(6):
            t2 = tr2[1, j] + bzh * tr2[2, j]
            t4 = tr2[3, j] - bzh * tr2[0, j]
            tr2[1, j] = t2
            tr2[3, j] = t4
    dpxi = np.zeros(6)
    dpxi[1] = 1.0
    dpxi[2] = bzhr0
    dpyi = np.zeros(6)
    dpyi[3] = 1.0
    dpyi[0] = -bzhr0
    ddpz = np.empty(6)
    dpz00 = np.empty(6)
    dpxr0 = np.empty(6)
    dpz0 = np.empty(6)
    for j in range(6):
        ddpz[j] = (tr2[5, j] * pr - tr2[1, j] * px - tr2[3, j] * py) / pz
        dpz00[j] = (-pxi * dpxi[j] - pyi * dpyi[j]) / pz00
    dpz00[5] += pr / pz00
    for j in range(6):
        dpxr0[j] = cp * dpxi[j] + sp * dpz00[j]
        dpz0[j] = (-pxr0 * dpxr0[j] - pyi * dpyi[j]) / pz0
    dpz0[5] += pr / pz0
    dh1r = p * p0 / h1 ** 2
    dxpa = np.empty(6)
    for j in range(6):
        dxpx = tr2[3, j] * pz0 + py * dpz0[j] - ddpz[j] * pyi - pz * dpyi[j]
        dxpy = ddpz[j] * pxr0 + pz * dpxr0[j] - tr2[1, j] * pz0 - px * dpz0[j]
        dxpz = tr2[1, j] * pyi + px * dpyi[j] - tr2[3, j] * pxr0 - py * dpxr0[j]
        if exactd:
            dxpa[j] = (xpx * dxpx + xpy * dxpy + xpz * dxpz) / xpa / pr ** 4      # v2: xpa = |xp|/pr² の微分は (xp·dxp)/(xpa pr⁴) (K64 は pr²)
        else:
            dxpa[j] = (xpx * dxpx + xpy * dxpy + xpz * dxpz) / xpa / pr ** 2
    dxpa[5] -= 2.0 * xpa / pr
    dal = np.empty(6)
    for j in range(6):
        dal[j] = -tr2[4, j]
    dal[4] += 1.0
    if gin is not None:
        for j in range(6):
            dal[j] -= gin[j]
            for i in range(6):
                dal[j] -= gout[i] * tr2[i, j]
    elif tin != 0.0 or tout != 0.0 or tiny != 0.0 or touty != 0.0:
        # v2: 偏向磁石の端の角による場の中の行路の変化 (al1 = al - dz - tin*x0 - tiny*y0 - tout*x - touty*y) の微分。
        # tiny・touty は回した座標系の MULT (v2next P3: x' = cos x − sin y)
        dal[0] -= tin
        dal[2] -= tiny
        for j in range(6):
            dal[j] -= tout * tr2[0, j] + touty * tr2[2, j]
    das = 1.0 / math.sqrt(1.0 - xpa ** 2)
    dtheta = dxpa * das
    danp = anrad * h1 * dtheta
    danp[5] += anp * (dh1r + dla * p * p0 / h1)
    duc = uct * dtheta - uc * dal / al1
    duc[5] += uc * (3.0 * dh1r + dlb * p * p0 / h1)
    ddg = -cave * (danp * uc + anp * duc)
    dddpx = np.empty(6)
    dddpy = np.empty(6)
    for j in range(6):
        if exactd:
            # v2: ddpx = rsgn·0.5·dpx·dg の微分は 0.5·rsgn·(d(dpx)·dg + dpx·d(dg))。K64 は 2 項目の dpx を ddpx (= 0.5·dpx·dg) に取り違えていて、dg 倍だけ小さい
            dddpx[j] = 0.5 * rsgn * ((tr2[1, j] - dpxr0[j]) * dg + (px - pxr0) * ddg[j])
            dddpy[j] = 0.5 * rsgn * ((tr2[3, j] - dpyi[j]) * dg + (py - pyi) * ddg[j])
        else:
            dddpx[j] = 0.5 * (rsgn * (tr2[1, j] - dpxr0[j]) * dg + ddpx * ddg[j])        # rsgn: v2 は ddpx の符号を反転 (ddpx は符号つきで渡る)
            dddpy[j] = 0.5 * (rsgn * (tr2[3, j] - dpyi[j]) * dg + ddpy * ddg[j])
    tr1 = np.zeros((6, 6))
    for j in range(6):
        tr1[0, j] = c1 * dddpx[j]
        tr1[2, j] = c1 * dddpy[j]
        if exactd:
            # v2: x' = x + ddpx·c1 の c1 = al1/(3 pr) の微分 ddpx·dc1 (dc1 = dal/(3 pr) - c1 dpr/pr)。K64 はこの項を落としている
            tr1[0, j] += ddpx * dal[j] / (3.0 * pr)
            tr1[2, j] += ddpy * dal[j] / (3.0 * pr)
        tr1[1, j] = (tr2[1, j] * dg + px * ddg[j]) / pr + dddpx[j]
        tr1[3, j] = (tr2[3, j] * dg + py * ddg[j]) / pr + dddpy[j]
        tr1[5, j] = ddg[j]
    if exactd:
        tr1[0, 5] -= ddpx * c1 / pr         # v2: ∂c1/∂δ = -c1/pr (K64 は c1 を掛けずに -ddpx/pr で、長さでなく角度の次元になっている)
        tr1[2, 5] -= ddpy * c1 / pr
        tr1[1, 5] -= px * dg / pr ** 2
        tr1[3, 5] -= py * dg / pr ** 2
    else:
        tr1[0, 5] -= ddpx / pr
        tr1[2, 5] -= ddpy / pr
        tr1[1, 5] -= px * dg / pr
        tr1[3, 5] -= py * dg / pr
    if exactd:
        # v2: 上の tr1 は区間の始点の座標 j に関する微分 (tr2 = 区間の写像の微分を使って求めた)。D・B の更新 (Dn = D + tr1 (T + D)、
        # Bn = (I + tr1) B (I + tr1)^T) は区間の終点の座標に関する微分を要る (T + D は終点までの行列)。K64 は始点に関する tr1 をそのまま使う
        # (スピンの更新 srot は始点の行列 transr を掛けていて正しい)。tr1 ← tr1 tr2^-1 (tr2 はシンプレクティックなので逆行列は -J tr2^T J)
        tr2i = np.empty((6, 6))
        for i in range(6):
            si = 1.0 if i % 2 == 0 else -1.0
            for j in range(6):
                l = j ^ 1
                sl = 1.0 if l % 2 == 0 else -1.0
                tr2i[i, j] = -(si * sl * tr2[l, i ^ 1])
        tr2i = np.linalg.inv(tr2)       # v2 (B224): 厳密な逆行列 (上の式は tr2 がシンプレクティックのときだけ正しい)
        tr1 = _mm(tr1, tr2i)
    if bzh != 0.0:
        for j in range(6):
            t2 = tr1[1, j] - bzh * tr1[2, j]
            t4 = tr1[3, j] + bzh * tr1[0, j]
            tr1[1, j] = t2
            tr1[3, j] = t4
    TD = T + D
    Dn = D + _mm(tr1, TD)
    Jm = tr1.copy()
    for i in range(6):
        Jm[i, i] += 1.0
    JB = _mm(Jm, B)
    Bn = _mm(JB, Jm.T.copy())
    dee = anp * uc ** 2 * cuu
    pxm = pxi + px
    pym = pyi + py
    if fixb39:
        b = bzh     # v2 (B106): 力学的運動量は px + bzh y なので ビームの広がりの項の b も bzh (K64: bzh*0.5)
    else:
        b = bzh * 0.5
    # K64 photons.f の dbeam(8) は beam(2)-beam(10)+b*beam(4) (= B12 - B44 + b B13): そのまま写す (B34 の誤記に見える。bzh=0 では効かない)
    d11 = (Bn[1, 1] + b * (2.0 * Bn[1, 2] + b * Bn[2, 2]) + (pxm ** 2 + pxi ** 2 + px ** 2) / 6.0) * dee
    if fixb39:
        # v2 (B39): <(px+b y)(py-b x)> = B[py,px] - b (B[x,px] - B[y,py] + b B[x,y])
        d31 = (Bn[3, 1] - b * (Bn[0, 1] - Bn[2, 3] + b * Bn[0, 2]) + (pxm * pym + pxi * pyi + px * py) / 6.0) * dee
    else:
        d31 = (Bn[3, 1] - b * (Bn[0, 1] - Bn[3, 3] + b * Bn[0, 2]) + (pxm * pym + pxi * pyi + px * py) / 6.0) * dee
    d33 = (Bn[3, 3] + b * (-2.0 * Bn[0, 3] + b * Bn[0, 0]) + (pym ** 2 + pyi ** 2 + py ** 2) / 6.0) * dee
    Bn[1, 1] += d11
    Bn[3, 1] += d31
    Bn[1, 3] += d31
    Bn[3, 3] += d33
    Bn[5, 1] += pxm * dee * 0.5
    Bn[1, 5] += pxm * dee * 0.5
    Bn[5, 3] += pym * dee * 0.5
    Bn[3, 5] += pym * dee * 0.5
    Bn[5, 5] += dee
    return True, Dn, Bn
