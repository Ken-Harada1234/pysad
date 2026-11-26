# -*- coding: utf-8 -*-
"""Common part of the quadrupole body integrators (B226): Strang-split integration of H = H_lin + H_nl as in K64 tsolqu and tsolque.
quadbody.py — 四極の本体の積分器の共通部 (B226)

K64 の四極の本体 (tsolqu・tsolque の bz=0 の枝) は、ハミルトニアン H = H_lin(x, y, px, py) + H_nl(px, py)
(H_lin = (px² + py²)/2 + k(x² − y²)/2 の厳密な線形写像 B、H_nl = −√(1 − p²) − (−1 + p²/2) の運動学的な残り =
「非線形ドリフト」D (運動量は変えない)) を、1 刻み (長さ aln) を D(aln/2) B(aln) D(aln/2) の Strang 分割 (2 次) で積分する。
刻みの数は 1 + int(|L·K1L|/ε) (CALC は ε = 0.1、追跡は 0.2) で、軌道のずれの大きさには依らない。
誤りは (aln²·k)·p² (p = 軌道の傾き) に比例し、軌道が 3 mm・2 mrad ずれた四極 (L = 0.3 m、K1 = −0.6) で 1e-7〜7e-7。

このモジュールは、その分割の構造を 1 か所に集めたもの (以前は 11 か所に同じ構造が重複していた):
- 積分の方式 (`scheme()`): 0 = K64 の Strang 分割 (k64 モード。従来とビット単位で同じ)、
  1 = 4 次の Yoshida の合成 (v2 モード): S(w1·h) S(w0·h) S(w1·h)、w1 = 1/(2 − 2^(1/3))、w0 = 1 − 2·w1 (< 0)。
  D と B はそれぞれ厳密な流れなので、合成は 4 次のシンプレクティックな写像になる。
- 刻みの数 (`ndiv()`・`ndiv_calc()`): 追跡の ε は k64 が 0.2 (EPS 倍)、v2 が 0.1 (四極 1 台の評価は numba で 2〜2.5 倍)。CALC の ε は k64 が 0.1、v2 が 0.05 (ヤコビアンの誤りを 1e-10 台にする)。
- 刻みの中の群 (`group()`): 「D → B」の組と、最後の D。放射の点 (tradk・tradke) は刻みの最後の B のあと、次の刻みの先頭の D の前。
  隣り合う刻みの D は 1 つにまとめる (K64 と同じ)。
- 追跡 (規格化した運動量・スカラー) の D と B (`make_track`): python 版 (track6) と numba 版 (track6n・track6n_ext) は、
  それぞれの sqrt1・tzsetparam0 (k64 のビット列を保つために別の実装) を渡して同じ式から作る。
CALC (正準運動量の COD 配列) の D と B (lattice の `_tqente_cod`・`_tsolque_body_cod`、elkern の `_tqente`・`_tsolque_body` と
ヤコビアン) は、numpy 版と numba 版で丸めが違うので従来のまま。それらの刻みの構造だけが `group()`・`ndiv()` を使う。

互換モード: `scheme()` は compat.v2() なら 1、k64 なら 0。`set_scheme(0|1|None)` で固定できる (試験・診断用)。
numba のカーネルは compat を読めないので、呼び出し側が scheme を引数 (rs[4]) で渡す (`track6n.make_rs`)。
"""
from __future__ import annotations

import math

from numba import njit

from . import compat as _compat

NJ = dict(cache=True, error_model='numpy', nogil=True)

# 4 次の Yoshida の係数
_W1 = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
_W0 = 1.0 - 2.0 * _W1
QB_M = (1, 3)                                    # 1 刻みあたりの本体 B の数 (方式 0・1)
QB_BW = ((1.0, 0.0, 0.0), (_W1, _W0, _W1))       # B の重み (刻みの長さに対する比)

_override = None


def set_scheme(s):
    """Fix the integration scheme to 0 (K64 Strang), 1 (4th-order Yoshida), or None (follow compat); return the previous setting.
    0 (K64 の Strang)・1 (4 次の Yoshida)・None (compat に従う) に固定する。戻り値: 以前の設定"""
    global _override
    old = _override
    _override = s
    return old


def scheme() -> int:
    if _override is not None:
        return int(_override)
    return 1 if _compat.v2() else 0


def ndiv(alak: float, eps: float, sch: int) -> int:
    """Return the number of tracking steps, 1 + int(|al*ak|/eps); v2 halves eps (B226).
    追跡の刻みの数 1 + int(|al·ak|/ε) (al = 長さ、ak = |K1|·長さ、ε = 0.2 (EPS 倍))。v2 は ε を半分 (0.1) にする: 4 次の合成で
    本体の評価が 3 倍、刻みが 1〜2 倍で、1 台の誤りは 1/16〜1/100 (docs/30 B226)。呼び出し側が従来どおり min(1000, ·) などをかける"""
    return 1 + int(abs(alak) / (eps if sch == 0 else eps * 0.5))


def ndiv_calc(alak: float, sch: int) -> int:
    """Return the number of steps for `CALC` (tsolque): eps = 0.1 in k64, and 0.05 in v2.
    CALC (tsolque) の刻みの数。k64 は ε = 0.1 (K64 のまま)、v2 は ε = 0.05 (刻みは 2 倍。3 mm・2 mrad ずれた四極のヤコビアンの誤りが
    7e-8 → 1.7e-10 になる)"""
    return 1 + int(abs(alak) / (0.1 if sch == 0 else 0.05))


def ngroup(ndiv_: int, sch: int) -> int:
    """Return the number of groups (m pairs of 'D then B' times ndiv, plus the final D).
    群の数 (m 個の「D → B」× ndiv と、最後の D)"""
    return QB_M[sch] * ndiv_ + 1


def group(g, aln1, ndiv_, sch):
    """Return the g-th group: lengths of D and B, B coefficient index jb, radiation point, and whether it is the final D.
    g 番目の群: (D の長さ, B の長さ, B の係数の番号 jb, 放射点 rad, 最後の D か)。
    刻みの先頭の D は、最初の刻みだけ半分、以降は前の刻みの最後の半分と合わせて全部 (K64 と同じ)。
    方式 0 では D = 0.5·aln1 (最初)・aln1 (以降)、B = aln1、最後の D = 0.5·aln1 で、従来の式とビット単位で同じ"""
    m = QB_M[sch]
    bw = QB_BW[sch]
    if g == m * ndiv_:
        return aln1 * bw[m - 1] * 0.5, 0.0, 0, True, True
    n = g // m
    j = g - n * m
    if j == 0:
        if n == 0:
            dl = aln1 * bw[0] * 0.5
        else:
            dl = aln1 * bw[0]
    else:
        dl = aln1 * 0.5 * (bw[j - 1] + bw[j])
    jb = 1 if (sch == 1 and j == 1) else 0
    return dl, aln1 * bw[j], jb, (j == m - 1 and n < ndiv_ - 1), False


group_nb = njit(**NJ)(group)         # numba 版 (elkern・track6n などの njit カーネルから呼ぶ。同じ式)
ndiv_nb = njit(**NJ)(ndiv)
ndiv_calc_nb = njit(**NJ)(ndiv_calc)


@njit(**NJ)
def sch_of(rs):
    """Read the integration scheme from the numba kernel state array rs (rs[4]; 0 = K64 if the array is short).
    njit カーネルの状態配列 rs から積分の方式を読む (rs[4]。配列が短ければ 0 = K64)"""
    return 1 if (rs.shape[0] > 4 and rs[4] == 1) else 0


def make_track(sqrt1, tzsetparam0, jit=None):
    """Build the D and B maps for tracking (normalized momentum px/(1+delta), scalar coordinates) as a translation of tsolqu.f.
    追跡 (規格化した運動量 px/(1+δ)、スカラーの座標) の D と B を作る。
    sqrt1・tzsetparam0 は呼び出し側のもの (python 版は track6、numba 版は track6n。k64 のビット列を保つため実装が違う)。
    jit は numba の njit (numba 版) か None (python 版)。式は track6._tsolqu・_tsolqur_ens (K64 の tsolqu.f) の直訳"""
    if jit is None:
        def jit(f):
            return f

    def drift(x, px, y, py, z, ln):
        """非線形ドリフト D (tqente 相当。運動量は変えない)。戻り: x, y, z"""
        ap = px**2 + py**2
        dpz = sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * ln
        x += px * r
        y += py * r
        z -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
        return x, y, z

    drift = jit(drift)
    _group = jit(group)
    _ndiv = jit(ndiv)

    def body(x, px, y, py, z, dv, par, ln, dx0, dy0, qum):
        """線形の本体 B (長さ ln、par = tzsetparam0(δ, ln, K1))。dx0・dy0 は MULT の K0 による軸のずれ。
        qum=True は tsolqum (z の dv の項が −0.25·(ap − dv)·ln。SAD のまま)。戻り: x, px, y, py, z"""
        w1, s1, xs1, dc1, sh2, xsh2, dch2 = par
        pxi, pyi = px, py
        ap = pxi**2 + pyi**2
        xi, yi = x + dx0, y + dy0
        u1 = xi * dc1 + pxi * s1 / w1
        u2 = -xi * w1 * s1 + pxi * dc1
        v1 = yi * dch2 + pyi * sh2 / w1
        v2 = yi * w1 * sh2 + pyi * dch2
        x += u1
        px = pxi + u2
        y += v1
        py = pyi + v2
        if qum:
            z = (z - 0.25 * (w1 * (xi**2 * xs1 - yi**2 * xsh2) + (ap - dv) * ln
                             + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2))
        else:
            z = (z - 0.25 * (w1 * (xi**2 * xs1 - yi**2 * xsh2) + ap * ln
                             + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2)
                 - dv * ln)
        return x, px, y, py, z

    body = jit(body)

    def nonrad(x, px, y, py, z, dp, dv, al, ak, eps0, ak0x, ak0y, qum, sch):
        """tsolqu / tsolqum (bz=0、放射なし、ak = |K1|·L > 0): 戻り x, px, y, py, z"""
        eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
        nd = _ndiv(al * ak, eps, sch)
        aln0 = al / nd
        dx0 = ak0x / ak
        dy0 = ak0y / ak
        akk = ak / al
        pa = tzsetparam0(dp, aln0 * QB_BW[sch][0], akk)
        pb = tzsetparam0(dp, aln0 * QB_BW[sch][1], akk) if sch == 1 else pa
        for g in range(QB_M[sch] * nd):
            dl, bl, jb, _rad, _fin = _group(g, aln0, nd, sch)
            x, y, z = drift(x, px, y, py, z, dl)
            x, px, y, py, z = body(x, px, y, py, z, dv, pb if jb == 1 else pa, bl, dx0, dy0, qum)
        dl, _bl, _jb, _rad, _fin = _group(QB_M[sch] * nd, aln0, nd, sch)
        x, y, z = drift(x, px, y, py, z, dl)
        return x, px, y, py, z

    nonrad = jit(nonrad)
    return drift, body, nonrad
