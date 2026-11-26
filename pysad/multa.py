"""
MULT (ANGLE≠0) の追跡。K64 の tmulta.f (tmulta)。tmulti が phia≠0 のとき呼ぶ。

曲げ磁石 (tbend の分割スライス: K1=0 は k_bend0、K1≠0 は k_bendi) と、曲率つきの多極キック (gkn の表) を、
ndiv 回繰り返す。1 スライスは「半スライス (入口) → [キック → 全スライス]… → キック → 半スライス (出口)」の順。

対応: Kn, SKn (n≤21)、ANGLE、E1 E2 AE1 AE2、FB1 FB2 FRINGE (FRMD)、DISFRIN、DX DY ROTATE/ROTATION DROTATE CHI2、RAD (放射、揺らぎ)、EPS。
未対応: BZ≠0 (K64 自体が abort)、VOLT (RF)、F1 (MULT の F1 は tsetfringep 用)、TAPER、光子 (PHOTONS)、
        スキュー双極 SK0 (K64 の ndiv 計算が fact(-1) を読む)。
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from . import track6n as _n
from .track6n import k_bend0, k_bendi, setr0_1, sqrt1, ndivrad, tbshift1, tbrot1

NJ = dict(cache=True, error_model='numpy')
NMULT = 21
AMPMAX = 0.05
EPS00 = 0.005
NDIVMAX = 2000

_FACT = np.array([math.factorial(i) for i in range(NMULT + 2)], dtype=np.float64)


def _make_gkn():
    g = np.zeros((NMULT + 1, NMULT + 1))
    for n in range(NMULT + 1):
        g[n, 0] = 1.0
        for k in range(1, NMULT - n + 1):
            g[n, k] = g[n, k - 1] * float((2 * k - 3) * (2 * k + 1)) / 8.0 / float(k * (k + n + 1))
    return g


GKN = _make_gkn()


def compile_multa(el, ph=None):
    """MULT ANGLE≠0 → (K_MULTA, params)。tffs.f tpara (icMULT) の前処理を写す"""
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if g('BZ') != 0.0:
        raise NotImplementedError('track6n: MULT の ANGLE≠0 と BZ≠0 は K64 も未対応 (abort)')
    if g('VOLT') + g('DVOLT') != 0.0:
        raise NotImplementedError('track6n: ANGLE≠0 の MULT の RF (VOLT) は未対応')
    if g('F1') != 0.0:
        raise NotImplementedError('track6n: ANGLE≠0 の MULT の F1 は未対応')
    if g('TAPER') != 0.0:
        raise NotImplementedError('track6n: MULT の TAPER は未対応')
    al = float(el.L)
    phi = g('ANGLE')
    psi1 = g('E1') * phi + g('AE1')
    psi2 = g('E2') * phi + g('AE2')
    fb1, fb2 = g('FB1'), g('FB2')
    frmd = int(round(g('FRINGE')))
    if frmd != 3 and frmd != 1:
        fb1 = 0.0
    if frmd != 3 and frmd != 2:
        fb2 = 0.0
    w = phi - psi1 - psi2
    if (fb1 != 0.0 or fb2 != 0.0) and al != 0.0 and phi != 0.0:
        al = al - ((phi * fb1)**2 + (phi * fb2)**2) / al / 48.0 * math.sin(0.5 * w) / math.sin(0.5 * phi)
    knl = {int(k): float(v) for k, v in (getattr(el, 'KnL', {}) or {}).items() if v}
    ksl = {int(k): float(v) for k, v in (getattr(el, 'KsL', {}) or {}).items() if v}
    if max(list(knl) + list(ksl) + [0]) > NMULT:
        raise NotImplementedError('track6n: MULT の次数が 21 を超える')
    ak = np.zeros(NMULT + 1, dtype=np.complex128)
    for n, v in knl.items():
        ak[n] += v
    for n, v in ksl.items():
        ak[n] += 1j * v
    if ak[0].imag != 0.0:
        raise NotImplementedError('track6n: ANGLE≠0 の MULT のスキュー双極 SK0 は未対応 (K64 の ndiv が fact(-1) を読む)')
    if al == 0.0:
        raise NotImplementedError('track6n: L=0 で ANGLE≠0 の MULT')
    from . import lattice as _lat
    theta = _lat._bend_theta(el)
    return _n.K_MULTA, (al, phi, psi1, psi2, ak, g('DX'), g('DY'), theta, g('DROTATE'), g('CHI2'),
                        0.5 * float(el.L), 0.5 * phi, g('EPS'), g('RAD') == 0.0, fb1, fb2, frmd, g('DISFRIN') == 0.0)


@njit(**NJ)
def _sw(phi0, psi1, psi2):
    w = phi0 - psi1 - psi2
    sinw = math.sin(w)
    cosw = math.cos(w)
    sqwh = sinw**2 / (1.0 + cosw) if cosw >= 0.0 else 1.0 - cosw
    return cosw, sinw, sqwh, math.sin(phi0 - psi2)


@njit(**NJ)
def _tbend(S, n, al, ak0, phi0, psi1, psi2, ak, fb1, fb2, mfr, fringe, krad, eps0, fluc, rs, p0, h0, anrad, cuc):
    """tbend (dx=dy=theta=dtheta=dchi2=alg=phig=0, ini=False): ak≠0 は tbendi、そうでなければ tbendcore"""
    cosp1, sinp1 = math.cos(psi1), math.sin(psi1)
    cosp2, sinp2 = math.cos(psi2), math.sin(psi2)
    phib = phi0 + ak0
    if ak != 0.0:
        k_bendi(S, n, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, fb1, fb2, mfr, krad, fringe, eps0, True, True, fluc, rs, p0, h0, anrad, cuc)
    else:
        cosw, sinw, sqwh, sinwp1 = _sw(phi0, psi1, psi2)
        k_bend0(S, n, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2, cosw, sinw, sqwh, sinwp1,
                fb1, fb2, mfr, fringe, krad, True, True, eps0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                fluc, rs, p0, h0, anrad, cuc, False)


@njit(**NJ)
def k_multa(S, n, al, phi, psi1, psi2, ak0, dx, dy, theta, dtheta, dchi2, alg, phig, eps0, radok,
            fb1, fb2, mfring, fringe, gkn, fact, rad, fluc, rs, p0, h0, anrad, cuc, b2=False):
    krad = rad and radok
    cost = math.cos(theta)
    sint = math.sin(theta)
    for i in range(n):
        tbshift1(i, S, dx, dy, phi, cost, sint, True)
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg, phig, dtheta, dchi2, True)
    eps = EPS00 if eps0 == 0.0 else EPS00 * eps0
    ak = ak0.copy()
    ak0r = ak[0].real
    ak1r = ak[1].real
    ak[0] = complex(0.0, ak[0].imag)
    ak[1] = complex(0.0, ak[1].imag)
    nmmax = 0
    for m in range(NMULT, -1, -1):
        if ak[m] != 0.0:
            nmmax = m
            break
    nmmin = nmmax
    for m in range(nmmax):
        if ak[m] != 0.0:
            nmmin = m
            break
    ndiv = 1
    for m in range(max(nmmin, 1), nmmax + 1):
        ndiv = max(ndiv, int(math.sqrt(AMPMAX**(m - 1) / 6.0 / fact[m - 1] / eps * abs(ak[m] * al))) + 1)
    phib = phi + ak0r
    if krad:
        ndiv = max(ndiv, ndivrad(phib, ak1r, 0.0, eps0, h0))
    ndiv = min(NDIVMAX, ndiv)
    aln = al / ndiv
    if not b2:          # v2 (B331): al は tffs.f:1368 で L' に補正済み。K64 はここでもう一度引く (二重)
        if fb1 != 0.0:
            aln = aln - (phi * fb1)**2 / al / 48.0 * math.sin(0.5 * (phi - psi1 - psi2)) / math.sin(0.5 * phi) / ndiv
        if fb2 != 0.0:
            aln = aln - (phi * fb2)**2 / al / 48.0 * math.sin(0.5 * (phi - psi1 - psi2)) / math.sin(0.5 * phi) / ndiv
    rho0 = al / phi
    phin = phi / ndiv
    ak0n = ak0r / ndiv
    ak1n = ak1r / ndiv
    akn = np.zeros(NMULT + 1, dtype=np.complex128)
    for m in range(nmmin, nmmax + 1):
        akn[m] = ak[m] / (fact[m + 1] * ndiv)
    if krad:
        for i in range(n):
            setr0_1(i, S)
    rk = np.zeros(NMULT + 1)
    for m_ in range(1, ndiv + 1):
        if m_ == 1:
            mfr = 0
            if mfring == 2:
                mfr = 0
            elif mfring != 0:
                mfr = -1
            # v2 (B340): tbendi (K1≠0) は mfr に関係なく fb1 を入口・fb2 を出口に入れる。K64 は半スライスにも両方渡すので、
            # MULT の途中に FB2 (前半の終わり)・FB1 (後半の始め) の端が入る。v2 は前半に fb2、後半に fb1 を渡さない
            _tbend(S, n, aln * 0.5, ak0n * 0.5, phin * 0.5, psi1, 0.0, ak1n * 0.5, fb1, 0.0 if b2 else fb2, mfr, fringe,
                   krad, eps0, fluc, rs, p0, h0, anrad, cuc)
        else:
            _tbend(S, n, aln, ak0n, phin, 0.0, 0.0, ak1n, 0.0, 0.0, 0, False, krad, eps0, fluc, rs, p0, h0,
                   anrad, cuc)
        for i in range(n):
            pr = 1.0 + S[i, 5]
            x = S[i, 0]
            cx1 = complex(x, S[i, 2])
            csl = complex(0.0, 0.0)
            csr = complex(0.0, 0.0)
            r = rho0 + x
            rk[0] = 1.0 + sqrt1(x / rho0)
            for k in range(1, NMULT - nmmin + 1):
                rk[k] = rk[k - 1] / r
            for m in range(NMULT, -1, -1):
                cl = complex(0.0, 0.0)
                cr = complex(0.0, 0.0)
                for k in range(max(0, m - nmmax), m - nmmin + 1):
                    m1 = m - k
                    cg = gkn[m1, k] * rk[k] * akn[m1]
                    cl = cl + (m + 1) * cg
                    cr = cr + (0.5 - k) * cg
                csl = csl * cx1 + cl
                csr = csr * cx1 + cr
            S[i, 1] = S[i, 1] - ((csr * cx1).real / r + csl.real) / pr
            S[i, 3] = S[i, 3] + csl.imag / pr
        if krad:
            for i in range(n):
                setr0_1(i, S)
    mfr = 0
    if mfring == 1:
        mfr = 0
    elif mfring != 0:
        mfr = -2
    _tbend(S, n, aln * 0.5, ak0n * 0.5, phin * 0.5, 0.0, psi2, ak1n * 0.5, 0.0 if b2 else fb1, fb2, mfr, fringe,
           krad, eps0, fluc, rs, p0, h0, anrad, cuc)
    for i in range(n):
        if dtheta != 0.0 or dchi2 != 0.0:
            tbrot1(i, S, alg - al, phig - phi, dtheta, dchi2, False)
        tbshift1(i, S, -dx, -dy, -phi, cost, -sint, False)
