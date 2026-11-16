"""
MULT (ANGLE≠0) の COD 写像と行列 (Twiss / Emittance[])。K64 の tmultae (tmulta.f 後半)。

分割スライスごとに「曲げ磁石 (tbende) → 曲率つきの多極キック」を繰り返す。曲げは emit6d._tbende_core / lattice._bend_core_cod_matrix
(誤差なしの BEND) に、スライスのパラメータをもった BEND 要素を渡して使う。キックは gkn の表 (multa.GKN) の解析式で、
行列は tmultae の trans1 (成分 (2,1) (2,3) (4,1) (4,3) だけ) をそのまま写す。

追跡 (multa.k_multa) との違い (K64 の tmultae そのまま): ndiv に放射の項がない、キックを (1+δ) で割らない、
rk(0)=sqrt(r/rho0)、フリンジ FB1/FB2 の tblfre は未対応 (NotImplementedError)、
放射あり (Emittance[]) の tradke は m ループの中にあり 1 スライスあたり nmult+1=22 回呼ばれる (そのまま写す。
キックによる角度変化の分しか放射せず、1 回でも結果は同じ (U0・ε とも 16 桁一致を確認))。
"""
from __future__ import annotations

import math

import numpy as np

from .multa import GKN, NMULT, AMPMAX, EPS00, NDIVMAX, _FACT

_NM1 = NMULT + 1
TRADKE_CALLS = _NM1        # K64 の tmultae は tradke を m ループ内 (nmult+1 回/スライス) で呼ぶ


def prep(el):
    """tmultae の前処理 (tpara の MULT と tmultae の先頭)。分割数・スライス長・係数など"""
    from .elements import Bend
    from . import lattice as _lat
    g = lambda k: float(getattr(el, k, 0.0) or 0.0)
    if g('BZ') != 0.0:
        raise NotImplementedError('MULT の ANGLE≠0 と BZ≠0 は K64 も未対応 (abort)')
    if g('VOLT') + g('DVOLT') != 0.0 or g('F1') != 0.0 or g('TAPER') != 0.0:
        raise NotImplementedError('ANGLE≠0 の MULT の RF (VOLT)・F1・TAPER は未対応')
    if g('CHI2') != 0.0:
        raise NotImplementedError('ANGLE≠0 の MULT の CHI2 は Twiss / Emittance[] では未対応')
    al = float(el.L)
    phi = g('ANGLE')
    if al == 0.0:
        raise NotImplementedError('L=0 で ANGLE≠0 の MULT')
    apsi1, apsi2 = g('AE1'), g('AE2')
    psi1r, psi2r = g('E1'), g('E2')               # 比 (tmultae は psi1, psi2 を ANGLE に対する比で受ける)
    fb1, fb2 = g('FB1'), g('FB2')
    mfring = int(round(g('FRINGE')))
    if mfring != 3 and mfring != 1:
        fb1 = 0.0
    if mfring != 3 and mfring != 2:
        fb2 = 0.0
    if fb1 != 0.0 or fb2 != 0.0:
        raise NotImplementedError('ANGLE≠0 の MULT のフリンジ FB1/FB2 は Twiss / Emittance[] では未対応 (tblfre)')
    ak = np.zeros(_NM1, dtype=np.complex128)
    for n, v in (getattr(el, 'KnL', {}) or {}).items():
        ak[int(n)] += float(v)
    for n, v in (getattr(el, 'KsL', {}) or {}).items():
        ak[int(n)] += 1j * float(v)
    if ak[0].imag != 0.0:
        raise NotImplementedError('ANGLE≠0 の MULT のスキュー双極 SK0 は未対応')
    eps = EPS00 * (g('EPS') or 1.0)
    ak0r, ak1r = ak[0].real, ak[1].real
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
        ndiv = max(ndiv, int(math.sqrt(AMPMAX**(m - 1) / 6.0 / _FACT[m - 1] / eps * abs(ak[m] * al))) + 1)
    ndiv = min(NDIVMAX, ndiv)
    aln = al / ndiv
    phib = phi + ak0r
    phin = phi / ndiv
    phibn = phib / ndiv
    ak1n = ak1r / ndiv
    akn = np.zeros(_NM1, dtype=np.complex128)
    for m in range(nmmin, nmmax + 1):
        akn[m] = ak[m] / (_FACT[m + 1] * ndiv)
    psi1 = psi1r * phi + apsi1
    psi2 = psi2r * phi + apsi2
    mf1 = 0 if (mfring == 2 or mfring == 0) else -1
    mf2 = 0 if (mfring == 1 or mfring == 0) else -2
    fringe = g('DISFRIN') == 0.0
    theta = _lat._bend_theta(el)
    dtheta = g('DROTATE')
    P = dict(al=al, phi=phi, ndiv=ndiv, aln=aln, phin=phin, phibn=phibn, ak1n=ak1n, rho0=aln / phin, akn=akn,
             nmmin=nmmin, nmmax=nmmax, psi1=psi1, psi2=psi2, fringe=fringe, mf1=mf1, mf2=mf2,
             dx=g('DX'), dy=g('DY'), theta=theta, dtheta=dtheta,
             alg=0.5 * float(el.L), phig=0.5 * phi)
    if dtheta != 0.0:
        P['dphix'] = phi * math.sin(0.5 * dtheta)**2
        P['dphiy'] = 0.5 * phi * math.sin(dtheta)
    else:
        P['dphix'] = P['dphiy'] = 0.0
    return P


def _bend(P, kind):
    """スライスの BEND: 'first' (半スライス、入口端 psi1)、'mid' (全スライス)、'last' (半スライス、出口端 psi2)"""
    from .elements import Bend
    if kind == 'mid':
        L, phi, phib, ak, e1, e2 = P['aln'], P['phin'], P['phibn'], P['ak1n'], 0.0, 0.0
    else:
        L, phi, phib, ak = P['aln'] * 0.5, P['phin'] * 0.5, P['phibn'] * 0.5, P['ak1n'] * 0.5
        e1, e2 = (P['psi1'], 0.0) if kind == 'first' else (0.0, P['psi2'])
    b = Bend(name='_multa', L=L, ANGLE=phi, K1=(ak / L), E1=e1, E2=e2)
    b.K0 = phib - phi
    b.DISFRIN = 0.0 if (P['fringe'] and kind != 'mid') else 1.0
    return b


def kick(c, P):
    """曲率つきの多極キック (tmultae): 出口 COD と trans1"""
    akn, nmmin, nmmax, rho0 = P['akn'], P['nmmin'], P['nmmax'], P['rho0']
    cx1 = complex(c[0], c[2])
    csl = csr = csxx = csxy = csyy = 0j
    r = rho0 + c[0]
    rk = np.zeros(_NM1)
    rk[0] = math.sqrt(r / rho0)
    for k in range(1, NMULT - nmmin + 1):
        rk[k] = rk[k - 1] / r
    for m in range(NMULT, -1, -1):
        cl = cr = cxx = cxy = cyy = 0j
        for k in range(max(0, m - nmmax), m - nmmin + 1):
            m1 = m - k
            cg = GKN[m1, k] * rk[k] * akn[m1]
            cl += (m + 1) * cg
            cr += (0.5 - k) * cg
            cxx += (0.25 * (4 * k**2 - 1) * cx1 - (2 * k - 1) * (m + 1) * r) * cg
            if m == 0:
                cxy += (1 - 2 * k) * cg
            else:
                cxy += (m + 1) * ((1 - 2 * k) * cx1 + 2 * m * r) * cg
                cyy += m * (m + 1) * cg
        csl = csl * cx1 + cl
        csr = csr * cx1 + cr
        csxx = csxx * cx1 + cxx
        if m != 0:
            csxy = csxy * cx1 + cxy
            csyy = csyy * cx1 + cyy
        else:
            csxy = csxy + cxy
    out = np.array(c, dtype=float)
    out[1] = c[1] - ((csr * cx1).real / r + csl.real)
    out[3] = c[3] + csl.imag
    J = np.eye(6)
    J[1, 0] = -(csxx.real / r**2 + csyy.real)
    J[1, 2] = csxy.imag * 0.5 / r
    J[3, 0] = J[1, 2]
    J[3, 2] = csyy.real
    return out, J


def _chg_in(P, ent):
    from . import lattice as _lat
    dx, dy, th, dth, alg, phig = P['dx'], P['dy'], P['theta'], P['dtheta'], P['alg'], P['phig']
    if ent:
        return (lambda c: _lat._tchge_cod(c, dx, dy, 0.0, th, dth, alg, phig, True),
                lambda c: _lat._tchge_jac(c, dx, dy, 0.0, th, dth, alg, phig, True))
    a, p = alg - P['al'], phig - P['phi']
    return (lambda c: _lat._tchge_cod(c, dx, dy, 0.0, th, dth, a, p, False),
            lambda c: _lat._tchge_jac(c, dx, dy, 0.0, th, dth, a, p, False))


def _shift(c, P):
    c = np.array(c, dtype=float)
    c[1] += P['dphix']
    c[3] += P['dphiy']
    return c


def lattice_step(el, cod, delta):
    """CALC (qtwiss): (M, 出口 COD)。M は各段の行列の積"""
    from . import lattice as _lat
    P = prep(el)
    M = np.eye(6)
    c = np.array(cod, dtype=float)

    def apply(Mn, cn):
        nonlocal M, c
        M = Mn @ M
        c = cn

    f, j = _chg_in(P, True)
    apply(j(c), f(c))
    if P['dtheta'] != 0.0:
        apply(np.eye(6), _shift(c, P))
    for n in range(1, P['ndiv'] + 1):
        c2, Mb = _lat._bend_core_cod_matrix(_bend(P, 'first' if n == 1 else 'mid'), c)
        apply(Mb, c2)
        c2, Jk = kick(c, P)
        apply(Jk, c2)
    c2, Mb = _lat._bend_core_cod_matrix(_bend(P, 'last'), c)
    apply(Mb, c2)
    if P['dtheta'] != 0.0:
        apply(np.eye(6), _shift(c, P))
    f, j = _chg_in(P, False)
    apply(j(c), f(c))
    return M, c


def emit(st, el, enarad=True):
    """tturne / tecalc 用 (emit6d._State に段ごとに積む)"""
    from . import emit6d as _e
    P = prep(el)
    f, j = _chg_in(P, True)
    st.piece(f, jac=j)
    if P['dtheta'] != 0.0:
        st.piece(lambda c: _shift(c, P))
    for n in range(1, P['ndiv'] + 1):
        _e._tbende_core(st, _bend(P, 'first' if n == 1 else 'mid'), enarad)
        st.piece(lambda c: kick(c, P)[0], jac=lambda c: kick(c, P)[1])
        if enarad:
            for _ in range(TRADKE_CALLS):          # K64: tradke は m ループの中 (22 回/スライス)
                st.tradke(P['aln'], 0.0)
    _e._tbende_core(st, _bend(P, 'last'), enarad)
    if P['dtheta'] != 0.0:
        st.piece(lambda c: _shift(c, P))
    f, j = _chg_in(P, False)
    st.piece(f, jac=j)
