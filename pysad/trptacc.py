"""
TRPT (輸送線) の加速: 基準運動量の分布 gammab と、加速のある CAVI / MULT(VOLT) の 4D Twiss 用の写像。

SAD:
- tfgeo.f tfgeo1 (acconly): trpt のとき CAVI/MULT ごとに tfacc で基準運動量 p (単位 amass) を進め gammab(i+1) に入れる。
  tfacc: v = |charge|·VOLT/amass、dh = max(oneev−h1, −v·sin(PHI·sign(charge)))、h2 = h1+dh、p2 = p1 + dh (h2+h1)/(p2'+p1)。
- qtwiss.f: CAVI は qcav → tcave (rfsw = rfsw .and. trpt)、MULT (VOLT≠0) は qmult1 → tmulte (acc = vc≠0 .and. rfsw)。
  どちらも gammab(l) → gammab(l+1) の加速: キックごとの dhg = (p1²−p0²)/(h1+h0)/ndiv で基準粒子の速度補正 (tgetdvh) を z に足し、
  要素の末尾で rg2 = gammab(l)/gammab(l+1) を px, py, δ の行と cod に掛ける (運動量の規格化を新しい基準に移す)。
- qcopymatg: 4×5 の行列 (x px y py; δ 列) は sqrt(gammab(l+1)/gammab(l)) 倍 (断熱減衰を外してユニモジュラにする)。
- qmat2twiss (非結合の枝): 入口の分散 (ex epx ey epy) に gr = gammab(l)/gammab(l+1) を掛けてから伝搬する。

ここでは tcave.f / tmulte.f の加速部分 (dx=dy=θ=0、V1/V20/V11/V02 キーワード 0、autophi=F、放射なし、wake なし) を
解析ヤコビアン (Fortran の trans1) ごと直訳する。MULT は K_n=0 (加速管) だけ (K≠0 かつ VOLT≠0 は未対応 → NotImplementedError)。
DX/DY/ROTATE: CAVI は tchge (入口は回してからずらす)、MULT は tsolrote (ずらしてから回す、行列は数値微分)。
"""
import math
from typing import List

import numpy as np

from . import lattice as _lat
from .emit6d import _Phys, _State, _p2h, _h2p, _tdrife0, CVELOC, ELMASS

ONEEV = 1.0 + 3.83e-12


# ---------------------------------------------------------------- 要素の RF パラメータ
def _g(el, k):
    return float(getattr(el, k, 0.0) or 0.0)


def rf_params(el):
    """(vc, harm, freq, phi): CAVI は VOLT (+DVOLT)、HARM か _harm、FREQ、PHI。MULT も同じキー"""
    tn = type(el).__name__
    if tn not in ('Cavity', 'Multipole'):
        return 0.0, 0.0, 0.0, 0.0
    vc = _g(el, 'VOLT') + _g(el, 'DVOLT')
    harm = _g(el, 'HARM') or _g(el, '_harm')
    return vc, harm, _g(el, 'FREQ'), _g(el, 'PHI')


def is_acc_element(el) -> bool:
    """TRPT で基準運動量を変える要素 (CAVI / MULT で VOLT≠0)"""
    return rf_params(el)[0] != 0.0


def tfacc(el, p1: float, h1: float, charge: float, amass: float):
    """tfgeo.f tfacc (dir=.true.): 要素 el の出口の基準運動量 (単位 amass) と h"""
    vc, harm, freq, phi = rf_params(el)
    v = abs(charge) * vc / amass
    if v == 0.0:
        return p1, h1
    phic = phi * (1.0 if charge >= 0 else -1.0)
    dh = max(ONEEV - h1, -v * math.sin(phic))
    h2 = h1 + dh
    p2 = _h2p(h2)
    p1 = p1 + dh * (h2 + h1) / (p2 + p1)
    return p1, h2


def gammab_profile(elements: List, p0: float, charge: float = 1.0, amass: float = ELMASS) -> np.ndarray:
    """tfgeo1(acconly): gammab(1..nlat)。gammab(1)=p0 (単位 amass)、CAVI/MULT のたびに tfacc"""
    n = len(elements)
    gb = np.empty(n + 1)
    gb[0] = p0
    p1 = p0
    h1 = _p2h(p1)
    for i, el in enumerate(elements):
        if is_acc_element(el):
            p1, h1 = tfacc(el, p1, h1, charge, amass)
        gb[i + 1] = p1
    return gb


# ---------------------------------------------------------------- tcavfrie (解析ヤコビアン付き)
def _cavfrie_piece(st: _State, al, v, w, dphis, s0, p0):
    """tcave.f tcavfrie (autophi=F)。出口は v → −v で呼ぶ"""
    vf = v / al / p0 * 0.5

    def f(c):
        c = c.copy()
        dp1r = c[5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        h1 = _p2h(p1)
        v1 = p1 / h1
        t = -c[4] / v1
        ph = w * t - dphis
        sph, cph = math.sin(ph), math.cos(ph)
        dpt = vf * (sph + s0)
        c[1] += c[0] * dpt
        c[3] += c[2] * dpt
        wc = -w * vf
        dha = wc * (c[0]**2 + c[2]**2) * 0.5
        dh = dha * cph * p0
        h2 = h1 + dh
        a = dh * (h1 + h2) / p1**2
        dpr = a / (1.0 + math.sqrt(1.0 + a))
        dp2r = dp1r + p1r * dpr
        p2r = 1.0 + dp2r
        v2 = p2r * p0 / h2
        c[5] = dp2r
        c[4] = c[4] * v2 / v1
        return c

    def jac(c):
        dp1r = c[5]
        p1r = 1.0 + dp1r
        p1 = p0 * p1r
        h1 = _p2h(p1)
        v1 = p1 / h1
        t = -c[4] / v1
        ph = w * t - dphis
        sph, cph = math.sin(ph), math.cos(ph)
        dpt = vf * (sph + s0)
        wc = -w * vf
        dha = wc * (c[0]**2 + c[2]**2) * 0.5
        dh = dha * cph * p0
        h2 = h1 + dh
        a = dh * (h1 + h2) / p1**2
        dpr = a / (1.0 + math.sqrt(1.0 + a))
        dp2r = dp1r + p1r * dpr
        p2r = 1.0 + dp2r
        v2 = p2r * p0 / h2
        J = np.eye(6)
        at = c[4] / p1r / h1**2
        J[1, 0] = dpt
        J[1, 4] = -c[0] * vf * w / v1 * cph
        J[1, 5] = -J[1, 4] * at
        J[3, 2] = dpt
        J[3, 4] = -c[2] * vf * w / v1 * cph
        J[3, 5] = -J[3, 4] * at
        J[5, 0] = wc * c[0] / v2 * cph
        J[5, 2] = wc * c[2] / v2 * cph
        J[5, 4] = w * dha * sph / v1 / v2
        J[5, 5] = v1 / v2 - J[5, 4] * at
        az = c[4] / v1 / h2**3 * p0
        J[4, 0] = az * J[5, 0]
        J[4, 2] = az * J[5, 2]
        J[4, 4] = v2 / v1 + az * J[5, 4]
        J[4, 5] = az * J[5, 5] - at * v2 / v1
        return J

    st.piece(f, symp_z=False, jac=jac)


# ---------------------------------------------------------------- 加速キック (tcave / tmulte の acc 部分)
def _acc_kick_piece(st: _State, p0, vn, w, wi, phic, phis, offset1, v20a, v02a):
    def parts(c):
        p1 = p0 * (1.0 + c[5])
        h1 = _p2h(p1)
        h1 = p1 + 1.0 / (h1 + p1)
        v1 = p1 / h1
        t = -c[4] / v1
        va = vn + v20a * c[0]**2 + v02a * c[2]**2
        phii = w * t + phic - phis
        sp, cp = math.sin(phii), math.cos(phii)
        dh = max(ONEEV - h1, -va * (sp + offset1))
        h2 = h1 + dh
        p2 = _h2p(h2)
        pf = (h2 + h1) / (p2 + p1) * dh
        p2 = p1 + pf
        v2 = p2 / h2
        av = -(cp * wi - offset1 * t) / p0
        return p1, h1, v1, t, va, sp, cp, dh, h2, p2, pf, v2, av

    def f(c):
        c = c.copy()
        p1, h1, v1, t, va, sp, cp, dh, h2, p2, pf, v2, av = parts(c)
        c[1] += 2.0 * v20a * c[0] * av
        c[3] += 2.0 * v02a * c[2] * av
        c[5] += pf / p0
        c[4] = -t * v2
        return c

    def jac(c):
        p1, h1, v1, t, va, sp, cp, dh, h2, p2, pf, v2, av = parts(c)
        a = -va * w * cp
        dpxa = 2.0 * v20a * c[0]
        dpya = 2.0 * v02a * c[2]
        dav = -(sp + offset1) / p0
        davdz = dav / v1
        davdp = dav / p1 / h1**2 * t * p0
        J = np.eye(6)
        J[1, 0] = 2.0 * v20a * av
        J[1, 4] = dpxa * davdz
        J[1, 5] = dpxa * davdp
        J[3, 2] = 2.0 * v02a * av
        J[3, 4] = dpya * davdz
        J[3, 5] = dpya * davdp
        ddhdx = -(2.0 * v20a * c[0]) * (sp + offset1)
        ddhdy = -(2.0 * v02a * c[2]) * (sp + offset1)
        ddhdz = -a / v1
        ddhdp = -a * t / p1 / h1**2
        J[5, 0] = ddhdx / v2 / p0
        J[5, 2] = ddhdy / v2 / p0
        J[5, 4] = ddhdz / v2 / p0
        J[5, 5] = (v1 + ddhdp) / v2
        J[4, 0] = -ddhdx * t / p2 / h2**2
        J[4, 2] = -ddhdy * t / p2 / h2**2
        J[4, 4] = v2 / v1 - ddhdz * t / p2 / h2**2
        J[4, 5] = -t * p0 * (ddhdp - dh * (h2 * (h2 + h1) + p1**2) / p1 / h1**2) / p2 / h2**2
        return J

    st.piece(f, symp_z=False, jac=jac)


def _tgetdvh(p0, h0, dh):
    """tturne.f tgetdvh: 基準粒子の加速 dh (単位 amass の h) による速度補正 dv"""
    if dh == 0.0:
        return 0.0
    h1 = h0 + dh
    p1 = _h2p(h1)
    return -(p0 + p1) / (p1 * h0 + p0 * h1) / p1 / p0 * dh


def acc_step(el, cod_in, gb0: float, gb1: float, charge: float = 1.0, amass: float = ELMASS,
             omega0: float = 0.0, trf0: float = 0.0, radcod: bool = False, rfsw: bool = True):
    """加速のある CAVI / MULT(VOLT, K=0) の 1 要素: (6×6 の cod まわりのヤコビアン, 出口 cod)。
    tcave.f (CAVI) / tmulte.f の acc 部分 (MULT) の直訳。gb0 = gammab(l)、gb1 = gammab(l+1) (単位 amass)。
    基準運動量の規格化は出口で rg2 = gb0/gb1 を px, py, δ に掛けて gb1 に移す (trans も同じ)。
    行列は qcopymatg の sqrt(gb1/gb0) 倍は掛けていない (呼び出し側 ins_twiss で掛ける)"""
    tn = type(el).__name__
    vc, harm, freq, phi = rf_params(el)
    al = float(el.L)
    is_mult = tn == 'Multipole'
    if is_mult:
        knl = {int(k): float(v) for k, v in (getattr(el, 'KnL', {}) or {}).items() if v}
        ksl = {int(k): float(v) for k, v in (getattr(el, 'KsL', {}) or {}).items() if v}
        if knl or ksl or _g(el, 'ANGLE') != 0.0 or _g(el, 'BZ') != 0.0:
            raise NotImplementedError('trptacc: MULT の VOLT と K_n/ANGLE/BZ の併用は未対応')
    cdx, cdy = _g(el, 'DX'), _g(el, 'DY')
    cth = _g(el, 'ROTATION') + (_g(el, 'DROTATE') if is_mult else 0.0)       # MULT: theta2 = ROTATE + DROTATE (+ akang, K1=0 なので 0)
    chg = cdx != 0.0 or cdy != 0.0 or cth != 0.0
    ph = _Phys(gb0 * amass, 0.0, amass=amass, charge=charge)
    ph.omega0 = omega0
    st = _State(np.asarray(cod_in, dtype=float), ph)
    st.want_db = False
    if chg:
        if is_mult:
            # tmulte.f: tsolrote (bz=0, dz=chi=0: ずらして theta2 で回す)。行列は数値微分
            ct_, st_ = math.cos(cth), math.sin(cth)
            st.piece(lambda c: _lat._mult_rot_in(c, cdx, cdy, ct_, st_), symp_z=False)
        else:
            # tcave.f: tchge(dx,dy,0,theta,…,.true.)
            st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True),
                     jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, True))
    st.trf0 = trf0
    st.rfsw = rfsw
    st.radcod = radcod
    p0 = ph.p0
    h0 = ph.h0
    p1g = gb1
    h1g = _p2h(p1g)
    c = CVELOC
    w = omega0 * harm / c if harm != 0.0 else 2.0 * math.pi * freq / c
    wi = 1.0 / w if w != 0.0 else 0.0
    v = vc / amass * abs(charge)
    fringe = _g(el, 'DISFRIN') == 0.0
    mfring = int(round(_g(el, 'FRINGE')))
    if getattr(el, '_rev', False):
        mfring = mfring * (11 + mfring * (2 * mfring - 9)) // 2        # 逆向き: 1↔2
    if is_mult:
        # tmulte: eps = eps00m (0.005)·EPS、ndiv ≤ 300 (nmmax=0 なので多極の分割は 1)
        eps = 0.005 * (_g(el, 'EPS') or 1.0)
        ndiv = 1
        if rfsw:
            ndiv = min(300, max(ndiv, 1 + int(min(abs(w * al), math.sqrt((v * (1.0 / h0 + 1.0 / h1g))**2 / 3.0 / eps)))))
    else:
        eps = 1e-2
        ndiv = 1 + int(min(abs(w * al), math.sqrt((v * (1.0 / h0 + 1.0 / h1g))**2 / 3.0 / eps)))
    vn = v / ndiv
    aln = al / ndiv
    phis = trf0 * w
    phic = phi * (1.0 if charge >= 0 else -1.0)
    v20a = vn * (w * (0.5 / p0 + 0.5 / p1g))**2 / 4.0
    v02a = v20a
    dhg = (p1g - p0) * (p1g + p0) / (h1g + h0) / ndiv if p1g != p0 else 0.0
    if rfsw:
        offset1 = 0.0
        s0 = 0.0
        if al != 0.0 and fringe and mfring >= 0 and mfring != 2:
            _cavfrie_piece(st, al, v, w, phis - phic, s0, p0)
        dgb = 0.0
        for n in range(ndiv):
            if al != 0.0:
                if n == 0:
                    _tdrife0(st, aln * 0.5)
                else:
                    _tdrife0(st, aln)
                    st.cod[4] += _tgetdvh(p0, h0, dgb) * aln
            _acc_kick_piece(st, p0, vn, w, wi, phic, phis, offset1, v20a, v02a)
            dgb += dhg
        if al != 0.0:
            _tdrife0(st, aln * 0.5)
            st.cod[4] += _tgetdvh(p0, h0, dgb) * aln * 0.5
            if fringe and mfring >= 0 and mfring != 1:
                _cavfrie_piece(st, al, -v, w, phis - phic, s0, p0)
    else:
        if al != 0.0:
            _tdrife0(st, al)
    if chg and is_mult:
        ct_, st_ = math.cos(cth), math.sin(cth)
        st.piece(lambda c: _lat._mult_rot_out(c, cdx, cdy, ct_, st_), symp_z=False)     # tmulte: tsolrote の戻し (rg2 の前)
    if dhg != 0.0:
        rg2 = p0 / p1g
        D = np.eye(6)
        D[1, 1] = rg2
        D[3, 3] = rg2
        D[5, 5] = rg2
        st.apply(D)
        st.cod[1] *= rg2
        st.cod[3] *= rg2
        st.cod[5] = (st.cod[5] + 1.0) * rg2 - 1.0
    if chg and not is_mult:
        st.piece(lambda c: _lat._tchge_cod(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False),
                 jac=lambda c: _lat._tchge_jac(c, cdx, cdy, 0.0, cth, 0.0, 0.0, 0.0, False))     # tcave: 出口 tchge (rg2 の後)
    return st.T.copy(), st.cod.copy()


def make_acc(elements, p0: float, charge: float = 1.0, amass: float = ELMASS, omega0: float = 0.0,
             trf0: float = 0.0, radcod: bool = False, rfsw: bool = True):
    """ins_twiss の acc 引数: {'gb': gammab, 'idx': 加速要素の添字の集合, 'step': (el, i, cod) → (M, cod_out)}。
    加速要素が無ければ None"""
    idx = {i for i, el in enumerate(elements) if is_acc_element(el)}
    if not idx:
        return None
    gb = gammab_profile(elements, p0, charge, amass)

    def step(el, i, cod):
        return acc_step(el, cod, gb[i], gb[i + 1], charge=charge, amass=amass, omega0=omega0,
                        trf0=trf0, radcod=radcod, rfsw=rfsw)
    return {'gb': gb, 'idx': idx, 'step': step}
