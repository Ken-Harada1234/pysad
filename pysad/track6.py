"""Particle tracking for `TrackParticles` (translation of tftrack.f / tturn.f `tturn`, `tturn1` in SAD).
track6.py — SAD TrackParticles (tftrack.f → tturn.f tturn / tturn1) の直訳による粒子追跡

SAD の追跡は Emittance[]/CALC (tturne: tquade, tbende, tthine, tcave, tradke) とは別のルーチンを使う:
  ドリフト  tdrift.f  tdrift_free
  四極      tquad.f   tquad → tsolrot (回転), ttfrin (nord=4), tsolqu / tsolqur (放射あり)
  6 極      tquad.f   tthin (α/β/γ の 4 次積分) + ttfrin (nord=6 フリンジ)
  BEND      tbend.f   tbend → tbendcore (入口エッジ・本体・出口エッジを 1 つの閉じた式で)。
                      放射ありは tbendr で ndivrad 分割 (各分割で tbendcore + tradk)
  空洞      tcav.f    tcav (RAD の有無で offset を切り替え)
  放射      photons.f kradlib tradk → tradkn (FLUC OFF)。tradkfn (FLUC ON, 乱数) は未移植
内部座標は (x, px/(1+δ), y, py/(1+δ), z, δ) と粒子ごとの dv (tconvm)。入出力は px/p0。

未対応 (NotImplementedError): K1≠0 の BEND、BEND/四極/6 極のずれ・回転 (四極 ROTATE は可)、
  長さのある空洞、四極の F1/FRINGE、ソレノイド、MULT など。FLUC (揺らぎ) は未対応。
"""
from __future__ import annotations

import copy
import math
import sys
import time
from typing import List

import numpy as np

from . import emit6d as _e
from . import compat as _compat
from . import lattice as _lat
from . import quadbody as _qb
from .elements import Drift, Bend, Quadrupole
from .sadrandom import SadRandom, tdusrn, tdusrnpl
from .perfkern import edge_len as _edge_len

CVELOC = 299792458.0
_GMIN = -0.999999        # photons.f module tspin の gmin (tradkn の max(gmin, ...))

# 乱数 (SAD 既定プラグイン "SAD")。FLUC の追跡はここから引く。track_particles(seed=...) で SeedRandom[seed] 相当
RNG = SadRandom(17)
_FLUC = False            # rfluct (RAD と FLUC が ON)
_POL = False             # calpol (POL と RAD が ON): 粒子ごとのスピンを追う (tradk の sprot、bsi)
_RADPOL = True           # radpol (放射による偏極 Sokolov-Ternov)
_GSPIN = 0.001159652181280002      # tspin: g-2 (spin.GSPIN と同じ)
_CL = 1.0 + _GSPIN
_PST = 8.0 * math.sqrt(3.0) / 15.0
_CFPD = 5.0 * math.sqrt(3.0) / 8.0


# ---------------------------------------------------------------------------
# mathfun
# ---------------------------------------------------------------------------

def _sqrt1(a):
    """mathfun sqrt1: sqrt(1+a) − 1"""
    return a / (1.0 + math.sqrt(1.0 + a))


def _pxy2dpz(px, py):
    """mathfun pxy2dpz: sqrt(1 − px² − py²) − 1"""
    return _sqrt1(-(px**2 + py**2))


def _xsincos(x):
    """mathfun xsincos: (sin, x − sin, cos, cos − 1)"""
    return _lat._xsincos(x)


def _xsincosh(x):
    """mathfun xsincosh: (sinh, x − sinh, cosh, cosh − 1)  (tfloor.f)"""
    sh = math.sinh(x)
    ch = math.cosh(x)
    return sh, x - sh, ch, sh**2 / (1.0 + ch)


def _g(el, k):
    return float(getattr(el, k, 0.0) or 0.0)


# ---------------------------------------------------------------------------
# 粒子の状態
# ---------------------------------------------------------------------------

class _Particle:
    """c = 内部座標 (x, px/(1+δ), y, py/(1+δ), z, δ)、dv、放射の基準 (pxr0, pyr0, zr0)。"""

    def __init__(self, cod, ph):
        c = np.array(cod, dtype=float)
        pr = 1.0 + c[5]
        c[1] /= pr
        c[3] /= pr
        self.c = c
        self.ph = ph
        self.dv = ph.dv(c[5])          # tconvm (dv = ... + dvfs: ph.dv が dvfs を含む)
        self.s = np.array([0.0, 0.0, 1.0])     # スピン (sx, sy, sz)。POL のときだけ使う
        self.bsi = 0.0                          # kradlib bsi (磁場の縦成分の積分: tradk のスピン歳差に使う)
        self.bpadd = 0.0                        # v2 (B165): 本体の運動量の向きの磁場の積分 ∫B∥dl/(Bρ) (粒子の運動量で規格化)。次の tradk の sprot が運動量の向きに足す
        self.bpar = None                        # v2 (C2): 運動量の向きの磁場の積分 ∫B∥dl/(Bρ)。None でなければ sprot は bsi の代わりにこれを使う
        self.setr0()

    def setr0(self):
        self.pxr0 = self.c[1]
        self.pyr0 = self.c[3]
        self.zr0 = self.c[4]

    def out(self):
        c = self.c.copy()
        pr = 1.0 + c[5]
        c[1] *= pr
        c[3] *= pr
        return c


_BPSIGN = 1.0


def _steer_edge_exact(p, rhoc, psi, entrance):
    """v2 (B221): ANGLE=0 の BEND (tsteer) の端の面の厳密な写像 (CALC の tsteee と同じ lattice._tbedge_cod)。追跡の座標 (px/(1+δ)) と正準運動量の変換つき"""
    c = p.c
    pr = 1.0 + c[5]
    cc = np.array([c[0], c[1] * pr, c[2], c[3] * pr, c[4], c[5]])
    cc = _lat._tbedge_cod(cc, rhoc, psi, entrance)
    pr = 1.0 + cc[5]
    c[0], c[1], c[2], c[3], c[4] = cc[0], cc[1] / pr, cc[2], cc[3] / pr, cc[4]


def _zr0_edge(p, tan, sgn):
    """v2 (B220): 端の面 (軸に対し角 e の平面、tan = tan e) の手前 (sgn = -1: 入口) か先 (+1: 出口) の自由な飛行の行程 ℓ。
    K64 の線形の式 tan·x の px の 1 次・p/pz の因子を入れた厳密な式 (perfkern.edge_free)。tradk の al1 = al − z + zr0 から引く"""
    c = p.c
    pr = 1.0 + c[5]
    return _edge_len(tan, 0.0, c[0], c[2], c[1] * pr, c[3] * pr, c[5], sgn)


def _sprot(p, pxm, pym, bx0, by0, bz0, bsi, h, pp, ucde, cphi0, sphi0):
    """photons.f sprot: 1 粒子のスピン歳差 (磁場の方向 (bx0, by0, bz0) = 運動量の変化、bsi = 縦の磁場)。
    K64 のまま: 最後の sx 更新の後の sz は更新後の sx を使う"""
    s = p.s
    sx, sy, sz = s[0], s[1], s[2]
    pzm = 1.0 + _pxy2dpz(pxm, pym)
    if _compat.v2():
        # v2 (C2): (bx0, by0, bz0) は刻みの前後の向き (単位ベクトル) の外積で、大きさは sin θ (θ は向きの変化の角)。K64 は
        # 中点の pz (≈ cos(θ/2)) で割るので歳差の角が (1+aγ)·2 sin(θ/2) になり、刻みが粗いと θ²/24 だけ小さい
        # (偏向角 0.26 rad の BEND を 3 刻みで 3e-4)。一様な磁場の刻みで厳密な θ にする
        ppa_ = math.sqrt(bx0 * bx0 + by0 * by0 + bz0 * bz0)
        if ppa_ > 1e-12:
            f_ = math.asin(min(1.0, ppa_)) / ppa_ * pzm
            bx0, by0, bz0 = bx0 * f_, by0 * f_, bz0 * f_
    bx = bx0 / pzm
    by = by0 / pzm
    if p.bpar is None and p.bsi != 0.0 and _compat.v2() and _compat.spin_bpar():
        # v2 (B242): 端の縦磁場の積分 bsi (p0 で規格化した ψ) は、粒子の運動量で規格化し直す (回転角 (1+a)ψ/(1+δ))。K64 は割らないので、本体の Δψ/(1+δ) (B165)
        # と合わせた和が δΔψ だけ 0 にならない (δ = 0.02 の QUAD で R01 に 1.2e-7)
        bsi = bsi / (1.0 + p.c[5])
    if p.bpar is not None:
        # v2 (C2、SOL 領域・UND): 向きの変化から出した (bx, by, bz) は運動量に垂直な磁場の積分 B⊥。運動量の向きの成分は
        # p.bpar = ∫B∥dl/(Bρ) (粒子の運動量で規格化) を運動量の向きに置く。
        # K64 の式 (下) は bsi = bz·al (p0 で規格化) を z 軸の向きに足すので、(1) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならず、
        # (2) 軌道が傾いていると bsi·ẑ の運動量に垂直な成分が B⊥ に二重に入る
        bz = bz0 / pzm
        bp = bx * pxm + by * pym + bz * pzm
        btx = bx - bp * pxm
        bty = by - bp * pym
        btz = bz - bp * pzm
        blx = p.bpar * pxm
        bly = p.bpar * pym
        blz = p.bpar * pzm
    else:
        bz = bz0 / pzm + bsi
        bp = bx * pxm + by * pym + bz * pzm
        blx = bp * pxm
        bly = bp * pym
        blz = bp * pzm
        btx = bx - blx
        bty = by - bly
        btz = bz - blz
    if p.bpadd != 0.0:
        blx += p.bpadd * pxm
        bly += p.bpadd * pym
        blz += p.bpadd * pzm
    ct = 1.0 + h * _GSPIN
    gx = ct * btx + _CL * blx
    gy = ct * bty + _CL * bly
    gz = ct * btz + _CL * blz
    if ucde > 0.0 and _RADPOL:
        bt = math.sqrt(btx**2 + bty**2 + btz**2)
        if bt != 0.0:
            st = (sx * btx + sy * bty + sz * btz) / bt
            dst = (_PST - st) * _CFPD * (pp * p.ph.p0 / h**2)**2 * ucde / bt
            sx = sx + dst * btx
            sy = sy + dst * bty
            sz = sz + dst * btz
    g = math.sqrt(gx**2 + gy**2 + gz**2)
    if g != 0.0:
        sinu, _xs, cosu, dcosu = _xsincos(g)
        sw = -(sx * gx + sy * gy + sz * gz) * dcosu / g**2
        sinu = sinu / g
        sux = sy * gz - sz * gy
        suy = sz * gx - sx * gz
        suz = sx * gy - sy * gx
        sx, sy, sz = cosu * sx + sinu * sux + sw * gx, cosu * sy + sinu * suy + sw * gy, cosu * sz + sinu * suz + sw * gz
    sx = sx * cphi0 + sz * sphi0
    sz = (sz - sx * sphi0) / cphi0
    s[0], s[1], s[2] = sx, sy, sz


# ---------------------------------------------------------------------------
# 光子の記録 (PHOTONS、photons.f の tphrec・tphotonlist、synradcl.f の radangle1・BK1323)。BEAMBEAM 以外の要素の光子 (v2next B3)
# ---------------------------------------------------------------------------

_PH = None        # PHOTONS: dict(sink=list, geo=GeoResult, l=要素の番号 (1 始まり), ctx=pcvt|None)。None なら光子を記録しない


def _bk1323(x):
    """synradcl.f BK1323: BesselK(1/3, x)/BesselK(2/3, x) の近似 (相対誤差 3.83e-5、x<0 は 0)"""
    if x >= 0.8579:
        x1 = 1.0 / x
        return (((0.010874519 * x1 - 0.042281332) * x1 + 0.083425502) * x1 - 0.16539043) * x1 + 1.0
    if x >= 0.3551:
        return (((-0.37018753 * x + 1.1671599) * x - 1.4970125) * x + 1.0262830) * x + 0.55852340
    if x > 0.0:
        x1 = x ** 0.333333333333333333
        x2 = x1 * x1
        return ((((0.81410785 * x2 - 1.6504717) * x2 + 1.7410516) * x2 - 1.4937488) * x2 + 1.5701681) * x1
    return 0.0


def _radangle1(gam, rho, u, compton):
    """synradcl.f radangle1 (2023-03-19 版): 光子 (エネルギー u、ビームエネルギーで規格化) の放射角 (thx, thy) と偏光 (xi3, xi2)。tgauss を 2 回引く"""
    upsilon = compton * gam**2 / rho
    a = math.log(u / upsilon / 1.5)
    sigg = math.exp(-0.462627729148672 + a * (-0.4883575801292511 + a * -0.020021831341729676))
    thyg = sigg * RNG.tgauss()
    thy = thyg / gam
    thx = sigg * RNG.tgauss() / gam
    thyg2 = math.sqrt(1.0 + thyg**2)
    v = 1.0 / math.sqrt(1.0 + (thx * gam)**2 + thyg**2)
    x = 1.5 * upsilon * u * v**3
    z = x / (3.0 * upsilon) * thyg2**3
    r = -thyg / thyg2 * _bk1323(z)
    return thx, thy, (1.0 - r * r) / (1.0 + r * r), 2.0 * r / (1.0 + r * r)


def _ph_set(dx, dy, theta, dtheta, phi0, phig, al, fr0=0.0):
    """tsetpcvt: 光子の位置を決める要素の情報 (ずれ・ロール・曲げの半角・長さ) と、今の刻みの要素内の位置 fr0 (0..1)"""
    if _PH is None:
        return
    if phi0 == 0.0:
        cost, sint = math.cos(theta + dtheta), math.sin(theta + dtheta)
    else:
        cost, sint = math.cos(theta), math.sin(theta)
    _PH['ctx'] = dict(dx=dx, dy=dy, theta=theta, dtheta=dtheta, phi0=phi0, phig=phig, al=al, cost=cost, sint=sint, fr0=fr0)


def _ph_adv(frac):
    """pcvt%fr0 += al1/al (tmulti.f・tsolqu.f: 全粒子の tradk の後に 1 回)"""
    if _PH is not None and _PH.get('ctx') is not None:
        _PH['ctx']['fr0'] += frac


def _ph_rec(p, xi, pxi, yi, pyi, dp, ppx, ppy, p1, h1, rho, ds):
    """photons.f tphrec: 光子 1 個 (エネルギー dp = ビームの運動量で規格化、xi, yi: 放射点の位置) を PhotonList の行に足す。
    行 = [E(eV), gx, gy, gz, ux, uy, uz, xi1, xi2, xi3, 粒子番号, 要素の番号]。
    K64 は偏光の回転で c1 = (cos·ppx1 + sin·ppy1)/pp と、すでに規格化した ppx1・ppy1 をさらに pp (≈ 曲げ角) で割る (xi1 が ±1 を超えて 1e3〜1e5 になる)。
    v2 は割らない (|xi1|, |xi3| ≤ 1)"""
    ctx = _PH['ctx']
    ph = p.ph
    l = _PH['l']
    if ctx is None:
        # この要素 (MULT・SOL 内など) は光子の位置の情報を持たない (未対応): 乱数だけ K64 と同じに引いて記録しない
        _radangle1(h1, rho * p1 / ph.p0, dp, ph.cuc / 1.5)
        return
    al = ctx['al']
    fr = ctx['fr0'] if al == 0.0 else ctx['fr0'] + ds / al
    gv = _PH['geo'].geo_frac(l, fr)
    compton = ph.cuc / 1.5                       # rclassic/rcratio = ħ/(mc) (簡約コンプトン波長)
    # radangle1(h1, rho*p1/p0, dp, thu, thv, xi30, xi2): 引数の順は (THX, THY, XI3, XI2) = (thu, thv, xi30, xi2)
    thu, thv, xi30, xi2 = _radangle1(h1, rho * p1 / ph.p0, dp, compton)
    ppn = math.hypot(ppx, ppy)
    ppx1 = ppx / ppn
    ppy1 = ppy / ppn
    pxir = pxi + (thv * ppx1 + thu * ppy1)
    pyir = pyi + (-thv * ppy1 + thu * ppx1)
    cod = _lat._tchge_cod([xi, pxir, yi, pyir, 0.0, 0.0], ctx['dx'], ctx['dy'], 0.0, ctx['theta'], ctx['dtheta'], 0.0, ctx['phig'], False)
    from . import geo as _geo
    gv1 = _geo.tforbitgeo(gv, cod)
    if _compat.v2():
        # v2: 低い γ・軟らかい光子の放射角 (σ/γ が 1 rad を超える) で px² + py² ≥ 1 になり pz が定義できない。単位球の上に載せる
        # (K64 の radangle1 は γ ≫ 1 の近似。数 GeV の電子では起きない)
        _pt2 = cod[1] ** 2 + cod[3] ** 2
        if _pt2 > 0.999999:
            _f = math.sqrt(0.999999 / _pt2)
            cod[1] *= _f
            cod[3] *= _f
    pzi = 1.0 + _pxy2dpz(cod[1], cod[3])
    dpa = dp * (gv1[:, 0:3] @ np.array([cod[1], cod[3], pzi]))
    if _compat.v2():
        c1 = ctx['cost'] * ppx1 + ctx['sint'] * ppy1
        s1 = -ctx['sint'] * ppx1 + ctx['cost'] * ppy1
    else:
        c1 = (ctx['cost'] * ppx1 + ctx['sint'] * ppy1) / ppn
        s1 = (-ctx['sint'] * ppx1 + ctx['cost'] * ppy1) / ppn
    c2 = (c1 - s1) * (c1 + s1)
    s2 = 2.0 * c1 * s1
    dn = float(np.linalg.norm(dpa))
    # K64 の tphotonlist は dp·amass·gammab(l) (K64 の gammab は p/(mc) = βγ で γ ではない: 相対 1/(2γ²)、3 GeV で 1.45e-8。D1)。v2 は γ
    _g_e = ph.h0 if _compat.v2() else ph.p0
    if _compat.v2() and ph.anrad < 0.0:
        _g_e = ph.p0 * ph.p0 / ph.h0          # v2 (B233): β < 1 の粒子: 光子のエネルギー = 粒子が失うエネルギー = β·Δp c (dE = v dp)。dp は運動量の損失 (δ 単位)
    _PH['sink'].append([dn * ph.amass * _g_e, float(gv1[0, 3]), float(gv1[1, 3]), float(gv1[2, 3]),
                        float(dpa[0] / dn), float(dpa[1] / dn), float(dpa[2] / dn),
                        float(c2 * xi30), float(xi2), float(s2 * xi30), float(p.idx + 1), float(l)])


def _tradk(p: _Particle, al, phi0):
    """photons.f kradlib tradk → tradkn (rfluct=F, calpol=F)。最後に基準点を更新。"""
    _tradk_core(p, al, phi0)
    p.setr0()
    p.bpadd = 0.0
    if _POL:
        p.bsi = 0.0                    # tradk: calpol なら bsi=0


def _tradk_core(p: _Particle, al, phi0):
    """photons.f tradkn / tradk1 (FLUC OFF) と tradkfn / tradkf1 (FLUC ON) の 1 粒子分 (基準点は更新しない)。
    FLUC ON では tdusrn で光子数とエネルギーを乱数で決める (粒子ごとに必ず乱数を 1 回以上引く)。"""
    if al != 0.0:
        ph = p.ph
        c = p.c
        cphi0 = math.cos(phi0)
        sphi0 = math.sin(phi0)
        dpz0 = _pxy2dpz(p.pxr0, p.pyr0)
        px0 = cphi0 * p.pxr0 + sphi0 * (1.0 + dpz0)
        dpz0 = cphi0 * dpz0 - sphi0 * p.pxr0
        dpx = c[1] - px0
        dpy = c[3] - p.pyr0
        dpz = _pxy2dpz(c[1], c[3])
        if _FLUC:
            dpz0 = _pxy2dpz(px0, p.pyr0)          # tradkfn: 回した後に計算し直す
        ppx = c[3] * dpz0 - dpz * p.pyr0 + dpy
        ppy = dpz * px0 - c[1] * dpz0 - dpx
        ppz = c[1] * p.pyr0 - c[3] * px0
        ppa = math.sqrt(ppx**2 + ppy**2 + ppz**2)
        theta = _e._asinz(ppa)
        pr = 1.0 + c[5]
        pp = ph.p0 * pr
        h1 = _e._p2h(pp)
        al1 = al - c[4] + p.zr0
        anp = abs(ph.anrad) * h1 * theta
        fb = 1.0
        if ph.anrad < 0.0:
            # v2 (B167): 厳密な放射のスペクトル (β < 1)。光子数に a、光子エネルギーに b (radbeta.ab)
            from .radbeta import ab as _rb_ab
            fa, fb, _dla, _dlb = _rb_ab(h1)
            anp *= fa
        if _FLUC:
            if _PH is not None:
                dph, r1, r2, an, dpr_, rph_ = tdusrnpl(RNG, anp)
            else:
                dph, r1, r2, an = tdusrn(RNG, anp)
                dpr_ = rph_ = None
            if an == 0.0:
                if _POL:
                    _sprot(p, px0 + dpx * 0.5, p.pyr0 + dpy * 0.5, ppx, ppy, ppz, p.bsi, h1, pp, 0.0, cphi0, sphi0)
                return
            rho = al1 / theta
            uc = ph.cuc * h1**3 / ph.p0 / rho * fb
            if dpr_ is not None:
                # tradkf1: 光子ごとに tphrec (放射点 xr = x − r (px − dpx r/2) al、ds = al − r al)
                for i_ in range(int(an)):
                    dgp = dpr_[i_] * uc
                    xr = c[0] - rph_[i_] * (c[1] - 0.5 * dpx * rph_[i_]) * al
                    yr = c[2] - rph_[i_] * (c[3] - 0.5 * dpy * rph_[i_]) * al
                    _ph_rec(p, xr, c[1], yr, c[3], dgp, ppx, ppy, pp, h1, rho, al - rph_[i_] * al)
            de = dph * uc
            dg = -de / (1.0 + 2.0 * de)
            c[5] = max(_GMIN, c[5] + dg)
            ddpx = -r1 * dpx * dg
            ddpy = -r1 * dpy * dg
            c[0] += r2 * ddpx * al1
            c[2] += r2 * ddpy * al1
            ucde = uc * de
        else:
            uc = ph.cuc * h1**3 / ph.p0 * theta / al1 * fb
            de = _e.CAVE * anp * uc
            dg = -de / (1.0 + 2.0 * de)
            c[5] = max(_GMIN, c[5] + dg)
            ddpx = -0.5 * dpx * dg
            ddpy = -0.5 * dpy * dg
            c[0] += ddpx * al1 / 3.0
            c[2] += ddpy * al1 / 3.0
            ucde = uc * de
        c[1] += ddpx
        c[3] += ddpy
        pr = 1.0 + c[5]
        p2 = ph.p0 * pr
        h2 = _e._p2h(p2)
        p.dv = -c[5] * (1.0 + pr) / h2 / (h2 + p2) + ph.dvfs       # SAD のまま (h2+p2)。+dvfs (tradkn)
        c[4] = c[4] * p2 / h2 * h1 / pp
        if _POL:
            _sprot(p, px0 + dpx * 0.5, p.pyr0 + dpy * 0.5, ppx, ppy, ppz, p.bsi, h1, pp, ucde, cphi0, sphi0)


# ---------------------------------------------------------------------------
# ドリフト
# ---------------------------------------------------------------------------

def _tdrift_free(p: _Particle, al):
    """tdrift.f tdrift_free"""
    c = p.c
    dpz = _pxy2dpz(c[1], c[3])
    al1 = al / (1.0 + dpz)
    c[0] += c[1] * al1
    c[2] += c[3] * al1
    c[4] += dpz * al1 - p.dv * al


# ---------------------------------------------------------------------------
# 四極 (tquad.f tquad, ttfrin nord=4, tsolqu.f tsolqu / tsolqur)
# ---------------------------------------------------------------------------

def _ttfrin4(p: _Particle, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    c = p.c
    aki = ak / al / 4.0 / (1.0 + c[5])
    xi, yi = c[0], c[2]
    a = aki * xi**2
    b = aki * yi**2
    ab = a - b
    t = ab**2 / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    c[0] = xi + dx1
    c[2] = yi + dy1
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h**2
    px1 = (c[1] * yy + c[3] * h) / f
    py1 = (c[3] * xx - c[1] * h) / f
    c[4] = c[4] - px1 * (dx1 + t * xi) - py1 * (dy1 + t * yi)
    c[1] = px1
    c[3] = py1


def _tzsetparam0(dp, aln0, akk):
    """tsolque.f tzsetparam0: (w1, s1, xs1, dc1, sh2, xsh2, dch2)"""
    pr = 1.0 + dp
    akkp = akk / pr
    w1 = math.sqrt(akkp)
    phi1 = min(100.0, max(-100.0, aln0 * w1))
    s1, xs1, c1, dc1 = _xsincos(phi1)
    sh2, xsh2, ch2, dch2 = _xsincosh(phi1)
    return w1, s1, xs1, dc1, sh2, xsh2, dch2


_qb_drift, _qb_body, _qb_nonrad = _qb.make_track(_sqrt1, _tzsetparam0)     # 四極の本体 (B226): K64 の Strang 分割 (k64)・4 次の合成 (v2)


def _tsolqu(p: _Particle, al, ak, eps0=0.0, ak0x=0.0, ak0y=0.0, qum=False, ibsi=-1):
    """tsolqu.f tsolqu / tsolqum (bz=0, 放射なし)。ak = |K1|L > 0。ak0x/ak0y は MULT の K0 (dx0 = ak0x/ak)。
    ak = 0 なら tdrift (ak0 付きの厳密ドリフト)。
    qum=True は tsolqum (tmulti の放射なしで使う版): z の dv の項が tsolqu と違い
    −0.25·(ap − dv)·aln0 (tsolqu は −0.25·ap·aln0 − dv·aln0)。SAD のまま (Level 11 で確認)。
    本体の積分は quadbody (k64: K64 の Strang 分割、v2: 4 次の合成。docs/30 B226)"""
    if ak == 0.0:
        _tdrift_kick(p, al, ak0x, ak0y, False)
        return
    akk = ak / al
    c = p.c
    x, px, y, py, z = c[0], c[1], c[2], c[3], c[4]
    if _POL:
        dx0 = ak0x / ak
        dy0 = ak0y / ak
        if ibsi == 1:
            p.bsi = akk * (x + dx0) * (y + dy0)
        elif ibsi >= 0:
            p.bsi = 0.0
    x, px, y, py, z = _qb_nonrad(x, px, y, py, z, c[5], p.dv, al, ak, eps0, ak0x, ak0y, qum, _qb.scheme())
    c[0], c[1], c[2], c[3], c[4] = x, px, y, py, z
    if _POL and ibsi == 2:
        p.bsi = -akk * (x + ak0x / ak) * (y + ak0y / ak)


def _tdrift_kick(p: _Particle, al, ak0x, ak0y, enarad, bz=0.0):
    """tdrift.f tdrift (ak0x/ak0y ≠ 0): 一定の双極キック付きの厳密ドリフト (tsolconv)。bz≠0 は SOL 領域内 (BEND の直線は tdrift(al, bzs, phiy, phix))。
    内部座標 px/(1+δ)。放射は tradk1 (基準 px0, py0, z0 = 入口)。"""
    if ak0x == 0.0 and ak0y == 0.0:
        if bz == 0.0:
            _tdrift_free(p, al)
        else:
            _tdrift_solenoid(p, al, bz, enarad)
        return
    c = p.c
    if enarad and _POL:
        p.bsi = 0.0
    b = math.sqrt(ak0x**2 + ak0y**2 + (bz * al)**2)
    phix = ak0y / b
    phiy = ak0x / b
    phiz = bz * al / b
    dphizsq = phix**2 + phiy**2
    pr = 1.0 + c[5]
    alb = al * pr / b
    bsi0 = bz * al + ak0x * c[2] + ak0y * c[0]
    bzp = bz / pr
    if bz != 0.0:
        c[1] = c[1] + bzp * c[2] * 0.5
        c[3] = c[3] - bzp * c[0] * 0.5
    px0, py0, z0 = c[1], c[3], c[4]
    dpz0 = _pxy2dpz(c[1], c[3])
    pz0 = 1.0 + dpz0
    dpl = c[1] * phix + c[3] * phiy + dpz0 * phiz
    pl = phiz + dpl
    plx, ply, plz = pl * phix, pl * phiy, pl * phiz
    ptx = c[1] - plx
    pty = c[3] - ply
    ptz = dpz0 - dpl * phiz + dphizsq
    pbx = pty * phiz - ptz * phiy
    pby = ptz * phix - ptx * phiz
    pbz = ptx * phiy - pty * phix
    bpr = b / pr
    db = bpr - pbz
    ap = math.hypot(pz0, pbz)
    dphi = -math.atan2(pbz, pz0)
    phi = _lat._asinz(db / ap) - dphi
    phi, sinphi, xsinphi, cosphi, dcosphi = _lat._tsolconv(pz0, plz, pbz, bpr, phi)
    c[0] += (plx * phi + ptx * sinphi - pbx * dcosphi) * alb
    c[2] += (ply * phi + pty * sinphi - pby * dcosphi) * alb
    c[4] += ((dpl * phiz - dphizsq) * xsinphi + dpz0 * sinphi - pbz * dcosphi) * alb - p.dv * al
    px1 = px0 + ptx * dcosphi + pbx * sinphi
    py1 = py0 + pty * dcosphi + pby * sinphi
    if enarad:
        if _POL:
            p.bsi += bsi0 - ak0x * c[2] - ak0y * c[0]          # tdrift: bsi += bsi0 − ak0x·y − ak0y·x (bz=0 なので bsi0 = ak0x·y0 + ak0y·x0)
        c[1], c[3] = px1, py1
        p.pxr0, p.pyr0, p.zr0 = px0, py0, z0
        if _PH is not None:
            _ph_set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, al)          # tdrift.f: tsetpcvt (D1)
        _tradk_core(p, al, 0.0)
        px1, py1 = c[1], c[3]
    c[1] = px1 - bzp * c[2] * 0.5
    c[3] = py1 + bzp * c[0] * 0.5


_TRPT = False        # TRPT フラグ (tffs.f ndivrad: TRPT のとき放射の分割が 10 倍。D1: track6 の経路でも効かせる)


def _ndivrad(phi, ak, bz, eps, h0, trpt=None):
    """tffs.f ndivrad"""
    eps1 = 1.0 if eps <= 0.0 else eps
    aka = (abs(phi) + math.hypot(ak, bz) * 0.02) / eps1
    if trpt is None:
        trpt = _TRPT
    n = int(1.0 + h0 / 100.0 * aka * 10.0) if trpt else int(1.0 + h0 / 100.0 * aka)
    if _POL and _compat.spin_bpar():
        # v2 (B165): スピン追跡は 1 刻みの偏向角が大きいと、座標系の回転とスピン歳差の分け方の 1 次の誤差 (刻みの角に比例) が残る
        # (0.5 GeV の 15° 扇形を 3 刻みで dR/dpy の R12・R21 が 1e-2。刻み 0.087 → 0.017 → 0.0035 rad で 1.2e-2 → 2.7e-3 → 5.4e-4)。偏向角 5e-3 rad 以下にする
        n = max(n, 1 + int(abs(phi) / 5.0e-3))
    return n


def _tsolqur_ens(parts, al, ak, eps0=0.0):
    """tsolqu.f tsolqur (bz=0, ak0x=ak0y=0): スライス外側・粒子内側で、各スライスの後に全粒子 tradk
    (SAD と同じ順で乱数を引くため)。本体の積分は quadbody (docs/30 B226)"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    sch = _qb.scheme()
    ndiv = _qb.ndiv(abs(al) * ak, eps, sch)
    ndiv = min(1000, max(ndiv, _ndivrad(0.0, ak, 0.0, eps0, parts[0].ph.h0)))
    aln = al / ndiv
    akk = ak / al
    bpar_on = _POL and _compat.spin_bpar()          # v2 (B165): 本体の B∥ の積分 ∫B·dr = Δψ (ψ = k x y) を刻みごとに足す
    m = _qb.QB_M[sch]
    psi0 = [0.0] * len(parts)
    for g in range(m * ndiv):
        n, j = divmod(g, m)
        dl, bl, jb, _rad, _fin = _qb.group(g, aln, ndiv, sch)
        for ip, p in enumerate(parts):
            c = p.c
            par = _tzsetparam0(c[5], bl, akk)
            if j == 0:
                if _POL:
                    p.bsi = akk * c[0] * c[2] if n == 0 else 0.0
                if bpar_on:
                    psi0[ip] = akk * c[0] * c[2]
            c[0], c[2], c[4] = _qb_drift(c[0], c[1], c[2], c[3], c[4], dl)
            c[0], c[1], c[2], c[3], c[4] = _qb_body(c[0], c[1], c[2], c[3], c[4], p.dv, par, bl, 0.0, 0.0, False)
            if j == m - 1:
                if _POL and n == ndiv - 1:
                    # v2: ndiv = 1 (弱い四極) では入口の +ψ を出口の −ψ が上書きして落ちる (K64 は =)。v2 は足す
                    p.bsi = (p.bsi if (ndiv == 1 and _compat.v2()) else 0.0) - akk * c[0] * c[2]
                if bpar_on:
                    p.bpadd = (akk * c[0] * c[2] - psi0[ip]) / (1.0 + c[5])
        if j == m - 1:
            for p in parts:
                _tradk(p, aln, 0.0)
            if _PH is not None and _PH.get('ctx') is not None:
                _PH['ctx']['fr0'] += aln / _PH['ctx']['al']      # tsolqu.f: pcvt%fr0 += alr/al
    dl, _bl, _jb, _rad, _fin = _qb.group(m * ndiv, aln, ndiv, sch)
    for p in parts:
        c = p.c
        c[0], c[2], c[4] = _qb_drift(c[0], c[1], c[2], c[3], c[4], dl)


def _tsolqur(p: _Particle, al, ak, eps0=0.0):
    _tsolqur_ens([p], al, ak, eps0)


def _rot(c, ct, st):
    """tsolrot (dx=dy=chi1=chi2=0) の回転: x' = c x − s y。SAD は chi3 /= 0 のときだけ回す
    (無条件に掛けると x=inf, y=inf の発散粒子で 0*inf = NaN になり SAD と変わる。track6n.rot1 と同じ)"""
    if st == 0.0 and ct == 1.0:
        return
    x, px, y, py = c[0], c[1], c[2], c[3]
    c[0] = ct * x - st * y
    c[2] = st * x + ct * y
    c[1] = ct * px - st * py
    c[3] = st * px + ct * py


def _tsetfringep(el, akk):
    """tffs.f tsetfringep (順方向要素): (f1in, f2in, f1out, f2out) の係数。
    k64 (K64 のまま): f1 = −|akk·F1²|/24 (F1 の符号を捨てる。F1 < 0 でも f1 ≤ 0)。
    v2 (B323): f1 = −|akk|·F1·|F1|/24 (CALC の tqlfre と同じ。F1 の符号が a の符号を決める)。
    f2 = |akk|·F2 (F1 = F1+F1K1F/B, F2 = F2+F2K1F/B)。逆向き要素は未対応。"""
    f1in = _g(el, 'F1') + _g(el, 'F1K1F')
    f1out = _g(el, 'F1') + _g(el, 'F1K1B')
    f2in = _g(el, 'F2') + _g(el, 'F2K1F')
    f2out = _g(el, 'F2') + _g(el, 'F2K1B')
    if _compat.v2():
        # copysign(|akk·F1²|, F1) = |akk|·F1·|F1| (F1 ≥ 0 では k64 の式とビット同一)
        return (-math.copysign(abs(akk * f1in * f1in), f1in) / 24.0, abs(akk) * f2in,
                -math.copysign(abs(akk * f1out * f1out), f1out) / 24.0, abs(akk) * f2out)
    return (-abs(akk * f1in * f1in) / 24.0, abs(akk) * f2in,
            -abs(akk * f1out * f1out) / 24.0, abs(akk) * f2out)


def _tqlfr(p: _Particle, a1, b1):
    """tquad.f の線形フリンジ (mfring=1,3 入口 / 2,3 出口)。a1 = ±f1/p の f1, b1 = f2"""
    c = p.c
    pr = 1.0 + c[5]
    a = a1 / pr
    ea = math.exp(a)
    b = b1 / pr
    pxf = c[1] / ea
    pyf = c[3] * ea
    c[4] = (c[4] - (a * c[0] + b * (1.0 + 0.5 * a) * pxf) * c[1]
            + (a * c[2] + b * (1.0 - 0.5 * a) * pyf) * c[3])
    c[0] = ea * c[0] + b * c[1]
    c[2] = c[2] / ea - b * c[3]
    c[1] = pxf
    c[3] = pyf


def _spin_rot_z(p, chi3, ent):
    """tsolrot の calpol (chi1=chi2=0): スピンを z 軸まわりに chi3 回す (入口)、逆に戻す (出口)。trot33"""
    if _POL and chi3 != 0.0:
        c3, s3 = math.cos(chi3), math.sin(chi3)
        sx, sy = p.s[0], p.s[1]
        if ent:
            p.s[0] = c3 * sx - s3 * sy
            p.s[1] = s3 * sx + c3 * sy
        else:
            p.s[0] = c3 * sx + s3 * sy
            p.s[1] = -s3 * sx + c3 * sy


def _tquad_ens(parts, el, krad):
    """tquad.f tquad (bz=0)。全粒子をまとめて (放射ありは tsolqur の順で)。
    ROTATE は tsolrot の chi3 (theta2)、F1/F2/FRINGE は tsetfringep の係数。"""
    L = el.L
    K1L = _g(el, 'K1') * L
    if L == 0.0:
        # tquad.f: al=0 → tthin(nord=4, ak=K1 (積分値), dx, dy, theta (ROT。akang なし), krad=.false., fringe=.false.)
        ak = _g(el, 'K1')
        if ak == 0.0:
            return
        from types import SimpleNamespace
        thin = SimpleNamespace(L=0.0, DX=_g(el, 'DX'), DY=_g(el, 'DY'), ROTATION=_g(el, 'ROTATION'), DISFRIN=1.0)
        for p in parts:
            _tthin_gen(p, thin, False, ak, 1)
        return
    if K1L == 0.0:
        for p in parts:
            _tdrift_free(p, L)
        return
    dx, dy = _g(el, 'DX'), _g(el, 'DY')
    ak = abs(K1L)
    th2 = _g(el, 'ROTATION') + (0.5 * math.pi if K1L < 0 else 0.0)   # akang
    fringe = _g(el, 'DISFRIN') == 0.0
    mfring = int(round(_g(el, 'FRINGE')))
    f1in, f2in, f1out, f2out = _tsetfringep(el, K1L / L)
    ct, st = math.cos(th2), math.sin(th2)
    for p in parts:
        c = p.c
        c[0] -= dx                       # tsolrot (ent): ずらす → 回す
        c[2] -= dy
        _rot(c, ct, st)
        _spin_rot_z(p, th2, True)
        if krad:
            p.setr0()
            if _POL:
                p.bsi = 0.0
        if fringe and mfring > -4 and mfring != 2:
            _ttfrin4(p, ak, L)
        if mfring in (1, 3):
            _tqlfr(p, f1in, f2in)
    if krad:
        if _PH is not None:
            _ph_set(dx, dy, th2, 0.0, 0.0, 0.0, L, -0.5 * -abs(f1in) / L)     # tquad.f: tsetpcvt(l_track, dx, dy, theta2, 0, 0, 0, al)、fr0 = −f1in/(2 al) (放射の引数は −|f1in|: B323 の前と同じ)
        if f1in != 0.0:
            for p in parts:
                _tradk(p, -abs(f1in), 0.0)     # SAD のまま: 長さの引数に係数 f1in を渡す (K64 の癖。符号は B323 の前 (f1 ≤ 0) のまま −|f1in|)
        for p in parts:
            p.setr0()
        _tsolqur_ens(parts, L, ak)
    else:
        for p in parts:
            _tsolqu(p, L, ak)
    for p in parts:
        if mfring in (2, 3):
            _tqlfr(p, -f1out, f2out)
        if fringe and mfring > -4 and mfring != 1:
            _ttfrin4(p, -ak, L)
    if krad and f1out != 0.0:
        if _PH is not None and _PH.get('ctx') is not None:
            _PH['ctx']['fr0'] = 1.0 - 0.5 * -abs(f1out) / L           # tquad.f:108 (放射の引数は −|f1out|)
        for p in parts:
            _tradk(p, -abs(f1out), 0.0)
    for p in parts:
        c = p.c
        _rot(c, ct, -st)                 # tsolrot (exit): 回し戻す → ずらし戻す
        _spin_rot_z(p, th2, False)
        c[0] += dx
        c[2] += dy


def _tquad(p: _Particle, el, krad):
    _tquad_ens([p], el, krad)


# ---------------------------------------------------------------------------
# 6 極 (tquad.f tthin, ttfrin nord=6)
# ---------------------------------------------------------------------------

def _tthin_gen(p: _Particle, el, krad, ak, kord):
    """tquad.f tthin の kord≠2 (OCT 以上) の分岐 (DX/DY/θ は TENT/TEXIT)。フリンジは ttfrin の一般次数 (nord=2·kord+2)"""
    L = el.L
    if ak == 0.0:
        _tdrift_free(p, L)
        return
    c = p.c
    dx, dy, theta = _g(el, 'DX'), _g(el, 'DY'), _g(el, 'ROTATION')
    cost, sint = _tent(c, dx, dy, theta)
    if krad:
        p.setr0()
        if _PH is not None:
            _ph_set(dx, dy, theta, 0.0, 0.0, 0.0, L)          # tquad.f tthin: tsetpcvt (D1: SEXT・OCT… の光子)
    fringe = _g(el, 'DISFRIN') == 0.0
    if fringe:
        _ttfrin_n(p, 2 * kord + 2, ak, L)
    akf = ak / math.factorial(kord)
    ala = L * _TT_ALPHA1
    alb = L * _TT_ALPHA
    if L != 0.0:
        _tdrift_free(p, ala)
        akf *= 0.5

    def kick():
        aki = akf / (1.0 + c[5])
        cx = complex(c[0], -c[2])**kord
        c[1] -= aki * cx.real
        c[3] -= aki * cx.imag
    kick()
    if L != 0.0:
        if kord == 1:
            f1 = _TT_BETA * L
            f3 = _TT_GAMMA * 2.0 * L**2
            f4 = 0.5 * f1
            f5 = f3
        else:
            f1 = _TT_BETA * kord * L
            f2 = 2.0 * _TT_GAMMA * kord * (kord - 1) * L**2
            f3 = _TT_GAMMA * kord * (kord + 1) * L**2
            f4 = 0.5 * _TT_BETA * L
            f5 = 2.0 * _TT_GAMMA * kord * L**2
        dpz = _pxy2dpz(c[1], c[3])
        al1 = alb / (1.0 + dpz)
        xi = c[0] + c[1] * al1
        yi = c[2] + c[3] * al1
        zi = c[4] + dpz * al1 - p.dv * alb
        aki = akf / (1.0 + c[5]) * 2.0
        r = xi**2 + yi**2
        if kord == 1:
            rcx1 = (xi - yi) * (xi + yi)
            c[1] += aki**2 * (f1 - aki * f3) * xi
            c[3] += aki**2 * (f1 + aki * f3) * yi
            zi += aki**2 * (f4 * r - f5 * aki * rcx1)
        else:
            rk1 = r**(kord - 2)
            rk = r * rk1
            cx = complex(xi, -yi)**kord
            rcx1 = xi * cx.real + yi * cx.imag
            c[1] += aki**2 * (f1 * rk * xi - aki * (f2 * rk1 * xi * rcx1 + f3 * rk * cx.real))
            c[3] += aki**2 * (f1 * rk * yi - aki * (f2 * rk1 * yi * rcx1 + f3 * rk * cx.imag))
            zi += aki**2 * rk * (f4 * r - f5 * aki * rcx1)
        dpz = _pxy2dpz(c[1], c[3])
        al1 = alb / (1.0 + dpz)
        c[0] = xi + c[1] * al1
        c[2] = yi + c[3] * al1
        c[4] = zi + dpz * al1 - p.dv * alb
        kick()
        _tdrift_free(p, ala)
    if fringe:
        _ttfrin_n(p, 2 * kord + 2, -ak, L)
    if krad:
        _tradk(p, L, 0.0)
    _texit(c, dx, dy, theta, cost, sint)


def _ttfrin6(p: _Particle, ak, al):
    if al == 0.0 or ak == 0.0:
        return
    c = p.c
    aki = ak / al / 24.0 / (1.0 + c[5])
    cx = complex(c[0], c[2])
    cp = complex(c[1], -c[3])
    cz1 = cx * cx
    cz = cz1 * cx
    a = aki * cz.imag * 2.0
    ca = -aki * cx * (cz / 4.0 - cz.conjugate())
    cx1 = cx + ca
    d = 1.0 + a**2 - (aki * 3.0)**2 * (cz.real**2 + cz.imag**2)
    cp1 = (complex(1.0, a) * cp - (aki * 3.0) * cz1 * (cx * cp).conjugate()) / d
    c[4] = c[4] - (ca * cp1).real
    c[0] = cx1.real
    c[2] = cx1.imag
    c[1] = cp1.real
    c[3] = -cp1.imag


def _tent(c, dx, dy, theta):
    """inc/TENT.inc: ずらす → 回す (|cos|≥|sin| では y を更新後の x から求める SAD の書き方のまま)"""
    if theta != 0.0:
        cost, sint = math.cos(theta), math.sin(theta)
    else:
        cost, sint = 1.0, 0.0
    if dx != 0.0:
        c[0] -= dx
    if dy != 0.0:
        c[2] -= dy
    if theta != 0.0:
        if abs(cost) >= abs(sint):
            c[0] = cost * c[0] - sint * c[2]
            c[2] = (sint * c[0] + c[2]) / cost
            c[1] = cost * c[1] - sint * c[3]
            c[3] = (sint * c[1] + c[3]) / cost
        else:
            xi = c[0]
            c[0] = xi * cost - c[2] * sint
            c[2] = xi * sint + c[2] * cost
            pxi = c[1]
            c[1] = pxi * cost - c[3] * sint
            c[3] = pxi * sint + c[3] * cost
    return cost, sint


def _texit(c, dx, dy, theta, cost, sint):
    """inc/TEXIT.inc: 回し戻す → ずらし戻す"""
    if theta != 0.0:
        if abs(cost) >= abs(sint):
            c[0] = cost * c[0] + sint * c[2]
            c[2] = (-sint * c[0] + c[2]) / cost
            c[1] = cost * c[1] + sint * c[3]
            c[3] = (-sint * c[1] + c[3]) / cost
        else:
            xi = c[0]
            c[0] = xi * cost + c[2] * sint
            c[2] = -xi * sint + c[2] * cost
            pxi = c[1]
            c[1] = pxi * cost + c[3] * sint
            c[3] = -pxi * sint + c[3] * cost
    if dy != 0.0:
        c[2] += dy
    if dx != 0.0:
        c[0] += dx


_TT_ALPHA = 2.88675134594812882e-1
_TT_BETA = 2.23290993692602255e-2
_TT_GAMMA = 9.43738783765593145e-4
_TT_ALPHA1 = 0.5 - _TT_ALPHA


_THIN_KEY = {'Sextupole': ('K2', 2), 'Octupole': ('K3', 3), 'Decapole': ('K4', 4), 'Dodecapole': ('K5', 5)}


def _tthin_sext(p: _Particle, el, krad):
    """tquad.f tthin (SEXT: nord=6, kord=2; OCT: kord=3; DECA: kord=4; DODECA: kord=5)"""
    L = el.L
    kn, kord = _THIN_KEY[type(el).__name__]
    K = _g(el, kn)
    ak = K * L if L != 0.0 else K
    if kord != 2:
        _tthin_gen(p, el, krad, ak, kord)
        return
    if ak == 0.0:
        _tdrift_free(p, L)
        return
    c = p.c
    dx, dy, theta = _g(el, 'DX'), _g(el, 'DY'), _g(el, 'ROTATION')
    cost, sint = _tent(c, dx, dy, theta)
    if krad:
        p.setr0()
        if _PH is not None:
            _ph_set(dx, dy, theta, 0.0, 0.0, 0.0, L)          # tquad.f tthin: tsetpcvt (D1: SEXT・OCT… の光子)
    fringe = _g(el, 'DISFRIN') == 0.0
    if fringe:
        _ttfrin6(p, ak, L)
    akf = ak / 2.0          # fact(kord=2)
    ala = L * _TT_ALPHA1
    alb = L * _TT_ALPHA
    if L != 0.0:
        _tdrift_free(p, ala)
        akf *= 0.5

    def kick():
        aki = akf / (1.0 + c[5])
        x, y = c[0], c[2]
        c[1] -= aki * (x - y) * (x + y)
        c[3] += 2.0 * aki * x * y
    kick()
    if L != 0.0:
        f1 = _TT_BETA * 2.0 * L
        f2 = 4.0 * _TT_GAMMA * L**2
        f3 = _TT_GAMMA * 6.0 * L**2
        f4 = 0.5 * _TT_BETA * L
        f5 = f2
        dpz = _pxy2dpz(c[1], c[3])
        al1 = alb / (1.0 + dpz)
        xi = c[0] + c[1] * al1
        yi = c[2] + c[3] * al1
        zi = c[4] + dpz * al1 - p.dv * alb
        aki = akf / (1.0 + c[5]) * 2.0
        r = xi**2 + yi**2
        cx = complex(xi, -yi)**2
        rcx1 = xi * cx.real + yi * cx.imag
        c[1] += aki**2 * (f1 * r * xi - aki * (f2 * xi * rcx1 + f3 * r * cx.real))
        c[3] += aki**2 * (f1 * r * yi - aki * (f2 * yi * rcx1 + f3 * r * cx.imag))
        zi += aki**2 * r * (f4 * r - f5 * aki * rcx1)
        dpz = _pxy2dpz(c[1], c[3])
        al1 = alb / (1.0 + dpz)
        c[0] = xi + c[1] * al1
        c[2] = yi + c[3] * al1
        c[4] = zi + dpz * al1 - p.dv * alb
        kick()
        _tdrift_free(p, ala)
    if fringe:
        _ttfrin6(p, -ak, L)
    if krad:
        _tradk(p, L, 0.0)
    _texit(c, dx, dy, theta, cost, sint)


# ---------------------------------------------------------------------------
# BEND (tbend.f tbend, tbendr, tbendcore)
# ---------------------------------------------------------------------------

def _tbendcore(p: _Particle, al, phi0, cosp1, sinp1, cosp2, sinp2, mfring, fringe,
               cosw, sinw, sqwh, sinwp1, krad, alr, bp, bsi1=0.0, bsi2=0.0):
    """tbend.f tbendcore。bp = (rhob, rho0, drhob, fb1, fb2) (tbendcom のモジュール変数)。
    bsi1, bsi2: POL の bsi の入口・出口の寄与 (最初・最後のスライスで 1)"""
    rhob, rho0, drhob, fb1, fb2 = bp
    if (mfring > 0 or mfring == -1) and fb1 != 0.0:
        dxfr1 = fb1**2 / rhob / 24.0
        dyfr1 = fb1 / rhob**2 / 6.0
        dzfr1 = dxfr1 * sinp1
        dyfra1 = 4.0 * dyfr1 / fb1**2 if fringe else 0.0
    else:
        dxfr1 = dyfr1 = dyfra1 = dzfr1 = 0.0
    if (mfring > 0 or mfring == -2) and fb2 != 0.0:
        dxfr2 = fb2**2 / rhob / 24.0
        dyfr2 = fb2 / rhob**2 / 6.0
        dzfr2 = dxfr2 * sinp2
        dyfra2 = 4.0 * dyfr2 / fb2**2 if fringe else 0.0
    else:
        dxfr2 = dyfr2 = dyfra2 = dzfr2 = 0.0
    if cosp1 * cosp2 > 0.0:
        dcosp = (sinp2 - sinp1) * (sinp2 + sinp1) / (cosp1 + cosp2)
    else:
        dcosp = cosp1 - cosp2
    sinp2wp1 = sinp2 + sinwp1
    sinp2p1 = sinp2 + sinp1
    cosp1p2 = cosp1 * cosp2
    smin = 1e-4
    c = p.c
    if _POL:
        p.bsi += bsi1 * c[2] / rhob
    y_in = c[2]
    dp = c[5]
    pp = 1.0 + dp
    pinv = 1.0 / pp
    rhoe = rhob * pp
    pxi = c[1]
    pyi = c[3]
    dpzi = _pxy2dpz(pxi, pyi)
    pzi = 1.0 + dpzi
    sp1 = sinp1 / pzi
    x1 = c[0] / (cosp1 - pxi * sp1)
    dz1 = x1 * sp1
    y1 = c[2] + pyi * dz1
    px1 = pxi * cosp1 + pzi * sinp1
    x1 = x1 + dxfr1 * dp * pinv
    py1 = pyi + (dyfr1 - dyfra1 * y1**2) * y1 * pinv**2
    z1 = c[4] - dz1 + (dxfr1 * px1 + (0.5 * dyfr1 - 0.25 * dyfra1 * y1**2) * y1**2 * pinv) * pinv - dzfr1
    pv1sqi = 1.0 / max(smin, 1.0 - px1**2)
    fa = y1 / rhoe * math.sqrt(pv1sqi)
    f = (1.0 - (y1 / rhob)**2 / 6.0) * fa
    ff = 0.25 * (f + fa) * y1 * pv1sqi
    x2 = x1 + ff
    py2 = py1 - px1 * f
    z2 = z1 - px1 * ff
    dph2 = _sqrt1(-py2**2)
    ph2 = 1.0 + dph2
    dpx2 = pxi * cosp1 + (dpzi - dph2) * sinp1
    pz2 = 1.0 + _pxy2dpz(px1, py2)
    drho1 = drhob + rhob * dp
    drho = drho1 + rhoe * dph2
    t2 = (px1 + ph2 * sinp1) / (pz2 + ph2 * cosp1)
    dpx3a = (x2 * sinw - drho * sinp2wp1) / rhoe
    dpx3 = dpx3a - dpx2 * (cosw - sinw * t2)
    px3 = ph2 * sinp2 + dpx3
    sqwht2 = sqwh + sinw * t2
    px1px3 = ph2 * sinp2p1 + dpx3a + dpx2 * sqwht2
    dpz3 = _pxy2dpz(px3, py2)
    pz3 = 1.0 + dpz3
    dpz32 = px1px3 * (px1 - px3) / (pz2 + pz3)
    pz3ph2 = pz3 + ph2 * cosp2
    pz2px1 = pz2 * cosp1 + px1 * sinp1
    t3 = (px3 + ph2 * sinp2) / pz3ph2
    t2t3 = ((ph2 * sinp2 + px1px3) / pz3ph2
            + (ph2 * sinp1 + px1 * (dpz32 - ph2 * dcosp) / pz3ph2) / (pz2 + ph2 * cosp1))
    t4 = (cosp2 + t3 * sinp2) * pz2px1
    x3 = (x2 * (cosw - rho0 / rhoe * t3 * sinw / ph2)
          + (rho0 * (cosw * t2t3 + sinw * (1.0 - t2 * t3)) * dpx2
             - drho * (-sinp2wp1 * rho0 / rhoe * t3 - dpz32 - sqwh * pz2 - sinw * px1)) / ph2)
    da = _e._asinz(
        (dpx2 * (sinp1 * (t2 * (pz3 * cosp2 + px3 * sinp2) - px1 * cosp2)
                 - sinp2 * (t3 * pz2px1 - px3 * cosp1)
                 + cosp1p2 * dpz32
                 + t4 * sqwht2)
         + dpx3a * t4) / ph2**2)
    phi0a = phi0 + da
    y3 = y1 + py2 * rhoe * phi0a
    z3 = z2 - phi0 * drho1 - da * rhoe - p.dv * al
    pv2sqi = 1.0 / max(smin, 1.0 - px3**2)
    fa = y3 / rhoe * math.sqrt(pv2sqi)
    f = (1.0 - (y3 / rhob)**2 / 6.0) * fa
    ff = 0.25 * (f + fa) * y3 * pv2sqi
    py4 = py2 - px3 * f
    z4 = z3 - px3 * ff
    x4 = x3 - ff - dxfr2 * dp * pinv
    py4 = py4 + (dyfr2 - dyfra2 * y3**2) * y3 * pinv**2
    z4 = z4 + (dxfr2 * px3 + (0.5 * dyfr2 - 0.25 * dyfra2 * y3**2) * y3**2 * pinv) * pinv - dzfr2
    dpz4 = _pxy2dpz(px3, py4)
    pxo = -cosp2 * dpx3 + sinp2 * (dpz4 - dpz3 - dpx3 * t3)
    dz4 = x4 * sinp2 / (cosp2 * (1.0 + dpz4) + sinp2 * px3)
    c[0] = x4 * cosp2 + pxo * dz4
    c[1] = pxo
    c[2] = y3 + py4 * dz4
    c[3] = py4
    c[4] = z4 - dz4
    if _POL and _compat.spin_bpar():
        p.bpadd = _BPSIGN * (c[2] - y_in) / rhob / (1.0 + c[5])      # v2 (B165): 本体の B∥ ∫B·dr = (y_出 − y_入)/ρ
    if _POL:
        p.bsi -= bsi2 * c[2] / rhob
    if krad:
        if sinp2 != 0.0 and _compat.v2():
            p.zr0 -= _zr0_edge(p, sinp2 / cosp2, 1.0)      # v2: 出口の端の角で場の中の行路が短い (B220: 面までの直進の行程を厳密に。tradk の al1 = al - z + zr0 に入れる)
        _tradk(p, alr, phi0)


def _spin_rot_cs(p, cost, sint):
    """v2 (P1): 座標の ROTATE (回転 θ: x' = cos θ x − sin θ y、y' = sin θ x + cos θ y) と同じ回転をスピンに掛ける (θ = cost, sint。出口は −θ)。
    K64 の tbshift・tent は座標だけ回してスピンを回さない (ROTATE ≠ 0 の BEND のスピンが回転後の系でなく元の系の鉛直軸のまわりに歳差する。
    厳密な Thomas-BMT の積分と Emittance[] (tchge の srot の回転) は回した系の軸のまわり。docs/history/v2next_P1.md)"""
    if _POL and (sint != 0.0 or cost != 1.0) and _compat.v2():
        sx, sy = p.s[0], p.s[1]
        p.s[0] = cost * sx - sint * sy
        p.s[1] = sint * sx + cost * sy


def _tbshift(p: _Particle, dx, dy, phi0, cost, sint, ent):
    """tbend.f tbshift"""
    c = p.c
    _spin_rot_cs(p, cost, sint)
    if dx != 0.0 or dy != 0.0:
        phih = phi0 * 0.5
        sph = math.sin(phih)
        dcph = -2.0 * math.sin(phih * 0.5)**2
        st1 = sint if ent else -sint
        dxa = dx * cost - dy * st1
        dx1 = dx + dcph * cost * dxa
        dy1 = dy - dcph * st1 * dxa
        ds = dxa * sph
        if ent:
            al = ds / (1.0 + _pxy2dpz(c[1], c[3]))
            x1 = c[0] + c[1] * al - dx1
            c[2] = c[2] + c[3] * al - dy1
            c[4] = c[4] - al
            c[0] = x1 * cost - c[2] * sint
            c[2] = x1 * sint + c[2] * cost
            px1 = c[1]
            c[1] = px1 * cost - c[3] * sint
            c[3] = px1 * sint + c[3] * cost
        else:
            x1 = c[0] * cost - c[2] * sint
            c[2] = c[0] * sint + c[2] * cost
            px1 = c[1]
            c[1] = px1 * cost - c[3] * sint
            c[3] = px1 * sint + c[3] * cost
            al = ds / (1.0 + _pxy2dpz(c[1], c[3]))
            c[0] = x1 + c[1] * al - dx1
            c[2] = c[2] + c[3] * al - dy1
            c[4] = c[4] - al
    elif sint != 0.0 or cost != 1.0:
        y1 = c[2]
        c[2] = c[0] * sint + y1 * cost
        c[0] = c[0] * cost - y1 * sint
        px1 = c[1]
        c[1] = px1 * cost - c[3] * sint
        c[3] = px1 * sint + c[3] * cost


def _tbrot(p: _Particle, alg, phig, dtheta, dchi2, ent):
    """tbend.f tbrot (BEND のロール誤差、alg/phig の点まわり)"""
    c0 = math.cos(phig)
    s0 = math.sin(phig)
    s2, xs2, c2, dc2 = _xsincos(dchi2)
    st, xst, ct, dct = _xsincos(dtheta)
    dr11 = dct * c0**2 + dc2 * s0**2 + c0 * s0 * s2 * st
    dr22 = c2 * dct + dc2
    dr33 = dc2 * c0**2 + dct * s0**2 - c0 * s0 * s2 * st
    r11, r22, r33 = dr11 + 1.0, dr22 + 1.0, dr33 + 1.0
    if ent:
        r12 = s0 * s2 - c0 * c2 * st
        r13 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r21 = -ct * s0 * s2 + c0 * st
        r23 = c0 * ct * s2 + s0 * st
        r31 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r32 = -c0 * s2 - c2 * s0 * st
    else:
        r21 = s0 * s2 - c0 * c2 * st
        r31 = c0 * ((dct - dc2) * s0 - c0 * s2 * st)
        r12 = -ct * s0 * s2 + c0 * st
        r32 = c0 * ct * s2 + s0 * st
        r13 = s0 * (c0 * (dct - dc2) + s0 * s2 * st)
        r23 = -c0 * s2 - c2 * s0 * st
    dl = alg if phig == 0.0 else alg * s0 / phig
    xl0 = dl * s0
    zi = -dl * c0
    dxz = r13 * zi
    dyz = r23 * zi
    dzz = dr33 * zi
    c = p.c
    xi = c[0] + xl0
    yi = c[2]
    pxi = c[1]
    pyi = c[3]
    pzi = 1.0 + _pxy2dpz(pxi, pyi)
    xf = dr11 * xi + r12 * yi + dxz + c[0]
    yf = r21 * xi + dr22 * yi + dyz + yi
    zf = r31 * xi + r32 * yi + dzz
    pxf = r11 * pxi + r12 * pyi + r13 * pzi
    pyf = r21 * pxi + r22 * pyi + r23 * pzi
    pzf = r31 * pxi + r32 * pyi + r33 * pzi
    c[1] = pxf
    c[3] = pyf
    c[0] = xf - pxf / pzf * zf
    c[2] = yf - pyf / pzf * zf
    c[4] = c[4] + zf / pzf


# --- tbendi (K1≠0 の BEND): tbend.f module bendib (tbendiinit, tbendicorr, tbendibody), tbendi.f

class _BendI:
    """bendib のモジュール変数 (粒子 1 個分)"""
    pass


def _tbendiinit(m, ak1, al):
    m.p = 1.0 + m.dp
    m.rhosq = m.rho0 * m.rhoe
    m.drhop = (m.rhoe - m.rho0) / m.rhosq
    m.akk = ak1 / al
    m.akxsq = -1.0 / m.rhosq - m.akk / m.p
    m.akysq = m.akk / m.p
    if m.akxsq > 0.0:
        akx = math.sqrt(m.akxsq)
        phix = akx * al
        m.dcx = 2.0 * math.sinh(0.5 * phix)**2
        spx = math.sinh(phix)
        xspx = phix - spx
        m.dcxkx = m.dcx / m.akxsq
        m.sxkx = spx / akx
        m.xsxkx = xspx / akx / m.akxsq
    elif m.akxsq < 0.0:
        akx = math.sqrt(-m.akxsq)
        phix = akx * al
        m.dcx = -2.0 * math.sin(0.5 * phix)**2
        spx = math.sin(phix)
        xspx = phix - spx
        m.dcxkx = m.dcx / m.akxsq
        m.sxkx = spx / akx
        m.xsxkx = (xspx / akx / m.akxsq) if _compat.v2() else (-xspx / akx / m.akxsq)       # B250 (v2 は連続な式)
    else:
        m.dcx = 0.0
        m.dcxkx = 0.5 * al**2
        m.sxkx = al
        m.xsxkx = -1.0 / 6.0 * al**3
    if m.akysq > 0.0:
        aky = math.sqrt(m.akysq)
        phiy = aky * al
        m.dcy = 2.0 * math.sinh(0.5 * phiy)**2
        spy = math.sinh(phiy)
        m.aksy = aky * spy
        m.syky = spy / aky
    elif m.akysq < 0.0:
        aky = math.sqrt(-m.akysq)
        phiy = aky * al
        m.dcy = -2.0 * math.sin(0.5 * phiy)**2
        spy = math.sin(phiy)
        m.aksy = -aky * spy
        m.syky = spy / aky
    else:
        m.dcy = 0.0
        m.aksy = 0.0
        m.syky = al


def _tbendicorr(m, ak1, al, phi):
    s = min(0.9999, m.pxi**2 + m.pyi**2)
    dpzi = _sqrt1(-s)
    pzi = 1.0 + dpzi
    m.zi = m.zi + (0.5 * s + (1.0 + m.xi / m.rho0) * dpzi / pzi) * al
    m.xi = (-al * dpzi * m.pxi + pzi * m.xi) / (pzi - phi * m.pxi)
    m.pxi = m.pxi + phi * dpzi
    m.yi = m.yi + (-al * dpzi + m.xi * phi) * m.pyi / pzi
    xr = m.xi / m.rho0
    if getattr(m, 'lin', False):
        # v2 (B160 のオプション compat.bend_k1_def() == 'lin'): B' 一様 = 'log' + 6 極 K2 = +hK1。kick のポテンシャル ak1 h (x³/3 − x y²/2) / p
        # (x の kick は ak1 h x² (級数なし)、y の項 +ak1 h y²/2 と py の kick +ak1 h x y が増える)
        hk = ak1 / m.rho0 / m.p
        m.pyi = m.pyi + hk * m.xi * m.yi
        m.pxi = m.pxi - hk * (m.xi * m.xi - 0.5 * m.yi * m.yi)
    else:
        m.pxi = m.pxi - ak1 * m.xi * xr * (0.5 - xr * (2.0 - xr) / 12.0) / m.p


def _tbendibody(m, al):
    xiksq = m.drhop + m.akxsq * m.xi
    dxf = m.dcxkx * xiksq + m.sxkx * m.pxi
    dpxf = m.sxkx * xiksq + m.dcx * m.pxi
    dyf = m.dcy * m.yi + m.syky * m.pyi
    dpyf = m.aksy * m.yi + m.dcy * m.pyi
    hi = (0.5 * (m.pxi**2 + m.pyi**2 - m.akxsq * m.xi**2 - m.akysq * m.yi**2)
          - m.drhop * m.xi)
    m.zi = (m.zi - 0.25 * (dxf * m.pxi + m.xi * dpxf + dxf * dpxf
                           + dyf * m.pyi + m.yi * dpyf + dyf * dpyf
                           - (5.0 / m.rho0 - 1.0 / m.rhoe)
                           * (m.drhop * m.xsxkx - m.sxkx * m.xi - m.dcxkx * m.pxi))
            - hi * al * 0.5)
    m.xi += dxf
    m.pxi += dpxf
    m.yi += dyf
    m.pyi += dpyf


def _tbendi(p: _Particle, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak, dx, dy, theta,
            dtheta, cost, sint, dchi2, alg, phig, fb1, fb2, mfring, enarad, fringe, eps0,
            first=True, last=True):
    """tbendi.f tbendi (K1≠0 の BEND の追跡)。ak = 積分 K1。"""
    _tbshift(p, dx, dy, phi0, cost, sint, True)
    if dtheta != 0.0 or dchi2 != 0.0:
        _tbrot(p, alg, phig, dtheta, dchi2, True)
    tanp1 = sinp1 / cosp1
    tanp2 = sinp2 / cosp2
    m = _BendI()
    m.rhob = al / phib
    m.rho0 = al / phi0
    _bpar_on = _compat.spin_bpar()
    m.lin = _compat.bend_k1_lin()                        # v2 のオプション: K1 の定義 'lin' (B' 一様)。k64 は常に False
    eps = 1e-3 if eps0 == 0.0 else 1e-3 * eps0          # epsbend
    ndiv = 1 + int(abs(phi0 / eps))
    krad = enarad and al != 0.0
    f1r = f2r = 0.0
    n1, n2 = 1, 0
    alc = al
    v2 = _compat.v2()
    if krad:
        p.setr0()
        if v2:
            p.zr0 -= _zr0_edge(p, tanp1, -1.0)            # v2: 入口の端の角で場の中の行路が短い (B220)
        if _POL:
            p.bsi = 0.0
        if first and fb1 != 0.0:
            f1r = 0.5 * fb1
            n1 = -1
        if last and fb2 != 0.0:
            f2r = 0.5 * fb2
            n2 = 2
        alc = al - f1r - f2r
        ndiv = max(ndiv, _ndivrad(phib * alc / al, ak * alc / al, 0.0, eps0, p.ph.h0))
    ndiv = min(1000, ndiv)
    n2 = ndiv + n2
    c = p.c
    if fringe and mfring > -4 and mfring != 2:
        _ttfrin4(p, ak, al)
    if fb1 != 0.0:
        dxfr1 = fb1**2 / m.rhob / 24.0
        dyfr1 = fb1 / m.rhob**2 / 6.0
        dyfra1 = 4.0 * dyfr1 / fb1**2
        dp = c[5]
        pr = 1.0 + dp
        c[0] = c[0] + dxfr1 * dp / pr
        c[3] = c[3] + (dyfr1 - dyfra1 * c[2]**2) * c[2] / pr**2
        c[4] = c[4] + (dxfr1 * c[1] + (0.5 * dyfr1 - 0.25 * dyfra1 * c[2]**2) * c[2]**2 / pr) / pr
    af = 1.0 if fringe else 0.0
    aln = alc / ndiv
    alx, alr, akxn, phixn = {}, {}, {}, {}
    for n in range(n1, n2 + 1):
        alx[n], alr[n] = _tbendal(n, ndiv, f1r, f2r, aln)
        akxn[n] = ak * alx[n] / al
        phixn[n] = phi0 * alx[n] / al
    m.dp = c[5]
    m.p = 1.0 + m.dp
    m.rhoe = m.rhob * m.p
    yi = c[2]
    v2edge = v2 and fringe          # v2 (v2next B1): sinp≠0 の端は CALC (tbende) と同じ厳密な tbedge (K64 の薄い端は x の 2 次で 1e-5 の差)
    if v2edge and sinp1 != 0.0:
        pe = 1.0 + c[5]            # 追跡の px・py は運動量で割ったもの。tbedge は正準運動量
        ce = _lat._tbedge_cod(np.array([c[0], c[1] * pe, c[2], c[3] * pe, c[4], c[5]], dtype=float), m.rhob, math.atan2(sinp1, cosp1), True)
        m.xi, m.pxi, m.yi, m.pyi, m.zi = ce[0], ce[1] / pe, ce[2], ce[3] / pe, ce[4]
        if _compat.bend_k1_edge() == 'face':
            kt = ak / al * tanp1 / pe              # v2 (P2 項目 5): 四極の磁場も端面から始まる模型 (半径方向の入口の面の楔形の領域の蹴り)
            m.pxi += kt * m.xi * m.xi
            m.pyi -= kt * m.xi * m.yi
    else:
        f = yi / m.rhoe
        fpx = af * c[1]
        m.pyi = c[3] - (tanp1 + fpx) * f
        ff = af * yi * f * 0.5
        m.xi = c[0] + ff
        m.zi = c[4] - ff * fpx
        m.pxi = c[1] + tanp1 * m.xi / m.rhoe
        m.yi = yi
    akx0 = alx0 = phix0 = 0.0
    if krad and _PH is not None:
        _ph_set(dx, dy, theta, dtheta, phi0, phig, al, 0.0)       # tbendi.f: 粒子ごとに tsetpcvt (fr0=0)
    for n in range(n1, n2 + 1):
        _tbendiinit(m, akxn[n], alx[n])
        _tbendicorr(m, (akxn[n] + akx0) * 0.5, (alx[n] + alx0) * 0.5, (phixn[n] + phix0) * 0.5)
        akx0, alx0, phix0 = akxn[n], alx[n], phixn[n]
        if krad and _POL and v2 and n == n1:
            # v2 (B58): 入口の端の縦磁場の積分 ∫Bs ds = y/ρ + k·x·y (K1=0 の tbend と同じ。K64 の tbendi は y/ρ がなく、k·x·y を全スライスで足す)
            p.bsi += m.yi / m.rhob + ak / al * m.xi * m.yi
        if _POL and _bpar_on:
            psi_a = m.yi / m.rhob + ak / al * m.xi * m.yi
        _tbendibody(m, alx0)
        if _POL and _bpar_on:
            # v2 (B165): 本体の B∥ の積分 ∫B·dr = ψ_終 − ψ_始 (ψ = y/ρ + k x y)
            p.bpadd = _BPSIGN * (m.yi / m.rhob + ak / al * m.xi * m.yi - psi_a) / (1.0 + m.dp)
        if krad and _POL and not v2:
            p.bsi += akxn[n] / alx[n] * m.xi * m.yi          # tbendi: bsi(i) += akxn(n)/alx(n)·xi·yi
        if krad and _POL and v2 and n == n2:
            p.bsi -= m.yi / m.rhob + ak / al * m.xi * m.yi   # v2 (B58): 出口の端
        if krad and n != n2:
            # tradk1 (cphi0 = cos(phixn(n)))。dp, dv(i) が更新される
            c[0], c[1], c[2], c[3], c[4], c[5] = m.xi, m.pxi, m.yi, m.pyi, m.zi, m.dp
            _tradk_core(p, alr[n], phixn[n])
            p.bpadd = 0.0
            m.xi, m.pxi, m.yi, m.pyi, m.zi, m.dp = c[0], c[1], c[2], c[3], c[4], c[5]
            p.setr0()
            if _PH is not None and _PH.get('ctx') is not None:
                _PH['ctx']['fr0'] += alx[n] / al       # tbendi.f:141 pcvt%fr0 += alx(n)/al
            if v2:
                m.rhoe = m.rhob * (1.0 + m.dp)     # v2: 放射で変わった運動量で曲率半径を更新する (K64 は入口の値のまま: 減衰の和則が崩れる)
            if _POL:
                p.bsi = 0.0
    _tbendicorr(m, akx0 * 0.5, alx0 * 0.5, phix0 * 0.5)
    m.zi = m.zi - p.dv * al
    if v2edge and sinp2 != 0.0:
        pe = 1.0 + m.dp
        if _compat.bend_k1_edge() == 'face':
            kt = ak / al * tanp2 / pe
            m.pxi += kt * m.xi * m.xi
            m.pyi -= kt * m.xi * m.yi
        ce = _lat._tbedge_cod(np.array([m.xi, m.pxi * pe, m.yi, m.pyi * pe, m.zi, m.dp]), m.rhob, math.atan2(sinp2, cosp2), False)
        c[0], c[1], c[2], c[3], c[4], c[5] = ce[0], ce[1] / pe, ce[2], ce[3] / pe, ce[4], ce[5]
    else:
        px = m.pxi + tanp2 * m.xi / m.rhoe
        f = m.yi / m.rhoe
        fpx = af * px
        c[1] = px
        c[2] = m.yi
        c[3] = m.pyi - (tanp2 - fpx) * f
        ff = af * m.yi * f * 0.5
        c[0] = m.xi - ff
        c[4] = m.zi + ff * fpx
        c[5] = m.dp
    if fb2 != 0.0:
        dxfr2 = fb2**2 / m.rhob / 24.0
        dyfr2 = fb2 / m.rhob**2 / 6.0
        dyfra2 = 4.0 * dyfr2 / fb2**2
        dp = c[5]
        pr = 1.0 + dp
        c[0] = c[0] - dxfr2 * dp / pr
        c[3] = c[3] + (dyfr2 - dyfra2 * c[2]**2) * c[2] / pr**2
        c[4] = c[4] - (dxfr2 * c[1] - (0.5 * dyfr2 - 0.25 * dyfra2 * c[2]**2) * c[2]**2 / pr) / pr
    if fringe and mfring > -4 and mfring != 1:
        _ttfrin4(p, -ak, al)
    if krad and v2:
        p.zr0 -= _zr0_edge(p, tanp2, 1.0)                  # v2: 出口の端の角の分 (B220)
    # 最後の tradk (全粒子分をまとめて: tbendi.f は粒子ループの後に call tradk(np,…,alr(n2),phixn(n2)))
    # と tbrot/tbshift の戻しは呼び出し側 (_tbend_ens) で行う
    return (alr[n2], phixn[n2], dtheta, dchi2, alg, phig, dx, dy, cost, sint)


_RBH = 0.5 + 1.0 / math.sqrt(12.0)
_RBL = 1.0 / 6.0 / _RBH


def _tbendal(n, ndiv, f1r, f2r, aln):
    """tbend.f tbendal: (alx, alr)"""
    if 0 < n <= ndiv:
        return aln, aln
    if n == -1:
        return _RBL * f1r, f1r
    if n == 0:
        return _RBH * f1r, f1r
    if n == ndiv + 1:
        return _RBH * f2r, f2r
    return _RBL * f2r, f2r


def _tsteer_ens(parts, el, krad):
    """tsteer.f tsteer (tbend.f から phib 引数 = −K0 で呼ばれる)。F1 フリンジ (fb1, fb2) は未対応。
    スライス外側・粒子内側で、各スライスの後に全粒子 tradk (SAD と同じ乱数の順)。
    入口エッジは tanp1 = tan(AE1) (最初のスライスだけ)、出口の tanp2 は fb2≠0 のときだけ入る (fb2=0 なら 0 のまま: SAD のまま)。"""
    al = el.L
    ak0 = _g(el, 'K0')
    dx, dy = _g(el, 'DX'), _g(el, 'DY')
    theta = _lat._bend_theta(el) + _g(el, 'DROTATE')      # tbend.f: tsteer には theta+dtheta (CALC の tsteee は theta だけ)
    if al == 0.0:
        if ak0 == 0.0:
            return
        for p in parts:
            c = p.c
            cost, sint = _tent(c, dx, dy, theta)
            c[1] -= ak0 / (1.0 + c[5])
            _texit(c, dx, dy, theta, cost, sint)
        return
    if ak0 == 0.0:
        for p in parts:
            _tdrift_free(p, al)
        return
    phib = -ak0
    rhob = al / phib
    psi1 = _g(el, 'E1')
    psi2 = _g(el, 'E2')
    tanp1 = math.sin(psi1) / math.cos(psi1)
    tanp2 = 0.0
    # F1・FB1・FB2 (tffs.f tpara: fb1 = F1+FB1、fb2 = F1+FB2。FRINGE ≤ 0 は −1 (入口だけ)・−2 (出口だけ) 以外 0。tsteer.f の fringe = DISFRIN==0)。
    # 入口 dxfr1 = −fb1²/ρb/24、dyfr1 = fb1/ρb²/6、断片 n = −1, 0 (長さ f1r = fb1/2 の放射)、出口 n = ndiv+1, ndiv+2、dxfr2 = −fb2²/ρb/24 (tsteer の符号)
    fb1, fb2, mfr, fringe_fl = _lat._bend_fringe_params(el)
    if mfr <= 0:
        if mfr != -1:
            fb1 = 0.0
        if mfr != -2:
            fb2 = 0.0
    f1r = f2r = 0.0
    dxfr1 = dyfr1 = dyfra1 = dxfr2 = dyfr2 = dyfra2 = 0.0
    n1, n2 = 1, 0
    if fb1 != 0.0:
        dxfr1 = -fb1**2 / rhob / 24.0
        dyfr1 = fb1 / rhob**2 / 6.0
        if fringe_fl:
            dyfra1 = 4.0 * dyfr1 / fb1**2
        if krad:
            f1r = 0.5 * fb1
            n1 = -1
    if fb2 != 0.0 and krad:
        f2r = 0.5 * fb2
        n2 = 2
    phic = phib * (1.0 - (f1r + f2r) / al)
    ctsn = []
    for p in parts:
        ctsn.append(_tent(p.c, dx, dy, theta))
        _spin_rot_cs(p, math.cos(theta), math.sin(theta))            # v2 (P1): ROTATE のスピンの座標系
    ndiv = 1
    if krad:
        for p in parts:
            p.setr0()
        ndiv = min(1000, max(1, _ndivrad(phic, 0.0, 0.0, _g(el, 'EPS'), parts[0].ph.h0)))
    alc = al * phic / phib
    aln = alc / ndiv
    v2 = _compat.v2()
    # v2 (B221): 端の面 (軸に対して AE 傾いた平面) を、CALC (tsteee) と同じ厳密な硬い端の鎖 (tbedgedrift → Maxwell → tbedgebody) で通す。
    # K64 の tsteer は x tan(e)/ρ の薄い蹴りだけで、K0³ tan²e の次数で厳密な写像と違い、ヤコビアンが K0² tan²e 違う (CALC との差の原因)。
    # F1・FB1・FB2 (tblfre の断片) があるときは K64 の薄い蹴りのまま (厳密な端の鎖と F1 の組み合わせは未確認)
    # v2 (B225): FRINGE = −2 (出口だけ) は入口の端、FRINGE = −1 (入口だけ) は出口の端を飛ばす (CALC の tsteee と同じ。K64 の追跡は両方の端を入れる)。
    # 磁石が隣と続く接合部の端 (場が続く) で縦収束・面の効果を入れないのが意図で、CALC の側が物理的に正しい
    sk1 = v2 and mfr == -2
    sk2 = v2 and mfr == -1
    ex1 = v2 and psi1 != 0.0 and not sk1
    ex2 = v2 and psi2 != 0.0 and not sk2
    if sk1:
        tanp1 = 0.0
    rhoc = -rhob          # 端の鎖に渡す曲率半径 (CALC の rhoe = L/K0)
    if ex1:
        if krad:
            for p in parts:
                p.zr0 -= _zr0_edge(p, tanp1, -1.0)      # 入口の面までの直進の行程 (B220)
        tanp1 = 0.0          # 薄い蹴りは使わない。鎖は最初の断片の tblfre の後に通す (CALC の tsteee と同じ順序)
    elif v2 and krad and tanp1 != 0.0:
        for p in parts:
            p.zr0 -= _zr0_edge(p, tanp1, -1.0)          # v2 (P3、B108、B220): 入口の端の角で場の中の行路が短い (tradk の al1 = al − z + zr0)
    for n in range(n1, ndiv + n2 + 1):
        alx, alr = _tbendal(n, ndiv, f1r, f2r, aln)
        if n == n1 + 1:
            dxfr1 = dyfr1 = dyfra1 = 0.0
            tanp1 = 0.0
        last_ex = ex2 and n == ndiv + n2
        last_skip = sk2 and n == ndiv + n2
        if n == ndiv + n2:
            if fb2 != 0.0:
                dxfr2 = -fb2**2 / rhob / 24.0
                dyfr2 = fb2 / rhob**2 / 6.0
                if fringe_fl:
                    dyfra2 = 4.0 * dyfr2 / fb2**2
            if (fb2 != 0.0 or v2) and not last_ex and not last_skip:
                # K64 の tsteer は出口の端の角 (AE2) を fb2≠0 のときだけ入れる (Emittance[]・CALC の tsteee は常に入れる: 不整合)。v2 は常に入れる (P3)
                tanp2 = math.sin(psi2) / math.cos(psi2)
        first_ex = ex1 and n == n1
        first_skip = sk1 and n == n1
        for p in parts:
            c = p.c
            dp = c[5]
            pr = 1.0 + dp
            rhoe = rhob * pr
            pxi = c[1] - tanp1 * c[0] / rhoe
            c[0] += dxfr1 * dp / pr
            c[3] += ((dyfr1 - dyfra1 * c[2]**2) / pr**2 + tanp1 / rhoe) * c[2]
            c[4] += (dxfr1 * pxi + (0.5 * dyfr1 - 0.25 * dyfra1 * c[2]**2) * c[2]**2 / pr) / pr
            if first_ex or first_skip:
                # 厳密な端の鎖 (Maxwell の項を含む) は tblfre の後 (CALC の tsteee と同じ順序)。飛ばす端 (FRINGE = −2) は何もしない
                if first_ex:
                    _steer_edge_exact(p, rhoc, psi1, True)
                pxi = c[1]
                pyi = c[3]
                xi = c[0]
            else:
                dpv1 = _sqrt1(-pxi**2)
                pv1 = 1.0 + dpv1
                fa = c[2] / rhoe / pv1
                f = (1.0 - (c[2] / rhob)**2 / 6.0) * fa
                ff = (0.5 - (c[2] / rhob)**2 / 24.0) * c[2] * fa / pv1**2
                pyi = c[3] + pxi * f
                xi = c[0] - ff
                c[4] += ff * pxi
            dpz1 = _pxy2dpz(pxi, pyi)
            pz1 = 1.0 + dpz1
            dpx = alx / rhoe
            pxf = pxi + dpx
            dpz2 = _pxy2dpz(pxf, pyi)
            pz2 = 1.0 + dpz2
            d = pxf * pz1 + pxi * pz2
            if d == 0.0 or pxi * pxf < 0.0:
                phsq = pxi**2 + pz1**2
                u = (dpx + pxf * dpz1 - pxi * dpz2) / phsq
                w = (pxf * dpz1 - pxi * dpz2 + dpx * pyi**2) / phsq
            else:
                u = dpx * (pxf + pxi) / d
                w = -dpx * (pxf * dpz1 + pxi * dpz2) / d
            dl = rhoe * (_lat._asinx(u) + w)
            c[0] = xi + alx * (pxf + pxi) / (pz2 + pz1)
            c[2] += pyi * (alx + dl)
            if last_ex or last_skip:
                # 出口の Maxwell の項と面は厳密な端の鎖 (下) の中にある。飛ばす端 (FRINGE = −1) は何もしない
                c[3] = pyi
                c[1] = pxf
                c[4] -= (dl + p.dv * alx)
                if last_ex:
                    _steer_edge_exact(p, rhoc, psi2, False)
            else:
                dpv2 = _sqrt1(-pxf**2)
                pv2 = 1.0 + dpv2
                fa = c[2] / rhoe / pv2
                f = (1.0 - (c[2] / rhob)**2 / 6.0) * fa
                ff = (0.5 - (c[2] / rhob)**2 / 24.0) * c[2] * fa / pv2**2
                c[3] = pyi - pxf * f
                c[0] += ff
                c[1] = pxf - c[0] * tanp2 / rhoe
                c[4] -= (dl + p.dv * alx + ff * pxf)
            c[0] -= dxfr2 * dp / pr
            c[3] += ((dyfr2 - dyfra2 * c[2]**2) / pr**2 + tanp2 / rhoe) * c[2]
            c[4] -= (dxfr2 * pxf - (0.5 * dyfr2 - 0.25 * dyfra2 * c[2]**2) * c[2]**2 / pr) / pr
        if krad:
            for p in parts:
                if last_ex:
                    p.zr0 -= _zr0_edge(p, math.tan(psi2), 1.0)  # 出口の面から z = L の面までの直進の行程 (B220)
                elif v2 and n == ndiv + n2 and tanp2 != 0.0:
                    p.zr0 -= _zr0_edge(p, tanp2, 1.0)  # v2 (P3、B108、B220): 出口の端の角の分
                _tradk(p, alr, 0.0)
    for p, (cost, sint) in zip(parts, ctsn):
        _texit(p.c, dx, dy, theta, cost, sint)
        _spin_rot_cs(p, math.cos(theta), -math.sin(theta))           # v2 (P1): 出口は逆回転


def _tbend_ens(parts, el, krad, first=True, last=True):
    """tbend.f tbend。全粒子をまとめて進める。K1≠0 は tbendi (粒子外側)。
    first/last = iprev/inext が 0 か (F1 の分割)。"""
    al0 = el.L
    phi0 = el.ANGLE
    ak0 = _g(el, 'K0')
    if phi0 == 0.0:
        # tbend.f: phi0=0, K1=0 → tsteer(.., −ak0, ..) (tsteer.f)。K0>0 は px を減らす (ANGLE≠0 の phib=ANGLE+K0 が曲げを強める向き)。
        # L=0 → tthin(nord=2, ak=K0): px −= K0/(1+δ) (TENT/TEXIT のずれ・ロール付き)、K0=0 → ドリフト
        if _g(el, 'K1') != 0.0:
            _tbend_k1_ens(parts, el, krad)
            return
        _tsteer_ens(parts, el, krad)
        return
    if al0 == 0.0:
        raise NotImplementedError('track6: L=0 で ANGLE≠0 の BEND (tbthin) は未対応')
    # tpara (tffs.f)
    psi1 = _g(el, 'E1')            # pysad の E1 は rad (= E1*ANGLE)
    psi2 = _g(el, 'E2')
    fb1, fb2, mfring, fringe = _lat._bend_fringe_params(el)
    if mfring <= 0:
        if mfring != -1:
            fb1 = 0.0
        if mfring != -2:
            fb2 = 0.0
    w = phi0 - psi1 - psi2
    al = al0
    if fb1 != 0.0 or fb2 != 0.0:
        al = al - ((phi0 * fb1)**2 + (phi0 * fb2)**2) / al / 48.0 \
            * math.sin(0.5 * w) / math.sin(0.5 * phi0)
    cosp1, sinp1 = math.cos(psi1), math.sin(psi1)
    cosp2, sinp2 = math.cos(psi2), math.sin(psi2)
    cosw, sinw = math.cos(w), math.sin(w)
    sqwh = sinw**2 / (1.0 + cosw) if cosw >= 0.0 else 1.0 - cosw
    sinwp1 = math.sin(phi0 - psi2)
    # tbend: ずれ (tbshift)、ロール誤差 DROTATE (tbrot, BEND の中点まわり: LGEO=L/2, ANGLGEO=ANGLE/2)
    dx, dy = _g(el, 'DX'), _g(el, 'DY')
    theta = _lat._bend_theta(el)        # ROT = 設計ロール ROTATE + ROTATION
    cost, sint = math.cos(theta), math.sin(theta)
    dtheta = _g(el, 'DROTATE')
    dchi2 = _g(el, 'CHI2')
    alg, phig = 0.5 * al0, 0.5 * phi0
    phib = phi0 + ak0
    ak1 = _g(el, 'K1') * al0            # pysad の Bend.K1 は長さあたり → SAD の積分値
    if ak1 != 0.0:
        tail = None
        for p in parts:                  # tbendi は粒子外側・スライス内側 (乱数の順もそのまま)
            tail = _tbendi(p, al, phib, phi0, cosp1, sinp1, cosp2, sinp2, ak1, dx, dy, theta,
                           dtheta, cost, sint, dchi2, alg, phig, fb1, fb2, mfring, krad, fringe,
                           _g(el, 'EPS'), first, last)
        alr_e, phi_e = tail[0], tail[1]
        if krad:
            if _PH is not None:
                _ph_set(dx, dy, theta, dtheta, phi0, phig, al, 1.0 - alr_e / al)      # tbendi.f: tsetpcvt, fr0 = 1 − alr(n2)/al
            for p in parts:              # tbendi.f: 粒子ループの後に全粒子の tradk
                _tradk(p, alr_e, phi_e)
        for p in parts:
            if dtheta != 0.0 or dchi2 != 0.0:
                _tbrot(p, alg - al, phig - phi0, dtheta, dchi2, False)
            _tbshift(p, -dx, -dy, -phi0, cost, -sint, False)
        return
    for p in parts:
        _tbshift(p, dx, dy, phi0, cost, sint, True)
        if dtheta != 0.0 or dchi2 != 0.0:
            _tbrot(p, alg, phig, dtheta, dchi2, True)
    if krad and _PH is not None:
        _ph_set(dx, dy, theta, dtheta, phi0, phig, al, 0.0)                         # tbend.f: tsetpcvt (光子の位置の情報)
    _tbend0_ens(parts, el, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2,
                cosw, sinw, sqwh, sinwp1, fb1, fb2, mfring, fringe, krad, first, last)
    for p in parts:
        if dtheta != 0.0 or dchi2 != 0.0:
            _tbrot(p, alg - al, phig - phi0, dtheta, dchi2, False)
        _tbshift(p, -dx, -dy, -phi0, cost, -sint, False)


def _bend_k1_proxy(el):
    """tbend.f:566-590 (phi0=0, K1≠0): K0=0 なら tquad(al, ak=K1·L, bz=0, dx, dy, theta+dtheta, theta+dtheta, krad, chro, fringe,
    f1in=f2in=f1out=f2out=0, mfring=0, eps) 相当の QUAD の代理、K0≠0 なら tmulti(akm(0)=K0, akm(1)=K1·L, theta+dtheta) 相当の
    MULT の代理 (lattice._bend_as_mult に ROTATION=theta+dtheta を入れる)。theta = ROTATE + ROTATION (_bend_theta)、dtheta = DROTATE"""
    from types import SimpleNamespace
    th = _lat._bend_theta(el) + _g(el, 'DROTATE')
    if _g(el, 'K0') == 0.0:
        # SAD K64 の癖 (29_sad_quirks A35): tbend.f:573 は tquad に theta2=theta+dtheta を渡し、akang (K1<0 の π/2) を含めない。
        # tquad は ak=sign(ak0,al) で符号を theta2 に任せるので、K1<0 の BEND は |K1| の集束四極として追跡される
        # (SAD で QUAD K1=+0.5 と BEND ANGLE=0 K1=-0.5 の追跡が全桁同じ、CALC は正しく -0.5)。それを再現する
        return 'quad', SimpleNamespace(L=el.L, K1=abs(_g(el, 'K1')), DX=_g(el, 'DX'), DY=_g(el, 'DY'), ROTATION=th,
                                       DISFRIN=_g(el, 'DISFRIN'), FRINGE=0.0, F1=0.0, F2=0.0, F1K1F=0.0, F1K1B=0.0,
                                       F2K1F=0.0, F2K1B=0.0, EPS=_g(el, 'EPS'), RAD=_g(el, 'RAD'))
    m = _lat._bend_as_mult(el)
    m.ROTATION = th
    m.DROTATE = 0.0
    m.RAD = _g(el, 'RAD')
    # tbend.f:584-591: tmulti に mfring・fb10 (= F1+FB1)・fb20 (= F1+FB2) を渡す (f1in 等 = 0、端の角 psi = 0。v2next P3)
    m.FB1 = _g(el, 'F1') + _g(el, 'FB1')
    m.FB2 = _g(el, 'F1') + _g(el, 'FB2')
    # 端の角: K64 の追跡 (tbend → tmulti) は psi1 = psi2 = 0 を渡す (AE1・AE2 を無視) が、Emittance[]・CALC (tbende → tmulte) は AE を使うので、
    # 同じ BEND の追跡と CALC が合わない。v2 は追跡も AE を使う (P3。gs はこの BEND (ANGLE=0、K0≠0、K1≠0) の USE で SIGSEGV するので K64 の確認は不可)
    m.AE1 = _g(el, 'AE1') if _compat.v2() else 0.0
    m.AE2 = _g(el, 'AE2') if _compat.v2() else 0.0
    return 'mult', m


def _tbend_k1_ens(parts, el, krad):
    """ANGLE=0, K1≠0 の BEND の追跡 (tbend.f: tquad / tmulti)"""
    kind, proxy = _bend_k1_proxy(el)
    if kind == 'quad':
        _tquad_ens(parts, proxy, krad)
    else:
        _tmulti_ens(parts, proxy, krad)


def _tbend(p: _Particle, el, krad, first=True, last=True):
    _tbend_ens([p], el, krad, first, last)


def _tbend0_ens(parts, el, al, phi0, ak0, phib, psi1, psi2, cosp1, sinp1, cosp2, sinp2,
                cosw, sinw, sqwh, sinwp1, fb1, fb2, mfring, fringe, krad, first, last):
    """tbend の tbshift/tbrot の間 (K1=0): tbendcore 1 回、または放射ありで tbendr
    (tbendr はスライス外側: 各スライスで全粒子の本体 → 全粒子の tradk)"""
    rhob = al / phib
    rho0 = al / phi0
    drhob = -rho0 * ak0 / phib
    bp = (rhob, rho0, drhob, fb1, fb2)
    n1, n2, ndiv = 1, 0, 1
    f1r = f2r = 0.0
    v2 = _compat.v2()
    if krad:
        for p in parts:
            p.setr0()
            if v2 and sinp1 != 0.0:
                p.zr0 -= _zr0_edge(p, sinp1 / cosp1, -1.0)    # v2: 入口の端の角で場の中の行路が短い (B220)
        if first and fb1 != 0.0:
            n1 = -1
            f1r = fb1 * 0.5
        if last and fb2 != 0.0:
            n2 = 2
            f2r = fb2 * 0.5
        phir = phib * (al - f1r - f2r) / al
        ndiv = min(1024, _ndivrad(phir, 0.0, 0.0, _g(el, 'EPS'), parts[0].ph.h0))
    n2 = ndiv + n2
    if n1 == n2:
        for p in parts:
            _tbendcore(p, al, phi0, cosp1, sinp1, cosp2, sinp2, mfring, fringe,
                       cosw, sinw, sqwh, sinwp1, krad, al, bp, 1.0, 1.0)
        return
    # tbendr
    aln = (al - f1r - f2r) / ndiv
    coswn = sinwn = sqwhn = sinwpn = 0.0
    for n in range(n1, n2 + 1):
        mfr1 = 0
        psi1n = psi2n = 0.0
        cosp1n, sinp1n, cosp2n, sinp2n = 1.0, 0.0, 1.0, 0.0
        alx, alr = _tbendal(n, ndiv, f1r, f2r, aln)
        if n == n1:
            if mfring > 0 or mfring == -1:
                mfr1 = -1
            psi1n, cosp1n, sinp1n = psi1, cosp1, sinp1
        elif n == n2:
            if mfring > 0 or mfring == -2:
                mfr1 = -2
            psi2n, cosp2n, sinp2n = psi2, cosp2, sinp2
        phi0n = alx / al * phi0
        if n <= 2 or n >= ndiv:
            wn = phi0n - psi1n - psi2n
            sinwn, xsinwn, coswn, sqwhn = _xsincos(wn)
            sinwpn = math.sin(phi0n - psi2n)
        for p in parts:
            _tbendcore(p, alx, phi0n, cosp1n, sinp1n, cosp2n, sinp2n, mfr1, fringe,
                       coswn, sinwn, -sqwhn, sinwpn, False, alr, bp,
                       1.0 if n == n1 else 0.0, 1.0 if (n == n2 and n != n1) else 0.0)
        for p in parts:
            if n == n2 and v2 and sinp2 != 0.0:
                p.zr0 -= _zr0_edge(p, sinp2 / cosp2, 1.0)    # v2: 出口の端の角の分 (B220)
            _tradk(p, alr, phi0n)
        if _PH is not None and _PH.get('ctx') is not None:
            _PH['ctx']['fr0'] += alx / _PH['ctx']['al']        # tbend.f:367 pcvt%fr0 += alx/al


# ---------------------------------------------------------------------------
# 空洞 (tcav.f tcav, L=0)
# ---------------------------------------------------------------------------

def _tcav(p: _Particle, el, trf0, rfsw, rad):
    """tcav.f tcav (V1=V20=V11=V02=0, DPHI=0, autophi=F, trpt=F → vnominal=0)。
    vnominal=0 では分割数 ndiv=1 (r=0 で sqrt(r²/3/eps)=0)。L≠0 では
    入口フリンジ tcavfrin → L/2 ドリフト → キック → L/2 ドリフト → 出口フリンジ (FRINGE≥0, DISFRIN=0)。
    rfsw=F でも tturn1 は空洞を呼ぶが、tlim=0 で t=0 になる。"""
    al = el.L
    ph = p.ph
    vc = _g(el, 'VOLT') + _g(el, 'DVOLT')
    harm = _g(el, 'HARM') or _g(el, '_harm')
    freq = _g(el, 'FREQ')
    w = 2.0 * math.pi * freq / CVELOC if harm == 0.0 else ph.omega0 * harm / CVELOC
    wi = 1.0 / w if w != 0.0 else 0.0
    charge = ph.charge
    phic = _g(el, 'PHI') * charge
    v = vc / ph.amass * abs(charge)
    phis = w * trf0
    p0 = ph.p0
    h0 = ph.h0
    pe = p0
    vcorr = v * (w * (0.5 / p0 + 0.5 / pe))**2 / 4.0
    v20a = vcorr
    v02a = vcorr
    qa = abs(charge) / ph.amass
    v1a, v11a = _g(el, 'V1') * qa, _g(el, 'V11') * qa            # tcav.f: v1a, v20a, v11a, v02a (V1 V20 V11 V02)
    vgen = any(q != 0.0 for q in (_g(el, 'V1'), _g(el, 'V20'), _g(el, 'V11'), _g(el, 'V02')))
    if vgen:
        v20a += 0.5 * _g(el, 'V20') * qa
        v02a += 0.5 * _g(el, 'V02') * qa
    dphis = phis - phic
    if rad:
        offset = math.sin(dphis)
        offset1 = 0.0
    else:
        offset = math.sin(dphis) - math.sin(phis)
        offset1 = math.sin(phis)
    tlim = 1e4 if rfsw else 0.0
    wsn = 1.0
    c = p.c
    fringe = _g(el, 'DISFRIN') == 0.0
    mfring = int(round(_g(el, 'FRINGE')))
    if al != 0.0 and fringe and mfring >= 0 and mfring != 2:
        _tcavfrin(p, al, v, w, dphis, offset)
    if al != 0.0:
        _tcav_drift(p, al * wsn * 0.5)
    dp1r = c[5]
    p1r = 1.0 + dp1r
    p1 = p0 * p1r
    if p.dv > 0.1:
        h1 = _e._p2h(p1)
    else:
        h1 = p1r * h0 / (1.0 - p.dv + p.ph.dvfs)
    t = min(tlim, max(-tlim, -c[4] * h1 / p1))
    phh = 0.5 * w * t
    vx = (v + (v1a + v20a * c[0] + v11a * c[2]) * c[0] + v02a * c[2]**2) if vgen else (v + v20a * c[0]**2 + v02a * c[2]**2)
    dh = max(1.0 + 3.83e-12 - h1, vx * (-2.0 * math.sin(phh) * math.cos(phh - dphis) + offset) * wsn)
    h2 = h1 + dh
    a = max(dh * (h1 + h2) / p1**2, -1.0)
    dpr = _sqrt1(a)
    dp2r = dp1r + p1r * dpr
    p2r = 1.0 + dp2r
    c[5] = dp2r
    p.dv = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + p.ph.dvfs
    av = -(math.cos(2.0 * phh - dphis) * wi - offset1 * t) * wsn / p0
    if vgen:
        dpx = (v1a + 2.0 * v20a * c[0] + v11a * c[2]) * av
        dpy = (v11a * c[0] + 2.0 * v02a * c[2]) * av
    else:
        dpx = 2.0 * v20a * c[0] * av
        dpy = 2.0 * v02a * c[2] * av
    c[1] = (c[1] * p1r + dpx) / p2r
    c[3] = (c[3] * p1r + dpy) / p2r
    c[4] = -t * p2r * p0 / h2
    if al != 0.0:
        _tcav_drift(p, al * wsn * 0.5)
        if fringe and mfring >= 0 and mfring != 1:
            _tcavfrin(p, al, -v, w, dphis, offset)


def _tcav_drift(p: _Particle, alx):
    """tcav.f のインラインのドリフト (do 10 / do 30)"""
    c = p.c
    a = c[1]**2 + c[3]**2
    dpz = _sqrt1(-a)
    al1 = alx / (1.0 + dpz)
    c[0] += c[1] * al1
    c[2] += c[3] * al1
    c[4] += dpz * al1 - p.dv * alx


def _tcavfrin(p: _Particle, al, v, w, dphis, offset):
    """tcav.f tcavfrin (空洞の入口/出口フリンジ。出口は v → −v で呼ぶ)。dv は dvfs を含む"""
    ph = p.ph
    p0, h0 = ph.p0, ph.h0
    vf = v / al / p0 * 0.5
    s0 = math.sin(dphis) - offset
    c = p.c
    dp1r = c[5]
    p1r = 1.0 + dp1r
    p1 = p0 * p1r
    if p.dv > 0.1:
        h1 = _e._p2h(p1)
    else:
        h1 = p1r * h0 / (1.0 - p.dv + ph.dvfs)
    t = -c[4] * h1 / p1
    phw = w * t
    sph = math.sin(phw - dphis) + s0
    if _compat.v2():
        # v2 (B330): 端の力は ∫F dt = ∫F ds/(β c) なので 1/β (粒子の β = p1/h1) がかかる。H_f ∝ 1/β(E) の整合のため、
        # エネルギー変化 (∂H_f/∂t) にも同じ 1/β、時間のずれ −∂H_f/∂E = −(x²+y²) κ0/(2 p1³) (κ0 = vf p0 sph) を入れる (シンプレクティック性)
        ib = h1 / p1
        dpt = vf * sph * ib / p1r
        c[1] += c[0] * dpt
        c[3] += c[2] * dpt
        r2 = c[0]**2 + c[2]**2
        dh = -w * p0 * vf * r2 * 0.5 * math.cos(phw - dphis) * ib
        dt = -0.5 * r2 * vf * p0 * sph / p1**3
        h2 = h1 + dh
        a = max(dh * (h1 + h2) / p1**2, -1.0)
        dpr = _sqrt1(a)
        dp2r = max(dp1r + p1r * dpr, 1e-10 - 1.0)
        p2r = 1.0 + dp2r
        c[5] = dp2r
        c[1] *= p1r / p2r            # v2 (B332): p が変わるので q = p_x/p を再スケール (正準運動量 p_x を保つ)
        c[3] *= p1r / p2r
        c[4] = -(t + dt) * p0 * p2r / h2
        p.dv = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + ph.dvfs
        return
    dpt = vf * sph / p1r
    c[1] += c[0] * dpt
    c[3] += c[2] * dpt
    dh = -w * p0 * vf * (c[0]**2 + c[2]**2) * 0.5 * math.cos(phw - dphis)
    h2 = h1 + dh
    a = max(dh * (h1 + h2) / p1**2, -1.0)
    dpr = _sqrt1(a)
    dp2r = max(dp1r + p1r * dpr, 1e-10 - 1.0)
    p2r = 1.0 + dp2r
    c[5] = dp2r
    c[4] = c[4] * p2r / p1r * h1 / h2
    p.dv = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + ph.dvfs


# ---------------------------------------------------------------------------
# MULT (tmulti.f tmulti; ANGLE=0, K0=K1=0, BZ=0, ずれなし)
# ---------------------------------------------------------------------------

_EPS00M = 0.005          # toplvl.f eps00m
_AMPMAXM = 0.05          # toplvl.f ampmaxm


def _ttfrin2(p: _Particle, ak, al):
    """tquad.f ttfrin (nord=2): 双極 (K0) のフリンジ"""
    c = p.c
    aki = ak / al / (1.0 + c[5])
    c[0] += 0.5 * aki * c[2]**2
    c[3] -= aki * c[1] * c[2]
    c[4] -= 0.5 * aki * c[1] * c[2]**2


def _ttfrin_gen(p: _Particle, nord, ak, al):
    """tquad.f ttfrin の一般次数 (nord ≥ 8, kord = nord/2) の分岐。an(k) = k"""
    kord = nord // 2
    c = p.c
    akk = ak / al / math.factorial(kord) / 4.0
    aki = akk / (1.0 + c[5])
    cx = complex(c[0], c[2])
    cp = complex(c[1], -c[3])
    cz1 = cx**(kord - 1)
    cz = cz1 * cx
    a = aki * cz.imag * 2.0
    ca = -aki * cx * (cz / (kord + 1) - cz.conjugate())
    cx1 = cx + ca
    d = 1.0 + a**2 - (aki * kord)**2 * (cz.real**2 + cz.imag**2)
    cp1 = (complex(1.0, a) * cp - (aki * kord) * cz1 * (cx * cp).conjugate()) / d
    c[4] = c[4] - (ca * cp1).real
    c[0] = cx1.real
    c[2] = cx1.imag
    c[1] = cp1.real
    c[3] = -cp1.imag


def _ttfrin_n(p: _Particle, nord, ak, al):
    """tquad.f ttfrin (bz=0, 実数 ak)。nord = 2, 4, 6, 8…"""
    if al == 0.0 or ak == 0.0:
        return
    if nord == 2:
        _ttfrin2(p, ak, al)
    elif nord == 4:
        _ttfrin4(p, ak, al)
    elif nord == 6:
        _ttfrin6(p, ak, al)
    else:
        _ttfrin_gen(p, nord, ak, al)


def _ttfrins(p: _Particle, nord, akc: complex, al):
    """tquad.f ttfrins: 複素 ak (回した座標系の係数) のフリンジ。主軸へ回して ttfrin、戻す"""
    if (akc.real == 0.0 and akc.imag == 0.0) or al == 0.0:
        return
    if akc.real == 0.0:
        theta = math.pi / nord
        aka = akc.imag
    elif akc.imag == 0.0:
        theta = 0.0
        aka = akc.real
    else:
        theta = math.atan2(akc.imag, akc.real) * 2.0 / nord
        aka = math.hypot(akc.real, akc.imag)
    c = p.c
    if theta != 0.0:
        ct, st = math.cos(theta), math.sin(theta)
        x0, px0 = c[0], c[1]
        c[0] = ct * x0 - st * c[2]
        c[2] = st * x0 + ct * c[2]
        c[1] = ct * px0 - st * c[3]
        c[3] = st * px0 + ct * c[3]
        _ttfrin_n(p, nord, aka, al)
        x0, px0 = c[0], c[1]
        c[0] = ct * x0 + st * c[2]
        c[2] = -st * x0 + ct * c[2]
        c[1] = ct * px0 + st * c[3]
        c[3] = -st * px0 + ct * c[3]
    else:
        _ttfrin_n(p, nord, aka, al)


def _tblfri(p: _Particle, al, ck0: complex, fb1):
    """tmulti.f tblfri: K0 (ck0.real = x 方向の曲げ、ck0.imag = y 方向) の端の非線形フリンジ (FB1/FB2)。呼ぶ側が入口 (+ck0, fb1)・出口 (−ck0, fb2) を決める"""
    dxfrx = dyfrx = dyfrax = dxfry = dyfry = dxfray = 0.0
    if ck0.real != 0.0:
        rhob = al / ck0.real
        dxfrx = fb1**2 / rhob / 24.0
        dyfrx = fb1 / rhob**2 / 6.0
        dyfrax = 4.0 * dyfrx / fb1**2
    if ck0.imag != 0.0:
        rhob = al / ck0.imag
        dyfry = fb1**2 / rhob / 24.0
        dxfry = fb1 / rhob**2 / 6.0
        dxfray = 4.0 * dxfry / fb1**2
    c = p.c
    pr = 1.0 + c[5]
    c[0] += dxfrx * c[5] / pr
    c[2] -= dyfry * c[5] / pr
    c[1] += (dxfry - dxfray * c[0]**2) / pr**2 * c[0]
    c[3] += (dyfrx - dyfrax * c[2]**2) / pr**2 * c[2]
    c[4] += (dxfrx * c[1] - dyfry * c[3] + ((0.5 * dyfrx - 0.25 * dyfrax * c[2]**2) * c[2]**2
                                           + (0.5 * dxfry - 0.25 * dxfray * c[0]**2) * c[0]**2) / pr) / pr


def _tmulk(p: _Particle, al, akr, apsi, nmmax):
    """tmulti.f tmulk: 端の角 apsi (AE1・AE2) による多極場の端の蹴り (fe = −tan(apsi)/al、(Σ K_n z^n/n!) z に比例)"""
    fe = -math.tan(apsi) / al
    c = p.c
    cx0 = complex(c[0], c[2])
    cx = 0j
    for kord in range(nmmax, 0, -1):
        cx = (cx + akr[kord]) * cx0 / kord
    cx = fe * (cx + akr[0]) * cx0
    pr = 1.0 + c[5]
    c[1] -= cx.real / pr
    c[3] += cx.imag / pr


def _mult_flow(p, a, al, akr, nmmax):
    """v2 (B223): MULT の本体の場 (長さあたりの強さは全長 al のものと同じ) を長さ a だけ通す Φ(a) (a < 0 は逆向き)。
    tsolqu の半分 + 薄い多極の蹴り (n ≥ 2) + tsolqu の半分 (tmulti の 1 刻みと同じ形)。a は端の面までの距離 (mm〜cm) なので 1 刻み"""
    if a == 0.0:
        return
    sc = a / al
    ak1 = (akr[1].real if nmmax >= 1 else 0.0) * sc * 0.5
    ak0 = akr[0] * sc * 0.5
    _tsolqu(p, a * 0.5, ak1, 0.0, ak0.real, ak0.imag, qum=False, ibsi=-1)
    if nmmax >= 2:
        c = p.c
        cx1 = complex(c[0], c[2])
        cx = akr[nmmax] * sc * cx1 / nmmax
        for n in range(nmmax - 1, 1, -1):
            cx = (cx + akr[n] * sc) * cx1 / n
        cx = cx * cx1 / (1.0 + c[5])
        c[1] -= cx.real
        c[3] += cx.imag
    _tsolqu(p, a * 0.5, ak1, 0.0, ak0.real, ak0.imag, qum=False, ibsi=-1)


def _mult_edge_ykick(p, apsi, al, akr, nmmax, sgn, th=None):
    """v2 (B224・B292): 端の面の縦方向の収束 (Maxwell の項)。面を横切る点での Δpy = −(tan(e + sgn·θ) − tan(sgn·θ))·Im ψ/(1+δ)、
    ψ = Σ_n (K_n/L) z^(n+1)/(n+1)! (横磁場のスカラーポテンシャル)、θ = atan(px/pz) は面を横切るときの水平方向の傾き (入口は面の外の入力の値で th に渡す、
    出口は面の内側 = 現在の値。sgn = +1 入口、−1 出口)。tan(sgn·θ) の分は本体 (tsolqu の双極の項) が持つ。lattice._mult_sequence の edge_ykick と同じ"""
    c = p.c
    if th is None:
        th = math.atan2(c[1], math.sqrt(max(1e-30, 1.0 - c[1] ** 2 - c[3] ** 2)))
    z = complex(c[0], c[2])
    g = 0j
    for n in range(nmmax, -1, -1):
        g = (g + akr[n] / al) * z / (n + 1)
    c[3] += -(math.tan(apsi + sgn * th) - math.tan(sgn * th)) * g.imag / (1.0 + c[5])


def _mult_edge_exact(p, entrance, apsi, al, akr, nmmax):
    """v2 (B223): MULT の端の面 (軸に対して AE 傾いた平面) の厳密な写像。入口は Φ(−z_f)∘Drift(z_f) (z_f = x tan e/(1 − (px/pz) tan e))、
    出口は Drift(ℓ)∘Φ(−ℓ) (ℓ = X_e tan e は Φ(−ℓ) で戻した軌道の x の不動点)。lattice._mult_sequence の exact_edge_in・exact_edge_out と同じ。
    K64 の tmulk は面を z = 0 の平面にして薄い蹴りだけを入れる"""
    t = math.tan(apsi)
    c = p.c
    if entrance:
        pzr = math.sqrt(max(1e-30, 1.0 - c[1] ** 2 - c[3] ** 2))
        zf = c[0] * t / (1.0 - c[1] / pzr * t)
        th_in = math.atan2(c[1], pzr)                       # 面を横切るときの傾き = 面の外 (入力) の値
        _tdrift_free(p, zf)
        _mult_edge_ykick(p, apsi, al, akr, nmmax, 1.0, th_in)      # 面の上 (Drift(z_f) の後、Φ(−z_f) の前) で縦の蹴り
        _mult_flow(p, -zf, al, akr, nmmax)
        return
    ell = t * c[0]
    for _ in range(40):
        q = copy.copy(p)
        q.c = p.c.copy()
        _mult_flow(q, -ell, al, akr, nmmax)
        new = t * q.c[0]
        done = abs(new - ell) < 1e-17 + 1e-15 * abs(ell)
        ell = new
        if done:
            break
    _mult_flow(p, -ell, al, akr, nmmax)
    _mult_edge_ykick(p, apsi, al, akr, nmmax, -1.0)         # 面の上 (Φ(−ℓ) の後、Drift(ℓ) の前)
    _tdrift_free(p, ell)


def _mult_fb_ae(el, pm):
    """tmulti の mfring・fb1・fb2・psi1・psi2: tpara (tffs.f) は FB1 を FRINGE = 1・3、FB2 を FRINGE = 2・3 のときだけ残す。
    端の角は psi = E1·ANGLE + AE1 (ANGLE=0 では AE1・AE2)"""
    mfr = pm['mfring']
    fb1 = pm['fb1'] if mfr in (1, 3) else 0.0
    fb2 = pm['fb2'] if mfr in (2, 3) else 0.0
    return mfr, fb1, fb2, _g(el, 'AE1'), _g(el, 'AE2')


def _tmulti_ens(parts, el, krad):
    """tmulti.f tmulti (ANGLE=0, BZ=0, SK なし, RF なし)。K0/K1/K2… とずれ・回転 (tsolrot) に対応。
    F1・F2 (四極の線形フリンジ、f1in 等)・FB1・FB2 (K0 の端の tblfri)・FRINGE = 1, 2, 3・AE1・AE2 (tmulk) に対応 (v2next P3)。
    tpara (tffs.f) の akr[n] = K_n·cr1^(n+1), theta2 = ROT+DROT+akang, ndiv0, dofr。
    各分割: tsolqu(al1, ak1, ak0) (+tradk 全粒子) → 多極キック (n≥2、(1+δ) で割る)。
    フリンジ: ttfrins で n=0 (双極), 1 (四極), 2 (6 極) 等 (dofr)。"""
    pm = _lat._mult_params(el)
    al, nmmax, akr, ndiv = pm['al'], pm['nmmax'], pm['akr'], pm['ndiv']
    if _g(el, 'ANGLE') or _g(el, 'VOLT'):
        raise NotImplementedError('track6: MULT の曲げ・RF は未対応')
    mfring, fb1, fb2, psi1, psi2 = _mult_fb_ae(el, pm)
    if nmmax < 0:
        for p in parts:
            _tdrift_free(p, al)
        return
    ct, st = math.cos(pm['theta2']), math.sin(pm['theta2'])
    dx, dy = pm['dx'], pm['dy']
    for p in parts:                      # tsolrot (ent, bz=0): ずらす → 回す
        p.c[0] -= dx
        p.c[2] -= dy
        _rot(p.c, ct, st)
        if _compat.v2():
            _spin_rot_z(p, pm['theta2'], True)       # v2 (B241): tsolrot の calpol はスピンも回す (QUAD の追跡・Emittance[] の tmulte は回す。MULT の追跡だけ欠けていた)
    akr1 = akr[1].real if nmmax >= 1 else 0.0
    akr0 = akr[0]
    eps = pm['eps']
    dofr = [akr[n] != 0 and (n <= 2 or abs(akr[n]) / al * 0.05**(n + 1) / (n + 1) > eps)
            for n in range(nmmax + 1)]
    fringe = pm['fringe']
    nz = al != 0.0
    # tpara: 追跡の四極の線形フリンジ (tsetfringep: f1 = −|K1/L·F1²|/24、f2 = |K1/L|·F2。tsetfringepe (Emittance[] 用の F1・F2 そのもの) とは別)
    f1in, f2in, f1out, f2out = _tsetfringep(el, akr1 / al) if nz else (0.0, 0.0, 0.0, 0.0)
    if nz:
        if krad:
            for p in parts:
                p.setr0()
            ndiv = max(ndiv, _ndivrad(abs(akr0), akr1, 0.0, _g(el, 'EPS'), parts[0].ph.h0))
            if _PH is not None:
                _ph_set(dx, dy, pm['theta2'], 0.0, 0.0, 0.0, al)      # tmulti.f: tsetpcvt (D1: MULT の光子)
        ndiv = min(1000, ndiv)
        if fringe and mfring != 2:
            for p in parts:
                for n in range(nmmax + 1):
                    if dofr[n]:
                        _ttfrins(p, n * 2 + 2, akr[n], al)
        if mfring == 1 or mfring == 3:
            for p in parts:
                if akr0 != 0 and fb1 != 0.0:
                    _tblfri(p, al, akr0, fb1)
                if f1in != 0.0 or f2in != 0.0:
                    _tqlfr(p, f1in, f2in)
    akr_k = list(akr)                    # tmulti: akr(1) = dble(akr(1)) (tmulk に渡す係数)
    if nmmax >= 1:
        akr_k[1] = akr1
    v2 = _compat.v2()
    if nz and psi1 != 0.0:
        for p in parts:
            if v2 and krad:
                p.zr0 -= _zr0_edge(p, math.tan(psi1), -1.0)       # v2 (P3、B108、B220): 入口の端の角で場の中の行路が短い (tradk の al1 = al − z + zr0)
            if v2:
                _mult_edge_exact(p, True, psi1, al, akr, nmmax)     # v2 (B223): 厳密な硬い端
            else:
                _tmulk(p, al, akr_k, psi1, nmmax)
    wi = 1.0 / ndiv
    ak1 = akr1 * wi * 0.5
    al1 = al * wi * 0.5
    ak01 = akr0 * wi * 0.5
    akrm = [a * wi for a in akr]
    ibsi = 1
    bpar_on = bool(krad and nz and _POL and v2 and _compat.spin_bpar() and nmmax >= 0)    # v2 (B165): 本体の B∥ = Δψ (ψ は横磁場のポテンシャル)

    def _psi(c):
        z = complex(c[0], c[2])
        g = 0j
        for n_ in range(nmmax, -1, -1):
            g = (g + akr[n_]) * z / (n_ + 1)
        return g.imag / al
    if bpar_on:
        psi_prev = [_psi(p.c) for p in parts]            # 入口の ψ (K64 の bsi の入口の点: 入口のフリンジと端の蹴りの後)
    for m in range(1, ndiv + 1):
        if nz:
            for p in parts:
                _tsolqu(p, al1, ak1, 0.0, ak01.real, ak01.imag, qum=(not krad) and not _compat.v2(), ibsi=ibsi)   # 放射なしは tsolqum (v2 (B14) は tsolqu の z の式)
            ibsi = 0
            if krad:
                if bpar_on:
                    # K64 の bsi (四極・二極の項は刻みの始め、n≥2 の項は最初の刻みの半分のあと) を、同じ点の ψ の全次数に置き換える
                    for ip_, p in enumerate(parts):
                        p.bsi = psi_prev[ip_] if m == 1 else 0.0
                        psi_c = _psi(p.c)
                        p.bpadd = (psi_c - psi_prev[ip_]) / (1.0 + p.c[5])
                        psi_prev[ip_] = psi_c
                elif m == 1 and _POL:
                    for p in parts:               # tmulti: bsi(i) += Im(.5·cx·cx1²)/al (cx は akr(n)·cx1·aninv(n+1) の和)
                        cx1 = complex(p.c[0], p.c[2])
                        cx = 0j
                        for n in range(nmmax, 1, -1):
                            cx = (cx + akr[n]) * cx1 / (n + 1)
                        p.bsi += (0.5 * cx * cx1**2).imag / al
                for p in parts:
                    _tradk(p, al1, 0.0)
                _ph_adv(al1 / al)
            wm = wi * 0.5 if m == ndiv else wi
            al1 = al * wm
            ak1 = akr1 * wm
            ak01 = akr0 * wm
            if nmmax >= 2:
                for p in parts:
                    c = p.c
                    cx1 = complex(c[0], c[2])
                    cx = akrm[nmmax] * cx1 / nmmax
                    for n in range(nmmax - 1, 1, -1):
                        cx = (cx + akrm[n]) * cx1 / n
                    cx = cx * cx1 / (1.0 + c[5])
                    c[1] -= cx.real
                    c[3] += cx.imag
        elif nmmax >= 1:
            for p in parts:
                c = p.c
                cx1 = complex(c[0], c[2])
                cx = akrm[nmmax] * cx1 / nmmax
                for n in range(nmmax - 1, 0, -1):
                    cx = (cx + akrm[n]) * cx1 / n
                cx = (cx + akrm[0]) / (1.0 + c[5])
                c[1] -= cx.real
                c[3] += cx.imag
    if nz:
        for p in parts:
            _tsolqu(p, al1, ak1, 0.0, ak01.real, ak01.imag, qum=(not krad) and not _compat.v2(), ibsi=2)
        if psi2 != 0.0:
            for p in parts:
                if v2:
                    _mult_edge_exact(p, False, psi2, al, akr, nmmax)    # v2 (B223): 厳密な硬い端 (出口の面は Drift(ℓ)∘Φ(−ℓ) の後の状態)
                    if krad:
                        p.zr0 -= _zr0_edge(p, math.tan(psi2), 1.0)   # v2 (P3、B108、B220): 出口の端の角の分
                else:
                    _tmulk(p, al, akr_k, psi2, nmmax)
        if mfring == 2 or mfring == 3:
            for p in parts:
                if f1out != 0.0 or f2out != 0.0:
                    _tqlfr(p, -f1out, f2out)
                if fb2 != 0.0 and akr0 != 0:
                    _tblfri(p, al, -akr0, fb2)
        if fringe and mfring != 1:
            for p in parts:
                for n in range(nmmax, -1, -1):
                    if dofr[n]:
                        _ttfrins(p, n * 2 + 2, -akr[n], al)
        if krad:
            if bpar_on:
                for ip_, p in enumerate(parts):
                    psi_c = _psi(p.c)
                    p.bsi = -psi_c
                    p.bpadd = (psi_c - psi_prev[ip_]) / (1.0 + p.c[5])
            elif _POL and nmmax >= 2:
                for p in parts:                   # tmulti: bsi(i) −= Im(.5·cx·cx1²)/al
                    cx1 = complex(p.c[0], p.c[2])
                    cx = akr[nmmax] * cx1 / (nmmax + 1)
                    for n in range(nmmax - 1, 1, -1):
                        cx = (cx + akr[n]) * cx1 / (n + 1)
                    p.bsi -= (0.5 * cx * cx1**2).imag / al
            for p in parts:
                _tradk(p, al1, 0.0)
            _ph_adv(al1 / al)
    for p in parts:                      # tsolrot (exit)
        _rot(p.c, ct, -st)
        if _compat.v2():
            _spin_rot_z(p, pm['theta2'], False)
        p.c[0] += dx
        p.c[2] += dy


def _tmulti(p: _Particle, el, krad):
    _tmulti_ens([p], el, krad)


# ---------------------------------------------------------------------------
# ソレノイド領域 (tsol.f tsol, tsfrin, trots; tdrift.f tdrift_solenoid; tsolqu.f の bz≠0 分岐)
#   SOL (BOUND=1) で始まり SOL (BOUND=1) で終わる。中の要素には bzs [1/m] = charge·BZ[T]/(p0/c) が掛かる
#   (tfinit.f tfbzs)。中に置けるのは DRIFT, 直線 BEND, QUAD, MULT, CAVI (ここでは DRIFT/QUAD のみ対応)。
#   CHI1/2/3=DX=DY=DZ=0 なら trots は恒等。F1=0 なら SOL の端で放射しない (tradks を呼ばない)。
# ---------------------------------------------------------------------------

class SolRegionError(NotImplementedError):
    """Raised when an element cannot be placed inside a SOL region (between SOL elements with BOUND nonzero).
    SOL 領域 (BOUND≠0 の SOL の間) に置けない要素がある (v2: K64 の tsol.f が出すはずの文で止める。K64 は SIGSEGV、docs/30 B35)"""


_SOL_OK = ('Drift', 'Bend', 'Quadrupole', 'Solenoid', 'Multipole', 'Cavity', 'Map', 'Marker', 'Monitor')


def _sol_check_type(el):
    """v2 (B35): tsol.f の select case が受け付ける型 (DRIFT BEND QUAD SOL MULT CAVI MAP APERT MARK MON) 以外は明確なエラー。
    k64 モードでは何もしない (今までの NotImplementedError・L=0 の素通し)"""
    from . import compat
    if compat.v2() and type(el).__name__ not in _SOL_OK:
        raise SolRegionError('Only DRIFT, BEND, QUAD, SOL, MULT, CAVI, MAP, APERT, MARK, MON are supported within SOL: '
                             f'{el.name} ({type(el).__name__})')


def _tfbzs(el, ph):
    """tfinit.f tfbzs: bzs = charge·(BZ+DBZ)·dir/(amass·(p0/amass)/c) = charge·BZ[T]/(p0[eV]/c) [1/m]。
    pysad の _Phys では brhoz = pgev/c (pgev は eV)。|bzs|<1e-20 は 0"""
    bz = (_g(el, 'BZ') + _g(el, 'DBZ')) * ph.charge / ph.brhoz
    return 0.0 if abs(bz) < 1e-20 else bz


def _tsfrin(p: _Particle, dbz):
    """tsol.f tsfrin: ソレノイドの端のフリンジ (bz の変化 dbz)"""
    if dbz == 0.0:
        return
    c = p.c
    x0, px0, y0, py0 = c[0], c[1], c[2], c[3]
    bq = -0.125 * dbz / (1.0 + c[5])
    bp = -2.0 * bq
    pr = x0 * px0 + y0 * py0
    pphi = x0 * py0 - y0 * px0
    phi0 = bq * pr
    w = math.exp(bq * pphi)
    sn, xs, cs, dc = _xsincos(phi0)
    c[0] = (cs * x0 - sn * y0) * w
    c[2] = (cs * y0 + sn * x0) * w
    c[1] = (cs * px0 - sn * py0) / w
    c[3] = (cs * py0 + sn * px0) / w
    c[4] = c[4] + bp * pr * pphi


def _tdrift_solenoid(p: _Particle, al, bz, enarad):
    """tdrift.f tdrift_solenoid"""
    c = p.c
    if enarad and _POL:
        p.bsi = 0.0
    pr = 1.0 + c[5]
    bzp = bz / pr
    pxi = c[1] + bzp * c[2] * 0.5
    pyi = c[3] - bzp * c[0] * 0.5
    zi = c[4]
    dpzi = _pxy2dpz(pxi, pyi)
    pzi = 1.0 + dpzi
    al1 = al / pzi
    phi = al1 * bzp
    a24, a12, a22, a14 = _xsincos(phi)
    c[0] += (a24 * pxi - a14 * pyi) / bzp
    c[2] += (a14 * pxi + a24 * pyi) / bzp
    c[4] += dpzi * al1 - p.dv * al
    px1 = a22 * pxi + a24 * pyi
    py1 = -a24 * pxi + a22 * pyi
    if enarad:
        # tradk1 (cphi0=1, sphi0=0) で基準は (pxi, pyi, zi)
        if _POL:
            p.bsi = bz * al
        c[1], c[3] = px1, py1
        p.pxr0, p.pyr0, p.zr0 = pxi, pyi, zi
        if _PH is not None:
            _ph_set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, al)          # tdrift.f: tsetpcvt (D1)
        _tradk_core(p, al, 0.0)
        px1, py1 = c[1], c[3]
    pr = 1.0 + c[5]
    bzp = bz / pr
    c[1] = px1 - bzp * c[2] * 0.5
    c[3] = py1 + bzp * c[0] * 0.5


def _ttfrin4_bz(p: _Particle, ak, al, bz):
    """tquad.f ttfrin (nord=4, bz≠0)"""
    if al == 0.0 or ak == 0.0:
        return
    c = p.c
    pr = 1.0 + c[5]
    aki = ak / al / 4.0 / pr
    bzph = 0.5 * bz / pr
    xi, yi = c[0], c[2]
    a = aki * xi**2
    b = aki * yi**2
    ab = a - b
    t = ab**2 / 6.0
    dx1 = xi * (a / 3.0 + b + t)
    dy1 = -yi * (a + b / 3.0 - t)
    c[0] = xi + dx1
    c[2] = yi + dy1
    d = a + b
    xx = 1.0 + d + ab * (5.0 * a - b) / 6.0
    yy = 1.0 - d + ab * (a - 5.0 * b) / 6.0
    h = 2.0 * aki * xi * yi * (1.0 - ab / 3.0)
    f = xx * yy + h**2
    if bz == 0.0:
        px1 = (c[1] * yy + c[3] * h) / f
        py1 = (c[3] * xx - c[1] * h) / f
        c[4] = c[4] - px1 * (dx1 + t * xi) - py1 * (dy1 + t * yi)
        c[1] = px1
        c[3] = py1
        return
    pxa = c[1] + bzph * c[2]
    pya = c[3] - bzph * c[0]
    c[1] = (pxa * yy + pya * h) / f - bzph * yi
    c[3] = (pya * xx - pxa * h) / f + bzph * xi
    c[4] = c[4] - c[1] * (dx1 + t * xi) - c[3] * (dy1 + t * yi)


class _TZ:
    """tsolque.f tzparams (bz≠0): tzsetparams"""

    def __init__(self, dp, aln0, akk, bz):
        self.aln = aln0
        pr = 1.0 + dp
        self.pr = pr
        bzp = bz / pr
        self.bzp = bzp
        akkp = akk / pr
        self.akkp = akkp
        w1 = math.sqrt((bzp**2 + math.hypot(bzp**2, 2.0 * akkp)) * 0.5)
        w2 = akkp / w1
        self.w1, self.w2 = w1, w2
        self.wss = 1.0 / (w1**2 + w2**2)
        self.ws = w1 + w2
        self.w12 = w1 - w2
        self.wd = bzp / self.ws
        self.phi1 = aln0 * w1
        self.phi2 = min(100.0, max(-100.0, aln0 * w2))
        self.s1, self.xs1, self.c1, self.dc1 = _xsincos(self.phi1)
        self.sh2, self.xsh2, self.ch2, self.dch2 = _xsincosh(self.phi2)
        self.g1 = (self.s1 * w2)**2 / akkp
        self.g2 = -(self.sh2 * w1)**2 / akkp
        self.wr1 = w1 / w2
        self.wr2 = w2 / w1
        self.cxs1 = self.dc1 * self.s1 - self.xs1
        self.cxs2 = self.dch2 * self.sh2 - self.xsh2
        self.aw1 = aln0 + (self.cxs2 * self.wr1) / self.ws
        self.aw2 = aln0 + (self.cxs1 * self.wr2) / self.ws
        self.cr2 = self.c1 * self.wr2
        self.cr3 = self.ch2 * self.wr1
        self.csw1 = akkp * aln0 * (-self.dc1 + self.dch2)
        self.csws = (self.ch2 * w1 * self.phi1 + self.c1 * w2 * self.phi2) * self.wss
        self.ca1 = -(self.ch2 * self.dc1) - self.dch2 - (self.s1 * self.sh2 * self.w12) / self.ws
        self.dxs = w2 * self.xs1 - w1 * self.xsh2
        self.dcw1 = self.dc1 * self.phi2
        self.dcw2 = self.dch2 * self.phi1

    def tztafs(self, zwu, pxt, pyt, zw1, wst, w12t, swss, g1t):
        """tsolque.f tztafs → dz"""
        w12xsh = w12t * self.xsh2
        wsxs = wst * self.xs1
        cwsd = self.c1 * w12xsh - self.ch2 * wsxs
        cswa = (self.csw1 + cwsd) * swss
        wsm = 4.0 * swss * self.akkp
        ws1 = (-1.0 + wsm) * wst
        ws2 = w12t * (1.0 + wsm)
        cr1 = self.bzp**2 * swss
        f1 = cswa + self.csws - zw1
        f2 = self.ca1 * cr1 + g1t
        f3 = (2.0 * cr1 * cwsd + self.cr3 * w12xsh
              + self.dcw1 * ws1 + self.dcw2 * ws2 - self.cr2 * wsxs)
        return zwu * (zwu * f3 / 2.0 + f2 * pxt + self.bzp * f1 * pyt) - (zw1 * pyt**2) / 2.0


def _tsolqu_bz_ens(parts, al, ak, bz, krad, eps0=0.0, fringe=True):
    """tsolqu.f tsolqu / tsolqur の bz≠0 分岐 (ak0x=ak0y=0)。ak = |K1|L > 0。全粒子をまとめて。
    放射あり (tsolqur) はスライス外側・粒子内側で、各スライスの後に全粒子 tradk。
    最後の半ステップ: tsolqur では最後に tzsetparams を呼んだ粒子 (最後のスライスの最後の粒子) の
    bzp を全粒子に使う (粒子の並びで結果が変わる SAD の挙動)。tsolqu は粒子ごとの値。"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    aka = math.hypot(ak, bz)
    ndiv = 1 + int(abs(al) * aka / eps)
    if krad:
        ndiv = min(1000, max(ndiv, _ndivrad(0.0, ak, bz, eps0, parts[0].ph.h0)))
        if _POL and _compat.v2():
            # v2 (P4): スピンの歳差は 1 刻みの前後の向きの変化 (弦) から出すので、らせんの回転角が大きいと不正確 (bz = 9 1/m で 9e-5)。
            # bz·al/ndiv ≤ 0.02 に刻む (領域内の DRIFT・BEND・Emittance[] の solemit._tquade と同じ)
            ndiv = min(1000, max(ndiv, 1 + int(abs(bz * al) / 0.02)))
    aln0 = al / ndiv
    akk = ak / al
    for p in parts:
        c = p.c
        bzp = bz / (1.0 + c[5])
        c[1] += bzp * c[2] * 0.5
        c[3] -= bzp * c[0] * 0.5
    tzs = [None] * len(parts)
    psi0_ = [0.0] * len(parts)
    for n in range(ndiv):
        for ip, p in enumerate(parts):
            c = p.c
            tz = _TZ(c[5], aln0, akk, bz)
            tzs[ip] = tz
            bzp = tz.bzp
            if krad and _POL:
                alr_ = aln0 * 0.5 if n == 0 else aln0
                p.bsi = (akk * c[0] * c[2] + bzp * alr_) if n == 0 else bzp * alr_
                if _compat.v2():
                    # v2 (C2): 縦の磁場の積分は 1 スライスに bz·aln0/(1+δ)、入口の端で + k x y/(1+δ) (出口で −)。
                    # K64 は最初のスライスを半分、最後のスライスを半分 (上書き) にするので、全体で 1 スライス分足りない
                    # v2 (P4): 端の ψ (硬い端の縦磁場) は端のフリンジを入れるとき (DISFRIN=0) だけ。DISFRIN=1 は縦磁場のない模型
                    p.bpar = bzp * aln0 + (akk * c[0] * c[2] / (1.0 + c[5]) if (n == 0 and fringe) else 0.0)
                    psi0_[ip] = akk * c[0] * c[2]
            ap = c[1]**2 + c[3]**2
            dpz = _sqrt1(-ap)
            r = -dpz / (1.0 + dpz) * aln0
            if n == 0:
                r *= 0.5
            a24, a12, a22, a14 = _xsincos(r * bzp)
            pxi = c[1]
            c[0] += (a24 * pxi - a14 * c[3]) / bzp
            c[2] += (a14 * pxi + a24 * c[3]) / bzp
            if _compat.v2():
                # v2 (B333): 運動量の回転は更新前の px で (行列式 1)。K64 (下) は更新後の pxi を使い、面積が cos²(wℓ) 倍になる
                pyi = -a24 * pxi + a22 * c[3]
                pxi = a22 * pxi + a24 * c[3]
            else:
                pxi = a22 * pxi + a24 * c[3]
                pyi = -a24 * pxi + a22 * c[3]          # SAD のまま: 更新後の pxi を使う
            c[4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
            xi, yi = c[0], c[2]
            w1, w2, ws, wss, wd = tz.w1, tz.w2, tz.ws, tz.wss, tz.wd
            a = (w2 * ws * xi - bzp * pyi) * wss
            bw = (ws * pxi - bzp * w2 * yi) * wss
            cc = (w1 * wd * xi + pyi) * wss
            dw = (-wd * pxi + w1 * yi) * wss
            u1w = a * tz.dc1 + bw * tz.s1
            u2w = -a * tz.s1 + bw * tz.dc1
            v1w = cc * tz.dch2 + dw * tz.sh2
            v2w = cc * tz.sh2 + dw * tz.dch2
            c[0] = c[0] + u1w + v1w * bzp
            c[1] = pxi + u2w * w1 + v2w * w2 * bzp
            c[2] = c[2] + wd * u2w + ws * v2w
            c[3] = pyi - wd * u1w * w1 + ws * v1w * w2
            awu = a / ws * w1
            dwu = dw * w2
            dz1 = tz.tztafs(awu, pxi, pyi, tz.aw1, ws, tz.w12, wss, tz.g1)
            dz2 = tz.tztafs(-dwu, -pyi, pxi, tz.aw2, -tz.w12, ws, -wss, tz.g2)
            c[4] = (c[4] + bzp * (-((awu * dwu * tz.dxs**2) / tz.akkp) + tz.ca1 * pxi * pyi * wss)
                    + dz1 + dz2 - aln0 * p.dv)
            if krad and _POL and _compat.v2() and _compat.spin_bpar():
                # v2 (B232): 四極の本体の運動量の向きの磁場 ∫B·dr = Δψ (ψ = k x y) を足す (K64 の bz·aln0 は螺旋の縦磁場だけで、本体の横磁場の B∥ が無い。
                # 端の ψ (DISFRIN=0) と本体の Δψ の和が 0 (Maxwell 場) になる。DISFRIN=1 は端なしの模型で本体の Δψ だけ残る)
                p.bpar += (akk * c[0] * c[2] - psi0_[ip]) / (1.0 + c[5])
            if krad and _POL and n == ndiv - 1:
                p.bsi = -akk * c[0] * c[2] + bzp * aln0 * 0.5
                if _compat.v2() and fringe:
                    p.bpar -= akk * c[0] * c[2] / (1.0 + c[5])
        if krad:
            for p in parts:
                _tradk(p, aln0, 0.0)
                p.bpar = None
            if _PH is not None and _PH.get('ctx') is not None and _compat.v2():
                _PH['ctx']['fr0'] += aln0 / _PH['ctx']['al']      # tsolqu.f: pcvt%fr0 += alr/al (v2 (P4): 領域内の四極の光子)
    if krad and _compat.v2():
        for p in parts:
            _tsolqu_bz_final(p, aln0, bz / (1.0 + p.c[5]))      # v2 (B08): 粒子ごとの bzp
    elif krad:
        for p in parts:
            _tsolqu_bz_final(p, aln0, tzs[-1].bzp)
    else:
        for ip, p in enumerate(parts):
            _tsolqu_bz_final(p, aln0, tzs[ip].bzp)


def _tsolqu_bz_final(p: _Particle, aln0, bzp):
    """tsolqu/tsolqur の bz≠0 分岐の最後の半ステップ。
    SAD の tsolqur (放射あり) では bzp は tz%tz0%bzp、すなわち最後に tzsetparams を呼んだ粒子
    (最後のスライスの最後の粒子) の値がそのまま全粒子に使われる (粒子の順番で結果が変わる SAD の挙動)。
    放射なしの tsolqu では粒子ごとの値。"""
    c = p.c
    ap = c[1]**2 + c[3]**2
    dpz = _sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0 * 0.5
    a24, a12, a22, a14 = _xsincos(r * bzp)
    pxi = c[1]
    c[0] += (a24 * pxi - a14 * c[3]) / bzp
    c[2] += (a14 * pxi + a24 * c[3]) / bzp
    c[1] = a22 * pxi + a24 * c[3] - bzp * c[2] * 0.5
    c[3] = -a24 * pxi + a22 * c[3] + bzp * c[0] * 0.5
    c[4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r


def _tquad_sol(parts, el, bzs, krad):
    """tsol.f の icQUAD → tquad (bz=bzs)。全粒子をまとめて進める。DX=DY=F1=FRINGE=0 のみ"""
    L = el.L
    K1L = _g(el, 'K1') * L
    if L == 0.0 or K1L == 0.0:
        raise NotImplementedError('track6: ソレノイド内の L=0 / K1=0 の四極')
    if any(_g(el, k) for k in ('DX', 'DY', 'F1', 'FRINGE')):
        raise NotImplementedError('track6: ソレノイド内の四極のずれ / F1 / FRINGE')
    ak = abs(K1L)
    th2 = _g(el, 'ROTATION') + (0.5 * math.pi if K1L < 0 else 0.0)
    fringe = _g(el, 'DISFRIN') == 0.0
    ct, st = math.cos(th2), math.sin(th2)
    if krad and _PH is not None and _compat.v2():
        _ph_set(0.0, 0.0, th2, 0.0, 0.0, 0.0, L, 0.0)        # tquad.f: tsetpcvt (v2 (P4): 領域内の四極の光子)
    for p in parts:
        _rot(p.c, ct, st)            # tsolrot (dx=dy=0): bzs=bz、回転のみ
        _spin_rot_z(p, th2, True)    # E: スピンも座標と一緒に回す (tsolrot の calpol。領域の外の _tquad_ens と同じ。K1<0 (th2 = ±π/2)・ROTATE≠0 で抜けていた)
        if krad:
            p.setr0()
        if fringe:
            _ttfrin4_bz(p, ak, L, bzs)
        if krad:
            p.setr0()
            if bzs != 0.0 and _compat.v2():
                # v2 (C2): 放射の基準の向きは力学的運動量 (track6n_ext.k_sol_quad と同じ)
                bzp_ = bzs / (1.0 + p.c[5])
                p.pxr0 = p.c[1] + bzp_ * p.c[2] * 0.5
                p.pyr0 = p.c[3] - bzp_ * p.c[0] * 0.5
    if bzs == 0.0:                   # tsolqu/tsolqur の bz=0 分岐
        if krad:
            _tsolqur_ens(parts, L, ak)
        else:
            for p in parts:
                _tsolqu(p, L, ak)
    else:
        _tsolqu_bz_ens(parts, L, ak, bzs, krad, fringe=fringe)
    for p in parts:
        if fringe:
            _ttfrin4_bz(p, -ak, L, bzs)
        _rot(p.c, ct, -st)
        _spin_rot_z(p, th2, False)


# ---------------------------------------------------------------------------
# v2 (C2): スピン追跡 (POL + RAD) の経路での、SOL 領域の端のずれ・回転と F1、領域内の BEND・MULT・CAVI、UND、PHSROT。
#   写像は高速追跡 (track6n・track6n_ext・und・phsrot) のカーネルを 1 回分ずつ呼び、放射とスピンの歳差 (sprot: Thomas-BMT) は
#   この module の tradk (Python) で行う。高速追跡と同じ軌道になる (FLUC OFF)。
# ---------------------------------------------------------------------------

def _state_of(parts):
    """粒子 → 高速追跡の状態 S (n, 11)。S[:, :6] は内部座標 (px は p で割った値)"""
    from . import track6n as _n
    ph = parts[0].ph
    S = _n.make_state(np.zeros((len(parts), 6)), ph)
    for k, p in enumerate(parts):
        S[k, :6] = p.c
        S[k, 6] = p.dv
        S[k, 7], S[k, 8], S[k, 9] = p.pxr0, p.pyr0, p.zr0
    return S


def _state_back(parts, S):
    for k, p in enumerate(parts):
        p.c[:] = S[k, :6]
        p.dv = float(S[k, 6])
        p.pxr0, p.pyr0, p.zr0 = float(S[k, 7]), float(S[k, 8]), float(S[k, 9])


def _kernel_rs():
    from . import track6n as _n
    return _n.make_rs(RNG.seed)


def _tradk_sol(p: _Particle, al, bz, bx=0.0, by=0.0, psi=0.0):
    """v2: SOL 領域の中の放射とスピンの歳差。p.c は力学的運動量、p.pxr0.. は刻みの前の力学的運動量。
    bz、(bx, by) はその刻みの磁場 (1/m、p0 で規格化)、psi は端の縦磁場の積分 (入口 +ψ、出口 −ψ。ψ = Im Σ k_n (x+iy)^(n+1)/(n+1)!)。
    運動量の向きの磁場の積分 ∫B∥dl/(Bρ) = [(bz pz + bx px + by py)·al/pz + psi]/(1+δ) を sprot に渡す"""
    if _POL:
        c = p.c
        pxm = 0.5 * (p.pxr0 + c[1])
        pym = 0.5 * (p.pyr0 + c[3])
        pzm = 1.0 + _pxy2dpz(pxm, pym)
        p.bpar = (bz * al + (bx * pxm + by * pym) * al / pzm + psi) / (1.0 + c[5])
    _tradk_core(p, al, 0.0)
    p.bpar = None


def _spin_only(p: _Particle, px0, py0, px1, py1):
    """v2: 長さのない蹴り (向き (px0, py0) → (px1, py1)) のスピンの歳差だけ (放射なし)。垂直な磁場による回転 (1 + aγ)θ"""
    if not _POL or (px0 == px1 and py0 == py1):
        return
    dpz0 = _pxy2dpz(px0, py0)
    dpz = _pxy2dpz(px1, py1)
    dpx = px1 - px0
    dpy = py1 - py0
    ppx = py1 * dpz0 - dpz * py0 + dpy
    ppy = dpz * px0 - px1 * dpz0 - dpx
    ppz = px1 * py0 - py1 * px0
    pp = p.ph.p0 * (1.0 + p.c[5])
    keep = p.bpar
    p.bpar = 0.0
    _sprot(p, px0 + dpx * 0.5, py0 + dpy * 0.5, ppx, ppy, ppz, 0.0, _e._p2h(pp), pp, 0.0, 1.0, 0.0)
    p.bpar = keep


def _sol_drift_v2(p: _Particle, al, bz, krad):
    """v2: 領域内の DRIFT。写像と放射は tdrift_solenoid と同じ (高速追跡と同じ軌道)、スピンは _tradk_sol"""
    if not (krad and _POL) or bz == 0.0:
        if bz != 0.0:
            _tdrift_solenoid(p, al, bz, krad)
        else:
            _tdrift_free(p, al)
        return
    c = p.c
    bzp = bz / (1.0 + c[5])
    pxi = c[1] + bzp * c[2] * 0.5
    pyi = c[3] - bzp * c[0] * 0.5
    zi = c[4]
    _tdrift_solenoid(p, al, bz, False)
    c[1] += bzp * c[2] * 0.5
    c[3] -= bzp * c[0] * 0.5
    p.pxr0, p.pyr0, p.zr0 = pxi, pyi, zi
    _tradk_sol(p, al, bz)
    bzp = bz / (1.0 + c[5])
    c[1] -= bzp * c[2] * 0.5
    c[3] += bzp * c[0] * 0.5


def _sol_bend_v2(p: _Particle, al, bz, kx, ky, krad):
    """v2: 領域内の BEND の 1 刻み (一様な双極 By = kx/al、Bx = ky/al とソレノイド)。写像と放射は tdrift (高速追跡の tdrift_bz1 と同じ)"""
    if not (krad and _POL):
        _tdrift_kick(p, al, kx, ky, krad, bz)
        return
    c = p.c
    bzp = bz / (1.0 + c[5])
    pxi = c[1] + bzp * c[2] * 0.5
    pyi = c[3] - bzp * c[0] * 0.5
    zi = c[4]
    _tdrift_kick(p, al, kx, ky, False, bz)
    c[1] += bzp * c[2] * 0.5
    c[3] -= bzp * c[0] * 0.5
    p.pxr0, p.pyr0, p.zr0 = pxi, pyi, zi
    _tradk_sol(p, al, bz, ky / al, kx / al)
    bzp = bz / (1.0 + c[5])
    c[1] -= bzp * c[2] * 0.5
    c[3] += bzp * c[0] * 0.5


class _RadHook:
    """track6n_ext.sol_mult_full の放射を、この module の tradk (スピンつき) で行う"""

    def __init__(self, parts, prm=None):
        self.parts = parts
        self.first = True
        self.last = False
        self.akr = None
        self.fringe = bool(prm[9]) if prm is not None else True       # v2 (P4): 端の ψ は端のフリンジ (DISFRIN=0) があるときだけ
        self.al_tot = float(prm[0]) if prm is not None else 0.0
        if prm is not None and prm[0] != 0.0:
            self.akr = [a / prm[0] for a in prm[2][:prm[1] + 1]]       # K_n/L (回した系)

    def _field(self, x, y):
        """(bx, by, ψ): By + i Bx = Σ k_n z^n/n!、ψ = Im Σ k_n z^(n+1)/(n+1)!"""
        if self.akr is None:
            return 0.0, 0.0, 0.0
        z = complex(x, y)
        f = 0j
        g = 0j
        zn = 1.0 + 0j
        fact = 1.0
        for n, k in enumerate(self.akr):
            f += k * zn / fact
            zn *= z
            fact *= (n + 1)
            g += k * zn / fact
        return f.imag, f.real, g.imag

    def setr0(self, S, n, bz):
        if _POL:
            for p in self.parts:
                p.bsi = 0.0

    def tradk(self, S, n, al, bz):
        for k, p in enumerate(self.parts):
            c = S[k, :6]                      # 行の view
            keep = p.c
            p.c = c
            p.dv = float(S[k, 6])
            p.pxr0, p.pyr0, p.zr0 = float(S[k, 7]), float(S[k, 8]), float(S[k, 9])
            bzp = bz / (1.0 + c[5])
            c[1] += bzp * c[2] * 0.5
            c[3] -= bzp * c[0] * 0.5
            bx, by, psi = self._field(c[0], c[2]) if _POL else (0.0, 0.0, 0.0)
            _tradk_sol(p, al, bz, bx, by, ((psi if self.first else 0.0) - (psi if self.last else 0.0)) if self.fringe else 0.0)
            p.setr0()
            S[k, 7], S[k, 8], S[k, 9] = p.pxr0, p.pyr0, p.zr0
            bzp = bz / (1.0 + c[5])
            c[1] -= bzp * c[2] * 0.5
            c[3] += bzp * c[0] * 0.5
            S[k, 6] = p.dv
            p.c = keep
        self.first = False
        if self.al_tot:
            _ph_adv(al / self.al_tot)             # tmulti.f: pcvt%fr0 += al1/al (全粒子の tradk の後に 1 回)


def _trots(p: _Particle, fr, ent):
    """tsol.f trots (track6n_ext.trots1 と同じ式) + スピンの回転 (rad かつ calpol: s ← R s)"""
    r, dx, dy, dz = fr
    c = p.c
    pxi, pyi = c[1], c[3]
    pzi = 1.0 + _pxy2dpz(pxi, pyi)
    xi, yi = (c[0], c[2]) if ent else (c[0] - dx, c[2] - dy)
    xf = r[0, 0] * xi + r[0, 1] * yi
    yf = r[1, 0] * xi + r[1, 1] * yi
    zf = r[2, 0] * xi + r[2, 1] * yi
    pxf = r[0, 0] * pxi + r[0, 1] * pyi + r[0, 2] * pzi
    pyf = r[1, 0] * pxi + r[1, 1] * pyi + r[1, 2] * pzi
    pzf = r[2, 0] * pxi + r[2, 1] * pyi + r[2, 2] * pzi
    c[1], c[3] = pxf, pyf
    if ent:
        c[0] = xf - pxf / pzf * zf + dx
        c[2] = yf - pyf / pzf * zf + dy
    else:
        c[0] = xf - pxf / pzf * zf
        c[2] = yf - pyf / pzf * zf
    c[4] = c[4] + zf / pzf + dz - p.dv * dz
    if _POL:
        p.s = r @ p.s


def _sol_edge(p: _Particle, pr, rad):
    """SOL の端 (track6n_ext.k_sol_edge と同じ): trots、tsfrin、F1 の放射 (tradks。スピンの歳差も)"""
    rot, r, dx, dy, dz, dbz, bz0, bz1, f1, radok = pr
    c = p.c
    krad = rad and radok
    if rot == 1:
        _trots(p, (r, dx, dy, dz), True)
    if krad:
        bzp = 0.5 * bz0 / (1.0 + c[5])
        p.pxr0 = c[1] + bzp * c[2]
        p.pyr0 = c[3] - bzp * c[0]
        p.zr0 = c[4]
        if _POL:
            p.bsi = 0.0
    _tsfrin(p, dbz)
    if krad:
        bzp = 0.5 * bz1 / (1.0 + c[5])
        c[1] += bzp * c[2]
        c[3] -= bzp * c[0]
        _tradk_sol(p, f1, 0.0)
        bzp = 0.5 * bz1 / (1.0 + c[5])
        c[1] -= bzp * c[2]
        c[3] += bzp * c[0]
        p.setr0()
    if rot == -1:
        _trots(p, (r, dx, dy, dz), False)


def _tsol_region_v2(parts, elements, i, rad, trf0, rfsw):
    """v2: SOL 領域を track6n_ext.compile_sol_region の表 (高速追跡と同じ) に沿って進める。戻り値は終点の添字"""
    from . import track6n as _n
    from . import track6n_ext as _x
    ph = parts[0].ph

    class _Col:
        def __init__(self):
            self.items = []

        def _add(self, k, p, name, cnt, sm):
            self.items.append((k, p))
    col = _Col()
    ke = _x.compile_sol_region(col, elements, i, ph)
    base_l = _PH['l'] if _PH is not None else 0
    for j, (kind, pr) in enumerate(col.items):
        el = elements[i + j]
        ph = parts[0].ph
        pref = getattr(ph, 'p0_ref', None)
        if _compat.v2() and pref is not None and pref != ph.p0 and kind in (
                _x.K_SOL_FRIN, _x.K_SOL_EDGE, _x.K_SOL_DRIFT, _x.K_SOL_BEND, _x.K_SOL_MULT, _x.K_SOL_QUAD):
            # v2 (B238): 領域の上流の CAVI が TRPT で基準運動量を変えた。bz はコンパイルのときの基準運動量 p0_ref で規格化してあるので p0_ref/p0 を掛ける
            # (track6n_ext.step_kind と同じ。K64 の tfbzs(l) は要素 l の入口の運動量で規格化)
            sc = pref / ph.p0
            if kind == _x.K_SOL_FRIN:
                pr = (pr[0] * sc,)
            elif kind == _x.K_SOL_EDGE:
                pr = pr[:5] + (pr[5] * sc, pr[6] * sc, pr[7] * sc) + pr[8:]
            elif kind == _x.K_SOL_DRIFT:
                pr = (pr[0], pr[1] * sc, pr[2])
            elif kind == _x.K_SOL_BEND:
                pr = (pr[0], pr[1] * sc) + pr[2:]
            elif kind == _x.K_SOL_MULT:
                pr = (pr[0], pr[1] * sc)
            else:
                pr = pr[:5] + (pr[5] * sc,) + pr[6:]
        if _PH is not None:
            _PH['l'] = base_l + j                        # l_track (v2 (P4): 領域内の要素の光子)
            _PH['ctx'] = None
        if kind == _n.K_NONE:
            continue
        if kind == _x.K_SOL_FRIN:
            for p in parts:
                _tsfrin(p, pr[0])
        elif kind == _x.K_SOL_EDGE:
            for p in parts:
                _sol_edge(p, pr, rad)
        elif kind == _x.K_SOL_DRIFT:
            al, bz, radok = pr
            # v2 (B1): 放射とスピンの歳差は 1 刻みの前後の向きの変化 (弦) から出すので、らせんの回転角 bz·al が大きいと不正確になる
            # (領域内の BEND と同じ。bz·al/ndiv ≤ 0.02)。K64 は 1 刻み (k64 モードは今まで通り)
            nd = (1 + int(abs(bz * al) / 0.02)) if (rad and radok and bz != 0.0 and _compat.v2()) else 1
            if _PH is not None and rad and radok and al != 0.0:
                _ph_set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, al)        # tdrift.f: tsetpcvt (v2 (P4): 領域内の DRIFT の光子)
            for _ in range(nd):
                for p in parts:
                    _sol_drift_v2(p, al / nd, bz, rad and radok)
                _ph_adv(1.0 / nd)
        elif kind == _x.K_SOL_QUAD:
            _tquad_sol(parts, el, pr[5], rad and pr[6])
        elif kind == _x.K_SOL_BEND:
            al, bz, kx, ky, radok = pr
            nd = (1 + int(abs(bz * al) / 0.02)) if (rad and radok) else 1        # k_sol_bend と同じ刻み
            for _ in range(nd):
                for p in parts:
                    _sol_bend_v2(p, al / nd, bz, kx / nd, ky / nd, rad and radok)
                _ph_adv(1.0 / nd)
        elif kind == _x.K_SOL_MULT:
            prm, bz = pr
            S = _state_of(parts)
            th = math.atan2(prm[8], prm[7])
            for p in parts:
                _spin_rot_z(p, th, True)
            if _PH is not None and rad and prm[11] and prm[0] != 0.0:
                _ph_set(prm[5], prm[6], th, 0.0, 0.0, 0.0, prm[0])   # tmulti.f: tsetpcvt (v2 (P4): 領域内の MULT の光子)
            ph2_ = _x.sol_mult_full(S, len(parts), ph, prm, bz, trf0, bool(_TRPT and _compat.v2()), rfsw, rad, _FLUC, _kernel_rs(), hook=_RadHook(parts, prm))
            _state_back(parts, S)
            if ph2_ is not None:
                for p in parts:
                    p.ph = ph2_              # v2 (B238): TRPT の加速。以降の放射・スピンの歳差 (aγ)・bz の規格化は空洞の後の基準運動量
            for p in parts:
                _spin_rot_z(p, th, False)
        elif kind == _x.K_MULT:
            _tmulti_ens(parts, el, rad and _g(el, 'RAD') == 0.0 and el.L != 0.0)
        elif kind == _n.K_CAV or kind == _n.K_DRIFT:
            for p in parts:
                _step_el(p, el, trf0, rfsw, rad)
        elif kind == _n.K_APERT:
            raise NotImplementedError('track6: ソレノイド領域内の APERT (スピン追跡の経路) は未対応')
        else:
            raise NotImplementedError(f'track6: ソレノイド内の要素 {type(el).__name__} ({el.name}) は未対応')
        for p in parts:
            if not np.all(np.isfinite(p.c)):
                raise FloatingPointError(f'track6: {el.name} で発散: {p.c}')
    return ke


def _und_b(x, y, s, F0, G0, dph, Kz, Kx, Ky, Qx, Qy):
    """アンジュレータの磁場 (Bx, By, Bs) = rot A、A = (ff, gg, 0) (und._und_a のポテンシャル。1/m、p0 で規格化)"""
    sn, cs = math.sin(Kz * s), math.cos(Kz * s)
    snp, csp = math.sin(Kz * s + dph), math.cos(Kz * s + dph)
    fy = (Ky * F0 / Kz * math.cos(Kx * x) * math.sinh(Ky * y) * sn
          - Qy * Qy * G0 / Qx / Kz * math.sinh(Qx * x) * math.cos(Qy * y) * snp)
    gx = (-G0 * Qx / Kz * math.cos(Qy * y) * math.sinh(Qx * x) * snp
          + Kx * Kx * F0 / Ky / Kz * math.sinh(Ky * y) * math.cos(Kx * x) * sn)
    by = F0 * math.cos(Kx * x) * math.cosh(Ky * y) * cs - Qy * G0 / Qx * math.sinh(Qx * x) * math.sin(Qy * y) * csp
    bx = G0 * math.cos(Qy * y) * math.cosh(Qx * x) * csp - Kx * F0 / Ky * math.sinh(Ky * y) * math.sin(Kx * x) * cs
    return bx, by, gx - fy


def _und_v2(parts, el, rad):
    """v2: UND (und.py) を track6 の経路で。磁場があれば刻みごとに [写像 (und_seg1) → 放射とスピンの歳差 (tradk)]"""
    from . import und as _un
    ph = parts[0].ph
    pr = _un.compile_und(el, ph.pgev)
    S = _state_of(parts)
    n = len(parts)
    rs = _kernel_rs()
    if pr[1] == 0.0 and pr[2] == 0.0:
        _un.k_und(S, n, *pr, ph.p0, ph.h0, _compat.v2())
        _state_back(parts, S)
        return
    (nslice, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds) = pr
    krad = rad and _g(el, 'RAD') == 0.0
    if not krad:
        _un.k_und_field(S, n, *pr, ph.p0, ph.h0, False, False, rs, ph.anrad, ph.cuc, _compat.v2())
        _state_back(parts, S)
        return
    for k, p in enumerate(parts):
        c = S[k, :6]
        keep = p.c
        p.c = c
        if _POL:
            p.bsi = 0.0
        for j in range(nslice):
            # 刻みの前の向き (力学的運動量 / p)。入口 (j=0) は und_seg1 が正準運動量からポテンシャルを引く
            if j == 0:
                a0x, a0y = _un._und_a(c[0], c[2], 0.0, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
                pn = 1.0 + c[5]
                vx0, vy0 = c[1] - a0x / pn, c[3] - a0y / pn
                _spin_only(p, c[1], c[3], vx0, vy0)          # 入口のポテンシャルの飛び (PHI≠0 のとき) は薄い蹴り
            else:
                vx0, vy0 = c[1], c[3]
            z0 = c[4]
            aox, aoy = _un._und_a(c[0], c[2], (j + 1) * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
            _un.und_seg1(k, S, j, j + 1, nslice + 1, F0, G0, dph, Kz, Kx, Ky, Qx, Qy, ds, False, False, rs, ph.p0, ph.anrad, ph.cuc)
            # 高速追跡 (und_seg1 の rad の枝) と同じ: 刻みの後の向きは、刻みの前の位置のポテンシャルで作った (px − ax)/p
            anx, any_ = _un._und_a(c[0], c[2], (j + 1) * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
            pn = 1.0 + c[5]
            c[1] += (anx - aox) / pn
            c[3] += (any_ - aoy) / pn
            p.pxr0, p.pyr0, p.zr0 = vx0, vy0, z0
            p.dv = float(S[k, 6])
            if _POL:
                bx_, by_, bs_ = _und_b(c[0], c[2], (j + 0.5) * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
                _tradk_sol(p, ds, bs_, bx_, by_)       # 運動量の向きの磁場 (軸を外れると Bs≠0) も歳差に入れる
            else:
                _tradk_sol(p, ds, 0.0)
            S[k, 6] = p.dv
            if _compat.v2():
                c[4] -= p.dv * ds          # v2 (B334): 速度の差による z の更新 −dv·ds (刻みごと。放射で p が変わるので)
        # 出口: 力学的運動量の向き → 正準運動量 (外側。出口の位置のポテンシャルを足す: und_seg1 の j1 = nslice と同じ)
        a1x, a1y = _un._und_a(c[0], c[2], nslice * ds, F0, G0, dph, Kz, Kx, Ky, Qx, Qy)
        pn = 1.0 + c[5]
        _spin_only(p, c[1], c[3], c[1] + a1x / pn, c[3] + a1y / pn)      # 出口のポテンシャルの飛び (周期の途中で終わるとき) は薄い蹴り
        c[1] += a1x / pn
        c[3] += a1y / pn
        h1 = math.sqrt(1.0 + (pn * ph.p0)**2)
        S[k, 6] = -c[5] * (1.0 + pn) / h1 / (h1 + pn * ph.h0) + S[k, 10]
        p.c = keep
    _state_back(parts, S)
    for p in parts:
        p.setr0()


def _phsrot_v2(parts, el, rad):
    """v2: PHSROT (phsrot.py) を track6 の経路で。位相空間の写像 (減衰・拡散つき) だけで、スピンは動かさない
    (PHSROT は磁場を持たないので歳差の大きさが決まらない)"""
    from . import phsrot as _prt
    ph = parts[0].ph
    pr = _prt.compile_phsrot(el)
    S = _state_of(parts)
    rs = _kernel_rs()
    _prt.k_phsrot(S, len(parts), *pr, rad, _FLUC, rs, ph.p0, ph.h0)
    RNG.seed = int(rs[0])
    _state_back(parts, S)
    for p in parts:
        p.setr0()


def _tsol_region(parts, elements, i, rad, trf0=0.0, rfsw=False):
    """tsol.f tsol: elements[i] が BOUND≠0 の SOL (始点)。終点 (次の BOUND≠0 の SOL) まで全粒子を追跡し、
    終点の添字を返す。GEO/CHI/DX/DY/DZ=0 (trots 恒等)、F1=0 (tradks なし) のみ。
    v2 は _tsol_region_v2 (ずれ・回転・F1、領域内の BEND・MULT・CAVI も)"""
    if _compat.v2():
        return _tsol_region_v2(parts, elements, i, rad, trf0, rfsw)
    ph = parts[0].ph
    sb = elements[i]
    if any(_g(sb, k) for k in ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')):
        raise NotImplementedError('track6: SOL のずれ・回転 (trots)')
    if _g(sb, 'F1') != 0.0:
        raise NotImplementedError('track6: SOL の F1 (tradks)')
    ke = None
    for j in range(i + 1, len(elements)):
        e = elements[j]
        if type(e).__name__ == 'Solenoid' and _g(e, 'BOUND') != 0.0:
            ke = j
            break
    if ke is None:
        raise ValueError('track6: ソレノイドの終端 (SOL BOUND=1) がない')
    bzs = _tfbzs(sb, ph)
    if _g(sb, 'DISFRIN') == 0.0:
        for p in parts:
            _tsfrin(p, bzs)
    for l in range(i + 1, ke + 1):
        el = elements[l]
        name = type(el).__name__
        _sol_check_type(el)
        if _is_apert(el):
            raise NotImplementedError('track6: ソレノイド領域内の APERT (tsol の tapert1) は未対応')
        if isinstance(el, Drift):
            for p in parts:
                if bzs != 0.0:
                    _tdrift_solenoid(p, el.L, bzs, rad and _g(el, 'RAD') == 0.0)
                else:
                    _tdrift_free(p, el.L)
        elif isinstance(el, Quadrupole):
            _tquad_sol(parts, el, bzs, rad and _g(el, 'RAD') == 0.0 and el.L != 0.0)
        elif name == 'Solenoid':
            bz1 = 0.0 if l == ke else _tfbzs(el, ph)
            if any(_g(el, k) for k in ('DX', 'DY', 'DZ', 'CHI1', 'CHI2', 'CHI3')) or _g(el, 'F1') != 0.0:
                raise NotImplementedError('track6: SOL のずれ・回転・F1')
            if _g(el, 'DISFRIN') == 0.0:
                for p in parts:
                    _tsfrin(p, bz1 - bzs)
            bzs = bz1
        elif name == 'BeamBeam':
            raise NotImplementedError('track6: ソレノイド領域内の BEAMBEAM は未対応 (K64 の tsol が落ちる)')
        elif el.L == 0.0:
            continue
        else:
            raise NotImplementedError(f'track6: ソレノイド内の要素 {name} は未対応 (DRIFT/QUAD のみ)')
        for p in parts:
            if not np.all(np.isfinite(p.c)):
                raise FloatingPointError(f'track6: {el.name} で発散: {p.c}')
    return ke


# ---------------------------------------------------------------------------
# 一周 (tturn1) と TrackParticles
#   SAD と同じく要素を外側、粒子を内側に回す (tsolqur の粒子間の持ち越しを再現するため)。
# ---------------------------------------------------------------------------

def _step_el(p: _Particle, el, trf0, rfsw, rad):
    """粒子ごとに進める要素 (ソレノイド領域の外、BEND/QUAD/MULT 以外)"""
    name = type(el).__name__
    if isinstance(el, Drift):
        _tdrift_free(p, el.L)
    elif name in _THIN_KEY:
        _tthin_sext(p, el, rad and _g(el, 'RAD') == 0.0 and el.L != 0.0)
    elif name == 'Cavity':
        cdx, cdy, cth = _g(el, 'DX'), _g(el, 'DY'), _g(el, 'ROTATION')
        ct_, st_ = _tent(p.c, cdx, cdy, cth)               # tcav.f: inc/TENT.inc (DX・DY・ROTATE)
        _tcav(p, el, trf0, rfsw, rad)
        _texit(p.c, cdx, cdy, cth, ct_, st_)
    elif name == 'BeamBeam':
        _step_beambeam(p, el)
    elif name == 'Coord':
        raise NotImplementedError('track6: COORD は track6n (高速追跡) だけが扱う')
    elif name == 'TCavity':
        raise NotImplementedError('track6: TCAVI は track6n (高速追跡) だけが扱う')
    elif name == 'Undulator':
        if _compat.v2():
            _und_v2([p], el, rad)
            return
        raise NotImplementedError('track6: UND は track6n (高速追跡) だけが扱う')
    elif name == 'PhsRot':
        if _compat.v2():
            _phsrot_v2([p], el, rad)
            return
        raise NotImplementedError('track6: PHSROT は K64 で動作しない (docs/30 B38)')
    elif el.L == 0.0:
        return          # MARK など
    else:
        raise NotImplementedError(f'track6: 要素 {name} は未対応')


def _cav_vnom(el, ph):
    """tffs.f:1455: vnominal = VOLT/amass |charge| sin(−PHI sign(charge)) (TRPT のときだけ使う。DVOLT は含まない)"""
    return _g(el, 'VOLT') / ph.amass * abs(ph.charge) * math.sin(-_g(el, 'PHI') * (1.0 if ph.charge >= 0 else -1.0))


def _tcav_acc_ens(parts, el, trf0, rfsw):
    """TRPT で基準運動量を変える CAVI (tcav.f の vnominal≠0: ndiv 分割と空洞の後の p0・h0 の更新)。高速追跡 (track6n.cav_full) と同じ関数を使い、
    全粒子の座標を S (n, 11) にまとめて通す。戻ったら全粒子の ph を空洞の後の基準運動量にする (以降の放射・スピンの歳差 aγ はこの p0・h0)。
    D1: これまで track6 の経路 (POL + RAD) は TRPT の加速を扱えず SADError にしていた。SOL 領域内の加速と RF つき MULT は未対応"""
    from . import track6n as _n6
    ph = parts[0].ph
    al = el.L
    vc = _g(el, 'VOLT') + _g(el, 'DVOLT')
    harm = _g(el, 'HARM') or _g(el, '_harm')
    freq = _g(el, 'FREQ')
    w = 2.0 * math.pi * freq / CVELOC if harm == 0.0 else ph.omega0 * harm / CVELOC
    wi = 1.0 / w if w != 0.0 else 0.0
    phic = (_g(el, 'PHI') + _g(el, 'DPHI')) * ph.charge
    v = vc / ph.amass * abs(ph.charge)
    vnom = _cav_vnom(el, ph)
    qa = abs(ph.charge) / ph.amass
    vq = (_g(el, 'V1') * qa, 0.5 * _g(el, 'V20') * qa, _g(el, 'V11') * qa, 0.5 * _g(el, 'V02') * qa)
    dphis = 0.0 - phic                     # tcav.f:50: TRPT は phis = 0
    offset = math.sin(dphis)
    offset1 = 0.0
    tlim = 1e4 if rfsw else 0.0
    cdx, cdy, cth = _g(el, 'DX'), _g(el, 'DY'), _g(el, 'ROTATION')
    n = len(parts)
    S = np.zeros((n, 11))
    for i, p in enumerate(parts):
        S[i, :6] = p.c
        S[i, 6] = p.dv
        S[i, 7], S[i, 8], S[i, 9] = p.pxr0, p.pyr0, p.zr0
        S[i, 10] = ph.dvfs
    if cdx != 0.0 or cdy != 0.0 or cth != 0.0:
        _n6.tent_all(S, n, cdx, cdy, cth)
    ph2 = _n6.cav_full(S, n, ph, al, v, w, wi, dphis, offset, offset1, tlim, _g(el, 'DISFRIN') == 0.0, int(round(_g(el, 'FRINGE'))),
                       vnom, vq, None, None)
    if cdx != 0.0 or cdy != 0.0 or cth != 0.0:
        _n6.texit_all(S, n, cdx, cdy, cth)
    for i, p in enumerate(parts):
        p.c[:] = S[i, :6]
        p.dv = float(S[i, 6])
        p.ph = ph2
        p.setr0()


def _step_beambeam(p: _Particle, el):
    """BEAMBEAM を track6 (Python の追跡) で通す (D1: PHOTONS で BEAMBEAM と BEND を同じラティスに置けるように)。
    beambeam.apply (track6n と同じ関数) に 1 粒子ぶんの S (n, 11) と、乱数の状態 (SadRandom の種) を渡し、戻す。
    光子の記録の粒子番号は apply の中で 1 になるので p.idx + 1 に直す"""
    from . import beambeam as _bb
    ph = p.ph
    key = (ph.pgev, id(getattr(el, '_photon', None)))
    cache = getattr(el, '_bb_cache', None)
    if cache is None or cache[0] != key:
        cache = (key, _bb.setup(el.params, ph.pgev, getattr(el, '_photon', None)))
        el._bb_cache = cache
    S = np.zeros((1, 11))
    S[0, :6] = p.c
    S[0, 6] = p.dv
    S[0, 7], S[0, 8], S[0, 9] = p.pxr0, p.pyr0, p.zr0
    S[0, 10] = ph.dvfs
    rs = np.array([RNG.seed, 1 if _compat.v2() else 0], dtype=np.int64)
    sink = _bb.PHOTON_SINK
    n0 = len(sink) if sink is not None else 0
    _bb.apply(S, 1, cache[1], rs)
    p.c[:] = S[0, :6]
    p.dv = float(S[0, 6])
    RNG.seed = int(rs[0])
    if sink is not None:
        for r in sink[n0:]:
            r[10] = float(p.idx + 1)


def _is_apert(el):
    """parser が APERT を Drift (L=0) + _apert=True で作る (パラメータは AX AY DX DY DX1.. の属性)"""
    return isinstance(el, Drift) and getattr(el, '_apert', False)


def _apert_pass(parts, nal, el, zlost):
    """tapert1 + tapert (tapert.f の直訳。numba 版は track6n.k_apert_gen)。parts[:nal] が生きている粒子。
    落ちた粒子は末尾の生きている粒子と入れ替えて外す (SAD と同じ並べ替え)。新しい nal を返す。
    LOSSAMPL による limittrack と最初の要素・15 要素ごとの tapert は track6 にはない (APERT 要素だけ)"""
    f = lambda k: float(getattr(el, k, 0.0) or 0.0)
    ax, ay, dx, dy, theta = f('AX'), f('AY'), f('DX'), f('DY'), f('ROTATION')
    xl, xh = min(f('DX1'), f('DX2')), max(f('DX1'), f('DX2'))
    yl, yh = min(f('DY1'), f('DY2')), max(f('DY1'), f('DY2'))
    pxj, pyj, dpj = f('JDPX'), f('JDPY'), f('DP')
    if pxj != 0.0 or pyj != 0.0:
        xa = min(abs(xh), abs(xl))
        ya = min(abs(yh), abs(yl))
        ax1 = 1.0 / xa if xa != 0.0 else 0.0
        ay1 = 1.0 / ya if ya != 0.0 else 0.0
        for p in parts[:nal]:
            c = p.c
            if (ax1 * c[0]) ** 2 + (ay1 * c[2]) ** 2 > 1.0:
                phi1 = 2.0 * math.pi * RNG.tran()
                phi2 = 2.0 * math.pi * RNG.tran()
                c[0] = xa * math.cos(phi1)
                c[2] = ya * math.sin(phi1)
                c[1] = pxj * math.cos(phi2)
                c[3] = pyj * math.sin(phi2)
                c[5] = c[5] + dpj
        return nal
    eli = ax != 0.0 and ay != 0.0
    if eli and (xh == xl or yh == yl):
        xl1, xh1, yl1, yh1 = -1e100, 1e100, -1e100, 1e100
    else:
        xl1, xh1, yl1, yh1 = xl, xh, yl, yh
    ax1 = 1.0 / ax if ax != 0.0 else 0.0
    ay1 = 1.0 / ay if ay != 0.0 else 0.0
    cost, sint = (math.cos(theta), math.sin(theta)) if theta != 0.0 else (1.0, 0.0)
    drop = [False] * nal
    for i, p in enumerate(parts[:nal]):
        x, y, z = p.c[0], p.c[2], p.c[4]
        if theta != 0.0:
            xa = cost * (x - dx) - sint * (y - dy)
            ya = sint * (x - dx) + cost * (y - dy)
        else:
            xa, ya = x - dx, y - dy
        ok = xl1 < xa < xh1 and yl1 < ya < yh1 and abs(z) <= zlost
        if ok and eli:
            ok = (ax1 * xa) ** 2 + (ay1 * ya) ** 2 <= 1.0
        drop[i] = not ok            # NaN は比較が偽になるので落ちる (tapert.f のコメント)
    i, np_ = 0, nal
    while i < np_:
        if drop[i]:
            while i < np_ - 1 and drop[np_ - 1]:
                np_ -= 1
            if not drop[np_ - 1] and np_ - 1 != i:
                parts[i], parts[np_ - 1] = parts[np_ - 1], parts[i]
                drop[i], drop[np_ - 1] = False, True
            np_ -= 1
        i += 1
    return np_


def _tapert_loss(full, nal, alost, zlost, mark, nodrop=False):
    """tturn の tapert (limittrack + 矩形 ±LOSSAMPL、|z| ≤ LOSSDZ)。track6n.k_apert と同じ判定。落ちた粒子は末尾へ寄せ、
    p.lost = mark (要素の番号 l_track, 周) を付ける。新しい nal を返す"""
    xlim = alost * 3.0
    drop = [False] * nal
    for i, p in enumerate(full[:nal]):
        c = p.c
        for k in (0, 2):
            v = c[k]
            if v != v:
                c[k] = xlim
            elif abs(v) > xlim:
                c[k] = math.copysign(xlim, v)
        for k in (1, 3):
            v = c[k]
            if v != v:
                c[k] = 0.9999
            elif abs(v) > 0.9999:
                c[k] = math.copysign(0.9999, v)
        v = c[4]
        if v != v:
            c[4] = 1e10
        elif abs(v) > 1e10:
            c[4] = math.copysign(1e10, v)
        ok = -alost < c[0] < alost and -alost < c[2] < alost and abs(c[4]) <= zlost and c[5] == c[5]
        drop[i] = (not ok) and not nodrop      # nodrop: K64 の最初の tapert (l_track=0) は落とす印が 0 になり落ちない (B33)
        if drop[i]:
            p.lost = mark
    i, np_ = 0, nal
    while i < np_:
        if drop[i]:
            while i < np_ - 1 and drop[np_ - 1]:
                np_ -= 1
            if not drop[np_ - 1] and np_ - 1 != i:
                full[i], full[np_ - 1] = full[np_ - 1], full[i]
                drop[i], drop[np_ - 1] = False, True
            np_ -= 1
        i += 1
    return np_


def _step_coord_tcav_ens(parts, el, rfsw, rad, kturn):
    """COORD (coord.k_coord) と TCAVI (tcav.apply_tcav) を track6 (Python の追跡) で通す。track6n と同じ kernel に、全粒子の S (n, 11) を渡す
    (K64 は tcoord・ttcave とも光子 (tphrec) を記録しない。BEND などが同じラティスにあるときの PHOTONS 用)。
    TCAVI の RANKICK・RANPHASE は周回ごとに 1 回引く (全粒子に同じ値) ので、粒子ごとに分けず全粒子まとめて呼ぶ。dv は COORD で更新しない (tcoord.f)"""
    name = type(el).__name__
    ph = parts[0].ph
    n = len(parts)
    S = np.zeros((n, 11))
    for i, p in enumerate(parts):
        S[i, :6] = p.c
        S[i, 6] = p.dv
        S[i, 7], S[i, 8], S[i, 9] = p.pxr0, p.pyr0, p.zr0
        S[i, 10] = ph.dvfs
    if name == 'Coord':
        from . import coord as _cd
        _cd.k_coord(S, n, *_cd.params(el))
    else:
        from . import tcav as _tc
        from . import track6n as _n6
        pr = _tc.compile_tcav(el)
        if rfsw:
            rs = np.array([RNG.seed, 1 if _compat.v2() else 0], dtype=np.int64)
            _tc.apply_tcav(S, n, pr, ph, bool(rad), bool(_FLUC), rs, kturn)
            RNG.seed = int(rs[0])
        else:
            _n6.k_drift(S, n, pr['al'])
    for i, p in enumerate(parts):
        p.c[:] = S[i, :6]
        p.dv = float(S[i, 6])
        if name != 'Coord':
            p.setr0()


def _tturn1(parts, elements, trf0, rfsw, rad, nal=None, zlost=100.0, alost=None, l0=0, kturn=1, lt_prev=0, nodrop_first=False,
            links=None, links_off=0):
    """tturn.f tturn1: 要素を外側、粒子を内側に回す。放射ありの BEND/QUAD/MULT はスライスごとに
    全粒子を進めてから全粒子の tradk (SAD と同じ順で乱数を引く)。
    alost (LOSSAMPL) を渡すと、最初の要素の前・15 要素ごと・区間の終わりに tapert (LOSSAMPL・LOSSDZ) を行い、発散した粒子も
    例外にせず落とす (スピン追跡など track6 の経路の損失)。l0: 区間の最初の要素の番号 − 1、kturn: 周の番号、lt_prev: 直前の l_track"""
    if nal is None:
        nal = len(parts)
    full = parts                     # 落ちた粒子は末尾へ寄せる (parts[:nal] が生きている)
    skip_to = -1
    la = 0
    for i, el in enumerate(elements):
        if i <= skip_to:
            continue
        name = type(el).__name__
        if alost is not None:
            if la <= 0:
                nal = _tapert_loss(full, nal, alost, zlost, (l0 + i if i >= 1 else lt_prev, kturn),
                                   nodrop=nodrop_first and i == 0)
                if nal <= 0:
                    break
                la = 15
            la -= 1
        parts = full[:nal]
        if _PH is not None:
            _PH['l'] = l0 + i + 1                        # l_track
            _PH['ctx'] = None
        if _is_apert(el):
            n0 = nal
            nal = _apert_pass(full, nal, el, zlost)
            for p in full[nal:n0]:
                p.lost = (l0 + i + 1, kturn)
            if nal <= 0:
                break
            la = 15
            continue
        if name == 'Solenoid':
            # BOUND=0 の SOL が領域の外にあるとき: K64 の tsol は BOUND≠0 でなくても、この SOL を領域の始点として次の BOUND≠0 の SOL まで
            # 領域として追う (tsgeo の幾何・光学は BOUND≠0 だけを始点とみなす。gs の追跡は BOUND=1 の始点と全桁一致: F4)
            skip_to = _tsol_region(parts, elements, i, rad, trf0, rfsw)
            continue
        krad = rad and _g(el, 'RAD') == 0.0 and el.L != 0.0
        if isinstance(el, Bend):
            bf_, bl_ = (links or {}).get(links_off + i, (True, True))        # tfgeo の iprev・inext (隣り合う BEND では F1 の分割を入れない)
            _tbend_ens(parts, el, krad, bf_, bl_)
        elif isinstance(el, Quadrupole):
            _tquad_ens(parts, el, krad)
        elif name == 'Multipole':
            _tmulti_ens(parts, el, krad)
        elif name == 'Cavity' and _TRPT and _g(el, 'VOLT') != 0.0 and _cav_vnom(el, parts[0].ph) != 0.0:
            _tcav_acc_ens(parts, el, trf0, rfsw)          # D1: TRPT の加速 (基準運動量が変わる。スピンの aγ は p.ph の新しい値を使う)
        elif name in ('Coord', 'TCavity'):
            _step_coord_tcav_ens(parts, el, rfsw, rad, kturn)      # 光子 (PHOTONS) と COORD・TCAVI の併用 (track6n の kernel を全粒子にまとめて通す)
        else:
            for p in parts:
                _step_el(p, el, trf0, rfsw, rad)
        if alost is not None:
            if any(not np.all(np.isfinite(p.c)) for p in parts):
                nal = _tapert_loss(full, nal, alost, zlost, (l0 + i + 1, kturn))      # 発散した粒子はここで落とす
                la = 15
                if nal <= 0:
                    break
            continue
        for p in parts:
            if not np.all(np.isfinite(p.c)):
                raise FloatingPointError(f'track6: {el.name} で発散: {p.c}')
    if alost is not None and nal > 0:
        nal = _tapert_loss(full, nal, alost, zlost, (l0 + len(elements), kturn))
    return nal


def track_particles(elements: List, coords, nturn: int, momentum_eV: float,
                    trf0: float = 0.0, rfsw: bool = True, rad: bool = False, fluc: bool = False,
                    seed=None, charge: float = 1.0, amass: float = None, progress: int = 0,
                    fshift: float = 0.0, circ: float = None, spins=None, radpol: bool = True,
                    alive=None, zlost: float = 100.0, alost=None, l0: int = 0, kturn0: int = 1, lt_prev: int = 0,
                    lost_info=None, nodrop_first: bool = False, links=None, links_off: int = 0, photons=None, photon_geo=None,
                    trpt: bool = False):
    """Track particles as `TrackParticles[{1, coords}, 1, 1, nturn]` does, from the lattice start for nturn turns.
    TrackParticles[{1, coords}, 1, 1, nturn] 相当 (始点 = 終点 = ラティスの先頭)。
    coords: (np, 6) の {x, px/p0, y, py/p0, z, δ}。nturn 周後の座標を返す。
    trf0: Emittance[] (emit) で決まった同期位置 (SAD の global trf0)。
    rad: RAD フラグ。fluc: FLUC フラグ (RAD と併用。乱数 RNG から引く)。seed: SeedRandom[seed] 相当
    (None なら RNG の現在の状態から続ける)。progress>0 なら progress 周ごとに経過を表示。
    SAD と同じく全粒子を要素ごとに進めるので、粒子の並びが結果に影響する場合 (ソレノイド内の
    放射あり四極) も、乱数を引く順番 (FLUC) も SAD と同じになる。"""
    global _FLUC, _POL, _RADPOL, _PH, _TRPT
    _TRPT = bool(trpt)
    if links is None and l0 == 0 and links_off == 0:
        from .track6n import bend_links as _bl, _bend_type as _bt
        links = _bl(elements) if any(_bt(e) for e in elements) else None
    if circ is None:                                        # 区間追跡 (TrackParticles の ls→ld) では全周長を渡す
        circ = sum(el.L for el in elements)
    kw = {} if amass is None else {'amass': amass}
    ph = _e._Phys(momentum_eV, circ, charge=charge, **kw)   # _Phys は eV 単位
    ph.dvfs = fshift if rfsw else 0.0                       # tsetdvfs: RFSW ON のとき FSHIFT (dv に足される = tconvm/tradk/tcav)
    if seed is not None:
        RNG.restore(seed)
    _FLUC = bool(rad and fluc)
    _PH = dict(sink=photons, geo=photon_geo, l=0, ctx=None) if (photons is not None and photon_geo is not None and rad and fluc) else None
    _POL = bool(spins is not None and rad)        # calpol: スピンは RAD が ON のときだけ動く (tradk の中で歳差)
    _RADPOL = bool(radpol)
    t0 = time.time()
    out = np.array(coords, dtype=float).copy()
    parts = [_Particle(out[ip], ph) for ip in range(out.shape[0])]
    for ip, p in enumerate(parts):
        p.idx = ip
    if spins is not None:
        for ip, p in enumerate(parts):
            p.s = np.array(spins[ip], dtype=float)
    # alive: 長さ np の bool 配列 (入出力)。APERT (tapert) で落ちた粒子は False になり、以降は追跡しない (座標は落ちた時点のまま)
    if alive is not None:
        parts.sort(key=lambda p: not alive[p.idx])        # 生きている粒子が前 (入力の順を保つ)
        nal = int(sum(1 for a in alive if a))
    else:
        nal = len(parts)
    try:
        for it in range(nturn):
            try:
                nal = _tturn1(parts, elements, trf0, rfsw, rad, nal, zlost, alost=alost, l0=l0, kturn=kturn0 + it,
                              lt_prev=(lt_prev if it == 0 else l0 + len(elements)),
                              nodrop_first=nodrop_first and it == 0, links=links, links_off=links_off)
            except FloatingPointError as e:
                raise FloatingPointError(f'{it + 1} 周目: {e}') from None
            if not rfsw:
                for p in parts[:nal]:
                    p.c[4] = 0.0          # tturn: rfsw=F なら z=0
            if progress and ((it + 1) % progress == 0 or it + 1 == nturn):
                for p in sorted(parts, key=lambda q: q.idx):
                    print(f'    [track6] P{p.idx + 1} turn {it + 1}/{nturn} ({time.time() - t0:.1f} s) '
                          + ' '.join(f'{v: .4e}' for v in p.out()), flush=True)
                sys.stdout.flush()
    finally:
        _FLUC = False
        _POL = False
        _PH = None
    for p in parts:
        out[p.idx] = p.out()
    if lost_info is not None:
        for p in parts[nal:]:
            lm = getattr(p, 'lost', None)
            if lm is not None and p.idx not in lost_info:
                lost_info[p.idx] = lm
    if alive is not None:
        for p in parts[nal:]:
            alive[p.idx] = False
    if spins is not None:
        sp = np.zeros((len(parts), 3))
        for p in parts:
            sp[p.idx] = p.s
        return out, sp
    return out


