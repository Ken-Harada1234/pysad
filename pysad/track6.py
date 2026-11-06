"""
track6.py — SAD TrackParticles (tftrack.f → tturn.f tturn / tturn1) の直訳による粒子追跡

SAD の追跡は Emittance[]/CALC (tturne: tquade, tbende, tthine, tcave, tradke) とは別のルーチンを使う:
  ドリフト  tdrift.f  tdrift_free
  四極      tquad.f   tquad → tsolrot (回転), ttfrin (nord=4), tsolqu / tsolqur (放射あり)
  六極      tquad.f   tthin (α/β/γ の 4 次積分) + ttfrin (nord=6 フリンジ)
  BEND      tbend.f   tbend → tbendcore (入口エッジ・本体・出口エッジを 1 つの閉じた式で)。
                      放射ありは tbendr で ndivrad 分割 (各分割で tbendcore + tradk)
  空洞      tcav.f    tcav (RAD の有無で offset を切り替え)
  放射      photons.f kradlib tradk → tradkn (FLUC OFF)。tradkfn (FLUC ON, 乱数) は未移植
内部座標は (x, px/(1+δ), y, py/(1+δ), z, δ) と粒子ごとの dv (tconvm)。入出力は px/p0。

未対応 (NotImplementedError): K1≠0 の BEND、BEND/四極/六極のずれ・回転 (四極 ROTATE は可)、
  長さのある空洞、四極の F1/FRINGE、ソレノイド、MULT など。FLUC (揺らぎ) は未対応。
"""
from __future__ import annotations

import math
import sys
import time
from typing import List

import numpy as np

from . import emit6d as _e
from . import lattice as _lat
from .elements import Drift, Bend, Quadrupole
from .sadrandom import SadRandom, tdusrn

CVELOC = 299792458.0
_GMIN = -0.999999        # photons.f module tspin の gmin (tradkn の max(gmin, ...))

# 乱数 (SAD 既定プラグイン "SAD")。FLUC の追跡はここから引く。track_particles(seed=...) で SeedRandom[seed] 相当
RNG = SadRandom(17)
_FLUC = False            # rfluct (RAD と FLUC が ON)


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


def _tradk(p: _Particle, al, phi0):
    """photons.f kradlib tradk → tradkn (rfluct=F, calpol=F)。最後に基準点を更新。"""
    _tradk_core(p, al, phi0)
    p.setr0()


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
        anp = ph.anrad * h1 * theta
        if _FLUC:
            dph, r1, r2, an = tdusrn(RNG, anp)
            if an == 0.0:
                return
            rho = al1 / theta
            uc = ph.cuc * h1**3 / ph.p0 / rho
            de = dph * uc
            dg = -de / (1.0 + 2.0 * de)
            c[5] = max(_GMIN, c[5] + dg)
            ddpx = -r1 * dpx * dg
            ddpy = -r1 * dpy * dg
            c[0] += r2 * ddpx * al1
            c[2] += r2 * ddpy * al1
        else:
            uc = ph.cuc * h1**3 / ph.p0 * theta / al1
            de = _e.CAVE * anp * uc
            dg = -de / (1.0 + 2.0 * de)
            c[5] = max(_GMIN, c[5] + dg)
            ddpx = -0.5 * dpx * dg
            ddpy = -0.5 * dpy * dg
            c[0] += ddpx * al1 / 3.0
            c[2] += ddpy * al1 / 3.0
        c[1] += ddpx
        c[3] += ddpy
        pr = 1.0 + c[5]
        p2 = ph.p0 * pr
        h2 = _e._p2h(p2)
        p.dv = -c[5] * (1.0 + pr) / h2 / (h2 + p2) + ph.dvfs       # SAD のまま (h2+p2)。+dvfs (tradkn)
        c[4] = c[4] * p2 / h2 * h1 / pp


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


def _tsolqu(p: _Particle, al, ak, eps0=0.0, ak0x=0.0, ak0y=0.0, qum=False):
    """tsolqu.f tsolqu / tsolqum (bz=0, 放射なし)。ak = |K1|L > 0。ak0x/ak0y は MULT の K0 (dx0 = ak0x/ak)。
    ak = 0 なら tdrift (ak0 付きの厳密ドリフト)。
    qum=True は tsolqum (tmulti の放射なしで使う版): z の dv の項が tsolqu と違い
    −0.25·(ap − dv)·aln0 (tsolqu は −0.25·ap·aln0 − dv·aln0)。SAD のまま (Level 11 で確認)。"""
    if ak == 0.0:
        _tdrift_kick(p, al, ak0x, ak0y, False)
        return
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    ndiv = 1 + int(abs(al * ak / eps))
    aln0 = al / ndiv
    dx0 = ak0x / ak
    dy0 = ak0y / ak
    akk = ak / al
    c = p.c
    w1, s1, xs1, dc1, sh2, xsh2, dch2 = _tzsetparam0(c[5], aln0, akk)
    x, px, y, py, z = c[0], c[1], c[2], c[3], c[4]
    dv = p.dv
    for n in range(ndiv):
        pxi, pyi = px, py
        ap = pxi**2 + pyi**2
        dpz = _sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * aln0
        if n == 0:
            r *= 0.5
        x += pxi * r
        y += pyi * r
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
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi**2 * xs1 - yi**2 * xsh2) + (ap - dv) * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2))
        else:
            z = (z - (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
                 - 0.25 * (w1 * (xi**2 * xs1 - yi**2 * xsh2) + ap * aln0
                           + u1 * px + xi * pxi * dc1 + v1 * py + yi * pyi * dch2)
                 - dv * aln0)
    ap = px**2 + py**2
    dpz = _sqrt1(-ap)
    r = -dpz / (1.0 + dpz) * aln0 * 0.5
    x += px * r
    y += py * r
    z -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
    c[0], c[1], c[2], c[3], c[4] = x, px, y, py, z


def _tdrift_kick(p: _Particle, al, ak0x, ak0y, enarad):
    """tdrift.f tdrift (bz=0, ak0x/ak0y ≠ 0): 一定の双極キック付きの厳密ドリフト (tsolconv)。
    内部座標 px/(1+δ)。放射は tradk1 (基準 px0, py0, z0 = 入口)。"""
    if ak0x == 0.0 and ak0y == 0.0:
        _tdrift_free(p, al)
        return
    c = p.c
    b = math.hypot(ak0x, ak0y)
    phix = ak0y / b
    phiy = ak0x / b
    phiz = 0.0
    dphizsq = phix**2 + phiy**2
    pr = 1.0 + c[5]
    alb = al * pr / b
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
        c[1], c[3] = px1, py1
        p.pxr0, p.pyr0, p.zr0 = px0, py0, z0
        _tradk_core(p, al, 0.0)
        px1, py1 = c[1], c[3]
    c[1] = px1
    c[3] = py1


def _ndivrad(phi, ak, bz, eps, h0, trpt=False):
    """tffs.f ndivrad"""
    eps1 = 1.0 if eps <= 0.0 else eps
    aka = (abs(phi) + math.hypot(ak, bz) * 0.02) / eps1
    return int(1.0 + h0 / 100.0 * aka * 10.0) if trpt else int(1.0 + h0 / 100.0 * aka)


def _tsolqur_ens(parts, al, ak, eps0=0.0):
    """tsolqu.f tsolqur (bz=0, ak0x=ak0y=0): スライス外側・粒子内側で、各スライスの後に全粒子 tradk
    (SAD と同じ順で乱数を引くため)"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    ndiv = 1 + int(abs(al) * ak / eps)
    ndiv = min(1000, max(ndiv, _ndivrad(0.0, ak, 0.0, eps0, parts[0].ph.h0)))
    aln = al / ndiv
    akk = ak / al
    for n in range(ndiv):
        alr = aln * 0.5 if n == 0 else aln
        for p in parts:
            c = p.c
            w1, s1, xs1, dc1, sh2, xsh2, dch2 = _tzsetparam0(c[5], aln, akk)
            ap = c[1]**2 + c[3]**2
            dpz = _sqrt1(-ap)
            r = -dpz / (1.0 + dpz) * alr
            c[0] += c[1] * r
            c[2] += c[3] * r
            c[4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r
            xi, yi = c[0], c[2]
            b, d = c[1], c[3]
            u1 = xi * dc1 + b * s1 / w1
            u2 = -xi * w1 * s1 + b * dc1
            v1 = yi * dch2 + d * sh2 / w1
            v2 = yi * w1 * sh2 + d * dch2
            c[0] += u1
            c[1] += u2
            c[2] += v1
            c[3] += v2
            c[4] = (c[4] - 0.25 * (w1 * (xi**2 * xs1 - yi**2 * xsh2) + (b**2 + d**2) * aln
                                   + u1 * (u2 + b) + xi * b * dc1 + v1 * (v2 + d) + yi * d * dch2)
                    - p.dv * aln)
        for p in parts:
            _tradk(p, aln, 0.0)
    for p in parts:
        c = p.c
        ap = c[1]**2 + c[3]**2
        dpz = _sqrt1(-ap)
        r = -dpz / (1.0 + dpz) * aln * 0.5
        c[0] += c[1] * r
        c[2] += c[3] * r
        c[4] -= (3.0 + dpz) * ap / 2.0 / (2.0 + dpz) * r


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
    f1 = −|akk·F1²|/24、f2 = |akk|·F2 (F1 = F1+F1K1F/B, F2 = F2+F2K1F/B)。逆向き要素は未対応。"""
    f1in = _g(el, 'F1') + _g(el, 'F1K1F')
    f1out = _g(el, 'F1') + _g(el, 'F1K1B')
    f2in = _g(el, 'F2') + _g(el, 'F2K1F')
    f2out = _g(el, 'F2') + _g(el, 'F2K1B')
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
        if krad:
            p.setr0()
        if fringe and mfring > -4 and mfring != 2:
            _ttfrin4(p, ak, L)
        if mfring in (1, 3):
            _tqlfr(p, f1in, f2in)
    if krad:
        if f1in != 0.0:
            for p in parts:
                _tradk(p, f1in, 0.0)     # SAD のまま: 長さの引数に係数 f1in を渡す
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
        for p in parts:
            _tradk(p, f1out, 0.0)
    for p in parts:
        c = p.c
        _rot(c, ct, -st)                 # tsolrot (exit): 回し戻す → ずらし戻す
        c[0] += dx
        c[2] += dy


def _tquad(p: _Particle, el, krad):
    _tquad_ens([p], el, krad)


# ---------------------------------------------------------------------------
# 六極 (tquad.f tthin, ttfrin nord=6)
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
               cosw, sinw, sqwh, sinwp1, krad, alr, bp):
    """tbend.f tbendcore。bp = (rhob, rho0, drhob, fb1, fb2) (tbendcom のモジュール変数)。"""
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
    if krad:
        _tradk(p, alr, phi0)


def _tbshift(p: _Particle, dx, dy, phi0, cost, sint, ent):
    """tbend.f tbshift"""
    c = p.c
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
        m.xsxkx = -xspx / akx / m.akxsq
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
    eps = 1e-3 if eps0 == 0.0 else 1e-3 * eps0          # epsbend
    ndiv = 1 + int(abs(phi0 / eps))
    krad = enarad and al != 0.0
    f1r = f2r = 0.0
    n1, n2 = 1, 0
    alc = al
    if krad:
        p.setr0()
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
    f = yi / m.rhoe
    fpx = af * c[1]
    m.pyi = c[3] - (tanp1 + fpx) * f
    ff = af * yi * f * 0.5
    m.xi = c[0] + ff
    m.zi = c[4] - ff * fpx
    m.pxi = c[1] + tanp1 * m.xi / m.rhoe
    m.yi = yi
    akx0 = alx0 = phix0 = 0.0
    for n in range(n1, n2 + 1):
        _tbendiinit(m, akxn[n], alx[n])
        _tbendicorr(m, (akxn[n] + akx0) * 0.5, (alx[n] + alx0) * 0.5, (phixn[n] + phix0) * 0.5)
        akx0, alx0, phix0 = akxn[n], alx[n], phixn[n]
        _tbendibody(m, alx0)
        if krad and n != n2:
            # tradk1 (cphi0 = cos(phixn(n)))。dp, dv(i) が更新される
            c[0], c[1], c[2], c[3], c[4], c[5] = m.xi, m.pxi, m.yi, m.pyi, m.zi, m.dp
            _tradk_core(p, alr[n], phixn[n])
            m.xi, m.pxi, m.yi, m.pyi, m.zi, m.dp = c[0], c[1], c[2], c[3], c[4], c[5]
            p.setr0()
    _tbendicorr(m, akx0 * 0.5, alx0 * 0.5, phix0 * 0.5)
    m.zi = m.zi - p.dv * al
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
    if _g(el, 'F1') != 0.0 or _g(el, 'FB1') != 0.0 or _g(el, 'FB2') != 0.0:
        raise NotImplementedError('track6: ANGLE=0 の BEND (tsteer) の F1/FB1/FB2 フリンジは未対応')
    phib = -ak0
    rhob = al / phib
    psi1 = _g(el, 'E1')
    tanp1 = math.sin(psi1) / math.cos(psi1)
    ctsn = []
    for p in parts:
        ctsn.append(_tent(p.c, dx, dy, theta))
    ndiv = 1
    if krad:
        for p in parts:
            p.setr0()
        ndiv = min(1000, max(1, _ndivrad(phib, 0.0, 0.0, _g(el, 'EPS'), parts[0].ph.h0)))
    aln = al / ndiv
    for n in range(1, ndiv + 1):
        for p in parts:
            c = p.c
            dp = c[5]
            pr = 1.0 + dp
            rhoe = rhob * pr
            pxi = c[1] - tanp1 * c[0] / rhoe
            c[3] += (tanp1 / rhoe) * c[2]
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
            dpx = aln / rhoe
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
            c[0] = xi + aln * (pxf + pxi) / (pz2 + pz1)
            c[2] += pyi * (aln + dl)
            dpv2 = _sqrt1(-pxf**2)
            pv2 = 1.0 + dpv2
            fa = c[2] / rhoe / pv2
            f = (1.0 - (c[2] / rhob)**2 / 6.0) * fa
            ff = (0.5 - (c[2] / rhob)**2 / 24.0) * c[2] * fa / pv2**2
            c[3] = pyi - pxf * f
            c[0] += ff
            c[1] = pxf
            c[4] -= (dl + p.dv * aln + ff * pxf)
        tanp1 = 0.0
        if krad:
            for p in parts:
                _tradk(p, aln, 0.0)
    for p, (cost, sint) in zip(parts, ctsn):
        _texit(p.c, dx, dy, theta, cost, sint)


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
        # SAD K64 の癖 (11_sad_quirks A35): tbend.f:573 は tquad に theta2=theta+dtheta を渡し、akang (K1<0 の π/2) を含めない。
        # tquad は ak=sign(ak0,al) で符号を theta2 に任せるので、K1<0 の BEND は |K1| の集束四極として追跡される
        # (SAD で QUAD K1=+0.5 と BEND ANGLE=0 K1=-0.5 の追跡が全桁同じ、CALC は正しく -0.5)。それを再現する
        return 'quad', SimpleNamespace(L=el.L, K1=abs(_g(el, 'K1')), DX=_g(el, 'DX'), DY=_g(el, 'DY'), ROTATION=th,
                                       DISFRIN=_g(el, 'DISFRIN'), FRINGE=0.0, F1=0.0, F2=0.0, F1K1F=0.0, F1K1B=0.0,
                                       F2K1F=0.0, F2K1B=0.0, EPS=_g(el, 'EPS'), RAD=_g(el, 'RAD'))
    m = _lat._bend_as_mult(el)
    m.ROTATION = th
    m.DROTATE = 0.0
    m.RAD = _g(el, 'RAD')
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
    if krad:
        for p in parts:
            p.setr0()
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
                       cosw, sinw, sqwh, sinwp1, krad, al, bp)
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
                       coswn, sinwn, -sqwhn, sinwpn, False, alr, bp)
        for p in parts:
            _tradk(p, alr, phi0n)


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
    dh = max(1.0 + 3.83e-12 - h1,
             (v + v20a * c[0]**2 + v02a * c[2]**2)
             * (-2.0 * math.sin(phh) * math.cos(phh - dphis) + offset) * wsn)
    h2 = h1 + dh
    a = max(dh * (h1 + h2) / p1**2, -1.0)
    dpr = _sqrt1(a)
    dp2r = dp1r + p1r * dpr
    p2r = 1.0 + dp2r
    c[5] = dp2r
    p.dv = -(1.0 + p2r) / h2 / (h2 + p2r * h0) * dp2r + p.ph.dvfs
    av = -(math.cos(2.0 * phh - dphis) * wi - offset1 * t) * wsn / p0
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
    dpt = vf * (math.sin(phw - dphis) + s0) / p1r
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


def _tmulti_ens(parts, el, krad):
    """tmulti.f tmulti (ANGLE=0, BZ=0, SK なし, RF なし, F1 なし)。K0/K1/K2… とずれ・回転 (tsolrot) に対応。
    tpara (tffs.f) の akr[n] = K_n·cr1^(n+1), theta2 = ROT+DROT+akang, ndiv0, dofr。
    各分割: tsolqu(al1, ak1, ak0) (+tradk 全粒子) → 多極キック (n≥2、(1+δ) で割る)。
    フリンジ: ttfrins で n=0 (双極), 1 (四極), 2 (六極) 等 (dofr)。"""
    pm = _lat._mult_params(el)
    al, nmmax, akr, ndiv = pm['al'], pm['nmmax'], pm['akr'], pm['ndiv']
    if _g(el, 'ANGLE') or _g(el, 'F1') or _g(el, 'FRINGE') or _g(el, 'VOLT'):
        raise NotImplementedError('track6: MULT の曲げ・F1・RF は未対応')
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
    akr1 = akr[1].real if nmmax >= 1 else 0.0
    akr0 = akr[0]
    eps = pm['eps']
    dofr = [akr[n] != 0 and (n <= 2 or abs(akr[n]) / al * 0.05**(n + 1) / (n + 1) > eps)
            for n in range(nmmax + 1)]
    fringe = pm['fringe']
    nz = al != 0.0
    if nz:
        if krad:
            for p in parts:
                p.setr0()
            ndiv = max(ndiv, _ndivrad(abs(akr0), akr1, 0.0, _g(el, 'EPS'), parts[0].ph.h0))
        ndiv = min(1000, ndiv)
        if fringe:
            for p in parts:
                for n in range(nmmax + 1):
                    if dofr[n]:
                        _ttfrins(p, n * 2 + 2, akr[n], al)
    wi = 1.0 / ndiv
    ak1 = akr1 * wi * 0.5
    al1 = al * wi * 0.5
    ak01 = akr0 * wi * 0.5
    akrm = [a * wi for a in akr]
    for m in range(1, ndiv + 1):
        if nz:
            for p in parts:
                _tsolqu(p, al1, ak1, 0.0, ak01.real, ak01.imag, qum=not krad)   # 放射なしは tsolqum
            if krad:
                for p in parts:
                    _tradk(p, al1, 0.0)
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
            _tsolqu(p, al1, ak1, 0.0, ak01.real, ak01.imag, qum=not krad)
        if fringe:
            for p in parts:
                for n in range(nmmax, -1, -1):
                    if dofr[n]:
                        _ttfrins(p, n * 2 + 2, -akr[n], al)
        if krad:
            for p in parts:
                _tradk(p, al1, 0.0)
    for p in parts:                      # tsolrot (exit)
        _rot(p.c, ct, -st)
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
        c[1], c[3] = px1, py1
        p.pxr0, p.pyr0, p.zr0 = pxi, pyi, zi
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


def _tsolqu_bz_ens(parts, al, ak, bz, krad, eps0=0.0):
    """tsolqu.f tsolqu / tsolqur の bz≠0 分岐 (ak0x=ak0y=0)。ak = |K1|L > 0。全粒子をまとめて。
    放射あり (tsolqur) はスライス外側・粒子内側で、各スライスの後に全粒子 tradk。
    最後の半ステップ: tsolqur では最後に tzsetparams を呼んだ粒子 (最後のスライスの最後の粒子) の
    bzp を全粒子に使う (粒子の並びで結果が変わる SAD の挙動)。tsolqu は粒子ごとの値。"""
    eps = 0.2 if eps0 == 0.0 else 0.2 * eps0
    aka = math.hypot(ak, bz)
    ndiv = 1 + int(abs(al) * aka / eps)
    if krad:
        ndiv = min(1000, max(ndiv, _ndivrad(0.0, ak, bz, eps0, parts[0].ph.h0)))
    aln0 = al / ndiv
    akk = ak / al
    for p in parts:
        c = p.c
        bzp = bz / (1.0 + c[5])
        c[1] += bzp * c[2] * 0.5
        c[3] -= bzp * c[0] * 0.5
    tzs = [None] * len(parts)
    for n in range(ndiv):
        for ip, p in enumerate(parts):
            c = p.c
            tz = _TZ(c[5], aln0, akk, bz)
            tzs[ip] = tz
            bzp = tz.bzp
            ap = c[1]**2 + c[3]**2
            dpz = _sqrt1(-ap)
            r = -dpz / (1.0 + dpz) * aln0
            if n == 0:
                r *= 0.5
            a24, a12, a22, a14 = _xsincos(r * bzp)
            pxi = c[1]
            c[0] += (a24 * pxi - a14 * c[3]) / bzp
            c[2] += (a14 * pxi + a24 * c[3]) / bzp
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
        if krad:
            for p in parts:
                _tradk(p, aln0, 0.0)
    if krad:
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
    for p in parts:
        _rot(p.c, ct, st)            # tsolrot (dx=dy=0): bzs=bz、回転のみ
        if krad:
            p.setr0()
        if fringe:
            _ttfrin4_bz(p, ak, L, bzs)
        if krad:
            p.setr0()
    if bzs == 0.0:                   # tsolqu/tsolqur の bz=0 分岐
        if krad:
            _tsolqur_ens(parts, L, ak)
        else:
            for p in parts:
                _tsolqu(p, L, ak)
    else:
        _tsolqu_bz_ens(parts, L, ak, bzs, krad)
    for p in parts:
        if fringe:
            _ttfrin4_bz(p, -ak, L, bzs)
        _rot(p.c, ct, -st)


def _tsol_region(parts, elements, i, rad):
    """tsol.f tsol: elements[i] が BOUND≠0 の SOL (始点)。終点 (次の BOUND≠0 の SOL) まで全粒子を追跡し、
    終点の添字を返す。GEO/CHI/DX/DY/DZ=0 (trots 恒等)、F1=0 (tradks なし) のみ。"""
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
        _tcav(p, el, trf0, rfsw, rad)
    elif el.L == 0.0:
        return          # MARK など
    else:
        raise NotImplementedError(f'track6: 要素 {name} は未対応')


def _tturn1(parts, elements, trf0, rfsw, rad):
    """tturn.f tturn1: 要素を外側、粒子を内側に回す。放射ありの BEND/QUAD/MULT はスライスごとに
    全粒子を進めてから全粒子の tradk (SAD と同じ順で乱数を引く)"""
    skip_to = -1
    for i, el in enumerate(elements):
        if i <= skip_to:
            continue
        name = type(el).__name__
        if name == 'Solenoid':
            if _g(el, 'BOUND') == 0.0:
                raise NotImplementedError('track6: BOUND=0 の SOL (領域外)')
            skip_to = _tsol_region(parts, elements, i, rad)
            continue
        krad = rad and _g(el, 'RAD') == 0.0 and el.L != 0.0
        if isinstance(el, Bend):
            _tbend_ens(parts, el, krad)
        elif isinstance(el, Quadrupole):
            _tquad_ens(parts, el, krad)
        elif name == 'Multipole':
            _tmulti_ens(parts, el, krad)
        else:
            for p in parts:
                _step_el(p, el, trf0, rfsw, rad)
        for p in parts:
            if not np.all(np.isfinite(p.c)):
                raise FloatingPointError(f'track6: {el.name} で発散: {p.c}')


def track_particles(elements: List, coords, nturn: int, momentum_eV: float,
                    trf0: float = 0.0, rfsw: bool = True, rad: bool = False, fluc: bool = False,
                    seed=None, charge: float = 1.0, amass: float = None, progress: int = 0,
                    fshift: float = 0.0, circ: float = None) -> np.ndarray:
    """TrackParticles[{1, coords}, 1, 1, nturn] 相当 (始点 = 終点 = ラティスの先頭)。
    coords: (np, 6) の {x, px/p0, y, py/p0, z, δ}。nturn 周後の座標を返す。
    trf0: Emittance[] (emit) で決まった同期位置 (SAD の global trf0)。
    rad: RAD フラグ。fluc: FLUC フラグ (RAD と併用。乱数 RNG から引く)。seed: SeedRandom[seed] 相当
    (None なら RNG の現在の状態から続ける)。progress>0 なら progress 周ごとに経過を表示。
    SAD と同じく全粒子を要素ごとに進めるので、粒子の並びが結果に影響する場合 (ソレノイド内の
    放射あり四極) も、乱数を引く順番 (FLUC) も SAD と同じになる。"""
    global _FLUC
    if circ is None:                                        # 区間追跡 (TrackParticles の ls→ld) では全周長を渡す
        circ = sum(el.L for el in elements)
    kw = {} if amass is None else {'amass': amass}
    ph = _e._Phys(momentum_eV, circ, charge=charge, **kw)   # _Phys は eV 単位
    ph.dvfs = fshift if rfsw else 0.0                       # tsetdvfs: RFSW ON のとき FSHIFT (dv に足される = tconvm/tradk/tcav)
    if seed is not None:
        RNG.restore(seed)
    _FLUC = bool(rad and fluc)
    t0 = time.time()
    out = np.array(coords, dtype=float).copy()
    parts = [_Particle(out[ip], ph) for ip in range(out.shape[0])]
    try:
        for it in range(nturn):
            try:
                _tturn1(parts, elements, trf0, rfsw, rad)
            except FloatingPointError as e:
                raise FloatingPointError(f'{it + 1} 周目: {e}') from None
            if not rfsw:
                for p in parts:
                    p.c[4] = 0.0          # tturn: rfsw=F なら z=0
            if progress and ((it + 1) % progress == 0 or it + 1 == nturn):
                for ip, p in enumerate(parts):
                    print(f'    [track6] P{ip + 1} turn {it + 1}/{nturn} ({time.time() - t0:.1f} s) '
                          + ' '.join(f'{v: .4e}' for v in p.out()), flush=True)
                sys.stdout.flush()
    finally:
        _FLUC = False
    for ip, p in enumerate(parts):
        out[ip] = p.out()
    return out


