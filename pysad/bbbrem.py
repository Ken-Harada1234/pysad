"""`BBBrem1`: generation of one beam-beam bremsstrahlung event (Kleiss-Burkhardt), translation of tbbbrem.f in K64.
BBBrem1[sqrt(s) (eV), ds/s min]: ビーム・ビーム制動放射 1 事象の生成 (Kleiss-Burkhardt)。K64 の tbbbrem.f (bcube1、tbcube1)。
戻り値は {散乱した陽電子の 4 元運動量 (eV, 符号反転と x y z の入れ替えなし: (-1,-1,-1,1)*p2), 近似断面積 (mbarn × 1e-31)}。
乱数は SAD の tran() (nran=3)。パッケージ bbbrem.n (BBBrem クラス) はこの関数を使う。
"""
from __future__ import annotations

import math

_MACPHYS_FINEST = 1.0 / 137.035999084
_ELMASS = 0.51099895000e6
_PLANKR = 6.62607015e-34 / (2.0 * math.pi)
_CVELOC = 299792458.0
_ELEMCH = 1.602176634e-19
_RME = _ELMASS / 1e9
_TOMB = (_PLANKR * _CVELOC / _ELEMCH) ** 2 / 1.0e18 / 1.0e-31
_RME2 = _RME ** 2
_TWOPI = 2.0 * math.pi


def bcube1(roots, rk0, random):
    """Compute the bcube1 event generation (roots in GeV), returning p2, qk, q2, weight, and sigapp.
    bcube1 (roots は GeV)。戻り: (p2[4], qk[4], q2[4], weight, sigapp)。random() は tran()"""
    alpha = _MACPHYS_FINEST
    rme2 = _RME2
    s = roots ** 2
    rme2s = rme2 / s
    rls = -math.log(rme2s)
    z0 = rk0 / (1 - rk0)
    a1 = math.log((1 + z0) / z0)
    a2 = math.log(1 + z0) / z0
    ac = a1 / (a1 + a2)
    sigapp = 8 * alpha ** 3 / rme2 * (-math.log(rme2s)) * (a1 + a2) * _TOMB
    eb = roots * 0.5
    pb = math.sqrt(eb * eb - rme2)
    rin2pb = 0.5 / pb
    q1 = [0.0, 0.0, pb, eb]
    # z
    if random() < ac:
        temp1 = random()
        z = 1.0 / (temp1 * (math.exp(a1 * random()) - 1))
    else:
        z = z0 / random()
    y = rme2s * z
    q0 = eb * y
    temp1 = pb * pb - eb * q0
    temp2 = temp1 * temp1 - rme2 * q0 * q0
    if temp2 < 0.0:
        return [0.0] * 4, [0.0] * 4, [0.0] * 4, 0.0, sigapp
    tmin = -2 * (temp1 + math.sqrt(temp2))
    tmax = rme2 * s * y * y / tmin
    sy = s * y
    w2 = sy + rme2
    temp1 = sy + tmax
    rlamx = math.sqrt(temp1 * temp1 - 4 * w2 * tmax)
    if temp1 <= 0.0:
        temp1 = rlamx - temp1
    else:
        temp1 = -4 * w2 * tmax / (rlamx + temp1)
    while True:
        b = math.exp(random() * math.log(1 + 2 * sy / temp1))
        t = -b * z * z * rme2 / ((b - 1) * (b * z + b - 1))
        if t >= tmin:
            break
    rlam = math.sqrt((sy - t) * (sy - t) - 4 * rme2 * t)
    eps = 4 * rme2 * w2 / (rlam * (rlam + w2 + rme2 - t))
    rl = math.log((2 + eps) / eps)
    vgam = eps * (math.exp(random() * rl) - 1)
    cgam = 1 - vgam
    sgam = math.sqrt(vgam * (2 - vgam))
    phi = _TWOPI * random()
    phig = _TWOPI * random()
    ql = (2 * eb * q0 - t) * rin2pb
    qt = math.sqrt((tmax - t) * (t - tmin)) * rin2pb
    q = [qt * math.sin(phi), qt * math.cos(phi), ql, q0]
    q2 = [q1[k] - q[k] for k in range(4)]
    r0 = eb + q0
    w = math.sqrt(w2)
    rin2w = 0.5 / w
    rinr0w = 1.0 / (r0 + w)
    eta = -(sy + 2 * w * q0 + t) * rin2w * rinr0w
    phat1 = -q[0] * (1 + eta)
    phat2 = -q[1] * (1 + eta)
    phat3 = pb * eta - ql * (1 + eta)
    phatl = rlam * rin2w
    phatt = math.sqrt(phat1 * phat1 + phat2 * phat2)
    sfhat = phat1 / phatt
    cfhat = phat2 / phatt
    sthat = phatt / phatl
    if phat3 > 0.0:
        vthat = sthat * sthat / (1 - math.sqrt(1 - sthat * sthat))
    else:
        vthat = sthat * sthat / (1 + math.sqrt(1 - sthat * sthat))
    cthat = vthat - 1
    sfg = math.sin(phig)
    cfg = math.cos(phig)
    temp1 = sgam * sfg
    temp2 = cthat * sgam * cfg + sthat * cgam
    veg = vthat + vgam - vthat * vgam - sthat * sgam * cfg
    qk4h = sy * rin2w
    qkh = [qk4h * (cfhat * temp1 + sfhat * temp2), qk4h * (-sfhat * temp1 + cfhat * temp2), qk4h * (veg - 1), qk4h]
    temp1 = pb * qkh[2]
    if temp1 > 0.0:
        temp2 = (rme2 * qkh[3] * qkh[3] + pb * pb * (qkh[0] * qkh[0] + qkh[1] * qkh[1])) / (eb * qkh[3] + temp1)
    else:
        temp2 = eb * qkh[3] - temp1
    qk = [0.0] * 4
    qk[3] = (temp2 + qkh[3] * q[3] + qkh[0] * q[0] + qkh[1] * q[1] + qkh[2] * q[2]) / w
    temp1 = (qk[3] + qkh[3]) * rinr0w
    qk[0] = qkh[0] + temp1 * q[0]
    qk[1] = qkh[1] + temp1 * q[1]
    qk[2] = qkh[2] + temp1 * (-pb + ql)
    p2 = [-q2[0] - qk[0], -q2[1] - qk[1], -q2[2] - qk[2], -q2[3] - qk[3] + roots]
    if qk[3] < eb * rk0:
        return p2, qk, q2, 0.0, sigapp
    c1 = math.log(1 + z) / math.log((2 + eps) / eps)
    temp1 = sy - tmax
    vnumx = math.sqrt(temp1 * temp1 - 4 * rme2 * tmax) + temp1
    temp1 = sy + tmax
    if temp1 < 0.0:
        vdenx = math.sqrt(temp1 * temp1 - 4 * w2 * tmax) - temp1
    else:
        vdenx = -4 * w2 * tmax / (math.sqrt(temp1 * temp1 - 4 * w2 * tmax) + temp1)
    temp1 = sy - tmin
    vnumn = math.sqrt(temp1 * temp1 - 4 * rme2 * tmin) + temp1
    temp1 = sy + tmin
    if temp1 < 0.0:
        vdenn = math.sqrt(temp1 * temp1 - 4 * w2 * tmin) - temp1
    else:
        vdenn = -4 * w2 * tmin / (math.sqrt(temp1 * temp1 - 4 * w2 * tmin) + temp1)
    c2 = 2 * rls / math.log((vnumx * vdenn) / (vdenx * vnumn))
    rlabl = (t - 2 * rme2 * y) * rin2pb
    rhat4 = -(2 * rme2 * y + (1 - y) * t) * rin2w
    etar = rhat4 * rinr0w
    rhat1 = -q[0] * (1 + etar)
    rhat2 = -q[1] * (1 + etar)
    rhat3 = rlabl + (pb - ql) * etar
    zz = s * (rhat4 * qkh[3] - rhat1 * qkh[0] - rhat2 * qkh[1] - rhat3 * qkh[2])
    s1 = 4 * eb * (eb - qk[3])
    d1 = sy * rlam * (eps + vgam) * rin2w * rin2w
    d2 = 0.5 * sy
    rind1 = 1.0 / d1
    rind12 = rind1 * rind1
    rind2 = 1.0 / d2
    rind22 = rind2 * rind2
    temp1 = s + t - 2 * d2
    temp2 = s1 + t + 2 * d1
    aa0 = (s * s + s1 * s1 + temp1 * temp1 + temp2 * temp2) * rind1 * rind2 * (-t)
    aa1 = -4 * rme2 * zz * zz * rind12 * rind22
    aa2 = -8 * rme2 * (d1 * d1 + d2 * d2) * rind1 * rind2
    aa3 = 16 * rme2 * rme2 * (d1 - d2) * zz * rind12 * rind22
    aa4 = -16 * rme2 * rme2 * rme2 * (d1 - d2) ** 2 * rind12 * rind22
    rmex = aa0 + aa1 + aa2 + aa3 + aa4
    rmap = 4 * s * s * rind1 * rind2 * (-t) * c1 * c2
    weight = rmex / rmap * sigapp
    return p2, qk, q2, weight, sigapp


def bbbrem1(ars, drs, random):
    """tbcube1: BBBrem1[sqrt(s) (eV), ds/s]。{(-1,-1,-1,1)*p2*GeV, sigapp*1e-31}"""
    p2, _qk, _q2, _w, sigapp = bcube1(ars / 1e9, drs, random)
    gev = 1.0e9
    return [-p2[0] * gev, -p2[1] * gev, -p2[2] * gev, p2[3] * gev], sigapp * 1.0e-31
