# -*- coding: utf-8 -*-
"""v2.5.0: (E) `tsolrot` の傾いた枝 (SOL 領域内の MULT・QUAD の CHI1・CHI2・DZ、Bz ≠ 0) の dcosphi の符号の不整合 (docs/30 B337)、
(F) MULT の CHI1・CHI2 の符号の確認。

(E) K64 `tsol.f:83-100` は `xsincos` の dcosphi = cos θ − 1 を、面の条件 (`s = plz*xsinphi + pz0*sinphi - pbz*dcosphi`) では正しい符号 (+ pbz (1 − cos θ)) で使い、
位置・運動量の更新 (`x += (plx*phi + ptx*sinphi + pbx*dcosphi)*alb`、`px -= ptx*dcosphi ...`) では 1 − cos θ の想定の符号で使う。同じらせん運動の式の中で符号が逆。
物理の基準 (独立): 長さ 0 の傾いた要素は恒等写像でなければならない。入口 (`ent`) と出口の写像を numpy に写し (下)、更新の符号を Fortran のまま (as) にすると恒等から 2e-7〜6e-6 ずれ、
反転 (fix) すると 1e-15 以内で恒等になる。
gs の保存した出力 (scripts/v250_tsolrot_tilt__{plain,tilt}.{sad,gs.txt}: 3 GeV/c、SOL (BZ = 13 T、bz = 1.3 /m、BOUND = 1) → MULT (L = 0、K2 = 1e-9。tilt は CHI1 = 0.03、CHI2 = −0.02、DZ = 0.01) → SOL (BZ = 0)、
4 粒子 (px は最大 0.05)): 傾きなし (plain) は入力に戻る (5e-16〜6e-15)。傾き (tilt) は **戻らず 2.2e-7、3.9e-7、1.4e-6、6.5e-6 ずれる**。numpy の as の恒等からのずれ (2.2e-7、3.9e-7、1.4e-6、6.5e-6) と
同じ大きさ。K64 の誤り (翻訳の誤りではない: pysad にはこの枝が無い)。pysad は未実装 (MULT の CHI1・CHI2・DZ、領域内の MULT) なので、v2 の修正は実装時に反映する (更新の dcosphi の符号を反転する)。

(F) gs の MULT (L = 1 m、K1 = 0.3) の軸上の粒子の出口の x・y・px・py: CHI1 = +0.05 で x = +1.1939e-3、px = −1.8396e-4、CHI2 = +0.05 で y = −1.3064e-3、py = −1.9141e-4 (符号を変えると符号が反転)。
ヘルプの 3×3 の式 (R₁、R₂) が与える要素の軸の向き (R の 3 行目: CHI1 では (−sin χ1, 0, cos χ1)、CHI2 では (0, −sin χ2, cos χ2)) の直線の周りの厚い四極 (x は収束、y は発散、K = 0.3 /m²) の ODE の
x = +1.1954e-3、px = −1.8393e-4、y = −1.3080e-3、py = −1.9144e-4 と 0.13% 以内で一致する。つまり K64 の MULT の CHI1・CHI2 の符号はヘルプの R₁・R₂ と同じ (`tmulti` が −chi を `tsolrot` に渡すのは内部の規約)。
残る 0.1% は、傾いた入口・出口の面の幾何の違い (K64 は面を傾けず座標を回して面まで直進させる)。
限界: (E) は 1 つの傾き (CHI1 = 0.03、CHI2 = −0.02、DZ = 0.01、領域の長さ 0) と bz = 1.3 /m。傾いた領域の中の MULT (L ≠ 0) の本体は見ていない。
(F) は L = 1 m の厚い四極 1 つ・CHI = ±0.05 で、DZ・CHI3 (ROTATE) の符号は見ていない。
"""
import math
import os
import re

import numpy as np
import pytest
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')
CL = 299792458.0
BRHO = 3e9 / CL
PTS = np.array([[1e-3, 2e-3, -5e-4, 1e-3, 0.0, 0.0], [3e-3, -4e-3, 2e-3, 5e-3, 0.0, 0.01],
                [-2e-3, 1e-2, 1e-3, -8e-3, 0.0, -0.02], [5e-4, 0.05, 5e-4, 0.03, 0.0, 0.0]])
CHI1, CHI2, DZ = 0.03, -0.02, 0.01


def parse(path):
    txt = open(os.path.join(SCRIPTS, path), encoding='utf-8').read()
    body = txt.split('R {', 1)[1]
    return np.array([[float(v.replace('D', 'E')) for v in n.split(',')] for n in re.findall(r'\{([^{}]*)\}', body)])


def xsincos(x):
    s = math.sin(x)
    return s, x - s, math.cos(x), math.cos(x) - 1.0


def pz_m1(px, py):
    s = px * px + py * py
    return -s / (1 + math.sqrt(1 - s))


def ent(c, chi1, chi2, dz, bz, alg, update_sign):
    """tsol.f tsolrot (ent、bz ≠ 0) の numpy への写し。update_sign = +1: 位置・運動量の更新で dcosphi = cos − 1 をそのまま使う (Fortran のまま)、−1: 1 − cos"""
    x, px, y, py, z, g = c
    schi1, _, cchi1, dcchi1 = xsincos(chi1)
    schi2, _, cchi2, dcchi2 = xsincos(chi2)
    s0s = -alg
    b = abs(bz)
    sg = math.copysign(1, bz)
    phix, phiy, phiz = -schi1 * sg, cchi1 * (-schi2) * sg, cchi1 * cchi2 * sg
    dphizsq = phix**2 + phiy**2
    bzs = phiz * b
    pr = 1 + g
    px = px + bz * y / pr * .5
    py = py - bz * x / pr * .5
    ds1 = schi1 * x + dcchi1 * s0s - cchi1 * dz
    x = cchi1 * x - schi1 * (s0s - dz)
    ds2 = -(schi2 * y + cchi2 * ds1 + dcchi2 * s0s)
    y = cchi2 * y - schi2 * (ds1 + s0s)
    pz0 = 1 + pz_m1(px, py)
    pz1 = schi1 * px + cchi1 * pz0
    px = cchi1 * px - schi1 * pz0
    py = cchi2 * py - schi2 * pz1
    bzp = bzs / pr
    alb = pr / b
    dpz0 = pz_m1(px, py)
    pz0 = 1 + dpz0
    dpl = px * phix + py * phiy + dpz0 * phiz
    pl = phiz + dpl
    plx, ply, plz = pl * phix, pl * phiy, pl * phiz
    ptx, pty = px - plx, py - ply
    ptz = dpz0 - dpl * phiz + dphizsq
    pbx, pby, pbz = pty * phiz - ptz * phiy, ptz * phix - ptx * phiz, ptx * phiy - pty * phix
    ps2 = ds2 / alb
    phi = math.asin(ps2 / pz0)
    for _ in range(10):
        sinphi, xsinphi, cosphi, dcos = xsincos(phi)
        s = plz * xsinphi + pz0 * sinphi - pbz * dcos
        dphi = (ps2 - s) / (pz0 + ptz * dcos + pbz * sinphi)
        phi0 = phi
        phi += dphi
        if phi0 == phi or abs(dphi / phi) < 3e-16:
            break
    sinphi, xsinphi, cosphi, dcos = xsincos(phi)
    dcu = update_sign * dcos
    x = x + (plx * phi + ptx * sinphi + pbx * dcu) * alb
    y = y + (ply * phi + pty * sinphi + pby * dcu) * alb
    z = z - phi * alb
    px = px - ptx * dcu + pbx * sinphi - bzp * y * .5
    py = py - pty * dcu + pby * sinphi + bzp * x * .5
    return np.array([x, px, y, py, z, g]), bzs


def exit_(c, chi1, chi2, dz, bz, bzs, alg):
    x, px, y, py, z, g = c
    schi1, _, cchi1, dcchi1 = xsincos(chi1)
    schi2, _, cchi2, dcchi2 = xsincos(chi2)
    s0s = alg
    pr = 1 + g
    px = px + bzs / pr * y * .5
    py = py - bzs / pr * x * .5
    pz0 = 1 + pz_m1(px, py)
    pz1 = -schi2 * py + cchi2 * pz0
    pyi = cchi2 * py + schi2 * pz0
    pz2 = -schi1 * px + cchi1 * pz1
    pxi = cchi1 * px + schi1 * pz1
    y1 = cchi2 * y - schi2 * s0s
    ds1 = -schi2 * y - dcchi2 * s0s
    ds2 = -schi1 * x + cchi1 * ds1 - dcchi1 * s0s + dz
    x1 = cchi1 * x + schi1 * (ds1 - s0s)
    bzp = bz / pr
    phi = -bzp * ds2 / pz2
    a24, a12, a22, a14 = xsincos(phi)
    xn = x1 + (a24 * pxi - a14 * pyi) / bzp
    yn = y1 + (a14 * pxi + a24 * pyi) / bzp
    return np.array([xn, a22 * pxi + a24 * pyi - bzp * yn * .5, yn, -a24 * pxi + a22 * pyi + bzp * xn * .5, z + ds2 / pz2, g])


def composite_dev(sign):
    bz = 13.0 / BRHO
    out = []
    for c in PTS:
        c1, bzs = ent(np.array(c), CHI1, CHI2, DZ, bz, 0.0, sign)
        out.append(np.abs(exit_(c1, CHI1, CHI2, DZ, bz, bzs, 0.0) - c).max())
    return np.array(out)


def test_numpy_composite_is_identity_only_with_the_sign_flipped():
    """numpy の写し: 長さ 0 の傾いた要素 (入口 + 出口) の恒等からのずれ。Fortran のまま (as) 2e-7〜6e-6、更新の dcosphi を反転 (fix) 1e-15 以内"""
    a = composite_dev(+1)
    f = composite_dev(-1)
    assert a == pytest.approx([2.15e-7, 3.87e-7, 1.45e-6, 6.48e-6], rel=0.01)
    assert f.max() < 1e-15


def test_gs_tilted_zero_length_mult_in_sol_is_not_identity():
    """gs: SOL 内の傾いた長さ 0 の MULT は恒等でない (2.2e-7、3.9e-7、1.4e-6、6.5e-6)。傾きなしは戻る。numpy の as のずれと 3% 以内で同じ大きさ"""
    plain = parse('v250_tsolrot_tilt__plain.gs.txt')[:6].T
    tilt = parse('v250_tsolrot_tilt__tilt.gs.txt')[:6].T
    assert np.abs(plain - PTS).max() < 1e-13
    dev = np.array([np.abs(tilt[i] - PTS[i]).max() for i in range(len(PTS))])
    assert dev == pytest.approx([2.19e-7, 3.88e-7, 1.43e-6, 6.53e-6], rel=0.01)
    assert dev == pytest.approx(composite_dev(+1), rel=0.03)
    assert np.all(dev > 1e3 * composite_dev(-1))


# ---------------------------------------------------------------- (F) MULT の CHI1・CHI2 の符号
def ode(chi, plane, K=0.3, L=1.0):
    """ヘルプの R₁・R₂ の要素の軸 (R の 3 行目) のまわりの厚い四極 (x 収束・y 発散) の出口の (座標, 運動量)。軸の傾き: CHI1 は x に −tan χ1、CHI2 は y に −tan χ2"""
    if plane == 'x':
        f = lambda s, u: [u[1], -K * (u[0] + math.tan(chi) * (s - L / 2))]
    else:
        f = lambda s, u: [u[1], +K * (u[0] + math.tan(chi) * (s - L / 2))]
    return solve_ivp(f, (0, L), [0, 0], rtol=1e-12, atol=1e-15).y[:, -1]


@pytest.mark.parametrize('name,chi,plane,col', [('chi1p', 0.05, 'x', (0, 1)), ('chi1m', -0.05, 'x', (0, 1)),
                                                ('chi2p', 0.05, 'y', (2, 3)), ('chi2m', -0.05, 'y', (2, 3))])
def test_mult_chi_sign_matches_help_rotation_matrices(name, chi, plane, col):
    """gs の MULT の CHI1・CHI2 (±0.05) の出口 (軸上の粒子) が、ヘルプの R₁・R₂ から決まる軸の傾きの ODE と 0.3% 以内 (符号まで一致)"""
    G = parse('v250_mult_chi_sign__%s.gs.txt' % name)
    got = np.array([G[col[0], 0], G[col[1], 0]])
    ref = ode(chi, plane)
    assert got == pytest.approx(ref, rel=3e-3)
    assert np.sign(got[0]) == np.sign(ref[0]) != 0


def test_mult_chi_flips_sign_and_none_is_zero():
    a, b = parse('v250_mult_chi_sign__chi1p.gs.txt'), parse('v250_mult_chi_sign__chi1m.gs.txt')
    assert np.abs(a[:4, 0] + b[:4, 0]).max() < 1e-15
    assert np.abs(parse('v250_mult_chi_sign__none.gs.txt')[:4, 0]).max() == 0.0
