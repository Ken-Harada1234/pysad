# -*- coding: utf-8 -*-
"""v2.5.0: CAVI の端の収束 (tcavfrin) の 1/β (docs/30 B330)。

K64 の端の収束は Δp_x/p0 = x·v/(2 L p0)·(sin(ωt − dphis) + s0) で、β = p/h の因子が無い。
軸の近くの Maxwell の場 (E_z = E0 g(z) f(t)、E_r = −(r/2)∂_z E_z、B_θ = (r/2c²)∂_t E_z、E_z の r² の項 −(r²/4)(∂_z² − ∂_t²)E_z) の中の
ローレンツ力の数値積分 (DOP853) では、端を通る間の Δp_x = x·v/(2L β)·(…) で 1/β がある。エネルギー変化 (x² の項) にも同じ 1/β がある。
v2 は追跡 (python 版・numba 版)・CALC・Emittance の写像に 1/β (粒子の β = p/h) と、H_f ∝ 1/β(E) の整合 (時間のずれ −∂H_f/∂E) を入れた。k64 は K64 のまま (gs と一致)。

独立な基準 (pysad のコードを使わない): 上の場の中の電子の運動方程式の積分。端の幅 d = 1e-5 m、L = 2e-4 m、E0 = 1e5 V/m、500 MHz。
gs の保存した出力 (scripts/v250_cavi_edge_beta__*.gs.txt。L = 1e-4 m、VOLT = 10 V、FREQ = 500 MHz、PHI = 0.3 の CAVI の FRINGE = 1 と DISFRIN = 1 の
TrackParticles 1 要素。x = 1e-5) と k64 を比べる。
限界: 端 1 つ (入口と出口を別々) の 1 次の収束とエネルギー変化の x² の項だけ。長い空洞の本体との組合せ (位相の進みの扱い) は見ない。
"""
import io
import math
import os
import re

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, track6, track6n_ext
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')
MC2 = 0.51099895e6
CLIGHT = 299792458.0
FREQ = 500e6
W = 2 * math.pi * FREQ / CLIGHT


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def mk(el, momentum):
    src = 'MOMENTUM = %s;\n%s\nMARK ST = ();\nLINE ONE = (ST CV);\nFFS;\nUSE ONE;\n' % (momentum, el)
    ffs, _ = run_sad_string(src, out=io.StringIO())
    return ffs


def trk(ffs, c, ext=False):
    a = np.zeros((1, 6))
    a[0, :] = c
    fn = track6n_ext.track_particles_ext if ext else track6.track_particles
    return np.asarray(fn(ffs.lat.elements, a, 1, ffs.momentum(), 0.0, True, False, False))[0]


def lorentz_edge(p0, delta, x0, v, L, phi, edge='in', d=1e-5, with_r2=False):
    """端 1 つを通る電子の運動 (c = 1、mc² = 1、長さ m)。入口: g = ½(1 + tanh(z/d))、出口: g = ½(1 − tanh(z/d))。
    qE_z(r, z, t) = −(v/L) g(z) sin(ωt − dphis) − (r²/4)(∂_z² − ∂_t²)(同)、qE_x = −(x/2)∂_z(qE_z)、qB_y = (x/2)∂_t(qE_z)。
    追跡の sin(ωt − dphis) の dphis は、phis = 0、電子 (charge = −1) で −PHI。粒子は z = 0 (端の中心) を t = 0 に通る。
    出口は一様な場の中 (z = −L) から出発して z = 40 d まで。戻り値: (px/p0、エネルギー変化 ∫ q E·v dt)"""
    dphis = -phi
    sgn = 1.0 if edge == 'in' else -1.0
    pz0 = p0 * (1.0 + delta)
    beta = pz0 / math.sqrt(1 + pz0 * pz0)
    A = -(v / L)

    def rhs(t, y):
        x, z, px, pz, E = y
        gam = math.sqrt(1 + px * px + pz * pz)
        s = math.tanh(z / d)
        g = 0.5 * (1 + sgn * s)
        gp = 0.5 * sgn * (1 - s * s) / d
        gpp = -sgn * s * (1 - s * s) / d**2
        psi = W * t - dphis
        sn, cs = math.sin(psi), math.cos(psi)
        Qz = A * g * sn
        if with_r2:
            Qz -= (x * x / 4) * A * (gpp * sn + g * W * W * sn)
        Qx = -(x / 2) * A * gp * sn
        By = (x / 2) * A * g * W * cs
        vx, vz = px / gam, pz / gam
        return [vx, vz, Qx - vz * By, Qz, Qz * vz + Qx * vx]

    if edge == 'in':
        zs, zend = -40 * d, L
    else:
        zs, zend = -L, 40 * d
    ev = lambda t, y: y[1] - zend
    ev.terminal = True
    ts = zs / beta
    r = solve_ivp(rhs, (ts, ts + 10), [x0, zs, 0.0, pz0, 0.0], method='DOP853', rtol=1e-13, atol=1e-19, events=ev)
    y = r.y[:, -1]
    return y[2] / p0, y[4]


MOMS = [('2MEV', 2.0e6), ('600KEV', 0.6e6)]
L = 2e-4
E0 = 1e5
V = E0 * L
VN = V / MC2


def cav(name, phi, fr):
    return mk('CAVI CV=(L=%g VOLT=%g FREQ=500E6 PHI=%g %s);' % (L, V, phi, fr), name)


# ---------------------------------------------------------------- (i) v2: 1/β のある収束が数値積分と合う
@pytest.mark.parametrize('name,mom', MOMS)
@pytest.mark.parametrize('phi', [0.3, 1.2, -0.7])
def test_v2_entrance_kick_matches_lorentz_integration(name, mom, phi):
    """v2: 入口の端 (FRINGE=1) の px/x が数値積分と 1% 以内 (残りは本体の Panofsky–Wenzel と薄い空洞の近似 O(ωL))。k64 は β 倍 (0.6 MeV/c で 0.76 倍)"""
    p0 = mom / MC2
    beta = p0 / math.sqrt(1 + p0 * p0)
    x0 = 1e-5
    ffs = cav(name, phi, 'FRINGE=1')
    px_num, _ = lorentz_edge(p0, 0.0, x0, VN, L, phi)
    c0 = np.array([x0, 0, 0, 0, 0, 0.0])
    compat.set_mode('v2')
    k2 = trk(ffs, c0)[1]
    compat.set_mode('k64')
    k6 = trk(ffs, c0)[1]
    assert abs(k2 / px_num - 1) < 0.01, (k2, px_num)
    assert abs(k6 / px_num - beta) < 0.012, (k6 / px_num, beta)


@pytest.mark.parametrize('name,mom', MOMS)
@pytest.mark.parametrize('phi', [0.3, 1.2])
def test_v2_exit_kick_matches_lorentz_integration(name, mom, phi):
    """v2: 出口の端 (FRINGE=2 と DISFRIN=1 の差) も 1/β。一様な場の中から出る粒子の px/x (数値積分は本体の寄与 O(ωL) = 1% を含む)"""
    p0 = mom / MC2
    beta = p0 / math.sqrt(1 + p0 * p0)
    x0 = 1e-5
    px_num, _ = lorentz_edge(p0, 0.0, x0, VN, L, phi, edge='out')
    ffs_o = cav(name, phi, 'FRINGE=2')
    ffs_n = cav(name, phi, 'DISFRIN=1')
    res = {}
    c0 = np.array([x0, 0, 0, 0, 0, 0.0])
    for m in ('v2', 'k64'):
        compat.set_mode(m)
        res[m] = trk(ffs_o, c0)[1] - trk(ffs_n, c0)[1]
    assert abs(res['v2'] / px_num - 1) < 0.02, (res['v2'], px_num)
    assert abs(res['k64'] / res['v2'] - beta) < 1e-3


@pytest.mark.parametrize('name,mom', MOMS)
def test_v2_energy_change_x2_matches_lorentz_integration(name, mom):
    """端のエネルギー変化の x² の項 (∂H_f/∂t) にも 1/β。E_z の r² の項まで入れた場の積分と 1% 以内、k64 は β 倍"""
    p0 = mom / MC2
    beta = p0 / math.sqrt(1 + p0 * p0)
    phi = 0.3
    x0 = 3e-4
    ffs = cav(name, phi, 'FRINGE=1')
    _, Ea = lorentz_edge(p0, 0.0, x0, VN, L, phi, with_r2=True)
    _, Eb = lorentz_edge(p0, 0.0, 0.0, VN, L, phi, with_r2=True)
    dE = Ea - Eb
    out = {}
    for m in ('v2', 'k64'):
        compat.set_mode(m)
        a = trk(ffs, np.array([x0, 0, 0, 0, 0, 0.0]))
        b = trk(ffs, np.array([0, 0, 0, 0, 0, 0.0]))
        out[m] = math.sqrt(1 + (p0 * (1 + a[5]))**2) - math.sqrt(1 + (p0 * (1 + b[5]))**2)
    assert abs(out['v2'] / dE - 1) < 0.01, (out['v2'], dE)
    assert abs(out['k64'] / dE - beta) < 0.012


@pytest.mark.parametrize('delta', [-0.05, 0.05, 0.15])
def test_v2_uses_the_particle_beta_not_the_design_beta(delta):
    """δ ≠ 0 の粒子の端の収束は、その粒子の β = p0(1+δ)/h の 1/β。数値積分 (初期の p_z = p0(1+δ)) と比べる (px は p0 規格化の正準運動量)。0.6 MeV/c、PHI = 0.3"""
    p0 = 0.6e6 / MC2
    phi = 0.3
    x0 = 1e-5
    ffs = cav('600KEV', phi, 'FRINGE=1')
    c = np.array([x0, 0, 0, 0, 0, delta])
    px_num, _ = lorentz_edge(p0, delta, x0, VN, L, phi)
    compat.set_mode('v2')
    o = trk(ffs, c)
    px2 = o[1]                             # 追跡の px は p0 規格化の正準運動量
    assert abs(px2 / px_num - 1) < 0.01, (px2, px_num)
    beta_p = p0 * (1 + delta) / math.sqrt(1 + (p0 * (1 + delta))**2)
    compat.set_mode('k64')
    ok = trk(ffs, c)
    pxk = ok[1]
    assert abs(pxk / px_num - beta_p) < 0.012, (pxk / px_num, beta_p)


# ---------------------------------------------------------------- (ii) k64: gs と同じ
def _parse(txt):
    out = {}
    for m in re.finditer(r'(ONE|TWO) PX \{\s*(\S+),\s*\S+\} DP \{\s*(\S+),\s*\S+\}', txt):
        out[m.group(1)] = (float(m.group(2).replace('D', 'E')), float(m.group(3).replace('D', 'E')))
    return out


@pytest.mark.parametrize('tag,name', [('3gev', '3GEV'), ('2mev', '2MEV'), ('600kev', '600KEV')])
def test_k64_matches_saved_gs(tag, name):
    """k64: gs の TrackParticles (L = 1e-4 m、VOLT = 10 V、FRINGE=1 / DISFRIN=1、x = 1e-5) と px・δ が相対 1e-9 以内で一致 (表示は 12〜13 桁)"""
    g = _parse(open(os.path.join(SCRIPTS, 'v250_cavi_edge_beta__%s.gs.txt' % tag), encoding='utf-8').read())
    assert set(g) == {'ONE', 'TWO'}
    compat.set_mode('k64')
    for key, fr in (('ONE', 'FRINGE=1'), ('TWO', 'DISFRIN=1')):
        ffs = mk('CAVI CV=(L=0.0001 VOLT=10 FREQ=500E6 PHI=0.3 %s);' % fr, name)
        o = trk(ffs, np.array([1e-5, 0, 0, 0, 0, 0.0]))
        assert o[1] == pytest.approx(g[key][0], rel=1e-9, abs=1e-19), (key, o[1], g[key][0])
        assert o[5] == pytest.approx(g[key][1], rel=1e-9, abs=1e-19)


def test_saved_gs_shows_no_beta_factor():
    """保存した gs の端の収束 (ONE − TWO) が x v/(2 L p0) sin(PHI) に 1e-6 で一致する (600 keV/c で 1/β = 1.31 倍にはなっていない)"""
    g = _parse(open(os.path.join(SCRIPTS, 'v250_cavi_edge_beta__600kev.gs.txt'), encoding='utf-8').read())
    kick = g['ONE'][0] - g['TWO'][0]
    p0 = 0.6e6 / MC2
    h0 = math.sqrt(1 + p0 * p0)
    formula = 1e-5 * (10.0 / MC2) / (2 * 1e-4 * p0) * math.sin(0.3)
    assert kick / formula == pytest.approx(1.0, abs=1e-6)
    assert abs(kick / formula - h0 / p0) > 0.3


# ---------------------------------------------------------------- (iii) 他の挙動
def test_v2_k64_numba_equals_python():
    """python 版と numba 版の追跡が、k64・v2 とも同じ (入口のみ・出口のみ・両方、δ ≠ 0)"""
    c = np.array([3e-3, 2e-4, -2e-3, 1e-4, -0.02, 0.03])
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        for fr in ('FRINGE=1', 'FRINGE=2', 'FRINGE=3'):
            ffs = mk('CAVI CV=(L=0.05 VOLT=2E4 FREQ=500E6 PHI=0.3 %s);' % fr, '600KEV')
            assert np.abs(trk(ffs, c, True) - trk(ffs, c)).max() < 1e-14, (m, fr)


def test_v2_ultrarelativistic_and_no_fringe_unchanged():
    """3 GeV/c (1/β − 1 = 1.5e-8) では v2 と k64 の差が 1e-7 (相対) 以下。DISFRIN=1 と L=0 (端なし) は v2 でも k64 と同じ"""
    c = np.array([1e-3, 2e-4, -5e-4, 1e-4, -0.02, 0.01])
    ffs = mk('CAVI CV=(L=0.1 VOLT=1E5 FREQ=500E6 PHI=0.3 FRINGE=3);', '3GEV')
    compat.set_mode('k64')
    a = trk(ffs, c)
    compat.set_mode('v2')
    b = trk(ffs, c)
    assert np.abs(a - b).max() < 1e-7 * np.abs(a).max() + 1e-13
    for el in ('CAVI CV=(L=0.1 VOLT=1E5 FREQ=500E6 PHI=0.3 DISFRIN=1);', 'CAVI CV=(L=0 VOLT=1E5 FREQ=500E6 PHI=0.3);'):
        ffs = mk(el, '600KEV')
        compat.set_mode('k64')
        a = trk(ffs, c)
        compat.set_mode('v2')
        assert np.array_equal(a, trk(ffs, c)), el


def _symp_err(ffs, c, h=1e-6):
    S = np.zeros((6, 6))
    for i in (0, 2, 4):
        S[i, i + 1] = 1
        S[i + 1, i] = -1
    J = np.zeros((6, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h
        J[:, j] = (trk(ffs, c + d) - trk(ffs, c - d)) / (2 * h)
    return np.abs(J.T @ S @ J - S).max()


def test_v2_edge_keeps_symplectic_structure():
    """端の 1/β を入れても、時間のずれ (−∂H_f/∂E) があるので、ヤコビアンのシンプレクティック性が k64 と同程度 (時間のずれが無いと 7e-4 で 1 桁悪い)。
    0.6 MeV/c、L = 1 mm、VOLT = 1 kV、x = 3 mm、δ = 0.03 (k64 自身が 6e-5 ずれる値。ずれの主因は端以外)"""
    c = np.array([3e-3, 2e-4, -2e-3, 1e-4, -0.02, 0.03])
    ffs = mk('CAVI CV=(L=1E-3 VOLT=1E3 FREQ=500E6 PHI=0.3 FRINGE=1);', '600KEV')
    compat.set_mode('k64')
    e6 = _symp_err(ffs, c)
    compat.set_mode('v2')
    e2 = _symp_err(ffs, c)
    assert e6 < 2e-4 and e2 < 2e-4, (e6, e2)
    assert e2 < 3 * e6


def test_v2_calc_edge_map_matches_tracking():
    """CALC・Emittance の端の写像 (emit6d._tcavfrie_map。dphis = −PHI、s0 = 0) の変化量 (px・δ・z) が、v2 で追跡の端 (FRINGE=1 と DISFRIN=1 の差) と 2e-3 (z は 6e-2) 以内。
    x = 3 mm、y = −2 mm、z = −2 mm、δ = 0 と 0.05、0.6 MeV/c、L = 1 mm、VOLT = 1 kV (差は本体との O(V²) の結合)。k64 は δ = 0 だけ (CALC の px は /p1r を持たない K64 の形)"""
    from pysad import emit6d
    p0 = 0.6e6 / MC2
    v = 1e3 / MC2
    ffs1 = mk('CAVI CV=(L=1E-3 VOLT=1E3 FREQ=500E6 PHI=0.3 FRINGE=1);', '600KEV')
    ffs0 = mk('CAVI CV=(L=1E-3 VOLT=1E3 FREQ=500E6 PHI=0.3 DISFRIN=1);', '600KEV')
    for m, deltas in (('v2', (0.0, 0.05)), ('k64', (0.0,))):
        compat.set_mode(m)
        f = emit6d._tcavfrie_map(1e-3, v, W, -0.3, 0.0, p0)
        for dl in deltas:
            c = np.array([3e-3, 0.0, -2e-3, 0.0, -0.002, dl])
            d_trk = trk(ffs1, c) - trk(ffs0, c)
            d_map = f(c) - c
            for k in ((1, 3, 4, 5) if m == 'v2' else (1, 3, 5)):
                ref = max(abs(d_trk[k]), 1e-300)
                tol = 6e-2 if k == 4 else 2e-3          # z は FRINGE=1 の本体のドリフトが δ の変化 (−9e-5) で z を 5e-8 ずらす分を含む
                assert abs(d_map[k] - d_trk[k]) < tol * ref + 1e-15, (m, dl, k, d_map[k], d_trk[k])
