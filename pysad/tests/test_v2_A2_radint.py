# -*- coding: utf-8 -*-
"""担当 A2 (3): RADINT・Emittance[]・追跡の整合を、大きなリングで自前の積分と比べる。

自前の計算は `_a2_indep.py` (要素表から 3x3 行列で I1..I5、U0・Jx・Je・τ・εx・σδ・α。SAD・pysad の光学を使わない)。
条件: 大きなリング A (kek-ls_v03_12.lat、3 GeV、周長 603 m、平面)、PF リング (2.5 GeV、187 m)、ASTRID2 (0.58 GeV、45.7 m) と MPW 入り。
結論 (docs/history/v2next_A2.md §3): RADINT = 自前の積分 (1e-6 以下)。Emittance[] は EMITDIV=20 で εx・σδ が 1e-4 台、Jx・Je・τ が U0/E の
オーダー (3 GeV で 1e-4、8 GeV で 1e-3) で合う (RADINT は U0/E の 1 次の理論なので、残りは Emittance[] と追跡の U0/E の 2 次)。
v2 の追加 (A2): ANGLE=0 で K0≠0 の BEND (挿入光源の極) を RADINT の積分に入れた (K64 は無視)。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, radint, track6n

from . import _a2_indep as A

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RINGS = {
    'A': (os.path.join(HERE, 'kek-ls_v03_12.lat'), 'RING'),
    'PF': (os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), 'RING'),
    'AST': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING'),
    'ASTMPW': (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad'), 'RINGMPW'),
}


def _load(name, extra='RFSW; RAD; NOFLUC;'):
    path, ring = RINGS[name]
    if not os.path.exists(path):
        pytest.skip(os.path.basename(path) + ' がない')
    with open(path, encoding='utf-8') as fh:
        src = fh.read()
    ffs, _ = ffscmd.run_sad_string(src + f'\nFFS;\nUSE {ring};\nCELL; CALC; {extra}\nExit;\n', out=io.StringIO())
    return ffs


def _radint_I(ffs):
    res, _m = radint.compute(ffs, ffs.ensure_rows())
    return [float(res['RI'][i]) for i in range(1, 12)]


@pytest.mark.parametrize('name', ['A', 'PF', 'AST'])
def test_radint_equals_own_integrals_v2(v2_mode, name):
    """RADINT の I1..I5 が、自前の積分 (3x3 行列、スライス + Richardson) と 1e-5 以内 (3 本の大きなリング)"""
    ffs = _load(name)
    els = ffs.lat.elements
    I, _C, _nu = A.integrals_rich(els)
    ri = _radint_I(ffs)[:5]
    for a, b in zip(I, ri):
        assert abs(b / a - 1) < 1e-5, (name, I, ri)


def _own(ffs):
    els = ffs.lat.elements
    I, C, _ = A.integrals_rich(els)
    return A.derived(I, C, ffs.momentum()), C


@pytest.mark.parametrize('name,tol', [
    ('A', dict(U0=1e-6, Jx=3e-4, Je=4e-4, tau=3e-4, ex=1e-4, sd=1e-4, alpha=2e-4)),
    ('PF', dict(U0=5e-5, Jx=1.2e-3, Je=1e-3, tau=1.2e-3, ex=5e-4, sd=5e-4, alpha=1e-3)),
])
def test_emittance_matches_integrals_with_emitdiv_v2(v2_mode, name, tol):
    """Emittance[] (EMITDIV=20) の U0・Jx・Je・減衰時間・εx・σδ・α が自前の積分と tol 以内 (残りは U0/E の 2 次: A 1e-4、PF 1.6e-4 の数倍)"""
    ffs = _load(name)
    d, C = _own(ffs)
    old = list(emit6d.EMITDIV)
    emit6d.EMITDIV[:] = [20.0] * 3
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
    finally:
        emit6d.EMITDIV[:] = old
    p = ffs.momentum()
    T0 = C / A.C_LIGHT
    dr = np.abs(e['damping_rate'])
    u0 = e['energy_loss_eV']
    assert abs(u0 / d['U0'] - 1) < tol['U0']
    assert abs(2 * p * dr[0] / u0 / d['Jx'] - 1) < tol['Jx']
    assert abs(2 * p * dr[2] / u0 / d['Je'] - 1) < tol['Je']
    assert abs(T0 / dr[0] / d['taux'] - 1) < tol['tau']
    assert abs(e['emittance_x'] / d['ex'] - 1) < tol['ex']
    assert abs(e['momentum_spread'] / d['sd'] - 1) < tol['sd']
    assert abs(e['momentum_compaction'] / d['alpha'] - 1) < tol['alpha']


@pytest.mark.parametrize('name', ['A', 'PF'])
def test_tracking_damping_matches_integrals_v2(v2_mode, name):
    """RAD つき 6 次元追跡の一周写像 (閉軌道の周りの数値微分) の 3 つのモードの減衰から求めた分配数が、自前の積分と 1e-3 以内
    (差は U0/E の 2 次の効果。A: Jx 1.56385 対 1.56395、PF: 1.00398 対 1.00422)"""
    ffs = _load(name)
    els = ffs.lat.elements
    p = ffs.momentum()
    circ = sum(x.L for x in els)
    e = emit6d.emittance6d(ffs.lat, p, rfsw=True, radcod=True)
    trf0 = e['trf0']

    def f(c):
        return np.asarray(track6n.track_particles_fast(els, np.array([c]), 1, p, trf0, True, True, False, circ=circ))[0]

    def jac(c, h=1e-7):
        J = np.zeros((6, 6))
        for i in range(6):
            d = np.zeros(6)
            d[i] = h
            J[:, i] = (f(c + d) - f(c - d)) / (2 * h)
        return J
    c0 = np.array(e['orbit_entrance'], float)
    for _ in range(4):
        c0 = c0 - np.linalg.solve(jac(c0) - np.eye(6), f(c0) - c0)
    ev = np.linalg.eigvals(jac(c0))
    fold = lambda x: min(x % 1.0, 1 - x % 1.0)
    mags = []
    for t in e['tunes']:
        z = min((z for z in ev if z.imag >= 0), key=lambda z: abs(fold(math.atan2(z.imag, z.real) / 2 / math.pi) - fold(t)))
        mags.append(abs(z))
    d, _C = _own(ffs)
    J = [2 * p * (1 - m) / d['U0'] for m in mags]
    assert abs(J[0] - d['Jx']) < 1e-3 and abs(J[1] - 1.0) < 1e-3 and abs(J[2] - d['Je']) < 1.5e-3, (J, d['Jx'], d['Je'])


def test_radint_includes_k0_poles_of_insertion_device_v2(v2_mode):
    """ASTRID2 + MPW: ANGLE=0 で K0≠0 の極 (磁場だけ、基準軌道は曲がらない) が RADINT の I2・I3 に入る。
    I2 = Σ (ANGLE+K0)²/L が独立な和と 1e-6、U0 が Emittance[] (追跡と同じ放射) と 2e-3 以内。K64 は極を無視し I2 = 3.8704 (U0 が 21% 小さい)"""
    ffs = _load('ASTMPW')
    els = ffs.lat.elements
    i2 = sum(((e.ANGLE + e.K0) ** 2 / e.L) for e in els if type(e).__name__ == 'Bend' and e.L > 0)
    i2_wig = sum((e.K0 ** 2 / e.L) for e in els if type(e).__name__ == 'Bend' and e.L > 0 and e.ANGLE == 0)
    ri = _radint_I(ffs)
    assert abs(ri[1] / i2 - 1) < 1e-6
    assert 1.0 < i2_wig < 1.05
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True)
    u0 = 2e9 * radint._coeff() * ri[1] * (ffs.momentum() / 1e9) ** 4
    assert abs(u0 / e['energy_loss_eV'] - 1) < 2e-3, (u0, e['energy_loss_eV'])
    assert abs(ri[2] / A_i3_plain_plus_wiggler(els) - 1) < 1e-6


def A_i3_plain_plus_wiggler(els):
    return sum((abs(e.ANGLE + e.K0) ** 3 / e.L ** 2) for e in els if type(e).__name__ == 'Bend' and e.L > 0)


def test_radint_ignores_k0_poles_k64(k64_mode):
    ffs = _load('ASTMPW')
    ri = _radint_I(ffs)
    assert abs(ri[1] - 3.870436) < 1e-4      # MPW なしの ASTRID2 と同じ


# ---------------------------------------------------------------------------------------------
# 結合のあるリング (傾いた四極): RADINT の I4V
# ---------------------------------------------------------------------------------------------

def _coupled(ks):
    return f'''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
QUAD  SQ = (L=0.2 K1={ks} ROTATE=0.7853981633974483);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=1.8E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE CELLS = (QF D2 B D2 SQ QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 4*(CELLS 5*CELLQ));
FFS;
USE RINGQ;
CELL; CALC; RFSW; RAD; NOFLUC;
Exit;
'''


@pytest.mark.parametrize('ks', [0.05, 0.15])
def test_radint_i4v_coupled_ring_v2(v2_mode, ks):
    """傾いた四極 (4 台、積分強度 ks) のあるリング: (a) RADINT の I4U + I4V = 独立な物理的 I4 (結合のある 5x5 行列の η_x から) と 1e-5、
    (b) RADINT の DU・DV・Je = 4 - DU - DV が、自前の 6x6 行列と放射の 1 次の摂動論 (固有ベクトル) による各モードの減衰分配と 3e-3 以内
    (自前の摂動論の基準の差 4e-4 を含む)。ks=0.15 は νy と νs が近い (0.047 と 0.028) 条件。条件: 3 GeV、RF 1.8 MV (νs=0.014)"""
    ffs, _ = ffscmd.run_sad_string(_coupled(ks), out=io.StringIO())
    ri = _radint_I(ffs)
    els = ffs.lat.elements
    i4, _i2 = A.coupled_i4(els)
    assert abs((ri[3] + ri[6]) / i4 - 1) < 1e-5, (ri[3], ri[6], i4)
    assert abs(ri[6]) > 1e-3                     # 結合で I4V が 0 でない
    modes, _u0 = A.coupled_partition(els, ffs.momentum())
    du, dv = 1 - ri[3] / ri[1], 1 - ri[6] / ri[1]
    dE = 4 - du - dv
    own_sorted = sorted(modes, key=lambda t: -t[0])          # ν の大きい順: x 様 (0.24→0.14)、y 様、z 様 (νs)
    jx, jz = own_sorted[0][1], own_sorted[2][1]
    jy = own_sorted[1][1]
    assert abs(du - jx) < 3e-3 and abs(dv - jy) < 3e-3 and abs(dE - jz) < 3e-3, ((du, dv, dE), (jx, jy, jz))
    assert abs(sum(m[1] for m in modes) - 4.0) < 1e-6        # 和則


def test_radint_i4v_is_unchanged_in_k64():
    """K64 でも I4V は同じ (v2 で変えていない)"""
    out = {}
    for mode in ('k64', 'v2'):
        old = compat._mode
        compat.set_mode(mode)
        try:
            ffs, _ = ffscmd.run_sad_string(_coupled(0.05), out=io.StringIO())
            out[mode] = _radint_I(ffs)
        finally:
            compat.set_mode(old)
    assert out['k64'][6] == pytest.approx(out['v2'][6], abs=1e-12)
