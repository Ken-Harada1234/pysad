# -*- coding: utf-8 -*-
"""担当 B1 (v2next 項目 6): SOL 領域内の POL (スピン)・加速 + RAD・bz ≥ 3 1/m の確認と、足りなかった DRIFT の刻み。

根拠 (gs は使わない。docs/46 §3 の「領域内の未確認」):
  - bz ≥ 3 1/m: 領域内の要素の写像がローレンツ力の数値積分と一致する範囲を bz = 20 1/m まで広げた (test_v2_C_sol の bz = 3 までの拡張)。放射損失は
    (2/3) r_e γ³ (K0/L)² L と bz = 20 で 3e-5。
  - 領域内の DRIFT のスピンと放射: 放射の 1 刻み (弦から損失を出す) が bz·al ≳ 1 で粗く、スピンの歳差が Thomas-BMT の積分と 1e-3〜1e-2 ずれた
    (領域内の BEND は bz·al/ndiv ≤ 0.02 に刻んでいたが DRIFT は 1 刻みだった)。v2 は DRIFT も刻んで 1e-6 以下にした (k64 はそのまま)。
  - POL の Emittance[]: 平面リングの SOL 領域 (DRIFT・MULT・CAVI) のスピンチューンが cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1+a) ∮Bz dl/(Bρ) と 1e-9 で一致する
    (aγ は同じリングのソレノイドなしの Emittance[] のスピンチューン)。bz = 4.5・9 1/m まで。
  - 加速 + RAD: 領域内の RF つき MULT (TRPT)。放射損失は (2/3) r_e mc² (bz γ0² θ0)² L で、加速に依らない (γ⁴κ² が一定: 横運動量は
    ソレノイドの縦磁場では変わらず曲率 κ ∝ p⊥/p²)。RAD ON−OFF の δ の差がこの式と 1e-3。
  - POL + TRPT + 加速 (track6 の経路) は基準運動量が変わらず黙って δ = 0.27 になっていた。明示のエラーにした。
"""
import io
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from pysad import compat, emit6d as _e, track6n as tn, track6n_ext as x

GS = 0.001159652181280002
ME = 0.51099895e6
CL = 299792458.0
RE = 2.8179403e-15


def _run(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    ffs, _ = run_sad_string(src, out=buf)
    return ffs, buf.getvalue()


def _elem(defs, mom='1GEV'):
    import pysad.ffscmd as F
    holder = {}
    orig = F.FFS.use

    def use(self, *a, **k):
        r = orig(self, *a, **k)
        holder['f'] = self
        return r
    F.FFS.use = use
    try:
        F.run_sad_string('MOMENTUM = %s;\n' % mom + defs + '\nMARK ST=();\nLINE BL=(ST E1);\nFFS;\nUSE BL;\nExit;\n', out=io.StringIO())
    finally:
        F.FFS.use = orig
    return [e for e in holder['f'].els if e.name == 'E1'][0]


class _Col:
    def __init__(self):
        self.items = []

    def _add(self, k, p, name, cnt, sm):
        self.items.append((k, p))


def _kin(v, bz):
    return np.array([v[0], v[1] + bz * v[2] / 2, v[2], v[3] - bz * v[0] / 2])


def _ode(c, L, bz, K):
    x0, px, y0, py = c[:4]
    pn = 1.0 + c[5]

    def f(s, u):
        X, PX, Y, PY = u
        pz = math.sqrt(pn * pn - PX * PX - PY * PY)
        w = complex(X, Y)
        b = 0j
        t = 1 + 0j
        fct = 1.0
        for n in range(len(K)):
            if n > 0:
                fct *= n
            b += K[n] / L * t / fct
            t *= w
        By, Bx = b.real, b.imag
        return [PX / pz, (PY * bz - pz * By) / pz, PY / pz, (pz * Bx - PX * bz) / pz]
    return solve_ivp(f, [0, L], [x0, px, y0, py], rtol=1e-13, atol=1e-16).y[:, -1]


C0 = np.array([1e-3, 2e-3, -1.5e-3, 1e-3, 0.0, 0.01])


def _trk(el, c, bz, p=1.0e9, rad=False):
    ph = _e._Phys(p, 100.0)
    if type(el).__name__ == 'Drift':
        k, pr = x.K_SOL_DRIFT, (float(el.L), bz, True)                # 領域内の DRIFT は compile ではなく呼び出し側が K_SOL_DRIFT にする
    else:
        col = _Col()
        x._compile_sol_v2(col, el, type(el).__name__, bz, ph, 2)
        k, pr = col.items[0]
    S = tn.make_state(np.array([c], float), ph)
    x.step_kind(k, pr, S, 1, ph, 0.0, False, rad, False, tn.make_rs(1))
    return tn.state_out(S, 1)[0], ph


_CASES = [
    ('DRIFT E1=(L=0.5);', [0.0]),
    ('BEND E1=(L=0.5 ANGLE=0.02);', [0.02]),
    ('MULT E1=(L=0.5 K0=0.02 DISFRIN=1);', [0.02]),
    ('MULT E1=(L=0.5 K1=0.3 DISFRIN=1);', [0, 0.3]),
    ('MULT E1=(L=0.5 K0=0.01 K1=0.3 K2=30 DISFRIN=1 EPS=0.01);', [0.01, 0.3, 30]),
]


@pytest.mark.parametrize('bz', [6.0, 10.0, 20.0])
@pytest.mark.parametrize('case', range(len(_CASES)))
def test_sol_region_element_matches_lorentz_force_at_large_bz(v2_mode, bz, case):
    """bz = 6・10・20 1/m (BZ = 60 T @ 3 GeV 相当の大きな値でカーネルの精度だけを見る): 領域内の DRIFT・BEND・MULT の写像 (力学的運動量)
    がローレンツ力の積分と 1e-8 以内。bz = 3 までの試験 (test_v2_C_sol) を広げたもの。残る差は多極の分割 (EPS)"""
    defs, K = _CASES[case]
    o, _ = _trk(_elem(defs), C0, bz)
    r = _ode(list(_kin(C0, bz)) + [0.0, C0[5]], 0.5, bz, [complex(k) for k in K])
    k1 = _kin(o, bz)
    assert max(abs(r[j] - k1[j]) for j in range(4)) < 1.2e-8


@pytest.mark.parametrize('bz', [6.0, 20.0])
def test_sol_region_radiation_loss_at_large_bz(v2_mode, bz):
    """放射損失 Δδ = −(2/3) r_e γ³ (K0/L)² L (BEND と K0 だけの MULT) が bz = 6・20 1/m でも 4e-5 以内"""
    for defs in ('BEND E1=(L=0.5 K0=0.02);', 'MULT E1=(L=0.5 K0=0.02 DISFRIN=1);'):
        o, ph = _trk(_elem(defs), np.zeros(6), bz, rad=True)
        exp = -(2.0 / 3.0) * RE * ph.h0**3 * (0.02 / 0.5)**2 * 0.5
        assert abs(o[5] / exp - 1.0) < 4e-5


@pytest.mark.parametrize('bz,th0', [(3.0, 1e-2), (9.0, 1e-2)])
def test_sol_region_drift_radiation_loss_helix(v2_mode, bz, th0):
    """領域内の DRIFT の放射 (v2: bz·al/ndiv ≤ 0.02 に刻む): 横運動量 p⊥ = p0 θ0 (正準 px = θ0 − ...) のらせんの放射損失 Δδ = −(2/3) r_e γ0³ (bz θ0)² L が
    bz = 3・9 1/m (bz L = 1.5・4.5) で 1e-4 以内 (K64 の 1 刻みは (bz L)²/12 だけ小さい: bz = 9 で 1.7)。3 GeV"""
    P = 3.0e9
    el = _elem('DRIFT E1=(L=0.5);', '3GEV')
    c = np.array([0.0, th0 - 0.0, 0.0, 0.0, 0.0, 0.0])
    # 力学的運動量 th0 を x 方向に持つ粒子 (正準 px = 力学的 px、y = 0)
    o, ph = _trk(el, c, bz, p=P, rad=True)
    exp = -(2.0 / 3.0) * RE * ph.h0**3 * (bz * th0)**2 * 0.5
    assert abs(o[5] / exp - 1.0) < 1e-4, (o[5], exp)


def test_sol_region_drift_radiation_k64_is_a_single_step(k64_mode):
    """k64: 領域内の DRIFT は 1 刻み (K64 の tdrift_solenoid)。bz L = 4.5 で放射損失が大きく外れる (記録: v2 の刻みの根拠)"""
    P = 3.0e9
    el = _elem('DRIFT E1=(L=0.5);', '3GEV')
    o, ph = _trk(el, np.array([0.0, 1e-2, 0.0, 0.0, 0.0, 0.0]), 9.0, p=P, rad=True)
    exp = -(2.0 / 3.0) * RE * ph.h0**3 * (9.0 * 1e-2)**2 * 0.5
    assert abs(o[5] / exp - 1.0) > 0.5


# ---------------------------------------------------------------- スピン (track6 の経路、POL + RAD)
def _tbmt(c, s, L, field, p0mc):
    """運動方程式と Thomas-BMT の数値積分 (test_v2_C2_spin と同じ)。c = (x, 力学的 px, y, 力学的 py, z, δ)、field(x, y, s) → B (1/m、p0 で規格化)"""
    pt = 1 + c[5]
    pv = np.array([c[1], c[3], math.sqrt(pt * pt - c[1] ** 2 - c[3] ** 2)])
    gam = math.sqrt(1 + (p0mc * pt) ** 2)

    def rhs(l, u):
        r, p, S = u[:3], u[3:6], u[6:]
        B = field(r[0], r[1], r[2])
        ph = p / pt
        Bl = np.dot(B, ph) * ph
        Om = ((1 + GS * gam) * (B - Bl) + (1 + GS) * Bl) / pt
        return np.concatenate([ph, np.cross(p, B) / pt, np.cross(S, Om)])
    ev = lambda l, u: u[2] - L
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2 * L + 1), np.concatenate([[c[0], c[2], 0.0], pv, s]), events=ev, rtol=1e-12, atol=1e-15, method='DOP853')
    return sol.y_events[0][0]


_BASE = '''MOMENTUM = %s GEV;
%s
SOL SA = (L=0 BZ=%s BOUND=1 GEO=1 DISFRIN=1) SB = (L=0 BZ=0 BOUND=1 DISFRIN=1);
MARK ST = (BX=5 BY=5);
LINE BL = (ST SA E1 SB);
FFS;
USE BL;
INS; NORAD; NOFLUC; CALC;
Exit;
'''
_P = np.array([[0, 0, 0, 0, 0, 0], [1e-3, 2e-4, -1e-3, 1e-4, 0, 1e-3], [-2e-3, -3e-4, 5e-4, 4e-4, 0, -5e-3]], float)
_S0 = ([0, 0, 1], [1, 0, 0], [0, 1, 0])


def _track(ffs, c, k0, k1, spin):
    c = np.asarray(c, float)
    n = len(c)
    cols = [list(c[:, j]) for j in range(6)]
    ffs.flags['RAD'] = True
    ffs.flags['FLUC'] = False
    ffs.flags['POL'] = True
    try:
        sp = np.asarray(spin, float)
        out = ffs.fn_track([[float(k0), cols + [[math.atan2(s[2], s[0]) for s in sp], [s[1] for s in sp], [1.0] * n]], float(k1)])
    finally:
        ffs.flags['POL'] = False
    r = out[1]
    so = [[math.sqrt(max(0.0, 1 - y * y)) * math.cos(a), y, math.sqrt(max(0.0, 1 - y * y)) * math.sin(a)] for a, y in zip(r[6], r[7])]
    return np.array(r[:6]).T, np.array(so)


@pytest.mark.parametrize('p,bzT', [(0.1, 1.5), (0.1, 3.0), (0.3, 3.0)])
def test_spin_in_sol_drift_at_large_bz_matches_thomas_bmt(v2_mode, p, bzT):
    """領域内の DRIFT (0.4 m) のスピン (POL + RAD): bz = 4.5・9 (0.1 GeV、BZ = 1.5・3 T)・3 (0.3 GeV) 1/m で Thomas-BMT の積分と 2e-6。
    軸を外れた粒子 (δ = 1e-3、−5e-3)、3 つの向きのスピン。旧版 (DRIFT が 1 刻み) は 1.5e-3・8.7e-3・1.3e-4 ずれた"""
    ffs, _ = _run(_BASE % (p, 'DRIFT E1 = (L=0.4);', bzT))
    bz = bzT * CL / (p * 1e9)
    fld = lambda x_, y_, s_: np.array([0.0, 0.0, bz])
    worst = 0.0
    for s0 in _S0:
        b, s = _track(ffs, _P, 1, len(ffs.els) + 1, [s0] * len(_P))
        for c, sv in zip(_P, s):
            ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
            u = _tbmt(ck, s0, 0.4, fld, p * 1e9 / ME)
            worst = max(worst, float(np.abs(sv - u[6:]).max()))
    assert worst < 2e-6, worst


def test_spin_in_sol_drift_at_large_bz_k64_is_coarse(k64_mode):
    """k64: 同じ DRIFT のスピンは 1 刻みで、bz = 9 1/m では 8.7e-3 ずれる (記録)"""
    p, bzT = 0.1, 3.0
    ffs, _ = _run(_BASE % (p, 'DRIFT E1 = (L=0.4);', bzT))
    bz = bzT * CL / (p * 1e9)
    fld = lambda x_, y_, s_: np.array([0.0, 0.0, bz])
    c = _P[1]
    b, s = _track(ffs, [c], 1, len(ffs.els) + 1, [[1, 0, 0]])
    ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
    u = _tbmt(ck, [1, 0, 0], 0.4, fld, p * 1e9 / ME)
    assert float(np.abs(s[0] - u[6:]).max()) > 1e-3


_RING = '''MOMENTUM = %(p)s GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.25);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=0.8e6 HARM=64);
MARK  ORIGIN = ();
SOL SA = (L=0 BZ=%(bz)s BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
%(extra)s
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA %(region)s SB);
LINE RINGQ = (ORIGIN 12*CELLQ CELLS 11*CELLQ RFC);
FFS;
USE RINGQ;
CELL; CALC;
NORFSW; NORAD; NOFLUC; POL;
e = Emittance[];
Print["SPIN ", SpinTune/.e];
Exit;
'''


def _spin_tune(region, extra, bzT, p):
    ffs, out = _run(_RING % dict(p=p, bz=bzT, region=region, extra=extra))
    line = [l for l in out.splitlines() if l.startswith('SPIN')][0]
    return float(line.split()[1])


@pytest.mark.parametrize('p,bzT', [(0.5, 1.0), (0.1, 1.5), (0.1, 3.0)])
@pytest.mark.parametrize('name,region,extra,L', [
    ('drift', 'DS DS', '', 0.5),
    ('mult', 'DS MM DS', 'MULT MM=(L=0.5 K1=0.0 DISFRIN=1);', 1.0),
    ('mult_k1', 'DS MM DS', 'MULT MM=(L=0.5 K1=0.05 DISFRIN=1);', 1.0),
    ('cavi', 'DS CC DS', 'CAVI CC=(L=0.5 VOLT=0 HARM=0 FREQ=5e8);', 1.0),
])
def test_emittance_pol_spin_tune_with_sol_region_equals_analytic(v2_mode, p, bzT, name, region, extra, L):
    """POL ON の Emittance[] の SpinTune (平面リング 24 台の偏向磁石、0.5 GeV: aγ = 1.135、0.1 GeV: aγ = 0.227): 領域 (DRIFT・MULT・K1 つき MULT・
    VOLT=0 の CAVI、長さ L の合計) が cos(μ/2) = cos(π ν0) cos(φ/2)、φ = (1+a) bz L と 1e-9 で一致する。
    ν0 は同じリングの BZ=0 の SpinTune (このリングの閉軌道 (周長の端数) が決める値で、aγ との差は 2.7e-5 (0.5 GeV)・3.9e-4 (0.1 GeV))。
    bz = 0.6・4.5・9 1/m。領域内の MULT・CAVI の POL (Emittance[]) の確認"""
    nu0 = _spin_tune(region, extra, 0.0, p)
    nu = _spin_tune(region, extra, bzT, p)
    bz = bzT * CL / (p * 1e9)
    phi = (1 + GS) * bz * L
    half = math.cos(math.pi * nu0) * math.cos(phi / 2)
    ana = math.acos(half) * 2 / (2 * math.pi)
    assert min(abs(abs(nu) - ana), abs(abs(nu) - (1 - ana))) < 1e-9, (nu, ana)


# ---------------------------------------------------------------- 加速 + RAD
@pytest.mark.parametrize('bz', [3.0, 9.0])
@pytest.mark.parametrize('volt', [0.0, 1e6, 2e7, 1e8])
def test_rf_mult_acceleration_does_not_change_radiation_loss_in_sol_region(v2_mode, bz, volt):
    """領域内の RF つき MULT (VOLT、1.3 GHz、位相 −1.2、基準運動量が +0 〜 +3.1% 上がる: TRPT) の放射損失 (RAD ON − OFF の δ の差) が
    −(2/3) r_e mc² (bz γ0² θ0)² L / E2 (θ0 = 1e-2 の横運動量、γ0 = 5871、3 GeV、L = 0.5 m) と 1e-3 で一致する。加速しても損失は変わらない
    (横運動量はソレノイドの縦磁場で変わらず、曲率 κ ∝ p⊥/p² なので γ⁴κ² = 一定)。E2 は出口の基準運動量"""
    P = 3.0e9
    el = _elem('MULT E1 = (L=0.5 VOLT=%g FREQ=1.3E9 PHI=-1.2 DISFRIN=1);' % volt, '3GEV')

    def run(rad):
        ph = _e._Phys(P, 100.0)
        col = _Col()
        x._compile_sol_v2(col, el, 'Multipole', bz, ph, 2)
        k, pr = col.items[0]
        S = tn.make_state(np.array([[0, 1e-2, 0, 0, 0, 0]], float), ph)
        ph2 = x.step_kind(k, pr, S, 1, ph, 0.0, True, rad, False, tn.make_rs(1), tm=10.0)
        return tn.state_out(S, 1)[0], (ph2 or ph)
    on, ph2 = run(True)
    off, _ = run(False)
    g0 = math.sqrt(1 + (P / ME) ** 2)
    U = (2.0 / 3.0) * RE * ME * (bz * g0**2 * 1e-2)**2 * 0.5
    exp = -U / (ph2.p0 * ME)
    assert abs((on[5] - off[5]) / exp - 1.0) < 1e-3, (on[5] - off[5], exp)


def test_pol_with_trpt_acceleration_in_sol_matches_fast_tracking(v2_mode):
    """POL (スピン追跡、track6 の経路) + TRPT + RFSW + SOL 領域内の CAVI (VOLT≠0) の加速。以前は基準運動量が変わらず δ が壊れ (δ = 0.27)、明示のエラーにしていた。
    v2.4.0 D3 (B238): track6 の経路も空洞の後の基準運動量・bz の規格化を更新する。POL なしの高速追跡 (領域内の加速も対応) と座標が 1e-9 以内で一致"""
    src = '''MOMENTUM = 1 GEV;
DRIFT D1 = (L=1.0) DS=(L=0.3);
CAVI CA = (L=0.5 VOLT=2E8 FREQ=1.3E9 PHI=-1.2);
SOL SA=(L=0 BZ=1 BOUND=1 GEO=1) SB=(L=0 BZ=0 BOUND=1);
MARK M0 = (BX=8 BY=3);
LINE TL = (M0 D1 SA DS CA DS SB D1);
FFS;
USE TL;
TRPT; RFSW; INS; CALC; RAD; NOFLUC;
Exit;
'''
    ffs, _ = _run(src)
    n = len(ffs.els) + 1
    c = np.array([[1e-3, 2e-4, -1e-3, 1e-4, 0, 1e-3]])
    b, _s = _track(ffs, c, 1, n, [[0, 0, 1]])
    ffs.flags['POL'] = False
    out = ffs.fn_track([[1.0, [[float(c[0, j])] for j in range(6)] + [[1.0]]], float(n)])
    fast = np.array([o[0] for o in out[1][:6]])
    assert abs(b[0, 5]) < 0.01 and np.abs(b[0] - fast).max() < 1e-9
