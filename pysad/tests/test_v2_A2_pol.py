# -*- coding: utf-8 -*-
"""担当 A2 (2): v2 の Polarization (線形のスピン軌道結合 + Sokolov-Ternov) の独立な検証。

根拠 (gs・追跡を使わない):
  1/τ_ST = (5√3/8) r_e ƛ_C c γ^5 <|ρ|^-3>、P_ST = (8/5√3) <ρ^-3>/<|ρ|^-3> (ρ は符号つきの曲率半径、<> は周長での平均)。
  偏向磁石の ANGLE と K0 から ∮ρ^-3 ds を自前で足す (光学・スピン追跡のコードを使わない)。
条件: 3 GeV の 24 セル FODO (周長 105.6 m、平面)、大きなリング A (kek-ls_v03_12.lat、3 GeV、周長 603 m)、PF リング (2.5 GeV、187 m、
example/pf/lattice/v9_16_asym.sad)、ASTRID2 (0.58 GeV、45.7 m) とその MPW 入り (example/astrid/lattice/ASTRID2_v01_02_id.sad)。
v2 の修正 (A2): 逆向きの磁場 (挿入光源の交互の極) のあるリングでは、K64 の sdamp = n·gintd は符号つきの Σ r_i (b_i·n) で、緩和率に使うと
P = 8/(5√3) (一定)・τ が実際より長く出る (ASTRID2 + MPW で τ が 31% 長い)。v2 は緩和率に Σ r_i (b_i·n)² (符号に依らない) を使う。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RE = 2.8179403262e-15
LC = 3.8615926796e-13
C = 299792458.0
ME = 0.51099895e6
PST = 8.0 / (5.0 * math.sqrt(3.0))


def _st_formula(els, p_eV):
    """偏向磁石 (ANGLE+K0 の曲げ) から ST の偏極時間と平衡偏極度 (独立な式)"""
    i3s = i3a = ctot = 0.0
    for e in els:
        ctot += e.L
        if type(e).__name__ == 'Bend' and e.L > 0:
            h = (e.ANGLE + getattr(e, 'K0', 0.0)) / e.L
            i3s += h ** 3 * e.L
            i3a += abs(h) ** 3 * e.L
    gam = math.sqrt(1 + (p_eV / ME) ** 2)
    rate = 5 * math.sqrt(3) / 8 * RE * LC * C * gam ** 5 * i3a / ctot
    return 1.0 / rate, PST * i3s / i3a


def _run(src_text, ring, extra=''):
    src = src_text + f'\nFFS;\nUSE {ring};\nCELL; CALC; RFSW; RAD; NOFLUC; POL;\n{extra}Exit;\n'
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    return ffs, e['spin']


def _read(path):
    if not os.path.exists(path):
        pytest.skip(os.path.basename(path) + ' がない')
    with open(path, encoding='utf-8') as fh:
        return fh.read()


_FODO = '''MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4.0E6 FREQ=5.0E8 PHI=0);
LINE CELL = (SF DS QF D2 B D2 QD DS SD D1);
LINE RING = (ORIGIN RF 24*CELL);
'''


def test_planar_fodo_ring_is_sokolov_ternov_v2(v2_mode):
    """平面リング: P = 8/(5√3) = 0.92376 (1e-12)、τ = 独立な式と 1e-6"""
    ffs, sp = _run(_FODO, 'RING')
    tau, pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert abs(sp["equpol"][0] - PST) < 1e-7
    assert abs(pst - PST) < 1e-12
    assert abs(sp['taup'] / tau - 1) < 1e-5, (sp['taup'], tau)
    assert sp['depol_rate'] == pytest.approx(0.0, abs=1e-20)


@pytest.mark.parametrize('name,path,ring,tol', [
    ('大きなリング A', os.path.join(HERE, 'kek-ls_v03_12.lat'), 'RING', 1e-4),
    ('PF リング', os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), 'RING', 1e-4),
    ('ASTRID2', os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RING', 1e-4),
])
def test_large_ring_polarization_time_matches_st_formula_v2(v2_mode, name, path, ring, tol):
    """大きなリング・PF・ASTRID2: 偏極時間が (5√3/8) r_e ƛ γ^5 c <ρ^-3> と 1e-4 以内 (実測 1e-6〜2e-5)、P = 92.376%"""
    ffs, sp = _run(_read(path), ring)
    tau, pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert abs(sp['taup'] / tau - 1) < tol, (name, sp['taup'], tau)
    assert abs(sp["equpol"][0] - PST) < 1e-7


def test_polarization_with_insertion_device_v2(v2_mode):
    """ASTRID2 + MPW (11 周期、逆向きの極): 偏極度と偏極時間が ∮ρ^-3 の符号つき・絶対値の比と 2e-3 以内 (P = 0.7019、τ = 21915 s)。
    K0 の極は ANGLE=0 の BEND (磁場だけ)"""
    path = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
    ffs, sp = _run(_read(path), 'RINGMPW')
    tau, pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert pst < 0.72 and abs(pst - 0.70163) < 1e-4          # 独立な式 (逆向きの極で下がる)
    assert abs(sp['equpol'][0] / pst - 1) < 2e-3, (sp['equpol'][0], pst)
    assert abs(sp['taup'] / tau - 1) < 2e-3, (sp['taup'], tau)


def test_polarization_with_insertion_device_k64(k64_mode):
    """k64: 偏極時間は符号つきの率 1/sdamp (28842 s) で、逆向きの極の分が相殺される (独立な式 21886 s より 31% 長い)"""
    path = os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_02_id.sad')
    ffs, sp = _run(_read(path), 'RINGMPW')
    tau, _pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert 1.25 < sp['taup'] / tau < 1.35


def test_large_ring_a_orbit_error_depolarization_scales_with_square_v2(v2_mode):
    """大きなリング A に鉛直方向の四極のずれ (rms 0.05 / 0.1 mm、60 台) を入れた脱偏極率: ずれの 2 乗に比例し (4 倍 ±10%)、P は 92.376% より下がる。
    (P1: rms 0.2 mm 以上は光学が不安定 (x の固有値 1.2) で、v2 は P = NaN を返す。A2・B3 の記録の 0.2・0.4 mm の値は意味がない)"""
    src = _read(os.path.join(HERE, 'kek-ls_v03_12.lat'))
    out = {}
    for amp in (0.5e-4, 1e-4):
        ffs, sp0 = _run(src, 'RING')
        rng = np.random.default_rng(3)
        quads = [q for q in ffs.lat.elements if type(q).__name__ == 'Quadrupole' and q.L > 0]
        for i in rng.choice(len(quads), 60, replace=False):
            quads[i].DY = float(rng.uniform(-amp, amp) * math.sqrt(3.0))
        ffs.calc(quiet=True)
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
        out[amp] = e['spin']
    r = out[1e-4]['depol_rate'] / out[0.5e-4]['depol_rate']
    assert 3.6 < r < 4.4, r
    assert out[0.5e-4]['equpol'][0] < PST and out[1e-4]['equpol'][0] < out[0.5e-4]['equpol'][0]
    assert out[1e-4]['equpol'][0] > 0.90
