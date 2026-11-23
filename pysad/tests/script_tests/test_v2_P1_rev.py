# -*- coding: utf-8 -*-
"""v2next 担当 P1 項目 3 (A2 の P1 の続き): 反転した BEND (曲率の符号が違う) のあるリングの P・τ が独立な式と合うこと、不安定なリングの扱い。

独立な式 (A2): τ = 1/[(5√3/8) r_e ƛ_C c γ⁵ ⟨|ρ|⁻³⟩]、P = (8/5√3) ⟨ρ⁻³⟩/⟨|ρ|⁻³⟩ (ρ は符号つきの曲率半径、⟨⟩ は周長での平均。偏向磁石の ANGLE と K0 から)。
リング: 小さな FODO (3 GeV) に ANGLE = +0.05, −0.10, +0.05 rad の組 (シケイン。曲げ角の和は 0) を入れたもの、大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)・
PF リング (2.5 GeV、187 m)・ASTRID2 (0.58 GeV、45.7 m) にも同じ組を入れたもの。
不安定なリング (1 周の写像の固有値 > 1) は v2 が P = NaN を返す (線形の式の周期解が不安定な軌道の上の値になるので、以前は 0.92 などを返していた)。
"""
import io
import math
import os

import numpy as np
import pytest

from pysad import emit6d, ffscmd, spin

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RE = 2.8179403262e-15
LC = 3.8615926796e-13
C = 299792458.0
ME = 0.51099895e6
PST = 8.0 / (5.0 * math.sqrt(3.0))


def _st_formula(els, p_eV):
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


_CHI = '''DRIFT DCHI = (L=0.1) DRCHI = (L=0.2);
BEND CB1 = (L=0.2 ANGLE=0.05) CB2 = (L=0.2 ANGLE=-0.10) CB3 = (L=0.2 ANGLE=0.05);
LINE CHI = (CB1 DCHI CB2 DCHI CB3 DRCHI);
'''

_FODO = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE CELLC = (QF D2 B D2 QD CHI);
'''


def _run(src, ring):
    full = src + f'\nFFS;\nUSE {ring};\nCELL; CALC; RFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(full, out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    return ffs, e['spin']


@pytest.mark.parametrize('line,tol', [('LINE RING = (ST 5*CELL1 CELLC 6*CELL1 RFC);', 1e-3)])
def test_small_ring_with_reversed_bends_matches_independent_formula_v2(v2_mode, line, tol):
    """シケインのある小さな FODO (3 GeV): P と τ が独立な式と 1e-3 (実測 1e-4 以下)。反転した曲げの ρ⁻³ が相殺しても、緩和率は |ρ|⁻³"""
    ffs, sp = _run(_FODO + _CHI + line, 'RING')
    tau, pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert pst < PST - 0.05
    assert abs(sp['equpol'][0] / pst - 1) < tol, (sp['equpol'][0], pst)
    assert abs(sp['taup'] / tau - 1) < tol, (sp['taup'], tau)
    assert sp['depol_rate'] == pytest.approx(0.0, abs=1e-20)


@pytest.mark.parametrize('path,ring,tol', [
    (os.path.join(HERE, '..', 'kek-ls_v03_12.lat'), 'RINGC', 2e-3),
    (os.path.join(ROOT, 'example', 'pf', 'lattice', 'v9_16_asym.sad'), 'RINGC', 2e-3),
    (os.path.join(ROOT, 'example', 'astrid', 'lattice', 'ASTRID2_v01_01.sad'), 'RINGC', 2e-3),
])
def test_large_rings_with_reversed_bend_chicane_match_independent_formula_v2(v2_mode, path, ring, tol):
    """大きなリング A・PF・ASTRID2 の先頭にシケインを入れたリングの P・τ が独立な式と 2e-3 (大きなリングでも 1e-3 台)"""
    if not os.path.exists(path):
        pytest.skip('lattice がない')
    src = open(path, encoding='utf-8').read() + '\n' + _CHI + '\nLINE RINGC = (CHI RING);\n'
    ffs, sp = _run(src, ring)
    tau, pst = _st_formula(ffs.lat.elements, ffs.momentum())
    assert abs(sp['equpol'][0] / pst - 1) < tol, (sp['equpol'][0], pst)
    assert abs(sp['taup'] / tau - 1) < tol, (sp['taup'], tau)


_VD = '''MOMENTUM = 0.5 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DV=(L=0.3);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
BEND VB1 = (L=0.4 ANGLE=0.15 ROTATE=90 DEG) VB2 = (L=0.4 ANGLE=-0.15 ROTATE=90 DEG);
CAVI RFC = (VOLT=0.2e6 HARM=64);
MARK ST = ();
LINE CELLQ = (QF D2 B D2 QD D1);
LINE VDOG = (VB1 DV VB2);
LINE RING = (ST 6*CELLQ VDOG 6*CELLQ RFC);
'''


def test_unstable_optics_with_finite_closed_orbit_gives_nan_v2(v2_mode):
    """垂直の偏向の組が鉛直方向の光学を不安定 (1 周の写像の固有値 1.8) にするリング: 閉軌道は見つかる (0) が、P は NaN で depol_error に理由が残る
    (v2.1.0 までは線形の式がそのまま P = 0.7 などを返していた)。安定な同じ構成 (QD = −0.55、±0.1 rad) は P が出る"""
    ffs, sp = _run(_VD, 'RING')
    assert np.isnan(sp['equpol'][0]) and 'unstable' in sp['depol_error']
    ffs, sp = _run(_VD.replace('K1=-0.6', 'K1=-0.55').replace('0.15', '0.1').replace('-0.15', '-0.1'), 'RING')
    assert np.isfinite(sp['equpol'][0]) and 'depol_error' not in sp


def test_rf_off_ring_uses_the_linear_formula_without_error_v2(v2_mode):
    """RF なし (NORFSW): 1 周の写像に縦の減衰が無い (νs = 0) ので側帯の式は使わず、線形の式 (depol_rate_linear) を使う。エラーも NaN も出ない"""
    src = _FODO + 'QUAD QFE = (L=0.3 K1=0.7 DY=0.003);\nLINE C2 = (QFE D2 B D2 QD D1);\nLINE RING = (ST 3*CELL1 C2 8*CELL1 RFC);\n'
    full = src + '\nFFS;\nUSE RING; CELL; CALC; NORFSW; RAD; NOFLUC; POL;\nExit;\n'
    ffs, _ = ffscmd.run_sad_string(full, out=io.StringIO())
    e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=False, radcod=True, pol=True)
    sp = e['spin']
    assert np.isfinite(sp['equpol'][0]) and 'depol_error' not in sp
    assert sp['depol_rate'] > 0 and not sp.get('depol_sideband')
