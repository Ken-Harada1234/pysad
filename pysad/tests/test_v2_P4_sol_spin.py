# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 1: SOL 領域内の QUAD・MULT・CAVI のスピン (POL + RAD)。

根拠 (gs は使わない): Thomas-BMT の数値積分 (`test_v2_B1_sol._tbmt`、DOP853、rtol 1e-12)。軸上と軸を外れた粒子 (δ = 1e-3、−5e-3)、3 つの向きのスピン。
直したもの (v2 だけ。k64 は K64 のまま):
  (a) 領域内の QUAD のスピン (track6 の POL 経路) の刻み: 放射・スピンの歳差は 1 刻みの前後の向きの変化 (弦) から出すので、bz·al/ndiv ≈ 0.2 では bz = 9 1/m で 9e-5
      ずれた (C2 の bpar の後)。bz·al/ndiv ≤ 0.02 に刻む (DRIFT・BEND・Emittance[] の solemit._tquade と同じ)。9e-5 → 1e-6。
  (b) 端の縦磁場の積分 ψ (QUAD の k x y、MULT の k0 y + k1 x y + …) は端のフリンジを入れるとき (DISFRIN=0) だけ。DISFRIN=1 (フリンジなし) のとき ψ を足すと、
      MULT の双極成分 (k0 y = 5e-5) の分だけ 1e-4 ずれた (bz に依らない)。Emittance[] (solemit) と追跡 (track6) の両方。1.6e-4 → 1e-6。
使った条件: 1 GeV (bz = 0.3 1/m)・0.3 GeV (3)・0.1 GeV (4.5・9)、長さ 0.4 m の QUAD (K1 = 0.3)・MULT (K0 = 0.02、K1 = 0.3)、DISFRIN = 1。
未確認: DISFRIN = 0 (端のフリンジあり) のとき ψ と軌道のフリンジ写像を合わせた厳密な基準 (硬い端の B_s の取り込み) は作っていない。ずれ・ROTATE の MULT、大きなリング。
"""
import math

import numpy as np
import pytest

from pysad.tests.test_v2_B1_sol import CL, ME, _BASE, _P, _S0, _run, _tbmt, _track

_CASES = {
    'quad': ('QUAD E1 = (L=0.4 K1=0.3 DISFRIN=1);', lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.75 * x, bz]))),
    'mult': ('MULT E1 = (L=0.4 K0=0.02 K1=0.3 DISFRIN=1);', lambda bz: (lambda x, y, s: np.array([0.75 * y, 0.05 + 0.75 * x, bz]))),
}
_ENERGY = [(1.0, 1.0), (0.3, 3.0), (0.1, 1.5), (0.1, 3.0)]            # (GeV, BZ [T]) → bz = 0.3・3.0・4.5・9.0 1/m


def _worst(case, p, bzT):
    el, fld = _CASES[case]
    ffs, _ = _run(_BASE % (p, el, bzT))
    bz = bzT * CL / (p * 1e9)
    f = fld(bz)
    worst = 0.0
    for s0 in _S0:
        b, s = _track(ffs, _P, 1, len(ffs.els) + 1, [s0] * len(_P))
        for c, sv in zip(_P, s):
            ck = [c[0], c[1] + 0.5 * bz * c[2], c[2], c[3] - 0.5 * bz * c[0], 0, c[5]]
            u = _tbmt(ck, s0, 0.4, f, p * 1e9 / ME)
            worst = max(worst, float(np.abs(sv - u[6:]).max()))
    return worst


@pytest.mark.parametrize('p,bzT', _ENERGY)
@pytest.mark.parametrize('case', sorted(_CASES))
def test_spin_in_sol_region_quad_and_mult_at_large_bz_matches_thomas_bmt_v2(v2_mode, case, p, bzT):
    """領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD): bz = 0.3・3・4.5・9 1/m で Thomas-BMT の積分と 2e-6 以内"""
    assert _worst(case, p, bzT) < 2e-6


def test_spin_in_sol_region_quad_k64_is_coarser(k64_mode):
    """k64: 領域内の QUAD のスピンは K64 の bsi の式 (最初・最後のスライスを半分にして全体で 1 スライス分足りない。C2 で v2 は bpar に直した) で
    Thomas-BMT の積分と 6e-2〜0.19 ずれる (記録。bz = 0.3〜9 1/m)。MULT は k64 では領域内が未対応 (NotImplementedError)"""
    assert _worst('quad', 0.1, 3.0) > 1e-2


_STUDY = ['quad', 'quad_rot', 'mult', 'mult_rot']


@pytest.mark.parametrize('name', _STUDY)
def test_spin_orbit_coupling_and_rotation_in_region_are_tighter_than_before_v2(v2_mode, name):
    """強いソレノイド (BZ = 2 T、0.5 GeV、bz = 1.2 1/m) の中の QUAD (ROTATE つき)・MULT (K0 + K1、ROTATE つき) のスピンの軌道座標についての微分と軸上の回転 R
    (`test_v2_E_spinsol` と同じ手順): Thomas-BMT の積分に対し、追跡 1e-8 以内 (以前は 1e-6)、Emittance[] の D 2e-6 以内、R 1e-6 以内 (以前の許容は 3e-5。
    MULT の K0 の ψ を DISFRIN=1 で足していた分が 1e-4 あった)"""
    from pysad.tests.test_v2_E_spinsol import _ELEMENTS, _study
    el = [e for e in _ELEMENTS if e[0] == name][0]
    rows, R, r0 = _study(el[1], el[3], el[2], 2.0, 0.5)
    for dx_, dt, dm in rows:
        assert np.abs(dt - dx_).max() < 1e-8
        assert np.abs(dm - dx_).max() < 2e-6
    assert np.abs(R @ np.array([1.0, 0, 0]) - r0).max() < 1e-6


@pytest.mark.parametrize('volt', [0.0, 2e5])
def test_emittance_pol_spin_tune_with_accelerating_cavi_in_sol_region_agrees_with_tracking_v2(v2_mode, volt):
    """リング (12 セル、0.5 GeV、RFSW + RAD、POL) の 12 セル目の SOL 領域 (BZ = 2 T、bz = 1.2 1/m) の中の CAVI (L = 0.4、HARM = 64、VOLT = 0 と 0.2 MV):
    Emittance[] の SpinTune が閉軌道から追った 1 周のスピン行列の回転角と 5e-7 以内、スピン行列 srot1 が 1e-5。
    弱い試験: 0.5 GeV の U0 は小さく空洞の位相は 0 なので、SpinTune はエネルギー変化にほとんど鈍い (VOLT = 0 と 0.2 MV で同じ値)。TRPT の加速そのものは対象外"""
    import math
    from pysad.tests.test_v2_E_spinsol import _run as _run_e, _track as _track_e
    src = f'''MOMENTUM = 0.5 GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS=(L=0.5);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK  ORIGIN = ();
CAVI CV = (L=0.4 VOLT={volt} HARM=64 PHI=0.0);
SOL SA = (L=0 BZ=2 BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
LINE CELLQ = (QF D2 B D2 QD D1);
LINE CELLS = (QF D2 B D2 QD SA DS CV DS SB);
LINE RINGQ = (ORIGIN 12*CELLQ CELLS 11*CELLQ);
FFS;
USE RINGQ;
CELL; CALC;
RFSW; RAD; NOFLUC; POL;
Exit;
'''
    ffs, _ = _run_e(src)
    ffs.run_emit(quiet=True)
    sp = ffs.emit['spin']
    n = len(ffs.els) + 1
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _b, s = _track_e(ffs, np.zeros((1, 6)), 1, n, spin_=[s0])
        M[:, j] = s[0]
    nu_tr = math.acos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))) / (2 * math.pi)
    assert abs(abs(sp['nup']) - nu_tr) < 5e-7, (abs(sp['nup']), nu_tr)
    assert np.abs(sp['srot1'] - M).max() < 1e-5


@pytest.mark.parametrize('region', ['DR QQ DR', 'DR MM DR', 'DR QQ DR MM DR'])
def test_spin_tune_with_quad_mult_in_sol_region_of_astrid2_agrees_with_tracking_v2(v2_mode, region):
    """ASTRID2 (0.58 GeV、周長 45.7 m、`example/astrid/lattice/ASTRID2_v01_01.sad`) の RF の前に SOL 領域 (BZ = 0.5 T、bz = 0.26 1/m) を置き、QUAD (K1 = 0.05)・
    MULT (K1 = 0.04、K2 = 0.5) を DRIFT と一緒に入れる (DISFRIN=1)。Emittance[] の SpinTune が追跡の 1 周のスピン行列の回転角と 3e-7 以内、srot1 が 1e-5。
    弱い条件 (閉軌道がほぼ 0 なので ψ や軌道との結合は効かない)。確認できる範囲: 領域の長さ 0.9 m までの ASTRID2"""
    import math
    from pathlib import Path
    from pysad.tests.test_v2_E_spinsol import _run as _run_e, _track as _track_e
    p = Path(__file__).resolve().parents[2] / 'example' / 'astrid' / 'lattice' / 'ASTRID2_v01_01.sad'
    if not p.exists():
        pytest.skip('ASTRID2 のラティスがない')
    src = p.read_text(encoding='utf-8', errors='replace') + f'''
DRIFT DR = (L=0.3);
QUAD QQ = (L=0.3 K1=0.05 DISFRIN=1);
MULT MM = (L=0.3 K1=0.04 K2=0.5 DISFRIN=1);
SOL SA = (L=0 BZ=0.5 BOUND=1 GEO=1) SB = (L=0 BZ=0 BOUND=1);
LINE RSOL = (MKS 3*NCELL SA {region} SB CAVI 3*NCELL);
FFS;
USE RSOL;
CELL; CALC;
NORFSW; NORAD; NOFLUC; POL;
Exit;
'''
    ffs, _ = _run_e(src)
    ffs.run_emit(quiet=True)
    sp = ffs.emit['spin']
    n = len(ffs.els) + 1
    M = np.zeros((3, 3))
    for j, s0 in enumerate(([1, 0, 0], [0, 1, 0], [0, 0, 1])):
        _b, s = _track_e(ffs, np.zeros((1, 6)), 1, n, spin_=[s0])
        M[:, j] = s[0]
    nu_tr = math.acos(max(-1.0, min(1.0, (np.trace(M) - 1) / 2))) / (2 * math.pi)
    assert abs(abs(sp['nup']) - nu_tr) < 3e-7, (abs(sp['nup']), nu_tr)
    assert np.abs(sp['srot1'] - M).max() < 1e-5
