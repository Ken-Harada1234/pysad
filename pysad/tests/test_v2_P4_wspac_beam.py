# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 4: 追跡の WSPAC のビーム行列に、空間電荷つきの自己無撞着な平衡値 (MINCOUP の下限つき) を使うオプション `WSPACSC = 1;` (v2)。

K64 の追跡 (tturn の twspac) は `beamsize(:,l)` (`CODPLOT; WSPAC; EMIT;` の最後の一周のビーム行列) を使う。gs で確かめた: 結合のあるリング (`example/05_tracking/flag_WSPAC_tune_shift.sad` の 0.5 GeV)
の追跡のチューンシフトは (−2.938e-5、−1.210e-3) で、ゼロ電流の行列 (v2 の既定: −2.963e-5、−1.196e-3) と 1% で合う。**MINCOUP の下限 (εy ≥ MINCOUP (εx + εy)) を入れた行列
(−2.79e-5、−3.81e-4) とは合わない**。垂直のエミッタンスが 0 の平らなリングは gs が NaN (ビームの高さ 0)。したがって K64 の追跡のビーム行列には MINCOUP の下限が入らない。
v2 の既定はゼロ電流の行列のまま (gs と 1%) で、**`WSPACSC = 1;` のとき Emittance[] と同じ反復 (`touschek_table(wspac=True)`) の最終パス (`BeamMatricesSC`、MINCOUP の下限つき) を
`ffs.beam_matrix_at` が返す** (強い空間電荷でビームが膨らむ場合と、垂直のエミッタンスが 0 のリングで WSPACBEAM なしに使うため。K64 にはないオプション)。k64 は従来どおり
(追跡に WSPAC を入れない)。

根拠: (1) 行列の εy の下限 MINCOUP (εx + εy)、Σyy(1)/βy(1) = EMITYC。(2) 試験粒子の追跡のチューンシフトが、同じ行列の各要素の入口に置いた
薄いレンズ (前の要素の中点からこの要素の中点までの長さ) の一次の摂動の式 Δν = −(1/4π) Σ β_u k_u と 2e-4 の相対精度で一致 (PBUNCH = 1e6、0.2 GeV)。
使った条件: 0.2 GeV、12 セルの FODO (周長 45.6 m、RF 0.2 MV、HARM 64)、PBUNCH = 1e6〜1e8、MINCOUP = 0.5。実機の大きなリングでは確かめていない
(強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9) は反復が収束しない)。ビームの大きさは追跡の間じゅう固定 (強-弱)。
"""
import io
import math

import numpy as np
import pytest

_SRC = """MOMENTUM = 0.2GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (L=0 VOLT=0.2e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
USE RING; CELL; CALC;
RFSW; RAD; NOFLUC;
PBUNCH = %r; MINCOUP = 0.5;
"""


def _setup(pb):
    from pysad.tests.test_v2_C_spch import setup
    return setup(_SRC % pb)


def test_tracking_wspac_uses_self_consistent_beam_matrix_with_mincoup_when_wspacsc_v2(v2_mode):
    ffs, I, out = _setup(1e8)
    zero = ffs.beam_matrix_at(1)                                  # WSPAC OFF: ゼロ電流の平衡値 (εy = 0)
    assert zero[2, 2] < 1e-30
    I.run("WSPACSC = 1; WSPAC; e = Emittance[];", 'ffs')
    g = I.globals
    exc, eyc = float(g['EMITXC']), float(g['EMITYC'])
    assert eyc / exc == pytest.approx(0.5, abs=1e-9)             # MINCOUP の下限 εy = 0.5 (εx + εy)
    b1 = ffs.beam_matrix_at(1)
    assert b1[2, 2] / ffs.row(1).bety == pytest.approx(eyc, rel=1e-2)          # 最終パスの行列と収束値 (PBUNCH = 1e8 で 0.5%、1e6 で 2e-5)
    assert b1[4, 4] > 0.0 and b1[0, 0] > zero[0, 0]               # 空間電荷つきの行列 (σz も)


def test_tracking_wspac_tune_shift_equals_first_order_formula_with_the_same_matrix_v2(v2_mode):
    from pysad import emit6d
    from pysad.tests.test_v2_C2_spch import tunes
    pb = 1e6
    ffs, I, out = _setup(pb)
    I.run("NORAD;", 'ffs')
    t0 = tunes(I, None, d=1e-6, det_tol=2e-3)
    I.run("WSPACSC = 1; WSPAC; NORAD;", 'ffs')
    t1 = tunes(I, None, d=1e-6, det_tol=2e-3)                    # WSPACBEAM なし・WSPACSC = 1: Emittance の空間電荷つきの行列 (εy の下限つき)
    els = ffs.els
    pos = np.zeros(len(els) + 1)
    for i, el in enumerate(els):
        pos[i + 1] = pos[i] + el.L
    ph = emit6d._Phys(0.2e9, pos[-1])
    d = [0.0, 0.0]
    sprev = 0.0
    for l in range(len(els)):
        s = 0.5 * (pos[l] + pos[l + 1])
        al, sprev = s - sprev, s
        B = ffs.beam_matrix_at(l + 1)
        r = ffs.row(l + 1)
        sx, sy, sz = math.sqrt(B[0, 0] - B[0, 4] ** 2 / B[4, 4]), math.sqrt(B[2, 2] - B[2, 4] ** 2 / B[4, 4]), math.sqrt(B[4, 4])
        f = 2 * pb * emit6d.ELRADI * al / (ph.p0 ** 2 * ph.h0 * math.sqrt(2 * math.pi) * sz * (sx + sy))
        d[0] += -(1 / (4 * math.pi)) * r.betx * f / sx
        d[1] += -(1 / (4 * math.pi)) * r.bety * f / sy
    for i in range(2):
        assert (t1[i] - t0[i]) == pytest.approx(-d[i], rel=2e-3)
    assert -d[0] > 3e-5 and -d[1] > 7e-5                          # シフトは有限 (εy = 0 の行列だと σy = 0 でレンズが発散する)


def test_tracking_wspac_k64_does_not_apply_space_charge(k64_mode):
    """k64: WSPAC の追跡は今までどおり入れない (gs の追跡の WSPAC は K64 の beamsize で、pysad の k64 には未実装)"""
    from pysad.tests.test_v2_C2_spch import tunes
    ffs, I, out = _setup(1e6)
    I.run("NORAD;", 'ffs')
    t0 = tunes(I, None, d=1e-6, det_tol=2e-3)
    I.run("WSPAC; NORAD;", 'ffs')
    assert tunes(I, None, d=1e-6, det_tol=2e-3) == t0


def test_wspacsc_default_is_zero_current_matrix_like_gs_v2(v2_mode):
    """WSPACSC なし (既定): ゼロ電流の行列で、εy = 0 のリングは「WSPACBEAM を与えること」のエラー (K64 の追跡は平らなリングで NaN。`test_v2_C2_spch.py`)。
    WSPACSC = 1 なら MINCOUP の下限つきの行列で走る"""
    ffs, I, out = _setup(1e6)
    I.run("WSPAC;", 'ffs')
    assert ffs.beam_matrix_at(1)[2, 2] < 1e-30
    I.run("WSPACSC = 1;", 'ffs')
    assert ffs.beam_matrix_at(1)[2, 2] > 1e-12
