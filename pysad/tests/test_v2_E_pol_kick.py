# -*- coding: utf-8 -*-
"""v2next E 項目 3: 高エネルギー・高 RF 電圧の Polarization の脱偏極率が FLUC つきのスピン追跡の約 2.3 倍小さい件の切り分けの道具 (δ の蹴りへの応答 G = ∂n/∂δ) の検証。
記録は docs/history/v2next_E.md。長い測定 (追跡 M2 の振幅別の表) は example/v2_next/E_pol_kick_response.py (PYSAD_SLOW=1 でここでも 1 本走る)。"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _e_spinkick as K
from pysad import naff


def test_spin_line_amplitude_is_recovered_from_a_synthetic_cone():
    """ζ = 軌道の線 2 本 + スピンチューンの線 (複素振幅 a)。NAFF で a を 1e-6 以内で取り出し、近い軌道の線 (間隔 0.0053) に影響されない"""
    n = np.arange(1500)
    a = 3.0e-4 * np.exp(0.7j)
    z = a * np.exp(2j * np.pi * 0.2359 * n) + 1.2e-2 * np.exp(2j * np.pi * 0.2306 * n + 0.3j) + 5e-3 * np.exp(2j * np.pi * 0.2906 * n)
    fr, am = K.line_amp(z, 0.2359, 1500)
    assert abs(fr - 0.2359) < 1e-6
    assert abs(am - a) < 2e-6, abs(am - a)


def test_linear_model_response_at_zero_amplitude_equals_slim_n_delta_v2(v2_mode):
    """線形の軌道モデル (M1) の閉軌道上の δ の蹴りへの応答が、Emittance[] の SLIM のスピン軌道結合 N0 の δ 列 (n に垂直な成分) と 0.5% 以内で一致する
    (方法の検算。4.2 GeV・12 MV の 12 セル FODO、四極の縦ずれ ±3 mm)。追跡 (M2) は同じ量が 13% 小さい (2.36 対 2.73。追跡の刻みを 12 倍にしても同じ。要判断)"""
    try:
        ring = K.Ring(4.2, 12, 0)
        secs = K.build_secs(ring)
        h = 3e-5
        r = K.track_pair_linear(ring, secs, ring.cod.copy(), 1500, [h, -h], abs(ring.sp['nup']))
        G = abs((r[1][1] - r[2][1]) / (2 * h))
        N0 = np.array(ring.sp['spin_orbit'])
        ref = np.linalg.norm(N0[:, 5] - ring.n0 * (N0[:, 5] @ ring.n0))
        assert abs(G / ref - 1.0) < 5e-3, (G, ref)
    finally:
        K.restore_cave()


@pytest.mark.skipif(not os.environ.get('PYSAD_SLOW'), reason='追跡 3 粒子 x 1500 周を数回 (数分)')
def test_tracking_response_is_within_a_factor_of_the_linear_model_at_one_sigma_v2(v2_mode):
    """追跡 (非線形な軌道 + 厳密なスピン) の δ の蹴りへの応答 <|G|²> は、1σ の振幅 (ベータトロン + シンクロトロン) で線形の軌道モデルの 0.6〜1.6 倍。
    FLUC 追跡の脱偏極率が式の 2.3 倍という差 (B3) は、この応答の違いでは説明できない (記録 §4)"""
    try:
        ring = K.Ring(4.2, 12, 0)
        secs = K.build_secs(ring)
        rng = np.random.default_rng(3)
        h = 3e-5
        g1, g2 = [], []
        for _ in range(4):
            u = ring.sample('both', 1.0, rng)
            r1 = K.track_pair_linear(ring, secs, u, 1500, [h, -h], abs(ring.sp['nup']))
            r2 = K.track_pair(ring, u, 1500, [h, -h], abs(ring.sp['nup']))
            g1.append(abs((r1[1][1] - r1[2][1]) / (2 * h)) ** 2)
            g2.append(abs((r2[1][1] - r2[2][1]) / (2 * h)) ** 2)
        assert 0.6 < np.mean(g2) / np.mean(g1) < 1.6
    finally:
        K.restore_cave()
