# -*- coding: utf-8 -*-
"""v2.5.0 (B341): Polarization の側帯の和 (spin.sideband_depolarization) の打ち切りの補正。

和 Σ_{k<K} e^{-iμk} G M^k u は K = min(kmax, ln tol / ln|λ|max) 項で切る。放射減衰が弱いリング (ASTRID2: |λ|max = 0.9999927、tol 1e-7 に 2.2e6 周が要るのに
kmax = 3e5) では、打ち切った和の残り (e^{-K/τ} ≈ 0.1 の大きさ) の位相が 1e-6 のチューンの変化で回り、変調なしの線形の和そのものが厳密な値から 6% ずれた
(B326 が 6 極の線形化を 1e-5 だけ動かして、傾いた MPW の P を 0.8% 動かして見えた)。補正: 打ち切られたとき同じ K・同じ試料の線形の和 (rate_lin_sameK) を返し、
slim_depolarization は「変調あり − 変調なし (同じ K)」の差だけを厳密な線形の脱偏極率に足す。
独立な基準: 同じ関数で kmax を収束する値 (tol を満たす K) まで上げた和。
"""
import math

import numpy as np

from pysad import spin


def _rot(nu, lam):
    c, s = math.cos(2 * math.pi * nu), math.sin(2 * math.pi * nu)
    return lam * np.array([[c, s], [-s, c]])


def _ring(nux, lam=1 - 2e-5, nuy=0.2664, nuz=0.0015, nus=0.3154):
    """1 要素の模型のリング: 3 つの減衰つきの振動 (x・y・z)、n = y のまわりのスピン回転 (チューン nus)、スピン軌道結合 D (δ の列は n 方向 = スピンチューンの変調を含む)"""
    T = np.zeros((6, 6))
    T[0:2, 0:2] = _rot(nux, lam)
    T[2:4, 2:4] = _rot(nuy, lam)
    T[4:6, 4:6] = _rot(nuz, lam)
    D = np.random.default_rng(3).normal(size=(3, 6)) * np.array([1, 1, 1, 1, 0.0, 0.3])
    D[1, :] = 0.0
    D[1, 5] = 0.8
    B = np.diag([1e-14, 1e-14, 1e-14, 1e-14, 0.0, 4e-14])
    c, s = math.cos(2 * math.pi * nus), math.sin(2 * math.pi * nus)
    R = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    return [(T, B, R, D)]


N0 = np.array([0.0, 1.0, 0.0])
S0 = np.diag([1e-8, 1e-8, 1e-8, 1e-8, 1e-9, 1e-8])


def test_truncated_sum_is_flagged_and_difference_is_stable_v2(v2_mode):
    """|λ|max = 1 − 2e-5 (K の必要数 8e5) を kmax = 1.6e5 で切ると、線形の和が厳密な値 (kmax = 1e6) から数 % 動き、チューンを 1e-6 変えるだけで
    その値が変わる (打ち切りの残りの位相)。補正した値 (厳密な線形 + (変調あり − 変調なし、同じ K)) は収束した変調ありの値と 0.1% で合う"""
    raw_err = []
    for nux in (0.1725, 0.172501, 0.172502):
        secs = _ring(nux)
        lin_t = spin.sideband_depolarization(secs, N0, S0, nsamp=16, kmax=160000, modulate=False)
        mod_t = spin.sideband_depolarization(secs, N0, S0, nsamp=16, kmax=160000, modulate=True)
        lin_c = spin.sideband_depolarization(secs, N0, S0, nsamp=16, kmax=1000000, modulate=False)
        mod_c = spin.sideband_depolarization(secs, N0, S0, nsamp=16, kmax=1000000, modulate=True)
        assert mod_t['truncated'] and not mod_c['truncated']
        assert 'rate_lin_sameK' not in mod_c
        assert abs(mod_t['rate_lin_sameK'] / lin_t['rate'] - 1) < 1e-12       # 同じ K・同じ試料の線形の和
        raw_err.append(mod_t['rate'] / mod_c['rate'] - 1)
        corrected = lin_c['rate'] + (mod_t['rate'] - mod_t['rate_lin_sameK'])
        assert abs(corrected / mod_c['rate'] - 1) < 1e-3, (nux, corrected, mod_c['rate'])
    assert max(abs(e) for e in raw_err) > 0.03, raw_err            # 補正なしの打ち切った和は数 % ずれる (この試験が B341 を見張る)
    assert max(raw_err) - min(raw_err) > 0.03, raw_err             # チューンの 1e-6 の変化で変わる
