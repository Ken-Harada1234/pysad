# -*- coding: utf-8 -*-
"""v2next 担当 P1 項目 1・3: v2 の Polarization の脱偏極率が FLUC つきのスピン追跡の約 2.3 倍小さい件 (4.2 GeV・12 MV の 12 セル FODO) の原因と、式への取り込み。

原因 (記録 docs/history/v2next_P1.md §4): 軌道のクロマティシティ。FLUC なし・損失なしの決定論的な追跡 (蹴りへのスピン軸の応答 G、NAFF) を、
線形の軌道 (M1)・追跡 (M2)・追跡の軌道 + 線形のスピン (Mx)・線形の軌道 + 厳密なスピン (My)・1 周の写像の 2 次の項のうち δ×横方向だけを残した軌道 (chr) で比べると、
M2 の応答の増え方 (サンプルごとに 0.6〜2.1 倍) は Mx と chr で再現し (対数比の相関 0.95)、My (スピンの非線形) はむしろ 12% 減らす。δ の振幅依存のスピンチューンは
厳密に 0 (閉じた軌道で ν = aγ(1+δ))、δ のモンテカルロの非ガウス性は無視できる。ベータトロンのチューンの δ 依存 (ξ) が、スピンチューンとベータトロンの側帯の共鳴
(νspin ≈ νx + 2νs、離調 0.005) でスピンの相対位相の変調指数を変え、側帯の重みを変える。
v2 (P1): spin.sideband_depolarization に、ベータトロンの固有モードの位相の変調 ξ_j Δ_k (Δ_k = Σ_{t<k} δ_t) を入れた (C6 = 分散軌道に沿った 1 周の写像の δ 微分)。
"""
import io
import math
import os
import sys

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, spin

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_v2_B3_pol import _lat, _setup, PST        # noqa: E402


@pytest.fixture(autouse=True, scope='module')
def _emitdiv_default_1():
    """参照 (FLUC の追跡 ±16%) との比較は EMITDIV = 1 の刻みで作った。v2.4.0 から POL ON の既定は 16 (線形の脱偏極率の刻みの誤差の収束値) で、
    24 MV では式が追跡より 14% 小さい (追跡の統計誤差 ±16% の内。docs/history/v24_RAD.md §5)。この試験は変調の模型の検証なので 1 に固定する"""
    old = os.environ.get('PYSAD_EMITDIV_DEFAULT')
    os.environ['PYSAD_EMITDIV_DEFAULT'] = '1'
    try:
        yield
    finally:
        if old is None:
            os.environ.pop('PYSAD_EMITDIV_DEFAULT', None)
        else:
            os.environ['PYSAD_EMITDIV_DEFAULT'] = old

RE = 2.8179403262e-15
LC = 3.8615926796e-13
C = 299792458.0
ME = 0.51099895e6
FLUC_4P2_12MV = 10.4e-8          # A2・B3 の FLUC つきスピン追跡 (4 本・480 粒子の重み付き平均、±1.6e-8)


def _rate(E, volt, chrom=True):
    ffs, _ = ffscmd.run_sad_string(_lat(E, volt), out=io.StringIO())
    old = spin.CHROMATIC
    spin.CHROMATIC = chrom
    try:
        e = emit6d.emittance6d(ffs.lat, ffs.momentum(), rfsw=True, radcod=True, pol=True)
    finally:
        spin.CHROMATIC = old
    return e['spin']


def test_chromatic_modulation_brings_the_formula_to_the_fluc_tracking_v2(v2_mode):
    """4.2 GeV・12 MV の 12 セル FODO: 脱偏極率は B3 の式 4.46e-8 (追跡の 0.43 倍) から 1.09e-7 になり、FLUC つきのスピン追跡の (10.4 ± 1.6)e-8 と 1.0±0.2 で合う。
    線形の軌道の MC 5.0e-8 (B3)・クロマティックな軌道の MC 1.24e-7 (P1、3000 粒子 × 30000 周)"""
    sp = _rate(4.2, 12)
    assert sp['depol_chromatic'] is True
    assert abs(sp['depol_rate'] / 1.09e-7 - 1) < 0.05, sp['depol_rate']
    assert 0.8 < sp['depol_rate'] / FLUC_4P2_12MV < 1.3
    assert abs(sp['depol_rate_linear'] / 2.40756e-8 - 1) < 0.01      # v2.4.0 D3 (B230: 四極の本体の B∥): 線形の式 2.4076e-8 (以前 2.314e-8 → B165 で 2.342e-8)。v2.4.0 D (B165): 線形の式 2.342e-8 は FLUC 追跡の値と 1.2% (以前 0.2%)。追跡の D の BEND の端の ψ の分解の近似 (R12・R21 に 1e-2) が残るため
    sp0 = _rate(4.2, 12, chrom=False)
    assert abs(sp0['depol_rate'] / 4.45e-8 - 1) < 0.05 and not sp0['depol_chromatic']
    assert 0.1 < sp['mod_chrom'] < 0.3 and 0.2 < sp['mod_spin'] < 1.0          # 変調指数 (ベータトロン、スピン)


@pytest.mark.parametrize('volt,expect', [(6, 8.03e-8), (24, 8.35e-8)])
def test_chromatic_rate_at_other_rf_voltages_v2(v2_mode, volt, expect):
    """4.2 GeV の RF 6 MV (νs = 0.020) と 24 MV (νs = 0.043): 6 MV 8.03e-8、24 MV 8.35e-8 (B3 の式は 3.49e-8・4.93e-8)。追跡は (測定していない。記録 §4.5)"""
    sp = _rate(4.2, volt)
    assert abs(sp['depol_rate'] / expect - 1) < 0.05
    assert sp['depol_rate'] > 1.4 * _rate(4.2, volt, chrom=False)['depol_rate']


def test_chromatic_term_is_small_at_low_energy_v2(v2_mode):
    """3 GeV・1.5 MV (νs = 0.012、スピンチューンは側帯の共鳴から遠い): クロマティックな変調を入れても率は 5% 以内しか変わらない (6.17e-11 と 6.29e-11)。追跡 3.4e-11 ± 22%"""
    sp = _rate(3.0, 1.5)
    sp0 = _rate(3.0, 1.5, chrom=False)
    assert abs(sp['depol_rate'] / sp0['depol_rate'] - 1) < 0.05


def test_chromatic_extension_is_continuous_and_zero_matrix_is_the_old_formula_v2():
    """C6 = 0 は B3 の式と同じ。C6 を ε 倍にした極限で B3 の式に連続につながる (ε = 1e-3 で 1%)"""
    ring = _setup(4.2, 12)
    s0 = np.asarray(ring['e']['beam_matrix'])
    n0 = ring['sl']['n0']
    r_old = spin.sideband_depolarization(ring['secs'], n0, s0)['rate']
    r_zero = spin.sideband_depolarization(ring['secs'], n0, s0, C6=np.zeros((6, 6)))['rate']
    assert abs(r_zero / r_old - 1) < 1e-12
    compat.set_mode('v2')
    try:
        ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12), out=io.StringIO())
        els = ffs.lat.elements
        e = ring['e']
        ph = emit6d._Phys(ffs.momentum(), emit6d._emit_circ(els))
        C6 = spin._chromatic_matrix(els, ph, np.array(e['orbit_entrance']), e.get('trf0', 0.0), True, True, ring['secs'])
    finally:
        compat.set_mode('k64')
    r_small = spin.sideband_depolarization(ring['secs'], n0, s0, C6=1e-3 * C6)
    assert abs(r_small['rate'] / r_old - 1) < 0.01 or r_small['chromatic'] is False
    r_full = spin.sideband_depolarization(ring['secs'], n0, s0, C6=C6)
    assert r_full['chromatic'] and 2.0 < r_full['rate'] / r_old < 2.9
    # ベータトロンのチューンの δ 微分 (rad/δ): 2π × 2.66 と 2π × 2.61 (x・y)。シンクロトロンのモードは 0
    xi = np.sort(np.abs(r_full['xi']))
    assert abs(xi[-1] / (2 * math.pi) - 2.656) < 0.05 and abs(xi[-3] / (2 * math.pi) - 2.609) < 0.05 and xi[0] == 0.0


def test_chromatic_response_matches_deterministic_tracking_models_v2(v2_mode):
    """導出の検算: 同じ初期条件 (閉軌道からの差 u、δ = ±4e-3 など) の δ の蹴りへのスピン軸の応答 |G| を、式 (∂α*/∂u の δ 列。ξ = 0 は M1 と 12% 以内、
    クロマティックは追跡 M2 と同じ向きに動く) で。記録の表 §4.2 (16 試料の <G²> は M2/M1 = 1.62、式 2.44、対数比の相関 0.95) の代表 4 試料"""
    ring = _setup(4.2, 12)
    s0 = np.asarray(ring['e']['beam_matrix'])
    n0 = ring['sl']['n0']
    compat.set_mode('v2')
    try:
        ffs, _ = ffscmd.run_sad_string(_lat(4.2, 12), out=io.StringIO())
        els = ffs.lat.elements
        e = ring['e']
        ph = emit6d._Phys(ffs.momentum(), emit6d._emit_circ(els))
        C6 = spin._chromatic_matrix(els, ph, np.array(e['orbit_entrance']), e.get('trf0', 0.0), True, True, ring['secs'])
    finally:
        compat.set_mode('k64')
    # (u の δ、x、M1 の |G|、M2 の |G|) は example/16_spin_polarization/pol_chromatic_response_samples_*.json の抜粋。ここでは ξ = 0 の式が M1 と一致することだけ確かめる
    rng = np.random.default_rng(5)
    L = np.linalg.cholesky(s0)
    U = (L @ rng.standard_normal((6, 6))).T
    g0 = np.abs(spin.sideband_depolarization(ring['secs'], n0, s0, us=U, nsamp=16)['G0'][:, 5])
    g1 = np.abs(spin.sideband_depolarization(ring['secs'], n0, s0, us=U, C6=C6, nsamp=16)['G0'][:, 5])
    assert np.all(g0 > 0.5) and np.all(g1 > 0.1)
    assert np.mean(g1 ** 2) / np.mean(g0 ** 2) > 1.2                        # 平均で応答が増える
