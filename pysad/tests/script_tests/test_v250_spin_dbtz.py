# -*- coding: utf-8 -*-
"""v2.5.0: Emittance[] (POL) のスピン軌道結合の z 成分の微分 `dbtz = dxpzb/pr − dblz` (`photons.f:1550`) が pr² でなく pr で割る (docs/30 B336)。

`tradke` は btx = xpx/pr² − blx、bty = xpy/pr² − bly、btz = xpzb/pr² − blz で、その微分は dbtx = dxpx/pr² − dblx (δ の列に −2 xpx/pr³)、dbty も同じ、
dbtz = dxpzb/pr − dblz (δ の列に −2 xpzb/pr³ = pr⁻² の微分の形)。z 成分だけ pr で割るのは、自分の btz の微分として不整合 (誤植)。δ ≠ 0 の軌道 (pr ≠ 1) で、
スピン軌道結合 D = srot(:,4:9) の z の行 (縦の歳差) が δ に比例して誤る。δ = 0 では pr = 1 で同じ。
独立な基準 (pysad の Emittance[] の解析式を使わない): 追跡のスピン回転 (`track6.track_particles` の spins) の中心差分から作る D (`test_v2_A_spin._spin_orbit_d` と同じ手順)。
結果 (3 GeV/c、1 要素 + ドリフト 1 m、軌道 (x, px, y, py) = (8e-3, 1e-3, 6e-3, −1e-3) の周りの D): 偏向磁石 (K1 = 0)・δ = 0.05 で、解析式と追跡の差が最大値の相対で
K64 の式 13.2%、直した式 1.20% (δ = 0.10 で 26.9% 対 2.24%)。K1 ≠ 0 の偏向磁石は 6.3% 対 0.18% (δ = 0.10 で 12.4% 対 0.15%)。δ = 0 は同じ (0.8%、0.2%)。
gs では確認できない: K64 の Polarization は NaN (B55)、srot の微分列は Emittance[] の出力に出ない。k64 は K64 のまま (pysad の k64 の式は誤植をそのまま写している)。
限界: 3 つの要素 (BEND の K1 = 0・≠ 0、QUAD) のうち、δ の効果が見えるのは偏向磁石だけ (四極は 0.2% で δ に依らない)。追跡 (TrackParticles) のスピンは別のコードで影響しない。
v2 の δ = 0.1 の偏向磁石の 2.2% の残りは別の項 (未調査)。
"""
import io

import numpy as np
import pytest

from pysad import compat, emit6d, ffscmd, spin, track6

pytestmark = pytest.mark.filterwarnings('ignore')

ELS = {
    'bend0': 'BEND E1 = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);',
    'bendk': 'BEND E1 = (L=1.2 ANGLE=0.2617993877991494 K1=0.3 E1=0.1308996938995747 E2=0.1308996938995747);',
    'quad': 'QUAD E1 = (L=0.3 K1=0.7);',
}
OFF = [8e-3, 1e-3, 6e-3, -1e-3, 0.0]


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def spin_orbit_d(kind, cod):
    src = 'MOMENTUM = 3GEV;\nDRIFT D1 = (L=1);\n%s\nMARK ST = ();\nLINE RING = (ST E1 D1);\nFFS;\nUSE RING;\n' % ELS[kind]
    ffs, _ = ffscmd.run_sad_string(src, out=io.StringIO())
    el = [e for e in ffs.lat.elements if e.name == 'E1'][0]
    p = ffs.momentum()
    ph = emit6d._Phys(p, 100.0)
    old = emit6d.POL_MODE
    emit6d.POL_MODE = True
    try:
        st = emit6d.tturne([el], np.array(cod, float), ph, rfsw=False, radcod=True)
    finally:
        emit6d.POL_MODE = old
    de = st.srot[:, 3:9].copy()

    def trk(c0):
        _o, sp = track6.track_particles([el], np.tile(c0, (3, 1)), 1, p, 0.0, False, True, False, spins=np.eye(3), radpol=False, circ=100.0)
        return np.asarray(sp).T
    c = np.array(cod, float)
    r0 = trk(c)
    h = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-4, 1e-5])
    dt = np.zeros((3, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h[j]
        a = (trk(c + d) - trk(c - d)) @ r0.T / (2 * h[j])
        dt[:, j] = -np.array([a[2, 1] - a[1, 2], a[0, 2] - a[2, 0], a[1, 0] - a[0, 1]]) * 0.5
    return de, dt


def rel_err(kind, delta, mode, pr2=True):
    compat.set_mode(mode)
    old = spin.DBTZ_PR2
    spin.DBTZ_PR2 = pr2
    try:
        de, dt = spin_orbit_d(kind, OFF + [delta])
    finally:
        spin.DBTZ_PR2 = old
    return float(np.abs(de - dt).max() / np.abs(dt).max())


@pytest.mark.parametrize('kind,delta,fixed_max,typo_min', [('bend0', 0.05, 0.02, 0.10), ('bend0', 0.10, 0.03, 0.20),
                                                           ('bendk', 0.05, 0.01, 0.05), ('bendk', 0.10, 0.01, 0.10)])
def test_spin_orbit_coupling_at_nonzero_delta(kind, delta, fixed_max, typo_min):
    """δ ≠ 0 の偏向磁石 (v2 の他の修正はそのまま): dbtz を pr² で割ると D が追跡の差分と fixed_max 以内、pr で割る (K64 の式) と typo_min 以上ずれる"""
    assert rel_err(kind, delta, 'v2', True) < fixed_max
    assert rel_err(kind, delta, 'v2', False) > typo_min


def test_dbtz_switch_is_ignored_in_k64():
    """k64 は DBTZ_PR2 によらず K64 の式 (pr)"""
    assert rel_err('bend0', 0.05, 'k64', True) == rel_err('bend0', 0.05, 'k64', False)


@pytest.mark.parametrize('kind', ['bend0', 'bendk', 'quad'])
def test_zero_delta_is_unchanged_by_the_dbtz_fix(kind):
    """δ = 0 では pr = 1 で dbtz の割り算が同じ: 両方の割り算で D の追跡との差が 3e-4 (絶対) 以内で同じ (放射のエネルギー損失で pr が 1e-6 ずれる分)、1% 以内"""
    a = rel_err(kind, 0.0, 'v2', True)
    assert a == pytest.approx(rel_err(kind, 0.0, 'v2', False), abs=3e-4) and a < 0.01


def test_quad_is_insensitive_to_delta():
    """四極は δ ≠ 0 でも v2 と追跡の差が 0.2% (btz への寄与が小さい)"""
    assert abs(rel_err('quad', 0.10, 'v2') - rel_err('quad', 0.0, 'v2')) < 5e-4
