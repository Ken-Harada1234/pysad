# -*- coding: utf-8 -*-
"""v2.5.0: ビームストラールングの光子エネルギーの逆引き表の累積 (docs/30 B335)。

K64 `bbstrhl.f` の `mkbesseltab` は、光子数のスペクトル dN/dξ ∝ Sw(ξ) = ∫_ξ^∞ K_{5/3} の累積 SwI を、格子 ξ_i = 1.1^40/1.02^i で
SwI(i) = SwI(i−1) + Sw(ξ_i)(ξ_{i−1} − ξ_i) (区間の下端の Sw = 矩形則) と作り、ξ_min = 1.1e-7 より下の光子 (全体の 0.6%) を捨てて SwI(1000) = 1 に規格化する。
乱数 u (一様) → ln ξ を SwI で線形補間して光子のエネルギー cubs·ρ·ξ を決める。表から (線形補間込みで) 求めたモーメントは、理論値 ⟨ξ⟩ = 0.30792、⟨ξ²⟩ = 0.40741 に対し
+1.60%・+2.62% (台形則にすると +0.62%・+0.65%、残りは捨てた裾による規格化)。v2 は台形則 + 裾 (小 x の漸近形 3·ξ_min·Sw(ξ_min)) を全体に入れ、+0.016%・+0.046%。
独立な基準: scipy の K_{5/3} (kv) の数値積分による理論値。gs の保存した出力 (scripts/v250_bstrl_table__big.{sad,gs.txt}: BEAMBEAM の BSTRL = 1、1.7 GeV、4000 粒子 × 7 周で 422 個の光子の PhotonList のエネルギー)
と k64 の光子のエネルギーが 1e-6 (相対) 以内で一致する (K64 の表と同じ表)。
限界: 表のモーメント (線形補間込みの厳密な積分) と、422 個の光子のエネルギーの gs との一致。gs 上の統計 (モーメントの標本平均) は取っていない
(表が gs と同じであることは、光子ごとのエネルギーの一致で確かめた)。スピンの追跡 (B336) とは別。
"""
import io
import math
import os

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import kv

from pysad import beambeam as B, compat
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def theory():
    sw = lambda x: quad(lambda t: kv(5 / 3., t), x, np.inf, limit=200)[0]
    den = quad(sw, 0, np.inf, limit=200)[0]
    m1 = quad(lambda x: x * sw(x), 0, np.inf, limit=200)[0] / den
    m2 = quad(lambda x: x * x * sw(x), 0, np.inf, limit=200)[0] / den
    return m1, m2


def table_moments(xi, swi, nbins):
    """K64 の逆引き (ln ξ が SwI で区間ごとに線形) から、u 一様の ⟨ξ⟩、⟨ξ²⟩ を区間ごとに厳密に積分する"""
    out = [0.0, 0.0]
    for j in range(1, nbins):
        du = swi[j + 1] - swi[j]
        a, b = xi[j], xi[j + 1]
        for k in (1, 2):
            out[k - 1] += (math.exp(k * b) - math.exp(k * a)) / (k * (b - a)) * du
    return out


def table(mode):
    compat.set_mode(mode)
    B._TABLE = None
    B._TABLE_V2 = None
    return B.besseltab()


def test_theory_values():
    m1, m2 = theory()
    assert m1 == pytest.approx(0.307920, abs=2e-6) and m2 == pytest.approx(0.407407, abs=2e-6)


def test_k64_table_moments_are_high_by_1p6_and_2p6_percent():
    """k64 (K64 の表): ⟨ξ⟩ が +1.60%、⟨ξ²⟩ が +2.62% (理論比)"""
    m1, m2 = theory()
    xi, swi = table('k64')
    t = table_moments(xi, swi, B.LBTAB)
    assert t[0] / m1 - 1 == pytest.approx(0.0160, abs=3e-4)
    assert t[1] / m2 - 1 == pytest.approx(0.0262, abs=3e-4)
    assert swi[B.LBTAB] == 1.0


def test_v2_table_moments_match_theory():
    """v2: ⟨ξ⟩・⟨ξ²⟩ が理論値と 0.1% 以内 (+0.016%・+0.046%)。SwI は単調で SwI(1000) < 1 (裾の 0.6%)、格子は K64 と同じ"""
    m1, m2 = theory()
    x6, s6 = table('k64')
    xi, swi = table('v2')
    t = table_moments(xi, swi, B.LBTAB + 1)
    assert abs(t[0] / m1 - 1) < 1e-3 and abs(t[1] / m2 - 1) < 1e-3
    assert np.all(np.diff(swi[1:B.LBTAB + 2]) >= 0.0)
    assert 0.99 < swi[B.LBTAB] < 0.998 and swi[B.LBTAB + 1] == 1.0
    assert np.array_equal(xi[1:B.LBTAB + 1], x6[1:B.LBTAB + 1])


def _energies(txt):
    return np.array([float(l.split()[1].replace('D', 'E')) for l in txt.splitlines() if l.startswith('PA ')])


def test_k64_photon_energies_match_saved_gs():
    """k64: 4000 粒子 × 7 周の BSTRL の 422 個の光子のエネルギーが gs と 1e-6 (相対) 以内 (K64 の表と同じ表)"""
    src = open(os.path.join(SCRIPTS, 'v250_bstrl_table__big.sad'), encoding='utf-8').read()
    g = _energies(open(os.path.join(SCRIPTS, 'v250_bstrl_table__big.gs.txt'), encoding='utf-8').read())
    compat.set_mode('k64')
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    p = _energies(buf.getvalue())
    assert len(p) == len(g) == 422
    assert (np.abs(p - g) / np.abs(g)).max() < 1e-6


def test_v2_photon_energy_sum_is_lower_by_the_table_bias():
    """v2: 同じ乱数で、光子の数は k64 と同じ (422)、エネルギーの総和は k64 の 1/1.016 倍 (0.984 ± 0.002)"""
    src = open(os.path.join(SCRIPTS, 'v250_bstrl_table__big.sad'), encoding='utf-8').read()
    out = {}
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        B._TABLE = None
        B._TABLE_V2 = None
        buf = io.StringIO()
        run_sad_string(src, out=buf)
        out[m] = _energies(buf.getvalue())
    assert len(out['v2']) == len(out['k64']) == 422
    assert out['v2'].sum() / out['k64'].sum() == pytest.approx(0.9838, abs=2e-3)
