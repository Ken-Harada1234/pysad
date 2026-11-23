# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 8: 陽子 (電子以外) の放射の β の因子 (B112) の確認。

条件: MASS = 938.272088 MeV の陽子、12 セルの FODO (周長 45.6 m、偏向 15 度 × 12 = π)、運動量 1.5・3・30 GeV/c (γ = 1.89・3.35・32.0、β = 0.848・0.955・0.9995)、RF 2 MV・h = 24。
- k64: Emittance[] の EnergyLossU0 は gs と 9 桁一致 (K64 は古典半径に電子の値を使うので、陽子の U0 は mp/me = 1836 倍)。
- v2: 古典半径を陽子の値にする (B112) と U0 の比は 1/mp·me になる。残る差は K64 の超相対論的な式の β の因子で、古典式 U0 = (2/3) r_p m c² β³ γ⁴ Σθ²/L に対し
  v2 の U0 はちょうど 1/β³ 倍 (6e-10 で一致)。1.5 GeV/c で 64%、3 GeV/c で 15%、30 GeV/c で 0.15% 大きい。**β の因子は v2 でも直していない** (限界: 下の定数が示す。
  放射が意味を持つ γ ≫ 1 の粒子では小さい。直すには tradkn の光子数・エネルギー分布の β 依存も入れる必要があり、独立な基準 (厳密な放射のスペクトル) が要るので要判断)。
"""
import math

import pytest

from . import _p5_util as U

ME, MP, RE = 0.51099895e6, 938.272088e6, 2.8179403262e-15
SUMTHETA2_OVER_L = 12 * (0.2617993877991494 ** 2 / 1.2)
CASES = [('p1p5', 1.5e9), ('p3', 3.0e9), ('p30', 30.0e9)]


def _u0(text):
    return U.floats([l for l in text.splitlines() if l.startswith('U0 ')][0][3:])[0]


def _classical(p):
    e = math.hypot(p, MP)
    g, b = e / MP, p / e
    return (2.0 / 3.0) * (RE * ME / MP) * MP * b ** 3 * g ** 4 * SUMTHETA2_OVER_L, b


@pytest.mark.parametrize('name,p', CASES)
def test_k64_u0_matches_gs_and_is_mp_over_me_times_too_large(name, p):
    src = U.example_src('proton_beta_' + name)
    gs = _u0(U.gs_text('proton_beta_' + name))
    k = _u0(U.run(src, 'k64'))
    assert k == pytest.approx(gs, rel=1e-8)
    v = _u0(U.run(src, 'v2'))
    b = _classical(p)[1]
    assert k / v == pytest.approx(MP / ME / b ** 3, rel=2e-6)      # v2.4.0 D (B167): v2 は古典半径 (mp/me) と β³ の両方を直した


@pytest.mark.parametrize('name,p', CASES)
def test_v2_u0_equals_classical_times_one_over_beta_cubed(name, p):
    """v2.4.0 D (B167): v2 の U0 は古典式 (2/3) r_p m c² β³ γ⁴ Σθ²/L と 2e-6 で一致する (以前は 1/β³ 倍だった。関数名は履歴のまま)"""
    v = _u0(U.run(U.example_src('proton_beta_' + name), 'v2'))
    u, b = _classical(p)
    assert v / u == pytest.approx(1.0, rel=2e-6)
