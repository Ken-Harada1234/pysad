# -*- coding: utf-8 -*-
"""v2.5.0: WSPAC (K64 の twspac) の縦の力の電位の基準 (翻訳の誤りの修正。v2 の WSPAC の追跡。k64 は追跡に WSPAC を入れない)。

K64 の縦のキックは Δδ = −(z − zc)/σz²·az·fu(x, y)。fu = `twspfu` の電位は軸上で 0 ではなく ln((1 + r)/2) (r = σ2/σ1。`twspu` の
log(.5 + .5 r) の項。遠方は丸いビームの形に合わせた基準)。pysad (v2) は基準を軸上 0 としていた (「K64 と同じ」と誤解)。
gs の保存した出力 (scripts/v250_wspac_axis__ring.{sad,gs.txt}: 0.5 GeV の電子リング (8 セル、周長 35.3 m、弱いスキュー四極、RF 0.2 MV、
`example/05_tracking/flag_WSPAC_tracking_beam.sad` と同じ)、`CODPLOT; WSPAC; EMIT;` のあと PBUNCH = 2e10・4e10 の WSPAC ON と OFF の
1 周の追跡の差 (8 粒子: 軸上の z = σz・2σz・−σz、x = 1e-3・y = 1e-4 の z = −σz・σz・2σz、δ = 0.01、y = 0)):
軸上 (x = y = 0) の粒子も縦のキックを受け (z = σz で +9.69e-10、2σz で 4.33e-10、−σz で −9.69e-10。比は 2e^(−3/2) = 0.4463、PBUNCH に比例)、
大きさは fu = ln((1 + r)/2) ≈ −0.64 (r ≈ 0.06) の −z/σz²·az·fu と 1% で合う。
v2 (修正後) は 8 粒子の Δδ が gs の最大 Δδ の 2.5% 以内 (最大 Δδ の 1.5% 以内、ビーム行列がゼロ電流の値であることなどの分)、Δx・Δy が 3% 以内。修正前の v2 は軸上で Δδ = 0。
限界: 0.5 GeV・結合のある 1 つのリング・1 周・PBUNCH = 2e10・4e10。ビーム行列は v2 の既定 (ゼロ電流の平衡値)。Δz (1 周後の z の差) は gs の
約 2 倍で v2 と合わない (未調査。Δδ とは別)。壁の半径 (WSPACRADIUS、K64 にない) と K64 の電位の関係は physics_guide 17 章。
"""
import io
import os
import re

import numpy as np
import pytest

from pysad import compat
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def parse(txt):
    txt = re.sub(r'\n {1,}(?=[ \d.\-])', ' ', txt)
    out = {}
    for l in txt.splitlines():
        m = re.match(r'^(DD|DX|DY|DZ) \{(.*)\}', l)
        if m:
            out.setdefault(m.group(1), []).append(np.array([float(v.replace('D', 'E')) for v in m.group(2).split(',')]))
    return out


def _src():
    return open(os.path.join(SCRIPTS, 'v250_wspac_axis__ring.sad'), encoding='utf-8').read()


GS = parse(open(os.path.join(SCRIPTS, 'v250_wspac_axis__ring.gs.txt'), encoding='utf-8', errors='replace').read())


def _v2():
    compat.set_mode('v2')
    buf = io.StringIO()
    run_sad_string(_src(), out=buf)
    return parse(buf.getvalue())


def test_gs_axis_particles_get_a_longitudinal_kick():
    """gs: 軸上の粒子の縦のキックが z に奇で (比 2e^(−3/2))、PBUNCH に比例 (基準が軸上 0 ではない)"""
    d0, d1 = GS['DD'][0], GS['DD'][1]
    assert d0[0] == pytest.approx(9.69e-10, rel=0.01) and d0[2] == pytest.approx(-d0[0], rel=1e-3)
    assert d0[1] / d0[0] == pytest.approx(2 * np.exp(-1.5), rel=2e-3)
    assert d1[:3] == pytest.approx(2 * d0[:3], rel=1e-3)


def test_v2_matches_gs_longitudinal_kick():
    """v2: 8 粒子の Δδ が gs の最大 Δδ の 2.5% 以内 (軸上の粒子を含む。修正前は軸上で 0)"""
    p = _v2()
    for i in range(2):
        assert np.abs(p['DD'][i] - GS['DD'][i]).max() < 0.025 * np.abs(GS["DD"][i]).max(), i
    assert abs(p['DD'][0][0]) > 8e-10


def test_v2_matches_gs_transverse_kick():
    """v2: x = 1e-3・y = 1e-4 の粒子の Δx・Δy が gs と 3% 以内"""
    p = _v2()
    for k in ('DX', 'DY'):
        for i in range(2):
            a, b = p[k][i][3:], GS[k][i][3:]
            assert np.abs(a - b).max() / np.abs(b).max() < 0.03, (k, i)
