# -*- coding: utf-8 -*-
"""v2.5.0: UND の z の更新に −dv·L が無い (docs/30 B334)。

K64 `undulator.f:109` の z の更新は −(fgx² + fgy²)dpds²/(2 ds) だけで、速度の差による −dv·ds の項が無い (DRIFT は持つ。`tdrift` の `z += dpz·al1 − dv·al`)。
磁場なしの UND (K64 の UND は磁場を設定できない: B37) は、DRIFT と同じはずなのに、δ ≠ 0 の粒子の z が dv·L だけ違う (R56 = L/γ² が欠ける)。
独立な基準: 磁場なしの自由飛行の z (SAD の符号: z > 0 は基準粒子より早い) は z = L (h0/p0 − h/p) (基準粒子と粒子の飛行時間の差 × c、p・h は mc 単位)。
gs の保存した出力 (scripts/v250_und_dv__{und,drift}_{3gev,10mev}.{sad,gs.txt}: UND (L = 1 m、KX = 0.1、QY = 0.2、POLE = 10、SLICE = 4、PHI = 30) と DRIFT (L = 1 m) の
TrackParticles の 5 粒子 (x = 1e-4、2e-4、px = 1e-4、2e-4、δ = 1e-3・1e-3・1e-2・−1e-2・0)): gs の UND − DRIFT の z は、3 GeV/c で δ = 1e-3 に −2.90e-11、δ = ±1e-2 に ∓2.9e-10、
10 MeV/c で δ = 1e-3 に −2.60e-6、δ = ±1e-2 に ∓2.6e-5 (= dv·L)。DRIFT の z は上の自由飛行の式と一致する。k64 は gs と全桁一致。v2 は UND の z が DRIFT と 1e-14 以内で一致する。
限界: 磁場なしの UND (gs では磁場を設定できない) だけを gs で確認。v2 の磁場つき UND は、磁場を 0 に近づけると DRIFT の z に一致することだけ確かめた
(磁場つきの z の更新の絶対値は、横の運動の経路の長さを含む K64 の式のまま)。放射つきの UND は −dv·ds を刻みごと (python 版) または最後に全長で (numba 版) 足す。
"""
import io
import math
import os
import re

import numpy as np
import pytest

from pysad import compat
from pysad.ffscmd import run_sad_string

pytestmark = pytest.mark.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, 'scripts')
MC2 = 0.51099895e6
DELTAS = np.array([1e-3, 1e-3, 1e-2, -1e-2, 0.0])


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def parse(txt):
    body = txt.split('R {', 1)[1]
    return np.array([[float(v.replace('D', 'E')) for v in n.split(',')] for n in re.findall(r'\{([^{}]*)\}', body)])


def gs(k, m):
    return parse(open(os.path.join(SCRIPTS, 'v250_und_dv__%s_%s.gs.txt' % (k, m)), encoding='utf-8').read())


def run(k, m, mode):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(open(os.path.join(SCRIPTS, 'v250_und_dv__%s_%s.sad' % (k, m)), encoding='utf-8').read(), out=buf)
    return parse(buf.getvalue())


MOMS = [('3gev', 3e9), ('10mev', 10e6)]


@pytest.mark.parametrize('m,mom', MOMS)
@pytest.mark.parametrize('k', ['und', 'drift'])
def test_k64_matches_saved_gs(k, m, mom):
    """k64 の 5 粒子の追跡が gs と 1e-15 以内"""
    assert np.abs(run(k, m, 'k64') - gs(k, m)).max() < 1e-15


@pytest.mark.parametrize('m,mom', MOMS)
def test_gs_und_lacks_dv_l_while_drift_has_it(m, mom):
    """gs: UND − DRIFT の z が、全粒子 (δ = 1e-3・1e-3・1e-2・−1e-2・0) で −L(h0/p0 − h/p) (自由飛行の式、−dv·L) と相対 2e-3 以内。δ = 0 の粒子は 0。
    DRIFT の x の傍軸近似の差 4e-12 は別 (UND は x' = px/pn)"""
    p0 = mom / MC2
    h0 = math.sqrt(1 + p0 * p0)
    p = p0 * (1 + DELTAS)
    h = np.sqrt(1 + p * p)
    exp = 1.0 * (h0 / p0 - h / p)
    gd, gu = gs('drift', m), gs('und', m)
    # 粒子ごとの δ は DELTAS (横の運動は UND と DRIFT で同じ大きさ。差を取ると横の項が消える)
    dz = gu[4] - gd[4]
    assert dz == pytest.approx(-exp, rel=2e-3, abs=1e-15)
    assert np.abs(dz).max() > 1e-11


@pytest.mark.parametrize('m,mom', MOMS)
def test_v2_und_equals_drift_in_z(m, mom):
    """v2: 磁場なしの UND の z が DRIFT と 1e-14 以内で一致 (δ = ±1e-2 を含む)。k64 は 3 GeV/c で 3e-10、10 MeV/c で 3e-5 違う"""
    zu = run('und', m, 'v2')[4]
    zd = run('drift', m, 'v2')[4]
    assert np.abs(zu - zd).max() < 1e-14
    assert np.abs(run('und', m, 'k64')[4] - run('drift', m, 'k64')[4]).max() > 2e-10


def test_v2_weak_field_und_approaches_drift_z():
    """v2: 磁場つき UND (FBX = 1e-9 T) の z が、磁場なしの DRIFT の z に 1e-12 以内 (磁場に比例する横の経路の分は 2 次の小ささ)。k64 は磁場を設定できない"""
    src = open(os.path.join(SCRIPTS, 'v250_und_dv__und_10mev.sad'), encoding='utf-8').read().replace('PHI=30', 'PHI=30 FBX=1e-9').replace('SLICE=4', 'SLICE=400')
    compat.set_mode('v2')
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    zf = parse(buf.getvalue())[4]
    zd = run('drift', '10mev', 'v2')[4]
    assert np.abs(zf - zd).max() < 1e-12
