# -*- coding: utf-8 -*-
"""v2.5.0: CAVI の端 (tcavfrin) が、エネルギー変化 dh で p を変えるのに q = p_x/p を再スケールしない (docs/30 B332)。

K64 の本体 (`tcav.f` の do 20) と `ttcav` は p1r → p2r の変化で px = px·p1r/p2r とするが、`tcavfrin` はしない (正準運動量 p_x が保たれない)。
端付きの空洞 (FRINGE ≥ 0、DISFRIN = 0) の追跡写像が、エネルギー変化の x² の項 ∝ x·q·ω p0 v_f のぶん、シンプレクティックでない。
gs の保存した出力 (scripts/v250_cavi_edge_symplectic__{3gev,600kev}.{sad,gs.txt}: L = 0.1 m、VOLT = 1e5 V、500 MHz、PHI = 0.5、FRINGE = 1、
点 (x, px, y, py, z, δ) = (1e-3, 1e-2, −5e-4, 5e-3, −0.02, 0.01) と 12 点の中心差分 (h = 1e-6) の 6×6 ヤコビアン) の |JᵀSJ − S| は 3 GeV/c で 1.3e-8、0.6 MeV/c で 8.1e-5。
k64 はこれと 12 桁一致。v2 は q の再スケールを入れ、0.6 MeV/c で 5.8e-8、3 GeV/c で ~1e-10 (中心差分の丸め) に下がる。
限界: 中心差分のヤコビアンのシンプレクティック性だけ (端付きの空洞の写像の正しさそのものは B330 の数値積分)。
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
H = 1e-6
SM = np.zeros((6, 6))
for _i in (0, 2, 4):
    SM[_i, _i + 1] = 1
    SM[_i + 1, _i] = -1


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def parse(txt):
    body = txt.split('R {', 1)[1]
    nums = re.findall(r'\{([^{}]*)\}', body)
    return np.array([[float(v.replace('D', 'E')) for v in n.split(',')] for n in nums])


def jac_err(out):
    J = np.zeros((6, 6))
    for j in range(6):
        J[:, j] = (out[:6, 1 + 2 * j] - out[:6, 2 + 2 * j]) / (2 * H)
    return float(np.abs(J.T @ SM @ J - SM).max()), J


def run(tag, mode):
    compat.set_mode(mode)
    src = open(os.path.join(SCRIPTS, 'v250_cavi_edge_symplectic__%s.sad' % tag), encoding='utf-8').read()
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return parse(buf.getvalue())


def gs(tag):
    return parse(open(os.path.join(SCRIPTS, 'v250_cavi_edge_symplectic__%s.gs.txt' % tag), encoding='utf-8').read())


@pytest.mark.parametrize('tag', ['3gev', '600kev'])
def test_k64_matches_saved_gs_bitwise(tag):
    """k64 の 13 点の追跡が gs の保存した出力と 1e-14 以内で一致 (表示は 12 桁以上)"""
    assert np.abs(run(tag, 'k64') - gs(tag)).max() < 1e-14


def test_gs_edge_is_not_symplectic_at_low_energy():
    """gs: 端付きの空洞のヤコビアンの欠陥 |JᵀSJ − S| が 0.6 MeV/c で 8.1e-5、3 GeV/c で 1.3e-8 (K64 の欠陥)。"""
    e6, _ = jac_err(gs('600kev'))
    e3, _ = jac_err(gs('3gev'))
    assert e6 == pytest.approx(8.09e-5, rel=0.01)
    assert e3 == pytest.approx(1.30e-8, rel=0.02)


@pytest.mark.parametrize('tag,lim', [('3gev', 1e-9), ('600kev', 1e-6)])
def test_v2_edge_jacobian_is_symplectic(tag, lim):
    """v2: 同じ点の |JᵀSJ − S| が 0.6 MeV/c で 1e-6 以下 (5.8e-8)、3 GeV/c で 1e-9 以下 (中心差分の丸め)。k64 の 1/20 以下 (0.6 MeV/c は 1/1000 以下)"""
    ev, _ = jac_err(run(tag, 'v2'))
    ek, _ = jac_err(run(tag, 'k64'))
    assert ev < lim and ev < ek / 20, (ev, ek)


def test_v2_no_edge_unchanged_and_body_still_symplectic_scale():
    """DISFRIN = 1 (端なし) は v2 でも k64 と全桁同じ"""
    src = open(os.path.join(SCRIPTS, 'v250_cavi_edge_symplectic__600kev.sad'), encoding='utf-8').read().replace('FRINGE=1', 'DISFRIN=1')
    outs = []
    for m in ('k64', 'v2'):
        compat.set_mode(m)
        buf = io.StringIO()
        run_sad_string(src, out=buf)
        outs.append(parse(buf.getvalue()))
    assert np.array_equal(outs[0], outs[1])
