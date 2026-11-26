# -*- coding: utf-8 -*-
"""v2.5.0: SOL 領域内の QUAD (tsolqu・tsolqum) の ΔH の運動量の回転が、更新後の px を py に使う (docs/30 B333)。

K64 `tsolqu.f:122-123` (と `tsolqum` の `:302-303`) は
    pxi = a22*pxi + a24*py(i)
    pyi = -a24*pxi + a22*py(i)
と書き、2 行目に更新後の pxi を使う。正しい回転は py' = −a24·px_old + a22·py で、K64 の写像は行列式 cos²(wℓ) (w = bz·(−dpz)/(1+dpz)、
ℓ = スライスの長さ) で、運動量の面積が縮む (シンプレクティックでない)。
gs の保存した出力 (scripts/v250_sol_quad_rotation__{p05,p10}.{sad,gs.txt}: 3 GeV/c、SOL (BZ = 13 T、bz = 1.3 /m) + QUAD (L = 0.4 m、K1 = 0.36 /m²) + SOL、NORAD、
点 (x, px, y, py) = (1e-3, 0.05, −5e-4, 0.02) と (1e-3, 0.1, −5e-4, 0.05) の周りの 13 点 (±1e-6) の中心差分の 6×6 ヤコビアン) の 4×4 ブロック (x, px, y, py) の |JᵀSJ − S| は
1.54e-7 (px = 0.05) と 3.39e-6 (px = 0.1)。det J − 1 = −1.50e-7、−3.23e-6。k64 は gs と全桁一致。v2 は 4×4 の欠陥が 5e-10 (中心差分の丸め)、det J − 1 = −6e-11、−5e-10。
限界: (y, δ)・(x, δ) の欠陥 (8.4e-7・1.0e-8 と 7.7e-6・1.7e-6) は v2 でも残る (別の原因。gs でスライスを細かくしても減らない。不確定。docs/history/V250_sol_quad_xdelta_defect.md)。
"""
import io
import os

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
    import re
    body = txt.split('R {', 1)[1]
    nums = re.findall(r'\{([^{}]*)\}', body)
    return np.array([[float(v.replace('D', 'E')) for v in n.split(',')] for n in nums])


def jac(out):
    J = np.zeros((6, 6))
    for j in range(6):
        J[:, j] = (out[:6, 1 + 2 * j] - out[:6, 2 + 2 * j]) / (2 * H)
    return J


def defect(out):
    J = jac(out)
    D = J.T @ SM @ J - SM
    return float(np.abs(D[:4, :4]).max()), float(np.linalg.det(J) - 1.0)


def run(tag, mode):
    compat.set_mode(mode)
    src = open(os.path.join(SCRIPTS, 'v250_sol_quad_rotation__%s.sad' % tag), encoding='utf-8').read()
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    return parse(buf.getvalue())


def gs(tag):
    return parse(open(os.path.join(SCRIPTS, 'v250_sol_quad_rotation__%s.gs.txt' % tag), encoding='utf-8').read())


@pytest.mark.parametrize('tag', ['p05', 'p10'])
def test_k64_matches_saved_gs_bitwise(tag):
    """k64 の 13 点の追跡が gs の保存した出力と 1e-15 以内で一致 (7 行 × 13 列)"""
    assert np.abs(run(tag, 'k64') - gs(tag)).max() < 1e-15


@pytest.mark.parametrize('tag,d4,det', [('p05', 1.54e-7, -1.50e-7), ('p10', 3.39e-6, -3.23e-6)])
def test_gs_jacobian_is_not_symplectic(tag, d4, det):
    """gs: 4×4 ブロックの欠陥と det J − 1 が、px = 0.05 で 1.5e-7、0.1 で 3.4e-6 (K64 の欠陥。中心差分の丸め 1e-10 より 3 桁以上大きい)"""
    e, dt = defect(gs(tag))
    assert e == pytest.approx(d4, rel=0.01) and dt == pytest.approx(det, rel=0.01)


@pytest.mark.parametrize('tag', ['p05', 'p10'])
def test_v2_jacobian_4x4_is_symplectic(tag):
    """v2: 4×4 ブロックの欠陥が 2e-9 以下、|det J − 1| が 2e-9 以下 (k64 の 1/100 以下)"""
    ek, dk = defect(run(tag, 'k64'))
    ev, dv = defect(run(tag, 'v2'))
    assert ev < 2e-9 and abs(dv) < 2e-9, (ev, dv)
    assert ev < ek / 100 and abs(dv) < abs(dk) / 100


def test_k64_det_deficit_has_the_size_of_the_cos2_rotation_error():
    """k64 の det J − 1 の大きさが、運動量の回転の角 θ = wℓ (w = bz(−dpz)/(1+dpz)、4 回) の θ² の和 (p10 で 3.0e-6 の見積り) の 0.8〜1.2 倍 (大きさの確認)"""
    _, dk = defect(run('p10', 'k64'))
    assert 2.4e-6 < abs(dk) < 3.6e-6
