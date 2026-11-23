# -*- coding: utf-8 -*-
"""v2.5.0: COORD の DIR ≠ 0 の行列とヘルプの式の食い違い (physics_guide 11 章 付録 B。K64 の誤りではなく、ヘルプの積の順序の誤記と判断。コードは変えない)。

ヘルプ: DIR = 0 は (x, y, z)₁ = R₃R₂R₁ (x − DX, y − DY, z − DZ)₀、DIR ≠ 0 は R₃R₂⁽⁻⁾R₁ (x, y, z)₀ + (DX, −DY, 0) (R₂⁽⁻⁾ は s₂ の符号を逆にしたもの)。
K64 (`tcoord.f`) の DIR ≠ 0 は D M₀ᵀ D = R₁R₂⁽⁻⁾R₃ (D = diag(1, −1, −1)) と、出口で (DX, −DY, −DZ) を足す = D·(真の逆変換)·D。ヘルプとは積の順序だけが逆。
gs で、K64 の DIR ≠ 0 は DIR = 0 の真の逆変換を作れる: COORD(DX, DY, DZ, CHI1, CHI2, CHI3; DIR = 0) の次に
COORD(DX, −DY, −DZ, −CHI1, CHI2, −CHI3; DIR = 1) を置くと、TrackParticles が入力に戻る (差 7e-18)。同じパラメータでは戻らない (差 4e-2)。
ヘルプの積の順序 (R₃R₂⁽⁻⁾R₁) で真の逆変換を作るパラメータは (−0.0212, −0.0292, −0.0406) (CHI = 0.02, −0.03, 0.04 に対し) で、単純な符号の反転ではない。
つまり K64 の形が、(DX, −DY, −DZ) の符号の反転 (ヘルプの式にも (DX, −DY) がある) と整合する。ヘルプの記述の誤記の可能性が高い (確定ではない: 座標系の定義の問題で、物理量として決まるものではない)。
限界: 3 GeV/c・1 つのパラメータの組 (CHI は 0.02〜0.04 rad、DX・DY・DZ は mm)。ヘルプの式の誤記か、K64 の意図かは、設計者に確かめられていない。
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
BEAM = np.array([[1e-3, -2e-3, 0.0, 0.0, 5e-4, 0.0], [2e-3, 0.0, 1e-3, 0.0, -1e-3, 3e-3], [-1e-3, 0.0, 2e-3, 1e-3, 0.0, 0.0],
                 [0.0, 4e-3, 0.0, -2e-3, 1e-3, 0.0], [1e-3, 0.0, -1e-3, 0.0, -2e-3, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.01]])


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def parse(txt):
    body = txt.split('R {', 1)[1]
    return np.array([[float(v.replace('D', 'E')) for v in n.split(',')] for n in re.findall(r'\{([^{}]*)\}', body)])


def gs(k):
    return parse(open(os.path.join(SCRIPTS, 'v250_coord_dir__%s.gs.txt' % k), encoding='utf-8').read())


def run(k, mode):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(open(os.path.join(SCRIPTS, 'v250_coord_dir__%s.sad' % k), encoding='utf-8').read(), out=buf)
    return parse(buf.getvalue())


@pytest.mark.parametrize('k', ['dir0', 'dir1', 'pair', 'samepair'])
def test_k64_matches_saved_gs(k):
    """pysad の k64 の COORD (DIR = 0・1、対) の 6 粒子の追跡が gs と 1e-15 以内"""
    assert np.abs(run(k, 'k64') - gs(k)).max() < 1e-15


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_pysad_pair_returns_to_input_like_gs(mode):
    """COORD(DIR = 0) + COORD(DX, −DY, −DZ, −CHI1, CHI2, −CHI3; DIR = 1) は入力に戻る (gs・pysad の k64・v2 とも 1e-15 以内)。同じパラメータでは戻らない"""
    assert np.abs(gs('pair')[:6] - BEAM).max() < 1e-15
    assert np.abs(gs('samepair')[:6] - BEAM).max() > 1e-2
    assert np.abs(run('pair', mode)[:6] - BEAM).max() < 1e-15


def R1(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def R2(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, s], [0, -s, c]])


def R3(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def test_k64_dir_matrix_is_flipped_inverse_and_help_order_cannot_invert_simply():
    """K64 の DIR ≠ 0 の回転 = D M₀ᵀ D = R₁R₂⁽⁻⁾R₃ (ヘルプの R₃R₂⁽⁻⁾R₁ は積の順序だけ逆)。真の逆変換 M₀ᵀ を作るパラメータは K64 の形で (−CHI1, CHI2, −CHI3)、
    ヘルプの順序では (−0.0212, −0.0292, −0.0406) (CHI = 0.02, −0.03, 0.04)"""
    from scipy.optimize import least_squares
    a, b, c = 0.02, -0.03, 0.04
    M0 = R3(c) @ R2(b) @ R1(a)
    D = np.diag([1.0, -1.0, -1.0])
    K = R1(a) @ R2(-b) @ R3(c)
    assert np.abs(D @ M0.T @ D - K).max() < 1e-15
    assert np.abs(R1(-a) @ R2(-b) @ R3(-c) - M0.T).max() < 1e-15          # K64 の形に (−a, b, −c) を入れる = 真の逆
    f = lambda x: (R3(x[2]) @ R2(-x[1]) @ R1(x[0]) - M0.T).ravel()
    r = least_squares(f, [0, 0, 0], xtol=1e-15, ftol=1e-15)
    assert np.abs(r.fun).max() < 1e-13
    assert np.abs(r.x - np.array([-a, b, -c])).max() > 1e-3
    assert r.x == pytest.approx([-0.02119253, -0.02916992, -0.04060916], abs=2e-7)
