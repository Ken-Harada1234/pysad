# -*- coding: utf-8 -*-
"""v2.5.0: B56 (`tbrote` の srot の行の逐次更新) の gs での確認: スピンの 1 周の回転行列は直交行列でなければならない。

K64 の `tbrote` (`tchge.f:197-202`) は rr (回転行列から単位行列を引いたもの) の srot への作用を、`do concurrent` の中で行を逐次更新する
(2 行目に更新後の 1 行目を使う) ので、(I + rr)·srot の行列積にならず、直交行列が直交行列に移らない。
gs の保存した出力 (scripts/v250_spin_tbrote__{base,droll}.{sad,gs.txt}: 3 GeV、24 セルの FODO リング (周長 105.6 m、矩形の偏向磁石)、
POL・RAD・RFSW、`EMIOUT; EMIT;` の「Spin precession matrix」(1 周のスピンの回転行列 srot(:,1:3)) の最後の 1 台の偏向磁石に DROTATE = 0.05 rad を与えたもの):
行列の直交性の欠陥 max|M Mᵀ − 1| は、ずれなし (base) で 1.26e-7 (表示の桁)、DROTATE = 0.05 で **9.76e-3** (gs)。pysad の k64 は gs と全桁一致、v2 は 8.4e-8。
限界: 1 つの要素・1 つの大きさ (DROTATE = 0.05 rad は大きめ。欠陥は回転角の 2 乗で効く)。DX・DY・DZ・ROTATE のずれは k64 が gs と一致 (test_v250_bend_chi2.py の記録)。
CHI2 は pysad が BEND の CHI2 を読んでいなかった (翻訳の抜け。直した。test_v250_bend_chi2.py)。
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
    rows = {}
    for l in txt.splitlines():
        m = re.match(r'\s+(sx|sy|sz)\s+(.*)$', l)
        if m:
            try:
                rows[m.group(1)] = [float(v.replace('D', 'E')) for v in m.group(2).split()]
            except ValueError:
                pass
    return np.array([rows['sx'][:3], rows['sy'][:3], rows['sz'][:3]])


def orth(M):
    return float(np.abs(M @ M.T - np.eye(3)).max())


def gs(tag):
    return parse(open(os.path.join(SCRIPTS, 'v250_spin_tbrote__%s.gs.txt' % tag), encoding='utf-8').read())


def run(tag, mode):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(open(os.path.join(SCRIPTS, 'v250_spin_tbrote__%s.sad' % tag), encoding='utf-8').read(), out=buf)
    return parse(buf.getvalue())


def test_gs_spin_matrix_is_not_orthogonal_with_a_rolled_bend():
    """gs: DROTATE = 0.05 の偏向磁石のあるリングの 1 周のスピン行列の直交性の欠陥が 9.8e-3 (ずれなしは表示の桁の 1.3e-7)"""
    assert orth(gs('base')) < 2e-7
    assert orth(gs('droll')) == pytest.approx(9.76e-3, rel=0.01)


@pytest.mark.parametrize('tag', ['base', 'droll'])
def test_k64_matches_saved_gs(tag):
    """k64 のスピン行列が gs の表示 (7 桁) と 1e-6 以内で一致"""
    assert np.abs(run(tag, 'k64') - gs(tag)).max() < 1e-6


def test_v2_spin_matrix_is_orthogonal():
    """v2: DROTATE = 0.05 でも直交 (欠陥 1e-6 以下。k64 の 1/1000 以下)"""
    assert orth(run('droll', 'v2')) < 1e-6
    assert orth(run('droll', 'v2')) < orth(run('droll', 'k64')) / 1000
