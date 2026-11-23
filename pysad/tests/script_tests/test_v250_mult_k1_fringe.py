# -*- coding: utf-8 -*-
"""v2.5.0: K1 ≠ 0・ANGLE ≠ 0 の MULT で、FB1・FB2 の端 (線形のフリンジ) が MULT の途中にも入る (docs/30 B340)。

K64 の `tmulta.f` は、MULT を「半スライス (入口) + 中間スライス + 半スライス (出口)」の `tbend` の列で作り、前半には (fb1, fb2, mfr = -1)、
後半には (fb1, fb2, mfr = -2) を渡す。K1 = 0 の `tbend` は mfr で fb を選ぶが、K1 ≠ 0 は `tbendi` に行き (`tbend.f:595`)、`tbendi.f:73,159` は
mfr に関係なく fb1 を入口に・fb2 を出口に入れる。結果、MULT の前半の終わり (FB2) と後半の始め (FB1) にも端のキックが入る。
BEND (F1) は 1 本の `tbendi` で、端は両端の 1 回ずつ。BEND は FRINGE = 1・2・3 のいずれでも F1 を両端に、-1 (入口のみ)・-2 (出口のみ) で片側に使う
(`tffs.f:1263-1277`)。MULT の FRINGE = 1・2・3 は入口のみ・出口のみ・両端 (`tffs.f:1346-1364`) で、仕様が違う (誤りではない)。
同じ磁場を比べるには BEND の FRINGE = -1・-2・3 と MULT の FRINGE = 1・2・3 を対にする。

gs の保存した出力 (scripts/v250_mult_k1_fringe.{sad,gs.txt}: 3 GeV/c、L = 1 m、ANGLE = 0.3、K1 = 0.3、E1 = 0.5、E2 = 0.3、F1 (FB1 = FB2) = 0.5、
NORAD・NORFSW、TrackParticles の 7 粒子 (x・px・y・py の 1 つずつと δ = 1e-2 を含む複数の同時励起)):
- F1 = 0 では BEND と MULT は 2.2e-11 で一致 (K1 の本体と MULT の K1 は同じ)。
- F1 ≠ 0 では、両端 (BEND 3・MULT 3) で最大 1.41e-4 (y・py)、入口のみ・出口のみで各 7.0e-5 (途中の余計な端 1 つぶん) 違う。
v2 は前半の `tbend` に fb2、後半に fb1 を渡さず、MULT が BEND と 1e-13 以内で一致する。k64 は gs と同じ。
限界: ANGLE ≠ 0・K1 ≠ 0・NORAD の追跡だけ (RAD ON の `tbendi` の first / last の扱い、CALC `tmultae`→`tbende` は未確認)。K1 = 0 の MULT は K64 でも正しい。
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
SRC = open(os.path.join(SCRIPTS, 'v250_mult_k1_fringe.sad'), encoding='utf-8').read()
GS = open(os.path.join(SCRIPTS, 'v250_mult_k1_fringe.gs.txt'), encoding='utf-8', errors='replace').read()


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def parse(txt):
    res = {}
    cur = None
    txt = re.sub(r'\n {2,}(?=[ \d.\-])', ' ', txt)
    for l in txt.splitlines():
        m = re.match(r'^T (\S+)', l)
        if m:
            cur = m.group(1)
            res[cur] = {}
            continue
        m = re.match(r'^([XPYQZD]) \{(.*)\}', l)
        if m and cur:
            res[cur][m.group(1)] = np.array([float(v.replace('D', 'E')) for v in m.group(2).split(',')])
    return res


def run(mode):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(SRC, out=buf)
    return parse(buf.getvalue())


def md(a, b):
    return max(float(np.abs(a[c] - b[c]).max()) for c in 'XPYQZD')


G = parse(GS)
NAMES = ('both', 'entr', 'exit', 'nof')


def test_gs_mult_differs_from_bend_only_with_f1_and_k1():
    """gs: F1 = 0 は BEND = MULT (2.2e-11)。F1 ≠ 0: 両端 1.41e-4、入口のみ・出口のみ各 7.0e-5 (途中の余計な端 1 つぶん)"""
    assert md(G['Bnof'], G['Mnof']) < 1e-10
    assert md(G['Bboth'], G['Mboth']) == pytest.approx(1.41e-4, rel=0.01)
    assert md(G['Bentr'], G['Mentr']) == pytest.approx(7.0e-5, rel=0.01)
    assert md(G['Bexit'], G['Mexit']) == pytest.approx(7.0e-5, rel=0.01)


@pytest.mark.parametrize('n', NAMES)
def test_k64_matches_saved_gs(n):
    """k64: BEND・MULT とも gs と 1e-15 以内 (表示の 16 桁)"""
    p = run('k64')
    assert md(p['B' + n], G['B' + n]) < 1e-15
    assert md(p['M' + n], G['M' + n]) < 1e-15


@pytest.mark.parametrize('n', NAMES)
def test_v2_mult_equals_bend(n):
    """v2: MULT が BEND と 1e-12 以内で一致 (k64 / gs は F1 ≠ 0 で 7e-5〜1.4e-4 違う)。v2 の BEND は gs と 1e-6 以内 (別の v2 修正 B1 などの分)"""
    p = run('v2')
    assert md(p['B' + n], p['M' + n]) < 1e-12
    assert md(p['B' + n], G['B' + n]) < 1e-6
