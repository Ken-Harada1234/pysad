# -*- coding: utf-8 -*-
"""v2.5.0: ANGLE ≠ 0 の MULT の追跡で、FB1・FB2 による長さ L' の補正が二重に入る (docs/30 B331)。

K64 は L' = L − ((φ FB1)² + (φ FB2)²)/(48 L)·sin(w/2)/sin(φ/2) を `tffs.f:1368` で p_L_MULT に入れ (BEND は `tffs.f:1282` で 1 回だけ)、
`tmulta.f:82-90` が受け取った al (すでに L') からスライス長 aln = al/ndiv をもう一度同じ式で引く。
ヘルプの BEND の F1 の定義 (L') を 1 回とし、MULT (ANGLE, FB1, FB2, FRINGE) は同じ磁場・同じフリンジの BEND (ANGLE, F1) と同じ追跡になるべき。
gs の保存した出力 (scripts/v250_mult_angle_fb_double__{bend,mult,multc}.{sad,gs.txt}: 3 GeV/c、L = 1 m、ANGLE = 0.3、E1 = 0.5、E2 = 0.3、F1 (FB1 = FB2) = 0.5、FRINGE = 3、
TrackParticles の 5 粒子 (x, px, y, py, δ の 1 つずつ)): BEND と MULT は x で 2.8e-7、z で 2.5e-8 食い違い、MULT の L を二重の補正の分だけ先に長くした multc
(L = 1.00018814134845: L − c(L) − c(L − c(L)) が BEND の L' に等しい) は BEND と 全成分で 0 (ビット単位) の差で一致する。
v2 は tmulta の 2 回目の補正を入れない (MULT が BEND と全桁一致)。k64 は gs と同じ。
限界: FRINGE = 1・2 の BEND との差は仕様の違い (BEND は F1 を両端に使う。BEND の -1・-2 と MULT の 1・2 を対にすると B331 の分だけ。test_v250_mult_k1_fringe.py)、試験は FRINGE = 3 だけ。ANGLE ≠ 0 の MULT の追跡だけ (gs で確認したのは K1 = 0、K2 以上なし、NORAD)。CALC (`tmultae`、K64 は同じ式を持つ) は pysad が FB1/FB2 を未対応。
K1 ≠ 0 の BEND と MULT の y の食い違いは別の誤り (B340、test_v250_mult_k1_fringe.py)。
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
BEAM = ('{{0.0, 1e-3, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 1e-3, 0.0}, '
        '{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 1e-2}, {1.0, 1.0, 1.0, 1.0, 1.0}}')


@pytest.fixture(autouse=True)
def _restore_mode():
    old = compat._mode
    yield
    compat.set_mode(old)


def script(elem, beam=BEAM):
    return '''MOMENTUM = 3GEV;
%s
MARK ST = ();
LINE ONE = (ST EL);
FFS; NORFSW; NORAD; NOFLUC;
$FORM = "22.16";
beam = %s;
USE ONE; CELL;
r = TrackParticles[{1, beam}, 1, 1, 1];
Print["X ", r[[2,1]]]; Print["P ", r[[2,2]]]; Print["Y ", r[[2,3]]]; Print["Q ", r[[2,4]]]; Print["Z ", r[[2,5]]]; Print["D ", r[[2,6]]];
''' % (elem, beam)


def parse(txt):
    d = {}
    for m in re.finditer(r'^([XPYQZD]) \{(.*)\}', txt, re.M):
        d[m.group(1)] = np.array([float(v.replace('D', 'E')) for v in m.group(2).split(',')])
    return d


def run(elem, mode, beam=BEAM):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(script(elem, beam), out=buf)
    return parse(buf.getvalue())


def maxdiff(a, b):
    return max(float(np.abs(a[c] - b[c]).max()) for c in 'XPYQZD')


def gs(tag):
    return parse(open(os.path.join(SCRIPTS, 'v250_mult_angle_fb_double__%s.gs.txt' % tag), encoding='utf-8').read())


def elem_of(tag):
    return [l for l in open(os.path.join(SCRIPTS, 'v250_mult_angle_fb_double__%s.sad' % tag), encoding='utf-8').read().splitlines() if l.startswith(('BEND', 'MULT'))][0]


# ---------------------------------------------------------------- gs
def test_gs_mult_differs_from_bend_but_precompensated_mult_equals_bend_bitwise():
    """gs: MULT は BEND と x で 2.8e-7・z で 2.5e-8 違うが、L を二重の補正の分だけ先に伸ばした MULT は BEND と全成分で同じ (差 0)"""
    b, m, mc = gs('bend'), gs('mult'), gs('multc')
    assert maxdiff(b, m) > 2e-7
    assert abs(b['X'] - m['X']).max() == pytest.approx(2.78e-7, rel=0.01)
    assert abs(b['Z'] - m['Z']).max() == pytest.approx(2.46e-8, rel=0.01)
    assert maxdiff(b, mc) == 0.0


@pytest.mark.parametrize('tag', ['bend', 'mult', 'multc'])
def test_k64_matches_saved_gs(tag):
    """k64: 3 つの台本の 5 粒子の x, px, y, py, z, δ が gs と 1e-15 以内で一致 (表示の 16 桁)"""
    assert maxdiff(run(elem_of(tag), 'k64'), gs(tag)) < 1e-15


# ---------------------------------------------------------------- v2: MULT = BEND
CASES = [
    ('L=1 ANGLE=0.3 E1=0.5 E2=0.3 FB1=0.5 FB2=0.5 FRINGE=3', 'L=1 ANGLE=0.3 E1=0.5 E2=0.3 F1=0.5 FRINGE=3'),
    ('L=1 ANGLE=0.1 FB1=0.2 FB2=0.2 FRINGE=3', 'L=1 ANGLE=0.1 F1=0.2 FRINGE=3'),
    ('L=2 ANGLE=0.2 E1=0.4 E2=0.1 FB1=0.4 FB2=0.4 FRINGE=3', 'L=2 ANGLE=0.2 E1=0.4 E2=0.1 F1=0.4 FRINGE=3'),
    ('L=1 ANGLE=0.3 FB1=0.5 FB2=0.5 FRINGE=3 DISFRIN=1', 'L=1 ANGLE=0.3 F1=0.5 FRINGE=3 DISFRIN=1'),
    ('L=1 ANGLE=0.3 E1=0.5 E2=0.3 FB1=0.5 FB2=0.5 FRINGE=3 DX=1e-3 ROTATE=0.01', 'L=1 ANGLE=0.3 E1=0.5 E2=0.3 F1=0.5 FRINGE=3 DX=1e-3 ROTATE=0.01'),
]


@pytest.mark.parametrize('mult,bend', CASES)
def test_v2_mult_equals_bend(mult, bend):
    """v2: ANGLE ≠ 0 の MULT (FB1 = FB2 = BEND の F1) の 5 粒子の追跡が BEND と 1e-15 以内で一致 (E1・E2・FRINGE 3・DISFRIN・DX・ROTATE)。k64 は 1e-9 以上違う"""
    b = run('BEND EL=(%s);' % bend, 'v2')
    m = run('MULT EL=(%s);' % mult, 'v2')
    assert maxdiff(b, m) < 1e-15, maxdiff(b, m)
    k = run('MULT EL=(%s);' % mult, 'k64')
    assert maxdiff(b, k) > 1e-9


def test_v2_unequal_fb1_fb2_follows_the_bend_formula():
    """FB1 ≠ FB2: BEND の F1 は 1 つなので、MULT の FB1 ≠ FB2 は BEND と比べられない。L' の式 ((φ FB1)² + (φ FB2)²)/(48 L) は同じなので、
    v2 の FB1 = 0.3・FB2 = 0.4 の MULT は、k64 の MULT の L を二重の補正の分だけ伸ばしたものと一致する (k64 の伸ばした L は ratio = sin(w/2)/sin(φ/2) を含めて解く)"""
    from scipy.optimize import brentq
    import math
    phi, e1, e2, fb1, fb2, Lb = 0.3, 0.5, 0.3, 0.3, 0.4, 1.0
    ratio = math.sin(0.5 * (phi - e1 * phi - e2 * phi)) / math.sin(0.5 * phi)
    c = lambda L: ((phi * fb1)**2 + (phi * fb2)**2) / 48.0 / L * ratio
    Lm = brentq(lambda L: (L - c(L)) - c(L - c(L)) - (Lb - c(Lb)), 0.5, 1.5, xtol=1e-16)
    el = 'MULT EL=(L=%.15g ANGLE=0.3 E1=0.5 E2=0.3 FB1=0.3 FB2=0.4 FRINGE=3);'
    # k64 の伸ばした L (二重の補正で L' = Lb − c(Lb)) と、v2 の L = Lb (1 回の補正で L' = Lb − c(Lb))
    assert maxdiff(run(el % Lb, 'v2'), run(el % Lm, 'k64')) < 1e-12


def test_v2_no_fringe_length_correction_unchanged():
    """FB1 = FB2 = 0 (または FRINGE = 0 で FB が効かない) の MULT は v2 でも k64 と全桁同じ (補正が無い)。K2 のキックつきも同じ"""
    for e in ('MULT EL=(L=1 ANGLE=0.3 E1=0.5 E2=0.3);',
              'MULT EL=(L=1 ANGLE=0.3 FB1=0.5 FB2=0.5 FRINGE=0);',
              'MULT EL=(L=1 ANGLE=0.3 E1=0.5 FRINGE=3 K2=0.5);'):
        assert maxdiff(run(e, 'v2'), run(e, 'k64')) == 0.0, e


def test_v2_k2_kick_with_fb_changes_only_through_the_length():
    """K2 ≠ 0・FB ≠ 0 の MULT (v2) は、K2 = 0 の差 (v2 − k64) と同じ長さの補正の分だけ k64 と違う: 差が 1e-5 以下で、K2 = 0 の差の 5 倍以内"""
    a = 'MULT EL=(L=1 ANGLE=0.3 E1=0.5 E2=0.3 FB1=0.5 FB2=0.5 FRINGE=3 K2=0.5);'
    b = 'MULT EL=(L=1 ANGLE=0.3 E1=0.5 E2=0.3 FB1=0.5 FB2=0.5 FRINGE=3);'
    d2 = maxdiff(run(a, 'v2'), run(a, 'k64'))
    d0 = maxdiff(run(b, 'v2'), run(b, 'k64'))
    assert 1e-8 < d2 < 1e-5 and d2 < 5 * d0, (d2, d0)
