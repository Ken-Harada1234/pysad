# -*- coding: utf-8 -*-
"""v2.5.0: BEND の CHI2 (tchge の tbrote の dchi2) の翻訳の抜けの補い (翻訳の誤りの修正。K64 の誤りではない)。

K64 の ANGLE ≠ 0 の BEND は CHI2 (BEND の中点の y 軸まわりのチルト) を tchge → tbrote (`tbende.f`・`tbend.f`) に渡す。pysad は BEND の CHI2 を
パーサが読まず (MULT の CHI2 だけ読む)、また ANGLE ≠ 0 の MULT の CHI2 は CALC・Emittance[] で NotImplementedError、tchge のスピンは dchi2 = 0 固定だった。
BEND の CHI2 は追跡・CALC・Emittance[]・1 周のスピン行列 (POL) のどれにも効かず、gs と食い違っていた (追跡で 4.5e-4、スピン行列で 2.8e-3)。
修正後、k64 が gs と一致する (追跡 1e-15、1 周のスピン行列 7 桁、Emittance[] の全量 5e-13、Twiss の最後の ν 6e-14)。

gs の保存した出力 (scripts/v250_bend_chi2__*): 3 GeV/c。
- track: 1 要素 (L = 1 m、ANGLE = 0.3、E1 = 0.5、E2 = 0.3) の TrackParticles の 7 粒子。BEND と MULT の CHI2 = 0.03 (と DROTATE・DX つき)。
- spin / spin_base: 24 セルの FODO リング (周長 105.6 m、RF・RAD・POL) の最後の 1 セルの偏向磁石に CHI2 = 0.03 を与え、`EMIT` の「Spin precession matrix」。
  gs では CHI2 = ±0.03 で直交性の欠陥 max|M Mᵀ − 1| が 1.9e-3 (K64 の `tbrote` の行の逐次更新 (B56) が dchi2 ≠ 0 で大きく出る)。v2 は 7.5e-8。
- emit_*: 同じリング (POL なし) の `Emittance[]` と最後の ν。BEND の CHI2 と、ANGLE ≠ 0 の MULT の CHI2。
限界: CHI2 = 0.03 rad・3 GeV/c・1 台。ANGLE = 0 の BEND の CHI2 は K64 の tsteee・tquade・tmulte に渡らないので使わない (gs でも効かない想定。未確認)。
ANGLE = 0 の MULT の CHI1・CHI2・DZ (tsolrot) は追跡・CALC・スピンとも未実装のまま (gs で CHI2 = 0.03 の L = 1 m の MULT が 7.8e-4 動く)。ANGLE ≠ 0 の MULT の
スピン (POL の Emittance[]) は未対応のまま。ampnf (TPSA) は BEND の CHI2 を入れない。
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


def _read(tag, ext):
    return open(os.path.join(SCRIPTS, 'v250_bend_chi2__%s.%s' % (tag, ext)), encoding='utf-8', errors='replace').read()


def _run(tag, mode):
    compat.set_mode(mode)
    buf = io.StringIO()
    run_sad_string(_read(tag, 'sad'), out=buf)
    return buf.getvalue()


def spin_matrix(txt):
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


def track(txt):
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


def tdiff(a, b):
    return max(float(np.abs(a[c] - b[c]).max()) for c in 'XPYQZD')


def nums(txt, tag):
    out = []
    for l in txt.splitlines():
        if l.startswith(tag):
            out += [float(x) for x in re.findall(r'[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?', l[len(tag):])]
    return np.array(out)


# ---------------------------------------------------------------- 追跡
def test_gs_bend_chi2_moves_the_particles():
    """gs: BEND の CHI2 = 0.03 は 7 粒子の追跡を 4.5e-4 動かす (ANGLE ≠ 0 の MULT の CHI2 と同じ大きさ)"""
    g = track(_read('track', 'gs.txt'))
    assert tdiff(g['bc'], g['b0']) == pytest.approx(4.54e-4, rel=0.01)
    assert tdiff(g['mc'], g['m0']) == pytest.approx(4.54e-4, rel=0.01)


@pytest.mark.parametrize('mode', ['k64', 'v2'])
@pytest.mark.parametrize('tag', ['b0', 'bc', 'bcr', 'm0', 'mc'])
def test_tracking_matches_gs(mode, tag):
    """k64・v2 とも BEND・MULT (ANGLE ≠ 0) の CHI2 (DROTATE・DX つきも) の追跡が gs と 1e-15 以内 (v2 は 1e-9)"""
    g = track(_read('track', 'gs.txt'))
    p = track(_run('track', mode))
    assert tdiff(p[tag], g[tag]) < (1e-15 if mode == 'k64' else 1e-9)


# ---------------------------------------------------------------- スピン行列
def test_gs_spin_matrix_is_not_orthogonal_with_a_tilted_bend():
    """gs: CHI2 = 0.03 の偏向磁石のあるリングの 1 周のスピン行列の直交性の欠陥が 1.9e-3 (CHI2 なしは表示の桁の 1.3e-7)"""
    assert orth(spin_matrix(_read('spin_base', 'gs.txt'))) < 2e-7
    assert orth(spin_matrix(_read('spin', 'gs.txt'))) == pytest.approx(1.9e-3, rel=0.05)


def test_k64_spin_matrix_matches_gs():
    """k64: 1 周のスピン行列が gs の表示 (7 桁) と 1e-6 以内 (修正前は 2.8e-3 違った)"""
    assert np.abs(spin_matrix(_run('spin', 'k64')) - spin_matrix(_read('spin', 'gs.txt'))).max() < 1e-6


def test_v2_spin_matrix_is_orthogonal():
    """v2: CHI2 = 0.03 でも直交 (欠陥 1e-6 以下。k64・gs の 1/1000 以下)"""
    assert orth(spin_matrix(_run('spin', 'v2'))) < 1e-6


# ---------------------------------------------------------------- CALC・Emittance[]
@pytest.mark.parametrize('tag', ['emit_none', 'emit_bend', 'emit_mult'])
def test_k64_emittance_and_tunes_match_gs(tag):
    """k64: Emittance[] の 41 個の数と Twiss の最後の ν が gs と 1e-9 (相対) 以内 (BEND の CHI2・ANGLE ≠ 0 の MULT の CHI2。修正前は ν が 1.3e-4 違った)"""
    g, p = _read(tag, 'gs.txt'), _run(tag, 'k64')
    for t in ('EM ', 'NU '):
        a, b = nums(g, t), nums(p, t)
        assert len(a) == len(b) > 0
        assert np.abs(a - b).max() / np.abs(a).max() < 1e-9, t


def test_chi2_changes_the_optics_in_gs():
    """gs: CHI2 = 0.03 は最後の ν を (41.15139、29.20493) から (41.15141、29.20506) に動かす (BEND と MULT で同じ)"""
    n0, nb, nm = (nums(_read(t, 'gs.txt'), 'NU ') for t in ('emit_none', 'emit_bend', 'emit_mult'))
    assert np.abs(nb - n0).max() > 1e-4
    assert np.abs(nb - nm).max() < 1e-9
