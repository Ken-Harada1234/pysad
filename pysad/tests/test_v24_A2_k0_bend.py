# -*- coding: utf-8 -*-
"""v2.4.0 A2: K0≠0 の偏向磁石 (挿入光源型) で ANGLE≠0・E1・E2≠0・K1≠0 がある場合の resdrive の解析式 (h・ADTS・ξ) を、独立な数値 (追跡の 1 周の写像のテイラー展開) と比べる。
リング: 小さな FODO の 1 セルに 5 台の偏向磁石のブロック (K0 と ANGLE の列 (1/2, -1, 1, -1, 1/2) 倍で閉軌道が局所のバンプ)。条件: k64、RF・放射なし。
ROTATE≠0 の偏向磁石は解析式に入れていない (警告が出る)。対応表は example/14_nonlinear_dynamics/normal_form/k0_bend_check.py の出力。"""
import os
import sys

import numpy as np
import pytest

from pysad import nform
from pysad import rdttrack as T
from pysad import resdrive as R

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'example', '14_nonlinear_dynamics', 'normal_form'))
import k0_bend_check as C  # noqa: E402


def _h_err(lat, rows):
    r = R.resonance_terms(lat, C.NAMES, rows=rows)
    rm = T.ring_map(lat.elements, C.P0, rows[0], amp=1e-5)
    h = T.rdt_from_map(rm)
    return max(abs(r.total[n] - h[n]) / abs(h[n]) for n in C.NAMES), r


@pytest.mark.parametrize('var,tol', [('ang0', 1e-8), ('ang', 5e-6), ('edge', 5e-6), ('edgeasym', 5e-6), ('k1', 5e-4)])
def test_k0_bend_block_h_matches_tracked_map(var, tol):
    """ANGLE=0 (MPW 型)・ANGLE≠0 かつ K0≠0・E1・E2≠0 (対称・非対称)・K1≠0 の K0 の偏向磁石の h (3 次の共鳴駆動項 5 種) が、追跡の 2 次のテイラー展開と tol 以内。警告なし"""
    lat = C.build(var)
    rows = lat.twiss_table()
    e, r = _h_err(lat, rows)
    assert e < tol
    assert not [w for w in r.warnings if 'K0' in w or '入れていない' in w]


def test_k0_bend_edge_adts_r6_matches_nform():
    """E1=E2≠0・ANGLE≠0・K0≠0 の磁石の ADTS (fourth='r6') が数値の正規形 (追跡の 3 次のテイラー展開) と 2e-4 以内 (K0=0 の FODO と同じ水準 1e-4)"""
    lat = C.build('edge')
    rows = lat.twiss_table()
    A = R.adts(lat, rows=rows, fourth='r6')
    rm3 = T.ring_map3(lat.elements, C.P0, rows[0], amp3=1e-4)
    an, _ = nform.adts_from_map(rm3)
    for k in ('xx', 'xy', 'yy'):
        assert abs(A.total[k] - an[k]) / abs(an[k]) < 2e-4


def test_natural_chromaticity_does_not_double_count_k0_bends():
    """K0≠0 の偏向磁石 (K1≠0) の自然な項は chromaticity_bend (要素の追跡から全体) が持つので chromaticity_natural は飛ばす。
    数え直す前は (自然 + 偏向磁石の分 + 6 極) − CALC が (−0.104, +0.108) だった → 0.02 以内"""
    lat = C.build('k1')
    rows = lat.twiss_table()
    ref = R.chromaticity_numeric(lat, order=1)['d1']
    nat = R.chromaticity_natural(lat, rows=rows)
    bend = R.chromaticity_bend(lat, rows=rows)
    sx = R.chromaticity_sextupole(lat, rows=rows)
    for i in (0, 1):
        assert abs(nat[i] + bend[i] + sx[i] - ref[i]) < 0.02


def test_rotated_k0_bend_is_included_in_coupled_path():
    """ROTATE≠0 の K0 の偏向磁石は、v2.4.0 G から結合した光学の模型 (ampnf) に入る (以前は解析式に入れず警告)。結合した経路 (`.coupled`) で、
    h は追跡の 2 次のテイラー展開と 1e-4 以内 (この ROTATE = 0.02 rad の小さな磁石では結合 3e-5)。警告は「結合した光学」(未対応でない)"""
    lat = C.build('rot')
    rows = lat.twiss_table()
    e, r = _h_err(lat, rows)
    assert getattr(r, 'coupled', False)
    assert not any('未対応' in w or '入れていない' in w for w in r.warnings)
    assert e < 1e-4
