# -*- coding: utf-8 -*-
"""K1≠0 かつ E1・E2≠0 の偏向磁石 (複合機能) の 3 次の閉形式 (bendnl.bend_analytic_h3: 本体 (h/2)x(px²+py²) + (K1 h/6)x³ の Gauss 積分 + 端の閉形式) と、
追跡の写像から取り出した原子の和 (bendnl.bend_atoms) の比較。以前は K1=0 の E≠0 と、K1≠0 の E=0 (ASTRID2) だけ確認していた (docs/25 §11.8)。
リング: FODO (3 GeV、105.6 m) の偏向磁石 B (L=1.2、ANGLE=0.2618) に K1・E1・E2 を足したもの (K1 の定義は 'log' = 既定)。条件: pysad の追跡、δ=0、閉軌道 0、nslice=16。
結果: v2 の追跡は閉形式と 1e-6 以内 (最大 6e-7。|K1| ≤ 0.2 [1/m²]、K1·L ≤ 0.24、E = -0.25〜0.4 rad。K1 に比例して増える差ではない)。
k64 の追跡は E≠0 かつ K1≠0 で 6e-4 (E2 だけ) 〜 6e-3 (E1) ずれる。k64 は sinp≠0 の端を薄い端 (K1 の磁場は半径方向の入口・出口の面から) で通し、閉形式は CALC (tbende) と
同じ厳密な端 (v2 の追跡) の式だから (track6.py の v2edge)。取り出しは追跡の写像から行うので、RDT・ADTS の値は k64 でも追跡と整合する。差を持つのは閉形式の比較だけ。"""
import os

import numpy as np
import pytest

from pysad import bendnl as B
from pysad import rdt2nd as Q
from pysad.parser import load_sad_string

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
FODO = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
BLANK = 'SF = (L=0.2 K2=0) SD = (L=0.2 K2=0)'
OLD = 'E1=0.1308996938995747 E2=0.1308996938995747'


def rel_diff(k1, e1, e2):
    s = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)').replace(OLD, 'K1=%r E1=%r E2=%r' % (k1, e1, e2))
    lat = load_sad_string(s).build_lattice('RING')
    rows = lat.twiss_table()
    i = [j for j, e in enumerate(lat.elements) if type(e).__name__ == 'Bend'][0]
    el = lat.elements[i]
    assert float(el.K1) != 0.0 and (e1 == 0.0 or float(el.E1) != 0.0)
    an = B.bend_analytic_h3(el, rows[i], rows[i + 1], Q.MONOS3)
    a, _ = B.bend_atoms(el, rows[i], rows[i + 1], Q.MONOS3, nslice=16)
    return np.abs(a.sum(axis=1) / 2 ** 1.5 - an).max() / np.abs(an).max()


@pytest.mark.parametrize('k1,e1,e2', [(0.1, 0.0, 0.0), (0.1, 0.3, 0.1), (0.1, -0.25, 0.4), (-0.1, 0.3, 0.1), (-0.2, -0.25, 0.4)])
def test_v2_combined_function_bend_closed_form_with_edges(v2_mode, k1, e1, e2):
    assert rel_diff(k1, e1, e2) < 1e-6


def test_k64_combined_function_bend_with_edges_differs_from_closed_form():
    """k64 の薄い端 (K1 の磁場は半径方向の面から) は、閉形式 (厳密な端) と E1 で 5e-3 台、E2 で 6e-4 だけ違う。E=0 では一致 (4e-7)"""
    assert rel_diff(0.1, 0.0, 0.0) < 1e-6
    assert 1e-3 < rel_diff(0.1, 0.3, 0.1) < 1e-2
    assert 1e-4 < rel_diff(0.1, 0.0, 0.1) < 2e-3


# --- 追跡 + NAFF + フィットの直接の ADTS との比較 (独立な基準)。FODO + 6 極 (K2=4・-6) の B に K1・E を足したもの。
# 精密版 adts(fourth='r6') と直接の ADTS (adtsdirect.direct_adts) の相対差: (K1, E1, E2) = (0.1, 0.3, 0.1) は v2 で −4e-5・−1e-5・−4e-5、k64 で −5e-5・−6e-6・−6e-5 (Jmax=1e-7、512 周)。
# K1=0.2・E=0 は νx ≈ −0.046 が整数に近く ADTS が 1.4e4 と大きいので、Jmax=1e-7 では追跡の J² の寄与で 6.9e-3 ずれる (解析の誤りではない)。Jmax を 3e-8・1e-8 にすると 5e-4・5e-5 に収束する。
def _adts_rel_diff(k1, e1, e2, jmax, nturn):
    from pysad import adtsdirect as D
    from pysad import resdrive as R
    s = FODO.replace(BLANK, 'SF = (L=0.2 K2=4) SD = (L=0.2 K2=-6)').replace(OLD, 'K1=%r E1=%r E2=%r' % (k1, e1, e2))
    lat = load_sad_string(s).build_lattice('RING')
    a = R.adts(lat, rows=lat.twiss_table(), fourth='r6').total
    t = D.direct_adts(lat, 3e9, jmax, order=2, nturn=nturn, nstep=3)
    return max(abs(a[k] - t[k]) / abs(t[k]) for k in ('xx', 'xy', 'yy'))


def test_v2_combined_function_bend_adts_matches_tracking(v2_mode):
    assert _adts_rel_diff(0.1, 0.3, 0.1, 1e-7, 512) < 2e-4
    assert _adts_rel_diff(0.2, 0.0, 0.0, 1e-8, 1024) < 2e-4


def test_k64_combined_function_bend_adts_matches_tracking():
    assert _adts_rel_diff(0.1, 0.3, 0.1, 1e-7, 512) < 2e-4
