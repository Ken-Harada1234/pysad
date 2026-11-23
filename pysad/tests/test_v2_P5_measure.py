# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 1: MEASURE・TRCOD を gs と比べる。

gs (K64、2026-11-19 の 10-19 ビルド) は、以前の記録 (B42) と違い、MEASURE を先頭の要素 (入口の MARK) の位置では出力する。ほかの位置 (既定の終端も) は
粒子の追跡 (tturn) で全粒子が落ちて np <= 1 になり、何も出力しない (FFS["MEASURE"] は {1, 0, 0})。TRCOD も先頭の位置だけ動き、ほかは SIGSEGV。
k64 モードは gs と同じにした (measure.py の `_k64_measure`・`_k64_trcod`): 初期分布 (ttinit、B41 の DPY の読み違いも)・乱数の消費 (MEASURE ごとに種を戻さない)・
統計の表示が、小さな FODO、大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)、PF リング (2.5 GeV、187 m)、ASTRID2 (0.58 GeV) で gs の表示と 4 桁で一致する。
gs の出力は example/02_optics_cod/fn_MEASURE_*.gs.txt。v2 は追跡を通せる (B42 の修正) ので、位置が先頭でなくても出力する (その値の検証は test_v2_C*_measure.py)。
"""
import re

import numpy as np
import pytest

from . import _p5_util as U

CASES = [('fodo', None), ('ring_a', 'A'), ('pf', 'PF'), ('astrid2', 'AST'), ('default_np', None)]


@pytest.mark.parametrize('name,lat', CASES)
def test_measure_and_trcod_statistics_match_gs_k64(name, lat):
    """MEASURE ORG・GAUSS・DP0・NOTRGAUSS・NP を変えた繰り返し・TRCOD ORG の統計ブロックが gs と一致 (gs の乱数も同じ順に消費される)"""
    src = U.example_src('measure_' + name, lat)
    got = U.stat_blocks(U.run(src, 'k64'))
    ref = U.stat_blocks(U.gs_text('measure_' + name))
    assert len(got) == len(ref) and len(ref) >= (1 if name == 'default_np' else 3)
    for i, (a, b) in enumerate(zip(got, ref)):
        assert U.blocks_close(a, b), (name, i, a[:5], b[:5])


def test_k64_py_centroid_is_the_b41_dp_jitter_not_zero():
    """gs の MEASURE の py の重心は設計値 0 でなく DP·gauss 程度 (4e-4: ttinit が軌道の揺らぎを DP から読む B41。x・px・y は 1e-5 以下)。k64 は gs と同じ"""
    ref = U.stat_blocks(U.gs_text('measure_fodo'))[0]
    com = [l for l in ref if l.startswith('   C of M')][0]
    f = U.floats(com.split(':', 1)[1])
    assert abs(f[3]) > 1e-4 and max(abs(f[0]), abs(f[1]), abs(f[2])) < 2e-5
    got = U.stat_blocks(U.run(U.example_src('measure_fodo'), 'k64'))[0]
    f2 = U.floats([l for l in got if l.startswith('   C of M')][0].split(':', 1)[1])
    assert np.allclose(f, f2, rtol=2e-3, atol=1e-19)


def test_k64_measure_at_other_positions_prints_nothing_like_gs():
    """終端 (`$$$`) など先頭以外の位置: gs は無出力で、FFS["MEASURE"] は np = 1 と 0 の重心・分散行列。k64 も同じ"""
    src = U.example_src('measure_fodo')
    gs = U.gs_text('measure_fodo')
    assert re.search(r'^R2 \{1,\{0,0,0,0,0,0\},\{\{0,0,0,0,0,0\}', gs, re.M)
    out = U.run(src, 'k64')
    assert re.search(r'^R2 \{1,\{0,0,0,0,0,0\},\{\{0,0,0,0,0,0\}', out, re.M)
    # `MEASURE $$$` の位置の統計ブロックは出ない: ブロック数は gs と同じ
    assert len(U.stat_blocks(out)) == len(U.stat_blocks(gs))


def test_k64_measure_consumes_the_random_seed_in_every_call_like_gs():
    """gs は FIXSEED でも MEASURE ごとに乱数の種を戻さない (NISTACK$FIXSEED@Push/Pop が効かない)。同じ位置の MEASURE を続けて呼ぶと分布が変わる (k64 も)。
    位置が先頭でない MEASURE (統計は出ない) も種を進める"""
    src = '''MOMENTUM = 3GEV;
DRIFT DS = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ORIGIN = ();
LINE RING = (ORIGIN 4*(QF DS QD DS));
FFS;
EMITX=1e-8; EMITY=1e-9; SIGZ=1e-3; DP=1e-3;
USE RING; CELL; CALC;
NP=300;
MEASURE ORIGIN;
MEASURE ORIGIN;
Exit;
'''
    b = U.stat_blocks(U.run(src, 'k64'))
    assert len(b) == 2 and b[0] != b[1]
    b2 = U.stat_blocks(U.run(src, 'v2'))
    assert len(b2) == 2
    # v2 は種を戻す (K64 の設計の意図): 2 回目も同じ分布
    assert b2[0] == b2[1]


def test_v2_initial_distribution_is_centered_on_the_design_orbit():
    """v2 (B41 の修正): py の重心は DP の揺らぎでなく設計値 (0) のまわりの統計誤差 (σpy/√N の数倍以内)"""
    out = U.run(U.example_src('measure_fodo'), 'v2')
    blocks = U.stat_blocks(out)
    for b in blocks[:3]:
        f = U.floats([l for l in b if l.startswith('   C of M')][0].split(':', 1)[1])
        sigpy = np.sqrt(U.floats([l for l in b if l.startswith('    py/p0')][0].split(':', 1)[1])[-1])
        assert abs(f[3]) < 6 * sigpy / np.sqrt(500)


def test_k64_trcod_beyond_first_element_is_refused_like_gs_sigsegv():
    """gs は先頭以外の位置の TRCOD が SIGSEGV (B42)。k64 は落とさず、理由を書いて何もしない (Twiss を変えない)"""
    src = '''MOMENTUM = 3GEV;
DRIFT DS = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ORIGIN = () M1 = ();
LINE RING = (ORIGIN 4*(QF DS QD DS) M1);
FFS;
USE RING; CELL; CALC;
TRCOD M1;
Print["DY ", Twiss["DPY", 1]];
Exit;
'''
    out = U.run(src, 'k64')
    assert 'SIGSEGV' in out
    dpy = U.floats([l for l in out.splitlines() if l.startswith('DY ')][0])[0]
    assert abs(dpy) < 1e-12                                  # TRCOD は Twiss の DPY を書かない
