# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 2: DETR >= 1 (x-y のモードが入れ替わる) の CALC・INS・Emittance[] を gs と比べる。

入れ替わりは、リングの両端に 90 度の座標回転 (COORD CHI3 = ±π/2) を入れて作る。回転の間 (リングの本体) は x モードが鉛直、y モードが水平になり、
K64 の DETR = 1 + xyth − R1 R4 + R2 R3 = 1.375 (R = 0) になる。大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m)・PF リング (2.5 GeV、187 m)・ASTRID2 (0.58 GeV)。
結果 (gs は example/02_optics_cod/flag_DETR_*.gs.txt):
- DETR・R・β・Emittance[] のチューンと減衰率・エミッタンス・U0: k64 は gs と一致 (K64 の誤りは見つからなかった)。
- 一周の位相 NX・NY の**整数部**: 純粋な 90 度回転は x < 0 で y が丸め誤差になる atan2 を通り、K64 の位相の跳びは丸め誤差の符号で +π か −π になる (K64 の誤りの候補 P5-3)。
  gs は 3 つのリングとも跳びが打ち消し (整数部が本体の物理的なチューン)、pysad の k64 は丸め誤差が違い ±1 ずれる。小数部は 8 桁で一致。
  v2 は入れ替わりに入る跳びを +π、出る跳びを −π に決めて対で打ち消し、gs と同じ整数部になる。
- v2 の Emittance[] は、回転なしのリングの Emittance[] とモードの並びを除いて同じ値 (入れ替わりが物理量に影響しない)。
"""
import re

import numpy as np
import pytest

from . import _p5_util as U

RINGS = [('astrid2', 'AST'), ('pf', 'PF'), ('ring_a', 'A')]


def _tag(text, tag):
    for l in text.splitlines():
        if l.startswith(tag + ' '):
            return U.floats(l[len(tag):])
    raise AssertionError(tag + ' がない')


def _tagblock(text, tag):
    """複数行にまたがる Print の出力 ('tag {...}')"""
    i = text.index('\n' + tag + ' ') + 1
    seg = text[i + len(tag):]
    m = re.search(r'\n(?=\s*[A-Za-z?])', seg)
    return U.floats(seg[:m.start()] if m else seg)


@pytest.mark.parametrize('name,lat', RINGS)
def test_detr_rows_and_beta_match_gs_in_both_modes(name, lat, ):
    """本体の DETR は 1.375、R はほぼ 0。β・位置は gs と一致。CELL CALC と INS CALC の両方 (k64・v2)"""
    ref = U.gs_text('detr_' + name)
    for mode in ('k64', 'v2'):
        out = U.run(U.example_src('detr_' + name, lat), mode)
        for tag in ('DETR', 'B', 'INS DETR', 'INS B'):
            g, p = _tagblock(ref, tag), _tagblock(out, tag)
            assert len(g) == len(p)
            assert np.allclose(g, p, rtol=1e-7, atol=1e-12), (mode, tag, g, p)
        d = _tagblock(out, 'DETR')
        assert abs(d[2] - 1.375) < 1e-9 and abs(d[3] - 1.375) < 1e-9 and abs(d[4] - 1.375) < 1e-9 and abs(d[0]) < 1e-9


@pytest.mark.parametrize('name,lat', RINGS)
def test_phase_advance_fraction_matches_gs_k64_and_v2_integer_too(name, lat):
    ref = U.gs_text('detr_' + name)
    g = _tag(ref, 'N')
    for tag in ('N', 'INS N'):
        gref = _tag(ref, tag)
        pk = _tag(U.run(U.example_src('detr_' + name, lat), 'k64'), tag)
        assert np.allclose(np.mod(pk, 1.0), np.mod(gref, 1.0), atol=1e-7) or np.allclose(np.mod(pk - gref, 1.0), 0, atol=1e-7)
        assert max(abs(np.array(pk) - np.array(gref))) <= 1.0 + 1e-7             # k64: 整数部は丸め誤差の符号で最大 1 ずれる (gs と同じとは限らない)
    out2 = U.run(U.example_src('detr_' + name, lat), 'v2')
    for tag in ('N', 'INS N'):
        assert np.allclose(_tag(out2, tag), _tag(ref, tag), atol=1e-7), (name, tag)


@pytest.mark.parametrize('name,lat', RINGS)
def test_emittance_with_swapped_modes_matches_gs_k64(name, lat):
    ref = U.gs_text('detr_' + name)
    out = U.run(U.example_src('detr_' + name, lat), 'k64')
    for tag, tol in (('Tunes', 1e-7), ('Emittances', 1e-6), ('Damp', 1e-6), ('U0', 1e-9)):
        g, p = _tagblock(ref, tag), _tagblock(out, tag)
        assert np.allclose(g, p, rtol=tol, atol=1e-30), (name, tag, g, p)


@pytest.mark.parametrize('name,lat', [('astrid2', 'AST'), ('pf', 'PF')])
def test_v2_emittance_is_independent_of_the_mode_exchange(name, lat):
    """v2: 回転で x・y を入れ替えたリングの Emittance[] (エミッタンス・減衰率・U0) は、回転なしのリングと、モードの並びを除いて同じ"""
    swapped = U.run(U.example_src('detr_' + name, lat), 'v2')
    plain_src = open(U.need(U.LAT[lat]), encoding='utf-8', errors='replace').read() + '''
FFS;
$FORM="14.8";
USE RING; CELL; CALC;
e = Emittance[];
Print["Tunes ", Tunes/.e];
Print["Emittances ", Emittances/.e];
Print["Damp ", DampingRate/.e];
Print["U0 ", EnergyLossU0/.e];
Exit;
'''
    plain = U.run(plain_src, 'v2')
    es, ep = sorted(_tagblock(swapped, 'Emittances')), sorted(_tagblock(plain, 'Emittances'))
    assert np.allclose(es[1:], ep[1:], rtol=1e-6)                                  # 0 に近い mode (縦ではない) の並びを除く
    assert np.allclose(sorted(_tagblock(swapped, 'Damp')), sorted(_tagblock(plain, 'Damp')), rtol=1e-6)
    assert np.isclose(_tagblock(swapped, 'U0')[0], _tagblock(plain, 'U0')[0], rtol=1e-9)
