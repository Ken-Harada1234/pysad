# -*- coding: utf-8 -*-
"""PF の履歴ラティス 4 本 (example/pf/lattice/pf_*.sad): cell calc emit が走り、周長・チューン・エミッタンスが記録した値と一致する。
gs の保存値 (example/pf/optics/FFS_cell_calc_emit_PF_*.gs.txt) とも k64 で一致する。
`-RQUAD` の中に `-RQUADM` がある二重の反転の展開 (pysad/parser.py) の試験を含む (1981・1986 のラティスが使う)。"""
import io
import math
import os
import re

import pytest

from pysad import compat, ffscmd

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), 'example', 'pf')
# (名前, 周長 m, νx, νy, εx k64 (gs と同じ), εx v2)
CASES = [
    ('1981_460nmrad', 187.0741280, 5.2986816, 4.1637558, 4.033741e-7, 4.5981e-7),
    ('1986_128nmrad', 187.0741280, 8.4540000, 3.2999997, 1.213028e-7, 1.2817e-7),
    ('1997_036nmrad', 187.0741284, 9.5971735, 4.2902488, 3.472886e-8, 3.5585e-8),
    ('2005_036nmrad', 187.0354880, 9.5975305, 5.2789602, 3.453798e-8, 3.5374e-8),
]


def _run(case):
    src = open(os.path.join(EX, 'optics', 'FFS_cell_calc_emit_PF_%s.sad' % case), encoding='utf-8').read()
    buf = io.StringIO()
    ffscmd.run_sad_string(src, out=buf)
    return buf.getvalue()


def _val(txt, pat):
    return [float(x) for x in re.search(pat, txt).groups()]


@pytest.mark.parametrize('case,circ,nux,nuy,ex,ex_v2', CASES)
def test_cell_calc_emit_k64(case, circ, nux, nuy, ex, ex_v2):
    compat.set_mode('k64')
    txt = _run(case)
    assert _val(txt, r'circumference \[m\] =\s*(\S+)')[0] == pytest.approx(circ, abs=1e-6)
    tx, ty = _val(txt, r'tune nx ny =\s*(\S+)\s+(\S+)')
    assert tx == pytest.approx(nux, abs=1e-6) and ty == pytest.approx(nuy, abs=1e-6)
    assert _val(txt, r'EMITX \[m rad\] =\s*(\S+)')[0] == pytest.approx(ex, rel=1e-5)
    # 偏向磁石 28 台、角度の和 = 2 pi (近似)
    nb, sa = _val(txt, r'bends =\s*(\S+)\s+sum of ANGLE \[rad\] =\s*(\S+)')
    assert nb == 28 and sa == pytest.approx(2 * math.pi, abs=1e-4)   # 1981〜1997 は 6.283115 (2π より 7e-5 小さい)


@pytest.mark.parametrize('case,circ,nux,nuy,ex,ex_v2', CASES)
def test_matches_gs_saved_output(case, circ, nux, nuy, ex, ex_v2):
    """保存した gs の画面表示と、Print の結果の行 (周長・チューン・EMITX・クロマティシティ) が表示桁で一致"""
    compat.set_mode('k64')
    txt = _run(case)
    gs = open(os.path.join(EX, 'optics', 'FFS_cell_calc_emit_PF_%s.gs.txt' % case), encoding='utf-8', errors='replace').read()
    num = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?'
    for pat in (r'circumference \[m\] =\s*(%s)' % num, r'tune nx ny =\s*(%s)\s+(%s)' % (num, num),
                r'EMITX \[m rad\] =\s*(%s)' % num, r'with the design sextupoles: \{\s*(%s),\s*(%s)\}' % (num, num),
                r'natural \(SF\* SD\* = 0\):\s*\{\s*(%s),\s*(%s)\}' % (num, num)):
        a = _val(gs, pat); b = _val(txt, pat)
        for x, y in zip(a, b):
            assert y == pytest.approx(x, rel=1e-8, abs=1e-12), pat


@pytest.mark.parametrize('case,circ,nux,nuy,ex,ex_v2', CASES)
def test_emit_v2(case, circ, nux, nuy, ex, ex_v2):
    """v2 モード: 光学は k64 と同じで、εx は Jx の端の項 (B75) のぶん違う (名前の 460・128・36 nm rad に近い)。v2.4.0 から EMITDIV の既定が 4 (εx は 1 の値より 8e-4〜2.1e-3 小さい)"""
    compat.set_mode('v2')
    try:
        txt = _run(case)
    finally:
        compat.set_mode('k64')
    tx, ty = _val(txt, r'tune nx ny =\s*(\S+)\s+(\S+)')
    assert tx == pytest.approx(nux, abs=1e-6) and ty == pytest.approx(nuy, abs=1e-6)
    assert _val(txt, r'EMITX \[m rad\] =\s*(\S+)')[0] == pytest.approx(ex_v2, rel=2e-4)


def test_double_reversal_of_line():
    """`-A` の中の `-B` は B を正順に戻す (gs と同じ)。A = (X -B Y)、B = (P Q) のとき、-A = (-Y B -X) = Y' P Q X' で順は Y P Q X"""
    compat.set_mode('k64')
    src = '''MOMENTUM = 1 GEV;
DRIFT P=(L=1) Q=(L=2) X=(L=4) Y=(L=8);
MARK M0=();
LINE B=(P Q) A=(X -B Y) R=(M0 -A);
FFS;
USE R;
Print["S ", LINE["S","$$$"], " ", LINE["S", LINE["POSITION","Q"]], " ", LINE["S", LINE["POSITION","P"]]];
Exit;
'''
    buf = io.StringIO(); ffscmd.run_sad_string(src, out=buf)
    # R = M0 + reverse(A) = M0 Y P Q X : P の入口は 8、Q の入口は 9 (逆転が正しければ P が Q より前)
    s = [float(x) for x in re.search(r'S\s+(\S+)\s+(\S+)\s+(\S+)', buf.getvalue()).groups()]
    assert s[0] == pytest.approx(15.0) and s[2] == pytest.approx(8.0) and s[1] == pytest.approx(9.0)
