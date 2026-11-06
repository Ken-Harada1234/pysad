# -*- coding: utf-8 -*-
"""item4 (docs/history/work_item4.md) の回帰テスト: SADScript のリスト関数、tsvd (SAD の SVD ソルバ)、制御構文、CalculateOptics の骨格。
SAD との数値照合は bench_level24*.sad (run_regression.py 24a..24f) で行う。ここは SAD なしで走る小さなテスト。"""
import io
import os

import numpy as np
import pytest

from pysad.sadscript import Interp
from pysad.tsvd import tsvd


def run(src):
    o = io.StringIO()
    I = Interp(out=o)
    I.run(src, 'script')
    return o.getvalue()


def test_list_functions_match_sad_output():
    # SAD (bench_level24a) の出力: $FORM 既定では 1.0 は "1"
    out = run('Print[Position[{1,2,3,2},2]]; Print[Cases[{1,a,2,b},_Real]]; Print[Take[{1,2,3,4,5},-2]];'
              'Print[Insert[{1,2,3},x,-1]]; Print[Delete[{{1,2},{3,4}},{2,1}]]; Print[Complement[{5,3,1,3},{9}]];'
              'Print[Partition[{1,2,3,4,5},3,1]]; Print[Apply[f,{{1,2},{3}},{1}]]')
    lines = out.strip().split('\n')
    assert lines[0] == '{{2},{4}}'
    assert lines[1] == '{1,2}'
    assert lines[2] == '{4,5}'
    assert lines[3] == '{1,2,x,3}'
    assert lines[4] == '{{1,2},{4}}'
    assert lines[5] == '{1,3,5}'
    assert lines[6] == '{{1,2,3},{2,3,4},{3,4,5}}'
    assert lines[7] == '{f[1,2],f[3]}'


def test_implicit_multiplication_and_symbolic_print():
    assert run('Print[2 Pi]').strip().startswith('6.2831853')
    assert run('x=.; Print[x+y+x]').strip() == 'x+x+y'
    assert run('x=.; Print[2 x]').strip() == '(2 x)'


def test_tsvd_matches_numpy_pinv_for_well_conditioned():
    rng = np.random.default_rng(3)
    for n, m in [(5, 5), (4, 7), (9, 11), (8, 3)]:
        A = rng.normal(size=(n, m))
        b = rng.normal(size=n)
        assert np.allclose(tsvd(A, b, 1e-6), np.linalg.pinv(A) @ b, atol=1e-10)


def test_tsvd_reproduces_sad_step_for_ill_conditioned_ins_fit():
    """pf2 matching_v30_12 の 2 段目 (INS、9 条件 × 11 変数、特異値の比 1e-7〜1): SAD の 1 歩目の K1 変化を再現する"""
    d = np.load(os.path.join(os.path.dirname(__file__), 'data_tsvd_sad.npz'))
    x = tsvd(d['qu'], d['df'], 1e-6)
    assert np.allclose(x, d['dx_sad'], atol=2e-9)
    # numpy の硬い打ち切りでは別の解になる (この問題が 2 段目の発散の原因だった)
    u, s, vt = np.linalg.svd(d['qu'], full_matrices=False)
    inv = np.array([1 / v if v > 1e-6 * s[0] else 0 for v in s])
    xn = vt.T @ (inv * (u.T @ d['df']))
    assert not np.allclose(xn, d['dx_sad'], atol=1e-4)


def test_repeat_until_if_else():
    from pysad.ffscmd import run_sad_string
    o = io.StringIO()
    run_sad_string('MOMENTUM=1GEV;\nDRIFT D=(L=1);\nLINE L=(D D);\nFFS;\n'
                   'k=0;\nREPEAT 5\n k=k+1;\nUNTIL k>=3\nPrint["k=",k];\n'
                   'IF k==3\n Print["then"];\nELSE\n Print["else"];\nENDIF\n'
                   'IF k==4\n Print["A"];\nELSEIF k==3\n Print["B"];\nENDIF\nstop;\n', out=o)
    t = o.getvalue()
    assert 'k=3' in t and 'then' in t and 'else' not in t and 'B' in t.split('then')[1] and 'A' not in t.split('then')[1]


def test_flag_and_element_forms():
    from pysad.ffscmd import run_sad_string
    o = io.StringIO()
    run_sad_string('MOMENTUM=1GEV;\nMARK M=();\nQUAD Q=(L=0.3 K1=0.5);\nDRIFT D=(L=1);\nLINE L=(M Q D Q D);\nFFS;\n'
                   'USE L;\nPrint[Element["NAME"]]; Print[Element["K1","Q"]]; Print[Element["K1",{"Q"}]];'
                   'Print[Element["POSITION","Q"]]; Print[FLAG["CELL"]];\nstop;\n', out=o)
    lines = [x for x in o.getvalue().split('\n') if x.strip()]
    assert any(x.strip() == '{"M","Q","D"}' for x in lines)
    assert any(x.strip() == '.5' for x in lines)
    assert any(x.strip() == '{.5}' for x in lines)
