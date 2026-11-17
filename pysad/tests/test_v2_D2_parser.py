# -*- coding: utf-8 -*-
"""D2 (1): MAIN 部の要素定義の式 (`ANGLE=PI/8`、`K1=SIN(0.1)`、`L=10 CM + 1 MM`) を gs と同じ文法・同じ値で読む。
参照値は gs (K64、2026-10-19 ビルド) の `Element["ANGLE","B1"]` の出力。pysad 自身の不足の修正なので互換モードによらない (k64 でも同じ)。
"""
import io

import pytest

# (式, gs の値)。None は gs が構文エラー (LALR(1) parser: syntax error) になる形
GS_CASES = [
    ('PI/8', .39269908169872414), ('2*PI', 6.283185307179586), ('-PI/8', -.39269908169872414),
    ('PI/8 DEG', 22.5), ('45 DEG', .7853981633974483), ('1 MM', .001), ('2 CM + 1 MM', .021),
    ('SIN(0.1)', .09983341664682815), ('sin(0.1)', .09983341664682815), ('SQRT(2)/10', .1414213562373095),
    ('COS(PI/3)', .5000000000000001), ('TAN(PI/8)', .41421356237309503), ('ATAN(1)', .7853981633974483),
    ('EXP(-1)', .36787944117144233), ('LOG(100)', 2.0), ('LN(10)', 2.302585092994046),
    ('SINH(0.1)', .10016675001984403), ('TANH(0.1)', .09966799462495582), ('COSH(0.1)', 1.0050041680558035),
    ('ATAN2(1+1,2)', .7853981633974483), ('ATAN2(PI,2 DEG)', 1.5596856728972892),
    ('2^3', 8.0), ('2^3^2', 512.0), ('-2^2', -4.0), ('-2*3', -6.0), ('PI*PI', 9.869604401089358),
    ('10 DEG/2', .08726646259971647), ('10/2 DEG', 286.4788975654116), ('1/0.5', 2.0), ('.5', .5), ('1E-1', .1),
    ('SIN(10 DEG)', .17364817766693033), ('SIN(0.1)^2', .009966711079379185), ('-SIN(0.1)', -.09983341664682815),
    ('SIN (0.1)', .09983341664682815), ('2*(1+2)', 6.0), ('SIN(PI)^2+COS(PI)^2', 1.0),
]
GS_ERRORS = ['pi/8', 'Pi/8', 'SQRT 2', 'ABS(-0.1)', 'ASIN(0.5)', 'LOG10(100)', '2**3', '2^-1', '2*-1', '--1', '1.5 PI',
             '(1+2) MM', '(PI/8) DEG', 'PI DEG', 'ATAN2(1 2)', 'ATAN2(1,2,3)', 'SIN(1,2)', '2 SIN(0.1)', 'SIN(0.1) DEG']


def _run(expr, key='ANGLE'):
    from pysad.ffscmd import run_sad_string
    src = f'''MOMENTUM=3 GEV;
MARK M0=();
BEND B1=(L=1 {key}={expr});
LINE R=(M0 B1);
FFS;
USE R;
Print["VAL ", Element["{key}","B1"]];
Exit;
'''
    buf = io.StringIO()
    import contextlib
    with contextlib.redirect_stdout(buf):
        run_sad_string(src, out=buf)
    return buf.getvalue()


def _val(out):
    rows = [l for l in out.splitlines() if l.startswith('VAL ')]
    assert rows, out[-400:]
    return float(rows[-1][4:])


@pytest.mark.parametrize('expr, ref', GS_CASES)
def test_element_value_expression_matches_gs(expr, ref):
    v = _val(_run(expr))
    assert v == pytest.approx(ref, rel=1e-15, abs=0)


@pytest.mark.parametrize('expr', GS_ERRORS)
def test_element_value_expression_rejected_like_gs(expr):
    out = _run(expr)
    assert 'syntax error in MAIN expression' in out
    assert 'Missing beamline' in out         # 要素が定義されないので LINE が作れない (gs は構文エラーで止まる)


def test_division_by_zero_and_domain_follow_gs():
    # gs: 0 での割り算は Infinity・NaN、定義域の外は NaN、LN(0) は -Infinity (エラーにならない)
    assert _val(_run('1/0')) == float('inf')
    assert _val(_run('-1/0')) == float('-inf')
    assert _val(_run('0/0')) != _val(_run('0/0'))
    assert _val(_run('LN(0)')) == float('-inf')
    assert _val(_run('SQRT(-1)')) != _val(_run('SQRT(-1)'))
    assert _val(_run('2^1000000')) == float('inf')


def test_main_variables_and_functions_in_element_values():
    from pysad.ffscmd import run_sad_string
    src = '''MOMENTUM=3 GEV;
A=PI/8; S=SINH(.1); C=1/0; D=ATAN2(1,2); F=A*2; G=10 DEG;
MARK M0=();
BEND B1=(L=1 ANGLE=A);
BEND B2=(L=1 ANGLE=S);
BEND B3=(L=1 ANGLE=C);
BEND B4=(L=1 ANGLE=D);
BEND B6=(L=1 ANGLE=F+A);
BEND B7=(L=1 ANGLE=G K1=2*PI/10 E1=-PI/100);
LINE R=(M0 B1 B2 B3 B4 B6 B7);
FFS;
USE R;
Print["VAL ", Element["ANGLE","B*"]];
Print["K1 ", Element["K1","B7"]];
Exit;
'''
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    out = buf.getvalue()
    # gs の出力 (A〜G の ANGLE)
    assert '{.39269908169872414,.10016675001984403,Infinity,.4636476090008061,1.1780972450961724,.17453292519943295}' in out
    assert 'K1 .6283185307179586' in out
