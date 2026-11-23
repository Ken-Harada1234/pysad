"""v2.4.0 担当 E: 組込の数学関数・特殊関数・リスト関数の監査 (docs/30 の B180〜B199)。

根拠は mpmath (高精度の直接計算) で作った参照値 `specfn_audit_ref.json` (式、実部、虚部。mpmath は環境に無いので参照値だけを持つ。
作り方: collective312 の venv の mpmath 1.4.1 で `mp.dps = 40`)。参照は主値の定義 (Mathematica と同じ切断線の取り方) で、gs の値ではない。
k64 は変えていないので、k64 の値 (gs の NaN・符号違い) が従来のままであることも固定する。
"""
import io
import json
import math
import os
import re

import pytest

_NUM = r'(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'
_REF = os.path.join(os.path.dirname(__file__), 'specfn_audit_ref.json')


def _run(src, mode):
    from pysad import compat
    from pysad.ffscmd import run_sad_string
    old = compat._mode
    compat.set_mode(mode)
    try:
        buf = io.StringIO()
        run_sad_string('FFS;\n$FORM="25.17";\nPageWidth=1999;\n' + src + '\nExit;\n', out=buf)
    finally:
        compat.set_mode(old)
    return buf.getvalue()


def _evals(exprs, mode='v2'):
    src = ''.join('Print["R%d= ",%s];\n' % (i, e) for i, e in enumerate(exprs))
    got = {}
    for line in _run(src, mode).splitlines():
        m = re.match(r'R(\d+)= (.*)', line)
        if m:
            got[int(m.group(1))] = m.group(2).strip()
    return [got.get(i) for i in range(len(exprs))]


def _parse(s):
    """出力の数 (複素は `( a+ b I)`・`( b I)`・`I`) を complex に"""
    if s is None:
        return None
    t = s.replace(' ', '')
    if t in ('NaN', 'Indeterminate'):
        return complex('nan')
    t = t.replace('-INF', '-1e999').replace('INF', '1e999').strip('()')
    if t.endswith('I'):
        body = t[:-1]
        m = re.fullmatch(r'(?:([-+]?%s)(?=[-+]))?([-+]?(?:%s)?)' % (_NUM, _NUM), body)
        if m is None:
            return None
        re_ = float(m.group(1)) if m.group(1) else 0.0
        im = m.group(2)
        return complex(re_, 1.0 if im in (None, '', '+') else (-1.0 if im == '-' else float(im)))
    try:
        return complex(float(t), 0.0)
    except ValueError:
        return None


def _s(e, mode='v2'):
    return _evals([e], mode)[0]


# ------------------------------------------------------------------------------ mpmath との比較
def test_v2_agrees_with_mpmath_reference():
    ref = json.load(open(_REF))
    got = _evals([r[0] for r in ref], 'v2')
    bad = []
    for (e, rr, ri), g in zip(ref, got):
        z = _parse(g)
        want = complex(rr, ri)
        if z is None or z != z or abs(z - want) > 2e-9 * abs(want) + 1e-300:
            bad.append((e, want, g))
    assert not bad, bad[:8]


# k64 は変えない: gs と同じ (誤った) 値のまま
@pytest.mark.parametrize('expr,expect', [
    ('ArcCosh[.5]', complex(0, -1.0471975511965979)),      # gs: -1.047 I (真値は +1.047 I。B180)
    ('ArcSin[2.]', complex(math.pi / 2, 1.3169578969248166)),
    ('ArcCos[2.]', complex(0, -1.3169578969248166)),
    ('Round[2.5]', 3.0), ('Round[-2.5]', -3.0), ('Round[0.5]', 1.0),
    ('Mod[7,-3]', -2.0),
])
def test_k64_unchanged_values(expr, expect):
    z = _parse(_s(expr, 'k64'))
    assert abs(z - complex(expect)) < 1e-14


@pytest.mark.parametrize('expr', ['ArcTanh[2.]', 'BesselY[0.,-1.5]', 'HermiteH[3,1.+1.*I]', 'LegendreP[2.5,1,.3]'])
def test_k64_gs_wrong_values_kept(expr):
    # gs (10-19 ビルド) が NaN を返すもの。LegendreP[2.5,1,.3] は gs が 0 (誤り)、pysad の k64 は NaN
    z = _parse(_s(expr, 'k64'))
    assert z is None or z != z or 'Hermite' in expr, (expr, z)
    assert 'wrongtype' in _run('Print[HermiteH[3,1.+1.*I]];', 'k64')


# gs (2026-11-22 の gs で実測) が真値を返すもの: k64 も v2 と同じ値 (従来の pysad の k64 は NaN・実部だけだった)
_GS_TRUE = [
    ('BesselJ[.5,-1.]', complex(0, 0.6713967071418031)), ('BesselI[.5,-1.]', complex(0, 0.9376748882454876)),
    ('BesselK[0.,-2.]', complex(0.11389387274953344, -7.161528439050255)),
    ('CoshIntegral[-2.5]', complex(3.5244254883541655, math.pi)),
    ('LogIntegral[-2.]', complex(0.03553227591356, 3.7351504825512546)),
    ('ExpIntegralE[1,-2.5]', complex(-7.0737658945786, -math.pi)),
    ('Hypergeometric2F1[1.,2.,3.,1.5]', complex(-0.7172025061689375, -2.792526803190927)),
    ('Gamma[0,1.]', 0.21938393439552029), ('Gamma[-.5,2.]', 0.030098757100186467), ('GammaRegularized[-.5,1.5]', -0.019522369872295782),
    ('ProductLog[-.5]', complex(-0.7940236323446894, 0.7701117505103791)),
    ('LerchPhi[-3.,2.,1.]', 0.646458473588903), ('PolyLog[-2.,2.]', -6.0),
    ('Hypergeometric1F1[.5,1.5,25.+3.*I]', complex(-1409411536.3294244, 380467084.496)),
    ('SinhIntegral[1.+1.*I]', complex(0.882453805007918, 1.104222658235582)),
]


@pytest.mark.parametrize('mode', ['k64', 'v2'])
def test_values_gs_computes_correctly_are_true_in_both_modes(mode):
    got = _evals([e for e, _ in _GS_TRUE], mode)
    for (e, want), g in zip(_GS_TRUE, got):
        z = _parse(g)
        assert z is not None and abs(z - complex(want)) < 2e-9 * abs(want) + 1e-12, (mode, e, g)
    assert _parse(_s('Binomial[-5.,-2.]', mode)).real == 0.0



# ------------------------------------------------------------------------------ 基本関数 (v2 と k64 の別・共通)
def test_arc_functions_branch_values():
    # 主値 (Mathematica): ArcSin[2.] = 1.5708 - 1.3170 I、ArcCos[2.] = 1.3170 I、ArcCosh[.5] = 1.0472 I、ArcTanh[2.] = .5493 - 1.5708 I
    z = _parse(_s('ArcSin[2.]'))
    assert abs(z - complex(math.pi / 2, -math.acosh(2.0))) < 1e-14
    z = _parse(_s('ArcSin[-2.]'))
    assert abs(z - complex(-math.pi / 2, math.acosh(2.0))) < 1e-14
    z = _parse(_s('ArcCos[2.]'))
    assert abs(z - complex(0, math.acosh(2.0))) < 1e-14
    z = _parse(_s('ArcCos[-2.]'))
    assert abs(z - complex(math.pi, -math.acosh(2.0))) < 1e-14
    z = _parse(_s('ArcCosh[.5]'))
    assert abs(z - 1j * math.acos(0.5)) < 1e-14
    z = _parse(_s('ArcTanh[2.]'))
    assert abs(z - complex(0.5 * math.log(3.0), -math.pi / 2)) < 1e-14
    # 恒等式 sin(ArcSin x) = x、cosh(ArcCosh x) = x が複素の分枝でも成り立つ
    for x in ('2.', '-3.5', '.3', '-.9'):
        r = _evals(['Sin[ArcSin[%s]]' % x, 'Cos[ArcCos[%s]]' % x, 'Cosh[ArcCosh[%s]]' % x, 'Tanh[ArcTanh[%s]]' % x if abs(float(x)) != 1 else '0'])
        for v in r[:3]:
            assert abs(_parse(v) - float(x)) < 1e-13


def test_both_modes_pysad_own_errors_fixed_like_gs():
    # gs が正しく、従来の pysad が誤っていたもの: k64 も v2 も gs に合わせる (B183)
    for mode in ('k64', 'v2'):
        r = _evals(['(-2.)^1501.', '(-2.)^1500.', 'Floor[Infinity]', 'Mod[Infinity,2]', 'Round[1.5+2.5*I]', 'Floor[1.5+2.5*I]', 'Sign[3.+4.*I]',
                    'Exp[1000.+1.*I]', 'Exp[1000.]'], mode)
        assert _parse(r[0]).real == -float('inf')
        assert _parse(r[1]).real == float('inf')
        assert _parse(r[2]).real == float('inf')
        assert _parse(r[3]) != _parse(r[3])
        assert _parse(r[4]) == (complex(2, 3) if mode == 'k64' else complex(2, 2))      # v2 は半整数を偶数へ (B184)
        assert _parse(r[5]) == complex(1, 2)
        assert abs(_parse(r[6]) - complex(.6, .8)) < 1e-15
        z = _parse(r[7])
        assert z.real == float('inf') and z.imag == float('inf')
        assert _parse(r[8]).real == float('inf')


def test_complex_with_zero_imag_becomes_real():
    # gs: (1+I)(1-I) = 2、I^2 = -1、Cos[I] は実数 (虚部が 0 の複素数は実数。B199)
    for mode in ('k64', 'v2'):
        r = _evals(['(1+I)*(1-I)', 'I^2', 'Cos[I]', 'I^I'], mode)
        assert all('I' not in t.replace('INF', '') for t in r), r
        assert _parse(r[0]) == 2.0 and _parse(r[1]) == -1.0


def test_sqrt_complex_exact():
    z = _parse(_s('Sqrt[-4.+0.*I]'))
    assert z == 2j or abs(z - 2j) == 0.0
    z = _parse(_s('Sqrt[3.+4.*I]'))
    assert abs(z - (2 + 1j)) < 1e-15
    # k64 は従来の複素数の冪 (実部に 1e-16 の雑音) のまま
    z = _parse(_s('Sqrt[-4.+0.*I]', 'k64'))
    assert abs(z - 2j) < 1e-14


def test_round_half_even_v2_and_mod():
    r = _evals(['Round[2.5]', 'Round[0.5]', 'Round[1.5]', 'Round[-2.5]', 'Round[-0.5]', 'Round[2.4]', 'Round[2.6]', 'Round[-2.6]'], 'v2')
    assert [_parse(t).real for t in r] == [2.0, 0.0, 2.0, -2.0, 0.0, 2.0, 3.0, -3.0]
    r = _evals(['Round[2.5]', 'Round[-2.5]', 'Round[0.5]'], 'k64')
    assert [_parse(t).real for t in r] == [3.0, -3.0, 1.0]
    # Mod は符号が除数と同じ (Mathematica)、大きな商でも正確 (v2)
    r = _evals(['Mod[-7,3]', 'Mod[7,-3]', 'Mod[5.5,2]', 'Mod[1.e17+1.,3.]', 'Mod[-1.e-20,1.]'], 'v2')
    assert [_parse(t).real for t in r[:3]] == [2.0, -2.0, 1.5]
    assert _parse(r[3]).real == (1e17 + 1.0) % 3.0
    assert _parse(r[4]).real == (-1e-20) % 1.0


# ------------------------------------------------------------------------------ 特殊関数の個別
def test_negative_x_bessel_and_identities():
    # J_ν(-x) = e^{iπν} J_ν(x) (主値)。独立な式
    from scipy import special as sp
    for nu in (.5, 1.5, 2.7):
        z = _parse(_s('BesselJ[%r,-3.]' % nu))
        assert abs(z - complex(math.cos(math.pi * nu), math.sin(math.pi * nu)) * sp.jv(nu, 3.0)) < 1e-13
        z = _parse(_s('BesselI[%r,-3.]' % nu))
        assert abs(z - complex(math.cos(math.pi * nu), math.sin(math.pi * nu)) * sp.iv(nu, 3.0)) < 1e-13
    # 整数次の J・I は実数
    assert _parse(_s('BesselJ[2.,-3.]')).imag == 0.0
    # Wronskian W[J, Y] = 2/(πz) が負の x でも成り立つ (主値の接続: J, Y とも同じ切断線の側)
    x = -2.5
    jn = lambda n: _parse(_s('BesselJ[%r,%r]' % (n, x)))
    yn = lambda n: _parse(_s('BesselY[%r,%r]' % (n, x)))
    nu = .7
    w = jn(nu + 1) * yn(nu) - jn(nu) * yn(nu + 1)
    assert abs(w - 2.0 / (math.pi * x)) < 1e-12


def test_gamma_incomplete_nonpositive_a_recurrence():
    # Γ(a+1, z) = a Γ(a, z) + z^a e^{-z} を独立に確かめる
    for a, z in ((-.5, 2.0), (-1.5, .3), (-2.0, .5), (0.0, 1.0), (-2.5, 4.0)):
        g0 = _parse(_s('Gamma[%r,%r]' % (a, z))).real
        g1 = _parse(_s('Gamma[%r,%r]' % (a + 1.0, z))).real
        assert abs(g1 - (a * g0 + z ** a * math.exp(-z))) < 1e-12 * abs(g1) + 1e-15, (a, z)
    # Γ(0, z) = E1(z)
    from scipy import special as sp
    assert abs(_parse(_s('Gamma[0,1.]')).real - sp.exp1(1.0)) < 1e-14
    # GammaRegularized[-n, z] = 0
    assert _parse(_s('GammaRegularized[-2.,3.]')).real == 0.0


def test_hyp1f1_complex_large_z_kummer_and_series_agree():
    # Kummer の変換 M(a,b,z) = e^z M(b-a,b,-z) (独立な経路)。従来の scipy は |z| = 25〜30 で破綻
    for a, b, z in ((.5, 1.5, '(0.1-30.*I)'), (1., 2., '(25.+3.*I)'), (.5, 1.5, '(12.-20.*I)')):
        m1 = _parse(_s('Hypergeometric1F1[%r,%r,%s]' % (a, b, z)))
        m2 = _parse(_s('Exp[%s]*Hypergeometric1F1[%r,%r,-%s]' % (z, b - a, b, z)))
        assert abs(m1 - m2) < 1e-9 * abs(m1), (a, b, z, m1, m2)


def test_legendre_nonint_degree_integer_order_matches_derivative_formula():
    # P_ν^m(x) = (-1)^m (1-x²)^{m/2} d^m/dx^m P_ν(x)。m = 1 を数値微分で
    nu, x, h = 2.5, .3, 1e-5
    p = lambda t: _parse(_s('LegendreP[%r,%r]' % (nu, t))).real
    d1 = (p(x + h) - p(x - h)) / (2 * h)
    want = -(1 - x * x) ** .5 * d1
    got = _parse(_s('LegendreP[2.5,1,.3]')).real
    assert abs(got - want) < 1e-8
    # k64 は NaN のまま
    assert _parse(_s('LegendreP[2.5,1,.3]', 'k64')) != _parse(_s('LegendreP[2.5,1,.3]', 'k64'))


def test_polylog_negative_integer_order_rational_closed_form():
    # Li_{-1}(z) = z/(1-z)²、Li_{-2}(z) = z(1+z)/(1-z)³、Li_{-3}(z) = z(1+4z+z²)/(1-z)⁴
    for z in (-.9, .5, 2.0, -7.0, 3.5):
        r = _evals(['PolyLog[-1.,%r]' % z, 'PolyLog[-2.,%r]' % z, 'PolyLog[-3.,%r]' % z])
        want = [z / (1 - z) ** 2, z * (1 + z) / (1 - z) ** 3, z * (1 + 4 * z + z * z) / (1 - z) ** 4]
        for g, w in zip(r, want):
            assert abs(_parse(g).real - w) < 1e-13 * abs(w)
    # Li_2(z) の小さい z は級数 z + z²/4 + ... (spence の桁落ちを避ける)
    z = 1e-10
    assert abs(_parse(_s('PolyLog[2.,1e-10]')).real - (z + z * z / 4)) < 1e-11 * z       # 出力は 12 桁 (指数が負で小さいとき)


def test_harmonic_number_small_argument():
    z = 1e-10
    want = math.pi ** 2 / 6 * z - 1.2020569031595942 * z * z
    assert abs(_parse(_s('HarmonicNumber[1e-10]')).real - want) < 1e-21


def test_productlog_branch_point_and_below():
    assert abs(_parse(_s('ProductLog[-.36787944117144233]')).real + 1.0) < 2e-8
    z = _parse(_s('ProductLog[-.5]'))
    assert abs(z * (2.718281828459045 ** z) - (-.5)) < 1e-14      # W e^W = x を独立に確かめる
    assert z.imag > 0


def test_lerch_phi_negative_z_integral_matches_series_and_functional_eq():
    # Φ(z,s,a) = a^{-s} + z Φ(z,s,a+1) (関数等式)。|z| ≥ 0.9 は積分表示、a+1 も同じ経路
    for z, s, a in ((-3., 2., 1.), (.95, 2., 1.5), (-.97, 2.5, 1.)):
        p0 = _parse(_s('LerchPhi[%r,%r,%r]' % (z, s, a))).real
        p1 = _parse(_s('LerchPhi[%r,%r,%r]' % (z, s, a + 1))).real
        assert abs(p0 - (a ** -s + z * p1)) < 1e-11 * abs(p0), (z, s, a)


def test_hyp2f1_above_cut_is_from_below():
    # 2F1(1/2,1/2;3/2;z) = ArcSin[√z]/√z。z = 2 は ArcSin[√2] = π/2 - i arccosh(√2) (切断線の下側)
    z = _parse(_s('Hypergeometric2F1[.5,.5,1.5,2.]'))
    want = complex(math.pi / 2, -math.acosh(math.sqrt(2.0))) / math.sqrt(2.0)
    assert abs(z - want) < 1e-12


def test_cosh_sinh_integral_complex_identities():
    # Shi(z) = -i Si(iz)、Chi は d/dz Chi = cosh(z)/z を数値微分
    z0 = 1.3 + .7j
    h = 1e-6
    ev = lambda f, z: _parse(_s('%s[(%r+%r*I)]' % (f, z.real, z.imag)))
    d = (ev('CoshIntegral', z0 + h) - ev('CoshIntegral', z0 - h)) / (2 * h)
    import cmath
    assert abs(d - cmath.cosh(z0) / z0) < 1e-8
    d = (ev('SinhIntegral', z0 + h) - ev('SinhIntegral', z0 - h)) / (2 * h)
    assert abs(d - cmath.sinh(z0) / z0) < 1e-8


# ------------------------------------------------------------------------------ 文法・評価器 (item 2)
def test_for_executes_body_v2_not_k64():
    src = 'For[i=0,i<3,i++,Print["b",i]];\nPrint["i=",i];'
    v2 = _run(src, 'v2')
    assert 'b 0' in v2.replace('b  ', 'b ') or 'b' in v2
    assert len([l for l in v2.splitlines() if l.startswith('b')]) == 3
    k = _run(src, 'k64')
    assert len([l for l in k.splitlines() if l.startswith('b')]) == 0      # gs は本体を実行しない (B71)


def test_trailing_operator_is_syntax_error_in_v2():
    assert 'unexpected' in _run('Print[1+];', 'v2') or 'syntax' in _run('Print[1+];', 'v2').lower()


def test_union_sametest_names_writebeamline_begin_context_v2():
    out = _run('Print[Union[{1,2,3,4,5,6},SameTest->(Abs[#-#2]<1.5&)]];', 'v2')
    nums = [float(t) for t in re.findall(_NUM, out.split('Print')[-1] if 'Print' in out else out)]
    assert nums[-3:] == [1.0, 3.0, 5.0]
    k = _run('Print[Union[{1,2,3,4,5,6},SameTest->(Abs[#-#2]<1.5&)]];', 'k64')
    assert [float(t) for t in re.findall(_NUM, k)][-1:] == [1.0]
    out = _run('Begin["Foo`"];xx=5;End[];Print[Names["Foo`*"]];', 'v2')
    assert 'Foo`xx' in out
    assert 'Context' in _run('Begin["Foo`"];', 'k64')


def test_polynomialfit_zero_mean_and_inverseerf_and_legendre_odd_m_v2():
    out = _run('Print[PolynomialFit[{{-1,2},{0,1},{1,2}},2]];', 'v2')
    nums = [float(t) for t in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', out.split('Print')[-1])][:3]
    assert abs(nums[0] - 1.0) < 1e-12 and abs(nums[1]) < 1e-12 and abs(nums[2] - 1.0) < 1e-12      # y = 1 + x² (k64 は {2, 2, 1})
    out = _run('Print[PolynomialFit[{{-1,2},{0,1},{1,2}},2]];', 'k64')
    assert '2.' in out
    assert abs(_parse(_s('InverseErf[.3]')).real - 0.2724627147267543) < 1e-15
    assert abs(_parse(_s('InverseErf[.3]', 'k64')).real - 0.272462711811) < 1e-12
    from scipy import special as sp
    # LegendreP[3,1,.4] = -(1-x²)^{1/2} P3'(x) = +0.2749545 (奇数 m の符号)
    assert abs(_parse(_s('LegendreP[3,1,.4]')).real - (-(1 - .16) ** .5 * (7.5 * .16 - 1.5))) < 1e-13
    assert sp is not None


# ------------------------------------------------------------------------------ リスト・文字列・パターン (item 3)
def test_insert_negative_position_v2_end_k64_gs():
    assert _evals(['Insert[{1,2,3},x,-1]'], 'v2')[0].replace(' ', '').endswith(',x}')       # Mathematica: 末尾
    assert _evals(['Insert[{1,2,3},x,-1]'], 'k64')[0].replace(' ', '').endswith(',x,3.000000000000000}')   # gs: 末尾の 1 つ手前
    assert _evals(['StringInsert["abc","X",-1]'], 'v2')[0] == 'abcX'
    assert _evals(['StringInsert["abc","X",-1]'], 'k64')[0] == 'abXc'            # gs
    assert _evals(['StringInsert["abc","X",2]'], 'v2')[0] == 'aXbc'


def test_table_over_explicit_list_v2():
    assert [_parse(t).real for t in _parse_list(_evals(['Table[i^2,{i,{1,3,7}}]'], 'v2')[0])] == [1.0, 9.0, 49.0]
    assert _parse(_s('Sum[i,{i,{1,3,7}}]')).real == 11.0
    assert 'wrongtype' in _run('Print[Table[i,{i,{1,3,7}}]];', 'k64')           # gs は受けない


def _parse_list(s):
    return [t for t in re.findall(_NUM, s)]


def test_take_drop_all_v2():
    assert _parse_list(_evals(['Take[{1,2,3},All]'], 'v2')[0]) == ['1.000000000000000', '2.000000000000000', '3.000000000000000']
    assert _evals(['Drop[{1,2,3},All]'], 'v2')[0] == '{}'
    assert 'wrongtype' in _run('Print[Take[{1,2,3},All]];', 'k64')


def test_inner_matches_mathematica_v2_and_gs_k64():
    assert _evals(['Inner[f,{1,2},{3,4},g]'], 'v2')[0].replace(' ', '').startswith('g[f[1.')
    assert _evals(['Inner[f,{1,2},{3,4},g]'], 'k64')[0].replace(' ', '').startswith('g[f[1.')
    # 行列 × ベクトル: Mathematica {17, 39}、gs は最上位の添字だけ対にするので {23, 34}
    assert [_parse(t).real for t in re.findall(_NUM, _evals(['Inner[Times,{{1,2},{3,4}},{5,6},Plus]'], 'v2')[0])] == [17.0, 39.0]
    assert [_parse(t).real for t in re.findall(_NUM, _evals(['Inner[Times,{{1,2},{3,4}},{5,6},Plus]'], 'k64')[0])] == [23.0, 34.0]
    assert [float(t) for t in re.findall(_NUM, _evals(['Inner[Times,{{1,2},{3,4}},{{1,0},{0,1}},Plus]'], 'v2')[0])] == [1.0, 2.0, 3.0, 4.0]
    assert _parse(_s('Inner[Times,{1,2},{3,4},Plus]')).real == 11.0


def test_mapindexed_level_spec_like_gs():
    for mode in ('k64', 'v2'):
        r = _evals(['MapIndexed[f,{{a},{b}},2]'], mode)[0].replace(' ', '')
        assert r.startswith('{f[{f[a,{1.')


def test_symbolic_logic_and_equal_stay_unevaluated_like_gs():
    for mode in ('k64', 'v2'):
        r = _evals(['x&&True', 'And[1,x]', 'Or[0,x]', 'x||False', 'If[x,1,2]', 'x==y', 'Head[x==y]', 'If[1,a,b]', 'And[]', 'Or[]'], mode)
        assert r[:4] == ['x', 'x', 'x', 'x']
        assert r[4].startswith('If[x')
        assert r[5] == 'x==y' or r[5].startswith('Equal[x,y')
        assert r[6] == 'Equal'
        assert r[7] == 'a'
    assert _evals(['x==x'], 'v2')[0] != 'x==x'            # v2: True
    assert _parse(_evals(['x==x'], 'v2')[0]).real == 1.0


def test_stringtake_zero_and_out_of_range():
    assert _evals(['StringTake["abc",0]'], 'v2')[0] == ''
    assert _evals(['StringTake["abc",2]'], 'v2')[0] == 'ab'
    assert _evals(['StringTake["abc",-2]'], 'v2')[0] == 'bc'
    assert _evals(['StringTake["abc",{2,3}]'], 'v2')[0] == 'bc'
    assert _evals(['StringTake["abc",5]'], 'v2')[0].startswith('StringTake[')


def test_complex_order_polylog_expintegrale_both_modes_like_gs():
    # gs: PolyLog[1.5+I,.3] = .327000552458529-.026234385397591 I、ExpIntegralE[1.5+I,2.] = .040081025470939-.010720907669603 I
    for mode in ('k64', 'v2'):
        z = _parse(_s('PolyLog[1.5+1.*I,.3]', mode))
        assert abs(z - complex(.327000552458529, -.026234385397591)) < 1e-13
        z = _parse(_s('ExpIntegralE[1.5+1.*I,2.]', mode))
        assert abs(z - complex(.040081025470939, -.010720907669603)) < 1e-13


def test_overflow_of_complex_trig_and_comparison_display_like_gs():
    for mode in ('k64', 'v2'):
        r = _evals(['Sin[1000.*I]', 'Cos[1000.*I]', 'Sinh[1000.+1.*I]', 'x<y', 'x>1&&y>1', 'x<y<z'], mode)
        assert r[0].replace(' ', '') == '(INFI)' and r[1].strip() == 'INF'
        z = _parse(r[2])
        assert z.real == float('inf') and z.imag == float('inf')
        assert r[3] == 'x<y' and r[4].replace(' ', '') == 'x>1.000000000000000&&y>1.000000000000000'.replace(' ', '') or r[4].startswith('x>')
        assert r[5] == '1.000000000000000&&x<y&&y<z' or r[5].replace(' ', '').startswith('1.0')
