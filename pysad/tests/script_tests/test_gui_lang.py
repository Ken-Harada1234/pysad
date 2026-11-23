"""Tkinter.n / KBFrame.n を動かすために足した SADScript の言語機能 (SAD K64 の実測に合わせたもの)。
参照値は SAD (WSL の gs) で perf_work の lst_probe*.sad / pt.sad を実行して取った。GUI そのもの (tkinter) の試験は test_tclbridge で、
こちらは GUI を使わない: Sym::tag、f@@[x, {1}]、Null[]、0x33、UpValue / SubValue / Hold 属性、Difference / SwitchCases / StringFill / Symbol。"""
import io

import pytest

from pysad.sadscript import Interp


def _run(src):
    I = Interp()
    out = io.StringIO()
    I.out = out
    I.run(src, 'script')
    return I, out.getvalue()


def test_list_and_string_functions_match_sad():
    I, _ = _run('''
d1 = Difference[{1,4,9,16}];
s1 = SwitchCases[{1,a,2,b,3.5,"s"},{_Real,_Symbol}];
s2 = SwitchCases[{a->1,b:>2,3,4},{_Rule|_RuleDelayed,_}];
s3 = SwitchCases[{1,2,3,4},{_?(#>2&),_?(#>1&),_}];
t1 = StringFill["ab","xyz",7];
t2 = StringFill["ab","xyz",-7];
t3 = StringFill["ab","xyz",0];
Symbol["abc"] = 5;
y3 = abc;
y4 = Symbol["x"//ToString[3]];
''')
    assert I.get('d1') == [3.0, 5.0, 7.0]
    assert len(I.get('s1')) == 2 and I.get('s1')[0] == [1.0, 2.0, 3.5]
    assert [len(x) for x in I.get('s2')] == [2, 2]
    assert [len(x) for x in I.get('s3')] == [2, 1, 1]
    assert I.get('t1') == 'abxyzxy' and I.get('t2') == 'yzxyzab' and I.get('t3') == ''
    assert I.get('y3') == 5.0


def test_difference_needs_two_elements():
    from pysad.sadscript import SADError
    with pytest.raises(SADError, match='wrongleng'):
        _run('x = Difference[{5}];')


def test_message_name_syntax_and_trailing_null_sequence():
    I, out = _run('''
KBF::main = "can be called only once";
Off[General::unprotconst];
m = KBF::main;
l = Join[{1,2}, {3}, Null[]];
n = {1, 2, Null[]};
h = 0x33 + 0xff;
''')
    assert I.get('m') == 'can be called only once'
    assert I.get('l') == [1.0, 2.0, 3.0] and I.get('n') == [1.0, 2.0]
    assert I.get('h') == 51.0 + 255.0


def test_apply_with_level_spec_and_empty_rule_rhs():
    I, _ = _run('''
a1 = f@@[{{1,2},{3}}, {1}];
r = {x->, y->3};
''')
    assert len(I.get('a1')) == 2


def test_chained_unsame():
    I, _ = _run('v1 = ("***" <=> (ids = "ab") <=> Null); v2 = ("***" <=> "***" <=> Null);')
    assert I.get('v1') == 1.0 and I.get('v2') == 0.0


def test_hold_attributes_and_held_matching():
    I, _ = _run('''
SetAttributes[hk, HoldAll];
hk[x_] := ToString[Hold[x]];
tag[x_] := x;
SetAttributes[pn, HoldAll];
pn[_] = "";
pn[root] = ".";
w1 = root;
r1 = pn[root];
r2 = pn[w1];
''')
    assert I.get('r1') == '.' and I.get('r2') == ''


def test_set_upvalue_intercepts_assignment():
    I, _ = _run('''
(w2_ = mk[opt___]) ^:= (w2 = made[opt]);
v = mk[1, 2];
''')
    assert str(I.get('v')).startswith('made[')


def test_subvalues():
    I, _ = _run('''
wd[b_, _][c_] := {b, c};
r = wd[1, 2][3];
ReturnTo[$W][Title] = "T";
t = ReturnTo[$W][Title];
''')
    assert I.get('r') == [1.0, 3.0]
    assert I.get('t') == 'T'


def test_baseform_and_string_substring():
    I, _ = _run('b = BaseForm[161, 36]; z = "abcdef"[2, 4];')
    assert I.get('b') == '4h' and I.get('z') == 'bcd'


def test_unset_removes_only_that_definition():
    I, _ = _run('''
f[1] = a; f[2] = b;
f[1] =.;
r = {f[1], f[2]};
''')
    r = I.get('r')
    assert str(r[1]) == 'b'
    assert str(r[0]).startswith('f[1')
