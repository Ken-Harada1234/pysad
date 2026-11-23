"""SAD の Class (pysad/sadclass.py) と、KBFrame / Tkinter.n を動かすために足した言語機能 (GUI を使わない部分)。
参照の挙動は SADHelp の Class の節と、WSL の SAD (gs) で確かめた値による。"""
import io

import pytest

from pysad.sadscript import Interp, SADError


def _run(src):
    I = Interp()
    out = io.StringIO()
    I.out = out
    I.run(src, 'script')
    return I, out.getvalue()


def test_class_instance_members_and_methods():
    I, _ = _run('''
cA = Class[{}, {cnt=0}, {v1=1, v2=0},
  Constructor[args__] := (v1 = Plus[args]; cnt++);
  sum[] := v1 + v2;
  setv2[x_] := (v2 = x; This);
  bump[] := (v1++; Helper[]);
  Helper[] := v1*10];
a1 = cA[1, 2, 3];
r1 = a1@v1;
a1@v2 = 5;
r2 = a1@sum[];
a2 = cA[10];
r3 = a2@v1;
r4 = a1@v1;
r5 = a1@bump[];
r6 = a1@cnt;
''')
    assert [I.get(k) for k in ('r1', 'r2', 'r3', 'r4', 'r5', 'r6')] == [6.0, 11.0, 10.0, 6.0, 70.0, 2.0]


def test_class_inheritance_and_method_chain():
    I, _ = _run('''
cA = Class[{}, {}, {v1=1, v2=0}, Constructor[x_] := (v1 = x); sum[] := v1 + v2; setv2[x_] := (v2 = x; This)];
cB = Class[{cA}, {}, {w=7}, sum[] := w + v1];
b1 = cB[4];
r1 = b1@sum[];
r2 = b1@setv2[3]@v2;
''')
    assert I.get('r1') == 11.0 and I.get('r2') == 3.0


def test_class_rules_set_instance_variables_and_options():
    I, _ = _run('''
cA = Class[{}, {}, {v1=1}, Constructor[] := (n = Length[Options]); getopt[] := Options];
a1 = cA[v1->9, Color->"red"];
r1 = a1@v1;
r2 = a1@getopt[];
''')
    assert I.get('r1') == 9.0
    assert len(I.get('r2')) == 1


def test_class_instances_do_not_share_instance_variables_and_destructor():
    I, out = _run('''
cA = Class[{}, {}, {v1=0}, Constructor[x_] := (v1 = x); Destructor[] := Print["bye ", v1]];
a1 = cA[1]; a2 = cA[2];
a1@v1 = 10;
r1 = {a1@v1, a2@v1};
a2 =.;
''')
    assert I.get('r1') == [10.0, 2.0]
    assert 'bye 2' in out


def test_at_operator_keeps_function_application_for_non_instances():
    I, out = _run('r1 = Print@Range[2]; r2 = f@g[x];')
    assert '{1,2}' in out.replace(' ', '')


def test_null_sequence_splice():
    I, _ = _run('r1 = f[a, Null@@{b, c}, d]; r2 = f[Null[x, y], z]; r3 = f[a, Null@@{}, d];')
    assert str(I.get('r1')) == 'f[a, b, c, d]' and str(I.get('r2')) == 'f[x, y, z]' and str(I.get('r3')) == 'f[a, d]'


def test_and_or_short_circuit():
    I, _ = _run('n = {}; a1 = (n<=>{}&&n[[1,2]]); a2 = (n===={}||n[[1]]);'.replace('====', '==='))
    assert I.get('a1') == 0.0 and I.get('a2') == 1.0


def test_with_semicolon_declaration_and_shorthand():
    I, _ = _run('''
x = 5; y = 6;
r1 = With[{a = x; b = y + 1, c = 3}, {a, b, c}];
r2 = With[{x, y}, x + y];
''')
    assert I.get('r1') == [5.0, 7.0, 3.0] and I.get('r2') == 11.0


def test_rule_object_patterns_and_replaceall_patterns():
    I, _ = _run('''
m1 = SwitchCases[{"a"->1, b:>2, 3}, {Rule[_String, _], RuleDelayed[_Symbol, _], _}];
t1 = {{F->Real,"*"},{S->String,"x"},{Slant->Q,"r"}} /. {{x_->Real,"*"}:>(x->"*"), {Slant->_, y_}:>(Slant->y), {x_->_,y_}:>(x->y)};
t2 = {1,2} /. n_Real -> n + 10;
''')
    assert [len(x) for x in I.get('m1')] == [1, 1, 1]
    assert str(I.get('t1')) == "[F->'*', S->'x', Slant->'r']"
    assert I.get('t2') == [11.0, 12.0]


def test_replaceall_does_not_evaluate_delayed_rule_rhs():
    I, out = _run('''
r = {Command:>Print["called"]} /. Text->Label;
''')
    assert 'called' not in out


def test_tagset_and_three_level_subvalues_and_unevaluated_dollar():
    I, _ = _run('''
a/:f[a, x_] := {x};
r1 = f[a, 3];
cf[1][c][Bg] = 5;
r2 = cf[1][c][Bg];
inc[i_] := (i++; i += 3; k[i] = 1; i++);
i = 5; inc[Unevaluated$[i]];
r3 = i; r4 = k[9];
''')
    assert I.get('r1') == [3.0] and I.get('r2') == 5.0 and I.get('r3') == 10.0 and I.get('r4') == 1.0


def test_hold_keeps_its_argument_and_union_of_same_heads():
    I, _ = _run('h = Hold[x + 1]; u = Union[Hold[], Hold[a]];')
    assert str(I.get('h')).startswith('Hold[') and str(I.get('u')) == 'Hold[a]'


def test_rgbcolor_and_date_string_format():
    I, _ = _run('c1 = RGBColor[1, 0.5, 0]; c2 = RGBColor[{{0,0,0},{1,1,1}}]; d = DateString[];')
    assert I.get('c1') == '#FF8000' and I.get('c2') == ['#000000', '#FFFFFF']
    import re
    assert re.fullmatch(r'\d\d/\d\d/\d{4} \d\d:\d\d:\d\d', I.get('d'))


def test_local_assignment_does_not_erase_global_definitions():
    I, _ = _run('''
w[1] = 10;
f[wt_] := Module[{w = wt}, w = w + 1; w];
r1 = f[5];
r2 = w[1];
''')
    assert I.get('r1') == 6.0 and I.get('r2') == 10.0


def test_set_with_switch_pattern_in_variable_and_typed_blank():
    I, _ = _run('''
types = _Real;
s1 = Switch[12, types, "real", _, "other"];
s2 = Switch["a", _Real, "real", _String, "string"];
''')
    assert I.get('s1') == 'real' and I.get('s2') == 'string'


def test_standard_form_ignores_script_form_setting():
    """init.n の StandardForm: 台本が $FORM = "6.2" としていても、StandardForm[...] の中では既定の書式 (Tcl へ送る命令に使う)"""
    I, _ = _run('$FORM = "6.2"; a = "x" // 1001; b = StandardForm["x" // 1001];')
    assert I.get('a') == 'x 1001.' and I.get('b') == 'x1001'
