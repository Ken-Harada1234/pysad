"""SADScript と入出力の、gs (K64) の実動作に合わせた箇所の試験。参照は gs で取った値
(Module の字句スコープ、Return の範囲、Switch/Which の一致なし、リストの ==、Out[n] と %、Close[{..}]、Write[$Output]、
DateString[date]、Read、文字列のエスケープ、StringMatchQ、Random[n,m]、SeedRandom)。"""
import io
import os
import re
import tempfile


def _run(body, echo=True):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=1999;\n' + body + '\nExit;\n', out=buf)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _lines(body):
    return [l for l in _run(body) if not l.startswith('Out[1]')]


def test_module_lexical_scope():
    r = _lines('''x = 10;
show[] := x;
m2[] := Module[{x = 7}, show[]];
Print["A ", m2[]];
Print["B ", Module[{x = 7, y = x}, y]];
Print["C ", Module[{x = 7}, Module[{x = x + 1}, x]]];
Print["D ", Module[{t}, t] === t, " ", Module[{x}, Head[x]]];
f[y_] := y + 1;
Print["E ", Module[{f}, f[x_] := x^2; f[4]], " ", f[4]];
mk[] := Module[{n = 0}, (n++ &)];
g = mk[]; g[]; Print["F ", g[]];
Print["G ", Module[{x = 7}, Map[# + x &, {1, 2}]]];
Print["H ", Module[{l}, l = {1, 2}; l[[2]] = 9; l]];
Print["I ", Module[{a, b}, {a, b} = {3, 4}; a + b]];
''')
    assert 'A 10' in r and 'B 10' in r and 'C 8' in r
    assert 'D 0 Symbol' in r
    assert 'E 16 5' in r
    assert 'F 1' in r or 'F 2' in r
    assert 'G {8,9}' in r and 'H {1,9}' in r and 'I 7' in r


def test_module_local_symbol_printed_with_suffix():
    r = _lines('Print["A ", Module[{t}, t]];')
    assert re.fullmatch(r'A t\$\d+', r[0])


def test_return_leaves_pure_function_and_loops():
    r = _lines('''f3[l_List] := Module[{r = 0}, Scan[If[# < 0, Return["negative found"]] &, l]; "none"];
Print["A ", f3[{1, -1}], " ", f3[{1, 2}]];
f4[] := (Map[If[# == 2, Return[#]] &, {1, 2, 3}]; "after");
Print["B ", f4[]];
f5[] := (Do[If[i > 2, Return[i]], {i, 5}]; "after");
Print["C ", f5[]];
f6[] := ((Return[3] &)[1]; 4);
Print["D ", f6[]];
g[] := Return[7];
f7[] := (g[]; 4);
Print["E ", f7[]];
h[] := (Module[{}, Return[1]]; 2);
Print["F ", h[]];
''')
    assert 'A negative found none' in r and 'B 2' in r and 'C 3' in r and 'D 3' in r and 'E 4' in r and 'F 1' in r


def test_return_outside_function_is_error():
    r = _lines('Print["a ", Return[3]];\nPrint["b"];')
    assert any(l.startswith('???General::unexpbreak') for l in r)
    assert 'b' in r and not any(l.startswith('a ') for l in r)


def test_switch_which_no_match_returns_expression():
    r = _lines('''Print["A ", Switch[9, 1, "one"]];
Print["B ", Switch[9, 1, "a", _String, "b"]];
a = 3;
Print["C ", Switch[a, 1, 2]];
Print["D ", Switch[{1, 2}, {1}, 1]];
Print["E ", Which[1 > 2, "a", 3 > 4, "b"]];
b = 5;
Print["F ", Which[b < 0, b]];
Print["G ", Switch[2, 1, "x", 2, "y"], Which[False, 1, True, 2]];
''')
    assert 'A Switch[9]' in r and 'B Switch[9]' in r
    assert 'C Switch[3,1,2]' in r and 'D Switch[{1,2}]' in r
    assert 'E Which[1>2,"a",3>4,"b"]' in r and 'F Which[b<0,b]' in r
    assert 'G y2' in r


def test_switch_which_wrong_arg_count():
    r = _lines('Print[Switch[3, 1]];\nPrint[Which[1]];')
    assert any('expected to be odd number in Switch[3,1]' in l for l in r)
    assert any('expected to be even number in Which[1]' in l for l in r)


def test_for_three_args_stays_unevaluated():
    r = _lines('j = 0; For[j = 0, j < 3, j++]; Print["A ", j];')
    assert 'A 0' in r


def test_list_equal_threads():
    r = _lines('''Print["A ", {1, 2} == {1, 3}];
Print["B ", Thread[{1, 2} == {1, 3}]];
Print["C ", {1, 2} == {1, 2}, {1, 2} == {1, 2, 3}];
Print["D ", {{1, 2}, {3}} == {{1, 2}, {4}}];
Print["E ", {1, 2} <> {0, 3}, {1, 2} >= 2, {1, 2} > {0, 3}];
Print["F ", {a, b} == {a, c}];
''')
    assert 'A {1,0}' in r and 'B {1,0}' in r
    assert 'C {1,1}0' in r and 'D {{1,1},{0}}' in r
    assert 'E {1,1}{0,1}{1,0}' in r and 'F {1,b==c}' in r


def test_out_and_percent():
    r = _run('''x = 5
x + 1;
2 * 3
Print["a ", %, " ", Out[1], " ", $Line];
{1, 2}
Print[Length[%], " ", %1, " ", Out[2]];
y = 4; z = 5
f[t_] := t
Print[$Line];
''')
    assert 'Out[1]:= 5' in r and 'Out[2]:= 6' in r
    assert 'a 6 5 2' in r
    assert 'Out[3]:= {1,2}' in r
    assert '2 5 6' in r
    assert 'Out[4]:= 5' in r and 'Out[5]:= t' in r and '5' in r


def test_close_list_and_write_output():
    with tempfile.TemporaryDirectory() as d:
        pa, pb = os.path.join(d, 'a.txt').replace('\\', '/'), os.path.join(d, 'b.txt').replace('\\', '/')
        r = _lines('''f = OpenWrite["%s"]; g = OpenWrite["%s"];
Write[f, "x"]; Write[g, "y"];
Print["A ", Close[{f, g}]];
Write[$Output, "w1 ", 3, {1, 2}];
WriteString[$Output, "ws\\n"];
Write[-1, "m1"];
Print["B ", Write[$Output, 1]];
Print[Write[$Output]];
''' % (pa, pb))
        assert open(pa).read().strip() == 'x' and open(pb).read().strip() == 'y'
    assert 'A Null' in r and 'w1 3{1,2}' in r and 'ws' in r and 'm1' in r
    assert '1' in r and 'B Null' in r
    assert any('expected to be 2 or more in Write[-1]' in l for l in r)


def test_datestring_with_argument():
    r = _lines('''Print[DateString[{2020, 1, 2, 3, 4, 5}]];
Print[DateString[1.5e9]];
Print[DateString[0]];
Print[DateString[{2020, 1, 2}]];
Print[DateString[1, 2]];
''')
    assert r[0] == '01/02/2020 03:04:05' and r[1] == '07/15/1947 02:40:00'
    assert r[2] == '01/01/1900 00:00:00' and r[3] == '01/02/2020 00:00:00'
    assert 'DateString$[1,2]' in r[4]


def test_read_words_lists_expressions():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, 'r.txt').replace('\\', '/')
        r = _lines('''a = OpenWrite["%s"];
Write[a, "alpha beta, gamma  12.5 7"];
Write[a, "x+y  3*4"];
Write[a, "1 2 3"];
Write[a, "4 5 6"];
Write[a, "end of file"];
Close[a];
c = OpenRead["%s"];
Print["A ", Read[c, Word]];
Print["B ", Read[c, Word]];
Print["C ", Read[c, {Word, Real, Real}]];
Print["D ", Read[c, Expression]];
Print["E ", Read[c, {Real, Real}]];
Print["F ", Read[c, {Real, Real, Real}]];
Print["G ", Read[c, Character]];
Print["H [", ReadString[c], "]"];
Print["I ", ReadString[c]];
Print["J ", ReadString[c]];
Close[c];
a = OpenWrite["%s"]; Write[a, "a:b:c"]; Write[a, "d:e"]; Close[a];
c = OpenRead["%s"];
Print["K ", Read[c, Word, WordSeparators -> ":"], Read[c, Word, WordSeparators -> ":"]];
Close[c];
c = OpenRead["%s"];
Print["L ", Read[c, Word, ReadNewRecord -> False, WordSeparators -> ":"], " ", Read[c, Word, ReadNewRecord -> False, WordSeparators -> ":"]];
Close[c];
c = StringToStream["10 20 30\\nnext line\\n3.5 abc"];
Print["M ", Read[c, 3*Real]];
Print["N ", ReadString[c]];
Print["O ", Read[c, {Real, Word}]];
Print["P ", Read[c, Word]];
c = StringToStream["{1, 2, x+1} f[3]"];
Print["Q ", Read[c, Expression], " ", Read[c, Expression]];
''' % (p, p, p, p, p))
    assert 'A alpha' in r and 'B beta' in r
    assert 'C {"gamma",12.5,7}' in r and 'D x+y' in r
    assert 'E {12,1}' in r and 'F {2,3,4}' in r and 'G 5' in r
    assert 'H [ 6]' in r and 'I end of file' in r and 'J EndOfFile' in r
    assert 'K ab' in r and 'L a b' in r
    assert 'M 102030' in r and 'N next line' in r and 'O {3.5,"abc"}' in r and 'P EndOfFile' in r
    assert 'Q {1,2,1+x} f[3]' in r


def test_string_escapes():
    r = _lines('''Print[ToCharacterCode["\\101"]];
Print[ToCharacterCode["\\101B\\60"]];
Print[ToCharacterCode["\\1234"]];
Print[ToCharacterCode["\\18"]];
Print[ToCharacterCode["\\8"]];
Print[ToCharacterCode["\\777"]];
Print[ToCharacterCode["\\x41"]];
Print[ToCharacterCode["\\xg"]];
Print[ToCharacterCode["a\\qb\\e\\f"]];
''')
    assert r[0] == '{65}' and r[1] == '{65,66,48}' and r[2] == '{83,52}' and r[3] == '{1,56}'
    assert r[4] == '{56}' and r[5] == '{255}' and r[6] == '{65}' and r[7] == '{0,103}'
    assert r[8] == '{97,113,98,27,12}'


def test_stringmatchq_wildcards():
    cases = [('Q1', 'Q%', 1), ('Q12', 'Q%', 0), ('Q12', 'Q%%', 1), ('Q', 'Q%', 0), ('QF', 'Q{FD}', 1), ('QX', 'Q{FD}', 0),
             ('QF', 'Q{^FD}', 0), ('QX', 'Q{^FD}', 1), ('QFD', 'Q{FD}*', 1), ('QFD', 'Q{FD}{FD}', 1), ('QF', 'Q?', 0),
             ('Q?', 'Q?', 1), ('QF1', '*F%', 1), ('QD', 'QF|QD', 1), ('QD', 'QF|Q*', 1), ('qf', 'QF', 0), ('', '*', 0),
             ('a', '{a-c}', 1), ('b', '{a-c}', 0), ('a{b', 'a{b', 0), ('a', '{a', 1), ('Q', 'Q*', 1), ('ab', 'a|b', 0)]
    body = ''.join('Print["%s|%s|", StringMatchQ["%s", "%s"]];\n' % (s, p, s, p) for s, p, _ in cases)
    r = _lines(body)
    for (s, p, v), l in zip(cases, r):
        assert l == '%s|%s|%d' % (s, p, v), (s, p, l)


def test_random_dimensions():
    r = _lines('''SeedRandom[1];
Print[Dimensions[Random[2, 3]], Dimensions[Random[2, 3, 4]], Random[0], Random[3, 0], Dimensions[Random[3, 2, 2]]];
Print[Dimensions[GaussRandom[2, 3]], Dimensions[ParabolaRandom[2, 3]], GaussRandom[0]];
SeedRandom[1];
a = Random[2, 3]; SeedRandom[1]; b = Random[6];
Print[Flatten[a] === b];
Print[Random[-1]];
Print[Random[{2, 3}]];
''')
    assert r[0] == '{2,3}{2,3,4}{}{{},{},{}}{3,2,2}'
    assert r[1] == '{2,3}{2,3}{}' and r[2] == '1'
    assert 'wrongnum' in r[3] and 'wrongtype' in r[5]


def test_seedrandom_state_forms():
    r = _lines('''SeedRandom[5];
Print[SeedRandom[]];
Print[SeedRandom["SAD"]];
SeedRandom["x"];
Print[SeedRandom[{"SAD", 5, 0, 0, 0}]];
SeedRandom[-5]; Print[SeedRandom[]];
SeedRandom[{7}]; Print[SeedRandom[]];
SeedRandom[3]; g = GaussRandom[]; Print[SeedRandom[]];
SeedRandom[1, 2];
''')
    assert r[0] == '{"SAD",5,0,0,0}' and r[1] == '{"SAD",5,0,0,0}'
    assert 'No such random plugin in SeedRandom["x"]' in r[2]
    assert r[3] == '{"SAD",5,0,0,0}' and r[4] == '{"SAD",4294967291,0,0,0}' and r[5] == '{"SAD",1,0,0,0}'
    assert r[6] == '{"SAD",2306130011,1,787469756,1071741574}'


def test_user_function_named_At_does_not_clash_with_member_ref():
    """ユーザー定義の関数名 At は、メンバー参照 obj@name の内部の頭と別 (gs: 101)。"""
    r = _lines('At[x_] := x + 100; Print["A ", At[1]];')
    assert 'A 101' in r
    r = _lines('Cc = Class[{}, {}, {x = 3}]; o = Cc[]; At[x_] := x + 100; Print["B ", o@x, " ", At[2]];')
    assert 'B 3 102' in r
