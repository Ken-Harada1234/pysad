# -*- coding: utf-8 -*-
"""担当 P5 (v2next) 項目 9: MAIN の要素定義の式の文法の細部 (docs/47 §(d) の D2 の 3 件) を gs と比べる。

gs (K64) の値は 2026-11-19 の 10-19 ビルドの実測 (`BEND B1=(L=1 ANGLE=式)` のあと `Element["ANGLE","B1"]`。gs はエラーのとき値を返さない)。
  D1 `45DEG` (数と単位の間に空白が無い)・`SIN[0.1]` (関数は丸括弧): gs は構文エラー。pysad は前者を受け (従来の読み取り)、後者は値 0 を入れていた。
     → `SIN[0.1]` は両モードで構文エラーにした (pysad の誤りの修正)。`45DEG` は k64 が gs と同じ構文エラー、v2 は受ける (入力の寛容さ。値は 45 度)。
  D2 MAIN の変数 `E`・`Pi`・`Degree` への代入 (`E=0.5;`): gs は要素定義の式でその値を使う。pysad は SADScript の定数 (e・π) のままだった → 両モードで直した
     (MAIN の代入は別の辞書に入れ、要素定義の式だけが見る。FFS の SADScript の `E`・`Pi` は変わらない)。
  D3 末尾の演算子 (`1+`・`2*3+`): gs は構文エラーにせず値を 1 にする (`2*3+` も 6 でなく 1)。k64 は同じ (1)、v2 は構文エラー (入力ミス)。
未解決: 未定義の名前 (`ANGLE=ZZZ`) は gs・pysad とも定義が失敗するがエラーの文言が違う。MAIN で代入していない `E` (`ANGLE=E`) は gs が構文エラー、pysad は定数 e を使う。
"""
import pytest

from . import _p5_util as U

ERR = None
# (名前, 事前の MAIN の代入, 式, gs の値 (None: エラーで値が出ない), v2 の値)
CASES = [
    ('e_var', 'E=0.5;', 'E', 0.5, 0.5),
    ('pi_var', 'Pi=3;', 'Pi/10', 0.3, 0.3),
    ('deg_var', 'Degree=0.1;', 'Degree', 0.1, 0.1),
    ('a_var', 'A=2;', 'A*0.1', 0.2, 0.2),
    ('trail_op', '', '1+', 1.0, ERR),
    ('trail_op2', '', '2*3+', 1.0, ERR),
    ('unit_nospace', '', '45DEG', ERR, 0.7853981633974483),
    # gs は D で始まる単位だけを D 指数 (`45D` + `EG`) と読んで構文エラーにする。D で始まらない単位は読む (2026-11-21 の gs の実測。
    # 数値は ANGLE に入れるので意味は無いが、読めるかどうかの確認)
    ('unit_nospace_mv', '', '1.7MV', 1700000.0, 1700000.0),
    ('unit_nospace_gev', '', '3GEV', 3000000000.0, 3000000000.0),
    ('unit_nospace_mm', '', '10MM', 0.01, 0.01),
    ('brackets', '', 'SIN[0.1]', ERR, ERR),
    ('ctrl_unit', '', '3 DEG', 0.05235987755982989, 0.05235987755982989),
    ('pi_upper', '', 'PI/10', 0.3141592653589793, 0.3141592653589793),
    ('lower_pi', '', 'pi/10', ERR, ERR),
    ('undef_expr', '', 'ZZZ*2', ERR, ERR),
]


def _angle(pre, expr, mode):
    src = (f'MOMENTUM=3 GEV;\n{pre}\nMARK M0=();\nBEND B1=(L=1 ANGLE={expr});\nLINE R=(M0 B1);\nFFS;\nPageWidth=1999;\nUSE R;\n'
           f'Print["RES ",Element["ANGLE","B1"]];\nExit;\n')
    out = U.run(src, mode)
    l = [x for x in out.splitlines() if x.startswith('RES ')]
    return float(l[0][4:]) if l else None


@pytest.mark.parametrize('name,pre,expr,gsv,v2v', CASES)
def test_main_expression_k64_equals_gs_and_v2_as_documented(name, pre, expr, gsv, v2v):
    k = _angle(pre, expr, 'k64')
    if gsv is None:
        assert k is None, (name, k)
    else:
        assert k == pytest.approx(gsv, rel=1e-12), (name, k)
    v = _angle(pre, expr, 'v2')
    if v2v is None:
        assert v is None, (name, v)
    else:
        assert v == pytest.approx(v2v, rel=1e-12), (name, v)


def test_main_user_value_of_e_does_not_change_the_sadscript_constant_in_ffs():
    """MAIN の `E=0.5;` は要素定義の式だけに効き、FFS の SADScript の E (自然対数の底) は変わらない (gs の実測: `ANG .5 EE 2.718281828459045`)"""
    src = '''MOMENTUM=3 GEV;
E=0.5;
MARK M0=();
BEND B1=(L=1 ANGLE=E);
LINE R=(M0 B1);
FFS;
USE R;
Print["ANG ",Element["ANGLE","B1"]," EE ",E];
Exit;
'''
    for mode in ('k64', 'v2'):
        out = U.run(src, mode)
        l = [x for x in out.splitlines() if x.startswith('ANG ')][0]
        v = U.floats(l[4:].replace('EE', ' '))
        assert v[0] == pytest.approx(0.5) and v[1] == pytest.approx(2.718281828459045)
