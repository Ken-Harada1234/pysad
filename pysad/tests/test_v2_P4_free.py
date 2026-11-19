# -*- coding: utf-8 -*-
"""v2next 担当 P4 項目 6: 要素ごとの FREE (`FREE QF.3 K1`) の「何もしない変数」。gs (K64、2026-10-19 ビルド、2026-11-19) の実測と全桁一致。

K64 の規則 (tffsfreefix の comp 分岐):
  - 既定のキー (QUAD の K1) の成分の変数は、要素全体の値の別名。先頭の成分 `QF.1` の変数は要素全体と同じ (ラベル QF.1)。先頭でない成分 `QF.3` の変数は
    結合先 `QF.1` を持つ何もしない変数で、フィットの変数の数 (反復の行の最後の数) と VAR の行に数えられ、値は先頭に追従する。
  - その変数だけで GO しても `?No variable.` は出ない (反復なし、残差の行だけ)。
  - 要素全体の変数がある所へ成分 (FREE QF K1 のあと FREE QF.3 K1) はエラー `Element already used as variable: QF`、成分のある所へ要素全体
    (FREE QF.3 K1 のあと FREE QF K1) はエラー `A component has been already used as variable: QF.3`。Info は変数を作るときだけ。
  - VAR の Saved の欄は設計値 (定義の値か SAVE の値)。GO では変わらず、RESET で Now が設計値に戻る。
使った条件: 3 GeV、周長 18.4 m の小さな FODO (QF・QD 各 4 台、BEND 8 台)。実機の大きなリングでは確かめていない。モードによらない (K64 の動作の再現)。
"""
import io
import re

import pytest

_LAT = '''MOMENTUM = 3 GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
MARK ST = ();
BEND B = (L=0.5 ANGLE=0.1);
LINE CELL1 = (QF D2 B D1 QD D2 B D1);
LINE RING = (ST 4*CELL1);
'''


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(_LAT + 'FFS;\nPageWidth=1999;\nUSE RING; CELL; CALC;\n' + body + '\nExit;\n', out=buf)
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _rows(lines):
    """VAR の表の行 (ラベル、キー、Now、Previous、Saved、結合先)"""
    out, cur = [], None
    for l in lines:
        if l.startswith('!Variable'):
            cur = []
            out.append(cur)
        elif cur is not None and ' !' in l and re.match(r' \S+\s+\S+\s+-?[.\d]', l):
            t = l.replace('!', ' ').split()
            cur.append((t[0], t[1], float(t[2]), float(t[3]), float(t[4]), t[7] if len(t) == 9 else ''))
    return out


def test_non_head_component_variable_is_a_noop_variable_counted_like_gs(k64_mode):
    # gs: VAR の行 `QF.3 K1 .7 ! .7 .7 … QF.1 1.0`。その変数だけの GO は ?No variable. を出さず、反復なしで残差 0.7344 のまま
    r = _run('FIT $$$ NX 1.30 NY 0.95;\nFREE QF.3 K1;\nVAR;\nGO;\nVAR;')
    assert ' ???-FFS-Info-Component QF.3 is coupled to QF.1 .' in r
    assert not any('No variable' in l for l in r)
    t = _rows(r)
    assert [x[0] for x in t[0]] == ['QF.3'] and t[0][0][5] == 'QF.1'
    assert t[0][0][2:5] == (0.7, 0.7, 0.7)
    res = [l for l in r if l.startswith(' Residual = ')]
    assert res and res[0].startswith(' Residual = 0.7344 ')
    assert not any('(NEWTON)' in l or '(DESCEND)' in l for l in r)
    assert [x[0] for x in t[1]] == ['QF.3']


def test_head_and_non_head_component_variables_fit_like_gs(k64_mode):
    # gs: FREE QF.3 K1; FREE QF.1 K1; FIT NX 1.30; GO → 変数 2 つ (後から FREE した QF.1 が先)。QF.1 = QF.3 = .743599255507、残差 5.3438E-02、反復の行の最後の数は 2
    # (gs の f2 の続き: 先に QD と QF.3 で NX 1.30・NY 0.95 に GO して QD を動かしてから、QF.3・QF.1 で NX 1.30 に GO する)
    r = _run('FREE QF.3 K1; FREE QD K1; FIT $$$ NX 1.30 NY 0.95; GO;\nFIX *;\nFREE QF.3 K1; FREE QF.1 K1; FIT $$$ NX 1.30;\nGO;\nVAR;')
    t = _rows(r)[-1]
    assert [x[0] for x in t] == ['QF.1', 'QF.3'] and t[1][5] == 'QF.1'
    assert t[0][2] == pytest.approx(0.743599255507, abs=2e-12) and t[1][2] == pytest.approx(0.743599255507, abs=2e-12)
    assert t[0][3] == 0.7 and t[0][4] == pytest.approx(0.7, abs=1e-9) and t[1][4] == pytest.approx(0.7, abs=1e-9)    # Saved は設計値のまま
    it = [l for l in r if '(NEWTON)' in l]
    assert it and it[-1].split()[-1] == '2' and it[-1].split()[2] == '5.3746E-02'          # 最後の GO の反復 (最初の GO は 0.5323)
    assert any(l.startswith(' Residual = 5.3438E-02') for l in r)


def test_element_and_component_variable_conflicts_like_gs(k64_mode):
    r = _run('FREE QF K1; FREE QF.3 K1; VAR;')
    assert ' ???-FFS-Error-Element already used as variable: QF' in r
    assert [x[0] for x in _rows(r)[0]] == ['QF']
    r = _run('FREE QF.3 K1; FREE QF K1; VAR;')
    assert ' ???-FFS-Error-A component has been already used as variable: QF.3' in r
    assert [x[0] for x in _rows(r)[0]] == ['QF.3']
    r = _run('FREE QF.1 K1; FREE QF K1; VAR;')
    assert ' ???-FFS-Error-A component has been already used as variable: QF.1' in r
    r = _run('FREE QF.1 K1; FREE QF.1 K1; VAR;')              # 2 回目は何も出さない (Info は変数を作るときだけ)
    assert len([l for l in r if 'Info-Component' in l]) == 3 and [x[0] for x in _rows(r)[0]] == ['QF.1']


def test_fix_removes_only_the_named_component_variable_like_gs(k64_mode):
    r = _run('FREE QF.3 K1; FREE QF.1 K1;\nFIX QF.1 K1;\nVAR;\nFIX QF.3 K1;\nVAR;')
    t = _rows(r)
    assert [x[0] for x in t[0]] == ['QF.3'] and t[1] == []
    r = _run('FREE QF.3 K1; FREE QD K1;\nFIX QF;\nVAR;\nFIX QF.*;\nVAR;')           # FIX QF (要素全体の語) は成分の変数に触れない。QF.* は消す
    t = _rows(r)
    assert [x[0] for x in t[0]] == ['QF.3', 'QD'] and [x[0] for x in t[1]] == ['QD']


def test_var_saved_column_is_the_design_value_and_reset_restores_now_like_gs(k64_mode):
    # gs: GO では Saved は設計値のまま (.7 と −.6)、SAVE で現在値になり、RESET で Now が Saved に戻る
    r = _run('FIT $$$ NX 1.30 NY 0.95;\nFREE QF QD K1;\nGO;\nVAR;\nSAVE;\nVAR;\nFIT $$$ NX 1.31;\nGO;\nRESET;\nVAR;')
    t = _rows(r)
    assert t[0][0][4] == pytest.approx(0.7, abs=1e-9) and t[0][1][4] == pytest.approx(-0.6, abs=1e-9)
    assert t[0][0][2] == pytest.approx(0.747774120771, abs=1e-9) and t[0][1][2] == pytest.approx(-0.642286021422, abs=1e-9)
    assert t[1][0][4] == pytest.approx(t[1][0][2], abs=1e-9)                           # SAVE のあと Saved = Now
    assert t[2][0][2] == pytest.approx(t[1][0][4], abs=1e-9) and t[2][1][2] == pytest.approx(t[1][1][4], abs=1e-9)    # RESET のあと Now = Saved


def test_nokeepexp_component_rule_is_unchanged(k64_mode):
    # 既定のキーの成分の変数 (QF.1・QF.3) だけでは展開を止めない (tfvcomp は既定でないキーの成分の変数があるとき真)
    r = _run('NOKEEPEXP;\nFIT $$$ NX 1.30;\nFREE QF.3 K1; FREE QF.1 K1;\nGO;')
    assert any(l.startswith(' Matched.') or l.startswith(' Residual') for l in r)
    assert not any('not expanded' in l for l in r)
