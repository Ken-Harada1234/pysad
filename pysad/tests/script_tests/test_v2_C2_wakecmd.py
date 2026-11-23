# -*- coding: utf-8 -*-
"""v2 (C2): 旧式の WAKE コマンド (tfwake.f の構文) の表を行をまたいで書ける。

K64: tfwake.f は呼ばれない (WAKECOMMAND[] を評価して失敗する)。tfwake.f 自身の読み方は、表が行をまたいでよく、ENDWAKE で終わる。
v2: `WAKE L|T 要素` の文は ENDWAKE (または ;) まで次の行に続く。1 行に何組でも、注釈 `!` 可、ENDWAKE の後の `;` は任意。
根拠: 同じ表を `WakeFunction[…] = …` で与えた追跡 (test_wake.py で gs と一致済み) と同じ結果。
"""
import io
import math

from pysad.tests.script_tests.test_wake import C3_SRC, C3_REF

MULTI = '''WAKE LONGITUDINAL RF
  -0.1   2e9      ! 注釈: 1 行に 1 組
  -0.02  1e9   0  5e8        ! 1 行に 2 組
  0.02 -1e9
  0.1 -2e9 0.2
  -1e9                       ! 組の途中で改行
ENDWAKE;
wake t RF        ! 小文字・省略形
  0 1e12  0.01 5e11
  0.05 1e11  0.2 0  0.3 0
ENDWAKE
'''


def _src(old, base=C3_SRC):
    lines = base.splitlines()
    out = [l for l in lines if not l.startswith('WakeFunction[')]
    k = out.index('PBUNCH = 1E12; NBUNCH = 1;')
    return '\n'.join(out[:k + 1]) + '\n' + old + '\n'.join(out[k + 1:]) + '\n'


def _rows(src):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string(src, out=buf)
    rows = []
    for l in buf.getvalue().splitlines():
        t = l.split()
        if len(t) == 5 and t[0].replace('.', '').isdigit():
            rows.append([float(v) for v in t[1:]])
    return rows, buf.getvalue()


def _same(rows, ref, tol):
    assert len(rows) == len(ref) == 6
    for r, q in zip(rows, ref):
        for a, b in zip(r, q):
            assert abs(a - b) <= tol * max(abs(b), 1e-3)


def _ref_same_mode():
    """C3_SRC (同じ表を WakeFunction で与えた台本) を今の互換モードで走らせた行。C3_REF は gs (K64 の 2 次の四極の積分器) の値で、
    v2 の既定 (4 次の合成。B226) の追跡とは 1e-9 台ずれるため、比べる相手は同じ積分器の WakeFunction の結果にする"""
    return _rows(C3_SRC)[0]


def test_multiline_table_equals_wakefunction(v2_mode):
    rows, out = _rows(_src(MULTI))
    assert '?' not in out
    _same(rows, _ref_same_mode(), 2e-9)


def _table(n, kind):
    if kind == 'L':
        return [(-0.1 + 0.4 * j / (n - 1), 2e9 * math.cos(9.0 * j / (n - 1)) * math.exp(-2.0 * j / (n - 1))) for j in range(n)]
    return [(0.4 * j / (n - 1), 1e12 * math.exp(-12.0 * j / (n - 1))) for j in range(n)]


def test_table_of_hundreds_of_points(v2_mode):
    n = 400
    tl, tt = _table(n, 'L'), _table(n, 'T')
    old = ('WAKE L RF\n' + ''.join('  %r %r\n' % p for p in tl) + 'ENDWAKE;\n'
           + 'WAKE T RF\n' + ''.join(' %r %r  %r %r\n' % (tt[j] + tt[j + 1]) for j in range(0, n, 2)) + 'ENDWAKE;\n')
    new = ('WakeFunction[Longitudinal, "RF"] = {' + ','.join('{%r,%r}' % p for p in tl) + '};\n'
           + 'WakeFunction[Transverse, "RF"] = {' + ','.join('{%r,%r}' % p for p in tt) + '};\n')
    r_old, out = _rows(_src(old))
    r_new, _ = _rows(_src(new))
    assert '?' not in out
    _same(r_old, r_new, 1e-12)
    # 6 点の表 (C3_REF) とは別の結果 = 400 点の表が実際に使われている
    assert any(abs(a - b) > 1e-6 * max(abs(b), 1e-3) for r, q in zip(r_old, C3_REF) for a, b in zip(r, q))


BASE2 = (C3_SRC.replace('CAVI RF = (VOLT=0 HARM=100 PHI=0.2 L=0.4);',
                        'CAVI RF1 = (VOLT=0 HARM=100 PHI=0.2 L=0.4) RF2 = (VOLT=0 HARM=100 PHI=0.2 L=0.2);')
         .replace('LINE BL = (ST QF DR RF DR QD DR);', 'LINE BL = (ST QF DR RF1 DR RF2 QD);'))


def test_wildcard_registers_all_cavities(v2_mode):
    assert 'RF1' in BASE2 and 'RF2 ' in BASE2
    old = MULTI.replace(' RF', ' RF*')
    new = ''.join('WakeFunction[Longitudinal, "%s"] = {{-0.1, 2e9}, {-0.02, 1e9}, {0, 5e8}, {0.02, -1e9}, {0.1, -2e9}, {0.2, -1e9}};\n'
                  'WakeFunction[Transverse, "%s"] = {{0, 1e12}, {0.01, 5e11}, {0.05, 1e11}, {0.2, 0}, {0.3, 0}};\n' % (c, c)
                  for c in ('RF1', 'RF2'))
    one = new.split('WakeFunction[Longitudinal, "RF2"]')[0]
    r_old, out = _rows(_src(old, BASE2))
    r_new, _ = _rows(_src(new, BASE2))
    r_one, _ = _rows(_src(one, BASE2))
    assert '?' not in out
    _same(r_old, r_new, 1e-12)
    assert any(abs(a - b) > 1e-6 * max(abs(b), 1e-3) for r, q in zip(r_old, r_one) for a, b in zip(r, q))     # RF2 にも入っている


def test_errors_are_explained(v2_mode):
    # ENDWAKE が無い: 次の ; までを表として読み、数でない語で止まる
    _, out = _rows(_src('WAKE L RF\n 0 1e9\n 0.1 0\nn0 = 3;\n'))
    assert '?Missing ENDWAKE.' in out and 'n0' in out
    # 表の途中の ;
    _, out = _rows(_src('WAKE L RF\n 0 1e9;\n'))
    assert '?Missing ENDWAKE.' in out and '表の途中に ; は書けない' in out
    # 値が奇数個
    _, out = _rows(_src('WAKE L RF\n 0 1e9\n 0.1\nENDWAKE;\n'))
    assert '?Missing value of wake.' in out and '3 個' in out
    # CAVI でない要素・無い要素・空の表
    _, out = _rows(_src('WAKE T QF\n 0 1e9\n 0.1 0\nENDWAKE;\n'))
    assert 'QF は CAVI ではない' in out
    _, out = _rows(_src('WAKE T NOSUCH\n 0 1e9\n 0.1 0\nENDWAKE;\n'))
    assert '?Undefined element NOSUCH' in out
    _, out = _rows(_src('WAKE T RF\nENDWAKE;\n'))
    assert '?No wake data.' in out
    # エラーの後も台本は続き、wake なしの追跡が出る
    rows, _ = _rows(_src('WAKE T QF\n 0 1e9\n 0.1 0\nENDWAKE;\n'))
    assert len(rows) == 6


def test_statement_after_endwake_on_same_line(v2_mode):
    rows, out = _rows(_src(MULTI.replace('ENDWAKE\n', 'ENDWAKE; Print["AFTER ", 7];\n')))
    assert 'AFTER 7' in out
    _same(rows, _ref_same_mode(), 2e-9)


def test_k64_mode_reads_line_by_line(k64_mode):
    from pysad.sadscript import Interp
    I = Interp(out=io.StringIO())
    src = 'WAKE L RF\n 0 1e9\n 0.1 0\nENDWAKE;\nx = 1;\n'
    assert len(I.split_statements(src)) == 5          # K64: 行ごと (従来どおり)
    from pysad import compat
    compat.set_mode('v2')
    try:
        st = I.split_statements(src)
    finally:
        compat.set_mode('k64')
    assert len(st) == 2 and st[0].split()[-1] == 'ENDWAKE'
    # WAKE 以外の文、ENDWAKE の無い 1 行の WAKE (WAKECOMMAND 用) はどちらのモードでも同じ
    assert len(I.split_statements('WAKE\nx = 1\ny = 2\n')) == 3
