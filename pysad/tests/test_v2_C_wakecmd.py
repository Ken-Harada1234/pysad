# -*- coding: utf-8 -*-
"""v2 の修正 B66: 旧式の WAKE コマンド (tfwake.f の構文)。K64 は tfwake.f を呼ばず WAKECOMMAND[] (配布に無い) を評価して失敗する。
v2 は `WAKE L|T 要素 s1 w1 ... ENDWAKE` を WakeFunction の表として登録する。根拠: 同じ表を WakeFunction で与えた追跡 (gs と一致済み) と同じ結果。"""
import io

from pysad.tests.test_wake import C3_SRC, C3_REF

OLD = '''WAKE LONGITUDINAL RF -0.1 2e9 -0.02 1e9 0 5e8 0.02 -1e9 0.1 -2e9 0.2 -1e9 ENDWAKE;
WAKE T RF 0 1e12 0.01 5e11 0.05 1e11 0.2 0 0.3 0 ENDWAKE;
'''


def _src():
    lines = C3_SRC.splitlines()
    out = [l for l in lines if not l.startswith('WakeFunction[')]
    k = out.index('PBUNCH = 1E12; NBUNCH = 1;')
    return '\n'.join(out[:k + 1]) + '\n' + OLD + '\n'.join(out[k + 1:]) + '\n'


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


def test_b66_old_wake_command_equals_wakefunction(v2_mode):
    rows, _ = _rows(_src())
    assert len(rows) == 6
    for r, q in zip(rows, C3_REF):
        for a, b in zip(r, q):
            assert abs(a - b) <= 2e-9 * max(abs(b), 1e-3)


def test_b66_missing_endwake_is_reported(v2_mode):
    _, out = _rows(_src().replace('0.2 0 0.3 0 ENDWAKE;', '0.2 0 0.3 0;'))
    assert '?Missing ENDWAKE.' in out


def test_b66_k64_mode_has_no_wake(k64_mode):
    rows, _ = _rows(_src())
    assert any(abs(a - b) > 1e-6 * max(abs(b), 1e-3) for r, q in zip(rows, C3_REF) for a, b in zip(r, q))
