# -*- coding: utf-8 -*-
"""notebook からの SAD の実行 (`pysad.notebook`: `Session`・`pysad.run`・`%%sad` セルのマジック) の試験。

ラティスとコマンドの文字列を貼って実行する。Session の試験は IPython なしで走る。マジックの試験は IPython、本物のカーネルの試験は nbclient・ipykernel があるときだけ走る (無ければ skip)。
"""
import io
import os

import pytest

import pysad
from pysad import notebook
from pysad.notebook import Session

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
LAT = open(os.path.join(ROOT, 'example', '00_lattice', 'FODO_ring.sad'), encoding='utf-8').read()
# ! で始まるコメント行と、"..." の文字列を含むコマンド (貼り付けたファイルをそのまま)
CMD = '''FFS;
! comment line (SAD)
USE RING; CELL; CALC;
Print["NX=", Twiss["NX","$$$"]];
stop;
'''
NX = 42.07427561351298        # FODO_ring.sad の水平チューン (CALC)


def _nx(text):
    for line in text.split('\n'):
        if line.startswith('NX='):
            return float(line[3:])
    raise AssertionError('NX= の行がありません: ' + text[-300:])


def test_session_runs_stored_lattice_plus_command():
    s = Session()
    s.add_lattice(LAT)
    buf = io.StringIO()
    s.run(CMD, out=buf)
    assert _nx(buf.getvalue()) == pytest.approx(NX, abs=1e-9)


def test_lattice_with_the_same_name_is_replaced_not_appended():
    s = Session()
    s.add_lattice(LAT)
    s.add_lattice(LAT)                       # セルを再実行した場合
    assert len(s.lattices) == 1
    buf = io.StringIO()
    s.run(CMD, out=buf)
    assert _nx(buf.getvalue()) == pytest.approx(NX, abs=1e-9)    # 要素の二重定義にならない
    s.add_lattice(LAT, 'other')
    assert list(s.lattices) == ['main', 'other']


def test_run_without_stored_lattice_takes_everything_from_the_text():
    s = Session()
    buf = io.StringIO()
    s.run(LAT + '\n' + CMD, out=buf)
    assert _nx(buf.getvalue()) == pytest.approx(NX, abs=1e-9)


def test_continue_after_stop_runs_ffs_statements_in_the_same_state():
    s = Session()
    s.add_lattice(LAT)
    s.run(CMD, out=io.StringIO())            # `stop;` で終わっている
    buf = io.StringIO()
    s.run('Print["again=", Twiss["NX","$$$"]];\nQF = 0;\n', cont=True, out=buf)
    assert 'again=%s' % repr(NX)[:8] in buf.getvalue().replace(' ', '')
    assert s.value('Twiss["NX","$$$"]') == pytest.approx(NX, abs=1e-9)


def test_continue_changes_the_state_that_later_cells_see():
    s = Session()
    s.add_lattice(LAT)
    s.run(CMD, out=io.StringIO())
    before = s.value('LINE["K1","QF"]')
    s.run('LINE["K1","QF"] = 0.5;', cont=True, out=io.StringIO())
    assert s.value('LINE["K1","QF"]') == pytest.approx(0.5)
    assert before != pytest.approx(0.5)


def test_continue_without_a_previous_run_is_an_error():
    s = Session()
    with pytest.raises(RuntimeError, match='まだありません'):
        s.run('CALC;', cont=True, out=io.StringIO())
    with pytest.raises(RuntimeError):
        s.value('1+1')


def test_value_returns_python_numbers_and_lists_and_does_not_print():
    s = Session()
    s.add_lattice(LAT)
    buf = io.StringIO()
    s.run(CMD, out=buf)
    shown = buf.getvalue()
    assert s.value('2+3') == 5.0
    assert s.value('Range[3]') == [1.0, 2.0, 3.0]
    assert isinstance(s.value('Twiss[{"BX","BY"},"QF"]'), list)
    assert buf.getvalue() == shown            # value() は標準の出力先に何も書かない


def test_reset_forgets_lattices_and_state():
    s = Session()
    s.add_lattice(LAT)
    s.run(CMD, out=io.StringIO())
    s.reset()
    assert not s.lattices and s.interp is None and s.ffs is None


def test_pysad_run_with_a_string_and_lattice_argument():
    notebook.session.reset()
    buf = io.StringIO()
    r = pysad.run(CMD, lattice=LAT, out=buf)
    assert r is notebook.session and _nx(buf.getvalue()) == pytest.approx(NX, abs=1e-9)
    buf2 = io.StringIO()
    pysad.run(LAT + '\n' + CMD, reset=True, out=buf2)           # 1 つの文字列に全部
    assert _nx(buf2.getvalue()) == pytest.approx(NX, abs=1e-9)
    notebook.session.reset()


def test_parse_line_options():
    assert notebook._parse_line('') == ('run', 'main', False, False)
    assert notebook._parse_line('lattice') == ('lattice', 'main', False, False)
    assert notebook._parse_line('lattice ring2') == ('lattice', 'ring2', False, False)
    assert notebook._parse_line('--continue') == ('run', 'main', True, False)
    assert notebook._parse_line('-r') == ('run', 'main', False, True)
    for bad in ('--nope', 'foo', 'lattice a b', 'lattice --continue'):
        with pytest.raises(ValueError):
            notebook._parse_line(bad)


# ---------------------------------------------------------------------- IPython のマジック
def _shell():
    pytest.importorskip('IPython')
    from IPython.core.interactiveshell import InteractiveShell
    ip = InteractiveShell.instance()
    ip.run_line_magic('load_ext', 'pysad.notebook')
    notebook.session.reset()
    return ip


def test_cell_magic_lattice_then_commands(capsys):
    ip = _shell()
    ip.run_cell_magic('sad', 'lattice', LAT)
    ip.run_cell_magic('sad', '', CMD)
    out = capsys.readouterr().out
    assert 'ラティス' in out and '`%%sad` のセル' in out and _nx(out) == pytest.approx(NX, abs=1e-9)
    assert ip.user_ns['sad'] is notebook.session
    assert notebook.session.value('Twiss["NX","$$$"]') == pytest.approx(NX, abs=1e-9)


def test_cell_magic_continue_and_reset_and_errors(capsys):
    ip = _shell()
    from IPython.core.error import UsageError
    with pytest.raises(UsageError):
        ip.run_cell_magic('sad', '--continue', 'CALC;')           # まだ実行していない
    ip.run_cell_magic('sad', 'lattice', LAT)
    ip.run_cell_magic('sad', '', CMD)
    capsys.readouterr()
    ip.run_cell_magic('sad', '--continue', 'Print["more=", Twiss["NX","$$$"]];')
    assert 'more=' in capsys.readouterr().out
    ip.run_line_magic('sad', '--reset')
    assert not notebook.session.lattices
    with pytest.raises(UsageError):
        ip.run_cell_magic('sad', '--bogus', 'x')
    with pytest.raises(UsageError):
        ip.run_line_magic('sad', '')


# ---------------------------------------------------------------------- 本物の Jupyter のカーネル
def test_notebook_executes_in_a_real_kernel():
    """nbclient で notebook を実行する。セルの本文に `!` のコメント行と `%` を含め、実際の経路で壊れないことを確かめる"""
    nbformat = pytest.importorskip('nbformat')
    nbclient = pytest.importorskip('nbclient')
    pytest.importorskip('ipykernel')
    nb = nbformat.v4.new_notebook()
    body = CMD.replace('! comment line (SAD)', '! comment line with a percent % sign and !bang\n!! another comment')
    nb.cells = [
        nbformat.v4.new_code_cell('%load_ext pysad.notebook'),
        nbformat.v4.new_code_cell('%%sad lattice\n' + LAT),
        nbformat.v4.new_code_cell('%%sad\n' + body),
        nbformat.v4.new_code_cell('%%sad --continue\nPrint["cont=", Twiss["NX","$$$"]];'),
        nbformat.v4.new_code_cell('sad.value(\'Twiss["NX","$$$"]\')'),
    ]
    env_old = os.environ.get('PYSAD_QUIET')
    os.environ['PYSAD_QUIET'] = '1'
    try:
        client = nbclient.NotebookClient(nb, timeout=300, kernel_name='python3', resources={'metadata': {'path': ROOT}})
        try:
            client.execute()
        except Exception as e:       # カーネルの仕様が見つからないなど (環境の問題) は skip。セルの失敗は失敗にする
            if 'No such kernel' in str(e) or type(e).__name__ == 'NoSuchKernel':
                pytest.skip('python3 のカーネルがありません: %s' % e)
            raise
    finally:
        if env_old is None:
            os.environ.pop('PYSAD_QUIET', None)
        else:
            os.environ['PYSAD_QUIET'] = env_old
    def text(cell):
        return ''.join(o.get('text', '') for o in cell.outputs if o.output_type == 'stream' and o.get('name') == 'stdout')
    assert _nx(text(nb.cells[2])) == pytest.approx(NX, abs=1e-9)
    assert 'cont=' in text(nb.cells[3])
    res = nb.cells[4].outputs[0]
    assert res.output_type == 'execute_result' and float(res['data']['text/plain']) == pytest.approx(NX, abs=1e-9)
