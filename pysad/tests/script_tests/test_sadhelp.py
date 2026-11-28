# -*- coding: utf-8 -*-
"""pysad ヘルプ (pysad/sadhelp.py・helpdata/sadhelp_index.json・コマンド `pysad-help`・`python -m pysad`) と、説明文・import の副作用の試験。

SADHelp の全文の取り出しは、SADHelp.html と同じ形式の小さな偽のファイルで試す (本物の SADHelp.html は使わない)。
"""
import inspect
import os
import re
import subprocess
import sys

import pytest

import pysad
from pysad import cli, sadhelp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
JA = re.compile(r'[぀-ヿ一-鿿]')

FAKE_HELP = '''<HEAD>
<TITLE>
</TITLE>
</HEAD>
<BODY>
<H3><A NAME=L1>
</A>
</H3>
<PRE>
</PRE>
<UL>
<LI>
<H3><A NAME=L2>
CALCULATE(CAL)</A>
</H3>
<PRE>Usage: (1) CALCULATE(CAL) [[NO]EXPAND]
calculates the optics &amp; matching functions.
</PRE>
<LI>
<H3><A NAME=L3>
elements</A>
</H3>
<PRE>An element.
</PRE>
<UL>
<LI>
<H3><A NAME=L4>
BEND</A>
</H3>
<PRE>A bending magnet (fake).
</PRE>
<UL>
<LI>
<H3><A NAME=L5>
AE1</A>
</H3>
<PRE>FAKE-AE1-UNDER-BEND
</PRE>
</UL>
<LI>
<H3><A NAME=L6>
MULT</A>
</H3>
<PRE>A multipole (fake).
</PRE>
<UL>
<LI>
<H3><A NAME=L7>
AE1</A>
</H3>
<PRE>FAKE-AE1-UNDER-MULT
</PRE>
</UL>
</UL>
</UL>
</BODY>
'''


# ---------------------------------------------------------------------- データと引き方
def test_index_is_bundled_and_has_the_main_kinds():
    idx = sadhelp.index()
    assert idx['sadhelp_version'] and len(idx['entries']) >= 500
    kinds = {e['kind'] for e in idx['entries']}
    assert {'ffs-command', 'function', 'special-variable', 'flag', 'element', 'element-keyword'} <= kinds
    n = {k: sum(1 for e in idx['entries'] if e['kind'] == k) for k in kinds}
    assert n['ffs-command'] == 56 and n['special-variable'] == 48 and n['flag'] == 60      # docs/34 の件数


def test_lookup_name_case_abbreviation_and_brackets():
    a = sadhelp.lookup('Twiss')
    assert len(a) == 1 and a[0]['kind'] == 'function' and a[0]['status'] == '実装'
    assert sadhelp.lookup('twiss') == a and sadhelp.lookup('Twiss[]') == a and sadhelp.lookup('`Twiss`') == a
    c = sadhelp.lookup('CAL')                      # 略形
    assert len(c) == 1 and c[0]['name'].startswith('CALCULATE') and sadhelp.lookup('calculate') == c
    assert sadhelp.lookup('NoSuchName123') == [] and sadhelp.lookup('') == []


def test_lookup_same_name_as_command_and_function_and_kind_filter():
    r = sadhelp.lookup('FIT')
    assert {e['kind'] for e in r} == {'ffs-command', 'function'}
    assert [e['kind'] for e in sadhelp.lookup('FIT', kind='function')] == ['function']


def test_lookup_element_keyword_with_parent():
    r = sadhelp.lookup('BEND.AE1')
    assert len(r) == 1 and r[0]['parent'] == 'BEND' and r[0]['kind'] == 'element-keyword'
    assert sadhelp.lookup('BEND AE1') == r and sadhelp.lookup('MULT.AE1')[0]['parent'] == 'MULT'
    assert len(sadhelp.lookup('AE1')) >= 2           # 親を付けなければ、AE1 を持つ要素すべて
    b = sadhelp.lookup('BEND', kind='element')[0]
    assert 'AE1' in b['keywords'] and b['pysad_type'] == 'あり'


def test_status_note_and_docs_link_come_from_docs35():
    f = sadhelp.lookup('Fit', kind='function')[0]
    assert f['status'] == '実装' and 'fitroot' in f['note'] and f['docs'] == ['06_sadscript_library.md']
    assert os.path.isfile(os.path.join(ROOT, 'docs', f['docs'][0])) or not os.path.isdir(os.path.join(ROOT, 'docs'))
    assert any(e.get('status') == '未対応' for e in sadhelp.entries())


def test_pysad_only_functions_have_english_and_japanese_lines():
    for n in ('FrequencyMap', 'MomentumAcceptance', 'TouschekLifetimeLMA', 'ResonanceDrivingTerms', 'AmplitudeTuneShift', 'ChromaticityTerms', 'NormalFormTune', 'NAFF'):
        e = sadhelp.lookup(n)[0]
        assert e['status'] == 'pysad のみ' and e['summary'] and not JA.search(e['summary']) and JA.search(e['summary_ja'])


def test_naff_is_a_method_entry_linked_to_the_primer():
    e = sadhelp.lookup('NAFF')[0]
    assert e['kind'] == 'method' and any('24_adts_rdt_primer' in d for d in e['docs'])


def test_search_finds_names_summaries_and_notes():
    names = [e['name'] for e in sadhelp.search('chromaticity')]
    assert 'ChromaticityTerms' in names and any(n.startswith('CHROMATICITY') for n in names)
    assert sadhelp.search('') == [] and sadhelp.search('zzzNoSuchKeyword') == []
    assert all(e['kind'] == 'flag' for e in sadhelp.search('on', kind='flag'))


# ---------------------------------------------------------------------- 表示
def test_text_for_a_name_shows_summary_status_source():
    t = sadhelp.text('Twiss')
    assert t.startswith('Twiss  [function') and 'SADHelp :' in t and 'pysad   : 実装' in t and 'SAD 側  :' in t


def test_text_shows_docs_link_for_fit():
    t = sadhelp.text('Fit', kind='function')
    assert 'https://ken-harada1234.github.io/pysad/docs/06_sadscript_library.html' in t


def test_text_for_unknown_name_suggests_and_search_lists():
    t = sadhelp.text('Twis')                       # 一部だけ
    assert '見つかりませんでした' in t and 'Twiss' in t
    assert '見つかりませんでした' in sadhelp.text('zzzNoSuch')
    s = sadhelp.text(search_for='tune', limit=3)
    assert '件' in s and '詳しくは' in s


def test_overview_and_sad_help_prints(capsys):
    ov = sadhelp.overview()
    assert 'pysad ヘルプ' in ov and str(len(sadhelp.entries())) in ov
    pysad.sad_help('CAL')
    out = capsys.readouterr().out
    assert 'CALCULATE(CAL)' in out and out.endswith('\n')


# ---------------------------------------------------------------------- 全文 (同梱の HelpMessages.n)
MINI_HELP = r"""!comment line
Help$Messages=(Print["Help"];{
{"",
 "top text",
  "",
 Sort[{
{"FOO(F)",
   "\
Usage: FOO [x]\n\
second line with \"quote\"",
  "BAR Baz",
  {}},
{"Bar",
   "bar text",
  "",
   Sort[{
     {"kw",
       "keyword of bar",
       "",
       {}}
   }]},
{"Table",
   "generated",
  "",
  Flatten[Table[{i}, {i, 3}]]}
}]}
});
"""


@pytest.fixture
def no_user_help(monkeypatch):
    monkeypatch.delenv('SADHELP_HTML', raising=False)
    monkeypatch.setattr(sadhelp, 'find_sadhelp_html', lambda: None)


def test_full_text_comes_from_the_bundled_helpmessages_without_any_user_file(no_user_help):
    """同梱の SAD のコピーの HelpMessages.n から、全文が出る (手元の SADHelp.html は要らない)"""
    bp = sadhelp.bundled_help_path()
    assert bp and bp.replace('\\', '/').endswith('sad_copy/Packages/HelpMessages.n')
    t = sadhelp.text('CAL', full=True)
    assert '全文 (SAD の HelpMessages.n (同梱の SAD のコピー)):' in t and 'matching-function1' in t
    assert 'Twiss[optical-function, component] returns' in sadhelp.text('Twiss', full=True)
    assert 'DP0' in sadhelp.text('DP0', full=True) and '全文 (SAD の' in sadhelp.text('RFSW', full=True)
    assert '関連 (See also):' in t                                      # 3 つ目の要素 (関連する名前)


def test_bundled_helpmessages_cover_commands_flags_variables_and_most_functions(no_user_help):
    miss = {}
    for k in ('ffs-command', 'flag', 'special-variable', 'function'):
        es = [e for e in sadhelp.entries() if e['kind'] == k]
        miss[k] = [e['name'] for e in es if sadhelp.full_text(e) is None]
        assert len(es) - len(miss[k]) >= {'ffs-command': 55, 'flag': 60, 'special-variable': 48, 'function': 94}[k], (k, miss[k])
    assert set(miss['ffs-command']) <= {'SHOW'}


def test_full_text_of_elements_is_partly_in_the_bundled_file(no_user_help):
    """要素のキーワードの全文は、同梱の HelpMessages.n にあるもの (BEND.AE1) と、実行時に作る部分で取れないもの (CAVI.DX) がある。
    取れないときは、そう表示し、SADHELP_HTML と CERN の URL を案内する"""
    ok = sadhelp.text('BEND.AE1', full=True)
    assert '全文 (SAD の HelpMessages.n' in ok and 'face angle at the entrance' in ok
    t = sadhelp.text('CAVI.DX', full=True)
    assert '見つかりませんでした' in t and 'SADHELP_HTML' in t and 'hep-project-sad.web.cern.ch' in t


def test_parser_of_helpmessages_format(tmp_path, monkeypatch, no_user_help):
    """`{名前, 本文, 関連, {下位項目}}` の入れ子: 行末の \\ は続き行、\\n・\\" は文字、Sort[..] は中身、実行時に作る部分 (Table など) は項目にしない"""
    (tmp_path / 'HelpMessages.n').write_text(MINI_HELP, encoding='utf-8')
    from pysad import pkgload
    monkeypatch.setattr(pkgload, 'packages_dir', lambda: str(tmp_path))
    monkeypatch.setattr(sadhelp, '_FULL', {})
    foo = {'name': 'FOO(F)', 'keys': ['FOO', 'F']}
    assert sadhelp.full_text(foo) == 'Usage: FOO [x]\nsecond line with "quote"\n\n関連 (See also): BAR Baz'
    assert sadhelp.full_text({'name': 'bar', 'keys': ['bar']}) == 'bar text'
    assert sadhelp.full_text({'name': 'kw', 'keys': ['kw'], 'parent': 'Bar'}) == 'keyword of bar'
    assert sadhelp.full_text({'name': 'kw', 'keys': ['kw'], 'parent': 'Nope'}) is None
    assert sadhelp.full_text({'name': 'Table', 'keys': ['Table']}) in (None, 'generated')


def test_parser_does_not_execute_code_in_helpmessages(tmp_path, monkeypatch, no_user_help):
    """HelpMessages.n は構文だけを解析する (評価しない): ファイルを開く・書くコードがあっても実行されない"""
    marker = tmp_path / 'written.txt'
    src = MINI_HELP.replace('Print["Help"]', 'Write[OpenWrite["%s"], "x"]' % str(marker).replace('\\', '/'))
    (tmp_path / 'HelpMessages.n').write_text(src, encoding='utf-8')
    from pysad import pkgload
    monkeypatch.setattr(pkgload, 'packages_dir', lambda: str(tmp_path))
    monkeypatch.setattr(sadhelp, '_FULL', {})
    assert sadhelp.full_text({'name': 'Bar', 'keys': ['Bar']}) == 'bar text'
    assert not marker.exists()


# ---------------------------------------------------------------------- 手元の SADHelp.html (全文。同梱より優先)
def test_full_text_from_a_user_file_in_the_sadhelp_format_has_priority(tmp_path, monkeypatch):
    p = tmp_path / 'SADHelp.html'
    p.write_text(FAKE_HELP, encoding='utf-8')
    monkeypatch.setenv('SADHELP_HTML', str(p))
    monkeypatch.setattr(sadhelp, '_FULL', {})
    assert sadhelp.find_sadhelp_html() == str(p)
    t = sadhelp.text('CAL', full=True)
    assert '全文 (あなたの SADHelp.html)' in t and 'calculates the optics & matching functions.' in t        # &amp; は戻す
    a = sadhelp.text('BEND.AE1', full=True)
    b = sadhelp.text('MULT.AE1', full=True)
    assert 'FAKE-AE1-UNDER-BEND' in a and 'FAKE-AE1-UNDER-MULT' not in a        # 親の要素で取り分ける
    assert 'FAKE-AE1-UNDER-MULT' in b and 'FAKE-AE1-UNDER-BEND' not in b
    t2 = sadhelp.text('Twiss', full=True)                                       # 手元のファイルに無い項目は、同梱の HelpMessages.n から
    assert '全文 (SAD の HelpMessages.n' in t2


def test_full_text_is_not_stored_in_the_index():
    """索引 (JSON) には全文を入れない: 説明は最初の行だけ (どの項目も改行を含まず、長くない)。全文は同梱の HelpMessages.n を実行時に読む"""
    for e in sadhelp.entries():
        assert '\n' not in e.get('summary', '') and len(e.get('summary', '')) <= 400
    assert os.path.getsize(os.path.join(ROOT, 'pysad', 'helpdata', 'sadhelp_index.json')) < 800_000


# ---------------------------------------------------------------------- 同梱の .n の関数・Tkinter
def test_n_package_functions_have_signature_and_file():
    r = sadhelp.lookup('CorrectOrbit')
    assert len(r) == 1 and r[0]['kind'] == 'package-function' and r[0]['n_files'] == ['correction.n']
    assert r[0]['signature'][0].startswith('CorrectOrbit[')
    t = sadhelp.text('CorrectOrbit')
    assert '.n      : Packages/correction.n' in t and '署名    : CorrectOrbit[' in t
    pk = [e for e in sadhelp.entries() if e['kind'] == 'package-function']
    assert len(pk) > 300 and all(e['n_files'] and e['signature'] for e in pk)
    assert not any('$' in e['name'] or e['name'][0].islower() for e in pk)       # 内部の名前は入れない


def test_n_functions_already_in_sadhelp_get_the_file_and_signature():
    f = sadhelp.lookup('Fit', kind='function')[0]
    assert f['n_files'] == ['Fit.n'] and f['signature']
    assert 'n_files' not in sadhelp.lookup('FIT', kind='ffs-command')[0]       # FFS コマンドの FIT は別物


def test_tkinter_widgets_are_listed_with_their_tk_source():
    w = sadhelp.lookup('Frame')
    assert len(w) == 1 and w[0]['kind'] == 'tk-widget' and 'tkFrame.c' in w[0]['summary']
    assert any(e['kind'] == 'tk-widget' for e in sadhelp.search('Tk widget'))
    assert sadhelp.lookup('After')[0]['n_files'] == ['Tkinter.n']
    assert 'tcltk_gui' in ' '.join(sadhelp.lookup('After')[0]['docs'])


# ---------------------------------------------------------------------- データが docs と同じか・同梱の設定
def test_index_is_up_to_date_with_docs():
    tool = os.path.join(ROOT, 'tools', 'build_help_data.py')
    if not (os.path.isfile(tool) and os.path.isdir(os.path.join(ROOT, 'docs'))):
        pytest.skip('tools/ と docs/ がある clone のときだけ')
    r = subprocess.run([sys.executable, tool, '--check'], capture_output=True, text=True, encoding='utf-8')
    assert r.returncode == 0 and 'same' in r.stdout, r.stdout + r.stderr


def test_helpdata_is_in_wheel_and_sdist_configuration_and_command_is_registered():
    import tomllib
    with open(os.path.join(ROOT, 'pyproject.toml'), 'rb') as f:
        cfg = tomllib.load(f)
    assert any('helpdata' in p for p in cfg['tool']['setuptools']['package-data']['pysad'])
    assert cfg['project']['scripts']['pysad-help'] == 'pysad.cli:help_main'
    man = [ln.strip() for ln in open(os.path.join(ROOT, 'MANIFEST.in'), encoding='utf-8').read().splitlines()]
    assert 'graft pysad/helpdata' in man
    assert os.path.isfile(os.path.join(ROOT, 'pysad', 'helpdata', 'sadhelp_index.json'))


# ---------------------------------------------------------------------- コマンド
def test_pysad_help_command(capsys):
    assert cli.help_main(['Twiss']) == 0
    assert 'Twiss  [function' in capsys.readouterr().out
    assert cli.help_main(['--search', 'chromaticity']) == 0
    assert 'ChromaticityTerms' in capsys.readouterr().out
    assert cli.help_main([]) == 0
    assert 'pysad ヘルプ' in capsys.readouterr().out
    assert cli.help_main(['BEND', 'AE1', '--kind', 'element-keyword']) == 0
    assert 'BEND.AE1' in capsys.readouterr().out


def test_pysad_help_command_help_option_is_bilingual(capsys):
    with pytest.raises(SystemExit) as e:
        cli.help_main(['--help'])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert 'Show help for' in out and 'SAD の FFS コマンド' in out and 'SADHELP_HTML' in out


def test_python_dash_m_pysad_dispatch(capsys):
    assert cli.module_main([]) == 2
    out = capsys.readouterr().out
    assert 'help' in out and 'run' in out and 'share' in out and 'warmup' in out
    assert cli.module_main(['help', 'CAL']) == 0
    assert 'CALCULATE(CAL)' in capsys.readouterr().out
    assert cli.module_main(['bogus']) == 2
    assert '不明なコマンド' in capsys.readouterr().out
    assert cli.module_main(['--help']) == 0


def test_python_dash_m_pysad_runs_as_a_module():
    env = dict(os.environ, PYSAD_QUIET='1', PYTHONIOENCODING='utf-8')
    r = subprocess.run([sys.executable, '-m', 'pysad', 'help', 'CAL'], capture_output=True, text=True, encoding='utf-8', cwd=ROOT, env=env)
    assert r.returncode == 0 and 'CALCULATE(CAL)' in r.stdout, r.stdout + r.stderr


def test_help_options_of_all_commands_are_bilingual_and_do_not_run_anything(capsys):
    assert cli.share_main(['--help']) == 0           # 以前は --help を無視して実行していた
    out = capsys.readouterr().out
    assert 'usage' in out and '使い方' in out and 'pysad-share' in out and 'GitHub' in out
    assert cli.run_main(['--help']) == 0
    out = capsys.readouterr().out
    assert 'usage' in out and '使い方' in out
    from pysad import warmup
    assert warmup.main(['--help']) == 0
    out = capsys.readouterr().out
    assert 'usage' in out and '使い方' in out and '--check' in out


# ---------------------------------------------------------------------- import の副作用・説明文
def test_importing_diag_modules_prints_nothing():
    """旧い診断 (pysad/diag_fit.py・diag_trace.py) は、import しただけでは何も出さない (以前は import で診断が走って表を出した)"""
    env = dict(os.environ, PYSAD_QUIET='1', PYTHONIOENCODING='utf-8')
    r = subprocess.run([sys.executable, '-c', 'import pysad.diag_fit, pysad.diag_trace'], capture_output=True, text=True, encoding='utf-8', cwd=ROOT, env=env)
    assert r.returncode == 0 and r.stdout == '' and r.stderr == '', r.stdout[:200] + r.stderr[:200]
    from pysad import diag_fit, diag_trace
    assert callable(diag_fit.main) and callable(diag_trace.main)


@pytest.mark.parametrize('obj', [pysad.run, pysad.sad_help, sadhelp.sad_help, sadhelp.lookup, sadhelp.search, sadhelp.text, sadhelp.full_text,
                                 cli.help_main, cli.module_main, cli.run_main, cli.share_main, cli.warmup_main])
def test_new_public_entry_points_have_english_and_japanese_docstrings(obj):
    d = inspect.getdoc(obj) or ''
    en_lines = [ln for ln in d.split('\n') if re.search(r'[A-Za-z]{4,}', ln) and not JA.search(ln)]
    assert d and JA.search(d) and en_lines, obj.__name__


def test_session_methods_and_magics_have_docstrings():
    from pysad import notebook
    for name in ('add_lattice', 'clear_lattices', 'reset', 'run', 'value', 'help'):
        d = inspect.getdoc(getattr(notebook.Session, name))
        assert d and JA.search(d) and re.search(r'[A-Za-z]{4,}', d), name
    assert 'continue' in notebook.sad_magic.__doc__ and 'sad.value' in notebook.sad_magic.__doc__
    assert '--search' in notebook.sad_help_magic.__doc__ and 'SADHELP_HTML' in notebook.sad_help_magic.__doc__


def test_package_docstring_describes_entry_points_and_its_example_runs():
    """`help(pysad)` の説明文: 入口を英語と日本語で書き、Quick start の例は実際に動いて、書いてある出力になる (古い API の例が残っていた)"""
    import io
    d = pysad.__doc__
    assert JA.search(d) and 'Entry points' in d and 'pysad.run' in d and 'sad_help' in d and '%load_ext pysad.notebook' in d and 'PYSAD_COMPAT=k64' in d
    for cmd in ('pysad-run', 'pysad-help', 'pysad-share', 'pysad-warmup'):
        assert cmd in d
    m = re.search(r">>> pysad\.run\(r\'\'\'(.*?)\'\'\'\)", d, re.S)
    assert m, 'Quick start の例がありません'
    src = '\n'.join(ln[4:] if ln.startswith('... ') else ln for ln in m.group(1).split('\n'))
    out = io.StringIO()
    pysad.run(src, out=out, reset=True)
    expected = re.search(r'prints (NX=[0-9.]+)', d).group(1)
    nx = [ln for ln in out.getvalue().split('\n') if ln.startswith('NX=')]
    assert nx and nx[0] == expected, (nx, expected)
