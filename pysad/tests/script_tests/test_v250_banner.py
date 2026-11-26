# -*- coding: utf-8 -*-
"""起動時の案内 (pysad/banner.py) の試験。

Python から SAD のインタプリタを呼んだとき、プロセスの中で最初の 1 回だけ、バージョン・元になっている SAD (K64 版)・ソースの場所を
標準エラー出力に出す。`import pysad` だけでは出さない。標準出力は変えない。`PYSAD_QUIET=1` で出さない。
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))


def _run(code, quiet=False):
    e = dict(os.environ)
    e.pop('PYSAD_QUIET', None)
    if quiet:
        e['PYSAD_QUIET'] = '1'
    e['PYTHONPATH'] = ROOT + os.pathsep + e.get('PYTHONPATH', '')
    e['PYTHONIOENCODING'] = 'utf-8'
    r = subprocess.run([sys.executable, '-c', code], env=e, capture_output=True, text=True, encoding='utf-8', cwd=ROOT, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout, r.stderr


def test_banner_text_has_version_origin_and_github():
    from pysad import __version__, banner
    t = banner.text()
    assert __version__ in t
    assert 'K64' in t and 'SAD' in t and 'Oide' in t
    assert 'https://github.com/Ken-Harada1234/pysad' in t
    assert 'https://ken-harada1234.github.io/pysad/' in t and 'not in the pip package' in t      # 例・文書・物理ガイド・試験は pip に入れない
    assert 2 <= len(t.strip().split('\n')) <= 6                  # 数行


def test_banner_is_not_shown_by_import_alone():
    out, err = _run('import pysad')
    assert 'pysad' not in err and err.strip() == '' and out == ''


def test_banner_is_shown_once_on_first_interpreter_call_to_stderr():
    out, err = _run('from pysad.sadscript import Interp\nI1 = Interp()\nI2 = Interp()\nprint("done")\n')
    assert out.strip() == 'done'                                  # 標準出力には出ない
    assert err.count('Source code: https://github.com/Ken-Harada1234/pysad') == 1
    assert 'K64' in err and 'Oide' in err


def test_banner_is_shown_for_run_sad_string_and_cli():
    out, err = _run('from pysad.ffscmd import run_sad_string\nimport io\nrun_sad_string("", out=io.StringIO())\n')
    assert 'github.com/Ken-Harada1234/pysad' in err


def test_banner_can_be_silenced_with_pysad_quiet():
    out, err = _run('from pysad.sadscript import Interp\nInterp()\n', quiet=True)
    assert err.strip() == ''
