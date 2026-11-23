"""Environment[名前, 値] が環境変数を設定する (Environment[名前] は値を返す)。

SAD のパッケージ (PYSAD_SAD_PACKAGES) を指定すると tclbridge.py が Environment を登録し直す。以前の登録は 1 引数形 (読み出し) だけで、
2 引数形を無視していた。example/02_optics_cod/fn_PUTMEASURE_GETMEASURE.sad が使う SADPOUT / SADPIN の設定が効かなくなる。
パッケージありの試験は、PYSAD_SAD_PACKAGES が指定されているときだけ走る (無ければ skip)。
"""
import io
import os

import pytest

PK = os.environ.get('PYSAD_SAD_PACKAGES')


def _run(body):
    from pysad.ffscmd import run_sad_string
    buf = io.StringIO()
    run_sad_string('FFS;\nPageWidth=1999;\n' + body + '\nExit;\n', out=buf, cwd='.')
    return [l.rstrip() for l in buf.getvalue().splitlines()]


def _check(name):
    os.environ.pop(name, None)
    try:
        r = _run('Environment["%s","hello"];\nPrint["V ",Environment["%s"]];' % (name, name))
        assert any(l.strip() == 'V hello' for l in r), r
        assert os.environ.get(name) == 'hello'
    finally:
        os.environ.pop(name, None)


def test_environment_two_args_sets_variable():
    _check('PYSAD_TEST_ENV_SET')


@pytest.mark.skipif(not PK or not os.path.isdir(PK), reason='環境変数 PYSAD_SAD_PACKAGES (SAD の Packages) が必要')
def test_environment_two_args_sets_variable_with_packages(monkeypatch):
    monkeypatch.setenv('PYSAD_SAD_PACKAGES', PK)
    _check('PYSAD_TEST_ENV_SET_PK')
