# -*- coding: utf-8 -*-
"""pysad に同梱した SAD のコピー (pysad/sad_copy) の試験。

pysad/sad_copy は KEK の SAD のソースからのコピーで (pysad/sad_copy/README.md)、環境変数 PYSAD_SAD_PACKAGES が無くても
SAD の Packages (CorrectOrbit・Steer・Monitor など) が使えるようにするためのもの。GUI (Tk) の接点は、画面を開かないよう、
環境変数 PYSAD_SAD_PACKAGES か PYSAD_GUI=1 が無いと付かない。
- 各ファイルが MANIFEST.sha256 と同じ (コピー後に変えていない)
- 環境変数が無いとき packages_dir() が同梱のコピーを返す。環境変数があればそちらが優先
- 環境変数が無い子プロセスで、Packages を使う台本 (CorrectOrbit・Steer・Monitor) が走り、元の SAD の Packages を使った結果と一致する (SAD が手元にあるとき)
- PYSAD_GUI=1 だけで、同梱のコピーから GUI の接点が動く (画面があるとき)
- LICENSE と、コピーであることを書いた README がある
"""
import hashlib
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
COPY = os.path.join(ROOT, 'pysad', 'sad_copy')
EXAMPLE = os.path.join(ROOT, 'example', '13_packages', 'fn_CorrectOrbit_Monitor_Steer.sad')
BS = chr(92)


def _env(without=True, **kw):
    e = dict(os.environ)
    if without:
        e.pop('PYSAD_SAD_PACKAGES', None)
    e.pop('PYSAD_GUI', None)
    e.update(kw)
    e['PYTHONPATH'] = ROOT + os.pathsep + e.get('PYTHONPATH', '')
    e['PYTHONIOENCODING'] = 'utf-8'
    return e


def _norm(p):
    return p.replace(BS, '/').rstrip('/')


def test_copy_files_match_the_manifest():
    lines = open(os.path.join(COPY, 'MANIFEST.sha256'), encoding='utf-8').read().split('\n')
    entries = [ln.split('  ', 1) for ln in lines if ln.strip()]
    assert len(entries) > 100
    for h, rel in entries:
        p = os.path.join(COPY, rel)
        assert os.path.isfile(p), rel
        assert hashlib.sha256(open(p, 'rb').read()).hexdigest() == h, rel
    listed = {rel for _, rel in entries}
    for sub in ('Packages', 'KBFrame'):
        for f in os.listdir(os.path.join(COPY, sub)):
            assert f'{sub}/{f}' in listed, f'{sub}/{f} が MANIFEST に無い'


def test_copy_has_license_and_provenance_readme():
    lic = open(os.path.join(COPY, 'LICENSE'), encoding='utf-8').read()
    assert 'KEK Accelerator Laboratory' in lic and 'Redistribution and use' in lic
    rd = open(os.path.join(COPY, 'README.md'), encoding='utf-8').read()
    assert 'SAD' in rd and 'コピー' in rd and 'pysad が書いたものではありません' in rd


def test_copy_does_not_include_machine_specific_files():
    names = os.listdir(os.path.join(COPY, 'Packages'))
    assert 'init.n' in names and 'correction.n' in names and 'Tkinter.n' in names
    assert 'init.local.n' not in names and not any(n.endswith('.bak') for n in names)


def test_packages_dir_prefers_env_then_bundled(monkeypatch, tmp_path):
    """pkgload (CorrectOrbit など) は、環境変数 > 同梱のコピー。tclbridge (GUI) は、環境変数 > PYSAD_GUI=1 のときの同梱のコピー > なし"""
    from pysad import pkgload, tclbridge
    monkeypatch.delenv('PYSAD_SAD_PACKAGES', raising=False)
    monkeypatch.delenv('PYSAD_GUI', raising=False)
    b = pkgload.packages_dir()
    assert b and _norm(b).endswith('pysad/sad_copy/Packages')
    assert tclbridge.packages_dir() is None                       # 画面を開かない (従来どおり)
    monkeypatch.setenv('PYSAD_GUI', '1')
    assert tclbridge.packages_dir() == b
    monkeypatch.delenv('PYSAD_GUI')
    monkeypatch.setenv('PYSAD_SAD_PACKAGES', str(tmp_path))
    assert _norm(pkgload.packages_dir()).endswith(tmp_path.name)
    assert _norm(tclbridge.packages_dir()).endswith(tmp_path.name)
    monkeypatch.setenv('PYSAD_SAD_PACKAGES', str(tmp_path / 'nonexistent'))      # 無効な指定は無視して同梱のコピー
    assert pkgload.packages_dir() == b


_RUN = '''import io, sys
from pysad.ffscmd import run_sad_string
src = open(sys.argv[1], encoding='utf-8').read()
buf = io.StringIO()
run_sad_string(src, out=buf)
sys.stdout.write(buf.getvalue())
'''


def _run(env):
    r = subprocess.run([sys.executable, '-c', _RUN, EXAMPLE], env=env, capture_output=True, text=True, encoding='utf-8', cwd=ROOT, timeout=900)
    assert r.returncode == 0, r.stderr[-2000:]
    return r.stdout


def test_correctorbit_script_runs_without_the_env_var():
    out = _run(_env(True))
    assert 'Matched.' in out and 'NX' in out and 'R:rmx1=' in out


def test_bundled_copy_gives_the_same_result_as_the_original_sad_packages():
    orig = os.environ.get('PYSAD_SAD_PACKAGES_ORIGINAL') or os.path.join(ROOT, '..', 'SAD-master', 'Packages')
    if not os.path.isfile(os.path.join(orig, 'correction.n')):
        pytest.skip('元の SAD の Packages が手元に無い (環境変数 PYSAD_SAD_PACKAGES_ORIGINAL で指定できる)')
    assert _run(_env(True)) == _run(_env(True, PYSAD_SAD_PACKAGES=orig))


_GUI = '''import io, os
os.environ['PYSAD_TK_NOLOOP'] = '1'
import tkinter
r = tkinter.Tk()
r.destroy()
from pysad.sadscript import Interp
I = Interp()
I.out = io.StringIO()
I.run('v = Tk$FullVersion; fs = $FontSystem;', 'script')
print(str(I.get('v'))[:3], I.get('fs'))
br = getattr(I, 'tclbridge', None)
print(br is not None and br.packages)
'''


def test_gui_bridge_works_from_the_bundled_copy_with_pysad_gui():
    """環境変数 PYSAD_SAD_PACKAGES が無く PYSAD_GUI=1 だけで、SAD の Tkinter.n・KBFrame.n (同梱のコピー) が動く (画面があるとき)"""
    r = subprocess.run([sys.executable, '-c', _GUI], env=_env(True, PYSAD_GUI='1'), capture_output=True, text=True, encoding='utf-8', cwd=ROOT, timeout=600)
    if r.returncode != 0 and 'TclError' in r.stderr:
        pytest.skip('tkinter が使えない')
    assert r.returncode == 0, r.stderr[-2000:]
    first, second = r.stdout.strip().split('\n')[:2]
    assert first.split()[0] in ('8.6', '9.0') and first.split()[1] == 'X11'
    assert _norm(second).endswith('pysad/sad_copy/Packages')


def test_no_gui_without_env_or_pysad_gui():
    code = "from pysad.sadscript import Interp\nI = Interp()\nprint(getattr(I, 'tclbridge', None) is None)\n"
    r = subprocess.run([sys.executable, '-c', code], env=_env(True), capture_output=True, text=True, encoding='utf-8', cwd=ROOT, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.strip() == 'True'


def test_nothing_outside_the_package_is_referenced():
    """別のフォルダから起動しても、Packages は pysad の中のコピーを指す"""
    code = "import os\nfrom pysad import pkgload\nb = pkgload.packages_dir()\nassert os.path.isdir(b)\nprint(b)\n"
    r = subprocess.run([sys.executable, '-c', code], env=_env(True), capture_output=True, text=True, encoding='utf-8', cwd=os.path.expanduser('~'), timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert 'pysad/sad_copy' in r.stdout.replace(BS, '/')
