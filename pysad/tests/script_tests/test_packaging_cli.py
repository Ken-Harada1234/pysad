# -*- coding: utf-8 -*-
"""配布の内容: pysad.cli (pysad-run・pysad-share・pysad-warmup)、setup.py の data_files (docs/・example/ を share/pysad/ に入れる)、版番号の一貫性。"""
import importlib.util
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def _load_setup():
    import setuptools
    cap = {}
    orig = setuptools.setup
    setuptools.setup = lambda **k: cap.update(k)
    try:
        spec = importlib.util.spec_from_file_location('pysad_setup', os.path.join(ROOT, 'setup.py'))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
    finally:
        setuptools.setup = orig
    return cap


def test_share_dir_in_clone_is_repo_root():
    from pysad.cli import share_dir
    d = share_dir()
    assert d is not None
    assert os.path.isdir(os.path.join(d, 'docs')) and os.path.isdir(os.path.join(d, 'example'))


def test_share_main_prints_dir(capsys):
    from pysad.cli import share_main
    assert share_main([]) == 0
    assert capsys.readouterr().out.strip()


def test_run_main_help(capsys):
    from pysad.cli import run_main
    assert run_main(['--help']) == 0
    assert 'pysad-run' in capsys.readouterr().out
    assert run_main([]) == 2


def test_data_files_cover_docs_and_example():
    cap = _load_setup()
    df = cap['data_files']
    targets = [t for t, _ in df]
    assert 'share/pysad/docs' in targets
    assert any(t.startswith('share/pysad/example/14_nonlinear_dynamics') for t in targets)
    allfiles = [f for _, fl in df for f in fl]
    assert 'docs/01_overview.md' in allfiles and 'run_sad.py' in allfiles
    assert not any('__pycache__' in f or f.endswith('.pyc') for f in allfiles)
    assert not any('/_deletable/' in f for f in allfiles)


def test_version_consistent():
    import pysad
    toml = open(os.path.join(ROOT, 'pyproject.toml'), encoding='utf-8').read()
    v = re.search(r'^version = "([^"]+)"', toml, re.M).group(1)
    assert v == pysad.__version__
    assert 'pysad-run = "pysad.cli:run_main"' in toml and 'pysad-share = "pysad.cli:share_main"' in toml
    assert 'pysad-warmup = "pysad.cli:warmup_main"' in toml


# ---- pysad-warmup (numba のカーネルの事前コンパイル)。ここでは実際の変換はしない (PYSAD_SLOW=1 のときだけ全部走らせる) ----

def test_warmup_main_help_and_args(capsys):
    from pysad.cli import warmup_main
    assert warmup_main(['--help']) == 0
    assert 'pysad-warmup' in capsys.readouterr().out
    assert warmup_main(['--nonsense']) == 2
    assert warmup_main(['--check', '--clear']) == 2
    assert 'pysad-warmup = "pysad.cli:warmup_main"' in open(os.path.join(ROOT, 'pyproject.toml'), encoding='utf-8').read()


def test_warmup_kernel_inventory():
    from pysad import warmup
    mods = warmup.kernel_modules()
    assert {'elkern', 'track6n', 'turnkern'} <= set(mods) and 'nbcache' not in mods
    ds = warmup.dispatchers()
    assert len(ds) > 100
    assert all(warmup._cache_file(d) is not None for _, _, d in ds)     # すべて cache=True
    assert warmup.cache_dirs(ds)


def test_warmup_check_does_not_modify_cache(capsys):
    import io
    from pysad import warmup
    dirs = warmup.cache_dirs()
    before = {f: os.stat(f).st_mtime_ns for f in warmup._cache_files(dirs)}
    out = io.StringIO()
    assert warmup.do_check(out=out) in (0, 1)
    txt = out.getvalue()
    assert 'numba のキャッシュ' in txt and ('判定' in txt or '記録がありません' in txt)
    assert warmup.main(['--check']) in (0, 1)
    assert {f: os.stat(f).st_mtime_ns for f in warmup._cache_files(dirs)} == before


def test_warmup_check_reports_missing_record(tmp_path, monkeypatch, capsys):
    from pysad import warmup
    monkeypatch.setattr(warmup, 'cache_dirs', lambda ds=None: [str(tmp_path)])
    assert warmup.do_check() == 1
    assert 'pysad-warmup' in capsys.readouterr().out


def test_warmup_clear_removes_only_pysad_kernel_cache(tmp_path, monkeypatch, capsys):
    from pysad import warmup
    names = ['elkern._tquad-100.py312.nbi', 'elkern._tquad-100.py312.0.nbc', 'track6n.k_bend-5.py312.nbi',
             warmup.RECORD_NAME, 'nbdeps.json', 'other_package.f-1.py312.nbi', 'notes.txt']
    for n in names:
        (tmp_path / n).write_bytes(b'x' * 10)
    monkeypatch.setattr(warmup, 'cache_dirs', lambda ds=None: [str(tmp_path)])
    assert warmup.main(['--clear']) == 0
    left = sorted(os.listdir(tmp_path))
    assert left == ['notes.txt', 'other_package.f-1.py312.nbi']
    assert '5 ファイル' in capsys.readouterr().out


@pytest.mark.skipif(os.environ.get('PYSAD_SLOW') != '1', reason='全カーネルの変換は長い (数十秒〜数分)。PYSAD_SLOW=1 で走らせる')
def test_warmup_full_run_then_check(capsys):
    from pysad import warmup
    assert warmup.do_warmup(verbose=False) == 0
    assert warmup.do_check() == 0
