# -*- coding: utf-8 -*-
"""配布の内容: pysad.cli (pysad-run・pysad-share)、setup.py の data_files (docs/・example/ を share/pysad/ に入れる)、版番号の一貫性。"""
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
