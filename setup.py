# -*- coding: utf-8 -*-
"""pyproject.toml の設定に加えて、文書 (docs/)、例 (example/)、物理ガイド (physics_guide/) を `<prefix>/share/pysad/` に入れる (data_files)。
pip install すると、パッケージ本体 (import pysad) のほかに、文書と例の全ファイル (HTML・図・出力を含む) が入る。
場所は `pysad-share` が表示する。"""
import os

from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_DIRS = {'__pycache__', '.pytest_cache', '_deletable'}
SKIP_SUFFIX = ('.pyc', '.pyo')


def _collect(top):
    out = []
    root = os.path.join(HERE, top)
    for dp, dn, fn in os.walk(root):
        dn[:] = sorted(d for d in dn if d not in SKIP_DIRS)
        files = sorted(os.path.join(os.path.relpath(dp, HERE), f).replace(os.sep, '/')
                       for f in fn if not f.endswith(SKIP_SUFFIX))
        if files:
            out.append((os.path.join('share', 'pysad', os.path.relpath(dp, HERE)).replace(os.sep, '/'), files))
    return out


data_files = []
for top in ('docs', 'example', 'physics_guide'):
    if os.path.isdir(os.path.join(HERE, top)):
        data_files += _collect(top)
for f in ('run_sad.py',):
    if os.path.isfile(os.path.join(HERE, f)):
        data_files.append(('share/pysad', [f]))
# 例の台本が相対パス pysad/tests/kek-ls_v03_12.lat で読むラティス (pip で入れたとき、share/pysad をカレントにして例が動くように、同じ相対パスに置く)
if os.path.isfile(os.path.join(HERE, 'pysad', 'tests', 'kek-ls_v03_12.lat')):
    data_files.append(('share/pysad/pysad/tests', ['pysad/tests/kek-ls_v03_12.lat']))

setup(data_files=data_files)
