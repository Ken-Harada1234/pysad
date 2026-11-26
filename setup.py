# -*- coding: utf-8 -*-
"""pyproject.toml の設定だけで、パッケージ本体 (import pysad、同梱の SAD のコピー pysad/sad_copy) を入れる。

文書 (docs/)・例 (example/)・物理ガイド (physics_guide/)・試験 (pysad/tests/) は、pip には入れない (GitHub と GitHub Pages にある。
場所は、起動時の案内と `pysad-share` が表示する)。"""
from setuptools import setup

setup()
