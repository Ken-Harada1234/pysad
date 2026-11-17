# -*- coding: utf-8 -*-
"""試験の共通設定。

既存の試験は gs (SAD K64) との数値比較なので、互換モード k64 (K64 の動作を再現) で走らせる。
v2 で修正した動作の試験は、``compat.set_mode('v2')`` か ``PYSAD_COMPAT=v2`` を、試験の中で明示して使う
(フィクスチャ ``v2_mode``)。環境変数は子プロセスにも引き継がれる。
"""
import os

os.environ.setdefault('PYSAD_COMPAT', 'k64')

import pytest


@pytest.fixture
def v2_mode():
    from pysad import compat
    old = compat._mode
    compat.set_mode('v2')
    try:
        yield
    finally:
        compat.set_mode(old)


@pytest.fixture
def k64_mode():
    from pysad import compat
    old = compat._mode
    compat.set_mode('k64')
    try:
        yield
    finally:
        compat.set_mode(old)
