"""Locations of the small data files bundled with pysad (pysad/data/), such as the KEK-LS lattice.
pysad に同梱した小さなデータファイル (pysad/data/) の場所を返す。

`data_path('kek-ls_v03_12.lat')`: KEK-LS (3 GeV、20 セル) のラティス。無いファイル名は FileNotFoundError。
`list_data()`: 同梱のデータファイル名の一覧。
"""
from __future__ import annotations

import os

_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def data_dir() -> str:
    """Return the pysad/data folder.
    pysad/data のフォルダ"""
    return _DIR


def list_data() -> list:
    """Return the names of the bundled data files (excluding README.md).
    同梱のデータファイル名 (README.md を除く)"""
    return sorted(f for f in os.listdir(_DIR) if os.path.isfile(os.path.join(_DIR, f)) and f != 'README.md')


def data_path(name: str) -> str:
    """Return the absolute path of a bundled data file; raise FileNotFoundError if it does not exist.
    同梱のデータファイルの絶対パス"""
    p = os.path.join(_DIR, name)
    if not os.path.isfile(p):
        raise FileNotFoundError(f'pysad に同梱のデータファイルが無い: {name} (あるもの: {", ".join(list_data())})')
    return p
