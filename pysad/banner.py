"""Startup banner shown once on standard error when the SAD interpreter is first called from Python.
起動時の案内 (バナー)。

Python から SAD のインタプリタ (`Interp`、`run_sad`、`pysad-run`、FFS・SADScript の実行) を呼んだとき、プロセスの中で最初の 1 回だけ、
数行を標準エラー出力に出す。バージョン、元になっている SAD (K. Oide 氏の K64 版) のこと、ソース・例・文書・物理ガイド・試験が GitHub (と GitHub Pages) にあること (pip には入れていない)。
`import pysad` だけでは出さない。標準出力には出さないので、台本の出力 (標準出力) は変わらない。

環境変数 `PYSAD_QUIET=1` で出さない (試験・自動実行向け)。`pysad.banner.show(force=True)` で、もう一度出せる。
"""
from __future__ import annotations

import os
import sys

URL = 'https://github.com/Ken-Harada1234/pysad'
PAGES = 'https://ken-harada1234.github.io/pysad/'
_shown = False


def text() -> str:
    """Return the banner text (a few lines).
    バナーの文字列 (数行)"""
    from . import __version__, compat
    mode = compat.mode() if hasattr(compat, 'mode') else getattr(compat, '_mode', 'v2')
    return (
        f'pysad {__version__} -- a Python translation of SAD (Strategic Accelerator Design), K64 version. compat mode: {mode}\n'
        'Based on SAD by K. Oide (KEK), K64 version: the physics and the SADScript language follow SAD.\n'
        'SAD packages (*.n) bundled in this package are copied from the SAD source (KEK Accelerator Laboratory, BSD license; see pysad/sad_copy/README.md).\n'
        f'Source code: {URL}\n'
        f'Examples, documentation, the physics guide and the tests are not in the pip package; see GitHub ({URL}) and {PAGES} (to hide this message: PYSAD_QUIET=1)\n'
    )


def show(force: bool = False, file=None) -> bool:
    """Print the banner to standard error only the first time, returning True if printed (suppressed when `PYSAD_QUIET` is set).
    最初の 1 回だけ、バナーを標準エラー出力に出す。出したとき True。`PYSAD_QUIET` が設定されているときは出さない"""
    global _shown
    if _shown and not force:
        return False
    _shown = True
    if os.environ.get('PYSAD_QUIET', '').strip().lower() in ('1', 'true', 'yes', 'on') and not force:
        return False
    f = file if file is not None else sys.stderr
    try:
        f.write(text())
        f.flush()
    except Exception:
        return False
    return True
