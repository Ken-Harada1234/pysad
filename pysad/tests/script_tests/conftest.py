# -*- coding: utf-8 -*-
"""試験の共通設定。

既存の試験は gs (SAD K64) との数値比較なので、互換モード k64 (K64 の動作を再現) で走らせる。
v2 で修正した動作の試験は、``compat.set_mode('v2')`` か ``PYSAD_COMPAT=v2`` を、試験の中で明示して使う
(フィクスチャ ``v2_mode``)。環境変数は子プロセスにも引き継がれる。
"""
import os

os.environ.setdefault('PYSAD_COMPAT', 'k64')
# SAD の Packages は、環境変数が無ければ pysad に同梱したコピー (pysad/sad_copy) を使う。GUI 系の試験が skip されないよう、
# 試験の中では環境変数の既定を同梱のコピーにする (元の SAD の Packages を試すときは、環境変数で上書きする)
os.environ.setdefault('PYSAD_SAD_PACKAGES', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'sad_copy', 'Packages'))

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


@pytest.fixture(autouse=True, scope='module')
def _close_leftover_tk_windows():
    """各試験ファイル (モジュール) のあと、残った Tk のルート窓 (tkinter.Tk()) を閉じる。GUI の試験の内容・結果は変えない。"""
    yield
    import gc
    import sys
    tk = sys.modules.get('tkinter')
    if tk is None:
        return
    roots = []
    try:
        roots = [o for o in gc.get_objects() if isinstance(o, tk.Tk)]
    except Exception:
        pass
    for r in roots:
        try:
            r.destroy()
        except Exception:
            pass
    try:
        tk._default_root = None
    except Exception:
        pass


@pytest.fixture(autouse=True, scope='module')
def _restore_global_state_per_module():
    """各試験ファイルのあと、グローバルの状態 (EMITDIV・互換モード・PYSAD_* の環境変数) を元に戻す。SAD 台本の `EMITDIV=20;` などが
    emit6d.EMITDIV を書き換えたまま後続のファイルの試験 (例: test_v2_P3_ae.py → test_emit_jac.py) に残る順序依存を防ぐ。
    ファイルの中の試験どうしの順序には触らない"""
    import os
    from pysad import compat, emit6d
    old_div = list(emit6d.EMITDIV)
    old_mode = compat._mode
    old_env = {k: v for k, v in os.environ.items() if k.startswith('PYSAD_')}
    yield
    emit6d.EMITDIV[:] = old_div
    compat.set_mode(old_mode)
    for k in [k for k in os.environ if k.startswith('PYSAD_') and k not in old_env]:
        del os.environ[k]
    os.environ.update(old_env)
