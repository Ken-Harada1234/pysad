# -*- coding: utf-8 -*-
"""pysad の互換モード (v2.0.0 から)。

v1 系は SAD (K64 版) のソースに忠実な翻訳で、K64 のバグや癖も gs の実動作どおりに再現した。
v2 系は、K64 のバグの修正・物理の誤りの修正・機能の追加を行う。結果が gs と変わる箇所は、次の 2 つのモードで切り替える。

- ``v2`` (既定): 修正した動作。
- ``k64``: v1 と同じ、K64 (gs) の動作を再現する。gs との数値比較 (試験・回帰ベンチ) はこのモードで行う。

切り替え: 環境変数 ``PYSAD_COMPAT=k64`` または ``PYSAD_COMPAT=v2``、または ``pysad.compat.set_mode('k64')``。
コードでは、修正する箇所に ``if compat.v2(): <修正した動作> else: <K64 の動作>`` と書く。修正の一覧は docs/47_v2_changes.md。
"""
import os

_mode = None


def mode():
    """現在のモード ('v2' または 'k64')。環境変数 PYSAD_COMPAT を毎回読む (set_mode で上書き可)。"""
    if _mode is not None:
        return _mode
    m = os.environ.get('PYSAD_COMPAT', 'v2').strip().lower()
    return 'k64' if m in ('k64', 'v1', 'compat') else 'v2'


def set_mode(m):
    """モードを設定する ('v2' | 'k64')。None で環境変数に戻す。"""
    global _mode
    if m is not None and m not in ('v2', 'k64'):
        raise ValueError("mode must be 'v2' or 'k64'")
    _mode = m


def v2():
    return mode() == 'v2'


def k64():
    return mode() == 'k64'
