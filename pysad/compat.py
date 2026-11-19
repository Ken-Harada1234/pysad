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


# ---------------------------------------------------------------------------------------------
# K1≠0 の BEND の端の模型 (v2next P2 項目 5)。既定は 'radial' (K64 の仕様: 双極磁場は端面から、四極 (K1) の磁場は半径方向の入口・出口の面から始まる)。
# 'face' は四極の磁場も端面から始まる硬い端の模型 (端の角 e ≠ 0 のとき、x tan(e) の幅の楔形の領域の四極の蹴りが違う)。v2 だけで効く (k64 は K64 のまま)。
# 効く範囲: CALC (4D・6D)・Emittance[]・EMIT・TrackParticles の K1≠0 かつ ANGLE≠0 の BEND (DISFRIN=0。ずれ・ロールつきも)。
# 切り替え: 環境変数 PYSAD_BEND_K1_EDGE=face、または pysad.compat.set_bend_k1_edge('face')。
# ---------------------------------------------------------------------------------------------
_bend_k1_edge = None


def bend_k1_edge():
    """K1≠0 の BEND の端の模型 ('radial' | 'face')。k64 モードでは常に 'radial'"""
    if mode() != 'v2':
        return 'radial'
    if _bend_k1_edge is not None:
        return _bend_k1_edge
    m = os.environ.get('PYSAD_BEND_K1_EDGE', 'radial').strip().lower()
    return 'face' if m == 'face' else 'radial'


def set_bend_k1_edge(m):
    """K1≠0 の BEND の端の模型を設定する ('radial' | 'face')。None で環境変数に戻す"""
    global _bend_k1_edge
    if m is not None and m not in ('radial', 'face'):
        raise ValueError("bend_k1_edge must be 'radial' or 'face'")
    _bend_k1_edge = m
