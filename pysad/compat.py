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


# ---------------------------------------------------------------------------------------------
# K1≠0 の BEND の K1 の定義 (v2.3.0 U2、B160・docs/11 A59)。K1 は中央面の勾配 g'(ρ0)/Bρ0 だけを決め、g''(ρ0) は決めない。
# 'log' (既定。K64 の定義): (1 + hx)∂b/∂x が一様。B_y = B0 + (K1/h) ln(r/ρ0)、B''(ρ0) = −h B'(ρ0)、x y² の項なし。
# 'lin': B' が x に一様 (B_y = B0 + B1 x、B_x = B1 y。elegant CSBEND の K2 = 0、MAD の SBEND と同じ)。'log' + 6 極 K2 = +h K1 の成分を本体に足したもの。
# v2 だけで効く (k64 は常に 'log')。効く範囲: δ≠0・オフセットのある粒子の BEND の本体 (追跡 TrackParticles・CALC・Emittance[]・
# MULT (ANGLE≠0) の K1 も同じ本体)。δ = 0 の設計軌道の光学・エミッタンス・放射積分は変わらない。
# 切り替え: 環境変数 PYSAD_BEND_K1_DEF=lin、または pysad.compat.set_bend_k1_def('lin')。
# ---------------------------------------------------------------------------------------------
_bend_k1_def = None


def bend_k1_def():
    """K1≠0 の BEND の K1 の定義 ('log' | 'lin')。k64 モードでは常に 'log'"""
    if mode() != 'v2':
        return 'log'
    if _bend_k1_def is not None:
        return _bend_k1_def
    m = os.environ.get('PYSAD_BEND_K1_DEF', 'log').strip().lower()
    return 'lin' if m == 'lin' else 'log'


def bend_k1_lin():
    """B' 一様の定義 ('lin') を使うか (v2 かつ bend_k1_def() == 'lin')"""
    return bend_k1_def() == 'lin'


def set_bend_k1_def(m):
    """K1≠0 の BEND の K1 の定義を設定する ('log' | 'lin')。None で環境変数に戻す"""
    global _bend_k1_def
    if m is not None and m not in ('log', 'lin'):
        raise ValueError("bend_k1_def must be 'log' or 'lin'")
    _bend_k1_def = m


# ---------------------------------------------------------------------------------------------
# B161 (v2.3.0 U2): CALC・Emittance[] の K1≠0 BEND の解析行列 (tbende.f tbendecorr) の dk (kick の x 微分) を、同じ関数の kick の導関数に直す。
# K64 は dk = −K1 xr (1 − xr (1 − 2 xr))、kick の導関数は −K1 xr (1 − xr/2 + xr²/3) (xr = x/ρ0)。v2 だけ直す (k64 は K64 のまま)。
# _bend_dk_fix は測定・試験用の上書き (None: v2 なら True、True/False: 強制。k64 では常に False)。
# ---------------------------------------------------------------------------------------------
_bend_dk_fix = None


def bend_dk_fix():
    if mode() != 'v2':
        return False
    return True if _bend_dk_fix is None else bool(_bend_dk_fix)


def set_bend_dk_fix(m):
    """測定・試験用: B161 の dk の修正を強制 (True/False)。None で既定 (v2 なら True) に戻す"""
    global _bend_dk_fix
    _bend_dk_fix = m


def bend_ecorr_flags():
    """CALC の BEND (K1≠0) の補正カーネル (perfkern_bend.bendecorr) に渡す整数の旗: 1 = K1 の定義 'lin' (B160)、2 = B161 の dk の修正"""
    return (1 if bend_k1_lin() else 0) | (2 if bend_dk_fix() else 0)
