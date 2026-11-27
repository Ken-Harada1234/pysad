# -*- coding: utf-8 -*-
"""Compatibility modes of pysad: `v2` (fixed behavior, default) and `k64` (reproduces K64 gs behavior).
pysad の互換モード (v2.0.0 から)。

v1 系は SAD (K64 版) のソースに忠実な翻訳で、K64 のバグや癖も gs の実動作どおりに再現した。
v2 系は、K64 のバグの修正・物理の誤りの修正・機能の追加を行う。結果が gs と変わる箇所は、次の 2 つのモードで切り替える。

- ``v2`` (既定): 修正した動作。
- ``k64``: v1 と同じ、K64 (gs) の動作を再現する。gs との数値比較 (試験・回帰ベンチ) はこのモードで行う。

切り替え: 環境変数 ``PYSAD_COMPAT=k64`` または ``PYSAD_COMPAT=v2``、または ``pysad.compat.set_mode('k64')``。
コードでは、修正する箇所に ``if compat.v2(): <修正した動作> else: <K64 の動作>`` と書く。修正の一覧は docs/28_v2_vs_k64.md。
"""
import os

_mode = None


def mode():
    """Return the current compatibility mode ('v2' or 'k64'), reading the PYSAD_COMPAT environment variable each time.
    現在のモード ('v2' または 'k64')。環境変数 PYSAD_COMPAT を毎回読む (set_mode で上書き可)。"""
    if _mode is not None:
        return _mode
    m = os.environ.get('PYSAD_COMPAT', 'v2').strip().lower()
    return 'k64' if m in ('k64', 'v1', 'compat') else 'v2'


def set_mode(m):
    """Set the compatibility mode ('v2' or 'k64'); None reverts to the environment variable.
    モードを設定する ('v2' | 'k64')。None で環境変数に戻す。"""
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
    """Return the edge model of a BEND with K1 nonzero ('radial' or 'face'); always 'radial' in k64 mode.
    K1≠0 の BEND の端の模型 ('radial' | 'face')。k64 モードでは常に 'radial'"""
    if mode() != 'v2':
        return 'radial'
    if _bend_k1_edge is not None:
        return _bend_k1_edge
    m = os.environ.get('PYSAD_BEND_K1_EDGE', 'radial').strip().lower()
    return 'face' if m == 'face' else 'radial'


def set_bend_k1_edge(m):
    """Set the edge model of a BEND with K1 nonzero ('radial' or 'face'); None reverts to the environment variable.
    K1≠0 の BEND の端の模型を設定する ('radial' | 'face')。None で環境変数に戻す"""
    global _bend_k1_edge
    if m is not None and m not in ('radial', 'face'):
        raise ValueError("bend_k1_edge must be 'radial' or 'face'")
    _bend_k1_edge = m


# ---------------------------------------------------------------------------------------------
# K1≠0 の BEND の K1 の定義 (v2.3.0 U2、B160・docs/29 A59)。K1 は中央面の勾配 g'(ρ0)/Bρ0 だけを決め、g''(ρ0) は決めない。
# 'log' (既定。K64 の定義): (1 + hx)∂b/∂x が一様。B_y = B0 + (K1/h) ln(r/ρ0)、B''(ρ0) = −h B'(ρ0)、x y² の項なし。
# 'lin': B' が x に一様 (B_y = B0 + B1 x、B_x = B1 y。elegant CSBEND の K2 = 0、MAD の SBEND と同じ)。'log' + 6 極 K2 = +h K1 の成分を本体に足したもの。
# v2 だけで効く (k64 は常に 'log')。効く範囲: δ≠0・オフセットのある粒子の BEND の本体 (追跡 TrackParticles・CALC・Emittance[]・
# MULT (ANGLE≠0) の K1 も同じ本体)。δ = 0 の設計軌道の光学・エミッタンス・放射積分は変わらない。
# 切り替え: 環境変数 PYSAD_BEND_K1_DEF=lin、または pysad.compat.set_bend_k1_def('lin')。
# ---------------------------------------------------------------------------------------------
_bend_k1_def = None


def bend_k1_def():
    """Return the definition of K1 for a BEND with K1 nonzero ('log' or 'lin'); always 'log' in k64 mode.
    K1≠0 の BEND の K1 の定義 ('log' | 'lin')。k64 モードでは常に 'log'"""
    if mode() != 'v2':
        return 'log'
    if _bend_k1_def is not None:
        return _bend_k1_def
    m = os.environ.get('PYSAD_BEND_K1_DEF', 'log').strip().lower()
    return 'lin' if m == 'lin' else 'log'


def bend_k1_lin():
    """Return whether the uniform-B' definition of K1 ('lin') is used (v2 and bend_k1_def() == 'lin').
    B' 一様の定義 ('lin') を使うか (v2 かつ bend_k1_def() == 'lin')"""
    return bend_k1_def() == 'lin'


def set_bend_k1_def(m):
    """Set the definition of K1 for a BEND with K1 nonzero ('log' or 'lin'); None reverts to the environment variable.
    K1≠0 の BEND の K1 の定義を設定する ('log' | 'lin')。None で環境変数に戻す"""
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
    """Force the B161 dk fix on or off for measurement and tests; None restores the default (True in v2).
    測定・試験用: B161 の dk の修正を強制 (True/False)。None で既定 (v2 なら True) に戻す"""
    global _bend_dk_fix
    _bend_dk_fix = m


def bend_ecorr_flags():
    """Return the integer flags passed to the BEND (K1 nonzero) correction kernel `perfkern_bend.bendecorr` of CALC.
    CALC の BEND (K1≠0) の補正カーネル (perfkern_bend.bendecorr) に渡す整数の旗: 1 = K1 の定義 'lin' (B160)、2 = B161 の dk の修正"""
    return (1 if bend_k1_lin() else 0) | (2 if bend_dk_fix() else 0) | (4 if v2() else 0)       # 4 = B250 (xsxkx の符号。bendebody だけが見る。v2 のとき)


# ---------------------------------------------------------------------------------------------
# B167 (v2.4.0): β < 1 の粒子の放射 (陽子・イオン、低い γ)。K64 の放射 (tradkn・tradke・RADINT) は β = 1 の式で、平均のエネルギー損失が古典式の 1/β³ 倍、
# 運動量の損失が Δp = ΔE/c (正しくは ΔE/(βc))、量子励起の 2 次モーメントが厳密なスペクトルと違う。v2 は厳密な放射のスペクトル (Schott の式) から決めた
# 係数 (pysad/radbeta.py) を掛ける。v2 だけ (k64 は K64 のまま)。γ ≥ 300 (電子の全部) には入れない (補正は 3e-5 以下)。
# 切り替え: 環境変数 PYSAD_RAD_BETA=0 で補正なし、または pysad.compat.set_radiation_beta(False)。
# ---------------------------------------------------------------------------------------------
_radiation_beta = None


def radiation_beta():
    if mode() != 'v2':
        return False
    if _radiation_beta is not None:
        return bool(_radiation_beta)
    return os.environ.get('PYSAD_RAD_BETA', '1').strip().lower() not in ('0', 'false', 'off', 'no')


def set_radiation_beta(m):
    """Force the beta factor of B167 on or off; None restores the default (True in v2).
    B167 の β の因子を強制 (True/False)。None で既定 (v2 なら True) に戻す。Emittance[]・追跡に使う _Phys は作った時点の値を持つ"""
    global _radiation_beta
    _radiation_beta = m


# ---------------------------------------------------------------------------------------------
# B165 (v2.4.0): 磁石の本体の運動量の向きの磁場の積分 (1 + a) B∥ ds (B136・B137)。K64 の追跡 (tradk の sprot)・Emittance[] (tradke) は、運動量の向きの変化
# (B⊥) と端の縦磁場 ψ (BEND・QUAD の入口 +ψ、出口 −ψ) だけでスピンを歳差させ、本体の中の B∥ = B·p̂ を入れない。静磁場の B∥ の経路積分は
# ∫B·dr = ψ_出口 − ψ_入口 (ψ = y g(r)、スカラーポテンシャルの差。BEND は y/ρ + k x y、QUAD は k x y) なので、y・x の傾きがあると
# 端の ψ だけが残って誤る (py = 5e-3 の 15° 扇形 BEND で R[0,1] が −4.05e-3、Maxwell 場つきの硬い端の Thomas-BMT の積分は −2.75e-3)。
# v2 は本体の各刻み (放射・スピンの刻み) に (ψ_終 − ψ_始)/(1+δ) を運動量の向きの磁場として足す。v2 だけ (k64 は K64 のまま)。
# 範囲: BEND (K1=0・K1≠0)・QUAD・MULT・SOL 領域内の四極の追跡と Emittance[] (v2.4.0 D3: B230・B231・B234)。SEXT・ステアリングは K64 が端の ψ も本体も入れず和が 0 で元から正しい。
# POL ON の追跡は 1 刻みの偏向角を 5e-3 rad 以下に刻む (B235。追跡の `_ndivrad`)。
# 切り替え: 環境変数 PYSAD_SPIN_BPAR=0 で無効、または pysad.compat.set_spin_bpar(False)。
# ---------------------------------------------------------------------------------------------
_spin_bpar = None


def spin_bpar():
    if mode() != 'v2':
        return False
    if _spin_bpar is not None:
        return bool(_spin_bpar)
    return os.environ.get('PYSAD_SPIN_BPAR', '1').strip().lower() not in ('0', 'false', 'off', 'no')


def set_spin_bpar(m):
    """Force the B-parallel treatment in the main part of B165 on or off; None restores the default (True in v2).
    B165 の本体の B∥ を強制 (True/False)。None で既定 (v2 なら True) に戻す"""
    global _spin_bpar
    _spin_bpar = m


# ---------------------------------------------------------------------------------------------
# B137 (v2): BEND の端の縦磁場 ψ の歳差を、刻みの本体の歳差と別の回転として順に掛ける (追跡)。K64 の tradk の sprot は入口の ψ を最初の刻みの、
# 出口の ψ を最後の刻みの回転ベクトルに足して 1 回で回す。回転は非可換なので、刻み 1 つ分の偏向角 θ_n に比例する誤差 (ψ θ_n の 1 次) が残り、
# スピン軌道結合 D (y・py の列) が刻みの数に反比例してずれる (0.5 GeV・15° 扇形・K1=0 で D[x,y] の誤差が −D/n)。v2 は入口の ψ を最初の刻みの本体の前に、
# 出口の ψ を最後の刻みの本体の後に、独立な回転として掛ける。v2 かつ B165 (spin_bpar) が有効なときだけ。k64 は K64 のまま。
# 切り替え: 環境変数 PYSAD_SPIN_EDGE_SPLIT=0 で無効、または pysad.compat.set_spin_edge_split(False)。
# ---------------------------------------------------------------------------------------------
_spin_edge_split = None


def spin_edge_split():
    if not spin_bpar():
        return False
    if _spin_edge_split is not None:
        return bool(_spin_edge_split)
    return os.environ.get('PYSAD_SPIN_EDGE_SPLIT', '1').strip().lower() not in ('0', 'false', 'off', 'no')


def set_spin_edge_split(m):
    """Force the separate edge rotation of B137 on or off; None restores the default (True in v2 with B165).
    B137 の端の ψ の独立な回転を強制 (True/False)。None で既定 (v2 かつ B165 なら True) に戻す"""
    global _spin_edge_split
    _spin_edge_split = m


# ---------------------------------------------------------------------------------------------
# EMITDIV の既定 (v2.4.0)。Emittance[]・EMIT の放射の刻みの倍率。K64 の既定は 1 (k64 は常に 1)。
# 1 では 1 刻みの放射の入れ方が粗く、εx が +1e-3〜+2e-2 大きい (大きなリング A +2.2%、PF +2.2e-3、矩形磁石の小リング +1.4e-3。
# EMITDIV = 4 で 2e-4 以下、実行時間の増加は 0.1 s 以下。docs/28 v2.4.0)。v2 の既定は 4 (大域変数 EMITDIV を台本で与えれば、その値)。
# 環境変数 PYSAD_EMITDIV_DEFAULT で変えられる (1 で K64 と同じ)。EMITDIVB・EMITDIVQ の既定は 1 のまま。
# ---------------------------------------------------------------------------------------------
def emitdiv_default(pol=False):
    """Return the default EMITDIV used when the script does not give EMITDIV (1 in k64, 4 in v2, 16 for POL ON spin).
    EMITDIV の既定 (台本が EMITDIV を与えないとき)。k64 は 1、v2 は 4 (POL ON のスピンの計算は 16: 線形の脱偏極率の刻みの誤差が EMITDIV = 1 で +24%、
    4 で +20%、16 で +4% (64 で収束 1.87e-8。4.2 GeV・12 MV の FODO。docs/history/v24_RAD.md §4))"""
    if mode() != 'v2':
        return 1.0
    d0 = 16.0 if pol else 4.0
    try:
        v = float(os.environ.get('PYSAD_EMITDIV_DEFAULT', str(d0)))
    except ValueError:
        v = d0
    return v if v > 0.0 else d0
