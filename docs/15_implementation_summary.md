# 15. pysad 実装まとめ

最終更新: 2026-11-06 (現状のスナップショット)。詳細は `01_overview.md` 〜 `14_fortran_status.md` と `work_*.md`。

## 1. 概要・所在

- 目的: SAD (**K64 fork**, Ver.1.4.4.2k64) を Python に直訳し、SAD の数値を再現する。SAD を物理モデルとして作り直さない。
  数値がずれたら SAD で実際に通るルーチンを特定し、行単位で訳す。SAD の癖 (バグに見えるもの) も再現する (→ `11_sad_quirks.md`)。
  SAD 自身のバグは再現せず `19_sad_bug_report.md` に集約する。
- 方針 (ユーザー指示): 1e-10 以下の差・表示の桁・丸めカオス後の乱数 seed は追わない。基準は「ユーザーのスクリプトが完走し、物理量が SAD と有意に一致する」。
  翻訳対象は Fortran/C の物理計算 (直接翻訳)。Packages/*.n は Python で書き直す。ファイル・リスト・GUI・描画は Python のライブラリで代替し翻訳しない (`15_sadhelp_purpose_policy.md`)。
- 場所: `<repo>\` (pysad 本体 `pysad\`、文書 `docs\`、作業ファイル `<work-perf>\`、ユーザースクリプトの照合 `<work-esrf>\` `<work-pf2>\`、配布 `pysad_20261106.tar.gz`)。
  SAD のソースは `<SAD>\`、参照値は WSL Ubuntu の `bin/gs` で取得 (`NPARA=1` を入れる)。
- 到達点: ユーザーの実スクリプト (`ユーザーの esrf ディレクトリ`, `ユーザーの pf2 ディレクトリ`) を `py -3.12 run_sad.py lattice script` でそのまま動かす
  (SADScript インタプリタ + FFS コマンド層)。86 本の現状は `12_user_scripts_status.md` (数値一致 62、一部一致 18、例外 1、SAD 側無効 4、対象外 1。<work-userscripts>/ の一体化ファイルで照合)。
- 実行環境: Python 3.12 (`py -3.12`) に numba / scipy / pytest。`PYTHONIOENCODING=utf-8`。長い計算は `<work-perf>/lowprio.py` (IDLE 優先度) で 1 本ずつ。
- 検証: `run_regression.py` (50 ベンチ、`QUICK=1` で短縮、FAIL 0)、pytest 277 件 pass (`--ignore` する 4 本は sys.exit で終わる旧式のスクリプト)。

## 2. モジュール構成と行数 (`wc -l pysad/*.py`、計 36821 行)

| ファイル | 行数 | 内容 (直訳元) |
|---|---:|---|
| lattice.py | 4038 | CALC: 4D 閉軌道 qcod、要素ごとの COD 写像 (tquade/tbende/tthine/tqlfre/tblfre/tchge/tbrote/tsteee/MULT)、twiss/ins_twiss (結合、不安定時の 1/cosμ、INS の初期軌道、TRPT の加速)、MARK OFFSET、解析行列と数値微分 (eps=1e-7) |
| ffscmd.py | 3568 | FFS コマンド層 (tffsa の解釈順)、Twiss[]/LINE[]/FFS[]/Element[]/Emittance[]/TrackParticles[]/CalculateOptics、EMIT 報告と結果変数、DISP/TYPE、誤差コマンド、RADINT、CHROMATICITY |
| sadscript.py | 2831 | SADScript インタプリタ (字句・Pratt 構文・評価、Module/Do/Table/純関数/パターン/REPEAT/UNTIL/IF、単位付き数値、`?FLAG`、autofg 書式) |
| emit6d.py | 2265 | Emittance[]/EMIT (temit/tturne/tradke/tecalc、RFSW × RADCOD、close_z、CAVI の DX/DY/ROTATE、EMITXE.. の下限) |
| track6.py | 2054 | TrackParticles の直訳 (純 Python。基準実装。tdrift/tquad/tthin/tbend/tbendi/tsteer/tcav/tmulti/tsol、tradk/tradkfn) |
| parser.py | 1876 | MAIN 部 (要素定義・LINE・単位・再定義・`-LINE`・`N*`) |
| track6n.py | 1826 | track6 の numba 版 (176 ns/粒子・要素)、粒子損失 (tapert、TrackRun) |
| elkern.py | 1723 | 要素カーネル (COD 写像とヤコビアンの numba 化) |
| match6d.py | 1468 | 解析微分 (qdcell/tffsqu)。`Matcher` の既定 |
| match6.py | 1421 | マッチング (tffsmatch/tffscalc/tdfun/twfit/tweigh/tfsolv/fmincube、FITP>1 と fam 点、フック、OpticsProlog/Epilog) |
| sadlist.py | 1167 | リスト・文字列関数 (functions.n, Misc.n 相当) |
| elements.py | 909 | 要素クラスと線形行列 |
| touschek.py | 849 | INTRA (IBS 付き emit)、TouschekTable/TouschekLifetime |
| track6n_ext.py | 696 | MULT と SOL 領域の numba カーネル |
| calc6d.py | 660 | CALC6D (qcell6d、CELL と INS)、etwiss2ri |
| radint.py | 630 | RADINT (INTGRL.f) |
| geo.py | 556 | 幾何 (tfgeo/geolib: GX.., GEO, ORG, DISP G/OG, fit の GX..) |
| teigen.py | 490 | temit の固有値分解 teigen/tnorm の直訳 (不安定な光学) |
| sadelem.py | 427 | Element/SetElement/BeamLine/WriteBeamLine |
| tsvd.py | 404 | tsolvg.f の SVD (2 段目フィット) |
| maintrack.py | 355 | MAIN 層の TRACK コマンド |
| perfkern.py / perfkern_bend.py | 354 / 273 | 4D 閉軌道の 1 周を numba 1 回にするカーネル |
| trptacc.py | 334 | TRPT の加速 (tfacc、tcave/tmulte の acc) |
| optics_fn.py | 315 | TransferMatrix/BeamMatrix/NormalCoordinates/OrbitGeo/GeoBase/SymplecticJ/SurvivedParticles/AccelerateParticles (.n の書き直し) |
| emit_extra.py | 305 | Emittance[] のオプション (OneTurnInformation/Matrix 系) |
| diag_trace.py / diag_fit.py | 266 / 244 | 照合用の診断 |
| sadrandom.py | 158 | SAD の乱数 (LCG, Box–Muller, tdusrn, ParabolaRandom) |
| dapert.py | 146 | DynamicApertureSurvey |
| turnkern.py | 137 | 1 周の Twiss 伝搬の numba カーネル |
| `__init__.py` / sadkw.py / constants.py | 121 / 34 / 10 | |
| **非推奨**: ffs.py / tracking.py / matching.py / optics.py | 1599 / 1239 / 602 / 471 | 独自実装 (SAD 未照合。tracking.py は放射係数に疑い)。新規コードは使わない。`__init__.py` の export は旧 API のまま |

テストは `pysad/tests/` (test_pysad.py ほか 28 ファイル)。

## 3. SAD 準拠状況

| 領域 | 状態 | 直訳元 |
|---|---|---|
| CALC (4D 閉軌道・Twiss・結合・オフモメンタム・INS・不安定ラインの続行) | 一致 (オンモメンタム 1e-11、COD/結合経由は相対 1e-7) | qcod, qcell1, qtwiss, qmdiag |
| 要素の COD 写像 (四極フリンジ、薄い多極、BEND の K0/K1/ロール/F1、ANGLE=0 BEND、MULT、ROTATE/ずれ付きの全要素) | 一致 | tquade, tthine, tchge, tbende, tsteee, tqlfre |
| CALC6D (CELL/INS)、幾何 (GX.., GEO, ORG)、RADINT | 一致 (1.7e-10、4.8e-13、7 桁) | qcell6d, tfgeo, INTGRL |
| TRPT の加速と MARK OFFSET (線形加速器) | 一致 (disp_bt_r の 3546 行の表が印字精度で同じ) | tfacc, tcave/tmulte acc, qcell |
| Emittance[]/EMIT (RFSW × RADCOD、空洞 L≠0、MULT、不安定な光学、報告の全行、結果変数) | 一致 (≤1e-11〜1e-9) | temit, tturne, tradke, tecalc, teigen |
| INTRA / Touschek、Emittance[] の拡張 | 一致 (≤7e-12、2.1e-10) | tintrb, emit.n |
| 粒子追跡 (NORAD/RAD/FLUC、K1≠0 BEND、MULT、SOL、粒子損失) | ≤5e-11、乱数の最終状態まで一致 (損失のない場合)。損失周は論理一致、実寸 DAsv の最終 seed はカオスで不一致 | tturn1 系, photons, tapert |
| マッチング FIT/FREE/FIX/GO (解析微分、DESCEND、fuzz、tsvd、2 段目、INS、フック) | 反復列・残差・K1 を SAD と同じ桁で再現 (esrf/matching で 12 桁) | tffsmatch, qdcell, tsolvg |
| FITP>1 (idp ループ、重み、fam 点) | 直訳のまま SAD と未照合 (SAD (WSL) が GO で再起動するため参照値なし)。1 段目の解は一致 | tffscalc, tffamsetup |
| DynamicApertureSurvey | NORAD は地図・Score 完全一致。RAD+FLUC の実寸は Score が 1〜3 % 違う (乱数のカオス) | trackd, tfdapert |
| 誤差コマンド (DELX.. の実現、SeedRandom) と CodCorr | 誤差の実現 (乱数列) は一致。補正量は SAD の teigen の縮退で丸め依存 | terror, tfadjst |
| SADScript / FFS 層 | ユーザースクリプト 86 本中 数値一致 62、一部一致 18 (`12_user_scripts_status.md`) | tfeval, tffsa |
| 光学関数 (TransferMatrix, BeamMatrix, NormalCoordinates, OrbitGeo, GeoBase, SymplecticJ, SurvivedParticles, AccelerateParticles) | 一致 (1e-9) | Optics.n, emit.n, dapert.n |

MON/MCOD/TCOD はこの SAD ビルドでは何もしない (式として評価されるだけ) ので pysad も受理のみ。

## 4. 主な仕様 (直訳で再現している SAD の扱い)

- **要素値**: `LINE["DX","QD.1"]=v` は CALC で全 QD に展開される (ファミリーの代表)。1 台だけの誤差は `DELX σ QD.1`。pysad API では `set_element` (ファミリー) と `set_error` (1 台) を分ける。
- **INS の初期軌道**: 先頭 MARK の DX/DPX/DY/DPY/DZ/DDP が初期座標 (誤差ではない)。
- **マッチング**: EMIT が失敗して EMITX/EMITY=NaN のときは fuzz の MAX/MIN が NaN を無視する gfortran の挙動を再現し、emminv 由来の ~1e-17 で EX/EPX の重みが 1e5 倍になる (esrf/matching が SAD と同じ解に着く条件)。fmincube の定数はオイラー定数 γ。
- **色収差**: SAD の CHRO は BEND 成分を含まない。ユーザースクリプトは DP0 スキャン + PolynomialFit を使う (Fit.n の式をそのまま評価)。
- **EMIT コマンド**は GAUSS フラグを立てる (Emittance[] は変えない)。SADScript の `1GEV` は SAD では `1*GEV`=1 だが pysad は単位として読む (意図的な拡張)。
- **FFS 内の ON/OFF の語**は無効 (`INTRA;` / `NOINTRA;` と書く)。MAIN 部の ON/OFF は有効。
- **乱数**: 既定の "SAD" 乱数 (LCG + Box–Muller) を int64 で写し、FLUC の引く順 (要素外側・粒子内側) まで同じ。`SeedRandom[]` の戻りは 5 要素。
- その他の癖は `11_sad_quirks.md` (A1〜A39, B1〜B6, C1〜C24, D1〜)。

## 5. 性能

- 1 回の cell calc 0.03 s (不変時) / 0.145 s (変更後)、Twiss[] 0.6 ms、追跡 176 ns/粒子・要素 (numba)。
- 600 s 内に終わらないスクリプト: 誤差入り 6D emit × 多数の実現 (print_all_error 300 実現 221 s、SAD 37 s)、DAsv 5000 周 × 969 粒子 (約 50 分、SAD も 300 s 超)、CodCorr の全実現。
  残る改善点は 6D の tturne 1 周 (0.03 s) と DAsv 系 (`history/work_perf2.md`, `history/work_perf3.md`)。

## 6. 未対応・制限 (優先順は `history/00_RESUME.md` に従う)

- 性能 (上記)。FITP>1 の最終解の照合 (SAD が再起動)。誤差系の要調査 2 件 (docs/12 §4)。
- 不安定ラティスの 1e199 級の値や NaN は再現しない。実寸 DAsv の最終 seed (カオス)。
- 要素: SOL の emit/幾何 (追跡のみ)、MULT の ANGLE≠0 (tmulta)、TCAVI、COORD、MAP、UND、BEAMBEAM、POL (偏極)、加速要素の K_n と VOLT の併用。
- 関数: SynchroBetaEmittance (SAD K64 に実体なし)、InitialOrbits/ElementValues (使用なし)、GaussianCoulomb/WakeFunction/RadiationField/ExternalMap (優先度 低)。
- 対象外: Tk/plot/EPICS/System、描画 (DRAW/GEO 描画)、ファイル入出力の失敗の扱い。
