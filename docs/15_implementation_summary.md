# 15. pysad 実装まとめ

最終更新: 2026-11-16 (現状のスナップショット)。詳細は `01_overview.md` 〜 `14_fortran_status.md`、要素・機能別の `20_*.md`〜`35_*.md`、大項目の解説 `36_*.md`〜`43_*.md`、`history/work_*.md`。

## 1. 概要・所在

- 目的: SAD (**K64 fork**, Ver.1.4.4.2k64) を Python に直訳し、SAD の数値を再現する。SAD を物理モデルとして作り直さない。
  数値がずれたら SAD で実際に通るルーチンを特定し、行単位で訳す。SAD の癖 (バグに見えるもの) も再現する (→ `11_sad_quirks.md`)。
  SAD 自身のバグは再現せず `19_sad_bug_report.md` に集約する。
- 方針 (ユーザー指示): 1e-10 以下の差・表示の桁・丸めカオス後の乱数 seed は追わない。基準は「計算が完走し、物理量が SAD と有意に一致する」。
  翻訳対象は Fortran/C の物理計算 (直接翻訳)。Packages/*.n は Python で書き直す。ファイル・リスト・GUI・描画は Python のライブラリで代替し翻訳しない (`15_sadhelp_purpose_policy.md`)。
- 場所: `<repo>\` (pysad 本体 `pysad\`、文書 `docs\`、作業ファイル `<work-perf>\`、配布 `pysad_20261106.tar.gz` は 11/06 の版で古い)。
  SAD のソースは `<SAD>\`、参照値は WSL Ubuntu の `bin/gs` で取得 (`NPARA=1` を入れる)。
- 到達点: SAD の台本 (SADScript + FFS) を `py -3.12 run_sad.py lattice script` でそのまま動かす (SADScript インタプリタ + FFS コマンド層)。
- 実行環境: Python 3.12 (`py -3.12`) に numba / scipy / pytest。`PYTHONIOENCODING=utf-8`。長い計算は `<work-perf>/lowprio.py` (IDLE 優先度) で 1 本ずつ。
- 検証 (2026-11-16 夜の記録 `local/wake/verify_log.md`。今回の文書更新では再実行していない): `PYTHONIOENCODING=utf-8 py -3.12 -u run_regression.py` (QUICK なし) で **50 ベンチをすべて実行して FAIL 0** (205 s。QUICK で SKIP される 7・12t・13a・16a も PASS)。
  `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` で **520 通過・3 skip (523 件収集)**。`--ignore` する 3 本は sys.exit で終わる旧式のスクリプト)。
  GUI 系 (test_plot、test_tclbridge) は全体実行で数回に 1 回失敗することがある (単独なら通る)。`PYSAD_SAD_PACKAGES=D:/work/SAD-master/SAD-master/Packages` が無いと GUI 系は skip。`PYSAD_SLOW=1` で遅い試験 (CSR MakeZL 約 7 分)。

## 2. モジュール構成と行数 (`wc -l pysad/*.py`、67 本・計 50169 行。2026-11-16 に数え直した)

| ファイル | 行数 | 内容 (直訳元) |
|---|---:|---|
| lattice.py | 4114 | CALC: 4D 閉軌道 qcod、要素ごとの COD 写像 (tquade/tbende/tthine/tqlfre/tblfre/tchge/tbrote/tsteee/MULT)、twiss/ins_twiss (結合、不安定時の 1/cosμ、INS の初期軌道、TRPT の加速)、MARK OFFSET、解析行列と数値微分 (eps=1e-7) |
| ffscmd.py | 4547 | FFS コマンド層 (tffsa の解釈順)、RADLIGHT (`_RadLight`)、旧 WAKE の経路、Twiss[]/LINE[]/FFS[]/Element[]/Emittance[]/TrackParticles[]/CalculateOptics、EMIT 報告と結果変数、DISP/TYPE、誤差コマンド、RADINT、CHROMATICITY |
| sadscript.py | 4248 | SADScript インタプリタ (字句・Pratt 構文・評価、Module/Do/Table/純関数/パターン/REPEAT/UNTIL/IF、単位付き数値、`?FLAG`、autofg 書式) |
| emit6d.py | 2422 | Emittance[]/EMIT (temit/tturne/tradke/tecalc、RFSW × RADCOD、close_z、CAVI の DX/DY/ROTATE、EMITXE.. の下限) |
| track6.py | 2209 | TrackParticles の直訳 (純 Python。基準実装。tdrift/tquad/tthin/tbend/tbendi/tsteer/tcav/tmulti/tsol、tradk/tradkfn) |
| parser.py | 1908 | MAIN 部 (要素定義・LINE・単位・再定義・`-LINE`・`N*`) |
| track6n.py | 2124 | track6 の numba 版 (176 ns/粒子・要素)、粒子損失 (tapert、TrackRun) |
| elkern.py | 1723 | 要素カーネル (COD 写像とヤコビアンの numba 化) |
| match6d.py | 1563 | 解析微分 (qdcell/tffsqu)。`Matcher` の既定 |
| match6.py | 1449 | マッチング (tffsmatch/tffscalc/tdfun/twfit/tweigh/tfsolv/fmincube、FITP>1 と fam 点、フック、OpticsProlog/Epilog) |
| sadlist.py | 1405 | リスト・文字列関数 (functions.n, Misc.n 相当) |
| elements.py | 884 | 要素クラスと線形行列 |
| touschek.py | 857 | INTRA (IBS 付き emit)、TouschekTable/TouschekLifetime |
| track6n_ext.py | 917 | MULT と SOL 領域の numba カーネル |
| calc6d.py | 671 | CALC6D (qcell6d、CELL と INS)、etwiss2ri |
| radint.py | 630 | RADINT (INTGRL.f) |
| geo.py | 573 | 幾何 (tfgeo/geolib: GX.., GEO, ORG, DISP G/OG, fit の GX..) |
| teigen.py | 490 | temit の固有値分解 teigen/tnorm の直訳 (不安定な光学) |
| sadelem.py | 435 | Element/SetElement/BeamLine/WriteBeamLine |
| tsvd.py | 404 | tsolvg.f の SVD (2 段目フィット) |
| maintrack.py | 355 | MAIN 層の TRACK コマンド |
| perfkern.py / perfkern_bend.py | 370 / 273 | 4D 閉軌道の 1 周を numba 1 回にするカーネル |
| trptacc.py | 334 | TRPT の加速 (tfacc、tcave/tmulte の acc) |
| optics_fn.py | 315 | TransferMatrix/BeamMatrix/NormalCoordinates/OrbitGeo/GeoBase/SymplecticJ/SurvivedParticles/AccelerateParticles (.n の書き直し) |
| emit_extra.py | 340 | Emittance[] のオプション (OneTurnInformation/Matrix 系) |
| diag_trace.py / diag_fit.py | 266 / 244 | 照合用の診断 |
| sadrandom.py | 158 | SAD の乱数 (LCG, Box–Muller, tdusrn, ParabolaRandom) |
| dapert.py | 146 | DynamicApertureSurvey |
| turnkern.py | 137 | 1 周の Twiss 伝搬の numba カーネル |
| csr.py | 1067 | CSR (csrtrack.f の CSRTrack・CSRHaissin、csroy.f の CSRSetupOY・CSROYMatrix、csrtest.f の CSRMatrix mode 1000 まで・CSRConvert・CSRInit、実 FFT)。docs/34、38 |
| synchrobeta.py | 651 | SynchroBetaEmittance (temits.f 一式) |
| beambeam.py | 731 | BEAMBEAM (bbstrhl.f、bbtool.f、matrix.f。BSTRL、光子、LuminosityMonitor)。docs/22、37 |
| fitroot.py | 618 | FindRoot / Fit / FitGaussian / FitEmit / NIntegrate (tffindroot.f、Fit.n)。docs/31 |
| specfn.py | 569 | 特殊関数・数学関数 (scipy.special で K64 と同じ関数名。K64 が誤値を返す式は真値)。docs/31 |
| gausscoul.py | 472 | GaussianCoulomb / U / Fitted (tspac.f、tintrb.f の twspu)。docs/30 |
| ffsextra.py | 463 | FFS の小さなコマンド (VARIABLES・ATTRIBUTE・MATRIX・SHOW・VARY・RECOVER・REVERSE・PRINT・EXECUTE ほか)。docs/33 |
| spin.py | 448 | POL のスピン (tradke のスピン部、spnorm・srequpol・spdepol、tchge の srot)。docs/34、41 |
| sadclass.py | 424 | SAD のクラス (Class・Constructor・`x@member`・This・継承) |
| tkcanvas.py | 423 | Tk の Canvas の組み込み関数 (TkCanvas*、CanvasClipLine)。docs/35、40 |
| miscfn.py | 414 | DownhillSimplex ほか optimize.n の移植、StringToStream・TemporaryName など |
| tclbridge.py | 322 | Tcl/Tk の接点 (TclCreateInterp・Tcl・TclArg・TkEventLoop を tkinter で)。docs/20、39 |
| wake.py | 308 | WakeFunction の追跡 (txwake.f)。docs/21、36 |
| corr.py | 285 | ResponseMatrix・LinearSolveConditioned (軌道補正の組み込み関数)。docs/32 |
| fileops.py | 275 | FileSystem.n のファイル・ディレクトリ関数 (Python の os・shutil) |
| sgeo.py | 237 | SOL 領域の幾何 (tsgeo.f)。docs/27 |
| multae.py / multa.py | 224 / 206 | MULT ANGLE≠0 の Twiss・Emittance[] / 追跡 (tmulta.f)。docs/26 |
| wakeopt.py | 200 | 光学の wake (tffswake.f、WAKEOPT)。docs/21 |
| tcav.py | 191 | TCAVI (ttcav.f、ttcave.f)。docs/24 |
| solemit.py / solmap.py | 168 / 84 | SOL 領域の Emittance[]・EMIT / Twiss (tsole1)。docs/27、42 |
| bbbrem.py | 162 | BBBrem1 (tbbbrem.f)。docs/30 |
| coord.py | 147 | COORD (tcoord.f、tcoorde.f、qcoord.f)。docs/23 |
| radfield.py | 139 | RadiationField / RadiationSpectrum (tlinit.f)。docs/43 |
| pkgload.py | 124 | SAD の Packages/*.n を SADScript のまま読み込む (CorrectOrbit・MakeBump・CSR.n の一部ほか)。docs/32 |
| und.py | 106 | UND (undulator.f。K64 は磁場なし)。docs/25 |
| measure.py | 66 | tlum・tluma (MEASURE の統計。ffscmd からは呼ばれない) |
| vecfn.py | 56 | Vector / ReleaseVector |
| `__init__.py` / sadkw.py / constants.py | 121 / 34 / 10 | |
| **非推奨**: ffs.py / tracking.py / matching.py / optics.py | 1599 / 1239 / 602 / 471 | 独自実装 (SAD 未照合。tracking.py は放射係数に疑い)。新規コードは使わない。`__init__.py` の export は旧 API のまま |

テストは `pysad/tests/` (`test_*.py` 60 ファイル、`def test` 303 個 (パラメータ化で 523 件に展開))。

## 3. SAD 準拠状況

| 領域 | 状態 | 直訳元 |
|---|---|---|
| CALC (4D 閉軌道・Twiss・結合・オフモメンタム・INS・不安定ラインの続行) | 一致 (オンモメンタム 1e-11、COD/結合経由は相対 1e-7) | qcod, qcell1, qtwiss, qmdiag |
| 要素の COD 写像 (四極フリンジ、薄い多極、BEND の K0/K1/ロール/F1、ANGLE=0 BEND、MULT、ROTATE/ずれ付きの全要素) | 一致 | tquade, tthine, tchge, tbende, tsteee, tqlfre |
| CALC6D (CELL/INS)、幾何 (GX.., GEO, ORG)、RADINT | 一致 (1.7e-10、4.8e-13、7 桁) | qcell6d, tfgeo, INTGRL |
| TRPT の加速と MARK OFFSET (線形加速器) | 一致 (3546 行の表が印字精度で同じ) | tfacc, tcave/tmulte acc, qcell |
| Emittance[]/EMIT (RFSW × RADCOD、空洞 L≠0、MULT、不安定な光学、報告の全行、結果変数) | 一致 (≤1e-11〜1e-9) | temit, tturne, tradke, tecalc, teigen |
| INTRA / Touschek、Emittance[] の拡張 | 一致 (≤7e-12、2.1e-10) | tintrb, emit.n |
| 粒子追跡 (NORAD/RAD/FLUC、K1≠0 BEND、MULT、SOL、粒子損失) | ≤5e-11、乱数の最終状態まで一致 (損失のない場合)。損失周は論理一致、実寸の DA の最終 seed はカオスで不一致 | tturn1 系, photons, tapert |
| マッチング FIT/FREE/FIX/GO (解析微分、DESCEND、fuzz、tsvd、2 段目、INS、フック) | 反復列・残差・K1 を SAD と同じ桁で再現 (12 桁) | tffsmatch, qdcell, tsolvg |
| FITP>1 (idp ループ、重み、fam 点) | 直訳のまま SAD と未照合 (SAD (WSL) が GO で再起動するため参照値なし)。1 段目の解は一致 | tffscalc, tffamsetup |
| DynamicApertureSurvey | NORAD は地図・Score 完全一致。RAD+FLUC の実寸は Score が 1〜3 % 違う (乱数のカオス) | trackd, tfdapert |
| 誤差コマンド (DELX.. の実現、SeedRandom) | 誤差の実現 (乱数列) は一致。補正量は SAD の teigen の縮退で丸め依存 | terror, tfadjst |
| 光学関数 (TransferMatrix, BeamMatrix, NormalCoordinates, OrbitGeo, GeoBase, SymplecticJ, SurvivedParticles, AccelerateParticles) | 一致 (1e-9) | Optics.n, emit.n, dapert.n |

| COORD・TCAVI・UND・MULT ANGLE≠0・MAP・DELCOR | gs と一致 (追跡 0〜10 桁、Twiss 1e-11〜1e-13、Emittance[] 1e-11〜5e-10)。TCAVI の Emittance[] は L=0 のみ、UND は磁場なし (K64 が設定不能)、DELCOR の CAVI・MAP の GEO は未対応 | tcoord, ttcav, undulator, tmulta, temap, tcorr |
| SOL (Twiss・Emittance[]・追跡・幾何) | Twiss 1e-11、Emittance[] εx 2e-12〜4e-9。追跡は SOL (BZ≠0) 内の四極 + RAD で現在の gs が保存参照と食い違う。現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い (B68。ソースの箇所は未特定)。領域内の BEND・MULT・CAVI は未対応。GEO=1、1.7〜3 GeV、BZ 0.3〜5 T の範囲だけで比較 | tsole1, tsolqu, tsgeo |
| POL・スピン (Emittance[]・EMIT・TrackParticles) | gs と 4e-12 以下 (追跡)、偏極時間 13 桁 (ずれ・傾きのあるリング)。`Polarization` は K64 が NaN。SOL 領域内の MULT・BEND・CAVI のスピンは未確認 | photons.f tradke・sprot, spnorm |
| CSR (CSRTrack・CSRHaissin・CSRSetupOY・CSRMatrix・CSRConvert・CSRInit、CSR.n の MakeZL ほか) | gs と 11 桁 (CSRTrack・CSRHaissin)〜1e-10 (CSRMatrix・CSRInit・MakeZL)。CSROYMatrix は未照合 (gs が wrongtype) | csrtrack.f, csroy.f, csrtest.f |
| WAKE (追跡・光学)・BEAMBEAM・RADLIGHT | gs と一致 (RADLIGHT の軌道は差 0)。BEAMBEAM は光学に効かず恒等。RADLIGHT は K64 が BEND の中を記録しない | txwake, tffswake, bbstrhl, tlinit |
| 特殊関数・Fit・FindRoot・SAD パッケージ・FFS の小コマンド・GaussianCoulomb・BBBrem・SynchroBetaEmittance | gs と一致 (特殊関数 458 式、Fit 5e-10〜1e-9、SynchroBeta emitx 4e-7〜1e-5) | gamma.f ほか, Fit.n, temits.f |
| GUI (Tcl/Tk)・描画 (ListPlot・Plot) | SAD 自身の `.n` を pysad のインタプリタで tkinter に表示。OpticsPlot・LatticePlot は枠・矩形・線が gs と一致。文字の数はフォント依存 | tfTkInter_.c |

MON/MCOD/TCOD はこの SAD ビルドでは何もしない (式として評価されるだけ) ので pysad も受理のみ。

## 4. 主な仕様 (直訳で再現している SAD の扱い)

- **要素値**: `LINE["DX","QD.1"]=v` は CALC で全 QD に展開される (ファミリーの代表)。1 台だけの誤差は `DELX σ QD.1`。pysad API では `set_element` (ファミリー) と `set_error` (1 台) を分ける。
- **INS の初期軌道**: 先頭 MARK の DX/DPX/DY/DPY/DZ/DDP が初期座標 (誤差ではない)。
- **マッチング**: EMIT が失敗して EMITX/EMITY=NaN のときは fuzz の MAX/MIN が NaN を無視する gfortran の挙動を再現し、emminv 由来の ~1e-17 で EX/EPX の重みが 1e5 倍になる (SAD と同じ解に着く条件)。fmincube の定数はオイラー定数 γ。
- **色収差**: SAD の CHRO は BEND 成分を含まない。色収差は DP0 スキャン + PolynomialFit で求める (Fit.n の式をそのまま評価)。
- **EMIT コマンド**は GAUSS フラグを立てる (Emittance[] は変えない)。SADScript の `1GEV` は SAD では `1*GEV`=1 だが pysad は単位として読む (意図的な拡張)。
- **FFS 内の ON/OFF の語**は無効 (`INTRA;` / `NOINTRA;` と書く)。MAIN 部の ON/OFF は有効。
- **乱数**: 既定の "SAD" 乱数 (LCG + Box–Muller) を int64 で写し、FLUC の引く順 (要素外側・粒子内側) まで同じ。`SeedRandom[]` の戻りは 5 要素。
- **周長と運動量**: K64 は SOL (BOUND≠0) の長さを周長に入れない (`emit6d._emit_circ`)。MOMENTUM の初期値は 30 GeV (K64 initb1.f。以前は 1 GeV と書いていた)。
- その他の癖は `11_sad_quirks.md` (A1〜A51、B1〜B7、C1〜C28、D1〜D5)。K64 自体の誤りは `19_sad_bug_report.md` (B01〜B62)。

## 5. 性能

- (2026-11-06 の測定値で、再測定していない) 1 回の cell calc 0.03 s (不変時) / 0.145 s (変更後)、Twiss[] 0.6 ms、追跡 176 ns/粒子・要素 (numba)。
- 600 s 内に終わらない計算: 誤差入り 6D emit × 多数の実現 (300 実現 221 s、SAD 37 s)、DA の 5000 周 × 969 粒子 (約 50 分、SAD も 300 s 超)。
  残る改善点は 6D の tturne 1 周 (0.03 s) と DA 系 (`history/work_perf2.md`, `history/work_perf3.md`)。

## 6. 未対応・制限 (優先順は `history/00_RESUME.md` に従う。一覧は `14_fortran_status.md` §9)

- 性能 (上記)。FITP>1 の最終解の照合 (SAD が再起動)。誤差系の要調査 2 件 (docs/12 §4)。
- 不安定ラティスの 1e199 級の値や NaN は再現しない。実寸の DA の最終 seed (カオス)。
- 要素: SOL 領域内の BEND・MULT・CAVI (BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げた)、SOL の DX DY CHI・F1 (Twiss・Emittance[]・追跡)、他の要素の光子 (PHOTONS の tphrec)、
  DELCOR の CAVI、MAP の GEO・CompiledMap、MULT ANGLE≠0 の Twiss の BZ・RF・F1・SK0、TCAVI の RANKICK。COORD・TCAVI・UND・MULT ANGLE・MAP・DELCOR・SOL・POL (Emittance[]・追跡のスピン) は実装済み (docs/23〜29、34、41)。PHSROT は K64 で動作しない。
- 未解決の差 (BEND・CAVI・MULT は原因未調査、`history/work_unimpl.md`): BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)、L≠0 の CAVI に DX・DY で 1e-3、
  SOL (BZ≠0) 内の四極 + RAD の追跡で現在の gs (2026-10-19 ビルド) が保存済みの参照 (bench_level9、10-14) と別の値を返す (現在の gs が RAD ON・BZ≠0 の領域内の四極の力を入れていない疑いが強い。2026-11-16 の単一粒子追跡で確認。B68)、
  K3 が強くずれの大きい MULT リングで CALC の閉軌道が gs と違う解になり Emittance が NaN。
- 関数: SynchroBetaEmittance (`synchrobeta.py`、12 セルと 24 セルの平面リングで照合。結合のあるラティスは gs が SIGSEGV で比較不能)、CSRInit・CSRMatrix mode 1000・CSRConvert・CSR.n (docs/34、38)、RADLIGHT の軌道の記録と RadiationField・RadiationSpectrum (docs/43)、GaussianCoulomb(U・Fitted)・BBBrem1・特殊関数・Fit/FindRoot・Packages/*.n は実装済み (docs/30〜33)。
  MessageName の関数形、複素 FindRoot、PolyLog の |z|≥150 の非整数次は未対応。
- 未翻訳: MEASURE/TRCOD (gs が動かない: B42)、MAIN の TRACK の Fourier 解析 (wjfit・wfres・wfbin・wiord・wjdraw)、`tspac.f` の追跡側。旧 WAKE コマンドは K64 が配布に無い WAKECOMMAND.n を呼ぶので移植不要。
- 対象外: Tk/plot/EPICS/System のうち GUI・描画は SAD 自身の `.n` を解釈して tkinter で表示 (docs/20、35、39、40)。ファイル入出力の失敗の扱い。

## 7. 11/15 夜の追加 (docs/31〜33、history/work_unimpl.md)
- 特殊関数・数学関数 (gs と 458 式が一致)、Fit/FindRoot/NIntegrate/FitGaussian、DownhillSimplex、ファイル関数、SingularValues、RadiationField/Spectrum。
- SAD の `Packages/*.n` を SADScript のまま実行 (CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・記号微分 D・rational・tableform)。ResponseMatrix・LinearSolveConditioned・Master・FFS["VAR"] を翻訳。
- FFS の VARIABLES・ATTRIBUTE・MATRIX・SHOW・VARY・RECOVER・REVERSE・STATUS・PRINT・EXECUTE・OUTPUT/APPEND/CLOSE・INPUT/READ。
- インタプリタを gs に合わせた箇所: Orderless 照合、Flat の積、`-x`、数の書式、`//.`、`Unevaluated$`、Head、Part、VectorQ。K64 の誤り (B43〜B51) は docs/19。
- 試験: 440 通過、6 skip。QUICK 回帰 (50 ベンチ) は FAIL なし。

## 8. 11/16 の追加 (docs/36〜43、history/work_unimpl.md)
- COORD (23)・TCAVI (24)・UND (25)・MULT ANGLE≠0 (26)・SOL の Twiss・Emittance[]・追跡 (27)・MAP (28)・DELCOR (29)・GaussianCoulomb・BBBrem (30)・POL/スピンと CSR (34)・描画 (35)・GUI/Tcl (20)・WAKE (21)・BEAMBEAM (22)。
- 本日: POL の Emittance[]・EMIT のスピン (MULT・空洞・SOL まで) と TrackParticles のスピン追跡 (`track6.py` の `sprot`、gs と 4e-12 以下)、RADLIGHT の軌道 (`ffscmd._RadLight`、gs と差 0)、
  SynchroBetaEmittance、CSRInit・CSRMatrix mode 1000・CSRConvert、tdgeo (DRIFT の L と BEND の ANGLE・L だけ解析微分)、Skip・Read の型、`measure.py` (tlum・tluma)。
- K64 の仕様の発見: SOL の長さは周長に入れない、MOMENTUM の既定は 30 GeV、tmulke の bsir0、B55〜B62 (docs/19)。
- 解説文書: `36_wake_guide.md`、`37_beambeam_guide.md`、`38_csr_guide.md`、`39_tcltk_guide.md`、`40_graphics_guide.md`、`41_pol_guide.md`、`42_sol_guide.md`、`43_radlight_guide.md`。
- 試験: 520 通過・3 skip (523 件収集)。回帰は 50 ベンチをすべて実行して FAIL 0。
- 2026-11-16 夜の修正 (gs との照合で見つけた差): FIT の UNIPOL・RELW (フラグを match6 につないだ。docs/11 A54)、DUMP の DK・DDK 列の書式、MATRIX (名前だけの位置・ラティスに無い名前の読み飛ばし)、LDSolve (`DiagonalMatrix` の複素数・行列演算の虚部 0 の扱い)、OpenRead (開けないファイルで入力行を打ち切る。A53)、tapert (セッション最初の TrackParticles の l_track=0。A52)。
- 未翻訳で確認した点: SPCH 要素は K64 の tturn.f:521 が PSPAC フラグで tpspac を呼ぶ空間電荷の要素で pysad は未翻訳 (PROF は K64 に無い)。LOSSMAP は gs が TrackParticles の戻りに 2 行 (落ちた要素と周) を足すが pysad は未実装。BIPOL の効果は FIT だけ。
