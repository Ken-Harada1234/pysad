# SADHelp (1.4.4.2k64 Updated: 10/19/2026) の項目の目的別分類と翻訳方針

ユーザー方針 (2026-11-06): **直接翻訳するのはソースが Fortran/C で内容が物理計算のもの**。**ソースが .n のものは Python で書き直す** (直訳ではない)。**ファイル操作・リスト操作・GUI・描画は Python のライブラリ (open, list/numpy, matplotlib 等) に置き換える**ので翻訳対象外 (必要なら Python で書き直す)。

分類と SAD 側・pysad の列は docs/14_pysad_vs_sadhelp.md と同じ機械照合 (生成: `<work-perf>/sadhelp_purpose.py`、分類は手で割り当て)。

## 方針の一覧

| 分類 | 方針 | 具体的に |
|---|---|---|
| 物理計算 (26) | 直接翻訳 | Fortran/C の式・分割・反復順序を写す (pysad の本体)。ソースが .n のものは Python で書き直す |
| マッチング・最適化の仕組み (21) | 直接翻訳 | tffsmatch/tffscalc/qdcell の直訳 (フックの呼び出しを含む)。.n のもの (OptimizeOptics, DownhillSimplex) は Python で書き直す |
| 要素・ビームラインの操作 (7) | Python で書き直し | pysad のデータ構造 (Element, Lattice, Env) の上で同じ意味になるよう書く。直訳ではない |
| 乱数 (6) | 直接翻訳 | FLUC の照合に SAD と同じ乱数列が要るので C (tfRandom_.c) を写す (sadrandom.py) |
| 数値計算・フィット (6) | Python ライブラリで代替 | numpy/scipy (curve_fit, integrate, interpolate, optimize) に置き換える。翻訳対象外、必要なら Python で書き直す |
| リスト・関数型 (言語機能) (21) | Python で書き直し (インタプリタ) | ユーザースクリプトを走らせるために sadscript.py が Python のリスト演算で実装する。pysad を Python から直接使うときは list/numpy |
| 文字列・式 (言語機能) (13) | Python で書き直し (インタプリタ) | sadscript.py が Python の str で実装。直接使うときは str |
| ファイル・入出力 (15) | Python ライブラリで代替 | open/print/csv/json に置き換える。翻訳対象外 (OpenRead の失敗の扱いなどは移植しない) |
| GUI・描画 (14) | Python ライブラリで代替 | matplotlib (plt.plot, scatter, contourf, hist) に置き換える。翻訳対象外 |
| システム・プロセス・並列 (7) | Python ライブラリで代替 | subprocess / multiprocessing / tempfile に置き換える。翻訳対象外 |
| ユーティリティ (時刻・計測・デバッグ) (5) | Python ライブラリで代替 | time / datetime / logging に置き換える。翻訳対象外 |
| 表示 (表・一覧) (9) | Python で書き直し | DISP/TYPE のように SAD と同じ文字列を出すものは照合のため書き直す。それ以外 (SHOW, STATUS 等) は必要になったら |
| 制御・入出力 (FFS コマンド) (19) | Python ライブラリで代替 | 入力の切り替え・停止・中断は Python の制御構造とファイル操作で置き換える。翻訳対象外 |
| 要素編集・変数 (FFS コマンド) (12) | Python で書き直し | FREE/FIX/SAVE/RESET/COUPLE 等は pysad の要素値の仕組みで書く |
| 誤差 (FFS コマンド) (1) | 直接翻訳 | terror.f の乱数の引き方まで写す (照合済み) |
| FFS フラグ・物理パラメータ (変数) (30) | 直接翻訳 | tfsetparam/tffsa が読む意味のとおりに使う |
| 表示・書式 (変数) (4) | Python で書き直し | $FORM/PageWidth/TITLE/CASE は出力の書式。照合に要る $FORM/PageWidth は書き直し済み、描画用は対象外 |
| フック・特殊 (変数) (6) | 直接翻訳 | FitFunction/MatchingAmplitude/InitialOrbits/OpticsProlog 等は SAD が評価する仕組みごと写す |

## 件数 (方針別)

| 方針 | 項目数 | うち pysad 実装・参照あり |
|---|---|---|
| Python ライブラリで代替 | 66 | 24 |
| Python で書き直し | 32 | 19 |
| 直接翻訳 | 77 | 64 |
| Python で書き直し (インタプリタ) | 34 | 27 |
| Python で書き直し (.n) | 13 | 9 |

## 物理計算 (26) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `BEAMSIZE(BEAM)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `CALCULATE(CAL)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `CHROMATICITY(CHRO)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `EMITTANCE(EMIT)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `MEASURE(MEA)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | 直接翻訳 |
| コマンド | `ORG` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `RADINT` | Fortran tffsa.f | 実装 | 直接翻訳 |
| 関数 | `AccelerateParticles` | .n dapert.n | 実装 | Python で書き直し (.n) |
| 関数 | `BeamMatrix` | .n Optics.n | 実装 | Python で書き直し (.n) |
| 関数 | `DynamicApertureSurvey` | .n dapert.n | 実装 | Python で書き直し (.n) |
| 関数 | `Emittance` | .n emit.n | 実装 | Python で書き直し (.n) |
| 関数 | `ExternalMap` | .n init.n | 未対応 | Python で書き直し (.n) |
| 関数 | `GaussianCoulomb` | Fortran tfinitn.f | 未対応 | 直接翻訳 |
| 関数 | `GeoBase` | .n Optics.n | 実装 | Python で書き直し (.n) |
| 関数 | `LINE` | Fortran tfinitn.f; .n OpticsPlot.n,OpticsPlot_230724a.n,correction.n,… | 実装 | 直接翻訳 |
| 関数 | `OrbitGeo` | .n Optics.n | 実装 | Python で書き直し (.n) |
| 関数 | `RadiationField` | Fortran tfinitn.f | 未対応 | 直接翻訳 |
| 関数 | `RadiationSpectrum` | Fortran tfinitn.f | 未対応 | 直接翻訳 |
| 関数 | `SurvivedParticles` | Fortran tfinitn.f | 実装 | 直接翻訳 |
| 関数 | `SymplecticJ` | .n emit.n | 実装 | Python で書き直し (.n) |
| 関数 | `SynchroBetaEmittance` | .n emit.n | 未対応 | Python で書き直し (.n) |
| 関数 | `TouschekLifetime` | .n emit.n | 実装 | Python で書き直し (.n) |
| 関数 | `TrackParticles` | Fortran tfinitn.f | 実装 | 直接翻訳 |
| 関数 | `TransferMatrix` | .n Optics.n | 実装 | Python で書き直し (.n) |
| 関数 | `Twiss` | Fortran tfinitn.f; .n correction.n,corutility.n,optimize.n | 実装 | 直接翻訳 |
| 関数 | `WakeFunction` | Fortran? tffswake.f | 未対応 | 直接翻訳 |

## マッチング・最適化の仕組み (21) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `FIT` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `FITPOINTS(FITP)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `FIX` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `FREE` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `GO` | Fortran tffsa.f | 実装 | 直接翻訳 |
| コマンド | `REJECT(REJ)` | Fortran tffsa.f | 実装 | 直接翻訳 |
| 関数 | `FFS` | Fortran tfinitn.f | 実装 | 直接翻訳 |
| 関数 | `FFS$SHOW` | Fortran tfinitn.f | 実装 | 直接翻訳 |
| 関数 | `FitValue` | Fortran? qdcell.f,tdfun.f,tffsa.f | 実装 (ユーザー定義のフック (tfgetfitval) を pysad が呼ぶ) | 直接翻訳 |
| 関数 | `FitWeight` | Fortran? tffscalc.f | 実装 (ユーザー定義のフック (twfit) を pysad が呼ぶ) | 直接翻訳 |
| 関数 | `OptimizeOptics` | .n optimize.n | 未対応 | Python で書き直し (.n) |
| 関数 | `VariableRange` | Fortran? tffsmatch.f,tfvars.f | 文字列参照あり | 直接翻訳 |
| 関数 | `VariableWeight` | Fortran? tffsmatch.f | 文字列参照あり | 直接翻訳 |
| 関数 | `DownhillSimplex` | .n optimize.n | 未対応 | Python で書き直し (.n) |
| 変数 | `CONVERGENCE` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `ExponentOfResidual` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `FFS$NumericalDerivative` | Fortran tffsmatch.f | 参照あり (ffscmd.py, sadscript.py) | 直接翻訳 |
| 変数 | `MatchingResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `NetResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `OffMomentumWeight` | Fortran tffsmatch.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `UnstableLevel` | Fortran tffscalc.f | 参照あり (ffscmd.py) | 直接翻訳 |

## 要素・ビームラインの操作 (7) — 方針: Python で書き直し

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `Element` | .n element.n | 実装 | Python で書き直し |
| 関数 | `SetElement` | .n element.n | 実装 | Python で書き直し |
| 関数 (言語機能の章) | `BeamLine` | Fortran (tfinitn.f) | 実装 | Python で書き直し |
| 関数 (言語機能の章) | `BeamLineName` | Fortran (tfinitn.f) | 実装 | Python で書き直し |
| 関数 (言語機能の章) | `ExtractBeamLine` | Fortran (tfinitn.f) | 実装 | Python で書き直し |
| 関数 (言語機能の章) | `PrintBeamLine` | Fortran (tfinitn.f) | 実装 | Python で書き直し |
| 関数 (言語機能の章) | `WriteBeamLine` | Fortran (tfinitn.f) | 実装 | Python で書き直し |

## 乱数 (6) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `SEED` | Fortran? initb1.f,tfinitn.f | 未対応 | 直接翻訳 |
| 関数 | `Random` | C tfRandom_.c | 実装 | 直接翻訳 |
| 関数 | `GaussRandom` | C tfRandom_.c | 実装 | 直接翻訳 |
| 関数 | `ParabolaRandom` | C tfRandom_.c | 実装 | 直接翻訳 |
| 関数 | `SeedRandom` | C tfRandom_.c | 実装 | 直接翻訳 |
| 関数 | `ListRandom` | C tfRandom_.c | 未対応 | 直接翻訳 |

## 数値計算・フィット (6) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `Fit` | .n Fit.n | 未対応 | Python ライブラリで代替 |
| 関数 | `FitEmit` | .n Fit.n | 未対応 | Python ライブラリで代替 |
| 関数 | `FitGaussian` | .n Fit.n | 未対応 | Python ライブラリで代替 |
| 関数 | `NIntegrate` | .n D.n,OpticsPlot_230724b.n | 未対応 | Python ライブラリで代替 |
| 関数 | `PolynomialFit` | .n Fit.n | 未対応 | Python ライブラリで代替 |
| 関数 | `Spline` | .n D.n,OpticsPlot_230724b.n | 実装 | Python ライブラリで代替 |

## リスト・関数型 (言語機能) (21) — 方針: Python で書き直し (インタプリタ)

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `Apply` | Fortran tfinitn.f | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `Cases` | Fortran tfinitn.f | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `DeleteCases` | Fortran tfinitn.f | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `Difference` | Fortran tfinitn.f | 未対応 | Python で書き直し (インタプリタ) |
| 関数 | `FixedPoint` | .n Misc.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `FixedPointList` | .n Misc.n | 未対応 | Python で書き直し (インタプリタ) |
| 関数 | `Map` | Fortran tfinitn.f; .n correction.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `MapThread` | Fortran tfinitn.f; .n Canvas.n,bump.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `Position` | Fortran tfinitn.f; .n functions.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `Scan` | Fortran tfinitn.f; .n Canvas.n,Canvas3D.n,CheckEntryList.n,… | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `ScanThread` | Fortran tfinitn.f; .n correction.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `SelectCases` | Fortran tfinitn.f | 未対応 | Python で書き直し (インタプリタ) |
| 関数 | `SwitchCases` | Fortran tfinitn.f | 未対応 | Python で書き直し (インタプリタ) |
| 関数 | `Thread` | Fortran tfinitn.f | 実装 | Python で書き直し (インタプリタ) |
| 関数 | `Class` | .n Class.n | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `Extract` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `Head` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `Length` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `List` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `Part` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `MatchQ` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |

## 文字列・式 (言語機能) (13) — 方針: Python で書き直し (インタプリタ)

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 (言語機能の章) | `FromCharacterCode` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringFill` | Fortran (tfinitn.f) | 未対応 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringJoin` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringLength` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringMatchQ` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringPart` | Fortran (tfinitn.f) | 未対応 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringPosition` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringReplace` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `StringTrim` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `Symbol` | Fortran (tfinitn.f) | 未対応 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `ToCharacterCode` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `ToExpression` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |
| 関数 (言語機能の章) | `ToString` | Fortran (tfinitn.f) | 実装 | Python で書き直し (インタプリタ) |

## ファイル・入出力 (15) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `$FORM` | .n KBFrame.n,KBMainFrame.n,Plots.n,… | 文字列参照あり | Python ライブラリで代替 |
| 関数 | `$Input` | .n init.n | 未対応 | Python ライブラリで代替 |
| 関数 | `$Output` | .n init.n | 未対応 | Python ライブラリで代替 |
| 関数 | `Close` | .n KBFrame.n,init.n | 実装 | Python ライブラリで代替 |
| 関数 | `OpenAppend` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `OpenRead` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `OpenWrite` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `PageWidth` | .n latex.n | 文字列参照あり | Python ライブラリで代替 |
| 関数 | `Print` | Fortran tfinitn.f; .n Class.n | 実装 | Python ライブラリで代替 |
| 関数 | `Read` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `ReadString` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `StandardForm` | Fortran tfinitn.f; .n KBMainFrame.n | 未対応 | Python ライブラリで代替 |
| 関数 | `StringToStream` | Fortran tfinitn.f | 未対応 | Python ライブラリで代替 |
| 関数 | `Write` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `WriteString` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |

## GUI・描画 (14) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `DRAW` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python ライブラリで代替 |
| 関数 | `BeamPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | 未対応 | Python ライブラリで代替 |
| 関数 | `ColumnPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `FitPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `GeometryPlot` | .n GeometryPlot.n | 未対応 | Python ライブラリで代替 |
| 関数 | `HistoPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `LaTeX` | .n latex.n | 未対応 | Python ライブラリで代替 |
| 関数 | `ListContourPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `ListDensityPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `ListPlot` | .n ListPlot.n | 未対応 | Python ライブラリで代替 |
| 関数 | `OpticsPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | 未対応 | Python ライブラリで代替 |
| 関数 | `ParametricPlot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `Plot` | .n Plots.n | 未対応 | Python ライブラリで代替 |
| 関数 | `PointSymbol` | .n? CSR.n,Canvas.n,ListPlot.n,… | 未対応 | Python ライブラリで代替 |

## システム・プロセス・並列 (7) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `Fork` | C tfProcess_.c | 未対応 | Python ライブラリで代替 |
| 関数 | `OpenShared` | .n shared.n | 未対応 | Python ライブラリで代替 |
| 関数 | `Shared` | .n CSR.n,shared.n | 未対応 | Python ライブラリで代替 |
| 関数 | `SharedSize` | Fortran tfinitn.f | 未対応 | Python ライブラリで代替 |
| 関数 | `System` | C tfProcess_.c | 実装 | Python ライブラリで代替 |
| 関数 | `TemporaryName` | Fortran tfinitn.f | 未対応 | Python ライブラリで代替 |
| 変数 | `NPARA` | Fortran initb1.f,tfefun1.f,tfsetparam.f,…; .n init.kekb.n,init.n,init.skekb.n | 参照あり (ffscmd.py, parser.py) | Python ライブラリで代替 |

## ユーティリティ (時刻・計測・デバッグ) (5) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 関数 | `DateString` | .n init.n | 実装 | Python ライブラリで代替 |
| 関数 | `MemoryCheck` | Fortran tfinitn.f | 未対応 | Python ライブラリで代替 |
| 関数 | `TimeUsed` | Fortran tfinitn.f | 実装 | Python ライブラリで代替 |
| 関数 | `Timing` | .n functions.n | 実装 | Python ライブラリで代替 |
| 関数 | `TracePrint` | Fortran tfinitn.f | 未対応 | Python ライブラリで代替 |

## 表示 (表・一覧) (9) — 方針: Python で書き直し

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `ATTRIBUTE(ATTR)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `DISPLAY(DISP)` | Fortran tffsa.f,tracke.f | 実装 | Python で書き直し |
| コマンド | `MATRIX(MAT)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `PRINT(PRI)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `REFERENCE(REF)` | Fortran tfdisp.f,tffsa.f | 実装 | Python で書き直し |
| コマンド | `SHOW` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `STATUS(STAT)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `TYPE(T)` | Fortran tffsa.f | 実装 | Python で書き直し |
| コマンド | `VARIABLES(VAR)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |

## 制御・入出力 (FFS コマンド) (19) — 方針: Python ライブラリで代替

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `ABORT` | Fortran tffsa.f | 実装 | Python ライブラリで代替 |
| コマンド | `APPEND(APP)` | Fortran phdrw.f,pstati.f,tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `CLOSE(CLO)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `ELSE` | Fortran tfif.f | 実装 | Python ライブラリで代替 |
| コマンド | `ELSEIF` | Fortran tfif.f | 未対応 | Python ライブラリで代替 |
| コマンド | `END` | Fortran tffile.f | 実装 | Python ライブラリで代替 |
| コマンド | `ENDIF` | Fortran tfif.f | 実装 | Python ライブラリで代替 |
| コマンド | `EXECUTE(EXEC)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `IF` | Fortran tfif.f | 実装 | Python ライブラリで代替 |
| コマンド | `INPUT(IN)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `OUTPUT(OUT)` | Fortran phdrw.f,pstati.f,tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `QUIT` | Fortran tffsa.f | 実装 | Python ライブラリで代替 |
| コマンド | `READ` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `REPEAT(REP)` | Fortran tffsa.f | 実装 | Python ライブラリで代替 |
| コマンド | `RESUME(RES)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `STOP` | Fortran pstati.f,tffsa.f | 実装 | Python ライブラリで代替 |
| コマンド | `SUSPEND(SUSP)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `TERMINATE(TERM)` | Fortran tffile.f | 未対応 | Python ライブラリで代替 |
| コマンド | `UNTIL` | Fortran tffsa.f | 実装 | Python ライブラリで代替 |

## 要素編集・変数 (FFS コマンド) (12) — 方針: Python で書き直し

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `BYE` | Fortran tffsa.f | 実装 | Python で書き直し |
| コマンド | `COUPLE(COUP)` | Fortran tffsa.f | 実装 | Python で書き直し |
| コマンド | `EXPAND` | Fortran tffsa.f,tftwiss.f | 実装 | Python で書き直し |
| コマンド | `RECOVER(REC)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `RENUMBER(RENUM)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `RESET` | Fortran tffsa.f | 実装 | Python で書き直し |
| コマンド | `REVERSE(REV)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `SAVE` | Fortran tffsa.f,tracke.f | 実装 | Python で書き直し |
| コマンド | `SPLIT` | Fortran tffsa.f | 未対応 | Python で書き直し |
| コマンド | `USE` | Fortran tffsa.f | 実装 | Python で書き直し |
| コマンド | `VARY` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) | Python で書き直し |
| コマンド | `VISIT` | Fortran tffsa.f | 実装 | Python で書き直し |

## 誤差 (FFS コマンド) (1) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| コマンド | `DUMP` | Fortran terror.f | 実装 | 直接翻訳 |

## FFS フラグ・物理パラメータ (変数) (30) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 変数 | `CHARGE` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `DP` | Fortran tffs.f,tfsetparam.f,tftwiss.f | 参照あり (emit_extra.py, ffscmd.py, match6.py…) | 直接翻訳 |
| 変数 | `DP0` | 不明 | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `DTSYNCH` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EFFRFFREQ` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EFFVC` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EFFVCRATIO` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EMITX` | Fortran initb1.f,temit.f,tffs.f,…; .n emit.n | 参照あり (ffscmd.py, maintrack.py, match6.py…) | 直接翻訳 |
| 変数 | `EMITXE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EMITY` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, match6.py…) | 直接翻訳 |
| 変数 | `EMITYE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `EMITZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, parser.py…) | 直接翻訳 |
| 変数 | `EMITZE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `FSHIFT` | Fortran ActTra.f,initb1.f,temitf.f,… | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `GCUT` | Fortran initb1.f; C tfRandom_.c,tgauss_.c | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `LOSSAMPL` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) | 直接翻訳 |
| 変数 | `LOSSDZ` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) | 直接翻訳 |
| 変数 | `MASS` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, optics_fn.py, parser.py…) | 直接翻訳 |
| 変数 | `MINCOUP` | Fortran initb1.f,tfbeam.f,tffs.f,… | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `MOMENTUM` | Fortran ActTra.f,bb.f,bbstrhl.f,…; .n dapert.n | 参照あり (ffs.py, ffscmd.py, optics_fn.py) | 直接翻訳 |
| 変数 | `NBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `NP` | Fortran initb1.f,tfdapert.f,tffsa.f,…; .n dapert.n | 参照あり (ffscmd.py, maintrack.py) | 直接翻訳 |
| 変数 | `OMEGA0` | Fortran ActTra.f,initb1.f,tfgeo.f,… | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `PBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `PhotonList` | Fortran photons.f | 未対応 (参照なし) | 直接翻訳 |
| 変数 | `PHICAV` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `SIGE` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, parser.py…) | 直接翻訳 |
| 変数 | `SIGZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, parser.py) | 直接翻訳 |
| 変数 | `SpeedOfLight` | Fortran tfinitn.f | 参照あり (__init__.py, parser.py, sadscript.py) | 直接翻訳 |
| 変数 | `Track$Excitation` | Fortran tffs.f | 未対応 (参照なし) | 直接翻訳 |

## 表示・書式 (変数) (4) — 方針: Python で書き直し

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 変数 | `$Line` | Fortran tfinitn.f; .n SADTerminal.n,init.n | 参照あり (sadscript.py) | Python で書き直し |
| 変数 | `CASE` | Fortran tffsa.f | 未対応 (参照なし) | Python で書き直し |
| 変数 | `DAPWIDTH` | Fortran initb1.f,tfdapert.f,trackd.f | 参照あり (ffscmd.py) | Python で書き直し |
| 変数 | `TITLE` | Fortran tffsa.f,twsdrw.f | 未対応 (参照なし) | Python で書き直し |

## フック・特殊 (変数) (6) — 方針: 直接翻訳

| 種別 | 名前 | SAD 側 | pysad | 方針 |
|---|---|---|---|---|
| 変数 | `ElementValues` | Fortran tffsa.f | 未対応 (参照なし) | 直接翻訳 |
| 変数 | `FitFunction` | Fortran qdcell.f,tdfun.f,tffsa.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `InitialOrbits` | Fortran tffsa.f | 未対応 (参照なし) | 直接翻訳 |
| 変数 | `MatchingAmplitude` | Fortran tffamsetup.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `OpticsEpilog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) | 直接翻訳 |
| 変数 | `OpticsProlog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) | 直接翻訳 |

