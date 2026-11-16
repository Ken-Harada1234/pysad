# pysad と SADHelp (1.4.4.2k64 Updated: 10/19/2026) の対応表 — SAD 側の定義場所つき

最終更新: 2026-11-16 (現状のスナップショット。SADHelp は 1.4.4.2k64 Updated: 10/19/2026 のまま。pysad 側の列は 2026-11-16 の Interp.builtins と FFS の表から再照合)

機械照合 (`perf_work/pysad_vs_sadhelp.py`)。「SAD 側」は名前が登録・定義されているソース: Fortran (src/*.f)、C (src/*.c の REGn)、.n (Packages/*.n の定義)。`?` 付きは登録の形が見つからず文字列として出現するだけのもの。「pysad」は名前の登録の有無 (実装の深さは docs/03・04・12 を参照): コマンドは FFS._COMMANDS/_abbrev、関数は Interp.builtins、変数はソース中の参照、フラグは _FLAGS_DEFAULT/_SINO、要素は parser/sadkw。

**読み方 (2026-11-16)**: 「実装」は名前の登録があるという意味で、gs (K64 の SAD) との一致の度合いは各要素の文書 (docs/20〜35 と 36〜43) にある。
中身が簡単な代用のものは行に「(限定: …)」と書いた。`Packages/*.n` の関数 (Plot・OptimizeOptics・CorrectOrbit など) は、組み込みの登録ではなく
SAD 自身の .n を pysad のインタプリタで動かす (環境変数 `PYSAD_SAD_PACKAGES` が要る。docs/32)。元の機械照合は組み込みの登録名と FFS の表だけを見るので、
そうした行とユーザー定義フックの行は手で直してある (この表を `perf_work/pysad_vs_sadhelp.py` で作り直すと手で直した行は消える。§6 を見て戻す)。

## 1. FFS コマンド (56): 実装 52、未対応 2、無視 2

| コマンド | SAD 側 | pysad |
|---|---|---|
| `ABORT` | Fortran tffsa.f | 実装 |
| `APPEND(APP)` | Fortran phdrw.f,pstati.f,tffile.f | 実装 |
| `ATTRIBUTE(ATTR)` | Fortran tffsa.f | 実装 |
| `BEAMSIZE(BEAM)` | Fortran tffsa.f | 実装 |
| `BYE` | Fortran tffsa.f | 実装 |
| `CALCULATE(CAL)` | Fortran tffsa.f | 実装 |
| `CHROMATICITY(CHRO)` | Fortran tffsa.f | 実装 |
| `CLOSE(CLO)` | Fortran tffile.f | 実装 |
| `COUPLE(COUP)` | Fortran tffsa.f | 実装 |
| `DISPLAY(DISP)` | Fortran tffsa.f,tracke.f | 実装 |
| `DRAW` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `DUMP` | Fortran terror.f | 実装 |
| `ELSE` | Fortran tfif.f | 実装 |
| `ELSEIF` | Fortran tfif.f | 実装 (sadscript.py の制御文 IF/ELSEIF/ELSE/ENDIF。対応する IF がなければ gs と同じ `???-FFS-Error-?ELSEIF without IF`) |
| `EMITTANCE(EMIT)` | Fortran tffsa.f | 実装 |
| `END` | Fortran tffile.f | 実装 |
| `ENDIF` | Fortran tfif.f | 実装 |
| `EXECUTE(EXEC)` | Fortran tffile.f | 実装 |
| `EXPAND` | Fortran tffsa.f,tftwiss.f | 実装 |
| `FIT` | Fortran tffsa.f | 実装 |
| `FITPOINTS(FITP)` | Fortran tffsa.f | 実装 |
| `FIX` | Fortran tffsa.f | 実装 |
| `FREE` | Fortran tffsa.f | 実装 |
| `GO` | Fortran tffsa.f | 実装 |
| `IF` | Fortran tfif.f | 実装 |
| `INPUT(IN)` | Fortran tffile.f | 実装 |
| `MATRIX(MAT)` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `MEASURE(MEA)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示)。tlum / tluma だけ `pysad/measure.py` に直訳。ほかの MEASURE/TRCOD の経路は未翻訳 (gs が動かない B42。docs/18) |
| `ORG` | Fortran tffsa.f | 実装 |
| `OUTPUT(OUT)` | Fortran phdrw.f,pstati.f,tffile.f | 実装 |
| `PRINT(PRI)` | Fortran tffsa.f | 実装 |
| `QUIT` | Fortran tffsa.f | 実装 |
| `RADINT` | Fortran tffsa.f | 実装 |
| `READ` | Fortran tffile.f | 実装 (ファイルの文を続けて実行。関数 `Read` の型は §2) |
| `RECOVER(REC)` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `REFERENCE(REF)` | Fortran tfdisp.f,tffsa.f | 実装 |
| `REJECT(REJ)` | Fortran tffsa.f | 実装 |
| `RENUMBER(RENUM)` | Fortran tffsa.f | 未対応 (何もしない。docs/33) |
| `REPEAT(REP)` | Fortran tffsa.f | 実装 |
| `RESET` | Fortran tffsa.f | 実装 |
| `RESUME(RES)` | Fortran tffile.f | 実装 |
| `REVERSE(REV)` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `SAVE` | Fortran tffsa.f,tracke.f | 実装 |
| `SEED` | Fortran? initb1.f,tfinitn.f | 未対応 (K64 でも obsolete。docs/33) |
| `SHOW` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `SPLIT` | Fortran tffsa.f | 実装 |
| `STATUS(STAT)` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `STOP` | Fortran pstati.f,tffsa.f | 実装 |
| `SUSPEND(SUSP)` | Fortran tffile.f | 実装 |
| `TERMINATE(TERM)` | Fortran tffile.f | 実装 |
| `TYPE(T)` | Fortran tffsa.f | 実装 |
| `UNTIL` | Fortran tffsa.f | 実装 |
| `USE` | Fortran tffsa.f | 実装 |
| `VARIABLES(VAR)` | Fortran tffsa.f | 実装 |
| `VARY` | Fortran tffsa.f | 実装 (pysad/ffsextra.py、docs/33。表の書式は gs と同じ) |
| `VISIT` | Fortran tffsa.f | 実装 |

**未対応**: `RENUMBER(RENUM)`、`SEED` (K64 でも obsolete)  
**無視 (読み飛ばし)**: `DRAW`、`MEASURE(MEA)` (SADHelp の表にない `GEO`・`TRCOD` なども警告して無視)

## 2. SADScript / FFS の関数 (94、pysad 実装 86)

### Data-Manipulation (6、pysad 実装 6)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Fit` | .n Fit.n | 実装 (pysad/fitroot.py。導関数は数値微分。docs/31) |
| `FitEmit` | .n Fit.n | 実装 |
| `FitGaussian` | .n Fit.n | 実装 |
| `NIntegrate` | .n D.n,OpticsPlot_230724b.n | 実装 |
| `PolynomialFit` | .n Fit.n | 実装 (SAD の Fit.n の定義を評価し、内部で `PolynomialFit$numpy`。xa=0 のとき xl=1 になる SAD の癖もそのまま) |
| `Spline` | .n D.n,OpticsPlot_230724b.n | 実装 (自然スプライン・`Spline[data,0,0]` の形が gs と全桁一致。後者は gs も pysad も未評価のまま。11-16) |

### functional-operations (14、pysad 実装 14)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Apply` | Fortran tfinitn.f | 実装 |
| `Cases` | Fortran tfinitn.f | 実装 |
| `DeleteCases` | Fortran tfinitn.f | 実装 |
| `Difference` | Fortran tfinitn.f | 実装 |
| `FixedPoint` | .n Misc.n | 実装 |
| `FixedPointList` | .n Misc.n | 実装 |
| `Map` | Fortran tfinitn.f; .n correction.n | 実装 |
| `MapThread` | Fortran tfinitn.f; .n Canvas.n,bump.n | 実装 |
| `Position` | Fortran tfinitn.f; .n functions.n | 実装 |
| `Scan` | Fortran tfinitn.f; .n Canvas.n,Canvas3D.n,CheckEntryList.n,… | 実装 (Sort[list, Greater] と共に gs と全桁一致。11-16) |
| `ScanThread` | Fortran tfinitn.f; .n correction.n | 実装 |
| `SelectCases` | Fortran tfinitn.f | 実装 |
| `SwitchCases` | Fortran tfinitn.f | 実装 |
| `Thread` | Fortran tfinitn.f | 実装 |

### FFS-dedicated-functions (28、pysad 実装 28)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `AccelerateParticles` | .n dapert.n | 実装 |
| `BeamMatrix` | .n Optics.n | 実装 |
| `DynamicApertureSurvey` | .n dapert.n | 実装 |
| `Element` | .n element.n | 実装 |
| `Emittance` | .n emit.n | 実装 |
| `ExternalMap` | .n init.n | 実装 (MAP 要素が呼ぶユーザー定義のフック `ExternalMap["TRACK"/"EMIT"/"OPTICS"]` を pysad が評価。docs/28) |
| `FFS` | Fortran tfinitn.f | 実装 |
| `FFS$SHOW` | Fortran tfinitn.f | 実装 |
| `FitValue` | Fortran? qdcell.f,tdfun.f,tffsa.f | 実装 (ユーザー定義のフック (tfgetfitval) を pysad が呼ぶ) |
| `FitWeight` | Fortran? tffscalc.f | 実装 (ユーザー定義のフック (twfit) を pysad が呼ぶ) |
| `GaussianCoulomb` | Fortran tfinitn.f | 実装 (`GaussianCoulombFitted`・`GaussianCoulombU` も。pysad/gausscoul.py、docs/30) |
| `GeoBase` | .n Optics.n | 実装 |
| `LINE` | Fortran tfinitn.f; .n OpticsPlot.n,OpticsPlot_230724a.n,correction.n,… | 実装 |
| `OptimizeOptics` | .n optimize.n | 実装 (SAD の optimize.n を pysad のインタプリタで動かす。K1・チューンは gs と 1e-10。docs/32) |
| `OrbitGeo` | .n Optics.n | 実装 |
| `RadiationField` | Fortran tfinitn.f | 実装 (pysad/radfield.py。BEND の中を記録しない K64 の仕様のまま。gs の RadiationField は SIGSEGV (docs/19 B49) なので値の照合は等速運動の場との比較だけ。docs/34) |
| `RadiationSpectrum` | Fortran tfinitn.f | 実装 (pysad/radfield.py。docs/34) |
| `SetElement` | .n element.n | 実装 |
| `SurvivedParticles` | Fortran tfinitn.f | 実装 |
| `SymplecticJ` | .n emit.n | 実装 |
| `SynchroBetaEmittance` | .n emit.n | 実装 (pysad/synchrobeta.py。gs と相対 1e-5、24 セルの平面リング + CAVI で emitx が 1.4e-7〜1.4e-6。結合のあるラティスは gs が SIGSEGV で比較不能。組み込みの登録ではなく emit.n の呼び出しを置き換える) |
| `TouschekLifetime` | .n emit.n | 実装 |
| `TrackParticles` | Fortran tfinitn.f | 実装 (POL かつ RAD ON でスピン追跡、RADLIGHT で `{結果, 軌道{t,x,y,z}}` を返す。docs/34) |
| `TransferMatrix` | .n Optics.n | 実装 (Optics.n の TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo が gs と表示桁で一致。gs は 1e-16 の雑音の要素を出し pysad は 0。11-16) |
| `Twiss` | Fortran tfinitn.f; .n correction.n,corutility.n,optimize.n | 実装 |
| `VariableRange` | Fortran? tffsmatch.f,tfvars.f | 実装 (ユーザー定義のフック (tffsvarfun) を pysad が呼ぶ) |
| `VariableWeight` | Fortran? tffsmatch.f | 実装 (ユーザー定義のフック (tffsvarfun) を pysad が呼ぶ) |
| `WakeFunction` | Fortran? tffswake.f | 実装 (ユーザー定義の表 `WakeFunction[Longitudinal/Transverse, 要素名]` を pysad が読み、ウェイク追跡と WAKEOPT に使う。docs/21) |

### Graphics (13、pysad 実装 9)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `BeamPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | .n は読み込めるが gs との比較はしていない (未確認。docs/35) |
| `ColumnPlot` | .n Plots.n | .n は読み込めるが動作は未確認 |
| `FitPlot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。43 項目。当てはめの値は 1e-8 相対。docs/35) |
| `GeometryPlot` | .n GeometryPlot.n | 実装 (SAD の .n を動かし tkinter に描く。線 74・多角形 37 が gs と一致。docs/35) |
| `HistoPlot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。38 項目が全部一致。docs/35) |
| `LaTeX` | .n latex.n | 未確認 (`$HaveLaTeX` は gs と同じく `""`) |
| `ListContourPlot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。種類と数は一致、等高線の始点・向きが違う。docs/35) |
| `ListDensityPlot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。全項目が gs と一致。docs/35) |
| `ListPlot` | .n ListPlot.n | 実装 (SAD の .n を動かし tkinter に描く。キャンバスの 43 項目が gs と全部一致。オプション (PlotRange・Scale・PlotLabel・GridLines・Frame・Legend・FontScale・Prolog・Epilog・Background・Plot->False・AspectRatio) を与えた戻りの Graphics 式も gs と文字列で完全に一致 (10 例。窓の見た目ではなく式の比較。`Options[g]` は gs が全オプションを返し pysad は {} の差)。docs/35) |
| `OpticsPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | 実装 (SAD の .n を動かし tkinter に描く。枠・矩形・線は一致、文字の数は Tk のフォントで違う。docs/35) |
| `ParametricPlot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。全項目が gs と一致。docs/35) |
| `Plot` | .n Plots.n | 実装 (SAD の .n を動かし tkinter に描く。全項目が gs と一致。ListPlot と同じオプションを与えた戻りの Graphics 式も gs と文字列で一致。docs/35) |
| `PointSymbol` | .n? CSR.n,Canvas.n,ListPlot.n,… | ListPlot が内部で使う (`CanvasSymbolDirect`)。単独では未確認 |

**未確認 (.n は読めるが gs と比べていないもの。2026-11-16 の確認ログ (local/wake/verify_log.md) にも答えが無く、未確認のまま)**: `BeamPlot`、`ColumnPlot`、`LaTeX`、`PointSymbol`。`CanvasSymbol` (Tcl 文字列を返す版)・`Canvas3D*` は未対応 (docs/35)。

### Input/Output (15、pysad 実装 11)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `$FORM` | .n KBFrame.n,KBMainFrame.n,Plots.n,… | 文字列参照あり |
| `$Input` | .n init.n | 文字列参照あり |
| `$Output` | .n init.n | 文字列参照あり |
| `Close` | .n KBFrame.n,init.n | 実装 |
| `OpenAppend` | Fortran tfinitn.f | 実装 |
| `OpenRead` | Fortran tfinitn.f | 実装 (開けないファイルでは gs と同じに入力行を打ち切る。11-16) |
| `OpenWrite` | Fortran tfinitn.f | 実装 |
| `PageWidth` | .n latex.n | 文字列参照あり |
| `Print` | Fortran tfinitn.f; .n Class.n | 実装 |
| `Read` | Fortran tfinitn.f | 実装 (型は gs に合わせた。docs/31・33) |
| `ReadString` | Fortran tfinitn.f | 実装 |
| `StandardForm` | Fortran tfinitn.f; .n KBMainFrame.n | 実装 (sadscript.py。HoldAll の表示形) |
| `StringToStream` | Fortran tfinitn.f | 実装 |
| `Write` | Fortran tfinitn.f | 実装 |
| `WriteString` | Fortran tfinitn.f | 実装 |

### Multiprocessing (4、pysad 実装 4)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Fork` | C tfProcess_.c | 実装 (限定: 本当には fork せず、自分のプロセス番号を返す) |
| `OpenShared` | .n shared.n | 実装 (限定: プロセス内の辞書。別プロセスとは共有しない) |
| `Shared` | .n CSR.n,shared.n | 実装 (限定: 同上) |
| `SharedSize` | Fortran tfinitn.f | 実装 (限定: 同上。辞書の項目数) |

### Object-oriented-programing (1、pysad 実装 1)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Class` | .n Class.n | 実装 |

### Random-number-functions (5、pysad 実装 5)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Random` | C tfRandom_.c | 実装 |
| `GaussRandom` | C tfRandom_.c | 実装 |
| `ParabolaRandom` | C tfRandom_.c | 実装 |
| `SeedRandom` | C tfRandom_.c | 実装 |
| `ListRandom` | C tfRandom_.c | 実装 (限定: `{"SAD"}` を返すだけ) |

### System-interface (2、pysad 実装 2)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `System` | C tfProcess_.c | 実装 |
| `TemporaryName` | Fortran tfinitn.f | 実装 |

### Utilities (5、pysad 実装 5)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `DateString` | .n init.n | 実装 |
| `MemoryCheck` | Fortran tfinitn.f | 実装 (限定: `{0,0,0,0}` を返すだけ) |
| `TimeUsed` | Fortran tfinitn.f | 実装 |
| `Timing` | .n functions.n | 実装 |
| `TracePrint` | Fortran tfinitn.f | 実装 (限定: 何もしない) |

### (直下) (1、pysad 実装 1)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `DownhillSimplex` | .n optimize.n | 実装 |

## 3. 特別変数 (48、pysad 参照あり 44)

| 変数 | SAD 側 | pysad |
|---|---|---|
| `$Line` | Fortran tfinitn.f; .n SADTerminal.n,init.n | 参照あり (sadscript.py) |
| `CASE` | Fortran tffsa.f | 未対応 (参照なし) |
| `CHARGE` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, ffsextra.py, parser.py…) |
| `CONVERGENCE` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py, ffsextra.py) |
| `DAPWIDTH` | Fortran initb1.f,tfdapert.f,trackd.f | 参照あり (ffscmd.py) |
| `DP` | Fortran tffs.f,tfsetparam.f,tftwiss.f | 参照あり (beambeam.py, emit_extra.py, ffscmd.py…) |
| `DP0` | 不明 | 参照あり (ffscmd.py, parser.py) |
| `DTSYNCH` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFRFFREQ` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFVC` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFVCRATIO` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) |
| `ElementValues` | Fortran tffsa.f | Couple.n (SetCoupledElements) を CALC の前に実行 (gs と一致)。VISIT の退避・BYE の復元は未対応 (docs/32) |
| `EMITX` | Fortran initb1.f,temit.f,tffs.f,…; .n emit.n | 参照あり (ffscmd.py, ffsextra.py, maintrack.py…) |
| `EMITXE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `EMITY` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, ffsextra.py, maintrack.py…) |
| `EMITYE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `EMITZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, ffsextra.py, maintrack.py…) |
| `EMITZE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `ExponentOfResidual` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py) |
| `FFS$NumericalDerivative` | Fortran tffsmatch.f | 参照あり (ffscmd.py, sadscript.py) |
| `FitFunction` | Fortran qdcell.f,tdfun.f,tffsa.f | 参照あり (ffscmd.py) |
| `FSHIFT` | Fortran ActTra.f,initb1.f,temitf.f,… | 参照あり (ffscmd.py, ffsextra.py, parser.py) |
| `GCUT` | Fortran initb1.f; C tfRandom_.c,tgauss_.c | 参照あり (csr.py, ffscmd.py, ffsextra.py) |
| `InitialOrbits` | Fortran tffsa.f | 参照あり (ffscmd.py) |
| `LOSSAMPL` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) |
| `LOSSDZ` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) |
| `MatchingAmplitude` | Fortran tffamsetup.f | 参照あり (ffscmd.py) |
| `MatchingResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) |
| `MASS` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, ffsextra.py, optics_fn.py…) |
| `MINCOUP` | Fortran initb1.f,tfbeam.f,tffs.f,… | 参照あり (ffscmd.py, ffsextra.py, parser.py) |
| `MOMENTUM` | Fortran ActTra.f,bb.f,bbstrhl.f,…; .n dapert.n | 参照あり (ffs.py, ffscmd.py, optics_fn.py)。初期値は K64 の initb1.f で 30 GeV (以前の文書の 1 GeV は誤り。docs/19) |
| `NBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, ffsextra.py, parser.py) |
| `NetResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) |
| `NP` | Fortran initb1.f,tfdapert.f,tffsa.f,…; .n dapert.n | 参照あり (beambeam.py, ffscmd.py, ffsextra.py…) |
| `NPARA` | Fortran initb1.f,tfefun1.f,tfsetparam.f,…; .n init.kekb.n,init.n,init.skekb.n | 参照あり (ffscmd.py, ffsextra.py, parser.py) |
| `OffMomentumWeight` | Fortran tffsmatch.f | 参照あり (ffscmd.py) |
| `OMEGA0` | Fortran ActTra.f,initb1.f,tfgeo.f,… | 参照あり (ffscmd.py) |
| `OpticsEpilog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) |
| `OpticsProlog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) |
| `PBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, ffsextra.py, parser.py) |
| `PhotonList` | Fortran photons.f | 参照あり (ffscmd.py) |
| `PHICAV` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `SIGE` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, ffsextra.py, maintrack.py…) |
| `SIGZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (beambeam.py, ffscmd.py, ffsextra.py…) |
| `SpeedOfLight` | Fortran tfinitn.f | 参照あり (__init__.py, parser.py, sadscript.py) |
| `TITLE` | Fortran tffsa.f,twsdrw.f | 未対応 (参照なし) |
| `Track$Excitation` | Fortran tffs.f | 未対応 (参照なし) |
| `UnstableLevel` | Fortran tffscalc.f | 参照あり (ffscmd.py) |

## 4. フラグ (60、pysad 登録 59)

| フラグ | SAD 側 | pysad |
|---|---|---|
| `ABSW` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `BIPOL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `CALC4D` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `CALC6D` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py, ffsextra.py, optics_fn.py) |
| `CALEXP` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `CELL` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `CMPLOT` | Fortran tffs.f (fname/sino) | 未対応 |
| `COD` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `CODPLOT` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `CONV` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `CONVCASE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `DAMPONLY` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `DAPERT` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `DIFFRES` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `ECHO` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `EMIOUT` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `FFSPRMPT` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `FIXSEED` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `FLUC` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `GAUSS` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py, optics_fn.py) |
| `GEOCAL` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `GEOFIX` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `HALFRES` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `IDEAL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `INS` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `INTRA` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `INTRES` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `JITTER` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `KEEPEXP` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `LOG` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `LOSSMAP` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `LWAKE` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `MOVESEED` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `PHOTONS` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `POL` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py。Emittance・EMIT のスピンと、RAD ON のときの TrackParticles のスピン追跡。docs/34) |
| `PRSVCASE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `PSPAC` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `QUIET` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `RAD` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `RADCOD` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `RADLIGHT` | Fortran tffs.f (fname/sino) | 実装 (TrackParticles が軌道 {t,x,y,z} を返す。K64 は BEND の中を記録せず横方向の座標も入れない。docs/34) |
| `RADPOL` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `RADTAPER` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `REAL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `RELW` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `RFSW` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py, maintrack.py) |
| `RING` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `SELFCOD` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `SORG` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `SPAC` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `STABLE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `SUMRES` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `SUS` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `TRPT` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `TWAKE` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `UNIFORM` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `UNIPOL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `UNSTABLE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `WAKEOPT` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `WSPAC` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |

## 5. 要素の型とキーワード (13 型)

SAD 側はいずれも Fortran (src/initb1.f の sethtb/kytbl、kyparam)。pysad の型は parser._SAD_ELEMENT_TYPES、キーワードは sadkw の表。

| 型 | pysad の型 | キーワード (SADHelp) | pysad に無いキーワード |
|---|---|---|---|
| `APERT` | あり | 0:  |  |
| `BEND` | あり | 18: AE1, AE2, ANGLE, DISFRIN, DISRAD, DROTATE, DX, DY, E1, E2, F1, FB1, FB2, FRINGE, K0, K1, L, ROTATE |  |
| `CAVI` | あり | 14: DISFRIN, DPHI, DVOLT, DX, DY, FREQ, HARM, L, PHI, ROTATE, V02, V1, V11, VOLT |  |
| `COORD` | あり | 0:  |  |
| `DECA` | あり | 7: DISFRIN, DISRAD, DX, DY, K4, L, ROTATE |  |
| `DODECA` | あり | 7: DISFRIN, DISRAD, DX, DY, K5, L, ROTATE |  |
| `DRIFT` | あり | 2: L, RADIUS |  |
| `MARK` | あり | 1: OFFSET |  |
| `MULT` | あり | 66: K0, SK0, K1, SK1, K2, SK2, K3, SK3, K4, SK4, K5, SK5, K6, SK6, K7, SK7, K8, SK8, K9, SK9, K10, SK10, K11, SK11, K12, SK12, K13, SK13, K14, SK14, K15, SK15, K16, SK16, K17, SK17, K18, SK18, K19, SK19, K20, SK20, K21, SK21, AE1, AE2, ANGLE, DISFRIN, DISRAD, DPHI, DVOLT, E1, E2, F1, F2, FB1, FB2, FREQ, FRINGE, HARM, L, misalignments, multipole_with_nonzero_ANGLE, PHI, RADIUS, VOLT | misalignments, multipole_with_nonzero_ANGLE |
| `OCT` | あり | 7: DISFRIN, DISRAD, DX, DY, K3, L, ROTATE |  |
| `QUAD` | あり | 10: DISFRIN, DISRAD, DX, DY, F1, F2, FRINGE, K1, L, ROTATE |  |
| `SEXT` | あり | 7: DISFRIN, DISRAD, DX, DY, K2, L, ROTATE |  |
| `SOL` | あり | 9: BOUND, BZ, DISFRIN, DPX, DPY, DX, DY, F1, GEO |  |

## 6. 手で直した行 (機械照合を作り直したら戻す)

- §1: MATRIX・RECOVER・REVERSE・SHOW・STATUS・VARY (ffsextra.py)、ELSEIF (sadscript.py)、MEASURE・RENUMBER・SEED・READ の注。
- §2: SynchroBetaEmittance・OptimizeOptics・PolynomialFit (.n を動かすもの)、ExternalMap・VariableRange・VariableWeight・WakeFunction (ユーザー定義のフック)、
  Graphics の行 (Plots.n などを動かす)、Fork・OpenShared・Shared・SharedSize・ListRandom・MemoryCheck・TracePrint (限定の代用)、
  Read・StandardForm・Fit・GaussianCoulomb・RadiationField・RadiationSpectrum・TrackParticles の注。
- §3・§4: MOMENTUM の注、POL・RADLIGHT の注。

## 7. 未確認・限界

- 「実装」は名前が登録されているという意味で、挙動が K64 と一致するかは要素ごとの文書による。Interp.builtins の登録名は 2026-11-16 に実行時に集めた (462 個。11/06 は 281 個)。
- 特別変数 (§3) は「ソース中に名前の文字列があるか」だけの照合。`CASE`・`ElementValues`・`TITLE`・`Track$Excitation` は参照が無い (未対応のまま)。
- `MessageName` は式のまま返す形だけ (関数形 `MessageName[f, "usage"]` の評価は見送り。docs/18)。`Skip` は gs に合わせた。
- `Polarization` (K64 が NaN: B55)、PHSROT (K64 で動作しない: B38) は pysad に無い。
- SADHelp は K64 版 `Documents/SADHelp.html` (1.4.4.2k64 Updated: 10/19/2026)。CERN の SADHelp は 1.4.0.0k64 2025-07-14 で同系統。KEK 公式の SAD ヘルプは別の SAD なので見ていない。

