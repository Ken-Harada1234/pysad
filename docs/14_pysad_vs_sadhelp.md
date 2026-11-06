# pysad と SADHelp (1.4.4.2k64 Updated: 10/19/2026) の対応表 — SAD 側の定義場所つき

機械照合 (`<work-perf>/pysad_vs_sadhelp.py`)。「SAD 側」は名前が登録・定義されているソース: Fortran (src/*.f)、C (src/*.c の REGn)、.n (Packages/*.n の定義)。`?` 付きは登録の形が見つからず文字列として出現するだけのもの。「pysad」は名前の登録の有無 (実装の深さは docs/03・04・12 を参照): コマンドは FFS._COMMANDS/_abbrev、関数は Interp.builtins、変数はソース中の参照、フラグは _FLAGS_DEFAULT/_SINO、要素は parser/sadkw。

## 1. FFS コマンド (56): 実装 32、未対応 12、無視 12

| コマンド | SAD 側 | pysad |
|---|---|---|
| `ABORT` | Fortran tffsa.f | 実装 |
| `APPEND(APP)` | Fortran phdrw.f,pstati.f,tffile.f | 未対応 |
| `ATTRIBUTE(ATTR)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `BEAMSIZE(BEAM)` | Fortran tffsa.f | 実装 |
| `BYE` | Fortran tffsa.f | 実装 |
| `CALCULATE(CAL)` | Fortran tffsa.f | 実装 |
| `CHROMATICITY(CHRO)` | Fortran tffsa.f | 実装 |
| `CLOSE(CLO)` | Fortran tffile.f | 未対応 |
| `COUPLE(COUP)` | Fortran tffsa.f | 実装 |
| `DISPLAY(DISP)` | Fortran tffsa.f,tracke.f | 実装 |
| `DRAW` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `DUMP` | Fortran terror.f | 実装 |
| `ELSE` | Fortran tfif.f | 実装 |
| `ELSEIF` | Fortran tfif.f | 未対応 |
| `EMITTANCE(EMIT)` | Fortran tffsa.f | 実装 |
| `END` | Fortran tffile.f | 実装 |
| `ENDIF` | Fortran tfif.f | 実装 |
| `EXECUTE(EXEC)` | Fortran tffile.f | 未対応 |
| `EXPAND` | Fortran tffsa.f,tftwiss.f | 実装 |
| `FIT` | Fortran tffsa.f | 実装 |
| `FITPOINTS(FITP)` | Fortran tffsa.f | 実装 |
| `FIX` | Fortran tffsa.f | 実装 |
| `FREE` | Fortran tffsa.f | 実装 |
| `GO` | Fortran tffsa.f | 実装 |
| `IF` | Fortran tfif.f | 実装 |
| `INPUT(IN)` | Fortran tffile.f | 未対応 |
| `MATRIX(MAT)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `MEASURE(MEA)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `ORG` | Fortran tffsa.f | 実装 |
| `OUTPUT(OUT)` | Fortran phdrw.f,pstati.f,tffile.f | 未対応 |
| `PRINT(PRI)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `QUIT` | Fortran tffsa.f | 実装 |
| `RADINT` | Fortran tffsa.f | 実装 |
| `READ` | Fortran tffile.f | 未対応 |
| `RECOVER(REC)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `REFERENCE(REF)` | Fortran tfdisp.f,tffsa.f | 実装 |
| `REJECT(REJ)` | Fortran tffsa.f | 実装 |
| `RENUMBER(RENUM)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `REPEAT(REP)` | Fortran tffsa.f | 実装 |
| `RESET` | Fortran tffsa.f | 実装 |
| `RESUME(RES)` | Fortran tffile.f | 未対応 |
| `REVERSE(REV)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `SAVE` | Fortran tffsa.f,tracke.f | 実装 |
| `SEED` | Fortran? initb1.f,tfinitn.f | 未対応 |
| `SHOW` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `SPLIT` | Fortran tffsa.f | 未対応 |
| `STATUS(STAT)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `STOP` | Fortran pstati.f,tffsa.f | 実装 |
| `SUSPEND(SUSP)` | Fortran tffile.f | 未対応 |
| `TERMINATE(TERM)` | Fortran tffile.f | 未対応 |
| `TYPE(T)` | Fortran tffsa.f | 実装 |
| `UNTIL` | Fortran tffsa.f | 実装 |
| `USE` | Fortran tffsa.f | 実装 |
| `VARIABLES(VAR)` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `VARY` | Fortran tffsa.f | 無視 (読み飛ばし、[pysad] FFS command … は未対応 と表示) |
| `VISIT` | Fortran tffsa.f | 実装 |

**未対応**: `APPEND(APP)`, `CLOSE(CLO)`, `ELSEIF`, `EXECUTE(EXEC)`, `INPUT(IN)`, `OUTPUT(OUT)`, `READ`, `RESUME(RES)`, `SEED`, `SPLIT`, `SUSPEND(SUSP)`, `TERMINATE(TERM)`  
**無視 (読み飛ばし)**: `ATTRIBUTE(ATTR)`, `DRAW`, `MATRIX(MAT)`, `MEASURE(MEA)`, `PRINT(PRI)`, `RECOVER(REC)`, `RENUMBER(RENUM)`, `REVERSE(REV)`, `SHOW`, `STATUS(STAT)`, `VARIABLES(VAR)`, `VARY`

## 2. SADScript / FFS の関数 (94、pysad 実装 46)

### Data-Manipulation (6、pysad 実装 1)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Fit` | .n Fit.n | 未対応 |
| `FitEmit` | .n Fit.n | 未対応 |
| `FitGaussian` | .n Fit.n | 未対応 |
| `NIntegrate` | .n D.n,OpticsPlot_230724b.n | 未対応 |
| `PolynomialFit` | .n Fit.n | 未対応 |
| `Spline` | .n D.n,OpticsPlot_230724b.n | 実装 |

**未対応**: `Fit`, `FitEmit`, `FitGaussian`, `NIntegrate`, `PolynomialFit`

### functional-operations (14、pysad 実装 10)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Apply` | Fortran tfinitn.f | 実装 |
| `Cases` | Fortran tfinitn.f | 実装 |
| `DeleteCases` | Fortran tfinitn.f | 実装 |
| `Difference` | Fortran tfinitn.f | 未対応 |
| `FixedPoint` | .n Misc.n | 実装 |
| `FixedPointList` | .n Misc.n | 未対応 |
| `Map` | Fortran tfinitn.f; .n correction.n | 実装 |
| `MapThread` | Fortran tfinitn.f; .n Canvas.n,bump.n | 実装 |
| `Position` | Fortran tfinitn.f; .n functions.n | 実装 |
| `Scan` | Fortran tfinitn.f; .n Canvas.n,Canvas3D.n,CheckEntryList.n,… | 実装 |
| `ScanThread` | Fortran tfinitn.f; .n correction.n | 実装 |
| `SelectCases` | Fortran tfinitn.f | 未対応 |
| `SwitchCases` | Fortran tfinitn.f | 未対応 |
| `Thread` | Fortran tfinitn.f | 実装 |

**未対応**: `Difference`, `FixedPointList`, `SelectCases`, `SwitchCases`

### FFS-dedicated-functions (28、pysad 実装 17)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `AccelerateParticles` | .n dapert.n | 実装 |
| `BeamMatrix` | .n Optics.n | 実装 |
| `DynamicApertureSurvey` | .n dapert.n | 実装 |
| `Element` | .n element.n | 実装 |
| `Emittance` | .n emit.n | 実装 |
| `ExternalMap` | .n init.n | 未対応 |
| `FFS` | Fortran tfinitn.f | 実装 |
| `FFS$SHOW` | Fortran tfinitn.f | 実装 |
| `FitValue` | Fortran? qdcell.f,tdfun.f,tffsa.f | 実装 (ユーザー定義のフック (tfgetfitval) を pysad が呼ぶ) |
| `FitWeight` | Fortran? tffscalc.f | 実装 (ユーザー定義のフック (twfit) を pysad が呼ぶ) |
| `GaussianCoulomb` | Fortran tfinitn.f | 未対応 |
| `GeoBase` | .n Optics.n | 実装 |
| `LINE` | Fortran tfinitn.f; .n OpticsPlot.n,OpticsPlot_230724a.n,correction.n,… | 実装 |
| `OptimizeOptics` | .n optimize.n | 未対応 |
| `OrbitGeo` | .n Optics.n | 実装 |
| `RadiationField` | Fortran tfinitn.f | 未対応 |
| `RadiationSpectrum` | Fortran tfinitn.f | 未対応 |
| `SetElement` | .n element.n | 実装 |
| `SurvivedParticles` | Fortran tfinitn.f | 実装 |
| `SymplecticJ` | .n emit.n | 実装 |
| `SynchroBetaEmittance` | .n emit.n | 未対応 |
| `TouschekLifetime` | .n emit.n | 実装 |
| `TrackParticles` | Fortran tfinitn.f | 実装 |
| `TransferMatrix` | .n Optics.n | 実装 |
| `Twiss` | Fortran tfinitn.f; .n correction.n,corutility.n,optimize.n | 実装 |
| `VariableRange` | Fortran? tffsmatch.f,tfvars.f | 文字列参照あり |
| `VariableWeight` | Fortran? tffsmatch.f | 文字列参照あり |
| `WakeFunction` | Fortran? tffswake.f | 未対応 |

**未対応**: `ExternalMap`, `GaussianCoulomb`, `OptimizeOptics`, `RadiationField`, `RadiationSpectrum`, `SynchroBetaEmittance`, `WakeFunction`

### Graphics (13、pysad 実装 0)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `BeamPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | 未対応 |
| `ColumnPlot` | .n Plots.n | 未対応 |
| `FitPlot` | .n Plots.n | 未対応 |
| `GeometryPlot` | .n GeometryPlot.n | 未対応 |
| `HistoPlot` | .n Plots.n | 未対応 |
| `LaTeX` | .n latex.n | 未対応 |
| `ListContourPlot` | .n Plots.n | 未対応 |
| `ListDensityPlot` | .n Plots.n | 未対応 |
| `ListPlot` | .n ListPlot.n | 未対応 |
| `OpticsPlot` | .n OpticsPlot.n,OpticsPlot_230724a.n | 未対応 |
| `ParametricPlot` | .n Plots.n | 未対応 |
| `Plot` | .n Plots.n | 未対応 |
| `PointSymbol` | .n? CSR.n,Canvas.n,ListPlot.n,… | 未対応 |

**未対応**: `BeamPlot`, `ColumnPlot`, `FitPlot`, `GeometryPlot`, `HistoPlot`, `LaTeX`, `ListContourPlot`, `ListDensityPlot`, `ListPlot`, `OpticsPlot`, `ParametricPlot`, `Plot`, `PointSymbol`

### Input/Output (15、pysad 実装 9)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `$FORM` | .n KBFrame.n,KBMainFrame.n,Plots.n,… | 文字列参照あり |
| `$Input` | .n init.n | 未対応 |
| `$Output` | .n init.n | 未対応 |
| `Close` | .n KBFrame.n,init.n | 実装 |
| `OpenAppend` | Fortran tfinitn.f | 実装 |
| `OpenRead` | Fortran tfinitn.f | 実装 |
| `OpenWrite` | Fortran tfinitn.f | 実装 |
| `PageWidth` | .n latex.n | 文字列参照あり |
| `Print` | Fortran tfinitn.f; .n Class.n | 実装 |
| `Read` | Fortran tfinitn.f | 実装 |
| `ReadString` | Fortran tfinitn.f | 実装 |
| `StandardForm` | Fortran tfinitn.f; .n KBMainFrame.n | 未対応 |
| `StringToStream` | Fortran tfinitn.f | 未対応 |
| `Write` | Fortran tfinitn.f | 実装 |
| `WriteString` | Fortran tfinitn.f | 実装 |

**未対応**: `$Input`, `$Output`, `StandardForm`, `StringToStream`

### Multiprocessing (4、pysad 実装 0)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Fork` | C tfProcess_.c | 未対応 |
| `OpenShared` | .n shared.n | 未対応 |
| `Shared` | .n CSR.n,shared.n | 未対応 |
| `SharedSize` | Fortran tfinitn.f | 未対応 |

**未対応**: `Fork`, `OpenShared`, `Shared`, `SharedSize`

### Object-oriented-programing (1、pysad 実装 1)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Class` | .n Class.n | 実装 |

### Random-number-functions (5、pysad 実装 4)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `Random` | C tfRandom_.c | 実装 |
| `GaussRandom` | C tfRandom_.c | 実装 |
| `ParabolaRandom` | C tfRandom_.c | 実装 |
| `SeedRandom` | C tfRandom_.c | 実装 |
| `ListRandom` | C tfRandom_.c | 未対応 |

**未対応**: `ListRandom`

### System-interface (2、pysad 実装 1)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `System` | C tfProcess_.c | 実装 |
| `TemporaryName` | Fortran tfinitn.f | 未対応 |

**未対応**: `TemporaryName`

### Utilities (5、pysad 実装 3)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `DateString` | .n init.n | 実装 |
| `MemoryCheck` | Fortran tfinitn.f | 未対応 |
| `TimeUsed` | Fortran tfinitn.f | 実装 |
| `Timing` | .n functions.n | 実装 |
| `TracePrint` | Fortran tfinitn.f | 未対応 |

**未対応**: `MemoryCheck`, `TracePrint`

### (直下) (1、pysad 実装 0)

| 関数 | SAD 側 | pysad |
|---|---|---|
| `DownhillSimplex` | .n optimize.n | 未対応 |

**未対応**: `DownhillSimplex`

## 3. 特別変数 (48、pysad 参照あり 42)

| 変数 | SAD 側 | pysad |
|---|---|---|
| `$Line` | Fortran tfinitn.f; .n SADTerminal.n,init.n | 参照あり (sadscript.py) |
| `CASE` | Fortran tffsa.f | 未対応 (参照なし) |
| `CHARGE` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) |
| `CONVERGENCE` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py) |
| `DAPWIDTH` | Fortran initb1.f,tfdapert.f,trackd.f | 参照あり (ffscmd.py) |
| `DP` | Fortran tffs.f,tfsetparam.f,tftwiss.f | 参照あり (emit_extra.py, ffscmd.py, match6.py…) |
| `DP0` | 不明 | 参照あり (ffscmd.py, parser.py) |
| `DTSYNCH` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFRFFREQ` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFVC` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `EFFVCRATIO` | Fortran initb1.f,track.f,tturne.f | 参照あり (ffscmd.py) |
| `ElementValues` | Fortran tffsa.f | 未対応 (参照なし) |
| `EMITX` | Fortran initb1.f,temit.f,tffs.f,…; .n emit.n | 参照あり (ffscmd.py, maintrack.py, match6.py…) |
| `EMITXE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `EMITY` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, match6.py…) |
| `EMITYE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `EMITZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, parser.py…) |
| `EMITZE` | Fortran initb1.f,temit.f | 参照あり (ffscmd.py) |
| `ExponentOfResidual` | Fortran tffsa.f,tffsmatch.f | 参照あり (ffscmd.py) |
| `FFS$NumericalDerivative` | Fortran tffsmatch.f | 参照あり (ffscmd.py, sadscript.py) |
| `FitFunction` | Fortran qdcell.f,tdfun.f,tffsa.f | 参照あり (ffscmd.py) |
| `FSHIFT` | Fortran ActTra.f,initb1.f,temitf.f,… | 参照あり (ffscmd.py, parser.py) |
| `GCUT` | Fortran initb1.f; C tfRandom_.c,tgauss_.c | 参照あり (ffscmd.py) |
| `InitialOrbits` | Fortran tffsa.f | 未対応 (参照なし) |
| `LOSSAMPL` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) |
| `LOSSDZ` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, maintrack.py, parser.py) |
| `MatchingAmplitude` | Fortran tffamsetup.f | 参照あり (ffscmd.py) |
| `MatchingResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) |
| `MASS` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, optics_fn.py, parser.py…) |
| `MINCOUP` | Fortran initb1.f,tfbeam.f,tffs.f,… | 参照あり (ffscmd.py, parser.py) |
| `MOMENTUM` | Fortran ActTra.f,bb.f,bbstrhl.f,…; .n dapert.n | 参照あり (ffs.py, ffscmd.py, optics_fn.py) |
| `NBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) |
| `NetResidual` | Fortran tffscalc.f | 参照あり (ffscmd.py) |
| `NP` | Fortran initb1.f,tfdapert.f,tffsa.f,…; .n dapert.n | 参照あり (ffscmd.py, maintrack.py) |
| `NPARA` | Fortran initb1.f,tfefun1.f,tfsetparam.f,…; .n init.kekb.n,init.n,init.skekb.n | 参照あり (ffscmd.py, parser.py) |
| `OffMomentumWeight` | Fortran tffsmatch.f | 参照あり (ffscmd.py) |
| `OMEGA0` | Fortran ActTra.f,initb1.f,tfgeo.f,… | 参照あり (ffscmd.py) |
| `OpticsEpilog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) |
| `OpticsProlog` | Fortran tffsa.f,tffscalc.f | 参照あり (ffscmd.py) |
| `PBUNCH` | Fortran initb1.f,tfsetparam.f,track.f | 参照あり (ffscmd.py, parser.py) |
| `PhotonList` | Fortran photons.f | 未対応 (参照なし) |
| `PHICAV` | Fortran initb1.f,tturne.f | 参照あり (ffscmd.py) |
| `SIGE` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, maintrack.py, parser.py…) |
| `SIGZ` | Fortran initb1.f,temit.f,tffs.f,… | 参照あり (ffscmd.py, parser.py) |
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
| `CALC6D` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py, optics_fn.py) |
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
| `EMIOUT` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
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
| `LWAKE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `MOVESEED` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `PHOTONS` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `POL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `PRSVCASE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `PSPAC` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `QUIET` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `RAD` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `RADCOD` | Fortran tffs.f (fname/sino) | 実装 (参照 ffscmd.py) |
| `RADLIGHT` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `RADPOL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
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
| `TWAKE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `UNIFORM` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `UNIPOL` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `UNSTABLE` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
| `WAKEOPT` | Fortran tffs.f (fname/sino) | 実装 (値を保持するだけ) |
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

