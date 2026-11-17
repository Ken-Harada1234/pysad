# 01. pysad 全体像 — SAD (K64 版) の Python 直訳

最終更新: 2026-11-17 (現状のスナップショット)

> 場所: コード・文書は `<repo>\` (pysad 本体は `pysad\`、文書は `docs\`、旧文書は `<work-olddocs>\`)。
> 対象の SAD は KEK 本流ではなく **K64 版** (https://github.com/KatsOide/SAD)。コマンド・要素キーワードの意味は CERN の SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html、1.4.0.0k64 2025-07-14) か、ローカルの `<SAD>\Documents\SADHelp.html` (1.4.4.2k64 2026-10-19、使っている WSL ビルドと同じ版) を先に見る。**KEK 公式サイト (acc-physics.kek.jp) は別系統の SAD なので参照しない** (ユーザー指示)。数値の最終根拠はローカルのソース `<SAD>\src` と `Packages\`。
> 再開時はまず `docs/history/00_RESUME.md` (状態・残課題) を読む。

## 1. 目的と方針

- 目的: SAD の光学計算・エミッタンス・マッチング・粒子追跡・DA を Python で再現し、SAD スクリプトを `run_sad.py` でそのまま動かす。
- 方針: **再導出せず直訳**。SAD の Fortran の式・分割数・反復順序・乱数の引き方を写し、SAD の数値 (印字精度 1e-10〜1e-11) を再現することを一致の基準にする。SAD の癖 (バグに見えるものも) は再現する (→ `11_sad_quirks.md`)。
- 進め方: SAD スクリプトを下書き → ユーザーが確認・改修 → WSL で実行して参照値 (`bench_level*.result`) → pysad で直訳して一致させる。長い計算は進捗表示 + 発散で停止。返答・コメントは日本語。
- 対象外: EPICS、`System[]` の実行、Python ブリッジ (SAD 側の PySad.c は非コンパイル)。Tcl/Tk の GUI と描画 (`ListPlot`・`Plot`・`OpticsPlot` など) は、SAD 自身の `Packages/*.n` を pysad のインタプリタで解釈して tkinter に表示する方式で対応した ([22_gui_tclbridge.md](22_gui_tclbridge.md)、[37_plot_canvas.md](37_plot_canvas.md))。ただし見た目 (色・大きさ・フォント)・操作の動作は SAD と比べておらず、`Canvas3D*` は未対応。11/16 夜の確認: `ListPlot`・`Plot` のオプション 10 例で戻りの Graphics 式は gs と文字列で完全一致 (見た目ではなく式の比較。`Options[g]` は gs が全オプションを返し pysad は `{}`)。**差あり**: `w=Window[]`・`e=Entry[w]` の直後、gs では値が `Widget[…]` になり `e[Insert]`・`e[Get]`・`e[Delete]` が働くが、pysad では w・e が裸のシンボルのままでメソッド呼び出しは何もしない (ウィジェット自体は作られ、親は toplevel でなく Root になる)。原因は `Tkinter.n` の `Widget[a$,w[master$]]` が値にならないところまで絞ったが未修正 (docs/22・39)。
- 要素: DRIFT/BEND/QUAD/SEXT/OCT/DECA/DODECA/MULT (ANGLE≠0 を含む)/SOL/CAVI/TCAVI/UND/COORD/MAP/BEAMBEAM/MARK/MONI/APERT/INS を扱う。PHSROT は K64 で動作しない (行列が 0、B38) ため追跡は NotImplementedError。SPCH・PROF は parser が読むだけで物理の翻訳はない (11/16 に K64 のソースで確認: tturn.f:521 の `case (icSPCH)` が PSPAC フラグのとき `tpspac` (粒子のメッシュによる空間電荷) を呼ぶが、pysad は未翻訳で PSPAC・WSPAC・SPAC フラグも効果なし。PROF は K64 のソースに現れない)。要素別の状況は §6 と各要素の文書 (23〜30)。

## 2. SAD の構造

```
MAIN 層 (toplvl.f)   要素定義 / LINE / 数値代入 / FFS, TRACK, EMIT の起動。SADScript は使えない
   └ FFS USE=…;  → track.f (tfinitn で SADScript 初期化、RFSW=ON, CALC6D=OFF 強制) → tffs → tffsa (コマンドループ)
        FFS 層 (tffsa.f)  1 語の解釈順: tfprint (式に見える語は SADScript として評価) → tflag (フラグ) → tgetfv (FIT 条件) → terror (誤差) → tffile → 要素名 (tfgetlineps→tfgetv) → コマンド群 → tfgetv (要素値、失敗なら ?Undefined command) (詳細は 03_ffs_commands §1)
        SADScript 層      tfetok (字句) / tfeval (構文、優先順位は tfstk.f iprior) / tfefun* (組み込み関数)、Packages/*.n
物理計算   CAL/GO: tffsmatch → tffscalc → qcell1 (qcod → qmdiag → qtwiss)
           Emittance[]: emit.n → temit → tturne (tquade/tbende/tthine/tcave + tradke) → tecalc
           TrackParticles[]: tftrack → tturn → tturn1 (tdrift/tquad/tthin/tbend/tcav/tmulti/tsol + tradk)
           DynamicApertureSurvey: dapert.n → tfdapert → trackd → tturn
```

### 主要ソースと役割

| ファイル | 役割 |
|---|---|
| tffsa.f | FFS コマンドループ (USE/CELL/INS/CAL/GO/FIT/FREE/SAVE…)、`FFS[]` (tfffs) |
| tffscalc.f | 1 回分の光学計算 + 残差 (運動量点ループ、重み twfit、nstab) |
| tffsmatch.f | GO の反復本体 (Newton/DESCEND 切替、fuzz、tfsolv/tgrad、fmincube、tweigh) |
| qdcell.f | tdfun (残差 df)、tffsqu/qdcell (解析微分) |
| qcell.f | qcell1: CELL の周期解 (qcod → qmdiag、結合込み分散) と INS、qtwiss 呼び出し |
| qtwiss.f | Twiss 伝搬 (qtwiss1/qmat2twiss)、**qcod** (4D 閉軌道、Newton + SVD)、qquad/qthin (ファイル qquad.f は非コンパイル) |
| temit.f | temit/tecalc: 6D 平衡ビーム行列 (固有値分解 + SVD)、synchm (RF あり) |
| tturne.f | emit 用の放射込み 6D 一周 (tquade/tbende/tthine/tcave/tmulte、trf0 更新) |
| tturn.f | 追跡の一周 (tturn1 の要素分岐、tapert 損失判定) |
| tbend.f / tbendi.f / tbende.f | 追跡 BEND (tbendcore/tbendr、K1≠0 は tbendi) / emit 用 BEND 行列 (qbend も内部) |
| tquad.f / tquade.f | 追跡の四極 tquad・薄い多極 tthin・ttfrin / emit 用四極行列 |
| tsolqu.f / tsolque.f | 四極本体の分割写像 (追跡 tsolqu/tsolqur/tsolqum、emit tsolque)、ソレノイド bz 分岐 |
| tthine.f | emit 用の薄い多極 (L/6-2L/3-L/6) |
| tcav.f / tcave.f | 空洞 (追跡 / emit)、tcavfrin フリンジ |
| tmulti.f / tmulte.f、tmulta.f | MULT (追跡 / emit、ANGLE=0 と ANGLE≠0) |
| photons.f | 放射: tradk→tradkn (平均) / tradkfn (FLUC、乱数)、tradke (emit 用、radcod) |
| trackd.f / tfdapert.f | DynamicAperture の粒子配置・境界探索・Score |
| tftwiss.f | Twiss[] / LINE[] / Element[] 関数 (値の取り出し。Twiss の計算本体ではない) |
| tgetfv.f / terror.f / tfgetv.f / tflag.f | FIT 条件 / 誤差コマンド DELX… / 要素値 "QF K1 v" / フラグ語 (tffs.f fname/sino) |
| autofg.f / tftype.f | 数値書式 (`$FORM`) / TYPE の書式 |
| tfstk.f | 中核モジュール (型コード icBEND=2 等、物理定数 macphys、演算子 opcode/iprior) |
| random_plugin_sad.c | 既定乱数 "SAD" (合同法、Box–Muller) |
| Packages/init.n, emit.n, Fit.n, dapert.n, Optics.n | 既定値 (NPARA=4 等)、Emittance[] の規則名、PolynomialFit、DynamicApertureSurvey、NormalCoordinates |

フラグ (RAD, RFSW, FLUC/DAMPONLY, CELL/INS, RADCOD, INTRES…) は tffs.f の `type flagset` と tflag.f。初期値は initb1.f の defflg と track.f のコピー。

## 3. pysad のモジュール地図

| モジュール | 内容 | 直訳元 |
|---|---|---|
| parser.py / sadkw.py / sadelem.py | MAIN 部 (要素定義、LINE、`-LINE`, `N*`, 再定義、LINE 名は書かれたとおり) の読み込み → `Lattice`。K1/K2/K3 は per-length に変換 (元の値は `_raw`)。Element/SetElement/BeamLine/WriteBeamLine | toplvl/doelem/doline/expln、tfelement |
| elements.py | 要素クラス (Drift, Bend, Quadrupole, Sextupole, Octupole, Multipole, Cavity, Solenoid, Marker…) と線形行列 | — |
| lattice.py / elkern.py / perfkern*.py / turnkern.py | CALC: 4D 閉軌道 `_qcod_periodic` (1 周を numba 1 回に)、要素ごとの COD 写像 (tquade/tbende/tthine/tqlfre/tblfre/tchge/tbrote/MULT/tsteee)、`twiss()`/`twiss_table()` (qcell1 結合版、不安定時の 1/cosμ 続行)、`ins_twiss()` (INS の初期軌道オフセット、TRPT の加速)、MARK OFFSET、`set_element`/`set_error`、要素カーネルの numba 化 | qtwiss.f, qcell.f, tquade.f, tbende.f, tthine.f, tchge.f |
| calc6d.py | CALC6D (qcell6d、CELL と INS)、NormalCoordinate6 (etwiss2ri) | qcell6d, temit |
| geo.py | 幾何 (GX/GY/GZ/CHI、GEO、ORG、DISP G/OG、fit の GX..) | tfgeo.f, geolib |
| radint.py | RADINT (放射積分 I1〜I5、結合・色収差) | INTGRL.f |
| wake.py / wakeopt.py | ウェイク場 (txwake、spsort。RF のある MULT は track6n_ext.mult_full)、光学の wake (WAKEOPT) と InitialOrbits | txwake.f, spkick.f, tffswake.f, tmulti.f tmultiacc |
| trptacc.py | TRPT の加速 (CAVI と VOLT 付き MULT の gammab 分布、tcave/tmulte の acc 部分) | tfgeo1 tfacc, tcave.f, tmulte.f |
| emit6d.py / teigen.py / emit_extra.py | Emittance[] / EMIT: 6D 一周 (T + D、量子励起 B)、tradke、tcave/trf0、tecalc (synchm、固有値 teigen/tnorm の直訳)、close_z、RADCOD、報告の全行、EMITXE.. の下限、Emittance[] のオプション (OneTurnInformation/Matrix) | temit.f, tturne.f, photons.f tradke, tcave.f, tmulte.f, emit.n |
| solmap.py | SOL 領域の COD 写像と行列 (CALC/INS の Twiss。追跡のカーネルを数値微分。docs/29) | tsole.f tsole1 |
| solemit.py | SOL 領域の放射込みの 1 周 (Emittance[]・EMIT: bz 付き tradke、tsolque の刻み。docs/29) | tsole.f, tsolque.f, tqente.f, photons.f tradke |
| multa.py / multae.py | MULT ANGLE≠0 (追跡 / COD 写像・行列。docs/28) | tmulta.f (tmulta, tmultae) |
| touschek.py | INTRA (IBS 付き emit)、TouschekTable / TouschekLifetime | temit.f tintrb, emit.n |
| match6.py / match6d.py / tsvd.py | FIT/FREE/FIX/GO: tffsmatch (Newton/DESCEND、fuzz、fmincube、tweigh、tfsolv SVD)、tdfun、twfit、FITP>1 と fam 点 (MatchingAmplitude/NFAMP)、FitValue/FitWeight/FitFunction/VariableRange/VariableWeight、OpticsProlog/Epilog。既定は解析微分 (match6d、qdcell/tffsqu の直訳)、`FFS$NumericalDerivative` で数値微分 | tffsmatch.f, tffscalc.f, qdcell.f, tgetfv.f, tsolvg.f |
| track6.py / track6n.py / track6n_ext.py | TrackParticles: 純 Python 版 (要素外側・粒子内側、tdrift/tquad/tthin/tbend/tbendi/tsteer/tcav/tmulti/tsol、tradk/tradkfn) と numba 版 (粒子配列 S(n,10)、同じ演算順序、乱数も int64)。粒子損失 (tapert、alive フラグ)、MULT/SOL 領域のカーネル。POL かつ RAD ON のときのスピン追跡 (`sprot`・`bsi`。純 Python 版のみで損失の処理なし、docs/36・43) | tftrack.f, tturn.f, tquad.f, tsolqu.f, tbend.f, tcav.f, tmulti.f, photons.f |
| dapert.py | DynamicApertureSurvey (スロット順・tapert・Score・62 進地図) | dapert.n, tfdapert.f, trackd.f |
| maintrack.py | MAIN 層の TRACK コマンド (`TRACK USE= NPART= TURNS= …`) | trackm (K64 では動かない機能) |
| optics_fn.py | TransferMatrix / BeamMatrix / NormalCoordinates / OrbitGeo / GeoBase / SymplecticJ / SurvivedParticles / AccelerateParticles | Optics.n, emit.n, dapert.n (Python で書き直し) |
| sgeo.py | SOL 領域の幾何 (tsgeo。領域内の BEND・MULT・CAVI と tdgeo の領域幾何は未対応。docs/29) | tsgeo.f |
| gausscoul.py / bbbrem.py | GaussianCoulomb(U・Fitted) (bbkick1・twspu・twspfu)、BBBrem1 (docs/32) | tspac.f, tintrb.f, bint.f, tbbbrem.f |
| beambeam.py | BEAMBEAM 要素 (強弱ビーム、ビームストラールング、ルミノシティ。追跡のみで Twiss には効かない。docs/24) | bbstrhl.f, bbtool.f, matrix.f |
| coord.py / tcav.py / und.py | COORD (追跡・6 次元写像・幾何、docs/25)、TCAVI (追跡 RAD FLUC まで、Emittance[] は L=0 のみ gs と一致、docs/26)、UND (追跡のみ。K64 は磁場を設定できず磁場なしのドリフト、docs/27) | tcoord.f, tcoorde.f, qcoord.f, ttcav.f, ttcave.f, undulator.f |
| specfn.py | 特殊関数・数学関数 (Bessel・超幾何・ガンマ族・ゼータ・多重対数・D.n の特殊関数など。K64 の InverseErf の誤差も再現。docs/33) | gamma.f, hg.f, Packages/D.n |
| fitroot.py / miscfn.py / fileops.py | Fit・FindRoot・NIntegrate・FitGaussian、DownhillSimplex、SingularValues、システム・ファイル関数、Parallelize (docs/33) | tffindroot.f, Fit.n, D.n, optimize.n, Matrix.n, FileSystem.n, shared.n |
| corr.py / pkgload.py | ResponseMatrix・LinearSolveConditioned・Master と、Packages/*.n (CorrectOrbit・MakeBump・Monitor・Steer・OptimizeOptics・D・rational など) の自動読み込み (docs/34) | pgrmat.f, pgsolvcond.f, pmicad.f, pgmast.f, correction.n ほか |
| ffsextra.py | FFS の VARIABLES・ATTRIBUTE・MATRIX・SHOW・VARY・RECOVER・REVERSE・STATUS・PRINT・EXECUTE・OUTPUT/APPEND/CLOSE・INPUT/READ (docs/35) | tfvars.f, tfattr.f, tfmat.f, tshow.f, tffsa.f |
| spin.py | POL (スピン): tradke のスピン部・tchge の srot 回転・spnorm・srequpol/spdepol・Emittance[]/EMIT の PolarizationTime/SpinTune/PolarizationVector (MULT・空洞・SOL まで。`Polarization` は K64 が NaN: B55。docs/36) | photons.f (tradke, tspin), tchge.f, temit.f |
| synchrobeta.py | SynchroBetaEmittance1 (シンクロ・ベータ結合のエミッタンス。運動量のずれごとの閉軌道と tturne を取り方位角モード展開を解く。gs と emitx 相対 4e-7〜1e-5) | tfemit.f tfemits, temits.f ほか |
| measure.py | tlum / tluma (ルミノシティの推定)。MEASURE が gs で動かない (B42) ので解析値との比較のみ | tlum.f |
| csr.py | CSRTrack・CSRHaissin・CSRSetupOY・CSROYMatrix・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRConvert (= tfcsrvec)・CSRInit と実FFT・spline1 (docs/36。`Packages/CSR.n` の MakeZL などは gs と一致、ScanTrack・OYStability は未照合) | csrtrack.f, csroy.f, csrtest.f, tcftr.f, spline.f |
| tkcanvas.py | Tk のキャンバスの組み込み関数 (TkCanvasCreateItem、CanvasClipLine、CanvasSymbolDirect …)。ListPlot・Plot の描画 (docs/37) | tfTkInter_.c, tfcanvasclip.f |
| vecfn.py | Vector・ReleaseVector (ベクトル評価。Plots.n の PlotInterpolate) | — |
| radfield.py | RadiationField・RadiationSpectrum。RADLIGHT の軌道記録 (`TrackParticles` が `{結果, 軌道{t,x,y,z}}` を返す) は `ffscmd._RadLight` (gs と差 0、docs/45) | tlinit.f, tlstore, tlresult |
| sadrandom.py | SAD の乱数 (LCG、Box–Muller、drndsr/tdusrn、ParabolaRandom) | random_plugin_sad.c, tfRandom_.c, drndsr.f |
| sadscript.py / sadlist.py | SADScript インタプリタ (字句・Pratt 構文・評価、Module/Do/Table/純関数/パターン/REPEAT/UNTIL/IF、単位付き数値、`?FLAG`、fmt_real=autofg) とリスト・文字列関数 | tfetok.f, tfeval.f, tfstk.f iprior, autofg.f, functions.n, Misc.n |
| wake.py | 追跡へのウェイク場の適用 (WakeFunction、LWAKE/TWAKE。TrackParticles の TRPT のみ) ([23_wake.md](23_wake.md)) | txwake.f, tffssetupwake (tffswake.f), tcav.f / tturn.f の呼び出し箇所 |
| tclbridge.py / sadclass.py | Tcl/Tk の接点 (tkinter) と SAD の Class。SAD 自身の Tkinter.n / KBFrame.n を SADScript として解釈して GUI を動かす ([22_gui_tclbridge.md](22_gui_tclbridge.md)) | tfTkInter_.c, tfTclArg.f, tfXlib_.c, Tkinter.n, KBFrame.n, KBMainFrame.n, TabFrame.n, SADHelp の Class |
| ffscmd.py | FFS コマンド層 (tffsa の解釈順)、USE/CELL/INS/CALC/CALC6D/GO/EMIT/TYPE/DISP/SAVE/RESET/RADINT/CHROMATICITY/誤差コマンド (DELX…)、Twiss[]/LINE[]/FFS[]/Element[]/Emittance[]/TrackParticles[]/CalculateOptics/SeedRandom[]、EMIT の結果変数、`FFS$SHOW` | tffsa.f, tftwiss.f, tflag.f, terror.f, tftype.f, tfdisp.f |
| diag_fit.py / diag_trace.py | 照合用の診断 (fit の反復表、追跡の途中値) | — |
| run_sad.py (<repo> 直下) | `py -3.12 run_sad.py lattice script …` = `gs` に複数ファイルを連結して渡すのと同じ | — |
| **非推奨** ffs.py / matching.py / tracking.py / optics.py | 独自実装 (SAD 未照合。tracking.py は放射係数に疑い)。新規コードは使わない。`__init__.py` の export は旧 API のまま残っている | — |

テスト: `pysad/tests/` (pytest 578 件を収集。11/17 の実行で **575 通過・3 skip**。回帰 `run_regression.py` は 50 ベンチをすべて実行して FAIL 0 (205 s、QUICK=1 は遅い 4 件 7・12t・13a・16a を SKIP して 46 件で、4 件も単独で実行して PASS)。`--ignore` する 3 本 `test_bend_feeddown.py`・`test_coord.py`・`test_sext_radiation.py` は sys.exit で終わる旧式のスクリプト。GUI 系 (`test_plot`・`test_tclbridge`) は全体実行で数回に 1 回失敗することがあり (単独なら通る)、環境変数 `PYSAD_SAD_PACKAGES=D:/work/SAD-master/SAD-master/Packages` が無いと skip。`PYSAD_SLOW=1` で遅い試験 (CSR の MakeZL、約 7 分) も走る)。試験が走らせる SAD スクリプトと gs・pysad の画面表示は `pysad/tests/scripts/` (`NAME.sad`・`NAME.gs.txt`・`NAME.pysad.txt`)、試験関数との対応は `pysad/tests/OUTPUTS.md`。
ベンチ: `bench_level1〜31*.sad/.result` と `*_pysad.py/.result`、`run_regression.py` (50 ベンチ。QUICK=1 は遅い 4 件 (7・12t・13a・16a) を SKIP して 46 件を実行、`QUICK=1` で短縮。11/16 に FAIL 0)。

## 4. docs/ の案内

| 文書 | 内容 |
|---|---|
| history/00_RESUME.md | 作業メモ (時系列、残課題)。再開時に最初に読む。履歴はここと work_*.md にだけ書く |
| 01_overview.md | 本書 |
| 02〜09 | SAD の構造と各計算の直訳の解説: MAIN/LINE/要素 (02)、FFS コマンド (03)、SADScript (04)、CALC/COD/Twiss (05)、マッチング (06)、Emittance (07)、追跡 (08)、DynamicApertureSurvey (09) |
| 10_benchmarks.md | ベンチ (bench_level*) の一覧と現在の一致度 |
| 11_sad_quirks.md | SAD の癖・仕様・ビルド依存の一覧 (現象 / 根拠 / pysad での扱い) |
| 12_pysad_api.md | pysad の Python API |
| 16_fortran_status.md | SAD の Fortran ファイル別の翻訳状況  |
| 17_implementation_summary.md | 実装まとめ (モジュール行数、SAD 準拠状況、課題) |
| 18_physics_guide.md / 19_source_analysis.md / 20_translation_candidates.md | 物理の式 / SAD ソースの解析 / SAD 機能の翻訳候補 |
| 21_sad_bug_report.md | SAD K64 のバグ報告 (B01〜) |
| 13_sad_help_index.md / 14_pysad_vs_sadhelp.md / 15_sadhelp_purpose_policy.md / sadhelp_summary.html | SADHelp の一覧、pysad との対応、目的別分類と翻訳方針 (<work-perf>/sadhelp_*.py で生成) |
| history/work_*.md | 項目ごとの作業メモ (経緯・照合の記録)。`history/work_unimpl.md` は未翻訳・未解決の差の最新の記録 |
| html/ | `make_html.py` で docs/*.md から生成した HTML |

### 4.1 要素・機能別の文書 (docs/22〜37、11/13〜11/16 に追加)

| 文書 | 1 行の説明 |
|---|---|
| [22_gui_tclbridge.md](22_gui_tclbridge.md) | Tcl/Tk の GUI (Tkinter.n・KBFrame.n) を、SAD 自身の `.n` を解釈して tkinter に出す仕組み。Class も含む |
| [23_wake.md](23_wake.md) | ウェイク場 (WakeFunction、LWAKE/TWAKE) の追跡と光学 (WAKEOPT・InitialOrbits) |
| [24_beambeam.md](24_beambeam.md) | BEAMBEAM 要素 (強弱ビームのビーム・ビーム相互作用、ビームストラールング) |
| [25_coord.md](25_coord.md) | COORD 要素 (座標変換: ずれと回転) の追跡・Twiss・幾何・Emittance[] |
| [26_tcav.md](26_tcav.md) | TCAVI (横偏向空洞・クラブ空洞)。追跡は RAD FLUC まで、Emittance[] は L=0 のみ |
| [27_und.md](27_und.md) | UND (アンジュレータ)。追跡のみ。K64 は磁場を設定できずドリフト (B37) |
| [28_mult_angle.md](28_mult_angle.md) | MULT の ANGLE≠0 (曲げつき多極) の追跡・Twiss・Emittance[] |
| [29_sol_optics.md](29_sol_optics.md) | SOL 領域の Twiss (CALC/INS)・Emittance[]/EMIT・幾何 (tsgeo)。領域内の MULT・BEND・CAVI は未対応 |
| [30_map.md](30_map.md) | MAP (ExternalMap) の追跡・Emittance[]・Twiss。GEO・CompiledMap は未対応 |
| [31_delcor.md](31_delcor.md) | DELCOR (相関のあるずれ)。tcorr・tvcorr と乱数列まで gs と一致。CAVI は未対応 |
| [32_gaussiancoulomb_bbbrem.md](32_gaussiancoulomb_bbbrem.md) | GaussianCoulomb(U・Fitted) と BBBrem1 |
| [33_specfn_fit_misc.md](33_specfn_fit_misc.md) | 特殊関数・数学関数、Fit・FindRoot・NIntegrate、DownhillSimplex、ファイル・文字列の小さな関数 |
| [34_sad_packages.md](34_sad_packages.md) | `Packages/*.n` (CorrectOrbit・MakeBump・OptimizeOptics・D など) を SADScript のまま使う仕組み |
| [35_ffs_small_commands.md](35_ffs_small_commands.md) | FFS の小さなコマンド (VARIABLES・ATTRIBUTE・MATRIX・SHOW・VARY・OUTPUT・INPUT・STATUS など) |
| [36_spin_csr.md](36_spin_csr.md) | POL (スピン、偏極時間・スピンチューン・偏極ベクトル) と CSR のウェイク追跡の関数 |
| [37_plot_canvas.md](37_plot_canvas.md) | 描画 (ListPlot・Plot・OpticsPlot と Tk のキャンバスの組み込み関数) |

### 4.2 解説文書 (docs/38〜45、初めて読む人向けに「できること・使い方・物理と K64 の式・実装・K64 の癖・gs との一致・限界」の順で書く)

| 文書 | 1 行の説明 |
|---|---|
| [38_wake_guide.md](38_wake_guide.md) | ウェイク場 (WAKE): LWAKE/TWAKE の追跡と光学 (WAKEOPT) の解説 |
| [39_beambeam_guide.md](39_beambeam_guide.md) | ビーム・ビーム (BEAMBEAM、GaussianCoulomb、ルミノシティ) の解説 |
| [40_csr_guide.md](40_csr_guide.md) | コヒーレント放射 (CSR): CSRTrack・CSRHaissin・CSROY 系・CSRMatrix の解説 |
| [41_tcltk_guide.md](41_tcltk_guide.md) | Tcl/Tk (GUI): Tkinter.n・KBFrame.n を pysad で動かす方法の解説 |
| [42_graphics_guide.md](42_graphics_guide.md) | グラフィック描画: ListPlot・Plot・OpticsPlot・Canvas の解説 |
| [43_pol_guide.md](43_pol_guide.md) | POL (スピン偏極): Emittance[]・EMIT・TrackParticles のスピンの解説 |
| [44_sol_guide.md](44_sol_guide.md) | SOL (ソレノイド): Twiss・Emittance[]・追跡・幾何の解説 |
| [45_radlight_guide.md](45_radlight_guide.md) | RADLIGHT (放射光の軌道記録) と RadiationField・RadiationSpectrum の解説 |

旧文書 (`<work-olddocs>/`) はソース照合済みの初期の解説で、pysad の状況は古い。現状は本書と上の各文書を見る。

## 5. 使い方

```
# SAD スクリプトをそのまま実行 (MAIN 部 + FFS 部)。numba 入りの Python 3.12 を使う
cd "<repo>"
py -3.12 run_sad.py <試験用ラティス>.lat bench_level13b.sad

# Python API (詳細は 12_pysad_api.md)
from pysad.parser import load_sad_file      # → SADScript (env に Lattice と LINE)
lat.twiss_table(delta=0.0)                   # CELL CALC (qcell1)。lat は Lattice
from pysad.emit6d import emittance6d; emittance6d(lat, 3e9, rfsw=True, radcod=False)
from pysad.match6 import Matcher; m = Matcher(lat); m.fit('NX', 47.6); m.free('Q*'); m.go()
from pysad.track6n import track_particles_fast   # TrackParticles の本体
from pysad.dapert import dynamic_aperture        # DynamicApertureSurvey の本体

# 回帰 (SAD 参照値と比較)
QUICK=1 py -3.12 run_regression.py
py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py

# SAD 本体 (参照値の取得): WSL Ubuntu
MSYS_NO_PATHCONV=1 wsl -d Ubuntu -- <SAD>/bin/gs file.sad > file.result   # NPARA=1; を入れる
```

## 6. 到達点と未対応

SAD と一致しているもの (詳細は 10_benchmarks.md):
- 線形光学 (CELL/INS、CALC6D、ROTATE/ずれ付きの全要素、ANGLE=0 BEND の K0/K1、TRPT の加速と MARK OFFSET、幾何 GX..、RADINT) ≤1e-10 (誤差・オフモメンタム付きは ≤1e-7)。
- Emittance[]/EMIT (RFSW × RADCOD の 4 通り、INTRA、不安定な光学、報告の全行) ≤1e-11〜1e-9、Touschek ≤7e-12。
- マッチング (解析微分、DESCEND、fuzz、tsvd、2 段目フィット、INS、フック) は反復列・残差・K1 を SAD と同じ桁で再現。
- 追跡 (NORAD/RAD/FLUC、粒子損失) ≤5e-11、DA の地図と Score、誤差の実現 (DELX.. の乱数列)。
- SADHelp の対応 (14_pysad_vs_sadhelp.md の 11/16 版): FFS コマンド 56 中 45 実装 (未対応 2、無視 9)、関数 94 中 46 (この数は表が数える SADHelp 記載の関数の数で、特殊関数・Fit・Packages の関数 (docs/33, 32) は SADHelp の外に増えている)、特別変数 48 中 44 を参照、フラグ 60 中 59 を登録。
- 要素・機能 (11/13〜11/16 に追加。各文書の「限界」を必ず見る): COORD・TCAVI・UND・MULT ANGLE≠0・SOL の Twiss/Emittance[]/追跡・MAP・DELCOR・GaussianCoulomb・BBBrem・BEAMBEAM・WAKE・CSR・POL のスピン (Emittance[]・EMIT・TrackParticles)・RADLIGHT・SynchroBetaEmittance・特殊関数・Fit/FindRoot・SAD パッケージ (CorrectOrbit・MakeBump ほか)・FFS の小さなコマンド・描画 (docs/22〜35)。gs との一致の桁は各文書の表にある。

| 未対応・制限 | 備考 |
|---|---|
| 性能 | 誤差の実現 ×多数・DynamicApertureSurvey 5000 周・チューンサーベイ 54×54 の本番サイズは長い。1 周の追跡は 176 ns/粒子・要素 |
| FITP>1 の最終解 | SAD (WSL) が GO で再起動するため参照値がない。1 段目・初回の解は一致、2 段目以降の多解は別解に着くことがある。fam 点は自己検証のみ |
| 不安定ラティス | SAD の 1e199 級の値や NaN は再現しない (安定域は一致) |
| 損失周の乱数 | 論理は SAD と同じだが、5000 周の DynamicApertureSurvey の最終 seed はカオスで一致しない (Score の差 1〜3 %) |
| 要素 | **SOL 領域 (BZ≠0) 内の BEND・MULT・CAVI** (追跡は NotImplementedError。BEND は K64 の tsol.f どおりに書いたが gs が偏向をほぼ無視して合わず取り下げた)、SOL の CALC6D・Touschek、SOL の BZ≠0 内の四極 + RAD の追跡 (現在の gs は保存済みの参照 bench_level9 と別の値を返す。11/16 に単一粒子で切り分けた: 領域内がドリフトだけなら gs = pysad、四極が入り RAD ON かつ BZ≠0 のとき gs は四極の収束・発散力を入れていないように振る舞う (四極なしと 8 桁で一致)。gs の 2 つのビルドは同じ値で、現在の gs が SOL 領域内の四極 + RAD を誤って扱う疑いが強い。ソースのどこかは未特定、保存済みの参照を出した旧版の実行ファイルは動かず確認できない。docs/10 §6)、MAP の GEO・CompiledMap、他の要素の光子 (PHOTONS の tphrec)、POL のスピンで SOL 領域内の MULT・BEND・CAVI、PHSROT (K64 で動作しない: B38) |
| コマンド・関数 | MEASURE/TRCOD (`tracke/ttinit/ttstat/trcoda/trackb`: gs で動かない (B42) うえ使うスクリプトも無く未翻訳。`tlum` は `measure.py`)、DRAW/GEO 描画コマンド、FFS の `SYNCHROB` 語・`INC`・`ALIGN`・`SCALE`・`RENUM`・`ZAP`・`DIMAD`・`REF`・`BEAM` は未対応 (警告して無視。関数 `SynchroBetaEmittance` は実装済み)、`Polarization` (K64 が NaN: B55)、旧式の MAIN の TRACK の Fourier 解析 (`wjfit/wfres/wfbin/wiord/wjdraw`。FFS から使えない)、`MessageName[f,"usage"]` の関数形 (gs と挙動が違う。`f::usage` の形は一致)、ElementValues の VISIT/BYE の退避・復元 (CALC の前の反映は docs/34)。WakeFunction は追跡も光学 (WAKEOPT、InitialOrbits) も実装済み: [23_wake.md](23_wake.md)。旧 WAKE コマンドは K64 でも `WAKECOMMAND.n` (配布に無い) を呼ぶので同じ経路 |
| 対象外 | EPICS・`System[]` の実行、ファイル入出力の失敗の扱い (Python ネイティブで行う)。Tk・描画は tkinter で対応済み (docs/22, 35。見た目と操作は SAD と未比較) |
| 意図的な差 | SADScript の `1GEV` は pysad では単位を掛ける (SAD は `1*GEV` = 1、11_sad_quirks C24)。MAIN の `MOMENTUM=` の初期値は K64 の initb1.f に合わせて 30 GeV |
| 未解決の差 | BEND に `F1=0.1 FRINGE=3` を付けて RAD ON で追跡すると gs と合わない (相対 0.7)。L≠0 の CAVI に DX・DY を付けた追跡が 1e-3 ずれる。K3 が強くずれの大きい MULT リングで pysad の CALC (qcod) の閉軌道が gs と別の解になり Emittance が NaN 。いずれも原因未調査 (`history/work_unimpl.md`。2026-11-16 夜の確認では取り上げていない)。`emit6d._tquade` (Python 版) の解析ヤコビアンは x がずれた軌道で δ の微分を落とし、`emit_extra`・`calc6d._turn`・`touschek` で 1e-5 ずれる (`tturne` の経路は 1e-12。`history/work_unimpl.md` の記録。直したかは今回コードで確かめていない = 未確認) |

