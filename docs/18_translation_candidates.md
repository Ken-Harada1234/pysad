# 18. pysad 翻訳候補と状況 (SAD の機能別)

最終更新: 2026-11-16

目的: SAD-master (`src/*.f`、`Packages/*.n`、SADHelp) の機能ごとに、pysad での対応状況 (済み・見送りと理由) を並べる。
難易度は SAD 側のルーチンの大きさと、pysad の既存モジュールからの距離による目安。照合は `bin/gs` (SAD の WSL ビルド) の出力との比較で、数値は各文書にある範囲 (小さなリング・試験の条件) のもの。実用のリングや強い非線形での確認は別である。

試験は `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本を除く) で 575 通過・3 skip (578 件収集。GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip)、回帰は 50 ベンチをすべて実行して FAIL 0 (2026-11-16 夜、`local/wake/verify_log.md`)。

## 1. FFS のコマンド (`tffsa.f` ほか)

| 機能 | SAD 側 | pysad の状況 |
|---|---|---|
| STOP / USE / VISIT / DISP / TYPE / SAVE / CALC / CAL / FIX / EXPAND / REJ / RESET / EXIT / ABORT / BYE | `tffsa.f` | 済み (`USE BeamLine[..]` も可。反転要素名 `-NAME` の差は docs/11) |
| FIT / FREE / GO / FITP / VariableRange / VariableWeight | `tffsmatch.f`、`qdcell.f` の tffsqu | 済み (解析微分 match6d.py が既定。多解問題の 2 段目以降は SAD と別解のことがある。FITP>1 は SAD 側が再起動を繰り返して落ちるので参照値なし) |
| DP / DP0 / GCUT / COUP / COUPLE | `tffsa.f` | 済み |
| DELX / DELY / DTHETA / DK / DDK / DL | `terror.f` | 済み (直訳、照合済み) |
| DUMP | `tffsa.f` | 済み (DK・DDK 列の書式 (Fortran の 1p g12.5) を `ffscmd._gwd` で合わせた。`test_dump_format_matches_gs`) |
| EMIT (コマンド/フラグ) | `temit.f` | 済み |
| 要素値設定 (`SDA K2 …`、`MAX`/`MIN`)・FITFUNCTION・入口 Twiss (`AXI` `BXI` … `EPXI`)・fit 条件 (`NX` `NY` `BX` …) | `tffsa.f` | 済み |
| UNTIL / REPEAT | `tffsa.f:255-290` (入力位置 `ipoint` を保存して巻き戻す) | 済み (最上位の FFS 文列) |
| IF / ELSEIF / ELSE / ENDIF | `tfif.f` | 済み (同上) |
| CHROMATICITY (CHRO) | `tfchro.f` (QUAD / SEXT / MULT ごとの薄肉近似) | 済み |
| MAIN の TRACK (`TRACK USE= NPART= TURNS= NX= NY= NZ=`) | `track.f`、`tracka.f`、`ttinit.f` | 済み (`pysad/maintrack.py`。K64 の MAIN の TRACK は SIGSEGV で動かず、SADScript で同じ手順を書いた版と一致)。PHSPLOTS の描画は無視、FOURIE ON は未対応 (docs/history/work_main_track.md) |
| PRINT / REF / ORG / SHOW / STATUS / MATRIX / BEAMSIZE / VARY / RECOVER / REVERSE | `tffsa.f` の各分岐、`tfmat.f` (MATRIX) | 済み (docs/33) |
| RENUM (要素番号の振り直し)・GEO・DRAW | `HelpMessages.n:332`、`tfgeo.f` | 警告して無視 (描画・表示系) |
| RADINT | `INTGRL.f` (693 行) | 済み (`radint.py` を FFS RADINT に接続。6 通りと小リング 3 通りが SAD と 7 桁一致) |
| APERTURE / ZAP | `tfaprt.f`、`tfzap.f` | tapert (追跡側の損失判定) は翻訳済み。コマンドは未対応 (表示系) |
| DELCOR | `tcorr.f`、`tvcorr.f`、`tffsa.f:982` | 済み (gs と一致。CAVI・OCT は未対応。docs/29) |
| MON / MCOD / TCOD | | 受理のみ (K64 ビルドでは SAD も何もしない) |
| MEASURE / PUTMEASURE / GETMEASURE / LOGMEASURE / ORBITJITTER / TRCOD / FIXCOD / MVCOD / STNCOD / LTRACK / PKSEED / PUSHSEED / POPSEED / EXCGSEED / DELWAVE | `tffsa.f:639-985`、`trcoda.f`、`trackb.f` | 受理のみ (TRCOD、MEASURE) または未定義。見送り (gs で TRCOD が SIGSEGV、MEASURE は何も出力しない: B42) |

### 1.1 ON/OFF のフラグ

| フラグ | pysad の状況 |
|---|---|
| LOG / COD / CELL / RADCOD / RFSW / FLUC / RAD / RING / INS / TRPT / CODPLOT | 済み (INS の初期軌道オフセットも一致。CODPLOT は表示なし) |
| DAPERT | 済み (dapert.py、DA の地図一致。カオス的な領域の最終 seed は不一致) |
| INTRA | 済み (EMIT が IBS 付き、TouschekLifetime の表。`temit.f` の INTRA 経路) |
| CALC6D | 済み (CELL と INS、1.7e-10。`qcell.f` の qcell6d) |
| CTIME / ECHO / EMIOUT / GAUSS / CONV | フラグとして受理。効果はソースで確認 (2026-11-16): GAUSS は MAIN の TRACK・MEASURE の初期分布、CTIME は CPU 時間の表示、CONV は収束の判定・表示、EMIOUT は EMIT 報告の追加行 (pysad は spin の行だけ対応) |
| BIPOL / UNIPOL / RELW / ABSW | BIPOL は FIT だけ (tffsmatch.f)。UNIPOL・RELW は FIT の結果を変え、pysad も match6 につないで gs と全桁一致 (test_item4.py) |
| WSPAC / SPAC / SELFCOD | 受理のみ。空間電荷の追跡側 (tspac.f。SPCH 要素は tturn.f:521 が PSPAC フラグで tpspac を呼ぶ) は未翻訳。tspac.f の bbkick1 と tintrb.f の twspu だけ `GaussianCoulomb(Fitted)` として翻訳済み (docs/30) |

## 2. SADScript の関数

「組」は 組み込み (B、`src/*.f` の `itfunaloc`) / パッケージ (P、`Packages/*.n`) / ラッパー (W: パッケージ経由で組み込みを呼ぶ)。

| 関数 | 組 | SAD 側 | pysad の状況 |
|---|---|---|---|
| Twiss | B | `tfetwiss` ほか | 済み。キー NX NY BX BY EX AX AY DX DY EY EPX DPX S `*`、SIGX/SIGY/SIGPX/SIGPY/SIGZ/SIGDP/SIGE/SX/SY (emit_extra.py、2.1e-10) |
| LINE | B | | 済み。キー NAME K1 LENGTH K2 POSITION K0 L ANGLE S ELEMENT DIR DX DY TYPE TYPENAME DK VOLT ROTATE E1 E2 GX GY GZ CHI1-3 OFFSET。RHO/RHOX/RHOY/KX/SK0 は要確認 |
| FFS | B | | 済み |
| Emittance | W (emit.n) | `tfemit.f` (226 行)、`temit.f` | 済み。オプション OneTurnInformation・Matrix・OneTurnTransferMatrix・SaveEMIT・Output は emit_extra.py (2.1e-10)。偏極・SpinTune は 0 固定 |
| DynamicApertureSurvey | P (dapert.n) | | 済み (遅い。既定の Output / ExpandElementValues は ffscmd.py:880 で再現) |
| TrackParticles | B | `tftrack.f` | 済み |
| CalculateOptics | B | `tfoptics.f` (117 行。初期条件 28 個を行 lb に書き込み qcell1 を実行) | 済み (4 次元の cell/INS・区間・結合・分散、idim=3 は全周) |
| Monitor / Steer | P (corutility.n) | `Master[]` (pgmast.f) と `LINE["type"/"element",…]` を呼ぶだけ | 済み (SAD の corutility.n を pysad のインタプリタで動かす。docs/32。FFS の MON とは別物) |
| ExtractBeamLine / BeamLine / PrintBeamLine / WriteBeamLine / ExpandBeamLine | B (tfbeamline.f:362-510) / P (beamline.n) | 逆向き要素は `-1*名前` | 済み (`USE BeamLine[..]` / `VISIT` / `BYE` も) |
| SetElement / Element | W (element.n → SetElement$ / Element$) | `tfbeamline.f:10-340` | 済み。`Element["KEYWORDS"/"DEFAULT"/"VALUE"/"POSITION"]`、型指定・新規定義とも。PROFILE, Saved->True, DK の連動を除く |
| FLAG | B (pgflag.f) | | 済み |
| GetMAIN | B (tffs.f tfmain) | | 済み (MAIN 形式ラティスの読込み) |
| FFSFlags | B (tfetok.f の tfflags) | | 登録名あり。内容は未対応 |
| CorrectOrbit / MakeBump | P (correction.n、bump.n) | `pgrmat.f` (応答行列)、`pgsolvcond.f` (条件付き SVD) | 済み (SAD の .n を動かす。ResponseMatrix・LinearSolveConditioned・Master は `corr.py` に直訳。補正量は gs と相対 1e-7。docs/32) |
| OptimizeOptics / DownhillSimplex | P (optimize.n) | | 済み (K1・チューンが gs と 1e-10。DownhillSimplex は `miscfn.py` で gs と 1e-13。docs/31・32) |
| TransferMatrix / BeamMatrix / NormalCoordinates / SymplecticInverse | P (Optics.n、emit.n) | | 済み (optics_fn.py、SAD と 1e-9) |
| PolynomialFit / Fit / FitEmit / FitGaussian / Spline / NIntegrate / Integrate / FindRoot / D | P (Fit.n、D.n) | | 済み (`fitroot.py`。gs と 5e-10〜1e-14。ベクトル独立変数と複素 FindRoot は未対応。D・Integrate は D.n を動かす。docs/31・32) |
| TouschekLifetime / TouschekTable | P (emit.n:129-) | `temit.f` の intra | 済み (15 ケース ≤7e-12) |
| SynchroBetaEmittance | B (`tfemit.f:144` の tfemits) | `temits.f` (1325 行) | 済み (`synchrobeta.py`、gs と相対 1e-5。24 セルの平面リング + CAVI。結合のあるラティスは gs が SIGSEGV で比較不能。不安定なチューンは未確認) |
| SeedRandom / GaussRandom / Random | B | | 済み (SAD と同じ乱数列) |
| LDSolve ほか (ExpMatrix、rational.n、extmap.n、Misc.n) | P | | SAD 自身の `.n` を解釈する方式で動く (docs/32)。LDSolve は複素数・行列演算の虚部の扱いを直して gs と一致。他は個別の照合なし |
| Polarization / SpinTune / PolarizationTime | | `temit.f:491,1188` の calpol、`polpar.f`、`photons.f`・`tturn.f` の spin 部 | 済み (Emittance[]・EMIT のスピン (MULT・空洞・SOL まで) と、POL かつ RAD ON の TrackParticles のスピン追跡。gs と 1e-11〜4e-12)。`Polarization` は K64 が NaN (B55) なので値の照合なし (docs/34・41) |
| PHOTONS (`PhotonList`) | | `tftrack.f:144,315`、`photons.f` の kphlist | 一部対応 (BEAMBEAM の光子だけ。他の要素の光子 tphrec は未対応。RADLIGHT の軌道と RadiationField・RadiationSpectrum は docs/34・43) |
| LMA | | | SAD に組み込みはない |

### 2.1 SADScript 言語系の組み込み

`pysad/sadlist.py` に MapThread, ScanThread, MapAt, Complement, DeleteCases, Cases, MatchQ, Extract, Intersection, Union, Position, Count, MemberQ, FreeQ, Select, Sort, Take/Drop, Insert/Delete/ReplacePart, Apply, Map, Fold, Nest, FixedPoint, Through ほかを実装済み (SAD 側は `tfinitn.f` の `itfunaloc` 登録と tfefun*.f / tfeeval.f のディスパッチ)。Total/Cross/Norm/Riffle/StringSplit は SAD にも無い (未評価のまま返る。gs で確認)。Sleep は no-op、Fourier / InverseFourier / BesselJ / BesselK は FFT / Bessel で実装。

### 2.2 パッケージ経由のものの注意

Emittance / Element / SetElement / Monitor / Steer / DynamicApertureSurvey / TouschekLifetime / PolynomialFit は SAD では SADScript パッケージ (`Packages/*.n`) の関数で、Fortran の組み込みではない。Python で個別に書き直したものは、SAD のパッケージの既定値・オプションとの差が出やすい。`pkgload.py` が SAD の `Packages/*.n` を最初に使うときに読む (環境変数 `PYSAD_SAD_PACKAGES` が要る。docs/32)。

## 3. 要素・物理モジュール

| 機能 | SAD 側 | pysad の状況 |
|---|---|---|
| BEND の K0 キック (steering) | `tsteer.f`、`tbende.f`、`tbend.f` | 済み (CALC の tsteee、追跡の tsteer) |
| 幾何 (LINE["GX"/"GY"/"GZ"/"CHI1".."CHI3"]、fit の GX) | `tfgeo.f` (411 行)、`tdgeo.f` | 済み (geo.py、4.8e-13。fit の幾何の微分は、DRIFT の L と BEND の ANGLE・L だけ解析微分 (tdgeo)、ほかは中心差分) |
| MULT ANGLE≠0 | `tmulta.f` (433 行) | 済み (docs/26) |
| SOL | `tsole.f`、`tsolque.f`、`tsgeo.f` | Twiss・Emittance[]・EMIT・追跡 (領域内は DRIFT・QUAD のみ)、幾何は gs と一致 (docs/27・42)。K64 は SOL (BOUND≠0) の長さを周長に入れない (B63) ので同じにした |
| TCAVI | `ttcav.f` (72 行) | 済み (docs/24)。L≠0 の Emittance[] が gs で壊れる: B36 |
| COORD | `tcoord.f` (82 行) | 済み (docs/23) |
| MAP | `temap.f` (625 行) | 済み (docs/28)。`ExternalMap["GEO"]` は未対応 (直進、長さ 0) |
| UND | `undulator.f` (165 行) | 済み (docs/25)。磁場は K64 に無い: B37 |
| WAKE | `tffswake.f`、`txwake.f` | 追跡と光学 (WAKEOPT) が gs と一致 (docs/21・36)。旧式の WAKE コマンド (tfwake.f) は K64 が WAKECOMMAND.n (配布に無い) を呼ぶので同じ経路で、移植不要 (B66) |
| BEAMBEAM | `bbstrhl.f` | 強弱ビームが gs と一致 (docs/22・37) |
| CSR | `csroy.f` ほか、CSR.n | CSRTrack・CSRHaissin・CSRSetupOY・CSRMatrix (mode 1000 = csrhelmat を含む)・CSRInit・CSRConvert、MakeZL の 3 周波数 (docs/34・38)。ScanTrack・OYStability は未照合 |
| RADLIGHT・放射場 | | TrackParticles が軌道 {t,x,y,z} を返す (gs と差 0)。RadiationField・RadiationSpectrum を実装 (RadiationField は gs が落ちる: B49) (docs/34・43) |
| GaussianCoulomb(Fitted)・BBBrem | `tspac.f` の bbkick1、`tintrb.f` の twspu、`tbbbrem.f` | 直訳 (docs/30) |
| ルミノシティ | `tlum.f`、`tluma` | `measure.py`。MEASURE が gs で動かないので解析値との比較だけ |
| 空間電荷の追跡側 (SPCH 要素) | `tspac.f`、`tturn.f:521` | 未翻訳 (SPCH は読むだけ)。PROF は K64 に無い |

## 4. GUI・描画

Tkinter.n・Plots.n・OpticsPlot.n などを SAD の .n のまま tkinter で動かす (docs/20・35・39・40)。BeamPlot は未比較。CanvasSymbol・Canvas3D* は未着手。

## 5. 見送り・未対応とその理由

| 項目 | 理由 |
|---|---|
| MEASURE / TRCOD (`tracke`・`ttinit`・`ttstat`・`trcoda`・`trackb`) | gs で TRCOD が SIGSEGV、MEASURE は何も出力せず参照が取れない (B42) |
| `wjfit`・`wfres`・`wfbin`・`wiord`・`wjdraw` | MAIN の TRACK の Fourier 解析で、FFS から使えない |
| SOL 領域内の MULT・BEND・CAVI と tsgeo / tdgeo の領域幾何 | BEND は K64 の tsol.f どおりに実装したが gs が偏向をほぼ無視して合わず取り下げ (B69)。MULT・CAVI は未着手 |
| `tspac.f` の追跡側 (空間電荷) | K64 の tturn.f:521 が PSPAC フラグで tpspac (粒子のメッシュによる空間電荷) を呼ぶ。pysad は未翻訳。bbkick1 と twspu だけ翻訳 |
| LOSSMAP | gs は TrackParticles の戻りに 2 行 (落ちた要素と周) を足すが、pysad は未実装 |
| `Polarization` | K64 が NaN を返す (B55) |
| PHSROT | K64 で動作しない (B38) |
| `MessageName` の関数形 | gs は f を無視する (B70)。`f::usage` の形は一致 |
| 複素引数の Bessel・超幾何、PolyLog の \|z\|≥150 (非整数次) | 未着手 |
| DELCOR の CAVI・OCT、MAP の GEO・CompiledMap、BeamPlot の比較 | 未着手 |
| FFSFlags の内容 | 未対応 |
| SAD の MON / MCOD / TCOD | K64 でも無効 (SAD 側で参照が取れない) |
| FITP>1 の GO | SAD が再起動を繰り返して落ちる (参照が取れない) |

## 6. 未解決の差 (pysad と gs が合わない。原因は確定していない)

BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡 (相対 0.7)、L≠0 の CAVI の DX・DY (相対 1e-3)、SOL (BZ≠0) 内の四極 + RAD の追跡で現在の gs が保存済みの参照 (10-14) と別の値を返す (B68)、
K3 が強くずれの大きい MULT リングで pysad の CALC (qcod) の閉軌道が gs と違う解になり Emittance が NaN (`local/wake/pm1.sad`)。詳細は `docs/history/work_unimpl.md`。
