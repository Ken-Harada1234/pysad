# 46. オリジナル SAD (K64) に対する pysad の未対応・未確認の一覧

最終更新: 2026-11-17 (現状のスナップショット)

この文書は、オリジナルの SAD (K64 版、`src/*.f`・`Packages/*.n`・SADHelp) に対して、pysad が**まだ対応していないもの**と**対応したが gs で確かめていないもの**を、1 か所にまとめたものです。

- 出所: [16_fortran_status.md](16_fortran_status.md) §9、[history/work_unimpl.md](history/work_unimpl.md)、[21_sad_bug_report.md](21_sad_bug_report.md)、`example/` の各例の「gs と pysad の比較」・「未確認・限界」節、[34_sad_packages.md](34_sad_packages.md) ほか各解説文書。
- この一覧のために、改めて pysad を走らせて確かめたものではありません。各項目の根拠は出所の文書にあります。
- 基準は gs (SAD K64) の実動作です。gs が落ちる・誤った値を返すものは「gs 側の問題」に分けています。
- `例` の欄は、その項目を示す `example/` の例です (パスは `example/` からの相対)。

凡例: 「未翻訳」は pysad に物理が入っていない。「一部未対応」は一部の入力・経路だけ未対応。「未確認」はコードはあるが gs と数値で比べていない。

## 1. 未翻訳

| 項目 | 理由・状況 | 例 |
|---|---|---|
| MEASURE / TRCOD (`trackb.f`、`tracke.f`、`trcoda.f`、`ttinit.f`、`ttstat.f`) | gs で TRCOD が SIGSEGV、MEASURE は出力なしで、参照が取れない。`tlum`・`tluma` だけ翻訳 (`measure.py`) | |
| `wjfit`、`wfres`、`wfbin`、`wiord`、`wjdraw` | MAIN の TRACK の Fourier 解析で、FFS から使えない | |
| SOL 領域内の MULT・BEND・CAVI (`tsole1`、`tmulte1`)、SOL の `DX DY DZ CHI`・`F1`、領域の `tdgeo` | gs が領域内の偏向をほぼ無視するため、比較できない | 07_spin_solenoid_special/elem_SOL |
| `SPCH` 要素と `tspac.f` の追跡、`PSPAC` フラグ | 読むだけで、空間電荷の物理は未翻訳。gs でも効果が確認できない例がある | 05_tracking/flag_SPAC |
| `LOSSMAP` | gs は `TrackParticles` の戻りに 2 行 (落ちた要素と周) を足すが、pysad は足さない | 05_tracking/flag_LOSSMAP |
| `RADTAPER` | フラグを保持するだけで、磁石の値が変わらない | 04_emittance_radiation/flag_RADTAPER |
| `System`、`Fork`、`Shared` | 登録だけの限定版。`System` は実行せず、`Fork` は子プロセスを作らない | 09_io_system/fn_System_TemporaryName、fn_Fork_Shared |
| `BeginPackage` / `Begin` の文脈 (名前空間) | gs の動作は調べたが、実装が `.n` パッケージの読み込みを壊す恐れがあるため見送り | 08_scripting/fn_Get_Needs_BeginPackage |
| `tphrec` (光子) | `BEAMBEAM` の分だけ対応。ほかの要素の光子は未対応 | |
| `tdgeo.f` の四極などの DX・DY の解析微分 (`qdtrans`・`tftmatu`) | 数値微分で代用 | |
| `PHSROT` | K64 で動作しない。追跡は NotImplementedError | |
| `Polarization` の有限値 | K64 が NaN を返す。pysad も NaN のまま写す | |
| `NOKEEPEXP`、`RENUMBER` の `INDX` 初期値、`Out[n]:=` の表示 | 未対応 | |

## 2. 一部未対応

| 機能 | 未対応の部分 | 例 |
|---|---|---|
| MULT `ANGLE`≠0 | Twiss の `BZ`・`RF`・`F1`・`SK0`・`FB1/FB2` (追跡は `FB1/FB2` まで対応) | 07_spin_solenoid_special/elem_MULT_ANGLE |
| MAP | `GEO`・`CompiledMap`。MAIN の TRACK・DA・Touschek は ExternalMap を呼ばない | 06_collective/fn_ExternalMap |
| DELCOR | `CAVI`、`OCT` | |
| POL | SOL 領域内の MULT・BEND・CAVI のスピン。スピン追跡は損失の処理なし | 07_spin_solenoid_special/flag_POL_RADPOL |
| TCAVI | `RANKICK`・`RANPHASE`。`Emittance[]` は L=0 のみ gs と一致 (L≠0 は gs が壊れる) | |
| UND | 磁場なしのドリフトのみ (K64 が磁場を設定できない) | |
| `APERT` | SOL 領域内は NotImplementedError。再入射 (`JDPX`・`JDPY`) は粒子 12 個以上で gs と乱数列が違う | 05_tracking/elem_APERT_LOSSAMPL |
| RF のある MULT (`tmultiacc`) | `F1`・`FRINGE`・`AUTOPHI` | |
| `ResponseMatrix` | TRPT (加速) の経路 | |
| `tsteer` | `F1/FB1/FB2` つきは NotImplementedError | |
| `FindRoot`・特殊関数 | 複素 `FindRoot`、複素引数の Bessel・超幾何 | 10_math_data/fn_Gamma_Erf |
| `ElementValues` (`Couple.n`) | `VISIT` / `BYE` の退避・復元 | 13_packages/fn_Couple |
| 描画 | `CanvasSymbol`、`Canvas3D*`、`UseBG->True` の文字背景。`OpticsPlot` の text 項目数 (gs 69、pysad 74。Tk のフォントによる `RemoveOverlap` の差)。`BeamPlot` (gs は `SIGXX=0` で楕円を描かない) | 11_graphics/fn_OpticsPlot、fn_BeamPlot |
| GUI | 最初の `Window[]` の親 (gs は Toplevel、pysad はルート窓)。起動時の窓の再現が要る | 12_gui_tkinter/README |
| `Packages/CSR.n` | 全体は読み込まない (必要な部分だけ) | 06_collective/fn_CSRInit_CSRTrack |

## 3. 未確認 (コードはあるが gs と数値で比べていない)

| 項目 | 状況 | 例 |
|---|---|---|
| `FITP>1` の最終解、`tffamsetup` (fam 点) | gs が `GO` で再起動することがある。自己検証のみ | 03_matching/FFS_FITPOINTS |
| `CSROYMatrix` | gs が wrongtype を返す。比較なし | |
| `tbrote` の `srot` 更新、`tbendebody` の `bsi2` | ソースの読みだけで、gs の個別確認なし | |
| `RadiationField` | gs が SIGSEGV。速度場との比較のみ | 07_spin_solenoid_special/flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum |
| `SynchroBetaEmittance` | 平面リングのみ照合。結合のあるラティスは gs が落ちる | 04_emittance_radiation/fn_SynchroBetaEmittance |
| `measure.py` (`tlum`・`tluma`) | gs で MEASURE が動かない。解析値のみ | |
| `tfloor.f` の `cfunc` (`xsincos` など) | 個別実装 | |
| `Optics.n` の `PeriodicR`、`CorrectOrbit` の Micado・`CorrectLocal` | 未確認 | 13_packages/fn_Optics_TransferMatrix、fn_CorrectOrbit_Monitor_Steer |
| `BeamPlot`・`ColumnPlot`・`LaTeX`・`PointSymbol` の見た目 | 動くことだけ確認。LaTeX は `$HaveLaTeX` が空の状態での一致 | 11_graphics |
| Tkinter の見た目・操作 | 色・大きさ・フォントは SAD と比べていない | 12_gui_tkinter |
| `dapert.n`、`element.n`、`emit.n`、`beamline.n`、`functions.n`、`Misc.n` | ファイル単位で gs と比べた記録がない | |
| 乱数列 (FLUC つきの追跡) | δ の rms が gs 4.21e-4、pysad 4.52e-4 (乱数列が違う) | 04_emittance_radiation/flag_RFSW_RAD_RADCOD_FLUC |

## 4. 既知の差 (未解決。原因は未調査)

| 差 | 状況 |
|---|---|
| BEND の `F1=0.1 FRINGE=3` + RAD の追跡 | 相対 0.7 の差 |
| L≠0 の CAVI + `DX`・`DY` | 1e-3 の差 |
| SOL 内の四極 + RAD | gs が RAD ON・BZ≠0 の領域内で四極の力を入れていない疑いが強い (gs 側の欠陥) |
| EMIT のあとに `USE` する | gs は EMIT の効果が追跡に残り、pysad は捨てる |
| `RADINT` と `EMIT` の Jx | 0.9415 と 0.9651 (2.4% の差)。I4 の端の項を 2 倍にすると近づく。pysad 側の原因は未調査 |
| `BUCKET HIGHT` (`RADINT`) と `BucketHeight` (`Emittance[]`) | 0.010045 と 0.010035。原因は不明 |
| `pm1.sad` の NaN、複数セグメントの `RadiationField` の失敗 | 未調査 |
| `LOSSAMPL`・RAD・FLUC つきの TRPT、MAIN の TRACK の損失報告、DAPERT での APERT | gs と比べていない |

## 5. gs 側の問題で比べられないもの

| 項目 | 状況 |
|---|---|
| `For` | gs は 4 番目の引数 (本体) を実行しない。K64 の癖と判断し、pysad は本体を実行する ([21_sad_bug_report.md](21_sad_bug_report.md) B71) |
| Bessel・HermiteH・Ei の一部 | gs が誤った値を返す (`BesselJ`(100) など)。pysad は触っていない |
| gs が落ちる | `RadiationField` の SIGSEGV、`TRCOD`、`TCAVI` の L≠0、`SynchroBeta` の結合ラティス、`SPAC` の小さい `RADIUS>0` |
| `draw` コマンド | gs は表示環境がないと出力せず、`FFS["draw ...",12]` は内部エラー。pysad は無視する |

## 6. 優先度 ([16_fortran_status.md](16_fortran_status.md) §9.4 から)

| 優先 | 項目 | 理由 |
|---|---|---|
| 中 | 誤差入り 6D emit と DA の反復 (性能) | 反復が 600 s を超える |
| 中 | SOL 領域内の MULT・BEND・CAVI、CALC6D・Touschek | 要素の完全対応 |
| 中 | §4 の未解決の差 | 原因の切り分け |
| 低 | `tdgeo` の四極の微分、`wfres` 系、`tspac` の追跡 | 必要になったら |
| 不要 | MON/MCOD/TCOD、Obsolete・非リンクのファイル、旧 WAKE | K64 で無効 |

## 7. 注意

- 条件: 例の数値は、サンプルのラティス (`example/00_lattice/FODO_ring.sad`: 3 GeV、周長 105.6 m) と、PF・ASTRID2 のラティス (`example/pf`、`example/astrid`) で得たものです。
- 表示や数字の桁の微妙な差は、ここには入れていません。完走することと、意味のある値が合うことを基準にしています。
- 一覧は変わります。項目を直したら、この文書と [16_fortran_status.md](16_fortran_status.md) §9 を一緒に更新してください。
