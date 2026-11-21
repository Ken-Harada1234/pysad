# 元のスクリプトとの対応 (v07_14)

新しい名前 ← 元のファイル名 (元のフォルダは読むだけで、書き換えていない)、変更点、元の設定と例の設定。元のスクリプトは拡張子が無く、ラティスの台本とスクリプトを 2 つ並べて `gs ラティス スクリプト` と走らせる使い方 (スクリプト単体ではラティスが無い) で、同じ名前の部分が版 (v07_12、v07_14、v07_16、v08_09 …) ごとに少しずつ違う複製になっている。ここでは v07_14 のラティスに合うものだけを選んだ。

## 元のファイルの特定 (v07_14 に関するもの)

| 役割 | 元のファイル | 判断 |
|---|---|---|
| ラティス (5 GeV) | `v07_14` | 6 極を 2 回定義していて、あとの定義 (20231205 の再調整版) が有効。`RING` を 3 回定義していて最後 (SP 輸送の並び) が有効 |
| ラティス (2.5 GeV) | `v07_14_2.5GeV` | `MOMENTUM = 2.5GEV`、RF 1.6 MV。6 極は**古い組だけ** (再調整版の定義が無い)。`RING` は最後の定義 (通常の入射用の並び)。6 極の再調整版は `v07_14_(人名)_20231205_2.5GeV` (同名の 20231205 版。個人名の部分は伏せる) に入っている |
| ラティスの複製 | `v07_14_(人名)_20231205`・`_5.0GeV`・`_2.5GeV`、`v07_14_SP_transport`、`v07_14.fin`、`v07_14_csr*` | 20231205 版は6 極の再調整版を整理した複製 (定義の順序を直したもの)。`v07_14_SP_transport` は再調整前の6 極。`v07_14.fin` は ELEGANT の出力 (SDDS) でラティスではない。`v07_14_csr*` は CSR (ELEGANT) 用。いずれも取り込まない |
| チューンのマッチング | `matching_v07_11` | `nx0 = 47.865`、`ny0 = 16.655`。v07_14 のチューン (47.865, 16.655) に一致する。`matching_v07_13` は `nx0 = 45.865`、`ny0 = 21.665` (チューンの打ち間違いの版。内容はほかは同じ)。`matching_v07_15_Chicane_ANG` は別のラティスの台本 (TCELLBP など)。取り込まない |
| 6 極の色収差補正 | `chrom_v07_12rev` | S1・S2 を 0 にして S3・S4 だけで合わせる。`chrom`・`chrom2`・`chrom3`・`chromA`・`chromN`・`chrom_ESSR` は別のリング用 |
| チューンサーベイ | **v07_14 用は無い** | `tune_sv_v06_10`・`_parm`・`_sx_err`・`_temp` は 28 セルの旧い版 v06_10 用 (`LCELLA`・`LCELL` を 28 セルとして `nx0/28` に合わせる)。`modules_tune_sv.n` は別の旧いリング (`BPI6`・`QFA1` など) 用で、v07_14 には適用できない。`s1_s2_sv`・`calc_ts_poly_sv*` は6 極 S1・S2 の振幅依存チューンシフトのサーベイ (旧い版の名前)。→ `survey/V0714_tune_survey` は、`tune_sv_v06_10_parm` の枠組みと `matching_v07_11` のマッチングを組み合わせた**新規の台本** |
| 誤差の計算 (5 GeV) | `correct_cod_v07_14_grad_res` | 出力ファイル名が `datas*_v07_16_1.txt` (v07_16 のコピーの名前のまま) |
| 誤差の計算 (2.5 GeV) | `correct_cod_v07_14_grad_res_2` | 5 GeV 版との違い: `ae = 5E-5/2`、出力ファイル名 `datas*_v07_14_2.txt` だけ |
| 別の版の誤差の計算 | `correct_cod_v07_12_grad_res(_2)`・`correct_cod_v07_16(_2)`・`correct_cod_v08_02_grad_res`・`correct_cod_v08_09_grad_res(_neg)` | v07_14 以外のラティス用 (手法は同じ)。`correct_cod01*`・`correct_cod02_RMS` は旧い版。取り込まない |
| COD 応答・補正の関数 | (各 `correct_cod_*` の中に直接) と `modules_write_correctdata.n` | `modules/modules_errors.n` にまとめた |
| 放射積分 | `calc_RadInt.sad` | |
| クロマティシティ | `calc_chromaticity.sad`・`amp_dep_ts_dp.sad` | `calc_chromaticity_v06` は別の版 |
| 振幅依存チューンシフト | `modules_amp_dep_ts.n`・`amp_dep_ts_x.sad`・`_y.sad` | `amp_dep_ts_xA`・`_x_octa`・`_y_octa`・`_x_poly`・`_y_poly`・`_dp_poly` は八極・多項式フィット用の別版 |
| 表示 | `disp`・`dispA` | `disp_error` は `NCELLA` (v07_14 に無いライン) を使う旧い版。`disp_*` の残りは別の用途 |
| 使われないモジュール | `modules_correctoptics.n`・`modules_correctoptics_pf.n`・`modules_matching_noID.n`・`modules_general.n` (一部) | `modules_correctoptics.n` (ETA・COD・β の応答行列と、固有値を選ぶ補正 `UseTevecoCST`・`UseDiscoCST`) は、v07_14 のどの台本からも読み込まれない (別の機会の分散・β補正用)。`modules_matching_noID.n` は ID 入りの別のリング用で、`matching_v07_11` は Get するが中の関数を使わない。`modules_general.n` は `DAsv[]`・`WriteGraphs[]` だけで、`DAsv[]` は誤差の計算の中に別の版を持つ |

## ラティス・モジュール

| 新しい名前 | 元 | 変更点 |
|---|---|---|
| `lattice/v07_14_5GeV.sad` | `v07_14` | 先頭に説明のヘッダを付けた。古い6 極の組 (あとで上書きされる) を除いた。使われない `RINGL`…の定義は残し、`RING = (RINGSP)` として元の有効な `RING` (SP 輸送の並び) と同じ並びにした。ほかのリングの並び (`RINGINJ` = 元の 2.5 GeV 版の `RING`、`RINGSYM` = `HRING RFC -HRING`) を名前を付けて足した。要素の定義の中身は元のまま。`cell calc emit` で元の台本と周長・チューン・εx・クロマティシティが全桁一致することを確認 (gs) |
| `lattice/v07_14_2.5GeV.sad` | `v07_14_2.5GeV` の光学 + `v07_14` の6 極 | `MOMENTUM = 2.5GEV`、`VOLT = 1.6e6`、`RING = (RINGINJ)` (元の 2.5 GeV 版の有効な並び)。**6 極は元の `v07_14_2.5GeV` ではなく、最終設定 (20231205 の再調整版) にした** (元の 2.5 GeV 版は古い6 極の組だけを持つ。再調整版は 2.5 GeV 用のファイル (`..._20231205_2.5GeV`) にもあり、その RF は 1.6 MV)。この違いで cell calc の光学は同じで、クロマティシティ・DA だけが変わる |
| `modules/modules_amp_dep_ts.n` | `modules_amp_dep_ts.n` | コメントアウトされた旧い行を除いた。式・順序は元のまま。`ResetInit[]` を追加 (使い終わった始点の初期座標を 0 に戻す。元は `amp_dep_ts_y.sad` の末尾に直接書いていた) |
| `modules/modules_errors.n` | 各 `correct_cod_v07_14_grad_res` の関数 + `modules_write_correctdata.n` の `ResultWrite`・`WriteKick`・`RMSWrite` | `CodCorr`・`CalculateCODResponse`・`RMSDisp`・`DAsv` は元のまま (`DAsv` は δ の点数と周回数をパラメータにし、結果を 1 行に並べた。`CalculateCODResponse` の `ICX` 回の無駄な `FFS["cell calc"]` を除いた。結果は変わらない)。`cell calc` を `FFS[cc]` (`cell calc noexpand`) に変えた (後述の修正 1)。`RMSWrite` の列「チューンのずれ」は、元は旧い設計のチューン 48.6・12.7 を引いていたので、設計のチューン (47.865・16.655) に直した。`BetaWrite`・`PrintResultBeta`・`MeasureBeta`・`ReadKick` は使われないので取り込まない |

## 設計・マッチング (`design/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 | 元の設定 | 例の設定 |
|---|---|---|---|---|---|
| `V0714_match_tune` | `matching_v07_11` | 8 | `Get[...]` 3 本 (使われない関数) とコメントアウトした旧い試行を除く。7 つのリングの δ スキャンの繰り返しを関数 `ChromFitRing[]` にした (画面表示は元と同じ)。`FitValue` の定義は、元と同じく各 fit 条件の直後 (先に定義すると収束しない)。パラメータ化: `nitA` (5)・`nitN` (20)・`doR56` (1)。元の typo (`MKMSA BY` の上限が `bymklsa` など) は挙動を変えないためそのまま | `nitA = 5`、`nitN = 20` | 同じ (gs 8 s、pysad 約 90 s で足りる) |
| `V0714_chrom_cell` | `chrom_v07_12rev` | 8 | 4 つのセルの同じ処理を関数 `CellChrom[]` にした (画面表示は元と同じ)。`datafile`・`stop` の後ろの到達しない部分を除く | — | 同じ |

## チューンサーベイ (`survey/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 | 元の設定 | 例の設定 |
|---|---|---|---|---|---|
| `V0714_tune_survey` | `tune_sv_v06_10_parm` の枠組み + `matching_v07_11` のマッチング (**新規**。v07_14 用の元の台本は無い) | 32 | 格子 (nx0, ny0) の二重ループ、設計値からの再マッチング、6 極の補正、ξ (中心差分)、ADTS、ファイルへの書き出しを、v07_14 の 2 セル単位 (TCELLA・TCELL、セルチューン nx0/15・ny0/15) に直した。ADTS は x・y に 5σ、位相の分割数 `thetan`。マッチングの結果は `SAVE ALL` で書き戻す | 元 (v06_10): nx0 = 43.11〜47.12、ny0 = 10.12〜14.13 (どちらも 0.1 刻み、41 × 41 = 1681 点)、位相 36 分割。v07_14 版の完全な設定: 設計のチューン ±2、0.1 刻み (41 × 41)、`nitA = 5`、`nitN = 20`、`thetan = 36` | 設計のチューン ±0.1、0.1 刻み (3 × 3 = 9 点)、`nitA = nitN = 2`、`thetan = 8` |

## 誤差の計算 (`errors/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 | 元の設定 | 例の設定 |
|---|---|---|---|---|---|
| `V0714_correct_cod` | `correct_cod_v07_14_grad_res` | 20 → 8 | 5 段階 (誤差を 1/5、2/5、…、1 倍) の同じ処理を段階の繰り返しにした (`SetPosErr[]`)。4 段階目のゲイン 0.6 (ほかは 0.5) は元のまま。コメントの「ae 20um」などの表示は段階の表示に直した。`codplot` は走らせない (結果は変わらない)。**修正 1** `USE RING NOEXPAND`・`cell calc noexpand`。**修正 2** VD の位置ずれを 0 に戻す (README の「元の台本の問題」)。**修正 3** 書き出し (COD・K0・統計) を DA の前に移した。ファイル名は `V0714_correct_cod_5GeV_*.txt` に。DA は x・y 方向、δ の点数・周回数をパラメータに | 種 100、段階 5、補正 5 回 (最後の段階は 50 回)、DA 5000 周、δ −5%〜+5% の 11 点 | 種 2、補正 3 回 (最後 10 回)、DA 200 周、δ 3 点 (−5%・0・+5%)、NPARA 8 |
| `V0714_correct_cod_2.5GeV` | `correct_cod_v07_14_grad_res_2` | 20 → 8 | 上と同じ。`ae = 2.5E-5`、ラティスは 2.5 GeV・RF 1.6 MV | 同上 | 同上 |
| `*_full.sad` | (上と同じ) | | `tools/mk_ddba_v07_14.py --full` で作る完全な設定の版。走らせてはいない (見込み時間は README) | | |

## パラメータ表示 (`display/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 | 元の設定 | 例の設定 |
|---|---|---|---|---|---|
| `V0714_ring_params` | `disp`・`dispA` | 32 | `draw`・`TkWait` (窓) と、v07_14 に無いライン (`RINGV`・`TESTV`) を使う後半を除き、MCELL・LCELL・MCELLA・LCELLA・TCELL・TCELLA・RINGSP・RINGINJ・RINGSYM の表示を足した | — | — |
| `V0714_radint` | `calc_RadInt.sad` | (4 をコメント) | コメント行と `datafile` を除く。式は元のまま | — | — |
| `V0714_chromaticity` | `calc_chromaticity.sad` と `amp_dep_ts_dp.sad` の δ スキャン | 8・64 | 七つのリングを関数 `ChromRing[]` にして、最終設定の6 極 ON と6 極 0 の両方で出す。**中心差分の ξ を足した** (元の PolynomialFit の係数 [1,2] は K64 のバグ B18 で ξ にならない)。元の末尾の `nx0 ny0 EMITX … emitring` の `Print` (未定義の変数) は除く | δ = ±1% (RING は ±0.5%) の 3 次フィット、δ スキャン 0〜±5% (0.5% 刻み) | δ スキャン 0〜±5% を 1% 刻み |
| `V0714_amp_dep_ts` | `amp_dep_ts_x.sad`・`amp_dep_ts_y.sad`・`amp_dep_ts_dp.sad` | 64 | 3 本を 1 本に。`ResetInit[]` を追加 | 振幅 1〜101σ を 10 刻み (11 点)、位相 36 (x)・10 (y) 分割 | 3 点 (1・51・101σ)、位相 8 分割 |

## 取り込まなかった元のスクリプト

全部の中身を読んだわけではない (約 280 ファイル。名前と、読んだ範囲から判断した)。
- `matching`・`matching_v06_*`・`matching_v07_02`〜`_09`・`_13`・`_15_Chicane_ANG`・`matching_v08_*`: 別の版・別のラティスの台本。`matching_v07_13` は v07_14 と同じ中身でチューンだけ誤り。
- `tune_sv_v06_*`・`sx_sv_*`・`s1_s2_sv`・`calc_ts_poly*`・`calc_da_file*`・`modules_tune_sv.n`: 旧い版 (v06・別のリング) のチューンサーベイ・6 極サーベイ・DA の集計。
- `correct_cod01*`・`correct_cod02_RMS`・`correct_cod_v07_12_*`・`correct_cod_v07_16*`・`correct_cod_v08_*`: 別のラティスの誤差の計算 (v07_14 と同じ手法)。
- `single_pass_cod*`・`design_example.sad`・`calc_R56_pfbt`・`calc_ave_beta`・`calc_zdif`・`touscheklife`・`sad2ele_*`・`csr_*`・`disp_*`・`draw_cad`・`drawmagvw3d.n`・`make_ID`・`make_sx_mult`: 別の用途 (単発の入射、IDの作成、Touschek、ELEGANT、描画)。
- `v05_and_older`・`DDBA` のサブフォルダ・`v06_*`・`v07_01`〜`v07_13`・`v07_15_*`・`v07_16_*`・`v08_*`・`v07_14_csr*`・`elegant`: ほかの版のラティス、ELEGANT 用の複製。
