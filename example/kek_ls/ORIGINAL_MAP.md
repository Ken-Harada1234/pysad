# 元のスクリプトとの対応

新しい名前 ← 元のファイル名 (元のフォルダは読むだけで、書き換えていない)、変更点、元の `NPARA` (並列数。例は 8 以下に直した)。元のスクリプトはファイル名に拡張子の無いものが多く、`USE` するライン名やモジュール名、要素の名前が版ごとに違う。ベースは設計版 `v03_68`。

## ユーザーの指定と実際のファイル名の違い

| ユーザーの指定 | 実際のファイル | 判断 |
|---|---|---|
| `matching_02_geo` (入射セルを含めてリングを閉じる) | `matching2_02_geo` (ほかに `matching16_geo`・`matching20_geo`。20_geo は `RINGI` を使う同じ内容) | ファイル名の誤り (`matching2_`)。→ `design/KEKLS_match_ring_geo` |
| `matching23_v03_35_mini` (ノーマルセルのチューンサーベイ) | `matching23_v03_35_mini` は全体を一様に縮小する台本。チューンサーベイではない | ノーマルセルのチューンサーベイは、ファイルの中身を読んで判断した: `modules_tune_sv_v03_68.n` (`InitMag`・`MatchingARCA`・`ChromCor` を持つ v03_68 用のモジュール) と、それを (v03_35 版 `modules_tune_sv_v03_35.n` の名前で) 使う `matching05` 系 (`matching05`・`_2`・`_3`・`_46`・`_46.5`・`_47`・`_47.5`・`_48`・`_48.5`・`_5`・`_6`・`_7`。セルのチューンを tunex, tuney の格子で振り、`MatchingARCA` → `ChromCor` → `DAsv` を行う)。→ `design/KEKLS_tune_survey_cell`。`modules_tune_sv_v03_68.n` を読み込む台本は元のフォルダに無い。`tune_survey*` (`tune_survey`・`_amp_dep_coef`(1〜3)・`_disp`・`_linear`・`_linear_from_zero`・`_momentum`・`_ncell`・`_write_QS`) は、別のリング (6 極 S1L・S2L…、チューン 15.x/5.x、`hsx_value.csv`) 用の `modules_tune_sv.n` を使い、v03_68 には適用できない |

## ラティス・モジュール

| 新しい名前 | 元 | 変更点 |
|---|---|---|
| `lattice/kek_ls_v03_68.sad` | `v03_68` (7127 B) | 先頭の 1 行 (`! 3.0 GEV RING, 20 CELLS`) をヘッダ (KEK-LS の説明、エネルギー・周長・セル数・ライン・チューン・εx) に替えた。改行を LF に。要素・ラインの定義は変えていない (先頭の 1 行以外が一致することを確認) |
| `modules/modules_general.n` | `modules_general.n` | `WriteGraphs[]` (`draw` コマンド) と使われない `GetA[]` を除いた。`DAsv[]` の中身は元のまま、追跡の周回数を変数 `DAturns` (既定 5000) にした |
| `modules/modules_amp_dep_ts.n` | `modules_amp_dep_ts.n` | コメントアウトされた旧い行を除いた。式・順序は元のまま。`ResetInit[]` を追加 (失敗した点の初期座標が次の点に残るのを防ぐ) |
| `modules/modules_cell_match.n` | `modules_tune_sv_v03_68.n` の `InitMag`・`MatchingARCA`・`MatchingARCA2`・`ChromCor` | v03_68 に合わせて直した (目標値、`MKSFA2` → `MKSFA1.2`、6 極 `SDA`・`SFA` → `SDH`・`SFH`、`BPI6/7` → `BPI6H/7H`)。`Chrom[]`・`GetKQ[]`・`PutKQ[]` を追加。元の `DAsv[]` は modules_general の `DAsv[]` とほぼ同じなので取り込まない |
| `modules/modules_errors.n` | (correct_cod.sad・print_*.sad の誤差を与える部分) | `SetErrors[]`・`ClearErrors[]`・`GetMagnetNames[]` にまとめた。位置ずれは `DELX`・`DELY` コマンドではなく成分ごとの `LINE["DX"/"DY"]` |
| (取り込まない) | `modules_matching_noID.n`・`modules_write_correctdata.n`・`modules_correctoptics.n`・`modules_correctoptics_pf.n`・`modules_tune_sv.n`・`modules_tune_sv_v03_35.n`・`drawmagvw3d.n` (元のフォルダに無い) | ID 入りラティス・別のリング・PF・v03_35 用、または描画。`correct_cod.sad` は `modules_correctoptics.n`・`modules_write_correctdata.n` を読み込まない |

## 設計・マッチング (`design/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 |
|---|---|---|---|
| `KEKLS_match_cell_tune` | `matching2_01` 前半 | 32 → 16 | `Get["modules_general.n"]`・`datafile`・描画・コメントの試行を除く。EMITX の FitFunction を除く。`free BPI6/7` は元の挙動 (効かない) を既定に残し `freebend` を足した。DAsv は含めない |
| `KEKLS_match_chrom` | `matching2_01` 後半 | (同上) | 前半から独立に、設計値の RINGC から始める。元の最後の `Print[...resx...]` は記号のまま出ていたので、ξ の差分表示に置き換えた |
| `KEKLS_match_injcell` | `matching2_03_INJ` | 16 | `WriteGraphs[]` を除く。50 回の繰り返し → 2 回。K1 を 1.2 倍にして始める手順を足した |
| `KEKLS_match_ring_geo` | `matching2_02_geo` (`matching20_geo`) | 64 | `RINGI` → `RING`。始めを 2.94 m にする手順を足した |
| `KEKLS_scale_mini` | `matching23_v03_35_mini` | 8 | v03_35 → v03_68 の名前、MKSFA1/2 → MKSFA1.1/.2、`free BPI6/7` → `BPI6H/7H`、factor の刻み 0.1、datas.txt への書き出しを除く |
| `KEKLS_scale_ncell` | `matching21_v03_35` (`matching21`・`matching22`・`matching21_v03_35_fix`) | 8 | v03_35 → v03_68 の名前、FitFunction (EMITX、`1^-9`) → `fit NX/NY`、j セルのリング `RINGJ..` を定義して εx を出す、j = 20,18,16,12,10 |
| `KEKLS_tune_survey_cell` | `matching05` 系 + `modules_tune_sv_v03_68.n` | 8 | モジュールは `modules/modules_cell_match.n` へ。格子 6×30 → 5×5。DA (5000 周) と `SAVE ALL`・`datas_46.txt` への書き出しを除く。RING → RINGC |
| `KEKLS_match_repeat` | (新規) | 8 | 再現性の診断 |

## パラメータ表示 (`display/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 |
|---|---|---|---|
| `KEKLS_ring_params` | `calc_emit.sad` 前半・`disp`・`disp_emit`・`chromaticity_plot.sad` | 64 (`disp` は 20 をコメント) | emit の表と ξ を RINGC・RING で出す。`TkWait`・`WriteGraphs`・`datas.txt` を除く |
| `KEKLS_optics_table` | `disp_optics`・`disp` (+ `draw_cad.sad` の入力) | 64 | 要素ごとの表と `disp geo $$$`。`draw_cad.sad` (`ring.vss`、`drawmagvw3d.n` が元のフォルダに無い) は取り込まない |
| `KEKLS_calc_H` | `calc_H.sad` (`calc_H_2.sad`・`calc_emit.sad` 後半) | 64 | γ を MOMENTUM から、Cq を m に、`DP` の上書きを `DJ` に |
| `KEKLS_chromaticity` | `chromaticity_plot.sad` (`amp_dep_ts_dp.sad`) | 64 | δ の刻み 0.5%、毎回の `DP0 = 0; cell calc` を省く、RINGC |
| `KEKLS_amp_dep_ts` | `amp_dep_ts.sad`・`amp_dep_ts_x.sad`・`amp_dep_ts_y.sad` | 64 | x・y の 2 本の表に、振幅の範囲を有効な範囲に、`ResetInit[]`、RINGC。後半の DAsv・旧い試行は除く |
| `KEKLS_sextupole_tune_scan` | `dp_sv_sx.sad` (`dp_sv_sx_dp_n0.03`・`_res`・`_v02_42` 系) | 8 | 6 極を SDH・SFH の 2 本に、基準を設計値に、5×5 |
| `KEKLS_scan_target` | `dp_sv.sad`・`dx_sv.sad`・`dx_sv_{ay,by,ex}.sad`・`dy_sv.sad`・`dy_sv_{ay,by,ex}.sad` | 8 | v02 の K1 の直接代入 → `InitMag[]`・`MatchingARCA[]`、AY・BY・EX を 3 点ずつ、RINGC |

## 誤差・補正 (`errors/`)

| 新しい名前 | 元 | 元の NPARA | 変更点 |
|---|---|---|---|
| `KEKLS_errors_cod` | `print_all_error.sad`・`print_dk_effect.sad`・`print_dt_error.sad`・`print_cod_magnitude.sad` | 8 | 誤差の種類ごとに表示、種 3 個。`DELX/DELY` → 成分ごとの `DX/DY` |
| `KEKLS_correct_cod` | `correct_cod.sad` (`_best`・`_worst`・`_noSX`・`_GCUT1`・`_v03_70`) | 64 | 種 2 個、補正 6 回、DA は任意。`MON ...`・`s!top;` を除く。`NOEXPAND`。ほかの `print_*.sad` の表示 (補正後の表など) の一部を統合 |
| `KEKLS_cod_response` | `write_cod_response.sad` | 64 | `BPM.*` → `PM.*`、書き出しを SAD の COD との照合に |

## 取り込まなかった元のスクリプト

全部の中身を読んだわけではない (約 290 ファイル。名前と、読んだ範囲から判断した)。

- `matching`・`matching01`〜`matching22`・`matching_v04_02` (ここで読んだのは 05・16_geo・20_geo・21・22・23 とその版違い): 設計の過程の試行。`matching05` 系は `design/KEKLS_tune_survey_cell`、`matching16_geo` は `STR`・`STRNINV` の別のライン (v03_68 に無い) を使い、`matching20_geo` は `matching2_02_geo` と同じ内容。
- `disp2`・`disp_arc_optics`・`disp_emit`・`disp_matching_tune`・`disp_matrix`: `modules_matching_noID.n`・`modules_write_correctdata.n` (ID 入りラティス用) を読む。
- `touscheklife.sad`・`touscheck*`: Touschek 寿命 (この例集の範囲外)。
- `da_xy`・`phase`・`test*`・`NG_script.sad`・`PF_v9_21_simple`・`SPring8`・`cod_pf.sad`: 別の用途・別のラティスの台本。
- `v01_*`・`v03_01`〜`v03_72` (v03_68 以外)・`v04_*`・`v05_*`・`v8_07_35nmrad`・`v03_35_*` (セル数や縮小率の変形)・`v03_44_*cells`・`~` つきのバックアップ・`#...#`: ほかの版のラティス。ベースは `v03_68` だけ。
- 版違い・断片の台本: `disp_b`・`disp_bunchlength`・`disp_eta`・`disp_NPARA_20/40`・`calc_H_2`・`dp_sv_2`・`dp_sv_R`・`dp_sv_v2_11`・`dx_sv_R`・`dx_sv_v2_11`・`dy_sv_R`・`dy_sv_v2_11`・`amp_dep_ts_{x,y,dp}`・`dp_sv_sx_res*`・`dp_sv_sx_v02_42`・`correct_cod_*` (設定違い)・`print_*.sad` (表示違い)。詳しい理由は各フォルダの README の「取り込んでいないもの」。
