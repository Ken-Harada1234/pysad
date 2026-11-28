# PF の例集 (2.5 GeV)

PF のラティスを使った、SAD の基本的な使い方と計算の例です。KEK の SAD チュートリアルの話題を、PF の実際のラティスで作り直しました。各例は `.sad`、`.gs.txt`、`.pysad.txt`、`.md` の 4 点です。式と数値の照合は各 `.md` にあります。

## 索引

| フォルダ | 例 | 内容 | gs と pysad |
|---|---|---|---|
| [`PF_in_notebook.ipynb`](PF_in_notebook.ipynb) | notebook 版 | ラティスを冒頭で読み込み、セルの計算・放射積分・分散の確認・振幅依存・運動量依存チューンを順に流す (出力つき) | ○ (`PYSAD_COMPAT=k64` で gs と一致) |
| [`lattice/`](lattice/README.md) | ラティス | PF のラティス (v9_16 の 2 本と、履歴ラティス 4 本。説明は [`PF_history.md`](lattice/PF_history.md)) | |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF`](optics/FFS_cell_calc_emit_PF.md) | 実行、`cell calc emit`、`SAVE`、周長・チューン・エミッタンス | ○ |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF_1981_460nmrad`](optics/FFS_cell_calc_emit_PF_1981_460nmrad.md) | 履歴ラティス (1981 年・建設時) の `cell calc emit`、直線部・標準セルの Twiss、クロマティシティ | ○ (k64) |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF_1986_128nmrad`](optics/FFS_cell_calc_emit_PF_1986_128nmrad.md) | 履歴ラティス (1986 年・最初の高輝度化) の `cell calc emit`、直線部・標準セルの Twiss、クロマティシティ | ○ (k64) |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF_1997_036nmrad`](optics/FFS_cell_calc_emit_PF_1997_036nmrad.md) | 履歴ラティス (1997 年・弧部改造) の `cell calc emit`、直線部・標準セルの Twiss、クロマティシティ | ○ (k64) |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF_2005_036nmrad`](optics/FFS_cell_calc_emit_PF_2005_036nmrad.md) | 履歴ラティス (2005 年・直線部増強) の `cell calc emit`、直線部・標準セルの Twiss、クロマティシティ | ○ (k64) |
| [`optics/`](optics/README.md) | [`FFS_disp_type_PF`](optics/FFS_disp_type_PF.md) | `disp`・`type` | ○ |
| [`optics/`](optics/README.md) | [`FFS_draw_OpticsPlot_PF`](optics/FFS_draw_OpticsPlot_PF.md) | `draw`、`OpticsPlot` (図の数値と PNG) | ○ (draw は gs 無出力・pysad 未対応) |
| [`optics/`](optics/README.md) | [`FFS_draw_topdrawer_PF`](optics/FFS_draw_topdrawer_PF.md) | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |
| [`optics/`](optics/README.md) | [`calc_chromaticity_PF`](optics/calc_chromaticity_PF.md) | `PolynomialFit` の 6 次、自然色収差、四極・6 極の和の式 | ○ |
| [`optics/`](optics/README.md) | [`calc_dispersion_PF`](optics/calc_dispersion_PF.md) | 分散の差分と $\alpha_c$ の確認 | ○ |
| [`optics/`](optics/README.md) | [`calc_respx_PF`](optics/calc_respx_PF.md) | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| [`optics/`](optics/README.md) | [`fn_PolynomialFit_PF`](optics/fn_PolynomialFit_PF.md) | チューン対運動量、`DP0`、`PolynomialFit` (SAD の係数の癖あり) | ○ (k64 モード。v2 の既定は x の平均が 0 のデータで真の係数を返し gs の壊れた係数と違う: B18) |
| [`optics/`](optics/README.md) | [`fn_Twiss_LINE_PF`](optics/fn_Twiss_LINE_PF.md) | `Twiss[]`・`LINE[]` | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_chrom_PF`](matching/FFS_FIT_chrom_PF.md) | 色収差補正の6 極フィット (設計値が戻る) | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_injbump_PF`](matching/FFS_FIT_injbump_PF.md) | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_modules_PF`](matching/FFS_FIT_modules_PF.md) | `ChromCorrect[]`・`MatchingNcellTune[]` のモジュール化 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`FFS_amp_dep_ts_PF`](tune_radiation/FFS_amp_dep_ts_PF.md) | 振幅依存チューンシフト (12 点平均) と図 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`FFS_mom_dep_ts_PF`](tune_radiation/FFS_mom_dep_ts_PF.md) | 運動量依存チューンシフト (6 極あり・なし) と図 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`calc_RadInt_PF`](tune_radiation/calc_RadInt_PF.md) | 放射積分、εx・U0・σE・Jx の式と数値の照合 | ○ (εx は SAD と 0.24%) |
| [`tracking/`](tracking/README.md) | [`fn_DynamicApertureSurvey_DAsv_PF`](tracking/fn_DynamicApertureSurvey_DAsv_PF.md) | DA サーベイ (x・y 方向、運動量 9 点、2000 周) | △ (y の 1 行だけ gs 11 / pysad 13) |
| [`tracking/`](tracking/README.md) | [`fn_DynamicApertureSurvey_PF`](tracking/fn_DynamicApertureSurvey_PF.md) | DA サーベイ (1000 周) | ○ (表示記号に 1〜2 文字の差) |
| [`tracking/`](tracking/README.md) | [`fn_DynamicApertureSurvey_da_xy_PF`](tracking/fn_DynamicApertureSurvey_da_xy_PF.md) | X-Y 面の DA (24 方向、1000 周) と図 | ○ |
| [`tracking/`](tracking/README.md) | [`fn_TrackParticles_PF`](tracking/fn_TrackParticles_PF.md) | 四極 1 個の転送行列 | ○ |
| [`tracking/`](tracking/README.md) | [`phase_PF`](tracking/phase_PF.md) | 位相空間 (60 周) と図 | ○ |
| [`injection/`](injection/README.md) | [`calc_traj_01_PF`](injection/calc_traj_01_PF.md) | 入射軌道 (キッカーの入射バンプと入射ビームの軌道) | ○ |
| [`injection/`](injection/README.md) | [`calc_traj_02_PF`](injection/calc_traj_02_PF.md) | 入射軌道 (計算の関数化と Courant-Snyder 不変量) | ○ |
| [`injection/`](injection/README.md) | [`calc_traj_03_PF`](injection/calc_traj_03_PF.md) | 入射軌道 (Tk の窓 (KBMainFrame・Canvas・OpticsPlot) と入力欄・ボタン) | ○ (計算も窓の項目数 114・115 も一致。2026-11-17 に修正) |
| [`injection/`](injection/README.md) | [`calc_traj_03b_PF`](injection/calc_traj_03b_PF.md) | 入射軌道 (pysad でも動く ListPlot の窓) | ○ (Canvas の項目数 382 → 384 が一致) |

## notebook 版

[`PF_in_notebook.ipynb`](PF_in_notebook.ipynb) は、このリングのラティスを冒頭で読み込み、このフォルダの例のコマンドを `%%sad` と `%%sad --continue` で順に流す notebook です (出力つきで保存してあります。`example/tools/make_ring_notebooks.py` が、元の `.sad` からコマンドを取り出して作ります。動かし方は [`../17_notebook/README.md`](../17_notebook/README.md))。

- 保存した出力は pysad の既定 (`v2`) の実行結果です。`PYSAD_COMPAT=k64` で実行すると、εx・σδ・ダンピング分配数・U0 は gs と全桁で一致しました (確認した量: εx 3.45380e-8 m、σδ 7.33656e-4、Jx 1.0313、U0 0.3991454 MV。出力の数値の行のうち、gs の出力に同じ値がある割合は 196 行中 188 行。残りは、書式の違いと 1e-12 以下の数値雑音と、gs にない出力 (pysad 独自の解析式など))。`v2` との差は、`v2` が K64 のバグを直した動作であることによります (docs/28)。
- 確認した範囲: Windows・Python 3.12・IPython 9.19.1・ipykernel 7.4.0・nbclient 0.11.0 で、本物のカーネルを起動して実行しました。実行時間は約 20 秒です。
