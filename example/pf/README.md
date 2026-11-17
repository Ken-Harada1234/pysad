# PF の例集 (2.5 GeV)

PF のラティスを使った、SAD の基本的な使い方と計算の例です。KEK の SAD チュートリアルの話題を、PF の実際のラティスで作り直しました。各例は `.sad`、`.gs.txt`、`.pysad.txt`、`.md` の 4 点です。式と数値の照合は各 `.md` にあります。

## 索引

| フォルダ | 例 | 内容 | gs と pysad |
|---|---|---|---|
| [`lattice/`](lattice/README.md) | ラティス | PF のラティス | |
| [`optics/`](optics/README.md) | [`FFS_cell_calc_emit_PF`](optics/FFS_cell_calc_emit_PF.md) | 実行、`cell calc emit`、`SAVE`、周長・チューン・エミッタンス | ○ |
| [`optics/`](optics/README.md) | [`FFS_disp_type_PF`](optics/FFS_disp_type_PF.md) | `disp`・`type` | ○ |
| [`optics/`](optics/README.md) | [`FFS_draw_OpticsPlot_PF`](optics/FFS_draw_OpticsPlot_PF.md) | `draw`、`OpticsPlot` (図の数値と PNG) | ○ (draw は gs 無出力・pysad 未対応) |
| [`optics/`](optics/README.md) | [`FFS_draw_topdrawer_PF`](optics/FFS_draw_topdrawer_PF.md) | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |
| [`optics/`](optics/README.md) | [`calc_chromaticity_PF`](optics/calc_chromaticity_PF.md) | `PolynomialFit` の 6 次、自然色収差、四極・六極の和の式 | ○ |
| [`optics/`](optics/README.md) | [`calc_dispersion_PF`](optics/calc_dispersion_PF.md) | 分散の差分と $\alpha_c$ の確認 | ○ |
| [`optics/`](optics/README.md) | [`calc_respx_PF`](optics/calc_respx_PF.md) | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| [`optics/`](optics/README.md) | [`fn_PolynomialFit_PF`](optics/fn_PolynomialFit_PF.md) | チューン対運動量、`DP0`、`PolynomialFit` (SAD の係数の癖あり) | ○ |
| [`optics/`](optics/README.md) | [`fn_Twiss_LINE_PF`](optics/fn_Twiss_LINE_PF.md) | `Twiss[]`・`LINE[]` | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_chrom_PF`](matching/FFS_FIT_chrom_PF.md) | 色収差補正の六極フィット (設計値が戻る) | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_injbump_PF`](matching/FFS_FIT_injbump_PF.md) | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| [`matching/`](matching/README.md) | [`FFS_FIT_modules_PF`](matching/FFS_FIT_modules_PF.md) | `ChromCorrect[]`・`MatchingNcellTune[]` のモジュール化 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`FFS_amp_dep_ts_PF`](tune_radiation/FFS_amp_dep_ts_PF.md) | 振幅依存チューンシフト (12 点平均) と図 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`FFS_mom_dep_ts_PF`](tune_radiation/FFS_mom_dep_ts_PF.md) | 運動量依存チューンシフト (六極あり・なし) と図 | ○ |
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
