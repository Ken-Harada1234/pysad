# PF / optics: 線形光学・分散・軌道応答・描画

各例は、`NAME.sad` (台本。ラティスを先頭に含む)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明) の 4 点です。一致の凡例は ○ 一致、△ 差あり、× 未対応です。

| 例 | 内容 | gs と pysad |
|---|---|---|
| [`FFS_cell_calc_emit_PF`](FFS_cell_calc_emit_PF.md) | 実行、`cell calc emit`、`SAVE`、周長・チューン・エミッタンス | ○ |
| [`FFS_disp_type_PF`](FFS_disp_type_PF.md) | `disp`・`type` | ○ |
| [`FFS_draw_OpticsPlot_PF`](FFS_draw_OpticsPlot_PF.md) | `draw`、`OpticsPlot` (図の数値と PNG) | ○ (draw は gs 無出力・pysad 未対応) |
| [`FFS_draw_topdrawer_PF`](FFS_draw_topdrawer_PF.md) | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |
| [`calc_chromaticity_PF`](calc_chromaticity_PF.md) | `PolynomialFit` の 6 次、自然色収差、四極・六極の和の式 | ○ |
| [`calc_dispersion_PF`](calc_dispersion_PF.md) | 分散の差分と $\alpha_c$ の確認 | ○ |
| [`calc_respx_PF`](calc_respx_PF.md) | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| [`fn_PolynomialFit_PF`](fn_PolynomialFit_PF.md) | チューン対運動量、`DP0`、`PolynomialFit` (SAD の係数の癖あり) | ○ (k64 モード。v2 の既定は x の平均が 0 のデータで真の係数を返し gs の壊れた係数と違う: B18) |
| [`fn_Twiss_LINE_PF`](fn_Twiss_LINE_PF.md) | `Twiss[]`・`LINE[]` | ○ |
