# example/14_nonlinear_dynamics (README・normal_form・primer) から取り除いた履歴の記述

旧パスは `example/v2_4/` (normal_form) と `example/v2_3/` (primer)。接頭辞の対応: A1_ → ampnf_vs_tracking、A2_ → k0_bend_check、A5_ → nu_of_J_fmap_grid、A7_ → ampnf_kekls、F1_ → rotated_bend、G_ → k0_*、H1_ → coupled_h4、I1_ → nu_of_J_orders、I2_ → island_vs_tracking。

## normal_form の台本の冒頭コメント
- 先頭の版・項目の符号を取り除いた: 「v2.4.0 A1:」「v2.4.0 A2:」「v2.4.0 A5:」「v2.4.0 A7:」「v2.4.0 F1:」「v2.4.0 G:」「v2.4.0 H1:」「v2.4.1 I1:」「v2.4.1 I2:」。
- `k0_ampnf.py`: 「(出力は G_k0_ampnf.out.txt、記録は docs/history/v24_G.md)」。`k0_ampnf.out.txt` の 2 行目: 「記録は docs/history/v24_G.md。」
- `nu_of_J_fmap_grid.py`: 「周波数マップの予想 (S3) が …」「保存済みの周波数マップ (Q1_fmap_*.table.txt …)」「# v2.4.0 A9: …」「# v2.4.0 A13: …」。使い方の行の「PYSAD_COMPAT=k64」(スクリプトは中で v2 に切り替える)。
- 周波数マップの表の付け替え: 保存した出力は、旧い 2 次の積分器で作った表 (旧 `Q1_fmap_astrid2.table.txt`・`Q1_fmap_kekls.table.txt`) を読んで作った。現行の読み込み先は `../frequency_map/fmap_astrid2.table.txt`・`fmap_kekls.table.txt` (4 次の合成の表)。出力の取り直しは README の「未確認・限界」に書いた範囲。
- `k0_bend_check.py` の使い方の行: 「(variant: base, ang0, ang, edge, rot)」(base・rot という名前の変種は無く、現在の名前に直した)。
- 出力ファイルの先頭の「# example/v2_4/G_k0_ampnf.py の出力」などの旧パスを、新しいパスに直した。

## primer の台本
- `primer_fmap_lma.py` の「(Q1_fmap_ring.py・Q1_lma_ring.py は 41×31 と 25 位置)」: 新しい名前 (frequency_map/fmap_ring.py・lma/lma_ring.py) に直した。
- `primer_breakdown.py` は、保存済みの `R3_scan_kekls.pysad.txt` (現 `../adts_rdt/adts_bend_scan_kekls.pysad.txt`) を読む。
