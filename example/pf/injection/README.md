# PF / injection: 入射軌道の計算と GUI

各例は、`NAME.sad` (台本。ラティスを先頭に含む)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明) の 4 点です。一致の凡例は ○ 一致、△ 差あり、× 未対応です。

| 例 | 内容 | gs と pysad |
|---|---|---|
| [`calc_traj_01_PF`](calc_traj_01_PF.md) | 入射軌道 (キッカーの入射バンプと入射ビームの軌道) | ○ |
| [`calc_traj_02_PF`](calc_traj_02_PF.md) | 入射軌道 (計算の関数化と Courant-Snyder 不変量) | ○ |
| [`calc_traj_03_PF`](calc_traj_03_PF.md) | 入射軌道 (Tk の窓 (KBMainFrame・Canvas・OpticsPlot) と入力欄・ボタン) | ○ (計算も窓の項目数 114・115 も一致。2026-11-17 に修正) |
| [`calc_traj_03b_PF`](calc_traj_03b_PF.md) | 入射軌道 (pysad でも動く ListPlot の窓) | ○ (Canvas の項目数 382 → 384 が一致) |
