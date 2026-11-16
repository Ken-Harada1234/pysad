# T5_injection: PF 入射軌道 (KEK SAD チュートリアル Tutorial5)

PF リング (2.5 GeV、周長 187.0355 m) の入射キッカー (KK1〜KK4) と入射ビームの軌道を、`TrackParticles` を 1 要素ずつ呼んで追う例。ラティスは `../lattice/v9_16_inj_noSF.sad` (`INJL` = 入射点から直線部の終わりまで 114 要素、27.7358 m)。放射・揺らぎは OFF。
各例は 4 ファイル (`NAME.sad`・`NAME.gs.txt`・`NAME.pysad.txt`・`NAME.md`)。走らせ方 (リポジトリの直下から): `py -3.12 example/tools/run_example.py example_PF_ASTRID/T5_injection/NAME.sad` (gs は WSL、pysad は `PYSAD_SAD_PACKAGES` つき)。凡例: ○ 一致 / △ 差あり (md に理由) / × gs 失敗 / × pysad 未対応。

| 例 | チュートリアルの節 | 内容 | gs と pysad |
|---|---|---|---|
| [calc_traj_01_PF](calc_traj_01_PF.md) | Tutorial5 計算部分の準備 | キッカーの入射バンプ (蓄積ビーム 11.2 mm@INJ)、入射ビーム 3 種の軌道。解析式 (Twiss) と 0.5〜1 % で一致 | ○ |
| [calc_traj_02_PF](calc_traj_02_PF.md) | Tutorial5 モジュール化 | 計算の関数化 (`CalcSBBump`・`DispCSInv` など)、Courant-Snyder 不変量、KK4 の掃引 | ○ |
| [calc_traj_03_PF](calc_traj_03_PF.md) | Tutorial5 GUI 入出力化 | KBMainFrame・Canvas・OpticsPlot・Entry・Button の窓 (教材の書き方) | ○ (計算の値も窓の Canvas の項目数 (114、再描画後 115) も一致。2026-11-17 に pysad の `Canvas$Widget=c1` と `OpticsPlot` の自前データを直した) |
| [calc_traj_03b_PF](calc_traj_03b_PF.md) | 同上 (ListPlot 版) | pysad でも動く `ListPlot` + `CanvasDrawer` の窓。pysad の窓の PNG | ○ (Canvas の項目 382 → 384 が一致。見た目は pysad のみ) |

図: `calc_traj_01_PF_traj.png` (gs の軌跡)、`calc_traj_03b_PF_window.png` (pysad が開いた窓)。
