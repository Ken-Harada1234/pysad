# 周波数マップの例

解説は `docs/26_frequency_map_lma.md`、試験は `pysad/tests/script_tests/test_v23_Q1_fmap_lma.py` (21 件)、コードは `pysad/trksetup.py`・`freqmap.py`・`lma.py`・`q1cmd.py`。
実行はリポジトリ直下から `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_ring.py astrid2` など (Windows では `PYTHONIOENCODING=utf-8` を付けると日本語の出力がファイルに書ける)。台本 `.sad` は `example/14_nonlinear_dynamics/frequency_map` で
`py -3.12 ../../../run_sad.py ../../astrid/lattice/ASTRID2_v01_01.sad fn_FrequencyMap_MomentumAcceptance.sad`。局所運動量アクセプタンス (LMA) の例は `example/15_lma_touschek/lma/` (その `README.md`)。すべて v2 モード。gs の出力は WSL の gs (K64、2026-10-19 ビルド、NPARA = 1) で取った
(timeout 1800 s、終了後に `ps` で sad1.exe の残りがないことを確かめた)。使ったラティスは `kek-ls_v03_12.lat` (大きなリング A、3 GeV、603 m、誤差なし) と ASTRID2 (0.58 GeV、45.7 m、誤差なし)。PF は使っていない。

| ファイル | 内容 |
|---|---|
| `fmap_ring.py` / `fmap_astrid2.{out.txt,table.txt,png}` / `fmap_kekls.{out.txt,table.txt,png}` | **周波数マップ**: 1 つの位置 (長直線部) の 41 × 31 の (x, y) の格子、1024 周、RFSW・RAD OFF。表 (x, y, δ, 前半・後半のチューン、D、落ちた周)、図 (左: チューン平面を D で色づけ + 次数 4 以下の共鳴線、右: 初期条件平面、失われた点は黒い ×)、出力 (作動点に近い共鳴線の表 = 整数・半整数を含む、D の大きい点の最も近い共鳴線、共鳴線の近くの点の統計) |
| `gs_fmap_compare_astrid2.py` / `gs_fmap_astrid2.{sad,gs.txt,pysad.txt,out.txt}` | **周波数マップの追跡の gs との比較**: gs で 1 周ずつ `TrackParticles` を回して 512 周の座標を印刷し、NAFF をかけて D を出す。k64 モードで座標の差 2.8e-7、D は 6 桁一致 (v2 の既定の追跡は四極の積分器が違い、一致しない) |
| `fn_FrequencyMap_MomentumAcceptance.sad` / `.pysad.txt`、`fn_FrequencyMap_MomentumAcceptance_fmap.{table.txt,png}`、`fn_FrequencyMap_MomentumAcceptance_lma.{table.txt,png}` | **SAD の台本から**: `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` (pysad のみ。gs では走らない) |

結果の要約 (条件は `docs/26_frequency_map_lma.md` §5・§6):

- 周波数マップ: 作動点 (ASTRID2: νx = 5.18、νy = 2.14。大きなリング A: νx = 47.58、νy = 17.62) の整数・半整数の共鳴は格子のチューンから離れている (0.08 以上)。D の大きい点は
  νx − νy = 整数 (差の共鳴)、3 次 (3νy、νx + 2νy、3νx)、5 次の線の上に載る。次数 4 以下の線に近い点は D > −5 の割合が遠い点の 1.9〜2.7 倍 (ASTRID2 2.7 倍、大きなリング A 1.9 倍)。
