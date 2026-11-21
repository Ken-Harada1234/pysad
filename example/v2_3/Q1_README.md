# example/v2_3: 担当 Q1 の例 (周波数マップ・局所運動量アクセプタンス、v2.3.0、2026-11-20)

解説は `docs/48_frequency_map_lma.md`、記録は `docs/history/v23_Q1.md`、試験は `pysad/tests/test_v23_Q1_fmap_lma.py` (21 件)、コードは `pysad/trksetup.py`・`freqmap.py`・`lma.py`・`q1cmd.py`。
実行はリポジトリ直下から `py -3.12 example/v2_3/Q1_xxx.py` (Windows では `PYTHONIOENCODING=utf-8` を付けると日本語の出力がファイルに書ける)。台本 `.sad` は `example/v2_3` で
`py -3.12 ../../run_sad.py ../astrid/lattice/ASTRID2_v01_01.sad Q1_script_commands.sad`。すべて v2 モード。gs の出力は WSL の gs (K64、2026-10-19 ビルド、NPARA = 1) で取った
(timeout 1800 s、終了後に `ps` で sad1.exe の残りがないことを確かめた)。使ったラティスは `kek-ls_v03_12.lat` (大きなリング A、3 GeV、603 m、誤差なし) と ASTRID2 (0.58 GeV、45.7 m、誤差なし)。PF は使っていない。

| ファイル | 内容 |
|---|---|
| `Q1_fmap_ring.py` / `Q1_fmap_astrid2.{out.txt,table.txt,png}` / `Q1_fmap_kekls.{out.txt,table.txt,png}` | **周波数マップ**: 1 つの位置 (長直線部) の 41 × 31 の (x, y) の格子、1024 周、RFSW・RAD OFF。表 (x, y, δ, 前半・後半のチューン、D、落ちた周)、図 (左: チューン平面を D で色づけ + 次数 4 以下の共鳴線、右: 初期条件平面、失われた点は黒い ×)、出力 (作動点に近い共鳴線の表 = 整数・半整数を含む、D の大きい点の最も近い共鳴線、共鳴線の近くの点の統計) |
| `Q1_lma_ring.py` / `Q1_lma_astrid2.{out.txt,table.txt,png}` / `Q1_lma_kekls.{out.txt,table.txt,png}` | **局所運動量アクセプタンス**: 周期 (ASTRID2 は 6 セル、大きなリング A は 20 セル) を確かめ、最初のセルの要素の境界で δ±(s) を求める。別のセル (3 つ) の同じ位置との比較、RF バケットの高さ・線形の見積もりとの比較、損失の判定の間隔 (check_every) の影響、Touschek 寿命 (局所の受容、一様な受容) |
| `Q1_gs_loss_compare.py` / `Q1_gs_loss_astrid2.{sad,gs.txt,pysad.txt,out.txt}` | **gs との損失判定の比較**: ASTRID2 (RF・放射つき) の 30 粒子を `TrackParticles[{loc, beam}, loc, 1, 1312]` で追跡。生き残りの旗が gs・pysad の台本・API で 30 / 30 一致。台本の中の二分法 (gs・pysad) の区間が同じで、LMA の区間がその内側。`.sad` が gs に渡した台本、`.gs.txt` が gs の出力、`.out.txt` が比較 |
| `Q1_gs_fmap_compare.py` / `Q1_gs_fmap_astrid2.{sad,gs.txt,pysad.txt,out.txt}` | **周波数マップの追跡の gs との比較**: gs で 1 周ずつ `TrackParticles` を回して 512 周の座標を印刷し、NAFF をかけて D を出す。座標の差 2.8e-7、D は 4 桁一致 |
| `Q1_script_commands.sad` / `.pysad.txt`、`Q1_script_fmap.{table.txt,png}`、`Q1_script_lma.{table.txt,png}` | **SAD の台本から**: `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` (pysad のみ。gs では走らない) |

結果の要約 (条件は `docs/48_frequency_map_lma.md` §5・§6):

- 周波数マップ: 作動点 (ASTRID2: νx = 5.18、νy = 2.14。大きなリング A: νx = 47.58、νy = 17.62) の整数・半整数の共鳴は格子のチューンから離れている (0.08 以上)。D の大きい点は
  νx − νy = 整数 (差の共鳴)、3 次 (3νy、νx + 2νy、3νx)、5 次の線の上に載る。次数 4 以下の線に近い点は D > −5 の割合が遠い点の約 2 倍。
- LMA: ASTRID2 (±10 mm、check_every = 1): δ+ = 0.83% (分散の大きい位置、線形の見積もり 0.865%)〜1.60% (δ_rf = 1.618%)。大きなリング A (±10 mm): 全位置が RF バケット (δ_rf = 2.99%) で決まり、
  δ+ = 3.01〜3.03%。別のセルの同じ位置との差は tol (1e-4) 以下 (最大 2.6e-5 と 4.9e-5)。
- 限界: 誤差・結合のないリングだけ。損失の判定は SAD の tapert (15 要素ごと) のままだと、アパーチャが限界のとき δ の限界を大きく見積もる (ASTRID2 の位置 17 で 1.05% 対 0.83%)。
