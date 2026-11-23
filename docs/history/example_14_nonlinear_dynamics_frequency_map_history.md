# example/14_nonlinear_dynamics/frequency_map の整理の履歴

`example/v2_3`・`v2_4` の作業フォルダから機能別のフォルダへ移したときに、md の本文・台本の冒頭コメントから取り除いた作業符号・経過の記述を、元の文脈のまま残す。

旧パス→新パスの対応 (主なもの): v2_3 の `Q1_fmap_*_v241`→`fmap_*`、`Q1_lma_*`→`lma_*`、`Q1_gs_*`→`gs_*`、`S1_*`→`fmap_*`、`S2_*`→`fmap_check_*` (`S2_data/`→`fmap_check_data/`)、`S3_*`→`fmap_predict*`、`R5_fmap_rdt2`→`fmap_rdt4`、`L2_*`→`lma_*` (`L2_elegant_isr`→`touschek_elegant/elegant_isr`)、`T1_*`→接頭辞なし、`T2_*`→`touschek_formula_*`、`U2_*`→`bend_k1_def_*`、v2_4 の `A3_*`→`touschek_bunch_*`、`A4_*`→`lma_support*`、`H2_*`→`lma_corrected_*`。旧い 2 次の積分器で作った `Q1_fmap_{astrid2,kekls}.{out.txt,table.txt,png}` (6 ファイル) は削除し、4 次の合成の `_v241` 版を `fmap_{astrid2,kekls}.*` にした。

## md から取り除いた文

### fmap_basic.md

元: # 周波数マップ・局所運動量アクセプタンスの例

後: # 周波数マップの例

### fmap_basic.md

元: 解説は `docs/48_frequency_map_lma.md`、記録は `docs/history/v23_Q1.md`、試験は

後: 解説は `docs/48_frequency_map_lma.md`、試験は

### fmap_basic.md

元: 実行はリポジトリ直下から `py -3.12 example/15_lma_touschek/touschek_elegant/Q1_xxx.py` (Windows

後: 実行はリポジトリ直下から `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_ring.py astrid2` など (Windows

### fmap_basic.md

元: 台本 `.sad` は `example/v2_3` で
> `py -3.12 ../../run_sad.py ../astrid/lattice/ASTRID2_v01_01.sad fn_FrequencyMap_MomentumAcceptance.sad`。

後: 台本 `.sad` は `example/14_nonlinear_dynamics/frequency_map` で
`py -3.12 ../../../run_sad.py ../../astrid/lattice/ASTRID2_v01_01.sad fn_FrequencyMap_MomentumAcceptance.sad`。局所運動量アクセプタンス (LMA) の例は `example/15_lma_touschek/lma/` (その `README.md`)。

### fmap_validation.md

元: 、記録は `docs/history/v23_S1.md`、試験は

後: 、試験は

### fmap_validation.md

元: 実行はリポジトリ直下から `py -3.12 example/15_lma_touschek/touschek_elegant/S1_xxx.py ...`

後: 実行はリポジトリ直下から `py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_variants.py ...` など

### fmap_validation.md

元: スクリプトは出力を `example/v2_3/` へ書く (保存済みの出力を残したいときは複製したフォルダで走らせる)。

後: スクリプトは出力を `example/14_nonlinear_dynamics/frequency_map/` へ書く (保存済みの出力を残したいときは複製したフォルダで走らせる)。

### fmap_resonance_lines.md

元: 、記録は `docs/history/v23_S3.md`、試験は

後: 、試験は

### fmap_reliability.md

元: 結果の全体は `docs/history/v23_S2.md`、解説は

後: 解説は

### fmap_reliability.md

元: `errors_astrid2.jsonl`・`fmap_check_errors_kekls.jsonl`

後: `errors_astrid2.jsonl`・`errors_kekls.jsonl`

### FMAP_summary.md

元: `example/v2_3/` の周波数マップ (frequency map、FMA) の台本

後: `example/14_nonlinear_dynamics/frequency_map/` の周波数マップ (frequency map、FMA) の台本

### FMAP_summary.md

元: 最後の §9 に簡潔にまとめました。ADTS・共鳴駆動項 (RDT) の同じ形のまとめは [`ADTS_summary.md`](ADTS_summary.md)。

後: 最後の §11 に簡潔にまとめました (LMA の例は `../../15_lma_touschek/lma/`)。ADTS・共鳴駆動項 (RDT) の同じ形のまとめは [`ADTS_summary.md`](../adts_rdt/ADTS_summary.md)。

### FMAP_summary.md

元: 解説は [`docs/48_frequency_map_lma.md`](../../docs/48_frequency_map_lma.md)、記録は [`docs/history/v23_Q1.md`](../../docs/history/v23_Q1.md) と [`docs/history/v23_S1.md`](../../docs/history/v23_S1.md)、
> 索引は

後: 解説は [`docs/48_frequency_map_lma.md`](../../../docs/48_frequency_map_lma.md)、
索引は

### FMAP_summary.md

元: **数値は `example/v2_3/` の出力ファイルから転記したもの。**

後: **数値はこのフォルダの出力ファイルから転記したもの。**

### FMAP_summary.md

元: 保存済みの出力 (`Q1_*`) を上書きしないよう

後: 保存済みの出力 (`fmap_astrid2.*` など) を上書きしないよう

### FMAP_summary.md

元: (スクリプトは出力を `example/v2_3/` へ書く)

後: (スクリプトは出力をこのフォルダへ書く)

### FMAP_summary.md

元: (cwd は `example/v2_3`)

後: (cwd は `example/14_nonlinear_dynamics/frequency_map`)

### FMAP_summary.md

元: cd example/v2_3
> py -3.12 ../../run_sad.py ../astrid/lattice/ASTRID2_v01_01.sad fn_FrequencyMap_MomentumAcceptance.sad

後: cd example/14_nonlinear_dynamics/frequency_map
py -3.12 ../../../run_sad.py ../../astrid/lattice/ASTRID2_v01_01.sad fn_FrequencyMap_MomentumAcceptance.sad

### FMAP_summary.md

元: 偽の ADTS を作る。docs/history/v23_R1.md)。

後: 偽の ADTS を作る)。

### FMAP_summary.md

元: MPW の損失の記録 (`docs/history/v23_S2.md`) と整合: νy

後: §14.2 の MPW の損失の切り分けと整合: νy

### FMAP_summary.md

元: 、記録は [`docs/history/v23_S2.md`](../../docs/history/v23_S2.md)、解説は [`docs/48`](../../docs/48_frequency_map_lma.md)

後: 、解説は [`docs/48`](../../../docs/48_frequency_map_lma.md)

### FMAP_summary.md

元: 図: [`lma_astrid2.png`](lma_astrid2.png)、[`lma_kekls.png`](lma_kekls.png)。

後: 図: [`lma_astrid2.png`](../../15_lma_touschek/lma/lma_astrid2.png)、[`lma_kekls.png`](../../15_lma_touschek/lma/lma_kekls.png)。

### FMAP_summary.md

元: ![](lma_astrid2.png)

後: ![](../../15_lma_touschek/lma/lma_astrid2.png)

## 台本の冒頭コメントから取り除いた・書き直した行

### fmap_check_coupling.py

元: `"""v2.3.0 S2 (1): 結合 (線形のスキュー四極) だけで D > −5 の割合が上がるかの確かめ。`

後: `"""結合 (線形のスキュー四極) だけで D > −5 の割合が上がるかの確かめ。`

### fmap_check_errors_scan.py

元: `"""v2.3.0 S2 (1): 誤差つき周波数マップの統計 (複数の乱数の種 × 複数の強さ × 誤差の種類の切り分け)。`

後: `"""誤差つき周波数マップの統計 (複数の乱数の種 × 複数の強さ × 誤差の種類の切り分け)。`

元: `KEK-LS (pysad/tests/kek-ls_v03_12.lat) は位置 3、x = 0.2〜6 mm・y = 0.2〜8 mm。格子は 21 × 16 (S1 の full の 1 つおき、336 点)。`

後: `KEK-LS (pysad/tests/kek-ls_v03_12.lat) は位置 3、x = 0.2〜6 mm・y = 0.2〜8 mm。格子は 21 × 16 (full の格子 41 × 31 の 1 つおき、336 点)。`

### fmap_check_errors_summary.py

元: `"""v2.3.0 S2 (1): fmap_check_errors_scan.py の結果 (fmap_check_data/errors_<ring>.jsonl) を (誤差の種類, 強さ) ごとに集計する。`

後: `"""fmap_check_errors_scan.py の結果 (fmap_check_data/errors_<ring>.jsonl) を (誤差の種類, 強さ) ごとに集計する。`

### fmap_check_grid.py

元: `"""v2.3.0 S2 (5): 格子の細かさへの依存の再確認。Q1 の 41 × 31 の表 (ASTRID2・KEK-LS) から 1 つおき・2 つおき・3 つおきの部分格子を取り出し、D の分布 (百分位・D > −5・D > −4・失われた点) を比べる`

後: `"""格子の細かさへの依存の再確認。fmap_<ring>.table.txt の 41 × 31 の表 (ASTRID2・KEK-LS) から 1 つおき・2 つおき・3 つおきの部分格子を取り出し、D の分布 (百分位・D > −5・D > −4・失われた点) を比べる`

### fmap_check_gs_errors.py

元: `"""v2.3.0 S2 (7): 誤差つき ASTRID2 の追跡を gs (SAD K64) と比べる。`

後: `"""誤差つき ASTRID2 の追跡を gs (SAD K64) と比べる。`

### fmap_check_gs_rfsw.py

元: `"""v2.3.0 S2 (3): RFSW ON (RAD OFF) の追跡 (ASTRID2、誤差なし、位置 2) を gs (SAD K64) と比べる。FFS に入ると RFSW は ON (R1) なので、gs の既定の台本はこの条件。
6 点 (Q1 と同じ x, y) を 512 周、1 周ずつ TrackParticles (RFSW・NORAD、LOSSAMPL = 0.1 m) で回し、(x, px, y, py, z, δ) を周ごとに印刷、pysad の同じ追跡
(TrackSetup.track(rfsw=True)) と台本版と比べる。RFSW OFF の同じ追跡 (Q1) との差 (z・δ の大きさ) も出す。`

後: `"""RFSW ON (RAD OFF) の追跡 (ASTRID2、誤差なし、位置 2) を gs (SAD K64) と比べる。FFS に入ると RFSW は ON なので、gs の既定の台本はこの条件。
6 点 (fmap_ring.py と同じ x, y) を 512 周、1 周ずつ TrackParticles (RFSW・NORAD、LOSSAMPL = 0.1 m) で回し、(x, px, y, py, z, δ) を周ごとに印刷、pysad の同じ追跡
(TrackSetup.track(rfsw=True)) と台本版と比べる。RFSW OFF の同じ追跡 (fmap_ring.py の既定) との差 (z・δ の大きさ) も出す。`

### fmap_check_k1opts.py

元: `"""v2.3.0 S2 (4): v2 と k64 で動作が違う設定での周波数マップ (ASTRID2、21 × 16、1024 周)。`

後: `"""v2 と k64 で動作が違う設定での周波数マップ (ASTRID2、21 × 16、1024 周)。`

元: `elegant・MAD と同じ定義、U2)、face_lin (両方)、k64 (互換モード k64)。設定は TrackSetup を作る前に切り替える (setup が edge・def の値を持ち、frequency_map が合わせる)。`

後: `elegant・MAD と同じ定義)、face_lin (両方)、k64 (互換モード k64)。設定は TrackSetup を作る前に切り替える (setup が edge・def の値を持ち、frequency_map が合わせる)。`

### fmap_check_kekls_orbit.py

元: `"""v2.3.0 S2 (1): KEK-LS (kek-ls_v03_12.lat) の誤差つき (補正なし) の線形の性質: 強さ ×0.03〜×3、種 1〜5 の閉軌道・線形の安定性・結合・フィードダウンの見積もり。`

後: `"""KEK-LS (kek-ls_v03_12.lat) の誤差つき (補正なし) の線形の性質: 強さ ×0.03〜×3、種 1〜5 の閉軌道・線形の安定性・結合・フィードダウンの見積もり。`

### fmap_check_mpw.py

元: `"""v2.3.0 S2 (2): 挿入光源 (ASTRID2_v01_02_id.sad の RINGMPW) で多くの点が失われる原因の切り分け。`

後: `"""挿入光源 (ASTRID2_v01_02_id.sad の RINGMPW) で多くの点が失われる原因の切り分け。`

### fmap_check_nturn.py

元: `"""v2.3.0 S2 (5): D の周回数 (512〜4096)・NAFF の窓・格子の細かさへの依存。`

後: `"""D の周回数 (512〜4096)・NAFF の窓・格子の細かさへの依存。`

元: `窓 (なし = 0、Hann の次数 1・2・3) も同じ追跡から作る。格子は S1 の full の格子 (41 × 31、Q1 の表) から 1 つおき・2 つおきなどの部分格子で、D の分布を比べる。`

後: `窓 (なし = 0、Hann の次数 1・2・3) も同じ追跡から作る。格子は full の格子 (41 × 31、fmap_<ring>.table.txt) から 1 つおき・2 つおきなどの部分格子で、D の分布を比べる。`

### fmap_check_rfsw.py

元: `"""v2.3.0 S2 (3): RFSW ON・RAD ON・FLUC ON の周波数マップ (ASTRID2、位置 2)。既定 (RFSW・RAD OFF) と比べる。`

後: `"""RFSW ON・RAD ON・FLUC ON の周波数マップ (ASTRID2、位置 2)。既定 (RFSW・RAD OFF) と比べる。`

### fmap_grid_summary.py

元: `"""v2.3.0 S1 (f): 格子の細かさへの依存。fmap_variants.py grid astrid2 81 61 の表 (fmap_grid_astrid2_81x61.table.txt、4941 点、1024 周) から、
1 つおき (step 2 = 41 × 31、Q1 と同じ格子)、3 つおき (step 4 = 21 × 16)、7 つおき (step 8 = 11 × 8) の部分格子を取り出し、D の分布の統計を比べる。`

後: `"""格子の細かさへの依存。fmap_variants.py grid astrid2 81 61 の表 (fmap_grid_astrid2_81x61.table.txt、4941 点、1024 周) から、
1 つおき (step 2 = 41 × 31、fmap_ring.py と同じ格子)、3 つおき (step 4 = 21 × 16)、7 つおき (step 8 = 11 × 8) の部分格子を取り出し、D の分布の統計を比べる。`

### fmap_predict.py

元: `予想: S1 (3 次の h)、R5 (+ 4 次の h の写像の係数)、S3 (+ 非共鳴の 3 次の摂動の補正 F4' = rdt4_normal_form)。ν(J) は解析の線形の ADTS (adts(fourth='r3', optics='calc'))。`

後: `予想: 3 次の h だけ (出力の表記 S1。fmap_rdt.py と同じ)、+ 4 次の h の写像の係数 (R5。fmap_rdt4.py と同じ)、+ 非共鳴の 3 次の摂動の補正 F4' = rdt4_normal_form (S3)。ν(J) は解析の線形の ADTS (adts(fourth='r3', optics='calc'))。`

### fmap_rdt.py

元: `"""v2.3.0 S1 (a): 周波数マップの D の大きい領域を、共鳴駆動項 (RDT) と ADTS から予想した共鳴の領域と定量的に比べる。`

後: `"""周波数マップの D の大きい領域を、共鳴駆動項 (RDT) と ADTS から予想した共鳴の領域と定量的に比べる。`

元: `  astrid2・kekls は Q1 の周波数マップの表 (fmap_<ring>.table.txt、1024 周) を読む (--recompute で追跡し直して S1_ の名前で保存)。`

後: `  astrid2・kekls は fmap_ring.py の周波数マップの表 (fmap_<ring>.table.txt、1024 周) を読む (--recompute で追跡し直して fmap_rdt_<ring>.table.txt に保存)。`

### fmap_rdt4.py

元: `"""周波数マップの D の大きい領域を、1 次の RDT だけ (S1) と 2 次の RDT (4 次の h_jklm、6 極の 2 次など) も入れた共鳴の領域の予想で比べる。`

後: `"""周波数マップの D の大きい領域を、1 次の RDT だけ (fmap_rdt.py。出力の表記 S1) と 2 次の RDT (4 次の h_jklm、6 極の 2 次など) も入れた共鳴の領域の予想で比べる。`

元: `周波数マップは S1・Q1 の保存済みの表を読む (astrid2・kekls: fmap_<ring>.table.txt、fodo: fmap_rdt_fodo.table.txt。1024 周、v2、RFSW・RAD OFF、LOSSAMPL 0.1 m)。
予想 (fmapres.ResonancePredictor): V = 2|F|、F = Σ h (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) (調和 a = (j−k, l−m) ごと)。S1 は 3 次の h だけ、R5 は 4 次の h (map の f4) も同じ調和に足す。
  振り子の式 (位相つき)・包絡の定義、ν(J) = 線形の ADTS (または追跡のチューン)、対照 = 領域を格子の上で巡回ずらし (400 回)、は S1 と同じ。`

後: `周波数マップは保存済みの表を読む (astrid2・kekls: fmap_<ring>.table.txt、fodo: fmap_rdt_fodo.table.txt。1024 周、v2、RFSW・RAD OFF、LOSSAMPL 0.1 m)。
予想 (fmapres.ResonancePredictor): V = 2|F|、F = Σ h (2Jx)^((j+k)/2) (2Jy)^((l+m)/2) (調和 a = (j−k, l−m) ごと)。1 次だけ (S1) は 3 次の h だけ、2 次も (R5) は 4 次の h (map の f4) も同じ調和に足す。
  振り子の式 (位相つき)・包絡の定義、ν(J) = 線形の ADTS (または追跡のチューン)、対照 = 領域を格子の上で巡回ずらし (400 回)、は fmap_rdt.py と同じ。`

### fmap_ring.py

元: `出力名の付け足し (例 _v241) を付けると fmap_<名前><付け足し>.* に書く (既存の表は上書きしない)。v2.3.0 の表 (付け足しなし) は旧い四極の積分器 (2 次)、
_v241 の表は v2.4.1 の 4 次の合成 (B226・B243) で作った。`

後: `出力名の付け足し (例 _test) を付けると fmap_<名前><付け足し>.* に書く (保存済みの表を上書きしない)。保存済みの表は 4 次の合成の四極の積分器で作った。`

### fmap_variants.py

元: `  k64 <ring>            (c) 互換モード k64 (PYSAD_COMPAT=k64 と同じ) で同じ周波数マップを走らせ、v2 (Q1 の表) と点ごとに比べる`

後: `  k64 <ring>            (c) 互換モード k64 (PYSAD_COMPAT=k64 と同じ) で同じ周波数マップを走らせ、v2 (fmap_<ring>.table.txt) と点ごとに比べる`

元: `  格子は full = 41 × 31 (Q1 と同じ)、coarse = 21 × 16 (full の 1 つおき)。周回数 1024 (前半・後半 512 周)。`

後: `  格子は full = 41 × 31 (fmap_ring.py と同じ)、coarse = 21 × 16 (full の 1 つおき)。周回数 1024 (前半・後半 512 周)。`

### fn_FrequencyMap_MomentumAcceptance.sad

元: `! v2.3.0 Q1: SAD の台本から使う FrequencyMap[]・MomentumAcceptance[]・TouschekLifetimeLMA[] (K64 にない pysad の追加機能。gs では走らない)。`

後: `! SAD の台本から使う FrequencyMap[]・MomentumAcceptance[]・TouschekLifetimeLMA[] (K64 にない pysad の追加機能。gs では走らない)。`

### gs_fmap_compare_kekls.py

元: `"""v2.3.0 S1 (g): 周波数マップの追跡の部分を gs (SAD K64) と比べる (gs_fmap_compare_astrid2.py の一般化。KEK-LS の kek-ls_v03_12.lat でも)。`

後: `"""周波数マップの追跡の部分を gs (SAD K64) と比べる (gs_fmap_compare_astrid2.py の一般化。KEK-LS の kek-ls_v03_12.lat でも)。`

元: `点は Q1 の表 (fmap_<ring>.table.txt) から D の値で選ぶ (生き残る 6 点 + 失われる 2 点)。512 周 (前半・後半 256 周)。`

後: `点は fmap_ring.py の表 (fmap_<ring>.table.txt) から D の値で選ぶ (生き残る 6 点 + 失われる 2 点)。512 周 (前半・後半 256 周)。`

### gs_fmap_kekls.sad

元: ``

後: `! KEK-LS is an ESRF-type HMBA designed with the help of Pantaleo Raimondi and Simone Liuzzo (ESRF).`

元: `beam = {{0.00107,0.0023749999999999999,0.0042599999999999999,0.00063500000000000004,0.001505,0.0036800000000000001,0.0046950000000000004,0.00020000000000000001}, {0,0,0,0,0,0,0,0}, {0.0017600000000000001,0.0041000000000000003,0.00046000000000000001,0.0051399999999999996,0.0074799999999999997,0.00332,0.0080000000000000002,0.00696}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {1,1,1,1,1,1,1,1}};`

後: `beam = {{0.00107,0.0029550000000000002,0.00165,0.0013600000000000001,0.00048999999999999998,0.00048999999999999998,0.0046950000000000004,0.0060000000000000001}, {0,0,0,0,0,0,0,0}, {0.0028,0.00097999999999999997,0.0061799999999999997,0.0074799999999999997,0.0061799999999999997,0.0064400000000000004,0.0080000000000000002,0.00046000000000000001}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {1,1,1,1,1,1,1,1}};`

## FMAP_summary.md の数値の更新 (四極の積分器を 4 次にしたあと)

ASTRID2・KEK-LS の x–y の周波数マップ (§1〜§4、§10、§12、§13) の数値を、現行の表 (`fmap_astrid2.table.txt`・`fmap_kekls.table.txt`) と、そこから `fmap_rdt.py`・`fmap_rdt4.py`・`fmap_predict.py` で出した値に直した。旧い値は四極が 2 次の積分器の表からの値で、旧い表からの再現は ASTRID2 で確認した (`fmap_grid_astrid2_81x61.table.txt` の 1 つおき)。主な旧 → 新 (ASTRID2 / KEK-LS):
- 失われた点 0 / 460 → 0 / 465 (落ちた周の最大 936 → 889)。D の範囲 −11.45〜−1.70 / −10.80〜−1.09 → −11.37〜−1.44 / −11.00〜−1.18。D の分布 566、613、71、20、1 / 211、358、180、52、10 (811 点) → 580、582、93、14、2 / 215、356、175、46、14 (806 点)。
- D > −5・D > −4 の点 48、21 / 124、62 → 49、16 / 112、60。近い / 遠い 379 / 892・507 / 304 → 379 / 892・501 / 305。D > −5 の割合 0.061 / 0.028・0.195 / 0.082 → 0.069 / 0.026・0.170 / 0.089。
- §4.5 (B): KEK-LS の線形の ADTS の cov 0.266 (対照 0.378、p = 0.83) → 0.295 (0.379、0.75)、追跡のチューン cov 0.734 (面積 0.438、対照 0.282 ± 0.131) → 0.696 (0.432、0.275 ± 0.133)。尾根と予想の中心の差 +1.29 ± 0.96・−1.35 ± 1.66 mm → +1.15 ± 1.00・−1.45 ± 1.82 mm。
- §12.5: ASTRID2 の 4 次の h も入れた予想 cov 0.479 / area 0.312 / 対照 0.315 ± 0.131 / p 0.20 → 0.551 / 0.312 / 0.317 ± 0.154 / 0.11。KEK-LS 追跡のチューンの 1 次の h だけ 0.758、4 次も 0.927 → 0.714、0.902。
- §13.4 の判別力の表 (AUC など) と §13.6・§13.7 の表も同様に更新した (ASTRID2 の AUC 0.675 → 0.682、KEK-LS 0.594 → 0.599)。
- §6 (k64): 「チューンの差 5e-11、失われた点の食い違い 0」は、2 次の積分器の表と k64 の表の比較。現行の v2 (4 次) との差はチューンの差の中央値 9e-5・3e-5、食い違い 0・2 点。
- §3・§10・§14.3・§14.7 の gs との比較 (座標 2.80e-7 m、KEK-LS の規則的な点 1.2e-11 m など) は k64 モードの結果として残した (k64 で再現できた)。v2 の既定の追跡では一致しない (ASTRID2 の 6 点で座標の最大の差 3.3e-2 m、KEK-LS の 8 点で 9.8e-3 m)。
