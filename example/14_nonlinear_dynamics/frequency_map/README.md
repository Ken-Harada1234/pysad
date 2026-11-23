# 周波数マップ (frequency map) の例

周波数マップ (チューン拡散 D = log10 √(Δνx² + Δνy²) の地図) の作り方、共鳴駆動項 (RDT)・振幅依存チューンシフト (ADTS) から予想した共鳴の領域との比較、gs (SAD K64) との比較、結果を信頼できる範囲の確認の例です。
実行はリポジトリ直下から `py -3.12 example/14_nonlinear_dynamics/frequency_map/<スクリプト> ...` (標準出力をファイルに落とすときは `PYTHONUTF8=1`)。すべて既定の v2 モードで (gs との比較のスクリプトも `compat.set_mode('v2')` のまま。gs と一致するのは `'k64'` に変えて走らせたとき。下の「使った条件・限界・未確認」)、周波数マップは RFSW・RAD OFF、LOSSAMPL = 0.1 m、1024 周 (前半・後半 512 周)、NAFF の窓は Hann。
gs を走らせるスクリプト (`gs_fmap_compare_*`・`fmap_check_gs_*`) は WSL の gs と環境変数 `PYSAD_SAD_PACKAGES` (K64 版 SAD の Packages のフォルダ) が要る。`--nogs` を付けると保存済みの gs の出力を使う。

説明の文書:

| 文書 | 内容 |
|---|---|
| [`fmap_basic.md`](fmap_basic.md) | 基本の例 (周波数マップ、gs との比較、SAD の台本から) |
| [`fmap_validation.md`](fmap_validation.md) | 検証 (RDT・ADTS との比較、δ・k64・誤差・MPW・格子) |
| [`fmap_reliability.md`](fmap_reliability.md) | 信頼性 (誤差の統計、MPW、RFSW・FLUC、K1 の設定、周回数、gs・誤差つき) |
| [`fmap_resonance_lines.md`](fmap_resonance_lines.md) | 共鳴線の予想 (3 次 + 4 次 + 非共鳴の補正) と判別力の評価 |
| [`FMAP_summary.md`](FMAP_summary.md) | 台本ごとのコマンド・式・結果・追跡の数値・比較・未確認のまとめ (全数値) |

## 例の一覧

名前の `fmap_*` は周波数マップの本体と検証、`gs_*` は gs との比較、`fmap_check_*` は信頼性の確認、`fn_*` は SAD の台本から呼ぶ関数の例です。各例は台本 (`.py`・`.sad`)・出力 (`.out.txt`・`.table.txt`・`.gs.txt`・`.pysad.txt`)・図 (`.png`) の 1 組です。

| 例 (台本) | 出力・図 | 何の例か・結果の要点 |
|---|---|---|
| `fmap_ring.py` | `fmap_astrid2.*`・`fmap_kekls.*` | ASTRID2・KEK-LS (kek-ls_v03_12.lat) の x–y の周波数マップ (41 × 31)。D の大きい点が次数 4 以下の共鳴線 (差の共鳴 νx − νy = 整数など) に載る |
| `fn_FrequencyMap_MomentumAcceptance.sad` | `.pysad.txt`・`_fmap.*`・`_lma.*` | SAD の台本から `FrequencyMap[]`・`MomentumAcceptance[]`・`TouschekLifetimeLMA[]` を呼ぶ (pysad のみ。gs では走らない)。周波数マップは 88 点、失われた点 0 |
| `gs_fmap_compare_astrid2.py`・`gs_fmap_astrid2.sad` | `gs_fmap_astrid2.*` | ASTRID2 の追跡 (512 周) を gs の `TrackParticles` と比較。k64 モードで座標の最大の差 2.8e-7 m、D は 6 桁一致 (v2 の既定では座標の差 3e-2 m) |
| `gs_fmap_compare_kekls.py`・`gs_fmap_kekls.sad` | `gs_fmap_kekls.*` | KEK-LS の追跡を gs と比較 (512 周、8 点)。k64 モードで規則的な点は座標の差 1.2e-11 m (v2 の既定では 6.1e-3 m) |
| `fmap_rdt.py` | `fmap_rdt_{fodo,astrid2,kekls}.*` | D の大きい領域を RDT・ADTS から予想した共鳴の領域と比較。FODO は 3νx の予想に D > −5 の点がすべて入る。ASTRID2 は 1 次の RDT では予想できない |
| `fmap_variants.py` | `fmap_{delta,k64,errors,mpw,grid,lambda}_*.*`、`fmap_variants_*.out.txt` | モード別の検証: x–δ・y–δ の地図、k64 との差 (チューンの差の中央値 9e-5・3e-5、四極の積分器の次数の違い)、位置ずれ・ロールの誤差つき、挿入光源 (MPW) つき、NAFF の窓・格子、6 極の強さ λ 倍 |
| `fmap_grid_summary.py` | `fmap_grid_summary.out.txt` | 81 × 61 の格子から部分格子を取り出し、D の分布の格子の細かさへの依存を見る |
| `fmap_rdt4.py` | `fmap_rdt4_{fodo,astrid2,kekls}.*` | 4 次の RDT も入れた共鳴の領域の予想との比較 |
| `fmap_predict.py` | `fmap_predict_{fodo,astrid2,kekls,mpw}.*` | 共鳴線の予想 (3 次 + 4 次 + 非共鳴の補正) の判別力の評価。FODO の 3νx のような孤立した低次の共鳴は予想できるが、ASTRID2・KEK-LS は改善しない |
| `fmap_check_errors_scan.py`・`fmap_check_errors_summary.py`・`fmap_check_coupling.py`・`fmap_check_kekls_orbit.py` | `fmap_check_data/{errors_*,coupling}.jsonl`、`fmap_check_errors_summary_*.out.txt`、`fmap_check_kekls_orbit.out.txt` | 誤差つき周波数マップの統計 (複数の種・強さ・要素の種類)。D > −5 の割合は誤差なしの 7〜12 倍、主因は差の共鳴を駆動する結合 |
| `fmap_check_mpw.py` | `fmap_check_data/mpw.jsonl`、`fmap_check_mpw.out.txt` | MPW で多くの点が失われる原因の切り分け。4νy = 9 の共鳴ではなく、MPW の内部の非線形 |
| `fmap_check_rfsw.py` | `fmap_check_data/rfsw.jsonl` | RFSW・RAD・FLUC を ON にした周波数マップ。D が拡散の指標でなくなる |
| `fmap_check_gs_rfsw.py`・`fmap_check_gs_rfsw.sad` | `fmap_check_gs_rfsw.*` | RFSW ON の追跡の gs との比較 |
| `fmap_check_k1opts.py` | `fmap_check_data/k1opts.jsonl`、`fmap_check_k1opts.out.txt` | K1 の設定 (v2 と k64 で動作が違う設定) ごとの周波数マップ |
| `fmap_check_nturn.py`・`fmap_check_grid.py` | `fmap_check_data/nturn_*.npz`、`fmap_check_nturn_*.out.txt`、`fmap_check_grid.out.txt` | D の周回数 (512〜4096)・NAFF の窓・格子の細かさへの依存 |
| `fmap_check_gs_errors.py`・`fmap_check_gs_errors.sad` | `fmap_check_gs_errors.*` | 誤差つき ASTRID2 の追跡の gs との比較 (pysad の誤差の値を gs の台本に明示して書く) |

`fmap_check_data/` は誤差・MPW・結合・K1・RFSW などの統計の保存データ (`.jsonl` は 1 行 = 1 つの周波数マップ、`.npz` は D_N の配列) です。`fmap_check_*.py` はこの保存データに追記します。

## 使った条件・限界・未確認

- ラティスは `example/astrid/lattice/ASTRID2_v01_01.sad` (0.58 GeV、周長 45.7 m)、`pysad/tests/kek-ls_v03_12.lat` (3 GeV、周長 603 m、20 セル)、`example/00_lattice/FODO_ring.sad` (24 セル、3 GeV、105.6 m)、`ASTRID2_v01_02_id.sad` の `RINGMPW`。磁場上限の制約はなく、既存のラティスをそのまま使っている。
- D の絶対値の閾値 (D > −5 など) は周回数・窓に依存する (周回数を 2 倍にすると D が約 1 下がる)。別の条件の結果とは比べられない。
- RDT・ADTS の予想は、1 次の RDT が駆動する共鳴だけの場合 ASTRID2 では当たらない。4 次の項・非共鳴の補正を入れても ASTRID2・KEK-LS は改善しなかった。
- 追跡の gs との比較は NORFSW・NORAD の座標だけで、pysad を k64 モード (四極 2 次の積分器。gs と同じ) にしたときに一致する。v2 の既定 (四極 4 次) は積分器が違うので一致しない。保存済みの `gs_fmap_kekls.*` の 8 点は、積分器を変える前の表から選んだ点 (`.pysad.txt` は v2 の既定で取り直した出力)。K64 に NAFF が無いので、gs の座標に pysad の NAFF をかけている。カオス的な点 (D > −3 など) は初期値の 1e-13 m の摂動の増幅と同じ大きさで食い違いうる。
- 補正つき・PF リング・MPW の RDT の予想は未確認。誤差つきの統計は ASTRID2 (補正なし) と KEK-LS (補正なし)。
- 保存済みの `fmap_astrid2.table.txt`・`fmap_kekls.table.txt` は 4 次の合成の四極の積分器で作ったもの。`gs_fmap_*.pysad.txt` など一部の保存済みの出力・台本の標準出力の見出しには、作業中の呼び名 (S1・R5 など) が残っている。
- 保存済みの `.pysad.txt` のうち `gs_fmap_astrid2`・`fmap_check_gs_rfsw`・`fn_FrequencyMap_MomentumAcceptance` は現行の v2 で取り直した (保存済みの出力との座標の差は相対 1e-5 程度 (絶対 3e-9 m)。gs との比較の出力 (`.out.txt`) は取り直していない)。他の台本の出力は取り直していない。
- 局所運動量アクセプタンス (LMA) の例は [`../../15_lma_touschek/`](../../15_lma_touschek/README.md)。ADTS・RDT の例は [`../adts_rdt/`](../adts_rdt/)。
