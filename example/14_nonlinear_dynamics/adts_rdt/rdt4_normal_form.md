# 4 次の共鳴駆動項 (6 極の 2 次を含む)・偏向磁石の生成関数の原理・数値の正規形・周波数マップの予想

**この文書で分かること**: 4 次の共鳴駆動項 (6 極の 2 次を含む)、偏向磁石の 3 次の生成関数、数値の正規形による ADTS、周波数マップの共鳴線の予想への 2 次の効果の入れ方と、追跡との一致度。
**前提**: `rdt_tracking.md` (3 次までの RDT) と ADTS の用語 (`rdt_adts_basic.md`)。
**読む順**: 結果の要約 → 限界・未確認。

台本ごとの式・解析の内訳・追跡 + FFT (NAFF) の数値・比較・コマンドは `ADTS_summary.md` §12 (rdt4_tracking.py・bend_generators.py・rdt4_gs_map.py・normal_form_numeric.py) と `../frequency_map/FMAP_summary.md` (fmap_rdt4.py)。解説は `docs/25` (§8.1 原理、§10 4 次の項)・`docs/24` (§3.4・§5.5.1)、試験は `pysad/tests/test_v23_R5_rdt4.py` (14 件)、コードは `pysad/rdt2nd.py` (4 次の項)・`pysad/nform.py` (数値の正規形)・`pysad/rdttrack.py` (追跡の 3 次のテイラー展開)・`pysad/resdrive.py` (`resonance_terms` の 4 次の名前)・`pysad/fmapres.py` (`order2=True`)。
ラティス: FODO (`example/00_lattice/FODO_ring.sad`、SF・SD の K2L = 4, −6、3 GeV、105.6 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、ASTRID2 MPW (`ASTRID2_v01_02_id.sad` の `RINGMPW`)、KEK-LS (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、603 m)。
条件 (すべて): pysad k64 (gs と同じ式) の追跡 (RF 空洞止め、放射なし、結合なし、閉軌道のまわり、δ=0)、リング先頭が位相の基準、3 次のテイラー展開の振幅 amp3 = FODO・ASTRID2 1e-4 √m、KEK-LS 1e-5 √m (正規化座標)、2 次の係数の振幅 1e-5 √m、偏向磁石の分割 16。gs は `NORFSW;`、timeout 1800 s、1 本ずつ (終了後 `sad1.exe` の残りなし)。実行は標準出力のリダイレクトに `PYTHONUTF8=1`、CPU 1 コア。**スクリプトは出力を `example/14_nonlinear_dynamics/adts_rdt/` へ上書きする**。

| ファイル | 内容 |
|---|---|
| `rdt4_tracking.py` | `table <ring>`: 解析 `rdt4` 対 追跡の $f_4$ (6 極 ×0・0.5・1・2)、`adts <ring>`: 対角 + 正規形の補正 = `adts()`、直接の追跡 (NAFF) との比較、`plot` |
| `rdt4_table_<ring>.pysad.{txt,json}`・`rdt4_adts_<ring>.pysad.{txt,json}`・`rdt4_table.png` | 上の出力 (fodo・fodo_rect・fodo_asym・astrid・astrid_mpw・kekls) |
| `bend_generators.py`・`.out.txt` | 偏向磁石の 3 次の生成関数: 追跡から取り出した係数 対 解析式 (B1)(B2)、E1/E2≠0 の端 (B3、解析式なし)、四極 1 個の 4 次の項 (運動学の台形則・フリンジ) |
| `rdt4_gs_map.py`・`rdt4_gs_map_{fodo,astrid}.gs.txt` | 検証用: gs の 1 周の追跡 360 点から取り出した $h_4$ 対 pysad・解析 |
| `normal_form_numeric.py`・`normal_form_numeric.out.txt` | 数値の正規形で $\partial\nu/\partial J$ (J の 1 次) |
| `../frequency_map/fmap_rdt4.py`・`fmap_rdt4_<ring>.{out.txt,png}` (周波数マップ側) | 周波数マップの共鳴線の予想: 1 次の RDT だけ 対 2 次 (4 次の h) も |

結果の要約:
- 解析 `rdt4` と追跡の $f_4$ の最大相対差 (35 項中 |h| の大きい 15〜19 項): FODO 1.3e-4、ASTRID2 1.2e-3、ASTRID2 MPW 1.8e-4、KEK-LS 5.3e-6 (6 極 ×1)。薄い 6 極の列の FODO は 2e-6。残りは K64 の厚い 6 極の積分器の分。6 極 ×0 は 3e-8 (MPW 8.7e-3)。λ の 2 次の依存: 解析は厳密、追跡は 1e-5〜4e-5。
- ADTS: 写像の対角 + 正規形の補正 = `adts()` (1e-15)。追跡 + NAFF との差は 1e-4 以下 (4 次の項を入れない解析 (`fourth='r3'` 相当) では −2.4e-3〜+2.6e-3) (KEK-LS 8e-6 以内、数値の正規形は 1e-7)。「6 極に依らない 4 次の定数」は、K64 の四極の積分器 (台形則、`kinematic_h4(quad_rule='k64')`) と偏向磁石の 4 次の項で解消。
- 周波数マップ: ASTRID2 は 1 次だけでは領域が格子に入らないが、2 次 (2νx−2νy = νx−νy = 3、4νy) で D > −5 の点の 48% が領域に入る (面積 31%、対照 30 ± 13%、p = 0.2): 有意には改善しない。原因の分析 (νx−νy ≈ 3 の帯が広い、次数 5 以上の共鳴) は `../frequency_map/FMAP_summary.md`。
- 数値の正規形 (J の 1 次): 直接の追跡と 1e-7〜2e-6 (MPW 3e-4)。J²・大振幅の改善は未実装 (設計は docs/25 §10.6)。

限界・未確認: docs/25 §10.8 (非対角の調和の正規形の補正、E1/E2≠0 の端の解析式、δ の 4 次の項、K0≠0 の偏向磁石の内部、J² の正規形、KEK-LS の大振幅)。
再現: `PYTHONUTF8=1 py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_tracking.py table fodo` (1 s)、`... table kekls` (30 s)、`... adts fodo` (数十秒)、`py -3.12 example/14_nonlinear_dynamics/adts_rdt/bend_generators.py`、`py -3.12 example/14_nonlinear_dynamics/adts_rdt/normal_form_numeric.py`、`py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_rdt4.py astrid2`、`py -3.12 example/14_nonlinear_dynamics/adts_rdt/rdt4_gs_map.py fodo` (gs、約 1 分)。
