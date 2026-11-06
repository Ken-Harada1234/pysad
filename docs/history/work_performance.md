# 作業ログ: 性能改善 (2026-10-16)

方針: **結果を変えずに速くする**。各段階で `QUICK=1 run_regression.py` FAIL 0、pytest (`--ignore` 3 本) 164 pass を確認。1e-10 以下の差・表示は無視。
変更ファイル: `pysad/ffscmd.py` (名前索引)、`pysad/lattice.py` (要素ごとのキャッシュ、twiss_table の直前結果キャッシュ、numba の Twiss 伝搬)、
`pysad/sadscript.py` (実数の四則の近道)、`pysad/emit6d.py` (tradke の D/B 更新を numba に。担当外だがホットスポットだったので最小限)、新規 `pysad/perfkern.py` (numba カーネル)。
`match6.py`, `match6d.py`, `dapert.py`, `track6n.py` は変更なし (下の「未対応」参照)。
CPU 制限 (ユーザー作業と並行) のため、numba の並列化 (prange) と multiprocessing は使っていない。長い実行は同時 1〜2 本。
作業用 (計測スクリプト・元コードの写し) は `<work-perf>/` (`<work-perf>/orig/pysad` = 作業前の pysad。`PYSAD_ORIG=1 python <work-perf>/prof_run.py 秒 出力.prof 格子 スクリプト` で中断付きプロファイル)。

## 1. プロファイル (手順 1) とホットスポット (作業前)

| 計測 | 結果 (作業前) | 主なホットスポット |
|---|---|---|
| (d) `Do[FFS["cell calc"],{20}]` (ESRF 3563 要素) | 0.49 s / 回 | `_transport_twiss` (Python 辞書で 1 要素 28 µs ×2 回)、`_one_turn_cod_matrix` (108 ms ×2)、`TwissRow` 生成 (27 ms ×2)。`cell calc` は `go` (残差計算で twiss) の後に `calc` でもう一度 twiss を計算する (計 2 回) |
| (d) `Do[Twiss["BX","Q*"],{200}]` | 534 ms / 回 | `positions()` のワイルドカード検索が O(n²): `_ordinal` が毎回先頭から数える (要素数 3563 → 540 万回の比較) |
| (d) `LINE["K1","Q*"]` | 520 ms / 回 | 同上 |
| (a) CodCorr 応答行列 (`FFS["cell calc"]` ×280) | 約 140 s 以上 | 上記 (毎回同じ状態の cell calc を計算し直す) |
| (a) CodCorr 補正 1 反復 (誤差入り) | 約 46 s (改善後も) | 誤差入りの閉軌道は要素ごとの Jacobian (BEND の tchge 数値微分、K1≠0 BEND の emit6d 直訳、薄レンズ) が Python で 1 要素 0.2〜1 ms。qcod の一周行列計算 10 回 ×1.5 s |
| (c) print_all_error | 1 実現 100 s 以上 | 誤差入りの `cell calc emit` (6D、tturne + tradke + BEND 本体) |
| (b) tune_survey 系 | MatchingWholeRing ×4 + DAsv | 4D の cell calc/マッチング (今回の改善対象) と DAsv (5000 周: SAD 785 s、pysad 1155 s) |

## 2. 対策と効果 (手順 2・3)

| # | 対策 | 内容 | 結果 |
|---|---|---|---|
| 1 | 要素名の索引 (`FFS._nidx`) | USE ごとに「名前→位置のリスト」「通し番号」「パターン→位置の結果」を作る。`positions/_family/_position/_ordinal/instances/elname` を O(1)〜O(一致数) に。ワイルドカードの正規表現はコンパイル結果を `lru_cache` | Twiss[...] 534 ms → **0.6 ms**、LINE[...] 520 ms → **0.8 ms** |
| 2 | 要素ごとの COD 写像キャッシュ (`lattice._el_cod_step_cached`) | 要素の属性 (`__dict__` の浅いコピー)・入口 COD・delta・RFSW が前回と同じなら前回の (M_el, 出口 COD) を返す (読み出し専用配列を共有。非スカラー属性を持つ要素はキャッシュしない) | 一周行列の計算 108 ms → 23 ms (設計軌道上・状態不変の場合)。cell calc 0.49 s → 0.25 s |
| 3 | numba の Twiss 伝搬 (`perfkern.twiss_chain`) | `_transport_twiss` (SAD qmat2twiss 直訳) を全要素に順に適用するカーネル。演算順は同じでビット単位で一致 (乱数行列 2000 個で最大差 0)。異常 (β≤0 / 非有限) は呼び出し側が従来どおり処理 | cell calc 0.25 s → **0.066 s** |
| 4 | `twiss_table` の直前結果キャッシュ (`Lattice._tt_*`) | 要素の属性・delta・閉軌道の初期値・安定判定の設定・MARK の値が同じなら、前回の TwissRow (複製)・副作用 (`_qstat`、`_cod_prev`、`_last_*` ほか) を復元して返す。`go` の後の `calc` の重複 twiss、状態不変の `cell calc` の繰り返しが対象 | cell calc **0.030 s** (状態不変) |
| 5 | sadscript: 実数どうしの Plus/Times/Subtract/Divide の近道 | `lmap2` を経由しない | 40000 回の式評価 2.34 s → 1.80 s |
| 6 | emit6d `tradke` の D・B 更新を numba に (`perfkern.tradke_db`) | numpy の小行列演算の連続を 1 カーネルに | 誤差入り `cell calc emit` 約 -10% (EMITX, COD は 1e-14 で一致) |

マイクロベンチ (ESRF 3563 要素、`<work-perf>/mb_calc.sad`):

| 項目 | 作業前 | 作業後 | 目標 |
|---|---|---|---|
| cell calc 1 回 (状態不変) | 0.49 s | 0.030 s | < 50 ms 達成 |
| cell calc 1 回 (QFA4 の K1 を変えた後。`go` の twiss と `calc` の twiss の 2 回分) | 0.42〜0.46 s | 0.145 s (twiss 1 回あたり約 95 ms: TwissRow 3563 個の生成 約 38 ms、一周行列 約 33 ms、memo の保存 約 20 ms) | < 50 ms 未達 |
| Twiss["BX","Q*"] 1 回 | 534 ms | 0.6 ms | < 1 ms 達成 |
| LINE["K1","Q*"] 1 回 | 520 ms | 0.8 ms | - |
| CodCorr 応答行列 (cell calc ×280) | 140 s 以上 | 7〜11 s | - |
| CodCorr 補正 1 反復 (誤差入り、2 calc + 280 回の K0 設定) | 約 46 s | 約 46 s (変化なし) | 未達 |

CodCorr の 1 実現 (40 反復) は **目標 (< 60 s) を達成できていない**: 誤差入りの閉軌道は毎反復ほぼ全要素の入口 COD が変わるのでキャッシュが効かず、
要素ごとの Jacobian (Python) が支配的 (BEND の `_tchge_jac` 数値微分 + `_bend_k1_cod_matrix` で一周 1.5 s ×10 回/反復)。
応答行列 (280 回の cell calc) の部分だけが 140 s → 約 10 s になった。

## 3. 未対応・見送り

- 誤差入りの要素ごとの Jacobian (BEND の tchge、K1≠0 BEND、薄レンズ) の numba 化。CodCorr (40 反復)、print_all_error、print_cod_magnitude の本命。解析式 (`_tchge_jac` の tbrote) が数値微分のまま。
- 6D (`calc6d`/`emit6d`) の要素ごとの計算 (tturne: 1 周数秒)。chromaticity_plot (6D→4D の off-momentum cell calc を 100 回)、disp_bunchlength (Emittance[] ×100) の本命。プロファイルが平らで (tradke、_tbendebody、_tbendecorr、piece、eye/copy が各 5〜10%) 小手先では 2 倍にならない。
- DAsv (`dapert`/`track6n`): 1 粒子 1 周 810 µs で SAD の 1.5 倍。FLUC ON は乱数列が粒子順に消費される (SAD の癖。docs/sad-tracking-particle-order-quirk) ので、粒子ごとの並列化は結果を変える。そもそも SAD も 5000 周の DA に 785 s かかる (tune_survey 系・disp_da・disp2・disp_test は SAD でも 300 s で打ち切り)。
- `match6.py`/`match6d.py` の累積行列の再利用: プロファイルでは Jacobian 計算より twiss/residual 計算が支配的で、要素キャッシュ (#2) の効果が既に入っている。

## 4. 検証

- `QUICK=1 run_regression.py`: 作業前 (索引だけ) / 作業後とも 50 ベンチ FAIL 0。pytest 164 passed (`--ignore` 3 本)。
- 出力の同一性: print_all_error の再実行出力は作業前の出力と、EMITY (1e-32 の雑音) 以外は同一。`mb_emiterr.sad` (誤差入り emit) は EMITX/EMITY/EMITZ/SIGE/COD が 1e-14 以内。
- キャッシュした配列は読み出し専用にしてあり、書き換えがあれば例外になる (回帰で例外なし)。

## 5. 打ち切り 10 本の再実行 (手順 4)


実行: `<work-perf>/rerun10.sh` (1 本ずつ順次、`timeout 600`、出力は `<work-recheck>/esrf.<name>.perf.out`、ログ `<work-recheck>/perf_log.txt`)。
比較は SAD 出力 (`esrf.<name>.sad.out`) の数値行と pysad の数値行を末尾から整列して 1e-6 (相対) で照合 (`<work-perf>/cmpnum.py` 相当)。

| スクリプト | 作業前 | 作業後 | SAD との数値照合 |
|---|---|---|---|
| chromaticity_plot | 600 s 打ち切り (148 行) | **完走 490 s** | 104 行 全て一致 |
| disp_bunchlength | 600 s 打ち切り (101 行) | **完走 545 s** | 100 行中 98 一致 (VOLT 0.1, 0.2 MV の BunchLength が 0.3% ずれる既知の差) |
| tune_survey_linear | 600 s 打ち切り (12471 行) | **完走 51 s** (38681 行) | 12790 行 全て一致 |
| tune_survey_disp | 600 s 打ち切り (1753 行) | **完走 542 s** (3613 行) | 722 行 全て一致 |
| tune_survey | 600 s 打ち切り (1457 行) | **完走 525 s** (2907 行) | 633 行 全て一致 |
| print_all_error | 600 s 打ち切り (60 行) | 打ち切り (61 行) | 先頭は作業前の出力と同一 (誤差入り 6D emit が支配的) |
| print_cod_magnitude | 600 s 打ち切り (78 行) | 打ち切り (129 行) | 同上 |
| disp2 | 600 s 打ち切り | 打ち切り (122 行) | DA (2000 周) は SAD も 296 s |
| disp_da | 600 s 打ち切り | 打ち切り (92 行) | DA (SAD も 300 s 内に終わらない) |
| disp_test | 600 s 打ち切り | 打ち切り (90 行) | DA 5000 周 (同上) |

打ち切り 10 本のうち **5 本が 600 s 内に完走** (目標の半分以上を達成)。tune_survey 3 本は MatchingWholeRing の繰り返し (4D の cell calc) が支配的で、
`Twiss`/`LINE` の索引と cell calc の高速化がそのまま効いた (いずれも点ごとに `daflag=1` で DA を飛ばす点が多い)。
残り 5 本は DAsv (5000 周の追跡: SAD でも 785 s) と誤差入りの 6D emit が支配的で、今回の対策では 600 s に収まらない。

QUICK 回帰は 50 ベンチ FAIL 0、pytest 164 passed (最終コードで再確認)。
