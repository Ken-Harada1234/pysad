# 作業ログ: 性能改善 (2026-10-16)

最終更新: 2026-11-16

方針: **結果を変えずに速くする**。各段階で `QUICK=1 run_regression.py` FAIL 0、pytest (`--ignore` 3 本) 164 pass を確認。1e-10 以下の差・表示は無視。
変更ファイル: `pysad/ffscmd.py` (名前索引)、`pysad/lattice.py` (要素ごとのキャッシュ、twiss_table の直前結果キャッシュ、numba の Twiss 伝搬)、
`pysad/sadscript.py` (実数の四則の近道)、`pysad/emit6d.py` (tradke の D/B 更新を numba に。担当外だがホットスポットだったので最小限)、新規 `pysad/perfkern.py` (numba カーネル)。
`match6.py`, `match6d.py`, `dapert.py`, `track6n.py` は変更なし (下の「未対応」参照)。
CPU 制限 (他の作業と並行) のため、numba の並列化 (prange) と multiprocessing は使っていない。長い実行は同時 1〜2 本。
作業用 (計測スクリプト・元コードの写し) は作業フォルダ (リポジトリに含まれない)。

## 1. プロファイル (手順 1) とホットスポット (作業前)

| 計測 | 結果 (作業前) | 主なホットスポット |
|---|---|---|
| (d) `Do[FFS["cell calc"],{20}]` (3563 要素のリング) | 0.49 s / 回 | `_transport_twiss` (Python 辞書で 1 要素 28 µs ×2 回)、`_one_turn_cod_matrix` (108 ms ×2)、`TwissRow` 生成 (27 ms ×2)。`cell calc` は `go` (残差計算で twiss) の後に `calc` でもう一度 twiss を計算する (計 2 回) |
| (d) `Do[Twiss["BX","Q*"],{200}]` | 534 ms / 回 | `positions()` のワイルドカード検索が O(n²): `_ordinal` が毎回先頭から数える (要素数 3563 → 540 万回の比較) |
| (d) `LINE["K1","Q*"]` | 520 ms / 回 | 同上 |
| `FFS["cell calc"]` を 280 回 (K1 や K0 の設定を挟む) | 約 140 s 以上 | 上記 (毎回同じ状態の cell calc を計算し直す) |
| 誤差入りのリングの閉軌道補正 1 反復 (2 calc + 280 回の K0 設定) | 約 46 s (改善後も) | 誤差入りの閉軌道は要素ごとの Jacobian (BEND の tchge 数値微分、K1≠0 BEND の emit6d 直訳、薄レンズ) が Python で 1 要素 0.2〜1 ms。qcod の一周行列計算 10 回 ×1.5 s |
| 誤差入りのリングの `cell calc emit` を繰り返す計算 | 1 実現 100 s 以上 | 誤差入りの `cell calc emit` (6D、tturne + tradke + BEND 本体) |
| 4D のマッチングの繰り返し + DynamicApertureSurvey | 4D の cell calc/マッチング (今回の改善対象) と DynamicApertureSurvey (5000 周: SAD 785 s、pysad 1155 s) |

## 2. 対策と効果 (手順 2・3)

| # | 対策 | 内容 | 結果 |
|---|---|---|---|
| 1 | 要素名の索引 (`FFS._nidx`) | USE ごとに「名前→位置のリスト」「通し番号」「パターン→位置の結果」を作る。`positions/_family/_position/_ordinal/instances/elname` を O(1)〜O(一致数) に。ワイルドカードの正規表現はコンパイル結果を `lru_cache` | Twiss[...] 534 ms → **0.6 ms**、LINE[...] 520 ms → **0.8 ms** |
| 2 | 要素ごとの COD 写像キャッシュ (`lattice._el_cod_step_cached`) | 要素の属性 (`__dict__` の浅いコピー)・入口 COD・delta・RFSW が前回と同じなら前回の (M_el, 出口 COD) を返す (読み出し専用配列を共有。非スカラー属性を持つ要素はキャッシュしない) | 一周行列の計算 108 ms → 23 ms (設計軌道上・状態不変の場合)。cell calc 0.49 s → 0.25 s |
| 3 | numba の Twiss 伝搬 (`perfkern.twiss_chain`) | `_transport_twiss` (SAD qmat2twiss 直訳) を全要素に順に適用するカーネル。演算順は同じでビット単位で一致 (乱数行列 2000 個で最大差 0)。異常 (β≤0 / 非有限) は呼び出し側が従来どおり処理 | cell calc 0.25 s → **0.066 s** |
| 4 | `twiss_table` の直前結果キャッシュ (`Lattice._tt_*`) | 要素の属性・delta・閉軌道の初期値・安定判定の設定・MARK の値が同じなら、前回の TwissRow (複製)・副作用 (`_qstat`、`_cod_prev`、`_last_*` ほか) を復元して返す。`go` の後の `calc` の重複 twiss、状態不変の `cell calc` の繰り返しが対象 | cell calc **0.030 s** (状態不変) |
| 5 | sadscript: 実数どうしの Plus/Times/Subtract/Divide の近道 | `lmap2` を経由しない | 40000 回の式評価 2.34 s → 1.80 s |
| 6 | emit6d `tradke` の D・B 更新を numba に (`perfkern.tradke_db`) | numpy の小行列演算の連続を 1 カーネルに | 誤差入り `cell calc emit` 約 -10% (EMITX, COD は 1e-14 で一致) |

マイクロベンチ (3563 要素のリング、`mb_calc.sad`):

| 項目 | 作業前 | 作業後 | 目標 |
|---|---|---|---|
| cell calc 1 回 (状態不変) | 0.49 s | 0.030 s | < 50 ms 達成 |
| cell calc 1 回 (四極 1 台の K1 を変えた後。`go` の twiss と `calc` の twiss の 2 回分) | 0.42〜0.46 s | 0.145 s (twiss 1 回あたり約 95 ms: TwissRow 3563 個の生成 約 38 ms、一周行列 約 33 ms、memo の保存 約 20 ms) | < 50 ms 未達 |
| Twiss["BX","Q*"] 1 回 | 534 ms | 0.6 ms | < 1 ms 達成 |
| LINE["K1","Q*"] 1 回 | 520 ms | 0.8 ms | - |
| cell calc ×280 (K1 や K0 の設定を挟む) | 140 s 以上 | 7〜11 s | - |
| 誤差入りの閉軌道補正 1 反復 (2 calc + 280 回の K0 設定) | 約 46 s | 約 46 s (変化なし) | 未達 |

誤差入りの閉軌道補正の 1 実現 (40 反復) は **目標 (< 60 s) を達成できていない**: 誤差入りの閉軌道は毎反復ほぼ全要素の入口 COD が変わるのでキャッシュが効かず、
要素ごとの Jacobian (Python) が支配的 (BEND の `_tchge_jac` 数値微分 + `_bend_k1_cod_matrix` で一周 1.5 s ×10 回/反復)。
280 回の cell calc の部分だけが 140 s → 約 10 s になった。

## 3. 未対応・見送り

- 誤差入りの要素ごとの Jacobian (BEND の tchge、K1≠0 BEND、薄レンズ) の numba 化。誤差入りの閉軌道補正・6D emit の繰り返しの本命。解析式 (`_tchge_jac` の tbrote) が数値微分のまま。
- 6D (`calc6d`/`emit6d`) の要素ごとの計算 (tturne: 1 周数秒)。off-momentum の cell calc を繰り返す計算や `Emittance[]` を繰り返す計算の本命。プロファイルが平らで (tradke、_tbendebody、_tbendecorr、piece、eye/copy が各 5〜10%) 小手先では 2 倍にならない。
- DynamicApertureSurvey (`dapert`/`track6n`): 1 粒子 1 周 810 µs で SAD の 1.5 倍。FLUC ON は乱数列が粒子順に消費される (SAD の癖。docs/sad-tracking-particle-order-quirk) ので、粒子ごとの並列化は結果を変える。そもそも SAD も 5000 周の DA に 785 s かかる。
- `match6.py`/`match6d.py` の累積行列の再利用: プロファイルでは Jacobian 計算より twiss/residual 計算が支配的で、要素キャッシュ (#2) の効果が既に入っている。

## 4. 検証

- `QUICK=1 run_regression.py`: 作業前 (索引だけ) / 作業後とも 50 ベンチ FAIL 0。pytest 164 passed (`--ignore` 3 本)。
- 出力の同一性: `mb_emiterr.sad` (誤差入り emit) は EMITX/EMITY/EMITZ/SIGE/COD が 1e-14 以内。
- キャッシュした配列は読み出し専用にしてあり、書き換えがあれば例外になる (回帰で例外なし)。

## 5. 最終確認

QUICK 回帰は 50 ベンチ FAIL 0、pytest 164 passed (最終コードで再確認)。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: この話題の記述は work_unimpl.md と docs/20〜35 に出てこない (grep で確認)、関係するモジュールに 11/06 以後この話題の変更の記録がない。この文書は 11/06 時点の記録のまま。
