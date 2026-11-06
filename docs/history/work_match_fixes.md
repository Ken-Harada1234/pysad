# 作業ログ: match6 / lattice の修正 — 現状記録 (2026-10-15)

担当範囲: `pysad/match6.py`, `lattice.py`。前担当 2 名が中断した作業の現状を記録する。コードは変更していない。
`pysad_20261013.tar.gz` (10/14 23:50 時点) を scratchpad の `old_pysad/` に展開して `diff -u` し、その上で CALC の時間と pf2 chrom を現行コードで測った。

## 1. diff の要点

- `match6.py`: **変更なし** (diff 空、mtime は 10/14 22:43 のまま)。
- `lattice.py` (+43 行): `Lattice.cell_initial` と `_twmov_to_marker` を追加。CELL CALC が安定なら、周期解の入口 Twiss (AX, BX, AY, BY, EX, EPX, EY, EPY, R1–R4, DETR) と閉軌道
  (DX, DPX, DY, DPY, DZ, DDP)、PSIX=PSIY=0 を先頭 MARK に書き戻す。SAD の `tffscalc` の `twmov(1,twiss,nlat,ndim,.false.)` の直訳。先頭が MARK でなければ何もしない。
- 関連: `ffscmd.py` の `row_frac` (qfraccomp、別文書)。`cell_initial` は **ffscmd.py / match6.py のどちらからも参照されていない** (grep で確認)。

## 2. 依頼項目ごとの判定

| 項目 | 判定 | 根拠 |
|---|---|---|
| (a) pf2 chrom の GO が NaN に発散 | 現行コードで解消 (match6.py は無変更なので、原因は間接的。どの変更が効いたかは未特定) | 下の 3 節 |
| (b) INS の初期軌道オフセットでチューンが違う (esrf amp_dep_ts 47.5924 vs 47.5800) | 未解決、途中 | `cell_initial` / MARK 書き戻しまでが入り、INS CALC 側がそれを入口条件に使う部分 (ffscmd の `ins_twiss` 呼び出し) が未着手。現行でも amp_dep_ts の 1 行目は 47.5924241431 (SAD 47.5800386899) |
| (c) CELL CALC の高速化 (3564 要素で 4 s → 0.5 s) | 差分なし (測定では既に目標以下で、「4 s」は再現せず) | 4 節 |

## 3. (a) pf2 chrom の現状

- 実行: `cd <work-pf2>; python -u ../run_sad.py chrom_r_full.sad` (ラティス込みのスクリプト。`pf_v10_16.lte` は Elegant 形式で SAD 用ではない)。2 分の打ち切りには掛からず **85 s で正常終了**、NaN の警告も `?No variable.` も出ない。
- 旧 pysad (tar 内) で同じ入力を 100 s 走らせると完走せず打ち切り。現行は完走する。
- 値の比較 (SAD は `chrom_r.sad.out`、`FIT; NX 11.6; NY 7.15; DP=0.001; FREE SD* SF* K2`):

| 量 | SAD | pysad | 相対差 |
|---|---|---|---|
| SD の K2 | -26.645487186 | -26.645683728 | 7e-6 |
| SFH の K2 | 23.916481408 | 23.916706565 | 9e-6 |
| 最終 NX / NY | 11.600007 / 7.150017 | 11.600007 / 7.150017 | 一致 |
| `Residual` (最終) | — | 7.9e-8 | — |
| tune (cell calc emit) | -0.3999939, 0.1500160, -0.0101703 | 同じ | 一致 |
| Equilibrium dz | 29.421510 mm | 29.423995 mm | 8e-5 (esrf disp と同じ傾向の差) |

- SAD の実行時間は 1.07 s。pysad は 85 s と約 80 倍遅い (GO の各反復が CALC 6 回分以上の評価を行うため)。
- 12:45〜12:57 の `chrom_dbg*`/`chrom_dbg3.pysad.out` は前担当のデバッグ記録。末尾の収束値は SD=-26.6456221, SFH=23.9166230 (残差 9.77e-8) で、上の値とほぼ同じ。

## 4. (c) CALC の所要時間 (pysad/tests/kek-ls_v03_12.lat、3564 要素、RING)

| 測定 | 旧 (tar) | 現行 |
|---|---|---|
| `twiss_table()` (`full=True`) | 0.072 s | 0.072 s |
| 同 (`full=False`) | — | 0.06 s |
| FFS 経由で `cell calc` を 3 回 (ラティス読込込み) | 未測定 (別の入力で失敗) | 0.5 s |
| `cell calc emit` 1 回 | — | 8〜12 s (測るたびにばらつく) |

`cProfile` では `_one_turn_cod_matrix` が 0.07 s、`_transport_twiss` が 0.03 s。つまり素の CALC は 0.5 s の目標を既に満たし、依頼にある「4 s」は
`CALC EMIT` のうち emit6d が占める時間を含んだ数字だった可能性が高い (確証はない)。時間が掛かるのは emit6d 側で、match6/lattice の範囲ではない。

## 5. 残り

1. (b): ffscmd の INS CALC が `lattice.cell_initial` (または先頭 MARK の書き戻し値) を入口条件に使うようにし、amp_dep_ts の 47.5800 を再現する。
   SAD では INS CALC の入口は MARK の値と `DX/DPX/DY/DPY` (LINE の値) から作られる (tffscalc)。
2. (a): 収束したが SD/SFH の K2 に 1e-5 の差が残る。dz の差 (8e-5) と同根かどうかは未調査。
3. (c): 測定上は不要。依頼の 4 s がどの操作を指すかを依頼者に確認する。
4. 回帰 (bench_level13a/13b、pytest 49 件) は確認済みとのことで再実行していない。

> 追記 (10/16、項目 6): `cell calc emit` の 8〜12 s は emit6d の数値ヤコビアンと閉軌道の反復が原因だった。解析ヤコビアン等で ESRF の emit は 1.6 s (RFSW+RADCOD)、`<work-esrf>/disp` 全体で 3.4 s になった (docs/history/work_item6.md §2)。
