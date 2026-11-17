# 作業ログ: match6 / lattice の修正 — 現状記録 (2026-10-15)

最終更新: 2026-11-16

担当範囲: `pysad/match6.py`, `lattice.py`。前担当 2 名が中断した作業の現状を記録する。コードは変更していない。
`pysad_20261013.tar.gz` (10/14 23:50 時点) を scratchpad の `old_pysad/` に展開して `diff -u` し、その上で CALC の時間を現行コードで測った。

## 1. diff の要点

- `match6.py`: **変更なし** (diff 空、mtime は 10/14 22:43 のまま)。
- `lattice.py` (+43 行): `Lattice.cell_initial` と `_twmov_to_marker` を追加。CELL CALC が安定なら、周期解の入口 Twiss (AX, BX, AY, BY, EX, EPX, EY, EPY, R1–R4, DETR) と閉軌道
  (DX, DPX, DY, DPY, DZ, DDP)、PSIX=PSIY=0 を先頭 MARK に書き戻す。SAD の `tffscalc` の `twmov(1,twiss,nlat,ndim,.false.)` の直訳。先頭が MARK でなければ何もしない。
- 関連: `ffscmd.py` の `row_frac` (qfraccomp、別文書)。`cell_initial` は **ffscmd.py / match6.py のどちらからも参照されていない** (grep で確認)。

## 2. 依頼項目ごとの判定

| 項目 | 判定 | 根拠 |
|---|---|---|
| INS の初期軌道オフセットでチューンが違う | 未解決、途中 | `cell_initial` / MARK 書き戻しまでが入り、INS CALC 側がそれを入口条件に使う部分 (ffscmd の `ins_twiss` 呼び出し) が未着手 |
| CELL CALC の高速化 (3564 要素で 4 s → 0.5 s) | 差分なし (測定では既に目標以下で、「4 s」は再現せず) | 4 節 |

## 3. CALC の所要時間 (3564 要素のリング)

| 測定 | 旧 (tar) | 現行 |
|---|---|---|
| `twiss_table()` (`full=True`) | 0.072 s | 0.072 s |
| 同 (`full=False`) | — | 0.06 s |
| FFS 経由で `cell calc` を 3 回 (ラティス読込込み) | 未測定 (別の入力で失敗) | 0.5 s |
| `cell calc emit` 1 回 | — | 8〜12 s (測るたびにばらつく) |


`cProfile` では `_one_turn_cod_matrix` が 0.07 s、`_transport_twiss` が 0.03 s。つまり素の CALC は 0.5 s の目標を既に満たし、依頼にある「4 s」は
`CALC EMIT` のうち emit6d が占める時間を含んだ数字だった可能性が高い (確証はない)。時間が掛かるのは emit6d 側で、match6/lattice の範囲ではない。

## 4. 残り

1. INS CALC が `lattice.cell_initial` (または先頭 MARK の書き戻し値) を入口条件に使うようにする。
   SAD では INS CALC の入口は MARK の値と `DX/DPX/DY/DPY` (LINE の値) から作られる (tffscalc)。
2. 測定上、CELL CALC の高速化は不要。依頼の 4 s がどの操作を指すかを依頼者に確認する。
3. 回帰 (bench_level13a/13b、pytest 49 件) は確認済みとのことで再実行していない。

> 追記 (10/16、項目 6): `cell calc emit` の 8〜12 s は emit6d の数値ヤコビアンと閉軌道の反復が原因だった。解析ヤコビアン等で 3564 要素のリングの emit は 1.6 s (RFSW+RADCOD) になった (docs/history/work_item6.md §2)。

## 2026-11-16 の確認・追記

- 2026-11-16 確認: 変化なし。根拠: この話題の記述は work_unimpl.md と docs/22〜35 に出てこない (grep で確認)、関係するモジュールに 11/06 以後この話題の変更の記録がない。この文書は 11/06 時点の記録のまま。
