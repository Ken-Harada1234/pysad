# FFS_RECOVER_SAVE_RESET (RECOVER・SAVE・RESET)
種別: FFS コマンド群。フォルダ: 03_matching。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
マッチングの結果を元に戻す・保存する操作。
- `RECOVER` (略 `REC`): FREE した変数を、**直前の `GO` の前の値と交換する** (もう一度 `RECOVER` すれば GO 後の値に戻る: トグル)。
- `SAVE [ALL] [パターン]`: 現在の値を「設計値」として保存する (`ALL` は全キーワード)。
- `RESET [ALL|TOTAL] [パターン]`: 値を保存してあった設計値に戻す。
`RECOVER` は GO の取り消し用、`SAVE`/`RESET` は任意の時点の値の保存・復帰用。

## スクリプト
```
FREE QF QD;   FIT $$$;   NX 6.7 NY 4.9;   GO;
RECOVER;   CALC;          ! GO 前の値に戻る
RECOVER;   CALC;          ! また GO 後の値
SAVE QF QD;   QF K1 0.5;   RESET QF QD;     ! 保存した値に戻る
```

## 結果の読み方
| 時点 | K1 (QF, QD) | NX, NY |
|---|---|---|
| GO の後 | (0.700221, -0.599636) | (6.7, 4.899998) |
| `RECOVER` | (0.7, -0.6) (元の値) | (6.696329, 4.906056) |
| `RECOVER` をもう一度 | (0.700221, -0.599636) | (6.7, 4.899998) |
| `SAVE` のあと `QF K1 0.5` | (0.5, -0.599636) | |
| `RESET QF QD` | (0.700221, -0.599636) (保存値) | |

## 物理・解析式
なし (値の保存・交換)。`RECOVER` で戻る値は FODO_ring の元の値で、チューンは (6.696329, 4.906056) の元に戻る (FFS_CALCULATE の値)。

## gs と pysad の比較
値はすべて一致。**pysad は `RESET QF QD` で `?Missing value for QF` / `?Missing value for QD` を出す** (結果の値は gs と同じに復帰するので、コマンドの読み方の一部が失敗している程度と思われる。△ 表示)。

## 未確認・限界
- `SAVE` なしで `RESET` した場合、`SAVE ALL`/`RESET ALL`/`RESET TOTAL`、MARK の保存値 (先頭 MARK の光学関数の保存) は未確認。
