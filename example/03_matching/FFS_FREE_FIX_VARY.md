# FFS_FREE_FIX_VARY (FREE・FIX・VARY・VARIABLES)
種別: FFS コマンド群。フォルダ: 03_matching。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
マッチングの**変数**の管理。
- `FREE パターン [キーワード] …`: 変数にする。キーワードを省くと既定キーワード (QUAD は K1、BEND は ANGLE、SEXT は K2、DRIFT は L)。同じ名前の成分はまとめて 1 変数 (ファミリー)。
- `FIX パターン [キーワード]`: 変数から外す。`FIX *` で全部外す。
- `VARY キーワード パターン`: その要素の既定キーワードを変える (例: QF の変数を K1 でなく ROTATE にする)。
- `VARIABLES` (略 `VAR`): 変数の表 (Variable・Keyword・Now・Previous・Saved・Minimum・Maximum・Couple・Coefficient)。

## スクリプト
```
USE RING; CELL; CALC;
FREE QF QD;     VARIABLES;
FIX QD;         VARIABLES;
FREE B;         VARIABLES;
VARY ROTATE QF;  FREE QF;   VARIABLES;
FIX *;          VARIABLES;
FREE QF K1;   FREE SF K2;   VARIABLES;
```

## 結果の読み方 (gs)
1. `FREE QF QD`: QF K1 = 0.7、QD K1 = -0.6 の 2 変数。最小・最大は ±1e10 (既定の限界)、`<--` は結合なしの印。
2. `FIX QD`: QF だけ残る。
3. `FREE B`: BEND は既定キーワード ANGLE (0.261799388) が変数に加わる。
4. `VARY ROTATE QF; FREE QF`: QF の既定キーワードが ROTATE になり、`QF ROTATE 0`、そして先に FREE していた QF K1 (非既定キーワードの変数。限界が ±1e99 と表示される) が並ぶ。
5. `FIX *`: 全部消える。
6. `FREE QF K1; FREE SF K2`: SF K2 と QF K1。

## 物理・解析式
なし (変数の管理)。

## gs と pysad の比較
変数の内容 (名前・キーワード・値) は一致。表示の差のみ: 4. で変数の**並び順**が違う (gs: ROTATE が先、pysad: K1 が先)。非既定キーワードの Minimum/Maximum の列は gs が `±1e99` で `<--` なし、pysad は `±1e10` と `<--` (K1 に既定の限界を付ける)。

## 未確認・限界
- `FREE` に `FREE QF.2` のような成分 (ファミリー先頭以外) を渡した場合 (docs/03 に `Info-Component … is coupled to`) は未確認。
- 実際に GO で非既定キーワード (ROTATE など) を動かす例は未確認 (VariableRange/Weight の例は fn_FitValue_... を参照)。
