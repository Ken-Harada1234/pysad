# var_Null_Out (Null, Out[n], %, `$Line`)
種別: 特殊変数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Null` (値なし)、`Out[n]` と `%` (n 番目と直前の出力)、`$Line`。FFS の最上位で `;` のない式を打つと `Out[n]:= 値` と表示され、`%` で参照できる。

## スクリプト
```
x = Null;  Head[Null];  {1, Null}   (* 表示は {1,} *)
1 + 1            (* Out[1]:= 2 *)
Print[Out[1]];  Print[%];  Print[$Line];
2 * 3;  Print[%%];
```

## 結果の読み方
gs: `Null` の `Head` は `Symbol`、`Length[{Null,Null}]`=2、`{1,Null}` の表示は `{1,}`。`Out[1]` は 2、`%` は直前の出力 2、`$Line`=1 (入力行の数ではなく Out の番号付き行の数)。`Out[-1]` と `%%` は評価されずそのまま出る (SAD に相対参照はない)。`BoundQ[x]` は値があれば 1。

## gs と pysad の比較
○ `%`・`%n`・`Out[n]`・`$Line` は一致 (2026-11-17 に修正: `%` は `Out[$Line]`、`%%` は gs でも記号のまま、未定義の `Out[-1]` は未評価)。差が残るのは `BoundQ` (gs は `BoundQ[zzz]` を 0、`BoundQ[x]` を 1 と評価、pysad は未評価): 今回の対象外。

## 未確認・限界
- `Out[n]` の番号は、`Print` の行や `;` で終わる行には付かない。どの行が番号を持つかは gs の画面で確認する。
