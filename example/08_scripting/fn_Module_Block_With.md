# fn_Module_Block_With (局所変数のスコープ)
種別: 関数 (モジュール化) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Module` (字句的スコープ: 局所変数)、`Block` (動的スコープ: 値の一時差し替え)、`With` (定数の置換)、入れ子、`Module` が作る局所記号。

## スクリプト
```
x = 10;
m1[] := Module[{x = 1, y}, y = x + 1; x = x + 100; {x, y}];
show[] := x;
b1[] := Block[{x = 7}, show[]];     m2[] := Module[{x = 7}, show[]];
With[{a = 3, b = 4}, Sqrt[a^2 + b^2]];   n1[] := Module[{a = 1}, Module[{a = 2}, a] + a];
```

## 結果の読み方
`m1[]`={101,2} で大域の `x` は 10 のまま。`b1[]`=7 (`Block` の中では `show[]` も差し替えた x=7 を見る)、`m2[]`=10 (`Module` の局所 `x` は `show[]` から見えない)。`With`=5。入れ子の `Module` は 2+1=3。`Module[{t}, t]` が返す局所記号は大域の `t` と別物 (gs で `c[]===t` が 0)。

## gs と pysad の比較
○ 一致する (2026-11-17 に修正: `Module` を字句スコープにした。局所変数を `名前$N` に改名して評価するので、`show[]` を `Module[{x = 7}, ...]` の中で呼んでも大域の 10 を返し、`Module[{t}, t] === t` は 0)。局所記号の表示は gs が `t$1073741824` のような番号、pysad は連番 (表示だけの差)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
