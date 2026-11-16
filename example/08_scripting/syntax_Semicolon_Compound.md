# syntax_Semicolon_Compound (括弧・関数の中の ;)
種別: 文法 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
括弧・関数の中の `;` は `CompoundExpression` の区切り。最後の式の値が返り、**末尾が `;` なら `Null` が返る**。

## スクリプト
```
f[x_] := (a = x; a^2);     g[x_] := (a = x; a^2;);
f[3]  (* 9 *)   g[3]  (* Null *)
Module[{t = 2}, t; t^2];   Module[{t = 2}, t; t^2;];
If[1 > 0, a = 1; a + 10, a = 2; a + 20];   Do[a = 1; a, {2}];
CompoundExpression[1, 2, 3];   Table[a = i; a^2, {i, 3}];
```

## 結果の読み方
`f[3]`=9、`g[3]`=Null (`Head` は `Symbol`)。`Module[..., t; t^2]`=4、末尾 `;` で Null。`If` の枝の中の `a=1; a+10` は 11。`Do[...]` の値は Null。`Table[a = i; a^2, ...]` は {1,4,9}。`(Print["step1"]; Print["step2"]; 42)` は 42。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。
- 最上位の `;` と括弧の中の `;` の違いは syntax_Semicolon_Toplevel.md。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
