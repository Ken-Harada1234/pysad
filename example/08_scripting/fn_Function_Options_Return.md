# fn_Function_Options_Return (オプション・既定値・Return・再帰)
種別: 関数 (ユーザー定義関数) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
ユーザー定義関数の書き方: オプション (`opt___Rule` と `Override[]`)、既定値つき引数、`Return`、`Break`、再帰 (階乗・最大公約数・ハノイ)、引数の個数による別定義、`SetAttributes[HoldAll]`、`Definition`。

## スクリプト
```
f1[x_, opt___Rule] := Module[{o = Override[opt, Scale->1, Offset->0]}, Scale*x + Offset /. o];
f2[x_, n_:2, m_:1] := x^n + m;
f3[l_List] := Module[{r = 0}, Scan[If[# < 0, Return["negative found"]] &, l]; Do[If[l[[i]] > 100, Return[i]], {i, Length[l]}]; "none"];
gcd[a_, 0] := a;  gcd[a_, b_] := gcd[b, Mod[a, b]];
SetAttributes[twice, HoldAll];  twice[e_] := (e; e);
```

## 結果の読み方
`f1[3]`=3、`f1[3, Scale->2]`=6、`f1[3, Offset->10, Scale->2]`=16。`f2[3]`=10、`f2[3,3]`=28、`f2[3,3,5]`=32。`gcd[48,18]`=6、`fact[10]`=3628800、`hanoi[6,1,2,3]`=63 ($2^6-1$)。`twice[k++]` は `k` を 2 回増やす (k=2)。引数の個数で定義を分けた例 (`f6`) は none / one / two / many を返す。

## gs と pysad の比較
△ (1) `Return` の範囲は 2026-11-17 に gs に合わせた (純関数・`Scan`・`Do`・`Module` を抜けて外側のユーザー定義関数まで届く。`f3[{1,-1}]` は "negative found")。(2) `Definition[fact]` は gs が `Hold[{fact:=fact, ...}]` を返し、pysad は `Null` (未対応のまま)。それ以外 (オプション、既定値、再帰、HoldAll、個数別定義) は一致。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
- `Override[opt, default...]` は `Packages/functions.n` 由来 (`emit.n` が使う書き方)。
