# fn_Get_Needs_BeginPackage (Get, BeginPackage, Begin, End, EndPackage, 文脈)
種別: 関数 (パッケージ化) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
補助の `.n` (パッケージ) をスクリプト自身が `OpenWrite`・`WriteString` で書き出し、`Get[]` で読み込む。`BeginPackage[Stat`]`・`Begin[Stat`]`・`End[]`・`EndPackage[]` で文脈 `Stat`` に関数・変数を定義し、`Stat`Mean2[..]` のように文脈つきで呼ぶ。

## スクリプト
```
fd = OpenWrite["ex_stat_pkg.n"];
WriteString[fd, "BeginPackage[Stat`]\n"]; WriteString[fd, "Begin[Stat`]\n"];
WriteString[fd, "Mean2[l_List] := (Counter++; Plus @@ l/Length[l]);\n"]; ...
WriteString[fd, "End[]\n"]; WriteString[fd, "EndPackage[]\n"]; Close[fd];
Get["ex_stat_pkg.n"];   Print[Stat`Mean2[{1,2,3,4}]];
```

## 結果の読み方
gs: `Stat`Mean2[{1,2,3,4}]`=2.5、`Stat`Var2`=1.25、`Stat`Counter`=2 (`Mean2` が 2 回呼ばれた)、`Context[Stat`sq]` は `Stat``。文脈なしの `Mean2[{1,2}]` は評価されない (`Stat`` は `$ContextPath` に入らない: この SAD の `BeginPackage` は `Class.n` の例のように文脈名で全て参照する)。`Context[]` は読み込み後も `System``。`Needs[]` は SADHelp にもなく、使っていない (別ファイルは `Get[]` で読む)。

## gs と pysad の比較
△ **pysad は文脈を無視する (未対応と判断)**: `BeginPackage`・`Begin` の文脈を作らないので、定義は大域に入り、`Stat`Mean2[...]` は評価されず (未定義の記号)、文脈なしの `Mean2[{1,2}]` が 1.5 と評価される。`Context[Stat`sq]` は `` ` ``、`Names["Stat`*"]` は {}。gs の動作は 2026-11-17 に調べた: 文脈の中で読んだ新しい記号は `Stat`名` になり (`Stat`Mean2[...]` で呼べる)、`End[]`・`EndPackage[]` のあとは文脈なしの `Mean2` は見えず (`$ContextPath` に `Stat`` が入らない)、`Context[]`・`$Context`・`$ContextPath`・`Names["Stat`*"]` (= {"Stat`"}) の値が独特。読み込み時の名前解決の模倣が要り、SAD の Packages の読み込みを壊しうるので実装を見送った (docs/04 §6)。

## 未確認・限界
- パッケージのファイルを example の中に置かず、スクリプトが一時的に書いて消す (`DeleteFile`)。
- `Needs[]` は試していない (SADHelp の関数一覧にない)。
- `$ContextPath` の操作 (`Begin` の入れ子など) は未確認。
