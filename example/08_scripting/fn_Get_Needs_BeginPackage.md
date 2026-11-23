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
△ **ほぼ一致 (2026-11-17)。** pysad は `BeginPackage`・`Begin` で文脈を作る。`Stat`Mean2[{1,2,3,4}]` = 2.5、`Stat`Var2` = 1.25、`Stat`Calls[]` = 2、`Stat`Counter` = 2、`Context[Stat`sq]` = `Stat``、文脈なしの `Mean2[{1,2}]` が未評価、`Context[]` が読み込み前後とも `System`` は gs と同じ。文脈の規則は `docs/28_v2_vs_k64.md` の「文脈」。
違いは `Names["Stat`*"]` だけ: gs は `{"Stat`"}`、pysad v2 (既定) は文脈の記号の一覧 (`Stat`Calls`・`Stat`Counter`・`Stat`Mean2` ほか)。`PYSAD_COMPAT=k64` では gs と同じ `{"Stat`"}` を返す (確認済み)。v2 は意図して一覧を返す (docs/28)。

## 未確認・限界
- パッケージのファイルを example の中に置かず、スクリプトが一時的に書いて消す (`DeleteFile`)。
- `Needs[]` は試していない (SADHelp の関数一覧にない)。
- `$ContextPath` の操作 (`Begin` の入れ子など) は未確認。
