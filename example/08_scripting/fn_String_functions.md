# fn_String_functions (character-string の項)
種別: 関数 (文字列) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
文字列関数: `StringLength` `StringJoin` と `//`、部分文字列 `s[n1,n2]`、`StringPosition` `StringReplace` `StringTrim` `StringFill` `StringMatchQ` (ワイルドカード) `ToUpperCase` `ToLowerCase` `Characters` `ToCharacterCode` `FromCharacterCode` `ToExpression` `Symbol` `SymbolName` `DigitQ` `LetterQ` `StringDrop` `StringInsert`、エスケープ。

## スクリプト
```
s = "Hello, SAD World";   StringLength[s];  s[1,5];  s[-5,-1];
StringPosition["abcbcbcbcb", "bcb"];  StringReplace["abcbcbcbc", {"bcb"->"xy", "cbc"->"pqrs"}];
StringFill["ab", "-", 7];  StringMatchQ["QF.1", "QF*"];  StringMatchQ["Q1", "Q%"];  StringMatchQ["QX", "Q{^FD}"];
ToCharacterCode["\101"];  "x=" // 1.5 // ", n=" // 3;
```

## 結果の読み方
`StringPosition` は重なりも数えて {{2,4},{4,6},{6,8},{8,10}} (SADHelp の例と同じ)、`StringReplace` は {"bcb"->"xy","cbc"->"pqrs"} で "axypqrsbc" (左から最初に合うもの)。`StringFill["ab","-",7]`="ab-----"、負なら "-----ab"。`//` は数も文字列にして連結 (`"x=1.5, n=3"`)。エスケープ `\t` `\"` `\\` `\101` (8 進) は gs で {9,34}、{92}、{65}。

## gs と pysad の比較
○ 一致する (2026-11-17 確認)。`StringMatchQ` は SADHelp (wildcards) の `*`、`%`、`{..}`、`{^..}`、`|` を gs と同じに扱い (`?` は普通の文字、大文字小文字を区別、空文字列は 0)、8 進エスケープ `"\101"` は "A" (65)。

## 未確認・限界
- `CharacterPosition` は gs でも評価されず (未実装か引数の書式が違う)、例から外した。
