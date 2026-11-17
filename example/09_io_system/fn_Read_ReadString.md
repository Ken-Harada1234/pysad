# fn_Read_ReadString
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`ReadString[f]` は次の 1 レコード (1 行) を文字列で返す。`Read[f, item]` は `Word`・`Real`・`Expression`・`Character` と、それらのリストを読む。オプション `WordSeparators`・`ReadNewRecord`。

## スクリプト
```
Read[c, Word];   Read[c, {Word, Real, Real}];   Read[c, Expression];   Read[c, {Real, Real}];
Read[c, Character];   ReadString[c];
Read[c, Word, WordSeparators -> ":"];   Read[c, Word, ReadNewRecord -> False, WordSeparators -> ":"];
```

## 結果の読み方
ファイルは 5 行 (`alpha beta, gamma  12.5 7` / `x+y  3*4` / `1 2 3` / `4 5 6` / `end of file`)。gs の読み方: `Word` は `alpha`、次は `beta` (既定の区切りは空白・`,`・タブ)、`{Word,Real,Real}`={"gamma",12.5,7}、`Expression`=`x+y`、`{Real,Real}`={12,1} (`Real` は式も読むので `3*4` が 12、続く 1 は次の行から)、`{Real,Real,Real}`={2,3,4}、`Character`=5、`ReadString` は読み残し `" 6"`、最後は `end of file`。`WordSeparators -> ":"` で `a:b:c` が `a`、`b` に分かれる。**オプション名は `WordSeparators` (複数形。SADHelp の表では `WordSeparator` と書かれているが、単数形は gs で `Wrong option`)**。

## gs と pysad の比較
○ 一致する (2026-11-17 に修正: gs の `tfreadstringf` を写した。区切り (`,` を含む既定の WordSeparators)、リスト形式 `{Word,Real,Real}`、式 `3*4` の評価、`WordSeparators`、`ReadNewRecord`、`ReadString` の読み残しが gs と同じ)。

## 未確認・限界
- この例では `NullWords`、`Read[$Input, ...]` は見ていない。`n*item` (`3*Real`) は `fn_StringToStream` と試験 (`pysad/tests/test_scriptgs.py`) で確認した。
- 式の切れ目 (`x+y  3*4` で切る) は gs の例 2 つに合わせた経験則で、未確認の形がありうる。
