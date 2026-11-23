# fn_OpenRead_OpenWrite_OpenAppend_Close
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
ファイルを書き込み (`OpenWrite`、既存ファイルは空になる)・追記 (`OpenAppend`)・読み込み (`OpenRead`) で開き、`Close` で閉じる。ファイル番号は実数。`ReadString` が `EndOfFile` を返すまで読む。`FileType`・`FileQ`・`DeleteFile`。

## スクリプト
```
a = OpenWrite[f]; Write[a, "line 1"]; Write[a, "line 2"]; Close[a];
b = OpenAppend[f]; Write[b, "line 3 (appended)"]; Close[b];
c = OpenRead[f]; While[(s = ReadString[c]) <=> EndOfFile, Print[s]]; Close[c];
CloseFileList: Close[{a, b}];   DeleteFile[f];   FileType[f];
```

## 結果の読み方
`OpenWrite` はファイル番号 (実数) を返す。書き込み後、3 行が読める (`line 1`、`line 2`、`line 3 (appended)`)。`OpenWrite` で開き直すと内容が置き換わる (`rewritten` だけ、次の `ReadString` は `EndOfFile`)。`FileType` は `File`、削除後は `None`、`FileQ` は 1→0。

## gs と pysad の比較
○ 一致する (`Close[{a, b}]` でリストの複数のファイルを閉じる)。

## 未確認・限界
- `OpenWrite` したままの `Close` 忘れ、読み書きの同時オープン、バイナリは見ていない。
- ファイル番号は実行ごとに変わりうるので表示していない。
