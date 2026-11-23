# fn_Print_Write_WriteString
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`Print` (画面 = `$Output` へ、末尾に改行)、`Write[f, ...]` (ファイルへ、改行つき)、`WriteString[f, ...]` (改行なし)。書いたファイルを `OpenRead`・`ReadString` で読み戻して確かめる。

## スクリプト
```
Print["...", 1, " ", 2.5, " ", {1,2}, " ", x + y];   Print[];
f = TemporaryName[];  fd = OpenWrite[f];
Write[fd, "Write adds a newline ", 1, " ", {1,2,3}];
WriteString[fd, "WriteString: no newline "];  WriteString[fd, "(same line) ", 3.5, "\n"];
Write[fd, a + b];  Close[fd];
Write[$Output, "Write to $Output"];
```

## 結果の読み方
`Print` は引数を文字列にして連結 (`1 2.5 str {1,2} x+y`)、実数は SAD の既定の書式 (`1/3`→`.3333333333333333`、`1.5e10`→`1.5e+10`)。ファイルから読み戻すと 1 行目 `Write adds a newline 1 {1,2,3}`、2 行目 `WriteString: no newline (same line) 3.5`、3 行目 `a+b`。`$Output`=-1 で、`Write[$Output, ...]` は画面に出る。

## gs と pysad の比較
△ (1) `Print[]` (引数なし): gs は `Null` を表示、pysad は空行。(2) 整数値の実数 `123456789012` は gs `123456789012.`、pysad `123456789012` (表示の差だけ)。(3) `Write[$Output, ...]` と `WriteString[$Output, ...]` は gs と一致 (`$Output` = -1 は標準出力。`Write[$Output]` は gs と同じ `General::narg`)。ファイルへの書き出し・読み戻しは一致。

## 未確認・限界
- 文字コードは UTF-8 で、日本語の出力は確認していない (`Print` の中は ASCII)。
