# fn_StringToStream
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
文字列をファイルのように読む (`StringToStream[s]` がファイル番号を返す)。`Read`・`ReadString` でテキストを解析できる。

## スクリプト
```
c = StringToStream["10 20 30\nnext line\n3.5 abc"];
Read[c, 3*Real];   ReadString[c];   Read[c, {Real, Word}];   ReadString[c];
c = StringToStream["a 1\nb 2\nc 3"];  While[...Length[r = Read[c, {Word, Real}]] == 2, ...];
c = StringToStream["{1, 2, x+1} f[3]"];  Read[c, Expression];
```

## 結果の読み方
gs: ファイル番号は実数。`Read[c, 3*Real]`=102030 (3 個の実数 10, 20, 30 が連結して表示。`3*Real` は 3 個読む `n*item`)、`ReadString`=`next line`、`Read[c,{Real,Word}]`={3.5,"abc"}、続きがないと `EndOfFile`。`{Word,Real}` を 1 行ずつ読んで key=a val=1 など。`Read[c, Expression]` は {1,2,1+x} と `f[3]` (式に評価されて `x+1` は `1+x`)。

## gs と pysad の比較
○ 一致する (fn_Read_ReadString と同じ実装で `n*item`・リスト形式・`Expression` が gs と同じ)。

## 未確認・限界
- 文字列のエスケープ `\n` は gs で改行として働く。
