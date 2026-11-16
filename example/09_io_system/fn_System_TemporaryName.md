# fn_System_TemporaryName (System, TemporaryName, Directory, Environment)
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`TemporaryName[]` は一時ファイル名 (呼ぶたびに違う)、`System[cmd]` はシェルコマンドを実行してリターンコードを返す。`Directory[]`、`Environment[]`。

## スクリプト
```
t1 = TemporaryName[];   r = System["echo from-shell > " // t1];
c = OpenRead[t1]; Print[ReadString[c]]; Close[c];   System["exit 3"] <> 0;   DeleteFile[t1];
```

## 結果の読み方
2 回呼ぶと別の名前 (`t1 <=> t2` が 1)、名前は文字列。gs で `System["echo from-shell > file"]` はファイルを作り (内容 `from-shell`)、リターンコードは 1 (**成功でも 0 ではなく 1 と返った**。SADHelp は "system return code" とだけ)、`System["exit 3"]` も 1 以外 (`<> 0` が真)。名前そのもの (`/tmp/...` の乱数部分) は実行ごとに変わるので表示していない。

## gs と pysad の比較
× **pysad は `System` を未対応 (意図的に実行しない)**: `[pysad] System[...] は実行しません` と表示してコード 0 を返し、ファイルは作られない (そのあとの `OpenRead` のエラーも `.pysad.txt` に残してある)。`TemporaryName`・`Directory`・`Environment` は一致。

## 未確認・限界
- `System` の結果は OS (gs は WSL の Ubuntu、パスは `/tmp/...`) に依る。
- 戻り値 1 の意味 (成功時) は未確認。
