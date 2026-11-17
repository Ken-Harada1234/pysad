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
△ **既定では意図して実行しない (v2.0.0、2026-11-17 の再実行)。** pysad v2 は `System` を既定で実行せず、`$Failed` を返して案内を出す (`.pysad.txt` はこの状態。ファイルは作られず、そのあとの `OpenRead` のエラーが残る)。環境変数 `PYSAD_ALLOW_SYSTEM=1` (起動前) か `pysad.miscfn.allow_system(True)` で許可すると、シェルで実行して終了コード (成功 0) を返す。許可して走らせると、ファイルが作られ内容 (`from-shell`。Windows の `echo` で末尾に空白が付く) も読め、`System["exit 3"] <> 0` は真で gs と同じ。戻り値は **v2 では成功が 0、gs は 1**: `PYSAD_COMPAT=k64` + 許可なら 1 で gs と合う (確認済み。k64 で許可なしは v1 と同じ 0 と「実行しません」の表示)。この 1 と 0 の違いは K64 の戻り値の癖で、成功時の意味は未確認 (SADHelp は "system return code" とだけ)。`TemporaryName`・`Directory`・`Environment` は一致。

## 未確認・限界
- `System` の結果は OS (gs は WSL の Ubuntu、パスは `/tmp/...`) に依る。
- 戻り値 1 の意味 (成功時) は未確認。
