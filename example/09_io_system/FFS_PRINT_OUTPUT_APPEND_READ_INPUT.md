# FFS PRINT・OUTPUT・APPEND・CLOSE・READ・INPUT
種別: FFS コマンド / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`PRINT 式` (式を評価して表示)、`OUTPUT "file"` (出力を file へ切り替え。ファイルは先頭から書く)、`APPEND "file"` (追記)、`CLOSE OUTPUT` (前の出力に戻る)、`READ "file"`・`INPUT "file"` (file の内容をコマンドとして読む。READ は巻き戻す、INPUT は巻き戻さない)。

## スクリプト
```
PRINT 1 + 2;   OUTPUT "ex_out.txt";   PRINT 3 + 4;   Print["..."];   CLOSE OUTPUT;
APPEND "ex_out.txt";   PRINT 5 * 6;   CLOSE OUTPUT;     (* ex_out.txt を読んで表示 *)
READ "ex_cmd.txt";   INPUT "ex_cmd.txt";     (* x = 11; Print[...] を書いたファイル *)
```

## 結果の読み方
`PRINT 1+2` は `Out[2]:= 3` を表示 (最初の `PRINT` の呼び出しだけ、gs は余分に `Out[1]:= PRINT` も出す)。`OUTPUT` のあいだの `PRINT 3+4` と `Print[...]` は画面に出ず `ex_out.txt` に入る (`Out[3]:= 7`、`Print inside OUTPUT is also redirected`)。`CLOSE OUTPUT` で画面に戻る。`APPEND` した `Out[5]:= 30` が続けて入る。`READ`・`INPUT` でファイルの `x = 11;` と `Print[...]` が実行される (gs は読んだ行をエコー表示)。

## gs と pysad の比較
○ 一致 (pysad は `Out[1]:= PRINT` の行を出さない、読んだ行のエコー表示がない、という表示の差だけ)。

## 未確認・限界
- **ファイル名は文字列 `"..."` で書く**: `OUTPUT ex_out.txt;` (引用なし) は gs が `Missing filename for OUT_PUT` で失敗し、続く `CLOSE OUTPUT` が画面の出力を閉じて以降の出力が消える。
- `READ` と `INPUT` の巻き戻しの違いは、ファイルを 2 度読んで確かめてはいない (SADHelp の記述のみ)。
- `TERMINATE` は見ていない。
