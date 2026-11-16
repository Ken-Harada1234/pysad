# var_FORM_Input_Output ($FORM, $Input, $Output)
種別: 特殊変数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`$FORM` (実数の書式)、`$Input` (入力ストリームのファイル番号。既定は -1 ではなく実行中のファイル)、`$Output` (出力のファイル番号。-1 が画面)。

## スクリプト
```
Print[$FORM, $Input, $Output];   $FORM = "8.3"; Print[Pi, " ", 1/3];   $FORM = "S17.15";
Write[$Output, "writing to the console via $Output"];
```

## 結果の読み方
`$Output`=-1、`$FORM` の初期値は空白 1 文字、`$Input` は gs・pysad とも 77 (SADHelp は既定 -1 だが、スクリプトを読み込んでいる間はそのファイルの番号 77 になる)。どちらも `Head` は `Real`。`$FORM="8.3"` で `Pi`→`   3.142` (幅 8、小数 3)、`1/3`→`    .333`。`Write[$Output, ...]` は画面に出る。

## gs と pysad の比較
○ 一致する (2026-11-17 に `Write[$Output, ...]` を修正)。`$FORM`・`$Input`・`$Output` の値は一致。

## 未確認・限界
- `$Input` の 77 は gs の読み込みの実装に依る値で、環境により違いうる。
- `$Input` を使った `Read[$Input, ...]` (標準入力から読む) は、標準入力が /dev/null なので見ていない。
