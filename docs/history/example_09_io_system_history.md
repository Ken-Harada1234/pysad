# example/09_io_system の履歴 (取り除いた記述)

md から、版・日付つきの修正・経過の記述を取り除いた。取り除いた元の行 (整理前) を、ファイルごとにそのまま残す。現在の説明は各 md を見る。

## example/09_io_system/README.md

整理前の 19 行目:

> 見つけた差は 2026-11-17 に `Read` (区切り・リスト・式・WordSeparators)、`Write[$Output]`、`Close[{..}]`、`DateString[date]` を gs に合わせて直した。残る差は `System` (v2 は既定で実行しない。`PYSAD_ALLOW_SYSTEM=1` で実行でき、戻り値は v2 が 0、gs と k64 モードが 1)、`Fork` (子を作らない。v2 は POSIX で子を作る (docs/47)。この例の範囲では未確認)。FFS のコマンド (OUTPUT・APPEND・READ・INPUT) と `$FORM` は一致した。

## example/09_io_system/fn_Fork_Shared.md

整理前の 19 行目:

> △ **Windows の pysad は子プロセスを作らない (2026-11-17、v2.0.0 の再実行。`.pysad.txt` はこの結果)。** `SharedSize[{1,2,3}]` は 32 になり gs と一致した (最初の版は 1。v2 で gs と同じ式にした)。`Fork[]` は Windows に fork が無いので、v2 は `???System::error: Fork[] is not available on this platform` を出し (k64 モードは従来どおり実数 1 を返す)、子の分岐を通らず、共有メモリは親のものだけで `{1,2,3}` のまま (gs は `{10,20,30}`)。Linux・WSL・macOS では v2 の `Fork[]` が `os.fork` で子を作る (`docs/47_v2_changes.md`。記録は WSL で確認)。**この例を WSL の pysad では再実行していない** (WSL の Python に numpy が無い)。要確認。

整理前の 24 行目:

> - 初版で `Wait[pid]` を使った例は gs が 120 秒止まったので外した。

## example/09_io_system/fn_System_TemporaryName.md

整理前の 17 行目:

> △ **既定では意図して実行しない (v2.0.0、2026-11-17 の再実行)。** pysad v2 は `System` を既定で実行せず、`$Failed` を返して案内を出す (`.pysad.txt` はこの状態。ファイルは作られず、そのあとの `OpenRead` のエラーが残る)。環境変数 `PYSAD_ALLOW_SYSTEM=1` (起動前) か `pysad.miscfn.allow_system(True)` で許可すると、シェルで実行して終了コード (成功 0) を返す。許可して走らせると、ファイルが作られ内容 (`from-shell`。Windows の `echo` で末尾に空白が付く) も読め、`System["exit 3"] <> 0` は真で gs と同じ。戻り値は **v2 では成功が 0、gs は 1**: `PYSAD_COMPAT=k64` + 許可なら 1 で gs と合う (確認済み。k64 で許可なしは v1 と同じ 0 と「実行しません」の表示)。この 1 と 0 の違いは K64 の戻り値の癖で、成功時の意味は未確認 (SADHelp は "system return code" とだけ)。`TemporaryName`・`Directory`・`Environment` は一致。

## example/09_io_system/fn_DateString_TimeUsed_Timing.md

整理前の 17 行目:

> ○ 一致する (2026-11-17 に修正: `DateString[date]` は 1900-01-01 からの秒、または {Y,m,d,H,M,S} を `mm/dd/yyyy HH:MM:SS` にする)。引数なし、`TimeUsed`、`Timing`、`Date` も一致 (`Timing` の合計は 16 桁目まで一致)。

## example/09_io_system/fn_OpenRead_OpenWrite_OpenAppend_Close.md

整理前の 19 行目:

> ○ 一致する (2026-11-17 に修正: `Close[{a, b}]` でリストの複数のファイルを閉じる)。

## example/09_io_system/fn_Read_ReadString.md

整理前の 18 行目:

> ○ 一致する (2026-11-17 に修正: gs の `tfreadstringf` を写した。区切り (`,` を含む既定の WordSeparators)、リスト形式 `{Word,Real,Real}`、式 `3*4` の評価、`WordSeparators`、`ReadNewRecord`、`ReadString` の読み残しが gs と同じ)。

## example/09_io_system/fn_StringToStream.md

整理前の 19 行目:

> ○ 一致する (2026-11-17 に fn_Read_ReadString と同じ修正で `n*item`・リスト形式・`Expression` が gs と同じ)。

## example/09_io_system/var_FORM_Input_Output.md

整理前の 17 行目:

> ○ 一致する (2026-11-17 に `Write[$Output, ...]` を修正)。`$FORM`・`$Input`・`$Output` の値は一致。


## 追加 (整理の 2 回目)

### example/09_io_system/fn_Print_Write_WriteString.md

整理前の 21 行目:

> △ (1) `Print[]` (引数なし): gs は `Null` を表示、pysad は空行。(2) 整数値の実数 `123456789012` は gs `123456789012.`、pysad `123456789012` (表示の差だけ)。(3) `Write[$Output, ...]` と `WriteString[$Output, ...]` は 2026-11-17 に修正 (`$Output` = -1 は標準出力。`Write[$Output]` は gs と同じ `General::narg`)。ファイルへの書き出し・読み戻しは一致。


## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### fn_Parallelize.sad (旧 v2_next/B2_parallelize.sad)
> (* v2next B2 項目 3: Parallelize[f, n, lshare, npara]。…
