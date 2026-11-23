# 09_io_system 一覧
入出力とシステムの例。ラティスは使わない。ファイルは作業ディレクトリ (gs は実行時の一時フォルダ) に作り、最後に削除する。実行ごとに変わる値 (日時、一時ファイル名、CPU 時間) は表示しない。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| FFS_PRINT_OUTPUT_APPEND_READ_INPUT | FFS PRINT・OUTPUT・APPEND・CLOSE・READ・INPUT | ○ |
| fn_DateString_TimeUsed_Timing | fn_DateString_TimeUsed_Timing (DateString, TimeUsed, Timing) | ○ |
| fn_FORM_PageWidth | fn_FORM_PageWidth (`$FORM`, PageWidth, StandardForm) | ○ |
| fn_Fork_Shared | fn_Fork_Shared (Fork, OpenShared, Shared, SharedSize) | △ (SharedSize は gs と一致。Windows の pysad は子を作らない。v2 は WSL で子を作る (この例では未再実行)) |
| fn_OpenRead_OpenWrite_OpenAppend_Close | fn_OpenRead_OpenWrite_OpenAppend_Close (OpenRead, OpenWrite, OpenAppend, Close) | ○ |
| fn_Print_Write_WriteString | fn_Print_Write_WriteString (Print, Write, WriteString) | △ |
| fn_Read_ReadString | fn_Read_ReadString (Read, ReadString) | ○ |
| fn_StringToStream | fn_StringToStream (StringToStream) | ○ |
| fn_System_TemporaryName | fn_System_TemporaryName (System, TemporaryName, Directory, Environment) | △ (v2 は既定で実行しない。許可すると実行でき、戻り値は v2 が 0、gs と k64 モードが 1) |
| var_FORM_Input_Output | var_FORM_Input_Output (`$FORM`, `$Input`, `$Output`) | ○ |
| fn_Parallelize | fn_Parallelize (`Parallelize[f, n, lshare, npara]`。POSIX (Linux・WSL) は fork した子で npara 個まで同時に計算する。Windows は逐次) | `.pysad.txt` は Windows の出力、`.pysad_wsl.txt` は WSL の出力 (gs は `.gs.txt`) |
| fn_SharedMemory_Parallelize_System_Fork | fn_SharedMemory_Parallelize_System_Fork (OpenShared・Shared、Parallelize、Fork、System。`Fork[]` は Linux・WSL・macOS だけ、Windows は v2 がエラーにする。`System[cmd]` は既定で実行せず、`PYSAD_ALLOW_SYSTEM=1` で許可したときだけ実行して終了コードを返す) | gs の出力なし。`.pysad.txt` は Windows の出力 (WSL では `child wrote: written by child`、許可すると `System の戻り値 = 0`) |

○一致 / △差あり / ×gs 失敗 / ×pysad 未対応。

`Read` (区切り・リスト・式・WordSeparators)、`Write[$Output]`、`Close[{..}]`、`DateString[date]` は gs と一致する。残る差は `System` (v2 は既定で実行しない。`PYSAD_ALLOW_SYSTEM=1` で実行でき、戻り値は v2 が 0、gs と k64 モードが 1)、`Fork` (子を作らない。v2 は POSIX で子を作る (docs/28)。この例の範囲では未確認)。FFS のコマンド (OUTPUT・APPEND・READ・INPUT) と `$FORM` は一致した。
