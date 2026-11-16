# 09_io_system 一覧
入出力とシステムの例。ラティスは使わない。ファイルは作業ディレクトリ (gs は実行時の一時フォルダ) に作り、最後に削除する。実行ごとに変わる値 (日時、一時ファイル名、CPU 時間) は表示しない。

| ファイル | 説明 | gs と pysad |
|---|---|---|
| FFS_PRINT_OUTPUT_APPEND_READ_INPUT | FFS PRINT・OUTPUT・APPEND・CLOSE・READ・INPUT | ○ |
| fn_DateString_TimeUsed_Timing | fn_DateString_TimeUsed_Timing | ○ |
| fn_FORM_PageWidth | fn_FORM_PageWidth ($FORM, PageWidth, StandardForm) | ○ |
| fn_Fork_Shared | fn_Fork_Shared (Fork, OpenShared, Shared, SharedSize) | × |
| fn_OpenRead_OpenWrite_OpenAppend_Close | fn_OpenRead_OpenWrite_OpenAppend_Close | ○ |
| fn_Print_Write_WriteString | fn_Print_Write_WriteString | △ |
| fn_Read_ReadString | fn_Read_ReadString | ○ |
| fn_StringToStream | fn_StringToStream | ○ |
| fn_System_TemporaryName | fn_System_TemporaryName (System, TemporaryName, Directory, Environment) | × |
| var_FORM_Input_Output | var_FORM_Input_Output ($FORM, $Input, $Output) | ○ |

○一致 / △差あり / ×gs 失敗 / ×pysad 未対応。

見つけた差は 2026-11-17 に `Read` (区切り・リスト・式・WordSeparators)、`Write[$Output]`、`Close[{..}]`、`DateString[date]` を gs に合わせて直した。残る差は `System` (意図的に実行しない)、`Fork` (子を作らない)。FFS のコマンド (OUTPUT・APPEND・READ・INPUT) と `$FORM` は一致した。
