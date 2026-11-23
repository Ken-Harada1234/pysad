# example/06_collective の md から取り除いた履歴の記述
例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきでそのまま残す (行番号は整理前のもの)。取り除いたあとの現在の記述は、各 md を参照。
整理日: 2026-11-23。

## flag_INTRA.md

元の 42 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## fn_TouschekLifetime.md

元の 40 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### flag_WAKE_command_table.sad (旧 v2_C/C2w_wake_table.sad)
> ! v2 (担当 C2): 旧式の WAKE コマンドの表を行をまたいで書く (K64 の tfwake.f の構文)。
> ! この例の C2w_wake_table.pysad.txt は pysad の v2 モード (既定) の出力。…
旧い作業符号「担当 C2」と旧い出力名 C2w_wake_table.pysad.txt を、現在のファイル名に直した。
