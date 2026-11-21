# PF / matching: マッチング (色収差補正・入射バンプ)

各例は、`NAME.sad` (台本。ラティスを先頭に含む)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明) の 4 点です。一致の凡例は ○ 一致、△ 差あり、× 未対応です。

| 例 | 内容 | gs と pysad |
|---|---|---|
| [`FFS_FIT_chrom_PF`](FFS_FIT_chrom_PF.md) | 色収差補正の6 極フィット (設計値が戻る) | ○ |
| [`FFS_FIT_injbump_PF`](FFS_FIT_injbump_PF.md) | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| [`FFS_FIT_modules_PF`](FFS_FIT_modules_PF.md) | `ChromCorrect[]`・`MatchingNcellTune[]` のモジュール化 | ○ |
