# fn_DateString_TimeUsed_Timing
種別: 関数 / フォルダ: 09_io_system / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。ファイルは `TemporaryName[]` か実行ディレクトリの作業ファイルで、最後に `DeleteFile` で消す。

## 何を見る例か
`DateString[]` は現在の日時 `"mm/dd/CCYY HH:MM:SS"` (19 文字)、`DateString[{Y,m,d,H,M,S}]` と `DateString[秒]` (1900-01-01 0:0:0 からの秒数) は日時の変換、`TimeUsed[]` は CPU 秒、`Timing[式]` は `{cpu 秒, 値}`。実行ごとに変わる値は表示せず、長さと型と大小関係だけを表示している。

## スクリプト
```
DateString[];  DateString[{2000, 1, 2, 3, 4, 5}];  DateString[0];  DateString[86400*365];
t0 = TimeUsed[];   r = Timing[s = 0; Do[s += Sin[i], {i, 20000}]; s];   TimeUsed[] >= t0;   Date[];
```

## 結果の読み方
`DateString[]` の長さ 19。`DateString[{2000,1,2,3,4,5}]`=`01/02/2000 03:04:05`、`DateString[0]`=`01/01/1900 00:00:00`、`DateString[86400*365]`=`01/01/1901 00:00:00` (365 日後)。`TimeUsed[]` は実数で、`Timing` は長さ 2 のリスト、値は 0.46196... ($\sum_{i=1}^{20000}\sin i$)。`Date[]` は長さ 6。

## gs と pysad の比較
○ 一致する (`DateString[date]` は 1900-01-01 からの秒、または {Y,m,d,H,M,S} を `mm/dd/yyyy HH:MM:SS` にする)。引数なし、`TimeUsed`、`Timing`、`Date` も一致 (`Timing` の合計は 16 桁目まで一致)。

## 未確認・限界
- 7 番目の要素 (offset) は無視している (未確認)。
- `Pause`、`Sleep` は見ていない。
