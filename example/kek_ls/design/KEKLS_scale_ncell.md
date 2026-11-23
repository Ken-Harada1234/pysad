# KEKLS_scale_ncell
種別: FFS の `fit`、`LINE[]` への代入、`LINE` の定義。フォルダ: `design`。元のスクリプト: matching21_v03_35 (v03_35 用。matching21・matching22・matching21_v03_35_fix は同系)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
標準セルの数を j (20, 18, 16, 12, 10) に変える。セル 1 個の偏向角を 20/j 倍し、四極 (と BPI6H・BPI7H) の K1 を調整してセルの位相進みを設計値に保ち、j セルのリング RINGJ の εx・チューン・周長を表示する。

## スクリプト
`KEKLS_scale_ncell.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
`cells j:` の行がセルの位相進みと βx(MKSFA1.1)・EX(MKC)、`ring RINGJ..:` の行がリングの周長・チューン・εx。εx は j = 20 の 134.3 pm rad から、18 で 184、16 で 261、12 で 615、10 で 1064 pm rad と、(20/j)³ で増える (偏向角の 3 乗)。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_scale_ncell.gs.txt`):
```
cells   20.0000000000: cell phase    2.4290011676     .8810039277  BX(MKSFA1.1)    8.3350203379  EX(MKC)     .0199999550
   ring RINGJ20: length  570.72  tune   48.5800233525   17.6200785537  EMITX [pm rad]  134.3162305767  SIGE  6.405157517E-4
cells   18.0000000000: cell phase    2.4290013056     .8810039784  BX(MKSFA1.1)    8.3214496772  EX(MKC)     .0221893396
   ring RINGJ18: length  513.648  tune   43.7220235008   15.8580716114  EMITX [pm rad]  184.3108969723  SIGE  6.752156659E-4
cells   16.0000000000: cell phase    2.4290019590     .8810042198  BX(MKSFA1.1)    8.3025747925  EX(MKC)     .0249114006
   ring RINGJ16: length  456.576  tune   38.8640313444   14.0960675171  EMITX [pm rad]  261.4223621971  SIGE  7.162601902E-4
cells   12.0000000000: cell phase    2.4290089958     .8810068619  BX(MKSFA1.1)    8.2336064284  EX(MKC)     .0329606386
   ring RINGJ12: length  342.432  tune   29.1481079502   10.5720823433  EMITX [pm rad]  614.8489290324  SIGE  8.274711234E-4
cells   10.0000000000: cell phase    2.4290011686     .8810039282  BX(MKSFA1.1)    8.1649326951  EX(MKC)     .0392413623
   ring RINGJ10: length  285.36  tune   24.2900116859    8.8100392818  EMITX [pm rad]  1063.993169439  SIGE  9.069090395E-4
```

## gs と pysad の比較
判定: **○** (比べた行 15、最大の相対差 0.0e+00、k64 モードの pysad)。全桁 (1e-10) 一致。

所要時間: gs 1 秒、pysad 3 秒 (この PC。k64 モード)。

## 未確認・限界
- 元は v03_35 用で、FitFunction が `(EMITX − 1*1^-9)*1E9` (`1^-9` は 1。10^-9 の誤り)、さらに ARCA では EMITX が NaN なので、元のままでは条件にならない (pysad は NaN の FitFunction で fit が動かない)。この例は、元の他の版でコメントにされていた `fit NX nxcell`・`fit NY nycell` (セルの位相進みを設計値に) を条件にした。
- USE し直すと ANGLE が元に戻るので、リングを USE した後に ANGLE と K1 を設定し直している (元の台本は ARCA だけで、εx は出さない)。
- j = 14 は、位相進みを保つとリングのチューンが 34.006 (整数に近い) になり `emit` が破綻するので外した。j = 20 → 2 の元の範囲では、整数・半整数共鳴に当たる j がほかにもある。
- RF 空洞は 1 台 (2 MV、ハーモニック数 952)。周長が変わるので RF 周波数・バケットは設計値と違う (εx には影響しない)。

