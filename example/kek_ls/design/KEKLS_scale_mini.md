# KEKLS_scale_mini
種別: FFS の `fit`・`FitFunction`・`LINE[]` への代入。フォルダ: `design`。元のスクリプト: matching23_v03_35_mini (v03_35 用)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
標準セル全体の長さを一様に factor 倍して縮小し (1 → 0.5)、K1 を 1/factor 倍にしてから、四極 (と、勾配つき偏向磁石 BPI6H・BPI7H) の K1 を再調整して、セルの位相進み (NX・NY) を元の値に保つ。

## スクリプト
`KEKLS_scale_mini.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
`design` が縮小前、`factor ...` の行が各縮小率でのセル長と位相進み、`BX*sin-0.3`・`sin_y` が FitFunction の 2 つの項 (6 極の対の位相条件)、`K1:` の行が再調整後の K1。factor 0.5 でセルは 14.268 m、四極の K1 は 2 倍前後。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_scale_mini.gs.txt`):
```
design: cell phase advance    2.4290011676     .8810039277  length   28.536
design: BX*sin(dmu_x)-0.3 =   -5.2438564772   sin(dmu_y) =     .0744900676
factor    1.0  length   28.536  phase    2.4290016591     .8810087818  BX*sin-0.3   -5.1314436899  sin_y  7.250357770E-9
factor     .9  length   25.6824  phase    2.4290011281     .8810083902  BX*sin-0.3   -4.6485049301  sin_y  -1.56265755E-8
factor     .8  length   22.8288  phase    2.4290010206     .8810085145  BX*sin-0.3   -4.1655287784  sin_y  2.022715614E-7
factor     .7  length   19.9752  phase    2.4290011812     .8810064067  BX*sin-0.3   -3.6824937063  sin_y  2.835676468E-7
factor     .6  length   17.1216  phase    2.4290011988     .8810056162  BX*sin-0.3   -3.1994084323  sin_y  1.524551278E-7
factor     .5  length   14.268  phase    2.4290012106     .8810049658  BX*sin-0.3   -2.7162749210  sin_y  6.201110970E-8
```

## gs と pysad の比較
判定: **○** (比べた行 14、最大の相対差 0.0e+00、k64 モードの pysad)。6 点とも位相進み (2.42900, 0.88100)・FitFunction の項・K1 が 2e-6 以内で一致。

所要時間: gs 1 秒、pysad 4 秒 (この PC。k64 モード)。

## 未確認・限界
- 元は v03_35 (1.5 GeV、別のセル) 用。v03_68 に合わせて、ライン名 `ARCA1` → `ARCA`、要素名 `BPI6`・`BPI7`・`QFA4A`・`QFA4B` → `BPI6H`・`BPI7H`・`QFA4`、`MKSFA1`・`MKSFA2` → `MKSFA1.1`・`MKSFA1.2` に直した。元の v03_35 は動かして確認していない (別のラティス)。
- 元の `free BPI6 K1` は名前が合わず効かないので `BPI6H`・`BPI7H` にした (`freebend = 0` で元の挙動)。
- 元の K1 の縮小は 1/factor 倍 (形のとおりなら 1/factor² 倍)。再調整で補う。
- 例は factor を 0.1 刻みの 6 点にした (元は 0.01 刻み 51 点)。元の範囲全部は `{factor,1,0.5,-0.01}`。

## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で 1e-8 (相対) 以内しか変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
