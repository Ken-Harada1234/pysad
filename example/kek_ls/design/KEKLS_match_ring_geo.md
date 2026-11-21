# KEKLS_match_ring_geo
種別: FFS の `fit GX`・`GY`、ジオメトリ。フォルダ: `design`。元のスクリプト: matching2_02_geo (ユーザーの指定は matching_02_geo)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
入射セルの長直線部の長さ `DRID1INJ` を変数にして、リング `RING` の周回後の位置 (GX, GY) を 0 にし、リングを幾何学的に閉じる。標準の 2.94 m から始めて、v03_68 の設計値 2.765460054 m に戻ることを見る。

## スクリプト
`KEKLS_match_ring_geo.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
`design` は設計値の閉じ方 (GX・GY ≈ 1e-7)、`start` は 2.94 m にしたとき (GX 0.282 m、GY 0.205 m)、`result` はマッチング後 (DRID1INJ 2.765460041 m、GX・GY ≈ 5e-8)、周長 570.72092 m。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_match_ring_geo.gs.txt`):
```
design : DRID1INJ    2.7654600539  closing GX GY  -2.83718165E-8  8.340120107E-8
start  : DRID1INJ    2.94  closing GX GY     .2824115368     .2051840958
result : DRID1INJ    2.7654600409  closing GX GY  -4.94618475E-8  6.807840382E-8  circumference  570.7209200818
```

## gs と pysad の比較
判定: **○** (比べた行 3、最大の相対差 0.0e+00、k64 モードの pysad)。全項目 5e-10 以内で一致。

所要時間: gs 1 秒、pysad 2 秒 (この PC。k64 モード)。

## 未確認・限界
- 元の `USE RINGI` は v03_68 にはなく (RING)、元の台本のまま走らせると gs が `Missing beamline in USE` で止まる (確認済み)。`RINGI` を `RING` に直した。
- 元の台本は設計値から始めるので、何も動かない。始めの値を 2.94 m にした手順は、この例で足した。
- 元の台本の `disp geo $$$` は、gs では `$$$` 1 行だけの表になる (要素ごとのジオメトリの表は出ない)。
