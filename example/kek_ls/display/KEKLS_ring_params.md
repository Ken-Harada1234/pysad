# KEKLS_ring_params
種別: FFS の `cell calc emit` の表、`Twiss[]`。フォルダ: `display`。元のスクリプト: calc_emit.sad (前半)、disp、disp_emit、chromaticity_plot.sad。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
リング RINGC・RING のパラメータ (εx、σδ、σz、減衰、αc、U0、ξ) の表示。`cell calc emit` の表 (U0、αc、減衰時間・分配、バケット高さ) と、変数と差分で出した ξ (6 極あり・自然)。

## スクリプト
`KEKLS_ring_params.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
各リングの `===== NAME` の後に emit の表 (gs・pysad とも出る)、その後の `circumference`・`tune`・`EMITX`・`SIGE`・`SIGZ`・`xi` の行。ξ の自然値 (K2 = 0) は (−57.56, −89.26)、6 極ありは (+0.119, −0.684)。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_ring_params.gs.txt`):
```
===== RINGC
circumference [m]  =  570.72
tune               =   48.5800233525   17.6200785537
EMITX [pm rad]     =  134.3162285432
SIGE               =  6.405157569E-4
SIGZ [mm]          =    2.7343242295
xi (with sextupoles)    =     .1186928406    -.6844283403
xi (natural, K2 = 0)    =  -57.5648538915  -89.2609293962
===== RING
circumference [m]  =  570.7209201079
tune               =   48.5800240787   17.6200831582
EMITX [pm rad]     =  135.6808820236
SIGE               =  6.423281594E-4
SIGZ [mm]          =    2.7312326662
xi (with sextupoles)    =    -.3458708868    -.6501105636
xi (natural, K2 = 0)    =  -58.0716501603  -89.2207641679
```

## gs と pysad の比較
判定: **○** (比べた行 14、最大の相対差 0.0e+00、k64 モードの pysad)。全桁一致 (6e-9 以内)。

所要時間: gs 1 秒、pysad 4 秒 (この PC。k64 モード)。

## 未確認・限界
- ξ は `DP0 = ±0.0005` の差分。作業メモの値 (RINGC の ξ = (+0.133, −0.658)、自然 (−57.58, −88.96)) とは違う (差分の刻み 1e-4〜1e-3 では 0.118〜0.120 で変わらない。メモの値の出し方は未確認)。
- pysad の v2 モードは εx・σδ が K64 と 0.04% 違う。表の値の違いは lattice/README.md。
