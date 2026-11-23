# KEKLS_cell_calc_emit
種別: FFS の `cell calc emit`。フォルダ: `lattice`。元のスクリプト: (全台本に共通の冒頭。個別の元の台本は無い)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
周回ライン `RINGC` (標準セル 20 個) と `RING` (入射セル込み) の周期解 (`cell calc emit`)。周長・チューン・εx・σδ・σz を、`LINE[]`・`Twiss[]`・変数 `EMITX`・`SIGE`・`SIGZ`・`MOMENTUM` で読む。

## スクリプト
`KEKLS_cell_calc_emit.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
`USE` するライン名を `FFS["USE " // nm]` で切り替えて 2 回回す。周長は `LINE["S","$$$"]` (終点の S)、チューンは `Twiss["NX","$$$"]/2/Pi`。RING は RINGC より 0.9 mm 長く (入射セルの `DRID1INJ`)、εx・σδ が少し大きい。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_cell_calc_emit.gs.txt`):
```
LINE
LINE RINGC
circumference [m]   =  570.72
tune nx ny          =   48.5800233525   17.6200785537
EMITX [m rad]       =  1.34316229E-10
SIGE                =  6.405157569E-4
SIGZ [m]            =     .0027343242
LINE RING
circumference [m]   =  570.7209201079
tune nx ny          =   48.5800240787   17.6200831582
EMITX [m rad]       =  1.35680882E-10
SIGE                =  6.423281594E-4
SIGZ [m]            =     .0027312327
```

## gs と pysad の比較
判定: **○** (比べた行 14、最大の相対差 0.0e+00、k64 モードの pysad)。全項目が全桁一致。

所要時間: gs 1 秒、pysad 3 秒 (この PC。k64 モード)。

## 未確認・限界
- pysad の v2 (既定) モードは εx・σδ が K64 と 0.04% 違う (RINGC εx 134.37 pm rad、σδ 6.4045e-4。K64 のバグ B75 の修正。`docs/28_v2_vs_k64.md`)。gs と同じ値は `PYSAD_COMPAT=k64`。
- エネルギー 3 GeV、周長 570.72 m、RF 2 MV、ハーモニック数 952 の設計値のラティス。誤差なし。
