# KEKLS_calc_H
種別: 解析式 (H 関数の積分によるエミッタンス)。フォルダ: `display`。元のスクリプト: calc_H.sad、calc_H_2.sad、calc_emit.sad (後半)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
標準セル (ARCA) の偏向磁石 (勾配つきを含む) の中点の Twiss から、各磁石の H 積分と分散の積分を解析式で求め、εx = Cq γ² ΣH / Σ(θ²/L) / Jx を、`cell calc emit` の RINGC の EMITX と比べる。

## スクリプト
`KEKLS_calc_H.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
磁石ごとの H 積分の行の後に、`Jx = 1.31104`、H 積分からの εx = 131.18 pm rad、SAD の εx = 134.32 pm rad、比 0.9767。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_calc_H.gs.txt`):
```
Jx =    1.3110389414
emittance from H integral [pm rad] =  131.1802005125
EMITX of RINGC (cell calc emit) [pm rad] =  134.3162285432
ratio =     .9766519052
```

## gs と pysad の比較
判定: **○** (比べた行 4、最大の相対差 0.0e+00、k64 モードの pysad)。全桁一致 (1e-10 以内)。

所要時間: gs 1 秒、pysad 3 秒 (この PC。k64 モード)。

## 未確認・限界
- H 積分は SAD の εx より 2.3% 小さい。差の原因は未確認 (中点の Twiss から前向き・後ろ向きに積分する近似、E1・E2 の端の効果なし、減衰分配 Jx の評価など)。解析式を SAD の結果に合わせる補正はしていない。
- 元の calc_emit.sad は勾配なしの式で、γ = 2500/0.511 (3 GeV に合わない)、途中の `stop;` で後半は動かない。calc_H.sad の Cq は `3.832/10000` (nm の単位)。
- 整理前後の比較: 元の calc_H.sad (ラティスを貼り、NPARA だけ 4 にして) の gs の Jx = 1.3110389414、εx = 0.1311802005 (nm rad) が、この例の 131.1802005 pm rad と一致。
