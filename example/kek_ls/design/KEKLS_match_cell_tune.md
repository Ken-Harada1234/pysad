# KEKLS_match_cell_tune
種別: FFS の `fit`・`free`・`go`。フォルダ: `design`。元のスクリプト: matching2_01 (前半)。ラティス: `../lattice/kek_ls_v03_68.sad` (3 GeV、周長 570.72 m、標準セル 20 個)。フラグ: `ON RAD FLUC COD RFSW RADCOD`。

## 何を見る例か
標準セル (`ARCA`) の位相進みを (48.58/20, 17.62/20) に合わせる、セルのチューンのマッチング。`MKC` の分散 EX の目標 exmkc を 0.004、0.002、0 と変えて 3 回行う。

## スクリプト
`KEKLS_match_cell_tune.sad` (先頭に `../lattice/kek_ls_v03_68.sad` を貼った自己完結のファイル。台本の本体は `FFS;` から。説明のコメントは台本の先頭にある)。

## 結果の読み方
条件 `fit NX`・`NY`・`MKC EX exmkc`・`MKS EX 0`、変数 `free Q*` (四極の K1)。各回の結果が `exmkc target ... result: ring tune ...` の行、そのときの四極の K1 が `K1:` の行。ring tune は セル×20 で (48.58, 17.62) に戻る。

gs の出力の抜粋 (指定の行だけ。全体は `KEKLS_match_cell_tune.gs.txt`):
```
design cell phase advance nx ny =    2.4290011676     .8810039277
target                          =    2.429     .881
SF-SF phase advance (x, y) [2pi] =    1.6010565067     .4881335462
exmkc target     .004  result: ring tune   48.58   17.6200000094  EX(MKC)     .004  EX(MKS)  -2.2556598E-14
exmkc target     .002  result: ring tune   48.5800004670   17.6200031580  EX(MKC)     .0020000119  EX(MKS)  -1.82118326E-8
exmkc target     .0  result: ring tune   48.5799997789   17.6200050225  EX(MKC)  1.366838534E-8  EX(MKS)  -1.06953239E-8
```

## gs と pysad の比較
判定: **△** (比べた行 9、最大の相対差 5.6e-03、k64 モードの pysad)。最初の 1 点 (exmkc = 0.004) は全桁一致。2 点目以降は残差 (チューン、EX) は同じ程度 (1e-6〜1e-7) に収束するが、変数の K1 は 1e-3 違う (変数が条件より多く、解が一意でない。最適化の経路で決まる)。

所要時間: gs 1 秒、pysad 2 秒 (この PC。k64 モード)。

## 未確認・限界
- 元の台本の `free BPI6 K1`・`free BPI7 K1` は、v03_68 の偏向磁石の名前 BPI6H・BPI7H に一致せず、K1 が動かない (gs・pysad とも)。元の挙動を残し、既定は `freebend = 0`。`freebend = 1` で BPI6H・BPI7H も変数にする。
- 元の台本の EMITX の FitFunction を除いた (閉じていない ARCA では EMITX は NaN で何も効かない)。ただし NaN の FitFunction があると pysad は fit が動かない (元の台本のまま走らせた比較で確認。gs は動く)。pysad の FitFunction の扱いの違いの疑い (README の要確認)。
- 整理前後の比較: 元の matching2_01 (ラティスを貼り、DAsv を除き、FitFunction を除いたもの) の gs と、この例の gs で、1 点目の値は一致、2・3 点目のチューンは 5e-6 以内で一致。

## 更新 (v2.4.0)
FitFunction・パーサ・USE の DP の修正のあとで `.pysad.txt` を再生成した。この例の結果は修正の前後で変わらず、上の比較と判定のままである (EX・K1 の行を gs と比べ直した)。
