# FFS_COUPLE (COUPLE / DECOUPLE: 要素値の結合)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`COUPLE slave master coef`: slave の既定キーワードの値を **coef × master の値** に結ぶ (電源を共有する 2 つの磁石の表現)。以後 master を変えると slave が追従する。`DECOUPLE slave` (または `INDEPENDENT`) で解除。
ここでは QD を QF に `COUPLE QD QF -0.857142857…` で結び (QD/QF の比を一定にする)、QF の K1 を 0.7→0.8 にしたときに QD が -0.6857 に追従することを見る。`ATTRIBUTE` の `Couple` 欄で結合先が見える。

## スクリプト
```
USE RING; CELL; CALC;
COUPLE QD QF -0.857142857142857;
ATTRIBUTE QD;
QF K1 0.8;   CALC;
DECOUPLE QD;   QF K1 0.7;   CALC;
```

## 結果の読み方
- 結合前: K1 (QF, QD) = (0.7, -0.6)、チューン 6.696329 / 4.906056。
- `ATTRIBUTE QD`: gs は `Couple` の欄に `QF.1` と係数 `-.857143` を出す。pysad は `<-- 1.000000` (結合を表示しない)。
- `QF K1 0.8` のあと: QF = 0.8、QD = -0.685714 (= -0.857143 × 0.8)、チューン 7.908561 / 5.859714 (gs・pysad 同じ)。

## 物理・解析式
QD は QF に従う: $K_{1,QD}=-0.857143\,K_{1,QF}$ (= -0.6/0.7 比の保存)。この K1 のペアのチューンを 2x2 行列で計算すると:

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| QD K1 (QF=0.8) | -0.685714 | -0.685714 | -0.685714 |
| NX (QF=0.8) | 7.90856 | 7.90856 | 7.90856 |
| NY (QF=0.8) | 5.85971 | 5.85971 | 5.85971 |


## gs と pysad の比較
- 結合の動作 (QD が QF に追従する K1 と、そのチューン) は**一致**。
- **△ 差あり (2 点)。** (1) `ATTRIBUTE QD` の `Couple` 欄: gs は結合を出し、pysad は出さない。(2) `DECOUPLE QD` のあと `QF K1 0.7; CALC`: pysad は `?Missing value for QD` を出し (DECOUPLE を読めていない)、チューンは NX=6.499372、NY=6.112523 (K1 = 0.7, -0.6857 の 2x2 行列の値: 6.499372、6.112523 と一致)。gs は同じ K1 の表示 (0.7, -0.685714) なのに `Vertical Unstable.` を出し、NX=6.688009、NY=4.975707 と別の値になる。(K1=0.7, -0.6857 で垂直は安定 (tr/2=-0.0295) なので、gs のこの結果は式と合わない。DECOUPLE のあとの gs の内部状態 (結合係数の再適用など) が表示値と食い違っている可能性があるが、原因は未調査。)

## 未確認・限界
- `COUPLE *` (すべての結合を再適用)、`INDEPENDENT`、`ElementValues` による一般の結合 (SADHelp の COUPLE の項) は未確認。
- gs の DECOUPLE 後の挙動は上のとおり式と合わない (K64 の癖か、このスクリプトの使い方の問題かは未確認)。
