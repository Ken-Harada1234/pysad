# FFS_REJECT (REJECT: 条件を外す)
種別: FFS コマンド。フォルダ: 03_matching。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`REJECT 関数名` は、いまの FIT 位置の該当する条件だけを外す (`REJECT TOTAL` は全位置の全条件を外す)。
NX=6.7 と NY=4.9 の 2 条件のうち NY を外して `GO` すると、条件 1 つに対して変数が 2 つ (QF、QD) なので解が一意に決まらず、**最小ノルムの解** (K1 の変化が最小) が選ばれる。

## スクリプト
```
FREE QF QD;     FIT $$$;     NX 6.7 NY 4.9;     SHOW;
REJECT NY;      SHOW;        GO;
Print["after GO (NY rejected): NX NY = ", Twiss["NX","$$$"]/2/Pi, " ", Twiss["NY","$$$"]/2/Pi, "  K1 QF QD = ", ...];
REJECT TOTAL;   SHOW;
```

## 結果の読み方
- 最初の `SHOW`: NX と NY の 2 行。`REJECT NY` のあとは NX だけ。`REJECT TOTAL` のあとは条件の表が空。
- `GO` のあと (gs): NX = 6.7 (目標どおり)、NY = 4.904387 (QF,QD の最小変化の副産物)、K1 (QF, QD) = (0.700277, -0.599951)。

## 物理・解析式
チューンの K1 に対する感度行列 $J=\partial(\nu_x,\nu_y)/\partial(K_{QF},K_{QD})$ (2x2 行列の模型から数値微分) を使い、条件が NX だけのときの最小ノルム解は
$$\Delta\mathbf k=\frac{\mathbf j_x\,\Delta\nu_x}{\mathbf j_x\!\cdot\!\mathbf j_x}\quad (\mathbf j_x=\text{J の第 1 行},\ \Delta\nu_x=6.7-6.696329)$$
(線形近似。重みが等しい (ABSW) 変数として)。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| K1(QF) | 0.700277 | 0.700277 | 0.700277 |
| K1(QD) | -0.599951 | -0.599951 | -0.599951 |
| NY | 4.90439 | 4.90439 | 4.90439 |

J = [[12.8511, 2.2976], [-3.4266, -14.553]]。線形の最小ノルム解と gs は K1 で 6 桁一致 (K1 の変化が 3e-4 と小さいので非線形の補正が見えない)。

## gs と pysad の比較
全桁一致。

## 未確認・限界
- `REJECT` の位置依存 (現在の FIT 位置の条件だけ外す)、`REJECT TOTALFIT` は未確認。
