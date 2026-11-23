# fn_Couple (Couple.n の ElementValues)

種別: パッケージ (Packages/Couple.n)・フォルダ 13_packages・ラティス FODO_ring (24 セル、3 GeV、周長 105.6 m)。pysad も動く

## 何を見る例か
`ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]}` で、要素の値同士を関係づける (QD の $K_1$ を QF の $-0.8$ 倍に結ぶ)。次の `CALC` で QD の $K_1$ が QF から決まる。QF/QD を同じ比で動かすなどの結合に使う。

## スクリプト
```
ElementValues={"K1"["QD"]->-0.8*"K1"["QF"]};
CALC;
Print["after: K1(QF.1)=",LINE["K1","QF.1"],"  K1(QD.1)=",LINE["K1","QD.1"]];
```

## 結果の読み方
`before:` は結合前 (QF 0.7、QD $-0.6$)、`after:` は結合後。`R:nx1`・`R:ny1` は結合後のチューン。

## gs と pysad の比較
| 量 | gs | pysad |
|---|---|---|
| $K_1$(QD) 結合後 | **-0.56** ($=-0.8\times0.7$) | **-0.56** |
| $\nu_x,\ \nu_y$ 結合後 | 6.788257, 4.306507 | 6.788257, 4.306507 |

- **○ 一致** (チューンは 6.78825745、4.3065066182 で全桁)。QD の $K_1$ が $-0.56$ に変わり、チューンも動く ($\nu_y$ は 4.906 から 4.307 へ。$K_1$ を 7% 弱めただけで大きく動くのは、$\nu_y$ が整数共鳴 ($\nu_y=4$) に近づくため)。
- pysad は `"K1"["QD"]` (文字列が頭部の式) を文字列の添字として評価すると `invalid literal for int() with base 10: 'QD'` で落ちるので、評価せず式のままにし、CALC の前に Couple.n の `SetCoupledElements` を SADScript のまま実行する (`With` が `Hold`・純関数の本体に値を埋め込む、`Hold[…][[1]]` の評価し直しなどを gs に合わせてある。docs/06)。
- `USE` で ElementValues が消えるのも gs と同じ (消えたあと CALC すると QD は -0.6 に戻る)。

## 未確認・限界
gs で `LINE["K1","QF"]=.75;` と QF を書き換えたあとに `CALC` すると QD が -0.6 に戻った (結合が保たれない。ElementValues の更新の仕組みは未確認)。`SetCoupledElements`、`ClearEV`、複数の結合は試していない。作動点の共鳴: この例の $\nu_y=4.31$ は半整数 4.5 からは離れているが、整数 4 には近いので、実際の設計では避ける値。
