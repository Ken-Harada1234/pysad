# fn_SpecialFunctions_D (D.n 系の特殊関数)
種別: 関数 (特殊関数。`Packages/D.n` と K64 の組み込み) / フォルダ: 10_math_data / ラティスなし

## 何を見る例か
多項式 (`LegendreP`、`HermiteH`)、指数・三角積分 (`ExpIntegralEi`、`ExpIntegralE`、`SinIntegral`、`CosIntegral`)、Fresnel、Airy、`Zeta`、`PolyLog`、`PolyGamma`、`ProductLog`、`Hypergeometric1F1` を 1 点ずつ評価して式の値と比べる。

## スクリプト
```
LegendreP[3,.5]  HermiteH[3,.5]  ExpIntegralEi[1.]  ExpIntegralE[2,1.]  SinIntegral[1.]  CosIntegral[1.]
FresnelS[1.]  FresnelC[1.]  AiryAi[1.]  AiryBi[1.]  Zeta[2.]  PolyLog[2,.5]  PolyGamma[1,1.]  ProductLog[1.]  Hypergeometric1F1[1.,2.,.5]
```

## 結果の読み方
`R:` の行が値。gs の `PolyLog` は虚部が 1e-15 の複素数 `(.58224...-7.2E-16 I)` で出る (数値の丸め)。`HermiteH` の前に gs は `hgu -1.500 0.5000 0.2500` という行を出す (デバッグ表示と思われる)。

## 物理・解析式
$P_3(x)=\frac12(5x^3-3x)$、$H_3(x)=8x^3-12x$ (物理学者のエルミート)、$\mathrm{Ei}(x)=-\mathrm{PV}\!\int_{-x}^\infty e^{-t}/t\,dt$、$E_n(x)=\int_1^\infty e^{-xt}t^{-n}dt$、$\mathrm{Si}(x)=\int_0^x\frac{\sin t}t dt$、
$\mathrm{Ci}(x)=-\int_x^\infty\frac{\cos t}t dt$、$S(x)=\int_0^x\sin\frac{\pi t^2}2dt$、$\zeta(2)=\psi_1(1)=\frac{\pi^2}6$、$\mathrm{Li}_2(1/2)=\frac{\pi^2}{12}-\frac{\ln^22}2$、
$W(1)e^{W(1)}=1$、${}_1F_1(1;2;z)=(e^z-1)/z$。

| 項目 | 式の値 (scipy.special / 式) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| LegendreP[3,.5] $=(5x^3-3x)/2$ | -0.4375 | -.4375 | -.4375 | 0.0e+00 | 0.0e+00 |
| HermiteH[3,.5] $=8x^3-12x$ | -5 | NaN | -5 | NaN | 0.0e+00 |
| ExpIntegralEi[1.] = Ei(1) | 1.89511781636 | (.256476306443282-2.71642432200215 I) | 1.895117816355937 | - | 0.0e+00 |
| ExpIntegralE[2,1.] = $E_2(1)$ | 0.148495506776 | .148495506775922 | .148495506775922 | 3.7e-16 | 3.7e-16 |
| SinIntegral[1.] = Si(1) | 0.946083070367 | .946083070367183 | .946083070367183 | 1.2e-16 | 1.2e-16 |
| CosIntegral[1.] = Ci(1) | 0.337403922901 | .337403922900968 | .337403922900968 | 4.9e-16 | 4.9e-16 |
| FresnelS[1.] | 0.43825914739 | .438259147390355 | .438259147390355 | 6.3e-16 | 6.3e-16 |
| FresnelC[1.] | 0.779893400377 | .779893400376823 | .779893400376823 | 0.0e+00 | 0.0e+00 |
| AiryAi[1.] | 0.135292416313 | .135292416312881 | .135292416312881 | 3.5e-15 | 3.5e-15 |
| AiryBi[1.] | 1.20742359495 | 1.207423594952871 | 1.207423594952872 | 3.7e-16 | 3.7e-16 |
| Zeta[2.] $=\pi^2/6$ | 1.64493406685 | 1.644934066848226 | 1.644934066848226 | 2.7e-16 | 2.7e-16 |
| PolyLog[2,.5] $=\frac{\pi^2}{12}-\frac{\ln^22}2$ | 0.582240526465 | (.582240526465012-7.2460401836E-16 I) | .582240526465014 | - | 2.5e-15 |
| PolyGamma[1,1.] $=\psi_1(1)=\pi^2/6$ | 1.64493406685 | 1.644934066848226 | 1.644934066848226 | 2.7e-16 | 2.7e-16 |
| ProductLog[1.] (Lambert W) | 0.56714329041 | .567143290409784 | .567143290409784 | 2.0e-16 | 2.0e-16 |
| Hypergeometric1F1[1,2,.5] $=2(e^{.5}-1)$ | 1.2974425414 | 1.297442541400256 | 1.297442541400256 | 3.4e-16 | 3.4e-16 |

## gs と pysad の比較
- 13 値は一致 (表の相対差が 1e-14 以下)。`AiryBi`・`PolyLog` は末尾 1〜2 桁の違い。
- **`HermiteH[3,.5]`: gs は `NaN`** (真値 -5)。pysad は -5。
- **`ExpIntegralEi[1.]`: gs は複素数 (.2565-2.7164 I)** (真値は実数 1.8951)。pysad は 1.895118。
- この 2 件は gs (K64 ビルド) の誤りで、pysad は真値を返す (docs/21 の方針)。△。

## 未確認・限界
1 点ずつの確認のみ。引数の範囲 (負、大きい)、次数の高い場合は試していない。複素引数: 最初の版では pysad 未対応と書いたが、v2 で Airy・Bessel・Hypergeometric0F1・1F1・2F1 (実数のパラメータ) が複素引数に対応した (`docs/47_v2_changes.md`。`AiryAi[1+I]` が値を返すことを 2026-11-17 に確認)。複素の PolyLog・HypergeometricU・Zeta は v2 でも未対応 (`Zeta[2+I]` は TypeError で確認)。この例の表には複素引数を入れておらず、gs との比較はしていない。`HermiteH`・`ExpIntegralEi` の gs の誤りの条件 (どの引数で誤るか) は調べていない。
