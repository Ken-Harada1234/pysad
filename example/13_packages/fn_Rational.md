# fn_Rational (rational.n の Rational・ContinuedFraction・SmallPrime)

種別: パッケージ (Packages/rational.n)・フォルダ 13_packages・ラティスなし

## 何を見る例か
分数 `Rational[n,d]` の演算 (和・積)、`Numerator`・`Denominator`・`FromRational`、連分数 `ContinuedFraction[x,n]`・`FromContinuedFraction`、素数表 `SmallPrime[i]`。連分数は、チューンの有理数近似に使える。

## スクリプト
```
Rational[1,3]+Rational[1,6]     Rational[1,3]*Rational[3,4]     Numerator[Rational[6,8]]     FromRational[Rational[7,4]]
ContinuedFraction[3.141592653589793,8]     FromContinuedFraction[{3,7,15,1}]     SmallPrime[i]    ContinuedFraction[0.696329,6]
```

## 結果の読み方
各行が `式=結果`。`R:cf`・`R:fcf` は連分数。

## gs と pysad の比較
- 連分数 $\pi=[3;7,15,1,292,1,1,1]$ は gs・pysad とも `{3,7,15,1,292,1,1,1}` (既知の値と一致)。$[3;7,15,1]=355/113=3.14159292\ldots$ (`Rational[355,113]`)。0.696329 の連分数 `{0,1,2,3,2,2}` も一致。
- `FromRational[Rational[7,4]]=1.75`、`Rational[2,-4]=Rational[-2,4]`、`Rational[5,1]=5`、`SmallPrime[1..10]` = 2,3,5,…,29 は一致。
- **○ 演算も一致 (2026-11-17 に修正)**: `Rational[1,3]+Rational[1,6]` は gs・pysad とも `Rational[3,6]` (約分は `ReducedRational` が $p^2>b$ で打ち切るため、和は約分されない。gs も同じ)、`Rational[1,3]*Rational[3,4]` は `Rational[1,4]`。以前の pysad は評価せず式のまま残した。`rational.n` の `Plus[Rational[…],…]^:=…` などの UpValue が、組み込みの Plus・Times・Minus・Divide・Subtract の定義として捨てられていたため。`-a`・`1/a`・`a-b`・`a/b`・`0.5+a` も gs と一致する (13 式、試験 `test_rational_arithmetic_uses_package_upvalues_like_gs`)。
- `Numerator[Rational[6,8]]` は 6 (約分されない。`Rational[6,8]` のまま)。

## 未確認・限界
`ReducedRational` の直接の呼び出し、大きな分母、負の数の連分数は試していない。
