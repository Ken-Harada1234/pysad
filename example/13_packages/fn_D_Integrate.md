# fn_D_Integrate (D.n の D と Integrate)

種別: パッケージ (Packages/D.n)・フォルダ 13_packages・ラティスなし

## 何を見る例か
記号微分 `D[f,x]`・`D[f,{x,n}]` (D.n) と、`Integrate` の動き。`D` は式を返し、そこに値を代入すれば数値が出る。

## スクリプト
```
D[x^3,x]   D[Sin[x]*Exp[x],x]   D[Log[x^2+1],x]   D[x^2*y,y]   D[Sin[x]^2,{x,2}]   D[Sqrt[1+x^2],x]
f=D[Sin[x]*Exp[x],x];   f/.x->0.7
Integrate[x^2,x]   Integrate[Sin[x],x]   Integrate[x^2,{x,0,3}]
```

## 結果の読み方
各行が `式=結果`。結果は SAD の表示形 (`(3 x^2)`、`(2 x (1+x^2)^-1)` など。積は括弧つきの `(a b)`、べきは `^`)。

## gs と pysad の比較
微分した式 $\cos x\,e^x+e^x\sin x$ に $x=0.7$ を代入した値と、元の関数の中心差分:

| 項目 | 解析式 $e^x(\sin x+\cos x)$ | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| $x=0.7$ での $\frac{d}{dx}(e^x\sin x)$ (D で作った式に代入) | 2.83749813731 | 2.837498137307049 | 2.837498137307049 | 1.6e-16 | 1.6e-16 |
| 同 (中心差分、刻み $10^{-6}$) | 2.83749813731 | 2.8374981373424646 | 2.8374981373424646 | 1.2e-11 | 1.2e-11 |

- `D` の 6 式の出力 (形) は gs と pysad で完全に一致 (○)。
- **`Integrate` は記号積分をしない**: `Integrate[x^2,x]`、`Integrate[Sin[x],x]`、`Integrate[Exp[2 x],x]` は評価されずそのまま返る (gs も pysad も同じ)。定積分 `Integrate[x^2,{x,0,3}]` も `Integrate$a[x^2,{x,0,3}]` という内部の式のまま。`Integrate` が評価されるのは `SplineData` の積分 (fn_Spline の `Integrate[sp[x],{x,0,Pi}]`) だけ。
- `Integrate[1/x,x]` の表示だけ違う (gs `Integrate[x^-1,x]`、pysad `Integrate[Divide[1,x],x]`、評価されない点は同じ)。

## 未確認・限界
`D` の `Derivative`・`ReduceExpression`・特殊関数の微分、`D` の出力の標準形 (docs/06 に「記号式の標準形は gs と同じではない」とある。59 式で一致を確認した範囲外は未確認) は試していない。
