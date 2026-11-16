# var_Constants (定数 Pi Degree GoldenRatio I INF Infinity NaN SpeedOfLight)

種別: 変数 (SADHelp functions の Constants)・フォルダ 10_math_data・ラティスなし

## 何を見る例か
`Degree GoldenRatio I INF* Infinity NaN* Pi SpeedOfLight` の値と、`INF`・`NaN` の算術・比較での振る舞い。

## スクリプト
```
Pi  Degree  GoldenRatio  SpeedOfLight  I*I  Infinity  -Infinity  INF  NaN  1/Infinity  NaN==NaN  NaNQ[NaN]  Max[1,Infinity]  Exp[1000.]
```

## 結果の読み方
`R:` の行は数値定数。他の行は `名前=値` (`Infinity=INF` は、`Infinity` が `INF` と同じものとして表示されることを示す)。

## 物理・解析式
$\mathrm{Degree}=\pi/180$、$\varphi=\frac{1+\sqrt5}2$ は $\varphi^2-\varphi-1=0$ を満たす。`SpeedOfLight` は定義値 299792458 m/s。

| 項目 | 式の値 | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| `Pi` | 3.14159265359 | 3.141592653589793 | 3.141592653589793 | 0.0e+00 | 0.0e+00 |
| `Degree` $=\pi/180$ | 0.0174532925199 | .017453292519943 | .017453292519943 | 1.7e-14 | 1.7e-14 |
| `GoldenRatio` $=(1+\sqrt5)/2$ | 1.61803398875 | 1.618033988749895 | 1.618033988749895 | 0.0e+00 | 0.0e+00 |
| `SpeedOfLight` | 299792458 | 299792458 | 299792458 | 0.0e+00 | 0.0e+00 |
| $\varphi^2-\varphi-1$ | 0 | 0 | 0 | 0.0e+00 | 0.0e+00 |

## gs と pysad の比較
- 数値定数は 15 桁一致。`INF`・`-INF`・`1/Infinity`=0・`Max[1,Infinity]`・`Exp[1000.]`=INF は一致。
- **`I*I`**: gs は `-1` (実数)、pysad は `(-1+0 I)` (複素数のまま)。表示が違う (△)。
- **`NaN==NaN`**: gs は 1 (真)、pysad は 0 (偽)。IEEE の比較なら偽。`NaNQ[NaN]` は両方 1。(△。K64 の挙動を確かめたが、pysad は IEEE のまま。)

## 未確認・限界
`NaN` を含む `Sort`・`Max`、`INF-INF`、`0*INF` は試していない。`I` の冪や `Sqrt[-1]` は fn_Elementary を参照。
