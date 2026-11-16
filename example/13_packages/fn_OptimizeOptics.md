# fn_OptimizeOptics (optimize.n の OptimizeOptics)

種別: パッケージ (Packages/optimize.n)・フォルダ 13_packages・ラティス FODO_ring (24 セル、3 GeV、周長 105.6 m、Qx=6.696、Qy=4.906)

## 何を見る例か
`FREE QF QD; FIT NX 6.75 NY 4.85;` で QF・QD の $K_1$ を変数、チューンを目標に決め、`OptimizeOptics[]` (DownhillSimplex を使うパッケージ) で合わせる。`FIT` の数値はチューンの整数部を含む値 (サイクル数)。

## スクリプト
```
FREE QF QD;
FIT NX 6.75 NY 4.85;
r=OptimizeOptics[Output->0,MaxIteration->60];
CALC;    Twiss["NX","$$$"]/(2Pi)
```

## 結果の読み方
`R:nx0`・`R:ny0` は最適化前のチューン (6.6963、4.9061)、`variables=2` は `FFS["VAR"]` の変数の数、`R:kqf`・`R:kqd` は結果の $K_1$、`R:nx1`・`R:ny1` は最適化後のチューン。

## gs と pysad の比較
| 量 | 初期 | 目標 | gs | pysad |
|---|---|---|---|---|
| $\nu_x$ | 6.696329 | 6.75 | 6.7499673 | 6.7499673 |
| $\nu_y$ | 4.906056 | 4.85 | 4.8504079 | 4.8504079 |
| $K_1$(QF) | 0.7 | - | 0.7036493 | 0.7036493 |
| $K_1$(QD) | -0.6 | - | -0.5970513 | -0.5970513 |

- gs と pysad は全桁一致 (○)。目標との残差は $3\times10^{-5}$ ($\nu_x$)、$4\times10^{-4}$ ($\nu_y$) 。止まった理由 (`Tolerance` の判定か `MaxIteration->60` か) は見ていない。
- 出力: gs・pysad とも `Output->0` で反復を出さない (2026-11-17 に修正。`OptimizeOptics` は選択肢を `Join[{opt},{VariableRange->limit}]` の規則のリストで `DownhillSimplex` に渡すので、pysad は入れ子のリストを読む。`MaxIteration` も効く)。ただし gs は最初に `{{`,System`},System`}` (BeginPackage の戻り値) を 1 行出し、pysad は `{$ContextPath,"`"}` を出す (文脈を持たないため。結果には影響しない)。

## 未確認・限界
`VariableRange`・`Weight`・他の FIT 項 (`BX`・`EX` など) は試していない。収束の判定の細かい値は見ていない。
