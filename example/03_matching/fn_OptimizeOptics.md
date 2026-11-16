# fn_OptimizeOptics (OptimizeOptics)
種別: SADScript 関数。フォルダ: 03_matching。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`OptimizeOptics[オプション]`: `GO` のように条件 (`FIT`) と変数 (`FREE`) を使い、**ダウンヒルシンプレックス法**で $1+\text{MatchingResidual}$ (または `OptimizeFunction->` の任意の関数) を最小化する。`GO` と違い、微分を使わないのでどんなキーワード (`L` など) も変数にできる。戻り値は最終のシンプレックス。変数は最小値の点に設定される。
`Tolerance->`、`SimplexSize->` (初期の相対的な大きさ) などの DownhillSimplex のオプションが使える。

## スクリプト
```
FREE QF QD;   FIT $$$;   NX 6.7 NY 4.9;
r = OptimizeOptics[Tolerance->1E-12, SimplexSize->0.1];
Print["OptimizeOptics: K1 QF QD = ", LINE["K1","QF.1"], " ", LINE["K1","QD.1"]];
CALC;     Print["NX NY = ", ...];
```

## 結果の読み方
`DownhillSimplex {limit, Option}: ...` のあと、シンプレックスの 3 頂点の関数値 (1+残差) が 1.00198… → 1.000000 に下がっていく行が並ぶ (gs の出力。pysad も同じ形)。最後の K1 と、`CALC` 後のチューン。

## 物理・解析式
目標 (6.7, 4.9) の厳密な解は 2x2 行列の模型で (0.700221, -0.599636) (FFS_FIT の解と同じ問題。この模型は FFS_FIT、fn_FitValue_... でも使用)。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| K1(QF) | 0.700221 | 0.700221 | 0.700225 |
| K1(QD) | -0.599636 | -0.599636 | -0.599635 |
| NX | 6.7 | 6.70001 | 6.70006 |
| NY | 4.9 | 4.9 | 4.89997 |

(`GO` の解 (0.700221, -0.599636) との差は 1e-5 以下 = シンプレックス法の収束の許容)

## gs と pysad の比較
最終値は一致 (K1 で 5e-6、チューンで 8e-6 以内)。**シンプレックスの反復の表示が違う**: gs は `{`,System`},System`}` や `DownhillSimplex` の限界・オプションの行 (デバッグ出力) と、最小値に達したあとの 1.000000 の行を 126 行ほど出す。pysad は反復を 17 行で打ち切って終わる (最終の `Residual = 1.9E-07` は pysad で、gs は `Matched. ( 3.55E-09 )`: gs のほうがやや深く収束)。

## 未確認・限界
- `L` など `GO` が扱えないキーワードを変数にする例 (SADHelp の例) は未確認。`OptimizeFunction->` に自作の関数を渡す例、`VariableRange` との組み合わせも未確認。
