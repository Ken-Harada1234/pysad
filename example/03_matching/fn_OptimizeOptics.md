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
**○ 一致 (2026-11-17)。** 最終値は gs と同じ: K1 (QF, QD) = (0.700221, -0.599636)、`Matched. ( 3.5504E-09)`、NX NY = 6.700009、4.900003 (6 桁)。シンプレックスの反復の表も gs と同じ行数 (最小値に達したあとの 1.000000 の行を含む) が出る。互換モードによらず同じ (`PYSAD_COMPAT=k64` でも gs と同じ値を確認)。表示の差が 2 行だけ残る: デバッグ出力 `{{`,System`},System`}` の引用符 (pysad は `"System`"`)、`DownhillSimplex {limit, Option}:` の行を pysad は `...` で省略する。

## 未確認・限界
- `L` など `GO` が扱えないキーワードを変数にする例 (SADHelp の例) は未確認。`OptimizeFunction->` に自作の関数を渡す例、`VariableRange` との組み合わせも未確認。
