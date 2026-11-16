# var_CONVERGENCE (特殊変数: CONVERGENCE・MAXITERATION・ExponentOfResidual・OffMomentumWeight・DP・DP0)
種別: 特殊変数。フォルダ: 03_matching。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
FFS の特殊変数 (代入で値を読み書きする)。マッチングに関係するもの:

| 変数 | 意味 | 既定 |
|---|---|---|
| `CONVERGENCE` | `GO` の収束判定 (残差がこれより小さければ収束) | 1e-9 |
| `MAXITERATION` (コマンド) | GO の反復上限 | 40 |
| `ExponentOfResidual` | 残差を作るときの指数 | 2 |
| `OffMomentumWeight` | オフモメンタム点の重み | 1 |
| `DP` | オフモメンタムの運動量幅 | 0.01 |
| `DP0` | CALC の設計運動量のずれ | 0 |

## スクリプト
```
Print["CONVERGENCE = ", CONVERGENCE, "  ExponentOfResidual = ", ExponentOfResidual, "  OffMomentumWeight = ", OffMomentumWeight, "  DP = ", DP, "  DP0 = ", DP0];
FREE QF QD;   FIT $$$;   NX 6.7 NY 4.9;
CONVERGENCE = 1E-3;   GO;      ...      CONVERGENCE = 1E-9;   GO;
MAXITERATION 2;   GO;
DP0 = 0.001;   CALC;
```

## 結果の読み方
- 既定値: CONVERGENCE = 1e-9、ExponentOfResidual = 2、OffMomentumWeight = 1、DP = 0.01、DP0 = 0。
- `CONVERGENCE = 1E-3`: `GO` は **1 回も反復せず** `Matched. ( 1.9799E-03 )` と出して止まる (K1 は元のまま (0.7, -0.6)、チューン 6.696329 / 4.906056)。残差 1.98e-3 が 1e-3 より大きいのに Matched と判定されるので、判定は表示される残差そのものとは別の量 (たぶん残差の平方など。`ExponentOfResidual` の定義) と比べている。正確な式は未確認。
- `CONVERGENCE = 1E-9` (既定): 収束して K1 = (0.700221, -0.599636)、NX,NY = (6.7, 4.9)。
- `MAXITERATION 2`: 結果は上と同じ (この問題は 2〜3 回の Newton 反復で収束するので上限が効かない。docs/03 には、この gs ビルドでは `MAXITERATION` が GO の反復数に効かないと書かれている)。
- `DP0 = 0.001`: チューンが (6.6925, 4.89345) に下がる (運動量偏差 +0.1 % のぶん)。

## 物理・解析式
`DP0` を変えたときのチューンの変化は chromaticity $\xi$ で: $\Delta\nu=\xi\,\delta$。FFS_CHROMATICITY の有限差分の $\xi$ を (K2=0、マッチ前) (-7.512, -6.565) とすると、マッチ後の K1 のリングに近い値で、δ=0.001 → $\Delta\nu\approx$ (-0.00751, -0.00656)。
| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| NX(DP0=0.001) − 6.7 | -0.00751214 | -0.007496 | -0.007496 |
| NY(DP0=0.001) − 4.9 | -0.00656452 | -0.006551 | -0.006551 |

(式の値は 1 次近似、マッチ後のリングの ξ はマッチ前とわずかに違うため、1e-4 程度の差は想定内。)

## gs と pysad の比較
全桁一致 (K1、チューン、`Matched` の残差とも)。

## 未確認・限界
- `ExponentOfResidual`、`OffMomentumWeight` を変えた GO は未確認 (オフモメンタム条件がある FFS_FITPOINTS の設定で変えると効くはず)。
- `MAXITERATION` の効き (反復が上限で止まるケース) は、反復が多く必要な問題で未確認。`TRPT`・`SIGE` など、マッチング以外の特殊変数は 04 以降の担当。
