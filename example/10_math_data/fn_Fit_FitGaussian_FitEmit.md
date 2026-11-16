# fn_Fit_FitGaussian_FitEmit (Fit・FitGaussian・FitEmit)

種別: 関数 (データ当てはめ)・フォルダ 10_math_data・ラティスなし (乱数なし、決定的なデータ)・`FFS;` の後に直接書く

## 何を見る例か
非線形最小二乗 `Fit`、1 次元データのガウス分布当てはめ `FitGaussian`、位相空間の粒子分布から $\alpha,\beta,\epsilon$ を求める `FitEmit` の 3 つ。
`FitGaussian`・`FitEmit` は `Fit` の仕組み (`Packages/Fit.n`、K64 の tffit1) の上に作られている。データは乱数を使わず、式で作った決定的な点を使う。

## スクリプト
```
dat=Table[{i,2*i+1+0.1*Sin[i*3]},{i,1,10}];
r=Fit[dat,a+b*x,x,{a,1},{b,1}];   (* {a->..,b->..,ChiSquare->..,GoodnessOfFit->..,ConfidenceInterval->..,CovarianceMatrix->..,DesignMatrix->..} *)
g=FitGaussian[d];                  (* {sigma, mean, {sigma_conf, mean_conf}, chisq} *)
e=FitEmit[x,px];                   (* {{xmean,pxmean,alpha,beta,emittance},{同じ並びの信頼区間}} *)
```
- (1) $y=2x+1+0.1\sin 3x$ (10 点) を $a+bx$ で当てはめる。`Fit` の引数は `Fit[data, expr, var, {par,初期値}, ...]`。初期値は必須。
- (2) `FitGaussian` の入力は、正規分布の分位点 $\sqrt2\,\mathrm{erf}^{-1}(2(i-0.5)/n-1)$ ($n=400$) に $\sigma=0.8$、平均 0.3、揺らぎ $0.05\sin 1.3i$ を足した 400 個の数。
- (3) `FitEmit` の入力は、(2) と同じ形の分位点を 2 組 (添字を `Mod[i*37,n]` で入れ替えて独立に近くした) 作り、$x=\sqrt{\epsilon\beta}\,u$、$p_x=\sqrt{\epsilon/\beta}\,(v-\alpha u)$ ($\epsilon=10^{-6}$、$\beta=5$、$\alpha=-1$) とした 300 粒子。

## 結果の読み方
- `fit=` の行が `Fit` の戻り値。`a->`、`b->` が当てはめ値、`ChiSquare` は残差の二乗和 (`dy` を与えないとき)、`GoodnessOfFit`、`ConfidenceInterval`、`CovarianceMatrix`、`DesignMatrix` が続く。
- `fg=` の行は `{sigma, mean, {sigma_conf, mean_conf}, chisq}`。
- `fe=` の行は 2 つのリスト。1 つ目が `{xmean, pxmean, alpha, beta, emittance}`、2 つ目が各量の信頼区間。
- `R:` で始まる行は、下の照合表のために数値だけを出したもの。

## 物理・解析式
(1) 線形最小二乗の解は $\hat\theta=(A^TA)^{-1}A^Ty$ ($A$ は `DesignMatrix`)。`ChiSquare` $=\sum_i(y_i-f_i)^2$、共分散行列は $\frac{\mathrm{ChiSquare}}{N-n_p}(A^TA)^{-1}$ ($N=10$、$n_p=2$。表示の CovarianceMatrix と一致することを確認した)。
`GoodnessOfFit` は、誤差を与えないとき分散の推定で規格化した $\chi^2=N-n_p=8$ を使い、`GammaRegularized[(N-n_p)/2, 8/2]` = `GammaRegularized[4,4]` = 0.4334701204 になる (gs・pysad とも .4334701204)。

| 項目 | 最小二乗 (numpy) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| $a$ | 1.0132993837 | 1.0132993837 | 1.0132993864 | 2.3e-12 | 2.7e-09 |
| $b$ | 1.99673824342 | 1.9967382434 | 1.9967382425 | 1.0e-11 | 4.6e-10 |
| ChiSquare (残差二乗和) | 0.0484421798622 | .0484421799 | .0484421799 | 7.8e-10 | 7.8e-10 |

(2) `FitGaussian` は分布の累積分布を逆誤差関数の直線で当てはめるので、$\sigma$ は標本の標準偏差とは少し違う値になる (真値 0.8、データに足した揺らぎのため 0.8032)。

| 項目 | 標本の統計 (np.std, np.mean) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| $\sigma$ (真値 0.8 に揺らぎ) | 0.798983946906 | .8032255747 | .8032255747 | 5.3e-03 | 5.3e-03 |
| 平均 (真値 0.3) | 0.30001439618 | .300007049 | .300007049 | 2.4e-05 | 2.4e-05 |

(3) 参照値は同じ粒子の 2 次モーメント (平均を引いた分散) から $\epsilon=\sqrt{\langle x^2\rangle\langle p^2\rangle-\langle xp\rangle^2}$、$\beta=\langle x^2\rangle/\epsilon$、$\alpha=-\langle xp\rangle/\epsilon$ とした値。
`FitEmit` は 1 次元の `FitGaussian` で分布の幅を求めるので、モーメント法とは一致しない。作った分布が理想の $\beta=5,\alpha=-1,\epsilon=10^{-6}$ から外れるのは、添字の入れ替えによる 2 変数の相関が残るため。

| 項目 | モーメント法 (下式) | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| $\alpha$ | -0.968499544999 | -1.0081380609 | -1.0081380613 | 4.1e-02 | 4.1e-02 |
| $\beta$ [m] | 5.00256138134 | 5.0646593893 | 5.0646593909 | 1.2e-02 | 1.2e-02 |
| $\epsilon$ [m] | 9.952015026e-07 | 8.778359965211E-7 | 8.778359962437E-7 | 1.2e-01 | 1.2e-01 |

## gs と pysad の比較
- `Fit`: ChiSquare、GoodnessOfFit、ConfidenceInterval、CovarianceMatrix は 9 桁以上一致。a と b は 3e-9 の差 (pysad は常に数値微分の経路、K64 は解析的な導関数を使う)。
  DesignMatrix の中身が違う: gs は {1,1},{1,2},... と正確な値だが、pysad は 1.0000000028、0.9999999984 のように数値微分の誤差の混じった値を返す。
- `FitGaussian`: 全桁一致。
- `FitEmit`: alpha、beta、emittance が 1e-9 以内で一致。
- 判定: ○ (値は一致。DesignMatrix の表示だけ pysad が数値微分の値)。

## 未確認・限界
- ConfidenceInterval の係数は式で確かめていない (標準誤差 0.0532 に対し 0.0805、比 1.515)。
- `Fit` の `D->False`・`Cutoff`・`MaxIterations`、$dy$ つきのデータ、非線形の式は例にしていない。
- `FitEmit` の信頼区間は表に入れていない。モーメント法との差の説明 (FitGaussian の仕様) は式の推定で、`FitEmit` のソースまでは確かめていない。
