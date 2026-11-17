# fn_Optics_TransferMatrix (Optics.n の TransferMatrix・BeamMatrix・OrbitGeo)

種別: パッケージ (Packages/Optics.n)・フォルダ 13_packages・ラティス FODO_ring (24 セル、3 GeV、周長 105.6 m、Qx=6.696、Qy=4.906)

## 何を見る例か
`TransferMatrix[i,j]` (要素 i から j への転送行列。`$$$` は終端)、`BeamMatrix[i]` (ビーム行列 $\Sigma$)、`OrbitGeo[i]` (軌道の幾何座標)。ここでは 1 周の行列 `TransferMatrix[1,$$$]` の水平 2×2 部分を調べる。

## スクリプト
```
USE RING; CELL; CALC;
t=TransferMatrix[1,$$$];          (* 1 周の 4x4 (RF なしの CALC なので 4 次元) *)
t[[1,1]]+t[[2,2]]    Cos[Twiss["NX","$$$"]]    t[[1,1]]*t[[2,2]]-t[[1,2]]*t[[2,1]]
BeamMatrix[1]    OrbitGeo[5]
```

## 結果の読み方
`dims={4,4}`、`R:m11`〜`m22` は水平の 2×2、`R:trace_x` は対角和、`R:cosmu` は $\cos(\mu_x)$ ($\mu_x$ = `Twiss["NX","$$$"]`、ラジアン)、`R:det_x` は行列式。`BeamMatrix[1] diag=` は $\sigma_{ii}$、`OrbitGeo[5]={.3,0,0}` は 5 番目の要素 (偏向磁石の入口?) の位置 $(x,y,z)$ (定義未確認)。

## gs と pysad の比較
1 周の行列のトレースは $2\cos\mu_x$、行列式は 1 (シンプレクティック)。

| 項目 | 式 | gs | pysad | gs の相対差 | pysad の相対差 (参照値が 0 の行は絶対差) |
|---|---|---|---|---|---|
| $\mathrm{tr}\,M_x=2\cos\mu_x$ ($\mu_x=2\pi\nu_x$) | -0.661739730794 | -.6617397308 | -.6617397308 | 9.0e-12 | 9.0e-12 |
| $\det M_x=1$ | 1 | 1 | 1 | 0.0e+00 | 0.0e+00 |

(gs・pysad とも `trace_x=-0.6617397308`、`cosmu=-0.3308698654` で、$2\cos\mu_x=-0.66174$ と一致。)

- 4×4 の転送行列の全要素・`OrbitGeo` は gs と pysad で表示 10 桁一致 (○)。
- `BeamMatrix[1]` の対角: 水平・垂直の $\sigma_{11},\sigma_{22}$ は一致。結合のない成分 ($\sigma_{33},\sigma_{44}$) は gs が 1e-37 の雑音、pysad は 0 (丸め誤差の差で、docs/34 と同じ)。

## 未確認・限界
`NormalCoordinates`・`PeriodicR` は試していない (`PeriodicR` は docs/34 でも呼び出し方が不明)。`OrbitGeo` の戻りの意味 (座標系) は確認していない。この例では RF なし (`RFSW` なし) で 4 次元。
