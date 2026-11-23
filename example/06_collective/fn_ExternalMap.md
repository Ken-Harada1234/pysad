# fn_ExternalMap
種別: 関数 (ユーザー定義の写像) / フォルダ: 06_collective / ラティス: `00_lattice/map_ring.sad` (新規。`FODO_ring.sad` の ORIGIN の直後に `MAP MP` を置いた `RINGM`)。3 GeV / フラグ: `NORFSW NORAD NOFLUC`

## 何を見る例か
`MAP` 要素は、SADScript の関数 `ExternalMap` で写像を与えるための要素 (定義がなければ恒等で長さ 0)。追跡 (`"TRACK"`) と光学 (`"OPTICS"`) の両方に、薄い四極 (集束 $px\to px-kx$、非集束 $py\to py+ky$) を与え、チューンが解析式どおり変わることを確かめる。

## スクリプト
```
MAP MP = (L=0);   LINE RINGM = (ORIGIN MP RF 24*CELL);
pm = LINE["POSITION", "MP"];        ! 2
k = 0.05;
ExternalMap["OPTICS", n_, cod_] := {cod, {{1,0,0,0,0,0},{-k,1,0,0,0,0},{0,0,1,0,0,0},{0,0,k,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1}}};
ExternalMap["TRACK", n_, nt_, x_] := {x[[1]], x[[2]] - k*x[[1]], x[[3]], x[[4]] + k*x[[3]], x[[5]], x[[6]], x[[7]]};
CALC;                                ! OPTICS の写像で Twiss
```
- `"OPTICS"`: 引数の軌道 `cod` に対し `{cod1, trans}` を返す (6×6 の行列 $(x,p_x,y,p_y,z,\delta)$)。
- `"TRACK"`: `x_` は (7, 粒子数) の表 $(x,p_x(1+\delta),y,p_y(1+\delta),z,\delta,\mathrm{flag})$。同じ形の表を返せば座標が置き換わる。
- 位置 `n` は `LINE["POSITION","MP"]` (= 2) と比べて使い分けてもよい (この例は全 MAP に同じ写像)。

## 結果の読み方
| 行 | gs = pysad |
|---|---|
| ring only の NX, NY | 6.69632894、4.90605604 |
| with MAP (CALC) NX, NY | 6.71594967、4.89631791 |
| with MAP (追跡のチューン小数部) x, y | 0.71594964、0.89631931 |
| 1 要素だけの追跡 | $x,p_x,y,p_y=10^{-6},\,9.5\times10^{-7},\,10^{-6},\,1.05\times10^{-6}$ |

## 物理・解析式
薄いレンズ $L_k=\begin{pmatrix}1&0\\-k&1\end{pmatrix}$ を 1 周の行列 $M$ に掛けた $L_kM$ から、$\cos\mu'=\cos\mu-\frac{k\beta}{2}\sin\mu$ ($\beta$ は MP の位置のリングのβ。垂直は $k\to-k$)。1 要素の追跡は $p_x=p_{x0}-kx_0=10^{-6}-0.05\times10^{-6}=9.5\times10^{-7}$、$p_y=p_{y0}+ky_0=1.05\times10^{-6}$ ($x_0=y_0=p_{x0}=p_{y0}=10^{-6}$)。

| 量 | 式 (行列) | gs | pysad |
|---|---|---|---|
| $\nu_x$ の小数部 | 0.7159497 | 0.71594967 (CALC)、0.71594964 (追跡) | 同左 |
| $\nu_y$ の小数部 | 0.8963179 | 0.89631791 (CALC)、0.89631931 (追跡) | 同左 |
| MAP の $p_x$ | $9.5\times10^{-7}$ | $9.5\times10^{-7}$ | 同左 |
| MAP の $p_y$ | $1.05\times10^{-6}$ | $1.05\times10^{-6}$ | 同左 |

($\nu_y$ の追跡値 0.8963193 の式との差 $1.4\times10^{-6}$ は、小振幅 $10^{-6}$ の 3 周の座標から求めた $\cos\mu$ の有限振幅の誤差。CALC は 8 桁で一致。)

## gs と pysad の比較
○ 一致 (全 4 行が印字桁まで同じ)。

## 未確認・限界
- 6×6 の恒等に近い単純な写像 (薄い四極) だけ。`"EMIT"` (Emittance[] 用) と `"GEO"` は見ていない。落ちた粒子の復活、`ExternalMap` が受け取る `nt` は見ていない。
- docs/19 のとおり、MAIN の `TRACK`・`DynamicApertureSurvey`・Touschek の経路は MAP を恒等として通す (ExternalMap を呼ばない)。
