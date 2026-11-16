# elem_QUAD (四極磁石)
種別: elem。フォルダ: 01_lattice_elements。ラティスは START-QUAD-FIN の輸送路 (INS、初期 BX=BY=10 m、AX=AY=0)。エネルギー 3 GeV (K1 は積分強度なので結果はエネルギーによらない)。

## 何を見る例か
`QUAD Q = (L=0.3 K1=0.7)`: **K1 は積分強度** ($K_1=k_1L$ [1/m]、正で水平収束)。1 台の転送行列と β の変換を、厚いレンズ (L=0.3) と、ほぼ薄いレンズ (L=0.001, K1=0.021) の両方で見る。

## スクリプト
```
QUAD Q = (L=0.3 K1=0.7) QT = (L=0.001 K1=0.021);
MARK START = (BX=10 BY=10 AX=0 AY=0) FIN = ();
LINE T = (START Q FIN);   LINE TT = (START QT FIN);
FFS;
USE T; INS; CALC;
m = TransferMatrix["START","FIN"];
Print["thick Q  M11 M12 M21 M22 = ", m[[1,1]], " ", m[[1,2]], " ", m[[2,1]], " ", m[[2,2]]];
```
`TransferMatrix[from, to]` は 4x4 (CALC4D) の行列 ({x, px, y, py})。全文は `elem_QUAD.sad`。

## 結果の読み方
`thick Q ...` の 2 行が x と y の 2x2 ブロック、`BX FIN`/`BY FIN` が出口の β、`NX FIN`/`NY FIN` が位相進み [rad]。最後の 1 行が薄い四極の $M_{21}$ と $M_{43}$。

## 物理・解析式
$k=K_1/L$、$\phi=\sqrt{k}\,L$ として、水平 (収束) と垂直 (発散):
$$M_x=\begin{pmatrix}\cos\phi&\sin\phi/\sqrt k\\-\sqrt k\sin\phi&\cos\phi\end{pmatrix},\quad
M_y=\begin{pmatrix}\cosh\phi&\sinh\phi/\sqrt k\\\sqrt k\sinh\phi&\cosh\phi\end{pmatrix},\quad
\beta_1=M_{11}^2\beta_0-2M_{11}M_{12}\alpha_0+M_{12}^2\gamma_0$$
薄いレンズの極限は $M_{21}=-K_1$ ($f=1/K_1$)。

| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| x: M11 = cos φ | 0.896825 | 0.896825 | 0.896825 |
| x: M12 = sin φ/√k | 0.28961 | 0.28961 | 0.28961 |
| x: M21 = -√k sin φ | -0.675756 | -0.675756 | -0.675756 |
| y: M33 = cosh φ | 1.10685 | 1.10685 | 1.10685 |
| y: M34 = sinh φ/√k | 0.310611 | 0.310611 | 0.310611 |
| y: M43 = √k sinh φ | 0.724759 | 0.724759 | 0.724759 |
| BX(FIN) | 8.05133 | 8.05133 | 8.05133 |
| BY(FIN) | 12.2608 | 12.2608 | 12.2608 |
| 薄レンズ M21 = -K1 | -0.021 | -0.021 | -0.021 |
| 薄レンズ M43 = +K1 | 0.021 | 43 | 43 |


## gs と pysad の比較
全桁一致。式との差は 1e-6 以下 (表示桁の範囲)。薄いレンズ (L=0.001) の $M_{21}$ は厳密な -K1 (表の値の差は L による 2 次の項)。

## 未確認・限界
- ずれ `DX DY`、回転 `ROTATE` (スキュー四極)、フリンジ (`DISFRIN`、`FRINGE`) は扱っていない。
- 非線形 (K1 に依らない高次) は追跡の例 (05) の担当。
