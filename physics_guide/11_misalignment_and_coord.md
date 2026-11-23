# 11. 取り付け誤差と COORD: 座標の正準変換

**この章で分かること**: 要素の取り付け誤差 (DX・DY・DZ・ROTATE・DROTATE・CHI1・CHI2) と COORD 要素が、力学ではなく**座標の正準変換**であること。(1) 粒子の軌跡を 3 次元の直線とみなし、剛体の運動 (平行移動と回転) を施して、新しい面で測り直す、という 1 つの構成で、COORD の追跡 (`tcoord`)・CALC の写像と 6×6 行列 (`tcoorde`)・偏向磁石のロール (`tbrot`・`tbrote`・`tbgrote`) が出ること。(2) 運動量 $p_x, p_y$・$z$・$\delta$ の扱い。(3) DX・DY・ROTATE の順序と符号。(4) MARK・APERT・MAP は力学の写像でないこと。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) (正準変数、$z$ の符号)、[02_drift.md](02_drift.md) (自由な直線運動)。

**読む順**: §1 (結論) → §3 (剛体運動と面への戻し) → §4 (COORD) → §5 (DX・DY・ROTATE) → §6 (偏向磁石のロール)。

## 1. 結論

- 取り付け誤差は、要素の前後で座標を取り替える操作で、要素の中の力学 (ハミルトニアン) は変えない。要素の入口で「ずれと回転を入れた座標」に移り、要素を通り、出口で元の座標に戻る。
- **一般の剛体運動** ($\vec r \to R(\vec r - \vec d)$、$\vec P \to R\vec P$) は、粒子の軌跡が自由な直線であるところ (要素の端の面) では、正準変換である。位置と運動量を同じ回転で回すのは 3 次元の点変換で、正準だからである。回転すると基準面が傾くので、**直線を新しい面 $z' = 0$ まで延ばして測り直す** (§3)。この「回転 + 面への戻し」の写像は、数値で $|J^\top SJ - S| \lesssim 10^{-12}$ (シンプレクティック)。
- 時間は保たれる: 粒子の位置を「基準粒子が面にいる時刻」のスナップショット $\vec r = (x, y, 0) + z\,\vec u$ ($\vec u = \vec P/p$ は進行方向の単位ベクトル。$z > 0$ は先) で表すと、剛体運動のあと新しい面から測った $z' = z'_{\rm 3D}/u'_z$ が出る (式 (11.4))。
- COORD の DIR = 0 は、ヘルプの行列 $R_3R_2R_1$ と一致する。**DIR ≠ 0 の K64 の行列は $R_1R_2^{(-)}R_3$ ($R_2^{(-)}$ は $R_2$ の $s_2$ の符号を逆にしたもの) で、ヘルプの式 $R_3R_2^{(-)}R_1$ とは積の順序が逆**。K64 の行列は、DIR = 0 の真の逆 $(R_3R_2R_1)^\top$ の $\chi_1$・$\chi_3$ の符号を反転したもの ($y$ と $z$ を反転した座標で書いた逆変換) とも読める (§4.1、付録 B。ヘルプの記述の誤記の可能性が高い)。DIR ≠ 0 は DIR = 0 の逆変換でもない (ヘルプが "NOT the inverse" と書く通り)。真の逆は $R^\top$ と平行移動の戻し。
- **ROTATE**: 入口で $(x, y) \to R(+\theta)(x - \mathrm{DX},\, y - \mathrm{DY})$ ($p_x, p_y$ も同じ $R(\theta)$ で回す)、出口で逆。DX・DY は回転の前に引く (QUAD で差 $10^{-19}$)。
- **MARK・APERT・MAP は座標を変える力学ではない**: MARK は恒等写像 (座標は変わらない)、APERT は通る粒子の座標を変えず、失う粒子の記録だけを作る (追跡のみ)、MAP はユーザー定義の写像 `ExternalMap` を呼ぶだけ (§8)。
- pysad の追跡は、直線の MULT の CHI1・CHI2・DZ を使わない (指定しても結果が変わらない。確認した)。K64 の `tsolrot` にはある (09 章 §8.2 で、$B_z = 0$ の枝が剛体運動 + 面への戻しと $10^{-17}$ で一致することを確認した)。BEND の CHI2・DROTATE (`tbrot`) は実装されていて、本章の構成と $10^{-17}$ で一致し、その回転行列は軸のまわりの回転の積として代数的に導ける (§6)。

## 2. 記号と単位

| 量 | SAD (ヘルプ) | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| 平行移動 | `DX`、`DY`、`DZ` [m] | `dx`、`dy`、`dz` (`tcoord`)、`dx0`・`dy0`・`dz0` (入口の引く量)、`dx1`・`dy1`・`dz1` (出口で足す量) | `coord._rot` の戻り値 |
| 回転角 | `CHI1` (y 軸まわり)、`CHI2` (x 軸まわり)、`CHI3` = `ROTATE` (z 軸まわり) | `chi1`・`chi2`・`chi3`、`cchi1 = cos chi1`・`schi1 = sin chi1` など | `coord.params` |
| 回転行列 $R$ | 3×3 | `r11`〜`r33` | `coord._rot` |
| 向き | `DIR` (0 で順変換) | `dir` (logical) | `direct` |
| 取り付け誤差 (CALC) | `ROTATE`、`DROTATE`、`DX`、`DY` | `tchge` の引数 `dx0`・`dy`・`dz`・`theta`・`dtheta`・`dchi2`、`ent` (入口なら真) | `lattice._tchge_cod`・`_tchge_jac` |
| 偏向磁石の支点まわりの量 | | `alg` (弦の半分の弧長)、`phig` (半分の偏向角)、`dl = alg sin(phig)/phig`、`c0 = cos phig`・`s0 = sin phig` | `track6._tbrot`、`elkern._tbrote` |
| 運動量 | $p_x, p_y$ (正準)、$\delta$ | CALC は `cod(2)`・`cod(4)`・`cod(6)`、`p = 1 + cod(6)`。追跡は $q = p_x/(1+\delta)$ | |

`xsincos(a)` は `(sin a, a - sin a の誤差の小さい形, cos a, cos a - 1)` を返す Fortran の関数 (`tfloor.f`)。コードの `dc` は $\cos a - 1$ (小さい角の桁落ちの回避)。

## 3. 剛体運動と面への戻し

### 3.1 なぜ正準変換か

SAD の変数 $(x, p_x, y, p_y, z, \delta)$ は、基準線に垂直な面 $s = $ 一定を粒子が横切るときの値である (01 章)。要素のずれ・回転は、**粒子ではなく要素 (座標系) を動かす**。要素の座標系で見た粒子の位置は、元の座標に剛体運動 $T$ を施したものである。

自由空間 (要素の端の面の外) では、粒子の軌跡は直線である。位置 $\vec r$ と運動量 $\vec P$ (大きさ $p = |\vec P|$、$p_0$ で規格化) について、直線は

$$
\vec r(\lambda) = \vec r_0 + \lambda\,\vec u,\qquad \vec u = \frac{\vec P}{p}
\tag{11.1}
$$

(11.1) で、剛体運動 $T:\ \vec r \to R(\vec r - \vec d)$、$\vec P \to R\vec P$ ($R$ は回転行列、$\vec d$ は平行移動) を施すと、直線は直線に写り、向きは $R\vec u$ になる。$(\vec r, \vec P)$ の変換 $(R\,\cdot,\, R\,\cdot)$ は、点変換 $\vec r' = R(\vec r - \vec d)$ とその正準な拡張 $\vec P' = (R^{-1})^\top \vec P = R\vec P$ ($R$ は直交) なので、**正準**である。母関数で書くと $F_2(\vec r, \vec P') = \vec P'\cdot R(\vec r - \vec d)$:

$$
\vec P = \frac{\partial F_2}{\partial \vec r} = R^\top \vec P',\qquad \vec r' = \frac{\partial F_2}{\partial \vec P'} = R(\vec r - \vec d)
\tag{11.2}
$$

($\vec P' = R\vec P$、$R^\top R = 1$ で整合する。)

### 3.2 面への戻し

回転すると、要素の座標系での基準面 $z' = 0$ は元の面 $s = 0$ と平行でない。粒子は元の面 $s = 0$ を $(x, y, 0)$ で横切った。新しい座標系では、$z' = 0$ の面を横切る点と時刻で測り直す必要がある。軌跡が直線 (11.1) なので、剛体運動のあと、直線を $z' = 0$ まで延ばすだけでよい。これは長さが不定の**ドリフト** (自由な直線運動) で、運動量は変わらない。

### 3.3 $z$ の扱い: スナップショット

SAD の $z$ は、基準粒子が面にいる時刻 $t_0$ に対して、$z = -v(t - t_0)$ ($v$ は粒子の速さ、先が正)。粒子が面を横切る時刻が $t$ なので、時刻 $t_0$ には粒子は面から進行方向に距離 $z$ だけ先にいる ($z > 0$ なら先に通過している)。時刻 $t_0$ の位置は

$$
\vec r = (x,\ y,\ 0) + z\,\vec u,\qquad
\vec u = \frac{(p_x,\ p_y,\ p_z)}{p},\quad p = 1 + \delta,\quad p_z = \sqrt{p^2 - p_x^2 - p_y^2}
\tag{11.3}
$$

剛体運動は、同じ時刻の位置を同時に動かす (静止した座標系の取り替え)。時刻 $t_0$ の新しい座標系での位置と運動量を $\vec r' = R(\vec r - \vec d)$、$\vec P' = R\vec P$ とする。新しい面 $z' = 0$ から測ると、粒子は進行方向 $\vec u' = \vec P'/p$ に距離 $z'$ だけ先にいる:

$$
z' = \frac{r'_z}{u'_z} = \frac{p\,r'_z}{P'_z},\qquad
x' = r'_x - u'_x\,z',\quad y' = r'_y - u'_y\,z',\quad p'_x = P'_x,\ p'_y = P'_y,\ \delta' = \delta
\tag{11.4}
$$

(11.4) の $x', y'$ は、粒子が面 $z' = 0$ を横切る点 ($r'_z - u'_z z' = 0$ になる) である。$\delta$ は大きさ $p$ だけで決まり、回転では変わらない。ここまでが `tcoord`・`tcoorde` の本体で、`z(i)=zf/pzf`、`x(i)=xf-px(i)*z(i)+dx1`、CALC の `cod(5)=p/pzf*zf` と同じ式。追跡の $q = p_x/p$ では $p$ が消えて `z = zf/pzf` になる。

**K64 の「$z$ に足す」形との同値**。K64 の `tsolrot`・`trots` は $z$ を $z \to z + z_f/p_{zf}$ のように**足す**形で書く。(11.3) の $z\vec u$ の項は、剛体運動のあとも同じ向き $\vec u'$ を向くので、$r' = R\bigl[(x, y, 0) + z\vec u\bigr] = R(x, y, 0) + z\,\vec u'$。したがって (11.4) の $z' = r'_z/u'_z$ は

$$
z' = \frac{[R(x, y, 0)]_z + z\,u'_z}{u'_z} = z + \frac{[R(x, y, 0)]_z}{u'_z} = z + \frac{z_f}{q_{zf}},\qquad
x' = [R(x, y, 0)]_x + z\,u'_x - u'_x\,z' = x_f - \frac{q_{xf}}{q_{zf}}z_f
$$

($z_f = [R(x, y, 0)]_z$。$z\,u'_x$ が打ち消し合う。) つまり、スナップショットの形 ((11.3)(11.4)) と、K64 の「$z$ に $z_f/q_{zf}$ を足し、$x$ から $q_x z_f/q_{zf}$ を引く」形は同じ写像で、09 章 §8.1 の `trots` の式になる。

**確認** ([mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) (b)): 7 組のずれ・角度 × 5 粒子 (角度は最大 0.5 rad) で、(11.3)(11.4) から numpy で独立に作った写像と、pysad の追跡 (`tcoord`) および CALC の COD の更新 (`tcoorde`) の差は $5.6 \times 10^{-17}$。

### 3.4 6×6 の行列 (`tcoorde`): 3 次元の変数を挟む連鎖律

CALC は $(x, p_x, y, p_y, z, \delta)$ の写像の 1 次の係数 $T = \partial(\mathrm{out})/\partial(\mathrm{in})$ が要る。写像は 3 段:

1. $(x, p_x, y, p_y, z, \delta) \to (\vec r, \vec P)$ (11.3)。
2. $(\vec r, \vec P) \to (R(\vec r - \vec d),\, R\vec P)$ (線形)。
3. $(\vec r', \vec P') \to (x', p_x', y', p_y', z', \delta')$ (11.4)。

**段 1** の偏微分 (変数は $x, p_x, p_y, z, \delta$、$p = 1+\delta$、$p_z = \sqrt{p^2 - p_x^2 - p_y^2}$)。位置の $x$ 成分 $r_x = x + p_x z/p$ は

$$
\frac{\partial r_x}{\partial p_x} = \frac{z}{p},\qquad
\frac{\partial r_x}{\partial z} = \frac{p_x}{p},\qquad
\frac{\partial r_x}{\partial \delta} = -\frac{p_x z}{p^2}
$$

($r_y$ も同じ形)。$r_z = p_z z/p$ の偏微分には $\partial p_z/\partial p_x = -p_x/p_z$、$\partial p_z/\partial\delta = p/p_z$ を使う:

$$
\frac{\partial r_z}{\partial p_x} = \frac{z}{p}\left(-\frac{p_x}{p_z}\right) = -\frac{p_x z}{p\,p_z},\qquad
\frac{\partial r_z}{\partial z} = \frac{p_z}{p},\qquad
\frac{\partial r_z}{\partial \delta} = z\left(\frac{1}{p_z} - \frac{p_z}{p^2}\right) = \frac{(p_x^2 + p_y^2)\,z}{p^2\,p_z}
$$

(最後の等号: $\dfrac{1}{p_z} - \dfrac{p_z}{p^2} = \dfrac{p^2 - p_z^2}{p^2 p_z}$ と $p^2 - p_z^2 = p_x^2 + p_y^2$。) 運動量の $z$ 成分は $\partial P_z/\partial p_x = -p_x/p_z$、$\partial P_z/\partial\delta = p/p_z$。Fortran の `trans1` (`trans1(1,2) = z/p`、`trans1(1,5) = pxi/p`、`trans1(5,2) = -pxi/p/pzi*z`、`trans1(5,6) = (pxi**2+pyi**2)/p**2/pzi*z`、`trans1(6,2) = -pxi/pzi`、`trans1(6,6) = p/pzi` など) と一致する。ここで `trans1` の第 5 行は $r_z$、第 6 行は $P_z$ を表す (中間の 6 変数は $(r_x, P_x, r_y, P_y, r_z, P_z)$ の並び)。

**段 2** は行列 $R$ を位置の行 (1, 3, 5) と運動量の行 (2, 4, 6) に別々に掛ける。Fortran の `do 110` のループ。

**段 3**: $z' = p\,r'_z/P'_z$、$p = |\vec P'|$ ($=\sqrt{P'^2_x + P'^2_y + P'^2_z}$ を $P'$ の関数とみなす。$\partial p/\partial P'_x = P'_x/p$ など) の偏微分:

$$
\frac{\partial x'}{\partial r'_z} = -\frac{P'_x}{P'_z},\quad
\frac{\partial x'}{\partial P'_x} = -\frac{r'_z}{P'_z},\quad
\frac{\partial x'}{\partial P'_z} = \frac{P'_x r'_z}{P'^2_z},\quad
\frac{\partial z'}{\partial r'_z} = \frac{p}{P'_z},
$$
$$
\frac{\partial z'}{\partial P'_x} = \frac{r'_z}{P'_z}\cdot\frac{P'_x}{p},\qquad
\frac{\partial z'}{\partial P'_z} = r'_z\left(\frac{P'_z/p}{P'_z} - \frac{p}{P'^2_z}\right) = -\frac{(P'^2_x + P'^2_y)\,r'_z}{p\,P'^2_z},\qquad
\frac{\partial \delta'}{\partial P'_{x}} = \frac{P'_x}{p},\ \frac{\partial \delta'}{\partial P'_{z}} = \frac{P'_z}{p}
\tag{11.5}
$$

(最後の $z'$ の $P'_z$ 微分の計算: $\dfrac{1}{p} - \dfrac{p}{P'^2_z} = \dfrac{P'^2_z - p^2}{p P'^2_z}$。) Fortran の `trans2` (`trans2(1,2) = -zf/pzf`、`trans2(1,6) = pxf/pzf**2*zf`、`trans2(5,5) = p/pzf`、`trans2(5,2) = pxf/p/pzf*zf`、`trans2(5,6) = -(pxf**2+pyf**2)/p/pzf**2*zf`、`trans2(6,2) = pxf/p`、`trans2(6,6) = pzf/p`) と一致する。最後に `trans1 * trans2` を全体の行列 `trans` の左から掛ける。

**確認** ((c)): 解析式の行列 (`coord.coord_map` の $T$) と、独立な写像の中心差分 (刻み $10^{-6}$) の最大の差は $5 \times 10^{-12}$。$|J^\top SJ - S|$ は解析式で $1.7 \times 10^{-18}$、数値微分で $5 \times 10^{-12}$ (シンプレクティック)。

**シンプレクティック性の理由**: 6 変数 $(x, p_x, \ldots)$ は、自由空間のハミルトニアン (経路を独立変数に取る) の流れの、ある面での値である。同じ流れの別の面での値への写像は、面の取り替えに関わらずシンプレクティックである (ポアンカレ断面間の写像の一般的な性質)。剛体運動 (11.2) は正準変換なので、その合成も正準。上の数値の確認はこの議論と整合する。数値で確かめたのは (b)(c) の範囲 (角度 $\le 0.5$ rad、$\delta = O(10^{-3})$、$z = O(0.05$ m$)$) であり、一般の証明ではない。

### 3.5 小さな角度の 1 次

(11.3)(11.4) の写像で、ずれ・角度のパラメータ ($\chi_1, \chi_2, \chi_3, \mathrm{DX}, \mathrm{DY}, \mathrm{DZ}$) の 1 次の係数を求める。$R = R_3R_2R_1$ を 1 次に展開すると $R \simeq 1 + \begin{pmatrix}0 & -\chi_3 & \chi_1\\ \chi_3 & 0 & \chi_2\\ -\chi_1 & -\chi_2 & 0\end{pmatrix}$ ($R_1 \simeq 1 + \chi_1\,(\hat x\hat z^\top - \hat z\hat x^\top)$ など)。位置は $\vec r = (x - \mathrm{DX},\ y - \mathrm{DY},\ 0) + z\vec u - (0, 0, \mathrm{DZ})$ とする ($\mathrm{DZ}$ は $r_z$ から引く)。

パラメータの 1 次で、$R \simeq 1 + A$ ($A$ は上の反対称行列)、$\vec r = \vec r_0 - \vec d$ ($\vec r_0 = (x, y, 0) + z\vec u$、$\vec d = (\mathrm{DX}, \mathrm{DY}, \mathrm{DZ})$)、$\vec P = (p_x, p_y, p_z)$ とする。0 次では $z'_0 = p\,r_{0z}/p_z = z$ ($r_{0z} = zp_z/p$)。

1. $\vec r' = (1 + A)(\vec r_0 - \vec d)$ の 1 次: $\delta r'_z = -\mathrm{DZ} + (A\vec r_0)_z = -\mathrm{DZ} - \chi_1r_{0x} - \chi_2r_{0y}$、$\delta r'_x = -\mathrm{DX} + (A\vec r_0)_x = -\mathrm{DX} - \chi_3r_{0y} + \chi_1r_{0z}$。
2. $\vec P' = (1 + A)\vec P$ の 1 次: $\delta P'_x = -\chi_3p_y + \chi_1p_z$、$\delta P'_y = \chi_3p_x + \chi_2p_z$、$\delta P'_z = -\chi_1p_x - \chi_2p_y$。
3. $z' = p\,r'_z/P'_z$ の 1 次: $\delta z' = \dfrac p{p_z}\,\delta r'_z - \dfrac{p\,r_{0z}}{p_z^2}\,\delta P'_z = \dfrac p{p_z}\bigl[-\mathrm{DZ} - \chi_1r_{0x} - \chi_2r_{0y}\bigr] + \dfrac{z}{p_z}(\chi_1p_x + \chi_2p_y)$ ($r_{0z} = zp_z/p$)。$r_{0x} = x + zp_x/p$、$r_{0y} = y + zp_y/p$ を入れると、$z$ を含む項 $-\chi_1zp_x/p_z - \chi_2zp_y/p_z$ が最後の項と打ち消し合い、
   $$z' \simeq z - \frac p{p_z}\,(\mathrm{DZ} + \chi_1x + \chi_2y)$$
4. $x' = r'_x - (P'_x/p)\,z'$ の 1 次: $\delta x' = \delta r'_x - \dfrac{\delta P'_x}{p}z - \dfrac{p_x}{p}\,\delta z'$。第 1 項 $-\mathrm{DX} - \chi_3(y + zp_y/p) + \chi_1zp_z/p$、第 2 項 $-(-\chi_3p_y + \chi_1p_z)z/p$ の $z$ を含む項 ($-\chi_3zp_y/p + \chi_1zp_z/p$ と $+\chi_3zp_y/p - \chi_1zp_z/p$) は打ち消し合う。第 3 項は $+(p_x/p_z)(\mathrm{DZ} + \chi_1x + \chi_2y)$ なので
   $$x' \simeq x - \mathrm{DX} - \chi_3\,y + \frac{p_x}{p_z}\bigl(\mathrm{DZ} + \chi_1x + \chi_2y\bigr)$$

($y'$ は、$x'$ の式で $\mathrm{DX} \to \mathrm{DY}$、$-\chi_3y \to +\chi_3x$、$p_x \to p_y$ としたもの。組 $(\chi_1x + \chi_2y + \mathrm{DZ})$ はそのまま。) この 1 次の係数は sympy で、(11.3)(11.4) の写像をパラメータで偏微分して確かめた ([mis_03_derivations.out.txt](outputs/mis_03_derivations.out.txt) (d)):

| 出力 | 1 次の項 |
|---|---|
| $x'$ | $-\mathrm{DX} - \chi_3\,y + \dfrac{p_x}{p_z}(\chi_1x + \chi_2y + \mathrm{DZ})$ |
| $p_x'$ | $+\chi_1p_z - \chi_3\,p_y$ |
| $y'$ | $-\mathrm{DY} + \chi_3\,x + \dfrac{p_y}{p_z}(\chi_1x + \chi_2y + \mathrm{DZ})$ |
| $p_y'$ | $+\chi_2p_z + \chi_3\,p_x$ |
| $z'$ | $-\dfrac{p}{p_z}(\chi_1x + \chi_2y + \mathrm{DZ})$ |
| $\delta'$ | 0 |

特に:

- **CHI1 = $\chi_1$** (x–z 面の回転): $p_x' \simeq p_x + \chi_1p_z$。設計軌道 ($p_x = 0$) は $p_x' = \chi_1$ (+1.0e-6 rad の入力で +1.0e-6 ([mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) (d)))。粒子は進行方向が $\chi_1$ だけ傾く (設計軌道が要素の軸に対して傾く)。
- **CHI2** も同じで $p_y' \simeq p_y + \chi_2p_z$ (数値: +1.0e-6)。
- **CHI3 (ROTATE)**: 設計軌道 ($x = y = p_x = p_y = 0$) は変わらない (係数がすべて $x, y, p_x, p_y$ に比例)。
- **DX**: $x' = x - \mathrm{DX}$、$p_x$ は変わらない。
- **DZ**: 設計軌道の位置は変わらないが $z' = z - \mathrm{DZ}$ ($p = p_z$ のとき。数値: $-1.0 \times 10^{-2}$ m)。$p_x \ne 0$ のとき $x' = x + (p_x/p_z)\,\mathrm{DZ}$ ($p_x = 10^{-3}$、DZ = 1 cm で $1.0 \times 10^{-5}$ m)。
- 2 次以上は、$R$ の非線形 ($\sin\chi$, $\cos\chi$) と $p_z$ の展開から出る。

## 4. COORD 要素

COORD は、座標を (11.3)(11.4) で取り替えるだけの要素で、長さがない。力学を持たない。

### 4.1 ヘルプの行列と K64

ヘルプの COORD (DIR = 0): $(x, y, z)_1 = R_3R_2R_1\,(x - \mathrm{DX},\ y - \mathrm{DY},\ z)_0$、

$$
R_3 = \begin{pmatrix} c_3 & -s_3 & 0\\ s_3 & c_3 & 0\\ 0&0&1\end{pmatrix},\quad
R_2 = \begin{pmatrix} 1&0&0\\ 0&c_2&s_2\\ 0&-s_2&c_2\end{pmatrix},\quad
R_1 = \begin{pmatrix} c_1&0&s_1\\ 0&1&0\\ -s_1&0&c_1\end{pmatrix}
\tag{11.6}
$$

行列の積を展開すると (sympy で確認)

$$
R_3R_2R_1 = \begin{pmatrix}
c_1c_3 + s_1s_2s_3 & -c_2s_3 & s_1c_3 - c_1s_2s_3\\
-s_1s_2c_3 + c_1s_3 & c_2c_3 & c_1s_2c_3 + s_1s_3\\
-s_1c_2 & -s_2 & c_1c_2
\end{pmatrix}
$$

で、`tcoord.f` の `r11`〜`r33` (`dir = .true.`) と一致する。DX・DY・DZ は `dx0 = dx`・`dy0 = dy`・`dz0 = dz` として、回転の前に引く (`zi = pzi*z - dz0`、`xi = x + pxi*z - dx0`)。ヘルプの式は DZ を書かない (キーワード表にも DZ が無い) が、K64 は DZ を持つ。

DIR ≠ 0 のヘルプの式 (原文で確認。付録 A) は、$(x, y, z)_1 = R_3\,R_2^{(-)}R_1\,(x, y, z)_0 + (\mathrm{DX}, -\mathrm{DY}, 0)$ で、$R_2^{(-)} = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}$ ($R_2$ の $s_2$ の符号を逆にしたもの)。一方、K64 の `dir = .false.` の行列は

$$
\begin{pmatrix}
c_1c_3 + s_1s_2s_3 & s_1s_2c_3 - c_1s_3 & s_1c_2\\
c_2s_3 & c_2c_3 & -s_2\\
-s_1c_3 + c_1s_2s_3 & c_1s_2c_3 + s_1s_3 & c_1c_2
\end{pmatrix},\qquad \text{出口で } (+\mathrm{DX},\ -\mathrm{DY},\ -\mathrm{DZ})\ \text{を足す}
\tag{11.7}
$$

である。**(11.7) はヘルプの式 $R_3R_2^{(-)}R_1$ とは一致しないが、$R_1R_2^{(-)}R_3$ (積の順序を逆にしたもの) とは一致する** (sympy: [mis_03_derivations.out.txt](outputs/mis_03_derivations.out.txt) (a)。$R_1R_2^{(-)}R_3 - (11.7) = 0$ は成り立ち、$R_3R_2^{(-)}R_1 - (11.7) = 0$ は成り立たない)。成分で見ると、$R_3R_2^{(-)}R_1$ の (1,1) 成分は $c_1c_3 - s_1s_2s_3$ で、(11.7) の $c_1c_3 + s_1s_2s_3$ と符号が違う。

**DIR = 0 の行列 $M_0 = R_3R_2R_1$ との関係**。真の逆は $M_0^\top = R_1^\top R_2^\top R_3^\top = R_1(-\chi_1)\,R_2^{(-)}\,R_3(-\chi_3)$ ($R_2^\top = R_2^{(-)}$、$R_1^\top$ は $s_1$ の符号を替えたもの)。(11.7) は、これの $R_1$ と $R_3$ の符号だけ ($s_1 \to -s_1$、$s_3 \to -s_3$) を元に戻した $R_1(\chi_1)\,R_2^{(-)}\,R_3(\chi_3)$ に等しい。$D = \mathrm{diag}(1, -1, -1)$ ($y$ と $z$ を反転) で挟むと $D\,M_0^\top D = R_1(\chi_1)R_2^{(-)}R_3(\chi_3)$ ($D$ による共役は、$z$ を含む回転 ($R_1$) と $y$ を含む回転 ($R_3$) の角の符号を反転し、$R_2$ は変えない): K64 の DIR ≠ 0 は、**「$y$ と $z$ を反転した座標で書いた逆変換」**。平行移動 $(+\mathrm{DX}, -\mathrm{DY}, -\mathrm{DZ})$ も $y$・$z$ の成分の符号が反転していて、整合する。

**どちらが意図か**。(i) ヘルプの平行移動 $(\mathrm{DX}, -\mathrm{DY}, 0)$ も $y$ 成分が反転しているので、ヘルプも「$y$ を反転した座標」の記述である。(ii) DIR ≠ 0 は「逆向きに使う」ための変換 (ヘルプ: "when you use a line with this element in the reverse direction") なので、回転の順序が DIR = 0 の逆 ($R_1$ が最後) になるのが自然で、ヘルプの $R_3R_2^{(-)}R_1$ (DIR = 0 と同じ順序) は**ヘルプの記述の誤記の可能性が高い**。(iii) ただし、gs で DIR ≠ 0 の COORD を試して確かめてはいない (**未確認**)。ヘルプの行列で作った写像は追跡と $x$ で $5 \times 10^{-5}$、$p_x$ で $6 \times 10^{-4}$ 違い、K64 の行列で作った写像は $10^{-18}$ で一致する ([mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) (a)(b2))。付録 B に書く。

### 4.2 逆変換ではないこと

同じパラメータの COORD (DIR = 0) の次に COORD (DIR = 1) を置いても元に戻らない (入力と出力の差の最大 $2.1 \times 10^{-2}$、同じ条件で真の逆 $R^\top$ では $9 \times 10^{-18}$。(e))。ヘルプが "NOT the inverse" と注意するとおり。真の逆変換は、$\vec r' = R^\top\vec r + \vec d$、$\vec P' = R^\top\vec P$ のあと (11.4)。

COORD が要素の対 (COORD、要素、COORD) の形で使われるときは、2 つ目のパラメータを自分で計算する (ヘルプの注意)。SOL は、この計算を自動でする。

## 5. DX・DY・ROTATE (要素の追跡と CALC)

### 5.1 平行移動と $s$ 軸まわりの回転

平行移動 $x \to x - d_x$、$y \to y - d_y$ は、$p_x, p_y, z, \delta$ を変えない点変換で、正準である (母関数 $F_2 = (x - d_x)p_x' + (y - d_y)p_y' + \ldots$)。$s$ 軸まわりの回転 $R(\theta)$ は、$(x, y)$ と $(p_x, p_y)$ を同じ $R(\theta)$ で回す。面は動かないので面への戻しは要らない (回転軸が面に垂直だから、$z'$ 成分は変わらない)。4×4 の行列 $\mathrm{diag}(R(\theta), R(\theta))$ は、$x$–$p_x$、$y$–$p_y$ の対をそのまま保つ。$|M^\top SM - S| = 2 \times 10^{-17}$ ((b))。

### 5.2 順序と符号

QUAD に DX、DY、ROTATE を付けた追跡を、「$(x - \mathrm{DX}, y - \mathrm{DY})$ を $R(\sigma\theta)$ で回し、誤差なしの QUAD を通し、$R(-\sigma\theta)$ で戻し、DX・DY を足す」と比べた ([mis_02_element_misalign.out.txt](outputs/mis_02_element_misalign.out.txt) (a)、$\theta = 0.1$、DX = 1 mm、DY = 2 mm):

| 手順 | 差 |
|---|---|
| $\sigma = +1$ (DX・DY を引いてから回す) | $2 \times 10^{-19}$ |
| $\sigma = -1$ | $6 \times 10^{-4}$ |
| 先に回してから DX・DY を引く | $1 \times 10^{-4}$ |

したがって入口は

$$
\begin{pmatrix} x\\ y\end{pmatrix}' = R(\theta)\begin{pmatrix} x - \mathrm{DX}\\ y - \mathrm{DY}\end{pmatrix},\quad
\begin{pmatrix} p_x\\ p_y\end{pmatrix}' = R(\theta)\begin{pmatrix} p_x\\ p_y\end{pmatrix},\qquad
R(\theta) = \begin{pmatrix} \cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}
\tag{11.8}
$$

で、出口は (11.8) の逆。ヘルプの「After displacing the magnet by DX and DY, rotate the magnet around the local s-axis by $-$(ROTATE)」は、磁石を $-\theta$ 回すこと、すなわち粒子の座標を $+\theta$ 回すこと。(11.8) と一致する。MULT のずれの 3×3 のヘルプの式 ([08_mult.md](08_mult.md) §5) の $R_3(\chi_3)$ も同じ向き ($\chi_3 = \theta$)。

### 5.3 CALC の `tchge`

CALC (`tchge.f`、入口 `ent = .true.`) は同じ操作を、閉軌道 `cod` と行列 `trans` に施す。

1. 角度 `th`: 直線の要素 (`phig = 0`) では `th = theta + dtheta` (ROTATE と DROTATE の和)、偏向磁石 (`phig ≠ 0`) では `th = theta` (DROTATE は §6 の `tbrote` で別に扱う)。
2. `th ≠ 0` のとき: `cod(1) = cos·x - sin·y`、`cod(3) = sin·x + cos·y`、`cod(2)`・`cod(4)` も同じ回転、`trans1` は $\mathrm{diag}(R, R)$ (4×4 の部分だけ。`trans1(1,1) = cost`、`trans1(1,3) = -sint`、`trans1(3,1) = sint`、`trans1(3,3) = cost` と、運動量の対)。
3. 変位の引き算: `dx = cost*dx0 - sint*dy`、`cod(3) -= sint*dx0 + cost*dy`。これは (11.8) $= R(\theta)(\vec x - \vec d)$ の分解で、$\vec x \to R\vec x$ の後に $R\vec d$ を引く ($R\vec d = (c\,d_x - s\,d_y,\ s\,d_x + c\,d_y)$) 形。$x$ 成分は `tbrote` の中の `cod(1) = cod(1) - dx1` で、$y$ 成分は上で引く。
4. 出口 (`ent = .false.`): 逆の順序 (`cod(3) += ...` のあとで $R^\top$)。

**K64 の既知の欠陥** ([docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md)): B09 (`tchge` の入口で `tbrote` などが `trans` を上書きし、一周行列に入口の回転が入らない)、B114 (`tbrote` の出口の `dz1 = dz*c0 - dx*s0` の符号が入口と逆で、出口のドリフトが逆向きになる)。この章の式は、回転 + 面への戻しの写像としての定義で、これらの欠陥を含む K64 の行列そのものではない。

**運動量の扱い**: CALC の `cod(2)`・`cod(4)` は正準運動量 ($p_0$ で規格化)。回転は $p_x, p_y$ の両方に同じ $R$ を掛けるだけで、$p = 1+\delta$ も $z$ も変えない。追跡 (`tsolrot` の $B_z = 0$ の枝、`chi3 /= 0`) も同じ式を、$q = p_x/p$ に対して書く ($R$ は線形なので $q$ にも同じ形)。

**ソレノイド中の場合**: 磁場 $B_z$ のあるとき、正準運動量 $p_x = \pi_x + a_x$ ($\pi_x$ は力学的運動量、$a_x = -B_z y/2$、規格化、ゲージ $\vec A = (-B_zy, B_zx, 0)/2$) は、原点をずらすとゲージが変わるので変わる。力学的運動量 $\pi$ は物理量で不変なので、原点を $(d_x, d_y)$ ずらして $y' = y - d_y$ にすると、$p_x' = \pi_x + a_x(y') = p_x + a_x(y') - a_x(y) = p_x + B_zd_y/2$ になる。`tsolrot` の `fx = bzs*dy*.5d0`、`px = px + fx/(1+g)` (追跡では $q$ に対して $1/(1+\delta)$ が付く) がこれである。$B_z = 0$ では $f_x = f_y = 0$ で、上の単純な式に戻る。この項は解析だけで、数値では確認していない ([09_solenoid.md](09_solenoid.md) §8)。

## 6. 偏向磁石のロール (DROTATE) と CHI2

偏向磁石の DROTATE と CHI2 (傾き) は、入口と出口の 2 つの面の間で、磁石を「弦の中点を支点にして」回すことに当たる。`tbrote` (CALC)・`tbrot` (追跡) は、§3 と同じ構成である。

面 (入口の面の座標系) で、支点は弦の中点にあり、面の中心の位置の逆向きに $(-d_l s_0,\ 0,\ d_lc_0)$ (支点の位置、$d_l = \mathtt{alg}\,\sin\phi_g/\phi_g$ は弦の半分の長さ、$s_0 = \sin\phi_g$、$c_0 = \cos\phi_g$、$\phi_g = \mathtt{phig}$ は半分の偏向角)。弦の方向は $\vec t = (-s_0, 0, c_0)$、弦に垂直な水平方向は $\vec n = (c_0, 0, s_0)$。入口の回転は

$$
R_{\rm ent} = R_{\vec t}(+\Theta)\,R_{\vec n}(-\chi_2),\qquad R_{\rm exit} = R_{\rm ent}^\top
\tag{11.9}
$$

($R_{\vec a}(\alpha)$ は軸 $\vec a$ のまわりに角 $\alpha$ だけ回す回転。$\Theta$ = DROTATE、$\chi_2$ = CHI2。) 粒子の位置は面の上の $(x, y, 0)$。支点まわりに回して $\vec r' = \vec r_{\rm piv} + R(\vec r - \vec r_{\rm piv})$、$\vec P' = R\vec P$ とし、(11.4) と同じ面への戻しをする。Fortran は、これを $(R - 1)$ を直接作って (`dr11 = dct*c0**2 + dc2*s0**2 + c0*s0*s2*st` など。小さな角で桁落ちしない) 書く:

$$
x_f = x + \mathtt{dr}_{11}\,(x + d_l s_0) + R_{12}\,y + R_{13}\,(-d_lc_0),\quad
z_f = R_{31}(x + d_ls_0) + R_{32}\,y + \mathtt{dr}_{33}\,(-d_lc_0)
\tag{11.10}
$$

(`xf = dr11*xi + r12*yi + dxz + x`、`zf = r31*xi + r32*yi + dzz`、$x_i = x + \mathtt{xl0}$、`xl0` $= d_ls_0$、$z_i = -d_lc_0$、`dxz = r13*zi`、`dzz = dr33*zi`。)

**(11.9) の代数的な確認**。軸 $\vec t = (-s_0, 0, c_0)$ のまわりの角 $\Theta$ の回転と、軸 $\vec n = (c_0, 0, s_0)$ のまわりの角 $-\chi_2$ の回転は、ロドリゲスの式 $R_{\vec a}(\alpha) = 1 + \sin\alpha\,[\vec a]_\times + (1 - \cos\alpha)\,[\vec a]_\times^2$ ($[\vec a]_\times$ は外積の行列) で書ける。この 2 つの積 $R_{\vec t}(\Theta)R_{\vec n}(-\chi_2)$ を成分ごとに展開すると、Fortran の行列 `rr` (入口) の 9 成分

$$
\begin{aligned}
R_{11} &= 1 + (\cos\Theta - 1)c_0^2 + (\cos\chi_2 - 1)s_0^2 + c_0s_0\sin\chi_2\sin\Theta,\quad
R_{22} = \cos\chi_2\cos\Theta,\\
R_{33} &= 1 + (\cos\chi_2 - 1)c_0^2 + (\cos\Theta - 1)s_0^2 - c_0s_0\sin\chi_2\sin\Theta,\\
R_{12} &= s_0\sin\chi_2 - c_0\cos\chi_2\sin\Theta,\quad R_{21} = -\cos\Theta\,s_0\sin\chi_2 + c_0\sin\Theta,\ \ \text{ほか}
\end{aligned}
$$

($R_{22} = 1 + \mathtt{dr}_{22} = 1 + c_2(c_\Theta - 1) + (c_2 - 1) = c_2c_\Theta$)と**すべて一致する** (sympy で差を簡約して 0: [mis_03_derivations.out.txt](outputs/mis_03_derivations.out.txt) (b))。順序を入れ替えた $R_{\vec n}(-\chi_2)R_{\vec t}(\Theta)$ は $2.7 \times 10^{-2}$、$R_{\vec t}(\Theta)R_{\vec n}(+\chi_2)$ は $0.24$、$R_{\vec t}(-\Theta)R_{\vec n}(-\chi_2)$ は $0.39$ ずれる ($\phi = 0.37$、$\Theta = 0.21$、$\chi_2 = -0.13$)。出口の `rr` は入口の転置 (Fortran の `ent = .false.` の代入で $(1,2) \leftrightarrow (2,1)$ などを入れ替えたもの) である。したがって、(11.9) の軸・符号・順序は、数値の総当たりではなく、Fortran の行列からの代数的な恒等式として確定した。$\Theta$ のみのとき $R_{11} = 1 + (\cos\Theta - 1)c_0^2$ は、軸 $\vec t$ の $x$ 成分が $-s_0$ なので $1 - t_x^2 = c_0^2$ から出る ($R_{11} = \cos\Theta + (1 - \cos\Theta)t_x^2 = 1 + (\cos\Theta - 1)(1 - t_x^2)$)。

**(11.10) の位置の式**。$\vec r = (x, y, 0)$、支点 $\vec r_p = d_l\vec t = d_l(-s_0, 0, c_0)$ として、回転した位置は $\vec r' = \vec r_p + R(\vec r - \vec r_p) = \vec r + (R - 1)(\vec r - \vec r_p)$。$\vec r - \vec r_p = (x + d_ls_0,\ y,\ -d_lc_0) = (x_i, y_i, z_i)$ (Fortran の `xi = x + xl0`、`zi = -dl*c0`) なので、$x_f = x + (R_{11} - 1)x_i + R_{12}y_i + R_{13}z_i$ ((11.10) の `xf = dr11*xi + r12*yi + dxz + x`) となる。3 成分とも、$(x_f, y_f, z_f) = \vec r_p + R(\vec r - \vec r_p)$ との差が $4 \times 10^{-17}$ (同 (b)。$z_f$ は元の面 $z = 0$ からの値で、$\vec r' - \vec r$ の $z$ 成分 $+\,\vec r$ の $z$ 成分 0)。向き $\vec P' = R\vec P$ のあと、(11.4) の面への戻しをする。

**確認** ([mis_02_element_misalign.out.txt](outputs/mis_02_element_misalign.out.txt) (c)): `track6._tbrot` (`ent` = 真と偽) を、(11.9) の形で独立に作った写像と比べた。DROTATE のみ、CHI2 のみ、両方の 3 通りで、差は $10^{-17}$ 程度 (入口は DROTATE が $\vec t$ まわり $+\Theta$、CHI2 が $\vec n$ まわり $-\chi_2$ で、$R_{\vec n}$ を先に掛ける。出口は逆)。符号と順序は、上のとおり Fortran の行列 `rr` からも代数的に確かめた。

**`tbgrote` の 6×6 行列** (CALC)。$(R - 1, \vec d, \text{面への戻し})$ の 6×6 の行列で、§3.4 と同じ連鎖律で書かれている (`dpzidpx = -pxi/pzi` など)。Fortran の `trans` の 36 成分のうち、$\partial z/\partial\delta$ の成分 $(5, 6)$ を除く 35 成分を、写像の中心差分 (刻み $10^{-6}$、5 つの乱数点。$x, y \sim 2$ mm、$p_x, p_y \sim 3$ mrad、$z \sim 1$ cm、$\delta \sim 2\times10^{-3}$、$\phi = 0.0654$、$\Theta = 0.02$、$\chi_2 = -0.015$) と比べると、差の最大は $2.7 \times 10^{-12}$ ([mis_03_derivations.out.txt](outputs/mis_03_derivations.out.txt) (c))。**$(5, 6)$ 成分だけは一致しない**: K64 の現在の式 `trans(5,6)=r33*(dpzf*prpzf+dpzi/pzi)-dr(3,3)` は、数値微分と $10^{-5}$ ずれる (値は $-8.5 \times 10^{-5}$ に対し数値微分 $-5.7 \times 10^{-9}$)。同じソースにコメントアウトされた式 `zfpzf*(1 - pr/pzf*r33*dpzidp)` は数値微分と 4 桁で一致する。回転がなくても (`dr = 0`) $dp_{zi}(p_r + 1)/p_{zi}$ が残る。これは K64 の欠陥として既に記録されている (docs/30 の B109。v2 は修正済み、k64 モードは K64 のまま)。pysad の対応関数 (`lattice._tbgrote_mat`) は、v2 で修正した式 (コメントアウトされた式) を持つ。**本章では、解析式 `trans` の各成分を 1 つずつ手で導出するところまでは行っていない** ($\partial z/\partial\delta$ は上のとおり)。35 成分が連鎖律の写像の微分に一致することを数値で確認した。

## 7. $z$ の符号と DZ

ヘルプの MULT の式 (原文で確認) は「$(x, y, \Delta s)_1 = R_3R_2^{(-)}R_1^{(-)}(x - \mathrm{DX}, y - \mathrm{DY}, -\mathrm{DZ})_0$」で、$R_2^{(-)} = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}$、$R_1^{(-)} = \begin{pmatrix}c_1&0&-s_1\\0&1&0\\s_1&0&c_1\end{pmatrix}$ (COORD の $R_2$・$R_1$ と、どちらも $s$ の符号が逆)。

K64 の `tsolrot` (入口) は、09 章 §8.2 で導いたとおり、$R_2'R_1'$ ($R_1' = R_1^{(-)}$ の形、$R_2' = R_2^{(-)}$ の形) を $(x - \mathrm{DX}, y - \mathrm{DY}, s_{0s} - \mathrm{DZ})$ に掛ける。`tmulti` は `-chi1, -chi2` を `tsolrot` に渡す (`tmulti.f:91`)。`tmulti` の `chi1`・`chi2` は `tturn.f:737` で `cmp%value(ky_CHI1_MULT)` などをそのまま受けるので、`tsolrot` に入る角は $-\mathrm{CHI1}$、$-\mathrm{CHI2}$ で、実際に掛かる回転は $R_2'(-\chi_2)R_1'(-\chi_1) = R_2R_1$ (COORD の $R_2$・$R_1$)。つまり**コードの MULT の回転は COORD の DIR = 0 の $R_3R_2R_1$ と同じ符号で、ヘルプの MULT の式 ($R_2^{(-)}$・$R_1^{(-)}$) とは CHI1・CHI2 の符号が逆になっている**。ヘルプの MULT の式の符号の誤記か、$\chi$ の定義の違い (回転の向きの規約) かは、gs で CHI1・CHI2 を指定した MULT を試していないので**未確認**。位置ベクトルの $z$ 成分 ($-\mathrm{DZ}$ に支点までの距離 $s_{0s}$ を足したもの) と、「出口は逆」は、09 章 §8.2 の導出と数値確認 (剛体運動との差 $10^{-17}$) のとおり。

**COORD と MULT の CHI1・CHI2・DZ の対応** は、DIR = 0 の COORD $= R_3R_2R_1(x - \mathrm{DX}, y - \mathrm{DY}, z)$ と MULT $= R_3R_2R_1(x - \mathrm{DX}, y - \mathrm{DY}, -\mathrm{DZ})$ で、行列は同じ (COORD は $z$ に DZ を引かないが、K64 は `zi = pzi*z - dz0` で DZ を引く。§4.1)。MULT は要素の中央を支点にするので、COORD より支点の分の平行移動が加わる (09 章 §8.2 の $\mathrm{ds1}$、$\mathrm{ds2}$ の $s_{0s}$ の項)。

pysad の追跡は直線の MULT の CHI1・CHI2・DZ を使わない ((d): 差 0.0)。したがって、これらの量の数値確認は、K64 の Fortran を numpy に写して剛体運動と比べたもの (09 章 §8.2) まで。

## 8. MARK・APERT・MAP

- **MARK**: 追跡では `case (icMARK)` で何もしない (`tturn.f`)。CALC では行列を作らない (`trans` は恒等)。光学関数の保存や、周期条件の位置の指定の役目だけを持つ (ヘルプ)。確認: 追跡の出力 − 入力は 0.0 ((e))。
- **APERT**: 追跡のみ (ヘルプ)。楕円 $(x - \mathrm{DX})^2/\mathrm{AX}^2 + (y - \mathrm{DY})^2/\mathrm{AY}^2 < 1$ と、箱 $\mathrm{DX1} < x - \mathrm{DX} < \mathrm{DX2}$ (DY も) を通った粒子だけが残る。`tapert` は、通る粒子の座標を変えない。外れた粒子には失う位置 (`kptbl(i,4)`) と周回数 (`kptbl(i,5)`) を書くだけ。確認: 通る粒子の出力 − 入力は 0.0。`PXJ`・`PYJ`・`DPJ` が非零のときは、外れた粒子を楕円の上に乱数で置き直す「ショートカット」の動作になる (読みだけ、未確認)。
- **MAP**: ユーザー定義の写像で、追跡では `ExternalMap("TRACK", l, n, {x, px, y, py, z, dp, ...})` を、CALC では `ExternalMap("OPTICS", l, cod)` を SAD の台本から呼ぶ。追跡へ渡す運動量は `px*(1+pr)` で、正準運動量 (正規化は $p_0$)。物理は台本の定義による。ガイドで扱う物理はない。

## 9. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `mis_01_coord_geometry.py` | (a) ヘルプの行列と `tcoord.f` の $r_{ij}$ (DIR = 0 は一致、DIR ≠ 0 は不一致、`D Mᵀ D` との関係、sympy)、(b) 独立な 3 次元の写像と追跡・CALC の写像の比較、(b2) DIR ≠ 0 のヘルプの式、(c) 6×6 の解析式と数値微分・シンプレクティック性、(d) 小さな角度の 1 次、(e) 逆変換でないこと | [mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) |
| `mis_03_derivations.py` | (a) DIR ≠ 0 の K64 の行列 = $R_1R_2^{(-)}R_3$ (sympy)、(b) `tbrot` の `rr` = $R_{\vec t}(\Theta)R_{\vec n}(-\chi_2)$ (ロドリゲスの式) と支点まわりの位置の式、(c) `tbgrote` の 6×6 行列と数値微分 ((5,6) 成分は K64 の式が合わない: B109)、(d) (11.3)(11.4) の 1 次の係数 | [mis_03_derivations.out.txt](outputs/mis_03_derivations.out.txt) |
| `mis_02_element_misalign.py` | (a) QUAD の DX・DY・ROTATE の順序と符号、(b) $\mathrm{diag}(R, R)$、(c) BEND の `tbrot` (DROTATE・CHI2) と剛体回転の比較、(d) 直線の MULT の CHI1・CHI2・DZ が追跡で無視されること、(e) MARK・APERT | [mis_02_element_misalign.out.txt](outputs/mis_02_element_misalign.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/mis_01_coord_geometry.py` (数秒)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、$\delta = O(10^{-3})$、$z = O(10^{-2})$ m、$x, y = O(10^{-3})$ m、$p_x, p_y = O(10^{-4})$。回転角は $\le 0.5$ rad (COORD)、0.1 rad (ROTATE)、0.02 rad (DROTATE・CHI2)。放射・スピン・ソレノイド磁場・RF は含まない。

**未確認**:
- `tsolrot` (CHI1・CHI2・DZ) は、pysad に実装がなく、K64 も動かしていない。K64 の式を numpy に写した範囲で、剛体運動 + 面への戻しと一致することは確認した (09 章 §8.2)。**ヘルプの MULT の式との CHI1・CHI2 の符号の関係は未確認** (§7)。
- 偏向磁石の (11.9) は、Fortran の行列から代数的に確かめた (§6)。CALC の `tbgrote` の 6×6 行列は、35 成分を数値微分と照合し、$(5, 6)$ 成分は K64 の誤り (B109) と確認した。各成分の手による導出はしていない。
- DIR ≠ 0 の COORD: K64 の行列がヘルプの式と積の順序が逆であること (確定)。どちらが意図かは (ヘルプの誤記の可能性が高いが) gs での確認が未了 (付録 B)。
- ソレノイド中の変位による $p_x$ のゲージのずれ ($f_x$) の数値確認。
- CALC の `tchge` の `trans1`・`trans2` の全成分の導出、`calpol` (スピン) の回転。
- 厳密な「直線で面に戻す」が成り立つのは、要素の端の面の近くに場がないとき。端の場 (フリンジ) が回転した軸に対して持つ構造は、この章では扱わない。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| COORD の DIR = 0 の行列 | **原文で確認** (TeX 版・平文版とも): $R_3R_2R_1$、(11.6)。K64 と一致。 |
| COORD の DIR ≠ 0 | **原文で確認**: $R_3 R_2^{(-)} R_1\,(x, y, z)_0 + (\mathrm{DX}, -\mathrm{DY}, 0)$。K64 の行列は $R_1R_2^{(-)}R_3$ で、積の順序が逆 (付録 B)。「these transformations are NOT the inverse to each other」は正しい (§4.2)。 |
| COORD のキーワード | ヘルプの用法 `COORD name=(DX DY CHI1 CHI2 CHI3 DIR)` とキーワード一覧に DZ が無い。K64 の `tcoord` は DZ を使う。 |
| MULT の misalignments | **原文で確認**: $(x, y, \Delta s)_1 = R_3R_2R_1^{(-)}(x - \mathrm{DX}, y - \mathrm{DY}, -\mathrm{DZ})_0$。「The inverse is applied at the exit」「The geometry of the design orbit is determined by the saved values of CHI1, CHI2, and DZ, while the current values are used for DX, DY, and ROTATE」。コード (`tsolrot` に渡す $-\chi_1$、$-\chi_2$) では、MULT の回転は COORD の DIR = 0 と同じ符号になり、この式と CHI1・CHI2 の符号が逆 (§7。未確認)。 |
| BEND の ROTATE・DROTATE | **原文で確認**: ROTATE は「rotate the magnet around the local s-axis by $-$(ROTATE)」。DROTATE は「does not affect the geometry of the ring」(§6)。DX は「applied before the rotation by ROTATE」(§5.2 の順序と一致)。 |
| APERT | **原文で確認**: 楕円と箱の条件 (TeX 版・平文版)。平文版の最後の式が `Max[DY2, DY2]` で、TeX 版の $\max(\mathrm{DY1}, \mathrm{DY2})$ の誤植とみられる。 |
| ヘルプに無いもの | (1) 面への戻しと $z$ の更新 (§3)。(2) 6×6 の行列 (§3.4)。(3) 偏向磁石の支点まわりの回転 (§6)。 |

## 付録 B. K64 の実装の癖・疑い (別枠。確定には SAD (gs) での確認が要る)

1. **COORD の DIR ≠ 0 の行列の積の順序がヘルプの式と逆** (§4.1)。K64 の行列は $R_1R_2^{(-)}R_3$ ($R_2^{(-)}$ は $s_2$ の符号を逆にした $R_2$)。ヘルプは $R_3R_2^{(-)}R_1$ と書く (DIR = 0 の $R_3R_2R_1$ の $s_2$ の符号だけを変えた形)。K64 の行列は、DIR = 0 の真の逆 $(R_3R_2R_1)^\top$ の $\chi_1$・$\chi_3$ の符号を反転したもので、$D\,M_0^\top D$ ($D = \mathrm{diag}(1, -1, -1)$) に等しい。ヘルプの平行移動 $(\mathrm{DX}, -\mathrm{DY}, 0)$ も $y$ の符号が反転した形なので、「$y$・$z$ を反転した座標での逆変換」が意図で、**ヘルプの積の順序の誤記の可能性が高い**。gs で DIR ≠ 0 の COORD を試して決めてはいない (未確認)。
2. **ヘルプの MULT の misalignments の式の CHI1・CHI2 の符号が、コードの実際の回転と逆** (§7。未確認): コードの MULT は COORD の DIR = 0 と同じ符号の回転。
3. **`tbgrote` の `trans(5,6)` の式の誤り** (B109、§6): コメントアウトされた式 (写像の微分に一致) の代わりに、現在の式は $O(q^2)$ の誤差を持つ (数値微分との比較で確認。v2 は修正済み)。`tchge` の入口が回転行列を上書きする (B09)、`tbrote` の出口の `dz1` の符号 (B114) も既知。
4. 平文版の APERT の式の `Max[DY2, DY2]` は誤植とみられる。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `dx0`・`dy0`・`dz0` | 入口で引く量 (DIR = 0 なら `dx`・`dy`・`dz`、DIR ≠ 0 なら 0) | `tcoord.f`、`tcoorde.f` | `coord._rot` |
| `dx1`・`dy1`・`dz1` | 出口で足す量 (DIR ≠ 0 のとき `dx`・`-dy`・`-dz`) | 同上 | 同上 |
| `r11`〜`r33` | 回転行列 (11.6)(11.7) | 同上 | 同上 |
| `xi`・`yi`・`zi` | スナップショットの位置 (11.3) を平行移動したもの。`zi = pzi*z - dz0`、`xi = x + pxi*z - dx0` | 同上 | `coord.k_coord` |
| `xf`・`yf`・`zf` | 回転後の位置 $\vec r'$ | 同上 | 同上 |
| `pzi`・`pzf` | $p_z$ (回転前・後)。追跡では $1 + \mathtt{sqrt1}(-a)$、CALC では $p\sqrt{1 - a/p^2}$ | 同上 | 同上 |
| `trans1`・`trans2` | 段 1・段 3 の偏微分 (§3.4) | `tcoorde.f` | `coord.coord_map` |
| `th` | `tchge` の $s$ 軸まわりの角 (直線: `theta + dtheta`、偏向: `theta`) | `tchge.f` | `lattice._tchge_cod` |
| `ent` | 入口 (真) か出口 (偽) か | `tchge`、`tbrote` | `ent` |
| `alg`・`phig` | 弦の半分の弧長・半分の偏向角 | `tbrote`、`tbrot` | `alg`・`phig` |
| `dl`・`xl0` | $d_l = \mathtt{alg}\,s_0/\phi_g$、$\mathtt{xl0} = d_ls_0$ | `tbrot`、`tbrote` | `dl`・`xl0` |
| `dr11`・`dr22`・`dr33` | $R_{ii} - 1$ (桁落ちの回避) | `tbrot`、`tbrote` | 同名 |
| `dc`・`dcchi1` など | $\cos a - 1$ (`xsincos` の 4 番目) | `tfloor.f`、`tsol.f` | `_xsincos` |
| `fx`・`fy` | ソレノイド中の変位によるゲージのずれ $\pm B_z d_{y,x}/2$ | `tsol.f` (`tsolrot`) | 未確認 |
