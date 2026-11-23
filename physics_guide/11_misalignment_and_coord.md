# 11. 取り付け誤差と COORD: 座標の正準変換

**この章で分かること**: 要素の取り付け誤差 (DX・DY・DZ・ROTATE・DROTATE・CHI1・CHI2) と COORD 要素が、力学ではなく**座標の正準変換**であること。(1) 粒子の軌跡を 3 次元の直線とみなし、剛体の運動 (平行移動と回転) を施して、新しい面で測り直す、という 1 つの構成で、COORD の追跡 (`tcoord`)・CALC の写像と 6×6 行列 (`tcoorde`)・偏向磁石のロール (`tbrot`・`tbrote`・`tbgrote`) が出ること。(2) 運動量 $p_x, p_y$・$z$・$\delta$ の扱い。(3) DX・DY・ROTATE の順序と符号。(4) MARK・APERT・MAP は力学の写像でないこと。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) (正準変数、$z$ の符号)、[02_drift.md](02_drift.md) (自由な直線運動)。

**読む順**: §1 (結論) → §3 (剛体運動と面への戻し) → §4 (COORD) → §5 (DX・DY・ROTATE) → §6 (偏向磁石のロール)。

## 1. 結論

- 取り付け誤差は、要素の前後で座標を取り替える操作で、要素の中の力学 (ハミルトニアン) は変えない。要素の入口で「ずれと回転を入れた座標」に移り、要素を通り、出口で元の座標に戻る。
- **一般の剛体運動** ($\vec r \to R(\vec r - \vec d)$、$\vec P \to R\vec P$) は、粒子の軌跡が自由な直線であるところ (要素の端の面) では、正準変換である。位置と運動量を同じ回転で回すのは 3 次元の点変換で、正準だからである。回転すると基準面が傾くので、**直線を新しい面 $z' = 0$ まで延ばして測り直す** (§3)。この「回転 + 面への戻し」の写像は、数値で $|J^\top SJ - S| \lesssim 10^{-12}$ (シンプレクティック)。
- 時間は保たれる: 粒子の位置を「基準粒子が面にいる時刻」のスナップショット $\vec r = (x, y, 0) + z\,\vec u$ ($\vec u = \vec P/p$ は進行方向の単位ベクトル。$z > 0$ は先) で表すと、剛体運動のあと新しい面から測った $z' = z'_{\rm 3D}/u'_z$ が出る (式 (11.4))。
- COORD の DIR = 0 は、ヘルプの行列 $R_3R_2R_1$ と一致する。**DIR ≠ 0 の K64 の行列は、ヘルプの式と一致しない** ($R_3R_2R_1$ に $\pm$ の符号を入れ替えた形。付録 B (疑い))。DIR ≠ 0 は DIR = 0 の逆変換でもない (ヘルプが "NOT the inverse" と書く通り)。真の逆は $R^\top$ と平行移動の戻し。
- **ROTATE**: 入口で $(x, y) \to R(+\theta)(x - \mathrm{DX},\, y - \mathrm{DY})$ ($p_x, p_y$ も同じ $R(\theta)$ で回す)、出口で逆。DX・DY は回転の前に引く (QUAD で差 $10^{-19}$)。
- **MARK・APERT・MAP は座標を変える力学ではない**: MARK は恒等写像 (座標は変わらない)、APERT は通る粒子の座標を変えず、失う粒子の記録だけを作る (追跡のみ)、MAP はユーザー定義の写像 `ExternalMap` を呼ぶだけ (§8)。
- pysad の追跡は、直線の MULT の CHI1・CHI2・DZ を使わない (指定しても結果が変わらない。確認した)。K64 の `tsolrot` にはある。BEND の CHI2・DROTATE (`tbrot`) は実装されていて、本章の構成と $10^{-17}$ で一致する (§6)。

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

(11.3)(11.4) を 1 次に展開する ([mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) (d))。$p_z \simeq p$:

- **CHI1 = $\chi_1$** (x–z 面の回転): $p_x' \simeq p_x + \chi_1 p_z$。設計軌道 ($p_x = 0$) は $p_x' = \chi_1$ (+1.0e-6 rad の入力で +1.0e-6)。粒子は進行方向が $\chi_1$ だけ傾く (設計軌道が要素の軸に対して傾く)。
- **CHI2** も同じで $p_y' \simeq p_y + \chi_2 p_z$ (数値: +1.0e-6)。
- **CHI3 (ROTATE)**: 設計軌道は変わらない。
- **DX**: $x' = x - \mathrm{DX}$、$p_x$ は変わらない。
- **DZ**: 設計軌道の位置は変わらないが $z' = z - \mathrm{DZ}$ (数値: $-1.0 \times 10^{-2}$ m)。$p_x \ne 0$ のとき $x' = x + (p_x/p_z)\,\mathrm{DZ}$ ($p_x = 10^{-3}$、DZ = 1 cm で $1.0 \times 10^{-5}$ m)。
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

DIR ≠ 0 のヘルプの式は、$(x, y, z)_1 = R_3\,R_2^{(-)}R_1\,(x, y, z)_0 + (\mathrm{DX}, -\mathrm{DY}, 0)$ で、$R_2^{(-)}$ は $R_2$ の $s_2$ の符号を逆にしたもの。一方、K64 の `dir = .false.` の行列は

$$
\begin{pmatrix}
c_1c_3 + s_1s_2s_3 & s_1s_2c_3 - c_1s_3 & s_1c_2\\
c_2s_3 & c_2c_3 & -s_2\\
-s_1c_3 + c_1s_2s_3 & c_1s_2c_3 + s_1s_3 & c_1c_2
\end{pmatrix},\qquad \text{出口で } (+\mathrm{DX},\ -\mathrm{DY},\ -\mathrm{DZ})\ \text{を足す}
\tag{11.7}
$$

である。(11.7) は、ヘルプの DIR = 0 の行列 $M_0$ について $D M_0^\top D$ ($D = \mathrm{diag}(1, -1, -1)$) に等しい (sympy で確認)。つまり $M_0^\top$ (真の逆回転) の $(1,2), (1,3), (2,1), (3,1)$ 成分の符号を反転したもの。ヘルプの DIR ≠ 0 の式の行列とは一致しない ([mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) (a)(b2): ヘルプの行列で作った写像は、追跡と $x$ で $5 \times 10^{-5}$、$p_x$ で $6 \times 10^{-4}$ 違う。K64 の行列で作った写像は $10^{-18}$ で一致)。どちらが意図かは、SAD (gs) で DIR ≠ 0 の COORD を試さないと決められない。**疑い**として付録 B に書く。

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

**確認** ([mis_02_element_misalign.out.txt](outputs/mis_02_element_misalign.out.txt) (c)): `track6._tbrot` (`ent` = 真と偽) を、(11.9) の形で独立に作った写像と比べた。DROTATE のみ、CHI2 のみ、両方の 3 通りで、差は $10^{-17}$ 程度 (軸の符号と合成の順序は総当たりで探した。入口は DROTATE が $\vec t$ まわり $+\Theta$、CHI2 が $\vec n$ まわり $-\chi_2$ で、$\vec R_{\vec n}$ を先に掛ける。出口は逆)。ただし、符号と順序は数値で当てたもので、Fortran の行列 `rr` からの代数的な導出ではない。軸が $\vec t$ と $\vec n$ になる理由は、$\Theta$ のみのとき $R_{11} = 1 + (\cos\Theta - 1)c_0^2$ ($\vec t$ の $x$ 成分の 2 乗が $s_0^2$ なので、$1 - a_x^2 = c_0^2$) から読める。

CALC の `tbgrote` は、$(R, \vec d, \text{面への戻し})$ の 6×6 の行列で、§3.4 と同じ連鎖律で書かれている (`dpzidpx = -pxi/pzi` など)。pysad の対応関数 (`lattice._tbgrote_mat`) の 1 つずつの導出と数値微分との照合は、この章では行っていない (**未確認**)。

## 7. $z$ の符号と DZ

ヘルプの MULT の式は「$(x, y, \Delta s)_1 = R_3R_2R_1^{(-)}(x - \mathrm{DX}, y - \mathrm{DY}, -\mathrm{DZ})_0$」で、$R_1^{(-)}$ の $s_1$ の符号が COORD の $R_1$ と逆である (ヘルプの MULT の式の $R_1$ は $(c_1, 0, -s_1; 0,1,0; s_1,0,c_1)$)。`tsolrot` は `-chi1`、`-chi2` を引数に渡して呼ぶ (`tmulti` の `call tsolrot(..., -chi1, -chi2, theta2, ...)`) ので、符号の違いはここで吸収されている、と読める。**COORD と MULT の CHI1・CHI2・DZ の符号と回転の支点の対応は、`tsolrot` を実行していないので未確認**。pysad の追跡は直線の MULT の CHI1・CHI2・DZ を使わない ((d): 差 0.0)。`tsolrot` の入口の枝は、`ds1 = schi1*x + dcchi1*s0s - cchi1*dz`、`ds2 = schi2*y + cchi2*ds1 + dcchi2*s0s`、$s_{0s} = -\mathtt{alg}$ (要素の半分の長さ。支点は要素の中央)、`a = ds2/pz2`、`x = x1 - a*px`、`z = z + a` で、§3 の「回転して直線を面に戻す」の形になっている (読み)。

## 8. MARK・APERT・MAP

- **MARK**: 追跡では `case (icMARK)` で何もしない (`tturn.f`)。CALC では行列を作らない (`trans` は恒等)。光学関数の保存や、周期条件の位置の指定の役目だけを持つ (ヘルプ)。確認: 追跡の出力 − 入力は 0.0 ((e))。
- **APERT**: 追跡のみ (ヘルプ)。楕円 $(x - \mathrm{DX})^2/\mathrm{AX}^2 + (y - \mathrm{DY})^2/\mathrm{AY}^2 < 1$ と、箱 $\mathrm{DX1} < x - \mathrm{DX} < \mathrm{DX2}$ (DY も) を通った粒子だけが残る。`tapert` は、通る粒子の座標を変えない。外れた粒子には失う位置 (`kptbl(i,4)`) と周回数 (`kptbl(i,5)`) を書くだけ。確認: 通る粒子の出力 − 入力は 0.0。`PXJ`・`PYJ`・`DPJ` が非零のときは、外れた粒子を楕円の上に乱数で置き直す「ショートカット」の動作になる (読みだけ、未確認)。
- **MAP**: ユーザー定義の写像で、追跡では `ExternalMap("TRACK", l, n, {x, px, y, py, z, dp, ...})` を、CALC では `ExternalMap("OPTICS", l, cod)` を SAD の台本から呼ぶ。追跡へ渡す運動量は `px*(1+pr)` で、正準運動量 (正規化は $p_0$)。物理は台本の定義による。ガイドで扱う物理はない。

## 9. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `mis_01_coord_geometry.py` | (a) ヘルプの行列と `tcoord.f` の $r_{ij}$ (DIR = 0 は一致、DIR ≠ 0 は不一致、`D Mᵀ D` との関係、sympy)、(b) 独立な 3 次元の写像と追跡・CALC の写像の比較、(b2) DIR ≠ 0 のヘルプの式、(c) 6×6 の解析式と数値微分・シンプレクティック性、(d) 小さな角度の 1 次、(e) 逆変換でないこと | [mis_01_coord_geometry.out.txt](outputs/mis_01_coord_geometry.out.txt) |
| `mis_02_element_misalign.py` | (a) QUAD の DX・DY・ROTATE の順序と符号、(b) $\mathrm{diag}(R, R)$、(c) BEND の `tbrot` (DROTATE・CHI2) と剛体回転の比較、(d) 直線の MULT の CHI1・CHI2・DZ が追跡で無視されること、(e) MARK・APERT | [mis_02_element_misalign.out.txt](outputs/mis_02_element_misalign.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/mis_01_coord_geometry.py` (数秒)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、$\delta = O(10^{-3})$、$z = O(10^{-2})$ m、$x, y = O(10^{-3})$ m、$p_x, p_y = O(10^{-4})$。回転角は $\le 0.5$ rad (COORD)、0.1 rad (ROTATE)、0.02 rad (DROTATE・CHI2)。放射・スピン・ソレノイド磁場・RF は含まない。

**未確認**:
- `tsolrot` (CHI1・CHI2・DZ、直線の MULT・SOL) の符号と支点: pysad に実装がなく、K64 は動かしていない。
- 偏向磁石の (11.9) は、軸と符号を数値で当てた (代数的な導出なし)。CALC の `tbgrote` の行列と pysad の対応関数の照合は未実施。
- DIR ≠ 0 の COORD: K64 の行列とヘルプの式の食い違いの決着 (付録 B)。
- ソレノイド中の変位による $p_x$ のゲージのずれ ($f_x$) の数値確認。
- CALC の `tchge` の `trans1`・`trans2` の全成分の導出、`calpol` (スピン) の回転。
- 厳密な「直線で面に戻す」が成り立つのは、要素の端の面の近くに場がないとき。端の場 (フリンジ) が回転した軸に対して持つ構造は、この章では扱わない。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| COORD の DIR = 0 の行列 | **原文で確認** (TeX 版・平文版とも): $R_3R_2R_1$、(11.6)。K64 と一致。 |
| COORD の DIR ≠ 0 | **原文で確認**: $R_3 R_2^{(-)} R_1\,(x, y, z)_0 + (\mathrm{DX}, -\mathrm{DY}, 0)$。K64 と一致しない (付録 B)。「these transformations are NOT the inverse to each other」は正しい (§4.2)。 |
| COORD のキーワード | ヘルプの用法 `COORD name=(DX DY CHI1 CHI2 CHI3 DIR)` とキーワード一覧に DZ が無い。K64 の `tcoord` は DZ を使う。 |
| MULT の misalignments | **原文で確認**: $(x, y, \Delta s)_1 = R_3R_2R_1^{(-)}(x - \mathrm{DX}, y - \mathrm{DY}, -\mathrm{DZ})_0$。「The inverse is applied at the exit」「The geometry of the design orbit is determined by the saved values of CHI1, CHI2, and DZ, while the current values are used for DX, DY, and ROTATE」。コードの符号との対応は §7 のとおり未確認。 |
| BEND の ROTATE・DROTATE | **原文で確認**: ROTATE は「rotate the magnet around the local s-axis by $-$(ROTATE)」。DROTATE は「does not affect the geometry of the ring」(§6)。DX は「applied before the rotation by ROTATE」(§5.2 の順序と一致)。 |
| APERT | **原文で確認**: 楕円と箱の条件 (TeX 版・平文版)。平文版の最後の式が `Max[DY2, DY2]` で、TeX 版の $\max(\mathrm{DY1}, \mathrm{DY2})$ の誤植とみられる。 |
| ヘルプに無いもの | (1) 面への戻しと $z$ の更新 (§3)。(2) 6×6 の行列 (§3.4)。(3) 偏向磁石の支点まわりの回転 (§6)。 |

## 付録 B. K64 の実装の癖・疑い (別枠。確定には SAD (gs) での確認が要る)

1. **COORD の DIR ≠ 0 の行列がヘルプの式と違う** (§4.1)。K64 の行列は $D M_0^\top D$。ヘルプは $R_3R_2^{(-)}R_1$ と書き、$s_2$ の符号だけを変えた形。どちらも直交行列だが一致せず、また逆変換 ($M_0^\top$) でもない。ヘルプの記述の誤りか、K64 の符号の誤りかは未決。出口側で使う変換として「$y$ と $z$ を反転した座標」で書いた逆変換、と読める可能性はある (確定していない)。
2. 平文版の APERT の式の `Max[DY2, DY2]` は誤植とみられる。

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
