# 09. ソレノイド (SOL) と、ソレノイド + 四極の本体

**この章で分かること**: 一様な縦磁場 $B_z$ (ソレノイド) の中の四極の本体の $H_2$ の**厳密解** ($w_1$、$w_2$、$w_+$)、運動学的な補正 $\Delta H$ の厳密解 (ラーモア回転)、$z$ の更新、K64 の `tsolqu` の分割 (Strang)、領域と座標の切替 (BOUND・GEO・CHI3)、端部の $H$ (DISFRIN)。ヘルプの式と Fortran (`tzparams` の変数) の対応表を付ける。

**前提**: [04_quad.md](04_quad.md) (四極の $H_2$、$\Delta H$、Strang)、[03_symplectic_integrators.md](03_symplectic_integrators.md)、[06_edge_fringe.md](06_edge_fringe.md) §4・§9.3 (端部)。ここでは $B_z = 0$ の四極は 04 章、ソレノイドの端の $H$ の導出は 06 章に任せる。

## 1. 結論

1. ソレノイド + 四極の $H_2$ (本体の 2 次の $H$) は線形で、運動方程式は $x'' = b_p y' - k_p x$、$y'' = -b_p x' + k_p y$ ($b_p = b_z/p$、$k_p = K/p$、$' = d/ds$)。特性方程式 $\lambda^4 + b_p^2\lambda^2 - k_p^2 = 0$ から、**振動数 $w_1$ の振動 (2 つ) と、増減率 $w_2$ の指数関数的な成長・減衰 (2 つ)**:
$$
w_1 = \sqrt{\tfrac12\bigl(b_p^2 + V\bigr)},\quad w_2 = \frac{k_p}{w_1},\quad V = \sqrt{b_p^4 + 4k_p^2} = w_1^2 + w_2^2,\quad w_+ = w_1 + w_2
\tag{9.1}
$$
(ヘルプの記号と同じ。sympy で確認: [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (a))。
2. 解は 2 つのモードの重ね合わせ。**振動モード**は $y/x = \pm i\,w_d$ ($w_d = b_p/w_+$)、**指数モード**は $y/x = \pm1/w_d$。初期条件から決まる係数が、ヘルプの $a, b, c, d$ (K64 のコードの `a`・`bw`・`c`・`dw`)。**ヘルプの式、K64 の式、4×4 の行列指数は $10^{-18}$ で一致する** (同 (b)。3 通りの $(B_z, K, \delta)$)。
3. $z$ の更新は $\Delta z = -\tfrac12\int(q_x^2 + q_y^2)\,ds$ ($q = P/p$、$P = p - a$ は運動学的運動量) で、K64 の `tztafs` はこれの閉じた形。ODE の積分と $10^{-19}$ で一致する (同 (c))。
4. $\Delta H = p - p_z - P^2/(2p)$ の解は**ラーモア回転**: 運動学的運動量 $P$ が角 $w\ell$ ($w = b_z(p - p_z)/(pp_z)$) だけ回り、位置はその間の速度 $P(1/p_z - 1/p)$ で動く。ヘルプの $\Delta x$、$\Delta y$ の式には $1/B_z$ の因子が要る (ヘルプの式のままでは ODE と $1.6 \times 10^{-6}$ ずれ、$1/B_z$ をつけると $4 \times 10^{-17}$。付録 A)。
5. K64 の本体は $\Delta H$ の半ステップ + ($H_2$ の厳密解、$\Delta H$) × $n$ の **Strang 分割 (2 次)**。SOL + QUAD + SOL (3 GeV、$b_z = 1.3$ /m、$K = 0.9$ /m²) の追跡は、独立な基準 (端の厳密な流れ + 全 $H$ の ODE) と $10^{-11}$ で一致 (分割の誤差)。DISFRIN = 1 でフリンジの効果 ($10^{-9}$) が消える ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt))。
6. `tztafs` の係数 (`aw1`、`ca1`、`cr2` など) は、(9.8) の積分を 4 つのモードで閉じた形にしたものと、10 個の単項式のすべてで代数的に一致する (§5。`aw1`・`aw2` は手で導出)。領域と座標の切替は、`trots` (領域の端) と `tsolrot` ($B_z = 0$ の枝) が剛体回転 + 面への戻しと一致することを確認した (§8)。確認できなかったこと: `tsgeo` の幾何の式、放射 (`tsolqur`)・スピンの扱い。§14。
7. **K64 の疑い** (付録 B): `tsolqu` の $\Delta H$ の運動量の回転が、更新後の `pxi` を使うため、行列式が $\cos^2 w\ell$ になりシンプレクティックでない (ソースの読みと、pysad のヤコビアンの数値で確かめた疑い。gs での確認は未了)。

## 2. 記号と単位

| 量 | 数学の記号 | K64 の Fortran | pysad | 注釈 |
|---|---|---|---|---|
| 軸上の縦磁場 ($B\rho$ で規格化) | $b_z = B_z/B\rho$ [1/m] | `bz` (SOL の `BZ` [T] を `tfbzs` で規格化。`bzs`) | `bz` | 荷電 × BZ[T]/(p0[eV]/c) |
| $p$ で割った $b_z$ | $b_p = b_z/p$ | `bzp` | `bzp` | 運動学的運動量 $q = P/p$ の方程式に現れる形 |
| 四極 | $K = \mathrm{K1}/L$ [1/m²] | `akk` $= \mathrm{ak}/\mathrm{al}$ | `akk` | `ak` は積分強度 K1 |
| $p$ で割った $K$ | $k_p = K/p$ | `akkp` | `akkp` | |
| 振動数・増減率 | $w_1$、$w_2$ | `w1`、`w2` | `w1`、`w2` | (9.1)。`w2 = akkp/w1` |
| 和・差・積の補助 | $w_+$、$w_-$、$1/V$ | `ws`、`w12`、`wss` | `ws`、`w12`、`wss` | `wss` $= 1/(w_1^2 + w_2^2)$ = ヘルプの $U$ |
| モードの比 | $w_d = b_p/w_+$ | `wd` | `wd` | $y/x = \pm iw_d$ (振動モード) |
| 位相 | $\phi_1 = w_1\ell$、$\phi_2 = w_2\ell$ | `phi1`、`phi2` | `phi1`、`phi2` | `phi2` は $\pm100$ に制限 (`phimax`) |
| 三角・双曲線関数 | $\sin\phi_1$、$\phi_1 - \sin\phi_1$、$\cos\phi_1$、$\cos\phi_1 - 1$ | `s1`、`xs1`、`c1`、`dc1` | 同名 | `xsincos` (`tfloor.f:260`) の戻り値 ($\mathtt{xs} = x - \sin x$、$\mathtt{dc} = \cos x - 1$)。小さい角で桁落ちしない形 |
| | $\sinh\phi_2$、$\phi_2 - \sinh\phi_2$、$\cosh\phi_2$、$\cosh\phi_2 - 1$ | `sh2`、`xsh2`、`ch2`、`dch2` | 同名 | `xsincosh` (`tfloor.f:330`) の戻り値 ($\mathtt{xsh} = x - \sinh x$ は負) |
| スライスの長さ | $\ell = L/n$ | `aln0`、`aln` | `aln` | $n = $ `ndiv` |
| 運動学的運動量 / $p$ | $q = P/p$ | `px`、`py` (追跡の内部変数) | `c[1]`、`c[3]` | ラーモア系では $q_x + b_py/2$ |
| ラーモア回転角 | $w\ell = b_z(p - p_z)\ell/(pp_z)$ | `r*bzp` (`r` $= -\mathtt{dpz}/(1 + \mathtt{dpz})\,\ell$) | | (9.9) |
| $p_z/p - 1$ | $\mathrm{dpz} = \sqrt{1 - q_x^2 - q_y^2} - 1$ | `dpz` (`sqrt1(-ap)`) | | `ap` $= q_x^2 + q_y^2$ |

座標の約束: 運動量は $p_0$ で規格化、$q = p_x/p$ で追跡する (04 章)。ベクトルポテンシャルは $\vec a = (-b_zy/2,\,b_zx/2)$ (一様な $B_z$)。運動学的運動量は $P_x = p_x + b_zy/2$、$P_y = p_y - b_zx/2$ (K64: `pxi = px + bzp*y/2`)。

## 3. ハミルトニアン

$s$ を独立変数とする $H$ (01a 章、曲率 $1/\rho = 0$。$E/(p_0v_0)$ の定数項を除く):

$$
H = -\sqrt{p^2 - P_x^2 - P_y^2} + \frac{K}{2}(x^2 - y^2),\qquad P_x = p_x + \frac{b_z}{2}y,\ P_y = p_y - \frac{b_z}{2}x
\tag{9.2}
$$

($K > 0$ は水平収束。垂直収束の四極は回転で $K > 0$ にそろえる。) 04 章と同じく $H = H_2 + \Delta H$ に分ける:

$$
H_2 = -p + \frac{P_x^2 + P_y^2}{2p} + \frac K2(x^2 - y^2),\qquad
\Delta H = p - p_z - \frac{P_x^2 + P_y^2}{2p}
\tag{9.3}
$$

$H_2$ は線形で厳密に解ける。$\Delta H$ は (§6) ラーモア回転になり厳密に解ける。

**$H_2$ の展開 (回転項と収束項)**。$P_x = p_x + \tfrac{b_z}2y$、$P_y = p_y - \tfrac{b_z}2x$ を 2 乗して足す:

$$
\begin{aligned}
P_x^2 &= p_x^2 + b_z\,y\,p_x + \tfrac{b_z^2}4y^2,\qquad
P_y^2 = p_y^2 - b_z\,x\,p_y + \tfrac{b_z^2}4x^2,\
P_x^2 + P_y^2 &= p_x^2 + p_y^2 - b_z\,(x\,p_y - y\,p_x) + \tfrac{b_z^2}4(x^2 + y^2)
\end{aligned}
$$

$p_\phi = xp_y - yp_x$ (軸まわりの正準角運動量) を使うと

$$
H_2 = -p + \frac{p_x^2 + p_y^2}{2p} - \frac{b_z}{2p}\,p_\phi + \frac{b_z^2}{8p}(x^2 + y^2) + \frac K2(x^2 - y^2)
\tag{9.3a}
$$

(sympy で恒等式を確認: [sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (a)。) 各項の意味: $p_\phi$ の項は、角速度 $-b_z/(2p)$ で位相空間を回す項 (ラーモア回転。$H = \omega\,p_\phi$ は座標を角度 $\omega s$ だけ回す)。$b_z^2$ の項は等方な収束力 (強さ $(b_z/2p)^2$)。

**$K = 0$ の場合の検算**。回転座標 (角速度 $b_p/2$ で回る系) では $p_\phi$ の項が消え、残りは等方な調和振動子 $H = \tfrac{p_x^2 + p_y^2}{2p} + \tfrac{p}{2}\,(b_p/2)^2(x^2 + y^2)$ (振動数 $b_p/2$)。したがって、固定系では (回転と振動の合成で) 振動数 $b_p/2 \pm b_p/2$、すなわち $b_p$ と $0$ の運動になる。これは (9.1) で $K = 0$ とした $w_1 = \sqrt{\tfrac12(b_p^2 + b_p^2)} = b_p$、$w_2 = 0$ と合う (固定系では $w_1 = b_p$。$b_p/2$ はラーモア系の振動数で、$w_1$ ではない)。$K \ne 0$ では四極が回転系で時間に依存するので、この単純な回転では解けず、次節の正規モードを使う。

## 4. $H_2$ の厳密解

### 4.1 運動方程式

正準方程式 $x' = \partial H_2/\partial p_x$、$p_x' = -\partial H_2/\partial x$ を作る。$P_x$、$P_y$ の微分は $\partial P_x/\partial p_x = 1$、$\partial P_y/\partial p_y = 1$、$\partial P_x/\partial y = b_z/2$、$\partial P_y/\partial x = -b_z/2$。$H_2 = -p + (P_x^2 + P_y^2)/(2p) + K(x^2 - y^2)/2$ から

$$
\begin{aligned}
x' &= \frac{P_x}p,\qquad y' = \frac{P_y}p,\
p_x' &= -\frac{\partial H_2}{\partial x} = -\frac{P_y}{p}\Bigl(-\frac{b_z}2\Bigr) - Kx = \frac{b_z}2\frac{P_y}p - Kx,\
p_y' &= -\frac{\partial H_2}{\partial y} = -\frac{P_x}p\cdot\frac{b_z}2 + Ky
\end{aligned}
$$

運動学的運動量 $P_x = p_x + b_zy/2$ の $s$ 微分は $P_x' = p_x' + \tfrac{b_z}2y'$。$y' = P_y/p$ を使うと

$$
P_x' = \frac{b_z}2\frac{P_y}p - Kx + \frac{b_z}2\frac{P_y}p = b_z\frac{P_y}{p} - Kx,\qquad
P_y' = p_y' - \frac{b_z}2x' = -\frac{b_z}2\frac{P_x}p + Ky - \frac{b_z}2\frac{P_x}p = -b_z\frac{P_x}p + Ky
$$

$p$ は $H_2$ の中で定数 ($H_2$ は $s$ によらず、$p$ は共役な $z$ を介してのみ動く) なので、$x'' = P_x'/p$、$y'' = P_y'/p$ から

$$
x'' = b_py' - k_px,\qquad y'' = -b_px' + k_py
\tag{9.4}
$$

($b_p = b_z/p$、$k_p = K/p$、$b_zP_y/p^2 = b_py'$。)

### 4.2 特性方程式と正規モード

$x = Xe^{\lambda s}$、$y = Ye^{\lambda s}$ を (9.4) に代入する。$x'' = \lambda^2X$ などから

$$
\lambda^2X = b_p\lambda Y - k_pX\ \Longrightarrow\ (\lambda^2 + k_p)X - b_p\lambda Y = 0,\qquad
\lambda^2Y = -b_p\lambda X + k_pY\ \Longrightarrow\ b_p\lambda X + (\lambda^2 - k_p)Y = 0
$$

$(X, Y) \ne 0$ の解があるための条件は係数行列の行列式 $= 0$:

$$
(\lambda^2 + k_p)(\lambda^2 - k_p) + b_p^2\lambda^2 = \lambda^4 + b_p^2\lambda^2 - k_p^2 = 0,\qquad
\lambda^2 = \frac{-b_p^2 \pm V}{2},\quad V = \sqrt{b_p^4 + 4k_p^2}
\tag{9.5}
$$

($\lambda^2$ の 2 次方程式の解の公式。) $V > b_p^2$ なので、$+$ の根は正、$-$ の根は負。負の根 $\lambda^2 = -w_1^2$ ($w_1^2 = (b_p^2 + V)/2$) は振動、正の根 $\lambda^2 = +w_2^2$ ($w_2^2 = (V - b_p^2)/2$) は指数関数。足し引きすると $w_1^2 + w_2^2 = V$、$w_1^2 - w_2^2 = b_p^2$、積 $w_1^2w_2^2 = (V^2 - b_p^4)/4 = k_p^2$ から $w_1w_2 = k_p$ (sympy で確認、同 (a))。$b_p^2 = w_1^2 - w_2^2 = w_+w_-$ ($w_\pm = w_1 \pm w_2$) なので

$$
w_d \equiv \frac{b_p}{w_+},\qquad w_d^2 = \frac{w_+w_-}{w_+^2} = \frac{w_-}{w_+}
$$

**固有ベクトル**: 第 1 式から $Y/X = (\lambda^2 + k_p)/(b_p\lambda)$。$k_p = w_1w_2$ を使う。
- $\lambda = \pm iw_1$: $\lambda^2 = -w_1^2$ なので $Y/X = (w_1w_2 - w_1^2)/(\pm ib_pw_1) = -w_-/(\pm ib_p) = \pm i\,w_-/b_p = \pm i\,w_-/(w_dw_+) = \pm i\,w_d$ ($w_-/w_+ = w_d^2$ を使った)。
- $\lambda = \pm w_2$: $Y/X = (w_2^2 + w_1w_2)/(\pm b_pw_2) = \pm w_+/b_p = \pm1/w_d$。

(sympy: 同 (e)。) 実の解にする。$e^{i\phi}$ の $(X, Y) = (1, iw_d)$ の実部・虚部をとると、$x = \cos\phi_1$ に $y = \mathrm{Re}(iw_d\,e^{i\phi_1}) = -w_d\sin\phi_1$、$x = \sin\phi_1$ に $y = \mathrm{Im}(iw_d\,e^{i\phi_1}) = w_d\cos\phi_1$。指数のモードは $e^{\pm w_2s}$ の $(X, Y) = (1, \pm1/w_d)$ なので、$x = \cosh\phi_2$ に $y = \sinh\phi_2/w_d$、$x = \sinh\phi_2$ に $y = \cosh\phi_2/w_d$ ($\cosh = (e^+ + e^-)/2$ に $Y/X$ の符号が $e^-$ で逆になることを使う)。したがって

$$
\begin{aligned}
x &= \bigl[A\cos\phi_1 + B\sin\phi_1\bigr] + \bigl[C\cosh\phi_2 + D\sinh\phi_2\bigr],\
y &= w_d\bigl[-A\sin\phi_1 + B\cos\phi_1\bigr] + \frac1{w_d}\bigl[C\sinh\phi_2 + D\cosh\phi_2\bigr]
\end{aligned}
\tag{9.6}
$$

($\phi_{1,2} = w_{1,2}s$。)

**係数**: 初期条件は、$s = 0$ の位置 $x_0$、$y_0$ と運動学的な傾き $x'(0) = P_x/p \equiv q_{xm}$、$y'(0) = P_y/p \equiv q_{ym}$ ($q_{xm} = q_{x0} + b_py_0/2$、$q_{ym} = q_{y0} - b_px_0/2$。$q_{x0} = p_{x0}/p$ は正準運動量の規格化値。K64 のコードでは `pxi`・`pyi` がこの $q_{xm}$、$q_{ym}$)。(9.6) とその微分を $s = 0$ で評価すると

$$
x_0 = A + C,\quad y_0 = w_dB + \frac D{w_d},\quad q_{xm} = w_1B + w_2D,\quad q_{ym} = -w_dw_1A + \frac{w_2}{w_d}C
$$

この 4 元連立 1 次方程式を $(A, B, C, D)$ について解く。第 1 式から $C = x_0 - A$ を第 4 式に入れ、$q_{ym} = -w_dw_1A + (w_2/w_d)(x_0 - A)$。$A$ の係数は $-(w_dw_1 + w_2/w_d) = -(w_d^2w_1 + w_2)/w_d$。$w_d^2w_1 + w_2 = [w_-w_1 + w_+w_2]/w_+ = (w_1^2 + w_2^2)/w_+ = V/w_+$ なので ($w_-w_1 + w_+w_2 = w_1^2 - w_1w_2 + w_1w_2 + w_2^2$)

$$
A = \frac{(w_2/w_d)\,x_0 - q_{ym}}{V/(w_dw_+)} = \frac{w_2w_+x_0 - w_dw_+q_{ym}}{V} = \frac{w_2w_+x_0 - b_pq_{ym}}{V}
$$

($w_dw_+ = b_p$。) 同様に ($B$ と $D$ は第 2、第 3 式の連立で、$C$ と $A$ は上の関係) 解くと、K64 の `a`、`bw`、`c`、`dw` と

$$
A = a,\quad B = bw,\quad C = b_p\,c,\quad D = b_p\,dw
\tag{9.7}
$$

の関係にある。

$$
a = \frac{w_2w_+x_0 - b_pq_{ym}}V,\quad
bw = \frac{w_+q_{xm} - b_pw_2y_0}V,\quad
c = \frac{w_1w_dx_0 + q_{ym}}V,\quad
dw = \frac{-w_dq_{xm} + w_1y_0}V
$$

(sympy で 4 元連立方程式を解いた結果と上の 4 式の差が恒等的に 0: [sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (b)。$w_d^2 = w_-/w_+$ を使って簡約した。4×4 を数値で解いた $(A, B, C, D)$ と K64 の式は $4 \times 10^{-19}$ で一致: [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (e)。)

### 4.3 ヘルプの式との対応

ヘルプの solution-H2 の記号と K64 のコードの対応 (変数は §2 の表):

| ヘルプ | K64 | 関係 |
|---|---|---|
| $U$ | `wss` | $1/V$ |
| $a$、$c$ | `a`、`c` | 同じ |
| $b$、$d$ | `bw`、`dw` | $b = w_1\cdot bw$、$d = w_2\cdot dw$ |
| $u_1, u_2$ | `u1w`、`u2w` | $u_{1,2} = w_1\,\mathtt{u1w,u2w}$ (`u1w` $= a\,\mathtt{dc1} + bw\,\mathtt{s1}$、`u2w` $= -a\,\mathtt{s1} + bw\,\mathtt{dc1}$) |
| $v_1, v_2$ | `v1w`、`v2w` | $v_{1,2} = w_2\,\mathtt{v1w,v2w}$ (`v1w` $= c\,\mathtt{dch2} + dw\,\mathtt{sh2}$、`v2w` $= c\,\mathtt{sh2} + dw\,\mathtt{dch2}$) |
| $\Delta x = u_1/w_1 + v_1B_z/(pw_2)$ | `x += u1w + bzp*v1w` | |
| $\Delta y = u_2B_z/(pw_1w_+) + w_+v_2/w_2$ | `y += wd*u2w + ws*v2w` | $w_d = b_p/w_+$ |
| $p_x$ の式 | `px = pxi + u2w*w1 + v2w*w2*bzp` (ラーモア系) | 正準運動量に直すと同じ |
| $p_y$ の式 | `py = pyi - wd*u1w*w1 + ws*v1w*w2` | |

(ヘルプの $p_x = p_{x0} + pu_2 + B_z(v_2 - \Delta y/2)$ は正準運動量 ($p_{x0}$ は正準) で、K64 は運動学的運動量 $q$ のラーモア系で更新して、最後に $-b_py/2$ を引いて正準 (の $q$) に戻す。)

**検証** ([sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (b)): $(B_z, K, \delta, \ell) = (1.3, 0.9, 0.01, 0.4)$、$(2.5, 1.2, -0.02, 0.4)$、$(0.3, 4.0, 0.02, 0.4)$ の 3 通りで、ヘルプの式と K64 の式 (pysad の `_TZ` で作った $w_1$ など + `tsolqu.f:117-142` の更新) は、4×4 の行列指数 $e^{A\ell}$ と $4 \times 10^{-19} \sim 2 \times 10^{-18}$ で一致 ($x_0 = 3$ mm、$p_x = 2$ mrad など)。

## 5. $z$ の更新

### 5.1 $-\tfrac12\int q^2ds$

$z$ は $p$ と共役。$H_2$ の $p$ 微分 (正準変数 $x, p_x$ などを固定。$P$ は $p$ によらない) は $\partial H_2/\partial p = -1 - \dfrac{P_x^2 + P_y^2}{2p^2}$。定数 $-1$ は基準軌道の経路 ($E/(p_0v_0)$ の項と合わせて $-dv$ の項) に入るので、残りが

$$
\Delta z_{H_2} = -\frac12\int_0^\ell\bigl(q_x^2 + q_y^2\bigr)\,ds,\qquad q = P/p
\tag{9.8}
$$

(04 章 §4.3 の導出と同じ形。) $q_x = P_x/p = x'$、$q_y = y'$ なので、(9.6) を $s$ で微分した

$$
x' = w_1(-A\sin\phi_1 + B\cos\phi_1) + w_2(C\sinh\phi_2 + D\cosh\phi_2),\quad
y' = w_dw_1(-A\cos\phi_1 - B\sin\phi_1) + \frac{w_2}{w_d}(C\cosh\phi_2 + D\sinh\phi_2)
$$

を (9.8) に入れる。$x'^2 + y'^2$ は、振動モードだけの項、指数モードだけの項、両者の交差項に分かれる:

$$
\begin{aligned}
x'^2 + y'^2 ={}& w_1^2\bigl[(-A\sin\phi_1 + B\cos\phi_1)^2 + w_d^2(A\cos\phi_1 + B\sin\phi_1)^2\bigr]\
&+ w_2^2\bigl[(C\sinh\phi_2 + D\cosh\phi_2)^2 + w_d^{-2}(C\cosh\phi_2 + D\sinh\phi_2)^2\bigr]\
&+ 2w_1w_2\bigl[(-A\sin\phi_1 + B\cos\phi_1)(C\sinh\phi_2 + D\cosh\phi_2) - (A\cos\phi_1 + B\sin\phi_1)(C\cosh\phi_2 + D\sinh\phi_2)\bigr]
\end{aligned}
$$

(交差項では $y'$ の $w_d\cdot w_d^{-1} = 1$ が相殺する。)

### 5.2 積分の表

$0 \to \ell$ の積分に出る 10 種類は、すべて初等関数。$s_1 = \sin\phi_1$、$c_1 = \cos\phi_1$、$\mathrm{sh}_2 = \sinh\phi_2$、$\mathrm{ch}_2 = \cosh\phi_2$、$\phi_{1,2} = w_{1,2}\ell$、$V = w_1^2 + w_2^2$ として (sympy の積分との差 $10^{-16}$: [sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (c3))

| 積分 ($\int_0^\ell\cdots ds$) | 値 | 導出 |
|---|---|---|
| $\cos^2w_1s$、$\sin^2w_1s$ | $\dfrac\ell2 \pm \dfrac{s_1c_1}{2w_1}$ | $\cos^2 = (1 + \cos2w_1s)/2$ |
| $\sin w_1s\cos w_1s$ | $\dfrac{s_1^2}{2w_1}$ | $\tfrac12\sin2w_1s$ の積分 |
| $\cosh^2w_2s$、$\sinh^2w_2s$ | $\pm\dfrac\ell2 + \dfrac{\mathrm{sh}_2\mathrm{ch}_2}{2w_2}$ | $\cosh^2 = (\cosh2w_2s + 1)/2$ |
| $\sinh w_2s\cosh w_2s$ | $\dfrac{\mathrm{sh}_2^2}{2w_2}$ | |
| $\cos w_1s\cosh w_2s$ | $\dfrac{w_1s_1\mathrm{ch}_2 + w_2c_1\mathrm{sh}_2}{V}$ | 部分積分を 2 回して $V$ で割る |
| $\cos w_1s\sinh w_2s$ | $\dfrac{w_1s_1\mathrm{sh}_2 + w_2c_1\mathrm{ch}_2 - w_2}{V}$ | 同 |
| $\sin w_1s\cosh w_2s$ | $\dfrac{-w_1c_1\mathrm{ch}_2 + w_2s_1\mathrm{sh}_2 + w_1}{V}$ | 同 |
| $\sin w_1s\sinh w_2s$ | $\dfrac{-w_1c_1\mathrm{sh}_2 + w_2s_1\mathrm{ch}_2}{V}$ | 同 |

(交差項の積分の分母 $w_1^2 + w_2^2 = V$ が、K64 の `wss` $= 1/V$ の因子として現れる。たとえば $I = \int_0^\ell\cos w_1s\cosh w_2s\,ds$ は、$\cos$ 側を 1 回、$\cosh$ 側を 1 回部分積分すると $I = [\sin w_1s\cosh w_2s/w_1]_0^\ell - (w_2/w_1)\int\sin w_1s\sinh w_2s\,ds$、さらに $\int\sin\sinh = [-\cos w_1s\sinh w_2s/w_1]_0^\ell + (w_2/w_1)I$ となり、代入して $I$ について解くと $(1 + w_2^2/w_1^2)\,I = [\sin\cosh/w_1 + w_2\cos\sinh/w_1^2]$、すなわち $(w_1^2 + w_2^2)\,I = [w_1\sin\phi_1\cosh\phi_2 + w_2\cos\phi_1\sinh\phi_2]$ の形に閉じる。) $\sin2\phi$、$\sinh2\phi$ は倍角で $s_1c_1$、$\mathrm{sh}_2\mathrm{ch}_2$ にした。

(9.8) の結果は $(A, B, C, D)$ の 2 次形式で、10 個の係数を持つ。sympy で積分した閉じた形は、数値積分と $10^{-17}$ で一致する (同 (c)、$w_1 = 0.9$、$w_2 = 0.35$、$\ell = 0.4$、$\Delta z = -1.2285626\times10^{-6}$)。

### 5.3 K64 の変数 $(\mathtt{awu}, \mathtt{dwu}, \mathtt{pxi}, \mathtt{pyi})$ への書き換え

K64 は $(A, B, C, D)$ の代わりに、`awu` $= a\,w_1/w_+$、`dwu` $= dw\cdot w_2$ と、入口のラーモア系の傾き `pxi` $= q_{xm}$、`pyi` $= q_{ym}$ の 4 つで $\Delta z$ を表す。(9.7) と §4.2 の $a, bw, c, dw$ の式から、$(x_0, y_0)$ を $(\mathtt{awu}, \mathtt{dwu})$ に取り替えると (sympy で恒等式を確認: 同 (c2))

$$
A = \frac{w_+}{w_1}\,\mathtt{awu},\qquad
B = \frac{\mathtt{pxi} - b_p\,\mathtt{dwu}}{w_1},\qquad
C = \frac{w_-}{w_2}\,\mathtt{awu} + \frac{w_d}{w_2}\,\mathtt{pyi},\qquad
D = \frac{b_p}{w_2}\,\mathtt{dwu}
\tag{9.8a}
$$

$\Delta z$ は $(\mathtt{awu}, \mathtt{dwu}, \mathtt{pxi}, \mathtt{pyi})$ の 2 次形式で、10 個の単項式を持つ。K64 の `tsolqu.f:137-141` は、これを

$$
\Delta z = \mathtt{bzp}\Bigl(-\frac{\mathtt{awu}\,\mathtt{dwu}\,\mathtt{dxs}^2}{\mathtt{akkp}} + \mathtt{ca1}\,\mathtt{pxi}\,\mathtt{pyi}\,\mathtt{wss}\Bigr) + \mathtt{dz1} + \mathtt{dz2}
$$

と書く。`dz1`、`dz2` は `tztafs` の戻り値で

$$
\mathtt{tztafs} = \mathtt{zwu}\Bigl(\tfrac12\mathtt{zwu}\,f_3 + f_2\,\mathtt{pxt} + \mathtt{bzp}\,f_1\,\mathtt{pyt}\Bigr) - \tfrac12\mathtt{zw1}\,\mathtt{pyt}^2
$$

を、第 1 の呼び出し $(\mathtt{zwu}, \mathtt{pxt}, \mathtt{pyt}, \mathtt{zw1}) = (\mathtt{awu}, \mathtt{pxi}, \mathtt{pyi}, \mathtt{aw1})$ と第 2 の呼び出し $(-\mathtt{dwu}, -\mathtt{pyi}, \mathtt{pxi}, \mathtt{aw2})$ (引数 `wst`・`w12t`・`swss` も入れ替えて符号を変える) で評価したもの。単項式への割り当ては

| 単項式 | K64 の項 |
|---|---|
| $\mathtt{pyi}^2$ | $-\tfrac12\mathtt{aw1}$ (第 1 呼び出し) |
| $\mathtt{pxi}^2$ | $-\tfrac12\mathtt{aw2}$ (第 2 呼び出し) |
| $\mathtt{awu}^2$、$\mathtt{dwu}^2$ | $\tfrac12f_3$ (第 1、第 2 呼び出し) |
| $\mathtt{awu}\,\mathtt{pxi}$、$\mathtt{dwu}\,\mathtt{pyi}$ | $f_2$ (第 1)、$+f_2$ (第 2) |
| $\mathtt{awu}\,\mathtt{pyi}$、$\mathtt{dwu}\,\mathtt{pxi}$ | $\mathtt{bzp}\,f_1$ (第 1)、$-\mathtt{bzp}\,f_1$ (第 2) |
| $\mathtt{awu}\,\mathtt{dwu}$ | $-\mathtt{bzp}\,\mathtt{dxs}^2/\mathtt{akkp}$ |
| $\mathtt{pxi}\,\mathtt{pyi}$ | $\mathtt{bzp}\,\mathtt{ca1}\,\mathtt{wss}$ |

**手で導ける 2 つの係数** ($\mathtt{pxi}^2$、$\mathtt{pyi}^2$)。$\mathtt{awu} = \mathtt{dwu} = \mathtt{pyi} = 0$ とすると (9.8a) から $(A, B, C, D) = (0, \mathtt{pxi}/w_1, 0, 0)$、つまり振動モードの 1 つだけ ($B$ モード)。このとき

$$
\Delta z = -\frac12\,\frac{\mathtt{pxi}^2}{w_1^2}\int_0^\ell\bigl[w_1^2\cos^2\phi_1 + w_d^2w_1^2\sin^2\phi_1\bigr]ds
= -\frac{\mathtt{pxi}^2}2\Bigl[\frac{1 + w_d^2}2\ell + \frac{1 - w_d^2}{2w_1}s_1c_1\Bigr]
$$

($\int\cos^2$、$\int\sin^2$ は §5.2 の表。) $1 + w_d^2 = 2w_1/w_+$、$1 - w_d^2 = 2w_2/w_+$ ($w_d^2 = w_-/w_+$) を入れると

$$
\Delta z\Big|_{\mathtt{pxi}^2} = -\frac{\mathtt{pxi}^2}{2}\cdot\frac{w_1^2\ell + w_2s_1c_1}{w_1w_+} = -\frac12\,\mathtt{aw2}\,\mathtt{pxi}^2,\qquad
\mathtt{aw2} = \ell + \frac{\mathtt{cxs1}\,w_2}{w_1w_+}
$$

と、`aw2` の定義 (§10。`cxs1` $= s_1c_1 - \phi_1$。$\ell + (s_1c_1 - w_1\ell)w_2/(w_1w_+) = (w_1^2\ell + w_2s_1c_1)/(w_1w_+)$ に等しい) と一致する。同様に、$\mathtt{pyi}$ だけの場合は $(A, B, C, D) = (0, 0, w_d\,\mathtt{pyi}/w_2, 0)$ ($C$ モードだけ) で、

$$
\Delta z = -\frac12\,\frac{w_d^2\mathtt{pyi}^2}{w_2^2}\int_0^\ell\bigl[w_2^2\sinh^2\phi_2 + w_2^2w_d^{-2}\cosh^2\phi_2\bigr]ds
= -\frac{\mathtt{pyi}^2}2\Bigl[w_d^2\Bigl(\frac{\mathrm{sh}_2\mathrm{ch}_2}{2w_2} - \frac\ell2\Bigr) + \frac\ell2 + \frac{\mathrm{sh}_2\mathrm{ch}_2}{2w_2}\Bigr]
= -\frac{\mathtt{pyi}^2}2\cdot\frac{w_2^2\ell + w_1\mathrm{sh}_2\mathrm{ch}_2}{w_2w_+} = -\frac12\,\mathtt{aw1}\,\mathtt{pyi}^2
$$

($\mathtt{aw1} = \ell + \mathtt{cxs2}\,w_1/(w_2w_+)$、`cxs2` $= \mathrm{sh}_2\mathrm{ch}_2 - \phi_2$)。これらの対応は sympy でも一致 (同 (c2))。

**残りの係数**。$\mathtt{awu}^2$ ($A$ と $C$ モードの混合)、$\mathtt{awu}\,\mathtt{dwu}$ などは、(9.8a) を §5.1 の 2 次形式に代入し、§5.2 の表の積分 (特に交差項の $1/V$) で整理した係数で、手で展開すると項数が多い。ここでは**代数的な恒等式として sympy で確認した**: $\Delta z$ の厳密な積分と、K64 の `dz` ($=\mathtt{bzp}(\ldots) + \mathtt{dz1} + \mathtt{dz2}$、§10 の定義をそのまま写したもの) の差の分子を、$\sin^2 + \cos^2 = 1$、$\cosh^2 - \sinh^2 = 1$、$w_d^2 = (w_1 - w_2)/(w_1 + w_2)$ で簡約すると、**全体で 0、10 個の単項式ごとでも 0** ([sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (c))。これは、`tztafs` の係数 (`aw1`、`aw2`、`g1`、`g2`、`ca1`、`cr2`、`cr3`、`csw1`、`csws`、`dcw1`、`dcw2`、`dxs`) が、(9.8) の閉じた形の係数を過不足なくまとめたものであることの確認になる (各係数が 10 個の単項式のどの部分に入るかを 1 つずつ手で追ってはいない)。

## 6. $\Delta H$ の解 (ラーモア回転)

$\Delta H = p - p_z - \dfrac{P^2}{2p}$ ($p_z = \sqrt{p^2 - P^2}$、$P^2 = P_x^2 + P_y^2$) は $P$ だけの関数 ($K = 0$)。$\partial p_z/\partial P_x = -P_x/p_z$ なので

$$
\frac{\partial\Delta H}{\partial P_x} = \frac{P_x}{p_z} - \frac{P_x}{p} = P_x\Bigl(\frac1{p_z} - \frac1p\Bigr) \equiv v_x,\qquad v_y\ \text{も同様}
$$

$\Delta H$ は $x, y$ に陽には依存せず、$P$ を通してだけ依存する。正準方程式は $x' = \partial\Delta H/\partial p_x = v_x$ ($\partial P_x/\partial p_x = 1$)、$y' = v_y$、$p_x' = -\partial\Delta H/\partial x = -v_y\,\partial P_y/\partial x = +\tfrac{b_z}2v_y$、$p_y' = -v_x\,\partial P_x/\partial y = -\tfrac{b_z}2v_x$:

$$
x' = v_x,\quad y' = v_y,\quad p_x' = \frac{b_z}{2}v_y,\quad p_y' = -\frac{b_z}{2}v_x
$$

運動学的運動量は $P_x' = p_x' + \tfrac{b_z}{2}y' = \tfrac{b_z}2v_y + \tfrac{b_z}2v_y = b_zv_y$、$P_y' = p_y' - \tfrac{b_z}2x' = -b_zv_x$。$v = P\,(1/p_z - 1/p)$ なので $b_zv = w\,P$ ($w \equiv b_z(1/p_z - 1/p) = b_z(p - p_z)/(pp_z)$)。したがって

$$
P_x' = w\,P_y,\qquad P_y' = -w\,P_x
$$

($|P|$ が一定なので $p_z$ も一定で、$w$ は定数。) これは角速度 $w$ の回転で、解は $P_x(s) = P_x\cos ws + P_y\sin ws$、$P_y(s) = P_y\cos ws - P_x\sin ws$。位置は $x' = v_x = (w/b_z)\,P_x(s)$ を積分する:

$$
\begin{aligned}
\Delta x &= \frac w{b_z}\int_0^\ell(P_x\cos ws + P_y\sin ws)\,ds = \frac1{b_z}\bigl[P_x\sin w\ell - P_y(\cos w\ell - 1)\bigr],\\
\Delta y &= \frac w{b_z}\int_0^\ell(P_y\cos ws - P_x\sin ws)\,ds = \frac1{b_z}\bigl[P_x(\cos w\ell - 1) + P_y\sin w\ell\bigr],\qquad w = \frac{b_z(p - p_z)}{pp_z}
\end{aligned}
\tag{9.9}
$$

($P_{x,y}$ は初期値の運動学的運動量。) 正準運動量は $p_x = p_{x0} + \tfrac{b_z}2\Delta y$、$p_y = p_{y0} - \tfrac{b_z}2\Delta x$ ($p_x' = \tfrac{b_z}2v_y$ を積分)。

**$z$ の更新**: $\partial\Delta H/\partial p = 1 - p/p_z + P^2/(2p^2)$ ($\partial p_z/\partial p = p/p_z$)。$P^2 = p^2 - p_z^2$ を入れると

$$
\Delta z = \Bigl(1 - \frac p{p_z} + \frac{p^2 - p_z^2}{2p^2}\Bigr)\ell = \Bigl(\frac32 - \frac p{p_z} - \frac{p_z^2}{2p^2}\Bigr)\ell
$$

**K64 の対応** (`tsolqu.f:107-123`): `dpz` $= p_z/p - 1$、`ap` $= q_x^2 + q_y^2 = 1 - (p_z/p)^2 = -\mathtt{dpz}\,(2 + \mathtt{dpz})$ ($q = P/p$)。`r` $= -\mathtt{dpz}/(1 + \mathtt{dpz})\,\ell$ とおくと回転角は

$$
r\,b_p = -\frac{\mathtt{dpz}}{1 + \mathtt{dpz}}\cdot\frac{b_z}p\,\ell = \frac{1 - p_z/p}{p_z/p}\cdot\frac{b_z}p\,\ell = \frac{b_z(p - p_z)}{p\,p_z}\,\ell = w\ell
$$

(第 2 の等号で $\mathtt{dpz} = p_z/p - 1$ を入れた。) `call xsincos(r*bzp,...)` が回転 (`a22` $= \cos w\ell$、`a24` $= \sin w\ell$、`a14` $= \cos w\ell - 1$)。位置の更新 `x += (a24*pxi - a14*py)/bzp`、`y += (a14*pxi + a24*py)/bzp` は、`pxi`、`py` が $q = P/p$、`bzp` $= b_z/p$ なので $P/b_z = q/b_p$ となり、(9.9) と同じ ($-\mathtt{a14} = 1 - \cos$)。$z$ の更新 $-(3 + \mathtt{dpz})\,\mathtt{ap}\,r/(2(2 + \mathtt{dpz}))$ に `ap`、`r` を代入すると、$d \equiv \mathtt{dpz}$ として

$$
-\frac{(3 + d)\cdot(-d(2 + d))}{2(2 + d)}\cdot\Bigl(-\frac{d}{1 + d}\ell\Bigr) = -\frac{(3 + d)\,d^2}{2(1 + d)}\,\ell
$$

一方、ヘルプの式で $p_z/p = 1 + d$ とすると $\dfrac32 - \dfrac1{1+d} - \dfrac{(1+d)^2}2 = \dfrac{3(1+d) - 2 - (1+d)^3}{2(1+d)} = -\dfrac{3d^2 + d^3}{2(1+d)}$ で、上と一致する (sympy でも確認: [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (d)、差 0)。数値の流れ ($\Delta H$ の正準方程式、DOP853、$B_z = 2$ /m、$\ell = 0.5$ m、$p = 1.02$) とは、(9.9) で $3.6 \times 10^{-17}$ で一致する。

## 7. 本体の分割 (`tsolqu`)

K64 の `tsolqu` (`tsolqu.f:1-166`) は $B_z \ne 0$ のとき:

1. 分割数 `ndiv = 1 + int(abs(al*hypot(ak, bz)/eps))` (`ak` は積分強度 K1、`al` は全長 $L$。$\mathrm{hypot}(K1, b_z)\,L/\varepsilon$ の整数部 + 1)。`eps` $= 0.2\times\mathrm{EPS}$ (EPS = 0 なら 0.2)。スライスの長さ $\ell = $ `aln0` $= L/n$。**この数の根拠**: Strang 分割の局所誤差は、$H_2$ と $\Delta H$ の交換子による $O(\ell^3)$ で、$H_2$ の固有の角 ($w_1\ell$、$w_2\ell$。$\ell\,\mathrm{hypot}(K, b_z)$ の程度) が 1 より十分小さいことが分割の条件になる、というところまでは読める。$\mathrm{hypot}(K1, b_z)\,L/n < \varepsilon$ の形はそれに合う。係数 0.2 は経験的な値で、ソースにコメントが無く、$\Delta H$ の大きさ (傾き $q^2$) を含む誤差の見積もりからは導けない (未導出)。
2. 入射でラーモア系に: `px += bzp*y/2`、`py -= bzp*x/2` ($q \to$ 運動学的運動量 / $p$)。
3. 各スライス: $\Delta H$ の流れ (回転 + ドリフト、最初と最後は半分、途中は隣り合う 2 つの半分を 1 つにまとめた長さ `aln0`) → $H_2$ の厳密解 (§4。長さ `aln0`) → …、最後に $\Delta H$ の半分。**$\Delta H$ の回転の運動量の更新 (`pyi = -a24*pxi + a22*py(i)`) は、更新後の `pxi` を使っている (付録 B。疑い)**。
4. 出射で元に戻す: `px = … - bzp*y/2`、`py = … + bzp*x/2`。

つまり $\Delta H(\tfrac\ell2)\,H_2(\ell)\,\Delta H(\ell)\,H_2(\ell)\cdots H_2(\ell)\,\Delta H(\tfrac\ell2)$ の **Strang 分割 (2 次)** ([03_symplectic_integrators.md](03_symplectic_integrators.md))。$B_z = 0$ の四極 (04 章) も同じ構造だが、$H_2$ は $\cos$、$\cosh$ の厳密解 (`tzsetparam0`: `w1 = sqrt(akkp)`、`phi1`、`phi2` が同じ $\phi_1$)。

**全体の確認** ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt)): SOL ($B_z = 13.0$ T = $b_z$ 1.3 /m、BOUND = 1) → QUAD ($K = 0.9$ /m²、$L = 0.4$ m、$n = 3$) → SOL ($B_z = 0$) の追跡と、独立な基準 (入口の $\exp(:f:)$ ((6.10) の厳密な流れ) → (9.2) の $H$ の正準方程式の DOP853 → 出口の $\exp(:f:)$) の差は、3 つの入射で $7 \times 10^{-12} \sim 1.3 \times 10^{-11}$ (3 mm・2 mrad、$\delta = +0.01, 0, -0.02$)。**この差は分割 (Strang、$n = 3$) の誤差で、SOL・QUAD・端のフリンジの符号や単位の規約を含めて、追跡が全 $H$ の流れと一致する**。

## 8. 領域と座標の切替 (BOUND・GEO・CHI3)

ヘルプの要約: SOL 要素は、ソレノイドの境界と、場が変わる所に置く。SOL の領域では、座標は**ソレノイドの軸**にとる (設計軌道が曲がっていても)。BOUND = 1 を領域の両端に置く。領域の中に DRIFT、BEND (直線)、QUAD、MULT を置いてよい (縦磁場が重なる)。領域の出口で座標が設計軌道に戻り、CHI3 は、出口の後の最初の GEO = 1 の MARK で 0 になるように決まる。DX・DY・DPX・DPY は GEO = 1 の SOL で設計軌道の軸に対するずれと傾きを与える。BZ だけを変える SOL (BOUND = 0) は座標を変えない。

実装は、役割で 2 つのサブルーチンに分かれる (`tsol.f`)。**領域の端での座標の取り替え**は `trots` (§8.1)、**領域の中の要素 (MULT・QUAD など) の取り付け誤差 (DX・DY・DZ・CHI1・CHI2・CHI3)** は `tsolrot` (§8.2、§8.3)。

### 8.1 領域の端: `trots` (`tsol.f:819`)

`tsol` (`tsol.f:565`、`:764`) は、領域の入口・出口の SOL 要素で `trots` を呼ぶ。引数は、`tpara` が CHI1・CHI2・CHI3 から作る回転行列 $R$ (`R11`〜`R33`。$R = R_3R_2R_1$。11 章 (11.6) と同じ形) と、SOL の DX・DY・DZ。この位置は領域の**外側** ($B_z = 0$) なので、正準運動量 $=$ 運動学的運動量で、11 章 §3 の「剛体回転 + 直線を面へ戻す」がそのまま使える。

**入口** (`ent`): 位置 $(x, y, 0)$ と向き $\vec q = (q_x, q_y, q_z)$ ($q_z = \sqrt{1 - q_x^2 - q_y^2}$。`pxy2dpz` は $q_z - 1$) を $R$ で回す:

$$
\begin{pmatrix}x_f\\y_f\\z_f\end{pmatrix} = R\begin{pmatrix}x\\y\\0\end{pmatrix},\qquad
\begin{pmatrix}q_{xf}\\q_{yf}\\q_{zf}\end{pmatrix} = R\begin{pmatrix}q_x\\q_y\\q_z\end{pmatrix}
$$

(コードの `xf = r11*xi + r12*yi`、`zf = r31*xi + r32*yi`、`pxf = r11*pxi + r12*pyi + r13*pzi`。位置の $z$ 成分が 0 なので位置には第 3 列を使わない。) 回した後の粒子は、直線 $\vec r(\lambda) = (x_f, y_f, z_f) + \lambda\vec q_f$ の上にある。新しい面 $z = 0$ との交点は $z_f + \lambda q_{zf} = 0$ より $\lambda = -z_f/q_{zf}$:

$$
x' = x_f - \frac{q_{xf}}{q_{zf}}z_f + \mathrm{DX},\quad
y' = y_f - \frac{q_{yf}}{q_{zf}}z_f + \mathrm{DY},\quad
z' = z + \frac{z_f}{q_{zf}} + \mathrm{DZ},\quad q'_{x,y} = q_{xf,yf}
$$

(`x(i) = xf - pxf/pzf*zf + dx`、`z(i) = z(i) + zf/pzf + dz`。$z$ は遅れで、面まで $\lambda$ だけ戻る分だけ $-\lambda = z_f/q_{zf}$ 変わる。) これは 11 章 (11.4) と同じ式で、(11.3) の「スナップショット」の $z\vec u$ の項は、回した後も同じ向きなので $z$ にそのまま加算されるだけになる。**出口**: 先に $(x - \mathrm{DX}, y - \mathrm{DY})$ を引いてから $R$ で回し、同じように面へ戻す (`xi = x(i) - dx`)。$z$ は $z + z_f/q_{zf} + \mathrm{DZ} - dv\,\mathrm{DZ}$ (出口の DZ だけ $dv$ をかける。K64。pysad の v2 は入口にも入れる)。

確認 ([sol_04_rotation_branches.out.txt](outputs/sol_04_rotation_branches.out.txt) (a)): Fortran の式を numpy に写し、直線と平面の交点を 3 次元で独立に解いたものと比べると、5 つの乱数入射 (角 0.05 rad 程度、$q \sim 0.03$) で差は $7 \times 10^{-18}$。

### 8.2 要素の取り付け誤差: `tsolrot` ($B_z = 0$ の枝)

`tsolrot` (`tsol.f:9`) は、`tmulti` (MULT)・`tquad` (QUAD) などが要素の前後で呼ぶ。引数 `alg` は支点までの距離 (要素の半分の長さ。入口の呼び出しで `alg`、出口の呼び出しで `alg - al` $= -$`alg`)、`dx, dy, dz`、`chi1, chi2, chi3`、`ent` (入口なら真)。MULT は `-chi1, -chi2` を渡す (`tmulti.f:91`、`:304`)。以下は直線の要素 ($B_z = 0$) の枝。

**入口** (`ent`): `s0s` $= -$`alg` (支点から見た入口の面の位置)。コードの位置の更新は次の 4 行 (と運動量の回転、最後の戻し) である:

- `ds1 = schi1*x + dcchi1*s0s - cchi1*dz`、`x = cchi1*x - schi1*(s0s - dz)` ($x$ は先に `x -= dx`)
- `ds2 = schi2*y + cchi2*ds1 + dcchi2*s0s`、`y = cchi2*y - schi2*(ds1 + s0s)` ($y$ は先に `y -= dy`)
- 運動量: `px = cchi1*px - schi1*pz0`、`py = cchi2*py - schi2*(schi1*px + cchi1*pz0)` (`px` は更新前の値、`pz0` $= q_z$)
- 戻し: `a = ds2/pz2`、`x -= a*px`、`y -= a*py`、`z += a`

位置ベクトルを、支点を原点として $\vec v = (x - \mathrm{DX},\ y - \mathrm{DY},\ s_{0s} - \mathrm{DZ})$ と書く (入口の面は $z = s_{0s}$。要素を DZ だけ前へずらすと、粒子は相対的に DZ だけ後ろ)。第 1 の回転 $R_1' = \begin{pmatrix}c_1&0&-s_1\\0&1&0\\s_1&0&c_1\end{pmatrix}$ ($x_1 = c_1x - s_1z$、$z_1 = s_1x + c_1z$) を $\vec v$ に施すと

$$
x_1 = c_1x - s_1(s_{0s} - \mathrm{DZ}),\qquad \mathrm{ds1} \equiv z_1 - s_{0s} = s_1x + (c_1 - 1)s_{0s} - c_1\mathrm{DZ}
$$

(`dcchi1` $= c_1 - 1$。) 第 2 の回転 $R_2' = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}$ を $(x_1, y, z_1)$ に施すと $y' = c_2y - s_2z_1 = c_2y - s_2(\mathrm{ds1} + s_{0s})$、$z_2 = s_2y + c_2z_1$、$\mathrm{ds2} \equiv z_2 - s_{0s} = s_2y + c_2\,\mathrm{ds1} + (c_2 - 1)s_{0s}$。これはコードの `y`、`ds1`、`ds2` と一致する。向きも同じ $R_2'R_1'$ で回る。回した後の粒子は面 $z = s_{0s}$ から $\mathrm{ds2}$ だけ前にいるので、向きに沿って $a = \mathrm{ds2}/q_{z2}$ だけ戻す。以上は、「要素の中央を支点にして要素を $(\mathrm{DX}, \mathrm{DY}, \mathrm{DZ})$ ずらし、$R_2'R_1'$ で回し、粒子の直線を新しい入口の面へ戻す」剛体運動 (11 章 §3) になっている。

MULT が渡す $\chi' = -\chi$ を代入すると $R_2'(-\chi_2)R_1'(-\chi_1) = \begin{pmatrix}1&0&0\\0&c_2&s_2\\0&-s_2&c_2\end{pmatrix}\begin{pmatrix}c_1&0&s_1\\0&1&0\\-s_1&0&c_1\end{pmatrix} = R_2R_1$ (11 章 (11.6) の $R_2$・$R_1$)。つまり MULT の CHI1・CHI2 は、COORD の $R_3R_2R_1$ ($R_3$ は `chi3` の回転で、コードでは最後に掛ける) と同じ符号の回転になっている。**ヘルプの MULT の式の $R_1^{(-)}$ (11 章 付録 A) との符号の関係は確認していない** (ヘルプの誤記か、`chi1` の符号の取り扱いを含めた読み方の違いか)。

**出口**: `s0s` $=$ (呼び出しの `alg`) $= -\mathtt{alg}_{\rm elem}$。コードは逆回転 $(R_2'R_1')^\top$ を $(x, y, -s_{0s})$ (出口の面は支点から前方 $+\mathtt{alg}_{\rm elem}$) に施し、$(+\mathrm{DX}, +\mathrm{DY}, +\mathrm{DZ})$ を足し、面へ戻す (`ds2 = ... + dz`、`x = ... + dx`)。

確認 (同 (b)): 入口と出口の Fortran の枝を numpy に写し、「支点まわりの剛体回転 + 直線を面へ戻す」の 3 次元の構成 ($R = R_2'R_1'$ を位置と向きに掛け、直線と平面の交点を求める) と比べると、差は $5 \times 10^{-17}$ (入口)、$7 \times 10^{-17}$ (出口)。

**CHI3**: `chi3` の回転 ($x, y, p_x, p_y$ を $R(\chi_3)$ で回す。11 章 (11.8)) は、入口では上の回転のあとに、出口では逆回転を先に行う。$B_z \ne 0$ のとき磁場の横成分 `bxs, bys` も同じ $R(\chi_3)$ で回す (`bxs = cchi3*bxs0 - schi3*bys`)。

**$\mathrm{CHI1} = \mathrm{CHI2} = \mathrm{DZ} = 0$ で $B_z \ne 0$ の場合** (DX・DY だけ): $x, y$ を DX、DY だけずらすと、正準運動量が $f_x = b_{zs}\,\mathrm{DY}/2$、$f_y = -b_{zs}\,\mathrm{DX}/2$ ずれる (11 章 §5.3 の導出)。

### 8.3 傾いたソレノイド ($B_z \ne 0$ かつ CHI1・CHI2・DZ ≠ 0) の入口

磁場の向き $\vec\phi$ は、軸 $\hat z = (0, 0, 1)$ に同じ回転 $R_2'R_1'$ を施したもの: 第 1 の回転で $(0, 0, 1) \to (-s_1, 0, c_1)$、第 2 の回転で $y = -s_2c_1$、$z = c_2c_1$ なので

$$
\vec\phi = \mathrm{sign}(b_z)\,(-\sin\chi_1,\ -\cos\chi_1\sin\chi_2,\ \cos\chi_1\cos\chi_2)
$$

(`phix, phiy, phiz`。$b_{xs}, b_{ys}, b_{zs} = |b_z|\vec\phi$。) 回転の前に、正準運動量を運動学的運動量に直す: `px += bz*y/(2 pr)`、`py -= bz*x/(2 pr)` (`pr` $= 1 + \delta$)。回した後の粒子は、傾いた一様磁場 $|b_z|\vec\phi$ の中で、面 $z = s_{0s}$ まで**らせん運動**する。単位ベクトル $\vec u = \vec q$ の運動方程式は、経路長 $l$ について $d\vec u/dl = b_p\,\vec u\times\vec\phi$ ($b_p = |b_z|/p$。$\vec\phi = \hat z$ なら $du_x/dl = b_pu_y$、$du_y/dl = -b_pu_x$ で、§6 の $P_x' = wP_y$、$P_y' = -wP_x$ と同じ向き)、$d\vec r/dl = \vec u$。$\vec u$ を $\vec\phi$ に平行な成分 $\vec u_l = (\vec u\cdot\vec\phi)\vec\phi$ と垂直な成分 $\vec u_t = \vec u - \vec u_l$ に分けると、$\vec u_l$ は一定、$d\vec u_t/dl = b_p\vec u_t\times\vec\phi$ は回転で、$\theta = b_pl$、$\vec u_b = \vec u_t\times\vec\phi$ として

$$
\vec u(\theta) = \vec u_l + \vec u_t\cos\theta + \vec u_b\sin\theta,\qquad
\vec r(\theta) = \vec r_0 + \frac1{b_p}\bigl[\vec u_l\,\theta + \vec u_t\sin\theta + \vec u_b\,(1 - \cos\theta)\bigr]
$$

(コードの `plx, ply, plz` $= \vec u_l$、`ptx, pty, ptz` $= \vec u_t$ ($\mathtt{ptz} = u_z - u_{lz}$)、`pbx, pby, pbz` $= \vec u_b$、`alb` $= 1/b_p$、`phi` $= \theta$。) 面の条件 $z(\theta) = s_{0s}$ は、$z$ 成分で

$$
u_{lz}\,(\theta - \sin\theta) + u_z\sin\theta + u_{bz}(1 - \cos\theta) = \mathrm{ps2} \equiv \mathrm{ds2}/\mathtt{alb}
$$

(左辺は $u_{lz}\theta + u_{tz}\sin\theta + u_{bz}(1 - \cos\theta)$ で、$u_{lz}\theta + u_{tz}\sin\theta = u_{lz}(\theta - \sin\theta) + u_z\sin\theta$ と書き換えた。コードの `s=plz*xsinphi+pz0*sinphi-pbz*dcosphi`。) をニュートン法で $\theta$ について解く (`itmax = 10`)。その後 $\vec r$、$\vec u$ の式で更新し、$z \mathrel{-}= \theta/b_p$、最後に正準運動量に戻す (`px -= bzp*y/2`、`bzp = bzs/pr`)。

**確認と疑い** ([sol_04_rotation_branches.out.txt](outputs/sol_04_rotation_branches.out.txt) (c)): Fortran の式を numpy に写し ($b_z = 1.3$ /m、$\chi_1 = 0.03$、$\chi_2 = -0.02$、DZ = 1 cm、半長 0.5 m)、一様磁場中のらせん軌道の ODE (DOP853) が面 $z = s_{0s}$ に着く点と比べた。面の条件を解いた経路長 $l$ は ODE と $2 \times 10^{-17}$ で一致する。一方、位置・運動量の更新は、`xsincos` の戻り値 `dcosphi` $= \cos\theta - 1$ をそのまま使うと ODE と $3 \times 10^{-7} \sim 1.4 \times 10^{-6}$ ずれる ($x$: $-3.2 \times 10^{-7}$、$y$: $+1.1 \times 10^{-6}$、$p_x$: $+1.4 \times 10^{-6}$)。更新だけ `dcosphi` の符号を反転して $1 - \cos\theta$ にすると $10^{-9}$ に下がる (ただし面の条件は合わなくなる)。つまり、面の条件の式と更新の式で `dcosphi` の符号の想定が逆になっている。`tsolrot` のコメントアウトされた古い式 (`dcosphi = 2.d0*sin(.5d0*phi)**2` $= 1 - \cos\theta$) が `xsincos` ($\mathtt{dc} = \cos - 1$) に置き換えられたときの符号の取り違えの可能性がある。**転記の誤りの可能性と、K64 を実行していないことから、「疑い」にとどめる** (差の大きさは $\theta^2 p_t$ の程度で、$\theta = b_pl$ が小さい通常の設定では $10^{-6}$ 以下)。pysad には `tsolrot` の傾いたソレノイド ($B_z \ne 0$) の枝は実装されていない。

### 8.4 GEO・CHI3 の決まり方 (未導出)

`tsgeo` (K64。pysad `sgeo.py`) が、GEO = 1 の SOL を基準に、領域の中の設計軌道を追い、各 SOL の CHI1・CHI2・CHI3・DX・DY・DZ を決める。この幾何は、軌道の追跡 (ドリフト・四極・SOL の写像) と座標の取り替えの繰り返しで、解析式の導出は行っていない (`docs/17_sol.md`)。**未導出**: 領域の出口の CHI3 を決める条件 (出口の後の最初の GEO = 1 の MARK で CHI3 = 0 になる) の式。

## 9. 端部 (DISFRIN)

SOL の境界で $b_z$ が $\mathrm{dbz}$ だけ変わるとき、`tsfrin` ($\exp(:f:)$、$f = -\dfrac{\mathrm{dbz}}{8p^2}p_\phi p_r$。[06_edge_fringe.md](06_edge_fringe.md) (6.10)) が、`ky_FRIN_SOL == 0` (DISFRIN = 0) のとき、`bz1 - bzs` の跳びごとに入る (`tsol.f:588`、`:753`)。線形な効果 (一様な $B_z$ の中の $p_x \mathrel{+}= b_zy/2$ の跳び) は、正準運動量の連続で自動的に入る。

**確認**: SOL + QUAD + SOL で、DISFRIN = 1 にすると、追跡の結果の変化 ($x$ で $-7.4 \times 10^{-10}$、$p_x$ で $5.2 \times 10^{-9}$) は、基準の $\exp(:f:)$ の効果と $3 \times 10^{-12}$ 以内で一致する ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt))。この確認では四極自身の DISFRIN は 1 (`ttfrin` の $B_z$ つきの枝は使っていない)。ヘルプの SOL の F1 (放射のための場の立ち上がりの長さ) は EMITTANCE の計算にだけ効く (この章の範囲外)。

## 10. ヘルプ ↔ Fortran: `tzparams` の全変数

`tzsetparams` (`tsolque.f:161-215`。$B_z = 0$ は `tzsetparam0`)。pysad の `track6._TZ`。$p = 1 + \delta$。

| 変数 | 定義 | 意味 |
|---|---|---|
| `aln` | $\ell$ | スライスの長さ |
| `pr` | $p$ | |
| `bzp`、`akkp` | $b_z/p$、$K/p$ | |
| `w1` | $\sqrt{(b_p^2 + \sqrt{b_p^4 + 4k_p^2})/2}$ (`hypot(bzp**2, 2*akkp)`) | 振動数 (9.1) |
| `w2` | $k_p/w_1$ | 増減率 |
| `wss` | $1/(w_1^2 + w_2^2)$ | $U$ |
| `ws`、`w12` | $w_1 + w_2$、$w_1 - w_2$ | $w_+$、$w_-$ |
| `wd` | $b_p/w_+$ | 振動モードの $y/x$ の比 (§4.2) |
| `phi1`、`phi2` | $\ell w_1$、$\ell w_2$ ($|{\rm phi2}| \le 100$) | |
| `s1`、`xs1`、`c1`、`dc1` | $\sin\phi_1$、$\phi_1 - \sin\phi_1$、$\cos\phi_1$、$\cos\phi_1 - 1$ | (`xsincos`。`xs1` は $\sin\phi_1/\phi_1$ ではない) |
| `sh2`、`xsh2`、`ch2`、`dch2` | $\sinh\phi_2$、$\phi_2 - \sinh\phi_2$、$\cosh\phi_2$、$\cosh\phi_2 - 1$ | (`xsincosh`) |
| `g1`、`g2` | $(s_1w_2)^2/k_p$、$-(\mathrm{sh}_2w_1)^2/k_p$ | `tztafs` の $\mathtt{g1t}$ (係数 $f_2$ の一部。$z$ の積分の項) |
| `wr1`、`wr2` | $w_1/w_2$、$w_2/w_1$ | |
| `cxs1`、`cxs2` | `dc1*s1 - xs1` $= \sin\phi_1\cos\phi_1 - \phi_1$、`dch2*sh2 - xsh2` $= \sinh\phi_2\cosh\phi_2 - \phi_2$ | $\int\sin^2$、$\int\sinh^2$ の定数でない部分 (§5.3) |
| `aw1`、`aw2` | $\ell + \mathrm{cxs2}\,\mathrm{wr1}/w_+$、$\ell + \mathrm{cxs1}\,\mathrm{wr2}/w_+$ | `tztafs` の `zw1` (積分の定数の部分) |
| `cr2`、`cr3` | $c_1w_2/w_1$、$\mathrm{ch}_2w_1/w_2$ | |
| `csw1` | $k_p\ell(-(\cos\phi_1 - 1) + (\cosh\phi_2 - 1))$ | |
| `csws` | $(\mathrm{ch}_2w_1\phi_1 + c_1w_2\phi_2)/V$ | |
| `ca1` | $-\mathrm{ch}_2\,\mathrm{dc}_1 - \mathrm{dch}_2 - s_1\mathrm{sh}_2w_-/w_+$ | |
| `dxs` | $w_2\,\mathrm{xs}_1 - w_1\,\mathrm{xsh}_2$ | |
| `dcw1`、`dcw2` | $(\cos\phi_1 - 1)\phi_2$、$(\cosh\phi_2 - 1)\phi_1$ | |

(`tzparamp` の $p$ 微分 (`w1p` など) は、`tsolque` (CALC) で $\partial z/\partial\delta$ の行を作るためのもの。この章では扱わない。) 各係数は定義を Fortran から写したもの。(9.8) の積分への割り当ては §5.3 (`aw1`・`aw2` は手で導出、全体は単項式ごとに sympy で恒等式を確認)。`g1`・`g2`・`cr2`・`cr3`・`csw1`・`csws`・`ca1`・`dcw1`・`dcw2`・`dxs` を単独の物理量 (どの積分のどの部分) に 1 つずつ対応づけることはしていない。

## 11. Fortran の対応表

| 内容 | ファイル:行 | サブルーチン | 式 |
|---|---|---|---|
| 本体 (追跡) | `tsolqu.f:1-166` | `tsolqu` | §7、(9.6)、(9.9) |
| パラメータ ($B_z \ne 0$) | `tsolque.f:161-215` | `tzsetparams` | (9.1)、§10 |
| パラメータ ($B_z = 0$) | `tsolque.f:217-239` | `tzsetparam0` | 04 章 |
| $z$ の更新 | `tsolque.f:39-86` | `tztafs` | (9.8) |
| 放射あり・本体の変形 | `tsolqu.f:167-345`、`:347-539` | `tsolqum`、`tsolqur` | 範囲外 |
| CALC の本体 | `tsolque.f:313-` | `tsolque` | |
| 座標の切替 (要素の取り付け誤差) | `tsol.f:9-252`、`:254-` | `tsolrot`、`tsolrote` | §8.2、§8.3 |
| 領域の端の座標の取り替え | `tsol.f:819-` | `trots` | §8.1 |
| ソレノイドの端 | `tsol.f:885-976`、`:978-` | `tsfrin`、`tsfrie` | 06 章 (6.10) |
| SOL 要素の領域 | `tsol.f:510-760` | `tsol` | §8、§9 |

## 12. pysad の対応表

| 内容 | 関数 | 対応する Fortran | 注意 |
|---|---|---|---|
| パラメータと $z$ の更新 | `track6._TZ` (`tztafs` はメソッド) | `tzsetparams`、`tztafs` | 同名の属性 |
| 本体 (追跡) | `track6._tsolqu_bz_ens`、`_tsolqu_bz_final`、`track6n_ext.tsolqu_bz_all` | `tsolqu` | $B_z \ne 0$。放射の有無で並びが変わる (docs/17) |
| ソレノイドの端 | `track6._tsfrin`、`track6n_ext.tsfrin1`、`sgeo._tsfrin` | `tsfrin` | |
| 領域の端の座標の取り替え | `track6._trots`、`track6n_ext.trots1` | `trots` | v2 は入口にも $dv\cdot\mathrm{DZ}$ を入れる。K64 は出口だけ (docs/17_sol.md) |
| 領域の幾何 (GEO) | `sgeo.py` (`tsgeo`、`region_frames`) | `tsgeo` | docs/17_sol.md |
| 要素の取り付け誤差 (`tsolrot`) | CHI3・DX・DY のみ (`tsolrote` 相当。`lattice.py`、`emit6d.py`) | `tsolrot`、`tsolrote` | **CHI1・CHI2・DZ は実装なし** (11 章 §7) |

## 13. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `sol_01_h2_exact.py` | 特性方程式と $w_1$、$w_2$ (sympy)、ヘルプの式・K64 (`_TZ`) と行列指数の一致、$z$ の更新、$\Delta H$ の解、固有モードと `a`、`bw`、`c`、`dw` | [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) |
| `sol_02_chain.py` | SOL + QUAD + SOL の追跡と、独立な基準 (端の流れ + 全 $H$ の ODE) の比較、DISFRIN | [sol_02_chain.out.txt](outputs/sol_02_chain.out.txt) |
| `sol_03_derivations.py` | $H_2$ の展開、$(A, B, C, D)$ と `a`・`bw`・`c`・`dw`、$z$ の積分 (積分の表、`tztafs` の恒等式)、`tsolqu` の $\Delta H$ の回転の行列式、SOL + QUAD + SOL のヤコビアンのシンプレクティック性の欠陥 | [sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) |
| `sol_04_rotation_branches.py` | `trots`、`tsolrot` の枝と 3 次元の構成・らせん軌道の ODE の比較 | [sol_04_rotation_branches.out.txt](outputs/sol_04_rotation_branches.out.txt) |

端のフリンジの確認は [edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (06 章)。

## 14. 未確認・限界・使った条件

- **使った条件**: $p_0$ = 3 GeV の電子 ($B\rho = 10.007$ T m)。$b_z = 0.3 \sim 2.5$ /m、$K = 0.9 \sim 4$ /m²、$\ell = 0.4$ m、振幅 mm、角 mrad。$\delta = \pm0.02$。$H_2$ の確認は 4×4 の行列指数 (正準変数)、$z$・$\Delta H$ は DOP853。
- **SAD (gs) との比較はしていない**。ヘルプの式は、K64 版同梱の原文 `Documents/SADHelp.html` の `solution-H2`・`solution-dH`・SOL の `DISFRIN` を直接読んで照合した (付録 A)。
- **確認できなかった**: (i) `tztafs` の係数のうち `aw1`・`aw2` 以外を、積分の部分ごとに 1 つずつ手で割り当てること (全体の代数的な一致は確認済み)。(ii) `tsgeo` の幾何 (領域の出口の CHI3 を決める条件など) の式 (§8.4)。`tsolrot` の傾いたソレノイドの枝の `dcosphi` の符号 (§8.3。疑い)。(iii) 放射あり (`tsolqur`)・スピン・CALC (`tsolque` の解析微分)。(iv) 四極自身のフリンジ (`ttfrin` の $B_z$ つきの枝) と SOL の端の合成。(v) $\Delta H$ の $K \ne 0$ との相互作用は分割の誤差 (Strang の 2 次) に任せている。分割数 $n$ の式の係数 ($\varepsilon = 0.2$) の根拠は、ソースに無く導けない (§7)。(vi) `tsolqu` の $\Delta H$ の回転の更新順序 (付録 B) の gs での確認。
- **範囲**: 曲率 $1/\rho = 0$ の領域。$b_z$ は一様。軸上の場と、軸対称の端の場の 3 次までを扱う。

## 15. この章で使う Fortran の変数名

§2、§10 の表を参照 (すべての変数に定義と式を付けた)。補足: `ap` $= q_x^2 + q_y^2$ (`px**2 + py**2`)、`dpz` $= \sqrt{1 - \mathtt{ap}} - 1$ (`sqrt1(-ap)` は桁落ちしない形)、`awu` $= a\,w_1/w_+$ ($= A\,w_1/w_+$、振動モードの振幅のうち位置から決まる分。(9.8a))、`dwu` $= dw\cdot w_2$ (指数モードの同じ量)、`a12`・`a14`・`a22`・`a24` は回転の `xsincos` の戻り値 (角 $r\,b_p$)。

## 付録 A. ヘルプの記述・K64 との違い (ソレノイド)

| 箇所 | 内容 |
|---|---|
| 原文との照合 (`SADHelp.html`) | `solution-H2` ($w_1$、$w_2$、$w_+$、$U$、$V$、$a, b, c, d$、$u_{1,2}$、$v_{1,2}$、$\Delta x$、$\Delta y$、$p_x$、$p_y$、$H_{2u}$、$z$) と `solution-dH` を原文で読み、$H_2$ の式は行列指数と一致を確認 (§4.3)。`solution-dH` は下の 2 行の食い違い。SOL の `DISFRIN` の $H$ は原文どおり |
| ヘルプ solution-dH の $\Delta x$、$\Delta y$ | $(p_{x0} + \tfrac{B_z}2y_0)\sin w\ell - (p_{y0} - \tfrac{B_z}2x_0)(\cos w\ell - 1)$ の式に **$1/B_z$ が無い** (**原文 `SADHelp.html` の `solution-dH` で確認**。TeX・平文とも)。(9.9) の導出と ODE ($3.6 \times 10^{-17}$)、K64 の `(a24*pxi - a14*py)/bzp` から、$1/B_z$ が要る (無いと $1.6 \times 10^{-6}$ ずれる)。ヘルプの誤記の疑い |
| ヘルプ solution-dH の $\Delta H$ | 式に `a_c` とあるが $a_x$ の誤記とみられる |
| ヘルプ solution-H2 | $K > 0$ にそろえる (回転) と書く。$H_2$ の式は行列指数と $10^{-18}$ で一致 |
| ヘルプの SOL の DISFRIN の $H$ | paraxial ($p^2$)。06 章 付録 A |
| ヘルプに無いもの | (1) 特性方程式の導出と $w_1$・$w_2$ の関係 (§4)。(2) 固有モード $y/x = \pm iw_d$、$\pm1/w_d$。(3) $z$ の更新が $-\tfrac12\int q^2$ であること (§5)。(4) ラーモア回転としての $\Delta H$ の解 (§6)。(5) `tzparams` の変数の定義 (§10) |

## 付録 B. K64 の誤りの疑い (別枠。確定には SAD (gs) での確認が要る)

この章の範囲 ($H_2$ の厳密解、$z$ の更新、$\Delta H$ の解、Strang 分割、SOL + QUAD + SOL の全体、DISFRIN) で、K64 の式は行列指数・ODE・全 $H$ の流れと、分割の誤差の範囲で一致する (§4〜§7)。ただし次の 2 点は**疑い**で、gs での確認はしていない。ヘルプの誤記の疑いは付録 A。

1. **`tsolqu` の $\Delta H$ の運動量の回転がシンプレクティックでない疑い** (ソースの読みと pysad (v2 既定、K64 のまま写してある部分) のヤコビアンの数値)。`tsolqu.f:122-123` (`tsolqum` の `:302-303`、`tsolqur` の `:487-488` も同じ) は

   - `pxi = a22*pxi + a24*py(i)`
   - `pyi = -a24*pxi + a22*py(i)`

   と書く。2 行目の `pxi` は、1 行目で更新した後の値である。正しい回転は $p_y' = -a_{24}\,p_{x,\rm old} + a_{22}p_y$。実際の写像は、$(p_{x,\rm old}, p_y) \to (p_x', p_y')$ の行列

   $$
   \begin{pmatrix}p_x'\\p_y'\end{pmatrix} = \begin{pmatrix}a_{22} & a_{24}\\ -a_{22}a_{24} & a_{22} - a_{24}^2\end{pmatrix}\begin{pmatrix}p_{x,\rm old}\\p_y\end{pmatrix},\qquad
   \det = a_{22}(a_{22} - a_{24}^2) + a_{22}a_{24}^2 = a_{22}^2 = \cos^2 w\ell
   $$

   で、正しい回転 (行列式 $a_{22}^2 + a_{24}^2 = 1$) と違い、面積を $\cos^2 w\ell$ 倍する ([sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (d)。sympy)。位置の更新 (`x`、`y`) は更新前の `pxi` を使っていて正しい。ループの後の最後の半ステップ (`tsolqu.f:151-156`) の運動量の更新も更新前の値を使っていて正しい。pysad は `track6.py:2517`・`track6n_ext.py:480` に「SAD のまま: 更新後の pxi を使う」と写している (k64・v2 とも。docs/30 に B 番号は無い)。

   **数値** (同 (e)。pysad v2 既定、SOL (BOUND = 1、$b_z$) + QUAD ($K = 0.36$ /m²、$L = 0.4$ m) + SOL ($B_z = 0$)、3 GeV、$\delta = 0$、$x = 1$ mm、$y = -0.5$ mm、中心差分の刻み $10^{-6}$ の 6×6 ヤコビアン $J$ で $\|J^\top SJ - S\|_{\max}$):

   | $b_z$ [/m] | $p_x = 0.003$ | $p_x = 0.05$ | $p_x = 0.1$ | $p_x = 0.1$、$p_y = 0.05$ |
   |---|---|---|---|---|
   | 1.3 | $1.9 \times 10^{-12}$ | $1.2 \times 10^{-7}$ | $1.9 \times 10^{-6}$ | $3.4 \times 10^{-6}$ |
   | 5 | $6.5 \times 10^{-11}$ | $2.4 \times 10^{-6}$ | $4.0 \times 10^{-5}$ | $3.2 \times 10^{-5}$ |
   | 10 | $7.3 \times 10^{-10}$ | $6.7 \times 10^{-6}$ | $1.1 \times 10^{-4}$ | $1.0 \times 10^{-4}$ |

   傾き $p_x = 3$ mrad の試験 (§7) の領域では $10^{-12} \sim 10^{-9}$ で見えないが、傾き 0.05 で $b_z = 1.3$ /m でも $1.2 \times 10^{-7}$ (ヤコビアンの行列式は $1 - 1.3 \times 10^{-7}$)、強い $B_z$ や大きな傾きで大きくなる。$\det J - 1 = -3.2 \times 10^{-6}$ ($b_z = 1.3$、$p_x = 0.1$、$p_y = 0.05$) は、予想される $\prod\cos^2(w\ell_i) - 1$ (各スライスの回転角の 2 乗の和程度) と同じ桁である。欠陥の大きさは $O\bigl((b_pq^2\ell)^2\bigr)$ ($w\ell \sim b_pq^2\ell/2$、$q = $ 傾き)。**確認の状況**: ソースの読み (確定) と pysad のヤコビアンの数値 (K64 の式を写した部分の挙動) まで。gs (K64) での確認と、v2 で更新前の `pxi` を使う修正をするかどうかは、別に進めていて結論が出ていないため、結論が出るまで「疑い」とする。

2. **`tsolrot` の傾いたソレノイドの枝 ($B_z \ne 0$、CHI1・CHI2・DZ ≠ 0) の `dcosphi` の符号の想定が、面の条件の式と更新の式で逆** (§8.3。転記の誤りの可能性があり、K64 を実行していない)。差は $\theta^2p_t$ の程度。
