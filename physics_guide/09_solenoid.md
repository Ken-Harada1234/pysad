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
4. $\Delta H = p - p_z - P^2/(2p)$ の解は**ラーモア回転**: 運動学的運動量 $P$ が角 $w\ell$ ($w = b_z(p - p_z)/(pp_z)$) だけ回り、位置はその間の速度 $P(1/p_z - 1/p)$ で動く。ヘルプの $\Delta x$、$\Delta y$ の式には $1/B_z$ の因子が要る (ヘルプの式のままでは ODE と $1.6 \times 10^{-6}$ ずれ、$1/B_z$ をつけると $3.6 \times 10^{-17}$。付録 A)。
5. K64 の本体は $\Delta H$ の半ステップ + ($H_2$ の厳密解、$\Delta H$) × $n$ の **Strang 分割 (2 次)**。SOL + QUAD + SOL (3 GeV、$b_z = 1.3$ /m、$K = 0.9$ /m²) の追跡は、独立な基準 (端の厳密な流れ + 全 $H$ の ODE) と $10^{-11}$ で一致 (分割の誤差)。DISFRIN = 1 でフリンジの効果 ($10^{-9}$) が消える ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt))。
6. `tztafs` の係数 (`aw1`、`ca1`、`cr2` など) は、(9.8) の積分を 4 つのモードで閉じた形にしたものと、10 個の単項式のすべてで、手で導いた閉じた形と一致する (§5.3。10 個とも手で導出し、sympy・数値積分で検算)。領域と座標の切替は、`trots` (領域の端) と `tsolrot` ($B_z = 0$ の枝) が剛体回転 + 面への戻しと一致することを確認した (§8)。確認できなかったこと: `tsgeo` の幾何の式、放射 (`tsolqur`)・スピンの扱い。§14。
7. **K64 の誤り** (付録 B): (a) `tsolqu` の $\Delta H$ の運動量の回転が、更新後の `pxi` を使うため、行列式が $\cos^2 w\ell$ になりシンプレクティックでない (gs のヤコビアンで確認、v2 は修正 = B333)。(b) `tsolrot` の傾いた枝 ($B_z \ne 0$、CHI1・CHI2・DZ ≠ 0) の更新の `dcosphi` の符号が、面の条件の式と逆 (gs で確認。pysad にはこの枝が無い = B337)。

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

**$K = 0$ の場合の検算** (回転系と固定系)。$H_2$ の流れで座標がどう動くかを、$\zeta \equiv x + iy$ (複素数) で書く。

1. $K = 0$ の (9.3a) は $H_2 = H_{\rm osc} + \omega_Lp_\phi$、$H_{\rm osc} = \dfrac{p_x^2 + p_y^2}{2p} + \dfrac{p}{2}\Bigl(\dfrac{b_p}{2}\Bigr)^2(x^2 + y^2)$、$\omega_L = -\dfrac{b_z}{2p} = -\dfrac{b_p}2$ ($-p$ の定数は省く。$b_z^2/(8p) = (p/2)(b_p/2)^2$ を使った)。
2. $p_\phi = xp_y - yp_x$ は $H_{\rm osc}$ (軸対称) と交換するので、$H_2$ の流れは「$H_{\rm osc}$ の流れ」と「$H = \omega_Lp_\phi$ の流れ」の積に分かれる。$H = \omega_Lp_\phi$ の流れは $\phi' = \partial H/\partial p_\phi = \omega_L$ ($\phi$ は $\zeta$ の偏角)、すなわち $\zeta \to e^{i\omega_Ls}\zeta = e^{-ib_ps/2}\zeta$ ($x, y$ を一定の角速度で回す)。
3. 回転系の座標 $\eta$ ($\zeta = e^{-ib_ps/2}\eta$) では $H_{\rm osc}$ だけが残り、これは等方な調和振動子 (質量 $p$、角振動数 $\Omega = b_p/2$) なので、$\eta'' = -(b_p/2)^2\eta$、解は $\eta = \alpha e^{ib_ps/2} + \beta e^{-ib_ps/2}$ ($\alpha, \beta$ は複素定数)。
4. 固定系に戻す: $\zeta = e^{-ib_ps/2}\eta = \alpha + \beta e^{-ib_ps}$。振動数は **$0$ と $b_p$** (1 つ目は一定、2 つ目は角速度 $b_p$ の回転)。「$b_p/2 \pm b_p/2$」の意味は、回転系の振動数 $b_p/2$ に、座標系の回転 $b_p/2$ を足し引きした $b_p$ と $0$ である。

**直接の確認**: (9.4) で $K = 0$ ($k_p = 0$) とすると $x'' = b_py'$、$y'' = -b_px'$ から $\zeta'' = b_py' - ib_px' = -ib_p(x' + iy') = -ib_p\zeta'$。$\zeta = \alpha + \beta e^{-ib_ps}$ はこの式を満たす ($\zeta' = -ib_p\beta e^{-ib_ps}$、$\zeta'' = -b_p^2\beta e^{-ib_ps} = -ib_p\zeta'$)。回転系でも $\eta = e^{ib_ps/2}\zeta$ が $\eta'' + (b_p/2)^2\eta = 0$ を満たす (sympy: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (B)。差 0)。これは (9.1) で $K = 0$ とした $w_1 = \sqrt{\tfrac12(b_p^2 + b_p^2)} = b_p$、$w_2 = 0$ と合う (固定系の振動数は $w_1 = b_p$ と 0。$b_p/2$ はラーモア系の振動数で、$w_1$ ではない)。$K \ne 0$ では四極が回転系で時間に依存するので、この単純な回転では解けず、次節の正規モードを使う。

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

($w_dw_+ = b_p$。) **$C$**: $C = x_0 - A$ に上の $A$ を入れ、通分する:

$$
C = \frac{V\,x_0 - w_2w_+x_0 + b_pq_{ym}}{V} = \frac{(V - w_2w_+)\,x_0 + b_pq_{ym}}{V}
$$

$V - w_2w_+ = w_1^2 + w_2^2 - w_1w_2 - w_2^2 = w_1(w_1 - w_2) = w_1w_-$。さらに $w_- = w_d^2w_+$ ($w_d^2 = w_-/w_+$) と $w_dw_+ = b_p$ から $w_1w_- = w_1w_d\,(w_dw_+) = b_p\,w_1w_d$ なので

$$
C = \frac{b_p\,(w_1w_d\,x_0 + q_{ym})}{V} \equiv b_p\,c,\qquad c = \frac{w_1w_d\,x_0 + q_{ym}}{V}
$$

**$B$ と $D$**: 第 2 式と第 3 式は $B$、$D$ だけの 2 元連立:

$$
\begin{pmatrix}w_d & 1/w_d\\ w_1 & w_2\end{pmatrix}\begin{pmatrix}B\\D\end{pmatrix} = \begin{pmatrix}y_0\\q_{xm}\end{pmatrix}
$$

係数行列の行列式は $\Delta = w_dw_2 - w_1/w_d = (w_d^2w_2 - w_1)/w_d$。$w_d^2w_2 - w_1 = (w_-w_2 - w_+w_1)/w_+ = (w_1w_2 - w_2^2 - w_1^2 - w_1w_2)/w_+ = -V/w_+$ なので $\Delta = -V/(w_dw_+) = -V/b_p$。Cramer の公式で

$$
\begin{aligned}
B &= \frac{y_0\,w_2 - q_{xm}/w_d}{\Delta} = \frac{b_p}{V}\Bigl(\frac{q_{xm}}{w_d} - w_2y_0\Bigr) = \frac{w_+q_{xm} - b_pw_2y_0}{V} \qquad (b_p/w_d = w_+)\\
D &= \frac{w_d\,q_{xm} - w_1y_0}{\Delta} = \frac{b_p\,(-w_dq_{xm} + w_1y_0)}{V}
\end{aligned}
$$

(Cramer の公式: $B$ は、係数行列の第 1 列 $(w_d, w_1)^\top$ を右辺 $(y_0, q_{xm})^\top$ に取り替えた行列式 $y_0w_2 - q_{xm}/w_d$ を $\Delta$ で割ったもの。$D$ は第 2 列を取り替えた行列式 $w_dq_{xm} - w_1y_0$ を $\Delta$ で割ったもの。$\Delta = -V/b_p$ なので $1/\Delta = -b_p/V$。) これらを K64 の `a`、`bw`、`c`、`dw` ($a = A$、$bw = B$、$c = C/b_p$、$dw = D/b_p$) と並べると

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

(上の中間式 (行列式 $-V/b_p$、$V - w_2w_+ = w_1w_-$、$w_d^2w_1 + w_2 = V/w_+$ など) も sympy で 0 を確認: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (A)。sympy で 4 元連立方程式を解いた結果と上の 4 式の差が恒等的に 0: [sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (b)。$w_d^2 = w_-/w_+$ を使って簡約した。4×4 を数値で解いた $(A, B, C, D)$ と K64 の式は $4 \times 10^{-19}$ で一致: [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (e)。)

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

**交差項 4 つの積分を、原始関数の微分で 1 つずつ確かめる** ($S = \sin w_1s$、$K = \cos w_1s$ (係数 $C$ と区別)、$\mathrm{Sh} = \sinh w_2s$、$\mathrm{Ch} = \cosh w_2s$。$S' = w_1K$、$K' = -w_1S$、$\mathrm{Sh}' = w_2\mathrm{Ch}$、$\mathrm{Ch}' = w_2\mathrm{Sh}$)。原始関数 $F(s)$ について $\int_0^\ell = F(\ell) - F(0)$ で、$F(0)$ が下限の定数項:

1. $F_{K\mathrm{Ch}} = (w_1S\,\mathrm{Ch} + w_2K\,\mathrm{Sh})/V$:
   $\dfrac{d}{ds}(w_1S\,\mathrm{Ch}) = w_1^2K\,\mathrm{Ch} + w_1w_2S\,\mathrm{Sh}$、$\dfrac{d}{ds}(w_2K\,\mathrm{Sh}) = -w_1w_2S\,\mathrm{Sh} + w_2^2K\,\mathrm{Ch}$。和は $(w_1^2 + w_2^2)\,K\,\mathrm{Ch} = V\,K\,\mathrm{Ch}$ ($S\,\mathrm{Sh}$ の項が打ち消す)。$F(0) = 0$ (表に定数項なし)。
2. $F_{K\mathrm{Sh}} = (w_1S\,\mathrm{Sh} + w_2K\,\mathrm{Ch})/V$:
   $\dfrac{d}{ds}(w_1S\,\mathrm{Sh}) = w_1^2K\,\mathrm{Sh} + w_1w_2S\,\mathrm{Ch}$、$\dfrac{d}{ds}(w_2K\,\mathrm{Ch}) = -w_1w_2S\,\mathrm{Ch} + w_2^2K\,\mathrm{Sh}$。和は $V\,K\,\mathrm{Sh}$。$F(0) = w_2/V$ なので $\int_0^\ell = F(\ell) - w_2/V$ (表の「$-w_2$」)。
3. $F_{S\mathrm{Ch}} = (-w_1K\,\mathrm{Ch} + w_2S\,\mathrm{Sh})/V$:
   $\dfrac{d}{ds}(-w_1K\,\mathrm{Ch}) = w_1^2S\,\mathrm{Ch} - w_1w_2K\,\mathrm{Sh}$、$\dfrac{d}{ds}(w_2S\,\mathrm{Sh}) = w_1w_2K\,\mathrm{Sh} + w_2^2S\,\mathrm{Ch}$。和は $V\,S\,\mathrm{Ch}$。$F(0) = -w_1/V$ なので $\int_0^\ell = F(\ell) + w_1/V$ (表の「$+w_1$」)。
4. $F_{S\mathrm{Sh}} = (-w_1K\,\mathrm{Sh} + w_2S\,\mathrm{Ch})/V$:
   $\dfrac{d}{ds}(-w_1K\,\mathrm{Sh}) = w_1^2S\,\mathrm{Sh} - w_1w_2K\,\mathrm{Ch}$、$\dfrac{d}{ds}(w_2S\,\mathrm{Ch}) = w_1w_2K\,\mathrm{Ch} + w_2^2S\,\mathrm{Sh}$。和は $V\,S\,\mathrm{Sh}$。$F(0) = 0$。

(sympy で 4 つとも $dF/ds - f = 0$、$F(0)$ は上のとおり: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (D)。)

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

**$\mathtt{pxi}^2$ と $\mathtt{pyi}^2$ の係数** ((9.8a) を直接入れる最も簡単な 2 つ)。$\mathtt{awu} = \mathtt{dwu} = \mathtt{pyi} = 0$ とすると (9.8a) から $(A, B, C, D) = (0, \mathtt{pxi}/w_1, 0, 0)$、つまり振動モードの 1 つだけ ($B$ モード)。このとき

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

**残りの 8 つの係数 (すべて手で導く)**。方針は、(1) $x'$、$y'$ を $(\mathtt{awu}, \mathtt{dwu}, \mathtt{pxi}, \mathtt{pyi})$ で書く、(2) $-\tfrac12\int(x'^2 + y'^2)\,ds$ の各単項式の係数を、§5.2 の積分の表で閉じた形にする、(3) K64 の `tztafs` の式を展開して同じ形になることを示す、の 3 段である。以下、$S = \sin w_1s$、$K = \cos w_1s$、$\mathrm{Sh} = \sinh w_2s$、$\mathrm{Ch} = \cosh w_2s$ ($s_1 = S(\ell)$、$c_1 = K(\ell)$、$\mathrm{sh}_2 = \mathrm{Sh}(\ell)$、$\mathrm{ch}_2 = \mathrm{Ch}(\ell)$)。$w_d$ の関係 $w_dw_+ = b_p$、$w_db_p = w_d^2w_+ = w_-$、$w_-/w_d = b_p$、$b_p/w_d = w_+$、$b_p^2 = w_+w_-$ を使う。

**(1) $x'$、$y'$ の書き換え**。§5.1 の $x' = w_1(-AS + BK) + w_2(C\,\mathrm{Sh} + D\,\mathrm{Ch})$ に (9.8a) を入れる。$w_1A = w_+\mathtt{awu}$、$w_1B = \mathtt{pxi} - b_p\mathtt{dwu}$、$w_2C = w_-\mathtt{awu} + w_d\mathtt{pyi}$、$w_2D = b_p\mathtt{dwu}$ なので

$$
x' = \mathtt{awu}\,(-w_+S + w_-\mathrm{Sh}) + b_p\,\mathtt{dwu}\,(\mathrm{Ch} - K) + \mathtt{pxi}\,K + w_d\,\mathtt{pyi}\,\mathrm{Sh}
$$

$y' = w_dw_1(-AK - BS) + (w_2/w_d)(C\,\mathrm{Ch} + D\,\mathrm{Sh})$ では、$w_dw_1A = w_dw_+\mathtt{awu} = b_p\mathtt{awu}$、$w_dw_1B = w_d(\mathtt{pxi} - b_p\mathtt{dwu})$、$(w_2/w_d)C = (w_-/w_d)\mathtt{awu} + \mathtt{pyi} = b_p\mathtt{awu} + \mathtt{pyi}$、$(w_2/w_d)D = (b_p/w_d)\mathtt{dwu} = w_+\mathtt{dwu}$。$w_db_p = w_-$ を使って

$$
y' = \mathtt{awu}\,b_p\,(\mathrm{Ch} - K) + \mathtt{dwu}\,(w_-S + w_+\mathrm{Sh}) - w_d\,\mathtt{pxi}\,S + \mathtt{pyi}\,\mathrm{Ch}
\tag{9.8b}
$$

($s = 0$ で $x' = \mathtt{pxi}$、$y' = \mathtt{pyi}$: $\mathrm{Ch} - K = 0$、$S = \mathrm{Sh} = 0$、$K = \mathrm{Ch} = 1$。sympy・数値で (9.8a) を (9.6) の微分に入れた値と一致: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (C)。)

**(2) 各単項式の被積分関数**。$\Delta z = -\tfrac12\int_0^\ell(x'^2 + y'^2)\,ds$ に (9.8b) を入れると、$x'^2 + y'^2$ の各単項式の係数は ($\mathtt{awu}$ と $\mathtt{dwu}$ に着目して整理):

| 単項式 | $x'^2 + y'^2$ での係数 |
|---|---|
| $\mathtt{awu}^2$ | $(-w_+S + w_-\mathrm{Sh})^2 + b_p^2(\mathrm{Ch} - K)^2$ |
| $\mathtt{dwu}^2$ | $b_p^2(\mathrm{Ch} - K)^2 + (w_-S + w_+\mathrm{Sh})^2$ |
| $\mathtt{awu}\,\mathtt{dwu}$ | $2b_p(\mathrm{Ch} - K)\bigl[(-w_+S + w_-\mathrm{Sh}) + (w_-S + w_+\mathrm{Sh})\bigr] = 4b_p(\mathrm{Ch} - K)(w_1\mathrm{Sh} - w_2S)$ |
| $\mathtt{awu}\,\mathtt{pxi}$ | $2K(-w_+S + w_-\mathrm{Sh}) - 2w_-(\mathrm{Ch} - K)S$ ($w_db_p = w_-$) |
| $\mathtt{awu}\,\mathtt{pyi}$ | $2w_d\mathrm{Sh}(-w_+S + w_-\mathrm{Sh}) + 2b_p(\mathrm{Ch} - K)\mathrm{Ch}$ |
| $\mathtt{dwu}\,\mathtt{pxi}$ | $2b_pK(\mathrm{Ch} - K) - 2w_dS(w_-S + w_+\mathrm{Sh})$ |
| $\mathtt{dwu}\,\mathtt{pyi}$ | $2w_-(\mathrm{Ch} - K)\mathrm{Sh} + 2\,\mathrm{Ch}\,(w_-S + w_+\mathrm{Sh})$ |
| $\mathtt{pxi}\,\mathtt{pyi}$ | $2w_dK\,\mathrm{Sh} - 2w_dS\,\mathrm{Ch} = 2w_d(K\,\mathrm{Sh} - S\,\mathrm{Ch})$ |
| $\mathtt{pxi}^2$、$\mathtt{pyi}^2$ | $K^2 + w_d^2S^2$、$w_d^2\mathrm{Sh}^2 + \mathrm{Ch}^2$ (上の $\mathtt{aw2}$、$\mathtt{aw1}$) |

($\mathtt{awu}\,\mathtt{dwu}$ の行: $w_- - w_+ = -2w_2$、$w_- + w_+ = 2w_1$ から $(-2w_2S + 2w_1\mathrm{Sh})$。)

**記法**: $I_{XY} \equiv \int_0^\ell XY\,ds$ ($X, Y \in \{S, K, \mathrm{Sh}, \mathrm{Ch}\}$。例: $I_{S^2} = \int S^2ds$)。値は §5.2 の表。

**補助の積分**。§5.2 の表から、2 つの組み合わせが何度も出る:

$$
J_1 \equiv \int_0^\ell(K\,\mathrm{Sh} - S\,\mathrm{Ch})\,ds = \frac{(w_1s_1\mathrm{sh}_2 + w_2c_1\mathrm{ch}_2 - w_2) - (-w_1c_1\mathrm{ch}_2 + w_2s_1\mathrm{sh}_2 + w_1)}{V} = \frac{w_-s_1\mathrm{sh}_2 + w_+(c_1\mathrm{ch}_2 - 1)}{V} \equiv \frac EV
$$

$$
J_2 \equiv \int_0^\ell(K\,\mathrm{Ch} - S\,\mathrm{Sh})\,ds = \frac{(w_1s_1\mathrm{ch}_2 + w_2c_1\mathrm{sh}_2) - (-w_1c_1\mathrm{sh}_2 + w_2s_1\mathrm{ch}_2)}{V} = \frac{w_-s_1\mathrm{ch}_2 + w_+c_1\mathrm{sh}_2}{V}
$$

($w_1 - w_2 = w_-$、$w_1 + w_2 = w_+$ でまとめた。) K64 の `ca1` $= -\mathrm{ch}_2\,\mathrm{dc}_1 - \mathrm{dch}_2 - s_1\mathrm{sh}_2w_-/w_+$ は、$\mathrm{dc}_1 = c_1 - 1$、$\mathrm{dch}_2 = \mathrm{ch}_2 - 1$ を入れると $-\mathrm{ch}_2c_1 + \mathrm{ch}_2 - \mathrm{ch}_2 + 1 - s_1\mathrm{sh}_2w_-/w_+ = 1 - c_1\mathrm{ch}_2 - s_1\mathrm{sh}_2w_-/w_+$ なので、$E = -w_+\,\mathtt{ca1}$。

**(3) K64 の `tztafs` の展開**。`tztafs` の中身 (`tsolque.f:59-71`。引数 $\mathtt{zwu}, \mathtt{pxt}, \mathtt{pyt}, \mathtt{zw1}, \mathtt{wst}, \mathtt{w12t}, \mathtt{swss}, \mathtt{g1t}$) は

$$
\begin{aligned}
&\mathtt{cwsd} = c_1\,\mathtt{w12t}\,\mathtt{xsh2} - \mathrm{ch}_2\,\mathtt{wst}\,\mathtt{xs1},\quad
\mathtt{cr1} = b_p^2\,\mathtt{swss},\quad \mathtt{wsm} = 4\,\mathtt{swss}\,k_p,\quad
\mathtt{ws1} = (-1 + \mathtt{wsm})\,\mathtt{wst},\quad \mathtt{ws2} = \mathtt{w12t}\,(1 + \mathtt{wsm}),\\
&f_1 = (\mathtt{csw1} + \mathtt{cwsd})\,\mathtt{swss} + \mathtt{csws} - \mathtt{zw1},\quad
f_2 = \mathtt{ca1}\,\mathtt{cr1} + \mathtt{g1t},\quad
f_3 = 2\,\mathtt{cr1}\,\mathtt{cwsd} + \mathtt{cr3}\,\mathtt{w12t}\,\mathtt{xsh2} + \mathtt{dcw1}\,\mathtt{ws1} + \mathtt{dcw2}\,\mathtt{ws2} - \mathtt{cr2}\,\mathtt{wst}\,\mathtt{xs1}
\end{aligned}
$$

(`xs1` $= \phi_1 - s_1$、`xsh2` $= \phi_2 - \mathrm{sh}_2$。`csw1`・`csws`・`cr2`・`cr3`・`dcw1`・`dcw2`・`g1` は §10 の表。) 第 1 の呼び出しは $(\mathtt{awu}, \mathtt{pxi}, \mathtt{pyi}, \mathtt{aw1}, w_+, w_-, 1/V, \mathtt g_1)$、第 2 の呼び出しは $(-\mathtt{dwu}, -\mathtt{pyi}, \mathtt{pxi}, \mathtt{aw2}, -w_-, w_+, -1/V, \mathtt g_2)$ なので、`tztafs` の戻り値 $\mathtt{zwu}\,(\mathtt{zwu}f_3/2 + f_2\mathtt{pxt} + b_pf_1\mathtt{pyt}) - \mathtt{zw1}\,\mathtt{pyt}^2/2$ から、単項式ごとの寄与は (第 2 の呼び出しは $\mathtt{zwu} = -\mathtt{dwu}$、$\mathtt{pxt} = -\mathtt{pyi}$、$\mathtt{pyt} = \mathtt{pxi}$)

$$
\begin{aligned}
&\text{第 1:}\quad \mathtt{awu}^2\cdot\tfrac12f_3^{(1)},\quad \mathtt{awu}\,\mathtt{pxi}\cdot f_2^{(1)},\quad \mathtt{awu}\,\mathtt{pyi}\cdot b_pf_1^{(1)},\quad \mathtt{pyi}^2\cdot(-\tfrac12\mathtt{aw1}),\\
&\text{第 2:}\quad \mathtt{dwu}^2\cdot\tfrac12f_3^{(2)},\quad \mathtt{dwu}\,\mathtt{pyi}\cdot f_2^{(2)},\quad \mathtt{dwu}\,\mathtt{pxi}\cdot(-b_pf_1^{(2)}),\quad \mathtt{pxi}^2\cdot(-\tfrac12\mathtt{aw2})
\end{aligned}
$$

以下、(2) の積分と (3) の展開を、単項式ごとに突き合わせる。

**$\mathtt{awu}\,\mathtt{dwu}$**。係数 $-\tfrac12\cdot4b_p\int(\mathrm{Ch} - K)(w_1\mathrm{Sh} - w_2S)\,ds = -2b_p\bigl[w_1I_{\mathrm{Ch\,Sh}} - w_2I_{S\mathrm{Ch}} - w_1I_{K\mathrm{Sh}} + w_2I_{KS}\bigr]$。§5.2 の表を入れる: $I_{\mathrm{Ch\,Sh}} = \mathrm{sh}_2^2/(2w_2)$、$I_{KS} = s_1^2/(2w_1)$、$I_{S\mathrm{Ch}} = (-w_1c_1\mathrm{ch}_2 + w_2s_1\mathrm{sh}_2 + w_1)/V$、$I_{K\mathrm{Sh}} = (w_1s_1\mathrm{sh}_2 + w_2c_1\mathrm{ch}_2 - w_2)/V$。$1/V$ の項の和は

$$
\frac{-w_2(-w_1c_1\mathrm{ch}_2 + w_2s_1\mathrm{sh}_2 + w_1) - w_1(w_1s_1\mathrm{sh}_2 + w_2c_1\mathrm{ch}_2 - w_2)}{V}
= \frac{-(w_1^2 + w_2^2)\,s_1\mathrm{sh}_2}{V} = -s_1\mathrm{sh}_2
$$

($w_1w_2c_1\mathrm{ch}_2$ と $\pm w_1w_2$ が打ち消す。) よって

$$
-2b_p\Bigl[\frac{w_1\mathrm{sh}_2^2}{2w_2} + \frac{w_2s_1^2}{2w_1} - s_1\mathrm{sh}_2\Bigr]
= -\frac{b_p}{w_1w_2}\bigl[w_1^2\mathrm{sh}_2^2 + w_2^2s_1^2 - 2w_1w_2s_1\mathrm{sh}_2\bigr]
= -\frac{b_p}{k_p}\,(w_1\mathrm{sh}_2 - w_2s_1)^2
$$

K64 の `dxs` $= w_2\,\mathrm{xs}_1 - w_1\,\mathrm{xsh}_2 = w_2(\phi_1 - s_1) - w_1(\phi_2 - \mathrm{sh}_2) = w_1\mathrm{sh}_2 - w_2s_1$ ($w_2\phi_1 = w_1\phi_2 = w_1w_2\ell$ が打ち消す) なので、これは `bzp*(-awu*dwu*dxs**2/akkp)` と**一致**する。

**$\mathtt{pxi}\,\mathtt{pyi}$**。係数 $-\tfrac12\cdot2w_d\int(K\,\mathrm{Sh} - S\,\mathrm{Ch})\,ds = -w_dJ_1 = -w_dE/V$。$E = -w_+\mathtt{ca1}$、$w_dw_+ = b_p$ より $-w_dE/V = b_p\,\mathtt{ca1}/V = \mathtt{bzp}\,\mathtt{ca1}\,\mathtt{wss}$ ($\mathtt{wss} = 1/V$) で、`bzp*ca1*pxi*pyi*wss` と**一致**。

**$\mathtt{awu}\,\mathtt{pxi}$**。係数 $-\tfrac12\int[2K(-w_+S + w_-\mathrm{Sh}) - 2w_-(\mathrm{Ch} - K)S]\,ds = \int[w_+SK - w_-K\,\mathrm{Sh} + w_-S\,\mathrm{Ch} - w_-SK] = 2w_2I_{SK} - w_-\int(K\,\mathrm{Sh} - S\,\mathrm{Ch})$ ($w_+ - w_- = 2w_2$)。$I_{SK} = s_1^2/(2w_1)$ と $J_1 = E/V$ から $\dfrac{w_2s_1^2}{w_1} - \dfrac{w_-E}{V}$。第 1 項は `g1` $= (s_1w_2)^2/k_p = w_2s_1^2/w_1$ ($k_p = w_1w_2$)。第 2 項は $E = -w_+\mathtt{ca1}$ から $-w_-E/V = w_+w_-\mathtt{ca1}/V = b_p^2\mathtt{ca1}/V = \mathtt{ca1}\cdot\mathtt{cr1}^{(1)}$ ($\mathtt{cr1}^{(1)} = b_p^2/V$)。したがって `f_2^{(1)} = ca1*cr1 + g1` に**一致**。

**$\mathtt{dwu}\,\mathtt{pyi}$**。係数 $-\tfrac12\int[2w_-(\mathrm{Ch} - K)\mathrm{Sh} + 2\,\mathrm{Ch}(w_-S + w_+\mathrm{Sh})]\,ds = -\int[w_-\mathrm{Ch\,Sh} - w_-K\,\mathrm{Sh} + w_-S\,\mathrm{Ch} + w_+\mathrm{Sh\,Ch}] = -2w_1I_{\mathrm{Ch\,Sh}} + w_-J_1$ ($w_- + w_+ = 2w_1$)。$I_{\mathrm{Ch\,Sh}} = \mathrm{sh}_2^2/(2w_2)$ から $-\dfrac{w_1\mathrm{sh}_2^2}{w_2} + \dfrac{w_-E}{V}$。K64 側は第 2 の呼び出しで $\mathtt{cr1}^{(2)} = b_p^2\cdot(-1/V)$ なので $f_2^{(2)} = -\mathtt{ca1}\,b_p^2/V + \mathtt g_2$、$\mathtt g_2 = -(\mathrm{sh}_2w_1)^2/k_p = -w_1\mathrm{sh}_2^2/w_2$。$w_-E/V = -b_p^2\mathtt{ca1}/V$ なので両者は**一致**。

**$\mathtt{awu}\,\mathtt{pyi}$**。係数 $-\tfrac12\int[2w_d\mathrm{Sh}(-w_+S + w_-\mathrm{Sh}) + 2b_p(\mathrm{Ch} - K)\mathrm{Ch}]\,ds = -\int[-b_pS\,\mathrm{Sh} + w_dw_-\mathrm{Sh}^2 + b_p\mathrm{Ch}^2 - b_pK\,\mathrm{Ch}]$ ($w_dw_+ = b_p$)。$b_p$ でくくると ($w_dw_-/b_p = w_-/w_+ = w_d^2$)

$$
-b_p\Bigl[w_d^2I_{\mathrm{Sh}^2} + I_{\mathrm{Ch}^2} - \bigl(I_{S\mathrm{Sh}} + I_{K\mathrm{Ch}}\bigr)\Bigr]
$$

第 1・2 項は $I_{\mathrm{Sh}^2} = -\ell/2 + \mathrm{sh}_2\mathrm{ch}_2/(2w_2)$、$I_{\mathrm{Ch}^2} = \ell/2 + \mathrm{sh}_2\mathrm{ch}_2/(2w_2)$ と $1 - w_d^2 = 2w_2/w_+$、$1 + w_d^2 = 2w_1/w_+$ から

$$
w_d^2I_{\mathrm{Sh}^2} + I_{\mathrm{Ch}^2} = \frac{1 - w_d^2}{2}\ell + \frac{1 + w_d^2}{2w_2}\mathrm{sh}_2\mathrm{ch}_2 = \frac{w_2\ell}{w_+} + \frac{w_1\mathrm{sh}_2\mathrm{ch}_2}{w_2w_+} = \frac{w_2^2\ell + w_1\mathrm{sh}_2\mathrm{ch}_2}{w_2w_+} = \mathtt{aw1}
$$

(§5.3 の $\mathtt{pyi}^2$ の計算と同じ組み合わせ。) 第 3 項は表の $I_{S\mathrm{Sh}} = (-w_1c_1\mathrm{sh}_2 + w_2s_1\mathrm{ch}_2)/V$、$I_{K\mathrm{Ch}} = (w_1s_1\mathrm{ch}_2 + w_2c_1\mathrm{sh}_2)/V$ の和で、$I_{S\mathrm{Sh}} + I_{K\mathrm{Ch}} = [w_+s_1\mathrm{ch}_2 - w_-c_1\mathrm{sh}_2]/V$ ($w_1 + w_2 = w_+$、$w_1 - w_2 = w_-$)。したがって係数は

$$
-b_p\Bigl[\mathtt{aw1} + \frac{w_-c_1\mathrm{sh}_2 - w_+s_1\mathrm{ch}_2}{V}\Bigr]
$$

K64 側は $b_pf_1^{(1)} = b_p\bigl[(\mathtt{csw1} + \mathtt{cwsd})/V + \mathtt{csws} - \mathtt{aw1}\bigr]$ で、$V\,[(\mathtt{csw1} + \mathtt{cwsd})/V + \mathtt{csws}]$ を展開する ($\mathtt{csw1} = k_p\ell(\mathrm{ch}_2 - c_1)$、$\mathtt{cwsd} = c_1w_-(\phi_2 - \mathrm{sh}_2) - \mathrm{ch}_2w_+(\phi_1 - s_1)$、$V\,\mathtt{csws} = \mathrm{ch}_2w_1\phi_1 + c_1w_2\phi_2$):

$$
\begin{aligned}
&k_p\ell\,(\mathrm{ch}_2 - c_1) + c_1w_-(\phi_2 - \mathrm{sh}_2) - \mathrm{ch}_2w_+(\phi_1 - s_1) + \mathrm{ch}_2w_1\phi_1 + c_1w_2\phi_2\\
&\quad=\ \mathrm{ch}_2\phi_1\,(w_2 - w_+ + w_1) + c_1\phi_2\,(-w_1 + w_- + w_2) - w_-c_1\mathrm{sh}_2 + w_+\mathrm{ch}_2s_1
= -w_-c_1\mathrm{sh}_2 + w_+s_1\mathrm{ch}_2
\end{aligned}
$$

($k_p\ell = w_1w_2\ell = w_2\phi_1 = w_1\phi_2$ なので $k_p\ell\,\mathrm{ch}_2 = w_2\phi_1\mathrm{ch}_2$、$k_p\ell\,c_1 = w_1\phi_2c_1$。$\phi$ を含む項は係数が 0 になって消える。) よって $b_pf_1^{(1)} = -b_p\bigl[\mathtt{aw1} + (w_-c_1\mathrm{sh}_2 - w_+s_1\mathrm{ch}_2)/V\bigr]$ で、上の係数と**一致**。

**$\mathtt{dwu}\,\mathtt{pxi}$**。係数 $-\tfrac12\int[2b_pK(\mathrm{Ch} - K) - 2w_dS(w_-S + w_+\mathrm{Sh})] = -b_p\int[K\,\mathrm{Ch} - K^2 - w_d^2S^2 - S\,\mathrm{Sh}]$ ($w_dw_- = w_d^2b_p$、$w_dw_+ = b_p$) $= -b_p\bigl[J_2 - (\tfrac{1 + w_d^2}2\ell + \tfrac{1 - w_d^2}{2w_1}s_1c_1)\bigr] = b_p\bigl[\mathtt{aw2} - J_2\bigr]$ ($\mathtt{aw2}$ は §5.3 の $\mathtt{pxi}^2$ の計算と同じ組み合わせ)、すなわち $b_p\bigl[\mathtt{aw2} - (w_-s_1\mathrm{ch}_2 + w_+c_1\mathrm{sh}_2)/V\bigr]$。K64 側は $-b_pf_1^{(2)}$、第 2 の呼び出し ($\mathtt{wst} = -w_-$、$\mathtt{w12t} = w_+$、$\mathtt{swss} = -1/V$) の $\mathtt{cwsd}^{(2)} = c_1w_+\mathtt{xsh2} + \mathrm{ch}_2w_-\mathtt{xs1}$ で

$$
-f_1^{(2)} = \frac{\mathtt{csw1} + \mathtt{cwsd}^{(2)}}{V} - \mathtt{csws} + \mathtt{aw2},\qquad
V\Bigl[\frac{\mathtt{csw1} + \mathtt{cwsd}^{(2)}}{V} - \mathtt{csws}\Bigr]
= \mathrm{ch}_2\phi_1(w_2 + w_- - w_1) + c_1\phi_2(-w_1 + w_+ - w_2) - w_+c_1\mathrm{sh}_2 - w_-\mathrm{ch}_2s_1
= -w_+c_1\mathrm{sh}_2 - w_-s_1\mathrm{ch}_2
$$

なので $-b_pf_1^{(2)} = b_p\bigl[\mathtt{aw2} - (w_+c_1\mathrm{sh}_2 + w_-s_1\mathrm{ch}_2)/V\bigr]$ で**一致**。

**$\mathtt{awu}^2$ と $\mathtt{dwu}^2$**。$\mathtt{awu}^2$ の係数は

$$
-\tfrac12\Sigma_a,\quad \Sigma_a = w_+^2I_{S^2} - 2w_+w_-I_{S\mathrm{Sh}} + w_-^2I_{\mathrm{Sh}^2} + w_+w_-\bigl(I_{\mathrm{Ch}^2} - 2I_{K\mathrm{Ch}} + I_{K^2}\bigr)
$$

($b_p^2 = w_+w_-$。) 積分の表を入れて、$\ell$、$s_1c_1$、$\mathrm{sh}_2\mathrm{ch}_2$、$1/V$ の項に分ける:

- $\ell$: $\tfrac{\ell}2\bigl[w_+^2 - w_-^2 + 2w_+w_-\bigr] = \ell\,(2w_1w_2 + w_1^2 - w_2^2)$ ($w_+^2 - w_-^2 = 4w_1w_2$、$w_+w_- = w_1^2 - w_2^2$)。
- $s_1c_1$: $\dfrac{s_1c_1}{2w_1}\bigl[-w_+^2 + w_+w_-\bigr] = -\dfrac{w_2w_+}{w_1}s_1c_1$。
- $\mathrm{sh}_2\mathrm{ch}_2$: $\dfrac{\mathrm{sh}_2\mathrm{ch}_2}{2w_2}\bigl[w_-^2 + w_+w_-\bigr] = \dfrac{w_1w_-}{w_2}\mathrm{sh}_2\mathrm{ch}_2$。
- $1/V$: $\dfrac{2w_+w_-}{V}\bigl[(w_1c_1\mathrm{sh}_2 - w_2s_1\mathrm{ch}_2) - (w_1s_1\mathrm{ch}_2 + w_2c_1\mathrm{sh}_2)\bigr] = \dfrac{2w_+w_-}{V}\bigl[w_-c_1\mathrm{sh}_2 - w_+s_1\mathrm{ch}_2\bigr]$。

K64 側は $\tfrac12f_3^{(1)}$。$f_3^{(1)}$ ($\mathtt{swss} = 1/V$、$\mathtt{wst} = w_+$、$\mathtt{w12t} = w_-$、$\mathtt{cr1} = w_+w_-/V$、$\mathtt{wsm} = 4w_1w_2/V$) を $\phi_{1,2}$ を含む項と含まない項に分けて展開する:

- $c_1\phi_2$: $\dfrac{2w_+w_-^2}{V} + w_+\Bigl(-1 + \dfrac{4w_1w_2}{V}\Bigr) = w_+\Bigl[\dfrac{2(w_1 - w_2)^2 + 4w_1w_2}{V} - 1\Bigr] = w_+$。これは $-\mathtt{cr2}\,w_+\mathtt{xs1}$ の $c_1\phi_1$ の項 $-\dfrac{w_2w_+}{w_1}c_1\phi_1 = -w_2w_+c_1\ell$ と打ち消す ($w_+c_1\phi_2 = w_2w_+c_1\ell$)。
- $\mathrm{ch}_2\phi_1$: $-\dfrac{2w_+^2w_-}{V} + w_-\Bigl(1 + \dfrac{4w_1w_2}{V}\Bigr) = \dfrac{w_-}{V}\bigl[-2w_+^2 + V + 4w_1w_2\bigr] = -w_-$。これは $\mathtt{cr3}\,w_-\mathtt{xsh2}$ の $\mathrm{ch}_2\phi_2$ の項 $\dfrac{w_1w_-}{w_2}\mathrm{ch}_2\phi_2 = w_1w_-\mathrm{ch}_2\ell$ と、$-w_-\mathrm{ch}_2\phi_1 = -w_1w_-\mathrm{ch}_2\ell$ で打ち消す。
- 三角関数を含まない $\phi$ の項: $-\phi_2w_+\bigl(-1 + \tfrac{4w_1w_2}{V}\bigr) - \phi_1w_-\bigl(1 + \tfrac{4w_1w_2}{V}\bigr)$ ($\mathtt{dcw1}\,\mathtt{ws1}$ と $\mathtt{dcw2}\,\mathtt{ws2}$ の $-1$ の分)。$\phi_{1,2} = w_{1,2}\ell$ を入れて $\dfrac{\ell}{V}\bigl[w_2w_+(V - 4w_1w_2) - w_1w_-(V + 4w_1w_2)\bigr]$。$w_2w_+ - w_1w_- = 2w_1w_2 - (w_1^2 - w_2^2)$、$w_2w_+ + w_1w_- = V$ なので $= \ell\bigl[2w_1w_2 - w_1^2 + w_2^2 - 4w_1w_2\bigr] = -\ell\,(2w_1w_2 + w_1^2 - w_2^2)$。
- $\phi$ を含まない項: $s_1$、$\mathrm{sh}_2$ だけの項は、$2\mathtt{cr1}\,\mathtt{cwsd}$ から $\dfrac{2w_+^2w_-}{V}\mathrm{ch}_2s_1 - \dfrac{2w_+w_-^2}{V}c_1\mathrm{sh}_2$、$\mathtt{cr3}$ の項から $-\dfrac{w_1w_-}{w_2}\mathrm{ch}_2\mathrm{sh}_2$、$\mathtt{cr2}$ の項から $+\dfrac{w_2w_+}{w_1}c_1s_1$。

以上をまとめると

$$
f_3^{(1)} = -(2w_1w_2 + w_1^2 - w_2^2)\ell + \frac{w_2w_+}{w_1}s_1c_1 - \frac{w_1w_-}{w_2}\mathrm{sh}_2\mathrm{ch}_2 - \frac{2w_+w_-}{V}\bigl(w_-c_1\mathrm{sh}_2 - w_+s_1\mathrm{ch}_2\bigr)
$$

で、各項が $-\Sigma_a$ の $\ell$、$s_1c_1$、$\mathrm{sh}_2\mathrm{ch}_2$、$1/V$ の項と**一致する** ($\tfrac12f_3^{(1)} = -\tfrac12\Sigma_a$)。

$\mathtt{dwu}^2$ も同じ手順で、係数は $-\tfrac12\Sigma_d$、$\Sigma_d = w_-^2I_{S^2} + 2w_-w_+I_{S\mathrm{Sh}} + w_+^2I_{\mathrm{Sh}^2} + w_+w_-(I_{\mathrm{Ch}^2} - 2I_{K\mathrm{Ch}} + I_{K^2})$。$\ell$: $\tfrac\ell2(w_-^2 - w_+^2 + 2w_+w_-) = \ell(-2w_1w_2 + w_1^2 - w_2^2)$。$s_1c_1$: $\dfrac{s_1c_1}{2w_1}(-w_-^2 + w_+w_-) = \dfrac{w_2w_-}{w_1}s_1c_1$。$\mathrm{sh}_2\mathrm{ch}_2$: $\dfrac{1}{2w_2}(w_+^2 + w_+w_-) = \dfrac{w_1w_+}{w_2}$。$1/V$: $\dfrac{2w_+w_-}{V}\bigl[(-w_1c_1\mathrm{sh}_2 + w_2s_1\mathrm{ch}_2) - (w_1s_1\mathrm{ch}_2 + w_2c_1\mathrm{sh}_2)\bigr] = -\dfrac{2w_+w_-}{V}\bigl[w_+c_1\mathrm{sh}_2 + w_-s_1\mathrm{ch}_2\bigr]$。K64 側 $f_3^{(2)}$ ($\mathtt{swss} = -1/V$、$\mathtt{wst} = -w_-$、$\mathtt{w12t} = w_+$、$\mathtt{cr1} = -w_+w_-/V$、$\mathtt{wsm} = -4w_1w_2/V$) を同じ要領で展開すると、$c_1\phi_2$ の係数 $-w_-$、$\mathrm{ch}_2\phi_1$ の係数 $-w_+$ がそれぞれ $\mathtt{cr2}$・$\mathtt{cr3}$ の項と打ち消し、$\ell$ の項 $\ell(-w_1^2 + 2w_1w_2 + w_2^2)$ ($= -\ell(w_1^2 - 2w_1w_2 - w_2^2)$)、$s_1c_1$ の項、$\mathrm{sh}_2\mathrm{ch}_2$ の項、$1/V$ の項が残り、

$$
f_3^{(2)} = (2w_1w_2 - w_1^2 + w_2^2)\ell - \frac{w_2w_-}{w_1}s_1c_1 - \frac{w_1w_+}{w_2}\mathrm{sh}_2\mathrm{ch}_2 + \frac{2w_+w_-}{V}\bigl(w_+c_1\mathrm{sh}_2 + w_-s_1\mathrm{ch}_2\bigr) = -\Sigma_d
$$

と**一致する**。

**数値の確認** ([sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (C))。$w_1 = 0.9$、$w_2 = 0.35$、$\ell = 0.4$ で、10 個の単項式の係数を (a) 数値積分 $-\tfrac12\int(x'^2 + y'^2)\,ds$ (30 桁の mpmath)、(b) 上で導いた閉じた形、(c) K64 の `tztafs` の式を写したもの、の 3 通りで計算すると、すべての差が $10^{-31}$ 以下: $\mathtt{pxi}^2$ $-0.195285$、$\mathtt{pyi}^2$ $-0.201889$、$\mathtt{awu}^2$ $-0.00913292$、$\mathtt{dwu}^2$ $-0.00931675$、$\mathtt{awu}\,\mathtt{pxi}$ $+0.0726969$、$\mathtt{dwu}\,\mathtt{pyi}$ $-0.0751671$、$\mathtt{awu}\,\mathtt{pyi}$ $-0.00369671$、$\mathtt{dwu}\,\mathtt{pxi}$ $+0.00376764$、$\mathtt{awu}\,\mathtt{dwu}$ $-2.5558\times10^{-5}$、$\mathtt{pxi}\,\mathtt{pyi}$ $+0.0294720$。さらに、$\Delta z$ の厳密な積分と K64 の `dz` の差の分子を $\sin^2 + \cos^2 = 1$ 等で簡約した sympy の恒等式 ([sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (c))も、全体と 10 単項式ごとに 0 である。

**係数の対応 (§10 の表との関係)**: $\mathtt{g1} = (s_1w_2)^2/k_p$ と $\mathtt{g2} = -(\mathrm{sh}_2w_1)^2/k_p$ は、$\mathtt{awu}\,\mathtt{pxi}$ の $2w_2I_{SK}$、$\mathtt{dwu}\,\mathtt{pyi}$ の $-2w_1I_{\mathrm{Ch\,Sh}}$ の部分。$\mathtt{ca1}\cdot\mathtt{cr1}$ はこの 2 つと $\mathtt{pxi}\,\mathtt{pyi}$ の $J_1$ の部分。$\mathtt{csw1}$・$\mathtt{csws}$・`cwsd` の組 ($f_1$) は $\mathtt{awu}\,\mathtt{pyi}$・$\mathtt{dwu}\,\mathtt{pxi}$ の $J_2$ の部分 ($\phi$ の項が打ち消した残り) と $\mathtt{aw1}$・$\mathtt{aw2}$。$f_3$ の中の `cr2`・`cr3`・`dcw1`・`dcw2` と `ws1`・`ws2` は、$\mathtt{awu}^2$・$\mathtt{dwu}^2$ の $\ell$、$s_1c_1$、$\mathrm{sh}_2\mathrm{ch}_2$ の項を、$\phi_{1,2}$ を含む形に書き換えて小さな角で桁落ちしないようにした K64 の書き方で、1 つの積分に 1 対 1 で対応するのではなく、上の展開で $\phi$ の項が消えて残る組み合わせである。

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

1. 分割数 `ndiv = 1 + int(abs(al*hypot(ak, bz)/eps))` (`ak` は積分強度 K1、`al` は全長 $L$。$\mathrm{hypot}(K1, b_z)\,L/\varepsilon$ の整数部 + 1)。`eps` $= 0.2\times\mathrm{EPS}$ (EPS = 0 なら 0.2)。スライスの長さ $\ell = $ `aln0` $= L/n$。**この数の根拠**: Strang 分割の局所誤差は、$H_2$ と $\Delta H$ の交換子による $O(\ell^3)$ で ([03_symplectic_integrators.md](03_symplectic_integrators.md) §4 の BCH の展開: $\log S_2(\ell) = \ell(A + B) + \ell^3\bigl(-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]\bigr) + O(\ell^5)$。$A$、$B$ は $\Delta H$ と $H_2$ の半ステップの生成子)、$H_2$ の固有の角 ($w_1\ell$、$w_2\ell$。$\ell\,\mathrm{hypot}(K, b_z)$ の程度) が 1 より十分小さいことが分割の条件になる、というところまでは読める。$\mathrm{hypot}(K1, b_z)\,L/n < \varepsilon$ の形はそれに合う。係数 0.2 は経験的な値で、ソースにコメントが無く、$\Delta H$ の大きさ (傾き $q^2$) を含む誤差の見積もりからは導けない (未導出)。
2. 入射でラーモア系に: `px += bzp*y/2`、`py -= bzp*x/2` ($q \to$ 運動学的運動量 / $p$)。
3. 各スライス: $\Delta H$ の流れ (回転 + ドリフト、最初と最後は半分、途中は隣り合う 2 つの半分を 1 つにまとめた長さ `aln0`) → $H_2$ の厳密解 (§4。長さ `aln0`) → …、最後に $\Delta H$ の半分。**$\Delta H$ の回転の運動量の更新 (`pyi = -a24*pxi + a22*py(i)`) は、更新後の `pxi` を使っている (K64 の誤り。gs で確認、v2 は修正 = B333。付録 B)**。
4. 出射で元に戻す: `px = … - bzp*y/2`、`py = … + bzp*x/2`。

つまり $\Delta H(\tfrac\ell2)\,H_2(\ell)\,\Delta H(\ell)\,H_2(\ell)\cdots H_2(\ell)\,\Delta H(\tfrac\ell2)$ の **Strang 分割 (2 次)** ([03_symplectic_integrators.md](03_symplectic_integrators.md))。$B_z = 0$ の四極 (04 章) も同じ構造だが、$H_2$ は $\cos$、$\cosh$ の厳密解 (`tzsetparam0`: `w1 = sqrt(akkp)`、`phi1`、`phi2` が同じ $\phi_1$)。

**全体の確認** ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt)): SOL ($B_z = 13.0$ T = $b_z$ 1.3 /m、BOUND = 1) → QUAD ($K = 0.9$ /m²、$L = 0.4$ m、$n = 3$) → SOL ($B_z = 0$) の追跡と、独立な基準 (入口の $\exp(:f:)$ ((6.10) の厳密な流れ。$:f:$ の約束は §9) → (9.2) の $H$ の正準方程式の DOP853 → 出口の $\exp(:f:)$) の差は、3 つの入射で $7 \times 10^{-12} \sim 1.3 \times 10^{-11}$ (3 mm・2 mrad、$\delta = +0.01, 0, -0.02$)。**この差は分割 (Strang、$n = 3$) の誤差で、SOL・QUAD・端のフリンジの符号や単位の規約を含めて、追跡が全 $H$ の流れと一致する**。

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
- 運動量: `px = cchi1*px - schi1*pz0`、`py = cchi2*py - schi2*(schi1*px + cchi1*pz0)` (`px` は更新前の値、`pz0` $= q_z$)。導出: 向き $\vec q = (q_x, q_y, q_z)$ に位置と同じ回転を順に施す。$R_1'$ で $q_{x}' = c_1q_x - s_1q_z$、$q_{z1} = s_1q_x + c_1q_z$ (`pz1 = schi1*px + cchi1*pz0`)。次に $R_2'$ (第 2 行) で $q_y' = c_2q_y - s_2q_{z1}$、すなわち `py` の式の括弧が $q_{z1}$ (回転で $z$ 成分に混ざった分)。$x$ の更新は新しい $q_{z1}$ に依らない。$q_{z2} = s_2q_y + c_2q_{z1}$ は後の戻し `pz2` ($= 1 + \mathrm{sqrt1}(-q_x'^2 - q_y'^2)$、同じ値) (sympy: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (E))
- 戻し: `a = ds2/pz2`、`x -= a*px`、`y -= a*py`、`z += a`

位置ベクトルを、支点を原点として $\vec v = (x - \mathrm{DX},\ y - \mathrm{DY},\ s_{0s} - \mathrm{DZ})$ と書く (入口の面は $z = s_{0s}$。要素を DZ だけ前へずらすと、粒子は相対的に DZ だけ後ろ)。第 1 の回転 $R_1' = \begin{pmatrix}c_1&0&-s_1\\0&1&0\\s_1&0&c_1\end{pmatrix}$ ($x_1 = c_1x - s_1z$、$z_1 = s_1x + c_1z$) を $\vec v$ に施すと

$$
x_1 = c_1x - s_1(s_{0s} - \mathrm{DZ}),\qquad \mathrm{ds1} \equiv z_1 - s_{0s} = s_1x + (c_1 - 1)s_{0s} - c_1\mathrm{DZ}
$$

(`dcchi1` $= c_1 - 1$。) 第 2 の回転 $R_2' = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}$ を $(x_1, y, z_1)$ に施すと $y' = c_2y - s_2z_1 = c_2y - s_2(\mathrm{ds1} + s_{0s})$、$z_2 = s_2y + c_2z_1$、$\mathrm{ds2} \equiv z_2 - s_{0s} = s_2y + c_2\,\mathrm{ds1} + (c_2 - 1)s_{0s}$。これはコードの `y`、`ds1`、`ds2` と一致する。向きも同じ $R_2'R_1'$ で回る。回した後の粒子は面 $z = s_{0s}$ から $\mathrm{ds2}$ だけ前にいるので、向きに沿って $a = \mathrm{ds2}/q_{z2}$ だけ戻す。以上は、「要素の中央を支点にして要素を $(\mathrm{DX}, \mathrm{DY}, \mathrm{DZ})$ ずらし、$R_2'R_1'$ で回し、粒子の直線を新しい入口の面へ戻す」剛体運動 (11 章 §3) になっている。

MULT が渡す $\chi' = -\chi$ を代入すると $R_2'(-\chi_2)R_1'(-\chi_1) = \begin{pmatrix}1&0&0\\0&c_2&s_2\\0&-s_2&c_2\end{pmatrix}\begin{pmatrix}c_1&0&s_1\\0&1&0\\-s_1&0&c_1\end{pmatrix} = R_2R_1$ (11 章 (11.6) の $R_2$・$R_1$)。積の成分は、$R_2'(-\chi_2)$ の各行を $R_1'(-\chi_1)$ の各列に掛けて (行ごとに $R_1'(-\chi_1)$ の行の線形結合):

$$
\text{第 1 行} = (1, 0, 0)\cdot R_1' = (c_1, 0, s_1),\quad
\text{第 2 行} = c_2(0, 1, 0) + s_2(-s_1, 0, c_1) = (-s_1s_2,\ c_2,\ c_1s_2),\quad
\text{第 3 行} = -s_2(0, 1, 0) + c_2(-s_1, 0, c_1) = (-s_1c_2,\ -s_2,\ c_1c_2)
$$

($R_1'(-\chi_1) = \begin{pmatrix}c_1&0&s_1\\0&1&0\\-s_1&0&c_1\end{pmatrix}$ の第 1・2・3 行 $= (c_1, 0, s_1)$、$(0, 1, 0)$、$(-s_1, 0, c_1)$。) これは 11 章 §4.1 で展開した $R_2R_1$ と同じ行列である (sympy: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (E))。つまり MULT の CHI1・CHI2 は、COORD の $R_3R_2R_1$ ($R_3$ は `chi3` の回転で、コードでは最後に掛ける) と同じ符号の回転になっている。**ヘルプの MULT の式の $R_1^{(-)}$ (11 章 付録 A) との符号の関係は確認していない** (ヘルプの誤記か、`chi1` の符号の取り扱いを含めた読み方の違いか)。

**出口**: `s0s` $=$ (呼び出しの `alg`) $= -\mathtt{alg}_{\rm elem}$ (支点から見た出口の面は前方 $+\mathtt{alg}_{\rm elem}$ にあるので、位置ベクトルの $z$ 成分は $-s_{0s}$)。入口の $R_2'R_1'$ の**逆回転** $(R_2'R_1')^\top = R_1'^\top R_2'^\top$ を、位置 $\vec v = (x,\ y,\ -s_{0s})$ と向き $\vec q$ に施す ($R_2'^\top$ を先に、$R_1'^\top$ を後に)。入口と同様に成分で書くと、$R_2'^\top = \begin{pmatrix}1&0&0\\0&c_2&s_2\\0&-s_2&c_2\end{pmatrix}$、$R_1'^\top = \begin{pmatrix}c_1&0&s_1\\0&1&0\\-s_1&0&c_1\end{pmatrix}$ なので:

1. $R_2'^\top$ を $(x, y, -s_{0s})$ に: $y_1 = c_2y + s_2(-s_{0s}) = c_2y - s_2s_{0s}$ (`y1 = cchi2*y - schi2*s0s`)、$z_1 = -s_2y + c_2(-s_{0s})$。面からの偏差 $\mathrm{ds1} \equiv z_1 + s_{0s} = -s_2y - (c_2 - 1)s_{0s}$ (`ds1 = -schi2*y - dcchi2*s0s`。$z_1 = \mathrm{ds1} - s_{0s}$)。
2. $R_1'^\top$ を $(x, y_1, z_1)$ に: $x_1 = c_1x + s_1z_1 = c_1x + s_1(\mathrm{ds1} - s_{0s})$ (`x1 = cchi1*x + schi1*(ds1 - s0s)`)、$z_2 = -s_1x + c_1z_1$。面からの偏差 $\mathrm{ds2} \equiv z_2 + s_{0s} = -s_1x + c_1\mathrm{ds1} - (c_1 - 1)s_{0s}$ ($c_1(\mathrm{ds1} - s_{0s}) + s_{0s} = c_1\mathrm{ds1} - (c_1 - 1)s_{0s}$)。ここに DZ を足す: `ds2 = -schi1*x + cchi1*ds1 - dcchi1*s0s + dz`。
3. 向き: $q_{y}' = c_2q_y + s_2q_z$ (`pyi`)、$q_{z1} = -s_2q_y + c_2q_z$ (`pz1`)、$q_x' = c_1q_x + s_1q_{z1}$ (`pxi`)、$q_{z2} = -s_1q_x + c_1q_{z1}$ (`pz2`)。
4. 面へ戻す: 回した後の粒子は出口の面 $z = $ (面の位置) から ds2 だけ前にいるので、$a = \mathrm{ds2}/q_{z2}$ だけ向きに沿って戻し、$(+\mathrm{DX}, +\mathrm{DY})$ を足す: `x = x1 - a*pxi + dx`、`y = y1 - a*pyi + dy`、`z = z + a`。

(入口の式 `ds1 = schi1*x + dcchi1*s0s - cchi1*dz`、`ds2 = schi2*y + cchi2*ds1 + dcchi2*s0s` と比べると、回転が転置になったので $\sin$ の符号が逆で、$z$ 成分が $+s_{0s}$ でなく $-s_{0s}$ で始まる。sympy で、上の成分が $(R_2'R_1')^\top$ を $(x, y, -s_{0s})$ に施した結果と一致することを確認: [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) (E)。)

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

**確認 (K64 の誤り。gs で確認 = B337)** ([sol_04_rotation_branches.out.txt](outputs/sol_04_rotation_branches.out.txt) (c)): Fortran の式を numpy に写し ($b_z = 1.3$ /m、$\chi_1 = 0.03$、$\chi_2 = -0.02$、DZ = 1 cm、半長 0.5 m)、一様磁場中のらせん軌道の ODE (DOP853) が面 $z = s_{0s}$ に着く点と比べた。面の条件を解いた経路長 $l$ は ODE と $2 \times 10^{-17}$ で一致する。一方、位置・運動量の更新は、`xsincos` の戻り値 `dcosphi` $= \cos\theta - 1$ をそのまま使うと ODE と $3 \times 10^{-7} \sim 1.4 \times 10^{-6}$ ずれる ($x$: $-3.2 \times 10^{-7}$、$y$: $+1.1 \times 10^{-6}$、$p_x$: $+1.4 \times 10^{-6}$)。更新だけ `dcosphi` の符号を反転して $1 - \cos\theta$ にすると $10^{-9}$ に下がる (ただし面の条件は合わなくなる)。つまり、面の条件の式と更新の式で `dcosphi` の符号の想定が逆になっている。`tsolrot` のコメントアウトされた古い式 (`dcosphi = 2.d0*sin(.5d0*phi)**2` $= 1 - \cos\theta$) が `xsincos` ($\mathtt{dc} = \cos - 1$) に置き換えられたときの符号の取り違えの可能性がある。**K64 の誤り (確定。gs で確認 = B337)**: 独立な基準 (長さ 0 の傾いた要素は恒等写像でなければならない) で、入口 (`ent`) と出口の写像を numpy に写して合成すると、更新を Fortran のまま ($\cos\theta - 1$) にしたものは恒等から $2.2\times10^{-7}$・$3.9\times10^{-7}$・$1.4\times10^{-6}$・$6.5\times10^{-6}$ (px 最大 0.05 の 4 粒子) ずれ、更新の符号を反転したものは $10^{-17}$ で恒等になる。**gs で確認**: 3 GeV/c、SOL (BZ = 13 T、$b_z$ = 1.3 /m、BOUND = 1) → MULT ($L$ = 0、K2 = $10^{-9}$、CHI1 = 0.03、CHI2 = $-0.02$、DZ = 0.01) → SOL (BZ = 0) の `TrackParticles` は入力に戻らず、同じ $2.2\times10^{-7}$・$3.9\times10^{-7}$・$1.4\times10^{-6}$・$6.5\times10^{-6}$ ずれる (numpy の「Fortran のまま」と 3% 以内。傾きなしは $5\times10^{-16}$〜$6\times10^{-15}$ で戻る)。ずれの大きさは $\theta^2 p_t$ の程度で、$\theta = b_pl$ が小さい通常の設定では $10^{-6}$ 以下。K64 自身の誤りで、pysad の翻訳の誤りではない。**pysad には `tsolrot` の傾いたソレノイド ($B_z \ne 0$) の枝 (SOL 領域内の MULT・QUAD の CHI1・CHI2・DZ) は未実装** (NotImplementedError。k64・v2 とも)。実装するときは、k64 は K64 のまま、v2 は更新の符号を反転する。限界: 1 つの傾き (CHI1 = 0.03、CHI2 = $-0.02$、DZ = 0.01、領域の長さ 0) と $b_z$ = 1.3 /m。傾いた領域の中の MULT ($L \ne 0$) の本体は見ていない。試験 `test_v250_tsolrot_tilt.py`、記録 `docs/history/B337_tsolrot_tilt.md`。

### 8.4 GEO・CHI3 の決まり方 (未導出)

`tsgeo` (K64。pysad `sgeo.py`) が、GEO = 1 の SOL を基準に、領域の中の設計軌道を追い、各 SOL の CHI1・CHI2・CHI3・DX・DY・DZ を決める。この幾何は、軌道の追跡 (ドリフト・四極・SOL の写像) と座標の取り替えの繰り返しで、解析式の導出は行っていない (`docs/17_sol.md`)。**未導出**: 領域の出口の CHI3 を決める条件 (出口の後の最初の GEO = 1 の MARK で CHI3 = 0 になる) の式。

## 9. 端部 (DISFRIN)

SOL の境界で $b_z$ が $\mathrm{dbz}$ だけ変わるとき、`tsfrin` ($\exp(:f:)$、$f = -\dfrac{\mathrm{dbz}}{8p^2}p_\phi p_r$。[06_edge_fringe.md](06_edge_fringe.md) (6.10)) が、`ky_FRIN_SOL == 0` (DISFRIN = 0) のとき、`bz1 - bzs` の跳びごとに入る (`tsol.f:588`、`:753`)。ここで `bzs` は SOL 要素の**手前**の $b_z$ (規格化した $B_z$。追跡の状態として持つ量)、`bz1` はその SOL 要素の**後**の $b_z$ (`tsol.f:740`: `bz1 = merge(0, tfbzs(l, kbz), l == ke)`。`tfbzs` は次の SOL 要素までの領域の $b_z$ を返し、領域の最後の要素 `l == ke` では 0)。したがって `bz1 - bzs` $= \mathrm{dbz}$ は、この要素の位置での $b_z$ の跳びで、領域の入口の要素では `tsfrin(…, bzs)` と、跳び (0 から新しい $b_z$ へ) の分だけが引数になる (`tsol.f:588`)。**$\exp(:f:)$ の約束**: SAD の $:f:$ の約束で ([C_poisson_lie.md](C_poisson_lie.md) §8.2。$f$ を長さをかけたハミルトニアンとみて、その時間 1 の流れ。括弧の符号が付録の約束と逆)、$f$ はヘルプの DISFRIN の $H$ と同じ符号である。06 章 (6.10) の $f$ ($+\dfrac{b_0}{8p^2}p_\phi p_r$。付録の約束) とは $f_{\rm SAD} = -f_{(6.10)}$ の関係 (06 章 §9.3 で、Fortran の `bq = -0.125*dbz/(1 + g)` と流れの向きを照合している)。線形な効果 (一様な $B_z$ の中の $p_x \mathrel{+}= b_zy/2$ の跳び) は、正準運動量の連続で自動的に入る。

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
| `g1`、`g2` | $(s_1w_2)^2/k_p = w_2s_1^2/w_1$、$-(\mathrm{sh}_2w_1)^2/k_p = -w_1\mathrm{sh}_2^2/w_2$ | `tztafs` の $\mathtt{g1t}$。$f_2$ の一部で、$\mathtt{awu}\,\mathtt{pxi}$ の $2w_2I_{SK}$、$\mathtt{dwu}\,\mathtt{pyi}$ の $-2w_1I_{\mathrm{Ch\,Sh}}$ の部分 (§5.3) |
| `wr1`、`wr2` | $w_1/w_2$、$w_2/w_1$ | |
| `cxs1`、`cxs2` | `dc1*s1 - xs1` $= \sin\phi_1\cos\phi_1 - \phi_1$、`dch2*sh2 - xsh2` $= \sinh\phi_2\cosh\phi_2 - \phi_2$ | $\int\sin^2$、$\int\sinh^2$ の定数でない部分 (§5.3) |
| `aw1`、`aw2` | $\ell + \mathrm{cxs2}\,\mathrm{wr1}/w_+$、$\ell + \mathrm{cxs1}\,\mathrm{wr2}/w_+$ | `tztafs` の `zw1` (積分の定数の部分) |
| `cr2`、`cr3` | $c_1w_2/w_1$、$\mathrm{ch}_2w_1/w_2$ | |
| `csw1` | $k_p\ell(-(\cos\phi_1 - 1) + (\cosh\phi_2 - 1))$ | |
| `csws` | $(\mathrm{ch}_2w_1\phi_1 + c_1w_2\phi_2)/V$ | |
| `ca1` | $-\mathrm{ch}_2\,\mathrm{dc}_1 - \mathrm{dch}_2 - s_1\mathrm{sh}_2w_-/w_+$ | |
| `dxs` | $w_2\,\mathrm{xs}_1 - w_1\,\mathrm{xsh}_2$ | |
| `dcw1`、`dcw2` | $(\cos\phi_1 - 1)\phi_2$、$(\cosh\phi_2 - 1)\phi_1$ | |

(`tzparamp` の $p$ 微分 (`w1p` など) は、`tsolque` (CALC) で $\partial z/\partial\delta$ の行を作るためのもの。この章では扱わない。) 各係数は定義を Fortran から写したもの。(9.8) の積分への割り当ては §5.3 (10 個の単項式をすべて手で導き、sympy・数値積分で検算)。`g1`・`g2`・`ca1`・`csw1`・`csws`・`dxs` などの部分の対応は §5.3 の最後に書いた。`cr2`・`cr3`・`dcw1`・`dcw2`・`ws1`・`ws2` は、$\ell$・$s_1c_1$・$\mathrm{sh}_2\mathrm{ch}_2$ の項を $\phi_{1,2}$ を含む形に書き換える K64 の桁落ち回避の書き方で、1 つの積分には 1 対 1 に対応しない。

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
| `sol_05_gaps.py` | (A) 4 元連立の Cramer の解、(B) $K = 0$ の検算、(C) (9.8b) と 10 個の単項式の係数 (数値積分・手で導いた閉じた形・K64 の `tztafs`)、(D) 交差項 4 つの積分の微分、(E) §8.2 の行列の積・出口の成分・運動量の回転 | [sol_05_gaps.out.txt](outputs/sol_05_gaps.out.txt) |
| `sol_04_rotation_branches.py` | `trots`、`tsolrot` の枝と 3 次元の構成・らせん軌道の ODE の比較 | [sol_04_rotation_branches.out.txt](outputs/sol_04_rotation_branches.out.txt) |

端のフリンジの確認は [edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (06 章)。

## 14. 未確認・限界・使った条件

- **使った条件**: $p_0$ = 3 GeV の電子 ($B\rho = 10.007$ T m)。$b_z = 0.3 \sim 2.5$ /m、$K = 0.9 \sim 4$ /m²、$\ell = 0.4$ m、振幅 mm、角 mrad。$\delta = \pm0.02$。$H_2$ の確認は 4×4 の行列指数 (正準変数)、$z$・$\Delta H$ は DOP853。
- **SAD (gs) との比較**: `tsolqu` の $\Delta H$ の回転の更新順序 (B333) と `tsolrot` の傾いた枝の `dcosphi` の符号 (B337) だけ、gs の `TrackParticles` と比べた (付録 B)。ほかの式は、行列指数・ODE との比較で、gs とは比べていない。ヘルプの式は、K64 版同梱の原文 `Documents/SADHelp.html` の `solution-H2`・`solution-dH`・SOL の `DISFRIN` を直接読んで照合した (付録 A)。
- **確認できなかった**: (i) `tztafs` の係数の `cr2`・`cr3`・`dcw1`・`dcw2` を単独の積分に対応づけること (K64 の桁落ち回避の書き方で、1 対 1 には対応しない。展開した結果は 10 個の単項式すべてで確認した: §5.3)。(ii) `tsgeo` の幾何 (領域の出口の CHI3 を決める条件など) の式 (§8.4)。`tsolrot` の傾いたソレノイドの枝の `dcosphi` の符号 (§8.3。K64 の誤りと確定 = B337。傾いた領域の中の MULT の本体は未確認)。(iii) 放射あり (`tsolqur`)・スピン・CALC (`tsolque` の解析微分)。(iv) 四極自身のフリンジ (`ttfrin` の $B_z$ つきの枝) と SOL の端の合成。(v) $\Delta H$ の $K \ne 0$ との相互作用は分割の誤差 (Strang の 2 次) に任せている。分割数 $n$ の式の係数 ($\varepsilon = 0.2$) の根拠は、ソースに無く導けない (§7)。(vi) `tsolqu` の $\Delta H$ の回転の更新順序 (付録 B) は gs で確認済み (B333)。SOL 領域内の QUAD の $(x, \delta)$・$(y, \delta)$ のヤコビアンの欠陥 (下) の原因は未調査。
- **範囲**: 曲率 $1/\rho = 0$ の領域。$b_z$ は一様。軸上の場と、軸対称の端の場の 3 次までを扱う。

## 15. この章で使う Fortran の変数名

§2、§10 の表を参照 (すべての変数に定義と式を付けた)。補足: `bzs` $=$ SOL 要素の手前の規格化した $b_z$ (追跡が持つ状態)、`bz1` $=$ その要素の後の $b_z$ (`tsol.f:740`)、`dbz` $= \mathtt{bz1} - \mathtt{bzs}$ (§9)。`ky_FRIN_SOL` は SOL の DISFRIN キーワードの値 (`cmp%value(ky_FRIN_SOL)`)、`ky_F1_SOL` は F1。`bzph` $= \tfrac12\mathtt{bzs}/(1 + g)$ (放射の計算の内部で使う、範囲外)。`ndiv`・`eps`・`aln0` は §7。 `ap` $= q_x^2 + q_y^2$ (`px**2 + py**2`)、`dpz` $= \sqrt{1 - \mathtt{ap}} - 1$ (`sqrt1(-ap)` は桁落ちしない形)、`awu` $= a\,w_1/w_+$ ($= A\,w_1/w_+$、振動モードの振幅のうち位置から決まる分。(9.8a))、`dwu` $= dw\cdot w_2$ (指数モードの同じ量)、`a12`・`a14`・`a22`・`a24` は回転の `xsincos` の戻り値 (角 $r\,b_p$)。

## 付録 A. ヘルプの記述・K64 との違い (ソレノイド)

| 箇所 | 内容 |
|---|---|
| 原文との照合 (`SADHelp.html`) | `solution-H2` ($w_1$、$w_2$、$w_+$、$U$、$V$、$a, b, c, d$、$u_{1,2}$、$v_{1,2}$、$\Delta x$、$\Delta y$、$p_x$、$p_y$、$H_{2u}$、$z$) と `solution-dH` を原文で読み、$H_2$ の式は行列指数と一致を確認 (§4.3)。`solution-dH` は下の 2 行の食い違い。SOL の `DISFRIN` の $H$ は原文どおり |
| ヘルプ solution-dH の $\Delta x$、$\Delta y$ | $(p_{x0} + \tfrac{B_z}2y_0)\sin w\ell - (p_{y0} - \tfrac{B_z}2x_0)(\cos w\ell - 1)$ の式に **$1/B_z$ が無い** (**原文 `SADHelp.html` の `solution-dH` で確認**。TeX・平文とも)。(9.9) の導出と ODE ($3.6 \times 10^{-17}$)、K64 の `(a24*pxi - a14*py)/bzp` から、$1/B_z$ が要る (無いと $1.6 \times 10^{-6}$ ずれる)。ヘルプの誤記の疑い |
| ヘルプ solution-dH の $\Delta H$ | 式に `a_c` とあるが $a_x$ の誤記とみられる |
| ヘルプ solution-H2 | $K > 0$ にそろえる (回転) と書く。$H_2$ の式は行列指数と $10^{-18}$ で一致 |
| ヘルプの SOL の DISFRIN の $H$ | paraxial ($p^2$)。06 章 付録 A |
| ヘルプに無いもの | (1) 特性方程式の導出と $w_1$・$w_2$ の関係 (§4)。(2) 固有モード $y/x = \pm iw_d$、$\pm1/w_d$。(3) $z$ の更新が $-\tfrac12\int q^2$ であること (§5)。(4) ラーモア回転としての $\Delta H$ の解 (§6)。(5) `tzparams` の変数の定義 (§10) |

## 付録 B. K64 の誤り (確定 = 1・2、疑い (未決) = 3、取り下げは無い)

この章の範囲 ($H_2$ の厳密解、$z$ の更新、$\Delta H$ の解、Strang 分割、SOL + QUAD + SOL の全体、DISFRIN) で、K64 の式は行列指数・ODE・全 $H$ の流れと、分割の誤差の範囲で一致する (§4〜§7)。ただし次の 2 点は **K64 の誤りと確定した** (gs で確認)。3 つ目は疑い (未決)。ヘルプの誤記の疑いは付録 A。

1. **`tsolqu` の $\Delta H$ の運動量の回転がシンプレクティックでない (確定。gs で確認。docs/30 の B333。v2 は修正)** (ソースの読みと gs・pysad のヤコビアンの数値)。`tsolqu.f:122-123` (`tsolqum` の `:302-303`、`tsolqur` の `:487-488` も同じ) は

   - `pxi = a22*pxi + a24*py(i)`
   - `pyi = -a24*pxi + a22*py(i)`

   と書く。2 行目の `pxi` は、1 行目で更新した後の値である。正しい回転は $p_y' = -a_{24}\,p_{x,\rm old} + a_{22}p_y$。実際の写像は、$(p_{x,\rm old}, p_y) \to (p_x', p_y')$ の行列

   $$
   \begin{pmatrix}p_x'\\p_y'\end{pmatrix} = \begin{pmatrix}a_{22} & a_{24}\\ -a_{22}a_{24} & a_{22} - a_{24}^2\end{pmatrix}\begin{pmatrix}p_{x,\rm old}\\p_y\end{pmatrix},\qquad
   \det = a_{22}(a_{22} - a_{24}^2) + a_{22}a_{24}^2 = a_{22}^2 = \cos^2 w\ell
   $$

   で、正しい回転 (行列式 $a_{22}^2 + a_{24}^2 = 1$) と違い、面積を $\cos^2 w\ell$ 倍する ([sol_03_derivations.out.txt](outputs/sol_03_derivations.out.txt) (d)。sympy)。位置の更新 (`x`、`y`) は更新前の `pxi` を使っていて正しい。ループの後の最後の半ステップ (`tsolqu.f:151-156`) の運動量の更新も更新前の値を使っていて正しい。pysad の k64 は `track6.py`・`track6n_ext.py` に「SAD のまま: 更新後の pxi を使う」と写している。v2 は更新前の `pxi` で回す (`track6._tsolqu`・`track6n_ext.tz_step1`。`solemit` の Emittance[] の写像も同じ)。

   **数値** (同 (e) の k64。SOL (BOUND = 1、$b_z$) + QUAD ($K_1 = 0.36$ (積分強度)、$L = 0.4$ m、すなわち $K = K_1/L = 0.9$ /m²) + SOL ($B_z = 0$)、3 GeV、$\delta = 0$、$x = 1$ mm、$y = -0.5$ mm、中心差分の刻み $10^{-6}$ の 6×6 ヤコビアン $J$ で $\|J^\top SJ - S\|_{\max}$):

   | $b_z$ [/m] | $p_x = 0.003$ | $p_x = 0.05$ | $p_x = 0.1$ | $p_x = 0.1$、$p_y = 0.05$ |
   |---|---|---|---|---|
   | 1.3 | $1.9 \times 10^{-12}$ | $1.2 \times 10^{-7}$ | $1.9 \times 10^{-6}$ | $3.4 \times 10^{-6}$ |
   | 5 | $6.5 \times 10^{-11}$ | $2.4 \times 10^{-6}$ | $4.0 \times 10^{-5}$ | $3.2 \times 10^{-5}$ |
   | 10 | $7.3 \times 10^{-10}$ | $6.7 \times 10^{-6}$ | $1.1 \times 10^{-4}$ | $1.0 \times 10^{-4}$ |

   傾き $p_x = 3$ mrad の試験 (§7) の領域では $10^{-12} \sim 10^{-9}$ で見えないが、傾き 0.05 で $b_z = 1.3$ /m でも $1.2 \times 10^{-7}$ (ヤコビアンの行列式は $1 - 1.3 \times 10^{-7}$)、強い $B_z$ や大きな傾きで大きくなる。$\det J - 1 = -3.2 \times 10^{-6}$ ($b_z = 1.3$、$p_x = 0.1$、$p_y = 0.05$) は、予想される $\prod\cos^2(w\ell_i) - 1$ (各スライスの回転角の 2 乗の和程度) と同じ桁である。欠陥の大きさは $O\bigl((b_pq^2\ell)^2\bigr)$ ($w\ell \sim b_pq^2\ell/2$、$q = $ 傾き)。

   **v2** (更新前の `pxi` を使う。同 (e) の v2): 上の 12 点すべてで $|J^\top SJ - S| \le 6.5 \times 10^{-11}$ (中心差分の丸め。$b_z = 1.3$、$p_x = 0.1$ で $1.2 \times 10^{-11}$、$\det J - 1 = +1.1 \times 10^{-11}$)。

   **gs で確認** (B333): 3 GeV/c、SOL (BZ = 13 T、$b_z$ = 1.3 /m、BOUND = 1) + QUAD ($L$ = 0.4 m、K1 = 0.36 /m²) + SOL の `TrackParticles` の 13 点 ($\pm10^{-6}$) の中心差分ヤコビアンの 4×4 ブロック $(x, p_x, y, p_y)$ の $\|J^\top SJ - S\|$ は、点 $(10^{-3}, 0.05, -5\times10^{-4}, 0.02)$ で $1.54\times10^{-7}$、$(10^{-3}, 0.1, -5\times10^{-4}, 0.05)$ で $3.39\times10^{-6}$、$\det J - 1$ は $-1.50\times10^{-7}$、$-3.23\times10^{-6}$ (回転の角 $\theta = w\ell$ の $\theta^2$ の和 $3.0\times10^{-6}$ の見積りと同じ大きさ)。pysad の k64 の追跡は gs と全桁一致。v2 は 4×4 の欠陥が $4\times10^{-10}$ 以下 (丸め)、$\det J - 1$ が $-6\times10^{-11}$・$-5\times10^{-10}$ に下がる。**限界**: 1 構成 (SOL + QUAD + SOL) ・2 点。$(y, \delta)$・$(x, \delta)$ のヤコビアンの欠陥 ($p_x$ = 0.1 で $7.7\times10^{-6}$・$1.7\times10^{-6}$。gs の点 $\delta = 0.01$) は v2 でも残る (別の原因。$z$ の更新の近似などと思われる。未調査。K64 の疑いは強いが厳密解との比較ができず、**未決**。この章の (e) の点 ($\delta = 0$) では v2 は上のとおり $10^{-11}$)。小さな傾き ($p_x \sim 3\times10^{-3}$) では $\theta^2 \sim 10^{-11}$ で見えない。試験 `test_v250_sol_quad_rotation.py`、記録 `docs/history/B333_sol_quad_rotation.md`。

2. **`tsolrot` の傾いたソレノイドの枝 ($B_z \ne 0$、CHI1・CHI2・DZ ≠ 0) の `dcosphi` の符号の想定が、面の条件の式と更新の式で逆 (確定。gs で確認。docs/30 の B337。pysad は未実装)** (§8.3)。差は $\theta^2p_t$ の程度 (px = 0.05・$b_z$ = 1.3 /m で最大 $6.5\times10^{-6}$)。

3. **SOL 領域内の QUAD の $(x, \delta)$・$(y, \delta)$ のヤコビアンの欠陥 (疑い。未決)**: 上の 1 を直した v2 でも、gs の点 ($\delta = 0.01$、$p_x$ = 0.1) で $7.7\times10^{-6}$・$1.7\times10^{-6}$ 残る。K64 の疑いは強いが、厳密解との比較ができず確定できない。
