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
6. 確認できなかったこと: `tztafs` の係数 (`aw1`、`ca1`、`cr2` など) の個別の導出、領域・座標の切替 (`tsolrot`) の数値確認、放射 (`tsolqur`)・スピンの扱い。§14。

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
| 三角・双曲線関数 | $\sin\phi_1$、$\sin\phi_1/\phi_1$、$\cos\phi_1$、$\cos\phi_1 - 1$ | `s1`、`xs1`、`c1`、`dc1` | 同名 | `xsincos` の戻り値。`dc1` は小さい角で桁落ちしない形 |
| | $\sinh\phi_2$、$\sinh\phi_2/\phi_2$、$\cosh\phi_2$、$\cosh\phi_2 - 1$ | `sh2`、`xsh2`、`ch2`、`dch2` | 同名 | |
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

$H_2$ は線形で厳密に解ける。$\Delta H$ は (§6) ラーモア回転になり厳密に解ける。$H_2$ を $P$ で展開すると $H_2 = -p + \dfrac{p_x^2 + p_y^2}{2p} - \dfrac{b_z}{2p}p_\phi + \dfrac{b_z^2}{8p}(x^2 + y^2) + \dfrac K2(x^2 - y^2)$ ($p_\phi = xp_y - yp_x$) で、$p_\phi$ の項が**ラーモア回転** (角速度 $b_z/2p$)、$b_z^2$ の項がソレノイドの収束。$K = 0$ ならラーモア系 (回転座標) で単なる等方振動子 ($w_1 = b_p/2$。下)。$K \ne 0$ では四極が回転系で時間に依存するので、この単純な回転では解けず、次節の正規モードを使う。

## 4. $H_2$ の厳密解

### 4.1 運動方程式

正準方程式 $x' = \partial H_2/\partial p_x = P_x/p$ など。$p_x' = -\partial H_2/\partial x$ は、$\partial P_y/\partial x = -b_z/2$ から $p_x' = \dfrac{b_z}{2}\dfrac{P_y}{p} - Kx$、同様に $p_y' = -\dfrac{b_z}{2}\dfrac{P_x}{p} + Ky$。運動学的運動量では $P_x' = p_x' + \tfrac{b_z}2y' = b_z\dfrac{P_y}p - Kx$、$P_y' = -b_z\dfrac{P_x}p + Ky$ (中間式の途中: $P_x' = p_x' + \tfrac{b_z}{2}y'$、$y' = P_y/p$)。$q = P/p$、$p$ は定数で

$$
x'' = b_py' - k_px,\qquad y'' = -b_px' + k_py
\tag{9.4}
$$

### 4.2 特性方程式と正規モード

$x = Xe^{\lambda s}$、$y = Ye^{\lambda s}$ を代入すると $(\lambda^2 + k_p)X = b_p\lambda Y$、$b_p\lambda X = -(\lambda^2 - k_p)Y$。行列式が 0 なので

$$
(\lambda^2 + k_p)(\lambda^2 - k_p) + b_p^2\lambda^2 = \lambda^4 + b_p^2\lambda^2 - k_p^2 = 0,\qquad
\lambda^2 = \frac{-b_p^2 \pm V}{2}
\tag{9.5}
$$

負の根 $\lambda^2 = -w_1^2$ ($w_1^2 = (b_p^2 + V)/2$) は振動、正の根 $\lambda^2 = +w_2^2$ ($w_2^2 = (V - b_p^2)/2$) は指数関数。関係 $w_1^2 + w_2^2 = V$、$w_1^2 - w_2^2 = b_p^2$、$w_1w_2 = k_p$ (sympy で確認、同 (a))。$b_p^2 = w_+w_-$ ($w_- = w_1 - w_2$) なので、$w_d \equiv b_p/w_+ = \sqrt{w_-/w_+}$。

**固有ベクトル**: $Y/X = (\lambda^2 + k_p)/(b_p\lambda)$。$k_p = w_1w_2$ を使うと、$\lambda = \pm iw_1$ で $Y/X = \pm i\,w_d$、$\lambda = \pm w_2$ で $Y/X = \pm1/w_d$ (sympy: 同 (e))。したがって

$$
\begin{aligned}
x &= \bigl[A\cos\phi_1 + B\sin\phi_1\bigr] + \bigl[C\cosh\phi_2 + D\sinh\phi_2\bigr],\\
y &= w_d\bigl[-A\sin\phi_1 + B\cos\phi_1\bigr] + \frac1{w_d}\bigl[C\sinh\phi_2 + D\cosh\phi_2\bigr]
\end{aligned}
\tag{9.6}
$$

($\phi_{1,2} = w_{1,2}s$)。係数は 4 つの初期条件 $x_0$、$y_0$、$q_{x0}$、$q_{y0}$ ($s = 0$) から決まる: $x_0 = A + C$、$y_0 = w_dB + D/w_d$、$q_{x0} = w_1B + w_2D$、$q_{y0} = -w_dw_1A + w_2C/w_d$。これを解いたものが、K64 の `a`、`bw`、`c`、`dw` と

$$
A = a,\quad B = bw,\quad C = b_p\,c,\quad D = b_p\,dw
\tag{9.7}
$$

の関係にある ($a = (w_2w_+x_0 - b_pq_{ym})/V$、$bw = (w_+q_{xm} - b_pw_2y_0)/V$、$c = (w_1w_dx_0 + q_{ym})/V$、$dw = (-w_dq_{xm} + w_1y_0)/V$。$q_{xm} = q_{x0} + b_py_0/2$、$q_{ym} = q_{y0} - b_px_0/2$)。4×4 を数値で解いた $(A, B, C, D)$ と K64 の式は $4 \times 10^{-19}$ で一致した (同 (e))。

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

$z$ は $p$ と共役。$H_2$ の $p$ 微分 (正準変数 $x, p_x$ などを固定) は $\partial H_2/\partial p = -1 - \dfrac{P_x^2 + P_y^2}{2p^2}$ ($P$ は $p$ によらない)。定数 $-1$ は基準軌道の経路 ($E/(p_0v_0)$ の項と合わせて $-dv$ の項) に入るので、残りが

$$
\Delta z_{H_2} = -\frac12\int_0^\ell\bigl(q_x^2 + q_y^2\bigr)\,ds,\qquad q = P/p
\tag{9.8}
$$

(04 章 §4.3 の導出と同じ形。) (9.6) の解を代入すると $\cos$、$\sin$、$\cosh$、$\sinh$ の積の積分になり、閉じた形になる。K64 の `tztafs` (`tsolque.f:39-86`) と、`tsolqu.f:137-141` の $z$ の更新 $\mathtt{bzp}(-\mathtt{awu}\,\mathtt{dwu}\,\mathtt{dxs}^2/\mathtt{akkp} + \mathtt{ca1}\,q_xq_y\,\mathtt{wss}) + \mathtt{dz1} + \mathtt{dz2}$ はこれを書いたもの。ヘルプは複素正規座標 $(u, p_u, v, p_v)$ で $H_{2u} = -p - iw_1up_u - w_2vp_v$ と書き、$z = z_0 + [-iup_u\partial w_1/\partial p - vp_v\partial w_2/\partial p + \Delta v]\ell$ ($up_u$、$vp_v$ は運動の定数)、$\partial w_1/\partial p = -w_1^3/(pV)$、$\partial w_2/\partial p = -w_2^3/(pV)$ (K64: `w1p = -w1**3/pr*wss`)。

**検証** (同 (c)): ODE の $-\int(q_x^2 + q_y^2)/2\,ds$ と K64 の $\Delta z_{H_2}$ (`dz1 + dz2 + …` の和) は、$(B_z, K, \delta) = (1.3, 0.9, 0.01)$ で $\Delta z = -2.0210463783 \times 10^{-7}$ と一致 (差 $8 \times 10^{-20}$)、ほかの 2 通りも差 $2 \times 10^{-19}$ 以下。

`tztafs` の係数 (`aw1`、`ca1`、`cr2`、…) は、(9.8) の積分を展開した各項のまとめで、**個別の導出はしていない** (全体は上の数値確認)。定義は §10 の表。

## 6. $\Delta H$ の解 (ラーモア回転)

$\Delta H = p - p_z - \dfrac{P^2}{2p}$ ($p_z = \sqrt{p^2 - P^2}$) は $P$ だけの関数 (K = 0)。運動方程式 ($\partial\Delta H/\partial P = P(1/p_z - 1/p) \equiv v$):

$$
x' = v_x,\quad y' = v_y,\quad p_x' = \frac{b_z}{2}v_y,\quad p_y' = -\frac{b_z}{2}v_x
$$

なので、運動学的運動量は $P_x' = p_x' + \tfrac{b_z}{2}y' = b_zv_y$、$P_y' = -b_zv_x$、すなわち $\vec P$ が**一定の大きさのまま角速度 $w = b_z(1/p_z - 1/p)$ で回転**する ($p_z$ も一定)。$x' = v_x = P_x\,w/b_z$ を積分して

$$
\Delta x = \frac1{b_z}\bigl[P_x\sin w\ell - P_y(\cos w\ell - 1)\bigr],\quad
\Delta y = \frac1{b_z}\bigl[P_x(\cos w\ell - 1) + P_y\sin w\ell\bigr],\quad w = \frac{b_z(p - p_z)}{pp_z}
\tag{9.9}
$$

$p_x = p_{x0} + \tfrac{b_z}2\Delta y$、$p_y = p_{y0} - \tfrac{b_z}2\Delta x$ ($P_{x,y}$ は初期値の運動学的運動量)。$z$ は $\partial\Delta H/\partial p = 1 - p/p_z + P^2/(2p^2)$ から $\Delta z = \bigl(\tfrac32 - \tfrac p{p_z} - \tfrac{p_z^2}{2p^2}\bigr)\ell$。

K64 の対応 (`tsolqu.f:107-123`): `r` $= -\mathtt{dpz}/(1 + \mathtt{dpz})\,\ell$ ($\mathtt{dpz} = p_z/p - 1$、`ap` $= q_x^2 + q_y^2$)、`call xsincos(r*bzp,...)` が回転 (`a22` $= \cos$、`a24` $= \sin$、`a14` $= \cos - 1$、`a12` $= \sin/\text{角}$)、`x += (a24*pxi - a14*py)/bzp` が (9.9)。回転角 $r\,b_p = w\ell$ と $z$ の更新 $-(3 + \mathtt{dpz})\,\mathtt{ap}\,r/(2(2 + \mathtt{dpz}))$ が、ヘルプの $w$ と $(3/2 - p/p_z - p_z^2/(2p^2))\ell$ に等しいことは sympy で確認した ([sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) (d)。差 0)。数値の流れ ($\Delta H$ の正準方程式、DOP853、$B_z = 2$ /m、$\ell = 0.5$ m、$p = 1.02$) とは、(9.9) で $3.6 \times 10^{-17}$ で一致する。

## 7. 本体の分割 (`tsolqu`)

K64 の `tsolqu` (`tsolqu.f:1-166`) は $B_z \ne 0$ のとき:

1. 分割数 `ndiv = 1 + int(abs(al*hypot(ak, bz)/eps))` (`ak` は積分強度 K1、`al` は全長 $L$。$\mathrm{hypot}(K1, b_z)\,L/\varepsilon$ の整数部 + 1)。`eps` $= 0.2\times\mathrm{EPS}$ (EPS = 0 なら 0.2)。スライスの長さ $\ell = $ `aln0` $= L/n$。
2. 入射でラーモア系に: `px += bzp*y/2`、`py -= bzp*x/2` ($q \to$ 運動学的運動量 / $p$)。
3. 各スライス: $\Delta H$ の流れ (回転 + ドリフト、最初と最後は半分、途中は隣り合う 2 つの半分を 1 つにまとめた長さ `aln0`) → $H_2$ の厳密解 (§4。長さ `aln0`) → …、最後に $\Delta H$ の半分。
4. 出射で元に戻す: `px = … - bzp*y/2`、`py = … + bzp*x/2`。

つまり $\Delta H(\tfrac\ell2)\,H_2(\ell)\,\Delta H(\ell)\,H_2(\ell)\cdots H_2(\ell)\,\Delta H(\tfrac\ell2)$ の **Strang 分割 (2 次)** ([03_symplectic_integrators.md](03_symplectic_integrators.md))。$B_z = 0$ の四極 (04 章) も同じ構造だが、$H_2$ は $\cos$、$\cosh$ の厳密解 (`tzsetparam0`: `w1 = sqrt(akkp)`、`phi1`、`phi2` が同じ $\phi_1$)。

**全体の確認** ([sol_02_chain.out.txt](outputs/sol_02_chain.out.txt)): SOL ($B_z = 13.0$ T = $b_z$ 1.3 /m、BOUND = 1) → QUAD ($K = 0.9$ /m²、$L = 0.4$ m、$n = 3$) → SOL ($B_z = 0$) の追跡と、独立な基準 (入口の $\exp(:f:)$ ((6.10) の厳密な流れ) → (9.2) の $H$ の正準方程式の DOP853 → 出口の $\exp(:f:)$) の差は、3 つの入射で $7 \times 10^{-12} \sim 1.3 \times 10^{-11}$ (3 mm・2 mrad、$\delta = +0.01, 0, -0.02$)。**この差は分割 (Strang、$n = 3$) の誤差で、SOL・QUAD・端のフリンジの符号や単位の規約を含めて、追跡が全 $H$ の流れと一致する**。

## 8. 領域と座標の切替 (BOUND・GEO・CHI3)

ヘルプの要約: SOL 要素は、ソレノイドの境界と、場が変わる所に置く。SOL の領域では、座標は**ソレノイドの軸**にとる (設計軌道が曲がっていても)。BOUND = 1 を領域の両端に置く。領域の中に DRIFT、BEND (直線)、QUAD、MULT を置いてよい (縦磁場が重なる)。領域の出口で座標が設計軌道に戻り、CHI3 は、出口の後の最初の GEO = 1 の MARK で 0 になるように決まる。DX・DY・DPX・DPY は GEO = 1 の SOL で設計軌道の軸に対するずれと傾きを与える。BZ だけを変える SOL (BOUND = 0) は座標を変えない。

実装: `tsolrot` (`tsol.f:9-252`、追跡) / `tsolrote` (`tsol.f:254-`、CALC) が、座標のずれ (`dx`、`dy`、`dz`) と回転 (`chi1`、`chi2`、`chi3`) を、**運動学的運動量に直してから**行い、正準運動量に戻す: $B_z \ne 0$ では `px += bz*y/(2 pr)`、`py -= bz*x/(2 pr)` (正準 → 運動学的) の後で座標を回し、最後に `px -= bzp*y/2`、`py += bzp*x/2` (`bzp = bzs/pr`、`bzs` は回転後の軸方向成分)。回転で縦磁場の成分 `bxs`、`bys`、`bzs` ($\phi_x = -\sin\chi_1$、$\phi_y = \cos\chi_1(-\sin\chi_2)$、$\phi_z = \cos\chi_1\cos\chi_2$ に $|b_z|$ をかけたもの) が作られ、傾いたソレノイドの場が横成分を持つ ($b_xs$、$b_ys$)。座標変換 (取り付け誤差) の一般論は `11_misalignment_and_coord.md` (予定)。

**この節は、ヘルプの記述とコードの読みで、`tsolrot` の式の導出と数値確認はしていない** (§14)。

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
| `s1`、`xs1`、`c1`、`dc1` | $\sin\phi_1$、$\sin\phi_1/\phi_1$、$\cos\phi_1$、$\cos\phi_1 - 1$ | |
| `sh2`、`xsh2`、`ch2`、`dch2` | 同じく双曲線関数 ($\phi_2$) | |
| `g1`、`g2` | $(s_1w_2)^2/k_p$、$-(\mathrm{sh}_2w_1)^2/k_p$ | `tztafs` の $\mathtt{g1t}$ (係数 $f_2$ の一部。$z$ の積分の項) |
| `wr1`、`wr2` | $w_1/w_2$、$w_2/w_1$ | |
| `cxs1`、`cxs2` | $(\cos\phi_1 - 1)\sin\phi_1 - \sin\phi_1/\phi_1$、$(\cosh\phi_2 - 1)\sinh\phi_2 - \sinh\phi_2/\phi_2$ (`dc1*s1 - xs1`、`dch2*sh2 - xsh2`) | 積分 $\int\cos\phi\sin\phi$ 型の項 |
| `aw1`、`aw2` | $\ell + \mathrm{cxs2}\,\mathrm{wr1}/w_+$、$\ell + \mathrm{cxs1}\,\mathrm{wr2}/w_+$ | `tztafs` の `zw1` (積分の定数の部分) |
| `cr2`、`cr3` | $c_1w_2/w_1$、$\mathrm{ch}_2w_1/w_2$ | |
| `csw1` | $k_p\ell(-(\cos\phi_1 - 1) + (\cosh\phi_2 - 1))$ | |
| `csws` | $(\mathrm{ch}_2w_1\phi_1 + c_1w_2\phi_2)/V$ | |
| `ca1` | $-\mathrm{ch}_2\,\mathrm{dc}_1 - \mathrm{dch}_2 - s_1\mathrm{sh}_2w_-/w_+$ | |
| `dxs` | $w_2\,\mathrm{xs}_1 - w_1\,\mathrm{xsh}_2$ | |
| `dcw1`、`dcw2` | $(\cos\phi_1 - 1)\phi_2$、$(\cosh\phi_2 - 1)\phi_1$ | |

(`tzparamp` の $p$ 微分 (`w1p` など) は、`tsolque` (CALC) で $\partial z/\partial\delta$ の行を作るためのもの。この章では扱わない。) 各係数は定義を Fortran から写したもので、(9.8) の積分の項への割り当ては導出していない。

## 11. Fortran の対応表

| 内容 | ファイル:行 | サブルーチン | 式 |
|---|---|---|---|
| 本体 (追跡) | `tsolqu.f:1-166` | `tsolqu` | §7、(9.6)、(9.9) |
| パラメータ ($B_z \ne 0$) | `tsolque.f:161-215` | `tzsetparams` | (9.1)、§10 |
| パラメータ ($B_z = 0$) | `tsolque.f:217-239` | `tzsetparam0` | 04 章 |
| $z$ の更新 | `tsolque.f:39-86` | `tztafs` | (9.8) |
| 放射あり・本体の変形 | `tsolqu.f:167-345`、`:347-539` | `tsolqum`、`tsolqur` | 範囲外 |
| CALC の本体 | `tsolque.f:313-` | `tsolque` | |
| 座標の切替 | `tsol.f:9-252`、`:254-` | `tsolrot`、`tsolrote` | §8 |
| ソレノイドの端 | `tsol.f:885-976`、`:978-` | `tsfrin`、`tsfrie` | 06 章 (6.10) |
| SOL 要素の領域 | `tsol.f:560-760` | | §8、§9 |

## 12. pysad の対応表

| 内容 | 関数 | 対応する Fortran | 注意 |
|---|---|---|---|
| パラメータと $z$ の更新 | `track6._TZ` (`tztafs` はメソッド) | `tzsetparams`、`tztafs` | 同名の属性 |
| 本体 (追跡) | `track6._tsolqu_bz_ens`、`_tsolqu_bz_final`、`track6n_ext.tsolqu_bz_all` | `tsolqu` | $B_z \ne 0$。放射の有無で並びが変わる (docs/17) |
| ソレノイドの端 | `track6._tsfrin`、`track6n_ext.tsfrin1`、`sgeo._tsfrin` | `tsfrin` | |
| 座標の切替・幾何 | `sgeo.py` | `tsolrot`、`tfgeo` ほか | docs/17_sol.md |

## 13. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `sol_01_h2_exact.py` | 特性方程式と $w_1$、$w_2$ (sympy)、ヘルプの式・K64 (`_TZ`) と行列指数の一致、$z$ の更新、$\Delta H$ の解、固有モードと `a`、`bw`、`c`、`dw` | [sol_01_h2_exact.out.txt](outputs/sol_01_h2_exact.out.txt) |
| `sol_02_chain.py` | SOL + QUAD + SOL の追跡と、独立な基準 (端の流れ + 全 $H$ の ODE) の比較、DISFRIN | [sol_02_chain.out.txt](outputs/sol_02_chain.out.txt) |

端のフリンジの確認は [edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (06 章)。

## 14. 未確認・限界・使った条件

- **使った条件**: $p_0$ = 3 GeV の電子 ($B\rho = 10.007$ T m)。$b_z = 0.3 \sim 2.5$ /m、$K = 0.9 \sim 4$ /m²、$\ell = 0.4$ m、振幅 mm、角 mrad。$\delta = \pm0.02$。$H_2$ の確認は 4×4 の行列指数 (正準変数)、$z$・$\Delta H$ は DOP853。
- **SAD (gs) との比較はしていない**。ヘルプの式は、K64 版同梱の原文 `Documents/SADHelp.html` の `solution-H2`・`solution-dH`・SOL の `DISFRIN` を直接読んで照合した (付録 A)。
- **確認できなかった**: (i) `tztafs` の係数の個別の導出。(ii) 領域・座標の切替 `tsolrot` の式の導出と数値確認 (§8。ヘルプとコードの読みだけ)。(iii) 放射あり (`tsolqur`)・スピン・CALC (`tsolque` の解析微分)。(iv) 四極自身のフリンジ (`ttfrin` の $B_z$ つきの枝) と SOL の端の合成。(v) $\Delta H$ の $K \ne 0$ との相互作用は分割の誤差 (Strang の 2 次) に任せていて、分割数 $n$ の式の根拠 ($\varepsilon = 0.2$) は確認していない。
- **範囲**: 曲率 $1/\rho = 0$ の領域。$b_z$ は一様。軸上の場と、軸対称の端の場の 3 次までを扱う。

## 15. この章で使う Fortran の変数名

§2、§10 の表を参照 (すべての変数に定義と式を付けた)。補足: `ap` $= q_x^2 + q_y^2$ (`px**2 + py**2`)、`dpz` $= \sqrt{1 - \mathtt{ap}} - 1$ (`sqrt1(-ap)` は桁落ちしない形)、`awu` $= a\,w_1/w_+$、`dwu` $= dw\cdot w_2$ (`tztafs` へ渡す量。モードの振幅)、`a12`・`a14`・`a22`・`a24` は回転の `xsincos` の戻り値 (角 $r\,b_p$)。

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

この章の範囲 ($H_2$、$z$ の更新、$\Delta H$、Strang 分割、SOL + QUAD + SOL の全体、DISFRIN) では、**K64 の誤りの疑いは見つからなかった**。K64 の式は、行列指数・ODE・全 $H$ の流れと、分割の誤差の範囲で一致する。ヘルプの誤記の疑い (付録 A) だけ。
