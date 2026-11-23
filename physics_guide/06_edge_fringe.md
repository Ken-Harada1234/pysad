# 06. 端部 (フリンジ)

**この章で分かること**: 磁石や空洞の端で、本体の $H$ (階段状の場を仮定したもの) に入らない効果を、Maxwell 方程式から導く。端の磁場の縦成分と高次の項 → ベクトルポテンシャル → $H$ の展開 → 階段の端の極限で $\delta$ 関数の項 → 生成関数 $f$ (リー変換 $\exp(:f:)$) → K64 が使う閉じた写像 (母関数)。偏向・四極・多極・ソレノイドの端が同じ手順で出ること、CAVI の端の $H_f$ の構造、FRINGE・DISFRIN の表を示す。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md)、[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)、付録 [C_poisson_lie.md](C_poisson_lie.md) ($\exp(:f:)$)。四極・偏向・多極の端の個別の式と検証は、[04_quad.md](04_quad.md) §7、[05_bend.md](05_bend.md) §6、[07_thin_multipole.md](07_thin_multipole.md) §9 にある。この章はそれらを**重複させず**、共通の導出と、残っていた穴 (母関数の意味、$z$ の更新、ソレノイド、CAVI・MULT、表) を埋める。

**読む順**: §1 (結論) → §3〜§5 (導出) → §6 (母関数) → §8〜§10 (要素ごと・表)。Fortran の対応は §11〜§12、変数名は §15。

## 1. 結論

1. 端の効果は、場が $s$ とともに変わるために現れる (i) **縦成分** $b_s$ と (ii) 横成分の**高次の項** ($\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ が要求するもの) である。どちらも、場の強さの $s$ 微分 ($\kappa'$、$\kappa''$) に比例する。
2. ベクトルポテンシャルを選ぶと ($\nabla\times\vec A = \vec b$)、端の効果は**横成分 $\vec a_\perp \propto \kappa'$** に集約される ($A_s$ は本体の場を作るだけ)。$H$ を $a_\perp$ で展開した 1 次の項は $H_1 = -(\vec p_\perp\cdot\vec a_\perp)/p_z$。階段の端 ($\kappa' = \kappa_0\,\delta(s)$) で積分すると、**端のリー生成関数** $f = \int H_1\,ds$ が出る (式 (6.8))。$\delta'$、$\delta''$ の項は、自由な運動に沿った微分になる (式 (6.9))。
3. この 1 つの手順から、(a) 偏向磁石の $f = \dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$ (05 章 (5.34))、(b) 四極の $f_{\rm in}$ ((4.27))、(c) $2m$ 極の位置のずれ $\Delta w$ ((7.17))、(d) ソレノイドの $f = -\dfrac{b_z}{8p^2}\,p_\phi p_r$ (ヘルプの DISFRIN の $H$) が、**すべて sympy で一致**した ([edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt))。
4. K64 の四極・多極の非線形フリンジは、$f$ の厳密な流れ $\exp(:f:)$ ではなく、**点変換 $x \to x + V + \tfrac12(V\cdot\nabla)V$ を母関数で正準化したもの**。$f = p\cdot V(x)$ の流れの 2 次までと一致し、**差は $K^3$ に比例**する ([edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (b): $K$ を 2 倍にして差が約 8 倍)。K64 の $z$ の更新は $\partial G/\partial p$ と一致する。
5. ソレノイドの端 (`tsfrin`) は $f$ の**厳密な流れ** (差 $4 \times 10^{-19}$)。Maxwell の端の場のローレンツ力の積分とも $b$ の 1 次の項で合う ($x$ で $0.9997$。他の成分は 3〜20%、§14)。
6. 偏向磁石の硬い端の Maxwell フリンジ ($\Delta x$、$\Delta z$) の因子は $p^2$ 型 (追跡・ヘルプ)。ローレンツ力の積分で独立に確認した ($\delta = 0, 0.3, -0.2$ で比 $1.0000$)。**CALC の `tbedgemaxwell` の COD の更新は $p^1$ 型で、$\delta \ne 0$ のとき $1 + \delta$ 倍ずれる (疑い。付録 B-1 の再確認)**。確定には SAD (gs) での確認が要る。
7. 端の効果が使える範囲は、場の展開が収束する $y \ll d$ ($d$ は端の幅)、1 次の摂動が効く $y^2/(2\rho_b d) \ll 1$。階段の極限 $d \to 0$ は**場の強さの 1 次の近似**で、それ以上の次数は意味を持たない (§4.4)。

## 2. 記号と単位

| 量 | 数学の記号 | K64 の Fortran | pysad |
|---|---|---|---|
| 正準運動量 ($p_0$ で規格化) | $p_x, p_y$ | `px`, `py` (追跡では $q = p_x/p$。`cod` は正準) | `c[1]`, `c[3]` |
| 全運動量 | $p = 1 + \delta$ | `pr` (`tbfrie.f`)、`p`、`pinv` $= 1/p$ | `pr` |
| 縦運動量 | $p_z = \sqrt{p^2 - p_x^2 - p_y^2}$ | `pzi`、`pv`・`pvi` $= \sqrt{p^2 - p_x^2}$ | `pvi` |
| 場の曲率半径 ($p_0$ に対して) | $\rho_b$ | `rhob` | `rhob` |
| 縦ベクトルポテンシャルなど | $a_x, a_y, A_s$ | (コードに現れない。導出用) | |
| 端での強さの変化 | $\kappa_0$ ($\kappa(s)$ の跳び) | `ak` (積分強度)、`akk` $= \mathrm{ak}/\mathrm{al}/(\ldots)$ | `ak` |
| 端の生成関数 | $f$ | (母関数 $G$ に直す。§6) | |
| 端の幅の長さ | $F1$、$F2$ | `f1in`、`fb1` など | |
| ソレノイドの軸上の場 | $b(s) = B_z/B\rho$ | `bz`、`bzs`、`dbz` (端での跳び) | `bz` |

$\kappa(s)$ は $2m$ 極の強さ ($m = 1$ 偏向、$m = 2$ 四極、…)。$\kappa_0$ は積分強度ではなく**長さあたり**の強さ (K64 の `ak/al`)。座標は直線の局所座標で、偏向磁石の曲率は端の幾何 ([05_bend.md](05_bend.md) §5) に入れ、端の磁場の導出では無視する。

## 3. 端の磁場 (Maxwell 方程式)

本体の磁場は、$s$ に依存しない横場 ($\kappa\,\mathrm{Im}\,w^m/m!$ 型) で、階段の端では $\kappa(s) = \kappa_0\Theta(s)$ と書く。しかし $\kappa$ が $s$ で変わると、$\nabla\cdot\vec b = 0$ が横場だけでは成り立たない。縦成分と高次の項を補う。

**3.1 偏向磁石** (端は $x$ に依らない直線。$g(s) = B_y(s)/B_0$ が 0 から 1 へ立ち上がる): [05_bend.md](05_bend.md) (5.30)〜(5.31)。$1/\rho_b$ を除いて

$$
b_x = 0,\quad b_y = g - \frac{y^2}{2}g'' + \frac{y^4}{24}g'''' - \cdots,\quad b_s = y\,g' - \frac{y^3}{6}g''' + \cdots
\tag{6.1}
$$

**3.2 四極と $2m$ 極**: [04_quad.md](04_quad.md) (4.24)〜(4.25)、[07_thin_multipole.md](07_thin_multipole.md) (7.16)。複素座標 $w = x + iy$ で、

$$
\Phi = \kappa\,\frac{\mathrm{Im}\,w^m}{m!} - \kappa''\,\frac{|w|^2\,\mathrm{Im}\,w^m}{4\,(m+1)!},\qquad \vec b = \nabla\Phi,\quad b_s = \partial_s\Phi = \kappa'\frac{\mathrm{Im}\,w^m}{m!} + \cdots
\tag{6.2}
$$

$\kappa''$ の項は、$m$ 回対称性を保つ形 ($|w|^2\mathrm{Im}\,w^m$) を選んだもの。$\nabla^2_\perp(|w|^2\mathrm{Im}\,w^m) = 4(m+1)\mathrm{Im}\,w^m$ から、$\partial_s^2\Phi$ と相殺する。**この選び方は、磁極が $m$ 回対称な磁石の場合の仮定**で、端の形が違えば高次の項は変わる (調和多項式の足し算の任意性)。偏向磁石 ($m = 1$) は、$x$ に依らない直線の端 (6.1) を仮定していて、(6.2) の $m = 1$ ($|w|^2 y$ を含み $x$ 依存) とは別の場。

**3.3 ソレノイド**: 軸上の $b(s) = B_z/B\rho$ (軸対称、$r = \sqrt{x^2 + y^2}$)。

$$
b_z = b - \frac{r^2}{4}b'' + \cdots,\qquad b_r = -\frac{r}{2}b' + \frac{r^3}{16}b''' - \cdots
\tag{6.3}
$$

($\nabla\cdot\vec b = r^{-1}\partial_r(rb_r) + \partial_sb_z$ が $b''$ まで 0。) 端で $b' \ne 0$ のとき、半径方向の場 $b_r$ が現れる。

**3.4 CAVI**: 軸上の電場 $E_z(s, t) = E(t)\,g(s)$ のとき、軸の近くで $E_r = -\tfrac{r}{2}E\,g'$、$B_\theta = \dfrac{r}{2c^2}\dot E\,g$ (Panofsky–Wenzell 型の関係。この章では導出せずに使う)。端 ($g' \ne 0$) で横方向の電場 $E_r$ が現れ、これが §9.1 の端の $H_f$。

## 4. ベクトルポテンシャルとハミルトニアン

### 4.1 ゲージ

横成分 $\vec a_\perp$ と縦成分 $A_s$ を、$\nabla\times\vec A = \vec b$ を満たすように選ぶ。成分で書くと

$$
b_x = \partial_yA_s - \partial_sa_y,\quad b_y = \partial_sa_x - \partial_xA_s,\quad b_s = \partial_xa_y - \partial_ya_x
\tag{6.4}
$$

選び方 (sympy で $\nabla\times\vec A = \vec b$ を確認した: [edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt) (a)〜(c)、(f)):

| 場 | $A_s$ (本体) | $\vec a_\perp$ (端。$\kappa'$ に比例) |
|---|---|---|
| 偏向 (1/ρ を除く) | $-x\,g$ | $a_x = -\dfrac{y^2}{2}g' + \dfrac{y^4}{24}g''' - \cdots$、$a_y = 0$ |
| 四極 | $-K\dfrac{x^2 - y^2}{2}$ | $a_x = -K'\dfrac{x^3 + 3xy^2}{12}$、$a_y = K'\dfrac{3x^2y + y^3}{12}$ |
| $2m$ 極 | $-\kappa\,\dfrac{\mathrm{Re}\,w^m}{m!}$ | $a_w \equiv a_x + ia_y = \dfrac{\kappa'}{4\,m!}\,w\left(\dfrac{w^m}{m+1} - \bar w^m\right)$ |
| ソレノイド | $0$ | $\vec a = (-y, x)\left(\dfrac{b}{2} - \dfrac{r^2}{16}b''\right)$ |

四極は $m = 2$ の特別な場合 ($a_w$ の式で一致)。$2m$ 極の $\nabla\times\vec A - \nabla\Phi$ は $\kappa'''$ の項だけが残る (これは $A$ の次の次数の補正で、この章では使わない)。ソレノイドの $b$ の 1 次の項 $\vec a = (-y, x)\,b/2$ は、本体のラーモア回転を作る (`09_solenoid.md`)。

### 4.2 $H$ の展開

$s$ を独立変数とする $H$ は (曲率と $E/(p_0v_0)$ の定数項を除く。01a 章)

$$
H = -\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2} - A_s
\tag{6.5}
$$

$\vec a_\perp$ で 1 次まで展開する。$P_\perp = p_\perp - a_\perp$ として $\partial H/\partial a_\perp = -P_\perp/p_z$ なので ($a_\perp = 0$ で)

$$
H = \underbrace{-p_z}_{\text{自由}} \;\underbrace{-\,A_s}_{\text{本体}} \;+\; H_1 + O(a_\perp^2),\qquad
H_1 = -\frac{p_xa_x + p_ya_y}{p_z},\quad p_z = \sqrt{p^2 - p_x^2 - p_y^2}
\tag{6.6}
$$

($a_\perp^2$ の項は $-(a_\perp^2)/(2p_z) - (p_\perp\cdot a_\perp)^2/(2p_z^3)$。$\kappa'$ の 2 乗で、階段の端では $\delta(s)^2$ になる。§4.4。)

「$y^4$ まで」とは、(6.1) の $b_y$ の $y^4 g''''/24$ ($a_x$ の $y^4g'''/24$) までを入れること。ヘルプの偏向の線形フリンジの生成関数 $V$ にある $y^4$ の項 ($y^4/(6\,\mathrm{F1}\,ho_b^2)$) は、$1/ho_b$ の 2 次の摂動の項で、05 章 §6.3 が扱う。この章の 1 次の $f$ に出る $y^4$ の項は、$g'''$ から来る $f_4$ (§4.3)。

### 4.3 階段の端の極限と生成関数

$f$ を、摂動 $H_1$ の自由な運動のまわりの 1 次 (Magnus 展開の最初の項) として

$$
f = \int H_1\bigl(x(s), p\bigr)\,ds,\qquad x(s) = x + \frac{p_x}{p_z}s,\ y(s) = y + \frac{p_y}{p_z}s,\ p_\perp\ \text{一定}
\tag{6.7}
$$

とする。流れ $\exp(:f:)$ は、$\Delta x = \partial f/\partial p_x$、$\Delta p_x = -\partial f/\partial x$、$\Delta z = \partial f/\partial p$ ([C_poisson_lie.md](C_poisson_lie.md) の約束)。階段の端 $\kappa' = \kappa_0\,\delta(s)$ では、$H_1$ が $\delta(s)$ に比例するので

$$
f = -\kappa_0\,\frac{p_x\tilde a_x + p_y\tilde a_y}{p_z}\quad\text{(入口。}\vec a_\perp = \kappa'\tilde a_\perp\text{)},\qquad f_{\rm out} = -f_{\rm in}
\tag{6.8}
$$

$\delta'$、$\delta''$ の項 (ソレノイドの $b'' = b_0\delta'$、偏向の $y^4$ の $g''' = \delta''$) は、部分積分で

$$
\int\delta^{(n)}(s)\,F(s)\,ds = (-1)^nF^{(n)}(0),\qquad F(s)\ \text{は自由な運動に沿った値}
\tag{6.9}
$$

**結果** ([edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt) (d)〜(f)):

- 偏向: $f = \dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$。$\Delta x = \dfrac{y^2}{2\rho_b}\dfrac{p^2 - p_y^2}{p_z^3}$、$\Delta p_y = -\dfrac{y}{\rho_b}\dfrac{p_x}{p_z}$、$\Delta z = -\dfrac{y^2}{2\rho_b}\dfrac{p\,p_x}{p_z^3}$。$p_y = 0$ で (5.35)。
- 四極: $f = \dfrac{K}{12p}(x^3p_x + 3xy^2p_x - 3x^2yp_y - y^3p_y)$ ((4.27)。差 0)。$2m$ 極: $\Delta w = \partial_{p_x}f + i\partial_{p_y}f = -\dfrac{\kappa_0}{p}\tilde a_w$ ((7.17)。$m = 2, 3, 4$ で差 0)。
- ソレノイド: $H_1 = -p_\phi\,(b/2 - r^2b''/16)/p_z$ ($p_\phi = xp_y - yp_x$)。$b/2$ の項は本体の回転。$b''$ の項に (6.9) ($n = 1$、$F = p_\phi r^2/(16p_z)$、$p_\phi$ は自由な運動で保存、$dr^2/ds = 2p_r/p_z$、$p_r = xp_x + yp_y$) を使い、

$$
f = -\frac{b_0}{8\,p_z^2}\,p_\phi\,p_r\ \ \xrightarrow{\ p_z\to p\ }\ -\frac{b_0}{8p^2}\,p_\phi p_r
\tag{6.10}
$$

ヘルプの DISFRIN の $H = -B_z/(8p^2B\rho)\,p_\phi p_r$ と一致する (K64 は paraxial の $p^2$)。

偏向の $y^4$ の項: $g''' = \delta''$ から (6.9) ($n = 2$) で $f_4 = -\dfrac{p_xp_y^2y^2}{2\rho_bp_z^3}$。$f_2 + f_4$ の $p_y^2$ までの展開は $\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_v}\left(1 - \dfrac{p_y^2}{2p_v^2}\right)$ ($p_v = \sqrt{p^2 - p_x^2}$)。**K64 の (5.36) は $p_y$ を $p_z$ に入れず ($p_v$)、この $p_y^2$ の補正を含まない**。補正の大きさは $(p_y/p_v)^2/2$ の相対値 (角の 2 乗)。**この補正は導出しただけで、ローレンツ力の積分では確かめていない** (展開が $y \ll d$ で成り立つ窓では $p_y$ を大きくできない。§14)。

### 4.4 適用範囲

(6.7) は場の強さの 1 次で、2 次の項 $O(a_\perp^2)$ は階段の端 ($d \to 0$) で $\delta(s)^2$ のように発散する。実際の端は有限の幅 $d$ を持ち、$\vec a_\perp \sim \kappa'/\ldots$ の大きさから、1 次で十分なのは**端の通過中の運動量のずれが小さい**とき: 偏向では $P_x$ のスパイクの高さ $\sim y^2/(2\rho_bd) \ll 1$ (数値で $y = 2$ mm、$\rho_b = 50$ m、$d = 1$ mm なら $4 \times 10^{-5}$)。また場の展開 (6.1)〜(6.3) は $y \ll d$ ($r \ll d$) でだけ収束する。$d \to 0$ の階段の極限で $y \gg d$ になる領域は、実際の磁石の端の形に依存する部分で、ヘルプ・K64 の式はその領域を表さない。

## 5. 要素ごとの生成関数 (結果の整理)

| 要素 | 端の項 | $f_{\rm in}$ (出口は符号反転) | 本章・他章 |
|---|---|---|---|
| 偏向 (K1 = 0) | $a_x = -y^2g'/(2\rho_b)$ | $\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$ | §4.3、05 章 (5.34)〜(5.36) |
| 四極 | $a_\perp = K'(\ldots)$ | (4.27) | 04 章 §7.2 |
| $2m$ 極 ($m \ge 3$) | $a_w = \kappa'(\ldots)$ | $\Delta w$: (7.17) | 07 章 §9 |
| ソレノイド | $a = (-y, x)(-r^2b''/16)$ | $-\dfrac{b_0}{8p^2}p_\phi p_r$ | (6.10)、09 章 |
| CAVI | $E_r = -\tfrac r2Eg'$ | $H_f$: §9.1 | |
| MULT ($n = 0$、双極子) | 偏向の式と同じ形 | `ttfrin` nord = 2: $\Delta x = \tfrac12 aki\,y^2$ | §9.2 |

## 6. 生成関数から閉じた写像へ (母関数の分解)

$f_{\rm in}$ (6.8) を四極で書くと $f = p_xV_x(x, y) + p_yV_y(x, y)$。**運動量について 1 次の生成子**は、位置の点変換 (ベクトル場 $V$ の流れ) と、運動量の逆転置による更新を作る。$\exp(:f:)$ の流れの位置の部分は

$$
\vec x_{\rm out} = \exp(V\cdot\nabla)\,\vec x = \vec x + V + \tfrac12(V\cdot\nabla)V + \tfrac16(V\cdot\nabla)^2V + \cdots
\tag{6.11}
$$

となり、閉じた式にならない (3 次以上が無限に続く)。K64 とヘルプの四極の母関数は、(6.11) の**2 次までで打ち切って正準にする**もの:

$$
G = H_0 + \tfrac12\left(\frac{\partial H_0}{\partial x_1}\frac{\partial H_0}{\partial p_{x2}} + \frac{\partial H_0}{\partial y_1}\frac{\partial H_0}{\partial p_{y2}}\right),\quad H_0 = p_{x2}V_x + p_{y2}V_y
\tag{6.12}
$$

型 2 の母関数 $F_2 = x_1p_{x2} + y_1p_{y2} + G$ から $\vec x_2 = \partial F_2/\partial\vec p_2 = \vec x_1 + V + \tfrac12(V\cdot\nabla)V$ が出る ($G$ の 2 項目の $\partial/\partial p_{x2}$ が $\tfrac12\sum_i V_i\partial_iV_x$)。**sympy で、四極の $V$ (4.27) について $x_2 - (x + V + \tfrac12 V\cdot\nabla V) = 0$ と、(4.29) の $x(1 + a/3 + b + t)$、$y(1 - a - b/3 + t)$ ($t = (a - b)^2/6$) との一致を確認した** ([edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (a))。つまり $t\,(x, y) = \tfrac12(V\cdot\nabla)V$ ($K$ の 2 次)。

- 点変換 (位置が運動量によらない) に直した写像は、運動量を $J^{-\top}$ で更新すれば、**構成からシンプレクティック**。
- $V \propto 1/p$ なので、$z_2 = z_1 + \partial G/\partial p$ は、$G$ の 1 次の部分 ($\propto 1/p$) と 2 次の部分 ($\propto 1/p^2$) で $-1/p$ 倍と $-2/p$ 倍。K64 の `dz1 = -(px(dx1 + t xi) + py(dy1 + t yi))/pr` は、`dx1` が $t$ を含むので 2 次の部分が 2 回入っていて、これと一致する ($\delta = 0, 0.01$ の 2 点で $\Delta z$ が $\partial G/\partial p$ と表示の 7 桁で一致し、「$-p_2\cdot(x_2 - x_1)/p$」とは $K$ の 2 次だけ違う: 同 (c))。$2m$ 極 ($m \ge 3$、`ttfrin` の一般の枝) は $V$ が 1 つの次数だけなので $z = -\mathrm{Re}(\Delta w\cdot\mathrm{cp1})$ (07 章 (7.18))。
- **厳密な $\exp(:f:)$ との差**: 数値積分 (DOP853、$f$ のハミルトン方程式、$\tau = 0 \to 1$) と K64 の式 (`ttfrin` 4 次) の差は、$K = 4 \times 10^3, 8 \times 10^3, 1.6 \times 10^4$ /m² で $x$ が $2.7 \times 10^{-9}, 2.3 \times 10^{-8}, 1.9 \times 10^{-7}$ (比 8.3、8.5)。$K^3$ に比例し、$K$ の 2 次までは一致する ([edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (b))。ただし**この 2 次の項 ($t$) 自体は、1 次の Maxwell の導出 (4.24) の範囲では物理的な意味を持たない** (母関数が 2 次を保つのは、点変換を厳密に正準にするため)。
- $m \ge 3$ の多極は (7.17) の 1 次の項だけ ($\tfrac12(V\cdot\nabla)V$ を足さない)。差 ($\kappa_0^2$ で減る) が出るのは四極だけ (07 章 §9.2)。

この「リー変換 → 点変換 + 逆転置」の形が、**母関数の分解**の中身: リー生成子 $f = p\cdot V$ の流れを、位置の点変換 (写像の 2 次までの級数) と運動量の逆転置の 2 段に分ける。BCH ([03_symplectic_integrators.md](03_symplectic_integrators.md)) の意味では、点変換の流れを 2 次で止めることは、生成子を $f + O(K^3)$ の別の生成子に置き換えることに当たる。

## 7. 線形フリンジ (F1・F2)

端が有限の幅で立ち上がるときの、硬い端からのずれの**線形な部分**。導出は、四極が [04_quad.md](04_quad.md) §7.4 ((4.31)〜(4.33): モーメント $I_1$、$I_2$ → $a = K_0I_1/p$、$b = K_0I_2/p$ → $\exp(a)$ 型の写像)、偏向が [05_bend.md](05_bend.md) §6.3 ((5.40)〜(5.42): $F1 = 6\int(g - g^2)ds$ から $\Delta x_{fr}$、$\Delta y_{fr}$、$\Delta y_{fra}$ と生成関数 $f_{\rm lin}$)。この章では繰り返さず、**統一的な見方**だけ書く。

- 本章の $f$ ((6.8)、1 次の摂動) は、端の幅 $d \to 0$ の極限 (硬い端)。$d$ が有限のときの補正は $d^1$、$d^2$ の項 (モーメント $I_n$) で、場の強さの 1 次では四極は線形な $(a, b)$ の写像、偏向では $\rho_b^{-2}$ (2 次の摂動) の $y$ の項になる。偏向の $F1$ は $y$ の集束 (2 次の摂動) の量で、四極の $F1$ ($x$ のずれ、1 次の摂動) と**定義が違う** (05 章 付録 A)。
- 多極 MULT の `tblfri` (§9.2) は、偏向の線形フリンジ ((5.42) の $\mathrm{dxfr}$、$\mathrm{dyfr}$、$\mathrm{dyfra}$) を、水平・垂直の両方の双極子成分 `ak0` (複素) に拡張した形。

## 8. 偏向磁石の端部 (`tbedge`)

### 8.1 構成

K1 = 0 の偏向磁石の端は、追跡の `tbendcore` (`tbend.f:372-527`) が 1 本の式にまとめ、CALC (`tbende`) は `tbedge` (`tbfrie.f:88-125`) の 3 段 ($\psi \ne 0$ のとき)

$$
\text{入口}:\ \texttt{tbedgedrift}\ (\text{面までの直進}) \to \texttt{tbedgemaxwell}\ (\text{Maxwell}) \to \texttt{tbedgebody}\ (\text{面から扇形の面へ})
$$

(出口は逆順) にする。`tbedgemaxwell` が (5.36) の更新 (COD) と行列、`tbedgedrift`・`tbedgebody` が面の幾何 (05 章 §5)。

### 8.2 ヘルプの式と `tbendcore` の同値

05 章 §6.2 に代数的な同値を示した ($ff = \tfrac14(f + fa)\,y\,\mathtt{pv1sqi}$ が、ヘルプの $y^2(1 - y^2/(12\rho_b^2))\,p^2/(2\rho_b(p^2 - p_x^2)^{3/2})$ に等しい)。ここでは数値で再確認した: ランダムな 2000 点 ($\delta \in \pm0.3$、$p_x/p \in \pm0.3$、$y \in \pm10$ mm) で、追跡の式 (`f`、`fa`、`ff` を $q$ 単位から正準に直したもの) とヘルプの式 (入口の $x_2$、$p_{y2}$、$z_2$) の最大差は $2 \times 10^{-20}$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (a))。なお、ヘルプの $p_{y2}$ の式は $p_1$ が 1 つ多く、$q$ の規格化と読めば一致する (05 章 付録 A)。

### 8.3 因子 $p$ の独立な再確認 (05 章 付録 B-1 の疑い)

**主張**: (5.35)〜(5.36) の $\Delta x = \dfrac{y^2}{2\rho_b}\dfrac{p^2}{p_v^3}$、$\Delta z = -\dfrac{y^2}{2\rho_b}\dfrac{p_xp}{p_v^3}$ は $p^2$、$p$ の因子を持つ。CALC の `tbedgemaxwell` (`tbfrie.f:209`、`:211`) の `cod(1) = cod(1) + bf*pr`、`cod(5) = cod(5) - bf*pxi` ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$) は、それぞれ $p^1$、$p^0$。

**(i) 導出** (sympy): (6.8) の偏向の $f$ から $\Delta x = \partial f/\partial p_x$ は $p_x/p_z$ の $p_x$ 微分 $= p^2/p_z^3$ (自由な運動の運動量 $p$ は $p_x$ を変えて $p_z$ を変える)。$\Delta z = \partial f/\partial p$ は $p_x/p_z$ の $p$ 微分 $-p\,p_x/p_z^3$。力が $\delta'$ なので、位置のずれの本質は「$P_x$ のスパイクの中で $dx/ds = P_x/p_z(P_x)$ が変わる」こと: $\Delta x = \epsilon\,\partial(P_x/p_z)/\partial P_x = \epsilon\,p^2/p_z^3$ ($\epsilon = y^2/(2\rho_b)$)。$\epsilon$ に $p$ は入らない ($b = eB/p_0$ が $p_0$ で規格化)。

**(ii) ローレンツ力の積分** (直線座標、Maxwell を $y^4$ まで満たす (6.1)、$g = (1 + \tanh(s/d))/2$、$d = 1, 0.5$ mm の差を $d \to 0$ に外挿、窓 $\pm8d$、$\rho_b = 50$ m、$y = 0.3$ mm、$p_y = 0$。$y \ll d$ の展開が使える範囲): $q = p_x/p = 0.1, 0.3$、$\delta = 0, +0.3, -0.2$ の 6 点で、**$\Delta x$、$\Delta p_y$、$\Delta z$ のどれも $p^2$ 型 (追跡・ヘルプ) との比が $0.9998 \sim 1.0001$**。CALC の更新との比は $\Delta x$、$\Delta z$ で **$1.3000$ ($\delta = 0.3$)、$0.8000$ ($\delta = -0.2$)** で、ちょうど $1 + \delta$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (b))。05 章の同じ確認 (bend_03 (b)) は $y = 2$ mm・$d = 1, 2$ mm だったが、そこでは $y > d$ で、(6.1) の展開の収束は保証されない窓だった。ここでは $y \ll d$ にしても同じ結論になる。

**(iii) CALC の内部の矛盾**: 解析行列 (`trans`、`pysad.lattice._tbedgemaxwell_mat`) は、$\Delta x$ は $p^2$ 型の、$\Delta z$ は $-bf\,p_x\,p$ ($p^1$ 型) の微分に、$\delta = 0, \pm0.1$ で $10^{-12}$ で一致する。COD の更新 ($p^1$・$p^0$) の中心差分とは、$\delta = 0.1$ で $1.7 \times 10^{-5}$、$\delta = -0.1$ で $2.5 \times 10^{-5}$ 違い ($\delta = 0$ でも $10^{-6}$。$\partial x/\partial\delta$ の成分 (1, 6) が 2 倍違う)。シンプレクティック性 $|M^	op SM - S|$ は、解析行列が 0 (丸め誤差以内)、COD の更新は $1.7 \times 10^{-5}$ ($\delta = 0.1$) で破れ、$p^2$・$p^1$ 型の更新は $1.5 \times 10^{-12}$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (c))。

**結論 (疑い)**: 追跡の式と解析行列は物理 ($p^2$ 型) に従い、**CALC の COD の更新だけが因子 $p$ を落としている**。ずれは $\delta\,y^2/(2\rho_b)$ 程度 ($y = 10$ mm、$\delta = 0.05$、$\rho_b = 10$ m で $2 \times 10^{-8}$ m)。**確定には SAD (gs) での COD の比較が要る** (この章では gs と比べていない)。pysad は K64 の直訳のまま (v2 も同じ)。

## 9. CAVI・SOL・MULT の端部

### 9.1 CAVI

ヘルプ: 入口の端の $\delta$ 関数の $H$ ($\omega$ は角周波数、$t$ は位相の時間、$V = $ VOLT + DVOLT)

$$
H_f = -\frac{eV}{L}\bigl(\sin(\omega t - d\varphi) + \sin d\varphi - \mathrm{offset}\bigr)\frac{x^2 + y^2}{4}\,\delta(s - s_0)\quad(\text{出口は符号反転})
\tag{6.13}
$$

$d\varphi$ と offset は、空洞の位相と放射による損失 (NORAD でだけ非 0) で決まる。構造は §3.4 の $E_r = -\tfrac r2E\,g'$ と同じ: 端で $E_z$ が $0 \to E\sin(\ldots)$ に跳ぶと、$E_r$ が $\delta(s)$ に比例する横方向の収束力 $F_r \propto -(r/2)E\,\delta(s)$ になり、$\int F_r\,dt$ が運動量のキック ($|\Delta p_r| \sim (eV/L)\,|\sin\ldots|\,r/2$) を与える。大きさは $H_f$ の $(x^2 + y^2)/4$ から $\Delta p_x = -\partial H_f/\partial x = (eV/L)(\sin\ldots)\,x/2$ と同じ ($1/v$ の因子と符号は未確認。下記)。**エネルギー**は $\partial H_f/\partial t$ から $\Delta E = -(\omega eV/L)\cos(\omega t - d\varphi)\,r^2/4$。

K64 `tcavfrin` (`tcav.f:209-249`): `vf = v/al/p0*.5`、`s0 = sin(dphis) - offset`、`dpt = vf*(sin(ph - dphis) + s0)/p1r`、`px += x*dpt`、`py += y*dpt` ($\Delta p_x = (V/2Lp_0)(\ldots)\,x/p$ の形)、`dh = -w*p0*vf*(x²+y²)/2*cos(ph - dphis)`、$h_2 = h_1 + dh$ から $\delta$ と $z$ ($z \to z\,p_2/p_1\cdot h_1/h_2$) と $dv$ を更新する。入口は `fringe` かつ `mfring >= 0` かつ `mfring /= 2` (`tcav.f:82`)、出口は `mfring /= 1` (`:185`)。

**未確認**: ここまでは、ヘルプの式とコードの**構造の対応を読んだだけ**で、(6.13) の導出 (Panofsky–Wenzell の関係と $E_z$ の端の形から)、$1/\beta$ の因子 ($\Delta p_x$ の物理的な式は $eEx/(2v)$ だが、コードは $v$ ではなく `1/p1r`。$q = p_x/p$ の単位と、$p_0$ の定義で説明がつくと推測)、`tcavfrin` の写像のシンプレクティック性は、数値で確認していない。詳しい扱いは `10_cavity.md` (予定)。

### 9.2 MULT

`tmulti` (`tmulti.f:125-139`、`:253-282`): 入口 `ttfrins(np,..., nord = 2n+2, akr(n), al, bzs)` を $n = 0, \ldots, nmmax$ (`dofr(n)` が真のもの) について、線形フリンジ ($\mathrm{F1}$・$\mathrm{F2}$: `f1in`、`f2in` の $\exp(a)$ 型 ((4.33) と同じ式)、`tblfri` ($\mathrm{FB1}$)) → 本体 → 出口 (逆順、$-\mathrm{akr}$)。`akr(n)` は複素 (実部が正規成分、虚部が skew)。`ttfrins` (`tquad.f:453-500`) は、成分を回転して正規成分にそろえ (角 $\theta = \tfrac{2}{\mathrm{nord}}\arg(\mathrm{ak}_1 + i\,\mathrm{ak}_2)$)、`ttfrin` を呼び、逆回転する。

- `ttfrin` の分岐: nord = 2 (双極子、$n = 0$): $x \mathrel{+}= \tfrac12\,\mathrm{aki}\,y^2$、$p_y \mathrel{-}= \mathrm{aki}\,p_xy$、$z \mathrel{-}= \tfrac12\,\mathrm{aki}\,p_xy^2$ ($\mathrm{aki} = \mathrm{ak}/(\mathrm{al}(1 + \delta))$)。これは (5.35) の $p_z \to p$、$p_x/p \to p_x$ の近似 ($\Delta x = y^2/(2\rho)$、$\Delta p_y = -y\,p_x/\rho$、$\Delta z = -y^2p_x/(2\rho)$) と同じ。双極子の端は (6.1) の直線の端を仮定していて、(7.17) を $m = 1$ に当てはめた式とは違う (§3.2)。nord = 4 は四極 (§6)、nord $\ge 6$ は (7.17)。
- `tblfri` (`tmulti.f:681-720`): 複素の `ck0` (水平の曲がりが実部、垂直が虚部) と `fb1` から、水平・垂直の両方の線形フリンジ: `dxfrx` $= \mathrm{fb1}^2/\rho_b/24$、`dyfrx` $= \mathrm{fb1}/\rho_b^2/6$、`dyfrax` $= 4\,\mathrm{dyfrx}/\mathrm{fb1}^2$ ((5.42) と同じ)。$x \mathrel{+}= \mathrm{dxfrx}\,\delta/p$、$p_y \mathrel{+}= (\mathrm{dyfrx} - \mathrm{dyfrax}\,y^2)\,y/p^2$ ($q$ 単位)。
- 以上は K64 のコードを読んだ範囲で、**MULT の端部を数値では確認していない**。ANGLE $\ne 0$ の MULT は `08_mult.md`。

### 9.3 ソレノイド

(6.10) の $f$ の流れが `tsfrin` (`tsol.f:885-976`)。Fortran のコメントの定義: $b = \tfrac12\,\mathrm{dbz}/p$、$bq = -b/4$、$bp = b/2$、$p_r = x p_x + y p_y$、$p_\phi = xp_y - yp_x$、$r_1 = r\,e^{bq\,p_\phi}$、$\phi_1 = \phi + bq\,p_r$、$z_1 = z + bp\,p_rp_\phi$ (コードの `bq = -0.125*dbz/(1 + g)`、`bp = -2*bq`)。これは、$f = -\tfrac{b_z}{8p^2}p_\phi p_r$ の座標 $(\phi, \ln r)$ (正準な共役は $(p_\phi, p_r)$) での流れ: $\Delta\ln r = \partial f/\partial p_r$、$\Delta\phi = \partial f/\partial p_\phi$、$\Delta z = \partial f/\partial p = +\tfrac{b_z}{4p^3}p_\phi p_r$ (積分)。

- **厳密な流れとの一致**: $f$ のハミルトン方程式を DOP853 で積分した流れと `tsfrin` は、$(x, p_x, y, p_y, z)$ で $4 \times 10^{-19}$ 以内 ($b = 1.7$ /m、3 mm・20 mrad) ([edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (a))。**$f$ は paraxial ($p_z \to p$) なので、閉じた式は $f$ の厳密な流れで、シンプレクティック性は構成から厳密**。
- **Maxwell の端の場との比較**: (6.3) の場 ($g = (1 + \tanh(s/d))/2$) のローレンツ力の積分 (kinematic momentum、窓 $\pm8d$)。モデルは「自由な運動 → 硬い端 (**正準運動量が連続**) → 一様な $B_z$ の中」で、これと ODE の差から**線形な端の集束 ($b$ の偶数次) を除くため、$b \to -b$ で符号が変わる部分 (奇数部分) を取り**、$d = 1, 0.5, 0.25$ mm から $d \to 0$ に外挿して `tsfrin` の効果と比べた ($b = 1$ /m、$x_0 = 0.3$ mm、$y_0 = -0.2$ mm、$p_x = 3$ mrad、$p_y = 2$ mrad。効果は $10^{-10}$ 程度): $x$: $0.9997$、$p_x$: $0.97$、$y$: $1.06$、$p_y$: $0.83$ ([edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (b))。効果が小さく ($10^{-11}\sim10^{-10}$)、$d$ への残りの依存 (2 つの外挿の差) が $p_y$ で 20% あるので、**$x$ 以外の成分は 3〜20% の確認**。
- 端部の線形な効果 (一様な $B_z$ の中の $p_x \mathrel{+}= b\,y/2$ の跳び) は正準運動量の連続で自動的に入る。**ヘルプの DISFRIN は、この非線形項 (6.10) だけを止める**。`tsol.f:588`、`:753`: `cmp%value(ky_FRIN_SOL) == 0` のとき `tsfrin` を $\mathrm{dbz}$ (端での $b_z$ の跳び) で呼ぶ。
- ヘルプの F1 (SOL) は EMITTANCE の計算にだけ効く (放射のため)。この章では扱わない。

## 10. FRINGE・DISFRIN の表

ヘルプの表 (QUAD・MULT・SEXT 以上・BEND で同じ形。SOL は DISFRIN が `tsfrin` だけ、CAVI は FRINGE と DISFRIN で端の $H_f$):

| FRINGE | DISFRIN = 0: 非線形 | 線形 (F1・F2) | DISFRIN ≠ 0: 非線形 | 線形 |
|---|---|---|---|---|
| 0 | 入口・出口 | なし | なし | なし |
| 1 | 入口 | 入口 | なし | 入口 |
| 2 | 出口 | 出口 | なし | 出口 |
| 3 | 入口・出口 | 入口・出口 | なし | 入口・出口 |

(ヘルプの BEND の表は、非線形の列だけを載せている。)

コードの判定 (K64 のソースを読んだ範囲):

| 要素 | 非線形フリンジ | 線形フリンジ |
|---|---|---|
| QUAD 追跡 `tquad` | 入口: `fringe .and. mfring > -4 .and. mfring /= 2`、出口: `... /= 1` (`tquad.f:49`、`:93`) | `mfring == 1 or 3` (入口)、`2 or 3` (出口) |
| QUAD CALC `tquade` | 入口: `fringe .and. mfring >= 0 .and. mfring /= 2` (`tquade.f:43`)、出口 `/= 1` (`:73`) | 同じ (`tqlfre`) |
| MULT 追跡 `tmulti` | 入口: `fringe`、`mfring /= 2`、$n$ ごとに `dofr(n)`、出口: `mfring /= 1` (`tmulti.f:125`、`:274`) | `mfring == 1 or 3` / `2 or 3`。F1・F2 (`f1in`..) と `tblfri` (FB1) |
| BEND 追跡 `tbendcore` (K1 = 0) | Maxwell は**常に入る**。DISFRIN は F1 の $y^3$ の項 (`dyfra`) だけを止める | F1 は `mfring > 0 or -1` (入口)、`> 0 or -2` (出口) |
| BEND K1 ≠ 0 追跡 | `ttfrin` (K1 の端) が DISFRIN に従う | F1 の項 |
| CAVI | `tcavfrin`: `fringe .and. mfring >= 0 .and. mfring /= 2` (入口)、`/= 1` (出口) | (なし。F1 は空洞にない) |
| SOL | `ky_FRIN_SOL == 0` のとき `tsfrin` ($b_z$ の跳びごと) | (なし) |

BEND の DISFRIN が Maxwell を止めない点は、ヘルプの記述と違う (05 章 付録 A、docs/30 の B134)。QUAD の FRINGE が $-1 \sim -3$ のときの追跡と CALC の差は 04 章 §7.1。

## 11. Fortran の対応表

| 内容 | ファイル:行 | サブルーチン | 式 |
|---|---|---|---|
| 四極・$2m$ 極の非線形フリンジ (追跡) | `tquad.f:325-452` | `ttfrin` | (6.12)、(4.29)〜(4.30)、(7.17)〜(7.18) |
| 四極の非線形フリンジ (CALC。$\delta$ の微分つき) | `tqfrie.f:1-129` | `tqfrie` | 同上 |
| MULT の回転つき非線形フリンジ | `tquad.f:453-500` | `ttfrins` | §9.2 |
| 線形フリンジ F1・F2 (追跡) | `tquad.f:501-528` | `tqlfri` (この版では未使用)。`tquad.f:53-67` 等にインライン | (4.33) |
| 偏向の端 (CALC) | `tbfrie.f:88-125` | `tbedge` | §8.1 |
| 偏向の Maxwell (更新と行列) | `tbfrie.f:190-229` | `tbedgemaxwell` | (5.36)。B-1 |
| 偏向の面の直進・面から扇形の面へ | `tbfrie.f:231-289`、`:127-188` | `tbedgedrift`、`tbedgebody` | 05 章 §5 |
| 偏向の線形フリンジ・Maxwell (追跡) | `tbend.f:372-527` | `tbendcore` | (5.36)、(5.42) |
| K1 ≠ 0 の端 (CALC のフリンジ。`ak` つき) | `tbfrie.f:1-86` | `tbfrie` | 05 章 §8 |
| 双極子成分の線形フリンジ (CALC) | `tbfrie.f:291-431` | `tbfrme`、`tbfrmle` | §9.2 |
| MULT の線形フリンジ (追跡) | `tmulti.f:681-720` | `tblfri` | §9.2 |
| CAVI の端 | `tcav.f:209-249` | `tcavfrin` | (6.13) |
| ソレノイドの端 (追跡・CALC) | `tsol.f:885-976`、`:978-` | `tsfrin`、`tsfrie` | (6.10) |

## 12. pysad の対応表

| 内容 | 関数 | 対応する Fortran | 注意 |
|---|---|---|---|
| 非線形フリンジ 4 次 / 一般 | `track6._ttfrin4`、`_ttfrin_n`、`_ttfrin6`、`_ttfrins` | `ttfrin`、`ttfrins` | 追跡は $q$ 単位。`track6n` (numba) に `ttfrin4_1` など |
| CALC の四極フリンジ (COD) | `lattice._tqfrie_cod`、`elkern._tqfrie` | `tqfrie` | 正準運動量 |
| 偏向の Maxwell (COD・行列) | `lattice._tbedgemaxwell_cod`、`_tbedgemaxwell_mat`、`elkern._tbedgemaxwell` | `tbedgemaxwell` | B-1 は直訳のまま |
| 偏向の端の合成 | `lattice._tbedge_cod`、`_tbedgedrift_cod`、`_tbedgebody_cod` | `tbedge`、`tbedgedrift`、`tbedgebody` | |
| 双極子の線形フリンジ (CALC) | `lattice._tbfrme_cod` | `tbfrme` | |
| MULT の線形フリンジ | `track6._tblfri` | `tblfri` | |
| CAVI の端 | `track6._tcavfrin`、`track6n.tcavfrin1` | `tcavfrin` | |
| ソレノイドの端 | `track6._tsfrin`、`sgeo._tsfrin`、`track6n_ext.tsfrin1` | `tsfrin` | |

## 13. 確認スクリプトと出力

`physics_guide/scripts/` (`py -3.12`、リポジトリのルートから。numpy・scipy・sympy。それぞれ数秒〜十数秒):

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `edge_01_gauge_generators.py` | $\nabla\times\vec A = \vec b$ (偏向・四極・$2m$ 極・ソレノイド)、(6.8)〜(6.10) から (4.27)・(7.17)・(5.35)・DISFRIN の $H$ が出る、偏向の $f_4$ ($p_y^2$ の補正) | [edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt) |
| `edge_02_generating_function.py` | 母関数 $G$ ↔ $x + V + \tfrac12 V\cdot\nabla V$ (sympy)、厳密な $\exp(:f:)$ との差 ($K^3$)、$z$ の更新 ($\partial G/\partial p$) | [edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) |
| `edge_03_bend_maxwell_recheck.py` | ヘルプ・追跡の同値、偏向の Maxwell の因子 $p$ (ローレンツ力の積分)、CALC の行列と更新 | [edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) |
| `edge_04_solenoid_fringe.py` | `tsfrin` = $f$ の流れ、Maxwell の端の場のローレンツ力の積分との比較 ($b$ の奇数部分) | [edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) |

## 14. 未確認・限界・使った条件

- **使った条件**: 偏向の ODE は $\rho_b = 50$ m、$y = 0.3$ mm、端の形 $\tanh$、$p_y = 0$、$q = 0.1, 0.3$、$\delta = 0, 0.3, -0.2$。エネルギーは規格化 ($p_0 = 1$) で、絶対値を使わない。四極の流れの比較は $K = 4 \times 10^3 \sim 1.6 \times 10^4$ /m² (大きな強さで 3 次の差を出した)。ソレノイドは $|b| = 1$ /m、振幅 0.3 mm。
- **SAD (gs) との比較はしていない** (ソース・ヘルプ・独立な数値計算だけ)。B-1 の確定には gs での COD の比較が要る。
- **確認できなかった**: (i) 偏向の $y^4/(12\rho_b^2)$ の係数 ((5.36) の $a$、$b$ の $y^2/\rho_b^2$ 項): 1 次の摂動 (6.8) では出ない ($1/\rho_b^3$ の 3 次の量)。導出していない。(ii) 偏向の $p_y$ 依存の補正 $f_4$ (§4.3): 導出のみで、ローレンツ力の積分では未確認 ($y \ll d$ の窓では $p_y$ を大きくして $y$ を動かせない)。K64 には無い。(iii) CAVI の端 (6.13) の導出と $1/\beta$ の因子、`tcavfrin` の数値確認。(iv) MULT の端部 (`ttfrins`、`tblfri`) の数値確認。(v) ソレノイドの $p_x$・$y$・$p_y$ は 3〜20% の確認 (§9.3)。(vi) 線形フリンジ ($F1$・$F2$) の導出は 04・05 章のまま。偏向の $L'$ の補正 (5.43) の導出は未。(vii) ヘルプは K64 版同梱の SADHelp.html の式 (テキスト形) を読んだ。画像として描画された式は参照していない。
- **範囲**: 場の強さの 1 次の摂動。$y \ll d$ の展開。階段の端 ($d \to 0$) の極限で、$d$ が有限のときの補正は線形フリンジ (§7) が担う。端の形が違えば高次の項 (特に (6.2) の調和多項式の足し算) は変わる。

## 15. この章で使う Fortran の変数名

| 変数 | 意味 (数学の標準記号ではない) | 式・なぜその形か |
|---|---|---|
| `pr` | $p = 1 + \delta$。`cod(6)` から | 規格化した全運動量 |
| `pvi`、`pvi2` | $p_v = \sqrt{p^2 - p_x^2}$、$p_v^2$ (`tbedgemaxwell`) | $p_y$ を入れない縦運動量 (`pr*sqrtl(pvi2/pr**2)`) |
| `af` | $a\,y/(\rho_bp_v)$、$a = 1 - y^2/(6\rho_b^2)$ | $\Delta p_y = -\mathtt{af}\,p_x$ の係数 |
| `bf` | $b\,y^2/(\rho_bp_v^3)$、$b = \tfrac12 - y^2/(24\rho_b^2)$ | $\Delta x = \mathtt{bf}\,p^2$ (正しい形)。CALC の更新は `bf*pr` |
| `yr` | $(y/\rho_b)^2$ | $a$、$b$、`c` の補正項 |
| `rhoe` | $\rho_e = \rho_bp$ ($\mathtt{tbendcore}$) | $q$ 単位の曲率半径 |
| `fa`、`f`、`ff` | `fa` $= y/(\rho_e\sqrt{1 - q^2})$、`f` $= a\cdot$`fa`、`ff` $= \tfrac14(f + fa)\,y/(1 - q^2)$ | $\Delta p_y$ と $\Delta x$ を同じ `fa` から作る (05 章 §6.2) |
| `pv1sqi` | $1/(1 - q_x^2)$ | `ff` の分母 |
| `dxfr1`、`dyfr1`、`dyfra1`、`dzfr1` | $F1^2/(24\rho_b)$、$F1/(6\rho_b^2)$、$4\,\mathrm{dyfr}/F1^2$、$\mathrm{dxfr}\sin\psi_1$ | 偏向の線形フリンジ (5.42)。`dyfra` は DISFRIN で 0 |
| `akk`、`aki` | `akk` $= \mathrm{ak}/\mathrm{al}/4$ ($m = 2$)。`aki` $= \mathrm{akk}/p$ | $K/(4p)$。$a = \mathrm{aki}\,x^2$、$b = \mathrm{aki}\,y^2$ |
| `a`、`b`、`d`、`ab`、`t` | $a = Kx^2/(4p)$、$b = Ky^2/(4p)$、$d = a + b$、$ab = a - b$、$t = ab^2/6$ | (4.28)〜(4.29)。$t = \tfrac12(V\cdot\nabla)V$ の係数 |
| `dx1`、`dy1` | $x(a/3 + b + t)$、$-y(a + b/3 - t)$ | 点変換のずれ ($t$ を含む)。`z` の更新は `dx1 + t*xi` ($t$ が 2 回) |
| `xx`、`yy`、`h`、`f` | $\partial x_2/\partial x_1$、$\partial y_2/\partial y_1$、$\partial x_2/\partial y_1$、$\det J$ | (4.30)。運動量の逆転置 |
| `kord`、`an` | 多極の次数 $m$、$\mathrm{an}(k) = k$ | (7.17) |
| `bq`、`bp`、`pphi`、`phi0`、`w` | ソレノイド: $bq = -\mathrm{dbz}/(8p)$、$bp = -2\,bq$、$p_\phi$、$bq\,p_r$、$e^{bq\,p_\phi}$ | (6.10) の流れ |
| `vf`、`s0`、`dpt`、`dh`、`offset`、`dphis` | CAVI: $V/(2Lp_0)$、$\sin\varphi_s - \mathrm{offset}$、キックの強さ、エネルギーの変化 | §9.1。(6.13) |
| `ck0` | MULT: 複素の $K_0$ (実部が水平、虚部が垂直) | `tblfri` |

## 付録 A. ヘルプの記述・K64 との違い (端部)

| 箇所 | 内容 |
|---|---|
| ヘルプ SOL の DISFRIN | $H = -B_z/(8p^2B\rho)\,p_\phi p_r$ は paraxial ($p_z \to p$)。(6.10) の厳密な形は $p_z^2$。K64 は paraxial |
| ヘルプ BEND の非線形フリンジ | $p_{y2}$ の $p_1$ が 1 つ多い (05 章 付録 A)。$p_y$ を $p_z$ に入れない (§4.3) |
| ヘルプ CAVI の $H_f$ | $\delta$ 関数の式のみ。導出・$1/\beta$ の因子は書かれていない |
| ヘルプの FRINGE・DISFRIN の表 | BEND の表は非線形の列だけ。BEND の DISFRIN は Maxwell を止めない (05 章 付録 A) |
| ヘルプに無いもの | (1) 端の磁場の縦成分と高次の項の導出 (§3)。(2) 母関数が 2 次までの打ち切りであること (§6)。(3) ソレノイドの端の $H$ の導出 (§4.3)。(4) CAVI の $H_f$ の導出 |

## 付録 B. K64 の誤りの疑い (別枠。確定には SAD (gs) での確認が要る)

**B-1 (再確認)**: CALC の `tbedgemaxwell` の COD の更新 (`tbfrie.f:209`、`:211`) が、$\Delta x$ で $p^1$、$\Delta z$ で $p^0$。正しい (追跡・ヘルプ・ローレンツ力の積分) のは $p^2$、$p^1$。ローレンツ力の積分で $\delta = 0.3$ で CALC は $\Delta x$、$\Delta z$ とも 1.3000 倍、$\delta = -0.2$ で 0.8000 倍 (§8.3)。解析行列は $p^2$ 型・$p^1$ 型の微分に一致し、COD の更新と内部で矛盾する。**疑い**。pysad は直していない。

**B-2 (05 章 B-2 と同じ)**: BEND の DISFRIN が Maxwell を止めない (ヘルプの記述との食い違い。追跡と CALC で一貫)。

**B-3 (規格化)**: ヘルプの偏向の式の $p_1$ の数が合わない (K64 のコードは正しい形。ヘルプの誤記とみられる)。

**この章で見つけた、K64 の誤りではないもの**: (i) K64 の四極の母関数 $G$ は $\exp(:f:)$ の 2 次までの打ち切り (仕様)。(ii) 偏向の Maxwell フリンジは $p_y$ の補正を持たない (4 次。仕様)。(iii) `tsetfringep` の $F1$ は $|F1|^2$ で符号を落とす (04 章の指摘)。
