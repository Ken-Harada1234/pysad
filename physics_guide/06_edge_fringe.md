# 06. 端部 (フリンジ)

**この章で分かること**: 磁石や空洞の端で、本体の $H$ (階段状の場を仮定したもの) に入らない効果を、Maxwell 方程式から導く。端の磁場の縦成分と高次の項 → ベクトルポテンシャル → $H$ の展開 → 階段の端の極限で $\delta$ 関数の項 → 生成関数 $f$ (リー変換 $\exp(:f:)$) → K64 が使う閉じた写像 (母関数)。偏向・四極・多極・ソレノイドの端が同じ手順で出ること、CAVI の端の $H_f$ の構造、FRINGE・DISFRIN の表を示す。途中式は省かず、各ステップを別の行に書く。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md)、[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)、付録 [C_poisson_lie.md](C_poisson_lie.md) ($\exp(:f:)$ の約束)。四極・偏向・多極の端の個別の式と検証は、[04_quad.md](04_quad.md) §7、[05_bend.md](05_bend.md) §6、[07_thin_multipole.md](07_thin_multipole.md) §9 にある。この章はそれらを**重複させず**、共通の導出と、母関数の意味、$z$ の更新、ソレノイド、CAVI・MULT、表を扱う。

**読む順**: §1 (結論) → §2 (符号の約束) → §3〜§5 (導出) → §6 (母関数) → §8〜§10 (要素ごと・表)。Fortran の対応は §11〜§12、変数名は §15。

## 要約と目次

**この章で分かること**

- **何を導くか**: 磁石・空洞の端で、本体の $H$ (階段状の場を仮定) に入らない効果を、Maxwell 方程式から導く。端の磁場の縦成分と高次の項 → ベクトルポテンシャル → $H$ の展開 → 階段の端の極限 → 生成関数 $f$ (リー変換 $\exp(:f:)$) → K64 が使う閉じた写像 (母関数)。偏向・四極・多極・ソレノイドの端が同じ手順で出る。
- **主な結果**: 端の磁場 (6.1)〜(6.3)、1 次の摂動の生成関数 (6.7) と階段の端の $f$ (6.8)、ソレノイドの端 (6.10)、リー変換から点変換 + 運動量の逆転置への分解 (6.11)(6.12)、K64 の写像が別の生成子の厳密な流れであること (6.12b)、CAVI の端の $H_f$ (6.13)、MULT の複素強度の回転 (6.13b)。
- **対応**: Fortran は `ttfrin`・`ttfrins`・`tqfrie` (四極・多極)、`tbedge` 系・`tbendcore` (偏向)、`tcavfrin` (CAVI)、`tsfrin` (ソレノイド)。pysad は `track6._ttfrin4`・`lattice._tbedgemaxwell_cod`・`track6._tcavfrin`・`track6._tsfrin` など (§11・§12 の表)。FRINGE・DISFRIN の表は §10。
- **確かめ方**: sympy、DOP853 によるハミルトン方程式とローレンツ力の数値積分 (`edge_01`〜`edge_06`、§13)、K64 のソースの読み。gs と比べたのは偏向の Maxwell フリンジの因子 $p$ と CAVI の端の $1/\beta$ だけ (§14)。
- **K64・ヘルプの誤りの確定・疑い** (付録 B と一致): **確定** = B-1 (CALC の `tbedgemaxwell` の COD の更新が因子 $p$ を落とす。gs で確認。05 章 付録 B-1、docs/30 の B324。v2 は修正)、B-4 (CAVI の端の $1/\beta$ の欠落 = B330 と、$q$ の再スケールの欠落 = B332。v2 は修正)。**疑い (未決)** = B-2 (BEND の DISFRIN が Maxwell を止めない。仕様の疑い)、B-3 (ヘルプの偏向の式の $p_1$ の数。誤植か規格化の省略か確定できない)。取り下げは無い。

**導出の流れ**

1. [§1 結論](#1-結論)・[§2 記号・単位・符号の約束](#2-記号単位符号の約束): 結果の要点と、$\exp(:f:)$ の符号の約束 (6.0)(6.0a)。ヘルプの $H$ とは $H_{\rm help} = -f$。
2. [§3 端の磁場 (Maxwell 方程式)](#3-端の磁場-maxwell-方程式): $\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ を級数で解き、縦成分と高次の項を作る。
   - 3.1 偏向 (6.1a)(6.1b) → (6.1)、3.2 四極と $2m$ 極 (6.2)、3.3 ソレノイド (6.3)、3.4 CAVI の $E_r$ (この 4 つは見出しでなく太字の段落)。
3. [§4 ベクトルポテンシャルとハミルトニアン](#4-ベクトルポテンシャルとハミルトニアン): 端の効果を $\vec a_\perp \propto \kappa'$ に集約し、$f$ を出す。
   - [§4.1 ゲージ](#41-ゲージ): (6.4) と表。$\nabla\times\vec A = \vec b$ の手計算。
   - §4.2 $H$ の展開: (6.5) → (6.6)。
   - [§4.3 階段の端の極限と生成関数](#43-階段の端の極限と生成関数): (6.7) の導出 (6 段) → (6.8)(6.9) → 偏向・四極・$2m$ 極・ソレノイド (6.10)・偏向の $f_4$。
   - [§4.4 適用範囲](#44-適用範囲): 1 次の摂動が使える範囲。
4. [§5 要素ごとの生成関数 (結果の整理)](#5-要素ごとの生成関数-結果の整理): 結果の表。
5. [§6 生成関数から閉じた写像へ (母関数の分解)](#6-生成関数から閉じた写像へ-母関数の分解): (6.11) の級数を 2 次で打ち切って母関数 (6.12) を作り、K64 の写像が $f' = f + \tfrac16\,p\cdot\mathcal D^2V$ (6.12b) の厳密な流れであることを示す。
6. [§7 線形フリンジ (F1・F2)](#7-線形フリンジ-f1f2): 04 章・05 章の導出の統一的な見方だけ。
7. [§8 偏向磁石の端部 (`tbedge`)](#8-偏向磁石の端部-tbedge): 構成と、因子 $p$ の再確認。
   - [§8.1 構成](#81-構成)、[§8.2 ヘルプの式と `tbendcore` の同値](#82-ヘルプの式と-tbendcore-の同値)、§8.3 因子 $p$ の独立な再確認 (B-1。力の積分と CALC の内部の矛盾)、§8.4 Maxwell フリンジの $y^2/\rho_b^2$ の補正 (導出できていない)。
8. [§9 CAVI・SOL・MULT の端部](#9-cavisolmult-の端部): [§9.1 CAVI](#91-cavi) ((6.13) を 5 段で導出。$1/\beta$)、[§9.2 MULT](#92-mult) ((6.13a)(6.13b)。回転の向き、垂直の線形フリンジ)、[§9.3 ソレノイド](#93-ソレノイド) ((6.10) の流れ、`tsfrin` との一致、Maxwell の場との比較)。
9. [§10 FRINGE・DISFRIN の表](#10-fringedisfrin-の表)・[§11 Fortran の対応表](#11-fortran-の対応表)・[§12 pysad の対応表](#12-pysad-の対応表)・[§13 確認スクリプトと出力](#13-確認スクリプトと出力)。
10. [§14 未確認・限界・使った条件](#14-未確認限界使った条件)・[§15 この章で使う Fortran の変数名](#15-この章で使う-fortran-の変数名)・[付録 A](#付録-a-ヘルプの記述k64-との違い-端部)・[付録 B](#付録-b-k64-の誤り-確定--b-1b-4疑い--b-2b-3取り下げは無い)。

主な式の連鎖: (6.1a) → (6.1b) → (6.1)、(6.5) → (6.6) → (6.7a) → (6.7) → (6.8) → (6.9) → (6.10)、(6.11) → (6.12) → (6.12a) → (6.12b)、(6.13a) → (6.13b)。

**読み方**: 先に [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)・[C_poisson_lie.md](C_poisson_lie.md) ($\exp(:f:)$ の約束)・[01_variables_and_conventions.md](01_variables_and_conventions.md) §8 を読む。四極・偏向・多極の端の個別の式は 04 章 §7、05 章 §6、07 章 §9 にあり、この章は共通の導出だけを扱う。§9.1〜§9.3 は要素ごとに独立で、必要なものだけ読めばよい。式を追うときは、$p_x$ (正準運動量) と追跡の $q = p_x/p$ の違い、そしてヘルプの $H$ との符号の違い (§2) に注意する。§8.4・§14・付録は確認や照合のときに読めばよい。

**未確認・限界の所在**: [§14 未確認・限界・使った条件](#14-未確認限界使った条件) ((i)〜(viii)。偏向の $y^2/\rho_b^2$ 項、$f_4$、CAVI、MULT、ソレノイドの確認の度合い)。ヘルプとの違いは付録 A、K64 の誤りの確定・疑いは付録 B。

## 1. 結論

1. 端の効果は、場が $s$ とともに変わるために現れる (i) **縦成分** $b_s$ と (ii) 横成分の**高次の項** ($\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ が要求するもの) である。どちらも、場の強さの $s$ 微分 ($\kappa'$、$\kappa''$) に比例する。
2. ベクトルポテンシャルを選ぶと ($\nabla\times\vec A = \vec b$)、端の効果は**横成分 $\vec a_\perp \propto \kappa'$** に集約される ($A_s$ は本体の場を作るだけ)。$H$ を $a_\perp$ で展開した 1 次の項は $H_1 = -(\vec p_\perp\cdot\vec a_\perp)/p_z$。階段の端 ($\kappa' = \kappa_0\,\delta(s)$) で積分すると、**端のリー生成関数** $f = -\int H_1\,ds$ が出る (式 (6.8)。符号は [C_poisson_lie.md](C_poisson_lie.md) の約束、§2)。$\delta'$、$\delta''$ の項は、自由な運動に沿った微分になる (式 (6.9))。
3. この 1 つの手順から、(a) 偏向磁石の $f = -\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$ (05 章 (5.34))、(b) 四極の $f_{\rm in}$ ((4.27))、(c) $2m$ 極の位置のずれ $\Delta w$ ((7.17))、(d) ソレノイドの $f = +\dfrac{b_z}{8p^2}\,p_\phi p_r$ (ヘルプの DISFRIN の $H = -\dfrac{B_z}{8p^2B\rho}p_\phi p_r$ が $H = -f$ に当たる) が出る。確認は sympy ([edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt)、途中式は [edge_05_derivation_steps.out.txt](outputs/edge_05_derivation_steps.out.txt))。**注**: (c) の「$\Delta w$ が (7.17) と一致」は、(7.17) の $\tilde a_w$ を同じ式から取っているので、それ自体は同語反復に近い (微分の符号と係数の確認)。独立な確認は、(a) $\nabla\times\vec A = \vec b$ の確認 (§4.1。手計算と sympy)、偏向の端のローレンツ力の積分 (§8.3)、ソレノイドの端のローレンツ力の積分 (§9.3)、07 章のローレンツ力の積分 (thin_03)。
4. K64 の四極・多極の非線形フリンジは、$f$ の厳密な流れ $\exp(:f:)$ ではなく、**点変換 $x \to x + V + \tfrac12(V\cdot\nabla)V$ を母関数で正準化したもの**。$f = -p\cdot V(x)$ の流れの 2 次までと一致し、**差は $K^3$ に比例**する (位置の差の主要項は $\tfrac16(V\cdot\nabla)^2V$。[edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (b)。$K \le 2000$ /m² で全成分の差が $K$ を 2 倍にして 7.7〜8.1 倍、$K \ge 8000$ では $p_x$・$z$ の増え方が小さくなる ($K^4$ 以上の項が混じる))。K64 の $z$ の更新は $\partial G/\partial p$ と一致する。厳密には、K64 の写像は $f' = f + \tfrac16\,p\cdot(V\cdot\nabla)^2V + O(K^4)$ の別の生成子の厳密な流れ ((6.12b))。
5. ソレノイドの端 (`tsfrin`) は $f$ の**厳密な流れ** (差 $4 \times 10^{-19}$。$f$ が $(p_\phi, p_u)$ だけの関数なので流れが閉じた形になる、§9.3)。Maxwell の端の場のローレンツ力の積分とも $b$ の 1 次の項で合う ($x$ で $0.9997$。他の成分は、$d \to 0$ の外挿の取り方でずれが変わり、$p_x$ で 2〜3%、$y$ で 3〜6%、$p_y$ で 11〜17%、§9.3・§14)。
6. 偏向磁石の硬い端の Maxwell フリンジ ($\Delta x$、$\Delta z$) の因子は $p^2$ 型 (追跡・ヘルプ)。ローレンツ力の積分で独立に確認した ($\delta = 0, 0.3, -0.2$ で比 $1.0000$)。**CALC の `tbedgemaxwell` の COD の更新は $p^1$ 型で、$\delta \ne 0$ のとき $1 + \delta$ 倍ずれる (K64 の誤り。gs で確認した (05 章 付録 B-1、docs/30 の B324)。pysad の v2 は直してある)**。
7. 端の効果が使える範囲は、場の展開が収束する $y \ll d$ ($d$ は端の幅)、1 次の摂動が効く $y^2/(2\rho_b d) \ll 1$。階段の極限 $d \to 0$ は**場の強さの 1 次の近似**で、それ以上の次数は意味を持たない (§4.4)。

## 2. 記号・単位・符号の約束

**符号の約束** ([C_poisson_lie.md](C_poisson_lie.md))。ポアソン括弧は $[f, g] = \sum_i\left(\dfrac{\partial f}{\partial q_i}\dfrac{\partial g}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i}\right)$、$:f:\,g = [f, g]$。ハミルトニアン $H_{\rm ham}$ (正準変数は $(x, p_x, y, p_y, z, \delta)$。$s$ が独立変数) による長さ $s$ の転送は $\exp(-s:H_{\rm ham}:)$、すなわち

$$
\exp(:f:),\qquad f = -\int H_{\rm ham}\,ds
\tag{6.0}
$$

$\exp(:f:)$ の流れは、座標・運動量それぞれ

$$
\Delta x = -\frac{\partial f}{\partial p_x},\qquad \Delta p_x = +\frac{\partial f}{\partial x},\qquad \Delta z = -\frac{\partial f}{\partial p}\ \ (\text{同様に } y, p_y)
\tag{6.0a}
$$

(出所: $:f:\,x = [f, x] = -\partial f/\partial p_x$、$:f:\,p_x = [f, p_x] = +\partial f/\partial x$)。**この章・04 章・05 章・07 章のすべての $\exp(:f:)$ は、この約束**。SAD のヘルプの $\exp(:H:)$ は括弧の符号が逆の流儀で、ヘルプの $H$ は「ハミルトニアン的な量 (長さをかけたもの)」を表し、この約束の $f$ との関係は $H_{\rm help} = -f$ ([C_poisson_lie.md](C_poisson_lie.md) §8.2)。ヘルプの DISFRIN の $H$ (ソレノイド) や CAVI の $H_f$ は、この符号の違いに注意して読む。

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

本体の磁場は、$s$ に依存しない横場 ($\kappa\,\mathrm{Im}\,w^m/m!$ 型) で、階段の端では $\kappa(s) = \kappa_0\Theta(s)$ と書く。しかし $\kappa$ が $s$ で変わると、$\nabla\cdot\vec b = 0$ が横場だけでは成り立たない。縦成分と高次の項を補う。以下 (§3.1〜§3.3)、真空中の静磁場で $\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ を級数で解く (確認: [edge_05_derivation_steps.out.txt](outputs/edge_05_derivation_steps.out.txt) (a)〜(c))。

**3.1 偏向磁石** (端は $x$ に依らない直線。$g(s) = B_y(s)/B_0$ が 0 から 1 へ立ち上がる): [05_bend.md](05_bend.md) (5.30)〜(5.31)。$1/\rho_b$ を除いて $b_x = 0$、$b_y(y, s)$、$b_s(y, s)$ と書く。$x$ に依らないので、$\nabla\cdot\vec b = 0$ と $(\nabla\times\vec b)_x = 0$ は

$$
\partial_yb_y + \partial_sb_s = 0,\qquad \partial_yb_s - \partial_sb_y = 0
\tag{6.1a}
$$

($(\nabla\times\vec b)_y$、$(\nabla\times\vec b)_s$ は $x$ 微分と $b_x = 0$ で自明に 0)。中央面 ($y = 0$) で $b_y = g(s)$、上下対称から $b_y$ は $y$ の偶関数、$b_s$ は奇関数。$b_y = \sum_nc_n(s)y^n$、$b_s = \sum_nd_n(s)y^n$ とおくと、(6.1a) の $y^n$ の係数から

$$
(n+1)\,d_{n+1} = c_n',\qquad (n+1)\,c_{n+1} = -d_n'
\tag{6.1b}
$$

$c_0 = g$、$d_0 = 0$ から順に:

| $n$ | 使う式 | 結果 |
|---|---|---|
| 0 | $d_1 = c_0'$ | $d_1 = g'$ |
| 1 | $2c_2 = -d_1'$ | $c_2 = -g''/2$ |
| 2 | $3d_3 = c_2'$ | $d_3 = -g'''/6$ |
| 3 | $4c_4 = -d_3'$ | $c_4 = g''''/24$ |

奇数の $c$ と偶数の $d$ は 0 のまま。よって

$$
b_x = 0,\quad b_y = g - \frac{y^2}{2}g'' + \frac{y^4}{24}g'''' - \cdots,\quad b_s = y\,g' - \frac{y^3}{6}g''' + \cdots
\tag{6.1}
$$

(sympy で漸化式を $y^4$ まで解いて一致。edge_05 (a)。)

**3.2 四極と $2m$ 極**: [04_quad.md](04_quad.md) (4.24)〜(4.25)、[07_thin_multipole.md](07_thin_multipole.md) (7.16)。複素座標 $w = x + iy$ で、横場はスカラーポテンシャル $\vec b = \nabla\Phi$ から出る ($\nabla\times\vec b = 0$ は自動)。$\Phi$ に課す条件は $\nabla^2\Phi = (\partial_s^2 + \nabla_\perp^2)\Phi = 0$ ($\nabla\cdot\vec b = 0$)。最低次は $\Phi_0 = \kappa\,\mathrm{Im}\,w^m/m!$ ($m = 2$: $\kappa\,xy$、$b_x = \kappa y$、$b_y = \kappa x$)。これは $\nabla_\perp^2\Phi_0 = 0$ ($\mathrm{Im}\,w^m$ は調和関数) を満たすが、$\kappa$ が $s$ に依ると $\partial_s^2\Phi_0 = \kappa''\,\mathrm{Im}\,w^m/m! \ne 0$。これを打ち消す次の項を $\Phi_2 = c\,\kappa''\,|w|^2\mathrm{Im}\,w^m$ とおく (「$m$ 回対称」を保つ形)。

$\nabla_\perp^2 = 4\partial_w\partial_{\bar w}$、$|w|^2\mathrm{Im}\,w^m = (w^{m+1}\bar w - w\bar w^{m+1})/(2i)$ なので

$$
\nabla_\perp^2\bigl(|w|^2\mathrm{Im}\,w^m\bigr) = \frac{4}{2i}\Bigl[(m+1)w^m - (m+1)\bar w^m\Bigr] = 4(m+1)\,\mathrm{Im}\,w^m
$$

($\partial_w\partial_{\bar w}(w^{m+1}\bar w) = (m+1)w^m$、$\partial_w\partial_{\bar w}(w\bar w^{m+1}) = (m+1)\bar w^m$)。$\kappa''$ の項の和を 0 にする: $\kappa''\,\mathrm{Im}\,w^m/m! + c\,\kappa''\,4(m+1)\,\mathrm{Im}\,w^m = 0$ より $c = -1/(4(m+1)\,m!) = -1/(4(m+1)!)$。よって

$$
\Phi = \kappa\,\frac{\mathrm{Im}\,w^m}{m!} - \kappa''\,\frac{|w|^2\,\mathrm{Im}\,w^m}{4\,(m+1)!},\qquad \vec b = \nabla\Phi,\quad b_s = \partial_s\Phi = \kappa'\frac{\mathrm{Im}\,w^m}{m!} + \cdots
\tag{6.2}
$$

($m = 1, 2, 3, 4$ で $c$ と恒等式を sympy で確認。edge_05 (b)。) 残るのは $\kappa''''$ の項 (次の次数)。**この選び方は、磁極が $m$ 回対称な磁石の場合の仮定**で、端の形が違えば高次の項は変わる (調和多項式の足し算の任意性)。偏向磁石 ($m = 1$) は、$x$ に依らない直線の端 (6.1) を仮定していて、(6.2) の $m = 1$ ($|w|^2 y$ を含み $x$ 依存) とは別の場。

**3.3 ソレノイド**: 軸上の $b(s) = B_z/B\rho$ (軸対称、$r = \sqrt{x^2 + y^2}$)。$b_z = b + \beta_2(s)r^2 + \cdots$、$b_r = \alpha_1(s)r + \alpha_3(s)r^3 + \cdots$ とおく。$(\nabla\times\vec b)_\phi = \partial_rb_z - \partial_sb_r = 0$ の $r^1$ の係数から $2\beta_2 = \alpha_1'$。$\nabla\cdot\vec b = r^{-1}\partial_r(rb_r) + \partial_sb_z = 0$ の $r^0$、$r^2$ の係数から

$$
2\alpha_1 + b' = 0,\qquad 4\alpha_3 + \beta_2' = 0
$$

よって $\alpha_1 = -b'/2$、$\beta_2 = \alpha_1'/2 = -b''/4$、$\alpha_3 = -\beta_2'/4 = b'''/16$:

$$
b_z = b - \frac{r^2}{4}b'' + \cdots,\qquad b_r = -\frac{r}{2}b' + \frac{r^3}{16}b''' - \cdots
\tag{6.3}
$$

(確認: $\nabla\cdot\vec b$ は $b^{(5)}$ まで 0、$(\nabla\times\vec b)_\phi$ は $r^5b^{(6)}/384$ が残る。edge_05 (c)。) 端で $b' \ne 0$ のとき、半径方向の場 $b_r$ が現れる。

**3.4 CAVI**: 軸上の電場 $E_z(s, t) = E(t)\,g(s)$ のとき、軸の近くで $E_r = r E_r^{(1)}(s, t)$ とおく。ガウスの法則 $\nabla\cdot\vec E = r^{-1}\partial_r(rE_r) + \partial_sE_z = 0$ の $r^0$ の係数から $2E_r^{(1)} + E g' = 0$、すなわち $E_r = -\dfrac{r}{2}E\,g'$。磁場は、アンペール・マクスウェルの法則の $s$ 成分 $r^{-1}\partial_r(rB_\theta) = c^{-2}\partial_tE_z$ から $B_\theta = \dfrac{r}{2c^2}\dot E\,g$ ($\dot E = \partial_tE$)。端 ($g' \ne 0$) で横方向の電場 $E_r$ が現れ、これが §9.1 の端の $H_f$。

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

四極は $m = 2$ の特別な場合 ($a_w$ の式で一致)。ソレノイドの $b$ の 1 次の項 $\vec a = (-y, x)\,b/2$ は、本体のソレノイドの回転 (磁場に巻きつく運動) を作る ([09_solenoid.md](09_solenoid.md))。

**手計算での確認** (各行で $\nabla\times\vec A = \vec b$ を (6.4) に代入する):

- 偏向: $b_x = \partial_y(-xg) - \partial_s(0) = 0$。$b_y = \partial_s(-\tfrac{y^2}{2}g' + \tfrac{y^4}{24}g''') - \partial_x(-xg) = -\tfrac{y^2}{2}g'' + \tfrac{y^4}{24}g'''' + g$。$b_s = \partial_x(0) - \partial_y(-\tfrac{y^2}{2}g' + \tfrac{y^4}{24}g''') = yg' - \tfrac{y^3}{6}g'''$。(6.1) に一致。
- 四極: $b_x = \partial_y\bigl(-K\tfrac{x^2-y^2}{2}\bigr) - \partial_s a_y = Ky - K''\tfrac{3x^2y + y^3}{12}$。$b_y = \partial_sa_x - \partial_xA_s = -K''\tfrac{x^3 + 3xy^2}{12} + Kx$。$b_s = \partial_xa_y - \partial_ya_x = \tfrac{K'}{12}(6xy) + \tfrac{K'}{12}(6xy) = K'xy$。(4.25) に一致。
- $2m$ 極 (複素表記): $\partial_w \equiv (\partial_x - i\partial_y)/2$、$\partial_{\bar w} \equiv (\partial_x + i\partial_y)/2$ を使う。縦成分は $b_s = \mathrm{Im}(2\partial_wa_w)$ ($= \partial_xa_y - \partial_ya_x$)。$\partial_w\bigl[w^{m+1}/(m+1) - w\bar w^m\bigr] = w^m - \bar w^m$ なので $2\partial_wa_w = \kappa'(w^m - \bar w^m)/(2\,m!)$、$b_s = \kappa'\,\mathrm{Im}\,w^m/m!$ ($= \partial_s\Phi$ の $\kappa'$ の項)。横成分は $b_x + ib_y = -2i\partial_{\bar w}A_s + i\,\partial_sa_w$ ($(6.4)$ の $b_x + ib_y = \partial_yA_s - i\partial_xA_s - \partial_sa_y + i\partial_sa_x$)。$A_s = -\kappa(w^m + \bar w^m)/(2m!)$ から $-2i\partial_{\bar w}A_s = i\kappa\bar w^{m-1}/(m-1)!$ で、これは $\Phi$ の第 1 項の $2\partial_{\bar w}\Phi_0 = \partial_x\Phi_0 + i\partial_y\Phi_0$ に等しい。$i\partial_sa_w = \dfrac{i\kappa''}{4m!}\Bigl(\dfrac{w^{m+1}}{m+1} - w\bar w^m\Bigr)$ は、$2\partial_{\bar w}\Phi_2 = \dfrac{i\kappa''}{4(m+1)!}\bigl(w^{m+1} - (m+1)w\bar w^m\bigr)$ に等しい。よって $\nabla\times\vec A = \nabla\Phi$ ($\kappa'$、$\kappa''$ の次数まで)。$b_s$ の $\kappa'''$ の項 ($\partial_s\Phi_2$) は、$A$ の次の次数の補正が要り、この章では使わない (残差は sympy で $k_3$ の項だけと確認。edge_01 (c))。
- ソレノイド: 円筒座標で $b_z = r^{-1}\partial_r(rA_\phi)$、$b_r = -\partial_sA_\phi$。$A_\phi = rb/2 - r^3b''/16$ から $b_z = r^{-1}(r^2b/2 - r^4b''/16)' = b - r^2b''/4$、$b_r = -rb'/2 + r^3b'''/16$。(6.3) に一致。

### 4.2 $H$ の展開

$s$ を独立変数とする $H$ は (曲率と $E/(p_0v_0)$ の定数項を除く。01a 章)

$$
H = -\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2} - A_s
\tag{6.5}
$$

$\vec a_\perp \to \epsilon\,\vec a_\perp$ として $\epsilon$ で展開する。根号の中は ($p_z^2 \equiv p^2 - p_x^2 - p_y^2$)

$$
u(\epsilon) = p^2 - (p_x - \epsilon a_x)^2 - (p_y - \epsilon a_y)^2 = p_z^2 + 2\epsilon\,(\vec p_\perp\cdot\vec a_\perp) - \epsilon^2a_\perp^2
$$

$-\sqrt u = -p_z\bigl[1 + X\bigr]^{1/2}$、$X = \bigl(2\epsilon\,\vec p_\perp\cdot\vec a_\perp - \epsilon^2a_\perp^2\bigr)/p_z^2$ に $(1 + X)^{1/2} = 1 + X/2 - X^2/8 + \cdots$ を入れる:

$$
-p_z\Bigl[1 + \frac{\epsilon\,p_\perp\cdot a_\perp}{p_z^2} - \frac{\epsilon^2a_\perp^2}{2p_z^2} - \frac{4\epsilon^2(p_\perp\cdot a_\perp)^2}{8p_z^4} + O(\epsilon^3)\Bigr]
= -p_z - \epsilon\frac{p_\perp\cdot a_\perp}{p_z} + \epsilon^2\left[\frac{a_\perp^2}{2p_z} + \frac{(p_\perp\cdot a_\perp)^2}{2p_z^3}\right] + O(\epsilon^3)
$$

($\epsilon^2$ の 2 つの項は、$X$ の 1 次の項の $\epsilon^2$ と、$X^2/8$ の $\epsilon^2$ から。sympy の 1 次・2 次の微分と一致。edge_05 (e)。) よって ($\epsilon = 1$)

$$
H = \underbrace{-p_z}_{\text{自由}} \;\underbrace{-\,A_s}_{\text{本体}} \;+\; H_1 + H_2 + O(a_\perp^3),\qquad
H_1 = -\frac{p_xa_x + p_ya_y}{p_z},\quad H_2 = \frac{a_\perp^2}{2p_z} + \frac{(\vec p_\perp\cdot\vec a_\perp)^2}{2p_z^3},\quad p_z = \sqrt{p^2 - p_x^2 - p_y^2}
\tag{6.6}
$$

($H_2$ は $\kappa'$ の 2 乗で、階段の端では $\delta(s)^2$ になる。§4.4。)

「$y^4$ まで」とは、(6.1) の $b_y$ の $y^4 g''''/24$ ($a_x$ の $y^4g'''/24$) までを入れること。ヘルプの偏向の線形フリンジの生成関数 $V$ にある $y^4$ の項 ($y^4/(6\,\mathrm{F1}\,\rho_b^2)$) は、$1/\rho_b$ の 2 次の摂動の項で、05 章 §6.3 が扱う。この章の 1 次の $f$ に出る $y^4$ の項は、$g'''$ から来る $f_4$ (§4.3)。

### 4.3 階段の端の極限と生成関数

**摂動の 1 次** ($H = H_0 + H_1(s)$、$H_0 = -p_z$ が自由な運動。ここで $H_1$ は小さい)。相互作用表示で書くと、始点 $s = 0$ から終点 $S$ までの転送は、自由な運動 $M_0(S)$ を最後に作用させた形

$$
\mathcal M(S) = M_0(S)\circ\exp(:f:) + O(H_1^2),\qquad f = -\int_0^SH_1\bigl(M_0(s)\,z_0\bigr)\,ds
\tag{6.7}
$$

になる。ここで $M_0(s)z_0$ は始点 $z_0$ から自由な運動で $s$ 進んだ位置と運動量 ($x(s) = x + \dfrac{p_x}{p_z}s$、$y(s) = y + \dfrac{p_y}{p_z}s$、$p_\perp$ 一定)。つまり、摂動の各点での値を自由な運動で始点に引き戻して足し、約束 (6.0) の $f = -\int H\,ds$ とする。

**(6.7) の導出** (一般の $H = H_0 + H_1(s)$。$H_1$ は小さいとする。確認: [edge_06_derivations.out.txt](outputs/edge_06_derivations.out.txt) (a))。次の 6 段で、各行に何をしたかを添える。

(1) **記号**。始点 $z_0$ から $s$ まで進んだ点を $\mathsf M(s)z_0$、関数 $g$ への作用を「終点で評価する」引き戻し $(\mathcal M(s)\,g)(z_0) = g\bigl(\mathsf M(s)z_0\bigr)$ と書く。$\exp(:f:)$ も同じ意味 ($g \mapsto g + [f, g] + \tfrac12[f, [f, g]] + \cdots$。$:f:$ が微分演算なので、点の流れ $z \mapsto z + [f, z] + \cdots$ を $g$ に代入したものに等しい)。

(2) **$\mathcal M$ の微分方程式**。正準方程式 $\dot z = [z, H]$ (約束 (6.0) と同じ: $x' = \partial H/\partial p_x = [x, H]$) から、粒子に沿った任意の関数 $g$ の変化は $\dfrac{d}{ds}g(z(s)) = [g, H] = -{:}H{:}\,g$。長さ $\mathrm d s$ の 1 歩は $\exp(-\mathrm d s\,{:}H(s){:})$ で、その前に $s$ まで進んでいるので

$$
\mathcal M(s + \mathrm d s)\,g = \mathcal M(s)\Bigl[e^{-\mathrm d s\,:H(s):}g\Bigr]
\quad\Longrightarrow\quad
\frac{d\mathcal M}{ds} = \mathcal M(s)\,\bigl(-{:}H(s){:}\bigr)
\tag{6.7a}
$$

(「点 $z(s)$ から 1 歩進めた値」$= \bigl(e^{-\mathrm ds:H:}g\bigr)(z(s)) = \mathcal M(s)\bigl[e^{-\mathrm ds:H:}g\bigr](z_0)$。右辺の $\mathrm d s$ の 1 次を取った。)

(3) **自由な運動を分ける**。$H_0$ だけの転送 $\mathcal M_0(s)$ は $d\mathcal M_0/ds = \mathcal M_0(-{:}H_0{:})$、$\mathcal M_0(0) = 1$。$\mathcal M = \mathcal N\,\mathcal M_0$ とおく ($\mathcal N(0) = 1$)。(6.7a) の左辺は $\mathcal N'\mathcal M_0 + \mathcal N\mathcal M_0(-{:}H_0{:})$、右辺は $\mathcal N\mathcal M_0\bigl(-{:}H_0{:} - {:}H_1{:}\bigr)$ なので、$-\mathcal N\mathcal M_0{:}H_0{:}$ が両辺で消えて

$$
\mathcal N'\,\mathcal M_0 = -\,\mathcal N\,\mathcal M_0\,{:}H_1{:}
\quad\Longrightarrow\quad
\mathcal N' = -\,\mathcal N\;\mathcal M_0{:}H_1{:}\mathcal M_0^{-1}
$$

(右から $\mathcal M_0^{-1}$ をかけた)。

(4) **$\mathcal M_0$ をポアソン括弧の中へ入れる**。$\mathsf M_0$ は正準変換 (自由な運動) なので、括弧が保たれる: $[\mathcal M_0a, \mathcal M_0b] = \mathcal M_0[a, b]$ (sympy で、自由な運動 $x \to x + p_xs$ に対して確認: edge_06 (a))。これを使うと、任意の $b$ について $\mathcal M_0{:}H_1{:}\mathcal M_0^{-1}b = \mathcal M_0[H_1, \mathcal M_0^{-1}b] = [\mathcal M_0H_1,\ b]$、すなわち $\mathcal M_0{:}H_1{:}\mathcal M_0^{-1} = {:}h{:}$。ここで $h(s; z_0) \equiv (\mathcal M_0H_1)(z_0) = H_1\bigl(\mathsf M_0(s)z_0,\ s\bigr)$ (自由な運動で $s$ 進んだ位置と運動量を $H_1$ に代入したもの。$x(s) = x + p_xs/p_z$ など)。よって $\mathcal N' = -\mathcal N\,{:}h(s){:}$。

(5) **1 次で積分**。$\mathcal N(0) = 1$ から、$H_1$ の 1 次まで $\mathcal N(S) = 1 - \displaystyle\int_0^S{:}h(s){:}\,ds + O(H_1^2) = \exp(:f:) + O(H_1^2)$、$f = -\displaystyle\int_0^Sh\,ds$ ($:f:$ は $f$ について線形なので積分が中に入る)。$O(H_1^2)$ は、異なる $s$ の ${:}h(s){:}{:}h(s'){:}$ の積 (一般に非可換。Magnus 展開の 2 次) で、$\delta$ 関数の端では $\delta(s)^2$ になる (§4.4)。

(6) **点の写像に戻す**。$\mathcal M g = \mathcal N(\mathcal M_0g)$ は、$(\mathcal Mg)(z_0) = (\mathcal M_0g)\bigl(\mathsf N z_0\bigr) = g\bigl(\mathsf M_0(S)\,\mathsf N z_0\bigr)$ ($\mathsf N$ は $\exp(:f:)$ の流れ)。つまり**点としては $\mathsf N$ が先、自由な運動 $\mathsf M_0(S)$ が後**: $\mathsf M = \mathsf M_0(S)\circ\exp(:f:)$。これが (6.7) の順序で、**$f$ の引数は始点の変数 $z_0$** ($h$ は始点から自由に進んだ点での値)。もし $\exp(:f:)$ を $\mathsf M_0$ の後に置く ($f$ を終点の変数で書く) と、$f$ が別の関数になる ($f_{\rm end} \ne f$、edge_06 (a) の参考行)。階段の端が始点 $s = 0$ にあるときは $\mathsf M_0(0) = 1$ なので順序は効かない。出口の端は、出口の面を始点にとった長さ 0 の要素 (自由な運動の長さが 0 なので $\mathsf M_0 = 1$) として扱い、本体 (と本体までの自由な運動) の写像のあとに合成する。(6.8) の $f_{\rm out}$ はその要素の生成関数で、出口の面での変数で書く。

**確認** (edge_06 (a)): 2 自由度 $H_0 = (p_x^2 + p_y^2)/2$、$H_1 = \epsilon\,x^2y\,(1 + s + s^2)$ で、厳密な 1 次の応答 ($\delta x' = \delta p_x$、$\delta p_x' = -\partial_xH_1(\text{自由軌道})$、…を $s$ で積分) と、(6.7) の $f$ から (6.0a) で作った $\mathsf M_0(S)\circ\exp(\epsilon:f:)$ の $\epsilon$ の 1 次が、$x, p_x, y, p_y$ の 4 成分とも **差 0** (sympy)。

**検算** (1 自由度のおもちゃ、$H_0 = p^2/2$、$H_1 = \epsilon x^2h(s)$、$h = s + s^2$、$0 \le s \le S$。edge_05 (f)): 厳密な 1 次の応答は $x_1'' = -2x_0h$、$x_0 = q_0 + p_0s$ から $\Delta x(S) = \int_0^S(S - s)(-2x_0h)\,ds$、$\Delta p(S) = \int_0^S(-2x_0h)\,ds$。(6.7) の $f = -\int_0^S(q_0 + p_0s)^2h\,ds$ から (6.0a) で $\Delta x_f = -\partial f/\partial p_0 = \int2(q_0 + p_0s)s\,h\,ds$、$\Delta p_f = \partial f/\partial q_0 = -\int2(q_0 + p_0s)h\,ds$、出力 $= (q_0 + \Delta x_f + S(p_0 + \Delta p_f),\ p_0 + \Delta p_f)$。自由な運動の部分を除くと $\Delta x_f + S\Delta p_f = \int2(q_0 + p_0s)h\,(s - S)\,ds$ で、厳密な 1 次の応答に一致 (sympy で差 0)。

階段の端 $\kappa' = \kappa_0\,\delta(s)$ では、$H_1$ が $\delta(s)$ に比例し、積分範囲の $s = 0$ だけが効くので ($M_0(0) = $ 恒等、$\vec a_\perp = \kappa'\tilde a_\perp$ とおく)

$$
f = +\kappa_0\,\frac{p_x\tilde a_x + p_y\tilde a_y}{p_z}\quad\text{(入口。}\vec a_\perp = \kappa'\tilde a_\perp\text{)},\qquad f_{\rm out} = -f_{\rm in}
\tag{6.8}
$$

(出口は $\kappa' = -\kappa_0\,\delta(s - s_{\rm out})$ なので符号が逆。$f = -\int H_1\,ds = +\int(p\cdot a)/p_z\,ds$ で、$\int\kappa'ds = \kappa_0$。)

$\delta'$、$\delta''$ の項 (ソレノイドの $b'' = b_0\delta'$、偏向の $y^4$ の $g''' = \delta''$) は、$n$ 回の部分積分 (端の外では $\delta$ が 0 なので境界項は消える) で

$$
\int\delta^{(n)}(s)\,F(s)\,ds = (-1)^nF^{(n)}(0),\qquad F(s)\ \text{は自由な運動に沿った値}
\tag{6.9}
$$

**結果** ([edge_01_gauge_generators.out.txt](outputs/edge_01_gauge_generators.out.txt) (d)〜(f)):

- **偏向**: $\tilde a_x = -y^2/(2\rho_b)$ ($g' = \delta$) なので $H_1 = -p_xa_x/p_z = +p_xy^2g'/(2\rho_bp_z)$、(6.8) から

$$
f = -\frac{y^2}{2\rho_b}\frac{p_x}{p_z}
$$

(6.0a) の微分を 1 行ずつ: $\partial p_z/\partial p_x = -p_x/p_z$、$\partial p_z/\partial p_y = -p_y/p_z$、$\partial p_z/\partial p = p/p_z$ から

$$
\Delta x = -\frac{\partial f}{\partial p_x} = \frac{y^2}{2\rho_b}\frac{\partial}{\partial p_x}\frac{p_x}{p_z} = \frac{y^2}{2\rho_b}\Bigl(\frac1{p_z} + \frac{p_x^2}{p_z^3}\Bigr) = \frac{y^2}{2\rho_b}\frac{p^2 - p_y^2}{p_z^3}
$$

$$
\Delta p_y = +\frac{\partial f}{\partial y} = -\frac{y}{\rho_b}\frac{p_x}{p_z},\qquad
\Delta z = -\frac{\partial f}{\partial p} = \frac{y^2}{2\rho_b}p_x\frac{\partial}{\partial p}\frac1{p_z} = -\frac{y^2}{2\rho_b}\frac{p\,p_x}{p_z^3},\qquad \Delta p_x = \frac{\partial f}{\partial x} = 0
$$

$p_y = 0$ で (5.35)。
- **四極**: $\tilde a_x = -(x^3 + 3xy^2)/12$、$\tilde a_y = (3x^2y + y^3)/12$、$\kappa_0 = K$ ($p_z \to p$ の近軸の式)。(6.8) から

$$
f = -\frac{K}{12p}\bigl[p_x(x^3 + 3xy^2) - p_y(3x^2y + y^3)\bigr] = -\bigl(p_xV_x + p_yV_y\bigr),\qquad V_x = \frac{K}{12p}(x^3 + 3xy^2),\ V_y = -\frac{K}{12p}(3x^2y + y^3)
$$

((4.27)。差 0。)$\Delta x = -\partial f/\partial p_x = V_x$、$\Delta y = V_y$ なので、位置のずれは $V$ そのもの。(6.8) を厳密な $p_z$ のまま使うと (4.27) の $p \to p_z$ に当たる ($p_z \to p$ の近軸が (4.27))。

- **$2m$ 極**: $f = \kappa_0\,\mathrm{Re}\bigl(p_w\overline{\tilde a_w}\bigr)/p$ ($p_w = p_x + ip_y$)。$\Delta w = -\partial_{p_x}f - i\,\partial_{p_y}f = -\dfrac{\kappa_0}{p}\tilde a_w$ ((7.17)。$m = 2, 3, 4$ で差 0)。
- **ソレノイド**: $H_1 = -p_\phi\,(b/2 - r^2b''/16)/p_z$ ($p_\phi = xp_y - yp_x$。$\vec p_\perp\cdot(-y, x) = p_\phi$)。$b/2$ の項は本体の回転。$b''$ の項は $H_1^{(b'')} = +p_\phi r^2b''/(16p_z)$ で、$b'' = b_0\delta'$ に (6.9) ($n = 1$) を使う: $F(s) = p_\phi r(s)^2/(16p_z)$。自由な運動で $p_\phi(s) = (x + p_xs/p_z)p_y - (y + p_ys/p_z)p_x = xp_y - yp_x$ は保存、$r(s)^2 = x^2 + y^2 + 2s\,p_r/p_z + (p_\perp s/p_z)^2$ ($p_r = xp_x + yp_y$) から $dr^2/ds|_0 = 2p_r/p_z$。よって $F'(0) = p_\phi p_r/(8p_z^2)$、$\int\delta'F\,ds = -F'(0)$、

$$
f = -\int H_1^{(b'')}ds = -b_0\int\delta'(s)\,F(s)\,ds = +b_0F'(0) = +\frac{b_0}{8\,p_z^2}\,p_\phi\,p_r\ \ \xrightarrow{\ p_z\to p\ }\ +\frac{b_0}{8p^2}\,p_\phi p_r
\tag{6.10}
$$

(ここで $H_1^{(b'')} = b''F(s)$、$F(s) = p_\phi r(s)^2/(16p_z)$ なので $\int H_1^{(b'')}ds = b_0\int\delta'F\,ds = -b_0F'(0)$。) ヘルプの DISFRIN の $H = -\dfrac{B_z}{8p^2B\rho}p_\phi p_r$ は、$H_{\rm help} = -f$ に当たる (§2。K64 は paraxial の $p^2$)。

- **偏向の $y^4$ の項** ($f_4$): $a_x$ の $y^4$ の項は $+y^4g'''/(24\rho_b)$ ($g''' = \delta''$)。$H_1^{(4)} = -p_xa_x^{(4)}/p_z = -p_xy^4g'''/(24\rho_bp_z)$ で、$f_4 = -\int H_1^{(4)}ds = +\dfrac{p_x}{24\rho_bp_z}\int\delta''(s)\,y(s)^4ds = +\dfrac{p_x}{24\rho_bp_z}\dfrac{d^2}{ds^2}\bigl[y(s)^4\bigr]_0$ ((6.9)、$n = 2$)。$y(s) = y + (p_y/p_z)s$ で $\dfrac{d^2}{ds^2}y(s)^4\big|_0 = 12y^2(p_y/p_z)^2$ なので

$$
f_4 = +\frac{p_xp_y^2y^2}{2\rho_bp_z^3}
$$

($f_2 = -\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$ に $y^2$ の項は $\delta$ 自身で $F(0)$ のみ。)$p_y$ の展開: $1/p_z = p_v^{-1}(1 - p_y^2/p_v^2)^{-1/2} = p_v^{-1}\bigl(1 + p_y^2/(2p_v^2) + \cdots\bigr)$ ($p_v = \sqrt{p^2 - p_x^2}$、$p_z^{-3} \approx p_v^{-3}$ ($p_y^2$ まで))。

$$
f_2 + f_4 = -\frac{y^2p_x}{2\rho_bp_v}\Bigl(1 + \frac{p_y^2}{2p_v^2}\Bigr) + \frac{y^2p_xp_y^2}{2\rho_bp_v^3} = -\frac{y^2}{2\rho_b}\frac{p_x}{p_v}\left(1 - \frac{p_y^2}{2p_v^2}\right) + O(p_y^4)
$$

(sympy で差 0。edge_01 (e)。)**K64 の (5.36) は $p_y$ を $p_z$ に入れず ($p_v$)、この $p_y^2$ の補正を含まない**。補正の大きさは $(p_y/p_v)^2/2$ の相対値 (角の 2 乗)。**この補正は導出しただけで、ローレンツ力の積分では確かめていない** (展開が $y \ll d$ で成り立つ窓では $p_y$ を大きくできない。§14)。

### 4.4 適用範囲

(6.7) は場の強さの 1 次で、2 次の項 $H_2$ ((6.6)) は階段の端 ($d \to 0$) で $\delta(s)^2$ のように発散する。実際の端は有限の幅 $d$ を持ち、$\vec a_\perp \sim \kappa'$ の大きさは $\kappa_0/d$ のオーダー。1 次で十分なのは**端の通過中の運動量のずれが小さい**とき: 偏向では $P_x$ のスパイクの高さ $\sim y^2/(2\rho_bd) \ll 1$ (数値で $y = 2$ mm、$\rho_b = 50$ m、$d = 1$ mm なら $y^2/(2\rho_bd) = 4\times10^{-6}/(100\times10^{-3}) = 4\times10^{-5}$)。また場の展開 (6.1)〜(6.3) は $y \ll d$ ($r \ll d$) でだけ収束する。$d \to 0$ の階段の極限で $y \gg d$ になる領域は、実際の磁石の端の形に依存する部分で、ヘルプ・K64 の式はその領域を表さない。

## 5. 要素ごとの生成関数 (結果の整理)

| 要素 | 端の項 | $f_{\rm in}$ (出口は符号反転。約束は §2) | 本章・他章 |
|---|---|---|---|
| 偏向 (K1 = 0) | $a_x = -y^2g'/(2\rho_b)$ | $-\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_z}$ | §4.3、05 章 (5.34)〜(5.36) |
| 四極 | $a_\perp = K'(\ldots)$ | (4.27) ($= -(p_xV_x + p_yV_y)$) | 04 章 §7.2 |
| $2m$ 極 ($m \ge 3$) | $a_w = \kappa'(\ldots)$ | $\Delta w$: (7.17) | 07 章 §9 |
| ソレノイド | $a = (-y, x)(-r^2b''/16)$ | $+\dfrac{b_0}{8p^2}p_\phi p_r$ | (6.10)、09 章 |
| CAVI | $E_r = -\tfrac r2Eg'$ | $f = -H_f$: §9.1 | |
| MULT ($n = 0$、双極子) | 偏向の式と同じ形 | `ttfrin` nord = 2: $\Delta x = \tfrac12 aki\,y^2$ | §9.2 |

## 6. 生成関数から閉じた写像へ (母関数の分解)

$f_{\rm in}$ (6.8) を四極で書くと $f = -(p_xV_x + p_yV_y)$ ((6.0a) で $\Delta x = V_x$、$\Delta y = V_y$)。**運動量について 1 次の生成子**は、位置の点変換 (ベクトル場 $V$ の流れ) と、運動量の逆転置による更新を作る。$\exp(:f:)$ の位置の部分を級数で書く: $\mathcal D \equiv V\cdot\nabla = V_x\partial_x + V_y\partial_y$ とおくと、

$$
:f:\,x = [f, x] = -\frac{\partial f}{\partial p_x} = V_x,\qquad :f:\,V_x = [f, V_x] = \sum_i\Bigl(\frac{\partial f}{\partial x_i}\frac{\partial V_x}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial V_x}{\partial x_i}\Bigr) = \sum_iV_i\,\partial_iV_x = \mathcal DV_x
$$

($V$ は $p_i$ に依らない ($p$ はパラメータ) ので第 1 項は 0)。同様に $:f:^n x = \mathcal D^{n-1}V_x$。よって

$$
\vec x_{\rm out} = \exp(:f:)\,\vec x = \vec x + V + \tfrac12(V\cdot\nabla)V + \tfrac16(V\cdot\nabla)^2V + \cdots
\tag{6.11}
$$

となり、閉じた式にならない (3 次以上が無限に続く)。K64 とヘルプの四極の母関数は、(6.11) の**2 次までで打ち切って正準にする**もの:

$$
G = H_0 + \tfrac12\left(\frac{\partial H_0}{\partial x_1}\frac{\partial H_0}{\partial p_{x2}} + \frac{\partial H_0}{\partial y_1}\frac{\partial H_0}{\partial p_{y2}}\right),\quad H_0 = p_{x2}V_x + p_{y2}V_y\ (= -f)
\tag{6.12}
$$

型 2 の母関数 $F_2 = x_1p_{x2} + y_1p_{y2} + G$ から、$\vec x_2 = \partial F_2/\partial\vec p_2$、$\vec p_1 = \partial F_2/\partial\vec x_1$。$H_0$ は $p_2$ について 1 次なので $\partial H_0/\partial p_{x2} = V_x(x_1)$、$\partial H_0/\partial x_1 = p_{x2}\partial_xV_x + p_{y2}\partial_xV_y$ で、$G$ の第 2 項は

$$
\tfrac12\Bigl[(p_{x2}\partial_xV_x + p_{y2}\partial_xV_y)V_x + (p_{x2}\partial_yV_x + p_{y2}\partial_yV_y)V_y\Bigr] = \tfrac12\,\vec p_2\cdot(V\cdot\nabla)V
$$

つまり $G = \vec p_2\cdot\vec W(\vec x_1)$、$\vec W = V + \tfrac12(V\cdot\nabla)V$。これから

- 位置: $\vec x_2 = \vec x_1 + \partial G/\partial\vec p_2 = \vec x_1 + \vec W$。$G$ の第 2 項の $\partial/\partial p_{x2}$ が $\tfrac12(V_x\partial_xV_x + V_y\partial_yV_x)$。
- 運動量: $\vec p_1 = \vec p_2 + \partial G/\partial\vec x_1 = \vec p_2 + (\partial_{x_1}\vec W)^\top\vec p_2 = J^\top\vec p_2$ ($J_{ij} = \partial(x_2)_i/\partial(x_1)_j$)。逆に解いて $\vec p_2 = J^{-\top}\vec p_1$ (運動量の逆転置)。

**sympy で、四極の $V$ (4.27) について $x_2 - (x + V + \tfrac12 V\cdot\nabla V) = 0$、$p_1 - J^\top p_2 = 0$ と、(4.29) の $x(1 + a/3 + b + t)$、$y(1 - a - b/3 + t)$ ($a = Kx^2/(4p)$、$b = Ky^2/(4p)$、$t = (a - b)^2/6$) との一致を確認した** ([edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (a))。手計算 ($c \equiv K/(12p)$、$V_x = c\,x(x^2 + 3y^2)$、$V_y = -c\,y(3x^2 + y^2)$): $V_x = x(a/3 + b)$ ($a/3 = Kx^2/(12p)$、$b = 3Ky^2/(12p)$)。$(V\cdot\nabla)V_x = V_x\,c(3x^2 + 3y^2) + V_y\,6c\,xy = c^2\bigl[(x^3 + 3xy^2)(3x^2 + 3y^2) - 6xy(3x^2y + y^3)\bigr] = c^2\,(3x^5 - 6x^3y^2 + 3xy^4) = 3c^2x(x^2 - y^2)^2$。$a - b = 3c\,(x^2 - y^2)$ なので $t\,x = \dfrac{9c^2(x^2 - y^2)^2}{6}\,x = \tfrac12\cdot3c^2x(x^2 - y^2)^2 = \tfrac12(V\cdot\nabla)V_x$。$y$ 成分も同じ手順で (edge_06 (b))。$\partial_xV_y = -6c\,xy$、$\partial_yV_y = -c(3x^2 + 3y^2)$ なので

$$
(V\cdot\nabla)V_y = V_x\,\partial_xV_y + V_y\,\partial_yV_y = c^2\bigl[-6x^2y(x^2 + 3y^2) + 3y(3x^2 + y^2)(x^2 + y^2)\bigr]
$$

(第 1 項は $c\,x(x^2 + 3y^2)\cdot(-6c\,xy)$、第 2 項は $\bigl(-c\,y(3x^2 + y^2)\bigr)\bigl(-c(3x^2 + 3y^2)\bigr)$)。展開して $(3x^2 + y^2)(x^2 + y^2) = 3x^4 + 4x^2y^2 + y^4$ を使うと $c^2y\bigl[(-6x^4 - 18x^2y^2) + (9x^4 + 12x^2y^2 + 3y^4)\bigr] = 3c^2y\,(x^4 - 2x^2y^2 + y^4) = 3c^2\,y\,(x^2 - y^2)^2$。$y$ の更新 $V_y = -c\,y(3x^2 + y^2) = -y(a + b/3)$ (4.29) に加わる $\tfrac12(V\cdot\nabla)V_y = \tfrac32c^2y(x^2 - y^2)^2$ は、$t\,y = \tfrac{9c^2(x^2 - y^2)^2}{6}\,y$ に等しい (符号は $x$ と同じ正)。つまり $t\,(x, y) = \tfrac12(V\cdot\nabla)V$ ($K$ の 2 次)。

- 点変換 (位置が運動量によらない) に直した写像は、運動量を $J^{-\top}$ で更新すれば、**構成からシンプレクティック**。
- $z$ の更新。$V \propto 1/p$、$\tfrac12(V\cdot\nabla)V \propto 1/p^2$ なので $G = \vec p_2\cdot\vec W_1/p + \vec p_2\cdot\vec W_2/p^2$ ($W_1 = pV$、$W_2 = \tfrac12p^2(V\cdot\nabla)V$ は $p$ に依らない)。(6.0a) の $\Delta z = -\partial f/\partial p$ に相当する母関数の形は $z_2 = z_1 + \partial G/\partial p$ で、$\vec p_2$ を固定して

$$
\frac{\partial G}{\partial p} = -\frac{\vec p_2\cdot\vec W_1}{p^2} - \frac{2\,\vec p_2\cdot\vec W_2}{p^3} = -\frac1p\,\vec p_2\cdot\Bigl[V + 2\cdot\tfrac12(V\cdot\nabla)V\Bigr]
$$

つまり 1 次の部分が $-1/p$ 倍、2 次の部分が $-2/p$ 倍。K64 の `dz1 = -(px(dx1 + t xi) + py(dy1 + t yi))/pr` は、`dx1` が $t$ を含む ($= V_x + t\,x$、つまり $W_x$) ので `dx1 + t*xi` $= V_x + 2t\,x$ となり、2 次の部分が 2 回入っていて、これと一致する ($\delta = 0, 0.01$ の 2 点で $\Delta z$ が $\partial G/\partial p$ と表示の 7 桁で一致し、「$-p_2\cdot(x_2 - x_1)/p$」とは $K$ の 2 次だけ違う: 同 (c))。$2m$ 極 ($m \ge 3$、`ttfrin` の一般の枝) は $V$ が 1 つの次数だけなので $z = -\mathrm{Re}(\Delta w\cdot\mathrm{cp1})$ (07 章 (7.18))。
- **厳密な $\exp(:f:)$ との差**: 位置は、$V$ が運動量に依らないので $\dot x = V(x)$ を厳密に解く流れで、(6.11) の無限級数そのもの。K64 の式との差の主要項は第 3 項 $\tfrac16(V\cdot\nabla)^2V$ (運動量の差は $J^{-\top}$ の第 3 次の項から)。数値積分 (DOP853、$f$ のハミルトン方程式 (6.0a)、$\tau = 0 \to 1$。入射 $(x, p_x, y, p_y) = (4, 3, -3, 2)$ mm・mrad、$\delta = 0$) と K64 の式 (`ttfrin` 4 次) の差と、$K$ を 2 倍にしたときの比 ([edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) (b)):

| $K$ [/m²] | $x$ の差 | $p_x$ の差 | $z$ の差 | 比 ($x$, $p_x$, $z$) |
|---|---|---|---|---|
| 1000 | $4.2\times10^{-11}$ | $1.6\times10^{-10}$ | $2.1\times10^{-13}$ | 8.03, 7.91, 7.88 |
| 2000 | $3.4\times10^{-10}$ | $1.2\times10^{-9}$ | $1.6\times10^{-12}$ | 8.07, 7.81, 7.76 |
| 4000 | $2.7\times10^{-9}$ | $9.4\times10^{-9}$ | $1.2\times10^{-11}$ | 8.13, 7.61, 7.50 |
| 8000 | $2.3\times10^{-8}$ | $6.8\times10^{-8}$ | $8.3\times10^{-11}$ | 8.27, 7.19, 6.93 |
| 16000 | $1.9\times10^{-7}$ | $4.2\times10^{-7}$ | $4.6\times10^{-10}$ | 8.53, 6.24, 5.52 |

$K \le 2000$ では $x$・$p_x$・$y$・$p_y$・$z$ の全成分の比が 7.7〜8.1 で、差は $K^3$ に比例する。$K \ge 8000$ で $p_x$・$z$ の比が 8 を下回るのは、$K^4$ 以上の項が $K^3$ の項と逆符号で混じるため ($K^4$ 以上の項と考えられる。個別の係数は確認していない)。$x$ の差の主要項 $\tfrac16(V\cdot\nabla)^2V_x$ との比は、$K = 1000$ で $1.008$、$4000$ で $1.034$、$16000$ で $1.14$ ($y$: 1.005、1.021、1.078)。ただし**この 2 次の項 ($t$) 自体は、1 次の Maxwell の導出 (4.24) の範囲では物理的な意味を持たない** (母関数が 2 次を保つのは、点変換を厳密に正準にするため)。
- $m \ge 3$ の多極は (7.17) の 1 次の項だけ ($\tfrac12(V\cdot\nabla)V$ を足さない)。差 ($\kappa_0^2$ で減る) が出るのは四極だけ (07 章 §9.2)。

この「リー変換 → 点変換 + 逆転置」の形が、**母関数の分解**の中身: リー生成子 $f = -p\cdot V$ の流れを、位置の点変換 (写像の 2 次までの級数) と運動量の逆転置の 2 段に分ける。この「2 次で止める」ことが、生成子を $f + O(K^3)$ の別の生成子に置き換えることに当たる、を式で示す (edge_06 (c))。

**K64 の写像は別の生成子の厳密な流れ**。K64 の写像 $\mathsf K$ は、位置が $x \mapsto \Phi(x) = x + V + \tfrac12(V\cdot\nabla)V$、運動量が $p \mapsto J^{-\top}p$ ($J = \partial\Phi/\partial x$) という点変換。$K$ ($V \propto K$) の級数で恒等写像から出発する点変換は、ベクトル場 $U = KU_1 + K^2U_2 + K^3U_3 + \cdots$ の時間 1 の流れとして一意に書ける (対数)。位置の流れの級数は、$\mathcal D_U = U\cdot\nabla$ として

$$
\Phi(x) = x + U + \tfrac12\mathcal D_UU + \tfrac16\mathcal D_U^2U + \cdots
\tag{6.12a}
$$

(前の (6.11) と同じ形で、$V$ を $U$ に替えたもの)。$f' = -\vec p\cdot U(\vec x)$ とおくと、(6.0a) から $\Delta x = -\partial f'/\partial p = U$ で位置は $\dot x = U(x)$ の流れ、運動量は $\dot p = +\partial f'/\partial x = -(\partial U)^{\top}p$。流れの $J(\tau)$ は $\dot J = (\partial U)J$ を満たすので $\tfrac{d}{d\tau}J^{-\top} = -(\partial U)^{\top}J^{-\top}$、したがって $p(\tau) = J(\tau)^{-\top}p_0$。時間 1 で、これが K64 の運動量の逆転置に一致する。すなわち**$\mathsf K = \exp(:f':)$ が厳密に成り立つ** (生成子 $f'$ は点変換から決まる)。

次に $f'$ を $f$ と比べる。(6.12a) の左辺を K64 の $\Phi = x + V + \tfrac12\mathcal DV$ ($\mathcal D = V\cdot\nabla$) とおいて、$K$ の次数ごとに比べる:

- $K^1$: $U_1 = V$。
- $K^2$: $U_2 + \tfrac12\mathcal DV = \tfrac12\mathcal DV$ より $U_2 = 0$。
- $K^3$: $U_3 + \tfrac12\bigl[(V\cdot\nabla)U_2 + (U_2\cdot\nabla)V\bigr] + \tfrac16\mathcal D^2V = 0$ (右辺は K64 の $K^3$ の項 = 0)。$U_2 = 0$ なので $U_3 = -\tfrac16\mathcal D^2V$。

よって

$$
f' = -\vec p\cdot U = f + \tfrac16\,\vec p\cdot(V\cdot\nabla)^2V + O(K^4)
\tag{6.12b}
$$

(sympy: 四極の $V$ で $U = V - \tfrac16\mathcal D^2V$ の流れを $K^4$ まで展開し、K64 の位置 $x + V + \tfrac12\mathcal DV$ との差が $K^1, K^2, K^3$ で 0、$K^4$ で 0 でない。$(V\cdot\nabla)^2V_x = 3c^3x(x^2 - y^2)(5x^4 + 26x^2y^2 + y^4)$、$(V\cdot\nabla)^2V_y = 3c^3y(x^2 - y^2)(x^4 + 26x^2y^2 + 5y^4)$)。つまり K64 の四極のフリンジは、**$f$ ではなく $f + \tfrac16\,p\cdot\mathcal D^2V$ の厳密な $\exp(:\cdot:)$** で、$f$ に対する差は $K^3$ から。厳密な $\exp(:f:)$ の位置が (6.11) の第 3 項 $+\tfrac16\mathcal D^2V$ を持ち、K64 がそれを持たない差 ($-\tfrac16\mathcal D^2V$) が、$f'$ の $+\tfrac16\,p\cdot\mathcal D^2V$ ($:f':x = -\partial f'/\partial p$ の符号を見ると $-\tfrac16\mathcal D^2V$) に当たる。

この置き換えの大きさは、BCH ([03_symplectic_integrators.md](03_symplectic_integrators.md)) の意味では、$\exp(:f + g:)$ ($g = \tfrac16\,p\cdot\mathcal D^2V$) が $\exp(:f:)\exp(:g:)$ と交換子 $\tfrac12[f, g] = O(K^4)$ だけ違う、という話で、$K^3$ の評価には影響しない。物理的には、$K^3$ の項は 1 次の Maxwell の導出 (4.24) の範囲に無いので、$f$ と $f'$ のどちらが正しいかを決める根拠は、この章の導出にない (場の強さの 1 次の摂動 (6.7) は生成子の $K^1$ の項だけを決める。$K^2$ 以上の項は、2 次以上の摂動 (§4.4。階段の端では発散する) を入れないと決まらず、$f$ の累乗 $\tfrac12{:}f{:}^2$ も $f'$ の $K^3$ の項も、その意味では同じく近似である)。

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

05 章 §6.2 に代数的な同値を示した。ここでも 3 行で追える。`fa` $= y/(\rho_e\sqrt{1 - q^2})$、`f` $= a\cdot$`fa`、`ff` $= \tfrac14(f + fa)\,y/(1 - q^2)$ ($a = 1 - y^2/(6\rho_b^2)$、$q = p_x/p$、$\rho_e = \rho_bp$、$1 - q^2 = p_v^2/p^2$) に対し

$$
ff = \frac{y^2}{4\rho_e}\frac{1 + a}{(1 - q^2)^{3/2}},\qquad \frac{1 + a}{4} = \frac12\Bigl(1 - \frac{y^2}{12\rho_b^2}\Bigr)
\ \Rightarrow\ ff = \frac{y^2}{2\rho_bp}\Bigl(1 - \frac{y^2}{12\rho_b^2}\Bigr)\frac{p^3}{p_v^3}
= \frac{y^2}{2\rho_b}\Bigl(1 - \frac{y^2}{12\rho_b^2}\Bigr)\frac{p^2}{p_v^3}
$$

で、ヘルプの $y^2(1 - y^2/(12\rho_b^2))\,p^2/(2\rho_b(p^2 - p_x^2)^{3/2})$ に等しい。数値でも再確認した: ランダムな 2000 点 ($\delta \in \pm0.3$、$p_x/p \in \pm0.3$、$y \in \pm10$ mm) で、追跡の式 (`f`、`fa`、`ff` を $q$ 単位から正準に直したもの) とヘルプの式 (入口の $x_2$、$p_{y2}$、$z_2$) の最大差は $2 \times 10^{-20}$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (a))。なお、ヘルプの $p_{y2}$ の式には $p_1$ が 1 つ多く見え、$p_y$ を ($q$ の単位の) $p_y/p$ と読めば一致する (誤植か規格化の省略か、原文からは確定できない。05 章 付録 A)。

### 8.3 因子 $p$ の独立な再確認 (05 章 付録 B-1。gs で確認済み)

**主張**: (5.35)〜(5.36) の $\Delta x = \dfrac{y^2}{2\rho_b}\dfrac{p^2}{p_v^3}$、$\Delta z = -\dfrac{y^2}{2\rho_b}\dfrac{p_xp}{p_v^3}$ は $p^2$、$p$ の因子を持つ。CALC の `tbedgemaxwell` (`tbfrie.f:209`、`:211`) の `cod(1) = cod(1) + bf*pr`、`cod(5) = cod(5) - bf*pxi` ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$) は、それぞれ $p^1$、$p^0$。

**(i) 導出**: §4.3 の $f$ から $\Delta x = -\partial f/\partial p_x$ は $p_x/p_z$ の $p_x$ 微分 $= (p^2 - p_y^2)/p_z^3$ (自由な運動の運動量 $p$ は $p_x$ を変えて $p_z$ を変える)。$\Delta z = -\partial f/\partial p$ は $p_x/p_z$ の $p$ 微分 $= -p\,p_x/p_z^3$。力の積分で独立に導く (edge_05 (g)): 直線座標で運動量 $\vec P$ (運動学的な運動量、$|\vec P| = p$ 一定) は $d\vec P/ds = \vec P\times\vec b/P_z$ に従い、$b_x = 0$、$b_y = \rho_b^{-1}(g - y^2g''/2)$、$b_s = \rho_b^{-1}y\,g'$ (6.1)。$b$ は $p_0$ で規格化した量で ($b = eB/p_0$)、$p$ で割らない。端の項 ($g'$、$g''$ を含む項。$g$ の項は本体) だけ取ると $P_y = 0$ から出て

$$
\frac{dP_x}{ds} = \frac{P_yb_s - P_zb_y}{P_z} = +\frac{y^2}{2\rho_b}g'',\qquad \frac{dP_y}{ds} = \frac{-P_xb_s}{P_z} = -\frac{P_x}{P_z}\frac{y}{\rho_b}g'
$$

$\int g'\,ds = 1$ (0 から 1 への立ち上がり) から $\Delta P_y = -(p_x/p_z)\,y/\rho_b$ ($= \partial f/\partial y$ に一致)。$P_x$ は $\Delta P_x(s) = \dfrac{y^2}{2\rho_b}g'(s)$ のスパイクになり、$\int\Delta P_x\,ds = \dfrac{y^2}{2\rho_b}\int g'ds = \dfrac{y^2}{2\rho_b}$。位置は $dx/ds = P_x/p_z(P_x)$ なので、スパイクの中で変わる分は 1 次で $\Delta x = \dfrac{\partial(P_x/p_z)}{\partial P_x}\int\Delta P_x\,ds = \dfrac{p^2 - p_y^2}{p_z^3}\cdot\dfrac{y^2}{2\rho_b}$ ($\epsilon = y^2/(2\rho_b)$ に $p$ は入らない)。§4.3 の $\Delta x$ と一致。

**(ii) ローレンツ力の積分** (直線座標、Maxwell を $y^4$ まで満たす (6.1)、$g = (1 + \tanh(s/d))/2$、$d = 1, 0.5$ mm の差を $d \to 0$ に外挿、窓 $\pm8d$、$\rho_b = 50$ m、$y = 0.3$ mm、$p_y = 0$。$y \ll d$ の展開が使える範囲): $q = p_x/p = 0.1, 0.3$、$\delta = 0, +0.3, -0.2$ の 6 点で、**$\Delta x$、$\Delta p_y$、$\Delta z$ のどれも $p^2$ 型 (追跡・ヘルプ) との比 (ODE の外挿値 ÷ 式) が $0.9998 \sim 1.0001$**。CALC の更新との比 (**ODE の外挿値 ÷ CALC の式**) は $\Delta x$、$\Delta z$ で **$1.3000$ ($\delta = 0.3$)、$0.8000$ ($\delta = -0.2$)** で、ちょうど $1 + \delta$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (b))。すなわち、CALC の式は正しい値の $1/(1 + \delta)$ 倍。窓は $y \ll d$ にしてあり、(6.1) の展開が収束する範囲での結論。

**(iii) CALC の内部の矛盾**: 解析行列 (`trans`、`pysad.lattice._tbedgemaxwell_mat`) は、$\Delta x$ は $p^2$ 型の、$\Delta z$ は $-bf\,p_x\,p$ ($p^1$ 型) の微分に、$\delta = 0, \pm0.1$ で $10^{-12}$ で一致する。COD の更新 ($p^1$・$p^0$) の中心差分とは、$\delta = 0.1$ で $1.7 \times 10^{-5}$、$\delta = -0.1$ で $2.5 \times 10^{-5}$ 違い ($\delta = 0$ でも $10^{-6}$。$\partial x/\partial\delta$ の成分 (1, 6) が 2 倍違う)。シンプレクティック性 $|M^\top SM - S|$ は、解析行列が 0 (丸め誤差以内)、COD の更新は $1.7 \times 10^{-5}$ ($\delta = 0.1$) で破れ、$p^2$・$p^1$ 型の更新は $1.5 \times 10^{-12}$ ([edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) (c))。

**結論 (確定)**: 追跡の式と解析行列は物理 ($p^2$ 型) に従い、**CALC の COD の更新だけが因子 $p$ を落としている**。1 つの端面での $\Delta x$ のずれは $\delta\,y^2/(2\rho_b)$ 程度で、$y = 10$ mm、$\delta = 0.05$、$\rho_b = 10$ m なら $0.05\times10^{-4}/20 = 2.5\times10^{-7}$ m ($y = 2$ mm、$\delta = 0.05$、$\rho_b = 50$ m なら $2\times10^{-9}$ m。edge_03 (d))。入口と出口は $f$ の符号が逆なので $\Delta x$ のずれも逆符号で、対称な磁石ではかなり相殺する。**gs で確認した**: 4 台の BEND ($|\rho_b| = 5$ m) と QUAD 2 台の 1 周 (3 GeV/c、$y = 4$〜$8$ mm、`NORFSW`) で、CALC の閉軌道と `TrackParticles` の閉軌道の差の $p_x$ は、$\delta = 0$ で $5 \times 10^{-13}$、$\delta = 0.05$ で $8.5 \times 10^{-8}$、$0.2$ で $3.8 \times 10^{-7}$、$-0.1$ で $-1.5 \times 10^{-7}$ ($\delta y^2/(2\rho_b)$ の大きさ)。pysad の k64 は gs と $10^{-15}$ 以下で一致し、v2 (COD の更新を $p^2$・$p^1$ 型に直した) では $\delta = 0.05$ で $3.8 \times 10^{-10}$ に下がる (05 章 付録 B-1、docs/30 の B324)。残りは端の Maxwell 以外の $\delta$ の項 (未特定)。

### 8.4 Maxwell フリンジの $y^2/\rho_b^2$ の補正 (a、b。導出できていない)

`tbedgemaxwell`・`tbendcore` の $a = 1 - y^2/(6\rho_b^2)$ ($\Delta p_y = -a\,p_xy/(\rho_bp_v)$)、$b = \tfrac12 - y^2/(24\rho_b^2)$ ($\Delta x = b\,y^2p^2/(\rho_bp_v^3)$) について、本章で確認できたことと足りないこと。

- **1 次の摂動 (6.8) では出ない**: $\Delta p_y$ の $y^3/\rho_b^3$ の項、$\Delta x$ の $y^4/\rho_b^3$ の項は、$\rho_b^{-1}$ の 3 次の量。(6.8) は $\epsilon = 1/\rho_b$ の 1 次なので $y/\rho_b$ と $y^2/\rho_b$ の項しか持たない。
- **2 次・3 次の摂動は階段の極限で定義できない**: $H_2$ ((6.6)) は $\delta(s)^2$、3 次の項は $\delta(s)^3$ を含み、$d \to 0$ で発散する (§4.4)。端の幅 $d$ が有限の模型が必要で、ヘルプ・K64 はその模型 (端の形) を書いていない。
- **確認できた一致** (導出ではない): $a$ と $b$ は、$\theta = y/\rho_b$ として $\dfrac{\sin\theta}{\theta} = 1 - \dfrac{\theta^2}{6} + \dfrac{\theta^4}{120} - \cdots$、$\dfrac{1 - \cos\theta}{\theta^2} = \dfrac12 - \dfrac{\theta^2}{24} + \dfrac{\theta^4}{720} - \cdots$ の $\theta^2$ までの展開に等しい ($\Delta p_y = -(p_x/p_v)\sin(y/\rho_b)$、$\Delta x = \rho_b(1 - \cos(y/\rho_b))\,p^2/p_v^3$ の $\theta^2$ まで。sympy: edge_05 (h))。K64 には $\theta^4$ 以上が無い。この三角関数の形が物理的に何を表すか (例えば有限の曲率半径の円弧上の幾何) は、確認していない。
- **場の $y$ のべき級数は閉じた形に足せるが、その収束域が限界を決める**。(6.1) は $b_y + i\,b_s = g(s + iy)$ ($g$ の複素平面への解析接続) と書ける (展開すると $g + iy\,g' - \tfrac{y^2}2g'' - \tfrac{i\,y^3}6g''' + \cdots$ で実部が $b_y$、虚部が $b_s$)。これは、$(s, y)$ 平面の調和関数 $\Phi$ ($\vec b = \nabla\Phi$) が中央面で $b_y = g(s)$、$b_s = 0$ を満たす唯一の解析的な場、という条件から出る ($F = \Phi + i\Psi$ が解析的で $F' = b_s - i\,b_y$、$F'(s) = -i\,g(s)$)。例えば $g = \tfrac12(1 + \tanh(s/d))$ は $s = \pm i\pi d/2$ に極を持つので、場は $|y| < \pi d/2$ でだけ有限。$y \gtrsim d$ の領域の場は、磁石の磁極の形 (境界条件) で決まり、端の関数 $g(s)$ だけからは決まらない。これが、$d \to 0$ の階段の極限で $y \gg d$ の粒子 ($y^2/\rho_b^2$ の補正の主な対象) を扱うときの根本的な限界である。
- **足りないもの**: (a) 端の形を指定した場合の、$\epsilon^3$ までの摂動計算 (または厳密な軌道) と、その $d \to 0$ の極限。ローレンツ力の積分で取り出すには、$y \ll d$ の窓 ($y = 0.3$ mm、$\rho_b = 50$ m) で補正が $(y/\rho_b)^2\sim4\times10^{-11}$ 倍と小さすぎ、ODE の精度 (相対 $10^{-4}$ 程度) に届かない。$y^3/\rho_b^3$ の項を数値で取り出すには、(i) $\rho_b \to -\rho_b$ (場の符号) で奇数次だけを残して $H_2$ (偶数次。階段の極限で $y^3/(\rho_b^2d)$ 型に大きくなる項) を除き、(ii) $\rho_b$ を小さく ($y/\rho_b \sim 10^{-2}$) して補正を $10^{-5}$ に上げ、(iii) ODE を相対 $10^{-10}$ 以下で解き、(iv) 端の幅の有限の補正 ($d/\rho_b$ の項) を $d$ の複数の値で外挿する、の 4 つが要る。さらに、窓が $y < \pi d/2$ (tanh の場合) に限られるので、得られる係数は「その端の形での $y \ll d$ の値」で、$y \gg d$ の階段の極限の値とは限らない (上の閉じた形)。この章ではやっていない。(b) $b$ の $\Delta x$ の補正の符号と係数の、独立な確認。この項の大きさは相対 $y^2/(6\rho_b^2)$ ($y = 10$ mm、$\rho_b = 10$ m で $1.7\times10^{-7}$) で、物理の結果への影響は小さい。

## 9. CAVI・SOL・MULT の端部

### 9.1 CAVI

ヘルプ: 入口の端の $\delta$ 関数の $H$ ($\omega$ は角周波数、$t$ は位相の時間、$V = $ VOLT + DVOLT)

$$
H_f = -\frac{eV}{L}\bigl(\sin(\omega t - d\varphi) + \sin d\varphi - \mathrm{offset}\bigr)\frac{x^2 + y^2}{4}\,\delta(s - s_0)\quad(\text{出口は符号反転})
\tag{6.13}
$$

($H_f$ は、ヘルプの流儀の $H$ ($H_{\rm help} = -f$、§2)。$p_0c$ で規格化する。) $d\varphi$ と offset は、空洞の位相と放射による損失 (NORAD でだけ非 0) で決まる。構造は §3.4 の $E_r = -\tfrac r2E\,g'$ と同じ。

**(6.13) の導出** (端の $E_z = E\,g(s)$ から。$c = 1$、$e$ は電荷。確認: [edge_06_derivations.out.txt](outputs/edge_06_derivations.out.txt) (d))。1 つだけ入力が要る。それは、$E(t)$ の符号と大きさを、ヘルプの $dE = -eV\sin(\ldots)$ (CAVI の本体のエネルギー利得) から読む約束 $eE(t) = -\dfrac{eV}{L}\bigl(\sin(\omega t - d\varphi) + \sin d\varphi - \mathrm{offset}\bigr)$ (本体の長さ $L$ にわたる $\int eE_z\,dz = eEL$) である。以下は、各行を式として追う。

**段 1: 端のスカラーポテンシャル**。§3.4 より $E_z = E\,g - \dfrac{r^2}4E\,g'' + \cdots$、$E_r = -\dfrac r2E\,g' + \cdots$ ($g$ は 0 → 1 に立ち上がる)。静電場として $\vec E = -\nabla\phi$ を満たす $\phi(r, s)$ を $r$ のべき級数で組む。$g = G'$ ($G$ は $g$ の積分) として

$$
\phi = -E\,G(s) + \frac{r^2}4E\,g'(s) + O(r^4)
$$

(確認: $-\partial_s\phi = E\,g - \tfrac{r^2}4E\,g''$、$-\partial_r\phi = -\tfrac r2E\,g'$。$\nabla\cdot\vec E = -\nabla^2\phi$ は $r^0$ の項が 0 で、残る $-\tfrac{r^2}4E\,g'''$ は $\phi$ の $r^4$ の項で打ち消す。sympy。) 第 1 項 $-E\,G$ は軸上の加速 (本体と同じ項) で、端に $\delta$ 関数を作るのは**第 2 項**、$\phi_\delta = \dfrac{r^2}4E\,g'(s)$。$E$ が時間変化するときの磁場 ($B_\theta = r\dot E\,g/(2c^2)$、§3.4) は $g$ に比例し、階段の端で段差を作るだけで $\delta$ 関数ではない (本体の項に入る)。

**段 2: $H$ への入り方**。位置依存の電位 $\phi$ のある $s$ 独立でない $H$ は $H = -\sqrt{(E_k - e\phi)^2 - m^2 - p_\perp^2}$ ($E_k$ はエネルギー)。$\phi$ の 1 次で $\partial H/\partial\phi\big|_0 = e\,E_k/p_z = e/\beta_z$ (edge_05 (d)、sympy)。したがって端の項は $H_f = \dfrac e{\beta}\,\phi_\delta$:

$$
H_f = \frac{e\,E(t)}{\beta}\,\frac{x^2 + y^2}4\,g'(s) \;\xrightarrow{\;g' = \delta(s - s_0)\;}\; -\frac{eV}{\beta L}\bigl(\sin(\omega t - d\varphi) + \sin d\varphi - \mathrm{offset}\bigr)\frac{x^2 + y^2}4\,\delta(s - s_0)
$$

($eE(t)$ に上の約束を代入。**$\beta = 1$ でヘルプの (6.13) に一致する**。ヘルプには $1/\beta$ が無い。)

**段 3: 横方向のキック**。$s$ が独立変数で $(x, p_x)$ が正準対なので (Hamilton の方程式 $dp_x/ds = -\partial H/\partial x$) $\delta$ の係数を積分して

$$
\Delta p_x = -\frac{\partial H_f}{\partial x} = -\frac{eE\,x}{2\beta} = +\frac{eV}{2\beta L}\bigl(\sin(\omega t - d\varphi) + \sin d\varphi - \mathrm{offset}\bigr)\,x
$$

$y$ も同じ ($\Delta p_y = -\partial H_f/\partial y$)。これは力の積分と同じ: 横方向の力の時間積分 $e\int E_r\,dt = \dfrac e\beta\int\bigl(-\tfrac x2\bigr)E\,g'\,ds = -\dfrac{eE\,x}{2\beta}$ ($dt = ds/\beta$、$\int g'\,ds = 1$)。静的な場 ($E$ 一定) の中の粒子を数値積分した結果は、$\beta = 0.3, 0.6, 0.95$ で $\Delta p_x\big/\bigl(-eE\,x/(2\beta)\bigr) = 1.0008$ (3 つとも同じ。$E = 10^{-3}$、$d = 10^{-3}$、$x = 10^{-4}$、$m = c = e = 1$。残りの $8\times10^{-4}$ は端の幅 $d$ の有限の補正)。つまり**キックは $1/\beta$ に比例する**。(約束 (6.0a): $f = -H_f$ から $\Delta p_x = +\partial f/\partial x$、同じ符号。)

**段 4: エネルギーの変化**。縦方向は、$t$ を座標、$-E_k$ をその運動量とする正準対にとる ($z$、$\delta$ への換算は段 5)。$s$ が独立変数で $H = -p_z$ のとき、Hamilton の方程式は $dt/ds = \partial H/\partial(-E_k) = (E_k - e\phi)/p_z = 1/\beta$ (時間が進む向きが正になる符号) と $d(-E_k)/ds = -\partial H/\partial t$、すなわち

$$
\frac{dE_k}{ds} = +\frac{\partial H}{\partial t}
$$

$H_f$ の $t$ 依存は $E(t)$ だけなので、$\delta$ 関数の層でのエネルギーの跳びは ($t$ は層の中で変わらない)

$$
\Delta E_k = \frac{\partial H_f}{\partial t} = -\frac{\omega eV}{\beta L}\cos(\omega t - d\varphi)\,\frac{x^2 + y^2}4
$$

(これが「$\partial H_f/\partial t$」の符号の約束: $t$ と $-E_k$ が共役で $dE_k/ds = +\partial H/\partial t$。)直接の確認として、$\phi_\delta$ のエネルギーへの寄与は力の仕事 $e\int E_z\,ds$ の $-\tfrac{r^2}4E\,g''$ の項で $e\cdot\bigl(-\tfrac{r^2}4\bigr)\int E(t)\,g''\,ds = +\dfrac{e\,r^2\dot E}{4\beta}$ ($g'' = \delta'$、(6.9) $n = 1$ で $-F'(0)$ に $F = E(t_0 + s/\beta)$、$F' = \dot E/\beta$)。$e\dot E = -(\omega eV/L)\cos(\ldots)$ で上と一致する。層の前後で $\phi$ は連続 ($G$ は連続) なので、$E_k$ の跳びは運動エネルギーの跳びに等しい。

**段 5: $z$ の更新**。SAD の縦座標は $t$ との関係 $t = -z/\beta = -z\,h/p$ ($h = E_k/(p_0c)$、$p$ は規格化した運動量、$\beta = p/h$。コードの `t = -z*h1/p1`) で決まる。層の中で $t$ は変わらず ($\partial H_f/\partial E_k = 0$ の 1 次まで)、$E_k$ だけが $h_1 \to h_2 = h_1 + dh$、運動量が $p_1 \to p_2$ に変わる。$z = -t\,p/h$ なので

$$
z_2 = -t\,\frac{p_2}{h_2} = z_1\,\frac{p_2/h_2}{p_1/h_1} = z_1\,\frac{p_2}{p_1}\,\frac{h_1}{h_2}
$$

($z\,p/h = -t$ が不変、という言い方もできる。sympy: edge_06 (d))。K64 の `z → z p2/p1 · h1/h2` と同じ。ここで $dh = \Delta E_k/(p_0c)$、コードの `dh = -w*p0*vf*(x²+y²)/2*cos(ph - dphis)` ($vf = V/(2Lp_0)$) は $-\omega V r^2\cos/(4L)$ で段 4 の $\beta = 1$ の値と一致する。

$d\varphi$ と offset は、空洞の位相と放射による損失 (NORAD でだけ非 0) で決まり、導出には入らない (上の約束の中)。構造は §3.4 の $E_r = -\tfrac r2E\,g'$ と同じ。

K64 `tcavfrin` (`tcav.f:209-249`): `vf = v/al/p0*.5`、`s0 = sin(dphis) - offset`、`dpt = vf*(sin(ph - dphis) + s0)/p1r`、`px += x*dpt`、`py += y*dpt`、`dh = -w*p0*vf*(x²+y²)/2*cos(ph - dphis)`、$h_2 = h_1 + dh$ から $\delta$ と $z$ ($z \to z\,p_2/p_1\cdot h_1/h_2$) と $dv$ を更新する。入口は `fringe` かつ `mfring >= 0` かつ `mfring /= 2` (`tcav.f:82`)、出口は `mfring /= 1` (`:185`)。コードを式に直すと ($v$ = 電圧 [V]、`p0` = 規準運動量 $p_0c$ [eV]、`p1r` $= p = 1 + \delta$): `px` は $q = p_x/p$ の変数なので、$\Delta q_x = \Delta P_x/(p_0p) = \dfrac{V}{2Lp_0}(\sin\ldots + s_0)\,\dfrac{x}{p}$ で、`1/p1r` は**この $q$ への変換**。この式は上の $\Delta P_x$ で $\beta = 1$ としたものに一致する。時間は `t = -z*h1/p1` ($= -z/\beta$、`h1/p1` $= 1/\beta$) とエネルギーと運動量から $\beta$ を使って計算しているので、コードは粒子の $\beta$ を持っているが、**`dpt` には $1/\beta$ が掛かっていない**。

**判断**: 「`1/p1r` が速度の代わりになる」という読みは成り立たない: `1/p1r` は $q$ への変換 ($p$ で割る) で、速度 ($v = \beta c$) の代わりにはならない。**$\beta \to 1$ でだけ正しい式** (K64 の誤り。ローレンツ力の数値積分と gs で確認、v2 は修正 = B330。[10 章](10_cavity.md) §7.1) で、低エネルギー (陽子・イオン) では $1/\beta$ だけ小さく見積もる。例えば $\beta = 0.9$ なら $1/\beta - 1 = 11\%$。上の導出 (段 2〜4) は、$\Delta p_x$ と $\Delta E_k$ のどちらにも $1/\beta$ が入ることを示す (`dpt` と `dh` の両方が $\beta = 1$ の式)。**導出で入れた仮定と確認**: (a) $E(t)$ の符号と大きさの読み (ヘルプの $dE = -eV\sin$ から。ヘルプの規格化の記述が $\beta$ を含むかは未確認)、端の層の中で $E(t)$ が変わらない、$B_\theta$ は $\delta$ 関数を作らない、$\phi$ の 1 次 (摂動の 1 次) の近似。(b) $1/\beta$ は、軸の近くの Maxwell の場のローレンツ力の数値積分 (0.6 MeV/c で k64 の収束が 24% 不足) と gs の `TrackParticles` で確認した。空洞全体の ODE との比較 (長い空洞の本体との組合せ) は未。(c) `tcavfrin` が $q = p_x/p$ を再スケールしない非シンプレクティックの欠陥は、10 章 §7.2 (B332。gs で確認、v2 は修正)。電子 ($\beta \approx 1$) では $1/\beta$ の差は出ない。詳しい扱いは [10_cavity.md](10_cavity.md)。

### 9.2 MULT

`tmulti` (`tmulti.f:125-139`、`:253-282`): 入口 `ttfrins(np,..., nord = 2n+2, akr(n), al, bzs)` を $n = 0, \ldots, nmmax$ (`dofr(n)` が真のもの) について、線形フリンジ ($\mathrm{F1}$・$\mathrm{F2}$: `f1in`、`f2in` の $\exp(a)$ 型 ((4.33) と同じ式)、`tblfri` ($\mathrm{FB1}$)) → 本体 → 出口 (逆順、$-\mathrm{akr}$)。`akr(n)` は複素 (実部が正規成分、虚部が skew)。`ttfrins` (`tquad.f:453-500`) は、成分を回転して正規成分にそろえ、`ttfrin` を呼び、逆回転する。この回転の向きを、座標の回転の定義から決める (確認: [edge_06_derivations.out.txt](outputs/edge_06_derivations.out.txt) (e))。

**(1) 複素強度とスカラーポテンシャル**。$\kappa_c = \mathrm{ak}/\mathrm{al} = \mathrm{ak}_1 + i\,\mathrm{ak}_2$ ($\mathrm{ak}_1$ = 正規、$\mathrm{ak}_2$ = skew) に対して、$m$ 極 ($m = \mathrm{nord}/2$) のポテンシャルを

$$
\Phi = \frac{1}{m!}\,\mathrm{Im}\bigl(\kappa_c\,w^m\bigr) = \frac{1}{m!}\bigl[\mathrm{ak}_1\,\mathrm{Im}\,w^m + \mathrm{ak}_2\,\mathrm{Re}\,w^m\bigr]
\tag{6.13a}
$$

と定義する (正規成分 $\mathrm{Im}\,w^m$ は §3.2 と同じ。$\mathrm{ak}_2$ が $+\mathrm{Re}\,w^m$ の係数)。これが `ttfrins` と、本体側の薄い蹴り `tmulk` の約束と合うことは、2 点で確かめる。(i) `ttfrins` の分岐 `ak(1) == 0` (純 skew) は $\theta = \pi/\mathrm{nord} = \pi/(2m)$ で、正規成分 $\mathrm{ak}_2$ を入れる。(ii) 薄い蹴りの `tmulk` (pysad の `_tmulk`) の Horner 形は、$p_x - ip_y$ の変化 $= -\sum_n\mathrm{akr}(n)\,w^n/n!$ ($p_x \mathrel{-}= \mathrm{Re}$、$p_y \mathrel{+}= \mathrm{Im}$)。(6.13a) では $F(w) = \kappa_cw^m/m!$ として $b_y + i\,b_x = F'(w) = \kappa_cw^{m-1}/(m-1)!$ ($\partial_x\mathrm{Im}F = \mathrm{Im}F'$、$\partial_y\mathrm{Im}F = \mathrm{Re}F'$)、力は $\Delta p_x = -b_y$、$\Delta p_y = +b_x$ なので $\Delta(p_x - ip_y) = -F'$。$n = m - 1$ で `tmulk` と一致する。

**(2) 回転の向き**。コード (Fortran・pysad とも同じ) は、入口で粒子の座標と運動量に

$$
x' = x\cos\theta - y\sin\theta,\quad y' = x\sin\theta + y\cos\theta
\quad\Longleftrightarrow\quad w' = e^{i\theta}w,\qquad p_{w}' = e^{i\theta}p_w
$$

($p_w = p_x + ip_y$ も同じ回転) を施し、$\theta = \dfrac2{\mathrm{nord}}\arg(\mathrm{ak}_1 + i\,\mathrm{ak}_2) = \dfrac{\arg\kappa_c}m$、強さ $|\kappa_c|$ の**正規**成分で `ttfrin` ((7.17)) を作用させ、$-\theta$ で戻す。回転した座標での作用は $\Delta w' = -\dfrac{|\kappa_c|}{4\,m!\,p}\,w'\Bigl(\dfrac{w'^m}{m + 1} - \bar w'^m\Bigr)$。戻すと ($w' = e^{i\theta}w$、$\Delta w = e^{-i\theta}\Delta w'$、$|\kappa_c|\,e^{im\theta} = \kappa_c$)

$$
\Delta w = -\frac{1}{4\,m!\,p}\;w\Bigl(\frac{\kappa_c\,w^m}{m + 1} - \bar\kappa_c\,\bar w^m\Bigr)
\tag{6.13b}
$$

(途中: $e^{-i\theta}\cdot e^{i\theta}w\cdot|\kappa_c|\bigl(e^{im\theta}w^m/(m + 1) - e^{-im\theta}\bar w^m\bigr) = w\bigl(\kappa_cw^m/(m + 1) - \bar\kappa_c\bar w^m\bigr)$。実の $\kappa_c$ なら (7.17)。) (6.13b) は、§4.1 の表の $a_w$ が強さについて実数線形であることを使い、$\kappa_c$ の実部 (正規) と虚部 (skew) に分けて足したものにも当たる (sympy で $m = 2, 3, 4$ について $b_s = \kappa'\,\mathrm{Im}(\kappa_cw^m)/m!$ を確認)。数値でも、pysad の `_ttfrins` ($m = 3$、ランダムな複素 `ak`、$\delta \neq 0$) の位置が (6.13b) と $4\times10^{-19}$ で一致する。**回転を逆向き ($w' = e^{-i\theta}w$) にすると $\bar\kappa_c$ になり**、同じ入力との差が $10^{-10}$ 出る。したがって (6.13a) の約束 ($\mathrm{ak}_2$ が $+\mathrm{Re}\,w^m$) のもとで、**向きは上の $+\theta$ に決まる**。残る未確認は、「SAD の skew 成分の符号が (6.13a) の約束と同じ向きであること」で、根拠は本体の `tmulk` との整合 (上の (ii)) だけである (SAD 本体 (gs) での直接の確認はしていない)。

- `ttfrin` の分岐: nord = 2 (双極子、$n = 0$): $x \mathrel{+}= \tfrac12\,\mathrm{aki}\,y^2$、$p_y \mathrel{-}= \mathrm{aki}\,p_xy$、$z \mathrel{-}= \tfrac12\,\mathrm{aki}\,p_xy^2$ ($\mathrm{aki} = \mathrm{ak}/(\mathrm{al}(1 + \delta))$)。これは (5.35) の $p_z \to p$、$p_x/p \to p_x$ の近似 ($\Delta x = y^2/(2\rho)$、$\Delta p_y = -y\,p_x/\rho$、$\Delta z = -y^2p_x/(2\rho)$) と同じ ((5.35) で $p_z = p_v \to p$、$p^2/p_z^3 \to 1/p$)。双極子の端は (6.1) の直線の端を仮定していて、(7.17) を $m = 1$ に当てはめた式とは違う (§3.2)。nord = 4 は四極 (§6)、nord $\ge 6$ は (7.17)。
- `tblfri` (`tmulti.f:681-720`): 複素の `ck0` (水平の曲がりが実部、垂直が虚部) と `fb1` から、水平・垂直の両方の線形フリンジを入れる。水平 (`ck0.real`、$\rho_b = \mathrm{al}/\mathrm{ck0.real}$) は (5.42) そのもの: `dxfrx` $= \mathrm{fb1}^2/\rho_b/24$ (= dxfr)、`dyfrx` $= \mathrm{fb1}/\rho_b^2/6$ (= dyfr)、`dyfrax` $= 4\,\mathrm{dyfrx}/\mathrm{fb1}^2$ (= dyfra)。$x \mathrel{+}= \mathrm{dxfrx}\,\delta/p$、$p_y \mathrel{+}= (\mathrm{dyfrx} - \mathrm{dyfrax}\,y^2)\,y/p^2$ ($q$ 単位)。

**垂直への拡張**は、水平の式を $+90°$ 回したものになる。座標の回転 $x' = -y$、$y' = x$、$p_x' = -p_y$、$p_y' = p_x$ ($w' = iw$。垂直の双極子 $\mathrm{Im}(i\kappa w) = \kappa\,x$ ($\mathrm{ak}_2$、$b_x = \kappa$) は、水平の双極子 $\mathrm{Im}(\kappa w') = \kappa\,y'$ を回したもの) を (5.42) の $f_{\rm lin}$ に代入すると

$$f_V = f_{\rm lin}(x' = -y,\ y' = x,\ p_x' = -p_y) = p_y\,\mathrm{dxfr}\Bigl(1 - \frac1p\Bigr) + \frac1p\Bigl(\frac{\mathrm{dyfr}}2x^2 - \frac{\mathrm{dyfra}}4x^4\Bigr)$$

で、(6.0a) から $\Delta y = -\partial f_V/\partial p_y = -\mathrm{dxfr}\,\delta/p$、$\Delta p_x = +\partial f_V/\partial x = (\mathrm{dyfr} - \mathrm{dyfra}\,x^2)\,x/p$、$\Delta z$ は $p_y$ の項が $-\mathrm{dxfr}\,p_y/p^2$、$x$ の項が $(\tfrac12\mathrm{dyfr} - \tfrac14\mathrm{dyfra}\,x^2)x^2/p^2$。これは、コードの `dyfry`・`dxfry`・`dxfray` (垂直の $\rho_b = \mathrm{al}/\mathrm{ck0.imag}$ から水平と同じ式で作る) による $y \mathrel{-}= \mathrm{dyfry}\,\delta/p$、$p_x \mathrel{+}= (\mathrm{dxfry} - \mathrm{dxfray}\,x^2)\,x/p^2$、$z \mathrel{+}= \bigl(-\mathrm{dyfry}\,p_y + (\tfrac12\mathrm{dxfry} - \tfrac14\mathrm{dxfray}\,x^2)x^2/p\bigr)/p$ と同じ ($q$ 単位)。数値でも、pysad の `_tblfri` の垂直の枝 (`ck0` が純虚) は、「座標・運動量を回して水平の枝を作用させ、戻す」と差 0 (edge_06 (e))。水平と垂直が両方あるときは、コードが 2 つの写像を 1 つの式の中で足し合わせる ($x$ を更新してから $p_x$ の更新に使う順序つき。$\mathrm{dyfry}\cdot\mathrm{dxfrx}$ などの交差項は 2 次の微小量で、導出の範囲外)。以上は 05 章 §6.3 の (5.42) の導出 (水平) を前提にした拡張で、Maxwell の場の数値積分は**垂直についてはしていない**。
- 以上は K64 のコードを読んだ範囲で、**MULT の端部を数値では確認していない**。ANGLE $\ne 0$ の MULT は [08_mult.md](08_mult.md)。

### 9.3 ソレノイド

(6.10) の $f$ の流れが `tsfrin` (`tsol.f:885-976`)。Fortran のコメントの定義: $b = \tfrac12\,\mathrm{dbz}/p$、$bq = -b/4$、$bp = b/2$、$p_r = x p_x + y p_y$、$p_\phi = xp_y - yp_x$、$r_1 = r\,e^{bq\,p_\phi}$、$\phi_1 = \phi + bq\,p_r$、$z_1 = z + bp\,p_rp_\phi$ (コードの `bq = -0.125*dbz/(1 + g)`、`bp = -2*bq`)。

**$f$ の流れの導出**: 座標 $(u, \phi)$、$u = \ln r$ と、共役な運動量 $(p_u, p_\phi)$、$p_u = xp_x + yp_y\ (= p_r)$、$p_\phi = xp_y - yp_x$ に正準変換する (ポアソン括弧 $[u, p_u] = [\phi, p_\phi] = 1$、他は 0。sympy で確認、edge_05 (i))。(6.10) の $f = k\,p_\phi p_u$ ($k = b_0/(8p^2)$) は $u$、$\phi$ を含まないので、$p_u$、$p_\phi$ は保存され、(6.0a) から

$$
\Delta u = -\frac{\partial f}{\partial p_u} = -k\,p_\phi,\qquad \Delta\phi = -\frac{\partial f}{\partial p_\phi} = -k\,p_u,\qquad \Delta z = -\frac{\partial f}{\partial p} = +\frac{2k}{p}\,p_\phi p_u = +\frac{b_0}{4p^3}p_\phi p_u
$$

($\Delta z$ は $p_\phi p_u$ が一定なので 1 ステップ。) **$(x, y, p_x, p_y)$ に戻す**。$u = \ln r$ なので $r_1 = r\,e^{\Delta u}$、角は $\phi_1 = \phi + \Delta\phi$:

$$
r_1 = r\,e^{-k\,p_\phi},\quad \phi_1 = \phi - k\,p_u,\qquad x_1 = r_1\cos\phi_1,\quad y_1 = r_1\sin\phi_1
$$

運動量は、$p_u$ と $p_\phi$ が保存される ($f$ が $u$、$\phi$ を含まない) ことと、定義 $p_u = x p_x + y p_y$、$p_\phi = x p_y - y p_x$ の逆から。$x$ 倍・$y$ 倍して足し引きすると: $x\,p_u - y\,p_\phi = x^2p_x + xy\,p_y - xy\,p_y + y^2p_x = r^2p_x$、$y\,p_u + x\,p_\phi = xy\,p_x + y^2p_y + x^2p_y - xy\,p_x = r^2p_y$。よって新しい位置 $(x_1, y_1)$、$r_1^2 = x_1^2 + y_1^2$ と、保存された $p_u$、$p_\phi$ から

$$
p_{x,1} = \frac{x_1\,p_u - y_1\,p_\phi}{r_1^2},\qquad p_{y,1} = \frac{y_1\,p_u + x_1\,p_\phi}{r_1^2},\qquad z_1 = z + \frac{2k}{p}\,p_\phi p_u
$$

(sympy: 定義から恒等式。数値: $b = 1.7$ /m、$\delta = 0.01$ で `tsfrin` と $1.7\times10^{-18}$ で一致。edge_06 (f)。$q$ 単位では $k \to b_0/(8p)$ で、$p_u$・$p_\phi$ も $q$ 単位。) コードの $bq\,p_\phi$ は、$q$ 単位の運動量 ($p_\phi/p$) を使うので $bq = -b_0/(8p)$ と $k\,p_\phi = (b_0/8p^2)p_\phi = (b_0/8p)(p_\phi/p)$ が対応し、$bp = -2bq$ が $\Delta z$ の $2k/p$ に対応する (向きは §2 の約束で一致)。

- **厳密な流れとの一致**: $f$ の $\exp(:f:)$ の流れ ((6.0a) のハミルトン方程式) を DOP853 で積分した流れと `tsfrin` は、$(x, p_x, y, p_y, z)$ で $4 \times 10^{-19}$ 以内 ($b = 1.7$ /m、3 mm・20 mrad) ([edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (a))。**$f$ は paraxial ($p_z \to p$) なので、閉じた式は $f$ の厳密な流れで、シンプレクティック性は構成から厳密**。
- **Maxwell の端の場との比較**: (6.3) の場 ($g = (1 + \tanh(s/d))/2$) のローレンツ力の積分 (kinematic momentum、窓 $\pm8d$)。モデルは「自由な運動 → 硬い端 (**正準運動量が連続**) → 一様な $B_z$ の中」で、これと ODE の差から**線形な端の集束 ($b$ の偶数次) を除くため、$b \to -b$ で符号が変わる部分 (奇数部分) を取り**、$d = 1, 0.5, 0.25$ mm から $d \to 0$ に外挿して `tsfrin` の効果と比べた ($b = 1$ /m、$x_0 = 0.3$ mm、$y_0 = -0.2$ mm、$p_x = 3$ mrad、$p_y = 2$ mrad。効果は $10^{-10}$ 程度): $x$: $0.9997$、$p_x$: $0.97$、$y$: $1.06$、$p_y$: $0.83$ ([edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) (b))。効果が小さく ($10^{-11}\sim10^{-10}$) 外挿の取り方に依るので、**$x$ 以外の成分は数 % 〜 17% の確認**にとどまる。具体的には、$d = 0.5/1$ mm の組から外挿した値 (外挿 1) と $d = 0.25/0.5$ mm の組から外挿した値 (外挿 2) の、`tsfrin` に対する比は

| 成分 | $x$ | $p_x$ | $y$ | $p_y$ |
|---|---|---|---|---|
| 外挿 1 ÷ `tsfrin` | 0.9988 | 1.016 | 0.972 | 1.110 |
| 外挿 2 ÷ `tsfrin` (上の $0.9997$、$0.97$、$1.06$、$0.83$) | 0.9997 | 0.974 | 1.064 | 0.835 |
| $\lvert$外挿 1 − 外挿 2$\rvert$ ÷ $\lvert$`tsfrin`$\rvert$ | 0.09% | 4.2% | 9.3% | **27.5%** |

最後の行が「$d$ への残りの依存」の定義で、2 つの外挿の差を `tsfrin` の効果の大きさで割った値 ($p_y$ の $1.1405\times10^{-10}$ と $8.58\times10^{-11}$ の差 $2.8\times10^{-11}$ を `tsfrin` の $1.03\times10^{-10}$ で割って約 27%。edge_06 (g))。$d$ の外挿が収束していないので、$p_y$ の一致の度合いは 17% (外挿 2 のずれ) よりよい保証はなく、27% の幅を持つ。$x$ は 0.1% 以内で確認でき、$p_x$・$y$ は数 %〜10% 以内、$p_y$ は 3 割程度の幅で `tsfrin` と符号・大きさが合う、という意味である。
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
| 線形フリンジ F1・F2 (追跡) | `tquad.f:501-528` | `tqlfri` (K64 では未使用。`tquad.f:53-67` 等にインライン) | (4.33) |
| 偏向の端 (CALC) | `tbfrie.f:88-125` | `tbedge` | §8.1 |
| 偏向の Maxwell (更新と行列) | `tbfrie.f:190-229` | `tbedgemaxwell` | (5.36)。05 章 付録 B-1 (COD の更新の因子 $p$ は K64 の誤り、gs で確認) |
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
| 偏向の Maxwell (COD・行列) | `lattice._tbedgemaxwell_cod`、`_tbedgemaxwell_mat`、`elkern._tbedgemaxwell` | `tbedgemaxwell` | k64 は K64 の直訳 (05 章 付録 B-1 の誤りを残す)、v2 は COD の更新を直した |
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
| `edge_02_generating_function.py` | 母関数 $G$ ↔ $x + V + \tfrac12 V\cdot\nabla V$ と $p_1 = J^\top p_2$ (sympy)、厳密な $\exp(:f:)$ との差 ($K^3$、$K$ の走査、$\tfrac16(V\cdot\nabla)^2V$ との比)、$z$ の更新 ($\partial G/\partial p$) | [edge_02_generating_function.out.txt](outputs/edge_02_generating_function.out.txt) |
| `edge_03_bend_maxwell_recheck.py` | ヘルプ・追跡の同値、偏向の Maxwell の因子 $p$ (ローレンツ力の積分)、CALC の行列と更新、1 つの端面のずれの大きさ | [edge_03_bend_maxwell_recheck.out.txt](outputs/edge_03_bend_maxwell_recheck.out.txt) |
| `edge_04_solenoid_fringe.py` | `tsfrin` = $f$ の流れ、Maxwell の端の場のローレンツ力の積分との比較 ($b$ の奇数部分) | [edge_04_solenoid_fringe.out.txt](outputs/edge_04_solenoid_fringe.out.txt) |
| `edge_06_derivations.py` | (a) (6.7) の一般の導出 (相互作用表示。2 自由度で厳密な 1 次と一致、括弧の保存)、(b) $(V\cdot\nabla)V_y$、(c) K64 の四極の写像 = $\exp(:f':)$、$U = V - \tfrac16\mathcal D^2V$、(d) CAVI: $\phi$・$H_f$・$\Delta p_x$・$\Delta E$・$z$ の更新・静的な場の ODE ($1/\beta$)、(e) MULT: 複素強度・`ttfrins` の回転の向き・`tblfri` の垂直、(f) ソレノイドの端の戻し式、(g) edge_04 (b) の外挿の差 | [edge_06_derivations.out.txt](outputs/edge_06_derivations.out.txt) |
| `edge_05_derivation_steps.py` | §3〜§4 の途中式: (6.1) の漸化式、(6.2) の $c$、(6.3) の発散・回転、CAVI の $E_r$ と $\partial H/\partial\phi = e/\beta_z$、(6.6) の 1 次・2 次、(6.7) の検算、偏向の力の積分、$a$・$b$ と三角関数の展開、ソレノイドの $(u, \phi)$ 座標 | [edge_05_derivation_steps.out.txt](outputs/edge_05_derivation_steps.out.txt) |

## 14. 未確認・限界・使った条件

- **使った条件**: 偏向の ODE は $\rho_b = 50$ m、$y = 0.3$ mm、端の形 $\tanh$、$p_y = 0$、$q = 0.1, 0.3$、$\delta = 0, 0.3, -0.2$。エネルギーは規格化 ($p_0 = 1$) で、絶対値を使わない。四極の流れの比較は $K = 250 \sim 1.6 \times 10^4$ /m²、入射 4 mm・3 mrad など (大きな強さで 3 次の差を出した)。ソレノイドは $|b| = 1$ /m、振幅 0.3 mm。
- **SAD (gs) との比較**: 偏向の Maxwell の因子 $p$ (05 章 付録 B-1) と、CAVI の端の $1/\beta$ ([10 章](10_cavity.md) §7.1。B330) だけ。偏向は gs の COD で確認した (4 次元、$\delta = 0.05$〜$0.2$・$-0.1$、$y \le 8$ mm、3 GeV/c、$|\rho_b| = 5$ m。$\Delta z$ の gs での確認はしていない)。ほかの式は、ソース・ヘルプ・独立な数値計算だけで、gs とは比べていない。
- **確認できなかった**: (i) 偏向の $a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$ ((5.36) の $y^2/\rho_b^2$ 項): 1 次の摂動 (6.8) では出ない ($\rho_b^{-3}$ の 3 次の量。2 次以上の摂動は階段の極限で定義できない)。$\sin\theta/\theta$、$(1 - \cos\theta)/\theta^2$ の展開に一致することは確かめたが、導出していない (足りないもの: 端の形を指定した $\epsilon^3$ の摂動計算または厳密な軌道、$\rho_b \to -\rho_b$ の奇数部分の数値分離、$y \gg d$ の磁極の形に依る領域の扱い。具体的な内訳は §8.4)。(ii) 偏向の $p_y$ 依存の補正 $f_4$ (§4.3): 導出と sympy の展開のみで、ローレンツ力の積分では未確認 ($y \ll d$ の窓では $p_y$ を大きくして $y$ を動かせない)。K64 には無い。(iii) CAVI の端: (6.13) は導出した (§9.1。$E(t)$ の符号と大きさの読み、層の中で $E$ が変わらないこと、$\phi$ の 1 次を仮定)。$1/\beta$ の因子は導出では入り、ヘルプ・`tcavfrin` には無い (K64 の誤り。ローレンツ力の数値積分と gs で確認、v2 は修正 = B330)。ヘルプの規格化、空洞全体の ODE との比較は未確認。`tcavfrin` の $q$ の再スケールの欠落 (非シンプレクティック) は B332 (gs で確認、v2 は修正)。(iv) MULT の端部: `ttfrins` の回転の向きは座標の回転から決めた (§9.2。SAD の skew の符号の約束の直接確認は未)、`tblfri` の垂直は水平の回転として導出した (数値の Maxwell 積分は垂直について未)。MULT 全体の端部を gs や磁場の数値積分と比べてはいない。(v) ソレノイドは $x$ で 0.1%、$p_x$・$y$・$p_y$ は数 %〜17% の確認で、2 つの外挿の差が $p_y$ で 27% ある (§9.3)。(vi) 線形フリンジ ($F1$・$F2$) の導出は 04・05 章。偏向の $L'$ の補正の導出は 05 章。(vii) ヘルプの式は、K64 版同梱の原文 `Documents/SADHelp.html` (式は TeX と平文の両方で埋め込まれている) を直接読んで照合した (付録 A)。ヘルプの図 (`SADHelp_img`) は端部の式に関係しない。(viii) $K^4$ 以上の個別の係数、点変換の 2 次の項 $t$ の物理的な意味 (§6。K64 の写像が $f' = f + \tfrac16\,p\cdot\mathcal D^2V + O(K^4)$ の厳密な流れであること (6.12b) は導出した)。
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
| `fa`、`f`、`ff` | `fa` $= y/(\rho_e\sqrt{1 - q^2})$、`f` $= a\cdot$`fa`、`ff` $= \tfrac14(f + fa)\,y/(1 - q^2)$ | $\Delta p_y$ と $\Delta x$ を同じ `fa` から作る (§8.2、05 章 §6.2) |
| `pv1sqi` | $1/(1 - q_x^2)$ | `ff` の分母 |
| `dxfr1`、`dyfr1`、`dyfra1`、`dzfr1` | $F1^2/(24\rho_b)$、$F1/(6\rho_b^2)$、$4\,\mathrm{dyfr}/F1^2$、$\mathrm{dxfr}\sin\psi_1$ | 偏向の線形フリンジ (5.42)。`dyfra` は DISFRIN で 0 |
| `akk`、`aki` | `akk` $= \mathrm{ak}/\mathrm{al}/4$ ($m = 2$)。`aki` $= \mathrm{akk}/p$ | $K/(4p)$。$a = \mathrm{aki}\,x^2$、$b = \mathrm{aki}\,y^2$ |
| `a`、`b`、`d`、`ab`、`t` | $a = Kx^2/(4p)$、$b = Ky^2/(4p)$、$d = a + b$、$ab = a - b$、$t = ab^2/6$ | (4.28)〜(4.29)。$t = \tfrac12(V\cdot\nabla)V$ の係数。**注意: Fortran の `d` は $a + b$ (四極の `ttfrin` の中の局所変数) で、本章 §1・§4.4・§7・§8 の端の幅 $d$ (磁場が立ち上がる長さ) とは別の量**。§6 の式 (6.11)〜(6.12b) には端の幅は出てこない |
| `dx1`、`dy1` | $x(a/3 + b + t)$、$-y(a + b/3 - t)$ | 点変換のずれ ($t$ を含む、$W - x$ に当たる)。`z` の更新は `dx1 + t*xi` ($t$ が 2 回) |
| `xx`、`yy`、`h`、`f` | $\partial x_2/\partial x_1$、$\partial y_2/\partial y_1$、$\partial x_2/\partial y_1$、$\det J$ | (4.30)。運動量の逆転置 |
| `kord`、`an` | 多極の次数 $m$、$\mathrm{an}(k) = k$ | (7.17) |
| `bq`、`bp`、`pphi`、`phi0`、`w` | ソレノイド: $bq = -\mathrm{dbz}/(8p)$、$bp = -2\,bq$、$p_\phi$、$bq\,p_r$、$e^{bq\,p_\phi}$ | (6.10) の流れ (§9.3) |
| `vf`、`s0`、`dpt`、`dh`、`offset`、`dphis` | CAVI: $V/(2Lp_0)$、$\sin\varphi_s - \mathrm{offset}$、キックの強さ ($q$ 単位)、エネルギーの変化 | §9.1。(6.13) |
| `ck0` | MULT: 複素の $K_0$ (実部が水平、虚部が垂直) | `tblfri` |

## 付録 A. ヘルプの記述・K64 との違い (端部)

| 箇所 | 内容 |
|---|---|
| 原文との照合 (`SADHelp.html`) | 次の箇所を原文で確認した。BEND: `F1`・`FRINGE`・`DISFRIN`・`transformation:BEND` (線形・非線形フリンジの式、$F1 = 6\int(B_y/B_0 - (B_y/B_0)^2)ds$、$V$ の $y^4/(6\,\mathrm{F1}\rho_b^2)$ の項)。QUAD: `transformation:QUAD` (母関数 $G$)、`F1`・`F2` ($I_n$ の定義)。CAVI: 本文の $H_f$ (`Hf = - (e (VOLT+DVOLT)/L)(Sin(omega t - dphi) + Sin(dphi) - offset) (x^2+y^2)/4 delta(s-s0)`)。SOL: `DISFRIN` の $H$。FRINGE・DISFRIN の表 (QUAD・MULT・BEND)。この章に書いた式は原文と一致 (食い違いは下の各行) |
| ヘルプ SOL の DISFRIN | $H = -B_z/(8p^2B\rho)\,p_\phi p_r$ は paraxial ($p_z \to p$)。(6.10) の厳密な形は $p_z^2$。K64 は paraxial。この章の $f = -H_{\rm help}$ (§2) |
| ヘルプ BEND の非線形フリンジ | $p_{y2}$ の式に $p_1$ が 1 つ多く見える (05 章 付録 A。$q$ の単位なら一致する。同じ式の $x_2$ は正準の $p_1$ で書かれ、規格化が明記されていないので、誤植か規格化の省略か確定できない)。$p_y$ を $p_z$ に入れない (§4.3) |
| ヘルプ CAVI の $H_f$ | $\delta$ 関数の式のみ。導出・$1/\beta$ の因子は書かれていない |
| ヘルプの FRINGE・DISFRIN の表 | BEND の表は非線形の列だけ。BEND の DISFRIN は Maxwell を止めない (05 章 付録 A) |
| ヘルプに無いもの | (1) 端の磁場の縦成分と高次の項の導出 (§3)。(2) 母関数が 2 次までの打ち切りであること (§6)。(3) ソレノイドの端の $H$ の導出 (§4.3)。(4) CAVI の $H_f$ の導出 |

## 付録 B. K64 の誤り (確定 = B-1・B-4、疑い = B-2・B-3。取り下げは無い)

**B-1 (確定。gs で確認。05 章 付録 B-1、docs/30 の B324。v2 は修正)**: CALC の `tbedgemaxwell` の COD の更新 (`tbfrie.f:209`、`:211`) が、$\Delta x$ で $p^1$、$\Delta z$ で $p^0$。正しい (追跡・ヘルプ・ローレンツ力の積分) のは $p^2$、$p^1$。ローレンツ力の積分 (ODE の外挿値 ÷ CALC の式) が $\delta = 0.3$ で $\Delta x$、$\Delta z$ とも 1.3000、$\delta = -0.2$ で 0.8000 (§8.3)。解析行列は $p^2$ 型・$p^1$ 型の微分に一致し、COD の更新と内部で矛盾する。gs の閉軌道の比較で確認 (§8.3)。pysad の v2 は直してある。

**B-2 (05 章 B-2 と同じ)**: BEND の DISFRIN が Maxwell を止めない (ヘルプの記述との食い違い。追跡と CALC で一貫)。

**B-3 (規格化)**: ヘルプの偏向の式の $p_1$ の数が合わない (K64 のコードは $q$ 単位で整合した形。ヘルプの誤植か規格化の省略かは確定できない)。

**B-4 (確定。10 章 §7.1・§7.2。v2 は修正)**: CAVI の端 (`tcavfrin`) の収束とエネルギー変化に $1/\beta$ が無い (B330。ローレンツ力の数値積分と gs で確認)。端のエネルギー変化のあと $q = p_x/p$ を再スケールしない (B332。gs のヤコビアンで確認)。

**この章で見つけた、K64 の誤りではないもの**: (i) K64 の四極の母関数 $G$ は $\exp(:f:)$ の 2 次までの打ち切り (仕様)。(ii) 偏向の Maxwell フリンジは $p_y$ の補正を持たない (4 次。仕様)。(iii) `tsetfringep` の $F1$ は $|F1|^2$ で符号を落とす (04 章の指摘)。
