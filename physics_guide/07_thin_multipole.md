# 07. 薄い多極 (SEXT・OCT・DECA・DODECA)

**この章で分かること**: 6 極 (SEXT)、8 極 (OCT)、10 極 (DECA)、12 極 (DODECA) の、(1) ハミルトニアンと磁場の関係 (積分強度 $K_k$ の定義)、(2) 薄いキックが**厳密に解ける** (リー演算子の級数が止まる) こと、(3) 長さのある多極を K64 が 4 次で追跡する構成 (ドリフト + キック + 補正キック。[03_symplectic_integrators.md](03_symplectic_integrators.md) の力勾配型を多極に当てはめる)、(4) CALC (行列) の 2 次の分割との関係、(5) 軌道がずれたとき (フィードダウン)、(6) 端部の非線形フリンジ (Maxwell 方程式から) を、途中の式を省かずに追う。ヘルプに結果の式だけが書いてある $\exp(:F_{\rm in}:)\exp(:aL:)\cdots$ の途中を埋める。式番号 $(7.n)$ から、§11 の Fortran の対応表と §12 の pysad の対応表を引ける。

**前提**: [03_symplectic_integrators.md](03_symplectic_integrators.md) (特に §8 の力勾配型の係数 $\alpha, \beta, \gamma$)、[C_poisson_lie.md](C_poisson_lie.md) ($\exp(:f:)$ の記法。ポアソン括弧の定義、合成の順序)、[02_drift.md](02_drift.md) (厳密なドリフト)、[04_quad.md](04_quad.md) §7 (四極の端部)。四極は薄い多極の最低次 ($k = 1$) だが、本体は別の構成 (04 章)。

**読む順**: §1 (結論) → §3〜§5 (キックと分割) → §6 (4 次) → §7〜§8 (行列・フィードダウン) → §9 (フリンジ)。コードを読みながら使うときは §11・§12 から入る。

## 1. 結論

- **ポテンシャル**: 2 次元の調和関数 $V_k = \dfrac{K}{(k+1)!}\,\mathrm{Re}\,(x - iy)^{k+1}$ ($k = 1$: 四極、2: 6 極、3: 8 極、4: DECA、5: DODECA)。$K$ は中面の磁場 $b_y(x, 0) = K\,x^k/k!$ の $k$ 階微分 ($K = \partial^kb_y/\partial x^k$ [1/m$^k$]。SAD の `K2`、`K3`、`K4`、`K5` は、これに長さをかけた**積分強度** $K_kL$)。ハミルトニアンは $H = T + V$ ($T$ はドリフトの運動エネルギーの項)。
- **薄いキックは厳密に解ける**: 長さ 0 の多極の流れは、座標を変えず運動量だけを

$$
\Delta(p_x + i\,p_y) = -\frac{K}{k!}\,(x - iy)^k
\tag{7.1}
$$

 だけ変える (K64 の `px -= aki*Re(cx)`、`py -= aki*Im(cx)`、`cx = (x - iy)^k`、`aki = ak/k!`)。リー演算子 $\exp(:f:)$ ($f = -V$) の級数は 2 項で止まる (§4)。この章の $\exp(:f:)$ は [C_poisson_lie.md](C_poisson_lie.md) の約束 ($[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg)$、転送は $f = -\int H\,ds$、$\Delta x = -\partial f/\partial p_x$、$\Delta p_x = +\partial f/\partial x$) で書く。SAD のヘルプの $\exp(:V:)$ は括弧の符号が逆で、$f = -V$ に当たる (C §8.2)。シンプレクティック。pysad の追跡と $9 \times 10^{-19}$ 以下で一致 ([thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (c))。
- **長さのある多極**は、ドリフト (厳密) とキックの積に分割する。K64 の**追跡** (`tthin`) は 4 次の力勾配型: ドリフト $\alpha_1L$、キック 1/2、ドリフト $\alpha L$、**補正キック**、ドリフト $\alpha L$、キック 1/2、ドリフト $\alpha_1L$ ($\alpha = 1/\sqrt{12}$、$\alpha_1 = \tfrac12 - \alpha$。補正キックは $|\nabla V|^2$ ($\beta$) と $V_jV_kV_{jk}$ ($\gamma$) の座標だけの関数。式 (7.9)〜(7.12))。**CALC** (`tthine`) は補正なしの 2 次の分割 ($L/6$、キック 1/2、$2L/3$、キック 1/2、$L/6$)。
- **4 次であることの確認**: $L$ を 2 倍にしたときの厳密な積分との差の比が、小さい $L$ で 32 ($L^5$)。6 極 ($\kappa = 400$ /m³、$x = 20$ mm) で $L = 0.05$ m: pysad の追跡 $2.1 \times 10^{-10}$、補正キックの $\gamma$ の項なし $2.7 \times 10^{-9}$、補正なし $2.2 \times 10^{-6}$、CALC の分割 $6.3 \times 10^{-8}$、1 つのキックを挟むだけの Strang $1.7 \times 10^{-5}$ ([thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) (b))。
- **軌道がずれたとき (フィードダウン)**: 中心から $(d_x, d_y)$ ずれた 6 極のキックは、(7.1) の $x - iy$ を $(x - d_x) - i(y - d_y)$ に替えたもので、展開すると 6 極 + 四極 ($K_1 = -K\,d_x$ と、強さ $K\,d_y$ のスキュー) + 二極のキック ($-K(d_x^2 - d_y^2)/2$、$K\,d_xd_y$) の和になる (§8)。
- **端部の非線形フリンジ** ($L \neq 0$、DISFRIN = 0): $2m$ 極 ($m = k + 1$) で、硬い端の座標のずれが、K の 1 次で

$$
\Delta(x + iy) = -\frac{\kappa}{4\,m!\,p}\,w\left(\frac{w^m}{m + 1} - \bar w^m\right),\qquad w = x + iy
\tag{7.2}
$$

 ($\kappa$ は長さあたりの強さ。四極 $m = 2$ で (4.27) の $\Delta x = (K/12p)(x^3 + 3xy^2)$)。Maxwell 方程式から導き (§9)、K64 の `ttfrin` の位置のずれと一致 ($m = 3, 4, 5, 6$ で $4 \times 10^{-19}$ 以下)、ローレンツ力の数値積分 (6 極・8 極) と 0.01〜0.5% で一致した ([thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt))。運動量は点変換の逆転置、$z$ は $-\Delta\vec r\cdot\vec p$ (正準変換)。
- **K64 の癖** (v2 は直してある): CALC・Emittance[] の `tthine` は、非線形フリンジを持たない (追跡の `tthin` は持つ。軌道が 6 極の中心からずれるとずれる。付録 A)。
- **使った条件と限界**は §14 に書く。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 多極の次数 $k$ ($V \propto (x - iy)^{k+1}$。四極 1、6 極 2、…) | `SEXT`・`OCT`・`DECA`・`DODECA` | `kord = nord/2 - 1` (`tthin`、`tthine`)。`nord` = 6 (6 極)、8、10、12 | `kord`。`_THIN_KEY` (`track6.py:895`) |
| 積分強度 $K_kL$ [1/m$^{k-1}$] | `K2`・`K3`・`K4`・`K5` | `ak` | `ak = K*L` (`L ≠ 0`)、`ak = K` (`L = 0`)。要素の属性 `K2` は長さあたり |
| 長さ $L$ [m] | `L` | `al` | `L` |
| $K_k/k!$ | | `aki = ak/fact(kord)`。`fact` は `multa` の階乗の表 | `aki` |
| キックの強さ | | `akf = ak/fact(kord)`。$L \neq 0$ では $\times\tfrac12$ (入口・出口の 2 回) | `akf` |
| ずれ・回転 | `DX`、`DY`、`ROTATE` | `dx`、`dy`、`theta` | `DX`、`DY`、`ROTATION` |
| フリンジ | `DISFRIN` (0 で非線形フリンジあり) | `fringe` | `DISFRIN` |
| 全運動量 $p = 1 + \delta$、規格化した運動量 $q = p_x/p$ | | `pr`、追跡の `px(i)` | `pr`、`c[1]` |
| $x - iy$ の $k$ 乗 | | `cx` ($= $ `dcmplx(x, -y)**kord`) | `cx` |

- 追跡は規格化した運動量 $q$ で書く (01 章 §5)。多極は薄いキックでは $p$ が一定で、正準運動量 $p_x$ のキックは $\delta$ によらないが、$q = p_x/p$ のキックは $1/p$ に比例する (`aki = akf/(1+g)`)。CALC は正準運動量 (`cod(2)`) をそのまま使い、$1/(1+\delta)$ は無い。
- 添字は、Fortran は 1 始まり、pysad は 0 始まり ([02_drift.md](02_drift.md) §2)。
- 数値の例は電子 3 GeV/$c$。

## 3. ハミルトニアンとポテンシャル

直線上 ($1/\rho = 0$)、$a_x = a_y = 0$ の多極の磁場は、ベクトルポテンシャル $a_s$ (01 章の式 (1.1)) で決まる:

$$
H = -\sqrt{p^2 - p_x^2 - p_y^2} - a_s + \frac{E}{p_0v_0},\qquad
a_s = -V_k(x, y)
\tag{7.3}
$$

($\vec b = e\vec B/p_0$ として、$b_x = \partial_ya_s$、$b_y = -\partial_xa_s$。01a の §9.3〜§9.5。)

**調和性**: 電流のない場 ($\nabla\times\vec b = 0$、$\nabla\cdot\vec b = 0$) で、$a_s$ は 2 次元ラプラス方程式を満たす (2 次元で $b_s = 0$ のとき)。$(x - iy)^n$ の実部・虚部は調和関数 ($\partial_x^2 + \partial_y^2 = 4\,\partial_z\partial_{\bar z}$ で、$\bar z^n$ は $z$ に依らない)。したがって

$$
V_k = \frac{K}{(k + 1)!}\,\mathrm{Re}\,(x - iy)^{k+1}
\tag{7.4}
$$

は $\nabla^2V_k = 0$ で、$k = 1$ なら $\tfrac K2(x^2 - y^2)$ (四極の式 (4.2) の $-a_s$)、$k = 2$ なら $\tfrac K6(x^3 - 3xy^2)$。

**規格化**: 偏微分は

$$
\frac{\partial V_k}{\partial x} + i\frac{\partial V_k}{\partial y} = \frac{K}{k!}\,(x - iy)^k
\tag{7.5}
$$

($\partial_x\mathrm{Re}\,\bar z^{k+1} = (k + 1)\mathrm{Re}\,\bar z^k$、$\partial_y\mathrm{Re}\,\bar z^{k+1} = (k + 1)\mathrm{Re}(-i\bar z^k) = (k + 1)\mathrm{Im}\,\bar z^k$ から。) $a_s = -V_k$ なので $b_x = \partial_ya_s = -\partial_yV_k$、$b_y = -\partial_xa_s = +\partial_xV_k$。中面 $y = 0$ では $b_y = \partial_xV_k = K\,x^k/k!$、$b_x = -\partial_yV_k = 0$ (中面で 0)。したがって

$$
K = \left.\frac{\partial^k b_y}{\partial x^k}\right|_{x=y=0}
\tag{7.6}
$$

(sympy で $k = 1, \ldots, 5$ について: $\nabla^2V_k = 0$、$b_y(y = 0) = Kx^k/k!$、$b_x(y = 0) = 0$、$\partial_x^kb_y = K$、(7.5) の差 0。[thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (a)。) これが、SAD の $K_k = B^{(k)}L/B\rho$ ($B^{(k)}$ は $k$ 階微分。K2 は $B''L/B\rho$) の定義と一致する。

$T$ は運動エネルギー (ドリフト) の項 $-\sqrt{p^2 - p_x^2 - p_y^2} + E/(p_0v_0)$ ([02_drift.md](02_drift.md) の式 (2.1))。**ハミルトニアン全体** $H = T(p_x, p_y, p) + V_k(x, y)$ は、運動量だけの関数 $T$ と、座標だけの関数 $V_k$ の和。この形が、分割 (§5〜§6) を可能にする。

## 4. 薄いキックは厳密に解ける ($\exp(:f:)$ が 2 項で止まる)

長さ $\ell \to 0$ で、強さ $K\ell$ を保つ極限 (薄いレンズ) を考える。$T$ の寄与は消え (ドリフトの長さ $\ell \to 0$)、$H\,\ell = V_k\ell$ だけが残る。ポアソン括弧は [C_poisson_lie.md](C_poisson_lie.md) の (C.1) ($[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg)$)。長さ $\ell$ の転送は $\exp(-\ell:H:)$ なので、$f = -\ell V_k$ ($K\ell \to K_{\rm int}$、積分強度) として $\exp(:f:)$。$f$ は**座標だけの関数** ($p_x, p_y$ を含まない)。

**級数**: (C.4)(C.5) から、$:f:\,g = [f, g]$。$f$ は $x, y$ だけの関数なので、$\partial f/\partial p_x = \partial f/\partial p_y = 0$。4 つの基本変数それぞれについて、括弧 $[f, g] = \sum_{q = x, y}(\partial_qf\,\partial_{p_q}g - \partial_{p_q}f\,\partial_qg)$ の第 1 項・第 2 項を書き下す:

| $g$ | 第 1 項 $\sum\partial_qf\,\partial_{p_q}g$ | 第 2 項 $-\sum\partial_{p_q}f\,\partial_qg$ | $[f, g]$ |
|---|---|---|---|
| $x$ | $\partial_xf\cdot0 + \partial_yf\cdot0 = 0$ ($x$ は $p$ を含まない) | $-0\cdot1 = 0$ ($f$ が $p$ を含まない) | $0$ |
| $y$ | $0$ | $0$ | $0$ |
| $p_x$ | $\partial_xf\cdot\partial_{p_x}p_x = \partial_xf$ | $-0\cdot\partial_xp_x = 0$ | $\partial f/\partial x$ |
| $p_y$ | $\partial_yf\cdot\partial_{p_y}p_y = \partial_yf$ | $0$ | $\partial f/\partial y$ |

つまり

$$
:f:\,x = 0,\qquad
:f:\,y = 0,\qquad
:f:\,p_x = \frac{\partial f}{\partial x},\qquad
:f:\,p_y = \frac{\partial f}{\partial y}
\tag{7.7}
$$

2 回目: $:f:^2p_x = :f:\,(\partial_xf) = [f, \partial_xf]$。$f$ も $\partial_xf$ も座標だけの関数なので、括弧の第 1 項は $\partial_qf\cdot\partial_{p_q}(\partial_xf) = 0$、第 2 項は $\partial_{p_q}f\cdot\partial_q(\partial_xf) = 0$。よって $:f:^2p_x = 0$ (同様に $:f:^2p_y = 0$、$:f:^2x = :f:^2y = 0$)。級数 $\exp(:f:) = 1 + :f: + \tfrac12:f:^2 + \cdots$ は 2 項で切れて

$$
\exp(:f:)\,x = x,\quad \exp(:f:)\,y = y,\quad
\exp(:f:)\,p_x = p_x + \frac{\partial f}{\partial x},\quad
\exp(:f:)\,p_y = p_y + \frac{\partial f}{\partial y}
\tag{7.8}
$$

**長さ 0 の極限**: 長さ $\ell$、長さあたりの強さ $\kappa$ ($V_k$ は $\kappa$ で書く) の多極の $f = -\ell V_k$ で、$\ell \to 0$、$\kappa\ell = K_{\rm int}$ を保つと、$f = -V_k|_{K\to K_{\rm int}}$ ((7.4) の $K$ を積分強度に置き換えたもの)。以下、$K$ はこの積分強度 ($L = 0$ の `ak`) とする。

$f = -V_k$ とおいて (7.5) を使うと (7.1): $\Delta p_x = -\partial V_k/\partial x = -\dfrac K{k!}\mathrm{Re}\,(x - iy)^k$、$\Delta p_y = -\partial V_k/\partial y = -\dfrac K{k!}\mathrm{Im}\,(x - iy)^k$。具体例:

- 6 極 ($k = 2$): $f = -\tfrac K6(x^3 - 3xy^2)$。$\partial_xf = -\tfrac K6(3x^2 - 3y^2) = -\tfrac K2(x^2 - y^2)$、$\partial_yf = -\tfrac K6(-6xy) = +K\,xy$。したがって $\Delta p_x = -\tfrac K2(x^2 - y^2)$、$\Delta p_y = +K\,xy$ (複素形: $(x - iy)^2 = x^2 - y^2 - 2ixy$ の実部・虚部に $-K/2!$)。
- 8 極 ($k = 3$): $f = -\tfrac K{24}\mathrm{Re}(x - iy)^4 = -\tfrac K{24}(x^4 - 6x^2y^2 + y^4)$。$\partial_xf = -\tfrac K{24}(4x^3 - 12xy^2) = -\tfrac K6(x^3 - 3xy^2)$、$\partial_yf = -\tfrac K{24}(-12x^2y + 4y^3) = \tfrac K6(3x^2y - y^3)$。

**シンプレクティック性** (新しい変数 $x' = x$、$y' = y$、$p_x' = p_x + f_x$、$p_y' = p_y + f_y$。$f_x \equiv \partial_xf$ など): 括弧の基本値 (C.2) から $[x, p_x'] = [x, p_x] + [x, f_x] = 1 + 0 = 1$ ($[x, f_x]$ は $f_x$ が $p$ を含まないので 0)、$[y, p_y'] = 1$、$[x, p_y'] = [x, p_y] + [x, f_y] = 0$、$[y, p_x'] = 0$、そして
$$[p_x', p_y'] = [p_x + f_x,\ p_y + f_y] = [p_x, f_y] + [f_x, p_y] = -\partial_xf_y + \partial_yf_x = 0$$
($[p_x, f_y] = \partial_xp_x\,\partial_{p_x}f_y - \partial_{p_x}p_x\,\partial_xf_y = -\partial_xf_y$、$[f_x, p_y] = \partial_yf_x\,\partial_{p_y}p_y = \partial_yf_x$。$[p_x, p_y] = [f_x, f_y] = 0$ も使った。最後は $\partial_x\partial_yf = \partial_y\partial_xf$。$V_k$ が調和でなくても成り立つ)。

([thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (b): sympy で $[f, x] = [f, y] = 0$、$[f, [f, p_x]] = [f, [f, p_y]] = 0$ (2 項で止まる)、新しい変数の括弧 $[x, p_x'] = [y, p_y'] = 1$、$[x, p_y'] = [y, p_x'] = [p_x', p_y'] = 0$。)

**pysad・K64 との対応**: $L = 0$ の `tthin` (`al == 0`) は、`akf = ak/fact(kord)` (全部の強さ) のキック 1 回だけ (§6)。pysad の薄い 6 極・8 極・DECA・DODECA の追跡 (`L = 0`、$\delta = 0, 0.01$) と (7.1) の最大差は、$4.3 \times 10^{-19}$ 以下 (同 (c))。正準運動量のキックは $\delta$ によらず、追跡の内部変数 $q$ では $1/p$ が付く。

## 5. 長さのある多極の分割

長さ $L$、長さあたりの強さ $\kappa = K_{\rm int}/L$ ($K_{\rm int}$ = SAD の `K2` など) の多極は、$H = T + V$ で、$T$ (厳密なドリフト、[02_drift.md](02_drift.md) の式 (2.6)) と $V$ (座標だけの関数) の流れ ($A$ と $B$。[03_symplectic_integrators.md](03_symplectic_integrators.md) §3) を交互に掛けて近似する。$B$ の長さ $\ell$ の流れは、(7.8) で $f = -\ell V$ とした厳密なキック。2 次の分割は 1 つの $L$ で

$$
S_2 = e^{(L/2)A}\,e^{LB}\,e^{(L/2)A}\quad(\text{ドリフト } L/2\text{、キック (全体)、ドリフト } L/2)
$$

誤差は [03_symplectic_integrators.md](03_symplectic_integrators.md) の (3.2) で $-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]$ に比例する。K64 の 2 つの実装は、これより良い 2 つの分割を使う。

| 実装 | 分割 ($L$ の割合) | 次数 | 消える項 |
|---|---|---|---|
| 追跡 `tthin` (`tquad.f:128`) | $\alpha_1$、キック 1/2、$\alpha$、補正、$\alpha$、キック 1/2、$\alpha_1$ | 4 次 | $[A,[A,B]]$、$[B,[B,A]]$ (3 次)、5 次の座標だけの項 |
| CALC `tthine` (`tthine.f:1`) | $1/6$、キック 1/2、$2/3$、キック 1/2、$1/6$ | 2 次 | $[B,[B,A]]$ (3 次の力の勾配の項)。$[A,[A,B]]$ が $\tfrac1{72}$ で残る |
| (比較) Strang | 1/2、キック、1/2 | 2 次 | なし |

**表の「消える項」の導出** (7 因子の BCH。確認: [thin_06_seven_stage_bch.out.txt](outputs/thin_06_seven_stage_bch.out.txt) (B))。1 刻み (長さ $h = L$) の積

$$
P = e^{aA}\,e^{B/2}\,e^{tA}\,e^{C}\,e^{tA}\,e^{B/2}\,e^{aA},\qquad a = \tfrac12 - t,\quad C = c_3\,[B,[A,B]]
$$

($A$、$B$ はドリフトとキックの生成の演算子に刻みの長さ $h$ を含めた 1 次の量、$C$ は 3 次の量。ドリフトの合計は $2a + 2t = 1$) の対数を、**内側の因子から順に、前後対称な 3 因子の積に分けて**求める。BCH の公式は [C_poisson_lie.md](C_poisson_lie.md) の (C.14) ($\log(e^Pe^Q) = P + Q + \tfrac12[P,Q] + \tfrac1{12}[P,[P,Q]] - \tfrac1{12}[Q,[P,Q]] + \cdots$。最後の項は $+\tfrac1{12}[Q,[Q,P]]$ と同じ) で、04 章 §6.1 の (4.19) の導出と同じ 2 回の BCH を使う。4 次以上は捨てる。

**(i) 前後対称な 3 因子の BCH**。$X$、$Z$ を 1 次の量として $\log(e^{X}e^{Z}e^{X})$ を 3 次まで求める。

1. $Y_1 \equiv \log(e^Xe^Z) = X + Z + \tfrac12[X,Z] + \tfrac1{12}[X,[X,Z]] - \tfrac1{12}[Z,[X,Z]]$ (BCH を $P = X$、$Q = Z$ で 1 回目)。
2. $\log(e^{Y_1}e^{X}) = Y_1 + X + \tfrac12[Y_1,X] + \tfrac1{12}[Y_1,[Y_1,X]] - \tfrac1{12}[X,[Y_1,X]]$ (2 回目。$P = Y_1$、$Q = X$)。右辺を 3 次まで展開する。
   - $\tfrac12[Y_1,X]$: $Y_1$ の 1 次と 2 次の部分を入れて $\tfrac12[Z,X] + \tfrac14[[X,Z],X]$ ($[X,X] = 0$)。$[Z,X] = -[X,Z]$、$[[X,Z],X] = -[X,[X,Z]]$ なので $-\tfrac12[X,Z] - \tfrac14[X,[X,Z]]$。
   - $\tfrac1{12}[Y_1,[Y_1,X]]$: 3 次では $Y_1 \to X + Z$ でよい。$[X + Z,\,[Z,X]] = -[X,[X,Z]] - [Z,[X,Z]]$。
   - $-\tfrac1{12}[X,[Y_1,X]]$: $Y_1 \to Z$ で $-\tfrac1{12}[X,[Z,X]] = +\tfrac1{12}[X,[X,Z]]$。
3. 次数ごとに集める。**2 次**: $\tfrac12[X,Z]$ ($Y_1$ から) $- \tfrac12[X,Z]$ (上の第 1 項) $= 0$ (前後対称なので偶数次は消える)。**3 次** の $[X,[X,Z]]$: $\tfrac1{12}$ ($Y_1$) $- \tfrac14 - \tfrac1{12} + \tfrac1{12} = -\tfrac16$。$[Z,[X,Z]]$: $-\tfrac1{12}$ ($Y_1$) $- \tfrac1{12} = -\tfrac16$。$[Z,[X,Z]] = -[Z,[Z,X]]$ なので $+\tfrac16[Z,[Z,X]]$。

$$
\log\bigl(e^{X}e^{Z}e^{X}\bigr) = 2X + Z + \tfrac16\,[Z,[Z,X]] - \tfrac16\,[X,[X,Z]] + O(5)
\tag{7.8a}
$$

($X = A/2$、$Z = B$ とすると $-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]$ で、03 章の (3.2) と同じ。非可換記号の自由代数で (7.8a) を 3 次まで確認した: 同 (B) (i)。)

**(ii) 内側の 3 因子** $Q = \log(e^{tA}e^{C}e^{tA})$。(7.8a) で $X = tA$、$Z = C$ (3 次の量)。$[Z,[Z,X]]$ は 7 次、$[X,[X,Z]] = t^2[A,[A,C]]$ は 5 次なので捨てて $Q = 2tA + C + O(5)$。

**(iii) 次の 3 因子** $W = \log(e^{B/2}e^{Q}e^{B/2})$。(7.8a) で $X = B/2$、$Z = Q = 2tA + C$。3 次の項は $Z$ の主要部 $2tA$ だけで書ける:

- $\tfrac16[Z,[Z,X]] = \tfrac16[2tA,[2tA,\tfrac B2]] = \tfrac16\cdot4t^2\cdot\tfrac12\,[A,[A,B]] = \tfrac{t^2}3\,[A,[A,B]]$。
- $-\tfrac16[X,[X,Z]] = -\tfrac16[\tfrac B2,[\tfrac B2,2tA]] = -\tfrac16\cdot\tfrac14\cdot2t\,[B,[B,A]] = -\tfrac t{12}\,[B,[B,A]]$。

$$
W = B + 2tA + C + \tfrac{t^2}{3}\,[A,[A,B]] - \tfrac{t}{12}\,[B,[B,A]] + O(5)
\tag{7.8b}
$$

**(iv) 外側** $\log P = \log(e^{aA}e^{W}e^{aA})$。(7.8a) で $X = aA$、$Z = W$。1 次: $2aA + (B + 2tA) = A + B$ ($2a + 2t = 1$)。3 次の新しい項は $Z$ の主要部 $B + 2tA$ で書く。$[A,A] = 0$ と $[A,[B,A]] = -[A,[A,B]]$ を使う:

- $\tfrac16[Z,[Z,X]] = \tfrac a6\,[B + 2tA,\,[B + 2tA,\,A]] = \tfrac a6\,[B + 2tA,\,[B,A]] = \tfrac a6\bigl([B,[B,A]] - 2t\,[A,[A,B]]\bigr)$。
- $-\tfrac16[X,[X,Z]] = -\tfrac{a^2}6\,[A,[A,B + 2tA]] = -\tfrac{a^2}6\,[A,[A,B]]$。

(iii) の 3 次の項と $C = c_3[B,[A,B]] = -c_3[B,[B,A]]$ を足して、**係数を集める**:

- $[A,[A,B]]$: $\dfrac{t^2}3 - \dfrac{at}3 - \dfrac{a^2}6$。$a = \tfrac12 - t$ を入れると $\dfrac{t^2}3 - \dfrac t6 + \dfrac{t^2}3 - \dfrac16\Bigl(\dfrac14 - t + t^2\Bigr) = \dfrac{t^2}2 - \dfrac1{24} = \dfrac{12t^2 - 1}{24}$。
- $[B,[B,A]]$: $-\dfrac t{12} + \dfrac a6 - c_3 = -\dfrac t{12} + \dfrac1{12} - \dfrac t6 - c_3 = \dfrac1{12} - \dfrac t4 - c_3$。

$$
\log P = (A + B) + \frac{12t^2 - 1}{24}\,[A,[A,B]] + \Bigl(\frac1{12} - \frac t4 - c_3\Bigr)[B,[B,A]] + O(5)
\tag{7.8c}
$$

(2 次の項は (i) のどの適用でも 0。$t$、$c_3$ を記号のまま、自由代数で 7 因子の積の対数を 3 次まで直接計算し、(7.8c) と一致 (差 0)、2 次の係数 0 を確認した: 同 (B) (ii)〜(vi)。)

**順序の規約** (リー演算子として掛ける順序が流れの順序と逆になること。[C_poisson_lie.md](C_poisson_lie.md) §6.1)。$P$ の因子の並びは前後対称 (回文) なので、順序を逆にしても同じ積。順序の規約の違いは、非可換な積の順序の反転 $\psi(X_1X_2) = X_2X_1$ に当たり、$\psi[X,Y] = -[X,Y]$ なので、2 重の交換子 $[A,[A,B]]$、$[B,[B,A]]$、$[B,[A,B]]$ は 2 回反転して符号が変わらない。したがって (7.8c) は順序の規約によらない (同 (B) (ix): 積が $\psi$ で不変、2 重の交換子が $\psi$ で不変、$\psi[A,B] = -[A,B]$ を確認)。

- 追跡 `tthin`: $t = 1/\sqrt{12}$ とすると $12t^2 = 1$ で第 1 係数が 0。$c_3 = \tfrac1{12} - \tfrac t4 = 0.0111645\ldots$ ($= \beta/2$) で第 2 係数が 0 (§6.2 の導出)。
- CALC `tthine`: $L/6$、キック 1/2、$2L/3$、キック 1/2、$L/6$ は $a = 1/6$、中央のドリフト $2tL = 2L/3$ すなわち $t = 1/3$、$C = 0$。第 1 係数は $(12/9 - 1)/24 = 1/72$ (残る)、第 2 係数は $\tfrac1{12} - \tfrac1{12} - 0 = 0$ (消える)。
- Strang: 上の標準形 $-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]$ (03 の (3.2))。両方が残る。


## 6. 追跡の 4 次の分割 (`tthin`)

### 6.1 構成

03 章の力勾配型 (3.9) を、多極 $V_k$ ((7.4)) に当てはめる。$A$ = ドリフト、$B$ = キック。$T = (p_x^2 + p_y^2)/(2p)$ (有効質量 $p$) で、補正は座標だけの関数 (3.15):

$$
f_C = c_3\,L^3\,\frac{|\nabla V|^2}{p} - 2c_5\,L^5\,\frac{V_jV_kV_{jk}}{p^2},\qquad c_3 = \frac\beta2,\quad c_5 = \frac\gamma2
\tag{7.9}
$$

(K64 では長さ $L$ の 1 つの刻み。$V = V_k$ は長さあたりの強さ $\kappa$ で書く。) 全体の並びは、ヘルプの並び (ヘルプの $\exp(:\cdot:)$ はハミルトニアンの時間 1 の流れで、この章の $\exp(:f:)$ とは $f = -H$ の関係。C 付録 §8.2) $\exp(:F_{\rm in}:)\exp(:aL:)\exp(:H/2:)\exp(:bL:)\exp(:V:)\exp(:aL:)\exp(:H/2:)\exp(:bL:)\exp(:F_{\rm out}:)$ のうち、2 番目の $\exp(:aL:)$ と最後の $\exp(:bL:)$ を $\exp(:bL:)$、$\exp(:aL:)$ に入れ替えた前後対称の形 (原文の並びは非対称で、そのままでは 2 次にしかならない。付録 A) に当たる:

1. 入口の非線形フリンジ (§9。`fringe` が真のとき)
2. ドリフト $\alpha_1L$ (`ala = al*alpha1`、$\alpha_1 = \tfrac12 - \tfrac1{\sqrt{12}}$)
3. キック (全体の強さの半分: `akf*0.5`)
4. ドリフト $\alpha L$ (`alb = al*alpha`)
5. **補正キック** (式 (7.12)。$\Delta q$ と $\Delta z$)
6. ドリフト $\alpha L$
7. キック (半分)
8. ドリフト $\alpha_1L$
9. 出口の非線形フリンジ (ak の符号が逆)
10. 放射の光子 (`tradk(…, al, 0)`)。1 回だけ ($L$ 全体を 1 刻みとして)

($L = 0$ なら 3 のキックが全体の強さ ($\mathrm{akf} \times 1$) で、2、4〜8 は無い。)

### 6.2 係数 $t$、$c_3$、$c_5$ と補正の母関数 (7.9) の導出

03 章 §8 の結果を、この章で使う形で一歩ずつ確かめる (sympy: [thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (B)(G))。

**(1) 演算子の対応**。C の約束では、転送 $\exp(:f:)$ の $f$ は $-\int H\,ds$ なので、ドリフトは $f_A = -T$、キックは $f_B = -V$ ($T = (p_x^2 + p_y^2)/(2m)$、有効質量 $m = p$)。演算子の交換子は括弧の演算子 $[:f:, :g:] = :[f, g]:$ ([C_poisson_lie.md](C_poisson_lie.md) §6.2) なので、$A$、$B$ の入れ子の交換子は、括弧の関数で次のようになる。

- $[A, B] \leftrightarrow [-T, -V] = [T, V] = \partial_{q_j}T\,\partial_{p_j}V - \partial_{p_j}T\,\partial_{q_j}V = 0 - \dfrac{P_j}{m}V_j = -\dfrac{\vec P\cdot\nabla V}{m}$ ($V$ は $\vec P$ を含まない)。
- $W = [B, [A, B]] \leftrightarrow \Bigl[-V,\ -\dfrac{P_jV_j}{m}\Bigr] = \dfrac1m\,[V, P_jV_j] = \dfrac1m\,\partial_{q_k}V\,\partial_{p_k}(P_jV_j) = \dfrac{V_kV_k}{m} = \dfrac{|\nabla V|^2}{m}$。
- $W_5 = [W, [A, B]] \leftrightarrow \Bigl[\dfrac{|\nabla V|^2}m,\ -\dfrac{P_jV_j}m\Bigr] = \partial_{q_k}\dfrac{|\nabla V|^2}m\cdot\partial_{p_k}\Bigl(-\dfrac{P_jV_j}m\Bigr) = \dfrac{2V_jV_{jk}}m\cdot\Bigl(-\dfrac{V_k}m\Bigr) = -\dfrac{2\,V_jV_kV_{jk}}{m^2}$ (第 2 項 $\partial_pW\cdot\partial_q(\cdot)$ は $W$ が $p$ を含まないので 0)。

($V_j = \partial_jV$、$V_{jk} = \partial_j\partial_kV$、同じ添字は和。sympy で一般の関数 $V(x, y)$ について差 0。)

**順序の規約によらない** ([C_poisson_lie.md](C_poisson_lie.md) §6.1 のように、リー演算子の積の順序は流れの順序と逆)。順序の規約を替えることは、積の順序の反転 $\psi(X_1X_2) = X_2X_1$ に当たり、交換子は $\psi[X,Y] = -[X,Y]$ と符号が反転する。したがって、交換子を $n$ 個含む入れ子の項は $(-1)^n$ 倍になる。この章で使う入れ子は $[A,[A,B]]$、$[B,[B,A]]$、$W = [B,[A,B]]$ ($n = 2$) と $W_5 = [W,[A,B]]$ ($n = 4$) で、すべて偶数個なので符号は変わらない。さらに、使う積 ((7.8c) の 7 因子) は前後対称で、順序を反転しても同じ積である。よって、補正の符号と係数は順序の規約によらない (04 章 §6.1 の (4.19) の注と同じ。確認は §5 の (7.8c) の後の段落)。

**(2) 補正の母関数**。補正 $C = c_3h^3W + c_5h^5W_5$ に当たる $f_C$ ($h$ は刻みの長さ $L$。$h$ のべきは交換子の段数に応じて 3 次、5 次) は、$m \to p$ として

$$
f_C = c_3\,L^3\,\frac{|\nabla V|^2}{p} - 2c_5\,L^5\,\frac{V_jV_kV_{jk}}{p^2}
$$

で、(7.9) そのもの。

**(3) 3 次の条件 ($t$ と $c_3$)**。1 自由度 ($H = p^2/2 + V(x)$、質量 1) で、7 段を**段ごとに** $h$ の級数で書き、厳密解と引く (sympy: [thin_06_seven_stage_bch.out.txt](outputs/thin_06_seven_stage_bch.out.txt) (A))。初期値は $x = 0$、$p = 0$、$V_n \equiv V^{(n)}(0)$。$a = \tfrac12 - t$ とし、$t$、$c_3$、$c_5$ は記号のまま残す。補正 (段 4) は (7.9) で $m \to 1$ としたもの ($|\nabla V|^2 = V'^2$、$V_jV_kV_{jk} = V'^2V''$): $\Delta p = \partial_xf_C = 2c_3h^3V'V'' - 2c_5h^5(2V'V''^2 + V'^2V''')$。

**厳密解** ($x'' = -V'(x)$ の Picard 反復)。$x(0) = p(0) = 0$ で 2 回積分する: $x(s) = -\int_0^s(s - u)\,V'(x(u))\,du$ ($\int_0^s\!du\int_0^u\!dv\,g(v) = \int_0^s(s - u)\,g(u)\,du$)。$x_0 = 0$ から反復 $x_{n+1}(s) = -\int_0^s(s - u)\,V'(x_n(u))\,du$ し、$\int_0^s(s - u)\,u^n\,du = s^{n+2}/\bigl((n + 1)(n + 2)\bigr)$ を使う。

- $x_1$: $V'(0) = V_1$ から $x_1 = -\tfrac12V_1s^2$。
- $x_2$: $V'(x_1) = V_1 + V_2x_1 + \tfrac12V_3x_1^2 = V_1 - \tfrac12V_1V_2u^2 + \tfrac18V_1^2V_3u^4$。$u^2 \to s^4/12$、$u^4 \to s^6/30$ なので $x_2 = -\tfrac12V_1s^2 + \tfrac1{24}V_1V_2s^4 - \tfrac1{240}V_1^2V_3s^6$。
- $x_3$: $V'(x_2)$ の $u^4$ の係数に $V_2\cdot\tfrac1{24}V_1V_2$ が加わる ($V_2x_2$ の項)。$-\tfrac1{24}V_1V_2^2\cdot s^6/30 = -\tfrac1{720}V_1V_2^2s^6$ が $x_2$ に加わる。$x_4$ は $s^6$ まで $x_3$ と同じ (収束)。
- $p = x'$ ($s \to h$):

$$
x_{\rm ex} = -\tfrac12V_1h^2 + \tfrac1{24}V_1V_2h^4 - \Bigl(\tfrac1{240}V_1^2V_3 + \tfrac1{720}V_1V_2^2\Bigr)h^6,\quad
p_{\rm ex} = -V_1h + \tfrac16V_1V_2h^3 - \Bigl(\tfrac1{40}V_1^2V_3 + \tfrac1{120}V_1V_2^2\Bigr)h^5
\tag{7.9a}
$$

**7 段**。各段で、前の段の $(x, p)$ から次を作る ($h^5$ まで残す。$x$ の $h^6$ 以上は不要):

| 段 | 写像と、何をしたか | 結果 |
|---|---|---|
| 0 | 初期値 | $x = 0$、$p = 0$ |
| 1 | ドリフト $ah$: $x \to x + ah\,p$ | $x = 0$、$p = 0$ ($p = 0$ なので動かない) |
| 2 | キック 1/2: $p \to p - \tfrac h2V'(x)$、$V'(0) = V_1$ | $p = -\tfrac12V_1h$ |
| 3 | ドリフト $th$: $x \to x + th\,p$ | $x = -\tfrac12V_1t\,h^2$ |
| 4 | 補正: $\Delta p = 2c_3h^3V'V'' - 2c_5h^5(2V'V''^2 + V'^2V''')$ を $x = -\tfrac12V_1th^2$ で評価。$V'V'' = V_1V_2 + (V_2^2 + V_1V_3)\,x + \cdots$ | $p = -\tfrac12V_1h + 2c_3V_1V_2h^3 - h^5\bigl[c_3t\,V_1(V_2^2 + V_1V_3) + 2c_5(2V_1V_2^2 + V_1^2V_3)\bigr]$ |
| 5 | ドリフト $th$: $x \to x + th\,p$ (段 3 の $x$ に $th\cdot(-\tfrac12V_1h + 2c_3V_1V_2h^3)$ を足す) | $x = -V_1t\,h^2 + 2c_3t\,V_1V_2h^4$ |
| 6 | キック 1/2: $p \to p - \tfrac h2V'(x)$。段 5 の $x$ を $V' = V_1 + V_2x + \tfrac12V_3x^2$ に入れて $V' = V_1 - V_1V_2t\,h^2 + (2c_3t\,V_1V_2^2 + \tfrac12t^2V_1^2V_3)h^4$ | $p = -V_1h + (2c_3 + \tfrac12t)V_1V_2h^3 - h^5\bigl[(c_3t + 2c_5 + \tfrac14t^2)V_1^2V_3 + (2c_3t + 4c_5)V_1V_2^2\bigr]$ |
| 7 | ドリフト $ah$: $x \to x + ah\,p$。$t + a = \tfrac12$、$2c_3(t + a) = c_3$、$\tfrac12at = \tfrac14t - \tfrac12t^2$ を使う | $x = -\tfrac12V_1h^2 + \bigl(c_3 + \tfrac14t - \tfrac12t^2\bigr)V_1V_2h^4$、$p$ は段 6 のまま |

(段 6 の $h^5$ の係数は、段 4 の $-c_3t(V_1V_2^2 + V_1^2V_3) - 2c_5(2V_1V_2^2 + V_1^2V_3)$ から、キックの $-\tfrac h2\cdot(2c_3tV_1V_2^2 + \tfrac12t^2V_1^2V_3)h^4$ を引いて作る。)

**差** (7 段 $-$ 厳密解 (7.9a))。
- $h^0$〜$h^2$: 7 段の $x = -\tfrac12V_1h^2 + O(h^4)$、$p = -V_1h + O(h^3)$ は (7.9a) の最初の項と同じで、**差は 0** ($t$、$c_3$ に依らない)。$x$ は $h$ の偶数次だけ、$p$ は奇数次だけ。
- $h^3$ ($p$):

$$
\Delta p\big|_{h^3} = V_1V_2\Bigl(2c_3 + \frac t2 - \frac16\Bigr)h^3
\tag{7.9b}
$$

- $h^4$ ($x$): 7 段の $(c_3 + \tfrac14t - \tfrac12t^2)$ から (7.9a) の $\tfrac1{24}$ を引いて

$$
\Delta x\big|_{h^4} = V_1V_2\,\frac{24c_3 + 6t - 12t^2 - 1}{24}\,h^4
\tag{7.9c}
$$

**条件**。(7.9b) $= 0$ から $c_3 = \tfrac1{12} - \tfrac t4$。これを (7.9c) に入れると $24c_3 = 2 - 6t$ なので $\Delta x|_{h^4} = V_1V_2\,\dfrac{1 - 12t^2}{24}\,h^4$。これが 0 になる条件は $t^2 = \tfrac1{12}$、すなわち $t = \dfrac1{\sqrt{12}}$ ($t > 0$、$a = \tfrac12 - t > 0$)。**この因子 $12t^2 - 1$ は、§5 の (7.8c) の $[A,[A,B]]$ の係数と同じ**で、$p = 0$ から出発する 1 自由度の展開の $x$ の $h^4$ にも見える。$c_3 = \tfrac1{12} - \tfrac1{4\sqrt{12}} = \tfrac1{12} - \tfrac{\sqrt3}{24} = 0.0111645\ldots$ で、(7.8c) の第 2 係数 $\tfrac1{12} - \tfrac t4 - c_3$ が 0 の条件と同じ。$\beta \equiv 2c_3 = \tfrac16 - \tfrac{\sqrt3}{12} = \tfrac16 - \tfrac1{\sqrt{48}} = 0.0223291$ (K64 の `beta`) なので $c_3 = \beta/2$。(上の 7 段に $t = 1/\sqrt{12}$、$c_3 = \beta/2$ を入れると、$h^0$〜$h^4$ の $x$ と $p$ の差は 0。同 (A)。)

**(4) 5 次の条件 ($c_5$)**。$t = 1/\sqrt{12}$、$c_3 = \tfrac1{12} - \tfrac t4$ を入れて、$h^5$ の $p$ の差 (段 6 の $h^5$ の係数 $-$ (7.9a) の $h^5$ の係数) を作る。$c_3t = \tfrac t{12} - \tfrac{t^2}4 = \tfrac{\sqrt3}{72} - \tfrac1{48}$ ($t = \sqrt3/6$)、$\tfrac14t^2 = \tfrac1{48}$ を使う:

- $V_1^2V_3$: $\tfrac1{40} - (c_3t + 2c_5 + \tfrac14t^2) = \tfrac1{40} - \tfrac{\sqrt3}{72} + \tfrac1{48} - \tfrac1{48} - 2c_5 = \tfrac1{40} - \tfrac{\sqrt3}{72} - 2c_5$。
- $V_1V_2^2$: $\tfrac1{120} - (2c_3t + 4c_5) = \tfrac1{120} - \tfrac{\sqrt3}{36} + \tfrac1{24} - 4c_5 = \tfrac1{20} - \tfrac{\sqrt3}{36} - 4c_5 = 2\bigl(\tfrac1{40} - \tfrac{\sqrt3}{72} - 2c_5\bigr)$。

$$
\Delta p\big|_{h^5} = \Bigl(\frac1{40} - \frac{\sqrt3}{72} - 2c_5\Bigr)\bigl(V_1^2V_3 + 2V_1V_2^2\bigr)h^5
\tag{7.9d}
$$

で、2 つの単項式 $V_1^2V_3$、$V_1V_2^2$ が同じ因子で消える (sympy で $t$、$c_3$ を入れた $h^5$ の差が (7.9d) と一致、$c_5$ を入れると 0: 同 (A))。$c_5 = \tfrac12\bigl(\tfrac1{40} - \tfrac{\sqrt3}{72}\bigr) = \tfrac1{80} - \tfrac{\sqrt3}{144} = \tfrac{9 - 5\sqrt3}{720} = 4.7187\times10^{-4}$。$\gamma \equiv 2c_5 = \tfrac1{40} - \tfrac{\sqrt3}{72} = \tfrac1{40} - \tfrac1{24\sqrt3} = 9.4374\times10^{-4}$ (K64 の `gamma`。$\sqrt3/72 = 1/(24\sqrt3)$)。この $\gamma$ の項は、5 次の誤差のうち運動量を含まない部分 ($p = 0$ で見える部分) だけを消す。$p^2$、$p^4$ に比例する部分は、03 (3.13) のように残る。($(V_1^2V_3 + 2V_1V_2^2) = \partial_x(V_1^2V_2)$: $f_C$ の第 2 項 $-2c_5h^5V_1^2V_2$ の $x$ 微分がちょうどこの組合せ。)

### 6.3 補正キックの具体形 (7.12) の導出

$|\nabla V|^2$ と $V_jV_kV_{jk}$ を、(7.5) から作る。$V_x + iV_y = \dfrac{\kappa}{k!}\bar z^k$ ($\bar z = x - iy$)。実部・虚部に分けると $V_x = \dfrac\kappa{k!}\mathrm{Re}\,\bar z^k$、$V_y = \dfrac\kappa{k!}\mathrm{Im}\,\bar z^k$ なので

$$
|\nabla V|^2 = V_x^2 + V_y^2 = \left(\frac\kappa{k!}\right)^2\bigl|\bar z^k\bigr|^2 = \left(\frac\kappa{k!}\right)^2 r^k,\qquad r \equiv x^2 + y^2
\tag{7.10}
$$

($|\bar z^k|^2 = |\bar z|^{2k} = r^k$)。次に、$V_jV_{jk} = \tfrac12\partial_k|\nabla V|^2$ ($|\nabla V|^2 = V_jV_j$ の $k$ 微分) から $V_jV_kV_{jk} = \tfrac12 V_k\,\partial_k|\nabla V|^2 = \tfrac12\nabla V\cdot\nabla|\nabla V|^2$。(7.10) から $\nabla|\nabla V|^2 = (\kappa/k!)^2\,k\,r^{k-1}\,(2x, 2y)$ ($\partial_xr = 2x$、$\partial_yr = 2y$)。内積は $\nabla V\cdot(2x, 2y) = 2(xV_x + yV_y) = 2\dfrac\kappa{k!}(x\,\mathrm{Re}\,\bar z^k + y\,\mathrm{Im}\,\bar z^k) = 2\dfrac\kappa{k!}\mathrm{rcx1}$ ($\mathrm{rcx1} \equiv x\,\mathrm{Re}\,\mathrm{cx} + y\,\mathrm{Im}\,\mathrm{cx}$、$\mathrm{cx} = \bar z^k$。K64 の変数名)。係数 $\tfrac12\cdot2 = 1$ をかけて

$$
V_jV_kV_{jk} = k\left(\frac\kappa{k!}\right)^3 r^{k-1}\,\mathrm{rcx1}
\tag{7.11}
$$

($k = 1..5$ で sympy が (7.10)(7.11) の差を 0 と確認。thin_05 (A)。)

(K64 では、補正キックは、中央の 1 つ目のドリフト ($\alpha L$) のあとの座標 $(x_i, y_i)$ (= `xi`、`yi`) で評価し、そのあとに 2 つ目の $\alpha L$ のドリフトが来る。) $\mathrm{aki} \equiv \dfrac{\mathrm{ak}}{k!\,p}$ ($\mathrm{ak} = \kappa L$。$p$ は全運動量。K64 の補正キックの中の `aki = akf/pr*2`) とすると、(7.9) は

$$
f_C = p\,L\,c_3\,\mathrm{aki}^2\,r^k - 2c_5\,k\,p\,L^2\,\mathrm{aki}^3\,r^{k-1}\,\mathrm{rcx1}
$$

(第 1 項: $c_3L^3(\kappa/k!)^2r^k/p = c_3L\,(\kappa L/k!)^2r^k/p$ で、$\mathrm{aki}^2 = (\kappa L/k!)^2/p^2$ なので $(\kappa L/k!)^2/p = p\,\mathrm{aki}^2$ から $c_3\,p\,L\,\mathrm{aki}^2r^k$。第 2 項: $2c_5L^5k(\kappa/k!)^3 r^{k-1}\mathrm{rcx1}/p^2 = 2c_5\,k\,L^2\,[(\kappa L/k!)^3/p^2]\,r^{k-1}\mathrm{rcx1}$ で、$(\kappa L/k!)^3/p^2 = p\,\mathrm{aki}^3$ から $2c_5\,k\,p\,L^2\mathrm{aki}^3\cdots$。)

**運動量の更新**: 規格化した運動量の更新は $\Delta q_x = (1/p)\,\partial f_C/\partial x$ ($\Delta p_x = \partial f/\partial x$ (7.8)、$q = p_x/p$。$p$ は $x$ によらない定数として微分する)。必要な偏微分は 3 つ。

- $\partial_xr^k = k\,r^{k-1}\,\partial_xr = 2k\,r^{k-1}x$。$\partial_x(r^{k-1}) = 2(k - 1)r^{k-2}x$。
- $\partial_x\mathrm{rcx1}$: $\mathrm{rcx1} = x\,\mathrm{Re}\,\bar z^k + y\,\mathrm{Im}\,\bar z^k$ を $x$ で微分する。$\partial_x\bar z^k = k\bar z^{k-1}$ ($\bar z = x - iy$ なので $\partial_x\bar z = 1$)。よって $\partial_x\mathrm{rcx1} = \mathrm{Re}\,\bar z^k + k\,(x\,\mathrm{Re}\,\bar z^{k-1} + y\,\mathrm{Im}\,\bar z^{k-1})$。括弧の中は $\mathrm{Re}\bigl[(x - iy)\bar z^{k-1}\bigr] = x\,\mathrm{Re}\,\bar z^{k-1} + y\,\mathrm{Im}\,\bar z^{k-1}$ ($\mathrm{Re}[(x - iy)(u_r + iu_i)] = xu_r + yu_i$) なので $\mathrm{Re}\,\bar z^k$ に等しい。したがって $\partial_x\mathrm{rcx1} = (k + 1)\,\mathrm{Re}\,\mathrm{cx}$。
- 同様に $\partial_y\mathrm{rcx1} = \mathrm{Im}\,\bar z^k + k\,[x\,\mathrm{Im}\,\bar z^{k-1} - y\,\mathrm{Re}\,\bar z^{k-1}]$ ($\partial_y\bar z^k = -ik\bar z^{k-1}$、$\mathrm{Re}(-iu) = u_i$、$\mathrm{Im}(-iu) = -u_r$)。括弧は $\mathrm{Im}[(x - iy)(u_r + iu_i)] = xu_i - yu_r$ なので $\mathrm{Im}\,\bar z^k$ に等しく、$\partial_y\mathrm{rcx1} = (k + 1)\,\mathrm{Im}\,\mathrm{cx}$。

$f_C$ の第 1 項の $x$ 微分: $\partial_x[p\,L\,c_3\,\mathrm{aki}^2r^k] = p\,L\,c_3\,\mathrm{aki}^2\cdot2k\,r^{k-1}x$。$1/p$ をかけて $2k\,c_3L\,\mathrm{aki}^2\,r^{k-1}x$。第 2 項: $\partial_x[r^{k-1}\mathrm{rcx1}] = 2(k - 1)r^{k-2}x\,\mathrm{rcx1} + (k + 1)\,r^{k-1}\mathrm{Re}\,\mathrm{cx}$ (積の微分)。$-2c_5kpL^2\mathrm{aki}^3/p$ をかけて $-2c_5kL^2\mathrm{aki}^3\bigl[2(k - 1)r^{k-2}x\,\mathrm{rcx1} + (k + 1)r^{k-1}\mathrm{Re}\,\mathrm{cx}\bigr]$。$y$ も同じ ($x \to y$、$\mathrm{Re} \to \mathrm{Im}$)。

$\Delta z = -\partial f_C/\partial p$ ($f_C$ の $p$ 依存は、$p\,\mathrm{aki}^2 \propto 1/p$、$p\,\mathrm{aki}^3 \propto 1/p^2$ から): $\partial_p(1/p) = -1/p^2$、$\partial_p(1/p^2) = -2/p^3$ なので

$$
-\frac{\partial f_C}{\partial p} = c_3L^3\Bigl(\frac\kappa{k!}\Bigr)^2\frac{r^k}{p^2} - 4c_5L^5k\Bigl(\frac\kappa{k!}\Bigr)^3\frac{r^{k-1}\mathrm{rcx1}}{p^3} = c_3L\,\mathrm{aki}^2r^k - 4c_5kL^2\mathrm{aki}^3r^{k-1}\mathrm{rcx1}
$$

以上をまとめて

$$
\begin{aligned}
\Delta q_x &= \mathrm{aki}^2\Bigl[\,f_1\,r^{k-1}x - \mathrm{aki}\,\bigl(f_2\,r^{k-2}\,x\,\mathrm{rcx1} + f_3\,r^{k-1}\,\mathrm{Re}\,\mathrm{cx}\bigr)\Bigr],\\
\Delta q_y &= \mathrm{aki}^2\Bigl[\,f_1\,r^{k-1}y - \mathrm{aki}\,\bigl(f_2\,r^{k-2}\,y\,\mathrm{rcx1} + f_3\,r^{k-1}\,\mathrm{Im}\,\mathrm{cx}\bigr)\Bigr],\\
\Delta z &= \mathrm{aki}^2\,r^{k-1}\Bigl(f_4\,r - f_5\,\mathrm{aki}\,\mathrm{rcx1}\Bigr)
\end{aligned}
\tag{7.12}
$$

係数は

$$
f_1 = \beta kL,\quad f_2 = 2\gamma k(k - 1)L^2,\quad f_3 = \gamma k(k + 1)L^2,\quad f_4 = \tfrac12\beta L,\quad f_5 = 2\gamma kL^2
$$

(対応: $f_1 = 2kc_3L = \beta kL$。$f_2 = 2c_5k\cdot2(k - 1)L^2 = 4c_5k(k - 1)L^2 = 2\gamma k(k - 1)L^2$。$f_3 = 2c_5k(k + 1)L^2 = \gamma k(k + 1)L^2$。$f_4 = c_3L = \beta L/2$ ($c_3L\,\mathrm{aki}^2r^k = \mathrm{aki}^2r^{k-1}\cdot c_3L\,r$)。$f_5 = 4c_5kL^2 = 2\gamma kL^2$。)

**具体例**。
- 6 極 ($k = 2$): $\mathrm{cx} = (x - iy)^2 = x^2 - y^2 - 2ixy$、$\mathrm{rcx1} = x(x^2 - y^2) + y(-2xy) = x^3 - 3xy^2$、$f_1 = 2\beta L$、$f_2 = 4\gamma L^2$、$f_3 = 6\gamma L^2$、$f_4 = \tfrac12\beta L$、$f_5 = 4\gamma L^2$。$r^{k-1} = r$、$r^{k-2} = 1$:
  $\Delta q_x = \mathrm{aki}^2[2\beta L\,r\,x - \mathrm{aki}(4\gamma L^2x(x^3 - 3xy^2) + 6\gamma L^2r(x^2 - y^2))]$、$\Delta q_y = \mathrm{aki}^2[2\beta L\,r\,y - \mathrm{aki}(4\gamma L^2y(x^3 - 3xy^2) - 12\gamma L^2r\,xy)]$、$\Delta z = \mathrm{aki}^2r(\tfrac12\beta L\,r - 4\gamma L^2\mathrm{aki}(x^3 - 3xy^2))$。
- **薄い四極 ($k = 1$)**: $V = \tfrac\kappa2(x^2 - y^2)$、$\nabla V = \kappa(x, -y)$、$|\nabla V|^2 = \kappa^2r$ ((7.10) で $k = 1$)、$V_jV_kV_{jk} = V_x^2V_{xx} + V_y^2V_{yy} = \kappa^2x^2\cdot\kappa + \kappa^2y^2\cdot(-\kappa) = \kappa^3(x^2 - y^2)$ ((7.11) で $k = 1$、$\mathrm{rcx1} = x^2 - y^2$)。$f_C = c_3L^3\kappa^2r/p - 2c_5L^5\kappa^3(x^2 - y^2)/p^2$ から $\Delta q_x = (1/p)\partial_xf_C = \mathrm{aki}^2(\beta L - 2\gamma L^2\mathrm{aki})\,x$、$\Delta q_y = \mathrm{aki}^2(\beta L + 2\gamma L^2\mathrm{aki})\,y$ ($\partial_y(x^2 - y^2) = -2y$ で第 2 項の符号が $y$ で逆)、$\Delta z = \mathrm{aki}^2(\tfrac12\beta L\,r - 2\gamma L^2\mathrm{aki}(x^2 - y^2))$。これは一般式 (7.12) で $k = 1$ ($r^{k-1} = 1$、$f_2$ の項は $k - 1 = 0$ で消える、$\mathrm{Re}\,\mathrm{cx} = x$、$\mathrm{Im}\,\mathrm{cx} = -y$、$f_3 = 2\gamma L^2$、$f_5 = 2\gamma L^2$) と同じ。K64 の `kord == 1` の分岐 (`tquad.f:238-246`: `f1 = beta*al`、`f3 = gamma*2*al**2`、`f4 = f1/2`、`f5 = f3`) とも同じ式 (sympy で差 0。この分岐の中では `aki` を設定する行が見当たらない)。**この分岐に到達するか**は、K64 のソースツリー全体 (`SAD-master` の全テキストファイル) で `tthin` を検索して確かめた。`tthin` が現れるのは `tquad.f`・`tquads.f`・`tsteer.f`・`tstrad.f`・`tturn.f` の 5 ファイルだけで、サブルーチンの定義 (`tquad.f:118`) 1 つと呼び出し 5 つ。呼び出しの第 1 引数 (`nord`) は次のとおり。

| 呼び出し | `nord` | `al` |
|---|---|---|
| `tquad.f:23` (`tquad` の `al == 0` の分岐の中) | リテラル 4 | リテラル `0.d0` |
| `tquads.f:26` (`tquads` の `al .eq. 0.d0` の分岐の中) | リテラル 4 | リテラル `0.d0` |
| `tsteer.f:28`、`tstrad.f:25` | リテラル 2 ($k = 0$、偏向) | 呼び出し側の長さ |
| `tturn.f:380` (`case (icSEXT, icOCTU, icDECA, icDODECA)` の中) | `lele` = 6、8、10、12 (`tfstk.f:278-279`、`sim/MACCODE.h:9-10` の `icSEXT = 6`、`icOCTU = 8`、`icDECA = 10`、`icDODECA = 12`) | 要素の長さ |

`kord == 1` には $\mathrm{nord} = 4$ が要るが、$\mathrm{nord} = 4$ の呼び出しは 2 つだけで、どちらも `al = 0` (`tthin` の `al == 0` の分岐は 1 回のキックだけ) なので、**ソースの静的な呼び出しの範囲では、`kord == 1` の分岐 (`tquad.f:238-246`) は実行されない**。限界: 検索は `tthin` という名前の文字列に対するもので、名前を実行時に組み立てて呼ぶような間接呼び出しがあれば見つからない (ソースツリーに `tthin` の他の出現は無い)。実行で確かめたものではなく、pysad の追跡にも対応する要素が無い (§14)。

K64 の `tquad.f:209-216` ($k = 2$: `f1 = beta*2*al`、`f2 = 4*gamma*al**2`、`f3 = gamma*6*al**2`、`f4 = .5*beta*al`、`f5 = f2`)、`:261-267` (一般の $k$: `f1 = beta*kord*al`、`f2 = 2*gamma*kord*(kord-1)*al**2`、`f3 = gamma*kord*(kord+1)*al**2`、`f4 = .5*beta*al`、`f5 = 2*gamma*kord*al**2`) の 3 つの分岐 (`kord == 2`、`kord == 1`、一般) を書き写した式は、$f_C$ から作った (7.12) と、**$k = 1, 2, 3, 4, 5$ のすべてで全成分 ($\Delta q_x$、$\Delta q_y$、$\Delta z$。$\kappa$、$p$ は記号のまま) 差 0**。[thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (A)。03 章の確認 ([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (d)) は $k = 2, 3$、この確認は $k = 1$〜$5$ (DECA・DODECA を含む)。pysad の追跡は DECA ($k = 4$)・DODECA ($k = 5$) も含めて (7.12) の式と $2 \times 10^{-18}$ 以下で一致 (同 (C)。補正キックを抜くと $10^{-4}$〜$10^{-3}$ 違う)。

### 6.4 数値確認 (厳密な積分との比較)

[thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) (a)(b)。厳密な運動 $H = -\sqrt{1 - p_x^2 - p_y^2} + V_k$ ($\delta = 0$、$p = 1$) を DOP853 で積分した結果との差 (最大の成分)。初期値 $(x, p_x, y, p_y) = (20\,\mathrm{mm}, 1\,\mathrm{mrad}, -10\,\mathrm{mm}, -0.5\,\mathrm{mrad})$、$\kappa = 400$ /m³ (6 極)、$4 \times 10^4$ /m⁴ (8 極) 一定、$L$ を変える (積分強度 $\kappa L$ を $L$ とともに変える):

| 6 極、方法 | $L = 0.0125$ | $0.025$ | $0.05$ | $0.1$ | $0.2$ |
|---|---|---|---|---|---|
| pysad の追跡 (`tthin`) | $4.5 \times 10^{-13}$ | $8.6 \times 10^{-12}$ | $2.1 \times 10^{-10}$ | $6.9 \times 10^{-9}$ | $3.1 \times 10^{-7}$ |
| 補正から $\gamma$ の項を除く | $2.6 \times 10^{-12}$ | $8.3 \times 10^{-11}$ | $2.7 \times 10^{-9}$ | $8.7 \times 10^{-8}$ | $3.0 \times 10^{-6}$ |
| 補正なし ($\alpha_1$、$\alpha$ の割り振りだけ) | $3.5 \times 10^{-8}$ | $2.8 \times 10^{-7}$ | $2.2 \times 10^{-6}$ | $1.8 \times 10^{-5}$ | $1.4 \times 10^{-4}$ |
| CALC の分割 (`tthine` 型) | $5.7 \times 10^{-10}$ | $5.7 \times 10^{-9}$ | $6.3 \times 10^{-8}$ | $9.4 \times 10^{-7}$ | $2.2 \times 10^{-5}$ |
| Strang (2 次) | $2.6 \times 10^{-7}$ | $2.1 \times 10^{-6}$ | $1.7 \times 10^{-5}$ | $1.3 \times 10^{-4}$ | $1.1 \times 10^{-3}$ |

- **$\gamma$ の項を除いた 4 次の方法**は、$L$ を 2 倍にして比 32.0〜32.5 ($L^5$ の 4 次)。**補正なし**と **Strang** は比 8.0 (2 次)。
- **追跡 (`tthin`)** の比は $L$ が小さい側から 19.1、24.2、32.8、45.4 (6 極。8 極は 19.6、35.6、41.6、49.7)。32 を下回る小さい $L$ と、32 を超える大きい $L$ の両方がある。理由は 2 つ。(i) $\gamma$ の項が 5 次の座標だけの項 ($V'^2V''$ 型) を消すので、傾きの小さい粒子では 5 次の残り (運動量を含む項) が小さく、$L$ の大きい側で 7 次の項が見え始める。5 次と 7 次の誤差は逆符号で、ある $L$ で打ち消し合う。1 自由度の 6 極 ([thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt) (C)、$\kappa = 40$、$x_0 = 0.3$、$p_{x0} = 0.1$) では、$x$ の誤差が $L = 0.15$ ($-2.4\times10^{-6}$) と $0.1$ ($+7.4\times10^{-8}$) の間で符号を変え、比は $L = 0.4\to0.2$、$0.2\to0.1$、$0.1\to0.05$、$0.05\to0.025$、$0.025\to0.0125$、$0.0125\to0.00625$ で 83.7、302、5.0、22.6、27.9、30.1 と、$L \to 0$ で 32 に下から近づく。この 1 自由度の系の対称な分割で、$\gamma$ ありの傾き (対数の傾き) が 6.4、7.6、2.9 と乱れて見えるのは、この符号反転のためで、5 次より高い次数が出ているわけではない。2 自由度の 6 極 ([thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (F)、小さい傾きの初期値) でも、$x$ の誤差は $L = 0.1$ ($-2.2\times10^{-9}$) と $0.05$ ($+7.3\times10^{-12}$) の間で符号を変える。(ii) $L = 0.0125$ の誤差 ($4.5\times10^{-13}$) は厳密な積分 (DOP853、`rtol = 1e-13`) の精度に近く、$0.0125 \to 0.025$ の比 19 にはその影響が混じる可能性がある (未確認: 厳密な積分の許容を変えての確認は、DOP853 の許容の下限が $2.2\times10^{-14}$ のため、この桁では分けられない)。$L = 0.05$ で、$\gamma$ の項を入れると誤差が 13 倍小さくなる ($2.7 \times 10^{-9} \to 2.1 \times 10^{-10}$)。
- **傾きが大きい粒子** ($p_x = 30$ mrad、$p_y = -20$ mrad): $\gamma$ の項を入れても誤差はほぼ同じ ($L = 0.05$ で $9.7 \times 10^{-9}$ 対 $8.1 \times 10^{-9}$)。5 次の誤差は $p^2$・$p^4$ の項 ((3.13)。運動量を含み、座標だけの補正では消せない) が支配する。したがって、$\gamma$ の項は、**傾きの小さい粒子 (リング内の大部分の粒子) の誤差を減らす**もの。
- **CALC の分割 (`tthine` 型)**: 3 次の項 $[A,[A,B]]$ は傾きの 2 乗に比例するので、傾きが小さいと小さく、5 次の項が見えはじめる。比は $L = 0.0125 \to 0.025,\ \ldots,\ 0.1 \to 0.2$ の 4 つの値の範囲で、次のとおり (出力 thin_02 (b))。

  | | 傾きが小さい (1 mrad、0.5 mrad) | 傾きが大きい (30 mrad、20 mrad) |
  |---|---|---|
  | 6 極 ($\kappa = 400$ /m³) | 9.9、11.1、15.0、22.9 (**9.9〜22.9**) | 8.1、8.2、8.5、9.1 (**8.1〜9.1**) |
  | 8 極 ($\kappa = 4\times10^4$ /m⁴) | 9.3、10.2、12.9、28.2 (**9.3〜28.2**) | 7.9、7.8、7.0、7.8 (**7.0〜7.9**) |

  傾きが小さいときの比は、$L$ が小さい側で 8 に近く (3 次の項が残る)、$L$ が大きい側で大きくなる (5 次の項が見えてくる、と読める)。傾きが大きいときは 3 次の項が支配して 2 次 (比 8 前後)。8 極の傾きが大きい場合の 7.0 のように 8 を少し下回る値の原因は、確かめていない (未確認)。
- 8 極も同じ傾向 ($\gamma$ なしの比 31.8〜33.5、補正なしの比 8.0〜8.2)。
- **K64 の式を書き直した関数** (この script の `tthin_k64`) と pysad の追跡の最大差は $6.9 \times 10^{-18}$ (6 極)、$1.3 \times 10^{-18}$ (8 極)。

## 7. CALC (行列) の `tthine`

CALC・Emittance[] の薄い多極は、COD $\xi$ を通る粒子の写像と、そのヤコビアン ([02_drift.md](02_drift.md) §8 の `tdrife1` と同じ考え方) を、2 次の分割 (§5 の表) で作る。ドリフトは `tdrife0` (厳密。$dv$ の項を含む)。キックは正準運動量のまま ($1/p$ なし): $\Delta p_x = -\mathrm{aki}\,\mathrm{Re}\,\mathrm{cx}$、$\Delta p_y = -\mathrm{aki}\,\mathrm{Im}\,\mathrm{cx}$ (`aki = ak/fact(kord)`、半分ずつの 2 回は `aki*.5`)。キックのヤコビアンは、$\mathrm{cx} = \bar z^k$ の微分 ($\partial_x\bar z^k = k\bar z^{k-1}$、$\partial_y\bar z^k = -ik\bar z^{k-1}$) から ($g_r + ig_i \equiv k\,\bar z^{k-1}$)。$\Delta p_x + i\Delta p_y = -\mathrm{aki}\,\bar z^k$ を $x$、$y$ で微分すると、$\partial_x(\Delta p_x + i\Delta p_y) = -\mathrm{aki}\,(g_r + ig_i)$ から $\partial_x\Delta p_x = -\mathrm{aki}\,g_r$、$\partial_x\Delta p_y = -\mathrm{aki}\,g_i$。$\partial_y(\Delta p_x + i\Delta p_y) = -\mathrm{aki}\,(-i)(g_r + ig_i) = -\mathrm{aki}\,(g_i - ig_r)$ から $\partial_y\Delta p_x = -\mathrm{aki}\,g_i$、$\partial_y\Delta p_y = +\mathrm{aki}\,g_r$。行列に並べると:

$$
\frac{\partial(\Delta p_x, \Delta p_y)}{\partial(x, y)} = -\,\mathrm{aki}\begin{pmatrix} g_r & g_i\\ g_i & -g_r\end{pmatrix}
\tag{7.13}
$$

(K64 `tthine.f` の `trans1(2,1) = -aki*dble(cx1)`、`trans1(2,3) = -aki*imag(cx1)`、`trans1(4,1) = -aki*imag(cx1)`、`trans1(4,3) = aki*dble(cx1)`、`cx1 = k*(x - iy)^(k-1)`。) 行列は、これと厳密なドリフトの行列 ((2.13)) を、分割の順に掛けたもの。

**具体例 (水平面の 2×2)**。COD が中面 ($y = 0$) で $p_x = p_y = 0$ のとき、$\bar z^{k-1} = x^{k-1}$ は実数なので $g_i = 0$、$g_r = k\,x^{k-1}$ で、(7.13) は $x$ 面と $y$ 面に分かれ、$x$ 面は $\begin{pmatrix}1 & 0\\ -h & 1\end{pmatrix}$ ($h = \tfrac{\mathrm{ak}}{k!}\cdot\tfrac12\cdot k\,x^{k-1}$。半分のキックなので `aki*.5`)、$y$ 面は $h \to -h$。傾きゼロのドリフトは $\begin{pmatrix}1 & l\\ 0 & 1\end{pmatrix}$。`tthine` の並び (ドリフト $l_1 = L/6$、キック 1、ドリフト $l_2 = 2L/3$、キック 2、ドリフト $l_1$。行列は右から順に掛ける) の $x$ 面の 2×2 を、キックの位置の $h_1$、$h_2$ を使って、**右から順に 4 回の積**で作る ([thin_06_seven_stage_bch.out.txt](outputs/thin_06_seven_stage_bch.out.txt) (C)で sympy が各段を確認)。$D(l) = \begin{pmatrix}1 & l\\ 0 & 1\end{pmatrix}$、$K(h) = \begin{pmatrix}1 & 0\\ -h & 1\end{pmatrix}$ (上の $x$ 面のキック)。

1. $N_1 = K(h_1)\,D(l_1) = \begin{pmatrix}1 & l_1\\ -h_1 & 1 - h_1l_1\end{pmatrix}$ (最初のドリフトのあとに 1 つめのキック。第 2 行 $=$ $-h_1\times$ 第 1 行 $+$ 第 2 行)。
2. $N_2 = D(l_2)\,N_1 = \begin{pmatrix}1 - h_1l_2 & l_1 + l_2 - h_1l_1l_2\\ -h_1 & 1 - h_1l_1\end{pmatrix}$ (第 1 行 $\to$ 第 1 行 $+\,l_2\times$ 第 2 行)。
3. $N_3 = K(h_2)\,N_2 = \begin{pmatrix}1 - h_1l_2 & l_1 + l_2 - h_1l_1l_2\\ -h_1 - h_2 + h_1h_2l_2 & 1 - h_1l_1 - h_2(l_1 + l_2) + h_1h_2l_1l_2\end{pmatrix}$ (第 2 行 $\to$ 第 2 行 $-\,h_2\times$ 第 1 行)。
4. $M_x = D(l_1)\,N_3$ (第 1 行 $\to$ 第 1 行 $+\,l_1\times$ 第 2 行)。成分は: $M_{11} = 1 - h_1l_2 + l_1(-h_1 - h_2 + h_1h_2l_2) = 1 - h_1(l_1 + l_2) - h_2l_1 + h_1h_2l_1l_2$。$M_{12} = l_1 + l_2 - h_1l_1l_2 + l_1\bigl(1 - h_1l_1 - h_2(l_1 + l_2) + h_1h_2l_1l_2\bigr) = 2l_1 + l_2 - h_1l_1(l_1 + l_2) - h_2l_1(l_1 + l_2) + h_1h_2l_1^2l_2$。第 2 行は $N_3$ のまま。

まとめると

$$
M_x = D(l_1)\,K(h_2)\,D(l_2)\,K(h_1)\,D(l_1) = \begin{pmatrix}
1 - h_1(l_1 + l_2) - h_2l_1 + h_1h_2l_1l_2 & 2l_1 + l_2 - h_1l_1(l_1 + l_2) - h_2l_1(l_1 + l_2) + h_1h_2l_1^2l_2\\[2pt]
-h_1 - h_2 + h_1h_2l_2 & 1 - h_1l_1 - h_2(l_1 + l_2) + h_1h_2l_1l_2
\end{pmatrix}
$$

($D(l)$ はドリフト、$K(h)$ はキックの行列。$h_1$ は 1 つめのキックの位置 $x_1$ ($=$ COD の $x$) での値、$h_2$ は 2 つめのキックの位置 $x_2 = x_1 + l_2p_{x1}$ での値で、$p_{x1} = -(\mathrm{ak}/k!)\,x_1^k/2$ は 1 つめのキック後の傾き。)$h_1 = h_2 = h$ のときは、$l_1 = L/6$、$l_2 = 2L/3$ を入れる。必要な組合せは $l_1 + l_2 = \tfrac56L$、$2l_1 + l_2 = L$、$l_1l_2 = \tfrac19L^2$、$l_1(l_1 + l_2) = \tfrac5{36}L^2$、$l_1^2l_2 = \tfrac1{54}L^3$。各成分は:
- $M_{11} = 1 - h\cdot\tfrac56L - h\cdot\tfrac16L + h^2\tfrac19L^2 = 1 - Lh + \tfrac19L^2h^2$。
- $M_{12} = L - 2h\cdot\tfrac5{36}L^2 + h^2\tfrac1{54}L^3 = L - \tfrac5{18}L^2h + \tfrac1{54}L^3h^2$。
- $M_{21} = -2h + h^2\cdot\tfrac23L$。
- $M_{22} = 1 - h\cdot\tfrac16L - h\cdot\tfrac56L + h^2\tfrac19L^2 = 1 - Lh + \tfrac19L^2h^2$ ($h_2(l_1 + l_2) + h_1l_1 = \tfrac56Lh + \tfrac16Lh$)。

すなわち

$$
M_x = \begin{pmatrix}1 - Lh + \dfrac{L^2h^2}{9} & L - \dfrac{5L^2h}{18} + \dfrac{L^3h^2}{54}\\[6pt] -2h + \dfrac{2Lh^2}3 & 1 - Lh + \dfrac{L^2h^2}{9}\end{pmatrix}
$$

(行列式は 1。$h$ の 1 次まででは、全体のキックの強さ $2h$ の薄いレンズを真ん中に置いた形に近い: $M_{21} = -2h$、$M_{12} = L - 5L^2h/18$。) $k = 2, 3$、$x = 4$ mm、$L = 0.1$ m で、`elkern._jac_thin` の ($x$, $p_x$) 2×2 ブロックと、上の式 ($h_2$ に $x_2$ を使う) の最大差は $5.1\times10^{-9}$ ($k = 2$)、$9.1\times10^{-11}$ ($k = 3$) ($h_2 = h_1$ とすると $2.1\times10^{-4}$、$1.1\times10^{-5}$ と大きい)。残りは厳密なドリフトの傾き依存 ($p_{x1}^2/2 \sim 10^{-8}$) と推定する。$y$ 面のブロックも $h \to -h$ で $1.7\times10^{-9}$、$3.0\times10^{-11}$ で一致する ([thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (E))。

- **設計軌道上** ($\xi = 0$): キックのヤコビアンは 0 ($k \geq 2$ で $\bar z^{k-1} = 0$) で、行列はドリフトそのもの。追跡と CALC の差は $1.3 \times 10^{-11}$ (L = 0.1 m、4×4 と $z$ の行。[thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) (c))。
- **軌道がずれたとき**: 追跡 (4 次) と CALC (2 次) は分割が違うので、ヤコビアンが違う。厳密な積分のヤコビアンとの最大差は、6 極 ($L = 0.1$ m、$\kappa = 400$ /m³、COD $x = 20$ mm など) で、追跡 $4.2 \times 10^{-6}$、CALC $2.3 \times 10^{-4}$。追跡と CALC の差は $2.3 \times 10^{-4}$ (同 (c))。**CALC の行列は、追跡のヤコビアンに対して、COD がずれるほど 2 次の分割の誤差だけ違う**。シンプレクティック性は、CALC の行列が $2 \times 10^{-16}$ (解析的)、追跡が $3 \times 10^{-12}$ (中心差分の精度)。
- 厚い多極の磁場の勾配が小さい (リングの 6 極は $K_2L \sim 10$〜$100$ /m²、$L \sim 0.1$ m) とき、この差は実用上小さい。

## 8. フィードダウン (ずれた多極)

多極が $(d_x, d_y)$ ずれて置かれると、粒子の座標 $(x, y)$ は磁石の中心から $(x - d_x, y - d_y)$ になる (pysad の `_tent` は K64 の `inc/TENT.inc` の翻訳で、`x -= dx`、`y -= dy` のあとに $\theta$ の回転。`TENT.inc` 自体は読んでいない)。薄いキック (7.1) は、$x - iy \to (x - d_x) - i(y - d_y) \equiv \bar z - \bar d$ ($\bar d = d_x - id_y$)。二項展開 $(\bar z - \bar d)^k = \sum_{j=0}^k\binom kj\bar z^{k-j}(-\bar d)^j$ から

$$
\Delta(p_x + ip_y) = -\frac K{k!}\,(\bar z - \bar d)^k = -\frac K{k!}\sum_{j=0}^k\binom kj\,\bar z^{k-j}(-\bar d)^j
\tag{7.14}
$$

$\binom kj/k! = 1/(j!\,(k - j)!)$ を使うと、$j$ 番目の項は
$$-\frac{K_{k-j}}{(k - j)!}\,\bar z^{k-j},\qquad K_{k-j} \equiv K\,\frac{(-\bar d)^j}{j!}$$
で、(7.1) の形の**次数 $k - j$ の多極** ($j = 0$ が元の多極、$j = 1$ が 1 つ低い次数の多極、$j = 2$ がもう 1 つ低い次数…) に等しい。強さ $K_{k-j}$ は**複素数**で、実部が通常 (法線) 成分、虚部がスキュー成分の強さを表す (虚部があると、(7.1) の $\mathrm{Re}$・$\mathrm{Im}$ が入れ替わった項になる)。例 ([thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (E)。$k = 2, 3, 4$ で (7.14) の二項展開の差 0):

| $k$ | $j$ | 次数 | 強さ $K_{k-j} = K(-\bar d)^j/j!$ (実部 + $i$ 虚部) |
|---|---|---|---|
| 2 (6 極) | 1 | 四極 | $-K d_x + i\,K d_y$ |
| 2 | 2 | 二極 | $K(d_x^2 - d_y^2)/2 - i\,K d_xd_y$ |
| 3 (8 極) | 1 | 6 極 | $-K d_x + i\,K d_y$ |
| 3 | 2 | 四極 | $K(d_x^2 - d_y^2)/2 - i\,K d_xd_y$ |
| 3 | 3 | 二極 | $-K(d_x^3 - 3d_xd_y^2)/6 + i\,K(3d_x^2d_y - d_y^3)/6$ |

6 極 ($k = 2$) の場合を、$\bar z = x - iy$ と $\bar d = d_x - id_y$ で順に展開する: $(\bar z - \bar d)^2 = \bar z^2 - 2\bar d\bar z + \bar d^2$ で、$\bar z^2 = x^2 - y^2 - 2ixy$、$\bar d\bar z = (d_x - id_y)(x - iy) = (d_xx - d_yy) - i(d_xy + d_yx)$、$\bar d^2 = d_x^2 - d_y^2 - 2id_xd_y$。(7.14) の $-\tfrac K2$ をかけて、実部と虚部に分ける:

$$
\begin{aligned}
\Delta p_x &= -\tfrac K2\bigl[(x^2 - y^2) - 2(d_xx - d_yy) + (d_x^2 - d_y^2)\bigr] = -\tfrac K2(x^2 - y^2) + K(d_x\,x - d_y\,y) - \tfrac K2(d_x^2 - d_y^2),\\
\Delta p_y &= -\tfrac K2\bigl[-2xy + 2(d_xy + d_yx) - 2d_xd_y\bigr] = K\,xy - K(d_y\,x + d_x\,y) + K\,d_xd_y
\end{aligned}
\tag{7.15}
$$

(元の式: $\Delta p_x = -\tfrac K2\bigl[(x - d_x)^2 - (y - d_y)^2\bigr]$、$\Delta p_y = K(x - d_x)(y - d_y)$ を直接展開しても同じ。)

- 第 1 項: 元の 6 極。
- 第 2 項: **四極 (フィードダウン)**。$d_y = 0$ なら $\Delta p_x = +Kd_x\,x$、$\Delta p_y = -Kd_x\,y$ で、四極 ($\Delta p_x = -K_1x$、$\Delta p_y = +K_1y$) の $K_1 = -K\,d_x$。$d_x$ が大きいほど強い四極として働く。$d_y \neq 0$ の項は、スキュー四極 ($\Delta p_x \propto y$、$\Delta p_y \propto x$) で、強さの大きさは $K\,d_y$。
- 第 3 項: **二極のキック** ($x, y$ によらない)。$-K(d_x^2 - d_y^2)/2$、$K\,d_xd_y$。

(sympy で (7.15) の展開との差 0。pysad の SEXT ($K_2 = 40$、DX = 1.5 mm、DY = −0.8 mm、薄い) の追跡と、$(x - \mathrm{DX}, y - \mathrm{DY})$ を使う閉じた式の最大差 $8.7 \times 10^{-19}$。軌道が $x = y = 0$ のときの二極のキックは $\Delta p_x = -3.220 \times 10^{-5}$、$\Delta p_y = -4.800 \times 10^{-5}$ で、(7.15) の第 3 項と一致 (数値: $-\tfrac{40}2(1.5^2 - 0.8^2)\times10^{-6} = -3.22\times10^{-5}$、$40\times1.5\times(-0.8)\times10^{-6} = -4.8\times10^{-5}$)。[thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (d)。)

**意味**: 取り付け誤差の座標変換 (ずらす→キック→戻す) と、分割した写像の合成だけで、フィードダウンは**自動的に**出る。ハミルトニアンに低次の項を足す必要はない ([11_misalignment_and_coord.md](11_misalignment_and_coord.md) で一般論)。長さのある多極でも同じ (4 次の分割の中の座標がずれている)。

## 9. 端部の非線形フリンジ

### 9.1 Maxwell 方程式から

長さ $L \neq 0$ の多極の端 ($s = 0$) で、強さが $\kappa(s)$ ($0 \to \kappa_0$) と変わるとき、磁場は Maxwell 方程式から縦成分と高次の項が決まる ([04_quad.md](04_quad.md) §7.2 の四極の場合の一般化)。$2m$ 極 ($m = k + 1$。四極 $m = 2$、6 極 $m = 3$) の、**法線成分の多極**のスカラーポテンシャル ($\vec b = \nabla\Phi$、$w = x + iy$) の主要部は $\Phi_0 = \kappa(s)\,\mathrm{Im}\,w^m/m!$。中面の磁場は $b_y(y = 0) = \partial_y\Phi_0 = \kappa\,\mathrm{Im}(i\,m\,x^{m-1})/m! = \kappa\,x^{m-1}/(m - 1)! = \kappa\,x^k/k!$ で、(7.6) の定義と一致する。

$\kappa$ が $s$ に依ると、$\Phi_0$ は 3 次元のラプラス方程式 $\nabla_\perp^2\Phi + \partial_s^2\Phi = 0$ ($\nabla\times\vec b = 0$、$\nabla\cdot\vec b = 0$ から。$\vec b = \nabla\Phi$ とおくと $\nabla\cdot\vec b = \nabla^2\Phi = 0$) を満たさない: $\nabla_\perp^2\Phi_0 = 0$ ($\mathrm{Im}\,w^m$ は調和) だが $\partial_s^2\Phi_0 = \kappa''\,\mathrm{Im}\,w^m/m! \neq 0$。これを打ち消す次の項 $\Phi_2$ を、$\nabla_\perp^2\Phi_2 = -\kappa''\,\mathrm{Im}\,w^m/m!$ から決める。補助の関数 $G \equiv |w|^2\,\mathrm{Im}\,w^m$ を取り、$\nabla_\perp^2 = 4\partial_w\partial_{\bar w}$ で計算する: $G = \dfrac{w\bar w\,(w^m - \bar w^m)}{2i} = \dfrac{w^{m+1}\bar w - w\bar w^{m+1}}{2i}$。

$$
\partial_{\bar w}G = \frac{w^{m+1} - (m + 1)\,w\bar w^m}{2i},\qquad
\partial_w\partial_{\bar w}G = \frac{(m + 1)\,w^m - (m + 1)\,\bar w^m}{2i} = (m + 1)\,\mathrm{Im}\,w^m
$$

よって $\nabla_\perp^2G = 4(m + 1)\,\mathrm{Im}\,w^m$。したがって $\Phi_2 = -\kappa''\,G/(4(m + 1)!)$ ($\nabla_\perp^2\Phi_2 = -\kappa''\cdot4(m + 1)\mathrm{Im}\,w^m/(4(m + 1)!) = -\kappa''\mathrm{Im}\,w^m/m!$ で、$\partial_s^2\Phi_0$ を消す)。まとめて

$$
\Phi = \kappa(s)\,\frac{\mathrm{Im}\,w^m}{m!} - \kappa''(s)\,\frac{|w|^2\,\mathrm{Im}\,w^m}{4\,(m + 1)!} + \cdots
\tag{7.16}
$$

($\cdots$ は $\kappa^{(4)}$ の項で、$\nabla^2\Phi$ の残りの最低次。sympy で $m = 2, 3, 4$ の残りが $\kappa^{(4)}$ に比例することを確認: [thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (D)。$\nabla_\perp^2(|w|^2\mathrm{Im}\,w^m) = 4(m + 1)\mathrm{Im}\,w^m$ も $m = 2$〜$6$ で sympy が 0。[thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (a)。) 四極 ($m = 2$) では $\mathrm{Im}\,w^2 = 2xy$ で、$\Phi = \kappa\,xy - \kappa''\,xy(x^2 + y^2)/12$ ((4.24) に一致)。磁場は $\vec b = (\partial_x\Phi, \partial_y\Phi, \partial_s\Phi)$ で、$b_s = \kappa'\,\mathrm{Im}\,w^m/m! + \cdots$。

力 ($d\vec p/ds = \vec p\times\vec b/p_z$、$\vec p \simeq (p_x, p_y, p)$ と近軸で): 成分は $(\vec p\times\vec b)_x = p_yb_s - p_zb_y$、$(\vec p\times\vec b)_y = p_zb_x - p_xb_s$。$p_z \simeq p$ として $b_s$ の項 ($p_y/p$、$p_x/p$ 倍) を落とすと

$$F_x = -b_y = -\partial_y\Phi,\qquad F_y = +b_x = +\partial_x\Phi$$

**$b_s$ の項を落とす理由と、その大きさ** ([thin_06_seven_stage_bch.out.txt](outputs/thin_06_seven_stage_bch.out.txt) (D))。$b_s = \kappa'\,\mathrm{Im}\,w^m/m!$ による力は $F_x \ni +(p_y/p)\,b_s$、$F_y \ni -(p_x/p)\,b_s$ で、傾き $\theta = p_\perp/p$ に比例する。

1. **局所の大きさの比**。$b_s$ の力は $\sim\theta\,\kappa'\,|w|^m/m!$、$\kappa''$ の項の力 ($-\partial\Phi_2$) は $\sim\kappa''\,(m + 2)|w|^{m+1}/\bigl(4\,(m + 1)!\bigr)$ なので、比は $\sim 4\,\theta\,(\kappa'/\kappa'')/|w|$ で、端の幅 $d$ の端 ($\kappa'/\kappa'' \sim d$) では $\sim 4\theta d/|w|$。$\theta = 3\times10^{-4}$、$d = 0.5$〜$2$ mm、$|w| = 4$ mm で $1.5\times10^{-4}$〜$6\times10^{-4}$。
2. **積分したときの形の違い**。$\kappa'' = \kappa_0\delta'(s)$ の力は運動量に $\delta(s)$ 型のスパイクを作り、位置を階段状にずらす (§9.2。$d$ によらない)。$\kappa' = \kappa_0\delta(s)$ の $b_s$ の力は、運動量に**有限の飛び** ($\sim\theta\,\kappa_0|w|^m/m!$) を作るだけで、位置は端の幅の中で $\sim d$ しか動かない。硬い端の極限 ($d \to 0$) で、位置のずれ (7.17) に $b_s$ は入らない。
3. **数値の確認** (tanh の端。$m = 3, 4$、$|w| = 4, 2$ mm、$d = 0.5, 2$ mm。ローレンツ力の積分で、$b_s$ を入れた場合と $b_s = 0$ にした場合の、端のすぐ後 ($s = 10d$) の位置のずれ $\Delta(x + iy)$ を (7.17) と比べる): $b_s$ を入れた場合の (7.17) からの相対差は 0.1〜0.6%、$b_s = 0$ の場合は 0.1〜0.7% で、**$b_s$ の有無で位置のずれはほぼ変わらない** ($|w| = 4$ mm、$d = 0.5$ mm で、6 極の $|\Delta w|$ は $b_s$ あり $1.260\times10^{-8}$ m、なし $1.259\times10^{-8}$ m、(7.17) $1.260\times10^{-8}$ m。8 極は $1.232\times10^{-9}$、$1.231\times10^{-9}$、$1.233\times10^{-9}$ m。$d$ を 4 倍にしても同じ程度)。
4. **運動量には効く**。同じ積分で、線形な場との運動量のずれ $\Delta p$ ($K$ の 1 次) のうち、$b_s$ を落とすと変わる分は $|\Delta p_{\rm 全て} - \Delta p_{b_s\text{ なし}}|/|\Delta p_{\rm 全て}| = 0.77$〜$0.86$ (6 極)、$0.54$ (8 極) で、**運動量のずれの半分以上が $b_s$ 由来**である ($b_s$ の運動量の飛び $\sim\theta\,\kappa_0|w|^m/m!$ は、位置のずれから出る運動量の変化と同じ $K$ の 1 次)。したがって、力から出すのは位置のずれ (7.17) だけにし、運動量は $b_s$ の寄与を含めて、次節 §9.3 の正準変換 (点変換の逆転置) で決める。この組合せが、§9.4 の積分 (**$b_s$ を入れた場の積分**) と 0.5% 以内で一致する。

(§9.4 の表の $\Delta$(ODE) は、$b_s$ の力を入れた積分。この節の確認は、$b_s$ を落とした積分でも位置のずれが変わらないことと、運動量は変わることの分離。)

### 9.2 位置のずれ

硬い端 $\kappa = \kappa_0\Theta(s)$ では $\kappa' = \kappa_0\delta(s)$、$\kappa'' = \kappa_0\delta'(s)$。$\kappa''$ の項 $\Phi_2$ の力は、複素数 $F_x + iF_y$ で

$$
F_x + iF_y = -\partial_y\Phi_2 + i\,\partial_x\Phi_2 = \frac{\kappa''}{4(m + 1)!}\,(\partial_y - i\,\partial_x)G
$$

($\Phi_2 = -\kappa''G/(4(m + 1)!)$ を代入した。符号は $-\partial_y(-\ldots) = +\ldots$、$i\partial_x(-\ldots) = -i\ldots$)。$\partial_y - i\partial_x = -i\,(\partial_x + i\partial_y) = -2i\,\partial_{\bar w}$ ($\partial_x + i\partial_y = 2\partial_{\bar w}$)。$2\partial_{\bar w}G = \dfrac{w^{m+1} - (m + 1)w\bar w^m}{i} = -i\,w\bigl[w^m - (m + 1)\bar w^m\bigr]$ ($1/i = -i$) なので

$$
F_x + iF_y = \frac{\kappa''}{4(m + 1)!}\cdot(-i)\cdot\bigl(-i\,w\bigl[w^m - (m + 1)\bar w^m\bigr]\bigr) = -\,\frac{\kappa''}{4(m + 1)!}\,w\bigl[w^m - (m + 1)\bar w^m\bigr]
$$

($(-i)(-i) = i^2 = -1$)。力は $\kappa'' = \kappa_0\delta'(s)$ に比例する。運動方程式 $dP/ds = F_x + iF_y$ ($P = p_x + ip_y$) を、端を挟んで積分する。$\int_{-\varepsilon}^{s}\delta'(s')\,g\bigl(w(s')\bigr)\,ds' = \delta(s)\,g\bigl(w(s)\bigr) - \int\delta(s')\,\dfrac{d}{ds'}g\bigl(w(s')\bigr)\,ds'$ (部分積分) の第 1 項が $\delta(s)$ 型のスパイク、第 2 項は有限の飛び (傾きに比例する小さな運動量の飛びで、位置には効かない)。**運動量が $\delta(s)$ 型のスパイク** $P \ni \kappa_0\,g\,\delta(s)$ ($g \equiv (F_x + iF_y)/\kappa''$ の係数部分) を持つと、$dw/ds = P/p$ ($p_z \simeq p$) から位置は $\delta(s)$ の積分で**階段状に**ずれる: $\Delta w = \kappa_0\,g(w)/p$。$g = -\dfrac{1}{4(m + 1)!}w\bigl[w^m - (m + 1)\bar w^m\bigr]$ なので

$$
\Delta w = -\frac{\kappa_0}{4\,(m + 1)!\,p}\,w\bigl[w^m - (m + 1)\bar w^m\bigr] = -\frac{\kappa_0}{4\,m!\,p}\,w\left(\frac{w^m}{m + 1} - \bar w^m\right)
\tag{7.17}
$$

(2 つ目の等号: $(m + 1)! = (m + 1)\,m!$ で、括弧の外に $(m + 1)$ を出した。式 (7.2)。スパイクの中で $w$ 自体が $\kappa_0$ のオーダーで動くことは、$\kappa_0$ の 2 次の項になり、ここでは入れていない: $K$ の 1 次。sympy で、(7.16) の $\kappa''$ の項の力から出した $\Delta w$ と (7.17) の差が $m = 2, 3, 4$ で 0。2 つの形の差も 0: thin_05 (D)。) 四極 ($m = 2$) では $\Delta w = -\dfrac{\kappa_0}{8p}\,w\left(\dfrac{w^2}{3} - \bar w^2\right)$ を $w = x + iy$ で展開する: $w^3 = (x^3 - 3xy^2) + i(3x^2y - y^3)$、$w\bar w^2 = |w|^2\bar w = (x^2 + y^2)(x - iy)$。実部は $\dfrac{x^3 - 3xy^2}3 - (x^2 + y^2)x = -\dfrac23(x^3 + 3xy^2)$、虚部は $\dfrac{3x^2y - y^3}3 + (x^2 + y^2)y = \dfrac23\,y(3x^2 + y^2)$。したがって

$$\Delta x = \frac{\kappa_0}{12p}(x^3 + 3xy^2),\qquad \Delta y = -\frac{\kappa_0}{12p}\,y\,(3x^2 + y^2)$$

で、$\Delta x$ は (4.27) (sympy: thin_05 (D))。出口は $\kappa' \to -\kappa_0\delta$ で符号が逆 ($\kappa_0 \to -\kappa_0$)。

**K64 の `ttfrin` (`nord = 2m`、$m \geq 3$、`tquad.f:325-452`) との対応**: `aki = ak/al/fact(kord)/4/(1+g)` ($\mathrm{kord} = m$) は $\kappa_0/(4\,m!\,p)$、`cx = dcmplx(x, y)` は $w$、`cz1 = cx**(kord-1)`、`cz = cz1*cx` $= w^m$、`ca = -aki*cx*(cz/an(kord+1) - conjg(cz))` ($\mathrm{an}(m+1) = m + 1$) は $-\mathrm{aki}\,w(w^m/(m + 1) - \bar w^m)$、すなわち $\Delta w$ そのもの。`cx1 = cx + ca` が新しい $(x, y)$ ($m = 3$ の `ttfrin` の分岐は同じ式で、数が `4`・`3` と直接書いてある)。pysad の `ttfrin` ($m = 2$〜6) の位置のずれと (7.17) は、$m = 3, 4, 5, 6$ で $4 \times 10^{-19}$ 以下で一致する (丸め誤差。$\kappa_0$ を半分にしても差が減らないので、(7.17) の $K$ の 1 次の項以外は無い。**多極 ($m \geq 3$) の点変換は (7.17) そのもの**)。四極 ($m = 2$) だけ、母関数 (4.28) から来る $K$ の 2 次の項 ($t$) があり、差が $\kappa_0$ の 2 乗で減る (比 4.0。[thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (b))。

### 9.3 運動量と $z$ (正準変換)

位置が $w \to w' = w + \Delta w(w, \bar w)$ という**点変換**で、運動量は、点変換のヤコビアン $J = \partial(x', y')/\partial(x, y)$ の逆転置 $J^{-\top}$ で更新する (点変換の母関数 $F_2 = \vec x'(\vec x)\cdot\vec p'$ から、$\vec p = \partial F_2/\partial\vec x = J^\top\vec p'$)。

**複素数での書き方**。$A \equiv \partial_w\Delta w$、$B \equiv \partial_{\bar w}\Delta w$ とすると $dw' = (1 + A)\,dw + B\,d\bar w$。(7.17) の $\mathrm{aki} \equiv \kappa_0/(4\,m!\,p)$ ($p$ を含めた K64 の変数名) で $\Delta w = -\mathrm{aki}\,\bigl(w^{m+1}/(m + 1) - w\bar w^m\bigr)$ なので

$$
A = -\mathrm{aki}\,(w^m - \bar w^m) = -2i\,\mathrm{aki}\,\mathrm{Im}\,w^m \equiv -i\,a,\qquad
B = m\,\mathrm{aki}\,w\,\bar w^{m-1},\qquad a \equiv 2\,\mathrm{aki}\,\mathrm{Im}\,w^m
$$

($\partial_w(w^{m+1}/(m + 1)) = w^m$、$\partial_w(w\bar w^m) = \bar w^m$、$\partial_{\bar w}(w\bar w^m) = m\,w\bar w^{m-1}$、$\partial_{\bar w}w^{m+1} = 0$。) 面積要素の変化 (ヤコビアンの行列式) は、$dw' = (1 + A)dw + B\,d\bar w$ の 2 次元実線形写像の行列式 $|1 + A|^2 - |B|^2$ ($\partial_w$・$\partial_{\bar w}$ で書いた複素ヤコビアン) で

$$
\det J = |1 - ia|^2 - m^2\mathrm{aki}^2|w|^{2m} = 1 + a^2 - m^2\mathrm{aki}^2|w|^{2m}
$$

($|B|^2 = m^2\mathrm{aki}^2|w|^2|w|^{2(m-1)}$)。

**運動量**。$\vec p\cdot d\vec r = p_xdx + p_ydy = \mathrm{Re}\,[\mathrm{cp}\,dw]$ ($\mathrm{cp} \equiv p_x - ip_y$。$\mathrm{Re}[(p_x - ip_y)(dx + i\,dy)] = p_xdx + p_ydy$)。点変換は $\vec p'\cdot d\vec r' = \vec p\cdot d\vec r$ を保つので $\mathrm{Re}[\mathrm{cp}_1\,dw'] = \mathrm{Re}[\mathrm{cp}\,dw]$ ($\mathrm{cp}_1$ = 新しい $p_x' - ip_y'$)。左辺に $dw' = (1 + A)dw + B\,d\bar w$ を入れ、$\mathrm{Re}[u\,d\bar w] = \mathrm{Re}[\bar u\,dw]$ を使うと $\mathrm{Re}\bigl[\{\mathrm{cp}_1(1 + A) + \overline{\mathrm{cp}_1}\,\bar B\}\,dw\bigr]$。任意の $dw$ で等しいので

$$
\mathrm{cp} = (1 + A)\,\mathrm{cp}_1 + \bar B\,\overline{\mathrm{cp}_1}
$$

この式とその複素共役 $\overline{\mathrm{cp}} = (1 + \bar A)\,\overline{\mathrm{cp}_1} + B\,\mathrm{cp}_1$ を、$(\mathrm{cp}_1, \overline{\mathrm{cp}_1})$ の連立 1 次方程式として解く (係数行列の行列式は $|1 + A|^2 - |B|^2 = \det J \equiv d$)。クラメルの公式で

$$
\mathrm{cp}_1 = \frac{(1 + \bar A)\,\mathrm{cp} - \bar B\,\overline{\mathrm{cp}}}{d} = \frac{(1 + ia)\,\mathrm{cp} - m\,\mathrm{aki}\,w^{m-1}\,\overline{w\,\mathrm{cp}}}{d}
$$

($1 + \bar A = 1 + ia$、$\bar B\,\overline{\mathrm{cp}} = m\,\mathrm{aki}\,\bar w w^{m-1}\overline{\mathrm{cp}} = m\,\mathrm{aki}\,w^{m-1}\,\overline{w\,\mathrm{cp}}$)。これが K64 の `d`、`a`、`cp1 = ((1 + i a) cp - m aki cz1 conj(cx cp))/d` ($\mathrm{cz1} = w^{m-1}$、$\mathrm{cx} = w$) そのもの。sympy で、$m = 3, 4, 5$ について、$d - \det J$ が 0、`cp1` と $J^{-\top}\vec p$ ($J$ は $x' = x + \mathrm{Re}\,\Delta w$、$y' = y + \mathrm{Im}\,\Delta w$ の 2×2 ヤコビ行列) の差が 0 ([thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) (D))。

**$z$**。位置のずれが運動量に依らず、$\mathrm{aki} \propto 1/p$ ($p = 1 + \delta$) を通して $\delta$ にだけ依るので、点変換の母関数に $\delta$ 依存を含めて $F = (\vec x + \Delta\vec r(\vec x, \delta))\cdot\vec p' + z\,\delta'$ ($z$ と $\delta$ が正準共役の対。$\delta' = \delta$) とおく。$z' = \partial F/\partial\delta' = z + (\partial\Delta\vec r/\partial\delta)\cdot\vec p'$。$\Delta\vec r \propto \mathrm{aki} \propto 1/(1 + \delta)$ なので $\partial\Delta\vec r/\partial\delta = -\Delta\vec r/(1 + \delta)$ から

$$
z \to z - \frac{\Delta\vec r\cdot\vec p_{\rm new}}{1 + \delta} = z - \Delta\vec r\cdot\vec q_{\rm new} = z - \mathrm{Re}\bigl(\Delta w\cdot\mathrm{cp1}\bigr)
\tag{7.18}
$$

($\Delta\vec r\cdot\vec q_{\rm new} = \Delta x\,q_x' + \Delta y\,q_y' = \mathrm{Re}[\Delta w\,(q_x' - iq_y')] = \mathrm{Re}[\Delta w\cdot\mathrm{cp1}]$。[04_quad.md](04_quad.md) §7.3 の四極と同じ構造。`z = z - dble(ca*cp1)`。$q$ 変数なので $q_{\rm new} = p_{\rm new}/p$。) 点変換から作る写像なので、シンプレクティック性は構成から厳密に成り立つ。符号の確認 (thin_05 (D)): $m = 3$、$\delta = 0.005$ の点で、$(x, p_x, y, p_y, z, \delta)$ の 6×6 ヤコビアンの $|M^\top\Omega M - \Omega|$ は (7.18) の符号で $1.1\times10^{-16}$、$z$ の符号を逆にすると $3.6\times10^{-7}$ (符号が逆だとシンプレクティックでない)。

**確認** ([thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (c)): 6 極 ($\kappa_0 = 500$ /m³、$\delta = 0.005$) で、pysad の `ttfrin` の運動量は、位置の写像の数値ヤコビアンの逆転置 $J^{-\top}p_{\rm old}$ と $1.1 \times 10^{-15}$ で一致。$z$ の変化は $-\Delta\vec r\cdot\vec p_{\rm new}$ と $3 \times 10^{-23}$ で一致。6×6 ヤコビアン ($\delta \neq 0$ を含む。正準運動量) の $|M^\top JM - J|$ は $3.8 \times 10^{-13}$ (中心差分の精度)。

### 9.4 ローレンツ力の数値積分との比較

[thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (d) ([04_quad.md](04_quad.md) §7.2 の検証と同じ構成): $\kappa(s) = \kappa_0(1 + \tanh(s/d))/2$ の入口 ($-0.06 < s < 0.06$ m)、$p = 1$。磁場は (7.16) の 3 次までの場 ($b_s = \kappa'\mathrm{Im}\,w^m/m!$)。$\Delta(\mathrm{ODE})$ = ODE (この磁場) $-$ ODE (線形な場 $\kappa\,\mathrm{Im}\,w^m/m!$ だけ。非線形フリンジの効果)、$\Delta$(モデル) = (ドリフト → `ttfrin` → 一様な多極) $-$ (ドリフト → 一様な多極)。

| 次数 | $\kappa_0$ | 振幅 $r$ [mm] | $d$ [mm] | $\Delta$(ODE) | $\Delta$(モデル) | 差/$\Delta$(モデル) |
|---|---|---|---|---|---|---|
| 6 極 | 400 /m³ | 4 | 0.5〜2 | $1.157 \times 10^{-8}$ | $1.156 \times 10^{-8}$ | 0.0004〜0.0052 |
| 6 極 | 400 /m³ | 2 | 0.5〜2 | $7.22 \times 10^{-10}$ | $7.22 \times 10^{-10}$ | 0.0001〜0.0027 |
| 8 極 | $4 \times 10^4$ /m⁴ | 4 | 0.5〜2 | $1.210 \times 10^{-9}$ | $1.210 \times 10^{-9}$ | 0.0002〜0.0010 |
| 8 極 | $4 \times 10^4$ /m⁴ | 2 | 0.5〜2 | $3.78 \times 10^{-11}$ | $3.78 \times 10^{-11}$ | 0.0000〜0.0003 |

**非線形フリンジの効果は、6 極・8 極でも (7.17) とその運動量・$z$ の更新で、0.5% 以内で再現される**。振幅を半分にすると効果は $2^{m+1}$ 倍小さくなる (6 極で 16 倍: $1.157 \times 10^{-8} \to 7.22 \times 10^{-10}$)。実際の大きさ: 6 極 $\kappa_0 = 400$ /m³ ($K_2L = 40$ /m² で $L = 0.1$ m 相当) で、振幅 4 mm でも効果は $10^{-8}$ m で、ふつうのリングでは微小 (最も効くのは、軌道が 6 極の中心から数 mm ずれたとき、または大きなリングの大振幅)。

## 10. 放射

放射を入れるとき ($\mathrm{krad}$)、`tthin` は `tradk(np, ..., al, 0)` を、本体全体 ($L$) を 1 刻みとして**最後に 1 回**呼ぶ (入口の位置 `pxr0`・`pyr0`・`zr0` を保存して、全体の変位から光子の放射を計算する。放射の式は [13_radiation.md](13_radiation.md)。`tsetpcvt` で光子の出る点を $L$ の範囲に設定)。K64 の `tthine` (CALC) は `enarad = .false.` で放射しない (Emittance[] の減衰・励起に多極が入らない。付録 A)。

## 11. Fortran の対応表

行番号は参照した K64 のソースのもの。

| 式・内容 | ファイル:行 | サブルーチン | 変数・メモ |
|---|---|---|---|
| 追跡の薄い多極 (全体) | `tquad.f:118-` | `tthin` | `nord`・`al`・`ak`・`dx`・`dy`・`theta`・`krad`・`fringe`。`ak == 0` → `tdrift_free` (`:152`)。`TENT.inc`・`TEXIT.inc` |
| 係数 $\alpha, \beta, \gamma, \alpha_1$ | `tquad.f:131-135` | (parameter) | 03 章 §2 |
| 入口・出口のフリンジ | `tquad.f:168` (入口)、`:316` (出口) | `ttfrin(…, nord, ±ak, al, 0)` | `fringe` が真のとき |
| 刻みのドリフト | `tquad.f:173-182` (入口側)、`:303-` (出口側) | | `ala = al*alpha1`、`alb = al*alpha`、`akf = akf*.5` (`al ≠ 0`) |
| キック (全体・半分) | `tquad.f:184-` (1 回目)、`:293-` (2 回目) | | $k = 2$: `px -= aki*(x-y)*(x+y)`、`py += 2*aki*x*y`。一般: `cx = dcmplx(x,-y)**kord`、`px -= aki*dble(cx)`、`py -= aki*imag(cx)` |
| 補正キックと $z$ (式 (7.12)) | `tquad.f:208-` ($k = 2$)、`:238-` ($k = 1$)、`:260-` (一般) | | `f1`〜`f5`。`aki = akf/pr*2`。`rk1 = r**(kord-2)`、`rk = r*rk1`、`rcx1 = xi*dble(cx)+yi*imag(cx)` |
| 放射 | `tquad.f:319` | `tradk` | 最後に 1 回 |
| 非線形フリンジ (7.17)〜(7.18) | `tquad.f:~325-452` | `ttfrin` | $m = 3$ (`nord == 6`)、一般 (`nord ≥ 8`): `akk = ak/al/fact(kord)/4`、`aki = akk/(1+g)`、`ca`、`cx1`、`d`、`cp1`、`z -= dble(ca*cp1)` |
| CALC の薄い多極 | `tthine.f:1-95` | `tthine` | ドリフト `tdrife0` (`:37`、`:64`、`:87`)。`ala = al/6`、`alb = al/1.5`、`aki = ak/fact(kord)`。`tchge` でずれ・回転 (`:25`、`:90`) |
| CALC のキックのヤコビアン (7.13) | `tthine.f:52-56`、`:76-80` | | `cx1 = dcmplx(cod(1),-cod(3))**(kord-1)*kord`、`trans1(2,1) = -aki*dble(cx1)` など |
| `tthine` の呼び出し | `tturne.f:592`、`qtwiss.f:1729` (`qthin`) | | フリンジは呼ばれない (確認した範囲。付録 A) |

## 12. pysad の対応表

| 内容 | pysad の関数 (ファイル:行) | 対応する Fortran | 注意 |
|---|---|---|---|
| 追跡の 6 極 ($k = 2$) | `track6._tthin_sext` (`track6.py:898`) | `tthin` | `ak = K*L` ($L \neq 0$)、$L = 0$ は `ak = K`。係数は `_TT_ALPHA` など (`:889-892`) |
| 追跡の 8 極以上 ($k = 3, 4, 5$) | `track6._tthin_gen` (`:746`) | `tthin` の一般 `kord` | `f1`〜`f5` の一般式 (`:776-786`。`kord == 1` の分岐も持つ) |
| 追跡の非線形フリンジ | `track6._ttfrin6` (`:821`、$m = 3$)、`_ttfrin_gen` (`:1938`、$m \geq 4$)、`_ttfrin_n` (`:1960`) | `ttfrin` | `a`、`ca`、`cx1`、`d`、`cp1`。`z -= (ca*cp1).real` |
| 追跡 (numba) | `track6n.k_thin` (`track6n.py:643`)、`ttfrin6_1` (`:582`)、`ttfrin_gen_1` (`:603`) | 同上 | |
| CALC の COD の写像 | `lattice._tthine_cod` (`lattice.py:1755`)、`elkern._tthine_map` (`elkern.py:360`) | `tthine` | 正準運動量のキック |
| CALC のヤコビアン | `elkern._jac_thin` (`elkern.py:330`)、`_jac_kick` (`:302`)、`emit6d._jac_thin` (`emit6d.py:469`) | `tthine` | 式 (7.13) |
| CALC (Emittance[]) の全体 | `emit6d._tthine` (`emit6d.py:1778`)、`_tthine_gen` (`:1794`)、`elkern.k_thin` (`:558`) | `tturne.f:592` | v2 は非線形フリンジ (`_lat._thin_fringe_flag`、`lattice.py:1783`) と、$L \neq 0$ の放射 (`krad`) を足す。k64 は K64 のまま (入れない) |
| 共鳴駆動項 (6 極の項) | `resdrive._thin_funcs` (`resdrive.py:277`) | — | [18_nonlinear_bridge.md](18_nonlinear_bridge.md) |

## 13. 確認スクリプトと出力

スクリプトはリポジトリのルートから `py -3.12 physics_guide/scripts/<名前>.py` で動かす (数秒〜十数秒。sympy・scipy が要る)。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `thin_01_potential_kick.py` | (a) 多極のポテンシャル (調和、規格化、複素形)。(b) $\exp(:f:)$ の級数が 2 項で止まる、シンプレクティック。(c) pysad の薄い多極 ($k = 2$〜5、$\delta \neq 0$) と (7.1)。(d) フィードダウン (7.15) と pysad の DX・DY | [thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) |
| `thin_02_fourth_order.py` | (a) pysad の `tthin` と K64 の式の書き直しの一致。(b) 厳密な積分との差と $L$ を 2 倍にしたときの比 (4 次・2 次、$\gamma$ の効果、傾きの大小)。(c) ヤコビアン (追跡・CALC・厳密) | [thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) |
| `thin_03_fringe.py` | (a) 磁場の調和性と位置のずれ (7.17) (sympy)。(b) pysad の `ttfrin` ($m = 2$〜6) との比較。(c) 運動量・$z$・シンプレクティック性。(d) Maxwell の場の積分との比較 (6 極・8 極) | [thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) |
| `thin_04_help_order.py` | ヘルプの並び (非対称) と補正 $V$ の符号の確認 (1 自由度の 6 極)。(C) 対称形・$\gamma$ ありの誤差の符号と、$L$ を細かく取ったときの比 | [thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt) |
| `thin_05_derivations.py` | (A) $f_C$ から (7.12) を $k = 1$〜$5$ で導き、K64 の 3 つの分岐と比較 (記号)。(B) $W$、$W_5$ の表現。(C) pysad の追跡 ($k = 2$〜$5$) と (7.12)。(D) (7.16) のラプラス方程式、(7.17)・`cp1`・(7.18) と 6×6 のシンプレクティック性。(E) フィードダウンの二項展開、(7.13)、CALC の 2×2 行列。(F) 6 極の誤差の符号。(G) $t$、$c_3$、$c_5$ の導出 | [thin_05_derivations.out.txt](outputs/thin_05_derivations.out.txt) |
| `thin_06_seven_stage_bch.py` | (A) 1 自由度の 7 段を段ごとに $h$ の級数で書き、厳密解 (Picard 反復) との差から $t$、$c_3$、$c_5$ (§6.2 の (7.9a)〜(7.9d))。(B) 7 因子の BCH を三つ組の BCH ((7.8a)〜(7.8c))、積の順序の反転で不変なこと。(C) CALC の $M_x$ の 2×2 の積と各成分。(D) 非線形フリンジの $b_s$ の項を入れた積分と落とした積分の比較 (§9.1) | [thin_06_seven_stage_bch.out.txt](outputs/thin_06_seven_stage_bch.out.txt) |
| `integ_01_bch_force_gradient.py` (03 章) | 補正キック (7.12) と K64 の `tthin` の式の記号的な一致 ((d)) | [integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) |

## 14. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV/$c$ ($\delta = 0$ が主)、6 極・8 極、$L \le 0.4$ m、振幅 数 mm〜20 mm、傾き 1〜30 mrad。多極は理想 (法線成分だけ、ポテンシャル (7.4))。

**導出の範囲**: (7.9) の係数 $c_3$、$c_5$ と $t$ は §6.2 で、1 自由度の $x = p = 0$ から出発する 7 段の展開 ((7.9a)〜(7.9d)) から導いた。$t = 1/\sqrt{12}$ と $c_3$ は、§5 の (7.8c) (7 因子の BCH を、この章の記号で 3 次まで導出) とも一致する。5 次の $p^2$・$p^4$ の項は、この章では導出していない (03 章 (3.13) の自由代数の計算に依る。$\gamma$ の項はそれを消さない、と §6.4 の数値で確認した範囲)。BCH は 3 次まで ((7.8c))、5 次の BCH は 03 章。(7.16) は $\kappa''$ までの 3 次元ラプラス方程式の導出で、$\kappa''''$ 以上の項は入れていない。(7.17) の運動量のスパイクの積分は、スパイクの中の $w$ の動き ($K$ の 2 次) を落とす近似で、この近似の大きさ ($K$ の 2 次の項) は多極 $m \geq 3$ では見積もっていない (四極は 04 章 §7)。

**この章の範囲外 (別の章)**: 多極 MULT (`tmulti`・`tmulte`。スキュー・複数の成分・ANGLE) は MULT の章 ([08_mult.md](08_mult.md))。端部の一般論 (磁場の立ち上がり F1・F2、BEND のフリンジ) は端部の章 ([06_edge_fringe.md](06_edge_fringe.md))。四極の本体は [04_quad.md](04_quad.md)。放射・スピンは別の章。

**未確認・限界**:
- **補正キックの式 (7.12) の確認範囲**: 記号計算で、$f_C$ から作った式が (7.12) と、$k = 1$〜$5$ で全成分差 0。K64 の `tthin` の 3 つの分岐 (`kord == 2`、`kord == 1`、一般) をソースから書き写した式とも $k = 1$〜$5$ で差 0。pysad の追跡は、SEXT・OCT・DECA・DODECA ($k = 2$〜$5$) で (7.12) の式と $2 \times 10^{-18}$ 以下 (thin_05 (A)(C))。未確認: (a) gs (K64 の実行) の DECA・DODECA の出力との数値比較 (ソースの式との一致までで、実行結果とは比べていない)。(b) 厳密な積分との比較 (§6.4) は 6 極・8 極だけ (DECA・DODECA は、(7.12) が $k$ の一般の式であること以外の確認が無い)。(c) 薄い四極 ($k = 1$、`kord == 1`) の分岐は、(7.12) と同じ式だが、K64 では `aki` が設定されない上に、ソースツリー全体で `tthin` の呼び出しを調べた範囲では到達しないコード (§6.3。間接呼び出しは検索できない)。pysad の追跡には対応する要素が無く、数値の確認はしていない。
- **補正キックの $\gamma$ の項は全体の誤差を必ず減らすとは限らない** (03 章 §8.3)。6 極の例では、傾きの小さい粒子で誤差を 13 倍減らし、傾きが大きい粒子ではほぼ効かない (§6.4)。一般の 1 自由度のポテンシャル ($x^4/4$) では、$\gamma$ の項を入れない方が誤差が小さかった (03 章 §9)。
- **K64 のコードとの比較は、pysad の Python 版と、K64 のソースの式を書き写した関数を通して**行った。gs (K64 の実行) の出力との数値比較は、この章ではしていない。
- **非線形フリンジ**: K の 1 次 (硬い端の近似)。端の磁場の展開は 3 次まで ($b_s = \kappa'\mathrm{Im}\,w^m/m!$ など。$\kappa''''$ の項は入れていない)。Maxwell の場の数値積分との比較は、$\kappa_0$、$r$、$d$ の限られた範囲 (表)。端の幅 $d$ が孔径に近い場合や、$K$ の 2 次の項 (四極で見えた残り。[04_quad.md](04_quad.md) §7.2) は確認していない。$m \geq 3$ の運動量・$z$ の更新を Maxwell の場から直接導いた (7.17 の位置のずれのように) 検証は、$\Delta(\mathrm{ODE})$ との一致 (表) による。
- **フリンジの線形部分** (多極の F1 に当たる、端の立ち上がりの形による効果。四極の §7.4) は、この章では扱わない (SEXT・OCT … は F1・F2 を持たない。端の立ち上がりの形で決まる $K$ の 1 次の項は、設計軌道上では $\propto x^{k}$ の非線形で、確認していない)。
- **フィードダウン**の式 (7.14)(7.15) は、薄いキックで確認した。長さのある多極での ROTATE・DROTATE・DZ との組合せ、`tchge` の癖 ([11_misalignment_and_coord.md](11_misalignment_and_coord.md)) は確認していない。
- **DISFRIN**: 追跡の `fringe` は `DISFRIN == 0` のとき。`tthin` の `fringe` の既定の値 (K64 の SEXT の `DISFRIN` の既定が 0 か) は、pysad の `_g(el, 'DISFRIN') == 0.0` で使うことしか確認していない。

## 15. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。03 章の表 (`alpha` など) は 03 章 §14。この表の定義は、実際のコード (Fortran: `SAD-master/src`、pysad: `pysad/`) で確認した。

| Fortran の名前 | 式 | 定義の場所 (Fortran) | pysad での名前 |
|---|---|---|---|
| `nord` | $2(k + 1)$ (6 極で 6、8 極で 8)。`tthin`・`tthine`・`ttfrin` の引数 | `tquad.f:118`、`tthine.f:1` | (呼び出し側の `kord`) |
| `kord` | `tthin`・`tthine`: $k = $ `nord/2 - 1`。**`ttfrin`**: $m = $ `nord/2` ($= k + 1$)。同じ名前で意味が違う | `tquad.f:170`・`:429` | `kord` (`tthin` 側)、`kord`・`nord` (`ttfrin`) |
| `fact(n)` | $n!$ (`multa` の表) | `multa` | `math.factorial` |
| `ak` | 積分強度 $K_kL$ | | `ak` |
| `akf` | $K/k!$ ($= $ `ak/fact(kord)`、`tquad.f:171`)。$L \neq 0$ では $\times\tfrac12$ (入口・出口の 2 回) | `tquad.f:171`、`:182` | `akf` |
| `aki` | キックの強さ。1 回のキック: `akf/pr` ($q$ 変数。`pr` $= p$)。補正キックの中: `akf/pr*2` $= K/(k!\,p)$。`ttfrin` の中: `akk/pr`、`akk = ak/al/fact(kord)/4` | `tquad.f` (キック・補正・`ttfrin`) | `aki` |
| `cx` | $\bar z^k = (x - iy)^k$ (キック・補正)。**`ttfrin` の中の `cx` は $w = x + iy$**。同じ名前で意味が違う | `tquad.f` (`tthin`・`ttfrin`) | `cx` |
| `r`、`rk1`、`rk` | $r = x^2 + y^2$、$r^{k-2}$、$r^{k-1}$ ($= r\cdot r^{k-2}$) | `tquad.f:260-` (一般 $k$) | `r`、`rk1`、`rk` |
| `rcx1` | $x\,\mathrm{Re}\,\mathrm{cx} + y\,\mathrm{Im}\,\mathrm{cx}$ ($= \nabla V\cdot\vec r\,k!/\kappa$。(7.11) の元) | `tquad.f:208-` ほか | `rcx1` |
| `f1`〜`f5` | (7.12) の係数。$f_1 = \beta kL$、$f_2 = 2\gamma k(k-1)L^2$、$f_3 = \gamma k(k+1)L^2$、$f_4 = \beta L/2$、$f_5 = 2\gamma kL^2$ | `tquad.f:209-216` など | `f1`〜`f5` |
| `ca` | $\Delta w$ (7.17)。`-aki*cx*(cz/an(kord+1) - conjg(cz))` | `tquad.f` (`ttfrin`) | `ca` |
| `cp1` | 新しい $p_x - ip_y$ (逆転置) | `ttfrin` | `cp1` |
| `d` | $\det J = 1 + a^2 - (m\,\mathrm{aki})^2|w|^{2m}$ (`ttfrin` の中。四極の `ttfrin` の `d` $= a + b$ とは別) | `ttfrin` | `d` |
| `ala`、`alb` | $\alpha_1L$、$\alpha L$ (`tthin`)、$L/6$、$L/1.5$ (`tthine`) | `tquad.f:173-174`、`tthine.f:36-37` | `ala`、`alb` |

## 付録 A. ヘルプの記述・K64 との違い (薄い多極)

**ヘルプの記述 (疑いを含む)**

| 箇所 | 内容 |
|---|---|
| ヘルプ SEXT・OCT・DECA・DODECA の transformation (`transformation:THIN`) の並び | **原文で確認** (4 要素とも同じ式。TeX 版・平文版とも)。$\exp(:F_{\rm in}:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:V_n:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:F_{\rm out}:)$、$a = \tfrac12 - \tfrac1{\sqrt{12}}$、$b = \tfrac12 - a$ ($L$ は長さ $L$ のドリフトのハミルトニアン、$H_n$ は薄いキック)。**確定した食い違い (誤記)**: 並びが前後対称でない ($a, H/2, b, V, a, H/2, b$)。K64 の追跡 (`tquad.f:173-182`)・03 章の力勾配型は $a, H/2, b, V, b, H/2, a$。1 自由度の 6 極 ($v = \kappa x^3/6$、$\kappa = 40$、$x_0 = 0.3$、$p_x0 = 0.1$) で厳密な流れと比べた誤差は、対称形が $L = 0.4 \to 0.05$ で $1.9 \times 10^{-3} \to 1.5 \times 10^{-8}$ (γ なしで 5 次)、原文の並びは $1.4 \times 10^{-2} \to 1.8 \times 10^{-4}$ (2 次。倍率の傾き 2.0) ([thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt))。2 番目の $\exp(:aL:)$ と最後の $\exp(:bL:)$ の入れ替えの誤記。 |
| ヘルプの SEXT・OCT・DECA・DODECA の説明の文面 | **原文で確認**。TeX 版は「thin sextupole / octupole / decapole / dodecapole kick」と正しいが、**平文版 (`EquAlign` の 2 つめの文字列) は「thin dipole kick」(SEXT)、「thin 3-pole kick」(OCT)、「thin quadrupole kick」(DECA)、「thin 5-pole kick」(DODECA)** と、極数でなく添字 $n$ を数えた誤記 (確定。式 $H_n = K_n/(n+1)!\,\mathrm{Re}(x - iy)^{n+1}$ は両方とも正しく、本章 (7.1) と同じ)。6 極のキックは二極ではなく (7.1) の $k = 2$、10 極は四極ではなく $k = 4$。 |
| ヘルプの補正 $V_n$ | **原文で確認**: $V_n = \sum_{j,k=(x,y)}\bigl[-\tfrac\beta2(H_{n,k})^2 + \gamma H_{n,j}H_{n,k}H_{n,jk}\bigr]$、$\beta = 1/6 - 1/\sqrt{48}$、$\gamma = 1/40 - 1/(24\sqrt3)$ (数値は K64 の `beta`・`gamma` と同じ。03 章 $c_3 = \beta/2$、$c_5 = \gamma/2$)。$\exp(:\cdot:)$ は「ハミルトニアンの流れ」の読み (C 付録 §8.2) で、$V_n = -f_C$ (本章 (7.9)。03 章は流れの生成元を $-H$ の側に書く) と**符号も係数も一致**する: 対称形で $V_n$ の符号を逆にすると 3 次にしかならず、正しい符号で 5 次以上 ([thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt))。**食い違い (略記の疑い)**: $H_n$ は積分強度なので、補正には長さ $L$ と $1/p$ の因子 ($-\tfrac\beta2 L(H_{n,k})^2/p$ など。(7.9) の $L^3|\nabla V|^2/p$、$L^5V_jV_kV_{jk}/p^2$ を積分強度で書いたもの) が要るが、原文の式には無い。 |
| 追跡 (4 次) の分岐の網羅 | K64 の `tthin` の補正キックは `kord == 2`・`kord == 1`・一般 ($k \geq 3$) の 3 つの分岐。いずれも (7.12) と同じ式 (§6.3、記号計算)。DECA・DODECA は一般の分岐を通り、$k = 4, 5$ で (7.12) と同じ。 |
| ヘルプに無いもの | (1) 補正キックの具体形 (7.12)。(2) 追跡 (4 次) と CALC (2 次) の分割の違い。(3) 非線形フリンジの Maxwell 方程式からの導出と、運動量・$z$ の更新。(4) フィードダウンの式。 |

**K64 の実装の癖・誤りの疑い**

| 箇所 | 内容 |
|---|---|
| CALC・Emittance[] の `tthine` にフリンジがない (B291) | 追跡の `tthin` は非線形フリンジ `ttfrin` を入れる (`DISFRIN = 0` が既定) が、K64 の `tthine` (`tturne.f:592`、`qtwiss.f:1729` の `qthin`) は入れない。フリンジの写像は座標の 3 次なので、軌道が 6 極の中心を通るとき (設計軌道) は影響がない。軌道が中心からずれる (ステアリング、ずれた 6 極、分散軌道) と、CALC の行列・COD が追跡の 1 周写像と合わない。**K64 の癖 (欠落)**。v2 は追跡と同じ `ttfrin` を CALC・Emittance[] の COD の写像とヤコビアンに入れる (`lattice._thin_fringe_flag`)。k64 は K64 のまま。docs/30 の B291。この章では、`tthin` と `tthine` の差をフリンジの有無で数値化していない (未確認)。 |
| CALC・Emittance[] の `tthine` に放射がない (B240) | 追跡の `tthin` は $L \neq 0$ の多極で 1 刻みの放射 (`tradk`) をするが、K64 の `tthine` は `enarad = .false.` (放射せず、Emittance[] の減衰・量子励起とスピンの歳差に多極が入らない)。v2 は追跡と同じ 1 刻みの `tradke` を入れる。docs/30 の B240。この章では扱っていない (放射の章で)。 |
| 追跡と CALC の分割の次数の違い | 追跡 4 次、CALC 2 次 (誤りではないが、軌道がずれると行列が追跡のヤコビアンからずれる。§7)。 |
| `ttfrin` の `cx` の意味 | `tthin` の `cx` は $\bar z^k$、`ttfrin` の `cx` は $w = x + iy$ (§15。同じ名前で別の量。誤りではない)。 |
| `tthin` の `kord == 1` の分岐 (`tquad.f:238-246`) | `aki` を設定する行が分岐の中に無く (一般の分岐は `aki=akf/pr*2` を設定する)、`tthin` を `nord = 4`・`al ≠ 0` で呼ぶ箇所はソースを読んだ範囲に無い (`tquad.f:23`、`tquads.f:26` は `al = 0`)。到達しないコード (K64 のソースツリー全体の `tthin` の検索で、呼び出しは 5 つで、`nord = 4` は `al = 0` の 2 つだけ。§6.3。名前を実行時に組み立てる間接呼び出しは検索できない。誤りの有無は未確認)。 |
| 薄い多極の `DISKIN` | pysad は `DISKIN` を使わない ([02_drift.md](02_drift.md) 付録 A)。薄い多極の `tthin` の中のドリフトは常に `pxy2dpz` の厳密な式 (`DISKIN` は多極には見ない。K64 の確認は `DRIFT`・`QUAD` の分岐だけ)。 |
