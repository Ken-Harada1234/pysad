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

 だけ変える (K64 の `px -= aki*Re(cx)`、`py -= aki*Im(cx)`、`cx = (x - iy)^k`、`aki = ak/k!`)。リー演算子 $\exp(:f:)$ ($f = -V$) の級数は 2 項で止まる (§4)。シンプレクティック。pysad の追跡と $9 \times 10^{-19}$ 以下で一致 ([thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (c))。
- **長さのある多極**は、ドリフト (厳密) とキックの積に分割する。K64 の**追跡** (`tthin`) は 4 次の力勾配型: ドリフト $\alpha_1L$、キック 1/2、ドリフト $\alpha L$、**補正キック**、ドリフト $\alpha L$、キック 1/2、ドリフト $\alpha_1L$ ($\alpha = 1/\sqrt{12}$、$\alpha_1 = \tfrac12 - \alpha$。補正キックは $|\nabla V|^2$ ($\beta$) と $V_jV_kV_{jk}$ ($\gamma$) の座標だけの関数。式 (7.9)〜(7.12))。**CALC** (`tthine`) は、k64 では補正なしの 2 次の分割 ($L/6$、キック 1/2、$2L/3$、キック 1/2、$L/6$)。**v2 は CALC も追跡と同じ 4 次** (docs/30 B326、§7)。
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
| 多極の次数 $k$ ($V \propto (x - iy)^{k+1}$。四極 1、6 極 2、…) | `SEXT`・`OCT`・`DECA`・`DODECA` | `kord = nord/2 - 1` (`tthin`、`tthine`)。`nord` = 6 (6 極)、8、10、12 | `kord`。`_THIN_KEY` (`track6.py:889`) |
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

**級数**: (C.4)(C.5) から、$:f:\,g = [f, g]$:

$$
:f:\,x = [f, x] = -\frac{\partial f}{\partial p_x} = 0,\qquad
:f:\,y = 0,\qquad
:f:\,p_x = [f, p_x] = \frac{\partial f}{\partial x},\qquad
:f:\,p_y = \frac{\partial f}{\partial y}
\tag{7.7}
$$

($[f, p_x] = \partial_xf\,\partial_{p_x}p_x - \partial_{p_x}f\,\partial_xp_x = \partial_xf$。) 2 回目は、$\partial_xf$ も座標だけの関数なので、$:f:\,\partial_xf = [f, \partial_xf] = 0$ ($f$ も $\partial_xf$ も $p$ を含まないので、括弧の各項が 0)。したがって

$$
\exp(:f:)\,x = x,\quad \exp(:f:)\,y = y,\quad
\exp(:f:)\,p_x = p_x + \frac{\partial f}{\partial x},\quad
\exp(:f:)\,p_y = p_y + \frac{\partial f}{\partial y}
\tag{7.8}
$$

$f = -V_k$ とおいて (7.5) を使うと (7.1): $\Delta p_x = -\partial V_k/\partial x = -\dfrac K{k!}\mathrm{Re}\,(x - iy)^k$、$\Delta p_y = -\partial V_k/\partial y = -\dfrac K{k!}\mathrm{Im}\,(x - iy)^k$。6 極 ($k = 2$) なら $\Delta p_x = -\tfrac K2(x^2 - y^2)$、$\Delta p_y = +K\,xy$。8 極なら $\Delta p_x = -\tfrac K6(x^3 - 3xy^2)$、$\Delta p_y = \tfrac K6(3x^2y - y^3)$。

([thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (b): sympy で $[f, x] = [f, y] = 0$、$[f, [f, p_x]] = [f, [f, p_y]] = 0$ (2 項で止まる)、新しい変数の括弧 $[x, p_x'] = [y, p_y'] = 1$、$[x, p_y'] = [y, p_x'] = [p_x', p_y'] = 0$ (シンプレクティック。最後の $[p_x', p_y'] = 0$ は $V_k$ が調和でなくても、$\partial_x\partial_yf = \partial_y\partial_xf$ から成り立つ)。)

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
| CALC `tthine` (`tthine.f:1`)。k64 (v2 は上の `tthin` と同じ 4 次) | $1/6$、キック 1/2、$2/3$、キック 1/2、$1/6$ | 2 次 | $[B,[B,A]]$ (3 次の力の勾配の項)。$[A,[A,B]]$ が $\tfrac1{72}$ で残る |
| (比較) Strang | 1/2、キック、1/2 | 2 次 | なし |

## 6. 追跡の 4 次の分割 (`tthin`)

### 6.1 構成

03 章の力勾配型 (3.9) を、多極 $V_k$ ((7.4)) に当てはめる。$A$ = ドリフト、$B$ = キック。$T = (p_x^2 + p_y^2)/(2p)$ (有効質量 $p$) で、補正は座標だけの関数 (3.15):

$$
f_C = c_3\,L^3\,\frac{|\nabla V|^2}{p} - 2c_5\,L^5\,\frac{V_jV_kV_{jk}}{p^2},\qquad c_3 = \frac\beta2,\quad c_5 = \frac\gamma2
\tag{7.9}
$$

(K64 では長さ $L$ の 1 つの刻み。$V = V_k$ は長さあたりの強さ $\kappa$ で書く。) 全体の並びは、ヘルプの $\exp(:F_{\rm in}:)\exp(:aL:)\exp(:H/2:)\exp(:bL:)\exp(:V:)\exp(:aL:)\exp(:H/2:)\exp(:bL:)\exp(:F_{\rm out}:)$ のうち、2 番目の $\exp(:aL:)$ と最後の $\exp(:bL:)$ を $\exp(:bL:)$、$\exp(:aL:)$ に入れ替えた前後対称の形 (原文の並びは非対称で、そのままでは 2 次にしかならない。付録 A) に当たる:

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

### 6.2 補正キックの具体形 (導出)

$|\nabla V|^2$ と $V_jV_kV_{jk}$ を、(7.5) から作る。$V_x + iV_y = \dfrac{\kappa}{k!}\bar z^k$ ($\bar z = x - iy$) なので

$$
|\nabla V|^2 = \left(\frac\kappa{k!}\right)^2 r^k,\qquad r \equiv x^2 + y^2
\tag{7.10}
$$

($|\bar z^k|^2 = r^k$)。$V_jV_kV_{jk} = \tfrac12\nabla V\cdot\nabla|\nabla V|^2$ ($\partial_k|\nabla V|^2 = 2V_jV_{jk}$) で、$\nabla|\nabla V|^2 = (\kappa/k!)^2\,k\,r^{k-1}\cdot2(x, y)$、$\nabla V\cdot(x, y) = \dfrac\kappa{k!}(x\,\mathrm{Re}\,\bar z^k + y\,\mathrm{Im}\,\bar z^k) = \dfrac\kappa{k!}\,\mathrm{rcx1}$ ($\mathrm{rcx1} \equiv x\,\mathrm{Re}\,\mathrm{cx} + y\,\mathrm{Im}\,\mathrm{cx}$、$\mathrm{cx} = \bar z^k$。K64 の変数名) から

$$
V_jV_kV_{jk} = k\left(\frac\kappa{k!}\right)^3 r^{k-1}\,\mathrm{rcx1}
\tag{7.11}
$$

(K64 では、補正キックは、中央の 1 つ目のドリフト ($\alpha L$) のあとの座標 $(x_i, y_i)$ (= `xi`、`yi`) で評価し、そのあとに 2 つ目の $\alpha L$ のドリフトが来る。) $\mathrm{aki} \equiv \dfrac{\mathrm{ak}}{k!\,p}$ ($\mathrm{ak} = \kappa L$。$p$ は全運動量。K64 の補正キックの中の `aki = akf/pr*2`) とすると、(7.9) は

$$
f_C = p\,L\,c_3\,\mathrm{aki}^2\,r^k - 2c_5\,k\,p\,L^2\,\mathrm{aki}^3\,r^{k-1}\,\mathrm{rcx1}
$$

(第 1 項: $c_3L^3(\kappa/k!)^2r^k/p = c_3L\,(\kappa L/k!)^2r^k/p = c_3\,p\,L\,\mathrm{aki}^2 r^k$。第 2 項: $2c_5L^5k(\kappa/k!)^3 r^{k-1}\mathrm{rcx1}/p^2 = 2c_5\,k\,p\,L^2\mathrm{aki}^3\cdots$)。

規格化した運動量の更新 $\Delta q_x = (1/p)\,\partial f_C/\partial x$ ($\Delta p_x = \partial f/\partial x$ (7.8)、$q = p_x/p$)。偏微分を 1 項ずつ: $\partial_xr^k = 2k\,r^{k-1}x$、$\partial_x(r^{k-1}) = 2(k - 1)r^{k-2}x$、$\partial_x\mathrm{rcx1} = (k + 1)\,\mathrm{Re}\,\mathrm{cx}$ ($\partial_x\mathrm{Re}\,\bar z^k = k\,\mathrm{Re}\,\bar z^{k-1}$、$\partial_x\mathrm{Im}\,\bar z^k = k\,\mathrm{Im}\,\bar z^{k-1}$ から、$x\,\mathrm{Re}\,\bar z^{k-1} + y\,\mathrm{Im}\,\bar z^{k-1} = \mathrm{Re}\,\bar z^k$ を使う)。

$$
\begin{aligned}
\Delta q_x &= \mathrm{aki}^2\Bigl[\,f_1\,r^{k-1}x - \mathrm{aki}\,\bigl(f_2\,r^{k-2}\,x\,\mathrm{rcx1} + f_3\,r^{k-1}\,\mathrm{Re}\,\mathrm{cx}\bigr)\Bigr],\\
\Delta q_y &= \mathrm{aki}^2\Bigl[\,f_1\,r^{k-1}y - \mathrm{aki}\,\bigl(f_2\,r^{k-2}\,y\,\mathrm{rcx1} + f_3\,r^{k-1}\,\mathrm{Im}\,\mathrm{cx}\bigr)\Bigr],\\
\Delta z &= -\frac{\partial f_C}{\partial p} = \mathrm{aki}^2\,r^{k-1}\Bigl(f_4\,r - f_5\,\mathrm{aki}\,\mathrm{rcx1}\Bigr)
\end{aligned}
\tag{7.12}
$$

係数は

$$
f_1 = \beta k L,\quad f_2 = 2\gamma k(k - 1)L^2,\quad f_3 = \gamma k(k + 1)L^2,\quad f_4 = \tfrac12\beta L,\quad f_5 = 2\gamma kL^2
$$

(第 1 項の係数: $2k\,c_3L\,\mathrm{aki}^2r^{k-1}x$ から $f_1 = 2kc_3L = \beta kL$。第 2 項: $-2c_5kL^2\mathrm{aki}^3[2(k-1)xr^{k-2}\mathrm{rcx1} + (k+1)r^{k-1}\mathrm{Re}\,\mathrm{cx}]$ から $f_2 = 4c_5k(k-1)L^2 = 2\gamma k(k-1)L^2$、$f_3 = 2c_5k(k+1)L^2 = \gamma k(k+1)L^2$。$\Delta z$: $\partial_p(1/p) = -1/p^2$ から $-\partial f_C/\partial p = c_3L\,\mathrm{aki}^2r^k - 4c_5kL^2\mathrm{aki}^3r^{k-1}\mathrm{rcx1}$、$f_4 = c_3L \cdot$ ($r^k = r\cdot r^{k-1}$) $= \beta L/2$、$f_5 = 4c_5kL^2 = 2\gamma kL^2$。)

これは K64 の `tquad.f:209-216` ($k = 2$: `f1 = beta*2*al`、`f2 = 4*gamma*al**2`、`f3 = gamma*6*al**2`、`f4 = .5*beta*al`、`f5 = f2`)、`:261-267` (一般の $k$: `f1 = beta*kord*al`、`f2 = 2*gamma*kord*(kord-1)*al**2`、`f3 = gamma*kord*(kord+1)*al**2`、`f4 = .5*beta*al`、`f5 = 2*gamma*kord*al**2`) と**同じ式**。$k = 2, 3$ で sympy で全成分 ($\Delta q_x$、$\Delta q_y$、$\Delta z$) の差が 0 ([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (d))。$k = 1$ (薄い四極、`kord == 1`) の式 (`tquad.f:239-245`: `f1 = beta*al`、`f3 = gamma*2*al**2`、`f4 = f1/2`、`f5 = f3`) は、$r^{k-1} = 1$、$f_2$ の項が無い ($\mathrm{rcx1} = x^2 - y^2$)。この章では $k \geq 2$ だけ数値で確認した。

### 6.3 数値確認 (厳密な積分との比較)

[thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) (a)(b)。厳密な運動 $H = -\sqrt{1 - p_x^2 - p_y^2} + V_k$ ($\delta = 0$、$p = 1$) を DOP853 で積分した結果との差 (最大の成分)。初期値 $(x, p_x, y, p_y) = (20\,\mathrm{mm}, 1\,\mathrm{mrad}, -10\,\mathrm{mm}, -0.5\,\mathrm{mrad})$、$\kappa = 400$ /m³ (6 極)、$4 \times 10^4$ /m⁴ (8 極) 一定、$L$ を変える (積分強度 $\kappa L$ を $L$ とともに変える):

| 6 極、方法 | $L = 0.0125$ | $0.025$ | $0.05$ | $0.1$ | $0.2$ |
|---|---|---|---|---|---|
| pysad の追跡 (`tthin`) | $4.5 \times 10^{-13}$ | $8.6 \times 10^{-12}$ | $2.1 \times 10^{-10}$ | $6.9 \times 10^{-9}$ | $3.1 \times 10^{-7}$ |
| 補正から $\gamma$ の項を除く | $2.6 \times 10^{-12}$ | $8.3 \times 10^{-11}$ | $2.7 \times 10^{-9}$ | $8.7 \times 10^{-8}$ | $3.0 \times 10^{-6}$ |
| 補正なし ($\alpha_1$、$\alpha$ の割り振りだけ) | $3.5 \times 10^{-8}$ | $2.8 \times 10^{-7}$ | $2.2 \times 10^{-6}$ | $1.8 \times 10^{-5}$ | $1.4 \times 10^{-4}$ |
| CALC の分割 (`tthine` 型) | $5.7 \times 10^{-10}$ | $5.7 \times 10^{-9}$ | $6.3 \times 10^{-8}$ | $9.4 \times 10^{-7}$ | $2.2 \times 10^{-5}$ |
| Strang (2 次) | $2.6 \times 10^{-7}$ | $2.1 \times 10^{-6}$ | $1.7 \times 10^{-5}$ | $1.3 \times 10^{-4}$ | $1.1 \times 10^{-3}$ |

- **$\gamma$ の項を除いた 4 次の方法**は、$L$ を 2 倍にして比 32.0〜32.5 ($L^5$ の 4 次)。**補正なし**と **Strang** は比 8.0 (2 次)。**追跡 (`tthin`)** は、$L \ge 0.05$ で比が 32 を超える (33〜45): $\gamma$ の項が 5 次の座標だけの項 ($V'^2V''$ 型) を消すので、小さい傾きの粒子では 5 次の残りが小さく、7 次の項が先に見える (小さい $L$ の 19〜24 は、厳密な積分の精度 $\sim 10^{-13}$ に近い)。$L = 0.05$ で、$\gamma$ の項を入れると誤差が 13 倍小さくなる ($2.7 \times 10^{-9} \to 2.1 \times 10^{-10}$)。
- **傾きが大きい粒子** ($p_x = 30$ mrad、$p_y = -20$ mrad): $\gamma$ の項を入れても誤差はほぼ同じ ($L = 0.05$ で $9.7 \times 10^{-9}$ 対 $8.1 \times 10^{-9}$)。5 次の誤差は $p^2$・$p^4$ の項 ((3.13)。運動量を含み、座標だけの補正では消せない) が支配する。したがって、$\gamma$ の項は、**傾きの小さい粒子 (リング内の大部分の粒子) の誤差を減らす**もの。
- **CALC の分割 (`tthine` 型)**: 傾きが小さいとき (1 mrad)、比 10〜23 (3 次の項 $[A,[A,B]]$ は傾きの 2 乗に比例して小さく、5 次の項が見える)。傾きが大きいとき (30 mrad)、比 8.1〜9.1 (3 次の項が見える。2 次)。
- 8 極も同じ傾向 ($\gamma$ なしの比 31.8〜33.5、補正なしの比 8.0〜8.2)。
- **K64 の式を書き直した関数** (この script の `tthin_k64`) と pysad の追跡の最大差は $6.9 \times 10^{-18}$ (6 極)、$1.3 \times 10^{-18}$ (8 極)。

## 7. CALC (行列) の `tthine`

CALC・Emittance[] の薄い多極は、COD $\xi$ を通る粒子の写像と、そのヤコビアン ([02_drift.md](02_drift.md) §8 の `tdrife1` と同じ考え方) を、k64 は 2 次の分割 (§5 の表) で作る (v2 の 4 次は下の段落)。ドリフトは `tdrife0` (厳密。$dv$ の項を含む)。キックは正準運動量のまま ($1/p$ なし): $\Delta p_x = -\mathrm{aki}\,\mathrm{Re}\,\mathrm{cx}$、$\Delta p_y = -\mathrm{aki}\,\mathrm{Im}\,\mathrm{cx}$ (`aki = ak/fact(kord)`、半分ずつの 2 回は `aki*.5`)。キックのヤコビアンは、$\mathrm{cx} = \bar z^k$ の微分 ($\partial_x\bar z^k = k\bar z^{k-1}$、$\partial_y\bar z^k = -ik\bar z^{k-1}$) から ($g_r + ig_i \equiv k\,\bar z^{k-1}$):

$$
\frac{\partial(\Delta p_x, \Delta p_y)}{\partial(x, y)} = -\,\mathrm{aki}\begin{pmatrix} g_r & g_i\\ g_i & -g_r\end{pmatrix}
\tag{7.13}
$$

(K64 `tthine.f` の `trans1(2,1) = -aki*dble(cx1)`、`trans1(2,3) = -aki*imag(cx1)`、`trans1(4,1) = -aki*imag(cx1)`、`trans1(4,3) = aki*dble(cx1)`、`cx1 = k*(x - iy)^(k-1)`。) 行列は、これと厳密なドリフトの行列 ((2.13)) を、分割の順に掛けたもの。

- **設計軌道上** ($\xi = 0$): キックのヤコビアンは 0 ($k \geq 2$ で $\bar z^{k-1} = 0$) で、行列はドリフトそのもの。追跡と CALC の差は $1.3 \times 10^{-11}$ (L = 0.1 m、4×4 と $z$ の行。[thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) (c))。
- **軌道がずれたとき**: 追跡 (4 次) と CALC (2 次) は分割が違うので、ヤコビアンが違う。厳密な積分のヤコビアンとの最大差は、6 極 ($L = 0.1$ m、$\kappa = 400$ /m³、COD $x = 20$ mm など) で、追跡 $4.2 \times 10^{-6}$、CALC $2.3 \times 10^{-4}$。追跡と CALC の差は $2.3 \times 10^{-4}$ (同 (c))。**CALC の行列は、追跡のヤコビアンに対して、COD がずれるほど 2 次の分割の誤差だけ違う**。シンプレクティック性は、CALC の行列が $2 \times 10^{-16}$ (解析的)、追跡が $3 \times 10^{-12}$ (中心差分の精度)。
- **v2 (B326)**: $L \neq 0$ の CALC・Emittance[] は、追跡と同じ 4 次の分割 (§6) にした。COD の写像は、ドリフト $\alpha_1L$ → キック 1/2 → ドリフト $\alpha L$ → 補正キック → ドリフト $\alpha L$ → キック 1/2 → ドリフト $\alpha_1L$。補正キックは座標 $(x, y)$ と $\delta$ だけの関数 (正準運動量では、追跡の規格化した運動量の補正に $(1+\delta)$ を掛けたもの: $\Delta p_x = \dfrac{a^2}{p}f_1 r^{k-1}x - \dfrac{a^3}{p^2}\bigl(f_2 r^{k-2}x\,R + f_3 r^{k-1}\mathrm{Re}\,\bar z^k\bigr)$、$\Delta z = \dfrac{a^2}{p^2}f_4 r^{k}- \dfrac{a^3}{p^3}f_5 r^{k-1}R$、$a = $ `ak/fact(k)`、$p = 1+\delta$、$r = x^2 + y^2$、$R = \mathrm{Re}\,\bar z^{k+1}$、$f_1$〜$f_5$ は (7.12) と同じ)。行列は、各段 (厳密なドリフト、キック (7.13)、補正キックの $(x, y, \delta)$ での解析微分) の積 (連鎖律)。$k = 2$〜5、$L = 0.1$〜$0.4$ m で、COD の写像は追跡の 1 刻みと $10^{-14}$ 以内で一致し、行列の厳密な積分との差は $L^5$ に比例 (k64 の 2 次は $L^3$。例: 6 極 $\kappa = 50$ /m³、$L = 0.2$ m で $5\times10^{-8}$ 対 $3.4\times10^{-5}$)。シンプレクティック性は $10^{-16}$。設計軌道上 ($\xi = 0$) では補正キックのヤコビアンも 0 で、行列は k64 と同じドリフト。MULT (`tmulke`) と BEND の 6 極成分は対象外 (従来のまま)。\n- 厚い多極の磁場の勾配が小さい (リングの 6 極は $K_2L \sim 10$〜$100$ /m²、$L \sim 0.1$ m) とき、この差は実用上小さい。

## 8. フィードダウン (ずれた多極)

多極が $(d_x, d_y)$ ずれて置かれると、粒子の座標 $(x, y)$ は磁石の中心から $(x - d_x, y - d_y)$ になる (pysad の `_tent` は K64 の `inc/TENT.inc` の翻訳で、`x -= dx`、`y -= dy` のあとに $\theta$ の回転。`TENT.inc` 自体は読んでいない)。薄いキック (7.1) は、$x - iy \to (x - d_x) - i(y - d_y) \equiv \bar z - \bar d$ ($\bar d = d_x - id_y$):

$$
\Delta(p_x + ip_y) = -\frac K{k!}\,(\bar z - \bar d)^k = -\frac K{k!}\sum_{j=0}^k\binom kj\,\bar z^{k-j}(-\bar d)^j
\tag{7.14}
$$

($j = 0$ が元の多極、$j = 1$ が 1 つ低い次数の多極、$j = 2$ がもう 1 つ低い次数…。) 6 極 ($k = 2$) の場合: $\Delta p_x = -\tfrac K2\bigl[(x - d_x)^2 - (y - d_y)^2\bigr]$、$\Delta p_y = K(x - d_x)(y - d_y)$ を展開して、

$$
\begin{aligned}
\Delta p_x &= -\tfrac K2(x^2 - y^2) + K(d_x\,x - d_y\,y) - \tfrac K2(d_x^2 - d_y^2),\\
\Delta p_y &= K\,xy - K(d_y\,x + d_x\,y) + K\,d_xd_y
\end{aligned}
\tag{7.15}
$$

- 第 1 項: 元の 6 極。
- 第 2 項: **四極 (フィードダウン)**。$d_y = 0$ なら $\Delta p_x = +Kd_x\,x$、$\Delta p_y = -Kd_x\,y$ で、四極 ($\Delta p_x = -K_1x$、$\Delta p_y = +K_1y$) の $K_1 = -K\,d_x$。$d_x$ が大きいほど強い四極として働く。$d_y \neq 0$ の項は、スキュー四極 ($\Delta p_x \propto y$、$\Delta p_y \propto x$) で、強さの大きさは $K\,d_y$。
- 第 3 項: **二極のキック** ($x, y$ によらない)。$-K(d_x^2 - d_y^2)/2$、$K\,d_xd_y$。

(sympy で (7.15) の展開との差 0。pysad の SEXT ($K_2 = 40$、DX = 1.5 mm、DY = −0.8 mm、薄い) の追跡と、$(x - \mathrm{DX}, y - \mathrm{DY})$ を使う閉じた式の最大差 $8.7 \times 10^{-19}$。軌道が $x = y = 0$ のときの二極のキックは $\Delta p_x = -3.220 \times 10^{-5}$、$\Delta p_y = -4.800 \times 10^{-5}$ で、(7.15) の第 3 項と一致。[thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) (d)。)

**意味**: 取り付け誤差の座標変換 (ずらす→キック→戻す) と、分割した写像の合成だけで、フィードダウンは**自動的に**出る。ハミルトニアンに低次の項を足す必要はない ([11_misalignment_and_coord.md] (予定) で一般論)。長さのある多極でも同じ (4 次の分割の中の座標がずれている)。

## 9. 端部の非線形フリンジ

### 9.1 Maxwell 方程式から

長さ $L \neq 0$ の多極の端 ($s = 0$) で、強さが $\kappa(s)$ ($0 \to \kappa_0$) と変わるとき、磁場は Maxwell 方程式から縦成分と高次の項が決まる ([04_quad.md](04_quad.md) §7.2 の四極の場合の一般化)。$2m$ 極 ($m = k + 1$。四極 $m = 2$、6 極 $m = 3$) の、**法線成分の多極**のスカラーポテンシャル ($\vec b = \nabla\Phi$、$w = x + iy$) は

$$
\Phi = \kappa(s)\,\frac{\mathrm{Im}\,w^m}{m!} - \kappa''(s)\,\frac{|w|^2\,\mathrm{Im}\,w^m}{4\,(m + 1)!} + \cdots
\tag{7.16}
$$

(**検証**: $\nabla_\perp^2(|w|^2\mathrm{Im}\,w^m) = 4(m + 1)\,\mathrm{Im}\,w^m$ ($\nabla_\perp^2 = 4\partial_w\partial_{\bar w}$ を $w\bar w\cdot(w^m - \bar w^m)/(2i)$ に作用させる)。したがって $\nabla_\perp^2\Phi = -\kappa''\mathrm{Im}\,w^m/m!$ で、$\partial_s^2\Phi = \kappa''\mathrm{Im}\,w^m/m!$ を足すと 0。$m = 2, 3, 4$ で sympy で確認: [thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (a)。) 四極 ($m = 2$) では $\mathrm{Im}\,w^2 = 2xy$ で、$\Phi = \kappa\,xy - \kappa''\,xy(x^2 + y^2)/12$ ((4.24) に一致)。磁場は $\vec b = (\partial_x\Phi, \partial_y\Phi, \partial_s\Phi)$ で、$b_s = \kappa'\,\mathrm{Im}\,w^m/m! + \cdots$。力は ($p_z \simeq p$) $F_x = -\partial_y\Phi$、$F_y = +\partial_x\Phi$ ($d\vec p/ds = \vec p\times\vec b/p_z$ から。$b_s$ の項は運動量の変化 (K の 1 次) に効く)。

### 9.2 位置のずれ

硬い端 $\kappa = \kappa_0\Theta(s)$ では $\kappa'' = \kappa_0\delta'(s)$。力が $\delta'(s)$ に比例するので、運動量は $\delta(s)$ 型のスパイクを持ち、位置は端で階段状にずれる ($\int\delta'(s)g(s)\,ds = -g'(0)$ から、運動量の $\delta(s)$ 項の係数が位置のずれの元)。$\kappa''$ の項の力は

$$
F_x + iF_y = \frac{\kappa''}{4(m + 1)!}\,(\partial_y - i\partial_x)\bigl(|w|^2\mathrm{Im}\,w^m\bigr) = -i\,\frac{\kappa''}{4(m + 1)!}\cdot2\,\partial_{\bar w}\bigl(|w|^2\mathrm{Im}\,w^m\bigr)
$$

($\partial_x + i\partial_y = 2\partial_{\bar w}$。) $2\partial_{\bar w}[w\bar w(w^m - \bar w^m)/(2i)] = (w/i)\bigl[w^m - (m + 1)\bar w^m\bigr] = -i\,w\bigl[w^m - (m + 1)\bar w^m\bigr]$ なので、位置のずれ (運動量の $\delta$ スパイクの係数を $p$ で割り、$\kappa'' \to \kappa_0\delta'$):

$$
\Delta w = -\frac{\kappa_0}{4\,(m + 1)!\,p}\,w\bigl[w^m - (m + 1)\bar w^m\bigr] = -\frac{\kappa_0}{4\,m!\,p}\,w\left(\frac{w^m}{m + 1} - \bar w^m\right)
\tag{7.17}
$$

(式 (7.2)。sympy で、(7.16) の $\kappa''$ の項の力から出した $\Delta w$ と (7.17) の差が $m = 2, 3, 4$ で 0。) 四極 ($m = 2$) で $\Delta w = -\dfrac{\kappa_0}{8p}\,w\left(\dfrac{w^2}{3} - \bar w^2\right)$ の実部は $\Delta x = \dfrac{\kappa_0}{12p}(x^3 + 3xy^2)$ で (4.27)。出口は $\kappa' \to -\kappa_0\delta$ で符号が逆 ($\kappa_0 \to -\kappa_0$)。

**K64 の `ttfrin` (`nord = 2m`、$m \geq 3$、`tquad.f:325-452`) との対応**: `aki = ak/al/fact(kord)/4/(1+g)` ($\mathrm{kord} = m$) は $\kappa_0/(4\,m!\,p)$、`cx = dcmplx(x, y)` は $w$、`cz1 = cx**(kord-1)`、`cz = cz1*cx` $= w^m$、`ca = -aki*cx*(cz/an(kord+1) - conjg(cz))` ($\mathrm{an}(m+1) = m + 1$) は $-\mathrm{aki}\,w(w^m/(m + 1) - \bar w^m)$、すなわち $\Delta w$ そのもの。`cx1 = cx + ca` が新しい $(x, y)$ ($m = 3$ の `ttfrin` の分岐は同じ式で、数が `4`・`3` と直接書いてある)。pysad の `ttfrin` ($m = 2$〜6) の位置のずれと (7.17) は、$m = 3, 4, 5, 6$ で $4 \times 10^{-19}$ 以下で一致する (丸め誤差。$\kappa_0$ を半分にしても差が減らないので、(7.17) の $K$ の 1 次の項以外は無い。**多極 ($m \geq 3$) の点変換は (7.17) そのもの**)。四極 ($m = 2$) だけ、母関数 (4.28) から来る $K$ の 2 次の項 ($t$) があり、差が $\kappa_0$ の 2 乗で減る (比 4.0。[thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) (b))。

### 9.3 運動量と $z$ (正準変換)

位置が $w \to w + \Delta w(w)$ という**点変換**で、運動量は、点変換のヤコビアン $J = \partial(x', y')/\partial(x, y)$ の逆転置 $J^{-\top}$ で更新する (点変換の母関数 $F_2 = \vec x'(\vec x)\cdot\vec p'$ から、$\vec p = J^\top\vec p'$)。複素数で書くと、$A \equiv \partial_w\Delta w = -\mathrm{aki}\,(w^m - \bar w^m) = -2i\,\mathrm{aki}\,\mathrm{Im}\,w^m$、$B \equiv \partial_{\bar w}\Delta w = m\,\mathrm{aki}\,w\bar w^{m-1}$ で、$\det J = |1 + A|^2 - |B|^2 = 1 + a^2 - m^2\mathrm{aki}^2|w|^{2m}$ ($a \equiv 2\,\mathrm{aki}\,\mathrm{Im}\,w^m$)。K64 の `d`、`a`、`cp1 = ((1 + i a) cp - m aki cz1 conj(cx cp))/d` ($\mathrm{cp} = p_x - ip_y$、`d` $= \det J$) は、これを書いたもの。$z$ は、位置のずれが運動量に依らないので、母関数の $p$ 依存性 ($\mathrm{aki} \propto 1/p$) から

$$
z \to z - \Delta\vec r\cdot\vec q_{\rm new} = z - \mathrm{Re}\bigl(\Delta w\cdot\mathrm{cp1}\bigr)
\tag{7.18}
$$

([04_quad.md](04_quad.md) §7.3 の四極と同じ構造。`z = z - dble(ca*cp1)`。$q$ 変数なので $q_{\rm new} = p_{\rm new}/p$。) 点変換から作る写像なので、シンプレクティック性は構成から厳密に成り立つ。

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

放射を入れるとき ($\mathrm{krad}$)、`tthin` は `tradk(np, ..., al, 0)` を、本体全体 ($L$) を 1 刻みとして**最後に 1 回**呼ぶ (入口の位置 `pxr0`・`pyr0`・`zr0` を保存して、全体の変位から光子の放射を計算する。放射の式は [13_radiation.md] (予定)。`tsetpcvt` で光子の出る点を $L$ の範囲に設定)。K64 の `tthine` (CALC) は `enarad = .false.` で放射しない (Emittance[] の減衰・励起に多極が入らない。付録 A)。

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
| 追跡の 6 極 ($k = 2$) | `track6._tthin_sext` (`track6.py:892`) | `tthin` | `ak = K*L` ($L \neq 0$)、$L = 0$ は `ak = K`。係数は `_TT_ALPHA` など (`:883-886`) |
| 追跡の 8 極以上 ($k = 3, 4, 5$) | `track6._tthin_gen` (`:740`) | `tthin` の一般 `kord` | `f1`〜`f5` の一般式 (`:771-776`) |
| 追跡の非線形フリンジ | `track6._ttfrin6` (`:815`、$m = 3$)、`_ttfrin_gen` (`:1932`、$m \geq 4$)、`_ttfrin_n` (`:1954`) | `ttfrin` | `a`、`ca`、`cx1`、`d`、`cp1`。`z -= (ca*cp1).real` |
| 追跡 (numba) | `track6n.k_thin` (`track6n.py:643`)、`ttfrin6_1` (`:582`)、`ttfrin_gen_1` (`:603`) | 同上 | |
| CALC の COD の写像 | `lattice._tthine_cod` (`lattice.py:1755`)、`elkern._tthine_map` (`elkern.py:360`) | `tthine` | 正準運動量のキック |
| CALC のヤコビアン | `elkern._jac_thin` (`elkern.py:330`)、`_jac_kick` (`:302`)、`emit6d._jac_thin` (`emit6d.py:469`) | `tthine` | 式 (7.13) |
| CALC (Emittance[]) の全体 | `emit6d._tthine` (`emit6d.py:1778`)、`_tthine_gen` (`:1794`)、`elkern.k_thin` (`:558`) | `tturne.f:592` | v2 は非線形フリンジ (`_lat._thin_fringe_flag`、`lattice.py:1783`) と、$L \neq 0$ の放射 (`krad`) を足す。k64 は K64 のまま (入れない) |
| 共鳴駆動項 (6 極の項) | `resdrive._thin_funcs` (`resdrive.py:277`) | — | [18_nonlinear_bridge.md] (予定) |

## 13. 確認スクリプトと出力

スクリプトはリポジトリのルートから `py -3.12 physics_guide/scripts/<名前>.py` で動かす (数秒〜十数秒。sympy・scipy が要る)。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `thin_01_potential_kick.py` | (a) 多極のポテンシャル (調和、規格化、複素形)。(b) $\exp(:f:)$ の級数が 2 項で止まる、シンプレクティック。(c) pysad の薄い多極 ($k = 2$〜5、$\delta \neq 0$) と (7.1)。(d) フィードダウン (7.15) と pysad の DX・DY | [thin_01_potential_kick.out.txt](outputs/thin_01_potential_kick.out.txt) |
| `thin_02_fourth_order.py` | (a) pysad の `tthin` と K64 の式の書き直しの一致。(b) 厳密な積分との差と $L$ を 2 倍にしたときの比 (4 次・2 次、$\gamma$ の効果、傾きの大小)。(c) ヤコビアン (追跡・CALC・厳密) | [thin_02_fourth_order.out.txt](outputs/thin_02_fourth_order.out.txt) |
| `thin_03_fringe.py` | (a) 磁場の調和性と位置のずれ (7.17) (sympy)。(b) pysad の `ttfrin` ($m = 2$〜6) との比較。(c) 運動量・$z$・シンプレクティック性。(d) Maxwell の場の積分との比較 (6 極・8 極) | [thin_03_fringe.out.txt](outputs/thin_03_fringe.out.txt) |
| `integ_01_bch_force_gradient.py` (03 章) | 補正キック (7.12) と K64 の `tthin` の式の記号的な一致 ((d)) | [integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) |

## 14. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV/$c$ ($\delta = 0$ が主)、6 極・8 極、$L \le 0.4$ m、振幅 数 mm〜20 mm、傾き 1〜30 mrad。多極は理想 (法線成分だけ、ポテンシャル (7.4))。

**この章の範囲外 (別の章)**: 多極 MULT (`tmulti`・`tmulte`。スキュー・複数の成分・ANGLE) は MULT の章 (予定。08)。端部の一般論 (磁場の立ち上がり F1・F2、BEND のフリンジ) は端部の章 (予定。06)。四極の本体は [04_quad.md](04_quad.md)。放射・スピンは別の章。

**未確認・限界**:
- **補正キックの式を確認した $k$**: 6 極 ($k = 2$)、8 極 ($k = 3$)。DECA・DODECA ($k = 4, 5$) は、式の一般形 ((7.12)) と pysad の一般 `kord` が同じ構造であることを確認していない (記号計算は $k = 2, 3$ だけ。数値の薄いキック (7.1) は $k = 2$〜5 で確認した)。薄い四極 ($k = 1$、`kord == 1` の分岐) の補正も確認していない。
- **補正キックの $\gamma$ の項は全体の誤差を必ず減らすとは限らない** (03 章 §8.3)。6 極の例では、傾きの小さい粒子で誤差を 13 倍減らし、傾きが大きい粒子ではほぼ効かない (§6.3)。一般の 1 自由度のポテンシャル ($x^4/4$) では、$\gamma$ の項を入れない方が誤差が小さかった (03 章 §9)。
- **K64 のコードとの比較は、pysad の Python 版と、K64 の式を書き直した関数を通して**行った。gs (K64 の実行) の出力との数値比較は、この章ではしていない。
- **非線形フリンジ**: K の 1 次 (硬い端の近似)。端の磁場の展開は 3 次まで ($b_s = \kappa'\mathrm{Im}\,w^m/m!$ など。$\kappa''''$ の項は入れていない)。Maxwell の場の数値積分との比較は、$\kappa_0$、$r$、$d$ の限られた範囲 (表)。端の幅 $d$ が孔径に近い場合や、$K$ の 2 次の項 (四極で見えた残り。[04_quad.md](04_quad.md) §7.2) は確認していない。$m \geq 3$ の運動量・$z$ の更新を Maxwell の場から直接導いた (7.17 の位置のずれのように) 検証は、$\Delta(\mathrm{ODE})$ との一致 (表) による。
- **フリンジの線形部分** (多極の F1 に当たる、端の立ち上がりの形による効果。四極の §7.4) は、この章では扱わない (SEXT・OCT … は F1・F2 を持たない。端の立ち上がりの形で決まる $K$ の 1 次の項は、設計軌道上では $\propto x^{k}$ の非線形で、確認していない)。
- **フィードダウン**の式 (7.15) は、薄いキックで確認した。長さのある多極での ROTATE・DROTATE・DZ との組合せ、`tchge` の癖 ([11_misalignment_and_coord.md] (予定)) は確認していない。
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
| ヘルプに無いもの | (1) 補正キックの具体形 (7.12)。(2) 追跡 (4 次) と CALC (2 次) の分割の違い。(3) 非線形フリンジの Maxwell 方程式からの導出と、運動量・$z$ の更新。(4) フィードダウンの式。 |

**K64 の実装の癖・誤りの疑い**

| 箇所 | 内容 |
|---|---|
| CALC・Emittance[] の `tthine` にフリンジがない (B291) | 追跡の `tthin` は非線形フリンジ `ttfrin` を入れる (`DISFRIN = 0` が既定) が、K64 の `tthine` (`tturne.f:592`、`qtwiss.f:1729` の `qthin`) は入れない。フリンジの写像は座標の 3 次なので、軌道が 6 極の中心を通るとき (設計軌道) は影響がない。軌道が中心からずれる (ステアリング、ずれた 6 極、分散軌道) と、CALC の行列・COD が追跡の 1 周写像と合わない。**K64 の癖 (欠落)**。v2 は追跡と同じ `ttfrin` を CALC・Emittance[] の COD の写像とヤコビアンに入れる (`lattice._thin_fringe_flag`)。k64 は K64 のまま。docs/30 の B291。この章では、`tthin` と `tthine` の差をフリンジの有無で数値化していない (未確認)。 |
| CALC・Emittance[] の `tthine` に放射がない (B240) | 追跡の `tthin` は $L \neq 0$ の多極で 1 刻みの放射 (`tradk`) をするが、K64 の `tthine` は `enarad = .false.` (放射せず、Emittance[] の減衰・量子励起とスピンの歳差に多極が入らない)。v2 は追跡と同じ 1 刻みの `tradke` を入れる。docs/30 の B240。この章では扱っていない (放射の章で)。 |
| 追跡と CALC の分割の次数の違い (B326) | 追跡 4 次、K64 の CALC 2 次 (軌道がずれると行列が追跡のヤコビアンからずれる。§7)。v2 は CALC も 4 次。 |
| `ttfrin` の `cx` の意味 | `tthin` の `cx` は $\bar z^k$、`ttfrin` の `cx` は $w = x + iy$ (§15。同じ名前で別の量。誤りではない)。 |
| 薄い多極の `DISKIN` | pysad は `DISKIN` を使わない ([02_drift.md](02_drift.md) 付録 A)。薄い多極の `tthin` の中のドリフトは常に `pxy2dpz` の厳密な式 (`DISKIN` は多極には見ない。K64 の確認は `DRIFT`・`QUAD` の分岐だけ)。 |
