# 12. 行列から光学へ: 閉軌道、周期解、x–y 結合、Twiss

**この章で分かること**: 要素ごとの写像と行列 (02〜11 章) から、リングの光学関数 (閉軌道、分散、x–y 結合、Twiss パラメータ、拡張 Twiss) が何をどう計算して出るか。(1) 閉軌道の Newton 法 (`tcod`) と収束の速さ。(2) 1 周の行列からの周期解と分散。(3) x–y 結合の分解 (`qmdiag`): Edwards–Teng 型の変換 $R$ を、ヘルプの定義と K64 の式で対応づける。(4) Courant–Snyder の行列と Twiss の伝搬。(5) 色収差。(6) 6 次元の拡張 Twiss (ヘルプの $V = PBR_6H$)。放射を含む平衡 (CALC6D の $D$、$\Sigma$) は [14_emittance.md](14_emittance.md)。

**前提**: 各要素の写像と、`trans(6,12)` の 1:6 列 (行列 $T$) ([02_drift.md](02_drift.md)、[04_quad.md](04_quad.md)、[05_bend.md](05_bend.md))。

**読む順**: §1 (結論) → §3 (閉軌道) → §5 (結合) → §6 (Twiss)。

## 1. 結論

- **閉軌道**は、1 周の写像 $f$ の不動点 $f(c) = c$。K64 (`tcod.f`) は Newton 法で解く: 更新 $(T - 1)\,\Delta = c - f(c)$。$T$ は 1 周の行列。収束判定は成分ごとの重み `codw` つきの残差 $r = \sum[(c - f(c))_i/\mathtt{codw}_i]^2 < \mathtt{conv}$。非線形のリング (6 極とずれた四極つき) で、残差は $9\times10^{-5} \to 2\times10^{-6} \to 1\times10^{-9} \to 2.5\times10^{-16}$ と下がり、毎回 $r_{k+1} \simeq C\,r_k^2$ ($C \simeq 250$) で 2 次収束した (§3.2)。
- **分散**は、同じ不動点の方程式を $\delta$ について 1 次に展開した解 $\eta = (1 - M_4)^{-1}m_6$ ($M_4$ は 1 周の $4\times4$、$m_6$ は $\delta$ の列)。追跡で $\delta = \pm10^{-5}$ の閉軌道の差分をとった値と $10^{-6}$ 以内で一致 (例のリングで $\eta_x = 2.28993$ m)。
- **x–y 結合の分解** $R M R^{-1} = \mathrm{diag}(A, B)$ の $R = \begin{pmatrix}c\,I & J r^\top J\\ r & c\,I\end{pmatrix}$ ($c^2 + \det r = 1$) は、K64 の `qmdiag` の閉じた式で求まる。結合した安定な行列で、$R$ はシンプレクティック、$R^{-1}$ は $r$ の符号を逆にしたもの (ヘルプの通り)、非対角ブロックは $10^{-16}$。pysad の `periodic_r` と一致。(12.6) が Riccati 式を満たすことは代数的に証明した (§5.3)。
- **Twiss**は、ブロック $A$ ($= $ 1 周の 2×2) に $B A B^{-1} = $ 回転 (§6.1 で行列の積を成分ごとに展開、sympy でも確認)。伝搬は $\beta_2 = t_{11}^2\beta - 2t_{11}t_{12}\alpha + t_{12}^2\gamma$、$\alpha_2$ も同様 (§6.2 で $A_2 = TA_1T^{-1}$ の成分から導出、ドリフトで pysad の Twiss と $9\times10^{-16}$)。
- **色収差**: 四極だけのリングで、追跡の Jacobian から出した $d\nu/d\delta$ は pysad の `lat.chromaticity()` と一致し、**厚い四極の積分 $-\frac1{4\pi}\sum K L\,\bar\beta$ ($\bar\beta$ は四極内の $\beta$ の平均) と 4 桁で一致**する。薄いレンズの式 (β は四極の両端の平均) との 4〜4.5% の差は、四極内の $\beta$ の曲がり ($\bar\beta$ と端の平均の差が $\phi^2/6$ 程度、$\phi = \sqrt K L$) が、QF と QD の項の差で拡大されたもの (§7)。**偏向磁石のあるリングでは、`lat.chromaticity()` だけが追跡・厳密なハミルトニアンの写像と合わない** (§7。`Lattice.chromaticity()` の計算の問題と切り分けた)。
- 6 次元の拡張 Twiss の行列 $H$ (ヘルプ) は、$a^2 + \det H_x + \det H_y = 1$ のとき**シンプレクティック**: (1,1)・(1,3)・(3,3) ブロックを §8 で手で展開して示し、全 6 ブロックを記号計算と有理数で厳密に、数値 (浮動小数点) で $2\times10^{-16}$ で確認した。

## 2. 記号と単位

| 量 | SAD (ヘルプ) | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| 閉軌道 | `DX`・`DPX`・`DY`・`DPY`・`DZ`・`DDP` | `cod(1:6)` (正準運動量) | `cod` |
| 1 周の行列 $T$ | | `trans(6,6)` (`trans(6,12)` の 1:6 列) | `lat.transfer_matrix()`、`calc6d` |
| 結合の行列 $r$ | `R1 R2 R3 R4` | `r1 r2 r3 r4` (`qmdiag` の出力)、`c` ($=\mu$) | `optics.periodic_r` の `R` |
| $\mu$ | (ヘルプの $\mu$) | `c` | `a` |
| $\det r$ | `DETR` ($\det r \ge 0.625$ で別の形、DETR は $1.375 - \det r$) | | `detr` |
| Twiss | `AX BX AY BY`、`NX NY` (位相 [rad]) | | `twiss` の `ax bx ay by nx ny` |
| 分散 | `EX EPX EY EPY` (拡張: `ZX ZPX` など) | | `ex epx ey epy` |
| 収束の重み | `CODCONV` | `codw(6) = [1e-6, 1e-7, 1e-7, 1e-8, 1e-5, 1e-6]` (RF なしは 0.01 倍)、`conv0 = 1e-4` | `emit6d._tcod_radcod` |

$J = \begin{pmatrix}0&1\\-1&0\end{pmatrix}$ (2×2)。バー記号 $\bar X \equiv J X^\top J$ (2×2 では $\bar X = -\mathrm{adj}\,X$、すなわち $\bar X = \begin{pmatrix}-X_{22} & X_{12}\\ X_{21} & -X_{11}\end{pmatrix}$。$X\bar X = -\det X\cdot I$)。

## 3. 閉軌道 (`tcod`)

### 3.1 不動点の方程式と Newton 法

閉軌道 $c$ は、1 周の写像 $f$ について $f(c) = c$。$c$ の近くで $f(c + \Delta) = f(c) + T\Delta + O(\Delta^2)$ ($T = \partial f/\partial c$ は 1 周の行列)。$f(c + \Delta) = c + \Delta$ とおいて 1 次で打ち切ると

$$
(T - 1)\,\Delta = c - f(c)
\tag{12.1}
$$

K64 では `dcod1 = codi - codf` ($c - f(c)$)、`trs = trans - 1`、`tsolvg(trs, dcod1, dcod)` が (12.1) を解き、`codi = codi + fact*dcod` (`fact` $\le 1$ の減衰) で更新する (`tcod.f` の `1 loop`)。`tturne` が 1 周を追跡して `codf` と `trans` (6×12) を同時に作る。

- **収束判定**: $r = \sum_i[(c - f(c))_i/\mathtt{codw}_i]^2$。$r <$ `conv` (= `conv0` $= 10^{-4}$、放射テーパーのときは 100 倍) で終了。$r$ が前回より増えたら `fact` を半分にして戻る (`fact < fmin` $= 10^{-4}$ で打ち切り)。
- **RF のないとき** (`rfsw = .false.`): $z$ と $\delta$ は閉じなくてよいので、`trs` の 5・6 列を 0 にして 5 行 (`im = 5`) を解く (成分 5 は 0)。$\delta$ は与えた値 (`cod(6)`) で固定。
- **RF のあるとき** (`rfsw = .true.`): まず RF なしで 4 次元の軌道を求め ($\delta$ 固定)、そのあと $z$ を、1 周のエネルギー損失 $U_0$ と空洞の電圧のつりあいから決める: `cod(5) = asin((u0*pgev - vcacc)/vceff)/wrfeff - trf0` (同期位相。導出は [14_emittance.md](14_emittance.md))。`radcod` のとき、放射の損失を含めて $\delta$ ごと不動点を解く (RADCOD。ずれ $dz$ を `trf0` に移す)。

### 3.2 収束の速さ

非線形な項があるとき (6 極と、ずれた四極 DX = 0.3 mm)、原点での線形解 $(T-1)^{-1}(c - f(c))$ は 1 回では真の解に届かず、2 次収束で数回で収束する ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (a))。残差 $\max|f(c) - c|$ の列:

| 反復 $k$ (更新の前の点 $c_k$ での残差。$c_1 = 0$) | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 残差 $r_k = \max\|f(c_k) - c_k\|$ | $8.9\times10^{-5}$ | $2.0\times10^{-6}$ | $1.0\times10^{-9}$ | $2.5\times10^{-16}$ | $2.7\times10^{-19}$ |

**2 次収束の理由**: 真の解を $c^*$、$k$ 回目の近似を $c_k = c^* + e_k$ とすると、$f(c_k) - c_k = (T^* - 1)e_k + \tfrac12 f''e_k^2 + \cdots$ ($T^* = \partial f/\partial c|_{c^*}$、$f''$ は $f$ の 2 階微分の項。$f''e_k^2$ は $\sum_{jl}\partial^2 f/\partial c_j\partial c_l\,(e_k)_j(e_k)_l$ の略記)。Newton の更新 (12.1) は $T_k$ ($c_k$ での行列。$T_k = T^* + f''e_k + \cdots$) を使うので、

$$
\begin{aligned}
e_{k+1} &= e_k - (T_k - 1)^{-1}\bigl(f(c_k) - c_k\bigr)\\
&= (T_k - 1)^{-1}\bigl[(T_k - 1)e_k - (T^* - 1)e_k - \tfrac12f''e_k^2\bigr] &&(\text{$(T_k-1)^{-1}$ をくくる})\\
&= (T_k - 1)^{-1}\bigl[f''e_k^2 - \tfrac12f''e_k^2\bigr] &&(T_k - T^* = f''e_k)\\
&= \tfrac12(T^* - 1)^{-1}f''\,e_k^2 + O(e_k^3)
\end{aligned}
$$

**係数 $C$**: 残差は $r_k \simeq |(T^* - 1)e_k|$ なので、$r_{k+1} \simeq C\,r_k^2$ の形になる ($C$ は $\tfrac12 f''$ と $(T^* - 1)^{-1}$ の大きさで決まる、問題ごとの定数)。上の表から

| $k \to k+1$ | $r_k^2$ | $r_{k+1}/r_k^2$ |
|---|---|---|
| 1 → 2 | $7.9\times10^{-9}$ | 252 |
| 2 → 3 | $4.0\times10^{-12}$ | 250 |
| 3 → 4 | $1.0\times10^{-18}$ | 250 |

(各行は $r_{k+1}$ を $r_k$ の 2 乗で割った値。$2.0\times10^{-6}/7.9\times10^{-9} = 252$ など。[opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (E)。) つまり、単に「$10^{-5} \to 10^{-6}$」のような純粋な 2 乗ではなく、$\log_{10}r_{k+1} \simeq 2\log_{10}r_k + \log_{10}250$ ($\log_{10}250 = 2.4$) で下がる。$C \simeq 250$ は、この例の $f''$ (6 極 $K_2 = 30$、ずれた四極) と $(T^* - 1)^{-1}$ (作動点が整数に近く、§4) で決まる値で、**$C$ を $f''$ と $(T^* - 1)^{-1}$ から独立に計算して比べてはいない** (未確認)。5 回目の残差 $2.7\times10^{-19}$ は、$250\,r_4^2 \sim 10^{-29}$ より大きく、丸め誤差の床で止まっている。閉軌道は $x = 3.764\times10^{-4}$ m。1 回目の更新 (線形解) は $3.847\times10^{-4}$ m。

## 4. 周期解と分散

小さな $\delta$ のもとで、$f(c, \delta) \simeq M_4\,c + m_6\,\delta$ ($c = (x, p_x, y, p_y)$、$M_4$ は 1 周の左上 4×4、$m_6$ は第 6 列の上 4 成分)。この線形化は、$f$ を $(c, \delta) = (0, 0)$ のまわりで展開し (原点が設計軌道 $f(0, 0) = 0$ であるとする。偏向磁石の傾きやずれのあるリングでは $f(0,0) \ne 0$ で、その場合は (12.1) の閉軌道の解が出発点)、**$O(c^2)$、$O(c\,\delta)$、$O(\delta^2)$ の項を捨てた**もの。$c\,\delta$ の項は、$\delta$ による $M_4$ の変化 (色収差、§7) で、$\eta$ への影響は $\delta$ の 2 次 ($\eta$ の $\delta$ 依存、2 次分散) なので、$\delta$ の 1 次の係数 $\eta$ には要らない。不動点 $c = M_4c + m_6\delta$ を整理すると $(1 - M_4)\,c = m_6\,\delta$ なので、$\delta$ の 1 次の係数は

$$
\eta \equiv \frac{dc}{d\delta} = (1 - M_4)^{-1}m_6
\tag{12.2}
$$

$(1 - M_4)$ は、整数共鳴 ($M_4$ が固有値 1 を持つ) で特異になる。例のリング ($\nu_x = 1.9403$、$\nu_y = 1.0280$) では (12.2) で $(\eta_x, \eta_x', \eta_y, \eta_y') = (2.28993, 0.413678, 0, 0)$。追跡で $\delta = \pm10^{-5}$ の閉軌道 (Newton 法) の差分は $(2.289930, 0.4136776, 5\times10^{-16}, -9\times10^{-17})$、pysad の Twiss の `ex`、`epx` は $(2.289930, 0.413678)$ ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (b))。

**作動点の確認**: $\nu_x = 1.9403$ は整数 2 の近く ($0.06$ 離れている)、$\nu_y = 1.0280$ は整数 1 の近く ($0.03$) で、(12.2) の逆行列は大きくなる (このリングは数値の例で、実機の作動点ではない)。

## 5. x–y 結合の分解 (`qmdiag`)

### 5.1 ヘルプの定義

ヘルプの x-y-coupling: 物理座標 $(x, p_x, y, p_y)$ から、x–y が分離した座標 $(X, P_x, Y, P_y)$ への変換

$$
R = \begin{pmatrix}\mu\,I & J r^\top J\\ r & \mu\,I\end{pmatrix},\quad
r = \begin{pmatrix}\mathtt{R1} & \mathtt{R2}\\ \mathtt{R3} & \mathtt{R4}\end{pmatrix},\quad
\mu^2 + \det r = 1,\qquad
R^{-1} = \begin{pmatrix}\mu\,I & -J r^\top J\\ -r & \mu\,I\end{pmatrix}
\tag{12.3}
$$

(原文は TeX 版・平文版とも同じ。成分で書くと上の行は $(\mu, 0, -\mathtt{R4}, \mathtt{R2})$ と $(0, \mu, \mathtt{R3}, -\mathtt{R1})$。) 場所 1 から 2 への物理の行列 $T$ について、$R_2TR_1^{-1} = \mathrm{diag}(T_X, T_Y)$。$\det r \ge 1$ では $\mu^2 = 1 - \det r$ が実数にならないので、別の形 $R = \begin{pmatrix}J r^\top J & \mu I\\ \mu I & r\end{pmatrix}$ ($\mu^2 + \det r = 1$) を使う (ヘルプ: この形は $\det r \ge 0.625$ で使われ、このとき `DETR` は $a - \det r$、$a = 1.375$ を返す)。K64 の `qmdiag` は通常の形だけを作る (別の形を作る箇所は未確認)。

### 5.2 分解の条件

**$R$ がシンプレクティックであること**。$\bar r \equiv Jr^\top J$ として、$R = \begin{pmatrix}\mu I & \bar r\\ r & \mu I\end{pmatrix}$ ($\mu^2 + \det r = 1$)。$R\,J_4$ ($J_4 = \mathrm{diag}(J, J)$) に右から $R^\top = \begin{pmatrix}\mu I & r^\top\\ \bar r^\top & \mu I\end{pmatrix}$ をかける。まず 2×2 の恒等式 $rJr^\top = \det r\cdot J$ を成分で示す ($r = \begin{pmatrix}a&b\\c&d\end{pmatrix}$):

$$
rJ = \begin{pmatrix}a&b\\c&d\end{pmatrix}\begin{pmatrix}0&1\\-1&0\end{pmatrix} = \begin{pmatrix}-b&a\\-d&c\end{pmatrix}\quad(\text{$J$ を右から掛けると列が入れ替わり、第 1 列の符号が変わる}),
$$
$$
rJr^\top = \begin{pmatrix}-b&a\\-d&c\end{pmatrix}\begin{pmatrix}a&c\\b&d\end{pmatrix} = \begin{pmatrix}-ab + ab & -bc + ad\\ -ad + bc & -cd + cd\end{pmatrix} = (ad - bc)\begin{pmatrix}0&1\\-1&0\end{pmatrix} = \det r\cdot J
$$

(同じ恒等式 $XJX^\top = \det X\cdot J$ を、以後すべての 2×2 行列 $X$ に使う。$X^\top JX = \det X\cdot J$ も同様。) 次に、$\bar r = Jr^\top J$ の性質: $\det\bar r = \det J^2\det r = \det r$、$\bar r^\top = J^\top rJ^\top = JrJ$ ($J^\top = -J$、符号が 2 回で打ち消す)、$\bar rJ = Jr^\top JJ = -Jr^\top$ ($JJ = -I$)、$J\bar r^\top = JJrJ = -rJ$。したがって $\bar rJ\bar r^\top = \det\bar r\cdot J = \det r\cdot J$、$Jr^\top + \bar rJ = 0$、$rJ + J\bar r^\top = 0$。$RJ_4$ は $\begin{pmatrix}\mu J & \bar rJ\\ rJ & \mu J\end{pmatrix}$ なので、ブロックごとの積は

$$
RJ_4R^\top = \begin{pmatrix}\mu^2J + \bar rJ\bar r^\top & \mu(Jr^\top + \bar rJ)\\ \mu(rJ + J\bar r^\top) & rJr^\top + \mu^2J\end{pmatrix}
= \begin{pmatrix}(\mu^2 + \det r)J & 0\\ 0 & (\det r + \mu^2)J\end{pmatrix} = J_4
$$

(非対角ブロックは $Jr^\top - Jr^\top = 0$ と、その転置。) したがって $R$ は $\mu^2 + \det r = 1$ のときシンプレクティック。

1 周の行列を $M = \begin{pmatrix}P & N\\ Q & K\end{pmatrix}$ (各 2×2) とする。$R M R^{-1} = \mathrm{diag}(A, B)$、すなわち $R M = \mathrm{diag}(A, B)\,R$ を、(12.3) の $R$ ($\mu = c$、$\bar r = J r^\top J$) について 4 つのブロックに分けて書くと

$$
\begin{aligned}
c P + \bar r Q &= A c,\qquad & c N + \bar r K &= A \bar r,\\
r P + c Q &= B r,\qquad & r N + c K &= B c
\end{aligned}
\tag{12.4}
$$

$RM = \begin{pmatrix}cP + \bar rQ & cN + \bar rK\\ rP + cQ & rN + cK\end{pmatrix}$、$\mathrm{diag}(A, B)\,R = \begin{pmatrix}Ac & A\bar r\\ Br & Bc\end{pmatrix}$ の成分を比べたのが (12.4) である。左上と右下から $A = P + \bar r Q/c$、$B = K + rN/c$。$\hat r \equiv r/c$ とおいて下左の式から $A, B$ を消去すると $r$ に対する Riccati 型の方程式

$$
Q + \hat r P - K\hat r - \hat r N\hat r = 0
\tag{12.5}
$$

になる (右上の式は、シンプレクティックの条件でこれと同値。理由: $R$ がシンプレクティックなので $D \equiv RMR^{-1}$ もシンプレクティック。$D = \begin{pmatrix}A & Y\\ 0 & B\end{pmatrix}$ ((12.4) の下左の式は $D$ の左下ブロックが 0 という意味) の $DJ_4D^\top = J_4$ の右上ブロックは $YJB^\top = 0$、$B$ は可逆なので $Y = 0$。つまり左下ブロックが 0 になれば、右上ブロックも自動的に 0 になる。)幾何の意味: $R^{-1}$ の最初の 2 列 $\binom{c\,I}{-r}$ は、モード 1 の不変な 2 次元面 (物理座標で $(y, p_y) = -\hat r\,(x, p_x)$ という形) を張る。(12.5) はこの面が $M$ で不変であるという条件。

### 5.3 閉じた解 (`qmdiag` の通常の場合)

K64 の `qmdiag.f` は、(12.5) の解を次の形で与える (Edwards–Teng の式):

$$
\hat r = \kappa\,(Q - \bar N),\qquad
\tau \equiv \mathrm{tr}(P - K),\quad
\mathrm{trab2} \equiv \tau^2 + 4\det(Q - \bar N),
$$
$$
c^2 = \frac12\left(1 + \frac{|\tau|}{\sqrt{\mathrm{trab2}}}\right),\qquad
\kappa = \frac{-1}{c^2\sqrt{\mathrm{trab2}}}\,\mathrm{sgn}(\tau) = -\frac{2}{\tau + \mathrm{sgn}(\tau)\sqrt{\mathrm{trab2}}}
\tag{12.6}
$$

(`qmdiag.f` の変数: `tmij` は 1 周の行列の $M_{ij}$、`trpq` $= M_{11} + M_{22} - M_{33} - M_{44} = \mathrm{tr}P - \mathrm{tr}K = \tau$ (コメントの `TrP-TrQ`。Fortran のコメントは右下のブロックを「Q」と呼ぶ。本章の $Q$ は左下ブロックなので、本章では $K$ と書く)、`trab2` $= \tau^2 + 4\det(Q - \bar N)$ ((12.6)。$\mathrm{tr}A - \mathrm{tr}B$ の 2 乗、§5.3 の手順 6)。Fortran の `trab2 = 4*((tm31+tm24)*(tm42+tm13)+(tm14-tm32)*(tm41-tm23)) + trpq**2` は、$Q - \bar N = \begin{pmatrix}M_{31} + M_{24} & M_{32} - M_{14}\\ M_{41} - M_{23} & M_{42} + M_{13}\end{pmatrix}$ の行列式が $(M_{31} + M_{24})(M_{42} + M_{13}) + (M_{14} - M_{32})(M_{41} - M_{23})$ に等しいことから、上の式と同じ。$\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$ で、$\ge 0$ が安定の条件。) 成分で書いた Fortran:

$$
\mathtt{aa} = \mp\frac{1}{c\sqrt{\mathrm{trab2}}},\quad
\begin{pmatrix}r_1 & r_2\\ r_3 & r_4\end{pmatrix} = \mathtt{aa}\begin{pmatrix}M_{31} + M_{24} & M_{32} - M_{14}\\ M_{41} - M_{23} & M_{42} + M_{13}\end{pmatrix}
\tag{12.7}
$$

(上の符号 $-$ は $\tau \ge 0$、$+$ は $\tau < 0$。行列は $Q - \bar N$ で、$\bar N = \begin{pmatrix}-M_{24} & M_{14}\\ M_{23} & -M_{13}\end{pmatrix}$。ここの `aa` は `qmdiag.f` のローカル変数で、$c\kappa$ ($= -\mathrm{sgn}(\tau)/(c\sqrt{\mathrm{trab2}})$、§5.3 の手順 5)。$r = c\hat r = c\kappa(Q - \bar N)$ の係数を 1 つにまとめたもの。後の `tffs.f` の `aa` とは別の変数 (§8)。)

**確認** ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (d)): 結合した安定なシンプレクティック行列 ($|M^\top J M - J| = 4\times10^{-16}$、$\nu = (0.31, 0.17)$ の対角形をランダムなシンプレクティック行列で結合させたもの) で、(12.6)(12.7) の $r$、$c$ を作ると
- $c^2 + \det r - 1 = 0$ ($c = 0.995759$、$\det r = 0.008465$)。
- $R$ はシンプレクティック ($|R^\top J_4R - J_4| = 1.1\times10^{-16}$、$J_4 = \mathrm{diag}(J, J)$)、ヘルプの $R^{-1}$ ($r$ の符号を逆) との積は $1.1\times10^{-16}$ で単位行列。
- $R M R^{-1}$ の非対角ブロックの最大 $8.9\times10^{-17}$ (対角ブロックの最大は 3.0)。
- (12.5) の残差 $5.3\times10^{-17}$、$\hat r = \kappa(Q - \bar N)$ との差 $1.4\times10^{-17}$。
- pysad の `optics.periodic_r` の $R$ と $0$ (差なし)、$A$ とは $2.2\times10^{-16}$。

**(12.6) が (12.5) を満たすことの代数的な証明**。記号: 2×2 行列 $A$ について $\bar A = JA^\top J = -\mathrm{adj}\,A$ (余因子行列の符号反転)。性質は (F1) $A\bar A = \bar AA = -\det A\cdot I$、(F2) $\overline{AB} = -\bar B\bar A$、(F3) $\bar{\bar A} = A$、(F4) $\mathrm{tr}\bar A = -\mathrm{tr}A$、(F5) $A - \bar A = (\mathrm{tr}A)\,I$ ($\mathrm{adj}\,A = (\mathrm{tr}A)I - A$ から)。証明を 1 行ずつ書く ($A = \begin{pmatrix}a&b\\c&d\end{pmatrix}$、$J^{-1} = J^\top = -J$、$JJ^\top = I$。[opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (C) で記号計算でも 0 を確認):
- 成分: $\bar A = JA^\top J$。$A^\top J = \begin{pmatrix}a&c\\b&d\end{pmatrix}\begin{pmatrix}0&1\\-1&0\end{pmatrix} = \begin{pmatrix}-c&a\\-d&b\end{pmatrix}$、$J\,(A^\top J) = \begin{pmatrix}0&1\\-1&0\end{pmatrix}\begin{pmatrix}-c&a\\-d&b\end{pmatrix} = \begin{pmatrix}-d&b\\c&-a\end{pmatrix}$。これは $-\mathrm{adj}\,A$ (§2 の式)。
- (F1): $A\bar A = (AJA^\top)J = (\det A\cdot J)J = -\det A\cdot I$ (前項の $AJA^\top = \det A\cdot J$ と $JJ = -I$)。成分では $\begin{pmatrix}a&b\\c&d\end{pmatrix}\begin{pmatrix}-d&b\\c&-a\end{pmatrix} = \begin{pmatrix}-ad + bc & ab - ab\\ -cd + dc & cb - da\end{pmatrix} = -\det A\cdot I$。$\bar AA$ も同じ ($\bar A = -\mathrm{adj}A$ は左逆元でもある)。
- (F2): $\overline{AB} = J(AB)^\top J = JB^\top A^\top J = (JB^\top J)(J^{-1}A^\top J)$ (間に $JJ^{-1} = I$ を入れた)。$JB^\top J = \bar B$、$J^{-1}A^\top J = -JA^\top J = -\bar A$ ($J^{-1} = -J$) なので $\overline{AB} = -\bar B\bar A$。
- (F3): $\bar{\bar A} = J(JA^\top J)^\top J = J\,J^\top A\,J^\top J = (JJ^\top)A(J^\top J) = A$。
- (F4): 成分から $\mathrm{tr}\bar A = -d - a = -\mathrm{tr}A$。
- (F5): $A - \bar A = \begin{pmatrix}a + d & 0\\0 & a + d\end{pmatrix} = (\mathrm{tr}A)I$ (成分から。$b - b = 0$、$c - c = 0$)。

1. **シンプレクティックな $M$ の逆行列**: $M^{-1} = -J_4M^\top J_4 = -\begin{pmatrix}\bar P & \bar Q\\ \bar N & \bar K\end{pmatrix}$。$MM^{-1} = M^{-1}M = 1$ の各ブロック (成分で $\bar{\cdot} = -\mathrm{adj}$) から、
   - $Q\bar P + K\bar N = 0$ (B)、$\bar NP + \bar KQ = 0$ (D)
   - $P\bar P + N\bar N = -1$ より $\det P + \det N = 1$、$\bar PP + \bar QQ = -1$ より $\det P + \det Q = 1$。したがって $\det N = \det Q$ (E)。
2. $X \equiv Q - \bar N$、$\tau \equiv \mathrm{tr}(P - K)$ とおく。**恒等式 (i)**: $XP - KX = \tau Q$。証明: $XP - KX = QP - \bar NP - KQ + K\bar N$。(D) から $\bar NP = -\bar KQ$、(B) から $K\bar N = -Q\bar P$。代入して $= QP + \bar KQ - KQ - Q\bar P = Q(P - \bar P) + (\bar K - K)Q = (\mathrm{tr}P)\,Q - (\mathrm{tr}K)\,Q$ ((F5) で $P - \bar P = \mathrm{tr}P\cdot I$、$\bar K - K = -\mathrm{tr}K\cdot I$)。
3. **恒等式 (ii)**: $XNX = \det X\cdot Q$。証明: $NX = NQ - N\bar N = NQ + \det N$ ((F1))。一方 $-\bar XQ$ は、$\bar X = \bar Q - \overline{\bar N} = \bar Q - N$ ((F3)) なので $-\bar XQ = -\bar QQ + NQ = \det Q + NQ$ ((F1))。(E) から $NX = -\bar XQ$。よって $XNX = X(NX) = -X\bar XQ = \det X\cdot Q$ ((F1) で $X\bar X = -\det X$)。
4. $\hat r = \kappa X$ を (12.5) の左辺に代入する。(i)、(ii) を使うと

$$
Q + \hat rP - K\hat r - \hat rN\hat r = Q + \kappa(XP - KX) - \kappa^2XNX = Q\,\bigl(1 + \tau\kappa - \det X\,\kappa^2\bigr)
$$

したがって (12.5) は**スカラーの 2 次方程式** $\det X\,\kappa^2 - \tau\kappa - 1 = 0$ と同値になる。解は $\kappa = [\tau \pm \sqrt{\tau^2 + 4\det X}]/(2\det X)$。ここで $\tau^2 + 4\det X = \mathrm{trab2}$ ($\det X = \det(Q - \bar N)$ は上で出た $\mathrm{trab2}$ の定義の第 2 項に等しい)。分母を有理化すると ($(\tau \pm\sqrt{\mathrm{trab2}})(\tau \mp\sqrt{\mathrm{trab2}}) = -4\det X$)

$$
\kappa = -\frac{2}{\tau \mp \sqrt{\mathrm{trab2}}}
$$

($\tau$ と同符号の項を足す方、$\kappa = -2/(\tau + \mathrm{sgn}(\tau)\sqrt{\mathrm{trab2}})$ を選ぶと桁落ちしない。これが (12.6) の 2 番目の等号。) 2 つの根は 2 つのモードの不変面に対応する。
5. **$c$ の式**: $c^2 + \det r = 1$、$r = c\hat r$ から $c^2(1 + \det\hat r) = 1$、$\det\hat r = \kappa^2\det X$。2 次方程式から $\kappa^2\det X = \tau\kappa + 1$ なので $c^2 = 1/(2 + \tau\kappa)$。$\kappa = -2/(\tau + s\sqrt{\mathrm{trab2}})$ ($s = \mathrm{sgn}\tau$) を入れると $2 + \tau\kappa = 2s\sqrt{\mathrm{trab2}}/(\tau + s\sqrt{\mathrm{trab2}})$、したがって $c^2 = \dfrac{\tau + s\sqrt{\mathrm{trab2}}}{2s\sqrt{\mathrm{trab2}}} = \dfrac12\Bigl(1 + \dfrac{|\tau|}{\sqrt{\mathrm{trab2}}}\Bigr)$ ((12.6) の $c^2$)。$\kappa = -s/(c^2\sqrt{\mathrm{trab2}})$ ((12.6) の 1 番目の等号) もこれから出る ($c^2\sqrt{\mathrm{trab2}} = (\sqrt{\mathrm{trab2}} + |\tau|)/2$)。$r = c\kappa X$ で $c\kappa = -s/(c\sqrt{\mathrm{trab2}})$ が (12.7) の `aa`。
6. **$\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$**: $A = P + \bar{\hat r}Q$、$B = K + \hat rN$ ($c$ が消える) なので $\mathrm{tr}A - \mathrm{tr}B = \tau + \kappa[\mathrm{tr}(\bar XQ) - \mathrm{tr}(XN)]$。$\bar X = \bar Q - N$ ((F3)) から $\mathrm{tr}(\bar XQ) = \mathrm{tr}(\bar QQ) - \mathrm{tr}(NQ) = -2\det Q - \mathrm{tr}(NQ)$、$\mathrm{tr}(XN) = \mathrm{tr}(QN) - \mathrm{tr}(\bar NN) = \mathrm{tr}(QN) + 2\det N$。したがって $\mathrm{tr}(\bar XQ) - \mathrm{tr}(XN) = -2[\det Q + \det N + \mathrm{tr}(NQ)]$。一方、2×2 の $\det(A + B) = \det A + \det B + \mathrm{tr}(\mathrm{adj}(A)B)$ から $\det X = \det Q + \det\bar N + \mathrm{tr}(\bar Q\bar N)$ ($\mathrm{adj}\,Q = -\bar Q$、$\det\bar N = \det N$)。$\mathrm{tr}(\bar Q\bar N) = -\mathrm{tr}\,\overline{NQ} = \mathrm{tr}(NQ)$ ((F2)(F4)) なので $\det X = \det Q + \det N + \mathrm{tr}(NQ)$。よって $\mathrm{tr}(\bar XQ) - \mathrm{tr}(XN) = -2\det X$ で、$\mathrm{tr}A - \mathrm{tr}B = \tau - 2\kappa\det X$。2 次方程式から $\kappa\det X = \tau + 1/\kappa$ なので $\mathrm{tr}A - \mathrm{tr}B = -\tau - 2/\kappa$。$\kappa = -2/(\tau + s\sqrt{\mathrm{trab2}})$ から $-2/\kappa = \tau + s\sqrt{\mathrm{trab2}}$ で、$\mathrm{tr}A - \mathrm{tr}B = s\sqrt{\mathrm{trab2}}$。すなわち $\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$ で、$\ge 0$ が安定の条件。

**数値確認** ([opt_03_qmdiag_etwiss.out.txt](outputs/opt_03_qmdiag_etwiss.out.txt) (a)): 200 個のランダムなシンプレクティック行列で、(i) $|XP - KX - \tau Q| \le 5 \times 10^{-16}$、(ii) $|XNX - \det X\,Q| \le 1.1 \times 10^{-16}$、$|\det N - \det Q| \le 1.9 \times 10^{-16}$、$|\det P + \det Q - 1| \le 4.4 \times 10^{-16}$、(12.6) の $\hat r$ での Riccati の残差 $\le 10^{-14}$、$|\det X\kappa^2 - \tau\kappa - 1| \le 4.4 \times 10^{-16}$、$c^2$ の式の差 $\le 1.8 \times 10^{-15}$。

**場合分け** (`qmdiag.f`): $\mathrm{trab2}$ の値で 3 つに分かれる ($\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$ で、Fortran では `trab2 ≥ 1e-14` が通常、`|trab2| ≤ 1e-14` が特別、`trab2 < -1e-14` が不安定)。
- **通常** ($\mathrm{trab2} \ge 10^{-14}$): (12.6)(12.7)。
- **特別** ($|\mathrm{trab2}| \le 10^{-14}$。2 つのモードのトレースが等しい。等しいチューンの結合、すなわち差共鳴の真上): 2 次方程式の判別式が 0 で重根 $\kappa = \tau/(2\det X) = -2/\tau$ ($\tau^2 = -4\det X$)、$c^2 = 1/(2 + \tau\kappa) = 1/0$ と発散するため (12.6) は使えない。K64 は $c = 1$ (または $1/\sqrt2$) とし、$r$ を次の量の場合分けで与える (すべて `qmdiag.f` のローカル変数。値は $M$ の成分 `tmij` $= M_{ij}$ から決まり、閾値で分岐する)。
  - `detn` $= \det N = M_{13}M_{24} - M_{14}M_{23}$ ($N$ は $M$ の右上 2×2 ブロック。閾値 $10^{-14}$)。$N$ は x–y の結合の強さを表す。
  - `trqk` $= M_{33} - M_{44}$、`trqj` $= M_{43} - M_{34}$、`trql` $= M_{34} + M_{43}$ (閾値 $10^{-7}$)。$M$ の右下ブロック $K = \begin{pmatrix}M_{33}&M_{34}\\M_{43}&M_{44}\end{pmatrix}$ について、$\mathrm{tr}(\sigma_zK) = M_{33} - M_{44}$ ($\sigma_z = \mathrm{diag}(1, -1)$)、$\mathrm{tr}(JK) = M_{43} - M_{34}$、$\mathrm{tr}(\sigma_xK) = M_{34} + M_{43}$ ($\sigma_x = \begin{pmatrix}0&1\\1&0\end{pmatrix}$) で、$K$ のトレースなし部分の 3 成分。Fortran のコメントは右下のブロックを「Q」、2×2 の基底 $\sigma_z$、$J$、$\sigma_x$ を「K」「J」「L」と呼んで `Tr(QK)`・`Tr(QJ)`・`Tr(QL)` と書く (本章の $Q$ は左下ブロックなので、本章ではこの記法を使わない)。`trqkj2` $= (\mathtt{trqk} + \mathtt{trqj})(\mathtt{trqk} - \mathtt{trqj})$ (閾値 $10^{-14}$)。
  - `aaa`、`bbb`、`ccc`: 小分岐ごとに決まる 3 つの数。$\det N \ne 0$ では $X_f = \begin{pmatrix}\mathtt{aaa} & \mathtt{bbb} + \mathtt{ccc}\\ \mathtt{ccc} - \mathtt{bbb} & -\mathtt{aaa}\end{pmatrix}$ (`xx11`〜`xx22`) を作り、$r = X_fN^{-1}$ とする ($N^{-1} = \mathrm{adj}N/\det N$ なので、Fortran の `r1 = (xx11*tm24 - xx12*tm23)/detn` などと同じ。[opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (G))。$\det N = 0$ では $Y_f = \begin{pmatrix}\mathtt{aaa} + \mathtt{bbb} & \mathtt{ccc}\\ \mathtt{ccc} & \mathtt{aaa} - \mathtt{bbb}\end{pmatrix}$ (`yy11`〜`yy22`) と $r = Y_f\begin{pmatrix}-M_{23} & M_{13}\\ -M_{24} & M_{14}\end{pmatrix}$。
  - $\det N \ne 0$ の 4 つの小分岐 ($(c;\ \mathtt{aaa}, \mathtt{bbb}, \mathtt{ccc})$): (i) $|\mathtt{trqk}| \ge 10^{-7}$ かつ $|\mathtt{trqkj2}| \ge 10^{-14}$: $c = 1$、$\mathtt{bbb} = -\mathtt{detn}/(\mathtt{trqj} + \mathtt{trqk})$、$\mathtt{ccc} = 0$、$\mathtt{aaa} = [\mathtt{detn}(1 - 2c^2)/c - \mathtt{bbb}\,\mathtt{trqj}]/\mathtt{trqk}$。(ii) $|\mathtt{trqk}| \ge 10^{-7}$ かつ $|\mathtt{trqkj2}| < 10^{-14}$: $c = 1$、$\mathtt{bbb} = -\mathtt{detn}/(2\,\mathtt{trqj})$、$\mathtt{ccc} = 0$、`aaa` は (i) と同じ式。(iii) $|\mathtt{trqk}| < 10^{-7}$、$|\mathtt{trqj}| \ge 10^{-7}$: $c = 1$、$\mathtt{aaa} = \mathtt{detn}/\mathtt{trqj}$、$\mathtt{bbb} = -\mathtt{detn}/\mathtt{trqj}$、$\mathtt{ccc} = 0$。(iv) $\mathtt{trqk} \approx \mathtt{trqj} \approx 0$: $|\mathtt{trql}| \ge 10^{-7}$ なら $c = 1$、$\mathtt{aaa} = 0$、$\mathtt{bbb} = \mathtt{detn}/\mathtt{trql}$、$\mathtt{ccc} = -\mathtt{detn}/\mathtt{trql}$。3 つとも 0 なら $c = 1/\sqrt2$ で、$\mathtt{detn} > 0$ は $\mathtt{bbb} = \sqrt{\mathtt{detn}/2}$ (他は 0)、$\mathtt{detn} < 0$ は $\mathtt{aaa} = \sqrt{-\mathtt{detn}/2}$ (他は 0)。
  - $\det N = 0$ の小分岐: $c = 1$ で、$|\mathtt{trqk}| \ge 10^{-7}$ なら $\mathtt{ccc} = 1/\mathtt{trqk}$ (他は 0)、そうでなく $|\mathtt{trqj}| \ge 10^{-7}$ なら $\mathtt{aaa} = -1/\mathtt{trqj}$、そうでなく $|\mathtt{trql}| \ge 10^{-7}$ なら $\mathtt{bbb} = -1/\mathtt{trql}$。3 つとも 0 の場合は $P$ の成分から: $p_3 = M_{21}$、$p_1 = M_{11} - (M_{11} + M_{22})/2$ として、$|p_3| \ge 10^{-7}$ なら $r = (0, M_{24}/p_3, 0, -M_{23}/p_3)$、そうでなく $|p_1| \ge 10^{-7}$ なら $r = (M_{24}/p_1, 0, -M_{23}/p_1, 0)$、どちらも 0 なら $c = 1/\sqrt2$、$r = (1, 0, 0, c)$。

    これらの式は `qmdiag.f` から読み取った表で、**各式がなぜ $RMR^{-1}$ のブロック対角化を与えるかの導出はしていない** (下の数値確認と §11)。
- **不安定** ($\mathrm{trab2} < -10^{-14}$): $\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$ が負は、$\mathrm{tr}A$、$\mathrm{tr}B$ が複素共役 (実 2×2 ブロックにならない) ことで、**和共鳴** ($\nu_x + \nu_y$ が整数の近く) の不安定。`stab = .false.`、"Sum resonance"。

**数値確認** (同 (b)): Fortran の分岐を numpy に写し、(1) $\mu_x = \mu_y = 0.9$ rad の対角形をシンプレクティックな変換で結合させた縮退した行列 ($\mathrm{trab2} = 4 \times 10^{-33}$) では、分岐は「$\det N \ne 0$、$|\mathtt{trqk}| \ge 10^{-7}$、$|\mathtt{trqkj2}| \ge 10^{-14}$」(上の (i)) で、$c = 1$、$R M R^{-1}$ の非対角ブロック $4.5 \times 10^{-17}$、$R$ のシンプレクティック性 $7 \times 10^{-20}$、$c^2 + \det r - 1 = 0$ と、特別な式が通常の式と同じ性質の $R$ を与える。(2) 結合のない等しいチューンの行列では $\det N = 0$ の分岐で $r = 0$、非対角ブロック 0。(3) $\mu_x = 0.9$、$\mu_y = -0.9$ に $x$–$y$ 結合 ($H$ に $k\,xy$) を足すと、$k = 0.1$、$0.3$ で $\mathrm{trab2} = -2.5 \times 10^{-2}$、$-0.22$、固有値の絶対値の最大が $1.05$、$1.16$ で、`qmdiag` は和共鳴の分岐を返す。**特別な場合の各小分岐の式の導出はしていない** (上の (1) の分岐が $R M R^{-1}$ をブロック対角にすることの数値確認まで。他の小分岐は縮退の特殊な行列でのみ通るため、一般の行列で網羅的に試してはいない)。

## 6. Twiss パラメータ

### 6.1 2×2 の安定な行列

結合を分解した後の 2×2 の 1 周行列 $A$ は、Courant–Snyder の形

$$
A = \begin{pmatrix}\cos\mu + \alpha\sin\mu & \beta\sin\mu\\ -\gamma\sin\mu & \cos\mu - \alpha\sin\mu\end{pmatrix},\qquad \gamma = \frac{1 + \alpha^2}{\beta}
\tag{12.8}
$$

に書ける。$\det A = (\cos\mu + \alpha\sin\mu)(\cos\mu - \alpha\sin\mu) + \beta\gamma\sin^2\mu = \cos^2\mu - \alpha^2\sin^2\mu + (1 + \alpha^2)\sin^2\mu = 1$ で、$\det A = 1$ の 2×2 行列は ($\mathrm{tr}A/2 = \cos\mu$、$\alpha$、$\beta$ の 3 つの数で) この形に書ける。逆に $\cos\mu = \mathrm{tr}A/2$、$\beta = A_{12}/\sin\mu$、$\alpha = (A_{11} - A_{22})/(2\sin\mu)$ (pysad: `a = (A11 - c)/s`)。$\sin\mu$ の符号は $A_{12}$ の符号で決める ($\beta > 0$)。$|\mathrm{tr}A| < 2$ が安定。

ヘルプの $B = \begin{pmatrix}1/\sqrt\beta & 0\\ \alpha/\sqrt\beta & \sqrt\beta\end{pmatrix}$ と回転 $P(\psi) = \begin{pmatrix}\cos\psi & \sin\psi\\ -\sin\psi & \cos\psi\end{pmatrix}$ について、$B A B^{-1} = P(\mu)$。$\det B = (1/\sqrt\beta)\sqrt\beta - 0 = 1$ で、2×2 では $\det = 1$ がシンプレクティックの条件なので $B$ はシンプレクティック。逆行列は ($\det B = 1$ の 2×2 の逆行列 $= $ 対角を入れ替え非対角の符号を変える) $B^{-1} = \begin{pmatrix}\sqrt\beta & 0\\ -\alpha/\sqrt\beta & 1/\sqrt\beta\end{pmatrix}$。**行列の積を成分ごとに**書く。$s \equiv \sin\mu$、$c_\mu \equiv \cos\mu$、$\gamma\beta = 1 + \alpha^2$ を使う。

1. $BA$ の第 1 行 ($B$ の第 1 行 $(1/\sqrt\beta,\ 0)$ を $A$ の各列に掛ける): $\dfrac{c_\mu + \alpha s}{\sqrt\beta}$ ($(1,1)$ 成分)、$\dfrac{\beta s}{\sqrt\beta} = \sqrt\beta\,s$ ($(1,2)$ 成分)。
2. $BA$ の第 2 行 ($B$ の第 2 行 $(\alpha/\sqrt\beta,\ \sqrt\beta)$):
   - $(2,1)$: $\dfrac{\alpha(c_\mu + \alpha s)}{\sqrt\beta} - \sqrt\beta\,\gamma s = \dfrac{\alpha c_\mu + \alpha^2 s - (1 + \alpha^2)s}{\sqrt\beta} = \dfrac{\alpha c_\mu - s}{\sqrt\beta}$ ($\beta\gamma = 1 + \alpha^2$)。
   - $(2,2)$: $\dfrac{\alpha\beta s}{\sqrt\beta} + \sqrt\beta(c_\mu - \alpha s) = \sqrt\beta\,\alpha s + \sqrt\beta\,c_\mu - \sqrt\beta\,\alpha s = \sqrt\beta\,c_\mu$。

   $$BA = \begin{pmatrix}\dfrac{c_\mu + \alpha s}{\sqrt\beta} & \sqrt\beta\,s\\[2mm] \dfrac{\alpha c_\mu - s}{\sqrt\beta} & \sqrt\beta\,c_\mu\end{pmatrix}$$
3. $(BA)B^{-1}$ ($BA$ の各行に $B^{-1}$ の各列 $(\sqrt\beta,\ -\alpha/\sqrt\beta)^\top$、$(0,\ 1/\sqrt\beta)^\top$ を掛ける):
   - $(1,1)$: $\dfrac{c_\mu + \alpha s}{\sqrt\beta}\sqrt\beta + \sqrt\beta\,s\left(-\dfrac{\alpha}{\sqrt\beta}\right) = c_\mu + \alpha s - \alpha s = c_\mu$。
   - $(1,2)$: $0 + \sqrt\beta\,s\cdot\dfrac{1}{\sqrt\beta} = s$。
   - $(2,1)$: $\dfrac{\alpha c_\mu - s}{\sqrt\beta}\sqrt\beta + \sqrt\beta\,c_\mu\left(-\dfrac{\alpha}{\sqrt\beta}\right) = \alpha c_\mu - s - \alpha c_\mu = -s$。
   - $(2,2)$: $\sqrt\beta\,c_\mu\cdot\dfrac{1}{\sqrt\beta} = c_\mu$。

   よって $BAB^{-1} = \begin{pmatrix}\cos\mu & \sin\mu\\ -\sin\mu & \cos\mu\end{pmatrix} = P(\mu)$ ([opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (A) で、成分ごとの $BA$ と $P(\mu)$ の sympy の差が 0。[opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (c) でも確認)。

したがって $B$ は、軌道を正規化座標 $(X, P_X)$ ( 位相空間で円) に直す。

### 6.2 伝搬

場所 1 から 2 への 2×2 の行列 $T$ ($\det T = 1$) について、$A_2 = TA_1T^{-1}$ (同じ位相の 1 周行列)。$T^{-1} = \begin{pmatrix}t_{22} & -t_{12}\\ -t_{21} & t_{11}\end{pmatrix}$ なので、$A_1 = (a_{ij})$ として $(TA_1)_{1\cdot} = (t_{11}a_{11} + t_{12}a_{21},\ t_{11}a_{12} + t_{12}a_{22})$ に $T^{-1}$ の第 2 列 $(-t_{12}, t_{11})^\top$ をかけると

$$
(A_2)_{12} = t_{11}^2a_{12} - t_{11}t_{12}(a_{11} - a_{22}) - t_{12}^2a_{21}
$$

(12.8) の $a_{12} = \beta_1\sin\mu$、$a_{11} - a_{22} = 2\alpha_1\sin\mu$、$a_{21} = -\gamma_1\sin\mu$ を入れると $(A_2)_{12} = \sin\mu\,(t_{11}^2\beta_1 - 2t_{11}t_{12}\alpha_1 + t_{12}^2\gamma_1)$。これが $\beta_2\sin\mu$ ((12.8) の (1,2) 成分) なので、$\mu$ が共通で消えて

$$
\beta_2 = t_{11}^2\beta_1 - 2t_{11}t_{12}\alpha_1 + t_{12}^2\gamma_1,\qquad
\alpha_2 = -t_{11}t_{21}\beta_1 + (t_{11}t_{22} + t_{12}t_{21})\alpha_1 - t_{12}t_{22}\gamma_1
\tag{12.9}
$$

**$\alpha_2$ の導出** ($A_2$ は (12.8) の形で $\alpha_2$ を含み、$(A_2)_{11} - (A_2)_{22} = 2\alpha_2\sin\mu$)。$(TA_1)$ の第 1 行 $= (t_{11}a_{11} + t_{12}a_{21},\ t_{11}a_{12} + t_{12}a_{22})$、第 2 行 $= (t_{21}a_{11} + t_{22}a_{21},\ t_{21}a_{12} + t_{22}a_{22})$。$T^{-1}$ の第 1 列は $(t_{22},\ -t_{21})^\top$、第 2 列は $(-t_{12},\ t_{11})^\top$ なので

$$
\begin{aligned}
(A_2)_{11} &= t_{22}(t_{11}a_{11} + t_{12}a_{21}) - t_{21}(t_{11}a_{12} + t_{12}a_{22}),\\
(A_2)_{22} &= -t_{12}(t_{21}a_{11} + t_{22}a_{21}) + t_{11}(t_{21}a_{12} + t_{22}a_{22})
\end{aligned}
$$

差を取ると、$a_{11}$、$a_{22}$、$a_{12}$、$a_{21}$ の係数は

$$
\begin{aligned}
a_{11}&: \ t_{11}t_{22} + t_{12}t_{21}, & a_{22}&: \ -t_{12}t_{21} - t_{11}t_{22},\\
a_{12}&: \ -t_{11}t_{21} - t_{11}t_{21} = -2t_{11}t_{21}, & a_{21}&: \ t_{22}t_{12} + t_{12}t_{22} = 2t_{12}t_{22}
\end{aligned}
$$

($a_{12}$ は $(A_2)_{11}$ の $-t_{21}t_{11}$ と $-(A_2)_{22}$ の $-t_{11}t_{21}$、$a_{21}$ は $t_{22}t_{12}$ と $+t_{12}t_{22}$)。したがって

$$
(A_2)_{11} - (A_2)_{22} = (t_{11}t_{22} + t_{12}t_{21})(a_{11} - a_{22}) - 2t_{11}t_{21}\,a_{12} + 2t_{12}t_{22}\,a_{21}
$$

(12.8) から $a_{11} - a_{22} = 2\alpha_1\sin\mu$、$a_{12} = \beta_1\sin\mu$、$a_{21} = -\gamma_1\sin\mu$ を入れる:

$$
(A_2)_{11} - (A_2)_{22} = \sin\mu\,\bigl[2(t_{11}t_{22} + t_{12}t_{21})\alpha_1 - 2t_{11}t_{21}\beta_1 - 2t_{12}t_{22}\gamma_1\bigr] = 2\alpha_2\sin\mu
$$

$2\sin\mu$ で割ると (12.9) の $\alpha_2$ の 3 項 ($-t_{11}t_{21}\beta_1$、$(t_{11}t_{22} + t_{12}t_{21})\alpha_1$、$-t_{12}t_{22}\gamma_1$) になる。($(A_2)_{11} + (A_2)_{22} = (t_{11}t_{22} - t_{12}t_{21})(a_{11} + a_{22}) = \mathrm{tr}A_1$ ($\det T = 1$) なので、$\cos\mu$ は場所によらない。) sympy で $a_{ij}$ の記号のまま各成分と $\alpha_2$、$\beta_2$ を展開し、(12.9) との差が 0: [opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (B) ([opt_03_qmdiag_etwiss.out.txt](outputs/opt_03_qmdiag_etwiss.out.txt) (d) でも $\det T = 1$ のもとで確認)。 位相の進みは、正規化座標 ($B$ で変換した座標。(12.8) の $BAB^{-1} = P(\mu)$) での回転角で、$B_2TB_1^{-1}$ の $(1,2)$ 成分と $(1,1)$ 成分の比: $B_1^{-1} = \begin{pmatrix}\sqrt{\beta_1} & 0\\ -\alpha_1/\sqrt{\beta_1} & 1/\sqrt{\beta_1}\end{pmatrix}$ から $(B_2TB_1^{-1})_{11} = (t_{11}\sqrt{\beta_1} - t_{12}\alpha_1/\sqrt{\beta_1})/\sqrt{\beta_2}$、$(B_2TB_1^{-1})_{12} = t_{12}/(\sqrt{\beta_1\beta_2})$。回転 $P(\Delta\psi)$ の $(1,2)/(1,1) = \tan\Delta\psi$ なので $\tan\Delta\psi = t_{12}/(\beta_1t_{11} - \alpha_1t_{12})$ (sympy でも確認)。ドリフト ($t_{11} = 1$、$t_{12} = L$、$t_{21} = 0$、$t_{22} = 1$) では $\beta_2 = \beta_1 - 2\alpha_1L + \gamma_1L^2$、$\alpha_2 = \alpha_1 - \gamma_1L$。pysad の Twiss の行 5 つのドリフトと、$9\times10^{-16}$ の差で一致。

結合があるときは、場所 1 と 2 のそれぞれの $R$ で $R_2TR_1^{-1} = \mathrm{diag}(T_X, T_Y)$ とし ((12.3) の下)、$T_X$、$T_Y$ に (12.9) を使う。

## 7. 色収差

$\xi = d\nu/d\delta$。1 周の行列から: $\cos\mu(\delta) = \mathrm{tr}A(\delta)/2$ の $\delta$ 微分と $\xi = -\dfrac{1}{2\pi\sin\mu}\dfrac{d\cos\mu}{d\delta}$ ($\nu = \mu/2\pi$、$d\cos\mu = -\sin\mu\,d\mu$)。ただし $\delta \ne 0$ では閉軌道も動く (分散)。四極の寄与は、薄いレンズ ($k = K_1L$、行列 $\begin{pmatrix}1&0\\-k&1\end{pmatrix}$) で、レンズの直後からの 1 周行列 $M_0$ に $\delta k$ だけ強さが変わると $M = M_0\begin{pmatrix}1&0\\-\delta k&1\end{pmatrix}$、$\mathrm{tr}M - \mathrm{tr}M_0 = -\delta k\,m_{12}$ (sympy で確認)。$\cos\mu = \mathrm{tr}M/2$ から $-\sin\mu\,\delta\mu = -\delta k\,m_{12}/2$、$m_{12} = \beta\sin\mu$ ((12.8)) より $\delta\mu = \tfrac12\beta\,\delta k$。運動量の偏差 $\delta$ で $K_1 \to K_1/(1+\delta)$ なので $\delta k = -k\delta$、$\delta\nu = \delta\mu/2\pi$ から $\xi_{\rm quad} = -\dfrac{1}{4\pi}\sum_i K_{1i}L_i\,\beta_i$ (垂直面は $K_1$ の符号が逆)。

**確認** ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (e2)): 偏向なしの 4 極だけのリング ($\nu = (1.1995, 1.0280)$) で、$d\nu/d\delta$ は
| 方法 | $d\nu_x/d\delta$ | $d\nu_y/d\delta$ |
|---|---|---|
| 薄いレンズの式 ($\beta$ は要素の両端の平均) | $-1.1664$ | $-1.1273$ |
| 追跡の Jacobian (差分 $\delta = \pm10^{-4}$) | $-1.2144$ | $-1.1799$ |
| `lat.chromaticity()` | $-1.2144$ | $-1.1799$ |

追跡と `lat.chromaticity()` は一致し、薄いレンズの式は水平で 4.0%、垂直で 4.5% 低い ($|{-1.1664}|/|{-1.2144}| = 0.960$、$|{-1.1273}|/|{-1.1799}| = 0.955$)。この差は四極の厚みによる。次に導く。

**厚い四極**。薄いレンズの導出は、$M_0$ が任意の「レンズの直後からの 1 周行列」でよいので、四極を $N$ 個の薄いスライス (強さ $k_i = K\Delta s$、$K$ は長さあたりの強さ) に分けて各スライスに適用でき、1 次の摂動は足し合わせられる。$N \to \infty$ で $\delta\mu = \tfrac12\sum_i\beta_i\,\delta k_i \to \tfrac12\int\beta(s)\,\delta K\,ds$、$\delta K = -K\delta$ から

$$
\xi = -\frac1{4\pi}\oint K(s)\,\beta(s)\,ds
\tag{12.11}
$$

薄いレンズの式は、この積分を「$\beta$ は四極の中で一定で、両端の値の平均」と近似したもの。(ここの $K$ は pysad の `el.K1` ($=$ SAD の入力の積分強度 `K1` を長さ $L$ で割ったもの [1/m²])。表の QF は SAD の `K1`=0.35 [1/m]、$L$=0.3 m で $K = 1.1667$ [1/m²]。垂直面は $K \to -K$。)

**四極の中の $\beta(s)$**。$\beta' = -2\alpha$、$\alpha' = K\beta - \gamma$ ($\gamma = (1 + \alpha^2)/\beta$) から $\beta'' = -2\alpha' = 2\gamma - 2K\beta$。入口 $(\beta_0, \alpha_0, \gamma_0)$ から、四極の行列 $\begin{pmatrix}C&S\\-KS&C\end{pmatrix}$、$C = \cos\sqrt K s$、$S = \sin(\sqrt K s)/\sqrt K$ ($K > 0$) で伝搬した $\beta$ は ((12.9) と同じ形。$t_{11} = C$、$t_{12} = S$)

$$
\beta(s) = C^2\beta_0 - 2CS\,\alpha_0 + S^2\gamma_0
$$

(上の微分方程式を満たすことを sympy で確認。(D)。) これを $0 \le s \le L$ で積分する。$\int_0^LC^2 ds = \tfrac L2 + \tfrac{\sin2\phi}{4\sqrt K}$、$\int_0^LS^2ds = \tfrac1K\bigl(\tfrac L2 - \tfrac{\sin2\phi}{4\sqrt K}\bigr)$、$\int_0^LCS\,ds = \tfrac{\sin^2\phi}{2K}$ ($\phi \equiv \sqrt K L$、四極の位相 [rad]。$\phi^2 = KL^2$ が薄さの目安) なので、四極内の平均 $\bar\beta \equiv \frac1L\int_0^L\beta\,ds$ は

$$
\bar\beta = \beta_0\Bigl(\frac12 + \frac{\sin2\phi}{4\phi}\Bigr) - \alpha_0L\Bigl(\frac{\sin\phi}{\phi}\Bigr)^2 + \gamma_0L^2\,\frac{1}{\phi^2}\Bigl(\frac12 - \frac{\sin2\phi}{4\phi}\Bigr)
\tag{12.12}
$$

($K < 0$ は $\phi \to i|\phi|$ で $\sin \to \sinh$。) 一方、薄いレンズの式の $(\beta_0 + \beta_1)/2$ ($\beta_1 = \beta(L)$)。$L$ で展開すると ($\sin2\phi/4\phi = \tfrac12 - \tfrac{\phi^2}{3} + \tfrac{2\phi^4}{15} - \cdots$ を使う。展開は sympy でも確認):

$$
\begin{aligned}
\bar\beta &= \beta_0 - \alpha_0L + \tfrac13(\gamma_0 - K\beta_0)L^2 + \tfrac13K\alpha_0L^3 + O(L^4),\\
\tfrac12(\beta_0 + \beta_1) &= \beta_0 - \alpha_0L + \tfrac12(\gamma_0 - K\beta_0)L^2 + \tfrac23K\alpha_0L^3 + O(L^4)
\end{aligned}
$$

(第 2 式は $\beta_1 = \beta_0 - 2\alpha_0L + (\gamma_0 - K\beta_0)L^2 + \tfrac43K\alpha_0L^3 + \cdots$ を $\beta'(0) = -2\alpha_0$、$\beta''(0) = 2(\gamma_0 - K\beta_0)$、$\beta'''(0) = 8K\alpha_0$ から作ったもの。)

$$
\bar\beta - \tfrac12(\beta_0 + \beta_1) = -\tfrac16(\gamma_0 - K\beta_0)\,L^2 - \tfrac13K\alpha_0\,L^3 + O(L^4)
\tag{12.13}
$$

**結論 (厚みの補正の大きさ)**: $O(L)$ の項 ($-\alpha_0L$) は両者で同じで消える。最初に残るのは $L^2$ の項で、$K\beta_0L^2/6 = \phi^2\beta_0/6$ を主要項として (γ0 の項は $\gamma_0/(K\beta_0)$ だけ小さい)、相対的な大きさ $\dfrac{\bar\beta - (\beta_0 + \beta_1)/2}{\beta_0} \simeq \dfrac{\phi^2}{6}\Bigl(1 - \dfrac{\gamma_0}{K\beta_0}\Bigr)$。小さいパラメータは $L/\beta$ でなく $\phi^2 = KL^2$ (四極の位相の 2 乗) で、1 次ではなく 2 次。例のリングの四極 ($L = 0.3$ m、$|K| \simeq 1.1$、$\phi^2 \simeq 0.1$) では 1 台あたり $\pm1.7$〜1.8% の補正で、符号は $K\beta$ の符号 (QF の水平 $\bar\beta$ は端の平均より大きく、QD の水平は小さい) による。

| 面・四極 | $K$ [1/m²] | $\phi^2$ | $\beta_0$ [m] | $\alpha_0$ | 端の平均 [m] | $\bar\beta$ (12.12) [m] | 比 $\bar\beta/$端の平均 | $1 + \frac{\phi^2}{6}(1 - \frac{\gamma_0}{K\beta_0})$ |
|---|---|---|---|---|---|---|---|---|
| x、QF | 1.1667 | 0.105 | 8.762 | −1.530 | 8.762 | 8.916 | 1.0176 | 1.0169 |
| x、QD | −1.1000 | −0.099 | 3.741 | 0.653 | 3.741 | 3.676 | 0.9827 | 0.9820 |
| y、QF | −1.1667 | −0.105 | 4.389 | 0.796 | 4.389 | 4.310 | 0.9820 | 0.9812 |
| y、QD | 1.1000 | 0.099 | 10.021 | −1.653 | 10.021 | 10.187 | 1.0166 | 1.0159 |

([opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (D)。$\beta_0$ は入口の値。このリングの四極は入口と出口で $\beta$ が等しい (中央に対称) ので、端の平均 $= \beta_0$。$\bar\beta$ は数値積分と (12.12) が 5 桁で一致。)

**4〜4.5% の理由**: $\sum K L\,\beta$ では QF の項が正、QD の項が負で、差を取っている。水平: 薄いレンズの式 QF $= 0.35 \times 8.7618 = 3.0666$、QD $= -0.33\times3.7408 = -1.2345$、和 $1.8321$ (1 セルあたり)。厚い四極は QF $= 0.35\times8.9159 = 3.1206$、QD $= -0.33\times3.6759 = -1.2130$、和 $1.9076$。QF は大きく、QD は小さい絶対値に補正される (どちらも和を増やす向き) うえ、補正が 1.7% でも 2 つの大きな項の差なので、和では $1.9076/1.8321 = 1.041$ ($4.1\%$) に拡大される。この 8 セルの和は薄いレンズの式 $-1.1664$ に対し、(12.11) を $\bar\beta$ で積分した値は $-1.2144$ (水平)、垂直は $-1.1273 \to -1.1799$ で、**追跡の Jacobian $(-1.2144,\ -1.1799)$ と 4 桁で一致する** (opt_04 (D))。したがって薄いレンズの式との差は、(12.11) の積分を薄い近似にした分で、全部説明できる。

**偏向磁石のあるリング** (例のリング): 3 通りで $d\nu/d\delta$ を比べた ([opt_02_chromaticity.out.txt](outputs/opt_02_chromaticity.out.txt))。

| 方法 | $d\nu_x/d\delta$ | $d\nu_y/d\delta$ |
|---|---|---|
| 厳密なハミルトニアン $H = -(1 + hx)\sqrt{p^2 - p_x^2 - p_y^2} + h(x + hx^2/2) + K(x^2 - y^2)/2$ の写像 (DOP853、フリンジなし。$\delta = \pm10^{-4}$ の閉軌道 (Newton 法) と 1 周のヤコビアンの位相の差分) | $-0.77973$ | $-0.83429$ |
| pysad の追跡 + Newton 法 + ヤコビアン (DISFRIN の指定によらず同じ) | $-0.77974$ | $-0.87868$ |
| `lat.chromaticity()` | $-1.45631$ | $-1.17990$ |

水平は、厳密な $H$ の写像と追跡が $10^{-5}$ で一致し ($-0.7797$)、`lat.chromaticity()` ($-1.456$) だけが大きく違う。垂直は、追跡 ($-0.879$) と厳密な $H$ の写像 ($-0.834$) が 0.045 違い、`lat.chromaticity()` ($-1.180$) は四極だけのリングの値 ($-1.1799$) と同じで、**偏向磁石の寄与を含まない**。単独の偏向磁石 (DISFRIN = 1) の追跡と厳密な $H$ の写像を比べると、$x$ のみ・$\delta$ のみでは差 0 ($10^{-16}$)、$y \ne 0$ で $p_y$ に差が出て ($y = 1$ mm、$\delta = 10^{-3}$ で $1.25 \times 10^{-7}$。$h^2L\delta y = 1.29 \times 10^{-7}$ と同じ大きさ)、これが垂直の 0.045 の差の原因と読めるが、BEND の追跡の $y$ 依存の項 (硬い端の Maxwell フリンジか本体か) の内訳は未確認。

したがって、**偏向磁石のあるリングで合わないのは pysad の `Lattice.chromaticity()` (要素の `transfer_matrix(δ)` を使う pysad 独自の計算。docstring に「SAD の直訳ではない」とある) で、K64 の問題でも追跡の問題でもない** (水平: 独立な厳密な $H$ の写像との一致で確定。垂直: 追跡が厳密な $H$ の写像と 0.045 違う理由は上のとおり未確認だが、`lat.chromaticity()` はそれらとも 0.3〜0.7 違う)。SAD の `CHROMATICITY` コマンドは四極・6 極だけの閉形式 (BEND を数えない) で、この例では TOTAL が $-30.6$ (SAD の位相 [rad]、$2\pi$ で割って $-4.87$) と別の量を出す (`tfchro` の式。何の量かはこの章では整理しなかった)。

## 8. 6 次元の拡張 Twiss

ヘルプの extended-Twiss-parameters と definitions: 6 次元の正規モードの行列 $V$ ($\bm U = V\bm u$、$\bm U = (X, P_x, Y, P_y, Z, P_z)$、$\bm u = (x, p_x, y, p_y, z, \delta)$) を

$$
V = P\,B\,R_6\,H,\qquad
R_6 = \begin{pmatrix}R & 0\\ 0 & I\end{pmatrix},\quad
H = \begin{pmatrix}
\left(1 - \dfrac{\det H_x}{1+a}\right)I & \dfrac{H_xJ_2H_y^\top J_2}{1+a} & -H_x\\[2mm]
\dfrac{H_yJ_2H_x^\top J_2}{1+a} & \left(1 - \dfrac{\det H_y}{1+a}\right)I & -H_y\\[2mm]
-J_2H_x^\top J_2 & -J_2H_y^\top J_2 & a\,I
\end{pmatrix}
\tag{12.10}
$$

($a^2 + \det H_x + \det H_y = 1$、$b^2 + \det R = 1$ ($b$ は (12.3) の $\mu$)、$P$、$B$ は各モードの (12.8) の回転と $B$。) 分散は $\begin{pmatrix}\mathtt{ZX}&\mathtt{EX}\\ \mathtt{ZPX}&\mathtt{EPX}\\ \mathtt{ZY}&\mathtt{EY}\\ \mathtt{ZPY}&\mathtt{EPY}\end{pmatrix} = R\begin{pmatrix}H_x\\ H_y\end{pmatrix}$ で定義される (EX・EPX・EY・EPY が分散、ZX 以下は $z$ モードとの結合 (synchro-beta))。同期ベータ結合がなければ、通常の Twiss と一致すると、ヘルプは書く。

**$H$ がシンプレクティックであることの代数的な導出**。記号: $u \equiv 1 + a$、$d_x \equiv \det H_x$、$d_y \equiv \det H_y$ ($a^2 + d_x + d_y = 1$)、$p_x \equiv 1 - d_x/u$、$p_y \equiv 1 - d_y/u$、2×2 のバー $\bar X = JX^\top J$ (§2。$\bar X^\top = JXJ$、$\det\bar X = \det X$、§5.2)。(12.10) の各ブロックは $H_{11} = p_xI$、$H_{12} = H_x\bar H_y/u$、$H_{13} = -H_x$、$H_{21} = H_y\bar H_x/u$、$H_{22} = p_yI$、$H_{23} = -H_y$、$H_{31} = -\bar H_x$、$H_{32} = -\bar H_y$、$H_{33} = aI$ ($H_xJ_2H_y^\top J_2 = H_x\bar H_y$ による)。$J_6 = \mathrm{diag}(J, J, J)$ に対して $(H^\top J_6H)_{ij} = \sum_kH_{ki}^\top JH_{kj}$ なので、$(H^\top J_6H)_{ij} = \delta_{ij}J$ の各ブロックを示す。2×2 の恒等式 $X^\top JX = \det X\cdot J$ (§5.2 と同じ) を使う。
- **(3,3) ブロック** ($k = 1, 2, 3$ の 3 項): $H_{13}^\top JH_{13} + H_{23}^\top JH_{23} + H_{33}^\top JH_{33} = H_x^\top JH_x + H_y^\top JH_y + a^2J = (d_x + d_y + a^2)J = J$。
- **(1,1) ブロック**: $H_{11}^\top JH_{11} = p_x^2J$。$H_{21}^\top JH_{21} = \dfrac{1}{u^2}(H_y\bar H_x)^\top J(H_y\bar H_x) = \dfrac{\det(H_y\bar H_x)}{u^2}J = \dfrac{d_xd_y}{u^2}J$。$H_{31}^\top JH_{31} = \bar H_x^\top J\bar H_x = d_xJ$。和の $J$ の係数に $u^2$ をかけると
  $$(u - d_x)^2 + d_xd_y + d_xu^2 = u^2 + d_x\,(d_x + d_y + u^2 - 2u) = u^2 + d_x\,(d_x + d_y + a^2 - 1) = u^2$$
  ($-2ud_x + d_x^2 + d_xd_y + d_xu^2$ を $d_x$ でくくり、$u^2 - 2u = (1 + a)^2 - 2(1 + a) = a^2 - 1$ を使った。最後は $a^2 + d_x + d_y = 1$)。したがって係数は 1。
- **(1,3) ブロック**: $H_{11}^\top JH_{13} = -p_x\,JH_x$。$H_{21}^\top JH_{23} = -\dfrac1u\bar H_x^\top H_y^\top JH_y$ で、$H_y^\top JH_y = d_yJ$、$\bar H_x^\top = JH_xJ$ なので $= -\dfrac{d_y}{u}(JH_xJ)J = +\dfrac{d_y}{u}JH_x$ ($JJ = -I$)。$H_{31}^\top JH_{33} = a(-\bar H_x^\top)J = -a(JH_xJ)J = +a\,JH_x$。和は $JH_x\,[-p_x + d_y/u + a] = JH_x\,\dfrac{-u + d_x + d_y + au}{u} = JH_x\,\dfrac{-(1 + a) + (1 - a^2) + a(1 + a)}{u} = 0$。
- **(1,2) ブロック**: $H_{11}^\top JH_{12} = \dfrac{p_x}{u}JH_x\bar H_y$。$H_{21}^\top JH_{22} = \dfrac{p_y}{u}(H_y\bar H_x)^\top J = \dfrac{p_y}{u}\bar H_x^\top H_y^\top J$。$H_y^\top J = J^{-1}(JH_y^\top J) = -J\bar H_y$ を使うと $= -\dfrac{p_y}{u}(JH_xJ)J\bar H_y = +\dfrac{p_y}{u}JH_x\bar H_y$。$H_{31}^\top JH_{32} = (-\bar H_x)^\top J(-\bar H_y) = (JH_xJ)J\bar H_y = -JH_x\bar H_y$。和は $JH_x\bar H_y\bigl[(p_x + p_y)/u - 1\bigr]$、$p_x + p_y = 2 - (d_x + d_y)/u = 2 - (1 - a^2)/(1 + a) = 2 - (1 - a) = u$ なので 0。
- **(2,2)、(2,3) ブロック**: 上の (1,1)、(1,3) で $H_x \leftrightarrow H_y$、ブロックの第 1 と第 2 を入れ替えたものと同じ ($H$ は、この入れ替えで $H_{12} \leftrightarrow H_{21}$ が互いに移る形)。(2,1)、(3,1)、(3,2) ブロックは $H^\top J_6H$ が反対称 ($J_6^\top = -J_6$) であることから従う。

以上で $H^\top J_6H = J_6$ ($a^2 + d_x + d_y = 1$ のもとで代数的に)。各中間式の sympy の確認 (中間式は差が 0) と、6 ブロックすべての記号計算 ($H_y$ の 1 成分を $a^2 + d_x + d_y = 1$ から消去)、有理数 30 組 ($a = 3/5$) での全 36 成分の厳密な確認: [opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) (F)。浮動小数点では、$a^2 + \det H_x + \det H_y = 1$ を満たす 5 組のランダムな $H_x, H_y$ で最大 $|H^\top J_6H - J_6| = 2.2\times10^{-16}$ ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (f))。**この導出は (12.10) の $H$ について**で、ヘルプの $V = PBR_6H$ 全体がシンプレクティックなことは、$P$、$B$、$R_6$ (§5.2) がそれぞれシンプレクティックなことから従う (未確認: ヘルプにある $V$ の他の主張 ($V$ が 6 次元の固有ベクトルとして 1 周行列を正規形にする、など) は導いていない)。

**K64 との対応** (`tffs.f` の `etwiss2ri`、`tfetwiss`)。`tecalc` は $V$ を (12.10) のように組み立てず、1 周の 6×6 行列の固有ベクトルを数値で求め (`teigenc`・`tnorm`・`tsymp`)、その行列 `ri` から `tfetwiss` が $(\beta, \alpha, \nu, R, H_x, H_y, \ldots)$ を分解して取り出す ((12.10) の逆)。逆に `etwiss2ri` は、パラメータから 6×6 行列を作る (位相 $P$ なし)。
- `etwiss2ri` (`tffs.f:2569`) の変数: `hx11`〜`hy22` は $H_x$、$H_y$ の成分 ($h_{x,11}$ など。分散 `ex` などから `hx11 = amu*zx - r2*zpy + r4*zy`、`hx12 = amu*ex - r2*epy + r4*ey` のように作る。$R\,[H_x; H_y]$ の逆)、`dethx` $= h_{x11}h_{x22} - h_{x21}h_{x12} = \det H_x$ (`dethy` も同様)、`azz` $= \sqrt{1 - \mathtt{dethx} - \mathtt{dethy}}$ ($= a$。(12.10) の $a$。$a^2 + \det H_x + \det H_y = 1$ から)、`aa` $= 1/(1 + \mathtt{azz})$ (`qmdiag.f` の `aa` とは別の変数。$1/(1 + a)$ を毎回割らずに済ませる)、`amux` $= 1 - \mathtt{dethx}\cdot\mathtt{aa}$ (`amuy` も同様)、`h(4,6)` は $H$ の上 4 行 (x・y モードの行)、`ria(6,6)` は関数 `etwiss2ri` の戻り値の 6×6 行列 (パラメータから作った行列。位相 $P$ なし。名前の `ri` の由来は未確認)。`ria(5:6, :)` が $z$ モードの行。`amu` $= \sqrt{1 - \mathtt{detr}}$ は (12.3) の $\mu$。
- `etwiss2ri` の `h(1,1)` $=$ `amux` $= 1 - \det H_x/(1 + a)$、`h(1,3)`・`h(1,4)` $= (H_xJ_2H_y^\top J_2)_{1\cdot}/(1 + a)$ (例: `h(1,3)` $= (h_{x12}h_{y21} - h_{x11}h_{y22})/(1 + a)$)、`h(1,5)`・`h(1,6)` $= -(h_{x11}, h_{x12})$、`ria(5,1:2)` $= (h_{x22}, -h_{x12})$ $= (-J_2H_x^\top J_2)_{1\cdot}$、`ria(5,5)` $= a$ と、(12.10) の $H$ の成分が 1 つずつ対応する。`br` は $B$ と $R$ の積 ($\mu/\sqrt{\beta_x}$ など)。
- `tfetwiss` (`tffs.f:2413`) の変数: 入力 `r(6,6)` は 1 周の行列の固有ベクトルの行列 (`tecalc` の結果)、`az` $= \sqrt{r_{55}r_{66} - r_{65}r_{56}}$ ($z$ モードの 2×2 ブロックの行列式の平方根。ブロックを行列式 1 にそろえる規格化の係数)、`uz11`・`uz12`・`uz21`・`uz22` $= r_{55}/\mathtt{az}$、$r_{56}/\mathtt{az}$、$r_{65}/\mathtt{az}$、$r_{66}/\mathtt{az}$ (行列式 1 に規格化した $z$ ブロック $U_z$。`az` は (12.10) の $a$ に当たる: $z$ モードのブロックは (12.10) で $a\times$(行列式 1 の 2×2) の形で、行列式が $a^2$ になるため。$H$ の $z$ 列の成分は、$x$、$y$ の行の $z$ 成分に $U_z^{-1}$ 相当の組合せをかけて取り出す)。続いて `f` $= 1/(1 + \mathtt{az})$、`ax` $= 1 - (h_{x11}h_{x22} - h_{x21}h_{x12})\,\mathtt{f}$ ($= $ `etwiss2ri` の `amux`、`ay` は `amuy`) と `hi` (上の `etwiss2ri` の `h` に当たる行列) を作る。
- `tfetwiss` は、$H_x$ の成分を `r(5..6, 1..4)` の行 (5・6 行目は $z$ モードの行) から `hx11 = r(6,2)*uz11 - r(5,2)*uz21` のように取り出す。分散は `twiss(mfitex) = axy*hx12 - r22*hy12 + r12*hy22` で、これは $R = \begin{pmatrix}\mu I & Jr^\top J\\ r & \mu I\end{pmatrix}$ の第 1 行と $[H_x; H_y]$ の第 2 列の積 $\mu h_{x12} + (-R_4h_{y12} + R_2h_{y22})$ である。すなわち **EX、EPX、EY、EPY は $R\,[H_x; H_y]$ の第 2 列、ZX、ZPX、ZY、ZPY は第 1 列**。

**確認** ([opt_03_qmdiag_etwiss.out.txt](outputs/opt_03_qmdiag_etwiss.out.txt) (c)): Fortran の `etwiss2ri` を numpy に写して作った 6×6 行列は、(12.10) の $B\,R_6\,H$ (ヘルプの式から作ったもの) と 5 組のランダムなパラメータで $1.1 \times 10^{-16}$ で一致し、シンプレクティック ($1.1 \times 10^{-16}$)。`tfetwiss` を写して `etwiss2ri` の出力に作用させると、$\beta_{x,y}$、$\alpha_{x,y}$、$R$ が $10^{-15}$ で戻り、EX、EPX、EY、EPY、ZX、ZPX、ZY、ZPY が $R\,[H_x; H_y]$ の第 2 列・第 1 列と $1.4 \times 10^{-17}$ で一致する。したがって、(12.10) と K64 の `etwiss2ri` は同じ行列で、$H_x$ の成分と分散の対応は上のとおり。`tecalc` が固有ベクトルからこの形の行列を得る手続き (正規化 `tnorm`、シンプレクティック化 `tsymp`) は、ここでは確認していない。

## 9. Fortran の対応

| 内容 | Fortran | 式 | pysad |
|---|---|---|---|
| 閉軌道の Newton 法 | `tcod.f` `tcod` (`1 loop`) | (12.1)、収束 `codw`、`conv0` | `emit6d._tcod_radcod`、`rdttrack.closed_orbit` (4 次元) |
| 一周の行列 | `tturne` (`trans`) | | `lattice._transport_6x6` |
| 結合の分解 | `qmdiag.f` `qmdiag` | (12.6)(12.7)、特別な場合 | `optics.periodic_r` |
| セルの周期解 | `qcell.f` `qcell1`・`qcell6d` | | `calc6d.calc6d` |
| Twiss の伝搬 | `qtwiss.f`、`qdtwis` など | (12.9) | `optics.transfer_matrix`、`lattice.twiss` |
| 色収差 | `tfchro.f` | §7 | `ffscmd.chro_cmd`、`lattice.chromaticity` |

## 10. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `opt_02_chromaticity.py` | 偏向磁石のあるリングの色収差: 厳密な $H$ の写像、pysad の追跡、`lat.chromaticity()` の比較、単独の偏向磁石の追跡と厳密な $H$ の差 | [opt_02_chromaticity.out.txt](outputs/opt_02_chromaticity.out.txt) |
| `opt_03_qmdiag_etwiss.py` | (12.6) の代数的な証明に使う恒等式、`qmdiag` の特別な場合、`etwiss2ri`・`tfetwiss` と (12.10)、(12.9) と位相の進み・薄いレンズの色収差 (sympy) | [opt_03_qmdiag_etwiss.out.txt](outputs/opt_03_qmdiag_etwiss.out.txt) |
| `opt_04_gaps.py` | 導出の途中式の確認: (A) $BAB^{-1}$ の行列の積 (§6.1)、(B) $A_2 = TA_1T^{-1}$ の成分と $\alpha_2$ (§6.2)、(C) 2×2 の恒等式と (F1)〜(F5) (§5.2・§5.3)、(D) 厚い四極の $\bar\beta$ (12.12) と端の平均との差 (12.13)・例のリングの数値 (§7)、(E) Newton 法の収束係数 $C \simeq 250$ (§3.2)、(F) $H$ のシンプレクティック性の手計算の中間式・全ブロックの記号計算と有理数での確認 (§8)、(G) `qmdiag` の $\det N 
e 0$、$= 0$ の分岐の式が $X_fN^{-1}$、$Y_f(\cdots)$ であること (§5.3) | [opt_04_gaps.out.txt](outputs/opt_04_gaps.out.txt) |
| `opt_01_linear_optics.py` | (a) 閉軌道の Newton 法 (2 次収束)、(b) 分散 (12.2) と追跡の差分・pysad の Twiss、(c) $BAB^{-1} = P(\mu)$ と (12.9) (sympy)・ドリフトでの伝搬、(d) `qmdiag` の式 (12.6)(12.7) の確認 (ブロック対角、Riccati 式)、(e) 色収差 (偏向あり)、(e2) 四極だけのリングの色収差、(f) ヘルプの $H$ のシンプレクティック性 | [opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/opt_01_linear_optics.py` (十数秒)。

## 11. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV。8 セルの FODO リング (QF・QD = 0.3 m、K1 = ±0.35・−0.33 [1/m]、偏向 1.2 m × 16 台の扇形、周長 41.6 m、$\nu = (1.9403, 1.0280)$)。放射・RF なし。偏向の端面なし、フリンジなし。結合のテストは数値の行列 ($\nu = (0.31, 0.17)$ の対角形をシンプレクティックな変換で結合)。

**未確認**:
- `qmdiag` の特別な場合 ($\mathrm{trab2} \approx 0$) の各小分岐の式の導出、$\det r \ge 0.625$ の別の形を K64 のどこで作るか。(12.6) が (12.5) を満たすことは代数的に証明した (§5.3)。
- 偏向磁石のあるリングの色収差: SAD の CHROMX (CALC で $\delta$ をずらした $\nu$ の差分) との比較。垂直の追跡と厳密な $H$ の写像の 0.045 の差の内訳 (§7)。`lat.chromaticity()` が偏向磁石の寄与を取り違えている点の具体的な原因 (コード側の調査) は未実施。
- `tecalc` が固有ベクトルから $V$ を得る手続き (`tnorm`・`tsymp`)。$V$ と $H$ の対応は `etwiss2ri`・`tfetwiss` で確認した (§8)。
- 閉軌道: RADCOD・RADTAPER の分岐 ($dz$ を `trf0` に移す部分、`dptaper`) の導出と確認。Newton 法の収束係数 $C \simeq 250$ (§3.2) を、$f''$ と $(T^* - 1)^{-1}$ から独立に計算して比べてはいない (2 次収束の式と、残差の列から読んだ値の一致のみ)。
- 色収差の厚い四極 (§7): (12.11)〜(12.13) と数値は偏向のない四極だけのリングの 1 例 (3 GeV、$L = 0.3$ m、$\phi^2 \simeq 0.1$) での確認。偏向磁石の端面・分散に由来する寄与 (偏向磁石のあるリング) は (12.11) の範囲外で、`lat.chromaticity()` との差はそこ (§7、未確認)。
- 拡張 Twiss (§8): (12.10) の $H$ のシンプレクティック性は代数的に示したが、ヘルプの $V = PBR_6H$ が固有ベクトルとして 1 周行列を正規形にするという主張、`tecalc` の手続きは導いていない。
- 結合の強いとき ($\det r$ が大きい) の整数共鳴・和共鳴の近くの挙動。作動点の整数・半整数共鳴は、(12.2) の逆行列で確認する (例のリングは整数共鳴から 0.03〜0.06 離れており、分散が大きい)。

## 12. この章で使う Fortran の変数名

[01_variables_and_conventions.md](01_variables_and_conventions.md) §8 の規約による一覧 (定義の確認は K64 の `tcod.f`・`qmdiag.f`・`tffs.f` の実際のコード)。行番号は K64 版ソース (`SAD-master/src/`) のもの。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `codi`、`codf` | $c$ (試行の閉軌道)、$f(c)$ (1 周追跡後の座標) | `tcod.f` `tcod` | `emit6d._tcod_radcod` の内部 |
| `dcod1` | $c - f(c)$ (12.1) の右辺 | `tcod.f:49` | 同 (内部) |
| `trs` | $T - 1$ (12.1) の左辺の行列 | `tcod.f:49` | 同 (内部) |
| `fact` | Newton の更新の減衰係数 (1 で全量、残差が増えたら半分、`fmin` $= 10^{-4}$ で打ち切り) | `tcod.f:83` | 同 (内部) |
| `codw`、`conv0` | 収束判定の重み (6 成分)、閾値 $10^{-4}$ (§2) | `tcod.f:3,5` | `emit6d._tcod_radcod` |
| `tmij` | 1 周の 6×6 行列の $M_{ij}$ (`qmdiag` の引数。4×4 で使う) | `qmdiag.f` の引数 | `optics.periodic_r` |
| `trpq` | $M_{11} + M_{22} - M_{33} - M_{44} = \tau = \mathrm{tr}P - \mathrm{tr}K$ | `qmdiag.f:34` | (内部) |
| `trab2` | $\tau^2 + 4\det(Q - \bar N) = (\mathrm{tr}A - \mathrm{tr}B)^2$、符号で通常・特別・不安定を分ける | `qmdiag.f:36` | (内部) |
| `aa` (`qmdiag.f`) | $c\kappa = \mp1/(c\sqrt{\mathrm{trab2}})$ ((12.7)) | `qmdiag.f:68` | (内部) |
| `detn` | $\det N = M_{13}M_{24} - M_{14}M_{23}$ (特別な場合の分岐の判定) | `qmdiag.f:84` | (内部) |
| `trqk`、`trqj`、`trql` | $M_{33} - M_{44}$、$M_{43} - M_{34}$、$M_{34} + M_{43}$ ($K$ のトレースなしの 3 成分。§5.3) | `qmdiag.f:89,92,132` (`detn = 0` では 169〜171) | (内部) |
| `trqkj2` | $(\mathtt{trqk} + \mathtt{trqj})(\mathtt{trqk} - \mathtt{trqj})$ | `qmdiag.f:93` | (内部) |
| `aaa`、`bbb`、`ccc` | 特別な場合の $r$ を作る 3 つの係数 ($X_f$、$Y_f$ の成分。§5.3) | `qmdiag.f:120` ほか | (内部) |
| `r1`〜`r4`、`c` | $r$ の成分 (R1〜R4)、$\mu$ | `qmdiag.f` の出力 | `optics.periodic_r` の `R`、`a` |
| `amu` | $\sqrt{1 - \det r}$ ((12.3) の $\mu$) | `tffs.f:2588` (`etwiss2ri`) | (内部) |
| `dethx`、`dethy` | $\det H_x$、$\det H_y$ | `tffs.f:2626` | (内部) |
| `azz` | $a = \sqrt{1 - \det H_x - \det H_y}$ ((12.10)) | `tffs.f:2628` | (内部) |
| `aa` (`etwiss2ri`) | $1/(1 + a)$ | `tffs.f:2629` | (内部) |
| `amux`、`amuy` | $1 - \det H_x/(1 + a)$、$1 - \det H_y/(1 + a)$ ((12.10) の対角) | `tffs.f:2630` | (内部) |
| `ria` | `etwiss2ri` の戻り値の 6×6 行列 (パラメータから作る。位相なし) | `tffs.f:2569` | (内部) |
| `az` | $z$ モードのブロックの行列式の平方根 ($= a$) | `tffs.f:2430` (`tfetwiss`) | (内部) |
| `uz11`〜`uz22` | 行列式 1 に規格化した $z$ ブロック | `tffs.f:2431` | (内部) |
| `hx11`〜`hy22` | $H_x$、$H_y$ の成分 | `tffs.f:2435` (`tfetwiss`)、`2618` (`etwiss2ri`) | (内部) |

「(内部)」は pysad 側で対応する名前を確認していない (K64 の変数は pysad では別の形で書かれている、または pysad に無い)。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| x-y-coupling の $R$、$R^{-1}$、$\mu^2 + \det r = 1$ | **原文で確認** (TeX 版・平文版とも同じ。成分表示も TeX 版に有る)。K64 の `qmdiag` と対応 (§5)。 |
| $\det r \ge 1$ の別の形と DETR | **原文で確認**: TeX 版は「$\det r \ge 1$」、「$\det r \ge 0.625$ で別の形」、「DETR は $a - \det r$、$a = 1.375$」と、同じ文の中で $1$ と $0.625$ の 2 つの閾値が出てくる。平文版も同じ。通常の形の $\mu$ が実数であるには $\det r \le 1$ でよいので、0.625 は DETR の見かけ上の切り替えの閾値とみられる (未確認)。 |
| extended-Twiss の定義 | **原文で確認**: (12.10)、TeX 版と平文版は同じ式 (記号のみ違う)。 |
| ヘルプに無いもの | (1) 閉軌道の Newton 法 (§3)。(2) 周期解と分散の式 (12.2)。(3) `qmdiag` の閉じた式 (12.6)。(4) Twiss の伝搬 (12.9)。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

1. 偏向磁石のあるリングで、pysad の `Lattice.chromaticity()` と追跡 (Jacobian) の $d\nu/d\delta$ が、水平 $-1.456$ 対 $-0.780$、垂直 $-1.180$ 対 $-0.879$ と大きく違う (§7)。四極だけのリングでは一致。**pysad 独自の `Lattice.chromaticity()` が偏向磁石の寄与を取り違えている (確定)**: 水平は、追跡 ($-0.7797$) が、独立な厳密なハミルトニアンの写像 ($-0.77973$。DOP853、フリンジなし、閉軌道を分散で動かす) と $10^{-5}$ で一致する。K64 の問題ではない。垂直の追跡値 ($-0.879$) は、厳密な $H$ の写像 ($-0.834$) と 0.045 違い、単独の BEND の追跡が $y \ne 0$ で厳密な $H$ の写像と $p_y$ で $h^2L\delta y$ の程度違うことと整合する (内訳は未確認)。
