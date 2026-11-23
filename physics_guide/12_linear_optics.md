# 12. 行列から光学へ: 閉軌道、周期解、x–y 結合、Twiss

**この章で分かること**: 要素ごとの写像と行列 (02〜11 章) から、リングの光学関数 (閉軌道、分散、x–y 結合、Twiss パラメータ、拡張 Twiss) が何をどう計算して出るか。(1) 閉軌道の Newton 法 (`tcod`) と収束の速さ。(2) 1 周の行列からの周期解と分散。(3) x–y 結合の分解 (`qmdiag`): Edwards–Teng 型の変換 $R$ を、ヘルプの定義と K64 の式で対応づける。(4) Courant–Snyder の行列と Twiss の伝搬。(5) 色収差。(6) 6 次元の拡張 Twiss (ヘルプの $V = PBR_6H$)。放射を含む平衡 (CALC6D の $D$、$\Sigma$) は [14_emittance.md](14_emittance.md)。

**前提**: 各要素の写像と、`trans(6,12)` の 1:6 列 (行列 $T$) ([02_drift.md](02_drift.md)、[04_quad.md](04_quad.md)、[05_bend.md](05_bend.md))。

**読む順**: §1 (結論) → §3 (閉軌道) → §5 (結合) → §6 (Twiss)。

## 1. 結論

- **閉軌道**は、1 周の写像 $f$ の不動点 $f(c) = c$。K64 (`tcod.f`) は Newton 法で解く: 更新 $(T - 1)\,\Delta = c - f(c)$。$T$ は 1 周の行列。収束判定は成分ごとの重み `codw` つきの残差 $r = \sum[(c - f(c))_i/\mathtt{codw}_i]^2 < \mathtt{conv}$。非線形のリング (6 極とずれた四極つき) で、残差は $9\times10^{-5} \to 2\times10^{-6} \to 1\times10^{-9} \to 2.5\times10^{-16}$ と 2 次収束した。
- **分散**は、同じ不動点の方程式を $\delta$ について 1 次に展開した解 $\eta = (1 - M_4)^{-1}m_6$ ($M_4$ は 1 周の $4\times4$、$m_6$ は $\delta$ の列)。追跡で $\delta = \pm10^{-5}$ の閉軌道の差分をとった値と $10^{-6}$ 以内で一致 (例のリングで $\eta_x = 2.28993$ m)。
- **x–y 結合の分解** $R M R^{-1} = \mathrm{diag}(A, B)$ の $R = \begin{pmatrix}c\,I & J r^\top J\\ r & c\,I\end{pmatrix}$ ($c^2 + \det r = 1$) は、K64 の `qmdiag` の閉じた式で求まる。結合した安定な行列で、$R$ はシンプレクティック、$R^{-1}$ は $r$ の符号を逆にしたもの (ヘルプの通り)、非対角ブロックは $10^{-16}$。pysad の `periodic_r` と一致。
- **Twiss**は、ブロック $A$ ($= $ 1 周の 2×2) に $B A B^{-1} = $ 回転 (sympy で確認)。伝搬は $\beta_2 = t_{11}^2\beta - 2t_{11}t_{12}\alpha + t_{12}^2\gamma$ (sympy で確認、ドリフトで pysad の Twiss と $9\times10^{-16}$)。
- **色収差**: 四極だけのリングで、追跡の Jacobian から出した $d\nu/d\delta$ は pysad の `lat.chromaticity()` と一致し、薄いレンズの式 $-\frac1{4\pi}\sum K_1L\,\beta$ の 3% 以内 (薄い近似の誤差)。**偏向磁石のあるリングでは、`lat.chromaticity()` と追跡で値が合わない** (§7。原因は未調査)。
- 6 次元の拡張 Twiss の行列 $H$ (ヘルプ) は、$a^2 + \det H_x + \det H_y = 1$ のとき数値でシンプレクティック ($2\times10^{-16}$)。

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

| 反復 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| 残差 | $8.9\times10^{-5}$ | $2.0\times10^{-6}$ | $1.0\times10^{-9}$ | $2.5\times10^{-16}$ | $2.7\times10^{-19}$ |

($(\Delta_{k+1}) \simeq \Delta_k^2\cdot O(\partial^2 f/\partial c^2)$ の 2 次収束。丸め誤差で $10^{-19}$ で止まる。) 閉軌道は $x = 3.764\times10^{-4}$ m。1 回目の更新 (線形解) は $3.847\times10^{-4}$ m。

## 4. 周期解と分散

小さな $\delta$ のもとで、$f(c, \delta) \simeq M_4\,c + m_6\,\delta$ ($c = (x, p_x, y, p_y)$、$M_4$ は 1 周の左上 4×4、$m_6$ は第 6 列の上 4 成分)。不動点 $c = M_4c + m_6\delta$ から

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

1 周の行列を $M = \begin{pmatrix}P & N\\ Q & K\end{pmatrix}$ (各 2×2) とする。$R M R^{-1} = \mathrm{diag}(A, B)$、すなわち $R M = \mathrm{diag}(A, B)\,R$ を、(12.3) の $R$ ($\mu = c$、$\bar r = J r^\top J$) について 4 つのブロックに分けて書くと

$$
\begin{aligned}
c P + \bar r Q &= A c,\qquad & c N + \bar r K &= A \bar r,\\
r P + c Q &= B r,\qquad & r N + c K &= B c
\end{aligned}
\tag{12.4}
$$

左上と右下から $A = P + \bar r Q/c$、$B = K + rN/c$。$\hat r \equiv r/c$ とおいて下左の式から $A, B$ を消去すると $r$ に対する Riccati 型の方程式

$$
Q + \hat r P - K\hat r - \hat r N\hat r = 0
\tag{12.5}
$$

になる (右上の式は、シンプレクティックの条件でこれと同値)。幾何の意味: $R^{-1}$ の最初の 2 列 $\binom{c\,I}{-r}$ は、モード 1 の不変な 2 次元面 (物理座標で $(y, p_y) = -\hat r\,(x, p_x)$ という形) を張る。(12.5) はこの面が $M$ で不変であるという条件。

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

(Fortran の `trab2 = 4*((tm31+tm24)*(tm42+tm13)+(tm14-tm32)*(tm41-tm23)) + trpq**2` は、$Q - \bar N = \begin{pmatrix}M_{31} + M_{24} & M_{32} - M_{14}\\ M_{41} - M_{23} & M_{42} + M_{13}\end{pmatrix}$ の行列式が $(M_{31} + M_{24})(M_{42} + M_{13}) + (M_{14} - M_{32})(M_{41} - M_{23})$ に等しいことから、上の式と同じ。$\mathrm{trab2} = (\mathrm{tr}A - \mathrm{tr}B)^2$ で、$\ge 0$ が安定の条件。) 成分で書いた Fortran:

$$
\mathtt{aa} = \mp\frac{1}{c\sqrt{\mathrm{trab2}}},\quad
\begin{pmatrix}r_1 & r_2\\ r_3 & r_4\end{pmatrix} = \mathtt{aa}\begin{pmatrix}M_{31} + M_{24} & M_{32} - M_{14}\\ M_{41} - M_{23} & M_{42} + M_{13}\end{pmatrix}
\tag{12.7}
$$

(上の符号 $-$ は $\tau \ge 0$、$+$ は $\tau < 0$。行列は $Q - \bar N$ で、$\bar N = \begin{pmatrix}-M_{24} & M_{14}\\ M_{23} & -M_{13}\end{pmatrix}$。)

**確認** ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (d)): 結合した安定なシンプレクティック行列 ($|M^\top J M - J| = 4\times10^{-16}$、$\nu = (0.31, 0.17)$ の対角形をランダムなシンプレクティック行列で結合させたもの) で、(12.6)(12.7) の $r$、$c$ を作ると
- $c^2 + \det r - 1 = 0$ ($c = 0.995759$、$\det r = 0.008465$)。
- $R$ はシンプレクティック ($|R^\top J_4R - J_4| = 1.1\times10^{-16}$、$J_4 = \mathrm{diag}(J, J)$)、ヘルプの $R^{-1}$ ($r$ の符号を逆) との積は $1.1\times10^{-16}$ で単位行列。
- $R M R^{-1}$ の非対角ブロックの最大 $8.9\times10^{-17}$ (対角ブロックの最大は 3.0)。
- (12.5) の残差 $5.3\times10^{-17}$、$\hat r = \kappa(Q - \bar N)$ との差 $1.4\times10^{-17}$。
- pysad の `optics.periodic_r` の $R$ と $0$ (差なし)、$A$ とは $2.2\times10^{-16}$。

(12.6) が (12.5) を満たすことの代数的な証明 ((12.5) に $\hat r = \kappa(Q - \bar N)$ を代入し、$M^{-1} = -J_4M^\top J_4$ から出る $P\bar K - N\bar Q = I$ などの関係式で整理して $\kappa$ の 2 次方程式 $\kappa^{-2}\ldots$ に帰着させる計算) は、この章では行っていない (**未確認**。数値では $10^{-17}$ で成立)。

**特別な場合** (`qmdiag.f`): $\mathrm{trab2} < 10^{-14}$ ($\tau = 0$ で結合が弱くない: 和共鳴ではなく差共鳴の近く) では、別の式 ($\det N$、$\mathrm{tr}(QK)$ などの場合分け)。$\mathrm{trab2} < 0$ は和共鳴 ($\nu_x + \nu_y$ = 整数) で、安定でない (`stab = .false.`、"Sum resonance")。この枝の導出はこの章では行わない (読みだけ)。

## 6. Twiss パラメータ

### 6.1 2×2 の安定な行列

結合を分解した後の 2×2 の 1 周行列 $A$ は、Courant–Snyder の形

$$
A = \begin{pmatrix}\cos\mu + \alpha\sin\mu & \beta\sin\mu\\ -\gamma\sin\mu & \cos\mu - \alpha\sin\mu\end{pmatrix},\qquad \gamma = \frac{1 + \alpha^2}{\beta}
\tag{12.8}
$$

に書ける ($\det A = 1$ のため)。逆に $\cos\mu = \mathrm{tr}A/2$、$\beta = A_{12}/\sin\mu$、$\alpha = (A_{11} - A_{22})/(2\sin\mu)$ (pysad: `a = (A11 - c)/s`)。$\sin\mu$ の符号は $A_{12}$ の符号で決める ($\beta > 0$)。$|\mathrm{tr}A| < 2$ が安定。

ヘルプの $B = \begin{pmatrix}1/\sqrt\beta & 0\\ \alpha/\sqrt\beta & \sqrt\beta\end{pmatrix}$ と回転 $P(\psi) = \begin{pmatrix}\cos\psi & \sin\psi\\ -\sin\psi & \cos\psi\end{pmatrix}$ について、$B A B^{-1} = P(\mu)$ (sympy で確認、[opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (c))。$\det B = 1$、$B$ はシンプレクティック。つまり $B$ は、軌道を正規化座標 $(X, P_X)$ ( 位相空間で円) に直す。

### 6.2 伝搬

場所 1 から 2 への 2×2 の行列 $T$ ($\det T = 1$) について、$A_2 = TA_1T^{-1}$ (同じ位相の 1 周行列)。(12.8) の $(1,2)$ 成分から

$$
\beta_2 = t_{11}^2\beta_1 - 2t_{11}t_{12}\alpha_1 + t_{12}^2\gamma_1,\qquad
\alpha_2 = -t_{11}t_{21}\beta_1 + (t_{11}t_{22} + t_{12}t_{21})\alpha_1 - t_{12}t_{22}\gamma_1
\tag{12.9}
$$

(sympy で $\beta_2$ の式を $\det T = 1$ のもとで確認。差 0。) 位相の進みは $\tan\Delta\psi = t_{12}/(\beta_1t_{11} - \alpha_1t_{12})$。ドリフト ($t_{11} = 1$、$t_{12} = L$、$t_{21} = 0$、$t_{22} = 1$) では $\beta_2 = \beta_1 - 2\alpha_1L + \gamma_1L^2$、$\alpha_2 = \alpha_1 - \gamma_1L$。pysad の Twiss の行 5 つのドリフトと、$9\times10^{-16}$ の差で一致。

結合があるときは、場所 1 と 2 のそれぞれの $R$ で $R_2TR_1^{-1} = \mathrm{diag}(T_X, T_Y)$ とし ((12.3) の下)、$T_X$、$T_Y$ に (12.9) を使う。

## 7. 色収差

$\xi = d\nu/d\delta$。1 周の行列から: $\cos\mu(\delta) = \mathrm{tr}A(\delta)/2$ の $\delta$ 微分と $\xi = -\dfrac{1}{2\pi\sin\mu}\dfrac{d\cos\mu}{d\delta}$。ただし $\delta \ne 0$ では閉軌道も動く (分散)。四極の寄与は、薄いレンズでは $\delta\mu = \tfrac12\beta\,\delta(K_1L)$、$K_1 \to K_1/(1+\delta)$ から $\xi_{\rm quad} = -\dfrac{1}{4\pi}\sum_i K_{1i}L_i\,\beta_i$ (垂直面は $K_1$ の符号が逆)。

**確認** ([opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (e2)): 偏向なしの 4 極だけのリング ($\nu = (1.1995, 1.0280)$) で、$d\nu/d\delta$ は
| 方法 | $d\nu_x/d\delta$ | $d\nu_y/d\delta$ |
|---|---|---|
| 薄いレンズの式 ($\beta$ は要素の両端の平均) | $-1.1664$ | $-1.1273$ |
| 追跡の Jacobian (差分 $\delta = \pm10^{-4}$) | $-1.2144$ | $-1.1799$ |
| `lat.chromaticity()` | $-1.2144$ | $-1.1799$ |

追跡と `lat.chromaticity()` は一致し、薄いレンズの式は 4% 低い (厚いレンズの補正 $\sim L/\beta$ の分)。

**偏向磁石のあるリング** (例のリング): 追跡 + Newton 法 (分散で動く閉軌道のまわり) の Jacobian からは $(-0.7797, -0.8787)$ ((e)。小数部の $\mu$ の差分)。`lat.chromaticity()` は $(-1.4563, -1.1799)$。垂直の値は四極だけのリングの値 ($-1.1799$) と同じで、**偏向磁石の寄与を含まない**ように見える。両者は大きく違い、**どちらが SAD (gs) の CALC に合うかは未確認**。SAD の `CHROMATICITY` コマンドは四極・6 極だけの閉形式 (BEND を数えない) で、この例では TOTAL が $-30.6$ (SAD の位相 [rad]、2π で割って $-4.87$) と別の量を出す (`tfchro` の式。何の量かはこの章では整理しなかった)。

この食い違いは、`Lattice.chromaticity` (pysad 独自の、要素の `transfer_matrix(δ)` を使う計算、docstring に「SAD の直訳ではない」とある) の問題か、追跡側の扱いか、未切り分け。付録 B に疑いとして書く。

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

**確認**: (12.10) の $H$ (TeX 版) は、$a^2 + \det H_x + \det H_y = 1$ を満たす 5 組のランダムな $H_x, H_y$ で数値的にシンプレクティック (最大 $|H^\top J_6H - J_6| = 2.2\times10^{-16}$、[opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) (f))。これはヘルプの式の整合性 (写しの誤りがないこと) の確認で、K64 の `tecalc` が実際にこの $V$ を作るか、$H_x$ の成分と $\eta$ の対応 (EX・EPX は $R\cdot[H_x; H_y]$ の第 2 列) は、この章では確認していない (**未確認**。`tecalc` は [14_emittance.md](14_emittance.md))。

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
| `opt_01_linear_optics.py` | (a) 閉軌道の Newton 法 (2 次収束)、(b) 分散 (12.2) と追跡の差分・pysad の Twiss、(c) $BAB^{-1} = P(\mu)$ と (12.9) (sympy)・ドリフトでの伝搬、(d) `qmdiag` の式 (12.6)(12.7) の確認 (ブロック対角、Riccati 式)、(e) 色収差 (偏向あり)、(e2) 四極だけのリングの色収差、(f) ヘルプの $H$ のシンプレクティック性 | [opt_01_linear_optics.out.txt](outputs/opt_01_linear_optics.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/opt_01_linear_optics.py` (十数秒)。

## 11. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV。8 セルの FODO リング (QF・QD = 0.3 m、K1 = ±0.35・−0.33 [1/m]、偏向 1.2 m × 16 台の扇形、周長 41.6 m、$\nu = (1.9403, 1.0280)$)。放射・RF なし。偏向の端面なし、フリンジなし。結合のテストは数値の行列 ($\nu = (0.31, 0.17)$ の対角形をシンプレクティックな変換で結合)。

**未確認**:
- (12.6) が (12.5) を満たすことの代数的な証明。`qmdiag` の特別な場合 ($\mathrm{trab2} \approx 0$ の枝)、$\det r \ge 0.625$ の別の形を K64 のどこで作るか。
- 偏向磁石のあるリングの色収差: `lat.chromaticity()` と追跡の差の原因、SAD の CHROMX (CALC で $\delta$ をずらした $\nu$ の差分) との比較。
- $H_x$ と分散 $\eta$ の関係 (§8)、`tecalc` による $V$ の作り方。
- 閉軌道: RADCOD・RADTAPER の分岐 ($dz$ を `trf0` に移す部分、`dptaper`) の導出と確認。
- 結合の強いとき ($\det r$ が大きい) の整数共鳴・和共鳴の近くの挙動。作動点の整数・半整数共鳴は、(12.2) の逆行列で確認する (例のリングは整数共鳴から 0.03〜0.06 離れており、分散が大きい)。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| x-y-coupling の $R$、$R^{-1}$、$\mu^2 + \det r = 1$ | **原文で確認** (TeX 版・平文版とも同じ。成分表示も TeX 版に有る)。K64 の `qmdiag` と対応 (§5)。 |
| $\det r \ge 1$ の別の形と DETR | **原文で確認**: TeX 版は「$\det r \ge 1$」、「$\det r \ge 0.625$ で別の形」、「DETR は $a - \det r$、$a = 1.375$」と、同じ文の中で $1$ と $0.625$ の 2 つの閾値が出てくる。平文版も同じ。通常の形の $\mu$ が実数であるには $\det r \le 1$ でよいので、0.625 は DETR の見かけ上の切り替えの閾値とみられる (未確認)。 |
| extended-Twiss の定義 | **原文で確認**: (12.10)、TeX 版と平文版は同じ式 (記号のみ違う)。 |
| ヘルプに無いもの | (1) 閉軌道の Newton 法 (§3)。(2) 周期解と分散の式 (12.2)。(3) `qmdiag` の閉じた式 (12.6)。(4) Twiss の伝搬 (12.9)。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

1. 偏向磁石のあるリングで、pysad の `Lattice.chromaticity()` と追跡 (Jacobian) の $d\nu/d\delta$ が、水平 $-1.456$ 対 $-0.780$、垂直 $-1.180$ 対 $-0.879$ と大きく違う (§7)。四極だけのリングでは一致。pysad 独自の計算 (K64 の直訳ではない) の問題か、K64 との比較で決まる (未決)。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `codi`・`codf` | 反復の入口の軌道、1 周後の軌道 | `tcod.f` | |
| `dcod1` | $c - f(c)$ = `codi - codf` | `tcod.f` | |
| `trs` | $T - 1$ (5・6 列は RF なしで 0) | `tcod.f` | |
| `fact` | Newton の減衰係数 ($\le 1$) | `tcod.f` | |
| `codw` | 成分ごとの収束の重み (`CODCONV` で倍率) | `tcod.f` (`codw0`) | |
| `conv`・`conv0` | 収束の閾値 $r <$ `conv` ($10^{-4}$) | `tcod.f` | |
| `tm11`…`tm44` | 1 周の 4×4 の行列要素 | `qmdiag.f` (引数) | `M` |
| `trpq` | $\mathrm{tr}(P - K)$ = `tm11 + tm22 - tm33 - tm44` | `qmdiag.f` | `pq` |
| `trab2` | $(\mathrm{tr}A - \mathrm{tr}B)^2$ | `qmdiag.f` | `d` (定義の形が違う) |
| `c` | $\mu = \sqrt{(1 + |\tau|/\sqrt{\mathrm{trab2}})/2}$ | `qmdiag.f` | `a` |
| `aa` | $\mp 1/(c\sqrt{\mathrm{trab2}})$ | `qmdiag.f` | `x` |
| `r1`〜`r4` | $r$ の成分 | `qmdiag.f` | `R` |
