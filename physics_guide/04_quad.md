# 04. 四極 (QUAD)

**この章で分かること**: 四極のハミルトニアンから、(1) 厳密に解ける部分と残りへの分割、(2) 追跡の写像 (K64 の Fortran と pysad の Python)、(3) 転送行列 (CALC・Emittance[] 用)、(4) 端部 (フリンジ) までを、途中の式を省かずに追う。ヘルプが書いていない途中を埋める。式番号 $(4.n)$ から、§12 の Fortran の対応表と §13 の pysad の対応表を引ける。

**前提**: [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) (ハミルトニアンと運動方程式の導出。§9.5 が四極) と、[01_variables_and_conventions.md](01_variables_and_conventions.md) の §1〜§2 と §8 (SAD の変数名の注釈の規約) (正準変数 $(x, p_x, y, p_y, z, \delta)$、追跡の内部変数 $q = p_x/(1+\delta)$、$dv$、K1 の約束)。解析力学と加速器の線形光学。SAD の内部は知らなくてよい。

**読む順**: §1 (結論) → §2 (記号) → §3〜§5 (H と、厳密に解ける 2 つの部分) → §6 (分割と刻みの数) → §7 (端部) → §8〜§10 (行列と追跡の関係)。Fortran・Python のコードを読みながら使うときは §12・§13 から入る。放射・ソレノイド内の四極は別の章 (§11、§15)。

## 1. 結論

- **四極の本体のハミルトニアンは 2 つに分かれ、どちらも厳密に解ける**: $H = H_2 + \Delta H$。$H_2$ は $(x, p_x, y, p_y)$ の 2 次式 (線形な運動方程式) で、解は $\cos$・$\sin$ (集束)・$\cosh$・$\sinh$ (発散)。$\Delta H$ は座標を含まず運動量だけの関数で、運動量が一定のまま座標だけが進む「非線形ドリフト」$D$ になる。$H_2$ の流れを $B$ と書く。
- **$B$ と $D$ は交換しない**ので、本体を $n$ 個の刻みに分け、1 刻みを $D(h/2)\,B(h)\,D(h/2)$ ($h = L/n$) の 2 次の Strang 分割 (K64 の動作) で積分する。v2 (既定) は、3 つの Strang をまとめた 4 次の Yoshida 合成 ($B$ の長さ $1.3512\,h$、$-1.7024\,h$、$1.3512\,h$)。どちらも $B$ と $D$ が厳密にシンプレクティックなので、合成も丸め誤差までシンプレクティック。刻みの数を 2 倍にしたときの誤差は、K64 が $1/4$、v2 が $1/16$ ([quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt))。
- **刻みの数**は $n = 1 + \lfloor L\cdot K1/\varepsilon \rfloor$ ($K1$ は積分強度、$\varepsilon = \varepsilon_0 \cdot$EPS)。軌道のずれの大きさには**よらない**。$\varepsilon_0$ は追跡で K64 が 0.2、v2 が 0.1、CALC で K64 が 0.1、v2 が 0.05。ヘルプの式 $1 + \lfloor 10|K1\,L|/\text{EPS}\rfloor$ は CALC の K64 の値で、追跡には当てはまらない (§6.3)。
- **端部**: 入口・出口に、(i) Maxwell 方程式から出る非線形フリンジ (磁場の縦成分と 3 次の項が作る、$K x^3/12$ 程度の座標のずれ) と、(ii) 磁場の立ち上がりの形 F1・F2 で決まる線形フリンジを、本体の外側で適用する (§7)。どちらも K の 1 次までの硬い端の近似で、$p_x$ などは正準変換 (点変換と指数写像) で更新する。
- **転送行列 (CALC・Emittance[])** は、追跡と同じ $D$・$B$・フリンジの各段を閉軌道のまわりで微分して掛けたもの。$z$ の行は、シンプレクティック条件から組む (§8)。設計軌道上では、$x$ の行列は $\begin{pmatrix}\cos\varphi & \sin\varphi/w\\ -w\sin\varphi & \cos\varphi\end{pmatrix}$ ($w = \sqrt{K1/L}$、$\varphi = wL$) になる。
- **数字の例** (電子 3 GeV、$L = 0.3$ m、K1 = 0.6 /m、$\gamma_0 = 5870.9$。[quad_01_h2_exact.out.txt](outputs/quad_01_h2_exact.out.txt)): 長さあたりの強さ $K = 2.0$ /m²、$\varphi = 0.4243$ rad (24.3°)。$x$ の行列は $\begin{pmatrix}0.911342 & 0.291081\\ -0.582161 & 0.911342\end{pmatrix}$、$y$ は $\begin{pmatrix}1.091358 & 0.309081\\ 0.618163 & 1.091358\end{pmatrix}$。$x$ の実効的な焦点距離は 1.718 m (薄肉近似の $1/K1 = 1.667$ m より 3% 長い、厚みの効果)。
- **精度の例**: 軌道が 3 mm・2 mrad ずれた四極 ($L = 0.3$ m、K1 = −0.6) の 4×4 ヤコビアンの、厳密な写像との差は、K64 が CALC で $1.8 \times 10^{-7}$、追跡で $7.1 \times 10^{-7}$、v2 が CALC で $5.0 \times 10^{-10}$、追跡で $7.9 \times 10^{-9}$ ([quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt))。設計軌道上 (COD が 0) では、追跡のヤコビアンと CALC の行列は $10^{-13}$ 台で一致する。
- **使った条件と限界**は §15 に書く。

![四極の本体の積分器の誤差](figures/quad_integrator_order.png)

図 (`quad_integrator_order.png`): 四極 ($L = 1$ m、K1 = 1.5 /m) を、10 mm・10 mrad の粒子が通るときの、厳密な写像との差。(a) 刻みの数 $n$ と誤差 (傾き $-2$ が Strang、$-4$ が Yoshida)。(b) 本体 $B$ の評価の数 (計算量) あたりの誤差。(c) 1 刻みの構成。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 積分強度 $K1 = B'L/B\rho$ [1/m]。正が水平収束 | `K1` | `ak` (`tsolqu`・`tsolque`・`tquade`)、`ak0` (`tquad`)。符号は回転で処理し、`tsolqu` に渡す `ak` は常に正 | 台本の K1。要素クラス `Quadrupole.K1` は長さあたり $K1/L$ [1/m²] |
| 長さあたりの強さ $K = K1/L$ [1/m²] ($K = eB'/p_0$、式 (4.1)) | | `akk = ak/al` | `akk` |
| 長さ $L$ [m] | `L` | `al` | `el.L`、`L` |
| 本体の刻みの数 $n$、1 刻みの長さ $h = L/n$ | | `ndiv`、`aln0`・`aln1` | `quadbody.ndiv()`、`aln1` |
| 精度パラメータ EPS | `EPS` (既定 1) | `eps0` | `eps0` |
| 全運動量 $p = 1+\delta$ | | `pr` | `pr` |
| 集束の波数 $w = \sqrt{K/p}$ [1/m]、位相 $\varphi = wL$ (刻みでは $w\,h$) | | `w1`、`phi1` (`tzsetparam0`) | `w1`、`phi1` (`track6._tzsetparam0`) |
| $\cos\varphi - 1$、$\varphi - \sin\varphi$、$\sin\varphi$ | | `dc1`、`xs1`、`s1` | `dc1`、`xs1`、`s1` |
| $\cosh\varphi - 1$、$\varphi - \sinh\varphi$、$\sinh\varphi$ | | `dch2`、`xsh2`、`sh2` | `dch2`、`xsh2`、`sh2` |
| 規格化した運動量 $q = p_x/p$ | | 追跡の `px(i)`、`py(i)`。CALC は `pxi = cod(2)/pr` | `track6` の `c[1]`、`c[3]` |
| 速度の差 $dv = 1 - v/v_0$、その $\delta$ 微分 | | `dv`、`dvdp` (`tgetdv`) | `p.dv`、`_dvz` |
| ずれ・回転 | `DX`、`DY`、`ROTATE` | `dx`、`dy`、`theta2`・`theta1` (K1 < 0 は $\pi/2$ 足す・引く) | `DX`、`DY`、`ROTATION` |
| フリンジ | `FRINGE`、`DISFRIN`、`F1`、`F2` | `mfring`、`fringe` (= DISFRIN が 0)、`f1in`・`f2in`・`f1out`・`f2out` | `FRINGE`、`DISFRIN`、`F1`、`F2`。係数は `track6._tsetfringep` |

$z$ の符号は先頭が正 (01 章)。以下、$p$ は一定の媒介変数として扱う ($H$ が $z$ を含まないので $\delta$ は保存する)。

## 3. ハミルトニアン

### 3.1 磁場とベクトルポテンシャル

理想的な四極の磁場は $B_y = B'x$、$B_x = B'y$ ($B'$ は勾配)。粒子の電荷を含めて強さを

$$
K \equiv \frac{eB'}{p_0}\quad [1/\mathrm{m}^2]
\tag{4.1}
$$

と規格化し、**$K > 0$ が水平収束**になるように符号を決める (SAD の K1 の定義。K1 $= K L$)。ベクトルポテンシャルは $\vec{A} = (0, 0, A_s)$ に取れる。$B_x = \partial_y A_s$、$B_y = -\partial_x A_s$ だから、$a_s = eA_s/p_0$ として

$$
a_s = -\frac{K}{2}\left(x^2 - y^2\right)
\tag{4.2}
$$

(これが水平収束 $p_x' = -Kx$ を与える符号は、次の式 (4.3) から $p_x' = -\partial H/\partial x$ を取って確かめられる。) 曲線座標での $A_s$ と磁場の関係 ($B_x = \partial_yA_s/(1 + x/\rho)$ など)、偏向磁石の $A_s$ との対比、(4.2) から運動方程式を出すまでの途中式は、[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §9.3〜§9.5。ラグランジアンから $H$ とハミルトンの方程式を出す手順 (§5)、$t \to s$ の入れ替え (§8) も、そこに書いてある。

### 3.2 $H$、$H_2$、$\Delta H$

01 章の式 (1.1) で、$1/\rho = 0$、$a_x = a_y = 0$、$a_s$ に (4.2) を入れる。**本体の磁場は (4.2) で厳密**なので、これが本体の厳密な $H$:

$$
H = -\sqrt{p^2 - p_x^2 - p_y^2} + \frac{K}{2}\left(x^2 - y^2\right) + \frac{E}{p_0 v_0}
\tag{4.3}
$$

根号を $p_x^2 + p_y^2$ の 2 次まで展開し、残りを $\Delta H$ とする:

$$
H_2 = -p + \frac{p_x^2 + p_y^2}{2p} + \frac{K}{2}\left(x^2 - y^2\right) + \frac{E}{p_0 v_0}
\tag{4.4}
$$

$$
\Delta H = H - H_2 = p - \sqrt{p^2 - p_x^2 - p_y^2} - \frac{p_x^2 + p_y^2}{2p}
\tag{4.5}
$$

(ヘルプの 2nd-order-Hamiltonian と QUAD の transformation の式。ヘルプの本体の式は $E/v_0$ のように $p_0$ の規格化を略して書かれている。) $\Delta H$ は運動量が小さいとき

$$
\Delta H \simeq \frac{(p_x^2 + p_y^2)^2}{8p^3}
\tag{4.6}
$$

で、横方向の運動量の 4 乗の大きさ。

**(4.6) の導出 (展開を 1 段ずつ)**。$Q \equiv p_x^2 + p_y^2$ ($\ll p^2$。条件: 傾き $\sqrt Q/p \ll 1$) とおく。

1. $\sqrt{p^2 - Q} = p\sqrt{1 - Q/p^2}$。$u \equiv Q/p^2$ (2 次の微小量) で $\sqrt{1 - u} = 1 - \dfrac u2 - \dfrac{u^2}{8} - \dfrac{u^3}{16} - \cdots$。
2. したがって $\sqrt{p^2 - Q} = p - \dfrac{Q}{2p} - \dfrac{Q^2}{8p^3} - \dfrac{Q^3}{16p^5} - \cdots$。
3. $\Delta H = p - \sqrt{p^2 - Q} - \dfrac{Q}{2p}$ に入れると、$p$ の項と $Q/(2p)$ の項が相殺して $\Delta H = \dfrac{Q^2}{8p^3} + \dfrac{Q^3}{16p^5} + \cdots$。
4. 最初の項 ($Q$ の 2 乗 = 運動量の 4 乗) を残し、次の $Q^3/(16p^5)$ (6 乗。最初の項との比が $Q/(2p^2)$) を捨てたのが (4.6)。

(sympy の `series` で確認: [ham_06_series.out.txt](outputs/ham_06_series.out.txt) (f)。ただし SAD は、$\Delta H$ を (4.5) のまま**厳密に**使う (§5)。(4.6) は大きさの見積もりだけに使う。)

### 3.3 なぜ分けるか

- $H_2$ は $(x, p_x, y, p_y)$ の 2 次式。運動方程式が線形で、**どんな強さ $KL$ でも厳密に解ける** (§4)。行列は強い四極でも (厚肉のまま) 正しい。ヘルプ (remarks-on-dynamics) の言う「設計軌道のまわりで厳密な線形の変換」はこれ。
- $\Delta H$ は座標 $x, y$ を含まない。運動量 $(p_x, p_y, p)$ が一定で、座標は一定の割合で進む。**厳密に解ける** (§5)。
- $H$ 全体は、$x$ を含む $K(x^2 - y^2)/2$ と、$p_x$ の根号が混ざるので、閉じた式では解けない。
- $\Delta H$ は 4 乗 (式 (4.6))。$H_2$ と $\Delta H$ は交換しない ($H_2$ が $x$ を含み $\Delta H$ が $p_x$ を含む) ので、$H_2$ の流れと $\Delta H$ の流れの積は $H$ の流れと**少し違う**。違いを刻みの数で抑える (§6)。
- どちらの流れも $H$ の形の流れ (正準変換) なので、シンプレクティック性は厳密に保たれる。

## 4. $H_2$ の厳密解

### 4.1 運動方程式と行列

$H_2$ から、$s$ を独立変数とする正準方程式:

$$
x' = \frac{p_x}{p},\quad p_x' = -Kx,\qquad
y' = \frac{p_y}{p},\quad p_y' = +Ky,\qquad
z' = \frac{\partial H_2}{\partial p} = -dv - \frac{p_x^2 + p_y^2}{2p^2},\quad p' = 0
\tag{4.7}
$$

($z'$ の第 1 項は $E/(p_0 v_0)$ の $p$ 微分で、01 章 §4。第 2 項は $\partial[(p_x^2 + p_y^2)/(2p)]/\partial p$。)

**偏微分を 1 項ずつ書く**。(4.4) の $H_2$ の各項を、$T_1 = -p$、$T_2 = (p_x^2 + p_y^2)/(2p)$、$T_3 = K(x^2 - y^2)/2$、$T_4 = E/(p_0v_0)$ として ($s$ が独立変数。ハミルトンの方程式は [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §8.3 の (1a.51a)):

| 項 | $\partial/\partial p_x$ | $\partial/\partial x$ | $\partial/\partial p_y$ | $\partial/\partial y$ | $\partial/\partial p$ |
|---|---|---|---|---|---|
| $T_1 = -p$ | 0 | 0 | 0 | 0 | $-1$ |
| $T_2 = (p_x^2 + p_y^2)/(2p)$ | $p_x/p$ | 0 | $p_y/p$ | 0 | $-(p_x^2 + p_y^2)/(2p^2)$ |
| $T_3 = K(x^2 - y^2)/2$ | 0 | $Kx$ | 0 | $-Ky$ | 0 |
| $T_4 = E/(p_0v_0)$ | 0 | 0 | 0 | 0 | $\partial[E/(p_0v_0)]/\partial p = v/v_0$ |

表の合計から $x' = \partial H_2/\partial p_x = p_x/p$、$p_x' = -\partial H_2/\partial x = -Kx$、$y' = p_y/p$、$p_y' = -\partial H_2/\partial y = +Ky$、$z' = \partial H_2/\partial p = -1 + v/v_0 - (p_x^2 + p_y^2)/(2p^2)$ ($= -dv - (p_x^2+p_y^2)/(2p^2)$。$dv = 1 - v/v_0$)。$p' = -\partial H_2/\partial z = 0$ ($H_2$ が $z$ を含まない)。これが (4.7)。偏微分の表は sympy でも確認した ([ham_06_series.out.txt](outputs/ham_06_series.out.txt) (d))。 規格化した運動量 $q_x = p_x/p$、$q_y = p_y/p$ と $w \equiv \sqrt{K/p}$ で書くと、

$$
x'' = -w^2 x,\qquad y'' = +w^2 y
\tag{4.8}
$$

調和振動子の解を、1 刻み (長さ $\ell$、$\varphi = w\ell$) の更新の形で書く (K64 の Fortran と同じ形)。**水平 (集束)**:

$$
\begin{aligned}
u_1 &= x\,(\cos\varphi - 1) + \frac{q_x}{w}\sin\varphi, & x &\to x + u_1,\\
u_2 &= -x\,w\sin\varphi + q_x\,(\cos\varphi - 1), & q_x &\to q_x + u_2
\end{aligned}
\tag{4.9}
$$

**垂直 (発散)**:

$$
\begin{aligned}
v_1 &= y\,(\cosh\varphi - 1) + \frac{q_y}{w}\sinh\varphi, & y &\to y + v_1,\\
v_2 &= y\,w\sinh\varphi + q_y\,(\cosh\varphi - 1), & q_y &\to q_y + v_2
\end{aligned}
\tag{4.10}
$$

**Fortran の変数名の注釈** (定義は §16 の表): 式 (4.9)(4.10) の $\cos\varphi - 1$ は `dc1`、$\sin\varphi$ は `s1`、$\cosh\varphi - 1$ は `dch2`、$\sinh\varphi$ は `sh2`、$w = \sqrt{K/p}$ は `w1`、$\varphi = w\ell$ は `phi1`、$u_1, u_2, v_1, v_2$ は `u1`、`u2`、`v1`、`v2` ($x$ の更新量、$q_x$ の更新量、$y$ の更新量、$q_y$ の更新量。`tsolqu.f` のローカル変数)。`dc1`、`dch2` のように「$-1$」を引いた形で持つのは、$\varphi$ が小さいとき $\cos\varphi \approx 1$ で桁落ちする引き算を避けるため。

「$x \to x + u_1$」と差の形で書き、$\cos\varphi - 1$ や $\varphi - \sin\varphi$ を級数で直接計算するのは、小さい $\varphi$ の桁落ちを避けるため (SAD の `xsincos`・`xsincosh`、`tfloor.f:260`・`:330`)。正準運動量 $p_x = p\,q_x$ の行列は、$p_x \to p\,q_x$ を代入して

$$
\begin{pmatrix}x\\ p_x\end{pmatrix}_{\rm out}
= \begin{pmatrix}\cos\varphi & \dfrac{\sin\varphi}{pw}\\[2mm] -pw\sin\varphi & \cos\varphi\end{pmatrix}
\begin{pmatrix}x\\ p_x\end{pmatrix},\qquad
\begin{pmatrix}y\\ p_y\end{pmatrix}_{\rm out}
= \begin{pmatrix}\cosh\varphi & \dfrac{\sinh\varphi}{pw}\\[2mm] pw\sinh\varphi & \cosh\varphi\end{pmatrix}
\begin{pmatrix}y\\ p_y\end{pmatrix}
\tag{4.11}
$$

行列式は 1 ($\cos^2 + \sin^2 = 1$、$\cosh^2 - \sinh^2 = 1$)。行列指数 $\exp(AL)$ ($A = \begin{pmatrix}0 & 1/p\\ \mp K & 0\end{pmatrix}$) との一致は [quad_01_h2_exact.out.txt](outputs/quad_01_h2_exact.out.txt) (a) で $10^{-16}$、pysad の `Quadrupole.transfer_matrix` との一致は (b) で厳密に 0 (同じ式)。

**ヘルプの solution-H2 との関係**: ヘルプはソレノイド ($B_z$) つきの一般の解を書く。$B_z = 0$ とすると、$V = \sqrt{(B_z/p)^4 + 4(K/p)^2} = 2K/p$、$w_1 = w_2 = w = \sqrt{K/p}$、$w_+ = 2w$、$U = 1/V = p/(2K)$ になり、係数は $a = x_0$、$b = q_{x0}$、$c = p_{y0}/(2K)$、$d = y_0/2$ ($p_{xm} = p_{x0}$、$p_{ym} = p_{y0}$)。ヘルプの $\Delta x = u_1/w_1 = a(\cos\varphi - 1) + (b/w)\sin\varphi = x_0(\cos\varphi - 1) + q_{x0}\sin\varphi/w$、$\Delta y = w_+v_2/w_2 = y_0(\cosh\varphi - 1) + q_{y0}\sinh\varphi/w$ は、式 (4.9)(4.10) の更新量と一致する。ソレノイドつきの 4×4 は `09_solenoid.md` (予定) に回す。

### 4.2 $K < 0$ と $K = 0$

- **$K < 0$ (水平発散)**: 運動方程式は $(x, p_x, y, p_y) \to (y, p_y, -x, -p_x)$ の 90 度の回転 (これは正準変換) で $K \to -K$ に写る。したがって、K1 < 0 の四極は「座標を 90 度回し、$|K|$ の集束四極として (4.9)(4.10) を適用し、戻す」。結果は、$x$ が $\cosh$・$\sinh$、$y$ が $\cos$・$\sin$ になる。Fortran は `tquade` の `theta1 = theta - m_pi_2`、`ak1 = -ak` (CALC)、追跡は呼び出し側が `theta2` に $\pi/2$ を足して $|K1|$ を渡す (`tquad.f:37` の `ak = sign(ak0, al)`)。回転の向きは追跡と CALC で逆だが、四極は $(x, y) \to (-x, -y)$ で不変なので、結果は同じ。
- **$K = 0$**: $\sin\varphi/w \to \ell$、$\cos\varphi \to 1$ の極限で、ドリフト $\begin{pmatrix}1 & \ell/p\\ 0 & 1\end{pmatrix}$ になる ($K = 10^{-12}$ で (4.11) は ドリフトと 12 桁一致、[quad_01_h2_exact.out.txt](outputs/quad_01_h2_exact.out.txt))。実装は `ak == 0` を別の枝 (`tdrift`・`tdrife`) にしている。
- **$\delta$ 依存**: $\delta$ は $w = \sqrt{K/(1+\delta)}$ と (4.11) の $p$ にだけ入る。エネルギーの高い粒子ほど集束が弱い、という自然クロマティシティはこれ。

### 4.3 $z$ の更新 (導出)

(4.7) から、本体 1 つぶん (長さ $\ell$) の $z$ の変化は

$$
\Delta z = -dv\,\ell - \frac{1}{2}\int_0^\ell \left(q_x^2 + q_y^2\right) ds
\tag{4.12}
$$

軌道に沿った積分は、保存量を使うと閉じた式になる。水平では $C_x \equiv q_x^2 + w^2x^2$ が一定、かつ $d(xq_x)/ds = q_x^2 - w^2x^2$ なので $q_x^2 = [\,C_x + d(xq_x)/ds\,]/2$。垂直では $C_y \equiv q_y^2 - w^2y^2$ が一定、$d(yq_y)/ds = q_y^2 + w^2y^2$ なので $q_y^2 = [\,C_y + d(yq_y)/ds\,]/2$。よって

$$
\Delta z = -dv\,\ell - \frac{1}{4}\Bigl[(C_x + C_y)\,\ell + (x q_x)\big|_{\rm in}^{\rm out} + (y q_y)\big|_{\rm in}^{\rm out}\Bigr]
\tag{4.13}
$$

K64 の Fortran (`tsolqu.f:83-88`) の書き方は、(4.13) の $w^2x^2\ell$ を $w\,x^2\varphi$ に直し、$\varphi$ を $(\varphi - \sin\varphi) + \sin\varphi$ に分けて桁落ちを避けたもの:

$$
\Delta z = -\frac{1}{4}\Bigl[\,w\bigl(x^2\,\mathrm{xs1} - y^2\,\mathrm{xsh2}\bigr) + (q_x^2 + q_y^2)\,\ell + u_1 q_x' + x\,q_x\,\mathrm{dc1} + v_1 q_y' + y\,q_y\,\mathrm{dch2}\Bigr] - dv\,\ell
\tag{4.14}
$$

($\mathrm{xs1} = \varphi - \sin\varphi$、$\mathrm{xsh2} = \varphi - \sinh\varphi$、$\mathrm{dc1} = \cos\varphi - 1$、$\mathrm{dch2} = \cosh\varphi - 1$。$q_x'$ は (4.9) の更新後。) 右辺の括弧は、$s_1$ の項が相殺して (4.13) の括弧に等しい (水平の分は $(x_f - x_i)q_{xf} + x_iq_{xi}\,\mathrm{dc1} = x_fq_{xf} - x_iq_{xi} + w x_i^2\sin\varphi$ を使う)。(4.13) の積分、(4.14)、`quadbody` の本体 $B$ が、$10^{-18}$ 以下で一致する ([quad_01_h2_exact.out.txt](outputs/quad_01_h2_exact.out.txt) (c)(d))。ヘルプには、この $z$ の更新に当たる式が (複素正規座標の形で) 1 行あるだけで、K64 の実座標の形との対応は書かれていない。

## 5. $\Delta H$ の解 (運動学的ドリフト)

$\Delta H(p_x, p_y, p)$ は $x, y, z$ を含まない。座標の速度は一定:

$$
x' = \frac{\partial\Delta H}{\partial p_x} = p_x\left(\frac{1}{p_z} - \frac{1}{p}\right),\quad
y' = p_y\left(\frac{1}{p_z} - \frac{1}{p}\right),\quad
z' = \frac{\partial\Delta H}{\partial p} = 1 - \frac{p}{p_z} + \frac{p_x^2 + p_y^2}{2p^2}
\tag{4.15}
$$

($p_z = \sqrt{p^2 - p_x^2 - p_y^2}$。複素ステップ微分で作った偏微分と手で導いたこの式は $10^{-16}$ で一致、[quad_02_dh_drift.out.txt](outputs/quad_02_dh_drift.out.txt) (a)。) 規格化した運動量 $q$ と $\zeta \equiv p_z/p = \sqrt{1 - q_x^2 - q_y^2}$ では、長さ $\ell$ ぶんの流れは

$$
x \to x + q_x\left(\frac{1}{\zeta} - 1\right)\ell,\quad
y \to y + q_y\left(\frac{1}{\zeta} - 1\right)\ell,\quad
z \to z - \frac{(2 + \zeta)(1 - \zeta)^2}{2\zeta}\,\ell
\tag{4.16}
$$

(最後の式は、(4.15) の $z' = 1 - 1/\zeta + q^2/2$ を $q^2 = 1 - \zeta^2$ で整理したもの。$q = 0$ で 0。) K64 の Fortran は (4.16) を `dpz = sqrt1(-ap)` ($= \zeta - 1$、`sqrt1(a) = a/(1+\sqrt{1+a})` で桁落ちなし)、`r = -dpz/(1+dpz)*aln` ($= (1/\zeta - 1)\ell$)、`x += px*r`、`z -= (3+dpz)*ap/2/(2+dpz)*r` の 4 行で書く (`tsolqu.f:65-66, 71, 83`)。$z$ の係数は、$(3 + {\rm dpz})/(2(2 + {\rm dpz}))\cdot q^2 r = (2 + \zeta)(1 - \zeta)^2/(2\zeta)\ell$ と一致する。DOP853 で (4.15) を積分した結果との差は $10^{-17}$ ([quad_02_dh_drift.out.txt](outputs/quad_02_dh_drift.out.txt) (b))。

小さい傾き $a = q$ で

$$
\Delta x \simeq \frac{a^3}{2}\,\ell,\qquad \Delta z \simeq -\frac{3a^4}{8}\,\ell
\tag{4.17}
$$

(1 m で、傾き 1 mrad なら $\Delta x = 5 \times 10^{-10}$ m、10 mrad なら $5 \times 10^{-7}$ m、30 mrad なら $1.35 \times 10^{-5}$ m、$\Delta z$ は 30 mrad で $-3.0 \times 10^{-7}$ m。[quad_02_dh_drift.out.txt](outputs/quad_02_dh_drift.out.txt) (c))。四極の中で $D$ が効くのは、傾きが数 mrad を超える大きな振幅の粒子のとき。

## 6. 本体の分割

### 6.1 Strang 分割 (2 次、K64)

$B(\tau)$ を $H_2$ の長さ $\tau$ の流れ、$D(\tau)$ を $\Delta H$ の流れとする。本体の長さ $L$ を $n$ 個の刻み ($h = L/n$) に分け、1 刻みを

$$
S_2(h) = D\!\left(\tfrac{h}{2}\right) B(h)\, D\!\left(\tfrac{h}{2}\right)
\tag{4.18}
$$

とする。隣り合う刻みの $D$ は 1 つにまとめる ($n$ 刻みで $D(h/2)\,B(h)\,D(h)\,B(h)\cdots B(h)\,D(h/2)$。`tsolqu.f:61-70` の `if(n == 1) r = r*.5` と最後の半分の $D$ がこの構造。$D$ は $n + 1$ 回、$B$ は $n$ 回)。$S_2$ は前後対称なので、$\log S_2(h)$ は $h$ の奇数次だけを含む:

$$
S_2(h) = \exp\!\Bigl(h\,(\mathcal{B} + \mathcal{D}) + h^3\,C_3 + O(h^5)\Bigr),
\qquad C_3 = \alpha\,[\mathcal{D},[\mathcal{D},\mathcal{B}]] + \beta\,[\mathcal{B},[\mathcal{B},\mathcal{D}]]
\tag{4.19}
$$

($\mathcal{B}, \mathcal{D}$ は $H_2$、$\Delta H$ のリー演算子。$\alpha, \beta$ は BCH の公式の標準の係数。)  $h$ の項が厳密な $H$ の流れで、誤りは $h^3$ から。1 刻みの誤り $O(h^3)$ が $n$ 回たまって、全体は $O(h^2) \propto n^{-2}$。刻みの数を 2 倍にすると誤りは 1/4 になる ([quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (b): $n \ge 4$ で比 3.8〜4.00)。

### 6.2 Yoshida の 4 次の合成 (v2)

$S_2$ を 3 つ、刻みの長さを $w_1h$、$w_0h$、$w_1h$ にして並べる:

$$
S_4(h) = S_2(w_1 h)\,S_2(w_0 h)\,S_2(w_1 h)
\tag{4.20}
$$

(4.19) から $\log S_4 = (2w_1 + w_0)\,h\,(\mathcal{B} + \mathcal{D}) + (2w_1^3 + w_0^3)\,h^3 C_3 + O(h^5)$ (3 つの対称な写像の積も対称なので $h^4$ は出ない)。$h(\mathcal B + \mathcal D)$ が正しい係数で、$h^3$ が消える条件は

$$
2w_1 + w_0 = 1,\qquad 2w_1^3 + w_0^3 = 0
\tag{4.21}
$$

$w_0 = 1 - 2w_1$ を入れると $(1 - 2w_1)^3 = -2w_1^3$、実数の立方根は $1 - 2w_1 = -2^{1/3}w_1$ なので

$$
w_1 = \frac{1}{2 - 2^{1/3}} = 1.351207\ldots,\qquad w_0 = 1 - 2w_1 = -1.702414\ldots
\tag{4.22}
$$

$w_0 < 0$ は、1 刻みの途中で**逆向きに**進む部分があること (図の (c) の 4 番の $B$ と 3・5 番の $D$)。それでも $B$ も $D$ も厳密な流れなので、逆向きの流れも厳密 (cosh/sinh の位相が負になるだけ)。隣り合う $D$ を合わせると、1 刻み (最初と最後以外) は $B$ が 3 回 ($1.3512\,h$、$-1.7024\,h$、$1.3512\,h$)、$D$ が 4 回 ($0.6756\,h$、$-0.1756\,h$、$-0.1756\,h$、$0.6756\,h$ で、前後の刻みとの境の $D$ は $2 \times 0.6756h = w_1 h$ に合わさる)。係数の確認は [quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (a) ($2w_1^3 + w_0^3 = -9 \times 10^{-16}$)。刻みの数を 2 倍にしたときの誤りの比は 15.6〜16.0 (4 次)。

**なぜ 4 次か**: 本体 $B$ の評価は 3 倍に増えるが、同じ計算量なら 4 次の方が誤りが小さい (図の (b))。v2 の既定の刻みの数 (K64 の 2 倍の細かさ) では、四極 ($L = 1$ m、K1 = 1.5、10 mm・10 mrad) の 1 台の誤りは、K64 の刻み 8 で $9.7 \times 10^{-9}$、v2 の刻み 16 ($B$ の評価 48 回) で $1.1 \times 10^{-11}$ ([quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (c))。1 刻みの位相が大きい ($K h^2 = 1.5$ の 1 刻み) と、4 次の方がかえって悪い (中振幅で $6.8\times 10^{-7}$ 対 $4.5\times 10^{-7}$)。既定の刻みの数の式 (§6.3) はそのような 1 刻みにならないが、EPS を大きくして強い四極を 1 刻みにすると起こる。

### 6.3 刻みの数

$$
n = 1 + \left\lfloor \frac{L\cdot K1}{\varepsilon_0\,\mathrm{EPS}} \right\rfloor
\tag{4.23}
$$

$K1$ は積分強度、$L K1 = KL^2 = \varphi^2$ は本体全体の位相の 2 乗、EPS は要素の精度パラメータ (既定 1。0 のときも 1 と同じ)。$n$ は設計軌道の強さだけで決まり、**粒子の振幅によらない**。したがって、振幅が大きい粒子ほど (傾き $q$ が大きく、$\Delta H$ が効くので) 誤りが大きい。

| 経路 | K64 の $\varepsilon_0$ | v2 の $\varepsilon_0$ | 場所 |
|---|---|---|---|
| 追跡 (`tsolqu`) | 0.2 | 0.1 | `tsolqu.f:47`、`quadbody.ndiv()` |
| CALC・Emittance[] (`tsolque`) | 0.1 | 0.05 | `tsolque.f:385`、`quadbody.ndiv_calc()` (CALC と放射つきの追跡は、刻みの数の上限が 1000) |
| ヘルプの式 $1 + \lfloor 10\,\lvert K1\,L\rvert/\mathrm{EPS} \rfloor$ | CALC の K64 の値と同じ ($\varepsilon_0 = 0.1$ のとき $1/0.1 = 10$) | | |

追跡の K64 は 0.2 (つまり $5\,L K1/\mathrm{EPS}$) で、ヘルプの式より刻みが少ない。同じ四極でも、追跡と CALC の刻みの数が違うので、軌道がずれたときの積分の誤りが違う (§1 の数字: 追跡 $7.1 \times 10^{-7}$、CALC $1.8 \times 10^{-7}$)。v2 はこれを、4 次の合成と刻みの半減で、どちらも $10^{-8}$ 以下にしてある。刻みの数の例 (EPS = 1。[quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (c)):

| $L$ [m] | K1 [/m] | $L\,K1$ | ヘルプ | 追跡 K64 | 追跡 v2 | CALC K64 | CALC v2 |
|---|---|---|---|---|---|---|---|
| 0.30 | 0.60 | 0.18 | 2 | 1 | 2 | 2 | 4 |
| 0.50 | 1.50 | 0.75 | 8 | 4 | 8 | 8 | 16 |
| 1.00 | 1.50 | 1.50 | 16 | 8 | 16 | 16 | 31 |
| 0.50 | 4.00 | 2.00 | 21 | 11 | 21 | 21 | 41 |
| 2.00 | 3.00 | 6.00 | 61 | 31 | 61 | 61 | 121 |

($n$ の最後の列は $1 + \lfloor L K1/0.05\rfloor$ で、1.5 / 0.05 は丸めで 30.0000…04 になるため 31。) 放射つきの追跡と CALC では、放射の分割数 (`ndivrad`、`itgetqraddiv`) との大きい方になる (§11)。

## 7. 端部 (フリンジ)

### 7.1 順序と選択

四極の入口・出口には、本体の $H$ (4.3) に含まれない効果が 2 つある。本体の磁場は $K\,\Theta(s)$ (階段) を仮定したが、実際の磁場は端で立ち上がる。この立ち上がりの効果を、ヘルプの remarks-on-dynamics が言う「ハミルトニアンの外の、正準変換による近似」で入れる。順序は (K64 の `tquad`、`tquade`、ヘルプの QUAD の transformation):

> [取り付け誤差・ロール] → **非線形フリンジ (入口)** → **線形フリンジ (入口)** → 本体 ($n$ 刻み) → **線形フリンジ (出口)** → **非線形フリンジ (出口)** → [逆の座標変換]

FRINGE と DISFRIN の組合せ (ヘルプの表): 非線形フリンジは DISFRIN = 0 のときだけ、FRINGE = 0 なら入口と出口の両方、1 なら入口、2 なら出口、3 なら両方。線形フリンジ (F1・F2) は FRINGE = 1、2、3 のとき、その入口・出口。実装の判定式は、追跡が `mfring > -4 && mfring /= 2` (入口)、CALC が `mfring >= 0 && mfring /= 2` (入口) で、FRINGE が $-1 \ldots -3$ のときだけ違う (`tquad.f:49`、`tquade.f:43`)。

### 7.2 非線形フリンジ: Maxwell 方程式から

端の磁場は、$K(s)$ が $s$ とともに変わるとき、$\nabla\times\vec{B} = 0$ と $\nabla\cdot\vec{B} = 0$ から決まる。磁気ポテンシャル $V$ ($\vec{B} = \nabla V$、$\nabla^2 V = 0$) を、$x, y$ の 4 次まで展開して取る:

$$
V = K(s)\,xy - \frac{K''(s)}{12}\,xy\,(x^2 + y^2) + \cdots
\tag{4.24}
$$

($\nabla^2 [xy(x^2 + y^2)] = 12xy$ なので、$\partial_s^2(Kxy) = K''xy$ と相殺する。) $\vec{b} = e\vec{B}/p_0$ として

$$
b_x = Ky - \frac{K''}{12}(3x^2y + y^3),\quad
b_y = Kx - \frac{K''}{12}(x^3 + 3xy^2),\quad
b_s = K'xy
\tag{4.25}
$$

縦成分 $b_s = K'xy$ は、端で $K$ が変わるときにだけ現れる (ヘルプの QUAD の非線形フリンジの出所)。力は $d(p_x, p_y)/ds = (\vec p\times\vec b)_{x,y}/p_z$ で、$p_z \simeq p$ とする。階段 $K(s) = K\,\Theta(s)$ では $K' = K\delta(s)$、$K'' = K\delta'(s)$。**K の 1 次**で、自由な運動 $x(s) = x_0 + x_0' s$ のまわりで端 ($s = 0$ のまわり) を積分する。$\int\delta'(s)\,g(s)\,ds = -g'(0)$ を使うと、$K''$ の項からの運動量のキックは、軌道に沿った微分 $g' = (x^2x' + x'y^2 + 2xyy')/4$ から

$$
\Delta p_x = \frac{K}{p}\left[-\frac{1}{4}(x^2 + y^2)p_x + \frac{1}{2}\,xy\,p_y\right]\quad(\text{入口})
\tag{4.26}
$$

($b_s$ の項 $K\,xy\,p_y/p$ が、$xy\,p_y$ の係数を $-\tfrac12 \to +\tfrac12$ にする。) $K''$ の力は $\delta'(s)$ なので、運動量は $\delta(s)$ 型の「スパイク」を持ち、座標は端で階段状に動く: $\Delta x = (K/12p)(x^3 + 3xy^2)$、$\Delta y = -(K/12p)(3x^2y + y^3)$。これらは次の**リー生成関数** $f$ の流れ ($\Delta x = \partial f/\partial p_x$、$\Delta p_x = -\partial f/\partial x$、…) と一致する:

$$
f_{\rm in} = \frac{K}{12p}\left(x^3 p_x + 3xy^2 p_x - 3x^2y\,p_y - y^3 p_y\right),\qquad f_{\rm out} = -f_{\rm in}
\tag{4.27}
$$

(出口は $K' \to -K\delta$ なので符号が逆。) これが、階段の端の極限 (端の幅が孔径などの長さに比べて十分短い極限) での非線形フリンジ。

**検証**: Maxwell 方程式を (4.24)〜(4.25) の 3 次まで満たす場 ($K(s) = K(1 \pm \tanh(s/d))/2$) の中の運動を、ローレンツ力で DOP853 で厳密に積分し、「線形な場だけ ($K'' = K' = 0$)」の積分との差 (= 非線形フリンジの効果) を、K64 の非線形フリンジの効果と比べた ([quad_05_fringe.out.txt](outputs/quad_05_fringe.out.txt) (b))。$K = 3$ /m²、振幅 $r = 1$〜2 mm で、効果は $8 \times 10^{-10}$〜$6 \times 10^{-9}$ m。**一致は、端の幅 $d = 1$ mm で 0.25〜0.5%、$d = 5$ mm で 3%** で、振幅を半分にしても変わらない (3 次の効果としては合っている)。残りは $d$ と $K$ に比例する小さな補正で、K の 2 次の項と推定する (未確認。フリンジの式は K の 1 次までの近似)。出口も同じ精度で一致する。

### 7.3 母関数と実際の写像

(4.27) は K の 1 次の結果である。指数写像 $\exp(:f:)$ を厳密に使うと、$x, y$ の 4 次以上の項が運動量にも及んで、閉じた式にならない。K64 とヘルプは、(4.27) の 1 次を保ち、シンプレクティックで閉じた式になる**母関数による正準変換**にする。入口 (ヘルプの QUAD の transformation):

$$
G(x_1, p_{x2}, y_1, p_{y2}, p) = H_0 + \frac{1}{2}\left(\frac{\partial H_0}{\partial x_1}\frac{\partial H_0}{\partial p_{x2}} + \frac{\partial H_0}{\partial y_1}\frac{\partial H_0}{\partial p_{y2}}\right),\quad
H_0 = p_{x2}\,\Delta x_1 + p_{y2}\,\Delta y_1
\tag{4.28}
$$

($\Delta x_1 = x_1(a/3 + b)$、$\Delta y_1 = -y_1(a + b/3)$、$a = Kx_1^2/(4p)$、$b = Ky_1^2/(4p)$。出口は $a, b$ の符号が逆。) 型 2 の母関数 $F = x_1 p_{x2} + y_1 p_{y2} + G$ から、$x_2 = \partial F/\partial p_{x2}$、$p_{x1} = \partial F/\partial x_1$ ($\mathbf{p}_1 = (D\Phi)^{\!\top}\mathbf p_2$) を計算すると、**座標は運動量によらない $(x,y)$ の点変換**になる:

$$
x_2 = x_1\left(1 + \frac{a}{3} + b + t\right),\quad
y_2 = y_1\left(1 - a - \frac{b}{3} + t\right),\qquad t = \frac{(a - b)^2}{6}
\tag{4.29}
$$

($t$ は、(4.28) の 2 項目の寄与。$x_1(a/3 + b)(a + b) + \ldots$ を整理して $t\,x_1$ になる。) 運動量は点変換のヤコビアン $J = \partial(x_2, y_2)/\partial(x_1, y_1) = \begin{pmatrix}x_x & h\\ -h & y_y\end{pmatrix}$ の逆転置で更新する ($d = a + b$、$\mathrm{ab} = a - b$):

$$
\begin{pmatrix}p_{x}\\ p_{y}\end{pmatrix}_{\rm out} = \frac{1}{f}\begin{pmatrix}y_y & h\\ -h & x_x\end{pmatrix}\begin{pmatrix}p_{x}\\ p_{y}\end{pmatrix}_{\rm in},\quad
\begin{aligned}
x_x &= 1 + d + \tfrac{1}{6}\,\mathrm{ab}\,(5a - b),\quad y_y = 1 - d + \tfrac{1}{6}\,\mathrm{ab}\,(a - 5b),\\
h &= \frac{K}{2p}\,x_1y_1\left(1 - \tfrac{\mathrm{ab}}{3}\right),\quad f = x_xy_y + h^2
\end{aligned}
\tag{4.30}
$$

(点変換 (4.29) の偏微分: $x_x = \partial x_2/\partial x_1$、$h = \partial x_2/\partial y_1 = -\partial y_2/\partial x_1$、$y_y = \partial y_2/\partial y_1$、$f = \det J$。) $z$ は $\partial G/\partial p$ で、$G$ が $p$ の 1 次 ($H_0$) と 2 次 ($t$ の項) の成分からなることから $z_2 = z_1 - \bigl[p_{x2}(x_2 - x_1) + p_{y2}(y_2 - y_1)\bigr]/p$。つまり**非線形フリンジは $(x, y)$ の点変換で、運動量は逆転置、$z$ はその $p$ 依存の補正**。点変換から作る写像なので、シンプレクティック性は構成から厳密に成り立つ。

[quad_05_fringe.out.txt](outputs/quad_05_fringe.out.txt) の (a): 母関数 (4.28) を複素ステップ微分で処理して作った写像と、K64 の式 (`ttfrin`、nord = 4、`tqfrie`) は、$x, p_x, y, p_y, z$ のすべてで $3 \times 10^{-18}$ の差で一致。1 次の項は (4.27) の流れ ($\Delta x$ は $7.300012 \times 10^{-7}$ に対し $7.300000 \times 10^{-7}$。差 $1.2 \times 10^{-12}$ は $t\,x$) 。(d): 6×6 ヤコビアンの $|M^\top J M - J| = 1.4 \times 10^{-12}$ (中心差分の精度)。

### 7.4 線形フリンジ: F1・F2

端で $K(s) = K_0\,g(s)$ が $0 \to 1$ に立ち上がるとき、硬い端 ($K_0\,\Theta(s - s_0)$) からのずれを $\delta K = K_0\,[g(s) - \Theta(s - s_0)]$ とし、その**モーメント**を定義する:

$$
I_n = \int (s - s_0)^n\,\bigl[g(s) - \Theta(s - s_0)\bigr]\,ds
\tag{4.31}
$$

(ヘルプの $I_n = \int_{-\infty}^{\infty}(s - s_0)^n K_1(s)/K_{10}\,ds$ は、$K_1(s)$ そのものでは積分が発散する。硬い端との差の積分、と読む。) $s_0$ は、有効な長さを決める端の位置で、$I_0 = 0$ になるように選ぶ。

**K の 1 次の摂動**: 水平の線形な運動 $x'' = -K(s)x/p$ の端の前後 (距離 $A$ 前、$B$ 後) の行列の硬い端からのずれは、ドリフト $D(\ell) = \begin{pmatrix}1 & \ell\\ 0 & 1\end{pmatrix}$ で挟んだ摂動の積分:

$$
\delta M = -\int D(B - u)\begin{pmatrix}0 & 0\\ \delta K/p & 0\end{pmatrix}D(u + A)\,du
= \frac{K_0}{p}\begin{pmatrix}I_1 & -(B - A)I_1 + I_2\\ 0 & -I_1\end{pmatrix}\quad(I_0 = 0)
$$

(ここで $u = s - s_0$。) これは、端に挿入した $F = \begin{pmatrix}1 + a & b\\ 0 & 1 - a\end{pmatrix}$ を、$D(B)(F - 1)D(A)$ として前後のドリフトで挟んだもの (中身は $\begin{pmatrix}a & b - a(B - A)\\ 0 & -a\end{pmatrix}$) と等しい。したがって

$$
a = \frac{K_0 I_1}{p},\qquad b = \frac{K_0 I_2}{p}
\tag{4.32}
$$

SAD の F1・F2 は、この $I_1$、$I_2$ を長さの次元に直したもの。ヘルプの定義 $\mathrm{F1} = \mathrm{sgn}(\alpha)\sqrt{|\alpha|}$、$\alpha = 24(I_0^2/2 - I_1)$、$\mathrm{F2} = I_2 - I_0^3/3$ で、$I_0 = 0$ なら $a = -K_0\,\mathrm{F1}\,|\mathrm{F1}|/(24p)$、$b = K_0\,\mathrm{F2}/p$ ($K_0 = K1/L$) で、これが (4.32) と同じ。**直線で立ち上がる (長さ $\ell$) 端は F1 = $\ell$、F2 = 0**、$\tanh(s/d)$ の端は F1 = $\pi d$、F2 = 0 ($I_1 = -\pi^2d^2/24$)。

**線形フリンジの写像** (入口。運動量は規格化した $q$ で書く):

$$
\begin{aligned}
x_2 &= e^{a}x_1 + b\,q_{x1}, & q_{x2} &= e^{-a}q_{x1},\\
y_2 &= e^{-a}y_1 - b\,q_{y1}, & q_{y2} &= e^{a}q_{y1},\\
z_2 &= z_1 - \bigl(a x_1 + b(1 + \tfrac{a}{2})\,q_{x2}\bigr)q_{x1} + \bigl(a y_1 + b(1 - \tfrac{a}{2})\,q_{y2}\bigr)q_{y1}
\end{aligned}
\tag{4.33}
$$

(出口は $a \to -a$、$b$ はそのまま、$z_2 = z_1 + (a x_1 - b(1 - a/2)q_{x2})q_{x1} - (a y_1 - b(1 + a/2)q_{y2})q_{y1}$。) 行列式は $e^{a}e^{-a} = 1$ で、1 次の (4.32) を指数に直して、厳密にシンプレクティックにしてある。$z$ は、$x$ 平面の母関数 $F_2 = e^{a}x_1 p_{x2} + \tfrac12 b'e^{a}p_{x2}^2$ ($b' = b/p$ は正準運動量での係数) の $p$ 微分 ($a \propto 1/p$、$b' \propto 1/p^2$ を使う) $z_2 - z_1 = \partial F_2/\partial p$ から出る。K64 の `tquad` の入口・出口の式と `tqlfre.f` は、(4.33) をそのまま ($q$ または正準運動量に直して) 書いたもの。`tquad.f:501` の `tqlfri` は同じ式に $B_z$ の項を足したサブルーチンだが、この版のソースでは呼ばれていない。

[quad_05_fringe.out.txt](outputs/quad_05_fringe.out.txt) (c): $K(s) = K_0 (s/\ell)^2$ ($0 \le s \le \ell$、$\ell = 0.1$ m) の立ち上がりは、$s_0 = 2\ell/3$ で $I_0 = 0$、$I_1 = -\ell^2/36$、$I_2 = \ell^3/405$ (数値積分で一致)、F1 $= 0.08165$、F2 $= 2.469 \times 10^{-6}$。この端を通る線形の運動を数値積分した 2×2 行列と、硬い端 + (4.32)(4.33) の差は、$K_0 = 10^{-3}$ /m² で $3.5 \times 10^{-12}$ (フリンジなしの差は $2.8 \times 10^{-7}$)、$K_0$ を 2 倍にすると 4 倍 (K の 2 次)。**線形フリンジの式は K の 1 次まで正しい**。Maxwell の場の端での線形部分 (tanh、$d = 1$ mm、$r = 1$ mm) の残りは $2.8 \times 10^{-12}$ (F1 なしでは $1.2 \times 10^{-9}$。[quad_05_fringe.out.txt](outputs/quad_05_fringe.out.txt) (b))。

## 8. 転送行列 (CALC・Emittance[])

CALC・Emittance[] の行列は、追跡と**同じ** $D$・$B$・フリンジの各段を、閉軌道 (COD) のまわりで微分して掛け合わせたもの。K64 の `tquade` は、COD $\xi$ (正準運動量の 6 成分) と行列 `trans` (6×12 の前半が $T$) を同時に更新する: 段 $k$ の写像 $\Phi_k$ と、その COD での微分 $M_k = \partial\Phi_k/\partial\xi\,|_{\xi_{k-1}}$ について

$$
\xi_k = \Phi_k(\xi_{k-1}),\qquad T = M_n M_{n-1}\cdots M_1
\tag{4.34}
$$

(`tmultr5(trans, trans1, irad)` が $T \leftarrow M_k T$)。$\Phi_k$ は上の各節の式そのもの (追跡の写像)。したがって、**設計軌道の近く (COD が 0 のとき) の行列は 4.1 節の行列 (4.11)**、COD が動くと、非線形ドリフト $D$ とフリンジの微分が加わる。各段:

**$D$ (`tqente`)**: 正準運動量 $p_x, p_y$ と $p$ で、$x \to x + r\,p_x$、$r = \ell\,(1/p_z - 1/p)$。偏微分は

$$
\frac{\partial x}{\partial p_x} = r + \frac{\ell\,p_x^2}{p_z^3},\quad
\frac{\partial x}{\partial p_y} = \frac{\ell\,p_xp_y}{p_z^3},\quad
\frac{\partial x}{\partial \delta} = p_x\,\frac{\partial r}{\partial p},\quad
\frac{\partial r}{\partial p} = \ell\left(-\frac{p}{p_z^3} + \frac{1}{p^2}\right)
\tag{4.35}
$$

($y$ は $x \leftrightarrow y$)。$z$ の行も解析的に持つ ($\partial z/\partial p_x \sim O(p^3)$。`tqente.f:32-40`)。

**$B$ (`tsolque` の本体)**: (4.9)(4.10) を、規格化した運動量 $(x, q_x, y, q_y)$ について微分し、$\delta$ の列は $w = \sqrt{K/p}$ の $\delta$ 依存だけ ($\partial w/\partial\delta = -w/(2p)$):

$$
\frac{\partial(x_f, q_{xf})}{\partial(x, q_x)} = \begin{pmatrix}\cos\varphi & \sin\varphi/w\\ -w\sin\varphi & \cos\varphi\end{pmatrix},\qquad
\frac{\partial x_f}{\partial\delta} = \left(-\ell\,x\sin\varphi + q_x\,\frac{\varphi\cos\varphi - \sin\varphi}{w^2}\right)\frac{\partial w}{\partial\delta}\ \ \text{(など)}
\tag{4.36}
$$

(Fortran の `trans1(1,6) = (-aln*xi*s1 + pxi*(phi1*dc1 + xs1)/w1**2)*w1p`、`tsolque.f:428`。`w1p` $= \partial w/\partial\delta$。)

**正準運動量に直す**: 追跡の内部変数 $q = p_x/p$ で作った微分 $M_n$ から、正準運動量の行列 $M$ に直す (K64 の `trans1(1:4,2)/pr`、`trans1(1:4,6)-...`、`trans1(2,1:6)*pr`、`tsolque.f:444-449, 575`):

$$
\frac{\partial F}{\partial p_x} = \frac{1}{p}\frac{\partial F}{\partial q_x},\qquad
\frac{\partial F}{\partial\delta}\bigg|_{p_x} = \frac{\partial F}{\partial\delta}\bigg|_{q_x} - \frac{q_x}{p}\frac{\partial F}{\partial q_x} - \frac{q_y}{p}\frac{\partial F}{\partial q_y}\qquad (F = x_f, y_f, \ldots)
\tag{4.37}
$$

出力の運動量 $p_{xf} = p\,q_{xf}(x, q, \delta)$ の行は、連鎖律から、上の $q_{xf}$ の行を $p$ 倍し、$\partial p_{xf}/\partial\delta$ にだけ $+q_{xf}$ を足す ($p_{yf}$ も同様)。

**$z$ の行: シンプレクティック条件で組む**: $z$ に $H$ が依存しないので、$z$ の行は $x, p_x, y, p_y$ の行・$\delta$ の列から決まる。$M^\top J M = J$ の $(j, \delta)$ の成分 ($j = x, p_x, y, p_y$) から

$$
T_{5j} = -T_{1j}T_{26} + T_{2j}T_{16} - T_{3j}T_{46} + T_{4j}T_{36},\qquad j = 1, 2, 3, 4
\tag{4.38}
$$

(Fortran: `tsolque.f:577`、`tqfrie.f:93`。$T_{66} = 1$、$T_{6j} = 0$ を使う。)  $T_{56}$ ($\partial z/\partial\delta$) は、(4.14) の $\delta$ 微分と $dv$ の $\delta$ 微分 `dvdp*aln`。(4.38) が追跡の数値微分と一致することは、[quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt) で、$1.3 \times 10^{-12}$ 以下 (中心差分の精度)。

pysad の v2 は、`emit6d._jac_tsolque_body` が (4.14) を $x, p_x, y, p_y$ で直接微分し、(4.38) と $10^{-19}$ で一致する ([quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt) (e))。

## 9. 追跡の写像 (`tquad` の流れ)

K64 の追跡は、座標の入口で $q = p_x/(1+\delta)$、$dv$ を作ってから (01 章 §5)、要素の写像を順に適用する。四極 (`tquad`) の流れは:

1. $L = 0$ なら薄肉 (`tthin`)。K1 = 0 ならドリフト (`tdrift_free`)。
2. K1 < 0 は、`ak = |K1L|` と回転角 `theta2`($\pi/2$ を足す) に直す (§4.2)。
3. `tsolrot`: 磁石のずれ (DX, DY) を引き、ROTATE の回転をかける。
4. 非線形フリンジ (入口): `ttfrin(…, 4, ak, al, bz)`。(4.28)〜(4.30)。条件: DISFRIN = 0 かつ `mfring > -4` かつ `mfring /= 2`。
5. 線形フリンジ (入口): `mfring = 1, 3` のとき (4.33)。
6. 本体: `tsolqu` (放射なし) か `tsolqur` (放射あり)。$n$ 刻みの $D$ と $B$ (§6)。
7. 線形フリンジ (出口): `mfring = 2, 3`。
8. 非線形フリンジ (出口): `ttfrin(…, -ak, …)`。
9. `tsolrot` で逆の座標変換。

追跡では、粒子ごとに $dv$、$w$、$\varphi$ などを 1 回計算して (`tzsetparam0`)、$n$ 刻みの中で使い回す。numba 版の追跡 (`track6n.k_quad`) と python 版 (`track6._tquad_ens`) は同じ式 (別実装。丸めが違う)。

## 10. 行列と追跡の関係

| 項目 | 追跡 (`tquad`) | CALC・Emittance[] (`tquade`) |
|---|---|---|
| 変数 | 規格化した運動量 $q = p_x/(1+\delta)$ | 正準運動量 (COD と行列) |
| 本体の分割 | $D\,B\,D$ (Strang)、v2 は Yoshida。`quadbody.group()` を共有 | 同じ (構造が共通) |
| 刻みの数 $\varepsilon_0$ | 0.2 (v2 0.1) | 0.1 (v2 0.05) |
| 本体の式 | (4.9)(4.10)(4.14) | 同じ式を微分 (4.36)〜(4.38) |
| 非線形フリンジ | `ttfrin` (nord = 4)。(4.28)〜(4.30) | `tqfrie`。同じ式の微分。判定は `mfring >= 0` |
| 線形フリンジ | (4.33)。係数は `tsetfringep` (`f1 = -\lvert K\,\mathrm{F1}^2\rvert/24`) | (4.33) の微分 (`tqlfre`。$\mathrm{F1}\,\lvert\mathrm{F1}\rvert$) |
| $z$ の更新 | 非線形の写像 | $z$ の行は (4.38) |
| 放射 | `tradk` を各刻みの $B$ のあと | `tradke`。COD が 0 なら寄与 0 |

- **同じ $H$ の同じ分割を、追跡は写像として、CALC は微分として使う**。したがって、設計軌道の近くでは追跡のヤコビアンと CALC の行列は一致する ($1.8 \times 10^{-13}$、[quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt))。軌道がずれると、刻みの数の違い ($\varepsilon_0$) による積分の誤りの違いだけ食い違う (K64 で $5.3 \times 10^{-7}$、v2 で $7.4 \times 10^{-9}$)。
- 次元の違いに注意: 追跡の `ttfrin` は $q$ で書いてあるが、CALC は正準運動量なので $p$ の因子が入る (たとえば `tqlfre.f` の `bp = b/p`)。
- 違う点: (1) 線形フリンジの係数で、追跡は $\mathrm{F1}$ の符号を捨てる ($f_1 = -\lvert K\,\mathrm{F1}^2\rvert/24$。`tffs.f:969`)。CALC は $-K\,\mathrm{F1}\,\lvert\mathrm{F1}\rvert/24$ で符号が残る (`tqlfre.f:12`)。F1 < 0 は、ヘルプの定義 $\mathrm{F1} = \mathrm{sgn}(\alpha)\sqrt{|\alpha|}$ では起こりうるが、F1 = −0.1 の四極 (K1 = 0.6、$L = 0.3$ m) では、追跡と CALC の 4×4 が $2.1 \times 10^{-3}$ 食い違う (F1 = +0.1 では $1.8 \times 10^{-13}$。[quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt) (f))。この違いは K64 のもので、v2 でも同じ。(2) フリンジの判定式 (§7.1)。

## 11. 放射 (概要)

放射を入れるとき、K64 の追跡 (`tsolqur`) は、刻みの数を $\max(n, n_{\rm rad})$ ($n_{\rm rad}$ は `ndivrad` で、1 刻みの光子数が $O(1)$ 以下になる数) にし、各刻みの $B$ のあとで放射のキック (`tradk`) を入れる (Strang の刻みは、$D(h/2)\,B(h)$ のあとで次の刻みの $D$ が続く)。四極の入口・出口の端も、`f1in`・`f1out` を長さの引数として `tradk` を呼ぶ (`tquad.f:74`、`:110`)。CALC の `tquade` は、COD から求めた分割数 (`itgetqraddiv`) で刻み、各刻みのあとに `tradke` (励起・減衰の行列) を入れる。臨界エネルギー・光子数・分散の式と、行列への入り方は `13_radiation.md` (予定) に書く。この章では、放射がない場合の分割と行列だけを扱う。

## 12. Fortran の対応表

行番号は参照した K64 のソースのもの。

| 式・内容 | ファイル:行 | サブルーチン | 変数・メモ |
|---|---|---|---|
| 追跡の全体 (順序) | `tquad.f:1-116` | `tquad` | `ak = sign(ak0, al)` (37)、`tsolrot` (38、112)、`ttfrin` (50、104)、線形入口 (52-)、`tsolqu`・`tsolqur` (79、83)、線形出口 (86-) |
| 本体の刻み $n$ (4.23) | `tsolqu.f:42-52` | `tsolqu` | `eps = 0.2*eps0` (43-45)、`ndiv = 1+int(abs(al*hypot(ak,bz)/eps))` (47)、`akk = ak/al` (52) |
| $w$、$\varphi$、`s1`・`xs1`・`dc1` 等 | `tsolque.f:217-239` | `tzsetparam0` | `akkp = akk/pr`、`w1 = sqrt(akkp)`、`phi1 = min(phimax, max(-phimax, aln*w1))`、`xsincos`・`xsincosh` |
| $D$ (4.16) | `tsolqu.f:64-72` (刻みの先頭)、`:90-95` (最後) | `tsolqu` | `dpz = sqrt1(-ap)` (65)、`r = -dpz/(1+dpz)*aln0` (66。最初だけ半分、67-69)、`x += px*r` (71)、`z -= (3+dpz)*ap/2/(2+dpz)*r` (83) |
| $B$ の $x, q_x, y, q_y$ (4.9)(4.10) | `tsolqu.f:75-82` | `tsolqu` | `u1, u2, v1, v2` |
| $B$ の $z$ (4.14) | `tsolqu.f:83-88` | `tsolqu` | `-0.25*(w1*(xi**2*xs1-yi**2*xsh2) + ap*aln0 + u1*px + xi*pxi*dc1 + ...) - dv*aln0` |
| 最後の半分の $D$ | `tsolqu.f:90-95` | `tsolqu` | |
| 非線形フリンジ (4.28)〜(4.30) | `tquad.f:325-370` | `ttfrin` (`nord == 4`) | `akk = ak/al/4`、`a, b, ab, t, dx1, dy1, xx, yy, h, f`、`px1 = (px*yy+py*h)/f`、`z -= px1*(dx1+t*xi)+py1*(dy1+t*yi)` (366) |
| 線形フリンジ (4.33) | `tquad.f:52-66` (入口)、`:86-100` (出口)、`:501-` (`tqlfri`。呼ばれていない) | `tquad` | `a = f1in/p`、`b = f2in/p`、`ea = exp(a)` |
| 線形フリンジの係数 | `tffs.f:952-980` | `tsetfringep` | `-abs(akk*f1in*f1in)/24`、`abs(akk)*f2in` |
| CALC の全体 | `tquade.f:1-82` | `tquade` | K1 < 0: `theta1 = theta - m_pi_2`、`ak1 = -ak` (30-35)。`tqfrie` (44、74)、`tqlfre` (47、71)、`tsolque` ($n$ 回、66-68) |
| CALC の本体 | `tsolque.f:313-596` | `tsolque` | `eps = 0.1*eps0` (381)、`ndiv` (385)、`tzsetparam` (395)、`tgetdv` (398)、`tqente` (401、593)、$B$ (407-438)、$\delta$ 列の補正 (444-448)、$z$ 行 (577)、`tmultr5` (584) |
| CALC の $D$ (4.35) | `tqente.f:1-40` | `tqente` | `r = -dpz/pz/pr*al` (20)、`trans2(1,2) = r + pxi*fx` (32)、`cod(5)` (24) |
| CALC の非線形フリンジ | `tqfrie.f:1-129` | `tqfrie` | `akk = ak/al/4` (15)、`dz1` (54)、`trans1(5,1)` (93) |
| CALC の線形フリンジ | `tqlfre.f:1-155` | `tqlfre` | `af1 = -ak/al*f1*abs(f1)/24` (12)、`af2 = ak/al*f2` (13) |
| $q = p_x/p$ と $dv$ | `tconvm.f:1-20` | `tconvm` | `px = px/pr` (14)、`dv = -g*(1+pr)/h1/(h1+pr*h0) + dvfs` (18) |
| $dv$ と $\partial dv/\partial\delta$ (CALC) | `tturne.f:869-882` | `tgetdv` | `dv`、`dvdp = h0/h1**3` |
| `xsincos`・`xsincosh`・`sqrt1`・`p2h` | `tfloor.f:260`・`:330`・`:505`・`:463` | | 桁落ち対策 |

## 13. pysad の対応表

| 内容 | pysad の関数 (ファイル:行) | 対応する Fortran | 注意 (添字・符号・単位・規格化) |
|---|---|---|---|
| 刻みの構造 ($D$・$B$ の並びと長さ) | `quadbody.group` (`quadbody.py:76`)、`ndiv`・`ndiv_calc` (`:59`・`:65`)、`QB_BW` (`:40`) | `tsolqu` の `n` のループ、`ndiv` | 追跡・CALC・numba 版が共有。方式は `scheme()`: k64 は 0 (Strang)、v2 は 1 (Yoshida)。`set_scheme(0/1/None)` で固定 |
| 追跡の $D$・$B$ (python と numba が同じ式) | `quadbody.make_track` (`:108`)。`drift`、`body`、`nonrad` | `tsolqu` | 規格化した運動量、スカラー。`sqrt1`・`tzsetparam0` は呼び出し側 (`track6._sqrt1`・`_tzsetparam0`、`track6n.sqrt1`・`tzsetparam0`) |
| 四極全体の追跡 (python) | `track6._tquad_ens` (`track6.py:662`)、`_tquad` (`:732`)、`_tsolqu` (`:461`)、`_tsolqur_ens` (`:563`) | `tquad`、`tsolqu`、`tsolqur` | 外から入る座標は正準運動量。`_Particle` が $p_x/(1+\delta)$ に直す。K1 < 0 は `th2 = ROTATION + π/2` |
| 四極全体の追跡 (numba) | `track6n.k_quad` (`track6n.py:499`)、`tsolqu1` (`:460`)、`ttfrin4_1` (`:334`)、`tqlfr1` (`:361`) | 同上 | `scheme` は `rs[4]` で渡す (`make_rs`) |
| 非線形フリンジ (追跡) | `track6._ttfrin4` (`:421`) | `ttfrin` (`nord == 4`) | `ak = abs(K1L)`、`al = L`。`c[1]`、`c[3]` は規格化した運動量 |
| 線形フリンジ (追跡) | `track6._tqlfr` (`:632`)、係数 `track6._tsetfringep` (`:621`) | 同上 | `f1 = -abs(akk*F1^2)/24`、`f2 = abs(akk)*F2`。F1 の符号は捨てる (§10) |
| CALC の COD の更新 | `lattice._tquade_cod` (`lattice.py:1656`)、`_tsolque_cod` (`:1080`)、`_tsolque_body_cod` (`:1045`)、`_tqente_cod` (`:1021`)、`_tqfrie_cod` (`:1595`)、`_tqlfre_cod` (`:1631`) | `tquade`、`tsolque`、`tqente`、`tqfrie`、`tqlfre` | 正準運動量。K1 < 0 は −90° 回転 (`(x,px,y,py)→(y,py,-x,-px)`)。`el.K1` も `_tquade_cod` の `K1` 引数も長さあたり (積分強度 `ak1 = abs(K1)*L` を内部で作る)。`_tsolque_cod`・`_tqfrie_cod` の `ak` は積分強度 |
| CALC の解析ヤコビアン (numpy) | `emit6d._jac_tqente` (`emit6d.py:505`)、`_jac_tqfrie` (`:539`)、`_jac_tsolque_body` (`:629`)、`_jac_tqlfre` (`:699`)、`_jac_quad_full` (`:748`)、`_tquade` (`:892`) | `tqente`、`tqfrie`、`tsolque`、`tqlfre` | $z$ の行: `_jac_tsolque_body` は v2 で $x, p_x, y, p_y$ での微分を持つ。`_jac_tqente` は $z$ の行の (5,2)・(5,4) を持たない (K64 の `trans2(5,2) = trans2(1,6)` に当たる $O(p^3)$) |
| CALC の解析ヤコビアン (numba) | `elkern._tqente` (`elkern.py:1093`)、`_jac_tqente` (`:1115`)、`_tsolque_body` (`:1147`)、`_jac_tsolque_body` (`:1182`)、`_tqfrie` (`:1236`)、`_jac_tqfrie` (`:1270`)、`_tqlfre` (`:1367`)、`_jac_tqlfre` (`:1389`)、`k_quad` (`:1569`) | 同上 | numba 版と numpy 版は丸めが違う (v2 で $10^{-12}$ 以内) |
| 1 要素の COD 更新と行列 | `lattice._el_cod_step` (`lattice.py:2616`) | `tquade` | 設計軌道上 (`cod = 0`、$\delta = 0$) は `el.transfer_matrix(0)` の行列 (`elements.Quadrupole.transfer_matrix`、`elements.py:133`)。4×4 と $T_{56}$ だけ。$z$ の行 (5,1..4) は持たない |
| 設計軌道の線形行列 | `elements.Quadrupole.transfer_matrix` | | 長さあたり K1、$\delta$ は `K1/(1+delta)` と px の係数 $1/pr$・$pr$ で入る (4.11)。フリンジなし。別系統 |

## 14. 確認スクリプトと出力

すべて、リポジトリのルートから `py -3.12 physics_guide/scripts/quad_xxx.py` で動く (数十秒以内)。出力は標準出力と `physics_guide/outputs/quad_xxx.out.txt` に出る。`physics_guide/scripts/_common.py` が共通部 (リポジトリのルートを `sys.path` に足す、出力の保存、日本語フォント)。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `quad_01_h2_exact.py` | (a) $H_2$ の行列指数と (4.11) の一致 ($10^{-16}$)、$K \to 0$ の極限。(b) `Quadrupole.transfer_matrix` との一致。(c) $z$ の積分と (4.14) の一致。(d) `quadbody` の本体 $B$ との一致 ($10^{-18}$)。(e) 数値例 | [quad_01_h2_exact.out.txt](outputs/quad_01_h2_exact.out.txt) |
| `quad_02_dh_drift.py` | $\Delta H$ の偏微分 (4.15)、DOP853 との一致 ($10^{-17}$)、$\Delta H$ の大きさ (4.17) | [quad_02_dh_drift.out.txt](outputs/quad_02_dh_drift.out.txt) |
| `quad_03_integrator_order.py` | Yoshida の係数 (4.21)(4.22)。刻みの数 $n$ と誤差 (比 4 と 16)。計算量あたりの誤差。刻みの数の表。図 | [quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt)、[図](figures/quad_integrator_order.png) |
| `quad_04_matrix_vs_map.py` | CALC の行列と追跡の写像の微分の一致。変分方程式との比較。$z$ の行 (4.38)。シンプレクティック性。K64 の動作と v2 の両方 | [quad_04_matrix_vs_map.out.txt](outputs/quad_04_matrix_vs_map.out.txt) |
| `quad_05_fringe.py` | (a) 母関数 (4.28) と K64 の式の一致 ($3 \times 10^{-18}$)。(b) Maxwell の場の端の積分との比較。(c) 線形フリンジ F1・F2。(d) シンプレクティック性 | [quad_05_fringe.out.txt](outputs/quad_05_fringe.out.txt) |
| `quad_06_invariants.py` | $B$ だけの $H_2$ の保存、分割した写像の $H$ の変化 ($n^{-2}$ と $n^{-4}$)、ヤコビアンのシンプレクティック性 | [quad_06_invariants.out.txt](outputs/quad_06_invariants.out.txt) |

主な出力 (抜粋)。四極 ($L = 1$ m、K1 = 1.5 /m) を、10 mm・10 mrad の粒子が通ったときの、厳密な写像 (DOP853) との最大の差 ($quad\_03$):

```
    nd    Strang 誤差    比      Yoshida 誤差   比
     1   4.483e-07        -   6.828e-07        -
     2   1.431e-07     3.13   4.370e-08    15.63
     4   3.806e-08     3.76   2.765e-09    15.80
     8   9.667e-09     3.94   1.735e-10    15.94
    16   2.426e-09     3.98   1.085e-11    15.98
    32   6.072e-10     4.00   6.786e-13    15.99
    64   1.518e-10     4.00   4.247e-14    15.98
```

CALC の行列と追跡のヤコビアン ($quad\_04$。K1 = −0.6、$L = 0.3$ m、3 GeV):

```
K64 (Strang):   設計軌道          CALC 4x4 と J_track の差 1.8e-13
                3 mm・2 mrad       CALC 4x4 と J_track の差 5.3e-07、CALC と変分方程式 1.8e-07、J_track と変分方程式 7.1e-07
v2 (Yoshida):   設計軌道          CALC 4x4 と J_track の差 1.7e-13
                3 mm・2 mrad       CALC 4x4 と J_track の差 7.4e-09、CALC と変分方程式 5.0e-10、J_track と変分方程式 7.9e-09
z の行 (4.38):  追跡の数値微分とシンプレクティック条件で組んだ値の差 5.3e-13 (K64)、1.3e-12 (v2)
```

非線形フリンジと Maxwell の場の端 ($quad\_05$ (b)、$K = 3$ /m²、振幅 $r$ = 1、2 mm): 非線形フリンジの効果 $\Delta$ (ODE と K64 の式) の相対差は、端の幅 $d = 0.5, 1, 2, 5$ mm で 0.13〜0.38%、0.25〜0.5%、1.0〜1.1%、2.9% (振幅に依らない)。

## 15. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV ($\gamma_0 = 5870.9$)、設計軌道とその近く (数 mm・数 mrad)。四極は理想的で、本体は $a_s = -K(x^2 - y^2)/2$ (多極成分なし)。ソレノイド ($B_z$)、加速、放射は含まない。行列の数値は線形・COD 0 (設計軌道) の値。

**この章の範囲外 (別の章)**:
- ソレノイドつきの四極 (`bz ≠ 0`、`tsolqu` の別の枝・`tsolque` の別の枝) と、ヘルプの solution-H2・solution-dH の記号の対応: `09_solenoid.md` (予定)。
- 放射の式 (臨界エネルギー・光子数・励起): `13_radiation.md` (予定)。
- MULT の本体 (K0・K2 以上の成分を持つ多極)・薄肉 (`L = 0`): `07_thin_multipole.md`・`08_mult.md` (予定)。
- 偏向磁石の複合機能 (K1 ≠ 0 の BEND)・端部の一般論: `05_bend.md`・`06_edge_fringe.md` (予定)。
- 取り付け誤差の座標変換の導出: `11_misalignment_and_coord.md` (予定)。

**未確認・限界**:
- **非線形フリンジ**: K の 1 次の硬い端の近似。Maxwell の場 (3 次まで) の端との相対差は、幅 $d = 1$ mm で 0.25〜0.5%、$d = 5$ mm で約 3% (§7.2)。残りが K の 2 次の項かは未確認。端の磁場の展開は 3 次までで、実際の端 ($d$ が孔径程度) では高次の項 ($K''''x^5$ など) が効く。
- **線形フリンジ**: K の 1 次まで。式 (4.31) の $I_n$ のとり方 (硬い端との差のモーメント) は、ヘルプの記述 (全域の積分) の読み替え。F1 < 0 では、追跡と CALC が食い違う (§10)。
- **ヘルプの $z$ の更新**: ヘルプの複素正規座標の式と、K64 の実座標の式 (4.14) の同値性は、$B_z = 0$ で数値的には比べていない。ここでの $z$ の導出は (4.7) の $\partial H_2/\partial p$ から直接行った。
- **刻みの数**: 4 次の方が悪くなる 1 刻みの強い四極 (§6.2)。v2 の追跡と CALC の刻みの基準 (0.1 と 0.05) の違いが残る。
- **CALC の $z$ の行**: pysad の `_el_cod_step` (4D の CALC 用) は $z$ の行 (5,1..4) を持たず、6 次元 (放射つき平衡エミッタンス) の `emit6d` が持つ。`emit6d._jac_tqente` の $z$ の行は (5,2)・(5,4) が欠けている (傾きの 3 乗の大きさ。$|J^\top J_0 J - J_0| = 2 \times 10^{-9}$。軌道が (3 mm, 2 mrad, 1 mm, −1 mrad、$\delta = 0.003$) の例)。
- **sympy**: この環境には sympy がないので、導出の検算は数値 (行列指数・数値積分・複素ステップ微分・DOP853) で行い、導出の代数は本文に手で書いた。
- ヘルプの誤植は [01_variables_and_conventions.md](01_variables_and_conventions.md) の付録 A。この章に関わる分は次の付録 A。

## 16. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。この表の定義は、実際のコード (Fortran: `SAD-master/src`、pysad: `pysad/`) で確認した。以下、「小さい $\varphi$ の桁落ち」とは、$\cos\varphi - 1$、$\varphi - \sin\varphi$ のような、近い数の差が、$\varphi \to 0$ で有効数字を失うこと。`xsincos`・`xsincosh` は、小さい $|\varphi|$ ではこれらを直接テイラー級数で計算して防ぐ。

| Fortran の名前 | 式 | 定義の場所 (Fortran) | pysad での名前 |
|---|---|---|---|
| `w1` | $w = \sqrt{K/p} = $ `sqrt(akkp)`、`akkp = akk/pr` ($p$ = `pr` $= 1 + \delta$、`akk` $= K1/L$) | `tsolque.f:217` `tzsetparam0` | `track6n.py:446` `tzsetparam0` の `w1` |
| `phi1` | $\varphi = w L$ (本体では刻みの長さ `aln0`)。`min(phimax, max(-phimax, aln*w1))` で $\pm$`phimax` に切る (`phimax = 100.d0`、`tsolque.f:2`。極端に強い四極で $\cosh$ がオーバーフローしないため。この理由は推定) | `tsolque.f:234` の `tzsetparam0` | `phi1` (同じ関数。100 で切る) |
| `s1` | $\sin\varphi$ | `xsincos(x,s,xs,c,dc)` の `s` (`tfloor.f:260`) | `geo.py:84` `xsincos` の戻り値 1 つ目 |
| `xs1` | $\varphi - \sin\varphi$ | 同 `xs` | 戻り値 2 つ目 |
| `c1` | $\cos\varphi$ | 同 `c` | (`tzsetparam0` は返さない) |
| `dc1` | $\cos\varphi - 1$ (小さい $\varphi$ では $-\varphi^2/2 + \varphi^4/24 - \cdots$ の級数。大きい $\varphi$ では $-\sin^2\varphi/(1 + \cos\varphi)$) | 同 `dc` | 戻り値 4 つ目 |
| `sh2` | $\sinh\varphi$ | `xsincosh(x,sh,xsh,ch,dch)` の `sh` (`tfloor.f:330`) | `xsincosh` の戻り値 1 つ目 |
| `xsh2` | $\varphi - \sinh\varphi$ (負の値。小さい $\varphi$ では $-\varphi^3/6 - \cdots$) | 同 `xsh` | 戻り値 2 つ目 |
| `ch2` | $\cosh\varphi$ | 同 `ch` | (返さない) |
| `dch2` | $\cosh\varphi - 1$ (大きい $\varphi$ では $\sinh^2\varphi/(1 + \cosh\varphi)$) | 同 `dch` | 戻り値 4 つ目 |
| `u1`、`u2`、`v1`、`v2` | 式 (4.9)(4.10) の $u_1 = x\,\mathrm{dc1} + q_x\,\mathrm{s1}/w$ など | `tsolqu.f` の本体ループ | `quadbody.py:133-150` `body` |
| `pxi`、`pyi`、`ap` | 刻みの入口の $q_x$、$q_y$ (規格化した運動量)、`ap` $= q_x^2 + q_y^2$ | `tsolqu.f` | 同 |
| `dpz`、`r` | `dpz = sqrt1(-ap)` ($= \sqrt{1 - q_x^2 - q_y^2} - 1$)、`r = -dpz/(1+dpz)*aln0` (ΔH の運動学的ドリフトの長さ。最初と最後の刻みは半分) | `tsolqu.f`、`sqrt1` は `tfloor.f:505` ($\sqrt{1+x} - 1$ を桁落ちなく) | `_qb_drift` |
| `dv` | $dv = 1 - v/v_0$ (01 章 §8) | `tgetdv`、`tconvm.f` | `dv` |
| `ndiv`、`aln0` | 刻みの数、1 刻みの長さ (§6.3) | `tsolqu.f` | `quadbody.ndiv()`、`aln1` |

未確認: `dpz`・`r` の式と §5 の式 (4.15)〜(4.17) との対応の行ごとの確認 (§12 の表にある範囲を超えるもの)。

## 付録 A. ヘルプの記述・K64 との違い (四極)

| 箇所 | 内容 |
|---|---|
| ヘルプ dynamics (Lagrangean・Hamiltonian) | 根号内の符号、母関数の $t_0$ の項などの誤植 ([01 章付録 A](01_variables_and_conventions.md#付録-a-ヘルプの記述との違い-この章に関わる分)) |
| ヘルプ QUAD > transformation の刻みの数 | $1 + \lfloor 10\lvert K1\,L\rvert/\mathrm{EPS}\rfloor$ は CALC の値 ($\varepsilon_0 = 0.1$)。追跡は 0.2 で、半分の刻み (§6.3) |
| ヘルプ QUAD > F1 の $I_n$ | $\int_{-\infty}^{\infty}$ は発散する。硬い端との差の積分と読む (§7.4)。$\mathrm{F1} = \mathrm{sgn}(\alpha)\sqrt{\alpha}$ の根号は $\sqrt{\lvert\alpha\rvert}$ |
| ヘルプ QUAD > transformation の $H_{2n}$・$\Delta H$ | $E/v_0$ は $E/(p_0v_0)$ ($p_0$ の規格化を略した書き方) |
| ヘルプに無いもの | (1) $H_2$ の解の導出 (ソレノイドつきの式だけ。$B_z = 0$ への縮約は §4.1)。(2) $\Delta H$ の流れ。(3) 分割の構成 ($D\,B\,D$ の順序、Yoshida)。(4) 非線形フリンジの Maxwell 方程式からの導出と、母関数から点変換への書き換え。(5) 線形フリンジの $a, b$ の導出。(6) 行列と追跡の関係。(7) $z$ の更新 (4.14) |
| K64 の実装の癖 | 追跡の線形フリンジが F1 の符号を捨てる (§10)。追跡と CALC で $\varepsilon_0$ が違う (§6.3) |
