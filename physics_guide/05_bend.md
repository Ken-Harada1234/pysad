# 05. 偏向磁石 (BEND)

**この章で分かること**: 偏向磁石 (BEND) の写像が、(1) 入口の端面 (面までの直進)、(2) 線形フリンジ (F1)、(3) 非線形の Maxwell フリンジ、(4) 本体 (円弧)、(5) 出口の側の同じ段、から成ることを、ハミルトニアンと一様磁場の幾何から途中の式を省かずに導き、K64 の Fortran (追跡の `tbendcore`、CALC の `tbende`) と pysad の関数に対応させる。ヘルプは各段の結果の式だけを書いている。その途中 (円弧の厳密解、端面の幾何、フリンジの導出、行列と追跡の関係) を埋める。式番号 $(5.n)$ から、§10 の Fortran の対応表と §11 の pysad の対応表を引ける。

**前提**: [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) (特に §9.4 の偏向磁石の $A_s$ と運動方程式)、[01_variables_and_conventions.md](01_variables_and_conventions.md) の §1〜§2 と §8 (SAD の変数名の注釈の規約)、[04_quad.md](04_quad.md) (章の書き方。非線形フリンジの母関数、行列と追跡の関係)。解析力学と加速器の線形光学。SAD の内部は知らなくてよい。

**読む順**: §1 (結論) → §2 (記号) → §3〜§4 (H と、K1 つきの本体の H2 の解) → §5 (K1 = 0 の円弧と端面の厳密解) → §6 (フリンジ) → §7 (`tbendcore` の読み方) → §8 (K1 ≠ 0) → §9 (行列と追跡)。コードを読みながら使うときは §10・§11・§14 から入る。K64 の誤りの疑いは付録 B に別にまとめた。

## 1. 結論

- **K1 = 0 の BEND の本体は、$H_2$ の近似ではなく、厳密解**を使う。一様な垂直磁場の中の粒子は、水平面内では円を描き、垂直方向には一定の速さで進む (らせん)。K64 の `tbendcore` は、この円と、端面を表す 2 本の直線の幾何を、閉じた式 (面内の角度は `asinz`) で書いたもの。実験室系の幾何だけで解いた基準と、端面角 (E1・E2・AE1・AE2)・ANGLE と K0 が違う場合 ($\rho_b \neq \rho_0$)・$\delta = 0, +0.02, -0.03$ のすべてで $1.1 \times 10^{-15}$ 以内で一致した ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (b))。
- 追跡の 1 要素は、次の順の合成: [取り付け誤差 (`tbshift`・`tbrot`)] → **入口の面までの直進** → **線形フリンジ (F1)** → **非線形の Maxwell フリンジ** → **円弧の本体** → Maxwell フリンジ (出口) → 線形フリンジ (出口) → **出口の面から基準面への直進** → [逆の取り付け誤差]。端面の 3 段 (直進・F1・Maxwell) と円弧が、`tbendcore` の 1 つのサブルーチンにまとまっている (§7)。
- **端面角の効果**は、面までの直進 (5.20) と円弧の幾何から出る。線形光学では、水平に $+\tan\psi/\rho$、垂直に $-\tan\psi/\rho$ の薄い蹴り (面の傾き $\psi = \mathrm{ANGLE}\cdot\mathrm{E} + \mathrm{AE}$)。垂直の収束は、Maxwell フリンジの $\Delta p_y = -(p_x/p_z)\,y/\rho$ に、面の座標系の $p_x \simeq \sin\psi$ を入れたもの。設計軌道の 2×2 行列は、端面の薄い蹴りと扇形の本体の積と $5 \times 10^{-13}$ で一致する ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (e))。
- **K1 ≠ 0 (複合機能)** は、本体を $H_2$ (曲率つきの四極。厳密に解ける) と、残りの非線形項に分け、スライスごとに交互に適用する (§8)。端面は薄い蹴りで近似する。$H_2$ の解の `xsxkx` の符号が $K1 > 0$ の枝で逆になっている (K64 の誤りの既知の例、付録 B)。
- **非線形の Maxwell フリンジ**は、端の磁場 ($\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ を $y$ の 3 次まで満たすもの) から導く生成関数 $f = \dfrac{y^2}{2\rho}\dfrac{p_x}{\sqrt{p^2 - p_x^2}}$ ($y^4$ の補正を除く) の流れで、閉じた式 (5.35) になる。磁場の端の幅 $d$ を持つ Maxwell の場の中の運動を積分すると、$d \to 0$ の外挿で、$\delta = 0$ と $\delta = 0.3$ のどちらも閉じた式と一致する ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。
- **線形フリンジ (F1)** は、ヘルプの定義 $\mathrm{F1} = 6\int (g - g^2)\,ds$ ($g = B_y(s)/B_0$) から導ける。垂直の収束は $F_1/(6\rho^2)$、$y^3$ の項は $\tfrac23\cdot\dfrac{1}{F_1\rho^2}$、水平の $\delta$ 依存のずれは $F_1^2\delta/(24\rho p)$。この 3 つが同じ F1 だけで書けるのは、端が直線の立ち上がり (と $\tanh$ の一部) のとき。一般の端の形では近似 (§6.3)。
- **転送行列 (CALC・Emittance[])** は、追跡と同じ段を COD のまわりで微分して掛けたもの (§9)。設計軌道の近くでは追跡のヤコビアンと $10^{-12}$ 台で一致する。軌道が大きく外れたとき ($y = 10$ mm、$\delta = 0.05$) は $2 \times 10^{-8}$。違いの主因は、CALC の `tbedgemaxwell` の COD の更新が $\Delta x$ と $\Delta z$ で因子 $p = 1 + \delta$ を落としていること (疑い。付録 B-1)。
- **数字の例** (電子 3 GeV、$\gamma_0 = 5870.85$、扇形、$L = 1$ m、ANGLE = 0.1 rad で $\rho = 10$ m。[bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (g)): $R_{16} = \rho(1 - \cos\theta) = 0.0499583$ m、$R_{26} = \sin\theta = 0.0998334$、$R_{56} = -(L - \rho\sin\theta) + L/\gamma_0^2 = -1.6658 \times 10^{-3}$ m (追跡の中心差分と 9 桁一致)。矩形 (E1 = E2 = 0.5) の垂直の行列要素 $R_{43} = -2\tan\psi/\rho + (\tan\psi/\rho)^2 L = -0.009983$。
- **使った条件と限界**は §13 に書く。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 偏向角 [rad]。正は $-x$ へ曲げる | `ANGLE` | `phi0` (`tbend` の引数。`tbendr` では分割ごとの `phi0n`) | `el.ANGLE` |
| 場の偏向角への補正 (K0)。場の角 $\mathrm{ANGLE} + \mathrm{K0}$ | `K0` | `ak0`、`phib = phi0 + ak0` (`tbend.f:565`) | `el.K0` |
| 長さ (円弧に沿う) | `L` | `al` (追跡は F1 の補正後の $L'$、`tffs.f:1283` の `p_L_BEND`)、CALC の `al0` (補正前) と `al` (補正後) | `el.L`、補正は `lattice._bend_cod_update` など |
| 基準軌道の曲率半径 $\rho_0 = L'/\mathrm{ANGLE}$、曲率 $h = 1/\rho_0$ | | `rho0 = al/phi0` (`tbend.f:613`) | `rho0` |
| 場の曲率半径 $\rho_b = L'/(\mathrm{ANGLE} + \mathrm{K0})$、$h_b = 1/\rho_b$ | | `rhob = al/phib` (`tbend.f:612`)、`drhob = rhob - rho0` (`:614`) | `rhob`、`drhob` |
| 運動量 $p = 1 + \delta$ に対する半径 $\rho_e = \rho_b p$ | | `rhoe = rhob*p` (`tbendcore` の局所変数) | `rhoe` |
| 端面の傾き $\psi_1 = \mathrm{ANGLE}\cdot\mathrm{E1} + \mathrm{AE1}$、$\psi_2$ も同様 [rad] | `E1`・`E2` (ANGLE との比)、`AE1`・`AE2` (角) | `psi1`・`psi2` (角)、`cosp1`・`sinp1`・`cosp2`・`sinp2` | `el.E1`・`el.E2` は**角** (= E1·ANGLE + AE1)。`cosp1` など |
| 面の間の角 $\omega = \mathrm{ANGLE} - \psi_1 - \psi_2$ | | `w`、`cosw`・`sinw`・`sqwh` ($= 1 - \cos\omega$)・`sinwp1` ($= \sin(\omega + \psi_1)$) | `cosw`・`sinw`・`sqwh`・`sinwp1` (`track6._tbend0_ens`) |
| 線形フリンジの長さ F1 (入口 `fb1 = F1 + FB1`、出口 `fb2 = F1 + FB2`) | `F1`・`FB1`・`FB2` | `fb1`・`fb2` | `lattice._bend_fringe_params` |
| FRINGE (線形フリンジの入口・出口の選択)、DISFRIN | `FRINGE`・`DISFRIN` | `mfring` (`FRMD`)、`fringe` (= DISFRIN が 0) | `mfring`・`fringe` |
| 全運動量 $p = 1 + \delta$、規格化した運動量 $q = p_x/p$ | | `pr`・`p`、追跡の `px(i)` は $q_x$ | 同 |
| $dv = 1 - v/v_0$ | | `dv` | `p.dv` |

追跡の内部の運動量は $q = p_x/p$ (01 章 §5)。この章で $p_x$ と書いたものは、断らない限り正準運動量 (CALC と同じ)。$z$ の符号は先頭が正。

## 3. ハミルトニアンと $H_2$

### 3.1 磁場とベクトルポテンシャル

曲線座標 ($h = 1/\rho_0$ は基準軌道の曲率) で、一様な垂直磁場 $B_x = 0$、$B_y = B_0$ を考える。01a §9.3 の式 (1a.60a) から $\partial A_s/\partial x = -(1 + hx)B_0$。$x$ で積分し、ゲージを $y$ と $s$ によらない形に選ぶと

$$
A_s = -\left(x + \frac{h x^2}{2}\right)B_0,\qquad
a_s \equiv \frac{eA_s}{P_0} = -\left(x + \frac{h x^2}{2}\right)h_b,\quad h_b \equiv \frac{eB_0}{P_0} = \frac{1}{\rho_b}
\tag{5.1}
$$

(01a の (1a.60b)・(1a.61) で $\rho_b = \rho_0$ としたもの。K0 ≠ 0 のとき場の曲率 $h_b$ は基準軌道の曲率 $h$ と違う。sympy で $B_x = 0$、$B_y = B_0$ を確認: [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (a)。)

### 3.2 厳密な $H$ と $H_2$

01 章の式 (1.1) で $a_x = a_y = 0$、$a_s$ に (5.1) を入れる ($E/(p_0v_0)$ の項は $z$ の方程式の $dv$ だけに効くので、以下では別に扱う):

$$
H = h_b\left(x + \frac{h x^2}{2}\right) - (1 + hx)\sqrt{p^2 - p_x^2 - p_y^2}
\tag{5.2}
$$

$Q \equiv p_x^2 + p_y^2$ とおき、$x$、$p_x$、$p_y$ を 1 次の微小量として、次数ごとに展開する (条件: $Q \ll p^2$、つまり傾きが 1 に比べて小さい。$x/\rho_0$ は 1 次の量の積としてそのまま残す)。

1. 根号: $\sqrt{p^2 - Q} = p\sqrt{1 - Q/p^2} = p - \dfrac{Q}{2p} - \dfrac{Q^2}{8p^3} - \cdots$ ($Q$ の 1 次が 2 次、2 乗が 4 次)。
2. $(1 + hx)$ をかける: $(1 + hx)\sqrt{p^2 - Q} = p + hpx - \dfrac{Q}{2p} - \dfrac{hxQ}{2p} - \dfrac{Q^2}{8p^3} - \cdots$。
3. $H$ に入れる: 次数ごとに
   - 0 次: $-p$ (定数。運動に効かない)。
   - 1 次: $(h_b - hp)\,x$ ($h_b x$ は $a_s$ から、$-hpx$ は $(1 + hx)\sqrt{\cdot}$ から)。
   - 2 次: $\dfrac{h h_b}{2}x^2 + \dfrac{Q}{2p}$。
   - 3 次: $\dfrac{hxQ}{2p}$。
   - 4 次: $\dfrac{Q^2}{8p^3}$。

   (sympy の `series` で一致: [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (a)。)

2 次までを残したものが $H_2$、残りが $H_{\rm rest}$:

$$
H_2 = -p + (h_b - hp)\,x + \frac{h h_b}{2}x^2 + \frac{p_x^2 + p_y^2}{2p},\qquad
H_{\rm rest} = H - H_2 = \frac{hxQ}{2p} + \frac{Q^2}{8p^3} + \cdots
\tag{5.3}
$$

$h_b = h$ ($\rho_b = \rho_0$、K0 = 0) では $(h_b - hp)x = -h\delta x$ と $hh_bx^2/2 = x^2/(2\rho^2)$ で、01a の式 (1a.62) と同じ。複合機能 (K1 ≠ 0) では、$a_s$ に四極の項が加わり、$H_2$ に $\dfrac{K}{2}(x^2 - y^2)$ ($K = K1/L$ [1/m²]) が足される (§8。$a_s$ の高次の項は $H_{\rm rest}$ 側)。

**$H_2$ が落とした項の大きさ**: 運動方程式で $H_{\rm rest}$ の 3 次の項は振幅の 2 次の力になる。厳密な (5.2) と $H_2$ の写像の差は、$L = 3$ m、$\rho = 10$ m で、初期振幅 $a$ = 1・2・4・8 mm に対して $2.6 \times 10^{-7}$、$1.0 \times 10^{-6}$、$4.2 \times 10^{-6}$、$1.7 \times 10^{-5}$ で、振幅を 2 倍にするごとに 4.01〜4.02 倍 ([bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (f))。

### 3.3 偏微分と運動方程式

$H_2$ の項を 1 つずつ微分する (定数項 $-p$ は落とす。$K$ は複合機能の四極成分):

| 項 | $\partial/\partial p_x$ | $\partial/\partial x$ | $\partial/\partial p_y$ | $\partial/\partial y$ |
|---|---|---|---|---|
| $(h_b - hp)\,x$ | 0 | $h_b - hp$ | 0 | 0 |
| $h h_b x^2/2$ | 0 | $h h_b x$ | 0 | 0 |
| $K(x^2 - y^2)/2$ | 0 | $Kx$ | 0 | $-Ky$ |
| $p_x^2/(2p)$ | $p_x/p$ | 0 | 0 | 0 |
| $p_y^2/(2p)$ | 0 | 0 | $p_y/p$ | 0 |

($s$ が独立変数: $x' = \partial H/\partial p_x$、$p_x' = -\partial H/\partial x$ など。) 表の合計から

$$
x' = \frac{p_x}{p},\quad p_x' = -(K + hh_b)\,x + (hp - h_b),\qquad
y' = \frac{p_y}{p},\quad p_y' = +Ky
\tag{5.4}
$$

追跡の内部変数 $q = p_x/p$ ($p$ は一定) で書き、K64 の名前を使うと ($\rho_e = \rho_b p$)

$$
q_x' = \texttt{drhop} + \texttt{akxsq}\; x,\qquad q_y' = \texttt{akysq}\; y,\qquad
\begin{cases}
\texttt{akxsq} = -\dfrac{hh_b + K}{p} = -\dfrac{1}{\texttt{rhosq}} - \dfrac{\texttt{akk}}{p}\\[2mm]
\texttt{akysq} = \dfrac{K}{p} = \dfrac{\texttt{akk}}{p}\\[2mm]
\texttt{drhop} = h - \dfrac{h_b}{p} = \dfrac{\rho_e - \rho_0}{\texttt{rhosq}}
\end{cases}
\tag{5.5}
$$

**Fortran の変数名の注釈**: `rhosq` ($= \rho_0\rho_e = \rho_0\rho_b p$。`tbend.f:58`。$1/\texttt{rhosq} = hh_b/p$)、`akk` ($= K1/\texttt{al} = K$。`tbend.f:60`。`al` はスライスの長さ)、`akxsq` ($-\omega_x^2$。水平の集束の波数の 2 乗の符号を反対にしたもの。振動する枝が `akxsq < 0`)、`akysq` ($K/p$。垂直。$K > 0$ の四極成分は垂直に発散するので `akysq > 0` が双曲型)、`drhop` (定数の力 $\delta/\rho$ 型の項。分散を作る)。定義は `tbend.f:41-116` の `tbendiinit`。sympy で、(5.5) の右辺が (5.4) から出ることを確認した ([bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (b): 差 0)。

## 4. $H_2$ の厳密解

### 4.1 水平と垂直

(5.5) の水平は、平衡点 $x_e$ のまわりの調和振動:

$$
x_e = \frac{\texttt{drhop}}{-\texttt{akxsq}} = \frac{hp - h_b}{K + hh_b}\ \ (K = 0 \text{ で } \rho_e - \rho_0),\qquad
\omega_x^2 = -\texttt{akxsq},\quad \varphi_x = \omega_x\ell
\tag{5.6}
$$

$$
x(\ell) - x_e = (x_0 - x_e)\cos\varphi_x + \frac{q_0}{\omega_x}\sin\varphi_x,\qquad
q(\ell) = -\omega_x (x_0 - x_e)\sin\varphi_x + q_0\cos\varphi_x
\tag{5.7}
$$

($\ell$ はスライスの長さ。) 1 刻みの更新の形 (K64 の `tbendibody`、`tbend.f:136-160`) では、$\texttt{xiksq} \equiv \texttt{drhop} + \texttt{akxsq}\,x = -\omega_x^2(x - x_e)$ として

$$
\Delta x = \texttt{dcxkx}\cdot\texttt{xiksq} + \texttt{sxkx}\cdot q_0,\qquad
\Delta q = \texttt{sxkx}\cdot\texttt{xiksq} + \texttt{dcx}\cdot q_0
\tag{5.8}
$$

`dcx` $= \cos\varphi_x - 1$、`sxkx` $= \sin\varphi_x/\omega_x$、`dcxkx` $= \texttt{dcx}/\texttt{akxsq} = (1 - \cos\varphi_x)/\omega_x^2$ ((5.7) に代入して確かめる: $\texttt{dcxkx}\cdot\texttt{xiksq} = (\cos\varphi_x - 1)(x_0 - x_e)$)。`dcx`・`spx` ($= \sin\varphi_x$ または $\sinh\varphi_x$) は小さい $\varphi_x$ の桁落ちを避けるため、`sxsin`・`sxsinh` ($\varphi - \sin\varphi$ などを級数で) で作る。`akxsq > 0` (K1 が強い負で水平に発散) のときは $\cos$・$\sin$ が $\cosh$・$\sinh$ になる。K1 = 0 でも $\omega_x^2 = hh_b/p > 0$ で、偏向磁石は水平に弱く集束する (01a の (1a.64))。垂直は `akysq` を波数の 2 乗として同様 ($\texttt{akysq} > 0$ で $\cosh$・$\sinh$、`dcy`・`syky`・`aksy`)。

**分散**: $x_0 = q_0 = 0$、$h_b = h$ ($p = 1 + \delta$) の 1 次の $\delta$ で

$$
R_{16} = \rho\left(1 - \cos\frac{L}{\rho}\right),\qquad R_{26} = \sin\frac{L}{\rho}
\tag{5.9}
$$

($x_e = \rho\delta$ と $\omega_x \simeq 1/\rho$ から。sympy: [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (c)。) ANGLE = 0.1 rad、$L = 1$ m で $R_{16} = 0.049958347$ m、$R_{26} = 0.099833417$ が、CALC の行列と 9 桁一致。

### 4.2 $z$ の更新

$H_2$ の $p$ 微分 (各項を 1 つずつ: $-p \to -1$、$(h_b - hp)x \to -hx$、$Q/(2p) \to -Q/(2p^2)$、$E/(p_0v_0) \to v/v_0 = 1 - dv$) から

$$
z' = -dv - h\,x - \tfrac12\left(q_x^2 + q_y^2\right)
\tag{5.10}
$$

$u = x - x_e$ の振動を使うと、軌道に沿う積分は閉じる。$\int x\,d\ell = x_e\ell + u_0\sin\varphi_x/\omega_x + q_0(1 - \cos\varphi_x)/\omega_x^2$。$q_x^2$ は、保存量 $C_x = q_x^2 + \omega_x^2u^2$ と $d(uq_x)/d\ell = q_x^2 - \omega_x^2u^2$ から $\int q_x^2 d\ell = [C_x\ell + (uq_x)|_{\rm in}^{\rm out}]/2$。垂直も同様 ($C_y = q_y^2 - \texttt{akysq}\,y^2$)。K64 の `hi` ($= \tfrac12(q_x^2 + q_y^2 - \texttt{akxsq}\,x^2 - \texttt{akysq}\,y^2) - \texttt{drhop}\,x$。$H_2$ の値から $-p$ を除いたもの。`tbend.f:149`) は $C_x/2 + C_y/2 - \omega_x^2x_e^2/2$ に等しい。整理すると ($T \equiv \texttt{drhop}\cdot\texttt{xsxkx} - \texttt{sxkx}\,x_0 - \texttt{dcxkx}\,q_0 = -\int x\,d\ell$ とおく。$\texttt{xsxkx} \equiv (\ell - \texttt{sxkx})/\texttt{akxsq}$):

$$
\Delta z = -\frac14\Bigl[\Delta(x q_x) + \Delta(y q_y)\Bigr] + \frac14\left(\frac{5}{\rho_0} - \frac{1}{\rho_e}\right)T - \frac{\texttt{hi}\cdot\ell}{2} - dv\,\ell
\tag{5.11}
$$

$\Delta(xq_x) = x_fq_{xf} - x_iq_{xi} = \texttt{dxf}\cdot q_0 + x_0\,\texttt{dpxf} + \texttt{dxf}\cdot\texttt{dpxf}$ (これが K64 の `dxf*pxi+xi*dpxf+dxf*dpxf`)。**係数 $5/\rho_0 - 1/\rho_e$ の出どころ**: $-h\int x = hT$ と、$-\tfrac12\int q_x^2$ から出る $-\tfrac14 x_e[\omega_x^2x_e\ell - \Delta q_x]$ が $+\tfrac14\,\texttt{drhop}\,T = \tfrac14(1/\rho_0 - 1/\rho_e)T$ になる ($\omega_x^2x_e\ell - \Delta q = \omega_x^2\int x = -\omega_x^2T$ と $\texttt{drhop} = \omega_x^2x_e$ から)。足すと $\left[1/\rho_0 + \tfrac14(1/\rho_0 - 1/\rho_e)\right]T = \tfrac14(5/\rho_0 - 1/\rho_e)T$。

(5.11) を $H_2$ の運動方程式の DOP853 の積分と比べると、$x, q_x, y, q_y$ は $10^{-18}$、$z$ は ($\texttt{xsxkx}$ を上の連続な式で取れば) $10^{-18}$ で一致する ($\rho_0 = \rho_b = 10$ m、$L = 1.3$ m、$\delta = 0.01$、2 mm・1 mrad。K1 = 0、集束・発散の 3 通り: [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (d))。**K64 の `xsxkx` は `akxsq < 0` の枝で符号が逆**で、$z$ に $\delta L^3/(3\rho^2)$ 程度 (この例で $7.3 \times 10^{-5}$ m) の誤差が出る (`tbend.f:82` の `xsxkx=-xspx/akx/akxsq`。`akxsq > 0` の枝 (`:72`) と $\texttt{akxsq} \to 0$ の極限 (`:91`) は上の連続な式と一致)。pysad は v2 で連続な式にしてあり、k64 は K64 のまま ([bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (e)。docs/30 の B250)。

## 5. K1 = 0 の円弧の厳密解と端面

### 5.1 一様磁場の円

K1 = 0 の本体の磁場は一様な $b_y = h_b$ ($\vec b = e\vec B/P_0$ で、(5.1) と同じ)。粒子の機械運動量 $\vec P$ ($|\vec P| = p$、$p_0$ で規格化) は、$d\tau = d\ell/p$ ($\ell$ は 3 次元の経路長) で、$d\vec r/d\tau = \vec P$、$d\vec P/d\tau = \vec P\times\vec b$。水平面 $(X, Z)$ (実験室系: $Z$ が入口の基準面の法線、$X$ が $x$ と同じ向き) で、$h_b$ は一様なので、$(P_X, P_Z)$ は一定の角速度 $h_b$ で回る ($y$ の成分 $P_y$ は一定)。$-X$ へ曲がる向き ($\mathrm{ANGLE} > 0$) を正に取ると、回転角 $\Theta = h_b\tau$ で

$$
\begin{aligned}
P_X(\Theta) &= P_{X0}\cos\Theta - P_{Z0}\sin\Theta,\qquad P_Z(\Theta) = P_{Z0}\cos\Theta + P_{X0}\sin\Theta,\\
X(\Theta) &= X_0 + \rho_b\bigl[P_{X0}\sin\Theta + P_{Z0}(\cos\Theta - 1)\bigr],\qquad
Z(\Theta) = Z_0 + \rho_b\bigl[P_{Z0}\sin\Theta + P_{X0}(1 - \cos\Theta)\bigr],\\
Y(\Theta) &= Y_0 + P_y\,\rho_b\Theta
\end{aligned}
\tag{5.12}
$$

($X' = P_X$ を $\tau$ で積分した。) 円の中心は $(X_c, Z_c) = (X_0 - \rho_bP_{Z0},\ Z_0 + \rho_bP_{X0})$ で、$(X - X_c,\ Z - Z_c) = \rho_b\,(P_Z,\ -P_X)$。3 次元の経路長は $\ell = p\tau = p\rho_b\Theta$。$z$ は $z_{\rm out} = z_{\rm in} - \ell + (v/v_0)L'$ ($= -\ell + (1 - dv)L'$)。これは $z' = -(1 + hx)p/P_z + v/v_0$ (01 章の式 (1.8)) を積分したものと同じ。

**曲線座標の $H$ の運動方程式との一致**: (5.2) の正準方程式 $x' = (1 + hx)p_x/P_z$、$p_x' = -h_b(1 + hx) + hP_z$、$y' = (1 + hx)p_y/P_z$、$p_y' = 0$、$z' = 1 - dv - (1 + hx)p/P_z$ ($P_z = \sqrt{p^2 - p_x^2 - p_y^2}$) を DOP853 で積分した結果と、(5.12) の幾何の差は、$\delta = 0, 0.02, -0.03$ で $1.1 \times 10^{-15}$ ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (f))。

### 5.2 入口の面までの直進

実験室系で、入口の基準面は $Z = 0$。入口の面は、基準軌道の入口の点 (原点) を通り、方向 $\hat e_{x1} = (\cos\psi_1, \sin\psi_1)$、法線 $\hat e_{z1} = (-\sin\psi_1, \cos\psi_1)$ の直線 (正の $\psi_1$ は、面が $x > 0$ の側で下流へ進む向き)。基準面の点 $(x_1, 0)$ から運動量 $(P_X, P_Z)$ の向きに直進して面に当たる点を $x_2\hat e_{x1}$ とする。直進のパラメータ $t$ ($d\vec r/dt = \vec P$) で、面の条件 $(x_1 + P_Xt)\sin\psi_1 = P_Zt\cos\psi_1$ から

$$
t = \frac{x_1\sin\psi_1}{P_Z\cos\psi_1 - P_X\sin\psi_1}
\ \Rightarrow\ x_2 = \frac{x_1}{\cos\psi_1 - (P_X/P_Z)\sin\psi_1}
\tag{5.20}
$$

面の座標系 ($\hat e_{x1}$ が $x$、$\hat e_{z1}$ が $z$) での運動量成分と他の変数 (経路は 3 次元の経路長 $pt$、$P_y$ は一定):

$$
\begin{aligned}
p_{x2} &= P_X\cos\psi_1 + P_Z\sin\psi_1,\qquad y_2 = y_1 + \frac{p_y}{P_Z}\,x_2\sin\psi_1,\\
z_2 &= z_1 - \frac{p}{P_Z}\,x_2\sin\psi_1,\qquad p_{y2} = p_y
\end{aligned}
\tag{5.21}
$$

(これが `tbendcore` の入口の `x1 = x/(cosp1 - pxi*sp1)` ($\texttt{sp1} = \sin\psi_1/p_z$)、`dz1 = x1*sp1`、`y1 = y + pyi*dz1`、`px1 = pxi*cosp1 + pzi*sinp1`。ヘルプの "drift to the entrance face" と同じ式。) 出口の面から出口の基準面への直進は、同じ式の逆 (面の座標系から基準面の座標系へ、$\psi \to \psi_2$ で $x_2 = x_1(\cos\psi_2 + (p_{x2}/p_{z2})\sin\psi_2)$ など) で、ヘルプの "drift from the exit face" に当たる。

**シンプレクティック性**: (5.20)(5.21) に $\delta$ (変えない) を加えた 6 変数 $(x, p_x, y, p_y, z, \delta)$ の写像のヤコビアン $M$ について、$M^\top SM - S$ を 40 桁で評価すると最大 $6 \times 10^{-42}$ ($\psi = 0.3$ rad、$(x, p_x, y, p_y, \delta) = (3\ \mathrm{mm}, 2\ \mathrm{mrad}, 2\ \mathrm{mm}, -1\ \mathrm{mrad}, 0.02)$ の 1 点。[bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (a))。面は固定された平面で、その上の位置と運動量の取り替え (自由粒子の直進) なので、正準変換になる。

### 5.3 面から面への円弧

入口の面の座標系で、位置 $(x_1, 0)$、運動量 $(p_{x1}, p_{z1})$ ($p_{z1} = \sqrt{p_h^2 - p_{x1}^2}$、$p_h^2 = p^2 - p_y^2$ は水平面内の運動量の大きさの 2 乗。ここでは $x_1$ は (5.20) の $x_2$、$p_{x1}$ は (5.21) の $p_{x2}$ に当たる) から始める。(5.12) の円の中心は、面の座標系 ($\hat e_{x1}$, $\hat e_{z1}$) で $(x_1 - \rho_bp_{z1},\ \rho_bp_{x1})$。出口の面の法線 $\hat N_2$ と方向 $\hat e_{f2}$ を実験室系で書く: 出口の基準面は基準軌道の出口 $E = (-\rho_0(1 - \cos\phi),\ \rho_0\sin\phi)$ ($\phi = \mathrm{ANGLE}$) を通り法線 $(-\sin\phi, \cos\phi)$、出口の面はそれを $\psi_2$ だけ傾けたもの:

$$
\hat e_{f2} = (\cos(\phi - \psi_2),\ \sin(\phi - \psi_2)),\qquad \hat N_2 = (-\sin(\phi - \psi_2),\ \cos(\phi - \psi_2))
$$

入口の面の座標系の基底との内積は $\hat e_{x1}\cdot\hat e_{f2} = \cos\omega$、$\hat e_{z1}\cdot\hat e_{f2} = \sin\omega$、$\hat e_{x1}\cdot\hat N_2 = -\sin\omega$、$\hat e_{z1}\cdot\hat N_2 = \cos\omega$ ($\omega = \phi - \psi_1 - \psi_2$)。基準軌道の出口の点 $E$ との内積は $E\cdot\hat N_2 = \rho_0[\sin(\omega + \psi_1) + \sin\psi_2]$、$E\cdot\hat e_{f2} = \rho_0[\cos\psi_2 - \cos(\omega + \psi_1)]$。

**出口の面に当たる条件は、出口の面の座標系の運動量 $p_{x2}$ について線形**: 出口の位置 $\vec r$ は $\vec r - \vec c = \rho_b(P_Z, -P_X)$ で、$(P_Z, -P_X)\cdot\hat N_2 = -\vec P\cdot\hat e_{f2} = -p_{x2}$ (出口の面の座標系の $x$ 運動量)。条件 $(\vec r - E)\cdot\hat N_2 = 0$ は $(\vec c - E)\cdot\hat N_2 - \rho_bp_{x2} = 0$ になり、$\vec c\cdot\hat N_2 = -(x_1 - \rho_bp_{z1})\sin\omega + \rho_bp_{x1}\cos\omega$ を入れて

$$
p_{x2} = -\frac{x_1}{\rho_b}\sin\omega + p_{z1}\sin\omega + p_{x1}\cos\omega - \frac{\rho_0}{\rho_b}\bigl[\sin\psi_2 + \sin(\omega + \psi_1)\bigr]
\tag{5.22}
$$

位置は $x_2 = (\vec r - E)\cdot\hat e_{f2} = (\vec c - E)\cdot\hat e_{f2} + \rho_bp_{z2}$ ($p_{z2} = \sqrt{p_h^2 - p_{x2}^2}$):

$$
x_2 = x_1\cos\omega + \rho_b\bigl(p_{z2} - p_{z1}\cos\omega + p_{x1}\sin\omega\bigr) + \rho_0\bigl[\cos(\omega + \psi_1) - \cos\psi_2\bigr]
\tag{5.23}
$$

運動量の回転角 $\Theta$ は、面内の向きの差から $\Theta = \arcsin(p_{x1}/p_h) - \arcsin(p_{x2}/p_h) + \omega$。垂直と $z$ は (5.12) から

$$
y_2 = y_1 + p_y\rho_b\Theta,\qquad
z_2 = z_1 - p\rho_b\Theta + \frac{v}{v_0}L'\ \ (\text{直進の分は (5.21) で入れてある})
\tag{5.24}
$$

(5.22)(5.23) は、ヘルプの "body of bend" の $p_{x2}$、$x_2$ の式と**同じ**。ヘルプの $s$ の式 ($s = \mathrm{ANGLE}\times\rho_b\times(\cdots)$) の $\mathrm{ANGLE}\times$ は、面内の経路長 $s = \rho_bp_h\Theta$ ($y_2 = y_1 + s\,p_y/p_h$ が (5.24) と一致する) の誤記と読む。(5.22)〜(5.24) とヘルプの式を、(5.12) の幾何 (円と面の交点を数値で解く) と比べると $2.4 \times 10^{-15}$ 以内で一致する ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (d))。

**K64 の `tbendcore` の対応** (§7 に行ごとの表): `ph2` $= p_h/p$ ($= \sqrt{1 - q_y^2}$)、`dph2` $= \texttt{ph2} - 1$、`drho1` $= \texttt{drhob} + \rho_b\delta = \rho_e - \rho_0$、`drho` $= \texttt{drho1} + \rho_e\,\texttt{dph2} = \rho_ep_h/p - \rho_0$、`dpx2` $= \texttt{px1} - \texttt{ph2}\sin\psi_1$ (面の座標系の $q_x$ から、設計の向きの粒子の $\texttt{ph2}\sin\psi_1$ を引いた量)、`dpx3a`・`dpx3`・`px3` は (5.22) を $q$ の単位 (÷$p$) で書いたもの。**`px3` は (5.22) の $-p_{x2}/p$**: `tbendcore` は出口の段を $x$ 運動量の符号を反対にした座標で組んでいる (出力側 `px(i) = -cosp2*dpx3 + ...` で戻す。$\psi_2 = 0$ なら $\texttt{px} = -\texttt{px3}$)。Fortran の本体の式を転記して (5.22) の符号つきの値と比べた (1 点で相対 $10^{-16}$)。`phi0a` $= \phi_0 + \texttt{da} = \Theta$: `y3 = y1 + py2*rhoe*phi0a` は (5.24) の $y_2$ ($\texttt{py2}\cdot\rho_bp = p_y\rho_b$)、`z3 = z2 - phi0*drho1 - da*rhoe - dv*al` は $-p\rho_b\Theta + L' - dv\,L'$ ($-\phi_0(\rho_e - \rho_0) - \texttt{da}\,\rho_e = -\rho_e(\phi_0 + \texttt{da}) + L'$。$\phi_0\rho_0 = L'$) に等しい。

### 5.4 端面角の効果 (線形光学)

設計軌道 ($x = p_x = 0$、$p = 1$、$p_y = 0$) のまわりで 1 次まで展開する。入口の面の座標系では、設計軌道は $p_{x1} = \sin\psi_1$、$p_{z1} = \cos\psi_1$ で面に当たる ((5.21))。この点のまわりの 1 次の変分は、(5.20)(5.21) から $\delta x_1 = \delta x/\cos\psi_1$、$\delta p_{x1} = \cos\psi_1\,\delta p_x$、$\delta p_{z1} = -\sin\psi_1\,\delta p_x$。これを (5.22)(5.23) と出口の面から基準面への直進に代入して 1 次を取る。結果は、標準的な端の行列と扇形の本体の積になる:

$$
\begin{pmatrix}x\\ p_x\end{pmatrix}_{\rm out} = E_x(\psi_2)\,S_x\,E_x(\psi_1)\begin{pmatrix}x\\ p_x\end{pmatrix},\quad
E_x(\psi) = \begin{pmatrix}1 & 0\\ \tan\psi/\rho & 1\end{pmatrix},\quad
S_x = \begin{pmatrix}\cos\theta & \rho\sin\theta\\ -\sin\theta/\rho & \cos\theta\end{pmatrix}
\tag{5.25}
$$

($\theta = \mathrm{ANGLE}$、$\rho_b = \rho_0 = \rho$ のとき。) この展開は sympy で行い、(5.20)〜(5.23) と出口の直進をつないだ写像の 2×2 ヤコビアンが (5.25) の積に 30 桁で一致する ($\psi_1 = 0.03$、$\psi_2 = 0.07$、$\phi = 0.1$、$\rho = 10$ m で差 $10^{-167}$。設計軌道は $x = p_x = 0$ に戻る): [bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (g)。水平の端の効果は発散 ($\Delta p_x = +\tan\psi\,x/\rho$、$\psi > 0$ は面を磁石の外側へ開く向き)。

**垂直**は、Maxwell フリンジの $\Delta p_y = -(p_x/p_z)\,y/\rho$ (5.33) に、面の座標系の $p_{x1} = \sin\psi$、$p_{z1} = \cos\psi$ を入れた $-\tan\psi\,y/\rho$ (収束)。垂直の 2×2 は $E_y(\psi_2)\,D(L)\,E_y(\psi_1)$、$E_y(\psi) = \begin{pmatrix}1 & 0\\ -\tan\psi/\rho & 1\end{pmatrix}$。扇形・矩形 (E1 = E2 = 0.5)・非対称 (E1 = 0.3、E2 = 0.7) で、追跡のヤコビアンの水平・垂直とも $5 \times 10^{-13}$ で一致する ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (e))。矩形では水平の $R_{21} = 0$、垂直の $R_{43} = -0.009983$ (ANGLE = 0.1、$\rho = 10$ m)。

## 6. 端部 (フリンジ)

### 6.1 端の磁場

端の近くの磁場 ($x$ 方向に一様な直線の端。$s$ は面の法線方向) を $g(s) = B_y(s)/B_0$ ($g$ は 0 から 1 へ立ち上がる) と、$\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ から作る。磁気スカラーポテンシャル $V$ ($\vec b = \nabla V$、$\nabla^2V = 0$) を

$$
V = \frac{1}{\rho_b}\left[y\,g - \frac{y^3}{6}g'' + \frac{y^5}{120}g'''' - \cdots\right]
\tag{5.30}
$$

に取ると ($x$ によらない: $\nabla^2V = \partial_y^2V + \partial_s^2V$ の各次数が相殺する)、

$$
b_x = 0,\qquad b_y = \frac{1}{\rho_b}\left(g - \frac{y^2}{2}g'' + \frac{y^4}{24}g'''' - \cdots\right),\qquad b_s = \frac{1}{\rho_b}\left(y\,g' - \frac{y^3}{6}g''' + \cdots\right)
\tag{5.31}
$$

$\nabla\cdot\vec b = \partial_yb_y + \partial_sb_s$ は、$y$ の高次 ($V$ の打ち切りの次数) を除いて 0 ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (a))。力は $d\vec P/ds = \vec P\times\vec b/P_z$ から

$$
\frac{dp_x}{ds} = \frac{p_yb_s}{p_z} - b_y,\qquad
\frac{dp_y}{ds} = -\frac{p_xb_s}{p_z}
\tag{5.32}
$$

以下は (5.31) の最低次、$1/\rho_b$ の 1 次で導く (曲率の効果は端面の幾何 (§5) に入れてある)。

### 6.2 硬い端の Maxwell フリンジ

階段 $g = \Theta(s)$ では $g' = \delta(s)$、$g'' = \delta'(s)$。端は幅がなく、粒子の $y$、$p_x$、$p_z$ は端を通る間に変わらないとすると、(5.32) の $p_y$ の式から

$$
\Delta p_y = -\frac{p_x}{p_z}\,\frac{y}{\rho_b}
\tag{5.33}
$$

($b_s = y\,\delta(s)/\rho_b$ の積分。) $p_x$ の式には $b_y$ の $-y^2g''/2$ の項が $\delta'(s)$ の力 (運動量の $\delta(s)$ 型のスパイク) として入り、座標は端で階段状に動く。これをまとめて、運動量 $p_x$ と座標 $y$ を保つ**正準変換**として書くために、生成関数 $f(y, p_x, p)$ を (5.33) から決める: $\Delta p_y = -\partial f/\partial y$ から

$$
f_{\rm in} = \frac{y^2}{2\rho_b}\,\frac{p_x}{p_z},\qquad p_z \equiv \sqrt{p^2 - p_x^2}
\tag{5.34}
$$

($y$ に依らない項 $g(p_x, p)$ は、$y = 0$ で座標が動かないことから 0)。$f$ は $x$、$p_y$、$z$ を含まないので、$\exp(:f:)$ の流れは**厳密に**閉じた形になる: $y$ と $p_x$ は変わらず ($\partial f/\partial p_y = 0$、$\partial f/\partial x = 0$)、

$$
\Delta x = \frac{\partial f}{\partial p_x} = \frac{y^2}{2\rho_b}\,\frac{p^2}{p_z^3},\qquad
\Delta p_y = -\frac{\partial f}{\partial y} = -\frac{y}{\rho_b}\,\frac{p_x}{p_z},\qquad
\Delta z = \frac{\partial f}{\partial p} = -\frac{y^2}{2\rho_b}\,\frac{p_xp}{p_z^3}
\tag{5.35}
$$

($\partial(p_x/p_z)/\partial p_x = p^2/p_z^3$、$\partial(p_x/p_z)/\partial p = -p_xp/p_z^3$。) 出口は $f_{\rm out} = -f_{\rm in}$ (出口では $\Delta x$、$\Delta p_y$、$\Delta z$ の符号が逆)。(5.34) は、$p_y$ を入れない $p_z$ ($\sqrt{p^2 - p_x^2}$) で書いた近似。ヘルプと K64 の式は、さらに $y^4$ の補正を入れる:

$$
f = \frac{1}{2\rho_b}\left(y^2 - \frac{y^4}{12\rho_b^2}\right)\frac{p_x}{\sqrt{p^2 - p_x^2}},\quad
\Delta x = \frac{b\,y^2}{\rho_b}\,\frac{p^2}{p_v^3},\ \ \Delta p_y = -\frac{a\,y}{\rho_b}\,\frac{p_x}{p_v},\ \ \Delta z = -\frac{b\,y^2}{\rho_b}\frac{p_xp}{p_v^3}
\tag{5.36}
$$

$a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$、$p_v = \sqrt{p^2 - p_x^2}$ (ヘルプの $(1 - y_1^2/(12\rho_b^2))$ と $(1 - y_1^2/(6\rho_b^2))$ の因子に当たる。$a = \partial(\cdot)/\partial y$ の関係から $f$ の $y^4$ の係数 $-1/(24\rho_b^3)$ と整合する)。**$y^4/(12\rho_b^2)$ の係数は、1 次の Maxwell の導出 (上) では出てこない ($1/\rho_b^3$ の 3 次の量) ので、この章では導いていない (§13)**。

**K64 のコードとの対応**: 追跡の `tbendcore` (`tbend.f:458-462`、$q$ の単位で `fa`、`f`、`ff`) と CALC の `tbedgemaxwell` (`tbfrie.f:190-229`、正準運動量) は、(5.36) を書いたもの。

| (5.36) | `tbendcore` (q の単位) | `tbedgemaxwell` (正準) |
|---|---|---|
| $\Delta p_y$ | `py2 = py1 - px1*f` ($f = a\,y/(\rho_ep_v/p)$、`fa = y1/rhoe*sqrt(pv1sqi)`、`pv1sqi = 1/(1-px1**2)`) | `cod(4) - af*pxi` ($\texttt{af} = a\,y/\rho_b/p_v$) |
| $\Delta x$ | `x2 = x1 + ff`、`ff = .25*(f+fa)*y1*pv1sqi` ($\simeq b\,y^2p^2/(\rho_bp_v^3)$) | `cod(1) + bf*pr` ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$) |
| $\Delta z$ | `z2 = z1 - px1*ff` | `cod(5) - bf*pxi` |

**ヘルプの式と `tbendcore` の式の形の違い**: ヘルプは (5.36) の形 (入口 $x_2 = x_1 + y_1^2\left(1 - \dfrac{y_1^2}{12\rho_b^2}\right)\dfrac{p_1^2}{2\rho_b(p_1^2 - p_{x1}^2)^{3/2}}$、$p_{y2} = p_{y1} - \cdots$、$z_2 = z_1 - \cdots$) で書く。`tbendcore` は $q = p_x/p$ の単位で、$\Delta p_y$ に使う `f` ($= a\cdot\texttt{fa}$) と、$\Delta x$ の `ff` ($= \tfrac14(f + fa)\,y\cdot\texttt{pv1sqi}$) を、同じ `fa` ($= y/(\rho_e\sqrt{1 - q_x^2})$) から作る。2 つは**代数的に等しい**: $f + fa = (2 - y^2/(6\rho_b^2))\,fa = 2\left(1 - \dfrac{y^2}{12\rho_b^2}\right)fa$ なので $ff = \tfrac12\left(1 - \dfrac{y^2}{12\rho_b^2}\right)fa\,y/(1 - q_x^2) = \dfrac{y^2}{2\rho_bp}\left(1 - \dfrac{y^2}{12\rho_b^2}\right)(1 - q_x^2)^{-3/2}$。$p^2/p_v^3 = 1/(p(1 - q_x^2)^{3/2})$ を使うと、ヘルプの $\Delta x$ に一致する (展開ではなく恒等式)。$\Delta p_y$ も $-q_xf\cdot p = -(p_x/p_v)\,a\,y/\rho_b$ で (5.36) と同じ。`(f + fa)` の対称な書き方は、ヘルプの $(1 - y^2/(12\rho_b^2))$ を $f$ の $(1 - y^2/(6\rho_b^2))$ から作るための式になっている。この形にした意図は、コードからは確認できない。

**検証 1: 端の磁場の中の運動** (直線座標、$1/\rho = 0.02$ /m、$y = 2$ mm、$q_x = p_x/p = 0.1$、端の形 $g = (1 + \tanh(s/d))/2$ で (5.31) の磁場、DOP853。硬い端との差を取り、$d = 2, 1$ mm から $d \to 0$ へ外挿): $\Delta x$ は $\delta = 0$ と $\delta = 0.3$ のどちらも (5.35) の $p^2/p_z^3$ の式と相対 $10^{-4}$ 以下で一致 (比 1.0000)。**CALC の $p$ の因子 (`bf*pr`) だけにした式とは $\delta = 0.3$ で比 1.3000 (= $p$) ずれる**。$\Delta p_y$ は (5.33) と $\delta = 0$・0.3 で相対 $1.5 \times 10^{-4}$、$1.2 \times 10^{-4}$ ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。

**検証 2: 追跡の式のシンプレクティック性**: 追跡の式 (`f`・`ff`) の 6×6 ヤコビアン (中心差分) の $|J^\top SJ - S|$ は $\delta = 0, \pm0.1$ で $10^{-12}$。CALC の `tbedgemaxwell` の COD の更新は $2.0 \times 10^{-5}$ ($\delta = 0$)、$1.7 \times 10^{-4}$ ($\delta = 0.1$)、$2.5 \times 10^{-4}$ ($\delta = -0.1$) で、シンプレクティックでない (付録 B-1。[bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (c))。

### 6.3 線形フリンジ (F1)

端の磁場が有限の幅で立ち上がるとき、硬い端 ($\Theta$) との違いは、端の形 $g(s) = B_y(s)/B_0$ のいくつかの積分で決まる。ヘルプと K64 は、これを F1 という**1 つの長さ** (ヘルプ: "Length of the slope of the field at the edge") で表す。ヘルプの厳密な定義は

$$
\mathrm{F1} = 6\int_{-\infty}^{\infty}\bigl(g(s) - g(s)^2\bigr)\,ds\qquad(\text{1 つのフリンジについて})
\tag{5.40}
$$

を使い、(5.31)〜(5.32) の入射 $p_x = 0$ の粒子で、硬い端との差を導く (端は直線座標、1 次の $1/\rho_b$ で $p_x$ は $p_x(s) = -G(s)/\rho_b$、$G(s) = \int^sg$。$p$ は全運動量で、$p_z \simeq p$):

1. **$y$ の項**: $dp_y/ds = -p_xb_s/p_z = (y/(\rho_b^2p))\,G\,g'$。$\int Gg'\,ds = [Gg] - \int g^2ds$ (積分は端の全体にわたる。$g$ が端の外側で 0 と 1 になるので、$[Gg]$ は $\int g\,ds$ に等しい) $= \int(g - g^2)ds = \mathrm{F1}/6$。$\Delta p_y = y\,\mathrm{F1}/(6\rho_b^2\,p)$。
2. **$y^3$ の項**: $b_s$ の $-y^3g'''/(6\rho_b)$ からの寄与は、$-\dfrac{y^3}{6\rho_b^2p}\int Gg'''\,ds = -\dfrac{y^3}{6\rho_b^2p}\int g'^2ds$ (部分積分 2 回)。さらに、$p_x$ が $b_y$ の $-y^2g''/(2\rho_b)$ の項で $+\dfrac{y^2}{2\rho_b}g'$ だけ変わり、$dp_y/ds = -p_xb_s/p_z$ に $-\dfrac{y^3}{2\rho_b^2p}g'^2$ を足す。合計 $-\left(\tfrac16 + \tfrac12\right)\dfrac{y^3}{\rho_b^2p}\int g'^2ds = -\dfrac23\,\dfrac{y^3}{\rho_b^2p}\int g'^2ds$。
3. **$x$ の項**: $dx/ds = p_x/p_z \simeq p_x/p$ から、硬い端との差 $\Delta X = -\dfrac{1}{\rho_bp}\int(G - G_{\rm hard})\,ds = \dfrac{I_1}{\rho_bp}$ ($I_1 = \int(s - s_0)(g - \Theta(s - s_0))\,ds$、$s_0$ は $I_0 = 0$ になる位置)。

端が**直線の立ち上がり (長さ $\ell$)** のとき、$\int g(1 - g) = \ell/6$、$\int g'^2 = 1/\ell$、$I_1 = -\ell^2/24$。つまり $\mathrm{F1} = \ell$ で、3 つの量が

$$
\Delta p_y = \frac{y}{p}\left(\frac{\mathrm{F1}}{6\rho_b^2} - y^2\,\frac{2}{3\,\mathrm{F1}\,\rho_b^2}\right),\qquad
\Delta X = -\frac{\mathrm{F1}^2}{24\rho_bp}
\tag{5.41}
$$

になる (ヘルプの $\Delta y_{fr} = \mathrm{F1}/(6\rho_b^2)$、$\Delta y_{fra} = 2/(3\,\mathrm{F1}\,\rho_b^2)$)。$\tanh(s/d)$ の端 ($\mathrm{F1} = 3d$) では、$\int g'^2 ds = 1/(3d) = 1/\mathrm{F1}$ で $y^3$ の項は同じ形になる。$x$ の項は $I_1 = -\pi^2d^2/24$ なので、(5.41) の $\mathrm{F1}^2/24 = 9d^2/24$ の $\pi^2/9 = 1.097$ 倍になる。誤差関数の端 ($\sigma$ = 2 mm) では、$\int g'^2 ds\cdot\mathrm{F1}$ が 0.955 (-4.5%) ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (d))。**ヘルプ・K64 の式は、1 つの F1 から 3 つの量を作る (直線の立ち上がりで厳密、一般の端の形では近似)**。

**設計軌道との関係**: $p = 1$ のときの $\Delta X = -\mathrm{F1}^2/(24\rho_b)$ は、$L'$ の補正 (§6.4) と $\rho_b = L'/(\mathrm{ANGLE} + \mathrm{K0})$ の定義で設計軌道に吸収される。$p \neq 1$ のずれの差 $\Delta X(p) - \Delta X(1) = \dfrac{\mathrm{F1}^2}{24\rho_b}\left(1 - \dfrac{1}{p}\right) = \mathrm{dxfr}\cdot\dfrac{\delta}{p}$ だけが写像に入る。

**生成関数と閉じた写像**: (5.41) の 3 つの効果は、$p_x$ と $y$ と $p$ だけの関数の生成関数で表せる。

$$
f_{\rm lin} = p_x\,\mathrm{dxfr}\left(1 - \frac1p\right) - \frac{1}{p}\left(\frac{\mathrm{dyfr}}{2}y^2 - \frac{\mathrm{dyfra}}{4}y^4\right),\quad
\mathrm{dxfr} = \frac{\mathrm{F1}^2}{24\rho_b},\ \mathrm{dyfr} = \frac{\mathrm{F1}}{6\rho_b^2},\ \mathrm{dyfra} = \frac{4\,\mathrm{dyfr}}{\mathrm{F1}^2}
\tag{5.42}
$$

($\mathrm{dyfra} = 2/(3\,\mathrm{F1}\,\rho_b^2)$)。$\Delta x = \partial f/\partial p_x = \mathrm{dxfr}\,\delta/p$、$\Delta p_y = -\partial f/\partial y = (\mathrm{dyfr} - \mathrm{dyfra}\,y^2)y/p$、$\Delta z = \partial f/\partial p = \bigl[\mathrm{dxfr}\,p_x + (\mathrm{dyfr}/2 - \mathrm{dyfra}\,y^2/4)y^2/p\bigr]/p$ が、`tbendcore` の入口の `x1 + dxfr1*dp*pinv`、`py1 + (dyfr1 - dyfra1*y1**2)*y1*pinv**2` ($q$ の単位で $p_y/p$ の更新)、`z1 + (dxfr1*px1 + (.5*dyfr1 - .25*dyfra1*y1**2)*y1**2*pinv)*pinv` ($q$ の単位) と一致する。$f$ が $x$、$p_y$、$z$ を含まないので流れは閉じた形で、**シンプレクティック性は構成から厳密**。出口は $\Delta x$ の符号が逆 ($\mathrm{dxfr2} = -\mathrm{F1}^2/(24\rho_b)$ とも書く) で、$\Delta p_y$ の符号は同じ。K64 は、これに端面の傾きによる $z$ の調整 `dzfr1 = dxfr1*sinp1` (入口)・`dzfr2 = dxfr2*sinp2` を足す。

**検証** (直線座標、$\tanh$ の端 $d = 2$ mm ($\mathrm{F1} = 6$ mm)、$\rho = 50$ m、$y = 0.5, 1$ mm の 2 点から $y$ と $y^3$ の係数を分離): $\delta = 0$ で $\Delta p_y$ の $y$ の係数は $\mathrm{F1}/(6\rho^2)$ の 1.0007 倍、$y^3$ の係数は $\tfrac23/(\mathrm{F1}\rho^2)$ の 1.031 倍 (残りは $y^5$ の項)。$\delta = 0.3$ では係数に $p$ をかけた値が $\delta = 0$ と一致 (すなわち $1/p$)。$x$ のずれ $\Delta X\cdot p/\mathrm{dxfr}$ は、直線の立ち上がりで $\delta = 0, 0.3, -0.2$ のどれも $-1.0000$、$\tanh$ では $-1.0966 = -\pi^2/9$。$\Delta X(\delta) - \Delta X(0)$ と $\mathrm{dxfr}\,\delta/p$ の比は、直線で 1.0000、$\tanh$ で 1.0966 ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (d))。

### 6.4 F1 による長さ補正と、フリンジ領域の分割

F1 があると、磁場が立ち上がる分、基準軌道の幾何を保つために、本体の円弧の長さを

$$
L' = L - \frac{(\phi_0\,\mathrm{F1})^2}{24\,L}\,\frac{\sin\bigl(\tfrac12(\phi_0(1 - \mathrm{E1} - \mathrm{E2}) - \mathrm{AE1} - \mathrm{AE2})\bigr)}{\sin(\phi_0/2)}
\tag{5.43}
$$

に取り直す (ヘルプ。K64 は入口と出口で別の F1 を許し、$((\phi_0\mathrm{fb1})^2 + (\phi_0\mathrm{fb2})^2)/(48L)$ と書く。`tffs.f:1282`、`tbende.f:384`)。同じ F1 では一致。**この式の導出は、この章では行っていない** (§13)。K0 = 0 で F1 ≠ 0 でも、設計軌道 (0, 0, 0, 0, 0, 0) が出口で 0 のまま保たれることは確認した ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (c)。$\rho_b = L'/\mathrm{ANGLE}$ にとるので、これは $L'$ の補正の物理的な正しさの確認にはならない)。

放射を入れる追跡 (`tbendr`) と CALC (`tbende`) は、本体の長さ $L'$ のうち入口と出口の $\mathrm{fb}/2$ ずつ (`f1r`・`f2r`) を、2 つの短い円弧 (長さ $\mathrm{rbl}\cdot f1r$ と $\mathrm{rbh}\cdot f1r$) に分ける (`tbendal`、`tbend.f:17-39`)。$\mathrm{rbh} = \tfrac12 + 1/\sqrt{12}$、$\mathrm{rbl} = 1/(6\,\mathrm{rbh})$ は 2 点ガウス・ルジャンドルの節 ($\tfrac12 \pm 1/(2\sqrt3)$、和 1、積 $1/6$) に一致する ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (e))。分割の目的 (放射の積分を磁場が立ち上がる領域で 2 点の求積にする、と推測) は確認していない。放射なしの追跡は、`tbendcore` を本体 $L'$ の 1 回で通す。CALC は、$\mathrm{ndiv} = 1$ でも、F1 があれば `tbendef1` (2 つの短い円弧) を使う (§9)。

## 7. `tbendcore` の読み方

追跡 (K1 = 0) の全体は `tbend` (`tbend.f:531-673`)。分岐: $\phi_0 = 0$ (ANGLE = 0) は `tsteer` (K1 = 0)、`tquad` (K0 = 0)、`tmulti` (どちらも ≠ 0)。$\phi_0 \neq 0$ で $K1 \neq 0$ は `tbendi` (§8)。K1 = 0 は、`tbshift` (DX・DY・ROTATE) → [`tbrot` (DROTATE)] → `phib = 0` なら `tbdrift` (磁場がなく、基準軌道だけ曲がる)、$L = 0$ なら `tbthin`、それ以外は `tbendcore` (放射のあるときは `tbendr` が `ndivrad` 分割で何度も呼ぶ) → 逆の取り付け誤差。

`tbendcore` (`tbend.f:372-527`) は、入口の面から出口の面までを 1 回で通す。変数は §14 の表。

| 段 | 行 | 内容 | 式 |
|---|---|---|---|
| 0 | 400-429 | F1 の係数 `dxfr1`、`dyfr1`、`dyfra1` (`fringe` のときだけ)、`dzfr1` (入口)、`dxfr2` ... (出口)。F1 の判定は `mfring > 0 or mfring == -1` (入口)、`> 0 or == -2` (出口) | (5.42) |
| 1 | 430-434 | `dcosp = cosp1 - cosp2` (桁落ちを避けて $(\sin^2\psi_2 - \sin^2\psi_1)/(\cos\psi_1 + \cos\psi_2)$) | |
| 2 | 440-448 | $p$、`rhoe`、`dpzi`、`pzi` ($= p_z/p$ で、`pxy2dpz` が $\sqrt{1 - q_x^2 - q_y^2} - 1$ を桁落ちなく) | |
| 3 | 449-453 | 入口の面までの直進と、面の座標系への回転 | (5.20)(5.21) |
| 4 | 454-457 | 線形フリンジ (入口) | (5.42) |
| 5 | 458-464 | Maxwell フリンジ (入口): `pv1sqi`、`fa`、`f`、`ff`、`x2`、`py2`、`z2` | (5.36) |
| 6 | 465-471 | 面内の運動量の大きさ `ph2` (面の座標系の $p_y$ が (5.36) で変わった後)、`dpx2`、`pz2`、`drho1`、`drho`、`t2` | §5.3 |
| 7 | 472-500 | 円弧: `dpx3a`、`dpx3`、`px3` ($= -p_{x2}/p$)、`sqwht2`、`px1px3`、`dpz3`、`dpz32`、`pz3ph2`、`pz2px1`、`t3`、`t2t3`、`t4`、`x3`、`da` | (5.22)(5.23)、$\Theta$ |
| 8 | 501-503 | `phi0a` ($= \Theta$)、`y3`、`z3` | (5.24) |
| 9 | 504-513 | Maxwell フリンジ (出口) と線形フリンジ (出口): `fa`、`f`、`ff`、`py4`、`z4`、`x4` | (5.36)(5.42) |
| 10 | 514-521 | 出口の面から出口の基準面への直進と、`px`・`x`・`y`・`py`・`z` の書き出し | (5.21) の逆 |

`t2`、`t3` ($= (p_{x1} + p_h\sin\psi_1)/(p_{z2} + p_h\cos\psi_1)$ など) は、$p_x = p_h\sin\theta$ と $p_z = p_h\cos\theta$ のとき恒等式 $(\sin\theta + \sin\psi)/(\cos\theta + \cos\psi) = \tan((\theta + \psi)/2)$ により半角の正接になる。これらと `t2t3`、`t4`、`x3`、`da` の式の各行の物理的な意味は、行ごとには追っていない。全体が、(5.12)〜(5.24) の幾何で $10^{-15}$ まで再現されることを根拠にしている (§13)。

放射を入れる場合は、`tbendr` が長さを `ndivrad` ($\mathrm{ndivmax} = 1024$) で刻み、各刻みの `tbendcore` のあとに `tradk` (光子の放射。この章の範囲外) を入れる。

## 8. K1 ≠ 0 (複合機能の偏向磁石)

K1 ≠ 0 の追跡は `tbendi` (`tbendi.f:1-200`) と `tbend.f` の `tbendiinit`・`tbendicorr`・`tbendibody`。本体を $n$ 個のスライス (長さ `alx`) に分け、各スライスで次の Strang 型の合成 ($\mathcal{C}$ は残りの非線形項、$\mathcal{B}$ は $H_2$ の厳密解 (5.8)):

$$
\mathcal{C}\!\left(\tfrac{\ell}{2}\right)\ \mathcal{B}(\ell)\ \mathcal{C}\!\left(\tfrac{\ell}{2}\right)
\tag{5.50}
$$

(隣り合う $\mathcal{C}$ は 1 つにまとめる: `tbendicorr` は `(akxn(n)+akx0)*.5`、`(alx(n)+alx0)*.5` で呼ぶ。`tbendi.f:115-121`、最後の半分は `:147`。) $\mathcal{B}$ が `tbendibody` で、(5.8)(5.11)。$\mathcal{C}$ が `tbendicorr` (`tbend.f:118-134`) で、$H_{\rm rest}$ (5.3) に当たる運動学的な非線形ドリフトと、K1 の非線形な蹴り (`pxi -= ak1*xi*xr*(0.5 - xr*(2-xr)/12)/p`、$x_r = x/\rho_0$) から成る。

- **端面**: K1 ≠ 0 の追跡は、端面を薄い蹴り (`pxi = px + tanp1*xi/rhoe`、`pyi = py - (tanp1 + fpx)*f`) で近似し、面までの直進と面から扇形の面までの弧 ($x\tan\psi$ の経路長の変化) を入れない。CALC (`tbende`) は K1 = 0 と同じ厳密な鎖 (`tbedge`) を使う。この違いは docs/30 の B117 (v2 は `DISFRIN = 0` で厳密な `tbedge` にそろえる)。DISFRIN = 0 のとき `ttfrin` (K1 の四極の非線形フリンジ) が入る (`mfring > -4`・`mfring /= 2`)。
- **分割数**: `tbendi.f:41` は $n = 1 + \mathrm{int}(|\phi_0/\varepsilon|)$、$\varepsilon = 10^{-3}\cdot\mathrm{EPS}$ (`epsbend = 1.d-3`、`tbend.f:12`)。放射があるときは `ndivrad` と大きい方 (上限 1000)。CALC (`tbende.f:455-460`) も $n = 1 + \mathrm{int}(|\phi_c/\varepsilon|)$ ($\phi_c = \phi_0\cdot r_{bc}$、F1 の分を除いた角)。**ヘルプの $1 + \lfloor\sqrt{|K1L'|/(12\cdot10^{-5}\mathrm{EPS})}\rfloor$ とは別の式** (§13・付録 A)。
- **K1 の定義**: K64 の K1 は、中央面の磁場 $B_y = B_0 + (K1/h)\ln(r/\rho_0)$ ($r = \rho_0 + x$) の勾配、すなわち曲線座標の勾配 $(1 + hx)\,\partial b/\partial x$ が一様、という暗黙の定義 (誤りではなく定義。docs/29 の A59、docs/30 の B160)。$\mathcal{C}$ の K1 の項はこの定義の $x^3$ までの級数 (この章では導出も数値の検証もしていない)。
- 水平・垂直の $H_2$ の解 (5.8)・(5.11) は、K1 の値によらず §4 で導いた。ただし `xsxkx` の符号は B250 (付録 B-4)。

pysad の対応: 追跡 `track6._tbendiinit`・`_tbendicorr`・`_tbendibody`・`_tbendi`、numba 版 `track6n.k_bendi`、CALC `lattice._bend_k1_cod_matrix`、`emit6d._tbendecorr`・`_tbendebody`、numba `perfkern_bend.bendecorr`・`bendebody`。K1 ≠ 0 の写像全体の厳密な基準との比較は、この章の確認 script には含めていない。

## 9. 行列 (CALC・Emittance[]) と追跡の関係

CALC の `tbende` (`tbende.f:351-561`) は、追跡と**同じ段**の各写像を COD ($\xi$) のまわりで微分して掛ける (式 (4.34) と同じ構造)。段ごとの解析行列 (`trans1`) は、K1 = 0 では次のサブルーチンに分かれている。

| 段 | 追跡 | CALC |
|---|---|---|
| 取り付け誤差 | `tbshift`・`tbrot` | `tchge` (`tbrote`) |
| F1 (入口) | `tbendcore` の中 (面までの直進の**後**) | `tblfre` (`tbende.f:650-678`)。**面までの直進の前** |
| 面までの直進・Maxwell・面から扇形の面 | `tbendcore` の中 | `tbedge` (`tbfrie.f:88-125`): `tbedgedrift` → `tbedgemaxwell` → `tbedgebody` (出口は逆順) |
| F1 の領域の本体 | 放射があるときだけ `tbendal` で 2 分割 | **常に** `tbendef1` (2 つの短い円弧。`tbende.f:319-347`、`ef1`)。F1 = 0 なら無い |
| 本体 | `tbendcore` (円弧の厳密解) | `tbendebody0` (`tbende.f:223-317`): 円弧を COD で微分。$\mathrm{ndiv} = 1$ (K1 = 0)、放射があるとき `ndivrad` |
| 分割 | 放射のあるとき `tbendr` | 放射のあるとき `nrad`、`emidiv` |
| $z$ の行 | 写像そのもの | シンプレクティック条件で組む ($T_{5j} = -T_{1j}T_{26} + T_{2j}T_{16} - \cdots$。`tbendebody0` の `trans1(5,1)` など、(4.38) と同じ形) |

**一致の程度** (電子 3 GeV、$L = 1$ m、ANGLE = 0.1 rad。追跡の中心差分ヤコビアンと、CALC の `lattice._el_cod_step` の 4×4 の差。[bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (a)): 設計軌道で $5 \times 10^{-13}$ (扇形・矩形・非対称とも)。小さいずれ ($\delta = 0.002$、0.3 mm、0.5 mrad) で $10^{-12}$ 〜 $3 \times 10^{-11}$。大きいずれ ($y = 10$ mm、$\delta = 0.05$) で、扇形 $2.2 \times 10^{-8}$、矩形 $3.0 \times 10^{-10}$、非対称 $5.7 \times 10^{-10}$。COD の更新 (出口の $x, p_x, y, p_y, z$) の差は小さいずれで $6 \times 10^{-11}$、大きいずれで $2.0 \times 10^{-8}$。F1 = 0.1 (FRINGE = 1、E1 = E2 = 0.5) では 4×4 が小さいずれで $10^{-7}$、大きいずれで $4.4 \times 10^{-5}$。

**差の出どころ**: (i) CALC の Maxwell フリンジの COD の更新が因子 $p$ を落としている (付録 B-1)。(ii) F1 ≠ 0 のとき、E = 0 でも、CALC の写像 (COD の更新) の微分が追跡のヤコビアンと食い違う (小さいずれ $x$ = 1 mm、$y$ = 0.3 mm、$\delta$ = 0.002、F1 = 0.1 で、$x$ の行の $y$ の列が $4 \times 10^{-8}$、$p_y$ の列が $1 \times 10^{-7}$。F1 = 0 では $10^{-12}$)。CALC の解析行列は自身の COD 更新の数値微分と $10^{-8}$ で一致するので、写像そのものが違う ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (f))。CALC は F1 の位置が追跡と違い (面までの直進の前。E = 0 では直進が無いので同じ位置)、F1 の領域を常に 2 つの円弧 (`tbendef1`) に分けるが、これらが原因かどうかは確認していない。物理的にどちらが正しいかも確認していない。(iii) 段ごとの解析行列は、同じ段の COD 更新の数値微分と、`tbendebody0` で $1.4 \times 10^{-9}$、`tbedgedrift` で $4 \times 10^{-10}$、`tbedgebody` で $6 \times 10^{-12}$ だが、**`tbedgemaxwell` では $4.5 \times 10^{-5}$** ($y = 10$ mm、$\delta = 0.05$): 解析行列の $\partial x/\partial y$ は COD 更新の数値微分の $1.05 = 1 + \delta$ 倍 ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (c2))。解析行列の方が (5.36) の $p^2$ の式と整合し、COD の更新が不整合、という内部の矛盾 (付録 B-1)。

**シンプレクティック性**: 追跡のヤコビアン (数値微分) は $10^{-11}$ 台 (数値微分の精度)、$z$ の行も上の条件と $10^{-11}$ で一致。CALC の 6×6 は (F1 = 0 で) $5 \times 10^{-17}$ ($z$ の行を条件で組むので)。F1 つきでは $4.4 \times 10^{-5}$ ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (d))。

## 10. Fortran の対応表

行番号は参照した K64 のソースのもの。

| 式・内容 | ファイル:行 | サブルーチン | 変数・メモ |
|---|---|---|---|
| 追跡の全体 (K1 = 0) | `tbend.f:531-673` | `tbend` | `phib = phi0 + ak0` (565)、`rhob`・`rho0`・`drhob` (612-614)、`ndiv` (619、640)、`tbendcore` (655)・`tbendr` (662) |
| 円弧と端面の 1 回の通過 (5.20)〜(5.24)、(5.36)、(5.42) | `tbend.f:372-527` | `tbendcore` | §7 の表 |
| 放射の分割 | `tbend.f:305-370` | `tbendr` | `aln = (al - f1r - f2r)/ndiv`、`tbendal` で長さを決める |
| F1 領域の分割の長さ | `tbend.f:17-39`、`:12-13` | `tbendal`、`rbh`・`rbl` | |
| 取り付け誤差 | `tbend.f:248-303`、`:168-246` | `tbshift`、`tbrot` | 章の範囲外 |
| 磁場のない BEND | `tbend.f:675-726` | `tbdrift` | `phib = 0` |
| 薄い BEND | `tbend.f:728-739` | `tbthin` | `al = 0` |
| $H_2$ の係数 (5.5)・(5.8) | `tbend.f:41-116` | `tbendiinit` | `akxsq`、`akysq`、`drhop`、`dcx`、`sxkx`、`dcxkx`、`xsxkx`、`aksx`、`dcy`、`syky`、`aksy` |
| $H_2$ の本体 (5.8)・(5.11) | `tbend.f:136-160` | `tbendibody` | `xiksq`、`dxf`、`dpxf`、`dyf`、`dpyf`、`hi`、`zi` |
| 残りの非線形項 (5.50) | `tbend.f:118-134` | `tbendicorr` | |
| K1 ≠ 0 の追跡 | `tbendi.f:1-200` | `tbendi` | 端面は薄い蹴り (`tanp1`、`tanp2`)、`ndiv = 1+int(abs(phi0/eps))` (41) |
| L' の補正 (5.43) | `tffs.f:1282-1285` (追跡)、`tbende.f:381-389` (CALC)、`tfgeo.f:229` (幾何) | `tpara` ほか | `al = al - ((phi*fb1)**2 + (phi*fb2)**2)/al/48 * sin(.5 w)/sin(.5 phi)` |
| CALC の全体 | `tbende.f:351-561` | `tbende` | K1 = 0 は `ndiv = 1` (460)、`tbedge` (472、548)、`tblfre` (441、553)、`tbendebody0`・`tbendef1` (494-509) |
| CALC の円弧 | `tbende.f:223-317` | `tbendebody0` | `dpxf`、`xf`、`da`、`trans1(1..5, ...)` |
| CALC の F1 領域の円弧 | `tbende.f:319-347` | `tbendef1` | `rb = rbl`、`1/6/rb` の 2 回 |
| CALC の端面 | `tbfrie.f:88-125` | `tbedge` | 入口 `tbedgedrift` → `tbedgemaxwell` → `tbedgebody`、出口は逆順 |
| 面までの直進 (5.20) | `tbfrie.f:231-289` | `tbedgedrift` | |
| 面から扇形の面へ | `tbfrie.f:127-188` | `tbedgebody` | |
| Maxwell フリンジ (5.36) | `tbfrie.f:190-229` | `tbedgemaxwell` | `a`・`af`・`b`・`bf`、`cod(1)`・`cod(4)`・`cod(5)` の更新は 209-211 行 |
| F1 の線形フリンジ (5.42) | `tbende.f:650-678` | `tblfre` | `dxfr`、`dyfr`、`dyfra` |
| ステアリング (ANGLE = 0) | `tsteee.f`、`tsteer.f` | `tsteee` ほか | 章の範囲外 |
| 薄い偏向 (CALC) | `tbende.f:680-704` | `tbthie` | |

## 11. pysad の対応表

| 内容 | pysad の関数 (ファイル:行) | 対応する Fortran | 注意 |
|---|---|---|---|
| 追跡の K1 = 0 の全体 | `track6._tbend_ens` (`track6.py:1614`)、`_tbend0_ens` (`:1729`) | `tbend` | `el.E1`・`el.E2` は角 (rad)。`psi1 = el.E1` |
| 円弧と端面 (python) | `track6._tbendcore` (`:963`) | `tbendcore` | 変数名は Fortran と同じ。`bp = (rhob, rho0, drhob, fb1, fb2)` |
| 円弧と端面 (numba) | `track6n.tbendcore1` (`track6n.py:738`)、`k_bend0` (`:939`)、`tbendal` (`:926`) | 同上 | |
| 取り付け誤差 | `track6._tbshift` (`:1087`)、`_tbrot` (`:1129`) | `tbshift`・`tbrot` | |
| F1 領域の分割 | `track6._tbendal` (`:1427`) | `tbendal` | |
| K1 ≠ 0 | `track6._tbendiinit` (`:1185`)、`_tbendicorr` (`:1235`)、`_tbendibody` (`:1254`)、`_tbendi` (`:1273`)、`track6n.k_bendi` (`:1098`) | `tbendiinit` ほか | `xsxkx` の符号は `compat.v2()` で分かれる (B250) |
| CALC の COD の更新 | `lattice._bend_cod_update` (`lattice.py:621`)、`_bend_fringe_params` (`:596`)、`_el_cod_step` (`:2616`) | `tbende` | `fringe`・`mfring`・`fb1`・`fb2`。F1 = 0 は `ndiv = 1` |
| CALC の端面の COD の更新 | `lattice._tbedgedrift_cod` (`:456`)、`_tbedgebody_cod` (`:476`)、`_tbedgemaxwell_cod` (`:512`)、`_tbedge_cod` (`:533`) | `tbedgedrift` ほか | |
| CALC の円弧と F1 | `lattice._tbendebody0_cod` (`:552`)、`_tblfre_cod` (`:606`) | `tbendebody0`・`tblfre` | |
| CALC の解析行列 (numpy) | `lattice._tbedgemaxwell_mat` (`:692`)、`_tbedgedrift_mat` (`:718`)、`_tbedgebody_mat` (`:756`)、`_tbendebody0_mat` (`:809`)、`_bend_edge_chain` (`:875`)、`_bend_cod_matrix` (`:894`) | 各 `trans` | |
| CALC (numba) | `elkern._tbedgedrift` (`elkern.py:798`)、`_tbedgemaxwell` (`:847`)、`_tbedgebody` (`:882`)、`_tbendebody0` (`:958`)、`_bend_core` (`:1756`)、`k_bend` (`:1838`) | 同上 | |
| Emittance[] | `emit6d._tbende` (`emit6d.py:1268`)、`_tbende_core` (`:1321`)、`_tbendecorr` (`:1160`)、`_tbendebody` (`:1167`)、`perfkern_bend.bendecorr` (`:65`)、`bendebody` (`:138`) | `tbende` ほか | |
| 設計軌道の線形行列 | `elements.Bend.transfer_matrix` (`elements.py:234`) | | 別系統。K0 ≠ 0 と ANGLE ≠ 0 の扱いに注意 |
| 3 次のリー生成関数の取り出し | `bendnl.bend_decomposition` (`bendnl.py:140`) | | 追跡の写像から。この章の範囲外 |

## 12. 確認スクリプトと出力

すべて、リポジトリのルートから `py -3.12 physics_guide/scripts/bend_0x_*.py` で動く (数秒)。出力は標準出力と `physics_guide/outputs/bend_0x_*.out.txt`。(a)〜(g) は各スクリプトの節の記号。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `bend_01_h2_body.py` | (a) $A_s$ から磁場、$H$ の 0〜3 次の展開 (sympy)。(b) 偏微分の表と運動方程式、K64 の `akxsq`・`drhop`・`akysq` との一致 (sympy で差 0)。(c) 水平の解が運動方程式を満たす、分散 $R_{16}$・$R_{26}$・$R_{56}$ (sympy)。(d) K64 の `tbendibody` の式と $H_2$ の ODE の比較 ($x, q, y, q_y$ は $10^{-18}$)、`xsxkx` の符号。(e) pysad の `_tbendiinit`・`_tbendibody` (k64 と v2)。(f) $H_2$ と厳密な $H$ の差の振幅依存 (比 4)。(g) 設計軌道の行列 | [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) |
| `bend_02_arc_faces.py` | (a) 面までの直進のシンプレクティック性 (sympy、40 桁)。(b) 実験室系の幾何 (面・円・フリンジ) と追跡の 5 成分が $10^{-15}$ で一致 (扇形・矩形・非対称・AE・K0 ≠ 0、$\delta = 0, 0.02, -0.03$)。(c) 設計軌道。(d) ヘルプの本体の式 (5.22)〜(5.24) と幾何。(e) DISFRIN = 1 でも Maxwell フリンジが入る。(f) 曲線座標の $H$ の ODE と幾何 | [bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) |
| `bend_03_fringe.py` | (a) 端の磁場 (sympy)。(b) 硬い端の Maxwell フリンジ (5.35) とローレンツ力の ODE。(c) CALC `tbedgemaxwell` と追跡の式のシンプレクティック性。(d) F1 の $y$・$y^3$・$x$ の項 (ODE)、F1 の定義と $\int g'^2$。(e) `rbh`・`rbl`、DISFRIN | [bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) |
| `bend_04_matrix_vs_map.py` | (a) CALC の 4×4・δ の列と COD の更新を追跡と比較。(c)(c2) CALC の解析行列と COD 更新の数値微分、段ごとの差。(d) シンプレクティック性。(e) 設計軌道の 2×2 と端面の行列。(f) F1 ≠ 0 の CALC と追跡の差 | [bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) |

主な出力 (抜粋)。K64 の `tbendcore` と実験室系の幾何の差 ($bend\_02$ (b)、最大):

```
扇形 (E1 = E2 = 0)                    δ=+0.00   6.3e-16
矩形 (E1 = E2 = 0.5)                  δ=+0.02   6.7e-16
非対称 (E1 = 0.3、E2 = 0.7)            δ=+0.00   9.1e-16
K0 = 0.005 (ρb ≠ ρ0)、E1 = E2 = 0.5   δ=+0.00   9.9e-16
最大 1.1e-15
```

Maxwell フリンジの $\Delta x$ ($bend\_03$ (b)。$d \to 0$ の外挿):

```
δ = 0.0: A = 4.0608e-08   閉じた式 (p^2) 4.0608e-08 (比 1.0000)、p の因子だけ 4.0608e-08
δ = 0.3: A = 3.1237e-08   閉じた式 (p^2) 3.1237e-08 (比 1.0000)、p の因子だけ 2.4028e-08 (比 1.3000)
```

## 13. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV ($\gamma_0 = 5870.85$)、$L = 1$ m、ANGLE = 0.1 rad ($\rho = 10$ m) が中心 ($H_2$ の例は $\rho = 10$ m、$L = 1.3$ m)。Maxwell フリンジと F1 の ODE は、直線座標、$1/\rho = 0.02$ /m ($\rho = 50$ m)、端の幅 $d = 1, 2$ mm、$y = 0.5$〜2 mm の $y$・$\delta = 0$・0.3 の範囲。放射、スピン、ソレノイドは含まない。

**この章の範囲外 (別の章)**: 取り付け誤差 (`tbshift`・`tbrot`・DROTATE の中点まわりの回転)、放射 (`tradk`・`ndivrad` の意味)、K1 ≠ 0 の複合機能の写像全体の厳密な基準との比較、ANGLE = 0 のステアリング (`tsteer`・`tsteee`)、端部の一般論 (CAVI・SOL・MULT の端)。

**未確認・限界**:
- **Maxwell フリンジの $y^4$ の係数** ($a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$): 1 次の導出 (§6.2) では出ない。ヘルプ・K64 の係数を使っただけで、独立な導出も数値の検証もしていない (検証した磁場では $y^2/\rho^2$ が $10^{-9}$ 以下で見えない)。また (5.34) は $p_y$ を入れない $p_z$ を使う近似で、$p_y \neq 0$ の高次の効果は確認していない。
- **F1 の式**: 1 つの F1 から $y$・$y^3$・$x$ の 3 つの効果を作るのは、直線の立ち上がりで厳密、$\tanh$ では $y^3$ の項が同じ形で $x$ の項が $\pi^2/9$ 倍、誤差関数では $y^3$ の項が -4.5% (§6.3)。K64 の F1 は端の形を 1 つの長さで表した近似。
- **$L'$ の補正 (5.43)** の導出は行っていない。設計軌道が保たれる ($\rho_b$ を $L'$ で決めるので当然) ことしか確認していない。
- **F1 の領域の 2 分割** (`rbh`・`rbl`) が 2 点ガウス・ルジャンドルの節であることは確認した。なぜその分割がよいか (放射の積分の求積と推測) は確認していない。
- **`tbendcore` の本体の各行** (`t2`、`t3`、`t2t3`、`t4`、`x3`、`da`) の式の物理的な意味は、行ごとには追っていない。全体の結果は、実験室系の幾何と $10^{-15}$ で一致する。
- **K1 ≠ 0**: `tbendicorr` の式の導出、K1 の非線形な蹴りの級数 (A59)、端面の薄い蹴りと厳密な端の差 (B117)、写像全体の厳密な基準との比較は、この章では行っていない。分割数の式 (§8) はコードを読んだもので、gs での実測はしていない。
- **CALC と追跡の違い** (§9): F1 ≠ 0 のとき、写像のヤコビアンが $10^{-7}$ の桁で違う。原因は特定していない。物理的にどちらが正しいかも確認していない。
- **ヘルプと Fortran の Maxwell フリンジの式**: ヘルプの $x_2$ の式と `tbendcore` の `ff` は代数的に等しい (§6.2)。$p_{y2}$ は、ヘルプが $1/p_1$ を余分に持つ (付録 A)。`f`・`fa`・`ff` の対称な書き方にした意図はコードから確認できない。
- **走らせていない検証**: gs (SAD 本体) との比較、K1 ≠ 0 の本体の厳密な ODE、放射つき、DROTATE つき、F1 つきの全体を (5.31) の端の磁場の ODE で直接通すこと。

## 14. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。下の「式」は数学的な定義。定義の場所は実際のコードで確認した。

| Fortran の名前 | 式 | 定義の場所 (Fortran) | pysad での名前 |
|---|---|---|---|
| `phi0` | ANGLE ($\phi_0$) | `tbend` の引数 | `phi0` |
| `ak0`・`phib` | K0、$\phi_b = \phi_0 + \mathrm{K0}$ | `tbend.f:565` | `ak0`・`phib` |
| `al` | $L'$ (F1 の補正後の長さ) | `tffs.f:1283` | `al` |
| `rho0`・`rhob`・`drhob`・`rhoe` | $\rho_0 = L'/\phi_0$、$\rho_b = L'/\phi_b$、$\rho_b - \rho_0$、$\rho_bp$ | `tbend.f:612-614`、`tbendcore` の局所変数 | 同名 |
| `psi1`・`psi2`・`cosp1`・`sinp1`・`cosp2`・`sinp2` | 端面の傾き $\psi_{1,2}$ (rad) と三角関数 | `tffs.f:1260-1262`、`tbend` の引数 | 同名 (pysad の `el.E1` が $\psi_1$) |
| `w`・`cosw`・`sinw`・`sqwh`・`sinwp1` | $\omega$、$\cos\omega$、$\sin\omega$、$1 - \cos\omega$ (`xsincos` の `dc` の符号反転)、$\sin(\phi_0 - \psi_2) = \sin(\omega + \psi_1)$ | `tbendr` (`tbend.f:356-360`)、`tffs.f:1280` | `cosw` など (`track6._tbend0_ens`) |
| `fb1`・`fb2` | $\mathrm{F1} + \mathrm{FB1}$、$\mathrm{F1} + \mathrm{FB2}$ | `tffs.f:1260-1262` | `fb1`・`fb2` |
| `mfring` | FRINGE の内部の値 (`FRMD`)。`FRMD ≤ 0` では `-1` 以外は `fb1 = 0`、`-2` 以外は `fb2 = 0` | `tffs.f:1271-1277` | `mfring` |
| `fringe` | DISFRIN = 0 のとき真。`tbendcore` では `dyfra` だけを切り替える | `tturn.f:342` (`ky_FRIN_BEND` = DISFRIN。`initb1.f:925`) | `fringe` |
| `dxfr1`・`dyfr1`・`dyfra1`・`dzfr1` | $\mathrm{F1}^2/(24\rho_b)$、$\mathrm{F1}/(6\rho_b^2)$、$4\,\mathrm{dyfr}/\mathrm{F1}^2$、$\mathrm{dxfr}\sin\psi_1$ | `tbend.f:400-429` | 同名 |
| `f1r`・`f2r` | $\mathrm{fb}/2$ (放射のとき本体から切り出す長さ) | `tbend.f:633-637` | 同名 |
| `rbh`・`rbl` | $\tfrac12 + 1/\sqrt{12}$、$1/(6\,\mathrm{rbh})$ | `tbend.f:12-13` | `_RBH`・`_RBL` (`track6.py:1423`、`lattice.py`) |
| `pv1sqi` | $1/\max(s_{\min}, 1 - q_{x1}^2)$ ($s_{\min} = 10^{-4}$)。(5.36) の $p^2/p_v^2$ | `tbend.f:458` | `pv1sqi` |
| `fa`・`f`・`ff` | $y/(\rho_e\sqrt{1 - q_x^2})$、$a\cdot\mathtt{fa}$、$\tfrac14(f + fa)\,y\cdot\mathtt{pv1sqi}$ ($\Delta x$) | `tbend.f:459-461`、`:505-507` | 同名 |
| `dpzi`・`pzi` | $p_z/p - 1$、$p_z/p$ (`pxy2dpz`) | `tbend.f:447-448` | `dpzi`・`pzi` |
| `ph2`・`dph2` | $p_h/p = \sqrt{1 - q_y^2}$ と、1 を引いたもの | `tbend.f:465-466` | `ph2`・`dph2` |
| `dpx2`・`pz2`・`drho1`・`drho`・`t2` | §5.3 | `tbend.f:467-471` | 同名 |
| `dpx3a`・`dpx3`・`px3`・`t3`・`da`・`phi0a` | (5.22) の $q$ 単位の各部分 ($\texttt{px3} = -p_{x2}/p$)、$\Theta - \phi_0$、$\Theta$ | `tbend.f:472-501` | 同名 |
| `akxsq`・`akysq`・`drhop` | (5.5) | `tbend.f:59-62` | `m.akxsq` など |
| `dcx`・`sxkx`・`dcxkx`・`xsxkx`・`aksx` | $\cos\varphi_x - 1$、$\sin\varphi_x/\omega_x$、$\texttt{dcx}/\texttt{akxsq}$、$(\ell - \texttt{sxkx})/\texttt{akxsq}$、$-\omega_x\sin\varphi_x$ ($x$ の $q$ の $x$ 係数) | `tbend.f:63-92` | `m.dcx` など |
| `dcy`・`syky`・`aksy` | 垂直の同様の量 | `tbend.f:93-114` | `m.dcy` など |
| `hi` | $\tfrac12(q_x^2 + q_y^2 - \texttt{akxsq}\,x^2 - \texttt{akysq}\,y^2) - \texttt{drhop}\,x$ | `tbend.f:149` | `hi` |
| `epsbend` | $10^{-3}$ (K1 ≠ 0 の分割の角) | `tbend.f:12` | 変数名なし (`track6._tbendi` の `eps = 1e-3 * eps0`、`track6.py:1287`) |
| `ndiv` | スライスの数 | `tbendi.f:41`、`tbend.f:640`、`tbende.f:460` | `ndiv` |

確認: `tbend.f` の `ky_FRIN_BEND` が DISFRIN に対応することは `initb1.f:925` (`kytbl(kwFRIN,0)=sethtb('DISFRIN ',...)`) で確認した。

## 付録 A. ヘルプの記述・K64 との違い (偏向磁石)

| 箇所 | 内容 |
|---|---|
| ヘルプ BEND > transformation、body | $s \equiv \mathrm{ANGLE}\times\rho_b\times(\arcsin\ldots)$ の $\mathrm{ANGLE}\times$ は誤記と読む。面内の経路長は $s = \rho_bp_h\Theta$ ($p_h = \sqrt{p^2 - p_y^2}$)。そう読めば $y_2 = y_1 + s\,p_y/\sqrt{p_1^2 - p_{y1}^2}$ が (5.24) と一致 (§5.3)。$p_{x2}$・$x_2$ の式は (5.22)(5.23) と同じ |
| ヘルプ BEND、nonlinear fringe の $p_{y2}$ | $p_{x1}\,y_1(1 - y_1^2/(6\rho_b^2))/(p_1\rho_b\sqrt{p_1^2 - p_{x1}^2})$ は、$p_1$ を 1 つ余分に持つ。正準運動量での導出 (§6.2) と、K64 の式 (`f` の $\rho_e = \rho_bp$ は $q$ の単位) は $-(p_x/p_v)\,y/\rho_b$。ヘルプの式の $p_1$ が (q の単位の) 規格化を指すなら一致する (ヘルプは運動量の規格化を明記していない) |
| ヘルプ BEND、linear fringe の $p_{y2}$ | $y_1(\Delta y_{fr} - \Delta y_{fra}y_1^2)/p_1^2$ の $1/p_1^2$ は、$q$ の単位 ($p_y/p$) の更新と読めば K64 と一致。正準運動量では $1/p$ (§6.3) |
| ヘルプ BEND の F1 | 厳密な定義は $6\int(g - g^2)\,ds$ (5.40)。四極の F1 (04 章の $I_1$ のモーメント) と**定義が違う**。偏向磁石の F1 は $y$ の収束の量で、$x$ のずれの $\mathrm{F1}^2/24$ は直線の立ち上がりでだけ厳密 (§6.3) |
| ヘルプ BEND、DISFRIN | 「DISFRIN が 0 でないと非線形の Maxwell フリンジを止める」。K1 = 0 の追跡 (`tbendcore`) では、DISFRIN は F1 の $y^3$ の項 (`dyfra`) だけを止め、硬い端の Maxwell フリンジ ((5.36) の `f`・`ff`) は止まらない ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (e)、[bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (e))。CALC も DISFRIN によらず Maxwell を入れる (docs/30 の B134)。K1 ≠ 0 の追跡は `ttfrin` (K1 の端) が DISFRIN に従う |
| ヘルプ BEND の K1 ≠ 0 の分割数 | $1 + \lfloor\sqrt{|K1L'|/(12\cdot10^{-5}\mathrm{EPS})}\rfloor$。K64 のコードは、追跡 (`tbendi.f:41`)・CALC (`tbende.f:460`) とも $1 + \mathrm{int}(|\phi/\varepsilon|)$、$\varepsilon = 10^{-3}\mathrm{EPS}$ (角で刻む)。放射があるときは `ndivrad` との大きい方 (gs での実測はしていない) |
| ヘルプ BEND、K1 ≠ 0 の端面 | 「E1・E2 の効果は薄い四極で近似する」。K64 のコードは、端面の薄い蹴りは偏向の蹴り (`tanp1*xi/rhoe`) と、垂直の Maxwell 型の蹴り (`(tanp1 + fpx)*f`) で、四極そのものではない (`tbendi.f:98-108`)。CALC は厳密な `tbedge` (B117) |
| ヘルプに無いもの | (1) 円弧の厳密解の導出 (§5)。(2) $\omega = \mathrm{ANGLE} - \psi_1 - \psi_2$ の意味と、面の座標系 (§5.3)。(3) Maxwell フリンジの導出 (§6.2)。(4) F1 の 3 つの係数の導出 (§6.3)。(5) $L'$ の導出 (この章でも未)。(6) K1 ≠ 0 の本体の $H_2$ と、`tbendicorr`・`tbendibody` の式 (§4・§8)。(7) 行列と追跡の関係 (§9) |

## 付録 B. K64 の誤りの疑い (別枠。確定には SAD (gs) での確認が要る)

この章は説明だけで、v2 (既定) が直す方針の動作は変えていない。ここに挙げるのは、この章の導出・数値の確認で見つけた、K64 のソースの食い違いの**疑い**。

**B-1. CALC の `tbedgemaxwell` の COD の更新が $\Delta x$ と $\Delta z$ で因子 $p = 1 + \delta$ を落としている (新しく見つけた疑い。docs/30 に無い)**

- **内容**: `tbfrie.f:209` の `cod(1)=cod(1)+bf*pr` と `:211` の `cod(5)=cod(5)-bf*pxi` は、(5.36) の $\Delta x = b\,y^2p^2/(\rho_bp_v^3)$、$\Delta z = -b\,y^2p_xp/(\rho_bp_v^3)$ より因子 $p$ だけ小さい ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$ に `pr` を 1 回しかかけない。(5.36) は `pr**2`、$z$ は `pxi*pr`)。$\Delta p_y$ は正しい。
- **根拠**: (i) 端の磁場 (5.31) の中のローレンツ力の積分で、$\Delta x$ は $\delta = 0.3$ で $p^2$ の式と一致 (比 1.0000)、`pr` の式とは 1.3000 倍違う ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。(ii) 追跡 (`tbendcore`) の式は $p^2$ の形で、6×6 ヤコビアンがシンプレクティック ($10^{-12}$)。CALC の更新はシンプレクティックでない ($\delta = 0.1$ で $1.7 \times 10^{-4}$。$\delta = 0$ でも、$\delta$ の列の微分がずれて $2.0 \times 10^{-5}$) ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (c))。(iii) CALC の解析行列 `trans(1,3) = af*pr**2/pvi2` は $p^2$ の式の微分で、COD の更新の数値微分と $1 + \delta$ 倍違う (内部の矛盾。[bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (c2))。
- **大きさ**: ずれは $\delta y^2/(2\rho)$ 程度。$y = 10$ mm、$\delta = 0.05$、$\rho = 10$ m の扇形で、CALC の COD の更新と追跡の写像が $2.0 \times 10^{-8}$ ($x$ と $z$) 違う。
- **状態**: 疑い (ソースの読みと数値)。pysad は k64・v2 とも同じ ($\mathtt{lattice.\_tbedgemaxwell\_cod}$ は Fortran の直訳)。直していない。

**B-2. DISFRIN が K1 = 0 の追跡の Maxwell フリンジを止めない (ヘルプの記述との食い違い。仕様の疑い)**: 付録 A。CALC と追跡の両方で止まらないので、両者の間では一貫している (B134 と同じ種類)。

**B-3. ヘルプの式の記述の食い違い (K64 のソースの誤りではなく、ヘルプの誤記・規格化の省略の疑い)**: 付録 A の $s$ の $\mathrm{ANGLE}\times$、nonlinear fringe の $p_{y2}$ の $1/p_1$。

**B-4. K1 ≠ 0 の本体の `xsxkx` の符号 (docs/30 の B250 の既知の誤り)**: (5.11) の導出から、`akxsq < 0` の枝の `xsxkx = -xspx/akx/akxsq` は符号が逆で、連続な式 $(\ell - \texttt{sxkx})/\texttt{akxsq}$ ($= \texttt{xspx}/\texttt{akx}/\texttt{akxsq}$) が正しい。$H_2$ の ODE との比で $z$ の誤差 $7.3 \times 10^{-5}$ m ($\rho_0 = \rho_b = 10$ m、$L = 1.3$ m、$\delta = 0.01$、K1 ≥ 0) が、符号を直すと $10^{-18}$ ([bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (d)(e))。v2 は直してある。
