# 05. 偏向磁石 (BEND)

**この章で分かること**: 偏向磁石 (BEND) の写像が、(1) 入口の端面 (面までの直進)、(2) 線形フリンジ (F1)、(3) 非線形の Maxwell フリンジ、(4) 本体 (円弧)、(5) 出口の側の同じ段、から成ることを、ハミルトニアンと一様磁場の幾何から途中の式を省かずに導き、K64 の Fortran (追跡の `tbendcore`、CALC の `tbende`) と pysad の関数に対応させる。ヘルプは各段の結果の式だけを書いている。その途中 (円弧の厳密解、端面の幾何、フリンジの導出、行列と追跡の関係) を埋める。式番号 $(5.n)$ から、§10 の Fortran の対応表と §11 の pysad の対応表を引ける。

**前提**: [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) (特に §9.4 の偏向磁石の $A_s$ と運動方程式)、[01_variables_and_conventions.md](01_variables_and_conventions.md) の §1〜§2 と §8 (SAD の変数名の注釈の規約)、[04_quad.md](04_quad.md) (章の書き方。非線形フリンジの母関数、行列と追跡の関係)。解析力学と加速器の線形光学。SAD の内部は知らなくてよい。

**読む順**: §1 (結論) → §2 (記号) → §3〜§4 (H と、K1 つきの本体の H2 の解) → §5 (K1 = 0 の円弧と端面の厳密解) → §6 (フリンジ) → §7 (`tbendcore` の読み方) → §8 (K1 ≠ 0) → §9 (行列と追跡)。コードを読みながら使うときは §10・§11・§14 から入る。K64 の誤りの疑いは付録 B に別にまとめた。

## 要約と目次

**この章で分かること**

- **何を導くか**: BEND の追跡の写像を、(入口の面までの直進) → (線形フリンジ F1) → (非線形の Maxwell フリンジ) → (円弧の本体) → (出口側の同じ段) の合成として、ハミルトニアンと一様磁場の幾何から導く。K1 ≠ 0 (複合機能) の本体と、CALC の行列と追跡の関係も扱う。
- **主な結果**: $H_2$ の運動方程式 (5.4)(5.5) と $z$ の更新 (5.11)、円弧と面の幾何 (5.20)〜(5.24)、端面角の効果の行列 (5.25)、硬い端の Maxwell フリンジ (5.35)(5.36)、F1 の効果 (5.41)(5.42)、F1 による長さ補正 $L'$ (5.43)(5.43c)、K1 ≠ 0 の分割 (5.50)。
- **対応**: 追跡は Fortran の `tbend`・`tbendcore` (K1 ≠ 0 は `tbendi`)、CALC は `tbende`・`tbedge`・`tblfre`。pysad は `track6._tbendcore`・`lattice._bend_cod_update`・`emit6d._tbende` など (§10・§11 の表)。
- **確かめ方**: sympy・mpmath (30〜40 桁)・DOP853 の数値積分・端の磁場の中のローレンツ力の積分・追跡のヤコビアンとの比較・K64 のソースの読み。一部 (付録 B-1) は gs の閉軌道で確かめた (§12 の `bend_01`〜`bend_09`)。
- **K64・ヘルプの誤りの確定・疑い** (付録 B と一致): **確定** = B-1 (CALC の `tbedgemaxwell` の COD の更新が因子 $p$ を落とす。gs で確認。docs/30 の B324。v2 は修正)、B-4 (K1 ≠ 0 の本体の `xsxkx` の符号。B250。v2 は修正)。**疑い (未決)** = B-2 (DISFRIN が K1 = 0 の追跡の Maxwell フリンジを止めない。仕様の疑い)、B-3 (ヘルプの式の誤記・規格化の省略の疑い)。取り下げは無い。

**導出の流れ**

1. [§1 結論](#1-結論)・[§2 記号と単位](#2-記号と単位): 結果の一覧と、SAD・K64 の Fortran・pysad の記号の対応。
2. §3 ハミルトニアンと $H_2$: 曲線座標のハミルトニアンを 2 次で切る。
   - [§3.1 磁場とベクトルポテンシャル](#31-磁場とベクトルポテンシャル): $A_s$ (5.1)。
   - §3.2 厳密な $H$ と $H_2$: (5.2) を次数ごとに展開して $H_2$ と $H_{\rm rest}$ (5.3) に分ける。
   - [§3.3 偏微分と運動方程式](#33-偏微分と運動方程式): (5.4) → (5.5)。K64 の `akxsq`・`akysq`・`drhop` との対応。
3. §4 $H_2$ の厳密解: 調和振動子の解 (5.6)〜(5.8) と分散 (5.9)。
   - [§4.1 水平と垂直](#41-水平と垂直)、§4.2 $z$ の更新: (5.10) の積分から (5.11a)〜(5.11c) を経て (5.11) を得る。K64 の `xsxkx` の符号の誤り (B-4) がここで分かる。
4. [§5 K1 = 0 の円弧の厳密解と端面](#5-k1--0-の円弧の厳密解と端面): K1 = 0 では $H_2$ の近似でなく、円と面の幾何で解く。
   - [§5.1 一様磁場の円](#51-一様磁場の円) (5.12)、[§5.2 入口の面までの直進](#52-入口の面までの直進) (5.20)(5.21)(5.21a)、[§5.3 面から面への円弧](#53-面から面への円弧) (5.22)〜(5.24)、[§5.4 端面角の効果 (線形光学)](#54-端面角の効果-線形光学) (5.25)。
5. [§6 端部 (フリンジ)](#6-端部-フリンジ): 本体の外側の立ち上がりを、Maxwell 方程式を満たす場から導く。
   - [§6.1 端の磁場](#61-端の磁場) (5.30)〜(5.32)、[§6.2 硬い端の Maxwell フリンジ](#62-硬い端の-maxwell-フリンジ) (5.33)〜(5.36) (生成関数 (5.34)。$y^4$ の係数は導けない)。
   - [§6.3 線形フリンジ (F1)](#63-線形フリンジ-f1) (5.40)〜(5.42)、[§6.4 F1 による長さ補正と、フリンジ領域の分割](#64-f1-による長さ補正とフリンジ領域の分割) ((5.43)、(5.43a)〜(5.43c)。放射用の 2 分割 `rbh`・`rbl`)。
6. [§7 `tbendcore` の読み方](#7-tbendcore-の読み方)・[§7.1 `t2t3`・`x3`・`da` の導出](#71-t2t3x3da-の導出-紙の上で-1-段ずつ): Fortran の長い式を、加法定理と和積公式で (5.22)〜(5.24) に結ぶ。
7. [§8 K1 ≠ 0 (複合機能の偏向磁石)](#8-k1--0-複合機能の偏向磁石): (5.50) の分割、K1 の級数、`tbendicorr` の各行。
8. [§9 行列 (CALC・Emittance[]) と追跡の関係](#9-行列-calcemittance-と追跡の関係): 同じ段を COD のまわりで微分する。差の出どころ (B-1)。
9. [§10 Fortran の対応表](#10-fortran-の対応表)・[§11 pysad の対応表](#11-pysad-の対応表)・[§12 確認スクリプトと出力](#12-確認スクリプトと出力)。
10. [§13 未確認・限界・使った条件](#13-未確認限界使った条件)・[§14 この章で使う Fortran の変数名](#14-この章で使う-fortran-の変数名)・[付録 A](#付録-a-ヘルプの記述k64-との違い-偏向磁石)・[付録 B](#付録-b-k64-の誤り-確定--b-1b-4疑い--b-2b-3)。

主な式の連鎖: (5.1) → (5.2) → (5.3) → (5.4) → (5.5) → (5.6)〜(5.8)、(5.10) → (5.11a)〜(5.11c) → (5.11)、(5.12) → (5.20) → (5.21) → (5.22) → (5.23) → (5.24) → (5.25)、(5.30) → (5.31) → (5.32) → (5.33) → (5.34) → (5.35) → (5.36)、(5.40) → (5.41) → (5.42)、(5.43a) → (5.43b) → (5.43c) → (5.43)。

**読み方**: 先に [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) の §9.4、[01_variables_and_conventions.md](01_variables_and_conventions.md) の §1〜§2・§8 (記号の規約)、[04_quad.md](04_quad.md) を読む。この章の $p_x$ は断らない限り正準運動量で、追跡の内部変数は $q = p_x/p$ (§2)。$\exp(:f:)$ の符号は [C_poisson_lie.md](C_poisson_lie.md) の約束に統一してあり、ヘルプの $H$ とは $H = -f$ の関係 (§1)。§7・§7.1 は Fortran の式を紙の上で追うための節で、飛ばしてよい。コードを読むだけなら §10・§11・§14 から入る。

**未確認・限界の所在**: [§13 未確認・限界・使った条件](#13-未確認限界使った条件) (Maxwell フリンジの $y^4$ の係数、F1 の式、$L'$ の補正、K1 ≠ 0、CALC と追跡の違い、走らせていない検証)。K64 の誤りの確定・疑いは付録 B、ヘルプとの違いは付録 A。

## 1. 結論

- **K1 = 0 の BEND の本体は、$H_2$ の近似ではなく、厳密解**を使う。一様な垂直磁場の中の粒子は、水平面内では円を描き、垂直方向には一定の速さで進む (らせん)。K64 の `tbendcore` は、この円と、端面を表す 2 本の直線の幾何を、閉じた式 (面内の角度は `asinz`) で書いたもの。円と面の幾何を実験室系で独立に解いた基準 (端面の Maxwell フリンジの式 (5.36) と F1 の閉じた式 (5.42) は K64・ヘルプと同じものを入れた。基準と本体が独立なのは、円と面の交点の幾何の部分) と、端面角 (E1・E2・AE1・AE2)・ANGLE と K0 が違う場合 ($\rho_b \neq \rho_0$)・$\delta = 0, +0.02, -0.03$ のすべてで $1.1 \times 10^{-15}$ 以内で一致した ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (b))。
- 追跡の 1 要素は、次の順の合成: [取り付け誤差 (`tbshift`・`tbrot`)] → **入口の面までの直進** → **線形フリンジ (F1)** → **非線形の Maxwell フリンジ** → **円弧の本体** → Maxwell フリンジ (出口) → 線形フリンジ (出口) → **出口の面から基準面への直進** → [逆の取り付け誤差]。端面の 3 段 (直進・F1・Maxwell) と円弧が、`tbendcore` の 1 つのサブルーチンにまとまっている (§7)。
- **端面角の効果**は、面までの直進 (5.20) と円弧の幾何から出る。線形光学では、水平に $+\tan\psi/\rho$、垂直に $-\tan\psi/\rho$ の薄い蹴り (面の傾き $\psi = \mathrm{ANGLE}\cdot\mathrm{E} + \mathrm{AE}$)。垂直の収束は、Maxwell フリンジの $\Delta p_y = -(p_x/p_v)\,y/\rho$ ($p_v = \sqrt{p^2 - p_x^2}$。(5.33)) に、面の座標系の $p_x \simeq \sin\psi$ を入れたもの。設計軌道の 2×2 行列は、端面の薄い蹴りと扇形の本体の積と $5 \times 10^{-13}$ で一致する ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (e))。
- **K1 ≠ 0 (複合機能)** は、本体を $H_2$ (曲率つきの四極。厳密に解ける) と、残りの非線形項に分け、スライスごとに交互に適用する (§8)。端面は薄い蹴りで近似する。$H_2$ の解の `xsxkx` の符号が $K1 > 0$ の枝で逆になっている (K64 の誤りの既知の例、付録 B)。
- **非線形の Maxwell フリンジ**は、端の磁場 ($\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ を $y$ の 3 次まで満たすもの) から導く生成関数 $f = -\dfrac{y^2}{2\rho}\dfrac{p_x}{\sqrt{p^2 - p_x^2}}$ ($y^4$ の補正を除く) の $\exp(:f:)$ で、閉じた式 (5.35) になる。この章の $\exp(:f:)$ は、[C_poisson_lie.md](C_poisson_lie.md) の約束 (括弧 $[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg)$。長さ $s$ の転送は $f = -sH$ で、$\Delta x = -\partial f/\partial p_x$、$\Delta p_x = +\partial f/\partial x$) に統一してある。SAD のヘルプの $\exp(:H:)$ は括弧の符号が逆なので、ヘルプの $H$ はこの章の $f$ の符号を反対にしたもの ($H = -f$。C §8.2)。磁場の端の幅 $d$ を持つ Maxwell の場の中の運動を積分すると、$d \to 0$ の外挿で、$\delta = 0$ と $\delta = 0.3$ のどちらも閉じた式と一致する ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。
- **線形フリンジ (F1)** は、ヘルプの定義 $\mathrm{F1} = 6\int (g - g^2)\,ds$ ($g = B_y(s)/B_0$) から導ける。垂直の収束は $F_1/(6\rho^2)$、$y^3$ の項は $\tfrac23\cdot\dfrac{1}{F_1\rho^2}$、水平の $\delta$ 依存のずれは $F_1^2\delta/(24\rho p)$。この 3 つが同じ F1 だけで書けるのは、端が直線の立ち上がり (と $\tanh$ の一部) のとき。一般の端の形では近似 (§6.3)。水平のずれは面に垂直な入射の式で、面が傾いて斜めに入射する場合は $1/\cos^3\psi$ 倍になる (§6.4 の (5.43a)。K64・ヘルプにこの因子は無い)。
- **転送行列 (CALC・Emittance[])** は、追跡と同じ段を COD のまわりで微分して掛けたもの (§9)。設計軌道の近くでは追跡のヤコビアンと $10^{-12}$ 台で一致する。軌道が大きく外れたとき ($y = 10$ mm、$\delta = 0.05$) は、COD の更新の差が $2 \times 10^{-8}$ (F1 = 0)。この差は、CALC の `tbedgemaxwell` の COD の更新が $\Delta x$ と $\Delta z$ で因子 $p = 1 + \delta$ を落としていること (K64 の誤り。gs で確認、v2 は直した。付録 B-1) で説明がつく: その更新だけを $p^2$ 型に差し替えると、$x$・$p_x$ の差は $10^{-16}$ 台、$z$ の差は $2 \times 10^{-8}$ から $1.4 \times 10^{-9}$ に下がる ([bend_07_b1_estimate.out.txt](outputs/bend_07_b1_estimate.out.txt) (b)。K64 の更新を差し替える試験。gs の COD との比較は付録 B-1)。F1 ≠ 0 の 4×4 の差は、k64 で $4 \times 10^{-5}$、v2 (B324 の修正を含む) で $2 \times 10^{-7}$ 以下に下がる (§9。下がった分の内訳は未分離。残りの原因は未特定)。
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

追跡の内部変数 $q = p_x/p$ ($p$ は一定) で書く。$q_x' = p_x'/p = -\dfrac{K + hh_b}{p}x + \left(h - \dfrac{h_b}{p}\right)$、$q_y' = p_y'/p = \dfrac{K}{p}y$ で、K64 の名前を使うと ($\rho_e = \rho_b p$)

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

(導出: $h_b = h = 1/\rho$ で、(5.6) の平衡点は $x_e = (hp - h)/(h\cdot h) = h\delta/h^2 = \rho\delta$、波数は $\omega_x^2 = h^2/p$ なので 0 次で $\omega_x = 1/\rho$。初期値 $x_0 = q_0 = 0$ を (5.7) に入れると $x(\ell) = x_e(1 - \cos\varphi_x)$、$q(\ell) = x_e\,\omega_x\sin\varphi_x$。$\delta$ の 1 次を取ると $x = \rho\delta(1 - \cos(L/\rho))$、$q = \delta\sin(L/\rho)$ ($\omega_x$ の $\delta$ 依存は 2 次の量にしか効かない)。$p_x = pq = q + O(\delta^2)$ なので $R_{26}$ も同じ。sympy: [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (c)、[bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (l)。)

$R_{56}$ は (5.10) から。$z' = -dv - hx - \tfrac12q^2$ で、$dv = -\delta/\gamma_0^2$ ($v/v_0 = 1 + \delta/\gamma_0^2$ の 1 次)、$x = \rho\delta(1 - \cos(\ell/\rho))$、$q^2 = O(\delta^2)$ なので、$\int_0^Lz'\,d\ell = \delta L/\gamma_0^2 - \delta\int_0^L(1 - \cos(\ell/\rho))\,d\ell = \delta\left[-(L - \rho\sin(L/\rho)) + L/\gamma_0^2\right]$。 ANGLE = 0.1 rad、$L = 1$ m で $R_{16} = 0.049958347$ m、$R_{26} = 0.099833417$ が、CALC の行列と 9 桁一致。

### 4.2 $z$ の更新

$z$ の方程式は、(5.3) の $H_2$ の $p$ 微分 ($z' = \partial H/\partial p$。$z$ の共役は $p$) で出る。各項を 1 つずつ微分する: $-p \to -1$、$(h_b - hp)\,x \to -hx$、$Q/(2p) \to -Q/(2p^2) = -\tfrac12(q_x^2 + q_y^2)$、$E/(p_0v_0) \to v/v_0 = 1 - dv$ (01 章)。$-1$ と $1 - dv$ をまとめて

$$
z' = -dv - h\,x - \tfrac12\left(q_x^2 + q_y^2\right)
\tag{5.10}
$$

この右辺を長さ $\ell$ にわたって積分する。必要な積分は $\int x\,d\ell$、$\int q_x^2\,d\ell$、$\int q_y^2\,d\ell$ の 3 つ。以下、$u \equiv x - x_e$、$\omega \equiv \omega_x$ とおく ((5.7) から $u' = q_x$、$q_x' = -\omega^2u$)。

**手順 1: $\int x\,d\ell$**。$x = x_e + u_0\cos\omega\ell + (q_0/\omega)\sin\omega\ell$ ((5.7)) の各項を積分する ($\int_0^\ell\cos\omega s\,ds = \sin\varphi_x/\omega$、$\int_0^\ell\sin\omega s\,ds = (1 - \cos\varphi_x)/\omega^2$):

$$
\int_0^\ell x\,ds = x_e\ell + u_0\frac{\sin\varphi_x}{\omega} + q_0\frac{1 - \cos\varphi_x}{\omega^2}
\tag{5.11a}
$$

K64 の記号 ($\texttt{sxkx} = \sin\varphi_x/\omega$、$\texttt{dcxkx} = (1 - \cos\varphi_x)/\omega^2$、$\texttt{xsxkx} = (\ell - \texttt{sxkx})/\texttt{akxsq}$、$\texttt{akxsq} = -\omega^2$、$\texttt{drhop} = \omega^2x_e$) で $T \equiv \texttt{drhop}\cdot\texttt{xsxkx} - \texttt{sxkx}\,x_0 - \texttt{dcxkx}\,q_0$ を書き下すと、$T = -x_e(\ell - \sin\varphi_x/\omega) - x_0\sin\varphi_x/\omega - q_0(1 - \cos\varphi_x)/\omega^2 = -x_e\ell - u_0\sin\varphi_x/\omega - q_0(1 - \cos\varphi_x)/\omega^2$ ($x_0 = x_e + u_0$ を使った)。これは (5.11a) の符号を反対にしたものなので

$$
T = -\int_0^\ell x\,ds
\tag{5.11b}
$$

**手順 2: $\int q_x^2\,d\ell$**。2 つの量を作る。$C_x \equiv q_x^2 + \omega^2u^2$ は保存量: $dC_x/d\ell = 2q_x(-\omega^2u) + 2\omega^2u\,q_x = 0$。また $d(uq_x)/d\ell = q_x\cdot q_x + u\cdot(-\omega^2u) = q_x^2 - \omega^2u^2$。2 式を足すと $C_x + d(uq_x)/d\ell = 2q_x^2$ なので

$$
\int_0^\ell q_x^2\,ds = \frac12\Bigl[C_x\,\ell + \Delta(u\,q_x)\Bigr],\qquad \Delta(uq_x) = (uq_x)|_{\rm out} - (uq_x)|_{\rm in}
\tag{5.11c}
$$

垂直は $q_y' = \texttt{akysq}\,y$ ($\texttt{akysq} = K/p$)。$C_y \equiv q_y^2 - \texttt{akysq}\,y^2$ ($dC_y/d\ell = 2q_y\,\texttt{akysq}\,y - 2\,\texttt{akysq}\,y\,q_y = 0$)、$d(yq_y)/d\ell = q_y^2 + \texttt{akysq}\,y^2$ なので、$C_y + d(yq_y)/d\ell = 2q_y^2$ で、$\int q_y^2\,d\ell = \tfrac12[C_y\ell + \Delta(yq_y)]$。

**手順 3: 組み立て**。(5.10) を積分し、手順 1・2 を使う ($-h\int x = hT$):

$$
\Delta z = -dv\,\ell + hT - \tfrac14\bigl[C_x\ell + \Delta(uq_x)\bigr] - \tfrac14\bigl[C_y\ell + \Delta(yq_y)\bigr]
$$

**手順 4: K64 の `hi` と $\Delta(xq_x)$ に書き直す**。K64 の `hi` ($= \tfrac12(q_x^2 + q_y^2 - \texttt{akxsq}\,x^2 - \texttt{akysq}\,y^2) - \texttt{drhop}\,x$。`tbend.f:149`) は、$H_2$ の値から $-p$ を除いて $p$ で割ったもの、つまり $\mathrm{hi} = (H_2 + p)/p$ ($H_2$ は $p_x$ で書かれていて $p_x^2/(2p) = p\,q_x^2/2$ なので、割る $p$ が要る。$\mathrm{hi}$ 自体は $q$ の単位の量)。$\texttt{akxsq} = -\omega^2$、$\texttt{drhop} = \omega^2x_e$ から $\tfrac12\omega^2x^2 - \omega^2x_ex = \tfrac12\omega^2(x - x_e)^2 - \tfrac12\omega^2x_e^2$ なので

$$
\mathrm{hi} = \tfrac12C_x + \tfrac12C_y - \tfrac12\omega^2x_e^2\quad\Rightarrow\quad \tfrac14(C_x + C_y)\ell = \tfrac12\,\mathrm{hi}\,\ell + \tfrac14\omega^2x_e^2\ell
$$

(sympy で $\mathrm{hi} = (H_2 + p)/p$ と上の関係を確認: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (a)。) また $\Delta(uq_x) = \Delta(xq_x) - x_e\Delta q_x$ ($x_e$ は定数)。手順 3 に入れると

$$
\Delta z = -dv\,\ell + hT - \tfrac14\Delta(xq_x) - \tfrac14\Delta(yq_y) - \tfrac12\,\mathrm{hi}\,\ell - \tfrac14x_e\bigl[\omega^2x_e\ell - \Delta q_x\bigr]
$$

**手順 5: $\Delta q_x$ を消す**。$\Delta q_x = \int q_x'\,d\ell = \int(\texttt{drhop} - \omega^2x)\,d\ell = \omega^2x_e\ell - \omega^2\int x\,d\ell = \omega^2x_e\ell + \omega^2T$ ((5.11b))。よって $\omega^2x_e\ell - \Delta q_x = -\omega^2T$ で、最後の項は $+\tfrac14x_e\omega^2T = +\tfrac14\texttt{drhop}\,T = \tfrac14\left(\tfrac1{\rho_0} - \tfrac1{\rho_e}\right)T$ ($\texttt{drhop} = h - h_b/p$)。$hT = T/\rho_0$ と足すと $\left[\tfrac1{\rho_0} + \tfrac14\left(\tfrac1{\rho_0} - \tfrac1{\rho_e}\right)\right]T = \tfrac14\left(\tfrac5{\rho_0} - \tfrac1{\rho_e}\right)T$ で、

$$
\Delta z = -\frac14\Bigl[\Delta(x q_x) + \Delta(y q_y)\Bigr] + \frac14\left(\frac{5}{\rho_0} - \frac{1}{\rho_e}\right)T - \frac{\mathrm{hi}\cdot\ell}{2} - dv\,\ell
\tag{5.11}
$$

$\Delta(xq_x) = x_fq_{xf} - x_iq_{xi} = \texttt{dxf}\cdot q_0 + x_0\,\texttt{dpxf} + \texttt{dxf}\cdot\texttt{dpxf}$ ($x_f = x_0 + \texttt{dxf}$、$q_{xf} = q_0 + \texttt{dpxf}$ を展開した。これが K64 の `dxf*pxi+xi*dpxf+dxf*dpxf`)。$\Delta(yq_y)$ も同じ形。**係数 $5/\rho_0 - 1/\rho_e$ の出どころ**は手順 5 の 2 項 ($T/\rho_0$ は $-h\int x$ から、$\tfrac14\texttt{drhop}\,T$ は $-\tfrac12\int q_x^2$ の中の $x_e\Delta q_x$ から)。(5.11a)〜(5.11c) と (5.11) は、数値積分 (mpmath、30 桁) と $10^{-43}$ で一致する ($K = 0$、$K > 0$、$K < 0$。$\rho_0 = 10$ m、$\rho_b = 9.7$ m、$p = 1.01$、$L = 1.3$ m。[bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (b))。

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

(導出: $\vec b = (0, h_b, 0)$ ($(X, Y, Z)$ の順) で $\vec P\times\vec b = (-P_Zh_b,\ 0,\ P_Xh_b)$ なので $dP_X/d\tau = -h_bP_Z$、$dP_Z/d\tau = +h_bP_X$ ($P_Y$ は一定)。これは角速度 $h_b$ の回転で、$\Theta = h_b\tau$ の関数として $P_X = P_{X0}\cos\Theta - P_{Z0}\sin\Theta$、$P_Z = P_{Z0}\cos\Theta + P_{X0}\sin\Theta$ (微分して $dP_X/d\Theta = -P_Z$、$dP_Z/d\Theta = P_X$ を確かめる)。位置は $d\tau = \rho_b\,d\Theta$ なので $X = X_0 + \rho_b\int_0^\Theta P_X\,d\Theta' = X_0 + \rho_b[P_{X0}\sin\Theta + P_{Z0}(\cos\Theta - 1)]$、$Z = Z_0 + \rho_b[P_{Z0}\sin\Theta + P_{X0}(1 - \cos\Theta)]$、$Y = Y_0 + P_Y\rho_b\Theta$。) 円の中心は $(X_c, Z_c) = (X_0 - \rho_bP_{Z0},\ Z_0 + \rho_bP_{X0})$ で、$(X - X_c,\ Z - Z_c) = \rho_b\,(P_Z,\ -P_X)$。3 次元の経路長は $\ell = p\tau = p\rho_b\Theta$。$z$ は $z_{\rm out} = z_{\rm in} - \ell + (v/v_0)L'$ ($= -\ell + (1 - dv)L'$)。これは $z' = -(1 + hx)p/P_z + v/v_0$ (01 章の式 (1.8)) を積分したものと同じ。

**曲線座標の $H$ の運動方程式との一致**: (5.2) の正準方程式 $x' = (1 + hx)p_x/P_z$、$p_x' = -h_b(1 + hx) + hP_z$、$y' = (1 + hx)p_y/P_z$、$p_y' = 0$、$z' = 1 - dv - (1 + hx)p/P_z$ ($P_z = \sqrt{p^2 - p_x^2 - p_y^2}$) を DOP853 で積分した結果と、(5.12) の幾何の差は、$\delta = 0, 0.02, -0.03$ で $1.1 \times 10^{-15}$ ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (f))。

### 5.2 入口の面までの直進

実験室系で、入口の基準面は $Z = 0$。入口の面は、基準軌道の入口の点 (原点) を通り、方向 $\hat e_{x1} = (\cos\psi_1, \sin\psi_1)$、法線 $\hat e_{z1} = (-\sin\psi_1, \cos\psi_1)$ の直線 (正の $\psi_1$ は、面が $x > 0$ の側で下流へ進む向き)。基準面の点 $(x_1, 0)$ から運動量 $(P_X, P_Z)$ の向きに直進して面に当たる点を $x_2\hat e_{x1}$ とする。直進のパラメータ $t$ ($d\vec r/dt = \vec P$) で、面の条件 $(x_1 + P_Xt)\sin\psi_1 = P_Zt\cos\psi_1$ から

$$
t = \frac{x_1\sin\psi_1}{P_Z\cos\psi_1 - P_X\sin\psi_1}
\ \Rightarrow\ x_2 = \frac{x_1}{\cos\psi_1 - (P_X/P_Z)\sin\psi_1}
\tag{5.20}
$$

(5.20) の導出の中身: 面は原点を通り方向 $\hat e_{x1}$ の直線なので、面の上の点は法線 $\hat e_{z1}$ との内積が 0。点 $\vec r = (x_1 + P_Xt,\ P_Zt)$ について $\hat e_{z1}\cdot\vec r = -(x_1 + P_Xt)\sin\psi_1 + P_Zt\cos\psi_1 = 0$ (上の面の条件と同じ)。これを $t$ について解くと $t\,(P_Z\cos\psi_1 - P_X\sin\psi_1) = x_1\sin\psi_1$ で、(5.20) の $t$ になる。面上の座標は $x_2 = \hat e_{x1}\cdot\vec r = (x_1 + P_Xt)\cos\psi_1 + P_Zt\sin\psi_1$。面の条件から $P_Zt\cos\psi_1 = (x_1 + P_Xt)\sin\psi_1$ なので、$(x_1 + P_Xt) = P_Zt\cos\psi_1/\sin\psi_1$ ($\psi_1 \neq 0$) を代入して $x_2 = P_Zt\,(\cos^2\psi_1 + \sin^2\psi_1)/\sin\psi_1 = P_Zt/\sin\psi_1$、すなわち $P_Zt = x_2\sin\psi_1$ (これは $\psi_1 = 0$ でも $t = 0$ で成り立つ)。$t$ を代入すると $x_2 = P_Zt/\sin\psi_1 = x_1P_Z/(P_Z\cos\psi_1 - P_X\sin\psi_1) = x_1/(\cos\psi_1 - (P_X/P_Z)\sin\psi_1)$ で、(5.20) の第 2 式。(sympy で面の条件の残差 0、$x_2$ の 2 通りの表現の差 0: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (c)。)

面の座標系 ($\hat e_{x1}$ が $x$、$\hat e_{z1}$ が $z$) での運動量成分と他の変数。$p_{x2}$ は $\vec P\cdot\hat e_{x1} = P_X\cos\psi_1 + P_Z\sin\psi_1$。$y$ は、直進の間 $P_y$ が一定なので $\Delta y = P_yt$ で、$t = x_2\sin\psi_1/P_Z$ を使うと $P_yt = (p_y/P_Z)x_2\sin\psi_1$。$z$ は、直進の間に 3 次元の経路長 $|\vec P|\,t = pt$ だけ進むので ($z' = -(\text{経路長}) + \ldots$ の経路長の部分。$v/v_0$ の項は、本体の $L'$ にまとめて (5.24) で入れる。§7 の表の段 8)、$\Delta z = -pt = -(p/P_Z)\,x_2\sin\psi_1$:

$$
\begin{aligned}
p_{x2} &= P_X\cos\psi_1 + P_Z\sin\psi_1,\qquad y_2 = y_1 + \frac{p_y}{P_Z}\,x_2\sin\psi_1,\\
z_2 &= z_1 - \frac{p}{P_Z}\,x_2\sin\psi_1,\qquad p_{y2} = p_y
\end{aligned}
\tag{5.21}
$$

(これが `tbendcore` の入口の `x1 = x/(cosp1 - pxi*sp1)` ($\texttt{sp1} = \sin\psi_1/p_z$)、`dz1 = x1*sp1`、`y1 = y + pyi*dz1`、`px1 = pxi*cosp1 + pzi*sinp1`。ヘルプの "drift to the entrance face" と同じ式 (原文で確認)。) 出口の面から出口の基準面への直進は、次の小節で導く (ヘルプの "drift from the exit face" に当たる)。

**出口の面から基準面への直進**。出口の面は基準軌道の出口 $E$ を通り、面の基底は $\hat e_{f2} = (\cos(\phi - \psi_2),\ \sin(\phi - \psi_2))$、法線 $\hat N_2 = (-\sin(\phi - \psi_2),\ \cos(\phi - \psi_2))$ (§5.3)。出口の基準面も $E$ を通り、基底は $\hat e_{\rm ref} = (\cos\phi,\ \sin\phi)$、法線 $\hat N_{\rm ref} = (-\sin\phi,\ \cos\phi)$ ($\psi_2 = 0$ に当たる面)。面の座標系 $(x_2, p_{x2}, p_{z2})$ から基準面の座標系 $(x_{\rm out}, p_{xr}, p_{zr})$ への変換を、入口の (5.20)(5.21) と同じ手順で作る。

1. **基底の内積** (加法定理): $\hat e_{\rm ref}\cdot\hat e_{f2} = \cos\phi\cos(\phi - \psi_2) + \sin\phi\sin(\phi - \psi_2) = \cos\psi_2$。$\hat e_{\rm ref}\cdot\hat N_2 = -\cos\phi\sin(\phi - \psi_2) + \sin\phi\cos(\phi - \psi_2) = \sin\psi_2$。$\hat N_{\rm ref}\cdot\hat e_{f2} = -\sin\phi\cos(\phi - \psi_2) + \cos\phi\sin(\phi - \psi_2) = -\sin\psi_2$。$\hat N_{\rm ref}\cdot\hat N_2 = \sin\phi\sin(\phi - \psi_2) + \cos\phi\cos(\phi - \psi_2) = \cos\psi_2$。
2. **運動量の成分**: $\vec P = p_{x2}\hat e_{f2} + p_{z2}\hat N_2$ の $\hat e_{\rm ref}$・$\hat N_{\rm ref}$ 成分を取ると、1 の内積から $p_{xr} = p_{x2}\cos\psi_2 + p_{z2}\sin\psi_2$、$\quad p_{zr} = -p_{x2}\sin\psi_2 + p_{z2}\cos\psi_2$。
3. **直進の時間**: 面の上の点 $\vec r = E + x_2\hat e_{f2}$ から、$d\vec r/dt = \vec P$ で進んで基準面 $(\vec r' - E)\cdot\hat N_{\rm ref} = 0$ に当たる条件は $x_2(\hat e_{f2}\cdot\hat N_{\rm ref}) + t\,p_{zr} = 0$。1 で $\hat e_{f2}\cdot\hat N_{\rm ref} = -\sin\psi_2$ なので $-x_2\sin\psi_2 + t\,p_{zr} = 0$、すなわち $t = x_2\sin\psi_2/p_{zr}$。
4. **基準面の座標**: $x_{\rm out} = (\vec r' - E)\cdot\hat e_{\rm ref} = x_2\cos\psi_2 + t\,p_{xr}$ ($\hat e_{f2}\cdot\hat e_{\rm ref} = \cos\psi_2$)。3 の $t$ を入れる。
5. **$y$ と $z$**: 直進の間 $p_y$ は一定なので $\Delta y = p_yt$、経路長は $pt$ で $z$ は $pt$ だけ遅れる (§5.2 の入口と同じ)。まとめて

$$
\begin{aligned}
&p_{xr} = p_{x2}\cos\psi_2 + p_{z2}\sin\psi_2,\qquad p_{zr} = p_{z2}\cos\psi_2 - p_{x2}\sin\psi_2,\qquad x_{\rm out} = x_2\left(\cos\psi_2 + \frac{p_{xr}}{p_{zr}}\sin\psi_2\right),\\
&y_{\rm out} = y_2 + \frac{p_y}{p_{zr}}\,x_2\sin\psi_2,\qquad z_{\rm out} = z_2 - \frac{p}{p_{zr}}\,x_2\sin\psi_2
\end{aligned}
\tag{5.21a}
$$

($p_{zr} > 0$ の向きに進むとき。$\psi_2 = 0$ なら恒等。`tbendcore` の最後の `dz4 = x4*sinp2/(cosp2*(1+dpz4) + sinp2*px3)` は、$q$ の単位で $p_{zr}$ を $\cos\psi_2\,p_z + \sin\psi_2\,(-q_{x2})$ と書いた分母 (`px3` $= -q_{x2}$、§5.3) で、3 の $t = x_2\sin\psi_2/p_{zr}$ に当たる。sympy で 1 の 4 つの内積、2 の成分、4 の式を確認し、$t$ の入った点が基準面に載ること (残差 $10^{-142}$) を数値で確認した: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (a)。) 入口の (5.20) は「基準面から面へ」の変換で $x_1$ に $1/(\cos\psi_1 - (p_x/p_z)\sin\psi_1)$ がかかるのに対し、これは「面から基準面へ」の変換で、係数は逆数の形 $\cos\psi_2 + (p_{xr}/p_{zr})\sin\psi_2$ になる。

**シンプレクティック性**: (5.20)(5.21) に $\delta$ (変えない) を加えた 6 変数 $(x, p_x, y, p_y, z, \delta)$ の写像のヤコビアン $M$ について、$M^\top SM - S$ を 40 桁で評価すると最大 $6 \times 10^{-42}$ ($\psi = 0.3$ rad、$(x, p_x, y, p_y, \delta) = (3\ \mathrm{mm}, 2\ \mathrm{mrad}, 2\ \mathrm{mm}, -1\ \mathrm{mrad}, 0.02)$ の 1 点。[bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (a))。面は固定された平面で、その上の位置と運動量の取り替え (自由粒子の直進) なので、正準変換になる。

### 5.3 面から面への円弧

入口の面の座標系で、位置 $(x_1, 0)$、運動量 $(p_{x1}, p_{z1})$ ($p_{z1} = \sqrt{p_h^2 - p_{x1}^2}$、$p_h^2 = p^2 - p_y^2$ は水平面内の運動量の大きさの 2 乗。ここでは $x_1$ は (5.20) の $x_2$、$p_{x1}$ は (5.21) の $p_{x2}$ に当たる) から始める。(5.12) の円の中心は、面の座標系 ($\hat e_{x1}$, $\hat e_{z1}$) で $(x_1 - \rho_bp_{z1},\ \rho_bp_{x1})$。出口の面の法線 $\hat N_2$ と方向 $\hat e_{f2}$ を実験室系で書く: 出口の基準面は基準軌道の出口 $E = (-\rho_0(1 - \cos\phi),\ \rho_0\sin\phi)$ ($\phi = \mathrm{ANGLE}$) を通り法線 $(-\sin\phi, \cos\phi)$、出口の面はそれを $\psi_2$ だけ傾けたもの:

$$
\hat e_{f2} = (\cos(\phi - \psi_2),\ \sin(\phi - \psi_2)),\qquad \hat N_2 = (-\sin(\phi - \psi_2),\ \cos(\phi - \psi_2))
$$

入口の面の座標系の基底との内積は $\hat e_{x1}\cdot\hat e_{f2} = \cos\omega$、$\hat e_{z1}\cdot\hat e_{f2} = \sin\omega$、$\hat e_{x1}\cdot\hat N_2 = -\sin\omega$、$\hat e_{z1}\cdot\hat N_2 = \cos\omega$ ($\omega = \phi - \psi_1 - \psi_2$)。基準軌道の出口の点 $E$ との内積は $E\cdot\hat N_2 = \rho_0[\sin(\omega + \psi_1) + \sin\psi_2]$、$E\cdot\hat e_{f2} = \rho_0[\cos\psi_2 - \cos(\omega + \psi_1)]$。

(内積の導出。$\hat e_{x1}\cdot\hat e_{f2} = \cos\psi_1\cos(\phi - \psi_2) + \sin\psi_1\sin(\phi - \psi_2) = \cos(\phi - \psi_2 - \psi_1) = \cos\omega$、$\hat e_{z1}\cdot\hat e_{f2} = -\sin\psi_1\cos(\phi - \psi_2) + \cos\psi_1\sin(\phi - \psi_2) = \sin(\phi - \psi_2 - \psi_1) = \sin\omega$。$\hat N_2$ との内積も同様。$E\cdot\hat N_2 = \rho_0[(1 - \cos\phi)\sin(\phi - \psi_2) + \sin\phi\cos(\phi - \psi_2)] = \rho_0[\sin(\phi - \psi_2) + (\sin\phi\cos(\phi - \psi_2) - \cos\phi\sin(\phi - \psi_2))] = \rho_0[\sin(\phi - \psi_2) + \sin\psi_2]$、ここで $\phi - \psi_2 = \omega + \psi_1$。$E\cdot\hat e_{f2} = \rho_0[-(1 - \cos\phi)\cos(\phi - \psi_2) + \sin\phi\sin(\phi - \psi_2)] = \rho_0[-\cos(\phi - \psi_2) + \cos(\phi - (\phi - \psi_2))] = \rho_0[\cos\psi_2 - \cos(\omega + \psi_1)]$。sympy で全部確認: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (c)。)

**出口の面に当たる条件は、出口の面の座標系の運動量 $p_{x2}$ について線形**: 出口の位置 $\vec r$ は $\vec r - \vec c = \rho_b(P_Z, -P_X)$ で、$(P_Z, -P_X)\cdot\hat N_2 = -\vec P\cdot\hat e_{f2} = -p_{x2}$ (出口の面の座標系の $x$ 運動量)。条件 $(\vec r - E)\cdot\hat N_2 = 0$ は $(\vec c - E)\cdot\hat N_2 - \rho_bp_{x2} = 0$ になり、$\vec c\cdot\hat N_2 = -(x_1 - \rho_bp_{z1})\sin\omega + \rho_bp_{x1}\cos\omega$ を入れて

$$
p_{x2} = -\frac{x_1}{\rho_b}\sin\omega + p_{z1}\sin\omega + p_{x1}\cos\omega - \frac{\rho_0}{\rho_b}\bigl[\sin\psi_2 + \sin(\omega + \psi_1)\bigr]
\tag{5.22}
$$

位置は $x_2 = (\vec r - E)\cdot\hat e_{f2} = (\vec c - E)\cdot\hat e_{f2} + \rho_bp_{z2}$ ($p_{z2} = \sqrt{p_h^2 - p_{x2}^2}$)。ここで $\vec r - \vec c = \rho_b(P_Z, -P_X)$ の $\hat e_{f2}$ 成分は $\rho_b(P_Z\cos(\phi - \psi_2) - P_X\sin(\phi - \psi_2)) = \rho_b\,\vec P\cdot\hat N_2 = \rho_bp_{z2}$ (出口の面の座標系の $z$ 運動量)。$\vec c\cdot\hat e_{f2}$ は、$\vec c$ の入口の面の座標系の成分 $(x_1 - \rho_bp_{z1},\ \rho_bp_{x1})$ と上の内積から $(x_1 - \rho_bp_{z1})\cos\omega + \rho_bp_{x1}\sin\omega$。$E\cdot\hat e_{f2}$ の式を引くと:

$$
x_2 = x_1\cos\omega + \rho_b\bigl(p_{z2} - p_{z1}\cos\omega + p_{x1}\sin\omega\bigr) + \rho_0\bigl[\cos(\omega + \psi_1) - \cos\psi_2\bigr]
\tag{5.23}
$$

運動量の回転角 $\Theta$ は、面内の向きの差から出る。入口の面の座標系で $p_{x1} = p_h\sin\theta_1$ ($\theta_1$ は面の法線 $\hat e_{z1}$ から $\hat e_{x1}$ の向きに測った運動量の角)、出口の面の座標系で $p_{x2} = p_h\sin\theta_2$ とおく。実験室系の $Z$ 軸から $+X$ へ測った運動量の角は、入口で $-\psi_1 + \theta_1$ (法線 $\hat e_{z1}$ が $Z$ 軸から $-\psi_1$)、出口で $-(\phi - \psi_2) + \theta_2$ (法線 $\hat N_2$ が $Z$ 軸から $-(\phi - \psi_2)$)。$-X$ の向きを正とする回転角は

$$
\Theta = \bigl[(\phi - \psi_2) - \theta_2\bigr] - \bigl[\psi_1 - \theta_1\bigr] = \omega + \theta_1 - \theta_2 = \arcsin\frac{p_{x1}}{p_h} - \arcsin\frac{p_{x2}}{p_h} + \omega
$$

(設計軌道は $\theta_1 = \psi_1$、$\theta_2 = -\psi_2$ なので $\Theta = \omega + \psi_1 + \psi_2 = \phi$ で、円弧の角が $\mathrm{ANGLE}$ に戻る。) 垂直と $z$ は (5.12) から

$$
y_2 = y_1 + p_y\rho_b\Theta,\qquad
z_2 = z_1 - p\rho_b\Theta + \frac{v}{v_0}L'\ \ (\text{直進の分は (5.21) で入れてある})
\tag{5.24}
$$

(5.22)(5.23) は、ヘルプの "body of bend" の $p_{x2}$、$x_2$ の式と**同じ** (原文で確認。ただし $p_{x2}$ の $x_1$ の項は、平文版が $-(x_1/\rho_b)\sin\omega$ なのに対し TeX 版は $-x_1/(\rho_b\sin\omega)$ で、$\omega \to 0$ で発散する誤植。平文版が (5.22) と同じ。[bend_05_help_checks.out.txt](outputs/bend_05_help_checks.out.txt) (e))。ヘルプの $s$ の式 ($s = \mathrm{ANGLE}\times\rho_b\times(\cdots)$) の $\mathrm{ANGLE}\times$ は、面内の経路長 $s = \rho_bp_h\Theta$ ($y_2 = y_1 + s\,p_y/p_h$ が (5.24) と一致する) の誤記と読む (原文は TeX 版・平文版とも $\mathrm{ANGLE}\times\rho_b\times(\ldots)$ で、$p_h$ も無い)。(5.22)〜(5.24) とヘルプの式を、(5.12) の幾何 (円と面の交点を数値で解く) と比べると $2.4 \times 10^{-15}$ 以内で一致する ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (d))。

**K64 の `tbendcore` の対応** (§7 に行ごとの表): `ph2` $= p_h/p$ ($= \sqrt{1 - q_y^2}$)、`dph2` $= \texttt{ph2} - 1$、`drho1` $= \texttt{drhob} + \rho_b\delta = \rho_e - \rho_0$、`drho` $= \texttt{drho1} + \rho_e\,\texttt{dph2} = \rho_ep_h/p - \rho_0$、`dpx2` $= \texttt{px1} - \texttt{ph2}\sin\psi_1$ (面の座標系の $q_x$ から、設計の向きの粒子の $\texttt{ph2}\sin\psi_1$ を引いた量)、`dpx3a`・`dpx3`・`px3` は (5.22) を $q$ の単位 (÷$p$) で書いたもの。**`px3` は (5.22) の $-p_{x2}/p$**: `tbendcore` は出口の段を $x$ 運動量の符号を反対にした座標で組んでいる (出力側 `px(i) = -cosp2*dpx3 + ...` で戻す。$\psi_2 = 0$ なら $\texttt{px} = -\texttt{px3}$)。Fortran の本体の式を転記して (5.22) の符号つきの値と比べた (1 点で相対 $10^{-16}$)。`phi0a` $= \phi_0 + \texttt{da} = \Theta$: `y3 = y1 + py2*rhoe*phi0a` は (5.24) の $y_2$ ($\texttt{py2}\cdot\rho_bp = p_y\rho_b$)、`z3 = z2 - phi0*drho1 - da*rhoe - dv*al` は $-p\rho_b\Theta + L' - dv\,L'$ ($-\phi_0(\rho_e - \rho_0) - \texttt{da}\,\rho_e = -\rho_e(\phi_0 + \texttt{da}) + L'$。$\phi_0\rho_0 = L'$) に等しい。

### 5.4 端面角の効果 (線形光学)

設計軌道 ($x = p_x = 0$、$p = 1$、$p_y = 0$) のまわりで 1 次まで展開する。入口の面の座標系では、設計軌道は $p_{x1} = \sin\psi_1$、$p_{z1} = \cos\psi_1$ で面に当たる ((5.21))。この点のまわりの 1 次の変分は、(5.20)(5.21) から $\delta x_1 = \delta x/\cos\psi_1$、$\delta p_{x1} = \cos\psi_1\,\delta p_x$、$\delta p_{z1} = -\sin\psi_1\,\delta p_x$。これを (5.22)(5.23) と出口の面から基準面への直進に代入して 1 次を取る。3 つの段に分けて行列を書く ($\rho_b = \rho_0 = \rho$、$\omega = \phi - \psi_1 - \psi_2$)。

1. **入口の直進** (基準面 $(x, p_x)$ → 入口の面の座標系 $(x_1, p_{x1})$): $\delta x_1 = \delta x/\cos\psi_1$、$\delta p_{x1} = \cos\psi_1\,\delta p_x$。行列 $J_{\rm in} = \mathrm{diag}(1/\cos\psi_1,\ \cos\psi_1)$。**対角になる理由**: (5.20) は $x_1 = x/(\cos\psi_1 - (p_x/p_z)\sin\psi_1)$ で、$x$ に比例する。$\partial x_1/\partial p_x = x\cdot\partial[\,\cdots]/\partial p_x$ は、設計軌道で $x = 0$ なので 0。$\partial x_1/\partial x = 1/(\cos\psi_1 - 0) = 1/\cos\psi_1$ ($p_x = 0$)。$p_{x1} = p_x\cos\psi_1 + p_z\sin\psi_1$ ((5.21)) は $x$ を含まないので $\partial p_{x1}/\partial x = 0$。$\partial p_{x1}/\partial p_x = \cos\psi_1 + (\partial p_z/\partial p_x)\sin\psi_1 = \cos\psi_1$ ($\partial p_z/\partial p_x = -p_x/p_z = 0$)。sympy で $x = 0$ の微分を確認: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (b)。
2. **円弧** ((5.22)(5.23) の 1 次、設計点 $x_1 = 0$、$p_{x1} = \sin\psi_1$、$p_{z1} = \cos\psi_1$。このとき $p_{x2} = -\sin\psi_2$、$p_{z2} = \cos\psi_2$): まず $\partial p_{z1}/\partial p_{x1} = -p_{x1}/p_{z1} = -\tan\psi_1$、$\partial p_{z2}/\partial p_{x2} = -p_{x2}/p_{z2} = +\tan\psi_2$。(5.22) から
   $\dfrac{\partial p_{x2}}{\partial x_1} = -\dfrac{\sin\omega}{\rho}$、$\dfrac{\partial p_{x2}}{\partial p_{x1}} = -\sin\omega\,\tan\psi_1 + \cos\omega = \dfrac{\cos(\omega + \psi_1)}{\cos\psi_1} = \dfrac{\cos(\phi - \psi_2)}{\cos\psi_1}$。
   (5.23) を微分する ($x_2 = x_1\cos\omega + \rho(p_{z2} - p_{z1}\cos\omega + p_{x1}\sin\omega) + \text{定数}$。$p_{z2}$ は $p_{x2}$ を通して、$p_{z1}$ は $p_{x1}$ を通して $x_1$・$p_{x1}$ に依る)。
   - $x_1$ で: $\dfrac{\partial x_2}{\partial x_1} = \cos\omega + \rho\tan\psi_2\,\dfrac{\partial p_{x2}}{\partial x_1} = \cos\omega - \tan\psi_2\sin\omega = \dfrac{\cos\omega\cos\psi_2 - \sin\omega\sin\psi_2}{\cos\psi_2} = \dfrac{\cos(\omega + \psi_2)}{\cos\psi_2} = \dfrac{\cos(\phi - \psi_1)}{\cos\psi_2}$ ($\omega + \psi_2 = \phi - \psi_1$)。
   - $p_{x1}$ で: $\dfrac{\partial x_2}{\partial p_{x1}} = \rho\left[\tan\psi_2\dfrac{\partial p_{x2}}{\partial p_{x1}} + \tan\psi_1\cos\omega + \sin\omega\right]$ ($\partial p_{z2}/\partial p_{x1} = \tan\psi_2\,\partial p_{x2}/\partial p_{x1}$、$-\cos\omega\,\partial p_{z1}/\partial p_{x1} = +\tan\psi_1\cos\omega$、$\partial(p_{x1}\sin\omega)/\partial p_{x1} = \sin\omega$)。上の $\partial p_{x2}/\partial p_{x1} = \cos(\omega + \psi_1)/\cos\psi_1$ を入れて、**通分** ($\cos\psi_1\cos\psi_2$ をかける) すると
     $\dfrac{\rho}{\cos\psi_1\cos\psi_2}\bigl[\sin\psi_2\cos(\omega + \psi_1) + \cos\psi_2(\sin\psi_1\cos\omega + \cos\psi_1\sin\omega)\bigr] = \dfrac{\rho}{\cos\psi_1\cos\psi_2}\bigl[\sin\psi_2\cos(\omega + \psi_1) + \cos\psi_2\sin(\omega + \psi_1)\bigr]$ ($\sin\psi_1\cos\omega + \cos\psi_1\sin\omega = \sin(\omega + \psi_1)$ を使った) $= \dfrac{\rho\,\sin(\omega + \psi_1 + \psi_2)}{\cos\psi_1\cos\psi_2} = \dfrac{\rho\sin\phi}{\cos\psi_1\cos\psi_2}$ (加法定理。$\omega + \psi_1 + \psi_2 = \phi$)。
   - $p_{x2}$ の $p_{x1}$ 微分 (上の再掲): $\partial p_{x2}/\partial p_{x1} = -\sin\omega\tan\psi_1 + \cos\omega = (\cos\omega\cos\psi_1 - \sin\omega\sin\psi_1)/\cos\psi_1 = \cos(\omega + \psi_1)/\cos\psi_1$。

   (sympy でここの各式の差 0: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (b)。) まとめて
   $J_{\rm arc} = \begin{pmatrix}\cos(\phi - \psi_1)/\cos\psi_2 & \rho\sin\phi/(\cos\psi_1\cos\psi_2)\\ -\sin\omega/\rho & \cos(\phi - \psi_2)/\cos\psi_1\end{pmatrix}$。
3. **出口の直進** (出口の面の座標系 → 基準面。(5.21a)): $x_{\rm out} = x_2(\cos\psi_2 + (p_{xr}/p_{zr})\sin\psi_2)$、$p_{xr} = p_{x2}\cos\psi_2 + p_{z2}\sin\psi_2$。$p_{x2} = \sin\theta_2$ ($p_{z2} = \cos\theta_2$、$p = 1$、$p_y = 0$) とおくと $p_{xr} = \sin(\theta_2 + \psi_2)$、$p_{zr} = \cos(\theta_2 + \psi_2)$ で、設計点は $\theta_2 = -\psi_2$ ($p_{xr} = 0$、$p_{zr} = 1$)。**対角になる理由**: 設計点で $x_2 = 0$ なので、$x_{\rm out} = x_2\cdot(\cdots)$ の $p_{x2}$ 微分は $x_2\,\partial(\cdots)/\partial p_{x2} = 0$。$\partial x_{\rm out}/\partial x_2 = \cos\psi_2 + (p_{xr}/p_{zr})\sin\psi_2 = \cos\psi_2$ ($p_{xr} = 0$)。$p_{xr}$ は $x_2$ を含まないので $\partial p_{xr}/\partial x_2 = 0$。$\partial p_{xr}/\partial p_{x2} = \cos\psi_2 + (\partial p_{z2}/\partial p_{x2})\sin\psi_2 = \cos\psi_2 + \tan\psi_2\sin\psi_2 = 1/\cos\psi_2$ ($\partial p_{z2}/\partial p_{x2} = -p_{x2}/p_{z2} = +\tan\psi_2$)。$J_{\rm out} = \mathrm{diag}(\cos\psi_2,\ 1/\cos\psi_2)$。

積 $J_{\rm out}J_{\rm arc}J_{\rm in}$ (対角行列を左右からかけるので、$(i,j)$ 成分に左の対角の $i$ 番目と右の対角の $j$ 番目がかかる) の 4 成分を陽に書く。
   - $M_{11} = \cos\psi_2\cdot\dfrac{\cos(\phi - \psi_1)}{\cos\psi_2}\cdot\dfrac1{\cos\psi_1} = \dfrac{\cos(\phi - \psi_1)}{\cos\psi_1}$、
   - $M_{12} = \cos\psi_2\cdot\dfrac{\rho\sin\phi}{\cos\psi_1\cos\psi_2}\cdot\cos\psi_1 = \rho\sin\phi$、
   - $M_{21} = \dfrac1{\cos\psi_2}\cdot\left(-\dfrac{\sin\omega}{\rho}\right)\cdot\dfrac1{\cos\psi_1} = -\dfrac{\sin\omega}{\rho\cos\psi_1\cos\psi_2}$、
   - $M_{22} = \dfrac1{\cos\psi_2}\cdot\dfrac{\cos(\phi - \psi_2)}{\cos\psi_1}\cdot\cos\psi_1 = \dfrac{\cos(\phi - \psi_2)}{\cos\psi_2}$。

これを、次の標準的な端の行列と扇形の本体の積 $E_x(\psi_2)S_xE_x(\psi_1)$ と比べる ($c \equiv \cos\phi$、$s \equiv \sin\phi$、$t_i \equiv \tan\psi_i$ と略す)。まず $S_xE_x(\psi_1) = \begin{pmatrix}c + s\,t_1 & \rho s\\ (-s + c\,t_1)/\rho & c\end{pmatrix}$。左から $E_x(\psi_2)$ をかけると 1 行目は変わらず、2 行目は ($t_2/\rho$)×1 行目 + 2 行目:
   - $(1,1)$: $c + s\,t_1 = (c\cos\psi_1 + s\sin\psi_1)/\cos\psi_1 = \cos(\phi - \psi_1)/\cos\psi_1$ ($= M_{11}$)。
   - $(1,2)$: $\rho s$ ($= M_{12}$)。
   - $(2,2)$: $(t_2/\rho)\,\rho s + c = c + s\,t_2 = \cos(\phi - \psi_2)/\cos\psi_2$ ($= M_{22}$)。
   - $(2,1)$: $\dfrac{t_2}{\rho}(c + s\,t_1) + \dfrac{-s + c\,t_1}{\rho} = \dfrac1\rho\bigl[c\,(t_1 + t_2) - s\,(1 - t_1t_2)\bigr]$。ここで $t_1 + t_2 = \dfrac{\sin(\psi_1 + \psi_2)}{\cos\psi_1\cos\psi_2}$、$1 - t_1t_2 = \dfrac{\cos(\psi_1 + \psi_2)}{\cos\psi_1\cos\psi_2}$ ($\tan A + \tan B = (\sin A\cos B + \cos A\sin B)/(\cos A\cos B)$、$1 - \tan A\tan B = (\cos A\cos B - \sin A\sin B)/(\cos A\cos B)$ の通分) なので、括弧 $= \dfrac{c\sin(\psi_1 + \psi_2) - s\cos(\psi_1 + \psi_2)}{\cos\psi_1\cos\psi_2} = \dfrac{\sin(\psi_1 + \psi_2 - \phi)}{\cos\psi_1\cos\psi_2} = -\dfrac{\sin\omega}{\cos\psi_1\cos\psi_2}$ ($= \rho M_{21}$)。

4 成分とも一致する (sympy で各式の差 0: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (b)、[bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (d)):

$$
\begin{pmatrix}x\\ p_x\end{pmatrix}_{\rm out} = E_x(\psi_2)\,S_x\,E_x(\psi_1)\begin{pmatrix}x\\ p_x\end{pmatrix},\quad
E_x(\psi) = \begin{pmatrix}1 & 0\\ \tan\psi/\rho & 1\end{pmatrix},\quad
S_x = \begin{pmatrix}\cos\theta & \rho\sin\theta\\ -\sin\theta/\rho & \cos\theta\end{pmatrix}
\tag{5.25}
$$

($\theta = \mathrm{ANGLE}$、$\rho_b = \rho_0 = \rho$ のとき。) この展開は sympy で行い、(5.20)〜(5.23) と出口の直進をつないだ写像の 2×2 ヤコビアンが (5.25) の積に 30 桁で一致する ($\psi_1 = 0.03$、$\psi_2 = 0.07$、$\phi = 0.1$、$\rho = 10$ m で差 $10^{-167}$。設計軌道は $x = p_x = 0$ に戻る): [bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (g)。水平の端の効果は発散 ($\Delta p_x = +\tan\psi\,x/\rho$、$\psi > 0$ は、面が $x > 0$ (曲がりの外側。ANGLE > 0 は $-x$ へ曲がる) の側で下流へ進む向き (§5.2)。矩形磁石 (E1 = E2 = 0.5) は両端とも $\psi > 0$)。

**垂直**は、Maxwell フリンジの $\Delta p_y = -(p_x/p_v)\,y/\rho$ (5.33) に、面の座標系の $p_{x1} = \sin\psi$、$p_v = \cos\psi$ ($p = 1$) を入れた $-\tan\psi\,y/\rho$ (収束)。垂直の 2×2 は $E_y(\psi_2)\,D(L)\,E_y(\psi_1)$、$E_y(\psi) = \begin{pmatrix}1 & 0\\ -\tan\psi/\rho & 1\end{pmatrix}$。扇形・矩形 (E1 = E2 = 0.5)・非対称 (E1 = 0.3、E2 = 0.7) で、追跡のヤコビアンの水平・垂直とも $5 \times 10^{-13}$ で一致する ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (e))。矩形では水平の $R_{21} = 0$、垂直の $R_{43} = -0.009983$ (ANGLE = 0.1、$\rho = 10$ m)。

## 6. 端部 (フリンジ)

### 6.1 端の磁場

端の近くの磁場 ($x$ 方向に一様な直線の端。$s$ は面の法線方向) を $g(s) = B_y(s)/B_0$ ($g$ は 0 から 1 へ立ち上がる) と、$\nabla\cdot\vec b = 0$、$\nabla\times\vec b = 0$ から作る。磁気スカラーポテンシャル $V$ ($\vec b = \nabla V$、$\nabla^2V = 0$) を

$$
V = \frac{1}{\rho_b}\left[y\,g - \frac{y^3}{6}g'' + \frac{y^5}{120}g'''' - \cdots\right]
\tag{5.30}
$$

に取る ($x$ によらない)。$\vec b = \nabla V$ から $b_x = 0$ で、$b_y = \partial V/\partial y$、$b_s = \partial V/\partial s$ を各項で微分すると

$$
b_x = 0,\qquad b_y = \frac{1}{\rho_b}\left(g - \frac{y^2}{2}g'' + \frac{y^4}{24}g'''' - \cdots\right),\qquad b_s = \frac{1}{\rho_b}\left(y\,g' - \frac{y^3}{6}g''' + \frac{y^5}{120}g^{(5)} - \cdots\right)
\tag{5.31}
$$

$\nabla\times\vec b = 0$ は $\vec b = \nabla V$ なので自動的に成り立つ。$\nabla^2V = 0$ の確認 (各次数が相殺する): (3 項で打ち切った (5.30) について) $\partial_y^2V = \dfrac1{\rho_b}\left[-y\,g'' + \dfrac{y^3}6g''''\right]$、$\partial_s^2V = \dfrac1{\rho_b}\left[y\,g'' - \dfrac{y^3}6g'''' + \dfrac{y^5}{120}g^{(6)}\right]$。$y$、$y^3$ の項は打ち消し合い、$y^5$ の項は $V$ に $-y^7g^{(6)}/5040$ を足せば ($\partial_y^2$ が $-y^5g^{(6)}/120$) 打ち消える (3 項で打ち切った $V$ では $\nabla\cdot\vec b = y^5g^{(6)}/(120\rho_b)$ が残る。sympy: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (f)。$y^5$ は 5 次なので、この章で使う $y^3$ までには効かない)。力は $d\vec P/ds = \vec P\times\vec b/P_z$ から。外積 $\vec P\times\vec b = (p_yb_s - p_zb_y,\ -p_xb_s,\ p_xb_y)$ ($\vec b = (0, b_y, b_s)$) の $x$・$y$ 成分を $P_z$ で割って

$$
\frac{dp_x}{ds} = \frac{p_yb_s}{p_z} - b_y,\qquad
\frac{dp_y}{ds} = -\frac{p_xb_s}{p_z}
\tag{5.32}
$$

($\nabla\cdot\vec b$ は [bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (a) でも確認。) 以下は (5.31) の最低次、$1/\rho_b$ の 1 次で導く (曲率の効果は端面の幾何 (§5) に入れてある)。

### 6.2 硬い端の Maxwell フリンジ

**約束**: 硬い端は、幅 $d$ の端 $g_d(s)$ ($g_d \to \Theta(s)$、$g_d' \to \delta(s)$、$g_d'' \to \delta'(s)$) の $d \to 0$ の極限。端を通る間の変化を $1/\rho_b$ の 1 次で求め、$d$ に比例して消える項 (端の中を進む長さの分) は落とし、$d$ によらない項を残す。$d$ とともに発散する項 ($\int g'^2ds \sim 1/d$ など、$1/\rho_b^2$ の量) は、有限の幅をもつ端の効果 (線形フリンジ F1、§6.3) に属する。

**(i) $p_y$ の変化**。(5.32) の $p_y$ の式に (5.31) の $b_s = (y\,g' - y^3g'''/6)/\rho_b$ を入れ、$p_x/p_z$ が端の中で一定とすると、$\int g'\,ds = 1$、$\int g'''\,ds = [g''] = 0$ なので

$$
\Delta p_y = -\frac{p_x}{p_v}\,\frac{y}{\rho_b},\qquad p_v \equiv \sqrt{p^2 - p_x^2}
\tag{5.33}
$$

($p_z$ のかわりに、$p_y$ を入れない $p_v$ で書いた。K64 の名前。**この近似の大きさ**: $p_z = \sqrt{p_v^2 - p_y^2}$ なので $1/p_z = p_v^{-1}(1 + p_y^2/(2p_v^2) + \cdots)$ で、(i) の $\Delta p_y$ は相対 $p_y^2/(2p_v^2)$ だけ小さく見積もる ($p_y = 1$ mrad で $5 \times 10^{-7}$、10 mrad で $5 \times 10^{-5}$。数値: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (e))。ただし端の磁場の $y^4$ の項から来る補正 $f_4$ ([06_edge_fringe.md](06_edge_fringe.md) §4.3) を一緒に入れると、$f = -\dfrac{y^2}{2\rho_b}\dfrac{p_x}{p_v}\left(1 - \dfrac{p_y^2}{2p_v^2}\right) + O(p_y^4)$ となり、$p_y$ の補正は符号が逆の $-p_y^2/(2p_v^2)$ (大きさは同じ) になる。K64 とヘルプの式 (5.36) は、この $p_y$ の補正を入れていない。)

**(ii) $p_x$ のスパイク**。(5.32) の $p_x$ の式で $p_y = 0$ とし、$b_y = (g - y^2g''/2)/\rho_b$ を入れると、$dp_x/ds = -g/\rho_b + y^2g''/(2\rho_b)$。積分すると

$$
p_x(s) = p_{x0} + \kappa\,g'(s) - \frac{G(s)}{\rho_b},\qquad \kappa \equiv \frac{y^2}{2\rho_b},\quad G(s) = \int^sg\,ds
\tag{5.33a}
$$

最後の項は曲がり (本体の円弧と F1) の効果、$\kappa\,g'(s)$ は端で $p_x$ が一瞬だけ大きく振れる**スパイク** (高さ $\sim\kappa/d$、幅 $d$、積分は $\kappa$)。

**(iii) 座標 $x$ のずれ**。$dx/ds = F(p_x)$、$F(p_x) \equiv p_x/\sqrt{p^2 - p_x^2}$ (= $p_x/p_z$)。(5.33a) のスパイクを $\kappa$ の 1 次で展開する: $F(p_{x0} + \kappa g') = F(p_{x0}) + F'(p_{x0})\,\kappa g' + O(\kappa^2g'^2)$。積分すると、第 1 項は $F(p_{x0})\times(\text{幅})\to0$、第 2 項は $F'(p_{x0})\,\kappa\int g'ds = \kappa F'(p_{x0})$、第 3 項は $\int g'^2ds \sim 1/d$ で発散する (上の約束により F1 の側)。したがって

$$
\Delta x = \kappa\frac{dF}{dp_x} = \frac{y^2}{2\rho_b}\frac{p^2}{p_v^3}
$$

($dF/dp_x = 1/p_v + p_x^2/p_v^3 = p^2/p_v^3$)。$z$ も同様: $z' = -p/p_z$ (経路長の部分。(5.21) と同じ理由) で、$d(1/p_v)/dp_x = p_x/p_v^3$ から $\Delta z = -p\,\kappa\,p_x/p_v^3 = -\dfrac{y^2}{2\rho_b}\dfrac{p_xp}{p_v^3}$。$y$ は、$p_y$ が端の中で変わる量が小さい ($\Delta p_y\sim y/\rho_b$) ので $\Delta y \sim (\Delta p_y/p)\times d \to 0$ で動かない。$p_x$ も、(5.33a) のスパイクは端を通りすぎると 0 に戻るので、硬い端の前後で変わらない。これらをまとめて、運動量 $p_x$ と座標 $y$ を保つ**正準変換**として書く。

**生成関数**。この章の $\exp(:f:)$ は C 付録の約束 ($\Delta x = -\partial f/\partial p_x$、$\Delta p_y = +\partial f/\partial y$、$\Delta z = -\partial f/\partial p$) で書く。(5.33) の $\Delta p_y = +\partial f/\partial y$ から

$$
f_{\rm in} = -\frac{y^2}{2\rho_b}\,\frac{p_x}{p_v}
\tag{5.34}
$$

($y$ に依らない項 $g(p_x, p)$ は、$y = 0$ で座標が動かないことから 0)。$f$ は $x$、$p_y$、$z$ を含まず、$p_x$、$y$、$p$ だけの関数なので、C 付録の括弧 ($[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg)$) で $:f:\,x = -\partial f/\partial p_x$、$:f:\,p_y = +\partial f/\partial y$、$:f:\,z = -\partial f/\partial p$ で、$:f:$ をもう一度かけると 0 (これらの量は $x$、$p_y$、$z$ を含まない)。つまり $\exp(:f:)$ の級数は 1 項で止まり、**厳密に**閉じた形になる ($y$ と $p_x$ は変わらない):

$$
\Delta x = -\frac{\partial f}{\partial p_x} = \frac{y^2}{2\rho_b}\,\frac{p^2}{p_v^3},\qquad
\Delta p_y = +\frac{\partial f}{\partial y} = -\frac{y}{\rho_b}\,\frac{p_x}{p_v},\qquad
\Delta z = -\frac{\partial f}{\partial p} = -\frac{y^2}{2\rho_b}\,\frac{p_xp}{p_v^3}
\tag{5.35}
$$

($\partial(p_x/p_v)/\partial p_x = p^2/p_v^3$、$\partial(p_x/p_v)/\partial p = -p_xp/p_v^3$。(iii) の $\Delta x$・$\Delta z$ と一致。C 付録の括弧での計算を sympy で確認: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (g)。) 出口は $g$ が 1 から 0 へ落ちるので $\int g'\,ds = -1$ となり、$f_{\rm out} = -f_{\rm in}$ (出口では $\Delta x$、$\Delta p_y$、$\Delta z$ の符号が逆)。(5.34) は、$p_y$ を入れない $p_v$ で書いた近似。ヘルプと K64 の式は、さらに $y^4$ の補正を入れる:

$$
f = -\frac{1}{2\rho_b}\left(y^2 - \frac{y^4}{12\rho_b^2}\right)\frac{p_x}{\sqrt{p^2 - p_x^2}},\quad
\Delta x = \frac{b\,y^2}{\rho_b}\,\frac{p^2}{p_v^3},\ \ \Delta p_y = -\frac{a\,y}{\rho_b}\,\frac{p_x}{p_v},\ \ \Delta z = -\frac{b\,y^2}{\rho_b}\frac{p_xp}{p_v^3}
\tag{5.36}
$$

$a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$、$p_v = \sqrt{p^2 - p_x^2}$ (ヘルプの $(1 - y_1^2/(12\rho_b^2))$ と $(1 - y_1^2/(6\rho_b^2))$ の因子に当たる)。

**$a$ と $b$ の整合 (確認済み)**: 上の $f$ から $\Delta p_y = \partial f/\partial y = -\dfrac1{2\rho_b}\left(2y - \dfrac{4y^3}{12\rho_b^2}\right)\dfrac{p_x}{p_v} = -\dfrac{y}{\rho_b}\left(1 - \dfrac{y^2}{6\rho_b^2}\right)\dfrac{p_x}{p_v}$ で $a$ が出る。$\Delta x = -\partial f/\partial p_x = \dfrac1{2\rho_b}\left(y^2 - \dfrac{y^4}{12\rho_b^2}\right)\dfrac{p^2}{p_v^3} = \dfrac{y^2}{\rho_b}\left(\dfrac12 - \dfrac{y^2}{24\rho_b^2}\right)\dfrac{p^2}{p_v^3}$ で $b$ が出る。$\Delta z$ も同様 (sympy で 0: (g))。つまり 3 つの式は 1 つの $f$ から出ていて、$y^4$ の係数 $-1/(24\rho_b^3)$ ($f$ の中の $-\tfrac{1}{2\rho_b}\cdot\tfrac{1}{12\rho_b^2}$) の 1 つで $a$ と $b$ が決まる。

**$y^4$ の係数の導出ができない理由** (何が足りないか)。
1. **1 次の摂動は $y$ について厳密**: (i) の $\int b_s\,ds = \dfrac1{\rho_b}\bigl[V\bigr]$ は、(5.30) の $V$ の端での値の差で、$g$ が $0 \to 1$ に立ち上がるとき $y\,g$ の項だけが残り ($g''$、$g''''$ は端の外で 0)、$y/\rho_b$ ちょうどになる。つまり (5.33) の $\Delta p_y = -(p_x/p_v)\,y/\rho_b$ は、$1/\rho_b$ の 1 次では $y$ の高次の補正を持たない。$a$・$b$ の補正 $y^2/\rho_b^2$ は、$\Delta p_y$ では $y^3/\rho_b^3$、$\Delta x$ では $y^4/\rho_b^3$ で、**$1/\rho_b$ の 3 次**の量。
2. **3 次の摂動に要るもの**: 3 次まで進めると、端の形 (幅 $d$) について $\int g'^2ds \sim 1/d$ のような $d$ とともに発散する積分と、端の中を進む長さ $\sim d$ の分の $y$・$p_y$ の変化 (消える量) の積が、$d$ によらない項として残る。さらに (5.31) は $y^2g''/g \sim y^2/d^2$ を展開パラメータとするべき級数なので、$d \to 0$ の極限は $y/d \to \infty$ で、この級数の収束域の外にある。硬い端の極限で $y^4/\rho_b^3$ の係数を一意に定めるには、(a) $d$ が有限の端の模型で $1/\rho_b^3$ までの摂動 (または厳密な軌道) を実行し、その $d \to 0$ の極限を取ること、(b) 発散する量と消える量の積のうち、どこまでを硬い端の側に、どこを F1 の側に割り振るかの約束 (ヘルプ・K64 が採った約束。資料に無い) が要る。(b) の約束が $d$ の形に依存するかどうかは、確認していない。
3. **確かめられた対応 (導出ではない)**: $a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$ は、$\theta = y/\rho_b$ として $\sin\theta/\theta = 1 - \theta^2/6 + \cdots$、$(1 - \cos\theta)/\theta^2 = \tfrac12 - \theta^2/24 + \cdots$ の $\theta^2$ までの展開に等しい ([06_edge_fringe.md](06_edge_fringe.md) §8.4。$\Delta p_y = -(p_x/p_v)\sin(y/\rho_b)$、$\Delta x = \rho_b(1 - \cos(y/\rho_b))\,p^2/p_v^3$ の形。sympy)。K64 には $\theta^4$ 以上が無い。この三角関数の形が、有限の曲率半径の円弧の幾何から出るのかどうかは確認していない。
4. **数値での確認も届かない**: 検証した端の磁場 (§13 の条件) では補正が $y^2/\rho^2 \lesssim 10^{-9}$ 倍で、ODE の精度 (相対 $10^{-4}$) に見えない。

ここではヘルプ・K64 の係数をそのまま使う。

**K64 のコードとの対応**: 追跡の `tbendcore` (`tbend.f:458-462`、$q$ の単位で `fa`、`f`、`ff`) と CALC の `tbedgemaxwell` (`tbfrie.f:190-229`、正準運動量) は、(5.36) を書いたもの。

| (5.36) | `tbendcore` (q の単位) | `tbedgemaxwell` (正準) |
|---|---|---|
| $\Delta p_y$ | `py2 = py1 - px1*f` ($f = a\,y/(\rho_ep_v/p)$、`fa = y1/rhoe*sqrt(pv1sqi)`、`pv1sqi = 1/(1-px1**2)`) | `cod(4) - af*pxi` ($\texttt{af} = a\,y/\rho_b/p_v$) |
| $\Delta x$ | `x2 = x1 + ff`、`ff = .25*(f+fa)*y1*pv1sqi` ($\simeq b\,y^2p^2/(\rho_bp_v^3)$) | `cod(1) + bf*pr` ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$) |
| $\Delta z$ | `z2 = z1 - px1*ff` | `cod(5) - bf*pxi` |

**ヘルプの式と `tbendcore` の式の形の違い**: ヘルプ (nonlinear fringe at entrance/exit。TeX 版・平文版とも同じ式、**原文で確認**) は (5.36) の形 (入口 $x_2 = x_1 + y_1^2\left(1 - \dfrac{y_1^2}{12\rho_b^2}\right)\dfrac{p_1^2}{2\rho_b(p_1^2 - p_{x1}^2)^{3/2}}$、$p_{y2} = p_{y1} - \cdots$、$z_2 = z_1 - \cdots$) で書く。`tbendcore` は $q = p_x/p$ の単位で、$\Delta p_y$ に使う `f` ($= a\cdot\texttt{fa}$) と、$\Delta x$ の `ff` ($= \tfrac14(f + fa)\,y\cdot\texttt{pv1sqi}$) を、同じ `fa` ($= y/(\rho_e\sqrt{1 - q_x^2})$) から作る。2 つは**代数的に等しい**: $f + fa = (2 - y^2/(6\rho_b^2))\,fa = 2\left(1 - \dfrac{y^2}{12\rho_b^2}\right)fa$ なので $ff = \tfrac12\left(1 - \dfrac{y^2}{12\rho_b^2}\right)fa\,y/(1 - q_x^2) = \dfrac{y^2}{2\rho_bp}\left(1 - \dfrac{y^2}{12\rho_b^2}\right)(1 - q_x^2)^{-3/2}$。$p^2/p_v^3 = 1/(p(1 - q_x^2)^{3/2})$ を使うと、ヘルプの $\Delta x$ に一致する (展開ではなく恒等式)。$\Delta p_y$ も $-q_xf\cdot p = -(p_x/p_v)\,a\,y/\rho_b$ で (5.36) と同じ。`(f + fa)` の対称な書き方は、ヘルプの $(1 - y^2/(12\rho_b^2))$ を $f$ の $(1 - y^2/(6\rho_b^2))$ から作るための式になっている。この形にした意図は、コードからは確認できない。

**検証 1: 端の磁場の中の運動** (直線座標、$1/\rho = 0.02$ /m、$y = 2$ mm、$q_x = p_x/p = 0.1$、端の形 $g = (1 + \tanh(s/d))/2$ で (5.31) の磁場、DOP853。硬い端との差を取り、$d = 2, 1$ mm から $d \to 0$ へ外挿): $\Delta x$ は $\delta = 0$ と $\delta = 0.3$ のどちらも (5.35) の $p^2/p_v^3$ の式と相対 $10^{-4}$ 以下で一致 (比 1.0000)。**CALC の $p$ の因子 (`bf*pr`) だけにした式とは $\delta = 0.3$ で比 1.3000 (= $p$) ずれる**。$\Delta p_y$ は (5.33) と $\delta = 0$・0.3 で相対 $1.5 \times 10^{-4}$、$1.2 \times 10^{-4}$ ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。

**検証 2: 追跡の式のシンプレクティック性**: 追跡の式 (`f`・`ff`) の 6×6 ヤコビアン (中心差分) の $|J^\top SJ - S|$ は $\delta = 0, \pm0.1$ で $10^{-12}$。CALC の `tbedgemaxwell` の COD の更新は $2.0 \times 10^{-5}$ ($\delta = 0$)、$1.7 \times 10^{-4}$ ($\delta = 0.1$)、$2.5 \times 10^{-4}$ ($\delta = -0.1$) で、シンプレクティックでない (K64 の誤り。v2 は直してあり、COD の更新の数値ヤコビアンが解析行列と $10^{-12}$ で一致し、シンプレクティック ($1.5 \times 10^{-12}$。[bend_03_fringe.v2.out.txt](outputs/bend_03_fringe.v2.out.txt) (c)、[edge_03_bend_maxwell_recheck.v2.out.txt](outputs/edge_03_bend_maxwell_recheck.v2.out.txt) (c))。付録 B-1。[bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (c) は k64)。

### 6.3 線形フリンジ (F1)

端の磁場が有限の幅で立ち上がるとき、硬い端 ($\Theta$) との違いは、端の形 $g(s) = B_y(s)/B_0$ のいくつかの積分で決まる。ヘルプと K64 は、これを F1 という**1 つの長さ** (ヘルプ: "Length of the slope of the field at the edge") で表す。ヘルプの厳密な定義は

$$
\mathrm{F1} = 6\int_{-\infty}^{\infty}\bigl(g(s) - g(s)^2\bigr)\,ds\qquad(\text{1 つのフリンジについて})
\tag{5.40}
$$

を使い、(5.31)〜(5.32) の入射 $p_x = 0$ の粒子で、硬い端との差を導く (端は直線座標、1 次の $1/\rho_b$ で $p_x$ は $p_x(s) = -G(s)/\rho_b$、$G(s) = \int^sg$。$p$ は全運動量で、$p_z \simeq p$):

まず (5.32) の $p_x$ の式 ($p_y \simeq 0$) を (5.33a) と同じ形に積分する: $p_x(s) = -G(s)/\rho_b + \dfrac{y^2}{2\rho_b}g'(s)$ ($G = \int^sg\,ds$。$y^2g''/2$ の項の積分がスパイク、$G$ の項は本体の曲がり)。ここでスパイクの項は硬い端 (§6.2) では $\Delta x$ などの $d$ によらない項を作ったが、有限の幅の端では、スパイクどうし・スパイクと $G$ の積が $1/\rho_b^2$ の量として残る。それを以下で求める ($p_z \simeq p$)。

1. **$y$ の項**: $dp_y/ds = -p_xb_s/p_z$ に $p_x = -G/\rho_b$、$b_s = y\,g'/\rho_b$ を入れて $dp_y/ds = (y/(\rho_b^2p))\,G\,g'$。端の領域を $[s_a, s_b]$ ($g(s_a) = 0$、$g(s_b) = 1$。区間の外では $g'$ が 0) に取り、入射で $p_x = 0$ なので $G(s_a) = 0$ とする。部分積分 $\int Gg'\,ds = [Gg]_{s_a}^{s_b} - \int g\cdot g\,ds$ ($G' = g$ を使った)。境界項は $G(s_b)\cdot1 - 0 = \int_{s_a}^{s_b}g\,ds$ なので、$\int Gg'\,ds = \int(g - g^2)ds = \mathrm{F1}/6$ ((5.40)。$g - g^2$ は端の外側で 0 なので積分区間は無限に取れる)。したがって $\Delta p_y = y\,\mathrm{F1}/(6\rho_b^2\,p)$。
2. **$y^3$ の項**: $b_s$ の $-y^3g'''/(6\rho_b)$ からの寄与は、$dp_y/ds \ni -\dfrac{p_x}{p}\left(-\dfrac{y^3g'''}{6\rho_b}\right) = -\dfrac{y^3}{6\rho_b^2p}\,G\,g'''$ で、$\int Gg'''\,ds = [Gg''] - \int g\,g''\,ds = -[g\,g'] + \int g'^2ds = \int g'^2ds$ (部分積分 2 回。$g'$、$g''$ は端の外側で 0 なので境界項は消える) から $-\dfrac{y^3}{6\rho_b^2p}\int g'^2ds$。さらに、$p_x$ のスパイクの項 $+\dfrac{y^2}{2\rho_b}g'$ が $dp_y/ds = -p_xb_s/p_z$ に $-\dfrac{y^2g'}{2\rho_b p}\cdot\dfrac{y\,g'}{\rho_b} = -\dfrac{y^3}{2\rho_b^2p}g'^2$ を足す。合計 $-\left(\tfrac16 + \tfrac12\right)\dfrac{y^3}{\rho_b^2p}\int g'^2ds = -\dfrac23\,\dfrac{y^3}{\rho_b^2p}\int g'^2ds$。
3. **$x$ の項**: $dx/ds = p_x/p_z \simeq p_x/p = -G/(\rho_bp)$ (スパイクの項は §6.2 で扱い済み)。硬い端の $G_{\rm hard}(s) = (s - s_0)\Theta(s - s_0)$ との差 $D(s) \equiv G - G_{\rm hard} = \int^s(g - \Theta)ds'$ ($s_0$ は $I_0 \equiv \int(g - \Theta)ds = 0$、つまり $D$ が端の外側で 0 になる位置) を使うと $\Delta X = -\dfrac{1}{\rho_bp}\int D\,ds$。部分積分 $\int D\,ds = [sD] - \int sD'\,ds = -\int s(g - \Theta)ds$ (境界項は $D \to 0$ で消える) より

$$
\Delta X = \frac{I_1}{\rho_bp},\qquad I_1 = \int(s - s_0)\bigl(g - \Theta(s - s_0)\bigr)\,ds
$$

**近似の次数と範囲**: 上の 1〜3 の $1/p$ 依存 ($\Delta p_y \propto 1/p$、$\Delta X \propto 1/p$) は、$p_z \simeq p$ の近似から来る。入射が $p_x = 0$ で、端の中の $p_x(s) = -G(s)/\rho_b$ の大きさは $|G|/\rho_b \lesssim \mathrm{F1}/\rho_b$ (例: F1 = 0.1 m、$\rho_b = 10$ m で $10^{-2}$) なので、$p_z = p\sqrt{1 - p_x^2/p^2 - \cdots}$ の相対誤差は $p_x^2/(2p^2) \lesssim 5 \times 10^{-5}$、すなわち $(\mathrm{F1}/\rho_b)^2$ の次数と、$p_y^2/(2p^2)$ (§6.2 の (5.33) と同じ)。**面が傾いていて斜めに入射する**とき ($p_x$ が $p\sin\psi$ で、$p_z \simeq p$ が成り立たない) は、$\Delta X$ に因子 $1/\cos^3\psi$ が付く ((5.43a)、§6.4。数値で確認)。$\Delta p_y$ の $y$・$y^3$ の項の斜め入射の因子は導出していない。K64 の `dxfr1`・`dyfr1`・`dyfra1` は、これらの因子を持たない。

端が**直線の立ち上がり (長さ $\ell$)** ($g = s/\ell + \tfrac12$、$-\ell/2 < s < \ell/2$) のとき: $\int g(1 - g)ds = \int_{-\ell/2}^{\ell/2}\left(\tfrac14 - \tfrac{s^2}{\ell^2}\right)ds = \tfrac\ell4 - \tfrac{\ell}{12} = \ell/6$、$\int g'^2ds = \ell\cdot\ell^{-2} = 1/\ell$、$I_0 = 0$ ($s_0 = 0$) で $I_1 = \int_{-\ell/2}^0s\left(\tfrac s\ell + \tfrac12\right)ds + \int_0^{\ell/2}s\left(\tfrac s\ell - \tfrac12\right)ds = -\ell^2/48 - \ell^2/48 = -\ell^2/24$ (sympy: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (h))。つまり $\mathrm{F1} = \ell$ で、3 つの量が

$$
\Delta p_y = \frac{y}{p}\left(\frac{\mathrm{F1}}{6\rho_b^2} - y^2\,\frac{2}{3\,\mathrm{F1}\,\rho_b^2}\right),\qquad
\Delta X = -\frac{\mathrm{F1}^2}{24\rho_bp}
\tag{5.41}
$$

になる (ヘルプの $\Delta y_{fr} = \mathrm{F1}/(6\rho_b^2)$、$\Delta y_{fra} = 2/(3\,\mathrm{F1}\,\rho_b^2)$)。**$\tanh$ の端** ($g = (1 + \tanh(s/d))/2$): $g(1 - g) = \tfrac14\,\mathrm{sech}^2(s/d)$、$\int\mathrm{sech}^2(s/d)\,ds = 2d$ より $\int g(1 - g)ds = d/2$ で $\mathrm{F1} = 6\cdot d/2 = 3d$。$g' = \mathrm{sech}^2(s/d)/(2d)$、$\int\mathrm{sech}^4(s/d)\,ds = 4d/3$ より $\int g'^2ds = \dfrac{1}{4d^2}\cdot\dfrac{4d}{3} = \dfrac1{3d} = \dfrac1{\mathrm{F1}}$ で、$y^3$ の項は直線の立ち上がりと同じ形になる。$x$ の項は、$g - \Theta$ が $s$ の奇関数 ($g(-s) = 1 - g(s)$) なので $s_0 = 0$、$s(g - \Theta)$ は偶関数で $I_1 = -2\int_0^\infty s\,(1 - g)\,ds = -d^2\int_0^\infty u\,(1 - \tanh u)\,du$ ($u = s/d$)。$1 - \tanh u = 2/(e^{2u} + 1)$ から $\int_0^\infty\dfrac{2u}{e^{2u} + 1}du = \tfrac12\int_0^\infty\dfrac{v}{e^v + 1}dv = \tfrac12\cdot\dfrac{\pi^2}{12} = \dfrac{\pi^2}{24}$ なので $I_1 = -\pi^2d^2/24$。(5.41) の $\mathrm{F1}^2/24 = 9d^2/24$ との比は $\pi^2/9 = 1.097$ (数値: $1.096623$)。**誤差関数の端** ($g = (1 + \mathrm{erf}(s/(\sqrt2\sigma)))/2$、$\sigma$ = 2 mm): $\int g'^2ds = 1/(2\sqrt\pi\,\sigma)$、$\int g(1 - g)ds = \sigma/\sqrt\pi$ (数値積分で確認) から $\mathrm{F1} = 6\sigma/\sqrt\pi$、$\int g'^2ds\cdot\mathrm{F1} = 3/\pi = 0.955$ (-4.5%) ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (d)、[bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (h))。**ヘルプ・K64 の式は、1 つの F1 から 3 つの量を作る (直線の立ち上がりで厳密、一般の端の形では近似)**。

**設計軌道との関係**: $p = 1$ のときの $\Delta X = -\mathrm{F1}^2/(24\rho_b)$ は、$L'$ の補正 (§6.4) と $\rho_b = L'/(\mathrm{ANGLE} + \mathrm{K0})$ の定義で設計軌道に吸収される。$p \neq 1$ のずれの差 $\Delta X(p) - \Delta X(1) = \dfrac{\mathrm{F1}^2}{24\rho_b}\left(1 - \dfrac{1}{p}\right) = \mathrm{dxfr}\cdot\dfrac{\delta}{p}$ だけが写像に入る。

**生成関数と閉じた写像**: (5.41) の 3 つの効果は、$p_x$ と $y$ と $p$ だけの関数の生成関数で表せる。

$$
f_{\rm lin} = -p_x\,\mathrm{dxfr}\left(1 - \frac1p\right) + \frac{1}{p}\left(\frac{\mathrm{dyfr}}{2}y^2 - \frac{\mathrm{dyfra}}{4}y^4\right),\quad
\mathrm{dxfr} = \frac{\mathrm{F1}^2}{24\rho_b},\ \mathrm{dyfr} = \frac{\mathrm{F1}}{6\rho_b^2},\ \mathrm{dyfra} = \frac{4\,\mathrm{dyfr}}{\mathrm{F1}^2}
\tag{5.42}
$$

($\mathrm{dyfra} = 2/(3\,\mathrm{F1}\,\rho_b^2)$。$f$ の符号は §6.2 と同じ C 付録の約束。) C 付録の約束 ($\Delta x = -\partial f/\partial p_x$、$\Delta p_y = +\partial f/\partial y$、$\Delta z = -\partial f/\partial p$) で、各項を 1 つずつ微分すると:

- $\Delta x = -\partial f/\partial p_x = \mathrm{dxfr}\,(1 - 1/p) = \mathrm{dxfr}\,\delta/p$ ($p_x$ の項だけ)。
- $\Delta p_y = +\partial f/\partial y = (\mathrm{dyfr}\,y - \mathrm{dyfra}\,y^3)/p = (\mathrm{dyfr} - \mathrm{dyfra}\,y^2)\,y/p$ ($\tfrac{\mathrm{dyfr}}2y^2 \to \mathrm{dyfr}\,y$、$\tfrac{\mathrm{dyfra}}4y^4 \to \mathrm{dyfra}\,y^3$)。
- $\Delta z = -\partial f/\partial p$: $\partial_p[-p_x\mathrm{dxfr}(1 - 1/p)] = -p_x\,\mathrm{dxfr}/p^2$、$\partial_p[(\cdots)/p] = -(\cdots)/p^2$ なので $\Delta z = \bigl[\mathrm{dxfr}\,p_x + (\mathrm{dyfr}/2 - \mathrm{dyfra}\,y^2/4)\,y^2\bigr]/p^2$ ($p_x$ は正準運動量。$q$ の単位 $q_x = p_x/p$ では $\bigl[\mathrm{dxfr}\,q_x + (\mathrm{dyfr}/2 - \mathrm{dyfra}\,y^2/4)\,y^2/p\bigr]/p$)。

(sympy で 3 式とも差 0: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (h)。) これが、`tbendcore` の入口の `x1 + dxfr1*dp*pinv`、`py1 + (dyfr1 - dyfra1*y1**2)*y1*pinv**2` ($q$ の単位で $p_y/p$ の更新)、`z1 + (dxfr1*px1 + (.5*dyfr1 - .25*dyfra1*y1**2)*y1**2*pinv)*pinv` ($q$ の単位。`px1` は $q_x$) と一致する。$f$ が $x$、$p_y$、$z$ を含まないので流れは閉じた形で、**シンプレクティック性は構成から厳密**。出口は $\Delta x$ の符号が逆 ($\mathrm{dxfr2} = -\mathrm{F1}^2/(24\rho_b)$ とも書く) で、$\Delta p_y$ の符号は同じ。K64 は、これに端面の傾きによる $z$ の調整 `dzfr1 = dxfr1*sinp1` (入口)・`dzfr2 = dxfr2*sinp2` を足す。

**ヘルプの $V$ との関係**: ヘルプ (平文版) の $V = -\dfrac1p\left(\dfrac{\mathrm{F1}^2}{24\rho_b}p_x + \dfrac{\mathrm{F1}}{12\rho_b^2}y^2 - \dfrac{y^4}{6\,\mathrm{F1}\rho_b^2}\right)$ は、ヘルプの流儀 ($\exp(:V:)$ が $V$ をハミルトニアンとする流れ) で書いてある。この章の約束では $f = -V$ に当たり、$-V = \dfrac{\mathrm{dxfr}\,p_x}{p} + \dfrac1p\left(\dfrac{\mathrm{dyfr}}2y^2 - \dfrac{\mathrm{dyfra}}4y^4\right)$。(5.42) の $f_{\rm lin}$ との差は $-V - f_{\rm lin} = \mathrm{dxfr}\,p_x$ で、$p$ によらない $x$ のずれ ($\Delta x = -\partial(\mathrm{dxfr}\,p_x)/\partial p_x = -\mathrm{dxfr} = -\mathrm{F1}^2/(24\rho_b)$。(5.41) の $p = 1$ の $\Delta X$ そのもの) だけ。この項が、設計軌道に吸収される部分 (上の「設計軌道との関係」、§6.4) に当たる (sympy: (h) の最後の行)。

**検証** (直線座標、$\tanh$ の端 $d = 2$ mm ($\mathrm{F1} = 6$ mm)、$\rho = 50$ m、$y = 0.5, 1$ mm の 2 点から $y$ と $y^3$ の係数を分離): $\delta = 0$ で $\Delta p_y$ の $y$ の係数は $\mathrm{F1}/(6\rho^2)$ の 1.0007 倍、$y^3$ の係数は $\tfrac23/(\mathrm{F1}\rho^2)$ の 1.031 倍 (残りは $y^5$ の項)。$\delta = 0.3$ では係数に $p$ をかけた値が $\delta = 0$ と一致 (すなわち $1/p$)。$x$ のずれ $\Delta X\cdot p/\mathrm{dxfr}$ は、直線の立ち上がりで $\delta = 0, 0.3, -0.2$ のどれも $-1.0000$、$\tanh$ では $-1.0966 = -\pi^2/9$。$\Delta X(\delta) - \Delta X(0)$ と $\mathrm{dxfr}\,\delta/p$ の比は、直線で 1.0000、$\tanh$ で 1.0966 ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (d))。

### 6.4 F1 による長さ補正と、フリンジ領域の分割

F1 があると、磁場が立ち上がる分、基準軌道の幾何を保つために、本体の円弧の長さを

$$
L' = L - \frac{(\phi_0\,\mathrm{F1})^2}{24\,L}\,\frac{\sin\bigl(\tfrac12(\phi_0(1 - \mathrm{E1} - \mathrm{E2}) - \mathrm{AE1} - \mathrm{AE2})\bigr)}{\sin(\phi_0/2)}
\tag{5.43}
$$

に取り直す (ヘルプ。**原文で確認**: $L'$ の式は TeX 版・平文版とも同じ。K64 は入口と出口で別の F1 を許し、$((\phi_0\mathrm{fb1})^2 + (\phi_0\mathrm{fb2})^2)/(48L)$ と書く。`tffs.f:1282`、`tbende.f:384`)。同じ F1 では一致。

**$L'$ の導出 (K0 = 0、$\mathrm{F1}^2$ の 1 次)**。方針: 端の磁場が幅をもって立ち上がるとき、設計軌道 (入口の基準面で $x = p_x = 0$) が硬い端の円弧からどれだけずれるかを (5.41) と同じ $1/\rho_b$ の 1 次で求め、そのずれを吸収する円の半径 ($= L'/\phi_0$) を、円と面の幾何から決める。ヘルプの (5.43) は $\psi_1 = \psi_2 \approx 0$ の場合に当たり、以下では $\psi_1 \neq \psi_2$ と斜め入射にまで広げて、ヘルプの式との関係を調べる。記号: $\phi_0 = \mathrm{ANGLE}$、$s_0 \equiv -I_1/\rho_b$ ($I_1 = \int(s - s_0')(g - \Theta)\,ds$。直線の立ち上がりで $s_0 = \mathrm{F1}^2/(24\rho_b)$、(5.41))、$\hat e_{x1}$・$\hat e_{f2}$・$\hat E$ は §5.3 と同じ ($\hat E = (-\sin\tfrac{\phi_0}2,\ \cos\tfrac{\phi_0}2)$ は弦の向き)。

**(A) 入口の端のずれ (斜め入射)**。入口の面の座標系 ($s$ が面の法線、$x$ が面に沿う向き) で、設計軌道は $p_x = p\sin\psi_1$、$p_z = p\cos\psi_1$ ($p = 1$ なら $\sin\psi_1$、$\cos\psi_1$) で入る ((5.21))。$y = p_y = 0$ とし、(5.32) を使う:
1. $s$ を独立変数にしているので ($ds = p_z\,d\tau$)、$dp_x/ds = -b_y = -g(s)/\rho_b$ に $p_z$ は現れない。$dx/ds = F(p_x)$、$F(p_x) \equiv p_x/\sqrt{p^2 - p_x^2}$。
2. 積分して $p_x(s) = p_{x0} - G(s)/\rho_b$ ($G = \int^sg$)。硬い端 ($G_{\rm hard} = (s - s_0')\Theta(s - s_0')$) との差を $D(s) \equiv G - G_{\rm hard}$ とおくと、$p_x = p_{x,{\rm hard}} - D/\rho_b$。$D$ は端の外で 0 (§6.3 の手順 3)。
3. $x$ のずれ: $\Delta X = \int\bigl[F(p_{x,{\rm hard}} - D/\rho_b) - F(p_{x,{\rm hard}})\bigr]ds$。$1/\rho_b$ の 1 次を取り、$D$ が非零の端の中 (幅 $\sim\mathrm{F1}$) では $p_{x,{\rm hard}} \simeq p_{x0}$ とみなすと ($p_{x,{\rm hard}}$ の変化は $\mathrm{F1}/\rho_b$ で、さらに 1 桁小さい)、$\Delta X = -\dfrac{F'(p_{x0})}{\rho_b}\displaystyle\int D\,ds$。
4. 部分積分 (§6.3 の手順 3 と同じ) で $\int D\,ds = -I_1$。$F'(p_x) = \dfrac1{\sqrt{p^2 - p_x^2}} + \dfrac{p_x^2}{(p^2 - p_x^2)^{3/2}} = \dfrac{p^2}{p_z^3}$。よって

$$
\Delta X_{\rm in} = \frac{I_1}{\rho_b}\,\frac{p^2}{p_{z}^3} = -\frac{s_0}{p\cos^3\psi_1}\qquad(\text{入口。}p_z = p\cos\psi_1)
\tag{5.43a}
$$

(面に垂直な入射 $\psi_1 = 0$ なら (5.41)。$p = 1$ で $-s_0/\cos^3\psi_1 \equiv -a_1$ と書く。)

5. これは、運動量が端の外で硬い端と同じ ($D \to 0$) で、位置だけが**面に沿って $\Delta X$** ずれた、という結果なので、実際の軌道は硬い端の軌道を $\Delta X\,\hat e_{x1}$ だけ平行移動したもの。一様磁場の円は平行移動に対して閉じているので、本体の円 $C'$ は、硬い端の円 $C_h$ (原点を通り、向き $\hat z = (0, 1)$) を $-a_1\hat e_{x1}$ だけ平行移動したもの。$C'$ は点 $A \equiv -a_1\hat e_{x1}$ を向き $\hat z$ で通る。

数値の確認 (端の形 $g = (1 + \tanh(s/d))/2$、$d = 10$ mm、$p = 1$ と 1.03、$\psi_1 = 0, 0.1, 0.3, 0.6$ rad。$\varepsilon = 10^{-6}$ で上の積分を mpmath で数値積分し、$\varepsilon I_1p^2/p_z^3$ との比を取る): **入口の比は 1.0000 (最大 $10^{-8}$ のずれ)**。1 次の $F'(p_{x0})$ の $1/\cos^3\psi$ (例えば $\psi = 0.6$ で 1.78 倍) がそのまま出る ([bend_08_lprime_general.out.txt](outputs/bend_08_lprime_general.out.txt) (b))。

**(B) 出口の端のずれ (ODE で通す)**。出口の面の座標系で、端の場は入口を裏返したもの: $g_{\rm out}(s) = g(-s)$ ($s$ が増える向きに 1 から 0 へ落ちる)。出口の面に入る向きは $p_x = -p\sin\psi_2$ ($\theta_2 = -\psi_2$。§5.3)。(A) と同じ式 $dp_x/ds = -g_{\rm out}/\rho_b$ を積分する。硬い端との差 $D_{\rm out}(s) \equiv \int^s(g_{\rm out} - \Theta(-s'))\,ds'$ は、$s' \to -s'$ の置き換えと $I_0 = \int(g - \Theta)ds = 0$ から $D_{\rm out}(s) = -D(-s)$。したがって $\int D_{\rm out}\,ds = -\int D\,ds = +I_1$ で、(A) の手順 3〜4 から

$$
\Delta X_{\rm out} = -\frac{F'(p_{x})}{\rho_b}\int D_{\rm out}\,ds = -\frac{I_1}{\rho_b}\frac{p^2}{p_z^3} = +\frac{s_0}{p\cos^3\psi_2}\ \ (\equiv +a_2/p)\qquad(\text{出口。}F'\text{ は }p_x\text{ の偶関数})
\tag{5.43b}
$$

**符号は入口と反対**: 入口では面の手前 ($s < 0$) にも磁場 ($g > 0$) がしみ出しているので、粒子は硬い端より先に曲がり始め、内側 ($-x$) に寄る ($D > 0$)。出口では、面の手前 ($s < 0$) で磁場が全強度に届かない分 ($g_{\rm out} < 1$)、硬い端より曲がり方が小さい期間があり ($D_{\rm out} < 0$)、面を出たあとの位置は外側 ($+x$) に寄る (面の外へしみ出した磁場は、積分した曲がり角を同じに戻すだけ)。つまり**出口を出る直線は、硬い端の直線を $+a_2\hat e_{f2}$ だけ平行移動したもの**で、符号は時間反転の対称性の議論ではなく、(5.32) の積分そのものから出る。数値 (同じ端の形): 出口の $-\Delta X/(\varepsilon I_1p^2/p_z^3)$ は 1.0000 ($\psi_2 = 0.6$ で $3 \times 10^{-7}$ のずれ)。([bend_08_lprime_general.out.txt](outputs/bend_08_lprime_general.out.txt) (b)。)

**(C) 円の幾何**。以下 $p = 1$。$C'$ が出口の**硬い**面と交わる点を $Q$ とする。実際の出口の直線は、$Q$ から出る直線を $+a_2\hat e_{f2}$ 平行移動したものなので、$Q + a_2\hat e_{f2}$ が (面の上の) 基準軌道の出口の点 $E$ に一致すれば、出口の位置は設計どおり。出口の向きが設計どおり ($p_{xr} = 0$) であるためには、$C'$ の円弧の角 $\Theta$ が $\phi_0$ に一致することが要る。$\Theta = \phi_0$ のとき、$A$ から $Q$ への弦は、(5.12) で $P_{X0} = 0$、$P_{Z0} = 1$ とおいた $(X, Z) = \rho'(\cos\phi_0 - 1,\ \sin\phi_0) = 2\rho'\sin\tfrac{\phi_0}2\,\hat E$ (半径 $\rho'$、$\sin\phi_0 = 2\sin\tfrac{\phi_0}2\cos\tfrac{\phi_0}2$、$\cos\phi_0 - 1 = -2\sin^2\tfrac{\phi_0}2$)。設計の $E = 2\rho_0\sin\tfrac{\phi_0}2\,\hat E$ との差を $\Delta C' \equiv 2(\rho' - \rho_0)\sin\tfrac{\phi_0}2$ とおく ($L' - L = \phi_0(\rho' - \rho_0)$ なので $L' - L = \phi_0\,\Delta C'/(2\sin\tfrac{\phi_0}2)$)。すると $Q = A + 2\rho'\sin\tfrac{\phi_0}2\hat E = E - a_1\hat e_{x1} + \Delta C'\,\hat E$。

**条件 (X): $Q$ が出口の面の上にある** ($\Theta = \phi_0$ で円が面に当たる条件): $(Q - E)\cdot\hat N_2 = 0$。内積は $\hat e_{x1}\cdot\hat N_2 = \cos\psi_1(-\sin(\phi_0 - \psi_2)) + \sin\psi_1\cos(\phi_0 - \psi_2) = -\sin\omega$ ($\omega = \phi_0 - \psi_1 - \psi_2$)、$\hat E\cdot\hat N_2 = \sin\tfrac{\phi_0}2\sin(\phi_0 - \psi_2) + \cos\tfrac{\phi_0}2\cos(\phi_0 - \psi_2) = \cos(\tfrac{\phi_0}2 - \psi_2)$。$a_1\sin\omega + \Delta C'\cos(\tfrac{\phi_0}2 - \psi_2) = 0$ から

$$
\Delta C' = -\,a_1\,\frac{\sin\omega}{\cos\left(\tfrac{\phi_0}2 - \psi_2\right)},\qquad
L' - L = \frac{\phi_0\,\Delta C'}{2\sin(\phi_0/2)} = -\,\frac{\phi_0\,s_0}{2\sin(\phi_0/2)}\,\frac{1}{\cos^3\psi_1}\,\frac{\sin\omega}{\cos(\tfrac{\phi_0}2 - \psi_2)}
\tag{5.43c}
$$

**(D) 特別な場合への帰着**。
- **対称な面 ($\psi_1 = \psi_2 = \psi$)**: $\tfrac{\phi_0}2 - \psi = \tfrac\omega2$ なので $\sin\omega/\cos(\omega/2) = 2\sin(\omega/2)$ で、$\Delta C' = -2a\sin(\omega/2)$、$L' - L = -\phi_0\,s_0\,\dfrac{\sin(\omega/2)}{\sin(\phi_0/2)}\dfrac1{\cos^3\psi}$。さらに $Q + a_2\hat e_{f2} - E$ の $\hat e_{\rm ref}$ 成分は 0 になり (条件 (X) で位置も合う。数値は下の表)、**向きと位置が同時に満たされる**。$\psi = 0$ (扇形: $\omega = \phi_0$) では $L' - L = -s_0\phi_0 = -\dfrac{(\phi_0\mathrm{F1})^2}{24L}$ ($s_0 = \mathrm{F1}^2/(24\rho_b)$、$\rho_b = L/\phi_0$)、一般の $\omega$ では $-s_0\phi_0\sin(\omega/2)/\sin(\phi_0/2)$ で、**ヘルプの (5.43) に一致する**。つまり (5.43) は、面に垂直な入射の $\Delta X$ ((5.41) の $-s_0$) を使った $\psi_1 = \psi_2$ の場合 (扇形・対称な面) で、$1/\cos^3\psi$ を除いて正しい。矩形 ($\omega = 0$) では $L' = L$。
- **斜め入射の因子**: ヘルプ・K64 は $1/\cos^3\psi$ を持たない。$\psi_1 = \psi_2 = 0.3$ rad で 15% (下の表の 0.495 と 0.568)、$\psi = 0.05$ rad で $3\psi^2/2 \approx 0.4\%$。
- **非対称な面 ($\psi_1 \neq \psi_2$)**: (5.43c) とヘルプの (5.43) ($\propto 2\sin(\omega/2)$) の比は ($1/\cos^3$ を除いて) $\dfrac{\cos(\omega/2)}{\cos(\omega/2 + (\psi_1 - \psi_2)/2)}$ ($\tfrac{\phi_0}2 - \psi_2 = \tfrac\omega2 + \tfrac{\psi_1 - \psi_2}2$)。$\psi_1 - \psi_2$ の 1 次の差 (例: $\psi_1 = -0.05$、$\psi_2 = 0.25$、$\omega = -0.1$ で 1.9%)。
- **弦の長さだけを合わせる見積り**: 向き (条件 (X)) を無視して、弦の $\hat E$ 方向の長さの 1 次の変化だけを $2(\rho' - \rho_0)\sin\tfrac{\phi_0}2$ に等しいとおくと、$\Delta C_\parallel = (a_1\hat e_{x1} - a_2\hat e_{f2})\cdot\hat E = -s_0\bigl[\sin(\tfrac{\phi_0}2 - \psi_1) + \sin(\tfrac{\phi_0}2 - \psi_2)\bigr] = -2s_0\sin\tfrac\omega2\cos\tfrac{\psi_1 - \psi_2}2$ (和積公式。$a_i \to s_0$)。ここに $\cos((\psi_1 - \psi_2)/2)$ の因子が出る。ヘルプの $2s_0\sin(\omega/2)$ ($= |\hat e_{x1} - \hat e_{f2}|\,s_0$。2 点の差のベクトルの大きさそのもの) との比は $\cos\tfrac{\psi_1 - \psi_2}2 \simeq 1 - (\psi_1 - \psi_2)^2/8$ (差 $(\psi_1 - \psi_2)^2/8$。$\psi_1 - \psi_2 = 0.04$ で $2 \times 10^{-4}$、0.2 で $5 \times 10^{-3}$。数値: [bend_08_lprime_general.out.txt](outputs/bend_08_lprime_general.out.txt) (a))。差のベクトル $\hat e_{x1} - \hat e_{f2}$ は弦の向き $\hat E$ に対して角 $(\psi_1 - \psi_2)/2$ だけ傾いていて、$\hat E$ 成分が $\cos$、直交する成分が $\sin\tfrac{\psi_1 - \psi_2}2$ (sympy で確認、(a))。ヘルプの式は**ベクトルの大きさ**を使った形に当たり、$\psi_1 = \psi_2$ でだけ弦の向きの成分と一致する。
- **条件 (X) で残る位置のずれ**: 一般には $x_{\rm out} = \bigl[-a_1\hat e_{x1} + \Delta C'\hat E + a_2\hat e_{f2}\bigr]\cdot\hat e_{\rm ref}$ が 0 にならない。半径 $\rho'$ という 1 つの量では、出口の向きと位置の 2 つを同時には満たせないため。$\psi_1 \neq \psi_2$ (または斜め入射で $a_1 \neq a_2$) の磁石では、物理的な端の磁場の設計軌道は出口で $x_{\rm out} \neq 0$ になる。

**ODE による確認** (実験室系の水平面内の運動。$p = 1$、$\phi_0 = 0.1$、$\rho_0 = 10$ m、$L = 1$ m。端の場 $b_y = g(s_1)\,g(-s_2)/\rho_b$ ($s_1$・$s_2$ は入口・出口の面の法線座標)、$g = (1 + \tanh(s/d))/2$、$F1 = 3d = 0.06, 0.03, 0.015$ m の 3 点から $F1 \to 0$ へ線形外挿。出口の向き $p_{xr} = 0$ になる $\rho_b$ を DOP853 の積分と割線法で求めて $L'$ とし、(5.43c) などと比べる。$L' - L$ と $x_{\rm out}$ は $s_I = -I_1/\rho_b = \pi^2d^2/(24\rho_b)$ を単位とする。tanh の端は $I_1 \neq -\mathrm{F1}^2/24$ ($\pi^2/9$ 倍。§6.3) なので、式の $s_0$ は $s_I$ に置き換えて比べる)。

| $\psi_1$ | $\psi_2$ | $\omega$ | $(L' - L)/s_I$: ODE | (5.43c) | (5.43c) で $a_i = s_0$ | ヘルプ (5.43) | $x_{\rm out}/s_I$ (ODE、$p_{xr} = 0$) |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 0.10 | $-0.10000$ | $-0.10000$ | $-0.10000$ | $-0.10000$ | 0 |
| 0.02 | 0.02 | 0.06 | $-0.06005$ | $-0.06005$ | $-0.06002$ | $-0.06002$ | 0 |
| 0.3 | 0.3 | $-0.50$ | $0.56774$ | $0.56774$ | $0.49501$ | $0.49501$ | 0 |
| 0.3 | 0 | $-0.20$ | $0.22824$ | $0.22824$ | $0.19900$ | $0.19975$ | $-0.11265$ |
| 0.3 | $-0.1$ | $-0.10$ | $0.11585$ | $0.11585$ | $0.10101$ | $0.10000$ | $-0.10820$ |
| $-0.05$ | 0.25 | $-0.10$ | $0.10229$ | $0.10229$ | $0.10191$ | $0.10000$ | $+0.07782$ |

読み方: 6 行すべてで**ODE は (5.43c) と 5 桁一致**する (1 次の外挿の精度)。$\psi = 0$ ではヘルプの式とも一致、斜め入射 ($\psi \neq 0$) の対称な面では $1/\cos^3\psi$ の分だけヘルプの式からずれ (3 行目)、非対称な面ではさらに $\cos$ の比の分がずれる。最後の列は、条件 (X) の $\rho_b$ で残る出口の位置のずれで、(D) の最後の式 $x_{\rm out} = [\cdots]\cdot\hat e_{\rm ref}$ の予測と同じ値 (5 桁)。ヘルプの $L'$ をそのまま使った $\rho_b$ では、$x_{\rm out}/s_I$ は 3〜6 行目で $-0.0036$、$-0.114$、$-0.109$、$+0.0777$、$p_{xr}/s_I$ は $-0.0074$、$-0.0028$、$-0.0016$、$-0.0002$ ([bend_08_lprime_general.out.txt](outputs/bend_08_lprime_general.out.txt) (c))。$s_I$ の大きさは、F1 = 0.1 m、$\rho_b = 10$ m の tanh の端で $4.6 \times 10^{-5}$ m なので、残りのずれは $0.1\,s_I \approx 5 \times 10^{-6}$ m の桁。F1 の有限性 ($F1/\rho_b \sim 0.006$) による 2 次の補正は、$(L' - L)/s_I$ の相対 $7 \times 10^{-4}$ 以下 (F1 = 0.06 の生の値と外挿値の差。$\psi = 0$ では 5 桁で変わらない)。

**この導出でまだ確かめられていないこと** (具体的に)。(1) **$K0 \neq 0$**: 場の角 $\phi_b = \phi_0 + K0$ が基準軌道の角 $\phi_0$ と違うと、$x = p_x = 0$ の軌道は出口で向きが $K0$ だけ曲がる (閉じない) ので、上の「出口の向き $\Theta = \phi_0$」の条件を置く根拠が無い。K64 の式 (`tffs.f:1282`) は $\phi_0$ だけを使い K0 を含まない。条件を何にするかの定義 (ヘルプにも書かれていない) が足りず、ODE でも確認していない。(2) **$\mathrm{F1}^2$ の 2 次以上**: 解析的には導出していない。ODE の外挿は、上のとおり $F1/\rho_b \sim 0.006$ で 2 次の補正が $10^{-3}$ 以下であることまで。(3) **端の場の模型**: ODE の場は、面に平行な等高線の水平面内の場 ($y = 0$)。実際の端の磁場の曲率と $y$ への依存 (Maxwell) は入れていない。(4) **K64 の設計軌道**: 上のとおり非対称・斜め入射では物理的な端の磁場で出口の位置が $0.1\,s_I$ 程度ずれるのに対し、K64 のコードは設計軌道を (x, p_x) = 0 に保つ (`dxfr1*dp/p` の $\delta$ 依存部分だけを入れ、$L'$ で閉じる) ので、K64 の設計軌道は物理的な端の磁場の設計軌道そのものではない (gs・実機の磁場との比較は未実施)。K0 = 0 で F1 ≠ 0 のとき設計軌道 (0, 0, 0, 0, 0, 0) が出口で 0 のまま保たれることは [bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (c) で確認した ($\rho_b = L'/\mathrm{ANGLE}$ にとるので、これは $L'$ の補正の物理的な正しさの確認にはならない。物理的な確認は上の ODE)。

**F1 領域の分割 (放射のある追跡 `tbendr` と CALC の `tbendef1`)**。放射を入れる追跡と CALC は、本体の長さ $L'$ のうち入口と出口の $\mathrm{fb}/2$ ずつ (`f1r`・`f2r`) を、2 つの短い円弧に分ける (`tbendal`、`tbend.f:17-39`): 入口は先に $\mathrm{rbl}\cdot f1r$、次に $\mathrm{rbh}\cdot f1r$ の長さ ($n = -1, 0$)、出口は $\mathrm{rbh}\cdot f2r$、$\mathrm{rbl}\cdot f2r$ ($n = \mathrm{ndiv} + 1, +2$)。$\mathrm{rbh} = \tfrac12 + 1/\sqrt{12}$、$\mathrm{rbl} = 1/(6\,\mathrm{rbh})$。

- **値**: $1/\sqrt{12} = 1/(2\sqrt3)$ なので $\mathrm{rbh} = \tfrac12 + \tfrac1{2\sqrt3}$、$\mathrm{rbl} = 1/(6\,\mathrm{rbh}) = \tfrac12 - \tfrac1{2\sqrt3}$ (積 $\mathrm{rbl}\cdot\mathrm{rbh} = \tfrac14 - \tfrac1{12} = \tfrac16$ より)。これは $[0, 1]$ の 2 点ガウス・ルジャンドルの節 $\tfrac12 \pm 1/(2\sqrt3)$ に一致し、和は 1、積は $1/6$ ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (e)、[bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (j))。
- **何をしているか (読み取り)**: 各円弧の偏向角は `phi0n = alx/al*phi0` ($\mathrm{alx} = \mathrm{rb}\cdot f1r$ なので $\mathrm{rb}\,f1r/L'\cdot\phi_0$) で、放射の呼び出し `tradk` に渡す長さは `alr = f1r` (両方の円弧で同じ)。放射は偏向角の変化から磁場を読むので、放射の長さ $f1r$ の間の平均の磁場は $\mathrm{rb}\,B_0$ になる。つまり、立ち上がりの全長 $\mathrm{F1} = 2f1r$ の領域を、「長さ $f1r$ ずつの 2 段の階段、高さ $\mathrm{rbl}\,B_0$ と $\mathrm{rbh}\,B_0$」で置き換えている (幾何は、磁場が全強度の長さ $f1r$ の円弧)。
- **この階段が合わせているもの**: 直線の立ち上がり $g(u) = u$ ($0 \le u \le 1$。全長 $\mathrm{F1}$) の $\int g^k\,du = 1/(k+1)$ と、2 段の階段の $\tfrac12(\mathrm{rbl}^k + \mathrm{rbh}^k)$ は、$k = 0, 1, 2, 3$ で一致する (2 点ガウス・ルジャンドルの精度は 3 次まで)。$k = 4$ は一致しない (差 $-5.6\times10^{-3}$)。$k = 1$ は偏向角 (磁場の積分)、$k = 2$ は放射のエネルギー損失 ($B^2$ の積分)、$k = 3$ は $B^3$ の積分に当たる。したがって、立ち上がり領域の放射は、磁場の積分・放射パワーの積分 ($B^2$)・3 乗の積分が、直線の立ち上がりと同じになる (sympy で確認: (j))。
- **未確認**: 上の「何をしているか」は、コードの引数 (`phi0n`・`alr`) からの読みで、`tradk` の内部 (偏向角から磁場をどう読むか) と、設計者の意図は確認していない。

放射なしの追跡は、`tbendcore` を本体 $L'$ の 1 回で通す。CALC は、$\mathrm{ndiv} = 1$ でも、F1 があれば `tbendef1` (2 つの短い円弧) を使う (§9)。

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

**段 6〜8 の各変数の意味** (`tbendcore` は $q = p_x/p$ の単位。$\texttt{ph2} = p_h/p$ は水平面内の運動量の大きさ)。入口の面の座標系で、$q_{x1}$ (`px1`) と面の法線方向 (`pz2` $= \sqrt{\texttt{ph2}^2 - q_{x1}^2}$) を、運動量の向きの角 $\theta_1$ で $q_{x1} = \texttt{ph2}\sin\theta_1$、$\texttt{pz2} = \texttt{ph2}\cos\theta_1$ と書く。出口の面の座標系では、`px3` $= -q_{x2}$ ((5.22) の符号を反対にした量) を $\texttt{px3} = \texttt{ph2}\sin\theta_3$、$\texttt{pz3} = \texttt{ph2}\cos\theta_3$ と書く ($\theta_3 = -\theta_2$)。設計軌道は $\theta_1 = \psi_1$、$\theta_3 = \psi_2$。以下の恒等式は、すべて 2 組の入力 ($p = 1.03$、$\psi = 0.03, 0.07$ など) で 40 桁まで確認した ([bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (e))。

| 変数 | 定義 (Fortran) | 意味 |
|---|---|---|
| `dpx2` | $q_{x1} - \texttt{ph2}\sin\psi_1$ | 面の座標系の運動量の、設計の向き ($\theta_1 = \psi_1$) からのずれ。$\texttt{ph2}(\sin\theta_1 - \sin\psi_1)$ |
| `t2` | $(q_{x1} + \texttt{ph2}\sin\psi_1)/(\texttt{pz2} + \texttt{ph2}\cos\psi_1)$ | $\tan\dfrac{\theta_1 + \psi_1}2$ (和積公式 $\dfrac{\sin\theta + \sin\psi}{\cos\theta + \cos\psi} = \tan\dfrac{\theta + \psi}2$)。これを使うと $\texttt{pz2} = \texttt{ph2}\cos\psi_1 - \texttt{dpx2}\cdot\texttt{t2}$ (桁落ちなく $\texttt{pz2} - \texttt{ph2}\cos\psi_1 = -\dfrac{q_{x1}^2 - \texttt{ph2}^2\sin^2\psi_1}{\texttt{pz2} + \texttt{ph2}\cos\psi_1}$ から) |
| `dpx3a` | $(x_2\sin\omega - \texttt{drho}\cdot\texttt{sinp2wp1})/\rho_e$、$\texttt{sinp2wp1} = \sin\psi_2 + \sin(\omega + \psi_1)$ | (5.22) の、位置 $x_2$ の項と、半径の差の項 ($\rho_0/\rho_e\cdot\texttt{sinp2wp1}$ を $\rho_e\,\texttt{ph2} - \rho_0 = \texttt{drho}$ で書き直したもの)。$\texttt{drho} = \rho_e\,\texttt{ph2} - \rho_0$ |
| `dpx3` | $\texttt{dpx3a} - \texttt{dpx2}(\cos\omega - \sin\omega\cdot\texttt{t2})$ | `px3` $- \texttt{ph2}\sin\psi_2$。導出: (5.22) の $-q_{x2}$ に $q_{x1} = \texttt{ph2}\sin\psi_1 + \texttt{dpx2}$、$\texttt{pz2} = \texttt{ph2}\cos\psi_1 - \texttt{dpx2}\,\texttt{t2}$ を入れると、$-\texttt{pz2}\sin\omega - q_{x1}\cos\omega = -\texttt{ph2}\sin(\omega + \psi_1) - \texttt{dpx2}(\cos\omega - \sin\omega\,\texttt{t2})$。残りの $\rho_0/\rho_e\cdot\texttt{sinp2wp1}$ と合わせて `px3` $= \texttt{ph2}\sin\psi_2 + \texttt{dpx3a} - \texttt{dpx2}(\ldots)$ |
| `px1px3` | $\texttt{ph2}(\sin\psi_2 + \sin\psi_1) + \texttt{dpx3a} + \texttt{dpx2}\cdot\texttt{sqwht2}$ | $q_{x1} + \texttt{px3}$ (和。恒等式で確認)。`sqwht2` $= (1 - \cos\omega) + \sin\omega\,\texttt{t2}$ |
| `dpz32` | $\texttt{px1px3}(q_{x1} - \texttt{px3})/(\texttt{pz2} + \texttt{pz3})$ | $\texttt{pz3} - \texttt{pz2}$ ($\texttt{pz}^2 = \texttt{ph2}^2 - \texttt{px}^2$ の差を和で割った、桁落ちのない形) |
| `t3` | $(\texttt{px3} + \texttt{ph2}\sin\psi_2)/(\texttt{pz3} + \texttt{ph2}\cos\psi_2)$ | $\tan\dfrac{\theta_3 + \psi_2}2$ (`t2` の出口版) |
| `t2t3` | 長い式 (`tbend.f` の段 7) | **積ではなく和 `t2` + `t3`**。分母の差 $D_3 - D_2 = \texttt{dpz32} - \texttt{ph2}\cdot\texttt{dcosp}$ ($D_2 = \texttt{pz2} + \texttt{ph2}\cos\psi_1$、$D_3 = \texttt{pz3} + \texttt{ph2}\cos\psi_2$) を使って通分すると余りが消える (§7.1 (1))。桁落ちを避けるための書き方 |
| `pz2px1` | $\texttt{pz2}\cos\psi_1 + q_{x1}\sin\psi_1$ | 入口の面の座標系の運動量 $(q_{x1}, \texttt{pz2})$ を基準面の座標系に戻した $z$ 成分 $= \texttt{ph2}\cos(\theta_1 - \psi_1)$ |
| `t4` | $(\cos\psi_2 + \texttt{t3}\sin\psi_2)\cdot\texttt{pz2px1}$ | $\texttt{ph2}\cos(\theta_1 - \psi_1)\,(\cos\psi_2 + \texttt{t3}\sin\psi_2)$ (定義から) |
| `x3` | `tbend.f` の段 7 | 出口の面の上の位置 (5.23)。構造: $\rho_e(\texttt{pz3} - \texttt{pz2}\cos\omega + q_{x1}\sin\omega) = \rho_e(\texttt{dpz32} + \texttt{sqwh}\,\texttt{pz2} + \sin\omega\,q_{x1})$ ($\texttt{pz3} - \texttt{pz2} + \texttt{pz2}(1 - \cos\omega)$ に分けた)。Fortran の式と (5.23) の差は $(\rho_0/\texttt{ph2})\,R$ で、$R$ は加法定理で 0 になる (§7.1 (2)。数値は $10^{-40}$ で一致) |
| `da` | $\sin(\texttt{da}) = [\ldots]/\texttt{ph2}^2$ を `asinz` したもの | $\texttt{da} = \Theta - \phi_0 = (\theta_1 - \psi_1) + (\theta_3 - \psi_2)$ ($\Theta$ は (5.12) の回転角。$\Theta = \theta_1 + \theta_3 + \omega$、$\phi_0 = \omega + \psi_1 + \psi_2$)。`asinz` の引数は $\sin(\theta_1 - \psi_1 + \theta_3 - \psi_2)$ に等しい (§7.1 (3)。$10^{-42}$)。$\theta_1 - \psi_1$、$\theta_3 - \psi_2$ が小さいので、角そのものでなく $\sin$ の和の形で書き、$\arcsin$ を最後に 1 回だけ取る |

(各変数の導出は、次の §7.1 で 1 段ずつ追う。全体が (5.12)〜(5.24) の幾何で $10^{-15}$ まで再現されることは §5.3・[bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (d) で確認してある。)

放射を入れる場合は、`tbendr` が長さを `ndivrad` ($\mathrm{ndivmax} = 1024$) で刻み、各刻みの `tbendcore` のあとに `tradk` (光子の放射。この章の範囲外) を入れる。

### 7.1 `t2t3`・`x3`・`da` の導出 (紙の上で 1 段ずつ)

`tbendcore` の段 7 の長い式を、中間変数を置いて 1 段ずつ (5.23) と $\Theta$ に結ぶ。以下、$q$ の単位で、$h \equiv \texttt{ph2}$ と書く。面の座標系の運動量の向きの角 $\theta_1$ ($q_{x1} = h\sin\theta_1$、$\texttt{pz2} = h\cos\theta_1$) と $\theta_3$ ($\texttt{px3} = h\sin\theta_3$、$\texttt{pz3} = h\cos\theta_3$) を使い、次の略記を置く。

$$
A \equiv \frac{\theta_1 + \psi_1}2,\ \ B \equiv \frac{\theta_3 + \psi_2}2,\ \ a_1 \equiv \frac{\theta_1 - \psi_1}2,\ \ a_3 \equiv \frac{\theta_3 - \psi_2}2,\ \ \Delta_1 \equiv \theta_1 - \psi_1 = 2a_1,\ \ \Delta_3 \equiv \theta_3 - \psi_2 = 2a_3
$$

使う道具は、和積公式 $\sin\theta - \sin\psi = 2\cos\tfrac{\theta + \psi}2\sin\tfrac{\theta - \psi}2$、$\cos\theta - \cos\psi = -2\sin\tfrac{\theta + \psi}2\sin\tfrac{\theta - \psi}2$ と、加法定理だけ。基本の関係を先に 4 つ。

- **(F1)** $\texttt{t2} = \tan A$、$\texttt{t3} = \tan B$ (§7 の表の和積公式。$\texttt{t2} = (q_{x1} + h\sin\psi_1)/(\texttt{pz2} + h\cos\psi_1) = (\sin\theta_1 + \sin\psi_1)/(\cos\theta_1 + \cos\psi_1)$)。
- **(F2)** $\texttt{dpx2} = q_{x1} - h\sin\psi_1 = h(\sin\theta_1 - \sin\psi_1) = 2h\cos A\sin a_1$。
- **(F3)** `dpx3a` の角の形。§7 の表の `dpx3` の行から $\texttt{dpx3a} = (\texttt{px3} - h\sin\psi_2) + \texttt{dpx2}(\cos\omega - \sin\omega\,\texttt{t2})$。第 1 項は $h(\sin\theta_3 - \sin\psi_2) = 2h\cos B\sin a_3$。第 2 項は $2h\cos A\sin a_1\cdot\dfrac{\cos\omega\cos A - \sin\omega\sin A}{\cos A} = 2h\sin a_1\cos(\omega + A)$。よって $\texttt{dpx3a} = 2h\bigl[\cos B\sin a_3 + \sin a_1\cos(\omega + A)\bigr]$。
- **(F4)** $\Theta = \theta_1 + \theta_3 + \omega$、$\phi_0 = \omega + \psi_1 + \psi_2$ (§5.3) から $A + B + \omega = (\Theta + \phi_0)/2$、$\Delta_1 + \Delta_3 = \Theta - \phi_0$。

**(1) `t2t3` は和 `t2` + `t3`**。Fortran の式 $\texttt{t2t3} = \dfrac{h\sin\psi_2 + \texttt{px1px3}}{D_3} + \dfrac{h\sin\psi_1 + q_{x1}(\texttt{dpz32} - h\,\texttt{dcosp})/D_3}{D_2}$ ($D_2 = \texttt{pz2} + h\cos\psi_1$、$D_3 = \texttt{pz3} + h\cos\psi_2$) で、`px1px3` $= q_{x1} + \texttt{px3}$ ($\texttt{px1px3}$ の表の式を $\theta$ で書いて確認)。
1. 第 1 項: $(h\sin\psi_2 + \texttt{px3} + q_{x1})/D_3 = \texttt{t3} + q_{x1}/D_3$ ($\texttt{t3} = (\texttt{px3} + h\sin\psi_2)/D_3$)。
2. 第 2 項の分子の $q_{x1}$ の項を分けて、$\dfrac{h\sin\psi_1}{D_2} + \dfrac{q_{x1}(\texttt{dpz32} - h\,\texttt{dcosp})}{D_2D_3}$。
3. 差 $D_3 - D_2 = \texttt{pz3} - \texttt{pz2} + h(\cos\psi_2 - \cos\psi_1) = \texttt{dpz32} - h\,\texttt{dcosp}$ ($\texttt{dcosp} = \cos\psi_1 - \cos\psi_2$)。これを 2 の第 2 項に入れると $\dfrac{q_{x1}(D_3 - D_2)}{D_2D_3} = q_{x1}\left(\dfrac1{D_2} - \dfrac1{D_3}\right)$。
4. 1〜3 を足すと $\texttt{t3} + \dfrac{q_{x1}}{D_3} + \dfrac{h\sin\psi_1}{D_2} + \dfrac{q_{x1}}{D_2} - \dfrac{q_{x1}}{D_3} = \texttt{t3} + \dfrac{q_{x1} + h\sin\psi_1}{D_2} = \texttt{t3} + \texttt{t2}$。

**(2) `x3` が (5.23) に一致する**。(5.23) を $q$ の単位で書く ($\rho_b p_{z} \to \rho_e\,\texttt{pz}$): $x_3 = x_2\cos\omega + \rho_eB_q + \rho_0K$、$B_q \equiv \texttt{pz3} - \texttt{pz2}\cos\omega + q_{x1}\sin\omega = \texttt{dpz32} + \texttt{sqwh}\,\texttt{pz2} + \sin\omega\,q_{x1}$、$K \equiv \cos(\omega + \psi_1) - \cos\psi_2$。
1. $\texttt{drho} = \rho_e h - \rho_0$ ($\texttt{drho} = \rho_e\,\texttt{ph2} - \rho_0$) から $\rho_e = (\texttt{drho} + \rho_0)/h$。よって $\rho_eB_q = \texttt{drho}\,B_q/h + \rho_0B_q/h$。
2. Fortran の `x3` の第 1 項 $x_2\bigl(\cos\omega - \tfrac{\rho_0}{\rho_e}\tfrac{\texttt{t3}\sin\omega}{h}\bigr)$ の $x_2\sin\omega$ を、`dpx3a` の定義 $\rho_e\,\texttt{dpx3a} = x_2\sin\omega - \texttt{drho}\,\texttt{sinp2wp1}$ で置き換える: $\dfrac{\rho_0\texttt{t3}}{\rho_eh}x_2\sin\omega = \dfrac{\rho_0\texttt{t3}}{h}\left(\texttt{dpx3a} + \dfrac{\texttt{drho}\,\texttt{sinp2wp1}}{\rho_e}\right)$。
3. Fortran の最後の項 $-\texttt{drho}\bigl(-\texttt{sinp2wp1}\tfrac{\rho_0}{\rho_e}\texttt{t3} - \texttt{dpz32} - \texttt{sqwh}\,\texttt{pz2} - \sin\omega\,q_{x1}\bigr)/h = \texttt{drho}\bigl(B_q + \texttt{sinp2wp1}\tfrac{\rho_0}{\rho_e}\texttt{t3}\bigr)/h$。2 の $\texttt{drho}\,\texttt{sinp2wp1}$ の項と**打ち消し合う** ($-\tfrac{\rho_0\texttt{t3}}{h}\tfrac{\texttt{drho}\,\texttt{sinp2wp1}}{\rho_e} + \tfrac{\texttt{drho}}{h}\texttt{sinp2wp1}\tfrac{\rho_0}{\rho_e}\texttt{t3} = 0$)。
4. 残ると Fortran の `x3` $= x_2\cos\omega + \texttt{drho}\,B_q/h + \dfrac{\rho_0}{h}\bigl[(\cos\omega\,(\texttt{t2} + \texttt{t3}) + \sin\omega\,(1 - \texttt{t2}\cdot\texttt{t3}))\,\texttt{dpx2} - \texttt{t3}\,\texttt{dpx3a}\bigr]$ (Fortran の変数 `t2t3` は (1) で和 $\texttt{t2} + \texttt{t3}$ と分かったので置き換え、$\texttt{t2}\cdot\texttt{t3}$ は積)。(5.23) との差は、1 の式から $(\rho_0/h)\,R$、
   $R \equiv \bigl[\cos\omega\,(\tan A + \tan B) + \sin\omega\,(1 - \tan A\tan B)\bigr]\texttt{dpx2} - \tan B\,\texttt{dpx3a} - B_q - hK$。
5. **$R = 0$**: ($\tan A + \tan B = \sin(A + B)/(\cos A\cos B)$、$1 - \tan A\tan B = \cos(A + B)/(\cos A\cos B)$ の通分で) 括弧 $= \dfrac{\cos\omega\sin(A + B) + \sin\omega\cos(A + B)}{\cos A\cos B} = \dfrac{\sin(A + B + \omega)}{\cos A\cos B}$。(F2) の $\texttt{dpx2} = 2h\cos A\sin a_1$ をかけて $2h\sin a_1\sin(A + B + \omega)/\cos B$。(F3) から $\tan B\,\texttt{dpx3a} = 2h\bigl[\sin B\sin a_3 + \sin a_1\sin B\cos(\omega + A)/\cos B\bigr]$。引くと、$\sin(A + B + \omega) - \sin B\cos(\omega + A) = \cos B\sin(\omega + A)$ ($\sin(X + B) = \sin X\cos B + \cos X\sin B$、$X = \omega + A$) から
   $\bigl[\cdots\bigr]\texttt{dpx2} - \tan B\,\texttt{dpx3a} = 2h\bigl[\sin a_1\sin(\omega + A) - \sin B\sin a_3\bigr]$。
   一方 $B_q + hK = h\bigl[\cos\theta_3 - \cos(\theta_1 + \omega) + \cos(\omega + \psi_1) - \cos\psi_2\bigr]$ ($B_q = h[\cos\theta_3 - \cos\theta_1\cos\omega + \sin\theta_1\sin\omega] = h[\cos\theta_3 - \cos(\theta_1 + \omega)]$)。和積公式で、$\cos(\omega + \psi_1) - \cos(\omega + \theta_1) = 2\sin(\omega + A)\sin a_1$ ($\tfrac{(\omega + \psi_1) + (\omega + \theta_1)}2 = \omega + A$、$\tfrac{(\omega + \psi_1) - (\omega + \theta_1)}2 = -a_1$)、$\cos\theta_3 - \cos\psi_2 = -2\sin B\sin a_3$。よって $B_q + hK = 2h[\sin(\omega + A)\sin a_1 - \sin B\sin a_3]$ で、$R = 0$。

**(3) `da` の `asinz` の引数が $\sin(\Delta_1 + \Delta_3)$**。引数は $N/h^2$、
$N = \texttt{dpx2}\bigl[\mathrm{Bk}_0 + \texttt{t4}\,\texttt{sqwht2}\bigr] + \texttt{dpx3a}\,\texttt{t4}$、
$\mathrm{Bk}_0 \equiv \sin\psi_1\bigl(\texttt{t2}(\texttt{pz3}\cos\psi_2 + \texttt{px3}\sin\psi_2) - q_{x1}\cos\psi_2\bigr) - \sin\psi_2\bigl(\texttt{t3}\,\texttt{pz2px1} - \texttt{px3}\cos\psi_1\bigr) + \cos\psi_1\cos\psi_2\,\texttt{dpz32}$。
1. **各部分を角で書く**: $\texttt{pz3}\cos\psi_2 + \texttt{px3}\sin\psi_2 = h\cos\Delta_3$、$\texttt{pz2px1} = \texttt{pz2}\cos\psi_1 + q_{x1}\sin\psi_1 = h\cos\Delta_1$ (加法定理)。$\texttt{t4} = (\cos\psi_2 + \tan B\sin\psi_2)\,h\cos\Delta_1 = h\cos\Delta_1\cos(B - \psi_2)/\cos B = h\cos\Delta_1\cos a_3/\cos B$ ($B - \psi_2 = a_3$)。$\texttt{sqwht2} = (1 - \cos\omega) + \sin\omega\tan A = (\cos A - \cos(\omega + A))/\cos A$。
2. **$\omega$ を含む項が消える**: $\texttt{dpx2}\,\texttt{t4}\,\texttt{sqwht2} = 2h\sin a_1\,\texttt{t4}\,(\cos A - \cos(\omega + A))$ と、(F3) の $\omega$ を含む項 $2h\sin a_1\cos(\omega + A)\,\texttt{t4}$ の和は $2h\sin a_1\cos A\,\texttt{t4}$ ($\cos(\omega + A)$ が打ち消し合う)。
3. (F3) の残りの項: $2h\cos B\sin a_3\cdot\texttt{t4} = h^2\cos\Delta_1\cdot2\sin a_3\cos a_3 = h^2\cos\Delta_1\sin\Delta_3$。
4. 2・3 と (F2) から $N = 2h\cos A\sin a_1\bigl(\mathrm{Bk}_0 + \texttt{t4}\bigr) + h^2\cos\Delta_1\sin\Delta_3$。これが $h^2\sin(\Delta_1 + \Delta_3) = h^2(\sin\Delta_1\cos\Delta_3 + \cos\Delta_1\sin\Delta_3)$ に等しいためには、$\sin\Delta_1 = 2\sin a_1\cos a_1$ から、(I2) $\ \mathrm{Bk}_0/h + \cos\Delta_1\cos a_3/\cos B = \cos a_1\cos\Delta_3/\cos A$ が要る。
5. **(I2) の証明**: $\mathrm{Bk}_0/h = \sin\psi_1\tan A\cos\Delta_3 - \sin\psi_1\sin\theta_1\cos\psi_2 - \sin\psi_2\tan B\cos\Delta_1 + \sin\psi_2\sin\theta_3\cos\psi_1 + \cos\psi_1\cos\psi_2(\cos\theta_3 - \cos\theta_1)$ (1 の置き換えと $q_{x1} = h\sin\theta_1$ 等)。(I2) の左辺 $-$ 右辺を 3 つのまとまりに分ける。
   - $\cos\Delta_3\bigl(\sin\psi_1\tan A - \cos a_1/\cos A\bigr) = \cos\Delta_3\,(\sin\psi_1\sin A - \cos a_1)/\cos A$。$a_1 = A - \psi_1$ なので $\cos a_1 = \cos A\cos\psi_1 + \sin A\sin\psi_1$、よって括弧 $= -\cos A\cos\psi_1$ で、このまとまりは $-\cos\psi_1\cos\Delta_3$。
   - $\cos\Delta_1\bigl(\cos a_3/\cos B - \sin\psi_2\tan B\bigr) = \cos\Delta_1(\cos a_3 - \sin\psi_2\sin B)/\cos B$。$a_3 = B - \psi_2$ なので $\cos a_3 = \cos B\cos\psi_2 + \sin B\sin\psi_2$、よって $= \cos\psi_2\cos\Delta_1$。
   - 残り: $-\sin\psi_1\sin\theta_1\cos\psi_2 + \sin\psi_2\sin\theta_3\cos\psi_1 + \cos\psi_1\cos\psi_2(\cos\theta_3 - \cos\theta_1)$。
   3 つを足し、$\cos\Delta_3 = \cos\theta_3\cos\psi_2 + \sin\theta_3\sin\psi_2$、$\cos\Delta_1 = \cos\theta_1\cos\psi_1 + \sin\theta_1\sin\psi_1$ を入れると、8 つの項が 2 つずつ打ち消し合って 0 になる ($-\cos\psi_1\cos\psi_2\cos\theta_3$ と $+\cos\psi_1\cos\psi_2\cos\theta_3$、$-\cos\psi_1\sin\psi_2\sin\theta_3$ と $+\sin\psi_2\sin\theta_3\cos\psi_1$、$+\cos\psi_2\cos\psi_1\cos\theta_1$ と $-\cos\psi_1\cos\psi_2\cos\theta_1$、$+\cos\psi_2\sin\psi_1\sin\theta_1$ と $-\sin\psi_1\sin\theta_1\cos\psi_2$)。
6. 以上から $N/h^2 = \sin(\Delta_1 + \Delta_3)$ で、`da` $= \arcsin(\cdots) = \Delta_1 + \Delta_3 = \Theta - \phi_0$ ((F4))。

各中間式 ((F1)〜(F3)、(1) の差、$R = 0$ の各行、2 の $\omega$ の消去、(I2)、$\texttt{pz2px1}$・$\texttt{t4}$・$\texttt{sqwht2}$ の角の形) を、3 組の入力 ($\psi_1 = 0.03, 0.3, -0.05$ など) で 40 桁で確認した (最大のずれ $1.4 \times 10^{-40}$): [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (c)。


## 8. K1 ≠ 0 (複合機能の偏向磁石)

K1 ≠ 0 の追跡は `tbendi` (`tbendi.f:1-200`) と `tbend.f` の `tbendiinit`・`tbendicorr`・`tbendibody`。本体を $n$ 個のスライス (長さ `alx`) に分け、各スライスで次の Strang 型の合成 ($\mathcal{C}$ は残りの非線形項、$\mathcal{B}$ は $H_2$ の厳密解 (5.8)):

$$
\mathcal{C}\!\left(\tfrac{\ell}{2}\right)\ \mathcal{B}(\ell)\ \mathcal{C}\!\left(\tfrac{\ell}{2}\right)
\tag{5.50}
$$

(隣り合う $\mathcal{C}$ は 1 つにまとめる: `tbendicorr` は `(akxn(n)+akx0)*.5`、`(alx(n)+alx0)*.5` で呼ぶ。`tbendi.f:115-121`、最後の半分は `:147`。) $\mathcal{B}$ が `tbendibody` で、(5.8)(5.11)。$\mathcal{C}$ が `tbendicorr` (`tbend.f:118-134`) で、$H_{\rm rest}$ (5.3) に当たる運動学的な非線形ドリフトと、K1 の非線形な蹴り (`pxi -= ak1*xi*xr*(0.5 - xr*(2-xr)/12)/p`、$x_r = x/\rho_0$) から成る。

- **端面**: K1 ≠ 0 の追跡は、端面を薄い蹴り (`pxi = px + tanp1*xi/rhoe`、`pyi = py - (tanp1 + fpx)*f`) で近似し、面までの直進と面から扇形の面までの弧 ($x\tan\psi$ の経路長の変化) を入れない。CALC (`tbende`) は K1 = 0 と同じ厳密な鎖 (`tbedge`) を使う。この違いは docs/30 の B117 (v2 は `DISFRIN = 0` で厳密な `tbedge` にそろえる)。DISFRIN = 0 のとき `ttfrin` (K1 の四極の非線形フリンジ) が入る (`mfring > -4`・`mfring /= 2`)。
- **分割数**: `tbendi.f:41` は $n = 1 + \mathrm{int}(|\phi_0/\varepsilon|)$、$\varepsilon = 10^{-3}\cdot\mathrm{EPS}$ (`epsbend = 1.d-3`、`tbend.f:12`)。放射があるときは `ndivrad` と大きい方 (上限 1000)。CALC (`tbende.f:455-460`) も $n = 1 + \mathrm{int}(|\phi_c/\varepsilon|)$ ($\phi_c = \phi_0\cdot r_{bc}$、F1 の分を除いた角)。**ヘルプの $1 + \lfloor\sqrt{|K1L'|/(12\cdot10^{-5}\mathrm{EPS})}\rfloor$ とは別の式** (§13・付録 A)。
- **K1 の定義**: K64 の K1 は、中央面の磁場 $B_y = B_0 + (K1/h)\ln(r/\rho_0)$ ($r = \rho_0 + x$) の勾配、すなわち曲線座標の勾配 $(1 + hx)\,\partial b/\partial x$ が一様、という暗黙の定義 (誤りではなく定義。docs/29 の A59、docs/30 の B160)。

  **$\mathcal{C}$ の K1 の項の導出** (確認済み。sympy: [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) (k)、[bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (d))。中央面 ($y = 0$) の磁場 $b_y(x)$ を、上の K1 の定義から 1 段ずつ作る ($K = K1/\mathrm{al}$ [1/m²]、$h = 1/\rho_0$)。
  1. **$b_y(x)$**: 定義 $(1 + hx)\,\partial b_y/\partial x = K$ (勾配が一様) から $\partial b_y/\partial x = K/(1 + hx)$。$x = 0$ で $b_y = h_b$ として $x$ で積分: $b_y(x) = h_b + K\int_0^x\dfrac{dx'}{1 + hx'} = h_b + \dfrac Kh\ln(1 + hx)$。小さい $x$ では $\dfrac Kh\cdot hx = Kx$。
  2. **力**: (5.2) の $a_s$ の $x$ 微分は 01a の (1a.60a) から $\partial a_s/\partial x = -(1 + hx)\,b_y$ (一様磁場 $b_y = h_b$ なら (5.1) の $-(1 + hx)h_b$)。$p_x' = -\partial H/\partial x$ で、$H = -a_s - (1 + hx)\sqrt{\cdots}$ の $a_s$ の項からの寄与は $\partial a_s/\partial x = -(1 + hx)\,b_y(x)$、根号の項からの寄与は $hP_z$ ($P_z = \sqrt{p^2 - p_x^2 - p_y^2}$)。$b_y = h_b + \Delta b$ ($\Delta b \equiv b_y - h_b = (K/h)\ln(1 + hx)$) とおくと $p_x' = -(1 + hx)h_b + hP_z - (1 + hx)\,\Delta b$。最初の 2 項は偏向の分 ((5.4)〜$H_2$ と $H_{\rm rest}$ に入っている)、**最後の項 $-(1 + hx)\Delta b$ が K1 の力**。
  3. **級数**: $u \equiv hx = x_r$ とおくと、K の力は $-(1 + u)\,\dfrac Kh\ln(1 + u)$。$x$ の 1 次は $-\dfrac Kh\cdot u = -Kx$ で、$H_2$ の $K(x^2 - y^2)/2$ の力 ((5.4) の $-Kx$) に入っている。級数展開すると

  $$
  (1 + u)\ln(1 + u) = u + \frac{u^2}2 - \frac{u^3}6 + \frac{u^4}{12} - \frac{u^5}{20} + \cdots
  $$

  ($(1 + u)(u - u^2/2 + u^3/3 - u^4/4 + u^5/5)$ の各次を足す: $u^2$: $-\tfrac12 + 1$、$u^3$: $\tfrac13 - \tfrac12$、$u^4$: $-\tfrac14 + \tfrac13$、$u^5$: $\tfrac15 - \tfrac14$)。最初の項 $(K/h)u = Kx$ が $H_2$ の線形な力で、残り ($u^2$ 以上) が $\mathcal{C}$ の非線形な蹴り: $\Delta p_x = -\dfrac{K\ell}{h}\left(\dfrac{u^2}{2} - \dfrac{u^3}6 + \dfrac{u^4}{12}\right) = -K\ell\,x\,x_r\left(\dfrac12 - \dfrac{x_r}6 + \dfrac{x_r^2}{12}\right)$ ($1/h = \rho_0$、$u^2/h = x\,x_r$)。K64 の `ak1*xi*xr*(0.5d0-xr*(2.d0-xr)/12.d0)` は $\tfrac12 - x_r(2 - x_r)/12 = \tfrac12 - \tfrac{x_r}6 + \tfrac{x_r^2}{12}$ で、**$u^4$ までの級数と一致** ($\mathrm{ak1}$ は $K\ell$。$q$ の単位への換算の $1/p$ が付く)。$u^5$ 以降 ($-u^5/20$) は捨てている (sympy: 差 $= -u^5/20$)。

- **`tbendicorr` のその他の行 (運動学の残り)**: $H_{\rm rest}$ ($H$ から $H_2$ を引いたもの) に対応する、厳密な運動方程式から $H_2$ の分を引いた量を、スライスの長さ $\ell$ だけ、運動量 ($\texttt{pxi}$、$\texttt{pyi}$) を初期値に固定して積分したもの。$s = \min(\texttt{psqmax}, q_x^2 + q_y^2)$、$\texttt{dpzi} = p_z/p - 1$ ($\texttt{sqrt1}(-s)$)、$\texttt{pzi} = 1 + \texttt{dpzi}$ とする。
  - `zi`: 厳密な $z' = (1 - dv) - (1 + hx)p/P_z$ (01 章 (1.8))。$H_2$ の分 (5.10) $z' = -dv - hx - \tfrac12q^2$ を引くと、$[1 - (1 + hx)/p_{zi}] - [-hx - q^2/2] = \tfrac12q^2 + (1 + hx)(1 - 1/p_{zi}) = \tfrac12s + (1 + x/\rho_0)\,\texttt{dpzi}/\texttt{pzi}$ ($1 - 1/p_{zi} = (p_{zi} - 1)/p_{zi}$)。これに $\ell$ をかけたものが `zi = zi + (.5*s + (1 + xi/rho0)*dpzi/pzi)*al`。
  - `pxi`: 厳密な $p_x' = hP_z - h_b(1 + hx)$ ((5.2) の $a_s$ と根号の項) のうち、$H_2$ の $p_x' \ni hp - h_b - hh_bx$ との差は $h(P_z - p)$。$q$ の単位では $q_x' = h\,\texttt{dpzi}$。$\ell$ をかけて $\texttt{pxi} = \texttt{pxi} + \texttt{phi}\cdot\texttt{dpzi}$ ($\texttt{phi} = h\ell$)。
  - `xi`: $x' = (1 + hx)\,q_x/p_{zi}$ から $H_2$ の $x' = q_x$ を引いた残り $= q_x(1/p_{zi} - 1) + hx\,q_x/p_{zi} = a + bx$、$a = -q_x\,\texttt{dpzi}/\texttt{pzi}$、$b = h\,q_x/\texttt{pzi}$ ($1/p_{zi} - 1 = (1 - p_{zi})/p_{zi} = -\texttt{dpzi}/\texttt{pzi}$)。$x$ の 1 次の線形方程式 ($q_x$ は固定) で、厳密解は $x\,e^{b\ell} + a(e^{b\ell} - 1)/b$。K64 は $(x + a\ell)/(1 - b\ell)$ ($= (-\ell\,\texttt{dpzi}\,q_x + \texttt{pzi}\,x)/(\texttt{pzi} - \texttt{phi}\,q_x)$。分母分子に $\texttt{pzi}$ をかけた形) を使う。これは $e^{b\ell} \to 1/(1 - b\ell)$ の置き換え (後退オイラー型)。**差の導出** ($\ell$ の 2 次まで展開): 厳密解 $= x(1 + b\ell + \tfrac12b^2\ell^2) + a\ell(1 + \tfrac12b\ell) + O(\ell^3)$ ($(e^{b\ell} - 1)/b = \ell(1 + b\ell/2 + \cdots)$)。K64 の式 $= (x + a\ell)(1 + b\ell + b^2\ell^2 + \cdots)$ ($1/(1 - b\ell)$ の級数) $= x(1 + b\ell + b^2\ell^2) + a\ell(1 + b\ell) + O(\ell^3)$。**厳密 $-$ K64** $= x\bigl(\tfrac12 - 1\bigr)b^2\ell^2 + a\ell\bigl(\tfrac12 - 1\bigr)b\ell = -\tfrac12b\ell^2(a + bx)$。$a + bx = (q_x/\texttt{pzi})(hx - \texttt{dpzi})$、$b = hq_x/\texttt{pzi}$ を入れると
    $$\text{厳密} - \text{K64} = \frac{\ell^2h\,q_x^2}{2\,\texttt{pzi}^2}\,(\texttt{dpzi} - hx) + O(\ell^3)$$
    (符号は**厳密 $-$ K64**。sympy で展開と代入の差 0: [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) (d)。)
  - `yi`: $y' = (1 + hx)q_y/p_{zi}$ から $H_2$ の $y' = q_y$ を引いた残り $= q_y(-\texttt{dpzi} + hx)/\texttt{pzi}$。積分して $\texttt{yi} + (-\ell\,\texttt{dpzi} + x\,\texttt{phi})\,q_y/\texttt{pzi}$。ここで $x$ は**更新後**の `xi` (中点でなく端の値)。
  - まとめると、$\mathcal{C}$ の各半ステップは、運動量を初期値に固定した 1 次精度の近似 (誤差は $\ell^2$ の項)。スライスの数 $n$ を増やすと全体が $O(1/n)$ で厳密解に近づくはずだが、Strang 型の合成 (5.50) での全体の次数は、この章では評価していない (未確認)。`zi` は更新前の `xi`、`yi` は更新後の `xi` を使うので、$z$ と $y$ の端の値の取り方が違う (どちらも $\ell^2$ の誤差)。
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

**一致の程度** (電子 3 GeV、$L = 1$ m、ANGLE = 0.1 rad。追跡の中心差分ヤコビアンと、CALC の `lattice._el_cod_step` の 4×4 の差。[bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (a)): 設計軌道で $5 \times 10^{-13}$ (扇形・矩形・非対称とも)。小さいずれ ($\delta = 0.002$、0.3 mm、0.5 mrad) で $10^{-12}$ 〜 $9 \times 10^{-11}$。大きいずれ ($y = 10$ mm、$\delta = 0.05$) で、扇形 $2.2 \times 10^{-8}$、矩形 $3.0 \times 10^{-10}$、非対称 $5.7 \times 10^{-10}$。COD の更新 (出口の $x, p_x, y, p_y, z$) の差は小さいずれで $6 \times 10^{-11}$、大きいずれで $2.0 \times 10^{-8}$。F1 = 0.1 (FRINGE = 1、E1 = E2 = 0.5) では 4×4 が小さいずれで $10^{-7}$、大きいずれで $4.4 \times 10^{-5}$ (以上は k64)。**v2** ([bend_04_matrix_vs_map.v2.out.txt](outputs/bend_04_matrix_vs_map.v2.out.txt)。COD の更新が $p^2$ 型): 大きいずれの 4×4 の差は扇形 $7 \times 10^{-13}$、矩形 $1 \times 10^{-12}$、非対称 $1 \times 10^{-12}$、COD の更新の差は $1.3 \times 10^{-9}$ ($z$)、F1 = 0.1 (FRINGE = 1、E1 = E2 = 0.5) の 4×4 は小さいずれで $1.4 \times 10^{-8}$、大きいずれで $1.9 \times 10^{-7}$。

**差の出どころ**: (i) CALC の Maxwell フリンジの COD の更新が因子 $p$ を落としている (K64 の誤り。gs で確認、v2 は修正。付録 B-1)。(ii) F1 ≠ 0 のとき、E = 0 でも、CALC の写像 (COD の更新) の微分が追跡のヤコビアンと食い違う (小さいずれ $x$ = 1 mm、$y$ = 0.3 mm、$\delta$ = 0.002、F1 = 0.1 で、$x$ の行の $y$ の列が $4 \times 10^{-8}$、$p_y$ の列が $1 \times 10^{-7}$。F1 = 0 では $10^{-12}$)。CALC の解析行列は自身の COD 更新の数値微分と $10^{-8}$ で一致するので、写像そのものが違う ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (f))。CALC は F1 の位置が追跡と違い (面までの直進の前。E = 0 では直進が無いので同じ位置)、F1 の領域を常に 2 つの円弧 (`tbendef1`) に分けるが、これらが原因かどうかは確認していない。物理的にどちらが正しいかも確認していない。(iii) 段ごとの解析行列は、同じ段の COD 更新の数値微分と、`tbendebody0` で $1.4 \times 10^{-9}$、`tbedgedrift` で $4 \times 10^{-10}$、`tbedgebody` で $6 \times 10^{-12}$ だが、**`tbedgemaxwell` では $4.5 \times 10^{-5}$** ($y = 10$ mm、$\delta = 0.05$): 解析行列の $\partial x/\partial y$ は COD 更新の数値微分の $1.05 = 1 + \delta$ 倍 ([bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (c2))。解析行列の方が (5.36) の $p^2$ の式と整合し、COD の更新が不整合、という内部の矛盾 (K64 の誤り。v2 は COD の更新を直して解消。付録 B-1)。

**(i) の大きさの確かめ** ([bend_07_b1_estimate.out.txt](outputs/bend_07_b1_estimate.out.txt))。
1. 1 つの面での $p^2$ 型 (5.36) と K64 の $p$ 型 (`bf*pr`) の差は、$\Delta x$ で $\dfrac{b\,y^2}{\rho_b}\dfrac{p^2 - p}{p_v^3}$ ($\simeq \delta\,y^2/(2\rho_b)$ 程度)、$\Delta z$ で $-\dfrac{b\,y^2}{\rho_b}\dfrac{p_x(p - 1)}{p_v^3}$。$y = 10$ mm、$\delta = 0.05$、$\rho_b = 10$ m で、$\delta y^2/(2\rho_b) = 2.5 \times 10^{-7}$、実際の式では $q_x = 0$ で $\Delta x$ の差 $2.27 \times 10^{-7}$、$q_x = 0.1$ で $2.30 \times 10^{-7}$ (と $\Delta z$ の差 $-2.3 \times 10^{-8}$)。**この $2.5 \times 10^{-7}$ は 1 つの面の差の大きさ**。
2. 出口の面では符号が逆 ($f_{\rm out} = -f_{\rm in}$) なので、入口と出口の差がほぼ打ち消し合い、残るのは入口と出口の $y$ の違い ($y^2$ に比例) などによる分。扇形で $y_{\rm out} = 10.2$ mm ($p_y = 0.2$ mrad が 1 m 進む) なので $y^2$ が 4% 大きく、差し引き $|2.27 - 2.35| \times 10^{-7} \approx 8.7 \times 10^{-9}$ で、COD の $x$ の差 $9.7 \times 10^{-9}$ と同じ桁。$z$ の差 ($2.0 \times 10^{-8}$) は、入口の $\Delta x$ の差 $\varepsilon = 2.3 \times 10^{-7}$ が本体で $z' = -hx$ を通して効く分 ($\varepsilon L/\rho_0 \approx 2.3 \times 10^{-8}$) と同じ桁。面そのものの $\Delta z$ の差は $q_x$ に比例し、扇形 ($q_x \sim 0.005$) では $10^{-9}$ 台。(§1・§9・06 章の「$\delta y^2/(2\rho_b)$ 程度」は、**1 つの面での差の大きさ**で、要素全体の COD の差はそれより 1〜2 桁小さい ($2.5 \times 10^{-7}$ に対して $2.0 \times 10^{-8}$)、と読む。)
3. **差し替えの試験**: `lattice._tbedgemaxwell_cod` を (5.36) の $p^2$ 型 ($\Delta x$ に $p^2$、$\Delta z$ に $p_xp$) に差し替えて ($\texttt{pysad}$ のファイルは変えず、試験のためだけ)、`_bend_cod_update` (Python 版の COD の更新) と追跡の写像の差を取ると、$y = 10$ mm、$\delta = 0.05$ で $x$: $9.7 \times 10^{-9} \to 7 \times 10^{-17}$、$p_x$: $2.3 \times 10^{-9} \to 5 \times 10^{-18}$、$z$: $2.0 \times 10^{-8} \to -1.4 \times 10^{-9}$。扇形・矩形・非対称、$y = 5, 10$ mm、$\delta = 0.02, \pm0.05$ のすべてで、$x$・$p_x$ の差が $10^{-15}$ 以下になる。残る $z$ の差 ($1.35 \times 10^{-9}$) は、$y$ によらず $\delta$ にほぼ比例する ($\delta = 0.05, 0.02, -0.05$ で $-1.35$、$-0.56$、$+1.57 \times 10^{-9}$) 別の効果で、原因は特定していない。
4. したがって、**COD の更新の差 ($x$・$p_x$ と $z$ の大部分) は $p$ の因子で説明がつく**。ただし、これは「更新を $p^2$ 型にすると追跡と合う」ということで、K64 (gs) がどちらを意図しているか、gs の出力との比較は未実施。

**F1 ≠ 0 の 4×4 の差 ($4.4 \times 10^{-5}$。k64) は、v2 で $2 \times 10^{-7}$ 以下に下がる** ([bend_04_matrix_vs_map.v2.out.txt](outputs/bend_04_matrix_vs_map.v2.out.txt) (a)・(f)、[bend_07_b1_estimate.v2.out.txt](outputs/bend_07_b1_estimate.v2.out.txt) (c): F1 = 0.1・E = 0 で $1 \times 10^{-8}$ 前後、E1 = E2 = 0.5 で $0.8$〜$9 \times 10^{-7}$ (F1 = 0.05〜0.2)、$y = 10$ mm・$\delta = 0.05$ で $1.8 \times 10^{-7}$ (F1 = 0.1)。k64 では $\partial x_{\rm out}/\partial p_y$ が $2 \times 10^{-5}$〜$6 \times 10^{-5}$、$y$ にほぼ比例、$\delta$ にほぼ比例、F1・E1・E2 によらない (bend_07.out (c)))。v2 の変更のうち、この差を下げたのは (B324 の) CALC の Maxwell フリンジの COD の更新だと考えられるが、F1 ≠ 0 の経路 (`tblfre`・`tbendef1`) との関係の内訳は分離していない。残りの $10^{-8}$〜$10^{-6}$ の差の原因も未特定。物理的にどちらが正しいかは確認していない。

**シンプレクティック性**: 追跡のヤコビアン (数値微分) は $10^{-11}$ 台 (数値微分の精度)、$z$ の行も上の条件と $10^{-11}$ で一致。CALC の 6×6 は (F1 = 0 で) $5.5 \times 10^{-17}$ (扇形)・$1.1 \times 10^{-16}$ (矩形)・$4.4 \times 10^{-16}$ (非対称。E1 = 0.3、E2 = 0.7) ($z$ の行を条件で組むので。3 つのうち最大は $4.4 \times 10^{-16}$)。F1 つきでは $4.4 \times 10^{-5}$ (k64。v2 は $1.0 \times 10^{-8}$。[bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt)・[bend_04_matrix_vs_map.v2.out.txt](outputs/bend_04_matrix_vs_map.v2.out.txt) (d))。

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
| 追跡の K1 = 0 の全体 | `track6._tbend_ens` (`track6.py:1620`)、`_tbend0_ens` (`:1735`) | `tbend` | `el.E1`・`el.E2` は角 (rad)。`psi1 = el.E1` |
| 円弧と端面 (python) | `track6._tbendcore` (`:969`) | `tbendcore` | 変数名は Fortran と同じ。`bp = (rhob, rho0, drhob, fb1, fb2)` |
| 円弧と端面 (numba) | `track6n.tbendcore1` (`track6n.py:738`)、`k_bend0` (`:939`)、`tbendal` (`:926`) | 同上 | |
| 取り付け誤差 | `track6._tbshift` (`:1093`)、`_tbrot` (`:1135`) | `tbshift`・`tbrot` | |
| F1 領域の分割 | `track6._tbendal` (`:1433`) | `tbendal` | |
| K1 ≠ 0 | `track6._tbendiinit` (`:1191`)、`_tbendicorr` (`:1241`)、`_tbendibody` (`:1260`)、`_tbendi` (`:1279`)、`track6n.k_bendi` (`:1098`) | `tbendiinit` ほか | `xsxkx` の符号は `compat.v2()` で分かれる (B250) |
| CALC の COD の更新 | `lattice._bend_cod_update` (`lattice.py:621`)、`_bend_fringe_params` (`:596`)、`_el_cod_step` (`:2616`) | `tbende` | `fringe`・`mfring`・`fb1`・`fb2`。F1 = 0 は `ndiv = 1` |
| CALC の端面の COD の更新 | `lattice._tbedgedrift_cod` (`:456`)、`_tbedgebody_cod` (`:476`)、`_tbedgemaxwell_cod` (`:512`)、`_tbedge_cod` (`:533`) | `tbedgedrift` ほか | |
| CALC の円弧と F1 | `lattice._tbendebody0_cod` (`:552`)、`_tblfre_cod` (`:603`) | `tbendebody0`・`tblfre` | |
| CALC の解析行列 (numpy) | `lattice._tbedgemaxwell_mat` (`:692`)、`_tbedgedrift_mat` (`:718`)、`_tbedgebody_mat` (`:756`)、`_tbendebody0_mat` (`:809`)、`_bend_edge_chain` (`:875`)、`_bend_cod_matrix` (`:894`) | 各 `trans` | |
| CALC (numba) | `elkern._tbedgedrift` (`elkern.py:798`)、`_tbedgemaxwell` (`:847`)、`_tbedgebody` (`:882`)、`_tbendebody0` (`:958`)、`_bend_core` (`:1756`)、`k_bend` (`:1838`) | 同上 | |
| Emittance[] | `emit6d._tbende` (`emit6d.py:1268`)、`_tbende_core` (`:1321`)、`_tbendecorr` (`:1160`)、`_tbendebody` (`:1167`)、`perfkern_bend.bendecorr` (`:65`)、`bendebody` (`:138`) | `tbende` ほか | |
| 設計軌道の線形行列 | `elements.Bend.transfer_matrix` (`elements.py:234`) | | 別系統。K0 ≠ 0 と ANGLE ≠ 0 の扱いに注意 |
| 3 次のリー生成関数の取り出し | `bendnl.bend_decomposition` (`bendnl.py:140`) | | 追跡の写像から。この章の範囲外 |

## 12. 確認スクリプトと出力

すべて、リポジトリのルートから `py -3.12 physics_guide/scripts/bend_0x_*.py` で動く (数秒。`bend_06` だけ 1〜2 分)。出力は標準出力と `physics_guide/outputs/bend_0x_*.out.txt`。(a)〜(g) は各スクリプトの節の記号。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `bend_01_h2_body.py` | (a) $A_s$ から磁場、$H$ の 0〜3 次の展開 (sympy)。(b) 偏微分の表と運動方程式、K64 の `akxsq`・`drhop`・`akysq` との一致 (sympy で差 0)。(c) 水平の解が運動方程式を満たす、分散 $R_{16}$・$R_{26}$・$R_{56}$ (sympy)。(d) K64 の `tbendibody` の式と $H_2$ の ODE の比較 ($x, q, y, q_y$ は $10^{-18}$)、`xsxkx` の符号。(e) pysad の `_tbendiinit`・`_tbendibody` (k64 と v2)。(f) $H_2$ と厳密な $H$ の差の振幅依存 (比 4)。(g) 設計軌道の行列 | [bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) |
| `bend_02_arc_faces.py` | (a) 面までの直進のシンプレクティック性 (sympy、40 桁)。(b) 実験室系の幾何 (面・円・フリンジ) と追跡の 5 成分が $10^{-15}$ で一致 (扇形・矩形・非対称・AE・K0 ≠ 0、$\delta = 0, 0.02, -0.03$)。(c) 設計軌道。(d) ヘルプの本体の式 (5.22)〜(5.24) と幾何。(e) DISFRIN = 1 でも Maxwell フリンジが入る。(f) 曲線座標の $H$ の ODE と幾何 | [bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) |
| `bend_03_fringe.py` | (a) 端の磁場 (sympy)。(b) 硬い端の Maxwell フリンジ (5.35) とローレンツ力の ODE。(c) CALC `tbedgemaxwell` と追跡の式のシンプレクティック性。(d) F1 の $y$・$y^3$・$x$ の項 (ODE)、F1 の定義と $\int g'^2$。(e) `rbh`・`rbl`、DISFRIN | [bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) |
| `bend_04_matrix_vs_map.py` | (a) CALC の 4×4・δ の列と COD の更新を追跡と比較。(c)(c2) CALC の解析行列と COD 更新の数値微分、段ごとの差。(d) シンプレクティック性。(e) 設計軌道の 2×2 と端面の行列。(f) F1 ≠ 0 の CALC と追跡の差 | [bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) |
| `bend_06_derivations.py` | 途中式を足した箇所の確認 (pysad は使わない。1〜2 分): (a) $\mathrm{hi} = (H_2 + p)/p$、(b) $z$ の更新 (5.11) と数値積分 (30 桁)、(c) 面の幾何・内積・$\Theta$、(d) 設計軌道の各段の 2×2 行列と積、(e) `tbendcore` の `t2`・`t3`・`t2t3`・`x3`・`da` の意味 (40 桁)、(f) 端の磁場、(g) 硬い端の Maxwell フリンジ (C 付録の括弧)、(h) F1 の積分と (5.42)、(i) $L'$ の導出、(j) `rbh`・`rbl`、(k) `tbendicorr`・K1 の級数、(l) 分散 | [bend_06_derivations.out.txt](outputs/bend_06_derivations.out.txt) |
| `bend_07_b1_estimate.py` | 付録 B-1 の大きさの見積りと、`_tbedgemaxwell_cod` を $p^2$ 型に差し替える試験 (a)(b)、F1 ≠ 0 の 4×4 の差の依存性 (c) | [bend_07_b1_estimate.out.txt](outputs/bend_07_b1_estimate.out.txt) |
| `bend_08_lprime_general.py` | $L'$ の導出の拡張 (§6.4): (a) 弦の幾何 (sympy。$\psi_1 \neq \psi_2$ の $\cos$ の因子、横ずれ、ヘルプの式との差 $(\psi_1 - \psi_2)^2/8$)、(b) 斜め入射の端の $x$ のずれ $\Delta X = I_1p^2/(\rho_bp_z^3)$ と出口の符号 (mpmath の求積)、(c) 実験室系の ODE (DOP853) で入口・出口の端の場を通し、出口の向きが設計どおりになる $\rho_b$ から $L'$ を求めて (5.43c)・ヘルプの式と比較 | [bend_08_lprime_general.out.txt](outputs/bend_08_lprime_general.out.txt) |
| `bend_09_algebra_steps.py` | 紙の上で追う導出の各行の確認 (pysad は使わない): (a) 出口の面から基準面への直進 (5.21a)、(b) §5.4 の微分・通分・$E_xS_xE_x$ の 4 成分・$J_{\rm in}$ が対角になる理由、(c) §7.1 の `t2t3`・`x3`・`da` の中間式 (40 桁)、(d) §8 の $b_y$・`xi` の厳密解との差、(e) (5.33) の $p_z \to p_v$ の近似の大きさ | [bend_09_algebra_steps.out.txt](outputs/bend_09_algebra_steps.out.txt) |
| `bend_05_help_checks.py` | ヘルプ原文の確認: (a) 線形フリンジの $V$ の TeX 版・平文版 (sympy)、(b) 本体の $s$、(d) 非線形フリンジの $p_{y2}$ の $1/p$、(e) 本体の $p_{x2}$ の TeX 版 | [bend_05_help_checks.out.txt](outputs/bend_05_help_checks.out.txt) |

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
- **Maxwell フリンジの $y^4$ の係数** ($a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$): $f$ の $y^4$ の係数 $-1/(24\rho_b^3)$ と $a$・$b$ が 1 つの $f$ から出る (整合) ことは §6.2 で確認した。係数そのものは $1/\rho_b^3$ の量で、$1/\rho_b$ の 1 次の導出には出ない (1 次は $y$ について厳密)。導出に足りないもの (§6.2 の 1〜4): (a) 端の幅 $d$ が有限の模型での $1/\rho_b^3$ までの摂動 (または厳密な軌道) とその $d \to 0$ の極限 ($y/d \to \infty$ は (5.31) の級数の収束域の外)、(b) $d$ とともに発散する量と消える量の積のうち、硬い端と F1 のどちらに割り振るかの約束 (ヘルプ・K64 が採った約束が資料に無い)。$a$・$b$ が $\sin\theta/\theta$、$(1 - \cos\theta)/\theta^2$ ($\theta = y/\rho_b$) の $\theta^2$ までの展開に等しいことは [06_edge_fringe.md](06_edge_fringe.md) §8.4 にある (導出ではない)。ヘルプ・K64 の係数を使っただけで、独立な導出も数値の検証もしていない (検証した磁場では $y^2/\rho^2$ が $10^{-9}$ 以下で見えない)。(5.34) は $p_y$ を入れない $p_v$ を使う近似で、$p_y$ の補正は相対 $-p_y^2/(2p_v^2)$ ([06_edge_fringe.md](06_edge_fringe.md) §4.3)。
- **F1 の式**: 1 つの F1 から $y$・$y^3$・$x$ の 3 つの効果を作るのは、直線の立ち上がりで厳密、$\tanh$ では $y^3$ の項が同じ形で $x$ の項が $\pi^2/9$ 倍、誤差関数では $y^3$ の項が -4.5% (§6.3。積分は導出して確認)。K64 の F1 は端の形を 1 つの長さで表した近似。斜めの入射での $\Delta X$ の因子 $1/\cos^3\psi$ は導出して ODE で確認したが ((5.43a)、§6.4)、K64 の `dxfr1` にこの因子は無い (gs との比較は未実施)。$\Delta p_y$ の $y$・$y^3$ の項の斜め入射の因子は未導出。
- **$L'$ の補正 (5.43)**: $\psi_1 = \psi_2$ (扇形・対称な面) の $\mathrm{F1}^2$ の 1 次を、端の $\Delta X$ (入口・出口を ODE で確認) と円の幾何から導出し、(5.43) に一致することを示した (§6.4)。一般の面では (5.43c) を導き、実験室系の ODE と 5 桁で一致した。(5.43) との違いは、斜め入射の因子 $1/\cos^3\psi$ ($\psi = 0.3$ rad で 15%) と、非対称な面の $\cos(\omega/2)/\cos(\omega/2 + (\psi_1 - \psi_2)/2)$ の比。非対称・斜め入射では、半径 1 つで出口の向きと位置を同時に満たせず、物理的な端の磁場の設計軌道は出口で $0.1\,s_I$ 程度 (F1 = 0.1 m、$\rho = 10$ m で $5 \times 10^{-6}$ m) ずれる。未確認: $K0 \neq 0$ (出口の向きの条件を置く根拠が無く、定義がヘルプに無い)、$\mathrm{F1}^2$ の解析的な 2 次 (ODE は $\mathrm{F1}/\rho_b = 0.006$ で $10^{-3}$ 以下)、実際の端の磁場の曲率・$y$ 依存、gs との比較。
- **F1 の領域の 2 分割** (`rbh`・`rbl`): 2 点ガウス・ルジャンドルの節であり、直線の立ち上がりの $\int g^k$ ($k \le 3$) を 2 段の階段で再現することを示した (§6.4)。放射 (`tradk`) が偏向角から磁場を読む、という前提と、設計者の意図は確認していない。
- **`tbendcore` の本体の各行**: `t2`・`t3`・`t2t3`・`px1px3`・`dpz32`・`t4`・`dpx3`・`x3`・`da` の意味と導出は §7 の表と §7.1 にある (恒等式は 40 桁で確認)。`x3` と `da` の長い括弧は §7.1 で紙の上で最後まで導出した (加法定理と和積公式だけ。40 桁の数値確認つき)。
- **K1 ≠ 0**: `tbendicorr` の各行 (運動学の残り) と K1 の非線形な蹴りの級数 (A59) は §8 で導出した。Strang 型の合成 (5.50) 全体の精度、端面の薄い蹴りと厳密な端の差 (B117)、写像全体の厳密な基準との比較は、この章では行っていない。分割数の式 (§8) はコードを読んだもので、gs での実測はしていない。
- **CALC と追跡の違い** (§9): F1 ≠ 0 のとき、写像のヤコビアンが k64 で $10^{-7}$ の桁、大きいずれでは $4 \times 10^{-5}$ ($\partial x_{\rm out}/\partial p_y$) で違う。v2 (B324 の修正を含む) では $2 \times 10^{-7}$ 以下に下がる。残り ($10^{-8}$〜$10^{-6}$) の原因と、下がった分の内訳は特定していない。物理的にどちらが正しいかも確認していない。
- **ヘルプと Fortran の Maxwell フリンジの式**: ヘルプの $x_2$ の式と `tbendcore` の `ff` は代数的に等しい (§6.2)。$p_{y2}$ は、ヘルプが $1/p_1$ を余分に持つ (付録 A。誤植か規格化の省略かは確定できない)。`f`・`fa`・`ff` の対称な書き方にした意図はコードから確認できない。
- **B-1**: `tbedgemaxwell` の因子 $p$ は gs の COD で確認済み (4 次元の閉軌道、$\delta = 0.05$〜$0.2$・$-0.1$、$y \le 8$ mm、3 GeV/c、$|\rho| = 5$ m)。$\Delta z$ は 4 次元の CALC に効かないので、gs の 6 次元の量との比較はしていない。更新を $p^2$ 型にすると追跡と合う (§9) が、残りの差 ($z$ の $1.4 \times 10^{-9}$ など) の原因は未特定。
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
| `rbh`・`rbl` | $\tfrac12 + 1/\sqrt{12}$、$1/(6\,\mathrm{rbh})$ | `tbend.f:12-13` | `_RBH`・`_RBL` (`track6.py:1429`、`lattice.py`) |
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
| `tanp1`・`tanp2` | $\tan\psi_1$、$\tan\psi_2$ (`sinp1/cosp1`)。数学の標準の記号ではなく、K64 の `tbendi` の局所変数。K1 ≠ 0 の端面の薄い蹴り $\Delta q_x = \tan\psi\cdot x/\rho_e$ に使う。形の理由: (5.25) の端の行列 $E_x(\psi) = \begin{pmatrix}1 & 0\\ \tan\psi/\rho & 1\end{pmatrix}$ の $p_x$ の行 ($\Delta p_x = \tan\psi\,x/\rho_b$) を $q = p_x/p$ に直したもの ($\rho_e = \rho_bp$) | `tbendi.f:98-108` | `tanp1`・`tanp2` |
| `fpx` | $\texttt{af}\cdot q_x$ ($\texttt{af}$ は DISFRIN = 0 のとき 1、そうでなければ 0。`fringe`)。端の Maxwell フリンジの垂直の蹴りの $p_x/p_v \simeq q_x$ ((5.33)。$q_x \ll 1$ の近似)。`pyi = py - (tanp1 + fpx)*f` ($f = y/\rho_e$) は、幾何の端面の収束 $-\tan\psi\,y/\rho_e$ ((5.25) 直後の $E_y$) と、Maxwell の $-q_xy/\rho_e$ の和。`ff = af*y*f/2` $= y^2/(2\rho_e)$ は (5.35) の $\Delta x$ ($p^2/p_v^3 \to 1/p$)、`z -= ff*fpx` は (5.35) の $\Delta z$ ($-y^2q_x/(2\rho_ep)$ 型) | `tbendi.f:98-108` | `fpx`・`ff` |
| `xr` | $x/\rho_0$ ($= hx = u$。§8 の $(1 + u)\ln(1 + u)$ の級数の変数)。K1 の非線形な蹴り `ak1*xi*xr*(0.5 - xr*(2 - xr)/12)/p` に使う | `tbend.f:118-134` (`tbendicorr`) | `xr` |
| `ak1` | `tbendicorr` の引数。隣り合う 2 つの半スライスをまとめた積分 K1 ($(\texttt{akxn}(n) + \texttt{akx0})/2$)。§8 の蹴りの係数 $K\ell$ に当たる | `tbend.f:118`、`tbendi.f:115-121` | `ak1` |
| `akxn`・`akx0` | `akxn(n)` $= \texttt{ak}\cdot\texttt{alx}(n)/\texttt{al}$ は $n$ 番目のスライスの積分 K1 ($K\ell_n$)。`akx0` は 1 つ前のスライスの値 (最初は 0、最後の半分は `akx0*.5`)。(5.50) の合成で、隣り合う $\mathcal{C}(\ell/2)$ を 1 つにまとめるために平均する | `tbendi.f:115-121`、`:147` | `akxn`・`akx0` |
| `psqmax` | $q_x^2 + q_y^2$ の上限。追跡では 0.9999 (`min(0.9999, pxi**2 + pyi**2)`)。$\sqrt{1 - s}$ ($\texttt{dpzi} = p_z/p - 1$) が負や 0 にならないための切り捨て。CALC の `tbendecorr` の引数 | `tbend.f:118-134`、`perfkern_bend.py:65` | `psqmax` |
| `phi` | $h\ell = \ell/\rho_0$ (スライスの偏向角。`phixn(n) = phi0*alx(n)/al` を隣と平均したもの)。§8 の `xi` の更新の分母 $\texttt{pzi} - \texttt{phi}\,q_x$ と、`pxi += phi*dpzi` に使う ($q_x' = h\,\texttt{dpzi}$ を $\ell$ 倍した量) | `tbend.f:118-134` | `phi` |
| `s` | `tbendicorr` の局所変数 $= \min(\texttt{psqmax}, q_x^2 + q_y^2)$。`dpzi = sqrt1(-s)` の引数 (§8 の「運動学の残り」)。§8 の本文の $a$・$b$ は `xi` の方程式 $x' = a + bx$ の係数を指す説明用の記号で、Fortran の変数ではない | `tbend.f:118-134` | `s` |
| `a`・`b` | `tbedgemaxwell` の局所変数: $a = 1 - y^2/(6\rho_b^2)$、$b = \tfrac12 - y^2/(24\rho_b^2)$ ((5.36)。`af = a*y/rhob/pvi`、`bf = b*y**2/rhob/pvi**3` として COD の更新に使う)。数学の標準の記号ではなく、$y^4$ の補正の係数 (§6.2。06 章 §8.4 の $\sin\theta/\theta$ との対応) | `tbfrie.f:190-229` | `a`・`b` (`lattice._tbedgemaxwell_cod`) |
| `epsbend` | $10^{-3}$ (K1 ≠ 0 の分割の角) | `tbend.f:12` | 変数名なし (`track6._tbendi` の `eps = 1e-3 * eps0`、`track6.py:1287`) |
| `ndiv` | スライスの数 | `tbendi.f:41`、`tbend.f:640`、`tbende.f:460` | `ndiv` |

確認: `tbend.f` の `ky_FRIN_BEND` が DISFRIN に対応することは `initb1.f:925` (`kytbl(kwFRIN,0)=sethtb('DISFRIN ',...)`) で確認した。

## 付録 A. ヘルプの記述・K64 との違い (偏向磁石)

| 箇所 | 内容 |
|---|---|
| ヘルプ BEND > transformation、body | **原文で確認、確定した食い違い**: $s \equiv \mathrm{ANGLE}\times\rho_b\times(\arcsin\ldots)$ の $\mathrm{ANGLE}\times$ は誤記と読む (TeX 版・平文版とも同じ式)。面内の経路長は $s = \rho_bp_h\Theta$ ($p_h = \sqrt{p^2 - p_y^2}$)。そう読めば $y_2 = y_1 + s\,p_y/\sqrt{p_1^2 - p_{y1}^2}$ が (5.24) と一致 (§5.3)。$p_{x2}$・$x_2$ の式は (5.22)(5.23) と同じ (実験室系の幾何と $2.4 \times 10^{-15}$)。$p_{x2}$ の $x_1$ の項は TeX 版だけ誤植 ($-x_1/(\rho_b\sin\omega)$。平文版 $-(x_1/\rho_b)\sin\omega$ が正しい) |
| ヘルプ BEND、nonlinear fringe の $p_{y2}$ | **原文で確認、食い違い (誤植または規格化の省略。どちらか確定できない)**: 入口 $p_{y2} = p_{y1} - p_{x1}(1 - y_1^2/(6\rho_b^2))\,y_1/(p_1\rho_b\sqrt{p_1^2 - p_{x1}^2})$ (出口は符号が逆。TeX・平文とも同じ) は、$p_1$ を 1 つ余分に持つように見える。正準運動量での導出 (§6.2) は $-(p_x/p_v)\,y/\rho_b\cdot(1 - y^2/(6\rho_b^2))$ で、原文の式との比は $1/p$ ([bend_05_help_checks.out.txt](outputs/bend_05_help_checks.out.txt) (d)、数値の ODE は [bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b)。ODE は正準運動量で、$\Delta p_y$ が $-(p_x/p_v)y/\rho_b$ と相対 $1.5\times10^{-4}$ で一致)。K64 の式 (`f` の $\rho_e = \rho_bp$ は $q$ の単位) は導出と同じ。ヘルプの $p_{y2}$ と $p_{y1}$ が ($q$ の単位の) $p_y/p$ なら、$p_1$ の 1 つ余分な因子は $q$ への換算として一致する。ただし同じ式の $x_2$ は正準の $p_1$ で書かれていて (正準なら $p^2/p_v^3$ の形が正しい)、ヘルプは規格化を明記していない。したがって「ヘルプの誤り」とは断定せず、$p_1 = 1$ では一致する、とだけ言える |
| ヘルプ BEND、linear fringe の変換 $V$ と $p_{y2}$ | **原文で確認**: 入口の $V$ は **TeX 版と平文版で式が違う**。平文版 $V = -f^2p_x/(24\rho_bp) - fy^2/(12\rho_b^2p) + y^4/(6\rho_b^2fp)$ は、ヘルプ自身の $x_2$・$p_{y2}$ と一致する (ヘルプの $\exp(:V:)$ は $V$ をハミルトニアンとする流れ。C 付録 §8.2)。(5.42) の $f_{\rm lin}$ とは、$-V = f_{\rm lin} + \mathrm{dxfr}\,p_x$ の関係 ($p$ に依らない $x$ のずれ $-\mathrm{F1}^2/(24\rho_b)$ の項だけ違う。§6.3)。TeX 版 $-\tfrac1p\bigl(\tfrac{F1^2}{24\rho_b}p_x - \tfrac{F1}{12\rho_b}y^2 + \tfrac{y^4}{6F1\rho_b^2}\bigr)$ は、$y^2$ の項が $\rho_b$ の 1 乗 (次元が合わない) で符号も逆、$y^4$ の項の符号も逆の誤植 ([bend_05_help_checks.out.txt](outputs/bend_05_help_checks.out.txt) (a))。$p_{y2} = p_{y1} + y_1(\Delta y_{fr} - \Delta y_{fra}y_1^2)/p_1^2$ の $1/p_1^2$ は、$q$ の単位 ($p_y/p$) の更新と読めば K64 と一致。正準運動量では $1/p$ (§6.3)。$z_2$ の $p_{x1}$ も $q$ の単位と読むと (5.42) の $\partial f/\partial p$ と一致 |
| ヘルプ BEND の F1 | **原文で確認**: 厳密な定義は $6\int(g - g^2)\,ds$ (5.40) (「1 つのフリンジについて積分」)。四極の F1 (04 章の $I_1$ のモーメント) と**定義が違う**。偏向磁石の F1 は $y$ の収束の量で、$x$ のずれの $\mathrm{F1}^2/24$ は直線の立ち上がりでだけ厳密 (§6.3)。「Only the effects up to $y^4$ in Hamiltonian」も原文どおり |
| ヘルプ BEND の末尾 | **原文で確認、実装との不一致の疑い (未確認)**: 「In FFS optics and Emittance calculations, or when the synchrotron radiation is turned on, the same algorithm as K1 <> 0 is applied.」。この章 §9 の記述 (CALC は K1 = 0 で厳密な `tbedge` の鎖と `ndiv = 1` を使う。`tbende.f:460`) とは合わない。K64 のコードでは K1 = 0 の枝 (`ak == 0`) が別にあり、ヘルプの記述が実装と違う可能性が高いが、gs で確かめていない |
| ヘルプ BEND、DISFRIN | **原文で確認**: 「If nonzero, the nonlinear Maxwellian fringe is suppressed.」 (DISFRIN の項)。K1 = 0 の追跡 (`tbendcore`) では、DISFRIN は F1 の $y^3$ の項 (`dyfra`) だけを止め、硬い端の Maxwell フリンジ ((5.36) の `f`・`ff`) は止まらない ([bend_02_arc_faces.out.txt](outputs/bend_02_arc_faces.out.txt) (e)、[bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (e))。CALC も DISFRIN によらず Maxwell を入れる (docs/30 の B134)。K1 ≠ 0 の追跡は `ttfrin` (K1 の端) が DISFRIN に従う |
| ヘルプ BEND の K1 ≠ 0 の分割数 | **原文で確認、確定した食い違い (ヘルプと実装が違う)**: 原文は 「$1 + \mathrm{floor}(\sqrt{|K1\,L'|/(12\cdot10^{-5}\,\mathrm{EPS})})$、EPS = 0 のときは EPS = 1」。コード `tbendi.f:41` は `ndiv = 1 + int(abs(phi0/eps))`、`eps = epsbend*eps0` (`epsbend = 1e-3`、EPS = 0 のとき `eps = epsbend`) で、式が違う (角で刻む)。追跡 (`tbendi.f:41`)・CALC (`tbende.f:460`) とも $1 + \mathrm{int}(|\phi/\varepsilon|)$、$\varepsilon = 10^{-3}\mathrm{EPS}$ (ソースを読んで確認)。放射があるときは `ndivrad` との大きい方 (gs での実測はしていない) |
| ヘルプ BEND、K1 ≠ 0 の端面 | **原文で確認** (「the effects from E1 and E2 are approximated by thin quadrupoles」)。「E1・E2 の効果は薄い四極で近似する」。K64 のコードは、端面の薄い蹴りは偏向の蹴り (`tanp1*xi/rhoe`) と、垂直の Maxwell 型の蹴り (`(tanp1 + fpx)*f`) で、四極そのものではない (`tbendi.f:98-108`)。CALC は厳密な `tbedge` (B117) |
| ヘルプ BEND > 出口の面からの直進の $\psi_2$ | **原文で確認、食い違い (省略)**: $\psi_2 \equiv \mathrm{ANGLE}\times\mathrm{E2}$ と定義していて $\mathrm{AE2}$ が無い。入口の $\psi_1 \equiv \mathrm{ANGLE}\times\mathrm{E1} + \mathrm{AE1}$ と、AE2 の説明 (有効な面の角は $\mathrm{E2}\times\mathrm{ANGLE} + \mathrm{AE2}$)、同じ節の $\Delta z_{fr}$ の式 ($\sin(\mathrm{ANGLE\,E2} + \mathrm{AE2})$) から、$\psi_2 = \mathrm{ANGLE}\times\mathrm{E2} + \mathrm{AE2}$ の書き落とし。 |
| ヘルプに無いもの | (1) 円弧の厳密解の導出 (§5)。(2) $\omega = \mathrm{ANGLE} - \psi_1 - \psi_2$ の意味と、面の座標系 (§5.3)。(3) Maxwell フリンジの導出 (§6.2)。(4) F1 の 3 つの係数の導出 (§6.3)。(5) $L'$ の導出 (§6.4。一般の面・斜め入射まで)。(6) K1 ≠ 0 の本体の $H_2$ と、`tbendicorr`・`tbendibody` の式 (§4・§8)。(7) 行列と追跡の関係 (§9) |

## 付録 B. K64 の誤り (確定 = B-1・B-4、疑い = B-2・B-3)

確定・疑い・取り下げの区別: **確定** = B-1 (gs で確認。docs/30 の B324)、B-4 (ODE で確認。B250)。**疑い (未決)** = B-2 (DISFRIN の仕様)、B-3 (ヘルプの誤記)。**取り下げ**: この章には無い。

この章は説明が中心で、v2 (既定) が直す動作は B-1 と B-4。B-2・B-3 は疑い (確定していない)。

**B-1. CALC の `tbedgemaxwell` の COD の更新が $\Delta x$ と $\Delta z$ で因子 $p = 1 + \delta$ を落としている (確定。gs で確認。docs/30 の B324。v2 は修正)**

- **内容**: `tbfrie.f:209` の `cod(1)=cod(1)+bf*pr` と `:211` の `cod(5)=cod(5)-bf*pxi` は、(5.36) の $\Delta x = b\,y^2p^2/(\rho_bp_v^3)$、$\Delta z = -b\,y^2p_xp/(\rho_bp_v^3)$ より因子 $p$ だけ小さい ($\texttt{bf} = b\,y^2/\rho_b/p_v^3$ に `pr` を 1 回しかかけない。(5.36) は `pr**2`、$z$ は `pxi*pr`)。$\Delta p_y$ は正しい。
- **根拠**: (i) 端の磁場 (5.31) の中のローレンツ力の積分で、$\Delta x$ は $\delta = 0.3$ で $p^2$ の式と一致 (比 1.0000)、`pr` の式とは 1.3000 倍違う ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (b))。(ii) 追跡 (`tbendcore`) の式は $p^2$ の形で、6×6 ヤコビアンがシンプレクティック ($10^{-12}$)。CALC の更新はシンプレクティックでない ($\delta = 0.1$ で $1.7 \times 10^{-4}$。$\delta = 0$ でも、$\delta$ の列の微分がずれて $2.0 \times 10^{-5}$) ([bend_03_fringe.out.txt](outputs/bend_03_fringe.out.txt) (c))。(iii) CALC の解析行列 `trans(1,3) = af*pr**2/pvi2` は $p^2$ の式の微分で、COD の更新の数値微分と $1 + \delta$ 倍違う (内部の矛盾。[bend_04_matrix_vs_map.out.txt](outputs/bend_04_matrix_vs_map.out.txt) (c2))。(iv) **gs で確認** (下の「gs での確認」)。(v) 更新だけを $p^2$ 型に差し替えると、CALC の COD の更新と追跡の写像の差のうち $x$・$p_x$ は $10^{-15}$ 以下に、$z$ は $2.0 \times 10^{-8}$ から $1.4 \times 10^{-9}$ に下がる ([bend_07_b1_estimate.out.txt](outputs/bend_07_b1_estimate.out.txt) (b))。
- **大きさ**: 1 つの面での $\Delta x$ の差は $\dfrac{b\,y^2}{\rho_b}\dfrac{p^2 - p}{p_v^3} \simeq \delta\,y^2/(2\rho_b)$ 程度で、$y = 10$ mm、$\delta = 0.05$、$\rho_b = 10$ m で $2.3$〜$2.5 \times 10^{-7}$ ($\delta y^2/(2\rho_b) = 2.5 \times 10^{-7}$)。入口と出口で符号が逆 ($f_{\rm out} = -f_{\rm in}$) なので、要素全体の COD の差は、入口と出口の $y^2$ の違い (扇形で 4%) と、入口の差が本体の $z' = -hx$ で $z$ に効く分 ($\varepsilon L/\rho_0 \approx 2 \times 10^{-8}$) が残り、扇形の $x$ で $9.7 \times 10^{-9}$、$z$ で $2.0 \times 10^{-8}$ になる。つまり $2.5 \times 10^{-7}$ は面ごとの値で、要素全体の差 ($2 \times 10^{-8}$) とは桁が違うのが正しい。F1 ≠ 0 の 4×4 の差 ($4.4 \times 10^{-5}$。k64) は、F1 = 0 の差とは桁が違い、v2 では $2 \times 10^{-7}$ 以下に下がる (§9。F1 = 0 では 4×4 の差は $3 \times 10^{-8}$ 以下)。(c2) の $4.5 \times 10^{-5}$ は 1 段の解析行列と更新の数値微分の食い違いで、数値が近いのは偶然とみなす。
- **gs での確認**: 4 台の BEND ($\mathrm{ANGLE} = \pm0.2$、$L = 1$ m、$\rho = \pm5$ m、FRINGE = 1、DISFRIN = 0、$E_1 = E_2 = 0$ と $0.05$) と QUAD 2 台 (DY $= \pm3$ mm で垂直の閉軌道、$y = 4$〜$8$ mm) の 1 周 (3 GeV/c、`NORFSW`)。DP0 を与えて `CALC` した閉軌道 $c$ から `TrackParticles` を 1 周し、追跡の写像の不動点を求めて CALC の閉軌道との差 $dc$ (4 次元) を出した。$\delta = 0$ は $dc_x, dc_{p_x} < 10^{-12}$。$\delta \ne 0$ では $dc_{p_x}$ が $\delta = 0.05$ で $8.5 \times 10^{-8}$、$\delta = 0.2$ で $3.8 \times 10^{-7}$ ($E_1 = E_2 = 0.05$ でも $3.8 \times 10^{-7}$)、$\delta = -0.1$ で $-1.5 \times 10^{-7}$、$\delta = 0.1$・$y = 8$ mm で $6.9 \times 10^{-7}$ で、$\delta y^2/(2\rho)$ の大きさ。pysad の k64 は gs と COD が $10^{-15}$、$dc$ が $10^{-16}$ 以下で一致する。COD の更新を $p^2$ 型にした v2 では $dc_{p_x}$ が $\delta = 0.05$ で $3.8 \times 10^{-10}$、$0.2$ で $3.2 \times 10^{-8}$、$0.1$・$y = 8$ mm で $5.3 \times 10^{-9}$ に下がる (k64 の 0.1〜9%。$\delta = 0$ は変わらない)。残りは端の Maxwell 以外の $\delta$ の項 (原因は未特定)。試験 `test_v250_bend_maxwell_calc_p.py`、記録 `docs/history/B324_bend_maxwell_calc_p.md`。
- **状態**: **確定** (gs の COD、ソースの読み、独立なローレンツ力の積分、解析行列との整合)。**v2 は COD の更新を $p^2$・$p^1$ 型に直した** (`lattice._tbedgemaxwell_cod`、`elkern._tbedgemaxwell`、追跡の端の `track6n._edge_maxwell1`)。k64 は K64 のまま (gs と同じ)。限界: 4 次元の閉軌道で確かめた。$\Delta z$ は 4 次元の CALC に効かないので、(5.36) の閉じた式と解析行列の数値微分で確かめた (gs の 6 次元の量との比較はしていない)。$\delta$ は 0.05〜0.2・$-0.1$、$y \le 8$ mm、$|\rho| = 5$ m、3 GeV/c。

**B-2. DISFRIN が K1 = 0 の追跡の Maxwell フリンジを止めない (ヘルプの記述との食い違い。仕様の疑い)**: 付録 A。CALC と追跡の両方で止まらないので、両者の間では一貫している (B134 と同じ種類)。

**B-3. ヘルプの式の記述の食い違い (K64 のソースの誤りではなく、ヘルプの誤記・規格化の省略の疑い)**: 付録 A の $s$ の $\mathrm{ANGLE}\times$、nonlinear fringe の $p_{y2}$ の $1/p_1$。

**B-4. K1 ≠ 0 の本体の `xsxkx` の符号 (docs/30 の B250 の既知の誤り)**: (5.11) の導出から、`akxsq < 0` の枝の `xsxkx = -xspx/akx/akxsq` は符号が逆で、連続な式 $(\ell - \texttt{sxkx})/\texttt{akxsq}$ ($= \texttt{xspx}/\texttt{akx}/\texttt{akxsq}$) が正しい。$H_2$ の ODE との比で $z$ の誤差 $7.3 \times 10^{-5}$ m ($\rho_0 = \rho_b = 10$ m、$L = 1.3$ m、$\delta = 0.01$、K1 ≥ 0) が、符号を直すと $10^{-18}$ ([bend_01_h2_body.out.txt](outputs/bend_01_h2_body.out.txt) (d)(e))。v2 は直してある。
