# 02. ドリフト (DRIFT)

**この章で分かること**: 磁場のない区間 (ドリフト) のハミルトニアンから、(1) 厳密な追跡の写像 ($x' = p_x/p_z$、$z$ の更新、$dv$ の項)、(2) 近軸近似 (DISKIN) との関係、(3) 閉軌道 (COD) のまわりで線形化した 6×6 の転送行列 (K64 の `tdrife1`) を、途中の式を省かずに追う。最後に、いちばん短い例として、ドリフトだけで起きる運動量依存の遅れ ($R_{56}$ に当たる量) を数値で示す。ヘルプが書いていない途中 (行列の各成分、$z$ の式の整理、$dv$ の扱い) を埋める。式番号 $(2.n)$ から、§11 の Fortran の対応表と §12 の pysad の対応表を引ける。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) (§1〜§2 と §4 の $dv$、§5 の $q = p_x/(1+\delta)$、§8 の変数名の注釈の規約)。偏微分と行列。ハミルトニアンの導出は [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) にあるが、この章では式 (1.1) の結果だけを使う。

**読む順**: §1 (結論) → §2 (記号) → §3〜§6 (写像の導出) → §7 (近軸との関係) → §8〜§9 (行列) → §10 (例)。コードを読みながら使うときは §11・§12 から入る。四極・偏向磁石・多極は、この章の写像を「刻みの途中の移動」として使う ([04_quad.md](04_quad.md)、[05_bend.md](05_bend.md)、[07_thin_multipole.md](07_thin_multipole.md))。

## 要約と目次

**この章で分かること**

- 磁場のない区間 (DRIFT) のハミルトニアン $H = -\sqrt{p^2 - p_x^2 - p_y^2} + E/(P_0v_0)$ (式 (2.1)) から、厳密な追跡の写像と、閉軌道 (COD) のまわりの 6×6 転送行列を、途中の式を省かずに導く。
- 結果の主な式: 厳密な写像 $x \to x + (p_x/p_z)L$、$z \to z - \bigl[(p_x^2 + p_y^2)/((p + p_z)p_z) + dv\bigr]L$ (式 (2.6)、$p_z = \sqrt{p^2 - p_x^2 - p_y^2}$)、厳密 = 近軸 × $D$ (式 (2.10))、行列の成分 (2.11)〜(2.13)、$z$ の行を作るシンプレクティック条件 (2.15)、設計軌道上の $R_{56} = L\,h_0/h_1^3 \to L/\gamma_0^2$ (式 (2.16)。$h_0 = \gamma_0$、$h_1$ はその粒子のローレンツ因子)。
- 対応: 追跡は K64 の `tdrift.f` の `tdrift_free` (pysad は `track6._tdrift_free`、numba 版 `track6n.drift1`)、行列は `toplvl.f` の `tdrife1` (pysad は `lattice._tdrife1_mat`、`elkern._jac_drift`)。`DISKIN` (近軸のドリフト) は pysad が読むだけで未実装。
- 確かめ方: sympy (偏微分・6×6 のヤコビアン)、DOP853 による正準方程式の積分 (差 $\le 1.8\times10^{-16}$)、pysad の追跡の数値微分 (行列との差 $\le 3.9\times10^{-11}$)。K64 の実行 (gs) との比較はしていない。
- K64 の誤り: なし (`dpz` の名前の二重の意味と、`DISKIN` が追跡だけに効くことは仕様で、B 番号なし)。ヘルプの誤り: **確定** (原文の式どうしの矛盾) として DRIFT の $E/v_0$ の符号 (QUAD・dynamics・BEND の式と逆。gs での確認はしていない)。疑い: なし。ヘルプに無いのは写像の式・行列・`DISKIN` (付録 A)。

**導出の流れ**

1. [§1 結論](#1-結論): 写像、$R_{56}$、近軸との関係、行列、`DISKIN` の結果を先に述べる。
2. [§2 記号と単位](#2-記号と単位): 長さ・運動量・$dv$ などの対応表。`dpz` が追跡と行列側で別の量であることに注意。
3. [§3 ハミルトニアン](#3-ハミルトニアン): 01 章の (1.1) を $1/\rho = 0$、$a = 0$ で書き (2.1)、$E/(P_0v_0) = h_0h_1/p_0^2$ (2.2) に直す。ヘルプの式との照合。
4. [§4 運動方程式](#4-運動方程式-偏微分を-1-項ずつ): $H$ の偏微分 (2.3) を 1 項ずつ書き、(2.4) を得る。$H$ が $x, y, z$ を含まないので右辺は定数。
5. [§5 $z$ の整理と $dv$](#5-zzmathzz195end-の整理と-zzmathzz196end): $-p/p_z + 1$ を有理化して桁落ちを避け、(2.5) にする。
6. [§6 厳密解と追跡の式](#6-厳密解と追跡の式): 積分して (2.6)、規格化した運動量での形 (2.7)、K64 の `tdrift_free` の形 (2.8)。`pxy2dpz` の級数と有理形の導出。数値確認。
7. [§7 近軸のドリフトと運動学的ドリフト、DISKIN](#7-近軸のドリフトと運動学的ドリフトdiskin)
    - [§7.1 厳密 = 近軸 × $D$](#71-厳密--近軸--zzmathzz283end): $H = H_2 + \Delta H$ の 2 つの流れが交換し、(2.9) → (2.10)。
    - [§7.2 DISKIN](#72-diskin): K64 が近軸の式だけを使う分岐と、pysad で未実装であること。
8. [§8 COD のまわりの 6×6 行列](#8-cod-のまわりの-66-行列)
    - [§8.1 偏微分を 1 項ずつ](#81-偏微分を-1-項ずつ): (2.6) を $p_x, p_y, \delta$ で微分して (2.11)・(2.12)。
    - [§8.2 K64 の `tdrife1`](#82-k64-の-tdrife1): (2.13) との一致、COD の更新 (2.14)、`tdrife` の掛け足し。
9. [§9 シンプレクティック性と $z$ の行](#9-シンプレクティック性と-zzmathzz442end-の行): $M^\top JM = J$ から (2.15) を導き、$T_{52} = T_{16}$ などを得る。前提 $T_{66} = 1$。
10. [§10 最短の例](#10-最短の例-ドリフトだけの運動量依存の遅れ): (2.16) と、3 GeV・10 MeV の電子の数値。
11. [§11 Fortran の対応表](#11-fortran-の対応表)・[§12 pysad の対応表](#12-pysad-の対応表)・[§13 確認スクリプトと出力](#13-確認スクリプトと出力)・[§14 未確認・限界・使った条件](#14-未確認限界使った条件)・[§15 この章で使う Fortran の変数名](#15-この章で使う-fortran-の変数名)・[付録 A](#付録-a-ヘルプの記述k64-との違い-ドリフト)。

式の連鎖: (2.1) → (2.2) → (2.3) → (2.4) → (2.5) → (2.6) → (2.7)・(2.8)。近軸との関係は (2.9) → (2.10)、行列は (2.6) → (2.11)・(2.12) → (2.13)、さらに (2.15) → (2.16)。

**読み方**: 前提は [01](01_variables_and_conventions.md) の §1〜§2、§4 ($dv$)、§5 ($q = p_x/(1+\delta)$)、§8 (記号の規約)。写像だけなら §1〜§6、行列は §8〜§9。コードから入るなら §11・§12。四極・偏向磁石・多極の章は、この章の写像を刻みの途中の移動として使う。

**未確認・限界の所在**: §14 (`DISKIN` の gs での未比較、行列の gs との未比較、$p_z \to 0$ の限界など) と付録 A (ヘルプの $E/v_0$ の符号、`DISKIN` の未実装)。

## 1. 結論

- **ハミルトニアン** (01 章の式 (1.1) で $1/\rho = 0$、$a_x = a_y = a_s = 0$): $H = -\sqrt{p^2 - p_x^2 - p_y^2} + E/(P_0 v_0)$。$H$ は $x, y, z$ を含まないので、$p_x, p_y, p$ は一定で、座標は一定の速さで進む。**厳密に解ける**。
- **厳密な写像** (長さ $L$。$p_z \equiv \sqrt{p^2 - p_x^2 - p_y^2}$):

$$
x \to x + \frac{p_x}{p_z}L,\qquad
y \to y + \frac{p_y}{p_z}L,\qquad
z \to z - \left[\frac{p_x^2 + p_y^2}{(p + p_z)\,p_z} + dv\right]L
$$

  ($p_x, p_y, \delta$ は変わらない。式 (2.6)。$z$ の第 1 項は、傾いて進む分の経路の伸び、第 2 項は速度の差 $dv = 1 - v/v_0$ による遅れ。)
- **設計軌道上の粒子 ($p_x = p_y = 0$) は $\delta$ に比例して $z$ が動く**: $\Delta z = -dv\,L \simeq L\,\delta/\gamma_0^2$ ($\gamma_0$ は設計のローレンツ因子)。電子 3 GeV ($\gamma_0 = 5870.85$) で、1 m あたり $\delta = 10^{-2}$ のとき $+2.858 \times 10^{-10}$ m (近似式 $L\delta/\gamma_0^2$ では $2.901 \times 10^{-10}$。差は $\delta$ の 2 次の項)。10 MeV の電子 ($\gamma_0 = 19.6$) では同じ $\delta$ で $+2.566 \times 10^{-5}$ m (近似式では $2.604 \times 10^{-5}$) ([drift_03_r56_example.out.txt](outputs/drift_03_r56_example.out.txt))。これが転送行列の $R_{56} = \partial z/\partial\delta = L\,h_0/h_1^3$ (式 (2.16)) で、「ドリフトの $R_{56}$」。
- **厳密なドリフトは、近軸のドリフトと、運動学的ドリフト $D$ (四極の章の $\Delta H$ の流れ) の積**で、両者は交換する (式 (2.10))。四極の本体の分割 ([04_quad.md](04_quad.md) §5〜§6) の $D$ は、$K = 0$ のときのこの差に当たる。
- **COD のまわりの 6×6 行列**は、(2.6) を正準変数で偏微分したもの (式 (2.11)〜(2.13))。$x, y, z$ の行だけが単位行列と違い、列は $p_x, p_y, \delta$ だけ。COD に傾き $p_x$ があると、$R_{16} = -p_x\,p\,L/p_z^3 \simeq -p_x L$ が現れる ($\delta$ が違うと $x$ が動く: 軌道の傾きによる経路・位置の運動量依存)。シンプレクティックで、$z$ の行はシンプレクティック条件からも組める (§9)。
- **DISKIN (近軸のドリフト)**: K64 の追跡は、`DISKIN` ≠ 0 のドリフトで近軸の式 ($x \to x + qL$、$z \to z - (q^2/2 + dv)L$) を使う。pysad の追跡は `DISKIN` を読むが使わず、厳密な式のまま (未実装。§7.2)。
- 数値の確認: 閉じた式は、正準方程式の数値積分 (DOP853) と $2 \times 10^{-16}$、pysad の追跡と $5.6 \times 10^{-17}$ で一致。行列は、sympy の偏微分と $6 \times 10^{-18}$、追跡の数値微分と $4 \times 10^{-11}$ 以下 (§13)。
- **使った条件と限界**は §14 に書く。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 長さ $L$ [m] | `L` | `al` | `L` |
| 全運動量 $p = 1 + \delta$ | | `pr` (`tdrife1`、`tdrife`)。追跡の $\delta$ は `g(i)` | `pr` |
| 縦方向の運動量 $p_z = \sqrt{p^2 - p_x^2 - p_y^2}$ (正準運動量で) | | `pzi` (`tdrife1`)。`pzi = pr*sqrt(1-a/pr**2)` | `pzi` |
| $p_x^2 + p_y^2$ | | `a` (`tdrife1`)、`pxisq`・`pyisq` (`qtwiss.f`) | `a` |
| 規格化した運動量 $q = p_x/p$、$r = p_y/p$ | | 追跡の `px(i)`、`py(i)` | `track6._Particle` の `c[1]`、`c[3]` |
| $\sqrt{1 - q^2 - r^2} - 1$ | | `dpz` (**追跡**の `tdrift_free`。`pxy2dpz(px,py)`) | `_pxy2dpz` |
| $p_z - p$ (**行列側**。canonical で、規格化していない) | | `dpz` (`tdrife1` の周辺、`qtwiss.f:212`。名前は同じで中身が違う) | `dpz` |
| 速度の差 $dv = 1 - v/v_0$、$-\partial dv/\partial\delta$ | | `dv`、`dvdp` (`tgetdv`) | `p.dv`、`_dvz` |
| 座標・運動量・$\delta$ | $\{x, p_x, y, p_y, z, \delta\}$ | `x(i)`, `px(i)`, … (追跡)。`cod(1..6)` (CALC) | `c[0..5]` |

- Fortran の `dpz` は 2 つの意味で使われる (**取り違えやすい**)。追跡の `tdrift_free` では $\zeta - 1 = \sqrt{1 - q^2 - r^2} - 1$ ($q = p_x/p$ は規格化した運動量)。CALC の `qtwiss.f` の `dpz` は $p_z - p$ ($p_x, p_y$ が正準運動量)。どちらも $-(\text{傾きの 2 乗})/2$ に近い小さい負の数。
- 添字は、Fortran は 1 始まり (`cod(1)` が $x$)、pysad は 0 始まり (`c[0]` が $x$)。
- 数値の例は電子、運動量 3 GeV/$c$ ($\gamma_0 = 5870.85$)。`MOMENTUM = 3GEV` から。

## 3. ハミルトニアン

01 章の式 (1.1) を、直線上 ($1/\rho = 0$)・磁場なし ($a_x = a_y = a_s = 0$) で書くと

$$
H = -\sqrt{p^2 - p_x^2 - p_y^2} + \frac{E}{P_0 v_0}
\tag{2.1}
$$

($E$ は粒子のエネルギー、$P_0 v_0$ は設計運動量と設計速度の積。$E/(P_0v_0)$ は $p$ だけの関数で、$z$ と共役な $p$ の側の量。) ヘルプ (K64 版同梱の `SADHelp.html`) の DRIFT > transformation には $H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ ($L$ をかけた形で、$\exp(:H:)$ に入れる。TeX 版・平文版とも同じ式) とある。**原文で確認**: $\sqrt{\ }$ の項と $L$ は (2.1) と同じ形で、定数 $1$ は運動に影響しない。$E/v_0$ は $P_0$ の規格化を略した書き方 (QUAD の $H_{2n}$ も同じ略記)。ただし**$E/v_0$ の符号は (2.1) と逆**で、これは誤植と判断した (ヘルプの式どうしの整合による判断で、SAD (gs) の実行での確認はしていない。付録 A。$\exp(:H:)$ を「$H$ を時間 1 のハミルトニアンとする流れ」と読む規約は [C_poisson_lie.md](C_poisson_lie.md) §8.2 に確定してある)。

$E$ を $p$ で書く。設計運動量を $mc$ の単位で $p_0$ ($= \beta_0\gamma_0$。無次元。次元のある設計運動量は $P_0 = p_0\,mc$。[01_variables_and_conventions.md](01_variables_and_conventions.md) §2)、$h_0 = \gamma_0 = \sqrt{1 + p_0^2}$、$h_1 = \sqrt{1 + p_0^2p^2}$ とする。

1. $E = \gamma mc^2$ で、粒子の運動量は $p_0p\;mc$ なので $\gamma = \sqrt{1 + (p_0p)^2} = h_1$、すなわち $E = h_1\,mc^2$。
2. $P_0v_0 = (p_0\,mc)(\beta_0c) = p_0\beta_0\,mc^2$。$\beta_0\gamma_0 = p_0$ から $\beta_0 = p_0/h_0$ なので $P_0v_0 = (p_0^2/h_0)\,mc^2$。
3. 比を取ると $E/(P_0v_0) = h_1mc^2\cdot h_0/(p_0^2mc^2)$:

$$
\frac{E}{P_0 v_0} = \frac{h_0\,h_1}{p_0^2}
\tag{2.2}
$$

($E/(P_0v_0)$ は無次元。$h_1$ は $p$ に依存。[drift_06_steps.out.txt](outputs/drift_06_steps.out.txt) (a) で sympy 確認。)

**$H$ が含む変数**: $p_x, p_y, p$ だけ。$x, y, z$ も $s$ も含まない。したがって $x' = \partial H/\partial p_x$ などの右辺は定数で (次の §4)、方程式は積分できる。

## 4. 運動方程式 (偏微分を 1 項ずつ)

正準方程式 ($s$ が独立変数。[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §8.3): $x' = \partial H/\partial p_x$、$p_x' = -\partial H/\partial x$、同様に $y$、$z' = \partial H/\partial p$、$p' = -\partial H/\partial z$。$p_z = (p^2 - p_x^2 - p_y^2)^{1/2}$ の偏微分は、鎖則で

$$
\frac{\partial p_z}{\partial p_x} = \frac{-p_x}{p_z},\qquad
\frac{\partial p_z}{\partial p_y} = \frac{-p_y}{p_z},\qquad
\frac{\partial p_z}{\partial p} = \frac{p}{p_z}
\tag{2.3}
$$

$H = T_1 + T_2$、$T_1 = -p_z$、$T_2 = E/(P_0v_0)$ として:

| 項 | $\partial/\partial p_x$ | $\partial/\partial p_y$ | $\partial/\partial x,\ \partial/\partial y,\ \partial/\partial z$ | $\partial/\partial p$ |
|---|---|---|---|---|
| $T_1 = -p_z$ | $+p_x/p_z$ | $+p_y/p_z$ | 0 | $-p/p_z$ |
| $T_2 = h_0h_1/p_0^2$ | 0 | 0 | 0 | $h_0\,p/h_1 = v/v_0$ |

($T_2$ の $p$ 微分は、$\partial h_1/\partial p = p_0^2p/h_1$ から $\partial T_2/\partial p = h_0p/h_1$。これは速度の比 $v/v_0$ そのもの: $v = c\,p_0p\,/h_1$、$v_0 = c\,p_0/h_0$ から $v/v_0 = p\,h_0/h_1$。[drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (a) の sympy で、$\partial T_2/\partial p - p\,h_0/h_1 = 0$。) 表の合計から

$$
x' = \frac{p_x}{p_z},\qquad y' = \frac{p_y}{p_z},\qquad
z' = -\frac{p}{p_z} + \frac{v}{v_0},\qquad
p_x' = p_y' = p' = 0
\tag{2.4}
$$

これが、運動方程式の全部。**$x' = p_x/p_z$ は、粒子の速度の $s$ 成分に対する $x$ 成分の比 (傾き)** で、$p_x/p$ ではなく $p_x/p_z$ (近軸近似では $p_z \to p$)。

## 5. $z$ の整理と $dv$

(2.4) の $z' = -p/p_z + v/v_0$ に、$dv = 1 - v/v_0$ (01 章の式 (1.9)) を入れる。$-p/p_z + 1$ の桁落ち (傾きが小さいと $p/p_z \to 1$) を避けるために、有理化する: $p - p_z = (p^2 - p_z^2)/(p + p_z) = (p_x^2 + p_y^2)/(p + p_z)$ なので

$$
z' = -\frac{p - p_z}{p_z} - dv = -\frac{p_x^2 + p_y^2}{(p + p_z)\,p_z} - dv
\tag{2.5}
$$

(sympy で、(2.4) の $z'$ と (2.5) の差が 0: 記号 0、数値 $1.4 \times 10^{-16}$。[drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (a)。)

- **第 1 項**は、傾いて進むと、設計軌道の $s$ 軸に沿った長さより長い道を通ること。小さい傾き $q = p_x/p$ で $-(q^2 + r^2)/2$ ($r = p_y/p$)。
- **第 2 項 $-dv$** は、速度が設計速度と違う分の遅れ ($dv > 0$ なら遅い粒子で $z$ は減る)。$dv$ は $\delta$ だけで決まる (式 (1.9))。$\delta \to 0$ で $dv \simeq -\delta/\gamma_0^2$。例: 3 GeV 電子で $\delta = 10^{-3}$ のとき $dv = -2.8970 \times 10^{-11}$ (近似 $-\delta/\gamma_0^2$ は $-2.9013 \times 10^{-11}$。差は $\delta$ の 2 次の項)。閉じた式 (1.9) と $1 - v/v_0$ の直接の式の差は $10^{-16}$ 以下 (同 (a))。

## 6. 厳密解と追跡の式

(2.4) の右辺は $s$ によらない定数 ($p_x, p_y, p$ が一定。$dv$ は $\delta$ の関数なので一定)。長さ $L$ ぶん積分して

$$
\begin{aligned}
x &\to x + \frac{p_x}{p_z}L, & y &\to y + \frac{p_y}{p_z}L,\\
z &\to z - \left[\frac{p_x^2 + p_y^2}{(p + p_z)\,p_z} + dv\right]L, & p_x, p_y, \delta &\to \text{不変}
\end{aligned}
\tag{2.6}
$$

**規格化した運動量 $q = p_x/p$、$r = p_y/p$、$\zeta \equiv p_z/p = \sqrt{1 - q^2 - r^2}$** で書くと、$p_x/p_z = q/\zeta$、$(p_x^2 + p_y^2)/((p + p_z)p_z) = (q^2 + r^2)/((1 + \zeta)\zeta) = (1 - \zeta)/\zeta$ ($q^2 + r^2 = 1 - \zeta^2$ から):

$$
x \to x + \frac{q}{\zeta}L,\qquad y \to y + \frac{r}{\zeta}L,\qquad
z \to z - \left[\frac{1 - \zeta}{\zeta} + dv\right]L
\tag{2.7}
$$

**K64 の追跡 (`tdrift_free`、`tdrift.f:10-23`)** は、これを次のように書く。`dpz = pxy2dpz(px, py)` $= \zeta - 1$ (追跡の `px`、`py` は $q$、$r$。`pxy2dpz` は `tfloor.f:491`。小さい引数 ($q^2 + r^2 < 8 \times 10^{-6}$) では級数 $-x(\tfrac12 + \tfrac18x + \tfrac1{16}x^2)$ ($x = q^2 + r^2$)、そうでなければ $-x/(1 + \sqrt{1 - x})$。$\sqrt{1 - x} - 1$ の桁落ちを避ける形。導出は下の段落)、`al1 = al/(1+dpz)` $= L/\zeta$:

```
x(i) = x(i) + px(i)*al1            ! x += q L/ζ
y(i) = y(i) + py(i)*al1
z(i) = z(i) + dpz*al1 - dv(i)*al   ! z += (ζ-1)/ζ L - dv L
```

(2.7) の $-(1 - \zeta)/\zeta = (\zeta - 1)/\zeta = $ `dpz*al1/al` なので一致。pysad の `track6._tdrift_free` (`track6.py:407`) と numba 版 `track6n.drift1` (`track6n.py:317`) は、これをそのまま書いてある。まとめて

$$
x \to x + q\,\mathrm{al1},\quad y \to y + r\,\mathrm{al1},\quad z \to z + \mathrm{dpz}\cdot\mathrm{al1} - dv\,L
\tag{2.8}
$$

(`al1` $= L/\zeta$、`dpz` $= \zeta - 1$。注: $q, r$ は変数 `px`、`py`。)

**`pxy2dpz` ($= \sqrt{1 - x} - 1$、$x = q^2 + r^2$) の 2 通りの形の導出** (`tfloor.f:491-502`。[drift_06_steps.out.txt](outputs/drift_06_steps.out.txt) (b) で確認)。

- **小さい $x$ の級数**: 二項級数 $\sqrt{1 - x} = \sum_k\binom{1/2}{k}(-x)^k$。$\binom{1/2}{1} = \tfrac12$、$\binom{1/2}{2} = \tfrac{(1/2)(-1/2)}{2!} = -\tfrac18$、$\binom{1/2}{3} = \tfrac{(1/2)(-1/2)(-3/2)}{3!} = \tfrac1{16}$、$\binom{1/2}{4} = -\tfrac5{128}$ なので、$(-x)^k$ の符号を入れて $\sqrt{1 - x} - 1 = -\tfrac x2 - \tfrac{x^2}8 - \tfrac{x^3}{16} - \tfrac{5x^4}{128} - \cdots = -x\bigl(\tfrac12 + \tfrac18x + \tfrac1{16}x^2\bigr) - \tfrac{5x^4}{128} - \cdots$。K64 は $x^3$ までの 3 項を取り、$x^4$ 以降を捨てる。捨てた項の相対的な大きさは $(5x^4/128)/(x/2) = 5x^3/64$ で、閾値 $x = 8\times10^{-6}$ では $4\times10^{-17}$ (絶対値で $1.6\times10^{-22}$)。倍精度の丸め ($1.1\times10^{-16}$) より小さい。
- **有理形**: 分子分母に $\sqrt{1 - x} + 1$ をかけて $\sqrt{1 - x} - 1 = \dfrac{(1 - x) - 1}{\sqrt{1 - x} + 1} = -\dfrac{x}{1 + \sqrt{1 - x}}$。引き算 $\sqrt{1 - x} - 1$ を直接行わないので、小さい $x$ でも桁落ちしない。直接の式は $x = 8\times10^{-6}$ で相対誤差が $1.1\times10^{-16}/(8\times10^{-6}) \approx 1.4\times10^{-11}$ 程度まで悪化する。
- **閾値 $8\times10^{-6}$ の根拠**: 級数の打ち切りの相対誤差 $5x^3/64$ が倍精度の丸めより小さくなる範囲の上限付近 ($5x^3/64 \sim 10^{-16}$ で $x \sim 10^{-5}$)。ソースには値だけが書いてあり、選び方の説明はない (推定)。

**確認** ([drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (b)): 長さ 1 m、3 GeV で、傾きと $\delta$ の違う 4 通り (最大 $q = 0.3$) の粒子について、(i) 閉じた式 (2.6)、(ii) 正準方程式 (2.4) を DOP853 ($rtol = 10^{-13}$) で積分した結果、(iii) pysad の追跡 (DRIFT 1 台) を比べた。(i) と (ii) の差は $1.8 \times 10^{-16}$ 以下、(i) と (iii) の差は $5.6 \times 10^{-17}$ 以下 (4 通りの差は $6.9 \times 10^{-18}$、0、0、$5.6 \times 10^{-17}$)。

## 7. 近軸のドリフトと運動学的ドリフト、DISKIN

### 7.1 厳密 = 近軸 × $D$

四極の章 ([04_quad.md](04_quad.md) §3.2〜§5) は、$H = H_2 + \Delta H$ に分けた。$K = 0$ のとき、

$$
H_2 = -p + \frac{p_x^2 + p_y^2}{2p} + \frac{E}{P_0v_0},\qquad
\Delta H = p - \sqrt{p^2 - p_x^2 - p_y^2} - \frac{p_x^2 + p_y^2}{2p}
\tag{2.9}
$$

で、和は (2.1) (定数 $-p + p = 0$ が相殺) に等しい。**どちらも $p_x, p_y, p$ だけの関数**なので、$x, y, z$ が動いても運動量は一定で、流れは (2.4) と同じ形の一定の速度の移動:

- $H_2$ の流れ (近軸のドリフト): $x \to x + qL$、$y \to y + rL$、$z \to z - [(q^2 + r^2)/2 + dv]L$。($x' = \partial H_2/\partial p_x = p_x/p = q$。$z' = \partial H_2/\partial p = -1 + v/v_0 - (p_x^2 + p_y^2)/(2p^2) = -dv - (q^2 + r^2)/2$。四極の章の式 (4.7) で $K = 0$。)
- $\Delta H$ の流れ ($D$): $x \to x + q(1/\zeta - 1)L$、$z \to z - (2 + \zeta)(1 - \zeta)^2/(2\zeta)\,L$ (四極の章の式 (4.16))。導出: $\Delta H = p - p_z - (p_x^2 + p_y^2)/(2p)$ を偏微分する ((2.3) を使う)。
  - $x' = \partial\Delta H/\partial p_x = p_x/p_z - p_x/p = q(1/\zeta - 1)$ ($p_x/p_z = q/\zeta$、$p_x/p = q$)。
  - $z' = \partial\Delta H/\partial p = 1 - p/p_z + (p_x^2 + p_y^2)/(2p^2) = 1 - 1/\zeta + (1 - \zeta^2)/2$ ($(p_x^2 + p_y^2)/p^2 = q^2 + r^2 = 1 - \zeta^2$)。
  - 通分: $[\,2\zeta - 2 + \zeta - \zeta^3\,]/(2\zeta) = -(\zeta^3 - 3\zeta + 2)/(2\zeta)$。$\zeta^3 - 3\zeta + 2 = (\zeta - 1)^2(\zeta + 2)$ (展開して確認できる: $(\zeta^2 - 2\zeta + 1)(\zeta + 2) = \zeta^3 - 3\zeta + 2$) なので $z' = -(2 + \zeta)(1 - \zeta)^2/(2\zeta)$。
  - 小さい傾きの展開 ($u = q^2 + r^2$、$\zeta = \sqrt{1 - u}$): $1/\zeta = (1 - u)^{-1/2} = \sum_k\binom{-1/2}{k}(-u)^k$ の係数は $1,\ \tfrac12,\ \tfrac38,\ \tfrac5{16},\ \tfrac{35}{128}$ ($\tfrac{(1/2)(3/2)}{2!} = \tfrac38$、$\tfrac{(1/2)(3/2)(5/2)}{3!} = \tfrac5{16}$)。よって $1/\zeta = 1 + \tfrac u2 + \tfrac{3u^2}8 + \tfrac{5u^3}{16} + O(u^4)$、$1/\zeta - 1 = \tfrac u2 + \tfrac{3u^2}8 + \cdots$。$z'$ は、上の $z' = 1 - 1/\zeta + (1 - \zeta^2)/2$ に $(1 - \zeta^2)/2 = u/2$ を入れて $z' = 1 - 1/\zeta + \tfrac u2 = -\tfrac{3u^2}8 - \tfrac{5u^3}{16} - O(u^4)$ (1 次の項は $-\tfrac u2 + \tfrac u2$ で相殺し 2 次から残る。捨てた次の項は $-\tfrac{35}{128}u^4$)。$x$ の進みは $q(1/\zeta - 1) = q(\tfrac u2 + \tfrac{3u^2}8 + \cdots)$。$x$ は $u$ の 1 次から、$z$ は 2 次から効く ([drift_05_expansions.out.txt](outputs/drift_05_expansions.out.txt) (a)(b)、[drift_06_steps.out.txt](outputs/drift_06_steps.out.txt) (c))。

2 つの流れは、運動量だけの関数の流れなので**交換する** (ポアソン括弧 $[H_2, \Delta H] = 0$: どちらも座標を含まず、括弧 $\sum(\partial_q f\,\partial_p g - \partial_p f\,\partial_q g)$ の各項が 0)。したがって、順序によらず、

$$
\text{厳密なドリフト}(L) = D(L) \circ \text{近軸のドリフト}(L) = \text{近軸のドリフト}(L) \circ D(L)
\tag{2.10}
$$

$x$ の和: $q L + q(1/\zeta - 1)L = qL/\zeta$ ✓。$z$ の和: $-(q^2 + r^2)L/2 - (2 + \zeta)(1 - \zeta)^2L/(2\zeta)$ が $-(1 - \zeta)L/\zeta$ に等しいことは、$q^2 + r^2 = 1 - \zeta^2 = (1 - \zeta)(1 + \zeta)$ を使って $-(1 - \zeta)\bigl[(1 + \zeta)\zeta + (2 + \zeta)(1 - \zeta)\bigr]/(2\zeta) = -(1 - \zeta)\cdot 2/(2\zeta)$ ($(1 + \zeta)\zeta + (2 + \zeta)(1 - \zeta) = \zeta + \zeta^2 + 2 - \zeta - \zeta^2 = 2$)。数値の確認は [drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (c): 近軸 → $D$ と厳密な式の差、$D$ → 近軸 との差は $1.4 \times 10^{-17}$ 以下。$D$ の寄与は、傾き $q = 0.3$、$r = 0.2$ の 1 m で $\Delta x = 2.2 \times 10^{-2}$ m、$\Delta z = -7.1 \times 10^{-3}$ m。小さい傾きで、上の展開 $1/\zeta - 1 \simeq (q^2 + r^2)/2$ から $\Delta x \simeq q(q^2 + r^2)L/2$ ($r = 0$ なら $q^3L/2$。四極の章の (4.17))。

したがって、**四極の本体の「$D$ を挟む」分割**は、$K = 0$ では厳密なドリフトに一致し、$K \neq 0$ では厳密なドリフトの代わりに、近軸の流れ ($H_2$ の厚肉の線形の流れ) と運動学的なずれ $D$ を別々に扱う構成になっている。

### 7.2 DISKIN

K64 の追跡 (`tturn.f:287-301`) は、`DISKIN` $\neq 0$ の DRIFT で、$D$ を落として近軸の式だけを使う:

```
x = x + px*al ; y = y + py*al ; z = z - ((px**2 + py**2)*0.5 + dv)*al   ! (px, py は規格化した運動量 q, r)
```

つまり $H_2$ ($K = 0$) の流れ。四極 (`tsolqu` の引数 `kin`) も `DISKIN` $\neq 0$ で $D$ を落とす (`tturn.f:368`、`tturne.f:576`、`qtwiss.f:309` が `cmp%value(ky_KIN_QUAD) == 0.d0` を渡す)。DRIFT の CALC (`qtwiss.f:205`、`tdrife.f`) は、`DISKIN` を見ない (常に厳密)。(`ky_KIN_DRFT`・`ky_KIN_QUAD` は、DRIFT・QUAD の `DISKIN` の値を指す要素の属性の名前、`kin` は四極の追跡ルーチンの論理引数。§15 の表。)

**pysad は `DISKIN` を読むが (`sadkw.py`)、追跡に使っていない**。`DRIFT` に `DISKIN=1` を付けた 1 m のドリフトで、$q = 0.3$ の粒子の追跡は、`DISKIN` なしと同じ (差 0)、厳密な式との差 $5.6 \times 10^{-17}$、近軸の式との差 $2.2 \times 10^{-2}$ ([drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (d))。`DISKIN` を使う台本では、傾きが大きい粒子 (数 mrad 以上。差は $q^3L/2$) で K64 と pysad が違う。傾きが小さければ ($q = 10^{-3}$ で $5 \times 10^{-10}$ m) 見えない。四極の `DISKIN` も pysad に無い (grep で、`DISKIN` を使う追跡のコードが無いことを確認。四極の追跡との比較はしていない)。付録 A に疑いではなく**未実装**として書く。

## 8. COD のまわりの 6×6 行列

CALC・Emittance[] の行列は、閉軌道 (COD) $\xi = (x, p_x, y, p_y, z, \delta)$ を通る粒子の写像 $\Phi$ のヤコビアン $M = \partial\Phi/\partial\xi$。ドリフトの写像 (2.6) を、$p_x, p_y, \delta$ (正準運動量) で偏微分する。$\partial\xi'/\partial\xi$ のうち、$p_x, p_y, \delta$ の行は恒等なので、$x, y, z$ の行だけ。

### 8.1 偏微分を 1 項ずつ

(2.3) を使う。$\partial(1/p_z)/\partial p_x = -p_z^{-2}\cdot(-p_x/p_z) = p_x/p_z^3$、$\partial(1/p_z)/\partial p_y = p_y/p_z^3$、$\partial(1/p_z)/\partial p = -p_z^{-2}\cdot p/p_z = -p/p_z^3$。

**$x$ の行** ($x_f = x + p_xL/p_z$):

$$
\begin{aligned}
\frac{\partial x_f}{\partial x} &= 1, &
\frac{\partial x_f}{\partial p_x} &= \frac{L}{p_z} + p_xL\cdot\frac{p_x}{p_z^3} = \frac{L}{p_z} + \frac{p_x^2L}{p_z^3},\\
\frac{\partial x_f}{\partial p_y} &= p_xL\cdot\frac{p_y}{p_z^3} = \frac{p_xp_yL}{p_z^3}, &
\frac{\partial x_f}{\partial\delta} &= p_xL\cdot\left(-\frac{p}{p_z^3}\right) = -\frac{p_x\,p\,L}{p_z^3}
\end{aligned}
\tag{2.11}
$$

($\partial/\partial\delta = \partial/\partial p$。$p = 1 + \delta$。) **$y$ の行**は $x \leftrightarrow y$ の入れ替え:
$\partial y_f/\partial p_y = L/p_z + p_y^2L/p_z^3$、$\partial y_f/\partial p_x = p_xp_yL/p_z^3$、$\partial y_f/\partial\delta = -p_y\,p\,L/p_z^3$。

**$z$ の行** ($z_f = z - L\,[\,(p/p_z) - 1\,] - dv\,L$。(2.4) の形 $z' = -p/p_z + 1 - dv$ を使う。$-p/p_z$ の微分):

$$
\begin{aligned}
\frac{\partial z_f}{\partial p_x} &= -L\,p\cdot\frac{p_x}{p_z^3} = -\frac{p_x\,p\,L}{p_z^3}, &
\frac{\partial z_f}{\partial p_y} &= -\frac{p_y\,p\,L}{p_z^3},\\
\frac{\partial z_f}{\partial\delta} &= -L\left(\frac{1}{p_z} - \frac{p^2}{p_z^3}\right) + L\,\frac{h_0}{h_1^3} = \frac{p_x^2 + p_y^2}{p_z^3}L + \frac{h_0}{h_1^3}L
\end{aligned}
\tag{2.12}
$$

($-L\,\partial(p/p_z)/\partial p = -L(1/p_z - p\cdot p/p_z^3) = -L(p_z^2 - p^2)/p_z^3 = +L(p_x^2 + p_y^2)/p_z^3$。$-L\,\partial dv/\partial\delta = +L\,\partial(v/v_0)/\partial\delta = L\,h_0/h_1^3$。後者は、$v/v_0 = p\,h_0/h_1$ の $p$ 微分: $h_0/h_1 - p\,h_0\,(p_0^2p/h_1)/h_1^2 = (h_0/h_1^3)(h_1^2 - p_0^2p^2) = h_0/h_1^3$ ($h_1^2 - p_0^2p^2 = 1$)。) $\partial z_f/\partial x = \partial z_f/\partial y = 0$、$\partial z_f/\partial z = 1$。

sympy で、(2.6) を $(x, p_x, y, p_y, z, \delta)$ で偏微分した 6×6 と、次の `tdrife1` の式との差は、記号で 0、数値で $6.4 \times 10^{-18}$ ([drift_02_linearization.out.txt](outputs/drift_02_linearization.out.txt) (a))。

### 8.2 K64 の `tdrife1`

`toplvl.f:1368` の `tdrife1(trans1, cod, dv, dvdp, al)` は、行列 `trans1` (単位行列から出発) と COD の更新を、次の名前で作る。`pr` $= p$、`pxi`、`pyi` $= p_x, p_y$ (COD の正準運動量)、`a` $= p_x^2 + p_y^2$、`pzi` $= p_z$ (`pr*sqrt(1-a/pr**2)`)、`ale` $= L/p_z$、`alz` $=$ `ale/pzi**2` $= L/p_z^3$ ($z$ の式に出る $1/p_z^3$ をまとめた名前):

$$
\begin{aligned}
&T_{12} = \mathrm{ale} + p_x^2\,\mathrm{alz},\quad T_{14} = p_xp_y\,\mathrm{alz},\quad T_{16} = -p_x\,p\,\mathrm{alz},\\
&T_{32} = T_{14},\quad T_{34} = \mathrm{ale} + p_y^2\,\mathrm{alz},\quad T_{36} = -p_y\,p\,\mathrm{alz},\\
&T_{52} = T_{16},\quad T_{54} = T_{36},\quad T_{56} = \mathrm{dvdp}\cdot L + a\,\mathrm{alz}
\end{aligned}
\tag{2.13}
$$

(式 (2.11)(2.12) の各成分に一致。$T_{56}$ の `dvdp` は $h_0/h_1^3$ (`tgetdv`、`tturne.f:880`。$dv$ の $\delta$ 微分の**符号を反転したもの**。上の (2.12) の最後の項)。) 添字は (行, 列) = ($x$, $p_x$, $y$, $p_y$, $z$, $\delta$) の 1〜6。COD も同じルーチンで進める:

$$
x \mathrel{+}= \frac{p_x}{p_z}L,\quad y \mathrel{+}= \frac{p_y}{p_z}L,\quad
z \mathrel{-}= \left(\frac{a}{(p + p_z)\,p_z} + dv\right)L
\tag{2.14}
$$

(`cod(5) = cod(5) - (a/(pr+pzi)/pzi + dv)*al`。(2.6) と同じ式。`tdrife1` の `dv` は引数で、`tgetdv` が返す値、または `tchge` から呼ぶときは 1 (座標の回転・ずれの幾何用)。) CALC の本体 (`tdrife.f`) は、`tdrife1` の結果を、すでにある行列 `trans` に掛け足す:

```
trans(1:5:2, 1:irad) = trans(1:5:2, 1:irad) + matmul(trans1(1:5:2, 2:6:2), trans(2:6:2, 1:irad))
```

$M = 1 + N$ で、$N$ は行 (1, 3, 5) と列 (2, 4, 6) にだけ成分を持つ。$T_{\rm new} = M T = T + N T$ の行 (1, 3, 5) の更新が、この 1 行 (`1:5:2` は 1, 3, 5。`2:6:2` は 2, 4, 6。`1:irad` の `irad` は、更新する列の数で、§15 の表)。

**$\delta$ の列の意味**: 正準運動量 $p_x$ が $p_x = q\,p$ なので、$\delta$ が違うと同じ $q$ でも $x$ の動き方が違う。$T_{16} = -p_x\,p\,L/p_z^3$ は、$p_x$ を固定して $\delta$ だけを動かしたときの $x$ の変化 (傾き $x' = p_x/p_z$ が、$p_z$ が大きいほど小さい)。設計軌道上 ($p_x = 0$) では 0。

## 9. シンプレクティック性と $z$ の行

(2.6) は、ハミルトニアン (2.1) の流れ (厳密な正準変換) なので、ヤコビアン $M$ は $M^\top J M = J$ ($J$ は $\begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$ の 3 つの直和) を満たす。$J$ の条件のうち、$z$ に関する成分から、次の関係が出る。**$z$ の行は、$x, y$ の行から組める**。導出を書く。$M^\top J M = J$ の $(j, 6)$ 成分 ($j = 1, \ldots, 4$。$J_{j6} = 0$) を、$J$ の非零成分 ($J_{12} = J_{34} = J_{56} = 1$、$J_{21} = J_{43} = J_{65} = -1$) を使って書き下す。$(M^\top J M)_{j6} = \sum_{a,b} M_{aj}J_{ab}M_{b6}$ で、$(a, b) = (1,2), (2,1), (3,4), (4,3), (5,6), (6,5)$ だけが残るから

$$
M_{1j}M_{26} - M_{2j}M_{16} + M_{3j}M_{46} - M_{4j}M_{36} + M_{5j}M_{66} - M_{6j}M_{56} = 0
$$

ここで、**$H$ が $z$ を含まないとき** ($z$ と共役な $p$ が保存するので、$z$ が出口の他の座標に影響しない) は $M_{66} = 1$、$M_{6j} = 0$ ($j = 1, \ldots, 5$)。$M_{66}$ と $M_{6j}$ の項を右辺に移して

$$
T_{5j} = -T_{1j}T_{26} + T_{2j}T_{16} - T_{3j}T_{46} + T_{4j}T_{36}\qquad (j = 1,\ldots,4)
\tag{2.15}
$$

(四極の章 (4.38) と同じ条件。**前提は $T_{66} = 1$、$T_{6j} = 0$**。一般のシンプレクティック行列では、上の一般形 $T_{5j}T_{66} - T_{6j}T_{56} = \ldots$ のままで、(2.15) は成り立たない: ランダムなシンプレクティック行列 ($T_{66} = 0.729$) では $j = 1$ の右辺 0.4064 に対し $T_{51} = 0.6465$ ([drift_05_expansions.out.txt](outputs/drift_05_expansions.out.txt) (d)。一般形の左辺は右辺と一致。(d) では上の $(j, 6)$ 成分の式も sympy で書き下した)。) ドリフトでは $T_{2j}$、$T_{4j}$ が単位行列の行 (第 2 行・第 4 行の対角が 1 でほかは 0。$T_{26} = T_{46} = 0$ も含む) なので、(2.15) の右辺は $T_{2j}T_{16} + T_{4j}T_{36}$ だけが残り、$j = 2$: $T_{52} = T_{16}$、$j = 4$: $T_{54} = T_{36}$、$j = 1, 3$: 0。これは (2.13) の $T_{52} = T_{16}$、$T_{54} = T_{36}$ と一致する。

数値の確認 ([drift_02_linearization.out.txt](outputs/drift_02_linearization.out.txt) (c)): COD が傾いた 2 通り (最大 $p_x = 0.3$、$\delta = 0.01$) で、$|M^\top JM - J|$ は丸め誤差以下 (0.0)、$z$ の行と (2.15) で組んだ値の差は 0。**$T_{16} = T_{52}$ という対称性**は、$x$ に対する $\delta$ の効きと、$z$ に対する $p_x$ の効きが同じ大きさであること (正準共役の関係 $x \leftrightarrow p_x$、$z \leftrightarrow \delta$)。

## 10. 最短の例: ドリフトだけの運動量依存の遅れ

設計軌道上の粒子 ($x = p_x = y = p_y = 0$) では (2.6) の $z$ は $z \to z - dv(\delta)\,L$、行列は $T_{12} = T_{34} = L/p$ (COD の $\delta$ が 0 なら $L$)、$T_{56} = L\,h_0/h_1^3$:

$$
R_{56} = L\frac{h_0}{h_1^3}\;\xrightarrow{\ \delta = 0\ }\;\frac{L}{\gamma_0^2}
\tag{2.16}
$$

- **$R_{56}$ は正**。$\delta > 0$ (速い) ほど $z$ が増える = 先頭に出る。ドリフトを通るごとに積算される (k 回通れば k 倍)。ヘルプの $z$ は、先頭が $z > 0$。
- 3 GeV の電子と 10 MeV の電子で、1 m のドリフトを 1〜4 回通った後の $z$ ([drift_03_r56_example.out.txt](outputs/drift_03_r56_example.out.txt)):

| 運動量 | $\gamma_0$ | $\delta$ | $dv$ (1.9) | 1 回後の $z$ | 4 回後の $z$ | $L\delta/\gamma_0^2$ (近似) |
|---|---|---|---|---|---|---|
| 3 GeV/$c$ | 5870.85 | $+10^{-3}$ | $-2.897 \times 10^{-11}$ | $+2.897 \times 10^{-11}$ | $+1.159 \times 10^{-10}$ | $+2.901 \times 10^{-11}$ |
| 3 GeV/$c$ | 5870.85 | $+10^{-2}$ | $-2.858 \times 10^{-10}$ | $+2.858 \times 10^{-10}$ | $+1.143 \times 10^{-9}$ | $+2.901 \times 10^{-10}$ |
| 3 GeV/$c$ | 5870.85 | $-10^{-2}$ | $+2.945 \times 10^{-10}$ | $-2.945 \times 10^{-10}$ | $-1.178 \times 10^{-9}$ | $-2.901 \times 10^{-10}$ |
| 10 MeV/$c$ | 19.595 | $+10^{-3}$ | $-2.601 \times 10^{-6}$ | $+2.601 \times 10^{-6}$ | $+1.040 \times 10^{-5}$ | $+2.604 \times 10^{-6}$ |
| 10 MeV/$c$ | 19.595 | $+10^{-2}$ | $-2.566 \times 10^{-5}$ | $+2.566 \times 10^{-5}$ | $+1.026 \times 10^{-4}$ | $+2.604 \times 10^{-5}$ |
| 10 MeV/$c$ | 19.595 | $-10^{-2}$ | $+2.644 \times 10^{-5}$ | $-2.644 \times 10^{-5}$ | $-1.058 \times 10^{-4}$ | $-2.604 \times 10^{-5}$ |

  追跡の $z$ と $k\cdot(-dv\,L)$ の差は $10^{-20}$ 以下 (表の 4 回後の列。丸め誤差)。$\delta$ の符号を変えると $z$ の符号が変わり、大きさは近似 $L\delta/\gamma_0^2$ と ($\delta$ の 1 次の補正の範囲で) 一致する (表の 2 行目と 3 行目、5 行目と 6 行目のように、$\delta = \pm10^{-2}$ で大きさが違う (3 GeV で 2.858 対 2.945、10 MeV で 2.566 対 2.644) のは、$\delta$ の 2 次の項)。
- 行列: 4 つのドリフト (4 m) の $T_{56} = 1.160533 \times 10^{-7}$ (= $4L/\gamma_0^2$)。$\delta = 10^{-3}$ の粒子の $z$ の変化は、行列で $1.1605 \times 10^{-10}$、追跡で $1.1588 \times 10^{-10}$ (差は $\delta$ の 2 次)。
- 傾いた COD ($p_x = 10^{-2}$) では、$T_{56} = L\,h_0/h_1^3 + p_x^2L/p_z^3$ に傾きの項が加わる (3 GeV で $1.0004 \times 10^{-4}$ 対 $\gamma_0^{-2}$ の項 $2.9 \times 10^{-8}$)。軌道の傾きが $10^{-2}$ なら、傾きの項が運動量の項 $L/\gamma_0^2$ より 3 桁以上大きい ([drift_02_linearization.out.txt](outputs/drift_02_linearization.out.txt) (d))。ビームラインの $R_{56}$ (等時性) の議論で、軌道の傾いた直線部では、このドリフトの項を忘れない。

## 11. Fortran の対応表

行番号は参照した K64 のソースのもの。

| 式・内容 | ファイル:行 | サブルーチン | 変数・メモ |
|---|---|---|---|
| 追跡のドリフト (2.8) | `tdrift.f:10-23` | `tdrift_free` | `dpz = pxy2dpz(px,py)` ($\zeta - 1$)、`al1 = al/(1+dpz)`、`x += px*al1`、`z += dpz*al1 - dv*al` |
| $\sqrt{1 - x} - 1$ | `tfloor.f:491` | `pxy2dpz` | 小さい $q^2 + r^2$ ($< 8 \times 10^{-6}$) は級数 |
| 追跡の DRIFT の呼び出し | `tturn.f:287-301` | `tturn` | `ky_KIN_DRFT == 0` なら `tdrift_free`、そうでなければ近軸式。`spac` なら `spdrift_free` |
| 行列 (2.13) と COD の更新 (2.14) | `toplvl.f:1368-1401` | `tdrife1` | `pzi`、`ale`、`alz`、`trans1(1,2)` など、`cod(5)` の更新 |
| CALC のドリフト (掛け足し) | `tdrife.f:6-` | `tdrife` | `tgetdv`、`tdrife1`、`trans(1:5:2,...) += matmul(...)`。$B_z$・ステアリング (`ak0x`、`ak0y`) があるときは `tsoldz` |
| COD が 0 のドリフトの行列 (Twiss 用のインライン) | `qtwiss.f:205-240` | `qtwiss` (`case (icDRFT)`) | `a12 = ale + pxisq*alz`、`a15 = -pr*pxi*alz`、`cod(5) += dpz*ale - dvemit*al` |
| $dv$ と $dvdp$ | `tturne.f:869-882` | `tgetdv` | `dv = -(pr+1)/h1/(h1+pr*h0)*dp + dvfs`、`dvdp = h0/h1**3` |
| 追跡の $dv$ | `tconvm.f:18` | `tconvm` | 01 章 §8 |
| Emittance[] のドリフト | `tturne.f:471-486` | `tturne` (`case (icDRFT)`) | `tdrife0` を呼ぶ (`DISKIN` を見ない)。IBS (`calint`) のときは半分ずつに分けて `tintrb` |
| 他の要素の中のドリフト | `tthine.f:37`、`tquade.f:27`、`tbende.f:579` など | `tdrife0` | `tdrife` から $B_z$・ステアリングを除いたもの (`tdrife.f` の `tdrife0`) |

## 12. pysad の対応表

| 内容 | pysad の関数 (ファイル:行) | 対応する Fortran | 注意 |
|---|---|---|---|
| 追跡のドリフト (python) | `track6._tdrift_free` (`track6.py:407`) | `tdrift_free` | `c[1]`、`c[3]` は規格化した運動量。`_pxy2dpz` は `_sqrt1(-(px²+py²))` |
| 追跡のドリフト (numba) | `track6n.drift1` (`track6n.py:317`)、`k_drift` (`:327`) | 同上 | `S[i,6]` が $dv$ |
| 行列 (2.13)・COD の更新 (CALC) | `lattice._tdrife1_mat` (`lattice.py:2103`)、`_tdrife1_cod` (`:2005`) | `tdrife1` | `dvdp` は引数。`pzi = pr*sqrt(1-a/pr²)` |
| 行列 (Emittance[]、numba) | `elkern._jac_drift` (`elkern.py:236`)、`k_drift` (`:457`)、`emit6d._jac_drift` (`emit6d.py:434`) | `tdrife1` を `tdrife0` から呼ぶ | `k = L/pz³`、`J[4,5] = k(px²+py²) + dvz·L`、`_dvz` $= h_0/h_1^3$ |
| COD が 0 のドリフトの 4×4 (旧 API の線形部) | `lattice._drift_matrix_at_cod` (`lattice.py:1697`) | `qtwiss.f:205-240` | $z$ の行を持たない (x, y の行だけ) |
| 旧 API (要素クラス) | `Drift.transfer_matrix` (`elements.py:103`)、`tracking.track_drift_exact` (`tracking.py:97`) | `tdrift_free` の一部 | `transfer_matrix` は $L/(1+\delta)$ の $x, y$ だけで、$z$ の行を持たない。`track_drift_exact` は $z$ に $dv$ の項が無い。K64 の `tdrift_free`・`tdrife1` とは別実装 (使用箇所は未確認) |
| $dv$ | `emit6d._Phys.dv` (`emit6d.py:193`)、`elkern._dv` (`elkern.py:103`)、`_dvz` (`:112`) | `tgetdv` | `dvdp` = `_dvz` = $h_0/h_1^3$ |
| `DISKIN` | `sadkw.py` の要素のキーワード表のみ | `tturn.f:296` | 追跡で使わない (未実装。§7.2) |

## 13. 確認スクリプトと出力

スクリプトはリポジトリのルートから `py -3.12 physics_guide/scripts/<名前>.py` で動かす (数秒)。出力は `physics_guide/outputs/<名前>.out.txt`。

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `drift_01_h_to_map.py` | (a) sympy: $H$ の偏微分 (2.4)、$z'$ の整理 (2.5)、$dv$ の閉じた式。(b) 閉じた式 (2.6)・DOP853・pysad の追跡の比較。(c) 厳密 = 近軸 × $D$ (2.10)。(d) DISKIN | [drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) |
| `drift_02_linearization.py` | (a) sympy: (2.6) の 6×6 偏微分と `tdrife1` の式 (2.13) の一致。(b) pysad の行列と追跡の数値微分。(c) シンプレクティック性と $z$ の行 (2.15)。(d) 設計軌道と傾いた COD の成分 | [drift_02_linearization.out.txt](outputs/drift_02_linearization.out.txt) |
| `drift_03_r56_example.py` | 最短の例: 3 GeV・10 MeV の電子で、ドリフトを 1〜4 回通った後の $z$。4 つのドリフトの行列 $T_{56}$ | [drift_03_r56_example.out.txt](outputs/drift_03_r56_example.out.txt) |
| `drift_04_help_convention.py` | ヘルプの $\exp(:H:)$ の向き (流れ) の確定と、DRIFT の $E/v_0$ の符号の食い違い (sympy。DRIFT・QUAD・SEXT・BEND の $V$ の 4 つの式を同じ読みで検証) | [drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt) |
| `drift_05_expansions.py` | 書き足した展開の確認 (sympy): (a) $D$ の流れ (§7.1)、(b) 小さい傾きの展開、(c) 近軸 + $D$ = 厳密、(d) (2.15) の導出と、$T_{66} \neq 1$ の行列で成り立たない例 (§9) | [drift_05_expansions.out.txt](outputs/drift_05_expansions.out.txt) |

主な数字: (2.6) と DOP853 の差 $\le 1.8 \times 10^{-16}$、pysad の追跡との差 $\le 5.6 \times 10^{-17}$ (傾き最大 0.3)。(2.10) の分解の差 $\le 1.4 \times 10^{-17}$。sympy の 6×6 と (2.13) の差 $6.4 \times 10^{-18}$。pysad の行列と追跡の数値微分 (中心差分 $10^{-6}$) の差 $\le 3.9 \times 10^{-11}$ (傾き 0.3 の COD で最大)。

## 14. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV/$c$ と 10 MeV/$c$、長さ 1 m、傾きは最大 $q = 0.3$。ドリフトは理想的 (磁場なし、孔径・散乱なし)。加速はない ($P_0$ 一定)。

**この章の範囲外 (別の章)**:
- 磁場のあるドリフト (ソレノイドの `tsoldz`・`tdrift_solenoid`、ステアリング $ak_{0x}$・$ak_{0y}$) は、ソレノイドの章 ([09_solenoid.md](09_solenoid.md))。`tdrife.f` の `bz` ≠ 0 の分岐は読んでいるが、この章では導出していない。
- 放射 (`tradke`、`tradk`)・空間電荷 (`spdrift_free`) は別の章 ([13_radiation.md](13_radiation.md)、[17_collective.md](17_collective.md))。ドリフトの中で光子は出ない ($B = 0$) ので、`tdrife` の `tradke` 呼び出しはステアリングや $B_z$ があるときだけ意味を持つ (読んだ範囲での理解。放射との組合せの確認はしていない)。
- 取り付け誤差 (DX など) は、ドリフトでは効かない (磁場がないので座標変換が相殺する。確認していない)。

**未確認・限界**:
- **DISKIN** は、Fortran の `tturn.f` の分岐 (近軸式) を読んだだけで、gs (K64 の実行) で数値を比べていない。pysad が `DISKIN` を使わないことは確認した (§7.2)。
- **ヘルプの $H$**: 原文 (DRIFT > transformation) を確認した。$E/v_0$ の規格化 ($P_0$ の略) は読み替え、符号は誤植と判断して (2.1) を使った (ヘルプの式どうしの整合による判断で、gs の実行での確認はしていない)。ヘルプにポアソン括弧の定義は無いが、$\exp(:H:)$ の向きは式どうしの整合から読んだ ([C_poisson_lie.md](C_poisson_lie.md) §8.2、[drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt))。
- **`tdrife.f` の `tdrife`** と `tdrife0` の違い (ステアリングの有無) は、引数と分岐を読んだだけ。
- **確認の範囲**: 閉じた式 (2.6) と DOP853 の比較は、同じ $H$ (2.1) から出した運動方程式 (2.4) を積分しており、「$H$ から写像への導出の正しさ」の確認である。$H$ そのもの (ヘルプの式、01 章の導出) の独立な確認ではない。$H$ の妥当性は [01_variables_and_conventions.md](01_variables_and_conventions.md)・[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) の導出による。
- 行列の確認は、pysad の追跡の数値微分と sympy の偏微分で行い、K64 (gs) の `trans` との比較は、この章ではしていない。
- 傾きが非常に大きい粒子 ($p_x^2 + p_y^2 \to p^2$) では $p_z \to 0$ で式が発散する。K64 の `tdrift_free` は `ampmax = 0.9999` ($(p_x/p_0)^2 + (p_y/p_0)^2$ の最大値。`tdrift.f:4`) を定義しているが、この章の範囲では使い方を追っていない。

## 15. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。この表の定義は、実際のコード (Fortran: `SAD-master/src`、pysad: `pysad/`) で確認した。

| Fortran の名前 | 式 | 定義の場所 (Fortran) | pysad での名前 |
|---|---|---|---|
| `pr` | $p = 1 + \delta$ | `tdrife1`、`tdrife` (`pr=1.d0+cod(6)`) | `pr` |
| `pxi`、`pyi` | COD の正準運動量 $p_x$、$p_y$ (`cod(2)`、`cod(4)`) | `toplvl.f:1375-1376` | `pxi`、`pyi` |
| `a` | $a = p_x^2 + p_y^2$ ($p_z^2 = p^2 - a$ に出る。傾きの 2 乗の目安) | `toplvl.f:1377` | `a` |
| `pzi` | $p_z = p\sqrt{1 - a/p^2}$ (`pr*sqrt(1-a/pr**2)`) | `toplvl.f:1378` | `pzi` |
| `ale` | $L/p_z$ (有効な長さ。傾いた粒子が $s$ 方向に進む長さあたりの経路) | `toplvl.f:1379` | `ale` |
| `alz` | $L/p_z^3$ = `ale/pzi**2` (偏微分で $1/p_z^3$ が出るのでまとめた名前) | `toplvl.f:1380` | `alz` |
| `dvdp` | $h_0/h_1^3 = -\partial dv/\partial\delta$ | `tturne.f:880` (`tgetdv`) | `dvdp`、`_dvz` |
| `dpz` (追跡) | $\sqrt{1 - q^2 - r^2} - 1 = \zeta - 1$ (小さい量。桁落ちなし) | `tdrift.f:19` (`pxy2dpz`) | `dpz` (`_pxy2dpz`) |
| `dpz` (行列側) | $p_z - p$ (canonical。$-a/(p(1 + \sqrt{1 - a/p^2}))$) | `qtwiss.f:212` | `dpz` |
| `al1` | $L/\zeta$ (追跡側。$x$ の進み $q\cdot$`al1`) | `tdrift.f:20` | `al1` |
| `dv` | $dv = 1 - v/v_0$ | `tgetdv`、`tconvm.f` | `dv` |
| `dvemit` | Emittance[]・Twiss 側が持つ $dv$。`tesetdv(dp)` が `dvemit = -(p1+p0)/h1emit/(p0*h1emit+p1*h0)*dp + dvfs` ($p_1 = p_0(1+\delta)$、`h1emit` $= h_1$) を作る (`tgetdv` と同じ式)。`qtwiss.f` の COD の $z$ の更新 (`cod(5)+dpz*ale-dvemit*al`)、`tbende.f` も使う | `tturne.f:847-856` (`tesetdv`) | `_Phys.dv()`・`elkern._dv` |
| `trans1` | 1 つの要素の 6×6 行列 (単位行列から出発) | `toplvl.f:1368` | `t`、`J` |
| `irad` | `tdrife` の引数。行列 `trans` の、更新する列の数。`qtwiss.f:83` と `INTGRL.f:355` では 6 (T の 6 列)。`photons.f:370` の呼び出しでは `ir0`。6 を超えるとき (放射の微分の列を含むと読める) の列の用途は未確認 | `tdrife.f:6`、`:14`、`:34-35` | 未確認 (pysad の行列は 6×6) |
| `kin` | 四極の追跡ルーチン `tquad` の論理引数。偽なら $\Delta H$ ($D$) を落とした近軸の式 (`DISKIN` ≠ 0 に当たる)。`tturn.f:368` が `cmp%value(ky_KIN_QUAD) == 0.d0` を渡す | `tturn.f:368` | 引数を渡す側まで確認 |
| `ky_KIN_DRFT`、`ky_KIN_QUAD` | DRIFT・QUAD の `DISKIN` の値 (`cmp%value(...)`) を指す要素の属性の識別子。0 なら厳密な式 | `tturn.f:293`、`:368` | 使い方まで確認 |
| `bz` | $B_z$ (ソレノイドの縦磁場。`tdrife` の引数)。0 でなければ `tsoldz` の分岐 | `tdrife.f:6` | 引数まで確認。分岐は 09 章 |
| `spac` | 空間電荷を入れる追跡かの論理変数。真なら DRIFT は `spdrift_free` | `tturn.f:287-` | 使い方まで確認。どこで真になるかは未確認 |
| `ampmax` | $0.9999$ (`tdrift.f:4`。$(p_x/p_0)^2 + (p_y/p_0)^2$ の最大値と書かれている。使い方は未確認) | `tdrift.f` | 未確認 |

## 付録 A. ヘルプの記述・K64 との違い (ドリフト)

**ヘルプの記述 (疑いを含む)**

| 箇所 | 内容 |
|---|---|
| ヘルプ DRIFT > transformation の $H$ | **原文で確認** (`SADHelp.html`、DRIFT の transformation。TeX 版・平文版で同じ)。$H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ は $L$ をかけた形で、$E/v_0$ は $E/(P_0v_0)$ ($P_0$ の規格化を略した書き方。QUAD の $H_{2n}$ と同じ略記)。**ヘルプの式どうしの食い違い (式どうしの整合から誤植と判断。gs 未確認)**: $E/v_0$ の符号が QUAD・dynamics・BEND の式と逆。これを誤植と判断したのは、式どうしの整合による (gs での確認はしていない)。$\exp(:H:)$ を「$H$ の時間 1 の流れ」と読むと ($\sqrt{\ }$ の項から $x' = p_x/p_z$ が出る向き)、$\Delta z = \partial H/\partial p = -(p/p_z)L - (v/v_0)L$ になり、正しい $\Delta z = (-p/p_z + v/v_0)L$ (BEND > body の $z_2 = z_1 - s\,p_1/\sqrt{\cdot} + (v_1/v_0)L'$ の直線の極限、QUAD の $H_{2n}$ の $+E/v_0$、dynamics の $H$ の $+E/(P_0v_0)$ と同じ) より $2(v/v_0)L$ ずれる ([drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt) (1))。$+1$ は $H$ を $-p + \ldots$ に近づけるための定数で運動に影響しない。 |
| ヘルプに無いもの | (1) 写像の具体式 (2.6)。(2) $dv$ の式 (1.9) と $z$ の整理 (2.5)。(3) COD のまわりの行列 (2.11)〜(2.13)。(4) 近軸の式との関係 (2.10)。(5) `DISKIN`: **原文 (`SADHelp.html`) に DISKIN は出てこない** (DRIFT のキーワードは `L` と `RADIUS` だけ。要素の一覧表にも無い)。意味は `tturn.f` から |

**K64 の実装の癖・誤りの疑い**

| 箇所 | 内容 |
|---|---|
| `dpz` の二重の意味 | 追跡の `tdrift_free` は $\zeta - 1$ (規格化)、CALC の `qtwiss.f` は $p_z - p$ (canonical)。誤りではなく、**同じ名前で中身が違う**ので注意 (§2)。 |
| `DISKIN` の効き方 | DRIFT の追跡は近軸、DRIFT の CALC は常に厳密。同じ DRIFT でも追跡と CALC の式が違う (誤りとは限らない。`DISKIN` ≠ 0 は近軸にする指定なので、CALC の行列は設計軌道上で一致する。COD が傾くと食い違う。数値では確認していない)。 |

pysad の側の未実装 (K64 と違う動作) は、`DISKIN` (§7.2)。
