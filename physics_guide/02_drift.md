# 02. ドリフト (DRIFT)

**この章で分かること**: 磁場のない区間 (ドリフト) のハミルトニアンから、(1) 厳密な追跡の写像 ($x' = p_x/p_z$、$z$ の更新、$dv$ の項)、(2) 近軸近似 (DISKIN) との関係、(3) 閉軌道 (COD) のまわりで線形化した 6×6 の転送行列 (K64 の `tdrife1`) を、途中の式を省かずに追う。最後に、いちばん短い例として、ドリフトだけで起きる運動量依存の遅れ ($R_{56}$ に当たる量) を数値で示す。ヘルプが書いていない途中 (行列の各成分、$z$ の式の整理、$dv$ の扱い) を埋める。式番号 $(2.n)$ から、§11 の Fortran の対応表と §12 の pysad の対応表を引ける。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) (§1〜§2 と §4 の $dv$、§5 の $q = p_x/(1+\delta)$、§8 の変数名の注釈の規約)。偏微分と行列。ハミルトニアンの導出は [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) にあるが、この章では式 (1.1) の結果だけを使う。

**読む順**: §1 (結論) → §2 (記号) → §3〜§6 (写像の導出) → §7 (近軸との関係) → §8〜§9 (行列) → §10 (例)。コードを読みながら使うときは §11・§12 から入る。四極・偏向磁石・多極は、この章の写像を「刻みの途中の移動」として使う ([04_quad.md](04_quad.md)、[05_bend.md](05_bend.md)、[07_thin_multipole.md](07_thin_multipole.md))。

## 1. 結論

- **ハミルトニアン** (01 章の式 (1.1) で $1/\rho = 0$、$a_x = a_y = a_s = 0$): $H = -\sqrt{p^2 - p_x^2 - p_y^2} + E/(p_0 v_0)$。$H$ は $x, y, z$ を含まないので、$p_x, p_y, p$ は一定で、座標は一定の速さで進む。**厳密に解ける**。
- **厳密な写像** (長さ $L$。$p_z \equiv \sqrt{p^2 - p_x^2 - p_y^2}$):

$$
x \to x + \frac{p_x}{p_z}L,\qquad
y \to y + \frac{p_y}{p_z}L,\qquad
z \to z - \left[\frac{p_x^2 + p_y^2}{(p + p_z)\,p_z} + dv\right]L
$$

  ($p_x, p_y, \delta$ は変わらない。式 (2.6)。$z$ の第 1 項は、傾いて進む分の経路の伸び、第 2 項は速度の差 $dv = 1 - v/v_0$ による遅れ。)
- **設計軌道上の粒子 ($p_x = p_y = 0$) は $\delta$ に比例して $z$ が動く**: $\Delta z = -dv\,L \simeq L\,\delta/\gamma_0^2$ ($\gamma_0$ は設計のローレンツ因子)。電子 3 GeV ($\gamma_0 = 5870.9$) で、1 m あたり $\delta = 10^{-2}$ のとき $+2.9 \times 10^{-10}$ m。10 MeV の電子 ($\gamma_0 = 19.6$) では同じ $\delta$ で $+2.6 \times 10^{-5}$ m ([drift_03_r56_example.out.txt](outputs/drift_03_r56_example.out.txt))。これが転送行列の $R_{56} = \partial z/\partial\delta = L\,h_0/h_1^3$ (式 (2.16)) で、「ドリフトの $R_{56}$」。
- **厳密なドリフトは、近軸のドリフトと、運動学的ドリフト $D$ (四極の章の $\Delta H$ の流れ) の積**で、両者は交換する (式 (2.10))。四極の本体の分割 ([04_quad.md](04_quad.md) §5〜§6) の $D$ は、$K = 0$ のときのこの差に当たる。
- **COD のまわりの 6×6 行列**は、(2.6) を正準変数で偏微分したもの (式 (2.11)〜(2.13))。$x, y, z$ の行だけが単位行列と違い、列は $p_x, p_y, \delta$ だけ。COD に傾き $p_x$ があると、$R_{16} = -p_x\,p\,L/p_z^3 \simeq -p_x L$ が現れる ($\delta$ が違うと $x$ が動く: 軌道の傾きによる経路・位置の運動量依存)。シンプレクティックで、$z$ の行はシンプレクティック条件からも組める (§9)。
- **DISKIN (近軸のドリフト)**: K64 の追跡は、`DISKIN` ≠ 0 のドリフトで近軸の式 ($x \to x + qL$、$z \to z - (q^2/2 + dv)L$) を使う。pysad の追跡は `DISKIN` を読むが使わず、厳密な式のまま (未実装。§7.2)。
- 数値の確認: 閉じた式は、正準方程式の数値積分 (DOP853) と $2 \times 10^{-16}$、pysad の追跡と $6 \times 10^{-17}$ で一致。行列は、sympy の偏微分と $6 \times 10^{-18}$、追跡の数値微分と $4 \times 10^{-11}$ 以下 (§13)。
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
| $p_z - p$ (**行列側**。canonical で、規格化していない) | | `dpz` (`tdrife1` の周辺、`qtwiss.f:219`。名前は同じで中身が違う) | `dpz` |
| 速度の差 $dv = 1 - v/v_0$、$-\partial dv/\partial\delta$ | | `dv`、`dvdp` (`tgetdv`) | `p.dv`、`_dvz` |
| 座標・運動量・$\delta$ | $\{x, p_x, y, p_y, z, \delta\}$ | `x(i)`, `px(i)`, … (追跡)。`cod(1..6)` (CALC) | `c[0..5]` |

- Fortran の `dpz` は 2 つの意味で使われる (**取り違えやすい**)。追跡の `tdrift_free` では $\zeta - 1 = \sqrt{1 - q^2 - r^2} - 1$ ($q = p_x/p$ は規格化した運動量)。CALC の `qtwiss.f` の `dpz` は $p_z - p$ ($p_x, p_y$ が正準運動量)。どちらも $-(\text{傾きの 2 乗})/2$ に近い小さい負の数。
- 添字は、Fortran は 1 始まり (`cod(1)` が $x$)、pysad は 0 始まり (`c[0]` が $x$)。
- 数値の例は電子、運動量 3 GeV/$c$ ($\gamma_0 = 5870.85$)。`MOMENTUM = 3GEV` から。

## 3. ハミルトニアン

01 章の式 (1.1) を、直線上 ($1/\rho = 0$)・磁場なし ($a_x = a_y = a_s = 0$) で書くと

$$
H = -\sqrt{p^2 - p_x^2 - p_y^2} + \frac{E}{p_0 v_0}
\tag{2.1}
$$

($E$ は粒子のエネルギー、$p_0 v_0$ は設計運動量と設計速度の積。$E/(p_0v_0)$ は $p$ だけの関数で、$z$ と共役な $p$ の側の量。) ヘルプの DRIFT の transformation には $H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ ($L$ をかけた形で、$\exp(:H:)$ に入れる) とある。定数 $1$ は運動に影響しない。$E/v_0$ の符号と、$p_0$ の規格化の略は、ヘルプの書き方 (付録 A。ヘルプの $H$ と (2.1) の関係は、[C_poisson_lie.md](C_poisson_lie.md) §8.2 の符号の注意による)。

$E$ を $p$ で書く: 設計運動量を $mc$ の単位で $p_0$ ($= \beta_0\gamma_0$)、$h_0 = \gamma_0 = \sqrt{1 + p_0^2}$、$h_1 = \sqrt{1 + p_0^2p^2}$ とすると、$E = h_1\,mc^2$、$p_0v_0 \to p_0\cdot\beta_0 c\cdot mc = p_0^2 mc^2/h_0$ で

$$
\frac{E}{p_0 v_0} = \frac{h_0\,h_1}{p_0^2}
\tag{2.2}
$$

($E/(p_0v_0)$ は無次元。$h_1$ は $p$ に依存。)

**$H$ が含む変数**: $p_x, p_y, p$ だけ。$x, y, z$ も $s$ も含まない。したがって $x' = \partial H/\partial p_x$ などの右辺は定数で (次の §4)、方程式は積分できる。

## 4. 運動方程式 (偏微分を 1 項ずつ)

正準方程式 ($s$ が独立変数。[01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §8.3): $x' = \partial H/\partial p_x$、$p_x' = -\partial H/\partial x$、同様に $y$、$z' = \partial H/\partial p$、$p' = -\partial H/\partial z$。$p_z = (p^2 - p_x^2 - p_y^2)^{1/2}$ の偏微分は、鎖則で

$$
\frac{\partial p_z}{\partial p_x} = \frac{-p_x}{p_z},\qquad
\frac{\partial p_z}{\partial p_y} = \frac{-p_y}{p_z},\qquad
\frac{\partial p_z}{\partial p} = \frac{p}{p_z}
\tag{2.3}
$$

$H = T_1 + T_2$、$T_1 = -p_z$、$T_2 = E/(p_0v_0)$ として:

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

**K64 の追跡 (`tdrift_free`、`tdrift.f:10-23`)** は、これを次のように書く。`dpz = pxy2dpz(px, py)` $= \zeta - 1$ (追跡の `px`、`py` は $q$、$r$。`pxy2dpz` は `tfloor.f:491`。小さい引数 ($q^2 + r^2 < 8 \times 10^{-6}$) では級数 $-x(\tfrac12 + \tfrac18x + \tfrac1{16}x^2)$ ($x = q^2 + r^2$)、そうでなければ $-x/(1 + \sqrt{1 - x})$。$\sqrt{1 - x} - 1$ の桁落ちを避ける形)、`al1 = al/(1+dpz)` $= L/\zeta$:

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

**確認** ([drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (b)): 長さ 1 m、3 GeV で、傾きと $\delta$ の違う 4 通り (最大 $q = 0.3$) の粒子について、(i) 閉じた式 (2.6)、(ii) 正準方程式 (2.4) を DOP853 ($rtol = 10^{-13}$) で積分した結果、(iii) pysad の追跡 (DRIFT 1 台) を比べた。(i) と (ii) の差は $1.8 \times 10^{-16}$ 以下、(i) と (iii) の差は $6.9 \times 10^{-17}$ 以下。

## 7. 近軸のドリフトと運動学的ドリフト、DISKIN

### 7.1 厳密 = 近軸 × $D$

四極の章 ([04_quad.md](04_quad.md) §3.2〜§5) は、$H = H_2 + \Delta H$ に分けた。$K = 0$ のとき、

$$
H_2 = -p + \frac{p_x^2 + p_y^2}{2p} + \frac{E}{p_0v_0},\qquad
\Delta H = p - \sqrt{p^2 - p_x^2 - p_y^2} - \frac{p_x^2 + p_y^2}{2p}
\tag{2.9}
$$

で、和は (2.1) (定数 $-p + p = 0$ が相殺) に等しい。**どちらも $p_x, p_y, p$ だけの関数**なので、$x, y, z$ が動いても運動量は一定で、流れは (2.4) と同じ形の一定の速度の移動:

- $H_2$ の流れ (近軸のドリフト): $x \to x + qL$、$y \to y + rL$、$z \to z - [(q^2 + r^2)/2 + dv]L$。($x' = \partial H_2/\partial p_x = p_x/p = q$。$z' = \partial H_2/\partial p = -1 + v/v_0 - (p_x^2 + p_y^2)/(2p^2) = -dv - (q^2 + r^2)/2$。四極の章の式 (4.7) で $K = 0$。)
- $\Delta H$ の流れ ($D$): $x \to x + q(1/\zeta - 1)L$、$z \to z - (2 + \zeta)(1 - \zeta)^2/(2\zeta)\,L$ (四極の章の式 (4.16))。

2 つの流れは、運動量だけの関数の流れなので**交換する** (ポアソン括弧 $[H_2, \Delta H] = 0$: どちらも座標を含まず、括弧 $\sum(\partial_q f\,\partial_p g - \partial_p f\,\partial_q g)$ の各項が 0)。したがって、順序によらず、

$$
\text{厳密なドリフト}(L) = D(L) \circ \text{近軸のドリフト}(L) = \text{近軸のドリフト}(L) \circ D(L)
\tag{2.10}
$$

$x$ の和: $q L + q(1/\zeta - 1)L = qL/\zeta$ ✓。$z$ の和: $-(q^2 + r^2)L/2 - (2 + \zeta)(1 - \zeta)^2L/(2\zeta)$ が $-(1 - \zeta)L/\zeta$ に等しいことは、$q^2 + r^2 = 1 - \zeta^2 = (1 - \zeta)(1 + \zeta)$ を使って $-(1 - \zeta)\bigl[(1 + \zeta)\zeta + (2 + \zeta)(1 - \zeta)\bigr]/(2\zeta) = -(1 - \zeta)\cdot 2/(2\zeta)$ ($(1 + \zeta)\zeta + (2 + \zeta)(1 - \zeta) = \zeta + \zeta^2 + 2 - \zeta - \zeta^2 = 2$)。数値の確認は [drift_01_h_to_map.out.txt](outputs/drift_01_h_to_map.out.txt) (c): 近軸 → $D$ と厳密な式の差、$D$ → 近軸 との差は $1.4 \times 10^{-17}$ 以下。$D$ の寄与は、傾き $q = 0.3$、$r = 0.2$ の 1 m で $\Delta x = 2.2 \times 10^{-2}$ m、$\Delta z = -7.1 \times 10^{-3}$ m。小さい傾きで $\Delta x \simeq q^3L/2$ (四極の章の (4.17))。

したがって、**四極の本体の「$D$ を挟む」分割**は、$K = 0$ では厳密なドリフトに一致し、$K \neq 0$ では厳密なドリフトの代わりに、近軸の流れ ($H_2$ の厚肉の線形の流れ) と運動学的なずれ $D$ を別々に扱う構成になっている。

### 7.2 DISKIN

K64 の追跡 (`tturn.f:287-301`) は、`DISKIN` $\neq 0$ の DRIFT で、$D$ を落として近軸の式だけを使う:

```
x = x + px*al ; y = y + py*al ; z = z - ((px**2 + py**2)*0.5 + dv)*al   ! (px, py は規格化した運動量 q, r)
```

つまり $H_2$ ($K = 0$) の流れ。四極 (`tsolqu` の引数 `kin`) も `DISKIN` $\neq 0$ で $D$ を落とす (`tturn.f:368`、`tturne.f:576`、`qtwiss.f:309` が `cmp%value(ky_KIN_QUAD) == 0.d0` を渡す)。DRIFT の CALC (`qtwiss.f:205`、`tdrife.f`) は、`DISKIN` を見ない (常に厳密)。

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

$M = 1 + N$ で、$N$ は行 (1, 3, 5) と列 (2, 4, 6) にだけ成分を持つ。$T_{\rm new} = M T = T + N T$ の行 (1, 3, 5) の更新が、この 1 行 (`1:5:2` は 1, 3, 5。`2:6:2` は 2, 4, 6)。

**$\delta$ の列の意味**: 正準運動量 $p_x$ が $p_x = q\,p$ なので、$\delta$ が違うと同じ $q$ でも $x$ の動き方が違う。$T_{16} = -p_x\,p\,L/p_z^3$ は、$p_x$ を固定して $\delta$ だけを動かしたときの $x$ の変化 (傾き $x' = p_x/p_z$ が、$p_z$ が大きいほど小さい)。設計軌道上 ($p_x = 0$) では 0。

## 9. シンプレクティック性と $z$ の行

(2.6) は、ハミルトニアン (2.1) の流れ (厳密な正準変換) なので、ヤコビアン $M$ は $M^\top J M = J$ ($J$ は $\begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$ の 3 つの直和) を満たす。$J$ の条件のうち、$z$ に関する成分から、次の関係が出る。**$z$ の行は、$x, y$ の行から組める**:

$$
T_{5j} = -T_{1j}T_{26} + T_{2j}T_{16} - T_{3j}T_{46} + T_{4j}T_{36}\qquad (j = 1,\ldots,4)
\tag{2.15}
$$

(四極の章 (4.38) と同じ条件。$J$ の条件のうち $(j, 6)$ 成分。) ドリフトでは $T_{2j}$、$T_{4j}$ が単位行列 (第 2 列・第 4 列の対角が 1 でほかは 0) なので、$j = 2$: $T_{52} = T_{16}$、$j = 4$: $T_{54} = T_{36}$、$j = 1, 3$: 0。これは (2.13) の $T_{52} = T_{16}$、$T_{54} = T_{36}$ と一致する。

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
| 3 GeV/$c$ | 5870.9 | $+10^{-3}$ | $-2.897 \times 10^{-11}$ | $+2.897 \times 10^{-11}$ | $+1.159 \times 10^{-10}$ | $+2.901 \times 10^{-11}$ |
| 3 GeV/$c$ | 5870.9 | $+10^{-2}$ | $-2.858 \times 10^{-10}$ | $+2.858 \times 10^{-10}$ | $+1.143 \times 10^{-9}$ | $+2.901 \times 10^{-10}$ |
| 10 MeV/$c$ | 19.6 | $+10^{-3}$ | $-2.601 \times 10^{-6}$ | $+2.601 \times 10^{-6}$ | $+1.040 \times 10^{-5}$ | $+2.604 \times 10^{-6}$ |
| 10 MeV/$c$ | 19.6 | $-10^{-2}$ | $+2.644 \times 10^{-5}$ | $-2.644 \times 10^{-5}$ | $-1.058 \times 10^{-4}$ | $-2.604 \times 10^{-5}$ |

  追跡の $z$ と $k\cdot(-dv\,L)$ の差は $10^{-20}$ 以下 (表の 4 回後の列。丸め誤差)。$\delta$ の符号を変えると $z$ の符号が変わり、大きさは近似 $L\delta/\gamma_0^2$ と ($\delta$ の 1 次の補正の範囲で) 一致する ($\delta = 10^{-2}$ の非対称 2.858 対 2.945 は、$\delta$ の 2 次の項)。
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

主な数字: (2.6) と DOP853 の差 $\le 1.8 \times 10^{-16}$、pysad の追跡との差 $\le 6.9 \times 10^{-17}$ (傾き最大 0.3)。(2.10) の分解の差 $\le 1.4 \times 10^{-17}$。sympy の 6×6 と (2.13) の差 $6.4 \times 10^{-18}$。pysad の行列と追跡の数値微分 (中心差分 $10^{-6}$) の差 $\le 3.9 \times 10^{-11}$ (傾き 0.3 の COD で最大)。

## 14. 未確認・限界・使った条件

**使った条件**: 数値例は電子 3 GeV/$c$ と 10 MeV/$c$、長さ 1 m、傾きは最大 $q = 0.3$。ドリフトは理想的 (磁場なし、孔径・散乱なし)。加速はない ($p_0$ 一定)。

**この章の範囲外 (別の章)**:
- 磁場のあるドリフト (ソレノイドの `tsoldz`・`tdrift_solenoid`、ステアリング $ak_{0x}$・$ak_{0y}$) は、ソレノイドの章 (予定。`09_solenoid.md`)。`tdrife.f` の `bz` ≠ 0 の分岐は読んでいるが、この章では導出していない。
- 放射 (`tradke`、`tradk`)・空間電荷 (`spdrift_free`) は別の章 (予定。`13_radiation.md`、`17_collective.md`)。ドリフトの中で光子は出ない ($B = 0$) ので、`tdrife` の `tradke` 呼び出しはステアリングや $B_z$ があるときだけ意味を持つ (読んだ範囲での理解。放射との組合せの確認はしていない)。
- 取り付け誤差 (DX など) は、ドリフトでは効かない (磁場がないので座標変換が相殺する。確認していない)。

**未確認・限界**:
- **DISKIN** は、Fortran の `tturn.f` の分岐 (近軸式) を読んだだけで、gs (K64 の実行) で数値を比べていない。pysad が `DISKIN` を使わないことは確認した (§7.2)。
- **ヘルプの $H$** の符号・$E/v_0$ の規格化は、(2.1) と比べて読み替えた。ヘルプのポアソン括弧の定義が書かれていないので、$\exp(:H:)$ の符号の対応は推定 ([C_poisson_lie.md](C_poisson_lie.md) §8.2)。
- **`tdrife.f` の `tdrife`** と `tdrife0` の違い (ステアリングの有無) は、引数と分岐を読んだだけ。
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
| `dpz` (行列側) | $p_z - p$ (canonical。$-a/(p(1 + \sqrt{1 - a/p^2}))$) | `qtwiss.f:219` | `dpz` |
| `al1` | $L/\zeta$ (追跡側。$x$ の進み $q\cdot$`al1`) | `tdrift.f:20` | `al1` |
| `dv` | $dv = 1 - v/v_0$ | `tgetdv`、`tconvm.f` | `dv` |
| `trans1` | 1 つの要素の 6×6 行列 (単位行列から出発) | `toplvl.f:1368` | `t`、`J` |
| `ampmax` | $0.9999$ (`tdrift.f:4`。$(p_x/p_0)^2 + (p_y/p_0)^2$ の最大値と書かれている。使い方は未確認) | `tdrift.f` | 未確認 |

## 付録 A. ヘルプの記述・K64 との違い (ドリフト)

**ヘルプの記述 (疑いを含む)**

| 箇所 | 内容 |
|---|---|
| ヘルプ DRIFT > transformation の $H$ | $H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ は、$L$ をかけた形で、$E/v_0$ は $E/(p_0v_0)$ ($p_0$ の規格化を略した書き方。四極の章の付録 A と同じ)。$\exp(:H:)$ の符号は [C_poisson_lie.md](C_poisson_lie.md) §8.2 の推定による (**疑い**: ヘルプにポアソン括弧の定義が無く、確認していない)。 |
| ヘルプに無いもの | (1) 写像の具体式 (2.6)。(2) $dv$ の式 (1.9) と $z$ の整理 (2.5)。(3) COD のまわりの行列 (2.11)〜(2.13)。(4) 近軸の式との関係 (2.10)。(5) `DISKIN` の意味の説明 (キーワードの名前のみ。内容は `tturn.f` から) |

**K64 の実装の癖・誤りの疑い**

| 箇所 | 内容 |
|---|---|
| `dpz` の二重の意味 | 追跡の `tdrift_free` は $\zeta - 1$ (規格化)、CALC の `qtwiss.f` は $p_z - p$ (canonical)。誤りではなく、**同じ名前で中身が違う**ので注意 (§2)。 |
| `DISKIN` の効き方 | DRIFT の追跡は近軸、DRIFT の CALC は常に厳密。同じ DRIFT でも追跡と CALC の式が違う (誤りとは限らない。`DISKIN` ≠ 0 は近軸にする指定なので、CALC の行列は設計軌道上で一致する。COD が傾くと食い違う。数値では確認していない)。 |

pysad の側の未実装 (K64 と違う動作) は、`DISKIN` (§7.2)。
