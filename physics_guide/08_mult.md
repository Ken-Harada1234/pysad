# 多極 MULT: ポテンシャルの展開、本体の分割、ずれ、ANGLE ≠ 0

**この章で分かること**: MULT は、1 つの要素に K0〜K21 (と SK0〜SK21) の多極を重ねられる要素である。この章は、(1) 多極のポテンシャル $A_s$ の展開と K・SK の意味 (ANGLE = 0 は直線の多極、ANGLE ≠ 0 は曲率つきの展開)、(2) ANGLE = 0 の追跡の構成 (厳密な四極 + 双極の本体 `tsolqu` と、$n \ge 2$ の薄いキック `tmulk` の対称な分割)、刻みの数 `ndiv`、端のフリンジ、(3) ずれ・回転 (DX、DY、DZ、CHI1、CHI2、ROTATE)、(4) ANGLE ≠ 0 の $A_s$ を Maxwell 方程式から導き、ヘルプの $g_{kn}$ との食い違いを確定する、までを、式を省かずに追う。四極の本体は [04_quad.md](04_quad.md)、薄い多極のキックと非線形フリンジは [07_thin_multipole.md](07_thin_multipole.md)、偏向磁石の本体と端は [05_bend.md](05_bend.md)、ソレノイド中の本体は [09_solenoid.md](09_solenoid.md)。RF (VOLT) つきの MULT は [10_cavity.md](10_cavity.md)。

**読む順**: §1 (結論) → §3 (ポテンシャル) → §4 (ANGLE = 0 の構成) → §6 (ANGLE ≠ 0) → §8 (確認)。

## 要約と目次

**この章で分かること**

- **何を導くか**: MULT (K0〜K21・SK0〜SK21 の多極を 1 つの要素に重ねる) の、(1) ポテンシャルの展開と K・SK の意味、(2) ANGLE = 0 の追跡の構成 (厳密な四極 + 双極の本体と、$n \ge 2$ の薄いキックの対称な分割)・刻みの数・端部、(3) ずれ・回転 (DX・DY・DZ・CHI1・CHI2・ROTATE)、(4) ANGLE ≠ 0 の $A_s$ の Maxwell 方程式からの導出。
- **主な結果**: キック (8.3)、SK$n$ は K$n$ の回転 (8.4)、分割 (8.7) と刻みの数 (8.8)・(8.8a)、端面の角の式 (8.9)、ずれの剛体変換 (8.10)、ANGLE ≠ 0 の $A_s$ の方程式 (8.13) と係数の漸化式 (8.14)、展開 (8.15)、キック (8.16)。
- **対応**: 追跡は Fortran の `tmulti` (ANGLE ≠ 0 は `tmulta`)、本体は `tsolqu`、キックは `tmulk`、ずれは `tsolrot`。CALC は `tmulte`・`tmultae`。pysad の対応は付録 B の表 (関数名つき)。
- **確かめ方**: sympy、DOP853 による厳密な運動との比較、K64 のソースの読み、gs の `TrackParticles` との比較 (§5・§6.4・§6.5。`mult_01`〜`mult_06`、§8)。
- **K64・ヘルプの誤りの確定・疑い** (付録 A と一致): K64 の**確定** = B331 (ANGLE ≠ 0 で円弧の長さ $L'$ の補正が 2 回入る。gs と ODE で確認。§6.4)、B340 (K1 ≠ 0 で FB1・FB2 の端のキックが途中にも入る。gs で確認。§6.5)。いずれも v2 は修正。ヘルプは、$g_{kn}$ の $(2k-1)!!$ が $(2k+1)!!$ の誤り、$z$ の冪 $n + k$ が $n + k + 1$ の誤り (どちらも確定。§6.3)。**取り下げ** = CHI1・CHI2 の符号 (gs で確認。誤りではない。§5)。疑い (未決) = (8.9) の $n \ge 1$ の項と pysad v2 の厳密な幾何との違い (どちらが厳密かは未確認)。

**導出の流れ**

1. [§1 結論](#1-結論)・[§2 記号と単位](#2-記号と単位): 結果の要点と、SAD・K64 の Fortran・pysad の記号の対応。
2. [§3 多極のポテンシャル (直線、ANGLE = 0)](#3-多極のポテンシャル-直線angle--0): ポテンシャル (8.1) からキック (8.2)(8.3) を作る。
   - [§3.1 K_n・SK_n のキック](#31-k_nsk_n-のキック)、[§3.2 SK_n は K_n の回転 (ヘルプの記述の確認)](#32-sk_n-は-k_n-の回転-ヘルプの記述の確認) (8.4)、[§3.3 基準面 (θ1)](#33-基準面-θ1)。
3. [§4 ANGLE = 0 の構成](#4-angle--0-の構成): 本体とキックの分割。
   - [§4.1 ハミルトニアンの分割](#41-ハミルトニアンの分割) (8.5)(8.6)、[§4.2 分割と順序 (`tmulti`)](#42-分割と順序-tmulti) (8.7。誤差が 2 次になること)、[§4.3 刻みの数 `ndiv`](#43-刻みの数-ndiv) ((8.8) → (8.8a) → (8.8b)。$\varepsilon$ の読みは推定)、[§4.4 フリンジ](#44-フリンジ) (非線形・線形・端面の角 (8.9) → (8.9a))。
4. [§5 ずれ・回転 (DX、DY、DZ、CHI1、CHI2、ROTATE)](#5-ずれ回転-dxdydzchi1chi2rotate): 剛体変換 (8.10) と `tsolrot` の対応。
   - [§5.1 `tsolrot` の式と (8.10) の対応 (導出)](#51-tsolrot-の式と-810-の対応-導出): 入口・出口が逆変換であること、回転の支点、CHI の符号 (gs で確認)。
5. §6 ANGLE ≠ 0 の $A_s$: 曲線座標の $A_s$ を導く。
   - [§6.1 ヘルプの定義](#61-ヘルプの定義) (8.11)、[§6.2 導出](#62-導出) ((8.12) → (8.13) → (8.14) → (8.15) → (8.16))、[§6.3 ヘルプの式との食い違い (確定)](#63-ヘルプの式との食い違い-確定)。
   - [§6.4 ANGLE ≠ 0 の追跡の構成 (`tmulta`)](#64-angle--0-の追跡の構成-tmulta) ($L'$ の二重補正 = B331)、[§6.5 K1 ≠ 0 の MULT の端 (FB1・FB2) が途中にも入る](#65-k1--0-の-mult-の端-fb1fb2-が途中にも入る-k64-の誤りgs-で確認v2-は修正--b340) (B340)。
6. [§7 行列 (CALC・Emittance[])](#7-行列-calcemittance)・[§8 確認](#8-確認): CALC の構成と、確認スクリプトの一覧。
7. [§9 未確認・限界・使った条件](#9-未確認限界使った条件)・[付録 A](#付録-a-ヘルプの記述k64-との違い-mult)・[付録 B](#付録-b-この章で使う-fortran-の変数名)。

主な式の連鎖: (8.1) → (8.2) → (8.3) → (8.4)、(8.5) → (8.6) → (8.7) → (8.8) → (8.8a) → (8.8b)、(8.9) → (8.9a)、(8.10)、(8.11) → (8.12) → (8.13) → (8.14) → (8.15) → (8.16)。

**読み方**: 先に [04_quad.md](04_quad.md) (四極の本体)、[07_thin_multipole.md](07_thin_multipole.md) (薄い多極のキックと非線形フリンジ)、ANGLE ≠ 0 は [05_bend.md](05_bend.md) を読む。四極・偏向の記号は 01 章 §8 の規約。§4.3 の $\varepsilon$ の読みと §4.4 の (8.9a) は、式の由来を追う節で、飛ばしてよい。§6.2 は $A_s$ の導出で、結果だけなら (8.14)(8.15) で足りる。コードの対応は付録 B から入る。K64 の `tmulti` の引数 $-\mathrm{CHI}$ は内部の規約 (§5.1)。

**未確認・限界の所在**: [§9 未確認・限界・使った条件](#9-未確認限界使った条件) (CHI1・CHI2・DZ、(8.9) と v2 の違い、(8.8) の $\varepsilon$ の読み、ANGLE ≠ 0 の限界、B331・B340 の未確認の範囲)。K64 の誤りの確定・取り下げ・疑いは付録 A の末尾。

## 1. 結論

- **K_n・SK_n の意味**: $\mathrm{K}n$ ($n = 0, 1, 2, \ldots$) は積分強度 $B^{(n)}L/B\rho$ [1/m$^n$] で、横方向のポテンシャル $U_n = \mathrm{Re}\,[(\mathrm{K}n + i\,\mathrm{SK}n)\,z^{n+1}]/(n+1)!$ ($z = x + iy$、長さあたりにするには $L$ で割る) を作る。キックは $\Delta(p_x, p_y) = (-\mathrm{Re}, +\mathrm{Im})\,[(\mathrm{K}n + i\,\mathrm{SK}n)\,z^n/n!]/(1+\delta)$ (追跡の内部変数 $q$ で)。SK$n$ > 0 は、K$n$ の磁石を ROTATE = $90°/(n+1)$ 回したものと**同じ** (ヘルプの記述を数値で確認: 差 $10^{-17}$)。
- **基準面の取り方**: K1 + i SK1 の位相の半分 $\theta_1 = \tfrac12\arg(\mathrm{K1} + i\,\mathrm{SK1})$ だけ座標を回して、スキュー四極の成分を 0 にする (ヘルプの「reference plane is defined so that the skew quadrupole component becomes zero」)。各成分は $(\mathrm{K}n + i\,\mathrm{SK}n)\,e^{-i(n+1)\theta_1}$ に変わる。
- **ANGLE = 0 の追跡** (`tmulti`): 本体は $B(h/2)\,K\,B(h)\,K\cdots K\,B(h/2)$ ($h = L/\mathrm{ndiv}$)。$B$ は K1 と K0 (と SK0) だけの**厳密な線形の本体** (`tsolqu`)、$K$ は $n \ge 2$ の多極の**薄いキック**。対称な 2 次の分割で、誤差は $\mathrm{ndiv}^{-2}$ に比例する (数値確認)。刻み数は $\mathrm{ndiv} = 1 + \bigl\lfloor\max_{n\ge2}\sqrt{a_{\max}^{\,n-1}\,|K_n L|/(6\,(n-1)!\,\varepsilon)}\,\bigr\rfloor$、$a_{\max} = 0.05$ m、$\varepsilon = 0.005\,\mathrm{EPS}$ (EPS = 0 は 1)、上限 1000。
- **ずれ**: DX、DY、DZ、CHI1、CHI2、ROTATE (= CHI3) は剛体の座標変換 (ヘルプの 3×3 の式)。DX・DY・ROTATE の向きは追跡で確認した。CHI1・CHI2・DZ は、K64 の `tsolrot` の式がヘルプの (8.10) の $R_2R_1$ と同じ剛体変換であること、回転の支点が**要素の中心** ($\mathtt{alg} = L/2$) であること、入口と出口が逆変換であることを sympy と数値で確認した (§5)。MULT が `tsolrot` に $-\mathrm{CHI}1$、$-\mathrm{CHI}2$ を渡す理由は、ANGLE ≠ 0 の `tbrot` と向きをそろえる内部の規約 (§5)。**CHI1・CHI2 の符号は gs で確認した** (要素の軸は CHI1 > 0 で $-x$ 側、CHI2 > 0 で $-y$ 側へ向く。(8.10) の $R$ を要素の向きを与える行列と読んだものと同じ。K64 の誤りではない。§5)。
- **端面の角 (AE1・AE2)**: K64 の `tmulk` の式 (8.9) は、中央面では「面までの余分な長さ × 力」に一致し、運動量だけのキックとして調和・勾配 (シンプレクティック) である (§4.4 で導出)。pysad の v2 は別の (面の厳密な幾何の) 式で、$n \ge 1$ の水平方向の項が違う (数値で比較)。
- **ANGLE ≠ 0**: 曲率 $1/\rho$ ($\rho = L/\mathrm{ANGLE}$) の座標で、$A_s$ は Maxwell 方程式から $\nabla^2F = F_x/(\rho + x)$ ($F \propto A_s$) を満たし、$A_s = \sum_{k,n}g_{kn}\,\mathrm{Re}\,[(\mathrm{K}n + i\,\mathrm{SK}n)\,(\rho + x)^{1/2-k}z^{n+k+1}]/((n+1)!\sqrt\rho)$、$g_{kn} = g_{k-1,n}\,(2k-3)(2k+1)/(8k(k+n+1))$ と解ける (§6.2。導出)。**ヘルプの式は、$g_{kn}$ の $(2k-1)!!$ (正しくは $(2k+1)!!$) と、$z$ の冪 $n + k$ (正しくは $n + k + 1$) が違う**。K64 (`gkninit`) の式は正しい。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 長さ | `L` | `al` (`tmulti`) | `el.L`、`pm['al']` |
| 多極の積分強度 $K_n$、$SK_n$ [1/m$^n$] | `K0`〜`K21`、`SK0`〜`SK21` | `ak(0:nmult)` (複素 $K_n + i\,SK_n$。`nmult = 21`)、回した後は `akr0`・`akr` | `pm['akr']` (複素の配列) |
| 曲げ角 | `ANGLE` | `phia` (`tmulti`)、`phi` (`tmulta`) | `el.ANGLE` |
| 端の角 | `E1`・`E2`・`AE1`・`AE2` | `psi1`・`psi2` ($= \mathrm{E1}\cdot\mathrm{ANGLE} + \mathrm{AE1}$) (`tffs.f:1344-1351`) | `_g(el, 'AE1')` ほか |
| ずれ・傾き | `DX`・`DY`・`DZ`・`CHI1`・`CHI2`・`ROTATE` (= `CHI3`)・`DROTATE` | `dx`・`dy`・`dz`・`chi1`・`chi2`・`theta`・`dtheta`、実際に回す角 `theta2` | `pm['dx']`・`pm['theta2']` |
| 刻み | `EPS` | `eps0`、`ndiv0` (`cmp%ivalue(1,p_NM_MULT)`) | `pm['ndiv']`、`pm['eps']` |
| 最高次 | (K・SK で 0 でない最大の $n$) | `nmmax` (`cmp%ivalue(2,p_NM_MULT)`) | `pm['nmmax']` |
| フリンジ | `FRINGE`・`DISFRIN`・`F1`・`F2`・`FB1`・`FB2` | `mfring`、`fringe` (= DISFRIN が 0)、`f1in`・`f2in`・`f1out`・`f2out`、`fb1`・`fb2`、`dofr(0:nmult)` | `pm['mfring']` ほか |
| 数 | (運動量 $p = 1 + \delta$) | `g(i)` ($= \delta$) | `c[5]` |

SAD の座標は $(x, p_x, y, p_y, z, \delta)$、$p_x, p_y$ は $p_0$ で規格化した正準運動量 ([01 章](01_variables_and_conventions.md))。追跡の内部では、キックを $1/(1+\delta)$ で割って $q = p/(1+\delta)$ の単位にする ([07 章](07_thin_multipole.md) の注意と同じ)。

## 3. 多極のポテンシャル (直線、ANGLE = 0)

### 3.1 K_n・SK_n のキック

薄い多極のポテンシャル ([07 章](07_thin_multipole.md) 式 (7.1) の複素の形) を、スキューまで含めて次のように置く:

$$
U = \sum_{n=0}^{n_{\max}}\mathrm{Re}\,\frac{(\mathrm{K}n + i\,\mathrm{SK}n)\,z^{n+1}}{(n+1)!},\qquad z = x + iy
\tag{8.1}
$$

($\mathrm{K}n$、$\mathrm{SK}n$ は積分強度。$H$ の中の項は $U/L$ を長さにわたって足したもの。) ポテンシャルの勾配は、$\partial_x\mathrm{Re}\,F(z) = \mathrm{Re}\,F'$、$\partial_y\mathrm{Re}\,F(z) = \mathrm{Re}\,(iF') = -\mathrm{Im}\,F'$ から

$$
\partial_xU = \mathrm{Re}\,C(z),\quad \partial_yU = -\mathrm{Im}\,C(z),\qquad C(z) \equiv \sum_n\frac{(\mathrm{K}n + i\,\mathrm{SK}n)\,z^n}{n!}
\tag{8.2}
$$

で、キック $\Delta p_x = -\partial_xU$、$\Delta p_y = -\partial_yU$ は

$$
\Delta p_x = -\mathrm{Re}\,C(z),\qquad \Delta p_y = +\mathrm{Im}\,C(z)
\tag{8.3}
$$

になる。(8.3) は、(8.1) の $U$ をハミルトニアンの項 $U/L$ の長さ $L$ の薄いレンズにしたとき $\Delta p = -\nabla U$ ということで、リー演算子の流れ $\exp(:f:)$ では $f = -U$ に当たる ([C_poisson_lie.md](C_poisson_lie.md) の約束: $\exp(:f:)\,p_x = p_x + \partial f/\partial x$。この章は $\exp(:f:)$ を直接は使わないが、$f$ に言い換えるときはこの約束)。

K64 の `tmulk`・`tmulti` の `px(i)=px(i)-dble(cx)`、`py(i)=py(i)+imag(cx)` ($\texttt{cx}$ = $C(z)$ の Horner 形) が (8.3) である。Horner 形は、`cx=(cx+akn(kord))*cx0*aninv(kord)` ($\texttt{cx0} = z$、$\texttt{aninv}(k) = 1/k$、`tffs.f:1162`) を $k = N, N-1, \ldots, 1$ の順に回したもので、$N = 3$ なら順に

$$
\texttt{cx}_3 = a_3\,\frac z3,\qquad
\texttt{cx}_2 = \bigl(a_3\tfrac z3 + a_2\bigr)\frac z2 = a_3\frac{z^2}{3\cdot2} + a_2\frac z2,\qquad
\texttt{cx}_1 = \bigl(\texttt{cx}_2 + a_1\bigr)\,z = a_3\frac{z^3}{3!} + a_2\frac{z^2}{2!} + a_1\,z
$$

(各行は、前の行に $a_{k-1}$ を足して $z/(k-1)$ をかけたもの。$a_n$ にかかる係数は $z^{n}\cdot\prod_{k=1}^{n}(1/k) = z^n/n!$)。最後に $a_0$ を足すと $\texttt{cx}_1 + a_0 = \sum_{n=0}^{3}a_nz^n/n! = C(z)$ になる (一般の $N$ で sympy が差 0、[mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (a))。(§8 の数値確認でも、薄い多極と K64 の追跡が $10^{-19}$ で一致。) `tmulk` (端面の角) は $k = 1$ までで止め、`fe*(cx+akn(0))*cx0` と最後に $a_0$ を足してから $z$ をかける (§4.4)。

### 3.2 SK_n は K_n の回転 (ヘルプの記述の確認)

まず ROTATE の意味を決める。ROTATE $= \theta$ は、`tsolrot` の入口で座標を $(x, y) \to (\cos\theta\,x - \sin\theta\,y,\ \sin\theta\,x + \cos\theta\,y)$、すなわち $z' = e^{i\theta}z$ と回し、要素の中の追跡 (係数が実数の K$n$) をして、出口で $z = e^{-i\theta}z'$ に戻す ([mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) (c): この向きの差 0、逆向きの差 $6 \times 10^{-4}$)。要素の中のポテンシャルは $z'$ で $\mathrm{Re}\,[\mathrm{K}n\,z'^{\,n+1}]/(n+1)!$ なので、実験室の座標 $z$ では $z'^{\,n+1} = e^{i(n+1)\theta}z^{n+1}$ から

$$
\mathrm{Re}\,\frac{\mathrm{K}n\,e^{i(n+1)\theta}\,z^{n+1}}{(n+1)!}
$$

つまり、実験室の複素係数は $\mathrm{K}n\,e^{i(n+1)\theta}$ になる。$\mathrm{K}n$ と $\mathrm{SK}n$ の組 $A = \mathrm{K}n + i\,\mathrm{SK}n$ が $A = \mathrm{K}n\,e^{i(n+1)\theta}$ に等しい ($\mathrm{K}n$ が実数の磁石を ROTATE したもの) ためには、逆に**要素の座標での係数**が

$$
(\mathrm{K}n + i\,\mathrm{SK}n)\;\to\;(\mathrm{K}n + i\,\mathrm{SK}n)\,e^{-i(n+1)\theta}
\tag{8.4}
$$

(新しい座標 $z' = e^{i\theta}z$ で書くと $z^{n+1} = e^{-i(n+1)\theta}z'^{\,n+1}$ だから、$A z^{n+1} = A\,e^{-i(n+1)\theta}z'^{\,n+1}$)。sympy で、$n = 1, 2, 3$ について両辺の実部の差が 0 であることを確認した ([mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (b))。

**SK$n$ との関係**。$\mathrm{K}n = 0$、$\mathrm{SK}n = s > 0$ のとき $A = i\,s = s\,e^{i\pi/2}$。(8.4) の因子が $e^{-i(n+1)\theta} = e^{-i\pi/2}$ になる $(n+1)\theta = \pi/2$、すなわち **$\theta = 90°/(n+1)$** で $A' = s\,e^{i\pi/2}e^{-i\pi/2} = s$ (実数) となる。逆に、実数の K$n$ ($= s$) の磁石を ROTATE $= 90°/(n+1)$ にすると、実験室の係数は $s\,e^{i(n+1)\theta} = i\,s$ で SK$n = s$ と同じ。ヘルプの「SK$n$ > 0 は、K$n$ の磁石を ROTATE = $90°/(n+1)$ にしたもの」と一致する (K1 なら 45°、K2 なら 30°、K3 なら 22.5°)。

### 3.3 基準面 (θ1)

MULT は、K1 の回転の位相を**要素の基準面**として取り込む。$a = \mathrm{K1} + i\,\mathrm{SK1} = |a|\,e^{2i\theta_1}$ と極形式に書き、$\theta_1 = \tfrac12\arg a$ だけ座標を回す。(8.4) の $n = 1$ ($e^{-i(n+1)\theta} = e^{-2i\theta_1}$) で、要素の座標の K1 + i SK1 は $a\,e^{-2i\theta_1} = |a|$ (実数) になる。つまり**スキュー四極の成分が 0 になる** (ヘルプの「reference plane is defined so that the skew quadrupole component becomes zero」)。`akang` (`tfloor.f:530`) は $\theta_1 = \tfrac12\,\mathrm{atan2}(\mathrm{SK1}, \mathrm{K1})$、$\texttt{cr1} = e^{-i\theta_1}$ とする。例: SK1 = 0 で K1 > 0 は $\theta_1 = 0$ ($\texttt{cr1} = 1$)、K1 < 0 は $\mathrm{atan2}(0, \mathrm{K1} < 0) = \pi$ だから $\theta_1 = \pi/2$、$\texttt{cr1} = -i$、K1 の係数は $\mathrm{K1}\cdot(-i)^2 = -\mathrm{K1} = |\mathrm{K1}|$ になる (座標を 90° 回すと、発散四極が収束四極として見える)。K0 は $n = 0$ なので $\mathrm{K0}\,e^{-i\theta_1}$ ($\theta_1 = \pi/2$ では $-i\,\mathrm{K0}$: 水平のキックが垂直のキックに入れ替わって見える)。`tffs.f:1398-1420` は、実際に回す角を `theta2` $= \mathrm{ROTATE} + \mathrm{DROTATE} + \theta_1$、各成分を `akr(n)` $= (\mathrm{K}n + i\,\mathrm{SK}n)\cdot\texttt{cr1}^{\,n+1}$ ($\texttt{cr}$ を $n$ のたびに `cr1` 倍して作る) と決める。これで、K1 を持つ MULT の本体 (§4) は四極の主軸の座標で解ける。

## 4. ANGLE = 0 の構成

### 4.1 ハミルトニアンの分割

回した座標 (K1 が実数) での $H$ は

$$
H = -\sqrt{p^2 - p_x^2 - p_y^2} + \frac{E}{p_0v_0} + \underbrace{\frac1L\Bigl[\tfrac{\mathrm{K1}}{2}(x^2 - y^2) + \mathrm{K0}\,x - \mathrm{SK0}\,y\Bigr]}_{H_{\rm body}} + \underbrace{\frac1L\sum_{n\ge2}\mathrm{Re}\,\frac{\mathrm{akr}(n)\,z^{n+1}}{(n+1)!}}_{U_{\rm kick}}
\tag{8.5}
$$

(第 1 項と第 2 項は 02 章 (2.1)。K0 = `akr(0)` の実部、SK0 は虚部。) $H_{\rm body}$ は、運動学の項 $-\sqrt{\cdot}$ と合わせて**四極 + 一様な双極の本体**で、[04 章](04_quad.md) の $H_2$ と $\Delta H$ で厳密に (運動学の部分は刻んで) 解ける。

**双極の項は中心のずれに吸収される**。$\mathrm{K1} \ne 0$ のとき、$\tfrac{\mathrm{K1}}{2L}(x^2 - y^2) + \tfrac{1}{L}(\mathrm{K0}\,x - \mathrm{SK0}\,y)$ を平方完成する:

$$
\frac{\mathrm{K1}}{2L}\Bigl[(x + d_x)^2 - (y + d_y)^2\Bigr] - \frac{\mathrm{K1}}{2L}\bigl(d_x^2 - d_y^2\bigr),\qquad d_x = \frac{\mathrm{K0}}{\mathrm{K1}},\ d_y = \frac{\mathrm{SK0}}{\mathrm{K1}}
\tag{8.6}
$$

(確認: $-\tfrac{\mathrm{K1}}{2L}\cdot 2d_yy = -\tfrac{\mathrm{SK0}}{L}y$、$\tfrac{\mathrm{K1}}{2L}\cdot 2d_xx = \tfrac{\mathrm{K0}}{L}x$。) 定数項は軌道に効かない (運動方程式に出ない) ので、本体の運動は、$(x, y)$ を $(x + d_x, y + d_y)$ に替えた四極の運動と同じ。`tsolqu` は `dx0=ak0x/ak`、`dy0=ak0y/ak`、`xi=x(i)+dx0` で本体の解 (04 章 (4.9)(4.10)) を $\texttt{xi}$、$\texttt{yi}$ に入れる (`tsolqu.f:55-56, 83-84`)。定数項そのものは $p$ (= $1 + \delta$) に依らない項なので $z$ にも直接は効かない。$z$ の更新式 (`w1*(xi**2*xs1-yi**2*xsh2)` など) は、振動の中心から測った座標 $\texttt{xi}$、$\texttt{yi}$ で書かれていて、中心のずれは $x \to \texttt{xi}$、$y \to \texttt{yi}$ の置き換えを通して入る (`tsolqu.f:83-84, 102-109`)。ここまでは K64 のコードの読みで、$z$ の式そのものの導出は [04_quad.md](04_quad.md) (4.12)〜(4.14) にある。K1 = 0 のときは `tdrift` に K0・SK0 の一様なキックを渡す (`tsolqu.f:34-40`)。

### 4.2 分割と順序 (`tmulti`)

$H = H_{\rm body}' + U_{\rm kick}$ ($H_{\rm body}'$ は運動学 + $H_{\rm body}$) とみて、2 次の対称な分割 ([03 章](03_symplectic_integrators.md) §4) をする。長さ $L$ を $\mathrm{ndiv}$ 個に分け、$h = L/\mathrm{ndiv}$:

$$
M = B\!\left(\tfrac h2\right)\,\mathcal K\,B(h)\,\mathcal K\,B(h)\cdots\mathcal K\,B\!\left(\tfrac h2\right)\qquad(\mathcal K\text{ は }\mathrm{ndiv}\text{ 個、}B(h)\text{ は }\mathrm{ndiv}-1\text{ 個})
\tag{8.7}
$$

$B(a)$ は長さ $a$ の分の $H_{\rm body}'$ の流れ (`tsolqu(al1, ak1, bzs, ak01x, ak01y, …)`、K1・K0 の強さは $a/L$ 倍)、$\mathcal K$ は強さが $1/\mathrm{ndiv}$ の $n \ge 2$ の薄いキック (8.3)。コードでは、`do m=1,ndiv` の中で本体を呼んだあとに `al1=al*wm` を更新し (`wm` は $m = \mathrm{ndiv}$ だけ $w_i/2$、それ以外 $w_i = 1/\mathrm{ndiv}$。`tmulti.f:205`)、キック (`akrm(n)=akr(n)*wi`) を入れる。最初の本体は `al1=al*wi*.5` (`tmulti.f:157-163`)、最後の半分はループの外 (`tmulti.f:244-248`)。キックの式は Horner 形 $\sum_{n\ge2}\mathrm{akrm}(n)\,z^n/n!\,/(1+\delta)$ で (`tmulti.f:209-221`)、07 章のキックと同じ。**3 次以上の力勾配型の補正 (07 章 §6) は MULT には無い**。K1 ≠ 0 のとき $B$ は $\mathcal K$ と**可換でない**ので、誤差は 2 次 (Strang) になる。

**誤差の次数の導出**。(8.7) の 1 刻み (長さ $h$) は、隣り合う本体 $B(h/2)B(h/2) = B(h)$ をまとめると $e^{hX/2}\,e^{hY}\,e^{hX/2}$ の形 ($X$ = 本体 (運動学 + K1・K0) の生成、$Y$ = キックの生成 (強さを長さあたりにしたもの。$e^{hY}$ が 1 刻みのキック)) である。対称な 3 つの積の BCH 展開 ([03_symplectic_integrators.md](03_symplectic_integrators.md) §4。係数 $-1/24$、$1/12$ の導出と、ランダムな行列での数値確認は [04_quad.md](04_quad.md) §6.1) は、$[X, Y]$ の入れ子の括弧を使って

$$
e^{hX/2}e^{hY}e^{hX/2} = \exp\Bigl\{h(X + Y) + h^3\Bigl(-\tfrac1{24}[X,[X,Y]] + \tfrac1{12}[Y,[Y,X]]\Bigr) + O(h^5)\Bigr\}
$$

(2 次の項 $h^2$ は対称性で消える。ランダムな行列で数値的に確かめた: $h$ を半分にすると残りの比が 1/4 になる、[mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (f2))。$\mathrm{ndiv}$ 個を重ねた誤差は $\mathrm{ndiv}\cdot h^3 = L^3/\mathrm{ndiv}^2$ に比例する。$X$ と $Y$ が可換 (K1 = 0 で運動学の項も無視できる) なら括弧が 0 で誤差も 0 のはずだが、運動学の項 $-\sqrt{p^2 - p_x^2 - p_y^2}$ はキックの $x$ と非可換なので、K1 = 0 でも誤差は 2 次で残る。1 自由度の調和振動子 (力の勾配 $\omega^2$) の場合は、1 刻みの行列が具体的に書ける:

$$
M = \begin{pmatrix}1 & h/2\\0 & 1\end{pmatrix}\begin{pmatrix}1 & 0\\-\omega^2h & 1\end{pmatrix}\begin{pmatrix}1 & h/2\\0 & 1\end{pmatrix}
= \begin{pmatrix}1 - \tfrac{(\omega h)^2}{2} & h\bigl(1 - \tfrac{(\omega h)^2}{4}\bigr)\\ -\omega^2h & 1 - \tfrac{(\omega h)^2}{2}\end{pmatrix}
$$

(2 つの行列の積を順に計算した。第 1 の積は $\begin{pmatrix}1 - \omega^2h^2/2 & h/2\\ -\omega^2h & 1\end{pmatrix}$、これに右から $\begin{pmatrix}1 & h/2\\0&1\end{pmatrix}$ をかけると上の形。) 厳密な行列は $\cos\omega h$、$\sin(\omega h)/\omega$、$-\omega\sin\omega h$ で、差は $M_{11}$ で $-(\omega h)^4/24$、$M_{12}$ で $-h(\omega h)^2/12$。位相は $\cos\theta = \tfrac12\mathrm{tr}M = 1 - (\omega h)^2/2$ から $\theta = \omega h\,[1 + (\omega h)^2/24 + \cdots]$ なので、$\mathrm{ndiv}$ 個で位相誤差は $\omega L\,(\omega h)^2/24 = (\omega L)^3/(24\,\mathrm{ndiv}^2)$ (sympy、同 (f))。

放射 ($\texttt{krad}$) があると、本体ごとに `tradk` (光子の放出) を呼び、`ndiv` を `ndivrad` との大きい方にする (`tmulti.f:119`)。SPAC (空間電荷) があると `ndiv` を $|L|/(0.05\,\mathrm{EPS})$ 以上にする (`:110-113`)。この章では扱わない。

### 4.3 刻みの数 `ndiv`

`tffs.f:1423-1434`: $\varepsilon = \varepsilon_{00m}\cdot\mathrm{EPS}$ (EPS = 0 は 1)、$\varepsilon_{00m} = 0.005$ (`toplvl.f:292`)、$a_{\max} = $`ampmaxm`$ = 0.05$。

$$
\mathrm{ndiv} = \min\Bigl(1000,\ \max\bigl(1,\ \max_{2\le n\le n_{\max}}\bigl\lfloor\sqrt{a_{\max}^{\,n-1}\,\|\mathrm{K}n + i\,\mathrm{SK}n\|\,L\,/(6\,(n-1)!\,\varepsilon)}\bigr\rfloor + 1\bigr)\Bigr)
\tag{8.8}
$$

(`fact(n-1)` は階乗の表、`hypot` は絶対値。) $n = 0, 1$ は厳密な本体で処理するので入らない。

**各因子の意味 (導出)**。$n$ 極 ($n \ge 2$) の力は、横方向の座標 $x$ で $F_x = -\mathrm{Re}\,C \to -\mathrm{K}n\,x^n/(n!\,L)$ (長さあたり、(8.2)(8.3) で $y = 0$)。ポテンシャルを $V(x)$ ($F_x = -V'$) とすると、その $x$ での勾配 (力の傾き) は $\omega^2(x) \equiv V''(x) = |\mathrm{K}n|\,x^{n-1}/((n-1)!\,L)$。振幅 $x = a_{\max} = 0.05$ m で、長さ $L$ をかけて無次元にした

$$
Q \equiv \omega^2L^2 = \frac{a_{\max}^{\,n-1}\,|\mathrm{K}n|\,L}{(n-1)!}
$$

が、(8.8) の根号の中の分子 ($Q$ は、$n$ 極の力の傾きの $a_{\max}$ での値 $\times L^2$)。(8.8) は $\mathrm{ndiv}^2 \ge Q/(6\varepsilon)$、すなわち

$$
\varepsilon \;\ge\; \frac{Q}{6\,\mathrm{ndiv}^2} = \frac{(\omega h)^2}{6}\qquad(h = L/\mathrm{ndiv})
\tag{8.8a}
$$

を要求する形。K64 のソースとヘルプには、$\varepsilon$ が何の誤差の上限かの説明が無い。以下では、**$\varepsilon$ を「$a_{\max}$ での $n$ 極の運動量の跳び (偏向) の相対誤差の上限」と読む (推定)** と、係数 $1/6$・$Q$ の冪 $Q^1$・$(n-1)!$ がすべて導出できることを示す (確認: [mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (a))。

**段 1: 分割の実効ハミルトニアン**。(8.7) の 1 刻みは $e^{hX/2}e^{hY}e^{hX/2}$。2 自由度 $(x, y)$ で、本体の運動学を放物近似した $X$ の生成関数を $f_X = -T$ ($T = p^2/2$)、キックを $f_Y = -V$ とおく (C 付録の約束 $f = -H$)。ポアソン括弧 $[f, g] = \sum_i(\partial_{q_i}f\,\partial_{p_i}g - \partial_{p_i}f\,\partial_{q_i}g)$ を 1 つずつ計算する:

- $[f_X, f_Y] = [T, V] = -\,p_i\,\partial_iV$ ($T$ は $q$ に依らず、$V$ は $p$ に依らないので、第 2 項だけ残る)。
- $[f_X, [f_X, f_Y]] = [-T,\ -p_j\partial_jV] = [T,\ p_j\partial_jV] = -\,p_ip_j\,\partial_i\partial_jV$ (添字は和)。
- $[f_Y, f_X] = [V, T] = +\,p_i\,\partial_iV$、$[f_Y, [f_Y, f_X]] = [-V,\ p_j\partial_jV] = -\,|\nabla V|^2$ ($p$ について微分されるのは $p_j\partial_jV$ だけで $\partial_{p_j}(p_k\partial_kV) = \partial_jV$)。

§4.2 の BCH の $h^3$ の項 $-\tfrac1{24}[f_X,[f_X,f_Y]] + \tfrac1{12}[f_Y,[f_Y,f_X]]$ に入れると $\tfrac1{24}p_ip_j\partial_i\partial_jV - \tfrac1{12}|\nabla V|^2$。$f = -h\tilde H$ なので 1 刻みの実効ハミルトニアン $\tilde H$ は

$$
\tilde H = T + V + h^2\Bigl[\frac1{12}\,|\nabla V|^2 - \frac1{24}\,p_ip_j\,\partial_i\partial_jV\Bigr] + O(h^4)
\tag{8.8b}
$$

(sympy で (8.8b) の BCH の計算と一致。**確認**: 6 極 ($n = 2$、$K2 = 40$ /m²、$L = 1$ m、初期値 $(4, 0.8, -3, -0.4)$ mm・mrad) の $\tilde H$ の厳密な流れと、厳密な運動方程式との差は、`ndiv` = 5 で $(-4.14, -6.35, -4.48, 5.75)\times10^{-7}$、分割の追跡と厳密解との差は $(-4.15, -6.35, -4.48, 5.74)\times10^{-7}$ (最大の差 $1.4\times10^{-9}$)、`ndiv` = 17 でも $3.6\times10^{-10}$ 以内で一致。つまり (8.8b) が分割の誤差の主要項を全部持つ。)

**段 2: 調和振動子 ($V = \omega^2x^2/2$)**。(8.8b) は $\tilde H = A\,p^2/2 + B\,x^2/2$、$A = 1 - (\omega h)^2/12$、$B = \omega^2\bigl[1 + (\omega h)^2/6\bigr]$ ($|\nabla V|^2/12 = \omega^4x^2/12$ が $x^2$ に、$-p^2V''/24$ が $p^2$ に入る)。この $\tilde H$ の運動は $\dot x = A\,p$、$\dot p = -B\,x$ なので、周波数は $\Omega = \sqrt{AB} = \omega\bigl[1 + (\omega h)^2/24 + \cdots\bigr]$ (位相誤差 $(\omega h)^2/24$ が、(§4.2) の $\theta = \omega h[1 + (\omega h)^2/24]$ と同じ)。長さ $L$ の後の運動量は $p(L) = -x_0\,\Omega\sin(\Omega L)/A$ (運動の解 $x = x_0\cos\Omega t$ と $p = \dot x/A$ から) で、$\Omega L \ll 1$ では $\sin(\Omega L) \simeq \Omega L$ なので $p(L) \simeq -x_0\,(\Omega^2/A)\,L = -B\,x_0L$ ($\Omega^2 = AB$)。厳密な (分割なしの) 跳びは $-\omega^2x_0L$ なので、**運動量の跳びの相対誤差 $= B/\omega^2 - 1 = (\omega h)^2/6$**。この $1/6$ は、$|\nabla V|^2/12 = \omega^4x^2/12$ が $\tfrac12Bx^2$ の $x^2$ の係数に入るので $B$ では 2 倍の $\omega^4h^2/6$ になる、という $\tfrac1{12}\times2$ (位相の $1/24$ だけではない)。この $(\omega h)^2/6$ が (8.8a) の右辺そのもの: 係数 $1/6$ は上の 3 行の括弧から出た $\tfrac1{12}\times2$ で、$Q^1$ なのは、跳び $\propto \omega^2L \propto Q$ の相対値だから ($Q^{3/2}$ は位相の絶対値の冪で、別の量)。

**一般の $n$ 極**。1 自由度の $V' = \mathrm{K}n\,x^n/(n!\,L)$ で、運動量の跳びは $-V'\,L$ の分割誤差を持つ。$-p^2V''$ の項は $\omega L\ll1$ で $O(p^2)$ なので落とし、(8.8b) の力は $-\partial_x\bigl[V + \tfrac{h^2}{12}V'^2\bigr] = -V'\bigl[1 + \tfrac{h^2}6V''\bigr]$。したがって跳びの相対誤差は $h^2V''(a_{\max})/6 = Q/(6\,\mathrm{ndiv}^2)$ で、**(8.8) の $a_{\max}^{n-1}$ と $(n-1)!$ ($V'' = K_na^{n-1}/((n-1)!L)$ から) を含めて (8.8a) と同じ**。2 自由度の 6 極 ($V = \tfrac{K}{6L}(x^3 - 3xy^2)$) では、$\Delta\vec p = -\tfrac{h^2L}{6}(\nabla\nabla V)(\nabla V)$ (静止の近似)。経路の中点で評価した値 ($\mathrm{ndiv} = 5$: $(-6.9, 5.1)\times10^{-7}$、17: $(-6.0, 4.4)\times10^{-8}$) が、測定 ($(-6.35, 5.74)\times10^{-7}$、$(-5.50, 4.97)\times10^{-8}$) と 10% 程度で合う (残りは経路上の勾配の変化と $p^2V''$ の項)。

**注意**: この読みが成り立つのは $\omega L \ll 1$ (軌道が 1 回の通過で大きく曲がらない) のとき。$a_{\max}$ では $Q = 2$、$\omega L = 1.4$ で、$\varepsilon = Q/(6\,\mathrm{ndiv}^2)$ は目安に過ぎない (例: EPS = 1 ($\varepsilon = 0.005$) で $\mathrm{ndiv} = 9$、$Q/(6\cdot81) = 0.0041 \le 0.005$)。この読みは、K64 の式の形と全部の数字 (係数、$Q$ の冪、階乗) と合うという意味での推定で、K64 の著者の定義を確認したものではない。(8.8) を「全位相誤差 $Q^{3/2}/(24\,\mathrm{ndiv}^2) \le \varepsilon$」と読むと $\mathrm{ndiv}^2 = Q^{3/2}/(24\varepsilon)$ で、(8.8) の $\mathrm{ndiv}^2$ との比は $4/\sqrt Q$ ($Q = 2$ で 2.8) になり、係数も冪も合わない。

数値例 (K2 = 40 /m²、$L$ = 1 m、$n = 2$): $Q = 0.05 \times 40 \times 1/1! = 2$。EPS = 4, 2, 1, 0.5, 0.25 ($\varepsilon = 0.005\,\mathrm{EPS} = 0.02, 0.01, 0.005, 0.0025, 0.00125$) に対し $\sqrt{Q/(6\varepsilon)} = 4.08, 5.77, 8.16, 11.5, 16.3$、$\mathrm{ndiv} = 1 + \lfloor\cdot\rfloor = 5, 6, 9, 12, 17$ (pysad の `ndiv` と一致、[mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (f3))。実際の振幅 4 mm では $(\omega L)^2 = 40 \times 0.004 \times 1 = 0.16$ で、運動量の跳びの見積り $a\,Q^2/(6\,\mathrm{ndiv}^2L)$ ($a = 4$ mm、$Q = 0.16$) は $\mathrm{ndiv} = 5, 17$ で $6.8\times10^{-7}$、$5.9\times10^{-8}$。下の厳密解との差の測定値 $6.4\times10^{-7}$、$5.5\times10^{-8}$ は、その 0.93 倍 (残りは上の 10%)。

**位相誤差だけの見積りが測定の約 1.5 倍になる理由**。位相誤差だけの見積り $a\,(\omega L)^3/(24\,\mathrm{ndiv}^2)$ ($4.3\times10^{-7}$、$3.7\times10^{-8}$) は、運動量の跳びの誤差のうち位相の $1/24$ 分だけを数え、運動量と座標の比 ($p = \dot x/A$ の $A$ など) の変化を数えていない (調和振動子で、運動量の跳びの相対誤差は位相の $1/24$ ではなく $1/6$)。2 つの見積りの比は $[Q^2/6]/[Q^{3/2}/24] = 4\sqrt Q = 4\,\omega L = 1.6$ ($Q = 0.16$) で、測定/位相の見積り $= 1.49$ (`ndiv` = 5 でも 17 でも) に近い (残り 7% は上の 10% と同じ原因)。つまり**約 1.5 倍は、非線形の力 ($x^2$) の 1 次の振幅依存ではなく、位相の見積りが誤差の一部しか数えていなかったため**。

**誤差の次数**: ANGLE = 0、K2 = 40 /m² (積分強度)、$L$ = 1 m、初期値 $(4, 0.8, -3, -0.4)\times$ (mm, mrad, mm, mrad) で、厳密な運動方程式 (DOP853) との差は、EPS が 4 → 0.25 のとき $6.4 \times 10^{-7} \to 5.5 \times 10^{-8}$。誤差の比は `ndiv` の比の 2 乗とほぼ等しい (1.44 対 1.44、2.25 対 2.25、1.78 対 1.78、2.01 対 2.01)。**2 次**である ([mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) (a))。K1 = 0.8 を加えた (b) でも、`ndiv` を 6 → 12 で $1.7 \times 10^{-6} \to 4.2 \times 10^{-7}$ (2 次)。

### 4.4 フリンジ

- **非線形フリンジ**: `fringe` (DISFRIN = 0) のとき、`mfring ≠ 2` で入口、`mfring ≠ 1` で出口に、$n$ ごとに `ttfrins(…, n*2+2, ±akr(n), al, bzs)` を呼ぶ (`tmulti.f:124-132, 273-281`。`n*2+2` は極の数 $2(n+1)$)。対象は `dofr(n)` が真の成分 (`tffs.f:1436-1443`): K$n$ ≠ 0 かつ ($n \le 2$ か、$\|K_n\|/L\cdot a_{\max}^{n+1}\,\texttt{aninv}(n+1) > \varepsilon$。`aninv(n+1)` は $1/(n+1)$ で、$1/(n+1)!$ ではない。判定量 $(\|K_n\|/L)\,a_{\max}^{n+1}/(n+1)$ は、$n$ 極の場 $K_n x^n/(n!\,L)$ を $x$ で積分した大きさ ($\sim K_na^{n+1}/((n+1)!\,L)$) に近いが $n!$ が無い見積りで、ソースに根拠は書かれていない)。判定量を $D_n \equiv (\|K_n\|/L)\,a_{\max}^{n+1}/(n+1)$ と書き、フリンジの位置のずれ ((7.17) の $m = n + 1$、$\kappa = K_n/L$。第 1 項の大きさ $|\Delta w|/a = \kappa\,a^{m}/(4\,m!\,(m+1))$) と比べると $D_n\big/(|\Delta w|/a) = 4\,n!\,(n + 2)$ (sympy: [mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (b))。つまり $D_n > \varepsilon$ は、「フリンジの位置のずれが $a_{\max}$ の $\varepsilon/(4\,n!\,(n+2))$ 倍より大きい」と読め、$\varepsilon = 0.005$ で $n = 2$: $1.6\times10^{-4}$、$n = 3$: $4.2\times10^{-5}$、$n = 4$: $8.7\times10^{-6}$、$n = 5$: $1.5\times10^{-6}$ (次数が上がると、より小さいずれまで入れる。$n!$ が無いことの効果)。**ソースの意図 ($n!$ を入れなかった理由) は確認できず、この読みは式の変形である**。

(足りないもの: K64 の著者が $\varepsilon$ をフリンジの何の誤差の上限としたかの定義。上の $n$ ごとの閾値が意図した値かどうかは、このソースからは決まらない。) 非線形フリンジの式は [07 章](07_thin_multipole.md) §9 (Maxwell から)。

- **線形フリンジ (F1・F2)**: 四極 (K1) の線形フリンジ ([04 章](04_quad.md) §7.4 と同じ式)。`f1in`・`f2in`・`f1out`・`f2out` は `tsetfringep` が $\mathrm{K1}/L$ から決める (`tffs.f:1380-1385`)。出口は $a \to -a$。`mfring = 1, 3` が入口、`2, 3` が出口。
- **K0 の端 (FB1・FB2)**: `tblfri(…, al, ±akr(0), fb)` (入口は `akr(0)`・`fb1`、出口は `-akr(0)`・`fb2`)。偏向磁石の線形フリンジ V ([05 章](05_bend.md) 式 (5.42)) を、$\rho_b = L/\mathrm{K0}$ ($x$ 方向)・$L/\mathrm{SK0}$ ($y$ 方向) で入れる。 垂直 (SK0) への拡張が (5.42) を $+90°$ 回したものになることは [06 章](06_edge_fringe.md) §9.2。
- **端面の角 (AE1・AE2)**: `tmulk(…, al, akr, psi, nmmax)` (`tmulti.f:4-31, 166-168, 250-252`。`fe` は 16 行目):

$$
C_\psi(z) = -\frac{\tan\psi}{L}\;z\sum_{k=0}^{n_{\max}}\frac{(\mathrm{K}k + i\,\mathrm{SK}k)\,z^k}{k!},\qquad
\Delta p_x = -\mathrm{Re}\,C_\psi/(1+\delta),\quad \Delta p_y = +\mathrm{Im}\,C_\psi/(1+\delta)
\tag{8.9}
$$

($\texttt{fe} = -\tan\psi/L$、`cx=fe*(cx+akn(0))*cx0`。§3.1 の Horner 形から、`tmulk` の `cx` は $-(\tan\psi/L)\,z\,C(z)$ に等しい。sympy で差 0、[mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (a)(c)。) **この式の高次の項 ($z\,C(z)$ の形) は、ソースに根拠が書かれていない**。次の手順で導ける (仮定つき)。以下 $t \equiv \tan\psi$。

1. **中央面 ($y = 0$) の値: 余分な長さ × 力**。入口の面は軌道に垂直な平面から $\psi$ だけ傾いていて、面上の点では $s = x\,t$ (矩形偏向磁石 $\psi = \mathrm{ANGLE}/2 > 0$ では、曲げの外側が $x > 0$ で、そこでは磁場の領域が短い。05 章 §5 の符号と同じ向き)。面の方程式を $s = s_f(x) = x\,t$ ($t \equiv \tan\psi$。面は $x$–$s$ 平面の直線で、$y$ に依らない) とおくと、この粒子が磁場の中を走る区間は $s_f \le s \le L$ で、その長さは $\int_{x t}^{L}ds = L - x\,t$。公称の全長 $L$ に対する差が**余分な長さ $-x\,t$** (出口も、`tmulk` に同じ符号の `psi2` を渡すので、$-x\,t$)。長さあたりの力は (8.3) を $L$ で割った $(-\mathrm{Re}\,C,\ +\mathrm{Im}\,C)/L$ なので、
$$
\Delta p_x = -(-x\,t)\,\frac{\mathrm{Re}\,C(x)}{L} = \frac{t\,x}{L}\,\mathrm{Re}\,C(x),\qquad
\Delta p_y = (-x\,t)\,\frac{\mathrm{Im}\,C(x)}{L} = -\frac{t\,x}{L}\,\mathrm{Im}\,C(x)
$$
(縦方向の磁場 $B_s$ による力は $\vec v\times\vec B_s$ の横成分で、横方向の速度 $v_\perp$ に比例する ($\propto v_\perp B_s$)。$v_\perp \ll 1$ の 1 次の近似では力に寄与しない。経路は $s$ 軸に平行な直線とした。) これが (8.9) の $y = 0$ での値に一致する (sympy で差 0、同 (c) の最終行)。
2. **$y \ne 0$ への拡張: 調和性**。面を横切る薄い層で積分した横方向の磁場がキックを決めるとして、キックの場 $(\Delta p_x, \Delta p_y)$ は、真空の Maxwell 方程式から、$\nabla\times\vec B = 0$ の $s$ 成分で $\partial_x\Delta p_y - \partial_y\Delta p_x = 0$ (勾配)、$\nabla\cdot\vec B = 0$ を $s$ で積分した式 ($B_s$ は領域の外で 0) から $\partial_x\Delta p_x + \partial_y\Delta p_y = 0$ (調和) を満たす。この 2 式は、複素量 $\Delta p_x - i\,\Delta p_y$ が $z$ の**解析関数**であるという条件 (コーシー・リーマンの式) と同じ。解析関数は実軸 ($y = 0$) の値で一意に決まるので、1. の値 $(t/L)\,x\,C(x)$ の $x \to z$ の拡張
$$
\Delta p_x - i\,\Delta p_y = \frac tL\,z\,C(z)\ \Longleftrightarrow\ C_\psi = -\frac tL\,z\,C(z)
$$
がただ一つの答えで、(8.9) と一致する。sympy で、(8.9) がコーシー・リーマンの式、$\mathrm{rot} = 0$、$\mathrm{div} = 0$ を満たすこと、ポテンシャル $\Phi = -(t/L)\,\mathrm{Re}\sum_n a_n z^{n+2}/((n+2)\,n!)$ ($\Delta p = -\nabla\Phi$、$\nabla^2\Phi = 0$) を持つことを確認した (同 (c))。C 付録の約束では、$\exp(:f:)\,p_x = p_x + \partial f/\partial x$ なので、この端面の薄いキックは $f = -\Phi = (t/L)\,\mathrm{Re}\sum_n a_n z^{n+2}/((n+2)\,n!)$。
3. **例**。K0 だけ: $\Delta p_x = +(t\,\mathrm{K0}/L)\,x$、$\Delta p_y = -(t\,\mathrm{K0}/L)\,y$ で、偏向磁石の端面の標準の収束 (水平 $\tan\psi/\rho$、垂直 $-\tan\psi/\rho$、[05 章](05_bend.md) §5) と同じ向き・大きさ。K1 だけ: $\Delta p_x = (t\,\mathrm{K1}/L)(x^2 - y^2)$、$\Delta p_y = -2(t\,\mathrm{K1}/L)\,xy$。
4. **仮定と限界**。(a) キックは横磁場の面通過の 1 次の積分 (経路は直線)。(b) 層の厚さは 0 の極限 (硬い端)。$n \ge 1$ の磁場では、この硬い端の $y \ne 0$ の値が 1. の「長さ × 力」だけでは決まらない。これを式で見る。「長さ × 力」を $y \ne 0$ へ素直に延ばす (余分な長さ $-x\,t$ は面が $x$ だけで決まるので $y$ に依らず、力は $z = x + iy$ での値) と

$$
\Delta p_x = \frac{t\,x}{L}\,\mathrm{Re}\,C(x + iy),\qquad \Delta p_y = -\frac{t\,x}{L}\,\mathrm{Im}\,C(x + iy)
\tag{8.9a}
$$

(複素の形では $\Delta p_x - i\,\Delta p_y = \dfrac{t\,x}{L}\,C(z)$)。この場の回転 $\mathrm{rot} \equiv \partial_x\Delta p_y - \partial_y\Delta p_x$ を、積の微分と、$C$ が解析的であることによるコーシー・リーマンの式 ($\partial_x\mathrm{Im}\,C = \mathrm{Im}\,C'$、$\partial_y\mathrm{Re}\,C = -\mathrm{Im}\,C'$) で計算する。$\partial_x\Delta p_y = -\dfrac tL\bigl[\mathrm{Im}\,C + x\,\mathrm{Im}\,C'\bigr]$、$\partial_y\Delta p_x = \dfrac{t\,x}L\,(-\mathrm{Im}\,C')$ なので

$$
\mathrm{rot} = -\frac tL\bigl[\mathrm{Im}\,C + x\,\mathrm{Im}\,C'\bigr] + \frac{t\,x}{L}\,\mathrm{Im}\,C' = -\frac tL\,\mathrm{Im}\,C(z),\qquad \mathrm{div} = \partial_x\Delta p_x + \partial_y\Delta p_y = \frac tL\,\mathrm{Re}\,C(z)
$$

K0 だけ ($C = \mathrm{K0}$ は実数) では $\mathrm{Im}\,C = 0$ で $\mathrm{rot} = 0$ (一様な力に長さの差をかけたもので、勾配)。K1 では $C = \mathrm{K0} + \mathrm{K1}\,z$、$\mathrm{Im}\,C = \mathrm{K1}\,y$ なので **$\mathrm{rot} = -t\,\mathrm{K1}\,y/L \ne 0$** (K1 の式 $\Delta p_x = t\,\mathrm{K1}\,x^2/L$、$\Delta p_y = -t\,\mathrm{K1}\,xy/L$ は勾配でなく、調和でもない)。(8.9) との差は、複素の形で $\dfrac tL\,z\,C - \dfrac{t\,x}L\,C = \dfrac{i\,t\,y}{L}\,C(z)$。K1 だけなら $\dfrac{i\,t\,y}L\,\mathrm{K1}\,(x + iy)$ の実部と (符号を替えた) 虚部から、$(\Delta p_x,\ \Delta p_y)$ に **$(-t\,\mathrm{K1}\,y^2/L,\ -t\,\mathrm{K1}\,xy/L)$** を足したものが (8.9) で、$\Delta p_x = t\,\mathrm{K1}(x^2 - y^2)/L$、$\Delta p_y = -2\,t\,\mathrm{K1}\,xy/L$ となり、rot・div が 0 の調和な場になる ([mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (c)。$C = \mathrm{K0} + \mathrm{K1}z + \mathrm{K2}z^2/2$ でも $\mathrm{rot} = -(t/L)\,\mathrm{Im}\,C$、(8.9) $-$ (8.9a) $= (i\,t\,y/L)\,C$ を sympy で確認)。調和性が必要とする補正を、(8.9) は 2 の拡張として全部入れている。

**pysad の v2 との違い**。pysad の v2 は、この端を面の厳密な幾何 ($z_f = x\tan\psi/(1 - (p_x/p_z)\tan\psi)$ の位置まで進み、本体の場で $-z_f$ だけ戻す。`_mult_edge_exact`、`track6.py:2040-2088`) と、縦方向の収束 $\Delta p_y = -(\tan(\psi \pm \theta) - \tan(\pm\theta))\,\mathrm{Im}\,\psi_{\rm pot}$ ($\psi_{\rm pot} = \sum_n(\mathrm{K}n/L)\,z^{n+1}/(n+1)!$、$\theta = \arctan(p_x/p_z)$。入口が $+$、出口が $-$) で置き換える (docs/30 の B223・B224・B292。k64 は `tmulk` のまま)。両者を小さな振幅で比べた ([mult_05_edge_v2_vs_k64.out.txt](outputs/mult_05_edge_v2_vs_k64.out.txt)。$L = 0.02$ m、AE1 = 0.05、$(x, y) = (a, a)$): K0 は $p_x$ の差が 0、$p_y$ の差 $< 3 \times 10^{-11}$ (v2 = K64 の収束)。K1 = 0.3 は $p_y$ の差が $2 \times 10^{-10}$ 以下で合い、**$p_x$ の差が $+t\,\mathrm{K1}\,y^2/L$** ($t\,\mathrm{K1}\,a^2/L = 7.51 \times 10^{-7}$ に対し、測定は $a = 1$ mm で $7.48 \times 10^{-7}$、$a = 2$ mm で $2.99 \times 10^{-6}$ (4 倍)) になる。つまり v2 は $\Delta p_x = t\,\mathrm{K1}\,x^2/L$ (余分な長さの項だけ) で、(8.9) の $-t\,\mathrm{K1}\,y^2/L$ を持たない。2 の調和性の議論からは、この項が運動量だけのキックを勾配 (シンプレクティック) にするために要る項に当たる。v2 の端の写像が $K1 \ne 0$ でシンプレクティックでないこと (docs/30 B224) とも整合するが、**どちらが硬い端の厳密な答えかは、この章の導出 (1 次・硬い端の仮定) からは確定できない** (v2 の基準は、B292 の独立な Maxwell 場の数値積分。$p_x$ の $y^2$ 項をその基準で検証したかは未確認)。

## 5. ずれ・回転 (DX、DY、DZ、CHI1、CHI2、ROTATE)

ヘルプ (`misalignments`) は、入口で粒子の座標を剛体の変換

$$
\begin{pmatrix}x\\y\\\Delta s\end{pmatrix}_1 = R_3(\chi_3)\,R_2(\chi_2)\,R_1(\chi_1)\begin{pmatrix}x - \mathrm{DX}\\ y - \mathrm{DY}\\ -\mathrm{DZ}\end{pmatrix}_0,\qquad
R_1 = \begin{pmatrix}c_1 & 0 & -s_1\\0&1&0\\s_1&0&c_1\end{pmatrix},\ R_2 = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}
\tag{8.10}
$$

($R_3$ は $z$ 軸まわりの回転。$c_i = \cos\mathrm{CHI}i$、$s_i = \sin\mathrm{CHI}i$。出口で逆変換) とする。粒子の線 (位置と向き) をこの変換で動かし、$\Delta s = 0$ の面まで直進させたのが入口の写像で、`tsolrot` の $B_z = 0$ の枝 (`tsol.f` の `ent = .true.`): `x-dx`、`y-dy` → $\mathrm{CHI}1$ の回転 (`cchi1*x - schi1*(s0s - dz)`、`px` も同様) → $\mathrm{CHI}2$ の回転 → 面までの直進 (`a=ds2/pz2`、`x-a*px`、`y-a*py`、`z+a`) → ROTATE (`theta2`) の回転 (`tsol.f` の `chi3`)、の順。MULT の `tmulti` は、`-chi1`、`-chi2` を渡して呼ぶ (符号は反転)。設計軌道の幾何は、保存された CHI1、CHI2、DZ で決まり、DX、DY、ROTATE は現在の値を使う (ヘルプ)。

### 5.1 `tsolrot` の式と (8.10) の対応 (導出)

$B_z = 0$ の入口の式を、コードの変数のまま書く。`s0s` $= -\mathtt{alg}$、`dcchi1` $= \cos\mathrm{CHI}1 - 1$、`dcchi2` $= \cos\mathrm{CHI}2 - 1$ ($\mathtt{xsincos}$ の第 5 出力は $\cos x - 1$。`tfloor.f:260-282`)。$\mathtt{alg}$ は、`tfgeo.f:213-216` で $L/2$ (設計軌道に沿った長さの半分。BEND・MULT の `p_LGEO_*`) が入る量で、`tmulti.f:305` の出口では `alg-al` $= -L/2$ が渡る。

位置ベクトルを、要素の中心を原点として $v_0 = (x - \mathrm{DX},\ y - \mathrm{DY},\ \mathtt{s0s} - \mathrm{DZ})$ と置く (第 3 成分 $= -L/2 - \mathrm{DZ}$: 入口面にいる粒子は、中心から $L/2$ 手前、中心が DZ だけ下流にずれていれば $L/2 + \mathrm{DZ}$ 手前)。$R_1v_0$、続いて $R_2(R_1v_0)$ を (8.10) の行列で計算すると

$$
R_1v_0 = \begin{pmatrix}c_1(x - \mathrm{DX}) - s_1(\mathtt{s0s} - \mathrm{DZ})\\ y - \mathrm{DY}\\ s_1(x - \mathrm{DX}) + c_1(\mathtt{s0s} - \mathrm{DZ})\end{pmatrix},\qquad
R_2R_1v_0 = \begin{pmatrix}(R_1v_0)_x\\ c_2(y - \mathrm{DY}) - s_2(R_1v_0)_s\\ s_2(y - \mathrm{DY}) + c_2(R_1v_0)_s\end{pmatrix}
$$

となる。コードは、第 1 成分を `x(i)=cchi1*x-schi1*(s0s-dz)` (一致)、第 3 成分を中心からの値 $(R_1v_0)_s$ でなく原点 (面の位置 `s0s`) からの変位 `ds1` $= (R_1v_0)_s - \mathtt{s0s}$ $= s_1x' + (c_1 - 1)\mathtt{s0s} - c_1\mathrm{DZ}$ (`schi1*x+dcchi1*s0s-cchi1*dz`、$x' = x - \mathrm{DX}$) として持ち、第 2 成分を `y=cchi2*y-schi2*(ds1+s0s)` $= c_2y' - s_2(R_1v_0)_s$ (一致)、`ds2=schi2*y+cchi2*ds1+dcchi2*s0s` $= (R_2R_1v_0)_s - \mathtt{s0s}$ とする。sympy で、`x`・`y`・`ds2` が $R_2R_1v_0$ の成分から $\mathtt{s0s}$ を引いたものに一致することを確認した ([mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (d))。

運動量は、$p_z = \sqrt{1 - p_x^2 - p_y^2}$ (`pz0`) として、方向ベクトル $(p_x, p_y, p_z)$ に同じ $R_2R_1$ をかける ($p_x' = c_1p_x - s_1p_z$、$p_{z1} = s_1p_x + c_1p_z$、$p_y' = c_2p_y - s_2p_{z1}$。コードの `px`・`pz1`・`py` と一致、同 (d))。回転後の位置は、粒子の面 (要素の入口面、回転後の座標で第 3 成分 $= \mathtt{s0s}$) の上にはなく、第 3 成分が $\mathtt{s0s} + \mathtt{ds2}$ にある。そこで、新しい方向 $(p_x', p_y', p_{z2})$ に沿って**面まで戻す**: 戻す距離 (第 3 成分の差 `ds2` を $p_{z2}$ で割った経路長) が $a = \mathtt{ds2}/p_{z2}$ で、$x \to x - a\,p_x'$、$y \to y - a\,p_y'$、$z \to z + a$。$z$ の符号は [01 章](01_variables_and_conventions.md) の約束 $z = -v\,(t - t_0)$ (先頭が $z > 0$) から決まる: 点を直線に沿って $a$ ($a > 0$ は面を越えている) だけ**後ろ**へ戻すと、その粒子が面を通る時刻は $a/v$ だけ**早い**ので、$z = -v(t - t_0)$ は $+a$ だけ増える。逆に前へ進めれば $z$ は減る。この「$-a\,(\text{位置の移動量})$ だけ $z$ が変わる」関係は、出口でも同じ (次の段)。

**出口** (`ent = .false.`、$\mathtt{s0s} = \mathtt{alg} - L = -L/2$。要素の中心を原点にした出口の面の第 3 成分は $-\mathtt{s0s} = +L/2$) は、(8.10) の逆変換で、次の 3 行になる (ここでは CHI3 = 0 とする。CHI3 を含めると第 1 行の先頭に $R_3^{-1}$ が付く)。入口の順序 (平行移動 → $R_1$ → $R_2$ → 面へ → CHI3) の逆順で、各段をその逆に置き換える。

1. **回転 (CHI3 → $R_1^{-1}R_2^{-1}$)**。要素の出口の面にいる粒子の位置 $v = (x, y, +L/2)$ と方向 $\vec p = (p_x, p_y, p_z)$ に ($R_3(\mathrm{CHI}3)^{-1}$ のあと) $(R_2R_1)^{-1} = R_1^{-1}R_2^{-1} = R_1^{\top}R_2^{\top}$ をかける: $v' = R_1^{\top}R_2^{\top}v$、$\vec p' = R_1^{\top}R_2^{\top}\vec p$。成分では $p_y' $ の中間量が `pz1=-schi2*py+cchi2*pz0`、$p_y^{\rm out} = c_2p_y + s_2p_z$ のように、入口の式の $s_i \to -s_i$ に当たる (入口が $R_i$、出口が $R_i^{\top}$)。
2. **実験室の出口の面まで**。回転した位置 $v'$ の第 3 成分は、実験室の面 (中心が DZ だけずれているので $v'_s + \mathrm{DZ} = +L/2$) からずれている。第 3 成分の残り $a = (v'_s + \mathrm{DZ} - L/2)/p'_z$ を、方向 $\vec p'$ に沿って戻す: $x \to v'_x - a\,p'_x$、$y \to v'_y - a\,p'_y$、$z \to z + a$ (符号は入口と同じ。上の約束から)。
3. **平行移動を戻す**: $x \to x + \mathrm{DX}$、$y \to y + \mathrm{DY}$ (入口の $x - \mathrm{DX}$、$y - \mathrm{DY}$ の逆)。運動量は第 1 行の $p'_x$、$p'_y$ のまま。

この 3 行を書き下した写像は、K64 の `tsolrot` の出口の枝を書き写したものと、ランダムな DX・DY・DZ・CHI1・CHI2 の 6 組で $1.8\times10^{-17}$ 以内で一致する ([mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (d))。**入口と出口が正確に逆変換であること**を、空の要素 (中は長さ $L$ の直進だけ) で調べた: 入口 → 要素内の直進 $L$ → 出口が、実験室での長さ $L$ の直進に $10^{-17}$ 以下で一致する ((DX, DY, DZ, CHI1, CHI2) = (0, 0, 0, 0.01, 0.02)、(1 mm, 2 mm, 3 mm, 0.01, 0.02)、(1 mm, 2 mm, 3 mm, −0.03, 0.05) など 5 組、同 (d))。これで、**回転の支点は要素の中心 ($\mathtt{alg} = L/2$) で、DZ の符号は ($v_0$ の第 3 成分の) 要素の中心を下流へずらす向き**、と確定した。

**CHI1・CHI2 の向き**。ヘルプの (8.10) の $R$ を、座標の変換 (実験室 → 要素) と文字どおりに読むと、$R_1^{T}$ の第 3 列 ($R_1$ の第 3 行) から、要素の進行方向の単位ベクトルは実験室で $(\sin\mathrm{CHI}1,\ 0,\ \cos\mathrm{CHI}1)$ (CHI1 > 0 は要素の軸が $+x$ 側へ向く)。さらに $R_2$ を合わせると $R_1^TR_2^T(0, 0, 1)^T = (s_1c_2,\ s_2,\ c_1c_2)$ (CHI2 > 0 は $+y$ 側)。MULT は `tsolrot` に $-\mathrm{CHI}1$、$-\mathrm{CHI}2$ を渡す (`tmulti.f:93, 305, 379, 661`) ので、**MULT の実際の向きは、この文字どおりの読みの逆** ($(-s_1c_2,\ -s_2,\ c_1c_2)$。$R$ の第 3 列、すなわち $R$ を要素の向きを与える行列 (要素 → 実験室) と読んだ向き)。この向きは **gs の追跡で確認した** (下の「gs で確認した CHI1・CHI2 の符号」)。ヘルプの記述の読み方 (受動/能動) の違いで、K64 の誤りではない。この符号 ($\mathrm{tsolrot}$ への $-\mathrm{CHI}$) の理由を、コードから次のように決めた ([mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (e))。ANGLE $\neq 0$ の MULT は `tmulta` で、CHI2 を `tbrot` (BEND 系の回転。`dchi2 = chi2`、符号反転なし) に渡す。$\phi_g = 0$、$\Delta\theta = 0$ の入口で `tbrot` を書き下すと $p_y' = c_2p_y + s_2p_z$、$y' = c_2y - s_2\,\mathtt{alg}$ で、これは `tsolrot` に $-\mathrm{CHI}2$ を渡したときの式 ($p_y' = c_2p_y - s_2^{(ts)}p_{z1}$、$s_2^{(ts)} = \sin(-\mathrm{CHI}2)$) と $5\times10^{-20}$ で一致する (+CHI2 を渡すと $0.11$ 違う)。つまり **$-\mathrm{CHI}$ は、同じ MULT の ANGLE = 0 の経路 (`tsolrot`) と ANGLE $\neq 0$ の経路 (`tbrot`) で CHI2 の向きを同じにするため** (BEND 系にそろえるため。SOL は `tsolrot` に $+$CHI を渡すので、MULT と SOL は逆向き)。CHI1 は `tmulta` が使わない (§6.4) ので、独立な対応物が無く、$-\mathrm{CHI}1$ も CHI2 と同じ呼び出しで同時に反転されていることしか確認できない。**未確認**: DZ・ROTATE (CHI3) の符号の gs での確認。pysad の追跡は CHI1・CHI2・DZ を実装していないので、pysad では確かめられない。

**確認したこと** (追跡): DX・DY は、$(x - \mathrm{DX}, y - \mathrm{DY})$ で通して $\mathrm{DX}$、$\mathrm{DY}$ を足すのと同じ (差 0)。ROTATE は、入口で座標を $+\theta$ 回し出口で戻すのと同じ (差 0、逆向きは $6 \times 10^{-4}$)。

**gs で確認した CHI1・CHI2 の符号** (docs/30 の B337 (F)。試験 `test_v250_tsolrot_tilt.py`): gs の MULT ($L$ = 1 m、K1 = 0.3、3 GeV/c) の軸上の粒子の出口は、CHI1 = +0.05 で $x = +1.1939\times10^{-3}$・$p_x = -1.8396\times10^{-4}$、CHI2 = +0.05 で $y = -1.3064\times10^{-3}$・$p_y = -1.9141\times10^{-4}$ (CHI の符号を変えると符号が反転)。要素の軸を $(-\sin\chi_1,\ 0,\ \cos\chi_1)$ (CHI1)、$(0,\ -\sin\chi_2,\ \cos\chi_2)$ (CHI2) (上の $R$ の第 3 列) とし、その軸のまわりの厚い四極 (水平収束・垂直発散) を ODE で解いた値 $+1.1954\times10^{-3}$・$-1.8393\times10^{-4}$・$-1.3080\times10^{-3}$・$-1.9144\times10^{-4}$ と 0.13% 以内で一致する (符号まで。残りは傾いた入口・出口の面の幾何の違い)。つまり MULT の CHI1・CHI2 の符号は、(8.10) の $R_1$・$R_2$ を要素の向きとして読んだものと同じで、誤りではない (取り下げ)。`tmulti` が $-\mathrm{chi}$ を `tsolrot` に渡すのは内部の規約。限界: $L$ = 1 m の四極 1 つ・CHI = ±0.05。DZ・CHI3 の符号、傾いた領域の中の MULT (SOL 領域内) の本体は見ていない。SOL 領域内の傾いた枝の更新の符号には K64 の誤りがある ([09 章](09_solenoid.md) §8.3、B337。MULT の CHI の符号の問題とは別)。

## 6. ANGLE ≠ 0 の $A_s$

### 6.1 ヘルプの定義

ヘルプ (`multipole_with_nonzero_ANGLE`) は、ANGLE ≠ 0 の多極を、$H$ の中の項 $A_s(x, y)$ として

$$
A_s = \sum_{k,n=0}^\infty g_{kn}\,\frac{\mathrm{K}n + i\,\mathrm{SK}n}{(n+1)!}\,(\rho + x)^{1/2-k}\,\frac{(x + iy)^{n+k}}{\sqrt\rho},\quad
g_{kn} = -\frac{(2k-1)!!\,(2k-3)!!\,(n+1)!}{8^k\,(n+k+1)!\,k!},\quad\rho = \frac{L}{\mathrm{ANGLE}}
\tag{8.11}
$$

と書く (原文。$n + k \le 21$ で打ち切る。「K0 と K1 は BEND のものと違う」。実部を取ることは明記されていない)。この式は、そのままでは**正しくない** (§6.3)。この節で、$A_s$ を Maxwell 方程式から導く。

### 6.2 導出

曲線座標 $(x, y, s)$ ($h = 1/\rho$、スケール因子 $1, 1, 1 + hx$) で、$s$ によらない場が $\vec A = a_s\hat s$ だけで書けるとする。$H$ の中の項は $-(1 + hx)a_s \equiv A_s$ ([01 章](01_variables_and_conventions.md) 式 (1.1))。$F \equiv (1 + hx)\,a_s$ とおく。直交曲線座標 (スケール因子 $h_x = h_y = 1$、$h_s = 1 + hx$) の回転は $(\nabla\times\vec A)_x = \frac{1}{h_yh_s}\bigl[\partial_y(h_sA_s) - \partial_s(h_yA_y)\bigr]$、$(\nabla\times\vec A)_y = \frac{1}{h_sh_x}\bigl[\partial_s(h_xA_x) - \partial_x(h_sA_s)\bigr]$、$(\nabla\times\vec A)_s = \frac{1}{h_xh_y}\bigl[\partial_x(h_yA_y) - \partial_y(h_xA_x)\bigr]$ で、$A_s = a_s$ だけ、$s$ によらない場では $A_x = A_y = 0$、$\partial_s = 0$ なので、$\vec B = \nabla\times\vec A$ の成分は ($h_sA_s = F$ を使って)

$$
B_x = \partial_ya_s = \frac{\partial_yF}{1 + hx},\qquad B_y = -\frac{1}{1 + hx}\,\partial_x\bigl((1 + hx)a_s\bigr) = -\frac{\partial_xF}{1 + hx}
\tag{8.12}
$$

($B_s = 0$。$\nabla\cdot\vec B = (1 + hx)^{-1}\bigl[\partial_x((1 + hx)B_x) + \partial_y((1 + hx)B_y)\bigr]$ で、$(1 + hx)B_x = \partial_yF$、$(1 + hx)B_y = -\partial_xF$ だから $(1 + hx)^{-1}(\partial_x\partial_yF - \partial_y\partial_xF) = 0$。恒等的に満たされる。) **場が真空の Maxwell 方程式を満たす条件**は $(\nabla\times\vec B)_s = \partial_xB_y - \partial_yB_x = 0$ で、(8.12) を入れて

$$
\partial_x\Bigl(\frac{\partial_xF}{1 + hx}\Bigr) + \frac{\partial_y^2F}{1 + hx} = 0\ \ \Longleftrightarrow\ \ \nabla^2F = \frac{\partial_xF}{\rho + x}
\tag{8.13}
$$

(第 2 式への変形: 第 1 式の第 1 項を積の微分で展開すると $\partial_x\bigl(\tfrac{F_x}{1 + hx}\bigr) = \tfrac{F_{xx}}{1 + hx} - \tfrac{h\,F_x}{(1 + hx)^2}$ ($\partial_x(1 + hx)^{-1} = -h(1 + hx)^{-2}$)。第 1 式全体は $\tfrac{1}{1 + hx}\bigl[F_{xx} + F_{yy} - \tfrac{h\,F_x}{1 + hx}\bigr] = 0$。$1 + hx$ をかけて $\nabla^2F = \tfrac{h}{1 + hx}F_x = \tfrac{F_x}{\rho + x}$ ($h = 1/\rho$ から $h/(1 + hx) = 1/(\rho + x)$)。$\nabla^2 = \partial_x^2 + \partial_y^2$。) $\rho \to \infty$ では $\nabla^2F = 0$ (調和関数 $\mathrm{Re}\,z^{n+1}$)。

**$\rho$ の冪の展開**。調和関数 $z^m$ ($m = n + 1$) を出発点に、$F = \sum_{k\ge0}g_k\,u^{a_k}z^{m+k}$、$u = \rho + x$、$a_k = \tfrac12 - k$ と置く。1 項ずつ計算する: $z^m$ は調和 ($\nabla^2z^m = 0$)、$\partial_x u^a = au^{a-1}$、$\partial_yu = 0$ だから

$$
\nabla^2(u^az^m) = a(a-1)\,u^{a-2}z^m + 2am\,u^{a-1}z^{m-1},\qquad
\frac{\partial_x(u^az^m)}{u} = a\,u^{a-2}z^m + m\,u^{a-1}z^{m-1}
$$

(第 1 式の第 2 項は $2\,\partial_xu^a\,\partial_xz^m = 2au^{a-1}mz^{m-1}$。) (8.13) の左辺 − 右辺 ($\nabla^2 - \partial_x/u$) に入れると、同じ型の項をまとめて

$$
\bigl[a(a - 1) - a\bigr]\,u^{a-2}z^m + \bigl[2am - m\bigr]\,u^{a-1}z^{m-1} = a(a-2)\,u^{a-2}z^m + m(2a - 1)\,u^{a-1}z^{m-1}
$$

(sympy で、$\bigl[\nabla^2 - \partial_x/u\bigr](u^az^m)/(u^{a-2}z^{m-1}) = a(a-2)\,z + m(2a-1)\,u$ と確認、[mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) (e)。$u = \rho + x$ なので $u^{a-2}z^m = u^{a-2}z^{m-1}\cdot z$)。$u^{a-2}z^m$ の係数が $a(a-2)$、$u^{a-1}z^{m-1}$ の係数が $m(2a - 1)$。ここで $F = \sum_k g_ku^{a_k}z^{m+k}$ ($a_k = \tfrac12 - k$) に適用する。$k$ 番目の項 ($a = a_k$、$m \to m + k$) の第 1 の型は $g_k\,a_k(a_k - 2)\,u^{a_k-2}z^{m+k}$、第 2 の型は $g_k\,(m + k)(2a_k - 1)\,u^{a_k-1}z^{m+k-1}$。$(k+1)$ 番目の項の第 2 の型は $g_{k+1}(m + k + 1)(2a_{k+1} - 1)\,u^{a_{k+1}-1}z^{m+k}$ で、$a_{k+1} - 1 = a_k - 2$ なので $k$ 番目の第 1 の型と同じ形 ($u^{a_k-2}z^{m+k}$) になる。全体が 0 であるためには、各形の係数の和が 0:

$$
g_k\,a_k(a_k - 2) + g_{k+1}\,(m + k + 1)(2a_{k+1} - 1) = 0
$$

(最も低い形 $u^{a_0-1}z^{m-1}$ は、$k = 0$ の項の第 2 の型 $g_0\,m(2a_0 - 1)$ だけで、これに釣り合う項が無い。したがって $m(2a_0 - 1) = 0$、**$a_0 = \tfrac12$**。これが $A_s$ の主項が $\sqrt{\rho + x}$ をかけた形になる理由で、「$a_k = \tfrac12 - k$」は、この要請と、$k$ が 1 つ増えるごとに $u$ の冪が 1 下がる (漸化式の型の一致) の 2 つから決まる。) $a_k(a_k - 2) = (\tfrac12 - k)(-\tfrac32 - k) = \tfrac14(2k - 1)(2k + 3)$、$2a_{k+1} - 1 = -2(k + 1)$、$m = n + 1$ を入れると

$$
g_{k+1} = g_k\,\frac{(2k - 1)(2k + 3)/4}{2(k + 1)(n + k + 2)} = g_k\,\frac{(2k - 1)(2k + 3)}{8\,(k + 1)(k + n + 2)}
$$

(sympy で同じ式を確認、同 (e))。$k + 1 \to k$ と読み替えると

$$
\boxed{\;g_{kn} = g_{k-1,n}\,\frac{(2k - 3)(2k + 1)}{8\,k\,(k + n + 1)},\qquad g_{0n} = 1\;}
\tag{8.14}
$$

($m = n + 1$ を使った。$g_{0n} = 1$ は、$\rho \to \infty$ で直線の多極 (8.1) に戻る規格化。) 閉じた形は、(8.14) を $g_{0n} = 1$ から $k$ 回かけた積

$$
g_{kn} = \prod_{j=1}^{k}\frac{(2j - 3)(2j + 1)}{8\,j\,(j + n + 1)}
$$

の各因子の積から出る: $\prod_{j=1}^{k}(2j-3) = (-1)\cdot1\cdot3\cdots(2k-3) = -(2k-3)!!$ ($j = 1$ の因子が $-1$、以後は正の奇数。$k = 1$ は $-1 = -(-1)!!$ で $(-1)!! = 1$ の約束)、$\prod_{j=1}^{k}(2j+1) = 3\cdot5\cdots(2k+1) = (2k+1)!!$、$\prod_{j=1}^{k}8j = 8^kk!$、$\prod_{j=1}^{k}(j + n + 1) = (n+2)(n+3)\cdots(n+k+1) = (n+k+1)!/(n+1)!$ なので

$$
g_{kn} = -\frac{(2k+1)!!\,(2k-3)!!\,(n+1)!}{8^k\,(n+k+1)!\,k!}\qquad(k \ge 1)
$$

($n, k \le 7$ で積と閉じた形が一致、同 (e))。最初の値 (同 (e2)): $n = 0$: $g_{0..3} = 1,\ -3/16,\ -5/256,\ -35/8192$。$n = 1$: $1,\ -1/8,\ -5/512,\ -7/4096$。$n = 2$: $1,\ -3/32,\ -3/512,\ -7/8192$。(8.14) は K64 の `gkninit` (`tffs.f:1188-1200`) と**同じ**。結局

$$
A_s = \sum_{k,n}g_{kn}\,\mathrm{Re}\,\frac{(\mathrm{K}n + i\,\mathrm{SK}n)\,z^{n+k+1}}{(n+1)!}\,\frac{(\rho + x)^{1/2-k}}{\sqrt\rho}
\tag{8.15}
$$

($\sqrt\rho$ で割る因子は、$\rho \to \infty$ で $k = 0$ の項 $u^{1/2}/\sqrt\rho = \sqrt{1 + x/\rho}$ が 1 になる規格化。$k = 0$ の項が $\sqrt{1 + x/\rho}$ を持つのは、上の $a_0 = \tfrac12$ の要請による。) $k$ の項は $u^{1/2-k}/\sqrt\rho = \rho^{-k}(1 + x/\rho)^{1/2-k}$ なので $(1/\rho)^k$ で小さくなる ($\rho$ についての展開)。K64 は $n + k \le 21$ で打ち切る (`nmult` = 21)。

**キック (`tmulta`)**。$\Delta p_x = -\partial_xA_s$、$\Delta p_y = -\partial_yA_s$。$\partial_x[u^{a_k}z^{n+k+1}] = u^{a_k}(n + k + 1)z^{n+k} + a_k\,u^{a_k-1}z^{n+k+1}$、$\partial_y[u^{a_k}z^{m}] = i\,m\,u^{a_k}z^{m-1}$ (sympy で確認、同 (e))。$A_s = \mathrm{Re}\sum c\,u^{a_k}z^{n+k+1}$ の係数 $c$ ($\mathrm{K}n + i\,\mathrm{SK}n$ を含む複素数) は $x, y$ に依らないので、微分は $\mathrm{Re}$ の中に入る: $\partial_xA_s = \mathrm{Re}\sum c\bigl[a_k\,u^{a_k-1}z^{n+k+1} + (n + k + 1)\,u^{a_k}z^{n+k}\bigr]$、$\partial_yA_s = \mathrm{Re}\sum c\,i\,(n + k + 1)\,u^{a_k}z^{n+k} = -\mathrm{Im}\sum c\,(n + k + 1)\,u^{a_k}z^{n+k}$。$\Delta p = -\nabla A_s$ なので、$u^{a_k-1}z^{n+k+1} = (u^{a_k}z^{n+k})\,z/u$ を使って

$$
\Delta p_x = -\frac{\mathrm{Re}\,[\,C_s\,z\,]/(\rho + x) + \mathrm{Re}\,C_l}{1 + \delta},\qquad
\Delta p_y = +\frac{\mathrm{Im}\,C_l}{1 + \delta},\qquad
C_l = \sum_{n,k}(n + k + 1)\,c_{kn}\,z^{n+k},\ \ C_s = \sum_{n,k}\bigl(\tfrac12 - k\bigr)\,c_{kn}\,z^{n+k}
\tag{8.16}
$$

($c_{kn} \equiv g_{kn}\,r_k\,\mathrm{akn}_n$、$r_k = (\rho + x)^{1/2-k}/\sqrt\rho$、$\mathrm{akn}_n = (\mathrm{K}n + i\,\mathrm{SK}n)/((n+1)!\,\mathrm{ndiv})$。`csl`=$C_l$、`csr`=$C_s$。`tmulta.f:170-200`: `rk(0)=1+sqrt1(x(i)/rho0)` ($= \sqrt{1 + x/\rho_0}$)、`rk(k)=rk(k-1)/r`、`px(i)=px(i)-(dble(csr*cx1)/r+dble(csl))/pr`、`py(i)=py(i)+imag(csl)/pr`。)

### 6.3 ヘルプの式との食い違い (確定)

(8.11) と (8.14)(8.15) を比べる ([mult_01_as_expansion.out.txt](outputs/mult_01_as_expansion.out.txt)):

- **$g_{kn}$ の因子**: 原文は $(2k-1)!!$、導出と K64 は $(2k+1)!!$。$n = 0$ で $k = 1$: 原文 $-1/16$、導出 $-3/16$。K64 の漸化式の閉じた形が $(2k+1)!!$ であることを $n, k \le 7$ で確認した。
- **$z$ の冪**: 原文は $(x + iy)^{n+k}$ (ポテンシャルの冪が $n + k$ になる)。$k = 0$ で直線の極限が $z^n$ で、K$n$ のポテンシャルの冪 $n + 1$ ($z^{n+1}/(n+1)!$) と合わない。(8.15) は $n + k + 1$。
- **検証**: (8.13) の左辺 − 右辺に、$k \le 3$ で打ち切った (8.15) を入れ、$\rho = 1$ で $x, y$ の冪に展開すると、K64 の $g_{kn}$ と冪 $n + k + 1$ の組だけが、打ち切りの次数まで 0 (残差の係数が見つからない)。原文の $g_{kn}$ や冪 $n + k$ では、低次から 0 でない残差が残る (例えば $n = 1$、冪 $n + k + 1$ で原文の $g$ は 2 次)。
- ヘルプの「実部を取る」の省略と、$n + k \le 21$ の打ち切り (原文で確認) は、K64 の `nmult` = 21、`gkn(n,k)` の表 ($k \le$ `nmult` − $n$) と同じ。

### 6.4 ANGLE ≠ 0 の追跡の構成 (`tmulta`)

`tmulti` は、ANGLE ($\texttt{phia}$) ≠ 0 のとき `tmulta` を呼ぶ (`tmulti.f:78-88`。$B_z \ne 0$ は `abort`)。

1. ずれ: `tbshift` (DX、DY、ROTATE)、`tbrot` (DROTATE、CHI2)。**CHI1・DZ は ANGLE ≠ 0 では使わない** (`tmulta` の引数に無い)。
2. `ak = ak0*rtaper`。K0 と K1 の実部は**偏向磁石の本体**に回し (`ak0r`、`ak1r`)、`ak(0)`・`ak(1)` からは実部を除く (`tmulta.f:50-55`)。つまり、K0・K1 の実部は、BEND のように厳密な円弧 (K1 ≠ 0 は複合機能の本体) で解き、スキュー成分 (SK0、SK1) と K2 以上の成分だけを $A_s$ (8.15) のキックで入れる。**K0 と K1 のキックの級数 ($k \ge 1$ の項) は入らない**。ヘルプの「K0 と K1 は BEND のものと違う」は、K0 は ANGLE の上の追加の偏向 ($\mathrm{ANGLE} + \mathrm{K0}$ が全体の偏向)、K1 は BEND の K1 と同じ暗黙の定義 ($\mathrm{B160}$) であることを指すと読める (未確認)。
3. 刻み: $n_{\min}$〜$n_{\max}$ (`ak` の 0 でない最小・最大の次数) について $\mathrm{ndiv} = 1 + \lfloor\max_n\sqrt{a_{\max}^{\,n-1}|\mathrm{ak}(n)L|/(6(n-1)!\,\varepsilon)}\rfloor$ ((8.8) と同じ形、`tmulta.f:69-78`。$n$ の下限が 0 からで、$n = 0$ では `fact(-1)` が要る。SK0 > 0 の MULT が pysad で未対応な理由)。上限 2000。
4. 各スライス: `tbend` (半分、入口の端面と Maxwell フリンジ込み) → [キック (8.16) → `tbend` (全体)] を $\mathrm{ndiv} - 1$ 回 → キック → `tbend` (半分、出口の端面)。BEND の $L'$ の補正 ([05 章](05_bend.md) 式 (5.43)) は、`aln` から $(\phi\,\mathrm{fb})^2/(48L)\,\sin(\ldots)/\sin(\phi/2)/\mathrm{ndiv}$ ずつ引く形 (`tmulta.f:82-90`)。**同じ補正が `tffs.f:1366-1371` でも `p_L_MULT` に入っている**。呼び出しの連鎖を追うと、`tturn.f:761` が `tmulti` に `cmp%value(p_L_MULT)` ($= L'$) を渡し、`tmulti.f:79` が `tmulta` にその `al` をそのまま渡す。つまり**追跡の `tmulta` は補正済みの $L'$ を受け取り、さらに引く**。1 回目 ($L' = L - c_1$、$c_1 = ((\phi\,\mathrm{fb}_1)^2 + (\phi\,\mathrm{fb}_2)^2)\,\sin(w/2)/(48\,L\sin(\phi/2))$) に続いて、2 回目は $c_2 = (\ldots)/(48L')\ldots$ (分母が $L'$) を、$\mathrm{ndiv}$ 個のスライスの $\mathtt{aln}$ から $c_2/\mathrm{ndiv}$ ずつ引くので、円弧の合計は $L' - c_2 = L - c_1 - c_2$。pysad の k64 は同じ構造 (`multa.compile_multa` が $L'$ を作り、`multa.k_multa` がさらに引く)。v2 は 2 回目を省く。数値例 ($L = 1$ m、ANGLE = 0.2、FB1 = FB2 = 0.2、E1 = E2 = 0): $c_1 = 6.667 \times 10^{-5}$、$c_2 = 6.667 \times 10^{-5}$、円弧の合計は $0.99986666$ ($L' = 0.99993333$ から $6.667 \times 10^{-5}$ 短い)。`compile_multa` が渡す `al` が $L'$ に一致することも確認した ([mult_03_Lprime.out.txt](outputs/mult_03_Lprime.out.txt) (1)(2))。**K64 の誤り。gs で確認、v2 は修正 = B331**。gs で、3 GeV/c、$L$ = 1 m、ANGLE = 0.3、E1 = 0.5、E2 = 0.3、F1 (FB1 = FB2) = 0.5、FRINGE = 3 の BEND と MULT の `TrackParticles` は $x$ で $2.8\times10^{-7}$、$z$ で $2.5\times10^{-8}$ 違う (ANGLE = 0.1・F1 = 0.2 でも $\delta = 10^{-3}$ の $x$ が $8.3\times10^{-10}$ 違い、二重の分の $L$ の短縮 $1.67\times10^{-5}$ m × 分散の予想と一致)。MULT の $L$ を二重の補正の分 ($L - c(L) - c(L - c(L)) = L'_{\rm BEND}$、$L$ = 1.00018814134845) だけ先に伸ばすと、全成分で BEND と差 0 (ビット単位) になる。F1 なしでは差 0。pysad の k64 は gs と $10^{-15}$ 以内で一致。**v2 は 2 回目の補正を省く** (MULT が BEND と全桁一致。試験 `test_v250_mult_angle_fb_double.py`)。限界: FRINGE = 3・K1 = 0 で確認。BEND と MULT の FRINGE = 1・2 の $x$ の $3\times10^{-5}$ 級の食い違いは仕様の違い (BEND は FRINGE ≠ 0 で F1 を両端に使い、MULT の FRINGE = 1・2 は片側。対にすれば差は $L'$ の二重補正の分だけ)。K1 ≠ 0 の $y$ の食い違いは別の誤り (§6.5、B340)。CALC `tmultae` も同じ式を持つが gs では確かめていない。以下の ODE との比較は、この二重補正が物理的にも過剰であることを、gs と独立に確かめたもの。

**ODE との比較で確定できたこと** ([mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) (f))。条件: 扇形 (E1 = E2 = 0)、ANGLE = 0.2、$L = 1$ m、FB1 = FB2 = 0.2 m (線形の立ち上がり、$\mathrm{F1} = 0.2$ m)、平面 ($y = 0$)。実験室 $(X, Z)$ で、入口の面からの距離 $s_{\rm in} = Z$ と、出口の面からの距離 $s_{\rm out}$ の、2 つの線形の立ち上がり $g_{\rm in}(s_{\rm in})g_{\rm in}(-s_{\rm out})/\rho'$ の磁場 (磁場の円の半径 $\rho' = L'/\phi$、$L'$ は (5.43) の 1 回だけの補正) の中で、$P' = P\times b$ を DOP853 で積分する ($y = 0$ の中央面では Maxwell フリンジの非線形項は 0)。設計粒子は出口で $x = 7\times10^{-10}$、$p_x = 1.5\times10^{-9}$ に出る (05 章の $L'$ と合う)。$\delta = \pm0.01$ の粒子の出口 $x$ (出口の面の座標) は

| | ODE | MULT | BEND | MULT − ODE | BEND − ODE |
|---|---|---|---|---|---|
| $\delta = +0.01$ | $9.8678\times10^{-4}$ | $9.8671\times10^{-4}$ | $9.8677\times10^{-4}$ | $-7.6\times10^{-8}$ | $-1.0\times10^{-8}$ |
| $\delta = -0.01$ | $-1.00651\times10^{-3}$ | $-1.00644\times10^{-3}$ | $-1.00650\times10^{-3}$ | $+7.3\times10^{-8}$ | $+5.4\times10^{-9}$ |

$\delta$ の奇数部分 ($x$ の分散 $R_{16}$ の誤差) は MULT − ODE $= -7.4\times10^{-6}$、BEND − ODE $= -7.9\times10^{-7}$、MULT − BEND $= -6.65\times10^{-6}$。一方、2 回目の補正 $c_2 = 6.667\times10^{-5}$ m を、角を変えずに円弧の長さから引くと、半径が $c_2/\phi$ だけ小さくなり $R_{16} = \rho(1 - \cos\phi)$ が $-c_2(1 - \cos\phi)/\phi = -6.645\times10^{-6}$ 変わる。**MULT と BEND の差は、この予想と 0.2% で一致する**。BEND (補正は 1 回) は ODE と $10^{-8}$ 以内 (残り $7.9\times10^{-7}$ は、$F1$ の 1 次の近似と、この ODE の違いで、未分離)、MULT は ODE から 10 倍離れる。したがって、**この条件では、`tmulta` の 2 回目の補正は物理的に過剰 (二重補正) で、円弧が $c_2$ 短すぎることによる $R_{16}$ の誤差 $-6.6\times10^{-6}$ を生む**。(意図かどうかは、ソースに説明が無いので決められないが、結果は ODE と合わない。pysad の k64 は K64 のまま、v2 は 2 回目の補正を省く。)

**確認していないこと**: (i) CALC 側 (`tmultae`、`tmulta.f:327-336`) の補正と、`qmult.f:28` が渡す `al` が $L$ か $L'$ か (pysad は FB1・FB2 を CALC で未対応)。(ii) $\psi \neq 0$ (矩形)、K2 $\neq 0$、$y \neq 0$、斜め入射 (gs の確認は E1 = 0.5・E2 = 0.3 を含む)。(iii) FB が小さい / 大きいときの $c_2$ の比例 (理論上 $\mathrm{FB}^2$)。(iv) $L'$ 自体の ODE 確認は 05 章 §6.4。

### 6.5 K1 ≠ 0 の MULT の端 (FB1・FB2) が途中にも入る (K64 の誤り。gs で確認、v2 は修正 = B340)

前半・後半の半スライスの `tbend` には (fb1, fb2) の両方が渡る (入口は mfr = -1、出口は mfr = -2。`tmulta.f:148,236`)。K1 = 0 の `tbend` は mfr で fb を選ぶ (`tbend.f:400,415`) ので正しいが、K1 ≠ 0 は `tbendi` に行き (`tbend.f:595`)、`tbendi` は mfr に関係なく fb1 を入口に・fb2 を出口に入れる (`tbendi.f:73,159`)。このため、MULT の前半の終わりに FB2、後半の始めに FB1 の端のキック ($x$ の $\delta$ 依存・$y$ の立方の項・$z$) が余計に入る (同じ磁場の BEND は端が両端の 1 回ずつ)。gs で、3 GeV/c、$L$ = 1 m、ANGLE = 0.3、K1 = 0.3、E1 = 0.5、E2 = 0.3、F1 (FB1 = FB2) = 0.5、NORAD・NORFSW の `TrackParticles` (7 粒子) を比べると、BEND と MULT は F1 = 0 で $2.2\times10^{-11}$、F1 ≠ 0 で両端 (BEND FRINGE = 3、MULT FRINGE = 3) が最大 $1.41\times10^{-4}$ ($y$・$p_y$)、入口のみ・出口のみで各 $7.0\times10^{-5}$ (途中の余計な端 1 つぶん) 違う (E1・E2 の有無、DISFRIN の 0・1 によらず同じ大きさ)。pysad の k64 は BEND・MULT とも gs と $10^{-15}$ 以内で一致する。**v2 は前半に fb2、後半に fb1 を渡さない**ので、MULT が BEND と $10^{-12}$ 以内で一致する (試験 `test_v250_mult_k1_fringe.py`、記録 `docs/history/B340_mult_k1_fringe.md`)。限界: NORAD の追跡だけ。RAD ON の `tbendi` の first・last の扱い、CALC (`tmultae` → `tbende`) は未確認。
5. 逆のずれ。

CALC (`tmultae`) は `tmulta` と同じ分割で行列を作る (`multae.py`)。キックを $(1+\delta)$ で割らないなど、追跡と違う点がある (読みだけ)。

## 7. 行列 (CALC・Emittance[])

`tmulte` (ANGLE = 0) は、追跡と同じ刻み数の式 (8.8) (`ampmax`・`eps00m`) と対称な分割 (8.7) を、COD のまわりで微分して行列にする (`tmulte.f:131-145`)。`tmulte` は、ずれ・回転に `tsolrote`、K0・K1 に `tsolque` を使う ([04 章](04_quad.md) §8、[09 章](09_solenoid.md))。この章では、CALC の MULT の行列と追跡の写像の比較はしていない。

## 8. 確認

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `mult_01_as_expansion.py` | (a) ヘルプと K64 の $g_{kn}$ の比較、(b) (8.13) の残差 (sympy。$k \le 3$ の打ち切り、冪 $n + k + 1$ と $n + k$)、(c) K64 の漸化式の閉じた形 | [mult_01_as_expansion.out.txt](outputs/mult_01_as_expansion.out.txt) |
| `mult_02_tracking.py` | (a) ANGLE = 0、K2: `ndiv` の式と誤差の次数 (2 次)。(b) K1 + K2。(c) SK$n$ = K$n$・ROTATE $90°/(n+1)$、DX・DY、ROTATE の向き。(d) ANGLE ≠ 0、K2: (8.15) から作った曲線座標の ODE との比較 | [mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) |
| `mult_03_Lprime.py` | ANGLE ≠ 0 の $L'$ の補正が 2 回入る構造 (§6.4 項 4): 補正量、`compile_multa` が渡す `al`、MULT と BEND の追跡の比較 | [mult_03_Lprime.out.txt](outputs/mult_03_Lprime.out.txt) |
| `mult_04_derivations.py` | 途中式の確認 (sympy・数値): (a) Horner 形 = $C(z)$ (b) 回転 (8.4) (c) 端面の角 (8.9): 解析性・ポテンシャル・K0・K1 の展開 (d) `tsolrot` = (8.10)、入口と出口が逆変換 (e) 漸化式 (8.14)・閉じた形・(8.16) の微分 (f) Strang の誤差 (調和振動子・BCH・ndiv の数値) | [mult_04_derivations.out.txt](outputs/mult_04_derivations.out.txt) |
| `mult_06_gaps.py` | (a) (8.8) の $\varepsilon$ の読み: 実効ハミルトニアン (8.8b) (BCH、sympy)・その厳密な流れと追跡の誤差の一致・調和振動子と 6 極の見積り。(b) `dofr(n)` の判定量と位置のずれの比 $4n!(n+2)$。(c) (8.9a) の rot・div と (8.9) との差。(d) 出口の 3 行と `tsolrot` の出口の一致。(e) `tbrot(+CHI2)` = `tsolrot(-CHI2)`。(f) ANGLE $\neq 0$ の $L'$ の二重補正: 平面 ODE と MULT・BEND の比較 | [mult_06_gaps.out.txt](outputs/mult_06_gaps.out.txt) |
| `mult_05_edge_v2_vs_k64.py` | 端面の角の kick: K64 の `tmulk` (8.9) と pysad v2 の厳密な面の幾何の比較 (K0 は一致、K1 は $p_x$ の $y^2$ 項が違う) | [mult_05_edge_v2_vs_k64.out.txt](outputs/mult_05_edge_v2_vs_k64.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/mult_01_as_expansion.py` (数十秒、sympy)、`py -3.12 physics_guide/scripts/mult_02_tracking.py` (数十秒)、`mult_03`・`mult_05` (数秒)、`mult_04_derivations.py` (数十秒、sympy)。

主な数字 (電子 3 GeV、$\delta = 0$): ANGLE = 0 で、EPS = 4 → 0.25 の `ndiv` 5 → 17 に対し、厳密解との差が $6.4 \times 10^{-7}$ → $5.5 \times 10^{-8}$ (2 次)。SK1・SK2・SK3 と、K$n$・ROTATE = 45°・30°・22.5° の差は $10^{-17}$ 以下。ANGLE = 0.2 ($\rho = 5$ m)、K2 = 40 /m²: 扇形の偏向だけの追跡と ODE の差が $4.8 \times 10^{-7}$ (硬い端の Maxwell フリンジ。ODE に入れていない)。K2 の効果 ($5.7 \times 10^{-4}$) を、K2 = 0 との差で比べると、追跡と ODE の差は EPS = 4 → 0.25 で $6.9 \times 10^{-7}$ → $2.3 \times 10^{-7}$ (相対 $4 \times 10^{-4}$)。この残りは、EPS を小さくしても $2 \times 10^{-7}$ に残る。K2 = 40 → 80 で EPS = 0.25 の残りは $2.25 \times 10^{-7}$ → $4.41 \times 10^{-7}$ (約 1.96 倍) と **K2 にほぼ比例する**ので、K2 の 2 乗の項ではなく、K2 の 1 次の項。K2 のキックで軌道が変わったとき、Maxwell フリンジ (K2 = 0 でも ODE との差 $4.8 \times 10^{-7}$ を生む項) の大きさが軌道に応じて変わる、K2 の 1 次の交差項と推測する (同じ桁。確認していない)。

**この比較が確かめること**。(d) の ODE は (8.15) から作った曲線座標の運動方程式なので、確かめられるのは「追跡 (`tmulta`) が (8.15) の展開を実装したものと合う」ことまで。(8.15) が Maxwell の解であることは、(8.13) の残差 ([mult_01_as_expansion.out.txt](outputs/mult_01_as_expansion.out.txt) (b)) と、§6.2 の導出で別に確認している。両方を合わせて、追跡が Maxwell の解を (打ち切りの次数まで) 実装している、と言える。

## 9. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、$\delta = 0$。ANGLE = 0 は $L$ = 1 m (確認 (a)(b))、0.3 m (確認 (c))、K2 = 40 /m²、初期値は mm・mrad の桁。ANGLE ≠ 0 は ANGLE = 0.2 rad、$L$ = 1 m、扇形 (E1 = E2 = 0)、K2 のみ (K0 = K1 = 0)、DISFRIN = 1。放射・スピン・ソレノイド・RF・ずれ (CHI1・CHI2・DZ) は含まない。周長・磁場の上限は使っていない。

**この章の範囲外**: ソレノイド中の MULT (`bzs ≠ 0`)、RF つきの MULT (`tmultiacc`。10 章)、放射 (`tradk`、`ndivrad`)、空間電荷 (`spkick`)、MULT の `tmulte` の行列と追跡の比較 (§7)。

**未確認・限界**:
- **CHI1・CHI2・DZ** (§5.1): `tsolrot` の式が (8.10) の $R_2R_1$ と一致すること、支点が要素の中心 ($L/2$)、入口と出口が (8.10) の逆変換であること (出口の 3 行)、MULT が $-\mathrm{CHI}$ を渡す理由 (ANGLE = 0 の `tsolrot` と ANGLE $\neq 0$ の `tbrot` の向きをそろえる。CHI2 で数値確認) は確認した。**CHI1・CHI2 の符号は gs で確認した** (§5。軸は $R$ の第 3 列の向き。K64 の誤りではない)。**未確認**は、DZ・ROTATE の符号と傾いた領域の中の MULT の gs での確認。pysad の追跡は CHI1・CHI2・DZ を実装していないので pysad の数値の確認ができない。
- **端面の角の式 (8.9)**: 中央面の値と調和性から導出した (§4.4。1 次・硬い端の仮定)。K0 の極限は標準の端面収束と一致。pysad v2 の式は K1 で $p_x$ の $y^2$ 項が違い、どちらが厳密かは確定できない。
- **刻み数の式 (8.8)**: $\varepsilon$ を「$a_{\max}$ での $n$ 極の運動量の跳びの相対誤差の上限」と読むと、係数 6、$Q^1$、$(n-1)!$ がすべて導出でき、実測 (位相誤差だけの見積りの約 1.5 倍) も位相の $1/24$ だけを数えた見積りの不足で説明できる (§4.3)。ただし $\varepsilon$ の読みは推定で、K64 の著者の定義は未確認。$\omega L \ll 1$ の近似なので、$a_{\max}$ ($\omega L = 1.4$) では目安。`dofr(n)` の判定量に $n!$ が無い理由も未確認 (閾値の読みは §4.4)。
- **ANGLE ≠ 0**: (8.15) は Maxwell 方程式の厳密な解 ($k$ について無限級数) の $n + k \le 21$ での打ち切り。K0・K1 の ($k \ge 1$ の) 項は入れず、本体 (BEND の円弧) で扱う。追跡との比較 (§8) は、Maxwell フリンジに起因する $\sim 2 \times 10^{-7}$ の残りまでで、これ以上は分離していない。複数の多極 ($n \ge 3$) や、スキュー (SK$n$)・K1 ≠ 0 の ANGLE ≠ 0 は、ODE との比較をしていない。
- **$L'$ の補正の二重適用** (§6.4): ソースの呼び出しの連鎖から追跡の `tmulta` が補正済みの $L'$ をさらに補正することを確認し、平面の ODE (扇形、線形の立ち上がり) との比較で、2 回目の補正が円弧を短くしすぎ ($R_{16}$ の誤差 $-6.6\times10^{-6}$。MULT と BEND の差と 0.2% で一致)、BEND (1 回) の方が ODE に合うことを確認した。**gs でも確認済み** (K64 の誤り、v2 は修正。B331)。FRINGE = 1・2 の BEND との食い違いは仕様の違い。未確認: CALC 側 (`tmultae`) の `al` が $L$ か $L'$ か。
- **K1 ≠ 0 の MULT の端の重複** (§6.5): gs で確認済み (K64 の誤り、v2 は修正。B340)。未確認: RAD ON、CALC。
- **`tmulk` の高次項** (§4.4): $n \ge 1$ の項は K64 の式と、pysad v2 の幾何厳密な式で水平方向が違う (上記)。F1・F2・FB の端の線形フリンジと (8.9) の重なり (二重に数えていないか) は調べていない。
- **ずれの CHI の ANGLE ≠ 0 版**: `tmulta` は CHI1・DZ を使わない (§6.4)。この場合の幾何は扱っていない。

## 付録 A. ヘルプの記述・K64 との違い (MULT)

原文 (`SADHelp.html` の `MULT`、`misalignments`、`multipole_with_nonzero_ANGLE`) で確認した。

| 箇所 | 内容 |
|---|---|
| MULT の説明 (reference plane) | **原文で確認**: 「the reference plane is defined so that the skew quadrupole component becomes zero」。§3.3 のとおり ($\theta_1$)。 |
| SK$n$ の説明 (「ROTATE = $90/(n+1)$ DEG」) | **原文で確認**、数値で一致 (§3.2)。ROTATE の値は $90/(n+1)$ 度 (K1 なら 45、K2 なら 30、…。K6 では 12.857…) と書かれている。「rotated around z-axis by $-90/(n+1)$ degree」は、座標を $+\theta$ 回すことに当たる。 |
| `misalignments` | **原文で確認**: (8.10) の 3×3 の式 (TeX 版・平文版とも)。「The inverse is applied at the exit」「the geometry of the design orbit is determined by the saved values of CHI1, CHI2, and DZ」。`tsolrot` の式が (8.10) の $R_2R_1$ と一致し、入口・出口が逆変換であることは確認した (§5.1)。MULT が $-\mathrm{CHI}1$、$-\mathrm{CHI}2$ を渡す理由は、ANGLE = 0 と ANGLE $\neq 0$ (`tbrot`) の向きをそろえるためと確認 (§5.1)。ヘルプの向きとの対応 (MULT の実際の向きは (8.10) の文字どおりの読みの逆) は未確認。 |
| `multipole_with_nonzero_ANGLE` の $g_{kn}$ | **原文で確認、確定した食い違い**: $(2k-1)!!$ は $(2k+1)!!$ の誤り (§6.3)。 |
| 同じ節の $(x + iy)^{n+k}$ | **原文で確認、確定した食い違い**: $z$ の冪は $n + k + 1$。 |
| 同じ節の「$n + k \le 21$」 | 原文で確認 (HTML では不等号が欠けて見える)。K64 の `nmult` = 21 と整合。 |
| `K0`・`K1` の説明 | 「Positive sign means horizontal focusing」は K0 にも付く (K0 は偏向であって収束ではない。コピーの誤記とみられる。確認していない)。SK$n$ の説明の「$L$ is the length」(K$n$ は effective length)。 |
| ヘルプに無いもの | (1) ANGLE = 0 の分割の構成 (8.7) と刻み数の式 (8.8)。(2) $A_s$ の導出と $g_{kn}$ の漸化式 (§6.2)。(3) 端面の角の式 (8.9)、FB1・FB2 の式。(4) ずれの `tsolrot` の式の詳細。 |

**K64 の誤り (確定。v2 は修正)**: §6.4 の $L'$ の二重補正 (gs と ODE で確認。B331)、§6.5 の K1 ≠ 0 の端の重複 (gs で確認。B340)。**取り下げ (誤りではない)**: MULT の CHI1・CHI2 の符号 (gs で確認。(8.10) の $R$ を要素の向きとして読んだものと同じ。§5)。**K64 の実装の癖・疑い**: (8.9) の $n \ge 1$ の項と pysad v2 の厳密な幾何の水平方向の違い (§4.4。どちらが厳密かは未確認)。K1 の暗黙の定義 (B160) は MULT の K1 にも同じ。pysad の k64 は `tmulk` (K64)、v2 は面の厳密な幾何 (B223・B224・B292)。

## 付録 B. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `akn`・`akr0`・`akr` | 多極の複素強度 $K_n + i\,SK_n$ (入力) / 基準面に回した後 ($\cdot\,\texttt{cr1}^{n+1}$) | `tmulti` の引数 / `tffs.f:1411-1418` (`p_K0R_MULT`) | `pm['akr']` |
| `cr1`・`theta1` | $e^{-i\theta_1}$、$\theta_1 = \tfrac12\arg(\mathrm{K1} + i\,\mathrm{SK1})$ | `akang` (`tfloor.f:530`) | `lattice._mult_params` |
| `theta2` | $\mathrm{ROTATE} + \mathrm{DROTATE} + \theta_1$ (入口で座標を回す角) | `tffs.f:1398-1402` | `pm['theta2']` |
| `aninv(k)`・`fact(k)` | $1/k$ ($k \ge 1$、`aninv(0)=aninv(1)=1`)、$k!$ | `tffs.f:1155-1186` (`multa` モジュール) | `_FACT`・`k!` |
| `ndiv`・`ndiv0`・`wi` | 刻みの数。`ndiv0` = `cmp%ivalue(1,p_NM_MULT)`。$w_i = 1/\mathrm{ndiv}$ | `tffs.f:1423-1434`、`tmulti.f:157` | `pm['ndiv']` |
| `nmmax` | K・SK が 0 でない最大の $n$ | `tffs.f:1405-1413` | `pm['nmmax']` |
| `eps00m`・`ampmaxm`・`ndivmaxm` | $0.005$、$0.05$、$1000$ | `toplvl.f:292-293` | `_EPS00M`・`_AMPMAXM` (`track6.py:1919`) |
| `ak1`・`al1`・`ak01` | 本体 1 回の K1・長さ・K0 ($\times w_i$ または $w_i/2$) | `tmulti.f:157-163` | (内部) |
| `akrm` | $\mathrm{akr}(n)\,w_i$ (1 回のキックの強さ) | `tmulti.f:163` | (内部) |
| `dx0`・`dy0` | $\mathrm{K0}/\mathrm{K1}$、$\mathrm{SK0}/\mathrm{K1}$ (本体の中心のずれ、(8.6)) | `tsolqu.f:55-56` | `track6._tsolqu` |
| `dofr(n)` | 非線形フリンジを入れる成分 | `tffs.f:1436-1443` | `dofr` |
| `fe` | $-\tan\psi/L$ (端面の角の係数) | `tmulk` (`tmulti.f:11`) | `track6._tmulk` |
| `gkn(n,k)` | $g_{kn}$ | `gkninit` (`tffs.f:1188`) | `multa.GKN` |
| `rho0`・`rhob`・`phib` | $L/\mathrm{ANGLE}$、$L/\phi_b$、$\mathrm{ANGLE} + \mathrm{K0}$ | `tmulta.f` | `multa.k_multa` |
| `rk(k)`・`csl`・`csr` | $(\rho + x)^{1/2-k}/\sqrt\rho$、$C_l$、$C_s$ ((8.16)) | `tmulta.f:170-200` | `multa.k_multa` |
| `ak0r`・`ak1r`・`ak0n`・`ak1n` | ANGLE ≠ 0 で、BEND の本体に回す $\mathrm{Re}\,\mathrm{K0}$、$\mathrm{Re}\,\mathrm{K1}$ (`ak(0)`・`ak(1)` からは実部を除く)。`ak0n`・`ak1n` は 1 スライスあたり (÷ndiv) | `tmulta.f:52-53, 92-93` | `multa.k_multa` |
| `phib`・`phin`・`aln`・`als` | $\phi_b = \mathrm{ANGLE} + \mathrm{Re}\,\mathrm{K0}$、$\phi/\mathrm{ndiv}$、スライスの円弧長 (補正後。$L'/\mathrm{ndiv}$ からさらに $c_2/\mathrm{ndiv}$ を引く、§6.4)、最初の半スライス | `tmulta.f:74-101` | `multa.k_multa` |
| `tmulk` | 端面の角 (AE1・AE2) の薄いキック (8.9)。`fe` $= -\tan\psi/L$、`cx` の Horner 形 | `tmulti.f:4-31` | `track6._tmulk` |
| `alg`・`phig` | 設計軌道に沿った長さの半分 ($L/2$)、曲げ角の半分 ($\phi/2$)。回転の支点 (要素の中心) までの距離 (§5.1)。出口では `alg-al` が渡る | `tfgeo.f:212-216`、`tmulti.f:93, 305` | `pm` の内部 |
| `s0s`・`cchi1`・`schi1`・`dcchi1` | $-\mathtt{alg}$ (入口) / $\mathtt{alg}$ (出口、引数 `alg-al` の符号反転を含む)、$\cos\mathrm{CHI}1$、$\sin\mathrm{CHI}1$、$\cos\mathrm{CHI}1 - 1$。`dcchi2`・`cchi2`・`schi2` も同様 | `tsol.f:30-165` (`tsolrot`)、`xsincos` (`tfloor.f:260`) | — (pysad の追跡は CHI1・CHI2・DZ 未実装) |
| `ds1`・`ds2` | 回転後の第 3 成分の面からの変位 (§5.1)。`a = ds2/pz2` が面まで戻す経路長 | `tsol.f` | — |
| `tsolrot`・`tsolqu`・`tdrift` | ずれ・回転 (入口 `ent = .true.`、出口 `.false.`)、四極 + 双極の厳密な本体、K1 = 0 のときのドリフト (K0・SK0 のキック付き) | `tsol.f`、`tsolqu.f`、`tmulti.f:93, 305` | (内部)。入口の回転・ずらしは `lattice._mult_rot_in` ほか |
| `ampmax`・`eps00`・`ndivmax` (`tmulta`) | (8.8) の $a_{\max} = 0.05$ m、$\varepsilon_{00} = 0.005$、上限 2000 (`parameter` 文で `tmulta` の中に定数として持つ。ANGLE = 0 は `tffs.f` の `ampmaxm`・`eps00m`・`ndivmaxm` (上限 1000))。$\varepsilon = \varepsilon_{00}\cdot\mathrm{EPS}$ (EPS = 0 は 1)。$\varepsilon$ の意味は §4.3 | `tmulta.f:15-16, 260` | `multa.AMPMAX`・`EPS00` |
| `rtaper` | テーパー比。放射と COD の補正 (`RAD`・`RADCOD`・`RADTAPER`) が有効なとき $1 + (\delta_{\rm COD}\text{ の要素の両端の平均}) + \mathtt{dptaper}$ (磁場を局所の運動量に比例させる)、そうでなければ 1。全成分に $\mathtt{akr}(n) \to \mathtt{akr}(n)\cdot\mathtt{rtaper}$、`tmulta` では $\mathrm{ak}(0) \mathrel{+}= (\mathtt{rtaper} - 1)\,\phi$ (ANGLE の分の場もテーパーする) | `tturn.f:389-392`、`tmulti.f:95, 102`、`tmulta.f:50-51` | (内部) |
| `psi1`・`psi2` | 入口・出口の端面の角 $\mathrm{E1}\cdot\phi + \mathrm{AE1}$、$\mathrm{E2}\cdot\phi + \mathrm{AE2}$ ($\psi$ は数学の標準記号と同じ向きの角。逆向きに置いた要素は入れ替わる)。`tmulk` の `apsi` に渡る | `tffs.f:1339-1351`、`tmulti.f:166, 250` | `pm['psi1']` ほか |
| `phia`・`phi` | ANGLE (曲げ角 [rad])。`tmulti` の引数が `phia`、`tmulta` の中では `phi` | `tmulti.f:78`、`tmulta.f` | `multa.k_multa` |
| `dtheta`・`theta`・`chi2` | DROTATE、ROTATE、CHI2 (ANGLE $\neq 0$ では `tbshift`・`tbrot` に渡す。CHI1・DZ は使わない) | `tmulta.f:40-43` | (内部) |
| `eps0` | 台本の EPS (0 は 1 として扱う)。(8.8) の $\varepsilon = \varepsilon_{00}\cdot\mathrm{EPS}$ に掛かる | `tmulta.f:43-47` | `eps0` |
| `fb1`・`fb2` | K0 の端 (FB1・FB2) の幅 (偏向磁石の F1 に当たる)。`tblfri` と $L'$ の補正に使う | `tffs.f:1344-1367` | `pm['fb1']` ほか |
| `mfring` | FRINGE の値 (1 入口、2 出口、3 両方 の判定に使う。§4.4) | `tmulti.f:126, 135, 253, 274` | `mfring` |
| `trans` | CALC の `tmulte`・`tmultae` の引数の転送行列 (6×12。§7)。数学の標準記号ではなく、Fortran の配列名 | `tmulte.f:66, 92` | `lattice` の行列 |
| `tmultr5` | K64 のソース (`tmulti.f`・`tmulte.f`・`tmulta.f`) に見つからず、この章では使わない (名前だけの指摘は未確認) | — | — |
| `fact(k)`・`aninv(k)` | $k!$、$1/k$ (`aninv(0)=aninv(1)=1`)。(8.8)・Horner 形で使う | `tffs.f:1155-1186` | `_FACT` |
