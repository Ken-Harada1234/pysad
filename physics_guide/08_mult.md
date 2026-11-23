# 多極 MULT: ポテンシャルの展開、本体の分割、ずれ、ANGLE ≠ 0

**この章で分かること**: MULT は、1 つの要素に K0〜K21 (と SK0〜SK21) の多極を重ねられる要素である。この章は、(1) 多極のポテンシャル $A_s$ の展開と K・SK の意味 (ANGLE = 0 は直線の多極、ANGLE ≠ 0 は曲率つきの展開)、(2) ANGLE = 0 の追跡の構成 (厳密な四極 + 双極の本体 `tsolqu` と、$n \ge 2$ の薄いキック `tmulk` の対称な分割)、刻みの数 `ndiv`、端のフリンジ、(3) ずれ・回転 (DX、DY、DZ、CHI1、CHI2、ROTATE)、(4) ANGLE ≠ 0 の $A_s$ を Maxwell 方程式から導き、ヘルプの $g_{kn}$ との食い違いを確定する、までを、式を省かずに追う。四極の本体は [04_quad.md](04_quad.md)、薄い多極のキックと非線形フリンジは [07_thin_multipole.md](07_thin_multipole.md)、偏向磁石の本体と端は [05_bend.md](05_bend.md)、ソレノイド中の本体は [09_solenoid.md](09_solenoid.md)。RF (VOLT) つきの MULT は [10_cavity.md](10_cavity.md)。

**読む順**: §1 (結論) → §3 (ポテンシャル) → §4 (ANGLE = 0 の構成) → §6 (ANGLE ≠ 0) → §8 (確認)。

## 1. 結論

- **K_n・SK_n の意味**: $\mathrm{K}n$ ($n = 0, 1, 2, \ldots$) は積分強度 $B^{(n)}L/B\rho$ [1/m$^n$] で、横方向のポテンシャル $U_n = \mathrm{Re}\,[(\mathrm{K}n + i\,\mathrm{SK}n)\,z^{n+1}]/(n+1)!$ ($z = x + iy$、長さあたりにするには $L$ で割る) を作る。キックは $\Delta(p_x, p_y) = (-\mathrm{Re}, +\mathrm{Im})\,[(\mathrm{K}n + i\,\mathrm{SK}n)\,z^n/n!]/(1+\delta)$ (追跡の内部変数 $q$ で)。SK$n$ > 0 は、K$n$ の磁石を ROTATE = $90°/(n+1)$ 回したものと**同じ** (ヘルプの記述を数値で確認: 差 $10^{-17}$)。
- **基準面の取り方**: K1 + i SK1 の位相の半分 $\theta_1 = \tfrac12\arg(\mathrm{K1} + i\,\mathrm{SK1})$ だけ座標を回して、スキュー四極の成分を 0 にする (ヘルプの「reference plane is defined so that the skew quadrupole component becomes zero」)。各成分は $(\mathrm{K}n + i\,\mathrm{SK}n)\,e^{-i(n+1)\theta_1}$ に変わる。
- **ANGLE = 0 の追跡** (`tmulti`): 本体は $B(h/2)\,K\,B(h)\,K\cdots K\,B(h/2)$ ($h = L/\mathrm{ndiv}$)。$B$ は K1 と K0 (と SK0) だけの**厳密な線形の本体** (`tsolqu`)、$K$ は $n \ge 2$ の多極の**薄いキック**。対称な 2 次の分割で、誤差は $\mathrm{ndiv}^{-2}$ に比例する (数値確認)。刻み数は $\mathrm{ndiv} = 1 + \bigl\lfloor\max_{n\ge2}\sqrt{a_{\max}^{\,n-1}\,|K_n L|/(6\,(n-1)!\,\varepsilon)}\,\bigr\rfloor$、$a_{\max} = 0.05$ m、$\varepsilon = 0.005\,\mathrm{EPS}$ (EPS = 0 は 1)、上限 1000。
- **ずれ**: DX、DY、DZ、CHI1、CHI2、ROTATE (= CHI3) は剛体の座標変換 (ヘルプの 3×3 の式)。この章では DX・DY・ROTATE の向きだけ数値で確認した。CHI1・CHI2・DZ の幾何は読みだけ。
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

になる。K64 の `tmulk`・`tmulti` の `px(i)=px(i)-dble(cx)`、`py(i)=py(i)+imag(cx)` ($\texttt{cx}$ = $C(z)$ の Horner 形) がこれ。$\texttt{cx}=\bigl(\cdots((a_N z/N + a_{N-1})z/(N-1) + \cdots\bigr)$ は、`cx=(cx+akn(kord))*cx0*aninv(kord)` ($\texttt{aninv}(k) = 1/k$、`tffs.f:1162`) を $k = N, \ldots, 1$ と回したもので、$a_n$ が通る係数は $\prod 1/k\cdot z^{n}$ なので $C(z)$ と等しい (§8 の数値確認で、薄い多極と K64 の追跡が $10^{-19}$ で一致)。

### 3.2 SK_n は K_n の回転 (ヘルプの記述の確認)

$K = 0$ で $\mathrm{SK}n > 0$ のとき (8.1) は $\mathrm{Re}\,[i\,\mathrm{SK}n\,z^{n+1}]/(n+1)!$。$i = e^{i\pi/2}$ なので、座標を $z \to e^{i\phi}z$ と回すと $e^{i(n+1)\phi}$ がかかる。$(n+1)\phi = \pi/2$ すなわち **$\phi = 90°/(n+1)$** の回転が、SK$n$ を同じ大きさの K$n$ にする。ヘルプの「SK$n$ > 0 は、K$n$ の磁石を ROTATE = $90°/(n+1)$ にしたもの」と一致する (K1 なら 45°、K2 なら 30°、K3 なら 22.5°)。ROTATE の向きは、`tsolrot` の入口で $(x, y) \to (\cos\theta\,x - \sin\theta\,y,\ \sin\theta\,x + \cos\theta\,y)$ (座標を $+\theta$ 回す) で、出口で逆に戻す ([mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) (c): この向きの差 0、逆向きの差 $6 \times 10^{-4}$)。この回転で (8.1) の係数は

$$
(\mathrm{K}n + i\,\mathrm{SK}n)\;\to\;(\mathrm{K}n + i\,\mathrm{SK}n)\,e^{-i(n+1)\theta}
\tag{8.4}
$$

になる (新しい座標 $z' = e^{i\theta}z$ で書くと $z^{n+1} = e^{-i(n+1)\theta}z'^{\,n+1}$)。

### 3.3 基準面 (θ1)

MULT は、K1 の回転の位相を**要素の基準面**として取り込む。$a = \mathrm{K1} + i\,\mathrm{SK1}$ の偏角の半分を $\theta_1$ として座標を $\theta_1$ 回し、(8.4) で K1 + i SK1 を実数 $|a|$ にする (`akang` (`tfloor.f:530`): $\theta_1 = \tfrac12\,\mathrm{atan2}(\mathrm{SK1}, \mathrm{K1})$、$\texttt{cr1} = e^{-i\theta_1}$。K1 < 0 の実数は $\theta_1 = \pi/2$、`cr1` $= -i$)。`tffs.f:1398-1420` は、実際に回す角を `theta2` $= \mathrm{ROTATE} + \mathrm{DROTATE} + \theta_1$、各成分を `akr(n)` $= (\mathrm{K}n + i\,\mathrm{SK}n)\cdot\texttt{cr1}^{\,n+1}$ ($\texttt{cr}$ を $n$ のたびに `cr1` 倍して作る) と決める。これで、K1 を持つ MULT の本体 (§4) は四極の主軸の座標で解ける。

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

(確認: $-\tfrac{\mathrm{K1}}{2L}\cdot 2d_yy = -\tfrac{\mathrm{SK0}}{L}y$、$\tfrac{\mathrm{K1}}{2L}\cdot 2d_xx = \tfrac{\mathrm{K0}}{L}x$。) 定数項は軌道に効かない (運動方程式に出ない) ので、本体の運動は、$(x, y)$ を $(x + d_x, y + d_y)$ に替えた四極の運動と同じ。`tsolqu` は `dx0=ak0x/ak`、`dy0=ak0y/ak`、`xi=x(i)+dx0` で本体の解 (04 章 (4.9)(4.10)) を $\texttt{xi}$、$\texttt{yi}$ に入れる (`tsolqu.f:55-56, 83-84`)。定数項が効くのは $z$ だけ ($z' = \partial H/\partial p$ は $\texttt{xi}$ を含む式 `w1*(xi**2*xs1-yi**2*xsh2)` から)。K1 = 0 のときは `tdrift` に K0・SK0 の一様なキックを渡す (`tsolqu.f:34-40`)。

### 4.2 分割と順序 (`tmulti`)

$H = H_{\rm body}' + U_{\rm kick}$ ($H_{\rm body}'$ は運動学 + $H_{\rm body}$) とみて、2 次の対称な分割 ([03 章](03_symplectic_integrators.md) §4) をする。長さ $L$ を $\mathrm{ndiv}$ 個に分け、$h = L/\mathrm{ndiv}$:

$$
M = B\!\left(\tfrac h2\right)\,\mathcal K\,B(h)\,\mathcal K\,B(h)\cdots\mathcal K\,B\!\left(\tfrac h2\right)\qquad(\mathcal K\text{ は }\mathrm{ndiv}\text{ 個、}B(h)\text{ は }\mathrm{ndiv}-1\text{ 個})
\tag{8.7}
$$

$B(a)$ は長さ $a$ の分の $H_{\rm body}'$ の流れ (`tsolqu(al1, ak1, bzs, ak01x, ak01y, …)`、K1・K0 の強さは $a/L$ 倍)、$\mathcal K$ は強さが $1/\mathrm{ndiv}$ の $n \ge 2$ の薄いキック (8.3)。コードでは、`do m=1,ndiv` の中で本体を呼んだあとに `al1=al*wm` を更新し (`wm` は $m = \mathrm{ndiv}$ だけ $w_i/2$、それ以外 $w_i = 1/\mathrm{ndiv}$。`tmulti.f:205`)、キック (`akrm(n)=akr(n)*wi`) を入れる。最初の本体は `al1=al*wi*.5` (`tmulti.f:157-163`)、最後の半分はループの外 (`tmulti.f:244-248`)。キックの式は Horner 形 $\sum_{n\ge2}\mathrm{akrm}(n)\,z^n/n!\,/(1+\delta)$ で (`tmulti.f:209-221`)、07 章のキックと同じ。**3 次以上の力勾配型の補正 (07 章 §6) は MULT には無い**。K1 ≠ 0 のとき $B$ は $\mathcal K$ と**可換でない**ので、誤差は 2 次 (Strang) になる。

放射 ($\texttt{krad}$) があると、本体ごとに `tradk` (光子の放出) を呼び、`ndiv` を `ndivrad` との大きい方にする (`tmulti.f:119`)。SPAC (空間電荷) があると `ndiv` を $|L|/(0.05\,\mathrm{EPS})$ 以上にする (`:110-113`)。この章では扱わない。

### 4.3 刻みの数 `ndiv`

`tffs.f:1423-1434`: $\varepsilon = \varepsilon_{00m}\cdot\mathrm{EPS}$ (EPS = 0 は 1)、$\varepsilon_{00m} = 0.005$ (`toplvl.f:292`)、$a_{\max} = $`ampmaxm`$ = 0.05$。

$$
\mathrm{ndiv} = \min\Bigl(1000,\ \max\bigl(1,\ \max_{2\le n\le n_{\max}}\bigl\lfloor\sqrt{a_{\max}^{\,n-1}\,\|\mathrm{K}n + i\,\mathrm{SK}n\|\,L\,/(6\,(n-1)!\,\varepsilon)}\bigr\rfloor + 1\bigr)\Bigr)
\tag{8.8}
$$

(`fact(n-1)` は階乗の表、`hypot` は絶対値。) $n = 0, 1$ は厳密な本体で処理するので入らない。$a_{\max}^{n-1}K_n/(n-1)!$ は、振幅 $a_{\max} = 50$ mm での $n$ 極の力の勾配 ($K_n a^{n-1}/(n-1)!$) で、それに $L$ をかけて $\varepsilon$ と比べる形。分母の 6 の根拠は**ソースに書かれていない** (Strang 分割の 3 次の誤差係数の見積りと推測。未確認)。pysad の `ndiv` は (8.8) と一致 (EPS = 4, 2, 1, 0.5, 0.25 で 5, 6, 9, 12, 17。K2 = 40 /m²、$L$ = 1 m)。

**誤差の次数**: ANGLE = 0、K2 = 40 /m² (積分強度)、$L$ = 1 m、初期値 $(4, 0.8, -3, -0.4)\times$ (mm, mrad, mm, mrad) で、厳密な運動方程式 (DOP853) との差は、EPS が 4 → 0.25 のとき $6.4 \times 10^{-7} \to 5.5 \times 10^{-8}$。誤差の比は `ndiv` の比の 2 乗とほぼ等しい (1.44 対 1.44、2.25 対 2.25、1.78 対 1.78、2.01 対 2.01)。**2 次**である ([mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) (a))。K1 = 0.8 を加えた (b) でも、`ndiv` を 6 → 12 で $1.7 \times 10^{-6} \to 4.2 \times 10^{-7}$ (2 次)。

### 4.4 フリンジ

- **非線形フリンジ**: `fringe` (DISFRIN = 0) のとき、`mfring ≠ 2` で入口、`mfring ≠ 1` で出口に、$n$ ごとに `ttfrins(…, n*2+2, ±akr(n), al, bzs)` を呼ぶ (`tmulti.f:124-132, 273-281`。`n*2+2` は極の数 $2(n+1)$)。対象は `dofr(n)` が真の成分 (`tffs.f:1436-1443`): K$n$ ≠ 0 かつ ($n \le 2$ か、$\|K_n\|/L\cdot a_{\max}^{n+1}\,\texttt{aninv}(n+1) > \varepsilon$。`aninv(n+1)` は $1/(n+1)$ で、$1/(n+1)!$ ではない)。式は [07 章](07_thin_multipole.md) §9 (Maxwell から)。
- **線形フリンジ (F1・F2)**: 四極 (K1) の線形フリンジ ([04 章](04_quad.md) §7.4 と同じ式)。`f1in`・`f2in`・`f1out`・`f2out` は `tsetfringep` が $\mathrm{K1}/L$ から決める (`tffs.f:1380-1385`)。出口は $a \to -a$。`mfring = 1, 3` が入口、`2, 3` が出口。
- **K0 の端 (FB1・FB2)**: `tblfri(…, al, ±akr(0), fb)` (入口は `akr(0)`・`fb1`、出口は `-akr(0)`・`fb2`)。偏向磁石の線形フリンジ V ([05 章](05_bend.md) 式 (5.42)) を、$\rho_b = L/\mathrm{K0}$ ($x$ 方向)・$L/\mathrm{SK0}$ ($y$ 方向) で入れる。
- **端面の角 (AE1・AE2)**: `tmulk(…, al, akr, psi, nmmax)` (`tmulti.f:4-31, 166-168, 250-252`。`fe` は 16 行目):

$$
C_\psi(z) = -\frac{\tan\psi}{L}\;z\sum_{k=0}^{n_{\max}}\frac{(\mathrm{K}k + i\,\mathrm{SK}k)\,z^k}{k!},\qquad
\Delta p_x = -\mathrm{Re}\,C_\psi/(1+\delta),\quad \Delta p_y = +\mathrm{Im}\,C_\psi/(1+\delta)
\tag{8.9}
$$

($\texttt{fe} = -\tan\psi/L$、`cx=fe*(cx+akn(0))*cx0`。) K0 だけの場合、$\Delta p_x = +(\tan\psi\,\mathrm{K0}/L)\,x$、$\Delta p_y = -(\tan\psi\,\mathrm{K0}/L)\,y$ で、偏向磁石の端面の標準の収束 (水平 $\tan\psi/\rho$、垂直 $-\tan\psi/\rho$、[05 章](05_bend.md) §5) と同じ向き・大きさ。**高次の項を $z\,C(z)$ とする根拠はソースに無い**。端面を斜めに切ったときの余分な長さ $x\tan\psi$ に、場 $C(z)$ をかけたものに、垂直方向の Maxwell の収束 ($y$ と $\mathrm{Re}\,C$ の積) が加わった形と読める (推測。導出していない)。pysad の v2 は、この端を面の厳密な幾何 (`_mult_edge_exact`) で置き換える (`track6.py:2040-2088`。docs/30 の B223・B224・B292)。k64 は `tmulk` のまま。

## 5. ずれ・回転 (DX、DY、DZ、CHI1、CHI2、ROTATE)

ヘルプ (`misalignments`) は、入口で粒子の座標を剛体の変換

$$
\begin{pmatrix}x\\y\\\Delta s\end{pmatrix}_1 = R_3(\chi_3)\,R_2(\chi_2)\,R_1(\chi_1)\begin{pmatrix}x - \mathrm{DX}\\ y - \mathrm{DY}\\ -\mathrm{DZ}\end{pmatrix}_0,\qquad
R_1 = \begin{pmatrix}c_1 & 0 & -s_1\\0&1&0\\s_1&0&c_1\end{pmatrix},\ R_2 = \begin{pmatrix}1&0&0\\0&c_2&-s_2\\0&s_2&c_2\end{pmatrix}
\tag{8.10}
$$

($R_3$ は $z$ 軸まわりの回転。$c_i = \cos\mathrm{CHI}i$、$s_i = \sin\mathrm{CHI}i$。出口で逆変換) とする。粒子の線 (位置と向き) をこの変換で動かし、$\Delta s = 0$ の面まで直進させたのが入口の写像で、`tsolrot` の $B_z = 0$ の枝 (`tsol.f` の `ent = .true.`): `x-dx`、`y-dy` → $\mathrm{CHI}1$ の回転 (`cchi1*x - schi1*(s0s - dz)`、`px` も同様) → $\mathrm{CHI}2$ の回転 → 面までの直進 (`a=ds2/pz2`、`x-a*px`、`y-a*py`、`z+a`) → ROTATE (`theta2`) の回転 (`tsol.f` の `chi3`)、の順。MULT の `tmulti` は、`-chi1`、`-chi2` を渡して呼ぶ (符号は反転)。設計軌道の幾何は、保存された CHI1、CHI2、DZ で決まり、DX、DY、ROTATE は現在の値を使う (ヘルプ)。

**確認したこと**: DX・DY は、$(x - \mathrm{DX}, y - \mathrm{DY})$ で通して $\mathrm{DX}$、$\mathrm{DY}$ を足すのと同じ (差 0)。ROTATE は、入口で座標を $+\theta$ 回し出口で戻すのと同じ (差 0、逆向きは $6 \times 10^{-4}$)。**確認していないこと**: CHI1・CHI2・DZ の幾何 (ヘルプの (8.10) の $\Delta s$ の符号とコードの `s0s`・`dz` の符号の対応、`alg` (設計軌道に沿った長さ) を回転の支点にとる意味)。`tsol.f` の式を読んだだけ。

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

曲線座標 $(x, y, s)$ ($h = 1/\rho$、スケール因子 $1, 1, 1 + hx$) で、$s$ によらない場が $\vec A = a_s\hat s$ だけで書けるとする。$H$ の中の項は $-(1 + hx)a_s \equiv A_s$ ([01 章](01_variables_and_conventions.md) 式 (1.1))。$F \equiv (1 + hx)\,a_s$ とおくと $\vec B = \nabla\times\vec A$ の成分は

$$
B_x = \partial_ya_s = \frac{\partial_yF}{1 + hx},\qquad B_y = -\frac{1}{1 + hx}\,\partial_x\bigl((1 + hx)a_s\bigr) = -\frac{\partial_xF}{1 + hx}
\tag{8.12}
$$

($B_s = 0$。$\nabla\cdot\vec B = (1 + hx)^{-1}\bigl[\partial_x((1 + hx)B_x) + \partial_y((1 + hx)B_y)\bigr]$ で、$(1 + hx)B_x = \partial_yF$、$(1 + hx)B_y = -\partial_xF$ だから $(1 + hx)^{-1}(\partial_x\partial_yF - \partial_y\partial_xF) = 0$。恒等的に満たされる。) **場が真空の Maxwell 方程式を満たす条件**は $(\nabla\times\vec B)_s = \partial_xB_y - \partial_yB_x = 0$ で、(8.12) を入れて

$$
\partial_x\Bigl(\frac{\partial_xF}{1 + hx}\Bigr) + \frac{\partial_y^2F}{1 + hx} = 0\ \ \Longleftrightarrow\ \ \nabla^2F = \frac{\partial_xF}{\rho + x}
\tag{8.13}
$$

(第 2 式は第 1 式に $1 + hx$ をかけ、$\partial_x(1/(1 + hx)) = -h/(1 + hx)^2$ を使った。$\nabla^2 = \partial_x^2 + \partial_y^2$。) $\rho \to \infty$ では $\nabla^2F = 0$ (調和関数 $\mathrm{Re}\,z^{n+1}$)。

**$\rho$ の冪の展開**。調和関数 $z^m$ ($m = n + 1$) を出発点に、$F = \sum_{k\ge0}g_k\,u^{a_k}z^{m+k}$、$u = \rho + x$、$a_k = \tfrac12 - k$ と置く。1 項ずつ計算する: $z^m$ は調和 ($\nabla^2z^m = 0$)、$\partial_x u^a = au^{a-1}$、$\partial_yu = 0$ だから

$$
\nabla^2(u^az^m) = a(a-1)\,u^{a-2}z^m + 2am\,u^{a-1}z^{m-1},\qquad
\frac{\partial_x(u^az^m)}{u} = a\,u^{a-2}z^m + m\,u^{a-1}z^{m-1}
$$

(第 1 式の第 2 項は $2\,\partial_xu^a\,\partial_xz^m = 2au^{a-1}mz^{m-1}$。) (8.13) の左辺 − 右辺は、$u^{a-2}z^m$ の係数が $a(a-2)$、$u^{a-1}z^{m-1}$ の係数が $m(2a - 1)$。$k$ 番目の項の第 1 の型 ($u^{a_k-2}z^{m+k}$) と、$(k+1)$ 番目の項の第 2 の型 ($u^{a_{k+1}-1}z^{m+k}$。$a_{k+1} - 1 = a_k - 2$) が同じ形になり、その和が 0:

$$
g_k\,a_k(a_k - 2) + g_{k+1}\,(m + k + 1)(2a_{k+1} - 1) = 0
$$

$a_k(a_k - 2) = (\tfrac12 - k)(-\tfrac32 - k) = \tfrac14(2k - 1)(2k + 3)$、$2a_{k+1} - 1 = -2(k + 1)$ を入れ、$k + 1 \to k$ と読み替えると

$$
\boxed{\;g_{kn} = g_{k-1,n}\,\frac{(2k - 3)(2k + 1)}{8\,k\,(k + n + 1)},\qquad g_{0n} = 1\;}
\tag{8.14}
$$

($m = n + 1$ を使った。$g_{0n} = 1$ は、$\rho \to \infty$ で直線の多極 (8.1) に戻る規格化。) 閉じた形は $g_{kn} = -(2k+1)!!\,(2k-3)!!\,(n+1)!/(8^k\,(n+k+1)!\,k!)$ ($(-1)!! = 1$、$(-3)!! = -1$)。(8.14) は K64 の `gkninit` (`tffs.f:1188-1200`) と**同じ**。結局

$$
A_s = \sum_{k,n}g_{kn}\,\mathrm{Re}\,\frac{(\mathrm{K}n + i\,\mathrm{SK}n)\,z^{n+k+1}}{(n+1)!}\,\frac{(\rho + x)^{1/2-k}}{\sqrt\rho}
\tag{8.15}
$$

($\rho$ で割る因子は $\rho \to \infty$ で $k = 0$ の項が 1 になるため。$\sqrt{1 + x/\rho}$ を $k = 0$ の項が持つのは、$(1 + hx)\,a_s$ の中の因子から来る。) $\rho$ について展開するので、$k$ の項は $(1/\rho)^k$ で小さくなり、K64 は $n + k \le 21$ で打ち切る (`nmult` = 21)。

**キック (`tmulta`)**。$\Delta p_x = -\partial_xA_s$、$\Delta p_y = -\partial_yA_s$。$\partial_x[u^{a_k}z^{n+k+1}] = u^{a_k}(n + k + 1)z^{n+k} + a_k\,u^{a_k-1}z^{n+k+1}$、$\partial_y[u^{a_k}z^{m}] = i\,m\,u^{a_k}z^{m-1}$ から

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
4. 各スライス: `tbend` (半分、入口の端面と Maxwell フリンジ込み) → [キック (8.16) → `tbend` (全体)] を $\mathrm{ndiv} - 1$ 回 → キック → `tbend` (半分、出口の端面)。BEND の $L'$ の補正 ([05 章](05_bend.md) 式 (5.43)) は、`aln` から $(\phi\,\mathrm{fb})^2/(48L)\,\sin(\ldots)/\sin(\phi/2)/\mathrm{ndiv}$ ずつ引く形 (`tmulta.f:82-90`)。**同じ補正が `tffs.f:1368` でも `p_L_MULT` に入っている** (`tmulta` に渡る `al` はすでに補正済み) ので、二重に補正している可能性がある (疑い。数値では確認していない。設計軌道の $x, p_x$ は 0 のままで、$z$ は比べていない)。
5. 逆のずれ。

CALC (`tmultae`) は `tmulta` と同じ分割で行列を作る (`multae.py`)。キックを $(1+\delta)$ で割らないなど、追跡と違う点がある (読みだけ)。

## 7. 行列 (CALC・Emittance[])

`tmulte` (ANGLE = 0) は、追跡と同じ刻み数の式 (8.8) (`ampmax`・`eps00m`) と対称な分割 (8.7) を、COD のまわりで微分して行列にする (`tmulte.f:131-145`)。`tmulte` は、ずれ・回転に `tsolrote`、K0・K1 に `tsolque` を使う ([04 章](04_quad.md) §8、[09 章](09_solenoid.md))。この章では、CALC の MULT の行列と追跡の写像の比較はしていない。

## 8. 確認

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `mult_01_as_expansion.py` | (a) ヘルプと K64 の $g_{kn}$ の比較、(b) (8.13) の残差 (sympy。$k \le 3$ の打ち切り、冪 $n + k + 1$ と $n + k$)、(c) K64 の漸化式の閉じた形 | [mult_01_as_expansion.out.txt](outputs/mult_01_as_expansion.out.txt) |
| `mult_02_tracking.py` | (a) ANGLE = 0、K2: `ndiv` の式と誤差の次数 (2 次)。(b) K1 + K2。(c) SK$n$ = K$n$・ROTATE $90°/(n+1)$、DX・DY、ROTATE の向き。(d) ANGLE ≠ 0、K2: (8.15) から作った曲線座標の ODE との比較 | [mult_02_tracking.out.txt](outputs/mult_02_tracking.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/mult_01_as_expansion.py` (数十秒、sympy)、`py -3.12 physics_guide/scripts/mult_02_tracking.py` (数十秒)。

主な数字 (電子 3 GeV、$\delta = 0$): ANGLE = 0 で、EPS = 4 → 0.25 の `ndiv` 5 → 17 に対し、厳密解との差が $6.4 \times 10^{-7}$ → $5.5 \times 10^{-8}$ (2 次)。SK1・SK2・SK3 と、K$n$・ROTATE = 45°・30°・22.5° の差は $10^{-17}$ 以下。ANGLE = 0.2 ($\rho = 5$ m)、K2 = 40 /m²: 扇形の偏向だけの追跡と ODE の差が $4.8 \times 10^{-7}$ (硬い端の Maxwell フリンジ。ODE に入れていない)。K2 の効果 ($5.7 \times 10^{-4}$) を、K2 = 0 との差で比べると、追跡と ODE の差は EPS = 4 → 0.25 で $6.9 \times 10^{-7}$ → $2.3 \times 10^{-7}$ (相対 $4 \times 10^{-4}$)。この残りは、EPS を小さくしても $2 \times 10^{-7}$ に残る。K2 のキックで軌道が変わったとき、Maxwell フリンジの大きさが変わる交差項 (同じ桁) と推測するが、確認していない。

## 9. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、$\delta = 0$。ANGLE = 0 は $L$ = 1 m (確認 (a)(b))、0.3 m (確認 (c))、K2 = 40 /m²、初期値は mm・mrad の桁。ANGLE ≠ 0 は ANGLE = 0.2 rad、$L$ = 1 m、扇形 (E1 = E2 = 0)、K2 のみ (K0 = K1 = 0)、DISFRIN = 1。放射・スピン・ソレノイド・RF・ずれ (CHI1・CHI2・DZ) は含まない。周長・磁場の上限は使っていない。

**この章の範囲外**: ソレノイド中の MULT (`bzs ≠ 0`)、RF つきの MULT (`tmultiacc`。10 章)、放射 (`tradk`、`ndivrad`)、空間電荷 (`spkick`)、MULT の `tmulte` の行列と追跡の比較 (§7)。

**未確認・限界**:
- **CHI1・CHI2・DZ** の符号と回転の支点 (§5)。DX・DY・ROTATE の向きだけ数値で確認した。
- **端面の角の式 (8.9)**: 高次の項を $z\,C(z)$ とする根拠は未導出 (§4.4)。K0 の極限だけ標準の端面収束と一致。pysad v2 は別の (幾何厳密な) 式で、K64 と違う。
- **刻み数の式 (8.8)** の係数 6 の根拠。2 次の次数は数値で確認した。
- **ANGLE ≠ 0**: (8.15) は Maxwell 方程式の厳密な解 ($k$ について無限級数) の $n + k \le 21$ での打ち切り。K0・K1 の ($k \ge 1$ の) 項は入れず、本体 (BEND の円弧) で扱う。追跡との比較 (§8) は、Maxwell フリンジに起因する $\sim 2 \times 10^{-7}$ の残りまでで、これ以上は分離していない。複数の多極 ($n \ge 3$) や、スキュー (SK$n$)・K1 ≠ 0 の ANGLE ≠ 0 は、ODE との比較をしていない。
- **$L'$ の補正の二重適用の疑い** (§6.4)。

## 付録 A. ヘルプの記述・K64 との違い (MULT)

原文 (`SADHelp.html` の `MULT`、`misalignments`、`multipole_with_nonzero_ANGLE`) で確認した。

| 箇所 | 内容 |
|---|---|
| MULT の説明 (reference plane) | **原文で確認**: 「the reference plane is defined so that the skew quadrupole component becomes zero」。§3.3 のとおり ($\theta_1$)。 |
| SK$n$ の説明 (「ROTATE = $90/(n+1)$ DEG」) | **原文で確認**、数値で一致 (§3.2)。ROTATE の値は $90/(n+1)$ 度 (K1 なら 45、K2 なら 30、…。K6 では 12.857…) と書かれている。「rotated around z-axis by $-90/(n+1)$ degree」は、座標を $+\theta$ 回すことに当たる。 |
| `misalignments` | **原文で確認**: (8.10) の 3×3 の式 (TeX 版・平文版とも)。「The inverse is applied at the exit」「the geometry of the design orbit is determined by the saved values of CHI1, CHI2, and DZ」。コードとの符号の対応は未確認 (§5)。 |
| `multipole_with_nonzero_ANGLE` の $g_{kn}$ | **原文で確認、確定した食い違い**: $(2k-1)!!$ は $(2k+1)!!$ の誤り (§6.3)。 |
| 同じ節の $(x + iy)^{n+k}$ | **原文で確認、確定した食い違い**: $z$ の冪は $n + k + 1$。 |
| 同じ節の「$n + k \le 21$」 | 原文で確認 (HTML では不等号が欠けて見える)。K64 の `nmult` = 21 と整合。 |
| `K0`・`K1` の説明 | 「Positive sign means horizontal focusing」は K0 にも付く (K0 は偏向であって収束ではない。コピーの誤記とみられる。確認していない)。SK$n$ の説明の「$L$ is the length」(K$n$ は effective length)。 |
| ヘルプに無いもの | (1) ANGLE = 0 の分割の構成 (8.7) と刻み数の式 (8.8)。(2) $A_s$ の導出と $g_{kn}$ の漸化式 (§6.2)。(3) 端面の角の式 (8.9)、FB1・FB2 の式。(4) ずれの `tsolrot` の式の詳細。 |

**K64 の実装の癖・疑い**: §6.4 の $L'$ の二重補正 (疑い)。K1 の暗黙の定義 (B160) は MULT の K1 にも同じ。pysad の k64 は `tmulk` (K64)、v2 は面の厳密な幾何 (B223・B224・B292)。

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
