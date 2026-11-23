# 16. アンジュレータ (UND) とビーム・ビーム (BEAMBEAM)

**この章で分かること**: SAD の `UND` と `BEAMBEAM` が追跡で行う写像を、ハミルトニアンと近似から組み立てる。どちらも**ヘルプに節がない** (ヘルプの要素の一覧に `UND`・`BEAMBEAM` の記述はない)。式はすべて K64 のソース (`undulator.f`、`bbstrhl.f`) から起こした。(1) UND: ベクトルポテンシャル $A_x, A_y$ の模型 (真空の Maxwell 方程式を満たす条件)、傍軸のハミルトニアン、スライスごとの陰的なシンプレクティック写像 (1 次)。(2) BEAMBEAM: 2 次元ガウス分布の電場 (Bassetti-Erskine)、スライスごとの 6 次元写像 (衝突点までの drift、キック、戻し)、強ビームの包絡とスライス、交差角・軌道・傾きの座標変換、ビームストラールング。2 次元ポテンシャルの数値積分と DOP853 を基準に確認する。

**範囲**: UND は追跡の写像 (K64 は磁場を設定できない: v2 だけ磁場つきの追跡と放射)。BEAMBEAM は弱-強モデルの追跡だけ (光学・Emittance[] には入らない)。**交差角の変換 (Lorentz boost に相当する 6 変数の閉じた写像) は式をソースから写し、シンプレクティック性と往復の恒等だけを確認した。導出はしていない**。光子リスト・ルミノシティ・制動放射 (BBBrem) は [docs/14_beambeam.md](../docs/14_beambeam.md) を参照。

**前提**: [13_radiation.md](13_radiation.md) (`anrad`・`cuc`・`cave`・`cuu`、光子の放射)、[03_symplectic_integrators.md](03_symplectic_integrators.md) (シンプレクティックな分割)、[02_drift.md](02_drift.md) (傍軸のドリフト)。

**読む順**: §1 (結論) → §3 (UND) → §6〜§8 (BEAMBEAM) → §12 (確認)。

## 1. 結論

- **UND の模型**: $A_x$、$A_y$ を (16.1) の形に与える。軸上で $B_y = F_0\cos(K_z s)$、$B_x = G_0\cos(K_z s + \varphi)$ ($F_0 = \mathtt{FBX}\,c/p$ [1/m]、$G_0 = \mathtt{FBY}\,c/p$、$K_z = 2\pi\,\mathtt{POLE}/L$)。$K_y^2 = K_z^2 + K_x^2$、$Q_x^2 = K_z^2 + Q_y^2$ のとき $\nabla\times\vec B = 0$ ($\nabla\cdot\vec B = 0$ は恒等)。[und_01_field_map.out.txt](outputs/und_01_field_map.out.txt) (a): sympy の $\nabla\times\vec B$ を数値代入して最大 $8\times10^{-15}$ (B ~ 0.3)。
- **UND の写像**: $H = [(p_x - A_x)^2 + (p_y - A_y)^2]/(2(1+\delta))$ (傍軸)。スライスごとに、2 型の母関数 $F_2 = q\cdot P + ds\,H(q, P; s_j)$ の陰的な写像 (記号は §3.3)。**1 次の精度**で厳密にシンプレクティック。確認 (b): 3 GeV、$B_y = 1$ T、$B_x = 0.4$ T、Kx・Qy・位相あり、10 周期・1 m で、DOP853 の積分との $p_x$ の誤差は slice = 1000・2000・4000・8000 で $1.07\times10^{-11}$、$7.07\times10^{-12}$、$3.97\times10^{-12}$、$2.09\times10^{-12}$ (比 1.51、1.78、1.90、→ 2)。1 次。
- **K64 の UND の癖**: 磁場のキーワードが設定できない (B37、v2 だけ `FBX`・`FBY`)。陰的な刻みの右辺に `dpds` を二重にかける (B78、v2 は 1 回)。**$-dv\,L$ の項 (速度の差による $z$ の更新) がない** (K64 の誤り。gs で確認、v2 は修正 = B334、付録 B): 磁場なしの UND と DRIFT で $z$ が $-2.897\times10^{-11}$ m ($= dv\,L$、3 GeV・$\delta = 10^{-3}$・$L = 1$ m) 違う。
- **BEAMBEAM の力**: 強ビーム (2 次元ガウス、$\sigma_x,\sigma_y$) が弱粒子に与えるキックは Bassetti-Erskine の式。[bb_01_kick_slice.out.txt](outputs/bb_01_kick_slice.out.txt) (a): 丸い・$\sigma_x > \sigma_y$・$\sigma_x < \sigma_y$ の 3 通り各 4 点で、$\Delta p_x$、$\Delta p_y$ が 2 次元ポテンシャルの数値積分の微分と $10^{-13}$ 以下の絶対差 (相対 $10^{-9}$) で一致。小振幅の傾き $-2c/(\sigma_x(\sigma_x+\sigma_y))$ ($c = N r_e/\gamma$) とも一致。
- **エネルギーのキック**: $\delta$ の変化 $= -g_{xy} + \tfrac14|\Delta p|^2$ ($g_{xy} = \sigma_{12x}\,\partial U/\partial\sigma_x^2 + \sigma_{12y}\,\partial U/\partial\sigma_y^2$) が、$\sigma_{12}\ne0$ (Twiss $\alpha\ne0$) の強ビームで $10^{-18}$ 以下の差で一致 (3 通り各 4 点)。
- **6 次元写像**: スライスごとに、衝突点まで drift ($x \mathrel{+}= p_x s^*$、$\delta \mathrel{-}=p_\perp^2/4$、$s^* = (z - z_s)/2$)、キック、同じ $s^*$ で戻す。$F = s^*(z)\,p_\perp^2/2$ の流れの共役で、厳密にシンプレクティック。確認 (b): SLICE = 5・交差角 10 mrad・軌道のずれ・傾いた強ビームで、$\|J^\top SJ - S\|$ (規格化) $\le 1.5\times10^{-11}$。キックなし (NP = 0) の往復は恒等 ($\le10^{-19}$)。
- **ビームストラールング**: 1 スライスの光子数 $N_\gamma = \mathtt{cnbs}\,\rho_i\,dsz = a_{\rm nrad}\gamma\theta$ (`cnbs` $= a_{\rm nrad}\gamma$ が比 $1.000000000$)、光子エネルギー $u = \mathtt{cubs}\,\rho_i\,\xi$ (`cubs` $= c_{\rm uc}\gamma^2$、比 $0.99999997$)。$\xi$ は光子数のスペクトルの逆引き表。表の $\langle\xi\rangle = 0.3128$、$\langle\xi^2\rangle = 0.4180$ は理論値 $8/(15\sqrt3) = 0.3079$、$11/27 = 0.4074$ に対して $+1.6\%$、$+2.6\%$ 高い (§9、疑い (小))。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| $A_x, A_y$ | UND のベクトルポテンシャル ($B\rho$ で規格化、1/m の次元の積分 ≙ 無次元) | `ff`、`gg` (`undulator.f:60-64`) | `und._und_a` |
| $F_0, G_0$ | 軸上の $B_y$・$B_x$ の振幅 / $B\rho$ [1/m] | `ulist(i_F0)`、`ulist(i_G0)` (`undinit`: `p_in(2)/pgev*c`) | `compile_und` の `f0`、`g0` |
| $K_z$ | $2\pi\,\mathtt{POLE}/L$ | `ulist(i_Kz)` (`len`、`npole`) | `kz` |
| $K_x, K_y, Q_x, Q_y$ | 横方向の波数 | `i_Kx`、`i_Ky`、`i_Qx`、`i_Qy` | `kx`、`ky`、`qx`、`qy` |
| $\varphi$ | $G_0$ の項の位相 (`PHI` [度]) | `i_dph` | `dph` |
| $ds$ | 1 スライスの長さ $L/\mathtt{SLICE}$ | `ulist(i_ds)`、`nslice` | `ds`、`nslice` |
| `dpds` | $ds/(1+\delta)$ | `dpds` | `dpds` |
| $\sigma_x^2, \sigma_{12x}$ ... | 強ビームのスライスの包絡 (5×5 の `Benv5`) | `colb%Benv5` | `bb['benv5']` |
| $z_s$ | スライスの重心 ($\sigma_z$ 単位) | `colb%zslice(1:nsli)` | `zc` |
| $c$ | $N r_e/\gamma$ | `colb%ce` ($\mathtt{rne}\cdot re/\gamma$) | `ce` |
| `nslice` | z のスライス数 (SLICE) | `colb%nslice` | `nslice` |
| $\varphi_c$ | 交差角 (XANGLE) | `colb%xangle(1:4)` = 角、$\sin$、$\cos$、$\tan$ | `xangle` |
| $\rho_i$、$dsz$ | 1/曲率半径、スライスの半長 | `rhoi`、`dsz` | |

粒子の変数は SAD の $(x, p_x, y, p_y, z, \delta)$。追跡の内部で `px = px*pn` ($p_n = 1+\delta$) と正準運動量にして計算し、最後に `px/pn` に戻す (UND も BEAMBEAM も同じ)。

## 3. UND: ハミルトニアンとスライスの写像

### 3.1 ベクトルポテンシャルの模型

直線の要素 ($1/\rho = 0$)、電場なし、縦のポテンシャル $A_s = 0$ のゲージ。K64 は

$$
\begin{aligned}
A_x &= \frac{F_0}{K_z}\cos(K_xx)\cosh(K_yy)\sin(K_zs) - \frac{Q_yG_0}{Q_xK_z}\sinh(Q_xx)\sin(Q_yy)\sin(K_zs+\varphi),\\
A_y &= -\frac{G_0}{K_z}\cos(Q_yy)\cosh(Q_xx)\sin(K_zs+\varphi) + \frac{K_xF_0}{K_yK_z}\sinh(K_yy)\sin(K_xx)\sin(K_zs)
\end{aligned}
\tag{16.1}
$$

(`ff` = $A_x$、`gg` = $A_y$。$B\rho$ で規格化済み。) 磁場は $\vec B = \nabla\times\vec A$ ($A_s = 0$): $B_x = -\partial_sA_y$、$B_y = \partial_sA_x$、$B_z = \partial_xA_y - \partial_yA_x$。第 1 項 ($F_0$) の寄与を計算すると

$$
\vec B_{F_0} = F_0\left(-\frac{K_x}{K_y}\sin(K_xx)\sinh(K_yy)\cos(K_zs),\ \cos(K_xx)\cosh(K_yy)\cos(K_zs),\ -\frac{K_z}{K_y}\cos(K_xx)\sinh(K_yy)\sin(K_zs)\right)
\tag{16.2}
$$

($B_z$ の計算: $\partial_xA_y - \partial_yA_x$ の $F_0$ の項は $\frac{F_0K_x^2}{K_yK_z}\sinh\cos\sin - \frac{F_0K_y}{K_z}\cos\sinh\sin = -\frac{F_0}{K_yK_z}(K_y^2 - K_x^2)\cos\sinh\sin = -\frac{F_0K_z}{K_y}\cos(K_xx)\sinh(K_yy)\sin(K_zs)$、$K_y^2 - K_x^2 = K_z^2$ を使った。) 軸上は $B_y = F_0\cos K_zs$ ($K_x = 0$ なら普通の平面アンジュレータ: $B_y = F_0\cosh(K_zy)\cos(K_zs)$、$B_z = -F_0\sinh(K_zy)\sin(K_zs)$)。第 2 項 ($G_0$) は $x \leftrightarrow y$ の役割を入れ替えた水平磁場の項で、軸上 $B_x = G_0\cos(K_zs+\varphi)$。

**Maxwell 方程式**: 各成分は調和関数でなければならない。(16.2) の $B_y$ の $\nabla^2$ は $(-K_x^2 + K_y^2 - K_z^2)B_y$ なので、$K_y^2 = K_z^2 + K_x^2$ のとき 0。第 2 項は $Q_x^2 = K_z^2 + Q_y^2$ が条件。**確認** (a): sympy で (16.1) から $\nabla\times\vec B$ を作り、$K_y, Q_x$ をこの関係にして、数値代入 10 点で最大 $8.2\times10^{-15}$ ($B \sim 0.3$)。$\nabla\cdot\vec B$ は恒等的に 0。軸上の $\vec B = (G_0\cos(K_zs+\varphi),\ F_0\cos K_zs,\ 0)$。

**単位**: $F_0 = \mathtt{FBX}\,c/p$ ($\mathtt{FBX}$ [T]、$p$ [eV]、$c$ [m/s]) $= B_y/B\rho$ [1/m]。軌道の最大の角は軸上の $p_x = -A_x/(1+\delta)$ から $F_0/K_z$ ($= K/\gamma$、偏向パラメータ $K = eB\lambda_u/2\pi mc$)。**確認** (c): $F_0 = 0.3$ /m、$K_z = 62.83$ /m で $F_0/K_z = 4.7746\times10^{-3}$、追跡の最大の $p_x$ (20 スライスおきの標本) $4.7741\times10^{-3}$。

**使える範囲**: 周期の数 (POLE) が整数、$\varphi = 0$ のとき $A(0) = A(L) = 0$ だが、磁場は端で最大 ($\cos$)。つまり**端の処理 (終端の極) が模型にない**。1 周期ごとに 1 次の積分 ($\int B_y\,ds$) は 0。$\varphi \ne 0$ のときは、入口と出口の $A_y \ne 0$ なので、**端でポテンシャルの飛びが蹴りになる** (確認 (b): $A_y = -1.9\times10^{-4}$)。

### 3.2 ハミルトニアン

$s$ を独立変数にした SAD のハミルトニアン ([01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)) は、$A_s = 0$ の直線で $H = -\sqrt{(1+\delta)^2 - (p_x-A_x)^2 - (p_y-A_y)^2}$ (定数の $E/(p_0v_0)$ の項は別)。運動学的な横運動量 $\pi_{x,y} = p_{x,y} - A_{x,y}$ で 2 次まで展開すると、$H = -(1+\delta) + \dfrac{\pi_x^2 + \pi_y^2}{2(1+\delta)} + O(\pi^4/p^3)$。定数部分を除いて

$$
H_{\rm UND} = \frac{(p_x - A_x)^2 + (p_y - A_y)^2}{2(1+\delta)}
\tag{16.3}
$$

(傍軸近似。落とした項は $\pi^4/(8p^3)$ 程度。) 運動方程式は正準方程式から: $x' = \pi_x/p_n$、$y' = \pi_y/p_n$、$p_x' = -\partial H/\partial x = (\pi_x\,\partial_xA_x + \pi_y\,\partial_xA_y)/p_n$、$p_y'$ も同様、$z' = \partial H/\partial\delta = -(\pi_x^2 + \pi_y^2)/(2p_n^2)$、$\delta' = 0$ ($p_n = 1+\delta$)。(SAD の $(z, \delta)$ は $(q, p)$ の対で、$z' = \partial H/\partial\delta$。)

### 3.3 スライスの写像

1 スライス ($ds = L/\mathtt{SLICE}$、$j$ 番目、$s_j = j\,ds$) で、K64 は次の写像を使う ($\vec q = (x,y)$、$\vec p = (p_x, p_y)$):

$$
\vec P = \vec p - ds\,\partial_{\vec q}H(\vec q, \vec P;\ s_j),\qquad
\vec q' = \vec q + ds\,\partial_{\vec P}H(\vec q, \vec P;\ s_j),\qquad
z' = z + ds\,\partial_\delta H(\vec q, \vec P;\ s_j),\qquad \vec p' = \vec P
\tag{16.4}
$$

($H$ の中の $A$ は古い位置 $\vec q$ と終点の $s_j$ で評価し、運動量は新しい $\vec P$。) これは 2 型の母関数 $F_2(\vec q, \vec P) = \vec q\cdot\vec P + ds\,H(\vec q, \vec P; s_j)$ から $\vec p = \partial F_2/\partial\vec q = \vec P + ds\,\partial_{\vec q}H$、$\vec q' = \partial F_2/\partial\vec P = \vec q + ds\,\partial_{\vec P}H$ と得られる、**正準変換** (ds によらずシンプレクティック。$(z, \delta)$ の項は $\delta' = \delta$ が保たれ、$z$ の更新は同じ $F_2$ の $\delta$ 微分)。精度は 1 次 (シンプレクティック Euler)。$H$ が $\vec P$ の 2 次なので陰的な式は線形で、2×2 の行列で解ける:

$$
\begin{pmatrix}1 - \hat f_x & -\hat g_x\\ -\hat f_y & 1 - \hat g_y\end{pmatrix}\begin{pmatrix}P_x\\P_y\end{pmatrix} = \begin{pmatrix}p_x - \widehat{fg}_x\\ p_y - \widehat{fg}_y\end{pmatrix},\quad
\hat f_x = \frac{ds}{p_n}\partial_xA_x,\ \hat g_x = \frac{ds}{p_n}\partial_xA_y,\ \ldots,\ \widehat{fg}_x = \hat f_xA_x + \hat g_xA_y
\tag{16.5}
$$

(最後の式は、$\partial_xH = -\frac1{p_n}[(P_x - A_x)\partial_xA_x + (P_y - A_y)\partial_xA_y]$ を $P_x(1-\hat f_x) - \hat g_xP_y = p_x - (\hat f_xA_x + \hat g_xA_y)$ に整理したもの。) 行列式 $\det_{1a} = (1 - \hat f_x)(1 - \hat g_y) - \hat f_y\hat g_x$。`undulator.f` の `fx, gx, fy, gy` に `dpds` を掛け、`det1a`、`fx = (1-gy)/det1a` などで解く (`undulator.f:84-103`)。更新: $x' = x + dpds\,(P_x - A_x)$、$y' = y + dpds\,(P_y - A_y)$、$z' = z - \tfrac12(\pi_x^2 + \pi_y^2)\,dpds^2/ds = z - \dfrac{ds}{2p_n^2}\pi^2$ (`undulator.f:105-109`)。

**入口・出口**: 外側 ($A = 0$) の運動量は正準運動量に等しい。入口では外の $p_x p_n$ を正準運動量として持ち込み、内側の力学的運動量は $p - A(0)$ になる。出口では最後の正準運動量 $/p_n$ を外の運動量として返す ($A(L) = 0$ でないとその分が蹴りになる)。**確認** (b): 入口の $A_y = -1.88\times10^{-4}$、出口も同じ値 (位相 $\varphi = 0.3$) のとき、DOP853 の積分も同じ約束 (入口で正準運動量 = 外の運動量) で行って、上の誤差になる。

**K64 の誤り (B78)**: `fgx = dpds*(ff*fx + gg*gx)` と、すでに `dpds` を掛けた `fx, gx` にもう 1 回 `dpds` を掛ける (`undulator.f:92-93`)。(16.5) の右辺の $\widehat{fg}$ が $ds/p_n$ 倍になり、収束力が刻み幅に比例して消える。v2 は 1 回だけ掛ける (本章の確認は v2 の式)。

### 3.4 刻み数

1 次の写像なので、誤差は刻み数に反比例。確認 (b) の例 (周期 0.1 m、$K \sim 9.4$〜) で、$p_x$ の誤差が $2\times10^{-12}$ (8000 スライス = 800 スライス/周期) に下がる。K64 の既定の目安は書かれていない ([docs/19_other_elements.md](../docs/19_other_elements.md) の例は 2000 スライス/20 極)。

## 4. UND の確認

(確認 (a)〜(d) は §1・§3.1・§3.3 の数字。) (d): 磁場なしの UND と DRIFT ($L = 1$ m、3 GeV、$x_0 = 10^{-4}$、$p_{x} = 2\times10^{-4}$、$\delta = 10^{-3}$): $x$ は $2.998001998\times10^{-4}$ (UND、傍軸 $x' = p_x/p_n$) と $2.998002038\times10^{-4}$ (DRIFT、$x' = p_x/p_z$)、差 $4\times10^{-12}$。$z$ は $-1.99600599\times10^{-8}$ と $-1.99310907\times10^{-8}$、差 $-2.897\times10^{-11}$ m。これは $dv\,L$ ($dv = -\delta(2+\delta)/[h_1(h_1 + p_nh_0)]$、$h_1 = \sqrt{1 + (p_np_0)^2}$、$L = 1$) の値と一致する: UND の $z$ の更新に $-dv\,L$ の項がない (DRIFT は持つ、[02_drift.md](02_drift.md))。$R_{56} = L/\gamma^2$ が欠ける ($\gamma = 5871$ で $\sim 3\times10^{-8}\,L$)。高エネルギーでは小さいが、低エネルギーでは大きな誤差になる (疑い、§13)。

## 5. BEAMBEAM の構成

BEAMBEAM は**弱-強モデル**: 強ビーム (固定、6 次元ガウス分布、`Benv` = 6×6 の包絡) が弱粒子に与えるキックだけを計算する。要素 1 個の内部は

1. 交差角の変換 (§8)、軌道のずれの補正、傾きの回転。
2. 強ビームを $z$ に沿って等確率の `SLICE` 個のスライスに分け (§7)、各スライスについて: 衝突点まで drift → キック → 戻す (§6)。
3. (ビームストラールング。§9)
4. 逆の変換。

光学 (Twiss、CALC、Emittance) には効かない (K64 の `qtwiss` に分岐がない。[docs/14_beambeam.md](../docs/14_beambeam.md) §4.3)。

## 6. 1 スライスの写像

### 6.1 キックの力: 2 次元ガウス分布の電場

強ビームのスライスは、電荷 $N/n_{\rm slice}$ の 2 次元のガウス分布 ($\sigma_x, \sigma_y$、中心は原点)。弱粒子 (電荷が反対) のキック ($p_0$ で規格化した運動量の変化) は、2 次元ポテンシャル $U$ の勾配で

$$
\Delta p_{x} = -c\,\frac{\partial U}{\partial x},\quad \Delta p_y = -c\,\frac{\partial U}{\partial y},\qquad c = \frac{Nr_e}{\gamma\,n_{\rm slice}},\quad
U = \int_0^\infty\frac{1 - \exp\!\left(-\frac{x^2}{2\sigma_x^2+t} - \frac{y^2}{2\sigma_y^2+t}\right)}{\sqrt{(2\sigma_x^2+t)(2\sigma_y^2+t)}}\,dt
\tag{16.6}
$$

($c$ は `ce/nslice`、$\mathtt{ce} = N r_e/\gamma$。) (16.6) の導出の要点: 2 次元の Poisson 方程式の解 (ガウス分布のポテンシャル) を、ガウス分布の積分表示 ($e^{-r^2/(2\sigma^2)}$ の Schwinger 表示) で書いたもの。**丸いビームの検算**: $\partial U/\partial x = \int_0^\infty\frac{2x\,e^{-r^2/(2\sigma^2+t)}}{(2\sigma^2+t)^2}dt = \frac{2x}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)$ ($u = 1/(2\sigma^2+t)$ と置換)。したがって

$$
\Delta\vec p_\perp = -\frac{2Nr_e}{\gamma n_{\rm slice}}\,\frac{\vec r}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)\qquad(\sigma_x = \sigma_y)
\tag{16.7}
$$

(コードの丸いビーム: `a2 = 2*(1-exp(hi))/(x**2+y**2)`、`fx0 = c1*a2*x`。`hi` $= -x^2/2\sigma_x^2 - y^2/2\sigma_y^2$。) 小振幅 $r \ll \sigma$ で $-\dfrac{N r_e}{\gamma\sigma^2}\vec r$ ($\sigma_x\ne\sigma_y$ でも $\Delta p_x = -2c\,x/(\sigma_x(\sigma_x+\sigma_y))$、ビーム・ビームパラメータ $\xi_x = N r_e\beta_x^*/(2\pi\gamma\sigma_x(\sigma_x+\sigma_y))$ の元)。

**一般の楕円 (Bassetti-Erskine)**: (16.6) の閉じた形は、複素誤差関数 $w(\zeta) = e^{-\zeta^2}\mathrm{erfc}(-i\zeta)$ で書かれる。$\sigma_x > \sigma_y$ のとき、$S = \sqrt{2(\sigma_x^2 - \sigma_y^2)}$、$\zeta_1 = (|x| + i|y|)/S$、$\zeta_2 = (|x|\sigma_y/\sigma_x + i|y|\sigma_x/\sigma_y)/S$ として

$$
f_y + i f_x = \frac{2\sqrt\pi\,c}{S}\left[w(\zeta_1) - e^{h}\,w(\zeta_2)\right],\qquad h = -\frac{x^2}{2\sigma_x^2} - \frac{y^2}{2\sigma_y^2},\quad \Delta p_x = -\mathrm{sgn}(x)\,f_x,\ \Delta p_y = -\mathrm{sgn}(y)\,f_y
\tag{16.8}
$$

(コード: `cerrf(d1x,d1y,w1y,w1x)` で $w_{1y} = \mathrm{Re}\,w$、$w_{1x} = \mathrm{Im}\,w$ を得て、`fx0 = c1*(w1x - a2x)`、`fy0 = c1*(w1y - a2y)`、`c1 = ce*2*sqrpi/S/nslice`。象限の符号は $|x|, |y|$ を使うので、最後に符号を戻す。) $\sigma_x < \sigma_y$ は $x \leftrightarrow y$ を入れ替えた式 (コードの 3 つめの枝)。(16.8) は標準の Bassetti-Erskine の式 (導出は引用)。**確認** (a): 丸い・$\sigma_x > \sigma_y$ ($0.2\ {\rm mm}$、$0.05\ {\rm mm}$)・$\sigma_x < \sigma_y$ (逆) の 3 通り、各 4 点 ($x$ が 0.2〜3 $\sigma$) で、(16.6) を scipy の 2 重の数値積分で微分した結果と、絶対差 $\le 1.4\times10^{-13}$ (キックの大きさ $\sim5\times10^{-5}$、相対 $3\times10^{-9}$) で一致。

### 6.2 $z$ に依存する包絡と $\delta$ のキック

強ビームのスライスは、衝突点の位置 $s^*$ により、衝突時の横の大きさが違う (包絡が $s^*$ の 2 次関数で、$\sigma_x^2(s) = \sigma_{11x} + 2\sigma_{12x}s + \sigma_{22x}s^2$ (`sig11x = Benv5(1) + 2*Benv5(2)*sz + Benv5(7)*sz**2`、$\sigma_{12x} = \mathtt{Benv5(2)} + \mathtt{Benv5(7)}\,s^*$))。したがってポテンシャルは $z$ にも依存し、運動量 $\delta$ が変わる。物理的なポテンシャルは、遠方で $\sigma$ に依らず $\ln r^2$ になるもの $U_{\rm phys} = U + 2\ln(\sigma_x + \sigma_y)$ ((16.6) の $U$ は遠方で $\ln r^2 - 2\ln(\sigma_x+\sigma_y)+$定数、定数部分は $\sigma$ の関数なので足す)。**$\delta$ のキックは $\delta \mathrel{-}= g_{xy}$**、

$$
g_{xy} = c\left(\sigma_{12x}\,\frac{\partial U_{\rm phys}}{\partial\sigma_x^2} + \sigma_{12y}\,\frac{\partial U_{\rm phys}}{\partial\sigma_y^2}\right)\quad\left(\text{丸いビーム: }\frac{\partial U_{\rm phys}}{\partial\sigma^2} = \frac{e^{h}}{\sigma^2},\ g_{xy} = c\,\frac{\sigma_{12x}}{\sigma_{11x}}e^{h}\right)
\tag{16.9}
$$

($\partial U_{\rm phys}/\partial\sigma^2$ は丸い場合 $E_1$ の微分 $e^{-w}/\sigma^2$。コード: `gxy = c1*sig12x/sig11x*exp(hi)`、楕円は `dux, duy` を使い `gxy = sig12x*dux + sig12y*duy`。)

**導出**: 相互作用 $H_{\rm int} = U_{\rm phys}(x, y; \sigma(s^*))\,\delta_D(s - s^*)$ の $\delta$ に関する運動方程式は $\delta' = -\partial H/\partial z$。衝突点 $s^* = (z - z_s)/2$ なので $\partial s^*/\partial z = \tfrac12$、$\partial\sigma_x^2/\partial s^* = 2\sigma_{12x}(s^*)$ ($\sigma_{12x}(s^*) = \sigma_{12x} + \sigma_{22x}s^*$、`sig12x`) から $\partial U/\partial s^* = 2\sigma_{12x}(s^*)\,\partial U/\partial\sigma_x^2 + 2\sigma_{12y}(s^*)\,\partial U/\partial\sigma_y^2$。よって $\delta \mathrel{-}= \tfrac12\,\partial U/\partial s^*\cdot c = g_{xy}$。**確認** (a): 3 通りの強ビーム ($\alpha_x = 0.3, 0.5, -0.2$ など、$\sigma_{12} = -\alpha\epsilon$) で、1 スライス・$z = 0$ ($s^* = 0$) の $\Delta\delta$ が $-g_{xy} + \tfrac14|\Delta p_\perp|^2$ ($g_{xy}$ は $U_{\rm phys}$ の $\sigma^2$ 微分を数値積分で作ったもの、$\tfrac14|\Delta p|^2$ は §6.3 の戻しの項) と絶対差 $\le5\times10^{-18}$ (値 $\sim10^{-9}$) で一致する。$2\ln(\sigma_x+\sigma_y)$ の項を入れないと一致しない (丸い場合に符号さえ違う)。

### 6.3 衝突点までの drift と戻し

粒子の位置 $z$ と、スライスの位置 $z_s$ ($\sigma_z$ 単位の重心 $\times\sigma_z$) の相手が、2 つのビームが対向して動くので、衝突点は $s^* = (z - z_s)/2$ (相対速度 $2c$)。粒子は衝突点まで drift してキックを受け、同じ距離だけ戻る:

$$
M_{\rm slice} = D^{-1}\circ K\circ D,\qquad D:\ x \mathrel{+}= p_x s^*,\ y \mathrel{+}= p_ys^*,\ \delta \mathrel{-}= \tfrac14(p_x^2 + p_y^2)
\tag{16.10}
$$

($D$ の中の $s^* = (z - z_s)/2$ は D の前の $z$。$D^{-1}$ は、同じ $F$ の流れの**時間 $-1$** を、キック後の $p$ で行う: $x \mathrel{-}= p_xs^*$、$\delta \mathrel{+}= \tfrac14p_\perp^2$ (同じ $s^*$。$K$ は $z$ を変えないので $s^*$ は変わらない)。$D$ の単純な逆写像 (古い $p$ を使う) ではなく、$p$ が $K$ で変わった分だけ違う。$K$: $p_{x,y} \mathrel{-}= f_{x,y}$、$\delta \mathrel{-}= g_{xy}$。)

**$D$ の導出**: 生成関数 $F(x, y, p_x, p_y, z) = s^*(z)\,H_0$、$H_0 = (p_x^2 + p_y^2)/2$ ($s^* = (z - z_s)/2$) の正準方程式を「時間 1」流す。$\dot x = \partial F/\partial p_x = s^*p_x$、$\dot p_x = -\partial F/\partial x = 0$、$\dot z = \partial F/\partial\delta = 0$ (傍軸の $H_0$ は $\delta$ に依らない)、$\dot\delta = -\partial F/\partial z = -\tfrac12H_0 = -\tfrac14p_\perp^2$。$F$ は自分の流れで保存するので、これは厳密に $x \mathrel{+}= s^*p_x$、$\delta \mathrel{-}= \tfrac14p_\perp^2$ ($s^*$ は一定)。したがって $D = e^{:F:}$ 型の写像で、$D^{-1}K D$ は共役なのでシンプレクティック。**確認** (b): SLICE = 5、交差角 10 mrad、軌道のずれ、$\alpha\ne0$ の強ビームで、6×6 ヤコビアン (中心差分) の $\|J^\top SJ - S\|$ (規格化) は $1.5\times10^{-11}$、傾いた (x-y 結合のある) 強ビームで $4.3\times10^{-12}$、SLICE = 1・交差角なしで $7.6\times10^{-12}$。行列式は $1 \pm 5\times10^{-10}$。

## 7. 強ビームの包絡とスライス

### 7.1 6 次元の包絡

強ビームの 6×6 の 2 次モーメント $\Sigma$ ($x, p_x, y, p_y, z, \delta$) を、Twiss パラメータから作る (`bbinit`、`beambeam.setup`): $\Sigma = (HRB)\,E\,(HRB)^\top$、$B$ は $(\alpha, \beta)$ ($x$、$y$、$z$)、$R$ は x-y 結合 (R1..R4)、$H$ は分散 (EX, EPX, EY, EPY, ZX, ZPX, ZY, ZPY)、$E = \mathrm{diag}(\epsilon_x, \epsilon_x, \epsilon_y, \epsilon_y, \sigma_\delta\sigma_z, \sigma_\delta\sigma_z)$。縦: $\epsilon_z = \sqrt{\sigma_\delta^2\sigma_z^2 - \langle z\delta\rangle^2}$、$\beta_z = \sigma_z^2/\epsilon_z$、$\alpha_z = -\langle z\delta\rangle/\epsilon_z$。`R11..R66` がすべて正なら、直接 $\Sigma$ を与える。交差角があれば $\Sigma \to C\Sigma C^\top$ (§8)。

### 7.2 スライスの条件付き包絡

$z$ を固定したときの $(x, p_x, y, p_y, \delta)$ の条件付きの分布 ($z$ 軸の方向に重みをつけたスライス) は、ガウス分布の性質から、**精度行列 $\Sigma^{-1}$ の対応する 5×5 のブロックの逆行列**が共分散 $\Sigma_5$、条件付き平均が $z\cdot(-\Sigma_5\,b)$ ($b$ = $\Sigma^{-1}$ の $z$ 行・5 元の列)。K64: `m_zslice_5d` が `Benv_inv` から `m_sub5`、`v_sub5` を作り、固有値分解で $\Sigma_5 = \mathrm{Benv5}$ を作る (固有値が 0 以下ならエラー)。`v_cen` は $\langle x_i\rangle = v_i\,z$ ($i = 1\ldots4$)。このずれ (`deform`) は、各スライスで `x -= zslice*sigz*v_cen(1)` などの座標の平行移動として入れる ((16.10) の $D$ の前後)。

x-y 結合 (`Benv5(3)` の非対角) は、傾きの回転 (`benv_norm_axis`、角 $t = \tfrac12\arctan2(2\Sigma_{13}, \Sigma_{11} - \Sigma_{33})$) で消す。**$\sigma_x < \sigma_y$ の強ビームは $t = \pi/2$ になる**: 包絡が 90 度回され、$\sigma_x > \sigma_y$ の枝で計算される。K64 は傾きの回転に `sint = xangle(5)` ($= t$ そのもの)、`cost = xangle(6)` ($= \sin t$) を使う (`bbstrhl.f:220-221`、`xangle(6)=sin`、`xangle(7)=cos` と書いてあるのに読む位置が 1 つずれる、B34)。**確認** (a): $\sigma_x < \sigma_y$ で k64 モードは $\Delta p$ が v2 と違う: $(x, y) = (2\times10^{-5}, -10^{-5})$ で $\Delta p_x$ = $-7.9\times10^{-5}$ (k64)、$-2.6\times10^{-5}$ (v2、(16.6) と一致)。$\sigma_x > \sigma_y$ のビーム (傾き 0) は影響を受けない。

### 7.3 z のスライス

等確率の $n$ 個のスライスの境界は $z_i = \sigma_z\,\Phi^{-1}(i/n)$ (`gauinv`。K64 は最大相対誤差 $6.4\times10^{-4}$ の近似式)。各スライスの重心 (σ_z 単位) は、ガウス分布の部分区間の平均

$$
\bar z_i = \frac{n}{\sqrt{2\pi}}\left(e^{-z_{i-1}^2/2} - e^{-z_i^2/2}\right)
\tag{16.11}
$$

(コードは符号を反転して「前から後ろ」の順に並べる: `colb%zslice(i) = -...`。両端のスライスは $z_0 = -\infty$、$z_n = +\infty$ に当たる。) ルミノシティの計算 `dLum` は各スライスの 2 次元ガウス密度 $\exp(-x^2/2\sigma_x^2 - y^2/2\sigma_y^2)/(2\pi\sigma_x\sigma_y)$ の和 ([docs/14_beambeam.md](../docs/14_beambeam.md))。

## 8. 交差角・軌道のずれ・傾き

**軌道のずれ**: 強ビームの軌道 $(x_c, p_{xc}, y_c, p_{yc}, z_c)$ を、弱粒子の座標から引いて (入口)、あとで足す (出口): `x -= cod(1)`、`px += cod(2)`、... (`bbstrhl.f:205-215`)。強ビームの軌道の $p$ は $\cos\varphi_c$ で割って入れる (`cod(2) = DPX/cos`)。

**交差角**: 弱ビームの座標から、強ビームが正面衝突する座標系への 6 変数の閉じた写像 ($\varphi_c$ = XANGLE、$p_n = 1+\delta$、$p_z = \sqrt{p_n^2 - p_x^2 - p_y^2}$、$h_1 = p_n - p_z$ (= `pxy2/(pn+pz)`))。

$$
\begin{aligned}
\delta' &= \delta - \tan\varphi_c\,p_x + \tan^2\varphi_c\,h_1,\qquad p_x' = \frac{p_x - \tan\varphi_c\,h_1}{\cos\varphi_c},\qquad p_y' = \frac{p_y}{\cos\varphi_c},\\
p_z' &= \left[(1+\delta')^2 - p_x'^2 - p_y'^2\right]^{-1/2}\ (\text{逆数}),\quad
x' = \tan\varphi_c\,z + \left(1 + \sin\varphi_c\,p_x'p_z'\right)x,\\
y' &= y + \sin\varphi_c\,p_y'p_z'\,x,\qquad z' = \frac{z}{\cos\varphi_c} - \frac{\tan\varphi_c}{\cos\varphi_c}\,h_1\,p_z'\,x
\end{aligned}
\tag{16.12}
$$

(`bbstrhl.f:184-201`。ここで `pz` は $1/p_z'$ を指す。出口で逆の写像 (`:482-498`)。) 小振幅・小角の極限で $x' \approx x + \varphi_cz$、$z' \approx z$、$\delta' \approx \delta - \varphi_cp_x$、$p_x' \approx p_x - \varphi_c\delta$ (座標系の回転 + Lorentz boost の合成に当たる)。この写像は**ソースから写したもので、導出していない**。**確認** (b): ここでの 6×6 ヤコビアンのシンプレクティック性 (上の数字)、キックなしの往復が恒等 ($|\text{出}-\text{入}| \le 10^{-19}$ (x, δ の 1e-19 は丸め)、$z = 0.01$ の粒子が $z = 0.01$ に戻る)。

**傾き**: 強ビームの x-y 結合を除く角 $t$ で、座標と運動量を回す: $x \to cx + sy$、$y \to cy - sx$ ($p$ も同様)。K64 は $c, s$ の読み位置を誤る (B34、k64 は $t \neq 0$ で誤り。v2 は $\cos t$、$\sin t$)。

## 9. ビームストラールング

BSTRL ≥ 1 のとき、各スライスのキックの大きさ $\theta = \sqrt{f_x^2 + f_y^2}$ ($= \mathtt{fx0}, \mathtt{fy0}$) から、曲率 $\rho_i = \theta/dsz$ ($dsz$ = スライスの半長 $= (z_{i-1} - z_i)\sigma_z/2$、`dsz`) とする。放射は [13_radiation.md](13_radiation.md) の kinematical method の考え方で、1 スライス 1 個以下の光子として扱う:

$$
N_\gamma = \mathtt{cnbs}\,\rho_i\,dsz = a_{\rm nrad}\,\gamma\,\theta,\qquad
\mathtt{cnbs} = \frac{5\sqrt3}{6}\alpha\gamma = a_{\rm nrad}\gamma,\qquad
u_c/E_0 = \mathtt{cubs}\,\rho_i,\quad \mathtt{cubs} = \frac{3}{2}\frac{\hbar c\gamma^2}{m_ec^2} = c_{\rm uc}\gamma^2
\tag{16.13}
$$

($5\sqrt3/6 = 5/(2\sqrt3)$ なので $\mathtt{cnbs} = a_{\rm nrad}\gamma$、$c_{\rm uc} = \tfrac32\bar\lambda_C$。**確認** (c): 比 $1.000000000$ と $0.99999997$。) 光子を出すかは乱数 $x_r\in(0,1)$ で $x_r < N_\gamma$ (`rnbs`)。$N_\gamma > 0.5$ は K64 が停止する ("Nbs is too large")。出すなら、$x_r/N_\gamma$ ($(0,1)$ に一様、乱数の再利用) を、光子数のスペクトル $dN/d\xi \propto \int_\xi^\infty K_{5/3}(x)dx$ ($\xi = u/u_c$) の累積分布の逆引き表 (`besseltab`: $\xi_i = 1.1^{40}/1.02^i$、累積 `SwI`) で $\xi_i$ に直し、$\delta \mathrel{-}= \mathtt{cubs}\,\rho_i\,\xi$。この分布の平均と 2 次のモーメントは $\langle\xi\rangle = c_{\rm ave}$、$\langle\xi^2\rangle = c_{\rm uu}$ に等しいはず ([13_radiation.md](13_radiation.md))。

**平均の損失**: $N_\gamma\,c_{\rm ave}\,u_c/E_0 = a_{\rm nrad}c_{\rm ave}c_{\rm uc}\gamma^3\theta^2/dsz = \tfrac23r_e\gamma^3\theta^2/dsz$ ($a_{\rm nrad}c_{\rm ave}c_{\rm uc} = \tfrac{8}{15\sqrt3}\cdot\tfrac{5\alpha}{2\sqrt3}\cdot\tfrac32\bar\lambda_C = \tfrac23r_e$)。古典の $\tfrac23r_e\gamma^3L/\rho^2$ ($\theta = L/\rho$) と一致 (`cpbs` = $\tfrac23r_e\gamma^3$ は定義されているが、平均の損失としては使わない)。

**確認** (c): 表を $x_r\in(0,1)$ の一様な 20 万点で逆引きして得たモーメントは $\langle\xi\rangle = 0.3128$、$\langle\xi^2\rangle = 0.4180$ (理論 $0.3079$、$0.4074$ に $+1.6\%$、$+2.6\%$)。同じ格子で累積を台形則で作り直して同じ手順にすると $0.3098$、$0.4098$ ($+0.6\%$、$+0.6\%$)。差は表の作り方 (矩形則に近い累積、2% 刻み) の誤差とみられる (§13 疑い (小)、影響は光子エネルギーの平均 ~1%)。

## 10. Fortran の対応

| 式 | ファイル:行 | サブルーチン | 変数・内容 |
|---|---|---|---|
| (16.1) | `undulator.f:54-64`、`:131-161` | `undulator`、`undinit` | `ff`、`gg`、`F0 = p_in(2)/pgev*c`、`Kz = 2π npole/len`、`Ky`、`Qx` |
| (16.4)、(16.5) | `undulator.f:66-109` | `undulator` | `fx, fy, gx, gy` (偏微分)、`det1a`、`fgx` (B78: `dpds` の二重)、`px/py/x/y/z` の更新 |
| 正準運動量 | `undulator.f:47-50`、`:116-122` | `undulator` | `px = px*pn`、最後に `px/pn`、`dv` |
| (16.6)〜(16.9) | `bbstrhl.f:285-361` | `beambeam` | `sig11x`、`sig12x`、`asp`、`cerrf`、`fx0`、`fy0`、`gxy`、`dux`、`duy` |
| (16.10) | `bbstrhl.f:250-257`、`:444-447` | `beambeam` | `sz = (z - zslice*sigz)/2`、`x += px*sz`、`g -= pxy2/4` (戻しは逆) |
| 強ビームの重心の移動 | `bbstrhl.f:262-267`、`:435-440` | `beambeam` | `v_cen` |
| (16.11) | `bbstrhl.f:776-803` | `bbinit` | `gauinv`、`zslice` |
| 包絡 (§7.1)、条件付き (§7.2) | `bbstrhl.f:629-663`、`:729-758` | `bbinit` | `Benv`、`m_zslice_5d`、`Benv5`、`v_cen`、`benv_norm_axis` |
| (16.12) | `bbstrhl.f:184-201`、`:482-498` | `beambeam` | `icross`、`tanx`、`sinx`、`cosx` |
| 傾き | `bbstrhl.f:219-230`、`:456-465` | `beambeam` | `sint = xangle(5)`、`cost = xangle(6)` (B34) |
| (16.13)、表 | `bbstrhl.f:365-430`、`:871-904` | `beambeam`、`mkbesseltab` | `cnbs`、`cpbs`、`cubs`、`rnbs`、`SwI`、`xi` |

(`bb.f` は同じ名前の旧い実装で、`tturn.f` が呼ぶのは `bbstrhl.f`。行番号は参照時点のもの。)

## 11. pysad の対応

| 関数 | 対応する Fortran | 注意 |
|---|---|---|
| `und.compile_und` | `undinit` | v2 は `FBX`・`FBY` を受ける。k64 は `F0 = G0 = 0` |
| `und.und_seg1`、`und.k_und_field` | `undulator` | v2: `dpds` は 1 回。`rad` のときは刻みごとに `tradk_core1` (K64 に放射はない)。`und.k_und` (k64) は磁場なしのドリフトの刻み |
| `beambeam.setup` | `bbinit` | 包絡、スライス、`_norm_axis` |
| `beambeam.bb_kick`、`apply` | `beambeam` | 傾きの読みは k64 は K64 のまま。v2 は B34 を直す。NaN は `limitnan` |
| `beambeam.cerrf` | `cerrf` (`bbtool.f`) | 複素誤差関数 (単精度の定数 3.2 が残る) |
| `beambeam.besseltab` | `mkbesseltab` | 添字は 1 始まりのまま (長さ 1002) |

## 12. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `und_01_field_map.py` | (a) $\nabla\times\vec B = 0$ (sympy)・軸上の磁場、(b) スライスの写像と DOP853 (誤差の次数)、(c) 最大の角 $F_0/K_z$、(d) UND と DRIFT の $z$ の差 $= dv\,L$ | [und_01_field_map.out.txt](outputs/und_01_field_map.out.txt) |
| `bb_01_kick_slice.py` | (a) Bassetti-Erskine のキックと $\delta$ のキック (ポテンシャルの数値積分)、傾きの k64・v2、(b) 6 次元写像のシンプレクティック性、往復の恒等、(c) ビームストラールングの定数と光子エネルギーの表 | [bb_01_kick_slice.out.txt](outputs/bb_01_kick_slice.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/und_01_field_map.py` (数秒)、`py -3.12 -W ignore physics_guide/scripts/bb_01_kick_slice.py` (数秒。scipy の積分の警告を抑える)。BEAMBEAM の gs との照合・光子リスト・ルミノシティ・ビーム・ビーム・チューンシフトの図は [docs/14_beambeam.md](../docs/14_beambeam.md)。

## 13. 未確認・限界・使った条件

**使った条件**: UND は 3 GeV、$L = 1$ m、10 周期、$B_y = 1$ T、$B_x = 0.4$ T (確認 (b))。BEAMBEAM は 1.7 GeV の電子 ($\gamma = 3326.8$)、$N = 10^{10}$、$\sigma = 0.05$〜0.2 mm、$\sigma_z = 5$〜10 mm、SLICE = 1〜5 (確認 (a)(b))。実機の衝突点 (β* が小さい、数百万粒子) の比較はしていない。

**未確認**:
- UND: 端の処理 (終端の極)、傍軸近似の誤差 ($\pi^4/8p^3$、$K$ が大きいとき)、放射 (v2 だけ)、スピンの結合 (v2 の POL)、`Pole` の数が整数でないとき。刻み数の選び方の指針。
- BEAMBEAM: **交差角の変換 (16.12) の導出** (Lorentz boost + 座標系の回転としての構成、(16.12) の各式の由来)。(16.8) の Bassetti-Erskine の閉じた形そのものの導出 (数値積分との一致だけ)。強-強の自己無撞着な計算 (弱-強だけ)。`cerrf` の精度 (単精度の定数 3.2)。スライス数の収束 (本章は SLICE = 1〜5 で式の確認のみ、ルミノシティ・チューンシフトの収束は未確認)。
- ビームストラールング: 光子エネルギーの表の精度 (§9)、$N_\gamma > 0.5$ の停止、光子 1 個の近似 (1 スライスに 2 個以上)、光子の放射角。
- SOL 領域内の BEAMBEAM は未対応。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| UND・BEAMBEAM | **ヘルプに記述がない** (SADHelp の要素の一覧・節に `UND`・`BEAMBEAM` の項目はなく、`UND` の文字列は `FOUND` などの一部だけ)。式・単位はすべて K64 のソースによる。 |
| ヘルプに無いもの | (16.1)〜(16.13) すべて。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

**確定 (独立な基準で確認済み。docs/30)**:
- **B37**: UND の磁場のキーワード (`FBX`・`FBY`・`BX`・`BY`) が K64 では設定できない (`kwFBX` が未定義)。そのため $F_0 = G_0 = 0$ の磁場なしの追跡になる。v2 は `FBX`・`FBY` を受ける。
- **B78**: `undulator.f:92-93` で `dpds` を二重に掛ける。磁場があっても収束力が刻み幅に比例して消える (200 と 2000 スライスで $p_y$ が 10 倍違う)。本章の確認 (b) は v2 の式で、DOP853 との一致 (1 次) を示す。
- **B34**: BEAMBEAM の傾きの回転が $c, s$ の読み位置を誤る (`sint = xangle(5)`、`cost = xangle(6) = sin t`)。**$\sigma_x < \sigma_y$ の強ビームは傾き $t = \pi/2$ になるので、k64 では $\Delta p$ が誤る** (本章確認 (a)、$-7.9\times10^{-5}$ 対 $-2.6\times10^{-5}$)。

- **B334**: UND の $z$ の更新に $-dv\,L$ (速度が $v_0$ と違うための遅れ、$R_{56} = L/\gamma^2$) がない (磁場なしの UND と DRIFT で $z$ が $dv\,L$ 違う、本章確認 (d))。DRIFT は持つ (02 章)。**gs で確認**: UND ($L$ = 1 m、KX = 0.1、QY = 0.2、POLE = 10、SLICE = 4、PHI = 30) と DRIFT の `TrackParticles` で、UND − DRIFT の $z$ が 3 GeV/c の $\delta = 10^{-3}$ で $-2.90\times10^{-11}$、$\delta = \pm10^{-2}$ で $\mp2.9\times10^{-10}$、10 MeV/c で $-2.60\times10^{-6}$、$\mp2.6\times10^{-5}$ ($= dv\,L$。自由飛行の式 $L(h_0/p_0 - h/p)$ と相対 $2\times10^{-3}$ で一致)。k64 は gs と全桁一致。v2 は $-dv\,L$ を足し、磁場なしの UND の $z$ が DRIFT と $10^{-14}$ で一致 (試験 `test_v250_und_dv.py`)。限界: 磁場なし (gs は磁場を設定できない) だけを gs で確認。v2 の磁場つきの $z$ の絶対値は未確認。

**疑い**:
- ビームストラールングの光子エネルギーの表: 累積の矩形則とみられる誤差で、平均が $+1.6\%$、2 次モーメントが $+2.6\%$ (理論比)。影響は小さい。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `ulist(1:13)` | UND のパラメータ (`nsli`、`npole`、`len`、`i_F0`、`i_G0`、`i_dph`、`i_Kz`、`i_Kx`、`i_Ky`、`i_Qx`、`i_Qy`、`lambda`、`i_ds`) | `undinit` | `compile_und` の戻り値 |
| `ff`、`gg` | $A_x$、$A_y$ | `undulator` | `_und_a` |
| `kx1..kx4` | $K_xx$、$K_yy$、$Q_xx$、$Q_yy$ | `undulator` | |
| `sinkz`、`sinkzp` | $\sin K_zs_j$、$\sin(K_zs_j + \varphi)$ | `undulator` | |
| `fx, fy, gx, gy` | $\partial_xA_x$、$\partial_yA_x$、$\partial_xA_y$、$\partial_yA_y$ (後に $\times\,dpds$、さらに逆行列の要素) | `undulator` | |
| `dpds` | $ds/p_n$ | `undulator` | `dpds` |
| `det1a` | $(1-\hat f_x)(1-\hat g_y) - \hat f_y\hat g_x$ | `undulator` | `det1a` |
| `fgx`、`fgy` | $\widehat{fg}_x$、$\widehat{fg}_y$ (B78: `dpds` が二重) | `undulator` | |
| `colb%Benv`、`Benv5` | 強ビームの 6×6、条件付き 5×5 の包絡 | `bbinit` | `bb['benv5']` |
| `colb%v_cen` | $\langle x_i\rangle = v_i z$ | `bbinit` | `vcen` |
| `colb%zslice` | スライスの重心・境界 ($\sigma_z$ 単位) | `bbinit` | `zc` |
| `colb%ce` | $Nr_e/\gamma$ | `bbinit` | `ce` |
| `colb%xangle(1:7)` | 交差角、$\sin$、$\cos$、$\tan$、傾き角、(6)(7) | `bbinit` | `xangle` |
| `colb%cod(1:5)` | 強ビームの軌道 | `bbinit` | `cod` |
| `sig11x`、`sig12x`、`sig11y`、`sig12y` | 衝突点でのスライスの包絡 | `beambeam` | |
| `asp` | $\sqrt{\sigma_{11y}/\sigma_{11x}}$ | `beambeam` | |
| `hi` | $h = -x^2/2\sigma_x^2 - y^2/2\sigma_y^2$ | `beambeam` | |
| `a2`、`a2x`、`a2y` | $e^h$、$e^h\,w_{2}$ | `beambeam` | |
| `w1x, w1y, w2x, w2y` | $\mathrm{Im}/\mathrm{Re}\,w(\zeta_{1,2})$ | `beambeam` | |
| `c1` | $2\sqrt\pi\,c/S$ (楕円)、$c$ (丸い) | `beambeam` | |
| `fx0`、`fy0` | キック $f_x, f_y$ | `beambeam` | |
| `gxy`、`dux`、`duy` | $g_{xy}$、$\partial U/\partial\sigma_x^2$、$\partial U/\partial\sigma_y^2$ 相当 | `beambeam` | |
| `sz` | $s^* = (z - z_s\sigma_z)/2$ | `beambeam` | |
| `rhoi`、`dsz`、`rnbs`、`ucE`、`delg` | $\rho_i$、スライスの半長、$N_\gamma$、$u_c$、放射する $\delta$ | `beambeam` | |
| `cnbs`、`cpbs`、`cubs` | $a_{\rm nrad}\gamma$、$\tfrac23r_e\gamma^3$、$c_{\rm uc}\gamma^2$ | `beambeam`、`bbinit` | `cnbs`、`cubs` |
| `SwI`、`xi` | 光子数分布の累積 (逆引き)、$\ln\xi$ | `besseltab` | `besseltab()` |
