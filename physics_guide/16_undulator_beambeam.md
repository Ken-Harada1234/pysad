# 16. アンジュレータ (UND) とビーム・ビーム (BEAMBEAM)

**この章で分かること**: SAD の `UND` と `BEAMBEAM` が追跡で行う写像を、ハミルトニアンと近似から組み立てる。どちらも**ヘルプに節がない** (ヘルプの要素の一覧に `UND`・`BEAMBEAM` の記述はない)。式は K64 のソース (`undulator.f`、`bbstrhl.f`) から起こし、途中式を本章で導く。(1) UND: ベクトルポテンシャル $A_x, A_y$ の模型と真空の Maxwell 方程式の条件、傍軸のハミルトニアン、スライスごとの陰的なシンプレクティック写像 (1 次)。(2) BEAMBEAM: 2 次元ガウス分布の電場 (ポテンシャルの Poisson 方程式から複素誤差関数の閉じた形まで)、スライスごとの 6 次元写像 (衝突点までの drift、キック、戻し)、強ビームの包絡とスライス、軌道・傾きの座標変換、**交差角の写像 (回転 + Lorentz boost) の導出**、ビームストラールング。2 次元ポテンシャルの数値積分と DOP853 を基準に確認する。

**範囲**: UND は追跡の写像 (K64 は磁場を設定できない: v2 だけ磁場つきの追跡と放射)。BEAMBEAM は弱-強モデルの追跡だけ (光学・Emittance[] には入らない)。交差角の写像 (16.12) は、座標軸の回転と Lorentz boost の合成から導き、$z$・$\delta$ の規格化まで含めて Fortran の式と一致することを確かめた (§8.2。この写像の幾何の読み — XANGLE が全交差角の半分に当たる — は式からの逆算で、ソースやヘルプに定義の記述はない)。光子リスト・ルミノシティ・制動放射 (BBBrem) は [docs/14_beambeam.md](../docs/14_beambeam.md) を参照。

**前提**: [13_radiation.md](13_radiation.md) (`anrad`・`cuc`・`cave`・`cuu`、光子の放射)、[03_symplectic_integrators.md](03_symplectic_integrators.md) (シンプレクティックな分割)、[02_drift.md](02_drift.md) (傍軸のドリフト)、[C_poisson_lie.md](C_poisson_lie.md) (ポアソン括弧とリー演算子)。

**読む順**: §1 (結論) → §3 (UND) → §6〜§8 (BEAMBEAM) → §12 (確認)。

## 1. 結論

- **UND の模型**: $A_x$、$A_y$ を (16.1) の形に与える。軸上で $B_y = F_0\cos(K_z s)$、$B_x = G_0\cos(K_z s + \varphi)$ ($F_0 = \mathtt{FBX}\,c/p$ [1/m]、$G_0 = \mathtt{FBY}\,c/p$、$K_z = 2\pi\,\mathtt{POLE}/L$)。$B = \nabla\times A$ は $\nabla\cdot B = 0$ を恒等的に満たし、$\nabla\times B = 0$ は $K_y^2 = K_z^2 + K_x^2$、$Q_x^2 = K_z^2 + Q_y^2$ のとき成り立つ (§3.1 で成分ごとに導く)。[und_01_field_map.out.txt](outputs/und_01_field_map.out.txt) (a): sympy の $\nabla\times\vec B$ を数値代入して最大 $8\times10^{-15}$ (B ~ 0.3)。成分の導出は [und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (a)。
- **UND の写像**: $H = [(p_x - A_x)^2 + (p_y - A_y)^2]/(2(1+\delta))$ (傍軸)。スライスごとに、2 型の母関数 $F_2 = q\cdot P + ds\,H(q, P; s_j)$ の陰的な写像 (記号は §3.3)。**1 次の精度**で厳密にシンプレクティック (任意の $ds$ で $\det J = 1$、$J^\top SJ = S$ が 40 桁で成り立つ: [und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (d))。確認 (b): 3 GeV、$B_y = 1$ T、$B_x = 0.4$ T、Kx・Qy・位相あり、10 周期・1 m で、DOP853 の積分との差は slice = 500・1000・2000・4000・8000 で、$|\Delta y|$ が $6.81\times10^{-10}$、$3.32\times10^{-10}$、$1.64\times10^{-10}$、$8.14\times10^{-11}$、$4.06\times10^{-11}$、$|\Delta p_y|$ が $1.30\times10^{-9}$、$6.35\times10^{-10}$、$3.14\times10^{-10}$、$1.56\times10^{-10}$、$7.76\times10^{-11}$。刻み数を 2 倍にしたときの比は $y$・$p_y$ とも 2.05、2.02、2.01、2.00 で、**1 次の収束**。$|\Delta x|$ の比も 2.24、2.13、2.06、2.03 で 2 に近づく。$|\Delta p_x|$ は $7.5\times10^{-12}$、$1.07\times10^{-11}$、$7.07\times10^{-12}$、$3.97\times10^{-12}$、$2.09\times10^{-12}$ で、500→1000 では増え (比 0.70)、その後 1.51、1.78、1.90 と 2 に近づく (単調でないので 1 次の証拠には $y$・$p_y$ を使う。$p_x$ の 1 次の係数が小さく、誤差が $10^{-11}$ の水準で高次の項や積分の許容誤差の影響を受けるためと読めるが、分けて確認していない)。
- **K64 の UND の癖**: 磁場のキーワードが設定できない (B37、gs で再現)。陰的な刻みの右辺に `dpds` を二重にかける (B78。**ソースの読み**。gs では磁場が設定できず確認できない。v2 は 1 回)。**$-dv\,L$ の項 (速度の差による $z$ の更新) がない** (疑い、§4、§13): 磁場なしの UND と DRIFT で $z$ が $-2.897\times10^{-11}$ m ($= dv\,L$、3 GeV・$\delta = 10^{-3}$・$L = 1$ m) 違う。
- **BEAMBEAM の力**: 強ビーム (2 次元ガウス、$\sigma_x,\sigma_y$) が弱粒子に与えるキックは、2 次元ガウス分布の電場の閉じた式 (複素誤差関数、(16.8))。ポテンシャル (16.6) が 2 次元の Poisson 方程式 $\nabla^2U = 4\pi\rho$ を満たすこと、遠方で $\ln r^2$ になること、閉じた形が (16.6) の積分から出ることを §6.1 で導く。[bb_01_kick_slice.out.txt](outputs/bb_01_kick_slice.out.txt) (a): 丸い・$\sigma_x > \sigma_y$・$\sigma_x < \sigma_y$ の 3 通り各 4 点で、$\Delta p_x$、$\Delta p_y$ が 2 次元ポテンシャルの数値積分の微分と $10^{-13}$ 以下の絶対差 (相対 $10^{-9}$) で一致。小振幅の傾き $-2c/(\sigma_x(\sigma_x+\sigma_y))$ ($c = N r_e/\gamma$) とも一致。
- **エネルギーのキック**: $\delta$ の変化 $= -g_{xy} + \tfrac14|\Delta p|^2$ ($g_{xy} = \sigma_{12x}\,\partial U/\partial\sigma_x^2 + \sigma_{12y}\,\partial U/\partial\sigma_y^2$) が、$\sigma_{12}\ne0$ (Twiss $\alpha\ne0$) の強ビームで $10^{-18}$ 以下の差で一致 (3 通り各 4 点)。
- **6 次元写像**: スライスごとに、衝突点まで drift ($x \mathrel{+}= p_x s^*$、$\delta \mathrel{-}=p_\perp^2/4$、$s^* = (z - z_s)/2$)、キック、同じ $s^*$ で戻す。$F = s^*(z)\,p_\perp^2/2$ の流れの共役で、厳密にシンプレクティック (流れの級数が 1 項・2 項で止まることを §6.3 で示す)。確認 (b): SLICE = 5・交差角 10 mrad・軌道のずれ・傾いた強ビームで、$\|J^\top SJ - S\|$ (規格化) $\le 1.5\times10^{-11}$。キックなし (NP = 0) の往復は恒等 ($\le10^{-19}$)。
- **交差角の写像 (16.12)**: 座標軸を $\varphi_c$ だけ回し、回した系の $x$ 方向に速度 $\beta = \sin\varphi_c$ で Lorentz boost し、運動量を新しい基準エネルギー $\cos\varphi_c$ で割る。運動量 ($\delta'$、$p_x'$、$p_y'$) と位置 ($x'$、$y'$、$z'$) の式がすべて、40 桁で一致する (差 $10^{-41}$) 形で再現される (§8.2、[bb_02_crossing_angle.out.txt](outputs/bb_02_crossing_angle.out.txt))。この入口の写像単独は $J^\top SJ = S/\cos\varphi_c$ (共形シンプレクティック。運動量の単位を $p_0\to p_0\cos\varphi_c$ に変えるため) で、出口の写像は入口の逆。小振幅の極限は $x' \approx x + \varphi_cz$、$z' \approx z$、$\delta' \approx \delta - \varphi_cp_x$、$p_x' \approx p_x$ ($p_\perp^2$ の補正だけ)。
- **ビームストラールング**: 1 スライスの光子数 $N_\gamma = \mathtt{cnbs}\,\rho_i\,dsz = a_{\rm nrad}\gamma\theta$ (`cnbs` $= a_{\rm nrad}\gamma$ が比 $1.000000000$)、光子エネルギー $u = \mathtt{cubs}\,\rho_i\,\xi$ (`cubs` $= c_{\rm uc}\gamma^2$、比 $0.99999997$)。$\xi$ は光子数のスペクトルの逆引き表。表の $\langle\xi\rangle = 0.3128$、$\langle\xi^2\rangle = 0.4180$ は理論値 $8/(15\sqrt3) = 0.3079$、$11/27 = 0.4074$ (§9 で Bessel 関数の積分から導き、数値でも確認) に対して $+1.6\%$、$+2.6\%$ 高い (§9、疑い (小))。

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

(`ff` = $A_x$、`gg` = $A_y$。$B\rho$ で規格化済み。) 磁場は $\vec B = \nabla\times\vec A$ ($A_s = 0$): $B_x = -\partial_sA_y$、$B_y = \partial_sA_x$、$B_z = \partial_xA_y - \partial_yA_x$。

**$F_0$ の項の磁場** ($K_y$ と $K_x$、$K_z$ の関係はまだ使わない)。$A_x^{(F)} = \dfrac{F_0}{K_z}\cos(K_xx)\cosh(K_yy)\sin(K_zs)$、$A_y^{(F)} = \dfrac{K_xF_0}{K_yK_z}\sinh(K_yy)\sin(K_xx)\sin(K_zs)$ について、1 項ずつ微分する:

$$
\begin{aligned}
B_x^{(F)} &= -\partial_sA_y^{(F)} = -\frac{K_x}{K_y}F_0\,\sin(K_xx)\sinh(K_yy)\cos(K_zs),\\
B_y^{(F)} &= \partial_sA_x^{(F)} = F_0\cos(K_xx)\cosh(K_yy)\cos(K_zs),\\
\partial_xA_y^{(F)} &= \frac{K_x^2F_0}{K_yK_z}\sinh(K_yy)\cos(K_xx)\sin(K_zs),\qquad
\partial_yA_x^{(F)} = \frac{K_yF_0}{K_z}\cos(K_xx)\sinh(K_yy)\sin(K_zs),\\
B_z^{(F)} &= \partial_xA_y^{(F)} - \partial_yA_x^{(F)} = -\frac{K_y^2 - K_x^2}{K_yK_z}F_0\cos(K_xx)\sinh(K_yy)\sin(K_zs)
\end{aligned}
\tag{16.2a}
$$

($B_z$ は $\frac{F_0}{K_yK_z}(K_x^2 - K_y^2)$ をくくり出した。sympy の出力と一致: [und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (a)。) 次の $\nabla\times\vec B$ の 3 成分を作る ($\vec B = \vec B^{(F)}$):

$$
\begin{aligned}
(\nabla\times\vec B)_x &= \partial_yB_z - \partial_sB_y
= -\frac{K_y^2 - K_x^2}{K_z}F_0\cos(K_xx)\cosh(K_yy)\sin(K_zs) + K_zF_0\cos(K_xx)\cosh(K_yy)\sin(K_zs)\\
&= \frac{F_0}{K_z}\left(K_z^2 + K_x^2 - K_y^2\right)\cos(K_xx)\cosh(K_yy)\sin(K_zs),\\
(\nabla\times\vec B)_y &= \partial_sB_x - \partial_xB_z = \frac{K_xK_z}{K_y}F_0\sin(K_xx)\sinh(K_yy)\sin(K_zs) - \frac{K_x(K_y^2 - K_x^2)}{K_yK_z}F_0\sin(K_xx)\sinh(K_yy)\sin(K_zs)\\
&= \frac{K_xF_0}{K_yK_z}\left(K_z^2 + K_x^2 - K_y^2\right)\sin(K_xx)\sinh(K_yy)\sin(K_zs),\\
(\nabla\times\vec B)_z &= \partial_xB_y - \partial_yB_x = -K_xF_0\sin(K_xx)\cosh(K_yy)\cos(K_zs) + K_xF_0\sin(K_xx)\cosh(K_yy)\cos(K_zs) = 0
\end{aligned}
\tag{16.2b}
$$

(sympy の出力と一致する: $x$ 成分 $\frac{F_0}{K_z}(K_x^2 - K_y^2 + K_z^2)\sin(K_zs)\cos(K_xx)\cosh(K_yy)$、$y$ 成分 $\frac{K_xF_0}{K_yK_z}(K_x^2 - K_y^2 + K_z^2)\sin(K_xx)\sin(K_zs)\sinh(K_yy)$。$z$ 成分は $\partial_xB_y = \partial_yB_x = -K_xF_0\sin(K_xx)\cosh(K_yy)\cos(K_zs)$ で打ち消し合う。) したがって $\nabla\times\vec B = 0$ の条件は

$$
K_y^2 = K_z^2 + K_x^2
\tag{16.2c}
$$

の 1 つだけで ($z$ 成分は恒等的に 0)、(16.2c) のもとで $B_z^{(F)} = -\dfrac{K_z}{K_y}F_0\cos(K_xx)\sinh(K_yy)\sin(K_zs)$ となり、

$$
\vec B_{F_0} = F_0\left(-\frac{K_x}{K_y}\sin(K_xx)\sinh(K_yy)\cos(K_zs),\ \cos(K_xx)\cosh(K_yy)\cos(K_zs),\ -\frac{K_z}{K_y}\cos(K_xx)\sinh(K_yy)\sin(K_zs)\right)
\tag{16.2}
$$

である。発散は条件なしで $\nabla\cdot\vec B = \partial_xB_x + \partial_yB_y + \partial_sB_z = \dfrac{F_0}{K_y}\left[-K_x^2 + K_y^2 - (K_y^2 - K_x^2)\right]\cos(K_xx)\sinh(K_yy)\cos(K_zs) = 0$ ($\vec B = \nabla\times\vec A$ なので恒等的に 0)。軸上 ($x = y = 0$) は $B_y = F_0\cos K_zs$、$B_x = B_z = 0$。$K_x = 0$ なら普通の平面アンジュレータ: $K_y = K_z$ で $B_y = F_0\cosh(K_zy)\cos(K_zs)$、$B_z = -F_0\sinh(K_zy)\sin(K_zs)$。

**$G_0$ の項の磁場**。$A_x^{(G)} = -\dfrac{Q_yG_0}{Q_xK_z}\sinh(Q_xx)\sin(Q_yy)\sin(K_zs+\varphi)$、$A_y^{(G)} = -\dfrac{G_0}{K_z}\cos(Q_yy)\cosh(Q_xx)\sin(K_zs+\varphi)$ について同じ手順で:

$$
\begin{aligned}
B_x^{(G)} &= -\partial_sA_y^{(G)} = G_0\cos(Q_yy)\cosh(Q_xx)\cos(K_zs+\varphi),\\
B_y^{(G)} &= \partial_sA_x^{(G)} = -\frac{Q_y}{Q_x}G_0\sinh(Q_xx)\sin(Q_yy)\cos(K_zs+\varphi),\\
B_z^{(G)} &= \partial_xA_y^{(G)} - \partial_yA_x^{(G)} = \left[-\frac{Q_x}{K_z} + \frac{Q_y^2}{Q_xK_z}\right]G_0\cos(Q_yy)\sinh(Q_xx)\sin(K_zs+\varphi) = -\frac{Q_x^2 - Q_y^2}{Q_xK_z}G_0\cos(Q_yy)\sinh(Q_xx)\sin(K_zs+\varphi)
\end{aligned}
\tag{16.2d}
$$

($\partial_xA_y^{(G)} = -\frac{G_0Q_x}{K_z}\cos(Q_yy)\sinh(Q_xx)\sin(\cdot)$、$\partial_yA_x^{(G)} = -\frac{Q_y^2G_0}{Q_xK_z}\sinh(Q_xx)\cos(Q_yy)\sin(\cdot)$。) $\nabla\times\vec B$ の $x$、$y$ 成分は、sympy の出力が $-\frac{G_0Q_y}{K_zQ_x}(K_z^2 - Q_x^2 + Q_y^2)\sin(Q_yy)\sin(K_zs+\varphi)\sinh(Q_xx)$ と $-\frac{G_0}{K_z}(K_z^2 - Q_x^2 + Q_y^2)\sin(K_zs+\varphi)\cos(Q_yy)\cosh(Q_xx)$ で、$z$ 成分は恒等的に 0。条件は

$$
Q_x^2 = K_z^2 + Q_y^2
\tag{16.2e}
$$

で、このとき $B_z^{(G)} = -\dfrac{K_z}{Q_x}G_0\cos(Q_yy)\sinh(Q_xx)\sin(K_zs+\varphi)$。軸上は $B_x = G_0\cos(K_zs+\varphi)$、$B_y = B_z = 0$。$G_0$ の項は、$F_0$ の項で $x\leftrightarrow y$、$K\leftrightarrow Q$ の役割を入れ替えた水平磁場の項である。

**Maxwell 方程式との関係**: 上の導出から、真空の Maxwell 方程式 ($\nabla\times\vec B = 0$、$\nabla\cdot\vec B = 0$) のうち、発散は $\vec B = \nabla\times\vec A$ から恒等的に、回転は (16.2c)(16.2e) のときに成り立つ。各成分が調和関数 ($\nabla^2 B = 0$) になる条件と同じである (例: (16.2) の $B_y$ の $\nabla^2$ は $(-K_x^2 + K_y^2 - K_z^2)B_y$)。**確認** (a): sympy で (16.1) から $\nabla\times\vec B$ を作り、$K_y, Q_x$ をこの関係にして、数値代入 10 点で最大 $8.2\times10^{-15}$ ($B \sim 0.3$)。条件を使わずに sympy で成分を出し、(16.2a)(16.2b)(16.2d) の因子 $(K_x^2 - K_y^2 + K_z^2)$、$(K_z^2 - Q_x^2 + Q_y^2)$ が現れ、条件のもとで 0 になることも確認した ([und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (a))。軸上の $\vec B = (G_0\cos(K_zs+\varphi),\ F_0\cos K_zs,\ 0)$。

**単位**: $F_0 = \mathtt{FBX}\,c/p$ ($\mathtt{FBX}$ [T]、$p$ [eV]、$c$ [m/s]) $= B_y/B\rho$ [1/m]。軌道の最大の角は軸上の $p_x = -A_x/(1+\delta)$ から $F_0/K_z$ ($= K/\gamma$、偏向パラメータ $K = eB\lambda_u/2\pi mc$)。**確認** (c): $F_0 = 0.3$ /m、$K_z = 62.83$ /m で $F_0/K_z = 4.7746\times10^{-3}$、追跡の最大の $p_x$ (20 スライスおきの標本) $4.7741\times10^{-3}$。

**使える範囲**: 周期の数 (POLE) が整数、$\varphi = 0$ のとき $A(0) = A(L) = 0$ だが、磁場は端で最大 ($\cos$)。つまり**端の処理 (終端の極) が模型にない**。1 周期の $\int B_y\,ds = F_0\int_0^{2\pi/K_z}\cos(K_zs)\,ds = F_0\,[\sin(K_zs)/K_z]_0^{2\pi/K_z} = 0$ なので、1 周期ごとに 1 次の積分は 0 (sympy でも確認、(b))。$\varphi \ne 0$ のときは、入口 ($s = 0$) の軸上の $A_y = -\dfrac{G_0\sin\varphi}{K_z} \ne 0$ で、出口 ($s = L$、整数周期) も同じ値なので、**端でポテンシャルの飛びが蹴りになる** (確認 (b): $A_y = -1.9\times10^{-4}$)。

### 3.2 ハミルトニアン

$s$ を独立変数にした SAD のハミルトニアン ([01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)) は、$A_s = 0$ の直線で $H = -\sqrt{(1+\delta)^2 - (p_x-A_x)^2 - (p_y-A_y)^2}$ (定数の $E/(p_0v_0)$ の項は別)。運動学的な横運動量 $\pi_{x,y} = p_{x,y} - A_{x,y}$、$\pi^2 = \pi_x^2 + \pi_y^2$ で、$p_n = 1+\delta$ として

$$
H = -p_n\sqrt{1 - \frac{\pi^2}{p_n^2}} = -p_n\left(1 - \frac{\pi^2}{2p_n^2} - \frac{\pi^4}{8p_n^4} - \cdots\right) = -p_n + \frac{\pi^2}{2p_n} + \frac{\pi^4}{8p_n^3} + O\!\left(\frac{\pi^6}{p_n^5}\right)
\tag{16.3a}
$$

($\sqrt{1-u} = 1 - u/2 - u^2/8 - \cdots$、$u = \pi^2/p_n^2$。sympy の級数展開と一致: [und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (c)。) 定数 $-p_n$ と $\pi^4$ 以上を落として

$$
H_{\rm UND} = \frac{(p_x - A_x)^2 + (p_y - A_y)^2}{2(1+\delta)}
\tag{16.3}
$$

(傍軸近似。落とした $\pi^4/(8p_n^3)$ は主項 $\pi^2/(2p_n)$ に対して $\pi^2/(4p_n^2)$ の比。最大の角 $\pi \approx F_0/K_z = 4.8\times10^{-3}$ ($K/\gamma$、§3.1 (c) の例) で $6\times10^{-6}$、確認 (b) の例 ($F_0/K_z\approx 1.6\times10^{-3}$) で $6\times10^{-7}$。落とした定数 $-p_n$ の $\delta$ 微分と $E/(p_0v_0)$ の項の組が、$z$ の更新の $-dv$ ([01_variables_and_conventions.md](01_variables_and_conventions.md) (1.9)) に当たる。UND はそれを入れていない: §4。)

運動方程式は正準方程式から。$H_{\rm UND}$ の偏微分を 1 つずつ書くと ($\pi_x = p_x - A_x$):

$$
x' = \frac{\partial H}{\partial p_x} = \frac{\pi_x}{p_n},\quad
y' = \frac{\pi_y}{p_n},\quad
p_x' = -\frac{\partial H}{\partial x} = \frac{\pi_x\,\partial_xA_x + \pi_y\,\partial_xA_y}{p_n},\quad
p_y' = \frac{\pi_x\,\partial_yA_x + \pi_y\,\partial_yA_y}{p_n},
$$
$$
z' = \frac{\partial H}{\partial\delta} = -\frac{\pi_x^2 + \pi_y^2}{2p_n^2},\qquad \delta' = -\frac{\partial H}{\partial z} = 0.
\tag{16.3b}
$$

($-\partial_x$ に $\pi_x = p_x - A_x$ の微分 $-\partial_xA_x$ が掛かり、符号が反転して $+$ になる。$H$ は $z$ に依らないので $\delta$ は保存。SAD の $(z, \delta)$ は $(q, p)$ の対で、$z' = \partial H/\partial\delta$。sympy の出力と一致、(c)。)

### 3.3 スライスの写像

1 スライス ($ds = L/\mathtt{SLICE}$、$j$ 番目、$s_j = j\,ds$) で、K64 は次の写像を使う ($\vec q = (x,y)$、$\vec p = (p_x, p_y)$):

$$
\vec P = \vec p - ds\,\partial_{\vec q}H(\vec q, \vec P;\ s_j),\qquad
\vec q' = \vec q + ds\,\partial_{\vec P}H(\vec q, \vec P;\ s_j),\qquad
z' = z + ds\,\partial_\delta H(\vec q, \vec P;\ s_j),\qquad \vec p' = \vec P
\tag{16.4}
$$

($H$ の中の $A$ は古い位置 $\vec q$ と終点の $s_j$ で評価し、運動量は新しい $\vec P$。) これは 2 型の母関数 $F_2(\vec q, z, \vec P, \delta') = \vec q\cdot\vec P + z\,\delta' + ds\,H(\vec q, \vec P, \delta'; s_j)$ から、2 型の変換の規則 $p_i = \partial F_2/\partial q_i$、$Q_i = \partial F_2/\partial P_i$ (新しい座標 $Q$ と運動量 $P$) で出る:

$$
\vec p = \frac{\partial F_2}{\partial\vec q} = \vec P + ds\,\partial_{\vec q}H,\quad
\vec q' = \frac{\partial F_2}{\partial\vec P} = \vec q + ds\,\partial_{\vec P}H,\quad
\delta = \frac{\partial F_2}{\partial z} = \delta',\quad
z' = \frac{\partial F_2}{\partial\delta'} = z + ds\,\partial_\delta H.
$$

(1 つめの式が (16.4) の $\vec P$ の式、3 つめは $\delta$ が保たれること。) 任意の $ds$ について**正準変換** (母関数から出る写像はシンプレクティック)。精度は 1 次 (シンプレクティック Euler)。数値の確認は、有限の $ds = 0.02$ で 1 スライスの $6\times6$ ヤコビアンを sympy で 40 桁で評価し、$\max|J^\top SJ - S| = 2.3\times10^{-41}$、$\det J = 1$ ([und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (d))。

$H$ が $\vec P$ の 2 次なので陰的な式は線形で、2×2 の行列で解ける。$\partial_xH = -\frac1{p_n}[(P_x - A_x)\partial_xA_x + (P_y - A_y)\partial_xA_y]$ なので、(16.4) の $\vec P$ の式を成分で書くと

$$
P_x = p_x + \frac{ds}{p_n}\left[(P_x - A_x)\partial_xA_x + (P_y - A_y)\partial_xA_y\right],\qquad
P_y = p_y + \frac{ds}{p_n}\left[(P_x - A_x)\partial_yA_x + (P_y - A_y)\partial_yA_y\right].
$$

$\hat f_x = \dfrac{ds}{p_n}\partial_xA_x$、$\hat g_x = \dfrac{ds}{p_n}\partial_xA_y$、$\hat f_y = \dfrac{ds}{p_n}\partial_yA_x$、$\hat g_y = \dfrac{ds}{p_n}\partial_yA_y$ とおいて $P$ の項を左辺に集めると ($A_x, A_y$ の項は右辺):

$$
\begin{pmatrix}1 - \hat f_x & -\hat g_x\\ -\hat f_y & 1 - \hat g_y\end{pmatrix}\begin{pmatrix}P_x\\P_y\end{pmatrix} = \begin{pmatrix}p_x - \widehat{fg}_x\\ p_y - \widehat{fg}_y\end{pmatrix},\quad
\widehat{fg}_x = \hat f_xA_x + \hat g_xA_y,\ \widehat{fg}_y = \hat f_yA_x + \hat g_yA_y
\tag{16.5}
$$

行列式 $\det_{1a} = (1 - \hat f_x)(1 - \hat g_y) - \hat f_y\hat g_x$ で、逆行列は余因子から $\dfrac{1}{\det_{1a}}\begin{pmatrix}1 - \hat g_y & \hat g_x\\ \hat f_y & 1 - \hat f_x\end{pmatrix}$。したがって

$$
P_x = \frac{(1-\hat g_y)\,r_x + \hat g_x\,r_y}{\det_{1a}},\quad
P_y = \frac{\hat f_y\,r_x + (1-\hat f_x)\,r_y}{\det_{1a}},\qquad r_x = p_x - \widehat{fg}_x,\ r_y = p_y - \widehat{fg}_y.
$$

`undulator.f` はこの逆行列の要素を `fx = (1-gy)/det1a`、`gx = gx/det1a`、`fy = fy/det1a`、`gy = (1-fx)/det1a` として持ち (`undulator.f:84-103`)、`px = fx*tmp + gx*(py - fgy)`、`py = fy*tmp + gy*(py - fgy)` と掛ける (`tmp = px - fgx`)。上の式と一致する。sympy で $\vec P$ を上の式で作り、(16.4) の陰的な式に代入した残差は $10^{-128}$ 以下 ((d))。位置と $z$ の更新は

$$
x' = x + dpds\,(P_x - A_x),\quad y' = y + dpds\,(P_y - A_y),\quad
z' = z - \frac{ds}{2p_n^2}\left[(P_x - A_x)^2 + (P_y - A_y)^2\right]
\tag{16.5a}
$$

(`undulator.f:105-109`。$z'$ は `z = z - 0.5*(…)*dpds**2/ds` で、$dpds^2/ds = ds/p_n^2$ なので (16.4) の $z' = z + ds\,\partial_\delta H$、$\partial_\delta H = -\pi^2/(2p_n^2)$ と一致する。) 更新に使う $A$ は古い位置 $\vec q$ のもの。

**入口・出口**: 外側 ($A = 0$) の運動量は正準運動量に等しい。入口では外の $p_x p_n$ を正準運動量として持ち込み、内側の力学的運動量は $p - A(0)$ になる。出口では最後の正準運動量 $/p_n$ を外の運動量として返す ($A(L) = 0$ でないとその分が蹴りになる)。**確認** (b): 入口の $A_y = -1.88\times10^{-4}$、出口も同じ値 (位相 $\varphi = 0.3$) のとき、DOP853 の積分も同じ約束 (入口で正準運動量 = 外の運動量) で行って、上の誤差になる。

**K64 の疑い (B78、ソースの読み)**: `fgx = dpds*(ff*fx + gg*gx)` と、すでに `dpds` を掛けた `fx, gx` にもう 1 回 `dpds` を掛ける (`undulator.f:92-93`。ソースで確認した)。(16.5) の右辺の $\widehat{fg}$ が $ds/p_n$ 倍に弱まる。$\widehat{fg}$ は力 $-\frac1{p_n}(A_x\partial_xA_x + A_y\partial_xA_y)$ (ワイグラー運動の収束の項、$A\,\partial A$) に当たるので、**収束の項が刻み幅に比例して小さくなる** ($ds\to0$ で消える)。v2 は 1 回だけ掛ける (本章の確認は v2 の式)。**根拠の程度**: ソースの読みと、刻み数を変えた数値実験 (docs/30)。gs では UND の磁場が設定できず (B37)、K64 で再現して確認できていない。K64 では UND に磁場が無いので、現状では結果に効かない。

### 3.4 刻み数

1 次の写像なので、誤差は刻み数に反比例 (§1 の確認: $y$、$p_y$ の比 2.0)。確認 (b) の例 (周期 0.1 m、$K\sim9.4$〜) で、$p_x$ の誤差が $2\times10^{-12}$ (8000 スライス = 800 スライス/周期) に下がる。K64 の既定の目安は書かれていない ([docs/19_other_elements.md](../docs/19_other_elements.md) の例は 2000 スライス/20 極)。

## 4. UND の確認

(確認 (a)〜(d) は §1・§3.1・§3.3 の数字。) (d): 磁場なしの UND と DRIFT ($L = 1$ m、3 GeV、$x_0 = 10^{-4}$、$p_{x} = 2\times10^{-4}$、$\delta = 10^{-3}$): $x$ は $2.998001998\times10^{-4}$ (UND、傍軸 $x' = p_x/p_n$) と $2.998002038\times10^{-4}$ (DRIFT、$x' = p_x/p_z$)、差 $4\times10^{-12}$ (傍軸の落とした項: $p_x/p_z - p_x/p_n \approx p_x\,p_x^2/(2p_n^3)\cdot L\approx 4\times10^{-12}$)。$z$ は $-1.99600599\times10^{-8}$ と $-1.99310907\times10^{-8}$、差 $-2.897\times10^{-11}$ m。これは $dv\,L$ ($dv = -\delta(2+\delta)/[h_1(h_1 + p_nh_0)]$、$h_1 = \sqrt{1 + (p_np_0)^2}$、$L = 1$) の値と一致する: UND の $z$ の更新に $-dv\,L$ の項がない (DRIFT は持つ、[02_drift.md](02_drift.md))。

**$R_{56}$ との関係**: $\delta$ が小さいとき $h_1\approx h_0 = \gamma$、$p_n\approx1$ なので $dv\approx-\delta\cdot2/(\gamma\cdot2\gamma) = -\delta/\gamma^2$。$\delta = 10^{-3}$、$\gamma = 5871$ で $-2.90\times10^{-11}$ (厳密な式との比は $0.9985$、$\delta$ の有限の分。0.5 GeV でも同じ比、[und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) (e))。したがって $z$ の更新の欠落は $R_{56} = L/\gamma^2$ の欠落に当たり、$\gamma = 5871$ で $\sim3\times10^{-8}\,L$。高エネルギーでは小さいが、低エネルギーでは大きな誤差になる (疑い、§13)。

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

($c$ は `ce/nslice`、$\mathtt{ce} = N r_e/\gamma$。)

**(16.6) の導出**。入力は、点状の (横の大きさのない) 薄いスライスの遠方のキック $\Delta\vec p_\perp = -\dfrac{2Nr_e}{\gamma}\dfrac{\vec r}{r^2}$ (2 次元のクーロンの法則。弱粒子の電荷が強ビームと反対で、収束) である。ガウス分布の 2 次元電場は、$\Delta\vec p = -c\,\nabla U$ と書いたとき、$U$ が次の 2 つで決まる: (i) 2 次元の Poisson 方程式 $\nabla^2U = 4\pi\rho$ ($\rho = \dfrac{1}{2\pi\sigma_x\sigma_y}e^{-x^2/2\sigma_x^2 - y^2/2\sigma_y^2}$ は規格化したガウス分布。点電荷なら $\rho = \delta_D(\vec r)$)、(ii) 遠方で $U\to\ln r^2$ (点電荷の $-c\nabla\ln r^2 = -2c\vec r/r^2$ に一致)。(16.6) の $U$ がこれを満たすことを示す。被積分関数の $x, y$ に依る部分を $g = -e^{-x^2/a - y^2/b}/\sqrt{ab}$ ($a = 2\sigma_x^2 + t$、$b = 2\sigma_y^2 + t$) とおくと、$E \equiv e^{-x^2/a - y^2/b}/\sqrt{ab}$ の $t$ 微分は $\partial_ta = \partial_tb = 1$ から

$$
\partial_tE = E\left[\frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac1{2a} - \frac1{2b}\right] = \frac14E\left[\left(\frac{4x^2}{a^2} - \frac2a\right) + \left(\frac{4y^2}{b^2} - \frac2b\right)\right],
$$

一方 $\nabla^2E$ ($x$ と $y$ の 2 階微分) $= E\left[\left(\dfrac{4x^2}{a^2} - \dfrac2a\right) + \left(\dfrac{4y^2}{b^2} - \dfrac2b\right)\right]$ なので、$\nabla^2g = -\nabla^2E = -4\,\partial_tE$ (sympy で確認、[bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (a))。(16.6) の $t$ 積分は $x, y$ に依らない部分 (定数に当たる "1" の項) の微分が 0 なので $\nabla^2U = \int_0^\infty\nabla^2g\,dt = -4[E]_0^\infty = 4E|_{t=0} = \dfrac{4}{2\sigma_x\sigma_y}e^{-x^2/2\sigma_x^2 - y^2/2\sigma_y^2} = 4\pi\rho$ ($E|_{t=\infty} = 0$、$a_0b_0 = 4\sigma_x^2\sigma_y^2$)。したがって (i) を満たす。遠方の挙動 (ii) は次の丸いビームで見る。

**丸いビームの検算** ($\sigma_x = \sigma_y = \sigma$、$u = 1/(2\sigma^2+t)$ と置換、$du = -u^2dt$): $\dfrac{\partial U}{\partial x} = \displaystyle\int_0^\infty\frac{2x\,e^{-r^2/(2\sigma^2+t)}}{(2\sigma^2+t)^2}dt = 2x\int_0^{1/2\sigma^2}e^{-r^2u}\,du = \frac{2x}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)$。したがって

$$
\Delta\vec p_\perp = -\frac{2Nr_e}{\gamma n_{\rm slice}}\,\frac{\vec r}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)\qquad(\sigma_x = \sigma_y)
\tag{16.7}
$$

(コードの丸いビーム: `a2 = 2*(1-exp(hi))/(x**2+y**2)`、`fx0 = c1*a2*x`。`hi` $= -x^2/2\sigma_x^2 - y^2/2\sigma_y^2$。) $r\gg\sigma$ で $-2c\vec r/r^2$ (点電荷、(ii) の係数) に、小振幅 $r \ll \sigma$ で $1 - e^{-r^2/2\sigma^2}\approx r^2/2\sigma^2$ から $-\dfrac{N r_e}{\gamma\sigma^2}\vec r$ になる ($\sigma_x\ne\sigma_y$ でも $\Delta p_x = -2c\,x/(\sigma_x(\sigma_x+\sigma_y))$、ビーム・ビームパラメータ $\xi_x = N r_e\beta_x^*/(2\pi\gamma\sigma_x(\sigma_x+\sigma_y))$ の元)。$U$ 自体は丸い場合 $w = r^2/2\sigma^2$ ($q = r^2/(2\sigma^2+t)$ と置換: $dt = -r^2dq/q^2$、$1/(2\sigma^2+t) = q/r^2$) で

$$
U_{\rm round} = \int_0^w\frac{1 - e^{-q}}{q}\,dq \equiv \mathrm{Ein}(w)\ \xrightarrow{w\to\infty}\ \ln w + \gamma_E = \ln r^2 - \ln2\sigma^2 + \gamma_E
\tag{16.7a}
$$

($\gamma_E$ はオイラー定数。$\mathrm{Ein}$ の数値積分との一致: [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (b)。) 遠方で $\ln r^2$ になる (ii)。$\ln r^2$ の係数が 1 で、$\sigma$ に依る定数 $-\ln2\sigma^2 + \gamma_E$ が付く。これが §6.2 の $2\ln(\sigma_x+\sigma_y)$ の項の理由になる。

**一般の楕円 (閉じた形)**: (16.6) の閉じた形は、複素誤差関数 $w(\zeta) = e^{-\zeta^2}\mathrm{erfc}(-i\zeta)$ で書かれる。$\sigma_x > \sigma_y$ のとき、$S = \sqrt{2(\sigma_x^2 - \sigma_y^2)}$、$\zeta_1 = (|x| + i|y|)/S$、$\zeta_2 = (|x|\sigma_y/\sigma_x + i|y|\sigma_x/\sigma_y)/S$ として

$$
f_y + i f_x = \frac{2\sqrt\pi\,c}{S}\left[w(\zeta_1) - e^{h}\,w(\zeta_2)\right],\qquad h = -\frac{x^2}{2\sigma_x^2} - \frac{y^2}{2\sigma_y^2},\quad \Delta p_x = -\mathrm{sgn}(x)\,f_x,\ \Delta p_y = -\mathrm{sgn}(y)\,f_y
\tag{16.8}
$$

(コード: `cerrf(d1x,d1y,w1y,w1x)` で $w_{1y} = \mathrm{Re}\,w$、$w_{1x} = \mathrm{Im}\,w$ を得て、`fx0 = c1*(w1x - a2x)`、`fy0 = c1*(w1y - a2y)`、`c1 = ce*2*sqrpi/S/nslice`。象限の符号は $|x|, |y|$ を使うので、最後に符号を戻す。) $\sigma_x < \sigma_y$ は $x \leftrightarrow y$ を入れ替えた式 (コードの 3 つめの枝)。

**(16.8) の導出** ($x, y > 0$、$c = 1$ としたときの $\partial U/\partial y + i\,\partial U/\partial x$ について。2 つの項の差が $a - b = 2(\sigma_x^2 - \sigma_y^2) = S^2$ (一定) であることを使う)。

1. (16.6) を微分する: $\partial_xU = 2x\displaystyle\int_0^\infty\frac{e^{-x^2/a - y^2/b}}{a^{3/2}b^{1/2}}dt$、$\partial_yU = 2y\displaystyle\int_0^\infty\frac{e^{-x^2/a - y^2/b}}{a^{1/2}b^{3/2}}dt$ ($\partial_x$ は指数の $-x^2/a$ から $-2x/a$ が出て、被積分関数の $1/\sqrt{ab}$ に掛かる。符号は $1 - e^{\ldots}$ の $-$ で $+$ になる)。
2. 積分変数を $\rho = \sqrt{b/a}$ に変える。$a - b = S^2$ から $a = S^2/(1-\rho^2)$、$b = \rho^2S^2/(1-\rho^2)$。$t = 0$ で $\rho = \rho_0 = \sigma_y/\sigma_x$、$t = \infty$ で $\rho = 1$。$dt = da = 2\rho S^2\,d\rho/(1-\rho^2)^2$、$a^{3/2}b^{1/2} = a^2\rho = \rho S^4/(1-\rho^2)^2$、$a^{1/2}b^{3/2} = a^2\rho^3$。指数は $-\dfrac{x^2}{a} - \dfrac{y^2}{b} = -\dfrac{(1-\rho^2)(x^2 + y^2/\rho^2)}{S^2}$。よって
   $$\partial_xU = \frac{4x}{S^2}\int_{\rho_0}^1e^{E(\rho)}d\rho,\qquad \partial_yU = \frac{4y}{S^2}\int_{\rho_0}^1\frac{e^{E(\rho)}}{\rho^2}d\rho,\qquad E(\rho) = -\frac{(1-\rho^2)(x^2 + y^2/\rho^2)}{S^2}.$$
   ($t$ の積分との数値の一致: [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (c)、相対 $10^{-12}$ 以下。)
3. 組み合わせ $-\partial_yU + i\,\partial_xU = \dfrac{4}{S^2}\displaystyle\int_{\rho_0}^1\left(-\frac{y}{\rho^2} + ix\right)e^{E(\rho)}d\rho$ を作る。複素変数 $w(\rho) = y/\rho + i\rho x$ の導関数が $dw/d\rho = -y/\rho^2 + ix$ で、括弧と一致する。指数は恒等式 $x^2 - y^2 + y^2/\rho^2 - \rho^2x^2 = (x-iy)^2 + w(\rho)^2$ (右辺の虚部 $-2xy + 2xy$ が打ち消し合う) により $E(\rho) = -[(x-iy)^2 + w(\rho)^2]/S^2$。したがって
   $$-\partial_yU + i\,\partial_xU = \frac{4}{S^2}e^{-(x-iy)^2/S^2}\int_{w(\rho_0)}^{w(1)}e^{-w^2/S^2}dw.$$
4. 積分は誤差関数: $\int e^{-w^2/S^2}dw = \dfrac{\sqrt\pi S}{2}\mathrm{erf}(w/S)$。端点は $w(1) = y + ix = i(x - iy)$、$w(\rho_0) = y\sigma_x/\sigma_y + i\sigma_yx/\sigma_x$。よって $-\partial_yU + i\,\partial_xU = \dfrac{2\sqrt\pi}{S}e^{-(x-iy)^2/S^2}\left[\mathrm{erf}\!\left(\dfrac{i(x-iy)}{S}\right) - \mathrm{erf}\!\left(\dfrac{w(\rho_0)}{S}\right)\right]$。
5. $w$ 関数に書き直す。$e^{u^2}\mathrm{erf}(u) = e^{u^2} - e^{u^2}\mathrm{erfc}(u) = e^{u^2} - w(iu)$ ($w(iu) = e^{u^2}\mathrm{erfc}(u)$)。1 つ目は $u = i(x-iy)/S$、$u^2 = -(x-iy)^2/S^2$ で $e^{-(x-iy)^2/S^2} - w(i u)$、$iu = -(x-iy)/S = -\bar\zeta_1$ ($\zeta_1 = (x+iy)/S$) なので $w(-\bar\zeta_1) = \overline{w(\zeta_1)}$ ($w(-\bar\zeta) = \overline{w(\zeta)}$ の性質)。2 つ目は $u_0 = w(\rho_0)/S$ として $e^{-(x-iy)^2/S^2}\mathrm{erf}(u_0) = e^{-(x-iy)^2/S^2} - e^{-(x-iy)^2/S^2 - u_0^2}\,w(iu_0)$ ($e^{-(x-iy)^2/S^2}\mathrm{erfc}(u_0) = e^{-(x-iy)^2/S^2 - u_0^2}\cdot e^{u_0^2}\mathrm{erfc}(u_0)$ と見る)。ここで $-(x-iy)^2/S^2 - u_0^2 = E(\rho_0) = h$ ($\rho_0$ での指数は $-x^2/2\sigma_x^2 - y^2/2\sigma_y^2$)、$iu_0 = -\bar\zeta_2$ ($\zeta_2 = (x\sigma_y/\sigma_x + iy\sigma_x/\sigma_y)/S$) なので $w(iu_0) = \overline{w(\zeta_2)}$。引き算 (1 つ目 − 2 つ目) で $e^{-(x-iy)^2/S^2}$ の項が打ち消し合い、
   $$-\partial_yU + i\,\partial_xU = \frac{2\sqrt\pi}{S}\left[e^{h}\,\overline{w(\zeta_2)} - \overline{w(\zeta_1)}\right].$$
6. 複素共役をとると $-\partial_yU - i\,\partial_xU = \dfrac{2\sqrt\pi}{S}[e^hw(\zeta_2) - w(\zeta_1)]$、すなわち $\partial_yU + i\,\partial_xU = \dfrac{2\sqrt\pi}{S}[w(\zeta_1) - e^hw(\zeta_2)]$。$f = c\,\partial U$ なので (16.8) の $f_y + if_x$ に一致する。$x, y$ の符号は、$U$ が $x, y$ の偶関数であることから $|x|, |y|$ で計算して $\partial_xU$、$\partial_yU$ にそれぞれ $\mathrm{sgn}(x)$、$\mathrm{sgn}(y)$ を掛ければよい。

(scipy の複素誤差関数 (`wofz`) で上の最終式を計算し、$t$ 積分の $\partial_yU + i\partial_xU$ と 12 桁で一致: (c)。) **確認** (a): 丸い・$\sigma_x > \sigma_y$ ($0.2\ {\rm mm}$、$0.05\ {\rm mm}$)・$\sigma_x < \sigma_y$ (逆) の 3 通り、各 4 点 ($x$ が 0.2〜3 $\sigma$) で、(16.6) を scipy の 2 重の数値積分で微分した結果と、絶対差 $\le 1.4\times10^{-13}$ (キックの大きさ $\sim5\times10^{-5}$、相対 $3\times10^{-9}$) で一致。

### 6.2 $z$ に依存する包絡と $\delta$ のキック

強ビームのスライスは、衝突点の位置 $s^*$ により、衝突時の横の大きさが違う (包絡が $s^*$ の 2 次関数で、$\sigma_x^2(s) = \sigma_{11x} + 2\sigma_{12x}s + \sigma_{22x}s^2$ (`sig11x = Benv5(1) + 2*Benv5(2)*sz + Benv5(7)*sz**2`)、$\sigma_{12x}(s) = \sigma_{12x} + \sigma_{22x}s$ (`sig12x = Benv5(2) + Benv5(7)*sz`、$\partial_s\sigma_x^2 = 2\sigma_{12x}(s)$)。したがってポテンシャルは $z$ にも依存し、運動量 $\delta$ が変わる。物理的なポテンシャルは、遠方で $\sigma$ に依らず $\ln r^2$ + 定数になるもの $U_{\rm phys} = U + 2\ln(\sigma_x + \sigma_y)$ である。理由: 丸い場合 (16.7a) で $U\to\ln r^2 - \ln2\sigma^2 + \gamma_E$ なので $\sigma$ に依る定数 $-\ln2\sigma^2$ が付き、$2\ln(2\sigma) = \ln4\sigma^2$ を足すと $\ln r^2 + \ln2 + \gamma_E$ で $\sigma$ に依らない。楕円でも、遠方で $U + 2\ln(\sigma_x+\sigma_y) - \ln r^2\to\ln2 + \gamma_E$ が $\sigma_x, \sigma_y$ に依らず、方向にも依らない (数値で確認: $\sigma = (1,1),(2,1),(3,0.5),(1,4)$ で $r = 300$、誤差 $\sim10^{-4}$ (残りは $O(\sigma^2/r^2)$)、[bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (b2))。遠方の粒子は、電荷の総量だけで決まる力を受けるはずで、ビームの大きさの $s^*$ 依存でエネルギーを変えてはいけない、というのが必要な理由である。

**$\delta$ のキックは $\delta \mathrel{-}= g_{xy}$**、

$$
g_{xy} = c\left(\sigma_{12x}\,\frac{\partial U_{\rm phys}}{\partial\sigma_x^2} + \sigma_{12y}\,\frac{\partial U_{\rm phys}}{\partial\sigma_y^2}\right)\quad\left(\text{丸いビーム: }\frac{\partial U_{\rm phys}}{\partial\sigma^2} = \frac{e^{h}}{\sigma^2},\ g_{xy} = c\,\frac{\sigma_{12x}}{\sigma_{11x}}e^{h}\right)
\tag{16.9}
$$

(コード: `gxy = c1*sig12x/sig11x*exp(hi)`、楕円は `dux, duy` を使い `gxy = sig12x*dux + sig12y*duy`。)

**導出**: キックを、位置・$s^*$ に依存する生成関数 $W = c\,U_{\rm phys}(x, y;\ \sigma(s^*))$ による変換 ($\vec p\to\vec p - \partial_{\vec q}W$、$\delta\to\delta - \partial_zW$) と見る ($\delta' = -\partial H/\partial z$ の積分)。衝突点 $s^* = (z - z_s)/2$ なので $\partial s^*/\partial z = \tfrac12$。$U_{\rm phys}$ が $z$ に依るのは $\sigma_x^2(s^*)$、$\sigma_y^2(s^*)$ を通してだけなので、連鎖律で

$$
\partial_zW = c\cdot\frac12\left[\frac{\partial\sigma_x^2}{\partial s^*}\frac{\partial U_{\rm phys}}{\partial\sigma_x^2} + \frac{\partial\sigma_y^2}{\partial s^*}\frac{\partial U_{\rm phys}}{\partial\sigma_y^2}\right] = c\left[\sigma_{12x}(s^*)\frac{\partial U_{\rm phys}}{\partial\sigma_x^2} + \sigma_{12y}(s^*)\frac{\partial U_{\rm phys}}{\partial\sigma_y^2}\right] = g_{xy}
$$

($\partial_s\sigma_x^2 = 2\sigma_{12x}(s)$ と $\frac12$ が約分)。丸い場合の微分は (16.7a) の $U_{\rm phys} = \mathrm{Ein}(w) + \ln4\sigma^2$ ($w = r^2/2\sigma^2$) から、$\mathrm{Ein}'(w) = (1 - e^{-w})/w$、$\partial w/\partial\sigma^2 = -w/\sigma^2$、$\partial_{\sigma^2}\ln4\sigma^2 = 1/\sigma^2$ を使って $\dfrac{1-e^{-w}}{w}\cdot\left(-\dfrac{w}{\sigma^2}\right) + \dfrac1{\sigma^2} = \dfrac{e^{-w}}{\sigma^2}$ (sympy で確認、(b))。この $e^{-w} = e^h$。**確認** (a): 3 通りの強ビーム ($\alpha_x = 0.3, 0.5, -0.2$ など、$\sigma_{12} = -\alpha\epsilon$) で、1 スライス・$z = 0$ ($s^* = 0$) の $\Delta\delta$ が $-g_{xy} + \tfrac14|\Delta p_\perp|^2$ ($g_{xy}$ は $U_{\rm phys}$ の $\sigma^2$ 微分を数値積分で作ったもの、$\tfrac14|\Delta p|^2$ は §6.3 の戻しの項) と絶対差 $\le5\times10^{-18}$ (値 $\sim10^{-9}$) で一致する。$2\ln(\sigma_x+\sigma_y)$ の項を入れないと一致しない (丸い場合に符号さえ違う)。

### 6.3 衝突点までの drift と戻し

粒子の位置 $z$ と、スライスの位置 $z_s$ ($\sigma_z$ 単位の重心 $\times\sigma_z$) の相手が、2 つのビームが対向して動くので、衝突点は $s^* = (z - z_s)/2$ (相対速度 $2c$)。粒子は衝突点まで drift してキックを受け、同じ距離だけ戻る:

$$
M_{\rm slice} = D^{-1}\circ K\circ D,\qquad D:\ x \mathrel{+}= p_x s^*,\ y \mathrel{+}= p_ys^*,\ \delta \mathrel{-}= \tfrac14(p_x^2 + p_y^2)
\tag{16.10}
$$

($D$ の中の $s^* = (z - z_s)/2$ は D の前の $z$。$D^{-1}$ は、同じ $F$ の流れの**時間 $-1$** を、キック後の $p$ で行う: $x \mathrel{-}= p_xs^*$、$\delta \mathrel{+}= \tfrac14p_\perp^2$ (同じ $s^*$。$K$ は $z$ を変えないので $s^*$ は変わらない)。$D$ の単純な逆写像 (古い $p$ を使う) ではなく、$p$ が $K$ で変わった分だけ違う。$K$: $p_{x,y} \mathrel{-}= f_{x,y}$、$\delta \mathrel{-}= g_{xy}$。)

**$D$ の導出**: 生成関数 $F(x, y, p_x, p_y, z) = s^*(z)\,H_0$、$H_0 = (p_x^2 + p_y^2)/2$ ($s^* = (z - z_s)/2$) の流れを時間 1 だけ進める。[C_poisson_lie.md](C_poisson_lie.md) の約束 (ポアソン括弧 $[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg)$、転送は $\exp(-:F:)$ で、$g\mapsto g + [g, F] + \frac12[[g, F], F] + \cdots$) で書くと、座標・運動量 $g$ について級数の各項は ($(z, \delta)$ を共役な対として):

| $g$ | $[g, F]$ | $[[g, F], F]$ |
|---|---|---|
| $x$ | $\partial F/\partial p_x = s^*p_x$ | $[s^*p_x, F] = s^*[p_x, F] + p_x[s^*, F] = 0$ ($[p_x, F] = -\partial_xF = 0$、$[s^*, F] = \partial_zs^*\,\partial_\delta F - 0 = 0$) |
| $y$ | $s^*p_y$ | 0 (同様) |
| $p_x$、$p_y$ | $-\partial F/\partial x = 0$ | 0 |
| $z$ | $\partial F/\partial\delta = 0$ ($H_0$ は $\delta$ に依らない) | 0 |
| $\delta$ | $-\partial F/\partial z = -\tfrac12H_0 = -\tfrac14p_\perp^2$ | $[-\tfrac12H_0, F] = 0$ ($H_0$ は $x, y$ に依らず、$p$ と $z$ の関数 $F$ とは括弧が 0) |

(sympy で同じ表を作った: [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (d)。) したがって級数は 1 項で止まり、$D$ は**厳密に** $x \mathrel{+}= s^*p_x$、$y \mathrel{+}= s^*p_y$、$\delta \mathrel{-}= \tfrac14p_\perp^2$ ($s^*$ は一定。$F$ 自身は流れで保存する)。$D = \exp(-:F:)$ 型の写像で、$D^{-1}KD$ は共役なのでシンプレクティック (シンプレクティック写像の合成と逆写像はシンプレクティック)。**確認** (b): SLICE = 5、交差角 10 mrad、軌道のずれ、$\alpha\ne0$ の強ビームで、6×6 ヤコビアン (中心差分) の $\|J^\top SJ - S\|$ (規格化) は $1.5\times10^{-11}$、傾いた (x-y 結合のある) 強ビームで $4.3\times10^{-12}$、SLICE = 1・交差角なしで $7.6\times10^{-12}$。行列式は $1 \pm 5\times10^{-10}$。

## 7. 強ビームの包絡とスライス

### 7.1 6 次元の包絡

強ビームの 6×6 の 2 次モーメント $\Sigma$ ($x, p_x, y, p_y, z, \delta$) を、Twiss パラメータから作る (`bbinit`、`beambeam.setup`): $\Sigma = (HRB)\,E\,(HRB)^\top$、$B$ は $(\alpha, \beta)$ ($x$、$y$、$z$) の 2×2 のブロック、$R$ は x-y 結合 (R1..R4)、$H$ は分散 (EX, EPX, EY, EPY, ZX, ZPX, ZY, ZPY)、$E = \mathrm{diag}(\epsilon_x, \epsilon_x, \epsilon_y, \epsilon_y, \sigma_\delta\sigma_z, \sigma_\delta\sigma_z)$。これは、座標が正規化した円分布の変数 $a$ ($\langle aa^\top\rangle = E$) から $X = HRB\,a$ と書けることによる ($\Sigma = \langle XX^\top\rangle = HRB\langle aa^\top\rangle(HRB)^\top$)。$x$ のブロックでは $B\,\mathrm{diag}(\epsilon_x,\epsilon_x)B^\top = \epsilon_x\begin{pmatrix}\beta & -\alpha\\-\alpha & (1+\alpha^2)/\beta\end{pmatrix}$ で、$\sigma_{11} = \epsilon_x\beta$、$\sigma_{12} = -\alpha\epsilon_x$ (§6.2 の確認で使った関係)。縦: $\epsilon_z = \sqrt{\sigma_\delta^2\sigma_z^2 - \langle z\delta\rangle^2}$、$\beta_z = \sigma_z^2/\epsilon_z$、$\alpha_z = -\langle z\delta\rangle/\epsilon_z$。`R11..R66` がすべて正なら、直接 $\Sigma$ を与える。交差角があれば $\Sigma \to C\Sigma C^\top$ ($C$ は (16.12) の原点まわりの 1 次の係数、§8.2 で導く)。

### 7.2 スライスの条件付き包絡

$z$ を固定したときの $(x, p_x, y, p_y, \delta)$ の条件付きの分布 ($z$ 軸の方向に重みをつけたスライス) は、ガウス分布の性質から求まる。$v = (v_a, z)$ ($v_a$ は $z$ 以外の 5 成分) として、分布 $\propto\exp(-\tfrac12v^\top\Sigma^{-1}v)$ の精度行列 $\Sigma^{-1}\equiv\Lambda$ を $z$ 以外のブロック $\Lambda_{aa}$、$\Lambda_{az}$ (5 成分のベクトル)、$\Lambda_{zz}$ に分ける。指数は

$$
-\tfrac12\left[v_a^\top\Lambda_{aa}v_a + 2z\,v_a^\top\Lambda_{az} + \Lambda_{zz}z^2\right]
= -\tfrac12\left(v_a + z\Lambda_{aa}^{-1}\Lambda_{az}\right)^\top\Lambda_{aa}\left(v_a + z\Lambda_{aa}^{-1}\Lambda_{az}\right) + (z\text{ だけの項})
$$

(平方完成)。$z$ を固定すると $v_a$ はガウス分布で、共分散は **$\Sigma_5 = \Lambda_{aa}^{-1}$** (精度行列の 5×5 のブロックの逆行列)、条件付き平均は $-z\,\Lambda_{aa}^{-1}\Lambda_{az} = z\cdot(-\Sigma_5\,b)$ ($b = \Lambda_{az}$)。これは $\Sigma$ のブロックで書いた $\Sigma_{aa} - \Sigma_{az}\Sigma_{zz}^{-1}\Sigma_{za}$ (Schur 補行列)、平均 $\Sigma_{az}\Sigma_{zz}^{-1}z$ と同じ (数値で $10^{-15}$ 以下: [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (e))。K64: `m_zslice_5d` が `Benv_inv` から `m_sub5`、`v_sub5` を作り、固有値分解で $\Sigma_5 = \mathrm{Benv5}$ を作る (固有値が 0 以下ならエラー)。`v_cen` は $\langle x_i\rangle = v_i\,z$ ($i = 1\ldots4$)。このずれ (`deform`) は、各スライスで `x -= zslice*sigz*v_cen(1)` などの座標の平行移動として入れる ((16.10) の $D$ の前後)。

x-y 結合 (`Benv5(3)` の非対角) は、傾きの回転 (`benv_norm_axis`) で消す。角 $t$ の回転 $R_t = \begin{pmatrix}\cos t&\sin t\\-\sin t&\cos t\end{pmatrix}$ で $(x, y)$ の 2×2 の共分散 $\begin{pmatrix}\Sigma_{11}&\Sigma_{13}\\\Sigma_{13}&\Sigma_{33}\end{pmatrix}$ の非対角を消す条件は、$R_t\Sigma R_t^\top$ の (1,2) 成分が $\Sigma_{13}(\cos^2t - \sin^2t) + (\Sigma_{33} - \Sigma_{11})\sin t\cos t = 0$、つまり $\tan2t = 2\Sigma_{13}/(\Sigma_{11} - \Sigma_{33})$。K64 は $t = \tfrac12\arctan2(2\Sigma_{13}, \Sigma_{11} - \Sigma_{33})$ とする (数値で確認 (e))。**$\sigma_x < \sigma_y$ の強ビームは $t = \pi/2$ になる**: $\Sigma_{13} = 0$、$\Sigma_{11} < \Sigma_{33}$ なら $\arctan2(0, \text{負}) = \pi$ なので $t = \pi/2$。包絡が 90 度回され、$\sigma_x > \sigma_y$ の枝で計算される。K64 は傾きの回転に `sint = xangle(5)` ($= t$ そのもの)、`cost = xangle(6)` ($= \sin t$) を使う (`bbstrhl.f:220-221`、`xangle(6)=sin`、`xangle(7)=cos` と書いてあるのに読む位置が 1 つずれる、B34)。**確認** (a): $\sigma_x < \sigma_y$ で k64 モードは $\Delta p$ が v2 と違う: $(x, y) = (2\times10^{-5}, -10^{-5})$ で $\Delta p_x$ = $-7.9\times10^{-5}$ (k64)、$-2.6\times10^{-5}$ (v2、(16.6) と一致)。$\sigma_x > \sigma_y$ のビーム (傾き 0) は影響を受けない。

### 7.3 z のスライス

等確率の $n$ 個のスライスの境界は、ガウス分布の累積分布 $\Phi$ の逆関数で $z_i = \sigma_z\,\Phi^{-1}(i/n)$ (`gauinv`。K64 は最大相対誤差 $6.4\times10^{-4}$ の近似式)。各スライス ($z_{i-1} < z < z_i$) の重心 (σ_z 単位) は、ガウス分布の部分区間の平均で、確率 $1/n$ で割る:

$$
\bar z_i = \frac{\int_{z_{i-1}}^{z_i}z\,\phi(z)\,dz}{1/n} = n\int_{z_{i-1}}^{z_i}\frac{z\,e^{-z^2/2}}{\sqrt{2\pi}}dz = \frac{n}{\sqrt{2\pi}}\left[-e^{-z^2/2}\right]_{z_{i-1}}^{z_i} = \frac{n}{\sqrt{2\pi}}\left(e^{-z_{i-1}^2/2} - e^{-z_i^2/2}\right)
\tag{16.11}
$$

(コードは符号を反転して「前から後ろ」の順に並べる: `colb%zslice(i) = -...`。両端のスライスは $z_0 = -\infty$、$z_n = +\infty$ に当たる。数値積分との一致 ($n = 5$): [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (f)。) ルミノシティの計算 `dLum` は各スライスの 2 次元ガウス密度 $\exp(-x^2/2\sigma_x^2 - y^2/2\sigma_y^2)/(2\pi\sigma_x\sigma_y)$ の和 ([docs/14_beambeam.md](../docs/14_beambeam.md))。

## 8. 交差角・軌道のずれ・傾き

### 8.1 軌道のずれと傾き

**軌道のずれ**: 強ビームの軌道 $(x_c, p_{xc}, y_c, p_{yc}, z_c)$ を、弱粒子の座標から引いて (入口)、あとで足す (出口): `x -= cod(1)`、`px += cod(2)`、`y -= cod(3)`、`py += cod(4)`、`z += cod(5)` (`bbstrhl.f:205-215`)。強ビームの軌道の $p$ は $\cos\varphi_c$ で割って入れる (`cod(2) = DPX/cos`)。位置は引き、運動量と $z$ は足す、という符号の非対称の理由は、ソースに書かれていない。次の読みで整合する (数値で確認): 軌道 $(X, P)$ はユーザーが**強ビームの進む向きの座標系**で与える。強ビームは弱粒子と逆向きに走るので、弱粒子の座標系 ($s$ の向きが逆) で見ると、傾き $dx/ds$ と $z$ の符号が反転し、位置 $x, y$ は反転しない。弱粒子から見た強ビームの位置は $s^*$ の関数として $X - P\,s^*$ で、衝突点での相対位置は $(x - X) + (p_x + P)s^*$ になる (弱粒子の drift $x\mathrel{+}= p_xs^*$ の前に `x -= X`、`px += P` を行うと、(16.10) の $D$ のあとで正確にこれになる)。確認: SLICE = 1 の粒子 ($z = 0.004$、$s^* = 0.002$、DPX $= P = 10^{-3}$) のキックは、$x = 1.5\times10^{-4}$ で $\Delta p_x = -7.63453134\times10^{-5}$、強ビームが軸上で粒子の $x$ を $+Ps^*$ ずらしたもの ($1.52\times10^{-4}$) で同じ値 $-7.63453134\times10^{-5}$、符号を逆にして $-Ps^*$ ずらすと $-7.61791380\times10^{-5}$ で一致しない ([bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (h))。つまり「$P$ は強ビームの進む向きの傾き」という読みと矛盾しない。この約束がヘルプ・docs に書かれているかは確認していない。

**傾き**: 強ビームの x-y 結合を除く角 $t$ (§7.2) で、座標と運動量を回す: $x \to cx + sy$、$y \to cy - sx$ ($p$ も同様)。K64 は $c, s$ の読み位置を誤る (B34、gs で再現。k64 は $t \neq 0$ で誤り。v2 は $\cos t$、$\sin t$)。

### 8.2 交差角の写像 (16.12) の導出

**写像の形**: 弱ビームの座標から、強ビームが正面衝突する座標系への 6 変数の閉じた写像 ($\varphi_c$ = XANGLE、$p_n = 1+\delta$、$p_z = \sqrt{p_n^2 - p_x^2 - p_y^2}$、$h_1 = p_n - p_z$ ($=\mathtt{pxy2/(pn+pz)}$)、$t = \tan\varphi_c$)。ここで $p_x, p_y$ は正準運動量 ($p_n$ 倍したもの)。

$$
\begin{aligned}
\delta' &= \delta - \tan\varphi_c\,p_x + \tan^2\varphi_c\,h_1,\qquad p_x' = \frac{p_x - \tan\varphi_c\,h_1}{\cos\varphi_c},\qquad p_y' = \frac{p_y}{\cos\varphi_c},\\
p_z' &= \left[(1+\delta')^2 - p_x'^2 - p_y'^2\right]^{-1/2}\ (\text{逆数}),\quad
x' = \tan\varphi_c\,z + \left(1 + \sin\varphi_c\,p_x'p_z'\right)x,\\
y' &= y + \sin\varphi_c\,p_y'p_z'\,x,\qquad z' = \frac{z}{\cos\varphi_c} - \frac{\tan\varphi_c}{\cos\varphi_c}\,h_1\,p_z'\,x
\end{aligned}
\tag{16.12}
$$

(`bbstrhl.f:184-201`。ここで `pz` は $1/p_z'$ を指す。$z'$ の式の $h_1$ は変換前の値。出口は逆の写像 (`:482-498`)。)

**導出の構成**。(16.12) は、(1) 座標軸の回転、(2) 回した系の $x$ 方向の Lorentz boost、(3) 運動量の規格化の変更、(4) 変換後の面 ($s' = 0$) への直進、を順に行ったものとして再現できる。粒子は質量を無視する ($E = |p| = p_n$、$\beta\to1$。(16.12) には $\gamma$ が現れないので、K64 の式もこの近似である)。単位は設計運動量 $p_0$ (長さの単位は m、$E$ は $c$ を掛けて運動量と同じ単位)。事象は $(x, y, s, ct)$、運動量は $(p_x, p_y, p_s, E)$。$c = \cos\varphi_c$、$s = \sin\varphi_c$ (ここでの $s$ は角の正弦で、座標 $s$ ではない。座標は $S$ と書く) とする。

*運動量*:

1. **座標軸の回転** (角 $\varphi_c$、$x$-$s$ 面): $p_{x1} = c\,p_x + s\,p_z$、$p_{s1} = c\,p_z - s\,p_x$、$p_y$ と $E = p_n$ は不変。
2. **boost** ($x_1$ 方向、$\beta = \sin\varphi_c$、$\gamma = 1/\cos\varphi_c$、$\gamma\beta = \tan\varphi_c$): $p_{x2} = \gamma(p_{x1} - \beta E)$、$E_2 = \gamma(E - \beta p_{x1})$、$p_{s2} = p_{s1}$、$p_{y2} = p_y$。代入して
   $$p_{x2} = \frac{c\,p_x + s\,p_z - s\,p_n}{c} = p_x - t\,(p_n - p_z) = p_x - t\,h_1,\qquad
   E_2 = \frac{p_n - s\,(c\,p_x + s\,p_z)}{c}.$$
   また $E_2 - p_{s2} = \dfrac{p_n - s\,c\,p_x - s^2p_z - c^2p_z + s\,c\,p_x}{c} = \dfrac{p_n - p_z}{c} = \dfrac{h_1}{c}$ ($s^2 + c^2 = 1$)。
3. **規格化の変更**: 設計粒子 ($p_x = p_y = 0$、$\delta = 0$、$p_z = p_n = 1$) は、回転後 $p_{x1} = s$、boost 後 $p_{x2} = \gamma(s - s) = 0$、$E_2 = \gamma(1 - s^2) = c$ になる。すなわち boost 後の系の設計運動量は $p_0' = p_0\cos\varphi_c$。新しい変数を $p_0'$ で規格化する: $1 + \delta' = E_2/c$、$p_x' = p_{x2}/c$、$p_y' = p_y/c$。
   $$1 + \delta' = \frac{p_n - s\,c\,p_x - s^2p_z}{c^2} = p_n(1 + t^2) - t\,p_x - t^2p_z = p_n - t\,p_x + t^2(p_n - p_z)\ \Rightarrow\ \delta' = \delta - t\,p_x + t^2h_1,$$
   $$p_x' = \frac{p_x - t\,h_1}{c},\qquad p_y' = \frac{p_y}{c}.$$
   (16.12) の $\delta'$、$p_x'$、$p_y'$ と一致する ($1/c^2 = 1 + t^2$ を使った)。縦運動量は $p_s' = p_{s2}/c = (c\,p_z - s\,p_x)/c = p_z - t\,p_x$ で、質量殻の条件 $E_2^2 = p_{x2}^2 + p_y^2 + p_{s2}^2$ から $(1+\delta')^2 - p_x'^2 - p_y'^2 = p_s'^2$ なので、(16.12) の $p_z' = 1/\sqrt{(1+\delta')^2 - p_x'^2 - p_y'^2}$ は $1/p_s'$ である。また $(1+\delta') - p_s' = (E_2 - p_{s2})/c = h_1/c^2$ で、新しい系での $h_1' = h_1/\cos^2\varphi_c$。

*位置*。変換前の面 $S = 0$ にいる粒子の事象は $(x, y, 0, ct)$ で、$ct = -z$ ($z>0$ は基準粒子より先行。この符号は、結果が (16.12) の $x'$ の $\tan\varphi_c\,z$ の符号と一致するように対応づけたもので、$z$ の定義から独立に導いたものではない)。

4. 回転: $x_1 = c\,x$、$S_1 = -s\,x$。boost: $x_2 = \gamma(x_1 - \beta\,ct) = x + t\,z$、$ct_2 = \gamma(ct - \beta x_1) = -z/c - s\,x$、$S_2 = S_1 = -s\,x$。この事象は新しい面 $S_2 = 0$ にはない。
5. **直進**: 粒子は新しい運動量の向き $(p_{x2}, p_y, p_{s2})$、エネルギー $E_2$ で直進するので、$S_2 = 0$ に達するまでの経路は、$\Delta S = s\,x$ に対し、運動量に比例して $\Delta x = \Delta S\,p_{x2}/p_{s2}$、$\Delta y = \Delta S\,p_y/p_{s2}$、$\Delta(ct) = \Delta S\,E_2/p_{s2}$。$p_{x2}/p_{s2} = p_x'/p_s' = p_x'p_z'$ (規格化の因子 $c$ は約分) を使って
   $$x' = x + t\,z + s\,x\,p_x'p_z' = t\,z + (1 + s\,p_x'p_z')\,x,\qquad y' = y + s\,p_y'p_z'\,x,$$
   $$ct' = -\frac zc - s\,x + s\,x\,\frac{1+\delta'}{p_s'}\ \Rightarrow\ z' = -ct' = \frac zc - s\,x\,\frac{(1+\delta') - p_s'}{p_s'} = \frac zc - s\,x\,\frac{h_1}{c^2}\,p_z' = \frac{z}{c} - \frac tc\,h_1p_z'\,x.$$
   (16.12) の $x'$、$y'$、$z'$ と一致する ($h_1$ は変換前の値、$s/c^2 = t/c$)。

**確認**: 上の 1〜5 をそのまま sympy で書き、(16.12) との差を記号計算で 0 にした (運動量・位置とも)。さらに、Fortran の式を逐語的に写した関数と、1〜5 の関数を、質量殻上の 5 点・40 桁で比べて差 $10^{-41}$ ([bb_02_crossing_angle.out.txt](outputs/bb_02_crossing_angle.out.txt) (1)〜(3))。

**出口の写像**: 入口の逆写像である。逆 boost ($\beta\to-\beta$) と逆回転から、運動量は ($h_1'$ は出口に入る側の $h_1' = (1+\delta') - p_s'$)

$$
p_x = c\,p_x' + s\,c\,h_1',\qquad p_y = c\,p_y',\qquad \delta = \delta' + s\,p_x'
$$

(逆 boost で $E_2 = c(1+\delta')$、$p_{x2} = c\,p_x'$ から $p_{x1} = \gamma(p_{x2} + \beta E_2) = p_x' + s(1+\delta')$、$E = \gamma(E_2 + \beta p_{x2}) = 1 + \delta' + s\,p_x'$、逆回転 $p_x = c\,p_{x1} - s\,p_{s1}$ ($p_{s1} = c\,p_s'$) $= c\,p_x' + s\,c\,[(1+\delta') - p_s']$)。位置は、入口の $x'$、$z'$ の式を $x, z$ について解く。$z' = z/c - s\,h_1'p_z'\,x$ ($h_1 = c^2h_1'$ を代入) から $z = c\,z' + s\,c\,h_1'p_z'\,x$、これを $x' = t\,z + x + s\,p_x'p_z'x$ に入れて $x' - s\,z' = x\,[1 + s\,p_z'(p_x' + s\,h_1')]$。したがって

$$
x = \frac{x' - s\,z'}{1 + s\,(p_x' + s\,h_1')\,p_z'},\qquad y = y' - s\,p_y'p_z'\,x,\qquad z = c\,z' + s\,c\,h_1'p_z'\,x.
$$

これが `bbstrhl.f:482-498` の出口の式 (`x = (x - sinx*z)/(1 + sinx*(px + sinx*h1)*pz)`、`y = y - sinx*py*pz*x`、`z = z*cosx + sinx*cosx*h1*pz*x`、`h1 = h1*cosx**2`、`px = xx*cosx + tanx*h1` ($= c\,p_x' + s\,c\,h_1'$)、`g = g + sinx*xx`) と一致する。出口 ∘ 入口が恒等であることを 40 桁で確認した (差 $10^{-43}$、(3))。

**シンプレクティック性**: Lorentz 変換と回転は、拡張した位相空間 $(x, y, S, ct;\ p_x, p_y, p_s, -E)$ の点変換 (正準) で、その上で面 $S = 0$ から面 $S_2 = 0$ への写像を作ったものが上の写像である (面の間の写像は、時間の代わりに $S$ を独立変数にしたハミルトニアン系の標準的な手順で、シンプレクティックに保たれる)。ただし運動量を $p_0' = p_0\cos\varphi_c$ で規格化し直すので、入口の写像単独は共形シンプレクティックで、$J^\top SJ = S/\cos\varphi_c$ になる。実際、$\varphi_c = 0.0137$ の点で 6×6 ヤコビアンを sympy で評価した $\max|J^\top SJ - S/\cos\varphi_c|$ は $6\times10^{-46}$、$S$ との差は $1/\cos\varphi_c - 1 = 9.4\times10^{-5}$ ((4))。出口の写像は逆なので $J^\top SJ = \cos\varphi_c\,S$ で、入口 → (中間のシンプレクティックな写像) → 出口の全体は $S$ を保つ。これが §6.3 の確認 (b) の $1.5\times10^{-11}$ に反映される。

**小振幅の極限**: 原点まわりの 1 次の係数 (ヤコビアン $J(0)$、行・列の順は $x, p_x, y, p_y, z, \delta$) は、sympy で $x' = x + \tan\varphi_c\,z$、$p_x' = p_x/\cos\varphi_c$、$y' = y$、$p_y' = p_y/\cos\varphi_c$、$z' = z/\cos\varphi_c$、$\delta' = \delta - \tan\varphi_c\,p_x$ ((5))。$h_1 = p_n - p_z = p_\perp^2/(p_n + p_z) = O(p_\perp^2)$ なので $h_1$ を含む項は 1 次に現れない。特に **$p_x'$ に $\delta$ の 1 次の項はなく、$p_x' = p_x/\cos\varphi_c + O(p_\perp^2)$、小角で $p_x' \approx p_x$** (数値: $p_x = 10^{-3}$ で $p_x' - p_x/\cos\varphi_c = -6.85\times10^{-9} = -\tan\varphi_c\,p_x^2/(2\cos\varphi_c)$。$\delta = 10^{-3}$、$p_x = 0$ では $p_x' = 0$ ((5)))。小角 $\varphi_c\to0$ で $x'\approx x + \varphi_cz$、$z'\approx z$、$\delta'\approx\delta - \varphi_cp_x$。この $J(0)$ は `beambeam.setup` が強ビームの包絡に掛ける $C$ (`C[0,4] = tanx`、`C[1,1] = C[3,3] = C[4,4] = 1/cos`、`C[5,1] = -tanx`) と一致する。強ビームの包絡も同じ線形写像で運ばれていることになる。

**未導出・未確認の点**: (i) 幾何の読み。上の構成で、設計粒子は回転後 $p_{x1} = s$ (回した系の軸から $\varphi_c$ 傾く) になり、boost で $p_x' = 0$ になる。対向する強ビームは同じ boost で $x$ 方向の速度が 0 になる必要があるので、回した系の軸 (2 つのビームの二等分線) に対し、弱ビームが $+\varphi_c$、強ビームが $-\varphi_c$ の向きに走る配置、すなわち**全交差角が $2\varphi_c$ で XANGLE は半角**、という読みになる。ソースにも docs にもこの定義の記述はなく、(16.12) の構造からの逆算である。(ii) キックの係数 $c = Nr_e/\gamma$ (`ce`) には $\cos\varphi_c$ の因子がなく (`bbinit`: `ce = NP*re/gamma`)、boost 後の系で電荷分布 (密度、長さの収縮) がキックにどう入るかは導いていない。$\varphi_c = 10$ mrad で $\cos\varphi_c - 1 = 5\times10^{-5}$ の相対差。(iii) $E = |p|$ ($\beta = 1$) の近似の誤差 ($\gamma^{-2}$ の項)。(iv) $ct = -z$ の対応は、結果の符号の一致から選んだ。

## 9. ビームストラールング

BSTRL ≥ 1 のとき、各スライスのキックの大きさ $\theta = \sqrt{f_x^2 + f_y^2}$ ($= \mathtt{fx0}, \mathtt{fy0}$) から、曲率 $\rho_i = \theta/dsz$ ($dsz$ = スライスの半長 $= (z_{i-1} - z_i)\sigma_z/2$、`dsz`) とする。放射は [13_radiation.md](13_radiation.md) の kinematical method の考え方で、1 スライス 1 個以下の光子として扱う。13 章の古典の式 (光子数 $N_\gamma = \tfrac{5}{2\sqrt3}\alpha\gamma\phi$ ($\phi$ は曲げ角)、臨界エネルギー $u_c = \tfrac32\bar\lambda_C\,mc^2\gamma^3/\rho$) で、曲げ角を $\phi = \theta$ (キックの角)、$1/\rho = \rho_i$ とおく:

$$
N_\gamma = \mathtt{cnbs}\,\rho_i\,dsz = a_{\rm nrad}\,\gamma\,\theta,\qquad
\mathtt{cnbs} = \frac{5\sqrt3}{6}\alpha\gamma = a_{\rm nrad}\gamma,\qquad
\frac{u_c}{E_0} = \mathtt{cubs}\,\rho_i,\quad \mathtt{cubs} = \frac{3}{2}\frac{\hbar c\gamma^2}{m_ec^2} = c_{\rm uc}\gamma^2
\tag{16.13}
$$

($\rho_i\,dsz = \theta$ なので $\mathtt{cnbs}\,\rho_i\,dsz = a_{\rm nrad}\gamma\theta$。$u_c/E_0 = \tfrac32\bar\lambda_C\,mc^2\gamma^3\rho_i/(\gamma mc^2) = \tfrac32\bar\lambda_C\gamma^2\rho_i$。$5\sqrt3/6 = 5/(2\sqrt3)$ なので $\mathtt{cnbs} = a_{\rm nrad}\gamma$、$c_{\rm uc} = \tfrac32\bar\lambda_C$。**確認** (c): 比 $1.000000000$ と $0.99999997$。) 光子を出すかは乱数 $x_r\in(0,1)$ で $x_r < N_\gamma$ (`rnbs`)。$N_\gamma > 0.5$ は K64 が停止する ("Nbs is too large")。出すなら、$x_r/N_\gamma$ ($(0,1)$ に一様、乱数の再利用) を、光子数のスペクトル $dN/d\xi \propto S(\xi) = \int_\xi^\infty K_{5/3}(x)dx$ ($\xi = u/u_c$) の累積分布の逆引き表 (`besseltab`: $\xi_i = 1.1^{40}/1.02^i$、累積 `SwI`) で $\xi_i$ に直し、$\delta \mathrel{-}= \mathtt{cubs}\,\rho_i\,\xi$。この分布の平均と 2 次のモーメントは $\langle\xi\rangle = c_{\rm ave}$、$\langle\xi^2\rangle = c_{\rm uu}$ に等しいはず。値は次のように Bessel 関数の積分から出る。$S(\xi)$ の部分積分 $\int_0^\infty\xi^\mu S(\xi)d\xi = \dfrac1{\mu+1}\int_0^\infty x^{\mu+1}K_{5/3}(x)dx$ と、積分公式 $\int_0^\infty x^\nu K_{5/3}(x)dx = 2^{\nu-1}\Gamma\!\left(\frac{\nu+8/3}2\right)\Gamma\!\left(\frac{\nu-2/3}2\right)$ から、$\langle\xi\rangle = \dfrac{\tfrac12M_2}{M_1}$、$\langle\xi^2\rangle = \dfrac{\tfrac13M_3}{M_1}$ ($M_\nu$ は上の積分) $= 8/(15\sqrt3) = 0.30792$、$11/27 = 0.40741$ (mpmath で $10^{-12}$ まで確認: [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) (g))。

**平均の損失**: $N_\gamma\,c_{\rm ave}\,u_c/E_0 = a_{\rm nrad}c_{\rm ave}c_{\rm uc}\gamma^3\theta^2/dsz = \tfrac23r_e\gamma^3\theta^2/dsz$。係数は $a_{\rm nrad}c_{\rm ave}c_{\rm uc} = \tfrac{5\alpha}{2\sqrt3}\cdot\tfrac{8}{15\sqrt3}\cdot\tfrac32\bar\lambda_C = \tfrac{5\cdot8\cdot3}{2\sqrt3\cdot15\sqrt3\cdot2}\alpha\bar\lambda_C = \tfrac{120}{180}\alpha\bar\lambda_C = \tfrac23r_e$ ($\alpha\bar\lambda_C = r_e$。数値 $0.666666666667$、(g))。古典の $\tfrac23r_e\gamma^3L/\rho^2$ ($\theta = L/\rho$、$L = dsz$ として $\theta^2/dsz = L/\rho^2$) と一致 (`cpbs` = $\tfrac23r_e\gamma^3$ は定義されているが、平均の損失としては使わない)。

**確認** (c): 表を $x_r\in(0,1)$ の一様な 20 万点で逆引きして得たモーメントは $\langle\xi\rangle = 0.3128$、$\langle\xi^2\rangle = 0.4180$ (理論 $0.3079$、$0.4074$ に $+1.6\%$、$+2.6\%$)。同じ格子で累積を台形則で作り直して同じ手順にすると $0.3098$、$0.4098$ ($+0.6\%$、$+0.6\%$)。差は表の作り方 (矩形則に近い累積、2% 刻み) の誤差とみられる (推定。§13 疑い (小)、影響は光子エネルギーの平均 ~1%)。

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
| 軌道のずれ | `bbstrhl.f:205-215`、`:470-478` | `beambeam` | `cod` (位置は引き、運動量と $z$ は足す) |
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
| `und_02_derivations.py` | (a) 条件を使わない $\vec B$ の成分・$\nabla\cdot\vec B$・$\nabla\times\vec B$ の因子、(b) 周期の積分・端の $A$、(c) $H$ の展開と運動方程式、(d) 陰的な式の 2×2 の解と 6×6 ヤコビアンのシンプレクティック性、(e) $-dv\,L$ と $L/\gamma^2$ | [und_02_derivations.out.txt](outputs/und_02_derivations.out.txt) |
| `bb_01_kick_slice.py` | (a) 閉じた形のキックと $\delta$ のキック (ポテンシャルの数値積分)、傾きの k64・v2、(b) 6 次元写像のシンプレクティック性、往復の恒等、(c) ビームストラールングの定数と光子エネルギーの表 | [bb_01_kick_slice.out.txt](outputs/bb_01_kick_slice.out.txt) |
| `bb_02_crossing_angle.py` | 交差角の写像 (16.12) の導出 (回転 + boost + 規格化 + 直進) と Fortran の式の一致 (記号と 40 桁)、出口 = 入口の逆、$J^\top SJ = S/\cos\varphi_c$、小振幅の極限 | [bb_02_crossing_angle.out.txt](outputs/bb_02_crossing_angle.out.txt) |
| `bb_03_kick_derivations.py` | (a) $\nabla^2U = 4\pi\rho$、(b) 丸いビームの $\mathrm{Ein}$・$\partial U_{\rm phys}/\partial\sigma^2$・楕円の遠方、(c) (16.8) の導出の各段、(d) $D$ のリー級数、(e) 条件付きガウス分布・傾き角、(f) (16.11)、(g) ビームストラールングの係数とモーメント、(h) 軌道のずれの符号 | [bb_03_kick_derivations.out.txt](outputs/bb_03_kick_derivations.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/und_01_field_map.py` (数秒)、`py -3.12 physics_guide/scripts/und_02_derivations.py` (20 秒程度)、`py -3.12 -W ignore physics_guide/scripts/bb_01_kick_slice.py` (数秒。scipy の積分の警告を抑える)、`py -3.12 physics_guide/scripts/bb_02_crossing_angle.py` (数秒)、`py -3.12 -W ignore physics_guide/scripts/bb_03_kick_derivations.py` (1 分程度)。BEAMBEAM の gs との照合・光子リスト・ルミノシティ・ビーム・ビーム・チューンシフトの図は [docs/14_beambeam.md](../docs/14_beambeam.md)。

## 13. 未確認・限界・使った条件

**使った条件**: UND は 3 GeV、$L = 1$ m、10 周期、$B_y = 1$ T、$B_x = 0.4$ T (確認 (b))。BEAMBEAM は 1.7 GeV の電子 ($\gamma = 3326.8$)、$N = 10^{10}$、$\sigma = 0.05$〜0.2 mm、$\sigma_z = 5$〜10 mm、SLICE = 1〜5 (確認 (a)(b))。交差角の導出の数値確認は $\varphi_c = 0.0137$ rad (質量殻上の 5 点)。実機の衝突点 (β* が小さい、数百万粒子) の比較はしていない。

**未確認**:
- UND: 端の処理 (終端の極)、傍軸近似の誤差 (§3.2 に見積り: $\pi^2/4p_n^2\sim10^{-6}$。$K$ がもっと大きいときは未確認)、放射 (v2 だけ)、スピンの結合 (v2 の POL)、`Pole` の数が整数でないとき。刻み数の選び方の指針。$p_x$ の誤差が刻み数に対して単調でない理由 (§1)。
- BEAMBEAM の交差角 (§8.2): 幾何の読み (XANGLE が全交差角の半分) はソースにもヘルプにも記述がなく、(16.12) の構造からの逆算。キックの係数 $c$ に $\cos\varphi_c$ の因子がなく、boost 後の電荷分布 (密度) の変換との整合は未導出。$\beta = 1$ の近似の誤差 ($\gamma^{-2}$) は評価していない。$ct = -z$ の対応は結果の符号から選んだ。
- BEAMBEAM の軌道のずれ (§8.1): 位置は引き、運動量と $z$ は足す符号の理由は「強ビームの軌道は強ビームの進む向きの座標で与える」という読みで、数値の整合 (§8.1) はとれたが、ソースに記述がなく、$z$ の符号 ($z$ に足す $\mathtt{cod(5)}$) の独立な確認はしていない。強-強の自己無撞着な計算 (弱-強だけ)。`cerrf` の精度 (単精度の定数 3.2)。スライス数の収束 (本章は SLICE = 1〜5 で式の確認のみ、ルミノシティ・チューンシフトの収束は未確認)。スライスの重心 `deform` の符号 ($x$ と $p_x$ ともに引く) の理由 (軌道のずれの符号と同じ流儀かどうか)。
- ビームストラールング: 光子エネルギーの表の精度の原因 (矩形則に近い累積という推定を表の作り方から確認していない)、$N_\gamma > 0.5$ の停止、光子 1 個の近似 (1 スライスに 2 個以上)、光子の放射角。
- SOL 領域内の BEAMBEAM は未対応。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| UND・BEAMBEAM | **ヘルプに記述がない** (SADHelp の要素の一覧・節に `UND`・`BEAMBEAM` の項目はなく、`UND` の文字列は `FOUND` などの一部だけ)。式・単位はすべて K64 のソースによる。 |
| ヘルプに無いもの | (16.1)〜(16.13) すべて。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

根拠の程度を次の 3 段に分ける: **gs で再現** (K64 の実行で確認)、**ソースの読み** (Fortran の式を読んで確認。gs では確認できていない)、**疑い** (意図かどうか未確認)。

**gs で再現 (docs/30)**:
- **B37**: UND の磁場のキーワード (`FBX`・`FBY`・`BX`・`BY`) が K64 では設定できない (`kwFBX` が未定義)。そのため $F_0 = G_0 = 0$ の磁場なしの追跡になる。v2 は `FBX`・`FBY` を受ける。
- **B34**: BEAMBEAM の傾きの回転が $c, s$ の読み位置を誤る (`sint = xangle(5)`、`cost = xangle(6) = sin t`)。**$\sigma_x < \sigma_y$ の強ビームは傾き $t = \pi/2$ になるので、k64 では $\Delta p$ が誤る** (本章確認 (a)、$-7.9\times10^{-5}$ 対 $-2.6\times10^{-5}$)。

**ソースの読み (gs では確認できていない)**:
- **B78**: `undulator.f:92-93` で `dpds` を二重に掛ける。(16.5) の右辺の $A\,\partial A$ の項が $ds/p_n$ 倍に弱まり、磁場があっても収束の項が刻み幅に比例して消える (docs/30 は、刻み数を変えた数値実験で、200 と 2000 スライスで $p_y$ が 10 倍違うとする)。gs では UND の磁場キーワードが無く (B37)、K64 の実行では確認できない。K64 では UND に磁場が無いので現状は効かない。本章の確認 (b) は v2 の式 (`dpds` を 1 回) で、DOP853 との一致 (1 次) を示す。

**疑い**:
- UND の $z$ の更新に $-dv\,L$ (速度が $v_0$ と違うための遅れ、$R_{56} = L/\gamma^2$) がない (磁場なしの UND と DRIFT で $z$ が $dv\,L$ 違う、本章確認 (d)、手計算の $\delta/\gamma^2\cdot L = 2.90\times10^{-11}$ とも一致)。DRIFT は持つ (02 章)。$H_{\rm UND}$ の定数項 $-(1+\delta)$ と $E/(p_0v_0)$ の項を落としているため。意図的かどうか未確認。
- ビームストラールングの光子エネルギーの表: 累積の矩形則とみられる (推定) 誤差で、平均が $+1.6\%$、2 次モーメントが $+2.6\%$ (理論比)。影響は小さい。

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
| `pxy2`、`h1`、`pz` | $p_\perp^2$、$h_1 = p_n - p_z$、(16.12) で $1/p_z'$ | `beambeam` | |
