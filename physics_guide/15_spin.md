# 15. スピンと偏極: Thomas-BMT から Sokolov-Ternov、脱偏極まで

**この章で分かること**: SAD の `POL` フラグ (K64 の `calpol`) が扱うスピンの物理を、出発点から組み立てる。(1) Thomas-BMT 方程式を、静止系の歳差と Thomas 歳差から導き、磁場の縦成分・横成分への分解 $(1+a\gamma)$、$(1+a)$ を出す。(2) 1 つの放射のスライス (`tradk1`・`tradke`) でのスピン回転 `sprot` ($g = (1+a\gamma)b_t + (1+a)b_l$、Rodrigues の式)。(3) 軌道との結合 $D$ (`srot(3,9)` の 4〜9 列) の積算。(4) 1 周の行列から偏極ベクトル $n_0$、スピンチューンを出す `spnorm`。(5) SLIM (線形のスピン軌道結合) の $N$ の式と脱偏極率。(6) Sokolov-Ternov の定数 (`pst`・`cfpd`・`sflc`) が `anrad`・`cuc`・`cave`・`cuu` から出ること、偏極時間、平衡偏極度、Derbenev-Kondratenko (DK) の式との対応。

**範囲**: 古典のスピン (ベクトル $S$) の運動と、放射による偏極の立ち上がり・脱偏極を、K64 の `photons.f` (`sprot`、`tradke`、`spnorm`、`srequpol`) の式に沿って扱う。Sokolov-Ternov の遷移率そのもの (QED) は導かず、入力として使い、定数が K64 の定数と一致することを確かめる。高次のスピン軌道共鳴、運動論的偏極の導出、側帯の項の詳細は、結果の要約と [docs/16_spin_pol.md](../docs/16_spin_pol.md) への参照にとどめる。

**前提**: [13_radiation.md](13_radiation.md) (`anrad`、`cuc`、`cave`、`cuu`、`tradke`)、[14_emittance.md](14_emittance.md) (`trans(6,12)`、`transr`)、[12_linear_optics.md](12_linear_optics.md) (1 周の行列)。

**読む順**: §1 (結論) → §3 (Thomas-BMT) → §4 → §7 (Sokolov-Ternov) → §11 (確認)。SLIM を知りたいときは §5〜§6。

## 1. 結論

- **Thomas-BMT** (電場なし): 実験室系の時間で $\dfrac{d\vec S}{dt} = \vec\Omega\times\vec S$、$\vec\Omega = -\dfrac{q}{\gamma m}\left[(1+a\gamma)\vec B_\perp + (1+a)\vec B_\parallel\right]$ ($a = (g-2)/2$、$\perp$・$\parallel$ は運動量に対して)。これは静止系の歳差と Thomas 歳差の和で、sympy で成分ごとに一致する ([spin_01_thomas_bmt.out.txt](outputs/spin_01_thomas_bmt.out.txt) (a))。経路長 $s$ で書くと、運動量の向きの回転ベクトル $\vec\theta'$ と縦磁場の積分 $\psi'$ から $d\vec S/ds = [(1+a\gamma)\vec\theta' + (1+a)\psi'\hat p]\times\vec S$。運動量に対する相対では、水平面内の偏向で $a\gamma\theta$ だけ回る (確認 (b2): 3 GeV で $\theta = 0.2618$ と $1.0$、差 $10^{-13}$)。
- **1 ステップの回転** (`sprot`): 運動量の変化 $\to$ $b_t$、縦磁場 $\to$ $b_l$ とし、$g = (1+a\gamma)\,b_t + (1+a)\,b_l$ のまわりに $|g|$ だけ回す (K64 の約束: $S' = \exp(-[g]_\times)S$)。刻みを 2 倍にすると誤差は 1/4 (2 次、確認 (b))。K64 は偏向角 $\theta$ の代わりに弦 $2\sin(\theta/2)$ を使い、相対 $-\theta^2/24$ (1 刻みで $-2.9\times10^{-3}$) だけ小さい (v2 は直す)。
- **軌道との結合 $D$**: 1 周の回転を座標の関数として $R(x_0) = \exp(-[D x_0]_\times)R_0$ と書く ($D$: 3×6、K64 の `srot(:,4:9)`)。1 スライスごとに $D \leftarrow R_{\rm step}\,(D + \partial g/\partial x\cdot M_{0\to i})$ と積算する (`tradke`)。
- **偏極ベクトルとスピンチューン**: $R_0 = \mathtt{srot(:,1:3)}$ の周期解の軸が $n_0$ (`spnorm`)、歳差角 $\mu_s$ から $\nu_s = \mu_s/2\pi$。平面リングでは $n_0 = \hat y$、$\nu_s = a\gamma\cdot\theta_{\rm tot}/2\pi$。確認 (e): 3 GeV・偏向の合計 $\pi$ のリングで SpinTune $= 0.40407411$ (= $a\gamma\pi/2\pi = 3.404074$ の端数 $0.404074$)。
- **SLIM**: スピン軸の座標依存 $n(x) = n_0 + Nx$ は $N M = R_0 N + [n_0]_\times D$ を満たし (導出 §6、確認 (d): 1 次で成り立つ)、放射の跳び $\xi$ ($\langle\xi\xi^\top\rangle = B_i$) による偏極の減りは 1 回あたり $\tfrac12|N\xi|^2$、脱偏極率 $= \tfrac12\sum_i\mathrm{tr}(N_iB_iN_i^\top)$ (1 周あたり)。
- **Sokolov-Ternov**: K64 の定数は $c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2 = r_e\bar\lambda_C$、$s_{\rm flc}\,a_{\rm nrad} = c_{\rm fpd}\,r_e\bar\lambda_C$ ($c_{\rm fpd} = 5\sqrt3/8$) で、相対 $10^{-12}$ で一致する (確認 (c))。偏極時間 $1/\tau_{ST} = \tfrac{5\sqrt3}{8}r_e\bar\lambda_Cc\gamma^5\langle\rho^{-3}\rangle$ は、平面リングで `PolarizationTime` $= 123.8153$ s と 6 桁以上一致。平衡偏極度は $8/(5\sqrt3) = 0.92376$ (v2 の `Polarization` $= 0.923760431$)。
- **DK との対応**: SLIM の $\delta$ の寄与 $\tfrac12 B_{\delta\delta}$ の係数 $\tfrac12 c_{\rm uu}\,a_{\rm nrad}c_{\rm uc}^2/(s_{\rm flc}\,a_{\rm nrad}) = 11/18$ は、DK の式の $\tfrac{11}{18}|\gamma\,\partial n/\partial\gamma|^2$ の係数に一致する (確認 (c))。
- **K64 の限界**: `Polarization` は k64 で NaN (`spdepol` の `cerfc` の桁あふれ、B55)。v2 は SLIM + 側帯の式で置き換えた。高エネルギー・大きな RF 電圧では v2 の脱偏極率が追跡 (FLUC) の 0.4 倍に過小評価される例があり ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §5.3)、確認は追跡で行う。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| $a$ | 異常磁気能率 $(g-2)/2 = 1.15965218\times10^{-3}$ | `gspin` (`tfstk.f:269`)、`cl` $= 1 + a$ | `spin.GSPIN`、`spin.CL` |
| $\gamma$ | ローレンツ因子 | `h1` (`p2h(p)`、$p = p_0(1+\delta)$)。`ct` $= 1 + h_1\,\mathtt{gspin}$ | |
| $S$ | スピンの単位ベクトル | `sx, sy, sz` | `_Particle.s` |
| $g$ | 1 ステップの回転ベクトル | `gx, gy, gz`、$|g|$ = `g` | |
| $b_t$、$b_l$ | 運動量に垂直・平行な磁場の積分 (運動量で規格化) | `btx..`、`blx..` | |
| $\psi$ | 縦磁場の積分 | `bsi`、`bsir0` (`temw`)、`bzh` (SOL) | `_Particle.bsi` |
| `ppx, ppy, ppz` | $p_{\rm out}\times p_{\rm in}$ (規格化前) | `ppx, ppy, ppz`、`xpx, xpy, xpz` | |
| $R$ | 1 周のスピン回転行列 | `srot(3,9)` の 1:3 列 | `st` (3×9) |
| $D$ | スピン軌道結合 (3×6) | `srot(3,9)` の 4:9 列 | `st[:, 3:9]` |
| $n_0$ | 偏極ベクトル | `sps(:,1)` (`spnorm`) | `PolarizationVector` |
| $\nu_s$ | スピンチューン | `smu` $= 2\pi\nu_s$ | `SpinTune` |
| 偏極の立ち上がり率 | Sokolov-Ternov の 1 周あたりの率 | `sdamp` $= n_0\cdot\mathtt{gintd}$。$\tau = 1/(\mathtt{sdamp}\,f_{\rm rev})$ | `PolarizationTime` |
| $P_{ST}$、$c_{\rm fpd}$、$s_{\rm flc}$ | $8\sqrt3/15$、$5\sqrt3/8$、$0.75\,(r_e/\alpha)^2$ | `pst`、`cfpd`、`sflc` (`photons.f:379-380`) | `PST`、`CFPD`、`sflc()` |
| $N$ | $\partial n/\partial x$ (3×6) | `drot` から (`srequpol`) | `slim_depolarization` |

SAD の座標は $(x, p_x, y, p_y, z, \delta)$。スピンは (sx, sy, sz) = 要素の局所座標 ($z$ = 進行方向) の成分で、偏向磁石では各スライスの偏向角 `phi0` ($\cos$ = `cphi0`、$\sin$ = `sphi0`) の分だけ座標系ごと回す (§4.3)。

## 3. Thomas-BMT 方程式

出発点は、静止系での古典のスピンの運動と、Lorentz 変換 (Thomas 歳差) だけである。単位は $c = 1$、電荷 $q$、質量 $m$、$\vec\beta = \vec v$。$g$ 因子の異常分は $a = (g-2)/2$、$g/2 = 1 + a$。

### 3.1 静止系の歳差

静止系 (記号に $'$) では、磁気モーメント $\vec\mu = \dfrac{gq}{2m}\vec S$ が磁場 $\vec B'$ から受けるトルクで $\dfrac{d\vec S}{dt'} = \vec\mu\times\vec B' = -\dfrac{gq}{2m}\vec B'\times\vec S$。回転角速度 $\vec\Omega' = -\dfrac{gq}{2m}\vec B'$ である。電場のない実験室系の磁場 $\vec B$ から静止系の磁場への変換は、

$$
\vec B' = \gamma\left[\vec B - \frac{\gamma}{\gamma+1}\vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.1}
$$

(縦成分は $B'_\parallel = B_\parallel$、横成分は $B'_\perp = \gamma B_\perp$。(15.1) の右辺で縦成分: $\gamma(1 - \gamma\beta^2/(\gamma+1))B_\parallel = B_\parallel$ ($\gamma\beta^2 = (\gamma^2-1)/\gamma$ を使う)。) 実験室系の時間 $dt = \gamma\,dt'$ で書くと、

$$
\vec\Omega_{\rm rest} = -\frac{(1+a)q}{m}\left[\vec B - \frac{\gamma}{\gamma+1}\vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.2}
$$

($\gamma$ が打ち消し合う。)

### 3.2 Thomas 歳差

加速される粒子の静止系の軸は、実験室系に対して、各瞬間の Lorentz ブーストの合成のために回転する。この角速度は $\vec\omega_T = \dfrac{\gamma^2}{\gamma+1}\,\vec w\times\vec v$ ($\vec w$ は加速度)。磁場だけの力 $\vec w = \dfrac{q}{\gamma m}\vec v\times\vec B$ を入れ、$(\vec v\times\vec B)\times\vec v = \beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)$ を使うと、

$$
\vec\omega_T = \frac{q\gamma}{m(\gamma+1)}\left[\beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.3}
$$

### 3.3 和: BMT

実験室系の歳差は (15.2) + (15.3):

$$
\vec\Omega = -\frac{q}{m}\left[\left(a + \frac1\gamma\right)\vec B - \frac{a\gamma}{\gamma+1}(\vec\beta\cdot\vec B)\vec\beta\right]
\tag{15.4}
$$

(整理: $\vec B$ の係数は $-(1+a) + \gamma\beta^2/(\gamma+1) = -(a + 1/\gamma)$ ($\gamma\beta^2/(\gamma+1) = (\gamma-1)/\gamma$ より)。$\vec\beta(\vec\beta\cdot\vec B)$ の係数は $(1+a)\gamma/(\gamma+1) - \gamma/(\gamma+1) = a\gamma/(\gamma+1)$。) これが Thomas-BMT 方程式 (電場なし) である。**確認** (a): sympy で $\vec\Omega_{\rm rest} + \vec\omega_T - \vec\Omega_{\rm BMT}$ の 3 成分が 0 になる。

### 3.4 縦・横への分解

磁場を運動量に平行 ($\parallel$) と垂直 ($\perp$) に分ける。$\vec\beta\cdot\vec B = \beta B_\parallel$ なので、平行成分の係数は $a + \dfrac1\gamma - \dfrac{a\gamma\beta^2}{\gamma+1} = a + \dfrac1\gamma - \dfrac{a(\gamma-1)}{\gamma} = \dfrac{1+a}{\gamma}$、垂直成分の係数は $a + \dfrac1\gamma = \dfrac{1+a\gamma}{\gamma}$:

$$
\vec\Omega = -\frac{q}{\gamma m}\left[(1+a\gamma)\vec B_\perp + (1+a)\vec B_\parallel\right]
\tag{15.5}
$$

(sympy の出力が同じ形。) 運動量の運動方程式 $d\vec p/dt = q\vec v\times\vec B$ から、運動量の向き $\hat p$ は $d\hat p/dt = -\dfrac{q}{\gamma m}\vec B_\perp\times\hat p$ で回る。(15.5) の $\perp$ の部分は、運動量の回転と同じ向きに $(1+a\gamma)$ 倍の速さで回ることを表す。運動量に対する相対では $a\gamma$ 倍。

### 3.5 経路長で書く

$dt = ds/v$、$p = \gamma mv$ で、(15.5) は $\vec\Omega\,dt = \vec\omega_s\,ds$:

$$
\frac{d\vec S}{ds} = \vec\omega_s\times\vec S,\qquad \vec\omega_s = (1+a\gamma)\,\vec\theta' + (1+a)\,\psi'\,\hat p,\qquad \vec\theta' \equiv -\frac{q}{p}\vec B_\perp,\quad \psi' \equiv -\frac{q}{p}B_\parallel
\tag{15.6}
$$

$\vec\theta'$ は運動量の向きの回転ベクトル (単位長さあたり、軸 × 角)、$\psi'$ は符号つきの縦磁場の積分 $\int B_\parallel ds/(B\rho)$ の密度 ($B\rho = p/q$)。SAD の $b_t$ が $\vec\theta'$ の積分、$b_l$ が $\psi'$ の積分に当たる (§4)。

**例 (平面リング)**: 垂直の一様磁場だけなら $\vec\theta'$ は一定の軸 $\hat y$ のまわり、運動量に対するスピンの相対歳差は $a\gamma\,\theta$ ($\theta$ は偏向角)。完全なリング ($\theta = 2\pi$) で $\nu_s = a\gamma$。$\gamma = 5870.85$ (3 GeV) で $a\gamma = 6.808$ (`NominalSpinTune` $= h_0\cdot g = 6.80814822$)。
**縦磁場**: $\hat p$ に平行な磁場だけなら運動量は回らず、スピンは $(1+a)\psi$ だけ回る (確認 (b2): $\psi = 0.5$ で回転角 $0.50057983 = (1+a)\psi$)。

## 4. 1 要素 (スライス) のスピン回転

### 4.1 K64 の `sprot`

放射の 1 スライス (`tradk1`) の前後の運動量 $p_{\rm in}$、$p_{\rm out}$ ($p_x, p_y$ と $p_z = 1 + \mathtt{pxy2dpz}$) から、

$$
(\mathtt{ppx}, \mathtt{ppy}, \mathtt{ppz}) = \vec p_{\rm out}\times\vec p_{\rm in}
\tag{15.7}
$$

(`tradk1` では `ppx = py*dpz0 - dpz*py0 + dpy` など。$p_z = 1 + dp_z$ を使った外積の成分で、大きさは $\sin\theta$、向きは運動量の回転軸の**逆**。) `sprot` は、中点の運動量 $\hat p_m = (\mathtt{pxm}, \mathtt{pym}, \mathtt{pzm})$ ($\mathtt{pzm} = 1 + \mathtt{pxy2dpz}(\mathtt{pxm}, \mathtt{pym})$) で割って

$$
\vec b = \frac{(\mathtt{ppx}, \mathtt{ppy}, \mathtt{ppz})}{\mathtt{pzm}} + (0, 0, \mathtt{bsi}),\quad
b_p = \vec b\cdot\hat p_m,\quad \vec b_l = b_p\,\hat p_m,\quad \vec b_t = \vec b - \vec b_l
\tag{15.8}
$$

として、回転ベクトルを

$$
\vec g = (1 + a\gamma)\,\vec b_t + (1 + a)\,\vec b_l\qquad (\texttt{ct} = 1 + h\,\mathtt{gspin},\ \texttt{cl} = 1 + \mathtt{gspin})
\tag{15.9}
$$

とする。(15.6) との対応: $\vec b_t = -\vec\theta$ (運動量の回転ベクトルの符号を反転したもの。ppx.. が $p_{\rm out}\times p_{\rm in}$ だから)、$\vec b_l = \hat p\,\mathtt{bsi}$。K64 は**スピンを $\exp(-[g]_\times)$ で回す**約束で、$-\vec g$ が右手系の回転ベクトルになり、(15.6) の $\vec\omega_s\,ds$ と一致する。$|g|$ の割り算: 偏向の 1 ステップで $\mathtt{ppx..}/\mathtt{pzm} = \sin\theta/\cos(\theta/2) = 2\sin(\theta/2)$ (弦の長さ) になる。

回転は Rodrigues の式。`xsincos(g, sinu, xsinu, cosu, dcosu)` が $\sin|g|$、$\cos|g|$、$\cos|g| - 1$ ($= -2\sin^2(|g|/2)$、桁落ちを避ける) を返し、

$$
\vec S' = \cos|g|\,\vec S + \frac{\sin|g|}{|g|}\,\vec S\times\vec g + \frac{1 - \cos|g|}{|g|^2}(\vec S\cdot\vec g)\,\vec g
\tag{15.10}
$$

(`sw = -(S·g) dcosu/g²`、`sinu = sin|g|/|g|`、`sux = S×g`)。(15.10) は $\exp(-[g]_\times)\vec S$ の閉じた式 ($-\vec g\times\vec S = \vec S\times\vec g$ だから)。

### 4.2 刻みの誤差

(15.9)〜(15.10) の列は、Thomas-BMT の厳密な積分に収束する。**確認** (b): $\hat B$ と運動量が $0.3$ rad の角をなす一様磁場 (らせん軌道、3 GeV、$\kappa L = 2$ rad、縦成分あり) で、(15.6) の DOP853 による積分を基準にすると、$\vec\omega_s$ の刻みごとの回転 $\mathrm{Rot}[(1+a\gamma)\vec\theta_{\rm step} + (1+a)\psi\,\hat p_m]$ の積は、刻み数 $N = 10, 20, 40, 80, 160$ で誤差 $1.08\times10^{-2}$、$2.69\times10^{-3}$、$6.73\times10^{-4}$、$1.68\times10^{-4}$、$4.21\times10^{-5}$、誤差の比はいつも $4.00$ で **2 次**。中点の $\hat p_m$ を縦成分の向きに使い、中点則になっているため。

**K64 の弦**: K64 は $\theta$ のかわりに弦 $2\sin(\theta/2)$ を使う。$\theta = 0.2618$ rad の偏向の 1・3・10 刻みで、歳差の角の差は $-5.83\times10^{-3}$、$-6.49\times10^{-4}$、$-5.84\times10^{-5}$ rad、相対 $-\theta_n^2/24$ ($\theta_n = \theta/n$) と一致する。v2 は $\sin\theta$ の外積から $\theta = \arcsin$ を取り直して厳密な $\theta$ にする (B83)。

### 4.3 座標系の回転と縦磁場 `bsi`

偏向磁石の 1 スライスでは、スピンを回転したあと、スライスの偏向角 $\phi_0$ の分だけ座標系を回す:

$$
s_x' = c\,s_x + s\,s_z,\qquad s_z' = (s_z - s\,s_x')/c\qquad (c = \mathtt{cphi0},\ s = \mathtt{sphi0})
\tag{15.11}
$$

(`sprot` 最後の 2 行。$s_z' = (s_z - s\,s_x')/c$ は、$s_x'$ の更新後の値を使うと $-s\,s_x'/c + s_z/c = (s_z - s(cs_x + s s_z))/c = c\,s_z - s\,s_x$ になり、正しい回転 $s_z' = -s\,s_x + c\,s_z$ と等しい。) 縦磁場の積分 `bsi` は、偏向磁石の端の縦磁場 ($\psi$ = `bsi1*yi/rhob` など)、MULT、SOL 領域の $b_z\,al$ で足される ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §3.2、§4)。ここでは、いずれも (15.9) の $\vec b_l$ に入る量だと理解すればよい。

## 5. 軌道との結合 $D$ と 1 周の行列

### 5.1 `tradke` の積算

Emittance[] 側 (`tradke`、`photons.f:1521-1616`) は、同じ $\vec g$ (15.9) を閉軌道のまわりで作り、1 周の回転と**座標への微分**を同時に積算する。`srot(3,9)` の 1:3 列が $R$ (3×3)、4:9 列が $D$ (3×6)。スライスごとに:

1. $\vec g$ と、スライスの入口の座標 $x_{\rm in}$ (6 成分) に関する偏微分 $\partial\vec g/\partial x_{\rm in}$ (`dgx, dgy, dgz`、6 成分ずつ) を作る。$\vec b_t$ は運動量の変化から、運動量の変化は 1 スライスの行列 `tr2` ($= \mathtt{trans}\cdot\mathtt{transr}^{-1}$) を通して $x_{\rm in}$ に依存する (`dxpx = tr2(4,:)*pz0 + py*dpz0 - ddpz*pyi - pz*dpyi` など。偏微分を 1 項ずつ書いた式)。$\gamma = h_1$ の $\delta$ 依存は `dct = ct*dh1r` ($\mathtt{dh1r} = pp_0/h_1^2$、$d\gamma/d\delta = pp_0/\gamma$ から)。
2. $x_{\rm in}$ から $x_0$ (リングの入口の座標) への写像 `transr` ($= M_{0\to i}$、前のスライスまでの累積行列) で $x_0$ に関する微分にして、$D$ に足す (`srot(1,4:9) += dgx(1:6)·transr`、`srot(2,..)`、`srot(3,..)` も同様)。
3. $\vec g$ のまわりの回転を `srot` の**全部の列**にかける (`srot(1,:) = cosu*sx + sinu*sux + sw*gx`、$\vec s_x$ は `srot(1,:)` の行)。
4. 座標系を $\phi_0$ だけ回す ((15.11) と同じ)。

**なぜこれで $D$ になるか**: スライス $i$ で座標が $x_0$ から出発した粒子の $\vec g$ が $\vec g_i + \delta\vec g_i$ ($\delta\vec g_i = (\partial\vec g_i/\partial x_{\rm in})M_{0\to i}x_0$) だとすると、$\mathrm{Rot}(g_i + \delta g_i) \approx \mathrm{Rot}(g_i)\exp(-[\delta g_i]_\times)$ (1 スライスの $|g|$ は小さいので、回転の非可換な補正は 2 次の微小量)。残りのスライスの回転 $R_N\cdots R_{i+1}$ を左からかけると、$\exp(-[R_N\cdots R_{i+1}\delta g_i]_\times)$ の形で全体の $R_{\rm tot}$ の左に出る:

$$
R(x_0) \approx \exp\!\left(-[D\,x_0]_\times\right)R_0,\qquad D = \sum_i R_N\cdots R_{i+1}\,\frac{\partial\vec g_i}{\partial x_{\rm in}}M_{0\to i}
\tag{15.12}
$$

(一般に $B\exp(-[v]_\times)B^{-1} = \exp(-[Bv]_\times)$ ($B$ は回転行列) による。) 手順 3 が「残りの回転を $D$ の列にかける」ことに、手順 2 が各スライスの寄与を足すことに当たる。したがって、**$D$ は 1 周の回転に対する座標の影響を表す 3×6 行列で、$R$ と同じ座標系 (1 周の終わり) で測る**。$\partial\vec g/\partial x$ の中の $\partial b_t/\partial x$ は偏向の位置の依存 (勾配・フリンジ)、$\partial/\partial\delta$ の項は $a\gamma$ の $\delta$ 依存と運動量 ($B\rho$) の変化から来る。

### 5.2 `spnorm`: 偏極ベクトルとスピンチューン

$R_0 = \mathtt{srot(:,1:3)}$ は、積算の丸め誤差のため直交からわずかにずれている。`spnorm` (`photons.f:767-839`) は次の順で処理する。

1. 列 2 を規格化、列 1 から列 2 の成分を引いて規格化、列 3 は外積 (直交化)。
2. 回転の軸の初期値: 反対称部分 $\mathtt{sps(:,1)} = (R_{32} - R_{23}, R_{13} - R_{31}, R_{21} - R_{12})$ を規格化 ($R = \exp(\mu\,[\hat n]_\times)$ のとき $\propto \sin\mu\,\hat n$)。大きさが $10^{-4}$ 未満なら、実固有値の固有ベクトル (`teigen`)。
3. $(R - 1)\,\delta\vec n = \vec n - R\vec n$ の最小ノルム解 (`tsolvg`) で軸を補正し、規格化。成分の符号が負側に偏っていれば反転。これが $n_0$ (`PolarizationVector`)。
4. $n_0$ に垂直な基底 $\vec e_2, \vec e_3$ を作る。$R\vec e_2$ の成分から $c_m = \vec e_2\cdot R\vec e_2$、$s_m = \vec e_3\cdot R\vec e_2$ を出し、$\mu_s = \mathrm{atan2}(-s_m, c_m)$、$\nu_s = \mu_s/2\pi$ (`SpinTune`)。
5. $\mathtt{sdamp} = n_0\cdot\mathtt{gintd}$ (Sokolov-Ternov の 1 周の率、§7)。

**確認** (e): 平面リング (3 GeV、12 セル、偏向 15 度 $\times$ 12) で $n_0 = (0, 1, 0)$、$\nu_s = 0.40407411$ ($a\gamma\theta_{\rm tot}/2\pi = 3.404074$ の端数 $0.404074$ と一致。放射損失と RF による $\gamma$ の補正は、この桁では見えない)。

## 6. SLIM: 軌道が変わると軸がどれだけ傾くか

### 6.1 $N$ の式

粒子の座標 $x_0$ (リング入口) に対するスピンの周期解の軸を $n(x_0)$ とする。軸は、1 周の回転と 1 周の写像 $x_1 = Mx_0$ に対して $R(x_0)\,n(x_0) = n(Mx_0)$ を満たす。1 次まで $n(x_0) = n_0 + Nx_0$ ($N$: 3×6、$n_0\cdot N = 0$) とおく。(15.12) より $R(x_0)\,n(x_0) = \exp(-[Dx_0]_\times)(R_0n_0 + R_0Nx_0) \approx n_0 + R_0Nx_0 - (Dx_0)\times n_0 = n_0 + R_0Nx_0 + n_0\times(Dx_0)$ ($R_0n_0 = n_0$)。右辺は $n_0 + NMx_0$。$x_0$ の係数を比べて

$$
N M = R_0\,N + [n_0]_\times D\qquad\Longleftrightarrow\qquad N = \left(R_0N + [n_0]_\times D\right)M^{-1}
\tag{15.13}
$$

($[n_0]_\times$ は外積の行列。) (15.13) は $N$ の線形方程式で、解は (チューンが $\nu_s\pm\nu_k$ = 整数 でない限り) 一意である ($\mathrm{vec}$ を取ると $(1 - M^{-\top}\otimes R_0)\mathrm{vec}(N) = \mathrm{vec}([n_0]_\times DM^{-1})$)。分母 $1 - e^{i(\mu_s\pm\mu_k)}$ が 0 に近づくのが、線形のスピン軌道共鳴 (チューン $\nu_s\pm\nu_x$、$\nu_s\pm\nu_z$ など)。要素ごとの形は $N_{i+1}T_i = R_iN_i + [n_{i+1}]_\times D_i$ (要素 $i$ の行列 $T_i$、回転 $R_i$、結合 $D_i$ (出口の座標系)。$N_{i+1}$ は出口の軸の座標依存)。

**確認** (d): $R_0$ を $n_0 = \hat y$ のまわりの回転 ($\nu_s = 0.403$)、$M$ を減衰 0.9 の 3 つの回転、$D$ を乱数にした玩具の 1 周写像で、(15.13) を解き ($n_0\cdot N = 2\times10^{-16}$)、$|R(x_0)\,n(x_0) - n(Mx_0)|$ が $|x_0| = 10^{-2}, 10^{-3}, 10^{-4}$ で $1.3\times10^{-5}$、$1.2\times10^{-7}$、$4.0\times10^{-10}$ ($|x_0|^2$ に比例、1 次の式として成り立つ)。

### 6.2 脱偏極率

放射が起こらなければ、スピン $\vec S$ は軸 $n(x)$ に追従し (断熱不変量 $\vec S\cdot n(x)$ が保存)、偏極は保たれる。1 個の光子の放射は、粒子の座標を $x\to x+\xi$ ($\xi$: 跳び、平均 0、共分散 $B_i$、`tradke` の `dbeam`) と変えるが、スピンは変えない (光子が運ぶスピンの寄与は Sokolov-Ternov 項として別に扱う)。軸は $n(x)\to n(x+\xi) = n(x) + N\xi$ になるので、

$$
\vec S\cdot n(x+\xi) = 1 - \tfrac12|N_\perp\xi|^2 + O(\xi^3)
\tag{15.14}
$$

($|n| = 1$ から $n\cdot N\xi = -\tfrac12|N\xi|^2$ を使った。) したがって、1 周あたりの偏極の減りの率 (脱偏極率) は

$$
\frac{1}{\tau_{\rm dep}\,f_{\rm rev}} = \frac12\sum_i\mathrm{tr}\left(N_iB_iN_i^\top\right)
\tag{15.15}
$$

(和は 1 周の放射のスライス $i$ について。$N_i$: そのスライスの出口での $N$。) $B_i$ のエネルギーの跳びの成分は $B_{\delta\delta} = N_\gamma\,u_c^2\,c_{\rm uu}$ ([13_radiation.md](13_radiation.md) の `dee`)。これが §7.3 で DK の係数 $11/18$ に対応する。

### 6.3 平衡偏極度

Sokolov-Ternov の偏極 (立ち上がりの率 $\tau_{ST}^{-1}$、到達値 $P_{ST}$) と脱偏極 (率 $\tau_{\rm dep}^{-1}$) の釣り合い $dP/dt = (P_{ST} - P)/\tau_{ST} - P/\tau_{\rm dep}$ の定常解は

$$
P_{\rm eq} = P_{ST}\,\frac{\tau_{\rm dep}}{\tau_{ST} + \tau_{\rm dep}} = P_{ST}\,\frac{\mathtt{sdamp}}{\mathtt{sdamp} + \mathtt{rate}}
\tag{15.16}
$$

(`rate` は (15.15) の 1 周あたりの率、`sdamp` は 1 周あたりの Sokolov-Ternov の率。) v2 の `Polarization` はこの形 (`slim_depolarization`)。K64 の `srequpol2` (`photons.f:861-907`) も同じ形の式 $P = \mathtt{pst}\,\mathtt{sdamp}/(\mathtt{rmd} + \mathtt{sdamp})$ を持つが、**どこからも呼ばれない** (呼ばれるのは `srequpol` だけ)。`rmd` は、正規座標に直した $D$ の $n_0$ に垂直な成分 `drot(2,:)`、`drot(3,:)` から作る SLIM の和で、分母 $|e^{-a_k + i(\mu_s + \mu_k)} - 1|^2$ が (15.13) の共鳴の分母に当たる。

**K64 の `srequpol`** は、`spdepol` (`photons.f:620`) で正規座標の 2 次の摂動展開 (6 次までの級数、`spintrm`) から `rm1` を作り、`equpol` $= \mathtt{pst}\,\mathtt{sdamp}/(\mathtt{sdamp} - \mathtt{rm1}(1,1))$ とする。`spdepol` の中の `exp(cb²−cφ)·cerfc(cb)` は、$|cb| > 18$ の漸近式で 0·∞ になり、`Polarization` は NaN (B55)。この級数の導出は確認していない。v2 はこの式を使わず、(15.15)〜(15.16) と側帯の項で置き換えた ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §5.3、B125)。

## 7. Sokolov-Ternov、定数の導出

### 7.1 QED の入力

電子が磁場の中で放射するとき、スピンの向きによって放射の確率がわずかに違う (スピンフリップ)。これを積分した結果 (QED、ここでは導かない): 一様な曲率半径 $\rho$ の磁場で、スピンが磁場の向きに沿う ($\uparrow\downarrow$) 偏極の立ち上がりは

$$
\frac{dP}{dt} = \frac{P_{ST} - P}{\tau_{ST}},\qquad
P_{ST} = \frac{8}{5\sqrt3} = 0.92376,\qquad
\frac1{\tau_{ST}} = \frac{5\sqrt3}{8}\,\frac{r_e\bar\lambda_C\,c\,\gamma^5}{|\rho|^3}
\tag{15.17}
$$

($\bar\lambda_C = \hbar/mc$)。磁場の向きが変わるリングでは $|\rho|^{-3}$ を $\langle|\rho|^{-3}\rangle$ (経路長の平均)、$P_{ST}$ の方向は $\hat b$ との内積で効く (§7.4)。

### 7.2 K64 の定数が (15.17) と一致すること

K64 は $P_{ST}$、$c_{\rm fpd}$、$s_{\rm flc}$ を

$$
\mathtt{pst} = \frac{8\sqrt3}{15} = \frac{8}{5\sqrt3},\qquad \mathtt{cfpd} = \frac{5\sqrt3}{8},\qquad \mathtt{sflc} = \frac34\left(\frac{r_e}{\alpha}\right)^2 = \frac34\bar\lambda_C^2
\tag{15.18}
$$

と持つ ($r_e/\alpha = \bar\lambda_C$)。`sprot` の項は、$\mathtt{ucde} = u_c\,\delta E = c_{\rm ave}\,N_\gamma\,u_c^2$ ($N_\gamma = a_{\rm nrad}\gamma\theta$、$u_c = c_{\rm uc}\gamma^3\theta/(p_0\,al)$、13 章) を使って $\mathtt{dst} = (\mathtt{pst} - s_t)\,c_{\rm fpd}\,(p\,p_0/\gamma^2)^2\,\mathtt{ucde}/|b_t|$。$\beta\to1$ ($pp_0/\gamma^2 \to \beta^2$) で、1 スライスの偏極の変化 ($s_t = 0$ のとき、$|b_t| = \theta$) は

$$
\Delta P = \mathtt{pst}\,\mathtt{cfpd}\,\beta^2\,c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2\;\frac{\gamma^5\theta^3}{al^2}\ \stackrel{\theta = al/\rho}{=}\ \mathtt{pst}\,\mathtt{cfpd}\,\beta^2\,(c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2)\,\frac{\gamma^5\,al}{\rho^3}
\tag{15.19}
$$

(15.17) の $dP = P_{ST}\,\tau_{ST}^{-1}\,dt$、$c\,dt = al/\beta$ は $\mathtt{pst}\,\mathtt{cfpd}\,r_e\bar\lambda_C\,\gamma^5 al/(\beta\rho^3)$。$\beta^3\to1$ で一致するためには、

$$
c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2 = \frac{8}{15\sqrt3}\cdot\frac{5\alpha}{2\sqrt3}\cdot\frac94\bar\lambda_C^2 = \alpha\bar\lambda_C^2 = r_e\bar\lambda_C
\tag{15.20}
$$

($\tfrac{8\cdot5}{15\cdot2\cdot3}\cdot\tfrac94 = \tfrac49\cdot\tfrac94 = 1$。$c_{\rm uc} = \tfrac32\bar\lambda_C$、$a_{\rm nrad} = \tfrac{5}{2\sqrt3}\alpha$。) **確認** (c): 比 $1.000000000000$。つまり `sprot` の 1 ステップの項は ST の式そのものである ($P_{ST}\,\mathtt{cfpd} = 1$)。

**Emittance[] 側の `gintd`**: `tradke` の `gintd += sw(1:3)*sflc*anp*(bt*h1*p/al1)**2` は、1 周の積算 ($\mathtt{sw}$ は $\hat g$ を (srot の列で) リング入口まで戻した成分、$\mathtt{sdamp} = n_0\cdot\mathtt{gintd}$ は (各スライスの $\hat g\cdot n$) の和)。係数 $s_{\rm flc}\,a_{\rm nrad} = \tfrac34\bar\lambda_C^2\cdot\tfrac{5\alpha}{2\sqrt3} = \tfrac{5\sqrt3}{8}\,r_e\bar\lambda_C = c_{\rm fpd}\,r_e\bar\lambda_C$ ($\tfrac{15}{8\sqrt3} = \tfrac{5\sqrt3}{8}$) で、$\mathtt{gintd}\cdot\hat n$ はスライスごとに $c_{\rm fpd}\,r_e\bar\lambda_C\,\gamma^5\beta^2\,\theta^3/al^2\cdot(\hat g\cdot\hat n)$ を足す。これは (15.17) の $\tau_{ST}^{-1}$ の経路長の積分 ($\hat g\cdot\hat n = 1$ のとき) に当たる。**確認** (c): $s_{\rm flc}a_{\rm nrad}/(c_{\rm fpd}r_e\bar\lambda_C) = 1.000000000000$。

**偏極時間**: $\tau = 1/(\mathtt{sdamp}\,f_{\rm rev})$ ($f_{\rm rev}$: 周回周波数。$\mathtt{sdamp}$ は 1 周あたり、$f_{\rm rev}$ は 1/s)。**確認** (e): 平面リング ($\rho = 4.5837$ m、偏向の長さ $14.4$ m、周長 $45.6$ m、$\langle\rho^{-3}\rangle = 3.2455\times10^{-3}$ m$^{-3}$、$\gamma = 5870.854$) で (15.17) の $\tau_{ST} = 123.8153$ s、Emittance[] の `PolarizationTime` $= 123.8153$ s (比 $1.000000$)。平衡偏極度 $P = 0.92376043$ (v2、$8/(5\sqrt3) = 0.923760$)。

### 7.3 SLIM の $\delta$ の項と DK の $11/18$

エネルギーの跳びの共分散 $B_{\delta\delta} = N_\gamma\,u_c^2\,c_{\rm uu}$ ($c_{\rm uu} = 11/27$。$N_\gamma = a_{\rm nrad}\gamma\theta$、$u_c = c_{\rm uc}\gamma^3\theta/(p_0\,al)$ から $B_{\delta\delta} = c_{\rm uu}a_{\rm nrad}c_{\rm uc}^2\gamma^5\theta^3/(\beta^2al^2)$) を (15.15) に入れる。$|\partial n/\partial\delta|^2$ の係数 $\tfrac12B_{\delta\delta}$ を、§7.2 の ST の 1 スライスの係数 $s_{\rm flc}\,a_{\rm nrad}\,\beta^2\gamma^5\theta^3/al^2$ ($\mathtt{gintd}$ の 1 項) で割ると ($\beta\to1$)

$$
\frac12\,\frac{c_{\rm uu}\,a_{\rm nrad}c_{\rm uc}^2}{s_{\rm flc}\,a_{\rm nrad}} = \frac12\cdot\frac{\dfrac{55}{24\sqrt3}}{\dfrac{15}{8\sqrt3}} = \frac{11}{18}
\tag{15.21}
$$

(確認 (c): 数値で $0.611111111 = 11/18$。式の形: $c_{\rm uu}a_{\rm nrad}c_{\rm uc}^2 = \tfrac{11}{27}\cdot\tfrac{5\alpha}{2\sqrt3}\cdot\tfrac94\bar\lambda_C^2 = \tfrac{55}{24\sqrt3}\alpha\bar\lambda_C^2$、$s_{\rm flc}a_{\rm nrad} = \tfrac{5\sqrt3}{8}\alpha\bar\lambda_C^2 = \tfrac{15}{8\sqrt3}\alpha\bar\lambda_C^2$。比 $\tfrac12\cdot\tfrac{55}{24\sqrt3}\cdot\tfrac{8\sqrt3}{15} = \tfrac{55\cdot8}{2\cdot24\cdot15} = \tfrac{11}{18}$。) したがって、$\delta$ の跳びによる脱偏極率は $\tfrac{5\sqrt3}{8}r_e\bar\lambda_Cc\gamma^5\langle|\rho|^{-3}\cdot\tfrac{11}{18}|\partial n/\partial\delta|^2\rangle$ の形になり、DK の式 (§7.4) の $\tfrac{11}{18}|\gamma\,\partial n/\partial\gamma|^2$ と同じ係数を持つ ($\beta\to1$ では $\gamma\,\partial/\partial\gamma = \partial/\partial\delta$)。SLIM の式 (15.15) は、DK の式の $\gamma\,\partial n/\partial\gamma$ の項を、全座標の $N$ に拡張した形である (水平・垂直の寄与は $B_{xx}$ などが入る)。

### 7.4 Derbenev-Kondratenko の式との対応

標準の DK の式 (導出は QED の摂動、ここでは引用):

$$
P_{\rm DK} = \frac{8}{5\sqrt3}\,\frac{\left\langle|\rho|^{-3}\,\hat b\cdot\left(\hat n - \gamma\,\dfrac{\partial\hat n}{\partial\gamma}\right)\right\rangle}{\left\langle|\rho|^{-3}\left[1 - \dfrac29(\hat n\cdot\hat v)^2 + \dfrac{11}{18}\left|\gamma\dfrac{\partial\hat n}{\partial\gamma}\right|^2\right]\right\rangle}
\tag{15.22}
$$

($\hat b$: 磁場の向き (進行方向と曲率の向き)、$\hat v$: 速度の向き、$\langle\rangle$: 経路長と位相空間の平均、$\gamma\,\partial\hat n/\partial\gamma$ はスピンチューン $\nu = a\gamma$ の変化に対する軸の応答。) 分子の第 2 項は運動論的偏極、分母の第 2 項が $\tfrac29$、第 3 項が (15.21) の $11/18$。(15.16) との対応: $\mathtt{sdamp}\leftrightarrow\langle|\rho|^{-3}[1 - \tfrac29(n\cdot v)^2]\rangle$ (v2 の `sdamp_rel`、`ST_EXACT`)、$\mathtt{rate}\leftrightarrow\tfrac{11}{18}$ の項 (SLIM の全座標への拡張)、運動論的偏極は v2 の `KINETIC`。K64 の `sdamp` は符号つきの $\langle|\rho|^{-3}\hat b\cdot\hat n\rangle$ ($\hat g\cdot n$ の和、逆向きの偏向で相殺する、B116) で、$(1-\tfrac29(n\cdot v)^2)$ と運動論の項は入らない。

## 8. 脱偏極の現状 (要約。詳細は docs/16)

- 側帯の項: スピンチューンが $a\gamma(1+\delta)$ のように $\delta$ に依存するため、シンクロトロン振動でスピンの位相が変調される。変調指数 $a\gamma\sigma_\delta/\nu_s\gtrsim0.3$ で線形の SLIM (15.15) は不十分になる。v2 は位相の変調を厳密に取る (`sideband_depolarization`)。
- 線形の軌道 + 厳密なスピン回転のモンテカルロとは $-10\%$ 以内で合う。FLUC つきのスピン追跡とは、4.2 GeV・12 MV の FODO で側帯の式が 0.43 倍、クロマティックな変調を入れて $1.05\pm0.16$、3 GeV 以下で 0.5〜1.6 倍 (統計 ±25%)。未分離の非線形 (色収差・運動学的な項・RF の非線形) が残っている。
- 逆向きの偏向を持つ挿入光源 (MPW) では、K64 の `sdamp` は符号つきで相殺するため $P = 92.376\%$ 一定となり不適切 (B116)。v2 は相殺しない緩和率を使う。

## 9. Fortran の対応

| 式 | ファイル:行 | サブルーチン | 変数・内容 |
|---|---|---|---|
| (15.1)〜(15.5) | `photons.f:1316-1319` | `sprot` | `ct = 1+h*gspin` (横)、`cl = 1+gspin` (縦) |
| (15.7)、(15.8) | `photons.f:1179-1181`、`1305-1316` | `tradk1`、`sprot` | `ppx, ppy, ppz`、`bx, by, bz`、`bp`、`blx..`、`btx..` |
| (15.9) | `photons.f:1317-1319` | `sprot` | `gx, gy, gz` |
| (15.10) | `photons.f:1345-1356` | `sprot` | `xsincos`、`sw`、`sux..` |
| (15.11) | `photons.f:1357-1358` | `sprot` | `sx = sx*cphi0+sz*sphi0`、`sz = (sz-sx*sphi0)/cphi0` |
| 追跡の呼び出し | `photons.f:1204-1208`、`1258-1262` | `tradk1`、`tradkn` | `pxm = px0+dpx/2`、`ucde = uc*de` |
| (15.12) | `photons.f:1521-1589`、`1614-1615` | `tradke` | `dgx..dgz` ($\partial g/\partial x$)、`transr`、`srot(:,4:9)` |
| $\partial b_t/\partial x$ | `photons.f:1546-1551` | `tradke` | `dbtx..dbtz`、`dxpx..dxpz`、`dh1r`、`dct` |
| $R$ の積算 | `photons.f:1590-1615` | `tradke` | `sinu, cosu, dcosu`、`sx..sz` ($srot$ の行) |
| `gintd` | `photons.f:1604-1605` | `tradke` | `gintd += sw(1:3)*sflc*anp*(bt*h1*p/al1)**2` |
| 偏向角の座標回転 | `photons.f:1614-1615` | `tradke` | `cp`、`sp` ($\phi_{r0}$ の半角から) |
| 初期化 | `photons.f:980-990` | `srotinit` | `srot = [I | 0]` |
| 1 周の処理 | `photons.f:767-839` | `spnorm` | `sps`、`smu`、`sdamp` |
| 平衡偏極度 | `photons.f:909-978`、`620`〜 | `srequpol`、`spdepol`、`spintrm` | `rm1`、`equpol` (NaN: B55)。`srequpol2` (`:861`) は未使用 |
| 呼び出し | `temit.f:493-496` | `tecalc` の末尾 | `call spnorm(srot,sps,spinmu,sdamp)`、`call srequpol(...)` |
| 定数 | `photons.f:379-382`、`tfstk.f:269` | module `tspin` | `pst`、`sflc`、`cfpd`、`cuu`、`cl`、`gspin` |

(行番号は参照時点のもの。)

## 10. pysad の対応

| 関数 | 対応する Fortran | 注意 |
|---|---|---|
| `track6._sprot` (`track6.py:134`) | `sprot` | 粒子ごとのスピン (`_Particle.s`)。K64 は弦 $2\sin(\theta/2)$。v2 は $\theta$ ($\arcsin$ で直す、B83)。`p.bpar` (SOL・UND) があるときは縦成分を運動量の向きに分けて足す (v2) |
| `spin.tradke_spin` (`spin.py:70`) | `tradke` のスピン部 | $\partial g/\partial x$ は解析式。添字は 0 始まり。`srot` は `st` (3×9) |
| `spin.tchge_rotate_spin`、`tbrote_spin`、`tchge_spin` | `tchge`、`tbrote` | ずれ・ROTATE・CHI2 の座標の回転。K64 は行を逐次更新 (B56)。v2 は行列積 |
| `spin.spnorm` (`spin.py:395`) | `spnorm` | `PolarizationVector`、`SpinTune`、`sdamp` |
| `spin.srequpol` (`:592`) | `srequpol` | k64: NaN。v2: `slim_depolarization` (`:793`) |
| `spin.sideband_depolarization` (`:623`) | (K64 にない) | 側帯の項。準モンテカルロ (Sobol 列) |
| `spin.spin_results` (`:897`) | `tecalc` の末尾 | `equpol`、`depol_rate`、`taudep`、`spin_orbit` (入口の $N$) |

SAD スクリプトでは、`POL` を有効にして `Emittance[]` の `PolarizationTime`、`SpinTune`、`NominalSpinTune`、`Polarization`、`PolarizationVector`、`SpinPrecessionMatrix` を見る。追跡は `TrackParticles` を 9〜11 列 (7 列目 = 歳差の角、8 列目 = $s_y/|s|$) で呼ぶ ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §2)。

## 11. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `spin_01_thomas_bmt.py` | (a) 静止系の歳差 + Thomas 歳差 = BMT (sympy)、縦・横の係数、(b) らせん軌道の BMT の厳密な積分 (DOP853) と 1 ステップ回転の積 (2 次)、(b2) 平面の偏向で $a\gamma\theta$、縦磁場で $(1+a)\psi$、K64 の弦 $-\theta^2/24$、(c) 定数 $c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2 = r_e\bar\lambda_C$、$s_{\rm flc}a_{\rm nrad} = c_{\rm fpd}r_e\bar\lambda_C$、DK の 11/18、(d) (15.13) が 1 次で $R(x)n(x) = n(Mx)$、(e) 平面リングの `SpinTune`・`PolarizationTime`・`Polarization` | [spin_01_thomas_bmt.out.txt](outputs/spin_01_thomas_bmt.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/spin_01_thomas_bmt.py` (1 分弱。(e) の Emittance[] が大部分)。K64 との照合 (gs と 4e-12 以下、Emittance[] の偏極時間 14 桁など)、DK の数値積分との比較 (3 GeV〜4.2 GeV、追跡 FLUC) は、[docs/16_spin_pol.md](../docs/16_spin_pol.md) と `pysad/tests/script_tests/test_spin*.py`・`test_v2_*spin*.py`・`_p1_dk.py` に既にある。この章では新しく取り直していない。

## 12. 未確認・限界・使った条件

**使った条件**: 電子 (a = 1.159652181280002e-3)、3 GeV (γ = 5870.85)、平面リング (12 セル、偏向 15 度 × 12 = π。周長 45.6 m、$\rho = 4.5837$ m、CAVI 1.5 MV・harmonic 64、`RFSW; RAD; NOFLUC; POL`)。(b) のらせん軌道は計算上の一様磁場 ($\kappa = 1$ rad/m、角 0.3 rad)。

**未確認**:
- Sokolov-Ternov の遷移率 (15.17) と DK の式 (15.22) の QED からの導出 (入力として使った。定数は K64 と式の間の一致だけ確認)。
- (15.12) の高次 ($\partial^2\vec g/\partial x^2$)、$|g|$ が大きい (非可換の補正が効く) 場合。K64 の $D$ は 1 次だけ。
- `spdepol`・`spintrm` の 6 次の級数の導出 (K64 は NaN で使えない)。v2 の式は独立な DK の数値積分・モンテカルロと比較してある範囲 (docs/16) でのみ確認済み。
- 運動論的偏極・高次のスピン軌道共鳴 ($N$ の分母が 0 に近い) の扱い。分母の 0 の近くは (15.13) が使えない。
- 領域内 (SOL) の縦磁場、傾いた磁場での $B_\parallel$ の扱い (v2 の `bpar`)。
- 本章の数値は 1 つの平面リングで確認。結合のあるリング・非平面の偏極ベクトルは docs/16 の試験に頼る。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| POL・RADPOL・TrackParticles・Emittance[] の説明 | **原文で確認**: 「POL が ON のとき TrackParticles はスピン追跡を行う。RADPOL は Sokolov-Ternov を入れる」。座標は `{x, px/p0, y, py/p0, z, dp/p0, sarg, sy}`。Emittance[] のキー: `Polarization`、`Polarization2`・`4`・`6`、`PolarizationVector`、`SpinTune`、`NominalSpinTune`。 |
| `Polarization2`・`4`・`6` | ヘルプにだけある (2・4・6 次)。K64 のソースに文字列がなく、gs の `Emittance[]` も返さない (docs/16)。 |
| ヘルプに無いもの | Thomas-BMT、`sprot` の式、$D$、$N$、Sokolov-Ternov の定数、`spnorm` の手順、`Polarization` の式。(15.1)〜(15.22) はすべて K64 のソースと標準の理論から。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

**確定 (独立な基準で確認済み。docs/16、docs/30)**:
- **`Polarization` が NaN** (`spdepol`、B55): 桁あふれ。さらに桁あふれだけ直した値は、追跡・独立な式に対して約 1000 倍の脱偏極率になる (B125。二重数え)。
- **符号つきの `sdamp`** (B116): 逆向きの偏向で相殺する。
- **弦 $2\sin(\theta/2)$** (B83): 刻みが粗いと歳差の角が $\theta^2/24$ だけ小さい (本章確認 (b2))。
- `tbrote`・`serot` の行の逐次更新 (B56)、`bsi2` の項の脱落 (B58): docs/16。

**疑い**:
- **確定 (K64 の誤り、v2 は修正 = B336、gs では確認できない)**: `tradke` の `dbtz = dxpzb/pr - dblz` (`photons.f:1550`): $x$ 成分・$y$ 成分は `dxpx/pr**2` だが、$z$ 成分だけ `pr` で割っている。自分の `btz = xpzb/pr**2 - blz` の微分として不整合 (誤植)。$\delta \ne 0$ の軌道で Emittance[] (POL) のスピン軌道結合 $D = \mathtt{srot}(:,4{:}9)$ の $z$ の行 (縦の歳差) が $\delta$ に比例して誤る。追跡のスピン回転の中心差分から作った $D$ との最大値の相対差は、3 GeV/c、15° の矩形の偏向磁石 ($K_1 = 0$)・軌道 $(x, p_x, y, p_y) = (8, 1, 6, -1)\times10^{-3}$ の周りで、$\delta$ = 0.05 に $\mathrm{pr}$ で割る式 13.2%、$\mathrm{pr}^2$ で割る式 1.20%、$\delta$ = 0.10 に 26.9% 対 2.24% ($K_1$ = 0.3 は 6.3% 対 0.18%、12.4% 対 0.15%)。$\delta$ = 0 は同じ。四極は 0.2% で $\delta$ に依らない。TrackParticles のスピン追跡には影響しない (別のコード)。gs では確認できない: K64 の Polarization は NaN (B55)、`srot` の微分列は Emittance[] の他の出力に出ない。限界: 1 要素 + ドリフトの比較 ($\delta$ = 0.1 の $K_1 = 0$ の偏向磁石の v2 の 2.2% の残りは別の項、未調査)。試験 `test_v250_spin_dbtz.py`。
- `tradke` の `call tadddsg(dsg,srot,dgx,dgy,dgx,dbeam)`: `dgz` のところに `dgx` を渡している (B57)。`dsg` は使われないので影響なし。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `gspin` | $a = (g-2)/2$ | `tfstk.f:269` | `GSPIN` |
| `cl` | $1 + a$ | `photons.f:382` | `CL` |
| `ct` | $1 + a\gamma$ (`1+h*gspin`) | `sprot`、`tradke` | |
| `pst` | $8\sqrt3/15$ | `photons.f:379` | `PST` |
| `cfpd` | $5\sqrt3/8$ | `photons.f:380` | `CFPD` |
| `sflc` | $0.75\,(r_e/\alpha)^2$ | `photons.f:380` | `sflc()` |
| `ucde` | $u_c\,\delta E$ (`uc*de`) | `tradk1` | |
| `ppx, ppy, ppz` | $\vec p_{\rm out}\times\vec p_{\rm in}$ | `tradk1` | |
| `pxm, pym, pzm` | 中点の運動量 | `sprot` | |
| `bx0..bz0`、`bx..bz` | 外積、中点の $p_z$ で割ったもの | `sprot` | |
| `bsi`、`bsir0` | 縦磁場の積分 | `sprot`、`temw` | `bsi` |
| `bp`、`blx..blz`、`btx..btz` | 縦・横の磁場成分 | `sprot`、`tradke` | |
| `gx, gy, gz`、`g` | 回転ベクトル、$|g|$ | `sprot`、`tradke` | |
| `xsincos` → `sinu, xsinu, cosu, dcosu` | $\sin|g|$、$x\cos x - \sin x$ の関数、$\cos|g|$、$\cos|g| - 1$ | `mathfun` | `xsincos` |
| `sw`、`sux..suz` | $(S\cdot g)/|g|$ (または $(1-\cos)$ の項)、$S\times g$ | `sprot`、`tradke` | |
| `cphi0, sphi0`、`cp, sp` | スライスの偏向角の余弦・正弦 | `tspin`、`tradke` | |
| `srot(3,9)` | $R$ (1:3 列) と $D$ (4:9 列) | `temw`、`tradke` | `st` |
| `transr` | $M_{0\to i}$ (前のスライスまでの累積行列) | `temw` | |
| `dgx, dgy, dgz`、`dct`、`dh1r` | $\partial g/\partial x_{\rm in}$ (6 成分)、$\gamma$ の $\delta$ 微分 | `tradke` | |
| `gintd(3)` | Sokolov-Ternov の傾きの積算 | `temw` | `gintd` |
| `sps(3,3)` | 偏極ベクトルと垂直な基底 | `spnorm` | |
| `smu`、`sdamp` | $2\pi\nu_s$、1 周の ST の率 | `spnorm` | `SpinTune`、 |
| `drot(3,6)` | `spsᵀ srot(:,4:9) r` (正規座標の $D$ を $n_0$ の基底で) | `srequpol` | |
| `rmd`、`rm1`、`equpol` | 脱偏極の和、脱偏極行列、平衡偏極度 | `srequpol2`、`srequpol` | `equpol` |
| `calpol`、`radpol` | `POL`、`RADPOL` のフラグ | `ffs_flag` | |
