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
- **SLIM**: スピン軸の座標依存 $n(x) = n_0 + Nx$ は $N M = R_0 N + [n_0]_\times D$ を満たし (導出 §6.1、確認 (d): 1 次の代数の整合)、放射の跳び $\xi$ ($\langle\xi\xi^\top\rangle = B_i$) による偏極の減りは 1 回あたり $\tfrac12|N\xi|^2$、脱偏極率 $= \tfrac12\sum_i\mathrm{tr}(N_iB_iN_i^\top)$ (1 周あたり)。この率は、確率的な軌道とスピンの玩具のモンテカルロで 5% 以内に確かめた ([spin_03_slim_mc.out.txt](outputs/spin_03_slim_mc.out.txt)、§6.2)。正規座標では $N$ の分母が $|e^{-a_k+i(\mu_s\pm\mu_k)}-1|$ になる (§6.4)。
- **Sokolov-Ternov**: K64 の定数は $c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2 = r_e\bar\lambda_C$、$s_{\rm flc}\,a_{\rm nrad} = c_{\rm fpd}\,r_e\bar\lambda_C$ ($c_{\rm fpd} = 5\sqrt3/8$) で、相対 $10^{-12}$ で一致する (確認 (c))。偏極時間 $1/\tau_{ST} = \tfrac{5\sqrt3}{8}r_e\bar\lambda_Cc\gamma^5\langle\rho^{-3}\rangle$ は、平面リングで `PolarizationTime` $= 123.8153$ s と 6 桁以上一致。平衡偏極度は $8/(5\sqrt3) = 0.92376$ (v2 の `Polarization` $= 0.923760431$)。
- **DK との対応**: SLIM の $\delta$ の寄与 $\tfrac12 B_{\delta\delta}$ の係数 $\tfrac12 c_{\rm uu}\,a_{\rm nrad}c_{\rm uc}^2/(s_{\rm flc}\,a_{\rm nrad}) = 11/18$ は、DK の式の $\tfrac{11}{18}|\gamma\,\partial n/\partial\gamma|^2$ の係数に一致する (確認 (c))。
- **K64 の限界**: `Polarization` は k64 で NaN (`spdepol` の `exp(cb²−cφ)·cerfc(cb)` の桁あふれ、B55。機構は §6.5)。v2 は SLIM + 側帯の式で置き換えた。高エネルギー・大きな RF 電圧では v2 の脱偏極率が追跡 (FLUC) の 0.4 倍に過小評価される例があり ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §5.3)、確認は追跡で行う。

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

静止系 (記号に $'$) では、磁気モーメント $\vec\mu = \dfrac{gq}{2m}\vec S$ が磁場 $\vec B'$ から受けるトルクで $\dfrac{d\vec S}{dt'} = \vec\mu\times\vec B' = -\dfrac{gq}{2m}\vec B'\times\vec S$。回転角速度 $\vec\Omega' = -\dfrac{gq}{2m}\vec B'$ である ($g/2 = 1+a$ なので $\vec\Omega' = -(1+a)\dfrac qm\vec B'$)。静止系には電場 $\vec E'$ も現れるが、静止した磁気モーメントは $\vec E'$ と結合しない (1 次のスピン方程式に入らない) ので落とす。

**静止系の磁場** (電場のない実験室系の磁場 $\vec B$ から)。電磁場テンソルの Lorentz 変換 (4×4 の速度 $\vec\beta$ のブースト $\Lambda$ による $F' = \Lambda F\Lambda^\top$) を、$\vec E = 0$ の $F$ について sympy で計算し、$\vec B'$ の 3 成分を取り出すと

$$
\vec B' = \gamma\left[\vec B - \frac{\gamma}{\gamma+1}\vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.1}
$$

に一致する ([spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (a)。差は 3 成分とも 0)。手で確かめるには、場の変換則 $B'_\parallel = B_\parallel$、$\vec B'_\perp = \gamma(\vec B - \vec\beta\times\vec E)_\perp = \gamma\vec B_\perp$ ($\vec E = 0$) から出発する。$\vec B_\parallel = \hat\beta(\hat\beta\cdot\vec B)$、$\vec B_\perp = \vec B - \vec B_\parallel$ なので

$$
\vec B' = \vec B_\parallel + \gamma(\vec B - \vec B_\parallel) = \gamma\vec B - (\gamma-1)\,\frac{\vec\beta(\vec\beta\cdot\vec B)}{\beta^2}
$$

で、$\gamma^2\beta^2 = \gamma^2 - 1 = (\gamma-1)(\gamma+1)$ から $(\gamma-1)/\beta^2 = \gamma^2/(\gamma+1)$。これを代入すると (15.1) になる。縦成分は $B'_\parallel = \gamma(1 - \gamma\beta^2/(\gamma+1))B_\parallel = B_\parallel$ ($\gamma\beta^2/(\gamma+1) = (\gamma-1)/\gamma$ を使う)、横成分は $\gamma B_\perp$ (sympy の確認 (a) の出力の 2 行目)。

**時間の変換**。静止系の固有時間 $t'$ と実験室系の時間 $t$ は $dt = \gamma\,dt'$。$\dfrac{d\vec S}{dt} = \dfrac1\gamma\dfrac{d\vec S}{dt'} = \dfrac1\gamma\vec\Omega'\times\vec S$ なので、実験室系の時間で見た角速度は $\vec\Omega_{\rm rest} = \vec\Omega'/\gamma$。(15.1) を代入すると、

$$
\vec\Omega_{\rm rest} = -\frac{(1+a)q}{m}\left[\vec B - \frac{\gamma}{\gamma+1}\vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.2}
$$

($\gamma$ が打ち消し合う。)

### 3.2 Thomas 歳差

加速される粒子の静止系の軸は、実験室系に対して、各瞬間の Lorentz ブーストの合成のために回転する。

**回転の大きさ**。時刻 $t$ の速度 $\vec v$ のブースト $\Lambda(\vec v)$ と、わずかに違う速度 $\vec v + \delta\vec v$ のブースト $\Lambda(\vec v+\delta\vec v)$ をつなぐ変換 $\Lambda(\vec v+\delta\vec v)\Lambda(\vec v)^{-1} = I + (\text{ブースト生成子}) + (\text{回転生成子})$ を 1 次で展開する。回転生成子は、この 4×4 行列の空間 3×3 部分の反対称部分である。ここで $\vec v = (\beta, 0, 0)$、$\delta\vec v = (\delta v_x, \delta v_y, \delta v_z)$ とすると、sympy ([spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (b)) で、反対称部分は回転ベクトル

$$
\vec w = \frac{\gamma^2}{\gamma+1}\,\delta\vec v\times\vec v
\tag{15.3a}
$$

に等しい。成分で書くと、$\delta\vec v\times\vec v = (0,\ \beta\,\delta v_z,\ -\beta\,\delta v_y)$ で、$\gamma^2\beta^2 = (\gamma-1)(\gamma+1)$ から $\gamma^2\beta/(\gamma+1) = (\gamma-1)/\beta$ なので、$\vec w = \dfrac{\gamma-1}{\beta}(0,\ \delta v_z,\ -\delta v_y)$。sympy の出力は、この 2 つの形の差が 3 成分とも 0 である。速度の変化 $\delta\vec v = \vec w_{\rm acc}\,dt$ ($\vec w_{\rm acc}$ は加速度) を入れると、静止系の軸の回転角速度は

$$
\vec\omega_T = \frac{\gamma^2}{\gamma+1}\,\vec w_{\rm acc}\times\vec v\qquad(\vec w_{\rm acc}:\text{加速度})
$$

(回転の向きの約束は、静止系で測ったスピンの歳差 (15.2) に足したとき、標準の BMT の式 (15.4) になるように合わせてある。確認 (a)。)

**磁場の場合**。磁場だけの力 $\vec w_{\rm acc} = \dfrac{q}{\gamma m}\vec v\times\vec B$ を入れ、$(\vec v\times\vec B)\times\vec v = \vec B(\vec v\cdot\vec v) - \vec v(\vec v\cdot\vec B) = \beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)$ ($\vec A\times\vec B\times\vec C$ の公式 $(\vec A\times\vec B)\times\vec C = \vec B(\vec A\cdot\vec C) - \vec A(\vec B\cdot\vec C)$ を $\vec A = \vec v$、$\vec C = \vec v$ に使った) を使うと、$\vec\omega_T = \dfrac{\gamma^2}{\gamma+1}\cdot\dfrac{q}{\gamma m}[\beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)]$ で、

$$
\vec\omega_T = \frac{q\gamma}{m(\gamma+1)}\left[\beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)\right]
\tag{15.3}
$$

### 3.3 和: BMT

実験室系の歳差は (15.2) + (15.3)。2 つを並べて書くと

$$
\vec\Omega = -\frac qm\Bigl[(1+a)\vec B - (1+a)\frac{\gamma}{\gamma+1}\vec\beta(\vec\beta\cdot\vec B)\Bigr] + \frac qm\frac{\gamma}{\gamma+1}\Bigl[\beta^2\vec B - \vec\beta(\vec\beta\cdot\vec B)\Bigr]
$$

で、$\vec B$ の項と $\vec\beta(\vec\beta\cdot\vec B)$ の項にまとめると:

$$
\vec\Omega = -\frac{q}{m}\left[\left(a + \frac1\gamma\right)\vec B - \frac{a\gamma}{\gamma+1}(\vec\beta\cdot\vec B)\vec\beta\right]
\tag{15.4}
$$

(整理: $q/m$ に対する $\vec B$ の係数は $-(1+a) + \gamma\beta^2/(\gamma+1)$。$\gamma\beta^2/(\gamma+1) = (\gamma^2-1)/(\gamma(\gamma+1)) = (\gamma-1)/\gamma = 1 - 1/\gamma$ なので、$-(1+a) + 1 - 1/\gamma = -(a + 1/\gamma)$。$q/m$ に対する $\vec\beta(\vec\beta\cdot\vec B)$ の係数は $(1+a)\gamma/(\gamma+1) - \gamma/(\gamma+1) = a\gamma/(\gamma+1)$。これが (15.4) の第 2 項の $-a\gamma/(\gamma+1)$ に当たる。) これが Thomas-BMT 方程式 (電場なし) である。**確認** (a): sympy で $\vec\Omega_{\rm rest} + \vec\omega_T - \vec\Omega_{\rm BMT}$ の 3 成分が 0 になる。

### 3.4 縦・横への分解

磁場を運動量に平行 ($\parallel$) と垂直 ($\perp$) に分ける。$\vec\beta\cdot\vec B = \beta B_\parallel$ なので、平行成分の係数は $a + \dfrac1\gamma - \dfrac{a\gamma\beta^2}{\gamma+1} = a + \dfrac1\gamma - \dfrac{a(\gamma-1)}{\gamma} = \dfrac{1+a}{\gamma}$、垂直成分の係数は $a + \dfrac1\gamma = \dfrac{1+a\gamma}{\gamma}$:

$$
\vec\Omega = -\frac{q}{\gamma m}\left[(1+a\gamma)\vec B_\perp + (1+a)\vec B_\parallel\right]
\tag{15.5}
$$

(sympy の出力が同じ形。) 運動量の運動方程式 $d\vec p/dt = q\vec v\times\vec B$ から、運動量の向き $\hat p$ は $d\hat p/dt = -\dfrac{q}{\gamma m}\vec B_\perp\times\hat p$ で回る。(15.5) の $\perp$ の部分は、運動量の回転と同じ向きに $(1+a\gamma)$ 倍の速さで回ることを表す。運動量に対する相対では $a\gamma$ 倍。

### 3.5 経路長で書く

$dt = ds/v$、$p = \gamma mv$ なので、$\dfrac{q}{\gamma m}dt = \dfrac{q}{\gamma mv}ds = \dfrac qp\,ds$。これを (15.5) の $\vec\Omega\,dt$ に入れると $\vec\Omega\,dt = -\dfrac qp\left[(1+a\gamma)\vec B_\perp + (1+a)\vec B_\parallel\right]ds$。運動量の向きの回転も同じ変数で書ける: $d\hat p/ds = -\dfrac qp\vec B_\perp\times\hat p$ ($d\hat p/dt$ の式を $v$ で割った)。そこで $\vec\theta' \equiv -\dfrac qp\vec B_\perp$ (運動量の向きの回転ベクトル) とおくと、

$$
\frac{d\vec S}{ds} = \vec\omega_s\times\vec S,\qquad \vec\omega_s = (1+a\gamma)\,\vec\theta' + (1+a)\,\psi'\,\hat p,\qquad \vec\theta' \equiv -\frac{q}{p}\vec B_\perp,\quad \psi' \equiv -\frac{q}{p}B_\parallel
\tag{15.6}
$$

$\vec\theta'$ は運動量の向きの回転ベクトル (単位長さあたり、軸 × 角)、$\psi'$ は符号つきの縦磁場の積分 $\int B_\parallel ds/(B\rho)$ の密度 ($B\rho = p/q$)。SAD の $b_t$ が $\vec\theta'$ の積分、$b_l$ が $\psi'$ の積分に当たる (§4)。

**例 (平面リング)**: 垂直の一様磁場だけなら $\vec\theta'$ は一定の軸 $\hat y$ のまわり ($\vec B_\parallel = 0$、$\psi' = 0$)。運動量の向き $\hat p$ は角速度 $\vec\theta'$ で回る座標系 (運動量に付いた座標系) を作る。この座標系から見たスピンの角速度は、$\vec\omega_s$ から座標系自身の角速度 $\vec\theta'$ を引いたもの: $\vec\omega_s - \vec\theta' = (1+a\gamma - 1)\vec\theta' = a\gamma\,\vec\theta'$。経路にそって積分すると、相対歳差は $a\gamma\,\theta$ ($\theta$ は偏向角)。完全なリング ($\theta = 2\pi$) で $\nu_s = a\gamma$。$\gamma = 5870.85$ (3 GeV) で $a\gamma = 6.808$ (`NominalSpinTune` $= h_0\cdot g = 6.80814822$)。
**縦磁場**: $\hat p$ に平行な磁場だけなら運動量は回らず、スピンは $(1+a)\psi$ だけ回る (確認 (b2): $\psi = 0.5$ で回転角 $0.50057983 = (1+a)\psi$)。

## 4. 1 要素 (スライス) のスピン回転

### 4.1 K64 の `sprot`

放射の 1 スライス (`tradk1`) の前後の運動量 $p_{\rm in}$、$p_{\rm out}$ ($p_x, p_y$ と $p_z = 1 + \mathtt{pxy2dpz}$) から、

$$
(\mathtt{ppx}, \mathtt{ppy}, \mathtt{ppz}) = \vec p_{\rm out}\times\vec p_{\rm in}
\tag{15.7}
$$

(`tradk1` では `ppx = py*dpz0 - dpz*py0 + dpy` など。$p_z = 1 + dp_z$ を使った外積の成分で、大きさは $\sin\theta$、向きは運動量の回転軸の**逆**。成分の導出: $\vec p_{\rm in} = (p_{x0}, p_{y0}, 1 + dp_{z0})$、$\vec p_{\rm out} = (p_x, p_y, 1 + dp_z)$ として、$x$ 成分は $p_{{\rm out},y}p_{{\rm in},z} - p_{{\rm out},z}p_{{\rm in},y} = p_y(1 + dp_{z0}) - (1 + dp_z)p_{y0} = p_y\,dp_{z0} - dp_z\,p_{y0} + (p_y - p_{y0})$ (最後の括弧が `dpy`)。$y$ 成分は $p_{{\rm out},z}p_{{\rm in},x} - p_{{\rm out},x}p_{{\rm in},z} = (1 + dp_z)p_{x0} - p_x(1 + dp_{z0}) = dp_z\,p_{x0} - p_x\,dp_{z0} - (p_x - p_{x0})$ (`ppy = dpz*px0 - px*dpz0 - dpx`)。$z$ 成分は $p_xp_{y0} - p_yp_{x0}$ (`ppz`)。$p_z$ を 1 からの差 $dp_z$ で持つのは、小さい角で $p_z\simeq1$ の桁落ちを避けるためと読める。) `sprot` は、中点の運動量 $\hat p_m = (\mathtt{pxm}, \mathtt{pym}, \mathtt{pzm})$ ($\mathtt{pzm} = 1 + \mathtt{pxy2dpz}(\mathtt{pxm}, \mathtt{pym})$) で割って

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

とする。(15.6) との対応: $\vec b_t = -\vec\theta$ (運動量の回転ベクトルの符号を反転したもの。ppx.. が $p_{\rm out}\times p_{\rm in}$ だから)、$\vec b_l = \hat p\,\mathtt{bsi}$。K64 は**スピンを $\exp(-[g]_\times)$ で回す**約束で、$-\vec g$ が右手系の回転ベクトルになり、(15.6) の $\vec\omega_s\,ds$ と一致する。**$\mathtt{pzm}$ で割る理由と弦の大きさ**。平面の偏向 (水平面内で、入口 $\hat p_{\rm in} = (0, 0, 1)$、出口 $\hat p_{\rm out} = (\sin\theta, 0, \cos\theta)$ 。$p_x$ と $p_y$ は $\hat p$ の成分) で成分を書き下すと、$\vec p_{\rm out}\times\vec p_{\rm in}$ の $y$ 成分は $p_{{\rm out},z}\,p_{{\rm in},x} - p_{{\rm out},x}\,p_{{\rm in},z} = \cos\theta\cdot0 - \sin\theta\cdot1 = -\sin\theta$、他の成分は 0 なので $(0, -\sin\theta, 0)$。大きさは $\sin\theta$ で、向きは運動量の回転軸 ($+\hat y$) の逆 (上の「逆」の内訳)。一方、中点の運動量は $p_{xm} = \sin\theta/2$、$\mathtt{pzm} = \sqrt{1 - p_{xm}^2} = \sqrt{1 - \sin^2\theta/4}$ で、これは $\cos(\theta/2)$ ($= \sqrt{1 - \sin^2(\theta/2)}$、入口と出口の二等分線の $z$ 成分) と $\sin^4(\theta/2)$ の差しかない。したがって $|\vec b_t| = \sin\theta/\mathtt{pzm} \simeq \sin\theta/\cos(\theta/2) = 2\sin(\theta/2)$ (弦の長さ、$\theta$ の代わり)。テイラー展開すると $2\sin(\theta/2) = \theta - \theta^3/24 + O(\theta^5)$ で、相対 $-\theta^2/24$ ($\theta^3$ の項は $\theta/2 - (\theta/2)^3/6 + \cdots$ の 2 倍から $2\cdot(\theta/2)^3/6 = \theta^3/24$)。$\mathtt{pzm}$ を使った実際の値は弦より相対 $-\theta^4/32$ だけさらに小さい (確認 [spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (j): $\theta = 0.2618$ で $-1.48\times10^{-4}$、$\theta^4/32 = 1.47\times10^{-4}$)。

回転は Rodrigues の式。導出: $[g]_\times$ ($[g]_\times\vec v = \vec g\times\vec v$) は $[g]_\times^2 = \vec g\vec g^\top - |g|^2I$、$[g]_\times^3 = -|g|^2[g]_\times$ を満たす (sympy で $K^3 + |g|^2K = 0$ と $K^2 - (gg^\top - |g|^2I) = 0$ を確認、[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (c))。べき級数 $\exp(-[g]_\times) = \sum_n(-[g]_\times)^n/n!$ を奇数べきと偶数べきに分けて、$[g]_\times^{2k+1} = (-|g|^2)^k[g]_\times$、$[g]_\times^{2k+2} = (-|g|^2)^k[g]_\times^2$ を入れると

$$
\exp(-[g]_\times) = I - \Bigl(\sum_{k}\frac{(-1)^k|g|^{2k}}{(2k+1)!}\Bigr)[g]_\times + \Bigl(\sum_k\frac{(-1)^k|g|^{2k}}{(2k+2)!}\Bigr)[g]_\times^2 = I - \frac{\sin|g|}{|g|}[g]_\times + \frac{1-\cos|g|}{|g|^2}[g]_\times^2
$$

(括弧の級数は $\sin|g|/|g|$、$(1-\cos|g|)/|g|^2$ のテイラー展開。) $\vec S$ にかけて $[g]_\times\vec S = \vec g\times\vec S = -\vec S\times\vec g$、$[g]_\times^2\vec S = \vec g(\vec g\cdot\vec S) - |g|^2\vec S$ を使うと、$\vec S - \dfrac{\sin|g|}{|g|}\vec g\times\vec S + \dfrac{1-\cos|g|}{|g|^2}[\vec g(\vec g\cdot\vec S) - |g|^2\vec S] = \cos|g|\,\vec S + \dfrac{\sin|g|}{|g|}\vec S\times\vec g + \dfrac{1-\cos|g|}{|g|^2}(\vec S\cdot\vec g)\vec g$ ($\vec S$ の項: $1 - (1-\cos|g|) = \cos|g|$)。この閉じた式を sympy と `scipy.linalg.expm` で確認した (差 $2\times10^{-15}$)。`xsincos(g, sinu, xsinu, cosu, dcosu)` が $\sin|g|$、$\cos|g|$、$\cos|g| - 1$ ($= -2\sin^2(|g|/2)$、桁落ちを避ける) を返し、

$$
\vec S' = \cos|g|\,\vec S + \frac{\sin|g|}{|g|}\,\vec S\times\vec g + \frac{1 - \cos|g|}{|g|^2}(\vec S\cdot\vec g)\,\vec g
\tag{15.10}
$$

(`sw = -(S·g) dcosu/g²`、`sinu = sin|g|/|g|`、`sux = S×g`)。(15.10) は $\exp(-[g]_\times)\vec S$ の閉じた式 ($-\vec g\times\vec S = \vec S\times\vec g$ だから)。

### 4.2 刻みの誤差

(15.9)〜(15.10) の列は、Thomas-BMT の厳密な積分に収束する。**確認** (b): $\hat B$ と運動量が $0.3$ rad の角をなす一様磁場 (らせん軌道、3 GeV、$\kappa L = 2$ rad、縦成分あり) で、(15.6) の DOP853 による積分を基準にすると、$\vec\omega_s$ の刻みごとの回転 $\mathrm{Rot}[(1+a\gamma)\vec\theta_{\rm step} + (1+a)\psi\,\hat p_m]$ の積は、刻み数 $N = 10, 20, 40, 80, 160$ で誤差 $1.08\times10^{-2}$、$2.69\times10^{-3}$、$6.73\times10^{-4}$、$1.68\times10^{-4}$、$4.21\times10^{-5}$、誤差の比はいつも $4.00$ で **2 次**。理由: 1 刻み ($\Delta s$) の厳密な転送は時間順序の指数関数で、刻みの中の角速度を $\vec\omega_s(s) = \vec\omega_m + \vec\omega_m'(s - s_m) + \cdots$ と中点のまわりに展開すると、1 次の項は 1 刻みの中で奇関数なので積分が 0 になり、2 次の項 (大きさ $\Delta s^3$ の $\vec\omega_m'\times\vec\omega_m$ 型の交換子) から誤差が出る (1 刻みの誤差 $O(\Delta s^3)$)。刻みの数は $1/\Delta s$ なので、全体の誤差は $O(\Delta s^2)$、$N$ を 2 倍にすると 1/4。K64 は縦成分の向きに中点の $\hat p_m$ を使い、この中点則の形になっている。

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

**なぜこれで $D$ になるか**。スライス $i$ ($i = 1, \dots, N$) の入口で座標が $x_{{\rm in},i} = M_{0\to i}x_0$ の粒子は、閉軌道上の粒子の $\vec g_i$ から $\delta\vec g_i = G_ix_0$、$G_i \equiv (\partial\vec g_i/\partial x_{\rm in})M_{0\to i}$ (3×6) だけずれた回転ベクトル $\vec g_i + \delta\vec g_i$ で回る。1 スライスごとの $|g_i|$ が小さいとして、指数写像の微分 (一般には $\exp(-[g+\delta g]_\times) = \exp(-[g]_\times)\exp(-[J(g)\,\delta g]_\times)$、$J(g) = I + O(|g|)$ の補正行列) の補正 $O(|g_i|)$ を落とすと

$$
\exp(-[g_i + \delta g_i]_\times)\ \simeq\ \exp(-[g_i]_\times)\,\exp(-[\delta g_i]_\times)
$$

(閉軌道の回転 $R_i = \exp(-[g_i]_\times)$ の前に、微小な追加の回転 $\exp(-[\delta g_i]_\times)$ がかかる形)。ここまでのスライスの回転が $\exp(-[D^{(i-1)}x_0]_\times)R^{(i-1)}$ ($R^{(i-1)} = R_{i-1}\cdots R_1$、$D^{(i-1)}$ は積算済みの $D$) の形だと仮定して、スライス $i$ を作用させると、微小な回転どうしは (2 次の微小量を除き) 回転ベクトルが足し算になるので

$$
R_i\exp(-[\delta g_i]_\times)\exp(-[D^{(i-1)}x_0]_\times)R^{(i-1)} \simeq R_i\exp(-[\delta g_i + D^{(i-1)}x_0]_\times)R^{(i-1)} = \exp(-[R_i(G_i + D^{(i-1)})x_0]_\times)\,R_iR^{(i-1)}
$$

(最後の等号は、一般に回転行列 $B$ について $B\exp(-[v]_\times)B^{-1} = \exp(-[Bv]_\times)$ であることによる。) したがって $D^{(i)} = R_i\,(D^{(i-1)} + G_i)$ ($D^{(0)} = 0$)。これを $N$ スライスまで繰り返すと

$$
R(x_0) \simeq \exp\!\left(-[D\,x_0]_\times\right)R_0,\qquad D = \sum_{i=1}^{N}R_N\cdots R_{i+1}R_i\,\frac{\partial\vec g_i}{\partial x_{\rm in}}M_{0\to i}
\tag{15.12}
$$

(和の各項は、スライス $i$ の追加の回転を、スライス $i$ とそれ以降の閉軌道の回転 $R_N\cdots R_i$ で出口の座標系へ運んだもの。) 手順 3 が「$R_i$ を $D$ (すでに積算した分と今の分の和) の全部の列にかける」ことに、手順 2 が各スライスの寄与 $G_i$ を足すことに当たる。この積算規則は、玩具の $N$ 個の等しいスライスで有限差分の $\partial R/\partial x_0$ と比べると、誤差がスライス数に比例して減る (4、8、16、32、64 スライスで $1.4\times10^{-1}$、$7.1\times10^{-2}$、$3.5\times10^{-2}$、$1.8\times10^{-2}$、$8.8\times10^{-3}$。[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (g))。これは落とした $O(|g_i|)$ の補正に当たり、スライスを細かくすると 1 次で消える (K64 の $D$ の刻みの誤差。$|g_i|$ が大きい刻みでは精度が落ちる)。したがって、**$D$ は 1 周の回転に対する座標の影響を表す 3×6 行列で、$R$ と同じ座標系 (1 周の終わり) で測る**。$\partial\vec g/\partial x$ の中の $\partial b_t/\partial x$ は偏向の位置の依存 (勾配・フリンジ)、$\partial/\partial\delta$ の項は $a\gamma$ の $\delta$ 依存と運動量 ($B\rho$) の変化から来る。

### 5.2 `spnorm`: 偏極ベクトルとスピンチューン

$R_0 = \mathtt{srot(:,1:3)}$ は、積算の丸め誤差のため直交からわずかにずれている。`spnorm` (`photons.f:767-839`) は次の順で処理する。

1. 列 2 を規格化、列 1 から列 2 の成分を引いて規格化、列 3 は外積 (直交化)。
2. 回転の軸の初期値: 反対称部分 $\mathtt{sps(:,1)} = (R_{23} - R_{32},\ R_{31} - R_{13},\ R_{12} - R_{21})$ (ソースの `sps(1,1)=srot(2,3)-srot(3,2)` など) を規格化。理由: K64 の約束の $R = \exp(-\mu[\hat n]_\times) = I - \sin\mu\,[\hat n]_\times + (1-\cos\mu)[\hat n]_\times^2$ (§4.1) で、$[\hat n]_\times^2$ は対称、$[\hat n]_\times$ は反対称で、$([\hat n]_\times)_{23} = -n_x$、$([\hat n]_\times)_{31} = -n_y$、$([\hat n]_\times)_{12} = -n_z$ なので、$R_{23} - R_{32} = -\sin\mu\,(-n_x) \cdot 2 = 2\sin\mu\,n_x$ (他の成分も同様) となり、$\mathtt{sps(:,1)} = 2\sin\mu\,\hat n$。大きさが $10^{-4}$ 未満 ($\sin\mu \approx 0$、$\mu \approx 0, \pi$) なら、実固有値の固有ベクトル (`teigen`)。
3. $(R - 1)\,\delta\vec n = \vec n - R\vec n$ の最小ノルム解 (`tsolvg`) で軸を補正し、規格化。理由: $\vec n + \delta\vec n$ が $R(\vec n + \delta\vec n) = \vec n + \delta\vec n$ を満たすには、$(R - 1)(\vec n + \delta\vec n) = 0$、つまり $(R-1)\delta\vec n = -(R-1)\vec n = \vec n - R\vec n$。成分の符号が負側に偏っていれば反転 (軸の向きは固有値問題では決まらないので、約束で正側にそろえる)。これが $n_0$ (`PolarizationVector`)。
4. $n_0$ に垂直な基底 $\vec e_2, \vec e_3 = n_0\times\vec e_2$ を作る (ソース: `sps(:,3)=outer(sps(:,1),sps(:,2))`)。$R\vec e_2$ の成分から $c_m = \vec e_2\cdot R\vec e_2$、$s_m = \vec e_3\cdot R\vec e_2$ を出す。$R\vec e_2 = \cos\mu\,\vec e_2 - \sin\mu\,\hat n\times\vec e_2 = \cos\mu\,\vec e_2 - \sin\mu\,\vec e_3$ なので $c_m = \cos\mu$、$s_m = -\sin\mu$。したがって $\mu_s = \mathrm{atan2}(-s_m, c_m)$ ($= \mu$)、$\nu_s = \mu_s/2\pi$ (`SpinTune`)。
5. $\mathtt{sdamp} = n_0\cdot\mathtt{gintd}$ (Sokolov-Ternov の 1 周の率、§7)。

**確認** (e): 平面リング (3 GeV、12 セル、偏向 15 度 $\times$ 12) で $n_0 = (0, 1, 0)$、$\nu_s = 0.40407411$ ($a\gamma\theta_{\rm tot}/2\pi = 3.404074$ の端数 $0.404074$ と一致。放射損失と RF による $\gamma$ の補正は、この桁では見えない)。

## 6. SLIM: 軌道が変わると軸がどれだけ傾くか

### 6.1 $N$ の式

粒子の座標 $x_0$ (リング入口) に対するスピンの周期解の軸を $n(x_0)$ とする。軸は、1 周の回転と 1 周の写像 $x_1 = Mx_0$ に対して $R(x_0)\,n(x_0) = n(Mx_0)$ を満たす。1 次まで $n(x_0) = n_0 + Nx_0$ ($N$: 3×6、$n_0\cdot N = 0$) とおく。(15.12) より $R(x_0)\,n(x_0) = \exp(-[Dx_0]_\times)(R_0n_0 + R_0Nx_0) \approx n_0 + R_0Nx_0 - (Dx_0)\times n_0 = n_0 + R_0Nx_0 + n_0\times(Dx_0)$ ($R_0n_0 = n_0$)。右辺は $n_0 + NMx_0$。$x_0$ の係数を比べて

$$
N M = R_0\,N + [n_0]_\times D\qquad\Longleftrightarrow\qquad N = \left(R_0N + [n_0]_\times D\right)M^{-1}
\tag{15.13}
$$

($[n_0]_\times$ は外積の行列、$[n_0]_\times v = n_0\times v$。)

**(15.13) の解の存在と一意性**。(15.13) の第 2 形は $N - R_0NM^{-1} = [n_0]_\times DM^{-1}$ で、$\mathrm{vec}(AXB) = (B^\top\otimes A)\,\mathrm{vec}(X)$ を使うと

$$
\left(1 - M^{-\top}\otimes R_0\right)\mathrm{vec}(N) = \mathrm{vec}\left([n_0]_\times DM^{-1}\right)
$$

の線形方程式 (18×18) になる。係数行列の固有値は $1 - \lambda_k^{-1}\,e^{\mp i\mu_s}$ と $1 - \lambda_k^{-1}$ ($\lambda_k$: $M$ の固有値、$R_0$ の固有値は $e^{\pm i\mu_s}$ と 1)。放射による減衰 ($|\lambda_k| < 1$) があれば $1 - \lambda_k^{-1}\ne0$、スピン共鳴 ($\lambda_k^{-1}e^{\mp i\mu_s} = 1$) でなければ前者も 0 にならないので、解は一意。さらに、(15.13) の左から $n_0^\top$ をかけると ($n_0^\top R_0 = n_0^\top$、$n_0^\top[n_0]_\times = 0$) $(n_0^\top N)M = n_0^\top N$ で、$M$ が固有値 1 を持たないので $n_0\cdot N = 0$ (最初の仮定と矛盾しない)。分母 $1 - e^{i(\mu_s\pm\mu_k)}$ が 0 に近づくのが、線形のスピン軌道共鳴 (チューン $\nu_s\pm\nu_x$、$\nu_s\pm\nu_z$ など)。

**要素ごとの形**。要素 $i$ を $x_{\rm out} = T_ix_{\rm in}$、$\vec S_{\rm out} = \exp(-[D_ix_{\rm in}]_\times)R_i\vec S_{\rm in}$ ($D_i$ は出口の座標系で測った結合、閉軌道上の軸の対応 $n_{i+1} = R_in_i$) と書く。入口の軸 $n_i(x_{\rm in}) = n_i + N_ix_{\rm in}$ に対して出口の軸は、$n_{i+1}(T_ix_{\rm in}) = \exp(-[D_ix_{\rm in}]_\times)R_i(n_i + N_ix_{\rm in}) \simeq n_{i+1} + R_iN_ix_{\rm in} + n_{i+1}\times(D_ix_{\rm in})$。左辺は $n_{i+1} + N_{i+1}T_ix_{\rm in}$ なので

$$
N_{i+1}T_i = R_iN_i + [n_{i+1}]_\times D_i
\tag{15.13'}
$$

を得る。これを 1 周でつなぐと (15.13) になる: 乱数の玩具の 6 要素で、$M = T_6\cdots T_1$、$R_0 = R_6\cdots R_1$、$D = \sum_iR_6\cdots R_{i+1}D_iT_{i-1}\cdots T_1$ (積算規則は (15.12) と同じ) を作ると、(15.13) の解から (15.13') で要素ごとに $N_i$ を送って 1 周後に元の $N$ に戻る (差 $9\times10^{-15}$、[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (h))。

**確認** (d): $R_0$ を $n_0 = \hat y$ のまわりの回転 ($\nu_s = 0.403$)、$M$ を減衰 0.9 の 3 つの回転、$D$ を乱数にした玩具の 1 周写像で、(15.13) を解き ($n_0\cdot N = 2\times10^{-16}$)、座標 $x_0$ の方向を固定して大きさを $|x_0| = 10^{-2}, 10^{-3}, 10^{-4}$ と変えると、残差 $|R(x_0)\,n(x_0) - n(Mx_0)|$ は $1.27\times10^{-6}$、$1.27\times10^{-8}$、$1.27\times10^{-10}$ (比はちょうど 100 で、$|x_0|^2$ に比例)。これは (15.13) の代数と 1 次の展開が自己無撞着であることの検査で、K64 の $D$・$N$ の検証ではない (玩具の $D$ は乱数)。確率過程としての意味 (脱偏極率) は §6.2 の確認で見る。

### 6.2 脱偏極率

放射が起こらなければ、スピン $\vec S$ は軸 $n(x)$ に追従し (断熱不変量 $\vec S\cdot n(x)$ が保存)、偏極は保たれる。1 個の光子の放射は、粒子の座標を $x\to x+\xi$ ($\xi$: 跳び、平均 0、共分散 $B_i$、`tradke` の `dbeam`) と変えるが、スピンは変えない (光子が運ぶスピンの寄与は Sokolov-Ternov 項として別に扱う)。軸は $n(x)\to n(x+\xi)$ に変わる。軌道が閉軌道に近い ($|Nx|\ll1$) とすると、$n(x) = n_0 + a$、$n(x+\xi) = n_0 + a + b$ ($a = Nx$、$b = N\xi$。$n_0\cdot a = n_0\cdot b = 0$) で、規格化した軸の内積は

$$
\hat n(x)\cdot\hat n(x+\xi) = \frac{1 + a^2 + a\cdot b}{\sqrt{1+a^2}\,\sqrt{1 + a^2 + 2a\cdot b + b^2}} = 1 - \tfrac12|b|^2 + O(|a||b|^2,\ |b|^3)
$$

(分母の展開は $\sqrt{1+u}\simeq1 + u/2 - u^2/8$ を $u = a^2$ と $u = a^2 + 2a\cdot b + b^2$ で使い、$|a|$ の 1 次までと $|b|$ の 2 次までを残した。$|a|\to0$ では分子 $= 1 + a\cdot b + a^2$ が 1 に、分母が $\sqrt{1+b^2} \simeq 1 + b^2/2$ になって直ちに $1 - \tfrac12|b|^2$。) スピンは跳びの直前の軸の向き ($\vec S = \hat n(x)$) にあるので、

$$
\vec S\cdot\hat n(x+\xi) = 1 - \tfrac12|N\xi|^2 + O(\xi^3)
\tag{15.14}
$$

(跳びのあと、スピンは新しい軸から角 $|N\xi|$ だけ離れ、軸のまわりを歳差するので、軸に垂直な成分は平均で 0 に位相がばらけて、軸の方向の成分だけが残る。つまり、偏極 $\langle\vec S\cdot\hat n\rangle$ は $1 - \tfrac12|N\xi|^2$ 倍に減ったまま続く。連続する跳びが独立 (相関しない) で、1 回の減りが小さい ($|N\xi|\ll1$) という近似。) $\langle\xi\xi^\top\rangle = B_i$ で平均すると $\langle|N\xi|^2\rangle = \mathrm{tr}(NB_iN^\top)$ ($|N\xi|^2 = \xi^\top N^\top N\xi$ の期待値)。したがって、1 周あたりの偏極の減りの率 (脱偏極率) は

$$
\frac{1}{\tau_{\rm dep}\,f_{\rm rev}} = \frac12\sum_i\mathrm{tr}\left(N_iB_iN_i^\top\right)
\tag{15.15}
$$

(和は 1 周の放射のスライス $i$ について。$N_i$: そのスライスの出口での $N$。) $B_i$ のエネルギーの跳びの成分は $B_{\delta\delta} = N_\gamma\,u_c^2\,c_{\rm uu}$ ([13_radiation.md](13_radiation.md) の `dee`)。これが §7.3 で DK の係数 $11/18$ に対応する。

**確認 (モンテカルロ)**: 玩具の軌道 $x_{n+1} = Mx_n + \xi_n$ ($M$: 弱い減衰 $\lambda$ を持つ 3 つの回転、$\xi\sim N(0, B)$)、スピン $\vec S_{n+1} = \exp(-[Dx_n]_\times)R_0\vec S_n$ (厳密な回転) で、2000 個のスピンを平衡の分布から出発させて追跡し、偏極 $P(n) = \langle\vec S_n\cdot\hat n(x_n)\rangle$ ($\hat n(x_n)$ は (15.13) の $N$ から作った軸) の減り方から率を読んだ ([spin_03_slim_mc.out.txt](outputs/spin_03_slim_mc.out.txt))。(15.15) の $r = \tfrac12\mathrm{tr}(NBN^\top)$ との比 (測定/予想) は、4 つの条件 ($\lambda = 10^{-2}$〜$3\times10^{-3}$、$|N|_{\max} = 0.6$〜2.7、$x$ の rms $\le0.05$、率 $10^{-5}$〜$6\times10^{-5}$ /周) で $0.977$、$1.024$、$1.045$、$1.037$。統計誤差は 1〜2% 程度なので、(15.15) は 5% 以内で成り立つ (残りは $|N x|$ の 2 次以降と統計)。K64 の $D$・$N$ を使った確認ではなく、(15.15) の式の確率過程としての検査。

### 6.3 平衡偏極度

Sokolov-Ternov の偏極 (立ち上がりの率 $\tau_{ST}^{-1}$、到達値 $P_{ST}$) と脱偏極 (率 $\tau_{\rm dep}^{-1}$) を合わせた偏極の運動方程式 $dP/dt = (P_{ST} - P)/\tau_{ST} - P/\tau_{\rm dep}$ で、$dP/dt = 0$ とおくと $P(\tau_{ST}^{-1} + \tau_{\rm dep}^{-1}) = P_{ST}\tau_{ST}^{-1}$ なので

$$
P_{\rm eq} = P_{ST}\,\frac{\tau_{\rm dep}}{\tau_{ST} + \tau_{\rm dep}} = P_{ST}\,\frac{\mathtt{sdamp}}{\mathtt{sdamp} + \mathtt{rate}}
\tag{15.16}
$$

($\tau^{-1} = (\text{1 周あたりの率})\times f_{\rm rev}$ なので、$f_{\rm rev}$ は分子分母で約分される。`rate` は (15.15) の 1 周あたりの率、`sdamp` は 1 周あたりの Sokolov-Ternov の率。) v2 の `Polarization` はこの形 (`slim_depolarization`)。K64 の `srequpol2` (`photons.f:861-907`) も同じ形の式 $P = \mathtt{pst}\,\mathtt{sdamp}/(\mathtt{rmd} + \mathtt{sdamp})$ を持つが、**どこからも呼ばれない** (呼ばれるのは `srequpol` だけ)。`rmd` の分母に $|e^{-a_k + i(\mu_s + \mu_k)} - 1|^2$ が現れ、(15.13) の共鳴の分母に当たる (`rmd` の他の係数は読み取っていない)。

### 6.4 正規座標での $N$ と共鳴の分母

$N$ を、$M$ の正規座標とスピンの基底 $(n_0, \vec e_2, \vec e_3 = n_0\times\vec e_2)$ で書く。$N$ の $n_0$ に垂直な 2 行を複素にまとめて $z = (\vec e_2 + i\vec e_3)^\top N$ ($1\times6$、複素)、同じく $d = (\vec e_2 + i\vec e_3)^\top D$ とおく。$R_0\vec e_2 = \cos\mu_s\vec e_2 - \sin\mu_s\vec e_3$、$R_0\vec e_3 = \sin\mu_s\vec e_2 + \cos\mu_s\vec e_3$ ($\mu_s$ は §5.2 の約束) から $(\vec e_2 + i\vec e_3)^\top R_0 = e^{-i\mu_s}(\vec e_2 + i\vec e_3)^\top$ (展開: $R_0^\top\vec e_2 = \cos\mu_s\vec e_2 + \sin\mu_s\vec e_3$、$R_0^\top\vec e_3 = -\sin\mu_s\vec e_2 + \cos\mu_s\vec e_3$ を足すと $(\cos\mu_s - i\sin\mu_s)\vec e_2 + (\sin\mu_s + i\cos\mu_s)\vec e_3 = e^{-i\mu_s}(\vec e_2 + i\vec e_3)$)。また $n_0\times(\vec e_2 + i\vec e_3) = \vec e_3 - i\vec e_2 = -i(\vec e_2 + i\vec e_3)$ から $(\vec e_2 + i\vec e_3)^\top[n_0]_\times = i(\vec e_2 + i\vec e_3)^\top$ ($a^\top[n_0]_\times v = a\cdot(n_0\times v) = v\cdot(a\times n_0) = -v\cdot(n_0\times a)$ を使った)。(15.13) の左から $(\vec e_2 + i\vec e_3)^\top$ をかけると

$$
z\,M = e^{-i\mu_s}\,z + i\,d
$$

$M$ を複素の固有値 $\lambda_{k\pm} = e^{-a_k\pm i\mu_k}$ と固有ベクトルで対角化し ($M = V\Lambda V^{-1}$)、$\tilde z = zV$、$\tilde d = dV$ とおくと、各成分で $\tilde z_{k\pm}(\lambda_{k\pm} - e^{-i\mu_s}) = i\tilde d_{k\pm}$ なので

$$
\tilde z_{k\pm} = \frac{i\,e^{i\mu_s}\,\tilde d_{k\pm}}{e^{-a_k\pm i\mu_k + i\mu_s} - 1}
\tag{15.13''}
$$

($\lambda_{k\pm} - e^{-i\mu_s}$ に $e^{i\mu_s}$ をかけて整理した。) 分母の絶対値が 0 に近づくのが、$\mu_s\pm\mu_k\simeq2\pi\times$整数 の線形のスピン軌道共鳴で、その幅は減衰率 $a_k$。(15.13'') は玩具で連立方程式の解と $3\times10^{-16}$ で一致する ([spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (k))。脱偏極率 (15.15) の $\tfrac12\mathrm{tr}(NBN^\top) = \tfrac12\langle|z\xi|^2\rangle$ は、$\xi$ を複素の正規振幅で展開すると、各モードの寄与 $\tfrac12|\tilde z_{k\pm}|^2\times(\text{そのモードの励起})$ の和になる。K64 の `srequpol2` (未使用) の分母 `cexp1(dcmplx(-damp, smu+amu))` は (15.13'') の分母と同じ形である。

### 6.5 K64 の `srequpol`・`spdepol` の構造と、導出できた範囲

**コードの構造 (ソースの読み)**。`spdepol` (`photons.f:620`) は、1 周のスピン回転の摂動を、モードごとの複素振幅 $\zeta_i,\bar\zeta_i$ ($i = x,y,z$) の単項式で展開した行列の表 (`rm(m)`、`scmat`。単項式の指数と、減衰・位相・スピン位相の個数を `ind(13)` に持つ) を使う。`mord = 6` で 6 次までの表を確保するが、実際に回すのは次だけ (ソースの読み):
1. 1 次の表 `rm(1)` を作る (各モードの $D$ の係数 `gx, gy, gz` から。`0.25*gzc*[...]` などの $3\times3$ 行列は、$n_0$、$\vec e_2$、$\vec e_3$ の基底での回転の生成子の複素成分)。
2. `spcalcres(..., k=2)` だけを回す (`do k=2,2`、他の次数の呼び出しは注釈の中)。`rm(2) = rm(1)·rm(1)/2` ($\exp(A)$ の 2 次の項 $A^2/2!$ の漸化式 $A^m/m! = (1/m)\,A\cdot A^{m-1}/(m-1)!$ の $m = 2$)、`rmi(2) = rm(2) + rm(1)·rmi(1)` (先の周回の項と今の周回の項の積の和)。
3. `spintrm` が `rmi` の各単項式の係数に、周回の和 (下の幾何級数) を掛ける。
4. 最後に `rmd(:,:,1)` だけが使われる (4 次・6 次の `rmd(:,:,2)`、`rmd(:,:,3)` は注釈)。`equpol` $= \mathtt{pst}\,\mathtt{sdamp}/(\mathtt{sdamp} - \mathtt{rm1}(1,1))$。

つまり K64 が実際に計算するのは**2 次の摂動の項だけ**である (6 次までの表は確保されるが、使われない)。

**2 次の相関和から SLIM が出ること (導出)**。周回ごとのスピン回転を $U_n = \exp(A_n)$ ($A_n$: 軌道の振幅に 1 次の小さな反対称行列) と書くと、$N$ 周後の平均は $\langle U_N\cdots U_1\rangle \simeq 1 + \sum_n\langle A_n\rangle + \sum_n\tfrac12\langle A_n^2\rangle + \sum_{n_2>n_1}\langle A_{n_2}A_{n_1}\rangle + \cdots$。$\langle A\rangle = 0$ なので、周回数に比例して積算される 2 次の項の 1 周あたりの量は $\tfrac12\langle A^2\rangle + \sum_{m\ge1}\langle A_{n+m}A_n\rangle$。軌道が $x_{n+1} = Mx_n + \xi_n$ の定常過程なら、相関は $\langle x_{n+m}x_n^\top\rangle = M^m\Sigma$ で、あるモードの寄与は $\varepsilon|d|^2\bigl[\tfrac12 + \sum_{m\ge1}q^m\bigr]$ の実部 ($q = e^{-a_k\pm i(\mu_s\pm\mu_k)}$ の型、$\varepsilon = \langle|x|^2\rangle$)。幾何級数 $\sum_{m\ge1}q^m = q/(1-q)$ から $\tfrac12 + q/(1-q) = \tfrac12(1+q)/(1-q)$、実部は $\tfrac12(1 - |q|^2)/|1-q|^2$ (sympy で確認、[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (l))。$|q|^2 = e^{-2a_k}$、平衡エミッタンスと励起の関係 $\varepsilon(1 - e^{-2a_k}) = \bar B_k$ (14 章 (14.5)) を使うと

$$
\text{脱偏極率(モードの寄与)} = \frac12\,\bar B_k\,\frac{|d|^2}{|1-q|^2}
$$

で、(15.15) と (15.13'') の組み合わせに一致する。つまり、2 次の摂動の結果は、平衡エミッタンス $\varepsilon$ ではなく、**1 周の励起 $\bar B_k$** で決まる (減衰率に依らない)。

**K64 の式との関係 (推定)**。K64 の `spintrm` は、各単項式に $\mathtt{cm0} = -\mathtt{mat}/(e^{c\varphi} - 1) = \mathtt{mat}\sum_{m\ge0}e^{mc\varphi}$ (幾何級数。$c\varphi$ は、減衰 $-\mathtt{ind}\times$(減衰率) と位相からなる複素数。ソースの `cphi=dcmplx(-ind(7)*dx-ind(8)*dy-ind(9)*dz, ind(10)*amx+...)`) を掛ける。単項式 $\zeta_i\bar\zeta_i$ (位相が打ち消し合う、`ind(10:13) = 0`) では、2 つの因子から減衰が $e^{-2a_k}$ になり、$c\varphi = -2a_k$、$\mathtt{cm0}\simeq\mathtt{mat}/(2a_k)$ となる。これに、`rm(1)` の係数に入っている (平衡の) エミッタンス $\varepsilon\simeq\bar B_k/(2a_k)$ が掛かるなら、上の導出の $\varepsilon\cdot\tfrac12(1-e^{-2a_k})/|1-q|^2$ に対して $1/(1-e^{-2a_k})\simeq1/(2a_k)$ 倍だけ過大になる (たとえば $a_k = 5\times10^{-4}$ で約 1000 倍。これは docs/30 の B125 の「`Polarization` の率が独立な値の約 1000 倍」と同じ桁)。docs/30 の B125 の位置の特定は、数値の依存性 (エミッタンスに比例、減衰率に反比例) による。**本章の再導出は、この二重の分母が原因である可能性を示す整合性の議論で、K64 の係数 (`rm(1)` の `0.25`、`0.5` の因子、`serot` の `se12`、`de12`) までは再導出していない。原因は推定**。

**ガウス広がりの項 (`a ≠ 0` の枝)**。$\mathtt{a} = \sqrt2\,|\mathtt{ind(13)}|\,\mathtt{ssprd}$ ($\mathtt{ssprd}$: スピンチューンの 1 周あたりの広がり、`drot(1,:)` と $\varepsilon$ から)。位相の広がりがガウス ($e^{-m^2a^2/4}$) のとき、周回の和を積分に置き換えると

$$
\int_0^\infty e^{-ct - (a/2)^2t^2}dt = \frac{\sqrt\pi}{a}\,e^{c_b^2}\,\mathrm{erfc}(c_b),\qquad c_b = \frac ca
\tag{15.23}
$$

($c = -c\varphi$ の実数の場合の簡単な形。微分方程式 $dI/dc = -(1 - cI)/(2a_g^2)$ ($a_g = a/2$)、$I(0) = \sqrt\pi/(2a_g)$ を満たすことを sympy で確認し、数値積分と一致 (10 桁)、(i)・(m))。K64 の $\mathtt{cb} = -(e^{-c\varphi}-1)/a$ は、$|c\varphi|\ll1$ で $\simeq c\varphi/a$ となり (15.23) の $c_b$ と同じ形になる (周回の和を離散のまま扱うための補正)。一方、前の係数は $\sqrt\pi/a$ でなく $2\pi\sqrt\pi/a$ ($\mathtt{cm} = -\mathtt{mat}\cdot2\pi\sqrt\pi\,e^{c_b^2-c\varphi}\mathrm{erfc}(c_b)/a$)。$a\to0$ の極限で $\mathtt{cm}/\mathtt{cm0}\to2\pi$ (数値: $a = 10^{-2}, 10^{-3}, 10^{-4}, 10^{-5}$ で $6.2866, 6.28322, 6.283186, 6.283185$) となり、$a = 0$ の枝 (`cm = cm0`) には $2\pi$ がない。したがって、`a ≠ 0` の枝の係数 $2\pi$ は、(15.23) に比べて $2\pi$ 倍大きい (読みと数値による。docs/30 の B125 も同じ指摘)。

**NaN の機構 (B55)**。`cerfc(z)` は $|z| > 18$ (`cezlim`) で漸近式 $\mathrm{erfc}(z)\simeq e^{-z^2}(1 - x + 3x^2 - \cdots)/(z\sqrt\pi)$ ($x = 1/(2z^2)$) に切り替わる。`spintrm` は $e^{c_b^2 - c\varphi}$ と `cerfc` を別々に評価して掛けるので、$c_b$ が実数で大きいとき、倍精度で $e^{c_b^2}$ は $c_b^2 > 709.78$ ($c_b > 26.64$) で桁あふれして ∞ になり、$e^{-c_b^2}$ は 0 に潰れ、積は ∞ または NaN になる (数値: $c_b = 26$ までは積は有限で $1/(c_b\sqrt\pi)$ に一致、$26.7$ と $27$ で ∞、$30$ で NaN。[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (i))。積自体は有限 ($\simeq1/(c_b\sqrt\pi)$) なので、$e^{c_b^2}\mathrm{erfc}(c_b)$ を直接 (Faddeeva 関数などで) 評価すれば有限である。docs/30 が書く「$|cb| > 18$」は漸近式に切り替わる境界で、積が桁あふれする境界は $|cb|\simeq26.6$ である (ここで確認した範囲)。B55 の NaN は gs で再現されている (確定)。原因を上の桁あふれとするのは推定で (docs/30: pysad で同じ式を写すと同じ NaN になり、桁あふれしない形に直すと有限になる)、実際のリングでの $c_b$ の値は本章では確かめていない。

**導出できていないこと**。(i) `rm(1)` の各係数 (`0.25`、`0.5`、$i$ の因子) の、$\exp(-[Dx]_\times)$ から複素振幅への展開による導出、(ii) `serot` が作る $\mathtt{se12}$、$\mathtt{de12}$ (光子の放射の跳びを、モードの複素振幅に写した励起) の導出、(iii) `ind` の添字の管理 (単項式の足し算) の正しさ、(iv) 4 次・6 次の項 (使われていない)。これらが足りないので、K64 の `Polarization` の式全体は再導出できておらず、原因の確定は gs での確認 (別途進めている) を待つ。v2 はこの式を使わず、(15.15)〜(15.16) と側帯の項で置き換えた ([docs/16_spin_pol.md](../docs/16_spin_pol.md) §5.3、B125)。

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

**QED の入力 (引用) から (15.17) を組み立てる**。QED の計算 (ここでは導かない) の結果は、磁場に沿った向き ($\uparrow$) から逆向き ($\downarrow$) へのスピン反転率と、その逆の率が

$$
w_{\uparrow\downarrow} = \frac{5\sqrt3}{16}\,w_0\left(1 + \frac{8\sqrt3}{15}\right),\qquad
w_{\downarrow\uparrow} = \frac{5\sqrt3}{16}\,w_0\left(1 - \frac{8\sqrt3}{15}\right),\qquad
w_0 = \frac{r_e\bar\lambda_C\,c\,\gamma^5}{|\rho|^3}
$$

であることを与える ($w_{\uparrow\downarrow}$ は、偏極の増える向き (到達する向き) へ反転する率、$w_{\downarrow\uparrow}$ はその逆向きへの率)。到達する向きの粒子数を $N_\uparrow$、逆向きを $N_\downarrow$ ($N = N_\uparrow + N_\downarrow$、$N_{\uparrow,\downarrow} = N(1\pm P)/2$) とすると、$dN_\uparrow/dt = w_{\uparrow\downarrow}N_\downarrow - w_{\downarrow\uparrow}N_\uparrow$ で、$P = (N_\uparrow - N_\downarrow)/N$ の変化は $dP/dt = 2\,dN_\uparrow/(N\,dt) = w_{\uparrow\downarrow}(1 - P) - w_{\downarrow\uparrow}(1 + P) = (w_{\uparrow\downarrow} - w_{\downarrow\uparrow}) - (w_{\uparrow\downarrow} + w_{\downarrow\uparrow})P$ なので、$\tau_{ST}^{-1} = w_{\uparrow\downarrow} + w_{\downarrow\uparrow}$、$P_{ST} = (w_{\uparrow\downarrow} - w_{\downarrow\uparrow})/(w_{\uparrow\downarrow} + w_{\downarrow\uparrow})$。上の率を入れると $w_{\uparrow\downarrow} + w_{\downarrow\uparrow} = \tfrac{5\sqrt3}{16}w_0\cdot2 = \tfrac{5\sqrt3}{8}w_0$ で $\tau_{ST}^{-1}$ の式、$P_{ST} = (8\sqrt3/15\cdot2)/2 = 8\sqrt3/15 = 8/(5\sqrt3)$ ($= 0.923760$) (sympy で確認、[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (e))。この章で QED として受け取るのは、上の $w_{\uparrow\downarrow}$、$w_{\downarrow\uparrow}$ の 2 つの数値係数 ($5\sqrt3/16$ と $8\sqrt3/15$) だけである。

### 7.2 K64 の定数が (15.17) と一致すること

K64 は $P_{ST}$、$c_{\rm fpd}$、$s_{\rm flc}$ を

$$
\mathtt{pst} = \frac{8\sqrt3}{15} = \frac{8}{5\sqrt3},\qquad \mathtt{cfpd} = \frac{5\sqrt3}{8},\qquad \mathtt{sflc} = \frac34\left(\frac{r_e}{\alpha}\right)^2 = \frac34\bar\lambda_C^2
\tag{15.18}
$$

と持つ ($r_e/\alpha = \bar\lambda_C$)。追跡 (`sprot`、`photons.f:1317-1340`) の Sokolov-Ternov の項は、1 スライスごとに、スピンを $\hat b_t$ の向き ($\vec b_t$ は運動量に垂直な磁場成分、磁場の向き) に次のように動かす。

```
st  = S·b_t / |b_t|                            (スピンの磁場方向の成分)
dst = (pst - st) * cfpd * (p*p0/h**2)**2 * ucde / |b_t|
S   = S + dst * b_t                            (成分ごと: sx = sx + dst*btx、...)
```

最後の行で $\vec b_t$ をかけるので、$|b_t|$ は約分され、磁場の向きの成分の増分は $|b_t|$ に依らず

$$
\Delta S_{\parallel b} = (\mathtt{pst} - s_t)\;\mathtt{cfpd}\;\Bigl(\frac{p\,p_0}{\gamma^2}\Bigr)^2\,\mathtt{ucde}
$$

である (磁場に沿う成分 $s_t$ が、到達値 `pst` に向かって、1 スライスごとに差 $(\mathtt{pst} - s_t)$ の一定割合ずつ近づく。この更新は 1 次の量で、スピンの長さの規格化は行わない。直後の回転は長さを変えない。)右辺の因子を 1 つずつ見る。(i) `p` は絶対運動量 ($p = p_0(1+\delta)$、$mc$ 単位)、`h` $= \gamma$ なので、$\delta = 0$ で $p\,p_0/\gamma^2 = p_0^2/\gamma^2 = \beta^2$ ($p_0 = \beta\gamma$)、その 2 乗が $\beta^4$。(ii) `ucde` $= u_c\,\delta E$ ($\mathtt{uc} = c_{\rm uc}\gamma^3\theta/(p_0\,al)$、$\mathtt{de} = c_{\rm ave}\,N_\gamma\,\mathtt{uc}$、$N_\gamma = \mathtt{anp} = a_{\rm nrad}\gamma\theta$。13 章) で、

$$
\mathtt{ucde} = c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2\;\frac{\gamma^7\theta^3}{p_0^2\,al^2}
$$

(代入: $c_{\rm ave}\cdot a_{\rm nrad}\gamma\theta\cdot(c_{\rm uc}\gamma^3\theta/(p_0\,al))^2$ の $\gamma$ のべきは $1 + 6 = 7$、$\theta$ のべきは $1 + 2 = 3$。) (i)(ii) を掛けて $p_0^2 = \beta^2\gamma^2$ を使うと、$\beta^4\cdot\gamma^7/(\beta^2\gamma^2) = \beta^2\gamma^5$。したがって、$s_t = 0$ のとき、1 スライスの磁場方向の増分は

$$
\Delta P = \mathtt{pst}\,\mathtt{cfpd}\,\beta^2\,(c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2)\,\frac{\gamma^5\,\theta^3}{al^2}\ \stackrel{\theta = al/\rho}{=}\ \mathtt{pst}\,\mathtt{cfpd}\,\beta^2\,(c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2)\,\frac{\gamma^5\,al}{\rho^3}
\tag{15.19}
$$

(偏向角 $\theta = al/\rho$ を入れて、$\theta^3/al^2 = al/\rho^3$。)

(15.17) の $dP = P_{ST}\,\tau_{ST}^{-1}\,dt$ ($s_t = 0$ で $P = 0$ から)、$c\,dt = al/\beta$ は $\mathtt{pst}\,\mathtt{cfpd}\,r_e\bar\lambda_C\,\gamma^5\,al/(\beta\rho^3)$。(15.19) と比べて、$\beta\to1$ (高エネルギー) で $\beta^2$ と $1/\beta$ の差は消えるので、一致するためには次が必要で、これが成り立つ:

$$
c_{\rm ave}\,a_{\rm nrad}\,c_{\rm uc}^2 = \frac{8}{15\sqrt3}\cdot\frac{5\alpha}{2\sqrt3}\cdot\frac94\bar\lambda_C^2 = \alpha\bar\lambda_C^2 = r_e\bar\lambda_C
\tag{15.20}
$$

($\tfrac{8\cdot5}{15\cdot2\cdot3}\cdot\tfrac94 = \tfrac49\cdot\tfrac94 = 1$。$c_{\rm uc} = \tfrac32\bar\lambda_C$、$a_{\rm nrad} = \tfrac{5}{2\sqrt3}\alpha$、$c_{\rm ave} = \tfrac{8}{15\sqrt3}$。) **確認** (c): 比 $1.000000000000$。つまり `sprot` の 1 ステップの項は ST の式そのものである ($P_{ST}\,\mathtt{cfpd} = 1$)。

**Emittance[] 側の `gintd`**: `tradke` の `gintd += sw(1:3)*sflc*anp*(bt*h1*p/al1)**2` は、1 周の積算 ($\mathtt{sw}$ は $\hat g$ を (srot の列で) リング入口まで戻した成分、$\mathtt{sdamp} = n_0\cdot\mathtt{gintd}$ は (各スライスの $\hat g\cdot n$) の和)。係数 $s_{\rm flc}\,a_{\rm nrad} = \tfrac34\bar\lambda_C^2\cdot\tfrac{5\alpha}{2\sqrt3} = \tfrac{5\sqrt3}{8}\,r_e\bar\lambda_C = c_{\rm fpd}\,r_e\bar\lambda_C$ ($\tfrac{15}{8\sqrt3} = \tfrac{5\sqrt3}{8}$) で、各スライスの項の中身は、$\mathtt{anp} = a_{\rm nrad}\gamma\theta$、$\mathtt{bt}\simeq\theta$ (磁場に垂直な成分の大きさ)、$\mathtt{h1} = \gamma$、$\mathtt{p} = p_0 = \beta\gamma$ ($\delta = 0$) から $\mathtt{anp}\,(\mathtt{bt}\,\mathtt{h1}\,\mathtt{p}/al)^2 = a_{\rm nrad}\gamma\theta\cdot(\theta\gamma\cdot\beta\gamma/al)^2 = a_{\rm nrad}\,\beta^2\gamma^5\theta^3/al^2$ ($\gamma$ のべき: $1 + 2 + 2 = 5$、$\theta$: $1 + 2 = 3$)。したがって $\mathtt{gintd}\cdot\hat n$ はスライスごとに $c_{\rm fpd}\,r_e\bar\lambda_C\,\gamma^5\beta^2\,\theta^3/al^2\cdot(\hat g\cdot\hat n)$ を足す。これは (15.17) の $\tau_{ST}^{-1}$ の経路長の積分 ($\hat g\cdot\hat n = 1$ のとき) に当たる。**確認** (c): $s_{\rm flc}a_{\rm nrad}/(c_{\rm fpd}r_e\bar\lambda_C) = 1.000000000000$。

**偏極時間**: $\tau = 1/(\mathtt{sdamp}\,f_{\rm rev})$ ($f_{\rm rev}$: 周回周波数。$\mathtt{sdamp}$ は 1 周あたり、$f_{\rm rev}$ は 1/s)。**確認** (e): 平面リング ($\rho = 4.5837$ m、偏向の長さ $14.4$ m、周長 $45.6$ m、$\langle\rho^{-3}\rangle = 3.2455\times10^{-3}$ m$^{-3}$、$\gamma = 5870.854$) で (15.17) の $\tau_{ST} = 123.8153$ s、Emittance[] の `PolarizationTime` $= 123.8153$ s (比 $1.000000$)。平衡偏極度 $P = 0.92376043$ (v2、$8/(5\sqrt3) = 0.923760$)。

### 7.3 SLIM の $\delta$ の項と DK の $11/18$

エネルギーの跳びの共分散 $B_{\delta\delta} = N_\gamma\,u_c^2\,c_{\rm uu}$ ($c_{\rm uu} = 11/27$。$N_\gamma = a_{\rm nrad}\gamma\theta$、$u_c = c_{\rm uc}\gamma^3\theta/(p_0\,al)$ から $B_{\delta\delta} = c_{\rm uu}\,a_{\rm nrad}\gamma\theta\,(c_{\rm uc}\gamma^3\theta/(p_0\,al))^2 = c_{\rm uu}a_{\rm nrad}c_{\rm uc}^2\gamma^7\theta^3/(p_0^2al^2) = c_{\rm uu}a_{\rm nrad}c_{\rm uc}^2\gamma^5\theta^3/(\beta^2al^2)$ ($p_0^2 = \beta^2\gamma^2$)) を (15.15) に入れる。スライスの $N$ の $\delta$ の列 ($\partial n/\partial\delta$) だけを取ると、(15.15) のスライスの項は $\tfrac12B_{\delta\delta}|\partial n/\partial\delta|^2$。$|\partial n/\partial\delta|^2$ の係数 $\tfrac12B_{\delta\delta}$ を、§7.2 の ST の 1 スライスの係数 $s_{\rm flc}\,a_{\rm nrad}\,\beta^2\gamma^5\theta^3/al^2$ ($\mathtt{gintd}$ の 1 項) で割ると ($\beta\to1$ で $\beta^{-2}/\beta^2\to1$)

$$
\frac12\,\frac{c_{\rm uu}\,a_{\rm nrad}c_{\rm uc}^2}{s_{\rm flc}\,a_{\rm nrad}} = \frac12\cdot\frac{\dfrac{55}{24\sqrt3}}{\dfrac{15}{8\sqrt3}} = \frac{11}{18}
\tag{15.21}
$$

(確認 (c): 数値で $0.611111111 = 11/18$。式の形: $c_{\rm uu}a_{\rm nrad}c_{\rm uc}^2 = \tfrac{11}{27}\cdot\tfrac{5\alpha}{2\sqrt3}\cdot\tfrac94\bar\lambda_C^2 = \tfrac{55}{24\sqrt3}\alpha\bar\lambda_C^2$、$s_{\rm flc}a_{\rm nrad} = \tfrac{5\sqrt3}{8}\alpha\bar\lambda_C^2 = \tfrac{15}{8\sqrt3}\alpha\bar\lambda_C^2$。比 $\tfrac12\cdot\tfrac{55}{24\sqrt3}\cdot\tfrac{8\sqrt3}{15} = \tfrac{55\cdot8}{2\cdot24\cdot15} = \tfrac{11}{18}$。) したがって、$\delta$ の跳びによる脱偏極率は $\tfrac{5\sqrt3}{8}r_e\bar\lambda_Cc\gamma^5\langle|\rho|^{-3}\cdot\tfrac{11}{18}|\partial n/\partial\delta|^2\rangle$ の形になり、DK の式 (§7.4) の $\tfrac{11}{18}|\gamma\,\partial n/\partial\gamma|^2$ と同じ係数を持つ。ここで $\gamma\,\partial/\partial\gamma$ と $\partial/\partial\delta$ の関係は、$p = p_0(1+\delta)$、$\gamma = \sqrt{1+p^2}$ から $d\gamma/d\delta = p_0\,p/\gamma$ なので $\gamma\,\partial/\partial\gamma = (\gamma/(d\gamma/d\delta))\,\partial/\partial\delta = (\gamma^2/(p_0p))\,\partial/\partial\delta$ で、$\delta = 0$ で $\gamma^2/p_0^2 = 1/\beta_0^2$、$\beta_0\to1$ で $\gamma\,\partial/\partial\gamma = \partial/\partial\delta$ (sympy で確認、[spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) (f))。SLIM の式 (15.15) は、DK の式の $\gamma\,\partial n/\partial\gamma$ の項を、全座標の $N$ に拡張した形である (水平・垂直の寄与は $B_{xx}$ などが入る)。

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
| `spin_01_thomas_bmt.py` | (a) 静止系の歳差 + Thomas 歳差 = BMT (sympy)、縦・横の係数、(b) らせん軌道の BMT の厳密な積分 (DOP853) と 1 ステップ回転の積 (2 次)、(b2) 平面の偏向で $a\gamma\theta$、縦磁場で $(1+a)\psi$、K64 の弦 $-\theta^2/24$、(c) 定数 $c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2 = r_e\bar\lambda_C$、$s_{\rm flc}a_{\rm nrad} = c_{\rm fpd}r_e\bar\lambda_C$、DK の 11/18、(d) (15.13) が 1 次で $R(x)n(x) = n(Mx)$ (代数の整合。方向を固定して残差が $|x_0|^2$ に比例)、(e) 平面リングの `SpinTune`・`PolarizationTime`・`Polarization` | [spin_01_thomas_bmt.out.txt](outputs/spin_01_thomas_bmt.out.txt) |
| `spin_02_derivations.py` | 本文の途中式の確認 (sympy・数値): (a) 磁場の Lorentz 変換 (15.1)、(b) Thomas 歳差の回転 (15.3a)、(c) Rodrigues の式 (15.10)、(d) (15.11) の恒等式、(e) スピン反転率から $P_{ST}$・$\tau_{ST}^{-1}$、(f) $\gamma\partial/\partial\gamma$ と $\partial/\partial\delta$、(g) $D$ の積算 (15.12) の有限差分との比較、(h) 要素ごとの $N$ (15.13') と 1 周の (15.13)、(i) ガウス広がりの和の積分 (15.23) と $e^{x^2}\mathrm{erfc}\,x$ の桁あふれ、(j) K64 の弦と $\mathtt{pzm}$、(k) 正規座標の $N$ (15.13'')、(l) 2 次の相関和と SLIM、(m) `spintrm` の $a\to0$ の極限の $2\pi$ | [spin_02_derivations.out.txt](outputs/spin_02_derivations.out.txt) |
| `spin_03_slim_mc.py` | SLIM の脱偏極率 (15.15) の、確率的な軌道 + 厳密なスピン回転の玩具のモンテカルロとの比較 | [spin_03_slim_mc.out.txt](outputs/spin_03_slim_mc.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/spin_01_thomas_bmt.py` (1 分弱。(e) の Emittance[] が大部分)、`spin_02_derivations.py` (数秒)、`spin_03_slim_mc.py` (1 分弱)。K64 との照合 (gs と 4e-12 以下、Emittance[] の偏極時間 14 桁など)、DK の数値積分との比較 (3 GeV〜4.2 GeV、追跡 FLUC) は、[docs/16_spin_pol.md](../docs/16_spin_pol.md) と `pysad/tests/script_tests/test_spin*.py`・`test_v2_*spin*.py`・`_p1_dk.py` に既にある。この章では新しく取り直していない。

## 12. 未確認・限界・使った条件

**使った条件**: 電子 (a = 1.159652181280002e-3)、3 GeV (γ = 5870.85)、平面リング (12 セル、偏向 15 度 × 12 = π。周長 45.6 m、$\rho = 4.5837$ m、CAVI 1.5 MV・harmonic 64、`RFSW; RAD; NOFLUC; POL`)。(b) のらせん軌道は計算上の一様磁場 ($\kappa = 1$ rad/m、角 0.3 rad)。

**未確認・導出できていないこと**:
- Sokolov-Ternov のスピン反転率 ($w_{\uparrow\downarrow}$、$w_{\downarrow\uparrow}$ の係数 $5\sqrt3/16$、$8\sqrt3/15$) と、DK の式 (15.22) の $\tfrac29$ の項・運動論的偏極の項の QED からの導出 (入力として使った。この章で確かめたのは、K64 の定数と式の間の一致 (15.20)・(15.21) と、反転率から (15.17) を組み立てる算術)。足りないもの: 磁場中の放射の確率のスピン依存の項 (Bessel 関数の積分) の計算。
- (15.12) の高次 ($\partial^2\vec g/\partial x^2$)、$|g_i|$ が大きい (指数写像の微分の補正 $J(g)$ が効く) 場合。K64 の $D$ は 1 次だけで、スライスを細かくして 1 次で収束する (玩具で確認、§5.1)。
- K64 の `spdepol`・`spintrm` の式全体: §6.5 で、実際に使われるのは 2 次の項だけであること、2 次の相関和が SLIM に一致すること (導出)、`a ≠ 0` の枝の $2\pi$ が $a\to0$ の極限で `a = 0` の枝と合わないこと、桁あふれの閾値を確かめた。残り: `rm(1)` の係数 (`0.25`、`0.5` など) の導出、`serot` の $\mathtt{se12}$・$\mathtt{de12}$、`ind` の添字の管理、B125 の「二重数え」が原因という推定の確定 (K64 の係数まで再導出するか、gs での確認が要る)。
- (15.15) のモンテカルロ確認は玩具の写像 (乱数の $D$) で行った。K64 の $D$・$N$ を持つ実リングの追跡 (FLUC つきのスピン追跡) との比較は docs/16 にあり、この章では取り直していない。
- 運動論的偏極・高次のスピン軌道共鳴 ($N$ の分母が 0 に近い) の扱い。分母の 0 の近くは (15.13) が使えない。
- 領域内 (SOL) の縦磁場、傾いた磁場での $B_\parallel$ の扱い (v2 の `bpar`)。
- 本章の数値は 1 つの平面リングで確認。結合のあるリング・非平面の偏極ベクトルは docs/16 の試験に頼る。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| POL・RADPOL・TrackParticles・Emittance[] の説明 | **原文で確認**: 「POL が ON のとき TrackParticles はスピン追跡を行う。RADPOL は Sokolov-Ternov を入れる」。座標は `{x, px/p0, y, py/p0, z, dp/p0, sarg, sy}`。Emittance[] のキー: `Polarization`、`Polarization2`・`4`・`6`、`PolarizationVector`、`SpinTune`、`NominalSpinTune`。 |
| `Polarization2`・`4`・`6` | ヘルプにだけある (2・4・6 次)。K64 のソースに文字列がなく、gs の `Emittance[]` も返さない (docs/16)。 |
| ヘルプに無いもの | Thomas-BMT、`sprot` の式、$D$、$N$、Sokolov-Ternov の定数、`spnorm` の手順、`Polarization` の式。(15.1)〜(15.23) はすべて K64 のソースと標準の理論から。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

各項目の根拠の程度を、docs/30 の記述に合わせて分けて書く。gs での確認が済んでいないものは、「疑い」または「ソースの読み」のままにしてある。

**独立な基準 (解析式・独立な数値計算・gs との比較) で確かめたもの**:
- **`Polarization` が NaN** (`spdepol`、B55): NaN になること自体は gs で再現済み。**原因は推定**: pysad で同じ式を写すと同じ NaN になり、`exp(cb²)·erfc(cb)` を桁あふれしない形で評価すると有限になる (docs/30)。本章の §6.5 は、倍精度での桁あふれの閾値 ($|c_b|\simeq26.6$) を確かめただけで、実際のリングの $c_b$ の値は未確認。
- **符号つきの `sdamp`** (B116): 逆向きの偏向で相殺する。独立な式 ($P = (8/5\sqrt3)\langle\rho^{-3}\rangle/\langle|\rho|^{-3}\rangle$、$\tau^{-1} = (5\sqrt3/8)r_e\bar\lambda_Cc\gamma^5\langle|\rho|^{-3}\rangle$) と K64 の追跡の緩和率の比で確認 (docs/30)。
- **弦 $2\sin(\theta/2)$** (B83): 刻みが粗いと歳差の角が $\theta^2/24$ だけ小さい。ソースの読みと解析式で、本章の確認 (b2)・§4.1 の (j) が一致する。docs/30 では、偏向磁石に `EPS=1E-5` を付けたときの gs と v2 の 1 周のスピン角が $10^{-8}$ で一致して、k64 の誤差は刻みを細かくすると消える。
- **磁石の本体の運動量の向きの磁場 $(1+a)B_\parallel\,ds$ が入らない** (B136・B165): K64 は、本体の中のスピンの歳差に $B_\parallel$ の積分を入れず、端の $\psi$ だけを入れる。傾いた軌道 ($p_y\ne0$) で厳密な値と食い違う。docs/30 では、スカラーポテンシャルから作った Maxwell 場の中の Thomas-BMT の数値積分と比べて確定している。v2 は本体の分を足す (切り替え `compat.set_spin_bpar`)。本章の (15.9) の $\vec b_l$ は、端の $\psi$ と、SOL・UND の縦磁場だけである。

**ソースの読み (gs での確認が部分的、または未確認)**:
- **B125**: `Polarization` の脱偏極率が、桁あふれを直しても、追跡・独立な式の約 1000 倍になる (docs/30 に数値)。**「二重数え」(位相が打ち消し合う 2 次の項の和 $1/(1-e^{-2a})$ が、すでに $1/(2a)$ を含む平衡エミッタンスにさらに掛かる、ガウスの項の $2\pi$ も余分) は推定**。docs/30 では、位置の特定は数値の依存性 (エミッタンスに比例、減衰率に反比例、`ssprd` に依らない) による。K64 の式の導出 (Mathematica 由来の展開) には戻っていない。本章 §6.5 の再導出は、2 次の相関和が (1 周の励起) で決まる SLIM に一致することと、$a\to0$ の極限で `a ≠ 0` の枝の $2\pi$ が `a = 0` の枝と合わないことで、この推定と整合するが、K64 の係数までは再導出していない。
- **B56**: `tbrote` (`tchge.f:197-202`)・`serot` (`photons.f:850-857`)・`tradke` の末尾 (`photons.f:1614-1615`) の `srot` の行を、更新後の値で逐次更新する。ソースの読み。gs での確認は `serot` と `tradke` の末尾だけで、`tbrote` は未確認 (`tradke` の末尾は平面リングの gs と一致したので意図かもしれない)。
- **B57**: `tradke` の `call tadddsg(dsg,srot,dgx,dgy,dgx,dbeam)` は `dgz` の引数に `dgx` を渡す (`photons.f:1589`)。`dsg` は使われない旧ルーチンでしか読まれないので影響なし。ソースの読み。
- **B58**: 偏向磁石の出口の `bsi2` の項が捨てられる。ソースの読みと、gs の追跡との比較 (出口の項だけを足すと SpinTune の差が $3.9\times10^{-4}\to2.7\times10^{-5}$。pysad を一時的に書き換えた試験で、恒久の変更ではない。docs/30)。

**疑い**:
- `tradke` の `dbtz = dxpzb/pr - dblz` (`photons.f:1550`): $x$ 成分・$y$ 成分は `dxpx/pr**2` だが、$z$ 成分だけ `pr` で割っている。$\delta\ne0$ での $\partial b_{tz}/\partial x$ の誤りになる (大きさは $\delta$ の程度)。誤植の可能性 (docs/16 は「写している」とだけ書く)。ソースの読みで、この章では数値で確かめていない。

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
