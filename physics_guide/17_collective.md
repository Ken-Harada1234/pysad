# 17. 集団効果: 空間電荷、ウェイク、CSR、Touschek・IBS

**範囲 (先に明記)**: この章は長いので、効果によって深さを変える。

| 効果 | 深さ | 節 |
|---|---|---|
| **空間電荷 (WSPAC)** | **導出する** (ローレンツ力 → $1/(\beta^2\gamma^3)$、ガウス分布のキック、$\delta$ のキック、線形のレンズ、数値で確認) | §3〜§4 |
| **ウェイク (LWAKE・TWAKE)** | **導出する** (ウェイク関数から $\delta$ と横のキック、自己項の 1/2、K64 の切片の方法、解析式で確認) | §5〜§6 |
| 空間電荷 (SPAC・PSPAC)、CSR、Touschek、IBS (INTRA) | 式とソースの対応表にとどめる (導出・数値確認は本章ではしていない) | §7 |

**ヘルプの状況**: ヘルプに式は無い。あるのは機能の説明だけ: `WSPAC` (「strong-weak モードの空間電荷。ビームの大きさは `EMIT` + `CODPLOT` で求める。すべての次元でガウス分布を仮定、バンチあたりの粒子数は `PBUNCH`」)、`LWAKE` (縦ウェイクつきの光学)、`INTRA` (バンチ内散乱)、`TouschekLifetime` の使い方、`MINCOUP` の式。式・単位は K64 のソースと標準の物理から起こした。

**前提**: [16_undulator_beambeam.md](16_undulator_beambeam.md) (2 次元ガウス分布の電場 = Bassetti-Erskine、(16.6)〜(16.8))、[14_emittance.md](14_emittance.md) (ビーム行列、MINCOUP)、[13_radiation.md](13_radiation.md)。

**読む順**: §1 (結論) → §3 → §5 → §8 (確認)。

## 1. 結論

- **空間電荷**: 長さ $L$ の区間で、粒子 ($z$、$(x, y)$) がバンチの電磁場から受けるキックは

$$
\Delta\vec p_\perp/p_0 = A(z)\,\nabla U(x,y),\qquad A(z) = \frac{N r_e L\,e^{-(z-z_c)^2/2\sigma_z^2}}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z},\qquad \Delta\delta = \frac{z-z_c}{\sigma_z^2}\,A(z)\,\bigl(U_{\rm ref} - U(x,y)\bigr)
$$

 ($U$ は 2 次元ガウス分布のポテンシャル (16.6)。式 (17.5)〜(17.7))。$1/(\beta^2\gamma^3)$ は、電場 ($\propto\gamma$ 倍) と磁場の引力の打ち消し ($1-\beta^2 = \gamma^{-2}$)、通過時間 $L/\beta c$、慣性 $\gamma$ から出る (§3)。K64 の `twspac` (追跡) の係数と一致し ([coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (b))、pysad の `WSpace.kick` は、ポテンシャルの数値積分から作った式と、横のキックで絶対差 $5\times10^{-12}$ 以下 (値 $\sim7\times10^{-6}$。内部変数 $p/(1+\delta)$ のため相対 $|\Delta\delta| \sim 7\times10^{-7}$ の差が残る)、$\delta$ のキックで $3\times10^{-16}$ 以下で一致。線形のレンズ $k_x = 2f/(\sigma_x(\sigma_x+\sigma_y))$ は式と $10^{-10}$ 以内。
- **ウェイク**: 縦ウェイク $W_L(s)$ [V/C] ($s$ = 試験粒子の前方の源までの距離) による $\delta$ の変化は $\Delta\delta(z) = -\dfrac{eN}{mc^2p_0}\displaystyle\int\lambda(z')\,W_L(z'-z)\,dz'$ (自己項は $W(0)$ の 1/2)。横ウェイク $W_T(s)$ [V/C/m] は源の位置 $x(z')$ に比例: $\Delta p_x = +\dfrac{eN}{mc^2p_0}\displaystyle\int\lambda(z')\,x(z')\,W_T(z'-z)\,dz'$。ガウス分布のバンチ + 指数関数のウェイク $W_0e^{-s/\ell}$ の平均の損失は $k = \tfrac12e^{\sigma^2/\ell^2}\mathrm{erfc}(\sigma/\ell)$ に比例し、`txwake` の直訳 (pysad の `WakeKick`) は $\ell = 2, 5, 20$ mm で解析式の $1.027$、$1.0115$、$1.0073$ 倍 (切片への粗視化の誤差。§6)。横ウェイクは $1.0115$ 倍。
- **対応表にとどめる効果** (§7): SPAC (3 次元 PIC)、PSPAC (2 次元 PIC)、CSR (縦のウェイクによる長さ・エネルギー幅の追跡と Haissinski 方程式)、Touschek・IBS (バンチ内散乱)。
- **K64 の疑い** (§12): WSPAC の縦の力の電位が壁の半径 $b$ に依る $\ln(b/\sigma)$ の項を持たない (基準は軸上で $\ln\frac{1+r}{2}$。未決)、WSPAC の Emittance[] の最後のスライスでレンズが 8 倍 (B127)、SPAC の力が $1/\gamma$ 倍 (B89)、など。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| $N$ | バンチあたりの粒子数 | `pbunch` (`PBUNCH`、既定 $10^{10}$) | `PBUNCH` |
| $r_e$ | 古典半径 | `rclassic` | `ELRADI` |
| $\gamma$、$\beta$ | ローレンツ因子、速度 | `h1 = p2h(p1)`、$p_1 = \beta\gamma$ (`p1`) | `_p2h` |
| $\sigma_z$ | バンチ長 | `sigzsq` = $\sigma_z^2$ (`beam(ia(5,5))`) | `sigz2` |
| $\sigma_x,\sigma_y$ | 横のビームの大きさ ($z$ を固定した条件付き。主軸) | `sigx`、`sigy` (`xx1`、`yy1`、`xy1` から) | `sig1`、`sig2` |
| $f$ | レンズの係数 $N r_eL/(\beta^2\gamma^3\sqrt{2\pi}\sigma_z)$ | `f` | `f` |
| $A(z)$ | $f\,e^{-(z-z_c)^2/2\sigma_z^2}$ | `az` | `az` |
| $U,\ \nabla U$ | 2 次元ガウス分布のポテンシャル・場 | `twspfu`: `fu`、`fx`、`fy` | `gauss_potential`、`gauss_field` |
| $W_L,\ W_T$ | 縦 [V/C]・横 [V/C/m] のウェイク関数の表 | `wakel(2,*)`、`waket(2,*)` | `wl`、`wt` |
| `fw` | $|{\rm charge}|\,e\,N\,n_{\rm bunch}/m c^2/n_{p}\times\tfrac12$ | `fw` | `kick_factor` |
| `ws(n)` | 切片 $n$ の粒子数 | `ws` | |
| `zs(n)` | 切片 $n$ の重心 | `zs` | |
| `wz(n)`、`wx(n)`、`wy(n)` | 切片 $n$ が受ける縦・横の和 | `wz`、`wx`、`wy` | |

## 3. 空間電荷の力

### 3.1 ローレンツ力

バンチの静止系 (記号に $'$) で、電荷分布は線密度 $\lambda' = \lambda/\gamma$ ($\lambda$: 実験室系の線密度、ローレンツ収縮) の静電場。横方向の電場 ($E'_\parallel$ は無視) は 2 次元ガウス分布の電場。丸い分布 ($\sigma_x = \sigma_y = \sigma$) なら、ガウスの法則から

$$
E'_r = \frac{mc^2}{e}\,\frac{2r_e\,\lambda'}{r}\left(1 - e^{-r^2/2\sigma^2}\right)
\tag{17.1}
$$

($e^2/(4\pi\varepsilon_0) = r_emc^2$)。実験室系に戻すと、横方向の電場は $E_\perp = \gamma E'_\perp$ (ローレンツ変換)、磁場は $B_\perp = (\beta/c)\gamma E'_\perp$ (バンチが速度 $\beta c$ で動く電流)。同じ速度で動く試験粒子が受ける力は

$$
F_\perp = e\,(E_\perp - \beta c\,B_\perp) = e\gamma E'(1 - \beta^2) = \frac{eE'}{\gamma}
\tag{17.2}
$$

(**確認** (a): sympy で $F\gamma/(eE') - 1 = 0$。電場の反発 ($\propto\gamma$) を磁場の引力 ($\propto\beta^2\gamma$) がほぼ打ち消す。) 長さ $L$ を通る時間は $\Delta t = L/\beta c$、試験粒子の運動量は $p_0 = \gamma m\beta c$。キックは

$$
\frac{\Delta p_\perp}{p_0} = \frac{F_\perp\Delta t}{p_0} = \frac{eE'}{\gamma}\,\frac{L}{\beta c\,\gamma m\beta c} = \frac{eE'L}{\gamma^2\beta^2\,mc^2}
\tag{17.3}
$$

(確認 (a): sympy で $(\Delta p/p_0)\,\gamma^2\beta^2mc^2/E' - 1 = 0$)。(17.1) の $E'$ と $\lambda' = \lambda/\gamma$ を入れて

$$
\frac{\Delta p_\perp}{p_0} = \frac{2r_e\lambda\,L}{\beta^2\gamma^3}\,\frac{\vec r}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)
\tag{17.4}
$$

**$1/(\beta^2\gamma^3)$ の由来**: $\gamma^{-2}$ (電場と磁場の打ち消し) $\times\ \gamma^{-1}$ (慣性) $\times\ \beta^{-2}$ ($\Delta t\propto1/\beta$、$p_0\propto\beta$)。バンチ内の位置 $z$ での線密度は $\lambda(z) = N e^{-z^2/2\sigma_z^2}/(\sqrt{2\pi}\sigma_z)$ なので、$A(z)$ は

$$
A(z) = \frac{Nr_eL}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z}\,e^{-(z-z_c)^2/2\sigma_z^2}
\tag{17.5}
$$

(丸い場合の (17.4) = $A(z)\,\nabla U$、$U = \mathrm{Ein}(r^2/2\sigma^2)$ ((16.6) の $c=1$ の丸い場合)。) 楕円の分布は、(16.6)〜(16.8) の Bassetti-Erskine のポテンシャル・場にそのまま置き換えればよい ($c \to N r_e L\,e^{-z^2/2\sigma_z^2}/(\beta^2\gamma^3\sqrt{2\pi}\sigma_z)$。ビーム・ビームでは電場と磁場が足し合って $c = N r_e/\gamma$ (対向)、ここでは打ち消し合って $1/(\beta^2\gamma^3)$ になる)。

### 3.2 強-弱モデルであること

WSPAC は、ビームの 2 次モーメント (ビーム行列) を固定したガウス分布が、1 個の粒子に与える力だけを計算する (ヘルプの「strong-weak mode」)。粒子どうしの力、ビーム行列の自己無撞着な変化は入らない。ビーム行列は `Emittance[]` の結果 (放射平衡) か、pysad では `WSPACBEAM = {εx, εy, σz, σδ}` で与える。要素の入口ごとに力を入れる (区間の長さ分を一括して)。

## 4. WSPAC のキックと行列

### 4.1 $z$ を固定した横の分散

ビーム行列 $\Sigma$ (6×6) から、バンチ内の位置 $z$ を固定した横の 2 次モーメント (条件付き) を作る: $\sigma_{xx}' = \Sigma_{11} - \Sigma_{15}^2/\Sigma_{55}$、$\sigma_{yy}' = \Sigma_{33} - \Sigma_{35}^2/\Sigma_{55}$、$\sigma_{xy}' = \Sigma_{13} - \Sigma_{15}\Sigma_{35}/\Sigma_{55}$ (`xx1`、`yy1`、`xy1`)。分散で $x$-$z$ の相関 (分散関数・クロマティックな項) の影響を除く。主軸の角 $\theta = \tfrac12\arctan2(2\sigma_{xy}', \sigma_{xx}' - \sigma_{yy}')$ ($c_1 = \cos\theta$、$s_1 = \sin\theta$)、主軸の大きさ $\sigma_{x,y} = \sqrt{\tfrac12(\sigma'_{xx} + \sigma'_{yy}) \pm \tfrac12\sqrt{(\sigma'_{xx} - \sigma'_{yy})^2 + 4\sigma_{xy}'^2}}$。16 章 §7.2 の「$z$ を固定した条件付き包絡」と同じ操作 (精度行列のブロックの逆行列と同じ結果: ガウス分布の条件付き分散)。

### 4.2 追跡のキック (`twspac`)

位置 $(x, y, z)$ の粒子について、主軸の座標 $(x_1, y_1) = R(\theta)(x - x_c, y - y_c)$ で

$$
\Delta p_{x_1,y_1} = A(z)\,\nabla_{(x_1,y_1)}U(x_1,y_1;\sigma_x,\sigma_y),\qquad
\Delta\delta = -\frac{z - z_c}{\sigma_z^2}\,A(z)\,U(x_1,y_1)
\tag{17.6}
$$

(`fx = -az*fx` の符号は `twspfu` の戻り値の符号の約束。`dg = -(z-zc)/sigzsq*az*fu`。$U$ は中心で 0 ($U(0,0) = 0$) にとる。)

**$\delta$ のキックの導出**: 空間電荷のポテンシャルエネルギー (長さ $L$ の分) は、力 $+A\nabla U$ (反発) の原始関数として $V = -A(z)\,(U - U_{\rm ref})$ ($U_{\rm ref}$: ポテンシャルの基準)。SAD の $(z, \delta)$ は $(q, p)$ の対で、$\delta' = -\partial V/\partial z = A'(z)\,(U - U_{\rm ref})$ ($A' = -\frac{z-z_c}{\sigma_z^2}A$) なので

$$
\Delta\delta = \frac{z - z_c}{\sigma_z^2}\,A(z)\,\bigl(U_{\rm ref} - U(x_1, y_1)\bigr)
\tag{17.7}
$$

K64 の `fu` (`twspfu`) は軸上で 0 ではなく、$-U_{\rm ref} = \ln\frac{1+r}{2}$ ($r = \sigma_2/\sigma_1$。`twspu` の $\log(0.5+0.5r)$ の項。遠方で丸いビームの形に合う基準): $\Delta\delta = -\frac{z-z_c}{\sigma_z^2}A\,[U + \ln\frac{1+r}{2}]$ ((17.6) の第 2 式の $U$ を $U + \ln\frac{1+r}{2}$ にしたもの)。丸いビーム ($r = 1$) では軸上 0。**gs で確認**: 0.5 GeV の結合のあるリング ($r \approx 0.06$) の軸上の粒子も縦のキックを受け ($z = \sigma_z$ で $+9.7\times10^{-10}$)、$-\frac{z}{\sigma_z^2}A\ln\frac{1+r}{2}$ と 1% で合う。pysad (v2) は以前この項を落としていた (軸上 0)。pysad (v2) は `WSPACRADIUS` ($> 0$) で、基準をダクトの半径 $b$ の電位 ($U_{\rm ref} = U(b)$) にする (物理的には、縦の電場はダクト壁の半径に依存する: $E_z\propto-\lambda'(z)\,[2\ln(b/r)+\ldots]$。K64 は壁なしで、軸上の縦の力は丸いビームで 0。§12)。**確認** (b): pysad の $\Delta\delta$ は、(17.7) を $U_{\rm ref} = 0$ として数値積分の $U$ から作った値と、**絶対差 $3\times10^{-16}$ 以下** (値 $\sim10^{-6}$) で一致する。

**運動量の更新** (`twspac` の末尾): $p_r = 1+\delta_{\rm old}$、$\delta_{\rm new} = \delta_{\rm old} + \Delta\delta$、$p_x \leftarrow (p_xp_r + \Delta p_x)/(1+\delta_{\rm new})$ (内部変数 $p_x/(1+\delta)$ ⇒ 正準運動量に戻して足し、新しい $\delta$ で割る)、$dv$ を更新。

### 4.3 線形のレンズ (`twspace`)

小振幅 ($r \ll \sigma$) で、(17.4)・(16.8) より $\nabla U \to 2(x_1/(\sigma_x(\sigma_x+\sigma_y)),\ y_1/(\sigma_y(\sigma_x+\sigma_y)))$ なので、中心 ($z = z_c$) のレンズは

$$
k_x = \frac{f'}{\sigma_x},\quad k_y = \frac{f'}{\sigma_y},\qquad f' = \frac{2Nr_eL}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z\,(\sigma_x+\sigma_y)}
\tag{17.8}
$$

(`akx = f/sx`、`aky = f/sy`、`f = 2 pbunch rclassic al/(p1**2 h1 (sx+sy) sqrt(2π sigzsq))`)。主軸が傾いていれば `trans(2,1) = aks + akd c2` などで回す ($aks = (k_x+k_y)/2$、$akd = (k_x-k_y)/2$、$c_2 = \cos2\theta$、$s_2 = \sin2\theta$)。**確認** (b): `WSpace.linear_k` は (17.8) と $10^{-10}$ 以内。

**Laslett のチューンシフト**: 各要素の入口にこのレンズを入れると、小振幅・バンチ中心の粒子のチューンシフトは $\Delta Q_{x} = -\dfrac{1}{4\pi}\sum\beta_xk_x$ ($\beta_x$: Twiss)。連続分布 ($L\to ds$) の極限は

$$
\Delta Q_x = -\frac{r_e\,\lambda_0}{2\pi\beta^2\gamma^3}\oint\frac{\beta_x\,ds}{\sigma_x(\sigma_x+\sigma_y)},\qquad \lambda_0 = \frac{N}{\sqrt{2\pi}\sigma_z}
\tag{17.9}
$$

(標準の Laslett の式。ガウス分布のバンチの中心の粒子)。pysad の WSPAC の Emittance[] のチューンシフトが (17.9) と 1〜5% で一致することは [docs/10_tracking.md](../docs/10_tracking.md) §9.7 (c) と docs/09 に記載 (本章では取り直していない)。K64 は、一周の最後のスライスで `gammab(nlat+1)` が配列の外で 0 になり、$p_1$ が半分、レンズが 8 倍強い (B127)。

## 5. ウェイク

### 5.1 ウェイク関数の定義と 1 個の源のキック

電荷 $q$ の点源が $z_s$ にあるとき、$z$ ($< z_s$、後方) の試験粒子 (電荷 $e$) が受けるエネルギーの変化は、ウェイク関数 $W_L$ [V/C] で

$$
\Delta E(z) = -e\,q\,W_L(z_s - z)
\tag{17.10}
$$

($W_L > 0$ は損失。SAD の約束は `dz = zs(m) - zs(n)` = 源 ($m$) の位置 $-$ 試験粒子 ($n$) の位置。$z > 0$ が先頭なので、因果的なウェイクは $dz > 0$ だけで 0 でない。) 横方向の力は、源の横のずれ $x_s$ に比例する (双極子ウェイク)。横ウェイク $W_T$ [V/C/m] で $\Delta p_xc = +e\,q\,W_T(z_s - z)\,x_s$ (これも $W_T > 0$ が「後続を源のずれの向きに引く」符号の約束)。

### 5.2 バンチ全体のキック

$N$ 個の粒子のバンチ (密度 $\lambda(z)$、$\int\lambda = 1$) は、各源のキックの和:

$$
\Delta\delta(z) = -\frac{eN}{mc^2\,p_0}\int\lambda(z')\,W_L(z'-z)\,dz',\qquad
\Delta p_x(z) = +\frac{eN}{mc^2\,p_0}\int\lambda(z')\,x(z')\,W_T(z'-z)\,dz'
\tag{17.11}
$$

($\Delta E/mc^2$ を $p_0$ ($= \beta\gamma$) で割ると、$\Delta\delta = \Delta p/p_0 = \Delta E/(\beta^2\gamma mc^2)$ なので、K64 の式は $\beta\to1$ の近似。$e$ は C、$W$ は V/C、$mc^2$ は eV。)

**自己項**: $z' = z$ の項 (粒子自身の源) は、後ろに寄与しないので、**$W_L(0)$ の 1/2** (beam loading の基本定理)。K64 は、同じ切片の内部 ($n = m$) で、表の最初の点が $z \ge 0$ のとき $\tfrac12W(z_1)\,ws(n)$ を使う。

**解析式** (確認 (c) の基準): $W_L(s) = W_0e^{-s/\ell}$ ($s \ge 0$) とガウス分布 ($\sigma_z = \sigma$) の平均の損失。2 粒子の距離 $s = z'-z$ は分散 $2\sigma^2$ のガウス分布、$s > 0$ の寄与だけ (源が先行):

$$
\langle W\rangle = W_0\int_0^\infty e^{-s/\ell}\,\frac{e^{-s^2/4\sigma^2}}{2\sigma\sqrt\pi}\,ds = \frac{W_0}{2}\,e^{\sigma^2/\ell^2}\,\mathrm{erfc}\!\left(\frac\sigma\ell\right)\equiv W_0\,k
\tag{17.12}
$$

($\int_0^\infty e^{-as-bs^2}ds = \tfrac12\sqrt{\pi/b}\,e^{a^2/4b}\mathrm{erfc}(a/2\sqrt b)$、$a = 1/\ell$、$b = 1/4\sigma^2$。$\sigma\to0$ で $k\to\tfrac12$ (自己項の 1/2 と一致)。) したがって平均の損失は $\langle\Delta\delta\rangle = -\dfrac{eNW_0k}{mc^2\,p_0}$ (横ウェイクで全体を $x_0$ ずらしたときの平均キックは $\langle\Delta p_x\rangle = +\dfrac{eNW_{0t}x_0k}{mc^2\,p_0}$)。

### 5.3 K64 の方法 (`txwake`): 切片と補間

粒子を $z$ で整列し (`spsort`)、$z$ の隣り合う間隔が閾値 $s_b$ ($\mathtt{nb} > 1$ のとき $(z_{\max}-z_{\min})/(2(n_b-1))$、そうでなければ全幅) を超えたところでバンチの境界とする (`txwdefslice`)。各バンチを、$\sqrt{n_m}$ 個ずつの粒子の**切片**に分ける ($n_m$ はバンチの粒子数。切片の数 $\approx\sqrt{n_m}$)。切片 $n$ の粒子数 `ws(n)`、座標の和 `xs(n)`、`ys(n)`、重心 `zs(n)`。切片の組 ($n$, $m$) について

$$
\mathtt{wz}(n) = \sum_mW_L(z_m - z_n)\,\mathtt{ws}(m),\qquad \mathtt{wx}(n) = \sum_mW_T(z_m - z_n)\,\mathtt{xs}(m)
\tag{17.13}
$$

($W$ は表の線形補間。$\mathtt{xs}(m)$ は切片の $x$ の**和**で、切片内の平均 $\times$ `ws`。(17.11) の $\int\lambda\,x\,W$ に当たる。) 各粒子 ($k$) のキックは、切片の重心の間で**線形補間**: $dwz = (1-w)\mathtt{wz}(n) + w\,\mathtt{wz}(n+1)$、$w = (z_k - z_n)/(z_{n+1} - z_n)$。

$$
\delta_k \leftarrow \delta_k - \mathtt{fwp}\cdot dwz,\qquad p_{x,k} \leftarrow \frac{p_a\,p_{x,k} + \mathtt{fwp}\cdot dwx}{p_b},\quad \mathtt{fwp} = \mathtt{fw}/p_0,\quad \mathtt{fw} = \frac{|{\rm charge}|\,e\,N\,n_{\rm bunch}}{mc^2\,n_p}\cdot\frac12
\tag{17.14}
$$

($p_a = 1+\delta_{\rm old}$、$p_b = 1+\delta_{\rm new}$。$n_p$ は `TrackParticles` に渡した粒子数 ($\mathtt{ws}$ の和 = マクロ粒子の数) で、$1/n_p$ が $\lambda$ の規格化。因子 $\tfrac12$ は、要素の前と後で 1 回ずつ同じ大きさのキックを入れるため (合計で (17.11) の 1 回分)。) **確認** (c): この方法を直訳した pysad の `WakeKick` (20000 粒子、$\sigma_z = 5$ mm、$N = 10^{10}$、1.7 GeV、$W_0 = 10^{12}$ V/C) の $\langle\Delta\delta\rangle$ は、(17.12) の解析式の $1.027$ ($\ell = 2$ mm)、$1.0115$ (5 mm)、$1.0073$ (20 mm) 倍。符号 (損失)・規格化 $eN/(mc^2p_0)$・自己項の 1/2 が合う。差 (+0.7〜2.7%) は、切片への粗視化: 同じ切片の内部の対は距離 0 とされ、表の $W(0)$ (全重み) が使われる (実際は半分が先行、半分が後続で平均 $W_0/2$)。切片の数 $\approx\sqrt{20000} = 141$ で、この過大評価は $\tfrac12\cdot\tfrac1{141}\cdot\tfrac{1}{k}$ ($= 3.3\%$、$1.7\%$、$0.9\%$ ($\ell$ = 2、5、20 mm))。観測値 (2.7%、1.15%、0.73%) と同程度だが、定量的には分離していない。

**光学 (`WAKEOPT`)**: `InitialOrbits` で与えた複数の軌道が互いにキックを受けるとして、軌道ごとの Twiss を作る (`tffswake.f`)。TRPT + INS の前提。本章では扱わない ([docs/13_wake.md](../docs/13_wake.md))。

## 6. ウェイクの確認

[coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (c): 縦ウェイクの規格化・符号・自己項は上のとおり。横ウェイクは、バンチ全体を $x_0 = 1$ mm ずらして $W_T = 10^{13}e^{-s/5\,\mathrm{mm}}$ V/C/m を入れたとき、$\langle\Delta p_x\rangle = +2.038\times10^{-9}$ (解析式の $1.0115$ 倍、縦と同じ粗視化の差)、$\Delta\delta = 0$ (縦ウェイクなし)。

## 7. 式とソースの対応表 (導出・数値確認は本章ではしていない)

| 効果 | 物理 (標準の式) | K64 | 状態 |
|---|---|---|---|
| **SPAC** (3 次元) | バンチの静止系の 3 次元 Poisson 方程式を、半径 $R$ の円筒の導体壁で解き、ローレンツ変換して実験室系に戻す。キック $\Delta p_\perp/p_0 = -C\partial\phi'/\partial(x,y)$、$C = r_cwL/(\beta_0^2\gamma_0^2)$、$\Delta\gamma = -r_cwL\,\partial\phi'/\partial\zeta$。一様な長い円柱で $x'' = Kx/a^2$、$K = 2r_c\lambda/(\beta^2\gamma^3)$ (ペルビアンス) | `spkick.f` (`spdrift_free`、`spfield`)、`tsol.f` (`spdrift_solenoid`)。DRIFT の `RADIUS > 0` | pysad v2: `spch.Solver3D`。**K64 は電場が $1/\gamma$ 倍に弱い (B89)** ほか (B90〜B92、B98) |
| **PSPAC** (2 次元 PIC) | 横 2 次元の Poisson 方程式 (矩形メッシュ、境界で電位 0) の電位 $\Phi(x,y)$ と縦の線密度 $\lambda(z)$ の積。$\Delta p_x/p_0 = -\frac{r_cwL}{\beta^2\gamma^3}\lambda(z)\partial\Phi/\partial x/N$ | `tspch_.c` (`tpspac`)、`spch.c`。`SPCH` 要素 | pysad v2: `spch.Solver2D`。K64 は既定 (`PSPACNZ = 1`) で落ちる (B93)、周の終わりのキックが入らない (B95) |
| **WSPAC** | 本章 §3〜§4 | `tintrb.f` (`twspac`、`twspace`、`qwspac`、`twspfu`、`twspu`) | 導出・確認済み。B127、B133、縦の力の基準 (§12) |
| **ウェイク** | 本章 §5 | `txwake.f`、`tffswake.f`、`tfwake.f` | 導出・確認済み |
| **CSR** (縦) | 1 周を `ndiv` 個の区間に分け、(i) $z \mathrel{+}= (r_{56}/n_{\rm div})\,\delta$、(ii) 密度 $\rho(z)$ (`csrrho`) とウェイクの畳み込み $w(z) = \int W(z-z')\rho(z')dz'$ (`csrwake`、インピーダンス `zl` を FFT)、(iii) $\delta \mathrel{+}= (aw/n_{\rm div})\,w(z)$ + RF (線形 $r_{65}z$、または $v_r\sin(\omega z+\phi_{\rm rf}) - dp_r$)、(iv) 放射減衰と量子励起 $\delta \to a_d\delta + b_d\,\xi$ ($a_d = 1 - 2d/n$、$b_d = \sqrt{4d/n\,(1-d/n)}\,\sigma_E$) を繰り返す。$b_d$ は、励起のない平衡で $\mathrm{Var}(\delta) = b_d^2/(1 - a_d^2) = \sigma_E^2$ になる値 | `csrtrack.f` (`csrtrack`、`csrwake`、`csrrho`)、`csroy.f` (Oide-Yokoya の安定性)、`Packages/CSR.n` (インピーダンス) | pysad: `csr.py`。K64 と照合済み ([docs/15_csr.md](../docs/15_csr.md))。ウェイクの導出と数値確認はしていない |
| **Haissinski** | ウェイクのある定常分布 $\rho(z) = K\exp\!\left[-\frac{z^2}{2\sigma_z^2} - aw\,V(z)\right]$、$V(z)$ = ウェイクの畳み込みの原始関数 (`csrwake` + `zi = zl/ak`)。Newton 法 (`tsolvg`) | `csrtrack.f` (`csrhaissin`) | 同上 |
| **IBS** (INTRA) | バンチの静止系の 3 次元 $(x, y, z)$ のビーム行列の固有値 $(e_1,e_2,e_3)$ と角度積分 $\int f(t)\,dt$、クーロン対数 $\ln(2b_{\max}/b_{\min})$ ($b_{\max} = \min(\sigma_{\min}, (V/N)^{1/3})$、$b_{\min} = \max(r_c/(p_t p_0)^2, \sqrt{V/\pi/(p_tp_0c/\gamma)/N/\tau_{\max}})$) から拡散行列 `bmi` を作り、励起 $\Delta$ に足す (Bjorken-Mtingwa 型) | `tintrb.f` (`tintrb`、`fintrb`、`fintrbu`) | pysad: `touschek.py`。MINCOUP の下限は `tintraconv` ([14_emittance.md](14_emittance.md) §5)。導出・数値確認なし |
| **Touschek** | 衝突の相対運動量が運動量アクセプタンスを超える散乱の損失率。K64 は平坦ビームの Bruck-Piwinski の関数の当てはめ `touckf(x)` $\propto\ln(1/x/\varepsilon_E + 1)e^{-x}(b + \varepsilon_Ex)/(b + x(a + 2x))/x$ ($a = 5.62966$、$b = 0.75159$、$\varepsilon_E = 7.98221$) を使い、損失率の表 `touckl`・`touckm`・`toucke` を積算 | `tintrb.f` (`touckf`、`tintrb` の `caltouck`)、`Packages/emit.n` (`TouschekLifetime`) | pysad v2: Piwinski 型の式が既定 (B162)。[docs/26](../docs/26_frequency_map_lma.md)。導出・数値確認なし |

## 8. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `coll_01_spacecharge_wake.py` | (a) 実験室系の力 $= eE'/\gamma$ (sympy)、(b) WSPAC のキック (横・$\delta$・線形のレンズ) と 2 次元ポテンシャルの数値積分、(c) ウェイク: ガウス + 指数関数の平均の損失 (縦)、全体のずれの平均キック (横) | [coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) |

実行: リポジトリのルートから `py -3.12 -W ignore physics_guide/scripts/coll_01_spacecharge_wake.py` (数秒)。

## 9. Fortran の対応

| 式 | ファイル:行 | サブルーチン | 変数・内容 |
|---|---|---|---|
| (17.5)〜(17.7) | `tintrb.f:378-455` | `twspac` | `f`、`az`、`twspfu`、`fx`、`fy`、`fu`、`dg`、`dpx`、`dpy`、`pr1`、`h1`、`dv` |
| §4.1 | `tintrb.f:395-407` | `twspac`、`twspace` | `xx1`、`yy1`、`xy1`、`a`、`c1`、`s1`、`sigx`、`sigy` |
| (17.8) | `tintrb.f:281-337` | `twspace` | `f`、`akx`、`aky`、`aks`、`akd`、`trans(2,1)` など |
| 光学への挿入 | `tintrb.f:339-376` | `qwspac` | `trs`、`coup` |
| (17.11)〜(17.14) | `txwake.f:5-170` | `txwake` | `ws`、`xs`、`ys`、`zs`、`wz`、`wx`、`wy`、`fwp`、`dwz`、`dwx`、`dwy` |
| 切片の定義 | `txwake.f:171-243` | `txwdefslice` | `izs`、`itab`、`mb` |
| `fw` | `tturn.f:259` | `tturn1` | `fw = (|charge| e pbunch anbunch/amass)/np0*.5` |
| CSR の追跡 | `csrtrack.f:229-407` | `csrtrack` | `ndiv`、`r56n`、`ad`、`bd`、`awn`、`rho`、`ddrho` |
| Haissinski | `csrtrack.f:131-227` | `csrhaissin` | `u`、`u0`、`f`、`g`、`wi` |
| IBS・Touschek | `tintrb.f:24-246`、`:6-19` | `tintrb`、`touckf` | `vol`、`eig`、`bmax`、`bmin`、`ci`、`touckl`、`touckm`、`toucke` |

(行番号は参照時点のもの。)

## 10. pysad の対応

| 関数 | 対応する Fortran | 注意 |
|---|---|---|
| `spch.WSpace` (`kick`、`linear_k`) | `twspac`、`twspace` | 縦の力の基準は K64 の $\ln\frac{1+r}{2}$ (軸上)。v2 は縦の力の基準 `rwall` (`WSPACRADIUS`) とバンチの速度での位置 `dz` を追加 (K64 は粒子の $z$ そのまま。差は $\delta/\gamma^2$)。`selfcod` |
| `spch.gauss_field`、`gauss_potential` | `twspfu`、`twspu` | Bassetti-Erskine の複素誤差関数と数値積分 |
| `spch.Solver3D`・`Solver2D` | `spkick.f`・`tspch_.c`+`spch.c` | v2 で正しい係数に書き直し (K64 の誤りは §12) |
| `wake.WakeKick`、`wake._kick` | `txwake` | K64 の癖をそのまま残す (§12)。切片の方法も直訳 |
| `wake.kick_factor` | `tturn.f:259` | $\mathtt{fw}$ |
| `csr.py` | `csrtrack.f`、`csroy.f` | K64 と照合済み |
| `touschek.py` | `tintrb.f` (IBS・Touschek) | Piwinski 型 (v2 の既定) と K64 の `touckf` (k64) |

## 11. 使った条件

- WSPAC (確認 (b)): 20 MeV の電子 ($\gamma = 40.139$、$\beta = 0.99969$)、$N = 10^{10}$、$\sigma_z = 5$ mm、区間の長さ 0.3 m、$\sigma = 1$ mm (丸い)、$\sigma_x = 1.5$ mm・$\sigma_y = 0.5$ mm (楕円)、$z$ は $0$、$\pm3$〜$6$ mm、$x, y$ は $0.4$〜$2.5$ mm。壁の半径は 0 (K64 と同じ基準: 軸上で $\ln\frac{1+r}{2}$、丸いビームでは 0)。
- ウェイク (確認 (c)): 1.7 GeV の電子、20000 個のガウス分布 ($\sigma_z = 5$ mm、種 5)、$N = 10^{10}$、$W_0 = 10^{12}$ V/C ($W_{0t} = 10^{13}$ V/C/m)、$\ell = 2, 5, 20$ mm。表は 4001 点 (0.05 mm 刻み、0〜0.2 m)。バンチ数 `nb` = 1。

## 12. 未確認・限界・K64 の疑い

**未確認**:
- SPAC・PSPAC・CSR・Haissinski・IBS・Touschek の導出と数値確認 (§7 は式とソースの対応のみ)。
- WSPAC: チューンシフトの Laslett の式との一致 (docs に記載の 1〜5% を取り直していない)、$\sigma_z$ や $\delta$ が大きいときの粗視化、ビーム行列の自己無撞着な収束 (K64 は `Emittance[]` の平衡値で、空間電荷を入れた平衡は反復で求める。docs/09)。
- ウェイク: 光学 (`WAKEOPT`)、切片の方法の誤差の完全な定量 (粗視化の差の見積りは大きさだけの一致)、複数バンチ (`nb > 1`)、実機のウェイク関数 (確認は人工的な指数関数だけ)。
- (17.11) の $\beta\to1$ の近似 ($\Delta\delta = \Delta E/(mc^2p_0)$、正しくは $/(\beta p_0)$)。

**疑い (K64。ソースの読みと一部の数値)**:
- **WSPAC の縦の力の電位の基準** (`fu` は中心で $\ln\frac{1+r}{2}$。丸いビームで 0。`tintrb.f:440` の `dg = -(z-zc)/sigzsq*az*fu`。pysad v2 は以前、軸上 0 としていたのを K64 の基準に直した: gs で軸上の粒子の縦のキックが $-\frac{z}{\sigma_z^2}A\ln\frac{1+r}{2}$ と 1% で合う): 縦の空間電荷の主要項 ($\ln(b/\sigma)$ に比例する、壁の半径に依る項) が入らない。pysad v2 は `WSPACRADIUS` で基準を壁の半径にできる。**仕様か誤りか未決 (コードは変えない)**。理由: (1) K64 の WSPAC はダクトの半径を持たない「壁のない自由空間」の 2 次元ガウス分布のモデルで、自由空間の線電荷のポテンシャルは対数発散し、基準 (定数) が定義できない。軸の電位に $\ln\frac{1+r}{2}$ を取るのは、その定数を決める 1 つの選択で、符号や式の誤りではない (横の位置に依る部分 $U(x_1, y_1) - U(0, 0)$ は、2 次元ポテンシャルの数値積分と絶対差 $3\times10^{-16}$ 以下で一致する: 確認 (b))。(2) ただし、落ちる定数は大きい: 丸いガウス分布 (rms $\sigma$) で、軸の電位は壁 (半径 $b$) に対して $(\lambda/4\pi\varepsilon_0)\,[\ln(b^2/2\sigma^2) + \gamma_E + E_1(b^2/2\sigma^2)]$ ($\gamma_E$: オイラー定数。私の導出を数値積分で確認)、軸から $\sigma$ までの電位差の $(\lambda/4\pi\varepsilon_0)\times0.444$ に対し、$b/\sigma = 10, 30, 100$ で 10 倍・15 倍・20 倍。したがって、壁のあるリング・輸送路の縦の空間電荷 (縦のインピーダンスの空間電荷の項、ポテンシャルウェルのゆがみ) は WSPAC (k64) に入らない。(3) K64 の作者の意図 (自由空間の横方向の力だけを使う、縦は補助) を示すものはヘルプにもソースのコメントにも無い。(4) gs で、軸上 ($x = y = 0$) の粒子の縦のキックを測った (台本の `Do` の中の `WSPAC;` は効かないので、`Do` を使わず上から並べる): 軸上の粒子も縦のキックを受け、基準が軸上 0 でないことが分かった (上記)。壁の半径を変える実験は、K64 に壁のパラメータが無いので gs ではできない。pysad の k64 の WSPAC の追跡は未実装 (v2 が実装)。結論: **未決**。解析は私の導出 (数値積分で確認済みの公式) の範囲。
- **B127**: WSPAC の `twspace` が最後のスライスで `gammab(nlat+1) = 0` を読み、レンズが 8 倍強い (docs/09 に記載、v2 は直した)。
- **B133**: WSPAC ON の CALC・CALC6D は入口以外の Twiss が NaN (gs で再現、原因は未調査。v2 は CALC6D と RFSW ON の CALC に空間電荷を入れた)。
- **B89〜B92・B98**: SPAC の力が $1/\gamma$ 倍に弱い (電場の $\gamma$ の抜け)、キック後の $p_x$ の二重の割り算、MULT の `RADIUS` の呼び出しの引数 (B91)、メッシュの整数の割り算、少数の粒子で落ちる (docs/10、docs/30。v2 は書き直し)。
- **B93〜B96**: PSPAC の既定で落ちる、線密度の配分、周の終わりのキックの欠落、static 変数 (docs/30)。
- **B97**: `tspac.f` の `tspac` は死んだコードで、係数 $1/(p_0(1+p^2)) = 1/(\beta\gamma\cdot\gamma^2)$ は $\beta$ が 1 つ足りない (呼ばれない)。
- **ウェイクの癖** (docs/13): (i) 追跡の表の最後の 1 行が使われない (`twxiwp` の `lwl0` は行数 − 1)、(ii) `tturn.f` で `dy = ROTATE` の代入、(iii) CAVI・MULT の分割で、wake の蹴りの合計が $(1.5n-0.5)/n$ 倍 (v2 は台形の重みで合計 1)、(iv) `txwdefslice` の配列の外への書き込みの可能性。
- **B162**: Touschek の損失率の式 `touckf` は平坦ビームの当てはめで、適用範囲に限界 (誤りというより近似。v2 は Piwinski 型が既定)。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| `WSPAC` のフラグの説明 | **原文で確認**: 「strong-weak モード。ビームの大きさは `EMIT` + `CODPLOT` で。ガウス分布を仮定、`PBUNCH`」。`SELFCOD` の説明 (力の中心を粒子の重心にするか閉軌道にするか)。`PSPAC` は PIC で、`SPAC`・`WSPAC` と別物 (ヘルプの注意書き)。 |
| `LWAKE`・`TWAKE`・`WAKEOPT` | **原文で確認**: 縦ウェイクつきの光学、TRPT の下で WAKEOPT が有効。`WakeFunction` は TRPT + INS + (TWAKE or LWAKE) のときだけ有効、追跡は `TrackParticles` だけ。 |
| `INTRA`・`MINCOUP`・`TouschekLifetime` | **原文で確認**: `MINCOUP` の式 $\varepsilon_{x,y} = \max[\varepsilon^0_{x,y}, \mathtt{MINCOUP}(\varepsilon^0_x+\varepsilon^0_y)]$、`TouschekLifetime[nx, ny, nz]` の 3 通りの使い方。 |
| ヘルプに無いもの | (17.1)〜(17.14)、IBS・Touschek・CSR の式。ウェイク関数の符号の約束 (`dz` = 源 − 試験粒子)、`fw` の規格化は K64 のソースから。 |

## 付録 B. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `pbunch` | $N$ (バンチあたりの粒子数) | `tmacro` (`PBUNCH`) | `PBUNCH` |
| `rclassic` | $r_e$ | `tmacro` | `ELRADI` |
| `gammab(l)` | 要素 $l$ の入口の $\beta\gamma$ (局所の設計運動量) | `ffs_pointer` | |
| `p1`、`h1` | $\beta\gamma(1+\delta)$、$\gamma$ | `twspac`、`twspace` | `p1`、`_p2h(p1)` |
| `f` | $N r_eL/(p_1^2h_1\sqrt{2\pi\sigma_z^2})$ | `twspac` | `f` |
| `az` | $A(z) = f\,e^{-(z-z_c)^2/2\sigma_z^2}$ | `twspac` | `az` |
| `xx1`、`yy1`、`xy1` | $z$ を固定した横の分散 | `twspac` | |
| `c1`、`s1` | 主軸の回転 $\cos\theta$、$\sin\theta$ | `twspac` | `c1`、`s1` |
| `sigx`、`sigy` | 主軸の $\sigma_x$、$\sigma_y$ | `twspac` | `sg1`、`sg2` |
| `fx`、`fy`、`fu` | $\nabla U$ と $U$ | `twspfu` | |
| `dg` | $\Delta\delta$ ((17.6)) | `twspac` | |
| `selfcod` | 力の中心: 重心 (真) か閉軌道 (偽) | `ffs_flag` (`SELFCOD`) | `selfcod` |
| `wakel(2,*)`、`waket(2,*)` | 縦・横のウェイク関数の表 (z, w) | `ffs_wake` | `wl`、`wt` |
| `lwl`、`lwt` | 表の行数 − 1 (追跡で使う行数) | `txwake` | |
| `itab(np)`、`izs(np)` | $z$ の昇順の並び、切片の番号 | `txwake`、`txwdefslice` | `itab`、`izs` |
| `ws`、`xs`、`ys`、`zs` | 切片の粒子数、座標の和、重心 | `txwake` | |
| `wz`、`wx`、`wy` | 切片 $n$ が受ける和 (17.13) | `txwake` | |
| `fw`、`fwp` | (17.14) の係数 | `tturn.f:259`、`txwake` | `fw` |
| `zw0l`、`zw0t` | 表の最初の $z$ と 0 の小さい方 (自己項の判定) | `txwake` | |
| `dzwr`、`dzlim` | バンチの幅、再整列の閾値 ($10^{-5}$) | `wakez` | `DZLIM` |
