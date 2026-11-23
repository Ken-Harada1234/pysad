# 17. 集団効果: 空間電荷、ウェイク、Touschek・IBS、CSR

**この章で分かること**: SAD の集団効果の機能 (空間電荷 `WSPAC`・`SPAC`・`PSPAC`、ウェイク `LWAKE`・`TWAKE`、バンチ内散乱 `INTRA` と Touschek 寿命、CSR と定常分布の方程式) が、何をどういう式で計算しているか。ヘルプには機能の説明しかなく、式は無い。そこで、標準的な出発点 (Lorentz 力、Gauss の法則、Poisson 方程式、クーロン散乱の断面積、シンクロトロン放射のスペクトル、Fokker-Planck 方程式) から式を導き、K64 のソースと対応させ、数値で確かめる。

**前提**: [16_undulator_beambeam.md](16_undulator_beambeam.md) (2 次元ガウス分布の電場 (16.6)〜(16.8))、[14_emittance.md](14_emittance.md) (ビーム行列、MINCOUP)、[13_radiation.md](13_radiation.md) (放射の減衰・励起)。

**範囲 (先に明記)**: 導出の深さは次のとおり。どの効果も、導出したのは式の構造と係数で、K64 の数値と一致するかの確認は、K64 のソースの式との照合 (読み) と、独立な数値 (数値積分・モンテカルロ・追跡) との比較で行った。gs を走らせた比較は本章では行っていない。

| 効果 | 導出 | 確認 | 節 |
|---|---|---|---|
| 空間電荷 (WSPAC) | Lorentz 力から $1/(\beta^2\gamma^3)$、ガウス分布のキック、$\delta$ のキックと壁の基準、条件付き分散、線形のレンズ、チューンシフト | 数値積分・sympy・追跡のチューンシフト | §3〜§4 |
| ウェイク (LWAKE・TWAKE) | ウェイク関数から $\delta$・横のキック、自己項の 1/2、K64 の切片の方法、解析式 | 追跡と解析式、切片の粗視化の粒子数依存 | §5〜§6 |
| 空間電荷 (SPAC・PSPAC) | 静止系の Poisson 方程式、実験室系のキック、円筒の壁、矩形の箱の級数解、長いバンチの極限 | 数値 (Poisson 方程式、箱と自由空間) | §7 |
| Touschek 寿命 | クーロン散乱の断面積から、損失の立体角、運動量の平均、寿命の式、平坦ビームの極限 | モンテカルロ、pysad の式との一致、K64 の当てはめ関数との比較 | §8 |
| IBS (INTRA) | 小角散乱の拡散と摩擦、拡散行列の構造 | モンテカルロと K64 の積分の比較 | §9 |
| CSR (自由空間の定常ウェイク) | シンクロトロン放射のスペクトルからインピーダンスと時間領域のウェイク | 数値積分 (周波数領域と時間領域) | §10 |
| 定常分布の方程式 (Haissinski) | Fokker-Planck の定常解 | K64 の Newton 解・追跡との比較 | §11 |

**読む順**: §1 (結論) → §3 → §5 → §8 (Touschek) → §13 (確認)。空間電荷だけなら §3、§4、§7。

**ヘルプの状況**: ヘルプに式は無い。あるのは機能の説明だけ: `WSPAC` (「strong-weak モードの空間電荷。ビームの大きさは `EMIT` + `CODPLOT` で求める。すべての次元でガウス分布を仮定、バンチあたりの粒子数は `PBUNCH`」)、`LWAKE` (縦ウェイクつきの光学)、`INTRA` (バンチ内散乱)、`TouschekLifetime` の使い方、`MINCOUP` の式。式・単位は K64 のソースと標準の物理から起こした。

## 1. 結論

- **空間電荷**: 長さ $L$ の区間で、粒子 ($z$、$(x, y)$) がバンチの電磁場から受けるキックは

$$
\Delta\vec p_\perp/p_0 = A(z)\,\nabla U(x,y),\qquad A(z) = \frac{N r_e L\,e^{-(z-z_c)^2/2\sigma_z^2}}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z},\qquad \Delta\delta = \frac{z-z_c}{\sigma_z^2}\,A(z)\,\bigl(U_{\rm ref} - U(x,y)\bigr)
$$

 ($U$ は 2 次元ガウス分布のポテンシャル (16.6)。式 (17.5)、(17.9))。$1/(\beta^2\gamma^3)$ は、電場 ($\propto\gamma$ 倍) と磁場の引力の打ち消し ($1-\beta^2 = \gamma^{-2}$)、通過時間 $L/\beta c$、慣性 $\gamma$ から出る (§3)。K64 の `twspac` の係数は、ソースの読み (`tintrb.f:420` の `f=pbunch*rclassic*al/(p1**2*h1*sqrt(2π sigzsq))`、$p_1^2h_1 = \beta^2\gamma^3$) と (17.5) が一致する。数値の確認 ([coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (b)) は、pysad の `WSpace.kick` を、ポテンシャルの数値積分から作った式と比べたもの (K64 の Fortran を走らせた比較ではない): 横のキックで絶対差 $5\times10^{-12}$ 以下 (値 $\sim7\times10^{-6}$。内部変数 $p/(1+\delta)$ のため相対 $|\Delta\delta| \sim 7\times10^{-7}$ の差が残る)、$\delta$ のキックで $3\times10^{-16}$ 以下。線形のレンズ $k_x = 2f/(\sigma_x(\sigma_x+\sigma_y))$ は式と $10^{-10}$ 以内。**縦のキック $\Delta\delta$ の基準** $U_{\rm ref}$ は、壁の半径 $b$ を入れると標準の縦の空間電荷の式 $2\ln(b/\sigma)$ の項を再現する。K64 は $U_{\rm ref}=0$ (§4.2)。
- **チューンシフト**: 各要素の入口のレンズの和 $\Delta Q = -\frac1{4\pi}\sum\beta k$ (§4.4)。pysad の追跡のチューンシフトは、0.5 GeV・周長 35.3 m のリングで、この和と $\Delta\nu_x$ で 0.9%、$\Delta\nu_y$ で 0.25% で一致した ([coll_03_laslett.out.txt](outputs/coll_03_laslett.out.txt))。
- **ウェイク**: 縦ウェイク $W_L(s)$ [V/C] ($s$ = 試験粒子の前方の源までの距離) による $\delta$ の変化は $\Delta\delta(z) = -\dfrac{eN}{mc^2p_0}\displaystyle\int\lambda(z')\,W_L(z'-z)\,dz'$ (自己項は $W(0)$ の 1/2)。横ウェイク $W_T(s)$ [V/C/m] は源の位置 $x(z')$ に比例: $\Delta p_x = +\dfrac{eN}{mc^2p_0}\displaystyle\int\lambda(z')\,x(z')\,W_T(z'-z)\,dz'$。ガウス分布のバンチ + 指数関数のウェイク $W_0e^{-s/\ell}$ の平均の損失は $k = \tfrac12e^{\sigma^2/\ell^2}\mathrm{erfc}(\sigma/\ell)$ に比例し、`txwake` の直訳 (pysad の `WakeKick`) は $\ell = 2, 5, 20$ mm で解析式の $1.027$、$1.0115$、$1.0073$ 倍。この差は切片の粗視化で、$1/(2n_sk)$ ($n_s \simeq \sqrt{n_p}$) と粒子数 2500〜160000 で合う (§6.2)。横ウェイクは $1.0115$ 倍。
- **SPAC・PSPAC** (§7): 静止系の Poisson 方程式の解を、実験室系の力 $\vec F_\perp = e\vec E'_\perp/\gamma$、$F_z = eE'_z$ に変換する。一様な円柱で $x'' = Kx/a^2$、$K = 2r_e\lambda/(\beta^2\gamma^3)$。長いバンチの 3 次元ガウス分布の場は、$\sigma_\zeta/\sigma_x = 3, 10, 30, 100$ で 2 次元の場 $\times\lambda(0)$ の $0.89, 0.979, 0.9965, 0.9996$ 倍に近づく (WSPAC・PSPAC の近似の条件)。接地した矩形の箱の級数解は、$a/\sigma \ge 20$ で自由空間の解と $r = 4\sigma$ でも 0.8% 以内。
- **Touschek 寿命** (§8): 静止系のクーロン散乱 (非相対論的な Møller の断面積) と、散乱後の縦運動量が運動量アクセプタンスを超える条件から、$1/\tau = r_e^2cN\,J/(\gamma\,V^*\sqrt{\lambda_1\lambda_2})$ ($J$ は 1 次元積分)。モンテカルロと 1.2σ 以内、pysad の損失率の式 (`touschek_piwinski`、v2 の既定) と係数まで一致 (比 1.000000)。横の運動量が 1 方向だけの極限は K64 の `touckf` (当てはめ関数) の関数形と $\varepsilon = 10^{-3}$〜$10$ で 1% 以内で一致する。
- **IBS** (§9): 小角クーロン散乱の拡散と摩擦から、静止系の運動量の共分散の増加率は $\propto\langle(q^2\delta_{ij}-3q_iq_j)/|q|^3\rangle$。K64 の `tintrb` の係数 $e_i = f_j+f_k-2f_i$ は、モンテカルロの平均と 4 組の固有値で $\pm0.3\%$ (統計誤差の範囲) で一致する。
- **CSR** (§10): 単一電子のシンクロトロン放射の低周波のスペクトルと因果律から、自由空間の定常ウェイクは $dE/ds = \dfrac{2Ne^2}{4\pi\epsilon_0\,3^{1/3}R^{2/3}}\displaystyle\int_0^\infty\lambda'(z+s)\,s^{-1/3}ds$ (符号は §10 の約束)。ガウス分布の平均の損失は $-\dfrac{\Gamma(2/3)}{3^{1/3}\Gamma(1/3)}\dfrac{Ne^2}{4\pi\epsilon_0R^{2/3}\sigma_z^{4/3}}$ ($=-0.3505\ldots$)。周波数領域と時間領域の数値積分が 8 桁で一致する。
- **定常分布の方程式** (§11): $\rho(z) = K\exp[-z^2/(2\sigma_{z0}^2) - a_wV(z)]$。K64 の `csrhaissin` の解は、独立に書いた Newton 解と $\sigma_z$ で $0.5\%$ 以内 (強い結合で $10^{-4}$)、K64 の追跡の平衡と、ウェイクの強さ $a_{wH} \le 1.25$ で $\sigma_z$ が $0.2\%$ 以内で一致する。$a_{wH} = 2.5$ では追跡の $\sigma_z$ が時間変動し (15%)、$\sigma_\delta$ も 24% 増えて、$\sigma_z$ は定常解より 3% 大きい。定常解が不安定で、追跡は揺らぐ状態になっている可能性がある (安定性は未解析。§11.3)。
- **K64 の疑い** (§17): WSPAC の縦の力の電位の基準が軸上 (壁の半径に依る $\ln(b/\sigma)$ の項が無い)、WSPAC の Emittance[] の最後のスライスでレンズが 8 倍 (B127)、SPAC の力が $1/\gamma$ 倍 (B89)、IBS の係数の $\gamma^4$ (導出は $(\beta\gamma)^4$、§9.3)。

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
| $\lambda_i$ | 静止系の運動量 (単位 $p_0$) の共分散 $P$ の固有値 (IBS)。Touschek では $(mc)^2$ 単位 | `eig` (`tintrb`) | `l1`、`l2` (`touschek.py`) |
| `ci`、`cintrb` | IBS の拡散行列の係数 | `tintrb`、`temit.f:801` | |
| `aw`、`inten` | CSR のキックの強さ (`inten`) と、定常分布の方程式の係数 `aw = inten/r65/sigz**2` | `csrtrack.f` | `aw` |

## 3. 空間電荷の力 (WSPAC の物理)

### 3.1 Lorentz 力

バンチの静止系 (記号に $'$) で、粒子の線密度 (単位長さあたりの個数) は $\lambda' = \lambda/\gamma$ ($\lambda$: 実験室系の線密度。長さが $\gamma$ 倍に伸びるので密度は $1/\gamma$ 倍)。長いバンチ ($\sigma_\zeta \gg \sigma_\perp$。条件は §7.3) では $\zeta$ 方向の電場 $E'_\parallel$ は横の電場より小さいので無視し、横方向は 2 次元の電場とする。丸い分布 ($\sigma_x = \sigma_y = \sigma$) の横の密度は $\dfrac{\lambda'}{2\pi\sigma^2}e^{-r^2/2\sigma^2}$。

**Gauss の法則** (半径 $r$、単位長さの円筒): 円筒の中の電荷は、$\int_0^r\dfrac{\lambda'}{2\pi\sigma^2}e^{-s^2/2\sigma^2}\,2\pi s\,ds\cdot e = e\lambda'\left(1 - e^{-r^2/2\sigma^2}\right)$ ($s^2/2\sigma^2 = t$ と置くと $\int_0^{r^2/2\sigma^2}e^{-t}dt$)。円筒の側面の電場の面積分は $2\pi r\,E'_r$ なので、

$$
E'_r = \frac{e\lambda'}{2\pi\epsilon_0r}\left(1 - e^{-r^2/2\sigma^2}\right),\qquad \frac{eE'_r}{mc^2} = \frac{2r_e\,\lambda'}{r}\left(1 - e^{-r^2/2\sigma^2}\right)
\tag{17.1}
$$

(2 つ目の等号は $e^2/(4\pi\epsilon_0) = r_emc^2$、すなわち $e^2/(2\pi\epsilon_0) = 2r_emc^2$ による。)

**実験室系に戻す**: 横方向の電場は $E_\perp = \gamma E'_\perp$ (横成分は $\gamma$ 倍)、磁場は、バンチが速度 $\beta c$ で動く電流の作る $B_\perp = \beta\gamma E'_\perp/c$ ($\vec B = \vec v\times\vec E/c^2$ の $\gamma$ 倍の変換)。同じ速度で動く試験粒子が受ける力 ($B$ の向きは、同じ向きに動く電荷が引き合う向き) は

$$
F_\perp = e\,(E_\perp - \beta c\,B_\perp) = e\,\gamma E'_\perp\,(1 - \beta^2) = \frac{eE'_\perp}{\gamma}
\tag{17.2}
$$

($1-\beta^2 = \gamma^{-2}$。**確認** ([coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (a)): sympy で $F\gamma/(eE') - 1 = 0$。電場の反発 ($\propto\gamma$) を磁場の引力 ($\propto\beta^2\gamma$) がほぼ打ち消す。) 長さ $L$ を通る時間は $\Delta t = L/\beta c$、試験粒子の運動量は $p_0 = \gamma m\beta c$。キックは

$$
\frac{\Delta p_\perp}{p_0} = \frac{F_\perp\Delta t}{p_0} = \frac{eE'}{\gamma}\,\frac{L}{\beta c\,\gamma m\beta c} = \frac{eE'L}{\gamma^2\beta^2\,mc^2}
\tag{17.3}
$$

(確認 (a): sympy で $(\Delta p/p_0)\,\gamma^2\beta^2mc^2/E' - 1 = 0$)。(17.1) の $eE'/mc^2$ と $\lambda' = \lambda/\gamma$ を入れて

$$
\frac{\Delta p_\perp}{p_0} = \frac{2r_e\lambda\,L}{\beta^2\gamma^3}\,\frac{\vec r}{r^2}\left(1 - e^{-r^2/2\sigma^2}\right)
\tag{17.4}
$$

**$1/(\beta^2\gamma^3)$ の由来**: $\gamma^{-2}$ (電場と磁場の打ち消し) $\times\ \gamma^{-1}$ ($\lambda' = \lambda/\gamma$、静止系の密度) $\times\ \beta^{-2}$ ($\Delta t\propto1/\beta$、$p_0\propto\beta$)。式 (17.3) の $1/(\gamma^2\beta^2)$ に、$\lambda' = \lambda/\gamma$ の $1/\gamma$ が掛かって $1/(\beta^2\gamma^3)$ になる。バンチ内の位置 $z$ での線密度は $\lambda(z) = N e^{-(z-z_c)^2/2\sigma_z^2}/(\sqrt{2\pi}\sigma_z)$ なので、$A(z)$ は

$$
A(z) = \frac{Nr_eL}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z}\,e^{-(z-z_c)^2/2\sigma_z^2}
\tag{17.5}
$$

(丸い場合の (17.4) = $A(z)\,\nabla U$、$U = \mathrm{Ein}(r^2/2\sigma^2)$、$\mathrm{Ein}(x) = \int_0^x(1-e^{-t})\,dt/t$。これは (16.6) の丸い場合: $\int_0^\infty dt\,[1 - e^{-r^2/(2\sigma^2+t)}]/(2\sigma^2+t)$ を $s = r^2/(2\sigma^2+t)$ で置換すると $\int_0^{r^2/2\sigma^2}(1-e^{-s})\,ds/s$ になる。$\nabla U = \dfrac{\vec r}{\sigma^2}\dfrac{1-e^{-x}}{x} = \dfrac{2\vec r}{r^2}(1-e^{-x})$、$x = r^2/2\sigma^2$ で (17.4) と一致する。) 楕円の分布は、(16.6)〜(16.8) の楕円ガウス分布の電位・場にそのまま置き換えればよい。ビーム・ビームでは電場と磁場が足し合って係数が $Nr_e/\gamma$ (対向)、ここでは打ち消し合って $1/(\beta^2\gamma^3)$ になる。

**この導出の条件**: (i) 長いバンチ ($\sigma_\zeta = \gamma\sigma_z \gg \sigma_\perp$。§7.3 で 3 次元の式との差を数値で示す)。(ii) 試験粒子の速度をバンチの速度と同じとする (試験粒子の $\delta$ による力の変化 $O(\delta/\gamma^2)$ を無視。K64 は `p1` に閉軌道の $\delta$ だけを使う)。(iii) 壁 (真空ダクト) は横のキックには効かない (中心のビームの場は壁の影響を受けない)。壁は縦のキックの基準に効く (§4.2)。

### 3.2 強-弱モデルであること

WSPAC は、ビームの 2 次モーメント (ビーム行列) を固定したガウス分布が、1 個の粒子に与える力だけを計算する (ヘルプの「strong-weak mode」)。粒子どうしの力、ビーム行列の自己無撞着な変化は入らない。ビーム行列は `Emittance[]` の結果 (放射平衡) か、pysad では `WSPACBEAM = {εx, εy, σz, σδ}` で与える。要素の入口ごとに力を入れる (区間の長さ分を一括して)。

## 4. WSPAC のキックと行列

### 4.1 $z$ を固定した横の分散

粒子が受けるキックは、粒子の位置 $z$ での横の分布で決まる。同じ $z$ にある粒子だけを集めた横の分布は、ビーム行列 $\Sigma$ (6×6) のガウス分布の**条件付き分布**。2 変数 $(x, z)$ の同時ガウス分布で、$x' = x - (\Sigma_{15}/\Sigma_{55})\,z$ を作ると、$\langle x'z\rangle = \Sigma_{15} - (\Sigma_{15}/\Sigma_{55})\Sigma_{55} = 0$。ガウス分布では無相関なら独立なので、$z$ を固定したときの $x$ の分散は $x'$ の分散に等しく

$$
\sigma_{xx}' = \Sigma_{11} - 2\frac{\Sigma_{15}}{\Sigma_{55}}\Sigma_{15} + \frac{\Sigma_{15}^2}{\Sigma_{55}^2}\Sigma_{55} = \Sigma_{11} - \frac{\Sigma_{15}^2}{\Sigma_{55}},\quad
\sigma_{yy}' = \Sigma_{33} - \frac{\Sigma_{35}^2}{\Sigma_{55}},\quad \sigma_{xy}' = \Sigma_{13} - \frac{\Sigma_{15}\Sigma_{35}}{\Sigma_{55}}
\tag{17.6}
$$

(`xx1`、`yy1`、`xy1`。$\sigma_{xy}'$ も同じ手順: $\langle x'y'\rangle$ を展開する。分散関数・クロマティックな項による $x$–$z$ の相関の影響を除く。$z$ に比例する平均のずれ $(\Sigma_{15}/\Sigma_{55})\,z$ は、K64 は入れない (中心は閉軌道か重心に固定)。平衡のバンチでは $\Sigma_{15} \approx 0$ なので影響は小さい。)

**主軸**: 2×2 の対称行列 $\begin{pmatrix}\sigma_{xx}'&\sigma_{xy}'\\\sigma_{xy}'&\sigma_{yy}'\end{pmatrix}$ を回転 $x_1 = x\cos\theta + y\sin\theta$、$y_1 = -x\sin\theta + y\cos\theta$ で対角にする。$\langle x_1y_1\rangle = \sigma_{xy}'\cos2\theta - \tfrac12(\sigma_{xx}'-\sigma_{yy}')\sin2\theta = 0$ から $\tan2\theta = 2\sigma_{xy}'/(\sigma_{xx}'-\sigma_{yy}')$ (K64: `theta = atan2(v,u)*.5`、$u = \sigma_{xx}'-\sigma_{yy}'$、$v = 2\sigma_{xy}'$)。このとき $\langle x_1^2\rangle = \tfrac12(\sigma_{xx}'+\sigma_{yy}') + \tfrac12(\sigma_{xx}'-\sigma_{yy}')\cos2\theta + \sigma_{xy}'\sin2\theta$ で、$(\cos2\theta,\sin2\theta) = (u,v)/a$、$a = \sqrt{u^2+v^2}$ を入れると

$$
\sigma_{x,y}^2 = \tfrac12(\sigma'_{xx}+\sigma'_{yy}) \pm\tfrac12\sqrt{(\sigma'_{xx}-\sigma'_{yy})^2 + 4\sigma_{xy}'^2}
\tag{17.7}
$$

($y_1$ は $-$ の側。固有値の公式と同じ。`sigx = (xx1+yy1)*.5`、`a = hypot(u,v)`。) 16 章 §7.2 の「$z$ を固定した条件付き包絡」と同じ操作。

### 4.2 追跡のキック (`twspac`) と縦のキックの基準

位置 $(x, y, z)$ の粒子について、主軸の座標 $(x_1, y_1) = R(\theta)(x - x_c, y - y_c)$ で、主軸のキックを (17.4) の楕円版で求める。実験室の座標に戻すには、回転を逆に ($\Delta p_x = \Delta p_{x_1}\cos\theta - \Delta p_{y_1}\sin\theta$、$\Delta p_y = \Delta p_{x_1}\sin\theta + \Delta p_{y_1}\cos\theta$。K64: `dpx = fx*c1 - fy*s1`、`dpy = fx*s1 + fy*c1`)。

$$
\Delta p_{x_1,y_1} = A(z)\,\nabla_{(x_1,y_1)}U(x_1,y_1;\sigma_x,\sigma_y)
\tag{17.8}
$$

(`fx = -az*fx` の符号は `twspfu` の戻り値の符号の約束。)

**$\delta$ のキックの導出**: (17.8) は力 $-\nabla V$ の形で、薄いキックの「相互作用エネルギー」$V(x_1,y_1,z) = -A(z)\,[U(x_1,y_1) - U_{\rm ref}]$ ($U_{\rm ref}$: ポテンシャルの基準 (定数)) から出る ($-\partial V/\partial x_1 = A\,\partial U/\partial x_1$)。$A$ は $z$ に依存するので、$V$ は $z$ にも依存する。SAD の $(z,\delta)$ は正準対 ($z$ が座標、$\delta$ が運動量) なので、薄いキックによる $\delta$ の変化は $\Delta\delta = -\partial V/\partial z = A'(z)\,[U - U_{\rm ref}]$。$A'(z) = -\frac{z-z_c}{\sigma_z^2}A(z)$ を入れて

$$
\Delta\delta = \frac{z - z_c}{\sigma_z^2}\,A(z)\,\bigl(U_{\rm ref} - U(x_1, y_1)\bigr)
\tag{17.9}
$$

(キック全体は、同じ $V$ から出るので、シンプレクティック。)

**基準 $U_{\rm ref}$ の物理**: $U_{\rm ref}$ は定数なので横の力には効かない。しかし $\Delta\delta$ には $x, y$ によらない項 $\frac{z-z_c}{\sigma_z^2}A\,U_{\rm ref}$ を足す。これは縦の電場の、横位置によらない部分で、**物理的には壁の半径で決まる**。丸いビームの電位を壁 ($r = b$、電位 0) から測ると、(17.1) を $r$ から $b$ まで積分して

$$
\phi'(r) = \int_r^b E'_rdr' = \frac{e\lambda'}{4\pi\epsilon_0}\bigl[\mathrm{Ein}(b^2/2\sigma^2) - \mathrm{Ein}(r^2/2\sigma^2)\bigr],\qquad
\mathrm{Ein}(b^2/2\sigma^2) \simeq \ln\frac{b^2}{2\sigma^2} + \gamma_E\ \ (b\gg\sigma)
\tag{17.10}
$$

($\gamma_E = 0.5772\ldots$ は Euler 定数。**確認** ([coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (b)): $\int_r^bE_rdr$ と $\tfrac12[\mathrm{Ein}(b^2/2\sigma^2)-\mathrm{Ein}(r^2/2\sigma^2)]$ が $r/\sigma = 0.5, 1.5, 3$、$b/\sigma = 20$ で 10 桁一致、$\mathrm{Ein}(b^2/2\sigma^2) - [\ln(b^2/2\sigma^2)+\gamma_E]$ は $10^{-15}$。) したがって $U_{\rm ref} = U(b)$ (壁の電位) と取るのが、物理的な縦の力に対応する。軸上 ($U = 0$) の縦の力は

$$
\Delta\delta\big|_{\rm axis} = -A'(z)\,U(b) = -\frac{r_eL}{\beta^2\gamma^3}\,\lambda'(z)\left[2\ln\frac b\sigma + \gamma_E - \ln2\right]
\tag{17.11}
$$

(ここで $A = r_eL\lambda(z)/(\beta^2\gamma^3)$ を使った。) これは標準の縦の空間電荷の式、$E_z = -\dfrac{g\,e}{4\pi\epsilon_0\gamma^2}\dfrac{d\lambda}{dz}$ ($g \simeq 2\ln(b/\sigma)$ + 定数) と $\Delta\delta = eE_zL/(\beta^2\gamma mc^2)$ で一致する。導出は §7.2 (静止系の縦の電場 $E'_z = -\partial\phi'/\partial\zeta$、$\partial/\partial\zeta = \gamma^{-1}\partial/\partial z$、$\lambda' = \lambda/\gamma$ から $\gamma^{-2}$)。**符号の物理**: $U_{\rm ref} = U(b)$ のとき、先頭 ($z > z_c$) の $\Delta\delta$ は正 (加速)、後尾は負 (減速) で、縦の空間電荷の向き (先頭が前へ押され、後尾が後ろへ押される) になる。K64 の $U_{\rm ref} = U(0,0) = 0$ では、軸外 ($U > 0$) で先頭が減速される向きになるが、これは基準の取り方の違い (定数 $U(b)$ の項を落とした) で、軸上の縦の力が 0 になる。K64 がこの定数を意図して落としたのか誤りかは、ソースだけでは決まらない (§17)。pysad (v2) は `WSPACRADIUS` ($> 0$) で $U_{\rm ref} = U(b)$ にする。**確認** (b): pysad の $\Delta\delta$ は、(17.9) を $U_{\rm ref} = 0$ として数値積分の $U$ から作った値と、**絶対差 $3\times10^{-16}$ 以下** (値 $\sim10^{-6}$) で一致する。

**運動量の更新** (`twspac` の末尾): SAD の内部変数 $p_x$ は $P_x/(p_0(1+\delta))$。キックの前の正準運動量は $p_xp_r$ ($p_r = 1+\delta_{\rm old}$)、キック後は $\delta_{\rm new} = \delta_{\rm old} + \Delta\delta$ なので

$$
p_x \leftarrow \frac{p_xp_r + \Delta p_x}{1+\delta_{\rm new}}
\tag{17.12}
$$

($\Delta p_x$ は $p_0$ で規格化した正準運動量の変化)。$dv$ も更新する (`dv = -dpr*(1+pr1)/h1/(h1+pr1*h0) + dvfs`、$h_1 = p2h(p_0\,pr_1)$)。

### 4.3 線形のレンズ (`twspace`)

小振幅 ($r \ll \sigma$) のキックは、(16.6) の $U$ を原点のまわりで展開して得る。被積分関数の $1 - \exp(-x^2/(2\sigma_x^2+t) - y^2/(2\sigma_y^2+t))$ は 2 次まで $x^2/(2\sigma_x^2+t) + y^2/(2\sigma_y^2+t)$ なので、

$$
\frac{\partial U}{\partial x} = 2x\int_0^\infty\frac{dt}{(2\sigma_x^2+t)^{3/2}(2\sigma_y^2+t)^{1/2}},\qquad
\int_0^\infty\frac{dt}{(A+t)^{3/2}(B+t)^{1/2}} = \frac{2}{\sqrt A(\sqrt A+\sqrt B)}
$$

($A = 2\sigma_x^2$、$B = 2\sigma_y^2$ とおくと、$\sqrt A\,(\sqrt A+\sqrt B) = 2\sigma_x(\sigma_x+\sigma_y)$ なので $\partial U/\partial x = 2x/(\sigma_x(\sigma_x+\sigma_y))$。積分公式は sympy で確認 ([coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (e))。丸い場合は $x/\sigma^2$ で、(17.4) の小振幅の極限 $2r_e\lambda L\cdot r/(2\sigma^2)$ と一致。) $y$ も同様 ($x\leftrightarrow y$)。中心 ($z = z_c$) のレンズは

$$
k_x = \frac{f'}{\sigma_x},\quad k_y = \frac{f'}{\sigma_y},\qquad f' = \frac{2Nr_eL}{\beta^2\gamma^3\sqrt{2\pi}\,\sigma_z\,(\sigma_x+\sigma_y)}
\tag{17.13}
$$

($\Delta p_x = k_xx$、$\Delta p_y = k_yy$。`akx = f/sx`、`aky = f/sy`、`f = 2 pbunch rclassic al/(p1**2 h1 (sx+sy) sqrt(2π sigzsq))`。) 主軸が傾いていれば、レンズの行列は $x_1,y_1$ で対角なので、回転して戻す: $k_{xx} = k_x\cos^2\theta + k_y\sin^2\theta$ など。`trans(2,1) = aks + akd c2` は $aks = (k_x+k_y)/2$、$akd = (k_x-k_y)/2$、$c_2 = \cos2\theta$ とした同じ式 ($\cos^2\theta = (1+c_2)/2$ から $k_{xx} = aks + akd\,c_2$)。**確認** (b): `WSpace.linear_k` は (17.13) と $10^{-10}$ 以内。

### 4.4 チューンシフト (標準の式)

薄いレンズ ($x'' $ に $+kx$ の発散。空間電荷は発散) をリングの 1 点に入れる。その点の Twiss パラメータで、1 周の行列は $M_0 = \begin{pmatrix}\cos\mu+\alpha\sin\mu & \beta\sin\mu\\-\gamma_T\sin\mu & \cos\mu-\alpha\sin\mu\end{pmatrix}$、レンズは $L = \begin{pmatrix}1&0\\k&1\end{pmatrix}$。積の跡は $\tfrac12\mathrm{tr}(M_0L) = \cos\mu + \tfrac12\beta k\sin\mu$ (sympy で確認、[coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (e))。レンズを入れたあとの位相 $\mu'$ は $\cos\mu' = \cos\mu + \tfrac12\beta k\sin\mu$。$k$ が小さいとき $\cos\mu' \simeq \cos\mu - \sin\mu\,\Delta\mu$ から $\Delta\mu = -\beta k/2$、

$$
\Delta Q = \frac{\Delta\mu}{2\pi} = -\frac{\beta k}{4\pi}
\tag{17.14}
$$

(要素の入口ごとのレンズの和 $\Delta Q_x = -\frac1{4\pi}\sum\beta_xk_x$ ($\beta_x$: Twiss)。$k$ が小さいので、位相の増分は足し合わせられる。) (17.13) の $k_x = f'/\sigma_x$ を入れ、レンズの長さ $L \to ds$ の連続分布の極限をとると

$$
\Delta Q_x = -\frac{1}{4\pi}\oint\beta_x\,\frac{2r_e\lambda_0}{\beta^2\gamma^3\,\sigma_x(\sigma_x+\sigma_y)}\,ds = -\frac{r_e\,\lambda_0}{2\pi\beta^2\gamma^3}\oint\frac{\beta_x\,ds}{\sigma_x(\sigma_x+\sigma_y)},\qquad \lambda_0 = \frac{N}{\sqrt{2\pi}\sigma_z}
\tag{17.15}
$$

(バンチ中心の粒子。$f'/L = 2r_e\lambda_0/(\beta^2\gamma^3(\sigma_x+\sigma_y))$ で、係数 $\frac{1}{4\pi}\times2 = \frac1{2\pi}$。標準の空間電荷チューンシフトの式。)

**確認** ([coll_03_laslett.out.txt](outputs/coll_03_laslett.out.txt)): 例 `example/05_tracking/flag_WSPAC_tune_shift.sad` の前半 (0.5 GeV の電子、8 セルのリング、周長 35.3 m、$N = 2\times10^{10}$、弱いスキュー四極 QS で結合を入れて $\varepsilon_y \neq 0$ にしてある。$\varepsilon_x = 1.98\times10^{-7}$、$\sigma_z = 18.4$ mm) を走らせた。WSPAC を入れた追跡の 1 周の行列から求めたチューンシフトは $\Delta\nu_x = -2.991\times10^{-5}$、$\Delta\nu_y = -1.203\times10^{-3}$。各要素の入口の $(17.13)$ のレンズ ($\Delta p = kx$ の $k$、(17.8) の主軸の傾き込み) と Twiss から作った (17.14) の和は $\Delta\nu_x = -3.018\times10^{-5}$、$\Delta\nu_y = -1.200\times10^{-3}$。比は $0.9911$ ($x$)、$1.0025$ ($y$)。$y$ は垂直のエミッタンスが小さく ($\sigma_y = 21\ \mu$m)、空間電荷が強い方向で、一致は 0.25%。$x$ は 0.9% で、$\Delta\nu_x$ の絶対値は $\Delta\nu_y$ の 1/40 しかなく、結合 (QS) による $x$–$y$ の混合や試験粒子の差分の誤差が相対的に効く可能性がある (原因は確認していない)。K64 は、一周の最後のスライスで `gammab(nlat+1)` が配列の外で 0 になり、$p_1$ が半分、レンズが 8 倍強い (B127。$p_1^2h_1$ が 1/8 になる)。

## 5. ウェイク

### 5.1 ウェイク関数の定義と 1 個の源のキック

電荷 $q$ の点源が $z_s$ にあるとき、$z$ ($< z_s$、後方) の試験粒子 (電荷 $e$) が受けるエネルギーの変化は、ウェイク関数 $W_L$ [V/C] で

$$
\Delta E(z) = -e\,q\,W_L(z_s - z)
\tag{17.16}
$$

($W_L > 0$ は損失。SAD の約束は `dz = zs(m) - zs(n)` = 源 ($m$) の位置 $-$ 試験粒子 ($n$) の位置。$z > 0$ が先頭なので、因果的なウェイクは $dz > 0$ だけで 0 でない。) 横方向の力は、源の横のずれ $x_s$ に比例する (双極子ウェイク)。横ウェイク $W_T$ [V/C/m] で $\Delta p_xc = +e\,q\,W_T(z_s - z)\,x_s$ (これも $W_T > 0$ が「後続を源のずれの向きに引く」符号の約束)。

### 5.2 バンチ全体のキック

$N$ 個の粒子のバンチ (密度 $\lambda(z)$、$\int\lambda = 1$) は、各源のキックの和。電荷 $q = e$ の源が $N\lambda(z')dz'$ 個あるので、試験粒子のエネルギーの変化は $\Delta E(z) = -e\cdot eN\int\lambda(z')\,W_L(z'-z)\,dz'$。運動量の変化は $\Delta E = \beta c\,\Delta p$ より $\Delta\delta = \Delta p/p_0 = \Delta E/(\beta c\,p_0\,mc) = \Delta E/(\beta\,p_0\,mc^2)$ なので、$\beta\to1$ とすると

$$
\Delta\delta(z) = -\frac{eN}{mc^2\,p_0}\int\lambda(z')\,W_L(z'-z)\,dz',\qquad
\Delta p_x(z) = +\frac{eN}{mc^2\,p_0}\int\lambda(z')\,x(z')\,W_T(z'-z)\,dz'
\tag{17.17}
$$

(2 つ目は、$\Delta p_xc = e\cdot eN\int\lambda\,x\,W_T$ を同様に $p_0\,mc$ で割ったもの。K64 の式は $\beta\to1$ の近似で、厳密には $\Delta\delta$ に $1/\beta$ が掛かる。$e$ は C、$W$ は V/C、$mc^2$ は eV で、$eW$ が eV になる。)

**自己項**: 粒子自身が作るウェイクを自分が受ける項 ($z' = z$) は、ウェイクが $s = 0$ で飛び ($s < 0$ で 0、$s = 0^+$ で $W(0^+)$)、**$W_L(0)$ の 1/2** で数える (beam loading の基本定理)。理由: 自分と十分近い、対称な分布の 2 粒子の間 ($s$ の分布が $s = 0$ について対称) で、$s > 0$ になる確率が 1/2 だから。粒子 1 個の幅を $\varepsilon$ の対称な分布として、$\int\!\!\int\rho\rho'\,W(z'-z) \to W(0^+)\,P(z'-z>0) = \tfrac12W(0^+)$ ($\varepsilon\to0$)。(17.18) の $\sigma\to0$ の極限 $k\to1/2$ も同じ。K64 は、同じ切片の内部の対 ($n = m$) で、表の最初の点が $z \ge 0$ のとき $\tfrac12W(z_1)\,ws(n)$ を使う。表が負の $z$ から始まると、$n = m$ も補間で $W(0)$ (全重み) が使われる (`txwake.f` の `if(n /= m .or. zw0l < 0.d0)` の読み。§5.3)。

**解析式** (確認 (c) の基準): $W_L(s) = W_0e^{-s/\ell}$ ($s \ge 0$) とガウス分布 ($\sigma_z = \sigma$) の平均の損失。2 粒子の距離 $s = z'-z$ は、独立な 2 つの $N(0,\sigma^2)$ の差なので、分散 $2\sigma^2$ のガウス分布 $e^{-s^2/4\sigma^2}/(2\sigma\sqrt\pi)$。$s > 0$ の寄与だけ (源が先行):

$$
\langle W\rangle = W_0\int_0^\infty e^{-s/\ell}\,\frac{e^{-s^2/4\sigma^2}}{2\sigma\sqrt\pi}\,ds = \frac{W_0}{2}\,e^{\sigma^2/\ell^2}\,\mathrm{erfc}\!\left(\frac\sigma\ell\right)\equiv W_0\,k
\tag{17.18}
$$

**積分の計算**: $s = 2\sigma u$ とおくと $\int_0^\infty e^{-s/\ell-s^2/4\sigma^2}ds = 2\sigma\int_0^\infty e^{-u^2-2bu}du$、$b = \sigma/\ell$。平方完成 $-u^2-2bu = -(u+b)^2+b^2$ から $\int_0^\infty e^{-u^2-2bu}du = e^{b^2}\int_b^\infty e^{-v^2}dv = e^{b^2}\frac{\sqrt\pi}2\mathrm{erfc}(b)$。したがって $\langle W\rangle = W_0\,\frac1{2\sigma\sqrt\pi}\,2\sigma\,e^{b^2}\frac{\sqrt\pi}2\mathrm{erfc}(b) = \frac{W_0}2e^{\sigma^2/\ell^2}\mathrm{erfc}(\sigma/\ell)$。検算: $\sigma/\ell\to0$ で $k\to\frac12$ (自己項の 1/2 と同じ。源が先行する確率が 1/2)、$\ell\to\infty$ でも $\frac12$ ($W = W_0$ の定数で、$s>0$ の確率が 1/2)。$\sigma = 5$ mm で $\ell = 2, 5, 20$ mm の $k$ は $0.10540$、$0.21379$、$0.38517$ ([coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (c))。したがって、(17.17) の平均 (バンチ内の全粒子の平均) は

$$
\langle\Delta\delta\rangle = -\frac{eNW_0\,k}{mc^2\,p_0},\qquad\langle\Delta p_x\rangle = +\frac{eNW_{0t}\,x_0\,k}{mc^2\,p_0}
\tag{17.19}
$$

(2 つ目は、横ウェイクで全体を $x_0$ ずらしたとき。$x(z') = x_0$ が積分の外に出る。)

### 5.3 K64 の方法 (`txwake`): 切片と補間

粒子を $z$ で整列し (`spsort`)、$z$ の隣り合う間隔が閾値 $s_b$ ($\mathtt{nb} > 1$ のとき $(z_{\max}-z_{\min})/(2(n_b-1))$、そうでなければ全幅) を超えたところでバンチの境界とする (`txwdefslice`)。各バンチ ($n_m$ 個) を、$m_1 = \lfloor\sqrt{n_m}\rfloor$ として 1 切片あたり $\lceil n_m/m_1\rceil$ 個ずつの粒子の**切片**に分ける (`mb(m)=(nm+m1-1)/m1`。$n_m = 20000$ なら $m_1 = 141$、1 切片 143 個、切片の数 $n_s \simeq \sqrt{n_m}$)。切片 $n$ の粒子数 `ws(n)`、座標の和 `xs(n)`、`ys(n)`、重心 `zs(n)`。切片の組 ($n$, $m$) について

$$
\mathtt{wz}(n) = \sum_mW_L(z_m - z_n)\,\mathtt{ws}(m),\qquad \mathtt{wx}(n) = \sum_mW_T(z_m - z_n)\,\mathtt{xs}(m)
\tag{17.20}
$$

($W$ は表の線形補間。$\mathtt{xs}(m)$ は切片の $x$ の**和**で、切片内の平均 $\times$ `ws`。(17.17) の $\int\lambda\,x\,W$ に当たる。) 各粒子 ($k$) のキックは、切片の重心の間で**線形補間**: $dwz = (1-w)\mathtt{wz}(n) + w\,\mathtt{wz}(n+1)$、$w = (z_k - z_n)/(z_{n+1} - z_n)$。

$$
\delta_k \leftarrow \delta_k - \mathtt{fwp}\cdot dwz,\qquad p_{x,k} \leftarrow \frac{p_a\,p_{x,k} + \mathtt{fwp}\cdot dwx}{p_b},\quad \mathtt{fwp} = \mathtt{fw}/p_0,\quad \mathtt{fw} = \frac{|{\rm charge}|\,e\,N\,n_{\rm bunch}}{mc^2\,n_p}\cdot\frac12
\tag{17.21}
$$

($p_a = 1+\delta_{\rm old}$、$p_b = 1+\delta_{\rm new}$。**`fw` の導出**: (17.17) の $\Delta\delta = -\frac{eN}{mc^2p_0}\int\lambda W$ で、$\int\lambda(z')W\,dz'$ をマクロ粒子の和 $\frac1{n_p}\sum_mW\,\mathtt{ws}(m)$ にすると、$\mathtt{fwp}\cdot\mathtt{wz} = \frac{eN}{mc^2\,p_0\,n_p}\sum_m W\,\mathtt{ws}(m)$。$n_p$ は `TrackParticles` に渡した粒子数 ($\mathtt{ws}$ の和 = マクロ粒子の数) で、$1/n_p$ が $\lambda$ の規格化。$n_{\rm bunch}$ は同じバンチの数 (複数バンチの周回)、$|{\rm charge}|$ は粒子の電荷数。因子 $\tfrac12$ は、要素の前と後で 1 回ずつ同じ大きさのキックを入れるため (合計で (17.17) の 1 回分)。K64 は `tturn.f:259` で `fw=(abs(charge)*e*pbunch*anbunch/amass)/np0*.5d0`。) **確認** (c): この方法を直訳した pysad の `WakeKick` (20000 粒子、$\sigma_z = 5$ mm、$N = 10^{10}$、1.7 GeV、$W_0 = 10^{12}$ V/C) の $\langle\Delta\delta\rangle$ は、(17.18) の解析式の $1.027$ ($\ell = 2$ mm)、$1.0115$ (5 mm)、$1.0073$ (20 mm) 倍。符号 (損失)・規格化 $eN/(mc^2p_0)$・自己項の 1/2 が合う。差の原因は §6.2。

**光学 (`WAKEOPT`)**: `InitialOrbits` で与えた複数の軌道が互いにキックを受けるとして、軌道ごとの Twiss を作る (`tffswake.f`)。TRPT + INS の前提。本章では扱わない ([docs/13_wake.md](../docs/13_wake.md))。

## 6. ウェイクの確認

### 6.1 横ウェイク

[coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) (c): 縦ウェイクの規格化・符号・自己項は上のとおり。横ウェイクは、バンチ全体を $x_0 = 1$ mm ずらして $W_T = 10^{13}e^{-s/5\,\mathrm{mm}}$ V/C/m を入れたとき、$\langle\Delta p_x\rangle = +2.038\times10^{-9}$ (解析式の $1.0115$ 倍、縦と同じ粗視化の差)、$\Delta\delta = 0$ (縦ウェイクなし)。

### 6.2 切片の粗視化による過大評価

切片の方法では、同じ切片の内部の対は、切片の重心の距離 0 で数えられる。表が負の $z$ から始まる (本確認の表。$-0.05$ m から) ので、$n = m$ の項も全重み $W(0) = W_0$ で入る。実際には、切片の中の 2 粒子の半分は先行、半分は後続で、平均は $W_0/2$ (自己項と同じ議論)。1 個の粒子にとって、同じ切片にいる源の割合は $m/n_p \simeq 1/n_s$ ($m$ は 1 切片の粒子数) なので、過大評価は $W_0(1 - \tfrac12)/n_s$。平均の損失 $W_0k$ との比で

$$
\frac{\langle\Delta\delta\rangle_{\rm K64}}{\langle\Delta\delta\rangle_{\rm exact}} - 1 \simeq \frac{1}{2\,n_s\,k},\qquad n_s\simeq\sqrt{n_p}
\tag{17.22}
$$

$n_p = 20000$ で $n_s = 141$: $\ell = 2, 5, 20$ mm で $3.4\%$、$1.7\%$、$0.9\%$ (観測の $2.7\%$、$1.15\%$、$0.73\%$ は、1 つの乱数の種の値で、有限個の標本の統計誤差 (下表の振れ、$\pm0.3$〜$0.4\%$) の範囲)。**粒子数を変えた確認** ([coll_07_wake_slices.out.txt](outputs/coll_07_wake_slices.out.txt)): $n_p = 2500$、10000、40000、160000 (3 つの種の平均) で、観測の過大評価は

| $\ell$ [mm] | $n_p = 2500$ | 10000 | 40000 | 160000 |
|---|---|---|---|---|
| 2 (観測 / (17.22)) | 0.101 / 0.095 | 0.046 / 0.047 | 0.026 / 0.024 | 0.0128 / 0.0119 |
| 5 | 0.049 / 0.047 | 0.022 / 0.023 | 0.013 / 0.012 | 0.0064 / 0.0058 |
| 20 | 0.027 / 0.026 | 0.012 / 0.013 | 0.0071 / 0.0065 | 0.0034 / 0.0032 |

と $1/\sqrt{n_p}$ で減り、(17.22) と 15% 以内で一致する (振れの範囲)。したがって差は切片の粗視化 (同じ切片の内部の対の重み) で説明される。切片の重心の間の線形補間の誤差は、$W$ が滑らか ($e^{-s/\ell}$ で切片の幅 $\sim0.1$ mm $\ll\ell$) なので、ここでは無視できる。**注意**: 表が $z = 0$ から始まるとき、$n = m$ の項は $\tfrac12W(z_1)\,ws$ になり、この過大評価は無い (K64 のソース `txwake.f` の `else` 枝の読み。数値では確認していない)。

## 7. 空間電荷 SPAC・PSPAC

K64 の SPAC (3 次元) と PSPAC (横 2 次元 × 縦の線密度) は、WSPAC と違って、粒子の分布そのものから場を作る (マクロ粒子の分布から Poisson 方程式を解く)。ここでは、場の定義とキックへの変換を導く。メッシュ上の離散化 (`spkick.f`、`spch.c` の差分と形状関数) は導出しない (§17)。

### 7.1 静止系の Poisson 方程式から実験室系のキックへ

バンチの静止系で、粒子の電荷分布 ($e\sum_j\delta^3(\vec r-\vec r_j)$) の静電ポテンシャル $\phi'$ は $\nabla'^2\phi' = -\rho'/\epsilon_0$。$e\phi'/mc^2 = r_e\sum_j1/|\vec r'-\vec r'_j|$ ($e^2/(4\pi\epsilon_0) = r_emc^2$)。実験室系の力は、(17.2) と同じ変換: 横の力は $e\vec E'_\perp/\gamma$、縦 ($z$) は電場の縦成分が不変なので $F_z = eE'_z$。バンチの座標 $\zeta = \gamma z$ (静止系の長さ)、$\vec E' = -\nabla'\phi'$ から、長さ $L$ のキックは

$$
\frac{\Delta p_{x}}{p_0} = -\frac{r_c\,w\,L}{\beta^2\gamma^2}\frac{\partial\varphi}{\partial x},\qquad
\Delta\gamma = -r_c\,w\,L\,\frac{\partial\varphi}{\partial\zeta},\qquad \varphi \equiv \sum_{j}\frac{1}{|\vec r'-\vec r'_j|}
\tag{17.23}
$$

($r_c$ = 粒子の古典半径、$w$ = マクロ粒子 1 個の実粒子数 ($\mathtt{PBUNCH}/n_{p0}$)。横は (17.3) の $eE'L/(\beta^2\gamma^2mc^2)$ に $E' = -\partial\phi'/\partial x$ を入れたもの。縦は $\Delta\gamma = eE_zL/(mc^2)$。$\Delta p/p_0 = \Delta\gamma/(\beta^2\gamma)$ なので縦の運動量の変化は $\Delta\delta = \Delta\gamma/(\beta^2\gamma)$。) $\varphi$ は静止系の座標の関数で、$\varphi\propto\lambda' = \lambda/\gamma$ を含むので、(17.4) の $1/(\beta^2\gamma^3)$ に戻る。pysad の `Solver3D` はこの式の係数 $C = r_cwL/(\beta_0^2\gamma_0^2)$ を使う (docstring)。**K64 の SPAC は電場が $\gamma$ 倍分、$1/\gamma$ 倍に弱い (B89)**: 例えば 5 MeV・$10^{11}$ 個・半径 1 mm の一様な円・長さ 0.1 m のバンチの 1 m のドリフトで、gs の rms $x$ 0.6226 mm に対し、包絡方程式 $a'' = K/a$ の式は 1.639 mm、pysad (v2) は 1.623 mm で、gs は $K/\gamma$ の式 (0.645 mm) に近い ([docs/10_tracking.md](../docs/10_tracking.md) §9.7 (c) の記載。本章では取り直していない)。

壁 (円筒、半径 $R$、電位 0) は、Green 関数に像の電荷を足して入れる。$\zeta$ 方向にフーリエ変換 ($e^{ik\zeta}$) すると、自由空間の Green 関数 $1/|\vec r-\vec r'|$ は 2 次元の遮蔽された Green 関数 $2K_0(k\,d)$ ($d$ は横の距離。$\int d\zeta\,e^{-ik\zeta}/\sqrt{d^2+\zeta^2} = 2K_0(|k|d)$)。円筒座標の加法定理 $K_0(k|\vec r-\vec r'|) = \sum_m\varepsilon_mK_m(kr_>)I_m(kr_<)\cos m(\theta-\theta')$ ($\varepsilon_0 = 1$、$\varepsilon_{m\ge1} = 2$、数値で確認、[coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (f)) に、$r = R$ で消える像の項 $-2\sum_m\varepsilon_m\frac{K_m(kR)}{I_m(kR)}I_m(kr)I_m(kr')\cos m(\theta-\theta')$ を足す (pysad の `Solver3D` の docstring の式。$r=R$ で $K_m(kR)I_m(kr')$ が項ごとに相殺する)。

### 7.2 一様な円柱とペルビアンス、縦の電場

**一様な円柱** (半径 $a$、線密度 $\lambda$): 静止系で $E'_r = e\lambda'r/(2\pi\epsilon_0a^2)$ (Gauss の法則)。実験室系の力 $F = e\gamma E'(1-\beta^2)$、キックは長さあたり $\Delta p/p_0 = F/(\beta c\,\gamma m\beta c)$。$e^2/(4\pi\epsilon_0) = r_emc^2$ を入れると

$$
x'' = \frac{K\,x}{a^2},\qquad K = \frac{2r_e\,\lambda}{\beta^2\gamma^3}\quad(\text{ペルビアンス})
\tag{17.24}
$$

(sympy で $K - 2r_e\lambda/(\beta^2\gamma^3) = 0$、$\gamma = 1/\sqrt{1-\beta^2}$ を入れて確認: [coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (a)。ガウス分布の小振幅の (17.4) $\Delta p/p_0 = 2r_e\lambda L\,r/(2\sigma^2\beta^2\gamma^3)$ と、$a^2\to2\sigma^2$ で同じ形。)

**縦の電場**: 丸いガウス分布の壁を基準にした電位 (17.10) を、$\zeta$ で微分する。静止系の縦の電場 $E'_z = -\partial\phi'/\partial\zeta$、$\phi'\propto\lambda'(\zeta)$ (長いバンチで $\lambda'$ の変化がゆるやかなとき)。実験室系の縦の電場は $E_z = E'_z$ (縦成分は変わらない)。$\partial/\partial\zeta = \gamma^{-1}\partial/\partial z$、$\lambda' = \lambda/\gamma$ なので $E_z = -\dfrac{e}{4\pi\epsilon_0\gamma^2}\dfrac{\partial\lambda}{\partial z}\bigl[\mathrm{Ein}(b^2/2\sigma^2)-\mathrm{Ein}(r^2/2\sigma^2)\bigr]$。軸上 ($r=0$) の $\Delta\delta = eE_zL/(\beta^2\gamma mc^2)$ が (17.11) と一致する ($e^2/4\pi\epsilon_0 = r_emc^2$)。$\gamma^{-2}$ なので、縦の空間電荷は高エネルギーでは小さい。

### 7.3 長いバンチの極限: 3 次元ガウス分布

3 次元ガウス分布 (標準偏差 $\sigma_x,\sigma_y,\sigma_\zeta$、全電荷 $Q$) の静電ポテンシャルは、$1/r = \frac{2}{\sqrt\pi}\int_0^\infty e^{-r^2\tau^2}d\tau$ とガウス分布の畳み込み ($\int e^{-(x-x')^2\tau^2}\mathcal N(x';\sigma^2)dx' = (1+2\sigma^2\tau^2)^{-1/2}\exp[-x^2\tau^2/(1+2\sigma^2\tau^2)]$) から

$$
\phi(x,y,\zeta) = \frac{Q}{4\pi\epsilon_0}\,\frac{2}{\sqrt\pi}\int_0^\infty d\tau\prod_{i=x,y,\zeta}\frac{\exp\!\bigl[-x_i^2\tau^2/(1+2\sigma_i^2\tau^2)\bigr]}{\sqrt{1+2\sigma_i^2\tau^2}}
\tag{17.25}
$$

**確認** ([coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (c)): (i) $\sigma\to0$ で $\phi\to1/r$ (点電荷: $r=3$ で $3\phi = 1.00000000$)。(ii) Poisson 方程式 $\nabla^2\phi = -4\pi\rho$ (単位 $Q/4\pi\epsilon_0 = 1$) を差分 ($h = 0.05\sigma$) で確認: 点 $(0.5, 0.3, 0.2)$ で $\nabla^2\phi = -0.60473$、$-4\pi\rho = -0.60501$ (差 0.05%、差分の誤差の範囲)。(iii) 長いバンチの極限: $\sigma_x = 1$、$\sigma_y = 0.6$、点 $(0.8, 0.4, 0)$ の $E_x$ が、(線密度 $\lambda(0) = 1/(\sqrt{2\pi}\sigma_\zeta)$) $\times$ (2 次元ガウス分布の場 (16.6)) に近づく割合は、$\sigma_\zeta/\sigma_x = 3, 10, 30, 100$ で $0.8901, 0.9794, 0.9965, 0.9996$。差 $1-$(比) は $0.110$、$0.0206$、$0.0035$、$0.00044$ と、$\sigma_\zeta$ が大きいほど急に小さくなる (べきの精密な評価はしていない)。WSPAC・PSPAC の長いバンチの近似の誤差は、静止系の長さ $\gamma\sigma_z$ と横の大きさの比で決まる。通常の電子のリング ($\gamma\sigma_z \gg\sigma_\perp$) では十分小さい。

### 7.4 PSPAC: 接地した矩形の箱の中の 2 次元 Poisson 方程式

PSPAC は、横 2 次元の Poisson 方程式を矩形の箱 (一辺 $a\times b$、境界で電位 0) で解き、縦の線密度 $\ell(z) = \lambda(z)/N$ ($\int\ell\,dz = 1$) を掛ける: すべての縦の位置の粒子の横の分布が同じ ($\rho(x,y,z) = \rho_\perp(x,y)\ell(z)$) という仮定。(17.4) の導出と同じ手順で ($\Lambda(z) = wN\ell(z)$ を実粒子の線密度として) $\Delta p_x/p_0 = -\dfrac{r_cwL}{\beta^2\gamma^3}\ell(z)\,\dfrac{\partial\Phi}{\partial x}$ ($\Phi = \sum_j(-2\ln d_j)$、$d_j$ は粒子 $j$ までの横の距離。$-\partial_x(-2\ln d) = 2(x-x_j)/d^2$ が単位線電荷の反発力)。

**箱の解**: 箱の固有関数 $\sin(m\pi x/a)\sin(n\pi y/b)$ で展開すると、$\nabla^2$ は $-\pi^2(m^2/a^2+n^2/b^2)$ を掛けることに当たる。したがって、電位 $\Phi = \sum\Phi_{mn}\sin\sin$、電荷分布 $\rho_{mn}$ の係数は $\Phi_{mn} = \rho_{mn}/[\epsilon_0\pi^2(m^2/a^2+n^2/b^2)]$ (単位電荷あたり $\epsilon_0 = 1$ なら $1/[\pi^2(m^2/a^2+n^2/b^2)]$)。中心 ($a/2$) の丸いガウス分布の係数は、各方向の射影から $g_m = \frac2a\sin\frac{m\pi}2\,e^{-m^2\pi^2\sigma^2/2a^2}$ (奇数の $m$ だけ)。**確認** ([coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) (d)): 正方形の箱 (一辺 $a$) で、中心から $r = \sigma, 2\sigma, 4\sigma$ の点の電場を、級数 (奇数の $m,n$ を 401 まで) と自由空間の解 $(1-e^{-r^2/2\sigma^2})/(2\pi\epsilon_0 r)$ で比べる:

| $a/\sigma$ | $r = \sigma$ | $r = 2\sigma$ | $r = 4\sigma$ |
|---|---|---|---|
| 20 | 1.00008 | 1.00055 | 1.00757 |
| 40 | 1.00000 | 1.00003 | 1.00047 |
| 80 | 1.00000 | 1.00000 | 1.00003 |

(箱 / 自由空間の比)。像の電荷による差は、正方形の 4 回対称のため $(r/a)^4$ に比例する (表の差は、$a$ を 2 倍にすると約 16 分の 1)。規格化と係数が合っている。K64・pysad (`Solver2D`) はこの連続の固有値の代わりに、5 点差分の固有値 $(2-2\cos(\pi p/n_x))/\Delta x^2$ ($\simeq(p\pi/a)^2[1 - (p\pi/n_x)^2/12]$) を使う (メッシュが粗いと差が出る: [docs/10_tracking.md](../docs/10_tracking.md) の PSPAC の警告)。K64 の PSPAC の誤り (B93〜B96) は本節の式とは別。

## 8. Touschek 効果

同じバンチの 2 粒子が散乱して、一方の縦運動量が運動量アクセプタンスを超え、粒子が失われる損失率を導く。実験室系ではなく、バンチの静止系 (記号に $*$) で考える。

### 8.1 静止系での散乱と損失の条件

静止系の運動量 ($mc$ 単位) は、横方向が $p^*_\perp = p_0x'$ ($p_0 = \beta\gamma$、$x'\sim10^{-4}$ で $p^*\sim0.6$。非相対論的 ($p^*\ll1$) は粗い近似になりうる: §17)。2 粒子の運動量 $\vec p_1,\vec p_2$ の相対運動量を $\vec q = \vec p_1 - \vec p_2$ とする。質量中心系では各粒子の運動量は $\pm\vec q/2$、大きさ $|q|/2$ で、弾性散乱は方向 $\hat n$ だけを変える。**縦運動量の変化**: 縦 ($z$) 成分は $\Delta p^*_z = \tfrac{|q|}{2}(n_z - \hat q_z)\simeq\tfrac{|q|}2n_z$ ($\hat q_z\simeq0$: 静止系の縦の運動量の広がりを横に比べて無視する)。実験室系では $\Delta p_s = \gamma\Delta p^*_z$ (静止系の運動量の縦成分は $\gamma$ 倍、エネルギーの変化は非相対論で小さい: $\Delta E^*\sim p\Delta p/m$ が $\gamma\beta\,\Delta p^*_z$ の項に比べて $p^*/mc$ 倍) なので、

$$
\Delta\delta = \frac{\Delta p_s}{p_0\,mc} = \frac{\gamma\,\Delta p_z^*}{\beta\gamma\,mc} = \frac{\Delta p_z^*}{\beta} = \frac{|q|\,n_z}{2\beta}
\tag{17.26}
$$

($mc$ 単位)。粒子が失われる条件 $\Delta\delta > a_+$ (正の運動量アクセプタンス) は $|q|\,n_z > c_+$、$c_\pm \equiv 2\beta a_\pm$ (負側 $\Delta\delta < -a_-$ も同様で $n_z < -c_-/|q|$)。$|q| < c_\pm$ の衝突は、どう散乱しても失われない。

### 8.2 クーロン散乱の断面積 (出発点)

質量中心系の各粒子の速度 $\beta^* = |q|/2$。距離の遠い散乱は Rutherford の断面積: 換算質量 $\mu = m/2$、相対速度 $v_{\rm rel} = 2\beta^*c$、ポテンシャルの強さ $k = r_emc^2$ から $\dfrac{d\sigma}{d\Omega} = \left(\dfrac{k}{2\mu v_{\rm rel}^2}\right)^2\dfrac1{\sin^4(\chi/2)} = \dfrac{r_e^2}{16\beta^{*4}\sin^4(\chi/2)} = \dfrac{r_e^2}{|q|^4\sin^4(\chi/2)}$ ($\chi$: 散乱角、$16\beta^{*4} = |q|^4$)。2 粒子は同種のフェルミ粒子 (電子) なので、交換項とスピン平均を入れる (**入力として使う標準の結果**: 非相対論的な同種の電子の散乱の断面積、Born 近似 $\alpha/\beta^*\ll1$ で位相の因子を 1 とする):

$$
\frac{d\sigma}{d\Omega} = \frac{r_e^2}{|q|^4}\left[\frac1{s^4}+\frac1{c^4}-\frac1{s^2c^2}\right] = \frac{4r_e^2}{|q|^4}\,K(\chi),\quad K(\chi) = \frac{4}{\sin^4\chi}-\frac{3}{\sin^2\chi},\quad s = \sin\tfrac\chi2,\ c = \cos\tfrac\chi2
\tag{17.27}
$$

(**整理**: $s^2c^2 = \sin^2\chi/4$、$s^4+c^4 = 1-2s^2c^2 = 1-\sin^2\chi/2$ から $\dfrac1{s^4}+\dfrac1{c^4} = \dfrac{s^4+c^4}{s^4c^4} = \dfrac{16(1-\sin^2\chi/2)}{\sin^4\chi} = \dfrac{16}{\sin^4\chi}-\dfrac8{\sin^2\chi}$。$\dfrac1{s^2c^2} = \dfrac4{\sin^2\chi}$ を引くと $\dfrac{16}{\sin^4\chi}-\dfrac{12}{\sin^2\chi} = 4K(\chi)$。数値で確認: [coll_02_touschek.out.txt](outputs/coll_02_touschek.out.txt) (a)、$\chi = 0.3, 1.1, 2.5$ で差 $4\times10^{-13}$ 以下。)

### 8.3 損失の立体角の積分 $W(\kappa)$

$\hat q$ は横の面 ($x$-$y$ 面) にあるとし ($q_z = 0$)、$\hat q = \hat x$ とする。散乱方向 $\hat n$ に、$z$ 軸を極軸とする極角 $\theta$・方位角 $\psi$ をとる: $n_z = \cos\theta$、$\hat n\cdot\hat q = \sin\theta\cos\psi = \cos\chi$。したがって $\sin^2\chi = 1-\sin^2\theta\cos^2\psi$。損失の条件 $n_z > c/|q| \equiv 1/\kappa$ ($\kappa = |q|/c$) の領域で $K(\chi)$ を積分する ($d\Omega = d(\cos\theta)\,d\psi$、$u = \cos\theta$)。$\psi$ の積分は、$a = \sin^2\theta = 1-u^2$ として

$$
\int_0^{2\pi}\frac{d\psi}{1-a\cos^2\psi} = \frac{2\pi}{\sqrt{1-a}} = \frac{2\pi}{u},\qquad
\int_0^{2\pi}\frac{d\psi}{(1-a\cos^2\psi)^2} = \frac{\pi(2-a)}{(1-a)^{3/2}} = \frac{\pi(1+u^2)}{u^3}
$$

(1 つ目は $\int_0^{2\pi}\frac{d\psi}{1-a\cos^2\psi}$ の標準の積分。2 つ目は、$I(a) = \int\frac{d\psi}{1-a\cos^2\psi}$ とすると $\int\frac{d\psi}{(1-a\cos^2\psi)^2} = I+aI'$ ($\frac{1}{(1-a\cos^2)^2} = \frac1{1-a\cos^2}+\frac{a\cos^2}{(1-a\cos^2)^2}$ と $I' = \int\frac{\cos^2\psi}{(1-a\cos^2\psi)^2}$ から)、$I+aI' = \frac{2\pi}{\sqrt{1-a}}+\frac{\pi a}{(1-a)^{3/2}}$。数値確認: [coll_02_touschek.out.txt](outputs/coll_02_touschek.out.txt) (b)、$u = 0.2, 0.5, 0.9$ で一致。) よって $K = 4(1-\ldots)^{-2}-3(1-\ldots)^{-1}$ の $\psi$ 積分は $4\cdot\frac{\pi(1+u^2)}{u^3}-3\cdot\frac{2\pi}{u} = \frac{4\pi}{u^3}-\frac{2\pi}{u}$ (sympy: $2\pi(2-u^2)/u^3$)。$u$ について $1/\kappa$ から 1 まで積分して (片側)

$$
W(\kappa) = \int_{1/\kappa}^1\left(\frac{4\pi}{u^3}-\frac{2\pi}{u}\right)du = 2\pi\left[\kappa^2-1-\ln\kappa\right]\quad(\kappa>1),\qquad W = 0\ (\kappa\le1)
\tag{17.28}
$$

($\int4\pi u^{-3}du = -2\pi u^{-2}$ から $2\pi(\kappa^2-1)$、$\int2\pi u^{-1}du = 2\pi\ln\kappa$。) **確認** (b): sympy の積分が $2\pi(\kappa^2-\ln\kappa-1)$ と一致 (差 0)。独立に、極軸を入射方向にとった $(\chi,\phi)$ の 2 重積分の数値 (被積分関数 $\sin\chi\,K(\chi)$、$\sin\chi\sin\phi>1/\kappa$ の領域) は、$\kappa = 1.5, 3, 10$ で $5.30636922$、$43.36269787$、$607.56777659$ で (17.28) と 8 桁一致する。損失の断面積は

$$
\sigma_{\rm loss}(\vec q) = \frac{4r_e^2}{|q|^4}\left[W\!\left(\frac{|q|}{c_+}\right)+W\!\left(\frac{|q|}{c_-}\right)\right]
\tag{17.29}
$$

(正側と負側。散乱後の 2 粒子は $\hat n$ と $-\hat n$ に出るので、それぞれの $n_z$ の符号で、損失する粒子を 1 個ずつ数える。) $\kappa^2$ の項は、小角散乱 ($1/\sin^4\chi$ の裾) を $n_z>1/\kappa$ で切ったことから来る。

### 8.4 損失率: 数え上げと運動量の平均

**数え上げ**: 静止系の密度 $n(\vec r,\vec p) = n_s(\vec r)f(\vec p)$、$\int n_s d^3r = N$。損失する粒子の単位体積・単位時間あたりの数は、順序つきの対について $\iint nn'\,v_{\rm rel}\,\sigma$ だが、(i) 対を重複して数えないので $\tfrac12$、(ii) 同種粒子の断面積 (17.27) は「$d\Omega$ に粒子が見つかる」断面積で、全立体角で積分すると衝突 1 回あたり 2 個の出射粒子 ($\hat n$ と $-\hat n$) をそれぞれ数える。(17.29) は損失する出射粒子を 1 個ずつ数えているので、上の $\tfrac12$ だけが残る: $R^* = \tfrac12\iint nn'v_{\rm rel}\sigma_{\rm loss}$。空間の積分は、ガウス分布 $n_s$ の $\int n_s^2d^3r = N^2/V^*$、$V^* = (4\pi)^{3/2}\sqrt{\det\Sigma^*}$ ($\Sigma^*$: 静止系の位置の共分散。$z^* = \gamma z$ なので $\Sigma^*_{zz} = \gamma^2\sigma_z^2$、$V^* = \gamma(4\pi)^{3/2}\sigma_x\sigma_y\sigma_z$)。粒子あたりの率は

$$
\frac1{\tau^*} = \frac{R^*}{N} = \frac{Nc}{2V^*}\Big\langle|q|\,\sigma_{\rm loss}(\vec q)\Big\rangle_{\vec q},\qquad v_{\rm rel} = c|q|\ (\text{非相対論})
\tag{17.30}
$$

**運動量の平均**: 2 粒子は同じ位置にいる。同時ガウス分布の位置を固定したときの運動量は、条件付き分布 (§4.1 と同じ): 共分散 $P = \Sigma_{pp}-\Sigma_{px}\Sigma_{xx}^{-1}\Sigma_{xp}$、平均 $\Sigma_{px}\Sigma_{xx}^{-1}\vec x$。同じ位置の 2 粒子は平均が同じなので、相対運動量 $\vec q = \vec p_1-\vec p_2$ では平均が相殺し、$\vec q\sim\mathcal N(0,2P)$。K64 の `pl` はこの $P$ ($p_0$ 単位) で、pysad は $p_0^2$ を掛けて $mc$ 単位にする (`P = pl[:2,:2]*p0**2`)。縦の広がりは無視し、$x$-$y$ の $2\times2$ のブロックの固有値を $\lambda_1\ge\lambda_2$、$\alpha_i = 1/(4\lambda_i)$ とする。$\vec q = k(\cos\varphi,\sin\varphi)$ の密度は $g = \dfrac{\exp[-k^2(\alpha_1\cos^2\varphi+\alpha_2\sin^2\varphi)]}{4\pi\sqrt{\lambda_1\lambda_2}}$ ($\det(2P) = 4\lambda_1\lambda_2$ から $2\pi\sqrt{4\lambda_1\lambda_2}$)。$\varphi$ の積分は、$\alpha_1\cos^2\varphi+\alpha_2\sin^2\varphi = \frac{\alpha_1+\alpha_2}{2}-\frac{\alpha_2-\alpha_1}{2}\cos2\varphi$ と $\int_0^{2\pi}e^{z\cos2\varphi}d\varphi = 2\pi I_0(z)$ から $2\pi e^{-k^2\alpha_1}I_{0e}(z)$、$z = k^2(\alpha_2-\alpha_1)/2$、$I_{0e}(z) = e^{-z}I_0(z)$ ($\frac{\alpha_1+\alpha_2}2 k^2 = \alpha_1k^2+z$)。(17.29) の $|q|\sigma_{\rm loss} = 4r_e^2k^{-3}[W+W]$ を入れ、$d^2q = k\,dk\,d\varphi$ で

$$
\Big\langle|q|\sigma_{\rm loss}\Big\rangle = \frac{4r_e^2}{4\pi\sqrt{\lambda_1\lambda_2}}\,2\pi\sum_\pm\int_{c_\pm}^\infty dk\,k^{-2}W(k/c_\pm)\,e^{-k^2\alpha_1}I_{0e}(z) = \frac{2r_e^2}{\sqrt{\lambda_1\lambda_2}}\sum_\pm J(c_\pm)
$$

$$
J(c) = \int_c^\infty dk\,k^{-2}\,W(k/c)\,e^{-k^2\alpha_1}\,I_{0e}\!\left(\frac{k^2(\alpha_2-\alpha_1)}2\right)
\tag{17.31}
$$

(積分の下限は $W = 0$ ($k<c$) による。) **確認** ([coll_02_touschek.out.txt](outputs/coll_02_touschek.out.txt) (c)): $\lambda_1 = 4\times10^{-6}$、$\lambda_2 = 10^{-6}$、$c = 1.2\times10^{-3}$ で、$\vec q$ を 400 万個サンプリングして $4W/|q|^3$ の平均をとると $2.87209\times10^9\pm7\times10^5$、(17.31) の $(2/\sqrt{\lambda_1\lambda_2})J$ は $2.87300\times10^9$ (差は 1.2σ)。

**実験室系の寿命**: (17.30) と (17.31) から、静止系の率は $1/\tau^* = r_e^2cN\sum_\pm J/(V^*\sqrt{\lambda_1\lambda_2})$。静止系の時間は実験室の時間の $1/\gamma$ (時間の遅れ) なので $1/\tau = (1/\gamma)(1/\tau^*)$:

$$
\frac1{\tau_{\rm inst}} = \frac{r_e^2\,c\,N}{\gamma\,V^*\sqrt{\lambda_1\lambda_2}}\sum_{\pm}J(2\beta a_\pm)
\tag{17.32}
$$

(1 点の瞬間の率。リング全体は、各区間 (長さ $\mathtt{al}$) の率を $\mathtt{al}/C$ で重みづけして足す: $1/\tau = \sum_{\rm slice}(\mathtt{al}/C)/\tau_{\rm inst}$。) $\lambda_1\lambda_2 = p_0^4\sigma_{x'}^2\sigma_{y'}^2$ (局所の角度の広がり) を入れると、$1/\tau\propto N/(\gamma^2\sigma_x\sigma_y\sigma_z\,p_0^2\sigma_{x'}\sigma_{y'})$ で、標準の式の $\gamma^{-2}$ の依存性。**確認** (d): 3 GeV の電子、$\sigma_x = 100\ \mu$m、$\sigma_y = 10\ \mu$m、$\sigma_z = 5$ mm、$\sigma_{x'}\ \sigma_{y'}$ = $10^{-4}$、$3\times10^{-5}$、$N = 10^{10}$、$a = 2\%$、長さ 1 m の 1 区間で、上の式を `scipy` の積分で組んだ値は $9.938313\times10^{-5}$ s$^{-1}$、pysad の `piwinski_slice_rate` は $9.938313\times10^{-5}$ s$^{-1}$ (比 $1.000000$)。(17.32) の係数 (数え上げの $\tfrac12$、$\tfrac1\gamma$、$V^*$、$I_{0e}$) まで、pysad の式と一致する。

### 8.5 平坦ビームの極限と K64 の `touckf`

横の運動量の広がりが 1 方向だけ ($\lambda_2\to0$、$z\to\infty$) の極限では $I_{0e}(z)\to1/\sqrt{2\pi z}\simeq2\sqrt{\lambda_2/\pi}/k$ なので、(17.32) は ($\sum_\pm$ の 1 つについて)

$$
\frac1{\tau_{\rm inst}} = \frac{2r_e^2cN}{\gamma V^*\sqrt{\pi\lambda_1}}\int_c^\infty dk\,k^{-3}W(k/c)\,e^{-k^2\alpha_1}
$$

$t = \alpha_1k^2$ (したがって $k^{-3}dk = \alpha_1dt/(2t^2)$、$\kappa = k/c = \sqrt{t/\varepsilon}$、$\varepsilon \equiv \alpha_1c^2 = \beta^2a^2/\lambda_1 = (a/\gamma\sigma')^2$) と置くと

$$
\frac1{\tau_{\rm inst}} = \frac{\sqrt\pi\,r_e^2\,c\,N}{2\gamma V^*\lambda_1^{3/2}}\,H(\varepsilon),\qquad
H(\varepsilon) = \int_\varepsilon^\infty\frac{dt}{t^2}\left[\frac t\varepsilon-1-\frac12\ln\frac t\varepsilon\right]e^{-t}
\tag{17.33}
$$

(1 方向の運動量の広がりの標準の式の形。$\varepsilon$ は、運動量アクセプタンスを静止系の運動量の広がりで割った量の 2 乗。) K64 の `tintrb`: `tf=al/vol/(ptrans*p0)**3` が $1/(V^*\lambda_1^{3/2})$ に当たり、`touckf(aez/eztrans)` が $H(\varepsilon)$ の当てはめ。**確認** (e): $H(\varepsilon)$ と `touckf`($=\ln(1/(x\varepsilon_E)+1)e^{-x}(b+\varepsilon_Ex)/(b+x(a+2x))/x$、$a = 5.62966$、$b = 0.75159$、$\varepsilon_E = 7.98221$) の比は、$\varepsilon = 10^{-3}, 10^{-2}, 0.03, 0.1, 0.3, 1, 3, 10$ で $1.0009$、$1.0018$、$0.9991$、$0.9957$、$1.0034$、$1.0070$、$0.9946$、$0.9904$ (1% 以内)。したがって K64 の `touckf` は (17.33) の $H(\varepsilon)$ の当てはめで、$\varepsilon\to0$ では $H\simeq\varepsilon^{-1}[\ln(1/\varepsilon)-\ln\varepsilon_E]$ の形になる (数値: $\varepsilon = 10^{-3}$ で $\varepsilon H(\varepsilon) = 4.853$、$\ln(1/\varepsilon)-\ln\varepsilon_E = 4.830$)。**この導出は、K64 の `touckf` の関数形が、静止系のクーロン散乱の損失断面積の平坦ビームの極限であることを示す**。K64 が `touckl`・`touckm` に積算した値を寿命に直す際の絶対値の係数 (`emit.n` の `TouschekLifetime`) は、本章では確かめていない (§17)。B162 (K64 の `touckf` の適用範囲): 平坦ビームの当てはめ (本節の極限) なので、2 方向の運動量の広がりがある場合 ((17.32) の $J$) は適用外。v2 は (17.32) の型が既定。

### 8.6 条件と限界

(i) 静止系で非相対論的 ($|q|\ll1$)。3 GeV で $p_0\sigma_{x'} = 0.6$ のように 1 に近い場合の相対論的な補正は入れていない。(ii) Born 近似 ($\alpha/\beta^*\ll1$)。しきい値付近 ($\beta^*\sim\beta a$) で $\alpha/\beta^*\sim1$ になりうる。(iii) 静止系の縦の運動量の広がりを無視 ($\sigma_\delta/\gamma\ll p_0\sigma_{x'}$)。(iv) 運動量アクセプタンス $a\gg\sigma_\delta$ (バンチ内の $\delta$ の分布の影響を無視。pysad の `touschek_piwinski` も同じ)。(v) 位置ごとの局所の値の $\mathtt{al}/C$ 平均 (局所運動量アクセプタンスは別: [docs/26](../docs/26_frequency_map_lma.md))。

## 9. バンチ内多重散乱 (IBS、`INTRA`)

IBS は、Touschek と同じクーロン散乱だが、損失ではなく、小角の散乱の累積による運動量の拡散と摩擦 (ビーム行列の成長) を扱う。

### 9.1 1 回の衝突の拡散と摩擦

静止系で、相対運動量 $\vec q$ の 2 粒子が、小さい角度 $\theta$ で散乱する。質量中心系で運動量が回転するので、$|q|$ は保存し、$\Delta\vec q = \Delta\vec q_\perp+\Delta q_\parallel\hat q$、$\Delta q_\parallel = -\frac12|\Delta\vec q_\perp|^2/|q|$ ($|q'|^2 = |q|^2$ から)。$|\Delta\vec q_\perp| = |q|\theta$ とし、$\Delta\vec q_\perp$ の方向は $\hat q$ に垂直な面で一様なので $\langle\Delta q_a\Delta q_b\rangle = w(\delta_{ab}-\hat q_a\hat q_b)$、$w = |q|^2\theta^2/2$、$\langle\Delta q_a\rangle = -(w/|q|)\hat q_a$。

粒子 1 の運動量は $\vec p_1 = \vec p_{\rm cm}+\vec q/2$、$\Delta\vec p_1 = \Delta\vec q/2$ ($\vec p_{\rm cm}$ は不変)。ビーム行列 $\langle p_ip_j\rangle$ の変化は、**拡散項** $\langle\Delta p_i\Delta p_j\rangle$ と**摩擦項** $\langle p_i\Delta p_j\rangle+\langle p_j\Delta p_i\rangle$ の和。$\vec p_{\rm cm}$ と $\vec q$ は独立で $\langle\vec p_{\rm cm}\rangle = 0$ なので、摩擦項は $\vec q/2$ だけが効いて $2\times\left(\tfrac{q_i}2\right)\left(-\tfrac{w}{2|q|}\hat q_j\right) = -\dfrac{w}{2}\dfrac{q_iq_j}{q^2}$。拡散項は $\tfrac w4(\delta_{ij}-\hat q_i\hat q_j)$。和は

$$
\frac{w}{4}\left(\delta_{ij}-\frac{q_iq_j}{q^2}\right)-\frac w2\frac{q_iq_j}{q^2} = \frac w4\left(\delta_{ij}-\frac{3q_iq_j}{q^2}\right)
\tag{17.34}
$$

(摩擦項は、2 つの粒子の運動量の差を縮める向き ($\Delta q_\parallel<0$) で、拡散項の $-2$ 倍の $q_iq_j$ を加える。)

### 9.2 小角の断面積と Coulomb 対数

小角の Rutherford (§8.2 の (17.27) 前の式): $d\sigma/d\Omega\simeq r_e^2/(|q|^4\sin^4(\theta/2))\simeq16r_e^2/(|q|^4\theta^4)$。$\int\theta^2d\sigma = 2\pi\int\theta^3\frac{16r_e^2}{|q|^4\theta^4}d\theta = \frac{32\pi r_e^2}{|q|^4}\ln\frac{\theta_{\max}}{\theta_{\min}}$ (対数発散を、最大・最小の角度で切る)。$\ln(\theta_{\max}/\theta_{\min})\equiv\ln\Lambda$。古典的な小角散乱は、衝突係数 $b$ で $\theta\simeq r_e/(\beta^{*2}b)$ ($\theta = 2k/(\mu v_{\rm rel}^2b)$) なので、$\ln\Lambda = \ln(b_{\max}/b_{\min})$。したがって $n$ を相手の密度として、$w$ の衝突率あたりの総和は $n\,v_{\rm rel}\int w\,d\sigma = nc|q|\cdot\frac{|q|^2}{2}\cdot\frac{32\pi r_e^2\ln\Lambda}{|q|^4} = \frac{16\pi r_e^2c\,n\ln\Lambda}{|q|}$ となり、(17.34) を掛けて

$$
\frac{d\langle p_ip_j\rangle}{dt^*} = 4\pi r_e^2\,c\,n\,\ln\Lambda\;\Big\langle\frac{q^2\delta_{ij}-3q_iq_j}{|q|^3}\Big\rangle_{\vec q}
\tag{17.35}
$$

($t^*$: 静止系の時間。$\langle\cdot\rangle$ は相対運動量 $\vec q\sim\mathcal N(0,2P)$ の平均。) 対角の和は $\sum_i(q^2-3q_i^2) = 0$ (トレースレス): 静止系の運動量空間の「エネルギー」$\sum\langle p_i^2\rangle$ はこの小角散乱で保存し、熱い方向 (大きな固有値) から冷たい方向へ移る。

**K64 の `bmin`・`bmax`** (`tintrb.f:223-231`): `bmax = min(bm, bn)`、`bm` = 最小のビームの大きさ (3 次元の主軸の最小)、`bn` = $(V/N)^{1/3}$ (粒子間の距離)。`bmin = max(rclassic/(ptrans*p0)**2, sqrt(vol/π/(ptrans*p0/h0*c)/pbunch/max(taurd)))`: 第 1 項は古典的な最小の衝突係数 $r_e/\beta^{*2}$、第 2 項は「減衰時間 (最大の $\tau_{\rm rd}$) の間に 1 回以上衝突する」断面積 $\pi b^2$ から決まる衝突係数 ($n\,v\,\tau\,\pi b^2 = 1$、$n = N/V$)。`log(2*bmax/bmin)` の定数 2 は K64 の取り方 (導出の $\ln(b_{\max}/b_{\min})$ との差は定数 $\ln2$ の取り方。未確認)。

### 9.3 K64 の拡散行列との対応

静止系の運動量の共分散 $P$ ($3\times3$、$p_0$ 単位。`pl`) の固有値を $\lambda_i$ (`eig`)、固有ベクトルを $R$ とする。(17.35) の平均を主軸で書くと、$\vec q$ の成分は $\mathcal N(0,2\lambda_i)$。K64 は

$$
f_i = 2\lambda_i\int_0^1\mathrm{fintrbu}(u;\lambda_i/\lambda_j,\lambda_i/\lambda_k)\,du,\qquad e_i = f_j+f_k-2f_i,\qquad \mathtt{pl}\leftarrow R\,\mathrm{diag}(e)\,R^\top
\tag{17.36}
$$

を作る (`fintrbu`: `u21 = 1-u²`、被積分関数 $4u\,u_{21}\sqrt{u_{21}^4(1-u_{21}^4)/[(1-u_{21}^4+a_1u_{21}^4)(1-u_{21}^4+a_2u_{21}^4)]}$)。$f_i$ は $\langle q_i^2/|q|^3\rangle$ に比例する量 ($\sum f_i\propto\langle1/|q|\rangle$) で、$e_i = \sum_jf_j-3f_i\propto\langle(q^2-3q_i^2)/|q|^3\rangle$ が (17.35) の主軸の対角成分に当たる。**確認** ([coll_06_ibs.out.txt](outputs/coll_06_ibs.out.txt)): 固有値 $(\lambda_1,\lambda_2,\lambda_3)$ = (1, 0.3, 0.1)、(1, 1, 0.2)、(4, 1, 0.05)、(1, 0.01, 0.001) の 4 組で、$\tau_{ii}\equiv(4\pi)^{3/2}\sqrt{\lambda_1\lambda_2\lambda_3}\,\langle(q^2-3q_i^2)/|q|^3\rangle$ (600 万サンプル) と $4\pi e_i$ が一致するかを調べた。比 $e_i/(\tau_{ii}/4\pi)$ は 12 個のうち 11 個が 0.996〜1.003 (統計誤差の範囲)、1 個が 1.027 ($e_2 = 2.15\times10^{-2}$ で相対的に小さい成分、$2\sigma$ の統計誤差)。$\sum e_i$ は $10^{-17}$ 以下 (トレースレス)。符号は、最大の固有値の方向が負 (冷える)、最小が正 (熱くなる) で、静止系で熱い横方向から冷たい縦方向へ運動量が移る (縦方向のエミッタンスの成長)。

**係数**: ビーム行列の増加 (静止系の運動量が $p_0$ 単位、長さ $\mathtt{al}$ の間、時間 $dt^* = \mathtt{al}/(\gamma\beta c) = \mathtt{al}/(p_0c)$) は、(17.35) に $n\to N/V^*$、運動量を $p_0$ 単位にして ($T\to T/p_0$、共分散 $/p_0^2$)

$$
\Delta\Sigma_{p_ip_j} = \frac{4\pi r_e^2N\ln\Lambda\,\mathtt{al}}{V^*p_0^4}\,T_{ij},\qquad T_{ij} = \Big\langle\frac{q^2\delta_{ij}-3q_iq_j}{|q|^3}\Big\rangle = \frac{\tau_{ij}}{(4\pi)^{3/2}\sqrt{\lambda_1\lambda_2\lambda_3}}
$$

$V^* = (4\pi)^{3/2}\sqrt{\det\Sigma^*_{xx}}$ ($=(4\pi)^{3/2}\mathtt{vol1}$) を入れて $\dfrac{r_e^2N\ln\Lambda\,\mathtt{al}}{(4\pi)^2\,\mathtt{vol1}\,\mathtt{pvol}\,p_0^4}\tau_{ij}$ ($\mathtt{pvol} = \sqrt{\lambda_1\lambda_2\lambda_3}$)。K64 は `ci = cintrb*al*log(2*bmax/bmin)/vol1/pvol/h0**4` (`cintrb = r_e²/(4π) * pbunch`、`temit.f:801`) に `pl` ($=\tau/4\pi$) を掛ける: $\dfrac{r_e^2N}{4\pi}\dfrac{\mathtt{al}\ln}{\mathtt{vol1\,pvol}\,h_0^4}\cdot\dfrac{\tau}{4\pi}$。**導出の式と一致し、違いは $p_0^4 = (\beta\gamma)^4$ が K64 では $h_0^4 = \gamma^4$ になっていること**。電子 ($\beta\approx1-1/(2\gamma^2)$、$\beta^4\approx1-2/\gamma^2$) では 3 GeV で $10^{-7}$ の差で、ガイドでは確認していない (K64 が $\beta\to1$ の近似にしているのか、導出の時間の換算の取り方が違うのか)。低エネルギー (陽子など $\beta<1$) では差が出る。この拡散の増加は、静止系 → 実験室系の座標変換 (`trans2` の逆) をして励起 $\Delta$ に足される (`beam(1:21) = bmi+beam(1:21)`、14 章の (14.2)) ので、IBS は $\Sigma = M\Sigma M^\top+\Delta(\Sigma)$ の反復 (`tintraconv`) になる。

### 9.4 IBS の限界

(i) 小角散乱のみ (大角の損失は Touschek: §8)。(ii) Coulomb 対数の取り方 ($b_{\max}$、$b_{\min}$、定数 2) は近似。(iii) 摩擦項は上の平均で導出 (運動量の分布がガウス分布で、3 次元の共分散で表せる場合)。(iv) 結合・分散のある場合の $P$ の取り方 (`beamsize` との関係) は導出していない (`tintrb` が行列で計算する部分のソースの読みだけ)。(v) 静止系の座標変換 (`trans1`、`trans2` のローレンツ変換) の詳細は導出していない。

## 10. CSR: 自由空間の定常ウェイク

ここでは、円軌道 (半径 $R$) の電子バンチが自由空間で出す、コヒーレントシンクロトロン放射 (CSR) の縦のウェイクを、標準のシンクロトロン放射のスペクトルから導く。真空ダクトの遮蔽・過渡的な効果は含まない (K64 の `Packages/CSR.n` はダクトを含むインピーダンスを数値で作る。§17)。

### 10.1 出発点: 単一電子のパワースペクトル

超相対論的な電子の 1 回の周回で、角周波数 $\omega$ あたりに出る平均のパワー ($\omega>0$ だけ) は、標準の結果として

$$
P_\omega = \frac{\sqrt3\,e^2\gamma}{8\pi^2\epsilon_0R}\,F\!\left(\frac\omega{\omega_c}\right),\quad F(y) = y\int_y^\infty K_{5/3}(t)\,dt,\quad \omega_c = \frac{3\gamma^3c}{2R}
\tag{17.37}
$$

(**規格化の確認** ([coll_04_csr_haissinski.out.txt](outputs/coll_04_csr_haissinski.out.txt) (a)(b)): $\int_0^\infty F\,dy = \frac12\int_0^\infty t^2K_{5/3}(t)\,dt = \frac{8\pi}{9\sqrt3} = 1.612266\ldots$ (数値積分で 12 桁一致)。全パワー $\int P_\omega d\omega = \frac{\sqrt3e^2\gamma}{8\pi^2\epsilon_0R}\omega_c\frac{8\pi}{9\sqrt3} = \frac{e^2c\gamma^4}{6\pi\epsilon_0R^2}$ ($\omega_c$ を入れた)。これは標準の全放射パワーで、比は 1.0 (数値)。) $\omega\ll\omega_c$ ($y\ll1$) の低周波では、小さい $t$ で $K_{5/3}(t)\simeq\frac12\Gamma(\tfrac53)2^{5/3}t^{-5/3}$ なので $\int_y^\infty K_{5/3}\simeq\frac12\Gamma(\tfrac53)2^{5/3}\cdot\frac32y^{-2/3}$、したがって

$$
F(y)\simeq\frac34\Gamma(\tfrac53)\,2^{5/3}\,y^{1/3} = \frac{4\pi}{\sqrt3\,\Gamma(1/3)}\left(\frac y2\right)^{1/3}
\tag{17.38}
$$

(2 つ目の等号は $\Gamma(5/3) = \tfrac23\Gamma(2/3)$、$\Gamma(1/3)\Gamma(2/3) = 2\pi/\sqrt3$ による。数値: 係数が両辺で $2.1495\ldots$ で一致、$y = 10^{-2}, 10^{-3}, 10^{-4}$ で $F(y)$ は漸近式の $0.961$、$0.992$、$0.998$ 倍。) (17.37) に代入すると、$\omega_c$ と $\gamma$ が消えて

$$
P_\omega\simeq\frac{e^2}{2\pi\epsilon_0\,\Gamma(1/3)\,3^{1/3}}\,\frac{k^{1/3}}{R^{2/3}},\qquad k = \omega/c
\tag{17.39}
$$

($y/2 = Rk/(3\gamma^3)$、$\sqrt3\gamma/(8\pi^2\epsilon_0R)\cdot4\pi/(\sqrt3\Gamma(1/3))\cdot(Rk/3\gamma^3)^{1/3} = \dfrac{\gamma}{2\pi\epsilon_0\Gamma(1/3)R}\dfrac{(Rk)^{1/3}}{3^{1/3}\gamma}$)。低周波のパワーは $\gamma$ によらない。

### 10.2 コヒーレントな放射とインピーダンス

剛体のバンチ ($N$ 個、規格化した線密度 $\lambda(z)$、$\int\lambda = 1$、$\tilde\lambda(k) = \int\lambda e^{ikz}dz$) を考える。粒子 $j$ の位置 $z_j$ からの放射の位相は $e^{ikz_j}$ で、全体の放射は $|\sum_je^{ikz_j}|^2$ に比例する。平均して $N+N(N-1)|\tilde\lambda|^2$。第 2 項がコヒーレントな放射で、パワーは

$$
P_{\rm coh} = N^2\int_0^\infty d\omega\,|\tilde\lambda(\omega/c)|^2P_\omega = N^2c\int_0^\infty dk\,|\tilde\lambda(k)|^2P_{ck}
\tag{17.40}
$$

**ウェイクの定義**: (17.16) と同じ約束で、粒子 1 個あたりの、長さあたりのエネルギー変化 (縦のウェイク関数 $w(s)$ [V/C/m]、源が試験粒子の前方 $s>0$): $dE/ds\,(z) = -Ne^2\int\lambda(z')\,w(z'-z)\,dz'$。全粒子で平均 (自己相関) すると、$\int dz\,\lambda(z)\lambda(z+s) = \int\frac{dk}{2\pi}|\tilde\lambda|^2e^{-iks}$ から

$$
\Big\langle\frac{dE}{ds}\Big\rangle = -Ne^2\int ds\,w(s)\int\frac{dk}{2\pi}|\tilde\lambda(k)|^2e^{-iks} = -Ne^2\int_0^\infty\frac{dk}\pi|\tilde\lambda(k)|^2\,\mathrm{Re}\,\hat w(k),\quad\hat w(k) = \int_0^\infty w(s)e^{iks}ds
\tag{17.41}
$$

(実数の $w$ では $\mathrm{Re}\,\hat w$ は $k$ について偶関数。) バンチ全体が失うエネルギーの時間あたりは $Nc|\langle dE/ds\rangle|$ で、これが放射のパワー (17.40) に等しい。任意の $\lambda$ について成り立つので

$$
\mathrm{Re}\,\hat w(k) = \frac{\pi}{e^2}P_{ck}\ \Longrightarrow\ \mathrm{Re}\,\hat w(k) = \frac{k^{1/3}}{2\epsilon_0\Gamma(1/3)3^{1/3}R^{2/3}}\quad(\text{(17.39) を入れた})
\tag{17.42}
$$

### 10.3 因果律と時間領域のウェイク

$w(s)$ は因果的 ($s<0$ で 0)。したがって $\hat w(k)$ は $k$ の上半面で解析的。自由空間・定常の問題では長さの尺度が $R$ だけなので、べき関数 $\hat w\propto k^{1/3}$ (17.42)。上半面で解析的なべき関数は $\hat w(k) = A(-ik)^{1/3}$ ($-ik$ の実部は上半面で正で、主値の枝)。実の $k>0$ で $(-i)^{1/3} = e^{-i\pi/6}$ なので $\mathrm{Re}\,\hat w = A\cos(\pi/6)k^{1/3} = \frac{\sqrt3}2Ak^{1/3}$。(17.42) と比べて

$$
A = \frac{1}{\sqrt3\,\epsilon_0\,\Gamma(1/3)\,3^{1/3}R^{2/3}}
\tag{17.43}
$$

(虚部は $\mathrm{Im}\,\hat w = -\frac12Ak^{1/3}$ と決まる。) 時間領域: $\int_0^\infty s^{\mu-1}e^{-(\varepsilon-ik)s}ds = \Gamma(\mu)(\varepsilon-ik)^{-\mu}$ (sympy で確認) より、$\mu = 2/3$ で $(-ik)^{-2/3}\leftrightarrow s^{-1/3}/\Gamma(2/3)$。$(-ik)^{1/3} = (-ik)\cdot(-ik)^{-2/3}$ で、$-ik$ を掛けることは $s$ での微分 ($\int w'e^{iks}ds = -ik\hat w$、$w(0)$ の境界項は超関数の意味で落とす) に当たるので

$$
w(s) = \frac{A}{\Gamma(2/3)}\frac{d}{ds}\left[s^{-1/3}\theta(s)\right] = -\frac{A}{3\Gamma(2/3)}\,s^{-4/3}\ (s>0)\ \text{(超関数)}
\tag{17.44}
$$

試験粒子が受けるエネルギー変化は、$\int\lambda(z')\,w(z'-z)\,dz'$ を部分積分して微分を $\lambda$ に移す ($s = z'-z$、$\frac{d}{ds}$ は $z'$ について):

$$
\frac{dE}{ds}(z) = +\frac{Ne^2A}{\Gamma(2/3)}\int_0^\infty\lambda'(z+s)\,s^{-1/3}ds = \frac{2Ne^2}{4\pi\epsilon_0\,3^{1/3}R^{2/3}}\int_0^\infty\lambda'(z+s)\,s^{-1/3}ds
\tag{17.45}
$$

(係数: $\dfrac{e^2A}{\Gamma(2/3)} = \dfrac{e^2}{\sqrt3\epsilon_0\Gamma(1/3)\Gamma(2/3)3^{1/3}R^{2/3}} = \dfrac{e^2}{2\pi\epsilon_0\,3^{1/3}R^{2/3}}$ ($\Gamma(1/3)\Gamma(2/3) = 2\pi/\sqrt3$。数値: $1/(\sqrt3\Gamma(1/3)\Gamma(2/3)) = 0.159154943092 = 1/2\pi$) $= \dfrac{2e^2}{4\pi\epsilon_0\,3^{1/3}R^{2/3}}$。符号: $z$ は先頭が正 (§5.1)。先頭 ($\lambda'<0$ の側) は損失、後尾は利得。) 単位: $e^2/(4\pi\epsilon_0)$ [J m] $\times R^{-2/3}\times\int\lambda'(z+s)s^{-1/3}ds$ [$\mathrm m^{-2}\cdot\mathrm m^{2/3} = \mathrm m^{-4/3}$] で J/m。

**ガウス分布の平均の損失**: $\lambda = e^{-z^2/2\sigma^2}/(\sqrt{2\pi}\sigma)$ で (17.41) を使う ($|\tilde\lambda|^2 = e^{-k^2\sigma^2}$、$\int_0^\infty k^{1/3}e^{-k^2\sigma^2}dk = \Gamma(2/3)/(2\sigma^{4/3})$):

$$
\Big\langle\frac{dE}{ds}\Big\rangle = -\frac{\Gamma(2/3)}{3^{1/3}\Gamma(1/3)}\,\frac{Ne^2}{4\pi\epsilon_0\,R^{2/3}\sigma_z^{4/3}},\qquad\frac{\Gamma(2/3)}{3^{1/3}\Gamma(1/3)} = 0.350472\ldots
\tag{17.46}
$$

**確認** ([coll_04_csr_haissinski.out.txt](outputs/coll_04_csr_haissinski.out.txt) (d)): 単位 $Ne^2/(4\pi\epsilon_0)$、$R = \sigma = 1$ で、周波数領域の式 (17.41) の数値積分 $-0.35047200$、時間領域の式 (17.45) を $\lambda$ で平均した数値積分 $-0.35047200$、閉じた式 (17.46) $-0.35047200$ が 8 桁で一致する。(17.45) の $dE/ds$ は、$z = +2\sigma$ (先頭側) で $-0.131$、$z=-2\sigma$ (後尾側) で $+0.180$ (同じ単位、$Ne^2/(4\pi\epsilon_0)$ を除く)。

### 10.4 K64 の CSR の追跡との関係と限界

`csrtrack` は 1 周を `ndiv` 個に分けて、(i) $z\mathrel{+}=(r_{56}/n_{\rm div})\delta$、(ii) 密度 $\rho(z)$ (`csrrho`、線形補間で格子に載せる) とウェイクの畳み込み $w(z) = \int W(z-z')\rho(z')\,dz'$ (`csrwake`、インピーダンス `zl` を FFT で畳み込む)、(iii) $\delta\mathrel{+}=(a_w/n_{\rm div})\,w(z)$ + RF (線形の $r_{65}z$、または $v_r\sin(\omega z+\phi_{\rm rf})-dp_r$)、(iv) 放射減衰と量子励起 $\delta\to a_d\delta+b_d\,\xi$ ($a_d = 1-2d/n$、$b_d = \sqrt{4d/n\,(1-d/n)}\,\sigma_E$、$\xi$ は標準正規乱数) を繰り返す。$b_d$ は、励起のない平衡で $\mathrm{Var}(\delta) = b_d^2/(1-a_d^2) = \sigma_E^2$ になる値 ($1-a_d^2 = \frac{4d}n(1-\frac dn)$)。K64 の `zl` (`Packages/CSR.n`) は、ダクトの形状を含む Maxwell 方程式の数値解のインピーダンスで、本節の自由空間の式 (17.45) の $R^{2/3}\sigma_z^{4/3}$ のスケーリング (`CSR.n` の `Enorm`、ソースの読み) と同じ形だが、係数と遮蔽の効果は本章では確かめていない。**限界**: (i) 自由空間・定常 (円軌道の長い偏向磁石の内部) のみ。偏向磁石の入口・出口の過渡的なウェイク、直線部のウェイク、ダクトによる遮蔽 (低周波のカットオフ) は導出していない。(ii) 因果律とべき則から複素インピーダンスを決めた (実部は放射のパワーから、虚部は解析性から)。時間領域の式 (17.45) を、リエナール・ヴィーヒェルトの場から直接導く方法 (遅延の条件) では確かめていない。(iii) `zl` の符号 ($\hat w$ と $e^{\pm ikz}$ の約束)・`csrwake` の畳み込みの向き (配列の添字の増える向きと $z$ の向き) は、`CSR.n` の側の約束で、確認していない。

## 11. 定常分布の方程式 (Haissinski 方程式)

### 11.1 Fokker-Planck の定常解

K64 の `csrtrack` の 1 周の写像は、$z\mathrel{+}=r_{56}\delta$ ($\mathtt{ndiv} = 1$ とする)、$\delta\to a_d[\delta+F(z)]+b_d\xi$、$F(z) = r_{65}z+a_w\,v(z)$ ($v(z) = \int W(z-z')\rho(z')dz'$ はウェイクの畳み込み)。$\delta$ の運動は、1 周の変化が小さい ($\mu_s\ll1$、$d\ll1$) 連続の極限で Langevin 方程式 $\dfrac{dz}{dn} = r_{56}\delta$、$\dfrac{d\delta}{dn} = F(z)-2d\,\delta+b_d\xi$ になり、位相空間の密度 $\psi(z,\delta)$ は Fokker-Planck 方程式 ($b_d^2\simeq4d\sigma_E^2$ より拡散係数 $D = b_d^2/2 = 2d\sigma_E^2$)

$$
\frac{\partial\psi}{\partial n} = -r_{56}\,\delta\,\frac{\partial\psi}{\partial z}-F(z)\frac{\partial\psi}{\partial\delta}+2d\,\frac\partial{\partial\delta}\Bigl(\delta\,\psi\Bigr)+2d\,\sigma_E^2\frac{\partial^2\psi}{\partial\delta^2}
\tag{17.47}
$$

に従う。$F = -U'(z)$ とする (ポテンシャル $U$)。**定常解**:

$$
\psi_\infty(z,\delta) \propto \exp\!\left[-\frac{\delta^2}{2\sigma_E^2}-\frac{U(z)}{r_{56}\,\sigma_E^2}\right]
\tag{17.48}
$$

(sympy で右辺が 0 になることを確認: [coll_04_csr_haissinski.out.txt](outputs/coll_04_csr_haissinski.out.txt) (e)。手で確かめるには、$\partial_z\psi = -\frac{U'}{r_{56}\sigma_E^2}\psi$、$\partial_\delta\psi = -\frac{\delta}{\sigma_E^2}\psi$ から、第 1 項 $-r_{56}\delta\,\partial_z\psi = +\frac{\delta U'}{\sigma_E^2}\psi$、第 2 項 $-F\partial_\delta\psi = -(-U')\left(-\frac{\delta}{\sigma_E^2}\right)\psi = -\frac{\delta U'}{\sigma_E^2}\psi$ で打ち消し合う。第 3 項 $2d\,\partial_\delta(\delta\psi) = 2d\left(1-\frac{\delta^2}{\sigma_E^2}\right)\psi$ と第 4 項 $2d\sigma_E^2\partial_\delta^2\psi = 2d\left(\frac{\delta^2}{\sigma_E^2}-1\right)\psi$ も打ち消し合う。) 密度 $\rho(z) = \int\psi d\delta$ は $\exp[-U/(r_{56}\sigma_E^2)]$ に比例。

### 11.2 Haissinski 方程式

$U = U_{\rm RF}+U_w$: 線形の RF の項 $F = r_{65}z$ は $U_{\rm RF} = -r_{65}z^2/2$、ウェイクの項は $U_w' = -a_w\,v(z)$。$-U_{\rm RF}/(r_{56}\sigma_E^2) = r_{65}z^2/(2r_{56}\sigma_E^2)\equiv-z^2/(2\sigma_{z0}^2)$、すなわち

$$
\sigma_{z0}^2 = \frac{r_{56}\,\sigma_E^2}{-r_{65}}
$$

(ウェイクなしの平衡のバンチ長。1 周の行列 $\begin{pmatrix}1&r_{56}\\r_{65}&1+r_{56}r_{65}\end{pmatrix}$ の固有構造から $\sigma_z^2/\sigma_\delta^2 = r_{56}/(-r_{65})$ とも一致する。)ウェイクの項は $-U_w/(r_{56}\sigma_E^2) = +\dfrac{a_w}{r_{56}\sigma_E^2}V(z)$、$V'(z) = v(z)$。$r_{65}\sigma_{z0}^2 = -r_{56}\sigma_E^2$ を使って $a_{wH}\equiv a_w/(r_{65}\sigma_{z0}^2) = -a_w/(r_{56}\sigma_E^2)$ とおくと

$$
\rho(z) = K\exp\!\left[-\frac{z^2}{2\sigma_{z0}^2}-a_{wH}\,V(z)\right],\qquad V'(z) = \int W(z-z')\,\rho(z')\,dz',\qquad\int\rho\,dz = 1
\tag{17.49}
$$

K64 の `tfcsrhaissin` は `aw = rbody(2)/r65/sigz**2` (`rbody(2)` = `inten`、`csrtrack` の `aw` と同じ量。`csrtrack` は `aw = inten/np`) で呼ぶので、$a_{wH} = \mathtt{inten}/(r_{65}\sigma_{z0}^2)$ と上の導出が一致する ($\mathtt{aw}_{\rm track}$ の $1/n_p$ は、`csrrho` が密度を個数 (規格化しない) で数えるため。`csrhaissin` の `rho` は $\sum\rho\,dz = 1$)。**定数の項**: ウェイクの直流成分 ($k=0$ の一定の力) は、同期位相のずれに吸収されるので除く (`zi(1) = zi(2) = 0`)。`csrwake` の畳み込み $v_n = \sum_mw[n-m]\rho_m$ は、格子の添字が増える向き ($z$ が増える向き) に因果的で、源は $z$ の小さい側。

**Newton 法** (`csrhaissin`): 未知数を $u_k = \ln\rho_k$ と $u_0 = \ln K$ とし、$f_k = u_k-u_0-u^{(0)}_k+a_{wH}V_k$ ($u^{(0)}_k = -z_k^2/(2\sigma_{z0}^2)$) $=0$ と規格化 $\sum\rho_k = 1$ を解く。ヤコビアン $g_{km} = \delta_{km}+a_{wH}\,\mathrm{wi}_{k-m}\rho_m$ ($\mathrm{wi}$: $V$ のウェイクの核。$V$ は $\rho_m = e^{u_m}$ について線形)、$g_{k0} = -1$、$g_{0m} = \rho_m$ の境界つきの線形方程式を `tsolvg`/`tsvd` で解き、ステップを半分にしながら $\sum|f_k|<10^{-5}$ まで反復する。

### 11.3 確認: K64 の追跡との比較

[coll_04_csr_haissinski.out.txt](outputs/coll_04_csr_haissinski.out.txt) (f)(g): $\mu_s = 0.04$ rad、$\sigma_E = 10^{-3}$、$\sigma_{z0} = 3$ mm ($r_{56} = 0.12$、$r_{65} = -0.0133$)、減衰率 $d = 4\times10^{-3}$/周、格子 $dz = 0.15$ mm ($n_r = 512$)、ウェイク $W(s) = e^{-s/\ell}$ ($\ell = 1$ mm、格子の向きは $z$ の小さい側が源) で、強さ `inten` を変える。`csrhaissin` の解の $\sigma_z$ と、独立に書いた Newton 型の求解 (実空間の同じ方程式)、`csrtrack` (K64 の直訳、20000 粒子、3000 周、最後の 1500 周の平均) の $\sigma_z$ を比べる。

| $a_{wH}$ | `csrhaissin` の $\sigma_z$ [mm] | 独立な Newton 解 [mm] | 追跡 $\sigma_z$ [mm] | 追跡の $\sigma_z$ の時間変動 | 追跡の $\sigma_\delta/\sigma_E$ | 追跡 / `csrhaissin` |
|---|---|---|---|---|---|---|
| 0.417 | 2.945 | 2.961 | 2.958 | 0.55% | 1.002 | 1.004 |
| 1.250 | 3.587 | 3.599 | 3.592 | 0.95% | 1.0015 | 1.001 |
| 2.500 | 4.704 | 4.705 | 4.847 | 14.9% | 1.244 | 1.030 |
| 3.750 | 5.585 | (求めていない) | 5.809 | 16.2% | 1.455 | 1.040 |

($a_{wH} = \mathtt{inten}/(r_{65}\sigma_{z0}^2)$。ウェイクなしは $\sigma_z = 3.000$ mm (追跡 3.007 mm)。表の最後の 1 行以外の追跡は $\mu_s = 0.04$・$d = 4\times10^{-3}$、$a_{wH} = 2.5$ では $\mu_s = 0.02$・$d = 2\times10^{-3}$ でも追跡 $4.841$ mm、変動 14.5%、比 1.029 で変わらない。)

**読み方**:
- `csrhaissin` の解 (K64 の Newton 法) は、独立に書いた Newton 型の求解と、$\sigma_z$ で $0.5\%$ ($a_{wH} = 0.42$)、$0.3\%$ (1.25)、$10^{-4}$ (2.5) の範囲で一致する。形 $\rho(z)$ の最大の差は最大値の 1%、3.6%、6.4% で、重心 $\langle z\rangle$ が $0.2$〜$0.9\%$ (格子の 0.03〜0.4 個分) 違うことによる (積分の離散化の違いと推測。未確認)。同じ方程式 (17.49) を解いていると読める。
- 追跡は、$a_{wH}\le1.25$ では $\sigma_z$ が `csrhaissin` と $0.4\%$ 以内で一致し ($\sigma_\delta$ はほぼ $\sigma_E$ のまま、$\sigma_z$ の時間変動は統計誤差の 1% 以内)、(17.49) が K64 の追跡の定常状態を表すことを示す。
- $a_{wH}\ge2.5$ では、追跡の $\sigma_z$ が時間とともに 15% 程度揺らぎ、$\sigma_\delta$ が $\sigma_E$ より 24% 以上大きくなり、平均の $\sigma_z$ は定常解より 3〜4% 大きい。$(\mu_s, d)$ を半分にしても変わらないので、有限の $\mu_s$・$d$ の効果ではない。定常解は存在するが不安定で、追跡は揺らぐ (乱流的な) 状態にある可能性がある。安定性の解析 (`CSRStability`) は本章では導出していない (§11.4)。

### 11.4 限界

(i) 連続極限 ($\mu_s\ll1$、$d\ll1$) の定常解で、1 周の写像の有限の $\mu_s$ の補正 ($O(\mu_s^2)$) は入れていない。(ii) 定常分布の存在 (方程式の解) のみで、その安定性 (CSR によるマイクロ波不安定性) は別の解析 (`CSRStability`、`csroy.f`) で、導出していない。(iii) ウェイクは $z$ だけの関数 (周回で記憶しない)。

## 12. 式とソースの対応表

| 効果 | 物理 (標準の式) | K64 | 状態 |
|---|---|---|---|
| **SPAC** (3 次元) | §7.1。静止系の 3 次元 Poisson 方程式 (円筒の導体壁) を解き、(17.23) で実験室系に戻す。一様な長い円柱でペルビアンス (17.24) | `spkick.f` (`spdrift_free`、`spfield`)、`tsol.f` (`spdrift_solenoid`)。DRIFT の `RADIUS > 0` | pysad v2: `spch.Solver3D`。係数と壁の Green 関数を導出。メッシュ上の離散化は未導出。**K64 は電場が $1/\gamma$ 倍に弱い (B89)** ほか (B90〜B92、B98) |
| **PSPAC** (2 次元 PIC) | §7.4。矩形の箱の正弦級数の解 | `tspch_.c` (`tpspac`)、`spch.c`。`SPCH` 要素 | pysad v2: `spch.Solver2D`。連続の解と自由空間の一致を確認、差分の離散化は未導出。K64 は既定 (`PSPACNZ = 1`) で落ちる (B93)、周の終わりのキックが入らない (B95) |
| **WSPAC** | §3〜§4 | `tintrb.f` (`twspac`、`twspace`、`qwspac`、`twspfu`、`twspu`) | 導出・確認済み。B127、B133、縦の力の基準 (§4.2) |
| **ウェイク** | §5〜§6 | `txwake.f`、`tffswake.f`、`tfwake.f` | 導出・確認済み |
| **Touschek** | §8 | `tintrb.f` (`touckf`、`tintrb` の `caltouck`)、`Packages/emit.n` (`TouschekLifetime`) | 損失率の式と `touckf` の関数形を導出・確認。`emit.n` の絶対値の換算は未確認。v2 は (17.32) の型が既定 (B162) |
| **IBS** (INTRA) | §9 | `tintrb.f` (`tintrb`、`fintrb`、`fintrbu`) | 拡散行列の構造と係数を導出、核を数値で確認。Coulomb 対数・係数 $h_0^4$ は未確認 |
| **CSR** (縦) | §10。自由空間の定常ウェイクを導出。追跡の写像は §10.4 | `csrtrack.f` (`csrtrack`、`csrwake`、`csrrho`)、`csroy.f` (CSR の安定性の解析)、`Packages/CSR.n` (インピーダンス) | pysad: `csr.py`。ウェイクの導出と数値確認 (自由空間・定常)。ダクトのインピーダンスは未導出。K64 と照合済み ([docs/15_csr.md](../docs/15_csr.md)) |
| **定常分布の方程式** | §11 | `csrtrack.f` (`csrhaissin`) | 導出し、K64 の Newton 解・追跡と比較 |

## 13. 確認スクリプトと出力

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `coll_01_spacecharge_wake.py` | (a) 実験室系の力 $= eE'/\gamma$ (sympy)、(b) WSPAC のキック (横・$\delta$・線形のレンズ) と 2 次元ポテンシャルの数値積分、(c) ウェイク: ガウス + 指数関数の平均の損失 (縦)、全体のずれの平均キック (横) | [coll_01_spacecharge_wake.out.txt](outputs/coll_01_spacecharge_wake.out.txt) |
| `coll_02_touschek.py` | (a) Møller の断面積の整理、(b) 損失の立体角 $W(\kappa)$、(c) 運動量の平均のモンテカルロと $J$、(d) pysad の `piwinski_slice_rate` との係数、(e) 平坦ビームの極限と `touckf` | [coll_02_touschek.out.txt](outputs/coll_02_touschek.out.txt) |
| `coll_03_laslett.py` | WSPAC のチューンシフトと (17.14) の和 | [coll_03_laslett.out.txt](outputs/coll_03_laslett.out.txt) |
| `coll_04_csr_haissinski.py` | (a)〜(d) CSR のスペクトル・インピーダンス・ウェイクの係数と平均の損失、(e) Fokker-Planck の定常解、(f)(g) Haissinski 方程式の解と追跡 | [coll_04_csr_haissinski.out.txt](outputs/coll_04_csr_haissinski.out.txt) |
| `coll_05_spac_pspac.py` | (a) 一様な円柱、(b) 壁の基準の電位、(c) 3 次元ガウス分布の電位、(d) 矩形の箱の級数解、(e) 線形のレンズの積分とチューンシフトの行列、(f) 円筒の Green 関数の加法定理 | [coll_05_spac_pspac.out.txt](outputs/coll_05_spac_pspac.out.txt) |
| `coll_06_ibs.py` | IBS の拡散の核 (17.36) と (17.35) の平均のモンテカルロ | [coll_06_ibs.out.txt](outputs/coll_06_ibs.out.txt) |
| `coll_07_wake_slices.py` | 切片の粗視化の誤差の粒子数依存 | [coll_07_wake_slices.out.txt](outputs/coll_07_wake_slices.out.txt) |

実行: リポジトリのルートから `py -3.12 -W ignore physics_guide/scripts/coll_01_spacecharge_wake.py` (数秒) など。`coll_04` は追跡を含み数分、`coll_02`・`coll_05`・`coll_06`・`coll_07` は 1 分前後。

## 14. Fortran の対応

| 式 | ファイル:行 | サブルーチン | 変数・内容 |
|---|---|---|---|
| (17.5)、(17.8)〜(17.9)、(17.12) | `tintrb.f:378-455` | `twspac` | `f`、`az`、`twspfu`、`fx`、`fy`、`fu`、`dg`、`dpx`、`dpy`、`pr1`、`h1`、`dv` |
| (17.6)、(17.7) | `tintrb.f:395-407` | `twspac`、`twspace` | `xx1`、`yy1`、`xy1`、`a`、`c1`、`s1`、`sigx`、`sigy` |
| (17.13) | `tintrb.f:281-337` | `twspace` | `f`、`akx`、`aky`、`aks`、`akd`、`trans(2,1)` など |
| 光学への挿入 | `tintrb.f:339-376` | `qwspac` | `trs`、`coup` |
| (17.17)〜(17.21) | `txwake.f:5-170` | `txwake` | `ws`、`xs`、`ys`、`zs`、`wz`、`wx`、`wy`、`fwp`、`dwz`、`dwx`、`dwy` |
| 切片の定義 | `txwake.f:171-243` | `txwdefslice` | `izs`、`itab`、`mb` |
| `fw` | `tturn.f:259` | `tturn1` | `fw = (|charge| e pbunch anbunch/amass)/np0*.5` |
| (17.32)、(17.33) | `tintrb.f:6-19`、`:135-190` | `touckf`、`tintrb` | `tf` (`:179`)、`touckl`・`toucke` (`:184-185`)、`touckm`、`ptrans`、`vol` |
| (17.35)、(17.36) | `tintrb.f:199-246` | `tintrb`、`fintrbu` | `eig`、`pl`、`f1`〜`e3` (`:206-213`)、`bmax`・`bmin`・`ci` (`:223-231`)、`bmi` |
| CSR の追跡 | `csrtrack.f:229-407` | `csrtrack` | `ndiv`、`r56n`、`ad`、`bd`、`awn`、`rho`、`ddrho` |
| (17.49) | `csrtrack.f:82-227` | `tfcsrhaissin`、`csrhaissin` | `aw`、`u`、`u0`、`f`、`g`、`wi`、`zi` |

(行番号は参照時点のもの。)

## 15. pysad の対応

| 関数 | 対応する Fortran | 注意 |
|---|---|---|
| `spch.WSpace` (`kick`、`linear_k`) | `twspac`、`twspace` | v2 は縦の力の基準 `rwall` (`WSPACRADIUS`) とバンチの速度での位置 `dz` を追加 (K64 は粒子の $z$ そのまま。差は $\delta/\gamma^2$)。`selfcod` |
| `spch.gauss_field`、`gauss_potential` | `twspfu`、`twspu` | 楕円ガウス分布の複素誤差関数と数値積分 |
| `spch.Solver3D`・`Solver2D` | `spkick.f`・`tspch_.c`+`spch.c` | v2 で正しい係数に書き直し (K64 の誤りは §17) |
| `wake.WakeKick`、`wake._kick` | `txwake` | K64 の癖をそのまま残す (§17)。切片の方法も直訳 |
| `wake.kick_factor` | `tturn.f:259` | $\mathtt{fw}$ |
| `csr.py` | `csrtrack.f`、`csroy.f` | K64 と照合済み |
| `touschek.py` | `tintrb.f` (IBS・Touschek) | (17.32) の型 (`piwinski_slice_rate`、v2 の既定) と K64 の `touckf` (k64) |

## 16. 使った条件

- WSPAC (確認 (b)): 20 MeV の電子 ($\gamma = 40.139$、$\beta = 0.99969$)、$N = 10^{10}$、$\sigma_z = 5$ mm、区間の長さ 0.3 m、$\sigma = 1$ mm (丸い)、$\sigma_x = 1.5$ mm・$\sigma_y = 0.5$ mm (楕円)、$z$ は $0$、$\pm3$〜$6$ mm、$x, y$ は $0.4$〜$2.5$ mm。壁の半径は 0 (K64 と同じく軸上で電位 0)。
- チューンシフト (§4.4): 0.5 GeV の電子、8 セル、周長 35.3 m、$N = 2\times10^{10}$、CAVI 0.2 MV・harmonic 20、RFSW・RAD・RADCOD、追跡は放射なしで、$d = 10^{-7}$ の試験粒子 4 本の 1 周の行列から $\cos$ の平均でチューン。
- ウェイク (確認 (c)): 1.7 GeV の電子、20000 個のガウス分布 ($\sigma_z = 5$ mm、種 5)、$N = 10^{10}$、$W_0 = 10^{12}$ V/C ($W_{0t} = 10^{13}$ V/C/m)、$\ell = 2, 5, 20$ mm。表は 4001 点 (0.05 mm 刻み、0〜0.2 m)。バンチ数 `nb` = 1。粗視化の粒子数依存 (§6.2) は同じ設定で、$n_p = 2500$〜160000、種 1〜3。
- 空間電荷の確認 (§7): 単位系は $\epsilon_0$ または $Q/4\pi\epsilon_0 = 1$。3 次元ガウス分布は $\sigma_x = 1$、$\sigma_y = 0.6$。箱は正方形、中心に丸いガウス分布。
- Touschek (§8): 運動量 $\lambda_1 = 4\times10^{-6}$、$\lambda_2 = 10^{-6}$ ($(mc)^2$)、$c = 1.2\times10^{-3}$ (モンテカルロ 400 万個)。係数の確認: 3 GeV の電子、$\sigma_x = 100\ \mu$m、$\sigma_y = 10\ \mu$m、$\sigma_z = 5$ mm、局所の角度の広がり $10^{-4}$・$3\times10^{-5}$、$N = 10^{10}$、$a = 0.02$、区間 1 m・周長 1 m。
- IBS (§9): 固有値は $p_0$ 単位の相対値 (1, 0.3, 0.1) など 4 組。サンプル 600 万個。
- CSR (§10): 単位 $Ne^2/(4\pi\epsilon_0)$、$R = \sigma = 1$。Haissinski (§11): $\mu_s = 0.04$、$\sigma_E = 10^{-3}$、$\sigma_{z0} = 3$ mm、$d = 4\times10^{-3}$、$\ell = 1$ mm、格子 0.15 mm・512 点、20000 粒子、3000 周。

## 17. 未確認・限界・K64 の疑い

**未確認**:
- SPAC・PSPAC: メッシュ上の離散化 (形状関数、`spkick.f` の円筒の Green 関数の離散化、`spch.c` の差分の次数)。導出は連続の式まで。
- WSPAC: $\sigma_z$ や $\delta$ が大きいときの粗視化、ビーム行列の自己無撞着な収束 (K64 は `Emittance[]` の平衡値で、空間電荷を入れた平衡は反復で求める。docs/09)、$z$ に比例する平均のずれ ($\Sigma_{15}/\Sigma_{55}$) を入れないことの影響、試験粒子の $\delta$ による力の変化 ($O(\delta/\gamma^2)$)。チューンシフトの $x$ の 0.9% の差の原因 (結合の影響の見積りは未確認)。
- ウェイク: 光学 (`WAKEOPT`)、複数バンチ (`nb > 1`)、実機のウェイク関数 (確認は人工的な指数関数だけ)、表が $z=0$ から始まるときの数値 (§6.2 の `else` 枝は読みだけ)。(17.17) の $\beta\to1$ の近似 ($\Delta\delta = \Delta E/(mc^2p_0)$、正しくは $/(\beta p_0)$)。
- Touschek: 相対論的な補正 (静止系の運動量が 1 に近い)、Born 近似の位相因子 (しきい値付近)、静止系の縦の運動量の広がり、バンチ内の $\delta$ の分布。`emit.n` の `TouschekLifetime` の絶対値の換算 (`touckl` 表から寿命へ)。K64 の `touckf` の定数 ($a$、$b$、$\varepsilon_E$) は当てはめで、導出していない。
- IBS: Coulomb 対数の定数 2 と $b_{\max}$・$b_{\min}$ の取り方、$h_0^4$ と $p_0^4$ の差、結合・分散がある場合の $P$ (`pl` の作り方)、ローレンツ変換 `trans1`・`trans2` の導出、`tintraconv` の収束。
- CSR: ダクトの遮蔽・過渡的なウェイクのインピーダンス (`CSR.n`)、遅延ポテンシャルからの直接の導出 (本章は放射のスペクトルと因果律)、`csroy.f` の安定性の解析 (作用角の格子の線形化)。
- Haissinski: 方程式の安定性、有限の $\mu_s$・$d$ の補正、強い結合での解の一意性。

**疑い (K64。ソースの読みと一部の数値)**:
- **WSPAC の縦の力の電位の基準が軸上** (`fu` が中心で 0): 縦の空間電荷の主要項 ($\ln(b/\sigma)$ に比例する、壁の半径に依る項。(17.11)) が入らない。K64 は軸上の縦の力が 0 になる。pysad v2 は `WSPACRADIUS` で基準を壁の半径にできる。K64 の仕様か誤りかは未決 (ヘルプに説明が無い。物理的には (17.10) の定数項が正しい基準)。
- **B127**: WSPAC の `twspace` が最後のスライスで `gammab(nlat+1) = 0` を読み、レンズが 8 倍強い (docs/09 に記載。v2 はこの誤りを入れない)。
- **B133**: WSPAC ON の CALC・CALC6D は入口以外の Twiss が NaN (gs で再現、原因は未調査。v2 は CALC6D と RFSW ON の CALC に空間電荷を入れた)。
- **B89〜B92・B98**: SPAC の力が $1/\gamma$ 倍に弱い (電場の $\gamma$ の抜け。gs での実測は docs/10 §9.7、本章では導出の係数 (17.23) と比べる形)、キック後の $p_x$ の二重の割り算、MULT の `RADIUS` の呼び出しの引数 (B91)、メッシュの整数の割り算、少数の粒子で落ちる (docs/10、docs/30。v2 は書き直し)。
- **B93〜B96**: PSPAC の既定で落ちる、線密度の配分、周の終わりのキックの欠落、static 変数 (docs/30)。
- **B97**: `tspac.f` の `tspac` は死んだコードで、係数 $1/(p_0(1+p^2)) = 1/(\beta\gamma\cdot\gamma^2)$ は $\beta$ が 1 つ足りない (呼ばれない)。
- **ウェイクの癖** (docs/13): (i) 追跡の表の最後の 1 行が使われない (`twxiwp` の `lwl0` は行数 − 1)、(ii) `tturn.f` で `dy = ROTATE` の代入、(iii) CAVI・MULT の分割で、wake の蹴りの合計が $(1.5n-0.5)/n$ 倍 (v2 は台形の重みで合計 1)、(iv) `txwdefslice` の配列の外への書き込みの可能性。
- **B162**: Touschek の損失率の式 `touckf` は平坦ビームの極限 (§8.4。導出で確認) で、2 方向の運動量の広がりがある場合は適用外 (誤りというより近似。v2 は (17.32) の型が既定)。
- **IBS の $h_0^4$**: 導出は $p_0^4 = (\beta\gamma)^4$ (§9.3)。電子では差が無視できる。疑いとして記録するだけで、K64 の意図は未確認。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| `WSPAC` のフラグの説明 | **原文で確認**: 「strong-weak モード。ビームの大きさは `EMIT` + `CODPLOT` で。ガウス分布を仮定、`PBUNCH`」。`SELFCOD` の説明 (力の中心を粒子の重心にするか閉軌道にするか)。`PSPAC` は PIC で、`SPAC`・`WSPAC` と別物 (ヘルプの注意書き)。 |
| `LWAKE`・`TWAKE`・`WAKEOPT` | **原文で確認**: 縦ウェイクつきの光学、TRPT の下で WAKEOPT が有効。`WakeFunction` は TRPT + INS + (TWAKE or LWAKE) のときだけ有効、追跡は `TrackParticles` だけ。 |
| `INTRA`・`MINCOUP`・`TouschekLifetime` | **原文で確認**: `MINCOUP` の式 $\varepsilon_{x,y} = \max[\varepsilon^0_{x,y}, \mathtt{MINCOUP}(\varepsilon^0_x+\varepsilon^0_y)]$、`TouschekLifetime[nx, ny, nz]` の 3 通りの使い方。 |
| ヘルプに無いもの | (17.1)〜(17.49)、IBS・Touschek・CSR・定常分布の方程式の式。ウェイク関数の符号の約束 (`dz` = 源 − 試験粒子)、`fw` の規格化は K64 のソースから。 |

**参考にした標準的な文献 (書名)**: 『Handbook of Accelerator Physics and Engineering』(空間電荷、ウェイク、Touschek 効果、バンチ内散乱、CSR の項目)、『Classical Electrodynamics』(Lorentz 変換、シンクロトロン放射のスペクトル、静電場)。本章の式は、これらの標準の結果を出発点として、本章で導出・数値確認したもので、文献の式との行ごとの照合は行っていない。

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
| `xx1`、`yy1`、`xy1` | $z$ を固定した横の分散 (17.6) | `twspac` | |
| `c1`、`s1` | 主軸の回転 $\cos\theta$、$\sin\theta$ | `twspac` | `c1`、`s1` |
| `sigx`、`sigy` | 主軸の $\sigma_x$、$\sigma_y$ (17.7) | `twspac` | `sg1`、`sg2` |
| `fx`、`fy`、`fu` | $\nabla U$ と $U$ | `twspfu` | |
| `dg` | $\Delta\delta$ ((17.9) で $U_{\rm ref}=0$) | `twspac` | |
| `selfcod` | 力の中心: 重心 (真) か閉軌道 (偽) | `ffs_flag` (`SELFCOD`) | `selfcod` |
| `wakel(2,*)`、`waket(2,*)` | 縦・横のウェイク関数の表 (z, w) | `ffs_wake` | `wl`、`wt` |
| `lwl`、`lwt` | 表の行数 − 1 (追跡で使う行数) | `txwake` | |
| `itab(np)`、`izs(np)` | $z$ の昇順の並び、切片の番号 | `txwake`、`txwdefslice` | `itab`、`izs` |
| `ws`、`xs`、`ys`、`zs` | 切片の粒子数、座標の和、重心 | `txwake` | |
| `wz`、`wx`、`wy` | 切片 $n$ が受ける和 (17.20) | `txwake` | |
| `fw`、`fwp` | (17.21) の係数 | `tturn.f:259`、`txwake` | `fw` |
| `zw0l`、`zw0t` | 表の最初の $z$ と 0 の小さい方 (自己項の判定) | `txwake` | |
| `dzwr`、`dzlim` | バンチの幅、再整列の閾値 ($10^{-5}$) | `wakez` | `DZLIM` |
| `touckf(x)` | $H(\varepsilon)$ の当てはめ (17.33) | `tintrb.f:6-19` | `touschek_f` |
| `tf`、`ptrans`、`vol` | $1/(V^*\lambda_1^{3/2})$ の係数、運動量の広がり、$V^*$ | `tintrb` | |
| `eig`、`pl`、`r` | 運動量の共分散の固有値、行列、固有ベクトル (17.36) | `tintrb` | |
| `bmax`、`bmin`、`ci`、`cintrb` | Coulomb 対数の上下限、拡散行列の係数 | `tintrb`、`temit.f:801` | |
| `fintrbu(u,a1,a2)` | 角度積分の被積分関数 | `tintrb.f` | `fintrbu` |
| `aw`、`inten`、`ndiv`、`r56`、`r65` | CSR のキックの強さ、分割数、$r_{56}$、$r_{65}$ | `csrtrack.f`、`CSR.n` | `aw`、`ndiv` |
