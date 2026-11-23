# 13. 放射: kinematical method

**この章で分かること**: SAD の放射 (追跡の `tradk`、CALC の `tradke`) が、要素の種類に依らず、粒子の向きの変化から曲率半径を決めて古典の放射の式を使う「kinematical method」であること。(1) 曲げ角 $\phi$ と曲率半径 $\rho_r$ の定義 (ヘルプ)、(2) 光子数 $N_\gamma$・臨界エネルギー $u_c$・平均損失・分散の古典の式と、K64 の定数 `anrad`・`cuc`・`cave`・`cuu` の対応、(3) 1 要素の分割 ($N_\gamma \lesssim 1$)、(4) 放射による運動量の縮み、横運動量・位置・$z$・$dv$ の更新、(5) FLUC のポアソン過程と乱数 `drndsr`、(6) CALC の励起 (`tradke` の $D$)、RADCOD、(7) 放射積分。平衡エミッタンスは [14_emittance.md](14_emittance.md)。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) ($p$、$\delta$、$dv$)、[10_cavity.md](10_cavity.md) ($z$ の更新で時間を保つ)。

**読む順**: §1 (結論) → §3 (古典の式) → §4 (追跡 1 スライス) → §6 (CALC の励起)。

## 1. 結論

- **kinematical method** (ヘルプ): 区間 (1, 2) で粒子の向き $\bm q = (p_x, p_y, p_z)/p$ が $\bm q_1 \to \bm q_2$ に変わったとき、一様な曲げと仮定して $\sin|\phi| = |\bm q_2 \times \bm q_1|$、$\rho_r = (L_{12} - z_2 + z_1)/|\phi|$ ($L_{12}$ は区間の名目の長さ、$z$ は経路の遅れ)。要素の種類 (多極、ソレノイド、フリンジ、電場でも) を問わず、同じ 1 つのルーチンで扱える。
- **古典の式** ($\gamma$: ローレンツ因子、$\lambda_C = \hbar/mc$ 換算コンプトン波長、$\alpha$ 微細構造定数): 臨界エネルギー $u_c = \tfrac32\lambda_C\,mc^2\gamma^3/\rho$、1 区間の光子数 $N_\gamma = \tfrac{5}{2\sqrt3}\alpha\gamma\phi$、平均エネルギー $\langle u\rangle = \tfrac{8}{15\sqrt3}u_c$、$\langle u^2\rangle = \tfrac{11}{27}u_c^2$。K64 の定数は `cuc` $= 1.5\lambda_C$ (単位 m)、`anrad` $= 5\alpha/(2\sqrt3)$、`cave` $= 8/(15\sqrt3)$、`cuu` $= 11/27$ (数値で確認した)。
- **確認** ([rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt)): 3 GeV の電子、偏向磁石 1 台 (ANGLE = $2\pi/16$、$L = 1.2$ m、$\rho = 3.056$ m、$u_c = 19.6$ keV、$N_\gamma = 24.3$) の追跡で、平均損失は $4.8850\times10^{-5}$ ($\delta$。古典の $\tfrac23 r_emc^2\gamma^4\theta/\rho$ は $4.8852\times10^{-5}$、差 $5\times10^{-5}$)。FLUC ON の標準偏差は $2.046\times10^{-5}$ (理論 $\sqrt{N_\gamma\cdot11/27}\,u_c = 2.055\times10^{-5}$。標本誤差 $10^{-7}$)。軸からずれた四極でも同じ式で $0.99993$ 倍 (§4.4)。16 台の和は $U_0 = 2.3448$ MeV で、$C_\gamma E^4 I_2/2\pi = 2.3449$ MeV。
- **1 スライスの運動量・位置の更新** (`tradk1`): 平均損失 $de = \mathtt{cave}\,N_\gamma\,u_c$、$dg = -de/(1+2de)$、$p_x \mathrel{+}= -\tfrac12 dp_x\,dg$、$x \mathrel{+}= (-\tfrac12dp_x\,dg)\,al/3$ ($dp_x$ はスライス内の $p_x$ の変化)。後ろの 2 式 (1/2 と 1/3) は、スライス内で一様に放射する連続過程の積分から導かれる (数値積分で確認、§4.3)。$dg$ の式は、曲率半径一定のときの連続過程の解 $(1 + 3de)^{-1/3} - 1$ と 2 次まで一致する (§4.2)。$dv$ の式は K64 では 01 章 (1.9) と分母の第 2 項が違い、低エネルギーで 1 に比べて小さくない差になる (§4.1 項目 5、付録 B)。CALC の `tradke` は追跡と式が一部違い、K64 の `tradke` そのままでは物理的に正しくない箇所がある (§6、付録 B)。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| $\bm q$ | 運動量の向き | `ppx ppy ppz` (外積の成分)、`px py` | `track6._tradk1` |
| $\phi$ (`theta`) | 曲げ角。$\sin\phi$ = `ppa` $= |\bm q_1\times\bm q_2|$ | `theta = asinz(ppa)` | |
| $\rho_r$ | 曲率半径。`rho = al1/theta` | `al1 = al - z + zr00` (区間の長さ $-$ $z$ の変化) | |
| $N_\gamma$ | 区間の光子数の期待値 | `anp = anrad*h1*theta` | |
| $u_c$ | 臨界エネルギー ($p_0c$ 単位) | `uc = cuc*h1**3/p0/rho`、または `uct = cuc*h1**3/p0/al1`、`uc = uct*theta` | |
| $\gamma$ | | `h1 = p2h(p)`、`p = p0*pr` | |
| `anrad` | $\tfrac{5}{2\sqrt3}\,r_c/\lambda_C = \tfrac{5}{2\sqrt3}\alpha$ | `tfsetparam.f`: `anrad = 5/2/sqrt(3)*rcratio`、`rcratio` = $r_e/(\hbar c/mc^2) = \alpha$ | |
| `cuc` | $\tfrac32 r_e/\mathtt{rcratio} = \tfrac32\lambda_C$ | `tfsetparam.f` | |
| `cave`・`cuu` | $\langle u\rangle/u_c$、$\langle u^2\rangle/u_c^2$ | `photons.f:381`: `cave = 8/15/sqrt(3)`、`cuu = 11/27` | |
| `gmin` | $g$ の下限 | `-0.999999` | |
| `urad`・`crad`・`erad` | 放射の定数 ($p_0$ で規格化したもの) | `tfsetparam.f` | 未確認 |

## 3. 古典の式

### 3.1 臨界エネルギーと光子数

曲率半径 $\rho$ の円運動 (向心加速度 $\gamma^2c^2/\rho$ の相対論的な電荷) の放射は、パワー (Liénard) $P = \tfrac23 r_emc^2\,\gamma^4c/\rho^2$、スペクトル $dP/du = (P/u_c)\,\dfrac{9\sqrt3}{8\pi}\,F(u/u_c)$ ($F(\xi) = \xi\int_\xi^\infty K_{5/3}(t)\,dt$、$K_{5/3}$ は変形ベッセル関数)、臨界エネルギー $u_c = \tfrac32\hbar c\gamma^3/\rho$ である。$F$ の積分は

$$
\int_0^\infty F\,d\xi = \frac{8\pi}{9\sqrt3},\qquad
\int_0^\infty\frac{F}{\xi}\,d\xi = \frac{5\pi}{3},\qquad
\int_0^\infty \xi F\,d\xi = \frac{55\pi}{81}
\tag{13.1}
$$

(数値積分で $1.61226610$、$5.23598776$、$2.13318020$。右辺と 8 桁で一致。)

**(13.1) の導出**。$G(\xi) = \int_\xi^\infty K_{5/3}(t)\,dt$ とおくと $F = \xi G$。部分積分で $\int_0^\infty\xi^{s-1}G\,d\xi = \bigl[\xi^sG/s\bigr]_0^\infty + \dfrac1s\int_0^\infty\xi^sK_{5/3}(\xi)\,d\xi = \dfrac1s\int_0^\infty t^sK_{5/3}(t)\,dt$ (境界項は、$\xi \to \infty$ で $G$ が指数的に減るので 0、$\xi \to 0$ では $K_{5/3}(t) \sim t^{-5/3}$ から $G \sim \xi^{-2/3}$ で、$\xi^sG \sim \xi^{s-2/3} \to 0$ ($s = \mu + 2 \ge 1$))。ベッセル関数のメリン変換 $\int_0^\infty t^sK_\nu(t)\,dt = 2^{s-1}\Gamma\bigl(\tfrac{1 + s + \nu}2\bigr)\Gamma\bigl(\tfrac{1 + s - \nu}2\bigr)$ を $\nu = 5/3$ で使う。$\int_0^\infty\xi^\mu F\,d\xi = \int\xi^{\mu+1}G\,d\xi$ ($s = \mu + 2$) なので

$$
\int_0^\infty\xi^\mu F\,d\xi = \frac{2^{\mu+1}}{\mu + 2}\,\Gamma\Bigl(\frac\mu2 + \frac73\Bigr)\Gamma\Bigl(\frac\mu2 + \frac23\Bigr)
$$

$\mu = 0$: $\Gamma(7/3)\Gamma(2/3) = \tfrac43\cdot\tfrac13\Gamma(\tfrac13)\Gamma(\tfrac23) = \tfrac49\cdot\dfrac{\pi}{\sin(\pi/3)} = \dfrac{8\pi}{9\sqrt3}$ ($\Gamma(x)\Gamma(1-x) = \pi/\sin\pi x$)。$\mu = -1$: $\Gamma(11/6)\Gamma(1/6) = \tfrac56\Gamma(\tfrac56)\Gamma(\tfrac16) = \tfrac56\cdot\dfrac\pi{\sin(\pi/6)} = \dfrac{5\pi}3$。$\mu = 1$: $\dfrac43\Gamma(\tfrac{17}6)\Gamma(\tfrac76) = \dfrac43\cdot\tfrac{11}6\cdot\tfrac56\cdot\tfrac16\,\Gamma(\tfrac56)\Gamma(\tfrac16) = \dfrac43\cdot\dfrac{55}{216}\cdot2\pi = \dfrac{55\pi}{81}$。数値積分とガンマ関数の式は 12 桁で一致する ([rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) (a))。

光子のエネルギー分布は $dN/du \propto (dP/du)/u$ ($dN = dP\,dt/u$)。規格化すると、1 個の光子のエネルギー ($\xi = u/u_c$) の平均と 2 次モーメントは

$$
\frac{\langle u\rangle}{u_c} = \frac{\int F\,d\xi}{\int (F/\xi)\,d\xi} = \frac{8\pi/(9\sqrt3)}{5\pi/3} = \frac{8}{15\sqrt3} = 0.30792,\qquad
\frac{\langle u^2\rangle}{u_c^2} = \frac{\int\xi F\,d\xi}{\int(F/\xi)\,d\xi} = \frac{55\pi/81}{5\pi/3} = \frac{11}{27}
\tag{13.2}
$$

(`cave`、`cuu` に等しい。数値: $0.30792014$、$0.40740741$。) 光子の発生率は、パワーを平均エネルギーで割ったもの:

$$
\dot N = \frac{P}{\langle u\rangle}
= \frac{\tfrac23 r_emc^2\gamma^4c/\rho^2}{\tfrac{8}{15\sqrt3}\cdot\tfrac32\hbar c\gamma^3/\rho}
= \frac{5}{2\sqrt3}\,\frac{r_emc^2}{\hbar c}\,\frac{\gamma c}{\rho}
= \frac{5}{2\sqrt3}\,\alpha\gamma\,\frac{c}{\rho}
\tag{13.3}
$$

($r_emc^2/\hbar c = \alpha$ を使った。係数: $\dfrac23\cdot\dfrac{15\sqrt3}{8}\cdot\dfrac23 = \dfrac{5\sqrt3}{6} = \dfrac{5}{2\sqrt3}$。) 曲げ角 $\phi$ の区間は時間 $\rho\phi/c$ かかるので $N_\gamma = \tfrac{5}{2\sqrt3}\alpha\gamma\phi$。これが `anp = anrad*h1*theta`。

**単位の確認**: $u_c$ を $p_0c$ (設計運動量 × $c$) で割ると $u_c/(p_0c) = \tfrac32\lambda_C\,\gamma^3/(\rho\,p_0)$ ($mc^2\gamma^3/\rho \cdot \lambda_C$ が $E$、$p_0$ は $mc$ 単位)。`uc = cuc*h1**3/p0/rho` と一致 (`cuc` $= \tfrac32\lambda_C$)。3 GeV、$\rho = 3.056$ m で $u_c/p_0c = 6.533\times10^{-6}$ ($19.60$ keV)。

### 3.2 平均損失と分散

1 区間で、光子数はポアソン分布 (平均 $N_\gamma$) に従い、各光子のエネルギーは独立にスペクトルから引く。放出エネルギー $\Delta u$ の平均と分散 (ポアソン過程の複合):

$$
\langle\Delta u\rangle = N_\gamma\langle u\rangle,\qquad
\mathrm{Var}(\Delta u) = N_\gamma\langle u^2\rangle
\tag{13.4}
$$

(平均 $N_\gamma$ のポアソン分布の合計 $\sum_{i=1}^{n}u_i$ の分散は $\langle n\rangle\langle u^2\rangle$。) $\gamma^4$ が出る確認: $\langle\Delta u\rangle = \tfrac{5}{2\sqrt3}\alpha\gamma\phi\cdot\tfrac{8}{15\sqrt3}\cdot\tfrac32\lambda_Cmc^2\gamma^3/\rho = \tfrac23\alpha\lambda_C\,mc^2\gamma^4\phi/\rho = \tfrac23 r_emc^2\gamma^4\phi/\rho$ (古典の 1 区間のエネルギー損失)。

## 4. 追跡 (`tradk1`、`tradk`、`tradkf1`)

### 4.1 1 スライスの流れ (FLUC OFF: `tradk1`)

1. スライスの入口と出口の向きから、外積 $\bm q_2\times\bm q_1$ の成分 (`ppx ppy ppz`) を作り、$\phi$ (`theta`)、$\rho_r$ (`al1 = al - z + zr00`、`rho = al1/theta`、$z$ はスライス内の経路の遅れの変化) を決める (§1)。
2. 平均損失 $de = \mathtt{cave}\cdot N_\gamma\cdot u_c$ ($u_c$ は `uct*theta`)。
3. 運動量: $dg = -de/(1 + 2\,de)$、$g \leftarrow \max(\mathtt{gmin},\,g + dg)$ ($g = \delta$)。
4. 横運動量・位置: `ddpx = -0.5*dpx*dg`、`px += ddpx`、`x += ddpx*al1/3` ($y$ も同様)。
5. $dv = -g(1 + p_r)/[h_2(h_2 + p_2)] + \mathtt{dvfs}$ ($p_r = 1 + g$、**$p_2 = p_0p_r$** (`p2=p0*pr`、$mc$ 単位の運動量))、$z \leftarrow z\,\dfrac{p_2/h_2}{p_1/h_1}$ (速度の比 $\beta_2/\beta_1$ をかける。時間を保つ: [10_cavity.md](10_cavity.md) の式と同じ)。

**$dv$ の式は 01 章 (1.9) と同じではない**: (1.9) と 10 章 §4 の $dv = -\delta(1 + p_r)/[h(h + p_rh_0)]$ は、$h = \sqrt{1 + p_r^2p_0^2}$ について厳密な $1 - v/v_0$ で、分母の第 2 項は $p_rh_0$ ($h_0 = \gamma_0$)。`tradk1` (`photons.f:1199-1202`) は $p_rh_0$ の代わりに $p_2 = p_0p_r$ を使う。差は $p_r(h_0 - p_0) = p_r/(h_0 + p_0)$ で、$h \simeq 2h_0$ に対して相対 $O(1/\gamma^2)$。数値 ([rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) (c)): $p_0 = 1\,mc$ ($0.51$ MeV/c)、$p_r = 1.1$ で、厳密な $1 - \beta/\beta_0 = -0.046433$、(1.9) は $-0.046433$ (一致)、`tradk1` の式は $-0.054613$ (相対 18% 違う)。$p_0 = 50\,mc$、$\delta = -0.01$ で相対 $10^{-4}$、3 GeV では (1.9) と `tradk1` の式の差は相対 $7 \times 10^{-9}$ で見えない。pysad は `track6.py:397` で「SAD のまま ($h_2 + p_2$)」と写す。この差は**低エネルギーで効く K64 の癖の疑い** (ソースの読みと数値例で確定。docs/30 に B 番号は無い。gs での確認はしていない)。

追跡の運動量変数は $q = p_x/p$ (01 章) なので、光子が速度の向きに出るとき $q$ は変わらず、$p_x = q\,p$ だけが $p$ に比例して縮む。

### 4.2 $dg = -de/(1+2de)$

$de$ は、損失をスライスの入口の $\gamma$ ($h_1$) で計算した量 ($de \propto \gamma^4$)。実際には、放射で $p$ が縮むと $\gamma$ も縮み、$\gamma^4$ 依存の分だけ損失は減る。**導出**: kinematical method では、スライスの間、曲率半径 $\rho_r$ (と曲げ角) は軌道の幾何から決まった一定値で、$p$ だけが変わると仮定する。古典の損失は $dp/ds \propto -\gamma^4/\rho_r^2 \propto -p^4$ ($\gamma \propto p$、$\gamma \gg 1$)。入口の値で規格化した $u = p/p_1$、$\tau = s/L \in [0, 1]$ として、入口の $\gamma$ での損失が $de$ なので

$$
\frac{du}{d\tau} = -de\,u^4,\qquad u(0) = 1\ \Longrightarrow\ u^{-3} = 1 + 3\,de\,\tau,\qquad dg = u(1) - 1 = (1 + 3de)^{-1/3} - 1
$$

($u^{-4}du = -de\,d\tau$ の積分。) 展開すると $(1 + 3de)^{-1/3} - 1 = -de + 2de^2 - \tfrac{14}3de^3 + \cdots$ (二項展開 $(1 + x)^{-1/3} = 1 - x/3 + \tfrac29x^2 - \tfrac{14}{81}x^3$、$x = 3de$)。K64 の $dg = -de/(1 + 2de) = -de(1 - 2de + 4de^2 - \cdots) = -de + 2de^2 - 4de^3 + \cdots$。**2 次まで一致し、3 次 ($-14/3$ 対 $-4$) から違う**。したがって `dg = -de/(1+2de)` は、$de^2$ までの近似で、$dg$ を $dg = -de(1 + 2dg)$ の陰的な形 ($dg(1 + 2de) = -de$) に書き換えたものと読める。数値 ([rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) (b)。$du/d\tau = -de\,u^4$ の ODE との比較): $de = 5 \times 10^{-5}$ で K64 との差 $8 \times 10^{-14}$、$de = 0.01$ で $6 \times 10^{-7}$、$0.05$ で $6 \times 10^{-5}$、$0.2$ で $2 \times 10^{-3}$ (3 次の差)。例の $|de| \sim 5 \times 10^{-5}$ では 3 次の差は見えない。ただし $\gamma \propto p$ の近似 ($1/\gamma^2$ の補正) は、この導出に含まれる。

### 4.3 横運動量の補正 ($-\tfrac12dp_x\,dg$) と位置の補正 ($al/3$)

スライス内で、力 (磁場の衝撃) が一様に加わり、運動量の損失も一様に ($g(s) = dg\cdot s/L$) 起こるとする。力による運動量の衝撃は損失に依らず ($\mathrm{d}p_x = F\,ds$、全体で $dp_x$)、運動量は $p(s) = 1 + g(s)$ なので、傾き $q = p_x/p$ の変化は

$$
\frac{dq}{ds} = \frac{F}{p(s)}\simeq F\,(1 - g(s)),\qquad F = \frac{dp_x}{L},\qquad g(s) = dg\,\frac{s}{L}
\tag{13.5}
$$

これを $s$ で積分すると $q(s) = dp_x\left[\dfrac{s}{L} - \dfrac12\,dg\left(\dfrac{s}{L}\right)^2\right]$。放射のないときの $dp_x\,s/L$ に対する補正は $-\tfrac12dp_x\,dg\,(s/L)^2$ である。したがって

$$
\Delta q(L) = -\tfrac12\,dp_x\,dg,\qquad
\Delta x(L) = \int_0^L\left(-\tfrac12dp_x\,dg\right)\left(\frac{s}{L}\right)^2ds = -\tfrac12\,dp_x\,dg\cdot\frac{L}{3}
\tag{13.6}
$$

(位置の補正の 1/3 は $\int_0^1t^2dt$。) `ddpx = -0.5*dpx*dg`、`x += ddpx*al1/3` と一致する。

**確認** ([rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt) (f)): 運動方程式を数値積分 (相対許容誤差 $10^{-12}$) すると、$dp_x = 2\times10^{-3}$、$dg = -10^{-3}$、$L = 1$ で $q(L) - dp_x = 1.0007\times10^{-6}$ (公式 $+1.0000\times10^{-6}$)、$x(L) - dp_xL/2 = 3.335\times10^{-7}$ (公式 $3.333\times10^{-7}$)。差は $dg$ の 2 次。

### 4.4 どんな要素でも同じ式

軸からずれた四極は、曲率 $1/\rho = K_1x(s)$ で、曲げ角は軌道の傾きの変化で決まる。ずれ $x_0$ の粒子の損失を、四極内の軌道 $x(s) = x_0\cos(\sqrt k\,s)$ ($k = K_1/L$) で古典の式 $\tfrac23r_emc^2\gamma^4\int(1/\rho)^2ds$ に積分した値と比べる ([rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt) (e)): $x_0 = 1, 2, 4$ mm で、追跡と古典の比は $0.99993$ (3 つとも同じ。薄い近似の誤差)。偏向磁石の場の 1 つの式で書く代わりに、粒子の向きの変化だけで決めるので、多極・ソレノイド・フリンジでも同じルーチンが使える (ヘルプ)。

### 4.5 FLUC (`tradk`、`tradkf1`)

FLUC では、スライスごとに光子数 $n$ をポアソン分布 (平均 `anp`) から引き (`tdusrn`)、$n$ 個の光子のエネルギーを `drndsr` で引く。

`tdusrn`: $P(n) = e^{-N}N^n/n!$ を累積して、一様乱数 $p$ に対し $p \le \sum_{k\le n}P(k)$ となる最小の $n$ を選ぶ (`p0 = exp(-anp)`、`x = x*anp/an`、`pn = pn + x*p0`)。累積が `cutoff` $= 0.9999$ を超えたら打ち切る (`npmax` $= 10^4$ 個まで)。出力は $\sum u_i/u_c$ (`dph`)、位置の重みつき平均 (`r1`、`r2`: 光子がスライス内のどこで出たか) と光子数 `an`。

`drndsr` (有理式による近似): 一様乱数 $p$ から、$P(y) = \int_0^y pp(y')dy'$ ($pp(y) = \dfrac{3}{5\pi}\int_y^\infty K_{5/3}(t)dt$、$y = u/u_c$) の逆関数 $y(p)$ を、3 つの区間の有理式で与える (近似の相対誤差 $< 2\times10^{-4}$)。

**確認** ((c)): 乱数 $2\times10^6$ 個の標本の平均と 2 次モーメントは $0.30791$、$0.40787$ (コメントの値 $0.30792$、$0.40741$。標本誤差 $4\times10^{-4}$ の範囲)。(d) の FLUC ON の標準偏差 $2.046\times10^{-5}$ は理論 (13.4) $2.055\times10^{-5}$ と標本誤差の範囲で一致。(h): 初期の $p_x = 3\times10^{-3}$ でも $\mathrm{Var}(\delta) = 4.206\times10^{-10}$ (理論 $4.223\times10^{-10}$)。

`photons` (光子を 1 個ずつ記録する) の `tphrec` は、光子のエネルギー、位置、向きを記録する付随の機能で、本章の対象外。

### 4.6 分割: $N_\gamma \lesssim 1$

放射を 1 つの点の損失で表すのは、スライスが十分に細かい ($N_\gamma \lesssim 1$) とき。K64 は `tbende.f:462-465` で、偏向磁石の分割数に次の項を含める:

$$
\mathrm{ndiv} = \max\Bigl(\mathrm{ndiv}_0,\ \mathrm{int}(\mathrm{nrad}\cdot\mathtt{emidiv}\cdot\mathtt{emidib}),\ 1 + \mathrm{int}\bigl(|\phi_bh_0\,\mathtt{anrad}|/\mathtt{epsr1}/10^6\cdot\mathtt{emidiv}\cdot\mathtt{emidib}\bigr)\Bigr)
$$

($\mathrm{nrad} = 1 + \mathrm{int}\bigl(|a_0\rho_{bc}/\mathtt{epsr1}\cdot\mathtt{crad}\,(h_0b)^2|\bigr)$。) `epsr1` $=$ `epsrad` $= 10^{-6}$ (`toplvl.f:291` の parameter。EPS の指定があればその倍)、`EMITDIV`・`EMITDIVB` の既定は 1 (`initb1.f:1703-1706`)。したがって**第 3 項は $1 + \mathrm{int}(N_\gamma)$** ($N_\gamma = \phi_bh_0\,\mathtt{anrad}$ は偏向 1 台の光子数 (13.3))、すなわち 1 刻みあたり $N_\gamma \lesssim 1$。第 2 項 `nrad` は、`crad` $= r_{\rm classic}(c/\mathtt{amass})^2/p_0/1.5$ ($\tfrac1{1.5} = \tfrac23$。Liénard の $\tfrac23r_e$ の係数) と $(h_0b)^2$ から、1 台の相対損失を `epsrad` で割った数に当たると読める (1 刻みの相対損失を $10^{-6}$ 以下にする)。例 (3 GeV、1.2 m の偏向、$N_\gamma = 24.28$): 第 3 項は 25 ([rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) (e))、追跡の相対損失 $4.885 \times 10^{-5}$ を $10^{-6}$ で割ると 48.9 で、第 2 項はこれに近い数 (約 50) になるはず。`crad` の単位まで含めた実際の値の計算、`tsolque.f:622-623`・`tsteee.f:79-82` の同じ形の式の確認、実際の分割数の確認は行っていない (**未確認**)。

## 5. 運動量の縮みと RADCOD

放射で $p$ が縮むと、リングでは RF が補う。`radcod` (RADCOD ON) は、閉軌道を求めるとき放射の損失を含めて固定点を探す (12 章 §3)。放射のないとき、閉軌道の $z$ は RF の同期位相で決まる。

放射テーパー (`radtaper`) は、放射による $p_0(s)$ の変化 (設計運動量が下流に向かって縮む) を取り入れる。この章では扱わない (読みだけ)。

## 6. CALC の励起 (`tradke`): $D$ と $B$

**先に注意** (K64 の `tradke` は追跡の `tradk1` と式が違う):
- 平均損失は `dg = -cave*anp*uc` (1 次。追跡は §4.2 の $dg = -de/(1 + 2de)$)。
- **横運動量の反作用の符号**: `tradke` は `ddpx = +.5*dpx*dg` (`photons.f:1417`)、追跡は `ddpx = -.5*dpx*dg` ($dg < 0$)。§4.3 の導出 (運動量が縮むと同じ磁場が運動量の小さい粒子をより強く曲げるので $-\tfrac12dp_x\,dg$) は追跡の符号を支持し、**`tradke` の符号が K64 の誤り** ([docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md) の B111。RADCOD ON の閉軌道が追跡で閉じない。v2 は `tradke` の符号を反転して修正、k64 モードは K64 のまま)。
- 微分 (減衰の行列) の節に 4 点の誤り・近似がある (B138、§6.3)。SOL 内の `dbeam(8)` の添字の誤記 (B39、§6.2)。

以下は式の構造の説明で、**K64 の `tradke` の式をそのまま写したものではなく、上の符号は追跡と §4.3 の導出に合わせた形**で書く。

### 6.1 平均損失

CALC (`tradke`) は、追跡 (`tradk1`) と同じ幾何で、$\theta$、$\rho_r$ ($al_1$)、$h_1$、$N_\gamma = \mathtt{anrad}\,h_1\theta$、$u_c = \mathtt{cuc}\,h_1^3\theta/(p_0\,al_1)$ を作り、平均損失 $dg = -\mathtt{cave}\cdot N_\gamma\cdot u_c$ (`dg = -cave*anp*uc`)、その位相空間の微分 (減衰の行列 $B$: `tr1`)、励起の行列 $D$ (`dee`、`dbeam`) を作る。

### 6.2 励起 (`dee`、`dbeam`)

**励起** (`dee = anp*uc**2*cuu`): 区間の $\delta$ の分散 $\mathrm{Var}(\Delta u) = N_\gamma\langle u^2\rangle$ ((13.4)、$p_0c$ 単位)。光子は速度の向きに出るので、$q$ は変わらず、$p_x = q\,p$ が $p$ に比例して変わる: 粒子 1 個について $\Delta p_x = p_x\,\Delta\delta/p_r \simeq p_x\Delta\delta$、$\Delta p_y \simeq p_y\Delta\delta$。したがって、粒子の集団についてスライス内で平均した

$$
\langle\Delta\delta^2\rangle = \mathtt{dee},\qquad
\langle\Delta p_x\,\Delta\delta\rangle = \langle p_x\rangle\,\mathtt{dee},\qquad
\langle\Delta p_x\Delta p_y\rangle = \langle p_xp_y\rangle\,\mathtt{dee},\qquad
\langle\Delta p_x^2\rangle = \langle p_x^2\rangle\,\mathtt{dee}
\tag{13.7}
$$

ここで $\langle\cdot\rangle$ は粒子の集団と、スライス内の位置 ($s$) の両方についての平均。$p_x$ は軌道 ($p_{xi}$ から $p_x$ まで線形に変わるとする) と、集団の広がり $\delta p_x$ の和なので、$\langle p_x\rangle = (p_{xi} + p_x)/2$、$\langle p_x^2\rangle = \dfrac{p_{xi}^2 + p_{xi}p_x + p_x^2}3 + \langle\delta p_x^2\rangle$ (軌道の $s$ 平均 $\tfrac1L\int p_x(s)^2ds$ に、ビーム行列の $\langle\delta p_x^2\rangle$ が加わる)。Fortran (`photons.f:1506-1516`): $pxm = p_{xi} + p_x$ として

- `dbeam(21) = dee` ($\langle\delta^2\rangle$)、`dbeam(17) = pxm*dee*0.5` ($\langle p_x\delta\rangle$)、`dbeam(19) = pym*dee*0.5` ($\langle p_y\delta\rangle$)。
- `dbeam(3) = (beam(3) + b*(2*beam(5) + b*beam(6)) + (pxm**2 + pxi**2 + px**2)/6)*dee`: $(pxm^2 + p_{xi}^2 + p_x^2)/6 = [(p_{xi} + p_x)^2 + p_{xi}^2 + p_x^2]/6 = (p_{xi}^2 + p_{xi}p_x + p_x^2)/3$ が軌道の $s$ 平均。**`beam(3)`** $= \langle p_x^2\rangle$ (ビーム行列の $\langle xx\rangle = $ `beam(1)`、$\langle xp_x\rangle = $ `beam(2)`、$\langle p_xp_x\rangle = $ `beam(3)`、$\langle xy\rangle = $ `beam(4)`、…の 21 成分の並び) は、集団の広がりの $\langle\delta p_x^2\rangle$。ソレノイドの磁場 $b_{zh}$ があるとき、力学的運動量は $p_x + b\,y$ ($b = b_{zh}/2$) なので $\langle(p_x + by)^2\rangle = \langle p_x^2\rangle + 2b\langle p_xy\rangle + b^2\langle y^2\rangle$ (`beam(3) + b*(2*beam(5) + b*beam(6))`)。
- `dbeam(10)` (py の対応する式)。`dbeam(8)` $= \langle\Delta p_x\Delta p_y\rangle$ は $\langle(p_x + by)(p_y - bx)\rangle = \langle p_xp_y\rangle - b\langle xp_x\rangle + b\langle yp_y\rangle - b^2\langle xy\rangle$ のはずだが、K64 は `beam(2) - beam(10) + b*beam(4)` と $\langle yp_y\rangle$ (`beam(9)`) の代わりに $\langle p_y^2\rangle$ (`beam(10)`) と書く (B39。$b_{zh} = 0$ では効かない)。

したがって励起が「軌道の $p_x$ だけで決まる」のではなく、**集団の広がり**にも依存する ($\langle\Delta p_x^2\rangle$ の $\langle\delta p_x^2\rangle\,\mathtt{dee}$ の項)。

この励起の行列を、1 周の $\Sigma = M\Sigma M^\top + \Delta$ に組む話は [14_emittance.md](14_emittance.md)。

**確認**: (13.7) の $\langle\Delta\delta^2\rangle$ の部分 (`dee`) は (d)(h) の標準偏差・分散で確認した。$\langle\Delta p_x\Delta\delta\rangle$、$\langle\Delta p_x^2\rangle$ の係数は、その後の輸送で位相空間の相関が混ざる (偏向磁石内の $R_{26}$ による $\delta$ と $p_x$ の相関: 初期 $p_x$ ゼロでも出る) ので、追跡の統計から直接は読めない。**数値では未確認**。上は (13.7) の導出 ($\Delta p_x = p_x\Delta\delta$ と軌道・集団の平均) と、Fortran の式の項ごとの読みである。

### 6.3 減衰 ($B$): 連鎖律

減衰の行列は、平均損失 $dg = -\mathtt{cave}\,N_\gamma\,u_c$ の位相空間の微分。$N_\gamma = \mathtt{anrad}\,h_1\theta$、$u_c = \mathtt{cuc}\,h_1^3\theta/(p_0\,al_1)$ なので $dg = -\mathtt{cave}\,\mathtt{anrad}\,\mathtt{cuc}\,h_1^4\theta^2/(p_0\,al_1)$ で、$\theta$ ($d\theta$ = `dtheta`)、$h_1$ ($p_r$ 経由。$h_1 = \sqrt{1 + p_0^2p_r^2}$ から $d\ln h_1/dp_r = p_0^2p_r/h_1^2 = $ `dh1r`)、$al_1$ ($d\,al_1$ = `dal`) に依存する。連鎖律で

$$
dN_\gamma = \mathtt{anrad}\,(h_1\,d\theta + \theta\,h_1\,\mathtt{dh1r}\,dp_r),\quad
du_c = u_{ct}\,d\theta - u_c\,\frac{d\,al_1}{al_1} + 3u_c\,\mathtt{dh1r}\,dp_r,\quad
d(dg) = -\mathtt{cave}\,(dN_\gamma\,u_c + N_\gamma\,du_c)
$$

($u_{ct} = u_c/\theta$。) コードの `danp = anrad*h1*dtheta`、`danp(6) += anp*dh1r`、`duc = uct*dtheta - uc*dal/al1`、`duc(6) += 3*uc*dh1r`、`ddg = -cave*(danp*uc + anp*duc)` は、この連鎖律そのままである (sympy で、$dg$ の $\theta$、$al_1$、$p_r$ での偏微分とコードの式の差が 0: [rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) (d))。$\theta = \arcsin(x_{pa})$、$x_{pa} = |\vec x_p|/p_r^2$ ($\vec x_p$ は外積 `xpx, xpy, xpz`) の微分は $d\theta = dx_{pa}/\sqrt{1 - x_{pa}^2}$ (`das`)、$dx_{pa} = \dfrac{\vec x_p\cdot d\vec x_p}{x_{pa}\,p_r^4} - \dfrac{2x_{pa}}{p_r}dp_r$ (sympy で確認、同 (d))。横運動量の更新 $\Delta p_x' = \tfrac12dp_x\,dg$ (符号は上のとおり) の微分は $d(\Delta p_x') = \tfrac12[(dp_x^{\rm end} - dp_{xr0})\,dg + dp_x\,d(dg)]$。

**K64 の `tradke` の微分の 4 点の誤り・近似** (B138): (1) 区間の始点の座標に関する微分 `tr1` を、終点の座標に関する微分を要する $D$・$B$ の更新にそのまま使う (正しくは $\mathtt{tr1}\cdot\mathtt{tr2}^{-1}$)。(2) `dddpx` の第 2 項が `ddpx*ddg` ($= \tfrac12dp_xdg\,d(dg)$) で、上の $dp_x\,d(dg)$ ではない。(3) `tr1(1,:)` の $x' = x + \Delta p_x'\,c_1$ ($c_1 = al_1/(3p_r)$) の $c_1$ の微分が無い。(4) `dxpa` の第 1 項が `/xpa/pr**2` で、上の $/(x_{pa}p_r^4)$ ではない ($\delta \ne 0$ の軌道で $2\delta$ の相対誤差)。(2)(4) は上の連鎖律との比較で確認した。(1)(3) は docs/30 の記述による (ここでは導出していない)。v2 は修正済み (`exactd`)、k64 モードは K64 のまま。**`tr1` の各成分を 1 つずつ手で導出し、数値微分と照合することはしていない** (**未確認**)。

## 7. 放射積分 (RADINT)

リングの 1 周の平均損失と、励起・減衰の積分は、偏向磁石の曲率 $1/\rho$ の積分で書かれる (標準):

$$
I_2 = \oint\frac{ds}{\rho^2},\qquad U_0 = \frac{C_\gamma E^4}{2\pi}I_2,\qquad C_\gamma = \frac{4\pi}{3}\frac{r_e}{(mc^2)^3} = 8.846\times10^{-5}\ \mathrm{m/GeV^3}
\tag{13.8}
$$

例のリング (16 台の扇形偏向、$I_2 = 16\phi^2/L = 2.0562$ /m) で $U_0 = 2.3449$ MeV。追跡 (16 台の和) は $2.3448$ MeV。pysad の `lat.radiation_integrals()` は $I_2 = 2.05617$ で一致する。他の積分 ($I_1 = 11.665$、$I_3 = 0.6729$、$I_4 = 1.2493$、$I_5 = 0.6486$ (単位は各積分による)) は、分散 $\eta$、$\mathcal{H}$ を含み、14 章で平衡エミッタンスに使う。ここでの確認は $I_2$ のみ。

SAD の RADINT は kinematical method (本体のずれた四極なども含める) の CALC (`tradke`) とは別の、偏向磁石の式による積分で、両者の違い (磁場のずれた四極の放射、K64 の B06・B07) は付録 B。

## 8. Fortran の対応

| 内容 | Fortran | 式 | pysad |
|---|---|---|---|
| 定数 | `tfsetparam.f` `tphyzp` (`cuc`、`anrad`、`urad`、`crad`)、`photons.f:381` (`cave`、`cuu`、`gmin`) | (13.1)〜(13.3) | `track6` の定数 |
| 平均の放射 (追跡) | `photons.f` `tradk1` | §4.1 | `track6._tradk1` |
| ポアソン + スペクトル (追跡) | `photons.f` `tradkf1`・`tradkfn`、`drndsr.f` `tdusrn`・`drndsr` | §4.5 | `track6._tradkf1` |
| 分割と呼び出し | `photons.f` `tradk` (`tradkn`) | §4.6 | |
| CALC の励起・減衰 | `photons.f` `tradke` (`dee`、`dbeam`、`tr1`) | (13.7) | `emit6d` |
| 放射積分 | `INTGRL.f` | (13.8) | `radint.compute`、`lat.radiation_integrals` |

## 9. 確認スクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `rad_02_derivations.py` | (a) (13.1) のメリン変換、(b) $dg$ の導出と ODE、(c) $dv$ の分母、(d) `tradke` の連鎖律 (sympy)、(e) 分割数 | [rad_02_derivations.out.txt](outputs/rad_02_derivations.out.txt) |
| `rad_01_kinematical.py` | (a) 定数の対応、(b) 放射スペクトルのモーメント (13.1)(13.2)、(c) `drndsr` のモーメント、(d) 偏向磁石 1 台の平均・分散、(e) ずれた四極の損失、(f) 横運動量・位置の補正 (13.6)、(g) $U_0 = C_\gamma E^4I_2/2\pi$、(h) $\delta$ の分散 | [rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/rad_01_kinematical.py` (数秒〜十数秒)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、偏向磁石 $L = 1.2$ m、ANGLE = $2\pi/16$ (扇形、端面なし、$\rho = 3.056$ m、$u_c = 19.6$ keV、$N_\gamma = 24.3$)、四極 $L = 0.3$ m、K1 = 0.35。RF なし (RFSW OFF)、$\delta = 0$ で初期条件を 0 (損失は 1 要素の通過)。FLUC は標本 $2\times10^4$ 個。

**未確認**:
- 分割数の `crad` を含む第 2 項の実際の値 (§4.6)。
- (13.7) の係数の数値確認、減衰の行列 `tr1` の各成分の手による導出と数値微分との照合 (§6)。$dg$ の式 (§4.2)、`dee`・`dbeam` の構造と連鎖律 (§6) は導出した。ソレノイドがあるときの項は読みだけ。
- RADCOD・RADTAPER の分岐、`tradkn` (多粒子の同時処理)、スピンの項 (`sprot`)、`tphrec` (光子の記録)。
- 軸からずれた四極の、$\gamma$ が大きい極限での古典放射の式の適用範囲 (臨界エネルギーに比べて十分短い、半古典 (量子) の補正は無視)。
- $\gamma$ が小さい (低エネルギー) とき、$N_\gamma$ の式、$u_c$ の式の適用。3 GeV 以外のエネルギーは確認していない。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| synchrotron-radiation の $\sin|\phi|$、$\rho_r$ | **原文で確認** (TeX 版・平文版とも)。`tradk1` の `ppa`・`al1 = al - z + zr00`、`rho = al1/theta` と一致。 |
| 「A component is sliced so that $N_\gamma\lesssim1$」 | **原文で確認**。分割数の式との対応は §4.6 のとおり未確認。 |
| 「derivatives ... obtained kinematically using the transfer matrix」 | **原文で確認**。`tradke` の `tr2`、`ddpz` など。 |
| 「In the region where the field is not uniform, such as the F1 region of a BEND, a special treatment for $\rho_r$」 | **原文で確認**。この章では扱わない (05 章、06 章)。 |
| ヘルプに無いもの | (1) $N_\gamma$・$u_c$・$\langle u\rangle$・$\langle u^2\rangle$ の式と定数 (本章 §3)。(2) 1 スライスの運動量・位置の更新 (§4)。(3) 励起の行列 (§6)。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

- **`tradke` の `ddpx` の符号が追跡と逆** (B111): `tradke` は $+\tfrac12dp_x\,dg$、追跡 (`tradk1`・`tradkn`) は $-\tfrac12dp_x\,dg$。§4.3 の導出は追跡の符号を支持する (v2 は `tradke` を反転。k64 モードは K64 のまま)。
- **`tradke` の微分の 4 点の誤り** (B138、§6.3) と、SOL 内の `dbeam(8)` の添字 (B39、§6.2)。
- **`tradk1` の $dv$ の分母が (1.9) と違う** ($h_2 + p_2$ 対 $h_2 + p_rh_0$。§4.1 項目 5): 相対 $O(1/\gamma^2)$ で低エネルギーで効く。ソースの読みと数値例で確定、B 番号なし、gs では未確認。
- 放射積分の計算は kinematical method とは独立な、偏向磁石の曲率による積分で、四極・6 極のずれた軌道の放射を含まない (pysad の `radiation_integrals` の docstring)。リングの平衡エミッタンスは `tradke` による 6 次元の $D$、$B$ (14 章) で、放射積分とは別の計算。K64 の B06・B07 (RADINT) は [docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md)。この章では確認していない。

## 付録 C. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `anrad` | $\tfrac{5}{2\sqrt3}\alpha$ | `tfsetparam.f` `tphyzp` | |
| `cuc` | $\tfrac32\lambda_C$ [m] | `tfsetparam.f` | |
| `rcratio` | $r_e/\lambda_C = \alpha$ | `tfsetparam.f` | |
| `cave`・`cuu` | $8/(15\sqrt3)$、$11/27$ | `photons.f:381` | |
| `gmin` | $g$ の下限 $-0.999999$ | `photons.f:381` | |
| `anp` | $N_\gamma$ = `anrad*h1*theta` | `photons.f` | |
| `uc`・`uct` | $u_c$ ($p_0c$ 単位)、$u_c/\phi$ | `photons.f` | |
| `theta` | $\phi$ = `asinz(ppa)` | `photons.f` | |
| `ppx ppy ppz` | $\bm q_2\times\bm q_1$ の成分 | `photons.f` | |
| `al1` | $L - z + z_{r0}$ (区間の長さ − $z$ の変化) | `photons.f` | |
| `de` | 平均損失 $\mathtt{cave}\cdot N_\gamma\cdot u_c$ | `tradk1` | |
| `dg` | $= -de/(1+2de)$ ($\delta$ の変化) | `tradk1` | |
| `dpx`・`dpy` | スライス内の $p_x$、$p_y$ の変化 (入口の $p_x$ を回して出口と引く) | `tradk1` | |
| `ddpx` | $-\tfrac12dp_x\,dg$ | `tradk1` | |
| `dee` | $N_\gamma\,\mathtt{cuu}\,u_c^2$ ($\delta$ の分散の増加) | `tradke` | |
| `dbeam` | 励起の増分 (ビーム行列 21 成分の並び) | `tradke` | |
| `tr1`・`tr2` | 減衰・伝搬の行列 | `tradke` | |
| `cutoff`・`npmax` | ポアソンの累積の打ち切り ($0.9999$)・光子数の上限 ($10^4$) | `drndsr.f` | |
