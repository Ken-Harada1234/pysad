# 13. 放射: kinematical method

**この章で分かること**: SAD の放射 (追跡の `tradk`、CALC の `tradke`) が、要素の種類に依らず、粒子の向きの変化から曲率半径を決めて古典の放射の式を使う「kinematical method」であること。(1) 曲げ角 $\phi$ と曲率半径 $\rho_r$ の定義 (ヘルプ)、(2) 光子数 $N_\gamma$・臨界エネルギー $u_c$・平均損失・分散の古典の式と、K64 の定数 `anrad`・`cuc`・`cave`・`cuu` の対応、(3) 1 要素の分割 ($N_\gamma \lesssim 1$)、(4) 放射による運動量の縮み、横運動量・位置・$z$・$dv$ の更新、(5) FLUC のポアソン過程と乱数 `drndsr`、(6) CALC の励起 (`tradke` の $D$)、RADCOD、(7) 放射積分。平衡エミッタンスは [14_emittance.md](14_emittance.md)。

**前提**: [01_variables_and_conventions.md](01_variables_and_conventions.md) ($p$、$\delta$、$dv$)、[10_cavity.md](10_cavity.md) ($z$ の更新で時間を保つ)。

**読む順**: §1 (結論) → §3 (古典の式) → §4 (追跡 1 スライス) → §6 (CALC の励起)。

## 1. 結論

- **kinematical method** (ヘルプ): 区間 (1, 2) で粒子の向き $\bm q = (p_x, p_y, p_z)/p$ が $\bm q_1 \to \bm q_2$ に変わったとき、一様な曲げと仮定して $\sin|\phi| = |\bm q_2 \times \bm q_1|$、$\rho_r = (L_{12} - z_2 + z_1)/|\phi|$ ($L_{12}$ は区間の名目の長さ、$z$ は経路の遅れ)。要素の種類 (多極、ソレノイド、フリンジ、電場でも) を問わず、同じ 1 つのルーチンで扱える。
- **古典の式** ($\gamma$: ローレンツ因子、$\lambda_C = \hbar/mc$ 換算コンプトン波長、$\alpha$ 微細構造定数): 臨界エネルギー $u_c = \tfrac32\lambda_C\,mc^2\gamma^3/\rho$、1 区間の光子数 $N_\gamma = \tfrac{5}{2\sqrt3}\alpha\gamma\phi$、平均エネルギー $\langle u\rangle = \tfrac{8}{15\sqrt3}u_c$、$\langle u^2\rangle = \tfrac{11}{27}u_c^2$。K64 の定数は `cuc` $= 1.5\lambda_C$ (単位 m)、`anrad` $= 5\alpha/(2\sqrt3)$、`cave` $= 8/(15\sqrt3)$、`cuu` $= 11/27$ (数値で確認した)。
- **確認** ([rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt)): 3 GeV の電子、偏向磁石 1 台 (ANGLE = $2\pi/16$、$L = 1.2$ m、$\rho = 3.056$ m、$u_c = 19.6$ keV、$N_\gamma = 24.3$) の追跡で、平均損失は $4.8850\times10^{-5}$ ($\delta$。古典の $\tfrac23 r_emc^2\gamma^4\theta/\rho$ は $4.8852\times10^{-5}$、差 $5\times10^{-5}$)。FLUC ON の標準偏差は $2.046\times10^{-5}$ (理論 $\sqrt{N_\gamma\cdot11/27}\,u_c = 2.055\times10^{-5}$。標本誤差 $10^{-7}$)。軸からずれた四極でも同じ式で $0.99993$ 倍 (§4.4)。16 台の和は $U_0 = 2.3448$ MeV で、$C_\gamma E^4 I_2/2\pi = 2.3449$ MeV。
- **1 スライスの運動量・位置の更新** (`tradk1`): 平均損失 $de = \mathtt{cave}\,N_\gamma\,u_c$、$dg = -de/(1+2de)$、$p_x \mathrel{+}= -\tfrac12 dp_x\,dg$、$x \mathrel{+}= (-\tfrac12dp_x\,dg)\,al/3$ ($dp_x$ はスライス内の $p_x$ の変化)。後ろの 2 式 (1/2 と 1/3) は、スライス内で一様に放射する連続過程の積分から導かれる (数値積分で確認、§4.3)。

## 2. 記号と単位

| 量 | 意味 | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| $\bm q$ | 運動量の向き | `ppx ppy ppz` (外積の成分)、`px py` | `track6._tradk1` |
| $\phi$ (`theta`) | 曲げ角。$\sin\phi$ = `ppa` $= |\bm q_1\times\bm q_2|$ | `theta = asinz(ppa)` | |
| $\rho_r$ | 曲率半径。`rho = al1/theta` | `al1 = al - z + zr00` (区間の長さ $-$ $z$ の変化) | |
| $N_\gamma$ | 区間の光子数の期待値 | `anp = anrad*h1*theta` | |
| $u_c$ | 臨界エネルギー ($p_0c$ 単位) | `uc = cuc*h1**3/p0/rho`、または `uct = cuc*h1**3/p0/al1`、`uc = uct*theta` | |
| $\gamma$ | | `h1 = p2h(p)`、`p = p0*pr` | |
| `anrad` | $\tfrac{5}{2\sqrt3}\,r_c/\lambda_C\cdot$ | `tfsetparam.f`: `anrad = 5/2/sqrt(3)*rcratio`、`rcratio` = $r_e/(\hbar c/mc^2) = \alpha$ | |
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
5. $dv = -g(1 + p_r)/[h_2(h_2 + p_2)] + \mathtt{dvfs}$ ($p_r = 1 + g$)、$z \leftarrow z\,\dfrac{p_2/h_2}{p_1/h_1}$ (速度の比 $\beta_2/\beta_1$ をかける。時間を保つ: [10_cavity.md](10_cavity.md) の式と同じ)。

追跡の運動量変数は $q = p_x/p$ (01 章) なので、光子が速度の向きに出るとき $q$ は変わらず、$p_x = q\,p$ だけが $p$ に比例して縮む。

### 4.2 $dg = -de/(1+2de)$

$de$ は損失をスライスの入口の $\gamma$ ($h_1$) で計算した量 ($de \propto \gamma^4$)。実際の損失は、放射で $p$ が縮むと、$\gamma$ も縮み、$\gamma^4$ 依存の分だけ減る。$dg$ を $dg = -de(1 + 2dg)$ の陰的な形に書いたもの ($dg\,(1+2de) = -de$ から $dg = -de/(1+2de)$) で、$de$ の 2 次の補正の取り入れ方 (終端の $p$ の寄与を 1 次で入れる) に見える。この形がなぜ正しいか (どの近似か) の導出は**行っていない** (未確認)。$|de| \sim 5\times10^{-5}$ 程度では補正は $10^{-4}$ の相対量で、(d) の一致 ($5\times10^{-5}$) には影響しない。

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

`drndsr` (Yokoya の近似式): 一様乱数 $p$ から、$P(y) = \int_0^y pp(y')dy'$ ($pp(y) = \dfrac{3}{5\pi}\int_y^\infty K_{5/3}(t)dt$、$y = u/u_c$) の逆関数 $y(p)$ を、3 つの区間の有理式で与える (近似の相対誤差 $< 2\times10^{-4}$)。

**確認** ((c)): 乱数 $2\times10^6$ 個の標本の平均と 2 次モーメントは $0.30791$、$0.40787$ (コメントの値 $0.30792$、$0.40741$。標本誤差 $4\times10^{-4}$ の範囲)。(d) の FLUC ON の標準偏差 $2.046\times10^{-5}$ は理論 (13.4) $2.055\times10^{-5}$ と標本誤差の範囲で一致。(h): 初期の $p_x = 3\times10^{-3}$ でも $\mathrm{Var}(\delta) = 4.206\times10^{-10}$ (理論 $4.223\times10^{-10}$)。

`photons` (光子を 1 個ずつ記録する) の `tphrec` は、光子のエネルギー、位置、向きを記録する付随の機能で、本章の対象外。

### 4.6 分割: $N_\gamma \lesssim 1$

放射を 1 つの点の損失で表すのは、スライスが十分に細かい ($N_\gamma \lesssim 1$) とき。K64 は `tbende.f` などで、分割数に $1 + \mathrm{int}(|\phi_bh_0\,\mathtt{anrad}|/\mathtt{epsr1}/10^6\cdot\mathtt{emidiv}\cdot\mathtt{emidib})$ のような項 ($\phi_bh_0\,\mathtt{anrad}$ は偏向 1 台の $N_\gamma$) を含める (`tbende.f:465`、`tsolque.f:623`、`tsteee.f:82`)。定数 `epsr1`、`emidiv`、`emidib` の意味と、この数が実際に $N_\gamma$ をどの程度に抑えるかは、この章では確認していない (**未確認**)。例 (3 GeV、1.2 m の偏向) では 1 台で $N_\gamma = 24.3$ 個で、1 スライス $\lesssim 1$ にするには 25 スライス以上が要る。

## 5. 運動量の縮みと RADCOD

放射で $p$ が縮むと、リングでは RF が補う。`radcod` (RADCOD ON) は、閉軌道を求めるとき放射の損失を含めて固定点を探す (12 章 §3)。放射のないとき、閉軌道の $z$ は RF の同期位相で決まる。

放射テーパー (`radtaper`) は、放射による $p_0(s)$ の変化 (設計運動量が下流に向かって縮む) を取り入れる。この章では扱わない (読みだけ)。

## 6. CALC の励起 (`tradke`): $D$ と $B$

CALC (`tradke`) は、追跡 (`tradk1`) と同じ幾何で、(i) 平均損失 $dg = -\mathtt{cave}\cdot\mathtt{anp}\cdot uc$ (`dg = -cave*anp*uc`)、(ii) その位相空間の微分 (減衰の行列 $B$: `tr1`、`ddg`、`danp`、`duc`、`dtheta`、`dal`)、(iii) 励起の行列 $D$ (`dee`、`dbeam`) を作る。

**励起** (`dee = anp*uc**2*cuu`): 区間の $\delta$ の分散 $\mathrm{Var}(\Delta u) = N_\gamma\langle u^2\rangle$ ((13.4)、$p_0c$ 単位)。光子は速度の向きに出るので、$q$ は変わらず、$p_x = q\,p$ が $p$ に比例して変わる ($\Delta p_x = q\,\Delta\delta$、$\Delta p_y = q_y\Delta\delta$)。したがって、スライス内の平均で

$$
\langle\Delta\delta^2\rangle = \mathtt{dee},\qquad
\langle\Delta p_x\,\Delta\delta\rangle = \langle p_x\rangle_s\,\mathtt{dee},\qquad
\langle\Delta p_x^2\rangle = \langle p_x^2\rangle_s\,\mathtt{dee}
\tag{13.7}
$$

($\langle\cdot\rangle_s$ はスライス内の一様平均。$p_x$ が入口 $p_{xi}$ から出口 $p_x$ まで線形に変わるとして、$\langle p_x\rangle_s = (p_{xi} + p_x)/2$、$\langle p_x^2\rangle_s = (p_{xi}^2 + p_{xi}p_x + p_x^2)/3$。) Fortran: `pxm = pxi + px` として `dbeam(17) = pxm*dee*0.5`、`dbeam(3) = (…+(pxm**2 + pxi**2 + px**2)/6)*dee` (括弧の中は $(p_{xi}+p_x)^2 + p_{xi}^2 + p_x^2 = 2(p_{xi}^2 + p_{xi}p_x + p_x^2)$ なので、$/6$ で (13.7) の $\langle p_x^2\rangle_s$ になる)。`dbeam(21) = dee`。ソレノイドの磁場 ($bzh$) があるときは、正準運動量と力学的運動量の差 $b = bzh/2$ の項が加わる。この励起の行列を、1 周の $\Sigma = M\Sigma M^\top + \Delta$ に組む話は [14_emittance.md](14_emittance.md)。

**確認**: (13.7) の $\langle\Delta\delta^2\rangle$ の部分 (`dee`) は (d)(h) の標準偏差・分散で確認した。$\langle\Delta p_x\Delta\delta\rangle$、$\langle\Delta p_x^2\rangle$ の係数の確認は、その後の輸送で位相空間の相関が混ざる (偏向磁石内の $R_{26}$ による $\delta$ と $p_x$ の相関: 初期 $p_x$ ゼロでも出る) ので、追跡の統計から直接は読めない。**未確認**。

**減衰** ($B$): 平均損失 $dg$ の位相空間の微分。$dg = -\mathtt{cave}\cdot N_\gamma(\phi)\,u_c$ で、$\phi$ ($\mathrm{d}\phi$ = `dtheta`)、$\gamma$ ($\propto p$、`dh1r`)、$\rho_r$ ($\mathrm{d}(\mathtt{al1})$ = `dal`) に依存する。連鎖律で、`danp = anrad*(h1*dtheta + theta*dh1r)`、`duc = uct*dtheta - uc*dal/al1 + 3 uc dh1r (δ 成分)`、`ddg = -cave*(danp*uc + anp*duc)`。これが 1 次の減衰の行列 `tr1` に入る。導出の確認は、この章では行っていない (**未確認**。`tests` の独立な放射の確認は pysad の試験にある)。

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
| `rad_01_kinematical.py` | (a) 定数の対応、(b) 放射スペクトルのモーメント (13.1)(13.2)、(c) `drndsr` のモーメント、(d) 偏向磁石 1 台の平均・分散、(e) ずれた四極の損失、(f) 横運動量・位置の補正 (13.6)、(g) $U_0 = C_\gamma E^4I_2/2\pi$、(h) $\delta$ の分散 | [rad_01_kinematical.out.txt](outputs/rad_01_kinematical.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/rad_01_kinematical.py` (数秒〜十数秒)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子 3 GeV ($p_0$ = 3 GeV/$c$)、偏向磁石 $L = 1.2$ m、ANGLE = $2\pi/16$ (扇形、端面なし、$\rho = 3.056$ m、$u_c = 19.6$ keV、$N_\gamma = 24.3$)、四極 $L = 0.3$ m、K1 = 0.35。RF なし (RFSW OFF)、$\delta = 0$ で初期条件を 0 (損失は 1 要素の通過)。FLUC は標本 $2\times10^4$ 個。

**未確認**:
- $dg = -de/(1+2de)$ の導出 (§4.2)、分割の定数 `epsr1`・`emidiv`・`emidib` (§4.6)。
- (13.7) の $\langle\Delta p_x\Delta\delta\rangle$・$\langle\Delta p_x^2\rangle$ の係数、減衰の行列 `tr1` の導出。ソレノイドがあるときの項。
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
