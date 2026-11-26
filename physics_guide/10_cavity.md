# 加速空洞 CAVI と TCAVI

**この章で分かること**: 高周波加速空洞 CAVI の追跡 (K64 の `tcav`) が、ΔE の式 (ヘルプ) から $\delta$・$z$・$p_x$ をどう更新するか、(1) 位相と時間 ($t = -z\,h_1/p_1$)、エネルギー利得と、RFSW・RADCOD・同期位相 ($\mathtt{phis}$、`offset`) の関係、(2) エネルギーが変わったあとの運動量・$z$・$dv$・横運動量の更新 (時間を保つ)、(3) 横方向の依存 (V1、V20、V11、V02) と Panofsky–Wenzel、非相対論の補正、(4) 長い空洞の分割 (重み `ws`) と、設計運動量 $p_0(s)$ が変わる場合 (TRPT)、(5) 端の効果 $H_f$、を式を省かずに追う。偏向キック型の空洞 TCAVI は、ヘルプに節が無いので、ソース `ttcav.f` から書く。MULT に RF を持たせる `tmultiacc` は、この章の CAVI と同じ手続きで、MULT の本体は [08_mult.md](08_mult.md)。

**読む順**: §1 (結論) → §3 → §4 (薄い空洞の写像) → §6 (長い空洞) → §8 (TCAVI) → §9 (確認)。

## 1. 結論

- **エネルギー利得** (ヘルプ): $\Delta E = -e\,(V + V_1x + \tfrac12V_{20}x^2 + \tfrac12V_{02}y^2 + V_{11}xy)\,\sin(\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$。追跡の中では $h = E/mc^2$ の変化 `dh`= $-v\,[\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{phis}]$ ($v = V/mc^2$、$\mathtt{phis} = \omega\,t_{rf0}$、$\mathtt{dphis} = \mathtt{phis} - \mathrm{PHI}$。RAD OFF のとき)。粒子の時間は位置から $t = -z\,h_1/p_1$ ($c = 1$) で決まる。
- **写像**: 薄い空洞は、時間 $t$ を保ったまま、(i) $h_2 = h_1 + \Delta h$、(ii) $p_2 = \sqrt{h_2^2 - 1}$ ($mc$ 単位)、(iii) $z_2 = -t\,p_2/h_2$、(iv) 正準運動量 $p_x, p_y$ は変えない (追跡の内部変数 $q = p/(1+\delta)$ は $p_1/p_2$ 倍)、で、追跡と一致する (差 $2 \times 10^{-13}$、$\delta$ と $z$ は桁落ちの範囲)。**端なしの薄い空洞はシンプレクティック** (数値微分のヤコビアンで $10^{-10}$。§4)。端 ($H_f$) つきの空洞は、K64 の `tcavfrin` が $q$ を再スケールしないため、**シンプレクティックでない (K64 の誤り。gs のヤコビアンで確認、v2 は修正 = B332。§7.2、付録 A)**。入口の端だけの写像の欠陥の大きさは $|\mathtt{dpr}| \simeq (h/p^2)\,|\Delta h|$ で、$p_x$ によらない (§7.2 の条件で、$x = 1$ mm、3 GeV/$c$ は $8 \times 10^{-10}$、0.6 MeV/$c$ は $6 \times 10^{-6}$)。
- **横方向**: $V_1, V_{20}, V_{11}, V_{02}$ は、電場の振幅の横依存。Panofsky–Wenzel の定理で、横のキックは $\Delta p_x = -\partial_xV\,\cos(\omega t - \mathtt{dphis})/(\omega p_0)$。非相対論の補正 $V_{20}, V_{02} \mathrel{+}= V(2\pi f/c)^2/(\gamma\beta)^2/4$ は、変形ベッセル関数 $I_0(\omega r/\beta\gamma c) = 1 + (\omega r/\beta\gamma c)^2/4 + \cdots$ の第 2 項 (数値確認: 2 MeV/c で ΔE の $x$ 依存が $1 + (\omega x/\beta\gamma c)^2/4$ に一致)。
- **長い空洞**: ドリフト−キック−ドリフトを、重みつきの長さ `al*ws(n)`・キックの割合 `ws(n)` で $\mathrm{ndiv}$ 回繰り返す (2 次の対称分割)。重みは $\sum\mathtt{ws} = 1$ で、設計運動量が空洞の中で増える (TRPT) とき、スライスごとに等比数列 $\mathtt{ws}(n) \propto (1+r)^{n/\mathrm{ndiv}}$ (数値確認: $e^{wl} = 1 + r$)。TRPT でない普通のリングでは $r = 0$ で一様 ($\mathtt{ws} = 1/\mathrm{ndiv}$)。
- **端**: 入口で $H_f = -(eV/L)\,[\sin(\omega t - \mathrm{dphi}) + \sin\mathrm{dphi} - \mathrm{offset}]\,(x^2 + y^2)/4\;\delta(s - s_0)$、出口は符号が逆 (ヘルプ)。Maxwell 方程式から導くと、$\beta = 1$ で同じ形の収束 (横の力 $F_r = -(e r/2)\,[\partial_sE_z + (\beta/c)\,\partial_tE_z]$) になる。**K64 の式は $1/\beta$ を略している** (K64 の誤り。低エネルギーで効く。K64 の式に $1/\beta$ が入っていないことは、gs と pysad の追跡で確認し、$1/\beta$ が物理的に正しいことは軸の近くの Maxwell の場のローレンツ力の数値積分と gs で確認した。v2 は直した = B330。§7.1)。
- **TCAVI**: 偏向キック $\Delta p_x = -K_0\sin(\omega t + \phi)$ と、Panofsky–Wenzel の相手の $\Delta h = K_0\,p_0\,\omega\,x\cos(\omega t + \phi)$ の組。シンプレクティック。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 電圧 [V] | `VOLT` ($V$)、`DVOLT` | `vc` ($= V + \mathrm{DVOLT}$、`tturn.f` の呼び出しで足す)。$v = \texttt{vc}/\texttt{amass}\cdot|\texttt{charge}|$ | `el.VOLT`、`v` |
| 周波数 | `FREQ` [Hz]、`HARM` (FREQ = 0 のとき $\omega = \omega_0\cdot\mathrm{HARM}/c$) | `w` ($= 2\pi\,\mathrm{FREQ}/c$ [rad/m]) | `w` |
| 位相 | `PHI`、`DPHI` | `phi`、`dphi`、`phic` ($= (\texttt{phi}+\texttt{dphi})\cdot\mathrm{sign}(\texttt{charge})$) | `phic` |
| 同期の時間 | (ヘルプの $t_s$) | `trf0` (global)、`phis` (TRPT のとき 0、autophi のとき `phic`、それ以外は $\texttt{w}\cdot\texttt{trf0}$。`tcav.f:35`)、`dphis` $= \texttt{phis} - \texttt{phic}$ | `trf0`、`phis`、`dphis` |
| 一定の項 | — | `offset`、`offset1` | `offset`、`offset1` |
| 横の依存 | `V1`、`V20`、`V11`、`V02` | `v1`、`v20`、`v11`、`v02` → `v1a`、`v20a`($= \tfrac12v_{20}/\texttt{amass}\ldots$)、`v11a`、`v02a` | `v1a` ほか |
| 設計運動量・エネルギー | (設計運動量 $p_0$、$h_0 = \gamma_0$) | `p0`、`h0` ($mc$ 単位) | `p0`、`h0` |
| 空洞の中の設計の変化 | (TRPT のとき $p_0(s)$) | `vnominal` $= v\sin(-\mathrm{PHI}\cdot\mathrm{sign})$、`he` $= h_0 + \texttt{vnominal}$、`pe` | `vnominal`、`he`、`pe` |
| 分割 | `L`、`FREQ` | `ndiv`、`ws(1:ndiv)`、`wsn` | `ndiv`、`ws` |
| 端 | `DISFRIN`、`FRINGE` | `fringe` (= DISFRIN が 0)、`mfring` | `fringe`、`mfring` |

**この表と本文に出る名前の注釈** (規約は [01 章](01_variables_and_conventions.md) §8。定義はソースで確認、章末の付録 B に表):

- `amass`: 粒子の静止エネルギー $mc^2$ [eV] (大域パラメータ `MASS`、`tfsetparam.f` の `amass=rgetgl1('MASS')`)。電圧 [V] を `amass` で割ると $mc^2$ 単位の無次元量になる ($e\,V/mc^2$、電荷の絶対値 `abs(charge)` が掛かる)。だから $v = \mathtt{vc}/\mathtt{amass}\cdot|\mathtt{charge}|$ (`tcav.f:34`) で、$h = E/mc^2$ と同じ単位で足せる。
- `charge`: 粒子の電荷 (電子 $-1$。大域パラメータ `CHARGE`)。符号は位相 `phic` に入る (`sign(charge)`)。
- `pgev`: 設計運動量 $p_0c$ [eV] (名前は gev だが単位は eV。`MOMENTUM` の値)。`tphyzp` (`tfsetparam.f:105`) は、`pgev` から $p_0$ = `pgev/amass`、$h_0$ = `p2h(p0)`、磁気剛性、放射の定数、周回の角周波数 $\omega_0$ など、設計運動量に従属する量を一括して作り直すサブルーチン。TRPT のとき、空洞を出たあとに `pgev=pe*amass` としてから `tphyzp` を呼ぶ (`tcav.f:202-203`) のは、設計運動量が空洞で変わったので、これらを新しい $p_e$ で更新するため。
- `rad`・`trpt`・`autophi`: 論理変数。`rad` は放射の有無 (大域の `RAD`、RADCOD。`tffs.f` で `fff%rad`)、`trpt` は輸送路 (TRPT。運動量が空洞で変わる宣言、`tffs.f:589`)、`autophi` は CAVI の `APHI` ≠ 0 (自動位相。`tturn.f:438`、`autophi=cmp%value(ky_APHI_CAVI) /= 0.d0`)。いずれも §3.3 で `offset` と `phis` の場合分けに使う。
- `mfring`・`frmd`: 端 (フリンジ) の有無を表す整数。`frmd` は入力の `FRINGE` (0, 1, 2, 3) で、`tffs.f:1447-1451` が `mfring` に直す: 要素が通常の向き (`cmp%ori`) なら `mfring = FRINGE`、逆向き (LINE の中で `-` をつけて逆向きに置いた要素) なら $\mathrm{frmd}\,(11 + \mathrm{frmd}\,(2\,\mathrm{frmd} - 9))/2$ で、FRINGE = 0, 1, 2, 3 がそれぞれ 0, 2, 1, 3 になる (代入で確かめる: 1 → $1\cdot(11 + 1\cdot(2 - 9))/2 = 2$、2 → $2\cdot(11 + 2\cdot(4 - 9))/2 = 1$、3 → $3\cdot(11 + 3\cdot(6 - 9))/2 = 3$)。逆向きの要素では入口と出口が入れ替わるので、1 (入口) と 2 (出口) を入れ替える。`tcav` の使い方: 入口の端は `mfring ≥ 0` かつ `mfring ≠ 2`、出口の端は `mfring ≥ 0` かつ `mfring ≠ 1` (`tcav.f:82`、`185`)。したがって `mfring` = 0 (FRINGE = 0) と 3 は両方、1 は入口だけ、2 は出口だけ。
- `dvfs`: RFSW での周波数のずらしによる $dv$ への足し込み ([01 章](01_variables_and_conventions.md) §8 の表)。`al`: 空洞の長さ $L$、`alx`: 刻みの中のドリフトの長さ (§6.1)。

SAD の座標は $(x, p_x, y, p_y, z, \delta)$ ([01 章](01_variables_and_conventions.md))。$z = -v(t - t_0(s))$、先頭が $z > 0$。追跡の内部では $p_x, p_y$ を $q = p/(1+\delta)$ の単位で持つ ([04 章](04_quad.md) 式 (4.37))。

## 3. エネルギー利得と時間

### 3.1 位置から時間へ

SAD の $z$ は $z = -v\,(t - t_0(s))$ ($v$ は粒子の速度、$t_0(s)$ は設計粒子の到着時間、[01 章](01_variables_and_conventions.md) 式 (1.7))。$c = 1$ の単位で $v = \beta = p/h$ ($p, h$ は粒子の運動量とエネルギーを $mc$、$mc^2$ で測ったもの) だから、設計粒子から見た粒子の時間のずれは

$$
t = -\frac{z}{\beta} = -z\,\frac{h_1}{p_1}\qquad(h_1 = \sqrt{1 + p_1^2},\ p_1 = p_0(1 + \delta))
\tag{10.1}
$$

(`t=-z(i)*h1/p1`、`tcav.f`)。$h_1$ は、$dv$ が小さいときは $dv$ から作る: $v/v_0 = p\,h_0/h_1$ ([01 章](01_variables_and_conventions.md) 式 (1.9)) から作る。導出 (1 行ずつ): $v/v_0 = (p/h_1)/(p_0/h_0) = p_{1r}\,h_0/h_1$ ($p = p_0p_{1r}$、$p_{1r} = 1 + \delta$)。一方、$dv$ の定義 $dv = 1 - v/v_0$ に周波数のずらし `dvfs` (RFSW での $dv$ への足し込み。`tcav.f` が `dvfs` を $dv$ に足す) を含めた式は $v/v_0 = 1 - dv + \mathtt{dvfs}$ で、両辺を等しくおく: $p_{1r}\,h_0/h_1 = 1 - dv + \mathtt{dvfs}$。これを $h_1$ について解くと $h_1 = p_{1r}\,h_0/(1 - dv + \mathtt{dvfs})$ (`dv` は粒子の持つ値。$\delta$ の小さい領域で $h_1 = \sqrt{1 + p_1^2}$ を直接引くより桁落ちが少ない)。$dv > 0.1$ では $h_1 = \sqrt{1 + p_1^2}$ (`p2h`)。RFSW OFF では `tlim=0` で、$t = \min(t_{\lim}, \max(-t_{\lim}, t)) = 0$ に潰れる (`tcav.f` の `tlim`)。つまり RF の位相に依存しない (§5.3)。

### 3.2 エネルギー利得

ヘルプの式を、1 つの粒子について $V(x, y) \equiv V + V_1x + \tfrac12V_{20}x^2 + \tfrac12V_{02}y^2 + V_{11}xy$ と書くと $\Delta E = -eV(x, y)\sin(\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$。$t_s$ は「放射の損失と加速のつりあいで決まる同期時間」(ヘルプ)。$\mathtt{phis} \equiv \omega t_s$ (運転の仕方で 3 通り。`tcav.f:35` の `phis = merge(0, merge(phic, w*trf0, autophi), trpt)`: TRPT のとき $\mathtt{phis} = 0$、autophi のとき $\mathtt{phis} = \mathtt{phic}$、それ以外は $\mathtt{phis} = \texttt{w}\cdot\texttt{trf0}$)、$\mathtt{phic} \equiv \mathrm{PHI} + \mathrm{DPHI}$ (電荷の符号つき) とおくと $\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI} = \omega t - \mathtt{dphis}$ ($\mathtt{dphis} = \mathtt{phis} - \mathtt{phic}$)。$h = E/mc^2$、$v = V/mc^2$ ($e|\mathtt{charge}|$ を含む) で、

$$
\Delta h = -v(x, y)\,\sin(\omega t - \mathtt{dphis})\ +\ \text{(一定の項。§3.3)}
\tag{10.2}
$$

K64 は、(10.2) を $\omega t$ が小さいときに桁落ちしない形で書く。$ph \equiv \omega t/2$ とおき、$2\sin A\cos B = \sin(A + B) + \sin(A - B)$ で ($A = ph$、$B = ph - \mathtt{dphis}$)

$$
-2\sin(ph)\cos(ph - \mathtt{dphis}) = -\sin(2\,ph - \mathtt{dphis}) - \sin(\mathtt{dphis}) = -\sin(\omega t - \mathtt{dphis}) - \sin\mathtt{dphis}
\tag{10.3}
$$

なので、`dh=(v+…)*(-2*sin(ph)*cos(ph-dphis)+offset)*wsn` は

$$
\Delta h = -v(x, y)\,\bigl[\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{dphis} - \mathtt{offset}\bigr]\cdot\mathtt{wsn}
\tag{10.4}
$$

(`wsn` はこの刻みに割り当てる割合、§6。薄い空洞は 1。) $t \to 0$ で第 1 項と第 2 項が消し合い、$\Delta h \to v\cdot\mathtt{offset}$ だけが残る形になる。

### 3.3 `offset` と同期位相

`offset` は、運転の仕方で決まる (`tcav.f`):

- **放射・TRPT・autophi** ($\mathtt{rad}$、$\mathtt{trpt}$、$\mathtt{autophi}$): $\mathtt{offset} = \sin\mathtt{dphis}$、$\mathtt{offset1} = 0$。(10.4) は $-v\sin(\omega t - \mathtt{dphis})$ で、ヘルプの $\Delta E$ そのもの。このとき $\mathtt{dphis} = \mathtt{phis} - \mathtt{phic}$ は、TRPT で $-\mathtt{phic}$、autophi で $0$ ($\mathtt{offset} = 0$)、放射だけのとき $\omega\,\mathtt{trf0} - \mathtt{phic}$ と、モードで中身が違う。基準の粒子 ($t = 0$) の利得は $v\sin\mathtt{dphis}$ で、放射の損失とつりあう (つりあいの位相 $t_s$ を `trf0` が持つ)。
- **それ以外 (RAD OFF の追跡)**: $\mathtt{offset} = \sin\mathtt{dphis} - \sin\mathtt{phis}$、$\mathtt{offset1} = \sin\mathtt{phis}$。(10.4) は $\Delta h = -v\,[\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{phis}]$。すなわち**基準の粒子の利得を引いた**形で、`trf0` = 0 ($\mathtt{phis} = 0$、FrequencyMap・LMA の既定) なら、ヘルプの $\Delta E$ と同じ。$\mathtt{phis} \neq 0$ では、$-v\sin\mathtt{phis}$ だけ一定の利得が加わる (放射がないので、基準の粒子の損失を考えず、つりあいの位相にすると基準の粒子が加速も減速もしない形に調整する。ソースからの読み。未確認)。

数値確認 ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (a)): 3 GeV の電子、VOLT = 2 MV、FREQ = 500 MHz、PHI = 0.3 rad で、$\mathtt{trf0} = 0$ と $0.4/\omega$ のとき、(10.1) と (10.4) から作った写像が、追跡と $2 \times 10^{-13}$ 以内で一致する (例: $\mathrm{trf0} = 0$、$z = -0.1$ m、$\delta = +0.002$ で追跡は $\delta$ を $+0.001349822$ にする。変化 $-6.5 \times 10^{-4}$ は、$t = -z/\beta = 0.1$ m、$\omega t = 1.048$、$\mathtt{dphis} = -0.3$ から $\Delta E/(p_0c) = -(2 \times 10^6/3 \times 10^9)\sin(1.048 + 0.3) = -6.5 \times 10^{-4}$ に合う。出力 (a)。$\omega = 2\pi\cdot500\,\mathrm{MHz}/c = 10.479$ /m)。

## 4. 薄い空洞の写像

(10.4) でエネルギーが $h_1 \to h_2 = h_1 + \Delta h$ に変わると、粒子の運動量・$z$・$dv$・横運動量を更新する。`tcav.f` の `do 20` (と `tcavfrin`、`ttcav.f`) の手順を式にする。

1. **運動量**: $h$ と $p$ の関係 $h^2 = 1 + p^2$ ($mc$ 単位) から $p_2^2 - p_1^2 = h_2^2 - h_1^2 = \Delta h\,(h_1 + h_2)$。$a \equiv \Delta h\,(h_1 + h_2)/p_1^2$ として

$$
\frac{p_2}{p_1} = \sqrt{1 + a},\qquad \mathtt{dpr} = \sqrt{1 + a} - 1\ (\texttt{sqrt1}(a)\text{、桁落ちなし}),\qquad p_{2r} = p_{1r}(1 + \mathtt{dpr}),\ \ \delta_2 = p_{2r} - 1
\tag{10.5}
$$

($\Delta h$ が小さいと $a$ も小さく、`sqrt1` が $a/(1 + \sqrt{1 + a})$ で桁落ちを避ける。$\sqrt{1+a} - 1 = [(1+a) - 1]/(\sqrt{1+a} + 1) = a/(1 + \sqrt{1+a})$。`max(a, -1)` は $p_2^2 \ge 0$ の下限。)

2. **$dv$**: 新しい $dv = 1 - v_2/v_0$ は、01 章の式 (1.9) の形で $dv = -(1 + p_{2r})\,\delta_2/[h_2(h_2 + p_{2r}h_0)] + \mathtt{dvfs}$ (`dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs`)。この式の導出: 速度比は $v/v_0 = (p_{2r}p_0/h_2)/(p_0/h_0) = p_{2r}h_0/h_2$ なので $dv = 1 - p_{2r}h_0/h_2 = (h_2 - p_{2r}h_0)/h_2$。分子は $h_2^2 - p_{2r}^2h_0^2 = (1 + p_{2r}^2p_0^2) - p_{2r}^2(1 + p_0^2) = 1 - p_{2r}^2 = -\delta_2(1 + p_{2r})$ ($h_2^2 = 1 + p_2^2 = 1 + p_{2r}^2p_0^2$、$h_0^2 = 1 + p_0^2$) なので、$h_2 - p_{2r}h_0 = -\delta_2(1 + p_{2r})/(h_2 + p_{2r}h_0)$。したがって $dv = -\delta_2(1 + p_{2r})/[h_2(h_2 + p_{2r}h_0)]$ ($\delta_2 = dp2r$)。01 章 (1.9) と同じ形で、分母の第 2 項は $p_{2r}h_0$ ($p_2 = p_{2r}p_0$ ではない。13 章 §4.1 の $dv$ との違い)。

3. **$z$**: 時間 $t$ は、空洞が薄い (粒子が動かない) 間は変わらない。新しい速度 $\beta_2 = p_2/h_2$ のもとで $z = -\beta t$ だから

$$
z_2 = -t\,\beta_2 = -t\,\frac{p_{2r}\,p_0}{h_2}
\tag{10.6}
$$

($p_2 = p_{2r}p_0$。`z(i)=-t*p2r*p0/h2-dzn`。`dzn` は設計の速度の変化の分、§6。)

4. **横運動量**: 縦方向の加速は、正準の $p_x, p_y$ を変えない ($H$ に $x, y$ が出ない薄い空洞の流れ。ただし V1 などがあるときは §5)。追跡の内部変数 $q = p_x/p$ は、$q_2 = q_1\,p_{1r}/p_{2r}$ (`px(i)=px(i)*p1r/p2r`)。つまり**加速で $q$ が減る (断熱減衰)**。

この写像は、(i) 運動量・エネルギーの関係が厳密、(ii) $z$ の更新が時間保存、(iii) 正準運動量が不変、なので、1 つの粒子についての写像として厳密である (空洞の長さが 0 の極限)。

**シンプレクティック性の理由** (端なし、V1 などなし)。座標を $(t, E)$ に取り替えて、行列を 1 段ずつ書く。以下 $E = h$ ($mc^2$ 単位)、$S = \begin{pmatrix}0&1\\-1&0\end{pmatrix}$。

(i) **エネルギーの足し算のヤコビ行列**。この写像は「時刻 $t$ を保ち、$E \to E + f(t)$ と進める」もの。$f(t) \equiv -v\sin(\omega t - \mathtt{dphis})$ は $t$ だけの関数 ((10.2))。新しい値は $t' = t$、$E' = E + f(t)$ なので

$$
J_{tE} = \begin{pmatrix}\partial t'/\partial t & \partial t'/\partial E\\ \partial E'/\partial t & \partial E'/\partial E\end{pmatrix} = \begin{pmatrix}1 & 0\\ f'(t) & 1\end{pmatrix},\qquad \det J_{tE} = 1
\tag{10.4a}
$$

2×2 の行列は $J^\top SJ = (\det J)\,S$ を満たす (成分で書く: $SJ = \begin{pmatrix}f' & 1\\ -1 & 0\end{pmatrix}$、$J^\top = \begin{pmatrix}1 & f'\\ 0 & 1\end{pmatrix}$ なので $J^\top SJ = \begin{pmatrix}f' - f' & 1\\ -1 & 0\end{pmatrix} = S$)。したがって $(t, E)$ の空間では、この写像は正準である。これは $\Delta E = \partial_tK$ の形のキック ($K(t) = \int^t f\,dt'$) の言い換えでもある: 正準の規則 $\Delta P = -\partial_tK$ ($P \equiv -E$) と同じ。

(ii) **$(t, E) \to (z, \delta)$ の座標変換のヤコビ行列**。$z = -\beta(E)\,t$、$\delta = p(E)/p_0 - 1$ で、$p = \sqrt{E^2 - 1}$、$\beta = p/E$。微分は $dp/dE = E/p = 1/\beta$、$d\beta/dE = [E\cdot E/p - p]/E^2 = 1/(E^2p)$ (商の微分: $\beta = p/E$、$\beta' = (p'E - p)/E^2$、$p'E - p = (E^2 - p^2)/p = 1/p$)。したがって

$$
M \equiv \frac{\partial(z, \delta)}{\partial(t, E)} = \begin{pmatrix}-\beta & -t/(E^2p)\\ 0 & 1/(\beta p_0)\end{pmatrix},\qquad \det M = -\beta\cdot\frac{1}{\beta p_0} = -\frac1{p_0}
\tag{10.4b}
$$

(上三角なので、行列式は対角の積。**定数**: $E$ にも $t$ にもよらない。) 括弧で書くと $[z, \delta]_{(t,P)} = \partial_tz\,\partial_P\delta - \partial_Pz\,\partial_t\delta = (-\beta)\cdot(-1/(\beta p_0)) - 0 = +1/p_0$ ($P = -E$ なので $\partial_P = -\partial_E$)。$x, y$ の対も、$p_x$ を $p_0$ で規格化した $p_x/p_0$ では $[x, p_x/p_0] = 1/p_0$。つまり SAD の 3 つの対 $(x, p_x/p_0)$、$(y, p_y/p_0)$、$(z, \delta)$ は、どれも標準の正準対の**共通の定数 $1/p_0$ 倍**で、座標変換は「正準変換 + 全体の定数倍」にあたる (付録 C の括弧の規約。この章の $(t, P)$ の対のキックは $\exp(-:K:)$、付録 C の記法)。

(iii) **$(z, \delta)$ の空間での写像**。$(z, \delta) = G(t, E)$ とおくと、写像は $G\circ F\circ G^{-1}$ ($F$: (i) の写像)。そのヤコビ行列は $J_{z\delta} = M(F(q))\,J_{tE}\,M(q)^{-1}$ (連鎖律。$q = (t, E)$)。行列式は

$$
\det J_{z\delta} = \det M\cdot\det J_{tE}\cdot\frac{1}{\det M} = 1
$$

($\det M$ が定数なので、入口と出口で打ち消し合う。)  2×2 で $\det = 1$ なら $J_{z\delta}^\top SJ_{z\delta} = S$ ((i) と同じ式)。$x, y, p_x, p_y$ は変わらない (単位行列) ので、6×6 のヤコビ行列は 2×2 のブロックの直和で、$J^\top SJ = S$ が全体で成り立つ。$z, \delta$ の定義 ([01 章](01_variables_and_conventions.md)) のうち、$dv$ は $E$ だけの関数 ((10.6) の $z_2 = -t\beta_2$) なので上の $z = -\beta(E)t$ に含まれる。

数値確認 ([cav_03_symplectic_proofs.out.txt](outputs/cav_03_symplectic_proofs.out.txt) (a)): sympy で $\det M = -1/p_0$、$\det J_{z\delta} = 1$、数値の $J_{z\delta}^\top SJ_{z\delta} - S$ の最大成分 $1 \times 10^{-16}$。

V1・V20 などで $\Delta E$ が $(x, y)$ に依存するときは、§5.1 の横のキックが相手を務める ($K(x, y, t)$ から $\Delta E = \partial_tK$、$\Delta p_x = -\partial_xK$ を同時に取るので正準。$x$-$t$ 成分の混合 $\partial_x\partial_tK$ が両方に入ることが要点。§5.1 に行列を書く)。

**確認**: 6×6 のヤコビアン (中心差分、刻み $10^{-7}$) で、V1・V20 なしの薄い空洞 $|J^\top SJ - S| = 1.2 \times 10^{-10}$ ($\det J - 1 = 1.2 \times 10^{-10}$)、V1・V20・V11・V02 ありで $2.1 \times 10^{-10}$ ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (b)。中心差分の精度)。**これは端 ($H_f$、§7) なしの空洞についての確認である。端つきの空洞の同じ量は (b) の 3 行目 ($L = 0.5$、フリンジあり) で $1.5 \times 10^{-10}$ (v2。k64 では $2.3 \times 10^{-10}$) で、フリンジなしの 4 行目 ($1.6 \times 10^{-10}$) と区別できない。この例 ($x = y = 1$ mm、3 GeV/$c$) の入口の端だけの欠陥の見積り (§7.2) は $|\mathtt{dpr}| \simeq (h/p^2)\,\omega p_0v_f\,(x^2 + y^2)\,|\cos\theta|/2 = 7 \times 10^{-9}\,|\cos\theta|$ ($\theta = \omega t - \mathtt{dphis}$) で、出口の端が $v \to -v$ で逆符号の寄与を持つこと (対称な配置で打ち消すと読めるが、その数値確認は未実施) と、中心差分の精度 ($10^{-10}$) から、この例では欠陥の有無が判定できない。端の欠陥の大きさは §7.2 の再実装で測る。**

## 5. 横方向の依存 (V1、V20、V11、V02) と非相対論の補正

### 5.1 Panofsky–Wenzel

エネルギー利得が $(x, y)$ に依存するとき、横のキックが同時に要る (Panofsky–Wenzel の定理)。SAD の座標 ($z$ と $\delta$ が共役、$z = -\beta t$) では、時間依存のキックのハミルトニアン $K(x, y, t)$ について、$\Delta\delta = -\partial K/\partial z = (1/\beta)\,\partial K/\partial t$ なので $\Delta E = \beta\,\Delta p = \partial K/\partial t$。(10.2) の $-eV(x, y)\sin(\omega t - \mathtt{dphis})$ を作るのは

$$
K(x, y, t) = \frac{e\,V(x, y)}{\omega}\,\cos(\omega t - \mathtt{dphis})
\tag{10.7}
$$

($\partial_tK = -eV\sin(\ldots)$)。横のキックは $\Delta p_x = -\partial K/\partial x$、$\Delta p_y = -\partial K/\partial y$ だから、$p_0$ で規格化して

$$
\Delta p_x = -\frac{\partial_xV}{\omega\,p_0}\cos(\omega t - \mathtt{dphis}),\qquad
\Delta p_y = -\frac{\partial_yV}{\omega\,p_0}\cos(\omega t - \mathtt{dphis})
\tag{10.8}
$$

(単位は $v = V/mc^2$ で測る)。$\partial_xV = V_1 + V_{20}x + V_{11}y$ ($=\mathtt{v1a} + 2\,\mathtt{v20a}\,x + \mathtt{v11a}\,y$、`v20a` $= \tfrac12V_{20}\ldots$)、$\partial_yV = V_{11}x + V_{02}y$。コードは `av=-(cos(2*ph-dphis)*wi-offset1*t)*wsn/p0`、`dpx=(v1a+2*v20a*x+v11a*y)*av`、`dpy=(v11a*x+2*v02a*y)*av`、`px=(px*p1r+dpx)/p2r` ($\mathtt{wi} = 1/\omega$)。$\cos(2\,ph - \mathtt{dphis}) = \cos(\omega t - \mathtt{dphis})$ で (10.8) と一致する。さらに $\mathtt{offset1}\ne0$ の項は、(10.4) の一定の利得 $-v(x, y)\,\mathtt{offset1}$ に対応するキックのハミルトニアン $K_{\rm extra} = -V(x, y)\,\mathtt{offset1}\,t$ ($\partial_tK_{\rm extra} = -V\,\mathtt{offset1}$) から $\Delta p_x = +\partial_xV\,\mathtt{offset1}\,t/p_0$ ($= -\partial_xK_{\rm extra}$)、コードの $-(\ldots - \mathtt{offset1}\,t)$ の第 2 項と同じ。つまり、$\mathtt{v1a}$ などのある空洞も、(10.7) と一定の項から導かれる**1 つの時間依存のハミルトニアンの流れ**に当たり、ヤコビアンの数値確認 (§4) と整合する。**横のキックの符号と正準性** (式を省かない形)。(10.7) の $K$ から、$(x, p_x, t, P)$ ($P = -E$) の写像は

$$
x' = x,\quad p_x' = p_x - s\,\partial_xK,\quad t' = t,\quad P' = P - \partial_tK\qquad(s = +1\ \text{が (10.8)、}s = -1\ \text{は符号を逆にした場合})
$$

(簡単のため $y$ は省く。$P' = P - \partial_tK$ は $\Delta E = +\partial_tK$ と同じ。) 偏導関数を $K_{xx} = \partial_x^2K$、$K_{xt} = \partial_x\partial_tK$、$K_{tt} = \partial_t^2K$ と書くと、ヤコビ行列と括弧の行列 $\Omega$ は

$$
J = \begin{pmatrix}1&0&0&0\\ -sK_{xx}&1&-sK_{xt}&0\\ 0&0&1&0\\ -K_{xt}&0&-K_{tt}&1\end{pmatrix},\qquad
\Omega = \begin{pmatrix}0&1&0&0\\ -1&0&0&0\\ 0&0&0&1\\ 0&0&-1&0\end{pmatrix}
$$

$D \equiv J^\top\Omega J - \Omega$ を成分で計算すると ($\Omega J$ を作り、$J^\top$ をかける。成分は sympy でも確認)、0 でない成分は $x$-$t$ の対だけで

$$
D_{xt} = -(s - 1)\,K_{xt},\qquad D_{tx} = +(s - 1)\,K_{xt}
$$

$s = +1$ なら $D = 0$ (勾配キックは正準)。$s = -1$ だと $D_{xt} = 2K_{xt} = -2\,\partial_xV\,\sin(\omega t - \mathtt{dphis})$ ($K_{xt} = \partial_x[\partial_tK] = -\partial_xV\sin(\ldots)$)。つまり**符号を逆にした欠陥は $2\,\partial_xV\sin(\omega t - \mathtt{dphis})$、$\partial_xV \sim V_1$ の大きさ**で。正しい符号 ($s = +1$) では、$J$ の $(p_x, t)$ 成分 $-K_{xt}$ と $(P, x)$ 成分 $-K_{xt}$ が、偏微分の順序の交換 $\partial_x\partial_t = \partial_t\partial_x$ で同じ値になり、$D_{xt}$ で打ち消し合う ($K_{xx}$・$K_{tt}$ は $D$ に入らない)。$(z, \delta)$、$p_x/p_0$ で測った値は、$D_{\rm SAD} = (1/p_0)\,A^{-\top}DA^{-1}$ ($A = \partial(x, p_x/p_0, z, \delta)/\partial(x, p_x, t, P)$。(10.4b) の $M$ を使う。$A^\top S_{\rm SAD}A = \Omega/p_0$ からの帰結) となる。

**数値** ([cav_03_symplectic_proofs.out.txt](outputs/cav_03_symplectic_proofs.out.txt) (b)(c)): 薄い空洞の写像 (10.1)(10.2)(10.5)(10.6)(10.8) を numpy で組み (オフセットなし。V1・V20・V02・V11 あり。$x = 2$ mm、$y = 1$ mm、$z = 0.05$ m)、中心差分のヤコビアンの $|J^\top SJ - S|$ は

| 運動量 | $s = +1$ (正しい符号) | $s = -1$ (符号を逆) | 解析の予想 ($s = -1$) |
|---|---|---|---|
| 3 GeV/$c$ | $6.5 \times 10^{-10}$ | $3.70 \times 10^{-4}$ | $3.70 \times 10^{-4}$ |
| 2 MeV/$c$ | $6.0 \times 10^{-10}$ | $0.615$ | $0.615$ |

(予想との差は $6 \times 10^{-10}$ 以内。条件: VOLT = 2 MV、500 MHz、dphis = $-0.3$、V1 = $2 \times 10^6$ V/m、V20 = $2 \times 10^8$、V02 = $-10^8$、V11 = $10^8$ (V/m²)。)  したがって (§4) の数値確認 (V1・V20・V11・V02 あり: $2 \times 10^{-10}$) は、横のキックの符号と、エネルギー利得との整合の確認にもなっている。ただし、この表の $s = \pm1$ は、いずれも式 (10.1)(10.2)(10.5)(10.6)(10.8) の numpy の再実装の値で、pysad の追跡そのもので符号を反転した試験は行っていない (pysad の追跡で符号が正しいことは、上の §4 の (b) の $2 \times 10^{-10}$ が示す)。

### 5.2 非相対論の補正

ヘルプ (VOLT の説明): 「(CAVI only) $(\mathrm{VOLT}+\mathrm{DVOLT})\,(2\pi\,\mathrm{FREQ}/c)^2/(\gamma\beta)^2/4$ が自動で V20・V02 に足される。ローレンツ因子は入口と出口の $1/(\beta\gamma)$ の平均の逆数」(原文で確認)。コードは `vcorr=v*(w*(.5d0/p0+.5d0/pe))**2/4`、`v20a=v20a+vcorr`、`v02a=v02a+vcorr` (`tcav.f`)。ここで $\mathtt{v20a} = \tfrac12V_{20}/\ldots$ が $x^2$ の係数 ($V_{20}x^2/2$ の $1/2$ を含む) なので、足される量は $\Delta V(x, y) = V\,(\omega/\beta\gamma c)^2\,(x^2 + y^2)/4$。

これは、$\omega$ で振動する軸対称の加速モード ($\mathrm{TM}_{010}$ 型) の $E_z$ の半径依存である。$\nabla^2E_z - (1/c^2)\partial_t^2E_z = 0$ に、粒子が速度 $\beta c$ で通るときの依存 $E_z \propto e^{i\omega(t - s/\beta c)}$ ($\partial_s^2 \to -\omega^2/\beta^2c^2$、$\partial_t^2 \to -\omega^2$) を代入すると、横方向の方程式は

$$
\frac1r\partial_r(r\,\partial_rE_z) - \Bigl(\frac{\omega^2}{\beta^2c^2} - \frac{\omega^2}{c^2}\Bigr)E_z = 0\ \ \Longleftrightarrow\ \ \nabla_\perp^2E_z = \Bigl(\frac{\omega}{\beta\gamma c}\Bigr)^2E_z
\tag{10.9}
$$

($1/\beta^2 - 1 = 1/(\beta\gamma)^2$ を使った。) 解は変形ベッセル関数 $E_z \propto I_0(\omega r/\beta\gamma c) = 1 + (\omega r/\beta\gamma c)^2/4 + O(r^4)$ ($I_0(u) = \sum(u/2)^{2k}/(k!)^2$)。第 2 項が `vcorr`。数値確認 ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (c)): 2 MeV/c の電子 ($\beta\gamma = 3.914$)、VOLT = 1 MV、FREQ = 500 MHz、PHI = $\pi/2$、$z = 0$ で、ΔE($x$)/ΔE(0) が $x = 0.01, 0.03, 0.05$ m で 1.000179、1.001613、1.004480、$1 + (\omega x/\beta\gamma c)^2/4$ と同じ ($I_0$ とは第 3 項 $u^4/64$ だけ違い、$x = 0.05$ で 1.004485)。

### 5.3 RFSW OFF

`rfsw` が偽 (RFSW OFF) の追跡では、$t_{\lim} = 0$ で $t = 0$ に潰れ、位相依存の項 $\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{dphis}$ が 0 になる (10.4) から、$\Delta h = v\cdot\mathtt{offset}$ (一定) だけが残る。RAD OFF では $\mathtt{offset} = \sin\mathtt{dphis} - \sin\mathtt{phis}$、$\mathtt{phis} = 0$ なら $\sin\mathtt{dphis} = -\sin\mathtt{phic}$ で、$\mathrm{PHI} \ne 0$ のときだけ一定のエネルギー利得がある (確認していない)。CALC の TCAVI (`ttcave.f:15-18`) は、RFSW OFF でドリフト (`tdrife0`) にする。CALC の CAVI (`tcave.f`) は `if(rfsw)then` の中で加速を入れる (RFSW OFF の分岐は読んでいない)。ヘルプの RFSW の説明: 「CAVI と TCAVI による加速を、粒子の追跡と emittance 計算で有効にする。FFS は起動時に必ず RFSW を ON にする」(原文で確認)。

## 6. 長い空洞の分割、$p_0(s)$ の変化

### 6.1 ドリフト–キック–ドリフト

$L \ne 0$ のとき、`tcav` は、入口の端 (§7) → 半分のドリフト → [キック → ドリフト] … → キック → 半分のドリフト → 出口の端、を作る。刻み $n$ のキックは強さの割合 `ws(n)`、ドリフトの長さは $L\cdot(\mathtt{ws}(n) + \mathtt{ws}(n-1))/2$ (最初と最後は $L\,\mathtt{ws}/2$)。各ドリフトは運動学的なドリフト ([02 章](02_drift.md))、$z$ の更新に `-dv*alx` を持つ。キックは §3〜§5。重みの和が 1 で全体の電圧に一致し ($\sum\mathtt{ws} = 1$)、各 $\mathtt{ws}(n)$ は隣り合う 2 つのドリフトに半分ずつ入るので、長さの合計は $L\sum\mathtt{ws} = L$ になる (2 次の対称分割)。

**刻みの数**: `tcav.f`: $\mathrm{ndiv} = 1 + \mathrm{int}\bigl(\min(|\omega L|,\ \sqrt{r^2/(3\varepsilon)})\bigr)$、$\varepsilon = 10^{-3}$、$r = |\mathtt{vnominal}\,(1/h_0 + 1/h_e)|$。$r = 0$ (TRPT でない普通のリング) のとき $\mathrm{ndiv} = 1$、キックは 1 回 (`ws = 1`)。$r > 10^{-6}$ のとき、$\mathrm{ndiv} = \min(1000, \max(\mathrm{ndiv}, \mathrm{nint}(10\,wl)))$。CALC (`tcave.f`) は $\varepsilon = 10^{-2}$ で、$\mathrm{ndiv} = 1 + \mathrm{int}\bigl(\min(|\omega L|,\ \sqrt{[v(1/h_0 + 1/h_1)]^2/(3\varepsilon)})\bigr)$ ($v$ は電圧そのもの。TRPT の `vnominal` ではない)。式の読み: $r = |\mathtt{vnominal}(1/h_0 + 1/h_e)|$ は、空洞の中の相対運動量変化の大きさの目安 ($\Delta p/p \simeq \Delta h\,h/p^2$ で、$h \gg 1$ なら $\Delta h/h$ の程度。$1/h_0 + 1/h_e \simeq 2/h$ なので $r \simeq 2\Delta h/h$)。$\mathrm{ndiv} \ge r/\sqrt{3\varepsilon}$ は、1 刻みあたりの相対変化 $r/\mathrm{ndiv}$ を $\sqrt{3\varepsilon} = 0.055$ ($\varepsilon = 10^{-3}$) 以下に抑える条件と読める。係数の 3 と $\varepsilon$ の値の根拠 (どの誤差を何に抑えるか) は**ソースに書かれていない**。

### 6.2 設計運動量が増える (TRPT) とき

`vnominal` $= v\sin(-\mathrm{PHI}\cdot\mathrm{sign}(\mathtt{charge}))$ ($\mathtt{trpt}$ のときだけ。TRPT は「ビームラインが蓄積リングではなく輸送路で、運動量が加速で変わる」ことを宣言、ヘルプの TRPT。DPHI は設計運動量 $p_0(s)$ に入れない、ヘルプの DPHI) は、空洞を設計粒子が通ったときのエネルギー利得 $h_e - h_0 = \mathtt{vnominal}$。空洞の中のエネルギーは一様に増えるとし、各刻みの利得を `vnominal*wsn` とする (`sv=sv+vnominal*wsn`、`h2=h0+sv`)。この間の設計の運動量・速度の変化を、$dv_n$ と $dz_n$ として $z$ に反映する (`dzn=-dvn*al*wsn*0.5`)、空洞を出たあと $\delta$ を新しい設計運動量 $p_e$ に対して測り直す: $\delta' = \delta\,p_0/p_e - \mathtt{dpepe}$ ($\mathtt{dpepe} = \mathtt{vnominal}\,(h_e + h_0)/((p_e + p_0)\,p_e)$。`g(i)=dp1`)、$dv$ も作り直し、`pgev=pe*amass`、`tphyzp` で設計運動量の物理量を更新する (`tcav.f` の最後)。

**空洞のあとの $\delta$ の測り直し**: 設計粒子のエネルギーは $h_e = h_0 + \mathtt{vnominal}$ になる。$p_e^2 - p_0^2 = h_e^2 - h_0^2 = \mathtt{vnominal}\,(h_e + h_0)$ ($h^2 = 1 + p^2$ の差) なので、$p_e - p_0 = \mathtt{vnominal}(h_e + h_0)/(p_e + p_0)$、したがって $\mathtt{dpepe} \equiv (p_e - p_0)/p_e = \mathtt{vnominal}(h_e + h_0)/[(p_e + p_0)p_e]$ (`tcav.f` の `dpepe` の式そのもの)。粒子の運動量は $p = p_0(1 + \delta)$、新しい設計運動量 $p_e$ で測った偏差 $\delta' = p/p_e - 1 = (p_0/p_e)(1 + \delta) - 1 = \delta\,p_0/p_e + (p_0 - p_e)/p_e = \delta\,p_0/p_e - \mathtt{dpepe}$ (`dp1=dp*p0/pe-dpepe`)。設計粒子 ($\delta = (p_e - p_0)/p_0$ で出る) は $\delta' = (p_e - p_0)/p_e - \mathtt{dpepe} = 0$ になる。`dv` は (1.9) の形 ($h_1$ は `p2h(pe*pr1)`、分母は $h_1 + p_{r1}h_e$) で作り直す。

**設計粒子の $z$ のずれ** (`dvn`、`dzn`): 各刻みで、ドリフトの `z -= dv*alx` の $dv$ は、空洞の入口の設計速度 $v_0$ に対する値 ($p_{2r}$ が入口の $p_0$ 規格)。しかし設計粒子自身の速度は、空洞の中で $v_0 \to v_n$ と変わる。`dvn` は、その刻みの設計粒子の $dv$ ($h_2 = h_0 + \mathtt{sv}$ の運動量 $p_2$、`dp2` $= \mathtt{sv}(h_2 + h_0)/((p_2 + p_0)p_0)$ で、上と同じ代数 ($p_2^2 - p_0^2 = \mathtt{sv}(h_2 + h_0)$) の $p_2/p_0 - 1$。(1.9) の形)。`dzn` $= -dvn\cdot\mathtt{al}\cdot\mathtt{wsn}/2$ を $z$ から引く操作は、各刻みの前後 (2 回) で、設計粒子に対してドリフトが余分に引く $dvn\cdot\mathtt{al}\cdot\mathtt{wsn}$ を打ち消して、設計粒子の $z$ が 0 のままになるようにする補正と読める。**この補正が設計粒子で $z = 0$ を保つことの数値確認は行っていない** (TRPT の追跡は pysad の通常の経路 (`track_particles`) では有効にならず、確認できなかった)。

**重み ws**: $r > 10^{-6}$ のとき

$$
wl = \mathrm{arcsinh}\,\frac{r(2 + r)}{2(1 + r)},\quad r_1 = \frac{wl}{2\,\mathrm{ndiv}},\quad \mathtt{ws}(1) = \frac{2e^{r_1}\sinh r_1}{r} = \frac{e^{2r_1} - 1}{r},\quad \mathtt{ws}(n) = \mathtt{ws}(n-1)\,e^{2r_1}
\tag{10.10}
$$

(`we=1+ws1`、$e^{2r_1} = 1 + 2e^{r_1}\sinh r_1$。) $\mathtt{ws}$ は公比 $e^{2r_1}$ の等比数列で、和は $\mathtt{ws}(1)\,(e^{2r_1\mathrm{ndiv}} - 1)/(e^{2r_1} - 1) = (e^{wl} - 1)/r$。さらに $x \equiv r(2 + r)/(2(1 + r))$ として $x^2 + 1 = (r^2 + 2r + 2)^2/(4(1 + r)^2)$ ($(r^2 + 2r)^2 + 4(1 + r)^2 = r^4 + 4r^3 + 8r^2 + 8r + 4 = (r^2 + 2r + 2)^2$) なので、$e^{wl} = x + \sqrt{x^2 + 1} = [(r^2 + 2r) + (r^2 + 2r + 2)]/(2(1 + r)) = 1 + r$。したがって**$\sum\mathtt{ws} = r/r = 1$、最後の刻みの重みは最初の $(1 + r)^{(\mathrm{ndiv}-1)/\mathrm{ndiv}}$ 倍** ($r = 0.5$、$\mathrm{ndiv} = 12$ で 1.4502、$r = 3$、$\mathrm{ndiv} = 20$ で 3.7321。[cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (d))。設計粒子のエネルギーが増えるにつれ、長さと電圧の割合を大きくする形で、運動量が大きいほどキックの長さが長い。**$wl$ の選び方**: 式 (10.10) の $wl$ は、$e^{wl} = 1 + r$、すなわち $wl = \ln(1 + r)$ と同じ (上の代数。数値: 出力 (d)、および [cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (e) で $\mathrm{arcsinh}\,x$ と $\ln(1 + r)$ が 15 桁一致)。$wl = \ln(1 + r)$ が**和が 1 になる条件から強制される**ことは導出できた (等比数列の和 $\sum\mathtt{ws} = (e^{wl} - 1)/r$ が 1 になるには $e^{wl} = 1 + r$)。**等比の重み (刻みが $(1 + r)^{1/\mathrm{ndiv}}$ ずつ長くなる) を選ぶ理由 (運動量が大きい側で粗くする意図と読めるが) は、ソースに書かれていない** (未確認)。

### 6.3 V20 の補正と平均

`vcorr` の運動量は、入口 $p_0$ と出口 $p_e$ の $1/p$ の平均で入れる ($0.5/p_0 + 0.5/p_e$、ヘルプ「入口と出口の $1/(\beta\gamma)$ の平均の逆数」)。

## 7. 端の効果 $H_f$

ヘルプ (CAVI の VOLT の説明、原文で確認): 「CAVI は最低次の端の効果を含む。入口 $s_0$ でのハミルトニアン $H_f = -(e(\mathrm{VOLT}+\mathrm{DVOLT})/L)\,(\sin(\omega t - \mathrm{dphi}) + \sin(\mathrm{dphi}) - \mathrm{offset})\,(x^2 + y^2)/4\,\delta(s - s_0)$。$\mathrm{dphi}$ と offset は空洞の位相と放射の損失で決まる (NORAD のときだけ 0 でない)。符号は出口で反転する。既知の導出と整合する」。

### 7.1 Maxwell からの導出

軸の近くの軸対称の場 $E_z(r, s, t) = E_0\,g(s)\cos(\omega t + \varphi_E)$ から。$\nabla\cdot\vec E = 0$ の $r$ 成分 $\tfrac1r\partial_r(rE_r) = -\partial_sE_z$ から $E_r = -\tfrac r2\partial_sE_z$。$\nabla\times\vec B = c^{-2}\partial_t\vec E$ の $z$ 成分 $\tfrac1r\partial_r(rB_\theta) = c^{-2}\partial_tE_z$ から $B_\theta = \tfrac{r}{2c^2}\partial_tE_z$。速度 $\beta c\,\hat s$ の粒子の半径方向の力は $F_r = e(E_r - \beta cB_\theta)$ ($(\vec v\times\vec B)_r = -v_sB_\theta$) なので

$$
F_r = -\frac{e\,r}{2}\Bigl[\partial_sE_z + \frac\beta c\,\partial_tE_z\Bigr]
\tag{10.11}
$$

入口の端では $g(s)$ が 0 から 1 に立ち上がり、$\partial_sE_z = E_0\,g'(s)\cos(\omega t + \varphi_E)$ が支配的である。第 2 項が小さいことを、幅 $w$ の端で積分して示す。$ds = \beta c\,dt$ で $\int F_r\,dt = \int F_r\,ds/(\beta c)$ なので、(10.11) の 2 つの項の積分は

$$
\text{第 1 項: }-\frac{e\,r}{2\beta c}\int\partial_sE_z\,ds = -\frac{e\,r}{2\beta c}\,E_0\cos(\omega t + \varphi_E)\int g'\,ds = -\frac{e\,r\,E_0}{2\beta c}\cos(\omega t + \varphi_E)\quad(\int g'\,ds = 1)
$$

$$
\Bigl|\text{第 2 項: }-\frac{e\,r}{2\beta c}\int\frac\beta c\,\partial_tE_z\,ds\Bigr| = \frac{e\,r}{2c^2}\Bigl|\int(-\omega E_0\,g(s)\sin(\omega t + \varphi_E))\,ds\Bigr| \le \frac{e\,r}{2c^2}\,\omega E_0\,w
$$

($\partial_tE_z = -\omega E_0\,g\sin(\ldots)$、$0 \le g \le 1$、端の幅 $w$ の積分区間。粒子が端の中で動く間の $t$ の変化は、$w \ll c/\omega$ では無視する。) 2 つの比は

$$
\frac{\text{第 2 項}}{\text{第 1 項}} \le \frac{\omega E_0\,w/c^2}{E_0/(\beta c)} = \beta\,\frac{\omega\,w}{c} = \beta\,\frac{2\pi w}{\lambda}\ \to\ 0\quad(w \to 0)
$$

($\lambda = 2\pi c/\omega$ は RF の波長、500 MHz で 0.6 m。) 第 2 項は、端の外側の一様な場 ($g = 1$) の区間にも出るが、それは長さ $L$ の本体の領域の力で、端の寄与ではない (本体の横のキックは §5.1)。したがって端の項は第 1 項だけで、端を通る間の運動量の変化は

$$
\Delta p_r = -\frac{e\,E_0\,r}{2\beta c}\cos(\omega t + \varphi_E)
\tag{10.12}
$$

($\int g'\,ds = 1$。) **ヘルプの記号への読み替えと単位換算**を 1 行ずつ書く。

1. 空洞を通ったときの電場による利得は、$g = 1$ の一様な場で $\Delta E = eE_0\,L\cos(\omega t + \varphi_E)$ ($eE_0L$ は電圧 $eV$ に当たる)。ヘルプの $\Delta E = -eV\sin(\omega t - \mathrm{dphi})$ と等しくおくと $E_0\cos(\omega t + \varphi_E) = -(V/L)\sin(\omega t - \mathrm{dphi})$ (すなわち $\varphi_E = \pi/2 - \mathrm{dphi}$、$E_0 = V/L$。$e$ は電荷の符号を含む)。
2. これを (10.12) に入れる: $\Delta p_r = -\dfrac{e\,r}{2\beta c}\cdot\Bigl[-\dfrac VL\sin(\omega t - \mathrm{dphi})\Bigr] = +\dfrac{e\,V\,r}{2L\,\beta c}\sin(\omega t - \mathrm{dphi})$ (運動量の次元)。
3. **単位を $mc$ に直す** ($c = 1$、運動量は $mc$ 単位、$v = eV/mc^2$): $\Delta p_r/(mc) = \dfrac{e\,V}{mc^2}\cdot\dfrac{r}{2L\,\beta} \sin(\ldots) = \dfrac{v\,r}{2L\,\beta}\sin(\omega t - \mathrm{dphi})$ ($eVr/(2L\beta c)$ を $mc$ で割ると $eV/(mc^2)$ が出て、残りが $r/(2L\beta)$)。
4. $p_0$ (設計運動量、$mc$ 単位) で規格化 (SAD の $p_x$ は $p_0$ で割った量): $\Delta p_r/(p_0\,mc) = \dfrac{v\,r}{2L\,\beta\,p_0}\sin(\omega t - \mathrm{dphi})$。

これは、ヘルプの $H_f$ から $\Delta p_x = -\partial_xH_f$ ($c = 1$ で $v$ を使う) $= (v/L)\,[\sin(\ldots) + \ldots]\,x/2$ ($p_0$ 規格化で $v\,x/(2Lp_0)[\ldots]$) と、$1/\beta$ を除いて同じ形 (定数の $\sin\mathrm{dphi} - \mathrm{offset}$ の項を除く。定数の項は、(10.4) の一定の利得 $v\,[\sin\mathtt{dphis} - \mathtt{offset}]$ の横方向の相手として、(10.7) 以下と同じ議論で出る)。$1/\beta$ の因子が、Maxwell から出る式 (手順 3、4) とヘルプ・K64 の式の違いである。出口では $g' \to -g'$ で符号が反転する (ヘルプ)。

**K64 は $1/\beta$ を略す**: 追跡 `tcavfrin` は `vf=v/al/p0*.5`、`dpt=vf*(sin(ph-dphis)+s0)/p1r`、`px=px+x*dpt` で、$\Delta q_x = x\,\dfrac{v}{2Lp_0}\bigl(\sin(\omega t - \mathtt{dphis}) + s_0\bigr)/p_{1r}$ ($s_0 = \sin\mathtt{dphis} - \mathtt{offset}$。正準運動量では $\Delta p_x/p_0 = x\,\dfrac{v}{2Lp_0}(\ldots)$)。(10.12) は $1/\beta$ を持つ。

**数値** ([cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (d)): 本体のキックが 0 になる位相 ($\omega t - \mathtt{dphis} = \pi/2$ に粒子を置く。本体は $\cos(\omega t - \mathtt{dphis})$ に比例) で、pysad の入口の端だけ (FRINGE = 1 と DISFRIN = 1 の差) を取り出すと、$\Delta p_x$ は K64 の式 $x\,v/(2Lp_0)$ と、3 GeV ($\beta = 1.0000$)、2 MeV/c ($\beta = 0.9689$)、0.6 MeV/c ($\beta = 0.7613$) のすべてで比 $1.000000$ で一致する。**K64 の式に $1/\beta$ は入っていない** (k64 の追跡の結果。v2 は $1/\beta$ を入れてあり、比は 2 MeV/c で 1.032、0.6 MeV/c で 1.314 ($= 1/\beta$)。[cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (d))。

**確認 (K64 の誤り。ローレンツ力の数値積分と gs で確認、v2 は修正 = B330。[cav_02_edge_lorentz.out.txt](outputs/cav_02_edge_lorentz.out.txt))**:

1. **gs**: $L = 10^{-4}$ m、VOLT = 10 V、500 MHz、PHI = 0.3 の CAVI 1 要素の `TrackParticles` の端の $p_x$ ($x = 10^{-5}$、FRINGE = 1 と DISFRIN = 1 の差) は、$x\,v\sin(\mathrm{PHI})/(2Lp_0)$ に 3 GeV/c・2 MeV/c・0.6 MeV/c の全部で比 1.000000 ($1/\beta$ は入っていない)。pysad の k64 と 12 桁一致 (出力 (d))。
2. **ローレンツ力の数値積分** (DOP853、軸の近くの場 $E_r = -\tfrac r2\partial_sE_z$、$B_\theta = \tfrac{r}{2c^2}\partial_tE_z$、$E_z$ の $r^2$ の項、端の幅 $10^{-5}$ m、$E_0 = 10^5$ V/m、電子): 追跡 / 数値積分の入口の $p_x/x$ の比は、k64 が 0.6 MeV/c で 0.754〜0.760 ($\beta = 0.7613$)、2 MeV/c で 0.962〜0.968 ($\beta = 0.9689$) で、**$1/\beta$ が正しい** (出力 (a))。端のエネルギー変化の $x^2$ の項 ($\Delta h = \partial_tH_f$) も同じ $1/\beta$ (k64 の比 0.762・0.970、v2 は 1.0009・1.0007。出力 (b))。$\delta \ne 0$ の粒子は設計の $\beta$ でなく**粒子の $\beta = p/h$** に追随する (0.6 MeV/c、$\delta = -0.05 \ldots +0.15$ で k64 の比 0.737〜0.795、v2 は 0.991 で一定。出力 (c))。
3. **v2 の修正**: キックとエネルギー変化に $1/\beta$ (粒子の $\beta$) をかけ、$H_f \propto 1/\beta(E)$ の整合の時間のずれ $\Delta t = -\partial H_f/\partial E = -(x^2 + y^2)\kappa_0/(2p^3)$ ($\kappa_0 = v_fp_0(\sin(\omega t - \mathtt{dphis}) + s_0)$、$\partial(1/\beta)/\partial h = -1/p^3$) を入れる。時間のずれが無いとヤコビアンのシンプレクティック条件の誤差が 2 桁悪化する (6.9e-4 対 1.7e-6。$q$ の再スケール (§7.2、B332) も入れた値。0.6 MeV/c、$x$ = 3 mm、$\delta$ = 0.03)。

入口と出口の端が同じ位相なら打ち消し合うので、効くのは端 1 つ・$\delta \ne 0$ や $z \ne 0$ で入口と出口の位相が違う場合。3 GeV/c では $1/\beta - 1 = 1.5 \times 10^{-8}$。

**限界**: 端 1 つの 1 次の収束とエネルギー変化の $x^2$ の項だけ (1% 以内)。長い空洞 ($\omega L/\beta$ が大きい) の本体との組合せ、定在波の位相の進み (K64 は基準粒子の位相が入口・中央・出口で同じ) は見ていない。時間のずれの項は $r^3$ までの場が要るので積分では確かめられず、シンプレクティック条件で決めた。

(1/β の有無の判定に使う測定は、cav_02 (d) の上の結果だけである。[cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (f) は、位相を選ばずに同じ量を測った参考の測定で、k64 の比が 0.9993 (2 MeV/c)、0.9734 (0.6 MeV/c) と 1 からずれる (v2 は 1.0314、1.2786)。これは端の式の性質ではなく、端のキックで動いた $x$ に、本体の $x$ 依存のキック ($V_{20}$・$V_{02}$ の補正) が反応する分が混ざった測定のずれで、上のように本体のキックが 0 の位相に置くと消える。(f) の出力の「$1/\beta$ をかけた式」の列は参考値で、判定には使わない。)

### 7.2 端のエネルギー変化と、$q$ の再スケールの欠落 (K64 の誤り。gs で確認、v2 は修正 = B332)

**端のエネルギー変化の符号と係数** (ヘルプの $H_f$ と整合する)。入口の端の項を、$s$ で積分したキックのハミルトニアン $K(x, y, t) = -(v/L)\,[\sin(\omega t - \mathtt{dphis}) + s_0]\,(x^2 + y^2)/4$ と書く ($c = 1$、エネルギーは $mc^2$ 単位。$v = V/mc^2$)。SAD のキックの規約 (§5.1: $\Delta E = \partial_tK$、$\Delta p_x = -\partial_xK$) を当てると

$$
\Delta E = \partial_tK = -\frac{\omega\,v\,(x^2 + y^2)}{4L}\cos(\omega t - \mathtt{dphis}),\qquad
\frac{\Delta p_x}{p_0} = -\frac{\partial_xK}{p_0} = \frac{v\,x}{2Lp_0}\bigl[\sin(\omega t - \mathtt{dphis}) + s_0\bigr]
$$

コードの `dh=-w*p0*vf*(x**2+y**2)*.5*cos(ph-dphis)` ($v_f = v/(2Lp_0)$ なので $-\omega p_0\,v_f\,(x^2 + y^2)/2 = -\omega v(x^2 + y^2)/(4L)$) と、`dpt` の分子 $v_f(\sin + s_0)$ は、この 2 式と一致する (sympy で差 0: [cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (a))。出口は $v \to -v$ (`tcavfrin(…, -v, …)`)。入口は `mfring ≥ 0` かつ `mfring ≠ 2` (`FRINGE` = 0, 1, 3)、出口は `mfring ≥ 0` かつ `mfring ≠ 1`。`fringe` は DISFRIN が 0 のとき。**ヘルプの FRINGE・DISFRIN の表** (DISFRIN ≠ 0 は端なし。FRINGE = 1 は入口、2 は出口、3 は両方、0 も両方): 原文で確認。K64 の `mfring` への変換 (`tffs.f:1447-1451` の `frmd`。§2) は、通常の向きの要素では FRINGE をそのまま `mfring` とし、逆向きの要素では 1 と 2 を入れ替える。FRINGE = 0 は `mfring` = 0 で入口・出口の両方の条件 (`mfring ≥ 0`) を満たし、ヘルプの表 (0 は両方) と整合する。

**`tcavfrin` の手続き** (`tcav.f:218-247`): (1) $h_1$ を作り、$t = -z\,h_1/p_1$、$\omega t$ を計算 (`ph=w*t`。本体と違い $\omega t$ そのもの)。(2) 傾きの kick $q \mathrel{+}= x\,\mathtt{dpt}$。(3) エネルギー変化 $dh$ から $p_{2r}$ を (10.5) と同じ形で計算、`g(i)=dp2r`。(4) $z \mathrel{*}= (p_{2r}/p_{1r})(h_1/h_2)$ ((10.6) と同じ。時間 $t$ を保つ)。(5) $dv$ を更新。

**欠落している操作 (K64 の誤り。確定。gs で確認、v2 は修正 = B332)**: 本体 `tcav` は、エネルギーが変わったあと `px=(px*p1r+dpx)/p2r` ($q_2 = (q_1p_{1r} + \Delta p_x/p_0)/p_{2r}$。正準運動量 $q\,p_r$ にキックを足してから新しい $p_{2r}$ で割る) とする。`tcavfrin` は $dh \ne 0$ で $p$ を変えるのに、**$q$ に $p_{1r}/p_{2r}$ をかけない** (キックだけ足す)。正しい式との差は

$$
q_2^{\rm correct} - q_2^{\rm K64} = \frac{q_1p_{1r} + X}{p_{2r}} - \Bigl(q_1 + \frac{X}{p_{1r}}\Bigr) = (q_1p_{1r} + X)\,\frac{p_{1r} - p_{2r}}{p_{1r}\,p_{2r}},\qquad X \equiv x\,v_f\bigl[\sin(\ldots) + s_0\bigr]
$$

($q_1 + X/p_{1r} = (q_1p_{1r} + X)/p_{1r}$ を使った。) この欠落が $J^\top SJ - S$ にどう入るかを、成分で書く。

**(1) 出力の正準 $p_x$ の差 $\eta$**。出力の正準運動量 ($p_0$ 規格化) は $q_2\,p_{2r}$ なので、再スケールありとなしの差は ($p_{x} = q_1p_{1r}$ は入力の正準 $p_x$、$X$ は上の定義の端のキック)
$$
\eta \equiv p_{2r}\,(q_2^{\rm correct} - q_2^{\rm K64}) = (q_1p_{1r} + X)\,\frac{p_{1r} - p_{2r}}{p_{1r}} = -(p_x + X)\,\mathtt{dpr}
\tag{10.13}
$$
($p_{2r} - p_{1r} = p_{1r}\,\mathtt{dpr}$、すなわち $p_{1r} - p_{2r} = -p_{1r}\,\mathtt{dpr}$ を使った。) 同じ式が $p_y$ にもある ($y$、$p_y$、$Y \equiv y\,v_f[\ldots]$)。$\mathtt{dpr} \simeq a/2 = \Delta h\,(h_1 + h_2)/(2p_1^2) \simeq (h/p^2)\,\Delta h$ ($h_2 \simeq h_1$)、端のエネルギー変化は $\Delta h = -A\cos\theta$、$A \equiv \omega p_0v_f\,(x^2 + y^2)/2$、$\theta \equiv \omega t - \mathtt{dphis}$ ($dh$ の式は上のとおり)。

**(2) ヤコビ行列の差**。再スケールありの写像のヤコビ行列を $J_c$ (本体と同じ扱い。シンプレクティック)、K64 の式のを $J_K$ とすると、$p_x$ の出力だけが $\eta$ だけ違うので $J_K = J_c - e_{p_x}\,g^\top$ ($g \equiv \nabla\eta$、$e_{p_x}$ は $p_x$ 成分の単位ベクトル)。$J_K^\top SJ_K - S$ を展開する ($J_c^\top SJ_c = S$ を使う):
$$
J_K^\top SJ_K - S = -J_c^\top S\,e_{p_x}\,g^\top - g\,e_{p_x}^\top SJ_c + g\,(e_{p_x}^\top Se_{p_x})\,g^\top
$$
$e_{p_x}^\top Se_{p_x} = 0$ (反対称)。$S$ の $(x, p_x)$ 成分は $S_{x\,p_x} = +1$、$S_{p_x\,x} = -1$ なので $Se_{p_x} = e_x$、$e_{p_x}^\top S = -e_x^\top$。端は $x$ を変えない ($x' = x$) ので $J_c$ の $x$ の行は $e_x^\top$ で、$J_c^\top e_x = e_x$、$e_x^\top J_c = e_x^\top$。したがって
$$
J_K^\top SJ_K - S = g\,e_x^\top - e_x\,g^\top\ (+\ y\ \text{方向の同じ項}\ g_y\,e_y^\top - e_y\,g_y^\top)
$$
(**反対称**。$g_x e_xe_x^\top$ は打ち消す。) すなわち欠陥は、$x$ の行・列に、$\eta$ の各座標での偏微分 $g_k = \partial\eta/\partial k$ ($k = p_x, z, \delta, y, p_y$) として入る。

**(3) 各成分の大きさ**。(10.13) から ($X$ の微分は $\mathtt{dpr}$ の分だけ小さく、高次):
$$
g_{p_x} = \frac{\partial\eta}{\partial p_x} = -\mathtt{dpr}\ \ (\text{ほぼ正確}),\qquad
g_z \simeq -(p_x + X)\,\partial_z\mathtt{dpr},\qquad \partial_z\mathtt{dpr} \simeq \frac h{p^2}\,\partial_z\Delta h = -\frac h{p^2}\,\frac{\omega A}{\beta}\sin\theta
$$
($\partial_z\theta = -\omega h_1/p_1 = -\omega/\beta$ から $\partial_z(-A\cos\theta) = A\sin\theta\cdot(-\omega/\beta)$。) $g_\delta$ も同様に $\partial_\delta\mathtt{dpr}$ に比例する。$g_{p_x}$ は $\mathtt{dpr} \sim (h/p^2)\,A$ そのもの (1 次) で、$g_z, g_\delta$ は $(p_x + X)\cdot(h/p^2)\,\omega A$ (さらに $p_x$ か $X$ と $\omega\times$長さの積) の大きさ。数値 (下の表の条件) では $g_{p_x}$ が最大で、$|J^\top SJ - S| \simeq |\mathtt{dpr}| = (h/p^2)\,A\,|\cos\theta|$、**傾き $p_x$ にはよらない** ($p_x = 0$ でも同じ値で残る。出力 (d))。$x^2$ に比例 (1 mm → 10 mm で 100 倍)、$h/p^2$ が大きい低エネルギーで大きくなる。**傾き $p_x$ の大きさで決まるのは $g_z, g_\delta$ だけで、これらは $|D|$ の主要部ではない。**

数値確認 ([cav_03_symplectic_proofs.out.txt](outputs/cav_03_symplectic_proofs.out.txt) (d)): (10.13) の $\eta$ は数値の差と一致 (例: 2 MeV/$c$、$x = 1$ mm、$p_x = 0.01$ で $1.306 \times 10^{-8}$、式の値も同じ)。$g_{p_x} = -\mathtt{dpr}$ (2 MeV/$c$、$x = 1$ mm で $1.297 \times 10^{-6}$)。$J_K^\top SJ_K - S$ の数値の行列と、$g$ の数値微分から組んだ $g\,e_x^\top - e_x\,g^\top + (y\ \text{項})$ の差は、全 12 条件 (3 つの運動量、$x = 1, 10$ mm、$p_x = 0.01, 0$) で $2 \times 10^{-11}$ 以下。

**数値** ([cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (b)(c)。CAVI、$L = 0.1$ m、VOLT = $10^5$ V、500 MHz、PHI = 0.5、$z = 0.02$ m、$\delta = 10^{-3}$、$(x, p_x) = (1\,\mathrm{mm}\ \text{または}\ 10\,\mathrm{mm},\ 0.01)$、中心差分): 入口の端だけを numpy で組み直した写像 ((b)。中心差分の精度は $10^{-9}$〜$10^{-10}$) の $|J^\top SJ - S|$ は、$q$ の再スケールなし (K64 の式) と、あり (本体と同じ扱い) で

| 運動量 | $x$ | 再スケールなし | あり |
|---|---|---|---|
| 3 GeV/c | 1 mm / 10 mm | $8.4 \times 10^{-10}$ / $8.8 \times 10^{-8}$ | $4.8 \times 10^{-10}$ / $4.0 \times 10^{-9}$ |
| 2 MeV/c | 1 mm / 10 mm | $1.3 \times 10^{-6}$ / $1.3 \times 10^{-4}$ | $2.4 \times 10^{-9}$ / $2.3 \times 10^{-8}$ |
| 0.6 MeV/c | 1 mm / 10 mm | $5.6 \times 10^{-6}$ / $5.6 \times 10^{-4}$ | $1.5 \times 10^{-8}$ / $2.0 \times 10^{-8}$ |

(再スケールを入れると、中心差分の精度の範囲の $10^{-8}$ 以下に下がる。) pysad の k64 の端つき CAVI 全体 ((c)。入口と出口の端、本体つき) では、フリンジあり $4.8 \times 10^{-8}$ (3 GeV、$x = 1$ mm)、$4.0 \times 10^{-6}$ (2 MeV/c)、$1.8 \times 10^{-5}$ (0.6 MeV/c)、$x = 10$ mm でそれぞれ $1.8 \times 10^{-8}$、$3.2 \times 10^{-5}$、$1.5 \times 10^{-4}$、DISFRIN = 1 では $10^{-8}$ 以下 ($3 \times 10^{-10}$〜$1.8 \times 10^{-7}$)。(これは K64 の式のまま、$q$ を再スケールしない値。)**v2 は $q$ の再スケールを入れ**、フリンジありは $x = 1$ mm で 3 GeV/c $4.8 \times 10^{-8}$、2 MeV/c $2.4 \times 10^{-8}$、0.6 MeV/c $5.4 \times 10^{-7}$、$x = 10$ mm で $5.6 \times 10^{-8}$、$5.9 \times 10^{-7}$、$1.4 \times 10^{-4}$ ([cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) (c) の v2。$x = 10$ mm の 0.6 MeV/c の $1.4 	imes 10^{-4}$ は k64 と同じ桁で残っており、原因は未調査)。

**gs で確認 (B332)**: gs の `TrackParticles` の 13 点の中心差分ヤコビアン ($L = 0.1$ m、VOLT = $10^5$ V、500 MHz、PHI = 0.5、FRINGE = 1、点 $(x, p_x, y, p_y, z, \delta) = (10^{-3}, 10^{-2}, -5 \times 10^{-4}, 5 \times 10^{-3}, -0.02, 0.01)$) の $\|J^\top SJ - S\|$ は、gs で 3 GeV/c が $1.3 \times 10^{-8}$、0.6 MeV/c が $8.1 \times 10^{-5}$ (pysad の k64 と全桁一致)。v2 は再スケールを入れ、0.6 MeV/c で $5.8 \times 10^{-8}$。限界: 1 点の中心差分のヤコビアンだけ。

**確認の状況**: K64 の `tcavfrin` に再スケールが無いことはソースの読み (確定)、非シンプレクティック性は gs のヤコビアンと pysad の再実装で確認した。(c) の数値は、DISFRIN = 1 (端なし) の値が $10^{-8}$ 以下〜$1.8 \times 10^{-7}$ と一様でなく、追跡側の中心差分の雑音を含む。(c) と、入口の端だけの (b) の欠陥 ((10.13) と $|\mathtt{dpr}|$ の見積り) との定量的な比較 (入口と出口の欠陥の重なり方を含む) は行っていない。

## 8. TCAVI (偏向空洞)

ヘルプに TCAVI の節は無い (キーワード一覧 `L ROTATE DX DY K0 V1 FREQ PHI HARM RANKICK RANPHASE` と RFSW の説明だけ。原文で確認)。`ttcav.f` (追跡):

1. $L \ne 0$ なら半分のドリフト。
2. $v = \mathrm{K0}\cdot p_0\,\omega$、$t = -z\,h_1/p_1$ ($|t| \le 10^4$)、位相 $ph = \omega t + \phi$ ($\phi = \mathrm{PHI}\cdot\mathrm{sign}(\mathtt{charge})$、周回ごとに HARM の端数 × $2\pi$ を足す `tturn.f`)。
3. **キック**: $\Delta p_x = -\mathrm{K0}\,\sin(ph)$ (内部 $q$ では $/p_{1r}$)、**エネルギー**: $\Delta h = v\,x\cos(ph) = \mathrm{K0}\,p_0\,\omega\,x\cos(ph)$。運動量・$dv$・$z$ は §4 と同じ (10.5)(10.6)、$q$ は $p_{1r}/p_{2r}$ 倍。
4. 半分のドリフト (放射があれば `tradk`)。

(1) キックと (3) のエネルギー変化は、ハミルトニアン $K(x, z, p) = \mathrm{K0}\,p_0\,x\,\sin(\phi - \omega z/\beta(p))$ ($\beta(p) = p/\sqrt{1 + p^2}$) の時間 1 の流れ: $\Delta p_x = -\partial_xK = -\mathrm{K0}\,p_0\sin(\ldots)$ ($p_0$ 規格化で $-\mathrm{K0}\sin(ph)$)、$\Delta p = -\partial_zK = \mathrm{K0}\,p_0\,x\,(\omega/\beta)\cos(\ldots)$ で、$\Delta h = \beta\,\Delta p = \mathrm{K0}\,p_0\,\omega\,x\cos(ph)$ と一致する。つまり **Panofsky–Wenzel の関係が成り立つ**。数値確認 ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (e)): 追跡は上の式と $6 \times 10^{-17}$ で一致、$|J^\top SJ - S| = 9.7 \times 10^{-11}$ (シンプレクティック)、追跡と $K$ の流れ (DOP853) の差は 0 ($\delta$、$p_x$ について。3 GeV では $z$ の更新が $3.7 \times 10^{-16}$ しか動かず、「時間を保つ」$z$ の更新と $\partial K/\partial p$ の違いは見えない。未確認)。`RANKICK`・`RANPHASE` は K0・PHI に Gauss 乱数を足す (`tturn.f:470-480`)。CALC (`ttcave`) は、同じ式の行列を作る (`ttcave.f`: `trans1(2,5) = ak*w/v1*cosp`、`trans1(6,1) = ak*w*cosp/v2` など)。

## 9. 確認

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `cav_02_derivations.py` | (a) 端の項 $H_f$ から $\Delta E = \partial_tK$・$\Delta p_x = -\partial_xK$ を取り `tcavfrin` の式と照合 (sympy)、(b)(c) `tcavfrin` の $q$ の再スケールの有無とヤコビアン (numpy の再実装と pysad の k64・v2)、(d) 本体のキックが 0 の位相での端の収束の $\beta$ 依存 (k64・v2)、(e) $wl = \ln(1 + r)$、`dpepe` | [cav_02_derivations.out.txt](outputs/cav_02_derivations.out.txt) |
| `cav_02_edge_lorentz.py` | 端の収束・端のエネルギー変化・$\delta \ne 0$ を、軸の近くの Maxwell の場のローレンツ力の数値積分と比べる (k64・v2)。gs の保存した出力との比 (§7.1。B330) | [cav_02_edge_lorentz.out.txt](outputs/cav_02_edge_lorentz.out.txt) |
| `cav_03_symplectic_proofs.py` | (a) $(t, E) \to (t, E + f(t))$ と $(t, E) \to (z, \delta)$ のヤコビ行列 ((10.4a)(10.4b)、sympy)、(b) 横のキックの符号 $s$ と $D = J^\top\Omega J - \Omega$ (sympy)、(c) $s = \pm1$ の薄い空洞の数値ヤコビアンと予想の比較、(d) `tcavfrin` の欠陥 $J^\top SJ - S = g\,e_x^\top - e_x\,g^\top + (y\text{ 項})$ と (10.13) の確認 | [cav_03_symplectic_proofs.out.txt](outputs/cav_03_symplectic_proofs.out.txt) |
| `cav_01_energy_map.py` | (a) 薄い CAVI の写像 (10.1)(10.4)〜(10.6) と追跡の一致、(b) シンプレクティック性 (V1・V20・フリンジつき)、(c) 非相対論の補正と $I_0$、(d) 重み `ws` の和と $e^{wl} = 1 + r$、(e) TCAVI の式・シンプレクティック性・$K$ の流れ、(f) 端の収束の測定 (位相を選ばない。参考。$1/\beta$ の判定は cav_02 (d)) | [cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/cav_01_energy_map.py` (数秒。numba のキャッシュは環境変数 `NUMBA_CACHE_DIR` で指定できる)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子、運動量 3 GeV/$c$ (§3〜§5、§8)、低エネルギーの確認は 0.6〜2 MeV/c (§5.2、§7)。VOLT = 2 MV (§3、§4)、FREQ = 500 MHz、PHI = 0.3 rad (§5.1 の符号の確認は V1 = $2 \times 10^6$ V/m などと dphis = $-0.3$、§7.2 の端の欠陥は VOLT = $10^5$ V、$L = 0.1$ m、PHI = 0.5)、$\delta = O(10^{-3})$、$z = O(0.1)$ m。放射・スピン・ウェイク (`txwake`)・TRPT (`vnominal ≠ 0`) の追跡は確認していない。

**この章の範囲外**: ウェイク場 (`txwake`、WAKE の章)、放射 (RADCOD、`tradk`、放射による損失と `trf0` の決め方)、`tmultiacc` (MULT + RF)、CALC (`tcave`) の行列の導出と追跡との比較、COORD・ずれ (DX・DY・ROTATE) の `TENT`/`TEXIT`。

**未確認・限界**:
- **TRPT (`vnominal ≠ 0`) の写像**は、読みだけ。$\mathtt{ws}$ の和と $e^{wl} = 1 + r$ は数値で確認したが、(10.10) の等比の重みを選ぶ理由、`ndiv` の式の係数 ($\sqrt{r^2/3\varepsilon}$ の 3 と $\varepsilon$。§6.1) の根拠、$dz_n$ の数値確認は行っていない。$wl = \ln(1 + r)$ と $dp_{epe}$ の式は導出した (§6.2)。
- **端の収束の $1/\beta$** (§7.1): ローレンツ力の数値積分と gs で確認済み (K64 の誤り、v2 は修正。B330。端の幅 $w \to 0$ の極限で第 2 項が落ちることは、上限の式で示した)。未確認: 長い空洞の本体との組合せ、定在波の位相の進み、時間のずれの項の積分による確認 (シンプレクティック条件で決めた)。
- **端のエネルギー変化**の符号と係数は、ヘルプの $H_f$ と整合する (§7.2、sympy)。**`tcavfrin` の $q$ の再スケールの欠落**は、gs のヤコビアンで確認済み (K64 の誤り、v2 は修正。B332。§7.2、付録 A)。欠陥の構造 ($J^\top SJ - S = g\,e_x^\top - e_x\,g^\top$) と大きさ ($|\mathtt{dpr}|$) は入口の端だけの再実装で確認したが、pysad の端つき CAVI 全体 (入口・出口・本体) での欠陥との定量的な比較は未実施。
- **設計粒子の $z$ の補正 `dzn`** (§6.2) の効果の数値確認、等比の重みを選ぶ理由。
- **`offset` の意味** (§3.3): ソースの読み。$\mathtt{phis} \ne 0$ のときの一定のエネルギー利得の物理的な意味は確認していない。
- **V1・V11・V02 は「Tracking only」** (ヘルプ) と書かれるが、`tcave.f` (CALC) は `v10`・`v20`・`v11`・`v02` を使う。CALC の V1 の行列の係数は確認していない。
- **ヘルプの V20・V02 の説明**: 「The non-relativistic corrections ... (VOLT+DVOLT)(2 Pi FREQ/c)^2/(gamma beta)^2/4 are automatically added to V20 and V02, respectively」は、足される量が V20 と V02 の係数 ($V_{20}x^2/2$ の $1/2$ を含む) の 2 倍の大きさ、すなわち $\Delta V = V(\omega/\beta\gamma c)^2(x^2 + y^2)/4$ を意味する。コードでは `v20a` ($= V_{20}/2$ の係数) に足す。文面は読み方で 2 通りあり、コードと数値 (§5.2) から後者。

## 付録 A. ヘルプの記述・K64 との違い (CAVI・TCAVI)

原文 (`SADHelp.html` の `CAVI`、`RFSW`、`RADCOD`、`TRPT`、キーワード一覧) で確認した。

| 箇所 | 内容 |
|---|---|
| CAVI の $\Delta E$ の式 | **原文で確認**: $\Delta E = -e(\mathrm{VOLT} + \mathrm{DVOLT} + V_1x + \tfrac12(V_{20}x^2 + V_{02}y^2) + V_{11}xy)\sin(2\pi\,\mathrm{FREQ}(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$ (PHI・DPHI・FREQ・VOLT の各項で同じ式)。§3 のとおり (10.4) で実装される。「$t_s$ は放射と加速のつりあいで決まる同期時間」。 |
| CAVI の VOLT の説明: 非相対論の補正 | **原文で確認、解釈を確定**: §5.2 (数値で一致)。 |
| CAVI の端の効果 $H_f$ | **原文で確認**: $H_f = -(e(\mathrm{VOLT} + \mathrm{DVOLT})/L)(\sin(\omega t - \mathrm{dphi}) + \sin\mathrm{dphi} - \mathrm{offset})(x^2 + y^2)/4\,\delta(s - s_0)$。**ヘルプは $\sin$ の中を $\omega t - \mathrm{dphi}$ と書き、DPHI の説明は $t_s$ の項を持つ形 ($\omega(t - t_s) + \ldots$)**。コードの $\sin(\omega t - \mathtt{dphis})$ ($\mathtt{dphis} = \mathtt{phis} - \mathtt{phic}$) と記号の定義が違う (ヘルプの dphi は `dphis` に当たると読める。確認していない)。Maxwell からの導出 (§7.1) と形は同じで、$1/\beta$ の因子が K64 に無い (§7.1。確認済み、v2 は修正)。 |
| CAVI の DISFRIN・FRINGE の表 | **原文で確認**: DISFRIN ≠ 0 は端なし。FRINGE = 1 入口、2 出口、3 入口と出口、0 は両方。他の要素 (QUAD) の表と違い、CAVI の表は「Linear」の列が無い。 |
| CAVI の DPHI の説明 | **原文で確認**: 「DPHI is not taken into account to determine the design momentum $p_0(s)$」。コードで `vnominal` は PHI だけを使う (`tffs.f` の `vnominal=…sin(-cmp%value(ky_PHI_CAVI)…)`)。一致。 |
| RADCOD・TRPT・RFSW | **原文で確認** (§5.3、§6.2)。TRPT の既定の運動量分布は「一様」(Emittance[] での扱い) は未確認。 |
| TCAVI | **原文に節が無い** (キーワード一覧と RFSW の言及のみ)。式はすべて `ttcav.f`・`ttcave.f` から。 |
| ヘルプに無いもの | (1) 薄い空洞の運動量・$z$・横運動量の更新 (§4)。(2) Panofsky–Wenzel の関係 (§5.1)。(3) 長い空洞の分割 (§6)、重み。(4) 端の導出 (§7)。(5) TCAVI の式。 |

**K64 の誤り (確定。v2 は修正)**:
1. 端の収束・エネルギー変化に $1/\beta$ が無い (§7.1。ローレンツ力の数値積分と gs で確認。B330)。低エネルギーの空洞では、ヘルプの「lowest order」の範囲を超える。
2. **`tcavfrin` が $q$ を再スケールしない (gs のヤコビアンで確認。B332)**: 端のエネルギー変化で $p$ が変わるのに、本体 `tcav` の `px=(px*p1r+dpx)/p2r` に相当する操作が無く、端つきの空洞が非シンプレクティックになる (§7.2。ソースの読み + gs と pysad のヤコビアンの数値。入口の端だけの写像で、欠陥は $|\mathtt{dpr}| \simeq (h/p^2)\,|\Delta h|$、$p_x$ によらず $x^2$ に比例: $x = 1$ mm で 3 GeV/$c$ は $8 \times 10^{-10}$、2 MeV/$c$ は $1.3 \times 10^{-6}$、0.6 MeV/$c$ は $5.6 \times 10^{-6}$)。入口の端だけの写像で、欠陥は $|\mathtt{dpr}| \simeq (h/p^2)\,|\Delta h|$、$p_x$ によらず $x^2$ に比例: $x = 1$ mm で 3 GeV/$c$ は $8 \times 10^{-10}$、2 MeV/$c$ は $1.3 \times 10^{-6}$、0.6 MeV/$c$ は $5.6 \times 10^{-6}$)。疑い・取り下げは無い。

## 付録 B. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `v` | $V/(mc^2)\cdot|\mathtt{charge}|$ (電圧を静止エネルギーで割った量) | `tcav.f`: `v=vc/amass*abs(charge)` | `v` |
| `w`・`wi` | $\omega = 2\pi f/c$ [rad/m]、$1/\omega$ | `tcav.f`: `w`、`wi=1.d0/w` | `w`・`wi` |
| `phic`・`phis`・`dphis` | $(\mathrm{PHI} + \mathrm{DPHI})\,\mathrm{sign}(\mathtt{charge})$、$\omega\,t_{rf0}$、$\mathtt{phis} - \mathtt{phic}$ | `tcav.f` | `phic`・`phis`・`dphis` |
| `trf0` | 同期時間 $t_s$ (global。RF の位相の基準) | `tmacro` | `trf0` |
| `offset`・`offset1` | 一定の利得の調整 ($\sin\mathtt{dphis}$ または $\sin\mathtt{dphis} - \sin\mathtt{phis}$、$0$ または $\sin\mathtt{phis}$) | `tcav.f` | `offset`・`offset1` |
| `amass` | 粒子の静止エネルギー $mc^2$ [eV] (電圧を割って $v$ を作る) | `tfsetparam.f:70` (`amass=rgetgl1('MASS')`) | `amass` |
| `charge` | 電荷 (電子 $-1$)。`abs(charge)` が $v$ に、`sign(charge)` が位相に入る | `tfsetparam.f` (`rgetgl1('CHARGE')`) | `charge` |
| `pgev`・`tphyzp` | $p_0c$ [eV] (名前は gev だが eV)。`tphyzp` は `pgev` から $p_0$ ($=$ `pgev/amass`)、$h_0$、磁気剛性、放射の定数、$\omega_0$ を作り直すサブルーチン (TRPT で空洞のあとに呼ぶ) | `tfsetparam.f:105` (`tphyzp`)、`tcav.f:202-203` | `p0`・`h0` ほか |
| `rad`・`trpt`・`autophi` | 論理変数: 放射 (RADCOD)、輸送路 (TRPT)、自動位相 (`APHI` ≠ 0)。`offset`・`phis` の場合分け (§3.3) | `tffs.f:585-589` (`fff%rad`、`fff%trpt`)、`tturn.f:438` (`autophi`) | `rad`・`trpt`・`autophi` |
| `mfring`・`frmd` | 端の有無の整数。`frmd` = `FRINGE`、逆向きの要素では 1 ↔ 2 を入れ替えて `mfring` にする。入口は `mfring ≥ 0` かつ `≠ 2`、出口は `≥ 0` かつ `≠ 1` | `tffs.f:1447-1451`、`tcav.f:82`・`185` | `mfring` |
| `dvfs`・`al`・`alx` | RFSW の周波数ずらしによる $dv$ への足し込み、空洞の長さ $L$、刻みの中のドリフトの長さ | `tturne.f:859-866` (`tsetdvfs`)、`tcav.f` | `dvfs`・`al` |
| `h1`・`p1`・`p1r` | 入口のエネルギー ($mc^2$)、運動量 ($mc$)、$1 + \delta$ | `tcav.f` | `amass` | 粒子の静止エネルギー $mc^2$ [eV] (電圧を割って $v$ を作る) | `tfsetparam.f:70` (`amass=rgetgl1('MASS')`) | `amass` |
| `charge` | 電荷 (電子 $-1$)。`abs(charge)` が $v$ に、`sign(charge)` が位相に入る | `tfsetparam.f` (`rgetgl1('CHARGE')`) | `charge` |
| `pgev`・`tphyzp` | $p_0c$ [eV] (名前は gev だが eV)。`tphyzp` は `pgev` から $p_0$ ($=$ `pgev/amass`)、$h_0$、磁気剛性、放射の定数、$\omega_0$ を作り直すサブルーチン (TRPT で空洞のあとに呼ぶ) | `tfsetparam.f:105` (`tphyzp`)、`tcav.f:202-203` | `p0`・`h0` ほか |
| `rad`・`trpt`・`autophi` | 論理変数: 放射 (RADCOD)、輸送路 (TRPT)、自動位相 (`APHI` ≠ 0)。`offset`・`phis` の場合分け (§3.3) | `tffs.f:585-589` (`fff%rad`、`fff%trpt`)、`tturn.f:438` (`autophi`) | `rad`・`trpt`・`autophi` |
| `mfring`・`frmd` | 端の有無の整数。`frmd` = `FRINGE`、逆向きの要素では 1 ↔ 2 を入れ替えて `mfring` にする。入口は `mfring ≥ 0` かつ `≠ 2`、出口は `≥ 0` かつ `≠ 1` | `tffs.f:1447-1451`、`tcav.f:82`・`185` | `mfring` |
| `dvfs`・`al`・`alx` | RFSW の周波数ずらしによる $dv$ への足し込み、空洞の長さ $L$、刻みの中のドリフトの長さ | `tturne.f:859-866` (`tsetdvfs`)、`tcav.f` | `dvfs`・`al` |
| `h1`・`p1`・`p1r` |
| `dh`・`h2`・`a`・`dpr`・`p2r` | $\Delta h$、$h_1 + \Delta h$、$\Delta h(h_1 + h_2)/p_1^2$、$\sqrt{1 + a} - 1$、$p_{1r}(1 + \mathtt{dpr})$ | `tcav.f` | 同名 |
| `t`・`ph`・`tlim` | $t = -z\,h_1/p_1$、$\omega t/2$、$|t|$ の上限 (RFSW OFF で 0) | `tcav.f` | 同名 |
| `v1a`・`v20a`・`v11a`・`v02a` | $V_1$、$V_{20}/2$、$V_{11}$、$V_{02}/2$ ($/\mathtt{amass}$) | `tcav.f` | 同名 |
| `vcorr` | $v\,(\omega\,(0.5/p_0 + 0.5/p_e))^2/4$ (非相対論の補正) | `tcav.f` | `vcorr` |
| `av`・`dpx`・`dpy` | $-\cos(2\,ph - \mathtt{dphis})/\omega\,\ldots/p_0$、横のキック | `tcav.f` | 同名 |
| `vnominal`・`he`・`pe`・`dpepe` | 設計の利得 $v\sin(-\mathrm{PHI})$、$h_0 + \mathtt{vnominal}$、$h_e$ の運動量、再基準化 | `tcav.f` | 同名 |
| `ndiv`・`ws`・`wsn`・`wl`・`r` | 刻みの数、重み、その刻みの重み、$wl$、$r = |\mathtt{vnominal}(1/h_0 + 1/h_e)|$ | `tcav.f` | 同名 |
| `sv`・`dzn`・`dvn` | 設計エネルギーの累積利得、設計の $z$ のずれ、設計の $dv$ | `tcav.f` | 同名 |
| `vf`・`s0`・`dpt` | 端の係数 $v/(2Lp_0)$、$\sin\mathtt{dphis} - \mathtt{offset}$、端の収束 | `tcavfrin` (`tcav.f`) | `vf`・`s0` |
