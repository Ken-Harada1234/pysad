# 加速空洞 CAVI と TCAVI

**この章で分かること**: 高周波加速空洞 CAVI の追跡 (K64 の `tcav`) が、ΔE の式 (ヘルプ) から $\delta$・$z$・$p_x$ をどう更新するか、(1) 位相と時間 ($t = -z\,h_1/p_1$)、エネルギー利得と、RFSW・RADCOD・同期位相 ($\mathtt{phis}$、`offset`) の関係、(2) エネルギーが変わったあとの運動量・$z$・$dv$・横運動量の更新 (時間を保つ)、(3) 横方向の依存 (V1、V20、V11、V02) と Panofsky–Wenzel、非相対論の補正、(4) 長い空洞の分割 (重み `ws`) と、設計運動量 $p_0(s)$ が変わる場合 (TRPT)、(5) 端の効果 $H_f$、を式を省かずに追う。偏向キック型の空洞 TCAVI は、ヘルプに節が無いので、ソース `ttcav.f` から書く。MULT に RF を持たせる `tmultiacc` は、この章の CAVI と同じ手続きで、MULT の本体は [08_mult.md](08_mult.md)。

**読む順**: §1 (結論) → §3 → §4 (薄い空洞の写像) → §6 (長い空洞) → §8 (TCAVI) → §9 (確認)。

## 1. 結論

- **エネルギー利得** (ヘルプ): $\Delta E = -e\,(V + V_1x + \tfrac12V_{20}x^2 + \tfrac12V_{02}y^2 + V_{11}xy)\,\sin(\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$。追跡の中では $h = E/mc^2$ の変化 `dh`= $-v\,[\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{phis}]$ ($v = V/mc^2$、$\mathtt{phis} = \omega\,t_{rf0}$、$\mathtt{dphis} = \mathtt{phis} - \mathrm{PHI}$。RAD OFF のとき)。粒子の時間は位置から $t = -z\,h_1/p_1$ ($c = 1$) で決まる。
- **写像**: 薄い空洞は、時間 $t$ を保ったまま、(i) $h_2 = h_1 + \Delta h$、(ii) $p_2 = \sqrt{h_2^2 - 1}$ ($mc$ 単位)、(iii) $z_2 = -t\,p_2/h_2$、(iv) 正準運動量 $p_x, p_y$ は変えない (追跡の内部変数 $q = p/(1+\delta)$ は $p_1/p_2$ 倍)、で、追跡と一致する (差 $2 \times 10^{-13}$、$\delta$ と $z$ は桁落ちの範囲)。**シンプレクティック** (数値微分のヤコビアンで $10^{-10}$)。
- **横方向**: $V_1, V_{20}, V_{11}, V_{02}$ は、電場の振幅の横依存。Panofsky–Wenzel の定理で、横のキックは $\Delta p_x = -\partial_xV\,\cos(\omega t - \mathtt{dphis})/(\omega p_0)$。非相対論の補正 $V_{20}, V_{02} \mathrel{+}= V(2\pi f/c)^2/(\gamma\beta)^2/4$ は、変形ベッセル関数 $I_0(\omega r/\beta\gamma c) = 1 + (\omega r/\beta\gamma c)^2/4 + \cdots$ の第 2 項 (数値確認: 2 MeV/c で ΔE の $x$ 依存が $1 + (\omega x/\beta\gamma c)^2/4$ に一致)。
- **長い空洞**: ドリフト−キック−ドリフトを、重みつきの長さ `al*ws(n)`・キックの割合 `ws(n)` で $\mathrm{ndiv}$ 回繰り返す (2 次の対称分割)。重みは $\sum\mathtt{ws} = 1$ で、設計運動量が空洞の中で増える (TRPT) とき、スライスごとに等比数列 $\mathtt{ws}(n) \propto (1+r)^{n/\mathrm{ndiv}}$ (数値確認: $e^{wl} = 1 + r$)。TRPT でない普通のリングでは $r = 0$ で一様 ($\mathtt{ws} = 1/\mathrm{ndiv}$)。
- **端**: 入口で $H_f = -(eV/L)\,[\sin(\omega t - \mathrm{dphi}) + \sin\mathrm{dphi} - \mathrm{offset}]\,(x^2 + y^2)/4\;\delta(s - s_0)$、出口は符号が逆 (ヘルプ)。Maxwell 方程式から導くと、$\beta = 1$ で同じ形の収束 (横の力 $F_r = -(e r/2)\,[\partial_sE_z + (\beta/c)\,\partial_tE_z]$) になる。K64 の式は $1/\beta$ を略している (低エネルギーで効く。解析だけで、数値の確認はしていない)。
- **TCAVI**: 偏向キック $\Delta p_x = -K_0\sin(\omega t + \phi)$ と、Panofsky–Wenzel の相手の $\Delta h = K_0\,p_0\,\omega\,x\cos(\omega t + \phi)$ の組。シンプレクティック。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran | pysad |
|---|---|---|---|
| 電圧 [V] | `VOLT` ($V$)、`DVOLT` | `vc` ($= V + \mathrm{DVOLT}$、`tturn.f` の呼び出しで足す)。$v = \texttt{vc}/\texttt{amass}\cdot|\texttt{charge}|$ | `el.VOLT`、`v` |
| 周波数 | `FREQ` [Hz]、`HARM` (FREQ = 0 のとき $\omega = \omega_0\cdot\mathrm{HARM}/c$) | `w` ($= 2\pi\,\mathrm{FREQ}/c$ [rad/m]) | `w` |
| 位相 | `PHI`、`DPHI` | `phi`、`dphi`、`phic` ($= (\texttt{phi}+\texttt{dphi})\cdot\mathrm{sign}(\texttt{charge})$) | `phic` |
| 同期の時間 | (ヘルプの $t_s$) | `trf0` (global)、`phis` $= \texttt{w}\cdot\texttt{trf0}$、`dphis` $= \texttt{phis} - \texttt{phic}$ | `trf0`、`phis`、`dphis` |
| 一定の項 | — | `offset`、`offset1` | `offset`、`offset1` |
| 横の依存 | `V1`、`V20`、`V11`、`V02` | `v1`、`v20`、`v11`、`v02` → `v1a`、`v20a`($= \tfrac12v_{20}/\texttt{amass}\ldots$)、`v11a`、`v02a` | `v1a` ほか |
| 設計運動量・エネルギー | (設計運動量 $p_0$、$h_0 = \gamma_0$) | `p0`、`h0` ($mc$ 単位) | `p0`、`h0` |
| 空洞の中の設計の変化 | (TRPT のとき $p_0(s)$) | `vnominal` $= v\sin(-\mathrm{PHI}\cdot\mathrm{sign})$、`he` $= h_0 + \texttt{vnominal}$、`pe` | `vnominal`、`he`、`pe` |
| 分割 | `L`、`FREQ` | `ndiv`、`ws(1:ndiv)`、`wsn` | `ndiv`、`ws` |
| 端 | `DISFRIN`、`FRINGE` | `fringe` (= DISFRIN が 0)、`mfring` | `fringe`、`mfring` |

SAD の座標は $(x, p_x, y, p_y, z, \delta)$ ([01 章](01_variables_and_conventions.md))。$z = -v(t - t_0(s))$、先頭が $z > 0$。追跡の内部では $p_x, p_y$ を $q = p/(1+\delta)$ の単位で持つ ([04 章](04_quad.md) 式 (4.37))。

## 3. エネルギー利得と時間

### 3.1 位置から時間へ

SAD の $z$ は $z = -v\,(t - t_0(s))$ ($v$ は粒子の速度、$t_0(s)$ は設計粒子の到着時間、[01 章](01_variables_and_conventions.md) 式 (1.7))。$c = 1$ の単位で $v = \beta = p/h$ ($p, h$ は粒子の運動量とエネルギーを $mc$、$mc^2$ で測ったもの) だから、設計粒子から見た粒子の時間のずれは

$$
t = -\frac{z}{\beta} = -z\,\frac{h_1}{p_1}\qquad(h_1 = \sqrt{1 + p_1^2},\ p_1 = p_0(1 + \delta))
\tag{10.1}
$$

(`t=-z(i)*h1/p1`、`tcav.f`)。$h_1$ は、$dv$ が小さいときは $dv$ から作る: $v/v_0 = p\,h_0/h_1$ ([01 章](01_variables_and_conventions.md) 式 (1.9)) から $h_1 = p_{1r}\,h_0/(1 - dv + \mathtt{dvfs})$ ($p_{1r} = 1 + \delta$、`dvfs` は RFSW での周波数のずらし)。$dv > 0.1$ では $h_1 = \sqrt{1 + p_1^2}$ (`p2h`)。RFSW OFF では `tlim=0` で、$t = \min(t_{\lim}, \max(-t_{\lim}, t)) = 0$ に潰れる (`tcav.f` の `tlim`)。つまり RF の位相に依存しない (§5.3)。

### 3.2 エネルギー利得

ヘルプの式を、1 つの粒子について $V(x, y) \equiv V + V_1x + \tfrac12V_{20}x^2 + \tfrac12V_{02}y^2 + V_{11}xy$ と書くと $\Delta E = -eV(x, y)\sin(\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$。$t_s$ は「放射の損失と加速のつりあいで決まる同期時間」(ヘルプ)。$\mathtt{phis} \equiv \omega t_s$ ($= \texttt{w}\cdot\texttt{trf0}$)、$\mathtt{phic} \equiv \mathrm{PHI} + \mathrm{DPHI}$ (電荷の符号つき) とおくと $\omega(t - t_s) + \mathrm{PHI} + \mathrm{DPHI} = \omega t - \mathtt{dphis}$ ($\mathtt{dphis} = \mathtt{phis} - \mathtt{phic}$)。$h = E/mc^2$、$v = V/mc^2$ ($e|\mathtt{charge}|$ を含む) で、

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

- **放射・TRPT・autophi** ($\mathtt{rad}$、$\mathtt{trpt}$、$\mathtt{autophi}$): $\mathtt{offset} = \sin\mathtt{dphis}$、$\mathtt{offset1} = 0$。(10.4) は $-v\sin(\omega t - \mathtt{dphis})$ で、ヘルプの $\Delta E$ そのもの。基準の粒子 ($t = 0$) の利得は $v\sin\mathtt{dphis}$ で、放射の損失とつりあう (つりあいの位相 $t_s$ を `trf0` が持つ)。
- **それ以外 (RAD OFF の追跡)**: $\mathtt{offset} = \sin\mathtt{dphis} - \sin\mathtt{phis}$、$\mathtt{offset1} = \sin\mathtt{phis}$。(10.4) は $\Delta h = -v\,[\sin(\omega t - \mathtt{dphis}) + \sin\mathtt{phis}]$。すなわち**基準の粒子の利得を引いた**形で、`trf0` = 0 ($\mathtt{phis} = 0$、FrequencyMap・LMA の既定) なら、ヘルプの $\Delta E$ と同じ。$\mathtt{phis} \neq 0$ では、$-v\sin\mathtt{phis}$ だけ一定の利得が加わる (放射がないので、基準の粒子の損失を考えず、つりあいの位相にすると基準の粒子が加速も減速もしない形に調整する。ソースからの読み。未確認)。

数値確認 ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (a)): 3 GeV の電子、VOLT = 2 MV、FREQ = 500 MHz、PHI = 0.3 rad で、$\mathtt{trf0} = 0$ と $0.4/\omega$ のとき、(10.1) と (10.4) から作った写像が、追跡と $2 \times 10^{-13}$ 以内で一致する ($\delta$ の変化 $-5.8 \times 10^{-4}$ が $\Delta E = -V\sin(\omega t - \mathtt{dphis})$ に合う)。

## 4. 薄い空洞の写像

(10.4) でエネルギーが $h_1 \to h_2 = h_1 + \Delta h$ に変わると、粒子の運動量・$z$・$dv$・横運動量を更新する。`tcav.f` の `do 20` (と `tcavfrin`、`ttcav.f`) の手順を式にする。

1. **運動量**: $h$ と $p$ の関係 $h^2 = 1 + p^2$ ($mc$ 単位) から $p_2^2 - p_1^2 = h_2^2 - h_1^2 = \Delta h\,(h_1 + h_2)$。$a \equiv \Delta h\,(h_1 + h_2)/p_1^2$ として

$$
\frac{p_2}{p_1} = \sqrt{1 + a},\qquad \mathtt{dpr} = \sqrt{1 + a} - 1\ (\texttt{sqrt1}(a)\text{、桁落ちなし}),\qquad p_{2r} = p_{1r}(1 + \mathtt{dpr}),\ \ \delta_2 = p_{2r} - 1
\tag{10.5}
$$

($\Delta h$ が小さいと $a$ も小さく、`sqrt1` が $a/(1 + \sqrt{1 + a})$ で桁落ちを避ける。`max(a, -1)` は $p_2^2 \ge 0$ の下限。)

2. **$dv$**: 新しい $dv = 1 - v_2/v_0$ は、01 章の式 (1.9) の形で $dv = -(1 + p_{2r})\,\delta_2/[h_2(h_2 + p_{2r}h_0)] + \mathtt{dvfs}$ (`dv(i)=-(1.d0+p2r)/h2/(h2+p2r*h0)*dp2r+dvfs`)。

3. **$z$**: 時間 $t$ は、空洞が薄い (粒子が動かない) 間は変わらない。新しい速度 $\beta_2 = p_2/h_2$ のもとで $z = -\beta t$ だから

$$
z_2 = -t\,\beta_2 = -t\,\frac{p_{2r}\,p_0}{h_2}
\tag{10.6}
$$

($p_2 = p_{2r}p_0$。`z(i)=-t*p2r*p0/h2-dzn`。`dzn` は設計の速度の変化の分、§6。)

4. **横運動量**: 縦方向の加速は、正準の $p_x, p_y$ を変えない ($H$ に $x, y$ が出ない薄い空洞の流れ。ただし V1 などがあるときは §5)。追跡の内部変数 $q = p_x/p$ は、$q_2 = q_1\,p_{1r}/p_{2r}$ (`px(i)=px(i)*p1r/p2r`)。つまり**加速で $q$ が減る (断熱減衰)**。

この写像は、(i) 運動量・エネルギーの関係が厳密、(ii) $z$ の更新が時間保存、(iii) 正準運動量が不変、なので、1 つの粒子についての写像として厳密である (空洞の長さが 0 の極限)。**シンプレクティック**であることを、6×6 のヤコビアン (中心差分、刻み $10^{-7}$) で確認した: V1・V20 なしの薄い空洞 $|J^\top SJ - S| = 1.2 \times 10^{-10}$ ($\det J - 1 = 1.2 \times 10^{-10}$)、V1・V20・V11・V02 ありで $2.1 \times 10^{-10}$ ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (b)。中心差分の精度)。

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

(単位は $v = V/mc^2$ で測る)。$\partial_xV = V_1 + V_{20}x + V_{11}y$ ($=\mathtt{v1a} + 2\,\mathtt{v20a}\,x + \mathtt{v11a}\,y$、`v20a` $= \tfrac12V_{20}\ldots$)、$\partial_yV = V_{11}x + V_{02}y$。コードは `av=-(cos(2*ph-dphis)*wi-offset1*t)*wsn/p0`、`dpx=(v1a+2*v20a*x+v11a*y)*av`、`dpy=(v11a*x+2*v02a*y)*av`、`px=(px*p1r+dpx)/p2r` ($\mathtt{wi} = 1/\omega$)。$\cos(2\,ph - \mathtt{dphis}) = \cos(\omega t - \mathtt{dphis})$ で (10.8) と一致する。さらに $\mathtt{offset1}\ne0$ の項は、(10.4) の一定の利得 $-v(x, y)\,\mathtt{offset1}$ に対応するキックのハミルトニアン $K_{\rm extra} = -V(x, y)\,\mathtt{offset1}\,t$ ($\partial_tK_{\rm extra} = -V\,\mathtt{offset1}$) から $\Delta p_x = +\partial_xV\,\mathtt{offset1}\,t/p_0$ ($= -\partial_xK_{\rm extra}$)、コードの $-(\ldots - \mathtt{offset1}\,t)$ の第 2 項と同じ。つまり、$\mathtt{v1a}$ などのある空洞も、(10.7) と一定の項から導かれる**1 つの時間依存のハミルトニアンの流れ**に当たり、ヤコビアンの数値確認 (§4) と整合する。

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

**刻みの数**: `tcav.f`: $\mathrm{ndiv} = 1 + \mathrm{int}\bigl(\min(|\omega L|,\ \sqrt{r^2/(3\varepsilon)})\bigr)$、$\varepsilon = 10^{-3}$、$r = |\mathtt{vnominal}\,(1/h_0 + 1/h_e)|$。$r = 0$ (TRPT でない普通のリング) のとき $\mathrm{ndiv} = 1$、キックは 1 回 (`ws = 1`)。$r > 10^{-6}$ のとき、$\mathrm{ndiv} = \min(1000, \max(\mathrm{ndiv}, \mathrm{nint}(10\,wl)))$。CALC (`tcave.f`) は $\varepsilon = 10^{-2}$ で、$\mathrm{ndiv} = 1 + \mathrm{int}\bigl(\min(|\omega L|,\ \sqrt{[v(1/h_0 + 1/h_1)]^2/(3\varepsilon)})\bigr)$ ($v$ は電圧そのもの。TRPT の `vnominal` ではない)。式の根拠 ($\sqrt{r^2/3\varepsilon}$ の 3) は**ソースに書かれていない**。

### 6.2 設計運動量が増える (TRPT) とき

`vnominal` $= v\sin(-\mathrm{PHI}\cdot\mathrm{sign}(\mathtt{charge}))$ ($\mathtt{trpt}$ のときだけ。TRPT は「ビームラインが蓄積リングではなく輸送路で、運動量が加速で変わる」ことを宣言、ヘルプの TRPT。DPHI は設計運動量 $p_0(s)$ に入れない、ヘルプの DPHI) は、空洞を設計粒子が通ったときのエネルギー利得 $h_e - h_0 = \mathtt{vnominal}$。空洞の中のエネルギーは一様に増えるとし、各刻みの利得を `vnominal*wsn` とする (`sv=sv+vnominal*wsn`、`h2=h0+sv`)。この間の設計の運動量・速度の変化を、$dv_n$ と $dz_n$ として $z$ に反映する (`dzn=-dvn*al*wsn*0.5`)、空洞を出たあと $\delta$ を新しい設計運動量 $p_e$ に対して測り直す: $\delta' = \delta\,p_0/p_e - \mathtt{dpepe}$ ($\mathtt{dpepe} = \mathtt{vnominal}\,(h_e + h_0)/((p_e + p_0)\,p_e)$。`g(i)=dp1`)、$dv$ も作り直し、`pgev=pe*amass`、`tphyzp` で設計運動量の物理量を更新する (`tcav.f` の最後)。

**重み ws**: $r > 10^{-6}$ のとき

$$
wl = \mathrm{arcsinh}\,\frac{r(2 + r)}{2(1 + r)},\quad r_1 = \frac{wl}{2\,\mathrm{ndiv}},\quad \mathtt{ws}(1) = \frac{2e^{r_1}\sinh r_1}{r} = \frac{e^{2r_1} - 1}{r},\quad \mathtt{ws}(n) = \mathtt{ws}(n-1)\,e^{2r_1}
\tag{10.10}
$$

(`we=1+ws1`、$e^{2r_1} = 1 + 2e^{r_1}\sinh r_1$。) $\mathtt{ws}$ は公比 $e^{2r_1}$ の等比数列で、和は $\mathtt{ws}(1)\,(e^{2r_1\mathrm{ndiv}} - 1)/(e^{2r_1} - 1) = (e^{wl} - 1)/r$。さらに $x \equiv r(2 + r)/(2(1 + r))$ として $x^2 + 1 = (r^2 + 2r + 2)^2/(4(1 + r)^2)$ ($(r^2 + 2r)^2 + 4(1 + r)^2 = r^4 + 4r^3 + 8r^2 + 8r + 4 = (r^2 + 2r + 2)^2$) なので、$e^{wl} = x + \sqrt{x^2 + 1} = [(r^2 + 2r) + (r^2 + 2r + 2)]/(2(1 + r)) = 1 + r$。したがって**$\sum\mathtt{ws} = r/r = 1$、最後の刻みの重みは最初の $(1 + r)^{(\mathrm{ndiv}-1)/\mathrm{ndiv}}$ 倍** ($r = 0.5$、$\mathrm{ndiv} = 12$ で 1.4502、$r = 3$、$\mathrm{ndiv} = 20$ で 3.7321。[cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (d))。設計粒子のエネルギーが増えるにつれ、長さと電圧の割合を大きくする形で、運動量が大きいほどキックの長さが長い。**この重みをこう選ぶ理由 (式 (10.10) の $wl$) は、ソースに書かれていない** (未確認)。

### 6.3 V20 の補正と平均

`vcorr` の運動量は、入口 $p_0$ と出口 $p_e$ の $1/p$ の平均で入れる ($0.5/p_0 + 0.5/p_e$、ヘルプ「入口と出口の $1/(\beta\gamma)$ の平均の逆数」)。

## 7. 端の効果 $H_f$

ヘルプ (CAVI の VOLT の説明、原文で確認): 「CAVI は最低次の端の効果を含む。入口 $s_0$ でのハミルトニアン $H_f = -(e(\mathrm{VOLT}+\mathrm{DVOLT})/L)\,(\sin(\omega t - \mathrm{dphi}) + \sin(\mathrm{dphi}) - \mathrm{offset})\,(x^2 + y^2)/4\,\delta(s - s_0)$。$\mathrm{dphi}$ と offset は空洞の位相と放射の損失で決まる (NORAD のときだけ 0 でない)。符号は出口で反転する。既知の導出と整合する」(導出者の名前が書かれている。この文書では省く)。

### 7.1 Maxwell からの導出

軸の近くの軸対称の場 $E_z(r, s, t) = E_0\,g(s)\cos(\omega t + \varphi_E)$ から。$\nabla\cdot\vec E = 0$ の $r$ 成分 $\tfrac1r\partial_r(rE_r) = -\partial_sE_z$ から $E_r = -\tfrac r2\partial_sE_z$。$\nabla\times\vec B = c^{-2}\partial_t\vec E$ の $z$ 成分 $\tfrac1r\partial_r(rB_\theta) = c^{-2}\partial_tE_z$ から $B_\theta = \tfrac{r}{2c^2}\partial_tE_z$。速度 $\beta c\,\hat s$ の粒子の半径方向の力は $F_r = e(E_r - \beta cB_\theta)$ ($(\vec v\times\vec B)_r = -v_sB_\theta$) なので

$$
F_r = -\frac{e\,r}{2}\Bigl[\partial_sE_z + \frac\beta c\,\partial_tE_z\Bigr]
\tag{10.11}
$$

入口の端では $g(s)$ が 0 から 1 に立ち上がり、$\partial_sE_z = E_0\,g'(s)\cos(\omega t + \varphi_E)$ が支配的 (第 2 項は端の外側の一様な場の区間でも出るが、端の幅 → 0 の極限で $g'$ の積分 1 に比べて小さい)。$ds = \beta c\,dt$ で $\int F_r\,dt = \int F_r\,ds/(\beta c)$ なので、端を通る間の運動量の変化は

$$
\Delta p_r = -\frac{e\,E_0\,r}{2\beta c}\cos(\omega t + \varphi_E)
\tag{10.12}
$$

($\int g'\,ds = 1$。) $\Delta E = -eV\sin(\omega t - \ldots)$ の振幅は $eE_0L = eV$ で、位相は $E_0\cos(\omega t + \varphi_E) \leftrightarrow -(V/L)\sin(\omega t - \mathrm{dphi})$ の読み替えになり、$\Delta p_r = +\dfrac{e\,V\,r}{2L\,\beta c}\sin(\omega t - \mathrm{dphi})$。これは、ヘルプの $H_f$ から $\Delta p_x = -\partial_xH_f = (V/L)(\sin(\ldots) + \ldots)\,x/2$ ($e$ は $\Delta E$ の式の $e$ に入れた電荷 $\to$ $v$ で測る) と同じ形 (定数の $\sin\mathrm{dphi} - \mathrm{offset}$ の項を除く。定数の項は、(10.4) の一定の利得 $v\,[\sin\mathtt{dphis} - \mathtt{offset}]$ の横方向の相手として、(10.7) 以下と同じ議論で出る)。出口では $g' \to -g'$ で符号が反転する (ヘルプ)。

**K64 は $1/\beta$ を略す**: 追跡 `tcavfrin` は `vf=v/al/p0*.5`、`dpt=vf*(sin(ph-dphis)+s0)/p1r`、`px=px+x*dpt` で、$\Delta p_x/p_0 = x\,\dfrac{v}{2Lp_0}\bigl(\sin(\omega t - \mathtt{dphis}) + s_0\bigr)$ ($s_0 = \sin\mathtt{dphis} - \mathtt{offset}$)。(10.12) は $1/\beta$ を持つ。数値で、追跡の入口の端 (FRINGE = 1 と DISFRIN = 1 の差) が $x\,v\sin(\ldots)/(2Lp_0)$ に一致する: 3 GeV ($\beta = 1.0000$) で比 1.0000、2 MeV/c ($\beta = 0.9689$) で 0.9993、0.6 MeV/c ($\beta = 0.7613$) で 0.9734 ([cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) (f))。つまり**K64 の式に $1/\beta$ は入っていない**。(10.12) の $1/\beta$ が正しいなら、低エネルギー (0.6 MeV/c で 31%) の端の収束が足りない。(10.12) の導出は**解析だけで、ローレンツ力の数値積分では確認していない**。

### 7.2 端のエネルギー変化

`tcavfrin` は、端の収束 (上) のほかに、端の中のエネルギー変化 $\Delta h = -\omega\,p_0\,v_f\,(x^2 + y^2)/2\cdot\cos(\omega t - \mathtt{dphis})$ ($v_f = v/(2Lp_0)$) を入れ、$p$、$z$ ($z \to z\,(p_{2r}/p_{1r})\,(h_1/h_2)$) を更新する (`tcavfrin`)。これは $H_f$ が時間 $t$ に依存することの相手 ($\Delta E = \partial H_f/\partial t$) で、(10.7) と同じ構造。**符号と係数の導出は、この章では行っていない** (ヘルプの $H_f$ から $\partial_tH_f$ を取ると $-\omega\,v\,(x^2 + y^2)\cos(\ldots)/(4L)$、符号は SAD の $\Delta E = \partial K/\partial t$ (§5.1) で決まる)。出口は $v \to -v$ (`tcavfrin(…, -v, …)`)。入口は `mfring ≥ 0` かつ `mfring ≠ 2` (`FRINGE` = 0, 1, 3)、出口は `mfring ≥ 0` かつ `mfring ≠ 1`。`fringe` は DISFRIN が 0 のとき。**ヘルプの FRINGE・DISFRIN の表** (DISFRIN ≠ 0 は端なし。FRINGE = 1 は入口、2 は出口、3 は両方、0 も両方): 原文で確認。K64 の `mfring` の変換 (`tffs.f` の `frmd`) は、$\mathrm{FRINGE}$ を 0, 1, 2, 3 から内部値に直す (FRINGE = 0 で両方になる表と整合するかは、ここでは確認していない)。

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
| `cav_01_energy_map.py` | (a) 薄い CAVI の写像 (10.1)(10.4)〜(10.6) と追跡の一致、(b) シンプレクティック性 (V1・V20・フリンジつき)、(c) 非相対論の補正と $I_0$、(d) 重み `ws` の和と $e^{wl} = 1 + r$、(e) TCAVI の式・シンプレクティック性・$K$ の流れ、(f) 端の収束の $1/\beta$ | [cav_01_energy_map.out.txt](outputs/cav_01_energy_map.out.txt) |

実行: リポジトリのルートから `py -3.12 physics_guide/scripts/cav_01_energy_map.py` (数秒。numba のキャッシュは環境変数 `NUMBA_CACHE_DIR` で指定できる)。

## 10. 未確認・限界・使った条件

**使った条件**: 電子、運動量 3 GeV/$c$ (§3〜§5、§8)、低エネルギーの確認は 0.6〜2 MeV/c (§5.2、§7)。VOLT = 2 MV (§3、§4)、FREQ = 500 MHz、PHI = 0.3 rad、$\delta = O(10^{-3})$、$z = O(0.1)$ m。放射・スピン・ウェイク (`txwake`)・TRPT (`vnominal ≠ 0`) の追跡は確認していない。

**この章の範囲外**: ウェイク場 (`txwake`、WAKE の章)、放射 (RADCOD、`tradk`、放射による損失と `trf0` の決め方)、`tmultiacc` (MULT + RF)、CALC (`tcave`) の行列の導出と追跡との比較、COORD・ずれ (DX・DY・ROTATE) の `TENT`/`TEXIT`。

**未確認・限界**:
- **TRPT (`vnominal ≠ 0`) の写像**は、読みだけ。$\mathtt{ws}$ の和と $e^{wl} = 1 + r$ は数値で確認したが、(10.10) の重みの選び方の理由、`ndiv` の式 ($\sqrt{r^2/3\varepsilon}$ の 3)、$dz_n$・$dp_{epe}$ の式の導出は行っていない。
- **端の収束の $1/\beta$** (§7.1): 解析の導出だけで、数値積分は未確認。K64 の式は $1/\beta$ を持たない (確認済み)。
- **端のエネルギー変化の符号と係数** (§7.2) は導出していない。
- **`offset` の意味** (§3.3): ソースの読み。$\mathtt{phis} \ne 0$ のときの一定のエネルギー利得の物理的な意味は確認していない。
- **V1・V11・V02 は「Tracking only」** (ヘルプ) と書かれるが、`tcave.f` (CALC) は `v10`・`v20`・`v11`・`v02` を使う。CALC の V1 の行列の係数は確認していない。
- **ヘルプの V20・V02 の説明**: 「The non-relativistic corrections ... (VOLT+DVOLT)(2 Pi FREQ/c)^2/(gamma beta)^2/4 are automatically added to V20 and V02, respectively」は、足される量が V20 と V02 の係数 ($V_{20}x^2/2$ の $1/2$ を含む) の 2 倍の大きさ、すなわち $\Delta V = V(\omega/\beta\gamma c)^2(x^2 + y^2)/4$ を意味する。コードでは `v20a` ($= V_{20}/2$ の係数) に足す。文面は読み方で 2 通りあり、コードと数値 (§5.2) から後者。

## 付録 A. ヘルプの記述・K64 との違い (CAVI・TCAVI)

原文 (`SADHelp.html` の `CAVI`、`RFSW`、`RADCOD`、`TRPT`、キーワード一覧) で確認した。

| 箇所 | 内容 |
|---|---|
| CAVI の $\Delta E$ の式 | **原文で確認**: $\Delta E = -e(\mathrm{VOLT} + \mathrm{DVOLT} + V_1x + \tfrac12(V_{20}x^2 + V_{02}y^2) + V_{11}xy)\sin(2\pi\,\mathrm{FREQ}(t - t_s) + \mathrm{PHI} + \mathrm{DPHI})$ (PHI・DPHI・FREQ・VOLT の各項で同じ式)。§3 のとおり (10.4) で実装される。「$t_s$ は放射と加速のつりあいで決まる同期時間」。 |
| CAVI の VOLT の説明: 非相対論の補正 | **原文で確認、解釈を確定**: §5.2 (数値で一致)。 |
| CAVI の端の効果 $H_f$ | **原文で確認**: $H_f = -(e(\mathrm{VOLT} + \mathrm{DVOLT})/L)(\sin(\omega t - \mathrm{dphi}) + \sin\mathrm{dphi} - \mathrm{offset})(x^2 + y^2)/4\,\delta(s - s_0)$。**ヘルプは $\sin$ の中を $\omega t - \mathrm{dphi}$ と書き、DPHI の説明は $t_s$ の項を持つ形 ($\omega(t - t_s) + \ldots$)**。コードの $\sin(\omega t - \mathtt{dphis})$ ($\mathtt{dphis} = \mathtt{phis} - \mathtt{phic}$) と記号の定義が違う (ヘルプの dphi は `dphis` に当たると読める。確認していない)。Maxwell からの導出 (§7.1) と形は同じで、$1/\beta$ の因子が K64 に無い (§7.1)。 |
| CAVI の DISFRIN・FRINGE の表 | **原文で確認**: DISFRIN ≠ 0 は端なし。FRINGE = 1 入口、2 出口、3 入口と出口、0 は両方。他の要素 (QUAD) の表と違い、CAVI の表は「Linear」の列が無い。 |
| CAVI の DPHI の説明 | **原文で確認**: 「DPHI is not taken into account to determine the design momentum $p_0(s)$」。コードで `vnominal` は PHI だけを使う (`tffs.f` の `vnominal=…sin(-cmp%value(ky_PHI_CAVI)…)`)。一致。 |
| RADCOD・TRPT・RFSW | **原文で確認** (§5.3、§6.2)。TRPT の既定の運動量分布は「一様」(Emittance[] での扱い) は未確認。 |
| TCAVI | **原文に節が無い** (キーワード一覧と RFSW の言及のみ)。式はすべて `ttcav.f`・`ttcave.f` から。 |
| ヘルプに無いもの | (1) 薄い空洞の運動量・$z$・横運動量の更新 (§4)。(2) Panofsky–Wenzel の関係 (§5.1)。(3) 長い空洞の分割 (§6)、重み。(4) 端の導出 (§7)。(5) TCAVI の式。 |

**K64 の実装の癖・疑い**: 端の収束に $1/\beta$ が無い (§7.1)。低エネルギーの空洞では、ヘルプの「lowest order」の範囲を超える。

## 付録 B. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `v` | $V/(mc^2)\cdot|\mathtt{charge}|$ (電圧を静止エネルギーで割った量) | `tcav.f`: `v=vc/amass*abs(charge)` | `v` |
| `w`・`wi` | $\omega = 2\pi f/c$ [rad/m]、$1/\omega$ | `tcav.f`: `w`、`wi=1.d0/w` | `w`・`wi` |
| `phic`・`phis`・`dphis` | $(\mathrm{PHI} + \mathrm{DPHI})\,\mathrm{sign}(\mathtt{charge})$、$\omega\,t_{rf0}$、$\mathtt{phis} - \mathtt{phic}$ | `tcav.f` | `phic`・`phis`・`dphis` |
| `trf0` | 同期時間 $t_s$ (global。RF の位相の基準) | `tmacro` | `trf0` |
| `offset`・`offset1` | 一定の利得の調整 ($\sin\mathtt{dphis}$ または $\sin\mathtt{dphis} - \sin\mathtt{phis}$、$0$ または $\sin\mathtt{phis}$) | `tcav.f` | `offset`・`offset1` |
| `h1`・`p1`・`p1r` | 入口のエネルギー ($mc^2$)、運動量 ($mc$)、$1 + \delta$ | `tcav.f` | `h1`・`p1`・`p1r` |
| `dh`・`h2`・`a`・`dpr`・`p2r` | $\Delta h$、$h_1 + \Delta h$、$\Delta h(h_1 + h_2)/p_1^2$、$\sqrt{1 + a} - 1$、$p_{1r}(1 + \mathtt{dpr})$ | `tcav.f` | 同名 |
| `t`・`ph`・`tlim` | $t = -z\,h_1/p_1$、$\omega t/2$、$|t|$ の上限 (RFSW OFF で 0) | `tcav.f` | 同名 |
| `v1a`・`v20a`・`v11a`・`v02a` | $V_1$、$V_{20}/2$、$V_{11}$、$V_{02}/2$ ($/\mathtt{amass}$) | `tcav.f` | 同名 |
| `vcorr` | $v\,(\omega\,(0.5/p_0 + 0.5/p_e))^2/4$ (非相対論の補正) | `tcav.f` | `vcorr` |
| `av`・`dpx`・`dpy` | $-\cos(2\,ph - \mathtt{dphis})/\omega\,\ldots/p_0$、横のキック | `tcav.f` | 同名 |
| `vnominal`・`he`・`pe`・`dpepe` | 設計の利得 $v\sin(-\mathrm{PHI})$、$h_0 + \mathtt{vnominal}$、$h_e$ の運動量、再基準化 | `tcav.f` | 同名 |
| `ndiv`・`ws`・`wsn`・`wl`・`r` | 刻みの数、重み、その刻みの重み、$wl$、$r = |\mathtt{vnominal}(1/h_0 + 1/h_e)|$ | `tcav.f` | 同名 |
| `sv`・`dzn`・`dvn` | 設計エネルギーの累積利得、設計の $z$ のずれ、設計の $dv$ | `tcav.f` | 同名 |
| `vf`・`s0`・`dpt` | 端の係数 $v/(2Lp_0)$、$\sin\mathtt{dphis} - \mathtt{offset}$、端の収束 | `tcavfrin` (`tcav.f`) | `vf`・`s0` |
