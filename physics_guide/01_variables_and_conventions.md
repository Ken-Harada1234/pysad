# 01. 変数と約束 (四極の章に必要な範囲)

**この章で分かること**: SAD の正準変数 $(x, p_x, y, p_y, z, \delta)$ が何で、ヘルプのラグランジアンからどう導かれるか。追跡の内部だけ運動量を $p_x/(1+\delta)$ に取り替えていること。$z$ の符号と、速度の差 $dv$。単位と K1 の約束。

**前提**: 古典力学 (ラグランジアン・ハミルトニアン・母関数) と、加速器の座標系 $(x, x', y, y', z, \delta)$。SAD の内部は知らなくてよい。

**読む順**: 先に §1 (結論) と §2 (記号の表) を見て、四極の章 ([04_quad.md](04_quad.md)) に進んでよい。導出 (§3〜§5) は、四極の章が「$p_x$ とは何か」「$z$ の更新の $dv$ はどこから来るか」と問うときに戻って読む。

この章は、ガイド全体の共通の前提のうち、四極の章が使う分だけを書く。残り (曲線座標の曲率 $1/\rho$ が効く導出、ソレノイドのベクトルポテンシャル、加速のある場合) は、各要素の章で足す。

数式の番号 $(1.n)$ は、この章の中と、他の章からの参照に使う。

## 1. 結論

- SAD は、基準線に沿った距離 $s$ を独立変数にし、正準変数 $(x, p_x, y, p_y, z, \delta)$ で粒子の運動を書く。$p_x, p_y$ は設計運動量 $p_0$ で規格化した**正準運動量**、$\delta = p - 1$ ($p = P/p_0$ は全運動量の規格化値)、$z = -v\,(t - t_0(s))$ はバンチ内の位置で、**先頭が $z > 0$**。
- ハミルトニアンは次の形 (ヘルプの dynamics > Hamiltonian。誤植を直した形):

$$
H = -\left(\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2} + a_s\right)\left(1 + \frac{x}{\rho}\right) + \frac{E}{p_0 v_0}
\tag{1.1}
$$

  ここで $a_x, a_y, a_s$ は規格化したベクトルポテンシャル $e\vec{A}/p_0$、$E$ は粒子のエネルギー、$v_0$ は設計速度。最後の項 $E/(p_0 v_0)$ が、$z$ の運動 (経路の遅れ) を作る。
- **追跡の内部変数は正準運動量そのものではない**: 追跡 (`tquad` など) は、座標の入口で $q_x = p_x/(1+\delta)$、$q_y = p_y/(1+\delta)$ に取り替えて計算し、出口で戻す (`tconvm`)。$q$ は粒子の進行方向と $s$ 軸のなす角に近い量 (傾き $x' = q_x/\sqrt{1 - q_x^2 - q_y^2}$)。転送行列 (CALC・Emittance[]) は、正準運動量のまま作る。
- 速度の差 $dv = 1 - v/v_0$ は、$\delta$ の関数として閉じた式 (式 (1.9)) で決まり、$z$ の自由な更新は $z' = -dv$。電子 3 GeV ($\gamma_0 \approx 5871$) なら $dv \approx -\delta/\gamma_0^2 \approx -2.9 \times 10^{-8}\,\delta$ で、$\delta = 10^{-3}$ のとき 1 m 進んで $z$ が $2.9 \times 10^{-11}$ m ずれる。
- K1 は**積分強度** ($K1 = B'L/B\rho$、単位 1/m)。K1 > 0 は水平収束 (§6)。

## 2. 記号と単位

| 量 | SAD (台本・ヘルプ) | K64 の Fortran の変数 | pysad |
|---|---|---|---|
| 水平・垂直の位置 $x, y$ [m] | 座標 `{x, y}` | `x(i)`, `y(i)` (追跡)、`cod(1)`, `cod(3)` (CALC) | `c[0]`, `c[2]` |
| 正準運動量 $p_x, p_y$ (無次元。$p_0$ で規格化) | `{px, py}` | `cod(2)`, `cod(4)` (CALC)。追跡の `px(i)`, `py(i)` は $p_x/(1+\delta)$ | 外の API は正準。`track6._Particle` の内部 `c[1]`, `c[3]` は $p_x/(1+\delta)$ |
| 経路の遅れ $z$ [m] (先頭が正) | `z` | `z(i)`, `cod(5)` | `c[4]` |
| 運動量偏差 $\delta = p - 1$ | `dp` | `g(i)`, `gp(i)` (追跡)、`cod(6)` (CALC) | `c[5]` |
| 全運動量の規格化値 $p = 1 + \delta$ | | `pr`, `p` | `pr` |
| 設計運動量 ($mc$ 単位 = $\beta_0\gamma_0$)・設計のローレンツ因子 $\gamma_0$ | `MOMENTUM` から | `p0`, `h0` (`tmacro`) | `emit6d._Phys.p0`, `.h0` |
| 速度の差 $dv = 1 - v/v_0$ | | `dv(i)` (追跡)、`dv` (`tgetdv`) | `_Particle.dv`、`_Phys.dv()` |
| K1 (積分強度) [1/m] | `K1` | `ak` (四極の `tsolqu` など)。長さあたりは `akk = ak/al` | 台本の `K1` は積分強度。`Quadrupole.K1` は長さあたり [1/m²] (パーサが `K1/L` に直す) |

- **添字**: Fortran は 1 始まり (`cod(1)` が $x$)、pysad は 0 始まり (`c[0]` が $x$)。この表の行は、変数の順 $(x, p_x, y, p_y, z, \delta)$ に対応する。
- **単位**: 長さは m。$p_x, p_y, \delta, q_x, q_y$ は無次元。ビームの運動量・エネルギーは台本の `MOMENTUM` で与える (eV)。

## 3. ラグランジアンからハミルトニアンまで

ヘルプの dynamics > Lagrangean・Hamiltonian の途中を、省かずに追う。(解析力学の基礎からの別の導き方 (ローレンツ力の一般化ポテンシャル、$t \to s$ の入れ替え、展開の途中式) は [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md)。)

**位置と時間の座標**。基準線に沿った距離 $s$ と、法線・従法線方向の変位 $(x, y)$ で位置を決める ($\vec{n}, \vec{b}, \vec{t}$ は右手系)。曲率 $1/\rho$ は水平面だけ (ヘルプの independent-variable)。微小な変位は

$$
d\vec{r} = \left(1 + \frac{x}{\rho}\right) ds\,\vec{t} + dx\,\vec{n} + dy\,\vec{b}
$$

粒子 (質量 $m$、電荷 $e$) の作用は、時間 $t$ を独立変数にすれば、規格化したベクトルポテンシャル $\vec{a} = e\vec{A}/p_0$ を使って

$$
S = \int L_t\, c\,dt, \qquad
L_t = -\frac{mc}{p_0}\sqrt{1 - \dot x^2 - \dot y^2 - \left(1 + \frac{x}{\rho}\right)^2 \dot s^2}
 + a_x \dot x + a_y \dot y + \left(1 + \frac{x}{\rho}\right) a_s \dot s
\tag{1.2}
$$

($\dot{\ } = d/d(ct)$)。**$s$ を独立変数に取り替える**には $L = L_t\, d(ct)/ds$ とする ($ct' \equiv d(ct)/ds$):

$$
L = -\frac{mc}{p_0}\sqrt{(ct')^2 - x'^2 - y'^2 - \left(1 + \frac{x}{\rho}\right)^2}
 + a_x x' + a_y y' + \left(1 + \frac{x}{\rho}\right) a_s
\tag{1.3}
$$

(式 (1.2)(1.3) の根号の中の符号は、ヘルプの印刷では $+\dot y^2$ などになっている誤植。次の正準運動量の式の分母の符号が正しい。付録 A)。

**正準運動量**。根号を $\Pi$ と書く ($\Pi = \sqrt{(ct')^2 - x'^2 - y'^2 - (1 + x/\rho)^2}$、$\hat m = mc/p_0$):

$$
p_x = \frac{\partial L}{\partial x'} = \frac{\hat m x'}{\Pi} + a_x,\quad
p_y = \frac{\hat m y'}{\Pi} + a_y,\quad
p_t = \frac{\partial L}{\partial t'} = -\frac{\hat m c^2 t'}{\Pi}
\tag{1.4}
$$

**ハミルトニアン** $H_t = x'p_x + y'p_y + t'p_t - L$。(1.4) を代入すると、$x'(p_x - a_x) + y'(p_y - a_y) + t'p_t = -\hat m\,[\Pi^2 + (1 + x/\rho)^2]/\Pi$ なので

$$
H_t = -\frac{\hat m (1 + x/\rho)^2}{\Pi} - \left(1 + \frac{x}{\rho}\right) a_s
$$

一方、(1.4) から $(p_t/c)^2 - (p_x - a_x)^2 - (p_y - a_y)^2 - \hat m^2 = \hat m^2 (1 + x/\rho)^2/\Pi^2$ なので、$\hat m(1 + x/\rho)/\Pi$ は左辺の平方根に等しい。結局

$$
H_t = -\left(\sqrt{\frac{p_t^2}{c^2} - \hat m^2 - (p_x - a_x)^2 - (p_y - a_y)^2} + a_s\right)\left(1 + \frac{x}{\rho}\right)
\tag{1.5}
$$

(ヘルプの式は最後の $(p_y - a_y)^2$ の符号が誤植。付録 A)。物理的な枝では $p_t = -E/p_0$ ($E$ は粒子のエネルギー) で、根号の中の $p_t^2/c^2 - \hat m^2$ は $p^2$ ($p = P/p_0$) に等しい ($E^2 - m^2c^4 = P^2c^2$ から)。つまり (1.5) は「$H_t = -$(縦方向の運動量) $\times (1 + x/\rho)$」の形で、これは曲線座標で $s$ を時間に取ったときの標準的な形。

**$(t, p_t)$ から $(z, p)$ への正準変換**。$t$ と $p_t$ は低エネルギー ($\gamma \sim 1$) で扱いにくいので、SAD は $(z, p)$ に変える。母関数を (古い運動量 $p_t$、新しい座標 $z$) の関数として

$$
G(p_t, z; s) = z\,p(p_t) + p_t\,t_0(s),
\qquad p(p_t) = \frac{1}{c}\sqrt{p_t^2 - \frac{m^2 c^4}{p_0^2}}
\tag{1.6}
$$

とおく ($t_0(s)$ は設計粒子が位置 $s$ に着く時刻)。ヘルプの規約では、変換の関係式は

$$
p = \frac{\partial G}{\partial z},\qquad t = \frac{\partial G}{\partial p_t},\qquad H = H_t - \frac{\partial G}{\partial s}
\tag{1.7}
$$

(通常の型 3 の母関数 $F_3 = -G$ の関係式 $q = -\partial F_3/\partial p_{\rm old}$、$P = -\partial F_3/\partial z$、$K = H + \partial F_3/\partial s$ と同じ)。計算すると:

- $p = \partial G/\partial z = p(p_t)$: 全運動量の規格化値。
- $t = z\,\partial p/\partial p_t + t_0$。$\partial p/\partial p_t = p_t/(c^2 p)$ で、$p_t = -E/p_0$ なので $= -E/(p_0 p c^2) = -1/v$ ($v = P c^2/E = p_0 p c^2/E$ は粒子の速さ)。よって $t = -z/v + t_0$、すなわち **$z = -v\,(t - t_0(s))$**。遅れて着く粒子 ($t > t_0$) が $z < 0$、先頭は $z > 0$。
- $\partial G/\partial s = p_t\,t_0'(s) = p_t/v_0 = -E/(p_0 v_0)$ ($t_0' = 1/v_0$)。よって

$$
H = H_t + \frac{E}{p_0 v_0} \quad\Rightarrow\quad (1.1)
$$

ヘルプの $G$ は、$t_0(s)$ の項が $-t_0(s)$ と印刷され、$p, t$ の式の $c$ の置き方もずれている (自然単位 $c = 1$ で読めば一致)。上の (1.6) は、(1.7) から $z = -v(t - t_0)$ と式 (1.1) の $H$ が出るように直した形 (付録 A)。

## 4. $z$ の運動と $dv$

(1.1) から、$z$ の運動方程式 ($z$ と $p$ は正準共役) は

$$
z' = \frac{\partial H}{\partial p} = -\frac{p}{\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2}}\left(1 + \frac{x}{\rho}\right) + \frac{v}{v_0}
\tag{1.8}
$$

(第 2 項は $\partial E/\partial p = p_0 v$ から)。磁場のない直線上 ($\rho = \infty$、$a = 0$) の軸上の粒子 ($p_x = p_y = 0$) では $z' = -1 + v/v_0 = (v - v_0)/v_0$ で、速い粒子が前に出る、という常識と合う。この差を

$$
dv \equiv 1 - \frac{v}{v_0}
$$

と書く。K64 の `dv` はこれ。$\delta$ から閉じた式で作る: 設計運動量を $mc$ 単位で $p_0 \,(= \beta_0\gamma_0)$、$h_0 = \gamma_0 = \sqrt{1 + p_0^2}$、$h_1 = \sqrt{1 + p_0^2 (1 + \delta)^2}$ とすると $v/v_0 = \beta/\beta_0 = (1 + \delta)\,h_0/h_1$ なので

$$
dv = 1 - \frac{(1+\delta)\,h_0}{h_1}
   = -\frac{\delta\,(2 + \delta)}{h_1\,\left(h_1 + (1 + \delta)\,h_0\right)}
\tag{1.9}
$$

(2 つ目の等号は $h_1^2 - (1 + \delta)^2 h_0^2 = 1 - (1 + \delta)^2 = -\delta(2 + \delta)$ を使った。差の形にして、$\delta \to 0$ の桁落ちを避けている。)

- `tconvm.f:18` (追跡の入口) の `dv=-g*(1.d0+pr)/h1/(h1+pr*h0)+dvfs` と `tturne.f:879` の `tgetdv` (CALC) は、この式 ($g = \delta$、$pr = 1 + \delta$) と $dv_{\rm fs}$ (RF 周波数のずらし。RFSW ON のときだけ入る) の和。
- 小さい $\delta$ で $dv \approx -\delta/\gamma_0^2$。したがって $\partial z/\partial\delta = +L/\gamma_0^2$ (長さ $L$ の自由な区間)。K64 の行列の `trans(5,6)` の `dvdp*aln` は、$dv$ の $\delta$ 微分 $= h_0/h_1^3$ (これは $\delta = 0$ で $1/\gamma_0^2$) に $L$ をかけたもの。
- $z' = -dv$ は、$z$ を「基準粒子との経路差」と見るとき、速度差だけから来る部分。傾きを持つ粒子の経路の伸び ($-(1/p_z - 1)$ に当たる項) は、各要素の中で別に数える (四極なら ΔH の流れ。[04_quad.md](04_quad.md) §5)。

## 5. 追跡の内部変数 $q = p_x/(1+\delta)$

ヘルプは $p_x$ を正準運動量として書くが、K64 の追跡は、座標の入口 (`tconvm`、`idir = -1`) で `px = px/pr`、`py = py/pr` ($pr = 1 + \delta$) に取り替え、出口 (`idir >= 0`) で `px = px + g*px` に戻す。内部の $q_x$ を使うと、式が次のように簡単になる。

- 直線上の傾き: $x' = p_x/p_z = q_x/\sqrt{1 - q_x^2 - q_y^2}$ ($p_z = \sqrt{p^2 - p_x^2 - p_y^2} = p\sqrt{1 - q_x^2 - q_y^2}$)。$\delta$ が式から消え、根号が 1 つの量 $q_x^2 + q_y^2$ だけで書ける。
- 四極の集束は、規格化した運動量で $q_x' = -(K/p)\,x$ ($K$ は長さあたりの強さ)。$\delta$ は $K/p$ の 1 か所だけに入る。

行列 (CALC・Emittance[]) は、解析的な微分の連鎖 (正準運動量のヤコビアン) を作る都合で、正準運動量のまま作る。この 2 つの変数の間の変換 (ヤコビアンの $\delta$ の列の補正) は、[04_quad.md](04_quad.md) 式 (4.37) に書く。

## 6. K1 の約束 (四極に必要な分)

- SAD の `K1` は積分強度 $K1 = B' L/B\rho$ [1/m] (`QUAD K1=...` に書く値)。**K1 > 0 は水平収束**、K1 < 0 は垂直収束 (ヘルプの QUAD > K1)。長さあたりの強さ (勾配) は $K = K1/L$ [1/m²]。Fortran の `tquad`・`tsolque` の引数 `ak` は積分強度、`akk = ak/al` が長さあたり (`tsolque.f:390`)。
- pysad のパーサは、台本の積分強度 K1 を要素クラスの `K1` (長さあたり $K1/L$) に直す。pysad の要素クラスの `transfer_matrix()` は、その長さあたりの値を使う別の系統 (線形・設計軌道のまわり) で、K64 の `tquade` とは別実装 ([04_quad.md](04_quad.md) §12)。

## 7. 未確認・限界

- この章は、曲率 $1/\rho$ が $0$ の場合の導出に絞った。$1/\rho \neq 0$ (偏向磁石) の $H$ の展開・ソレノイドの $a_x, a_y$ の取り扱いは、各要素の章で書く。
- 加速のある要素 (設計運動量 $p_0$ が $s$ とともに変わる場合) は、(1.7) の $t_0(s)$ と $E/(p_0 v_0)$ の取り扱いを変える。ヘルプも「特別な扱い」としか書いていない。この章では扱わない。
- SAD のヘルプの式と実装の関係 (「ヘルプの式は実装そのままではない」) は、各要素の章の「Fortran の対応」に書く。

## 8. 記号の凡例 (SAD の変数名の注釈の規約と、共通の名前の表)

物理ガイドは、SAD の Fortran の変数名をそのまま使う (`xs1`、`dc1`、`w1` など)。数学の標準の記号ではないので、次の規約で注釈する。

1. **書き方**: SAD・Fortran・pysad の変数名は、コード体 (`xs1`) で書く。数学の記号 ($\varphi - \sin\varphi$) とは区別する。
2. **初出で説明する**: 初めて出てくる場所で、「`xs1` ($= \varphi - \sin\varphi$、Fortran の `xsincos` が返す)」のように、(i) 中身の式、(ii) 定義の場所 (ファイル:行、サブルーチン)、(iii) なぜその形で持つか (例: $\varphi$ が小さいときの桁落ちの回避)、を書く。
3. **以後も都度注釈する**: 式の直後や表の脚注で、「`xs1` ($= \varphi - \sin\varphi$)」のように、その式に出る名前を毎回書く。
4. **章末の表**: 各章の末尾に「この章で使う Fortran の変数名」の表 (名前、式、定義の場所、pysad での名前) を置く。四極の章は [04_quad.md](04_quad.md) §16。
5. **確認の基準**: 定義は実際のコード (SAD の Fortran と pysad) で確認する。確認できない名前は「未確認」と書く。

この章で出てくる共通の名前:

| Fortran の名前 | 式・意味 | 定義の場所 | pysad での名前 |
|---|---|---|---|
| `pr` | $p = 1 + \delta$ (全運動量の規格化値) | `tconvm.f` (`pr = 1.d0 + g`)、`tgetdv` (`tturne.f:869`、`pr = 1.d0 + dp`) | `pr` |
| `g(i)`、`dp` | $\delta$ (運動量の偏差)。追跡の配列 `g(i)`、`tgetdv` の引数 `dp` | `tconvm.f`、`tturne.f:869` | `c[5]`・`dp` |
| `p0` | 設計運動量 ($mc$ 単位。$\beta_0\gamma_0$) | `tmacro` (モジュール) | `emit6d._Phys.p0` |
| `h0` | $\gamma_0 = \sqrt{1 + p_0^2}$ | `tmacro` | `emit6d._Phys.h0` |
| `p1`、`h1` | $p_1 = p_0\,p$、$h_1 = \sqrt{1 + p_1^2}$ (= `p2h(p1)`。`p2h` は `tfloor.f:463`。$p_1 > 10^3$ では級数) | `tconvm.f`、`tturne.f:869` | (式の内部) |
| `dv` | $dv = 1 - v/v_0$ (式 (1.9))。`tconvm.f` では `-g*(1.d0+pr)/h1/(h1+pr*h0) + dvfs`、`tgetdv` では `-(pr+1.d0)/h1/(h1+pr*h0)*dp + dvfs` | `tconvm.f`、`tturne.f:879` | `_Particle.dv`、`_Phys.dv()` |
| `dvfs` | RF 周波数のずらしによる $dv$ への足し込み (RFSW ON のときだけ) | `tmacro` (定義は未確認。使用は `tconvm.f`、`tturne.f:879`) | 未確認 |
| `dvdp` | $\partial dv/\partial\delta$ = `h0/h1**3` (`tturne.f:880`)。$\delta = 0$ で $1/\gamma_0^2$ | `tgetdv` | `_dvz` (`elkern.py:112`、`emit6d.py:427`) |
| `cod(1..6)` | CALC の閉軌道 ($x, p_x, y, p_y, z, \delta$) の配列 | `tturne.f` など | `c[0..5]` (0 始まり) |
| `akk` | $K1/L$ (長さあたりの四極強度、K1 は積分強度 `ak`) | `tsolqu.f:52` (`akk = ak/al`) | `akk` |

数学関数 (`xsincos`・`xsincosh`・`sqrt1`・`p2h`) は `tfloor.f` の `mathfun` モジュール。小さい引数での桁落ちを避けるための関数で、定義は [04_quad.md](04_quad.md) §16 の表。

## 付録 A. ヘルプの記述との違い (この章に関わる分)

ヘルプ (SADHelp、dynamics の節) の式を実際に導いてみて見つけた、記述の誤りと省略。内容は上の本文で正しい形に直してある。

| 箇所 | ヘルプの記述 | 正しい形 (本文の式) |
|---|---|---|
| Lagrangean の $L_t$ (式 1.2) の根号 | $\sqrt{1 - \dot x^2 + \dot y^2 + (1 + x/\rho)^2\dot s^2}$ | $\sqrt{1 - \dot x^2 - \dot y^2 - (1 + x/\rho)^2\dot s^2}$ |
| 同じ節の $L$ ($s$ を独立変数にしたもの、式 1.3) の根号 | $\sqrt{c^2t'^2 - x'^2 + y'^2 + (1 + x/\rho)^2}$ | $\sqrt{c^2t'^2 - x'^2 - y'^2 - (1 + x/\rho)^2}$ (次の $p_x, p_y$ の式の分母は正しい形) |
| Hamiltonian の $H_t$ (式 1.5) | 根号の中に $+(p_y - a_y)^2$ | $-(p_y - a_y)^2$ |
| 母関数 $G$ (式 1.6) | $G = \frac{z}{c}\sqrt{\cdots} - t_0(s)$ (これでは $t = \partial G/\partial p_t$ に $t_0$ が出ない) | $G = z\,p(p_t) + p_t\,t_0(s)$ |
| $p, t$ の式 | $p = \sqrt{p_t^2 p_0^2 - m^2c^4}/p_0$ ($c$ が足りない)、$t = -z\sqrt{p^2p_0^2 - m^2c^2}/(cpp_0) + t_0$ (根号の中が $-m^2c^2$。$E$ を作るなら符号は + ) | 式 (1.6)(1.7) とその下の計算 |
| $E$ の定義 | $E = \sqrt{m^2c^4 + p_0^2 p^2}$ | $E = \sqrt{m^2c^4 + p_0^2p^2c^2}$ ($c = 1$ なら同じ) |

省略されているもの: (1) 導出の途中全部 (本章の §3)。(2) 追跡の内部変数 $q = p_x/(1+\delta)$ (§5)。(3) $dv$ の式 (§4)。(4) K1 と $a_s$ の関係 ([04_quad.md](04_quad.md) §3)。
