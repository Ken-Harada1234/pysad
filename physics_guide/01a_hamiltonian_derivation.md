# 01a. ハミルトニアンから運動方程式まで (ラグランジアン、正準変換、$t \to s$ の入れ替え)

**この章で分かること**: 荷電粒子のハミルトニアンが、(1) 速度に依存するポテンシャルをもつラグランジアンから、(2) 曲線座標 (Frenet–Serret 座標) への正準変換と、(3) 独立変数 $t \to s$ の入れ替えを経て、(4) 加速器でおなじみの 2 次の運動方程式 (偏向磁石の $x'' = -x/\rho^2 + \delta/\rho$、四極の $x'' = -Kx$) になるまでを、途中の式を省かずに追う。とくに次の 2 点を丁寧に書く。

- **速度を含むポテンシャル**なので、教科書でよく見る形 $\dfrac{d}{dt}\dfrac{\partial T}{\partial \dot q}-\dfrac{\partial T}{\partial q}=-\dfrac{\partial U}{\partial q}$ にならず、**追加の項**が要ること (§3、§4)。
- $t \to s$ は、一般の正準変換 (母関数 $F_1$〜$F_4$) とも、単なる独立変数の付け替えとも違う**特殊な変換**であること (§8)。

**前提**: 大学の力学 (ニュートンの運動方程式、ベクトル解析)。解析力学は知らなくてよい (§2 で要点から書く)。SAD は知らなくてよい。

**読む順**: 先に §1 (結論と記号) を見る。ラグランジアンの導出から読むなら §2 → §3 → §4 → §5。解析力学に慣れている人は §3 (追加項) と §8 ($t \to s$) だけでも読める。四極の章 ([04_quad.md](04_quad.md)) は、この章の §9 の結果 (式 (1a.63)〜(1a.67)) を出発点にしている。SAD のヘルプとの対応は §10。

数式の番号 $(1a.n)$ は、この章の中と、他の章からの参照に使う。確認用の sympy スクリプトは `physics_guide/scripts/ham_*.py`、出力は `physics_guide/outputs/ham_*.out.txt` (§11)。

## 1. 結論と記号

- 質量 $m_0$、電荷 $e$ (符号を含む) の粒子が、静的な電磁場の中を運動する。**ラグランジアン**は

$$
L = -\frac{m_0 c^2}{\gamma} - e\,(\phi - \vec A\cdot\vec v)
\tag{1a.1}
$$

  で、ポテンシャル $U = e(\phi - \vec A\cdot\vec v)$ が**速度 $\vec v$ に依存する** (一般化ポテンシャル)。このとき一般化力は $Q_k = -\partial U/\partial q_k + \dfrac{d}{dt}\dfrac{\partial U}{\partial \dot q_k}$ で、第 2 項が**追加の項** (§3)。
- **一般化運動量 (正準運動量)** は $\vec p = \partial L/\partial\vec v = \gamma m_0 \vec v + e\vec A$。機械的な運動量 $\gamma m_0\vec v$ とは $e\vec A$ だけ違う。SAD の $p_x, p_y$ は正準運動量の方 (§3.3)。
- **ハミルトニアン**は $H = c\sqrt{(\vec p - e\vec A)^2 + m_0^2c^2} + e\phi$ (式 (1a.34))。
- **Frenet–Serret 座標**への変換は、母関数 $F_3(\vec p, \vec Q) = -\vec p\cdot\vec r(\vec Q)$ の正準変換で、$\vec P = (p_x,\,p_y,\,(1 + X/\rho)\,p_z)$ (§7)。
- **独立変数を $t$ から $s$ に替える**には、$s$ を座標から独立変数に、$t$ を独立変数から座標に、**入れ替える**。$t$ に共役な運動量は $-H$ で、新しいハミルトニアンは $H' = -P_s$ (§8)。結果は

$$
H' = -eA_s - \left(1 + \frac{X}{\rho}\right)\sqrt{p_{\rm m}^2 - (P_X - eA_X)^2 - (P_Y - eA_Y)^2}
\tag{1a.2}
$$

  ($p_{\rm m}$ は機械的運動量の大きさ)。
- 設計運動量 $P_0$ で割り、$p_{\rm m} = P_0(1 + \delta)$ とおき、横方向の 2 次まで展開すると $H = -\dfrac{eA_s}{P_0} + \dfrac{p_x^2 + p_y^2}{2(1+\delta)} - \dfrac{x}{\rho}(1 + \delta)$。偏向磁石は $x'' = -x/\rho^2 + \delta/\rho$、四極は $x'' = -Kx$、$y'' = +Ky$ ($\delta = 0$ のとき。§9)。

| 量 | 記号 | 約束 |
|---|---|---|
| 質量、電荷 | $m_0$、$e$ | $e$ は符号つき。電子なら $e < 0$ |
| ローレンツ因子、速度 | $\gamma = 1/\sqrt{1 - v^2/c^2}$、$\vec v$ | $\beta = v/c$ |
| 電磁ポテンシャル | $\phi$、$\vec A$ | $\vec E = -\nabla\phi - \partial\vec A/\partial t$、$\vec B = \nabla\times\vec A$ |
| 正準運動量 | $\vec p$ | $\gamma m_0\vec v + e\vec A$ |
| 機械的運動量 | $\vec\pi = \vec p - e\vec A = \gamma m_0 \vec v$ | その大きさ $p_{\rm m} = \lvert\vec\pi\rvert$ |
| Frenet–Serret 座標 | $(X, Y, s)$、$h = 1 + X/\rho$ | $X$ は水平、$Y$ は垂直、$s$ は軌道に沿う距離。$\rho$ は曲率半径 |
| 新しい座標に共役な運動量 | $(P_X, P_Y, P_s)$ | $P_s$ は $s$ 方向の正準運動量 (母関数 $F_3$ で決まる) |
| 設計運動量、偏差 | $P_0$、$\delta = (p_{\rm m} - P_0)/P_0$ | |
| 規格化した運動量 | $p_x = P_X/P_0$、$p_y = P_Y/P_0$、$a_s = eA_s/P_0$ | SAD の $(p_x, p_y)$ と同じ約束 |

以下、章の中で「SAD のヘルプ」は CERN 版の SADHelp を指す。

## 2. ラグランジュ形式

加速器の粒子は自由粒子 (束縛条件がない) なので、仮想変位の原理や束縛条件の取り扱いは要らない。それでも、ラグランジュの運動方程式がどんな前提の上にあるかを知っておくと、§3 の追加の項の意味が見える。そこで、要点を押さえる。

### 2.1 仮想仕事とダランベールの原理 (要点)

質点 $i$ ($i = 1\ldots N$) の運動量を $\vec P_i$ (この節だけの記号)、力を $\vec F_i$ とする。力は、仕事をする**有効力** $\vec F^a_i$ と、仕事をしない**束縛力** $\vec f_i$ の和 ($\vec f_i\cdot\delta\vec x_i = 0$) で、$\delta\vec x_i$ は束縛条件に従う仮想変位。静力学の釣り合いは $\sum_i\vec F_i\cdot\delta\vec x_i = 0$ (仮想仕事の原理)。運動方程式 $\dot{\vec P}_i = \vec F_i$ を「$\vec F_i - \dot{\vec P}_i = 0$ という力の釣り合い」と読むと、$-\dot{\vec P}_i$ を慣性抵抗力として、動力学にも同じ原理が使える (ダランベールの原理):

$$
\sum_{i=1}^{N}\left(\vec F^a_i - \dot{\vec P}_i\right)\cdot\delta\vec x_i = 0
\tag{1a.3}
$$

### 2.2 ラグランジュの運動方程式 (要点)

束縛条件を使って、座標を自由度の数 $K$ の一般化座標 $q_1\ldots q_K$ で書く: $\vec x_i = \vec x_i(q_1,\ldots,q_K,t)$。すると

$$
\dot{\vec x}_i = \sum_k \frac{\partial\vec x_i}{\partial q_k}\dot q_k + \frac{\partial\vec x_i}{\partial t},\qquad
\frac{\partial\dot{\vec x}_i}{\partial\dot q_k} = \frac{\partial\vec x_i}{\partial q_k},\qquad
\delta\vec x_i = \sum_k\frac{\partial\vec x_i}{\partial q_k}\delta q_k
\tag{1a.4}
$$

($\vec x_i$ は $\dot q_k$ に陽に依存しないので、$\dot{\vec x}_i$ は $\dot q_k$ の 1 次式で、その係数が $\partial\vec x_i/\partial q_k$ であることを使った。) $\dot{\vec P}_i = m_i\ddot{\vec x}_i$ の項の変形がこの導出の中心で、

$$
\ddot{\vec x}_i\cdot\frac{\partial\vec x_i}{\partial q_k}
= \frac{d}{dt}\!\left(\dot{\vec x}_i\cdot\frac{\partial\vec x_i}{\partial q_k}\right) - \dot{\vec x}_i\cdot\frac{d}{dt}\frac{\partial\vec x_i}{\partial q_k}
= \frac{d}{dt}\frac{\partial}{\partial\dot q_k}\frac{\dot{\vec x}_i^{\,2}}{2} - \frac{\partial}{\partial q_k}\frac{\dot{\vec x}_i^{\,2}}{2}
\tag{1a.5}
$$

(積の微分の逆。2 つ目の等号は (1a.4) の $\partial\vec x_i/\partial q_k = \partial\dot{\vec x}_i/\partial\dot q_k$ と、$d/dt$ と $\partial/\partial q_k$ が交換すること $\frac{d}{dt}\frac{\partial\vec x_i}{\partial q_k} = \frac{\partial\dot{\vec x}_i}{\partial q_k}$ を使う。) 運動エネルギー $T = \sum_i \tfrac12 m_i\dot{\vec x}_i^{\,2}$ と、**一般化力**

$$
Q_k \equiv \sum_{i=1}^{N}\vec F^a_i\cdot\frac{\partial\vec x_i}{\partial q_k}
\tag{1a.6}
$$

を使うと、(1a.3) は $\sum_k\left[Q_k - \dfrac{d}{dt}\dfrac{\partial T}{\partial\dot q_k} + \dfrac{\partial T}{\partial q_k}\right]\delta q_k = 0$ になる。$\delta q_k$ が独立なので、各 $k$ について

$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k} - \frac{\partial T}{\partial q_k} = Q_k
\tag{1a.7}
$$

力が**保存力**で、ポテンシャルが位置と時間だけの関数 $U(q, t)$ ($\partial U/\partial\dot q_k = 0$) なら、$Q_k = -\partial U/\partial q_k$。そこで $L = T - U$ とおくと、$\partial L/\partial\dot q_k = \partial T/\partial\dot q_k$ なので

$$
\frac{d}{dt}\frac{\partial L}{\partial\dot q_k} - \frac{\partial L}{\partial q_k} = 0
\tag{1a.8}
$$

これが、教科書でよく見るラグランジュの運動方程式。**ここまでの導出には、$U$ が速度に依存しないという前提が入っている**。

### 2.3 変分法でのオイラー・ラグランジュ方程式

同じ方程式を、作用の極値として導く。始点 $x_1$ と終点 $x_2$ で値を固定した関数 $y(x)$ に対し、

$$
J = \int_{x_1}^{x_2} f(x, y, y')\,dx
\tag{1a.9}
$$

の極値を与える $y$ を探す。$y$ に、端点で $0$ になる任意の関数 $\eta(x)$ ($\eta(x_1) = \eta(x_2) = 0$) を足して $Y(x,\alpha) = y(x) + \alpha\eta(x)$ とし、$J(\alpha) = \int f(x, Y, Y')dx$ を $\alpha$ の関数と見る。$\alpha = 0$ が本来の $y$ なので、極値の条件は $\left.\partial J/\partial\alpha\right|_{\alpha=0} = 0$。$\partial Y/\partial\alpha = \eta$、$\partial Y'/\partial\alpha = \eta'$ から

$$
\frac{\partial J}{\partial\alpha} = \int_{x_1}^{x_2}\left(\frac{\partial f}{\partial Y}\eta + \frac{\partial f}{\partial Y'}\eta'\right)dx
= \int_{x_1}^{x_2}\frac{\partial f}{\partial Y}\eta\,dx + \left[\frac{\partial f}{\partial Y'}\eta\right]_{x_1}^{x_2} - \int_{x_1}^{x_2}\frac{d}{dx}\!\left(\frac{\partial f}{\partial Y'}\right)\eta\,dx
\tag{1a.10}
$$

(第 2 項を部分積分した。) 端点の項は $\eta = 0$ で消える。$\alpha = 0$ とすると

$$
\int_{x_1}^{x_2}\left[\frac{\partial f}{\partial y} - \frac{d}{dx}\frac{\partial f}{\partial y'}\right]\eta\,dx = 0
\quad\Rightarrow\quad
\frac{\partial f}{\partial y} - \frac{d}{dx}\frac{\partial f}{\partial y'} = 0
\tag{1a.11}
$$

($\eta$ が任意なので括弧の中が $0$。) $f \to L$、$x \to t$、$y \to q_k$ とすれば (1a.8) に一致する。**この導出も、$L$ の形は問わない**。$L(q, \dot q, t)$ が速度に依存してもそのまま使える。したがって、速度に依存するポテンシャルを扱う道は 2 つある:

1. (1a.7) の一般化力 $Q_k$ を、速度依存の場合に正しく書く (§3.1)。
2. 作用 $\int L\,dt$ の極値の原理を出発点にし、$L$ を適切に選ぶ (§3.2〜§4)。ローレンツ力が (1a.1) の $L$ で表せることを §4 で示す。

## 3. 速度に依存するポテンシャルと、追加の項

### 3.1 一般化力の追加項

ポテンシャル $U(q, \dot q, t)$ が速度に依存するとき、(1a.7) の右辺の一般化力を

$$
Q_k = -\frac{\partial U}{\partial q_k} + \frac{d}{dt}\frac{\partial U}{\partial\dot q_k}
\tag{1a.12}
$$

と書ける力に対しては、$L = T - U$ でラグランジュの運動方程式 (1a.8) がそのまま使える。実際、$L = T - U$ を (1a.8) に入れて整理すると

$$
\frac{d}{dt}\frac{\partial T}{\partial\dot q_k} - \frac{\partial T}{\partial q_k}
= -\frac{\partial U}{\partial q_k} + \frac{d}{dt}\frac{\partial U}{\partial\dot q_k}
\tag{1a.13}
$$

で、右辺が (1a.12) の $Q_k$ そのもの。つまり、

> $U$ が速度を含むと、運動方程式の右辺は $-\partial U/\partial q_k$ **だけでは済まず**、$\dfrac{d}{dt}\dfrac{\partial U}{\partial\dot q_k}$ が加わる。$U$ が速度を含まないときは第 2 項が 0 で、よく見る形に戻る。

逆に、力が (1a.12) の形に書けるかどうかは、個々の力で確かめるしかない。ローレンツ力が書けることは §4 で示す。

注意: 速度に依存する力のすべてが (1a.12) の形になるわけではない。例えば、速度に比例する摩擦力 $f = -\eta\dot x$ は、$D = \eta\dot x^2/2$ を使って $-\partial D/\partial\dot x$ と書き、運動方程式の右辺に $-\partial D/\partial\dot q$ を加える形 (散逸関数) になる。$L = T - U$ の枠を越えて拡張が必要になる別の例で、エネルギーを失う力である。磁場は速度に垂直な力しか与えず、仕事をしないので、(1a.12) の形に収まる。

### 3.2 一般化運動量

ラグランジアン $L(q, \dot q, t)$ から、一般化運動量を

$$
p_k \equiv \frac{\partial L}{\partial\dot q_k}
\tag{1a.14}
$$

と定義する。ラグランジュの方程式 (1a.8) は $\dot p_k = \partial L/\partial q_k$ となる。$U$ が速度を含まなければ $p_k = \partial T/\partial\dot q_k$ で、通常の運動量 (質量 × 速度) になる。**$U$ が速度を含むと、$p_k$ に $-\partial U/\partial\dot q_k$ の分が加わる**。

### 3.3 ビームの場合: 正準運動量と機械的運動量

(1a.1) の $L$ を使う。デカルト座標 $\vec r = (x, y, z)$、速度 $\vec v = \dot{\vec r}$ で、速度は $v^2 = v_x^2 + v_y^2 + v_z^2$。$p_x = \partial L/\partial v_x$ を計算する。$\partial v^2/\partial v_x = 2v_x$ と、$\phi$ が速度を含まないこと、$\vec A\cdot\vec v = A_xv_x + A_yv_y + A_zv_z$ から

$$
p_x = \frac{\partial}{\partial v_x}\!\left[-m_0c^2\sqrt{1 - \frac{v^2}{c^2}}\right] + eA_x
= -m_0c^2\cdot\frac12\left(1 - \frac{v^2}{c^2}\right)^{-1/2}\left(-\frac{2v_x}{c^2}\right) + eA_x
= \gamma m_0 v_x + eA_x
\tag{1a.15}
$$

3 成分をまとめて

$$
\vec p = \gamma m_0\vec v + e\vec A,\qquad
\vec\pi \equiv \vec p - e\vec A = \gamma m_0\vec v
\tag{1a.16}
$$

- $\vec p$ が**正準運動量** (一般化運動量)。ハミルトニアンの $\vec p$ はこれ。SAD の $p_x, p_y$ も、$P_0$ で規格化した正準運動量である。
- $\gamma m_0\vec v$ が**機械的な運動量** (ニュートン力学の「質量 × 速度」で、質量に相対論的質量 $\gamma m_0$ を使ったもの)。ベクトルポテンシャル $\vec A$ のあるところでは、正準運動量と値が違う。軌道の傾き (速度の向き) を決めるのは機械的運動量の方。

確認: sympy で $\partial L/\partial\vec v - (\gamma m_0\vec v + e\vec A) = 0$ ([ham_01_lorentz_potential.out.txt](outputs/ham_01_lorentz_potential.out.txt) (d))。

## 4. ローレンツ力が一般化ポテンシャルで書けること

ローレンツ力 $\vec F = e(\vec E + \vec v\times\vec B)$ が、$U = e(\phi - \vec A\cdot\vec v)$ を使った (1a.12) の形 $\vec F = -\partial U/\partial\vec r + \dfrac{d}{dt}\dfrac{\partial U}{\partial\vec v}$ で書けることを示す。

### 4.1 場をポテンシャルで書く

Maxwell 方程式の $\nabla\cdot\vec B = 0$ から、磁場はベクトルポテンシャルで書ける:

$$
\vec B = \nabla\times\vec A
\tag{1a.17}
$$

(発散のないベクトル場は、あるベクトル場の回転で書ける。) ファラデーの法則 $\nabla\times\vec E + \partial\vec B/\partial t = 0$ に (1a.17) を入れ、$\partial/\partial t$ と $\nabla\times$ が交換することを使うと

$$
\nabla\times\left(\vec E + \frac{\partial\vec A}{\partial t}\right) = 0
\tag{1a.18}
$$

回転が $0$ のベクトル場は、スカラー関数の勾配で書ける。そこで $\vec E + \partial\vec A/\partial t = -\nabla\phi$ とおくと

$$
\vec E = -\nabla\phi - \frac{\partial\vec A}{\partial t}
\tag{1a.19}
$$

ローレンツ力は

$$
\vec F = e\left[-\nabla\phi - \frac{\partial\vec A}{\partial t} + \vec v\times(\nabla\times\vec A)\right]
\tag{1a.20}
$$

### 4.2 ベクトル恒等式

3 つのベクトル場 $\vec a$、$\vec b$ について、次の恒等式が成り立つ:

$$
\nabla(\vec a\cdot\vec b) = (\vec a\cdot\nabla)\vec b + (\vec b\cdot\nabla)\vec a + \vec a\times(\nabla\times\vec b) + \vec b\times(\nabla\times\vec a)
\tag{1a.21}
$$

注意: (i) $(\vec a\cdot\nabla)$ は**スカラーの微分演算子** ($a_x\partial_x + a_y\partial_y + a_z\partial_z$) で、$(\vec a\cdot\nabla)\vec b$ は $\vec b$ の各成分にこの演算子を作用させたベクトル。内積 $\vec a\cdot\vec b$ は順序を替えてよいが、$\vec A\cdot\nabla$ と $\nabla\cdot\vec A$ は別物 (前者は演算子、後者は値)。(ii) $\vec a\times(\vec b\times\vec c)$ と $(\vec a\times\vec b)\times\vec c$ は一般に違う。(iii) どこに微分が掛かっているかを追うこと。

ここで $\vec a = \vec v$ (粒子の速度)、$\vec b = \vec A$ とおく。**$\vec v$ は、粒子の位置 $\vec r$ とは独立な変数**である (位置の関数としての速度場ではない。$\nabla$ は $\vec r$ について微分し、$\vec v$ は一定に保つ)。したがって $\nabla$ が $\vec v$ に掛かる項、すなわち (1a.21) の第 2 項 $(\vec A\cdot\nabla)\vec v$ と第 4 項 $\vec A\times(\nabla\times\vec v)$ は $0$ になる。残りは

$$
\nabla(\vec v\cdot\vec A) = (\vec v\cdot\nabla)\vec A + \vec v\times(\nabla\times\vec A)
\;\;\Rightarrow\;\;
\vec v\times(\nabla\times\vec A) = \nabla(\vec v\cdot\vec A) - (\vec v\cdot\nabla)\vec A
\tag{1a.22}
$$

**成分で確認する (補足)**。$x$ 成分で、右辺の各項を書き下す。左辺 $[\vec v\times(\nabla\times\vec A)]_x = v_y(\partial_xA_y - \partial_yA_x) + v_z(\partial_xA_z - \partial_zA_x)$ ($\nabla\times\vec A$ の $z$ 成分は $\partial_xA_y - \partial_yA_x$、$y$ 成分は $\partial_zA_x - \partial_xA_z$ で、$[\vec v\times\vec C]_x = v_yC_z - v_zC_y$)。右辺は $\partial_x(\vec v\cdot\vec A) = v_x\partial_xA_x + v_y\partial_xA_y + v_z\partial_xA_z$ ($\vec v$ は定数として出る) から $(\vec v\cdot\nabla)A_x = v_x\partial_xA_x + v_y\partial_yA_x + v_z\partial_zA_x$ を引くと $v_y\partial_xA_y + v_z\partial_xA_z - v_y\partial_yA_x - v_z\partial_zA_x$。左辺と一致する ($v_x\partial_xA_x$ が相殺)。$y$、$z$ 成分も同じ。sympy でも、一般の関数 $\vec a(\vec r)$、$\vec b(\vec r)$ に対する (1a.21) の全成分と、$\vec v$ が定数のときの (1a.22) を確かめた ([ham_01_lorentz_potential.out.txt](outputs/ham_01_lorentz_potential.out.txt) (a))。

### 4.3 $\vec A$ の全微分

粒子が動くと、$\vec A$ が感じる位置が $\vec r(t)$ のように変わる。$\vec A = \vec A(\vec r(t), t)$ の時間変化は、鎖則で

$$
\frac{d\vec A}{dt} = \frac{\partial\vec A}{\partial t} + \frac{\partial\vec A}{\partial x}\frac{dx}{dt} + \frac{\partial\vec A}{\partial y}\frac{dy}{dt} + \frac{\partial\vec A}{\partial z}\frac{dz}{dt}
= \frac{\partial\vec A}{\partial t} + (\vec v\cdot\nabla)\vec A
\tag{1a.23}
$$

(成分ごとに $\dfrac{dA_i}{dt} = \partial_tA_i + v_x\partial_xA_i + v_y\partial_yA_i + v_z\partial_zA_i$。$\partial_t\vec A$ は固定した場所での時間変化、$d\vec A/dt$ は粒子が見る時間変化。) 変形して

$$
\frac{\partial\vec A}{\partial t} = \frac{d\vec A}{dt} - (\vec v\cdot\nabla)\vec A
\tag{1a.24}
$$

sympy で、成分ごとに異なる形の関数 $A_i(x, y, z, t)$ に対し、$\vec r(\tau)$ を代入して $\tau$ で微分した結果が $\partial_tA_i + \vec v\cdot\nabla A_i$ に一致することを確かめた ((b))。

### 4.4 ローレンツ力を勾配と全微分の形に

(1a.22) と (1a.24) を (1a.20) に代入する:

$$
\vec F = e\left[-\nabla\phi - \frac{d\vec A}{dt} + (\vec v\cdot\nabla)\vec A + \nabla(\vec v\cdot\vec A) - (\vec v\cdot\nabla)\vec A\right]
= e\,\nabla(-\phi + \vec v\cdot\vec A) - e\frac{d\vec A}{dt}
\tag{1a.25}
$$

($(\vec v\cdot\nabla)\vec A$ が相殺した。)

### 4.5 一般化ポテンシャル

$U = e(\phi - \vec A\cdot\vec v)$ とおいて (1a.12) の右辺を作る。

$$
\frac{\partial U}{\partial\vec v} = -e\vec A
\;\;\Rightarrow\;\;
\frac{d}{dt}\frac{\partial U}{\partial\vec v} = -e\frac{d\vec A}{dt},
\qquad
-\frac{\partial U}{\partial\vec r} = e\nabla(-\phi + \vec A\cdot\vec v)
\tag{1a.26}
$$

($\partial/\partial\vec r$ は $\vec v$ を一定にした位置微分。) 両者の和は (1a.25) に一致する:

$$
\vec F = -\frac{\partial U}{\partial\vec r} + \frac{d}{dt}\frac{\partial U}{\partial\vec v},
\qquad U = e\,(\phi - \vec A\cdot\vec v)
\tag{1a.27}
$$

したがって、ラグランジアンは $L = T - U$ (自由粒子の相対論的な $T$ 相当の部分は $-m_0c^2/\gamma$、次節) で、(1a.1) になる。

**成分で確認する (補足)**。$x$ 成分で書く。ローレンツ力は

$$
F_x = -e\partial_x\phi - e\partial_tA_x + e\left[v_y(\partial_xA_y - \partial_yA_x) + v_z(\partial_xA_z - \partial_zA_x)\right]
$$

一方 (1a.27) の右辺は、$-\partial_xU = -e\partial_x\phi + e(v_x\partial_xA_x + v_y\partial_xA_y + v_z\partial_xA_z)$ と、$\dfrac{d}{dt}\dfrac{\partial U}{\partial v_x} = -e\dfrac{dA_x}{dt} = -e(\partial_tA_x + v_x\partial_xA_x + v_y\partial_yA_x + v_z\partial_zA_x)$ の和。$v_x\partial_xA_x$ が相殺し、$-e\partial_x\phi - e\partial_tA_x + e[v_y\partial_xA_y + v_z\partial_xA_z - v_y\partial_yA_x - v_z\partial_zA_x]$ で、$F_x$ と一致する。**追加の項 $\frac{d}{dt}\frac{\partial U}{\partial v}$ を落とすと、$x$ 成分に $e\,dA_x/dt$ の差が残る** (sympy の (c) の最後の行)。3 成分とも sympy で $Q_i - F_i = 0$ を確認した。

### 4.6 運動方程式

(1a.1) のラグランジュの方程式 (1a.8) を $\vec r$ の成分に使うと、$\partial L/\partial\vec v = \vec p$ (式 (1a.16)) から

$$
\frac{d\vec p}{dt} = \frac{\partial L}{\partial\vec r}
\;\Leftrightarrow\;
\frac{d}{dt}(\gamma m_0\vec v) + e\frac{d\vec A}{dt} = -e\nabla\phi + e\nabla(\vec A\cdot\vec v)
\;\Leftrightarrow\;
\frac{d}{dt}(\gamma m_0\vec v) = e(\vec E + \vec v\times\vec B)
\tag{1a.28}
$$

(右辺の最後の変形は (1a.25) の逆。) つまり、(1a.1) の $L$ から、相対論的なローレンツ力の運動方程式が出る。sympy で、$\left[\dfrac{d}{dt}\dfrac{\partial L}{\partial\vec v} - \dfrac{\partial L}{\partial\vec r}\right] - \left[\dfrac{d}{dt}(\gamma m_0\vec v) - e(\vec E + \vec v\times\vec B)\right] = 0$ を 3 成分で確認した ((d))。

## 5. ルジャンドル変換とハミルトニアン

### 5.1 ルジャンドル変換

2 変数の関数 $f(x, y)$ の全微分を $df = u\,dx + w\,dy$ ($u = \partial f/\partial x$、$w = \partial f/\partial y$) とする。$g \equiv f - ux$ をつくると

$$
dg = df - d(ux) = u\,dx + w\,dy - u\,dx - x\,du = w\,dy - x\,du
\;\Rightarrow\;
\frac{\partial g}{\partial u} = -x,\quad\frac{\partial g}{\partial y} = w
\tag{1a.29}
$$

$g$ は変数が $(x, y)$ から $(u, y)$ に替わった関数になる。これが**ルジャンドル変換**。1 変数なら、$y(x)$ の接線 ($y = ux + z$、$u$ は傾き) の切片 $z = y - ux$ を傾き $u$ の関数として表すことに当たる。

### 5.2 ラグランジアンからハミルトニアンへ (定義から順に)

1 自由度で書く (多自由度は和を取るだけ)。ラグランジアン $L(q, \dot q, t)$ から始める。

**(1) 一般化運動量の定義**

$$
p \equiv \frac{\partial L}{\partial\dot q}
\tag{1a.30a}
$$

**(2) ラグランジュの方程式**は (1a.8) から、(1) を使って $\dfrac{d}{dt}\dfrac{\partial L}{\partial\dot q} = \dfrac{\partial L}{\partial q}$、つまり

$$
\dot p = \frac{\partial L}{\partial q}
\tag{1a.30b}
$$

**(3) $L$ の全微分**。$L$ は $q$、$\dot q$、$t$ の関数なので

$$
dL = \frac{\partial L}{\partial q}dq + \frac{\partial L}{\partial\dot q}d\dot q + \frac{\partial L}{\partial t}dt
= \dot p\,dq + p\,d\dot q + \frac{\partial L}{\partial t}dt
\tag{1a.30c}
$$

(第 1 項は (2)、第 2 項は (1) を使った。)

**(4) ハミルトニアンの定義** (ルジャンドル変換: $\dot q$ を $p$ に取り替える):

$$
H \equiv p\,\dot q - L
\tag{1a.30d}
$$

**(5) $H$ の全微分**。積の微分 $d(p\dot q) = \dot q\,dp + p\,d\dot q$ と (3) から

$$
dH = d(p\dot q) - dL = \dot q\,dp + p\,d\dot q - \dot p\,dq - p\,d\dot q - \frac{\partial L}{\partial t}dt
= \dot q\,dp - \dot p\,dq - \frac{\partial L}{\partial t}dt
\tag{1a.30e}
$$

($p\,d\dot q$ が相殺して、$d\dot q$ が消える。$H$ は $(q, p, t)$ の関数になった。)

**(6) 係数の比較**。$H = H(q, p, t)$ の全微分は $dH = \dfrac{\partial H}{\partial q}dq + \dfrac{\partial H}{\partial p}dp + \dfrac{\partial H}{\partial t}dt$。(5) と $dq$、$dp$、$dt$ の係数を比べて

$$
\frac{\partial H}{\partial p} = \dot q,\qquad
\frac{\partial H}{\partial q} = -\dot p,\qquad
\frac{\partial H}{\partial t} = -\frac{\partial L}{\partial t}
\tag{1a.30f}
$$

最初の 2 つが**ハミルトンの正準方程式** $\dot q = \partial H/\partial p$、$\dot p = -\partial H/\partial q$。$L$ に $t$ が陽に入らないとき ($\partial L/\partial t = 0$)、$\partial H/\partial t = 0$ で、$dH/dt = \frac{\partial H}{\partial q}\dot q + \frac{\partial H}{\partial p}\dot p = -\dot p\dot q + \dot q\dot p = 0$、すなわち $H$ は保存する。多自由度では、(1)〜(6) を各 $k$ で行う ($H = \sum_kp_k\dot q_k - L$)。

### 5.3 荷電粒子のハミルトニアン

$L$ は (1a.1)、$\vec p$ は (1a.16)。$H = \vec p\cdot\vec v - L$:

$$
H = (\gamma m_0\vec v + e\vec A)\cdot\vec v + \frac{m_0c^2}{\gamma} + e(\phi - \vec A\cdot\vec v)
= \gamma m_0v^2 + \frac{m_0c^2}{\gamma} + e\phi
\tag{1a.31}
$$

($e\vec A\cdot\vec v$ が相殺した。) 右辺の前 2 項は

$$
\gamma m_0v^2 + \frac{m_0c^2}{\gamma}
= \gamma m_0\left[v^2 + c^2\left(1 - \beta^2\right)\right] = \gamma m_0c^2
\tag{1a.32}
$$

($1/\gamma = \gamma(1 - \beta^2)$ を使った。) さらに $\gamma m_0c^2$ を機械的運動量 $\vec\pi = \gamma m_0\vec v$ で書く:

$$
(\gamma m_0c^2)^2 = \gamma^2m_0^2c^4 = m_0^2c^4\left[(1 - \beta^2)\gamma^2 + \beta^2\gamma^2\right] = m_0^2c^4 + c^2\gamma^2m_0^2v^2 = m_0^2c^4 + c^2\pi^2
\tag{1a.33}
$$

($\gamma^2(1 - \beta^2) = 1$、$\gamma^2\beta^2c^4m_0^2 = c^2\gamma^2m_0^2v^2$。) 以上から、**ハミルトニアン**

$$
H = \sqrt{m_0^2c^4 + c^2(\vec p - e\vec A)^2} + e\phi
= c\sqrt{(\vec p - e\vec A)^2 + m_0^2c^2} + e\phi
\tag{1a.34}
$$

$\phi = 0$ で静的な場なら、$H$ は粒子の全エネルギー $E$ に等しく、$E^2 = m_0^2c^4 + c^2p_{\rm m}^2$ ($p_{\rm m}$ は機械的運動量の大きさ) を満たす。

確認 (sympy、[ham_02_legendre.out.txt](outputs/ham_02_legendre.out.txt)): $\vec v = (\vec p - e\vec A)/(\gamma m_0)$、$\gamma = \sqrt{1 + (\vec p - e\vec A)^2/(m_0c)^2}$ が (1a.16) を満たし、$\vec p\cdot\vec v - L = c\sqrt{\ldots} + e\phi$ が恒等的に成り立つ。$\partial H/\partial\vec p = \vec v$、および $-\partial H/\partial\vec r|_{\vec p} = \partial L/\partial\vec r|_{\vec v}$ ($\vec v \to \vec v(\vec p)$ を入れた後) も確認した。

## 6. 正準変換

ハミルトニアン系 $(q, p, H)$ から、同じ形のハミルトン方程式を満たす別の系 $(Q, P, H')$ への変換を**正準変換**という。作用 $\delta\int(p\dot q - H)\,dt = 0$ が、新しい変数でも $\delta\int(P\dot Q - H')\,dt = 0$ となるように決める。両者の被積分関数の差が、ある関数 $F$ の全微分 ($dF/dt$) で書ければ、端点を固定した変分でその差は $[\delta F]$ の境界項になって効かない:

$$
(p\dot q - H) - (P\dot Q - H') = \frac{dF}{dt}
\;\Rightarrow\;
dF = p\,dq - P\,dQ - (H - H')\,dt
\tag{1a.35}
$$

($\dot q\,dt = dq$ を使った。) $F$ を**母関数**という。どの 2 つの変数を独立に選ぶかで、4 種類ある。$F_1$ から、全微分の $d(PQ)$、$d(pq)$ の項を足し引きして次の形に移す (ルジャンドル変換)。

| 母関数 | 独立変数 | 定義 | 関係式 | 新しいハミルトニアン |
|---|---|---|---|---|
| $F_1$ | $(q, Q, t)$ | | $p = \dfrac{\partial F_1}{\partial q}$、$P = -\dfrac{\partial F_1}{\partial Q}$ | $H' = H + \dfrac{\partial F_1}{\partial t}$ |
| $F_2$ | $(q, P, t)$ | $F_2 = F_1 + PQ$ | $p = \dfrac{\partial F_2}{\partial q}$、$Q = \dfrac{\partial F_2}{\partial P}$ | $H' = H + \dfrac{\partial F_2}{\partial t}$ |
| $F_3$ | $(p, Q, t)$ | $F_3 = F_1 - pq$ | $q = -\dfrac{\partial F_3}{\partial p}$、$P = -\dfrac{\partial F_3}{\partial Q}$ | $H' = H + \dfrac{\partial F_3}{\partial t}$ |
| $F_4$ | $(p, P, t)$ | $F_4 = F_1 - pq + PQ$ | $q = -\dfrac{\partial F_4}{\partial p}$、$Q = \dfrac{\partial F_4}{\partial P}$ | $H' = H + \dfrac{\partial F_4}{\partial t}$ |

(1a.36)

導き方の例 ($F_2$): $dF_2 = dF_1 + d(PQ) = p\,dq - P\,dQ - (H - H')dt + P\,dQ + Q\,dP = p\,dq + Q\,dP - (H - H')dt$ なので、$F_2(q, P, t)$ の偏微分が $p = \partial F_2/\partial q$、$Q = \partial F_2/\partial P$、$H' - H = \partial F_2/\partial t$。$F_3$、$F_4$ も同様 ($-d(pq)$ で $dq$ の項が $-q\,dp$ に変わる)。使い方 (例 $F_1$): $p = \partial F_1/\partial q$ を $Q$ について解いて $Q(q, p, t)$ を作り、$P = -\partial F_1/\partial Q$ に代入して $P(q, p, t)$ を得る。

(この表の $F_1$〜$F_4$ は一般の母関数で、SAD の端部の係数 F1、F2 ([04_quad.md](04_quad.md) §7) とは別のもの。)

**重要な前提**: この種類の正準変換では、**独立変数 $t$ は変わらない**。変換は、各時刻 $t$ の相空間 $(q, p) \to (Q, P)$ の写像で、時間は媒介変数のまま。§8 の $t \to s$ は、この枠に入らない。

## 7. Frenet–Serret 座標への正準変換

### 7.1 座標系

加速器の基準軌道 $\vec r_0(s)$ ($s$ は軌道に沿う距離) のまわりの座標を取る。単位ベクトルは

$$
\vec t = \frac{d\vec r_0}{ds},\qquad
\vec n = -\rho\frac{d^2\vec r_0}{ds^2},\qquad
\vec b = \vec t\times\vec n
\tag{1a.37}
$$

($\vec t$ は接線、$\vec n$ は主法線 (曲率の中心とは**逆向き**。この章と SAD の規約は、$x$ が外向きに正)、$\vec b$ は陪法線。$\rho$ が曲率半径。) 水平面内の偏向だけ (ねじれ $\tau = 0$) とすると

$$
\frac{d\vec t}{ds} = -\frac{1}{\rho}\vec n,\qquad
\frac{d\vec n}{ds} = \frac{1}{\rho}\vec t,\qquad
\frac{d\vec b}{ds} = 0
\tag{1a.38}
$$

(ねじれがあるときは $d\vec n/ds = \vec t/\rho - \tau\vec b$、$d\vec b/ds = +\tau\vec n$ ($\vec n\cdot\vec b = 0$ を $s$ で微分した $\vec n'\cdot\vec b + \vec n\cdot\vec b' = 0$ から、$\vec b'$ の $\vec n$ 成分は $\vec n'$ の $\vec b$ 成分と逆符号)。この章では使わない。) 粒子の位置は、基準軌道の点 $\vec r_0(s)$ から水平に $X$、垂直に $Y$ 離れた点として

$$
\vec r(X, Y, s) = \vec r_0(s) + X\,\vec n(s) + Y\,\vec b(s)
\tag{1a.39}
$$

(進行方向の変位は考えない。$s$ が進行方向の座標。) 偏向のない直線部では $\rho \to \infty$。

**独立変数に関する注意**: この時点では、$s$ は粒子の座標のひとつ (時刻 $t$ の関数 $s(t)$) で、独立変数は時間 $t$ のまま。独立変数としての $s$ (§8) と、座標としての $s$ を混同しないこと。

### 7.2 母関数 $F_3$

変換前の座標を小文字 $\vec q = (x, y, z)$ (デカルト座標)、運動量を $\vec p = (p_x, p_y, p_z)$ (この節では**デカルト成分**。§1 の記号 $p_x$ とは別。) 変換後の座標を $\vec Q = (X, Y, s)$、運動量を $\vec P = (P_X, P_Y, P_s)$ とする。母関数を

$$
F_3(\vec p, \vec Q) = -\vec p\cdot\vec r(\vec Q)
\tag{1a.40}
$$

($\vec r(\vec Q)$ は (1a.39))。$F_3$ は時間を陽に含まないので $H' = H$。(1a.36) の関係式から

$$
\vec q = -\frac{\partial F_3}{\partial\vec p} = \vec r(\vec Q),
\qquad
\vec P = -\frac{\partial F_3}{\partial\vec Q} = \vec p\cdot\frac{\partial\vec r}{\partial\vec Q}
\tag{1a.41}
$$

第 1 式は、座標が $\vec Q$ の関数として決まる、つまり点変換であることを言う。第 2 式の各成分は

$$
\frac{\partial\vec r}{\partial X} = \vec n,\qquad
\frac{\partial\vec r}{\partial Y} = \vec b,\qquad
\frac{\partial\vec r}{\partial s} = \vec t + X\frac{d\vec n}{ds} = \left(1 + \frac{X}{\rho}\right)\vec t
\tag{1a.42}
$$

(最後は (1a.38) の $d\vec n/ds = \vec t/\rho$ を使った。ねじれがあるなら $Y\,d\vec b/ds$ と $-X\tau\vec b$ の項が加わる。) これを (1a.41) に入れて

$$
\vec P = (P_X, P_Y, P_s) = \left(\vec p\cdot\vec n,\;\vec p\cdot\vec b,\;\left(1 + \frac{X}{\rho}\right)\vec p\cdot\vec t\right)
\tag{1a.43}
$$

$\vec p\cdot\vec n$ などは、運動量の各方向 (水平・垂直・進行方向) の成分で、これを $(p_x, p_y, p_z)$ (Frenet 基底の成分) と書けば $\vec P = (p_x, p_y, (1 + X/\rho)p_z)$。**$s$ 方向の正準運動量だけ、因子 $(1 + X/\rho)$ が付く** (内側と外側で、同じ $ds$ の弧の長さが違うため)。

ベクトルポテンシャルも、$\vec p$ と同じ変換則 (共変ベクトル) に従うので、

$$
\vec A \to (A_X, A_Y, A_s) = \left(A_x,\;A_y,\;\left(1 + \frac{X}{\rho}\right)A_z\right)
\tag{1a.44}
$$

($A_z$ は進行方向の成分。)

### 7.3 変換後のハミルトニアン

$H' = H$ なので、(1a.34) の $(\vec p - e\vec A)^2$ を新しい変数で書き直すだけでよい。$\vec p - e\vec A$ の Frenet 基底の成分が $(P_X - eA_X,\;P_Y - eA_Y,\;(P_s - eA_s)/(1 + X/\rho))$ なので

$$
H = c\sqrt{(P_X - eA_X)^2 + (P_Y - eA_Y)^2 + \frac{(P_s - eA_s)^2}{(1 + X/\rho)^2} + m_0^2c^2} + e\phi
\tag{1a.45}
$$

(座標の形が変わるのは $s$ 成分だけ。)

確認 (sympy、[ham_03_frenet_f3.out.txt](outputs/ham_03_frenet_f3.out.txt)): 円弧 $\vec r_0(s) = \rho(\cos(s/\rho), 0, \sin(s/\rho))$ で、(1a.38)、(1a.42)、(1a.43)、$|\vec p|^2 = P_X^2 + P_Y^2 + P_s^2/(1 + X/\rho)^2$、変換のポアソン括弧 ($\{Q_i, P_j\} = \delta_{ij}$、$\{P_i, P_j\} = \{Q_i, Q_j\} = 0$、つまり正準変換であること)。

## 8. 独立変数 $t \to s$ の入れ替え

### 8.1 なぜ特別な変換か

§6 の正準変換は、独立変数 $t$ を固定し、座標と運動量の組を取り替える。一方、加速器では、粒子が**どの位置 $s$ にあるとき座標と運動量が何か**を知りたい (同じ位置を毎周通る)。したがって独立変数を $s$ にしたい。ところが $s$ は、(1a.45) のハミルトニアンでは**座標のひとつ**で、独立変数だったのは時間 $t$ である。

**単なる付け替えだけで済まない理由**。独立変数を替えるだけなら、$dt = t'\,ds$ ($t' = dt/ds$)、$\dot q = dq/dt = q'/t'$ で、作用は

$$
\delta\int L\,dt = \delta\int\left(\sum p\dot q - H\right)dt = \delta\int\left(\sum p\,\frac{q'}{t'} - H\right)t'\,ds = \delta\int\left(\sum pq' - Ht'\right)ds
\tag{1a.46}
$$

となる。これは形式的には正しい。しかし、ここには 2 つの問題が残る。

1. 和 $\sum$ には座標 $s$ も含まれ、$s' = ds/ds = 1$ という自明な項 $p_s\cdot 1$ が出る。$s$ はもう座標として動かない (独立変数そのもの) のに、運動量 $p_s$ が残る。
2. 時間 $t$ は、独立変数だった立場から、$s$ の関数 $t(s)$ という「従属変数」に変わる。その運動を、座標と同じ正準方程式で追う形になっていない。

**必要なのは入れ替え**: $s$ を座標から独立変数に移し、$t$ を独立変数から**座標に戻す**。「座標」が減る代わりに「座標」が増え、位相空間の次元は同じ (6 次元) になる。これは、同じ $t$ のもとの相空間 $(q, p) \to (Q, P)$ の写像 (§6) ではなく、座標とパラメータの役割を交換する変換 (一般の正準変換の枠に入らない) で、次のように組む。この手続きは Bell (1955) の Appendix 1 に書かれている手続きと同じで、記号の対応は §8.6。

### 8.2 拡張位相空間と新しいラグランジアン

時間 $t$ を座標と見て、$(x, p_x, y, p_y, s, p_s, t, -H)$ で張る空間を**拡張位相空間**と呼ぶ。以下で、$t$ に共役な運動量が $-H$ になることを示す。

作用を、独立変数 $s$ (パラメータ) で書く。$s$ の関数として $x(s)$、$y(s)$、$t(s)$ を座標とみなす。独立変数の替えかたは (1a.46) と同じ、$\int L\,dt = \int L\,t'\,ds$ だが、**$L\,t'$ を新しいラグランジアン**と見る。

$$
\tilde L \equiv L\,t',\qquad
\dot x = \frac{x'}{t'},\;\;\dot y = \frac{y'}{t'},\;\;\dot s = \frac{1}{t'}\;(\because s' = 1)
\tag{1a.47}
$$

$\tilde L$ は、座標 $(x, y, t)$ とその $s$ による微分 $(x', y', t')$ の関数 ($t'$ が増えると速度 $\dot q = q'/t'$ が小さくなる)。座標に共役な運動量は、ラグランジアンを座標の微分で偏微分したもの。$t$ に共役な運動量 $P_t$ は

$$
P_t = \frac{\partial\tilde L}{\partial t'} = L + t'\sum_{i = x, y, s}\frac{\partial L}{\partial\dot q_i}\frac{\partial\dot q_i}{\partial t'}
= L + t'\sum_i p_i\left(-\frac{q_i'}{t'^2}\right) = L - \sum_ip_i\dot q_i = -H
\tag{1a.48}
$$

($\partial\dot q_i/\partial t' = -q_i'/t'^2 = -\dot q_i/t'$ を使い、$p_i = \partial L/\partial\dot q_i$、$H = \sum p_i\dot q_i - L$。和は $i = x, y, s$。)

$x$、$y$ に共役な運動量は

$$
P_{x} = \frac{\partial\tilde L}{\partial x'} = t'\frac{\partial L}{\partial\dot x}\frac{\partial\dot x}{\partial x'} = t'\,p_x\cdot\frac{1}{t'} = p_x,
\qquad P_y = p_y
\tag{1a.49}
$$

つまり $x, y, p_x, p_y$ は恒等 (単なる名前の付け替え)。$t$ だけが、共役運動量 $P_t = -H$ をもつ新しい座標になった。

### 8.3 新しいハミルトニアン

新しい座標 $(x, y, t)$ と運動量 $(P_x, P_y, P_t)$ に対するハミルトニアンは、(1a.30b) と同じ手順で

$$
\tilde H = P_xx' + P_yy' + P_tt' - \tilde L
\tag{1a.50}
$$

(和は $x, y, t$ だけ。$s$ は座標でなくなった。) $\tilde L = Lt' = t'\left(p_x\dot x + p_y\dot y + p_s\dot s - H\right) = p_xx' + p_yy' + p_s - Ht'$ ($t'\dot s = 1$、$t'\dot x = x'$ を使った) と、$P_t = -H$ を入れると

$$
\tilde H = p_xx' + p_yy' - Ht' - \left(p_xx' + p_yy' + p_s - Ht'\right) = -p_s
\tag{1a.51}
$$

**新しいハミルトニアンは $\tilde H = -p_s$。** (サニティチェック: $\tilde L$ を $\int(p_xx' + p_yy' + P_tt' - \tilde H)\,ds = \int(p_xx' + p_yy' - Ht' + p_s)\,ds$ に戻せて、(1a.46) の被積分関数 ($p_s$ の項を含めて) と一致する。つまり同じ作用を、正準形 $\sum PQ' - \tilde H$ に書き直した。)

**$s$ を独立変数とする正準方程式** (§5.2 と同じ手続きを、$t \to s$、$q \to Q = (x, y, t)$、$p \to P = (P_x, P_y, P_t)$、$L \to \tilde L$、$H \to \tilde H$ として繰り返す)。

1. 定義: $P_i = \partial\tilde L/\partial Q_i'$ ((1a.48)(1a.49))。
2. $\tilde L$ は $\int\tilde L\,ds$ の被積分関数なので、(1a.11) より $\dfrac{d}{ds}\dfrac{\partial\tilde L}{\partial Q_i'} = \dfrac{\partial\tilde L}{\partial Q_i}$、すなわち $P_i' = \partial\tilde L/\partial Q_i$。
3. 全微分: $d\tilde L = P_i'\,dQ_i + P_i\,dQ_i' + \dfrac{\partial\tilde L}{\partial s}ds$ (添字 $i$ は和)。
4. $\tilde H = \sum_iP_iQ_i' - \tilde L$ ((1a.50))。
5. 全微分: $d\tilde H = Q_i'\,dP_i + P_i\,dQ_i' - P_i'\,dQ_i - P_i\,dQ_i' - \dfrac{\partial\tilde L}{\partial s}ds = Q_i'\,dP_i - P_i'\,dQ_i - \dfrac{\partial\tilde L}{\partial s}ds$。
6. $dP_i$、$dQ_i$ の係数を比べて

$$
Q_i' = \frac{\partial\tilde H}{\partial P_i},\qquad P_i' = -\frac{\partial\tilde H}{\partial Q_i}
\tag{1a.51a}
$$

つまり $x' = \partial\tilde H/\partial P_x$、$P_x' = -\partial\tilde H/\partial x$、$y$ も同様、$t' = \partial\tilde H/\partial P_t$、$P_t' = -\partial\tilde H/\partial t$。$\tilde H = -p_s$ なので、$x' = -\partial p_s/\partial p_x$ など、**$s$ の運動方程式は $-p_s$ の偏微分で書ける** (ここで $p_s$ は、$P_t = -H$ を解いて $(x, y, P_x, P_y, t, P_t)$ の関数にしたもの。下)。

**これを使える形にする**には、$-p_s$ を新しい変数 $(x, y, P_x, P_y, t, P_t; s)$ だけで書く必要がある。そのために、$P_t = -H$ を**$p_s$ について解く** (逆に言えば、$P_t = -H(\ldots, p_s, \ldots)$ を、新しい運動量 $P_t$ と古い運動量 $p_s$ の関係式として使う)。(1a.45) で $\phi = 0$ (または $\phi$ を含む場合 $H - e\phi = -P_t - e\phi$) として、$E \equiv -P_t$:

$$
\left(\frac{E - e\phi}{c}\right)^2 = (P_X - eA_X)^2 + (P_Y - eA_Y)^2 + \frac{(P_s - eA_s)^2}{(1 + X/\rho)^2} + m_0^2c^2
$$

$$
\Rightarrow\;
P_s = eA_s + \left(1 + \frac{X}{\rho}\right)\sqrt{\left(\frac{E - e\phi}{c}\right)^2 - m_0^2c^2 - (P_X - eA_X)^2 - (P_Y - eA_Y)^2}
\tag{1a.52}
$$

(粒子が $+s$ 向きに進む枝の平方根を選ぶ。) 根号の中の $\left(\frac{E - e\phi}{c}\right)^2 - m_0^2c^2$ は、機械的運動量の大きさの 2 乗 $p_{\rm m}^2$ である ($E - e\phi = \gamma m_0c^2$ と (1a.33) から)。したがって

$$
\tilde H = -P_s = -eA_s - \left(1 + \frac{X}{\rho}\right)\sqrt{p_{\rm m}^2 - (P_X - eA_X)^2 - (P_Y - eA_Y)^2}
\tag{1a.53}
$$

これが (1a.2)。**静的な場** ($\phi = 0$、$\vec A$ が時間によらない) では、$\tilde H$ が $t$ を含まないので $P_t = -E$ は保存し、$p_{\rm m} = \sqrt{(E/c)^2 - m_0^2c^2}$ は一定の定数 ($P_0(1 + \delta)$) として扱える。

### 8.4 元の運動方程式と同じか (確認)

$\tilde H$ のハミルトン方程式 $dx/ds = \partial\tilde H/\partial P_X$ などが、元の方程式 ($dx/dt$ など) を $ds/dt = \dot s$ で割ったものになることを、陰関数の微分で確かめる。$H(X, P_X, P_s) = -P_t$ から $P_s(X, P_X, P_t)$ を定めると、$\partial P_s/\partial P_X = -(\partial H/\partial P_X)/(\partial H/\partial P_s)$、同様に $X$ と $P_t$ について。したがって

$$
\frac{dX}{ds} = \frac{\partial\tilde H}{\partial P_X} = -\frac{\partial P_s}{\partial P_X} = \frac{\partial H/\partial P_X}{\partial H/\partial P_s} = \frac{\dot X}{\dot s},
\qquad
\frac{dP_X}{ds} = -\frac{\partial\tilde H}{\partial X} = -\frac{\partial H/\partial X}{\partial H/\partial P_s} = \frac{\dot P_X}{\dot s},
\qquad
\frac{dt}{ds} = \frac{\partial\tilde H}{\partial P_t} = \frac{1}{\partial H/\partial P_s} = \frac{1}{\dot s}
\tag{1a.54}
$$

($\dot X = \partial H/\partial P_X$、$\dot P_X = -\partial H/\partial X$、$\dot s = \partial H/\partial P_s$。$H = -P_t$ が一定の面の上での微分。) sympy で、(1a.45) の $H$ (磁場のみ、$\vec A$ が $X, Y$ の関数) について、(1a.52) の解を使った $\tilde H$ から、(1a.54) の 5 本の式が一致することを確かめた ([ham_04_t_to_s.out.txt](outputs/ham_04_t_to_s.out.txt) (c))。同じスクリプトの (a)(b) では、Frenet 座標のビームのラグランジアン $L = -m_0c^2\sqrt{1 - (\dot x^2 + \dot y^2 + h^2\dot s^2)/c^2} + e(A_x\dot x + A_y\dot y + A_s\dot s)$ から、$P_t + H = 0$、$P_x = p_x$、$\tilde H + p_s = 0$ を確かめた。

### 8.5 $t$ と $z$

$t$ の運動 ($dt/ds = \partial\tilde H/\partial P_t$) を、粒子の「バンチ内の位置」$z$ に変える変換は別にある。SAD は、母関数 $G(p_t, z; s) = z\,p(p_t) + p_t\,t_0(s)$ ([01_variables_and_conventions.md](01_variables_and_conventions.md) §3 の式 (1.6)) で $(t, p_t) \to (z, p)$ とし、ハミルトニアンに $E/(p_0v_0)$ の項が加わる (式 (1.1))。この章は、そこに行く前の、横方向の運動だけを扱う ($\delta$ は一定の定数)。

### 8.6 Bell (1955) Appendix 1 との対応

M. R. Bell, "Non-linear equations of motion in the synchrotron", A.E.R.E. T/M 125 (1955) の Appendix 1 (原本 9〜10 ページ。原本の画像で確認した) は、$t$ 以外の座標を独立変数に取り替えたときの新しいハミルトニアンを、次の流れで導いている。この節は、その流れと §8.2〜§8.3 の対応を示す。

1. **作用の書き換え**。ハミルトンの原理 $\delta\int L\,dt = 0$ (Bell の (14)) で、独立変数を $t$ から別の変数 $u$ に替えると $\delta\int L\,\dfrac{dt}{du}\,du = 0$、すなわち $\delta\int M\,du = 0$ ((14)')。$M = L\,dt/du$ は $(q_n, t, \partial q_n/\partial u, \partial t/\partial u)$ の関数。$L = \sum p_n\dot q_n - H$ を使うと、(14)' は $\delta\int\left(\sum p_n\dfrac{dq_n}{du} - H\dfrac{dt}{du}\right)du = 0$。(本章の (1a.46) に当たる。)
2. **対称性の議論**。この式は、$-H$ を運動量、$t$ を座標とみなすと、$\sum p\,dq$ の形と対称になる。この対称性から、「$t$ が独立変数のときのハミルトニアンが $H$ なら、$q_n$ が独立変数のときのハミルトニアンは $-p_n$」と Bell は言う。
3. **計算による確認**。独立変数を $q_r$ にとり ($u = q_r$)、$M = L\,t'$ ($' = d/dq_r$)。$t$ に共役な運動量 $P = \partial M/\partial t' = L - \sum_n(\partial L/\partial\dot q_n)\,q_n'/t' = L - \sum p_n\dot q_n = -H$。ほかの座標 ($q_r$ 以外) の運動量は $P_j = \partial M/\partial q_j' = p_j$ (古い運動量のまま)。新しいハミルトニアンは $R = \sum_jP_jq_j' - M = -t'p_r\dot q_r = -p_r$ ($p_r$ は $q_r$ に共役な、古い運動量。新しい独立変数 $q_r$ に対応する運動量で、**符号が逆になる**)。

| Bell (1955) | この章 | 備考 |
|---|---|---|
| $L$、$H$ (原文は花文字)、$p_n$、$q_n$ | $L$、$H$、$p_i$、$q_i$ | 時間 $t$ を独立変数とする古い系 |
| 新しい独立変数 $u$、特に $q_r$ | $s$ ($q_r = s$) | |
| $M = L\,dt/dq_r = L\,t'$ | $\tilde L = L\,t'$ (式 (1a.47)) | |
| $P = \partial M/\partial t' = -H$ | $P_t = -H$ (式 (1a.48)) | $t$ が座標になり、共役運動量が $-H$ |
| $P_j = \partial M/\partial q_j' = p_j$ | $P_x = p_x$、$P_y = p_y$ (式 (1a.49)) | $j$ は $q_r$ 以外の座標 |
| $R = \sum P_jq_j' - M = -p_r$ | $\tilde H = -p_s$ (式 (1a.50)(1a.51)) | 新しいハミルトニアン。符号が逆 |

**Bell の Appendix 1 と本章の違い**: (i) 本章は、$p_s$ を $P_t = -H$ から解いて $(x, y, P_x, P_y, t, P_t)$ の関数にする段階 ((1a.52)(1a.53)) と、元の方程式との一致の確認 ((1a.54)) を加えている。Appendix 1 の 9〜10 ページには、この部分は見当たらなかった。(ii) 本章の Frenet–Serret 座標、規格化、展開 (§7、§9) に当たる記述は、Appendix 1 にはない。**Bell の本文 (Appendix 1 以外のページ) は、本章では内容を確認していない (未確認)**。

## 9. 規格化と展開、偏向磁石と四極の運動方程式

### 9.1 $P_0$ と $\delta$ による規格化

加速器では、設計運動量 $P_0$ で割った量を使う。$H'$ ($= \tilde H$) を $P_0$ で割り、運動量を $p_x = P_X/P_0$、$p_y = P_Y/P_0$、$a_s = eA_s/P_0$ とする。$A_X = A_Y = 0$ (進行方向の磁場はない) として

$$
H = \frac{\tilde H}{P_0} = -a_s - \left(1 + \frac{X}{\rho}\right)\sqrt{\left(\frac{p_{\rm m}}{P_0}\right)^2 - p_x^2 - p_y^2}
\tag{1a.55}
$$

この規格化は、運動量を定数 $1/P_0$ 倍し、ハミルトニアンも同じ倍率にしたもの (ハミルトンの方程式の形は変わらない: $dx/ds = \partial(\tilde H/P_0)/\partial(P_X/P_0)$ 等。スケール変換)。運動量の偏差を

$$
\delta \equiv \frac{p_{\rm m} - P_0}{P_0},\qquad \frac{p_{\rm m}}{P_0} = 1 + \delta
\tag{1a.56}
$$

とおくと

$$
H = -a_s - \left(1 + \frac{x}{\rho}\right)\sqrt{(1 + \delta)^2 - p_x^2 - p_y^2}
\tag{1a.57}
$$

($X \to x$、$Y \to y$ と小文字に戻した。) **$p_x$ は規格化した正準運動量**で、角度 (傾き) とは $x' = p_x/(1 + \delta)$ の関係 (後の式 (1a.60))。傾きの意味の「$x'$」と $p_x$ を区別して書く。

### 9.2 2 次の展開 (1 段ずつ)

**成り立つ条件**: (i) $x/\rho \ll 1$ (偏向磁石の $\rho$ が数 m〜数十 m、$x$ が mm 程度なら $10^{-3}$〜$10^{-4}$)。(ii) $p_x, p_y \ll 1 + \delta$ (角度が mrad 程度なら $10^{-3}$)。$\delta$ は小さくなくてもよい ($\delta > -1$)。**次数の数え方**: $x$、$y$、$p_x$、$p_y$ を 1 次の微小量とし、$\delta$ は数えず (厳密に残す)、2 次までを残す。

**手順 1: 根号の中身を因数分解する**。

$$
u \equiv \frac{p_x^2 + p_y^2}{(1 + \delta)^2},\qquad
\sqrt{(1 + \delta)^2 - p_x^2 - p_y^2} = (1 + \delta)\sqrt{1 - u}
\tag{1a.58a}
$$

$u$ は $(p_x, p_y)$ の 2 次の量 ($\ll 1$)。

**手順 2: 根号を $u$ で展開する** ($|u| < 1$ で成り立つテイラー展開):

$$
\sqrt{1 - u} = 1 - \frac{u}{2} - \frac{u^2}{8} - \frac{u^3}{16} - \cdots
\tag{1a.58b}
$$

$u$ は 2 次、$u^2$ は 4 次、$u^3$ は 6 次。**$u$ の 1 次まで (つまり 2 次まで) 残し**、$u^2$ 以降 (4 次以上) を捨てる。捨てた最初の項は $-u^2/8$ で、残した $u/2$ に対する比は $u/4$ (角度が 1 mrad なら $\sim 2.5\times10^{-7}$)。

**手順 3: $(1+\delta)$ をかける**:

$$
\sqrt{(1 + \delta)^2 - p_x^2 - p_y^2} \approx (1 + \delta)\left(1 - \frac{u}{2}\right) = (1 + \delta) - \frac{p_x^2 + p_y^2}{2(1 + \delta)}
\tag{1a.58c}
$$

(第 2 項: $(1+\delta)\cdot\dfrac{u}{2} = (1 + \delta)\dfrac{p_x^2 + p_y^2}{2(1 + \delta)^2} = \dfrac{p_x^2 + p_y^2}{2(1 + \delta)}$。) **$1/(1+\delta)$ は展開せず、そのまま残す**。小さい $\delta$ で展開すれば $1/(1 + \delta) = 1 - \delta + \delta^2 - \cdots$ だが、これを 1 次で切ると $\delta$ の 2 次以上を捨てることになる。ここでは行わない。

**手順 4: $(1 + x/\rho)$ をかける**。$(1 + x/\rho)$ は 1 次までで終わっている (展開不要)。掛け算の各項を、次数ごとに分ける:

$$
\left(1 + \frac{x}{\rho}\right)\left[(1 + \delta) - \frac{p_x^2 + p_y^2}{2(1 + \delta)}\right]
= \underbrace{(1 + \delta)}_{0\text{ 次}} + \underbrace{\frac{x}{\rho}(1 + \delta)}_{1\text{ 次}} - \underbrace{\frac{p_x^2 + p_y^2}{2(1 + \delta)}}_{2\text{ 次}} - \underbrace{\frac{x}{\rho}\,\frac{p_x^2 + p_y^2}{2(1 + \delta)}}_{3\text{ 次}}
\tag{1a.58d}
$$

**手順 5: $H$ に代入して符号を直す** ($H = -a_s - [\text{上}]$):

$$
H = -a_s - (1 + \delta) - \frac{x}{\rho}(1 + \delta) + \frac{p_x^2 + p_y^2}{2(1 + \delta)} + \frac{x}{\rho}\,\frac{p_x^2 + p_y^2}{2(1 + \delta)} + O(p^4)
\tag{1a.58e}
$$

**手順 6: 捨てる項**。(a) 定数項 $-(1 + \delta)$: $x$、$p$ を含まないので、$\partial H/\partial x$ も $\partial H/\partial p$ も変えず、運動に影響しない。(b) 3 次の項 $x(p_x^2 + p_y^2)/(2\rho(1 + \delta))$: 2 次の項 $p^2/(2(1 + \delta))$ に対する比が $x/\rho$ ($\ll 1$)。(c) $O(p^4)$ (手順 2 の $u^2$ の項)。残りは

$$
H = -a_s + \frac{p_x^2 + p_y^2}{2(1 + \delta)} - \frac{x}{\rho}(1 + \delta)
\tag{1a.59}
$$

($a_s$ の中の $x$、$y$ の展開は、磁場の形に応じて各要素で行う (§9.4、§9.5)。) 以上の各展開は sympy の `series` で確認した ([ham_06_series.out.txt](outputs/ham_06_series.out.txt) (a)(b))。$\varepsilon$ 展開で、$\varepsilon^0$ が $-(1 + \delta)$、$\varepsilon^1$ が $-(x/\rho)(1 + \delta)$、$\varepsilon^2$ が $(p_x^2 + p_y^2)/(2(1 + \delta))$、$\varepsilon^3$ が $x(p_x^2 + p_y^2)/(2\rho(1 + \delta))$。別のスクリプト ([ham_05_expansion.out.txt](outputs/ham_05_expansion.out.txt) (a)) でも同じ。

### 9.3 ベクトルポテンシャルと磁場の関係

$a_s$ に入れる $A_s$ は、**Frenet 座標のベクトルポテンシャルの $s$ 成分** ((1a.44)、$A_s = (1 + x/\rho)A_z$) である。磁場との関係は、直交曲線座標 (スケール因子が $1$、$1$、$1 + X/\rho$) での回転の公式から

$$
B_X = \frac{1}{1 + X/\rho}\frac{\partial A_s}{\partial Y},\qquad
B_Y = -\frac{1}{1 + X/\rho}\frac{\partial A_s}{\partial X}
\tag{1a.60a}
$$

($A_X = A_Y = 0$ のとき。$B_s = 0$。) デカルト座標のスカラーポテンシャルによる磁場の解は、この関係でベクトルポテンシャルに直してから使う。sympy で、$A_t$ を多項式にして、デカルト座標で $\nabla\times\vec A$ を直接計算した値が (1a.60a) と一致することを、いくつかの点で確かめた ([ham_03_frenet_f3.out.txt](outputs/ham_03_frenet_f3.out.txt) (e))。

### 9.4 偏向磁石

一様な垂直磁場 $B_X = 0$、$B_Y = B_0$ (一定) とする。(1a.60a) から $\partial A_s/\partial X = -(1 + X/\rho)B_0$ なので、$X$ について積分し、ゲージを $Y$、$s$ の関数がない形に選ぶと

$$
A_s = -\left(x + \frac{x^2}{2\rho}\right)B_0
\tag{1a.60b}
$$

(確認: $B_X = 0$、$B_Y = B_0$、[ham_05_expansion.out.txt](outputs/ham_05_expansion.out.txt) (b)。) 設計軌道が基準軌道 (曲率半径 $\rho$) であること、つまり設計軌道に沿って運動量 $P_0$ の粒子が曲がる条件は

$$
\frac{eB_0}{P_0} = \frac{1}{\rho}
\tag{1a.61}
$$

($e$ の符号を含めて $\rho$ の符号が決まる。)

**ステップ 1: $a_s$ を求める**。$a_s = eA_s/P_0 = -\left(x + \dfrac{x^2}{2\rho}\right)\dfrac{eB_0}{P_0} = -\dfrac{1}{\rho}\left(x + \dfrac{x^2}{2\rho}\right)$ ((1a.61) を使った)。この $a_s$ は $x$ の 2 次までで終わっている (1 次と 2 次の項だけ)。

**ステップ 2: $H$ に代入** ((1a.59) の $-a_s$ を使う):

$$
H = \frac{1}{\rho}\left(x + \frac{x^2}{2\rho}\right) + \frac{p_x^2 + p_y^2}{2(1 + \delta)} - \frac{x}{\rho}(1 + \delta)
$$

**ステップ 3: 項を整理する**。$x$ の 1 次の項: $\dfrac{x}{\rho} - \dfrac{x}{\rho}(1 + \delta) = -\dfrac{\delta\,x}{\rho}$ ($x/\rho$ が相殺し、$\delta$ に比例する項だけ残る)。

$$
H = \frac{x^2}{2\rho^2} - \frac{\delta\,x}{\rho} + \frac{p_x^2}{2(1 + \delta)} + \frac{p_y^2}{2(1 + \delta)}
\tag{1a.62}
$$

**ステップ 4: 偏微分を 1 項ずつ** (4 つの項を $T_1 = \dfrac{x^2}{2\rho^2}$、$T_2 = -\dfrac{\delta x}{\rho}$、$T_3 = \dfrac{p_x^2}{2(1+\delta)}$、$T_4 = \dfrac{p_y^2}{2(1+\delta)}$ とする。$\delta$ は定数):

| 項 | $\partial/\partial p_x$ | $\partial/\partial x$ | $\partial/\partial p_y$ | $\partial/\partial y$ |
|---|---|---|---|---|
| $T_1 = x^2/(2\rho^2)$ | $0$ | $x/\rho^2$ | $0$ | $0$ |
| $T_2 = -\delta x/\rho$ | $0$ | $-\delta/\rho$ | $0$ | $0$ |
| $T_3 = p_x^2/(2(1+\delta))$ | $p_x/(1+\delta)$ | $0$ | $0$ | $0$ |
| $T_4 = p_y^2/(2(1+\delta))$ | $0$ | $0$ | $p_y/(1+\delta)$ | $0$ |
| 合計 | $p_x/(1+\delta)$ | $x/\rho^2 - \delta/\rho$ | $p_y/(1+\delta)$ | $0$ |

(例: $\partial T_3/\partial p_x = \dfrac{2p_x}{2(1+\delta)} = \dfrac{p_x}{1+\delta}$。)

**ステップ 5: ハミルトンの方程式** ((1a.51a)。$s$ が独立変数): $x' = \partial H/\partial p_x$、$p_x' = -\partial H/\partial x$ など。表の合計から

$$
x' = \frac{p_x}{1 + \delta},\qquad
p_x' = -\frac{x}{\rho^2} + \frac{\delta}{\rho},\qquad
y' = \frac{p_y}{1 + \delta},\qquad
p_y' = 0
\tag{1a.63}
$$

**ステップ 6: 2 階の方程式にする**。$\delta$ は定数なので $x'' = p_x'/(1 + \delta)$。$1/(1+\delta) = 1 - \delta + \delta^2 - \cdots$ を $p_x'$ にかけて展開すると $x'' = (-x/\rho^2 + \delta/\rho)(1 - \delta + \cdots) = -x/\rho^2 + \delta/\rho + \delta x/\rho^2 - \delta^2/\rho + \cdots$。$x$、$\delta$ について 1 次までを残し、$\delta x/\rho^2$、$\delta^2/\rho$ (2 次) を捨てると

$$
x'' = -\frac{x}{\rho^2} + \frac{\delta}{\rho},\qquad y'' = 0
\tag{1a.64}
$$

(sympy の `series` で $\delta$ の 2 次までを確認: [ham_06_series.out.txt](outputs/ham_06_series.out.txt) (c)(e)。) $p_x$ を $x'$ とみなして ($\delta \to 0$ の分母で) も同じ式になる。

この符号は、偏向磁石の弱い集束 ($-x/\rho^2$) と、運動量の大きい粒子 ($\delta > 0$) が偏向されにくく外側 ($x > 0$ は曲率の中心と逆向き) へずれる効果 ($+\delta/\rho$) と読める。分散関数は、$x'' + x/\rho^2 = \delta/\rho$ の周期解として $\eta(s) = \rho\,[1 - \cos(s/\rho)]$ 型になる。


### 9.5 四極磁石

直線部 ($\rho \to \infty$) で、$B_X = B'y$、$B_Y = B'x$ ($B'$ は勾配)。(1a.60a) から $\partial A_s/\partial y = B'y$、$-\partial A_s/\partial x = B'x$ なので

$$
A_s = -\frac{B'}{2}\left(x^2 - y^2\right),\qquad
a_s = -\frac{K}{2}\left(x^2 - y^2\right),\quad K \equiv \frac{eB'}{P_0}
\tag{1a.65}
$$

(確認: $\partial_yA_s = B'y$、$-\partial_xA_s = B'x$。これが [04_quad.md](04_quad.md) の式 (4.1)(4.2)。) $1/\rho = 0$ なので (1a.59) の $-(x/\rho)(1 + \delta)$ は消え、$-a_s = +K(x^2 - y^2)/2$ を入れる:

$$
H = \frac{K}{2}\left(x^2 - y^2\right) + \frac{p_x^2 + p_y^2}{2(1 + \delta)}
\tag{1a.66}
$$

**偏微分を 1 項ずつ**:

| 項 | $\partial/\partial p_x$ | $\partial/\partial x$ | $\partial/\partial p_y$ | $\partial/\partial y$ |
|---|---|---|---|---|
| $Kx^2/2$ | $0$ | $Kx$ | $0$ | $0$ |
| $-Ky^2/2$ | $0$ | $0$ | $0$ | $-Ky$ |
| $p_x^2/(2(1+\delta))$ | $p_x/(1+\delta)$ | $0$ | $0$ | $0$ |
| $p_y^2/(2(1+\delta))$ | $0$ | $0$ | $p_y/(1+\delta)$ | $0$ |

ハミルトンの方程式 ($x' = \partial H/\partial p_x$、$p_x' = -\partial H/\partial x$) から

$$
x' = \frac{p_x}{1 + \delta},\;\; p_x' = -Kx,\qquad
y' = \frac{p_y}{1 + \delta},\;\; p_y' = +Ky
\;\;\Rightarrow\;\;
x'' = -\frac{Kx}{1 + \delta},\;\; y'' = +\frac{Ky}{1 + \delta}
\tag{1a.67}
$$

($x'' = p_x'/(1 + \delta)$)。$\delta = 0$ では $x'' = -Kx$、$y'' = +Ky$。$\delta \ne 0$ の形は、四極の章の式 (4.8) の $w^2 = K/p$ ($p = 1 + \delta$)。$K > 0$ が水平収束。sympy で偏微分を確認 ([ham_06_series.out.txt](outputs/ham_06_series.out.txt) (d))。

### 9.6 展開の次数の目安 (数値)

(1a.59) は、$H$ の 3 次の項 $x(p_x^2 + p_y^2)/(2\rho(1 + \delta))$ などを落としている。偏向磁石 ($\rho = 10$ m、$\delta = 0$) で、厳密な (1a.57) の運動方程式を DOP853 で 3 m 積分した結果と、(1a.62) の線形の運動方程式の差は、初期振幅 $a$ ($x = a$、$p_x = a/2$) が 1・2・4・8 mm のとき $1.9\times10^{-7}$、$7.6\times10^{-7}$、$3.1\times10^{-6}$、$1.2\times10^{-5}$ で、振幅が 2 倍になるごとに 4 倍になる ([ham_05_expansion.out.txt](outputs/ham_05_expansion.out.txt) (e))。落とした項は $H$ の 3 次で、運動方程式への寄与 (力) は振幅の 2 次、という予想どおり。

## 10. SAD のヘルプ、pysad・K64 との対応

### 10.1 ヘルプのハミルトニアンとの対応

SAD のヘルプ (dynamics > Hamiltonian) の式は、[01_variables_and_conventions.md](01_variables_and_conventions.md) の式 (1.1):

$$
H_{\rm SAD} = -\left(\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2} + a_s\right)\left(1 + \frac{x}{\rho}\right) + \frac{E}{p_0v_0}
$$

この章の (1a.55) $H = -a_s^{\rm (本章)} - (1 + x/\rho)\sqrt{p^2 - p_x^2 - p_y^2}$ ($p = p_{\rm m}/P_0 = 1 + \delta$、$a_x = a_y = 0$) との対応は次のとおり。

| 項目 | この章 | SAD のヘルプ |
|---|---|---|
| 運動量 | $P_X, P_Y$ を $P_0$ で割った $p_x, p_y$ | $p_x, p_y$ (同じ) |
| 全運動量の規格化値 | $p = p_{\rm m}/P_0 = 1 + \delta$ | $p = 1 + \delta$ (同じ) |
| ベクトルポテンシャルの $s$ 成分 | **共変成分** $A_s = (1 + x/\rho)A_t$。式 (1a.55) の $-a_s^{\rm (本章)}$ | **接線方向の物理成分** $A_t$ ($a_s = eA_t/p_0$)。式 (1.1) は $-a_s(1 + x/\rho)$ |
| 同じ量の書き方 | $-a_s^{\rm (本章)} = -(1 + x/\rho)\,a_s^{\rm (SAD)}$ | |
| 縦方向 | $(t, P_t = -E)$ をそのまま座標・運動量にする (§8)。静的な場では $E$ は定数 | $(z, p)$ に正準変換 ([01_variables_and_conventions.md](01_variables_and_conventions.md) §3)。そのとき $+E/(p_0v_0)$ の項が加わる |

したがって、**この章の式 (1a.55) は、SAD のヘルプの式から、$E/(p_0v_0)$ の項 (経路の遅れ $z$ の運動を作る項。$\delta$ が一定のとき横方向の運動に影響しない) を除いたもの**になる (記号 $a_s$ の意味の違いを揃えれば)。横方向の運動方程式は、どちらから出しても同じ。01 章は、ヘルプのラグランジアンから $(t, p_t) \to (z, p)$ の母関数まで書く。この章は、それより基本的な解析力学の部分 (速度依存のポテンシャル、ローレンツ力の一般化ポテンシャル、$t \to s$ の入れ替え) と、偏向・四極の運動方程式への展開を詳しく書く。

### 10.2 pysad・K64 の対応

- 四極の本体のハミルトニアン (04 章の式 (4.3)) は、この章の (1a.57) と (1a.65) で $1/\rho = 0$ としたもの。$H_2$ と $\Delta H$ への分割、厳密解、追跡の写像と転送行列は [04_quad.md](04_quad.md)。
- SAD は、$H$ を上のように 2 次で打ち切らず、根号を厳密なまま残し、2 次の部分 $H_2$ (厳密に解ける) と残り $\Delta H$ (運動学的なドリフトとして解ける) に分けて、交互に適用する。(1a.59) は $H_2$ に相当する (ただし $\delta$ に関する扱いと、$x/\rho$ の項を含む偏向磁石は [05_bend.md] (予定))。
- ヘルプの転送は $\exp(:H:)$ の積で書かれている。この記法 (ポアソン括弧とリー演算子) の意味と、ヘルプの $H$ との符号の対応 (推定) は [C_poisson_lie.md](C_poisson_lie.md)。
- この章で使う SAD の記号: 全運動量の規格化値 $p = 1 + \delta$ は Fortran の `pr` (追跡の `g(i)` が $\delta$、CALC の `cod(6)` が $\delta$)、四極の強さは `akk` ($= K1/L$) ([01_variables_and_conventions.md](01_variables_and_conventions.md) §2 と §8 の表)。この章の導出自体は、Fortran の変数 (`xs1`、`dc1`、`w1` など) を使わない。それらは [04_quad.md](04_quad.md) の §16 の表。
- 偏向磁石 (BEND) の K64 の実装 (Fortran と pysad) の対応は、この章には書かない (偏向磁石の章で)。この章は物理の導出に限る。

## 11. 確認スクリプトと出力

sympy (1.14) で式を確かめる。`py -3.12 physics_guide/scripts/ham_01_lorentz_potential.py` のように、リポジトリのルートから動かす (数秒)。結果は `physics_guide/outputs/<スクリプト名>.out.txt`。

| スクリプト | 節 | 確認すること |
|---|---|---|
| `ham_01_lorentz_potential.py` | §3.3、§4 | ベクトル恒等式 (4 項、$\vec v$ が独立のとき 2 項)、$dA/dt$ の鎖則、一般化力 $Q = -\partial U/\partial r + d(\partial U/\partial v)/dt$ が $e(\vec E + \vec v\times\vec B)$ と一致、$\vec p = \gamma m_0\vec v + e\vec A$、ラグランジュ方程式がローレンツ力の運動方程式 |
| `ham_02_legendre.py` | §5 | $\vec v(\vec p)$、$\vec p\cdot\vec v - L = c\sqrt{\ldots} + e\phi$、$\partial H/\partial\vec p = \vec v$、$-\partial H/\partial\vec r = \partial L/\partial\vec r$ |
| `ham_03_frenet_f3.py` | §7、§9.3 | $\partial\vec r/\partial Q$、$\vec P = \vec p\cdot\partial\vec r/\partial\vec Q$、ポアソン括弧、$\lvert\vec p\rvert^2$ の形、$A_s$ と磁場 (回転を直接計算) |
| `ham_04_t_to_s.py` | §8 | $P_t = -H$、$P_{x,y} = p_{x,y}$、$\tilde H = -p_s$、$P_s$ の解、ハミルトンの方程式が元の方程式 $/\dot s$ に一致 |
| `ham_05_expansion.py` | §9 | $\varepsilon$ 展開 (各次の係数)、$A_s \to B$、偏向磁石と四極の運動方程式、厳密な $\tilde H$ と 2 次の展開の差の振幅依存性 |
| `ham_06_series.py` | §9.2、§9.4、§9.5 | $\sqrt{1 \pm u}$ と $1/(1+\delta)$ の展開、$H$ の因子ごとの展開と次数、偏向磁石・四極の偏微分の 1 項ずつの表、$x'' = p_x'/(1+\delta)$ の展開、$\Delta H$ の展開 (四極の章の (4.6)) |
| `ham_07_lie_operator.py` | 付録 C | ポアソン括弧、$\exp(:f:)$ のドリフト・集束四極・6 極の例、合成の順序、BCH、2 次の $f$ と行列 ([C_poisson_lie.md](C_poisson_lie.md)) |

確認できたこと: 上の各行。確認できなかったこと (または確認が間接的なもの): (i) ベクトル恒等式 (1a.21) は成分の展開で確認し、一般の証明は書いていない。(ii) §8 の「同じ作用の別の書き方」(1a.46)(1a.51) は、$P_t = -H$ と $\tilde H = -p_s$ の代数の確認で、作用の変分 (境界条件を含む) の厳密な議論は行っていない。(iii) Bell (1955) は Appendix 1 だけを確認した (§8.6)。小磯の OHO テキストは I-1〜I-7 ページだけを目で照合した (参考文献)。Stupakov と Penn は未入手。

## 12. 未確認・限界・使った条件

- **相対論的な扱い**: 運動方程式は厳密に相対論的 (§3〜§5、§8)。展開 (§9) は、横方向の運動量が小さいこと ($p_x^2 + p_y^2 \ll (1 + \delta)^2$) に基づく 2 次の近似。$\delta$ は展開していない。
- **静的な場の仮定**: $\phi = 0$ (または時間によらない)、$\vec A$ が時間によらない場で、$H$ が保存し $\delta$ が一定になる。加速空洞のような時間変化する場、ビーム自身の電磁場 (空間電荷、ウェイク)、放射 (放射の反作用は、ハミルトン系ではない) は扱っていない。
- **座標系の仮定**: 水平面内の曲率だけ (ねじれ $\tau = 0$)、進行方向の磁場なし ($A_X = A_Y = 0$)。ソレノイド ($A_X, A_Y \neq 0$) や、垂直偏向・skew の曲線座標は扱っていない。
- **ゲージ**: 偏向磁石の $A_s$ (式 (1a.60b)) と四極の $A_s$ (式 (1a.65)) は、$Y$ と $s$ だけの関数を足す自由度 (ゲージ) の中の 1 つの選び方。別のゲージは、同じ磁場を与え、運動方程式 ($x'$、$p_x'$ の関係) は変わらない ($p_x$ の値は、正準運動量なのでゲージ変換で変わる)。
- **展開の次数**: (1a.59) の 3 次の項の影響は §9.6 の数値例のとおり。SAD は 2 次の打ち切りの代わりに $H_2 + \Delta H$ の分割を使う (04 章)。
- **偏向磁石の端部・多極**: 端部 (フリンジ)、複合機能、6 極以上の多極の $A_s$ の展開は、各要素の章で書く。
- **$z$ の運動と $\delta$**: 縦方向の運動 ($z$、$dv$) は [01_variables_and_conventions.md](01_variables_and_conventions.md) と四極の章。この章では $\delta$ を定数として扱う。

## 参考文献

- SAD のヘルプ (CERN 版 SADHelp) の dynamics の節: Independent variable、Lagrangean、Hamiltonian ([01_variables_and_conventions.md](01_variables_and_conventions.md))。
- M. R. Bell, "Non-linear equations of motion in the synchrotron", A.E.R.E. T/M 125 (1955)。**Appendix 1 (原本 9〜10 ページ) を確認した**: 独立変数を $t$ から別の座標に替えたときの新しいハミルトニアン ($M = Lt'$、$P = -H$、$R = -p_r$)。本章 §8 と同じ手続きで、対応は §8.6。Appendix 1 以外は未確認。
- 小磯晴代、OHO (高エネルギー加速器セミナー) 1991 年のテキスト「ビーム貯蔵リングの基礎」(I 部、全 22 ページ)。**I-1〜I-7 ページを確認した** (8 ページ以降は、本章では未確認):
  - I-1 (§2.1): ラグランジアンとハミルトニアン、正準運動量 $p_i = \gamma m\,dx_i/dt + eA_i$ が力学的運動量と異なること、正準方程式。本章 §3〜§5 と同じ内容。
  - I-2〜I-3 (§2.2): Frenet–Serret 座標 (主法線 $\mathbf e_x$ はリングの外向きが正)、微分の公式、母関数 $G = -\mathbf p\cdot\mathbf r$、新しい運動量、$A_s$ の定義、$H$ ($s$ 成分に $1/(1+\kappa_xx+\kappa_yy)^2$)。本章 §7 と同じ内容。
  - I-3 (§2.3): 座標 $q_j$ ($dq_j/dt \ne 0$、つまり $t$ の単調関数) を独立変数にして作用を書き直し、$H_a = -p_s$、$H = E$ を $p_s$ について解いた形、$s$ を独立変数とする正準方程式 (式 (2-3-8)。$dt/ds = \partial H_a/\partial(-E)$ を含む)。本章 §8 と同じ内容。単調性の条件は本章 §8 では書いていなかった。
  - I-4 (§2.3〜§2.4): 運動量とハミルトニアンの規格化 (力学的全運動量 $p$ で割る。$p = 1$ のとき以外は正準変換ではないと書かれている)、$A_s$ から磁場を作る式 ($B_x = (1/h)(\partial A_s/\partial y - \partial A_y/\partial s)$ など。$h = 1 + \kappa_xx + \kappa_yy$)。本章の (1a.60a) と同じ。
  - I-5 (§2.4): 偏向・四極・スキュー四極・六極・ソレノイドの $A_s$ (四極は $\frac{e}{p_0}A_s = \frac12K(y^2 - x^2)$。本章の (1a.65) と同じ)。
  - I-6〜I-7 (§2.5〜§2.7): $H$ の 3 次までの展開と、ディスパージョンの方程式 ($D_x' = D_{px}$、$D_{px}' = -(\kappa_x^2 + K)D_x + \kappa_x$、$L = N = 0$)。本章 (1a.64) と同じ符号。
  - **本章との違い**: OHO は $\delta \equiv (p - p_0)/p$ (本章・SAD は $(p - P_0)/P_0$) で、運動量を設計運動量 $P_0$ ではなく粒子自身の全運動量 $p$ で規格化する。$\tau$ の符号の約束 (OHO の (2-2-3) は $\mathbf e_x' = -\tau\mathbf e_y + \kappa\mathbf e_s$、$\mathbf e_y' = \tau\mathbf e_x$) は、本章 (1a.38) の下の注と同じ。
  - 本章の式は、これらのページと目で照合した (sympy ではない)。
- G. Stupakov and G. Penn, "Classical Mechanics and Electromagnetism in Accelerator Physics": 入手しておらず、内容は書いていない (未確認)。
