# 18. 非線形への橋渡し: 要素の写像 → 1 周の写像 → 正規形 → ADTS・RDT

**この章で分かること**: 02〜11 章で作った**要素ごとのハミルトニアンと写像**が、リング全体の**非線形な共鳴・振幅依存チューンシフト (ADTS) の解析**にどうつながるかを、途中式つきで導く。(1) 要素の kick・運動学的な項のハミルトニアン (§3.1)、(2) 正規化座標・複素変数 $\zeta^\pm$・ポアソン括弧 (§3.2〜3.3)、(3) 要素の写像の合成が 1 周の写像になること (順序つき指数の展開) と、共鳴駆動項 (RDT) $h_{jklmp}$ の導出 (§3.4〜3.6)、(4) 1 次の摂動 (8 極・運動学・端部フリンジ) と 2 次の摂動 (6 極の対) による ADTS の導出 (§3.7〜3.8)、(5) クロマティシティ、1 本の共鳴の近くのハミルトニアン、強さと振幅の次数の勘定 (§3.9〜3.11)、(6) 周波数マップ・局所運動量アクセプタンス (LMA) の式 (§3.12)、(7) 追跡の写像から係数を取り出す方法 (§4)。

**範囲**: §3 と §4 の式は本章で導出し、sympy ([nl_01_lie_rdt_adts_derivations.py](scripts/nl_01_lie_rdt_adts_derivations.py)、出力 [nl_01_lie_rdt_adts_derivations.out.txt](outputs/nl_01_lie_rdt_adts_derivations.out.txt)) と、薄い 6 極・8 極の列の写像の追跡 ([nl_02_adts_formula_vs_tracking.py](scripts/nl_02_adts_formula_vs_tracking.py)、出力 [nl_02_adts_formula_vs_tracking.out.txt](outputs/nl_02_adts_formula_vs_tracking.out.txt)) で確かめた。途中式の補足 (ヤコビ恒等式と BCH、起点の独立性、共鳴近傍の写像、Hann 窓の裾) は [nl_03_gap_fills.py](scripts/nl_03_gap_fills.py) (出力 [nl_03_gap_fills.out.txt](outputs/nl_03_gap_fills.out.txt))。nl_02 は pysad を使わない自己完結の模型で、**pysad の `resdrive`・`rdt2nd` の数値との突き合わせは本章では取り直していない**。pysad の解析と追跡の比較 (§6) は [docs/24](../docs/24_adts_rdt_primer.md)・[docs/25](../docs/25_resonance_terms_adts.md)・[docs/26](../docs/26_frequency_map_lma.md) からの転記。導出していないもの (偏向磁石の 3 次の生成関数の解析式、厚い要素の閉形式、$\delta\ne0$ の ADTS、$J^2$ の項など) は §8 に分けて書く。

**前提**: 要素の章 (02〜11)、[C_poisson_lie.md](C_poisson_lie.md) (ポアソン括弧とリー演算子)、[03_symplectic_integrators.md](03_symplectic_integrators.md) (BCH)。

**読む順**: §1 (結論) → §2 (流れ) → §3 (導出。§3.4 が中心) → §6 (比較の要約) → §8 (限界)。

**記法**: ポアソン括弧は C_poisson_lie.md の約束 $[f, g] = \sum_i(\partial_{q_i}f\,\partial_{p_i}g - \partial_{p_i}f\,\partial_{q_i}g)$、$:f:\,g = [f, g]$。要素 (長さ $ds$、ハミルトニアン $H$) の写像は $\exp(-ds:H:)$、つまり $f = -H\,ds$ の $\exp(:f:)$ ($z\circ\Phi = \exp(:f:)\,z$ の形)。docs/24・25 は括弧を $\{\ ,\ \}$ と書くが、同じ約束 ($\{\zeta, \bar\zeta\} = 2i$) である。

## 1. 結論

- **要素の写像の合成が 1 周の写像**: リングの 1 周は、02〜11 章の各要素の写像を順に合成したもの $\mathcal M = \mathcal E_N\circ\cdots\circ\mathcal E_1$。どの $\mathcal E_i$ もシンプレクティックなので $\mathcal M$ もシンプレクティック。**閉軌道のまわりの 1 次** (微分) が 12 章の 1 周の行列 $M$。**2 次以上**が非線形な部分で、これが共鳴と ADTS を決める。
- **1 周の写像のリー表現** (§3.4): 線形部分の流れで座標を取り直す (相互作用表示) と、$g\circ\Phi_{\rm int} = e^{:F:}g$、$F = -h_3 - h_4 - \cdots$、$h_3 = \int H_3\,ds$、$h_4 = \int H_4\,ds + \tfrac12\sum_{s_l > s_e}[H_3(s_l), H_3(s_e)]\,ds_l\,ds_e$ ((18.18))。1 次は素子ごとの足し算、2 次は並びに依る対の和。
- **共鳴駆動項 (RDT)** $h_{jklmp}$ (§3.5) は、$H_3$ を複素変数 $\zeta^\pm = \hat x \mp i\hat p = \sqrt{2J}\,e^{\pm i\phi}$ の単項式 $(\zeta_x^+)^j(\zeta_x^-)^k(\zeta_y^+)^l(\zeta_y^-)^m\delta^p$ に展開した係数。薄い 6 極の式は (18.19)〜(18.22)。単項式の位相 $e^{i[(j-k)\phi_x + (l-m)\phi_y]}$ が 1 周ごとにそろう条件 $(j-k)\nu_x + (l-m)\nu_y = $ 整数が共鳴 (§3.6)。
- **ADTS** (§3.7〜3.8): 平均化したハミルトニアン $K(J)$ から $\nu = \nu_0 + \frac1{2\pi}\partial K/\partial J$。対角の 4 次の項 $h_{22000}$ などから $\partial\nu_x/\partial J_x = 4h_{22000}/\pi$ ((18.24))。**1 次の摂動** (8 極・運動学的な項・端部フリンジ) は $h_4$ の対角にそのまま入る (§3.7)。**6 極は 3 次なので 1 次では ADTS を作らず**、2 次 (6 極 × 6 極) で作る: 正規形の変換 $g_q = F_q/(e^{i\mathbf m\cdot\boldsymbol\mu} - 1)$ の後に $K = h_4 - \frac12[g, Rg]$ の対角が残り ((18.32))、1 次元では (18.35)。1 個の薄い 6 極なら $\partial\nu_x/\partial J_x = -\frac{a^2}{64\pi}[\cot3\pi\nu + 3\cot\pi\nu]$ ($a = K_2L\,\beta^{3/2}$、(18.36))。強さの 2 乗に比例し、3$\nu_x$・$\nu_x$ が整数に近いと発散する。
- **追跡の写像から取り出せる** (§4): 1 周の追跡の写像のテイラー展開の 2 次の係数から $f_3$、3 次の係数から $f_4$ が一意に決まる ((18.46)・(18.47))。pysad は `rdttrack.py`・`nform.py`・`bendnl.py` でこれを使う。基準値の ADTS は追跡 + NAFF + 多項式フィット (`adtsdirect.py`)。
- **確認** (§6): 本章の式 (18.32) と 3 個の 6 極の列の追跡は $10^{-4}$ 前後で一致 (nl_02)。pysad の解析と追跡の比較は、6 極 2 次・8 極 1 次が $10^{-3}$ 前後、全部入りの合計が追跡の直接のフィットと 0.01〜0.26% (6 極 ×1)。残りは 6 極に依らない 4 次の定数 (追跡 − 解析: 大きなリング A で $7.0, 6.9, 17$ [1/(m rad)]。原因は四極の本体の積分器の端数と推定、未解決)。

## 2. 流れ

| 段 | 内容 | 章・節 |
|---|---|---|
| 1 | 要素ごとの $H$ と厳密な流れ | 02〜11 章、§3.1 |
| 2 | 要素の写像を 1 周に合成。COD のまわりのテイラー展開 | 12 章 (1 次)、§3.4、§4 (2 次以上) |
| 3 | 線形部分 (Twiss、チューン) の正規化 $\hat x = x/\sqrt\beta$、$\hat p = (\alpha x + \beta p_x)/\sqrt\beta$。複素変数 $\zeta^\pm$ | 12 章、§3.2 |
| 4 | 各素子の非線形項 ($H_3$、$H_4$) を相互作用表示で書く。6 極: $H_3 = \frac{K_2}{6}(x^3 - 3xy^2)$、8 極: $H_4 = \frac{K_3}{24}(x^4 - 6x^2y^2 + y^4)$ ([07_thin_multipole.md](07_thin_multipole.md))。偏向磁石・四極の端部・運動学的な項 $(p_x^2 + p_y^2)^2/8$ も同じ形の多項式の項 | 05・06・07 章、§3.1 |
| 5 | **1 次**: $h_3 = \int H_3\,ds$ を 1 周で足し (厚い 6 極は閉じた式)、$h_{jklmp}$ にする | §3.5 |
| 6 | **2 次**: 対の和で $h_4$。対角の $h_{22000}$ など → ADTS。非対角 → 4 次の共鳴 | §3.4、§3.7〜3.8 |
| 7 | **正規形**: 非共鳴項を $e^{:g:}$ で消し、平均化した $K(J)$ | §3.8 |
| 8 | 追跡 + FFT (NAFF) で $\partial\nu/\partial J$ を直接求め、解析式と比べる | §3.12、`adtsdirect.py` |

## 3. 要素の $H$ から $h_{jklmp}$・ADTS へ

### 3.1 要素の非線形ハミルトニアン

設計軌道のまわりの横方向 $(x, p_x, y, p_y)$ ($\delta = 0$) の、$s$ を独立変数とするハミルトニアンの非線形な部分を、次数ごとに書く。

**多極磁石** (07 章): 薄い多極の kick を生む 2 次元のポテンシャルは、$c_n = K_n$ (通常の多極) として

$$
H_{n+1} = \frac{1}{(n+1)!}\,\mathrm{Re}\left[c_n\,(x + iy)^{n+1}\right]
\tag{18.1}
$$

$n = 2$ (6 極) では $\mathrm{Re}(x + iy)^3 = x^3 - 3xy^2$ なので $H_3 = \frac{K_2}{6}(x^3 - 3xy^2)$、$n = 3$ (8 極) では $\mathrm{Re}(x + iy)^4 = x^4 - 6x^2y^2 + y^4$ なので $H_4 = \frac{K_3}{24}(x^4 - 6x^2y^2 + y^4)$ (sympy で $\mathrm{Re}$ の展開との差 0)。正準方程式 $p' = -\partial H/\partial q$ から kick は

$$
\Delta p_x = -\frac{K_2}{2}(x^2 - y^2),\quad \Delta p_y = K_2\,xy\qquad(6\text{ 極、薄いレンズでは }K_2 \to K_2L)
\tag{18.2}
$$

$$
\Delta p_x = -\frac{K_3}{6}(x^3 - 3xy^2),\quad \Delta p_y = \frac{K_3}{6}(3x^2y - y^3)\qquad(8\text{ 極})
$$

((18.2) は pysad の追跡の kick の符号と同じ約束 (docs/25 §2)。sympy で $-\partial H_3/\partial x = -\frac{K_2}{2}(x - y)(x + y)$、$-\partial H_3/\partial y = K_2xy$ を確認。)

**運動学的な項**: ドリフトのハミルトニアン (02 章) は $\delta = 0$ で $-\sqrt{1 - p_x^2 - p_y^2}$。$p^2 = p_x^2 + p_y^2$ として

$$
-\sqrt{1 - p^2} = -1 + \frac{p^2}{2} + \frac{p^4}{8} + \frac{p^6}{16} + \cdots
\tag{18.3}
$$

($p^2/2$ が線形のドリフト、$p^4/8$ が最初の非線形項。sympy で級数を確認。) したがって 4 次の項は $H_4^{\rm kin} = (p_x^2 + p_y^2)^2/8$。これは全周にわたって存在する (6 極の位置に限らない)。

**四極の端部フリンジ**: 06 章の生成関数 $f = -\frac{K}{12p}[p_x(x^3 + 3xy^2) - p_y(3x^2y + y^3)]$ ($p \to 1$) は、端での 1 回の kick $e^{:f:} = e^{-:H_f:}$ に対応し、

$$
H_f = \sigma\frac{K_1}{12}\left[(x^3 + 3xy^2)p_x - (y^3 + 3x^2y)p_y\right]\delta(s - s_{\rm edge}),\quad\sigma = \begin{cases}+1 & \text{入口}\\ -1 & \text{出口}\end{cases}
\tag{18.4}
$$

(4 次の項: 座標の 3 乗 × 運動量。$\sigma$ の符号は docs/25 §2 の記述による。ここでは (18.4) を受け取り、ADTS への寄与だけを §3.7 で導く。)

**分散のある場所**: $\delta \ne 0$ では閉軌道が $x_\delta = \eta\delta$ だけずれるので、$x \to x + \eta\delta$ と置いた展開が $\delta$ の項を作る (§3.5、§3.9)。

### 3.2 正規化座標と複素変数

Twiss パラメータ ($\alpha$、$\beta$、$\gamma = (1 + \alpha^2)/\beta$) で正規化した座標を

$$
\hat x = \frac{x}{\sqrt\beta},\qquad \hat p = \frac{\alpha x + \beta p_x}{\sqrt\beta}
\tag{18.5}
$$

とする。$J = (\hat x^2 + \hat p^2)/2 = (\gamma x^2 + 2\alpha xp_x + \beta p_x^2)/2$ (単位 m rad)、$\hat x = \sqrt{2J}\cos\phi$、$\hat p = -\sqrt{2J}\sin\phi$ (pysad と SAD の向きの規則)。複素変数

$$
\zeta^+ \equiv \zeta = \hat x - i\hat p = \sqrt{2J}\,e^{i\phi},\qquad \zeta^- \equiv \bar\zeta = \hat x + i\hat p = \sqrt{2J}\,e^{-i\phi}
\tag{18.6}
$$

($y$ 方向は同様に $v = \hat y - i\hat q$。) 線形な運動では $\phi$ が 1 周で $2\pi\nu$ 進むので、$\zeta \to e^{2\pi i\nu}\zeta$。逆に $\hat x = (\zeta + \bar\zeta)/2$、$\hat p = i(\zeta - \bar\zeta)/2$。

**位置 $s$ での座標**: リングの起点 ($s = 0$) の $\zeta$ を基準にする。位置 $s$ までの位相進み $\mu(s)$ ($\mu' = 1/\beta$) で、その位置の局所の正規化座標は $\zeta_{\rm loc} = e^{i\mu}\zeta$ ($\hat x_{\rm loc} = \hat x\cos\mu + \hat p\sin\mu$、$\hat p_{\rm loc} = -\hat x\sin\mu + \hat p\cos\mu$)。実座標は $x = \sqrt\beta\,\mathrm{Re}\,\zeta_{\rm loc}$ なので、$w(s) \equiv \sqrt{\beta(s)}\,e^{i\mu(s)}$ として

$$
x = \tfrac12\left(w\,\zeta + \bar w\,\bar\zeta\right),\qquad w = \sqrt\beta\,e^{i\mu}
\tag{18.7}
$$

運動量は (18.5) を逆に解いて $p_x = (\hat p_{\rm loc} - \alpha\hat x_{\rm loc})/\sqrt\beta$。$w' = dw/ds$ を計算する ($\beta' = -2\alpha$、$\mu' = 1/\beta$):

$$
w' = \left(\frac{\beta'}{2\sqrt\beta} + i\sqrt\beta\,\mu'\right)e^{i\mu} = \frac{-\alpha + i}{\sqrt\beta}\,e^{i\mu}
$$

$\mathrm{Re}(w'\zeta) = \frac1{\sqrt\beta}\mathrm{Re}[(-\alpha + i)\zeta_{\rm loc}]$。$\zeta_{\rm loc} = \hat x_{\rm loc} - i\hat p_{\rm loc}$ なので $\mathrm{Re}[(-\alpha + i)(\hat x_{\rm loc} - i\hat p_{\rm loc})] = -\alpha\hat x_{\rm loc} + \hat p_{\rm loc}$。したがって

$$
p_x = \tfrac12\left(w'\,\zeta + \bar w'\,\bar\zeta\right),\qquad |w'|^2 = \frac{1 + \alpha^2}{\beta} = \gamma
\tag{18.8}
$$

((18.7)(18.8) と $dw/ds$ を sympy で確認: [nl_01](outputs/nl_01_lie_rdt_adts_derivations.out.txt) (c)。) $|w|^2 = \beta$ も使う: $w\bar w = \beta$。

### 3.3 ポアソン括弧: $[\zeta, \bar\zeta] = 2i$ と作用・角変数

(18.6) から $[\zeta, \bar\zeta] = [\hat x - i\hat p,\ \hat x + i\hat p] = i[\hat x, \hat p] - i[\hat p, \hat x] = 2i[\hat x, \hat p] = 2i$ ($[\hat x, \hat p] = 1$)、$[\zeta, \zeta] = [\bar\zeta, \bar\zeta] = 0$。$\zeta$ と $\bar\zeta$ の関数 $f$、$g$ には連鎖律で

$$
[f, g] = \frac{\partial f}{\partial\zeta}\frac{\partial g}{\partial\bar\zeta}[\zeta, \bar\zeta] + \frac{\partial f}{\partial\bar\zeta}\frac{\partial g}{\partial\zeta}[\bar\zeta, \zeta] = 2i\left(\frac{\partial f}{\partial\zeta}\frac{\partial g}{\partial\bar\zeta} - \frac{\partial f}{\partial\bar\zeta}\frac{\partial g}{\partial\zeta}\right)
\tag{18.9}
$$

($y$ 方向の $v$、$\bar v$ があれば、同じ形の項を足す。) 単項式 $f = \zeta^a\bar\zeta^b$、$g = \zeta^c\bar\zeta^d$ では $f_\zeta = a\zeta^{a-1}\bar\zeta^b$、$f_{\bar\zeta} = b\zeta^a\bar\zeta^{b-1}$、$g_\zeta = c\zeta^{c-1}\bar\zeta^d$、$g_{\bar\zeta} = d\zeta^c\bar\zeta^{d-1}$ なので、$f_\zeta g_{\bar\zeta} = ad\,\zeta^{a+c-1}\bar\zeta^{b+d-1}$、$f_{\bar\zeta}g_\zeta = bc\,\zeta^{a+c-1}\bar\zeta^{b+d-1}$。よって

$$
[\zeta^a\bar\zeta^b,\ \zeta^c\bar\zeta^d] = 2i\,(ad - bc)\,\zeta^{a+c-1}\bar\zeta^{b+d-1}
\tag{18.10}
$$

($a, b, c, d = 0\ldots3$ の 256 通りを sympy で $\hat x, \hat p$ の括弧から直接計算して確認。符号は $[\hat x, \hat p] = +1$ の約束 ((C.1)) から決まる。**括弧の符号の約束が逆の文献の式とは、$\zeta^\pm$ の定義 (あるいは $\hat p$ の向き) が変わるので直接は比べられない**: 比べるときは $[\zeta, \bar\zeta]$ の符号を最初に合わせる。)

**作用・角変数**: $\zeta = \sqrt{2J}e^{i\phi}$ とすると $[\zeta, \bar\zeta]_{(\phi,J)} = \partial_\phi\zeta\,\partial_J\bar\zeta - \partial_J\zeta\,\partial_\phi\bar\zeta = i\zeta\cdot\frac{\bar\zeta}{2J} - \frac{\zeta}{2J}\cdot(-i\bar\zeta) = 2i$ で (18.9) と同じ値なので、$(\phi, J)$ は ($\phi$ を座標、$J$ を運動量とする) 正準対 ($[\hat x, \hat p]_{(\phi,J)} = 1$ も sympy で確認)。したがって

$$
[F, J] = \frac{\partial F}{\partial\phi},\qquad [F, \phi] = -\frac{\partial F}{\partial J}
\tag{18.11}
$$

$e^{:F:}$ を 1 回作用させた (1 次の) 変化は $\Delta J = [F, J] = \partial F/\partial\phi$、$\Delta\phi = [F, \phi] = -\partial F/\partial J$。$F = -h$ なら

$$
\Delta J = -\frac{\partial h}{\partial\phi},\qquad \Delta\phi = +\frac{\partial h}{\partial J}
\tag{18.12}
$$

(例: 薄い 6 極 1 個の kick $\Delta\hat p = -\frac a2\hat x^2$ ($a = K_2L\,\beta^{3/2}$) の 1 次の $\Delta J = \hat p\,\Delta\hat p$ は $\frac a2(2J)^{3/2}\sin\phi\cos^2\phi$、$h = \frac a6\hat x^3$ から $-\partial h/\partial\phi$ も同じ。sympy (d) の最終行。)

**回転**: 線形な 1 周 ($\zeta \to e^{i\mu}\zeta$、$\mu = 2\pi\nu$) を関数に作用させる演算子を $R$ ($(Rg)(\zeta) = g(e^{i\mu_x}\zeta, e^{i\mu_y}v)$) と書く。単項式は $R$ の固有関数で

$$
R\left(\zeta^j\bar\zeta^kv^l\bar v^m\right) = e^{i[(j-k)\mu_x + (l-m)\mu_y]}\,\zeta^j\bar\zeta^kv^l\bar v^m
\tag{18.13}
$$

($[J, \zeta] = J_\phi\zeta_J - J_J\zeta_\phi = -i\zeta$ から $R = e^{-\mu:J:}$、つまり $H_0 = \nu J$ ($\theta = 2\pi s/C$ あたり) の流れ $e^{-2\pi:H_0:}$。) $R$ はシンプレクティックな写像による合成なので $R[a, b] = [Ra, Rb]$、すなわち $R\,:g:\,R^{-1} = :Rg:$ が成り立つ (§3.8 で使う)。

### 3.4 要素の写像の合成と 1 周の写像

**相互作用表示**: 位置 $s$ までの線形な流れ $\rho(s)$ (12 章の行列) で座標を取り直し、$\tilde z = \rho(s)^{-1}z$ とする。線形部分が消え、非線形部分のハミルトニアンは、$\tilde z$ (すなわち起点の $\zeta$) の関数

$$
H_{\rm int}(\zeta, s) = H_{\rm nl}\bigl(x(s), p_x(s), \ldots\bigr)\Big|_{(18.7),\,(18.8)}
\tag{18.14}
$$

になる (非線形項の $x$、$p_x$ に (18.7)(18.8) を代入する。例: 6 極は $H_3 = \frac{K_2}{6}(x^3 - 3xy^2)$ の $x$ に $\frac12(w\zeta + \bar w\bar\zeta)$ を入れる)。1 周の写像は、この相互作用表示の写像 $\Phi_{\rm int}$ のあとに全周の線形な回転 $\rho(C)$ をかけたもの。

**合成の順序**: 要素 $i$ ($i = 1$ が最初に通る) の相互作用表示の写像を $\Phi_i$、$z\circ\Phi_i = e^{:F_i:}z$、$F_i = -H_i\,ds_i$ とする。任意の関数 $g$ について $g\circ\Phi_i = e^{:F_i:}g$。2 つの要素 (1 が先、2 が後) の合成は、$u \equiv g\circ\Phi_2 = e^{:F_2:}g$ に $\Phi_1$ を合成して

$$
g\circ\Phi_2\circ\Phi_1 = u\circ\Phi_1 = e^{:F_1:}u = e^{:F_1:}e^{:F_2:}g
\tag{18.15}
$$

**早い要素が左**に来る。具体例で確かめる: $H_1 = c_1q^3$ (kick、先)、$H_2 = c_2p^3$ (後) の点の写像は $(q, p)\to(q,\ p - 3c_1q^2)\to(q + 3c_2(p - 3c_1q^2)^2,\ p - 3c_1q^2)$ ($e^{-:H:}$ の級数が 2 項で止まる ((C.12) と同じ形))。次の $F$ の $e^{:F:}q$、$e^{:F:}p$ と 3 次まで一致する (sympy (f)):

$$
F = F_1 + F_2 + \tfrac12[F_1, F_2] = -(H_1 + H_2) + \tfrac12[H_1, H_2]
\tag{18.16}
$$

括弧の符号を逆にした $F$ では、$q$ の 3 次の項が $-18c_1c_2\,p\,q^2$、$p$ の 3 次の項が $+18c_1c_2\,p^2q$ だけ違ってしまう (sympy)。

**(18.16) の導出**: 2 つの指数の積の BCH 公式 $e^Ae^B = e^{A + B + \frac12[A, B] + \cdots}$ ($[A, B] = AB - BA$ は演算子の交換子) を、$A = :F_1:$、$B = :F_2:$ に使う。まず演算子の交換子を、括弧 $[F_1, F_2]$ の演算子 $:[F_1, F_2]:$ に結びつける (手順 1〜2)。

1. 演算子 $:a:$ は $:a:\,g = [a, g]$。交換子を関数 $g$ に作用させると、$[:a:, :b:]\,g = [a, [b, g]] - [b, [a, g]]$。
2. ヤコビ恒等式 $[a, [b, g]] + [b, [g, a]] + [g, [a, b]] = 0$ に、$[g, a] = -[a, g]$、$[g, [a, b]] = -[[a, b], g]$ を入れると $[a, [b, g]] - [b, [a, g]] - [[a, b], g] = 0$。したがって $[:a:, :b:]\,g = [[a, b], g]$、すなわち
   $$[:a:, :b:] = :[a, b]:$$
   (一般の 3 次式 $a, b$ と $g$ に (18.9) の括弧を使って sympy で確認: [nl_03_gap_fills.out.txt](outputs/nl_03_gap_fills.out.txt) (a))。
3. 手順 2 を BCH の $\frac12[A, B]$ に使うと、2 つの指数の積は $e^{:F_1:}e^{:F_2:} = \exp\{:F_1 + F_2 + \frac12[F_1, F_2]:\ +\cdots\}$ で、(18.16) の形になる (残りの項は 3 重の括弧 $[F_1, [F_1, F_2]]$ など。$F_i$ が 3 次式ならこれは 5 次式)。
4. 3 つの指数の積は、BCH を 2 回使う。1 回目: $e^ae^b = e^{a + b + \frac12[a, b]}$ (3 重の括弧は 5 次式なので落とす)。2 回目: $e^{a + b + \frac12[a, b]}\,e^c = \exp\bigl\{(a + b + \tfrac12[a, b]) + c + \tfrac12[a + b + \tfrac12[a, b],\ c]\bigr\}$。最後の括弧を展開すると $\frac12[a, c] + \frac12[b, c] + \frac14[[a, b], c]$ で、最後の項は 5 次式。したがって
   $$e^ae^be^c = \exp\Bigl\{a + b + c + \tfrac12\bigl([a,b] + [a,c] + [b,c]\bigr) + (5\text{ 次以上})\Bigr\}$$
   (3 次式 $a, b, c$ について、$\zeta$ への作用を 3 次まで比べて確認: nl_03 (b))。

**$N$ 個の要素**: 3 次式どうしの括弧は 4 次式 ($3 + 3 - 2 = 4$)、3 重の括弧は 5 次式なので、4 次までは

$$
F = \sum_iF_i + \frac12\sum_{e < l}[F_e, F_l]
\tag{18.17}
$$

$F_i = -H_i\,ds_i$ ($H_i = H_{3,i} + H_{4,i} + \cdots$) を代入し、$[F_e, F_l] = [H_e, H_l]\,ds_e\,ds_l = -[H_l, H_e]\,ds_e\,ds_l$ を使うと、3 次: $F_3 = -\sum_iH_{3,i}ds_i$、4 次: $F_4 = -\sum_iH_{4,i}ds_i - \frac12\sum_{l > e}[H_{3,l}, H_{3,e}]\,ds_l\,ds_e$。$F = -h_3 - h_4 - \cdots$ と置くと

$$
h_3 = \int H_3\,ds,\qquad
h_4 = \int H_4\,ds + \tfrac12\sum_{s_l > s_e}[H_3(s_l),\,H_3(s_e)]\,ds_l\,ds_e
\tag{18.18}
$$

連続極限では二重積分 $\iint_{s_l > s_e}$。**同じ薄い要素の中の対は 0** ($[H, H] = 0$)。厚い要素の内部の対は $\iint_{s_1 > s_2\in{\rm elem}}$ で、ドリフト (6 極の内部) では閉じた式で積分できる (docs/25 §10.2。本章では導出していない)。

最後に、$g\circ\mathcal M = (g\circ\rho(C))\circ\Phi_{\rm int} = e^{:F:}(Rg)$ ($Rg = g\circ\rho(C)$ は (18.13) の $R$)。以下、非線形部分は $F$、線形部分は $R$ で表す。

**$h_4$ の構成のまとめ** (対の項は並びに依る): 1 次の項 $h_3$ は並びによらず足すだけだが、$h_4$ の対の項は「後の要素 − 前の要素」の順序つきの和なので、**6 極の 2 次の ADTS は 6 極の並び (位相) に依る**。ただし (18.35) のとおり、起点を移したときの $h_4$ の変化は正規形の補正と打ち消し合い、ADTS は起点に依らない (§3.8)。

### 3.5 1 次の共鳴駆動項 $h_{jklmp}$

$H_3$ を $\zeta^\pm$、$v^\pm$、$\delta$ の単項式に展開する。(18.7) と同様に $y = \frac12(w_yv + \bar w_y\bar v)$、$w_y = \sqrt{\beta_y}e^{i\mu_y}$、分散 $x \to x + \eta\delta$。薄い 6 極 (積分強度 $K_2L$) では、$H_3 = \frac{K_2L}{6}\left[(X + \eta\delta)^3 - 3(X + \eta\delta)Y^2\right]\delta(s - s_n)$ ($X = \frac12(w_x\zeta + \bar w_x\bar\zeta)$、$Y = \frac12(w_yv + \bar w_y\bar v)$)。

**$x^3$ の展開** (二項定理、$w_x\bar w_x = \beta_x$):

$$
X^3 = \frac18\left[w_x^3\zeta^3 + 3w_x^2\bar w_x\,\zeta^2\bar\zeta + 3w_x\bar w_x^2\,\zeta\bar\zeta^2 + \bar w_x^3\bar\zeta^3\right]
= \frac18\left[w_x^3\zeta^3 + 3\beta_xw_x\,\zeta^2\bar\zeta + \text{c.c.}\right]
$$

$H_3$ の係数 $K_2L/6$ をかけて、単項式の係数を読む:

$$
h_{30000} = \frac{K_2L}{6}\cdot\frac18\,w_x^3 = \frac{K_2L}{48}\,\beta_x^{3/2}e^{3i\mu_x},\qquad
h_{21000} = \frac{K_2L}{6}\cdot\frac38\,\beta_xw_x = \frac{K_2L}{16}\,\beta_x^{3/2}e^{i\mu_x}
\tag{18.19}
$$

**$xy^2$ の展開**: $Y^2 = \frac14\left(w_y^2v^2 + 2\beta_yv\bar v + \bar w_y^2\bar v^2\right)$ ($w_y\bar w_y = \beta_y$)。$-3XY^2\cdot\frac{K_2L}{6} = -\frac{K_2L}{2}XY^2$ で、$XY^2 = \frac18(w_x\zeta + \bar w_x\bar\zeta)(w_y^2v^2 + 2\beta_yv\bar v + \bar w_y^2\bar v^2)$ なので ($-\frac{K_2L}{2}\cdot\frac18 = -\frac{K_2L}{16}$)

$$
h_{10110} = -\frac{K_2L}{16}\cdot2\beta_yw_x = -\frac{K_2L}{8}\sqrt{\beta_x}\,\beta_y\,e^{i\mu_x},\quad
h_{10020} = -\frac{K_2L}{16}\,w_x\bar w_y^2 = -\frac{K_2L}{16}\sqrt{\beta_x}\,\beta_y\,e^{i(\mu_x - 2\mu_y)},\quad
h_{10200} = -\frac{K_2L}{16}\sqrt{\beta_x}\,\beta_y\,e^{i(\mu_x + 2\mu_y)}
\tag{18.20}
$$

($h_{10020}$: $\zeta\,\bar v^2$ すなわち $(j, k, l, m) = (1, 0, 0, 2)$、$h_{10200}$: $\zeta v^2$ すなわち $(1, 0, 2, 0)$。)

**分散の項**: $(X + \eta\delta)^3$ の $\delta^1$ の項は $3X^2\eta\delta$、$-3(X + \eta\delta)Y^2$ の $\delta^1$ の項は $-3\eta\delta\,Y^2$。$\frac{K_2L}{6}$ をかけて $\frac{K_2L}{2}\eta\delta\,(X^2 - Y^2)$。$X^2 = \frac14(w_x^2\zeta^2 + 2\beta_x\zeta\bar\zeta + \bar w_x^2\bar\zeta^2)$ から

$$
h_{11001} = \frac{K_2L}{2}\eta\cdot\frac{2\beta_x}{4} = \frac{K_2L}{4}\,\eta\beta_x,\quad
h_{00111} = -\frac{K_2L}{4}\,\eta\beta_y,\quad
h_{20001} = \frac{K_2L}{8}\,\eta\beta_x e^{2i\mu_x},\quad
h_{00201} = -\frac{K_2L}{8}\,\eta\beta_y e^{2i\mu_y}
\tag{18.21}
$$

$\delta^2$ の項は $\frac{K_2L}{6}\cdot3X\eta^2\delta^2$ で $h_{10002} = \frac{K_2L}{4}\eta^2\sqrt{\beta_x}\,e^{i\mu_x}$ ((18.19)〜(18.21) すべてを sympy で展開して確認: nl_01 (d))。

リング全体では、**要素 (厚い要素は積分) について足す**:

$$
h_{jklmp} = \sum_{\rm elem}c_{jklmp}\int w_x^j\bar w_x^kw_y^l\bar w_y^m\,\eta^p\,ds,\qquad h_{kjmlp} = \bar h_{jklmp}
\tag{18.22}
$$

($c_{jklmp}$ は $H_3$ の多項式の係数に二項係数と $2^{-(j+k+l+m)}$ をかけたもの。薄い要素は積分を位置での値で置き換える。)

**意味**: (18.12) から、1 周の $J_x$ の変化は $\Delta J_x = -\partial h_3/\partial\phi_x$。例えば $h_{30000}(2J_x)^{3/2}e^{3i\phi_x}$ + c.c. $= 2|h_{30000}|(2J_x)^{3/2}\cos(3\phi_x + \chi)$ から $\Delta J_x = 6|h_{30000}|(2J_x)^{3/2}\sin(3\phi_x + \chi)$。位相空間のゆがみ (三角形) を作る。$h_{11001}$ ($\zeta\bar\zeta\,\delta = 2J_x\delta$) は $J_x$ に比例する項なので**チューンのずれ** (§3.9)、$h_{20001}$ ($\zeta^2\delta$) は $\cos2\phi_x$ に比例する項で、$\delta$ に比例した位相空間の楕円の変形 (ベータ関数の $\delta$ 依存)。

**6 極の寄与の打ち消し**: $h$ は 6 極 1 本ずつの寄与 (複素数、位相 $e^{3i\mu}$ など) の和。位相が回って打ち消し合うので、和は個々の寄与の絶対値の和より桁で小さくなる (FODO の 48 本の 6 極で $|\sum h_{30000}| = 0.7035$ 対 $\sum|h| = 31.861$、docs/24 §4.2)。「6 極の組」の考え方はこれを使う。

### 3.6 共鳴の条件

(18.22) の単項式の位相は $e^{i[(j-k)\phi_x + (l-m)\phi_y]}$。1 周ごとに $\phi_{x,y}$ は $2\pi\nu_{x,y}$ 進むので、$n$ 周後の位相は $e^{in\cdot2\pi[(j-k)\nu_x + (l-m)\nu_y]}$。

- $\mathbf m\cdot\boldsymbol\nu \equiv (j-k)\nu_x + (l-m)\nu_y$ が整数でなければ、周回とともに位相が回り、$n$ について足すと平均 0 (非共鳴)。
- **整数なら**位相が毎周そろい、同じ向きの kick が積算する (共鳴)。正規形では、(18.13) の固有値 $e^{2\pi i\,\mathbf m\cdot\boldsymbol\nu} = 1$ となり、(18.30) の変換 $g_q = F_q/(e^{2\pi i\,\mathbf m\cdot\boldsymbol\nu} - 1)$ の分母が 0 になる (§3.8)。摂動論が破れる。

3 次の項 ($j + k + l + m = 3$) の調和 $\mathbf m = (j-k, l-m)$:

| 単項式 | $\mathbf m$ | 条件 |
|---|---|---|
| $\zeta^3$ ($h_{30000}$) | $(3, 0)$ | $3\nu_x = N$ |
| $\zeta^2\bar\zeta$ ($h_{21000}$)、$\zeta v\bar v$ ($h_{10110}$) | $(1, 0)$ | $\nu_x = N$ (整数共鳴。1 次の位相空間のゆがみは非共鳴でも出る) |
| $\zeta\bar v^2$ ($h_{10020}$) | $(1, -2)$ | $\nu_x - 2\nu_y = N$ |
| $\zeta v^2$ ($h_{10200}$) | $(1, 2)$ | $\nu_x + 2\nu_y = N$ |

4 次の項: $\mathbf m = (2, 0)$ ($\zeta^3\bar\zeta$、$2\nu_x$)、$(4, 0)$ ($\zeta^4$)、$(2, \pm2)$、$(0, \pm2)$、$(0, 4)$ など。対角 $\mathbf m = (0, 0)$ ($\zeta^2\bar\zeta^2$、$\zeta\bar\zeta v\bar v$、$v^2\bar v^2$) は共鳴ではなく ADTS を作る。

**$\nu_x - \nu_y$ そのものが 6 極の 2 次の共鳴で駆動されない理由**: 上下対称な磁石の $H$ は $(y, p_y)$ について偶数次 ($y \to -y$、$p_y \to -p_y$ で不変)。($l + m$) が偶数になる。$\nu_x - \nu_y$ の共鳴は $\mathbf m = (1, -1)$ で、$l - m = -1$ は $l + m$ が奇数でないと作れない。したがって上下対称な磁石の列は $\nu_x - \nu_y$ を駆動せず、4 次の $2\nu_x - 2\nu_y$ ($\mathbf m = (2, -2)$) などが最低次になる (docs/24 §3.4 の表と同じ結論)。

### 3.7 1 次の摂動による ADTS: 8 極・運動学的な項・端部フリンジ

(18.12) の $\Delta\phi = \partial h/\partial J$ から、角度に依らない項 ($j = k$、$l = m$) だけが平均化後に残り、チューンを変える。4 次の対角項を

$$
K = h_{22000}(\zeta\bar\zeta)^2 + h_{11110}\,\zeta\bar\zeta\,v\bar v + h_{00220}(v\bar v)^2
= h_{22000}(2J_x)^2 + h_{11110}(2J_x)(2J_y) + h_{00220}(2J_y)^2
\tag{18.23}
$$

($\zeta\bar\zeta = 2J_x$) と書くと、1 周のチューン ($\nu = \nu_0 + \frac1{2\pi}\partial K/\partial J$) は

$$
\nu_x = \nu_{x0} + \frac1{2\pi}\left(8h_{22000}J_x + 4h_{11110}J_y\right),\qquad
\nu_y = \nu_{y0} + \frac1{2\pi}\left(4h_{11110}J_x + 8h_{00220}J_y\right)
$$

$$
\frac{\partial\nu_x}{\partial J_x} = \frac{4h_{22000}}{\pi},\qquad
\frac{\partial\nu_x}{\partial J_y} = \frac{\partial\nu_y}{\partial J_x} = \frac{2h_{11110}}{\pi},\qquad
\frac{\partial\nu_y}{\partial J_y} = \frac{4h_{00220}}{\pi}
\tag{18.24}
$$

($\nu = \nu_0 + \frac1{2\pi}\partial K/\partial J$ は、(18.12) の $\Delta\phi = -\partial F/\partial J = \partial K/\partial J$ ($F = -K$) を 1 周の位相の進み $2\pi\Delta\nu$ とおいたもの。)

**8 極**: $h_{22000}$ は $\zeta^2\bar\zeta^2$ の係数。$X^4 = \frac1{16}(w\zeta + \bar w\bar\zeta)^4$ の $\zeta^2\bar\zeta^2$ の係数は $\frac{\binom42}{16}w^2\bar w^2 = \frac{3}{8}\beta_x^2$、$X^2Y^2$ の $\zeta\bar\zeta\,v\bar v$ の係数は $\frac12\beta_x\cdot\frac12\beta_y = \frac14\beta_x\beta_y$ (各 $X^2$、$Y^2$ の $\zeta\bar\zeta$ の係数が $\frac{2\beta}{4}$)。$H_4 = \frac{K_3L}{24}(X^4 - 6X^2Y^2 + Y^4)$ から

$$
h_{22000} = \frac{K_3L}{24}\cdot\frac38\beta_x^2 = \frac{K_3L\beta_x^2}{64},\quad
h_{11110} = -\frac{K_3L}{24}\cdot6\cdot\frac{\beta_x\beta_y}{4} = -\frac{K_3L\beta_x\beta_y}{16},\quad
h_{00220} = \frac{K_3L\beta_y^2}{64}
$$

(18.24) に入れて ($4/64 = 1/16$、$2/16 = 1/8$)

$$
\frac{\partial\nu_x}{\partial J_x} = \sum\frac{K_3L\,\beta_x^2}{16\pi},\qquad
\frac{\partial\nu_x}{\partial J_y} = -\sum\frac{K_3L\,\beta_x\beta_y}{8\pi},\qquad
\frac{\partial\nu_y}{\partial J_y} = \sum\frac{K_3L\,\beta_y^2}{16\pi}
\tag{18.25}
$$

(docs/24 §5.3 と同じ。sympy で $h$ の係数と (18.25) を確認。nl_02 (C): 8 極 1 個 ($K_3L = 40$) を 6 極の列に足すと、式の変化は上の 3 式と一致し、追跡との相対差は $2\times10^{-4}$、$4\times10^{-4}$、$7\times10^{-6}$。)

**運動学的な項**: $H_4^{\rm kin} = (p_x^2 + p_y^2)^2/8$、(18.8) の $p_x = \frac12(w'\zeta + \bar w'\bar\zeta)$ ($|w'|^2 = \gamma$)。$p_x^4$ の $\zeta^2\bar\zeta^2$ の係数は $\frac{\binom42}{16}|w'|^4 = \frac38\gamma_x^2$、$p_x^2p_y^2$ の $\zeta\bar\zeta v\bar v$ の係数は $\frac{\gamma_x}{2}\cdot\frac{\gamma_y}{2} = \frac{\gamma_x\gamma_y}{4}$。$(p_x^2 + p_y^2)^2 = p_x^4 + 2p_x^2p_y^2 + p_y^4$ なので

$$
h_{22000} = \frac18\cdot\frac38\gamma_x^2 = \frac{3\gamma_x^2}{64},\quad
h_{11110} = \frac18\cdot2\cdot\frac{\gamma_x\gamma_y}{4} = \frac{\gamma_x\gamma_y}{16}
$$

(要素ごとに $\int ds$)。(18.24) から

$$
\frac{\partial\nu_x}{\partial J_x} = \frac{3}{16\pi}\oint\gamma_x^2\,ds,\quad
\frac{\partial\nu_x}{\partial J_y} = \frac{1}{8\pi}\oint\gamma_x\gamma_y\,ds,\quad
\frac{\partial\nu_y}{\partial J_y} = \frac{3}{16\pi}\oint\gamma_y^2\,ds
\tag{18.26}
$$

(正の値: 振幅とともにチューンが上がる。6 極も 8 極もない裸のラティスでも、運動学的な項と端のフリンジで ADTS が出る。FODO の例で 6 極なしの J = 1e-6 m rad のチューンが 7.5e-5 動く: docs/24 §2。)

**端部のフリンジ** (18.4): $x^3p_x$ の $\zeta^2\bar\zeta^2$ の係数を求める。$x^3 = \frac18[w^3\zeta^3 + 3w^2\bar w\zeta^2\bar\zeta + 3w\bar w^2\zeta\bar\zeta^2 + \bar w^3\bar\zeta^3]$、$p_x = \frac12(w'\zeta + \bar w'\bar\zeta)$ から、$\zeta^2\bar\zeta^2$ を作る組は $(3w\bar w^2\zeta\bar\zeta^2)\times(w'\zeta/2)$ と $(3w^2\bar w\zeta^2\bar\zeta)\times(\bar w'\bar\zeta/2)$:

$$
\frac{3}{16}\left(w\bar w^2w' + w^2\bar w\bar w'\right) = \frac{3\beta}{16}\left(\bar ww' + w\bar w'\right) = \frac{3\beta}{16}\cdot2\,\mathrm{Re}(\bar ww')
$$

$\bar ww' = \sqrt\beta\,e^{-i\mu}\cdot\frac{-\alpha + i}{\sqrt\beta}e^{i\mu} = -\alpha + i$ なので $2\mathrm{Re} = -2\alpha$、係数は $-\frac38\alpha\beta$。$H_f$ の $\sigma\frac{K_1}{12}$ をかけて $h_{22000} = -\sigma\frac{K_1\alpha_x\beta_x}{32}$、(18.24) から

$$
\frac{\partial\nu_x}{\partial J_x} = -\frac{1}{8\pi}\sum\sigma K_1\alpha_x\beta_x
\tag{18.27}
$$

同様に $\partial\nu_x/\partial J_y = \frac1{8\pi}\sum\sigma K_1(\alpha_y\beta_x - \alpha_x\beta_y)$、$\partial\nu_y/\partial J_y = \frac1{8\pi}\sum\sigma K_1\alpha_y\beta_y$ ($y$ の項は (18.4) の $-(y^3 + \ldots)p_y$ の符号を含む。sympy で (18.4) から 3 式すべてを展開して確認。$\alpha$、$\beta$ は端での値。docs/25 §2 と同じ。)

**なぜ 6 極の 1 次ではチューンシフトが出ないか**: $\zeta^j\bar\zeta^kv^l\bar v^m$ の角度に依らない項は $j = k$ かつ $l = m$、つまり次数 $j + k + l + m = 2(j + l)$ は偶数。3 次の単項式は奇数次なので $j = k$、$l = m$ を満たせず、**$h_3$ の対角項は無い**。物理の言葉では、6 極の力 $\Delta p \propto x^2$ は $x \to -x$ で向きが変わらず、チューンを変える「$x$ に比例する力」を平均では作らない ((18.12) の $\Delta J = -\partial h_3/\partial\phi$ は位相空間の三角形のゆがみで、$\phi$ について平均すると 0)。ADTS が出るのは、4 次の項 (8 極の 1 次、運動学的な項、端のフリンジ) と、3 次の項の 2 次 (6 極 × 6 極、6 極 × 偏向磁石の 3 次の項)。

### 3.8 2 次の摂動: 正規形と 6 極の ADTS

3 次の項は強さの 1 乗で位相空間をゆがめるだけだが、**ゆがんだ後にもう一度同じ種類の kick を受ける** (2 次) と、$J^2$ に比例するハミルトニアンの項が残る。これを、1 周の写像の正規形で導く。

**1 周の写像**: §3.4 から $g\circ\mathcal M = e^{:F:}(Rg)$、$F = F_3 + F_4 + \cdots$ ($F_3 = -h_3$、$F_4 = -h_4$)。この写像を、シンプレクティックな座標変換 $e^{:g:}$ で共役にする ($M' = e^{:g:}Me^{-:g:}$。共役変換でチューンと ADTS は変わらない):

$$
M' = e^{:g:}\,e^{:F:}\,R\,e^{-:g:}
\tag{18.28}
$$

(演算子の合成。) $R\,:g:\,R^{-1} = :Rg:$ ((18.13) の後) から $R\,e^{-:g:} = e^{-:Rg:}R$。よって

$$
M' = e^{:g:}\,e^{:F:}\,e^{-:Rg:}\,R
\tag{18.29}
$$

**3 次の項を消す**: 3 つの指数の積の BCH 公式 (§3.4。$g$、$F_3$、$Rg$ は 3 次式) を 4 次まで使うと、積は $e^{:d:}R$ で

$$
d = \underbrace{(g + F_3 - Rg)}_{3\text{ 次}} + \underbrace{F_4 + \tfrac12[g, F_3] - \tfrac12[g, Rg] - \tfrac12[F_3, Rg]}_{4\text{ 次}} + \cdots
$$

($g + F + (-Rg)$ の和と、3 組の括弧 $\frac12([g, F_3] + [g, -Rg] + [F_3, -Rg])$。) 3 次の項を 0 にする条件は $g + F_3 - Rg = 0$。単項式 $\zeta^j\bar\zeta^kv^l\bar v^m$ ごとに、(18.13) の固有値 $e^{i\mathbf m\cdot\boldsymbol\mu}$ ($\mathbf m\cdot\boldsymbol\mu = (j-k)\mu_x + (l-m)\mu_y$) を使うと $g_q(1 - e^{i\mathbf m\cdot\boldsymbol\mu}) = -F_q$ から

$$
g_q = \frac{F_q}{e^{i\mathbf m\cdot\boldsymbol\mu} - 1},\qquad (Rg)_q = e^{i\mathbf m\cdot\boldsymbol\mu}g_q = \frac{e^{i\mathbf m\cdot\boldsymbol\mu}F_q}{e^{i\mathbf m\cdot\boldsymbol\mu} - 1}
\tag{18.30}
$$

($F_q = -h_q$。3 次の単項式は $\mathbf m = 0$ を作れないので分母は 0 にならない。**共鳴 $\mathbf m\cdot\boldsymbol\nu \in \mathbb Z$ では分母が 0** で変換が定義できない (§3.6)。)

**4 次の項**: $F_3 = Rg - g$ を代入すると、$\frac12[g, F_3] - \frac12[F_3, Rg] = \frac12[g, Rg - g] - \frac12[Rg - g, Rg] = \frac12[g, Rg] + \frac12[g, Rg] = [g, Rg]$ ($[g, g] = [Rg, Rg] = 0$)。よって $d_4 = F_4 + [g, Rg] - \frac12[g, Rg]$:

$$
d_4 = F_4 + \tfrac12\,[g,\ Rg]
\tag{18.31}
$$

(sympy で一般の 3 次式について (18.30) の $F_3 + g - Rg = 0$ と $d_4 - F_4 = \frac12[g, Rg]$ を確認。) $d_4$ の**対角の項** (角度に依らない) が $-K$ で、$K = -d_4|_{\rm diag}$。$F_4 = -h_4$ なので

$$
K = h_4\big|_{\rm diag} - \tfrac12[g,\ Rg]\big|_{\rm diag}
\tag{18.32}
$$

(非対角の項は 4 次の共鳴を作る。次の次数の変換で消す。ADTS には影響しない。)

**1 次元 ($\zeta, \bar\zeta$ だけ、$\zeta^2\bar\zeta^2$ の係数)**: 3 次の項は調和 $\pm3$ ($\zeta^3$、$\bar\zeta^3$) と $\pm1$ ($\zeta^2\bar\zeta$、$\zeta\bar\zeta^2$)。$g$、$Rg$ は (18.30) から

$$
g_{30} = \frac{F_{30}}{e^{3i\mu} - 1},\ g_{03} = \bar g_{30}\ (= \tfrac{\bar F_{30}}{e^{-3i\mu} - 1}),\ g_{21} = \frac{F_{21}}{e^{i\mu} - 1},\ g_{12} = \bar g_{21},\quad Rg_{30} = e^{3i\mu}g_{30},\ Rg_{03} = e^{-3i\mu}g_{03},\ \ldots
$$

$[g, Rg]$ のうち $\zeta^2\bar\zeta^2$ を作る対を (18.10) で数える。条件は $a + c = b + d = 2$ (4 次の対角): $(3,0)\times(0,3)$: $2i(3\cdot3 - 0) = 18i$、$(0,3)\times(3,0)$: $2i(0 - 3\cdot3) = -18i$、$(2,1)\times(1,2)$: $2i(2\cdot2 - 1\cdot1) = 6i$、$(1,2)\times(2,1)$: $2i(1\cdot1 - 2\cdot2) = -6i$。他の組 (例: $(3,0)\times(2,1)$ は $\zeta^4\bar\zeta$) は非対角。よって

$$
[g, Rg]\big|_{\zeta^2\bar\zeta^2} = 18i\left(g_{30}Rg_{03} - g_{03}Rg_{30}\right) + 6i\left(g_{21}Rg_{12} - g_{12}Rg_{21}\right)
$$

$g_{30}Rg_{03} = |g_{30}|^2e^{-3i\mu}$、$g_{03}Rg_{30} = |g_{30}|^2e^{3i\mu}$ なので括弧内は $-2i\sin3\mu\,|g_{30}|^2$。同様に $g_{21}Rg_{12} - g_{12}Rg_{21} = -2i\sin\mu\,|g_{21}|^2$。したがって $[g, Rg]|_{\zeta^2\bar\zeta^2} = 36\sin3\mu\,|g_{30}|^2 + 12\sin\mu\,|g_{21}|^2$。$|g_{jk}|^2 = |F_{jk}|^2/|e^{im\mu} - 1|^2 = |F_{jk}|^2/(4\sin^2(m\mu/2))$ を使い、$\sin3\mu = 2\sin\frac{3\mu}2\cos\frac{3\mu}2$ で約分して ($36\cdot2/4 = 18$、$12\cdot2/4 = 6$)

$$
\tfrac12[g, Rg]\big|_{\zeta^2\bar\zeta^2} = 9|F_{30}|^2\cot\frac{3\mu}{2} + 3|F_{21}|^2\cot\frac{\mu}{2} = 9|h_{30000}|^2\cot3\pi\nu + 3|h_{21000}|^2\cot\pi\nu
\tag{18.33}
$$

($\mu = 2\pi\nu$。$F_{jk} = -h_{jk}$。数値代入で sympy の結果と $5\times10^{-16}$ で一致。) 1 周の写像の 4 次の対角 $h_4|_{\rm diag}$ は、8 極などの 4 次の項 $\bigl(\int H_4\,ds\bigr)$ と、(18.18) の対の項の対角:

**対の項の対角** (1 次元): $H_3(s)$ の単項式の係数を $c_{30}(s)$ ($\zeta^3$、$= \frac{K_2L}{48}\beta^{3/2}e^{3i\mu}$ の要素ごとの値)、$c_{21}(s)$ とすると、$[H_l, H_e]$ の $\zeta^2\bar\zeta^2$ の係数は、同じ 4 組の対から $18i\left(c_{30}(l)\bar c_{30}(e) - \bar c_{30}(l)c_{30}(e)\right) + 6i\left(c_{21}(l)\bar c_{21}(e) - \bar c_{21}(l)c_{21}(e)\right)$。$z - \bar z = 2i\,\mathrm{Im}\,z$ なので $18i\cdot2i\,\mathrm{Im}[c_{30}(l)\bar c_{30}(e)] = -36\,\mathrm{Im}[\ldots]$、同様に $-12\,\mathrm{Im}[\ldots]$。$\frac12$ をかけて

$$
\tfrac12\sum_{l > e}[H_l, H_e]\Big|_{\zeta^2\bar\zeta^2} = -18\,\mathrm{Im}\,P_3 - 6\,\mathrm{Im}\,P_1,\qquad
P_3 = \sum_{l > e}c_{30}(l)\bar c_{30}(e),\quad P_1 = \sum_{l > e}c_{21}(l)\bar c_{21}(e)
\tag{18.34}
$$

(18.32)〜(18.34) をまとめて、**1 次元の $\zeta^2\bar\zeta^2$ の係数** $K_{2200}$ ($K = K_{2200}(2J_x)^2$) は

$$
K_{2200} = h^{(4)}_{22000} - 18\,\mathrm{Im}\,P_3 - 6\,\mathrm{Im}\,P_1 - 9|h_{30000}|^2\cot3\pi\nu_x - 3|h_{21000}|^2\cot\pi\nu_x,\qquad
\frac{\partial\nu_x}{\partial J_x} = \frac{4K_{2200}}{\pi}
\tag{18.35}
$$

($h^{(4)}_{22000}$: 8 極・運動学・端部フリンジなどの 1 次の項 (§3.7)。$h_{30000} = \sum c_{30}$、$h_{21000} = \sum c_{21}$ は全周の和。)

**起点の取り方に依らないこと**。1 周の写像の起点 (周期の切れ目) を、要素 1 個ぶん先に動かしても ADTS が変わらないことを、(18.35) の $\mathrm{Im}\,P_3$ の項と $\cot$ の項の変化が打ち消し合うことで示す。$c_e \equiv c_{30}(e)$ (要素 $e$ の $\zeta^3$ の係数) は、起点から測った位相 $\mu_e$ を使って $c_e = |c_e|\,e^{3i\mu_e}$ の形で、$c_{21}$ は $e^{i\mu_e}$ の形である ((18.7) で $x\propto w\zeta$、$w\propto e^{i\mu}$ だから、$x^3$ の $\zeta^3$ の係数は $w^3\propto e^{3i\mu}$、$\zeta^2\bar\zeta$ は $w^2\bar w\propto e^{i\mu}$)。

1. **起点を動かしたときの係数**: 起点を最初の要素 (1 とする) の直後に動かすと、すべての位相が $\mu_1$ だけ引かれ ($c_e\to c_e e^{-3i\mu_1}$、全要素に共通の位相因子)、要素 1 は最後に回って、さらに 1 周ぶん $\mu$ ($= 2\pi\nu$) が加わる ($c_1\to c_1e^{-3i\mu_1}e^{3i\mu}$)。$P_3$ ($c$ の双線形で $c_l\bar c_e$ の形)、$|h_{30000}|^2$ は共通の位相因子で変わらないので、**要素 1 の係数だけに $e^{i\theta}$ ($\theta = 3\mu$) を掛け、要素 1 を和の最後に回す**、と考えてよい。
2. **$P_3$ の変化**: $S' \equiv \sum_{e\ne1}c_e$、$P_3' \equiv$ (要素 1 以外のあいだの対の和) とする。動かす前は、要素 1 が先頭なので、$l > e = 1$ の対は $\sum_{l\ne1}c_l\bar c_1 = S'\bar c_1$:
   $$P_3 = P_3' + S'\,\bar c_1.$$
   動かしたあとは、要素 1 が最後 ($l = 1$、$e$ はその他すべて) で、係数が $c_1e^{i\theta}$:
   $$P_3^{\rm new} = P_3' + c_1e^{i\theta}\,\bar S'.$$
   $X \equiv c_1\bar S'$ とおくと、$S'\bar c_1 = \bar X$ なので
   $$\Delta P_3 = Xe^{i\theta} - \bar X,\qquad \mathrm{Im}\,\Delta P_3 = \mathrm{Im}(Xe^{i\theta}) + \mathrm{Im}\,X = \mathrm{Im}\bigl[X(e^{i\theta} + 1)\bigr]\quad(\mathrm{Im}\,\bar X = -\mathrm{Im}\,X).$$
3. **$|h_{30000}|^2$ の変化**: $h_{30000} = S' + c_1$ が $S' + c_1e^{i\theta}$ になる。$|S' + c|^2 = |S'|^2 + |c|^2 + 2\,\mathrm{Re}(\bar S'c)$ なので
   $$\Delta|h_{30000}|^2 = 2\,\mathrm{Re}\bigl[\bar S'c_1(e^{i\theta} - 1)\bigr] = 2\,\mathrm{Re}\bigl[X(e^{i\theta} - 1)\bigr].$$
4. **$\theta/2$ で書き直す**: $e^{i\theta} + 1 = e^{i\theta/2}(e^{i\theta/2} + e^{-i\theta/2}) = 2\cos\frac\theta2\,e^{i\theta/2}$、$e^{i\theta} - 1 = 2i\sin\frac\theta2\,e^{i\theta/2}$。$\mathrm{Re}(iw) = -\mathrm{Im}\,w$ を使うと
   $$\mathrm{Im}\,\Delta P_3 = 2\cos\tfrac\theta2\,\mathrm{Im}\bigl(Xe^{i\theta/2}\bigr),\qquad \Delta|h_{30000}|^2 = -4\sin\tfrac\theta2\,\mathrm{Im}\bigl(Xe^{i\theta/2}\bigr).$$
5. **$K_{2200}$ の変化**: (18.35) の係数 $-18$ (Im $P_3$) と $-9\cot3\pi\nu_x = -9\cot\frac\theta2$ ($\theta = 3\mu = 6\pi\nu$ なので $3\pi\nu = \theta/2$) を使い、
   $$\Delta K_{2200} = -18\cdot2\cos\tfrac\theta2\,\mathrm{Im}(\cdots) - 9\cot\tfrac\theta2\cdot\bigl(-4\sin\tfrac\theta2\bigr)\mathrm{Im}(\cdots) = (-36 + 36)\cos\tfrac\theta2\,\mathrm{Im}\bigl(Xe^{i\theta/2}\bigr) = 0.$$
   調和 1 ($\zeta^2\bar\zeta$) は、$\theta\to\mu$、係数が $-6$ と $-3\cot\frac\mu2$ で、同じ計算が $(-12 + 12) = 0$ になる。8 極などの $h^{(4)}_{22000}$ は実数で位相に依らない。
6. したがって **$h_4$ の対の項の起点依存は、正規形の補正の $\cot$ の項が打ち消す**。数値で確認: ランダムな 5 要素 (位相・強さ・$c_{21}$ も任意、$\nu_x = 0.2347$) で起点を 1 つずつ動かして、(18.35) の $P$ と $\cot$ の項の和は 15 桁以上一致する ($-0.87037331765626$。要素 1 個の移動の $\Delta P_3$、$\Delta|h|^2$ の式との差も $10^{-15}$ 程度): [nl_03_gap_fills.out.txt](outputs/nl_03_gap_fills.out.txt) (c)。nl_02 (D) では、起点を最初の 6 極の後ろに移した式が小数 6 桁で一致する。

**1 個の薄い 6 極** ($P_3 = P_1 = 0$、$h_4^{(4)} = 0$): $|h_{30000}| = a/48$、$|h_{21000}| = a/16$ ($a = K_2L\beta^{3/2}$) を (18.35) に入れる。$9/48^2 = 1/256$、$3/16^2 = 3/256$:

$$
\frac{\partial\nu_x}{\partial J_x} = -\frac4\pi\cdot\frac{a^2}{256}\left[\cot3\pi\nu + 3\cot\pi\nu\right] = -\frac{(K_2L)^2\beta_x^3}{64\pi}\left[\cot3\pi\nu_x + 3\cot\pi\nu_x\right]
\tag{18.36}
$$

(sympy で (18.35) との差が 0。**追跡との比較** (nl_02 (A)): $a = 0.6\times5^{3/2}$ で $\nu_x = 0.2347$ のとき式 $-0.5722$、追跡 ($J$ の格子でフィット、NAFF) $-0.5722$ (相対差 $-7\times10^{-7}$)、$\nu_x = 0.41$ で $-0.4489$ (相対差 $1.7\times10^{-5}$)。)

**複数の要素 ($x$-$y$ を含む一般の場合)**: 4 変数 ($\zeta, \bar\zeta, v, \bar v$) の多項式の括弧 (18.9) を $y$ の項も入れて使い、対角の単項式 $\zeta^2\bar\zeta^2$、$\zeta\bar\zeta v\bar v$、$v^2\bar v^2$ の係数 $K_{2200}$、$K_{1111}$、$K_{0022}$ を (18.32) で求める。(18.24) から $\partial\nu_x/\partial J_x = 4K_{2200}/\pi$、$\partial\nu_x/\partial J_y = 2K_{1111}/\pi$、$\partial\nu_y/\partial J_y = 4K_{0022}/\pi$。nl_02 では、位相・ベータ関数・強さの違う 3 個の薄い 6 極の列 ($\nu = (0.2347, 0.3141)$。$3\nu_x$、$\nu_x\pm2\nu_y$ は整数から離れている) で:

| | $\partial\nu_x/\partial J_x$ | $\partial\nu_x/\partial J_y$ | $\partial\nu_y/\partial J_y$ |
|---|---|---|---|
| 式 (18.32) (多項式の括弧と正規形) | 5.6856 | −0.9613 | 0.6391 |
| 追跡 (NAFF + 2 次フィット、$J \le 3\times10^{-5}$) | 5.6853 | −0.9611 | 0.6392 |
| 相対差 | $-5\times10^{-5}$ | $-1.6\times10^{-4}$ | $5\times10^{-5}$ |

1 次元の明示の式 (18.35) の $\partial\nu_x/\partial J_x$ は、一般の式と 6 桁一致 (5.685598)。強さを ×0.5、×0.25 にすると式は 1.4214、0.3553 ($\propto$ 強さの 2 乗)、追跡との相対差は強さによらず同じ ($J$ の高次の寄与で、強さの高次ではない)。**6 極の ADTS は強さの 2 乗に比例する**: $g \propto F_3 \propto K_2L$ なので (18.32) の $[g, Rg] \propto (K_2L)^2$。

**docs との関係**: docs/24 §3.4・docs/25 §10.4 の「写像の対角部分 ($h_{22000}$ など) + 3 次の非共鳴の摂動 (正規形の補正、$\cot\pi x_k\,|h_{3,k}|^2$ の型)」は、(18.32) の第 1 項 ($h_4|_{\rm diag}$) と第 2 項 ($\cot$ の項) に対応する構造である (pysad の `rdt2nd.nf_correction` との係数の一致は本章では取り直していない)。

### 3.9 クロマティシティ: $\xi = h_{11001}/\pi$

$h_{11001}\,\zeta\bar\zeta\,\delta = h_{11001}(2J_x)\delta$ は $J_x$ に比例する項で、$K = 2h_{11001}J_x\delta$。(18.24) と同じ導出で $\nu_x = \nu_{x0} + \frac1{2\pi}\partial K/\partial J_x = \nu_{x0} + \frac{h_{11001}}{\pi}\delta$:

$$
\xi_x = \frac{h_{11001}}{\pi} = \frac1{4\pi}\sum K_2L\,\eta_x\beta_x,\qquad
\xi_y = \frac{h_{00111}}{\pi} = -\frac1{4\pi}\sum K_2L\,\eta_x\beta_y
\tag{18.37}
$$

(6 極の 1 次の寄与。四極の自然クロマティシティは別。)

**独立な導出**: 分散のある位置の 6 極は、閉軌道のずれ $\eta\delta$ のまわりで $-\frac{K_2L}{2}(x + \eta\delta)^2 \ni -K_2L\,\eta\delta\,x$ の線形な kick、つまり強さ $k = K_2L\,\eta\delta$ の薄い四極 ($\Delta p = -kx$) として働く。正規化座標で 1 周の行列は回転 $R(\mu)$ と kick $\begin{pmatrix}1 & 0\\ -k\beta & 1\end{pmatrix}$ の積で、トレースの半分 $\cos\mu' = \cos\mu - \frac{k\beta}{2}\sin\mu$。$\mu' = \mu + \Delta\mu$ の 1 次: $-\sin\mu\,\Delta\mu = -\frac{k\beta}2\sin\mu$ から $\Delta\mu = k\beta/2$、すなわち $\Delta\nu = k\beta/(4\pi) = K_2L\,\eta\beta\,\delta/(4\pi)$。$\xi_x$ の (18.37) と一致 (sympy (i))。

### 3.10 1 本の共鳴の近く

共鳴 $\mathbf m\cdot\boldsymbol\nu = N$ の近くでは、その調和は (18.30) で消せない (分母 $e^{i\mathbf m\cdot\boldsymbol\mu} - 1\to0$)。例: $3\nu_x = N$ ($\nu_x = N/3 + \delta$、$|\delta|\ll1$)。1 周の写像から、共鳴角 $3\phi$ だけを含む**1 自由度のハミルトニアン**を次の順で作る。

1. **共鳴の調和だけを残す**: §3.8 の変換 $e^{:g:}$ で、$g$ を共鳴でない調和に対してだけ (18.30) で作る。$\zeta^3$、$\bar\zeta^3$ の調和は分母が 0 に近いので $g$ に入れない。変換後の写像は $M' = e^{:F':}R$ ((18.29) の $e^{:d:}R$)、
   $$F' = -\bigl(h_{30000}\,\zeta^3 + \bar h_{30000}\,\bar\zeta^3\bigr) - c\,J^2$$
   (3 次: $F_3 = -h_3$ のうち共鳴の調和の項だけが残る。4 次: 対角項 $-K = -cJ^2$)。ここで $cJ^2$ は (18.32) の対角項から、共鳴の調和 $(3,0)$、$(0,3)$ を作る対 ($[g, Rg]$ の $18i$、$-18i$ の組、(18.33) の第 1 項) を**除いた**もの。$g_{30}$ が存在しないので、その対はない。$h_4$ の側の対 (18.34) と 8 極などの項は残る。$(2J)^2 = 4J^2$ から
   $$c = 4\left[h^{(4)}_{22000} - 18\,\mathrm{Im}\,P_3 - 6\,\mathrm{Im}\,P_1 - 3|h_{21000}|^2\cot\pi\nu\right].$$
   本節の例 (1 個の薄い 6 極、$P_3 = P_1 = 0$、$h^{(4)} = 0$) では $c = -12|h_{21000}|^2\cot\pi\nu$。
2. **回転座標に移す**: $R = e^{-\mu:J:}$ で $\mu = 2\pi N/3 + 2\pi\delta$ なので $R = R_{N/3}R_\delta$ ($R_{N/3} = e^{-(2\pi N/3):J:}$、$R_\delta = e^{-2\pi\delta:J:}$。どちらも $J$ の関数の流れで可換)。$R_{N/3}$ は $\zeta^p\bar\zeta^q$ に $e^{i(p-q)2\pi N/3}$ を掛けるので、$\zeta^3$、$\bar\zeta^3$ ($p - q = \pm3$ で $e^{\pm2\pi iN} = 1$) と $J$ の関数は**変えない**: $R_{N/3}F' = F'$。したがって $R_{N/3}\,:F':\,R_{N/3}^{-1} = :R_{N/3}F': = :F':$ で、$R_{N/3}$ は $e^{:F':}$ と可換:
   $$M' = R_{N/3}\;e^{:F':}R_\delta,\qquad M'^3 = R_{N/3}^3\bigl(e^{:F':}R_\delta\bigr)^3 = \bigl(e^{:F':}R_\delta\bigr)^3\qquad(R_{N/3}^3 = e^{-2\pi N:J:} = 1).$$
   つまり、1 周あたり $2\pi N/3$ ずつ回る座標 ($\phi\to\phi + 2\pi N/3$ を 1 周の回転として取り除いた系) で見ると、残りの回転は $2\pi\delta$ だけで、$\zeta^3$ の位相 $3\phi$ は 1 周に $6\pi\delta$ しか動かない。元の座標の周期 3 の不動点 (3 周で戻る 3 点) は、回転座標では $e^{:F':}R_\delta$ の不動点に対応する。
3. **1 周の写像を 1 つのハミルトニアンの流れで近似する**: BCH で $e^{:F':}e^{-2\pi\delta:J:} = \exp\{:F' - 2\pi\delta J + \tfrac12[F', -2\pi\delta J] + \cdots:\}$。共鳴の近傍では $\delta J$ と $F'\sim h(2J)^{3/2}$ が同じ大きさ ($\epsilon$ とおく) で、交換子の項は $\epsilon\,\delta$ (小さいパラメータ $\delta$ の分だけ小さい)。$O(\epsilon\delta)$ を落とすと、1 周の写像は時間 1 の流れ $e^{-:K:}$、$K = -F' + 2\pi\delta J$:
   $$K(J, \phi) = 2\pi\delta\,J + c\,J^2 + 2|h_{30000}|\,(2J)^{3/2}\cos(3\phi + \chi)
   \tag{18.38}$$
   (共鳴の項: $\zeta = \sqrt{2J}e^{i\phi}$ ((18.9) の次の段落) から $h\zeta^3 + \bar h\bar\zeta^3 = 2|h|(2J)^{3/2}\cos(3\phi + \chi)$、$\chi = \arg h_{30000}$。**「3 周を連続の時間とみなす」平均化の近似**は、この BCH の $O(\epsilon\delta)$ を落とすことに当たり、$\delta$ が小さいほど良い。)
   数値の確認 (nl_03 (d)): $h = 0.02$、$\chi = 0.4$、$c = -0.5$、$\delta = 10^{-3}$ ($\nu = 1/3 + \delta$) の写像 $e^{:F':}R$ (流れは数値積分) の周期 3 の不動点の $J = 9.77494\times10^{-4}$ が、(18.38) の不動点 ($3\phi + \chi = \pi$、$\partial K/\partial J = 0$) の $9.77457\times10^{-4}$ と相対 $4\times10^{-5}$ で一致する ($c$ を無視した (18.39) は $1.37\times10^{-3}$ で 40% 違う)。

**固定点**: $\partial K/\partial\phi = 0$ から $\sin(3\phi + \chi) = 0$。$\partial K/\partial J = 2\pi\delta + 2cJ + 6\sqrt2|h|\sqrt J\cos(3\phi + \chi) = 0$ ($(2J)^{3/2}$ の $J$ 微分は $3\sqrt2\sqrt J$)。$\cos = -1$ の枝で $c \to 0$ なら

$$
\sqrt{J_{\rm fp}} = \frac{\pi\delta}{3\sqrt2\,|h_{30000}|},\qquad J_{\rm fp} = \frac{\pi^2\delta^2}{18|h_{30000}|^2}
\tag{18.39}
$$

($c$ を含めると $\sqrt J = \frac{\sqrt2}{2c}\left[3|h| - \sqrt{9|h|^2 - 2\pi c\delta}\right]$ など 2 根。sympy。) **追跡との比較** (nl_02 (E): 1 個の薄い 6 極、$a = 0.5\times5^{3/2}$、$h_{30000} = a/48$): $\delta = 0.004$ で理論の $J_{\rm fp} = 5.973\times10^{-4}$ ($c$ を無視すると $6.468\times10^{-4}$)。1 周の写像の周期 3 の不動点 ($M^3z = z$) の $J$ は、3 点で $5.850\times10^{-4}$ と $6.087\times10^{-4}$、平均 $5.969\times10^{-4}$ (理論との相対差 $-7\times10^{-4}$)。$\delta = 0.008$ で理論 $2.234\times10^{-3}$、追跡の平均 $2.254\times10^{-3}$ (相対差 $9\times10^{-3}$)。**非共鳴の調和の項 $c$ を入れないと 8〜16% ずれる** ($\delta = 0.004$ で $6.468	imes10^{-4}/5.973	imes10^{-4}$ の 8.3%、$\delta = 0.008$ で $2.587	imes10^{-3}/2.234	imes10^{-3}$ の 15.8%) ので、(18.38) の $c J^2$ は必要である。

**島の幅 (一般の共鳴、2 自由度)**: 作用 $(J_x, J_y)$ から共鳴角 $\theta = m_x\phi_x + m_y\phi_y$ とその共役な作用 $I$ に変数を取り替える ($J_a = J_a^{(0)} + m_aI$ と書ける: 母関数 $\theta I + \phi_yI_2$)。非共鳴の 2 次の項 $\frac12\sum\partial^2K/\partial J_a\partial J_b\,\Delta J_a\Delta J_b = \pi\sum A_{ab}\Delta J_a\Delta J_b$ ($A_{ab} = \partial\nu_a/\partial J_b$、$K = \pi A_{ab}J_aJ_b$ から $\nu = \frac1{2\pi}\partial K/\partial J$ で整合) に $\Delta J_a = m_aI$ を入れると $\pi\,G\,I^2$、$G \equiv \sum m_aA_{ab}m_b$。1 周のハミルトニアンは

$$
K(I, \theta) = \pi G\,I^2 + V\cos\theta,\qquad V = 2\Bigl|\textstyle\sum_q h_q\,(2J_x)^{(j+k)/2}(2J_y)^{(l+m)/2}\Bigr|\ \ (q\text{ はその共鳴の調和に属する単項式})
\tag{18.40}
$$

($I$ の 1 次の項 (離調) は $I$ の原点のずらしで消す。) $G > 0$、$V > 0$ なら不安定固定点は $(I, \theta) = (0, 0)$ ($K = V$)、安定は $\theta = \pi$。分離線 $K = V$: $\pi GI^2 = V(1 - \cos\theta)$、最大は $\theta = \pi$ で $I_{\max} = \sqrt{2V/(\pi G)}$。チューン平面の離調 $\Delta_0 = \mathbf m\cdot\boldsymbol\nu(\mathbf J) - N$ は $\Delta_0 = G\,\Delta I$ なので、島の半幅

$$
W = G\,I_{\max} = \sqrt{\frac{2\,G\,V}{\pi}}
\tag{18.41}
$$

(docs/26 §3.5 の式と同じ。sympy (j)。島に捕らえられる条件は $\pi\Delta_0^2/|G| < V(1 - {\rm sgn}(G)\cos\theta)$。**この式を追跡 (周波数マップの D の大きい領域) と比べた結果は docs/26 §3.5・§9.1**: 予想できるのは小振幅・孤立した低次共鳴だけ。)

### 3.11 強さと振幅の次数の勘定

| 量 | 強さへの依存 | 振幅への依存 | 理由 |
|---|---|---|---|
| $h_3$ (6 極の RDT) | $K_2L$ の 1 乗 | $J^{3/2}$ | (18.22)、$\zeta^{(j+k)}\ldots$ の次数 3 |
| $h_4$ の 8 極・運動学・端部 | $K_3L$ の 1 乗 (運動学は強さに依らない) | $J^2$ | (18.25)〜(18.27) |
| $h_4$ の 6 極の対の項と正規形の補正 | $(K_2L)^2$ | $J^2$ | $[H_3, H_3]$、$[g, Rg]$ ($g \propto F_3 \propto K_2L$) |
| 6 極 × 偏向磁石の 3 次の項 | $K_2L$ の 1 乗 (偏向磁石の項は 6 極に依らない) | $J^2$ | 同じ括弧の交差項 |
| ADTS ($\partial\nu/\partial J$) | $J^2$ の項の係数 ($\partial/\partial J$) | $J^0$ | $K \propto J^2$ |

**偶数・奇数の分離**: 6 極の強さを $\lambda$ 倍すると ADTS は $\lambda^2$ の項 (6 極の 2 次) と $\lambda$ の項 (交差項) と $\lambda^0$ の項 (運動学など) の和なので、$f(\lambda) = A\lambda^2 + B\lambda + C$。$f(+1) + f(-1) - 2f(0) = 2A$、$(f(+1) - f(-1))/2 = B$。これが docs/24 §5.1 の「偶数部が 6 極 2 次」の比較の仕組み。なお (18.35) の導出は 2 次の摂動までで、強さの 4 乗以上 ($g$ の高次の括弧) は含まない。

### 3.12 周波数マップ・局所運動量アクセプタンスの式

**NAFF** (Numerical Analysis of Fundamental Frequencies): 信号 $z_n = A\,e^{2\pi i\nu n}$ ($n = 0\ldots N-1$、$z = \hat x - i\hat p$ なら $\nu > 0$) に窓 $w_n$ をかけ、試験周波数 $\nu'$ との内積 $S(\nu') = \frac AN\sum_nw_n\,e^{2\pi i(\nu - \nu')n}$ の絶対値が最大になる $\nu'$ を探す。$\Delta = \nu - \nu'$ として、矩形窓 ($w_n = 1$) の核は

$$
D(\Delta) = \frac1N\sum_{n=0}^{N-1}e^{2\pi i\Delta n} = e^{i\pi\Delta(N-1)}\frac{\sin\pi N\Delta}{N\sin\pi\Delta},\qquad |D| \to \frac{1}{\pi N|\Delta|}\ (N\Delta \gg 1)
\tag{18.42}
$$

(等比級数の和。) 1 本の線だけなら $|S|$ の最大は $\Delta = 0$ で、窓に依らず $\nu' = \nu$ (窓は正で対称)。2 本目 ($\nu_1$、振幅 $\epsilon$) があると、裾 $S\sim\epsilon\,D(\nu_1 - \nu')$ が 1 本目のピークを引き寄せる。**余弦 2 乗窓** (Hann 窓、$w_n = 1 - \cos(2\pi n/N)$) は $e^{\pm2\pi in/N}$ の和なので、核は $D(\Delta)$ の 3 つの和

$$
S_{\rm Hann}(\Delta) = D(\Delta) - \tfrac12D\left(\Delta - \tfrac1N\right) - \tfrac12D\left(\Delta + \tfrac1N\right)
\tag{18.43}
$$

($\Delta = 3.37/N$ で直接の和との差 $6\times10^{-17}$。nl_02 (F)。) (18.43) の括弧を、(18.42) の形に戻して展開する。

1. **ずらした核**: $D(\Delta\mp1/N) = e^{i\pi(\Delta\mp1/N)(N-1)}\dfrac{\sin(\pi N\Delta\mp\pi)}{N\sin\pi(\Delta\mp1/N)}$。位相因子は $e^{i\pi\Delta(N-1)}\,e^{\mp i\pi(N-1)/N} = e^{i\pi\Delta(N-1)}\cdot(-1)\,e^{\pm i\pi/N}$ ($e^{\mp i\pi} = -1$、$-(N-1)/N = -1 + 1/N$)、分子は $\sin(\pi N\Delta\mp\pi) = -\sin\pi N\Delta$。2 つの $-1$ が打ち消し合うので
   $$D\!\left(\Delta\mp\tfrac1N\right) = e^{i\pi\Delta(N-1)}\,\frac{\sin\pi N\Delta}{N}\cdot\frac{e^{\pm i\pi/N}}{\sin\pi(\Delta\mp1/N)}.$$
   (振動因子 $e^{i\pi\Delta(N-1)}\sin\pi N\Delta$ は 3 項に共通で、符号は変わらない。位相は $e^{\pm i\pi/N}$ が残る。)
2. **括弧**: (18.43) に代入すると
   $$S_{\rm Hann}(\Delta) = \frac{e^{i\pi\Delta(N-1)}\sin\pi N\Delta}{N}\left[\frac1{\sin\pi\Delta} - \frac12\frac{e^{i\pi/N}}{\sin\pi(\Delta - 1/N)} - \frac12\frac{e^{-i\pi/N}}{\sin\pi(\Delta + 1/N)}\right].$$
3. **括弧の展開**: $\varphi(\eta) \equiv e^{-i\pi\eta}/\sin\pi(\Delta + \eta)$ と置くと、括弧は $\varphi(0) - \tfrac12[\varphi(1/N) + \varphi(-1/N)]$ (第 2 項は $\eta = -1/N$ で $e^{-i\pi\eta} = e^{i\pi/N}$、第 3 項は $\eta = +1/N$)。$\eta$ の 2 次まで展開すると、1 次の項は打ち消し合って
   $$\varphi(0) - \tfrac12\left[\varphi\!\left(\tfrac1N\right) + \varphi\!\left(-\tfrac1N\right)\right] = -\frac1{2N^2}\varphi''(0) + O(N^{-4}).$$
   $\Delta\ll1$ ($\sin\pi\Delta\simeq\pi\Delta$) で、$\varphi\simeq e^{-i\pi\eta}/[\pi(\Delta + \eta)]$ の 2 階微分は $\varphi''(0) = \dfrac{2}{\pi\Delta^3} + \dfrac{2i}{\Delta^2} - \dfrac{\pi}{\Delta}$ (積の微分: $(1/\pi(\Delta+\eta))'' = 2/\pi\Delta^3$、交差項 $2\cdot(-i\pi)\cdot(-1/\pi\Delta^2) = 2i/\Delta^2$、$(e^{-i\pi\eta})'' = -\pi^2$ から $-\pi/\Delta$)。$N\Delta\gg1$ では第 1 項が主要部で、括弧 $\simeq -\dfrac{1}{\pi N^2\Delta^3}$。
4. したがって $|S_{\rm Hann}|\simeq\dfrac{|\sin\pi N\Delta|}{N}\cdot\dfrac1{\pi N^2\Delta^3} \le \dfrac{1}{\pi(N\Delta)^3}$ で、**裾は $1/\Delta^3$** に減る (矩形窓の $1/(\pi N\Delta)$ に対し、差分で $1/\Delta$ と $1/\Delta^2$ の項が消える)。数値 ([nl_03_gap_fills.out.txt](outputs/nl_03_gap_fills.out.txt) (e)): $N\Delta = 5.5$ ($|\sin| = 1$) で $|S| = 1.979\times10^{-3}$、式 $1/(\pi\cdot5.5^3) = 1.913\times10^{-3}$、$N\Delta = 20.5$ で $3.7026\times10^{-5}$ と $3.6948\times10^{-5}$。$N = 1024$ では $N\Delta = 80.5$ でも $6.1013\times10^{-7}$ と $6.1019\times10^{-7}$。$\Delta$ が小さくない場合 ($N = 256$、$N\Delta = 80.5$ で $\Delta = 0.31$) は $\sin\pi\Delta\simeq\pi\Delta$ が崩れて式から外れる。括弧の式と (18.43) の直接の和は 5 桁以上一致する


nl_02 (F) の確認 ($|S|/|S(0)| = 2.4\times10^{-2}$ ($\Delta = 2.5/N$)、$2.0\times10^{-3}$ (5.5)、$3.7\times10^{-5}$ (20.5)、$5.6\times10^{-7}$ (80.5)。矩形窓の $|D|$ は 0.127、0.058、0.0157、0.0047)。2 本の線 (振幅比 $10^{-3}$、距離 0.0245) の下での $\nu_0$ の誤差は、$N = 256$、1024、4096 で矩形窓 $-8\times10^{-8}$、$-1\times10^{-8}$、$-2\times10^{-10}$、Hann 窓 $4\times10^{-9}$、$3\times10^{-11}$、$5\times10^{-14}$。(精度の $N$ 依存の指数は、本章では測っていない。docs/24 §6 は数値で確認している。)

**周波数マップ**: 追跡の前半・後半の各面の信号に NAFF をかけてチューン $\nu_1$、$\nu_2$ を求め、$D = \log_{10}\sqrt{(\nu_{x2} - \nu_{x1})^2 + (\nu_{y2} - \nu_{y1})^2}$。規則的な運動 ($\nu$ が一定) では $\Delta\nu$ は NAFF の誤差で、$D$ は小さい (FODO の線形な運動で $D\approx-11$、docs/26 §3.1)。共鳴・カオスで $\nu$ が周回とともに動くと $D$ が上がる。ADTS (§3.8) は $\nu$ が $J$ に依ることを意味し、RDT (§3.5) で $J$ が周回とともに変動すると、ADTS を通じてチューンが周回ごとに動き、前後半のずれ ($D$) になる (定性的な説明。定量的な関係は docs/26 §9.1)。「周波数マップが結果の地図、RDT と ADTS が原因の係数」という位置づけである。

**局所運動量アクセプタンス (LMA) の RF バケットの限界**: 縦方向を、1 周あたり $\Delta\psi = \kappa\,\delta$、$\Delta\delta = -a\left[\sin(\varphi_s + \psi) - \sin\varphi_s\right]$ の写像 ($\psi$ は同期位相からの位相のずれ、$\kappa = 2\pi h\alpha_p$ ($= kC\alpha_p$、14 章の $2\pi h\alpha_p$ と同じ大きさ。符号の約束は下の段落)、$a = V/(pc/e)$ ($V$ は 1 周の電圧、$\sin\varphi_s = U_0/V$)) とし、1 周を連続の時間とみなすと、ハミルトニアン

$$
H(\psi, \delta) = \tfrac12\kappa\,\delta^2 + a\left[-\cos(\varphi_s + \psi) - \psi\sin\varphi_s\right]
\tag{18.44}
$$

($-\partial H/\partial\psi = -a[\sin(\varphi_s + \psi) - \sin\varphi_s]$ が $\Delta\delta$ を与える。$\kappa a\cos\varphi_s > 0$ で $\psi = 0$ が安定。$\kappa > 0$ を仮定。) 力が 0 の点は $\sin(\varphi_s + \psi) = \sin\varphi_s$ から $\psi = 0$ (安定) と $\psi_u = \pi - 2\varphi_s$ (不安定)。分離線は不安定固定点を通る等高線で、$\psi = 0$ での $\delta_{\rm rf}$ は

$$
\tfrac12\kappa\,\delta_{\rm rf}^2 = a\left[2\cos\varphi_s - (\pi - 2\varphi_s)\sin\varphi_s\right]
\ \Rightarrow\ \delta_{\rm rf}^2 = \frac{2a}{\kappa}\left[2\cos\varphi_s - (\pi - 2\varphi_s)\sin\varphi_s\right]
\tag{18.45}
$$

($H(\psi_u, 0) - H(0, 0) = a[-\cos(\pi - \varphi_s) - (\pi - 2\varphi_s)\sin\varphi_s + \cos\varphi_s] = a[2\cos\varphi_s - (\pi - 2\varphi_s)\sin\varphi_s]$。sympy で一致。) pysad の `lmabucket.rf_bucket_analytic` の式 ($A = V/E_0$ ($E_0 = p_0c$)、$\kappa = kC\alpha$ ($kC = 2\pi h$) と置いた形。ソースの docstring) と同じ形で、14 章 §5 の `bh` とも同じ式 ($a = V/p$。符号の約束の対応は下の段落)。減衰・$\alpha_p$ の 2 次を含まない理想化した式で、docs/26 の LMA の値は追跡で求める。docs/26 の `zlimit` $= (\pi + 2|\varphi_s|)/k$ は、力が 0 になる 2 つの不安定点 $\psi = \pi - 2\varphi_s$、$-\pi - 2\varphi_s$ の外側の値で、分離線が $\delta = 0$ と交わる反対側の位相 ($\varphi_s > 0$ ではこれより内側。(18.44) の線形項 $-a\psi\sin\varphi_s$ のため) とは別である。判定の端としての妥当性は本章では確認していない。

**14 章 §5 との符号の約束の対応**: (18.44) は、SAD・pysad の約束 ($\psi = \varphi_s + kz$、$z$ は先行が正、`lmabucket` の docstring) の写像 $\Delta\psi = -\kappa\delta$、$\Delta\delta = +a[\sin(\varphi_s + \psi) - \sin\varphi_s]$ ($H_{\rm SAD} = -\tfrac12\kappa\delta^2 + a[\cos(\varphi_s + \psi) + \psi\sin\varphi_s]$) の**両式の右辺の符号を同時に反転した** ($H = -H_{\rm SAD}$、時間を逆向きにたどる) ものである。軌道の曲線 (等高線)・固定点・安定性・分離線は同じで、(18.45) の $\delta_{\rm rf}$ も同じ ($H_{\rm SAD}$ では、分離線の条件 $-\tfrac12\kappa\delta_{\rm rf}^2 + a\cos\varphi_s = a[-\cos\varphi_s + (\pi - 2\varphi_s)\sin\varphi_s]$ から (18.45) が出る)。14 章 §5 は位相を時間の向き ($\varphi = \omega_{\rm RF}t$、遅れて着くと正) で書き、$d\delta/dn = +(V/E)[\sin(\varphi_s + \varphi) - \sin\varphi_s]$、$d\varphi/dn = 2\pi h\alpha_p\delta$ とする。これは SAD の約束から $\varphi = -\psi$、$\varphi_{s,14} = \pi - \varphi_s$ と置いた形である ($\sin(\pi - \varphi_s + \varphi) = \sin(\varphi_s - \varphi) = \sin(\varphi_s + \psi)$、$\sin\varphi_{s,14} = \sin\varphi_s = U_0/V$ は同じ)。したがって 14 章の $\cos\varphi_{s,14}$ は本章の $\cos\varphi_s$ と符号が逆で、安定性の条件は 14 章の $\alpha_pV\cos\varphi_{s,14} < 0$ と本章の $\kappa a\cos\varphi_s > 0$ で、同じ内容の別の書き方になる (14 章 §5 の式は別途見直されるので、ここに書いた対応は上の形の 14 章を前提とする)。

## 4. 追跡の写像から取り出す

**生成関数とテイラー係数**: 1 周の写像 (または 1 要素の写像) を $z\circ\Phi = e^{:f_3:}e^{:f_4:}z$ ($f_3$、$f_4$ は 3 次式、4 次式。$[f_3, f_4]$ は 5 次なので 4 次までは $e^{:f_3 + f_4:}$ と同じ) と書くと、$z$ のテイラー展開は

$$
z_i\circ\Phi = z_i + [f_3, z_i] + \underbrace{[f_4, z_i] + \tfrac12[f_3, [f_3, z_i]]}_{3\text{ 次}} + \cdots
\tag{18.46}
$$

($e^{:f_3:}e^{:f_4:}z = e^{:f_3:}(z + [f_4, z] + \cdots) = z + [f_3, z] + \frac12[f_3, [f_3, z]] + [f_4, z] + (4\text{ 次以上})$。) **2 次の項**: $u_i \equiv [f_3, z_i]$ は 1 自由度では $(u_q, u_p) = (-\partial f_3/\partial p,\ \partial f_3/\partial q)$ ($[f, q] = -\partial_pf$、$[f, p] = \partial_qf$)。したがって写像の 2 次のテイラー係数 $T_2$ は $(\partial_pf_3, \partial_qf_3)$ を直接与え、**$f_3$ は一意に決まる** (3 次の同次式は $f_3 = \frac13z\cdot\nabla f_3$ で復元できる: オイラーの関係。sympy で 3 次式について確認)。例: 6 極の kick $\Delta p = -\frac{K_2L}{2}q^2$ から $\partial f_3/\partial q = -\frac{K_2L}{2}q^2$、$f_3 = -\frac{K_2L}{6}q^3 = -H_3L$。$F = -h_3$ の符号どおり。

**3 次の項**: 一般の 3 次式 $f$ について $[f, [f, z_i]] = [f, u_i] = (u\cdot\nabla)u_i$ ($[f, g] = \sum(\partial_qf\,\partial_pg - \partial_pf\,\partial_qg) = u_q\partial_qg + u_p\partial_pg$ で、$u = (-\partial_pf, \partial_qf)$ を使った。sympy で一般の 3 次式について差 0)。3 次のテイラー係数を $\tau_3$ とすると (18.46) から

$$
[f_4, z_i] = \tau_{3,i} - \tfrac12(u\cdot\nabla)u_i,\qquad f_4 = \tfrac14\,z\cdot\nabla f_4
\tag{18.47}
$$

(4 次の同次式のオイラーの関係)。これが pysad の `rdttrack.f4_from_map3` の原理 (docs/25 §10.5 の「$v = \tau_3 - \frac12(u\cdot\nabla)u$ から $\nabla h_4 = -Sv$」と同じ。符号は $f_4 = -h_4$)。$\tau_3$ は方向を変えた奇関数の差分から取り、誤差は 7 次の項で決まる (docs/25 §10.5。本章では導出しない)。

**偏向磁石の 3 次の生成関数** (`bendnl.py`): 偏向磁石を「入口の端・$n$ スライス (既定 16)・出口の端」に分けると、各スライスは恒等に近いので、生成関数 = **スライスのハミルトニアン × 長さ** (刻みの 1 次まで)。取り出した生成関数は、解析式 $\tfrac h2ds\,x(p_x^2 + p_y^2)$、端 $\pm\tfrac h2p_xy^2$、複合機能の本体 $\tfrac{K_1h}{6}ds\,x^3$ と一致する (docs/24 §5.5.1)。解析式 ((18.18)) で積分するには、6 極との交差項を閉じた式にするのが難しいため、1 要素の写像のテイラー展開から係数を取る。

解析式は、**入れた効果の範囲と摂動論の次数で原理的に作る**。追跡に合わせるための補正項・当てはめは足さない。例外はこの偏向磁石の 3 次の生成関数。

| 道具 | 何をするか |
|---|---|
| `resdrive.py` (`resonance_terms`、`adts`、`chromaticity_*`) | 解析式: $h_{jklmp}$、ADTS、クロマティシティ。6 極 + 偏向磁石の 3 次の項 (`bend_cubic=True`)。4 次 (`order2=True`) |
| `rdt2nd.py` | 4 次の項 $h_{jklm}$ (6 極の 2 次) |
| `rdttrack.py` | 追跡の 1 周の写像のテイラー展開から $h$ を取り出す (解析との照合) |
| `nform.py` | 数値の正規形 (追跡の写像から) |
| `adtsdirect.py` | 実際の強さの追跡 + NAFF + 多項式フィット (ADTS の基準値) |
| `bendnl.py` | 偏向磁石の 3 次のリー生成関数 (追跡の写像から取り出す) |
| `freqmap.py`、`lma.py` | 周波数マップ、局所運動量アクセプタンス |

K64 (gs) には RDT・ADTS・周波数マップ・LMA の機能がない (すべて pysad のみ)。追跡自体は K64 と一致する (docs/26 §5)。

## 5. この章の位置づけ (ほかの章との関係)

- 要素の章の「行列と追跡の関係」: 1 周の写像の 1 次は CALC の行列 (設計軌道上で追跡のヤコビアンと一致)。この章はその 2 次以上。
- 分割の次数 (03・04・07 章) は、1 周の写像の 4 次の項に**数値積分器の端数**を残す: 四極の本体が $\int p^4/8$ の台形則 (K64) であること (`kinematic_h4(quad_rule="k64")`) が、解析との 6 極に依らない 4 次の定数の差の原因と推定されている (docs/25、未解決)。
- 放射・スピン (13〜15 章) は、1 周の写像がシンプレクティックでなくなる (放射) が、ここでの解析は放射なし。

## 6. 解析と追跡の比較の要約

**この章の確認** (§3 の導出を確かめたもの。pysad は使わない):

| 比較 | 結果 |
|---|---|
| 括弧・展開・合成の公式 ((18.1)〜(18.22)、(18.31)、(18.46)) | sympy で差 0 ([nl_01](outputs/nl_01_lie_rdt_adts_derivations.out.txt)) |
| 1 個の薄い 6 極の ADTS (18.36) vs 追跡 | $\nu_x = 0.2347$: 相対差 $-7\times10^{-7}$、$\nu_x = 0.41$: $1.7\times10^{-5}$ ([nl_02](outputs/nl_02_adts_formula_vs_tracking.out.txt) (A)) |
| 3 個の 6 極の列 (18.32) vs 追跡 | $-5\times10^{-5}$、$-1.6\times10^{-4}$、$5\times10^{-5}$ (xx、xy、yy) (B) |
| 同、8 極 1 個を足す | $2\times10^{-4}$、$4\times10^{-4}$、$-7\times10^{-6}$ (C) |
| 起点の取り方 | 式が 6 桁で一致 (D) |
| 3$\nu_x$ の近くの周期 3 の不動点の $J$ (18.39) | 相対差 $-7\times10^{-4}$ ($\delta = 0.004$)、$9\times10^{-3}$ ($\delta = 0.008$) (E) |

**pysad の解析と追跡の比較** ([docs/24](../docs/24_adts_rdt_primer.md) §7 から転記。数字は取り直していない):

| 比較 | 一致 |
|---|---|
| 6 極 2 次 (強さ ×(+1, −1, 0) の追跡の偶数部) | FODO で $4\times10^{-5}$・$7\times10^{-6}$・$2\times10^{-4}$ 以内、ASTRID2 で $1.5\times10^{-3}$・$1.1\times10^{-3}$・$1.7\times10^{-3}$、大きなリング A で 6 極 $2\times10^{-5}$・$3\times10^{-5}$・$10^{-3}$、8 極 $4\times10^{-6}$・$8\times10^{-6}$・$2\times10^{-5}$ |
| 合計 (全部入り) vs 追跡 + NAFF + フィット (6 極 ×1) | FODO $-0.04\%$、ASTRID2 $-0.24\%$・$-0.08\%$・$-0.19\%$、大きなリング A $-0.014\%$・$-0.017\%$・$-0.26\%$ |
| 6 極の強さを 0.5〜2 倍 | $0.003$〜$0.7\%$ |
| 厚い 6 極 vs 薄レンズ列 | $1/n^2$ で収束 |
| RDT: 符号・規格化・位相 (FODO) | pysad の $h$ は $H$ の係数そのもの ((18.22))。docs/24 §4.1 は、文献の薄い 6 極の式との比 $-1$ を報告している。符号の差は $h$ の定義 ($\zeta^\pm$ の向き、$F = \pm h$) の違いによると考えられるが、文献との符号の照合は未了。1 周の $\Delta J$ と 0.5% 以内、3$\nu_x$ の近くの固定点の振幅と 0.15〜1.3% |
| RDT: 追跡の 1 周の写像の 2 次のテイラー展開 | 6 極の分は 6 リングで $10^{-6}$ 以内、偏向磁石つきの全体は $10^{-8}$ (FODO・ASTRID2)〜$4\times10^{-7}$ (大きなリング A)・$5\times10^{-5}$ (MPW) |
| ADTS: 4 次の項 (`adts(fourth='r6')`) と写像の対角 + 正規形の補正 | $10^{-15}$ ([docs/25](../docs/25_resonance_terms_adts.md) §10.4) |

大きなリング A の $yy$ で、解析に何を足すかによる差の推移: 6 極 2 次 + 8 極 1 次だけ $-34\%$ → 運動学的な項 + K1 の端のフリンジを足して $+3.6\%$ → 偏向磁石の 3 次どうしと 6 極との交差項を足して $-0.26\%$ (docs/24 §7.3)。

## 7. 使った条件

nl_02 の模型: 正規化座標の線形回転と薄い 6 極・8 極 (kick の式は §3.1)、チューン $(0.2347, 0.3141)$、$\beta_{x,y} = (6, 2), (2, 7), (9, 3)$ m、$K_2L = 0.8, -1.1, 0.9$ [1/m²]、1024〜4096 周、$J \le 3\times10^{-5}$ m rad (比較の J の格子は $\{0, 1, 2, 3\}\times10^{-5}$)。エネルギー・周長・磁場は模型に入らない (正規化座標)。

pysad の比較 (§6 の転記): エネルギー・周長: FODO 3 GeV、105.6 m (24 セル、$\nu_x = 6.696$、$\nu_y = 4.906$)。ASTRID2 0.58 GeV、45.7 m。大きなリング A 3 GeV、603 m (HMBA 型、$\nu = 47.58, 17.62$)。追跡の条件の共通部分: k64 モード (gs と同じ追跡)、RF 空洞なし (`rfsw=False`)、放射なし、結合・誤差なし、$\delta = 0$。磁場の上限は各ラティスの設計値 (docs/24 の冒頭)。**gs の台本で ADTS を求めるときは `NORFSW;` が必要** (gs は FFS で RFSW が ON になり、RF つきのラティスで自然クロマティシティ × $\delta(J)$ の偽の ADTS が入る。docs/25 §5)。

## 8. 未確認・限界・導出の範囲

**本章で導出したもの**: (18.1)〜(18.47) のうち、多極・運動学・フリンジの $H$、$\zeta^\pm$ と括弧、順序つき指数の展開 ($h_3$、$h_4$)、薄い 6 極の $h_{jklmp}$、1 次の ADTS (8 極・運動学・端部)、2 次の ADTS (正規形、1 次元の閉じた式、一般の式の構造)、クロマティシティ、1 本の共鳴、NAFF の窓、RF バケットの高さ、テイラー係数からの生成関数。

**導出していないもの** (docs/24・25・26 にある。何が足りないか):
- 厚い要素の積分 $\int w^j\bar w^k\cdots ds$ と要素内の対 $\iint_{s_1 > s_2}$ の閉形式 (ドリフトでは $w$ が $s$ の 1 次式、K1≠0 では $\cos/\sin$ か $\cosh/\sinh$ の和: 積分の実行が必要。docs/25 §3・§10.2)。
- 偏向磁石の 3 次の生成関数の解析式・$E_1$、$E_2\ne0$ の端の閉形式 (docs/25 §8・§11.3。05・06 章の写像の展開が必要)。
- $\delta\ne0$ の ADTS、$\delta$ 微分 (docs/25 §11.4)、非対角の調和の正規形の補正 (docs/25 §11.5。(18.30) の $g$ は全調和に使えるが、4 次の非対角の項を消す次の次数の変換と、その共鳴の幅への効果は取り直していない)。
- $J^2$ の項 ($\partial^2\nu/\partial J^2$): 5・6 次の生成関数と、3 次の $g_3$ の反復 (3 回の括弧) が要る。docs/25 §13。
- 結合・ROTATE した磁石の光学 (docs/25 §14)。
- pysad の `resdrive`・`rdt2nd` の数値と本章の式・nl_02 の直接の突き合わせ (nl_02 は pysad を使わない)。

**限界**:
- 解析式は、誤差なし・結合なしを前提 (ROTATE した磁石は結合した光学の模型に自動で切り替える。docs/25 §14)。ASTRID2 の追跡 1 例で、結合・閉軌道 mm 級の誤差で 10〜40% 動く。
- 解析に入っていない項: 6 極に依らない 4 次の定数 (上、原因は積分器の端数と推定。未解決)、端の 4 次の項、6 極 × 8 極・8 極の 2 次、10 極以上。
- $\delta\ne0$ の ADTS の $\delta$ 微分は精密版 (`fourth='r6'`・`optics='tracked'`) で $x$-$x$・$x$-$y$ が 0.3% 以内 ($y$-$y$ は 28% 外れ、原因未特定)、既定は $10^{-3}$ 台の違いが微分で数十% に増える。$\partial^2\nu/\partial\delta^2$ の解析式はない。
- RDT の追跡との照合: 1 周の $\Delta J$・共鳴線の近くの検証は FODO のみ (3$\nu_x$、$y = 0$ の 1 本)。$\delta$ の項 ($p > 0$) は $h_{11001}$・$h_{00111}$ とクロマティシティの比較まで。スペクトルの線から $h$ を求める方法は未実装。
- v2 モードの追跡での偏向磁石の 3 次の項の検証は、k64 の追跡とだけ。
- 文献との照合 (符号・規格化) は未了 (docs/24 §4.1 が薄い 6 極の式との比を報告しているが、文献の定義との照合はしていない)。本章の (18.19)〜(18.21) は定義 ($H_3$ の単項式の係数をそのまま読む) から導いたもので、文献の式には依らない。
- 共鳴の近くの (18.38) は、1 周の写像を時間 1 の流れで置き換える近似 (§3.10 手順 3 の BCH で $O(\epsilon\delta)$ を落とす)。$\delta$ が小さいほど良く、$\delta = 0.008$ で 0.9% (nl_02 (E)、$J^{3/2}$ の項の高次は無視)。
- (18.45) は減衰・放射・$\alpha_p$ の 2 次を含まない理想化した式。docs/26 の LMA の値は追跡で求め、この式は上限の見積りに使う。
- ヘルプには RDT・ADTS の式は一切ない。

## 付録 A. ヘルプの記述・K64 との違い

| 箇所 | 内容 |
|---|---|
| ヘルプ | RDT・ADTS・正規形・周波数マップ・LMA の記述はない (K64 の機能にもない)。この章の式は標準の理論 (リー代数の正規形) と pysad の実装による。 |
| K64 の仕様 | 四極の本体の $\int p^4/8$ が台形則 (積分器の端数) で、ADTS の解析との 6 極に依らない 4 次の定数の差の原因と推定される (docs/25)。誤りではなく積分器の仕様。 |

## 付録 B. K64 の実装の癖・疑い (別枠)

この章で新しく見つけた疑いはない。関連する K64 の事項: 追跡と CALC の分割の次数の違い (薄い多極 4 次と 2 次、四極の刻み。04・07 章)、CALC の `tthine` に非線形フリンジがない (B291、07 章)。これらは 1 周の写像の 2 次以上と、CALC の行列の違いに関わる。

## 付録 C. この章で使う名前

(この章は Fortran の変数を使わない。)

| 名前 | 意味 | ファイル |
|---|---|---|
| $h_{jklmp}$ | 共鳴駆動項 | `resdrive.resonance_terms` |
| $J$、$\phi$、$\zeta^\pm$ | 作用・角度・複素変数 | docs/24 §2〜3、本章 §3.2 |
| `adts().total`・`.sext`・`.octu` | ADTS の合計・6 極 2 次・8 極 1 次 | `resdrive.adts` |
| `bend_cubic`、`fourth`、`optics` | 偏向磁石の 3 次の項、4 次の項の精密版、光学の取り方 | `resdrive` |
| `direct_adts` | 追跡 + NAFF + フィットの ADTS | `adtsdirect` |
| `nf_correction` | 正規形の補正 ((18.33) の $\cot$ の項に対応する構造) | `rdt2nd` |
| `f4_from_map3` | 3 次のテイラー係数から $f_4$ ((18.47)) | `rdttrack` |
| `adts_formula`・`pbracket` | 本章の確認用: (18.32) の一般の式・多項式の括弧 | `scripts/nl_02_adts_formula_vs_tracking.py` |
| `pb`・`kfun`・`flow` | 本章の確認用 (補足): (18.9) の括弧、(18.35) の非 $h_4$ 部、共鳴近傍の流れ | `scripts/nl_03_gap_fills.py` |
