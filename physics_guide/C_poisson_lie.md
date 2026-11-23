# 付録 C. ポアソン括弧とリー演算子 $\exp(:f:)$

**この付録で分かること**: SAD のヘルプには、要素の転送が $\exp(:H:)$ や $\exp(:V:)$ の積で書いてある。この記法の意味を、ポアソン括弧から順に、最も簡単な例 (ドリフト、集束四極、6 極のキック) で途中式つきで説明する。

**前提**: ハミルトンの正準方程式 ([01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §5、§8)。リー代数は知らなくてよい。

**読む順**: §1 (結論) → §2〜§4 (記法の意味) → §5 (例。ここが中心) → §6〜§7 (合成の順序、線形写像) → §8 (SAD での対応)。

数式の番号 $(C.n)$ は、この付録の中と他の章からの参照に使う。確認の sympy スクリプトは `physics_guide/scripts/ham_07_lie_operator.py`、出力は [ham_07_lie_operator.out.txt](outputs/ham_07_lie_operator.out.txt)。

## 1. 結論

- **ポアソン括弧** $[f, g]$ は、相空間の関数 2 つから関数を作る演算で、正準方程式は $q' = [q, H]$、$p' = [p, H]$ と書ける。
- **リー演算子** $:f:$ は「$f$ とのポアソン括弧を取る」演算子 ($:f:\,g = [f, g]$)。$\exp(:f:) = \sum_n :f:^n/n!$ は、$f$ を生成関数とする相空間の変換 (リー変換) を表す。
- ハミルトニアン $H$ (位置 $s$ によらない) で、長さ $s$ だけ進む転送は $\exp(-s:H:)$、つまり $f = -sH$ の $\exp(:f:)$。ドリフトは $x \to x + s p$ に、集束四極は $\cos$・$\sin$ の行列になり、6 極のキックは、級数が 2 項で止まって $p \to p - k_2x^2$ になる (§5)。
- 2 次の $f$ の $\exp(:f:)$ は線形写像 (行列 $\exp(F)$)、3 次以上の $f$ は非線形な写像。いずれもシンプレクティック (§7)。
- この付録の規約: $[f, g] = \sum_i\left(\dfrac{\partial f}{\partial q_i}\dfrac{\partial g}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i}\right)$。**文献によって括弧の符号が逆**のものがある。SAD のヘルプの $\exp(:H:)$ ($H$ は長さをかけたハミルトニアン) は、この付録と括弧の符号が逆の流儀で、ハミルトニアン (長さをかけたもの) の時間 1 の流れ (§8.2。ヘルプの式どうしの整合で確定)。

## 2. ポアソン括弧と正準方程式

自由度 $n$ の系の座標 $q_i$、運動量 $p_i$ ($i = 1\ldots n$) の関数 $f$、$g$ に対し、ポアソン括弧を

$$
[f, g] \equiv \sum_{i=1}^{n}\left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i} - \frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right)
\tag{C.1}
$$

と定義する。基本的な値は

$$
[q_i, p_j] = \delta_{ij},\qquad [q_i, q_j] = [p_i, p_j] = 0,\qquad [f, g] = -[g, f],\qquad [f, f] = 0
\tag{C.2}
$$

(例: $[x, p] = \partial_xx\,\partial_pp - \partial_px\,\partial_xp = 1$。)

**正準方程式との関係**。独立変数を $s$ とし、$q' = \partial H/\partial p$、$p' = -\partial H/\partial q$ (§5、§8 の導出)。相空間の関数 $g(q, p)$ が、粒子が動くにつれてどう変わるかを、鎖則で書くと

$$
\frac{dg}{ds} = \sum_i\left(\frac{\partial g}{\partial q_i}q_i' + \frac{\partial g}{\partial p_i}p_i'\right)
= \sum_i\left(\frac{\partial g}{\partial q_i}\frac{\partial H}{\partial p_i} - \frac{\partial g}{\partial p_i}\frac{\partial H}{\partial q_i}\right)
= [g, H]
\tag{C.3}
$$

$g = q_i$ とすれば $q_i' = [q_i, H] = \partial H/\partial p_i$、$g = p_i$ とすれば $p_i' = [p_i, H] = -\partial H/\partial q_i$ (つまり、正準方程式 2 本が、1 つの式 $g' = [g, H]$ にまとまる)。例: $H = p^2/2 + V(x)$ で、$[x, H] = p$、$[p, H] = -V'(x)$。確認は sympy (a)。

## 3. リー演算子 $:f:$

関数 $f$ を 1 つ決めて、関数 $g$ に $[f, g]$ を対応させる演算子を $:f:$ と書く:

$$
:f:\,g \equiv [f, g]
\tag{C.4}
$$

繰り返すと、

$$
:f:^2 g = [f, [f, g]],\qquad :f:^3 g = [f, [f, [f, g]]],\qquad \ldots
\tag{C.5}
$$

(左から順に $f$ とのポアソン括弧を取る。$:f:^0 g = g$。) (C.3) は $g' = [g, H] = -[H, g] = -:H:\,g$ と読める。

## 4. $\exp(:f:)$ と転送

指数関数の級数と同じように

$$
\exp(:f:) \equiv \sum_{n=0}^{\infty}\frac{:f:^n}{n!} = 1 + :f: + \frac{:f:^2}{2} + \frac{:f:^3}{6} + \cdots
\tag{C.6}
$$

と定義する ($\exp(:f:)\,g = g + [f, g] + \frac12[f, [f, g]] + \cdots$)。

**転送との関係**。ハミルトニアン $H$ が $s$ によらないとする。(C.3) を繰り返し使う ($g$ の $s$ 微分が $-:H:g$、その $s$ 微分が $(-:H:)^2g$、...)。出発点 $z_0$ から長さ $s$ 進んだ点での $g$ の値を、$s$ のテイラー級数で書くと

$$
g\big(z(s)\big) = \sum_{n=0}^{\infty}\frac{s^n}{n!}\left(\frac{d^ng}{ds^n}\right)_{z_0}
= \sum_{n=0}^{\infty}\frac{s^n}{n!}\left(-:H:\right)^n g\Big|_{z_0}
= \bigl(\exp(-s:H:)\,g\bigr)(z_0)
\tag{C.7}
$$

特に $g$ を座標 $z_i$ ($x, p_x$ など) にとると、**長さ $s$ の要素の転送は $\exp(-s:H:)$**、つまり $f = -sH$ として $\exp(:f:)$。(右辺は「初期値の関数」。出発点の座標から出口の座標への写像を与える。) 級数の収束は、ここでは形式的に扱う (§9)。

## 5. 最も簡単な例

以下、1 自由度 $(x, p)$。$f$ が $s$ を含むのは、$f = -sH$ だから。確認は、スクリプトの (b)(c)(d) の節。

### 5.1 例 1: ドリフト $f = -\dfrac{s\,p^2}{2}$ (スクリプトの (b))

$H = p^2/2$ (小さい角度の近似のドリフト)。(C.1) から $:f:\,x = [f, x] = \partial_xf\,\partial_px - \partial_pf\,\partial_xx = 0 - (-sp)\cdot1 = s\,p$。次に、$:f:^2x = [f, sp] = \partial_xf\,\partial_p(sp) - \partial_pf\,\partial_x(sp) = 0 - (-sp)\cdot0 = 0$ (以降も 0)。$:f:\,p = [f, p] = \partial_xf\cdot1 - \partial_pf\cdot0 = 0$。したがって、級数は

$$
\exp(:f:)\,x = x + s\,p,\qquad \exp(:f:)\,p = p
\tag{C.8}
$$

ドリフトの行列 $\begin{pmatrix}1 & s\\ 0 & 1\end{pmatrix}$ そのもの。

### 5.2 例 2: 集束四極 (線形) $f = -\dfrac{s}{2}\left(k\,x^2 + p^2\right)$ (スクリプトの (c))

$H = (kx^2 + p^2)/2$、$k > 0$ (集束)。$:f:\,x = [f, x] = -\partial_pf = s\,p$、$:f:\,p = [f, p] = \partial_xf = -s\,k\,x$。繰り返すと

$$
\begin{aligned}
:f:\,x &= s\,p, & :f:^2x &= s\,(:f:p) = -k s^2x, & :f:^3x &= -ks^2(:f:x) = -k s^3p, & :f:^4x &= k^2s^4x,\ \ldots\\
:f:\,p &= -ks\,x, & :f:^2p &= -ks\,(:f:x) = -ks^2p, & :f:^3p &= k^2s^3x, &&
\end{aligned}
\tag{C.9}
$$

(線形なので、$:f:$ は $(x, p)$ の中でしか動かない。) (C.6) に入れ、$w \equiv \sqrt k$ とすると

$$
\exp(:f:)\,x = x\left(1 - \frac{w^2s^2}{2} + \frac{w^4s^4}{24} - \cdots\right) + \frac{p}{w}\left(ws - \frac{w^3s^3}{6} + \cdots\right)
= x\cos ws + \frac{p}{w}\sin ws
\tag{C.10}
$$

$$
\exp(:f:)\,p = -w\,x\sin ws + p\cos ws
\tag{C.11}
$$

(級数がそれぞれ $\cos$、$\sin$ のテイラー級数になる。) これは、四極の章 ([04_quad.md](04_quad.md) §4.1) の式 (4.11) の水平の行列で $p = 1$ とおいたものに一致する。$k < 0$ では $\cos \to \cosh$、$\sin \to \sinh$。$k \to 0$ の極限は、(C.8) のドリフト。スクリプトの (c) で、$s$ の 12 次未満の級数が $\cos$、$\sin$ の級数と一致することを確認した。

### 5.3 例 3: 6 極のキック $f = -\dfrac{k_2}{3}x^3$ (スクリプトの (d))

6 極の薄いキックを、位置 $x$ だけの関数のハミルトニアン $H = (k_2/3)x^3$ に長さ 1 (薄いレンズ) の流れとして書いた例。$k_2$ はここでは単なる係数で、SAD の K2 ($K2 = B''L/B\rho$) との関係 (係数の $1/2$ など) は確認していない。

$:f:\,x = [f, x] = -\partial_pf = 0$、$:f:\,p = [f, p] = \partial_xf = -k_2x^2$。次に $:f:^2p = [f, -k_2x^2] = \partial_xf\,\partial_p(-k_2x^2) - \partial_pf\,\partial_x(-k_2x^2) = 0$ ($f$ は $p$ を含まず、$-k_2x^2$ も $p$ を含まないので、どちらの項も 0)。したがって級数は 2 項で止まる:

$$
\exp(:f:)\,x = x,\qquad \exp(:f:)\,p = p - k_2\,x^2
\tag{C.12}
$$

ハミルトンの方程式から直接作る ($\Delta p = -\partial H/\partial x$) と、$\Delta p = -k_2x^2$ で一致。**$f$ が座標だけ (または運動量だけ) の多項式なら、$\exp(:f:)$ は有限の項で止まる**。逆に、(b) のように座標と運動量の両方を含む $f$ では、級数が無限に続く ($\cos$、$\sin$)。

(参考: 薄いレンズの変換 (C.12) のポアソン括弧は $[x, p - k_2x^2] = 1$ で、シンプレクティック。スクリプトの (d) で確認した。)

## 6. 合成と作用の順序

### 6.1 順序の規約

2 つの変換を続けて行うとき、$\exp(:f:)\exp(:g:)$ が表す変換の順序に注意が必要。リー演算子は**関数**に作用する。座標の関数 $z_i$ に作用させると

$$
\exp(:f:)\,\exp(:g:)\,z_i = \exp(:f:)\big[\exp(:g:)\,z_i\big]
\tag{C.13}
$$

右の演算子が、先に関数 $z_i$ に作用する。ところが、座標の関数への作用は、粒子の座標を代入する操作と順序が逆になる。具体例 (ドリフト $f = -sp^2/2$ と、薄いキック $g = -kx^2/2$ ($p \to p - kx$)):

- $\exp(:g:)x = x$、$\exp(:g:)p = p - kx$。
- $\exp(:f:)\exp(:g:)x = \exp(:f:)\,x = x + sp$。
- $\exp(:f:)\exp(:g:)p = \exp(:f:)(p - kx) = p - k(x + sp)$。

これは、粒子に**先にドリフト ($f$)、次にキック ($g$)** を作用させた結果 ($x \to x + sp$、$p \to p - k(x + sp)$) と同じである。逆の積 $\exp(:g:)\exp(:f:)$ は $x \to x + sp - ksx$、$p \to p - kx$ で、違う (sympy (e))。

**本ガイドの規約**: 式 $\exp(:f_1:)\exp(:f_2:)\cdots\exp(:f_n:)$ は、**左から順に、粒子が通る順** ($f_1$ が先、$f_n$ が最後) の変換を表す。ビームラインの要素の並びと式の並びが同じ向きになる。

注意: 文献には、積を右から作用する演算子の順 (上と逆)、または、写像の合成 $M_n\circ\cdots\circ M_1$ の順に書くものもある。他の文献や、SAD のヘルプの式を読むときは、まず順序の規約を確かめる。

### 6.2 BCH 公式

$\exp(:f:)\exp(:g:)$ は、一般に $\exp(:f + g:)$ に**等しくない** ($:f:$ と $:g:$ が可換でないため)。リー演算子はヤコビ恒等式から $[:f:, :g:] = :[f, g]:$ を満たす。ベイカー・キャンベル・ハウスドルフ (BCH) の公式を使うと

$$
\exp(:f:)\exp(:g:) = \exp\left(:\,f + g + \frac12[f, g] + \frac1{12}\Big([f, [f, g]] + [g, [g, f]]\Big) + \cdots:\right)
\tag{C.14}
$$

($f$ と $g$ が可換 $[f, g] = 0$ なら、$\exp(:f + g:)$ に等しい。) 上のドリフトとキックの例 ($f = -\varepsilon p^2/2$、$g = -\varepsilon x^2/2$) で、$[f, g]/2$ までの $\exp$ が $\varepsilon^3$ の差を残し、3 次の項 $\frac{1}{12}([f,[f,g]] + [g,[g,f]])$ を足すとその差が消えることを、sympy で確認した ((f))。$[:f:, :g:] = :[f, g]:$ も、いくつかの多項式の $f$、$g$、$h$ で確認した。

**注意**: BCH の無限級数が実際に収束するか、どの次数まで足すか、は別の問題 (§9)。SAD のヘルプは、要素の変換を $\exp(:\cdot:)$ の積 (例えば DECA) で書き、次数をそろえるための補正項 $\exp(:V_4:)$ を足している (§8)。

### 6.3 $H$ が $s$ に依存する場合

(C.7) は $H$ が $s$ によらないとして導いた。$H$ が $s$ に依存すると、$\exp(-\int H\,ds\,:\,:)$ とは書けない。区間を細かく分けて各区間の積 $\exp(-\Delta s_1:H(s_1):)\exp(-\Delta s_2:H(s_2):)\cdots$ (§6.1 の順) を取り、分割を細かくした極限が正しい転送 (時間順序つき指数)。1 つの $\exp(:\Omega:)$ にまとめた形が**マグヌス展開**で、$\Omega = -\int H\,ds + \frac12\int\!\!\int[\ldots] + \cdots$ のように、ポアソン括弧の項が続く。この付録では、詳細は扱わない。区間ごとに $H$ が一定な (要素ごとの) 場合は、区間ごとの $\exp(-\Delta s:H:)$ の積 (§6.1 の順) になる。

## 7. 線形写像と 2 次の $f$、シンプレクティック性

$f$ が $z = (q_1, p_1, \ldots)$ の**2 次の斉次多項式**なら、$:f:\,z_k$ は $z$ の 1 次式なので、$:f:$ は行列で表せる。$H = \frac12\mathbf z^TA\mathbf z$ ($A$ は対称行列) で、(C.3) は $z_k' = \sum_l J_{kl}\,\partial H/\partial z_l$ ($J$ は $\begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$ の直和) だから、$\mathbf z' = JA\mathbf z$。したがって

$$
\exp(-s:H:)\,\mathbf z = \exp(sJA)\,\mathbf z = M\,\mathbf z,\qquad M = \exp(F),\;\;F = sJA
\tag{C.15}
$$

5.2 の例 ($A = \mathrm{diag}(k, 1)$) では $F = \begin{pmatrix}0 & s\\ -ks & 0\end{pmatrix}$ で、行列指数は (C.10)(C.11) の行列に一致する (sympy (g))。**3 次以上の $f$ のリー変換は非線形写像**になる ($:f:$ が次数を保たず、(C.12) のように、高次の項を作る)。

**シンプレクティック性**。$M = \exp(F)$、$F = sJA$ とすると、$A$ が対称、$J^TJ = 1$、$J^2 = -1$ から $F^TJ + JF = s(AJ^TJ + J\,JA) = s(A - A) = 0$。したがって $\frac{d}{ds}\left(M^TJM\right) = M^T(F^TJ + JF)M = 0$ で、$s = 0$ では $M = 1$ なので $M^TJM = J$ ($M$ はシンプレクティック行列)。(5.2 の例で、$\det M = 1$、$M^TJM = J$ を確認。) 一般の (非線形の) $\exp(:f:)$ も、座標のポアソン括弧 $[Z_i, Z_j] = [z_i, z_j]$ を保つ (正準変換)。(C.12) で、$[x, p - k_2x^2] = 1$ を確認した。これが、リー変換の積で作る写像が、丸め誤差までシンプレクティックに保たれる理由。

## 8. SAD のヘルプ、pysad・Fortran との対応

### 8.1 ヘルプに出てくる式

K64 版同梱のヘルプ原文 (`Documents/SADHelp.html`) で、$\exp(:\cdot:)$ が出てくるのは次の式 (原文で確認。DRIFT・QUAD・BEND・SEXT・OCT・DECA・DODECA。CAVI・MULT・SOL の節の式は 06・08・09・10 章):

| ヘルプの節 | 式 (ヘルプの記述の骨子) |
|---|---|
| DRIFT の transformation | $\exp(:H:)$、$H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ ($E/v_0$ の符号は誤植とみられる。02 章付録 A) |
| QUAD の本体 | スライスごとに $\exp(:H_{2n}:)$ ($H_{2n} = \frac1n\{(-p + \frac{p_x^2 + p_y^2}{2p} + \frac E{v_0})L + \frac{K1}{2}(x^2 - y^2)\}$)。スライスの間に $\exp(:\Delta H:)$ ($\Delta H = \frac1n(p - \sqrt{p^2 - p_x^2 - p_y^2} - \frac{p_x^2 + p_y^2}{2p})L$) |
| BEND の線形フリンジ | $\exp(:V:)$、平文版 $V = -\tfrac1p\bigl(\tfrac{F1^2}{24\rho_b}p_x + \tfrac{F1}{12\rho_b^2}y^2 - \tfrac{y^4}{6F1\rho_b^2}\bigr)$ (TeX 版は $y^2$・$y^4$ の項が違う誤植。05 章付録 A)。QUAD のフリンジは母関数 (04 章) |
| SEXT・OCT・DECA・DODECA の transformation | $\exp(:F_{\rm in}:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:V_n:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:F_{\rm out}:)$ ($L$ はドリフト、$H_n$ は薄い多極のキック $K_n/(n+1)!\,\mathrm{Re}(x - iy)^{n+1}$、$a = \frac12 - \frac1{\sqrt{12}}$、$b = \frac12 - a$、$V_n$ は 3 次を揃える補正)。原文の並びは前後対称でない (誤記。07 章付録 A) |

ヘルプの $H$ は、位置 $s$ の長さ $L$ (またはスライスの長さ) をかけた量 ($H = H_{\rm ham}\cdot L$) になっていて、そのうえで $\exp(:H:)$ と書かれている。

### 8.2 括弧の符号の食い違い (重要)

この付録の規約では、転送は $\exp(-s:H_{\rm ham}:)$、つまり $\exp(:f:)$ の $f$ は $-sH_{\rm ham}$ ((C.7))。ヘルプは、$f = +L\,H_{\rm ham}$ を $\exp(:f:)$ に入れている (例: DRIFT の $H$ は、ハミルトニアン $-\sqrt{p^2 - p_x^2 - p_y^2} + \ldots$ に $L$ をかけたもの)。ヘルプの $:f:$ は、この付録の $:f:$ と符号が逆、$:f:_{\rm SAD}\,g = -[f, g] = [g, f]$、すなわち **$\exp(:H:)$ は「$H$ を時間 1 のハミルトニアンとする流れ」** ($x \to x + \partial H/\partial p_x$、$p_x \to p_x - \partial H/\partial x$)。

**確認** (ヘルプ原文には括弧の定義が無い。式どうしの整合で確定した。[drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt) (1)): (i) DRIFT の $\sqrt{\ }$ の項から $x' = p_x/p_z$ が正しい向きで出る。(ii) QUAD の $H_{2n}$ の $K1(x^2 - y^2)/2$ から $\Delta p_x = -K1\,x$ (K1 > 0 で水平収束。ヘルプの K1 の説明と一致)。(iii) SEXT の $H_2 = K2/3!\,\mathrm{Re}(x - iy)^3$ から $\Delta p_x = -(K2/2)(x^2 - y^2)$。(iv) BEND の線形フリンジ $V$ (平文版) から、ヘルプ自身の $x_2$、$p_{y2}$ の式が出る。(v) 04 章・07 章・03 章の数値 ($V_n$ の符号) とも整合 ([thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt))。標準の括弧で $\exp(:H:)$ と読むと、これらはすべて逆向きになる。唯一の例外が DRIFT の $E/v_0$ の符号 (QUAD・dynamics・BEND の body と逆、誤植とみられる)。

式を $\exp(:H_{\rm help}:)$ から本付録の記法に直すときは、$H_{\rm help} \to -H_{\rm help}$ とする。

順序の規約 (§6.1) については、ヘルプの DECA などの式が、入口のフリンジ $F_{\rm in}$ を最左、出口の $F_{\rm out}$ を最右に書いている (入口から出口の順)。これは §6.1 の規約 (左が先) と同じ向きに読める。ヘルプは順序の規約を明記していない (原文で確認)。ここでは、DECA などの式の並び (入口から出口の順) と整合する読みとして使う。

### 8.3 pysad・Fortran との対応

- **pysad `bendnl.py`** (偏向磁石の 3 次のリー生成関数を、追跡の写像から取り出す)。コードの説明に「薄い kick の生成関数 $f$ は、kick 写像 $z \to z + \tau(z)$ から $S\nabla f = \tau$ ($S = \begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$) で決まる」とある。$S\nabla f = (\partial_pf, -\partial_xf)$ なので、これはハミルトニアン $f$ の流れ ($\Delta x = \partial f/\partial p$、$\Delta p = -\partial f/\partial x$)、すなわちこの付録の記法では $\exp(:-f:)$ (ヘルプの $\exp(:f:)$ と同じ符号の向き)。`docs/25_resonance_terms_adts.md` の $M = R\,e^{:f_3:}e^{:f_4:}\cdots$ との符号の対応は、この付録では確かめていない。
- **K64 の Fortran** (四極の本体 `tsolqu`、薄い多極 など)。実装は、$\exp(:H_2:)$ などを級数ではなく**閉じた式** (四極の本体なら (4.9)(4.10)、キックなら $p \to p - \ldots$) で書いている。ヘルプの $\exp(:H:)$ を、コードの更新式にどう対応させたかは、各要素の章で書く。薄い多極のコード (OCT・DECA の $a$、$b$、$V$ の積の実装) との対応はこの付録では**未確認**。

## 9. 未確認・限界

- 級数 (C.6) (C.7) は**形式的なもの**として扱った (収束の議論はしていない)。ドリフト (C.8) や 6 極 (C.12) のように級数が止まる場合と、(C.10) のようにテイラー級数が既知の関数になる場合だけを例にした。
- BCH (C.14) は 3 次の項までを確認した。それ以上の項と、実際に使う分割の精度の議論は、[03_symplectic_integrators.md] (予定) に回す。
- ヘルプの $\exp(:\cdot:)$ の向きは、ヘルプに括弧の定義が無いので、式どうしの整合 (DRIFT・QUAD・SEXT・BEND の 4 つ) からの確定。定義の明記は原文に無い (§8.2)。順序の規約も明記が無い。
- マグヌス展開 ($H$ が $s$ に依存する場合) は、名前の紹介だけで、式を導いていない。
- 多自由度 (結合、$x$-$y$) の例は書いていない。2 次の $f$ の一般的な対応 (C.15) は、1 自由度の例だけ確認した。
- この付録の 6 極の係数 $k_2$ と SAD の K2 の関係は確認していない。

## 参考

- SAD のヘルプ (K64 版同梱の `Documents/SADHelp.html`): DRIFT、QUAD、BEND、SEXT・OCT・DECA・DODECA などの transformation の節。
- [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md): 正準方程式の導出 (§5、§8)。
- [04_quad.md](04_quad.md): 四極の線形写像、$H_2$ と $\Delta H$ の分割 (ヘルプの $\exp(:H_{2n}:)$、$\exp(:\Delta H:)$ に当たる)。
- リー演算子の一般論の文献は、この付録では参照していない (内容は、ポアソン括弧と指数関数の定義から導いたもの)。
