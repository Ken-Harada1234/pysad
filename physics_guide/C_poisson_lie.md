# 付録 C. ポアソン括弧とリー演算子 $\exp(:f:)$

**この付録で分かること**: SAD のヘルプには、要素の転送が $\exp(:H:)$ や $\exp(:V:)$ の積で書いてある。この記法の意味を、ポアソン括弧から順に、最も簡単な例 (ドリフト、集束四極、6 極のキック) で途中式つきで説明する。

**前提**: ハミルトンの正準方程式 ([01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) §5、§8)。リー代数は知らなくてよい。

**読む順**: §1 (結論) → §2〜§4 (記法の意味) → §5 (例。ここが中心) → §6〜§7 (合成の順序、線形写像) → §8 (SAD での対応)。

数式の番号 $(C.n)$ は、この付録の中と他の章からの参照に使う。確認の sympy スクリプトは `physics_guide/scripts/ham_07_lie_operator.py` (出力は [ham_07_lie_operator.out.txt](outputs/ham_07_lie_operator.out.txt)) と、途中式の確認を足した `physics_guide/scripts/ham_08_bch_jacobi.py` (出力は [ham_08_bch_jacobi.out.txt](outputs/ham_08_bch_jacobi.out.txt))、BCH の 3 次の語の表・ライプニッツ則・(C.15) の帰納法・§8.2 の偏微分を足した `physics_guide/scripts/ham_10_audit_steps.py` (出力は [ham_10_audit_steps.out.txt](outputs/ham_10_audit_steps.out.txt))。

## 要約と目次

**この付録で分かること**

- ヘルプの式 $\exp(:H:)$ の記法の意味を、ポアソン括弧 (C.1) から順に、ドリフト・集束四極・6 極のキックの 3 つの例で途中式つきで導く。
- 結果の主な式: 正準方程式の 1 本化 $dg/ds = [g, H]$ (C.3)、長さ $s$ の要素の転送 $\exp(-s:H:)$ (C.7)、ドリフト (C.8)・集束四極 (C.10)(C.11)・6 極のキック (C.12) の閉じた形、BCH の公式 (C.14)、2 次の $f$ が作る行列 $M = \exp(sJA)$ (C.15)、非線形な $\exp(:f:)$ も括弧を保つこと (C.16)(C.17)。
- 対応: ヘルプの DRIFT・QUAD・BEND の線形フリンジ・SEXT〜DODECA の `exp(:...:)` の式の読み方 (§8)。Fortran・pysad のコードそのものの導出ではない (K64 は閉じた式で書いてあり、薄い多極の $a$、$b$、$V$ の積との対応は未確認。pysad の `bendnl.py` の符号の向きだけ §8.3 に書いた)。
- 確かめ方: sympy (`ham_07`・`ham_08`・`ham_10`)。ヘルプの $\exp(:H:)$ の向きは、ヘルプの式どうしの整合から読んだもので、gs では確認していない。
- K64 の誤り: なし。ヘルプの誤り: この付録では新たな確定・疑いなし。§8.2 で、$\exp(:H:)$ を「$H$ を時間 1 のハミルトニアンとする流れ」と読むと DRIFT・QUAD・SEXT・BEND の式が整合し、唯一の例外が DRIFT の $E/v_0$ の符号 (02 章付録 A の確定した誤植と同じ。gs での確認なし)。B 番号なし。

**導出の流れ**

1. [§1 結論](#1-結論): 括弧・演算子・転送・例・括弧の符号の規約を先に述べる。
2. [§2 ポアソン括弧と正準方程式](#2-ポアソン括弧と正準方程式): 括弧の定義 (C.1)・基本の値 (C.2) から、正準方程式が $g' = [g, H]$ (C.3) になることを示す。
3. [§3 リー演算子 $:f:$](#3-リー演算子-f): $:f:\,g = [f, g]$ (C.4)、繰り返し (C.5)。
4. [§4 $\exp(:f:)$ と転送](#4-zzmathzz85end-と転送): 定義 (C.6)。$d^ng/ds^n = (-:H:)^ng$ の帰納法とテイラー展開で (C.7)。
5. [§5 最も簡単な例](#5-最も簡単な例)
    - [§5.1 ドリフト](#51-例-1-ドリフト-zzmathzz121end-スクリプトの-b): 級数が止まって (C.8)。
    - [§5.2 集束四極](#52-例-2-集束四極-線形-zzmathzz128end-スクリプトの-c): 級数 (C.9) が $\cos$・$\sin$ になって (C.10)(C.11)。
    - [§5.3 6 極のキック](#53-例-3-6-極のキック-zzmathzz167end-スクリプトの-d): 級数が 2 項で止まって (C.12)。
6. [§6 合成と作用の順序](#6-合成と作用の順序)
    - [§6.1 順序の規約](#61-順序の規約): (C.13)。本ガイドは左から粒子が通る順。
    - [§6.2 BCH 公式](#62-bch-公式): 3 次までの導出で (C.14)。
    - [§6.3 $H$ が $s$ に依存する場合](#63-zzmathzz336end-が-zzmathzz337end-に依存する場合): 時間順序つき指数とマグヌス展開の紹介 (導出なし)。
7. [§7 線形写像と 2 次の $f$、シンプレクティック性](#7-線形写像と-2-次の-zzmathzz348endシンプレクティック性): (C.15) を導き、$M$ がシンプレクティックであることを示す。
    - [§7.1 非線形な $\exp(:f:)$ も括弧を保つ](#71-非線形な-zzmathzz391end-も括弧を保つ-正準変換): ライプニッツ則と二項展開で (C.16) → (C.17)。
8. [§8 SAD のヘルプ、pysad・Fortran との対応](#8-sad-のヘルプpysadfortran-との対応)
    - [§8.1 ヘルプに出てくる式](#81-ヘルプに出てくる式): $\exp(:\cdot:)$ が出る式の一覧。
    - [§8.2 括弧の符号の食い違い](#82-括弧の符号の食い違い-重要): ヘルプの $:f:$ はこの付録と符号が逆。4 つの式で向きを確かめる。
    - [§8.3 pysad・Fortran との対応](#83-pysadfortran-との対応)
9. [§9 未確認・限界](#9-未確認限界)、[参考](#参考)。

式の連鎖: (C.1) → (C.3) → (C.6) → (C.7) → (C.8)・(C.10)(C.11)・(C.12)。合成は (C.13) → (C.14)、線形・非線形の写像の性質は (C.15) → (C.16) → (C.17)。

**読み方**: 前提は [01a](01a_hamiltonian_derivation.md) の §5・§8 (正準方程式)。§5 の例が中心。ヘルプの式を読むだけなら §1 と §8.2 の結論 (向きの規約、$H_{\rm help} \to -H_{\rm help}$) で足りる。BCH の導出 (§6.2) と非線形の括弧の保存 (§7.1) は、[03](03_symplectic_integrators.md) などで使うまで飛ばしてよい。文献によって括弧の符号が逆なので、他の文献と式を比べるときは §1 の規約を確かめる。

**未確認・限界の所在**: §9 (級数の収束、5 次以降の BCH、マグヌス展開、多自由度の例、ヘルプの向きが gs で未確認であること)。

## 1. 結論

- **ポアソン括弧** $[f, g]$ は、相空間の関数 2 つから関数を作る演算で、正準方程式は $q' = [q, H]$、$p' = [p, H]$ と書ける。
- **リー演算子** $:f:$ は「$f$ とのポアソン括弧を取る」演算子 ($:f:\,g = [f, g]$)。$\exp(:f:) = \sum_n :f:^n/n!$ は、$f$ を生成関数とする相空間の変換 (リー変換) を表す。
- ハミルトニアン $H$ (位置 $s$ によらない) で、長さ $s$ だけ進む転送は $\exp(-s:H:)$、つまり $f = -sH$ の $\exp(:f:)$。ドリフトは $x \to x + s p$ に、集束四極は $\cos$・$\sin$ の行列になり、6 極のキックは、級数が 2 項で止まって $p \to p - k_2x^2$ になる (§5)。
- 2 次の $f$ の $\exp(:f:)$ は線形写像 (行列 $\exp(F)$)、3 次以上の $f$ は非線形な写像。いずれもシンプレクティック (§7)。
- この付録の規約: $[f, g] = \sum_i\left(\dfrac{\partial f}{\partial q_i}\dfrac{\partial g}{\partial p_i} - \dfrac{\partial f}{\partial p_i}\dfrac{\partial g}{\partial q_i}\right)$。**文献によって括弧の符号が逆**のものがある。SAD のヘルプの $\exp(:H:)$ ($H$ は長さをかけたハミルトニアン) は、この付録と括弧の符号が逆の流儀で、ハミルトニアン (長さをかけたもの) の時間 1 の流れ (§8.2。ヘルプの式どうしの整合から読んだもので、SAD (gs) の実行での確認はしていない)。

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

**転送との関係**。ハミルトニアン $H$ が $s$ によらないとする。(C.3) を繰り返し使う。(C.3) は任意の関数 $h$ について $dh/ds = [h, H] = -:H:h$ ($H$ が $s$ を含まないので、右辺の $:H:$ は $s$ によらない演算子)。$g$ の 1 階微分は $dg/ds = -:H:g$。これも相空間の関数なので、$h = -:H:g$ に (C.3) を使うと $d^2g/ds^2 = -:H:(-:H:g) = (-:H:)^2g$。一般に $d^{n-1}g/ds^{n-1} = (-:H:)^{n-1}g$ なら、$h = (-:H:)^{n-1}g$ に使って $d^ng/ds^n = (-:H:)^ng$ (帰納法。例: $H = p^2/2$、$g = x^2$ で $n = 0\ldots3$ を厳密解のテイラー展開と比べた。[ham_08_bch_jacobi.out.txt](outputs/ham_08_bch_jacobi.out.txt) (c))。出発点 $z_0$ から長さ $s$ 進んだ点での $g$ の値を、$s$ のテイラー級数で書くと

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

(線形なので、$:f:$ は $(x, p)$ の中でしか動かない。) 規則は、$:f:^{2m}x = (-ks^2)^m x$、$:f:^{2m+1}x = (-ks^2)^m\,s\,p$、$:f:^{2m}p = (-ks^2)^m p$、$:f:^{2m+1}p = (-ks^2)^m(-ks)\,x$ ($m = 0, 1, 2, \ldots$。$n = 0\ldots6$ で sympy 確認、[ham_08_bch_jacobi.out.txt](outputs/ham_08_bch_jacobi.out.txt) (d))。(C.6) に入れ、$w \equiv \sqrt k$ とすると、$x$ の係数は偶数次だけ ($n = 2m$): $\sum_m\dfrac{(-1)^m(ws)^{2m}}{(2m)!} = \cos ws$ ($(-ks^2)^m = (-1)^m(ws)^{2m}$)。$p$ の係数は奇数次だけ ($n = 2m + 1$): $\sum_m\dfrac{(-1)^m s\,(ws)^{2m}}{(2m+1)!} = \dfrac1w\sum_m\dfrac{(-1)^m(ws)^{2m+1}}{(2m+1)!} = \dfrac{\sin ws}{w}$。つまり

$$
\exp(:f:)\,x = x\left(1 - \frac{w^2s^2}{2} + \frac{w^4s^4}{24} - \cdots\right) + \frac{p}{w}\left(ws - \frac{w^3s^3}{6} + \cdots\right)
= x\cos ws + \frac{p}{w}\sin ws
\tag{C.10}
$$

$$
\exp(:f:)\,p = -w\,x\sin ws + p\cos ws
\tag{C.11}
$$

($p$ の式は、$x$ の係数が $-w\sin ws$ ($:f:^{2m+1}p$ の $x$ の項 $-ks(-ks^2)^m\,x$ の和)、$p$ の係数が $\cos ws$。) これは、四極の章 ([04_quad.md](04_quad.md) §4.1) の式 (4.11) の水平の行列で $p = 1$ とおいたものに一致する。$k < 0$ では $\cos \to \cosh$、$\sin \to \sinh$。$k \to 0$ の極限は、(C.8) のドリフト。スクリプトの (c) で、$s$ の 12 次未満の級数が $\cos$、$\sin$ の級数と一致することを確認した。

### 5.3 例 3: 6 極のキック $f = -\dfrac{k_2}{3}x^3$ (スクリプトの (d))

6 極の薄いキックを、位置 $x$ だけの関数のハミルトニアン $H = (k_2/3)x^3$ に長さ 1 (薄いレンズ) の流れとして書いた例。$k_2$ はここでは単なる係数。SAD の K2 ($K2 = B''L/B\rho$) とは、§8.2 (iii) の式との比較で $k_2 = K2/2$ ($y = 0$) の関係になる (§9)。

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

($f$ と $g$ が可換 $[f, g] = 0$ なら、$\exp(:f + g:)$ に等しい。)

**[:f:, :g:] = :[f, g]: の導出**: 任意の関数 $h$ に作用させる。$[:f:, :g:]h = :f:(:g:h) - :g:(:f:h) = [f,[g,h]] - [g,[f,h]]$。ヤコビ恒等式 $[f,[g,h]] + [g,[h,f]] + [h,[f,g]] = 0$ ([ham_08_bch_jacobi.out.txt](outputs/ham_08_bch_jacobi.out.txt) (a) で、一般の関数 $f, g, h$ について sympy が 0 になることを確認) から $[f,[g,h]] - [g,[f,h]] = [f,[g,h]] + [g,[h,f]] = -[h,[f,g]] = [[f,g],h]$。したがって $[:f:, :g:]h = :[f,g]:\,h$。

**BCH の 3 次までの導出**。$X$、$Y$ を非可換な演算子として $\log(e^Xe^Y)$ を展開する。$e^Xe^Y = 1 + u_1 + u_2 + u_3 + \cdots$ ($u_n$ は $n$ 次)、

$$
u_1 = X + Y,\quad u_2 = \tfrac12X^2 + XY + \tfrac12Y^2,\quad u_3 = \tfrac16X^3 + \tfrac12X^2Y + \tfrac12XY^2 + \tfrac16Y^3
$$

($e^Xe^Y$ の積で、$X$ の項は必ず左、$Y$ の項は必ず右)。$\log(1 + u) = u - u^2/2 + u^3/3 - \cdots$ で次数ごとに集める。

- **2 次**: $u_2 - \tfrac12u_1^2 = \tfrac12X^2 + XY + \tfrac12Y^2 - \tfrac12(X^2 + XY + YX + Y^2) = \tfrac12(XY - YX) = \tfrac12[X, Y]$。
- **3 次**: $u_3 - \tfrac12(u_1u_2 + u_2u_1) + \tfrac13u_1^3$。まず積を展開する ($u_1u_2$ は、$(X + Y)$ の $X$ と $Y$ を、$u_2 = \tfrac12XX + XY + \tfrac12YY$ の各項の左にかけたもの。$u_2u_1$ は各項の右にかけたもの):

$$
\begin{aligned}
u_1u_2 &= \tfrac12XXX + XXY + \tfrac12XYY + \tfrac12YXX + YXY + \tfrac12YYY,\\
u_2u_1 &= \tfrac12XXX + \tfrac12XXY + XYX + XYY + \tfrac12YYX + \tfrac12YYY
\end{aligned}
$$

  $u_1^3 = (X + Y)^3$ は、8 つの語がすべて係数 1。語 (非可換な積。$XXY$ は $X$ を 2 回かけてから $Y$) ごとに、$u_3$、$u_1u_2$、$u_2u_1$、$u_1^3$ の係数を拾って $u_3 - \tfrac12(u_1u_2 + u_2u_1) + \tfrac13u_1^3$ を作る:

| 語 | $u_3$ | $u_1u_2$ | $u_2u_1$ | $u_1^3$ | 3 次の係数 |
|---|---|---|---|---|---|
| $XXX$ | $\tfrac16$ | $\tfrac12$ | $\tfrac12$ | 1 | $\tfrac16 - \tfrac12(\tfrac12 + \tfrac12) + \tfrac13 = 0$ |
| $XXY$ | $\tfrac12$ | 1 | $\tfrac12$ | 1 | $\tfrac12 - \tfrac12(1 + \tfrac12) + \tfrac13 = \tfrac1{12}$ |
| $XYX$ | 0 | 0 | 1 | 1 | $0 - \tfrac12(0 + 1) + \tfrac13 = -\tfrac16$ |
| $YXX$ | 0 | $\tfrac12$ | 0 | 1 | $0 - \tfrac12(\tfrac12 + 0) + \tfrac13 = \tfrac1{12}$ |
| $XYY$ | $\tfrac12$ | $\tfrac12$ | 1 | 1 | $\tfrac12 - \tfrac12(\tfrac12 + 1) + \tfrac13 = \tfrac1{12}$ |
| $YXY$ | 0 | 1 | 0 | 1 | $0 - \tfrac12(1 + 0) + \tfrac13 = -\tfrac16$ |
| $YYX$ | 0 | 0 | $\tfrac12$ | 1 | $0 - \tfrac12(0 + \tfrac12) + \tfrac13 = \tfrac1{12}$ |
| $YYY$ | $\tfrac16$ | $\tfrac12$ | $\tfrac12$ | 1 | $\tfrac16 - \tfrac12(\tfrac12 + \tfrac12) + \tfrac13 = 0$ |

  純粋な冪 $XXX$、$YYY$ の係数は 0 (3 次に $X^3$、$Y^3$ は残らない。$Y = 0$ なら $\log e^X = X$ だから当然)。語ごとの係数は sympy の非可換記号でも数えた ([ham_08_bch_jacobi.out.txt](outputs/ham_08_bch_jacobi.out.txt) (b))。表は [ham_10_audit_steps.out.txt](outputs/ham_10_audit_steps.out.txt) (a) で再計算した。
- 交換子で書く: $[X,[X,Y]] = X(XY - YX) - (XY - YX)X = XXY - 2XYX + YXX$、$[Y,[Y,X]] = YYX - 2YXY + XYY$。表の係数 ($XXY$、$YXX$、$YYX$、$XYY$ が $\tfrac1{12}$、$XYX$、$YXY$ が $-\tfrac16$、$XXX$、$YYY$ が 0) は、$\tfrac1{12}\bigl([X,[X,Y]] + [Y,[Y,X]]\bigr)$ の各語の係数 ($XXY$: $\tfrac1{12}$、$XYX$: $-\tfrac2{12}$、$YXX$: $\tfrac1{12}$、$YYX$: $\tfrac1{12}$、$YXY$: $-\tfrac2{12}$、$XYY$: $\tfrac1{12}$、$XXX$・$YYY$: 0) に一致する。
よって $\log(e^Xe^Y) = X + Y + \tfrac12[X,Y] + \tfrac1{12}([X,[X,Y]] + [Y,[Y,X]]) + \cdots$ (これが (C.14) の $X = :f:$、$Y = :g:$ の場合。上の交換子の関係で、演算子の交換子 $\to$ 関数の括弧 $:[f,g]:$ に直す)。

上のドリフトとキックの例 ($f = -\varepsilon p^2/2$、$g = -\varepsilon x^2/2$) で、$[f, g]/2$ までの $\exp$ が $\varepsilon^3$ の差を残し、3 次の項 $\frac{1}{12}([f,[f,g]] + [g,[g,f]])$ を足すとその差が消えることを、sympy で確認した ((f))。$[:f:, :g:] = :[f, g]:$ も、いくつかの多項式の $f$、$g$、$h$ で確認した。

**注意**: BCH の無限級数が実際に収束するか、どの次数まで足すか、は別の問題 (§9)。SAD のヘルプは、要素の変換を $\exp(:\cdot:)$ の積 (例えば DECA) で書き、次数をそろえるための補正項 $\exp(:V_4:)$ を足している (§8)。

### 6.3 $H$ が $s$ に依存する場合

(C.7) は $H$ が $s$ によらないとして導いた。$H$ が $s$ に依存すると、$\exp(-\int H\,ds\,:\,:)$ とは書けない。区間を細かく分けて各区間の積 $\exp(-\Delta s_1:H(s_1):)\exp(-\Delta s_2:H(s_2):)\cdots$ (§6.1 の順) を取り、分割を細かくした極限が正しい転送 (時間順序つき指数)。1 つの $\exp(:\Omega:)$ にまとめた形が**マグヌス展開**で、$\Omega = -\int H\,ds + \frac12\int\!\!\int[\ldots] + \cdots$ のように、ポアソン括弧の項が続く。この付録では、詳細は扱わない。区間ごとに $H$ が一定な (要素ごとの) 場合は、区間ごとの $\exp(-\Delta s:H:)$ の積 (§6.1 の順) になる。

## 7. 線形写像と 2 次の $f$、シンプレクティック性

$f$ が $z = (q_1, p_1, \ldots)$ の**2 次の斉次多項式**なら、$:f:\,z_k$ は $z$ の 1 次式なので、$:f:$ は行列で表せる。$H = \frac12\mathbf z^TA\mathbf z$ ($A$ は対称行列) で、(C.3) は $z_k' = \sum_l J_{kl}\,\partial H/\partial z_l$ ($J$ は $\begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$ の直和) だから、$\mathbf z' = JA\mathbf z$。次の 3 段で、$(-:H:)^n$ が行列のべきになることを示す。

1. (C.3) を $g = z_k$ に使うと、$-:H:\,z_k = z_k' = \sum_lJ_{kl}\,\partial H/\partial z_l = \sum_l(JA)_{kl}\,z_l$ ($\partial H/\partial z_l = (A\mathbf z)_l$)。
2. $-:H:$ は、関数の線形結合に線形に作用する (定数係数は外に出る)。1. の右辺は $z_l$ の 1 次式なので、$(-:H:)^2z_k = \sum_l(JA)_{kl}\,(-:H:\,z_l) = \sum_l(JA)_{kl}\sum_m(JA)_{lm}\,z_m = \sum_m[(JA)^2]_{km}\,z_m$。
3. $(-:H:)^nz_k = \sum_m[(JA)^n]_{km}\,z_m$ と仮定する。両辺に $-:H:$ を作用させて、$(-:H:)^{n+1}z_k = \sum_m[(JA)^n]_{km}\,(-:H:\,z_m) = \sum_m\sum_l[(JA)^n]_{km}(JA)_{ml}\,z_l = \sum_l[(JA)^{n+1}]_{kl}\,z_l$ (帰納法。$n = 1$ は 1.)。

(C.7) の級数に入れると、$\exp(-s:H:)\,z_k = \sum_n\frac{s^n}{n!}\sum_m[(JA)^n]_{km}z_m = \sum_m[\exp(sJA)]_{km}\,z_m$。したがって

$$
\exp(-s:H:)\,\mathbf z = \exp(sJA)\,\mathbf z = M\,\mathbf z,\qquad M = \exp(F),\;\;F = sJA
\tag{C.15}
$$

(確認: 2 自由度の一般の対称行列 $A$ で $n = 1\ldots5$ を sympy で確認。[ham_10_audit_steps.out.txt](outputs/ham_10_audit_steps.out.txt) (c)。)

5.2 の例 ($A = \mathrm{diag}(k, 1)$) では $F = \begin{pmatrix}0 & s\\ -ks & 0\end{pmatrix}$ で、行列指数は (C.10)(C.11) の行列に一致する (sympy (g))。**3 次以上の $f$ のリー変換は非線形写像**になる ($:f:$ が次数を保たず、(C.12) のように、高次の項を作る)。

**シンプレクティック性**。$M = \exp(F)$、$F = sJA$ とすると、$A$ が対称、$J^TJ = 1$、$J^2 = -1$ から $F^TJ + JF = s(AJ^TJ + J\,JA) = s(A - A) = 0$。したがって $\frac{d}{ds}\left(M^TJM\right) = M^T(F^TJ + JF)M = 0$ で、$s = 0$ では $M = 1$ なので $M^TJM = J$ ($M$ はシンプレクティック行列)。(5.2 の例で、$\det M = 1$、$M^TJM = J$ を確認。) 一般の (非線形の) $\exp(:f:)$ については、次の小節で証明する。

### 7.1 非線形な $\exp(:f:)$ も括弧を保つ (正準変換)

$f$ を任意の (形式的な) 関数とし、座標の関数 $a$、$b$ について次を示す。

$$
\exp(:f:)\,[a, b] = \bigl[\exp(:f:)\,a,\ \exp(:f:)\,b\bigr]
\tag{C.16}
$$

1. **$:f:$ は括弧のライプニッツ則を満たす**。ヤコビ恒等式 $[f,[a,b]] + [a,[b,f]] + [b,[f,a]] = 0$ で、$[b,f] = -[f,b]$、$[b,[f,a]] = -[[f,a],b]$ を使うと、$[f,[a,b]] = [a,[f,b]] + [[f,a],b]$。すなわち $:f:[a,b] = [:f:a,\ b] + [a,\ :f:b]$。
2. **$n$ 回の作用**: $:f:^n[a,b] = \sum_{k=0}^n\binom nk\,[:f:^ka,\ :f:^{n-k}b]$。$n = 0$ は自明。$n$ で成り立つと仮定して両辺に $:f:$ を作用させ、各項に 1. を使うと、$\sum_k\binom nk\bigl([:f:^{k+1}a,\ :f:^{n-k}b] + [:f:^ka,\ :f:^{n-k+1}b]\bigr)$。$[:f:^ka,\ :f:^{n+1-k}b]$ の係数は $\binom n{k-1} + \binom nk = \binom{n+1}k$ (パスカルの関係) で、$n + 1$ の式になる (帰納法)。
3. **指数に直す**: $\exp(:f:)[a,b] = \sum_n\frac1{n!}\sum_k\binom nk[:f:^ka,\ :f:^{n-k}b]$。$\frac1{n!}\binom nk = \frac1{k!\,(n-k)!}$ なので、$m = n - k$ として、$\sum_{k,m}\bigl[\frac{:f:^ka}{k!},\ \frac{:f:^mb}{m!}\bigr] = [\exp(:f:)a,\ \exp(:f:)b]$ (括弧は双線形)。これが (C.16)。
4. **座標に適用する**: $a = z_i$、$b = z_j$ とする。右辺は、新しい座標 $Z_i = \exp(:f:)z_i$ の括弧 $[Z_i, Z_j]$。左辺の $[z_i, z_j]$ は定数 ($0$ か $\pm1$) で、定数 $c$ には $:f:\,c = [f, c] = 0$ なので $\exp(:f:)c = c$。したがって

$$
[Z_i, Z_j] = [z_i, z_j]
\tag{C.17}
$$

ヤコビ行列を $M_{ik} = \partial Z_i/\partial z_k$ とすると、(C.1) から $[Z_i, Z_j] = (MJM^{\mathsf T})_{ij}$ ($J$ は $[z_i, z_j]$ を並べた行列)。(C.17) は $MJM^{\mathsf T} = J$ と同じで、$J^{-1} = -J$ を使って転置・逆行列を取ると $M^{\mathsf T}JM = J$ と同値。つまり**非線形の $\exp(:f:)$ のヤコビ行列は、各点でシンプレクティック行列**。リー変換の積 (§6.1 の順) も同じ性質を持つので、積で作る写像は、丸め誤差を除いてシンプレクティックに保たれる。(C.12) の $[x, p - k_2x^2] = 1$ は、この (C.17) の例。多項式の $f$、$a$、$b$ で、ライプニッツ則、二項展開、(C.16) の $\varepsilon^6$ までの展開、$[Z_x, Z_p] = 1$ を sympy で確認した ([ham_10_audit_steps.out.txt](outputs/ham_10_audit_steps.out.txt) (b))。

## 8. SAD のヘルプ、pysad・Fortran との対応

### 8.1 ヘルプに出てくる式

K64 版同梱のヘルプ原文 (`Documents/SADHelp.html`) で、$\exp(:\cdot:)$ が出てくるのは次の式 (原文で確認。DRIFT・QUAD・BEND・SEXT・OCT・DECA・DODECA。CAVI・MULT・SOL の節の式は 06・08・09・10 章):

| ヘルプの節 | 式 (ヘルプの記述の骨子) |
|---|---|
| DRIFT の transformation | $\exp(:H:)$、$H = (-\sqrt{p^2 - p_x^2 - p_y^2} + 1 - E/v_0)L$ ($E/v_0$ の符号は、ヘルプの式どうしの整合から誤植と判断した。gs での確認なし。02 章付録 A) |
| QUAD の本体 | スライスごとに $\exp(:H_{2n}:)$ ($H_{2n} = \frac1n\{(-p + \frac{p_x^2 + p_y^2}{2p} + \frac E{v_0})L + \frac{K1}{2}(x^2 - y^2)\}$)。スライスの間に $\exp(:\Delta H:)$ ($\Delta H = \frac1n(p - \sqrt{p^2 - p_x^2 - p_y^2} - \frac{p_x^2 + p_y^2}{2p})L$) |
| BEND の線形フリンジ | $\exp(:V:)$、平文版 $V = -\tfrac1p\bigl(\tfrac{F1^2}{24\rho_b}p_x + \tfrac{F1}{12\rho_b^2}y^2 - \tfrac{y^4}{6F1\rho_b^2}\bigr)$ (TeX 版は $y^2$・$y^4$ の項が違う誤植。05 章付録 A)。QUAD のフリンジは母関数 (04 章) |
| SEXT・OCT・DECA・DODECA の transformation | $\exp(:F_{\rm in}:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:V_n:)\exp(:aL:)\exp(:H_n/2:)\exp(:bL:)\exp(:F_{\rm out}:)$ ($L$ はドリフト、$H_n$ は薄い多極のキック $K_n/(n+1)!\,\mathrm{Re}(x - iy)^{n+1}$、$a = \frac12 - \frac1{\sqrt{12}}$、$b = \frac12 - a$、$V_n$ は 3 次を揃える補正)。原文の並びは前後対称でない (誤記。07 章付録 A) |

ヘルプの $H$ は、位置 $s$ の長さ $L$ (またはスライスの長さ) をかけた量 ($H = H_{\rm ham}\cdot L$) になっていて、そのうえで $\exp(:H:)$ と書かれている。

### 8.2 括弧の符号の食い違い (重要)

この付録の規約では、転送は $\exp(-s:H_{\rm ham}:)$、つまり $\exp(:f:)$ の $f$ は $-sH_{\rm ham}$ ((C.7))。ヘルプは、$f = +L\,H_{\rm ham}$ を $\exp(:f:)$ に入れている (例: DRIFT の $H$ は、ハミルトニアン $-\sqrt{p^2 - p_x^2 - p_y^2} + \ldots$ に $L$ をかけたもの)。ヘルプの $:f:$ は、この付録の $:f:$ と符号が逆、$:f:_{\rm SAD}\,g = -[f, g] = [g, f]$、すなわち **$\exp(:H:)$ は「$H$ を時間 1 のハミルトニアンとする流れ」** ($x \to x + \partial H/\partial p_x$、$p_x \to p_x - \partial H/\partial x$)。

**確認** (ヘルプ原文には括弧の定義が無い。式どうしの整合から読んだもので、gs での確認はしていない。[drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt) (1)、[ham_10_audit_steps.out.txt](outputs/ham_10_audit_steps.out.txt) (d))。ヘルプの $\exp(:H:)$ を「$H$ を時間 1 のハミルトニアンとする流れ」($\Delta x_i = \partial H/\partial p_i$、$\Delta p_i = -\partial H/\partial x_i$) と読んで、各式の偏微分を取る:

- (i) **DRIFT**: $H = \bigl(-\sqrt{p^2 - p_x^2 - p_y^2} + \cdots\bigr)L$ の $p_x$ 微分は、$\partial_{p_x}\bigl[-(p^2 - p_x^2 - p_y^2)^{1/2}\bigr] = -\tfrac12(\cdots)^{-1/2}\cdot(-2p_x) = p_x/\sqrt{p^2 - p_x^2 - p_y^2}$。よって $\Delta x = L\,p_x/p_z$ ($p_z = \sqrt{p^2 - p_x^2 - p_y^2}$) で、$x' = p_x/p_z$ が正しい向きに進む。逆の読み ($\Delta x = -\partial H/\partial p_x$) では $x$ が戻る。
- (ii) **QUAD**: $H_{2n}$ の座標を含む部分は $\dfrac{K1}{2n}(x^2 - y^2)$。$\Delta p_x = -\partial_xH = -\dfrac{K1}{n}\,x$、$\Delta p_y = -\partial_yH = +\dfrac{K1}{n}\,y$ ($n$ 個のスライスの合計で $\Delta p_x = -K1\,x$)。$K1 > 0$ で水平収束、垂直発散 (ヘルプの K1 の説明と一致)。
- (iii) **SEXT**: $\mathrm{Re}(x - iy)^3 = x^3 - 3xy^2$ ($(x - iy)^3 = x^3 - 3ix^2y - 3xy^2 + iy^3$ の実部)。$H_2 = \dfrac{K2}{3!}(x^3 - 3xy^2)$ から、$\Delta p_x = -\partial_xH_2 = -\dfrac{K2}6(3x^2 - 3y^2) = -\dfrac{K2}2(x^2 - y^2)$、$\Delta p_y = -\partial_yH_2 = -\dfrac{K2}6(-6xy) = K2\,xy$。
- (iv) **BEND の線形フリンジ** (平文版): $V = -\dfrac1p\Bigl(\dfrac{F1^2}{24\rho_b}\,p_x + \dfrac{F1}{12\rho_b^2}\,y^2 - \dfrac{y^4}{6F1\rho_b^2}\Bigr)$。$\Delta x = \partial V/\partial p_x = -\dfrac{F1^2}{24\,p\,\rho_b}$。$\Delta p_y = -\partial V/\partial y = \dfrac1p\Bigl(\dfrac{F1}{6\rho_b^2}\,y - \dfrac{2y^3}{3F1\rho_b^2}\Bigr) = \dfrac{y\,(F1^2 - 4y^2)}{6\,F1\,p\,\rho_b^2}$ (通分: $\dfrac{F1}{6} - \dfrac{2y^2}{3F1} = \dfrac{F1^2 - 4y^2}{6F1}$)。これがヘルプ自身の線形フリンジの式 ($x_2 = x_1 + \mathtt{dxfr}\,(p - p_0)/p$ の $p$ 依存部分 $-\mathtt{dxfr}/p$、$p_{y2} = p_{y1} + y\,(\mathtt{dyfr} - \mathtt{dyfra}\,y^2)/p_1^2$ の $y(\mathtt{dyfr} - \mathtt{dyfra}\,y^2)/p$ の形) と一致する ([drift_04_help_convention.out.txt](outputs/drift_04_help_convention.out.txt) (1))。
- (v) 04 章・07 章・03 章の数値 ($V_n$ の符号) とも整合 ([thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt))。

標準の括弧で $\exp(:H:)$ と読む ($\Delta x = [H, x] = -\partial H/\partial p_x$) と、上の (i)〜(iv) はすべて逆向きになる。唯一の例外が DRIFT の $E/v_0$ の符号 (QUAD・dynamics・BEND の body と逆。式どうしの整合から誤植と判断した。gs での確認なし)。

式を $\exp(:H_{\rm help}:)$ から本付録の記法に直すときは、$H_{\rm help} \to -H_{\rm help}$ とする。

順序の規約 (§6.1) については、ヘルプの DECA などの式が、入口のフリンジ $F_{\rm in}$ を最左、出口の $F_{\rm out}$ を最右に書いている (入口から出口の順)。これは §6.1 の規約 (左が先) と同じ向きに読める。ヘルプは順序の規約を明記していない (原文で確認)。ここでは、DECA などの式の並び (入口から出口の順) と整合する読みとして使う。

### 8.3 pysad・Fortran との対応

- **pysad `bendnl.py`** (偏向磁石の 3 次のリー生成関数を、追跡の写像から取り出す)。コードの説明に「薄い kick の生成関数 $f$ は、kick 写像 $z \to z + \tau(z)$ から $S\nabla f = \tau$ ($S = \begin{pmatrix}0 & 1\\ -1 & 0\end{pmatrix}$) で決まる」とある。$S\nabla f = (\partial_pf, -\partial_xf)$ なので、これはハミルトニアン $f$ の流れ ($\Delta x = \partial f/\partial p$、$\Delta p = -\partial f/\partial x$)、すなわちこの付録の記法では $\exp(:-f:)$ (ヘルプの $\exp(:f:)$ と同じ符号の向き)。`docs/25_resonance_terms_adts.md` の $M = R\,e^{:f_3:}e^{:f_4:}\cdots$ との符号の対応は、この付録では確かめていない。
- **K64 の Fortran** (四極の本体 `tsolqu`、薄い多極 など)。実装は、$\exp(:H_2:)$ などを級数ではなく**閉じた式** (四極の本体なら (4.9)(4.10)、キックなら $p \to p - \ldots$) で書いている。ヘルプの $\exp(:H:)$ を、コードの更新式にどう対応させたかは、各要素の章で書く。薄い多極のコード (OCT・DECA の $a$、$b$、$V$ の積の実装) との対応はこの付録では**未確認**。

## 9. 未確認・限界

- 級数 (C.6) (C.7) は**形式的なもの**として扱った (収束の議論はしていない)。ドリフト (C.8) や 6 極 (C.12) のように級数が止まる場合と、(C.10) のようにテイラー級数が既知の関数になる場合だけを例にした。
- BCH (C.14) は 3 次の項までをこの付録で導いた。5 次までの項 (自由リー代数の基底とヤコビ恒等式を使った展開) と、実際に使う分割の精度の議論は、[03_symplectic_integrators.md](03_symplectic_integrators.md) §4.1・§5・§8.3 に書いてある。
- ヘルプの $\exp(:\cdot:)$ の向きは、ヘルプに括弧の定義が無いので、式どうしの整合 (DRIFT・QUAD・SEXT・BEND の 4 つ) から読んだもの (gs での確認なし)。定義の明記は原文に無い (§8.2)。順序の規約も明記が無い。
- マグヌス展開 ($H$ が $s$ に依存する場合) は、名前の紹介だけで、式を導いていない。
- 多自由度 (結合、$x$-$y$) の例は書いていない。2 次の $f$ の一般的な対応 (C.15) は、1 自由度の例だけ確認した。
- この付録の 6 極の係数 $k_2$ は、§8.2 (iii) (ヘルプの SEXT の $H_2 = K2/3!\,\mathrm{Re}(x - iy)^3$ から $\Delta p_x = -(K2/2)(x^2 - y^2)$) と (C.12) ($\Delta p = -k_2x^2$、$y = 0$) を比べて $k_2 = K2/2$ に当たる。ヘルプの `K2` の定義の文面そのものとの照合はしていない。

## 参考

- SAD のヘルプ (K64 版同梱の `Documents/SADHelp.html`): DRIFT、QUAD、BEND、SEXT・OCT・DECA・DODECA などの transformation の節。
- [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md): 正準方程式の導出 (§5、§8)。
- [04_quad.md](04_quad.md): 四極の線形写像、$H_2$ と $\Delta H$ の分割 (ヘルプの $\exp(:H_{2n}:)$、$\exp(:\Delta H:)$ に当たる)。
- リー演算子の一般論の文献は、この付録では参照していない (内容は、ポアソン括弧と指数関数の定義から導いたもの)。
