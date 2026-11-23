# 03. シンプレクティックな分割の一般論

**この章で分かること**: K64 は、四極の本体も薄い多極も、「厳密に解ける 2 つの流れ」を交互に掛けて、ハミルトニアン全体の流れを近似する。この章は、その分割の一般論 — (1) ハミルトニアンを 2 つに分ける考え方、(2) 2 次の対称分割 (Strang) の誤差の導出、(3) 4 次への合成 (Yoshida)、(4) シンプレクティック性と長時間の振る舞い、(5) 薄い多極が使う**力勾配型**の 4 次の分割 (係数 $a, b, \beta, \gamma$ を導く)、(6) 刻みの数の決め方 — を、要素に依らない形で導く。四極と薄い多極の具体的な式は [04_quad.md](04_quad.md)、[07_thin_multipole.md](07_thin_multipole.md) に書き、この章はそれらが共通に使う道具を与える。式番号 $(3.n)$ から、§10 の Fortran の対応表と §11 の pysad の対応表を引ける。

**前提**: [C_poisson_lie.md](C_poisson_lie.md) (ポアソン括弧、$\exp(:f:)$、合成の順序、BCH 公式 (C.14))。ほかに [04_quad.md](04_quad.md) §6 (四極の Strang と Yoshida の結果) を読んでおくと、一般論と具体例の対応が分かる。

**読む順**: §1 (結論) → §3〜§5 (2 次・4 次) → §6〜§7 (性質) → §8 (力勾配型。薄い多極に必要) → §9 (刻みの数)。四極の章の図 ([図](figures/quad_integrator_order.png)) と数値 ([quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt)) は、この章の一般論の 1 つの例。この章の数値確認は、それとは別の簡単な系 (1 自由度) と、記号計算で行う。

## 1. 結論

- **分割の前提**: $H = H_a + H_b$ で、$H_a$ の流れ $A$ も $H_b$ の流れ $B$ も、閉じた式で厳密に解ける (K64 では、四極の本体なら $B$ = $H_2$ (2 次式) の線形な流れ、$A$ = $\Delta H$ (運動量だけ) の運動学的ドリフト。薄い多極なら $A$ = ドリフト、$B$ = 座標だけのポテンシャルのキック)。$A$ と $B$ は一般に**交換しない**ので、$e^{h(A+B)}$ は $e^{hA}e^{hB}$ に等しくなく、積の順序と長さの割り振りで近似の次数が決まる。
- **対称分割 (Strang)** $S_2(h) = e^{hA/2}\,e^{hB}\,e^{hA/2}$ は
$$
\log S_2(h) = h(A + B) + h^3\Bigl(-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]\Bigr) + O(h^5)
$$
 (式 (3.4))。偶数次の項が出ず、誤差は 1 刻み $h^3$、全体 $h^2$ (2 次)。
- **Yoshida の 4 次**: $S_4(h) = S_2(w_1h)\,S_2(w_0h)\,S_2(w_1h)$、$w_1 = 1/(2 - 2^{1/3}) = 1.351207$、$w_0 = 1 - 2w_1 = -1.702414$。$3$ 次が消えて誤差は $O(h^5)$ (全体 $h^4$)。$w_0 < 0$ (逆向きに進む部分がある)。
- **どの分割も、$A$・$B$ が厳密な流れ (正準変換) ならシンプレクティック**で、丸め誤差までそのまま保たれる (式 (3.8))。さらに、分割した写像は、元の $H$ に近い**別のハミルトニアン** $\tilde H = H + O(h^2)$ の厳密な流れになる (後方誤差解析。§6)。$H$ が $s$ によらないときの $H$ の誤差は、時間とともに増えずに振動する ([integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (c): 5 万刻みで Strang・Yoshida・力勾配型は比 1.0、RK4 は 7.7 倍に増え続ける)。
- **力勾配型の 4 次** (K64 の薄い多極 `tthin`、ヘルプの OCT・DECA・DODECA の transformation): $e^{aA}\,e^{B/2}\,e^{tA}\,e^{C}\,e^{tA}\,e^{B/2}\,e^{aA}$ ($a = \tfrac12 - t$、中央に**補正キック** $C$ — 座標だけの関数の流れ)。3 次の 2 つの条件から $t = 1/\sqrt{12}$ と補正の強さ $c_3 = \beta/2$ ($\beta = \tfrac16 - \tfrac{1}{2\sqrt{12}}$) が決まり、5 次の「運動量を含まない項」を消す条件から $c_5 = \gamma/2$ ($\gamma = \tfrac1{40} - \tfrac1{24\sqrt3}$) が決まる (§8)。これが K64 の `alpha`、`alpha1`、`beta`、`gamma` で、補正キックの関数形も K64 の式 (多極の 6 極、8 極) と記号計算で一致した ([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (c)(d))。力の評価は 1 刻みで 2 回 + 勾配 1 回 (Yoshida は 3 回)。
- **刻みの数**は要素ごとの経験式で、軌道の大きさに**よらず**設計軌道の強さだけで決まる (§9)。
- **使った条件と限界**は §12 に書く。

## 2. 記号と単位

| 量 | この章 | K64 の Fortran | pysad |
|---|---|---|---|
| 刻みの長さ | $h$ | `aln0`・`aln1`・`aln` (四極の 1 刻み)、`al` (薄い多極は 1 刻みで全体) | `aln1`、`L` |
| 流れの生成元 (リー演算子) | $A, B$ (§3) | | |
| 2 次の対称分割 | $S_2(h)$ | `tsolqu` の $D\,B\,D$ | `quadbody` の `group` (方式 0) |
| 4 次の合成の係数 | $w_1, w_0$ (式 (3.9)) | (K64 になし) | `_W1`、`_W0` (`quadbody.py`) |
| 力勾配型のドリフトの長さ | $a\,h$、$t\,h$ ($t = \alpha$) | `ala = al*alpha1`、`alb = al*alpha` | `ala`、`alb` (`_TT_ALPHA1`、`_TT_ALPHA`) |
| 補正キックの強さ | $c_3 = \beta/2$、$c_5 = \gamma/2$ | `beta`、`gamma` (`tquad.f:131-135` の `parameter`) | `_TT_BETA`、`_TT_GAMMA` |
| 精度パラメータ | $\varepsilon$ | `eps`・`eps0`・`epsdef` | `eps0` |
| 刻みの数 | $n$ | `ndiv` | `quadbody.ndiv()` |

**$\alpha, \beta, \gamma$ の注釈** (K64 `tquad.f:131-135`): `alpha` $= 1/\sqrt{12} = 0.288675\ldots$ (中央のドリフト 2 本の長さの比。数学の標準の記号ではない)、`alpha1` $= \tfrac12 - $ `alpha` $= 0.211325\ldots$ (両端のドリフトの比)、`beta` $= \tfrac16 - \alpha/2 = 0.0223291\ldots$ (補正キックの 3 次の強さ)、`gamma` $= \tfrac1{40} - \tfrac{1}{24\sqrt3} = 9.43739 \times 10^{-4}$ (補正キックの 5 次の強さ)。ソースのコメントは「`alpha=1/sqrt(12), beta=1/6-alpha/2, gamma=1/40-1/24/sqrt(3)`」。それぞれが何を消す係数かは §8 で導く。

## 3. ハミルトニアンの分割と、リー演算子の書き方

流れを「演算子の指数」で書く。[C_poisson_lie.md](C_poisson_lie.md) の規約で、ハミルトニアン $H$ (位置 $s$ によらない) の長さ $h$ の転送は $\exp(-h:H:)$。$A \equiv -:H_a:$、$B \equiv -:H_b:$ (演算子。$A, B$ は関数への作用で、$:f:\,g = [f, g]$)。$e^{hA}$ は $H_a$ の長さ $h$ の流れ、$e^{hB}$ は $H_b$ の流れ。積 $e^{hA}e^{hB}$ は「$e^{hA}$ が先、$e^{hB}$ が後」に粒子が通る順 (C.13 の規約)。回路の順は、対称な積では結果に影響しないので、以下の積の順序の議論はそのまま使える。

リー演算子は交換子で閉じる: $[:f:, :g:] = :[f, g]:$ (C.14 の前)。以下の $[A, B]$ は演算子の交換子 $AB - BA$ で、関数のポアソン括弧 $:[f_a, f_b]:$ に対応する。

K64 の分割の例:

| 要素 | $A$ | $B$ | 1 刻み | 次数 |
|---|---|---|---|---|
| 四極の本体 (追跡・CALC、[04_quad.md](04_quad.md)) | $\Delta H$ の流れ $D$ (運動学的ドリフト) | $H_2$ の厳密な線形の流れ | $D(h/2)\,B(h)\,D(h/2)$。v2 は 3 つを Yoshida で合成 | 2 次 (K64)、4 次 (v2) |
| 薄い多極 (追跡 `tthin`、[07_thin_multipole.md](07_thin_multipole.md)) | 厳密なドリフト | 多極のキック | 力勾配型 (§8) | 4 次 |
| 薄い多極 (CALC `tthine`) | 厳密なドリフト | 多極のキック | $L/6$、キック 1/2、$2L/3$、キック 1/2、$L/6$ (補正なし) | 2 次 (§8.4) |

## 4. 2 次の対称分割 (Strang)

$$
S_2(h) = e^{hA/2}\,e^{hB}\,e^{hA/2}
\tag{3.1}
$$

**BCH 公式** (C.14): $e^{X}e^{Y} = \exp\bigl(X + Y + \tfrac12[X, Y] + \tfrac1{12}([X,[X,Y]] + [Y,[Y,X]]) + \cdots\bigr)$。2 回使って (3.1) の $\log$ を 3 次まで出す。以下 $h = 1$ と置き、次数 (演算子 $A, B$ の個数) で整理する ($h$ は各項に次数の分だけ付く)。

**手順 1**: $X = A/2$、$Y = B$ とする。$Z_1 = \log\bigl(e^{A/2}e^{B}\bigr)$:

$$
Z_1 = \left(\tfrac A2 + B\right) + \tfrac12\left[\tfrac A2, B\right] + \tfrac1{12}\left(\left[\tfrac A2,\left[\tfrac A2, B\right]\right] + \left[B,\left[B, \tfrac A2\right]\right]\right) + \cdots
= \underbrace{\left(\tfrac A2 + B\right)}_{Z_1^{(1)}} + \underbrace{\tfrac14[A, B]}_{Z_1^{(2)}} + \underbrace{\tfrac1{48}[A,[A,B]] + \tfrac1{24}[B,[B,A]]}_{Z_1^{(3)}}
$$

**手順 2**: $e^{Z_1}e^{A/2}$ に BCH を使う ($X = Z_1$、$Y = A/2$)。$Z = Z_1 + \tfrac A2 + \tfrac12[Z_1, \tfrac A2] + \tfrac1{12}\bigl([Z_1,[Z_1, \tfrac A2]] + [\tfrac A2,[\tfrac A2, Z_1]]\bigr) + \cdots$。次数ごとに集める:

- **1 次**: $Z_1^{(1)} + \tfrac A2 = A + B$。
- **2 次**: $Z_1^{(2)} + \tfrac12[Z_1^{(1)}, \tfrac A2] = \tfrac14[A, B] + \tfrac12[B, \tfrac A2] = \tfrac14[A,B] - \tfrac14[A,B] = 0$ (前後対称なので偶数次が消える)。
- **3 次**: (i) $Z_1^{(3)}$ そのもの。(ii) $\tfrac12[Z_1^{(2)}, \tfrac A2] = \tfrac12\left[\tfrac14[A,B], \tfrac A2\right] = \tfrac1{16}[[A,B], A] = -\tfrac1{16}[A,[A,B]]$。(iii) $\tfrac1{12}[Z_1^{(1)},[Z_1^{(1)}, \tfrac A2]]$: $[Z_1^{(1)}, \tfrac A2] = [B, \tfrac A2] = \tfrac12[B,A]$、$[Z_1^{(1)}, \tfrac12[B,A]] = \tfrac12\left[\tfrac A2 + B, [B,A]\right] = \tfrac14[A,[B,A]] + \tfrac12[B,[B,A]] = -\tfrac14[A,[A,B]] + \tfrac12[B,[B,A]]$。$\tfrac1{12}$ をかけて $-\tfrac1{48}[A,[A,B]] + \tfrac1{24}[B,[B,A]]$。(iv) $\tfrac1{12}[\tfrac A2,[\tfrac A2, Z_1^{(1)}]] = \tfrac1{12}\cdot\tfrac14[A,[A, B]]\cdot$ ($[A/2, A/2] = 0$ なので $B$ の項だけ) $= \tfrac1{48}[A,[A,B]]$。

 $[A,[A,B]]$ の係数: $\tfrac1{48} - \tfrac1{16} - \tfrac1{48} + \tfrac1{48} = \tfrac{1 - 3 - 1 + 1}{48} = -\tfrac{2}{48} = -\tfrac1{24}$。$[B,[B,A]]$ の係数: $\tfrac1{24} + \tfrac1{24} = \tfrac1{12}$。したがって

$$
\log S_2(h) = h(A + B) + h^3\,C_3 + O(h^5),\qquad
C_3 = -\frac1{24}[A,[A,B]] + \frac1{12}[B,[B,A]]
\tag{3.2}
$$

(sympy で、非可換な自由代数を使って $\log S_2$ を直接計算し、$c_1 = -1/24$、$c_2 = +1/12$ を得た。2 次が消えることも確認: [integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (a)。)

**対称性から 5 次まで偶数次が出ない理由**: $S_2(h)S_2(-h) = 1$ (前後を入れ替えても同じ形)。$\log S_2(h) = \Omega(h)$ とすると $\Omega(-h) = -\Omega(h)$ だから、$\Omega$ は $h$ の奇数次だけ。

**誤差の見積もり**: $n = L/h$ 刻みで、$S_2(h)^n = \exp\bigl(L(A + B) + L h^2 C_3 + \cdots\bigr)$ (同じ形の積を $n$ 回かけると指数が $n$ 倍になる)。したがって**誤差は $L h^2 C_3$** で、刻みの数 $n$ を 2 倍にすると $1/4$ ([integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (a): Strang の比 4.0)。

四極の本体で $A = D$ (運動学的ドリフト)、$B = H_2$ の流れの場合は、$[B,[B,A]]$ が (4.19) の $\beta$ に、$[A,[A,B]]$ が $\alpha$ に当たる。

## 5. 4 次への合成 (Yoshida)

対称な 2 次の方法 $S_2$ を 3 つ、長さ $w_1h$、$w_0h$、$w_1h$ で並べる:

$$
S_4(h) = S_2(w_1h)\,S_2(w_0h)\,S_2(w_1h)
\tag{3.3}
$$

(3.2) を 3 回使って、$\log S_4$ を求める。3 つの指数の和の形になる (これら 3 つの $\Omega(w_ih) = w_ih(A+B) + w_i^3h^3C_3 + O(h^5)$ の積の BCH では、交換子 $[\Omega_i, \Omega_j]$ のうち $(A + B)$ どうしは 0 (同じ $A + B$ の定数倍どうしは交換)。$(A + B)$ と $h^3C_3$ の交換子 ($h^4$) は残りうるが、積全体が前後対称なので偶数次は消え、残りは $h^5$ 以上):

$$
\log S_4 = (2w_1 + w_0)\,h(A+B) + (2w_1^3 + w_0^3)\,h^3 C_3 + O(h^5)
\tag{3.4}
$$

(対称な積なので偶数次は出ない。) 正しい $h(A+B)$ の係数 1 と、$h^3$ の消去から

$$
2w_1 + w_0 = 1,\qquad 2w_1^3 + w_0^3 = 0
\tag{3.5}
$$

$w_0 = 1 - 2w_1$ を入れ、$(1 - 2w_1)^3 = -2w_1^3$ の実数の立方根 $1 - 2w_1 = -2^{1/3}w_1$ から

$$
w_1 = \frac{1}{2 - 2^{1/3}} = 1.3512071919\ldots,\qquad w_0 = 1 - 2w_1 = -1.7024143839\ldots
\tag{3.6}
$$

(四極の章の (4.20)〜(4.22) と同じ。)

**一般化 (高次への飛び越し)**: 対称な $p$ 次 ($p$ は偶数) の方法 $S$ の誤差が $c\,h^{p+1}$ のとき、$S(w_1h)S(w_0h)S(w_1h)$ の誤差は $(2w_1^{p+1} + w_0^{p+1})\,c\,h^{p+1}$。$2w_1 + w_0 = 1$ と $2w_1^{p+1} + w_0^{p+1} = 0$ から $w_1 = 1/(2 - 2^{1/(p+1)})$。$p = 2$ が (3.6)、$p = 4$ なら $1/(2 - 2^{1/5})$ で 6 次の方法ができる。K64・pysad は 4 次までしか使わない。

**確認**: 自由代数で、$w_1$ を入れた $S_4$ の $\log$ の 1 次は $A + B$ ($2 w_1 + w_0 = 1$)、2・3・4 次は 0 (係数 $10^{-41}$ 以下)、5 次は残る (語の数 32、最大の係数 0.10)。5 次の関数は、具体的な表現 ($A = -p^2/2$、$B = -V(x)$) で

$$
H_5^{\rm Yoshida} = 0.05209\,V'^2V''\,p^0 + \left(-0.01199\,V'V''' + 0.01405\,V''^2\right)p^2 + 0.000414\,V''''\,p^4
\tag{3.7}
$$

([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (b)。$V' = V_x$ など)。**4 次の誤差は、これらの項の和で決まる**。

## 6. シンプレクティック性と長時間の振る舞い

- **シンプレクティック性**: $e^{hA}$、$e^{hB}$ はそれぞれ正準変換 ($A$、$B$ の流れが厳密に解けるので、ヤコビアン $M$ が $M^\top JM = J$ を厳密に満たす。[C_poisson_lie.md](C_poisson_lie.md) §7)。正準変換の積も正準変換なので、**分割した写像は、刻みの数・次数・係数によらず、丸め誤差までシンプレクティック**:

$$
M_{\rm split} = M_n\cdots M_1,\qquad M_i^\top J M_i = J\ \Rightarrow\ M_{\rm split}^\top J M_{\rm split} = J
\tag{3.8}
$$

 (係数が負 ($w_0 < 0$) でも、逆向きの流れ $e^{-|w|hA}$ も正準変換なので同じ。) 比較の RK4 は、1 刻みで $\det J - 1 = -1.4 \times 10^{-4}$ (1 自由度で $h = 0.3$)、Strang・Yoshida・力勾配型は $10^{-11}$ の程度 (中心差分の精度 $10^{-6}$ による。[integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (b))。
- **後方誤差解析 (modified Hamiltonian)**: $S_2(h) = \exp\bigl(h(A + B) + h^3C_3 + \cdots\bigr)$ は、ある $\tilde H(h)$ の厳密な流れ $\exp(-h:\tilde H:)$ (形式的に。$\tilde H = H + h^2H_3 + \cdots$) になる。$\tilde H$ は $s$ によらないので、$\tilde H$ は分割した写像で**保存**される。元の $H$ は、$\tilde H - H = O(h^2)$ だけ、位相に応じて振動する。したがって、$s$ によらない $H$ の流れを分割で長く追跡したときの $H$ の誤差は、**時間とともに増えず、$O(h^2)$ (4 次では $O(h^4)$) の幅で振動する**。(厳密には、級数 $\tilde H$ は一般に発散する漸近級数で、「指数的に長い時間」の範囲の話。)
- **数値**: [integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (c) ($H = p^2/2 + x^2/2 + x^4/4$、$h = 0.25$、5 万刻み): 最初の 5000 刻みと全体の $\max|\Delta H|$ の比が、Strang で 1.0 ($2.6 \times 10^{-3}$)、Yoshida で 1.0 ($1.3 \times 10^{-4}$)、力勾配型で 1.0 ($3.6 \times 10^{-6}$)。RK4 は 7.2×10⁻³ から 5.5×10⁻² へ 7.7 倍に増える (系統的なずれ。非シンプレクティック)。**多数周回する蓄積リングの追跡では、この性質が効く** (周回を重ねても、ずれが蓄積しない)。
- 注意: $s$ に陽に依る $H$ (たとえば加速、RF の時間依存) では、$H$ は保存しない。シンプレクティック性そのもの ((3.8)) は保たれる。

## 7. 一般論のまとめ: 何が分割の「良さ」を決めるか

誤差の次数は、BCH の交換子の何次の項まで消せるかで決まる。2 つの生成元の 3 次の交換子は 2 つ ($[A,[A,B]]$、$[B,[B,A]]$)、5 次は (自由代数で) 6 つ。分割の自由度 (ドリフトの長さ・キックの強さの割り振り) を使って、低い次数から消す。座標だけの関数の補正キック $C$ を使うと、$[B,[B,A]]$ の項 (**力の勾配** $\propto |\nabla V|^2$ の項) をキックで打ち消せる。これが**力勾配型**で、次の §8。

## 8. 力勾配型の 4 次の分割

### 8.1 構成

$A$ = ドリフトの流れ ($H_a = T = P^2/(2m)$)、$B$ = 座標だけのキック ($H_b = V(\vec x)$)。K64 の薄い多極の 1 刻みは、(中央の補正キックを $e^{C}$ として)

$$
S_{\rm FG}(h) = e^{a h A}\,e^{hB/2}\,e^{t h A}\,e^{C}\,e^{t h A}\,e^{hB/2}\,e^{a h A},\qquad a = \tfrac12 - t
\tag{3.9}
$$

(ドリフトの長さの合計 $2ah + 2th = h$、キックの合計 $h/2 + h/2 = h$。前後対称。)

補正キックは、$B$ が座標だけの関数なので、座標だけの関数の流れになる。自由代数の元として、$C = c_3\,W + c_5\,W_5$ ($h$ のべきを含める)、

$$
W \equiv [B,[A,B]]\ (3\text{ 次}),\qquad W_5 \equiv [W, [A,B]]\ (5\text{ 次})
\tag{3.10}
$$

を取る。具体的な表現 ($A = -:T:$、$B = -:V:$) では、$[A, B]$ は $-(P/m)\cdot\nabla V$ の流れ、$W$ は $|\nabla V|^2/m$ (座標だけの関数)、$W_5$ は $-2\,V_jV_kV_{jk}/m^2$ (座標だけの関数。$V_j = \partial V/\partial x_j$、$V_{jk}$ は 2 階微分、重複する添字は和) に比例する。**座標だけの関数** $f(\vec x)$ の流れは、運動量を $\Delta p_j = \partial f/\partial x_j$ だけ変える (座標は不変) キックで、閉じた式。

### 8.2 3 次の条件 (手で)

(3.9) を、Strang の公式 (3.2) を 2 回使って組む。

**内側** $e^{hB/2}\,e^{tA + C/h}\cdots$: 中央の $e^{tA}e^{C}e^{tA} = \exp(2tA + C + O(h^5))$ ($C$ は 3 次、$[A, C]$ は 4 次の項になるが、前後対称な積なので偶数次が消えて、5 次以上に回る)。これを $e^{B/2}\,e^{2tA + C}\,e^{B/2}$ ((3.1) で $A \to B$、$B \to 2tA + C$、$h = 1$ の形) に入れ、(3.2) を使う:

$$
Y = B + (2tA + C) - \tfrac1{24}[B,[B, 2tA]] + \tfrac1{12}[2tA,[2tA, B]] + O(5)
$$

($C$ は 3 次なので、(3.2) の 3 次の交換子の中では 5 次以上。) $[B,[B, 2tA]] = 2t[B,[B,A]]$、$[2tA,[2tA, B]] = 4t^2[A,[A,B]]$、$C = c_3[B,[A,B]] = -c_3[B,[B,A]]$ なので

$$
Y = B + 2tA + \Bigl(\tfrac{t^2}3\,[A,[A,B]] - \bigl(\tfrac{t}{12} + c_3\bigr)[B,[B,A]]\Bigr) + O(5)
$$

**外側** $e^{aA}\,e^{Y}\,e^{aA}$ ((3.1) で $A \to 2aA$、$B \to Y$): $X \equiv 2aA$ として、$\exp\bigl(X + Y - \tfrac1{24}[X,[X,Y]] + \tfrac1{12}[Y,[Y,X]]\bigr)$。3 次の項だけを集める ($Y$ の 1 次部分 $B + 2tA$ を交換子の中に入れる)。$[X,[X,Y]] = 4a^2[A,[A,B]]$ ($2tA$ は $A$ と交換)。$[Y,[Y,X]] = 2a\,[B + 2tA,[B, A]] = 2a[B,[B,A]] - 4at[A,[A,B]]$。したがって、$\log S_{\rm FG}$ の 3 次の項の係数は

$$
\begin{aligned}
[A,[A,B]]:&\quad \tfrac{t^2}3 - \tfrac{a^2}6 - \tfrac{at}3 = \tfrac{12t^2 - 1}{24}\quad(a = \tfrac12 - t\text{ を代入}),\\
[B,[B,A]]:&\quad -\tfrac t{12} - c_3 + \tfrac a6 = \tfrac{1}{12} - \tfrac t4 - c_3
\end{aligned}
\tag{3.11}
$$

(sympy の自由代数による直接計算で、同じ式 $c_1 = (12t^2 - 1)/24$、$c_2 = -(12c_3 + 3t - 1)/12$ を得た。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (c)。) **3 次が消える条件**は

$$
12t^2 = 1\ \Rightarrow\ t = \frac1{\sqrt{12}} = 0.288675\ldots,\qquad
c_3 = \frac1{12} - \frac t4 = \frac1{12} - \frac1{4\sqrt{12}} = 0.0111645\ldots
\tag{3.12}
$$

($t > 0$ の解。$t = -1/\sqrt{12}$ は負のドリフトで採らない。) これは、K64 の `alpha = 1/sqrt(12)`、`alpha1 = 1/2 - alpha` ($= a$)、そして `beta` $= \tfrac16 - \tfrac\alpha2 = 0.0223291$ の**半分**: $c_3 = \beta/2$ (数値で $c_3 = 0.011164549685$、$\beta/2 = 0.011164549685$)。

**意味**: ドリフトの長さの割り振り ($t = 1/\sqrt{12}$) は、$[A,[A,B]]$ (運動量の項。$p^2V'''$ 型) を消す。補正キック $c_3 W$ は、$[B,[B,A]]$ ($|\nabla V|^2$ 型の力の勾配の項) を消す。補正キックは、座標だけの関数の流れで、**$B$ が座標だけの関数だから閉じた式で作れる**。

### 8.3 5 次の条件 ($c_5$)

3 次を消した後の 5 次の項 ($t$、$c_3$ を代入) は、具体的な表現 ($A = -p^2/2$、$B = -V$。1 自由度) で

$$
H_5^{\rm FG} = \frac{-720\,c_5 + 9 - 5\sqrt3}{360}\,V'^2V''\,p^0 + \left[\Bigl(\tfrac{\sqrt3}{288} - \tfrac1{120}\Bigr)V'V''' + \Bigl(\tfrac1{240} - \tfrac{\sqrt3}{288}\Bigr)V''^2\right]p^2 + \frac{V''''}{4320}\,p^4
\tag{3.13}
$$

(自由代数の $\log S_{\rm FG}$ の 5 次の項を、Dynkin の写像で関数に直したもの。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (c)。)

- $p^2$、$p^4$ の項は、$C$ に含まれる座標だけの関数ではどうにもならない (運動量を含むので)。これは 4 次の誤差として残る。
- **$p^0$ の項 ($V'^2V''$) は、座標だけの関数なので、$c_5$ を使って消せる**: $-720c_5 + 9 - 5\sqrt3 = 0$ より

$$
c_5 = \frac{9 - 5\sqrt3}{720} = \frac1{80} - \frac{\sqrt3}{144} = 4.71869 \times 10^{-4} = \frac{\gamma}{2}
\tag{3.14}
$$

 (K64 の `gamma` $= \tfrac1{40} - \tfrac1{24\sqrt3} = 9.43739 \times 10^{-4}$ の半分。$\tfrac1{24\sqrt3} = \tfrac{\sqrt3}{72}$。比 $c_5/\gamma = 0.500000000$。)

多自由度では、$p^0$ の項の座標だけの関数は $G \equiv V_jV_kV_{jk}$ (5 次のただ 1 つの座標だけの交換子 $W_5$。$[B,[B,[A,[A,B]]]]$ も、ヤコビ恒等式で $W_5$ に比例する) だけ。1 自由度の式 ($V'^2V''$) はその特別な場合。**この $\gamma$ の項は 5 次の誤差の一部 ($p^0$) だけを消す**ので、全体の誤差が小さくなるとは限らない: この章の例 (1 自由度の $x^4/4$) では、$c_5$ を入れない方 (c3 だけ) が誤差が小さい ([integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (a)(d)。$n = 64$ で $1.5 \times 10^{-9}$ 対 $6.5 \times 10^{-9}$)。多極 ($V \propto \mathrm{Re}\,(x - iy)^{n+1}$) では、07 章で確認する。

### 8.4 補正キックの具体形と K64 との一致

質量に当たるものが全運動量 $p$ ($T = (p_x^2 + p_y^2)/(2p)$。$A$ の $H_a$ は $p$ に依る) のとき、補正キックの母関数 ($\exp(:f:)$ の $f$。$\Delta p_j = \partial f/\partial x_j$、$\Delta z = -\partial f/\partial p$。[C_poisson_lie.md](C_poisson_lie.md) の規約) は

$$
f_C = c_3\,h^3\,\frac{|\nabla V|^2}{p} - 2c_5\,h^5\,\frac{V_jV_kV_{jk}}{p^2}
\tag{3.15}
$$

(1 自由度で $m = p$ とすると、$W = V'^2/m$、$W_5 = -2V'^2V''/m^2$ で $f_C = c_3 W + c_5 W_5$)。規格化した運動量 $q = p_x/p$ の更新は $\Delta q = \Delta p_x/p$。$V = \dfrac{K}{(k+1)!}\mathrm{Re}\,(x - iy)^{k+1}$ (多極のポテンシャル。07 章) として、(3.15) から作った $\Delta q_x$、$\Delta q_y$、$\Delta z$ は、K64 の `tthin` の補正キックの式 (`tquad.f:208-250` の $f_1$〜$f_5$ を使う式) と、$k = 2$ (6 極)、$k = 3$ (8 極) の**記号計算で完全に一致した** (差 0。$K$、$p$ は記号のまま。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (d))。具体形は [07_thin_multipole.md](07_thin_multipole.md) で書く。

**CALC の `tthine` の分割 (2 次)**: $L/6$、キック 1/2、$2L/3$ (= $t = 1/3$ の 2 本)、キック 1/2、$L/6$ で、補正 $C = 0$。(3.11) に $t = 1/3$、$c_3 = 0$ を入れると、$[B,[B,A]]$ の係数は $\tfrac1{12} - \tfrac1{12} = 0$ (消える)、$[A,[A,B]]$ の係数は $\tfrac{12/9 - 1}{24} = \tfrac1{72}$ (残る)。つまり**補正なしで $[B,[B,A]]$ だけを消した 2 次の分割**で、Strang の誤差の係数 $|-\tfrac1{24}|$ と $\tfrac1{12}$ より小さい ([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (e))。

## 9. 数値確認 (一般論の系で)

[integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt)。$H = p^2/2 + x^2/2 + x^4/4$ (1 自由度)、$(x, p) = (0.5, 0.3)$ から $s = 2$。誤差は DOP853 ($rtol = 10^{-13}$) との差:

| 方法 | $n = 8$ | $n = 16$ | $n = 32$ | $n = 64$ | 比 (最後) | 力の評価/刻み |
|---|---|---|---|---|---|---|
| Strang | $8.3 \times 10^{-3}$ | $2.1 \times 10^{-3}$ | $5.1 \times 10^{-4}$ | $1.3 \times 10^{-4}$ | 4.0 | 1 |
| Yoshida | $1.3 \times 10^{-3}$ | $7.7 \times 10^{-5}$ | $4.8 \times 10^{-6}$ | $3.0 \times 10^{-7}$ | 16.0 | 3 |
| 力勾配型 ($c_3$, $c_5$) | $2.7 \times 10^{-5}$ | $1.7 \times 10^{-6}$ | $1.0 \times 10^{-7}$ | $6.5 \times 10^{-9}$ | 16.0 | 2 + 勾配 1 |
| 力勾配型 ($c_3$ だけ) | $6.3 \times 10^{-6}$ | $3.9 \times 10^{-7}$ | $2.5 \times 10^{-8}$ | $1.5 \times 10^{-9}$ | 16.0 | 2 + 勾配 1 |
| RK4 (非シンプレクティック) | $8.5 \times 10^{-5}$ | $5.3 \times 10^{-6}$ | $3.3 \times 10^{-7}$ | $2.1 \times 10^{-8}$ | 16.1 | 4 |

力の評価の総数を約 64 回にそろえた誤差 (同 (d)): Strang $1.3 \times 10^{-4}$、Yoshida $2.6 \times 10^{-5}$、力勾配型 $5.6 \times 10^{-7}$ ($c_3$ だけは $1.3 \times 10^{-7}$)、RK4 $5.3 \times 10^{-6}$。この系では、力勾配型は力の評価あたりの誤差が一番小さい (勾配の評価を 1 回の力の評価と数えた)。**力勾配型には力の勾配 $\nabla|\nabla V|^2$ の式が要る**ので、$V$ が多項式 (多極) のように勾配が簡単に書ける場合の方法。

## 10. 刻みの数の決め方 (K64 の式の一覧)

刻みの数は、要素ごとの経験式。**いずれも設計軌道の強さ (と放射の強さ) だけで決まり、粒子の振幅には依らない**。

| 場所 | 式 | $\varepsilon$ | 備考 |
|---|---|---|---|
| 四極の追跡 (`tsolqu`、`tsolqu.f:47`) | $n = 1 + \lfloor L\,\lVert(K, B_z)\rVert/\varepsilon\rfloor$ (`ak` = K1 の積分強度、`bz` = ソレノイド) | $\varepsilon = 0.2\,\mathrm{EPS}$ (`epsdef = 0.2d0`、`tsolqu.f:14`) | 04 章 §6.3 |
| 四極の CALC (`tsolque`、`tsolque.f:385`) | 同じ式。上限 `ndivmax` | $0.1\,\mathrm{EPS}$ (`tsolque.f:381-383`) | 追跡の半分の $\varepsilon$ |
| 放射つきの四極の追跡 (`tsolqu.f:400-402`) | $n = \max\bigl(1 + \lfloor L\,K/\varepsilon\rfloor,\ \mathrm{ndivrad}(\lVert(a_{k0x}, a_{k0y})\rVert, K, B_z, \mathrm{EPS})\bigr)$、上限 `ndivmax` | $0.2\,\mathrm{EPS}$ | `ak0x`、`ak0y` はステアリング (偏向の積分強度) |
| 偏向磁石 (K1 ≠ 0) の CALC (`tbende.f:460`) | $n = 1 + \lfloor\lvert\phi_c/\varepsilon\rvert\rfloor$ ($K1 = 0$ なら 1) | `epsbend` $\times$ EPS (定義は未確認) | $\phi_c$ は偏向角に当たる量 (定義は未確認)。放射・スピン偏極があれば、ほかの条件との大きい方 |
| 放射つきの偏向磁石の追跡 (`tbend.f:640`) | $n = \min(\mathrm{ndivmax}, \mathrm{ndivrad}(\phi_r, 0, 0, \varepsilon))$ | EPS | $\phi_r = \phi_b\,(L - f_{1r} - f_{2r})/L$ (端の分を除いた偏向角) |
| 多極 MULT の追跡 (`tmulti.f:107-124`) | $n = \max(\mathrm{ndiv0}, \mathrm{ndivrad}(\ldots), \ldots)$、上限 `ndivmax` | EPS | `ndiv0` は要素のパラメータ (名前の対応は未確認) |
| 薄い多極 (`tthin`) | $n = 1$ (本体は 1 刻み) | — | `call setndivelm(l_track, 1)` |
| `ndivrad(phi, ak, bz, eps)` (`tffs.f:641-653`) | $\lfloor 1 + (h_0/100)\,[(\lvert\phi\rvert + \lVert(K, B_z)\rVert\cdot 0.02)/\varepsilon_1]\rfloor$ ($\varepsilon_1 = $ EPS (EPS ≤ 0 なら 1)。`trpt` が真なら $\times 10$) | EPS | $h_0 = \gamma_0$ (`tmacro` の `h0`)。係数 0.02 は `arad`。放射の光子数が 1 刻みで多くなりすぎないための分割数と思われる (推定) |

**式の読み方**: 四極なら $n \propto L\,K1$。$L\,K1 = KL^2 = \varphi^2$ は本体全体の位相の 2 乗 ($\varphi = wL$、$K$ は長さあたりの強さ)。1 刻みの位相は $\varphi/n \approx \sqrt{\varepsilon/n}$ (近似)。誤差は 1 刻みの位相の何乗かに比例するので、$\varepsilon$ を決めることが 1 刻みの位相を決めることになる。

数値例は [04_quad.md](04_quad.md) §6.3 の表 ($L = 0.5$ m、K1 = 1.5 /m で追跡 K64 は $n = 4$、v2 は 8)。

## 11. 数値確認で使うスクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `integ_01_bch_force_gradient.py` | (a) Strang の BCH (非可換な級数、sympy)。(b) Yoshida の係数と 3 次の消去、5 次の残り。(c) 力勾配型: 3 次の係数 ((3.11))、$t, c_3$ ((3.12))、$c_5$ ((3.14))、残る 5 次の項 ((3.13))。(d) 2 自由度で、導いた補正キックと K64 の `tthin` の式 ($\Delta q$、$\Delta z$) が 6 極・8 極で記号的に一致。(e) CALC の `tthine` の分割の係数 | [integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) |
| `integ_02_order_symplectic.py` | 誤差の次数 (Strang・Yoshida・力勾配型・RK4)、シンプレクティック性、長時間の $H$ の誤差、力の評価の数あたりの誤差 | [integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) |

実行はリポジトリのルートから `py -3.12 physics_guide/scripts/integ_01_bch_force_gradient.py` (sympy。約 12 秒)、`py -3.12 physics_guide/scripts/integ_02_order_symplectic.py` (約 1 秒)。四極の積分器の次数の確認は [quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (この章とは重複しない)。

## 12. Fortran・pysad の対応

| 内容 | Fortran (ファイル:行) | pysad |
|---|---|---|
| 四極の Strang の構造 ($D\,B\,D$、隣り合う $D$ を 1 つに) | `tsolqu.f:61-95` (本章の範囲では刻みの構造だけ読んだ。式は 04 章) | `quadbody.group` (`quadbody.py:76`)、`ngroup` |
| Yoshida の係数 | (K64 に無い) | `_W1`、`_W0` (`quadbody.py:41-42`)、`QB_BW` (`:44`) |
| 方式の選択 (k64 は Strang、v2 は Yoshida) | — | `quadbody.scheme()`、`set_scheme` (`quadbody.py:55-66`) |
| 力勾配型の係数 | `tquad.f:131-135` (`alpha`、`beta`、`gamma`、`alpha1`) | `track6._TT_ALPHA`・`_TT_ALPHA1`・`_TT_BETA`・`_TT_GAMMA` |
| 力勾配型の本体 (追跡) | `tquad.f:118-` (`tthin`)。`ala = al*alpha1`、`alb = al*alpha` (`:173-174`)。補正キックは `kord = 2` が `:208-`、`kord = 1` が `:239-`、一般の `kord` が `:261-` | `track6._tthin_sext` (`track6.py:892`)、`_tthin_gen` (`:740`) |
| CALC の 2 次の分割 | `tthine.f:36-87` (`ala = al/6`、`alb = al/1.5`) | `lattice._tthine_cod` (`lattice.py:1755`)、`elkern._tthine_map` (`:360`) |
| 刻みの数の式 | §10 | `quadbody.ndiv`・`ndiv_calc` (`quadbody.py:59-71`) |

(行番号は参照した K64 のソース、pysad のファイルのもの。)

## 13. 未確認・限界・使った条件

**使った条件**: 一般論の数値確認は、1 自由度 $H = p^2/2 + x^2/2 + x^4/4$ と、2 自由度の多極のポテンシャル ($k = 2, 3$) で、記号計算 (sympy) と数値積分 (DOP853)。K64 のソースとの比較は、`tthin` の補正キックの式を記号計算の結果と比べたもの (pysad の追跡を数値で動かした比較は 07 章)。

**未確認・限界**:
- **BCH 級数の収束**は議論していない。誤差の次数は、形式的な級数の係数の次数で決まる (数値で確認)。$\tilde H$ の「指数的に長い時間」の保存は、解析的な定理を引用しただけで、この章では確かめていない (数値は 5 万刻みまで)。
- **刻みの数の経験式**は、K64 のソースの式を読んだもの。$\varepsilon$ の値の物理的な根拠 (なぜ 0.2、0.1 か) は不明 (K64 の定数)。`epsbend`、`epsrad`、`ndivmax`、`ndiv0`、$\phi_c$ の定義は未確認。
- **力勾配型の 5 次の全誤差**の評価は、1 自由度の $p^0$ の項だけを示した。$p^2$・$p^4$ の項の大きさ (3.13) と、多極での実際の誤差は 07 章で扱う。
- **5 次で消せる項が「$p^0$ の項だけ」**であること (多自由度で $G$ のみ) は、1 自由度の自由代数の計算と Jacobi 恒等式の議論による。多自由度の完全な基底の分解は計算していない。
- **Strang の誤差を $A, B$ のどちらを外側にするかで比較する議論**は、この章では扱っていない (四極の章の構成は $A = D$ が外側)。
- 加速のある要素・時間依存の $H$ への拡張は扱わない (シンプレクティック性だけが保たれる)。

## 14. この章で使う Fortran の変数名

SAD の変数名の注釈の規約は [01_variables_and_conventions.md](01_variables_and_conventions.md) §8。

| Fortran の名前 | 式 | 定義の場所 (Fortran) | pysad での名前 |
|---|---|---|---|
| `alpha` | $\alpha = 1/\sqrt{12} = 2.88675134594812882 \times 10^{-1}$。(3.12) の $t$。中央のドリフト 2 本のそれぞれの長さ $= \alpha\,L$ | `tquad.f:132` (`parameter`) | `_TT_ALPHA` |
| `alpha1` | $\tfrac12 - \alpha$。両端のドリフトの長さ $= \alpha_1 L$ ((3.12) の $a$) | `tquad.f:135` | `_TT_ALPHA1` |
| `beta` | $\beta = \tfrac16 - \alpha/2 = 2.23290993692602255 \times 10^{-2}$。補正キックの 3 次の強さ ($c_3 = \beta/2$) | `tquad.f:133` | `_TT_BETA` |
| `gamma` | $\gamma = \tfrac1{40} - \tfrac1{24\sqrt3} = 9.43738783765593145 \times 10^{-4}$。補正キックの 5 次の強さ ($c_5 = \gamma/2$) | `tquad.f:134` | `_TT_GAMMA` |
| `ala`、`alb` | $\alpha_1 L$、$\alpha L$ (`tthin` の両端、中央の各ドリフトの長さ) | `tquad.f:173-174` | `ala`、`alb` |
| `f1`〜`f5` | 補正キックの係数。$k = 2$ (6 極) で $f_1 = 2\beta L$、$f_2 = 4\gamma L^2$、$f_3 = 6\gamma L^2$、$f_4 = \tfrac12\beta L$、$f_5 = f_2$。一般の $k$ では (07 章) | `tquad.f:209-216`、`:239-245`、`:261-267` | `_tthin_sext` の `f1`〜`f5` |
| `aki` | $\mathrm{ak}/(k!\,p)$ (補正キックの中。`akf/pr*2`)。1 回のキックでは `akf/pr` ($\mathrm{akf} = \mathrm{ak}/k!/2$) | `tquad.f` | `aki` |
| `eps`、`eps0`、`epsdef` | $\varepsilon = \mathrm{epsdef}\times\mathrm{eps0}$ (`eps0` = 要素の EPS。0 なら `epsdef`) | `tsolqu.f:14、43-46` | `eps0`、`quadbody.ndiv` の `eps` |
| `ndiv` | 刻みの数 | `tsolqu.f:47` ほか | `quadbody.ndiv()` |
| `ndivrad` | 放射による分割数 (§10) | `tffs.f:641` | `track6._ndivrad` (`track6.py:549`) |

## 付録 A. ヘルプの記述・K64 との違い (分割の一般論)

| 箇所 | 内容 |
|---|---|
| ヘルプ OCT・DECA・DODECA の transformation の $V$ | $V_n = \sum\bigl(-\tfrac\beta2 H_{,k}^2 + \gamma H_{,j}H_{,k}H_{,jk}\bigr)$ の形 (ヘルプの $H$ は多極のポテンシャルに当たる。符号と規格化は本章 (3.15) と対応するが、ヘルプの記号の定義は未確認)。**疑い**: 係数の $\beta$、$\gamma$ の前の $1/2$ の有無は、ヘルプの $\exp(:\cdot:)$ の符号・規格化に依るので、本章の導出 ($c_3 = \beta/2$、$c_5 = \gamma/2$。K64 のコードと記号的に一致) で読み替えた。ヘルプの式の単独の検証はしていない。 |
| ヘルプに無いもの | (1) 力勾配型の係数 $\alpha, \beta, \gamma$ の導出。(2) 補正キックが力の勾配の項を消す理由。(3) 分割の次数の議論 |

K64 の側の注意 (誤りではない): CALC の `tthine` は 2 次、追跡の `tthin` は 4 次で、同じ薄い多極でも分割が違う。したがって、軌道が大きくずれたとき、CALC の行列と追跡のヤコビアンは 2 次の分割の誤差だけ違う (07 章)。
