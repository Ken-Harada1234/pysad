# 03. シンプレクティックな分割の一般論

**この章で分かること**: K64 は、四極の本体も薄い多極も、「厳密に解ける 2 つの流れ」を交互に掛けて、ハミルトニアン全体の流れを近似する。この章は、その分割の一般論 — (1) ハミルトニアンを 2 つに分ける考え方、(2) 2 次の対称分割 (Strang) の誤差の導出、(3) 4 次への合成 (Yoshida)、(4) シンプレクティック性と長時間の振る舞い、(5) 薄い多極が使う**力勾配型**の 4 次の分割 (係数 $a, b, \beta, \gamma$ を導く)、(6) 刻みの数の決め方 — を、要素に依らない形で導く。四極と薄い多極の具体的な式は [04_quad.md](04_quad.md)、[07_thin_multipole.md](07_thin_multipole.md) に書き、この章はそれらが共通に使う道具を与える。式番号 $(3.n)$ から、§10 の Fortran の対応表と §11 の pysad の対応表を引ける。

**前提**: [C_poisson_lie.md](C_poisson_lie.md) (ポアソン括弧、$\exp(:f:)$、合成の順序、BCH 公式 (C.14) と 3 次までの導出)。5 次の項は、この章の §4.1 で、自由リー代数の基底とヤコビ恒等式から手で追える形に導く。ほかに [04_quad.md](04_quad.md) §6 (四極の Strang と Yoshida の結果) を読んでおくと、一般論と具体例の対応が分かる。

**読む順**: §1 (結論) → §3〜§5 (2 次・4 次) → §6〜§7 (性質) → §8 (力勾配型。薄い多極に必要) → §9 (刻みの数)。四極の章の図 ([図](figures/quad_integrator_order.png)) と数値 ([quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt)) は、この章の一般論の 1 つの例。この章の数値確認は、それとは別の簡単な系 (1 自由度) と、記号計算で行う。

## 1. 結論

- **分割の前提**: $H = H_a + H_b$ で、$H_a$ の流れ $A$ も $H_b$ の流れ $B$ も、閉じた式で厳密に解ける (K64 では、四極の本体なら $B$ = $H_2$ (2 次式) の線形な流れ、$A$ = $\Delta H$ (運動量だけ) の運動学的ドリフト。薄い多極なら $A$ = ドリフト、$B$ = 座標だけのポテンシャルのキック)。$A$ と $B$ は一般に**交換しない**ので、$e^{h(A+B)}$ は $e^{hA}e^{hB}$ に等しくなく、積の順序と長さの割り振りで近似の次数が決まる。
- **対称分割 (Strang)** $S_2(h) = e^{hA/2}\,e^{hB}\,e^{hA/2}$ は
$$
\log S_2(h) = h(A + B) + h^3\Bigl(-\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]\Bigr) + O(h^5)
$$
 (式 (3.2))。偶数次の項が出ず、誤差は 1 刻み $h^3$、全体 $h^2$ (2 次)。
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
| 4 次の合成の係数 | $w_1, w_0$ (式 (3.6)) | (K64 になし) | `_W1`、`_W0` (`quadbody.py`) |
| 力勾配型のドリフトの長さ | $a\,h$、$t\,h$ ($t = \alpha$) | `ala = al*alpha1`、`alb = al*alpha` | `ala`、`alb` (`_TT_ALPHA1`、`_TT_ALPHA`) |
| 補正キックの強さ | $c_3 = \beta/2$、$c_5 = \gamma/2$ | `beta`、`gamma` (`tquad.f:131-135` の `parameter`) | `_TT_BETA`、`_TT_GAMMA` |
| 精度パラメータ | $\varepsilon$ | `eps`・`eps0`・`epsdef` | `eps0` |
| 刻みの数 | $n$ | `ndiv` | `quadbody.ndiv()` |

**$\alpha, \beta, \gamma$ の注釈** (K64 `tquad.f:131-135`): `alpha` $= 1/\sqrt{12} = 0.288675\ldots$ (中央のドリフト 2 本の長さの比。数学の標準の記号ではない)、`alpha1` $= \tfrac12 - $ `alpha` $= 0.211325\ldots$ (両端のドリフトの比)、`beta` $= \tfrac16 - \alpha/2 = 0.0223291\ldots$ (補正キックの 3 次の強さ)、`gamma` $= \tfrac1{40} - \tfrac{1}{24\sqrt3} = 9.43739 \times 10^{-4}$ (補正キックの 5 次の強さ)。ソースのコメントは「`alpha=1/sqrt(12), beta=1/6-alpha/2, gamma=1/40-1/24/sqrt(3)`」。それぞれが何を消す係数かは §8 で導く。

**§10 の表に出る名前** (確認できた範囲。定義の場所はソース):

| 名前 | 意味 | 定義の場所 | 確認 |
|---|---|---|---|
| `aln0` | 四極の追跡の 1 刻みの長さ $= $ `al/ndiv` | `tsolqu.f:49` | 確認 |
| `ndiv` (四極) | 刻みの数 $= 1 + \lfloor L\lVert(K1, B_z)\rVert/\varepsilon\rfloor$ (`ndiv=1+int(abs(al*hypot(ak,bz)/eps))`。`ak` が積分強度 K1 なので、式の $K1$ は積分強度) | `tsolqu.f:47` | 確認 |
| `ak`、`bz` | 四極の積分強度 K1、ソレノイドの $B_z$ (`tsolqu` の引数) | `tsolqu.f:13` (宣言) | 確認 |
| `ak0x`、`ak0y` | ステアリング (偏向の積分強度) | `tsolqu.f:13` (宣言)。意味は §10 の表の注 | 宣言まで確認 |
| `eps`、`eps0`、`epsdef` | `eps0` = 台本の EPS、`epsdef = 0.2d0`。`eps0 == 0` なら `eps = epsdef`、そうでなければ `eps = epsdef*eps0` | `tsolqu.f:14`、`:41-46` | 確認 |
| `ndivmax` | 刻みの数の上限。`tbend.f:549` では 1024、`tbendi.f:16` では 1000 (ファイルごとに別の定数。ほかのファイルは未確認) | `tbend.f:549`、`tbendi.f:16` | 一部確認 |
| `phic` | 偏向磁石 (CALC) の $\phi_c$ = `phi0*rbc` (`phi0` は偏向角。`rbc` の定義は未確認) | `tbende.f:454` | 式まで確認 |
| `epsbend` | `eps = merge(epsbend, epsbend*eps0, eps0 == 0.d0)` に出る量。`epsbend` の値の定義は未確認 | `tbende.f:455` | 未確認 |
| `arad` | 放射の分割 `ndivrad` の係数 0.02 (`real*8, parameter :: arad=0.02d0`) | `tffs.f:644` | 確認 |
| `trpt` | `ndivrad` に出る真偽値 (`fff%trpt` を指す。意味は未確認) | `tffs.f:589` | 未確認 |
| `phir`、`f1r`、`f2r` | 放射つき偏向磁石の追跡で使う偏向角と端の長さ (本文 §10 の式のとおり。定義は未確認) | `tbend.f` | 未確認 |
| `ndiv0` | MULT の要素のパラメータ (`tmulti` の引数。`ndiv=ndiv0` が `tmulti.f:107`) | `tmulti.f:62`、`:107` | 一部確認 |
| `kord` | 薄い多極の次数の指標 `kord = nord/2 - 1` (`nord` = 4 なら 1、6 なら 2 (6 極)、8 なら 3 (8 極))。補正キックの式の場合分け (§14 の `f1`〜`f5`) に使う | `tquad.f:170` (`tthin`) | 確認 |
| `akf` | `akf = ak/fact(kord)` ($= \mathrm{ak}/k!$)。長さのある多極 (`al /= 0`) では、その直後に `akf = akf*0.5d0` で $\mathrm{ak}/k!/2$ になる (§14 の `aki` の定義と同じ) | `tquad.f:171`、`:182` | 確認 |
| `setndivelm(ia, ndiv)` | 要素番号 `ia` の刻みの数 `ndiv` を、要素のパラメータの記録に登録するサブルーチン。`tthin` は `call setndivelm(l_track,1)` (刻み 1) | `tffs.f:1537` (定義)、`tquad.f:150` | 定義の冒頭まで確認。記録の使い道は未確認 |
| `l_track` | `setndivelm` の第 1 引数に渡される変数。追跡中の要素の番号と読めるが、宣言の場所は未確認 | `tquad.f:27、70、82、150` (使用) | 未確認 |
| `epsrad` | 放射による分割数の精度の定数 $10^{-6}$ (`real*8, parameter :: epsrad=1.d-6`)。偏向磁石の CALC では `epsr1 = merge(epsrad, epsrad*eps0, eps0 == 0.d0)` として `nrad` (放射の分割数) の式に入る | `toplvl.f:291`、`tbende.f:462`、`tsolque.f:621` | 確認 |

## 3. ハミルトニアンの分割と、リー演算子の書き方

流れを「演算子の指数」で書く。[C_poisson_lie.md](C_poisson_lie.md) の規約で、ハミルトニアン $H$ (位置 $s$ によらない) の長さ $h$ の転送は $\exp(-h:H:)$。$A \equiv -:H_a:$、$B \equiv -:H_b:$ (演算子。$A, B$ は関数への作用で、$:f:\,g = [f, g]$)。$e^{hA}$ は $H_a$ の長さ $h$ の流れ、$e^{hB}$ は $H_b$ の流れ。積 $e^{hA}e^{hB}$ は「$e^{hA}$ が先、$e^{hB}$ が後」に粒子が通る順 (C.13 の規約)。以下で主に使う積 (3.1)・(3.3)・(3.9) は前後対称な並びなので、積の順序の規約を逆にしても同じ式になる (非対称な積が出るのは、§4 の途中の $e^{A/2}e^{B}$ だけで、そこでは $\log$ の係数をそのまま BCH 公式 (C.14) から読む)。

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

四極の本体で $A = D$ (運動学的ドリフト)、$B = H_2$ の流れの場合は、$[B,[B,A]]$ は (4.19) の $C_3$ の $[\mathcal B,[\mathcal B,\mathcal D]]$ の項 (係数 $1/12$) に、$[A,[A,B]]$ は $[\mathcal D,[\mathcal D,\mathcal B]]$ の項 (係数 $-1/24$) に当たる。

### 4.1 5 次の項 (補題 S): 自由リー代数の基底とヤコビ恒等式

(3.2) の次の項 (5 次) を導く。§5 の Yoshida の 5 次の残り (3.7) と、§8.3 の力勾配型の 5 次 (3.13) が、この結果を使う。数値は後でまとめて照合するが、ここでは**係数を 1 つずつ手で導く**。確認は [integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (スクリプト `integ_04_bch5_steps.py`。自由代数での直接計算と、各項の値)。

**(a) 次数 5 の基底と、ヤコビ恒等式による簡約**。$X$、$Y$ を非可換な生成元とし、括弧の名前を決める:

$$
U \equiv [X, Y],\qquad V \equiv [X, U] = [X,[X,Y]],\qquad W \equiv [Y, U] = [Y,[X,Y]]
\tag{3.2a}
$$

3 次の項は $V$ と $W$ ((3.2) の $[A,[A,B]] = V$、$[B,[B,A]] = -W$)。4 次の括弧は $[X,V]$、$[Y,V]$、$[X,W]$、$[Y,W]$ だが、ヤコビ恒等式 $[X,[Y,U]] = [[X,Y],U] + [Y,[X,U]]$ ($[U,U] = 0$ を使うと) から

- (J1) $[X, W] = [X,[Y,U]] = [U, U] + [Y,[X,U]] = [Y, V]$

なので、独立なのは $[X,V]$、$[Y,V]$、$[Y,W]$ の 3 つ。5 次 (生成元 5 個の括弧) の独立な元の数は、自由リー代数の次元の公式 (Witt の公式。2 生成元、次数 5 で $\tfrac15(2^5 - 2) = 6$) から 6 個。次の 6 個を基底に選ぶ ($X$ と $Y$ の個数ごとに並べた):

| 名前 | 定義 | ($X$ の数, $Y$ の数) |
|---|---|---|
| $c_1$ | $[X,[X,V]]$ | (4, 1) |
| $c_2$ | $[Y,[X,V]]$ | (3, 2) |
| $c_3$ | $[U, V]$ | (3, 2) |
| $c_4$ | $[X,[Y,W]]$ | (2, 3) |
| $c_5$ | $[U, W]$ | (2, 3) |
| $c_6$ | $[Y,[Y,W]]$ | (1, 4) |

(2 つの生成元の入れ替えで $c_1 \leftrightarrow c_6$ などが対応する。$c_1 \ldots c_6$ が 1 次独立で、生成元 5 個の右にネストした括弧 $[x_1,[x_2,[x_3,[x_4,x_5]]]]$ ($x_i \in \{X, Y\}$、32 通り) がすべてこの 6 個の 1 次結合になることは、自由代数で確認した ([integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (a))。) 基底以外の括弧は、次のヤコビ恒等式で基底に直す。以下の (J2)〜(J5) は、(J1) $[X,W] = [Y,V]$ と、$[Y,X] = -U$、$[U,\cdot]$ の定義だけから出る:

- (J2) $[X,[Y,V]] = [[X,Y],V] + [Y,[X,V]] = [U,V] + [Y,[X,V]] = c_3 + c_2$ (ヤコビ恒等式 $[X,[Y,V]] = [[X,Y],V] + [Y,[X,V]]$)。
- (J3) $[Y,[Y,V]] = [Y,[X,W]]$ ((J1) で $[Y,V] = [X,W]$) $= [[Y,X],W] + [X,[Y,W]] = -[U,W] + [X,[Y,W]] = c_4 - c_5$。
- (J4) $[X,[X,W]] = [X,[Y,V]]$ ((J1)) $= c_3 + c_2$ ((J2))。
- (J5) $[Y,[X,W]] = [Y,[Y,V]]$ ((J1)) $= c_4 - c_5$ ((J3))。

**(b) 補題 S**。前後対称な積の 5 次までの対数は

$$
\log\bigl(e^{X/2}e^{Y}e^{X/2}\bigr) = X + Y - \tfrac1{24}V - \tfrac1{12}W + \tfrac7{5760}c_1 + \tfrac7{1440}c_2 + \tfrac1{480}c_3 + \tfrac1{180}c_4 - \tfrac1{120}c_5 + \tfrac1{720}c_6 + O(7)
\tag{3.2b}
$$

($-\tfrac1{24}V - \tfrac1{12}W$ は (3.2) の $C_3$。偶数次の項は 0。$O(7)$ は 7 次以上。**$X$、$Y$ は任意のリー代数の元に置き換えてよい** (自由代数の恒等式だから)。§5、§8.3 は、$X$、$Y$ に和や括弧を入れて使う。) 以下、6 個の係数 ($c_1$ から $c_6$) を、$X$、$Y$ の個数ごとに導く。

**(c-0) 準備: 指数写像の微分**。$S(\varepsilon) = e^{Z(\varepsilon)}$ とする。$\frac{d}{d\varepsilon}e^{Z} = \int_0^1e^{\sigma Z}\,Z'\,e^{(1-\sigma)Z}d\sigma$ ($Z'$ は $\varepsilon$ 微分。べき級数の項ごとに確認できる標準の公式) の右から $e^{-Z}$ をかけると、$e^{\sigma Z}Z'e^{-\sigma Z} = e^{\sigma\,\mathrm{ad}_Z}Z'$ ($\mathrm{ad}_Z\,g \equiv [Z, g]$) なので

$$
\frac{dS}{d\varepsilon}\,S^{-1} = \int_0^1e^{\sigma\,\mathrm{ad}_Z}Z'\,d\sigma = \frac{e^{\mathrm{ad}_Z} - 1}{\mathrm{ad}_Z}\,Z'
\quad\Rightarrow\quad
Z' = F(\mathrm{ad}_Z)\left(\frac{dS}{d\varepsilon}S^{-1}\right),\qquad F(u) \equiv \frac{u}{e^u - 1} = 1 - \frac u2 + \frac{u^2}{12} - \frac{u^4}{720} + \cdots
\tag{3.2c}
$$

($F$ の係数 $f_0 = 1$、$f_1 = -\tfrac12$、$f_2 = \tfrac1{12}$、$f_3 = 0$、$f_4 = -\tfrac1{720}$ はベルヌーイ数 ($B_n/n!$)。$F(u) + u/2 = \tfrac u2\coth\tfrac u2$ が偶関数なので $f_3 = 0$。)

**(c-1) $(X, Y) = (4, 1)$ ($Y$ が 1 個)**。$S(\varepsilon) = e^{X/2}e^{\varepsilon Y}e^{X/2}$、$Z(0) = X$。$dS/d\varepsilon = e^{X/2}Ye^{\varepsilon Y}e^{X/2}$ に右から $S^{-1} = e^{-X/2}e^{-\varepsilon Y}e^{-X/2}$ をかけると $e^{X/2}Ye^{-X/2} = e^{\mathrm{ad}_X/2}Y$ ($\varepsilon$ によらない)。$\varepsilon = 0$ の (3.2c) から、$Z = X + \varepsilon Z_1 + \cdots$ の $Z_1$ は

$$
Z_1 = F(\mathrm{ad}_X)\,e^{\mathrm{ad}_X/2}Y = \frac{\mathrm{ad}_X/2}{\sinh(\mathrm{ad}_X/2)}\,Y = \Bigl(1 - \tfrac1{24}\mathrm{ad}_X^2 + \tfrac7{5760}\mathrm{ad}_X^4 - \cdots\Bigr)Y
$$

($F(u)e^{u/2} = u/(e^{u/2} - e^{-u/2}) = (u/2)/\sinh(u/2)$。$\sinh y = y(1 + y^2/6 + y^4/120 + \cdots)$ から $y/\sinh y = 1 - y^2/6 + 7y^4/360$、$y = u/2$ とおく。) したがって $-\tfrac1{24}[X,[X,Y]] = -\tfrac1{24}V$ ((3.2) の係数に一致) と、5 次の $\tfrac7{5760}\mathrm{ad}_X^4Y = \tfrac7{5760}[X,[X,[X,[X,Y]]]] = \tfrac7{5760}c_1$ ($[X,[X,[X,[X,Y]]]] = [X,[X,V]] = c_1$)。

**(c-2) $(1, 4)$ ($X$ が 1 個)**。$S(\varepsilon) = e^{\varepsilon X/2}e^{Y}e^{\varepsilon X/2}$、$Z(0) = Y$。$dS/d\varepsilon = \tfrac X2\,S + e^{\varepsilon X/2}e^{Y}\tfrac X2\,e^{\varepsilon X/2}$ (積の微分) に右から $S^{-1} = e^{-\varepsilon X/2}e^{-Y}e^{-\varepsilon X/2}$ をかけると、第 1 項は $\tfrac X2$、第 2 項は $e^{\varepsilon X/2}e^{Y}\tfrac X2\,e^{-Y}e^{-\varepsilon X/2} = e^{\varepsilon\,\mathrm{ad}_X/2}e^{\mathrm{ad}_Y}\tfrac X2$ ($e^{A}Be^{-A} = e^{\mathrm{ad}_A}B$):

$$
\frac{dS}{d\varepsilon}S^{-1} = \frac X2 + e^{\varepsilon\,\mathrm{ad}_X/2}e^{\mathrm{ad}_Y}\,\frac X2
\tag{3.2d}
$$

$\varepsilon = 0$ で $h_0 \equiv \tfrac12(1 + e^{\mathrm{ad}_Y})X$。$Z_1 = F(\mathrm{ad}_Y)h_0 = \dfrac{u}{e^u - 1}\cdot\dfrac{1 + e^u}2X = \dfrac u2\coth\dfrac u2\,X$ ($u = \mathrm{ad}_Y$) $= \bigl(1 + \tfrac{u^2}{12} - \tfrac{u^4}{720} + \cdots\bigr)X$。3 次: $\tfrac1{12}[Y,[Y,X]] = -\tfrac1{12}W$ ((3.2) の係数に一致)。5 次: $-\tfrac1{720}\mathrm{ad}_Y^4X$。$\mathrm{ad}_Y^4X = [Y,[Y,[Y,[Y,X]]]] = -[Y,[Y,[Y,U]]] = -[Y,[Y,W]] = -c_6$ ($[Y,X] = -U$) なので、$+\tfrac1{720}c_6$。

**(c-3) $(3, 2)$ ($Y$ が 2 個)**。(c-1) の続き ($\varepsilon^2$ の項)。(3.2c) を $\varepsilon$ で微分して $\varepsilon = 0$ とすると ($g \equiv e^{\mathrm{ad}_X/2}Y$ は $\varepsilon$ によらない)、$2Z_2 = \bigl[\tfrac{d}{d\varepsilon}F(\mathrm{ad}_Z)\bigr]g$。$F(\mathrm{ad}_Z) = \sum_nf_n\,\mathrm{ad}_Z^n$、$\mathrm{ad}_Z = \mathrm{ad}_X + \varepsilon\,\mathrm{ad}_{Z_1}$ なので ($\mathrm{ad}_X$ と $\mathrm{ad}_{Z_1}$ は可換でないから、$\mathrm{ad}_{Z_1}$ を入れる位置で $n$ 通りの項になる)

$$
2Z_2 = \sum_nf_n\sum_{i+j=n-1}\mathrm{ad}_X^i\,\mathrm{ad}_{Z_1}\,\mathrm{ad}_X^j\,g,\qquad
Z_1 = Y - \tfrac1{24}\mathrm{ad}_X^2Y + \cdots,\quad g = Y + \tfrac12\mathrm{ad}_XY + \tfrac18\mathrm{ad}_X^2Y + \tfrac1{48}\mathrm{ad}_X^3Y + \cdots
$$

($g = \sum_k\mathrm{ad}_X^kY/(2^kk!)$)。$X$ が 3 個、$Y$ が 2 個の項を拾う ($\mathrm{ad}_X$ の総数 $i + j$、$Z_1$ の中の $\mathrm{ad}_X$ の数 $m$、$g$ の中の数 $k$ で $i + j + m + k = 3$)。$\mathrm{ad}_XY = U$、$\mathrm{ad}_X^2Y = V$、$\mathrm{ad}_X^3Y = [X,V]$ を使い、(J2)、(J4) で $c_2$、$c_3$ に直す。各項 (係数 $f_n$ を含めた値):

| $n$ | 項 | 括弧の値 | $c_2$、$c_3$ での値 |
|---|---|---|---|
| 1 ($f_1 = -\tfrac12$) | $[Y,\ \tfrac1{48}[X,V]]$ ($m = 0$、$k = 3$) | $\tfrac1{48}c_2$ | $-\tfrac1{96}c_2$ |
| 1 | $[-\tfrac1{24}V,\ \tfrac12U]$ ($m = 2$、$k = 1$) | $-\tfrac1{48}[V,U] = \tfrac1{48}c_3$ | $-\tfrac1{96}c_3$ |
| 2 ($f_2 = \tfrac1{12}$) | $[Y,\ \mathrm{ad}_X(\tfrac18V)]$ ($m = 0$、$k = 2$、$j = 1$) | $\tfrac18[Y,[X,V]] = \tfrac18c_2$ | $\tfrac1{96}c_2$ |
| 2 | $\mathrm{ad}_X[Y,\ \tfrac18V]$ ($m = 0$、$k = 2$、$i = 1$) | $\tfrac18[X,[Y,V]] = \tfrac18(c_2 + c_3)$ ((J2)) | $\tfrac1{96}(c_2 + c_3)$ |
| 2 | $[-\tfrac1{24}V,\ \mathrm{ad}_XY]$ ($m = 2$、$k = 0$、$j = 1$) | $-\tfrac1{24}[V,U] = \tfrac1{24}c_3$ | $\tfrac1{288}c_3$ |
| 2 | $\mathrm{ad}_X[-\tfrac1{24}V,\ Y]$ ($m = 2$、$k = 0$、$i = 1$) | $\tfrac1{24}[X,[Y,V]] = \tfrac1{24}(c_2 + c_3)$ | $\tfrac1{288}(c_2 + c_3)$ |
| 4 ($f_4 = -\tfrac1{720}$) | $[Y,\ \mathrm{ad}_X^3Y]$ ($i = 0$、$j = 3$) | $[Y,[X,V]] = c_2$ | $-\tfrac1{720}c_2$ |
| 4 | $\mathrm{ad}_X[Y,\ \mathrm{ad}_X^2Y]$ ($i = 1$、$j = 2$) | $[X,[Y,V]] = c_2 + c_3$ | $-\tfrac1{720}(c_2 + c_3)$ |
| 4 | $\mathrm{ad}_X^2[Y,\ \mathrm{ad}_XY]$ ($i = 2$、$j = 1$) | $[X,[X,[Y,U]]] = [X,[X,W]] = c_2 + c_3$ ((J4)) | $-\tfrac1{720}(c_2 + c_3)$ |

($n = 3$ は $f_3 = 0$。$[Y,Y] = 0$ なので $j = 0$ の項は 0。) $c_2$ の係数の和: $-\tfrac1{96} + \tfrac1{96} + \tfrac1{96} + \tfrac1{288} - \tfrac3{720} = \tfrac1{72} - \tfrac1{240} = \tfrac7{720}$。$c_3$ の係数の和: $-\tfrac1{96} + \tfrac1{96} + \tfrac2{288} - \tfrac2{720} = \tfrac1{144} - \tfrac1{360} = \tfrac1{240}$。$Z_2$ はこの半分で、**$\tfrac7{1440}c_2 + \tfrac1{480}c_3$**。

**(c-4) $(2, 3)$ ($X$ が 2 個)**。(c-2) の続き ($\varepsilon^2$ の項)。(3.2d) を $\varepsilon$ で微分して $\varepsilon = 0$: $h_1 = \tfrac14\mathrm{ad}_Xe^{\mathrm{ad}_Y}X = \tfrac14\sum_k\tfrac1{k!}[X,\mathrm{ad}_Y^kX]$。(3.2c) を微分して

$$
2Z_2 = \sum_nf_n\sum_{i+j=n-1}\mathrm{ad}_Y^i\,\mathrm{ad}_{Z_1}\,\mathrm{ad}_Y^j\,h_0 + F(\mathrm{ad}_Y)\,h_1,\qquad
Z_1 = X + \tfrac1{12}\mathrm{ad}_Y^2X - \cdots,\quad h_0 = X + \tfrac12\mathrm{ad}_YX + \tfrac14\mathrm{ad}_Y^2X + \tfrac1{12}\mathrm{ad}_Y^3X + \cdots
$$

$Y$ が 3 個の項を拾う。$\mathrm{ad}_YX = -U$、$\mathrm{ad}_Y^2X = -W$、$\mathrm{ad}_Y^3X = -[Y,W]$、$[X,\mathrm{ad}_Y^3X] = -c_4$ などと、(J3)、(J5) を使う。第 1 部 (各項は係数 $f_n$ 込み):

| $n$ | 項 | 括弧の値 | 値 |
|---|---|---|---|
| 1 ($-\tfrac12$) | $[X,\ \tfrac1{12}\mathrm{ad}_Y^3X]$ | $-\tfrac1{12}c_4$ | $+\tfrac1{24}c_4$ |
| 1 | $[\tfrac1{12}\mathrm{ad}_Y^2X,\ \tfrac12\mathrm{ad}_YX]$ | $\tfrac1{24}[-W,-U] = \tfrac1{24}[W,U] = -\tfrac1{24}c_5$ | $+\tfrac1{48}c_5$ |
| 2 ($\tfrac1{12}$) | $[X,\ \mathrm{ad}_Y(\tfrac14\mathrm{ad}_Y^2X)]$ ($j = 1$) | $-\tfrac14c_4$ | $-\tfrac1{48}c_4$ |
| 2 | $\mathrm{ad}_Y[X,\ \tfrac14\mathrm{ad}_Y^2X]$ ($i = 1$) | $-\tfrac14[Y,[X,W]] = -\tfrac14(c_4 - c_5)$ ((J5)) | $-\tfrac1{48}c_4 + \tfrac1{48}c_5$ |
| 2 | $[\tfrac1{12}\mathrm{ad}_Y^2X,\ \mathrm{ad}_YX]$ ($j = 1$) | $[-\tfrac1{12}W,\ -U] = \tfrac1{12}[W,U] = -\tfrac1{12}c_5$ | $-\tfrac1{144}c_5$ |
| 2 | $\mathrm{ad}_Y[\tfrac1{12}\mathrm{ad}_Y^2X,\ X]$ ($i = 1$) | $\tfrac1{12}[Y,[X,W]] = \tfrac1{12}(c_4 - c_5)$ | $\tfrac1{144}(c_4 - c_5)$ |
| 4 ($-\tfrac1{720}$) | $[X,\ \mathrm{ad}_Y^3X]$ ($i = 0$、$j = 3$) | $-c_4$ | $+\tfrac1{720}c_4$ |
| 4 | $\mathrm{ad}_Y[X,\ \mathrm{ad}_Y^2X]$ ($i = 1$、$j = 2$) | $[Y,[X,-W]] = -(c_4 - c_5)$ | $\tfrac1{720}(c_4 - c_5)$ |
| 4 | $\mathrm{ad}_Y^2[X,\ \mathrm{ad}_YX]$ ($i = 2$、$j = 1$) | $[Y,[Y,[X,-U]]] = -[Y,[Y,V]] = -(c_4 - c_5)$ ((J3)) | $\tfrac1{720}(c_4 - c_5)$ |

第 2 部 ($F(\mathrm{ad}_Y)h_1$。$Y$ の総数 3):

| $n$ | 項 | 括弧の値 | 値 |
|---|---|---|---|
| 0 ($1$) | $\tfrac14\cdot\tfrac1{3!}[X,\ \mathrm{ad}_Y^3X]$ | $-\tfrac1{24}c_4$ | $-\tfrac1{24}c_4$ |
| 1 ($-\tfrac12$) | $\mathrm{ad}_Y\bigl(\tfrac14\cdot\tfrac1{2!}[X,\ \mathrm{ad}_Y^2X]\bigr)$ | $-\tfrac18[Y,[X,W]] = -\tfrac18(c_4 - c_5)$ | $+\tfrac1{16}(c_4 - c_5)$ |
| 2 ($\tfrac1{12}$) | $\mathrm{ad}_Y^2\bigl(\tfrac14[X,\ \mathrm{ad}_YX]\bigr)$ | $-\tfrac14[Y,[Y,V]] = -\tfrac14(c_4 - c_5)$ | $-\tfrac1{48}(c_4 - c_5)$ |

$c_4$ の係数: 第 1 部 $\tfrac1{24} - \tfrac2{48} + \tfrac1{144} + \tfrac3{720} = \tfrac1{90}$、第 2 部 $-\tfrac1{24} + \tfrac1{16} - \tfrac1{48} = 0$、合計 $\tfrac1{90}$。$c_5$ の係数: 第 1 部 $\tfrac1{48} + \tfrac1{48} - \tfrac1{144} - \tfrac1{144} - \tfrac2{720} = \tfrac1{40}$、第 2 部 $-\tfrac1{16} + \tfrac1{48} = -\tfrac1{24}$、合計 $-\tfrac1{60}$。$Z_2$ はこの半分で、**$\tfrac1{180}c_4 - \tfrac1{120}c_5$**。

(c-1)〜(c-4) で、(3.2b) の 6 個の係数がそろった。自由代数で $\log(e^{X/2}e^Ye^{X/2})$ を直接計算した 5 次の係数 (7/5760、7/1440、1/480、1/180、$-$1/120、1/720) と一致する ([integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (b)。(c-3)、(c-4) の表の各項の値も同出力に全部ある)。$X$、$Y$ に一般のリー代数の元 (例: $2A + 3[A,B]$ と $B$) を入れても (3.2b) が 5 次まで成り立つことも確認した。

**(d) 関数表現**。具体的な表現 ($A = -:T:$、$B = -:V:$、$T = p^2/2$。1 自由度。$[f,g] = f_xg_p - f_pg_x$) で、$X = A$、$Y = B$ とおいた括弧の関数 ($:f_{\ldots}:$ の $f$。リー演算子の交換子は関数の括弧に対応、$[:f:, :g:] = :[f,g]:$) を求める。$A \mapsto -p^2/2$、$B \mapsto -V(x)$:

- $U = [A,B] \mapsto [-p^2/2,\ -V] = [p^2/2,\ V] = -pV'$、$V_A \equiv [A,U] \mapsto [-p^2/2,\ -pV'] = -p^2V''$ ($V_A$ は括弧 $V$ ((3.2a)) の関数。ポテンシャル $V(x)$ との混同を避けて $V_A$ と書く)、$W = [B,U] \mapsto [-V,\ -pV'] = V'^2$。
- $[A,V_A] \mapsto [-p^2/2,\ -p^2V''] = -p^3V'''$、$c_1 = [A,[A,V_A]] \mapsto [-p^2/2,\ -p^3V'''] = -p^4V''''$。
- $c_2 = [B,[A,V_A]] \mapsto [-V,\ -p^3V'''] = 3p^2V'V'''$。$c_3 = [U,V_A] \mapsto [-pV',\ -p^2V''] = 2p^2V''^2 - p^2V'V'''$ ($f_xg_p - f_pg_x = (-pV'')(-2pV'') - (-V')(-p^2V''')$)。
- $[B,W] \mapsto [-V,\ V'^2] = 0$ (どちらも $x$ だけの関数)。したがって $c_4 = [A,[B,W]] \mapsto 0$、$c_6 = [B,[B,W]] \mapsto 0$。
- $c_5 = [U,W] \mapsto [-pV',\ V'^2] = 0 - (-V')(2V'V'') = 2V'^2V''$。

まとめ ($V_j$ ではなく 1 自由度の $V'$ など):

| 括弧 | 関数 |
|---|---|
| $U$、$V$、$W$ ((3.2a) で $X = A$、$Y = B$) | $-pV'$、$-p^2V''$、$V'^2$ |
| $c_1$、$c_2$、$c_3$ | $-p^4V''''$、$3p^2V'V'''$、$p^2(2V''^2 - V'V''')$ |
| $c_4$、$c_5$、$c_6$ | $0$、$2V'^2V''$、$0$ |

(確認: [integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (c)。) Strang の $C_5$ ($X = A$、$Y = B$) の関数は、(3.2b) の係数を掛けて

$$
C_5 \mapsto F_5 \equiv -\tfrac1{60}V'^2V'' + \Bigl(\tfrac1{80}V'V''' + \tfrac1{240}V''^2\Bigr)p^2 - \tfrac7{5760}V''''p^4
\tag{3.2e}
$$

($p^2V'V'''$: $\tfrac7{1440}\cdot3 - \tfrac1{480} = \tfrac{21 - 3}{1440} = \tfrac1{80}$。$p^2V''^2$: $\tfrac2{480} = \tfrac1{240}$。$V'^2V''$: $-\tfrac1{120}\cdot2 = -\tfrac1{60}$。$p^4$: $\tfrac7{5760}\cdot(-1)$。)

## 5. 4 次への合成 (Yoshida)

対称な 2 次の方法 $S_2$ を 3 つ、長さ $w_1h$、$w_0h$、$w_1h$ で並べる:

$$
S_4(h) = S_2(w_1h)\,S_2(w_0h)\,S_2(w_1h)
\tag{3.3}
$$

$S_2(w_ih)$ は (3.1) の $A$、$B$ を $w_iA$、$w_iB$ に置き換えたものなので、(3.2)・(3.2b) から ($Z_1 \equiv A + B$)

$$
\Omega_i \equiv \log S_2(w_ih) = w_ih\,Z_1 + (w_ih)^3C_3 + (w_ih)^5C_5 + O(h^7)
$$

(各項は $h$ と $w_i$ について、その次数のべき。偶数次はない。$C_3 = -\tfrac1{24}[A,[A,B]] + \tfrac1{12}[B,[B,A]]$、$C_5$ は (3.2b) の 5 次の項で $X = A$、$Y = B$ としたもの。)

$S_4 = e^{\Omega_1}e^{\Omega_0}e^{\Omega_1}$ は前後対称な 3 つの積で、補題 S (3.2b) の形 $e^{X'/2}e^{Y'}e^{X'/2}$ に、**$X' = 2\Omega_1$、$Y' = \Omega_0$ としてそのまま当てはまる** (BCH を 2 回重ねる式を別に作る必要はない)。(3.2b) から

$$
\log S_4 = X' + Y' + C_3(X', Y') + C_5(X', Y') + O(h^7),\qquad
C_3(X', Y') = -\tfrac1{24}[X',[X',Y']] + \tfrac1{12}[Y',[Y',X']]
\tag{3.3a}
$$

(補題 S の $C_5$ を $C_5(X', Y')$ と書いた。$X'$、$Y'$ について 5 次の括弧の和。) 項ごとに次数を数える。

- **線形の部分**: $X' + Y' = (2w_1 + w_0)\,hZ_1 + (2w_1^3 + w_0^3)\,h^3C_3 + (2w_1^5 + w_0^5)\,h^5C_5 + O(h^7)$。
- **2 次の交換子 $[X', Y']$ は、(3.3a) に現れない** (前後対称な積だから。$e^{\Omega_1}e^{\Omega_0}$ を先に作って BCH を 2 回重ねると、$\tfrac12([\Omega_1,\Omega_0] + [\Omega_1,\Omega_1] + [\Omega_0,\Omega_1])$ のような 2 次の項が途中に出るが、$[\Omega_1,\Omega_1] = 0$ と、$[\Omega_1,\Omega_0] + [\Omega_0,\Omega_1] = 0$ (反対称) で相殺する。補題 S は、その相殺を済ませた形)。
- **3 重の括弧 $C_3(X', Y')$**: $X' = a_1hZ_1 + a_3h^3C_3 + \cdots$、$Y' = b_1hZ_1 + b_3h^3C_3 + \cdots$ ($a_1 = 2w_1$、$a_3 = 2w_1^3$、$b_1 = w_0$、$b_3 = w_0^3$)。$[Z_1, Z_1] = 0$ なので、$X$ だけの入れ子は消え、$[X', Y'] = (a_1b_3 - a_3b_1)\,h^4[Z_1, C_3] + O(h^6)$。したがって、$[X',[X',Y']] = a_1\Delta\,h^5[Z_1,[Z_1,C_3]] + O(h^7)$、$[Y',[Y',X']] = -b_1\Delta\,h^5[Z_1,[Z_1,C_3]] + O(h^7)$ ($\Delta \equiv a_1b_3 - a_3b_1$。$[Y',X'] = -[X',Y']$)。$h^3$ 以下に交換子の項は出ず、$h^5$ に
$$
C_3(X',Y')\Big|_{h^5} = \rho\,h^5\,[Z_1,[Z_1,C_3]],\qquad \rho \equiv -\Delta\Bigl(\frac{a_1}{24} + \frac{b_1}{12}\Bigr)
$$
 が出る ($h \cdot h \cdot h^3 = h^5$)。
- **5 重の括弧 $C_5(X', Y')$**: $X'$ と $Y'$ を 5 個含み、$h^5$ の項は $Z_1$ だけの括弧 ($=0$)。残りは $h^7$ 以上。

よって 3 次までと、5 次の形は

$$
\log S_4 = (2w_1 + w_0)\,h(A+B) + (2w_1^3 + w_0^3)\,h^3 C_3 + h^5\Bigl[(2w_1^5 + w_0^5)\,C_5 + \rho\,[Z_1,[Z_1,C_3]]\Bigr] + O(h^7)
\tag{3.4}
$$

(上のとおり、偶数次の項は相殺する。) 正しい $h(A+B)$ の係数 1 と、$h^3$ の消去から

$$
2w_1 + w_0 = 1,\qquad 2w_1^3 + w_0^3 = 0
\tag{3.5}
$$

$w_0 = 1 - 2w_1$ を入れると 2 番目の式は $(1 - 2w_1)^3 = -2w_1^3$。実数の範囲で両辺の立方根を取ると $1 - 2w_1 = -2^{1/3}w_1$、これを $w_1$ について解く: $1 = (2 - 2^{1/3})\,w_1$。数値は $2^{1/3} = 1.2599210$、$2 - 2^{1/3} = 0.7400790$、$w_1 = 1.3512072$、$w_0 = 1 - 2\times1.3512072 = -1.7024144$ ([integ_03_expansions.out.txt](outputs/integ_03_expansions.out.txt) (a)。$2w_1^3 + w_0^3 = 2\times10^{-122}$)。

$$
w_1 = \frac{1}{2 - 2^{1/3}} = 1.3512071919\ldots,\qquad w_0 = 1 - 2w_1 = -1.7024143839\ldots
\tag{3.6}
$$

(四極の章の (4.20)〜(4.22) と同じ。)

**一般化 (高次への飛び越し)**: 対称な $p$ 次 ($p$ は偶数) の方法 $S$ は、$\log S(h) = hX + c\,h^{p+1} + O(h^{p+3})$ ($c$ は演算子の定数。偶数次は出ない)。(3.4) と同じ議論 (2 次の交換子は補題 S の対称形で最初から無く、3 重以上は $h^{p+2}$ 以上) で、$S(w_1h)S(w_0h)S(w_1h)$ の $\log$ は $(2w_1 + w_0)hX + (2w_1^{p+1} + w_0^{p+1})\,c\,h^{p+1} + \cdots$。$2w_1 + w_0 = 1$ と $2w_1^{p+1} + w_0^{p+1} = 0$ から、(3.5) と同じ手順で $(1 - 2w_1)^{p+1} = -2w_1^{p+1}$、$p + 1$ は奇数なので実数の $(p+1)$ 乗根が取れて $1 - 2w_1 = -2^{1/(p+1)}w_1$、したがって $w_1 = 1/(2 - 2^{1/(p+1)})$。$p = 2$ が (3.6)、$p = 4$ なら $1/(2 - 2^{1/5}) = 1.1746718$ で 6 次の方法ができる (数値は同 (a))。K64・pysad は 4 次までしか使わない。

**5 次の残り (3.7) の導出**。(3.4) の $h^5$ の項を、(3.5) を使って整理する。

1. $\Delta = a_1b_3 - a_3b_1 = 2w_1w_0^3 - 2w_1^3w_0$。(3.5) の $w_0^3 = -2w_1^3$ を第 1 項に入れ、$\Delta = -4w_1^4 - 2w_1^3w_0 = -2w_1^3(2w_1 + w_0) = -2w_1^3$ ($2w_1 + w_0 = 1$)。
2. $\dfrac{a_1}{24} + \dfrac{b_1}{12} = \dfrac{2w_1}{24} + \dfrac{w_0}{12} = \dfrac{w_1 + w_0}{12}$。
3. $\rho = -\Delta\cdot\dfrac{w_1 + w_0}{12} = \dfrac{w_1^3(w_1 + w_0)}{6} = \dfrac{w_1^3(1 - w_1)}{6} = -0.1444036$。$s_5 \equiv 2w_1^5 + w_0^5 = (2 - 2^{5/3})\,w_1^5 = -5.2914471$ ($w_0 = -2^{1/3}w_1$ から $w_0^5 = -2^{5/3}w_1^5$)。
4. 具体的な表現 ((d)) での関数: $C_5 \mapsto F_5$ ((3.2e))。$C_3 = -\tfrac1{24}V - \tfrac1{12}W \mapsto \tfrac1{24}p^2V'' - \tfrac1{12}V'^2$ ($V \mapsto -p^2V''$、$W \mapsto V'^2$)、$Z_1 = A + B \mapsto -p^2/2 - V$ (ここの $V$ は括弧ではなくポテンシャル)。$[Z_1, C_3]$ の関数 $G_1 = f_xg_p - f_pg_x$ ($f = -p^2/2 - V$、$g = \tfrac1{24}p^2V'' - \tfrac1{12}V'^2$): $f_x = -V'$、$f_p = -p$、$g_p = \tfrac1{12}pV''$、$g_x = \tfrac1{24}p^2V''' - \tfrac16V'V''$ なので
 $G_1 = -\tfrac1{12}pV'V'' + p\bigl(\tfrac1{24}p^2V''' - \tfrac16V'V''\bigr) = \tfrac1{24}p^3V''' - \tfrac14pV'V''$。
 さらに $[Z_1, G_1]$: $(G_1)_p = \tfrac18p^2V''' - \tfrac14V'V''$、$(G_1)_x = \tfrac1{24}p^3V'''' - \tfrac14p(V''^2 + V'V''')$ なので $G_5 \equiv [Z_1,[Z_1,C_3]] \mapsto (-V')(\tfrac18p^2V''' - \tfrac14V'V'') + p\bigl(\tfrac1{24}p^3V'''' - \tfrac14p(V''^2 + V'V''')\bigr)$。整理して
$$
G_5 = \tfrac14V'^2V'' - \Bigl(\tfrac38V'V''' + \tfrac14V''^2\Bigr)p^2 + \tfrac1{24}V''''p^4
$$
5. $H_5 = s_5F_5 + \rho\,G_5$ の $p$ のべきごとの係数 (係数は $F_5$、$G_5$ の係数に $s_5$、$\rho$ を掛けて足す):
 $p^0$: $V'^2V''$ の係数 $-\tfrac{s_5}{60} + \tfrac\rho4 = 0.0881908 - 0.0361009 = 0.0520899$。
 $p^2$: $V'V'''$ の係数 $\tfrac{s_5}{80} - \tfrac{3\rho}8 = -0.0661431 + 0.0541513 = -0.0119918$。$V''^2$ の係数 $\tfrac{s_5}{240} - \tfrac\rho4 = -0.0220477 + 0.0361009 = 0.0140532$。
 $p^4$: $V''''$ の係数 $-\tfrac{7s_5}{5760} + \tfrac\rho{24} = 0.0064305 - 0.0060168 = 0.0004137$。

(係数は $2^{1/3}$ を含む代数的な数 (有理数係数の $w_1 = 1/(2 - 2^{1/3})$ の多項式) で、小数は有効 6 桁。たとえば $s_5 = 2(1 - 2^{2/3})/(2 - 2^{1/3})^5$。)

**確認**: 自由代数で、$w_1$ を入れた $S_4$ の $\log$ の 1 次は $A + B$ ($2 w_1 + w_0 = 1$)、2・3・4 次は 0 (係数 $10^{-41}$ 以下)、5 次は残る (語の数 32、最大の係数 0.10)。5 次の項を、上の $(2w_1^5 + w_0^5)C_5 + \rho[Z_1,[Z_1,C_3]]$ と比べると、差は $3\times10^{-41}$ 以下 ([integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (c))。5 次の関数は、具体的な表現 ($A = -p^2/2$、$B = -V(x)$) で

$$
H_5^{\rm Yoshida} = 0.05209\,V'^2V''\,p^0 + \left(-0.01199\,V'V''' + 0.01405\,V''^2\right)p^2 + 0.000414\,V''''\,p^4
\tag{3.7}
$$

($V' = V_x$ など。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (b) の直接計算と、上の 5. の手計算が一致。) **4 次の誤差は、これらの項の和で決まる**。

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

を取る。具体的な表現 ($A = -:T:$、$B = -:V:$) では、$[A, B]$ は $-(P/m)\cdot\nabla V$ の流れ、$W$ は $|\nabla V|^2/m$ (座標だけの関数)、$W_5$ は $-2\,V_jV_kV_{jk}/m^2$ (座標だけの関数。$V_j = \partial V/\partial x_j$、$V_{jk}$ は 2 階微分、重複する添字は和) に比例する。

導出 (1 自由度、$T = p^2/2$、$m = 1$。ポアソン括弧 $[f, g] = f_xg_p - f_pg_x$。リー演算子は $[:f:, :g:] = :[f,g]:$ と、$[-:f:, -:g:] = [:f:, :g:]$ から):

1. $[A, B] = [:T:, :V:] = :[T, V]:$。$[T, V] = T_xV_p - T_pV_x = 0 - p\,V' = -pV'$。よって $[A, B] = :-pV':$。
2. $W = [B, [A, B]] = [-:V:, :-pV':] = -:[V, -pV']: = :[V, pV']:$。$[V, pV'] = V_x\,\partial_p(pV') - V_p\,\partial_x(pV') = V'\cdot V' - 0 = V'^2$。よって $W = :V'^2:$ (座標だけの関数 $V'^2 = |\nabla V|^2$)。
3. $W_5 = [W, [A, B]] = [:V'^2:, :-pV':] = :[V'^2, -pV']:$。$[V'^2, -pV'] = (V'^2)_x\cdot(-V') - 0 = 2V'V''\cdot(-V') = -2V'^2V''$。よって $W_5 = :-2V'^2V'':$。

(sympy: [integ_03_expansions.out.txt](outputs/integ_03_expansions.out.txt) (b)。)

**多自由度** ($T = \sum_jp_j^2/(2m)$、$V(x_1, x_2, \ldots)$、$V_j = \partial V/\partial x_j$、$V_{jk} = \partial^2V/\partial x_j\partial x_k$) でも、括弧に和が付くだけで、同じ 3 ステップになる:

1. $[T, V] = \sum_j\bigl(T_{x_j}V_{p_j} - T_{p_j}V_{x_j}\bigr) = -\sum_j(p_j/m)\,V_j$ ($T_{x_j} = 0$、$V_{p_j} = 0$ なので第 1 項は 0。$T_{p_j} = p_j/m$)。よって $[A,B] = :-\sum_jp_jV_j/m:$。
2. $W = [B,[A,B]] = -:[V,[T,V]]:$ (1 自由度と同じ符号の計算)。$[V, g] = \sum_jV_j\,\partial g/\partial p_j$ ($V$ が座標だけなので) に $g = -\sum_kp_kV_k/m$ ($\partial g/\partial p_j = -V_j/m$) を入れて $[V,[T,V]] = -\sum_jV_j^2/m$。したがって $W = :\sum_jV_j^2/m:$。
3. $W_5 = [W,[A,B]] = :[\sum_jV_j^2/m,\ -\sum_kp_kV_k/m]:$。$[f, g] = \sum_l(f_{x_l}g_{p_l} - f_{p_l}g_{x_l})$ で、$f = \sum_jV_j^2/m$ は座標だけなので $f_{p_l} = 0$、$f_{x_l} = 2\sum_jV_jV_{jl}/m$、$g_{p_l} = -V_l/m$。よって $W_5 = :-2\sum_{jl}V_jV_lV_{jl}/m^2:$。

1 自由度の結果は $V_j \to V'$ とした特別な場合。2 自由度の関数表現は [integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (e) で確認した。**座標だけの関数** $f(\vec x)$ の流れは、運動量を $\Delta p_j = \partial f/\partial x_j$ だけ変える (座標は不変) キックで、閉じた式。

### 8.2 3 次の条件 (手で)

(3.9) を、Strang の公式 (3.2) を 2 回使って組む。

**内側**: 中央の $e^{tA}e^{C}e^{tA}$ は、(3.1) で $A \to 2tA$、$B \to C$ とした形なので、(3.2) を使うと $\exp\bigl(2tA + C - \tfrac1{24}[2tA,[2tA, C]] + \tfrac1{12}[C,[C, 2tA]] + \cdots\bigr)$。$A$ は 1 次、$C$ は 3 次なので、$[A,[A,C]]$ は 5 次、$[C,[C,A]]$ は 7 次。したがって $e^{tA}e^{C}e^{tA} = \exp(2tA + C + O(5))$ ($[A, C]$ のような 4 次の項は、前後対称なので出ない)。これを $e^{B/2}\,e^{2tA + C}\,e^{B/2}$ ((3.1) で $A \to B$、$B \to 2tA + C$、$h = 1$ の形) に入れ、(3.2) を使う:

$$
Y = B + (2tA + C) - \tfrac1{24}[B,[B, 2tA]] + \tfrac1{12}[2tA,[2tA, B]] + O(5)
$$

($C$ は 3 次なので、(3.2) の 3 次の交換子の中では 5 次以上。) $[B,[B, 2tA]] = 2t[B,[B,A]]$、$[2tA,[2tA, B]] = 4t^2[A,[A,B]]$、$C = c_3[B,[A,B]] = -c_3[B,[B,A]]$ ($[A,B] = -[B,A]$)。係数をまとめる: $-\tfrac1{24}\cdot 2t = -\tfrac t{12}$、$\tfrac1{12}\cdot4t^2 = \tfrac{t^2}3$。よって

$$
Y = B + 2tA + \Bigl(\tfrac{t^2}3\,[A,[A,B]] - \bigl(\tfrac{t}{12} + c_3\bigr)[B,[B,A]]\Bigr) + O(5)
$$

**外側** $e^{aA}\,e^{Y}\,e^{aA}$ ((3.1) で $A \to 2aA$、$B \to Y$): $X \equiv 2aA$ として、$\exp\bigl(X + Y - \tfrac1{24}[X,[X,Y]] + \tfrac1{12}[Y,[Y,X]]\bigr)$。3 次の項だけを集める ($Y$ の 1 次部分 $B + 2tA$ を交換子の中に入れる)。$[X,[X,Y]] = 4a^2[A,[A,B]]$ ($2tA$ は $A$ と交換)。$[Y,[Y,X]] = 2a\,[B + 2tA,[B, A]] = 2a[B,[B,A]] - 4at[A,[A,B]]$。各項の係数は、$[A,[A,B]]$: $\tfrac{t^2}3$ (内側 $Y$ から)、$-\tfrac1{24}\cdot4a^2 = -\tfrac{a^2}6$、$\tfrac1{12}\cdot(-4at) = -\tfrac{at}3$。$[B,[B,A]]$: $-(\tfrac t{12} + c_3)$ (内側から)、$\tfrac1{12}\cdot2a = \tfrac a6$。$a = \tfrac12 - t$ を代入すると、$[A,[A,B]]$: $\tfrac{t^2}3 - \tfrac{(1/2 - t)^2}6 - \tfrac{(1/2 - t)t}3 = \tfrac{t^2}3 - \tfrac1{24} + \tfrac t6 - \tfrac{t^2}6 - \tfrac t6 + \tfrac{t^2}3 = \tfrac{t^2}2 - \tfrac1{24}$。$[B,[B,A]]$: $-\tfrac t{12} - c_3 + \tfrac1{12} - \tfrac t6 = \tfrac1{12} - \tfrac t4 - c_3$。したがって、$\log S_{\rm FG}$ の 3 次の項の係数は

$$
\begin{aligned}
[A,[A,B]]:&\quad \tfrac{t^2}3 - \tfrac{a^2}6 - \tfrac{at}3 = \tfrac{12t^2 - 1}{24}\quad(a = \tfrac12 - t\text{ を代入。上の計算}),\\
[B,[B,A]]:&\quad -\tfrac t{12} - c_3 + \tfrac a6 = \tfrac{1}{12} - \tfrac t4 - c_3
\end{aligned}
\tag{3.11}
$$

(sympy の自由代数による直接計算で、同じ式 $c_1 = (12t^2 - 1)/24$、$c_2 = -(12c_3 + 3t - 1)/12$ を得た。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (c)。係数の代数は [integ_03_expansions.out.txt](outputs/integ_03_expansions.out.txt) (c) でも確認。) **3 次が消える条件**は

$$
12t^2 = 1\ \Rightarrow\ t = \frac1{\sqrt{12}} = 0.288675\ldots,\qquad
c_3 = \frac1{12} - \frac t4 = \frac1{12} - \frac1{4\sqrt{12}} = 0.0111645\ldots
\tag{3.12}
$$

($t > 0$ の解。$t = -1/\sqrt{12}$ は負のドリフトで採らない。) これは、K64 の `alpha = 1/sqrt(12)`、`alpha1 = 1/2 - alpha` ($= a$)、そして `beta` $= \tfrac16 - \tfrac\alpha2 = 0.0223291$ の**半分**: $c_3 = \beta/2$ (数値で $c_3 = 0.011164549685$、$\beta/2 = 0.011164549685$)。

**意味**: ドリフトの長さの割り振り ($t = 1/\sqrt{12}$) は、$[A,[A,B]]$ (運動量の項。$p^2V'''$ 型) を消す。補正キック $c_3 W$ は、$[B,[B,A]]$ ($|\nabla V|^2$ 型の力の勾配の項) を消す。補正キックは、座標だけの関数の流れで、**$B$ が座標だけの関数だから閉じた式で作れる**。

### 8.3 5 次の条件 ($c_5$)

3 次を消した後の 5 次の項を、§4.1 の補題 S (3.2b) を 3 回使って組む (内側 → 中間 → 外側。§8.2 の 3 次と同じ構成)。6 個の係数を 1 つずつ追う。確認は [integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (d)。

**記号**。補題 S の基底 $c_1\ldots c_6$ ((3.2a) の表) で $X = A$、$Y = B$ としたものを $b_1\ldots b_6$ と書く:

$$
U = [A,B],\quad V = [A,U],\quad W = [B,U]\ (\text{(3.10) の }W),\qquad
b_1 = [A,[A,V]],\ b_2 = [B,[A,V]],\ b_3 = [U,V],\ b_4 = [A,[B,W]],\ b_5 = [U,W],\ b_6 = [B,[B,W]]
$$

(以下 $V$ は括弧 $[A,[A,B]]$ で、ポテンシャルではない。) ヤコビ恒等式は §4.1 の (J1)〜(J5) で $X = A$、$Y = B$ とおいたもの: $[A,W] = [B,V]$、$[A,[B,V]] = b_3 + b_2$、$[B,[B,V]] = b_4 - b_5$、$[A,[A,W]] = b_2 + b_3$、$[B,[A,W]] = b_4 - b_5$。補正は $C = c_3W + c_5W_5$、$W_5 = [W,U] = -b_5$。

**入れ替えの規則**: 中間の段で補題 S を $(X, Y) = (B, A)$ として使う ($e^{B/2}\ldots e^{B/2}$ の外側が $B$)。このとき補題の $U$、$V$、$W$ は $[B,A] = -U$、$[B,-U] = -W$、$[A,-U] = -V$ で、基底は符号つきで入れ替わる:

$$
c_1 \to -b_6,\quad c_2 \to -b_4,\quad c_3 \to b_5,\quad c_4 \to -b_2,\quad c_5 \to b_3,\quad c_6 \to -b_1
$$

(たとえば $c_1 = [X,[X,V_{XY}]] \to [B,[B,-W]] = -b_6$、$c_3 = [U_{XY},V_{XY}] \to [-U,-W] = b_5$。)

**内側** ($e^{tA}e^{C}e^{tA}$。補題 S で $X = 2tA$、$Y = C$): $C$ は 3 次以上なので、$-\tfrac1{24}[X,[X,C]]$ だけが 5 次に残る ($[C,[C,X]]$ と $C_5$ は 7 次以上)。

$$
Y_1 = 2tA + C - \tfrac{t^2}6c_3[A,[A,W]] + O(7) = 2tA + c_3W - c_5\,b_5 - \tfrac{t^2c_3}6\,(b_2 + b_3) + O(7)
$$

($[A,[A,W]] = b_2 + b_3$。)

**中間** ($e^{B/2}e^{Y_1}e^{B/2}$。補題 S で $X = B$、$Y = Y_1$)。$Y_1$ の 1 次部分 $2tA$、3 次部分 $c_3W$、5 次部分 $Y_1^{(5)}$ で、$Y_2 = B + Y_1 + C_3(B, Y_1) + C_5(B, Y_1)$ の次数ごとの項:

- 1 次: $B + 2tA$。
- 3 次: $c_3W + C_3(B, 2tA)$。$C_3(B, 2tA) = -\tfrac1{24}[B,[B,2tA]] + \tfrac1{12}[2tA,[2tA,B]]$、$[B,[B,A]] = -[B,U] = -W$、$[A,[A,B]] = V$ から $\tfrac t{12}W + \tfrac{t^2}3V$。したがって $Y_2^{(3)} = \tfrac{t^2}3V + \bigl(c_3 + \tfrac t{12}\bigr)W$。
- 5 次: 3 つの寄与。
 (i) $Y_1^{(5)} = -c_5b_5 - \tfrac{t^2c_3}6(b_2 + b_3)$。
 (ii) $C_3$ の中の $Y_1$ の 3 次部分 $c_3W$ を 1 つ含む項: $-\tfrac1{24}[B,[B,c_3W]] + \tfrac1{12}\bigl([c_3W,[2tA,B]] + [2tA,[c_3W,B]]\bigr)$。$[B,[B,W]] = b_6$、$[W,[A,B]] = [W,U] = -b_5$、$[A,[W,B]] = -[A,[B,W]] = -b_4$ なので、$-\tfrac{c_3}{24}b_6 - \tfrac{tc_3}6(b_4 + b_5)$。
 (iii) $C_5(B, 2tA)$。補題 S の $c_i$ を、入れ替えの規則で $b$ に直し、$Y$ に当たる $2tA$ の個数 ($c_1$: 1 個、$c_2$、$c_3$: 2 個、$c_4$、$c_5$: 3 個、$c_6$: 4 個) だけ $2t$ のべきをかける: $-\tfrac7{5760}(2t)\,b_6 - \tfrac7{1440}(2t)^2\,b_4 + \tfrac1{480}(2t)^2\,b_5 - \tfrac1{180}(2t)^3\,b_2 + \bigl(-\tfrac1{120}\bigr)(2t)^3\,b_3 - \tfrac1{720}(2t)^4\,b_1$。
 3 つを足すと
$$
Y_2^{(5)} = -\tfrac{t^4}{45}b_1 - \Bigl(\tfrac{2t^3}{45} + \tfrac{t^2c_3}6\Bigr)b_2 - \Bigl(\tfrac{t^3}{15} + \tfrac{t^2c_3}6\Bigr)b_3 - \Bigl(\tfrac{7t^2}{360} + \tfrac{tc_3}6\Bigr)b_4 + \Bigl(\tfrac{t^2}{120} - \tfrac{tc_3}6 - c_5\Bigr)b_5 - \Bigl(\tfrac{7t}{2880} + \tfrac{c_3}{24}\Bigr)b_6
$$
 (たとえば $b_2$: (iii) の $-\tfrac1{180}\cdot8t^3 = -\tfrac{2t^3}{45}$ と (i) の $-\tfrac{t^2c_3}6$)。

**外側** ($e^{aA}e^{Y_2}e^{aA}$。補題 S で $X = 2aA$、$Y = Y_2$、$a = \tfrac12 - t$)。$Y_2$ の 1 次部分 $B + 2tA$、3 次部分 $Y_2^{(3)}$、5 次部分 $Y_2^{(5)}$ を使う。$X$ は $A$ の定数倍なので、$[X, 2tA] = 0$。

- 1 次: $2aA + B + 2tA = A + B$ ($a + t = \tfrac12$)。
- 3 次: $Y_2^{(3)} + C_3(2aA, B + 2tA)$。$-\tfrac1{24}[X,[X,B]] = -\tfrac{a^2}6V$。$[B + 2tA, X] = -2aU$ なので $[Y,[Y,X]] = -2a[B + 2tA, U] = -2a(W + 2tV)$、$\tfrac1{12}$ をかけて $-\tfrac a6W - \tfrac{at}3V$。よって $V$ の係数は $\tfrac{t^2}3 - \tfrac{a^2}6 - \tfrac{at}3 = \tfrac{12t^2 - 1}{24}$、$W$ の係数は $c_3 + \tfrac t{12} - \tfrac a6 = c_3 + \tfrac t4 - \tfrac1{12}$ ((3.11) と同じ。$[B,[B,A]] = -W$)。
- 5 次: 3 つの寄与。
 (i) $Y_2^{(5)}$ (上の式)。
 (ii) $C_3$ の中の $Y_2^{(3)}$ を 1 つ含む項: $-\tfrac1{24}[X,[X,Y_2^{(3)}]] + \tfrac1{12}\bigl([Y_2^{(3)},[Y_2^{(1)},X]] + [Y_2^{(1)},[Y_2^{(3)},X]]\bigr)$。$\kappa \equiv c_3 + \tfrac t{12}$ とおくと $Y_2^{(3)} = \tfrac{t^2}3V + \kappa W$。$[A,[A,V]] = b_1$、$[A,[A,W]] = b_2 + b_3$、$[B,[A,V]] = b_2$、$[B,[A,W]] = b_4 - b_5$、$[V,U] = -b_3$、$[W,U] = -b_5$ から、
  $[X,[X,Y_2^{(3)}]] = 4a^2\bigl[\tfrac{t^2}3b_1 + \kappa(b_2 + b_3)\bigr]$、
  $[Y_2^{(3)},[Y_2^{(1)},X]] = 2a\bigl[\tfrac{t^2}3b_3 + \kappa b_5\bigr]$ ($[Y_2^{(1)},X] = -2aU$ なので $-2a[Y_2^{(3)},U]$)、
  $[Y_2^{(1)},[Y_2^{(3)},X]] = -2a\Bigl[\tfrac{t^2}3(b_2 + 2tb_1) + \kappa\bigl(b_4 - b_5 + 2t(b_2 + b_3)\bigr)\Bigr]$ ($[Y_2^{(3)},X] = -2a[A,Y_2^{(3)}]$、$[B + 2tA,\cdot]$ を作用)。
 (iii) $C_5(2aA, B + 2tA)$。補題 S の括弧を $X = 2aA$、$Y = B + 2tA$ で作る: $U_o = [X,Y] = 2aU$、$V_o = [X,U_o] = 4a^2V$、$W_o = [Y,U_o] = 2a(W + 2tV)$ ($[A,A] = 0$ による)。
  $c_1^o = 16a^4\,b_1$、$c_2^o = 8a^3(b_2 + 2tb_1)$、$c_3^o = 8a^3\,b_3$、$c_4^o = 4a^2\bigl(b_4 + 4t(b_2 + b_3) + 4t^2b_1\bigr)$、$c_5^o = 4a^2(b_5 + 2tb_3)$、$c_6^o = 2a\bigl(b_6 + 6tb_4 - 4tb_5 + 12t^2b_2 + 8t^2b_3 + 8t^3b_1\bigr)$
  (たとえば $c_4^o = [X,[Y,W_o]]$ は、$[B + 2tA,\,W + 2tV] = [B,W] + 4t[B,V] + 4t^2[A,V]$ ($[A,W] = [B,V]$) に $2a\cdot2a[A,\cdot]$ を作用: $[A,[B,W]] = b_4$、$[A,[B,V]] = b_3 + b_2$、$[A,[A,V]] = b_1$)。$C_5 = \tfrac7{5760}c_1^o + \tfrac7{1440}c_2^o + \tfrac1{480}c_3^o + \tfrac1{180}c_4^o - \tfrac1{120}c_5^o + \tfrac1{720}c_6^o$。

$t = 1/\sqrt{12}$、$c_3 = \tfrac1{12} - \tfrac t4$、$a = \tfrac12 - t$ を入れた 3 つの寄与 ($b_1 \ldots b_6$ の係数。$s \equiv \sqrt3$):

| 寄与 | $b_1$ | $b_2$ | $b_3$ | $b_4$ | $b_5$ | $b_6$ |
|---|---|---|---|---|---|---|
| (i) $Y_2^{(5)}$ | $-\tfrac1{6480}$ | $-\tfrac1{864} - \tfrac{s}{25920}$ | $-\tfrac1{864} - \tfrac s{2880}$ | $\tfrac1{540} - \tfrac s{432}$ | $\tfrac1{240} - \tfrac s{432} - c_5$ | $-\tfrac1{288} + \tfrac{23s}{17280}$ |
| (ii) | $-\tfrac1{1296}$ | $-\tfrac1{216} + \tfrac s{648}$ | $0$ | $-\tfrac1{108} + \tfrac s{216}$ | $\tfrac1{54} - \tfrac s{108}$ | $0$ |
| (iii) | $\tfrac1{1440}$ | $\tfrac1{270} - \tfrac s{1080}$ | $\tfrac7{2160} - \tfrac s{720}$ | $\tfrac{13}{2160} - \tfrac s{432}$ | $-\tfrac{11}{1080} + \tfrac s{216}$ | $\tfrac1{720} - \tfrac s{2160}$ |
| **合計** | $-\tfrac1{4320}$ | $-\tfrac1{480} + \tfrac s{1728}$ | $\tfrac1{480} - \tfrac s{576}$ | $-\tfrac1{720}$ | $\tfrac1{80} - \tfrac s{144} - c_5$ | $-\tfrac1{480} + \tfrac s{1152}$ |

(例 $b_1$: $-\tfrac1{6480} - \tfrac1{1296} + \tfrac1{1440} = \tfrac{-2 - 10 + 9}{12960} = -\tfrac1{4320}$。$b_4$ の $s$ の項: $-\tfrac1{432} + \tfrac1{216} - \tfrac1{432} = 0$。)

一般の $t$、$c_3$、$c_5$ での合計は

$$
\begin{aligned}
b_1&: \tfrac{240t^4 - 120t^2 + 7}{5760},\qquad b_2: \tfrac{-60c_3 + 60t^3 - 60t^2 - 15t + 7}{1440},\qquad b_3: -\tfrac{20c_3 + 60t^3 - 20t^2 + 5t - 1}{480},\\
b_4&: -\tfrac{60c_3 + 15t - 4}{720},\qquad b_5: -\tfrac{60c_3t - 20c_3 + 120c_5 + 5t^2 - 5t + 1}{120},\qquad b_6: -\tfrac{120c_3 + 15t - 4}{2880}
\end{aligned}
$$

(自由代数で $\log S_{\rm FG}$ を直接計算した 5 次の係数と一致。[integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (d))。

**関数表現**。§4.1 (d) の表から、$b_1 \mapsto -p^4V''''$、$b_2 \mapsto 3p^2V'V'''$、$b_3 \mapsto p^2(2V''^2 - V'V''')$、$b_4 \mapsto 0$、$b_5 \mapsto 2V'^2V''$、$b_6 \mapsto 0$ (1 自由度。$b_i \mapsto$ は $c_i$ の $X = A$、$Y = B$ の場合)。上の合計を掛けて足すと

- $p^4$: $-b_1V'''' = \tfrac{V''''}{4320}$。
- $p^2$: $V'V'''$ の係数 $3b_2 - b_3 = -\tfrac3{480} + \tfrac{3s}{1728} - \tfrac1{480} + \tfrac s{576} = -\tfrac1{120} + \tfrac{\sqrt3}{288}$、$V''^2$ の係数 $2b_3 = \tfrac1{240} - \tfrac{\sqrt3}{288}$。
- $p^0$: $2b_5\,V'^2V'' = 2\bigl(\tfrac1{80} - \tfrac{\sqrt3}{144} - c_5\bigr)V'^2V'' = \tfrac{-720c_5 + 9 - 5\sqrt3}{360}V'^2V''$ ($\tfrac2{80} = \tfrac9{360}$、$\tfrac{2\sqrt3}{144} = \tfrac{5\sqrt3}{360}$)。

すなわち、3 次を消した後の 5 次の項は、

$$
H_5^{\rm FG} = \frac{-720\,c_5 + 9 - 5\sqrt3}{360}\,V'^2V''\,p^0 + \left[\Bigl(\tfrac{\sqrt3}{288} - \tfrac1{120}\Bigr)V'V''' + \Bigl(\tfrac1{240} - \tfrac{\sqrt3}{288}\Bigr)V''^2\right]p^2 + \frac{V''''}{4320}\,p^4
\tag{3.13}
$$

(sympy の自由代数の直接計算 ([integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (c)) と、この手計算の結果の差は 0。**独立な数値検証**も足した: $H = p^2/2 + x^4/4$ ($V' = x^3$、$V'' = 3x^2$、$V''' = 6x$、$V'''' = 6$) で、(3.9) の 1 刻み ($c_3 = \tfrac1{12} - \tfrac t4$、$c_5$ は下の (3.14) の値、補正キックは $f_C = c_3h^3V'^2 - 2c_5h^5V'^2V''$ で $\Delta p = \partial f_C/\partial x$) と厳密な流れ (mpmath 50 桁の ODE) の差を $h^5$ で割ると、$h = 0.04, 0.02, 0.01$ で $x$ 成分 $0.001017, 0.001006, 0.001001$、$p$ 成分 $-0.001416, -0.001395, -0.001385$ と、$h \to 0$ で $(0.000995,\ -0.001374)$ に収束する ($h$ の 1 次の補正を含む)。これは (3.13) の $H_5$ を $x_0 = 0.5$、$p_0 = 0.3$ で評価した $(\partial H_5/\partial p,\ -\partial H_5/\partial x) = (-0.000995,\ +0.001374)$ の**符号を反転した値**。すなわち、1 刻みの誤差は $-h^5\times(H_5\text{ の流れ})$ で、大きさと成分の比は (3.13) と一致する (分割した写像は、$H - h^4H_5$ の流れに近い)。$c_5$ を 0 にしたときも同じ検証をして、(3.13) の $c_5$ 依存 ($p$ 成分 $-0.001204$ 対 $-0.001197$) が合う ([integ_03_expansions.out.txt](outputs/integ_03_expansions.out.txt) (d))。)

- $p^2$、$p^4$ の項は、$C$ に含まれる座標だけの関数ではどうにもならない (運動量を含むので)。これは 4 次の誤差として残る。
- $c_5$ が $H_5$ に入る係数は、$C \ni c_5h^5W_5 = :-2c_5h^5V'^2V'':$ が $-2c_5\,V'^2V''\,p^0$ として入るだけ ($c_5$ は上の $Y_1$ の $-c_5b_5$ にだけ入り、$b_5 \mapsto 2V'^2V''$)。(3.13) の $V'^2V''p^0$ の係数 $(-720c_5 + 9 - 5\sqrt3)/360$ の $c_5$ の部分 $-720/360 = -2$ と一致する。
- **$p^0$ の項 ($V'^2V''$) は、座標だけの関数なので、$c_5$ を使って消せる**: $-720c_5 + 9 - 5\sqrt3 = 0$ より

$$
c_5 = \frac{9 - 5\sqrt3}{720} = \frac1{80} - \frac{\sqrt3}{144} = 4.71869 \times 10^{-4} = \frac{\gamma}{2}
\tag{3.14}
$$

 (K64 の `gamma` $= \tfrac1{40} - \tfrac1{24\sqrt3} = 9.43739 \times 10^{-4}$ の半分。$\tfrac1{24\sqrt3} = \tfrac{\sqrt3}{72}$。比 $c_5/\gamma = 0.500000000$。)

**多自由度で、5 次の座標だけの項が $W_5$ だけであることの証明**。(a) 6 個の $b_i$ の関数表現 (2 自由度で計算): $b_4 = [A,[B,W]]$、$b_6 = [B,[B,W]]$ は、$[B,W]$ が座標だけの関数 ($B \mapsto -V$、$W \mapsto \sum_jV_j^2/m$) どうしの括弧で 0 なので、関数が 0。$b_1$、$b_2$、$b_3$ は運動量の 4 次、2 次、2 次の項だけを含み、$b_5 = -W_5 \mapsto +2\sum_{jk}V_jV_kV_{jk}/m^2$ は運動量を含まない (2 自由度で確認、[integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (e))。(b) ほかの 5 次の括弧は、ヤコビ恒等式で基底に直る。たとえば $[B,[B,[A,[A,B]]]] = [B,[B,V]]$ は、(J3) で $b_4 - b_5$ になる ($[B,[B,V]] = [B,[A,W]]$ ((J1)) $= [[B,A],W] + [A,[B,W]] = -[U,W] + [A,[B,W]] = b_4 - b_5$)。関数表現で $b_4 \mapsto 0$ なので、$-b_5 = [W,U] = W_5$ と等しい (1 自由度で、両方とも $-2V'^2V''$、比 1)。よって 5 次の運動量を含まない項は、基底では $b_5$ のみ。1 自由度の式 $V'^2V''$ は、多自由度の $G \equiv \sum_{jk}V_jV_kV_{jk}$ の特別な場合 ($W_5 \mapsto -2G/m^2$)。

**この $\gamma$ の項は 5 次の誤差の一部 ($p^0$) だけを消す**ので、全体の誤差が小さくなるとは限らない: この章の例 (1 自由度の $x^4/4$) では、$c_5$ を入れない方 (c3 だけ) が誤差が小さい ([integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) (a)(d)。$n = 64$ で $1.5 \times 10^{-9}$ 対 $6.5 \times 10^{-9}$)。多極 ($V \propto \mathrm{Re}\,(x - iy)^{n+1}$) では、07 章で確認する。

### 8.4 補正キックの具体形と K64 との一致

質量に当たるものが全運動量 $p$ ($T = (p_x^2 + p_y^2)/(2p)$。$A$ の $H_a$ は $p$ に依る) のとき、補正キックの母関数 ($\exp(:f:)$ の $f$。$\Delta p_j = \partial f/\partial x_j$、$\Delta z = -\partial f/\partial p$。[C_poisson_lie.md](C_poisson_lie.md) の規約) は

$$
f_C = c_3\,h^3\,\frac{|\nabla V|^2}{p} - 2c_5\,h^5\,\frac{V_jV_kV_{jk}}{p^2}
\tag{3.15}
$$

(1 自由度で $m = p$ とすると、$W = V'^2/m$、$W_5 = -2V'^2V''/m^2$ で $f_C = c_3 W + c_5 W_5$)。規格化した運動量 $q = p_x/p$ の更新は $\Delta q = \Delta p_x/p$。$V = \dfrac{K}{(k+1)!}\mathrm{Re}\,(x - iy)^{k+1}$ (多極のポテンシャル。07 章) として、(3.15) から作った $\Delta q_x$、$\Delta q_y$、$\Delta z$ は、K64 の `tthin` の補正キックの式 (`tquad.f:208-250` の $f_1$〜$f_5$ を使う式) と、$k = 2$ (6 極)、$k = 3$ (8 極) の**記号計算で完全に一致した** (差 0。$K$、$p$ は記号のまま。[integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) (d))。具体形は [07_thin_multipole.md](07_thin_multipole.md) で書く。

**近軸の質量を使うことの限界** ([integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) (f))。(3.15) は、$A$ を近軸のドリフト $T = (p_x^2 + p_y^2)/(2p)$ とみて $W = |\nabla V|^2/p$ を使った形で、K64 の補正キックの式と一致する (上)。$A$ を**厳密なドリフト** $H_a = -p_z$ ($p_z = \sqrt{P^2 - p_x^2 - p_y^2}$、全運動量 $P$ は定数) にすると、1 自由度 ($p_y = 0$) で $[H_a, V] = -(p_x/p_z)V'$ ($V$ の $p_x$ 微分は 0 で $H_a$ の $x$ 微分も 0。$\partial H_a/\partial p_x = p_x/p_z$)、$W = -[V,[H_a,V]] = V'^2\,\partial_{p_x}(p_x/p_z) = V'^2P^2/p_z^3$ ($\partial_{p_x}(p_x/p_z) = 1/p_z + p_x^2/p_z^3 = P^2/p_z^3$)。近軸の $V'^2/P$ との比は $P^3/p_z^3 = (1 - u)^{-3/2} \approx 1 + \tfrac32u$ ($u = p_x^2/P^2$ は傾きの 2 乗)。補正キックは座標だけの関数の流れ ($p_x$ に依る $W$ は、座標も動かす) なので、K64 の補正キックは近軸の形にとどまり、厳密なドリフトに対する $[B,[B,A]]$ の消去は相対 $\tfrac32u$ だけ不完全で、$c_3h^3\cdot\tfrac32u\,V'^2/P$ 程度の 3 次の残りが出る (傾きの小さい粒子では小さい)。この残りが実際の追跡誤差にどれだけ効くかは、この章では評価していない (**未確認**。07 章の数値は傾きの小さい粒子と大きい粒子 ($p_x = 30$ mrad、$p_y = -20$ mrad) で確認している)。

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
| 四極の追跡 (`tsolqu`、`tsolqu.f:47`) | $n = 1 + \lfloor L\,\lVert(K1, B_z)\rVert/\varepsilon\rfloor$ (`ak` = K1 の積分強度、`bz` = ソレノイド) | $\varepsilon = 0.2\,\mathrm{EPS}$ (`epsdef = 0.2d0`、`tsolqu.f:14`) | 04 章 §6.3 |
| 四極の CALC (`tsolque`、`tsolque.f:385`) | 同じ式。上限 `ndivmax` | $0.1\,\mathrm{EPS}$ (`tsolque.f:381-383`) | 追跡の半分の $\varepsilon$ |
| 放射つきの四極の追跡 (`tsolqu.f:400-402`) | $n = \max\bigl(1 + \lfloor L\,K/\varepsilon\rfloor,\ \mathrm{ndivrad}(\lVert(a_{k0x}, a_{k0y})\rVert, K, B_z, \mathrm{EPS})\bigr)$、上限 `ndivmax` | $0.2\,\mathrm{EPS}$ | `ak0x`、`ak0y` はステアリング (偏向の積分強度) |
| 偏向磁石 (K1 ≠ 0) の CALC (`tbende.f:460`) | $n = 1 + \lfloor\lvert\phi_c/\varepsilon\rvert\rfloor$ ($K1 = 0$ なら 1) | `epsbend` $\times$ EPS (`eps0 == 0` なら `epsbend` そのまま。`epsbend` の値は未確認。確認できたのは `tbende.f:454-460` の式まで) | $\phi_c$ = `phi0*rbc` (偏向角に当たる量。`rbc` の定義は未確認)。放射・スピン偏極があれば、ほかの条件との大きい方 |
| 放射つきの偏向磁石の追跡 (`tbend.f:640`) | $n = \min(\mathrm{ndivmax}, \mathrm{ndivrad}(\phi_r, 0, 0, \varepsilon))$ | EPS | $\phi_r = \phi_b\,(L - f_{1r} - f_{2r})/L$ (端の分を除いた偏向角) |
| 多極 MULT の追跡 (`tmulti.f:107-124`) | $n = \max(\mathrm{ndiv0}, \mathrm{ndivrad}(\ldots), \ldots)$、上限 `ndivmax` | EPS | `ndiv0` は要素のパラメータ (名前の対応は未確認) |
| 薄い多極 (`tthin`) | $n = 1$ (本体は 1 刻み) | — | `call setndivelm(l_track, 1)` |
| `ndivrad(phi, ak, bz, eps)` (`tffs.f:641-653`) | $\lfloor 1 + (h_0/100)\,[(\lvert\phi\rvert + \lVert(K, B_z)\rVert\cdot 0.02)/\varepsilon_1]\rfloor$ ($\varepsilon_1 = $ EPS (EPS ≤ 0 なら 1)。`trpt` が真なら $\times 10$) | EPS | $h_0 = \gamma_0$ (`tmacro` の `h0`)。係数 0.02 は `arad`。放射の光子数が 1 刻みで多くなりすぎないための分割数と思われる (推定) |

**式の読み方**: 四極なら $n \propto L\,K1$ ($K1$ は積分強度)。$L\,K1 = KL^2 = \varphi^2$ は本体全体の位相の 2 乗 ($\varphi = wL$、$K = K1/L$ は長さあたりの強さ)。1 刻みの位相は $\varphi/n \approx \sqrt{\varepsilon/n}$ (近似)。誤差は 1 刻みの位相の何乗かに比例するので、$\varepsilon$ を決めることが 1 刻みの位相を決めることになる。

数値例は [04_quad.md](04_quad.md) §6.3 の表 ($L = 0.5$ m、K1 = 1.5 /m で追跡 K64 は $n = 4$、v2 は 8)。

## 11. 数値確認で使うスクリプト

| スクリプト | 確認すること | 出力 |
|---|---|---|
| `integ_01_bch_force_gradient.py` | (a) Strang の BCH (非可換な級数、sympy)。(b) Yoshida の係数と 3 次の消去、5 次の残り。(c) 力勾配型: 3 次の係数 ((3.11))、$t, c_3$ ((3.12))、$c_5$ ((3.14))、残る 5 次の項 ((3.13))。(d) 2 自由度で、導いた補正キックと K64 の `tthin` の式 ($\Delta q$、$\Delta z$) が 6 極・8 極で記号的に一致。(e) CALC の `tthine` の分割の係数 | [integ_01_bch_force_gradient.out.txt](outputs/integ_01_bch_force_gradient.out.txt) |
| `integ_04_bch5_steps.py` | 5 次の BCH の手計算の確認 (§4.1、§5、§8.3): (a) 基底 $c_1 \ldots c_6$ とヤコビ恒等式 (J1)〜(J5)、(b) 補題 S の 6 個の係数 ((c-3)、(c-4) の表の各項を含む)、(c) Yoshida の 5 次 (3.7)、(d) 力勾配型の 3 段の組み立てと (3.13)、(e) $W_5$ と 2 自由度の関数、(f) 厳密なドリフトの $W$ | [integ_04_bch5_steps.out.txt](outputs/integ_04_bch5_steps.out.txt) |
| `integ_03_expansions.py` | 書き足した展開の確認: (a) Yoshida の係数と一般化、(b) 関数表現 ($[A,B]$、$W$、$W_5$)、(c) (3.11) の係数の代数、(d) 力勾配型の 1 刻みの誤差と (3.13) の比較 (50 桁) | [integ_03_expansions.out.txt](outputs/integ_03_expansions.out.txt) |
| `integ_02_order_symplectic.py` | 誤差の次数 (Strang・Yoshida・力勾配型・RK4)、シンプレクティック性、長時間の $H$ の誤差、力の評価の数あたりの誤差 | [integ_02_order_symplectic.out.txt](outputs/integ_02_order_symplectic.out.txt) |

実行はリポジトリのルートから `py -3.12 physics_guide/scripts/integ_01_bch_force_gradient.py` (sympy。約 12 秒)、`py -3.12 physics_guide/scripts/integ_02_order_symplectic.py` (約 1 秒)。四極の積分器の次数の確認は [quad_03_integrator_order.out.txt](outputs/quad_03_integrator_order.out.txt) (この章とは重複しない)。

## 12. Fortran・pysad の対応

| 内容 | Fortran (ファイル:行) | pysad |
|---|---|---|
| 四極の Strang の構造 ($D\,B\,D$、隣り合う $D$ を 1 つに) | `tsolqu.f:61-95` (本章の範囲では刻みの構造だけ読んだ。式は 04 章) | `quadbody.group` (`quadbody.py:76`)、`ngroup` |
| Yoshida の係数 | (K64 に無い) | `_W1`、`_W0` (`quadbody.py:37-38`)、`QB_BW` (`:40`) |
| 方式の選択 (k64 は Strang、v2 は Yoshida) | — | `quadbody.scheme()`、`set_scheme` (`quadbody.py:45`)、`scheme` (`:53`) |
| 力勾配型の係数 | `tquad.f:131-135` (`alpha`、`beta`、`gamma`、`alpha1`) | `track6._TT_ALPHA`・`_TT_ALPHA1`・`_TT_BETA`・`_TT_GAMMA` |
| 力勾配型の本体 (追跡) | `tquad.f:118-` (`tthin`)。`ala = al*alpha1`、`alb = al*alpha` (`:173-174`)。補正キックは `kord = 2` が `:208-`、`kord = 1` が `:239-`、一般の `kord` が `:261-` | `track6._tthin_sext` (`track6.py:898`)、`_tthin_gen` (`:746`) |
| CALC の 2 次の分割 | `tthine.f:36-87` (`ala = al/6`、`alb = al/1.5`) | `lattice._tthine_cod` (`lattice.py:1755`)、`elkern._tthine_map` (`:360`) |
| 刻みの数の式 | §10 | `quadbody.ndiv`・`ndiv_calc` (`quadbody.py:59`・`:65`) |

(行番号は参照した K64 のソース、pysad のファイルのもの。)

## 13. 未確認・限界・使った条件

**使った条件**: 一般論の数値確認は、1 自由度 $H = p^2/2 + x^2/2 + x^4/4$ と、2 自由度の多極のポテンシャル ($k = 2, 3$) で、記号計算 (sympy) と数値積分 (DOP853)。K64 のソースとの比較は、`tthin` の補正キックの式を記号計算の結果と比べたもの (pysad の追跡を数値で動かした比較は 07 章)。

**未確認・限界**:
- **BCH 級数の収束**は議論していない。誤差の次数は、形式的な級数の係数の次数で決まる (数値で確認)。$\tilde H$ の「指数的に長い時間」の保存は、解析的な定理を引用しただけで、この章では確かめていない (数値は 5 万刻みまで)。
- **刻みの数の経験式**は、K64 のソースの式を読んだもの。$\varepsilon$ の値の物理的な根拠 (なぜ 0.2、0.1 か) は不明 (K64 の定数)。`epsbend`、`ndivmax`、`ndiv0`、$\phi_c$ の定義は未確認 (`epsrad` は $10^{-6}$ の定数、§2 の表)。
- **力勾配型の 5 次の全誤差**: $p^2$・$p^4$ の項の係数は (3.13) に導いたが、その大きさが実際の誤差に効く程度と、多極での実際の誤差は 07 章で扱う。
- **5 次の係数の導出の範囲**: 補題 S (3.2b) の 6 個の係数は §4.1 で手で導き、自由代数 (語の長さ 5 まで、sympy) の直接計算と一致する。7 次以上の項は計算していない。(3.2b) は形式的な級数の係数で、収束は議論していない。
- **5 次で消せる項が「$p^0$ の項だけ」**であること (多自由度で $G$ のみ) は、ヤコビ恒等式 (J3) と関数表現 ($[B,W] \mapsto 0$) による (§8.3)。2 自由度で 6 個の基底の関数の運動量の次数を確認した。3 自由度以上も同じ構造だが、別には計算していない。
- **厳密なドリフトに対する補正キックの近軸近似** ((3.15) の質量 $= p$) による 3 次の残り (相対 $\tfrac32u$) の大きさの評価は、この章ではしていない (§8.4)。
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
| ヘルプ SEXT・OCT・DECA・DODECA の transformation の $V_n$ | **原文で確認** (`SADHelp.html` の `transformation:THIN`)。$V_n = \sum_{j,k}\bigl(-\tfrac\beta2 H_{n,k}^2 + \gamma H_{n,j}H_{n,k}H_{n,jk}\bigr)$、$\beta = 1/6 - 1/\sqrt{48} = 0.0223291$、$\gamma = 1/40 - 1/(24\sqrt3) = 9.4374 \times 10^{-4}$。$H_n$ は多極のキックのハミルトニアン (積分強度)。ヘルプの $\exp(:\cdot:)$ は「ハミルトニアンの流れ」(C 付録 §8.2) なので、本章 (3.15) の $f_C$ (流れの生成元を $-H$ の側に書く) とは $V_n = -f_C$ の関係で、$c_3 = \beta/2$、$c_5 = \gamma/2$ と符号・係数が一致する (対称形での数値確認: [thin_04_help_order.out.txt](outputs/thin_04_help_order.out.txt))。ヘルプの式には $L$ と $1/p$ の因子が無い (略記の疑い。07 章付録 A)。ドリフトと H の並びは原文が非対称 (誤記。07 章付録 A)。 |
| ヘルプに無いもの | (1) 力勾配型の係数 $\alpha, \beta, \gamma$ の導出。(2) 補正キックが力の勾配の項を消す理由。(3) 分割の次数の議論 |

K64 の側の注意 (誤りではない): CALC の `tthine` は 2 次、追跡の `tthin` は 4 次で、同じ薄い多極でも分割が違う。したがって、軌道が大きくずれたとき、CALC の行列と追跡のヤコビアンは 2 次の分割の誤差だけ違う (07 章)。
