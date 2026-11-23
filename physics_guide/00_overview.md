# 00. 全体像: ハミルトニアンから Fortran・Python まで

**この章で分かること**: 物理ガイド全体の地図。(1) ハミルトニアンから追跡の写像・転送行列・K64 の Fortran・pysad の Python に至る**一本の流れ**。(2) K64 の**3 系統の実装** (追跡 `t*.f`、行列 `t*e.f`、解析微分 `q*`) と、要素ルーチンの命名則。(3) 各章の位置づけと読み方。(4) **ガイド全体の確認の一覧** (どの章で、どの式を、どの独立な基準で確かめたか)。(5) **ガイド全体の未確認・K64 とヘルプの誤りの疑いの一覧** (各章の付録からの集約。確定と疑いを分ける)。

**前提**: 加速器物理の基礎。SAD の内部は知らなくてよい。

**読む順**: §1 (流れ) → §2 (3 系統) → §4 (章の読み方)。確認の一覧 (§6)、疑いの一覧 (§7) は、必要なときに引く。

## 1. 流れ: ハミルトニアンから Python まで

| 段 | 何をするか | 結果 | 章 |
|---|---|---|---|
| 0 | **変数と約束**: $s$ を独立変数にし、$(x, p_x, y, p_y, z, \delta)$ ($p_{x,y}$ は設計運動量で規格化した正準運動量、先頭が $z > 0$) | SAD の相空間 | [01](01_variables_and_conventions.md) |
| 1 | **ハミルトニアン**: ローレンツ力のラグランジアン → 正準運動量 → ルジャンドル変換 → Frenet–Serret 座標 → $t \to s$ の入れ替え。結果 $H = -\bigl(\sqrt{p^2 - (p_x - a_x)^2 - (p_y - a_y)^2} + a_s\bigr)(1 + x/\rho) + E/(p_0v_0)$ | 要素ごとの $H$ | [01a](01a_hamiltonian_derivation.md) |
| 2 | **記法**: ポアソン括弧、$\exp(:H:)$ (ヘルプの指数写像の読み方) | ヘルプの式を読む道具 | [C](C_poisson_lie.md) |
| 3 | **運動方程式と厳密解**: $H$ が閉じた形で解けるもの (ドリフト、偏向磁石の円弧、四極 + ソレノイドの 2 次の $H_2$、薄い多極のキック) は厳密解。解けないもの ($\Delta H$ など) は運動量だけの関数や座標だけの関数に分け、それぞれ厳密な流れにする | $x \to x + \frac{p_x}{p_z}L$ など | [02](02_drift.md)・[04](04_quad.md)・[05](05_bend.md)・[07](07_thin_multipole.md)・[09](09_solenoid.md) |
| 4 | **分割**: $e^{h(A+B)}$ を $e^{hA}$、$e^{hB}$ の積で近似 (Strang 2 次、Yoshida 4 次、力勾配型 4 次)。各段が正準変換なので、合成もシンプレクティック。刻みの数は要素ごとの経験式 | 本体の写像 | [03](03_symplectic_integrators.md)・[04](04_quad.md)〜[10](10_cavity.md) |
| 5 | **端部**: 磁場の立ち上がり (Maxwell 方程式から縦成分と高次の項) をリー生成関数 $f = \int H_1ds$ にして、母関数による閉じた写像にする | 入口・出口の写像 | [06](06_edge_fringe.md) |
| 6 | **座標の取り替え**: 取り付け誤差・COORD は、剛体運動 + 直線の面への戻し (正準変換) | 要素の前後の写像 | [11](11_misalignment_and_coord.md) |
| 7 | **追跡の写像と転送行列**: 追跡は上の写像そのもの。転送行列 (CALC・Emittance[]) は、**同じ段を閉軌道のまわりで微分して掛けたもの**。$z$ の行はシンプレクティック条件からも組める | 6×6 の $T$ | 各要素の章の「行列と追跡の関係」 |
| 8 | **行列から光学へ**: 閉軌道 (Newton)、分散、x–y 結合の分解、Twiss、拡張 Twiss | $\beta$、$\eta$、$\nu$ | [12](12_linear_optics.md) |
| 9 | **放射と平衡**: kinematical method で放射 (平均損失と量子励起)。$\Sigma = M\Sigma M^\top + \Delta$ でエミッタンス | $\varepsilon$、$\sigma_\delta$、減衰率 | [13](13_radiation.md)・[14](14_emittance.md) |
| 10 | **スピン**: Thomas-BMT → 1 スライスの回転 → 軌道との結合 $D$ → SLIM・Sokolov-Ternov | 偏極時間、偏極度 | [15](15_spin.md) |
| 11 | **特殊な要素**: アンジュレータ (UND)、ビーム・ビーム (BEAMBEAM) | 写像 | [16](16_undulator_beambeam.md) |
| 12 | **集団効果** (空間電荷、ウェイク、CSR、Touschek・IBS) と**非線形への橋渡し** (1 周の写像、正規形、ADTS・RDT) | | [17](17_collective.md)・[18](18_nonlinear_bridge.md) |
| 13 | **Fortran と Python**: 各式が K64 のどのファイルのどのサブルーチンの、どの変数か。pysad のどの関数か | 各章の対応表 | 各章の「Fortran の対応」「pysad の対応」 |

**一枚の流れ図** (矢印は「から」「へ」の変換。左が物理、右が実装):

- 物理: $H$ (01・01a) → [解ける部分 = 厳密解 / 解けない部分 = 厳密な流れに分割] (02〜10) → 端部の母関数 (06) → 座標の取り替え (11)
- 実装 (K64): 要素の写像 `t<要素>.f` (追跡) ↔ その微分 `t<要素>e.f` (行列) → リングの 1 周 `tturn` / `tturne` → 光学 `qtwiss` ・平衡 `tecalc`
- 実装 (pysad): 追跡 `track6.py`・`track6n*.py`・`elkern.py` (numba) ・行列 `lattice._transport_6x6`・`emit6d.py`・`calc6d.py` ・光学 `optics.py`

## 2. 3 系統の実装

同じ要素の物理が、K64 では 3 通りに実装されている。「同じ $H$ の同じ分割」を使うが、別のコードなので、**次数が違う箇所**がある (結果が変わるものだけを各章に書いた)。

| 系統 | 何をするか | 呼び出しの入口 | 要素ルーチン (例) | 配列 |
|---|---|---|---|---|
| **追跡** (`t<要素>.f`) | 粒子の座標を 1 要素ぶん進める。非線形、粒子ごと。`TrackParticles`、ダイナミックアパーチャ、周波数マップ | `tftrack` → `tturn` → `tturn1` (`tturn.f`) | `tdrift`、`tquad` (本体は `tsolqu`、`tsolqur`)、`tbend` (`tbendcore`、K1 ≠ 0 は `tbendi`)、`tthin`、`tmulti` (ANGLE ≠ 0 は `tmulta`)、`tcav`、`ttcav`、`tsol`、`tcoord`、`tchge` (座標の取り替え)、`trad`/`tradk` (放射) | `x, px, y, py, z, g, dv, sx, sy, sz` (粒子の配列) |
| **行列** (`t<要素>e.f`) | 閉軌道のまわりの 6×6 の行列 (`trans(6,12)`)、閉軌道 (`cod`)、ビーム行列 (`beam(42)`)、スピン (`srot(3,9)`) を、1 要素ぶん更新。`Emittance[]`・`CALC6D`・RADCOD | `temit` → `tcod` → `tturne` → `tturne1` (`tturne.f`) | `tdrife`、`tquade` (本体は `tsolque`)、`tbende` (`tbendebody`)、`tthine`、`tmulte`、`tcave`、`ttcave`、`tsole`、`tcoorde`、`tchge`、`trade`/`tradke` | `trans(6,12)` (1:6 が $T$、7:12 が放射の微分)、`cod(6)`、`beam(42)`、`srot(3,9)` |
| **解析微分** (`q*`) | 4 次元の光学 (`CALC`): `qtwiss1` が各要素の `q<要素>` (`qquad`、`qthin`、`qmult`、`qcav`...) を呼ぶ。これらは**同じ `t*e` を呼んで 4×5 に写し取る** (`qquad`: `tquade` → `qcopymat`) 薄いラッパー。さらに FIT の勾配のため、パラメータに関する Twiss の**解析微分** `qd*` (`qdquad`、`qdthin`、`qdbend`、`qddrif`、`qdtwis`) | `tffscalc` → `qcell1` → `qtwiss` (`qtwiss.f`)、`qdtwis` | `qdrift`、`qquad`、`qthin`、`qmult`、`qcav`、`qcoord`、`qchg`、`qdquad`... | `trans(4,5)`、`twiss` |

- **追跡 ↔ 行列の関係**: 設計軌道の近くでは、追跡のヤコビアンと CALC の行列は $10^{-12}$〜$10^{-13}$ で一致する (四極: [04](04_quad.md) §10、偏向: [05](05_bend.md) §9、ドリフト: [02](02_drift.md) §13)。軌道がずれると、分割の次数や刻みの数の違いだけ食い違う (薄い多極: 追跡 4 次・CALC 2 次、[07](07_thin_multipole.md) §7。四極の刻み: 追跡 0.2、CALC 0.1)。
- **追跡の内部変数**: `tquad` などは、入口で $q = p_x/(1+\delta)$ に換えて計算し、出口で戻す (`tconvm`)。行列は正準運動量のまま ([01](01_variables_and_conventions.md) §5)。
- **pysad の対応**: 追跡は `track6.py` (Python、遅い、放射・スピン用) と、numba の `track6n*.py`・`elkern.py`・`perfkern*.py`・`turnkern.py` (高速、スクリーニング向け)。行列は `lattice._transport_6x6`・`emit6d.py`・`calc6d.py`・`solemit.py`・`solmap.py`、光学は `optics.py`・`optics_fn.py`。FIT の微分は `match6.py`・`match6d.py`。pysad の**要素クラス (`elements.py`) の `transfer_matrix()` は別系統**で、K1 などは**長さあたり**の値 (SAD の積分強度 K1 との違いは [01](01_variables_and_conventions.md) §6)。

## 3. 要素ルーチンの命名則

K64 のサブルーチンの名前は、頭の文字と接尾の文字で系統と役割が決まる (主な例)。

| 名前 | 意味 | 例 |
|---|---|---|
| `t<要素>` | 追跡 (粒子の配列)。`t` は tracking | `tdrift`、`tquad`、`tbend`、`tthin`、`tmulti`、`tcav`、`tsol` |
| `t<要素>e` | 行列・閉軌道・ビーム行列 (`trans`、`cod`、`beam`、`srot`)。`e` は emittance | `tdrife`、`tquade`、`tbende`、`tthine`、`tmulte`、`tcave`、`tsole` |
| `t<要素>i` | 複合機能 (K1 ≠ 0 の偏向) | `tbendi` |
| `t<要素>a`・`tmultae` | ANGLE ≠ 0 の MULT | `tmulta`、`tmulte` |
| `t<要素>r` | 放射つきの本体 | `tsolqur` (`tsolqu` の放射つき) |
| `tsolqu`・`tsolque` | 四極 + ソレノイドの本体 (追跡・行列)。四極も `tquad` から `tsolqu` を呼ぶ | |
| `tradk`・`tradk1`・`tradkn`・`tradke` | 放射 (kinematical method): `tradk` 追跡 (1 粒子/配列)、`tradke` 行列 | [13](13_radiation.md) |
| `q<要素>` | CALC の 4×5 の行列 (`t*e` を呼んで写す) | `qquad`、`qthin`、`qmult` |
| `qd<要素>` | CALC の解析微分 (FIT の勾配) | `qdquad`、`qdthin`、`qdbend` |
| `tbfrie`・`tqfrie`・`tqlfre` | 端部 (偏向の Maxwell フリンジ、四極のフリンジ、線形フリンジ) | [06](06_edge_fringe.md) |
| `tcoord`・`tchge`・`tbrot` | 座標の取り替え (COORD、ずれ・ROTATE)、偏向磁石のロール | [11](11_misalignment_and_coord.md) |
| `tturn`・`tturne` | リングの 1 周 (追跡・行列) | |
| `tecalc` | 平衡エミッタンス (`beam` と `trans` から) | [14](14_emittance.md) |
| `tcod` | 閉軌道 (Newton 法) | [12](12_linear_optics.md) |
| `spnorm`・`srequpol` | スピン | [15](15_spin.md) |

(規則から外れるものもある: `tdrife1` (`toplvl.f`) は行列の 1 つの更新、`tradke` は放射の行列。名前は上の表の「基本形」と考える。)

## 4. 各章と読み方

| 章 | 内容 | ヘルプでの $H$・式の所在 | K64 の Fortran | pysad |
|---|---|---|---|---|
| [01](01_variables_and_conventions.md) | 変数と約束 | dynamics (Lagrangean・Hamiltonian。誤植あり) | `tconvm`、`tfstk.f` | `track6.py` |
| [01a](01a_hamiltonian_derivation.md) | ハミルトニアンの導出 | 同上 | — | — |
| [02](02_drift.md) | ドリフト | DRIFT の transformation | `tdrift.f`、`tdrife.f`、`toplvl.f` (`tdrife1`) | `elkern.py` |
| [03](03_symplectic_integrators.md) | シンプレクティックな分割 | THIN (OCT〜DODECA) の transformation | `tthin`、`tsolqu`、`tbendi` | `quadbody.py` |
| [04](04_quad.md) | 四極 | QUAD の transformation | `tquad`、`tsolqu`、`tquade`、`tsolque`、`tqfrie` | `quadbody.py` |
| [05](05_bend.md) | 偏向磁石 | BEND の transformation | `tbend`、`tbendi`、`tbende`、`tbfrie` | `lattice.py` |
| [06](06_edge_fringe.md) | 端部 | FRINGE・DISFRIN、各要素 | `tbfrie`、`tqfrie`、`ttfrin`、`tsfrin` | |
| [07](07_thin_multipole.md) | 薄い多極 | SEXT〜DODECA の transformation | `tthin`、`tthine`、`ttfrin` | `elkern.py` |
| [08](08_mult.md) | 多極 MULT | MULT、`multipole_with_nonzero_ANGLE` | `tmulti`、`tmulte`、`tmulta` | `multa.py`、`multae.py` |
| [09](09_solenoid.md) | ソレノイド | solution-H2・solution-dH | `tsol`、`tsolqu`、`tzparams` | `solmap.py`、`solemit.py` |
| [10](10_cavity.md) | 加速空洞 | CAVI (ΔE の式、$H_f$) | `tcav`、`tcave`、`ttcav` | `tcav.py` |
| [11](11_misalignment_and_coord.md) | 取り付け誤差・COORD | COORD、misalignments | `tcoord`、`tcoorde`、`tchge`、`tbrot` | `lattice.py` |
| [12](12_linear_optics.md) | 行列から光学へ | x-y-coupling、extended-Twiss | `tcod`、`qmdiag`、`qtwiss` | `optics.py`、`lattice.py` |
| [13](13_radiation.md) | 放射 | synchrotron-radiation | `tradk`、`tradke`、`drndsr` | `track6n.py`、`emit6d.py` |
| [14](14_emittance.md) | 平衡エミッタンス | equilibrium-beam-envelope | `temit`、`tecalc`、`tmulbs` | `emit6d.py` |
| [15](15_spin.md) | スピン・偏極 | POL・RADPOL (式なし) | `sprot`、`tradke`、`spnorm`、`srequpol` | `spin.py` |
| [16](16_undulator_beambeam.md) | UND・BEAMBEAM | **節なし** | `undulator.f`、`bbstrhl.f` | `und.py`、`beambeam.py` |
| [17](17_collective.md) | 集団効果 (空間電荷・ウェイクを導出、CSR・Touschek・IBS は対応表) | WSPAC・LWAKE・TWAKE・INTRA・Touschek (機能の説明だけ。式なし) | `tintrb`(`twspac`)、`txwake`、`tffswake`、`spkick`、`csrtrack` | `spch.py`、`wake.py`、`csr.py`、`touschek.py` |
| [18](18_nonlinear_bridge.md) | 非線形への橋渡し (要約と参照) | — (式なし) | — | `resdrive.py`、`nform.py`、`rdttrack.py`、`adtsdirect.py` |
| [C](C_poisson_lie.md) | ポアソン括弧とリー演算子 | — | — | — |

**読み方の目安**:

- 最短: **00 → 01 → 02 → 04 → 05**。四極の章 ([04](04_quad.md)) が書き方の手本 (結論 → 記号 → $H$ → 導出 → 行列 → 追跡 → 行列と追跡の関係 → Fortran → pysad → 確認 → 未確認)。
- 導出を詳しく: 01 → **01a** (ラグランジアンから $H$) → C (リー演算子)。
- 要素別に拾い読み: 各章の冒頭の「この章で分かること・前提」から入る。
- 放射・エミッタンスに進む: 12 → 13 → 14 → (スピン 15)。
- ヘルプの式を疑うとき: 各章の付録 A (原文で確認した式と、確定した食い違い)。

## 5. 全体に共通する約束

- **記号**: [01](01_variables_and_conventions.md) §8 の規約。SAD の Fortran の変数名はそのまま使い、使うたびに注釈し、章末の「この章で使う Fortran の変数名」に表にする。
- **積分強度**: SAD の K1、K2... は積分強度 $B^{(n)}L/B\rho$ [1/m$^n$]。pysad の要素クラスは長さあたり。
- **ヘルプの原文** (K64 版同梱の `SADHelp.html`): 式は TeX 版と平文版の 2 つ併記で、食い違う箇所がある。各章の付録 A に原文で確認した結果を書いた。
- **K64 と v2 (既定) の違い**: 結果が変わるものだけを現在の仕様として書く。K64 の誤りの疑いは各章の付録 B (別枠)。
- **確認スクリプト**: `scripts/<章>_0n_*.py`、出力 `outputs/*.out.txt`。独立な基準 (DOP853 の積分、複素ステップや中心差分、行列指数、sympy) で確かめる ([README.md](README.md) の表)。

## 6. ガイド全体の確認の一覧

どの章で、どの式を、どの独立な基準で確かめたか。数字は各章の出力から ([README.md](README.md) の表にスクリプトの一覧)。

| 章 | 確かめたこと | 基準 | 結果 |
|---|---|---|---|
| 01a | ラグランジアン → 一般化力・正準運動量 → $H = c\sqrt{(\vec p - e\vec A)^2 + m_0^2c^2} + e\phi$、$F_3$、$t \to s$、2 次の展開 | sympy (`ham_01`〜`ham_06`) | 恒等式・展開の係数が一致。展開の誤差は振幅の 4 次 |
| C | ポアソン括弧、$\exp(:f:)$ の例、合成の順序、BCH | sympy (`ham_07`) | 一致。ヘルプの $\exp(:H:)$ の向きは式どうしの整合で確定 |
| 02 | $H \to$ 厳密写像、$dv$、COD のまわりの 6×6 | DOP853、sympy の偏微分、追跡の数値微分 | $2\times10^{-16}$ (積分)、$6\times10^{-18}$ (sympy)、$4\times10^{-11}$ (数値微分)。最短の例 $R_{56} = L/\gamma_0^2$ |
| 03 | BCH の係数 ($S_2$、$S_4$、力勾配型の $\beta, \gamma$)、K64 の補正キックの式 | 自由代数の BCH (sympy)、1 自由度の系での次数 | 記号的に一致。誤差の次数 2・4、$H$ の誤差は有界 (RK4 は 7.7 倍に増加) |
| 04 | $H_2$ の厳密解、$z$ の更新、$\Delta H$、積分器の次数、行列と追跡、フリンジ、保存量 | 行列指数、DOP853、複素/中心差分、Maxwell の場の積分 | $H_2$ は $10^{-18}$ 台。次数 2・4 (誤差の比 4・16)。CALC と追跡のヤコビアンは設計軌道で $10^{-13}$、ずれた軌道で K64 $1.8\times10^{-7}$ (CALC)・$7.1\times10^{-7}$ (追跡) |
| 05 | K1 = 0 の円弧と端面、Maxwell フリンジ、F1、K1 ≠ 0 の本体 | 実験室系の幾何、ローレンツ力の積分、$H_2$ の ODE、中心差分 | 円弧は $1.1\times10^{-15}$、Maxwell フリンジの因子 $p^2$ は比 $1.0000$ |
| 06 | 端のリー生成関数 $f = \int H_1ds$ から (4.27)・(7.17)・(5.35)・DISFRIN が出る、母関数と $\exp(:f:)$ の差 ($K^3$)、`tsfrin` | sympy、ローレンツ力の積分 | 4 つとも sympy で一致。ソレノイドの端は厳密な流れ (差 $4\times10^{-19}$) |
| 07 | 薄いキック、4 次の分割 (補正キック)、フィードダウン、非線形フリンジ | $\exp(:f:)$ の閉じた式、厳密な積分、ローレンツ力の積分 | キックは $9\times10^{-19}$。次数 $L^5$ (比 32)。フリンジの位置のずれは $4\times10^{-19}$ |
| 08 | $K_n$・SK$n$ の意味、基準面、ANGLE = 0 の刻み、$A_s$ の $g_{kn}$ | 数値 (ROTATE との等価)、Maxwell の残差 (sympy)、厳密な ODE | 差 $10^{-17}$。2 次の誤差。ヘルプの $g_{kn}$ は誤り (確定) |
| 09 | $H_2$ の $w_1, w_2$、$z$ の更新、$\Delta H$ のラーモア回転、SOL + QUAD + SOL | sympy、行列指数、ODE | $10^{-18}$ (行列指数)、$10^{-11}$ (全体、分割の誤差) |
| 10 | CAVI の写像、非相対論の補正、重み `ws`、端の収束、TCAVI | 追跡、数値、Maxwell から | 差 $2\times10^{-13}$。シンプレクティック $10^{-10}$ |
| 11 | 剛体運動 + 面への戻し、COORD の行列、ROTATE の順序、`tbrot` | 独立な剛体運動の計算、sympy | シンプレクティック $10^{-12}$。DIR ≠ 0 はヘルプと K64 が不一致 |
| 12 | 閉軌道の 2 次収束、分散 $(1-M_4)^{-1}m_6$、`qmdiag`、Twiss | 追跡の差分、sympy | 残差 $9\times10^{-5} \to 2.5\times10^{-16}$、分散 $10^{-6}$、行列 $10^{-16}$ |
| 13 | `anrad`・`cuc`・`cave`・`cuu`、1 スライスの更新、$U_0$ | 放射スペクトルのモーメント、数値積分、放射積分 | $U_0 = 2.3448$ 対 $2.3449$ MeV、分散 $2.046$ 対 $2.055\times10^{-5}$ |
| 14 | $\Sigma = M\Sigma M^\top + \Delta$、$\Delta$ の組立、各モードの $\varepsilon_k$ | 放射積分の式、13 章の励起、正規座標 | $\varepsilon_x$ 2.2%、$\sigma_\delta$ 0.5%、$J_x$ 1.8%、$\varepsilon_k$ の式 0.03% |
| 15 | Thomas + 静止系 = BMT、1 ステップ回転の次数、定数 $c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2 = r_e\bar\lambda_C$、SLIM の $N$、偏極時間 | sympy、DOP853 (らせん軌道)、玩具の 1 周写像、Sokolov-Ternov の式 | BMT 厳密一致 (sympy)、2 次収束 (比 4.00)、定数 $10^{-12}$、$N$ は 1 次、$\tau_{ST} = 123.8153$ s |
| 16 | UND の Maxwell と写像、BE のキックと $\delta$ のキック、6 次元写像のシンプレクティック性 | sympy、DOP853、2 次元ポテンシャルの数値積分 | $\nabla\times\vec B$ $8\times10^{-15}$、1 次収束、キック $10^{-13}$ (絶対)、$\delta$ のキック $10^{-18}$、シンプレクティック $10^{-11}$ |
| 17 | 空間電荷の係数 $1/(\beta^2\gamma^3)$、WSPAC のキック・線形のレンズ、ウェイクの規格化・符号・自己項の 1/2 | sympy、2 次元ポテンシャルの数値積分、ガウス + 指数関数の解析式 | 実験室系の力 $= eE'/\gamma$。WSPAC のキック $5\times10^{-12}$ (絶対)、$\delta$ のキック $3\times10^{-16}$。ウェイクの平均の損失は解析式の 1.007〜1.027 倍 (切片への粗視化) |
| 18 | (参照の章。docs/24 §7 の比較を転記) | 追跡 + NAFF + フィット、1 周の写像のテイラー展開 | 6 極 2 次 $10^{-3}$ 前後、全部入り 0.01〜0.26% (6 極 ×1) |

**確認の方針**: (a) 行列の $M^\top JM = J$ (シンプレクティック性)、(b) 追跡の数値微分と行列の一致、(c) 独立な数値積分 (DOP853 など) との一致、(d) 刻み数・振幅への収束 (次数)、(e) sympy による係数・恒等式、の 5 種類。**gs (K64 の実機) との比較は主役にしない**。

## 7. ガイド全体の未確認・K64 とヘルプの誤りの疑い

(各章の付録 A・B と「未確認・限界」からの集約。**確定** = 独立な基準 (積分・恒等式・原文) で確認、**疑い** = ソースの読みや一部の数値だけで、SAD (gs) での確定が要る。)

### 7.1 ヘルプの誤りの疑い

**確定 (ヘルプの原文の式を導き直して食い違いを確認)**:

| 章 | 箇所 | 内容 |
|---|---|---|
| 01 | Lagrangean・Hamiltonian | $L_t$ の根号内の符号、$H_t$ の $(p_y - a_y)^2$ の符号、母関数 $G$ の $t_0$ の項、$E$ の $c$ の抜け |
| 02 | DRIFT の transformation | $E/v_0$ の符号 (QUAD は $+E/v_0$ で正しい) |
| 04 | QUAD の線形フリンジ | $b$ に $1/p_1$ が無い ($a$ には有る) |
| 04 | QUAD の刻みの数 | ヘルプの式 $1 + \lfloor10\lvert K1L\rvert/\mathrm{EPS}\rfloor$ は CALC の値。追跡は 0.2 で半分 |
| 05 | BEND の本体 | $s \equiv \mathrm{ANGLE}\times\rho_b\ldots$ の $\mathrm{ANGLE}\times$ は誤記 |
| 05 | 非線形フリンジの $p_{y2}$ | $p_1$ が 1 つ多い |
| 05 | 線形フリンジの $V$ | TeX 版と平文版で式が違う (平文版がヘルプ自身の他の式と一致) |
| 05 | K1 ≠ 0 の分割数 | ヘルプの式とコード (`tbendi.f`) が違う |
| 05 | $\psi_2$ の定義 | AE2 が無い |
| 07 | SEXT〜DODECA の説明の文面 | 平文版は「thin dipole kick」(SEXT)、「thin quadrupole kick」(DECA) など (TeX 版は正しい) |
| 08 | `multipole_with_nonzero_ANGLE` | $g_{kn}$ の $(2k-1)!!$ は $(2k+1)!!$、$z$ の冪は $n + k + 1$ |
| 09 | solution-dH | $\Delta x$、$\Delta y$ の式に $1/B_z$ が無い (ODE で確認、ヘルプのままだと $1.6\times10^{-6}$ ずれ) |
| 11 | APERT (平文版) | `Max[DY2, DY2]` は誤植 (TeX 版は $\max(\mathrm{DY1}, \mathrm{DY2})$) |

**疑い (誤植とみられるが、他の解釈の余地が残る)**:

| 章 | 箇所 | 内容 |
|---|---|---|
| 04 | solution-H2 の $z$ (複素正規座標) | 4 行 3 列の $i$ が誤植か。原因は未特定 |
| 05 | BEND の DISFRIN | 「nonlinear Maxwellian fringe is suppressed」だが、K1 = 0 の追跡は $y^3$ の項だけを止め、Maxwell の硬い端は止めない (仕様の疑い) |
| 09 | solution-dH の $\Delta H$ | `a_c` は $a_x$ の誤記か |
| 10 | CAVI の $H_f$ | $\sin$ の中の位相の書き方 (章の付録 A) |
| 11 | COORD の DIR ≠ 0 | ヘルプの行列と K64 が一致しない (どちらが誤りかは未決) |
| 12 | DETR の閾値 | 同じ文で $\det r \ge 1$ と $0.625$ の 2 つが出る |
| 07 | SEXT〜DODECA の並び | 原文は非対称で 2 次。対称形は 4 次以上 (`thin_04_help_order.out.txt`) |

**ヘルプに無い**: DRIFT の `DISKIN`、TCAVI、UND、BEAMBEAM、スピンの式 (`POL`・`RADPOL` の文だけ)、`Polarization2`・`4`・`6` は K64 のソースに無い。ヘルプの記述の無い式はすべてソースと標準の理論から起こした。

### 7.2 K64 の誤りの疑い

**確定 (独立な式・数値で確認。番号は [docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md))**:

| 章 | 内容 | v2 |
|---|---|---|
| 04 | 追跡の四極の線形フリンジが F1 の符号を捨てる (CALC と不一致、B323) | 直した |
| 08 | MULT (ANGLE ≠ 0) の FB1・FB2 による $L'$ の補正が二重 (`tffs.f:1368` と `tmulta.f:82-90`)。gs の BEND と MULT が $x$ で $2.8\times10^{-7}$ 違い、L を先に伸ばした MULT は BEND と差 0 (B331) | 直した |
| 09 | SOL 領域内の QUAD の $\Delta H$ の運動量の回転が更新後の $p_x$ を使う (`tsolqu.f:122-123`。行列式 $\cos^2 w\ell$、gs のヤコビアンの欠陥 $p_x$ = 0.1 で $3.4\times10^{-6}$、B333) | 直した |
| 10 | CAVI の端 (`tcavfrin`) がエネルギー変化のあと $q = p_x/p$ を再スケールしない (端付き空洞が非シンプレクティック。gs のヤコビアンの欠陥 0.6 MeV/c で $8.1\times10^{-5}$、B332) | 直した |
| 10 | CAVI の端の収束・エネルギー変化に $1/\beta$ が無い (gs と軸の近くの Maxwell の場の数値積分で確認。0.6 MeV/c で収束が 24% 不足、B330) | 直した |
| 05 | K1 > 0 の本体の `xsxkx` の符号 (B250)。$z$ の誤差 $7.3\times10^{-5}$ m (10 m の偏向) | 直した |
| 07 | CALC の `tthine` に非線形フリンジも放射も無い (B291、B240) | 直した |
| 15 | `Polarization` が NaN (`spdepol`、B55)。桁あふれを直しても約 1000 倍の脱偏極率 (B125)。符号つきの `sdamp` (B116)。刻みの角が弦 (B83)。`tbrote` などの行の逐次更新 (B56)、`bsi2` の脱落 (B58) | 直した |
| 16 | UND の磁場が設定できない (B37)、`dpds` の二重 (B78)、BEAMBEAM の傾きの読み位置 (B34、$\sigma_x < \sigma_y$ で $\Delta p$ が誤る) | 直した |
| 15 | `tradke` の `dbtz = dxpzb/pr - dblz` (`photons.f:1550`) が $x$・$y$ 成分 ($\mathrm{/pr^2}$) と不整合 (追跡のスピン回転の差分との差が、δ = 0.05 の偏向磁石で 13% 対 1.2%。gs では確認できない、B336) | 直した |
| 16 | ビームストラールングの光子エネルギーの表の累積が矩形則 + 裾を捨てて規格化 (gs の光子 422 個のエネルギーが再実装と 1.6e-7 で一致、表のモーメントが理論値に対し $+1.60\%$・$+2.62\%$、B335) | 直した |
| 16 | UND の $z$ の更新に $-dv\,L$ が無い ($R_{56} = L/\gamma^2$ の欠落。gs で UND − DRIFT の $z$ が 3 GeV/c・$\delta = 10^{-3}$ で $2.9\times10^{-11}$ m、10 MeV/c で $2.6\times10^{-6}$ m、B334) | 直した |
| 17 | WSPAC の最後のスライスで `gammab(nlat+1) = 0` を読みレンズが 8 倍 (B127)、SPAC の力が $1/\gamma$ 倍 (B89) など (docs/09・docs/10・docs/30 の記録。本章では係数の導出で $1/(\beta^2\gamma^3)$ を確認したのみ) | 直した |

**疑い (ソースの読み・一部の数値。確定には gs での確認が要る。pysad は直していない)**:

| 章 | 内容 | 大きさ・状態 |
|---|---|---|
| 05・06 | CALC の `tbedgemaxwell` の COD の更新 (`tbfrie.f:209`・`:211`) が $\Delta x$ で $p^1$、$\Delta z$ で $p^0$ (正しくは $p^2$、$p^1$)。解析行列とも矛盾 | $\delta y^2/(2\rho)$。$y = 10$ mm・$\delta = 0.05$・$\rho = 10$ m で $2\times10^{-8}$ |
| 17 | WSPAC の縦の力の電位の基準が軸上 (壁の半径に依る縦の空間電荷の主要項が入らない) | 軸上の縦の力が 0。落ちる定数は $b/\sigma$ = 30 で軸から $\sigma$ までの電位差の 15 倍。v2 は `WSPACRADIUS` で基準を壁にできる。**仕様か誤りか未決** (自由空間のモデルでは基準が定義できない。gs の確認は台本が動かず未達) |
| 12 | (pysad 側) 偏向磁石のあるリングで `lat.chromaticity()` と追跡の $d\nu/d\delta$ が合わない | $-1.456$ 対 $-0.780$ など。K64 の誤りではなく pysad 独自の計算 |

**取り下げ (K64 の誤りではない)**: COORD の DIR ≠ 0 の行列が、ヘルプの式と積の順序だけ逆 (11 章)。gs で K64 の形が真の逆変換を符号の反転で表せることを確認。ヘルプの誤記の可能性が高い (コードは変えない)。

**K64 の誤りではないもの (仕様・定義の違い)**: 四極の母関数 $G$ は $\exp(:f:)$ の 2 次までの打ち切り (06)。偏向の Maxwell フリンジは $p_y$ の補正を持たない (06)。K1 の暗黙の定義 (B160、05)。追跡と CALC の分割の次数の違い (04・07)。`dpz` の名前が追跡と CALC で別の量 (02)。

### 7.3 各章の未確認 (主なもの)

| 章 | 未確認 |
|---|---|
| 01 | 曲率 $1/\rho \ne 0$ の $H$ の展開とソレノイドの $a_x, a_y$ は各要素の章。加速のある要素 ($p_0(s)$ の変化) の $t_0(s)$ の扱い |
| 01a | 静的な場の仮定 (時間変化する場・空間電荷・放射の反作用は扱わない)、ねじれ $\tau = 0$ と進行方向の磁場なしの座標系、3 次以上の展開の影響 (§9.6 の数値例のみ) |
| 02 | `DISKIN` (K64 を gs で実行して比べていない。pysad は未実装)、`tdrife` と `tdrife0` の違い、$p_z \to 0$ の限界 (`ampmax`) |
| 03 | BCH 級数の収束、$\tilde H$ の「指数的に長い時間」の保存 (定理の引用、数値は 5 万刻みまで)、刻みの数の経験式の根拠 ($\varepsilon$ = 0.2、0.1)、多自由度での 5 次の完全な分解 |
| 04 | 非線形フリンジの K の 2 次の項と実際の端の高次、F1 < 0 の k64 の食い違い (v2 は直した)、ヘルプの $z$ の更新の食い違いの原因、**pysad の `emit6d._jac_tqente` の $z$ の行 (5,2)・(5,4) の欠け** ($2\times10^{-9}$、pysad 側) |
| 05 | Maxwell フリンジの $y^4$ の係数 (独立な導出なし)、$L'$ の補正の導出、K1 ≠ 0 の `tbendicorr`・端面、**F1 ≠ 0 のとき CALC と追跡のヤコビアンが $10^{-7}$ の桁で違う (原因未特定)** |
| 06 | 偏向の $y^4$ の係数、$p_y$ 依存の補正 $f_4$、CAVI の端の導出 ($1/\beta$ は 10 章で確認済み)、MULT の端 (`ttfrins`・`tblfri`)、ソレノイドの $p_x, y, p_y$ の 3〜20% の確認。SAD (gs) との比較はしていない |
| 07 | 補正キックの確認は $k = 2, 3$ (DECA・DODECA・四極の `kord == 1` は一般形の構造を確認していない)、$\gamma$ の項が誤差を減らすとは限らないこと、$K$ の 2 次のフリンジ、ROTATE・DROTATE・DZ との組合せ |
| 08 | CHI1・CHI2・DZ、端面の角の式 (8.9) の導出、刻み数の式の係数 6 の根拠、複数の多極・SK$n$・K1 ≠ 0 の ANGLE ≠ 0 の ODE との比較、$L'$ の二重適用 (確認済み。FRINGE = 1・2 と K1 ≠ 0 の BEND との食い違いは未調査) |
| 09 | `tztafs` の係数の個別の導出、`tsolrot` の数値確認、放射 (`tsolqur`)・スピン・CALC の解析微分、四極のフリンジと SOL の端の合成 |
| 10 | TRPT の写像・重みの選び方・`ndiv` の式、端の収束の $1/\beta$ (確認済み。長い空洞の本体との組合せ・定在波の位相の進みは未確認)、端のエネルギー変化の符号と係数 ($1/\beta$ の因子は確認済み)、`offset` の物理的な意味、V1・V11・V02 の CALC の係数 |
| 11 | `tsolrot` (CHI1・CHI2・DZ)、偏向磁石の回転 (11.9) の代数的な導出、DIR ≠ 0 のヘルプの誤記の確認 (設計者へ)、ソレノイド中の変位による $p_x$ のゲージ、CALC の `tchge` の全成分とスピン、端のフリンジが回転した軸に対して持つ構造 |
| 12 | `qmdiag` の特別な場合、RADCOD、`tecalc` の $V$、偏向ありの色収差 (pysad 側の不一致) |
| 13 | RADCOD・RADTAPER、`tr1` の導出、分割の定数 |
| 14 | `tnorm`・`btr` の導出、IBS、SynchroBeta、結合のあるリング・MINCOUP の数値 |
| 15 | Sokolov-Ternov と DK の QED の導出、`spdepol` の級数、高次のスピン軌道共鳴、運動論的偏極 |
| 16 | 交差角の変換 (16.12) の導出、BE の閉じた形の導出、強-強、スライス数の収束、UND の端の処理 |
| 17 | SPAC・PSPAC・CSR・Haissinski・IBS・Touschek の導出と数値確認 (対応表のみ)、WSPAC の Laslett 式との一致の取り直し、ウェイクの光学 (`WAKEOPT`)・複数バンチ |
| 18 | 参照の章。式の導出・数値は docs/24・25・26 |

## 8. 付録・ほかの文書との関係

- **付録 C** ([C_poisson_lie.md](C_poisson_lie.md)): ポアソン括弧とリー演算子。ヘルプの $\exp(:H:)$ の読み方。
- **付録 A** (ヘルプの式と K64 の食い違いの総覧): 各章の付録 A に書いた。総覧は予定 (§7.1 が現状の集約)。
- **付録 B** (確認スクリプトの一覧): [README.md](README.md) の表。
- **ほかの文書**: K64 の式の早見表は [docs/23_physics_guide.md](../docs/23_physics_guide.md)、K64 のバグは [docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md)、v2 と k64 の違いは [docs/28_v2_vs_k64.md](../docs/28_v2_vs_k64.md)、使い方は [docs/01_overview.md](../docs/01_overview.md)、スピンの使い方は [docs/16_spin_pol.md](../docs/16_spin_pol.md)。

## 9. 限界

この章は地図であり、式の導出は各章にある。§6・§7 の数字と判定は、2026 年の各章の付録・出力からの転記で、各章の「使った条件」(エネルギー、磁場、周長) が前提。§7.2 の「疑い」は、確定ではない。
