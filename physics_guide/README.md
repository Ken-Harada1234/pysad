# 物理ガイド: ハミルトニアンから転送行列・追跡、Fortran・Python まで

SAD のヘルプには、各要素の転送や計算の元になるハミルトニアンが書いてある。しかし、ハミルトニアンから実際の転送行列・追跡の写像にどう至るのか、その途中は書かれていない。この「物理ガイド」は、その隙間を埋め、**ハミルトニアン → 運動方程式 → 厳密解または分割 → 追跡の写像と転送行列 → K64 の Fortran → pysad の Python** を、要素ごとに 1 章ずつ説明する。

- **読み手**: 加速器物理は知っているが、SAD の内部は知らない人。
- **各章の約束**: 結論を先に書く。記号と単位を表にする (SAD・K64 の Fortran・pysad の 3 列)。導出は途中の式を省かず、式に番号を付ける。Fortran と Python の対応表は、式番号で引ける。確認スクリプトと出力を付ける。未確認・限界・使った条件 (エネルギー、磁場、周長など) を書く。
- **物理が主役**: ヘルプの式を導き直し、K64 の実装との違い (癖・バグの可能性) は、各章の「Fortran の対応」「限界」に分けて書く。v2 (既定) と k64 (K64 と同じ動作) の違いは、結果が変わるものだけを現在の仕様として書く。
- ほかの文書との分担: 使い方と SAD の機能の仕様は `docs/`、物理の導出はこのガイド。K64 の式の早見表は [docs/23_physics_guide.md](../docs/23_physics_guide.md)、K64 のバグは [docs/30_sad_bug_report.md](../docs/30_sad_bug_report.md)、v2 と k64 の違いは [docs/28_v2_vs_k64.md](../docs/28_v2_vs_k64.md)。

## 目次

✓ は書けているもの。「予定」は、まだ書いていないもの。

| 章 | 内容 | 状態 |
|---|---|---|
| [01_variables_and_conventions.md](01_variables_and_conventions.md) | 変数と約束。正準変数、ラグランジアンから母関数まで、$dv$、$p_x/(1+\delta)$、K1 の約束 (四極の章に必要な範囲) | ✓ (四極に必要な範囲。ほかは予定) |
| [01a_hamiltonian_derivation.md](01a_hamiltonian_derivation.md) | **ハミルトニアンの導出**。ラグランジアン、速度依存のポテンシャルと一般化力・一般化運動量、ルジャンドル変換、正準変換、Frenet–Serret 座標、$t \to s$ の入れ替え、2 次の展開と偏向・四極の運動方程式 | ✓ |
| 00_overview.md | 全体像。ハミルトニアンから Fortran・Python までの流れ、3 系統の実装 (追跡・行列・解析微分) | 予定 |
| [02_drift.md](02_drift.md) | **ドリフト**。$H$ から $x' = p_x/p_z$、$z$ の更新と $dv$ の項、厳密 = 近軸 × $D$、近軸 (DISKIN)、COD のまわりの 6×6 (`tdrife1`)、最短の例 ($R_{56} = L/\gamma_0^2$) | ✓ |
| [03_symplectic_integrators.md](03_symplectic_integrators.md) | **シンプレクティックな分割の一般論**。BCH と Strang (2 次)、Yoshida (4 次)、後方誤差解析、力勾配型の係数 $\alpha, \beta, \gamma$ の導出、K64 の刻みの数の式の一覧 | ✓ |
| [04_quad.md](04_quad.md) | **四極**。$H_2$ の厳密解、運動学的な $\Delta H$、Strang と Yoshida、刻みの数、フリンジ、行列と追跡の関係 | ✓ |
| [05_bend.md](05_bend.md) | **偏向磁石**。$H_2$ の厳密解、K1 = 0 の円弧と端面角の厳密解、Maxwell フリンジ・線形フリンジ (F1)、複合機能 (K1 ≠ 0)、行列と追跡の関係 | ✓ |
| [06_edge_fringe.md](06_edge_fringe.md) | **端部 (フリンジ)**。Maxwell 方程式から縦磁場と高次の項、ベクトルポテンシャルと $H$ の展開、階段の端のリー生成関数 $f$ (偏向・四極・多極・ソレノイドが同じ手順)、母関数による閉じた写像 (リー変換の 2 次までの打ち切り)、偏向の Maxwell フリンジの因子 $p$ の再確認 (CALC の疑い)、CAVI・SOL・MULT の端、FRINGE・DISFRIN の表 | ✓ |
| [07_thin_multipole.md](07_thin_multipole.md) | **薄い多極 (SEXT・OCT・DECA・DODECA)**。ポテンシャルと積分強度、薄いキックの厳密解 ($\exp(:f:)$ の級数が止まる)、4 次の力勾配型 (`tthin`) の補正キックの導出、CALC の 2 次の分割 (`tthine`)、フィードダウン、非線形フリンジ (Maxwell 方程式から) | ✓ |
| [08_mult.md](08_mult.md) | **多極 MULT**。ポテンシャルの展開と K・SK の意味、基準面 ($\theta_1$)、ANGLE = 0 の分割 (厳密な四極 + 双極の本体と薄いキックの対称な 2 次分割)、刻みの数、フリンジ、ずれ、ANGLE ≠ 0 の $A_s$ を Maxwell 方程式から導く (ヘルプの $g_{kn}$ との食い違い) | ✓ (CHI1・CHI2・DZ、RF・放射・ソレノイドつきは読みだけ・範囲外) |
| [09_solenoid.md](09_solenoid.md) | **ソレノイド**。$H_2$ の厳密解 ($w_1$、$w_2$、固有モード)、$z$ の更新、$\Delta H$ のラーモア回転、Strang 分割、端部 (DISFRIN)、ヘルプ ↔ `tzparams` の対応表 | ✓ (領域・座標の切替 `tsolrot` は読みだけ) |
| [10_cavity.md](10_cavity.md) | **加速空洞 (CAVI)・TCAVI**。$\Delta E$ から $\delta$・$z$・$p_x$ の更新 (時間を保つ)、同期位相と `offset`、横方向の依存と Panofsky–Wenzel、非相対論の補正 ($I_0$)、長い空洞の分割 (重み `ws`)、$p_0(s)$ の変化、端の効果 $H_f$ (Maxwell からの導出)、TCAVI | ✓ (TRPT の写像・CALC は読みだけ) |
| [11_misalignment_and_coord.md](11_misalignment_and_coord.md) | **取り付け誤差と COORD**。剛体運動 + 直線の面への戻し (スナップショットと $z$)、`tcoord`・`tcoorde` の 6×6、COORD の DIR (K64 とヘルプの食い違いの疑い)、DX・DY・ROTATE の順序と符号、偏向磁石のロール (`tbrot`)、MARK・APERT・MAP は力学でないこと | ✓ (CHI1・CHI2・DZ の `tsolrot` は読みだけ) |
| [12_linear_optics.md](12_linear_optics.md) | **行列から光学へ**。閉軌道の Newton 法 (`tcod`、2 次収束)、周期解と分散 $(1-M_4)^{-1}m_6$、x–y 結合の分解 (`qmdiag`: ヘルプの $R$ と K64 の閉じた式)、Twiss と伝搬、色収差、拡張 Twiss の $H$ | ✓ (`qmdiag` の特別な場合・RADCOD・`tecalc` の $V$ は読みだけ。偏向ありの色収差は未決) |
| [13_radiation.md](13_radiation.md) | **放射 (kinematical method)**。曲げ角 $\phi$ と $\rho_r$ (ヘルプ)、$N_\gamma$・$u_c$・$\langle u\rangle$・$\langle u^2\rangle$ の古典の式と `anrad`・`cuc`・`cave`・`cuu`、追跡の 1 スライス (`tradk1`: 運動量の縮み、$-\tfrac12dp_x\,dg$ と $al/3$ の導出)、FLUC のポアソン過程と `drndsr`、CALC の励起 (`tradke` の `dee`)、放射積分 | ✓ (RADCOD・RADTAPER・`tr1` の導出・分割の定数は未確認) |
| [14_emittance.md](14_emittance.md) | **平衡エミッタンス**。$\Sigma = M\Sigma M^\top + \Delta$ の導出、$\Delta$ の組み立て (`tmulbs`)、`tecalc` (固有値・正規座標・減衰率・分配数)、各モードの $\varepsilon_k = \bar\Delta_k/(1-e^{-2a_k})$、同期運動と RFSW、MINCOUP | ✓ (`tnorm`・`btr` の導出・IBS・SynchroBeta は未確認) |
| [15_spin.md](15_spin.md) | **スピンと偏極**。Thomas-BMT (静止系の歳差 + Thomas 歳差、縦・横の係数 $(1+a)$・$(1+a\gamma)$)、1 スライスの回転 `sprot` (Rodrigues、刻みの 2 次、弦の誤差)、軌道との結合 $D$ (`tradke` の積算)、`spnorm` (偏極ベクトル・スピンチューン)、SLIM の $N$ の式と脱偏極率、Sokolov-Ternov の定数 (`pst`・`cfpd`・`sflc` が `anrad`・`cuc`・`cave`・`cuu` から出る)、DK との対応、平衡偏極度 | ✓ (ST・DK の QED の導出、`spdepol` の級数、高次共鳴は未確認) |
| [16_undulator_beambeam.md](16_undulator_beambeam.md) | **アンジュレータ (UND) とビーム・ビーム (BEAMBEAM)**。ヘルプに節がないのでソースから。UND: $A_x, A_y$ の模型 (Maxwell 方程式)、傍軸の $H$、陰的な 1 次のシンプレクティック写像 (DOP853 と比較)。BEAMBEAM: Bassetti-Erskine のキック (ポテンシャルの数値積分と一致)、$\delta$ のキック、スライスの 6 次元写像 ($s^* = (z-z_s)/2$ への drift と戻し、$F = s^*H_0$ の流れの共役)、強ビームの包絡とスライス、交差角・傾き、ビームストラールング | ✓ (交差角の変換 (16.12) の導出は未確認。UND の $-dv\,L$ の欠落は疑い) |
| 17_collective.md | 集団効果 (空間電荷、ウェイク、CSR、Touschek・IBS) | 予定 |
| 18_nonlinear_bridge.md | 非線形への橋渡し (正規形、ADTS・RDT) | 予定 |
| [C_poisson_lie.md](C_poisson_lie.md) | 付録: ポアソン括弧とリー演算子 $\exp(:f:)$ (ドリフト、四極、6 極の例、合成の順序、SAD のヘルプの記法との対応) | ✓ |
| A_help_vs_source.md | 付録: ヘルプの式と K64 のソースの食い違い・誤植の一覧 | 予定 (01 章と 04 章の付録 A に先に書いてある) |
| B_script_index.md | 付録: 確認スクリプトの一覧 | 予定 (下の表) |

読む順の目安: **01 → 01a (導出を詳しく追うとき) → 04 → 05**。ヘルプの $\exp(:H:)$ の記法は付録 C。**01 → 04** (四極を例に、ガイドの書き方を掴む)。そのあと、要素別に拾い読みしてよい。

## 共通の記号の表 (最小限)

詳しくは各章の「記号と単位」の表。

| 量 | 意味 | 単位・約束 |
|---|---|---|
| $(x, p_x, y, p_y, z, \delta)$ | SAD の正準変数 | $x, y, z$ は m。$p_x, p_y$ は $p_0$ で規格化した正準運動量。$\delta = p - 1$。先頭が $z > 0$ |
| $p = 1 + \delta$ | 全運動量の規格化値 (設計運動量 $p_0$ に対して) | 無次元 |
| $q_x = p_x/(1+\delta)$ | 追跡の内部変数 (角度に近い) | 追跡 (`tquad` など) だけ。行列 (CALC) は正準運動量 |
| $dv = 1 - v/v_0$ | 速度の差 | $z' = -dv$ (01 章 §4) |
| K1、K2 … | 積分強度 ($K1 = B'L/B\rho$ など) | 1/m、1/m² …。K1 > 0 は水平収束 |
| $\rho$、$\rho_b$ | 曲率半径 | 偏向磁石の章で定義 |
| cod | 閉軌道 | 正準運動量の 6 成分 |
| `trans(6,12)` | CALC・Emittance[] の行列。1:6 が $T$、7:12 が励起 $D$ (放射の章) | |

## ヘルプの原文

各章が「ヘルプ」と言うのは、K64 版に同梱の `Documents/SADHelp.html` (画像は同じフォルダの `SADHelp_img`。CERN 版と同系統) の原文で、式・係数・記号を読んで照合してある。ヘルプの式は、TeX 版 (`EquAlign` の 1 つめ) と平文版 (2 つめ) が併記されていて、**両者が食い違う箇所がある** (BEND の線形フリンジ $V$、本体の $p_{x2}$、SEXT〜DODECA の説明文)。各章の付録 A に、原文の式を引いて「原文で確認」「確定した食い違い」「疑い」に分けて書いた。ヘルプには、ポアソン括弧・$\exp(:H:)$ の向き・積の順序の定義が無い。向きは式どうしの整合から確定した ([C_poisson_lie.md](C_poisson_lie.md) §8.2。$\exp(:H:)$ は $H$ (長さをかけた量) を時間 1 のハミルトニアンとする流れ)。

## 確認スクリプトの実行方法

**前提**: 四極の線形フリンジの F1 の符号の修正 (v2 の追跡が CALC と一致する) は v2.4.4 に入る。この物理ガイド (dev/v2.5.0) は v2.4.4 の後に載せる前提で書いてあり、`quad_04_matrix_vs_map.py` の (f) の v2 の結果は、その修正の入った pysad で得た。この作業ツリーの pysad がまだ修正前のときは、環境変数 `PYSAD_SRC` に v2.4.4 以降の pysad のルートを指定して実行する (`scripts/_common.py` が読む)。

スクリプトは `physics_guide/scripts/` にあり、リポジトリのルートから実行する (Python 3.12。pysad は同じフォルダのものを使うので、インストールや `PYTHONPATH` は要らない。numpy・scipy・matplotlib が要る):

```
py -3.12 physics_guide/scripts/quad_01_h2_exact.py
py -3.12 physics_guide/scripts/quad_02_dh_drift.py
py -3.12 physics_guide/scripts/quad_03_integrator_order.py      # 図 physics_guide/figures/quad_integrator_order.png も作る
py -3.12 physics_guide/scripts/quad_04_matrix_vs_map.py
py -3.12 physics_guide/scripts/quad_05_fringe.py
py -3.12 physics_guide/scripts/quad_06_invariants.py
py -3.12 physics_guide/scripts/drift_01_h_to_map.py             # 02 章 (ドリフト) の確認 (drift_01〜drift_03)
py -3.12 physics_guide/scripts/integ_01_bch_force_gradient.py   # 03 章 (分割の一般論) の確認 (integ_01、integ_02)
py -3.12 physics_guide/scripts/thin_01_potential_kick.py        # 07 章 (薄い多極) の確認 (thin_01〜thin_03)
py -3.12 physics_guide/scripts/bend_01_h2_body.py               # 05 章 (偏向磁石) の確認 (bend_01〜bend_04)
py -3.12 physics_guide/scripts/edge_01_gauge_generators.py       # 06 章 (端部) の確認 (edge_01〜edge_04)
py -3.12 physics_guide/scripts/sol_01_h2_exact.py               # 09 章 (ソレノイド) の確認 (sol_01、sol_02)
py -3.12 physics_guide/scripts/drift_04_help_convention.py     # ヘルプ原文の確認 (drift_04、quad_07、bend_05、thin_04。各章の付録 A)
py -3.12 physics_guide/scripts/mult_01_as_expansion.py          # 08 章 (MULT) の確認 (mult_01、mult_02)
py -3.12 physics_guide/scripts/cav_01_energy_map.py            # 10 章 (空洞) の確認
py -3.12 physics_guide/scripts/mis_01_coord_geometry.py        # 11 章 (取り付け誤差・COORD) の確認 (mis_01、mis_02)
py -3.12 physics_guide/scripts/opt_01_linear_optics.py          # 12 章 (行列から光学へ) の確認
py -3.12 physics_guide/scripts/rad_01_kinematical.py            # 13 章 (放射) の確認
py -3.12 physics_guide/scripts/emit_01_equilibrium.py           # 14 章 (平衡エミッタンス) の確認
py -3.12 physics_guide/scripts/spin_01_thomas_bmt.py            # 15 章 (スピン) の確認
py -3.12 physics_guide/scripts/und_01_field_map.py              # 16 章 (UND・BEAMBEAM) の確認 (und_01、bb_01。bb_01 は -W ignore を付けると scipy の警告が出ない)
py -3.12 -W ignore physics_guide/scripts/bb_01_kick_slice.py
py -3.12 physics_guide/scripts/ham_01_lorentz_potential.py      # 01a、付録 C の sympy の確認 (ham_01〜ham_07)
```

結果は標準出力と `physics_guide/outputs/<スクリプト名>.out.txt` に出る (UTF-8)。本文に引用した数字は、この出力のもの。スクリプトは、どれも数十秒以内。`physics_guide/scripts/_common.py` が共通部 (リポジトリのルートを `sys.path` に足す、出力の保存、日本語フォントの設定、四極 1 台のラティスの作成)。

| スクリプト | 章 | 確認すること |
|---|---|---|
| `quad_01_h2_exact.py` | 04 | $H_2$ の厳密解 (行列、$z$ の更新) と pysad の本体・`transfer_matrix` との一致 |
| `quad_02_dh_drift.py` | 04 | $\Delta H$ の解 (運動学的ドリフト) |
| `quad_03_integrator_order.py` | 04 | Strang (2 次)・Yoshida (4 次) の誤差の次数、刻みの数の表、図 |
| `quad_04_matrix_vs_map.py` | 04 | CALC の行列と追跡の写像の微分の一致、$z$ の行、シンプレクティック性 |
| `quad_05_fringe.py` | 04 | フリンジ (母関数、Maxwell の場の端の積分、F1・F2、シンプレクティック性) |
| `quad_06_invariants.py` | 04 | 保存量と分割の誤差 |
| `quad_07_help_checks.py` | 04 | ヘルプ原文: 線形フリンジの $b$ の $1/p$、solution-H2 の $z$ の式、F1 の符号 |
| `drift_01_h_to_map.py` | 02 | $H$ の偏微分・$z'$ の整理・$dv$ (sympy)、閉じた式と DOP853・pysad の追跡の比較、厳密 = 近軸 × $D$、DISKIN (pysad は未実装) |
| `drift_02_linearization.py` | 02 | 6×6 の偏微分 (sympy) と `tdrife1` の式、pysad の行列と追跡の数値微分、シンプレクティック性と $z$ の行 |
| `drift_03_r56_example.py` | 02 | 最短の例: ドリフトを 1〜4 回通った後の $z$ (3 GeV・10 MeV)、4 つのドリフトの $T_{56}$ |
| `drift_04_help_convention.py` | 02、付録 C | ヘルプの $\exp(:H:)$ の向きの確定、DRIFT の $E/v_0$ の符号の食い違い (sympy) |
| `integ_01_bch_force_gradient.py` | 03 | 自由代数の BCH (sympy): Strang・Yoshida の 3 次の係数、力勾配型の $t$・$c_3$・$c_5$ ($\alpha, \beta, \gamma$)、K64 の `tthin` の補正キックと記号的に一致、`tthine` の分割 |
| `integ_02_order_symplectic.py` | 03 | 1 自由度の系で、Strang・Yoshida・力勾配型・RK4 の誤差の次数、シンプレクティック性、長時間の $H$ の誤差 |
| `thin_01_potential_kick.py` | 07 | 多極のポテンシャル (調和・規格化)、$\exp(:f:)$ の級数が 2 項で止まる、pysad の薄い多極と閉じた式、フィードダウン (DX・DY) |
| `thin_02_fourth_order.py` | 07 | 長さのある多極の `tthin` (4 次) を厳密な積分と比較: 誤差の次数 (L を 2 倍にした比)、$\gamma$ の項の効果、CALC の `tthine` (2 次) とのヤコビアンの違い |
| `thin_03_fringe.py` | 07 | 多極 ($m = 2$〜6) の非線形フリンジの位置のずれ (Maxwell から導いた式と `ttfrin`)、運動量・$z$・シンプレクティック性、6 極・8 極のローレンツ力の数値積分との比較 |
| `thin_04_help_order.py` | 03、07 | ヘルプ原文: SEXT〜DODECA の並び (原文は非対称で 2 次、対称形は 4 次以上)、$V_n$ の符号 |
| `mult_01_as_expansion.py` | 08 | ANGLE ≠ 0 の $A_s$: ヘルプと K64 の $g_{kn}$、Maxwell 方程式の残差 (sympy) |
| `mult_02_tracking.py` | 08 | MULT の追跡と厳密な ODE の比較 (`ndiv` の式・2 次の誤差、SK$n$ と ROTATE、DX・DY、ANGLE ≠ 0 の K2) |
| `cav_01_energy_map.py` | 10 | CAVI の薄い写像 (追跡との一致)、シンプレクティック性、非相対論の補正と $I_0$、重み `ws`、TCAVI、端の収束の $1/\beta$ |
| `mis_01_coord_geometry.py` | 11 | ヘルプの回転行列と `tcoord.f` の $r_{ij}$ (DIR ≠ 0 は不一致)、独立な剛体運動 + 面への戻しと追跡・CALC の一致、6×6 の解析式と数値微分、シンプレクティック性、DIR の非可逆 |
| `mis_02_element_misalign.py` | 11 | QUAD の DX・DY・ROTATE の順序と符号、BEND の `tbrot` (DROTATE・CHI2) と支点まわりの剛体回転、直線の MULT の CHI1・CHI2・DZ が追跡で無視されること、MARK・APERT |
| `opt_01_linear_optics.py` | 12 | 閉軌道の Newton 法の収束、分散 $(1-M_4)^{-1}m_6$ と追跡の差分、$BAB^{-1}=P(\mu)$ と Twiss の伝搬 (sympy)、`qmdiag` の式で $RMR^{-1}$ がブロック対角、色収差、ヘルプの $H$ のシンプレクティック性 |
| `rad_01_kinematical.py` | 13 | 定数 `anrad`・`cuc`・`cave`・`cuu`、放射スペクトルのモーメント、`drndsr` のモーメント、偏向磁石 1 台の平均・分散、ずれた四極の損失、$-\tfrac12dp_x\,dg$ と $al/3$ の数値積分、$U_0 = C_\gamma E^4I_2/2\pi$ |
| `emit_01_equilibrium.py` | 14 | `calc6d` (`tecalc`) の結果と放射積分の式 ($U_0$・$J_x$・$\varepsilon_x$・$\sigma_\delta$・$\alpha_p$・$\nu_s$)、$\Delta = \Sigma - M\Sigma M^\top$ と 13 章の励起、正規座標の $\Sigma_n$、各モードの $\varepsilon_k$ の式 |
| `spin_01_thomas_bmt.py` | 15 | Thomas + 静止系 = BMT (sympy)、らせん軌道の BMT の DOP853 と 1 ステップ回転の積 (2 次)、平面の偏向 $a\gamma\theta$・縦磁場 $(1+a)\psi$、K64 の弦、定数 $c_{\rm ave}a_{\rm nrad}c_{\rm uc}^2 = r_e\bar\lambda_C$ と DK の 11/18、$NM = R_0N + [n_0]_\times D$、平面リングの SpinTune・PolarizationTime・Polarization |
| `und_01_field_map.py` | 16 | $\nabla\times\vec B = 0$ (sympy)・軸上の磁場、スライスの写像と DOP853 (1 次)、最大角 $F_0/K_z$、UND と DRIFT の $z$ の差 $= dv\,L$ |
| `bb_01_kick_slice.py` | 16 | Bassetti-Erskine のキックと $\delta$ のキック (ポテンシャルの数値積分、3 通りの楕円)、k64 の傾きの誤り (B34)、6 次元写像のシンプレクティック性・往復の恒等、ビームストラールングの定数と光子エネルギーの表 |
| `bend_01_h2_body.py` | 05 | $H_2$ の展開・運動方程式・厳密解、$z$ の更新 (K64 の `xsxkx` の符号)、分散・$R_{56}$ |
| `bend_02_arc_faces.py` | 05 | K1 = 0 の円弧と端面角の厳密解 (実験室系の幾何・曲線座標の $H$ の ODE との比較)、端面の行列 |
| `bend_03_fringe.py` | 05 | 端の磁場、Maxwell フリンジ、F1 (線形フリンジ)、CALC の `tbedgemaxwell` のシンプレクティック性 |
| `bend_04_matrix_vs_map.py` | 05 | CALC の行列と追跡の写像の関係、段ごとの比較、端面の 2×2 行列 |
| `bend_05_help_checks.py` | 05 | ヘルプ原文: 線形フリンジ $V$ の TeX 版・平文版、本体の $s$・$p_{x2}$、非線形フリンジの $p_{y2}$ |
| `edge_01_gauge_generators.py` | 06 | $\nabla\times\vec A = \vec b$ (偏向・四極・$2m$ 極・ソレノイド)、端の生成関数 $f = \int H_1ds$ から (4.27)・(7.17)・(5.35)・DISFRIN の $H$ が出る (sympy)、偏向の $y^4$ の項 |
| `edge_02_generating_function.py` | 06 | 母関数 $G$ = 点変換 $x + V + \frac12V\cdot\nabla V$ (sympy)、厳密な $\exp(:f:)$ との差 ($K^3$)、$z$ の更新 ($\partial G/\partial p$) |
| `edge_03_bend_maxwell_recheck.py` | 06 | 偏向の Maxwell フリンジの因子 $p$ (ローレンツ力の積分、CALC の COD の更新と解析行列の比較) |
| `edge_04_solenoid_fringe.py` | 06 | `tsfrin` = $f$ の厳密な流れ、Maxwell の端の場のローレンツ力の積分との比較 |
| `sol_01_h2_exact.py` | 09 | 特性方程式と $w_1$・$w_2$ (sympy)、ヘルプの式・K64 (`_TZ`) と行列指数の一致、$z$ の更新、$\Delta H$ の解 (ヘルプの $1/B_z$ の欠落)、固有モード |
| `sol_02_chain.py` | 09 | SOL + QUAD + SOL の追跡と、端の流れ + 全 $H$ の ODE の基準の比較、DISFRIN |
| `ham_01_lorentz_potential.py` | 01a | ベクトル恒等式、$dA/dt$、一般化力 (追加項)、$\vec p = \gamma m_0\vec v + e\vec A$、ラグランジュ方程式 |
| `ham_02_legendre.py` | 01a | ルジャンドル変換、$H = c\sqrt{(\vec p - e\vec A)^2 + m_0^2c^2} + e\phi$ |
| `ham_03_frenet_f3.py` | 01a | 母関数 $F_3$ の $\vec P$、ポアソン括弧、$A_s$ と磁場 |
| `ham_04_t_to_s.py` | 01a | $t \to s$: $P_t = -H$、$\tilde H = -p_s$、陰関数の微分 |
| `ham_05_expansion.py` | 01a | 2 次の展開、偏向磁石・四極の運動方程式、展開の誤差の振幅依存性 |
| `ham_06_series.py` | 01a、04 | 展開の途中式 (series)、偏微分の 1 項ずつの表、$\Delta H$ の展開 |
| `ham_07_lie_operator.py` | 付録 C | ポアソン括弧、$\exp(:f:)$ の例、合成の順序、BCH、線形写像 |

**環境の注意**: 四極の章の検算は数値 (numpy・scipy) で行う。導出の章 (01a と付録 C) の検算は sympy (`ham_*.py`。`py -3.12` に sympy 1.14 が入っていること。数秒) で行う。numba のキャッシュが古いときは `py -3.12 -m pysad.nbcache --force`。台本の実行 (`run_sad_string`) を使うので、SAD の Packages は要らない。

## 数式の書き方と HTML の注意

- 数式は md の中に TeX で書く: インラインは `$x_1$`、独立した行は `$$ ... $$`、式番号は `\tag{4.7}`。文字の `$` は `\$`。インラインの `$` の内側の前後は空白にしない (閉じる `$` の直後に数字を置かない)。コードブロックと、バッククォートの中の `$` は数式にならない。
- HTML は `py -3.12 tools/build_html.py --include-untracked physics_guide` で作る (md と同じ名前の `.html`)。**このリポジトリの `tools/build_html.py` は数式 (TeX) を変換しない** (`tools/README.md` の「限界」)。数式は TeX のまま HTML に入るので、数式として読むときは md を、TeX を描画できるビューアで読む。KaTeX による描画 (Node.js が要る) の変換器は、この作業ツリーの `tools/` にはない。
- 図は `physics_guide/figures/` に置き、md から相対パスで参照する (300 KB 以下の png。日本語のラベルは、スクリプトの `_common.font()` の設定で描く)。
- 文書間のリンクは `[04_quad.md](04_quad.md)` のような相対リンク。docs への参照は、`docs/<名前>` の形で統一する (docs の番号が変わったときの一括置換のため)。

## 書き方の約束 (章を足すときの手本)

四極の章 ([04_quad.md](04_quad.md)) の構成: 冒頭に「この章で分かること・前提・読む順」→ 1 結論 (数字つき) → 2 記号と単位 (SAD・K64・pysad の 3 列) → 3 ハミルトニアン → 導出 (4〜7) → 行列 (8) → 追跡の写像 (9) → 行列と追跡の関係の表 (10) → 放射などほかの章への参照 (11) → Fortran の対応表 (12。ファイル:行、サブルーチン、変数) → pysad の対応表 (13。関数、対応する Fortran、添字・符号・単位の注意) → 確認スクリプトと出力 (14) → 未確認・限界・使った条件 (15) → 付録 A (ヘルプとの違い)。

- 数式の番号は章ごとに (章番号.通し番号)。本文・Fortran の対応表・pysad の対応表の列に式番号を書く。
- 導出で、本質でない代数は確認スクリプトに回し、スクリプトの出力 (誤差の数字) を本文に引用する。独立な基準 (DOP853 の数値積分、複素ステップ微分、行列指数) で確かめる。
- ヘルプの誤植・省略は、本文では正しい形で書き、付録 A に集める。
- 履歴・版番号は書かない。K64 と v2 の違いは、結果が変わるものだけ。
