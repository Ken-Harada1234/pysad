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
| 02_drift.md | ドリフトと $dv$、COD のまわりの線形化 | 予定 |
| 03_symplectic_integrators.md | 分割の一般論 (Strang、Yoshida、力勾配型)、刻みの数の一覧 | 予定 |
| [04_quad.md](04_quad.md) | **四極**。$H_2$ の厳密解、運動学的な $\Delta H$、Strang と Yoshida、刻みの数、フリンジ、行列と追跡の関係 | ✓ |
| [05_bend.md](05_bend.md) | **偏向磁石**。$H_2$ の厳密解、K1 = 0 の円弧と端面角の厳密解、Maxwell フリンジ・線形フリンジ (F1)、複合機能 (K1 ≠ 0)、行列と追跡の関係 | ✓ |
| 06_edge_fringe.md | 端部の一般論 (Maxwell、F1・F2、非線形フリンジ) | 予定 (四極の端部は 04 章 §7) |
| 07_thin_multipole.md | 薄い 6 極・8 極など。4 次の分割 | 予定 |
| 08_mult.md | MULT (多極) | 予定 |
| 09_solenoid.md | ソレノイド。$H_2$ の 4×4 の厳密解、記号の対応 | 予定 |
| 10_cavity.md | 加速空洞 (CAVI)・TCAVI | 予定 |
| 11_misalignment_and_coord.md | 取り付け誤差と COORD (座標の正準変換) | 予定 |
| 12_linear_optics.md | 行列から光学へ (閉軌道、周期解、結合、Twiss) | 予定 |
| 13_radiation.md | 放射 (kinematical method、光子数、励起と減衰) | 予定 |
| 14_emittance.md | 平衡エミッタンス | 予定 |
| 15_spin.md | スピンと偏極 | 予定 |
| 16_undulator_beambeam.md | アンジュレータ、ビーム・ビーム | 予定 |
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
py -3.12 physics_guide/scripts/bend_01_h2_body.py               # 05 章 (偏向磁石) の確認 (bend_01〜bend_04)
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
| `bend_01_h2_body.py` | 05 | $H_2$ の展開・運動方程式・厳密解、$z$ の更新 (K64 の `xsxkx` の符号)、分散・$R_{56}$ |
| `bend_02_arc_faces.py` | 05 | K1 = 0 の円弧と端面角の厳密解 (実験室系の幾何・曲線座標の $H$ の ODE との比較)、端面の行列 |
| `bend_03_fringe.py` | 05 | 端の磁場、Maxwell フリンジ、F1 (線形フリンジ)、CALC の `tbedgemaxwell` のシンプレクティック性 |
| `bend_04_matrix_vs_map.py` | 05 | CALC の行列と追跡の写像の関係、段ごとの比較、端面の 2×2 行列 |
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
