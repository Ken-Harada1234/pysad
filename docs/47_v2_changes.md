# 47. pysad v2 の変更の一覧 (v2.0.0・v2.1.0・v2.2.0・v2.3.0。v1 = K64 に忠実な翻訳との違い)

最終更新: 2026-11-22 (v2.4.1 開発中の DRAW・ampnf の限界の拡張を追記)。節は新しい版が先: v2.4.1 (開発中)、v2.3.0 (下の (a)〜(e))、v2.2.0 ((a)〜(f))、v2.1.0、v2.0.0。

- v1 系は SAD (K64 版) のソースに忠実な翻訳で、K64 のバグや癖も gs の実動作に合わせて再現した。v2 系は、K64 のバグの修正、物理の誤りの修正、機能の追加を行う ([`../CHANGELOG.md`](../CHANGELOG.md))。
- この文書は一覧。根拠の数値・条件・使い方の詳細は作業の記録 (v2.0.0 は `history/v2_changes_*.md`、v2.1.0・v2.2.0 は `history/v2next_*.md`、v2.3.0 は `history/v23_*.md`) にある。記録と例の対応は文書の末尾の「記録と例」。
- K64 の側のバグの番号 (B01〜B162) は [`21_sad_bug_report.md`](21_sad_bug_report.md)。B160〜B162 は v2.3.0 (B160 は T1 で見つけ U1 で「K1 の暗黙の定義」に再分類、B161 は U1 の疑いを U2 が確定して v2 で直した、B162 は T2 の Touschek の式の適用範囲)。B72〜B110 は v2.0.0、B111〜B134 は v2.1.0、B135〜B159 は v2.2.0 の作業で新しく見つけたもの。
- 使った条件: 電子。3 GeV の FODO リング (24 セル 105.6 m / 91.2 m、12 セル 45.6 m / 62.4 m)、1 GeV・3 GeV の短い輸送路、0.5 GeV・2 GeV のリング、ソレノイドは 0.3〜3 T/(Bρ) 相当、
  RADLIGHT は 10.2 MeV、空間電荷は 5〜50 MeV。個別の条件は各表と記録に書いた。**実機の大きなリング・挿入光源・強い結合では確かめていない。**
- v2 の値の根拠は、解析式、独立な数値計算 (scipy の solve_ivp によるローレンツ力・Thomas-BMT の積分など)、保存則、別のコードとの一致で、gs ではない。
  gs (SAD 1.4.4.2k64、2026-10-19 ビルド) は「K64 の動作の確認」と「K64 が正しく動く機能の参照」に使った。

## v2.4.1 (開発中)

- **6 極の最適化の項の指定 (`pysad/sxtopt.py` の拡張。新機能。互換モードによらず同じ動作)**: 項ごとに `objective`/`constraint`/`ignore`・重み・目標値・制約の幅を指定して K2L (ファミリー) を計算する `SxtSpec`・`SxtSpecObjective`・`design_k2l` (1 回の SLSQP。反復・探索はしない)。項は 3 次の RDT 5 項・δ の 1 次の駆動項 3 項・4 次の個別の項 (`h31000` など)・ADTS 3 成分・ξ・誤差感度 (`err_h21000` など)・グループ (`w3`・`chrom`・`w4`・`adts`・`xi`・`err3`)。辞書・文字列・JSON・規則の列で指定でき、SAD の台本のオプションから同じ形で渡せる。既定 `legacy_spec()` は従来の `SxtObjective` と同じ値・勾配 (1e-12)。作動点の周りの共鳴の一覧 `near_resonances` (次数・距離・駆動項。b = l−m が奇数の線は水平面対称の磁石では何次でも駆動されない)・警告 `resonance_warnings`。例 [`example/ddba/v07_14/sextupole_v2/`](../example/ddba/v07_14/sextupole_v2/README.md) (v07_14、5 GeV、周長 749.48 m、20 ファミリー、|K2L| ≤ 32.65 1/m²): K2L のランダム誤差 3%・10% の 5 シードの DA で評価し、gs の `TrackParticles` で確認。**ADTS だけの最適化は完成版の x-y DA の 0.62 倍 (誤差なし)・0.64 倍 (3%)・0.81 倍 (10%)、解析式の項を全部小さくしても届かない (gs、5 シードの平均)。完成版のそばの 1 項ずつの走査で、δ の駆動項 h00201 を 0 にした候補 r6b が完成版の 1.30 倍 (誤差なし)・1.17 倍 (3%)、10% では同じ** (r6a は 3% で 1.09 倍で揺らぎの範囲)。誤差で生じる項はファミリーの強さでは 3% 以下しか変わらない。LMA (RF 6.5 MV・放射あり、4 点、誤差なし) は h00201 を 0 にしても有意には下がらなかった (完成版 3.17%、r6a 3.19%、r6b 3.06%)。未確認: 誤差つき LMA・Touschek 寿命・位置ずれ・補正。試験 `test_v241_sxtopt_spec.py` (10 件)。
- **四極の本体の積分器の共通化と v2 の 4 次の合成 (B226・B243。v2 のみ。k64 は従来とビット単位で同じ)**: K64 の `tsolqu`・`tsolque` (bz = 0) は、非線形ドリフト D と厳密な線形の本体 B を D(h/2) B(h) D(h/2) の Strang 分割 (2 次) で積分し、刻みの数は軌道のずれに依らない。pysad では同じ構造が 11 か所に重複していたので、`pysad/quadbody.py` に集約した (刻みの群 `group()`、刻みの数 `ndiv()`・`ndiv_calc()`、追跡の D・B `make_track()`。追跡の python 版・numba 版・MULT・放射つき、CALC の python 版・numba 版 (`_quad_full`・`k_quad`)・解析ヤコビアン・MULT の `tsolque_pieces`・放射の刻み)。v2 は **4 次の Yoshida の合成** (3 本体/刻み) にし、刻みの基準を追跡 ε = 0.2 → 0.1、CALC ε = 0.1 → 0.05 にした。numba のカーネルへは `rs[4]` (追跡)・`PH[13]` (CALC) で方式を渡す。
- 結果が変わるもの (v2): 3 mm・2 mrad ずれた四極 (L = 0.3、K1L = −0.6) のヤコビアンの厳密な写像との差 CALC 1.8e-7 → 5e-10、追跡 7.1e-7 → 7.9e-9。CALC・Emittance[] の四極の行列は軌道がずれたとき 1e-7〜1e-6 変わる (ステアリングの閉軌道を含むリングの CALC のチューンが 3e-10〜2.5e-7 動く)。追跡の座標は 1e-9〜2e-8 (10 mm・数 mrad の軌道)。numba の四極 1 台の評価は 2〜2.5 倍。ADTS の精密版 (`fourth='r6'`) の `kinematic_h4(quad_rule='k64')` は v2 では厳密な積分と同じ (台形則の誤差が無いため)。
- 限界・未確認: 追跡と CALC のチューンの差はステアリング K0 = ±0.006 で +2.3e-6 → −6.3e-8 (0.012 で +9.2e-6 → −2.5e-7) に減るが 0 にならない (符号が変わる。原因は未確認)。リング全体の追跡時間・DA の結果の変化は測っていない。EPS で刻みを細かくすると 4 次で減る (k·aln² > 1 の粗い刻みでは Yoshida の方が悪いので EPS を大きくしない)。SOL 内の四極 (bz ≠ 0)・MULT ANGLE ≠ 0・偏向磁石の本体は別の積分器で対象外。`pysad/tracking.py` の `track_thick_quad` は厳密な行列 + z だけの別の簡易版で、変えていない。試験 `test_v241_quad_integrator.py` (独立な基準は運動方程式と変分方程式の DOP853 積分)。B227 (放射の刻みの既定) と B238 の残り (低エネルギーの加速 + POL) は確認のみで、直していない (docs/46)。
- **周波数マップの例の表の再作成 (B226・B243 の結果としての更新。例の表のみで、コードの誤りではない)**: v2.3.0 で作った周波数マップの例の表 (`example/v2_3/Q1_fmap_astrid2.table.txt`・`Q1_fmap_kekls.table.txt`。41 × 31、1024 周、RFSW・RAD OFF) は旧い 2 次の四極の積分器の値で、v2 の既定が 4 次の合成になったため、同じ点を追跡し直すとチューンが 1e-7 以上ずれて試験 `test_v23_Q1_fmap_lma.py::test_example_frequency_map_tables` が落ちた (`quadbody.set_scheme(0)` なら通る)。ずれの原因は四極の積分の誤りの減少: ASTRID2 の表の全点で、チューンの差は中央値 5e-5〜7e-5、95% 点 1.5e-4〜1.6e-4 (一部のカオス的な点は最大 0.02〜0.04)。ASTRID2 の 4 点 (x = 2〜10 mm・y = 2〜6 mm、256 周) で、K64 の方式の誤りは最大 5e-5、四極の刻みを 4 倍・16 倍細かくすると 3.6e-6・2.3e-7 に減り、4 次の合成の結果 (刻み 4 倍細かい 4 次との差 1.2e-6) に収束することを確認した (旧い値の方が粗い)。そこで表は v2.4.1 の既定で再作成して新しい名前 `Q1_fmap_astrid2_v241.{table.txt,out.txt,png}`・`Q1_fmap_kekls_v241.{…}` で保存した (`Q1_fmap_ring.py <名前> <workers> _v241`。条件は同じ。v2.3.0 の表・図は上書きせず残してあり、docs/48・50・S1 の例・A5 は従来のまま v2.3.0 の表を使う) 。試験の参照は `_v241` の表に向けた (許容は変えていない: 追跡し直しの一致 1e-8)。再作成の表の結果: ASTRID2 は失われた点 0 のまま (D の範囲 −11.37〜−1.44、旧 −11.45〜−1.70)。KEK-LS は失われた点 460 → 465 (落ちる集合が 9 点違う。DA の境界の点)、次数 4 以下の共鳴線の近くで D > −5 の割合は近い 0.195・遠い 0.082 (共鳴線に載る傾向は同じ)。もう 1 件 `test_v23_S1_fmapres.py::test_k64_mode_runs_same_map` は、k64 と v2 の周波数マップが 1e-6 以内で一致するという前提が、四極の積分器の違いで成り立たなくなった (k64 は quadbody の方式 0 と v2 方式 0 とビット単位で同じ。差は 5.4e-5 までの積分の誤り): v2 の積分を方式 0 にした場合と 1e-6 以内で一致することと、v2 の既定との差が 1e-6〜1e-4 であることを確かめる形に直した。
- **MULT のスキュー成分の gs との比較の試験の許容 (B226 の結果としての更新。試験のみで、コードの誤りではない)**: `test_v2_C_sol.py::test_mult_skew_components_match_gs[v2]` が、v2 の MULT (L = 0.5、K1 = 0.1、SK0〜SK3・K2 を持つ) の追跡 `TR` と gs の差 2.7e-12〜1.1e-11 (許容 2e-13) で落ちた。原因は四極の本体の 4 次の合成 (B226) だけ: v2 の積分を `quadbody.set_scheme(0)` にすると gs と 0 (ビット単位) で一致し、B290 の K0K1 BEND・B291 の 6 極フリンジ (この MULT は DISFRIN=1)・B292 の MULT 端の縦収束 (AE = 0) は関与しない。差は K64 の Strang 分割 (2 次) の誤りで、MULT の EPS を 1 → 0.1 → 0.01 と小さくすると 1.1e-11 → 1.5e-12 → 1.6e-13 と刻みの 2 乗に比例して減る (TW の BX・BY の差 2.5e-9・1.2e-8 も同様に EPS に比例して減る。相対 3.3e-9)。MULT 全体の刻みの誤り (EPS = 1 で収束値との差 7e-8〜1.3e-7) よりずっと小さい。試験は、k64 と v2 の方式 0 は従来どおり gs と 2e-13 (TW は 1e-9 の相対)、v2 の既定は TR 2e-11・TW 1e-8 の相対、かつ方式 0 との差が 1e-12 を超える (4 次の合成が効いている) ことを確かめる形にした (許容を数値だけ緩めたのではなく、期待を K64 の積分誤りの量に更新)。本体 (`lattice.py` など) は変更なし。
- **v2 の既定が B226 (四極の 4 次の合成)・B291 (6 極のフリンジの CALC への追加) で動いたことによる試験の期待の更新 (8 件。試験のみで、コードの誤りではない)**: v2.4.1 の変更のあとの全試験 (2506 件) で 8 件が落ちた。いずれも HEAD (v2.4.0) では通る。原因は 2 つで、`quadbody.set_scheme(0)` (四極を K64 の方式にする) と、CALC の 6 極のフリンジを外す (`lattice._thin_fringe_flag` を常に False にする。B291) で従来の値に戻る。評価器 (`sadscript.py` の RuleList など)・FOURIER・光子の変更は関与しない。(1) **B226 だけ (5 件)**: `test_v2_B_ffs.py::test_b30_main_off_rfsw_reaches_ffs` は v2 の追跡を k64 の NORFSW の追跡と全桁一致で比べていた (互換モードをまたいだ比較が、四極の積分器の違いで 1e-9 台ずれた)。v2 の追跡は v2 の NORFSW・RFSW と全桁一致、k64 との差は 1e-7 以内、方式 0 なら k64 と全桁一致、の形にした。`test_v2_C_wakecmd.py::test_b66_old_wake_command_equals_wakefunction`・`test_v2_C2_wakecmd.py` の 2 件 (`test_multiline_table_equals_wakefunction`・`test_statement_after_endwake_on_same_line`) は、旧式の `WAKE` コマンドの追跡を gs の値 (`C3_REF`) と相対 2e-9 で比べていたが、v2 の既定とは 4.4e-8 違う (K64 の 2 次の積分の誤り)。比べる相手を、同じ表を `WakeFunction` で与えた台本の同じ互換モードの結果 (差 0) にし、方式 0 では gs の値と従来の許容で一致することを確かめる形にした (WAKE コマンドの読み方の試験で、積分器の試験ではない)。`test_v2_C_track.py::test_spin_tracking_loss_matches_gs[v2]` は POL + RAD の追跡で、12 mm の粒子 (7 周) の x が gs と 6.5e-7 違い (許容 2e-7)、方式 0 にすると 2e-7 以内に戻る: v2 の x の許容を 2e-6 にし、方式 0 では従来の 2e-7 で gs と一致、かつ既定と方式 0 の差が 1e-7 を超えることを確かめる形にした。(2) **B226 と B291 の両方 (3 件)**: `test_v24_H3_qcod_robust.py::test_v2_robust_qcod_finds_cod_where_k64_fails` は KEK-LS (3 GeV、周長 570.72 m、誤差 1.5 倍で COD が最大 1.9・3.4 mm) の v2 のチューン NY が 17.6923195779 → 17.6923143372 (5.2e-6)、1 倍の NX が 48.5760358865 → 48.5760353019 (5.8e-7)、COD 1.3085727e-3 → 1.3085717e-3 動いた。方式 0 でフリンジありは 17.6923227887 (B291 の分 3e-6)、方式 0 ・フリンジなしで従来の値に戻る (v2 の新しい値を許容 1e-8 で、戻した値が従来の値と 1e-8 で一致することを試験に入れた)。`test_v2_P5_micado_fitp_oy.py::test_fitp_final_solution_matches_gs_k64[five]`・`[beta]` は FITP の最終解を gs と比べる試験で、v2 の最終の残差 (five: 1.8068e-5、gs 1.8079e-5。方式 0 でフリンジあり (B291 だけ) は 1.8089e-5、方式 0 でフリンジなしは k64 と同じ 1.8076e-5) が許容 2e-4 を超えた (6e-4。オフモメンタム δ = ±2% の COD が数 mm〜cm の 6 極で、チューンが 1e-7〜1e-6 動いた結果の合わせ残り)。K1 は 6 桁一致、K2 は 2e-4 以内のまま。許容を、v2 の残差 2e-3 (K1・K2 は従来のまま)、v2 と k64 の最終解が残差 2e-3・K 1e-3 以内、と更新し、`five` では方式 0・フリンジなしの v2 が k64 と全桁同じ最終解になることを確かめる形にした (従来の「FITP は両モードで同じ」の全桁一致は、四極の積分器の違いで成り立たない)。k64 と v2 (方式 0・フリンジなし) は従来の許容で gs と一致し、v2 の既定の期待は原因と大きさを書いて物理的な量 (積分の誤りの減少と、追跡と CALC の整合) に更新した。本体 (`lattice.py` など) は変更なし。

- **FFS の `DRAW` を matplotlib で実装** (新機能。K64 では `$DisplayFunction` が CanvasDrawer (Tk) のときだけ動き、gs は表示環境がないと何も出さない。v1 の pysad は未対応として無視していた)。互換モード (v2・k64) によらず同じ動作。結果が変わる既存の機能はない (DRAW が「無視」から「描く」に変わる。`$DisplayFunction` が Identity で保存先がなければ従来どおり何も出さない。メッセージ `[pysad] FFS command DRAW は未対応 (無視)` は出なくなった)。
- 内容: `pysad/draw_mpl.py` (構文解析は CANVASDRAW.n、点列・単位・ラティス図・InfoLabel は OpticsPlot.n の移植)。保存先は `PYSAD_DRAW_FILE`・`Draw$Option` の `File`・Python API の `save`。`pip install pysad-k64[plot]`。詳細・限界は [docs/42 §7b](42_graphics_guide.md)。
- 試験: `test_v25_draw_mpl.py` (33 件。構文解析、OpticsPlot の Graphics と全点一致、gs の Graphics との一致、保存、`matplotlib` が無いときの警告)。例: `example/11_graphics/fn_DRAW.*`。
- 未解決: SADScript の `OpticsPlot[…]` の Graphics を matplotlib に流す選択肢 (`PYSAD_DISPLAY=mpl` で OpticsPlot を差し替える形) は未実装 (環境変数 `PYSAD_DISPLAY=mpl` は今のところ DRAW を matplotlib にするだけ)。窓での見た目・クリックは未確認。

### ANGLE=0 の K1≠0・K0≠0 の BEND・6 極のフリンジ・MULT の端の縦収束 (B290〜B292)

- **B290 (K64 の誤り。v2・k64 とも直した)**: ANGLE=0・K1≠0・K0≠0 の BEND は、gs で追跡が長さ L のただのドリフト (`tbend` が `tmulti` に rtaper = 0 を渡し、nmmax の判定も逆)、CALC・Emittance[] が SIGSEGV (`tbende` が `tmulte` に引数を 2 個渡さない。B02 の原因) になる。pysad は gs の誤りを再現せず、`tbend` が渡すはずだった値で同じ K0・K1・L の MULT と同じに計算する (追跡・CALC・Emittance[]。F1 は MULT の FB1・FB2 = F1 + FB)。**gs との照合は MULT で行い、12 セル FODO の追跡 (NORAD・RAD)・CALC・Emittance[] が全桁一致** (F1・FB・FRINGE=3 つきも)。numba の追跡 (K_MULT) に対応 (F1・FB・AE が無いとき。numba と Python は 1e-13。以前は常に Python の追跡)。独立な基準: 一様な K0・K1 の運動方程式 (k64 4e-9、v2 6e-10)、v2 の追跡の一周写像のチューンと CALC のチューン 1e-7 以内 (AE・F1・FB・DX・ROTATE つき)。k64 の追跡は端の角 AE を無視する (B142)。
- **B291 (K64 の不整合。v2 で直した。k64 は K64 のまま)**: L≠0 の SEXT・OCT・DECA・DODECA と MULT の n ≥ 2 の端の非線形フリンジ ttfrin (Lee-Whiting 型、DISFRIN=0 が既定) を、追跡は入れるが CALC・Emittance[] は入れない (`tthine`・`tmulte` にフリンジが無い)。座標の 3 次なので理想の軌道では効果 0、軌道がずれると一周写像が追跡とずれる (K2 = ±6 の FODO にステアリング K0 = 0.012 を入れてチューン 1.9e-3、MULT の K3 + ROTATE で 6.8e-3)。v2 は追跡と同じ式を COD の写像と解析ヤコビアン (`elkern._fringe_step`: 実数の式を複素ステップ微分) に入れた (numba の `k_thin` に `fringe`、Python の経路・turnkern・emit6d・MULT の `_mult_sequence`)。1.9e-3 → 4e-5 (フリンジなしと同じ。残りは本体の積分の次数の違い)。**結果が変わるもの (v2)**: 軌道が 6 極の中心を外れた CALC・Emittance[]・LMA・周波数マップの前提の光学 (理想の軌道では変わらない)。`lattice._el_cod_step_cached` のキーに互換モードを加えた。
- **B292 (K64 の薄い蹴りを v2 で厳密にしたものの誤り。v2 で直した。k64 は K64 のまま)**: MULT の硬い端 (AE1・AE2) の縦方向の収束 (docs/46 の「厳密な Maxwell と 2e-5」) の原因は 2 つ。(1) 面を横切る水平の傾き θ への依存が無い (厳密には tan(e ± θ) − tan(±θ)、tan(±θ) は本体の双極の項が持つ)。(2) K1 以上は z B(z) でなく横磁場のポテンシャル ψ (水平方向の運動を厳密な流れが担うので、K1 の項が 2 倍に数えられていた)。独立な基準 (滑らかな階段の Maxwell 場 ∇(ψ f) の中のローレンツ力の数値積分、w の外挿。ANGLE=0 の BEND の端の鎖とも px 1e-11 で一致) との差が py 5e-7 → 1.6e-8、K1 = 0.4 で 4e-7 → 2e-10 (追跡・CALC・δ = 0.01 とも)。
- 実装しない・直さない (理由は docs/46 の「v2.4.1 A〜E」): B239 (θγ ≲ 1 の陽子の端部放射。局所の光子の模型では書けず、基準が無い。目安をスケーリングで 1/(θγ) 倍の桁に見直した)、B238 の残り (RF の向きの変化の分離。効果 8e-6〜1.4e-5、大きな変更)。
- 試験: `test_v241_A_bend_k1_k0.py` (36 件)、`test_v241_B_sext_fringe.py` (14 件。k64 が gs の追跡 (フリンジあり) と CALC (フリンジなし) を再現する確認を含む)、`test_v241_C_mult_edge_ykick.py` (12 件)。回帰 (QUICK=1、k64) は 50 ベンチ FAIL 0 (290 s)。関連の既存試験 (SEXT・MULT・AE・追跡・Emittance[] を含む 87 ファイル) は、test_apert〜test_v23_S1 (約 420 件) のうち 3 件 (`test_v23_Q1_fmap_lma.py::test_example_frequency_map_tables` の 2 件と `test_v23_S1_fmapres.py::test_k64_mode_runs_same_map`) が落ちる。**この 3 件は本項の変更ではなく v2.4.1 の四極の積分器 (B226・B243) の v2 の 4 次の合成による (`quadbody.set_scheme(0)` で 3 件とも通る)。例の表 (v2.3.0 で作成) を再作成 (`Q1_fmap_*_v241.*`)、k64 の試験を直して、3 件とも通るようにした (上の「周波数マップの例の表の再作成」)**。それ以降のファイルは全部は走り切っていない (計算機の負荷が高く)。

### ampnf の限界の拡張 (J の高次・5・6 次の駆動項・1 本の共鳴の島・K0≠0 の磁石のオプション)

- **`pysad/ampnf_ext.py` と `ampnf.py` の追加 (新機能。互換モードによらず同じ動作。既定の呼び出しの結果は変わらない)**: (1) J² は `ampnf` の D ≥ 5 で既に入っていたので (D ≥ 7 で J³)、API だけ足した: `AmpTuneResult.nu(plane, jx, jy, δ, jorder=)` (J の次数で切った級数)、`d2nu(plane)` (∂²ν/∂Ji∂Jj)。5・6 次の共鳴の駆動項は `ampnf_ext.drive_table(res)` (純粋な n = 3〜6 次の共鳴の写像の単項式の係数。`ampnf.normal_form(keep=)` で共鳴の単項式を消さずに残した正規形から)。(2) 1 本の共鳴の島: `resonance_normal_form`・`island_theory` (固定点の作用と振り子の島の半幅)・`track_island` (追跡の周期 n の固定点 O・X と分離線の作用の範囲。検証用)。(3) `k0_nonlinear=True|False` (`amplitude_tune`・`one_turn_series`・`build_ops`): K0≠0 の磁石の非線形を他の素子のオプションと独立に入れる・落とす。`None` (既定) は従来の規則。
- 確認 (docs/49 §15): ν(J) の J の次数ごとの誤差を FODO・ASTRID2・KEK-LS で追跡 (pysad k64、基準の gs の TrackParticles は 2.6e-11 以下で一致) と比べた (J² で J¹ のみの 1/3〜1/20、J に比例して残る誤差は入れていない効果 (K64 の積分器など)、KEK-LS は差の共鳴の上 (J = 6.4e-7) で崩れる)。島は、24 セルのうち 12 セルに 6 極を置いた非周期の FODO の 6νx = 5 で、固定点の作用が追跡と 0.04〜1.5%、島の端が 0.3〜1.3%。5 次・4 次は島が広く振り子の式が外れる (`pendulum_ok` = False)。K0≠0 の磁石 (3 極の挿入光源、K0 = 0.06 rad) 自身の非線形の J¹ への寄与は 0.02%。試験 `test_v241_ampnf_ext.py` (5 件)。例 `example/v2_4/I1_nu_of_J_orders.py`・`I2_island_vs_tracking.py`。
- **誤りの指摘・K64 の修正はない** (B 番号なし)。注意すべき限界: **リングの対称性で禁止された共鳴にも `island_theory` は有限の島を返す** (24 セルの FODO の 4νx = 3: 追跡は島なし、理論は半幅 J_r の 14%)。K0≠0 の磁石の `bend_cubic`・`kin`・`fringe` の個別の選択は閉じた式が分けられないためできない (docs/49 §15.4)。
- 未確認・未解決 (docs/46 の v2.4.1 ampnf の節): 垂直面・和差の共鳴の島の追跡との比較、ASTRID2・KEK-LS の島、`drive_table` の h_jklm との対応、`sxtopt` の目的関数への J² の導入 (`sxtopt.py` は他の作業で編集中のため触っていない)。

### DXF (ラティスの平面図の CAD 図)

- **`pysad/dxf_lattice.py` と SADScript `DrawMagnetsDXF[]` (新機能。K64 にない。互換モードによらず同じ動作)**: ラティスの平面図を ASCII の DXF (R12 相当) に書き出す。基準軌道・偏向磁石・四極などの外形・MARK・要素名をレイヤー別に出す。幾何は `LINE["GEO"]` と同じ値。参考にした SAD の `drawmagvw3d.n` (VectorScript) との違いと、磁石の寸法の既定・限界は [docs/51](51_lattice_dxf.md)。結果が変わる既存の機能はない。
- 確認: gs と pysad の台本の出力が 8 本とも全行一致、Python 関数とも一致、閉じたリングで始点と終点が合う、偏向磁石の外形の独立な幾何の確認。試験 `test_dxf_lattice.py`。例 `example/dxf_lattice/`。
- 追加 (既定の出力は変わらない): `MagnetSpec[MagnetType[要素名]]` の読み込み (`UseMagnetSpec`・`magnet_spec`。優先順は量ごとに `Spec` の要素名 > 型名 > `MagnetSpec` > 既定値) と、要素名の重なり消し (`Overlap->"remove"`・`overlap='remove'`。既定 `"keep"`)。gs・pysad の台本・Python 関数が KEK-LS の標準セル (全周も) で全行一致 ([docs/51](51_lattice_dxf.md) §2b・§2c)。
- 未解決: 立体・コイル・端板は未対応。`MagnetSpec`・`Spec` なしの磁石の寸法は仮の既定値。`MagnetSpec` の書式は呼ぶ側の使い方から読んだもので実際の定義ファイルでは未確認 (`Height` は読むだけ)。重なり消しの外接矩形は見積もり (CAD の字体とは違う)。CAD ソフトで開く確認は未実施 (自前のパーサと matplotlib の描き直しのみ)。

## v2.3.0 (2026-11-21)

v2.2.0 の後に、K64 にない解析機能 (周波数マップ・局所運動量アクセプタンス (LMA)・共鳴駆動項 (RDT)・振幅依存チューンシフト (ADTS)・Touschek の独立な式・数値の正規形) を足し、それぞれの信頼性を調べた版。解説は [`48_frequency_map_lma.md`](48_frequency_map_lma.md)・[`49_resonance_terms_adts.md`](49_resonance_terms_adts.md)・入門 [`50_adts_rdt_primer.md`](50_adts_rdt_primer.md)、記録は `history/v23_*.md` (Q1・Q2・R1〜R6・S1〜S3・T1・T2・U1・U2・L2)、例は [`example/v2_3/`](../example/v2_3/)。**v2.2.0・v2.1.0・v2.0.0 の節は変えずに残す**。この節の (a)〜(e) は v2.3.0 だけの一覧。

- **新機能は K64 に無い追加で、互換モード (v2・k64) に関係なく同じ動作**。numba のカーネルは編集していない。結果が変わる (a) は、新機能の既定の選択と、既存の v2 の 2 つの機能 (`TouschekLifetime[INF, INF, nz]`・CALC の BEND の行列) だけ。
- **解析式の位置づけ**: RDT・ADTS・共鳴線の予想などの解析式は、摂動論の定義された次数と、入れた効果の範囲で原理的に作る。追跡 (全非線形) とは入れている効果が違うので、厳密な一致は目的ではなく、追跡に合わせるための補正項・当てはめは足していない (S3 では、いったん入れた追跡からの当てはめ `fit_tunes` を削除した)。追跡との一致の数値は、入れた効果と次数の検証の結果として docs/49・[`example/v2_3/ADTS_summary.md`](../example/v2_3/ADTS_summary.md) に載せ、差の主因は docs/49 §11.9 にある (原因が分からないものは「未特定」)。追跡から係数を取り出す唯一の部分は偏向磁石の生成関数 (`bendnl.py`) で、原理は「写像の 2 次のテイラー係数が生成関数を一意に決める」(docs/49 §8.1)。
- K64 の B 番号は B160 (K1 の暗黙の定義。誤りではない)・B161 (CALC の `tbendecorr` の dk。v2 で直した)・B162 (Touschek の損失率の式の適用範囲。v2 は Piwinski 型を既定にした)。
- 使ったリング: ASTRID2 ([`example/astrid/lattice/ASTRID2_v01_01.sad`](../example/astrid/lattice/ASTRID2_v01_01.sad)、0.58 GeV、45.7 m)、大きなリング A (KEK-LS。[`pysad/tests/kek-ls_v03_12.lat`](../pysad/tests/kek-ls_v03_12.lat)、3 GeV、周長 603 m、3563 要素)、小さな FODO (3 GeV、24 セル)、MPW つき ASTRID2 (`RINGMPW`)。gs は 2026-11-19 の 10-19 ビルド (WSL)。elegant は 2026.4.0 (T1・T2・L2)。
- 新しい example: PF の履歴ラティス ([`example/pf/lattice/`](../example/pf/lattice/)・`PF_history.md`、1981・1986・1997・2005 の 4 つ)、KEK-LS の例集 ([`example/kek_ls/`](../example/kek_ls/)、v03_68、3 GeV、周長 570.72 m)、DDBA の例集 ([`example/ddba/`](../example/ddba/)、DDBA_01〜28 の版、5 GeV・2.5 GeV、周長 731〜757 m)、DDBA の最終版 v07_14 ([`example/ddba/v07_14/`](../example/ddba/v07_14/)、周長 749.48 m、30 セル)。

### (a) v2.3.0 で結果が変わるもの

| 項目 | 変わり方 | モード |
|---|---|---|
| LMA の `check_every` の既定 (L2) | 0 (SAD の tapert と同じ、15 要素ごと) から **1 (要素ごと)** へ。0 は警告つきで残る。0 は ASTRID2 の 25 位置で δ+ を最大 1.46 倍に見積もる。時間は ASTRID2 で +7%、大きなリング A で ×2.4 | 両 (新機能) |
| `FrequencyMap[]` の既定 (S2) | RFSW・RAD・FLUC OFF が既定。RF ON・RAD OFF の同期位相は 0 (gs と同じ)。RFSW ON では D が拡散の指標でなくなる | 両 (新機能) |
| `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` (T2) | **Piwinski 型の式 (`touschek.touschek_piwinski`) が v2 の既定**。K64 の式 (`touckf`) とは結合 1% で +0.1〜+6%、10% で +1〜+14%、丸いビームで −14〜+21% 違う (B162)。`Formula->"K64"`・k64 モードは K64 の表。`TouschekTable` (EMIT の報告・横の振幅の受容) は K64 のまま | v2 |
| CALC・`Emittance[]` の BEND (K1≠0) の行列の dk (B161、U2) | kick の導関数に直した (`compat.bend_dk_fix()`)。実際のリングでは ξ で 4e-6、∂²ν/∂δ² で 0.01、Emittance[] の EMITX・SIGE で 1e-9 (相対) 以内。K1 の強い扇形磁石では δ = 2% でチューン 1e-4 | v2 (k64 は K64 のまま) |
| `adts()` (R2〜R6) | 新機能の中での既定の選択: **既定は物理的な項だけ (`fourth='r3'`・`optics='calc'`)**。四極の積分器の台形則と偏向磁石の 4 次の項 (数値積分器の 0.1〜0.3% の端数) は精密版のオプション `fourth='r6'`、追跡から作る光学は `optics='tracked'` (V2 で R6 の既定から戻した)。δ≠0 の kick の強さの 1/(1+δ) 倍の誤り (R2 から) は修正して残した | 両 (新機能) |
| `resonance_terms()` (R4) | 既定で偏向磁石の 3 次の項を足す (`bend_cubic=True`。`.total_sext` が 6 極のみ) | 両 (新機能) |

K64 の動作 (k64 モード) は、新機能のほかは v2.2.0 と同じ。

### (b) 新しく使える機能

モードは**両** (k64・v2 で同じ)。SAD の台本の関数は pysad だけ (gs には無い)。使える範囲・条件・未確認は docs/46 の v2.3.0 の行。

| 機能 | 使える範囲と条件 | 使えない・未確認 | 例・記録 |
|---|---|---|---|
| **周波数マップ** (`freqmap.frequency_map`・`fmapcheck.py`、台本 `FrequencyMap[]`) | 1 つの位置から (x, y)・(x, δ)・(y, δ) の格子を追跡し、NAFF で前半・後半のチューンを出す。D = log10 √(Δνx² + Δνy²)。損失点は nan と損失周回数。表・図・共鳴線 (整数・半整数を必ず含む)。追跡は gs と座標 2.8e-7・D 4 桁。S2: 周回数の換算 `D_norm = D + log2(nturn/1024)`、誤差の統計 (種を複数)、警告 12 種、gs と誤差つき・RFSW ON が一致 | FLUC つきの D の意味、PF リング、SP8 (大きなリング B)。D の閾値は周回数で動く。結合が強いと νx・νy の対応が崩れる | `example/v2_3/Q1_*`・`S1_*`・`S2_*`、docs/48 §3・§10 |
| **周波数マップの共鳴線の予想** (`fmapres.py`・`fmappred.py`) | RDT (3 次・4 次、非共鳴の正規形の補正)・線形の ADTS・振り子の島の幅から共鳴の領域を予想し、D の大きい点と比べる (連続のスコア・信頼度)。FODO の 3νx の島では AUC 0.94 | **ASTRID2 (0.68)・KEK-LS (0.59) では振幅だけのスコア (0.78・0.78) に及ばない**。原因: 次数 5 以上の共鳴、島の中は規則的で D の大きい点は縁と重なり、大振幅の ν(J) の J² の項は v2.4.0 A1 の `ampnf` (docs/49 §13) で入れたが、予想は改善しない (追跡のチューンを与える診断でも AUC は上がらない。ν(J) ではなく 5 次以上の共鳴・島の重なりが原因、`docs/history/v24_A.md` §A6) | `example/v2_3/S3_*`、docs/49 §12、`FMAP_summary.md` §13 |
| **局所運動量アクセプタンス** (`lma.momentum_acceptance`・`lmabucket.py`、台本 `MomentumAcceptance[]`) | 位置ごとの δ±。粗い走査 + 二分法、放射・RF つき、RF バケットの外も損失、周回数 = シンクロトロン振動 2 周期 (既定)、並列。損失判定は gs と 30 / 30 一致、elegant と δ± が 0.5% 以内。L2: `check_every = 1`、RF バケットとアパーチャを区別 (`limit_types`)、`verify`、`turncheck`、周期性の確認、誤差つき・MPW は全周で走る、警告 9 種 | 判定は要素の境界 (BEND・QUAD の中は見ない)。k64 モード、PF リング、SP8 (大きなリング B)。FLUC は 2 周期では二分の幅以下だが 1 つのリングだけ | `example/v2_3/Q1_lma_*`・`L2_*`、docs/48 §4・§4.5 |
| **Touschek への受け渡しと Piwinski 型の式** (`lma.touschek_lifetime_local`、`touschek.touschek_piwinski`、台本 `TouschekLifetimeLMA[]`) | 要素ごとの損失率から 1/τ = Σ ½[R_l(δ+) + R_l(\|δ−\|)]。T2: 独立な式は 3 次元の一般の積分・実験室系の対のモンテカルロと 0.3% 以内、elegant と −1.6〜+0.3% | バンチ内の δ 分布が一様でない効果 (一様な受容だけ)、elegant の残り −0.3〜−1.3% の原因、IBS 込み (1e9 粒子/バンチ) だけ | `example/v2_3/T1_*`・`T2_*`、docs/48 §4.4・§4.6 |
| **追跡の設定の共通部分** (`trksetup`) | 閉軌道の Newton 法での求め直し、周期性の検出、RF バケット、優先度を下げた並列。RAD OFF の同期位相 0 (S2) | 空洞 (L = 0) とマーカーを除いて周期性を比べる。誤差のあるリングでは周期が無い | docs/48 §2 |
| **共鳴駆動項** (`resdrive.resonance_terms`・`rdt2nd.py`、台本 `ResonanceDrivingTerms[]`) | 6 極の 1 次の h (3 次、10 個の既定) + 偏向磁石の 3 次の項 + (R5) 4 次の h_jklm (6 極の 2 次・8 極・運動学・フリンジ・偏向磁石の 4 次。`order2=True` が既定)。追跡の 1 周の写像の 2 次・3 次のテイラー展開 (`rdttrack.py`) と、6 極の分は 6 リングで 1e-6 以内、偏向磁石つきの全体は 1e-8 (FODO・ASTRID2)〜4e-7 (大きなリング A)。R6: K0≠0 の偏向磁石 (MPW) は 5e-5 → 2.5e-9 | 結合・垂直分散・DX/DY のある素子 (警告)、K0≠0 の磁石の E・ROTATE、10 極以上。文献の δ の項の式との照合は記憶に基づく | `example/v2_3/Q2_*`・`R4_*`・`R5_*`・`R6_*`、docs/49 |
| **ADTS** (`resdrive.adts`、台本 `AmplitudeTuneShift[]`) | 6 極の 2 次 + 8 極の 1 次 + 運動学的な項 p⁴/8 + 四極・K1 つき偏向磁石の端部のフリンジ + 偏向磁石の 3 次の項の対 (`.total`。既定 `fourth='r3'`)。追跡 (pysad k64 = gs、70 粒子の NAFF + 多項式フィット) との差 (6 極 ×1): FODO 4e-4、ASTRID2 2.4e-3、KEK-LS yy 2.6e-3。**精密版 `fourth='r6'` は 1e-4 以下** (数値積分器の端数を含む)。差の主因は docs/49 §11.9 | 結合、垂直分散 (誤差で 10〜40% 動く)、端の 4 次、δ 微分 (既定は数十%、精密版も ASTRID2 の yy が 28% 外れる。数値の正規形を使う)、∂²ν/∂δ² の解析式 | `example/v2_3/Q2_adts_*`・`R2_*`〜`R6_*`、`ADTS_summary.md` |
| **数値の正規形** (`nform.py`) | 追跡の 1 周の写像の 3 次のテイラー展開から h3・h4 を取り出して ADTS を出す。δ 微分は `adts_delta_derivative(method='nform')` (追跡と 1e-3 以内) | J の 1 次まで。J² と大振幅 (KEK-LS) は未対応 | docs/49 §10 |
| **クロマティシティ** (`chromaticity_numeric`・`chromaticity_sextupole`・`chromaticity_natural`・`chromaticity_bend`、台本 `ChromaticityTerms[]`) | CALC の δ 依存の多項式フィットと 6 極・偏向磁石の 1 次の寄与。MPW: CALC − 解析 = (+6.2e-4, +3.5e-4) | ∂²ν/∂δ² の解析式なし | docs/49 §4 |
| **6 極の検討 (試作)** (`sxtopt.py`、例 [`example/ddba/v07_14/sextupole/`](../example/ddba/v07_14/sextupole/)、R7) | v07_14 (5 GeV、周長 749.48 m) の6 極 20 ファミリーの K2L の多項式モデル (1 次の RDT・δ の項は線形、4 次の RDT・ADTS は二次形式) + 解析的な勾配 + SLSQP。解析式の目的関数と追跡の DA 面積の順位相関 ρ = −0.89 (6 極の組 32 通り)。検証 (小〜中振幅の DA の健全性チェック) には条件つきで使える。**最適化の試作は最終設定より DA の広い解を出せなかった** (DA 面積は最終設定の 0.60〜0.71 倍。簡略化した追跡: 512 周、誤差なし)。上限 |K2L| ≤ 最終設定の最大 (仮定) | LMA・x-δ の DA の予測 (最終設定の近くで |ρ| ≲ 0.2)。J² の ADTS・δ 依存・5 次以上の駆動項は目的関数にない。完全な設定の DA・LMA での比較・軌道補正つきの誤差の比較は未実施 | [`example/ddba/v07_14/sextupole/README.md`](../example/ddba/v07_14/sextupole/README.md)、[`pysad/tests/test_ddba_sextupole.py`](../pysad/tests/test_ddba_sextupole.py) (8 件) |
| **B160 のオプション** (`compat.set_bend_k1_def('log'\|'lin')`、環境変数 `PYSAD_BEND_K1_DEF`) | 'lin' = B' 一様 (elegant CSBEND K2 = 0・MAD と同じ)。既定 'log' = K64 (K1 の暗黙の定義 (1 + hx)∂b/∂x 一様、docs/11 A59)。v2 だけ | k64 では無視。K1 の端のフリンジの比較なし | docs/history/v23_U2.md |

### (c) pysad 自身の不具合の修正

K64 のバグではなく、pysad のパーサの誤り (B 番号は付けていない)。

- 逆向きの文脈 (`-LINE` を逆向きに展開している中) の `-LINE` が二重に反転していた (`parser.py` `_expand_refs`)。PF の履歴ラティスの例で見つけ、gs と合うようにした。
- k64 モードで、数に隣接する単位 (`1.7MV`・`3GEV`・`10MM`) を構文エラーにしていた。gs は読む (D で始まる単位 `45DEG` だけが Fortran 流の D 指数 `45D` と読まれて構文エラー)。DDBA_06 の `CAVI VOLT=1.7MV` で見つけた。
- LMA・周波数マップ側: RF ON・RAD OFF の同期位相 (`trksetup.from_ffs`) を gs と同じ 0 にした (S2)。`TrackSetup` の K1 の設定 (`set_bend_k1_def`・`set_bend_k1_edge`) の引き継ぎ (S2)。

### (d) 要判断・未解決

各記録の「要判断」「未確認」を集めた。**解決したもの**: LMA の `check_every` の既定 (L2: 1 に変更)、RF バケットの判定の独立な確認 (L2)、二分法の見落とし (L2: `verify`。窓は 0 個)、周波数マップの規格化の既定 (S2: `norm='data'` のまま。`D_norm` の換算を足した)、Touschek 寿命の式のどちらが正確か (T2: Piwinski 型が正確、K64 の式は適用範囲に限界)、`adts()` の既定 (V2: 物理的な項だけ。精密版はオプション)、B160 の既定・オプション (U1・U2: 既定は K64 の定義のまま、'lin' はオプション)、ADTS の gs との不一致 (R1: 原因は gs の台本の RFSW)。

**判断が要るもの (既定は担当の選択のまま)**
- LMA の `check_every = 0` (SAD と同じ判定) を gs との比較用にだけ残すか (現状は警告つきで可)。`verify` の既定 0 (6 倍遅いので既定にしていない。窓は 0 個だった)。
- Touschek の Piwinski 型を v2 の既定にした (T2)。K64 の `TouschekTable` を使う EMIT の報告・横の振幅の受容は K64 のまま。
- 周波数マップの共鳴線の予想 (S3) は、判別力が振幅だけのスコアに及ばないリングがある。予想として使うか、診断に留めるか。

**未解決・未確認**
- 周波数マップ・LMA: k64 モードの LMA、PF リング、SP8 (大きなリング B)、FLUC つきの周波数マップの D の意味、誤差つき・結合のあるリングの LMA の統計 (全周で走る・警告が出る。統計は取っていない)。D の閾値は周回数で動く。周波数マップの全周の追跡は小さな格子 (336 点) と 1024 周が中心。LMA の判定は要素の境界で、BEND・QUAD の中は見ない。
- ADTS・RDT: 精密版の δ 微分の ASTRID2 の yy が 28% 外れる (原因は未特定)、K1≠0 の偏向磁石の本体の閉形式の 1e-4〜3e-4 (K64 の `tbendi` の積分器と推定、照合は未実施)、K0≠0 の偏向磁石に E1・E2≠0・ROTATE がある場合、共鳴線の上での非対角の補正、J² の正規形・大振幅 (KEK-LS)、∂²ν/∂δ² の解析式、結合・垂直分散・DX/DY、端の 4 次、10 極以上。MPW の ADTS の λ=2 の yy −0.5% は追跡側のフィットの偏り (Jmax・order で解析に収束)。
- Touschek: バンチ内の δ 分布が一様でない効果 (T2)、elegant の損失率が少し大きい残り (−0.3〜−1.3%) の原因、IBS 込み (電流 1 点) だけ。
- K1≠0 の BEND (B160・U1・U2): 「定義であって誤りではない」は中〜高の確度、K1 の端のフリンジ・E1/E2≠0・RAD の本体の比較なし、1 つのリング (n = 0.4 の 16 台)・δ = ±0.4% だけ。
- **ユーザー定義の関数名 `At` が pysad の `@` (内部の頭 `At`) と衝突して呼べない** (確認した: gs は `At[1]` = 101、pysad は警告を出して None。docs/46)。
- **`FitFunction` を `fit` と併用すると pysad が gs と別の解に落ちる** (KEK-LS の例。NaN の `FitFunction` では pysad は fit が動かない。原因は未確認)。
- DDBA_15 (等時セル + 逆曲げ BEND) の εx が gs と +0.33% (理由は調べていない)。DDBA v07_14 の 5 GeV の誤差計算の DA は gs が 0〜1σ、pysad が 27〜38σ (gs が補正後の状態を追跡していない疑い、原因は追っていない)。
- 互換モード: LMA・周波数マップは v2 で確かめた (k64 は未確認)。RDT・ADTS の追跡の比較は k64。

### (e) 試験と回帰の状態 (2026-11-21)

| 項目 | 結果 |
|---|---|
| 全試験 | `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本を `--ignore`、`PYSAD_SAD_PACKAGES` つき、互換モードは k64 が既定、`py -3.12 -m pysad.nbcache` の後): **1899 通過・8 skip・失敗 0 (1907 件収集。R7 の試験 [`test_ddba_sextupole.py`](../pysad/tests/test_ddba_sextupole.py) 8 件を含む)** (V3 の最終の通し。1759 通過・8 skip・失敗 0 は 2026-11-20 の V1 の時点。R6・S1〜S3・L2・T2・U2・X1〜X4 の試験が加わった) |
| 回帰 | `run_regression.py` (`QUICK=1`): **50 ベンチ FAIL 0 (`QUICK=1`、133 s、k64)** |
| skip | PYSAD_SLOW=1 が要る 4 件 ([`test_csr.py`](../pysad/tests/test_csr.py)・[`test_v2_C3_measure.py`](../pysad/tests/test_v2_C3_measure.py)・[`test_v2_D1_emit.py`](../pysad/tests/test_v2_D1_emit.py)・[`test_v2_E_pol_kick.py`](../pysad/tests/test_v2_E_pol_kick.py))、`PYSAD_ERR_LATTICE` が要る 2 件 ([`test_elkern.py`](../pysad/tests/test_elkern.py))、POSIX だけ・fork が要る 2 件 ([`test_v2_B_script.py`](../pysad/tests/test_v2_B_script.py)・[`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py)) |
| v2.3.0 の追加分 | [`test_v23_Q1_fmap_lma.py`](../pysad/tests/test_v23_Q1_fmap_lma.py)・`Q2_resdrive.py`・`R1`〜`R6`・`S1`〜`S3`・`L2`・`T1`・`T2`・`U1`・`U2`、PF 履歴・KEK-LS・DDBA の例の試験 ([`test_pf_history_lattices.py`](../pysad/tests/test_pf_history_lattices.py) ほか) |
| 窓の後始末 | [`pysad/tests/conftest.py`](../pysad/tests/conftest.py) の autouse の module fixture が、各試験ファイルのあと残った Tk のルート窓を閉じる (GUI の試験の内容・結果は変えない。全試験のあと窓が残らないことを確認) |
| 限界 | v2.3.0 の機能の gs との比較は、追跡と損失判定と ADTS の一部だけ。PF リング・SP8 (大きなリング B) の周波数マップ・LMA は使っていない |

## v2.2.0 (2026-11-20)

v2.1.0 の「要判断・未解決」を担当 P1〜P5 が追求した版 (記録は `history/v2next_P1.md`〜`P5.md`、例は `example/v2_next/P*_*`)。互換モード ([`pysad/compat.py`](../pysad/compat.py)、既定 v2、`PYSAD_COMPAT=k64` で K64 の動作) の使い方は下の v2.0.0 の (a) と同じ。**v2.1.0・v2.0.0 の節は変えずに残す**。この節の (a)〜(f) は v2.2.0 だけの一覧。

- 担当: P1 = 偏極の残り (脱偏極率、n が傾く・非平面のリング、反転した BEND)、P2 = 放射・エミッタンス (SynchroBeta、tradke、RADCOD の微分、K1≠0 の BEND の端)、P3 = TrackParticles の F1 と MULT の AE、P4 = SOL・空間電荷・DAPWIDTH・FREE、P5 = gs との検証の拡充 (MEASURE、DETR ≥ 1、複素の特殊関数、Micado、陽子の β、文法)。
- 使ったリング: 小さなリング (12 セル・24 セルの FODO、周長 18〜106 m、0.2〜4.2 GeV)、大きなリング A ([`pysad/tests/kek-ls_v03_12.lat`](../pysad/tests/kek-ls_v03_12.lat)、3 GeV、周長 603 m、3563 要素。KEK-LS は ESRF 型の HMBA (Hybrid Multi-Bend Achromat) で、ESRF の Pantaleo Raimondi 氏と Simone Liuzzo 氏の協力で設計した。)、PF リング ([`example/pf/lattice/v9_16_asym.sad`](../example/pf/lattice/v9_16_asym.sad)、2.5 GeV、187 m)、ASTRID2 ([`example/astrid/lattice/ASTRID2_v01_01.sad`](../example/astrid/lattice/ASTRID2_v01_01.sad)、MPW 入りは `ASTRID2_v01_02_id.sad`、0.58 GeV、45.7 m)。SP8 (大きなリング B、8 GeV、周長 1436 m、矩形の偏向磁石 88 台) は P2 が数値だけ確認に使った (SynchroBeta・EMITDIV・RADTAPER の閉じ方)。
- K64 の側のバグの番号は **B135〜B159** (25 件。[`21_sad_bug_report.md`](21_sad_bug_report.md))。B42・B108・B111・B112・B116・B128・B129 は更新。B128 (SynchroBeta の 40 倍) は P2 が原因を特定したので新しい番号は付けていない。
- gs は 2026-11-19 の 10-19 ビルド (WSL)。v2 の値の根拠は v2.1.0 と同じ (解析式、独立な数値計算、gs)。
- **k64 モードの結果を変えたのは、K64 の動作の追加 (これまで未対応・無視だったものを gs と同じに実装した) と pysad 自身の不具合の修正だけ**: MEASURE・TRCOD (先頭の位置)、PSPAC、TrackParticles の F1・MULT の AE、要素ごとの FREE の何もしない変数、DAPWIDTH の窓、MAIN の文法 (`SIN[0.1]`・`45DEG`・`1+`)。K64 の誤りを直した動作は v2 だけ。

### (a) v2.2.0 で結果が変わるもの (K64 = k64 モードの動作との違い)

「K64」は gs の動作で、k64 モードの動作でもある。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B128 | SynchroBetaEmittance (リングの始点の分散 η ≠ 0) | 源泉に閉軌道 c(δ) (1 次の分散を含む) をそのまま使い、x_b と δ の相関を入れない。始点に依って εx が数倍〜40 倍過大 (24 セル FODO の 3 通りの始点で 6.54・6.54・4.74 倍、大きなリング A 40 倍) | 源泉から η δ を引き、軌道の δ 依存を解析的な 2 項で入れる (`synchrobeta.tesumb_v2`)。Emittance[] との比: FODO 1.0003〜1.0005、大きなリング A 1.0084 (AzimuthalModes 3。32 で 1.0011)、SP8 (大きなリング B) 1.0045、PF 1.005、ASTRID2 1.006 | 合成のトイ (軌道 c = a δ の回転 + 減衰 + 拡散の解析式 1.02499 に v2 が 1e-3、MC 1.025〜1.029。K64 の式は μ に依り 1.027〜1.086) と、始点を変えた同じリングの値の比較 (gs でも再現) | `test_v2_P2_rad.py::test_synchrobeta_*` (11 件) |
| B111 | tradke の符号 (確定)・RADTAPER の追跡の BEND の重み | tradke の `ddpx` が追跡と逆符号。1 台の偏向磁石 (3 GeV・1.2 m・0.2618 rad) で既定の EMITDIV = 1 の Δpx が厳密値の 0.910 倍。RADTAPER は Emittance[] 3/4・追跡 1/2 | v2 の tradke は放射反作用つきのローレンツ力の数値積分と 6 桁一致。RADTAPER の BEND の係数は (δ_l + δ_l+1)/2 (調和平均) に統一 (v2.1.0 は 3/4) | 放射反作用つきのローレンツ力の積分 (`_p2_indep.bend_with_radiation`)。SP8 (大きなリング B) で係数の違いは軌道 2e-10 m・チューン・U0 は 1e-6 の桁 (Emittance[] と追跡の係数を揃えれば閉軌道は閉じる) | `test_tradke_sign_agrees_with_exact_radiation_reaction_in_v2_not_k64`、`test_radtaper_coefficients_v2_*` |
| B138 | tradke の微分 (RADCOD の解析 T + D) の 4 点の誤り | 4 点 (tr1 の始点・終点の取り違え、dddpx の 2 項目、c1 の微分、dxpa のべき)。BEND 1 台 5.1e-7、1 セル 1.9e-6、24 セル 4.3e-6 の数値微分との差。WSPAC の Emittance[] のチューンシフトで 2e-3 (相対) | 4 点を直す (`perfkern.tradke_db` の `exactd`)。BEND 1 台 1.7e-10、QF 1 台 6.7e-16、1 セル 6.4e-9、24 セル 1.9e-7 (残りは B139 = BEND の T[z, δ] の近似) | 写像の 4 次の中心差分 (刻み 5e-6〜1e-4 で同じ) | `test_radcod_analytic_matrix_matches_numeric_derivative_{v2,k64}` |
| B116 | Polarization・PolarizationTime の緩和率 (n が傾く・非平面のリング) | 逆向きの極は v2.1.0 (A2) の Σ r (b·n)² で数えるが、n ⊥ b (垂直の偏向・水平磁場の MPW) で 0。1 − (b·n)² の分だけ過小 | Sokolov-Ternov の厳密式 Σ r (1 − 2/9 (n·v)²) (`spin.ST_EXACT`)。分子に運動論的偏極 b·∂n/∂δ (`spin.KINETIC`)、脱偏極率を放射の刻みごとのスピン軸の場で積分 (`spin.STEP_RATE`)。v2.1.0 (A2) の式は水平磁場の MPW (ROTATE = 90°) で τ が 32% 短い | 独立な DK の数値積分 (刻み 0.05 m の追跡から N の周期解と刻みごとの ∂n/∂δ を作る: `_p1_dk.py`)。τ の比 (v2 / DK) 1.000001〜1.0016 (垂直の偏向、ソレノイド 1 T、縦ずれの四極、ASTRID2 + MPW の ROTATE = 0・0.3 rad・90°)、P の比 0.976〜1.0053 (0.5 GeV の FODO。ソレノイドと縦ずれのリングは脱偏極率の差が残って 2.2〜2.4%) | [`test_v2_P1_st.py`](../pysad/tests/test_v2_P1_st.py) (9 件) |
| B125 の更新 | Polarization の脱偏極率: ベータトロンのチューンの δ 依存 (クロマティシティ) | (K64 の `spdepol` は B125 のまま) | クロマティックな変調 (チューンの δ 依存 ξ = dν/dδ) を側帯の式に入れる (`spin.CHROMATIC`、`sideband_depolarization(..., C6=...)`)。4.2 GeV・12 MV の 12 セル FODO で式 / 追跡 (FLUC) = 0.43 → 1.05 ± 0.16。**原因の切り分け**: 振幅依存のスピンチューン (厳密に 0)・δ の非ガウス性 (無視できる) ではなく、ベータトロンの位相の変調 (変調指数 0.16) とスピンの位相の変調 (0.30) の差が、共鳴に近いリング (スピンチューンが νx + 2νs に離調 0.005) の側帯の重みを決める。軌道の非線形 (チューンの δ 依存) が応答を 2.7 倍に増やし、スピンの非線形は 0.75 倍で増やす向きではない | 決定論的な応答の測定 (NAFF)・FLUC の追跡・MC (`P1_pol_chromatic.py`)。表は `history/v2next_P1.md` §4.4: 2.903 GeV・1.5 MV 1.22、ASTRID2 0.99、3 GeV・1.5 MV 1.8 (安全側)、PF・大きなリング A は式の値だけ (クロマティックな変調は入れない) | [`test_v2_P1_pol_chrom.py`](../pysad/tests/test_v2_P1_pol_chrom.py) (6 件) |
| B135 | 追跡 (track6) のスピンが BEND の ROTATE で回らない | `tbshift`・`tent` が座標だけ回す。ROTATE = 90°・K0 = 0.064・L = 0.0875 m の BEND 1 台で追跡のスピン回転行列 R[0,2] = −0.1479、厳密 R[1,2] = +0.1479 | 座標と同じ回転をスピンに掛ける (`track6._spin_rot_cs`: K1=0・K1≠0 の BEND と ANGLE=0 の BEND の極)。厳密な Thomas-BMT の積分と 1e-9 | 厳密な Thomas-BMT の積分。Emittance[] の srot は一致 (gs の追跡のスピンで同じかは未確認) | [`test_v2_P1_st.py`](../pysad/tests/test_v2_P1_st.py) (2 件) |
| 不安定リングの P | Polarization (不安定なリング) | (K64 は P が NaN: B55) | **P = NaN** (`depol_error = 'unstable optics: …'`): 1 周の写像 T + D の固有値の大きさが 1 + 1e-6 を超える。v2.1.0 は不安定な軌道の上で線形の式を解いて 0.92 などを返していた (pysad 自身の不具合: (c)) | 大きなリング A の縦ずれ 0.2・0.4 mm (x の固有値 1.2)、垂直の偏向の組 (y が 1.8)。反ダンピングのシケインも NaN | [`test_v2_P1_rev.py`](../pysad/tests/test_v2_P1_rev.py) (6 件) |
| B140 | 追跡 `tsteer` の出口の AE2 | fb2 ≠ 0 のときだけ入れる (CALC・Emittance[] は常に) | 常に入れる。追跡の一周写像の固有値と CALC のチューンが大きなリング A (補正磁石 AE2 = 0.3) で k64 1.4e-4 → v2 6e-6 | gs で再現 (AE2 だけのステアリングの追跡のチューンが AE なしと同じ) | [`test_v2_P3_ae.py`](../pysad/tests/test_v2_P3_ae.py) |
| B141 | `tmulte` (CALC・Emittance[] の MULT) の出口の端の蹴りの位置 | 最後の半スライスの前 (追跡は後) | 後 (追跡と同じ)。MULT (K1 = 0.7・−0.6、K0 = ±0.01、AE2 = 0.3) の追跡と CALC のチューンの差 k64 3e-5 → v2 3e-8 | gs の CALC・追跡の差 3e-5〜1.2e-4 | [`test_v2_P3_ae.py`](../pysad/tests/test_v2_P3_ae.py) |
| B142 | K1≠0・K0≠0 の ANGLE=0 BEND の追跡の AE | `tmulti` に psi = 0 を渡し AE を無視 (ソースの読み) | AE を使う。同じ BEND の CALC のチューンと 2e-6 (k64 は 5e-4 以上) | gs はこの BEND で SIGSEGV (B02)。K64 の値の確認は無い | `test_bend_k1_k0_edge_angle_tracking_uses_ae_in_v2_only` |
| B143・B108 の拡張 | MULT・追跡 (`tsteer`・`tmulti`) の端の角の放射 | 行路を変えない (B75・B108 と同じ。端の項が約 2 倍: ΔJx の比 2.0〜2.1) | 行路を x tan(e) 短くする (追跡は `zr0`、Emittance[] は `tradke` の tin・tout)。厳密な幾何との ΔJx の比 1.03〜1.04 (MULT K1 = 0・0.06 の勾配つき、BEND の ANGLE=0、AE = (0.3, 0.3)・(−0.2, 0.4)) | 硬い端の極の厳密な幾何 + ロビンソンの定理 (`P3_b108_exact.py`)。~~残る 2.7% は原因未特定 ((d))~~ → 基準側の差と判明 (v2.4.0 D2、B222)。MULT の端の薄い蹴りは B223 で厳密にした | `test_edge_angle_term_v2_matches_exact_geometry_emit_and_tracking` ほか |
| B144 | 補正磁石の ROTATE と DX・DY を同時に付けた CALC | gs の CALC のチューンが追跡と違う (12 セル FODO で νy 0.1258 対 0.1337) | v2 の CALC は追跡と 5e-7 (B09 の修正が入る) | gs の CALC と追跡。k64 の CALC は gs と 1e-4〜7e-4 違い K64 を再現していない (未調査) | [`test_v2_P3_f1.py`](../pysad/tests/test_v2_P3_f1.py) |
| B146 | SOL 領域内の端の縦磁場の積分 ψ と QUAD の刻み | ψ が DISFRIN に関係なく入る。QUAD の刻み bz·al/ndiv ≈ 0.2 | ψ は DISFRIN=0 のときだけ。POL のとき QUAD を bz·al/ndiv ≤ 0.02 に刻む。領域内の QUAD・MULT のスピンが Thomas-BMT の積分と 9.3e-7 以下 (以前 QUAD 9e-5、MULT 1.6e-4) | Thomas-BMT の積分 (bz = 0.3〜9 1/m、0.4 m) と ASTRID2 の Emittance[] の SpinTune が追跡と 1.2e-7 | [`test_v2_P4_sol_spin.py`](../pysad/tests/test_v2_P4_sol_spin.py) (18 件) |
| B150 | SOL 領域内の DRIFT・QUAD・MULT の光子 (RAD + FLUC + PHOTONS) | 領域内の QUAD は記録しない (MULT・DRIFT は領域の始点または要素の範囲) | DRIFT・QUAD・MULT (ANGLE=0) で軌道上の位置に記録。エネルギーの合計が古典値と統計の範囲 (MULT で 309 eV/粒子、古典 285 eV。1500 粒子) | MULT K0 = 0.01 の放物線と 5e-7、RAD NOFLUC の損失が領域内・外 (bz = 1e-9) で 1e-17 | [`test_v2_P4_sol_photons.py`](../pysad/tests/test_v2_P4_sol_photons.py) (4 件) |
| B153 | 入れ替わる (DETR ≥ 1) INS の一周の位相の整数部 | 純粋な 90 度回転の位相の跳び ±π が丸め誤差の符号で決まる。pysad の k64 は ASTRID2・PF の INS・大きなリング A で gs と整数部が +1・+1・−1 ずれる | 入る要素の跳びを +π、出る要素を −π と決めて打ち消す (`lattice._exch_phase`、`perfkern._step` の `v2` 引数)。3 つのリングで gs と同じ整数部 | gs (小数部は 8 桁で一致) と、回転なしのリングの Emittance[] との 1e-6 の一致 | [`test_v2_P5_detr.py`](../pysad/tests/test_v2_P5_detr.py) (11 件) |
| 文法 | MAIN の `SIN[0.1]`・`45DEG`・`1+`・`E`・`Pi`・`Degree` | gs: `SIN[0.1]` 構文エラー、`45DEG` 構文エラー、`1+` は値 1 (B156)、`E=0.5;` の代入を要素定義の式が使う | `SIN[0.1]` は両モードで構文エラー (pysad は値 0 を入れていた誤りを直した)、`E`・`Pi`・`Degree` の代入を使う (両モードで直した)。`45DEG` は k64 だけ構文エラー (v2 は受ける)、`1+` は k64 だけ値 1 (v2 は構文エラー) | gs と 1 式ずつ比較 (`P5_README.md`) | [`test_v2_P5_syntax.py`](../pysad/tests/test_v2_P5_syntax.py) (13 件) |
| (P2 5) | K1≠0 の BEND の端の模型 | 四極の磁場を半径方向の面から始める (`radial`) | オプション (下の (b)) の `face` では四極も端面から始まる硬い端の磁石。**既定は `radial` のまま** (結果は変わらない) | ローレンツ力の厳密な写像 (`_p2_indep.exact_map`) と CALC 3e-10・追跡 1e-9 | `test_k1_bend_face_model_*` (5 件) |

- v2.1.0 の B111 の行 (RADTAPER の追跡の重み (3δ_l + δ_l+1)/4) は、v2.2.0 で (δ_l + δ_l+1)/2 に変わった (上の B111 の行)。v2.1.0 の節の記述は履歴として残してある。
- B108 の ANGLE=0 の BEND (補正磁石) の「傾き 30 mrad で 0.02〜0.03」は、30 mrad に特有ではなく K0 に依らない 2.7% の系統差だった (AE に比例、極の長さ・間隔に依る)。
- 数値の効果が小さいもの (放射の式の 6 桁の一致、チューンの変化 1e-5 の桁) は記録の表にある。

#### 互換モードごとの違い (v2.2.0 で追加・変更したオプションと既定)

| 項目 | k64 | v2 (既定) | 切り替え |
|---|---|---|---|
| Polarization の緩和率・運動論的偏極・脱偏極率の刻み・クロマティックな変調 | K64 のまま (Σ r (b·n) など) | `spin.ST_EXACT`・`KINETIC`・`STEP_RATE`・`CHROMATIC` = True | モジュール変数を False にすると v2.1.0 の式 (`CHROMATIC = False` で B3 の式) |
| K1≠0 の BEND の端の模型 | 常に `radial` | 既定 `radial`。`face` を選べる | `compat.set_bend_k1_edge('face')` または `PYSAD_BEND_K1_EDGE=face` (k64 は設定しても `radial`) |
| 追跡の WSPAC のビーム行列 | 追跡に WSPAC を入れない | 既定はゼロ電流の平衡行列。`WSPACSC = 1;` で空間電荷つき (MINCOUP の下限つき) | FFS の大域変数 `WSPACSC` |
| DAPWIDTH の既定 | 7 (gs と同じ) | 7 (v2.1.0 は 51) | `DAPWIDTH = 51;` で全点 |
| PSPAC | K64 の形状関数 (v2.2.0 から。`spch_k64.py`) | 2 次の B スプライン | モード |
| MEASURE・TRCOD | 先頭の位置だけ (gs と同じ。種を戻さない) | 追跡を通す (種を戻す) | モード |
| RADTAPER の BEND の係数 | Emittance[] 3/4・追跡 1/2 (K64) | 両方 1/2 | モード |
| 複素の特殊関数 (追加した 15 関数) | 真値 (gs が誤る点も) | 同じ | — |

### (b) モードによらない追加・修正と、新しく使える機能

「モード」の列: **v2** = v2 モードだけ、**両** = 両モードで同じ、**k64** = k64 モードだけ。例は [`example/v2_next/`](../example/v2_next/) (`.gs.txt` = gs、`.pysad.txt` = v2、`.pysad_k64.txt` = k64)。

| 機能 | モード | 使える範囲と条件 | 使えない・未確認 | 例・記録 |
|---|---|---|---|---|
| **n が傾く・非平面のリングの Polarization** (`spin.ST_EXACT`・`KINETIC`・`STEP_RATE`) | v2 | 垂直の偏向、ソレノイド (1 T・0.8 m、n が傾く)、水平磁場の挿入光源 (ASTRID2 + MPW の ROTATE = 90°)、垂直の分散のあるリング。0.5〜3 GeV の FODO、ASTRID2。τ は DK の積分と 1e-5〜2e-3、P は ±1% (脱偏極率が分母の 10% 以下) | P の脱偏極率の部分は Emittance[] の D (追跡の D と違う: B136・B137)。分母の 15% を超える (傾いた MPW・強いソレノイド) と 2〜3% 違う。傾いた n の大きなリングの独立な計算は未実施。運動論的偏極は ∂n/∂δ (位置・角度を固定) の δ 列 | `P1_pol_vertical_bend.sad`・`P1_pol_mpw_rotated.sad`・`P1_pol_reversed_bend.sad`・`P1_st_dk_check.py` |
| **クロマティックな側帯の式** (`spin.CHROMATIC`、`CHROM_KMAX = 6e4`・`CHROM_NEL = 1200`) | v2 | 変調指数 (ベータトロン) ≥ 0.03 かつスピンチューンが νx + kνs の近く (離調 ≲ 2νs) のリング。4.2 GeV・12 MV の 12 セル FODO で式 / 追跡 1.05 ± 0.16。そうでないリングは B3 の式と ±5% | 追跡の測定は 4.2 GeV・12 MV が主 (6 MV・24 MV・PF・大きなリング A は式の値だけ)。K > 6e4・要素数 > 1200 (大きなリング A) では入れない。6 極のあるリングの振幅依存のチューンシフトと小さな ξ の模型は未検証 | `P1_pol_chromatic.py` |
| **反転した BEND のシケインのリング** (A2 の確認の拡張) | v2 | 曲げ角の和が 0 のシケイン (+0.05, −0.10, +0.05 rad) を 12 セル FODO (3 GeV)・大きなリング A・PF・ASTRID2 に入れて P・τ が独立な式と 1e-4 以内 (大きなリング A の P は −0.105699 で負) | 反ダンピングの大きなシケイン (固有値 > 1) は P = NaN | [`test_v2_P1_rev.py`](../pysad/tests/test_v2_P1_rev.py) |
| **`compat.set_bend_k1_edge('face')`** (K1≠0 の BEND の端の模型 `face`) | v2 | K1≠0 かつ ANGLE≠0 の BEND、E1・E2 ≠ 0、DISFRIN=0。ローレンツ力の厳密な写像と CALC 3e-10・追跡 1e-9。軌道が中心のリングでは `radial` と完全に同じ。違いは軌道がずれた BEND の x0 tan e の項 (x0 = 1 mm で 1.1e-8、K1 = 0.3 で 4.5e-7) | DISFRIN=1、ANGLE=0 の K1≠0 の BEND・MULT の端、SOL 領域内、RADINT、F1・FB の端。放射 (tradke) は楔形の蹴りを偏向角に入れない。`face` では numba の高速な経路 (`_bend_step_fast`・`_bende_fast`) を使わない | `P2_bend_k1_face.py`、`P2_README.md` |
| **TrackParticles の F1・FB1・FB2・FRINGE** (`tsteer` (ANGLE=0 の BEND))、**MULT の F1・F2・F1K1F/B・FB1・FB2・FRINGE = 1・2・3** (`tmulti`) | 両 (K64 の動作の追加) | k64 は gs と 10 桁表示の全桁一致 (12 セル FODO: ステアリング 6 通り、MULT 9 通り、大きなリング A の補正磁石 2 台)。numba と Python は 1e-13。独立な基準: 軸上の偏向 1e-14、放射の積分 1e-5、Maxwell 場の端の運動 (fb1 = fb2 で x が 1.3e-11)。F1 つきの MULT を含むリングは Python の経路に落ちる (1 周 1 粒子 0.4 s 程度) | 逆向きの要素 (`direlc < 0`) の入れ替え、RF つき MULT の F1、POL の MULT の F1 の端の項。fb1 ≠ fb2 で K64 の模型は厳密な Maxwell 場と (K0/al)(fb1² − fb2²)/24 違う (B145。直さない) | `P3_steer_f1_track.sad`・`P3_mult_f1_ae_track.sad`・`P3_maxwell_ramp.py` |
| **MULT の AE1・AE2** (CALC・Emittance[]・追跡) と K1≠0・K0≠0 の ANGLE=0 BEND の AE (CALC・Emittance[]) | 両 (これまで黙って無視していた) | k64 は gs と全桁一致 (5 通り、大きなリング A)。大きなリング A (補正磁石 2 台、AE ≤ 0.5 rad、K0 = 5 mrad〜30 mrad の極) で追跡と CALC のチューンが 2e-7 以内 (v2) | K0 の 2 次の追跡 − CALC の差 (K0 = 0.003 で 1.6e-5。K64 の端の蹴りの 2 次の違い、直さない)。放射の残る 2.7%。実験室系で y に傾く端の面 (K1 < 0 の MULT・ROTATE) は厳密な幾何との比較なし | `P3_mult_ae_emit.sad`・`P3_steer_ae2_tune.sad`・`P3_ring_a.py`・`P3_b108_exact.py` |
| **PSPAC の K64 の形状関数** ([`pysad/spch_k64.py`](../pysad/spch_k64.py)) | k64 (新しい機能) | 横 CIC・縦「ほぼ最近接」・5 点差分 + 双 3 次補間・3 次スプライン。gs と rms が 12 桁以上一致 (5 MeV・1e11 個・5 条件、多周回・凍結 4 条件、ASTRID2 に SPCH 3 か所・PBUNCH = 1e14・10 周で rms x を 14 桁)。**NPARA = 1 で比べた** | gs の NPARA > 1 (子プロセスごとの電位) の再現、TrackParticles ごとに新しく始める (gs は 2 回目に PSPAC が働かない)、PSPACNZ < 3 はエラー (gs は SIGSEGV)、SPAC・WSPAC・ウェーク・放射との併用、PF・大きなリング A | `P4_pspac_k64_shape.sad`・`P4_pspac_astrid2.sad` |
| **`WSPACSC = 1;`** (追跡の WSPAC に空間電荷つきのビーム行列) | v2 (K64 にない) | 0.2 GeV の 12 セル FODO、PBUNCH = 1e6〜1e8、MINCOUP = 0.5、RF あり。チューンシフトが同じ行列の一次の式と 2e-4 | gs の追跡と合わない (MINCOUP の下限が入るため)。強い空間電荷 (PBUNCH ≳ 1e9)・RF なし・大きなリング | `P4_wspac_tracking_beam.sad` |
| **DAPWIDTH の窓** (`dapert.dynamic_aperture` が `trackd.f` のスロット管理の直訳) | 両 (K64 の動作の再現。DAPWIDTH の既定を 51 → 7) | NPARA = 1、`Axes->"ZX"`、`Output->6` か 0。NORAD の地図・Score が gs と一致 (24 セル FODO で DAPWIDTH 既定・3・12・51、大きなリング A の 60 周)。k64 は RAD + FLUC の最終 seed まで一致 (DAPWIDTH 1〜6) | NPARA > 1、`Damping->True`、v2 は放射の式が違うので RAD + FLUC の seed は合わない。損失の境界の粒子は 1〜数周違う | `P4_dapwidth_window.sad`・`P4_dapwidth_ringA.out.txt` |
| **MEASURE・TRCOD・LOGMEASURE** (k64) | k64 (K64 の動作の追加) | 先頭の MARK の位置だけ gs と同じ統計 (重心・分散行列・エミッタンス・β・R1..R4・σ が 4 桁一致。FODO + skew 四極、大きなリング A・PF・ASTRID2)。それ以外の位置は無出力 (`FFS["MEASURE"]` は np = 1)、TRCOD は警告 (`SIGSEGV (B42)`)。種を戻さない (B152)、TRCOD の粒子数は最後の `tfsetparam` の NP (B154)、`NP` の既定 1000 | v2 の追跡を通した値は gs と比べられない (gs が追跡で落ちる: B151)。TRCOD FIXCOD 等・GETMEASURE・FIFO は未確認 | `P5_measure_*.sad`、`P5_README.md` |
| **要素ごとの FREE の何もしない変数** | 両 (gs と全桁一致) | `FREE QF.3 K1` は結合先 QF.1 を持つ何もしない変数として VAR に出る。変数の数・VAR の行・競合のエラー (`Element already used as variable`)・FIX・Saved の欄・RESET が gs と一致 | `COUPLE` で結合した要素の成分、`ElementValues` との連動、`BYE ALL`、反復の途中の overflow 行の数 | `P4_free_noop_variable.sad` |
| **SOL 領域内の DRIFT・QUAD・MULT の光子** | v2 | 3 GeV、bz = 0.1 1/m、傾き・ずれのない領域 | 領域内の BEND・SEXT・OCT、傾いた・ずれた領域の位置 (`geo_frac`)、偏光 ξ の値、gs との位置・数の比較 (gs は QUAD を記録せず MULT の偏向を無視する) | `P4_sol_photons.sad` |
| **複素の特殊関数 15 関数** ([`pysad/specfn_cx2.py`](../pysad/specfn_cx2.py): ExpIntegralE・ExpIntegralEi・PolyGamma・WhittakerM/W・ProductLog・CosIntegral・SinIntegral・FresnelS/C・HarmonicNumber・GammaRegularized・Beta・Factorial・`Gamma[a,z]`・LegendreP・複素の次数の BesselJ) | 両 | 複素の引数が 1 つでもあれば複素の実装 (実数は従来のまま)。732 式の走査で gs と 1e-9 が 436 → 641 式。残り 91 式は gs が誤る・不正確 (B158) で、pysad が独立な基準 (scipy・級数・恒等式・微分方程式の残差) と合う。\|z\| ≲ 30 | 複素の次数の ExpIntegralE・`HarmonicNumber[z,r]`・`LegendreP[ν,m,z]` (m ≠ 0)・複素の a の GammaRegularized・BesselY/I/K の複素の次数 (SADError か未確認)。k64 が gs の誤りを再現しない (要判断) | `P5_specfn_cx_scan.sad`、[`test_v2_P5_specfn_cx.py`](../pysad/tests/test_v2_P5_specfn_cx.py) (23 件) |
| **gs との照合 (コードの変更なし)** | — | `OYStability` の後半 (CSROYMatrix → Eigensystem → 成長率) が gs と 1e-12。Micado の水平・2 平面同時・`SVDmode->1`・`MakeBump` が gs の k64 と 5e-6 (`SVDmode->1` は n を無視: B157)。FITP>1 (3 例) の最終の残差が gs と 2e-4 以内。陽子 (1.5・3・30 GeV/c) の U0 の k64 が gs と全桁。Python 3.14.9 + Tk 9.0 の GUI が動き窓の大きさが Tk 8.6 と同じ (PNG 4 枚) | `OYStability`・`OYScan` の通し (gs が動かない)、`CorrectOrbit` の `Bump->` (B155)、FITP>1 の反復の経路、大きなリングの FITP>1、gs の窓の撮影 | `P5_oy_pipeline*.sad`・`P5_micado_*.sad`・`P5_fitp_*.sad`・`P5_proton_beta_*.sad`・`P5_tk90_*.png` |

### (c) pysad 自身の不具合の修正 (K64 のバグではない。B 番号を付けていない)

| 項目 | 不具合 | 修正 | モード・根拠 |
|---|---|---|---|
| `DynamicApertureSurvey` が常に全 51 点を追う | DAPWIDTH の既定 51、`Maximum number of particles` も違う (gs 42、pysad 306)。窓の動きで追跡する点が変わり RAD + FLUC の乱数の列が変わる | trackd.f の直訳 (上の (b)) | 両。gs と地図・Score・seed が一致 |
| 不安定なリングの P が 0.92 などを返す | 閉軌道は見つかるが光学が不安定なとき、不安定な軌道の上で線形の式を解いて値を返していた (v2.1.0 A2 の表の 2 行は意味がなかった) | P = NaN と `depol_error = 'unstable optics: …'` | v2。垂直の偏向・ソレノイドのリングでも以前は黙って返していた |
| MULT の AE1・AE2 を読まない | `_mult_params` が読まず、CALC・Emittance[]・追跡のどれでも黙って無視 (`tbende` から tmulte に回る ANGLE=0 の K1≠0・K0≠0 BEND の AE も) | 読む (上の (b)) | 両。k64 は gs と全桁 |
| 複素引数で虚部を黙って捨てる | `ExpIntegralEi`・`SinIntegral`・`FresnelS/C`・`HarmonicNumber`・`LegendreP` が実数を返した。`ExpIntegralE`・`PolyGamma`・`WhittakerM/W`・`ProductLog`・`CosIntegral`・`GammaRegularized`・`Beta`・`Factorial`・`Gamma[a,z]`・複素の次数の `BesselJ` は何も返さなかった | 複素の実装 (上の (b)) | 両 |
| `SIN[0.1]` が値 0 | 角括弧の関数呼び出しを SAD の文法で構文エラーにすべきなのに値 0 を入れていた | 構文エラー | 両。gs と同じ |
| MAIN の `E`・`Pi`・`Degree` の代入を要素定義の式が使わない | SADScript の定数 e・π を使っていた | 代入を `env._main_user_vars` に入れ、要素定義の式だけが見る (FFS の SADScript の `E` は変わらない) | 両。gs と同じ |
| 要素ごとの FREE の何もしない変数 | 先頭でない成分の変数を作らず、gs と変数の数・VAR の行が 1 つ違った | 作る (上の (b)) | 両。gs と全桁 |
| `track6` の領域内の QUAD の刻み | 弦の誤差 (bz·al/ndiv ≈ 0.2) | POL のとき bz·al/ndiv ≤ 0.02 | v2 (k64 は K64 のまま) |

### (d) 要判断・未解決

各記録の「要判断」「残した項目」を、v2.2.0 で解決したものを除いて集めた。**要判断の項目は、各担当の選んだ既定のままにしてある。判断が要ることと選択肢を書く (既定は変えていない)。**

**判断が要るもの (既定は担当の選択のまま)**

| 項目 | 現在の既定 (担当の選択) | 判断が要ること・選択肢 |
|---|---|---|
| EMITDIV の既定 (P2) | 1 (K64 と同じ)。v2 の Emittance[] は SP8 (大きなリング B、矩形磁石) で +1.4e-3 の離散化の誤差を含む | v2 の既定を 2〜4 にするか。変えると v2 の多くの試験の基準値が動く。当面は RADINT との比較に EMITDIV ≥ 2 (できれば 4) を使う |
| DAPWIDTH の既定 (P4) | 7 (gs と同じ。v2.1.0 は 51) | 51 を前提にした台本 (DAPWIDTH を書かないもの) は窓の外が `*` になり Score が変わる。全点を追うには `DAPWIDTH = 51;`。v2 の既定を 51 に戻すか |
| `WSPACSC` (P4) | 0 (追跡のビーム行列はゼロ電流の平衡行列。gs の追跡と 1〜2%) | 既定を空間電荷つきの行列に変えるか (K64 が `EMIT` の beamsize を使う意図)。`WSPACSC = 1` は K64 にないオプションで gs の追跡とは合わない (MINCOUP の下限が入る) |
| K1≠0 の BEND の端の模型 (P2) | `radial` (K64 の仕様。k64 は常に `radial`) | v2 の既定を `face` (厳密な硬い端の磁石) にするか。軌道が中心のリングでは結果は同じで、違いは軌道がずれた BEND の x0 tan e の項 (1e-8〜1e-7) |
| 放射の β の因子 (B112、P5) | 直していない (v2 の U0 は古典式の 1/β³ 倍: 陽子 1.5・3・30 GeV/c で 64%・15%・0.15%) | 光子数と光子のエネルギー分布の β 依存を厳密な放射のスペクトル (β < 1) で入れ、独立な基準 (厳密な放射の式の積分) と比べて直すか。`crad` に β³ を掛ける近道は量子励起を無視するので入れなかった。当面は低い γ の陽子で放射を使わない |
| SOL 領域内の端の縦磁場 ψ (B146、P4) | DISFRIN=1 で入れない (v2)。DISFRIN=0 では入れる | DISFRIN=0 の ψ が物理的に正しいかは、硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準が要る (未着手) |
| 偏向磁石の本体の B∥ (B136・B137、P1) | 直していない (k64・v2 とも K64 のまま) | tradke・tradk・tbend の全ての要素の本体に (1+a) B∥ ds を足しフレームの扱いを揃えるか。B58・B64・B130 のスピンの結果が全部動く。直すと脱偏極率の D の追跡との食い違い (鉛直の軌道・分散の大きいリング) も減る見込み |
| 複素の特殊関数の k64 の値 (B158・B118〜B120、P5) | k64 も真値 (gs の誤りを再現しない) | K64 の `gamma.f`・D.n の直訳 (規則が単純な LogGamma の分枝・ArcCosh の符号・Erf/Erfc の Re z < 0・ExpIntegralEi から) で k64 に gs の誤りを再現するか |
| `CSROYMatrix` の gs の wrongtype (B60、P5) | 再現しない (pysad は `CSRSetupOY` の出力を受ける) | 再現すると k64 で `OYStability`・`OYScan` が使えなくなる |
| 文法 (P5) | `45DEG` は v2 が受ける・k64 は構文エラー。`1+` は v2 が構文エラー・k64 は値 1 | v2 の寛容さ (前者) と入力ミスの検出 (後者) を gs に合わせるか |
| BEND の T[z, δ] の近似 (B139、P2) | 直していない (K64 の仕様。1 周で 1e-7) | v2 で `tbendebody0` の (5,6) の px² の項を厳密にするか |
| `tsteer` の F1 の fb1 ≠ fb2 (B145、P3) | 直していない (K64 の模型の近似) | 磁場の立ち上がりの中心を要素の端に置く取り方の判断が要る |
| K1≠0・K0≠0 の ANGLE=0 BEND の追跡・CALC の AE (B142、P3) | v2 は AE を使う (ソースの読み) | gs が SIGSEGV で K64 の値を確認できないので、読みが正しいかは未確認 |
| SOL 領域内の加速 (POL + TRPT、P4) | 実装せず明示のエラー | RF の電場を含む運動方程式の積分の独立な基準を作って、`_RadHook` の基準運動量の更新と `_tradk_sol` の p(s) を実装するか |

**未解決・未確認 (記録ごとの限界)**
- P1: 脱偏極率の追跡との式 / 追跡の比は **4.2 GeV・12 MV の 1 点** (6 MV・24 MV は FLUC のスピン追跡が 1 条件数時間で未実施、PF・大きなリング A は 1 粒子 1 周 0.1〜0.4 s)。式のモード分解は応答を追跡の 1.5 倍に見積もり、結合項 b(u) の 2 次の欠落 (−10%) と打ち消している可能性。大きなリング A・PF はクロマティックな変調を入れない (要素数・K の上限)。6 極の振幅依存のチューンシフトと、6 極で補正した小さな ξ の模型は未検証。n が傾く大きなリングの独立な計算は未実施。P の脱偏極率の部分は Emittance[] の D (追跡の D と違う) で、分母の 15% 以上では 2〜3%。**v2.2.0 の spin 系の既存の試験の許容を緩めた箇所**は (e)。
- P2: SB のシンクロ・ベータ共鳴の効果の長時間の追跡との比較は未取得 (CPU が必要)。SB と Emittance[] の差 0.3〜0.8% の原因 (物理か離散化か) は区別していない。RF が分散のある場所にある矩形リングの 6D の εx が独立な 6x6 の平衡と +1.7e-3 違う件は**原因を特定できていない** (どちらが誤りか決められない)。`face` の RADINT・DISFRIN=1・ANGLE=0。(4) の微分は大きなリング A・B・PF・ASTRID2 では未実施 (要素単体で確認)。
- P3: ~~B108 の残る 2.7% の系統差 (AE に比例、極の長さ・間隔に依る)~~ (v2.4.0 D2 で基準側の差と判明、B222) と AE のある面の和則の破れ 2.7e-4 は原因未特定。逆向きの要素の F1・AE の入れ替え。F1 つきの MULT の追跡の高速化。RF つき MULT の F1。実験室系で y に傾く端の面 (K1 < 0・ROTATE の MULT) は K64 の約束 (面が要素と一緒に回る) が意図かの判断。POL の MULT の AE の bsi の端の項。
- P4: 大きなリング A での確認は DAPWIDTH だけ。PSPAC (3 SPCH・10 周)・SOL 領域内のスピンは ASTRID2 (0.58 GeV、45.7 m) で代えた。603 m のリングでの PSPAC・領域内の光子・WSPAC の追跡は未確認。PF リングは未使用。DAPWIDTH のスロット管理は Python のループ (100 行 × 51 点 × 100 周で数秒)。領域内の CAVI のスピンの確認は弱い (位相 0)。
- P5: v2 の追跡経由の MEASURE・TRCOD の値は gs と比べられない。MARK の生の DETR ≥ 1 の INS (B159)・MEASURE の入口が DETR ≥ 1。`OYStability`・`OYScan` の通しと実用の物理の独立な基準 (Vlasov 方程式の既知の解など)。`CorrectOrbitBump` の直接呼び (gs が 300 s で終わらず)。FITP>1 の反復の経路 (gs は 1〜2 反復早く止まる) と大きなリング・PF・ASTRID2 の FITP>1。Tk 9.0 の `cget` の単位つきの値 (`4p`) を数として使う SAD の `.n`、gs の窓の撮影、高 DPI。陽子の Touschek・Intrabeam・RADINT の古典半径。
- 従来から: RADINT の K0 の極の εx (3.5%)・Jx (0.7%) の差、SEXT・BEND の軌道のずれ (B126)、空洞の位相。WSPAC の CALC (4D) の直訳の検証。`TrackParticles` の POL の領域内の四極のずれ。MAIN の FOURIER。`Parallelize` の Windows の並列。numba のキャッシュの自動無効化 (`pysad.nbcache` を手で走らせる)。SP8 (大きなリング B) の POL (四極 60 台に 0.1 mm で不安定)。

### (e) 試験と回帰の状態 (2026-11-20)

| 項目 | 結果 |
|---|---|
| 全試験 | `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本 [`test_bend_feeddown.py`](../pysad/tests/test_bend_feeddown.py)・[`test_coord.py`](../pysad/tests/test_coord.py)・[`test_sext_radiation.py`](../pysad/tests/test_sext_radiation.py) を `--ignore`、`PYSAD_SAD_PACKAGES` つき、互換モードは k64 が既定、`py -3.12 -m pysad.nbcache` の後): **1646 通過・9 skip・失敗 0** (1655 件収集、26 分 51 秒、他の重い計算なし)。このあと [`test_v2_A_spin.py`](../pysad/tests/test_v2_A_spin.py) に 1 件足して (下)、そのファイルを再実行して 22 件通過 → **1647 通過・9 skip・失敗 0** (1656 件収集)。v2.1.0 は 1345 通過・8 skip |
| v2.2.0 の追加分 | `test_v2_P1_{st,pol_chrom,rev}.py` 21 件、[`test_v2_P2_rad.py`](../pysad/tests/test_v2_P2_rad.py) 22 件、`test_v2_P3_{f1,ae}.py` 124 件、`test_v2_P4_{free,dapwidth,pspac_k64,wspac_beam,sol_spin,sol_photons}.py` 57 件、`test_v2_P5_{measure,detr,micado_fitp_oy,specfn_cx,proton_beta,syntax,tk90}.py` 78 件 (計 302 件。単独で走らせて 302 通過・skip なし、11 分 43 秒) と、[`test_v2_A_spin.py`](../pysad/tests/test_v2_A_spin.py) の 1 件。ほかに独立な計算の補助 `_p1_dk.py`・`_p1_models.py`・`_p2_indep.py`・`_p3_indep.py`・`_p3_gsdata.py`・`_p5_util.py`・`_p5_gscx.py` |
| 回帰 | `run_regression.py` (`QUICK=1`): **50 ベンチ FAIL 0** (k64。138 s) |
| 既存の試験の変更 (動作が変わったための追従) | [`test_v2_B3_pol.py`](../pysad/tests/test_v2_B3_pol.py) 3 か所 (変調なしの側帯の式との比較を要素の出口の N の率 `rate_el` に、4.2 GeV・12 MV の脱偏極率の期待値 4.46e-8 → 1.09e-7 (クロマティック込み。`spin.CHROMATIC = False` で 4.46e-8)、P の恒等式に運動論的偏極の項)、`test_v2_A2_pol.py::test_large_ring_a_orbit_error_depolarization_scales_with_square_v2` (ずれの rms を 0.1・0.2 mm → 0.05・0.1 mm。0.2 mm は光学が不安定で NaN)、`test_v2_A_spin.py::test_polarization_with_orbit_errors_v2` (下)、`test_v2_C_track.py::test_radtaper_matches_gs[v2]` (v2 のチューンの許容 2e-5 → 5e-5)、`test_v2_C2_optics.py::test_radtaper_in_calc6d` (v2 の軌道の許容 5e-3 → 1.5e-2)、[`test_v2_D1_emit.py`](../pysad/tests/test_v2_D1_emit.py) の WSPAC のチューンシフトの期待値 (−4.1386e-3・−8.4669e-3 → −4.1295e-3・−8.4361e-3) と SB の警告の試験 (大きなリング A で収束して警告なしに)、`test_v2_B2_free.py::test_mixed_variables_fit_and_table_order_like_gs` (VAR の並び)、`test_v2_C2_spin.py::test_spin_path_follows_the_same_orbit_as_fast_tracking` (許容 1e-10 → 5e-9)、[`test_v2_C_measure.py`](../pysad/tests/test_v2_C_measure.py) の k64 の MEASURE・TRCOD の試験 (gs の実動作に) |
| **spin 系の緩めた許容の確認** (`test_v2_A_spin.py::test_polarization_with_orbit_errors_v2`) | P1 と P2 の両方の変更の後で**通る**。P1 は許容を 3e-3 → 1.5e-2 に緩めていた (3 GeV・12 セル FODO、四極の縦ずれ 3 mm、期待値 0.8798)。実際の P は **0.8695** (−1.0e-2) で、内訳を `spin.KINETIC`・`STEP_RATE`・`ST_EXACT` を 1 つずつ切って確かめた: 運動論的偏極 −7.1e-3、厳密な緩和率 −3.0e-3 (τ 124.18 → 123.73 s。平面リングの 123.8 s に近づく)、刻みごとの脱偏極率 −2.5e-4、`CHROMATIC` は P に効かない。3 つを切ると 0.87983 で v2.1.0 の 0.8798 に戻り、P2 の変更 (tradke の微分など) はこの値を 4 桁で動かさない。したがって緩めた許容は妥当 (P1 の変更だけで説明できる) だが、5 倍緩いのは緩すぎるので、許容を **0.8695 ± 3e-3** に戻し、3 つの項を切ると 0.8798 ± 3e-3 に戻ることを確かめる試験 `test_p1_terms_account_for_the_change_of_polarization_v2` を足した。独立な基準: P1 の記録の 0.5 GeV の縦ずれのリングで、運動論的偏極を入れた P は DK の積分と 2.2% で、入れないと 1.7% 悪い (大きな 3 GeV のこのリングの DK の積分は未実施) |
| skip | GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip。`PYSAD_SLOW=1` で遅い試験が走る (SB の収束の警告、POL の kick の応答、大きなリング A の MEASURE)。9 skip の内訳は個別に調べていない (v2.1.0 の 8 + 1。P の試験は skip なし) |
| Python 3.14 | v2.1.0 の確認 (GUI 系 25 件と D2 の 90 件) のまま。P5 が Python 3.14.9 + Tk 9.0 の GUI の起動と窓の大きさ (Tk 8.6 と同じ) を確認した ([`test_v2_P5_tk90.py`](../pysad/tests/test_v2_P5_tk90.py) は専用の venv があれば走り、無ければ skip)。全試験は 3.14 では走らせていない |
| 限界 | v2 の値は gs では確かめていない (gs が K64 のバグを持つため)。SP8 (大きなリング B、8 GeV) は数値だけ (試験・例に使わない)。実機の強い結合・挿入光源のずれは未確認。脱偏極率の追跡との比は 4.2 GeV・12 MV の 1 点 |

### (f) v2.1.0 の (d) 「要判断・未解決」のうち v2.2.0 で解決したものの扱い

下の v2.1.0 の (d) は変えていない。次の項目は v2.2.0 で解決・更新した (残りは上の (d) に引き継いだものか、v2.1.0 のまま)。

| v2.1.0 の (d) の項目 | v2.2.0 での扱い |
|---|---|
| tradke の符号 (B111)・RADTAPER の追跡の重み | **確定**: 符号は放射反作用つきのローレンツ力の積分と 6 桁一致 (k64 が誤り)。RADTAPER は (δ_l + δ_l+1)/2 が物理的 (v2 は Emittance[] も追跡もこれ)。v2 の Emittance[] が gs とずれるのは K64 の誤差。SOL 領域内の tradke (bzh ≠ 0) の追跡との整合・SP8 (大きなリング B) の U0/E の確認は B で RADTAPER の閉じ方だけ確認 |
| v2 の Emittance[] と RADINT の差 (矩形磁石 1.4e-3) | **解決**: v2 の誤りではなく EMITDIV の既定の離散化 (+1.4e-3)。収束値は RADINT と 1.3e-6。k64 の「1e-4 の整合」は +1.3e-3 と −1.3e-3 (端の項 B75) の偶然の相殺。(d) の EMITDIV の既定が残る |
| v2 の Polarization の脱偏極率の過小評価 (0.43 倍) | **解決**: 原因はベータトロンのチューンの δ 依存。4.2 GeV・12 MV で 1.05 ± 0.16。追跡の測定が 1 点なのが残る |
| n が傾く・非平面のリングの緩和率 (B116) | **解決**: 厳密式・運動論的偏極・脱偏極率の刻みごとの積分 |
| B108 の残り (K1≠0・K0≠0 の ANGLE=0 BEND・MULT の AE、追跡の端の放射、傾き 30 mrad の 0.02〜0.03) | **解決・切り分け**: AE を実装、追跡・MULT の放射を修正 (B143)。30 mrad は 2.7% の系統差と判明 (~~原因未特定~~ → v2.4.0 D2: 基準側の四極のフリンジ・U0 の定義の差と考えられる、B222) |
| 大きなリング A に鉛直の四極の誤差 0.2・0.4 mm で光学が不安定 | **扱いを直した**: 不安定なリングは P = NaN (以前は 0.92 などを返した) |
| SynchroBetaEmittance (B128、40 倍) | **解決**: 原因を特定して v2 で修正 (始点の分散)。共鳴の効果の追跡との比較が残る |
| WSPAC の追跡のビーム行列 | **更新**: gs の追跡のビーム行列には MINCOUP の下限が入らない。v2 のオプション `WSPACSC = 1` |
| RADCOD の 1 周の解析 T と数値微分の差 2e-3 | **解決 (一部)**: B138 の 4 点を直し 2e-3 → BEND 1 台 1.7e-10・24 セル 1.9e-7。残りは B139 (K64 の仕様) |
| K1≠0 の BEND の端の模型 | **オプションを実装**: `face` (既定は `radial`) |
| 放射の β の因子 (B112)・陽子の TouschekLifetime 等 | **確認 (直さない)**: k64 は gs と全桁。1/β³。要判断 |
| 領域内の QUAD・MULT のスピン (bz ≥ 3)・領域内の光子・領域内の CAVI の POL | **解決**: スピン 9e-7 以下 (B146)、光子は DRIFT・QUAD・MULT、CAVI は VOLT = 0・0.2 MV で追跡と 1.6e-7。領域内の加速は実装しない (判断) |
| DETR ≥ 1 の INS | **確認**: 3 つのリングで β・DETR・Emittance[]・U0 が gs と一致、位相の整数部は B153。MARK の生の DETR は B159 (未解決) |
| MEASURE・TRCOD の gs との比較 | **更新**: gs は先頭の位置だけ動く。k64 を gs と同じ統計で実装。v2 の追跡経由は比べられない |
| TrackParticles の F1 (`tsteer`・`tmulti`)・MULT の AE | **解決**: 実装 (k64 は gs と全桁) |
| PSPAC (K64 の形状関数)・`DAPWIDTH` の窓・先頭でない成分の FREE・複素の特殊関数の未対応の関数 | **解決**: k64 の PSPAC、trackd の直訳、gs と全桁の何もしない変数、15 関数。複素の特殊関数の k64 の値は (d) に残る |
| `45DEG`・`SIN[0.1]`・`1+`・MAIN の `E`・`Pi` | **解決**: (a) の文法の行。v2 の `45DEG`・`1+` は (d) |
| `OYStability`・`OYScan`・Micado の変種・FITP>1・Tk 9.0 | **確認**: (b) の最後の行 |

## v2.1.0 (2026-11-18)

v2.0.0 の「未対応・未確認」を追求して直した版 (記録は `history/v2next_A1.md`・`A2`・`B1`・`B2`・`B3`・`C`・`D1`・`D2`・`E`、例は [`example/v2_next/`](../example/v2_next/))。互換モード ([`pysad/compat.py`](../pysad/compat.py)、既定 v2、`PYSAD_COMPAT=k64` で K64 の動作) の使い方は下の v2.0.0 の (a) と同じ。**v2.0.0 の節は変えずに残す**。この節の (a)〜(f) は v2.1.0 だけの一覧で、v2.0.0 の (a)〜(f) とは別。

- 使ったリング: 小さなリング (12 セル・24 セルの FODO、周長 41.8〜105.6 m、0.2〜4.2 GeV)、大きなリング A ([`pysad/tests/kek-ls_v03_12.lat`](../pysad/tests/kek-ls_v03_12.lat)、3 GeV、周長 603 m、3563 要素)、PF リング ([`example/pf/lattice/v9_16_asym.sad`](../example/pf/lattice/v9_16_asym.sad)、2.5 GeV、187 m)、ASTRID2 ([`example/astrid/lattice/ASTRID2_v01_01.sad`](../example/astrid/lattice/ASTRID2_v01_01.sad)、MPW 入りは `ASTRID2_v01_02_id.sad`、0.58 GeV、45.7 m)。SP8 (大きなリング B、8 GeV、周長 1436 m、4 回対称、矩形の偏向磁石 88 台) は数値だけ確認に使った (A2 の RADINT・Polarization、D1・D2 の NORFSW・COD・Emittance[])。
- K64 の側のバグの番号は B111〜B134 ([`21_sad_bug_report.md`](21_sad_bug_report.md))。記録の項目名との対応は `local/wake/v21_bnum.md`。
- 根拠の種類は v2.0.0 と同じ (解析式、独立な数値計算、gs)。gs は 2026-11-17〜18 の 10-19 ビルド。gs で確かめられない v2 の値は「v2 の値は独立な根拠だけ」と書いた。

### (a) v2.1.0 で結果が変わるもの (K64 = k64 モードの動作との違い)

「K64」は gs の動作で、k64 モードの動作でもある。k64 モードを変えないのが原則 (v2.1.0 で k64 モードの結果を変えたのは、pysad 自身の不具合の修正 (下の (c)) だけ)。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B111 | 放射の反作用の角度の符号 (tradke と tradkn)・RADTAPER の追跡の BEND の重み | tradke (Emittance[]・閉軌道) の `ddpx` の符号が tradkn (追跡) と逆。Emittance[] の閉軌道から追跡を 1 周すると x −4.3e-7 m、z 2.4e-6 m 動く (小リング)。RADTAPER の追跡の BEND の重みは (δ_l + δ_l+1)/2 | 符号を反転。RADTAPER の追跡の重みを (3δ_l + δ_l+1)/4 (tturne と同じ) | 閉軌道の 1 周の追跡の動き (v2 / k64): 小リング x 1.4e-9 / −4.35e-7 m、大きなリング A (RADCOD) −4.5e-11 / −2.0e-8 m、PF 2.7e-11 / −6.2e-9 m、ASTRID2 (MPW 入り) −5.4e-9 / −1.17e-7 m。**v2 の Emittance[] (RADCOD) は gs と相対 1e-4〜1e-3 ずれる** (k64 は全桁一致)。符号の根拠は導出と「閉軌道が閉じる」ことで、K64 の作者の意図の裏づけは無い | `test_v2_A1_next.py::test_radcod_closed_orbit_closes_in_tracking_v2`・`…_does_not_close_k64`・`test_tradke_python_and_numba_agree_v2` |
| B112 | 放射の古典半径 (MASS≠電子) | `re0 = elradi` (質量によらない)。30 GeV/c の陽子の U0 が 1.2651 eV | `re0 = elradi·elmass/amass`。U0 = 6.88e-4 eV (古典式) | 古典式 U0 = (2/3) r_p mc² β³γ⁴ Σθ²/L と追跡の 1 周の δ の変化が 1e-6 (質量 1×・2× 陽子)。**電子の放射は両モードでビット単位で同じ** | `test_v2_A1_next.py::test_proton_radiation_loss_v2`・`…_k64`・`test_electron_radiation_same_in_both_modes` |
| B114 | ずれた BEND の CALC (tbrote の出口) | 出口のドリフトが逆向き。矩形磁石 1 台の DX = 2 mm でチューンが 2.7e-5 動く (追跡は動かない) | 出口も入口と同じ回し方 | 厳密な幾何 (弦に垂直なずれの帯の中の円弧) と CALC が 2e-8。矩形は光学が不変 | [`test_v2_A2_dx.py`](../pysad/tests/test_v2_A2_dx.py) (10 件) |
| B115 | RADINT の ANGLE=0・K0≠0 の BEND (挿入光源の極) | 無視 (ASTRID2 + MPW で U0 6.17 keV) | K0 を曲げ角として積分に入れる (U0 7.81 keV) | I2 = Σ(ANGLE+K0)²/L と 1e-6、U0 は Emittance[] と 1.1e-3。**Jx 0.7%・εx 3.5% は Emittance[] とずれる** | `test_v2_A2_radint.py::test_radint_includes_k0_poles_of_insertion_device_v2` |
| B116 | Polarization・PolarizationTime (逆向きの BEND = 挿入光源の極) | 符号つきの率 Σ r (b·n) で極が相殺 (ASTRID2 + MPW で P = 92.376%、τ = 28842 s) | 符号に依らない Σ r (b·n)² (P = 0.70190、τ = 21915 s) | 独立な式 P = 0.70163、τ = 21886 s と 2e-3。K64 の追跡の緩和率の比 (MPW あり/なし 1.311) が v2 1.315、式 1.317 | [`test_v2_A2_pol.py`](../pysad/tests/test_v2_A2_pol.py) (7 件) |
| B117 | K1≠0 の BEND の追跡の端 (tbendi) | 薄い端の写像。CALC と x-px の行列が 1.1e-5 違う (x0 = 1 mm) | `DISFRIN=0` かつ端の角 ≠ 0 で CALC と同じ厳密な tbedge | 厳密な写像 (ローレンツ力の数値積分) と 1e-10 (K64 は 1e-10〜9e-9)。CALC との行列の差 4e-8 | `test_v2_B1_optics.py::test_k1_bend_*` |
| B121・B122 | `Union[…, SameTest->f]`・`WriteBeamLine[…, Name->"L2"]` (文字列) | 先頭の要素だけを返す・`Insufficient number of Names` | 真値 (`{1,3,5}`)・文字列の Name も受ける | gs の実測と定義 | `test_v2_B2_script.py::test_union_samplest_true_value_in_v2`・`test_writebeamline_name_string_k64_vs_v2` |
| B123 | 光子の偏光 ξ1 (tphrec) | pp で 2 回割り \|ξ1\| が 1e3〜1e6 | pp で割らない (\|ξ1\|, \|ξ3\| ≤ 1) | ξ1² + ξ2² + ξ3² ≈ 1。k64 は gs と全列一致 | [`test_v2_B3_photons.py`](../pysad/tests/test_v2_B3_photons.py) (`test_stokes_*`) |
| B125 | `Polarization` の脱偏極率 (spdepol) | emittance/減衰率に比例 (4.2 GeV・12 MV の FODO で P = 1.27e-4。独立な率の約 1000 倍) | 側帯 (シンクロトロン振動によるスピンの位相の変調) を厳密に取った独立な式。`depol_rate` (側帯込み)・`depol_rate_linear` | 線形の軌道 + 厳密なスピン回転の MC と 10% 以内。FLUC つきの追跡には 0.43 倍 (過小評価。(d) 参照) | [`test_v2_B3_pol.py`](../pysad/tests/test_v2_B3_pol.py) (7 件)、`test_v2_D1_emit.py::test_k64_spdepol_rate_scales_as_emittance_over_damping_v2` |
| B107・B126・B109・B106 | RADINT の BMAG の CKL1Y・QUADRU の軌道のずれ、tbgrote の (5,6)、tradke の b | CKL1Y = −cos θ·DKL1、QUADRU は四極自身の DX・DY を引かず K1<0 で符号が逆、(5,6) に px² + py² が残る、b = bzh·0.5 | CKL1Y = +sin θ·DKL1、(dx, dy) = 閉軌道 − 要素のずれ、コメントアウトされた式、b = bzh | I4U + I4V が物理的な I4 と 1e-7 (K64 は 1.6〜2.7% ずれる)。ΔJx = +0.0078 (追跡 +0.00781、K64 +0.0028)。(5,6) は有限差分と 4e-12。**B106 の数値の効果は 1e-13 以下** | [`test_v2_C3_radint.py`](../pysad/tests/test_v2_C3_radint.py) (21 件)、[`test_v2_C3_misc.py`](../pysad/tests/test_v2_C3_misc.py) |
| B108 | ANGLE=0 の BEND (補正磁石・挿入光源の極) の端の角 AE1・AE2 の放射 | 端の項が約 2 倍 (B75 と同じ) | 行路を x tan(e) 短くする (ANGLE≠0 の `tbende` と同じ) | 厳密な幾何 + ロビンソンの定理: ΔJx が 厳密 : v2 : k64 = 1 : 1.03 : 2.08 (傾き 5 mrad の極)。傾き 30 mrad では v2 は 0.02〜0.03 ずれる | `test_v2_E_steer.py::test_steering_edge_angle_*` |
| B127 | WSPAC の Emittance[]・EMIT のチューンシフト | 最後のスライスの `gammab(nlat+1) = 0` でレンズが 8 倍強く、シフトが 3〜4% 強い | p1 = (1+δ) p0 | 一次のラスレットの積分と 0.2〜1.8% (大きなリング A・PF・ASTRID2。K64 は 0.9〜7%)。k64 は gs と全桁一致 | `test_v2_D1_emit.py::test_wspac_*` (4 件) |
| B128 | SynchroBetaEmittance の収束の警告 | 警告なし (大きなリング A で 40 倍違う値を返す) | **値は変えず**、AzimuthalModes m と 1.5 倍の m の emitx が 20% 以上違えば警告 | 大きなリング A で m = 9 が 5.1e-9、m = 14 が 2.5e-9 → 警告。B・ASTRID2 は警告なし | [`test_v2_D1_emit.py`](../pysad/tests/test_v2_D1_emit.py) (警告 2 件。1 件は `PYSAD_SLOW=1`) |
| B129 | CorrectOrbit の `Solver->{"Micado",n}` | 全補正量の最小二乗の先頭 nx 個を使い、補正が悪化する | 選んだ nx 本だけを解く | 選んだ列の `lstsq` と 1e-9。24 セルで n = 1 から 24 まで単調に下がる (4.74e-5 → 2.06e-5 m) | `test_v2_D2_packages.py::test_micado_*` (4 件) |
| B130 | Emittance[]・EMIT の POL の量 (傾いた・ずれたソレノイド) | 縦磁場を z 軸の向きに足す (SpinTune が追跡と 1.5e-4〜1.1e-3、n0 が 2.5e-3〜1.5e-2 rad 違う) | B⊥ と B∥ を運動量の向きに分けて足す。領域内の DRIFT・QUAD の放射・スピンの刻みも bz·L/n ≤ 0.02。向きの変化の大きさを sin θ でなく θ に (領域の外にも効く: N8、平面リングで SpinTune が 2.6e-5 変わる) | 回転した縦磁場の歳差の閉じた式と 5e-7。SpinTune の差 1e-13〜5e-7、n0 の差 0〜1e-5 rad (大きなリング A・PF・ASTRID2・ASTRID2 + MPW でも) | [`test_v2_E_spinsol.py`](../pysad/tests/test_v2_E_spinsol.py) (20 件) |
| B131 | SOL 領域内の CAVI のスピン (追跡) | 歳差しない | 放射・スピンの刻みに通す | Thomas-BMT と 1e-9 台 | [`test_v2_E_spinsol.py`](../pysad/tests/test_v2_E_spinsol.py) |
| B132 | PhotonList のエネルギー列 | γ でなく βγ を掛ける (3 GeV で 1.45e-8) | γ を掛ける | gs と E 列だけ 1.4e-8 違う (k64 は全列一致) | `test_v2_D1_emit.py::test_photons_*` |
| B133 | WSPAC ON の CALC・CALC6D | NaN (gs)。k64 は従来どおり空間電荷を入れない | 空間電荷を入れる (6D の経路。RFSW ON の 4D の CALC も) | 一次のラスレットの積分と 0.15%、Emittance[] のチューンと 2e-6。K64 の `qwspac` の直訳ではない (gs で比べられない) | `test_v2_B1_optics.py::test_wspac_*` |
| B83 (N8) | 平面リングのスピン角 (example/07 の 500 周) | 刻みの誤差が積算 (0.8794 rad) | 厳密な θ (1.1516 rad) | 刻みを EPS=1E-5 に細かくすると gs と v2 が 1e-8 で一致。**保存値 0.8794 rad は K64 の刻みの誤差で物理の値ではない** | [`test_v2_A2_spin.py`](../pysad/tests/test_v2_A2_spin.py) (9 件) |
| (B1) | SOL 領域内の DRIFT のスピン・放射の刻み | 1 刻み (bz = 9 で放射が 1.7 倍外れ、スピンが 1e-2) | ndiv = 1 + int(\|bz·al\|/0.02) | Thomas-BMT の積分と 2.3e-7〜8.9e-7、らせんの放射損失と 1e-4 (3 GeV、bz = 3・9) | `test_v2_B1_sol.py::test_spin_in_sol_drift_at_large_bz_matches_thomas_bmt`・`test_sol_region_drift_radiation_*` |
| (A1) | PSPAC のメッシュ | (PSPAC は v2 だけ) | 結果は変えず、粗いメッシュの警告 (`???-PSPAC-Warning-…`) | 横の rms が 1.5 セル未満・縦の rms が 1.5 セル未満・粒子の 1% 以上がメッシュの外。メッシュが十分細かいと解析式と 1〜2%、gs と 1% | `test_v2_A1_next.py::test_pspac_resolution_vs_analytic_and_gs` |

- B99 (DELCOR の規格化)・B40 (CAVI の符号) は v2.0.0 のまま。規格化の根拠 (K64 は複素数の ⟨\|c\|²⟩ = 1) と、CAVI の符号は追跡に効かず観測できないことが確定した (`history/v2next_C.md` §3.4)。
- B83 の行の根拠の数値は 24 セルの FODO (3 GeV、RF 4 MV)。k64 の 1 周のスピン角 −1.1940907 rad、v2 −1.1935299 rad、EPS=1E-5 で gs −1.1935395510、v2 −1.1935395636。

### (b) モードによらない追加・修正と、新しく使える機能

「モード」の列: **v2** = v2 モードだけ、**両** = 両モードで同じ (pysad の追加・修正、または gs と同じ結果が出る機能)。例は [`example/v2_next/`](../example/v2_next/) (各担当の README・`.gs.txt` が gs、`.pysad.txt` が v2、`.pysad_k64.txt` が k64)。

| 機能 | モード | 使える範囲と条件 | 使えない・未確認 | 例・記録 |
|---|---|---|---|---|
| **NAFF** (チューンの抽出。[`pysad/naff.py`](../pysad/naff.py): `naff`・`naff_real`・`tune_xp`・`tunes_from_coords`) | 両 (新しい関数) | 放射を切った追跡の座標列 (1/周回の単位)。周波数の間隔が約 2/N 以上の成分。既知の信号で周波数 1e-11。追跡の座標列 (1024 周) のチューンが CALC と 12 セル FODO 3e-9、大きなリング A 6e-10・8e-10、PF 1.4e-10、ASTRID2 5e-11。MAIN の TRACK の `wtunem` (FFT の放物線補間) は 4e-5〜2e-4 | 放射減衰のある列 (包絡 exp(−λn) で周波数がずれる)、共鳴線の近傍 (間隔 < 1/N)、6 極つきで振幅が大きいと振幅依存のチューン。**MAIN の TRACK の FOURIER (wjfit・wfres・wfbin・wiord) は未翻訳** (gs は SIGSEGV で比べられない。pysad は FOURIER を無視する) | `B3_naff_tunes.py`、B3 §6 |
| **光子 (PhotonList)** (RAD・FLUC・PHOTONS ON の TrackParticles) | 両 (xi1 の誤りと E 列は v2 だけ直す) | BEND (K1=0・K1≠0・ROTATE・DX・DY・DROTATE)・QUAD・MULT (K0・K1・K2・DX・DY・ROTATE)・SEXT・OCT (DX・ROTATE)・磁場つきドリフト・BEAMBEAM + BEND (TRPT つきでも)。k64 は gs と光子 358 個・325 個・71 個が全列一致 (E 列は B132)。v2 は光子数が解析式の 4σ 以内、⟨u⟩/u_c = 8/(15√3)、スペクトルの KS 検定。条件: 3 GeV の 12・24 セル FODO (1 粒子 1 周 1 ms 程度)、1.7 GeV (BEAMBEAM) | **SOL 領域内の要素の光子は記録しない** (乱数だけ K64 と同じに引く)。COORD・TCAVI・MAP があるラティスは従来の高速追跡 (BEND がなければ光子なし)。POL の追跡との併用は未確認。QUAD の光子の位置は gs と未照合 (BEND は全桁)。BEND があると Python の追跡になり遅い | `B3_photons_bend.sad`、`D1_photons_mult.sad`、`D1_photons_bb_bend.sad`、docs/45 |
| **WSPAC の Emittance[]・EMIT・CALC6D・RFSW ON の CALC** | Emittance[]・EMIT は両、CALC・CALC6D は v2 | 反復して空間電荷つきの一周の写像からチューン・エミッタンス・Twiss を作る (INTRA とも併用可)。0.2 GeV の 12 セル FODO (PBUNCH 1e5〜1e8、MINCOUP 0.5・0.2)。**k64 は gs と全桁一致** (小リング、大きなリング A・PF・ASTRID2 で、PBUNCH 1e9〜1e10)。v2 は一次のラスレットと 0.2〜1.8% | 追跡の WSPAC のビーム行列は空間電荷なし・MINCOUP なしの平衡値のまま (`LINE["SIG…"]`)。強い空間電荷 (0.2 GeV で PBUNCH ≳ 1e9 は反復が収束しない)、RF なしの CALC、TouschekLifetime との併用、実機の大きなリング | `B1_wspac_emittance.sad`、[`test_v2_B1_optics.py`](../pysad/tests/test_v2_B1_optics.py)、[`test_v2_D1_emit.py`](../pysad/tests/test_v2_D1_emit.py)、B1 項目 3、D1 項目 2 |
| **成分ごとの FREE** (`FREE QF.3 K1`・`FIX B.3` など) | 両 (pysad の修正) | 語に `.` があれば成分の変数。既定のキー (QUAD の K1 など) は先頭の成分に結合 (`???-FFS-Info-Component QF.2 is coupled to QF.1 .`)、既定でないキー (BEND の K1・ROTATE) はその成分だけ。CALC・GO・VAR・FIX・NOKEEPEXP。小さな FODO で反復の残差・NX が gs と全桁。大きなリング A (`FREE BPI6.2 K1; FREE BPI7.5 K1`) でも反復の行・チューン・K1 が gs と一致 | **先頭でない成分の既定のキーの変数は作らない** (gs は値を持たない変数として数える: 反復の行の変数の数が 1 つ違い、その変数だけで GO すると pysad は `?No variable.`)。COUPLE で結合した要素の成分、`USE NOEXPAND` の展開済みの成分は未確認 | `B2_free_component.sad`、[`test_v2_B2_free.py`](../pysad/tests/test_v2_B2_free.py) |
| **複素の PolyLog・HypergeometricU・Zeta・HurwitzZeta** ([`pysad/specfn_cx.py`](../pysad/specfn_cx.py)) | 両 (**k64 も真値**) | 複素の z・s。Zeta・HurwitzZeta は \|s\| ≤ 500、PolyLog は実数の次数・\|z\| ≤ 1e8、HypergeometricU は \|z\| ≤ 30 (Re a > 0 は積分表示、a が負の整数は多項式、ほかは漸化式と Kummer の変換)。独立な参照 (Borwein、Euler–Maclaurin、spence、積分表示、接続公式) と 1e-10〜1e-11。gs の誤り (B118〜B120) の点でも真値。負の実数の z は上側 | **k64 が gs の誤った値を再現しない** (再現には K64 の複素の不完全ガンマなどの直訳が要る。要判断)。複素の次数の PolyLog、Re a ≤ 0 の複素の a・b の HypergeometricU (SADError)、z = 0 の U、\|z\| が範囲外 | `B2_special_functions_complex.sad`、[`test_v2_B2_specfn_cx.py`](../pysad/tests/test_v2_B2_specfn_cx.py) |
| **文脈・`Symbol`・`Message`・VISIT/BYE・`Ffs;`** | 両 (gs に合わせた) | `$Context` への文脈の記号 (名前の末尾にバッククォートを付けたもの) の代入と `$ContextPath` への代入が `Begin` と同じ規則で働く、`WithContext`、入れ子の文脈 (台本の 19 行が gs と全行一致)。大域の名前 51 個を最初から持つ。`Symbol["x"]` は値を返す、`Message[f::tag, 値…]` は `???f::tag:  本文`、`General::tag` へ。VISIT・BYE・USE が FitFunction・FitValue・ElementValues・OpticsProlog・OpticsEpilog を退避・復元 (8 行が gs と一致)。`FFS` と `ffs` だけが FFS の語 (`Ffs;` は構文エラー)。LINE 名は書かれたとおり。**`Begin` に文脈を文字列で渡す形は k64 が拒否、v2 は受ける** | 文脈つきの `Names` (規則が決められない)、init.n が大域に作る名前の全部 (51 個しか調べていない)、`Off[f::tag]`、`BYE ALL` の gs との比較 (gs は SIGSEGV)、VISIT の未定義の LINE のメッセージ | `B2_context.sad`・`B2_symbol_message.sad`・`B2_visit_bye.sad`、[`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py) |
| **`Parallelize`** ([`pysad/fileops.py`](../pysad/fileops.py)) | 両 | POSIX (Linux・WSL・macOS) は npara 個まで fork で同時に計算 (結果の並びは gs と同じ。子の代入は親に戻らない。lshare 超過は `Shared::toolarge` で `$Failed`、Timeout 秒で子を殺す)。WSL で 4 つの `Sleep[1]` が 1.02 秒。`PYSAD_PARALLELIZE=serial` で逐次。親でカーネルを一度走らせてコンパイル・キャッシュを済ませてから呼ぶ | **Windows は逐次** (spawn では Interp・FFS・ラティスを pickle できない)。numba と fork の併用は実際には試せていない (WSL に numba が無い) | `B2_parallelize.sad`、[`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py) |
| **`TrackParticles` の MASS・CHARGE** | 両 (放射の古典半径は v2 だけ: B112) | FFS の中で与えた陽子・重い粒子・多価イオンの TrackParticles・Emittance[]・CALC6D・RF の位相・ソレノイド。k64 は gs と全桁一致 (30 GeV/c の陽子、大きなリング A の 3 GeV/c の陽子)。ドリフトの z、ソレノイドのラーモア回転 (電荷 1・−1・2)、Qs が解析式と 1e-3〜1e-4 | `CHARGE` は FFS に入るたびに 1 に戻る (gs と同じ)。放射の β の因子 (β < 0.99 では U0 の 1/β² 倍)、TouschekLifetime・Intrabeam・RADINT の古典半径は電子の値のまま (陽子での確認なし) | `A1_proton_ring.sad`、[`test_v2_A1_next.py`](../pysad/tests/test_v2_A1_next.py) |
| **Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE**、**ANGLE=0 の BEND (補正磁石) と MULT の F1** | 両 (K64 の動作の追加) | k64 は gs と全桁一致: BEND (K1 の有無、FRINGE = 1・2・3・−1・−2、F1・FB1・FB2、K0、DX・DY・ROTATE) を小リング 9 通り、大きなリング A・PF・ASTRID2・MPW 入りの全偏向磁石に F1 を付けて。ステアリング (L = 0.2 m) と MULT (K1 = 0.7・−0.6) の F1 は FODO 14 通りと、大きなリング A・PF・ASTRID2 に 2 台挿入して gs と全桁。CALC のチューン = Emittance[] のチューン (RF・RAD なし、1e-9)。F1 のある BEND は numba の `k_bend` を使わず Python の経路 | **TrackParticles の F1 (ANGLE=0 の BEND の `tsteer`、MULT の `tmulti`) は未対応**。MULT の AE1・AE2 は読まない。`tbfrmle` の行列は数値微分。逆向きの要素 (`direlc<0`) の F1 の入れ替え。ステアリングに ROTATE と DX/DY を同時に付けると gs と CALC のチューンが 1e-4 違う (F1 に関係なく既存、未調査)。MULT の FB1・FB2 の `tblfre` (`multae.py`) は NotImplementedError のまま | `B1_bend_f1_emittance.sad`、`E_steer_f1.sad`、`E_mult_f1.sad`、[`test_v2_B1_optics.py`](../pysad/tests/test_v2_B1_optics.py)・[`test_v2_E_steer.py`](../pysad/tests/test_v2_E_steer.py) |
| **TrackParticles の BEND の F1・FRINGE (B113)・CAVI の DX・DY と V1 V20 V11 V02** | 両 (gate なし) | `track6n.bend_links` が tfgeo の inext・iprev を再現。CAVI は TENT・TEXIT と V1 V20 V11 V02。gs と全桁一致 (NORAD・RAD、POL の Python の追跡も)。大きなリング A (空洞に DX = 2e-4・DY = −3e-4、QDA2 を 0.03 rad 回して RAD の 5 周追跡) で差 0 | WakeFunction のある CAVI・SOL 領域内の CAVI の V1 V20 V11 V02 (NotImplementedError)、Emittance[] の CAVI の V1 V20 V11 V02 (`tcave` は V = 0 を前提) | `A1_bend_f1_adjacent.sad`・`A1_cavity_offset.sad`、[`test_v2_A1_next.py`](../pysad/tests/test_v2_A1_next.py) |
| **EMIT のあとの `USE`** | 両 (K64 の仕様) | `trf0` を `USE` で持ち越す (EMIT・Emittance[]・CALC6D のあとの追跡の RF の位相が EMIT の値のまま)。gs と全桁 | `emit`・Touschek の結果は USE で捨てる (K64 での残り方は未確認) | `A1_use_after_emit.sad`、`test_v2_A1_next.py::test_use_after_emit_keeps_trf0` |
| **CALC6D・RADCOD つきの CALC・INTRA と SOL 領域**、**`DP0` の別名** | 両 (gate なし) | SOL 領域のあるリングの CALC (4D・6D・RADCOD)・Emittance[]・EMIT・INTRA。k64 は gs と 1e-9 (SOL 領域のあるリング、3 GeV、BZ = 2 T)。CELL CALC6D のあと `DP0` を先頭の `DDP` にする。CALC6D + RADTAPER は k64 もテーパーを入れる (gs と 8e-8) | TouschekLifetime は SOL のあるリングで未確認 (gs が NaN)。BOUND=0 の SOL (BZ≠0) の CALC は NotImplementedError。RADTAPER + SOL・WSPAC・INTRA + SOL の併用は走るが値の確認なし。**4 次元の CALC は RADTAPER を見ない (K64 の仕様。実装しない)** | `B1_sol_calc6d_intra.sad`、[`test_v2_B1_optics.py`](../pysad/tests/test_v2_B1_optics.py) |
| **傾いたソレノイドの POL** (B130 の修正の機能面) | v2 | v2 の Emittance[]・EMIT の POL (SpinTune・n0・Polarization・脱偏極ベクトル) が、ずれ・傾き・F1 つきの SOL 領域 (領域内は DRIFT・QUAD・BEND・MULT・CAVI・RF つき MULT・SOL) で追跡と一致。1〜3 T・0.3〜3 GeV | bz ≥ 3 1/m の QUAD・MULT の D は 1e-4 台。領域内の QUAD の端の ψ は hard edge の Bs なし (軌道が数 mm の範囲で確認、cm のずれは未確認)。領域の外の B∥ (垂直な軌道の傾きの中の偏向磁石の B_y y') は K64 も追跡も入れていない | `E_sol_tilt_spin.py`、[`test_v2_E_spinsol.py`](../pysad/tests/test_v2_E_spinsol.py) |
| **POL + TRPT + 加速** (スピン追跡 track6 の VOLT≠0 の CAVI) | v2 | SOL 領域の外の CAVI。DX・DY・ROTATE も。断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 と 5e-7 (1 GeV、扇形の偏向磁石の間に 0・0.2・0.5 GV)。track6 の `_ndivrad` が TRPT の 10 倍を見ていなかったのを直した (POL・RADPOL・光子の追跡の TRPT の結果が gs に近づく) | SOL 領域内の加速、RF つき MULT・TCAVI の加速、リングの何周もの TRPT (K64 の TRPT は輸送路用)、スピン追跡の損失の処理 (これらは SADError)。k64 は従来どおり明示のエラー | `D1_spin_accel.py`、`test_v2_D1_emit.py::test_spin_tracking_with_trpt_*` |
| **領域内の四極の DX・DY と幾何の fit** | 両 (領域内の四極のずれは k64 も gs と同じ) | 追跡・`solemit`・`sgeo` で領域内の四極の DX・DY (gs と LINE["GX"…] が全桁一致)。幾何の fit (QF の DX、QD の DY) が gs と同じ 2 回の Newton で収束 (変数 1e-12)。解析微分は幾何の有限差分と 2e-3 | SEXT・OCT・MULT の領域内のずれは未確認。四極の F1・FRINGE と領域内の BEND・MULT のずれは未対応。領域内の四極のずれの POL の追跡 (track6) は NotImplementedError。k64 は SOL 自体のずれが未対応。CHI3 の fit は K64 でも収束しない (B124) | `B3_geofit_quad_dx.sad`、[`test_v2_B3_geofit.py`](../pysad/tests/test_v2_B3_geofit.py) |
| **`NPARA=n;`** (FLUC の追跡の子プロセスへの分け方) | 両 (既定は NPARA=1 のまま) | gs の既定 NPARA=4 と同じ分け方 (`TrackParticles` は np ≥ 13 で粒子を分け、子の乱数の種は親の種 + 2k)。RAD + FLUC の追跡だけ (POL・RADLIGHT・PHOTONS・WAKE は分けない)。NPARA=1 と 4 の rms dp・最終の乱数の種が gs と全桁一致 (FODO、np = 200、1・10・100 周。PF リング、np = 40) | v2 は A1・A2 の放射の式の修正で rms が 6e-6 (1 周)〜1e-3 (100 周) 違う (乱数の引き方・最終の種は同じ) | `C3_fluc_nparallel.sad`、`test_v2_C3_misc.py::test_fluc_*` |
| **MAIN の要素定義の式** (`ANGLE=PI/8`・`SIN(...)`・`10 DEG`) | 両 (pysad の不足) | gs の文法 (約 100 通りの実測): `PI`・`SIN COS TAN ATAN SQRT EXP LOG LN SINH COSH TANH ATAN2`・`+ - * / ^`・単位は数の直後だけ。大きなリング A・B、PF、ASTRID2 の全要素の値が v2.0.0 の読み取りとビット一致 | gs が構文エラーにするが pysad が受ける形 (`45DEG`、`SIN[0.1]`、未定義の名前)。gs が受けるが pysad が構文エラーにする形 (末尾の `1+`)。MAIN の変数を `E`・`Pi`・`Degree` で代入しても要素定義の式の `E` は SADScript の定数 | `D2_main_expr.sad`、[`test_v2_D2_parser.py`](../pysad/tests/test_v2_D2_parser.py) |
| **Packages の関数**: `Log2$`・`FunctionDefinition`・`Fill$String`・`Remove$Paren`・`Emittance->False`・`TURNS` の既定・`Element[..]+=` | 両 (pysad の不足) | 小さな台本を gs と行ごとに突き合わせた (functions.n・Misc.n・element.n・beamline.n・emit.n・dapert.n)。`Element["K1","QF"]+=0.1`・`LINE[..]+=` がフックを通らない不具合を直した。`SetElement` の戻り値、`Touschek` のメッセージ、`DynamicApertureSurvey` の TURNS の既定 −1 | `DAPWIDTH` の既定 (gs は 7 で窓の追跡、pysad は 51 で全点)、`Axes->"XZ"` の z が 3 点のとき、`Names`、`PROFILE` の代入、`Through` の文言、除算 `a/b` が gs では `a*(1/b)` (最後の桁の 1 ulp) | `B2_packages_compare.sad`、[`test_v2_B2_script.py`](../pysad/tests/test_v2_B2_script.py) |
| **PeriodicR・CorrectLocal** | 両 | 4x4 の一周行列 6 通りで R・A・B の 12 成分が gs と 3e-16。CorrectLocal (12 セル FODO、領域 {BA3, BA8}、垂直) の補正量・BPM の DY が 1e-14 | 水平・2 平面同時の Micado、`Bump`・`SVDmode->1` つきの Micado は未確認。PeriodicR は `PeriodicOptics` を先に呼ぶ (gs も同じ) | `D2_periodicr.sad`・`D2_correctlocal.sad`・`D2_micado.sad`、[`test_v2_D2_packages.py`](../pysad/tests/test_v2_D2_packages.py) |
| **`pysad.nbcache`** (numba のキャッシュの無効化) | 両 | `py -3.12 -m pysad.nbcache`。numba の索引は関数を定義したファイルの内容でだけ無効になり、呼ばれる側 (別のファイル) の変更で無効にならない。`__pycache__/nbdeps.json` に依存先のハッシュを覚えて、変わったモジュールの `.nbi`・`.nbc` を捨てる。**kernel (`elkern.py`・`perfkern*.py`・`turnkern.py`・`track6n*.py`) を編集したら試験の前に 1 回走らせる** | 自動では呼ばない。最初の 1 回は記録だけ (`--force`) | `test_nbcache_purges_caller_when_callee_file_changes` |
| **Python 3.14 (Tk 9.0) の GUI** | 両 | [`test_tclbridge.py`](../pysad/tests/test_tclbridge.py)・[`test_plot.py`](../pysad/tests/test_plot.py)・[`test_gui_lang.py`](../pysad/tests/test_gui_lang.py) の 25 件と [`example/12_gui_tkinter`](../example/12_gui_tkinter) の 6 本が通る。Windows の Tk 9.0 で `$DefaultFixedFont` を `Courier New` にする。3.12 と 3.14 で [`test_spin.py`](../pysad/tests/test_spin.py) の値が同じ (複素数の inf の扱いの差を吸収) | 実際の画面の撮影、マウス操作、macOS・Linux。Tk 9.0 の `cget` が長さを pt で返す (`4p`) のを `.n` が数として使うと壊れる (未確認) | `D2_README.md`、D2 §2.2 |

### (c) pysad 自身の不具合の修正 (K64 のバグではない。B 番号を付けていない)

| 項目 | 不具合 | 修正 | モード・根拠 |
|---|---|---|---|
| `Twiss["DETR"]` が INS の CALC の後で 0 (N12) | `ffscmd._rows_from_tws` が `detr` を渡していなかった | 渡す | 両。小さな結合リング・大きなリング A・PF リングで gs と全桁 (DETR < 1)。**DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認** |
| `TrackParticles` が MASS・CHARGE を渡さない | 電子の質量・電荷 1 のまま追跡 | 渡す (`ffscmd._mass_charge()`) | 両。gs と全桁 |
| BEND の F1 の放射の断片 (A1) | `first`・`last` を常に True にして、隣り合う BEND と 1 個の BEND の輸送路で F1 の断片を余計に入れた | tfgeo の inext・iprev を再現 (B113) | 両。RAD の追跡が gs と 1e-6 (8 桁表示の丸め) |
| CAVI の DX・DY・ROTATE と V1 V20 V11 V02 | 追跡の CAVI に TENT・TEXIT が無く、パーサが V1 V20 V11 V02 を捨てていた | `tent_all`・`texit_all`、`k_cav_q`、パーサ | 両。gs と全桁 |
| EMIT のあとの `USE` で `trf0` が 0 に戻る | `USE` で `self.trf0 = 0.0` | 持ち越す | 両。gs と全桁 |
| NORFSW; RAD の Emittance[] の閉軌道 (D1 項目 1) | 空洞があって RFSW OFF のとき、z の戻りを δ で閉じていた (K64 の tcod は δ も z も動かさない) | 「空洞があって RFSW OFF」のときだけ z を δ で閉じない | 両。SP8 (大きなリング B、QF5 に DX = 5e-5): gs 7.1529048e-9、k64 7.1529045e-9 (以前 6.8685e-9)。大きなリング A・PF も gs と同じ。空洞のないリングは従来どおり (回帰 3c・3d・26b) |
| `_tquade` の軌道ゼロの四極のスピン軌道結合 D が 0 (POL ON、E の (a)) | tradke を通らない | `calpol` のときは通す | 両。斜め四極つきのリングで D の 18 成分が gs と全桁 |
| 領域内の四極 (K1<0 の −90°、ROTATE) でスピンの座標系を回さない (E の (b)) | `track6._tquad_sol`・`solemit._tquade` | `_spin_rot_z`・`solemit._spin_rot` | 両 (k64 の Emittance[] の srot が gs と全桁)。Thomas-BMT と比べて K1<0 で 2.4e-3 ずれていた |
| `RadiationField` の複数区間・複数周 (C) | 区間の始めの点と直前の要素の後の点が同じで `x must be strictly increasing` | 完全に同じ点を除く (時刻が戻る点は SADError) | 両。Liénard-Wiechert の厳密な数値と振幅・位相・スペクトル・全エネルギーで一致 (gs は RadiationField で SIGSEGV: B49) |
| 不安定なリングの POL の Emittance[] が 20 分以上 (D1 項目 8) | NaN の閉軌道の 1 周 (POL では 21〜65 秒) を 30 回以上繰り返す | NaN で打ち切り、`depol_error` つきの NaN を返す | 両。SP8 (大きなリング B、四極 60 台に DY 0.1 mm) 20 分以上 → 2 秒 |
| `FREE QF.3 K1` がファミリー全体の変数になる | gs は成分の変数 | 成分ごとの FREE (上の (b)) | 両。gs と全桁 |
| `Element[..]+=`・`LINE[..]+=` が利用者の規則を作る | 複合代入がフックを通らない | 包み直した `Set` を呼ぶ | 両。gs と一致 |
| 要素定義の式 `ANGLE=PI/8` で要素を捨てる | `PI`・関数を文字列として扱い警告 | gs の文法で評価 (上の (b)) | 両。gs の 37 通りの値と 19 通りの構文エラー。**v2.0.0 の `C2` の試験の台本 `PHI=-Pi/2` は gs で構文エラーの形** で、`PHI=-PI/2` に直した |
| Python 3.14 で `test_spin.py::test_mult_ring_matches_gs` が失敗 | 不安定なリングの `_spintrm` で inf を含む複素数の演算が例外 (3.14 の複素数の扱いの変更) | `cb` が有限でなければ `nan+nanj`、`_cerfc` の `z ** 2` を `z * z` | 両。3.12 の結果と全桁同じ (k64 の値は変わらない) |
| Tk 9.0 (Windows) の固定幅フォント | `courier` が Courier New に置き換わらない | `$DefaultFixedFont` を `Courier New` に | Tk 9.0 の Windows だけ |
| `emit6d._tquade` の解析ヤコビアンの z の行 | 関数単体として z の行 (列 x px y py) が 0 | v2 は完全にした (結果は変わらない: `State.piece` が組み直していた) | v2。数値微分と 1e-10。**旧記述 (T の δ 列が 1e-5 ずれる) は再現しなかった** |
| 改行 (LF/CRLF) の整理 | 作業コピーが HEAD と違う改行の 7 ファイル | HEAD の改行に戻した | 内容は変えていない |

### (d) 要判断・未解決

各記録の「要判断」「残した項目」を、v2.1.0 で解決したものを除いて集めた。

**放射・エミッタンス・スピン**
- **tradke の符号 (B111) と RADTAPER の追跡の重み**: どちらも「閉軌道が追跡で閉じる」を根拠にした選択で、K64 の作者の意図の裏づけは無い。v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 ずれる。RADTAPER の重みのどちらが物理的かは未判断。SOL 領域内の `tradke` (bzh ≠ 0) も符号が変わるが、追跡 (`tsolqur` など) との整合は未確認。U0/E が大きいリング (8 GeV、周長 1.4 km: SP8 (大きなリング B) 相当) は未確認。
- **v2 の Emittance[] と RADINT の差 (矩形磁石)**: SP8 (大きなリング B) で v2 の εx 6.9651e-9 は RADINT の 6.9554e-9 より 1.4e-3 大きい (k64 は 1e-4 で整合)。v2 の修正 (A1・A2・C の矩形の BEND の E1・E2 と放射) の側の問題かもしれない (D2 が指摘、未調査)。
- **v2 の Polarization は脱偏極率を過小評価する**: 4.2 GeV・12 MV の 12 セル FODO で FLUC つきの追跡の 0.43 倍 (4 本・480 粒子で 3.7σ)。残りの 2 倍は式の導出ではなく軌道の非線形の効果 (確定的な蹴りへの応答で、追跡のスピン軸の応答が線形の軌道モデルの 2.07 ± 0.24 倍。δ の応答 1.9 倍、水平 1.3 倍、垂直 1.0 倍。振幅 2σ 以上、シンクロトロン + ベータトロンとともに増える)。低エネルギー (ASTRID2) では 1.0。どの非線形かは未分離。**式には入れられず**、限界の表 (E §4.4) にした。使える範囲の目安: 変調指数 aγσδ/νs ≲ 0.3。≳ 0.5 の高エネルギー・高 RF 電圧のリングでは追跡 (FLUC) で確かめる。その他: 結合項 b の u の 2 次の項を入れていない (MC に対し −10%)、全座標の応答の比は 12 試料・位置 1 か所・4.2 GeV・12 MV のみ (6 MV・24 MV・PF は未測定)、ゼロ振幅の線形応答が追跡で 13% 小さい (領域の外の偏向磁石の B∥ = B_y y' を入れると変わるか未確認)。
- **大きなリング A に鉛直の四極の誤差 rms 0.2・0.4 mm を入れると光学が不安定** (x の固有値 1.2)。側帯の項は不安定なリングでは線形の式に戻す。Polarization は 0.1 mm だけ意味がある。SP8 (大きなリング B) の POL は四極 60 台に 0.1 mm で不安定 (NaN)。
- **n が傾く・非平面のリングの緩和率** (B116: Σ r (b·n)² が最適か。ST の厳密式の (1 − 2/9 (n·v)²) の項を入れていない)。ずれのある MPW のリングの P・τ は未確認。
- **B108 の残り**: K1≠0・K0≠0 の ANGLE=0 BEND と MULT の端の角 (pysad が AE を読まない)。追跡 (`tsteer`) の端の放射は未確認。軌道の傾き 30 mrad の極で独立な値と 0.02〜0.03 違う。
- **RADINT**: K0 の極の εx は Emittance[] と 3.5%、Jx は 0.7% 違う (扇形の行列による近似)。軌道のずれ (B126) は SEXT・BEND を入れていない。結合リングで νy ≈ νs のときの J は RADINT と追跡・Emittance[] が 0.015 違う (モードが混ざり定義できない条件)。K1≠0 の BEND の傾き θ ≠ 0 の CKL1Y は数値の確認が無い。RADINT の空洞の位相は VOLT を単純に足す (v2.0.0 のまま)。
- **SynchroBetaEmittance (B128)**: 大きなリング A で 40 倍違う (K64 の式の限界か、gs・pysad k64 が共通に持つ未特定の誤りか)。原因は未特定。確かめて原因ではないもの: 閉軌道の探索、二次の分散、クロマティシティ、RF の位置の分散。SB の式の導出 (temits の `amus`・`dispp`) か、6 次元追跡の長時間の平衡 (4τ ≈ 50000 周) との比較が要る。
- **WSPAC**: 追跡の WSPAC のビーム行列は空間電荷なし・MINCOUP なしのまま。CALC (4D) の WSPAC は CALC6D と同じ 6D の経路 (K64 が NaN を返す機能なので直訳の検証ができない)。
- **その他の K64 の疑い**: CHI3 を使う幾何の fit (B124)、CALC の DISFRIN=1 の Maxwell のフリンジ (B134)、PhotonList 以外の光子 (MULT・SEXT・SOL 内・ドリフトの SOL 領域内は記録しない)、POL の追跡 (track6) の領域内の四極のずれ、K1=0 の BEND の追跡のスピン (tbend + tradk) が Emittance[]・tbendi (v2) と 2e-5 rad 違う (放射の分割が粗い K64 の精度の範囲)。
- **MAIN の TRACK の FOURIER** (wjfit の maxord > 0): gs が SIGSEGV で比べられない。NAFF だけ追加。

**追跡・要素・SOL**
- **DELCOR (B99)**: 規格化の根拠は確定 (複素数の ⟨\|c\|²⟩ = 1) だが、「rms = delx」が作者の意図かは推定のまま。CAVI の符号 (B40) は追跡に効かず観測できない。
- 領域内の QUAD・MULT のスピン (bz ≥ 3) は 1e-4 台、領域内の CAVI の POL (Emittance[]) は VOLT=0 だけ確認、SOL 領域内の光子・加速 (POL)、RF つき MULT・TCAVI の加速 (POL)。大きなリング A への SOL の確認は小リングのみ。
- **K1≠0 の BEND の端の模型** (四極が半径方向の面で始まる): CALC・追跡とも同じ K64 の仕様。端面で始める厳密な硬い端の磁石とは x0 tan e に比例する差 (x0 = 1 mm、E=0.5 で 2.3e-5)。CALC の解析行列は y ≠ 0 の軌道で自身の写像の数値微分と y = 2 mm で 1.5e-6 違う。
- **RADCOD の 1 周の解析 T と写像の数値微分の差** 2e-3 (K64 の `tradke` の微分が近似。k64・v2 とも)。BEND の T[z, δ] の近似 (`tbendebody0` の (5,6)) による 2e-5 (大きなリング A。K64 の仕様)。
- 放射の β の因子 (B112)、TouschekLifetime・Intrabeam・RADINT の古典半径 (陽子)。
- **DETR ≥ 1 の INS と MEASURE** (モードが入れ替わる結合リング) は未確認。
- 4 次元の CALC の RADTAPER は実装しない (K64 の仕様)。SOL 領域 + RADTAPER は未確認。

**集団効果**
- PSPAC: gs の x と y の 16% の非対称は再現しない (元の例の乱数のビームか SPC8 の状態の持ち越しの可能性)。粗いメッシュの定量値は gs と違う (K64 の形状関数 (CIC・双 3 次・「ほぼ最近接」) を入れていない。k64 モードに K64 の形状関数の PSPAC を入れる案)。gs は PSPAC の状態を C の static 変数に持つので、1 つのセッションで TrackParticles を繰り返すと 2 回目から PSPAC が働かない (SPC8)。
- WSPAC の Emittance[] が gs と 4% (x)・3% (y) 違う件は B127 で解決。
- 空間電荷の k64 モード (K64 の動作 = 空間電荷なし)・空間電荷の未確認項目 (ウェーク・放射との組み合わせ、SOL 領域内の四極・MULT など) は v2.0.0 のまま。

**スクリプト・FFS・特殊関数・パッケージ**
- **複素の特殊関数の k64 モードの値**: gs が誤る点でも真値を返す (B118〜B120)。k64 は「gs と同じ」の規則を満たしていない。再現には K64 の gamma.f (複素の不完全ガンマ `cgamma2`・`czconvb`・`cpolylogl`) の直訳が要る。
- 先頭でない成分の既定のキーの変数は作らない。`DAPWIDTH` の窓 (既定 51 のまま。trackd.f の追跡順の直訳が要る)。評価器の除算の逆数の積 (1 ulp。全試験の基準値に影響)。文脈の `Names`。init.n が大域に作る名前の全部 (51 個)。`Message` の `Tkinter::undef`。
- `45DEG`・`SIN[0.1]`・未定義の名前を pysad が受ける、`1+` を受けない、MAIN の `E`・`Pi` の代入。Tk 9.0 の `cget` の単位つきの値 (`4p`)。
- `Parallelize` の Windows の並列、numba と fork の併用。
- CSROYMatrix (B60) の gs の wrongtype を k64 で再現するか (今は再現しない)。`OYStability`・`OYScan` の gs との比較 (gs が完走しない)。tffamsetup (fam 点): 独立な基準を 2 つ作った (Courant-Snyder の不変量、8 極の薄いキックのヤコビアン) が、有限振幅の周期解の条件そのもの・振幅依存チューンシフトを目標にした FIT の収束・結合や6 極のあるラティスの fam 点は未確認 (追跡の NAFF か6 極の解析的な振幅依存チューンシフトが要る)。
- FITP>1 の最終解: pysad は独立な最小値と 5e-8 で一致、gs は最小値の手前で止まる (残差が相対 5e-5 大きい)。
- Micado: 水平・2 平面同時・`Bump`・`SVDmode->1` は未確認。
- numba のキャッシュの自動無効化は入れていない (`pysad.nbcache` を手で走らせる)。

**MEASURE・TRCOD** (v2.0.0 の (e) のまま): gs との比較は 1 件も無い。C の確認 (閉軌道の式、固有エミッタンス、ビームサイズ、FIXCOD・ALIGN) は独立な基準で誤りが見つからなかった。周期解を matched 分布として 1 周追う試験は、δ が大きいと色収差で崩れる (DP = 1e-3 で 40%)。

### (e) 試験と回帰の状態 (2026-11-18)

| 項目 | 結果 |
|---|---|
| 全試験 | `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本 [`test_bend_feeddown.py`](../pysad/tests/test_bend_feeddown.py)・[`test_coord.py`](../pysad/tests/test_coord.py)・[`test_sext_radiation.py`](../pysad/tests/test_sext_radiation.py) を `--ignore`、`PYSAD_SAD_PACKAGES` つき、互換モードは k64 が既定): **1345 通過・8 skip・失敗 0** (担当 E の報告、2026-11-18)。DOC が収集だけを確かめた件数は 1353 件 (E の実行より 6 件多い。**要確認: 最終の通過数は再実行して更新する**) |
| v2.1.0 の追加分 | `test_v2_{A1_next,A2_*,B1_*,B2_*,B3_*,C3_*,D1_emit,D2_*,E_*}.py` の 412 件 (A1 40、A2 38、B1 87、B2 30、B3 32、C 53、D1 22、D2 67、E 43)。ほかに独立な計算の補助 `_a2_indep.py`・`_c3_indep.py`・`_e_spinkick.py` |
| 回帰 | `run_regression.py` (`QUICK=1`): **50 ベンチ FAIL 0** (k64。E の実行で 159 s) |
| 既存の試験の変更 | 動作が変わったための更新 (期待値の追従): [`test_v2_A_spin.py`](../pysad/tests/test_v2_A_spin.py) の 4 か所 (側帯の `depol_rate_linear`、b = bzh、B の他の成分を相対 1e-10)、`test_v2_C_track.py::test_radtaper_matches_gs`・`test_v2_C2_optics.py::test_radtaper_in_calc6d`・[`test_v2_C2_spch.py`](../pysad/tests/test_v2_C2_spch.py) の 3 件 (RADTAPER の符号反転で v2 の Emittance[] が gs と 1e-6〜3e-3 ずれるため許容誤差を v2 だけ緩めた、k64 は同じ)、`test_sol_optics.py::test_sol_calc6d_is_refused` → `test_sol_calc6d_runs`、A1 の警告の試験の置き換え、B1 の `test_emittance_f1_unsupported_…` → `…_implemented_by_e`、B3 の光子の試験 (BEAMBEAM + BEND)、[`test_geo_fit_sig.py`](../pysad/tests/test_geo_fit_sig.py) (絶対許容 1e-25 を相対 1e-14 に)、`test_tclbridge.py::test_tk_version_and_font_system` (Tk 8.6 か 9.0)、[`test_v2_C2_script.py`](../pysad/tests/test_v2_C2_script.py) (`1/0`・`SQRT(-1)` を gs の値に) |
| skip | GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip。`PYSAD_SLOW=1` で遅い試験が走る (SB の収束の警告、POL の kick の応答、大きなリング A の MEASURE)。5 skip の内訳は未確認 |
| Python 3.14 | `py314` の venv で GUI 系 25 件と D2 の 90 件 (`test_v2_D2_*`・[`test_spin.py`](../pysad/tests/test_spin.py)・[`test_corr.py`](../pysad/tests/test_corr.py)・[`test_csr.py`](../pysad/tests/test_csr.py)・[`test_fam_points.py`](../pysad/tests/test_fam_points.py)) が通る。全試験は 3.14 では走らせていない |
| 限界 | v2 の値は gs では確かめていない (gs が K64 のバグを持つため)。SP8 (大きなリング B、8 GeV) は数値だけ (試験・例に使わない)。実機の強い結合・挿入光源のずれは未確認 |

### (f) v2.0.0 の (e) 「要判断・未解決」のうち v2.1.0 で解決したものの扱い

下の v2.0.0 の (e) は変えていない。次の項目は v2.1.0 で解決・更新した (残りは (d) に引き継いだものか、v2.0.0 のまま)。

| v2.0.0 の (e) の項目 | v2.1.0 での扱い |
|---|---|
| スピン追跡の 500 周後の `spin angle` (example/07) | **解決**: 原因は N8 (B83) の刻みの誤差。K64 の 0.8794 rad は物理の値ではない。EPS=1E-5 で gs と v2 が 1e-8 |
| tradke の b (B106)、tbgrote (B109)、RADINT の CKL1Y (B107)、ANGLE=0 の BEND の端の角 (B108) | **解決・修正**: 全部 v2 で修正。B107 は QUADRU も (B126)。B108 は K64 の誤りと確定 (端の項が約 2 倍) |
| K64 の spdepol・v2 の Polarization の限界 | **更新**: 側帯の項を入れた (B125)。追跡との差は 0.25 → 0.43。残りは (d) |
| CALC (tbende) と追跡 (tbend) の、ずれた BEND | **解決**: B114 (出口の符号)。K1≠0 は B117 |
| B21 の tcod.f (rfsw=.false. は δ を動かさない) | **解決**: D1 項目 1 で pysad の不具合を直した (空洞があって RFSW OFF) |
| RADINT の BEND の K0 (B74) | **解決**: ANGLE=0 の K0 は B115 で v2 が入れる。ANGLE≠0 の BEND の K0 は変えていない |
| DELCOR の規格化 (B99)・CAVI の符号 (B40) | **確定 (一部)**: 規格化の根拠は確定、意図は推定のまま。CAVI の符号は観測できない |
| RADTAPER の追跡の閉軌道 | **解決**: B111 (RADTAPER の追跡の重みも) |
| CALC (4 次元) の RADTAPER | **解決 (実装しない)**: K64 の 4D の CALC は RADTAPER を見ない (gs で確認) |
| Emittance[]・EMIT の POL の量 (SOL 領域内) | **解決**: B130 |
| N12: `Twiss["DETR"]` | **解決**: (c)。DETR ≥ 1 は未確認 |
| gs との差の未解決 (BEND の F1 + RAD、L≠0 の CAVI の DX・DY、EMIT のあとの USE) | **解決**: gs と全桁一致 |
| 未対応 (SOL・UND まわり): track6 の TRPT の加速 | **解決 (一部)**: SOL 領域の外の CAVI は実装 (POL + TRPT + 加速)。SOL 領域内・RF つき MULT・TCAVI は SADError のまま |
| 領域内の未確認 (領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3、加速 + RAD) | **確認**: B1 項目 6、E (bz = 20 まで、SpinTune の解析式 1e-9、放射損失 1e-3)。QUAD・MULT のスピンは 1e-4 台 |
| `TrackParticles` が MASS・CHARGE を渡さない | **解決**: (c) |
| 光学・Emittance への WSPAC | **解決**: (b)。gs との差の 4% は B127 |
| PSPAC の gs との 16% の差 | **原因を特定**: メッシュの分解能 (16% は再現しない)。警告を足した |
| NOKEEPEXP の `tfvcomp` | **解決**: 成分ごとの FREE |
| 複素引数の PolyLog・HypergeometricU・Zeta | **解決**: (b) |
| `Parallelize` の本当の並列 | **解決 (POSIX)**: Windows は逐次 |
| `Symbol["x"]`・`Message[f::usage]`・文脈の未確認・`Ffs;` | **解決**: (b) |
| B60 (CSROYMatrix)・B68 | B60 は原因を特定 (zj の head の上書き)。B68 は未解決のまま |
| MEASURE・TRCOD の gs との比較・MEASURE の初期分布 | v2.0.0 のまま (gs との比較は無い) |

## v2.0.0 (2026-11-17)

以下は v2.0.0 の一覧 (変更していない。v2.1.0 で解決したものの扱いは上の v2.1.0 の (f))。互換モードの使い方 (a) は v2.1.0 でも同じ。

### (a) 互換モードの使い方

| 項目 | 内容 |
|---|---|
| モード | `v2` (既定): K64 のバグを直した動作。`k64`: v1 と同じ、K64 (gs) の動作 |
| 切り替え | 環境変数 `PYSAD_COMPAT=k64` (または `v2`)。Python からは `from pysad import compat; compat.set_mode('k64')` (`None` で環境変数に戻す)。`compat.mode()` が現在のモード、`compat.v2()`・`compat.k64()` が真偽 |
| 試験・回帰 | 既存の試験は gs との数値比較なので k64 モードで走る ([`pysad/tests/conftest.py`](../pysad/tests/conftest.py) が `PYSAD_COMPAT` の既定を `k64` にする)。v2 の試験は `test_v2_*.py` で、フィクスチャ `v2_mode` で v2 を明示する。回帰 (`run_regression.py`) も k64 で走らせる |
| gs と同じ結果が要るとき | `PYSAD_COMPAT=k64` を付ける。`py -3.12 run_sad.py lattice.sad script.sad` は既定が v2 なので、gs と数値を比べるときは必ず付ける |
| 例 | [`example/v2_B/`](../example/v2_B/)・[`example/v2_C/`](../example/v2_C/) は既定 (v2) で走る。K64 と同じ動作は `PYSAD_COMPAT=k64` (v2 だけの機能は k64 では動かない・無視される) |

- **モードで切り替えない変更がある**: pysad 自身の不具合の修正 (§(d))、K64 に無い機能・動かない機能のうち gs と同じ結果が出るもの (LOSSMAP、RADTAPER など)、pysad が gs に合わせて足した機能 (NOKEEPEXP、INDEX など)。表の「k64 モード」の列でモードごとの動作を書く。
- **k64 モードは「v1 と同じ」ではない点が 2 つある** (gs に合わせて変えた): (1) `For` の本体を実行しない (B71。v1 の pysad は実行していた)。(2) 利用者のコードの `BeginPackage`・`Begin` が名前空間を作る (v1 は無視して全部大域。SAD の `Packages` の読み込みは変えていない)。
- 内部の約束: 追跡のカーネルが `rs[1]` (v2 の印) を読むので、`rs` は必ず `track6n.make_rs()` で作る。

### (b) v2 で結果が変わるもの (K64 = k64 モードの動作との違い)

表の「K64」は K64 (gs) の動作で、k64 モードの動作でもある (例外は注記)。ID は [`21_sad_bug_report.md`](21_sad_bug_report.md) の B 番号。
k64 モードを持たず「従来どおり未対応」のものは「K64」の欄にそう書いた。試験名は [`pysad/tests/`](../pysad/tests/) の下。

#### 放射・エミッタンス・光学 (記録 A)

試験は `test_v2_A_*.py` (67 件)。使った条件: 3 GeV の FODO リング、偏向磁石 L=1.2 m・15 度 (ρ=4.58 m)、RF 500 MHz。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B06 | RADINT の BEND の L/ANGLE | 全 BEND に最後の BEND の L/ANGLE | その BEND 自身の L/ANGLE | 曲率の違う 2 種の BEND のリングで自前の 3x3 行列の積分と I1〜I5 が 2e-5 で一致 (K64 は I1・I4・I5 が 1e-3 以上ずれる) | `test_v2_A_radint.py::test_per_bend_radius_v2` |
| B07 | RADINT の6 極 | 直前の DL を使う。L=0 は 0.2 m。先頭が SEXT だと NaN | 自分の L と K2。L=0 は寄与なし | 軌道 0 なら6 極は I1〜I5 を変えない (K2=3 と 0 で 12 桁同じ) | `…::test_sextupole_first_and_thin` |
| B12 | RADINT の定数 | 単精度 (1e27、1.E-6、0.2、2.3026) | 倍精度 | ELOSS = Cγ E⁴ I2/(2π) と 1e-13 で一致 (K64 は相対 1.15e-8 小さい) | `…::test_double_precision_constants` |
| B13 | RADINT の CAVI | 最後の空洞の値 | VOLT は総和。FREQ・HARM は電圧のある最初の空洞 (違う周波数が混ざると注意) | 2.5 MV + 1.5 MV → 4000 kV。νs = √(hαV cosφ/2πE) と 2e-3 | `…::test_cavity_voltage_is_summed` |
| B72 | RADINT の `QUANTUM LIFETIME(hours)` | exp(2.3026·τ[秒]) | τ[秒]/3600 | τq = τE/2·exp(ξ)/ξ、ξ = (δmax/σδ)²/2 | `…::test_quantum_lifetime_hours` |
| B73 | RADINT の `LINEAR CHROMATICITY` | −I9/4π (磁石の中だけの ∮γ ds で色収差にならない。例: −2.36、実際は −7.30) | 運動量を ±1e-4 ずらした光学のチューンの差分 dν/dδ (6 極・偏向の寄与を含む) | 同じスクリプトの `DP0` をずらした CALC の差分と 2e-5 で一致 | `…::test_chromaticity_is_tune_derivative` |
| B74 | RADINT の MULT・SOL | 黙って無視する | 場のある MULT・SOL があれば「入っていない」と表示 (積分には入れない)。k64 モードは表示なし | — | `…::test_mult_warning` |
| B75 | Jx の端の項 (docs/46 の「RADINT と EMIT の Jx の差」) | Emittance[]・EMIT・追跡は、端の角 e の蹴りを偏向角 θ に入れるが行路 al1 を変えない。I4 の端の項が 2 倍になる。矩形磁石のリングで Jx = 1.097 | 端の角のある面で場の中の行路 al1 を x·tan(e) だけ短くする (微分も D 行列に入れる) | (a) 解析式 Jx = 1 − I4/I2 と (b) 厳密な幾何 + Robinson の Je = 2 + dln l/dδ が 2e-6 で一致 (E1=E2=0: 0.92082、0.1309: 0.94147、0.5: 1.00575)。**RADINT が正しく K64 の EMIT が誤り**。v2 の Emittance[] は 1.00586 (EMITDIV=20 で 1.00580)、追跡の一周写像の固有値は 1.0033。gs の Emittance[] の Jx は 1.0974、同じ gs の RADINT は 1.00575 | [`test_v2_A_damping.py`](../pysad/tests/test_v2_A_damping.py) (7 件) |
| B76 | K1 つき BEND の追跡の曲率半径 | tbendi は放射で運動量が変わっても rhoe を入口の値のまま。減衰の和が 4 にならない (K1=0.1 で 4.05、Jx=0.7385) | スライスごとの放射のあとに rhoe = rhob·(1+δ) を更新 | 和が 4.000、Jx = 0.6875 (解析 0.68762、RADINT 0.68762、Emittance[] 0.6883) | `…::test_track_edge_and_sum_rule_v2` |
| B77 | MomentumCompaction (Emittance[]) | RFSW OFF は −R56/C だけ (Σ R5j η_j が抜け、FODO の例で 0.022898、正しくは 0.022137: 3.4% 大きい)。RFSW ON は μ² から (μ²/12 だけ大きい: 0.25%) | RFSW OFF は −(R56 + Σ R5j η_j)/C。RFSW ON は 4 sin²(μ/2) | 自前の I1/C と 2e-4 (FODO の例は 6 桁)。BucketHeight 0.009877 → 0.010045 (RADINT 0.010045、RFSW ON 0.010047)。docs/46 の「BUCKET HIGHT の差」はこの α の差 | `test_v2_A_radint.py::test_momentum_compaction_is_i1_over_c_v2` |
| B09 | tchge の入口のロール + ずれ | tbrote が回転行列を上書きし、一周行列に回転が入らない | 入口の回転を行列に掛ける | `_tchge_jac` が `_tchge_cod` の数値微分と 1e-8 で一致。ロール + ずれのある6 極・BEND のリングの CALC のチューンが追跡の一周写像の固有値と一致 (6 極 1e-9、BEND 4e-6。K64 は 0.02〜0.18 ずれる) | [`test_v2_A_tchge.py`](../pysad/tests/test_v2_A_tchge.py) (5 件) |
| B56 | tbrote のスピン回転 srot | 行の逐次更新で行列積にならない | (I+rr)·srot | 直交性が 1e-14 で保たれ、入口と出口が互いに逆 (K64 は 1e-5 以上崩れる) | [`test_v2_A_tchge.py`](../pysad/tests/test_v2_A_tchge.py) (3 件) |
| B11 | 4D の CALC の縦方向の Twiss 列 | 作り直さず 6D の値が残る | 既定値 (AZ=0、BZ=1、他 0) | 4D では縦方向を計算しない | `test_v2_A_ffs.py::test_calc4d_resets_longitudinal_columns_v2` |
| B15 | `UnstableLevel` | 常に 0 | 不安定な点の数 nstab | MatchingResidual = 10·nstab + NetResidual (例: 2 と 20) | `…::test_unstable_level_v2` |
| B16 | CalculateOptics の stabz・tracez | 未初期化 (k64 モードは 1 と 0) | idim=3 は tracez = 2cos(μz) と stabz、idim=2 は stabz=1・tracez=0 | Emittance[] のシンクロトロンチューンと 2e-3 | `…::test_calculate_optics_idim3_returns_z_trace_v2` |
| B39 | tradke (SOL 内) の dbeam(8) | B12 − B44 + b·B13 | B12 − B34 + b·B13 | 力学的運動量 (px + b y, py − b x) の共分散 A B Aᵀ | `test_v2_A_spin.py::test_kinetic_momentum_covariance_formula`、`…::test_tradke_px_py_diffusion_uses_b_y_py_v2` |
| B58 | BEND の端の縦磁場の積分 bsir0 (Emittance[]・追跡の POL) | 出口の項が常に捨てられる。K1 の項の次元が合わない。スピン軌道結合に端の項の y・x 依存が入らない。追跡の tbendi は y/ρ がない | 入口で +(y/ρ + k x y)、出口で −(…) (k = K1/L)。軌道についての微分もスピン軌道結合に入れる。追跡の tbendi も同じ | ∫Bs ds/(Bρ) = ψ(x,y)。軌道のずれた 1 要素のスピン回転が Emittance[] と追跡で一致: K1=0.3 の BEND 1.7e-6 (K64 2.0e-3)、K1=0 の BEND 7.5e-5 (K64 1.0e-3)。スピン軌道結合は追跡の差分と 2% (K64 は 1/6) | [`test_v2_A_spin.py`](../pysad/tests/test_v2_A_spin.py) (3 件) |
| B64 | MULT の bsir0 (Emittance[] の POL) | imag(cx)/aln (次元 1/m)。四極は入れない。分割は放射の大きさだけで決まる | ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し出口で引く。四極は k x y。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 になる数に分割 | MULT (K2=30) 3.1e-7 (K64 9.4e-5)、MULT (K1・K2・SK2) 3.8e-6 (3.1e-4)。四極の分割が 1 個だと歳差の微分が追跡の差分と 50% 違い、細かくすると 2% で合う | 同上 |
| B55 | `Polarization` (POL の平衡偏極度) | NaN (spintrm の exp(cb²−cφ)·cerfc(cb) が 0·∞)。k64 モードも NaN | (1) exp(cb²)·erfc(cb) を Faddeeva 関数 w(i·cb) で。(2) 平衡偏極度は線形のスピン軌道結合 (DK の (11/18)\|γ∂n/∂γ\|² の項) の脱偏極率 rate から P = (8/5√3)·sdamp/(sdamp+rate) | 平面リングで 92.376% (Sokolov–Ternov)。誤差のあるリングは FLUC つきのスピン追跡で測った脱偏極率と 0.5〜1.6 倍で合う (3 例: 3 GeV 1.6 倍、2.903 GeV 1.2 倍、4.2 GeV 0.53 倍。16 粒子、統計 ±25%)。K64 の式 (桁あふれだけ直す) は 3 桁以上大きい | [`test_v2_A_spin.py`](../pysad/tests/test_v2_A_spin.py) (5 件) |

#### 追跡・要素 (記録 A・C・C2・C_measure)

試験は [`test_v2_A_track.py`](../pysad/tests/test_v2_A_track.py)、`test_v2_C_*.py`、`test_v2_C2_*.py`。使った条件: 電子、0.5・1・3 GeV、周長 91 m の FODO リングと 4〜6 m の輸送路、ソレノイド 1〜3 T (bz = 0.3〜3 1/m)。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B08 | tsolqur (RAD、SOL 領域内の四極) | 最後の半ステップで最後の粒子の bzp を全粒子に使う | 粒子ごとの BZ/(1+δ) | 1 粒子の結果が、一緒に追う粒子に依らない (差 0。K64 は 2.6e-6) | `test_v2_A_track.py::test_solenoid_quad_rad_is_independent_of_other_particles_v2` |
| B14 | tsolqum (放射なしの MULT) の z の dv の項 | −0.25·(ap−dv)·aln0 | tsolqu と同じ −dv·aln0 | 同じ K1・L の QUAD と MULT の z が 1e-12 で一致 (K64 は 20 MeV・L=0.5 m で 1.9e-5 m ずれる) | `…::test_mult_z_equals_quad_and_drift_v2` |
| B10 | tdusrn の NaN の anp | npmax = 10000 回まわり乱数を 20001 個引く | NaN・±Inf なら光子を出さず乱数は 1 個 (有限の anp は K64 と同じ乱数列) | 同じ種で同じ戻り値と種 | `test_v2_C_track.py::test_b10_*` (3) |
| B17 | NaN と MAX/MIN | gfortran の MAX/MIN が NaN を無視する前提 | 乱数の消費だけ B10 で NaN を明示的に扱う。asinz・gfloor・sqrt1 は変えていない | 要判断 1 (§(e)) | (B10 の試験) |
| B33 | セッション最初の TrackParticles の最初の tapert | l_track = 0 で、LOSSAMPL の外の粒子が落ちず追跡が続く | 最初の tapert で落ちる (座標は 3·LOSSAMPL に頭打ち、落ちた周 1、落ちた要素は始点の 1 つ前) | 2 回目以降の K64 の動作と同じ | `test_b33_*` (2) |
| B35 | SOL 領域内に DRIFT BEND QUAD SOL MULT CAVI MAP APERT MARK MON 以外の要素 | SIGSEGV (k64 モードは NotImplementedError。L=0 は素通し) | `SolRegionError` (K64 の tsol.f が出すはずの文)。長さ 0 の SEXT なども黙って通さない | 仕様 (tsol.f の case default) | `test_b35_*` (1) |
| B37・B78 | UND の磁場 | 磁場のキーワードが無く、ドリフト。さらに undulator.f の陰的刻みの右辺で dpds を二重に掛け (B78)、磁場があっても収束力が刻み幅に比例して消える。k64 モードは従来どおり (FBX は NotImplementedError、光学は恒等) | `FBX`・`FBY` [T] を受ける (追跡・CALC・Emittance[]・EMIT。RAD で放射も)。磁場つきカーネルは dpds を 1 回だけ掛ける。**v2 では磁場なしでも UND は長さ L のドリフトとして光学に入る** | 平面アンジュレータの解析式: 軸上の軌道 (出口で x = x' = 0)、行路長 (1e-6)、垂直の収束 k = F0²/2 (5e-4)、シンプレクティック性 1e-7、放射損失 (2/3) r_e γ³ (F0²/2) L (1e-4)。リング: Δνy = ∫βk ds/4π (1%)、ΔU0 (2e-3)、減衰率の和 = 2U0/E (1e-3) | [`test_v2_C_und.py`](../pysad/tests/test_v2_C_und.py) (9) |
| B38・B79 | PHSROT | 行列がすべて 0 で粒子の座標が 0 になる。PX PY PZ は指定できない。拡散の乱数が非対角の U で共分散を壊す (B79)。k64 モードは追跡が NotImplementedError、光学は恒等 | 意図された写像 M = N·Rot(2πPX, 2πPY, 2πPZ)·(1−D)·N⁻¹ + 拡散 (RAD・FLUC)。拡散は N U (√λ ⊙ ξ) | M0 のシンプレクティック性 1e-13、固有値 = exp(±2πi P)、リング (ST PHSROT) の Twiss・チューン = 入力 (1e-9)、平衡エミッタンス b/(2d−d²) (Emittance[] 1e-6、追跡の統計 8%)、減衰率 = D | [`test_v2_C_phsrot.py`](../pysad/tests/test_v2_C_phsrot.py) (5) |
| B40・B62・B99 | DELCOR | MARK には DETR・DPX、CAVI には HARM・PHI に (符号も逆に) 書く (B40)。corlength≠0 は波数 −1 の成分 cv(n) が未設定で、負の波数の重みが 1 つずれる (B62)。corlength≠0 の規格化が複素 c の絶対値の 2 乗平均 = 1 で DX の rms が delx/√2 (B99) | MARK・CAVI も DX・DY に他の磁石と同じ符号で入れる。cv(n) を対称に置く (乱数の順は同じ)。`s = 1/sqrt(2Σw²+1)` で DX の rms = delx、DY の rms = dely | 相関長 ≫ 周長で全要素のずれが同じ (QF・RFC・ST が 5% 以内、5 km で全桁同じ)。相関関数 exp(−Δ²/(2 x0²)) (300 標本で 5%、虚部 2%)。場の rms (x0 = 15 m): v2 0.987・1.009、k64 0.708・0.703。DELCOR コマンド (62.4 m のリング、150 回): rms(DX)/delx = 1.021 (k64 0.713)。B99 は作者の意図の裏づけが無い (推定) | [`test_v2_C_delcor.py`](../pysad/tests/test_v2_C_delcor.py)、[`test_v2_C2_delcor.py`](../pysad/tests/test_v2_C2_delcor.py) (7) |
| B65 | RADLIGHT の軌道 | 横方向の座標が入らず、BEND で向きが変わらず、BEND の中を記録しない | 横の軸を幾何から取り、BEND は円弧で進め、中も 1/γ の角度あたり `RADLIGHTDIV` 点 (既定 16) 記録。曲げのある MULT は出口だけ、CAVI も記録 | 振幅 0 の粒子の点 = 幾何 `LINE["GX"/"GY"/"GZ"]` (1e-11)、時刻 = 行路長/(βc) (1e-9)。偏向磁石の放射場のピーク E = e β²/(4πε0 ρ R (1−β)²) と RadiationField が 0.25% (RADLIGHTDIV 16)、5e-5 (64) | [`test_v2_C_radlight.py`](../pysad/tests/test_v2_C_radlight.py) (3) |
| B69・B80 | SOL 領域内の BEND・MULT・CAVI・APERT | gs が領域内の BEND の偏向をほぼ無視する。k64 モードは NotImplementedError。放射 (tdrift) は 1 刻みの弦から損失を出すので、らせんの回転角 bz·L が大きいと小さく出る (相対 (bz·L)²/12。bz·L = 1 で 8%) (B80) | BEND = 直線の座標のままの一様な双極磁場 (tsol.f の `tdrift(al, bzs, phiy, phix)`)、MULT = tmulti1、CAVI = tmultiacc を bz つきで。領域内の BEND・MULT は bz·L/ndiv ≤ 0.02 に刻む。APERT で落ちる | ローレンツ力 dp/ds = (p × B)/pz の数値積分 (solve_ivp、rtol 1e-13): BEND・K0 だけの MULT 1e-14、K1 つき 1e-10、6 極・8 極・スキューは 1e-8。放射損失 (2/3) r_e γ³ (K0/L)² L が bz = 0・0.4・2 で 2e-4。bz → 0 の領域内の要素 = 領域の外 (リングの Twiss・Emittance[] が 2e-7)。CALC の軌道 = 追跡 (1e-12)。bz = 3 T のリングで減衰率の和 = 2U0/E (2e-3) | [`test_v2_C_sol.py`](../pysad/tests/test_v2_C_sol.py) (30) |
| B81 (N6) | スピン追跡 sprot の縦磁場 bsi | bz = bz0/pzm + bsi (p0 で規格化したまま z 軸の向きに足す)。δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)。傾いた軌道で bsi·ẑ の垂直成分が B⊥ に二重に入る。双極 + ソレノイドの bsi が T-BMT と 3.6e-3 rad 合わない | SOL 領域と UND では ∫B∥dl/(Bρ) = [(bz pz + Bx px + By py)·al/pz + ψ]/(1+δ) を運動量の向きに置く | Thomas-BMT の積分: v2 1e-7〜1e-5、K64 の式 6e-4〜4e-3 | [`test_v2_C2_spin.py`](../pysad/tests/test_v2_C2_spin.py) (17 件中) |
| B82 (N7) | tsolqur の縦磁場と放射の基準 (SOL 領域内の四極、RAD) | bsi を最初・最後のスライスは半分で上書き (1 スライス分足りない)。放射の基準 pxr0 が正準運動量のままで、最初のスライスに (bz/2)(y, −x) の見かけの偏向 | スライスごとに bz·aln/(1+δ)。放射の基準は力学的運動量 (高速追跡も) | 軸上の粒子のスピンの回転 (1+a) bz L: K64 は半分、v2 は 5e-15。軸外の放射損失: v2 の四極と MULT が 3e-4 で一致、K64 は 2.2 倍 | 同上 (`[quad_sol]`) |
| B83 (N8) | sprot の刻みの歳差の角 (すべての要素の追跡のスピン) | (1+aγ)·2 sin(θ/2)。刻みが粗いと θ²/24 小さい (0.26 rad の BEND を 3 刻みで 3e-4) | 厳密な θ | 平面リングの 1 周の回転角 = 2π aγ (K64 の式は 4e-3 rad ずれる) | `test_ring_spin_tune_with_solenoid[0.0]` |
| B85 (N11) | AUTOPHI | TRPT でないとき phis = phic で dphis = 0 になり、DZ·w も PHI も効かない。光学 (phii = phic) と合わない。全空洞が AUTOPHI だと実効の RF の傾きが 0。k64 モードは NotImplementedError | 「位相の基準 = 設計軌道の到着時刻」(TRPT のときの形 phis = 0)。追跡 (領域の内外)・CALC (TRPT)・Emittance[] に実装。傾きは軌道の位相での値を足す | 設計軌道が DZ = 3 mm 早く着く輸送路: ありは設計軌道の粒子の δ が 0 (CALC 1e-12、追跡 1e-9)。gs (TRPT) と CALC の DZ・DDP・追跡が全桁一致 | [`test_v2_C2_solrf.py`](../pysad/tests/test_v2_C2_solrf.py) (AUTOPHI の 2、Emittance[] の 1) |
| B86 (N13) | SOL の trots の z の DZ | 入口は z に DZ だけ、出口は DZ·(1 − dv) | 入口・出口とも (1 − dv)。k64 モードは SOL のずれ・回転が未対応 | 差は dv·DZ (3 GeV で 1e-13 m 以下) | — |
| B41・B87・B88 (MEASURE の初期分布) | MEASURE の初期分布 (ttinit) | 軌道の揺らぎを MARK の添字 26..29 (ZPX ZY ZPY DP) から読む (B41)。QUIET の重心合わせの基準が分布の作り方と合わない (B87)。TWAKE/LWAKE ON で全粒子 z = 0 (B88)。MEASURE は無出力 (B42)。k64 モードは MEASURE を無視 | JDX JDPX JDY JDPY (34..37) を読む。QUIET は分布と同じ写像で求めた重心に合わせる。z の分布を作る (`ttinit(b41=True, center_fix=False, wake=True)` で K64 の動作) | n=200000 の JITTER ON の平均が v2 の基準と 4σ 以内 (K64 の基準から 20σ 以上)。MK の例で px の差 2.5e-5 (統計誤差 3e-7)。B41 は MARK のキーワード表の添字。**ソースの読みのみ (gs は MEASURE が無出力)** | `test_v2_C_measure.py::test_b41_*`、`…::test_ttinit_quiet_centroid_*` |
| B42 (MEASURE・TRCOD) | MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` | MEASURE は無出力、TRCOD は SIGSEGV、MAIN の TRACK は落ちる。k64 モードは MEASURE・TRCOD を「未対応 (無視)」、TURNS=-1 は 1 周のリングの追跡 | K64 のソースの意図どおりに動く (§(c))。MEASURE の間だけ TRPT を ON にして戻す | 初期分布の 2 次モーメント (n=100000) が統計誤差 1.3% 以内、tbdecoup の非結合 emittance が 4×4 行列の固有値と 1e-10、MEASURE の出口の分散行列 M Σ0 Mᵀ と 1.9e-7、TRCOD の軌道が CALC (INS) と 5.4e-11。**gs との比較は無い** | [`test_v2_C_measure.py`](../pysad/tests/test_v2_C_measure.py) (19) |
| B100〜B105 | TRCOD の FIXCOD・MVCOD・STNCOD・ALIGN、ORBITJITTER、PUTMEASURE・GETMEASURE | ソースの誤り・旧い番号の名残 (docs/21 の B100〜B105)。K64 の FFS からは基準を読めない。k64 モードは語として認識しない / 未対応 | §(c) の機能として実装 (基準は REFERENCE と同じ参照光学に置き `Twiss["RDX"]`・`Twiss["DDX"]` で読む) | FIXCOD の差 Δx が解析式と相対 1e-5 (輸送路)・1e-4 (リング)、ALIGN が平行移動と 1e-13 (BEND を含むと 4e-10)、PUTMEASURE → GETMEASURE が 7 桁。**gs との比較は無い** | [`test_v2_C2_measure.py`](../pysad/tests/test_v2_C2_measure.py) (20) |

#### 集団効果 (記録 C・C2・C_spch・C2_spch)

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B34 | BEAMBEAM の傾き回転 | cost = sin(傾き角)、sint = 傾き角 (回転にならず座標が潰れる) | (cos t, sin t) の回転。強ビームに x-y 結合があるときだけ結果が変わる | 回転に対する共変性: t = 0.3 で 1e-9 (蹴りに対する相対) | [`test_v2_C_beambeam.py`](../pysad/tests/test_v2_C_beambeam.py) (2) |
| B84 (N10) | 空洞の wake の重み (tcav.f・tmulti.f の `fw=fw0*(wsn+ws(n-1)*.5d0)`) | 空洞が 2 つ以上に分割されると蹴りの合計が (1.5n − 0.5)/n 倍。k64 モードは再現 | 台形の重み fw0·(ws(n) + ws(n−1))/2 (合計 1)。領域の外の空洞・RF つき MULT も | 位置によらない表で、エネルギーの変化が分割数によらない: 20 MV (1 分割) と 200 MV (4 分割以上) で 5e-3。K64 は 1.25〜1.5 倍 | `test_v2_C2_solrf.py::test_wake_weight_sum_is_one_in_v2_and_one_and_half_in_k64` |
| B66 | 旧式の `WAKE` コマンド | 配布に無い WAKECOMMAND.n を読もうとして失敗 (tfwake.f は死んだコード) | `WAKE L 要素 s1 w1 … ENDWAKE` (横は `T`) を WakeFunction の表として登録 (CAVI だけ。行をまたいで書ける §(c)) | 同じ表を `WakeFunction[…] = …` で与えた追跡と 2e-9 (gs の参照値と一致済みの追跡) | [`test_v2_C_wakecmd.py`](../pysad/tests/test_v2_C_wakecmd.py) (3)、[`test_v2_C2_wakecmd.py`](../pysad/tests/test_v2_C2_wakecmd.py) (6) |
| B59 | CSRSetupOY の戻り | dpjn が zjn と同じ記憶域 (zj に運動量が入り、dpj は 0)。点が無いと i1 が未初期化 | zj は位置、dpj は運動量。点が無ければ失敗 (CSR::irreg) | 蛙跳びの関係 z_{k+1} − z_k = p_k r56 dt1 (1e-9) | [`test_v2_C_csr.py`](../pysad/tests/test_v2_C_csr.py) (2) |
| B89〜B92 (SPC1〜SPC4) | SPAC (DRIFT の `RADIUS>0`) の空間電荷 | 力が 1/γ 倍に弱い (B89)。px を二重に (1+δ) で割る (B90)。メッシュが整数の割り算で np<1000 では縦 10 枚、np=2000 でも横 8×8 (B92)。k64 モードは SPAC・PSPAC・WSPAC を追跡に入れない | バンチの静止系の 3 次元のポアソン方程式 (円筒の壁は像の電荷) で正しい係数。24〜64 ノードの 2 次 B スプライン | 一様な長い円柱 x'' = K x/a² (K = 2 r_c λ/(β²γ³))、KV 包絡方程式 (5 MeV・1e11 個で a(1 m) = 3.996 mm vs 追跡 2·rms 5e-3)、Bassetti-Erskine 5e-3、リングのチューンシフト (Laslett) 追跡/モデル 0.97〜1.00。gs は 5 MeV で rms 0.6226 mm (解析 K: 1.639、K/γ: 0.645、v2 1.623) | [`test_v2_C_spch.py`](../pysad/tests/test_v2_C_spch.py) (22)、[`test_v2_C2_spch.py`](../pysad/tests/test_v2_C2_spch.py) (35) |
| B91 (SPC3) | MULT の `RADIUS` の空間電荷 | spkick の呼び出しの引数が合っていない。最後の回だけ長さが半分 | ANGLE=0・VOLT=0・SOL 領域の外の MULT に入れる。四極を MULT で書くと全周に入る | 四極の中の KV ビームの包絡 2·rms 0.07%・0.2%、リングのチューンシフト 追跡/モデル 0.94〜0.97 (SPAC)・Laslett の式と 1% | [`test_v2_C2_spch.py`](../pysad/tests/test_v2_C2_spch.py) |
| B93〜B98 (SPC5〜SPC10) | PSPAC (SPCH 要素) と SPAC の小粒子数 | `PSPACNZ=1` (既定) で gs が SIGSEGV (B93)。λ(z) の配分の不連続 (B94)。周の終わりのキックの長さ (B95)。状態が static (B96)。tspac は死んだコード (B97)。SPAC の np が少ないと SIGSEGV・Bus error・NaN (B98) | `PSPACNZ<4` はバンチに合わせた 64 ノードを自動で張る。2 次 B スプライン。区間の始めから数え終わりに残りを入れる。呼び出しごとに作り直す。粒子数によらず走る (2・5・15・33・81 個) | PSPAC の 3 周後の rms が gs と 4 桁目まで (1.777/1.042 mm と 1.776/1.042 mm)、Laslett の式 (PSPAC: 追跡/モデル 0.98・1.004)。B98 は gs で np = 2〜300 を走らせて原因を特定 | 同上 |

#### 特殊関数・数値 (記録 B。pysad 自身の不具合も含む)

試験は [`test_v2_B_specfn.py`](../pysad/tests/test_v2_B_specfn.py) (13 件)。実数の倍精度。mpmath は無く、scipy.special と scipy に依らない独立な計算 (閉じた式、級数、AGM、積分表示、Wronskian) の 2 つを突き合わせた。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B43 | `Hypergeometric2F1`・`BesselJ` | `Hypergeometric2F1[.5,.5,1,.99]` = 2.348414、`BesselJ[.5,30.]` = −.143924709693 (打ち切り) | 2.352715816780、−.143929653370 (pysad は以前から真値。モードによらない) | 2F1(½,½;1;m) = 2K(m)/π を AGM で、J_{1/2}(x) = sqrt(2/πx) sin x。相対 1e-13 | `test_b43_hyp2f1_besselj_precision` |
| B44 | `ExpIntegralEi` | `Ei[1.5]` = .164−2.985 I、`Ei[-.1]` = −2.2265 | 3.301285449130、−1.822923958419 (実数。モードによらない) | 級数 γ + ln\|x\| + Σ x^k/(k k!)、負の側は −E1(x)、大きい x は漸近級数。1e-10〜2e-12 | `test_b44_b54_expintegralei` |
| B45 | `HermiteH[3,.7]`・`LegendreP`・`Binomial[-3,2]` | NaN | −5.656、+.274954541697 (`LegendreP[3,1,.4]`)、2.73 (`LegendreP[2,2,.3]`)、6 | H3 = 8x³−12x、P_3^1 = −(1−z²)^½ (15z²−3)/2、P_2^2 = 3(1−z²)、P_4^3 = −105 z (1−z²)^{3/2} | `test_b45_*` |
| B46 | `HypergeometricU`・`BesselY`・`BesselK`・`BesselJ[-1,…]`・`BesselJZero` | INF・NaN・無応答 | 2、−1.534238651350、1.806573612779、−.557936507910、2.404825557696 (モードによらない) | U(1,2,z) = 1/z、Y_0 の級数、K_{3/2} の閉じた式、Wronskian J_1 Y_0 − J_0 Y_1 = 2/(πx) | `test_b46_*` |
| B47 | `InverseErf` | `InverseErf[.3]` = .272462711811 (Newton の係数が sqrt(π/2)、相対 1e-8 の誤差) | .2724627147267543 (scipy の erfinv)。`InverseErfc` は小さい z で桁落ちしない erfcinv。**モードで変わる** (k64 は K64 の値) | erf(InverseErf(x)) = x が 4e-16 以内 | `test_b47_*`、`test_fitgaussian_uses_true_inverseerf_in_v2` |
| B53 | 大きい引数の `BesselJ/K/Y`・`FresnelS`・`AiryAi` | `BesselJ[0,33.]` 相対 3e-3、`FresnelS[6.1]` = 1.3e7、`AiryAi[15.5]` = −16 | 真値 (モードによらない) | J_n・I_0 の積分表示、K_ν の積分、Fresnel の Simpson、Ai = (1/π) sqrt(x/3) K_{1/3}(2/3 x^{3/2})、Ai Bi' − Ai' Bi = 1/π | `test_b53_*` |
| B54 | `PolyLog` | 複素数 (虚部 1e-15 の雑音つき) | 実数 (z ≤ 1)。z > 1 は切断線の下側 (Li_2(2) = π²/4 − iπ ln 2)。非整数次は \|z\| の上限なし (反転公式・積分表示。s ≤ 7.9、\|z\| ≤ 1e12 で相対 1e-11) | Li_2(½) = π²/12 − ln²2/2、Fermi–Dirac 積分、Li_{s−1} = z dLi_s/dz | `test_b54_polylog_*` |
| B18 | `PolynomialFit` (x の平均が 0) | `xl=1` のまま二項係数の和を取り係数が壊れる (`{{-1,1},{0,0},{1,1},{2,4.1},{-2,4}}` の 2 次で {1.0257, 2.0486, 1.0143}) | {−.0085714, .02, 1.0142857}。**モードで変わる** (k64 は K64 の値。平均が 0 でないデータは両モード同じ) | numpy.polyfit と 1e-13、3 次の多項式の係数の回復 | `test_b18_polynomialfit_mean_zero_v2` |
| B52 | `Fit[data,fun,{x1,x2},…]` | narg を出して最小にならない値 | 最小二乗解 (モードによらない) | numpy.linalg.lstsq と一致 (係数 1e-6、ChiSquare 1e-8) | `test_b52_*` |
| B61 | 複素の初期値の `FindRoot` | 壊れた値 (`{z->Restrict[(.25+I),…]}`)。k64 モードは TypeError | 複素の Newton 法 (§(c)) | z² = −1 → i、e^z = 2 → ln 2、z³ = 1、連立 zw = i・z + w = 2 の残差 1e-12 以下 | `test_b61_findroot_complex_v2` |

#### スクリプト・FFS (記録 B)

試験は [`test_v2_B_script.py`](../pysad/tests/test_v2_B_script.py) (25 件、うち 1 件は POSIX だけ)、[`test_v2_B_ffs.py`](../pysad/tests/test_v2_B_ffs.py) (14 件)。使った条件: 1 GeV・周長 6〜7 m の四極だけのリング、3 GeV・周長 105.6 m の FODO 24 セル。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B71 | `For[ini,cond,incr,body]` | 本体を実行せず incr だけ繰り返す。3 引数は未評価。**k64 モードは gs と同じにした** (v1 の pysad は本体を実行していた) | 本体を実行 (Break・Continue も効く)。3 引数は本体なしのループ | t = Σ i² = 30 | `test_b71_for_body_v2`、`test_b71_for_body_k64` |
| B70 | `MessageName[f,"usage"]` (関数形) | f を無視して `Undefined MessageName General::"usage"` | f の本文を返す (無ければ `General::tag` の本文)。モードによらない | `General::mytag` の継承は gs の `h::mytag` の実測と同じ | `test_b70_messagename` |
| B22 | `MON`・`MCOD`・`TCOD` | コマンドでなく、後ろの語が式として読まれて `Missing value for PM*` で文が止まる | 「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない。k64 モードも同じ) | K64 の tffsa.f にこれらのコマンドは無い | `test_b22_mon_mcod_tcod` |
| B24 | `LINE["LENGTH","$$$"]` | narg エラー | エラーのまま。メッセージに位置の取り方 (`LINE["POSITION",loc]`、`LINE["LENG",loc]`・`LINE["S",loc]`) を足した | K64 のソース (tftwiss.f:508) が `LENGTH` を引数 1 個に限っている (仕様) | `test_b24_line_length_with_location` |
| B26 | `ORG B3 1.5 -2.5 0.5 10 20 30` | `1.5-2.5` と読まれ値が 1 つ足りない | 空白で区切った `-2.5` は別の値。`3/2`、`2 *5` などの式は書ける | K64 で曖昧さのない書き方 (`1.5 0-2.5 …`) と同じ ORG の値 | `test_b26_org_negative_numbers` |
| B30 | MAIN の `OFF RFSW` | FFS に届かない (FFS は RFSW ON で始まる) | 届く (`ON RFSW` も)。FFS の中の `RFSW;`・`NORFSW;` が後から優先 | FFS の `NORFSW;` と全桁同じ追跡結果 | `test_b30_main_off_rfsw_reaches_ffs` |
| B67 | MOMENTUM の既定 | 30 GeV で、書き忘れても何も出ない | 既定は同じ 30 GeV。MOMENTUM を MAIN でも FFS でも設定しないまま運動量を使う計算 (追跡・EMIT・Emittance・RADINT・CALC6D・DA ほか) に入ると 1 度だけ警告 | 明示の `MOMENTUM=30e9` では警告しない | `test_b67_momentum_warning` |

### (c) 新しく使える機能と、モードによらない追加・修正

「モード」の列: **v2** = v2 モードだけ (k64 モードは従来どおり未対応・無視・恒等)、**両** = 両モードで同じ (pysad の追加・修正、または gs と同じ結果が出る機能)。
例は [`example/v2_B/`](../example/v2_B/)・[`example/v2_C/`](../example/v2_C/) (README に表がある)。使える範囲は「確かめた条件」で、それ以外は未確認。

#### 放射・光学・要素

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **UND の磁場** (`FBX`・`FBY` [T]。B37) | v2 | 追跡 (`TrackParticles`、MAIN の TRACK、DA)、CALC の Twiss、Emittance[]・EMIT・CALC6D・Touschek。RAD で放射 (FLUC も)。`SLICE` は `POLE` の 20 倍以上 (刻みが周期の整数倍だと磁場が見えない)。`KX`・`QY` は K64 と同じ。試験は 3 GeV、1 T・1.8 T、λ = 5 cm、20 周期。磁場なしでも UND は長さ L のドリフトとして光学に入る。`SLICE=0` は恒等 | `BX`・`BY` という名前 (K64 にも無い。ValueError)、端の磁極の整形。PHI≠0 の入口はベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま) | [`example/v2_C/C1_undulator.sad`](../example/v2_C/C1_undulator.sad)、`history/v2_changes_C.md` §2.1 |
| **PHSROT** (B38) | v2 | PX PY PZ はチューンの単位の位相、Dij は正規座標での減衰行列 (RAD で効く)、Bij は 1 回あたりの拡散行列 (RAD かつ FLUC)。R1..R4、EX EPX EY EPY ZX ZPX ZY ZPY も使える。追跡 (track6n)、CALC、Emittance[]・EMIT・CALC6D・Touschek。(z, δ) の回転は向きが逆 (matrix.f の m_prot_set_6d。縦のチューンは −PZ と出る)。分散つきの Twiss を CALC で読むと BX・EX は入力と 1e-3 ほど違う | k64 モード。拡散行列に負の固有値があると ValueError。スピン追跡の経路 (track6) では位相空間の写像だけで、PHSROT はスピンを動かさない (リング全体を 1 個で表すとスピンチューンは 0。偏極の計算には使えない) | `C2_phsrot.sad`、C §2.2、C2 §2.3 |
| **SOL 領域内の BEND・MULT・CAVI・APERT** (B69) | v2 | BEND: `L`、`ANGLE + K0`、`ROTATE + DROTATE`、`RAD`。直線の座標のままの一様な双極磁場。MULT: K0..Kn と SK0..SKn、`DX DY ROTATE`、端のフリンジ、RAD、`VOLT`。CAVI: `L`、`VOLT`、`PHI`、`DX DY ROTATE`、端のフリンジ。追跡、CALC (INS/CELL) の Twiss・軌道、Emittance[]・EMIT。条件: bz = 0.3〜3 1/m で確かめた | BEND の `K1`・`DX`・`DY` (NotImplementedError。四極成分は QUAD か MULT で)、ANGLE≠0 の MULT、`V1 V20 V11 V02`、SEXT・COORD・BEAMBEAM など (SolRegionError)、CALC6D・RADCOD つきの CALC (SADError)。未確認: 領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3 1/m を超える場 | `C3_sol_region.sad`、C §2.3 |
| **ずれて傾いたソレノイド** (`GEO=1` の SOL の `DX DY DPX DPY`、`DX DY DZ CHI1..3`、端の放射 `F1`) | v2 | 追跡 (高速追跡と track6)、CALC の Twiss と軌道、Emittance[]・EMIT、幾何 (`LINE["GX"…]`) と幾何の fit。領域の中の要素の軸はソレノイドの軸 (領域の中の Twiss・軌道はソレノイドの系)。`F1≠0` は RAD で端の蹴りを長さ F1 の一様な偏向として放射させる。領域内の BEND・MULT・CAVI・K1=0 の四極も軌道を追う | CALC6D・RADCOD つきの CALC、追跡を領域の途中から始める・終える、領域より上流に軸まわりの回転 (GCHI3≠0) があるラティスで追跡を先頭以外から始めた場合、領域内の APERT のスピン追跡。未確認: 軸が設計軌道と 10 mrad 以上ずれる場合、GEO のない領域に手で書いた CHI の符号の約束。gs との数値比較は行っていない | `C2a_sol_tilt.sad`、C2 §2.1 |
| **領域内の空洞の TRPT の加速・wake・AUTOPHI** | v2 | `TRPT; RFSW;` で領域内の CAVI・`VOLT≠0` の MULT が加速する (追跡、CALC (INS))。`LINE["GAMMABETA", k]` が各位置の基準運動量。wake は領域内の CAVI・RF つき MULT に付けられる (TRPT のとき)。`AUTOPHI=1` は空洞の位相を設計軌道が着く時刻 (Twiss の DZ) を基準に測る。RF つき MULT の Emittance[] は bz≠0 の領域の中だけ | スピン追跡の経路 (track6) での TRPT の加速、`V1 V20 V11 V02`、ANGLE≠0 の MULT。AUTOPHI の意味は K64 のソースの読み (SADHelp に説明が無い)。未確認: 領域内の空洞の wake の横の蹴りを gs と比べていない、加速 + RAD | `C2b_sol_rf.sad`、C2 §2.2 |
| **スピン追跡 (POL + RAD、track6) の SOL 領域・UND・PHSROT** | v2 | ビームは 9 列 (7 列目 = atan(sz/sx)、8 列目 = sy)。スピンは RAD が ON のときだけ動く。通る要素: SOL 領域 (ずれ・回転・F1、領域内の DRIFT・QUAD・BEND・MULT・CAVI・SOL)、UND (磁場つきも。縦磁場と周期の途中で終わる端のポテンシャルの飛びの歳差も入る)、PHSROT (位相空間の写像だけでスピンは動かさない)。精度: 一様磁場の刻みは厳密、四極・多極の中は刻みの粗さ (1e-6〜5e-5。`EPS` を小さくすると細かくなる)。リングのスピンチューン: 平面リング (0.5 GeV) + 長さ 1 m のソレノイド (0・1・2.5 T) で cos(μ/2) = cos(π aγ) cos(φ/2) と 2e-5 | **Emittance[]・EMIT の POL の量 (SpinTune・Polarization など) は SOL 領域の中が K64 の式 (N6) のまま**: ずれて傾いたソレノイドではスピンチューンで 2.3e-4、n0 の向きで 4e-3 rad 違う。領域の外の Emittance[] のスピンチューンは N8 の分だけ追跡と違う (平面リングで 2.8e-5)。領域内の APERT、k64 モード | `C2c_spin.sad`、C2 §2.3 |
| **Polarization (POL) の平衡偏極度** (B55) | v2 | 平面リング 92.4% (Sokolov–Ternov)。誤差のあるリング: 脱偏極率が FLUC つきの追跡と 0.5〜1.6 倍 (3 例)。`emit6d.emittance6d(…, pol=True)['spin']` に `equpol`・`depol_rate`・`taudep`・`spin_spread`・`spin_orbit`・`equpol_k64` (K64 の式の参考値) | 線形 (SLIM 相当) の見積もり。高次の側帯共鳴、運動論的偏極 (b·γ∂n/∂γ の項) なし。エネルギー幅の大きいリングでは側帯の寄与が増えるはず。確かめたのは小さな 1 つのリング。POL のときは四極を細かく分割するので Emittance[] が遅くなる | A §3・§7 |
| **RADTAPER** | 両 (gs と同じ結果) | `RFSW; RAD; RADCOD; RADTAPER;` が ON のとき、Emittance[]・EMIT・TrackParticles・CALC6D が磁石の強さを閉軌道の局所の運動量に合わせる。`Twiss["DDP"]` に閉軌道の δ。閉軌道の z は 0、位相の基準は `EquilibriumPosition` (trf0)。gs の例: 軌道 −7.887e-5 → −2.198e-7 (pysad −2.173e-7: gs は収束判定が緩く、`CODCONV = 1e-4` で −2.17319735e-7 と 6 桁一致)、U0 8 桁一致、εx 1.7372e-7、εz 1.9077e-5 | CALC (4 次元の光学、K64 の qtwiss の radtaper は未翻訳)、SOL 領域内の要素、RFSW OFF。見つからないとき `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける。追跡の閉軌道は Emittance[] と 1e-6 m 程度違う (K64 も同じ) | `C4_loss_and_taper.sad`、`C2d_taper_optics.sad`、C §2.5、C2 §2.4 |
| **LOSSMAP と損失** | 両 (gs と同じ結果) | `LOSSMAP;` で TrackParticles の戻りの flag の後に 2 行 (落ちた要素の l_track、落ちた周。生存は 0 と nend+1。POL のときは 10・11 行目)。POL + RAD の TrackParticles も LOSSAMPL・LOSSDZ で落ちる。gs の例と一致 (要素 {0,0,242,242,242}、周 {31,31,4,1,1}) | MAP 要素 (ExternalMap) が落とした粒子の要素番号 (0)、空間電荷が落とした粒子の要素番号。POL なしで高速追跡に乗らない要素の並びは SADError | `C4_loss_and_taper.sad`、C §2.4 |
| **MULT のスキュー成分** (`SK0 SK1 SK2 …`) | 両 (pysad の修正) | CALC・Emittance[]・追跡。gs と追跡 4 成分が 2e-13、Twiss (BX BY R1 R2 DX DY) が 1e-9 | — | §(d) N2、`C6_misc.sad` |
| **ResponseMatrix の TRPT** | 両 (gs と同じ結果) | 応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。gs と全桁一致 (1 → 1.2 GeV) | — | `C6_misc.sad` |
| **RADLIGHT の軌道** (B65) | v2 | `RADLIGHTDIV = n;` で BEND の中の点数 (1/γ の角度あたり。既定 16)。3 GeV・0.26 rad の BEND で 2.5 万点 | SOL 領域を含む区間 (SADError)、COORD による幾何の変化、曲げのある MULT の中の点。DX・DY・F1・FRINGE のある BEND は中の点が 1e-6 m 程度ずれる | `C5_radlight.sad`、C §2.6 |
| **DELCOR (v2 の規格化)** | v2 | corlength > 0 の DX の rms = delx、DY の rms = dely (K64 は 1/√2)。CAVI・MARK にも入る | rms は多数回の平均としての値 (GCUT を掛ければ小さい)。相関長が格子間隔より短いと線形補間で rms が下がる。corlength=0 の枝は別の量 (分散 = delx²·s(1−s/周長))。[`example/v2_C/C6_misc.pysad.txt`](../example/v2_C/C6_misc.pysad.txt) の DELCOR の行の値は √2 倍になる | `C2w_delcor_rms.sad`、C2_misc §2 |
| **旧式の WAKE コマンド** (B66) | v2 | `WAKE LONGITUDINAL 名前 s1 w1 … ENDWAKE` (横は `TRANSVERSE`)。ENDWAKE (または `;`) まで次の行に続く。CAVI だけ。名前はワイルドカード可。点数に上限なし (400 点で確認)。追跡は WakeFunction と同じ (TRPT と LWAKE/TWAKE) | 表の途中の `;`、CAVI 以外。ENDWAKE を書き忘れると次の `;` までを表として読んでエラー | `C2w_wake_table.sad`、C2_misc §1 |
| **RADINT の警告** (B74) | v2 | 場のある MULT・SOL があれば「入っていない」と表示 | — | A §1 |
| **EMITDIV・EMITDIVB・EMITDIVQ** | 両 (pysad の修正) | Emittance[] の放射の分割の倍率。MAIN の `EMITDIV = 20;` も効く。gs: `EMITDIV=20` で DampingRate が −.000286831320204 (pysad k64 と 12 桁一致) | — | §(d) |

#### MEASURE・TRCOD と周辺

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1`** (B41・B42) | v2 (K64 の gs は MEASURE が無出力、TRCOD が SIGSEGV) | 先頭が MARK のビームライン (輸送路、リング。リングは `CELL; CALC;` の後)。NP ≥ 2。MOMENTUM を設定する。輸送路は先頭 MARK に BX BY AX AY (と EX.. DX..) を書き EMITX EMITY DP SIGZ を MARK か FFS の変数で与える。**MARK に DP を書かないと DP は 0.01 (既定) で運動量の広がり 1% (一様) の分布**。追跡は TrackParticles と同じ経路 (RAD・FLUC・RFSW・APERT の損失)。`FFS["MEASURE 位置"]` は `{np, {重心 6}, {分散行列 6×6}}`。MEASURE の TRPT は追跡の間だけ ON (K64 は戻さない) | **gs との比較は無い**。`DRAW` の中の粒子追跡 (K64 でも死んだコード)。TRGAUSS OFF の分布と TWAKE/LWAKE ON の MEASURE は動くが独立な検証は薄い (LWAKE・TWAKE は一様・階段状のウェイクで重心と広がりを確認)。**強い結合のリング (モードが入れ替わる形) の MEASURE は未確認** (`Twiss["DETR"]` が常に 0: §(e))。MAIN の TRACK TURNS=-1 は粒子種の語 (PROTON など) を見ない | `C2m_*.sad`、C_measure、C2_measure |
| **TRCOD FIXCOD / MVCOD / STNCOD / ALIGN、ORBITJITTER、PUTMEASURE / GETMEASURE** | v2 | FIXCOD: EX EPX EY EPY DX DPX DY DPY を基準として保存し `Twiss["RDX",k]` (基準)・`Twiss["DDX",k]` (差) で読む (REFERENCE と同じ参照光学に置く。`REF;` は基準も上書きするので戻すには `TRCOD STNCOD;`)。ALIGN: 対象は BEND QUAD SEXT OCT DECA DODECA MULT の DX・DY だけ。ORBITJITTER: 今の Twiss の DX・DY に位置ごとの正規乱数。PUTMEASURE・GETMEASURE: 環境変数 SADPOUT・SADPIN のファイル (Windows は普通のファイル、Linux・WSL は FIFO を想定) | ALIGN の SOL・MONI。GETMEASURE の対象 (DX・DY) は判断。FIFO の待ち合わせは未確認。区間に BEND を含む ALIGN は平行移動と 4e-10 m まで | `C2m_trcod_reference.sad`、`C2m_putget_measure.sad`、C2_measure |

#### 集団効果

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **空間電荷 SPAC・PSPAC (+ SPCH 要素)** (B89〜B98) | v2 | `TrackParticles` の中。SPAC は `DRIFT` (`RADIUS>0`) と SOL 領域内のドリフト、MULT の `RADIUS` (ANGLE=0・VOLT=0・領域外。`SPACSTEP`)、PSPAC は `SPCH` 要素。電子 5 MeV〜50 MeV、1e8〜1e12 個/バンチ、ビーム半径 0.1〜2 mm、ダクト 3 mm〜3 cm。マクロ粒子は 1000 個以上 (チューンシフトを数 % で見るなら 4000 個以上)。条件: 粒子に z の広がり、`RFSW` ON (`NORFSW` はエラー。`RADIUS<0` の絞りだけなら走る)、ビームの外縁が 0.9 R 以内。FFS の変数 `PBUNCH`・`SPACSTEP`・`SPACNX/NY/NZ`・`PSPACNX/NY/DX/DY/NZ/DZ/NTURNCALC`・`WSPACBEAM`・`WSPACRADIUS`。フラグは FFS の最上位に書く (`Do[...]` の中では切り替わらない) | **k64 モードは空間電荷を追跡に入れない**。MULT の `ANGLE≠0`・`VOLT≠0`・SOL 領域内の MULT (メッセージを出して空間電荷なし)。四極・偏向・6 極の中 (偏向のあるリングでは SPAC は全周に入らない。全周は PSPAC + SPCH か WSPAC)。光学・Emittance への空間電荷。MAIN の TRACK、DA、スピン追跡。複数バンチ、壁の抵抗、円でないダクト。未確認: 空間電荷とウェーク・放射の組み合わせ、バンチの静止系で相対論的になる場合 (γ·x' ≳ 1)、電子以外 (式は確かめたが SAD スクリプトからは使えない: `TrackParticles` が MASS・CHARGE を追跡に渡さない) | `C2s_mult_radius.sad`、`example/05_tracking/flag_SPAC` (作り直した)、C_spch、C2_spch |
| **WSPAC・SELFCOD** (強-弱のモデル) | v2 | 各要素の入口でガウス分布の力を蹴る。ビーム行列は既定で Emittance の平衡値、`WSPACBEAM = {εx, εy, σz[, σδ]}` で Twiss から作る。`WSPACRADIUS` で縦の力の電位の基準。電子 5 MeV (WSPACBEAM)・0.5 GeV (Emittance)、Δν = −3e-5〜−0.025。gs と 1% (0.5 GeV のリング、2e10 個) | 光学 (CALC)・Emittance への空間電荷 (ビーム行列は空間電荷なしの平衡値)。垂直のエミッタンスが 0 のリングは `WSPACBEAM` を与える。SOL 領域内・SELFCOD と閉軌道のあるリング・放射との組み合わせは値の確認なし | `C2s_wspac.sad`、C2_spch §2 |

#### スクリプト・FFS・特殊関数 (モードによらない追加・修正)

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **文脈** (`BeginPackage`・`Begin`・`End`・`EndPackage`) | 両 (k64 も名前空間を作る) | K64 (init.n と gs の実測) の規則: 名前は「現在の文脈 → `$ContextPath`」の順に探し、無ければ現在の文脈に作る。`$Context`・`$ContextPath`・`Context[sym]`。大域に公開する関数は `` `Name[…] := … ``。組み込み関数・FFS の関数はそのまま使える。`Names["Stat`*"]` は v2 が文脈の記号の一覧、k64 が gs と同じ `{"Stat`"}`。SAD 自身の Packages は従来どおり文脈を無視して大域で読む。K64 との差 (意図的): pysad の System` は K64 より広い、表示は文脈名を `$ContextPath` にあるものだけ省く、`Symbol["Stat`x"]` は記号を返す | `$Context = …`・`$ContextPath = …` への代入、`WithContext`。入れ子の文脈は確かめていない | `B2_package_context.sad`、B §2.4 |
| **Fork・Shared・Wait・Kill・Sleep** | 両 (Windows の Fork は v2 でエラー、k64 は従来どおり pid を返す) | `OpenShared`・`Shared[s]` (数・文字列・記号・リスト・式・規則)・`SharedSize` (K64 と同じ式。`Shared::toolarge`)。POSIX (Linux・WSL・macOS) の `Fork[]` は `os.fork` で分ける。`Wait[]`・`Wait4`・`Kill[pid]` (この pysad が作った子だけ)・`Sleep[t]`。WSL で example と 3 個の子の試験 | Windows の Fork (`System::error` で打ち切り。並列は `Parallelize`、pysad は逐次)。ラティスを使う重い計算を子で走らせる使い方 (numba のスレッドと fork の併用) は未確認。子の `Print` は出力先が StringIO のとき親に届かない | `B3_shared_system.sad`、B §2.5 |
| **System** | v2 (許可なし) | 既定では実行しない (`$Failed` を返し案内を出す)。許可: 環境変数 `PYSAD_ALLOW_SYSTEM=1` (起動前)、または `pysad.miscfn.allow_system(True)`。SADScript からは変えられない。許可があるとシェルで実行し終了コード (成功 0) を返す。`[ -e path ]`・`[ -d path ]`・`[ -f path ]`・`[ -r path ]`・`[ -s path ]`・`test …` は許可なしでも判定する。`PYSAD_SYSTEM_SHELL` で sh・bash を指定 | k64 モードは許可なしで 0 (「実行しません」)、許可ありで戻り値 1。SAD の Packages が呼ぶ個々のコマンド (`import`、`latex` など) は未確認 | `B3_shared_system.sad`、B §2.6 |
| **複素引数** (Bessel J/Y/I/K (実数の次数)、Hypergeometric0F1・1F1・2F1 (実数のパラメータ)、Airy) | 両 | 積分表示・級数・Wronskian と 1e-12〜1e-13 | 複素の PolyLog・HypergeometricU・Zeta (scipy に実装が無い・分岐の向きを確かめる参照が無い)。`PolyLog` は非整数次で負の s が大きい (s < −10 程度) 精度は未確認 | `B1_special_functions.sad`、B §2.2 |
| **複素の FindRoot** (B61) | v2 | 式が変数について正則と仮定した Newton 法 (前進差分のヤコビアン + 最小二乗、減速つき)。連立可。`MaxIterations`・`AccuracyGoal` が効く | 変数の範囲 `{var, ini, {min, max}}` は複素では無視。Abs・Conjugate・Re・Im を含む式は収束を保証しない。実数の初期値で複素の値になる式は従来どおり `General::wrongval` | B §2.3 |
| **MOMENTUM の警告** (B67)、**`MON`・`MCOD`・`TCOD` の警告** (B22) | v2 | MOMENTUM 未設定で運動量を使う計算に入ると 1 度だけ警告 (CALC と幾何では出さない。SADScript の式の中の `MOMENTUM` の直接の使用では出ない) | — | `B4_ffs_v2.sad` |
| **MAIN の代入と `FFS` の文の読み方** | 両 (pysad の修正) | gs の動作 (小さな台本 15 本) に合わせた: MAIN の `NAME = 式;` (数 [単位]、`+ - * /`、`^` (右結合)、関数 `SIN COS TAN ATAN SQRT EXP LOG LN COSH`、`PI`、`MOMENTUM`) が要素定義と FFS に届く。`FFS USE RING;`、行の途中の `FFS`、`FFS;NPARA=3;` の読み。32 個の式を gs の値と 2e-14 | gs で構文エラーになる形は代入しない (`Sqrt[16]`、`ABS(…)`、`2*-A` など。`(1+2)*3` だけは読む)。残した違い: `Ffs;` (大小まじり)・`ffs use ring;` は pysad は読む、`FFS` の次の行が USE 以外のとき、gs が FFS の終了後に MAIN の文として読む残りは pysad は読まない、`A = 3 B = 4;` は代入しない、`LINE[{"EMITY",…},1]` は MARK の値 (0) を返す、要素定義の中の MAIN の関数 (`L=2*SIN(A)`) は未確認 | `C2w_main_script.sad`、C2_misc §3・§4 |
| **NOKEEPEXP・RENUMBER の INDEX・`LINE["INDEX",…]`** | 両 (gs と同じ結果) | `NOKEEPEXP` のとき EXPAND/NOEXPAND を付けない CALC・GO は `calexp = .not. tfvcomp()`。INDEX 属性が 2 以上の要素の番号はその値から。gs の実測と一致 | `tfvcomp` は常に偽 (成分ごとの FREE が無い)。RENUMBER 後の `LINE["POSITION",…]` の gs の癖は写していない | `B4_ffs_v2.sad`、B §1 |

### (d) pysad v1 の不具合の修正

K64 のバグではなく、pysad 自身の誤り・未対応だったもの (K64 の番号 B72 以降には入れていない: [`21_sad_bug_report.md`](21_sad_bug_report.md) §6)。記載は「モードによらない」ものを含む。

| 項目 | v1 の不具合 | 修正 | 根拠・試験 |
|---|---|---|---|
| `LegendreP[n,m,z]` の奇数 m の符号 | scipy の `lpmv` に (−1)^m を掛けていた (`lpmv` は Condon–Shortley の位相つき。`LegendreP[3,1,.4]` が −.2749…) | 掛けない。**+.2749545416973507** | 式の m = 1 ± 1e-6 の値と閉じた式 (`test_b45_legendre_integer_m_is_limit_of_formula`)。gs は NaN (B45) で参照が無く気づけなかった。モードによらない |
| MULT の `SK0 SK1 SK2 …` (N2) | `lattice._mult_params` が `KsL` を読まず、スキュー成分を黙って無視 (CALC・Emittance[]・追跡) | `akr(n) = cmplx(K_n, SK_n)·cr1^(n+1)`、cr1 は複素 K1 の主軸 (tpara) | gs: 追跡 4 成分が 2e-13、Twiss が 1e-9 (`test_v2_C_sol.py::test_mult_skew_components_match_gs`)。モードによらない |
| `CalculateOptics` の TRPT の加速 (N3) | CALC は入れるが CalculateOptics は入れない | TRPT かつ RFSW の加速を区間に入れる (区間の途中からでも) | gs: 終端の BX AX BY AY EX EPX NX NY が全桁、CalculateOptics − CALC = 0 (`test_v2_C2_optics.py::test_calculate_optics_includes_trpt_acceleration`) |
| 加速区間の分散 (N9) | qmat2twiss の gr = gammab(出口)/gammab(入口) の向きが逆 (入口に分散のある加速区間の EX EPX が gs と合わなかった: 0.1109 / 0.0734) | K64 のとおり gammab(出口)/gammab(入口) を入口の分散に掛ける | gs と全桁 (薄い空洞・厚い空洞)。`test_dispersion_through_accelerating_cavity_matches_gs` |
| EMITDIV・EMITDIVB・EMITDIVQ | pysad は無視していた | Emittance[] に効く (MAIN の `EMITDIV = 20;` も) | `EMITDIV=20` で DampingRate が gs と 12 桁 (`test_emitdiv_is_used`、`…_in_main_section`) |
| SOL の L と S 座標 (B63) | pysad は S 座標 (LENG) にだけ SOL の L を入れていた (46.9) | 入れない (gs の LENG 46.6)。SOL に L のキーワードは無い (SADHelp) | `test_v2_A_ffs.py::test_sol_has_no_length`。モードによらない |
| 特殊関数 | 非整数次の `PolyLog` が z < −12・z > 60 で精度が落ちる (z = −140、s = −1.5 で相対誤差 1。\|z\| ≥ 150 は NaN)。複素引数の Bessel・超幾何・Airy が TypeError | `_polylog_far` (反転公式・積分表示)。複素引数の対応 (§(c)) | s ≤ 7.9、\|z\| ≤ 1e12 で相対 1e-11 |
| `SharedSize`・`LINE["INDEX",…]`・INDEX 属性の番号・`h::tag` の `General::tag` への継承・POSIX の `Fork`・`Sleep` | `SharedSize` が 1 (gs は 32)、`LINE["INDEX",…]` が 0、INDEX 属性が効かない、`Fork` が分けない、`Sleep` 未定義 | gs の動作に合わせた (§(c)) | gs の実測 |
| `fn_track` (B33 の修正) | POL + RAD のとき後ろの `_compat.v2()` が UnboundLocalError (関数内の import が原因。k64 の [`test_spin_track.py`](../pysad/tests/test_spin_track.py) が 5 件落ちた) | モジュール先頭の import だけにした (`measure.py`・`_exec_word` も同様) | 試験 |
| `Element["RADIUS", ..]` (DRIFT・MULT) | 0 を返していた | 定義の値を返す (gs と同じ。追跡・光学の値は変わらない) | `test_mult_radius_readable`、`test_element_radius_readable`。モードによらない |
| 描画の試験 ([`test_plot.py`](../pysad/tests/test_plot.py)) | Tk がマウスの下の項目に付ける `current` タグで判定が変わる | 読み取りで `current` を外す (判定は同じ) | — |
| 例 `example/05_tracking/flag_SPAC` | `NOSPAC;`・`SPAC;` を `Do[...]` の中に書き、フラグが切り替わらず「gs でも SPAC の効果が出ない」と誤記 | 最上位に書き、解析式と比べる例に作り直した (gs は力が 1/γ: B89) | — |

### (e) 要判断・未解決

**v2.1.0 で解決したものは上の v2.1.0 の (f) に対応を書いた (この節の文面は v2.0.0 のときのまま)。**

各記録の「要判断」「残した違い」「使えない範囲」を、解決済みのものを除いて集めた。「決定済み」は親の判断で変更していないもの。

#### 放射・エミッタンス・スピン (記録 A)

- **スピン追跡の 500 周後の `spin angle` (`example/07_spin_solenoid_special/flag_POL_RADPOL`)**: v2 の既定では 1.1516 rad で、gs と k64 モードの 0.8794 rad と違う。原因は特定していない (A の B58・B56・K1 つき BEND の追跡の修正のいずれかによる意図した違いの可能性が高いが、独立な根拠で 500 周後の値を確かめていない)。小さなリングの 1 周あたりのスピンの回転角を解析式 (スピンチューン aγ と、ソレノイド・K1 つき BEND の端の項) で確かめる試験が要る。

| 項目 | 見つけたこと | 直さなかった理由・状態 |
|---|---|---|
| tradke の dbeam の b = bzh·0.5 (B106) | 力学的運動量は px + bzh·y なので b = bzh のはず | bzh の定義の別の約束がある可能性を排除できない。SOL 内の放射のビーム行列の「ビームの広がり」の項だけに効く |
| tbgrote の trans(5,6) (B109) | コメントアウトされた微分の式ではなく別の式。dr=0 でも dpzi·(pr+1)/pzi が残る | 作者が意図して置き換えた式。2 次の小さい量 |
| RADINT の CKL1Y (B107) | sin の誤記に見える (結合のあるときの I4V) | 結合のあるリングの独立な I4V の計算を持っていない |
| RADINT の BEND の K0 (B74) | 無視される (ANGLE だけ) | 仕様かもしれない |
| RADINT の空洞の位相 | VOLT を単純に足す (PHI を見ない) | 位相の違う空洞のベクトル和にするかは仕様の判断 |
| K64 の spdepol (2 次の摂動、シンクロトロン側帯) | 大きさが追跡と 3〜5 桁合わない (側帯の項に 2π が余分、1/減衰率の因子)。v2 は線形の式に置き換えた | K64 の式を直すには導出が要る |
| v2 の Polarization の限界 | 軌道の振幅について 1 次 (SLIM 相当)。追跡との一致は 3 例で 0.5〜1.6 倍 (16 粒子、統計 ±25%)。高次の側帯共鳴、運動論的偏極なし | 確かめたのは小さな 1 つのリングだけ |
| CALC (tbende) と追跡 (tbend) の、ずれた BEND | 矩形磁石を DX=2 mm ずらすと CALC のチューンが 2.9e-5 変わる (追跡は変わらない。追跡が正しい) | tbende の 2 次の幾何の項。場所を特定していない |
| K1=0 の BEND の追跡のスピン (tbend + tradk) | 1 個の BEND のスピン回転角が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う | 放射の分割が粗い K64 の精度の範囲 |
| B21 の K64 のソース | tcod.f の rfsw=.false. の枝は δ を動かさない。z を閉じる δ は temit の呼び出し側で決まるはず | 呼び出し側を追っていない。pysad の動作は gs の実測と合い、物理としても妥当 |
| ANGLE=0 の BEND・K0 つき MULT の端の角 (B108)、tthine の bsir0 (B110) | B75 と同じ作りのはず / 次元が合わない | 独立な値が無い / 影響なし (直していない: ANGLE≠0 の BEND だけ B75 を直した) |
| 小さな Jx の残り | 既定の分割 (EMITDIV=1) では Emittance[] の Jx が 1.5e-3 大きい (K64 の精度)。Jx+Jy+Je = 3.9990 (U0/2E の 2 次) | `EMITDIV=20;` で 1e-4 まで合う |
| 範囲の限界 | 大きなリング・挿入光源・強い結合では確かめていない。RADINT の V 面 (I4V) は確かめていない。SOL 領域内の MULT・BEND のスピンは (Emittance[]) 未対応 | — |

#### 追跡・要素 (記録 C・C2)

| 項目 | 状態 |
|---|---|
| B17 (NaN の返し方) | 決定済み: 現状のまま。asinz・gfloor・sqrt1 の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効く。「NaN になった粒子をその場で落とす」に変えると落ちる要素・周が K64 と変わる (LOSSMAP の値)。マッチングの fuzz (tffsmatch) は未確認 |
| DELCOR の CAVI の符号 (B40) | 決定済み: v2 の +delx (K64 は設計の DX − delx)。SOL の式を写した意図があるなら逆 |
| DELCOR の corlength≠0 の規格化 (B99) | v2 は rms = delx に直した (作者の文書による裏づけは無い。推定)。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける |
| UND の光学 | 決定済み: v2 は磁場なしでも光学にドリフトとして入れる (追跡と揃える)。K64 の「光学では恒等」を仕様とみるなら磁場つきのときだけにする選択もある |
| B33 の落ちた要素 | v2 は「始点の 1 つ前の要素 (周回で数える)」。K64 の 2 回目以降は「前の追跡の最後の要素」 |
| スピン追跡の損失の k64 の最初の呼び出し | K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs 要素 60、pysad 15)。高速追跡 (track6n) は全桁一致 |
| 未対応のまま (SOL・UND まわり) | スピン追跡の経路 (track6) での TRPT の加速 (RF は入るが基準運動量は変えない)、領域内の APERT のスピン追跡、UND の端の磁極の整形、CALC (4 次元) の RADTAPER (K64 も無し)、bz = 0 の領域と領域の外の RF つき MULT の Emittance[] (C2 §2.2。C2 §1 の項目 2d は「bz = 0 の領域の中の CAVI・MULT は領域の外と同じ写像」とも書いており、**記録内の表現に食い違いがある: 要確認**) |
| Emittance[]・EMIT の POL の量 (SOL 領域内) | K64 の式 (N6 の bsi) のまま。`spin.tradke_spin` の `xpzb = xpz + (bsir0 + bzh·2·al)·pr²` を直すには 6 次元の微分 dxpzb も一緒に直す必要がある (未着手) |
| N12: `Twiss["DETR"]` が CALC の後で常に 0 | **v2.0.0 の時点では未修正の不具合 (v2.1.0 で修正済み、上の N12 の行)** (`ffscmd._rows_from_tws` が detr を渡していなかった。k64・v2 とも)。結合の強いリングの MEASURE に効く |
| N13 (trots の入口の DZ)、AUTOPHI の意味 (N11) | K64 のソースの読み。gs で確かめられたのは TRPT の AUTOPHI だけ。リング (TRPT でない) の追跡の AUTOPHI は K64 と違う |
| RADTAPER の追跡の閉軌道 | Emittance[] の閉軌道から追跡を始めても、1 周で x が 4e-7 m、z が 2e-6 m 動く (テーパーなし)。K64 の追跡と tturne の放射の入れ方の違い (pysad の k64 の追跡は gs と 1e-12 で合う)。直していない |
| 領域内の未確認 | 領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3 1/m を超える場、領域内の空洞の wake の横の蹴り (gs と未比較)、加速 + RAD。gs は SOL + 四極 + RAD で保存値と食い違う (B68) |
| gs との差の未解決 (docs/46) | BEND の `F1 FRINGE=3` + RAD の追跡、L≠0 の CAVI + `DX`・`DY`、EMIT のあとに `USE` したときの追跡 は v2 でも未解決 |

#### MEASURE・TRCOD (記録 C_measure・C2_measure)

| 項目 | 状態 |
|---|---|
| gs との比較 | 1 件も無い (gs は MEASURE が無出力、TRCOD が SIGSEGV)。実装は K64 のソースの読み。乱数を引く順はソースどおりだが粒子ごとの値を gs と比べたことはない |
| MEASURE の TRPT の副作用 | 決定済み: v2 は追跡の間だけ ON にして戻す (K64 は戻さない)。K64 のとおりにしたければ MEASURE の後に `TRPT;` と書く (または `cmd_measure` の `finally` の 1 行を外す) |
| QUIET の重心合わせ (B87)・TWAKE/LWAKE の z = 0 (B88) | v2 は直した (K64 の式の誤りと判断。名残と判断)。TRCOD の途中の統計 (K64 は装置番号 0) は出さない |
| GETMEASURE の対象 (B103) | DX・DY の条件と判断。旧い DDX.. (基準との差) だった可能性があり、その場合の目標値は基準の分だけ違う |
| FIXCOD の基準の置き場 (B100) | REFERENCE と同じ参照光学に置いた (`measure._ref_rows` の 1 か所で替えられる)。ALIGN の対象に DECA DODECA MULT を足した (`measure._ALIGN_V2` を `_ALIGN_K64` にすれば K64 の BEND QUAD SEXT OCT だけ) |
| 前任の試験の古くなった 1 行 | C2 の記録では `test_v2_C_measure.py::test_measure_particle_loss_and_errors` の `'?TRCOD FIXCOD: not supported'` の 1 行 (FIXCOD の実装で出なくなる) と [`test_v2_C_delcor.py`](../pysad/tests/test_v2_C_delcor.py) の ⟨\|c\|²⟩ の期待値が古くなっていた。現在の試験は `'  trcodf'` と 2.0 を見ていて、全試験は通る (解決済み) |
| 未翻訳 | `DRAW` の中の粒子追跡 (K64 でも死んだコード)、PUTMEASURE・GETMEASURE を FIFO で使う待ち合わせ (未確認)、ORBITJITTER のパイプ |

#### 集団効果 (記録 C・C_spch・C2_spch)

| 項目 | 状態 |
|---|---|
| 空間電荷の k64 モード | 決定済み: K64 の動作 = 空間電荷なし。PSPAC だけ v2 の実装を k64 でも使うかは未決定 (gs に 3〜4 桁で近づくが「k64 = 空間電荷なし」が変わる) |
| `TrackParticles` が `MASS`・`CHARGE` を追跡に渡さない | 電子以外は追えない (空間電荷と関係のない pysad の既存の制限。他の追跡にも影響するので親の判断) |
| 光学・Emittance への WSPAC (K64 の `twspace` の反復) | 未対応。入れるとビーム行列が自己無撞着になり Twiss にも空間電荷のチューンシフトが出る (`emit6d`・`calc6d` の変更が要る) |
| 空間電荷の未確認 | ウェーク・放射との組み合わせ、SOL 領域内の四極・MULT、`ANGLE≠0` の MULT、静止系で相対論的な運動、放射 + 空間電荷の MULT (走るが PBUNCH=0 の一致だけ)。PSPAC の gs との 16% の差 (直線の例で x と y が違う) は原因を調べていない |
| 文書の更新 | `docs/46`・`docs/01` などの「SPAC は読むだけ」「gs でも効果が確認できない」の記述 (別の担当が更新) |

#### 特殊関数・スクリプト・FFS (記録 B・C2_misc)

| 項目 | 判断 | 理由・調べたこと |
|---|---|---|
| B24 `LINE["LENGTH",loc]` | 値は返さず案内だけ足した | K64 のソースが `LENGTH` の引数を 1 個に限る。位置を付けたときの意味 (位置か、軌道長か) がソースからもヘルプからも決まらない |
| B22 MON・MCOD・TCOD の機能 | 実装せず警告して読み飛ばす | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料が無い。軌道補正は CorrectOrbit パッケージ (docs/34) |
| `Out[n]:=` の gs の表示 | 写さない | gs は未定義の記号 1 つの文で前の `Out[n]` を繰り返す。K64 の表示の癖 |
| NOKEEPEXP の `tfvcomp` | 常に偽 | pysad の `FREE` はファミリー単位で成分ごとの変数が無い。成分ごとの FREE を実装するときに `FFS._tfvcomp` を直す |
| RENUMBER 後の名前の引き方 | pysad は `LINE["NAME"]` の表示と一致する位置を返す | gs は RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う (内部の表の更新の癖。INDEX の有無によらず以前からの差) |
| 複素引数の PolyLog・HypergeometricU・Zeta | 未対応 (SADError・TypeError) | scipy に複素の実装が無い、分岐の向きを確かめる参照が無い |
| `Parallelize` の本当の並列 | 逐次のまま | POSIX では Fork で実装できるが、numba のスレッドと fork の併用を確かめていない |
| `Symbol["x"]`・`Message[f::usage]` の gs との差 | そのまま | gs は `Symbol` が値を返し、`Message[f::usage]` が `???f::usage:  本文` を表示する。pysad は記号を返し、何も表示しない |
| 文脈の未確認 | `$Context = …` への代入、`WithContext`、入れ子の文脈 | §(c) |
| `Ffs;`・`ffs use ring;` など | pysad が gs より寛容 | §(c) の「残した違い」 |

#### K64 側の未解決 (docs/21)

B68 (gs の 10-19 ビルドが SOL + 四極 + RAD で保存した参照と食い違う: 現在の gs が領域内の四極の力を入れていない疑い)、B60 (CSROYMatrix が gs で wrongtype)、B49・B51 (RadiationField・RenumberElement が gs で SIGSEGV) などは v2 でも比べられない。B72〜B110 のうち確度が「疑い」のもの (B88、B99、B106〜B109) は K64 の作者の意図を確かめていない。

### (f) 試験と回帰の状態 (2026-11-17 20:30 に実行)

| 項目 | 結果 |
|---|---|
| 全試験 | `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本 [`test_bend_feeddown.py`](../pysad/tests/test_bend_feeddown.py)・[`test_coord.py`](../pysad/tests/test_coord.py)・[`test_sext_radiation.py`](../pysad/tests/test_sext_radiation.py) を `--ignore`、環境変数 `PYSAD_SAD_PACKAGES` つき、互換モードは k64 が既定): **937 通過・4 skip・失敗 0** (941 件収集) |
| v2 の試験 | `test_v2_A_*.py` 67 件、`test_v2_B_*.py` 52 件、`test_v2_C_*.py` 114 件、`test_v2_C2_*.py` 130 件 (計 363 件、29 ファイル)。v2 の動作は `compat.set_mode('v2')` で明示して試験する |
| 回帰 | `run_regression.py` (`QUICK=1`): **50 ベンチ FAIL 0**。k64 モードで走る (gs との比較) |
| skip | GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip (環境変数を付けると走る)。`PYSAD_SLOW=1` で遅い試験が走る。4 skip の内訳は未確認 |
| 条件 | Windows (Python 3.12)。`Fork` の試験は WSL (Ubuntu、Python 3.14、numba は代用の空モジュール) で確認 |
| 限界 | 試験の参照は解析式・独立な数値計算・gs。v2 の値は gs では確かめていない (gs が K64 のバグを持つため)。MEASURE・TRCOD は gs との比較が無い。実機の大きなリング・挿入光源・強い結合では未確認 |

## 記録と例

| 記録 (`history/`) | 内容 |
|---|---|
| [`v2_changes_A.md`](history/v2_changes_A.md) | 放射・エミッタンス・スピン・RADINT・tchge (B06〜B16、B39、B55・B56・B58・B64 ほか、Jx の 2.4% の差の結論、平衡偏極度) |
| [`v2_changes_B.md`](history/v2_changes_B.md) | 特殊関数・数値、SADScript、FFS の小項目 (文脈、Fork・Shared、System、複素の FindRoot) |
| [`v2_changes_C.md`](history/v2_changes_C.md)、[`C_measure`](history/v2_changes_C_measure.md)、[`C_spch`](history/v2_changes_C_spch.md) | 追跡・集団効果・要素 (UND、PHSROT、SOL 領域内、LOSSMAP、RADTAPER、RADLIGHT、旧 WAKE)、MEASURE・TRCOD、空間電荷 |
| [`v2next_A1.md`](history/v2next_A1.md)・[`A2`](history/v2next_A2.md)・[`B1`](history/v2next_B1.md)・[`B2`](history/v2next_B2.md)・[`B3`](history/v2next_B3.md)・[`C`](history/v2next_C.md)・[`D1`](history/v2next_D1.md)・[`D2`](history/v2next_D2.md)・[`E`](history/v2next_E.md) | v2.1.0。A1: 放射 (反作用の符号、古典半径、F1 の追跡)・CAVI・RADTAPER・PSPAC、A2: スピン角・Polarization の検証・RADINT・ずれた BEND、B1: Emittance[] の F1・K1≠0 の BEND の端・WSPAC・SOL の CALC6D、B2: FREE・複素の特殊関数・Parallelize・文脈・Packages、B3: 偏極の側帯・光子・幾何の fit・NAFF、C: MEASURE・RadiationField・SynchroBeta・DELCOR・B106〜B109・FLUC の乱数、D1: NORFSW・WSPAC の 8 倍・POL + 加速・光子・spdepol、D2: 要素定義の式・Tk 9.0・Micado・CSROYMatrix・nbcache、E: 傾いたソレノイドの POL・ANGLE=0 の BEND の F1・B108・脱偏極率の応答 |
| [`v23_R5.md`](history/v23_R5.md) | v2.3.0。R5: 4 次の共鳴駆動項 (6 極の 2 次)・追跡の 3 次のテイラー展開・数値の正規形・周波数マップの予想 (K64 の誤りなし) |
| [`v23_S1.md`](history/v23_S1.md) | v2.3.0。S1: 周波数マップの検証 (RDT・ADTS からの共鳴の予想、δ・誤差・k64・MPW・NAFF・格子、KEK-LS の gs との比較。K64 の新しい誤りなし) |
| [`v23_T1.md`](history/v23_T1.md) | v2.3.0。T1: elegant との比較 (LMA・Touschek 寿命。SAD → elegant の変換、光学の確認、ASTRID2・kek-ls_v03_12.lat の LMA の比較、一様な δ と局所の LMA の寿命の比較、式の差の切り分け)。B160 (K1≠0 の BEND のオフモメンタムの光学。U1 で K1 の暗黙の定義に再分類。直していない) |

| [`v23_U1.md`](history/v23_U1.md) | v2.3.0。U1: B160 の切り分け (K64 のソース、SADHelp、導出、独立なローレンツ力の積分。K1 は (1 + hx)∂b/∂x 一様の定義、B'' = −hB')。B161 (疑い) |
| [`v23_R6.md`](history/v23_R6.md)・[`S2`](history/v23_S2.md)・[`S3`](history/v23_S3.md)・[`L2`](history/v23_L2.md)・[`T2`](history/v23_T2.md)・[`U2`](history/v23_U2.md) | v2.3.0。R6: K0≠0 (MPW)・δ≠0 の ADTS (kick の強さの誤りの修正)・端の閉形式・非対角の正規形の補正・精密版 `fourth='r6'` (K64 の誤りなし)。S2: 周波数マップの信頼性 (誤差の統計・MPW・RFSW ON・`D_norm`・警告・gs と誤差つき)。S3: 共鳴線の予想の作り直し。L2: LMA の信頼性 (`check_every` の既定 1、RF バケットとアパーチャの区別、`verify`、周期性)。T2: Touschek の式 (B162、Piwinski 型が v2 の既定)。U2: B160 のオプション `set_bend_k1_def('lin')` と B161 の修正 |
| [`v23_Q1.md`](history/v23_Q1.md)・[`Q2`](history/v23_Q2.md) | v2.3.0。Q1: 周波数マップ・LMA・Touschek への受け渡し・追跡の設定の共通部分 (設計の判断、損失の判定の間隔、RF バケット、gs との比較)。Q2: 6 極の共鳴駆動項・ADTS・クロマティシティ (厚いレンズの解析式、追跡との照合、gs の ADTS の不一致)。まとめは上の v2.3.0 の節、解説は [`48_frequency_map_lma.md`](48_frequency_map_lma.md)・[`49_resonance_terms_adts.md`](49_resonance_terms_adts.md) |
| [`v2next_P1.md`](history/v2next_P1.md)・[`P2`](history/v2next_P2.md)・[`P3`](history/v2next_P3.md)・[`P4`](history/v2next_P4.md)・[`P5`](history/v2next_P5.md) | v2.2.0。P1: 偏極 (n が傾く・非平面のリングの緩和率、クロマティックな側帯、反転した BEND、追跡のスピンの ROTATE)、P2: 放射・エミッタンス (SynchroBeta の始点の分散、EMITDIV、tradke の符号と微分、RADTAPER、K1≠0 の BEND の端 `face`)、P3: TrackParticles の F1 と MULT の AE・B108 の残り、P4: SOL 領域内のスピンと光子・k64 の PSPAC・WSPACSC・DAPWIDTH の窓・要素ごとの FREE、P5: gs との検証 (MEASURE・DETR ≥ 1・複素の特殊関数・Micado・FITP・陽子の β・文法・Tk 9.0)。K64 の新しい誤りは B135〜B159 |
| [`v2_changes_C2.md`](history/v2_changes_C2.md)、[`C2_measure`](history/v2_changes_C2_measure.md)、[`C2_misc`](history/v2_changes_C2_misc.md)、[`C2_spch`](history/v2_changes_C2_spch.md) | 残りの仕上げ (ずれたソレノイド、領域内の空洞、スピン追跡、RADTAPER の収束、TRCOD の残り、旧 WAKE の表・DELCOR・台本の読み方、MULT の RADIUS・WSPAC) |

| 例 | 内容 |
|---|---|
| [`example/v2_B/`](../example/v2_B/README.md) | `B1_special_functions`・`B2_package_context`・`B3_shared_system`・`B4_ffs_v2` (特殊関数、文脈、共有メモリ・System・Fork、FFS の小項目) |
| [`example/v2_C/`](../example/v2_C/README.md) | `C1_undulator`、`C2_phsrot`、`C3_sol_region`、`C4_loss_and_taper`、`C5_radlight`、`C6_misc`、`C2a_sol_tilt`、`C2b_sol_rf`、`C2c_spin`、`C2d_taper_optics`、`C2s_mult_radius`、`C2s_wspac`、`C2m_trcod_reference`、`C2m_putget_measure`、`C2w_wake_table`、`C2w_delcor_rms`、`C2w_main_script` |
| [`example/v2_3/`](../example/v2_3/README.md) | v2.3.0 の例 (Q1: 周波数マップ・LMA・gs との比較・台本の関数、Q2: 共鳴駆動項・ADTS・厚い ⇔ 薄い・追跡と gs の比較。Q1 は `Q1_README.md`、Q2 は `Q2_README.md`) |
| [`example/v2_next/`](../example/v2_next/) | v2.1.0・v2.2.0 の例 (v2.2.0 は `P1_*`〜`P5_*`。P2・P3・P4・P5 は README (`P2_README.md` ほか)、P1 は台本と出力。`P1_response_samples_*.json` は決定論的な測定のデータ)。v2.1.0 の例。A1・B1・B3・C3・D2 は README (`A1_README.md`・`B1_README.md`・`B3_README.md`・`C3_README.md`・`D2_README.md`)、そのほかは台本 (`A2_*`・`B2_*`・`D1_*`・`E_*`) と出力 (`.gs.txt` = gs、`.pysad.txt` = v2、`.pysad_k64.txt` = k64、`.out.txt` = Python の例) |

## v2.4.0 (放射・エミッタンス)

記録 [`history/v24_RAD.md`](history/v24_RAD.md)、試験 `test_k64_tcod_notfound.py`・`test_v24_rad_pole.py` (独立な計算 `_rad_indep.py`)。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| tcod が収束しないときの Emittance[] (B163) | gs と同じに、収束していない軌道から続行して `???-Emittance[]-closed orbit not found.` (Stable = 0)。DDBA_15 の εx が gs と全桁一致 (これまで +0.33% 違った) | 従来どおり 6D 固定点を Newton 法で解く (メッセージなし) | gs の出力、tcod のトレース。収束した値は設計メモの値と 0.06% で合う |
| RADINT の ANGLE=0・K0≠0 の極の εx・Jx (B115 の残り) | 極を積分に入れない (K64 のまま) | 極の中を直線の基準軌道の場として積分 (扇形の h² の収束・I4 の h³ の項なし)。CALC の分散が極の分散の発生を含まない分を周期解で補正 | I1〜I5 が独立な積分と 5e-5 (ASTRID2 + MPW)。RF なしの Emittance[] と εx 1e-4・Jx 2e-4 |
| **EMITDIV の既定** | 1 (K64) | **4** (`compat.emitdiv_default()`、環境変数 `PYSAD_EMITDIV_DEFAULT`。台本が EMITDIV を与えればその値。EMITDIVB・EMITDIVQ は 1) | EMITDIV = 1 の εx の誤差: 大きなリング A +2.2%、PF +2.2e-3、矩形磁石の 16 セルの小リング +1.4e-3 (RF の位置が分散のある所でも同じ)。4 で 2.2e-4 以下。実行時間の増加 0.1 s 以下 |

**EMITDIV の既定の例外**: POL ON のスピンの計算では既定 16 (v2 の線形の脱偏極率は偏向磁石の刻みで収束し、EMITDIV = 1 は収束値 (約 1.87e-8) より +24%、4 は +20%、16 は +4%: 4.2 GeV・12 MV の FODO。`test_v24_pol_emitdiv.py`)。RADTAPER ON の閉軌道は k64 でも収束した軌道 (tcod の直訳は RADTAPER OFF)。試験の参照値が EMITDIV = 1 のもの (`test_v2_D1_emit.py` の WSPAC のチューンシフト、`test_v2_P2_rad.py` の SynchroBeta の Emittance[]、`test_pf_history_lattices.py` の v2 の εx) は更新または明示した。

Emittance[] と RADINT の差の整理 (ASTRID2 + MPW・16 セルの小リング + 3 極の挿入光源): RF なしの Emittance[] は RADINT と一致する。RF あり Emittance[] の Jx が +0.7〜2%、εx が −1.4% 違うのは、RF の位置の分散が効く物理的な効果で、独立な 4 次元の放射つき写像 (RK4 + 空洞) の固有値が再現する (Jx 0.7038 対 0.7040)。矩形磁石 (E1 = E2 = 0.5·ANGLE) で RF が分散のある場所にある 16 セルのリングの +1.7e-3 の差 (v2.2.0 の未決) は **v2 の誤りではない**: EMITDIV = 20 では独立な 4 次元の写像のビーム行列と 3e-5 で一致 (扇形も矩形も)、EMITDIV = 1 で +1.4e-3。


## v2.4.0 (パーサ・評価器・FIT・追跡の周辺: gs との照合)

記録は docs/46 §2・§4 の該当行、試験 `test_v24_A_parser_fit.py`・`test_v24_A2_mult_angle_twiss.py`・`test_v24_A3_das_damping.py`・`test_apert.py`・`test_map.py`。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| 実数の `a/b` | `a*(1/b)` (gs: `Log[100]/Log[10]` = 1.9999999999999998) | IEEE の `a/b` | gs と全桁 |
| FitFunction の行 | 不等式型 (gs の `tfsolv` が範囲外の添字を読む結果) | 等式の条件 | KEKLS_match_repeat が gs と 7〜8 桁 |
| NaN の FitFunction | 条件にしない (gs と反復の履歴まで同じ) | 同じ | gs |
| USE のたびの DP・EMITX 等 | MARK の正の値、DP は 0 以下なら 0.01 (`tffsinitparam`) | 同じ | gs |
| `Emittance[]` の tcod の初期値 | qcod が見つからなければ 0 から (gs の qcod は初期値を返す) | 同じ | pm1 が gs と 1e-9 |
| ユーザー関数 `At`、`A = 3 B = 4;`、`Off[f::tag]`、位置のリスト (`LINE["POSITION",{{1,2},{3,4}}]`)、`{f,g}[1,2]` のエラー | gs と同じ | 同じ | gs |
| `ElementValues` (Couple.n) とフィット変数の連動 | gs と同じ (結合を更新しながら数値微分) | 同じ | gs 1e-9 |
| TCAVI の `RANKICK`・`RANPHASE` (追跡) | 実装 (1 周ごとに tgauss 2 回) | 同じ | gs の 8 粒子・6 周と 1e-14 |
| APERT 再入射 (JDPX/JDPY)・PSPAC の NPARA > 1 | `NPARA=n;` で gs の子プロセスの分け方を再現 (粒子 12 個以上の乱数列の差の原因) | 同じ | gs と全桁 |
| LOSSMAP の要素番号 (MAP・空間電荷が落とした粒子) | MAP は gs と同じ (temap の `kptbl(j,4) = l`)。空間電荷は v2 のみ (gs と比べられない) | 同じ | gs (MAP) |
| `DynamicApertureSurvey` の `Damping->True`・`DampEnough` | 実装 (tpdamp) | 同じ | gs (6 極なしは全点一致。カオス的な縁は数周〜数百周違う) |

実装しないと決めたもの (理由は docs/46): MULT の SK0 (gs の応答が物理でない)、SOL の BOUND=0・BZ≠0 と領域内の SEXT・COORD (gs が SIGSEGV)、`Hold` の部分評価、UND の端の磁極 (基準が無い)、DAS の NPARA > 1 の RAD + FLUC (子プロセスの競合)。

追記 (v2.4.0 担当 A の第 3・4 弾): 単独の `NOEXPAND;`・`EXPAND;` は gs と同じく黙って受ける (以前は `?Undefined command or element` を出した。結果は変わらない)。FFS の `COUPLE QD QF -0.8; FREE QF; FIT; GO;` は gs と 1e-8 以内 (試験あり)。`DynamicApertureSurvey` は gs の NPARA = 4 と 1 で 510 点の ntloss が全点一致するので、pysad は NPARA を無視する (NORAD の地図)。KEK-LS の例 (`example/kek_ls/design/*`・`display/KEKLS_scan_target`) の `.pysad.txt` と 3 つの `.md` (`match_repeat`・`tune_survey_cell`・`scan_target`) を再生成後の値に更新した (`match_repeat` は gs と 2.7e-5 以内)。k64 の振る舞いの変更はすべて上の表のとおり (除算 `a*(1/b)`、FitFunction の行、USE での DP、tcod の初期値、NPARA の分け方、TCAVI の乱数、LOSSMAP の MAP の要素番号)。

## v2.4.0 (E: 組込の数学関数・特殊関数・リスト関数の監査、B180〜B199)

方針 (ユーザー指示 2026-11-22): v2 は K64 の特殊関数・数学の誤りを、誤りと指摘した上で v2 で直す。k64 (gs の再現) は、gs が誤る点・NaN を返す点は変えない (gs と 2026-11-22 に照合し、**gs が真値を返すのに pysad の k64 が NaN・実部だけだった点は pysad の誤りとして k64 も直した**。B185〜B196 の表と docs/21 の注記)。基準は mpmath (collective312 の venv、`mp.dps = 40`) の主値で、v2 の切断線は Mathematica と同じ。参照値は [`pysad/tests/specfn_audit_ref.json`](../pysad/tests/specfn_audit_ref.json)、試験は [`test_v24_E_audit.py`](../pysad/tests/test_v24_E_audit.py) (49 件)。記録は [`history/v24_E.md`](history/v24_E.md)。

**走査の範囲**: 約 4000 点 (実数: 負・小 (1e-300)・大 (1e300)・境界・整数次・非整数次、複素: 20 点の格子 (|z| = 1e-3〜35、負の実軸の上下)) × 約 110 関数。v2 の値が mpmath と相対 1e-9 で一致しないものを洗い出し、全部直した (残りは 12 桁の表示の限界・零点近傍・オーバーフロー・`Pochhammer[-.5,-.5]` = ComplexInfinity のみ)。B43〜B47・B53・B54・B118〜B122・B158 の v2 の値は再確認して全点が真値 (B47: `InverseErf`・`InverseErfc` は 1e-300〜1−1e-15 で一致、B121・B122・B18・LegendreP の奇数 m の符号・`For` の本体は試験で固定)。

| B | 内容 | v2 | k64 |
|---|---|---|---|
| B180〜B182 | `ArcCosh[x]` (|x| < 1 の符号)・`ArcTanh[x]` (|x| > 1 が NaN)・`ArcSin`/`ArcCos` (x > 1 の符号) | Mathematica の主値 | gs の値 |
| B183 | pysad 自身の誤り (gs は正しい): 奇数乗のオーバーフローの符号・`Floor[Infinity]`・`Mod[Infinity,2]`・複素数の `Round`/`Floor`・`Sign`・`Exp` のオーバーフロー | 修正 | **修正 (gs に合わせた)** |
| B184 | `Round` の半整数 (K64 は 0 から遠い方、Mathematica は偶数へ。**台本の結果が変わる**)、`Mod` の桁落ち | 偶数へ・正確な剰余 | 従来 |
| B185〜B189 | x < 0 の Bessel・負/複素の Shi・Chi・li・E_n (z < 0)・2F1 (z > 1) が NaN・INF・実部だけ。**gs (10-19 ビルド) は BesselJ/I/K・Chi・li・E_n・2F1 を正しく返す (pysad の k64 の誤り)** | 複素数の真値 | **gs と同じ真値に直した** (gs が NaN の `BesselY[x<0]` だけ k64 は NaN) |
| B190〜B196 | (gs が真値: 複素の 1F1・不完全ガンマ・ProductLog・Gamma[a,z]・LerchPhi・PolyLog[−n,z]・Binomial[−5,−2] → k64 も直した。gs が誤り・NaN: HermiteH の複素・LegendreP[ν,m,x]・U の大きい z・Li2 の小さい z・HarmonicNumber の小さい x・Binomial[3.5,−1] → v2 のみ) 複素の 1F1 (|z| ≳ 12 で scipy が相対 1e-4〜1)・U の精度・HermiteH の複素・LegendreP[ν,m,x]・ProductLog (z < −1/e・分岐点)・不完全ガンマ (a ≤ 0、z < 0)・LerchPhi (|z| ≥ 0.9)・PolyLog の小さい z と負の整数次・HarmonicNumber の小さい x・`Binomial[n,−k]` | 真値 | 上のとおり |
| B197 | `Insert`・`StringInsert` の負の位置 (gs は 1 つずれる) | Mathematica (末尾) | gs (pysad の StringInsert の不一致は直した) |
| B198 | `Table`・`Do`・`Sum`・`Product` の `{i,{a,b,c}}`、`Take[l,All]`・`Drop[l,All]` | 受ける | gs と同じエラー |
| B199 | pysad 自身の誤り (gs は正しい): `Inner`・`MapIndexed` のレベル・決まらない条件の `If`/`And`/`Or`/`Equal`・虚部 0 の複素数・`StringTake[s,0]` | 修正 | **修正 (gs に合わせた)**。`Inner` の行列は gs の「最上位の添字だけ対」を再現 |

**結果が変わる既定 (v2)**: `Round[x]` の半整数、`Insert[..,-1]`・`StringInsert[..,-1]`、`Equal` で決まらない比較 (従来 0 → 式のまま。`x==x` は 1)、`If[x,..]` の決まらない条件 (従来 Null → 式のまま)。k64 の変更は B183・B199 の pysad 自身の誤り (gs と合わせたので gs との一致は増えた)。

**複素の次数の `PolyLog`・`ExpIntegralE`** (v2.4.0 E で実装。k64・v2 とも。gs は対応している): `PolyLog[s, z]` は |z| ≤ 0.5 が級数、|Log z| ≤ 5 が μ = Log z の展開、それ以外 (|z| ≥ 2) は反転公式 Li_s(z) = −e^{iπs} Li_s(1/z) + (2π)^s e^{iπs/2}/Γ(s) ζ(1−s, ½ + ln(−z)/(2πi)) (`specfn_cx2.polylog_cs`。z > 1 の実軸は ln(−z) = ln z + iπ)。`ExpIntegralE[ν, z]` は z^{ν−1} Γ(1−ν, z): |z| ≤ 8 は級数、|z| > 8 かつ |arg z| < 2.6 は連分数、負の実軸の近く (|arg z| ≥ 2.6) は |z| ≤ 60 が Kummer 型の同符号級数 γ(a,z) = z^a Σ(−z)^k/(k!(a+k))、|z| > 60 が漸近展開 (`expn_cn`・`_gamma_upper_ca`)。mpmath と 340 点 (|z| 〜 1e4、負の実軸の上下を含む) で 1e-9、gs と `PolyLog[1.5+I,.3]`・`ExpIntegralE[1.5+I,2.]` が 1e-14。**SADError で止まる領域は無い** (非整数でない負の整数の次数 ν = 0, −1, … の ExpIntegralE は整数次の `expn_c`、PolyLog の整数次は従来の経路)。
**|z| > 35 の複素の走査 (613 点、|z| = 36〜1000)**: 誤りは無かった (Sin[1000 I]・Cos[1000 I] の桁あふれのエラーだけ。gs と同じ INF 成分の複素数にした)。mpmath の多倍長のフォールバックは要らないので入れていない (mpmath を依存に加えない)。
**表示**: 決まらない比較は `x<y`・`x>1&&y>1`・`x<y<z` (= `1&&x<y&&y<z`) と gs と同じに出す (従来は `$Cmp[x,y]`。k64・v2 とも)。

### v2 で台本の結果が変わるもの (E。影響の出る書き方と対処)

| 変更 | 影響が出る台本の書き方 | 対処 (k64 の値が要るとき) |
|---|---|---|
| `Round` の半整数を偶数へ (B184) | `Round[2.5]`・`Round[x+.5]`・`Round[n/2]` で小数部がちょうど .5 になる値 (丸め先が 3 → 2、0.5 → 0、−2.5 → −2)。間隔を丸める `Round[x,d]` も半整数のとき | `PYSAD_COMPAT=k64`、または `Floor[x+.5]` と書く |
| `Insert[list,x,-n]`・`StringInsert[s,x,-n]` (B197) | 負の位置。`Insert[l,x,-1]` が末尾 (k64 は末尾の 1 つ手前) | 正の位置 `Length[l]+1` と書く |
| `x==y`・`If[x,..]` が決まらない (B199 は k64 も) | シンボル・未定義の値を `==`・`If` の条件にする台本。従来は 0・Null (If の else) と評価されたが、gs と同じく式のまま返る。`x==x` は v2 で 1 (k64 は式のまま) | 値が決まってから比較する。式の同一判定には `===` を使う |
| `And[1,x]`・`x&&True` が x (B199 は k64 も) | 未定義のシンボルを `&&`・`||` に使う台本 (従来は 0) | 同上 |
| 虚部 0 の複素数が実数 (B199 は k64 も) | `(1+I)(1-I)`・`I^2` を `Re`・`Im`・`NumberQ[..]`・`Head[..]` で調べる台本 (従来は Complex、いまは Real) | — (gs と同じ) |
| `Table`・`Do`・`Sum`・`Product` の `{i,{a,b,c}}` (B198) | k64 ではエラーだった書式が通る (結果が変わるのではなく通るようになる) | — |
| `Mod[x,y]` が正確な剰余 (B184) | x/y が 1e15 を超える・x が y の倍数に近い値 (最後の桁と、0 か y に近い方の判定が変わる) | — |
| 数学関数の定義域外の値 (B180〜B182, B185〜B196) | `ArcCosh[x]` (|x|<1)・`ArcSin`/`ArcCos` (x>1) の虚部の符号が反転、NaN だった `ArcTanh[|x|>1]`・x<0 の `BesselY`・`HermiteH` の複素・非整数次の `LegendreP[ν,m,x]` が数になる | 符号を使う台本は `Conjugate[..]` で戻せる。k64 は gs の値 |

**v2 でも直さないもの (理由つき)**: (1) `Hold` の部分評価 (gs は引数の算術を評価する。規則が数個の例から決められない。pysad は標準の HoldAll)、(2) `LINE["POSITION",0.5]`・範囲外の位置 (gs が範囲外の位置を作る。数学的な定義が無い)、(3) `Take[{1,2,3},4]`・`Drop[l,5]` が `{}` を返す (Mathematica は `Take::take` のエラー。gs・pysad は `{}`。直すと SAD のパッケージが依拠する可能性があり、独立な基準は無い)、(4) `FractionalPart[-2.5]` = .5 (SAD の Packages/functions.n の定義 x − Floor[x]。Mathematica は −.5)、`FixedPointList` が初期値を含まない (SAD の HelpMessages.n の定義)、(5) `Pochhammer[-.5,-.5]`・`Factorial[-1.]` の ComplexInfinity (pysad は NaN)、(6) (実装した: 複素の次数の PolyLog・ExpIntegralE)、(7) (直した: `$Cmp` の表示)、(8) 末尾の `1+` は v2 が構文エラーのまま (gs は値 1。入力ミスを隠す癖)。

**未確認**: Mathematica 本体は無く、「Mathematica の定義」は mpmath の主値と公開されている定義 (切断線の取り方) による。複素数の 1F1 の十進の級数は |z| ≲ 700 まで (項数 20000 の上限)。複素の特殊関数の |z| > 35 の領域は走査していない。

## v2.4.0 要判断項目の一覧 (v2.4.0 の初期の状態。その後の変更は各節 (D・D2・D3・E・F・G・H) が最新)

docs/47 v2.2.0 の (d) の「要判断」と、v2.4.0 で足した判断を 1 項目ずつ並べた。この表を作った時点の推奨はすべて「今の既定を据え置く」だった。**その後、v2 は K64 の誤りを直す方針に決めたので、B∥ (行 8)・T[z, δ] (行 11)・β の因子 (行 5)・Touschek のバンチ内 δ (行 15)・ExpandElementValues (行 16)・`Bend.transfer_matrix` (行 17) などは v2 で修正した** (表の各行と、D・D2・D3・F の節を参照)。「影響」は変えたときに動くもの。

| # | 項目 | 現状 (既定) | 選択肢 | 影響 (変えたとき) | 推奨 | 根拠・参照 |
|---|---|---|---|---|---|---|
| 1 | 複素の特殊関数の k64 の値 (B158・B118〜B120。**v2.4.0 E: v2 の真値を mpmath の 1871 点で再確認、複素の 1F1・Shi・Chi を B186・B190 で追加修正**) | k64 も真値 (gs の誤りを再現しない)。732 式の走査で gs と 1e-9 で一致 641 式、残り 91 式は gs が誤る・不正確 | (a) 据え置き (b) k64 だけ gs の誤りを再現: 規則が単純な群 (LogGamma の分枝・ArcCosh の符号・Erf/Erfc の Re z < 0・ExpIntegralEi) から K64 の `gamma.f`・D.n を直訳 | (b) は k64 で gs と 91 式のうち規則の単純な分が一致 (回帰のベンチは複素特殊関数を使わないので既存試験は動かない)。v2 は変わらない | **据え置き**。k64 の目的は gs との一致だが、誤った値を基準にしても使う人がいない。必要な式が出てから (b) を式ごとに | docs/46 の特殊関数の行、`v2next_P5` |
| 2 | `DAPWIDTH` の既定 | 7 (gs と同じ。v2.1.0 は 51) | (a) 7 (b) 51 に戻す | DAPWIDTH を書かない台本の DA の Score と窓の外の `*` が変わる。51 だと全点を追うが遅く、gs の既定の出力と違う | **7 のまま**。全点を追う台本には `DAPWIDTH = 51;` を書く (DDBA の例の `DAsv` は 51 を設定済み) | `v2next_P4`、[`test_v2_P4_dapwidth.py`](../pysad/tests/test_v2_P4_dapwidth.py) |
| 3 | `WSPACSC` | 0 (追跡のビーム行列はゼロ電流の平衡行列。gs の追跡と 1〜2%)。**v2.4.0 F2 で検討し据え置き**: `WSPACSC = 1` の行列は MINCOUP の下限で垂直のビームサイズを人為的に大きくする (B149 の例: MINCOUP の下限が自然な結合より大きいとき垂直のチューンシフトが 1/3 になる。−1.2e-3 対 −3.8e-4)。解析式の空間電荷チューンシフト (ラスレット) はビームの大きさを与えたときの式なので、どちらの行列が実際のビームかは判定できず (どちらの行列でも一次の式と 2e-4 で一致)、独立な基準は gs の追跡だけで、それはゼロ電流の行列と 1〜2% で合う。MINCOUP の下限は数値の正則化であって物理ではない | (a) 0 (b) 既定を空間電荷つきの行列にする (K64 が `EMIT` の beamsize を使う意図) | 追跡の WSPAC の値が 1〜2% 動く。`WSPACSC = 1` は K64 にないオプションで、gs の追跡とは合わない (MINCOUP の下限が入る) | **0 のまま**。gs との比較ができる値を既定にする | [`test_v2_D1_emit.py`](../pysad/tests/test_v2_D1_emit.py)、`v2next_P4` |
| 4 | EMITDIV の既定 | v2.4.0 で v2 の既定を 4 に変更済み (k64 は 1)。k64 の既定は変えていない | (a) 4 のまま (b) 1 に戻す (c) 2 | 4 → 1 で矩形磁石・低エミッタンスのリングの Emittance[] の εx が +1.4e-3〜+2.1% (SP8 で 1.4e-3、大きなリング A で 2.1%)。計算時間は EMITDIV に比例 | **4 のまま** (解決済みの項目。確認だけ)。k64 は gs との一致のため 1 | docs/46 の RADINT の行、`v2next_P2` |
| 5 | 放射の β の因子 (B112) | **v2 で修正 (v2.4.0 D、B167。以下は修正前の記述)**: 直していない (v2 の U0 は古典式の 1/β³ 倍: 陽子 1.5・3・30 GeV/c で 64%・15%・0.15%) | (a) 据え置き (b) 光子数・エネルギー分布の β 依存を厳密な放射のスペクトル (β < 1) で入れ、独立な基準 (厳密な放射の式の積分) と比べて直す | 低い γ の陽子の放射 (U0・εx・σδ・Touschek) が変わる。電子では β = 1 で変わらない | **据え置き**。電子のラティスしか使わない間は影響がない。陽子の放射が要るときに (b) を実装。それまでは低い γ の陽子で放射を使わない | docs/46、`v2next_P5` |
| 6 | K1≠0 の BEND の端の模型 (`face`) | `radial` (K64 の仕様。k64 は常に `radial`) | (a) `radial` (b) v2 の既定を `face` (厳密な硬い端の磁石) にする | 軌道が中心のリングでは同じ。軌道がずれた BEND の x0 tan e の項 (1e-8〜1e-7) | **`radial` のまま**。`face` は DISFRIN=1・ANGLE=0・SOL 内・RADINT・F1/FB の端で使えないので、既定にすると使えない場合のエラー分岐が増える。必要な利用者が `compat.set_bend_k1_edge('face')` を選ぶ | docs/46 の RADINT の行、`v2next_P2` |
| 7 | SOL 領域内の端の縦磁場 ψ (B146) | DISFRIN=1 で入れない (v2)。DISFRIN=0 では入れる | (a) 据え置き (b) DISFRIN=0 の ψ が物理的に正しいかを、硬い端の B_s と軌道のフリンジ写像を合わせた厳密な基準で確かめてから決める | SOL を含む RADINT・POL・追跡のフリンジ項 | **据え置き**。独立な基準が未着手で、どちらが正しいか決められない。DISFRIN=1 で ψ を使わないのは K64 の仕様 | `v2next_P4`、docs/46 の SOL | **[v2.4.0 D3] 解決: DISFRIN=0 の ψ は物理的に正しい (本体の Δψ が欠けていた = B234。本体を足して Maxwell 場の積分と一致、DISFRIN=1 は硬い端の積分と一致)。**
| 8 | 偏向磁石の本体の B∥ (B136・B137) | **v2 で修正 (BEND の本体のみ、v2.4.0 D、B165。QUAD・MULT 等は未対応。以下は修正前の記述)**: 直していない (k64・v2 とも K64 のまま) | (a) 据え置き (b) tradke・tradk・tbend の全ての要素の本体に (1+a) B∥ ds を足しフレームの扱いを揃える | B58・B64・B130 のスピンの結果が全部動く。脱偏極率の追跡との食い違い (鉛直の軌道・分散の大きいリング) が減る見込み。**スピン・偏極は B の担当で、この表の他の項目と違い影響が広い** | **据え置き**。基準 (PF・大きなリング A の FLUC のスピン追跡、1 条件数時間) が揃うまで動かさない。揃ったら (b) を v2 だけで | docs/46、`v2next_P1` | **[v2.4.0 D3] QUAD・MULT・SOL 内の四極も v2 で対応 (B230・B231・B234)。SEXT・ステアリングは元から正しい。追跡の BEND の R12・R21 の 1e-2 は刻みの 1 次誤差 (B235)。**
| 9 | `CSROYMatrix` の gs の wrongtype (B60) | 再現しない (pysad は `CSRSetupOY` の出力を受ける) | (a) 据え置き (b) k64 で再現 | 再現すると k64 で `OYStability`・`OYScan` が使えなくなる | **据え置き** (使えなくなる側を選ぶ理由がない) | `v2next_P5` |
| 10 | 文法の寛容さ (**v2.4.0 E: 確認済み。v2 のまま**) | `45DEG` は v2 が受け k64 は構文エラー。`1+` は v2 が構文エラー・k64 は値 1 | (a) このまま (b) v2 も gs に合わせる | 台本の書き方の許容。(b) は v2 で動いていた入力が止まる | **このまま** (v2 の寛容さは害が少ない。末尾の `1+` を値 1 にするのは gs の癖で、入力ミスを隠す) | docs/46 の台本の読み方 |
| 11 | BEND の T[z, δ] の近似 (B139) | **v2 で修正 (v2.4.0 D、B164。以下は修正前の記述)**: 直していない (K64 の仕様。1 周で 1e-7) | (a) 据え置き (b) v2 で `tbendebody0` の (5,6) の px² の項を厳密にする | 6D の z–δ の行列が 1e-7 動く (RADCOD の 1 周の解析 T の残り)。k64 は変えない | **据え置き**。影響が 1e-7 で、直すと v2 の Emittance[]・RADCOD の基準値が全部動く | `v2next_P2` |
| 12 | `tsteer` の F1 の fb1 ≠ fb2 (B145) | 直していない (K64 の模型の近似) | (a) 据え置き (b) 磁場の立ち上がりの中心を要素の端に置く取り方を決めて直す | F1 つきのステアリングの追跡と CALC の微小な差 | **据え置き**。取り方の判断に独立な基準 (厳密なフリンジ場) が要る | `v2next_P3` | **[v2.4.0 D3] 再分類: 誤りではなく基準軌道の取り方 (設計軌道 δ = 0 をフリンジが動かさない)。独立な厳密な運動との差は δ に依らない定数 (K0/al)(fb1² − fb2²)/24 だけ (B236)。変えない。**
| 13 | K1≠0・K0≠0 の ANGLE=0 BEND の追跡・CALC の AE (B142) | (v2 はすでに修正済み。k64 は変えない方針) v2 は AE を使う (ソースの読み)。k64 の追跡は ψ = 0 を渡し AE を無視 (CALC と 1e-4 以上違う) | (a) 据え置き (b) k64 の追跡も AE を使う | k64 の追跡が CALC と整合するが gs と違う | **据え置き**。gs が SIGSEGV で K64 の値を確認できないので、読みが正しいかが未確認のまま | [`test_v2_P3_ae.py`](../pysad/tests/test_v2_P3_ae.py) |
| 14 | SOL 領域内の加速 (POL + TRPT) | 実装せず明示のエラー | (a) エラーのまま (b) RF の電場を含む運動方程式の積分の独立な基準を作り、`_RadHook` の基準運動量の更新と `_tradk_sol` の p(s) を実装 | SOL を含むリングの加速の追跡が可能になる | **エラーのまま** (必要な利用者が出てから) | docs/46 の POL の行 |
| 15 | (v2.4.0 追加) Touschek のバンチ内の δ の分布 (`TouschekLifetimeLMA[BunchDelta->True, BunchQs->True]`) | **解決 (v2.4.0 F2、B200): v2 の既定を「入れる」に変更** (k64 は入れない = δ_old = 0。K64・elegant・Piwinski と同じ基準)。`->False` で切れる | (a) 入れない (b) 既定で入れる → **(b) に決定** | 損失率が増え寿命が短くなる: ASTRID2 a/σδ = 9 で +4.2〜4.5% (δ_c + q_s)、14 で +1.9%、23 で +0.7%。elegant との比較 (T1・T2) は「δ = 0 の基準との差」になる | **(b)** (誤りの指摘は B200: δ_old = 0 の基準はバンチ内の δ の分布を無視する近似で、独立な基準 (実験室系の対のモンテカルロ) が入れた値と 0.3〜2% で一致)。受容が σδ の 8 倍を超える範囲でだけ有効 (警告あり) | docs/48 §4.4、[`docs/history/v24_A.md`](history/v24_A.md) A3・A10 |
| 16 | (v2.4.0 追加) `DynamicApertureSurvey` の `ExpandElementValues` の既定 | **解決 (v2.4.0 F2、B201): v2 の既定を False に変更** (k64 は True = gs と同じ。成分ごとの K0 などを代表の値に揃える) | (a) True (b) v2 の既定を False → **(b) に決定** | 補正した状態 (成分ごとの K0・位置ずれ) を持つリングの DA が、補正のない状態の値から補正後の値に変わる (DDBA v07_14: 0〜38σ の偶然の値 → 29〜36σ。小さな FODO で 15 → 51)。gs の既定と違う | **(b)** (誤りの指摘は B201: 既定 True は補正磁石の成分を消す)。k64 は gs と同じ。`->True` で従来どおり揃える。台本に書く必要は無くなった (DDBA の例の `DAsv` は明示のまま有効) | [`docs/history/v24_A.md`](history/v24_A.md) A5、[`test_v24_A5_dapert_options.py`](../pysad/tests/test_v24_A5_dapert_options.py) |
| 17 | (v2.4.0 追加) `Bend.transfer_matrix` の K0≠0・ANGLE≠0 | 直していない (CALC・追跡・resdrive の経路では使われない API。K0 = 0 では 5e-16、K0≠0 で追跡のヤコビアンと 5e-5〜2e-4) | (a) 据え置き (b) 場の角 (ANGLE + K0) と端・K1 の組み合わせを合わせて直す | `transfer_matrix` を CALC の外から直接呼ぶ場合 (K0≠0 かつ ANGLE≠0) の値 | **据え置き** (呼ぶ経路が無い。emit/radint の B のファイル) | [`docs/history/v24_A.md`](history/v24_A.md) A12 |

追記 (v2.4.0 担当 A の第 5 弾): USE のたびに DP0 を 0 に戻す (gs と照合。k64・v2 とも)。SADScript の評価器を、結果を変えずに約 15% 速くした (DDBA v07_14 の誤差計算の例 93 s → 80 s。出力は全行一致): `_copy`・`_has_complex` を型の集合 (C の速さ) で判定し、最も多い式ノードを `eval` の型判定を通さず `eval_expr` へ直接渡す。残りの時間は追跡 (numba、約 1/3) と評価器の再帰 (式 1 個あたり約 100 μs) で、これ以上は評価器の構造 (束縛・パターン・保持引数の判定が 1 回の eval ごとに入る) を変えないと縮まらない。

## v2.4.0 F (ROTATE した偏向磁石の RDT・ADTS、DA 関係の v2 の既定の見直し)

記録は `docs/history/v24_F.md`、試験 `test_v24_F1_rotated_bend.py`・`test_v24_F2_defaults.py`、解説は docs/49 §14。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| ROTATE した BEND・QUAD の RDT・ADTS・ν(J, δ) の級数 | 同じ (新機能。互換モードによらない) | `resonance_terms`・`adts` が結合した光学 (4×4 の固有モード) の模型で通す (全次数の θ)。以前はその磁石を飛ばして警告 (BEND)。ampnf は以前 SEXT・OCT の ROTATION を黙って無視していた (pysad 自身の不具合を修正) | 同じ模型の追跡 3e-8、実際の追跡 3e-6 (ν)、追跡の 1 周の写像の 2 次のテイラー展開 6e-4、ROTATE = π = 逆曲げ 1e-9 (docs/49 §14) |
| Touschek のバンチ内の δ の分布 (B200) | 入れない (δ_old = 0) | **既定で入れる** (`TouschekLifetimeLMA[]`・`TouschekLifetime[INF, INF, nz]` の Piwinski 型。BunchDelta・BunchQs。`->False` で切る) | 実験室系の対のモンテカルロ (δ_c 2%、q_s 込み 0.3%) |
| `DynamicApertureSurvey` の `ExpandElementValues` の既定 (B201) | True (gs と同じ) | **False** (補正磁石の成分を消さない。`->True` で従来) | 小さな FODO で DA 15 → 51 (消える)、DDBA v07_14 で 29〜36σ → 偶然の値 |
| FrequencyMap・MomentumAcceptance・NormalFormTune の TrackSetup への受け渡し (pysad 自身) | 同じ | NOCALEXP (`NOEXPAND`) のとき成分ごとの値を揃えない (以前は無条件に揃えて補正後の状態が消えた)。CALEXP (既定) は従来どおり | `test_v24_F2_defaults.py` |
| `WSPACSC` の既定 | 0 | 0 のまま (要判断 3 の理由) | docs/47 の表の行 3 |

**既定を変えたことで動く既存の参照値**: v2 モードの Touschek 寿命 (`TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の Piwinski 型) は δ_c つきの値 (ASTRID2 a/σδ = 9 で −2.5% (δ_c のみ)・−4.1〜4.3% (q_s 込み)) になる。`test_v23_T2_touschek.py` の v2 の既定の試験は `BunchDelta->0` を明示して δ_c = 0 の値を確認する形に変えた。elegant との比較 (T1・T2) の値は δ_c = 0 の基準のまま (`example/v2_3/T1_*`・`T2_*` は `touschek_piwinski` を直接呼ぶ)。k64 モード・gs との比較は変わらない。

**誤りの指摘**: B200 (Touschek の δ_old = 0 はバンチ内の δ の分布を無視する近似)、B201 (DAS の既定は補正磁石の成分を消す)。K64 の誤りと言えるが gs の既定を再現する k64 は変えていない。

**要判断・未解決**: `WSPACSC` の既定 (据え置き)、K0≠0 の偏向磁石に ROTATE がある場合 (ampnf が K0≠0 を扱えない。docs/49 §14.4)、結合した経路の 4 次の h・δ≠0・ref≠None。

### v2.4.0 G (K0≠0 の偏向磁石 = 挿入光源・逆曲げ・MPW・傾けた MPW の ν(J, δ)・RDT・ADTS)

記録は `docs/history/v24_G.md`、試験 `test_v24_G_k0_ampnf.py`、解説は docs/49 §14.4、例 `example/v2_4/G_k0_ampnf.py`。新機能 (互換モードによらない)。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| ampnf (`amplitude_tune`・`one_turn_series`・`rdt_series`・`model_tune`・`tracked_tunes`) の K0≠0 の BEND | 同じ (新機能) | 入る (以前は `NotImplementedError`)。K64 の追跡の閉じた式 (`tbendcore1`・`k_steer`) と積分器 (`tbendi`、K1≠0) を磁石内の軌道のまわりで TPSA 評価 (`pysad/ampnf_k0.py`)。ANGLE = 0・≠0、E1・E2、逆曲げ、ROTATE (傾けた MPW)。閉軌道が原点にない場合は 1 周の写像の不動点を求めて展開 (`res['orbit']`) | 素子の写像が追跡と 1e-13、ν(J, δ) が実際の追跡 + NAFF と 4e-6 以内 (3 極の挿入光源)、h が追跡の 2 次のテイラー展開と 1.5e-3 以内 |
| `resonance_terms`・`adts` の傾けた K0≠0 の BEND | 同じ | 結合した経路・`.rotate` に自動で入る (以前は飛ばして警告)。ロールなしの K0≠0 は A2 の解析式のまま | 傾けた (0.2 rad) 3 極の `adts().total` が追跡のフィットと 3e-3 以内 |

**pysad 自身の変更**: A2 の試験 `test_rotated_k0_bend_is_not_included_and_warns` は、傾けた K0 の磁石が入るようになったので `test_rotated_k0_bend_is_included_in_coupled_path` に変えた。

**G の続き (K1≠0・F1・据え付け・独立な基準・実機規模)**: K1≠0 かつ K0≠0 は K64 の `tbendi` の数値積分器そのものを TPSA で評価、F1・FB1・FB2、DX・DY・DROTATE・CHI2、場が 0 の ANGLE + K0 = 0 (K64 は phib = 0 で割って追跡が NaN。phib → 0 の極限で入れ、独立な自由な直進と 1.4e-8) も入れた。**独立な基準** (実験室の直交座標での運動方程式の数値積分、Maxwell を満たす端): K64 の `tbendcore1`・`k_steer` が水平面 1e-12、Maxwell の項も w → 0 で一致 (K1 = 0)。見つけた K64 の近似: `tsteer` (ANGLE = 0) の楔は薄い蹴りだけ (v2 で 3e-7。放射担当の B221 で v2 は厳密な鎖になり 5e-13)、k64 は出口の楔を入れない (B140、9e-5)。ANGLE≠0 は厳密。`bend_cubic`・`kin`・`fringe` は、どれかが True なら K0≠0 の磁石は全部入り、全部 False なら軌道のまわりの 1 次だけ。実機規模: ASTRID2 + 2 T MPW 23 極 (D = 5 で 10 s)、KEK-LS + MPW 23 極 (144 s) で ν が追跡と 9e-7 以内。h10002 の 1.3e-3 は D = 3 の打ち切り (運動学 4 次のフィードダウン。D = 5 で 3.5e-5)。

**pysad 自身の不具合 (H の報告)**: `ampnf.amplitude_tune` が素子の DX・DY を黙って無視していた → DX・DY・DZ・DROTATE のある lat は H1 の `ampnf_coupled.amplitude_tune_ext` に回す。K0≠0 の磁石は結合した経路でも G の `tmap` で通る (`ampnf_coupled.build_ops_ext`・`propagate`・`numeric_track` に最小の変更)。

**G の残りの項目 (2026-11-22 夜)**: (1) ANGLE = 0 で K1≠0 かつ K0≠0 (K64 の tmulti。gs は SIGSEGV) は、厳密な Hamiltonian の対称な 4 次の分割の TPSA (Lee–Whiting + 双極子の Maxwell の端つき) で入れた。独立な運動方程式の積分と w² で収束し、pysad の `tmulti` の追跡は 4e-7〜1.6e-5 外れる。ANGLE=0 の CHI2 は追跡に合わせて無視 (警告)。(2) K1≠0 の曲がった磁石の本体 (水平面) の独立な基準: K64 の `tbendi`・K1 の 'log'/'lin' の定義が 1e-9 以内で正しい。(3) **`rdttrack` の修正** (pysad 自身の不具合): `tracked_rdt_table` が展開の中心を原点にしていて、閉軌道 ≠ 0 のリング先頭で h10002 が壊れていた (6.6e3 対 0.55) → δ = 0 の閉軌道を中心にした。(4) ASTRID2 の δ の差 1e-5 は偏向磁石の分割数 (nslice_bend)。(5) **B250**: tbendi の xsxkx の符号 (akxsq < 0 の枝が極限と逆。行路 z だけ ~1e-9 m) を v2 で直した (追跡・Emittance[] の python 経路・CALC の numba 版 (G7 で追加)。k64 は不変)。

**G7 (2026-11-22 夜)**: B250 を CALC の numba 版 (`perfkern_bend.bendebody`、`compat.bend_ecorr_flags()` の bit2) にも入れた。K1≠0 の曲がった磁石の垂直方向・端の独立な基準 (多極展開の端 + Maxwell を満たす場の運動方程式の積分) を作り、K64 の `tbendi` (本体・Maxwell の端・Lee–Whiting・楔) が E なしで w² で一致、E≠0 は v2 の 'face' で一致、K64 既定の 'radial' は E≠0・K1≠0・K0≠0 で 3.5e-7 (面の取り方の違い)。ANGLE≠0 の K1≠0 かつ K0≠0 の F1・FB を `tbendi` の fb の項で入れた (追跡と 1e-16)。

**未対応・要判断**: ANGLE=0 で K1≠0 かつ K0≠0 の F1・FB1・FB2 は `NotImplementedError` (tmulti。gs が SIGSEGV で基準なし)。

## v2.4.0 H (結合した経路の RDT の拡張: 4 次の h・δ≠0・ref・DX/DY/DROTATE、KEK-LS の補正つき LMA)

記録は `docs/history/v24_H.md`、試験 `test_v24_H1_coupled_h4.py`、解説は docs/49 §14.5・docs/48 §11、コードは新しい `pysad/ampnf_coupled.py` (`ampnf.py` は編集せず、`resdrive.py` は呼び出しの追加)。F の「要判断・未解決」の最後の項目 (結合した経路の 4 次の h・δ≠0・ref≠None) は H で入った。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| ロール・ずれのある lat の 4 次の h・δ≠0・ref | 同じ (新機能。互換モードによらない) | `resonance_terms` が結合した経路で 4 次の名前・δ・ref を受ける (以前は NotImplementedError → 従来の式でロールを飛ばす)。h4 = 写像 e^{−:h3 + h4:} の h4 | 追跡の写像のテイラー展開と h3 7e-6〜3e-4・h4 2e-4〜8e-4 (結合 0.01〜0.24)、ADTS が ampnf の正規形と 1e-13 |
| DX・DY のある lat の `resonance_terms` (pysad 自身の不具合) | 同じ | 以前は DX・DY を黙って無視していた (解析式は閉軌道 0・ずれなしの前提) → 結合した経路 (閉軌道 ≠ 0 のまわりの展開) に回す (警告つき) | 閉軌道が追跡と 1e-10〜1e-8、h3・h4 は上の水準 |
| DX・DY・DZ・DROTATE の TPSA | 同じ (出口の符号は幾何として正しい符号に固定。K64 の追跡 tbshift も同じ符号で、k64 モードの追跡と 1e-13。K64 の CALC・Emittance の tchge の出口の符号は逆 = A2 が v2 で直した誤り) | SAD の tchge (入口・出口、tbrote の 3 次元の回転) をそのまま TPSA にした | 要素 1 個の写像が追跡と 2e-10 (v2・k64 とも) |
| `adts().rotate` の δ≠0・D | 同じ | δ≠0 は J¹ の δ 依存 (ν の級数 D=5、p ≤ 2) で評価 (以前は δ=0 の値)。傾きの項の D を 3 → 4 に上げた (D=3 では運動学・フリンジが TPSA の次数に入らず落ちていた)。DX・DY・DROTATE のある lat は閉軌道のまわりの正規形 | δ = 1e-3 で δ=0 の値から 4〜13% 動き、δ の閉軌道のまわりの写像の h の ADTS と 1e-3 以内 |
| `MomentumAcceptance[]` に補正後の状態 (F2 の修正) の KEK-LS での実走 | 同じ | NOEXPAND の台本で補正後の K0 が TrackSetup に渡る (KEK-LS v03_68・v03_12 で確認)。既定の CALEXP は補正が消える | docs/48 §11 |

**結果が変わるもの**: (1) ロール・ずれのある lat の `resonance_terms(ref=, delta=)` と 4 次の名前が、以前の例外・従来の式から結合した経路の値になる。(2) DX・DY だけのある lat の `resonance_terms` の値 (以前は無視)。(3) `adts().rotate` の値が D=4 の模型の J¹ の差になる (傾きと運動学・フリンジの交差が入る。差は ADTS の 1e-3 以下)。

**v2 の既定の変更 (H4、B260)**: v2 の周期閉軌道の探索 (CALC・Emittance・追跡の初期軌道) は、K64 の qcod が収束しなかったときだけ堅牢な解法 (`lattice._qcod_robust`) にフォールバックする。収束していた計算は変わらない (ビット単位)。収束しなかった計算 (誤差の大きい低エミッタンスのリング) は、以前は不動点でない軌道や設計の光学で続行していたが、真の不動点 (安定なものを優先) を返す。k64 は不変。`lattice._QCOD_ROBUST = False` で v2 でも切れる。

**要判断・未解決**: 結合が強い (0.7 以上) ASTRID2 の DROTATE の h の 2e-3〜4e-3、誤差つき LMA の周回数 (2000 周以上で見る。1 万周以上・FLUC は未確認。docs/48 §11)、`_qcod_robust` の収束の保証 (docs/46 の H の節)。閉軌道 ≠ 0 の ADTS の NAFF のフィットの J の格子は整理済み (定数項つき `tracked_adts_fit`、H4-2)。

### v2.4.0 D (放射・スピン・偏向磁石の物理)

- **B164** (B139): `tbendebody0` の dω/dp を厳密にした (v2)。1 セルの解析 T と写像の数値微分の差 6.4e-9 → 4.5e-11。k64 不変。
- **B165** (B136・B137): BEND (K1=0・K1≠0) の本体に (ψ_終 − ψ_始)/(1+δ) を運動量の向きの磁場として足した (v2)。Maxwell 場つきの硬い端の Thomas-BMT の積分が基準 (dR/dpy の誤差 0.26 → 8e-4 (追跡 R01)・4e-4 (Emittance[])。B136 の「(2+a)θ py」のうち θ py は端の ψ だった)。QUAD・MULT・SOL 内の四極は v2.4.0 D3 で対応 (B230・B234)、SEXT・ステアリングは元から正しい (和が 0)。
- **B167** (B112 の β): 陽子の U0 は古典式と 1e-6、εx・減衰率・σδ は r2/β⁴・β² の因子 (`pysad/radbeta.py`)。RADINT も陽子で電子の定数を使わなくなった (v2)。γ ≥ 300 (電子の全部) には入れないので電子の結果は動かない。(v2.4.0 D3 で確認: 独立な (x, px, δ, z) 写像と εx 4e-5・σδ 8e-6、RF の同期位相 B232、PhotonList の β B233)。
- **v2.4.0 D2 で閉じた**: 項目 4 (B143/B140 の端の幾何 → B220・B221・B223)、5 (CALC の分散の K0 の発生 → 誤りではない、B222)。項目 6 (B145・`Bend.transfer_matrix`・B146 ψ・SOL 内の加速) は **v2.4.0 D3 で閉じた (下)**。理由と次の手は docs/46。



### v2.4.0 D2 (硬い端の極の端の幾何・CALC の分散)

試験 `pysad/tests/test_v2_D2_edge.py` (17 件)、独立な基準 `_d2_indep.py`・`rad_indep_d2.py` (硬い端の面を実験室の直線座標で厳密に積分する写像)。バグ表は docs/21 の B220〜B222。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| B220: 端の角の放射の「場の中の行程」(`tin`・`tout`、追跡の `zr0`) | tan(e)·x の 1 次 (v2.2.0 で v2 に入れた形) | 面までの直進の行程を厳密に: ℓ = x tan e/(1 ∓ (px/pz) tan e)·p/pz (CALC・Emittance[]・追跡の全ての端の角。`perfkern.edge_free`・`edge_len`) | 1 つの極の放射つき写像のヤコビアンの δ の行が厳密解 (RK4) と 1e-12 以下 (以前 1% = 8e-9)。リングの減衰率が AE = 0.1 で 5e-5 → 4e-6 |
| B221: 追跡 `tsteer` の端 | 薄い蹴り x tan(e)/ρ (K0³ tan²e で外れる) | CALC の `tsteee` と同じ厳密な硬い端の鎖 (F1・FB1・FB2 なしのとき)。python 版・numba 版 (`track6.py`・`track6n.py`) とも | 1 要素の写像が厳密解と 1e-13 以内。FODO + ステアリング 2 台 (AE = 0.3) で、追跡と CALC のチューンの AE に依る差が 2e-7 以下 (以前 K0 = 0.0015・0.003・0.006 で 6e-6・2.5e-5・1e-4) |
| B222: CALC の分散 | 変更なし | 変更なし (誤りではない) | EX・EPX が追跡の閉軌道の δ 微分と一致、δ 列が追跡のヤコビアンと一致 (docs/21 B222) |

**結果が変わるもの (v2 のみ)**: 端の角 (E1・E2・AE1・AE2) のある BEND・MULT の放射 (U0・εx・減衰率) が px tan²e の次数 (典型的な蓄積リングで 1e-6 以下) で動く。K0 の大きい極 (0.06 rad) では減衰率が 5e-5 動く。ANGLE=0 の BEND の追跡 (tsteer、AE ≠ 0) の写像が K0³ tan²e の次数で動く。

**訂正 (B143 の「残る 2.7%」と v2.4.0 D の「K64 の端の模型 tbedge の AE の 2 次の違い」)**: 原因は端の模型でなく基準側だった。(1) 基準のリングに四極のハードエッジの非線形フリンジ (tqfrie、SAD の既定) が無く、AE で軌道が四極の中心から外れるので AE² で効く (チューンの差 8e-5・9e-4)。(2) 基準の U0 が RF なしの軌道上の値で、閉軌道のエネルギーの沈み込み (δ = −8.8e-4) の分 0.2% 違う。damping rate を直接比べると、四極のフリンジを切った厳密解と AE = 0・0.1・0.3 で 3e-6・4e-6・1.2e-5 で一致する。CALC の端の鎖 (tbedge) は硬い端の厳密な写像と 1e-16 で一致している。

**第 2 弾 (D2 の続き)**: (1) B223: MULT の端の角 (tmulke の薄い蹴り) を、入口 Φ(−z_f)∘Drift(z_f)・出口 Drift(ℓ)∘Φ(−ℓ) の厳密な硬い端にした (v2。CALC・Emittance[]・追跡 tmulti)。K1 = 0 は厳密解と 1e-15、K1 ≠ 0 は 3e-7。MULT 3 極のリング (K0 = 0.06・−0.12・0.06、四極のフリンジなし) は厳密解とチューン 2.5e-6、減衰 1e-5 (AE = 0.1)・5.5e-5 (AE = 0.3) (以前 2.9e-4・2.5e-3、2.6e-4・3.2e-3)。(2) B221 の続き: tsteer の F1・FB1・FB2 つきの端も厳密な鎖 (CALC と同じ順序 tblfre → 鎖 → 断片)。FRINGE = −1・−2 (片側) は除く。(3) 追跡と CALC のチューンの AE に依らない K0² の差の原因: **四極の本体の積分器** (K64 の tsolqu と tsolque)。軌道が 3 mm・2 mrad ずれた四極 1 台のヤコビアンが追跡と CALC で 5e-7 違う (厳密な写像との差: CALC 1.8e-7、追跡 7e-7)。ステアリング・扇形の偏向磁石・ドリフトは 1e-9 以下で一致。K64 の積分器の刻みの限界で、v2 でも直していない。

**結果が変わるもの (追加)**: 端の角のある MULT の v2 の CALC・Emittance[]・追跡 (K0 の大きい極では 1e-5〜1e-3、典型的な補正磁石・K0 ≤ 1e-3 では 1e-10 以下)。

**第 3 弾 (D2 の続き)**: (1) B224: K1 ≠ 0 の MULT の端 (K1 < 0 は面の傾きが実験室系の y): 原因は 2 つ。(a) 厳密な硬い端の写像はシンプレクティックでなく (K1 ≠ 0 で 5e-4〜3e-2)、`tradke` の D・B の更新のシンプレクティックの逆行列の式が効かない → v2 は厳密な逆行列 (`perfkern.tradke_db`)。(b) MULT の厳密な端に縦収束 (Maxwell) が無かった (BEND の鎖と y・py が 8e-4 違った) → K64 の薄い蹴りの縦成分を足した (2e-5 以内)。δ の行が数値微分と一致 (1e-17)、K1 < 0 の MULT の Emittance[] の ΔJx −0.011 → 8e-6 (追跡と一致)。(2) B225: FRINGE = −1・−2 のステアリングの追跡が、CALC と同じく反対側の端を飛ばす (v2。接合部の端で縦収束・面を入れないのが意図)。追跡と CALC のチューンの差は AE に依らず 6e-7 (四極の積分器) だけ。(3) B226: 四極の本体の積分器 (tsolqu・tsolque の Strang 分割の 2 次の誤り (aln²·k)·p²) は D2 では直していない (**v2.4.1 で 4 次の合成にして直した**。上の v2.4.1 の節)。厳密にするには刻み 14 倍か 4 次の合成が要り、11 か所の重複した実装の共通化なしには性能を落とす (docs/21 B226)。

**結果が変わるもの (追加)**: K1 ≠ 0 の MULT の端 (K0 の大きい極で 1e-5 以下)。K1 < 0 の MULT の AE のある Emittance[]。FRINGE = −1・−2 のステアリングの追跡 (AE の有無によらず)。

**要判断・未解決 (D2)**: (a) 四極の本体の積分器 (B226) を v2 で厳密にするか → v2.4.1 で 4 次の合成を共通化して入れた (性能と DA の結果が動く)。(b) MULT の端の縦収束は薄い蹴りの縦成分で、厳密な Maxwell (ANGLE=0 の BEND の鎖) と 2e-5 違う (K0 = 0.06・AE = 0.3)。(c) 追跡の放射の刻み (K64 の `ndivrad` の既定の EPS) が粗く、追跡の Jx が Emittance[] の収束値より 0.27% 小さい (B227。K1 < 0 の MULT に限らず BEND の極でも同じ。原因は特定済み、v2 でも未修正)。(d) ユーザー指示の「CALC の δ 列に K0 を足す」は B222 のとおり行っていない (足すと追跡・6D 行列と不整合)。RADINT の極の補正 (B115) は残す。

### v2.4.0 D3 (スピン・放射・SOL の残り: B230〜B242)

記録 `docs/history/v24_D3.md`、バグ表は docs/21 の B230〜B241、試験 `pysad/tests/test_v24_D3_spin_bpar_magnets.py` (17 件)・`test_v24_D3_rad_rf_sol.py` (12 件)、独立な基準 `_spin_maxwell_indep2.py` (直線の磁石の Maxwell 場の Thomas-BMT。ずれ・ROTATE つき)・`_rad_beta_indep.py` (Ring4B にエネルギー kick の空洞と modal エミッタンスを追加)。k64 は変えていない (k64 の試験・回帰は動かない)。

| 項目 | k64 | v2 | 根拠 |
|---|---|---|---|
| B230・B231: QUAD・MULT の本体の運動量の向きの磁場 (B165 の続き。ずれ・ROTATE・K1 < 0・F1/FRINGE つきも確認) | 端の ψ だけ (本体なし)。ndiv = 1 の四極は入口の +ψ を出口の −ψ が上書き | 刻みごとに Δψ/(1+δ) を足す (追跡 `p.bpadd`・MULT は全次数の ψ、Emittance[] `st.bpot`)。ndiv = 1 は足し算 | 直線の磁石の Maxwell 場の Thomas-BMT の積分: z 軸まわりの回転 R01 の誤差 1e-5 → 3e-9 (QUAD)、6.5e-6 → 3e-9 (MULT)。SEXT・ステアリングは元から 2e-8 以内 (端も本体も入れないので和が 0) |
| B232: RF の同期位相 (エネルギー損失・電荷) | Δp c を電圧と比べる (β = 1 の前提)、電荷を割らない | β·Δp c/\|z\| (`emit6d._u0_energy`。γ < 300 の β の係数が有効なとき) | 同期位相の追跡 1 周の δ の残り/(U0/pc): 3 GeV/c の陽子 (z = 1000) で 1e-6 以下 (β を入れないと 0.05、電荷を入れないと 1e3) |
| B233: PhotonList の光子のエネルギー | dp m βγ (β = 1 の電子では γ と同じ) | 陽子は dp·β p0 m (エネルギー保存)。低 γ の軟らかい光子の方向は単位球に載せる | 光子の総和/粒子の損失 = 1.000 (6 GeV/c)・0.998 (3 GeV/c)。以前は 1/β² 倍 |
| B167: 独立な (x, px, δ, z) 写像との εx・σδ | — | 変更なし (確認のみ) | 3 GeV/c の陽子 (z = 1000)、24 セルの扇形 FODO、RF あり: Emittance[] と独立な写像 (z_I = 1.7e5、刻み 0.0125 m) が εx 4e-5、σδ 8e-6 で一致 (刻み 0.02 m で 1e-3・3e-4)。1.5 GeV/c (γ = 1.89、β = 0.848) も εx 6e-5・σδ 4e-4。r2/β⁴ の因子の独立な確認 |
| B234: SOL 領域内の四極の本体の B∥ (B146 の ψ の確定) | 螺旋の bz·al + DISFRIN=0 の端の ψ | 本体の Δψ/(1+δ) を足す (追跡 `p.bpar`、Emittance[] `solemit._quad_trk`) | 一様なソレノイド + Maxwell の四極 (DISFRIN=0) と硬い端 (DISFRIN=1) の Thomas-BMT: 以前 5.5e-6 → 追跡 1e-6・Emittance[] 5e-7 以内。**B146 の DISFRIN=0 の ψ は物理的に正しい** |
| B235: 追跡の BEND のスピンの刻み | 放射の刻み (0.5 GeV の 15° 扇形で 5° の 3 刻み) | POL ON の追跡は 1 刻みの偏向角 5e-3 rad 以下 (`_ndivrad`) | dR/dpy の R12・R21: 1.2e-2 → 6e-4 (刻みに比例する 1 次誤差だった。端の ψ の分解の近似ではなかった) |
| B236: F1 の fb1 ≠ fb2 (B145) | — | 変更なし (誤りではない) | 厳密な運動との差は δ に依らない定数 (K0/al)(fb1² − fb2²)/24 = 設計軌道の変位。px・δ 依存は 1e-10。基準軌道の取り方 |
| B237: `Bend.transfer_matrix` の K0≠0・ANGLE≠0 | 旧式 (h = ANGLE/L) | CALC の `_bend_element_cod_matrix` (追跡の COD = 0 での微分) | 追跡のヤコビアンと 1e-12 (K1 = 0)。K1 ≠ 0 は 8e-6 (追跡の数値積分の刻み) |
| B238: SOL 内の CAVI の TRPT の加速 + POL | SADError | `sol_mult_full` に TRPT を渡し、新しい基準運動量を全粒子に、後続の bz を p0_ref/p0 倍 | E_z + B_z の運動方程式・Thomas-BMT (E×β を含む): γ が 15% 増える加速で 8e-6 (V = 0 で 4e-8)。残りの原因は、加速による向きの変化を磁場の曲げと同じ (1+aγ) で回す K64 の扱い (E 場では aγ + γ/(γ+1)、差 θ/(γ+1)): 基準の E×β の係数をその比倍にすると 57〜78% 減る |
| B239: 短い磁石・硬い端の放射 (低 γ) | — | 変更なし (限界) | 基準が取れない。誤差の目安は docs/46 |
| B240: Emittance[] の 6 極 (L ≠ 0) の放射・スピン | tthine は放射しない (spin・D・B に入らない) | 追跡と同じ 1 刻みの tradke (`_tthine(enarad=True)`、numba の k_thin は L ≠ 0 で使わない) | Maxwell 場の Thomas-BMT で spin 回転行列 3.5e-4 → 2e-6。U0 は追跡と 4e-6 |
| B242: 端の ψ の運動量 (1+δ) の規格化 | bsi・bsir0 は p0 規格化のまま (本体 Δψ/(1+δ) との和が δΔψ 残る) | (1+δ) で割る (`_sprot`・`tradke_spin`) | δ = 0.02・−0.03 の QUAD で R01 が Maxwell 場の積分と 1.2e-7 → 3e-9 (追跡・Emittance[]) |
| B241: 追跡の MULT の ROTATE のスピン | 座標だけ回す (pysad の欠け) | 入口で回し出口で戻す | MULT (ROTATE = 0.2) のスピン回転行列が基準と 5.7e-4 → 4.5e-7 |

**結果が変わるもの (v2 のみ)**: POL ON の追跡 (QUAD・MULT のスピン、BEND の刻み)、POL の Emittance[] の QUAD・MULT・SOL 内の四極・6 極 (L ≠ 0。放射も: 軌道が 6 極の中心を通らないとき U0・減衰・εx が動く)、陽子・イオン (γ < 300) の RF の同期位相 (trf0・RADCOD OFF の閉軌道の z)・PhotonList のエネルギー、`Bend.transfer_matrix` (K0≠0・ANGLE≠0)。電子の放射・RF は変わらない。参照値・許容の更新 (QUAD の本体の B∥・追跡の刻み・MULT の ROTATE で動いた POL の値): `test_v24_pol_emitdiv` の r1 (2.34228e-8 → 2.40756e-8)、`test_v2_B3_pol`・`test_v2_P1_pol_chrom` の線形の脱偏極率 (→ 2.40756e-8)、`test_v2_A_spin` の P (0.8695 → 0.8663)、`test_v2_B3_radpol` の丸め (1e-14 → 2e-13)、`test_v2_B1_sol` の POL + TRPT (SADError の試験 → 高速追跡との一致)、`test_v2_P1_st` の P_eq の許容 (2e-3 → 8e-3。docs/46 の未解決)。

**未解決 (D3)**: docs/46 の「v2.4.0 D3」。

## v2.4.1 LMA の誤差つき・長時間追跡の gs との比較 (確認のみ。コード変更なし)

**誤差つき LMA と 10⁴ 周級の長時間追跡 (放射・FLUC) の gs との比較 (確認。pysad のコードは変えていない。互換モードの k64 は gs と一致、v2 は K64 の誤りの修正で違う)**:
- **誤差つき LMA** (ASTRID2、0.58 GeV、周長 45.7 m、RF 50 kV、RFSW・RAD ON、FLUC なし、アパーチャ ±10 mm、LOSSDZ = バケット、2624 周 = シンクロトロン振動 4 周期): 位置ずれ・ロール (四極 50 μm・0.2 mrad、6 極 50 μm、偏向磁石 50 μm・0.1 mrad の rms) と K の相対誤差 (四極 K1・6 極 K2 とも 1e-3 rms) の 4 つの種を pysad で作り、**誤差の値を要素定義として台本へ書いて gs と pysad に同じ台本を渡した** (乱数列は同じにできないため)。位置 6 か所 (2・17・40・75・110・150)、台本の二分法 (9 回、[0, 0.02]、幅 3.9e-5) で δ± の 48 区間を比べた。**gs と pysad k64 は 48 / 48 が完全に同じ区間**。pysad v2 は 38 / 48 が同じ、6 が 1 刻み違い、4 が 2 刻み以上 (最大 6.6e-4 = 0.07%、アパーチャが限界の位置の δ− と 150 番の δ+)。δ+ の平均・標準偏差 (gs 1.435% ± 0.304%、v2 1.435% ± 0.303%)、δ− (1.446% ± 0.268%、1.450% ± 0.263%) は同じ。v2 の差は四極の 4 次の積分器 (B226・B243) と放射の修正による軌道の違いで、LMA の統計には出ない。1312 周 (2 周期) と 2624 周 (4 周期) の比較 (種 1・2): 位置 17 (δ−)・150 (δ±) の限界が 1e-4〜6e-4 (絶対) 動く (誤差なしの長直線部・バケットで決まる位置は動かない)。**誤差つきは 2000 周以上** の既知事項と整合。二分法は損失の窓 (狭い領域) を見落とすので、粗い走査つきの `momentum_acceptance` (check_every = 0) は種 1・位置 150 の δ+ で二分法より小さい 1.135% を返した (二分法 1.164%、gs も同じ 1.164%)。これは gs と pysad の差でなくアルゴリズムの差。
- **10⁴ 周の全周追跡 (放射 RAD・RF ON)**: KEK-LS 全周 (3 GeV、603.12 m、3563 要素、`kek-ls_v03_12.lat`) と ASTRID2 全周で、NOFLUC (10 粒子、初期振幅 x 0.5 mm・y 0.3 mm・δ 3e-3 [ASTRID2 は 1 mm・0.5 mm・3e-3]) と FLUC (40 粒子、振幅 0 の閉軌道から)。250 周ごとに記録。**gs と pysad k64 の座標は粒子ごとに一致**: NOFLUC 6e-14 (KEK-LS)・7e-14 (ASTRID2)、**FLUC でも 1e-14 (KEK-LS)・1.7e-14 (ASTRID2)、どちらも 40 粒子 10⁴ 周。K64 の乱数列まで pysad が再現している**。したがって統計誤差を気にせず、K64 と pysad k64 の同一性は粒子ごとに確認できた (統計つきの比較が要るのは v2 だけ)。
- **v2 は gs と違う** (K64 の修正の結果): KEK-LS NOFLUC の座標の差は 250 周で 3.8e-6、10⁴ 周で 2.7e-4 (ASTRID2 は 1.7e-5、1.7e-4)。FLUC (40 粒子) は粒子ごとには一致せず、座標の差が 250 周で 2.3e-3 (KEK-LS)・7.7e-3 (ASTRID2)、10⁴ 周で 1.4e-2・3.0e-2 (統計で比べる。乱数の使い方と放射の修正のどちらの寄与かは未分離)。原因の切り分け (KEK-LS 250 周): RAD OFF では v2 − k64 = 2.1e-7 で、`quadbody.set_scheme(0)` (Strang) にすると 3e-16 (**四極の 4 次の積分器 B226・B243 だけ**)。RAD ON では Strang にしても 3.7e-6 残る (**v2 の放射の修正**。個別の項目への切り分けは未実施、B76 の rhoe の更新が候補)。
- **減衰率 (NOFLUC、10 粒子の J の対数の傾き。Emittance[] の α [振幅/周] との比)**: 追跡の x: gs (= k64) は KEK-LS 0.945 ± 0.006・ASTRID2 0.762 ± 0.05、v2 は 0.991 ± 0.007・1.042 ± 0.05。y: 両方 0.97〜1.00。z (KEK-LS の包絡): 0.999。Emittance[] は放射積分の和則 Jx + Jy + Jz = 4 を満たす (ASTRID2 Jx = 1.375、KEK-LS 1.57) ので、**K64 の追跡の x の減衰が放射積分から計算した値より遅い** (K1≠0 の BEND の放射の追跡側の誤り B76 (v2 で修正済) の大きさの確認。ASTRID2 の複合機能偏向磁石 (K1L = −1.68、30°) で 24%、KEK-LS で 5.5%)。**新しい K64 のバグではなく、B76 の実リングでの大きさの確認**で、B76 に原因を切り分けたわけではない (K64 の追跡の放射だけを個別に切り替える手段がない)。ASTRID2 の x の傾きの誤差は 6% (10⁴ 周が τx = 1.4×10⁵ 周の 7% のため)。
- **FLUC の統計 (KEK-LS・ASTRID2)**: Emittance[] の平衡値 (εx = 1.549e-10 m、σδ = 6.947e-4、σz = 3.77 mm) と減衰率から作る予想 ε(n) = ε_eq (1 − e^(−2αn)) に対する比 (t ≥ 1500 周の平均、標準誤差は 1500 周おきの 6 点の見積もりで粗い): gs・k64 N = 40: σδ² 0.86 ± 0.04、εx 0.90 ± 0.06、σz² 0.96。v2 N = 40: 1.18 ± 0.07、0.92 ± 0.07、1.08 (これだけでは 0.86 と 1.18 は N = 40 の標本の揺らぎ。1 点ごとの統計は σδ² で 22%、εx で 16%)。**k64 N = 120: σδ² 0.99 ± 0.03、εx 0.95 ± 0.03、v2 N = 120: 1.04 ± 0.05、1.02 ± 0.02、σz² 0.90** (N = 40 の揺らぎが N = 120 で 1 に近づいた)。垂直 εy は gs・pysad とも 1e-39 (結合なし、光子の角度による励起は入れていない)。10⁴ 周は KEK-LS で τx の 0.77 倍 (εx の予想の 79%) で、**平衡ではない**。ASTRID2 は τx = 1.4×10⁵ 周 (10⁴ 周で 14%) で平衡への近さの試験にならない。**ASTRID2 N = 40 (予想は gs の Emittance[] の値 εx 1.268e-8 m、σδ 4.33e-4、σz 22 mm): gs = k64 は σδ² 0.95 ± 0.10、εx 1.20 ± 0.04、σz² 0.95。v2 は 0.80 ± 0.05、0.94 ± 0.02、0.83**。上の ± は 6 点の見積もりだが、6 点は同じ 40 粒子の追跡で互いに独立でなく過小。N = 40 の 1 点の統計は σδ² 22%・εx 16% なので、gs の εx 1.20 も v2 の σδ² 0.80 も 1 と区別できず、gs と v2 の違い (εx 1.20 対 0.94) も統計の範囲 (K64 の追跡の x の減衰が遅い (B76 候補) ことの効果は、この N では見えない)。
- 台本・参照出力: `pysad/tests/scripts/v241_lma_err_s{1..4}_t2624.{sad,gs.txt,pysad.txt,pysadk64.txt}`、`v241_lma_err_s{1,2}_t1312.*`、`v241_lma_long_{kekls,astrid2}_{nofluc,fluc}_n*_t10000.*`、`v241_lma_long_summary.out.txt`、`v241_lma_err_summary.out.txt`。生成・解析は `pysad/tests/_v241_lma_*.py`。試験 `pysad/tests/test_v241_lma_long.py` (16 件。15 件が既定で通り、live の二分法 1 件は `PYSAD_SLOW=1`)。
- **未確認・限界**: (1) 誤差つき LMA は ASTRID2 の 4 種・6 位置・2624 周の二分法だけ (KEK-LS の誤差つきの gs 比較・補正つき・FLUC つきは未実施。位置の合計 48 区間で gs = k64 なので追加の標本の価値は小さい)。gs の乱数は使っていない (誤差の値を共通にした)。(2) 長時間追跡は 10 粒子 (NOFLUC)・40 粒子 (FLUC) の 10⁴ 周で、**平衡エミッタンスには届いていない** (KEK-LS 0.77 τx、ASTRID2 0.07 τx)。平衡の比較は立ち上がりの傾き (予想との比) まで。(3) v2 の放射の差の項目別の切り分け (B76・B227 など) は未実施。(4) 1 セルの周期的な追跡・衝突・集団効果・RAD なしの長時間は見ていない。(5) 損失の起こる長時間 (LMA の限界付近の 10⁴ 周) は見ていない。(6) 計算はこの PC の 1 台 (gs: WSL、KEK-LS FLUC 40 粒子 10⁴ 周が 2912 s、NOFLUC 10 粒子 47 s)。


## v2.4.1 F (docs/46 の未実装項目の整理: FOURIER・光子と COORD・TCAVI・FIFO・SOL BOUND=0・qwspac)

条件: 互換モードによらず同じ動作 (追加は K64 の動作の再現で、v2 の物理の変更は無い)。K64 のバグは docs/21 の B320〜B322。

| 項目 | 結果 | 基準 |
| --- | --- | --- |
| MAIN の FOURIER 解析 (wjfit・wfres・wfbin・wiord) | **実装** (`pysad/fourier_fit.py`、`maintrack.run_track` から。`ON SMEAR FOURIER;`、MAXORDER・ADDTERMS・ADDDENSE)。`wjdraw` (描画) は実装しない | gs は MAIN の TRACK が SIGSEGV。K64 の Fortran 自体を gfortran で組んで比較: 7 条件で表示は全行一致、係数は相対 1e-7、チューンは 1e-12 |
| 光子 (PHOTONS) と COORD・TCAVI | **対応** (track6 に COORD・TCAVI を通す `_step_coord_tcav_ens`。ffscmd の光子の条件から COORD・TCAVI の除外を外した)。MAP は実装しない (理由は docs/46) | gs と k64: 光子 278 個が個数・12 列 (gs の出力桁) で一致 |
| TRCOD の FIFO (PUTMEASURE・GETMEASURE) | **確認** (コード変更なし)。WSL の mkfifo の FIFO で書き手・読み手が待ち合わせて同じ値を受け渡す。試験は Windows で skip | nfgetm.c の書式。gs の比較は無い (B151) |
| SOL BOUND=0 で BZ≠0 の領域外 | **追跡だけ対応** (K64 の tsol は BOUND=0 でも次の BOUND≠0 まで領域として追う: B321)。Twiss・CALC・幾何は未対応のまま (gs が SIGSEGV) | gs と追跡 18 座標が全桁一致 (両モード) |
| qwspac (4 次元 CALC の空間電荷) | **実装しない** (理由は docs/46・B320: gs は εy = 0 で NaN、結合があっても 4D CALC は Unstable) | gs の NaN・Unstable の原因調査 |

- 変更したコード: 新規 `pysad/fourier_fit.py`。`pysad/maintrack.py` (FOURIER ON のとき `fourier_report` を呼ぶ 8 行)、`pysad/track6.py` (`_step_coord_tcav_ens` の追加、BOUND=0 の SOL を領域外で拒否する例外の削除)、`pysad/ffscmd.py` (光子の条件 1 か所)。numba の kernel は編集していない。k64 のビット同一: 回帰 QUICK=1 の 50 ベンチ FAIL 0 (下記の確認)。
- 試験: `test_fourier_fit.py` (6 件)、`test_v2_F2_photons_coord_tcavi.py` (2 件)、`test_v2_F3_measure_fifo.py` (POSIX だけ)、`test_v2_F4_sol_bound0.py` (2 件)。例: `example/v2_next/F2_photons_coord_tcavi.sad`・`F4_sol_bound0_track.sad` (gs の出力 `*.gs.txt` つき)。
- 未確認: FOURIER を実際の追跡データで Fortran と比べていない (合成データ)、TCAVI の RANKICK・RANPHASE と光子の併用、COORD の幾何の上の光子の位置 (gs と一致した条件は DX DY DZ CHI1 CHI2 CHI3 の対)、FIFO の ORBITJITTER。

## v2.4.1 (偏光の FLUC つき追跡と gs の比較、B270)

記録 `docs/history/v241_pol_fluc.md`、試験 `pysad/tests/test_v241_pol_fluc.py`、台本・参照出力 `pysad/tests/scripts/v241_pol_*`。

| 項目 | 結果 |
|---|---|
| 同一性 | pysad k64 の FLUC つきスピン追跡は 12 セル FODO (6・24 MV)・PF・ASTRID2 で gs と粒子ごとのスピンまで全桁一致 (相対差 0〜7e-9) |
| 脱偏極率 (式 / gs の追跡) | 24 MV 1.05 ± 0.04、PF 0.98 ± 0.10、6 MV 0.49 ± 0.12 (追跡の平均が少数の粒子に支配される重い裾。上位 1% を除くと 0.82)。ASTRID2 は 40 粒子 × 18000 周で率 (1.04 ± 0.68)e-10 /周 (D。A・L は負) が式 5.4e-14 の 1σ の感度 (700〜1250 倍) に届かず評価不能 |
| スピンチューン | k64 の Emittance[] は gs の追跡と 3.8e-4 (FODO)・9.2e-4 (PF)・4.8e-3 (ASTRID2) 違う。主因は K64 の B58 (出口の端の縦磁場を捨てる)。v2 は 1e-5 台 |
| **B270 (pysad 自身の誤り、v2.4.1 で修正)** | 台本の `Emittance[]` (`emit_extra._turn_record`) が L ≠ 0 の 6 極に B240 を渡さなかった。6 極つきリングで軌道がずれると POL の量 (SpinTune・PolarizationVector・Polarization・PolarizationTime) が `emittance6d` と違った (ASTRID2 3.6e-3、PF 1e-3)。k64 は不変 |

**結果が変わるもの (v2 のみ)**: 6 極 (L ≠ 0) のあるリングの POL ON の台本の `Emittance[]` (軌道が 6 極の中心からずれるとき)。
**未解決**: ASTRID2 の追跡の統計、6 MV の裾 (数万粒子が要る)、v2 の物理の追跡の率 (pysad の Python 追跡が遅い)、ASTRID2 の Emittance[] の SpinTune が gs と pysad k64 で 4e-7 違う (閉軌道の収束の差と推定)。
