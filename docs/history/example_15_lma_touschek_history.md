# example/15_lma_touschek の整理の履歴

`example/v2_3`・`v2_4` の作業フォルダから機能別のフォルダへ移したときに、md の本文・台本の冒頭コメントから取り除いた作業符号・経過の記述を、元の文脈のまま残す。

旧パス→新パスの対応 (主なもの): v2_3 の `Q1_fmap_*_v241`→`fmap_*`、`Q1_lma_*`→`lma_*`、`Q1_gs_*`→`gs_*`、`S1_*`→`fmap_*`、`S2_*`→`fmap_check_*` (`S2_data/`→`fmap_check_data/`)、`S3_*`→`fmap_predict*`、`R5_fmap_rdt2`→`fmap_rdt4`、`L2_*`→`lma_*` (`L2_elegant_isr`→`touschek_elegant/elegant_isr`)、`T1_*`→接頭辞なし、`T2_*`→`touschek_formula_*`、`U2_*`→`bend_k1_def_*`、v2_4 の `A3_*`→`touschek_bunch_*`、`A4_*`→`lma_support*`、`H2_*`→`lma_corrected_*`。旧い 2 次の積分器で作った `Q1_fmap_{astrid2,kekls}.{out.txt,table.txt,png}` (6 ファイル) は削除し、4 次の合成の `_v241` 版を `fmap_{astrid2,kekls}.*` にした。

## md から取り除いた文

### lma_reliability.md

元: 記録は `docs/history/v23_L2.md`、解説は

後: 解説は

### lma_reliability.md

元: | `lma_common.py` | ラティスの読み込み、T1 の elegant の表の読み込み | — |

後: | `lma_common.py` | ラティスの読み込み、`../touschek_elegant/` の elegant の表の読み込み | — |

### lma_reliability.md

元: ## 結果の要約 (数値は `docs/history/v23_L2.md`)

後: ## 結果の要約

### lma_reliability.md

元: - **`check_every` の既定を 1 (要素ごとの判定) に変えた**。SAD の tapert

後: - **`check_every` の既定は 1 (要素ごとの判定)**。SAD の tapert

### lma_reliability.md

元: はすべてアパーチャで決まる (T1 が正しい)。

後: はすべてアパーチャで決まる (`../touschek_elegant/` の elegant の LMA と一致)。

### bend_k1_definition.md

元: # K1≠0 の BEND の K1 の定義のオプションと B161 の例

後: # K1≠0 の BEND の K1 の定義のオプションと CALC の dk の修正の例

### bend_k1_definition.md

元: 記録は `docs/history/v23_U2.md`、試験は

後: 試験は

### bend_k1_definition.md

元: (3) B161 の影響 (n = 0.4 のリング

後: (3) CALC の dk の修正の影響 (n = 0.4 のリング

### bend_k1_definition.md

元: T1 の `astrid2_all.lte` の K2

後: `sad2ele.py` の `astrid2_all.lte` の K2

### bend_k1_definition.md

元: elegant との比較は T1 の保存値 (`bend_k1_elegant.out.txt`

後: elegant との比較は保存値 (`bend_k1_elegant.out.txt`

### bend_k1_definition.md

元: ## 結果の要約 (条件・限界は `docs/history/v23_U2.md`)

後: ## 結果の要約

### bend_k1_definition.md

元: 相対 1.6e-3)。B160 に鈍感。

後: 相対 1.6e-3)。K1 の定義 (`log`・`lin`) に鈍感。

### bend_k1_definition.md

元: B161 (CALC の dk): v2 は

後: CALC の dk (`compat.set_bend_dk_fix`): v2 は

### elegant_comparison.md

元: 記録は `docs/history/v23_T1.md`、試験は `pysad/tests/test_v23_T1_elegant.py`。この例は

後: 試験は `pysad/tests/test_v23_T1_elegant.py`。この例は

### elegant_comparison.md

元: sad2ele.py astrid2 example/v2_3
> 

後: sad2ele.py astrid2 example/15_lma_touschek/touschek_elegant


### elegant_comparison.md

元: ## 結果の要約 (条件は `docs/history/v23_T1.md` §2)

後: ## 結果の要約

### elegant_comparison.md

元: - LMA (kek-ls): `docs/history/v23_T1.md` §4.2。
> 

後: - LMA (kek-ls): 14 位置 (1500 周) の elegant と、同じ位置の pysad の LMA を `lma_compare_kekls.out.txt` に載せた (周回数・RF バケットの扱いは `../lma/lma_reliability.md`)。


### elegant_comparison.md

元: 。B160 の根拠 |

後: 。B160 (K1 の暗黙の定義) の根拠 |

### elegant_comparison.md

元: である (docs/11 A59、`docs/history/v23_U1.md`)。elegant には

後: である (docs/11 A59)。elegant には

### touschek_formula.md

元: 記録は `docs/history/v23_T2.md`、試験は `pysad/tests/test_v23_T2_touschek.py` (8 件)。

後: 試験は `pysad/tests/test_v23_T2_touschek.py` (8 件)。

### touschek_formula.md

元: T1 の lte (`<名前>_all.lte`) を使う。

後: `sad2ele.py` の lte (`<名前>_all.lte`、`elegant_comparison.md`) を使う。

### touschek_formula.md

元: ## 結果の要約 (条件は `docs/history/v23_T2.md`)

後: ## 結果の要約

## 台本の冒頭コメントから取り除いた・書き直した行

### lma/gs_loss_astrid2.sad

元: `LOSSAMPL = 0.01; LOSSDZ = 1.5406914179745479;`

後: `LOSSAMPL = 0.01; LOSSDZ = 1.5406914179745532;`

元: `beam = {{-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11,-3.014645389754925e-11},{-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.2120939087369584e-07,-2.21209`

後: `beam = {{-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11,-3.0146460663953359e-11},{-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.2120939087526666e-07,-2.21209390875266`

元: `beam = {{-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07,-5.0856794153997456e-07},{2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.1481157368861496e-07,2.14811`

後: `beam = {{-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07,-5.0856794153663778e-07},{2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.1481157369118038e-07,2.14811`

元: `cod = {-3.014645389754925e-11,-2.2120939087369584e-07,8.5362930388775991e-35,-4.682044508848058e-35,-2.039658521037257e-07,-1.0379519492506624e-11}; lo = 0; hi = 0.02;`

後: `cod = {-3.0146460663953359e-11,-2.2120939087526666e-07,1.066290925630897e-34,-5.9571858257343504e-35,-2.0396585261861938e-07,-1.0379517242241593e-11}; lo = 0; hi = 0.02;`

元: `cod = {-3.014645389754925e-11,-2.2120939087369584e-07,8.5362930388775991e-35,-4.682044508848058e-35,-2.039658521037257e-07,-1.0379519492506624e-11}; lo = 0; hi = 0.02;`

後: `cod = {-3.0146460663953359e-11,-2.2120939087526666e-07,1.066290925630897e-34,-5.9571858257343504e-35,-2.0396585261861938e-07,-1.0379517242241593e-11}; lo = 0; hi = 0.02;`

元: `cod = {-5.0856794153997456e-07,2.1481157368861496e-07,-4.7339140179979798e-35,-4.4722591934182531e-35,-5.298653076911784e-08,-8.8603150216192072e-07}; lo = 0; hi = 0.02;`

後: `cod = {-5.0856794153663778e-07,2.1481157369118038e-07,-6.153834394270825e-35,-5.6176506379420251e-35,-5.2986531282534043e-08,-8.8603150215967089e-07}; lo = 0; hi = 0.02;`

元: `cod = {-5.0856794153997456e-07,2.1481157368861496e-07,-4.7339140179979798e-35,-4.4722591934182531e-35,-5.298653076911784e-08,-8.8603150216192072e-07}; lo = 0; hi = 0.02;`

後: `cod = {-5.0856794153663778e-07,2.1481157369118038e-07,-6.153834394270825e-35,-5.6176506379420251e-35,-5.2986531282534043e-08,-8.8603150215967089e-07}; lo = 0; hi = 0.02;`

### lma/gs_loss_compare.py

元: `"""v2.3.0 Q1: pysad の追跡 (trksetup.TrackSetup.track = LMA・周波数マップが使う追跡) と gs (SAD K64) の TrackParticles が、同じ損失判定になることの確認。`

後: `"""pysad の追跡 (trksetup.TrackSetup.track = LMA・周波数マップが使う追跡) と gs (SAD K64) の TrackParticles が、同じ損失判定になることの確認。`

### lma/lma_cell.py

元: `"""v2.3.0 L2 項目 5: 1 セルの結果を全周に延ばす影響。RF 空洞は 1 か所なので、全セルの同じ位置の LMA は同じにならない (周期性が崩れる)。`

後: `"""1 セルの結果を全周に延ばす影響。RF 空洞は 1 か所なので、全セルの同じ位置の LMA は同じにならない (周期性が崩れる)。`

### lma/lma_check_every.py

元: `"""v2.3.0 L2 項目 1: check_every (損失の判定の間隔) の影響と計算時間。elegant (要素ごとの MAXAMP = check_every 1 と同じ) と比べる。`

後: `"""check_every (損失の判定の間隔) の影響と計算時間。elegant (要素ごとの MAXAMP = check_every 1 と同じ) と比べる。`

元: `        locs = [l for l in locs if l in ele]          # T1 の elegant と同じ位置 (最初のセルの 11 点)`

後: `        locs = [l for l in locs if l in ele]          # elegant_lma.py の elegant と同じ位置 (最初のセルの 11 点)`

### lma/lma_common.py

元: `"""v2.3.0 L2 の共通部分: ラティスの読み込みと位置の選び方 (ASTRID2・KEK-LS・MPW 版 ASTRID2)。リポジトリ直下から実行する。"""`

後: `"""LMA の例の共通部分: ラティスの読み込みと位置の選び方 (ASTRID2・KEK-LS・MPW 版 ASTRID2)。リポジトリ直下から実行する。"""`

元: `    """T1 の elegant (kek-ls、1500 周) の位置 (pysad の loc) と δ± の表: {loc: (s, dp+, dp-)}"""`

後: `    """elegant_lma.py の elegant (kek-ls、1500 周) の位置 (pysad の loc) と δ± の表: {loc: (s, dp+, dp-)}"""`

### lma/lma_corrected_kekls_v03_12.template.sad

元: `! v2.4.0 H2: KEK-LS v03_12 (kek-ls_v03_12.lat、3 GeV、603 m、3563 要素。docs/48 §4 の A4・L2 の「大きなリング A」) 版。lma_corrected_kekls_v03_68.template.sad のラティスを v03_12 に、`

後: `! KEK-LS v03_12 (kek-ls_v03_12.lat、3 GeV、603 m、3563 要素。「大きなリング A」) 版。lma_corrected_kekls_v03_68.template.sad のラティスを v03_12 に、`

### lma/lma_corrected_kekls_v03_68.template.sad

元: `! v2.4.0 H2: KEKLS_correct_cod.sad (誤差・COD の補正) の補正のあとに MomentumAcceptance[] を足したもの (pysad のみ。gs には MomentumAcceptance[] が無い)。`

後: `! KEKLS_correct_cod.sad (誤差・COD の補正) の補正のあとに MomentumAcceptance[] を足したもの (pysad のみ。gs には MomentumAcceptance[] が無い)。`

元: `! 変更点: seeds = {SEEDVAL}、niter = NITERVAL、ae・fe・te、ctrl = 1 のとき LMA の前に 'cell calc expand' (既定の CALEXP の対照)、補正のあとの DA を LMA に置き換え (Distance = Dist m ごと、`

後: `! 元の台本との違い: seeds = {SEEDVAL}、niter = NITERVAL、ae・fe・te、ctrl = 1 のとき LMA の前に 'cell calc expand' (既定の CALEXP の対照)、補正のあとの DA を LMA に置き換え (Distance = Dist m ごと、`

### lma/lma_corrected_run.py

元: `"""v2.4.0 H2: KEK-LS (v03_68、3 GeV、570.72 m、RING = 標準セル 19 + 入射セル) に誤差を与え、COD を補正し、補正後の状態で MomentumAcceptance[] を走らせる`

後: `"""KEK-LS (v03_68、3 GeV、570.72 m、RING = 標準セル 19 + 入射セル) に誤差を与え、COD を補正し、補正後の状態で MomentumAcceptance[] を走らせる`

### lma/lma_corrected_summary.py

元: `"""H2: v03_12 (RAD ON、982 周、9 位置の 8 つ = s = 0, 75.4, …, 527.8 m。最後の位置は最初と同じ) の lma_corrected_kekls_v03_12_rad_t982_p9_*.out.txt をまとめる。`

後: `"""v03_12 (RAD ON、982 周、9 位置の 8 つ = s = 0, 75.4, …, 527.8 m。最後の位置は最初と同じ) の lma_corrected_kekls_v03_12_rad_t982_p9_*.out.txt をまとめる。`

### lma/lma_errors.py

元: `"""v2.3.0 L2 項目 7: 誤差つきのリングと MPW つきのリングの LMA。
  errors: ASTRID2 に位置のずれ・ロール (S1 の fmapres.apply_alignment_errors: 四極 50 µm・0.2 mrad、6 極 50 µm、偏向 50 µm・0.1 mrad の rms ×level) を、`

後: `"""誤差つきのリングと MPW つきのリングの LMA。
  errors: ASTRID2 に位置のずれ・ロール (fmapres.apply_alignment_errors: 四極 50 µm・0.2 mrad、6 極 50 µm、偏向 50 µm・0.1 mrad の rms ×level) を、`

### lma/lma_fluc.py

元: `"""v2.3.0 L2 項目 6: FLUC (量子励起) つきの LMA。複数の乱数の種 (既定 6) で pysad の δ± の平均・標準偏差を求め、FLUC なしと、elegant の ISR つき`

後: `"""FLUC (量子励起) つきの LMA。複数の乱数の種 (既定 6) で pysad の δ± の平均・標準偏差を求め、FLUC なしと、elegant の ISR つき`

### lma/lma_limit_types.py

元: `"""v2.3.0 L2 項目 2: RF バケットの限界とアパーチャの限界の区別 (pysad/lmabucket.py)。`

後: `"""RF バケットの限界とアパーチャの限界の区別 (pysad/lmabucket.py)。`

### lma/lma_ring.py

元: `"""v2.3.0 Q1: 局所運動量アクセプタンス (LMA) の例 (pysad/lma.py)。1 セル (1 周期) の位置を計算し、全周の別のセルの数点で 1 セルの結果と合うことを確かめ、`

後: `"""局所運動量アクセプタンス (LMA) の例 (pysad/lma.py)。1 セル (1 周期) の位置を計算し、全周の別のセルの数点で 1 セルの結果と合うことを確かめ、`

### lma/lma_support.py

元: `"""v2.4.0 A4: LMA (局所運動量アクセプタンス) の FLUC・誤差・結合・k64 の対応範囲の確認 (L2 で未確認だった項目)。`

後: `"""LMA (局所運動量アクセプタンス) の FLUC・誤差・結合・k64 の対応範囲の確認。`

### lma/lma_turns.py

元: `"""v2.3.0 L2 項目 4: LMA の周回数への依存 (kek-ls)。T1 の elegant (1500 周、14 位置) と同じ位置で、周回数 491・1000・1500・2000 の δ± を求める。`

後: `"""LMA の周回数への依存 (kek-ls)。elegant_lma.py の elegant (1500 周、14 位置) と同じ位置で、周回数 491・1000・1500・2000 の δ± を求める。`

### lma/lma_turns_short.py

元: `"""v2.3.0 L2 項目 4 (続き): シンクロトロン振動 1 周期未満・数周期での収束。周回数をシンクロトロン周期の 0.25〜3 倍にして、バケットで決まる δ+ と、`

後: `"""シンクロトロン振動 1 周期未満・数周期での収束。周回数をシンクロトロン周期の 0.25〜3 倍にして、バケットで決まる δ+ と、`

### lma/lma_verify.py

元: `"""v2.3.0 L2 項目 3: 二分法が最初の損失しか見ないことの影響。'verify=n' (限界の手前 [0, δ] を n 点で再追跡) で、粗い走査の点のあいだの狭い損失の窓を数える。`

後: `"""二分法が最初の損失しか見ないことの影響。'verify=n' (限界の手前 [0, δ] を n 点で再追跡) で、粗い走査の点のあいだの狭い損失の窓を数える。`

### touschek_elegant/bend_k1_def_astrid2_chrom.py

元: `"""v2.3.0 U2: ASTRID2 (30 度の複合機能偏向磁石 K1 = -1.678 /m^2) の色収差 ξ・∂²ν/∂δ² を、pysad の K1 の定義 'log' (K64)・'lin' (B' 一様) と
elegant (CSBEND K2 = -h K1 と K2 = 0、EDGE_ORDER = 1・2) で比べる。B160 の「∂²ν/∂δ² の y の 3 倍の差が解消するか」の確認。`

後: `"""ASTRID2 (30 度の複合機能偏向磁石 K1 = -1.678 /m^2) の色収差 ξ・∂²ν/∂δ² を、pysad の K1 の定義 'log' (K64)・'lin' (B' 一様) と
elegant (CSBEND K2 = -h K1 と K2 = 0、EDGE_ORDER = 1・2) で比べる。「∂²ν/∂δ² の y の 3 倍の差が解消するか」の確認 (compat.set_bend_k1_def)。`

### touschek_elegant/bend_k1_def_lma_astrid2.py

元: `"""v2.3.0 U2: ASTRID2 の局所運動量アクセプタンス (1 セル、pysad の LMA) を K1 の定義 'log' (K64) と 'lin' (B' 一様) で比べる。`

後: `"""ASTRID2 の局所運動量アクセプタンス (1 セル、pysad の LMA) を K1 の定義 'log' (K64) と 'lin' (B' 一様) で比べる。`

### touschek_elegant/bend_k1_def_verify.py

元: `"""v2.3.0 U2: K1 の定義のオプション ('log'|'lin'、B160) と B161 (CALC の dk の修正) の検証の表を作る。すべて v2 モード。`

後: `"""K1 の定義のオプション ('log'|'lin'、compat.set_bend_k1_def) と CALC の dk の修正 (compat.set_bend_dk_fix) の検証の表を作る。すべて v2 モード。`

元: ` 1. n = 0.4 の扇形磁石 16 台のリング (3 GeV 相当): dν/dδ (CALC・追跡・独立な積分・解析式・elegant (T1 の保存値))。追跡のチューンと積分の差 (5 つの δ)`

後: ` 1. n = 0.4 の扇形磁石 16 台のリング (3 GeV 相当): dν/dδ (CALC・追跡・独立な積分・解析式・elegant (bend_k1_elegant.out.txt の保存値))。追跡のチューンと積分の差 (5 つの δ)`

元: ` 3. B161: 同じリングの CALC のチューンの δ 依存 (dk 修正の有無) と ASTRID2・大きなリング A の ξ・∂²ν/∂δ²・Emittance[] への影響`

後: ` 3. dk の修正: 同じリングの CALC のチューンの δ 依存 (dk 修正の有無) と ASTRID2・大きなリング A の ξ・∂²ν/∂δ²・Emittance[] への影響`

元: `import _u2_indep as U2  # noqa: E402`

後: `import _u2_indep as indep  # noqa: E402`

元: `        t = np.array([U2.tunes(rr, x, fields=fl, edge=True) for x in ds])`

後: `        t = np.array([indep.tunes(rr, x, fields=fl, edge=True) for x in ds])`

### touschek_elegant/bend_k1_elegant.py

元: `"""v2.3.0 T1: gs_bend_k1.py の一様な磁場指数のリング (n = 0.4、16 台の扇形磁石) を elegant (CSBEND) で計算し、dν/dδ を出す。`

後: `"""gs_bend_k1.py の一様な磁場指数のリング (n = 0.4、16 台の扇形磁石) を elegant (CSBEND) で計算し、dν/dδ を出す。`

### touschek_elegant/elegant_isr.py

元: `"""v2.3.0 L2 項目 6: elegant の momentum_aperture (ISR つき = 量子励起あり、SYNCH_RAD = 1・ISR = 1) を、位置ごと・種ごとに serial の elegant で走らせる (WSL、最大 4 本同時)。`

後: `"""elegant の momentum_aperture (ISR つき = 量子励起あり、SYNCH_RAD = 1・ISR = 1) を、位置ごと・種ごとに serial の elegant で走らせる (WSL、最大 4 本同時)。`

### touschek_elegant/elegant_lma.py

元: `"""v2.3.0 T1: elegant の momentum_aperture (追跡、RF・放射つき) で局所運動量アクセプタンスを求める (WSL の Pelegant)。`

後: `"""elegant の momentum_aperture (追跡、RF・放射つき) で局所運動量アクセプタンスを求める (WSL の Pelegant)。`

元: `    astrid2: 全周。kekls: 最初のセルの 0.5 m おきの位置 (pysad の Q1 と同じ select_locations) + 別のセル (3・9・15 番目) の同じ位置の 3 点ずつ"""`

後: `    astrid2: 全周。kekls: 最初のセルの 0.5 m おきの位置 (lma_ring.py と同じ select_locations) + 別のセル (3・9・15 番目) の同じ位置の 3 点ずつ"""`

### touschek_elegant/elegant_summary.py

元: `"""v2.3.0 T1: 保存済みの結果 (pysad_*.ref.json、elegant_*.touschek.txt、elegant_*.lma.txt) から、寿命と LMA の比較の表を作る (計算はしない)。`

後: `"""保存済みの結果 (pysad_*.ref.json、elegant_*.touschek.txt、elegant_*.lma.txt) から、寿命と LMA の比較の表を作る (計算はしない)。`

### touschek_elegant/elegant_touschek.py

元: `"""v2.3.0 T1: elegant の touschekLifetime (Piwinski の式。Xiao・Borland の積分) による Touschek 寿命。pysad と同じ入力にして比べる。`

後: `"""elegant の touschekLifetime (Piwinski の式。Xiao・Borland の積分) による Touschek 寿命。pysad と同じ入力にして比べる。`

### touschek_elegant/elegant_touschek_scan.py

元: `"""v2.3.0 T1: Touschek 寿命の式 (pysad = SAD K64 の TouschekTable と、elegant の touschekLifetime = Piwinski の式) の違いを切り分ける。`

後: `"""Touschek 寿命の式 (pysad = SAD K64 の TouschekTable と、elegant の touschekLifetime = Piwinski の式) の違いを切り分ける。`

### touschek_elegant/gs_bend_k1.py

元: `"""v2.3.0 T1: K1≠0 の BEND (複合機能) のオフモメンタムの光学を gs (SAD K64) と pysad で比べ、解析式・elegant と比べる。`

後: `"""K1≠0 の BEND (複合機能) のオフモメンタムの光学を gs (SAD K64) と pysad で比べ、解析式・elegant と比べる。`

### touschek_elegant/lma_compare.py

元: `"""v2.3.0 T1: pysad の LMA (lma.py) と elegant の momentum_aperture の局所運動量アクセプタンスを同じ位置で比べる。`

後: `"""pysad の LMA (lma.py) と elegant の momentum_aperture の局所運動量アクセプタンスを同じ位置で比べる。`

### touschek_elegant/pysad_ref.py

元: `"""v2.3.0 T1: pysad 側の基準値 (光学、IBS 込みのビーム、Touschek 寿命) を作る。elegant と同じ条件を揃えるための入力にもなる。`

後: `"""pysad 側の基準値 (光学、IBS 込みのビーム、Touschek 寿命) を作る。elegant と同じ条件を揃えるための入力にもなる。`

### touschek_elegant/sad2ele.py

元: `"""v2.3.0 T1: SAD のラティス (リング) を elegant の lte に変換する (pysad で展開した要素の並びから書く)。`

後: `"""SAD のラティス (リング) を elegant の lte に変換する (pysad で展開した要素の並びから書く)。`

元: `                # (一様な磁場指数のリングで 4 桁一致: docs/history/v23_T1.md)。sad_k1=False なら K2 を足さない (厳密な勾配)。`

後: `                # (一様な磁場指数のリングで 4 桁一致: gs_bend_k1.py・bend_k1_elegant.py)。sad_k1=False なら K2 を足さない (厳密な勾配)。`

### touschek_elegant/touschek_bunch_mc.py

元: `"""v2.4.0 A3: バンチ内の δ の分布 (touschek_bunch.py) の独立な検証: 実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学、touschek_formula_mc.py と同じ)。`

後: `"""バンチ内の δ の分布 (touschek_bunch.py) の独立な検証: 実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学、touschek_formula_mc.py と同じ)。`

### touschek_elegant/touschek_bunch_qs_3d.py

元: `"""v2.4.0 A10: 相対運動量の縦成分 q_s の効果を入れた 3 次元の積分 (touschek_bunch.qs_ratio) と、実験室系の対のモンテカルロ (touschek_bunch_mc.py の出力) の比較。`

後: `"""相対運動量の縦成分 q_s の効果を入れた 3 次元の積分 (touschek_bunch.qs_ratio) と、実験室系の対のモンテカルロ (touschek_bunch_mc.py の出力) の比較。`

### touschek_elegant/touschek_bunch_qs_ring.py

元: `"""v2.4.0 A10: ASTRID2 全周の Touschek の損失率 (一様な受容 a) の、δ_c = 0 (elegant・K64・Piwinski と同じ基準) に対する増え: δ_c だけ (1 次元) と δ_c + q_s (3 次元を 6 / 12 / 24 スライスで補間)。`

後: `"""ASTRID2 全周の Touschek の損失率 (一様な受容 a) の、δ_c = 0 (elegant・K64・Piwinski と同じ基準) に対する増え: δ_c だけ (1 次元) と δ_c + q_s (3 次元を 6 / 12 / 24 スライスで補間)。`

### touschek_elegant/touschek_formula_cause.py

元: `"""v2.3.0 T2: K64 の Touschek の式と独立な式の差の原因を、近似を 1 つずつ置き換えて切り分ける。`

後: `"""K64 の Touschek の式と独立な式の差の原因を、近似を 1 つずつ置き換えて切り分ける。`

### touschek_elegant/touschek_formula_check3d.py

元: `"""v2.3.0 T2: 独立な検証 (b-1): 3 次元の一般の積分。touschek_piwinski の 1 次元の式 (q_s = 0、x-y ブロックの 2 固有値) の近似を外して評価する。`

後: `"""独立な検証 (b-1): 3 次元の一般の積分。touschek_piwinski の 1 次元の式 (q_s = 0、x-y ブロックの 2 固有値) の近似を外して評価する。`

### touschek_elegant/touschek_formula_limits.py

元: `"""v2.3.0 T2: 解析式との比較 (極限) と、近似の大きさの見積もり。`

後: `"""解析式との比較 (極限) と、近似の大きさの見積もり。`

### touschek_elegant/touschek_formula_lma_lifetime.py

元: `"""v2.3.0 T2: 局所の運動量アクセプタンス (LMA) の Touschek 寿命を、K64 の表と Piwinski 型の式 (touschek_lifetime_local_piwinski) で出す。
LMA の表は Q1・T1 の保存済みの表 (ASTRID2: lma_astrid2.table.txt、KEK-LS: pysad_kekls.lma.txt。1 セルを周期的に延ばす)。IBS 込み PBUNCH = 1e9、MINCOUP = 0.01。`

後: `"""局所の運動量アクセプタンス (LMA) の Touschek 寿命を、K64 の表と Piwinski 型の式 (touschek_lifetime_local_piwinski) で出す。
LMA の表は lma_ring.py・lma_compare.py の保存済みの表 (ASTRID2: lma_astrid2.table.txt、KEK-LS: pysad_kekls.lma.txt。1 セルを周期的に延ばす)。IBS 込み PBUNCH = 1e9、MINCOUP = 0.01。`

### touschek_elegant/touschek_formula_mc.py

元: `"""v2.3.0 T2: 独立な検証 (b-2): 実験室系の 6 次元ガウス分布からの 2 粒子の対のモンテカルロ (厳密な相対論の運動学)。`

後: `"""独立な検証 (b-2): 実験室系の 6 次元ガウス分布からの 2 粒子の対のモンテカルロ (厳密な相対論の運動学)。`

### touschek_elegant/touschek_formula_scan.py

元: `"""v2.3.0 T2: Touschek 寿命の 3 通りの比較 (pysad = K64 の式 / elegant touschekLifetime / 独立な式 touschek_piwinski)。`

後: `"""Touschek 寿命の 3 通りの比較 (pysad = K64 の式 / elegant touschekLifetime / 独立な式 touschek_piwinski)。`

元: ` 入力 (T1 の lte '<名前>_all.lte' を使う): -emitInput = EMITXC、-emityInput = EMITYC、-deltaInput = SIGE、-length = SIGZ、-particles = PBUNCH。`

後: ` 入力 (sad2ele.py の lte '<名前>_all.lte' を使う): -emitInput = EMITXC、-emityInput = EMITYC、-deltaInput = SIGE、-length = SIGZ、-particles = PBUNCH。`

元: `    """全要素を 8 分割した elegant の twiss を ~/t2/<名前>/ に作る (T1 の lte を使う)"""`

後: `    """全要素を 8 分割した elegant の twiss を ~/t2/<名前>/ に作る (sad2ele.py の lte を使う)"""`

### touschek_elegant/touschek_formula_validate.py

元: `"""v2.3.0 T2: touschek_piwinski (1 次元の式) を、独立な 2 つの数値積分 (touschek_formula_check3d: 3 次元の一般の積分、touschek_formula_mc: 実験室系の対のモンテカルロ、厳密な相対論) と、スライスごとに比べる。`

後: `"""touschek_piwinski (1 次元の式) を、独立な 2 つの数値積分 (touschek_formula_check3d: 3 次元の一般の積分、touschek_formula_mc: 実験室系の対のモンテカルロ、厳密な相対論) と、スライスごとに比べる。`
