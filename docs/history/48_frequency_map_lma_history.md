# 48_frequency_map_lma の履歴 (整理で本文から取り除いた記述)

本文 (`docs/48_frequency_map_lma.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 3 行 (節: 48. 周波数マップと局所運動量アクセプタンス (LMA) の解説)

> 最終更新: 2026-11-20 (v2.3.0、S1 の検証を §9 に追加)

## 元 5 行 (節: 48. 周波数マップと局所運動量アクセプタンス (LMA) の解説)

> 関連: [docs/49 (共鳴駆動項・ADTS・クロマティシティ。v2.4.0 の J²・δ 依存は §13)](49_resonance_terms_adts.md)。**入門 (何をどう求めて何と比べたか) は [docs/50](50_adts_rdt_primer.md)**。

## 元 10 行 (節: 48. 周波数マップと局所運動量アクセプタンス (LMA) の解説)

> FFS からは使えない (docs/history/work_main_track.md)。チューンの取り出しは [`pysad/naff.py`](../pysad/naff.py) (NAFF、v2.1.0) を使う。

## 元 22 行 (節: 1. 何ができるか)

> | FLUC (量子励起) つきの LMA | 対応 (`TrackSetup.fluc = True`・`seed`。複数の種の平均・分散で使う。2 周期では効果は二分の幅以下、§4.5 (f))。周波数マップは S2 §10.3 | `lma.momentum_acceptance` |

## 元 23 行 (節: 1. 何ができるか)

> | 誤差つきのリング (位置のずれ・ロール、全周の追跡) | 周波数マップは走る (S1、ASTRID2、§9.4)。LMA も走る (L2、ASTRID2 level 1・3 種。v2.4.0 A4: level 5 = 結合の見積もり |c−| 0.107 が差の共鳴までの距離 0.04 を超える強い結合・閉軌道 4.8 mm でも走る。大きなリング A は `docs/history/v24_A.md` §A4)。補正つき (COD 補正後) は閉軌道を Newton 法で求めるので補正の状態をそのまま使えるが、台本から取り込む経路は未検証 | `fmapres.apply_alignment_errors` |

## 元 24 行 (節: 1. 何ができるか)

> | δ を含む周波数マップ (x–δ、y–δ) | 対応 (既存の `frequency_map` の dp の格子。S1 で ASTRID2・KEK-LS を確認、§9.3) | `freqmap.frequency_map` |

## 元 25 行 (節: 1. 何ができるか)

> | 共鳴駆動項 (RDT)・ADTS からの共鳴の領域の予想 | 対応 (1 次の RDT が駆動する共鳴だけ。S1、§3.5・§9.1) | `fmapres.ResonancePredictor` |

## 元 82 行 (節: 3.2 Python の API)

> - v2.3.0 S2 の追加 (§3.6、§10): `frequency_map(..., scaling=True, warn=True)`。`r.D_half` (最初の nturn/2 周の D)・`r.slope` (周回数への依存の傾き)、`r.D_norm` (1024 周に換算した D = D + log2(nturn/1024))、

## 元 98 行 (節: 3.3 SAD の台本から)

> 戻りは `{X->…, Y->…, DP->…, NuX->…, NuY->…, NuX2->…, NuY2->…, Diffusion->…, LostTurn->…}` (点ごとのリスト。失われた点の Diffusion は −99)。追跡の条件: **既定は RFSW OFF・RAD OFF・FLUC OFF** (v2.3.0 S2。現在のフラグは使わない。FFS に入ると RFSW は ON なので、`NORFSW;` を書き忘れると D > −5 の割合が 0.03 から 0.57 になるため。現在のフラグを使うときは `UseFlags->True`)、LOSSAMPL・LOSSDZ は現在の値。警告は ` FrequencyMap: WARNING [コード] …` で出力される (§3.6)。

## 元 125 行 (節: 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1))

> ### 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1)

## 元 138 行 (節: 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1))

> normal 6 極では調和は 3νx、νx、νx ± 2νy の 4 つだけ。**νx − νy などの 2 次の共鳴 (6 極の 2 次の項)、4 次以上の共鳴の強さの解析式はない**ので、予想できない (S1、§9.1)。

## 元 141 行 (節: 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1))

> **S3 の作り直し (v2.3.0、[`pysad/fmappred.py`](../pysad/fmappred.py)、`ResonancePredictor.from_lattice(..., model='s3')`、導出・限界は docs/49 §12)**: 駆動項に 4 次の h (6 極の 2 次ほか、`rdt2nd.rdt4`) と非共鳴の 3 次の補正 (`rdt4_normal_form` の F4')、

## 元 142 行 (節: 3.5 RDT・ADTS から予想する共鳴の領域 (`pysad/fmapres.py`、v2.3.0 S1))

> ν(J) は解析の線形の ADTS (R6)、幅は上の振り子の式 (J ごとの G0)。点ごとに連続のスコア log10(W/|Δ0|)、包絡に入る調和の数、信頼度 (0/1/2) と全体の警告を返す。**追跡に合わせた当てはめの項は足していない。**

## 元 147 行 (節: 3.6 信頼して使うための条件と警告 (`pysad/fmapcheck.py`、v2.3.0 S2))

> ### 3.6 信頼して使うための条件と警告 (`pysad/fmapcheck.py`、v2.3.0 S2)

## 元 237 行 (節: 4.4 Touschek への受け渡し)

> **Piwinski 型の式 (v2.3.0 T2、v2 の既定)**: K64 の表の損失率 (`touckf`) は、平坦ビームの関数に ptrans² = tr(P) を入れた近似で、結合が小さくないリングで寿命が ±20% ずれる (B162)。

## 元 248 行 (節: 4.4 Touschek への受け渡し)

> **バンチ内の δ の分布 (v2.4.0 A3。v2.4.0 F2 から v2 の既定。B200)**: **v2 の `TouschekLifetimeLMA[]`・`TouschekLifetime[INF, INF, nz]` (Piwinski 型) と `lma.touschek_lifetime_local_piwinski` は既定で BunchDelta・BunchQs つき (`->False` で切る)。k64 は従来の δ_c = 0 の基準。δ_old = 0 の基準 (K64・Piwinski・elegant) はバンチ内の δ の分布を無視する近似なので、elegant との比較 (T1・T2) の値は「δ_c = 0 の基準との差」である (比較の値自体は `touschek_piwinski` を直接呼ぶので変わらない)。以下の「既定は入れない」は v2.4.0 F2 より前の記述で、`BunchDelta` を省略したときの k64 の動作にあたる。** 台本 `TouschekLifetimeLMA[BunchDelta->True]`・`lma.touschek_lifetime_local_piwinski(bunch_delta=True)`・`touschek_bunch.touschek_piwinski_bunch`。散乱後の粒子の δ は δ_c + Δ (δ_c = 対の重心の δ ~ N(0, σδ²/2)、Δ は相対運動量の散乱だけで決まる) なので、受容を a ∓ δ_c にして平均する。受容が σδ の 8〜14 倍のとき損失率は ASTRID2 で +1.3〜+2.7%、KEK-LS で +1.8〜+4.3% 増える (寿命が短くなる)。実験室系の 2 粒子の対のモンテカルロ (厳密な相対論の運動学) で δ_c の効果が 2% 以内で再現される。(v2.4.0 A3 の時点の既定は入れない): 相対運動量の縦成分 q_s の効果 (MC で通常の kick の判定より損失率が ASTRID2 で 2%・KEK-LS で 0.5% 小さい、つまり絶対値の判定との差がその分大きい) は `BunchDelta->True` だけでは入らない。**v2.4.0 A10: `BunchQs->True` を足すと q_s の効果も入る** (3 次元の積分 `touschek_bunch.qs_ratio`。絶対値 / kick の損失率の比が MC と ASTRID2 +4.81・+2.22% (MC +4.84・+2.22%)、KEK-LS +4.5・+1.7% (MC +4.5・+1.8%) で 0.3% 以内。12 スライスで求めて補間するので全周で 70 秒。ASTRID2 全周で a/σδ = 9・14・23 の寿命が δ_c = 0 に対し δ_c だけ −2.5・−1.2・−0.5%、q_s も入れて −4.1〜−4.3・−1.9・−0.7% (補間の不確かさ ±0.3%))。elegant・K64・Piwinski の式は kick の基準を δ_old = 0 とし、elegant の `touschekLifetime` には δ の分布を入れる機能が無い。a ≲ 3σδ では使えない (警告)。詳細・限界は `pysad/touschek_bunch.py` の冒頭と `docs/history/v24_A.md` §A3、検証は `example/v2_4/A3_mc_bunch.py`。

## 元 250 行 (節: 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2))

> ### 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2)

## 元 262 行 (節: 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2))

>   理想化した解析式 (減衰なし・α 一定) との差 (ASTRID2 で δ+ が +1.5%、δ− が +0.4%) は、放射の減衰・α2 (a2/a1 = −0.76) と離散写像で説明がつき、**LMA の δ+ が δ_rf を 1% 超える (Q1・T1) のは有限の周回数と、この効果**。

## 元 267 行 (節: 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2))

>   **以前の「位置 2 は RF」は check_every = 0 の値 (δ+ = 1.635%、理由 'rf') の取り違え**。`reason` は最初に失われた粒子を落とした直接の原因 ('rf' は |z| の判定が先に効いたとき) で、限界の種類とは別。バケットの限界でも、バケットから出てすべった粒子はアパーチャに当たって 'aper' と記録されうる。

## 元 282 行 (節: 4.5 信頼して使うための条件と警告 (`pysad/lma.py`・`pysad/lmabucket.py`、v2.3.0 L2))

> 振幅 (0.8〜1.6%) が σδ (4.3e-4) の 20〜40 倍なので、量子励起による振幅の変化 ≈ σδ²/δ が 2e-6〜2e-5 と二分の幅より小さいため。大きなリング A (σδ = 6.9e-4、位置 16・89、6 種、tol 1e-4): pysad の FLUC あり (491 周と 1500 周で同じ) は δ+ = 3.027%・3.029% (σ 1.1〜1.3e-4)、δ− = −3.017%・−3.005% (σ 4e-5〜5e-5)、FLUC なしとの差 +0.8e-5〜+8.1e-5 (σ 以内)。elegant の ISR あり 1500 周 (6 種) は δ+ = 3.008%・3.010% (σ 1.0e-4・1.2e-4)、δ− = −3.022%・−3.009% (σ 5.9e-5・7.7e-5)、pysad との差は δ+ で +1.9e-4 (σ の 1.5〜2 倍)、δ− で +0.5e-4〜+0.4e-4 (T1 の ISR なしの差 +5.7e-5〜+1.05e-4 に近い)。elegant の ISR あり 491 周は位置 89 で δ± = +3.19%・−3.25% (σ 7e-4・1.4e-3) と収束しない (バケットから出た粒子がアパーチャに届かないため。T1 と同じ理由。pysad はバケットの判定があるので 491 周でよい)。

## 元 303 行 (節: 5. 検証)

> | L2 (§4.5) | check_every・バケットとアパーチャの区別・窓・周回数・1 セルから全周・FLUC・誤差・MPW | [`example/v2_3/L2_README.md`](../example/v2_3/L2_README.md)、[`docs/history/v23_L2.md`](history/v23_L2.md)。「RF バケット」の行の δ+ の超過の原因と、位置 2 の分類は §4.5 (b) |

## 元 307 行 (節: 5. 検証)

> | elegant との比較 (T1) | ASTRID2 (全周 199 位置の elegant `momentum_aperture`、pysad は 1 セル 25 位置、check_every = 1、±10 mm、1312 周)、kek-ls_v03_12.lat (elegant 14 位置・1500 周、pysad 52 位置・491 周)。RF・放射つき (ISR なし)、同じアパーチャ・同じ周回数 (ASTRID2)。Touschek は elegant の `touschekLifetime` (Piwinski の式) に pysad の IBS 込みのビームを渡す | **LMA: δ± の相対差 0.5% 以内** (ASTRID2: アパーチャが限界の δ = 0.83% も RF バケットで決まる 1.60% も、最大の差 7e-5。kek-ls: 全位置が RF バケット、0.35% 以内)。**Touschek (一様な δ、表の点): elegant / pysad = −4.3〜+2.4% (ASTRID2)、−3.4〜+1.2% (kek-ls)**。局所の LMA の寿命: ASTRID2 pysad 9502 s・elegant 9069 s (−4.6%)、kek-ls 26876 s・27117 s (+0.9%)。差は LMA でなく寿命の式 (pysad は SAD の `touschek_f` の当てはめと等方な運動量分散の近似) による。elegant の `momentum_aperture` はバケットから出ただけでは損失にしないので、周回数が足りないと δ+ が見つからない (kek-ls は 1500 周)。[`example/v2_3/T1_README.md`](../example/v2_3/T1_README.md)、[`docs/history/v23_T1.md`](history/v23_T1.md)。**T2: 独立な式 (Piwinski 型、3 次元積分・厳密な運動学のモンテカルロで確認) では elegant (標本化を細かくしたとき) と −1.3〜+0.1%、K64 の式が不正確 (B162)。局所の LMA の寿命 (Piwinski 型): ASTRID2 9327 s、大きなリング A 26862 s** |

## 元 327 行 (節: 7. 未確認・限界)

> - 周波数マップは、誤差なしのほかに誤差つき (全周の追跡、§9.4・§10.1)・MPW つき (§9.5・§10.2)・k64 (§9.6)・RFSW ON・RAD ON・FLUC ON (§10.3) でも走ることを確かめた (S2 で複数の種・強さ・種類を統計、KEK-LS も)。**LMA は、誤差つき ASTRID2 (3 種) と MPW 版 ASTRID2 でも全周で走ることを確かめた (§4.5 (g))。補正つき・KEK-LS の誤差つきは v2.4.0 H2 で実走した (§11)。結合の強いリングの補正つきは未確認**。

## 元 334 行 (節: 7. 未確認・限界)

> - **Touschek 寿命の式 (T1・T2)**: K64 の表の式 (`touckf`) は、同じビーム・同じ一様な δ で独立な式 (Piwinski 型) と −14〜+21% 違う (結合 1% で +0〜+6%、丸いビームで ±20%。原因は x-y 面の運動量の分布を平坦ビームの 1 次元の関数と tr(P) で近似していること。B162)。elegant の式は独立な式と −1.3〜+0.1% (elegant は要素の出口だけの標本化で、粗いリングで +1〜+5% 長く出る)。v2 の既定は Piwinski 型。限界: 受容 a が σδ の 10 倍以下では、バンチ内の δ の分布 (a − δ_old) を入れていない (4 つの式とも共通) ことが 5〜20% 効きうる (概算。[`docs/history/v23_T2.md`](history/v23_T2.md) §7)。K1≠0 の BEND の K1 は B' 一様ではなく (1 + hx)∂b/∂x 一様の定義 (B160・docs/11 A59。誤りではない) で、elegant の CSBEND には K2 = −hK1 を足して合わせた。ASTRID2 の LMA には影響が無かった。

## 元 340 行 (節: 7. 未確認・限界)

> - 互換モードは v2 で確かめた。**k64 の周波数マップは S1 で確かめた** (ASTRID2・KEK-LS で v2 と同じ、チューンの差 5e-11、§9.6)。**LMA の k64 は v2.4.0 A4 で確かめた** (ASTRID2 の 8 位置で v2 と δ± が 6e-8 以内。check_every = 1。k64 の最初の tapert の扱い (B33) の違いは check_every = 0 のときだけで、見ていない)。

## 元 341 行 (節: 7. 未確認・限界)

> - 周波数マップの D の絶対値の閾値 (D > −5 など) は周回数・窓で動く (周回数 2 倍で約 −1、§9.2)。S2 で 1024 周への換算 `D_norm` を足した (§3.6、§10.5。512〜4096 周)。RDT・ADTS の予想は 1 次の共鳴だけで、ASTRID2 の D の大きい領域 (νx − νy = 3 など 2 次の共鳴) は予想できない (§9.1)。

## 元 351 行 (節: 8. 実装の構成)

> | [`pysad/fmapcheck.py`](../pysad/fmapcheck.py) | (S2) 誤差の統計・フィードダウンの見積もり (`errors_summary`・`feeddown_estimate`・`center_sextupoles`)、`map_stats`、周回数の換算 (`threshold_for_nturn`)、警告 (`check_setup`・`check_result`) |

## 元 352 行 (節: 8. 実装の構成)

> | [`pysad/fmapres.py`](../pysad/fmapres.py) | (S1) `ResonancePredictor` (RDT・ADTS から共鳴の領域を予想)、`coverage_stats`、`lines_distance`・`other_mn`、`apply_alignment_errors`、`plot_with_prediction` |

## 元 353 行 (節: 8. 実装の構成)

> | [`pysad/lmabucket.py`](../pysad/lmabucket.py) | L2: `rf_bucket_analytic`・`longitudinal_params`・`bucket_limits_1d`・`bucket_limits_6d`・`aperture_limits_4d`・`limit_types` (RF バケットの限界とアパーチャの限界の区別、§4.5) |

## 元 359 行 (節: 9. S1 の検証 (v2.3.0、2026-11-20))

> ## 9. S1 の検証 (v2.3.0、2026-11-20)

## 元 375 行 (節: 9.1 共鳴線と RDT・ADTS の定量的な接続 (a))

>   原因は、6 極の強さを λ 倍した周波数マップ (D > −5 の点: λ = 0、0.5 で 0 点、1 で 10 点、1.5 で 37 点) の通り、2 次以上の共鳴。**6 極の 2 次の共鳴駆動項の解析式が `resdrive.py` に無い** (R4 の担当範囲)。

## 元 405 行 (節: 9.8 S1 の未確認)

> ### 9.8 S1 の未確認

## 元 409 行 (節: 10. S2 の検証: 周波数マップを信頼して使う (v2.3.0、2026-11-21))

> ## 10. S2 の検証: 周波数マップを信頼して使う (v2.3.0、2026-11-21)

## 元 448 行 (節: 10.3 RFSW ON・RAD ON・FLUC ON)

> 直したもの: RAD OFF のときの同期位相を 0 に (放射損失がないので。`Emittance[]` の φs では正味の加速が入り、gs と 1% ずれた。0 にして gs とビット一致)。台本 `FrequencyMap[]` の既定を RFSW・RAD OFF に。

## 元 465 行 (節: 10.7 S2 の未確認)

> ### 10.7 S2 の未確認

## 元 469 行 (節: 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`))

> ## 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`)

## 元 472 行 (節: 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`))

> `MomentumAcceptance[]` (・`FrequencyMap[]`・`NormalFormTune[]`) は、台本が `USE RING NOEXPAND;` と `cell calc noexpand` (NOCALEXP) のときは**補正した成分ごとの K0・位置ずれのまま** TrackSetup を作る (v2.4.0 F2 の修正。既定の CALEXP は同名の素子の代表の値に揃える)。この経路を KEK-LS の 2 つのラティスで確かめた。誤差 = `KEKLS_correct_cod` の誤差 (四極・6 極の位置ずれ 40〜50 μm、強さの誤差 4e-4、回転 100 μrad)、補正 = COD 応答行列の特異値分解の疑似逆行列を 0.5 倍ずつ 10 回 (補正磁石 HD・VD の K0)。LMA は v2、RFSW・RAD ON、FLUC なし、アパーチャ ±10 mm、周回数は自動 (シンクロトロン振動 2 周期)。

## 元 476 行 (節: 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`))

> | v03_12 (3563 要素、603 m、1.4 MV、バケット 2.994%、5 位置) | 誤差なし | 0 | 3.021 (バケット) | −3.016 (A4 と同じ値) |

## 元 486 行 (節: 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`))

> - 補正した状態は TrackSetup に渡る: NOEXPAND の台本では HD の K0 が成分ごとの値のまま (LMA の後も保たれる)。誤差なしの v03_12 の LMA が A4 の値と一致する (台本の経路と `TrackSetup.from_source` の経路が同じ結果)。対照 (既定の CALEXP) は補正した K0 が代表の値に揃い、位置ずれ・強さの誤差も同名の素子で同じ値になるので、補正も誤差も別の状態になる (v03_68 は補正なしに近い値、v03_12 は誤差なしに近い値)。

## 元 488 行 (節: 11. 誤差つき・補正つきの KEK-LS の LMA (v2.4.0 H2。台本の経路: 誤差 → COD の補正 → 同じ FFS の中で `MomentumAcceptance[]`))

> - **誤差を 1.5 倍 (v03_68) にすると k64 (gs) は補正が始まらない**: K64 の qcod (反復 10 回の Newton 法) に収束の保証がなく、gs でも 1.5 倍は `cell calc` の COD がゼロ (設計のチューンのまま。docs/21 の B260)。**v2 は堅牢な解法にフォールバックして COD を見つけ (1.5 倍で 1.9・3.4 mm、2 倍で別の安定な不動点)、補正が走る** (v2.4.0 H4。`docs/history/v24_H.md` H4-1)。この節の LMA の表は 1 倍の誤差 (K64 が収束する) なので影響を受けない。

