# 変更履歴

バージョンは `vメジャー.マイナー.パッチ` の形です。メジャーの区切りに GitHub のリリースを作ります。

**バージョンの方針**
- **v1 系: ソースに忠実な翻訳。** SAD (K64 版) のソースの動作を、そのまま Python に写すことを目的とします。K64 のバグや癖も、gs の実動作に合わせて再現します (癖は [`docs/11_sad_quirks.md`](docs/11_sad_quirks.md)、K64 側の不具合の記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md))。v1.x のマイナー・パッチは、翻訳の漏れの追加と、翻訳の誤りの修正です。
- **v2 系: ソースのバグの修正と機能の追加。** K64 のバグを直し、ソースにない機能を加えた版です。gs と同じ結果にならない箇所が出てきます。v1 との違いは、リリースごとにここに書きます。

## v2.3.0 (2026-11-21)

v2.2.0 の後に、K64 にない解析機能 (周波数マップ・局所運動量アクセプタンス (LMA)・共鳴駆動項 (RDT)・振幅依存チューンシフト (ADTS)・Touschek の独立な式・数値の正規形) を足し、それぞれの信頼性を調べた版です (PyPI は v2.2.0 まで公開。v2.3.0 は未公開)。K64 の動作・互換モード (`pysad/compat.py`、既定 `v2`、`PYSAD_COMPAT=k64`) は、下の「既定の変更」のほかは変えていません。**新機能は K64 に無い追加なので、互換モードに関係なく同じ動作で使えます**。numba のカーネルは変えていません。K64 について B 番号を 3 つ足しました (B160: K1≠0 の BEND の K1 の暗黙の定義。誤りではない。B161: CALC の `tbendecorr` の行列の dk の不整合。v2 で直した。B162: Touschek の損失率の式の適用範囲。v2 は Piwinski 型の式を既定にした)。K64 のバグの記録は B01〜B162。

**解析式の位置づけ**: RDT・ADTS・共鳴線の予想などの解析式は、摂動論の定義された次数と、入れた効果の範囲で原理的に作ります。追跡 (全非線形) とは入れている効果が違うので、厳密な一致は目的ではなく、追跡に合わせるための補正項・当てはめは足していません。追跡との一致の数値は、入れた効果と次数の検証の結果として載せ、入れていない効果と差の主因を [`docs/49`](docs/49_resonance_terms_adts.md) §11.9 に書きました (原因が分からないものは「未特定」)。追跡から係数を取り出す唯一の部分は偏向磁石の生成関数 (`bendnl.py`) で、原理は「写像の 2 次のテイラー係数が生成関数を一意に決める」です (docs/49 §8.1)。
変更の一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.3.0 の節、解説は [`docs/48_frequency_map_lma.md`](docs/48_frequency_map_lma.md) と [`docs/49_resonance_terms_adts.md`](docs/49_resonance_terms_adts.md)、作業の記録は `docs/history/v23_Q1.md`・`v23_Q2.md`、例は [`example/v2_3/`](example/v2_3/) です。

**新機能** (使える範囲・条件は docs/48・49 と docs/47 の v2.3.0 の (b))
- **周波数マップ** (`pysad/freqmap.py`): 1 つの位置から (x, y)・(x, δ)・(y, δ) の格子を追跡し、NAFF で前半・後半のチューンを出して D = log10 √(Δνx² + Δνy²) の地図を作る。損失点は D = nan と損失周回数。表・図 (チューン平面 + 共鳴線、初期条件平面)、共鳴線の表 (整数・半整数を必ず含む)。
- **局所運動量アクセプタンス (LMA)** (`pysad/lma.py`): 位置ごとの δ+(s)・δ−(s) を、粗い走査 + 二分法、放射・RF つき、RF バケットの外も損失、周回数の既定はシンクロトロン振動 2 周期、並列 (既定 2) で求める。1 セルだけの計算と周期への延長 (`trksetup.find_periodicity`)。Touschek への受け渡し (`touschek_lifetime_local`)。
- **共鳴駆動項・ADTS・クロマティシティ** (`pysad/resdrive.py`): 6 極 (SEXT・K2≠0 の MULT) の 1 次の共鳴駆動項 h_jklm (10 個。要素ごとの寄与つき)、振幅依存チューンシフト (6 極の 2 次 + 8 極の 1 次)、クロマティシティ (数値微分と 6 極の 1 次の寄与)。厚いレンズの解析式 (既定) と数値積分。
- **追跡の設定の共通部分** (`pysad/trksetup.py`): 閉軌道の求め直し、周期性の検出、位置の選び方、RF バケットの限界、優先度を下げた並列。

**使える関数**
- SAD の台本 (pysad のみ。gs には無い): `FrequencyMap[位置, nturn, X->…, Y->…, DP->…, Output->…, Plot->…]`、`MomentumAcceptance[位置, Cell->…, CheckEvery->…, Workers->…, …]`、`TouschekLifetimeLMA[]`、`ResonanceDrivingTerms[Terms->…, Method->…, Slices->…]`、`AmplitudeTuneShift[…]`、`ChromaticityTerms[Order->…, …]`。
- Python: `pysad.freqmap.frequency_map`・`resonance_check`・`resonance_lines`、`pysad.lma.momentum_acceptance`・`touschek_lifetime_local`、`pysad.trksetup.TrackSetup`・`find_periodicity`・`select_locations`、`pysad.resdrive.resonance_terms`・`adts`・`chromaticity_numeric`・`chromaticity_sextupole` (API は docs/12 §16.8)。

**確かめたこと**: 周波数マップの追跡は gs の `TrackParticles` と座標が 2.8e-7 以内・D が 4 桁一致。LMA の損失判定は gs と 30 / 30 一致。1 セルの LMA と別のセルの差は最大 5e-5。LMA の δ+ は RF バケット高さと ±3% 以内、ASTRID2 の線形の見積もりと 0.96 倍。Touschek の局所の式は一様な受容で既存の表と 1e-10。ADTS は小さな FODO で解析式・pysad の追跡・gs の追跡が 4e-5〜2e-4 以内、ASTRID2・大きなリング A で pysad の追跡と 1e-3 以内。共鳴駆動項は 1 周の ΔJ の追跡と 0.5%。8 極 1 次は薄い式と 1e-9。厚いレンズ ⇔ 薄レンズ列は 1/n²。

**R2 の追加 (ADTS の直接の検証)**: 実際の強さの追跡 (pysad k64・gs) + NAFF + J の多項式フィットで ADTS を求める `pysad/adtsdirect.py` (`direct_adts`) を足し、解析の ADTS と比べた。Q2 の C の比較は交差項を消していたため、全部入りの追跡の ADTS が 6 極 2 次 + 8 極 1 次と最大 34% 違うことに気づかなかった。原因の運動学的な項 p⁴/8 と四極 (K1 つき偏向磁石) の端部の非線形フリンジの 1 次を `resdrive.adts` に足した (`.kin`・`.fringe`・`.total`、台本の `Kinematic->`・`QuadFringe->`・`Total->`)。追跡との差は FODO +2.8%・−1.7%・−8.6%、ASTRID2 −9.9%・−13%・+7.0%、大きなリング A −0.19%・−1.2%・−3.7%、ASTRID2 MPW −9.9%・−13%・+8.5%。残りは 6 極 × (偏向磁石の曲率・硬い端) の交差項 (原因特定済み、未実装)。gs と pysad の追跡のフィットは 1e-10。文献 (Bengtsson 型) の h との符号 (pysad = −文献)・規格化の照合、共鳴線の近くの周期点 (0.15〜1.3%)、8 極の走査 (2e-6〜2e-5)、δ 依存・結合・MPW も確認 (docs/49 §7、`docs/history/v23_R2.md`、`example/v2_3/R2_README.md`)。

**R3 の追加 (偏向磁石の 3 次の項)**: 偏向磁石の曲率 (h/2)x(px²+py²)・硬い端 ±(h/2)px y²・複合機能の本体の K1 h x³/6・E1/E2≠0 の端の 3 次の項を `resdrive.adts` に足した (`.bend` = 3 次どうしの対の項、`.cross` = 6 極との交差項、`.total` に入る。`bend_cubic=True` が既定、台本 `AmplitudeTuneShift[BendCubic->…, BendSlices->…]`)。係数は K64 の追跡の写像の 2 次のテイラー展開から取り出す新しい `pysad/bendnl.py` (偏向磁石を入口の端・16 個のスライス・出口の端に分け、硬いフリンジは隣り合うスライスで打ち消す)。追跡 (pysad k64 = gs の直接のフィット) との差 (6 極 ×1、解析 − 追跡): FODO −0.04%・−0.04%・−0.03% (R2 −2.7%・+1.6%・+9.5%)、ASTRID2 −0.24%・−0.08%・−0.19% (R2 +10.9%・+15.3%・−6.5%)、ASTRID2 MPW −0.25%・−0.08%・−0.22%、大きなリング A −0.014%・−0.017%・−0.26% (R2 +0.19%・+1.2%・+3.6%)、長方形・非対称の端の FODO も 0.1% 以内。6 極の強さ ×0〜×2 の走査で奇数次・偶数次の両方が合う。残りは 6 極に依らない 4 次の定数 (四極・K1 つき偏向磁石の本体、未解決。ASTRID2 で 1.1・1.1・2.6、大きなリング A で 7.0・6.9・17 [1/(m rad)])。gs 3 例 (ASTRID2 ×0.5、FODO ×2、大きなリング A ×2) で gs = pysad k64 (1e-9)。K64 の誤りなし (B 番号の追加なし)。K64 の K1≠0 の偏向磁石の本体の 3 次の項が Maxwell 方程式から導いた形と違う点は要判断として記録 (docs/49 §8、`docs/history/v23_R3.md`、`example/v2_3/R3_README.md`)。

**R4 の追加 (共鳴駆動項の追跡による検証と偏向磁石の項)**: 共鳴駆動項 h を追跡の 1 周の写像の 2 次のテイラー展開から取り出す `pysad/rdttrack.py` を足し (kick の和なので 1 次の摂動で厳密。δ の項は閉軌道のヤコビ行列の δ 微分、3νx の近くの不動点も)、`resonance_terms` を FODO・FODO 長方形・非対称・ASTRID2・ASTRID2 MPW・大きなリング A で検証した。6 極の分は 1e-6 以内で正しい (δ の項 h11001・h00111・h20001・h00201・h10002 も 1e-7 以内)。偏向磁石の 3 次の項は h の 1.4〜3.4% (FODO)、0.8〜11% (ASTRID2)、最大 22% (MPW)、0.005〜0.3% (大きなリング A) あったので `resonance_terms(bend_cubic=True)` (既定) に足し、追跡と 1e-8 (FODO・ASTRID2)・4e-7 (大きなリング A)・5e-5 (MPW。K0≠0 の偏向磁石の内部) 以内。**既定の変更**: 6 極のみの値は `bend_cubic=False`。δ の項は偏向磁石の分散による軌道のずれの展開で、自然クロマティシティの解析式 `chromaticity_natural` + 偏向磁石の分 `chromaticity_bend` が CALC の ξ (6 極なし) と 3e-3 (FODO)・6e-4 (ASTRID2)・1.4e-5 (大きなリング A) 以内で一致 (台本 `ChromaticityTerms[]` の `NaturalPart->`・`BendPart->`、`ResonanceDrivingTerms[BendCubic->…]`)。3νx の近くの不動点 (対称なリングでは曲線) の J(φ) の 3φ 成分から求めた |h30000| は解析と 0.1〜0.6%、gs の 1 周の写像から取り出した h は pysad k64 と 6e-10 以内。K64 の誤りなし (B 番号の追加なし)。互換モードの分岐は不要 (新機能)。未確認: K0≠0 の偏向磁石の内部、δ≠0 の h の偏向磁石の項、4 次以上の h、結合・誤差 (docs/49 §9、`docs/history/v23_R4.md`、`example/v2_3/R4_README.md`)。

**R5 の追加 (4 次の共鳴駆動項・6 極の 2 次)**: `pysad/rdt2nd.py` (4 次の h_jklm (35 項、p=0) の解析式: 6 極の 2 次 (Magnus の第 2 項の順序つきの対の和、厚い 6 極の内部は閉形式)・8 極 1 次・運動学的な項・四極端部のフリンジ・偏向磁石の 3 次どうし/6 極との交差項・偏向磁石の 4 次 (要素の写像から)。`resonance_terms(["h22000","h31000"], order2=True)` (既定)、台本 `ResonanceDrivingTerms[Terms->{"h22000"}, Order2->True]`)、`rdttrack.ring_map3`・`rdt4_from_map` (追跡の 1 周の写像の 3 次のテイラー展開から Dragt–Finn の f4)、`pysad/nform.py` (数値の正規形、J の 1 次)、`fmapres` の `order2=True`。ADTS との整合: 写像の対角部分 + 正規形の補正 = `adts()` が 1e-15。追跡の f4 との最大相対差は FODO 1.3e-4・ASTRID2 1.2e-3・MPW 1.8e-4・KEK-LS 5.3e-6 (薄い 6 極の列の FODO は 2e-6。残りは K64 の厚い 6 極の積分器)。**R3 の「未解決」の 6 極に依らない 4 次の定数の原因は K64 の四極の積分器 (台形則、`kinematic_h4(quad_rule="k64")`)** で、入れると ADTS の追跡 (NAFF) との差が FODO −3.6e-4 → −8.8e-5、ASTRID2 −2.4e-3 → −6e-5〜−2e-4、KEK-LS yy +2.6e-3 → +8e-6 (`adts()` の既定は R3 のまま、`rdt2nd.adts_from_rdt4` で出る)。周波数マップ: ASTRID2 は 2 次で領域が現れる (2νx−2νy = νx−νy=3、4νy) が D > −5 の点の 48% が領域に入るだけで対照 (30±13%) と有意に違わない (p=0.2。原因: νx−νy≈3 の帯が格子の 30%、次数 5 以上の共鳴、非対角の調和の正規形の補正なし)。K64 の誤りなし (B 番号の追加なし)。未確認: 非対角の調和の正規形の補正、E1/E2≠0 の端の 3 次の生成関数の解析式 (取り出しのみ)、δ の 4 次の項、K0≠0 の偏向磁石の内部、J² の正規形・大振幅 (docs/49 §8.1・§10、`docs/history/v23_R5.md`、`example/v2_3/R5_README.md`)。

**S1 の追加 (周波数マップの検証)**: `pysad/fmapres.py` (1 次の RDT と ADTS から共鳴の領域を予想する `ResonancePredictor`、統計 `coverage_stats`、誤差の与え方 `apply_alignment_errors`)。FODO の 3νx では D > −5 の 44 点がすべて予想の領域に入り、共鳴の中心が D の尾根と 0.5 mm 以内。ASTRID2 は D の大きい点が 2 次の共鳴 (νx − νy など) で予想できない。δ を含む格子・誤差つき・k64 (v2 と同じ)・挿入光源つきで周波数マップが走る。D の閾値は周回数で動く (2 倍で約 −1)。KEK-LS の gs との比較 (規則的な点は 1e-11 m)。台本ごとのまとめは `example/v2_3/FMAP_summary.md`、記録は `docs/history/v23_S1.md`、docs/48 §9。

**T1 の追加 (elegant との比較)**: pysad の LMA と Touschek 寿命を elegant (2026.4.0、`momentum_aperture`・`touschekLifetime`) と同じ条件で比べた (`example/v2_3/T1_*`、`docs/history/v23_T1.md`、`pysad/tests/test_v23_T1_elegant.py`。pysad のコードは変えていない)。SAD → elegant の変換 (`T1_sad2ele.py`): 光学 (チューン・αc・U0・自然エミッタンス・σδ) が合い、色収差は K1≠0 の BEND に K2 = −h K1 と EDGE_ORDER = 2 を足して合わせた。**LMA**: ASTRID2 (全周 199 位置の elegant と pysad の 1 セル) は δ± の相対差 0.5% 以内 (アパーチャが限界の位置も RF バケットで決まる位置も)、kek-ls_v03_12.lat は全位置が RF バケットで 0.35% 以内 (elegant はバケットから出ただけでは損失にしないので 1500 周が要った)。**Touschek 寿命**: 同じビーム・同じ一様な δ で elegant / pysad = −4.3〜+2.4% (ASTRID2)、−3.4〜+1.2% (kek-ls)。差は LMA でなく式 (SAD は `touschek_f` の当てはめ + 等方な運動量分散の近似、elegant は Piwinski の式)。どちらが正確かは未確認。**B160**: K1≠0 の BEND のオフモメンタムの光学 (一様な磁場指数のリングで gs・pysad の dν/dδ = (−0.172、+0.211)、B' 一様の解析式・elegant (−0.602、+0.738))。U1 で「誤りではなく K1 の暗黙の定義」に再分類 (下)。

**U1 の追加 (B160 の切り分け)**: K1≠0 の BEND の本体は、中央面の磁場 B_y = B0 + (K1/h) ln(r/ρ0)、すなわち曲線座標の勾配 (1 + hx)∂b/∂x が一様な定義 (B'' = −hB') を実装している。独立なローレンツ力の積分と 16 台の非線形の写像が 1e-10 で一致し、B' 一様とは 6 極 K2 = +hK1 の差。**誤りではなく定義**で、SAD → elegant・MAD の変換には K2 = −hK1 を足す (docs/11 A59)。pysad のコードは変えていない (`docs/history/v23_U1.md`、`test_v23_U1_bend_k1_field.py` 7 件)。要判断: 既定を変えない (推奨)、`$B'$` 一様のオプションを足すか。B161 (疑い、CALC・Emittance[] の行列だけ)。

**R6 の追加 (入れていなかった物理の項と誤りの修正)**: K0≠0 の偏向磁石 (挿入光源 = MPW) の項 (閉軌道のまわりのテイラー展開。h の 3 次 4.6e-5 → 2.5e-9)、δ≠0 の ADTS で **kick の強さを 1/(1+δ) 倍していた誤り** (R2 から。追跡の kick は p0 で規格化した場の強さで δ に依らない。修正して残した)、E1・E2≠0 の端の生成関数の閉形式 (K1=0 で 1.6e-8)、非対角の調和の正規形の補正 (独立な数値のリー変換と 1.2e-8)。`adts()` の精密版 `fourth='r6'` (四極の積分器の台形則 + 偏向磁石の 4 次の項) と `optics='tracked'` (δ≠0 の光学を追跡の閉軌道から) を足した。

**R7 の追加 (6 極の検討、試作)**: `pysad/sxtopt.py` (試作。20 ファミリーの K2L に対する 1 次の RDT・δ の項 (線形)・4 次の RDT と ADTS (二次形式) の多項式モデル、解析的な勾配、SLSQP の最適化。`SxtModel.check` で解析式と 1e-9 以内)、例 `example/ddba/v07_14/sextupole/` ([README](example/ddba/v07_14/sextupole/README.md)。v07_14 (5 GeV、周長 749.48 m) の6 極 20 ファミリー、別の研究者による最終設定との比較)、試験 `test_ddba_sextupole.py` 8 件。結論: 6 極の組 32 通りで、解析式の目的関数と追跡の x–y DA 面積の順位相関は ρ = −0.89。検証 (小〜中振幅の DA の健全性チェック) には条件つきで使える。LMA の予測 (最終設定の近くで |ρ| ≲ 0.2) には使えない。**最適化の試作は、最終設定より DA の広い解を出せなかった**: 解析式の目的関数 (物理的な重み。DA に合わせた調整なし) を最小化した解は、目的関数が最終設定の約 1/8 なのに、DA 面積は最終設定の 0.60〜0.71 倍 (劣る。数値の揺らぎ 0.3〜1.1% より大きい。簡略化した追跡: 512 周、誤差なし)。不足は、J² の ADTS (5・6 次の写像)、δ 依存の項、5 次以上の駆動項。 試作の最適化は走り、制約 (ξ、|K2L| ≤ 最終設定の最大) を満たし、40 通りの初期値が同じ底に収束するので、障害は最適化の収束でなく、目的関数と DA の関係の弱さ。R7 の追跡は簡略化した設定 (512〜1024 周、光線の走査、6 極の位置ずれ rms 50 μm の 3 種、補正なし) で、完全な設定の DA・LMA での比較は未実施 (README §7 の次の課題 6 件)。

**V2 での `adts()` の既定の変更**: R6 は `fourth='r6'` を既定にしたが、これは数値積分器の端数 (0.1〜0.3%) を含むので、**既定は物理的な項だけ (`fourth='r3'`・`optics='calc'`) に戻した** (R6 は精密版のオプション)。既定の追跡との差は FODO 4e-4、ASTRID2 2.4e-3、KEK-LS yy 2.6e-3 (主因は入れていない 4 次の積分器の分。`fourth='r6'` で 1e-4 以下)。あわせて、r3 では偏向磁石の 3 次の項を原点のまわりの分解で取る (K0≠0 の磁石だけ閉軌道のまわり。`bendnl.ORBIT_TOL`) ことで、CALC の呼び出しの履歴に残る 1e-10 の雑音で ADTS が 1.4% 変わる不安定さを除いた。

**S2 の追加 (周波数マップの信頼性)**: `pysad/fmapcheck.py`。誤差つき (ASTRID2 の D > −5 の割合は誤差 ×1 で 0.47、誤差なし 0.03。原因は差の共鳴 νx − νy = 3 を誤差の結合が駆動すること)、MPW の失われる点の原因 (内部の非線形による y の力学口径の縮小。4νy = 9 ではない)、RFSW ON (D は拡散の指標でなくなる)、D の周回数の換算 `D_norm`、警告 12 種、gs と誤差つき・RFSW ON の追跡の一致。**既定の変更**: `FrequencyMap[]` は RFSW・RAD・FLUC OFF が既定、RF ON・RAD OFF の同期位相は 0 (gs と同じ)。

**S3 の追加 (共鳴線の予想の作り直し)**: `pysad/fmappred.py`。3 次・4 次の h と非共鳴の正規形の補正、線形の ADTS、振り子の島の幅から連続のスコアと信頼度を出す。追跡から J² の項を当てはめる版は、補正項を足していく方針に反するので削除した。FODO の 3νx の島では AUC 0.94 だが、ASTRID2 (0.68)・KEK-LS (0.59) は振幅だけのスコア (0.78・0.78) に及ばない (次数 5 以上・島の中は規則的で D の大きい点は縁と重なり・大振幅で J² の項が効く)。

**L2 の追加 (LMA の信頼性)**: `pysad/lmabucket.py`。**`check_every` の既定を 1 (要素ごと) に変更** (0 は SAD の tapert と同じ 15 要素ごとの判定で、ASTRID2 の 25 位置で δ+ を最大 1.46 倍に見積もる。警告つきで残る)。RF バケットの限界を解析式・1 次元の写像・6 次元の追跡で確認 (3e-5〜1e-4)、RF バケットとアパーチャの区別 (`limit_types`)、`verify`・`turncheck`、周期性の確認、警告 9 種。

**T2 の追加 (Touschek 寿命の式)**: T1 の「どちらが正確か」を解決した。独立な Piwinski 型の式 (`touschek.touschek_piwinski`) は 3 次元の一般の積分・実験室系の対のモンテカルロと 0.3% 以内、elegant と −1.6〜+0.3% で一致し、K64 の式 (`touckf`) は結合 1% で +0.1〜+6%、丸いビームで −14〜+21% 違う (B162)。**v2 の `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定を Piwinski 型にした** (`Formula->"K64"`・k64 モードは K64 の表)。

**U2 の追加 (B160 のオプションと B161 の修正)**: `compat.set_bend_k1_def('lin')` (B' 一様。elegant・MAD と同じ。既定 'log' = K64。v2 だけ)。**B161 は v2 で直した** (`compat.bend_dk_fix()`。k64 は K64 のまま。実際のリングへの影響は小さい)。

**pysad 自身の不具合の修正 (パーサ)**: 逆向きの文脈の `-LINE` の二重反転 (PF の履歴ラティスで発見)、k64 モードの数に隣接する単位 (`1.7MV`・`3GEV`・`10MM`) の構文エラー (gs は読む。D で始まる単位 `45DEG` だけが構文エラー。DDBA_06 で発見)。

**新しい example**: PF の履歴ラティス (`example/pf/lattice/`、1981〜2005 の 4 つ、`PF_history.md`)、KEK-LS の例集 (`example/kek_ls/`、v03_68、3 GeV、周長 570.72 m。設計・マッチング・パラメータ表示・誤差と COD 補正)、DDBA の例集 (`example/ddba/`、DDBA_01〜28 の版の設計手順、5 GeV・2.5 GeV、周長 731〜757 m)、DDBA の最終版 v07_14 (`example/ddba/v07_14/`、チューンのマッチング・チューンサーベイ・誤差を段階的に大きくしながら COD を補正する誤差の計算・パラメータ表示)。参考文献: KEK-LS CDR (2017 年 5 月 22 日)、Hybrid Ring の論文 (DOI 10.1107/S1600577521012753)、ASTRID2 のデザインレポート。個人名・内部のフォルダ名は例に残さない。

**テスト基盤**: `pysad/tests/conftest.py` に、各試験ファイルのあと残った Tk のルート窓を閉じる autouse の module fixture を足した (GUI の試験の内容・結果は変えない)。

**要判断・未確認の経緯 (R5・T1)**: R5 の `adts()` の既定に四極の台形則・偏向磁石の 4 次を入れるかは、V2 で「入れない (精密版のオプション)」にした。T1 の Touschek 寿命の式の差は T2 で解決。R5 の未確認 (非対角の調和の正規形の補正・E1/E2≠0 の端・K0≠0 の偏向磁石の内部) は R6 で解消。一覧は docs/46・47 の v2.3.0 の (d)。

**試験の状態** (2026-11-21、V3 の最終の通し): 全試験 **1899 通過・8 skip・失敗 0 (1907 件収集。R7 の試験 `test_ddba_sextupole.py` 8 件を含む)**、回帰 `QUICK=1` 50 ベンチ FAIL 0 (133 s)。skip の内訳: PYSAD_SLOW=1 が要る 4 件・PYSAD_ERR_LATTICE が要る 2 件・POSIX/fork が要る 2 件。以前の記録: V1 (2026-11-20) は 1759 通過・8 skip・失敗 0、R5 の後は 1752 通過、R4 の後は 1718 通過、v2.2.0 は 1647 通過・9 skip、v2.1.0 は 1345 通過・8 skip、v2.0.0 は 937 通過・4 skip。`test_plot.py::test_fitgaussian_plot_matches_gs` は回帰と同時だと PermissionError で落ちることがあり、単独で通る。

**条件と限界**: 使ったリングは ASTRID2 (0.58 GeV、45.7 m)、大きなリング A (3 GeV、周長 603 m。KEK-LS は ESRF 型の HMBA (Hybrid Multi-Bend Achromat) で、ESRF の Pantaleo Raimondi 氏と Simone Liuzzo 氏の協力で設計した。)、小さな FODO (3 GeV)。いずれも**誤差・結合なし**。誤差・結合のあるリング、FLUC つきの LMA、PF リング、挿入光源つきのリングは未確認。k64 モードでの LMA・周波数マップは未確認。LMA の損失判定は要素の境界 (既定は SAD の tapert の 15 要素ごと。アパーチャが限界だと δ の限界を 1.2〜1.3 倍に見積もるので `check_every = 1〜5` を勧める)。共鳴駆動項・ADTS は結合・垂直分散・6 極 × 8 極・8 極の 2 次・10 極以上・6 極に依らない 4 次の項 (四極・K1 つき偏向磁石の本体、ADTS の 0.02〜1.3%) を入れない (運動学的な項・四極の端部のフリンジは R2、偏向磁石の 3 次の項は R3 で入れた)。文献の h_jklm の符号・規格化との照合は未実施。

**未解決・未確認** ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.3.0 の (d)、[`docs/46`](docs/46_unimplemented_and_unverified.md)): 精密版の δ 微分の ASTRID2 の yy が 28% 外れる (原因は未特定)、K1≠0 の偏向磁石の本体の閉形式の 1e-4〜3e-4 (K64 の積分器と推定、未照合)、共鳴線の予想 (S3) が ASTRID2・KEK-LS で振幅だけのスコアに及ばない、k64 モードの LMA・PF リング・SP8 (大きなリング B)、Touschek のバンチ内の δ 分布の効果、**ユーザー定義の関数名 `At` が pysad の `@` と衝突して呼べない** (gs は呼べる)、**`FitFunction` を `fit` と併用すると pysad が gs と別の解に落ちる** (KEK-LS の例)、DDBA_15 の εx (gs と +0.33%)、DDBA v07_14 の 5 GeV の誤差計算の DA (gs 0〜1σ、pysad 27〜38σ)。

## v2.2.0 (2026-11-20)

v2.1.0 の「要判断・未解決」を追求した版です。PyPI: [`pysad-k64` 2.2.0](https://pypi.org/project/pysad-k64/2.2.0/) (2026-11-20 公開。`pip install pysad-k64`)。GitHub のリリースは v2.2.0 (Latest)。互換モード (`pysad/compat.py`、既定 `v2`、`PYSAD_COMPAT=k64` で K64 (gs) の動作) は v2.1.0 と同じで、k64 モードの結果を変えたのは、K64 の動作の追加 (これまで未対応・無視だったものを gs と同じに実装) と pysad 自身の不具合の修正だけです。
変更の一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.2.0 の節、K64 側のバグの記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (v2.2.0 で新しく見つけたものは B135〜B159)、作業の記録は `docs/history/v2next_P1.md`〜`P5.md`、例は [`example/v2_next/`](example/v2_next/) (`P1_*`〜`P5_*`) です。

**概要**: 偏極 (n が傾く・非平面のリング、脱偏極率の追跡との差の原因)、放射 (SynchroBeta の 40 倍の原因、tradke の符号の確定と微分の 4 点の誤り)、TrackParticles の F1 と MULT の AE、空間電荷 (k64 の PSPAC)、DAPWIDTH の窓、MEASURE の k64、複素の特殊関数を追求し、K64 の誤りを 25 件新しく見つけた (B135〜B159)。

**主な修正 (v2 で結果が変わるもの。K64 の動作は k64 モードで再現)**
- 偏極: 緩和率を厳密式 Σ r (1 − 2/9 (n·v)²) に (n が傾く・垂直の偏向・水平磁場の挿入光源。B116)、運動論的偏極と放射の刻みごとの脱偏極率。脱偏極率の追跡との差 (v2.1.0 は 0.43 倍) の原因はベータトロンのチューンの δ 依存で、クロマティックな変調を側帯の式に入れて 4.2 GeV・12 MV で式 / 追跡 = 1.05 ± 0.16。不安定なリングの P は NaN。追跡のスピンが BEND の ROTATE で回らない (B135)。
- 放射・エミッタンス: SynchroBetaEmittance の始点の分散による数倍〜40 倍の過大を修正 (B128 の原因の特定)。tradke の符号は放射反作用つきのローレンツ力の積分で確定 (B111)、RADTAPER の BEND の係数を Emittance[]・追跡とも (δ_l + δ_l+1)/2 に統一、tradke の微分の 4 点の誤りを修正 (B138)。
- 端の角: 追跡 `tsteer` の出口の AE2 (B140)、MULT の出口の蹴りの位置 (B141)、K1≠0・K0≠0 の ANGLE=0 BEND の追跡の AE (B142)、MULT・追跡の端の角の放射の行路 (B143。B108 の拡張)。
- SOL: 領域内の QUAD・MULT のスピンを Thomas-BMT の積分と 9e-7 以下に (B146)、領域内の DRIFT・QUAD・MULT の光子 (B150)。
- 光学: 入れ替わる (DETR ≥ 1) INS の一周の位相の整数部 (B153)、補正磁石の ROTATE + DX・DY の CALC (B144)。

**新しく使える機能・オプション** (使える範囲・条件は docs/47 の v2.2.0 の (b))
- `spin.ST_EXACT`・`KINETIC`・`STEP_RATE`・`CHROMATIC` (v2、Polarization)。`compat.set_bend_k1_edge('face')` (K1≠0 の BEND の端の模型。環境変数 `PYSAD_BEND_K1_EDGE=face`。既定は `radial`)。`WSPACSC = 1;` (追跡の WSPAC に空間電荷つきのビーム行列。既定は 0)。
- TrackParticles の F1・FB1・FB2・FRINGE (ANGLE=0 の BEND・MULT)、MULT の AE1・AE2 (k64 は gs と全桁一致)。k64 の PSPAC (K64 の形状関数。`pysad/spch_k64.py`。NPARA = 1 で gs と rms が 12 桁以上)。k64 の MEASURE・TRCOD (先頭の位置。gs と統計が 4 桁一致)。DAPWIDTH の窓 (`trackd.f` の直訳。**既定を 51 から 7 に変更**)。要素ごとの FREE の何もしない変数 (gs と全桁)。
- 複素の特殊関数 15 関数 (`pysad/specfn_cx2.py`。732 式の走査で gs と 1e-9 が 436 → 641 式)。MAIN の文法 (`SIN[0.1]` を構文エラーに、`E`・`Pi`・`Degree` の代入)。gs との照合: `OYStability` の後半・Micado の変種・FITP>1・陽子の U0・Python 3.14 + Tk 9.0 の GUI。

**pysad 自身の不具合の修正**: DAPWIDTH を無視して全点を追っていた、MULT の AE を黙って無視、複素引数で虚部を黙って捨てる特殊関数、不安定なリングで P を返す、`SIN[0.1]` が値 0 など ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.2.0 の (c))。

**試験の状態** (2026-11-20): 全試験 1647 通過・9 skip・失敗 0 (1656 件収集。`PYSAD_SAD_PACKAGES` つき)、回帰 (`QUICK=1`) は 50 ベンチ FAIL 0 (138 s)。v2.1.0 は 1345 通過・8 skip。v2.2.0 の追加は担当 P1〜P5 の `test_v2_P*.py` 19 ファイル。詳細は docs/47 の v2.2.0 の (e)

**条件と限界**: v2 の値の根拠は解析式・独立な数値計算・保存則で、gs ではない。確かめたのは 0.2〜4.2 GeV の FODO リング、大きなリング A (3 GeV、周長 603 m)・PF リング・ASTRID2 (MPW 入りを含む) などで、SP8 (大きなリング B) は数値だけ。SP8 (大きなリング B) は SPring-8 型で、デザインレポートより再現した (8 GeV、周長 1436 m、矩形の偏向磁石 88 台)。脱偏極率の式の追跡との比は 4.2 GeV・12 MV の 1 点。実機の強い結合・挿入光源のずれのあるリングでは未確認。

**未解決**: [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.2.0 の (d)。**判断が要るもの** (既定は担当の選択のまま): EMITDIV の既定 (1)、DAPWIDTH の既定 (7)、`WSPACSC` の既定 (0)、K1≠0 の BEND の端の模型の既定 (`radial`)、放射の β の因子 (B112)、領域内の ψ の DISFRIN、偏向磁石の本体の B∥ (B136)、複素の特殊関数の k64 の値。ほかに B108 の残る 2.7%、RF が分散のある場所にある矩形リングの εx の +1.7e-3、SOL 領域内の加速 (POL + TRPT)、SB のシンクロ・ベータ共鳴の追跡との比較。

## v2.1.0 (2026-11-18)

v2.0.0 の「未対応・未確認」の項目を追求して直した版です。PyPI: [`pysad-k64` 2.1.0](https://pypi.org/project/pysad-k64/2.1.0/) (2026-11-18 公開。`pip install pysad-k64`)。互換モード (`pysad/compat.py`、既定 `v2`、`PYSAD_COMPAT=k64` で K64 (gs) の動作) は v2.0.0 と同じで、k64 モードの結果は pysad 自身の不具合の修正を除いて変えていません。
変更の一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.1.0 の節、K64 側のバグの記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (v2.1.0 で新しく見つけたものは B111〜B134)、作業の記録は `docs/history/v2next_*.md` (A1・A2・B1・B2・B3・C・D1・D2・E)、例は [`example/v2_next/`](example/v2_next/) です。

**概要**: 放射・スピン・光学・空間電荷・パッケージの未対応を gs と全桁一致 (k64 モード) または独立な式との一致 (v2 モード) で埋め、K64 の誤りを 24 件新しく見つけて v2 で直した。電子以外の粒子 (`MASS`・`CHARGE`)、光子、NAFF、複素の特殊関数、成分ごとの `FREE`、`Parallelize` などが使えるようになった。

**主な修正 (v2 で結果が変わるもの。K64 の動作は k64 モードで再現)**
- 放射: tradke と tradkn の反作用の角度の符号が逆 (RADCOD の閉軌道が追跡で閉じない: B111)、放射の古典半径が質量によらない (B112)、WSPAC の最後のスライスのレンズが 8 倍強い (B127)、光子の偏光 ξ1 が ±1 を超える (B123)。
- RADINT: ANGLE=0 の K0 の極 (B115)、K1≠0 の BEND の勾配の y 成分 (B107)、軌道がずれた四極 (B126)。tbgrote の (5,6) (B109)、tradke の b (B106)。
- 光学・追跡: ずれた BEND の CALC の出口の符号 (B114)、K1≠0 の BEND の追跡の端 (B117)、ANGLE=0 の BEND の端の角 (B108)。
- 偏極: `spdepol` を側帯の項つきの独立な式に置き換え (B125。K64 は約 1000 倍大きい)、逆向きの BEND の緩和率 (B116)、傾いたソレノイドの Emittance[] の POL の量 (B130)、領域内の CAVI のスピン (B131)。
- パッケージ: CorrectOrbit の Micado が選んだ補正量だけを使う (B129)。`Union` の `SameTest`・`WriteBeamLine` の `Name` (B121・B122)。複素の特殊関数は gs が誤る点でも真値 (B118〜B120。k64 も真値)。

**新しく使える機能** (使える範囲・条件は docs/47 の v2.1.0 の (b))
- `TrackParticles`・Emittance[]・CALC6D に `MASS`・`CHARGE` (陽子・イオン)、光子 (PhotonList) の BEND・QUAD・MULT・SEXT・OCT・ドリフト・BEAMBEAM + BEND、`NPARA=n;` (gs の子プロセスの乱数の分け方の再現)。
- NAFF (`pysad/naff.py`: チューンの抽出)、複素の `Zeta`・`HurwitzZeta`・`PolyLog`・`HypergeometricU`、成分ごとの `FREE`、`$Context` の代入・`WithContext`・入れ子の文脈、`Symbol`・`Message`、VISIT/BYE、`Parallelize` (POSIX は fork で並列)。
- Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE と ANGLE=0 の BEND・MULT の F1、光学・Emittance への空間電荷 (WSPAC)、CALC6D・RADCOD つきの CALC・INTRA と SOL 領域、領域内の四極の DX・DY と幾何の fit、POL + TRPT + 加速 (SOL 領域の外)、PSPAC の粗いメッシュの警告、SynchroBeta の収束の警告。
- MAIN の要素定義の式 (`ANGLE=PI/8`)、Python 3.14 (Tk 9.0) での GUI、`pysad.nbcache` (numba のキャッシュの無効化)。

**pysad 自身の不具合の修正** (モードによらない): `Twiss["DETR"]` が INS の CALC の後で 0、`TrackParticles` が `MASS`・`CHARGE` を渡さない、BEND の F1 の放射の断片、CAVI の DX・DY・ROTATE と `V1 V20 V11 V02`、EMIT のあとの `USE`、NORFSW の Emittance[] の閉軌道、複数周の `RadiationField`、不安定なリングの POL の Emittance[] (20 分以上 → 2 秒)、成分ごとの FREE、`Element[..]+=` のフックなど ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.1.0 の (c))。

**試験の状態** (2026-11-18): 全試験 1345 通過・8 skip・失敗 0 (`PYSAD_SAD_PACKAGES` つき、k64 モードが既定。収集は 1353 件)。うち v2 の試験は 775 件 (v2.1.0 の追加は 412 件、`test_v2_*.py` 54 ファイル)。回帰 (`QUICK=1`) は 50 ベンチ FAIL 0。

**条件と限界**: v2 の値の根拠は解析式・独立な数値計算・保存則で、gs ではない。確かめたのは 0.2〜4.2 GeV の FODO リング、大きなリング A (3 GeV、周長 603 m)・PF リング・ASTRID2 (MPW 入りを含む) などで、SP8 (大きなリング B、8 GeV) は数値だけ。実機の強い結合・挿入光源のずれのあるリングでは未確認。v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 違う (tradke の符号の修正による意図した差)。

**未解決**: [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.1.0 の (d)。主なもの: `Polarization` が高エネルギー・高 RF 電圧のリングで脱偏極率を追跡の約 0.4 倍に過小評価する (軌道の非線形)、`SynchroBetaEmittance` が大きなリング A で 40 倍違う (B128。警告のみ)、v2 の Emittance[] と RADINT の矩形磁石の 1.4e-3 の差、tradke の符号と RADTAPER の重みの判断、複素の特殊関数の k64 の値 (gs の誤りを再現しない)、TrackParticles の F1 (ANGLE=0 の BEND・MULT)、SOL 領域内の光子、`DAPWIDTH` の窓、`Parallelize` の Windows の並列。

## v2.0.0 (2026-11-17)

pysad の 2 つ目の区切りで、**K64 のバグの修正と機能の追加**を行った版です。v1.0.0 (ソースに忠実な翻訳) と同じ日に出します。結果が gs と変わる箇所は、互換モードで切り替えられます。
変更の一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md)、K64 側のバグの記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (v2 で新しく見つけたものは B72〜B110)、作業の記録は `docs/history/v2_changes_*.md`。

**互換モード**
- 既定は `v2` (修正した動作)。`PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`) で v1 と同じ K64 (gs) の動作になる。gs と同じ結果が要るときは `k64` を使う。
- 既存の試験と回帰ベンチは gs との比較なので k64 モードで走る (`pysad/tests/conftest.py`、`run_regression.py`)。v2 の動作は `test_v2_*.py` が `v2` を明示して試験する。
- k64 モードは gs に合わせて 2 点だけ以前の pysad と違う: `For` の本体を実行しない (B71)、利用者のコードの `BeginPackage`・`Begin` が名前空間を作る。

**主な修正 (v2 で結果が変わるもの。K64 の動作は k64 モードで再現)**
- 放射・エミッタンス: RADINT の BEND ごとの曲率・6 極・倍精度・複数空洞 (B06・B07・B12・B13)、量子寿命とクロマティシティ (B72・B73)。
  偏向磁石の端の角のある Emittance[]・追跡の Jx (K64 は矩形磁石のリングで 1.097、正しくは 1.006。RADINT は正しい: B75)、K1 つき BEND の追跡の減衰の和 (B76)、MomentumCompaction (B77)。
  tchge のロール + ずれ (B09)、スピン (B56・B58・B64)、平衡偏極度 `Polarization` (K64 は NaN。v2 は線形のスピン軌道結合で平面リング 92.4%: B55)。
- 追跡・要素: UND の磁場 (B37)、PHSROT (B38)、SOL 領域内の BEND・MULT・CAVI (B69)、ソレノイド内のスピンの縦磁場・刻みの角 (B81〜B83)、wake の重み (B84)、BEAMBEAM の傾き (B34)、
  DELCOR (B40・B62・B99)、RADLIGHT の軌道 (B65)、tsolqur・tsolqum・tdusrn (B08・B14・B10)。
- 空間電荷: K64 の SPAC は力が 1/γ 倍に弱い (B89)。v2 は正しい係数で SPAC・PSPAC・WSPAC を追跡に入れる。
- 特殊関数: `Hypergeometric2F1`・Bessel・`ExpIntegralEi`・`HermiteH`・`LegendreP`・`InverseErf`・`PolyLog` ほか (B43〜B47・B53・B54)、`PolynomialFit` (B18)。
- スクリプト・FFS: `For` の本体 (B71)、`MessageName` の関数形 (B70)、MAIN の `OFF RFSW` (B30)、`ORG` の負数 (B26)、複素の `FindRoot` (B61)。

**新しく使える機能** (使える範囲・条件は docs/47 の (c))
- UND の磁場 (`FBX`・`FBY`)、PHSROT、SOL 領域内の BEND・MULT・CAVI・APERT、ずれて傾いたソレノイド (`DX DY DZ CHI1..3`、`F1`)、領域内の空洞の TRPT の加速・wake・AUTOPHI。
- スピン追跡 (POL + RAD) の SOL 領域・UND・PHSROT、`RADTAPER`、`LOSSMAP`、ResponseMatrix の TRPT、MULT のスキュー成分。
- MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1`、TRCOD FIXCOD / MVCOD / STNCOD / ALIGN、ORBITJITTER、PUTMEASURE / GETMEASURE (K64 のソースの読みで実装。gs では動かず比較できていない)。
- 空間電荷 SPAC・PSPAC (SPCH)・WSPAC・SELFCOD (MULT の `RADIUS` を含む)、旧式の `WAKE` コマンド、RADLIGHT の 3 次元の軌道。
- 文脈 (`BeginPackage`・`Begin`)、`Fork`・`Shared`・`Wait` (POSIX)、`System` (既定は実行しない。`PYSAD_ALLOW_SYSTEM=1` で許可)、複素引数の Bessel・超幾何・Airy。
- 走る例: [`example/v2_B/`](example/v2_B/README.md)、[`example/v2_C/`](example/v2_C/README.md)。

**pysad v1 の不具合の修正** (モードによらない): `LegendreP` の奇数 m の符号、MULT のスキュー成分、`CalculateOptics` の TRPT、加速区間の分散、EMITDIV など ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の (d))。

**試験の状態** (2026-11-17 20:30 に実行): 全試験 937 通過・4 skip・失敗 0 (941 件収集。`PYSAD_SAD_PACKAGES` つき、k64 モードが既定)。うち v2 の試験は 363 件 (`test_v2_*.py` 29 ファイル)。回帰 (`QUICK=1`) は 50 ベンチ FAIL 0。

**条件と限界**: v2 の値の根拠は解析式・独立な数値計算・保存則で、gs ではない (gs が K64 のバグを持つため)。確かめたのは 3 GeV 前後の FODO リングなど小さな系で、実機の大きなリング・挿入光源・強い結合では未確認。
MEASURE・TRCOD は gs との比較が無い。

**未解決**: [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の (e)。主なもの: `Twiss["DETR"]` が CALC の後で常に 0 (pysad の不具合。v2.1.0 で修正)、Emittance[] の POL の量の SOL 領域内 (K64 の式のまま)、空間電荷の光学・Emittance への組み込みと k64 モードの扱い、
`TrackParticles` が `MASS`・`CHARGE` を渡さない (電子以外が追えない)、残る gs との差 (BEND の `F1 FRINGE=3` + RAD の追跡、L≠0 の CAVI + `DX`・`DY`、EMIT のあとに `USE` したときの追跡)。

## v1.0.0 (2026-11-17)

PyPI: [`pysad-k64` 1.0.0](https://pypi.org/project/pysad-k64/1.0.0/) (2026-11-17 公開。`pip install pysad-k64==1.0.0`)。GitHub のタグは `v1.0.0`、PyPI 用のビルドは `release/pypi-v1.0.0` ブランチ。パッケージ名は、PyPI にある異常検知ライブラリ `pysad` と別物のため `pysad-k64` です (import 名は `pysad`)。

pysad の最初の区切りで、**ソースに忠実な翻訳**の版です。SAD (K64 版) の Python への直訳で、SAD スクリプトを gs と同じ引数で実行できます (`py -3.12 run_sad.py lattice.sad script.sad`)。

**この版に含まれるもの**
- `pysad/` 本体: 67 モジュール。SADScript の評価器、FFS のコマンド層、光学 (CELL/INS、CALC)、マッチング (FIT)、エミッタンス (EMIT、`Emittance[]`)、粒子追跡 (`TrackParticles`)、ダイナミックアパーチャ、WAKE・BEAMBEAM・CSR・POL・SOL・RADLIGHT、特殊関数・Fit、SAD の `Packages/*.n` の読み込み、Tk の GUI・描画。
- 試験: `pysad/tests/` に 578 件 (575 通過・3 skip)、回帰ベンチ 50 本 (FAIL 0)。試験が走らせる SAD 台本と、その gs・pysad の画面表示は `pysad/tests/scripts/` ([`OUTPUTS.md`](pysad/tests/OUTPUTS.md))。
- 例集: [`example/`](example/README.md) (SAD の機能・関数ごとの短い例。台本、gs の画面表示、pysad の画面表示、説明の 4 点)、PF・ASTRID2 のラティスでの例 ([`example/pf`](example/pf/README.md)、[`example/astrid`](example/astrid/README.md))。
- 文書: [`docs/`](docs/01_overview.md) に解説文書 46 本 (番号 01〜46)、[`docs/history/`](docs/history/00_RESUME.md) に作業履歴。

**基準と条件**
- 検証の基準は SAD (K64 版) の gs の出力です。K64 のソースは GitHub `KatsOide/SAD` の `master` (2026-02-18) と同じ内容です。
- ヘルプの参照は CERN の SADHelp (1.4.0.0k64) です。
- 動作確認は、基本的に自分で書いた SAD スクリプトと、サンプルのラティス、PF・ASTRID2 のラティスに対して行っています。

**未対応・未確認 (主なもの)**。最新の一覧は [`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md)。
- 未翻訳: MEASURE / TRCOD、`LOSSMAP`、`SPCH` (空間電荷)、`RADTAPER`、SOL 領域内の MULT・BEND・CAVI、`BeginPackage` / `Begin` の文脈、`System`・`Fork`。
- 未解決の差: BEND の `F1 FRINGE=3` + RAD の追跡、L≠0 の CAVI + `DX`・`DY`、EMIT のあとに `USE` したときの追跡、`RADINT` と EMIT の Jx の 2.4% の差。
- gs 側の問題で比べられないもの: `For` の本体、一部の特殊関数の値、`RadiationField` などの SIGSEGV。
- 表示や数字の桁の微妙な差は追っていません。基準は、完走することと、意味のある値が合うことです。

**この版を作った過程**: SAD の Fortran / C の物理計算の直訳、SADScript の評価器、GUI と描画を順に進め、例集と文書を gs の実動作に合わせて整えました。作業の経緯は [`docs/history/`](docs/history/00_RESUME.md)。

## 次の版 (v2 系の続き)

v2.1.0 の続きは `main` で進める。v1 系の保守 (翻訳の漏れ・誤りの修正) が必要になったら、`v1.0.0` から `v1` ブランチを切る。

- 今後の候補 ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の v2.1.0 の (d) の要判断・未解決と、[`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md) の残り):
  - 要判断: tradke の符号 (B111) と RADTAPER の追跡の重み (K64 の作者の意図の裏づけが無い)、複素の特殊関数の k64 モードの値 (gs の誤った値を再現するか)、`Polarization` の追跡との差 (軌道の非線形をどう式に入れるか)、`SynchroBetaEmittance` の大きなリング A の不一致 (B128)、v2 の Emittance[] と RADINT の矩形磁石の差、DELCOR の規格化の意図 (B99)、MEASURE・TRCOD と gs の比較 (gs が動かないので別の参照が要る)、CSROYMatrix の gs の wrongtype を k64 で再現するか。
  - 未対応: TrackParticles の F1 (ANGLE=0 の BEND・MULT)、MULT の AE1・AE2・`FB1/FB2` の `tblfre`、SOL 領域内の光子・領域内の四極のずれの POL の追跡、POL + TRPT + 加速の SOL 領域内・RF つき MULT・TCAVI、`DAPWIDTH` の窓、MAIN の TRACK の FOURIER (NAFF で代用)、`Parallelize` の Windows の並列、UND の端の磁極の整形、BOUND=0 の SOL の CALC。
  - 残る gs との差と未確認 (docs/46): gs 側の問題で比べられないもの (B49・B51・B60・B68・B133 ほか)、DETR ≥ 1 の INS と MEASURE、SP8 (大きなリング B、8 GeV) の確認、K3 が強い MULT リングの CALC の閉軌道、PSPAC の粗いメッシュの gs との差。
- v1 との差は、v2 のリリースノートに「v1 と結果が変わるもの」として一覧にする (v2.0.0 は上の節と docs/47 の v2.0.0 の節、v2.1.0 は docs/47 の v2.1.0 の節)。
