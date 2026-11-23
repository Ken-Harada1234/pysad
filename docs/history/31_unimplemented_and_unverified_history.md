# 31 未対応・未確認の一覧: 取り除いた経過の語り (元の文脈つき)

旧 46 の §6「実装・解消して一覧から外したもの」の全文 (変化の語りを含む)。docs/31 では現在形に直している。現在の状態は docs/31 の §6 にある。

## 6. 実装・解消して一覧から外したもの (現状)

記録 (`history/v2_changes_*.md`・`history/v2next_*.md`) と試験で確かめた。「両」はモードによらない。

| 項目 | モード | 現状 |
|---|---|---|
| `LOSSMAP` | 両 | `TrackParticles` の戻りに 2 行 (落ちた要素・落ちた周) を足す。gs の例と一致 |
| `RADTAPER` | 両 | Emittance[]・EMIT・追跡・CALC6D に入る。gs (CODCONV=1e-4) と x 6 桁、U0 10 桁、チューン全桁。v2 は追跡の閉軌道が 1 周で閉じる (B111) |
| MULT の `SK0 SK1 SK2 …` (スキュー成分) | 両 | 以前は黙って無視していた (pysad の不具合)。gs と追跡 2e-13、Twiss 1e-9 |
| `ResponseMatrix` の TRPT、`CalculateOptics` の TRPT の加速、加速区間の分散 | 両 | gs と全桁 |
| `NOKEEPEXP`、`RENUMBER` の INDEX 属性、`LINE["INDEX",…]` | 両 | gs の実測と一致。成分ごとの `FREE` で `tfvcomp` も gs と同じ |
| `BeginPackage`・`Begin` の文脈、`$Context`・`$ContextPath` への代入、`WithContext`、入れ子 | 両 | 利用者のコードで実装 (名前空間を作る)。gs の 19 行の台本と全行一致 |
| `Shared`・`SharedSize`・`OpenShared`、`Parallelize` (POSIX) | 両 | fork した子と共有 (無名の mmap)。`SharedSize[{1,2,3}]` = 32 (gs と同じ)。Parallelize は fork で並列 (Windows は逐次) |
| 複素引数の Bessel・超幾何・Airy | 両 | 実数のパラメータで実装 (1e-12〜1e-13) |
| 複素引数の `PolyLog`・`HypergeometricU`・`Zeta`・`HurwitzZeta` | 両 | 実装 (`specfn_cx.py`。独立な参照と 1e-10〜1e-11。gs の誤り B118〜B120 の点でも真値) |
| PolyLog の非整数次の大きい \|z\| (z < -12・z > 60 の精度、\|z\| ≥ 150 の NaN) | 両 | `_polylog_far` で相対 1e-11 (s ≤ 7.9、\|z\| ≤ 1e12) |
| `LegendreP` の奇数 m の符号 | 両 | 以前は符号が逆 (pysad の不具合)。+.2749545416973507 |
| `EMITDIV`・`EMITDIVB`・`EMITDIVQ` | 両 | 以前は無視していた。Emittance[] に効く (MAIN の `EMITDIV = 20;` も) |
| MAIN の `NAME = 式;` の代入、`FFS;` の書き方、要素定義の式 (`ANGLE=PI/8`) | 両 | gs の動作に合わせた (代入 32 個の式を gs と 2e-14。要素定義は 37 通りの値と 19 通りの構文エラー) |
| `For` の本体 | k64 / v2 | 上の §5 |
| RADINT と EMIT の Jx の 2.4% の差 (0.9415 と 0.9651) | v2 | 原因は K64 の Emittance[]・追跡が偏向磁石の端の角で放射の行路を変えないこと (B75)。**RADINT が正しく K64 の EMIT が誤り**。v2 は Emittance[] を直し、E1=E2=0.5 で解析値 1.005751 に対し 1.00586 (`EMITDIV=20` で 1e-4)。k64 モードは 0.9651 のまま |
| `BUCKET HIGHT` (`RADINT`) と `BucketHeight` (`Emittance[]`) の 0.010045 と 0.010035 | v2 | 原因は K64 の MomentumCompaction の誤り (B77)。v2 は 0.010047 と 0.010045。k64 モードは K64 のまま |
| RADINT の量子寿命・線形クロマティシティ・BEND ごとの半径・6 極・単精度の定数・CAVI | v2 | 直した (B06・B07・B12・B13・B72・B73)。k64 モードは K64 のまま |
| DELCOR の MARK・CAVI の書き込み先、波数 −1 の成分 (B40・B62) | v2 | 直した。k64 モードは K64 のまま (CAVI は HARM・PHI に入る) |
| `PolynomialFit` の平均 0 のデータ (B18)、`InverseErf` (B47) | v2 | 直した。k64 モードは K64 の値 |
| MOMENTUM を設定しない計算の警告 (B67) | v2 | 1 度だけ警告 (値は変えない) |
| BEND の `F1=0.1 FRINGE=3` + RAD の追跡の gs との差 (相対 0.7) | 両 | 解決: tfgeo の inext・iprev の再現 (B113)。gs と 1e-6 (8 桁表示の丸め)。大きなリング A・PF・ASTRID2 も一致 |
| L≠0 の CAVI + `DX`・`DY` の差 (1e-3) | 両 | 解決: CAVI の TENT・TEXIT と V1 V20 V11 V02 を追跡に足した。gs と全桁 |
| EMIT のあとに `USE` すると追跡の RF の位相が変わる差 | 両 | 解決: `trf0` を持ち越す (K64 の仕様)。gs と全桁 |
| `TrackParticles` が `MASS`・`CHARGE` を追跡に渡さない | 両 | 解決: 渡す。陽子・イオンの追跡・Emittance[]・CALC6D・RF・ソレノイドが gs と全桁一致。放射の古典半径は v2 だけ質量に比例 (B112) |
| 乱数列 (FLUC つきの追跡) の δ の rms が gs と違う (gs 4.21e-4、pysad 4.52e-4) | 両 | 解決: 原因は gs の既定 NPARA=4 (子プロセスに分ける)。NPARA=1 で pysad (k64) は gs と全桁一致、`NPARA=n;` で分け方を再現 (FODO 200 粒子 1・10・100 周、PF リング 40 粒子) |
| 光学・Emittance への空間電荷 (WSPAC) | Emittance[]・EMIT は両、CALC・CALC6D は v2 | 実装 (反復して空間電荷つきの一周の写像から)。k64 は gs と全桁一致 (K64 の最後のスライスの 8 倍のレンズ B127 を再現)、v2 は一次のラスレットの積分と 0.2〜1.8% |
| Emittance[]・EMIT・CALC6D の BEND の F1・FRINGE、ANGLE=0 の BEND (補正磁石)・MULT の F1 | 両 | 実装。k64 は gs と全桁 (FODO と大きなリング A・PF・ASTRID2・MPW 入り)。CALC のチューン = Emittance[] のチューン |
| 光子 (`tphrec`) の BEND・QUAD・MULT・SEXT・OCT・磁場つきドリフト・BEAMBEAM + BEND | 両 | 実装。k64 は gs と光子 358・325・71 個が全列一致 (E 列は 1.4e-8 違うのは B132 の仕様) |
| tdgeo の領域内の四極の DX・DY の解析微分と幾何の fit | 両 (領域内の四極のずれ) / v2 (解析微分) | 実装 (幾何が gs と全桁、fit が 2 回の Newton で収束) |
| CALC6D・RADCOD つきの CALC・Touschek (INTRA) と SOL 領域 | 両 | 実装。k64 は gs と 1e-9 |
| 傾いてずれたソレノイドの Emittance[] の POL の量 | v2 | 修正 (B130)。SpinTune が追跡と 1e-13〜5e-7 |
| POL + TRPT + 加速 (SOL 領域の外の CAVI) | v2 | 実装 (断熱加速の解析式と 5e-7) |
| ステアリング・MULT の K0・端の角 (B108) | v2 | 修正 (厳密な幾何と 3〜4%) |
| NORFSW; RAD の Emittance[] で閉軌道がずれたときの gs との差 (1.8〜4%) | 両 | 解決: pysad の不具合 (z を δ で閉じていた)。SP8 (大きなリング B、QF5 に DX = 5e-5) で gs 7.1529048e-9、k64 7.1529045e-9 |
| 不安定なリングの POL の Emittance[] が終わらない (SP8 (大きなリング B) で 20 分以上) | 両 | 解決: NaN で打ち切る (2 秒、`depol_error` つき) |
| `emit6d._tquade` の解析ヤコビアンが δ の微分を落とす | 両 | **再現しなかった** (Python 版・numba 版・数値微分が 1e-12 で一致)。z の行の関数単体の欠落を補った |
| 複数セグメント・複数周の `RadiationField` の失敗 | 両 | 解決 (重なった点を除く) |
| `Twiss["DETR"]` が INS の CALC の後で 0 | 両 | 解決 (gs と全桁。DETR ≥ 1 は未確認: §3) |
| `Parallelize`・複素の特殊関数・`$Context` の代入・`Symbol`・`Message`・VISIT/BYE・成分ごとの FREE | 両 | 実装・gs に合わせた (上) |
| `Optics.n` の `PeriodicR`・`CorrectLocal`・Micado | 両 | gs と比較 (Micado の誤りは B129 で v2 が修正) |
| `dapert.n`・`element.n`・`emit.n`・`beamline.n`・`functions.n`・`Misc.n` の gs との比較 | 両 | 小さな台本で行ごとに比べた。差の一部は直した (`Log2$`・`FunctionDefinition`・`Fill$String`・`Remove$Paren`・`Emittance->False`・TURNS の既定・`Element[..]+=` のフック) |
| Python 3.14 (Tk 9.0) の試験の失敗 2 件 | 両 | 解決 (`test_tk_version_and_font_system`、不安定なリングの `_spintrm`) |
| TrackParticles の F1・FB1・FB2・FRINGE (`tsteer` (ANGLE=0 の BEND)・`tmulti` (MULT: F1・F2・F1K1F/B)) | 両 | 実装 (K64 の動作の追加)。k64 は gs と 10 桁表示の全桁一致 (12 セルの FODO・ステアリング 6 通り・MULT 9 通り、大きなリング A の補正磁石 2 台)。numba の高速追跡と Python は 1e-13。独立な基準 (軸上の偏向 1e-14、放射の積分 1e-5、Maxwell 場の端の運動 1e-3〜1e-10) |
| MULT の AE1・AE2 (CALC・Emittance[]・追跡)・K1≠0・K0≠0 の ANGLE=0 BEND の AE (CALC・Emittance[]) | 両 | 実装 (これまで黙って無視)。k64 は gs と全桁一致 (5 通り + 大きなリング A。**K1≠0・K0≠0 の ANGLE=0 BEND は gs の CALC が落ちるので (B290)、この組み合わせは同じ値の MULT で比べた**)。v2 は B140・B141・B142・B143 を直す。大きなリング A で追跡と CALC のチューンが 2e-7 以内 |
| TrackParticles の F1・AE の放射の行路 (B108 の MULT・追跡への拡張) | v2 | 修正 (厳密な幾何 + ロビンソンの定理と ΔJx が 1.03〜1.04 倍。K64 は約 2 倍) |
| 要素ごとの FREE の何もしない変数 | 両 | gs と全桁一致 (`FREE QF.3 K1` は結合先 QF.1 を持つ何もしない変数として VAR に出る、変数の数・競合のエラー・FIX・Saved の欄・RESET) |
| `DynamicApertureSurvey` の DAPWIDTH の窓 (既定 7) | 両 | `trackd.f` のスロット管理を直訳。NORAD の地図・Score が gs と一致 (24 セル FODO 4 通り、大きなリング A 2 通り)。k64 は RAD + FLUC の最終 seed まで一致。既定を 51 から 7 に変えた (51 を前提にした台本は `DAPWIDTH = 51;`) |
| SynchroBetaEmittance の大きなリング A での 40 倍 (B128) | v2 | 原因を特定して修正 (上の §2) |
| 追跡のスピンが BEND の ROTATE で回らない (B135) | v2 | 修正 (厳密な Thomas-BMT の積分と 1e-9) |
| tradke の微分 (B138)・tradke の符号 (B111)・RADTAPER の重み | v2 | 確定・修正 (B111 の符号は放射反作用つきのローレンツ力の積分と 6 桁。RADTAPER の BEND の係数は (δ_l + δ_l+1)/2 に統一) |
| n が傾く・非平面のリングの緩和率 (B116) | v2 | 厳密式 Σ r (1 − 2/9 (n·v)²) と運動論的偏極を入れた (独立な DK の数値積分と τ が 1e-5〜2e-3) |
| 領域内の QUAD・MULT のスピン (bz ≥ 3)・領域内の光子 | v2 | QUAD・MULT のスピンは 9e-7 以下 (B146)、DRIFT・QUAD・MULT の光子を実装 (上の §1・§2) |
| PSPAC の K64 の形状関数 (k64 モード) | k64 | 実装 (gs と rms が 12 桁以上。NPARA = 1) |
| WSPAC の追跡のビーム行列 | v2 | オプション `WSPACSC = 1;` を実装 (既定はゼロ電流の行列) |
| MEASURE・TRCOD の k64 (先頭の位置)・DETR ≥ 1 の INS・複素の特殊関数 15 関数・`OYStability` の後半・Micado の変種・FITP>1・Tk 9.0 の見た目・文法の細部 (`SIN[0.1]`・`45DEG`・MAIN の `E`・`Pi`) | 両 | gs と比較・実装 (上の §2・§3 と docs/history/v2_changes_all.md の (b)) |


## 追記: 表のセルに残っていた変化の語り (現在形に直した文の旧い表現)

- 「v2 の既定は変えていない (要判断)。」(EMITDIV の既定の節。現在は v2 の既定が EMITDIV = 4 なので削除)
- 「v2 の既定は EMITDIV = 4 に変えた。」「(…v2 の既定を EMITDIV = 4 に変更)」「**対応** (既定を 4 に変更、+1.7e-3 は EMITDIV = 1 の誤差)」「v2 の EMITDIV の既定 (1)・RF が分散のある場所にある矩形リングの εx の +1.7e-3」
- 「`check_every` の既定を 1 (要素ごと) に変えた」「既定を 51 から 7 に変えた」(DAPWIDTH)
- 「`ALIGN` の対象に DECA DODECA MULT を足した」「`D_norm = …` を足している」「`compat.set_bend_k1_def('lin')` を足している」「`_step_coord_tcav_ens` を足している」「警告を足した」「CAVI の TENT・TEXIT と V1 V20 V11 V02 を追跡に足した」
- 「メンバー参照の内部の頭を `At` から `SAD$Member` に変えた」「pysad は False を読まなかったので入れた」「v2 は大きなリング A で出さなくなった」「追跡の MULT の ROTATE はスピンを回していなかった (B241、v2 で修正)」
- 「解決 (重なった点を除く)」(`RadiationField`)、「(確認。pysad のコードは変えていない。…)」、「`sxtopt.py` は他の作業で編集中のため触っていない」、「出力は変更前と完全一致」(SADScript 評価器の速度)、「SAD 自身の Packages の読み込みは従来どおり文脈を無視」、「決定済み: 現状のまま (再確認、変更なし)」
