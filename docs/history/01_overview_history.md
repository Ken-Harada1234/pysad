# 01_overview.md から取り除いた履歴の記述

01_overview.md (現在の仕様の解説) から、版・経過に結びついた記述を取り除いた。下は取り除く前の元の行を、節の見出しつきでそのまま並べたもの (行番号は整理前のもの)。

## # 01. pysad 全体像 — SAD (K64 版) の Python 直訳

整理前の 3 行目:

````text
最終更新: 2026-11-20 (v2.3.0 の現状のスナップショット)
````

## ## 1. 目的と方針

整理前の 14 行目:

````text
- **v1 と v2**: v1 系は K64 に忠実な翻訳 (上の方針)。v2 系 (v2.0.0 から) は K64 のバグの修正・物理の誤りの修正・機能の追加を含み、互換モード [`pysad/compat.py`](../pysad/compat.py) で切り替える。既定は `v2` (修正した動作)。`PYSAD_COMPAT=k64` か `pysad.compat.set_mode('k64')` で v1 と同じ K64 (gs) の動作になる。**gs と数値を比べるときは k64 モード**にする (既存の試験と回帰は k64 で走る)。v2 の変更の一覧は [47_v2_changes.md](47_v2_changes.md)。この本文の「再導出せず直訳」「癖も再現」は k64 モードの方針で、v2 モードで変わった箇所は各文書の「v2 モード」の節に分けて書いた。
````

整理前の 16 行目:

````text
- 要素: DRIFT/BEND/QUAD/SEXT/OCT/DECA/DODECA/MULT (ANGLE≠0 を含む)/SOL/CAVI/TCAVI/UND/COORD/MAP/BEAMBEAM/MARK/MONI/APERT/INS を扱う。PHSROT は K64 で動作しない (行列が 0、B38) ため k64 モードでは追跡が NotImplementedError、v2 モードは意図された写像を実装した (docs/08 §9.7)。SPCH は K64 の tturn.f:521 の `case (icSPCH)` が PSPAC フラグのとき `tpspac` (粒子のメッシュによる空間電荷) を呼ぶ要素で、**v2 モードでは SPAC・PSPAC・WSPAC が `TrackParticles` に入る** (k64 モードは従来どおり入らない。docs/08 §9.7)。PROF は K64 のソースに現れず、parser が読むだけ。要素別の状況は §6 と各要素の文書 (23〜30)。
````

## ## 3. pysad のモジュール地図

整理前の 77 行目:

````text
| touschek.py | INTRA (IBS 付き emit)、TouschekTable / TouschekLifetime。v2.3.0 T2: `touschek_piwinski` (独立な Piwinski 型の式。v2 の `TouschekLifetime[INF, INF, nz]`・`TouschekLifetimeLMA[]` の既定。`Formula->"K64"`・k64 モードは K64 の表。docs/48 §4.6、B162) | temit.f tintrb, emit.n |
````

整理前の 87 行目:

````text
| compat.py | **互換モード** (v2.0.0 から)。`mode()`・`set_mode()`・`v2()`・`k64()`。環境変数 `PYSAD_COMPAT=k64\|v2`。docs/12 §16.6、docs/47 | — (pysad 独自) |
````

整理前の 90 行目:

````text
| spch.py | 空間電荷 (SPAC・PSPAC・WSPAC。v2 モードだけ。`TrackParticles` に入る。docs/08 §9.7。Emittance[]・EMIT・CALC6D の WSPAC は v2.1.0 で `touschek.py`・`calc6d.py` に入った: docs/07) | spkick.f, tspch_.c・spch.c・psn.c, tintrb.f twspac |
````

整理前の 92 行目:

````text
| naff.py | NAFF 風のチューン抽出 (`naff`・`naff_real`・`tune_xp`・`tunes_from_coords`。v2.1.0。放射を切った追跡の座標列からチューンを取る。MAIN の TRACK の FOURIER の代用。docs/12 §16.7、47) | K64 の直訳ではない (wtunem・wjfit に代わる独立な実装) |
````

整理前の 93 行目:

````text
| specfn_cx.py | 複素の Zeta・HurwitzZeta・PolyLog・HypergeometricU (v2.1.0。両モード。gs が誤る点でも真値。docs/33) | gamma.f の複素の部分に代わる独立な実装 |
````

整理前の 94 行目:

````text
| specfn_cx2.py | 複素の ExpIntegralE・ExpIntegralEi・PolyGamma・WhittakerM/W・ProductLog・CosIntegral・FresnelS/C・HarmonicNumber・GammaRegularized・Beta・Factorial・`Gamma[a,z]`・LegendreP・複素の次数の BesselJ (v2.2.0。両モード。複素の引数が 1 つでもあれば `specfn.install` が切り替える。docs/33) | scipy の複素対応と級数・漸化式 |
````

整理前の 95 行目:

````text
| spch_k64.py | PSPAC の K64 の形状関数 (v2.2.0。k64 モードだけ。CIC・「ほぼ最近接」・5 点差分 + 双 3 次補間・3 次スプライン。NPARA = 1 で gs と rms が 12 桁以上) | spch.c・psn.c・tspch_.c |
````

整理前の 96 行目:

````text
| trksetup.py | 追跡の設定の共通部分 (v2.3.0。`TrackSetup` = 閉軌道の求め直し・`track`・RF バケット、`find_periodicity`、`select_locations`、優先度を下げた並列。docs/48 §2) | K64 にない追加 (`track6n.py` を使う) |
````

整理前の 97 行目:

````text
| freqmap.py | 周波数マップ (v2.3.0。1 つの位置の格子から NAFF でチューンの拡散 D、表・図、共鳴線の表。台本 `FrequencyMap[]`。docs/48 §3) | K64 にない追加 (naff.py・track6n.py を使う) |
````

整理前の 98 行目:

````text
| lma.py | 局所運動量アクセプタンス (v2.3.0。位置ごとの δ±、Touschek への受け渡し `touschek_lifetime_local`。台本 `MomentumAcceptance[]`・`TouschekLifetimeLMA[]`。docs/48 §4) | K64 にない追加 (touschek.py の表を使う) |
````

整理前の 99 行目:

````text
| lmabucket.py | RF バケットの限界とアパーチャの限界の区別 (v2.3.0 L2。解析式の δ_rf、縦方向の 1 次元の写像、`limit_types`。docs/48 §4.5) | K64 にない追加 |
````

整理前の 100 行目:

````text
| fmapcheck.py | 周波数マップの信頼性の確認 (v2.3.0 S2。誤差の統計、警告 12 種、D の周回数の換算、6 極のフィードダウンの見積もり。docs/48 §3.6・§10) | K64 にない追加 |
````

整理前の 101 行目:

````text
| q1cmd.py | 周波数マップ・LMA・共鳴駆動項などの SAD の台本の関数の登録 (v2.3.0。`ffscmd.py` の入口から `install`) | K64 にない追加 |
````

整理前の 102 行目:

````text
| resdrive.py | 6 極の共鳴駆動項・ADTS・クロマティシティ (v2.3.0。厚いレンズの解析式と数値積分。台本 `ResonanceDrivingTerms[]`・`AmplitudeTuneShift[]`・`ChromaticityTerms[]`。docs/49) | K64 にない追加 (R2: 運動学・端部のフリンジの ADTS。追跡 + NAFF + フィットの `adtsdirect.py`。R3: 偏向磁石の 3 次の項を追跡の写像から取り出す `bendnl.py`。R4: 共鳴駆動項にも偏向磁石の 3 次の項 (δ の項は分散による軌道のずれ) を入れ、追跡の 1 周の写像から RDT を取り出す `rdttrack.py` で検証) |
````

整理前の 103 行目:

````text
| fmapres.py | 周波数マップの D の大きい領域を、RDT・ADTS から予想した共鳴の領域と比べる (v2.3.0 S1。誤差の与え方・統計。docs/48 §3.5) | K64 にない追加 |
````

整理前の 104 行目:

````text
| fmappred.py | 周波数マップの共鳴線の予想の作り直し (v2.3.0 S3。3 次・4 次の h と非共鳴の正規形の補正、線形の ADTS、振り子の島の幅、連続のスコア、信頼度。docs/49 §12) | K64 にない追加 |
````

整理前の 105 行目:

````text
| rdt2nd.py・nform.py・rdttrack.py・bendnl.py・adtsdirect.py | 4 次の共鳴駆動項 (6 極の 2 次)・数値の正規形・追跡の 1 周の写像からの RDT の取り出し・偏向磁石の生成関数 (追跡の写像の 2 次のテイラー係数から)・追跡 + NAFF の ADTS (v2.3.0 R2〜R6。docs/49) | K64 にない追加 |
````

整理前の 106 行目:

````text
| sxtopt.py | [試作] 6 極ファミリーの K2L の多項式モデルと SLSQP の最適化 (v2.3.0 R7。RDT・ADTS の解析式を k の線形・二次形式で組む。検証用。台本の関数なし。[`example/ddba/v07_14/sextupole/`](../example/ddba/v07_14/sextupole/)) | K64 にない追加 |
````

整理前の 108 行目:

````text
| nbcache.py | numba のキャッシュの無効化 (v2.1.0。`py -3.12 -m pysad.nbcache`。kernel を編集したら試験の前に 1 回走らせる。docs/12 §16.7) | — |
````

整理前の 128 行目:

````text
テスト: [`pysad/tests/`](../pysad/tests/) (pytest 1907 件を収集 (2026-11-21、v2.3.0 の最終)。**1899 通過・8 skip・失敗 0** (docs/47 v2.3.0)。v2.1.0 の 11/18 は 1345 通過・8 skip (v2.0.0 の 11/17 は 937 通過・4 skip)。既存の試験と回帰は gs との比較なので k64 モードで走る (`conftest.py`)。v2 の試験は `test_v2_*.py` 73 本 (1081 件。v2.1.0 の追加は 412 件) で v2 を明示する。回帰 `run_regression.py` (`QUICK=1`) は 50 ベンチ FAIL 0。`--ignore` する 3 本 [`test_bend_feeddown.py`](../pysad/tests/test_bend_feeddown.py)・[`test_coord.py`](../pysad/tests/test_coord.py)・[`test_sext_radiation.py`](../pysad/tests/test_sext_radiation.py) は sys.exit で終わる旧式のスクリプト。GUI 系 (`test_plot`・`test_tclbridge`) は全体実行で数回に 1 回失敗することがあり (単独なら通る)、環境変数 `PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>` が無いと skip。`PYSAD_SLOW=1` で遅い試験 (CSR の MakeZL、約 7 分) も走る)。試験が走らせる SAD スクリプトと gs・pysad の画面表示は [`pysad/tests/scripts/`](../pysad/tests/scripts/) (`NAME.sad`・`NAME.gs.txt`・`NAME.pysad.txt`)、試験関数との対応は [`pysad/tests/OUTPUTS.md`](../pysad/tests/OUTPUTS.md)。
````

## ### 4.1 要素・機能別の文書 (docs/22〜37、11/13〜11/16 に追加)

整理前の 149 行目:

````text
### 4.1 要素・機能別の文書 (docs/22〜37、11/13〜11/16 に追加)
````

整理前の 162 行目:

````text
| [31_delcor.md](31_delcor.md) | DELCOR (相関のあるずれ)。k64 モードは tcorr・tvcorr と乱数列まで gs と一致 (CAVI は HARM・PHI に入る K64 のまま)。v2 モードは MARK・CAVI も DX・DY、規格化を直した |
````

## ## 6. 到達点と未対応

整理前の 216 行目:

````text
- 要素・機能 (11/13〜11/16 に追加。各文書の「限界」を必ず見る): COORD・TCAVI・UND・MULT ANGLE≠0・SOL の Twiss/Emittance[]/追跡・MAP・DELCOR・GaussianCoulomb・BBBrem・BEAMBEAM・WAKE・CSR・POL のスピン (Emittance[]・EMIT・TrackParticles)・RADLIGHT・SynchroBetaEmittance・特殊関数・Fit/FindRoot・SAD パッケージ (CorrectOrbit・MakeBump ほか)・FFS の小さなコマンド・描画 (docs/22〜35)。gs との一致の桁は各文書の表にある。
````

整理前の 224 行目:

````text
| 要素 | **SOL 領域 (BZ≠0) 内の BEND・MULT・CAVI** (k64 モードの追跡は NotImplementedError。BEND は K64 の tsol.f どおりに書いたが gs が偏向をほぼ無視して合わず取り下げた。v2 モードは実装した: ローレンツ力の数値積分との一致で、gs との比較ではない。docs/44 §7a)、SOL の BZ≠0 内の四極 + RAD の追跡 (現在の gs は保存済みの参照 bench_level9 と別の値を返す。11/16 に単一粒子で切り分けた: 領域内がドリフトだけなら gs = pysad、四極が入り RAD ON かつ BZ≠0 のとき gs は四極の収束・発散力を入れていないように振る舞う (四極なしと 8 桁で一致)。gs の 2 つのビルドは同じ値で、現在の gs が SOL 領域内の四極 + RAD を誤って扱う疑いが強い。ソースのどこかは未特定、保存済みの参照を出した旧版の実行ファイルは動かず確認できない。docs/10 §6)、MAP の GEO・CompiledMap、SOL 領域内の光子 (PHOTONS の tphrec。BEND・QUAD・MULT・SEXT・OCT・ドリフト・BEAMBEAM + BEND は v2.1.0 で対応。docs/45)、TrackParticles の F1 (`tsteer`・`tmulti`)。(SOL の CALC6D・Touschek と、v2 の Emittance[] の POL の量の SOL 領域内は v2.1.0 で対応: docs/44・43) |
````

整理前の 225 行目:

````text
| コマンド・関数 | MEASURE/TRCOD は v2 モードで実装した (K64 のソースの読み。gs で動かない (B42) ので gs との比較は無い。k64 モードは未対応で無視)、DRAW/GEO 描画コマンド、FFS の `SYNCHROB` 語・`INC`・`ALIGN` (単独)・`SCALE`・`ZAP`・`DIMAD`・`BEAM` は未対応 (警告して無視。関数 `SynchroBetaEmittance` は実装済み)、`Polarization` は v2 モードで有限値 (k64 モードは K64 が NaN: B55)、旧式の MAIN の TRACK の Fourier 解析 (`wjfit/wfres/wfbin/wiord/wjdraw`。FFS から使えない)、ElementValues のフィット変数との連動 (VISIT/BYE の退避・復元は v2.1.0 で実装: docs/04)。WakeFunction は追跡も光学 (WAKEOPT、InitialOrbits) も実装済み: [23_wake.md](23_wake.md)。旧 WAKE コマンドは K64 では `WAKECOMMAND.n` (配布に無い) を呼ぶので k64 モードは同じ経路、v2 モードは `WAKE … ENDWAKE` の表を読む。複素の次数の PolyLog、Re a ≤ 0 の複素の HypergeometricU は未対応 (複素引数の PolyLog・HypergeometricU・Zeta・HurwitzZeta と光学・Emittance への WSPAC の空間電荷は v2.1.0 で実装)。一覧は [46_unimplemented_and_unverified.md](46_unimplemented_and_unverified.md) |
````

整理前の 228 行目:

````text
| 未解決の差 | K3 が強くずれの大きい MULT リングで pysad の CALC (qcod) の閉軌道が gs と別の解になり Emittance が NaN (原因未調査。`history/work_unimpl.md`)。SOL (BZ≠0) 内の四極 + RAD (現在の gs が領域内の四極の力を入れていない疑い: B68)。v2 の Emittance[] (RADCOD) は gs と 1e-4〜1e-3 違う (意図した差: B111)、v2 の Emittance[] と RADINT が矩形磁石で 1.4e-3 違う (未調査)。v2.0.0 の BEND の F1 + RAD、L≠0 の CAVI の DX・DY、EMIT のあとの USE の差は v2.1.0 で解決 (B113・CAVI の TENT・TEXIT・`trf0`)。`emit6d._tquade` の解析ヤコビアンの δ の微分の欠落は再現しなかった (46 §6) |
````

