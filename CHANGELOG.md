# 変更履歴

バージョンは `vメジャー.マイナー.パッチ` の形です。メジャーの区切りに GitHub のリリースを作ります。

**バージョンの方針**
- **v1 系: ソースに忠実な翻訳。** SAD (K64 版) のソースの動作を、そのまま Python に写すことを目的とします。K64 のバグや癖も、gs の実動作に合わせて再現します (癖は [`docs/11_sad_quirks.md`](docs/11_sad_quirks.md)、K64 側の不具合の記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md))。v1.x のマイナー・パッチは、翻訳の漏れの追加と、翻訳の誤りの修正です。
- **v2 系: ソースのバグの修正と機能の追加。** K64 のバグを直し、ソースにない機能を加えた版です。gs と同じ結果にならない箇所が出てきます。v1 との違いは、リリースごとにここに書きます。

## v2.0.0 (2026-11-17)

pysad の 2 つ目の区切りで、**K64 のバグの修正と機能の追加**を行った版です。v1.0.0 (ソースに忠実な翻訳) と同じ日に出します。結果が gs と変わる箇所は、互換モードで切り替えられます。
変更の一覧は [`docs/47_v2_changes.md`](docs/47_v2_changes.md)、K64 側のバグの記録は [`docs/21_sad_bug_report.md`](docs/21_sad_bug_report.md) (v2 で新しく見つけたものは B72〜B110)、作業の記録は `docs/history/v2_changes_*.md`。

**互換モード**
- 既定は `v2` (修正した動作)。`PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`) で v1 と同じ K64 (gs) の動作になる。gs と同じ結果が要るときは `k64` を使う。
- 既存の試験と回帰ベンチは gs との比較なので k64 モードで走る (`pysad/tests/conftest.py`、`run_regression.py`)。v2 の動作は `test_v2_*.py` が `v2` を明示して試験する。
- k64 モードは gs に合わせて 2 点だけ以前の pysad と違う: `For` の本体を実行しない (B71)、利用者のコードの `BeginPackage`・`Begin` が名前空間を作る。

**主な修正 (v2 で結果が変わるもの。K64 の動作は k64 モードで再現)**
- 放射・エミッタンス: RADINT の BEND ごとの曲率・六極・倍精度・複数空洞 (B06・B07・B12・B13)、量子寿命とクロマティシティ (B72・B73)。
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

**未解決**: [`docs/47_v2_changes.md`](docs/47_v2_changes.md) の (e)。主なもの: `Twiss["DETR"]` が CALC の後で常に 0 (pysad の不具合)、Emittance[] の POL の量の SOL 領域内 (K64 の式のまま)、空間電荷の光学・Emittance への組み込みと k64 モードの扱い、
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

v2.0.0 の続きは `main` で進める。v1 系の保守 (翻訳の漏れ・誤りの修正) が必要になったら、`v1.0.0` から `v1` ブランチを切る。

- 今後の候補 ([`docs/47_v2_changes.md`](docs/47_v2_changes.md) の (e) の要判断と、[`docs/46_unimplemented_and_unverified.md`](docs/46_unimplemented_and_unverified.md) の残り):
  - 要判断: 空間電荷の k64 モードの扱い、光学・Emittance への WSPAC、`TrackParticles` への `MASS`・`CHARGE` の受け渡し、K64 の spdepol の導出、`Twiss["DETR"]` の修正、
    疑いのまま残した K64 の項目 (B106〜B109) の確認、MEASURE・TRCOD と gs の比較 (gs が動かないので別の参照が要る)。
  - 未対応: スピン追跡の経路での TRPT の加速、Emittance[] の POL の量の SOL 領域内、UND の端の磁極の整形、CALC (4 次元) の RADTAPER、複素引数の PolyLog・HypergeometricU・Zeta、`Parallelize` の並列化。
  - 残る gs との差と未確認 (docs/46): BEND の `F1 FRINGE=3` + RAD の追跡、L≠0 の CAVI + `DX`・`DY`、EMIT のあとに `USE` したときの追跡、gs 側の問題で比べられないもの (B49・B51・B60・B68 ほか)。
- v1 との差は、v2 のリリースノートに「v1 と結果が変わるもの」として一覧にする (v2.0.0 は上の節と docs/47)。
