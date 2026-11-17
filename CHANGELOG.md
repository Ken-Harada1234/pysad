# 変更履歴

バージョンは `vメジャー.マイナー.パッチ` の形です。メジャーの区切りに GitHub のリリースを作ります。

## v1.0.0 (2026-11-17)

pysad の最初の区切りです。SAD (K64 版) の Python への直訳で、SAD スクリプトを gs と同じ引数で実行できます (`py -3.12 run_sad.py lattice.sad script.sad`)。

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

## 次の版 (v2 系)

`main` で続けます。未対応・未確認の一覧 (docs/46) の項目を減らす作業が中心です。
