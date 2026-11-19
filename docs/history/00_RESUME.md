# 引き継ぎメモ (pull した人と AI アシスタント向け)

最終更新: 2026-11-20 (§5a を追加。§1〜§8 の本文は 2026-11-16 のまま残し、現状は §5a と docs/46・47 を見る)

このリポジトリを clone / pull して作業を始めるときに最初に読む文書。時系列の詳細ログは置かない (作業者のローカル用メモは git に入れていない)。
経緯と照合の記録は項目ごとに `docs/history/work_*.md` (一覧は §10)、現状の説明は `docs/01_overview.md` 〜 `docs/45_radlight_guide.md` (地図は §3)、SAD の癖は `docs/11_sad_quirks.md`、SAD 自身の誤りは `docs/21_sad_bug_report.md`。
AI アシスタント (Claude Code など) は、リポジトリ直下の `CLAUDE.md` からここに来る。
未対応・未翻訳・未解決の差の最新の一覧は `docs/history/work_unimpl.md`。この文書の §6 はその要約。

## 0. 文書中の記号

文書 (特に `docs/history/`) の `<repo>` はこのリポジトリのルート、`<SAD>` は K64 版 SAD のソースツリー (`src/`、`Packages/`、`bin/gs`、`Documents/` がある所)、`<work-…>` は作業者がローカルに持つ作業フォルダ (ベンチ、デバッグ用など) で、このリポジトリには含まれない。`bench_level*` などの名前は SAD の参照値付きベンチを指すが、リポジトリには入っていない。`local/wake/` なども作業者のローカルで、リポジトリには入っていない。

## 1. これは何か

pysad は加速器設計コード SAD (K64 版、https://github.com/KatsOide/SAD) の Python への**直訳**。目的は、SAD のスクリプト (MAIN のラティス + FFS + SADScript) を
そのまま実行して、SAD と同じ数値を出すこと。一致の基準は SAD の印字精度 (1e-10〜1e-11)。物理を再導出せず、Fortran / C の式・分割数・反復順序・乱数の引き方を写す。
SAD の癖 (バグに見えるものも) は再現する (`docs/11_sad_quirks.md`)。SAD 自身のバグは再現せず `docs/21_sad_bug_report.md` (B1〜B62) に集める。

## 2. 方針 (作業者の指示。守ること)

- **翻訳の対象**: Fortran / C で書かれた物理計算は直接翻訳する。SADScript (`Packages/*.n`) は Python で書き直す (直訳ではない)。ただし GUI と、SAD の `Packages/*.n` のうち SADScript のまま動かせるもの (CorrectOrbit・MakeBump・D など) は、pysad のインタプリタでそのまま実行する (`docs/22_gui_tclbridge.md`、`docs/34_sad_packages.md`)。
  ファイル操作・リスト操作・描画は Python のライブラリ (open、list/numpy、matplotlib、tkinter) に任せ、翻訳しない (`docs/15_sadhelp_purpose_policy.md`)。
- **照合の基準**: 「物理量が SAD と有意に一致する」。1e-10 以下の雑音、表示の桁・空白、丸めカオス後の乱数 seed は追わない。検証の基準は SAD (WSL の gs) の出力。
- **SAD の参照先**: コマンドや関数の意味は CERN の SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html) か K64 版ソース同梱の `Documents/SADHelp.html`。
  KEK 公式サイトの SAD は別系統なので、ヘルプとしては参照しない (ライセンスの出典だけは KEK のページ)。数値の最終根拠は K64 版のソース (`src/*.f`、`src/*.c`、`Packages/*.n`)。
- **文書の書き方**: 解説文書 (`docs/*.md`) は履歴なしの現状のスナップショット。「最初 X だったが直した」のような経緯・日付付きの追記は `docs/history/` にだけ書く。
  未確認・限界・使った制約条件 (エネルギー、磁場、周長) は必ず書く。数値だけ引用しない。
- **版**: 成果物・html・ZIP は加算式で上書きしない。
- **計算資源**: 他の計算と同居することが多いので、長い実行は 1 本ずつ・最低優先度で。動作テストでは誤差の実現数・チューンサーベイの格子・DA の周回数を間引く。

## 3. 構成

| 場所 | 内容 |
|---|---|
| `pysad/` | 本体。モジュールの地図は `docs/01_overview.md` §3、行数と準拠状況は `docs/17_implementation_summary.md` |
| `pysad/tests/` | pytest 578 件収集 (64 ファイル)。SAD (gs) から取った参照値を埋め込み、または `data_*.npz` に持つ。ラティスファイルを読むテストあり (誤差入りラティスの 2 件は環境変数 `PYSAD_ERR_LATTICE` を指定したときだけ走る) |
| `run_sad.py` | `py -3.12 run_sad.py lattice.sad script.sad` (SAD の `gs` と同じ引数、複数ファイルは連結) |
| `docs/` | 解説 (現状)。`01_overview.md` から。API は `12_pysad_api.md` |
| `docs/history/` | この文書と `work_*.md` (項目ごとの経緯・照合の記録)、`review_notes.md`、`sad_patch/` (元の SAD の SIGSEGV を直す差分) |

`docs/` の一覧 (番号順):

| 番号 | 内容 |
|---|---|
| 01〜12 | 全体像、ラティスと MAIN、FFS、SADScript、CALC、マッチング、エミッタンス、追跡、DA、ベンチ、SAD の癖 (11) |
| 13〜15 | API (12_pysad_api)、SADHelp の索引 (13_sad_help_index)、Fortran の状況 (16_fortran_status)、SADHelp との対応 (14_pysad_vs_sadhelp)、実装の要約 (17_implementation_summary)、SADHelp の目的と方針 (15_sadhelp_purpose_policy) |
| 16〜19 | 物理ガイド、ソース解析、翻訳候補、SAD のバグ報告 (B1〜B62) |
| 20〜35 | 要素・機能別: GUI/Tcl (20)、WAKE (21)、BEAMBEAM (22)、COORD (23)、TCAVI (24)、UND (25)、MULT ANGLE≠0 (26)、SOL の光学 (27)、MAP (28)、DELCOR (29)、GaussianCoulomb・BBBrem (30)、特殊関数・Fit (31)、SAD パッケージ (32)、FFS の小コマンド (33)、POL/スピンと CSR (34)、描画 (35) |
| 36〜43 | 大項目の解説文書 (何ができるか、使い方の例、物理と K64 の式、実装の構成、K64 の癖、gs との一致、限界): `38_wake_guide.md`、`39_beambeam_guide.md`、`40_csr_guide.md`、`41_tcltk_guide.md`、`42_graphics_guide.md`、`43_pol_guide.md`、`44_sol_guide.md`、`45_radlight_guide.md` |
| `sadhelp_summary.html`、`html/` | SADHelp の要約と、文書の html 版 (`make_html.py` で作る生成物。手で編集しない) |

使わない旧モジュール: `pysad/ffs.py`, `matching.py`, `tracking.py`, `optics.py` (SAD 未照合の独自実装、非推奨)。新しいコードでは使わない。

## 4. 環境と実行

- Python 3.12 に numpy、scipy、numba、pytest。Windows では `py -3.12`。出力に日本語や記号が入るので `PYTHONIOENCODING=utf-8`。
- テスト: `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py`
  (無視する 3 本は sys.exit で終わる旧式のスクリプト)。**575 通過・3 skip** (578 件収集、2026-11-17)。**2026-11-20 (v2.2.0) は 1647 通過・9 skip・失敗 0 (1656 件収集、約 27 分)**。
  - GUI 系 (`test_plot.py`、`test_tclbridge.py`) は環境変数 `PYSAD_SAD_PACKAGES=D:/work/SAD-master/SAD-master/Packages` (K64 の `Packages` の場所) が無いと skip。評価器 (`sadscript.py` など) を変えたら、これを設定して全試験と回帰を回す。
  - 全体実行では GUI 系が数回に 1 回、窓の更新待ちで失敗することがある (単独なら通る)。
  - 試験が走らせる SAD スクリプトと gs・pysad の画面表示は `pysad/tests/scripts/`、対応表は `pysad/tests/OUTPUTS.md` (`py -3.12 pysad/tests/make_outputs_md.py` で作る。取り出し方は `docs/history/work_test_outputs.md`)。
  - `PYSAD_SLOW=1` で遅い試験 (CSR の MakeZL、約 7 分) も走る。
- 回帰 (ローカルの `run_regression.py`、リポジトリには含まない): `QUICK=1 PYTHONIOENCODING=utf-8 py -3.12 -u run_regression.py` は遅い 4 件 (7・12t・13a・16a) を SKIP して 46 件。夜に 4 件も実行し、**50 ベンチすべて実行して FAIL 0** (非 QUICK で 205 s、2026-11-16)。
- SAD 本体 (参照値を新しく取るとき): K64 版をビルドした `bin/gs` が必要 (このリポジトリには含まない)。WSL で `gs file.sad > file.out`。ラティスとスクリプトは**1 ファイルに連結**して渡す
  (`gs lat script` の 2 引数は標準入力待ちで止まる)。`NPARA = 1;` を入れる。WSL の SAD は `draw` / `TkWait` で落ちる・止まるので、描画の行は外す。lattice 台本の先頭に `FFS;` を付けない (SADScript の台本は逆に必要)。
- 参照値は多くがテストに埋め込み済みなので、SAD 本体が無くても `pytest` は通る。SAD 本体が必要なのは、新しい機能を SAD と照合するときだけ。
- 元の SAD は SIGSEGV で落ちる場合が 3 つある (wake・SOL・tapert)。直す差分は `docs/history/sad_patch/wake_sigsegv_fixes.patch` (経緯は `work_wake.md`)。

## 5. 現状 (SAD と一致しているもの、していないもの) — 2026-11-16

- 一致: CALC (4D/6D、CELL/INS、結合、オフモメンタム、不安定ラインの続行)、EMIT / Emittance[] (RFSW × RADCOD、INTRA、報告の全行)、FIT/GO (解析微分、DESCEND、fuzz、2 段目、INS、フック)、
  追跡 (NORAD/RAD/FLUC、粒子損失)、DynamicApertureSurvey (NORAD)、RADINT、幾何 (GX…、GEO。SOL 領域の tsgeo は翻訳済みで docs/29 に記述、gs との一致の度合いは docs/29 を参照)、TRPT の加速、Touschek、誤差の実現、光学関数 (TransferMatrix 等)。詳細は `docs/17_implementation_summary.md` §3。
- 試験: 575 通過・3 skip (578 件収集)、回帰は 50 ベンチすべて実行して FAIL 0 (QUICK=1 は遅い 4 件を SKIP して 46 件。§4)。
- 11/13〜11/16 に実装して gs と照合した要素・機能 (各文書に、一致の度合い・使った条件・未確認がある):
  COORD (23)、TCAVI (24。Emittance[] は L=0 のみ gs と一致、L≠0 は gs が壊れる B36)、UND (25。K64 は磁場を設定できずドリフト B37)、MULT ANGLE≠0 (26)、SOL の Twiss・Emittance[]・追跡 (27)、MAP/ExternalMap (28)、DELCOR (29)、
  GaussianCoulomb(Fitted)・BBBrem (30)、特殊関数 458 式・Fit・FindRoot (31)、SAD の Packages (CorrectOrbit・MakeBump ほか、32)、FFS の小コマンド (33)、POL/スピンと CSR (34)、描画 (ListPlot・OpticsPlot ほか、35)、GUI/Tcl (20)、WAKE (21)、BEAMBEAM (22)。
- 11/17: APERT 要素の追跡 (楕円・矩形・DX DY・ROTATE・JDPX、gs と一致。`work_apert.md`)。
- 11/17: SADScript・入出力を gs に合わせた (Module の字句スコープ、Return、Switch/Which、リストの ==、Out[n]・%、Read、StringMatchQ、エスケープ、Random・SeedRandom ほか。For の本体と BeginPackage の文脈は合わせず。`work_scriptgs.md`)。
- 11/16 の追加 (詳細は `work_session_20261116.md`): POL の Emittance・EMIT のスピン (MULT・空洞・SOL まで) と TrackParticles のスピン追跡 (gs と 4e-12 以下)、RADLIGHT の軌道 (gs と差 0)、SynchroBetaEmittance (gs と相対 1e-5〜6e-6)、
  CSRInit・CSRMatrix mode 1000・CSRConvert・MakeZL、tdgeo (DRIFT の L と BEND の ANGLE・L だけ解析微分)、Skip・Read の型、`measure.py` (tlum/tluma)、SOL の長さを周長に入れない (K64 の癖)、初期の MOMENTUM は 30 GeV。
- GUI: SAD 自身の Tkinter.n / KBFrame.n / KBMainFrame.n を SADScript として解釈し、tkinter で表示する (`pysad/tclbridge.py`、`pysad/sadclass.py`、`docs/22_gui_tclbridge.md`)。
  KBMainFrame と demo-LabelFrame のウィジェット木は SAD (WSL の gs) と一致。入れ子の部品 (KBFFrameToggled / KBFFrameTabed)、見た目と操作の一致は未確認。`PYSAD_SAD_PACKAGES` が無ければ無効。
- WAKE: 追跡 (TRPT、LWAKE / TWAKE) と光学 (WAKEOPT、InitialOrbits) を実装 (`pysad/wake.py`、`wakeopt.py`、`docs/23_wake.md`)。gs と 10 桁一致。旧式の WAKE コマンドは K64 では `WAKECOMMAND.n` (配布に無い) を呼ぶので、pysad も同じ経路 (tfwake.f は移植しない)。
- BEAMBEAM: 追跡・BSTRL・LuminosityMonitor・PHOTONS (BEAMBEAM の光子だけ) が gs と一致 (`pysad/beambeam.py`、`docs/24_beambeam.md`)。

## 5a. 現状 (2026-11-20、v2.2.0)

- **版**: v2.2.0 (`pysad.__version__`・`pyproject.toml`)。PyPI は v2.1.0 まで公開 (v2.2.0 は未公開。push・タグ・PyPI はしていない)。v2.0.0 → v2.1.0 → v2.2.0 の変更の一覧は `docs/47_v2_changes.md`、K64 側のバグは `docs/21_sad_bug_report.md` (B01〜B159。B135〜B159 が v2.2.0)、未対応・未確認は `docs/46_unimplemented_and_unverified.md`。経緯と根拠は `docs/history/v2next_*.md` (v2.2.0 は P1〜P5)、例は `example/v2_next/`。
- **試験**: `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本を `--ignore`、`PYSAD_SAD_PACKAGES` つき、k64 が既定) は **1647 通過・9 skip・失敗 0** (1656 件収集)。回帰 (`QUICK=1 PYTHONIOENCODING=utf-8 py -3.12 -u run_regression.py`) は **50 ベンチ FAIL 0** (138 s)。v2.1.0 は 1345 通過・8 skip。numba のカーネル (`elkern.py`・`perfkern*.py`・`turnkern.py`・`track6n*.py`) を編集したら試験の前に `py -3.12 -m pysad.nbcache`。
- **互換モード**: 既定 `v2`、`PYSAD_COMPAT=k64` で K64 (gs) の動作。k64 モードは K64 の動作の追加 (MEASURE・TRCOD の先頭の位置、PSPAC の K64 の形状関数、TrackParticles の F1・MULT の AE、DAPWIDTH の窓、要素ごとの FREE の何もしない変数) と pysad 自身の不具合の修正だけ変えた。新しいオプション (v2): `spin.ST_EXACT`・`KINETIC`・`STEP_RATE`・`CHROMATIC`、`compat.set_bend_k1_edge('face')`、FFS の `WSPACSC`。
- **spin 系の試験の許容 (v2.2.0 の確認)**: `test_v2_A_spin.py::test_polarization_with_orbit_errors_v2` は P1・P2 の変更の後で通る。P1 が許容を 3e-3 → 1.5e-2 に緩めていたが、P は 0.8798 → 0.8695 で、運動論的偏極 −7.1e-3・厳密な緩和率 −3.0e-3・刻みごとの脱偏極率 −2.5e-4 の内訳で説明できる (3 つを切ると 0.8798)。許容を 0.8695 ± 3e-3 に戻し、内訳の試験を足した (docs/47 の v2.2.0 の (e))。
- **残課題 (詳細は docs/47 の v2.2.0 の (d))**: 判断が要るもの (既定は担当の選択のまま): EMITDIV の既定 (1)・DAPWIDTH の既定 (7)・`WSPACSC` の既定 (0)・K1≠0 の BEND の端の模型の既定 (`radial`)・放射の β の因子 (B112)・SOL 領域内の ψ の DISFRIN・偏向磁石の本体の B∥ (B136)・複素の特殊関数の k64 の値・`CSROYMatrix` の wrongtype。未解決: 脱偏極率の追跡との比は 4.2 GeV・12 MV の 1 点、B108 の残る 2.7% の系統差、RF が分散のある場所にある矩形リングの εx の +1.7e-3、SOL 領域内の加速 (POL + TRPT)、SB のシンクロ・ベータ共鳴の追跡との比較、MARK の生の DETR ≥ 1 の INS (B159)、逆向きの要素の F1・AE の入れ替え。
- **文書・例に出してよいラティス**: `pysad/tests/kek-ls_v03_12.lat` (大きなリング A: 3 GeV、周長 603 m)、`example/pf/lattice/*` (PF)、`example/astrid/lattice/*` (ASTRID2)。ほかの大きなリングは「大きなリング B (8 GeV、1436 m)」のようにエネルギーと周長で書き、名前・パスを書かない。
- §6 以降は 2026-11-16 の時点のまま (多くは解決済み。現状の未解決は上と docs/46・47)。

## 6. 残課題 (優先順。最新の詳細は `docs/history/work_unimpl.md`)

### 6.1 見つけた未解決の差 (pysad と gs が合わない。原因は未確定のものを含む)

1. BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)。DX・DY・ROTATE だけなら 1e-13 で一致。
2. L≠0 の CAVI に DX・DY を付けた RAD ON の追跡で 1e-3 ずれる。
3. SOL (BZ≠0) 内の四極 + RAD の追跡で、現在の gs (2026-10-19 ビルド) が、保存済みの参照 (bench_level9、2026-10-14 に取得。pysad はこれと一致) と別の値を返す。**11/16 夜に絞り込み: gs は RAD ON かつ BZ≠0 の領域内で四極の力を入れていない疑いが強く (四極を外した場合と 8 桁一致、BZ=0.01 と 0.2 でほぼ同じ値で BZ に連続でない)、gs の疑い**。ただし K64 ソースのどこかは未特定 (落とし穴 §8、`work_session_20261116.md` §11.3)。
3b. (11/17 に解決) Tkinter の Widget 値: `e=Entry[w]`・`w=Window[]` が gs と同じ `Widget[…]` になり、メソッド (Insert・Get・Delete) が働く。原因は `With` が入れ子の規則の呼び出し元の引数名を見ていたこと (`work_example_fixes_20261117.md`)。残る差は最初の `Window[]` が Toplevel でなく `.` 自身であること。
4. SOL 領域内の BEND・MULT・CAVI は未対応 (NotImplementedError)。BEND は K64 の tsol.f どおりに書いたが、gs が偏向をほぼ無視して合わず取り下げた。原因未確認。
5. K3 が強くずれの大きい MULT リング (`local/wake/pm1.sad`: K3=300、ずれ 1e-4) で、pysad の CALC (qcod) の閉軌道が gs と別の解 (y=-0.033) になり Emittance が NaN。gs は小さい軌道で通る。qcod の初期値・後退の違いと疑うが未確認。
6. `emit6d._tquade` (Python 版) の解析ヤコビアンが、x がずれた軌道で δ の微分を落とす (DX/DY のある QUAD を含むリングで T の δ 列が 1e-5、ε_y で 1e-4)。`emit_extra`・`calc6d._turn`・`touschek` が使う。`emit6d.tturne` (elkern) は 1e-12 で一致。SOL つきの Emittance[] を emit_extra に入れられないのもこのため (`work_unimpl.md` 後で直す一覧 1・2)。
7. 上以外の未対応・未確認: DELCOR の CAVI・OCT、SOL 領域の CALC (数値微分の 1e-10 差)、GEO のない SOL の S 座標の途切れ (K64 が周長を狂わせる。pysad は再現しない)、ExternalMap["GEO"]・CompiledMap。

### 6.2 未翻訳 (理由つき)

- SPCH 要素 (`tturn.f:521` が PSPAC のとき `tpspac` を呼ぶ。pysad は読むだけで物理なし。PSPAC・WSPAC・SPAC フラグも無効)、LOSSMAP。
- 結合リングの SynchroBetaEmittance: gs が SIGSEGV で比較不能、pysad は `ValueError: cannot convert float NaN to integer` で失敗 (エラーの出し方が粗い)。RadiationField の複数区間 (複数周) の軌道は区間の始めの点が重なり `x must be strictly increasing sequence` で失敗 (gs は RadiationField が SIGSEGV)。
- MEASURE / TRCOD (`tracke`・`ttinit`・`ttstat`・`trcoda`・`trackb`): gs が動かない (B42) うえ、使うスクリプトも見つからない。
- `wjfit` / `wfres` / `wfbin` / `wiord` / `wjdraw`: MAIN の TRACK の Fourier 解析で、FFS から使えない。
- SOL 領域内の MULT / BEND / CAVI と、`tspac.f` の追跡側。tdgeo の四極などの変数 (`qdtrans`・`tftmatu` が要る。中心差分のまま)。
- `Polarization` (K64 が NaN: B55)、PHSROT (K64 で動作しない: B38)。
- MessageName の関数形 `MessageName[f,"usage"]` (gs は "Undefined MessageName" を出す。見送り)。
- PolyLog の非整数次数で |z|≥150、複素引数の Bessel・超幾何。ResponseMatrix の TRPT (加速) 経路。
- BEAMBEAM 以外の要素の光子 (tphrec。PHOTONS が gs と乱数列でずれる)。BeamPlot の描画の gs との比較、CanvasSymbol、Canvas3D*。

### 6.3 古い課題 (変化なし)

- C6: 補正量 (SAD の teigen が縮退固有値で丸めに依存する)。
- fit の多解・FITP>1 (SAD が再起動するので参照値が取れない)、不安定域の閉軌道探索の初期値の差は、SAD 側の仕様なので記録のみ。
- (11/16 夜に解決したので課題から外した: B33 の出力座標の差 = セッション最初の l_track=0、UNIPOL・RELW、DUMP の書式、MATRIX コマンド、LDSolve の複素行列。`work_session_20261116.md` §11.2)
- GUI の残り: 最初の Window の親 (gs は Toplevel)、入れ子の部品 (demo-KBFrame1)、Canvas / TabFrame の描画、見た目と操作の比較 (`docs/22_gui_tclbridge.md`、`work_gui_tclbridge.md`)。
- 未確認のまま: SynchroBetaEmittance の実用の大きなリング・不安定なチューン・結合・RFSW の組み合わせ。CSR の ScanTrack・OYStability などの長い流れ。tdgeo の動く例の試験。

## 7. 新しい機能を足す・直すときの手順

1. 新しく書く前に、該当の `docs/NN_*.md` を読み、既存のコードをコピーして直す。
2. SAD のソース (`src/*.f`、`src/*.c`、`Packages/*.n`) で、実際に通るルーチンを特定する (itfunaloc、`word=='…'`、`REGn("Name")` で登録を辿る。ビルドに使われるファイルは `mk/sad.obj.mk` で確かめる。例: BEAMBEAM は `bbstrhl.f` で `bb.f` は未使用)。
3. 式・分割数・反復順序・乱数の引き方を Python に写す。数値微分や厳密な式への置き換えをしない (SAD の解析行列は厳密微分と 2e-8 ずれる)。
4. 小さなラティス (7 セル + RFC 程度) と短い SAD スクリプトを書き、SAD の `gs` で実行して pysad の出力と照合する。差が出たら SAD のソースに戻る。gs が SIGSEGV で落ちる・何も出さない場合 (B35、B42) は原因を特定して記録する。
5. 照合したケースを `pysad/tests/` に参照値付きで足す。GUI 系の試験は `PYSAD_SAD_PACKAGES` を設定して回す。
6. 文書: 解説文書の主題の節に現状として書く (未確認と限界、使った条件も)。SAD の癖は `docs/11_sad_quirks.md` に、SAD の誤りは `docs/21_sad_bug_report.md` に番号付きで足す。経緯と照合の記録は `docs/history/work_*.md` に書く。
7. 未実装のものを黙ってドリフトや 0 にしない。`SADError` / `NotImplementedError` で拒否する (過去に COORD・TCAVI・UND・BEAMBEAM は黙ってドリフトに、SOL は無視されて結合 0 になっていた)。

## 8. 落とし穴

- SAD の FFS の `ON` / `OFF` は FFS 内では無効 (`INTRA;` / `NOINTRA;` と書く)。EMIT コマンドは GAUSS を立てる。SADScript の `1GEV` は SAD では `1*GEV` (=1)、pysad は単位として読む (意図的な拡張)。
- FFS[] の中のコマンドがエラーを出すと最上位の文 (Do など) が打ち切られる。LINE[…]= の未定義要素も同じ (`docs/11_sad_quirks.md` C25)。
- qcod が溢れても SAD は止まらず壊れた値で続行する。pysad は直前の値か NaN (C26)。
- 乱数: 既定の "SAD" 乱数 (合同法 + Box–Muller) を int64 で写してある。FLUC は要素外側・粒子内側の順に引く。損失周の最終 seed はカオスで一致しない。
- INTRA ON の EMIT 報告は 2 段 (零電流 → Touschek の表 → IBS 込み)、PBUNCH の既定は 1e10 (A40, A41)。
- 4D の CALC では空洞はドリフト、TRPT かつ RFSW のときだけ加速する。TRPT の放射は ndivrad が 10 倍 (tffs.f。四極の分割 gs 4)。
- 数値の丸めに強く依存する比較 (ゼロ同然の値、1e-17 など) は一致を求めない。
- **gs のビルドと保存済み参照の食い違い (2026-11-16)**: 現在の `bin/gs` は 2026-10-19 ビルドで、10-14 に保存した参照 (bench_level9 など) と同じ入力で別の値を返すことがある (SOL BZ≠0 内の四極 + RAD の追跡: S R T1 P1 の x が -.0027867898 から -.0014539780)。pysad は保存した参照と一致し、現在の gs とは合わない。保存した参照と現在の gs のどちらが正しいかは、11/16 夜の調査で現在の gs の疑いが強い (SOL 領域内の四極の力が入らない。未確定)。参照を新しく取るときは gs のビルド日付を記録し、保存済みの参照と食い違ったら、NORAD など差の出ない条件で両方が一致するかをまず見る。
- 初期の MOMENTUM は 30 GeV (K64 initb1.f)。以前の文書は 1 GeV と書いていた。エネルギーを書かないスクリプトでは効く。
- SOL (BOUND≠0) の長さは K64 では周長 (周回周波数・空洞の位相・MomentumCompaction・偏極時間) に入らない。`emit6d._emit_circ` と追跡の `circ` が同じ扱いにしてある。
- GEO のない SOL では K64 の S 座標が途切れ、周長が狂う (MomentumCompaction が 100 倍、偏極時間が約半分)。pysad は GEO を見ないので影響を受けず、再現しない。
- クラスのメンバーでない名前の `obj@name[...]` は式のまま (gs に合わせた)。`MessageName` の `f::usage` の形は一致、関数形は不一致。
- RADLIGHT の軌道: K64 は BEND の中を記録せず、横方向の座標も入れない。`example.oide.sad` の流れ (B を ANGLE 0 の STEER にする) でないと意味のある場にならない (`docs/45_radlight_guide.md`)。
- **gs は開けない OpenRead で入力行を打ち切る**: ファイルが無いと gs はメッセージを出して、その行の残りを捨てる (`Read[h,…]`・`Close[h]` は wrongtype)。pysad も同じにしてある (FileNotFoundError で落とさない)。
- **l_track=0 の最初の tapert**: K64 の tapert は落とす印に `kptbl(i,4)=l_track` を書く。l_track はセッション最初の TrackParticles では 0 (CALC をしても変わらない) なので、最初の要素では粒子が落ちず (座標は 3·LOSSAMPL に頭打ち)、後の tapert で落ちる。2 回目の TrackParticles からは最初の要素で落ちる。pysad は `TrackRun.l_track0` で再現。同じセッションで比べるときは、gs の何回目の呼び出しかを揃える (B33)。
- RAD=1 は「その要素は放射しない」(K64 の仕様、B29 はバグではない: `krad = (RAD_QUAD==0) .and. al/=0`)。MON・MCOD・TCOD のコマンドは K64 に無い。
- 元の gs は、SOL 領域始点の wake、光学の片方だけの表、TRPT で最初から LOSSAMPL の外の粒子、で SIGSEGV (`sad_patch`)。「SOL 領域内で粒子が落ちると SIGSEGV」は誤りだった (LOSSAMPL が初期振幅より小さかっただけ)。
- TrackParticles の第 2 引数は終点の要素番号 (周回数ではない)。呼び出し形式: `TrackParticles[{loc, {{x..},{px..},{y..},{py..},{z..},{dp..},{flg..}}}, 終点, 始点, 終点]`。POL のスピン追跡では 9〜11 列。
- 文字化け: WSL の heredoc 内の python は `\` などが壊れる。ファイルに書く。gs に `Total` は無い (`Plus@@`)。
- 生成物の html (`docs/html/`) は親が `make_html.py` で作り直す。手で編集しない。

## 9. リポジトリの扱い

- 追跡するのは `pysad/`、`docs/`、`run_sad.py`、`README.md`、`LICENSE`、`CLAUDE.md`、`.gitignore`。ベンチ、回帰スクリプト、出力、`local/` は含めない。
- ライセンスは SAD 同梱の LICENSE と同じ文面 (Copyright (c) 1986-2016, KEK Accelerator Laboratory)。pysad は SAD の派生物としてこの条件に従う。

## 10. 履歴文書 (`docs/history/`) の索引

- v2 系: `v2_changes_*.md` (v2.0.0)、`v2next_A1`・`A2`・`B1`・`B2`・`B3`・`C`・`D1`・`D2`・`E` (v2.1.0)、`v2next_P1`〜`P5` (v2.2.0)。
- 最新: `work_example_fixes_20261117.md` (example 10〜13・PF T5 の gs との差の修正、11/17)、`work_unimpl.md` (未実装要素の実装と未解決の差、11/14〜)、`work_session_20261116.md` (11/16 の作業を話題別に)。
- 要素別 (11/14〜): `work_wake.md`、`work_beambeam.md`、`work_gui_tclbridge.md`。
- 10 月〜11/06 の作業 (ほぼ確定。各ファイル末尾に 2026-11-16 の確認がある): 解析微分 (`work_analytic_deriv`)、CALC6D (`work_calc6d`)、BEND K0 (`work_bend_k0`)、DESCEND (`work_descend`)、emit (`work_emit_extra`・`work_emit_report`・`work_emit_unstable`)、
  幾何 (`work_geo`・`work_geo_sig`)、INS (`work_ins_offset`)、統合 (`work_integration`)、item4・item6、K0 fit (`work_k0fit`)、損失の乱数 (`work_loss_rng`・`work_loss_rng2`)、MAIN の TRACK (`work_main_track`)、match (`work_match_fixes`)、性能 (`work_performance`・`work_perf2`・`work_perf3`)、
  ROTATE (`work_rotate`)、小項目 (`work_small_items`)、Touschek (`work_touschek`)、追跡の拡張 (`work_track6n_ext`)、TRPT の加速 (`work_trpt_acc`)、不安定ライン (`work_unstable`)、z/dp (`work_zdp`)、レビュー (`review_notes`)。
