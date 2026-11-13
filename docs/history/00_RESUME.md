# 引き継ぎメモ (pull した人と AI アシスタント向け)

このリポジトリを clone / pull して作業を始めるときに最初に読む文書。時系列の詳細ログは置かない (作業者のローカル用メモは git に入れていない)。
経緯と照合の記録は項目ごとに `docs/history/work_*.md`、現状の説明は `docs/01_overview.md` 〜 `docs/19_sad_bug_report.md`、SAD の癖は `docs/11_sad_quirks.md`。
AI アシスタント (Claude Code など) は、リポジトリ直下の `CLAUDE.md` からここに来る。

## 0. 文書中の記号

文書 (特に `docs/history/`) の `<repo>` はこのリポジトリのルート、`<SAD>` は K64 版 SAD のソースツリー (`src/`、`Packages/`、`bin/gs`、`Documents/` がある所)、`<work-…>` は作業者がローカルに持つ作業フォルダ (esrf/pf2 のスクリプトの写し、ベンチ、デバッグ用など) で、このリポジトリには含まれない。`bench_level*` などの名前は SAD の参照値付きベンチを指すが、リポジトリには入っていない。

## 1. これは何か

pysad は加速器設計コード SAD (K64 版、https://github.com/KatsOide/SAD) の Python への**直訳**。目的は、SAD のスクリプト (MAIN のラティス + FFS + SADScript) を
そのまま実行して、SAD と同じ数値を出すこと。一致の基準は SAD の印字精度 (1e-10〜1e-11)。物理を再導出せず、Fortran / C の式・分割数・反復順序・乱数の引き方を写す。
SAD の癖 (バグに見えるものも) は再現する (`docs/11_sad_quirks.md`)。SAD 自身のバグは再現せず `docs/19_sad_bug_report.md` に集める。

## 2. 方針 (作業者の指示。守ること)

- **翻訳の対象**: Fortran / C で書かれた物理計算は直接翻訳する。SADScript (`Packages/*.n`) は Python で書き直す (直訳ではない)。
  ファイル操作・リスト操作・GUI・描画は Python のライブラリ (open、list/numpy、matplotlib) に任せ、翻訳しない (`docs/15_sadhelp_purpose_policy.md`)。
- **照合の基準**: 「ユーザーのスクリプトが完走し、物理量が SAD と有意に一致する」。1e-10 以下の雑音、表示の桁・空白、丸めカオス後の乱数 seed は追わない。
- **SAD の参照先**: コマンドや関数の意味は CERN の SADHelp (https://hep-project-sad.web.cern.ch/SADHelp/SADHelp.html) か K64 版ソース同梱の `Documents/SADHelp.html`。
  KEK 公式サイトの SAD は別系統なので、ヘルプとしては参照しない (ライセンスの出典だけは KEK のページ)。数値の最終根拠は K64 版のソース (`src/*.f`、`src/*.c`、`Packages/*.n`)。
- **文書の書き方**: 解説文書 (`docs/*.md`) は履歴なしの現状のスナップショット。「最初 X だったが直した」のような経緯・日付付きの追記は `docs/history/` にだけ書く。
- **計算資源**: 他の計算と同居することが多いので、長い実行は 1 本ずつ・最低優先度で。動作テストでは誤差の実現数・チューンサーベイの格子・DA の周回数を間引く。

## 3. 構成

| 場所 | 内容 |
|---|---|
| `pysad/` | 本体。モジュールの地図は `docs/01_overview.md` §3、行数と準拠状況は `docs/15_implementation_summary.md` |
| `pysad/tests/` | pytest 277 件。SAD から取った参照値を埋め込み、または `data_*.npz` に持つ。ラティス `kek-ls_v03_12.lat` を読むテストあり (誤差入りラティスの 2 件は環境変数 `PYSAD_ERR_LATTICE` を指定したときだけ走る) |
| `run_sad.py` | `py -3.12 run_sad.py lattice.sad script.sad` (SAD の `gs` と同じ引数、複数ファイルは連結) |
| `docs/` | 解説 (現状)。`01_overview.md` から。API は `13_pysad_api.md` |
| `docs/history/` | この文書と `work_*.md` (項目ごとの経緯・照合の記録)、`review_notes.md` |

使わない旧モジュール: `pysad/ffs.py`, `matching.py`, `tracking.py`, `optics.py` (SAD 未照合の独自実装、非推奨)。新しいコードでは使わない。

## 4. 環境と実行

- Python 3.12 に numpy、scipy、numba、pytest。Windows では `py -3.12`。出力に日本語や記号が入るので `PYTHONIOENCODING=utf-8`。
- テスト: `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py --ignore=pysad/tests/test_track6n_ext.py`
  (無視する 4 本は sys.exit で終わる旧式のスクリプト)。約 45 秒。
- SAD 本体 (参照値を新しく取るとき): K64 版をビルドした `bin/gs` が必要 (このリポジトリには含まない)。WSL で `gs file.sad > file.out`。ラティスとスクリプトは**1 ファイルに連結**して渡す
  (`gs lat script` の 2 引数は標準入力待ちで止まる)。`NPARA = 1;` を入れる。WSL の SAD は `draw` / `TkWait` で落ちる・止まるので、描画の行は外す。
- 参照値は多くがテストに埋め込み済みなので、SAD 本体が無くても `pytest` は通る。SAD 本体が必要なのは、新しい機能を SAD と照合するときだけ。

## 5. 現状 (SAD と一致しているもの、していないもの)

- 一致: CALC (4D/6D、CELL/INS、結合、オフモメンタム、不安定ラインの続行)、EMIT / Emittance[] (RFSW × RADCOD、INTRA、報告の全行)、FIT/GO (解析微分、DESCEND、fuzz、2 段目、INS、フック)、
  追跡 (NORAD/RAD/FLUC、粒子損失)、DynamicApertureSurvey (NORAD)、RADINT、幾何 (GX…、GEO)、TRPT の加速、Touschek、誤差の実現、光学関数 (TransferMatrix 等)。詳細は `docs/15_implementation_summary.md` §3。
- ユーザーの実スクリプト 86 本 (esrf 52、pf2 34) を SAD と同条件で走らせた結果: 数値一致 62、一部一致 18、例外 1、SAD 側無効 4、対象外 1 (`docs/12_user_scripts_status.md`)。
  一部一致の原因は SAD 側の仕様 (FITP>1 の GO で SAD が再起動、不安定域の壊れた値、RAD+FLUC の乱数のカオス) が主で、pysad の不具合ではないと確認済み。
- GUI: SAD 自身の Tkinter.n / KBFrame.n / KBMainFrame.n を SADScript として解釈し、tkinter で表示する (`pysad/tclbridge.py`、`pysad/sadclass.py`、`docs/20_gui_tclbridge.md`)。
  KBMainFrame と demo-LabelFrame のウィジェット木は SAD (WSL の gs) と一致。入れ子の部品 (KBFFrameToggled / KBFFrameTabed)、Canvas の描画、見た目と操作の一致は未確認。環境変数 `PYSAD_SAD_PACKAGES` が無ければ無効。
- 未対応・制限は `docs/15_implementation_summary.md` §6 と `docs/01_overview.md` §6: SOL の emit・幾何、MULT の ANGLE≠0、TCAVI、COORD、MAP、UND、BEAMBEAM、POL (偏極)、
  SADHelp の関数のうち GaussianCoulomb・WakeFunction・RadiationField・ExternalMap など (優先度 低)、SynchroBetaEmittance (K64 に実体なし)。

## 6. 残課題 (優先順)

1. C6: CodCorr 系の補正量 (SAD の teigen が縮退固有値で丸めに依存する)。
2. C8: 小さな数値差 (esrf の amp_dep_ts_dz で 1.3e-4、原因未調査)。
3. fit の多解・FITP>1 (SAD が再起動するので参照値が取れない)、不安定域の閉軌道探索の初期値の差は、SAD 側の仕様なので記録のみ。
4. GUI の残り: 入れ子の部品 (demo-KBFrame1)、Canvas / TabFrame の描画、見た目と操作の比較 (`docs/20_gui_tclbridge.md` §6、経緯は `work_gui_tclbridge.md`)。
5. 次の実装 (利用者の指示、実例の台本なし): wake (WakeFunction ほか) と beam-beam (BEAMBEAM)。SADHelp と K64 のソースを根拠にする。
6. 仕上げ (文書の整理、テストの整備)。必要になったら SADHelp の優先度 低の関数。

## 7. 新しい機能を足す・直すときの手順

1. SAD のソース (`src/*.f`、`src/*.c`、`Packages/*.n`) で、実際に通るルーチンを特定する (itfunaloc、`word=='…'`、`REGn("Name")` で登録を辿る)。
2. 式・分割数・反復順序・乱数の引き方を Python に写す。数値微分や厳密な式への置き換えをしない (SAD の解析行列は厳密微分と 2e-8 ずれる)。
3. 小さなラティス (7 セル + RFC 程度) と短い SAD スクリプトを書き、SAD の `gs` で実行して pysad の出力と照合する。差が出たら SAD のソースに戻る。
4. 照合したケースを `pysad/tests/` に参照値付きで足す。
5. 文書: 解説文書の主題の節に現状として書く。SAD の癖は `docs/11_sad_quirks.md` に番号付きで足す。経緯と照合の記録は `docs/history/work_*.md` に書く。

## 8. 落とし穴

- SAD の FFS の `ON` / `OFF` は FFS 内では無効 (`INTRA;` / `NOINTRA;` と書く)。EMIT コマンドは GAUSS を立てる。SADScript の `1GEV` は SAD では `1*GEV` (=1)、pysad は単位として読む (意図的な拡張)。
- FFS[] の中のコマンドがエラーを出すと最上位の文 (Do など) が打ち切られる。LINE[…]= の未定義要素も同じ (`docs/11_sad_quirks.md` C25)。
- qcod が溢れても SAD は止まらず壊れた値で続行する。pysad は直前の値か NaN (C26)。
- 乱数: 既定の "SAD" 乱数 (合同法 + Box–Muller) を int64 で写してある。FLUC は要素外側・粒子内側の順に引く。損失周の最終 seed はカオスで一致しない。
- INTRA ON の EMIT 報告は 2 段 (零電流 → Touschek の表 → IBS 込み)、PBUNCH の既定は 1e10 (A40, A41)。
- 4D の CALC では空洞はドリフト、TRPT かつ RFSW のときだけ加速する。
- 数値の丸めに強く依存する比較 (ゼロ同然の値、1e-17 など) は一致を求めない。

## 9. リポジトリの扱い

- 追跡するのは `pysad/`、`docs/`、`run_sad.py`、`README.md`、`LICENSE`、`CLAUDE.md`、`.gitignore`。ベンチ、回帰スクリプト、ユーザーのラティスとスクリプトの写し、出力は含めない。
- ライセンスは SAD 同梱の LICENSE と同じ文面 (Copyright (c) 1986-2016, KEK Accelerator Laboratory)。pysad は SAD の派生物としてこの条件に従う。
