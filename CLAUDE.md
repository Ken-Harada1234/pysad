# CLAUDE.md

pysad は加速器設計コード SAD (K64 版) の Python への直訳 (v1 系)。v2 系 (v2.0.0 から) は K64 のバグの修正と機能の追加を含み、互換モードで切り替える。作業の前に次を読むこと。

1. `docs/history/00_RESUME.md`: 引き継ぎメモ (方針、構成、環境、現状、残課題、手順、落とし穴)。未解決の差の最新は `docs/history/work_unimpl.md`。
2. `docs/01_overview.md`: 全体像とモジュールの地図。
3. 必要に応じて `docs/11_sad_quirks.md` (SAD の癖)、`docs/12_pysad_api.md` (API)、`docs/46_unimplemented_and_unverified.md` (未対応・未確認の一覧)、`docs/47_v2_changes.md` (v2 の変更の一覧: v1 との違い、互換モード、要判断)。K64 側のバグは `docs/21_sad_bug_report.md` (B01〜B110)。
4. 要素・機能別: `docs/22_*.md`〜`docs/37_*.md` (GUI/Tcl、WAKE、BEAMBEAM、COORD、TCAVI、UND、MULT ANGLE、SOL、MAP、DELCOR、GaussianCoulomb・BBBrem、特殊関数・Fit、SAD パッケージ、FFS の小コマンド、POL・CSR、描画)。
   大項目の解説文書: `docs/38_wake_guide.md`、`39_beambeam_guide.md`、`40_csr_guide.md`、`41_tcltk_guide.md`、`42_graphics_guide.md`、`43_pol_guide.md`、`44_sol_guide.md`、`45_radlight_guide.md`。

要点:
- **互換モード**: 既定は `v2` (K64 のバグを直した動作)。gs (K64) と同じ結果が要るときは `PYSAD_COMPAT=k64` (または `pysad.compat.set_mode('k64')`)。既存の試験と回帰 (`run_regression.py`) は gs との比較なので k64 で走る (`pysad/tests/conftest.py`)。v2 の動作の試験は `test_v2_*.py` で `v2` を明示する。結果が変わる箇所を足すときは `compat.v2()` で分け、`docs/47_v2_changes.md` と `docs/21_sad_bug_report.md` に書く。
- 再導出せず直訳する (v1 の方針。v2 で変えたものは docs/47 の一覧)。数値の根拠は K64 版のソース (v2 の修正の根拠は解析式・独立な数値計算)。ヘルプは CERN の SADHelp か K64 版同梱の SADHelp.html (KEK 公式サイトは別系統なので参照しない)。
- 翻訳の対象は Fortran / C の物理計算。`Packages/*.n` は Python で書き直し、ファイル・リスト・GUI・描画は Python のライブラリに任せる。
  例外: GUI (Tkinter.n / KBFrame.n など) と SAD のパッケージ (CorrectOrbit など) は書き直さず、SAD 自身の `.n` を pysad のインタプリタで解釈する (`docs/22_gui_tclbridge.md`、`docs/34_sad_packages.md`)。
- テストは `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本を `--ignore`、手順は 00_RESUME.md §4)。937 通過・4 skip (941 件収集)、回帰 (`QUICK=1`) は 50 ベンチ FAIL 0 (2026-11-17。`PYSAD_SAD_PACKAGES` つき、k64 が既定)。GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip。SAD 本体が無くても通る。試験が走らせる SAD スクリプトと gs・pysad の画面表示は `pysad/tests/scripts/`、対応表は `pysad/tests/OUTPUTS.md` (`make_outputs_md.py` で作る)。
- 解説文書 (`docs/*.md`) は履歴なしの現状。経緯は `docs/history/` に書く。未確認・限界・使った条件 (エネルギー、磁場、周長) を必ず書く。
- gs の 10-19 ビルドは、保存済みの参照 (10-14) と値が違うことがある (00_RESUME.md §8)。
- コミットしてよいのは追跡対象 (00_RESUME.md §9) だけ。出力や中間ファイルを入れない。
