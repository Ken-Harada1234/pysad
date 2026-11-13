# CLAUDE.md

pysad は加速器設計コード SAD (K64 版) の Python への直訳。作業の前に次を読むこと。

1. `docs/history/00_RESUME.md`: 引き継ぎメモ (方針、構成、環境、現状、残課題、手順、落とし穴)。
2. `docs/01_overview.md`: 全体像とモジュールの地図。
3. 必要に応じて `docs/11_sad_quirks.md` (SAD の癖)、`docs/13_pysad_api.md` (API)、`docs/12_user_scripts_status.md` (実スクリプトの状況)。

要点:
- 再導出せず直訳する。数値の根拠は K64 版のソース。ヘルプは CERN の SADHelp か K64 版同梱の SADHelp.html (KEK 公式サイトは別系統なので参照しない)。
- 翻訳の対象は Fortran / C の物理計算。`Packages/*.n` は Python で書き直し、ファイル・リスト・GUI・描画は Python のライブラリに任せる。
  例外: GUI (Tkinter.n / KBFrame.n など) は書き直さず、SAD 自身の `.n` を pysad のインタプリタで解釈して tkinter に表示する (`docs/20_gui_tclbridge.md`)。
- テストは `py -3.12 -m pytest pysad/tests -q` (旧式の 4 本を `--ignore`、手順は 00_RESUME.md §4)。SAD 本体が無くても通る。
- 解説文書 (`docs/*.md`) は履歴なしの現状。経緯は `docs/history/` に書く。
- コミットしてよいのは追跡対象 (00_RESUME.md §9) だけ。出力や中間ファイルを入れない。
