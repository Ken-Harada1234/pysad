# 作業ログ: pysad の試験が走らせる SAD スクリプトと gs・pysad の画面表示 (2026-11-17)

## 依頼と方針

pysad の試験にも、gs の結果のリダイレクトと pysad の結果のリダイレクトを (画像出力があるものは画像も込みで) ファイル名で分かるように置く。
`example/` と pysad の試験は別のスクリプトなので、それぞれのフォルダにそのフォルダのスクリプトの結果を置く。`example/` の出力を試験で読むことはしない (最初に作りかけた `test_examples.py`・`examples_util.py` は削除)。

## したこと

- 試験の実行中に `pysad.ffscmd.run_sad_string` を差し替えるプラグインで、台本の文字列と画面出力 (標準出力 + 標準エラー) を全試験から集めた (試験のコードは変えていない。判定も同じ)。
  重複 (同じ試験の同じ台本) は除いた。328 本の台本、46 本の試験ファイル。
- `pysad/tests/scripts/NAME.sad` (台本)、`NAME.pysad.txt` (pysad の画面出力)、`NAME.gs.txt` (同じ台本を WSL の gs に `gs NAME.sad < /dev/null` で 40 秒打ち切り、3 並列、低優先度で走らせた画面表示) を置いた。
  NAME は `<試験ファイル名 (test_ を除く)>__<試験関数 (クラス名を含む)>[_2, _3, …]`。
- `pysad/tests/make_outputs_md.py` で `pysad/tests/OUTPUTS.md` を作る (試験関数とファイルの相対パスの表、対象外の試験、`pysad/tests/` の入力ファイル、回帰ベンチの場所)。

## 限界・未確認

- 試験の台本は pysad 用の整形で、gs でそのまま走るとは限らない。gs の画面表示と pysad の画面表示の一致は判定していない (試験自身の判定は埋め込みの SAD 参照値)。gs が SIGSEGV で落ちたものは OUTPUTS.md の列に書いた (6 本)。
- `scripts/` は取り出した写しで、試験は読まない。台本を直したら取り出し直す (プラグインは作業用で、リポジトリには含めていない)。
- `run_sad_string` を通らない試験 (`run_sad`・`load_sad_string`・Python 関数の直接呼び出し) は取り出していない。OUTPUTS.md の「対象外」に挙げた。
- 画像を書き出す試験は無い (描画・GUI の図は `example/` の PNG)。
- 回帰ベンチ (`bench_level*.sad`・`.result`、`run_regression.py`) はコピーを作らず、場所と対応を OUTPUTS.md に書いた。
- `.gs.txt` の `Today:` の時刻、`.pysad.txt` の一時名は実行ごとに変わる。
- gs のビルドは 2026-10-19 (保存済みの参照 10-14 と値が違うことがある。00_RESUME.md §8)。

## 反映した文書

`CLAUDE.md`、`README.md`、`docs/01_overview.md`、`docs/history/00_RESUME.md` §4。試験の数 (575 通過・3 skip) は変わらない。
