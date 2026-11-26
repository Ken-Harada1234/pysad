# pysad/tests/ (試験)

pysad の試験は、目的の違う 2 つのサブフォルダに分けてある。

| フォルダ | 目的 | 中身 | 実行 |
|---|---|---|---|
| [`script_tests/`](script_tests/README.md) | 機能ごとの自動試験 (スクリプトテスト) | 198 本の試験ファイル、収集 2534 件。保存した SAD (K64) の出力、解析式、独立な数値計算と比べる | pytest。全部で 1 時間 45 分ほど |
| [`regression_tests/`](regression_tests/README.md) | SAD (K64) との照合ベンチ (回帰テスト) | 50 ベンチ (Level 1〜14・16・22・24・26〜31)。保存済みの基準結果 (`.result`) と比べる。SAD の実行は不要 | `run_regression.py`。全部で約 17 分、`QUICK=1` で約 7 分 |

`pysad/tests/` の直下には、この README と、両方が共有するラティス `kek-ls_v03_12.lat` (大きなリング A: 3 GeV、周長 603 m、3563 要素) がある。

## 目次

1. [2 つの違い](#2-つの違い)
2. [全体の実行方法](#全体の実行方法)
3. [所要時間の目安](#所要時間の目安)
4. [結果の読み方](#結果の読み方)

## 2 つの違い

- **script_tests**: 要素・コマンド・関数ごとに、小さなラティスや SAD スクリプトを pysad で走らせて判定する。基準は、(a) 保存した SAD (K64) の出力、(b) 解析式、(c) SAD・pysad のコードを使わない独立な数値計算。既定で互換モード k64 で走る。v2 の動作の試験は試験の中で `v2` を明示する。
- **regression_tests**: SAD (K64) で取って保存した出力 (`bench_levelN.result`) と、pysad の出力を、ベンチごとの許容 (追跡は座標 1e-9 など) で突き合わせる。ベンチの内容はラティス全体の光学・エミッタンス・追跡・マッチング・DA・SADScript の関数など。
- 一部の試験 (`test_emit_*`・`test_rotate`・`test_small_items`・`test_unstable`・`test_main_track`・`test_fourier_fit`) は、`regression_tests/` の基準結果と判定の部品 (`cmp_*.py`) を読む。

## 全体の実行方法

リポジトリのルートから (Python 3.12)。

```
(環境変数は要らない。別の SAD の Packages を試すときだけ: set PYSAD_SAD_PACKAGES=<SAD の Packages フォルダ>。PowerShell: $env:PYSAD_SAD_PACKAGES="...")
py -3.12 -m pysad.nbcache --force                                numba のキャッシュを作り直す (初回・カーネルを直したあと)

# script_tests (pytest)
py -3.12 -m pytest pysad/tests -q -p no:cacheprovider `
    --ignore=pysad/tests/script_tests/test_bend_feeddown.py `
    --ignore=pysad/tests/script_tests/test_coord.py `
    --ignore=pysad/tests/script_tests/test_sext_radiation.py

# regression_tests
set QUICK=1
py -3.12 -u pysad/tests/regression_tests/run_regression.py
```

- `script_tests/conftest.py` が `PYSAD_SAD_PACKAGES` の既定を、pysad に同梱した SAD のコピー (`pysad/sad_copy/Packages`) にするので、環境変数なしでも GUI 系 (`test_plot.py`・`test_tclbridge.py` ほか) が走る (画面が無く `tkinter.Tk()` が作れないときは skip)。別の SAD の `Packages` で試すときは `PYSAD_SAD_PACKAGES` で指定する (指定が優先)。同梱のコピーの試験は `script_tests/test_v250_bundled_sad_copy.py`。SAD 本体 (gs) が無くても通る。
- `--ignore` する旧式の 3 本 (`test_bend_feeddown.py`・`test_coord.py`・`test_sext_radiation.py`) は、スクリプト形式で import 時に `sys.exit` を呼ぶため、pytest の収集を止める。
- `PYSAD_SLOW=1` を設定すると、遅い試験 (CSR MakeZL など) も走る。ほかの環境変数 (`PYSAD_ERR_LATTICE`・`PYSAD_PY314`・`V241_POL_CACHE`) は [`script_tests/README.md`](script_tests/README.md) にある。
- 試験と回帰を同時に走らせると、`test_plot.py::test_fitgaussian_plot_matches_gs` が PermissionError で落ちることがある (単独では通る)。続けて走らせる。
- 回帰は `QUICK=1` なしで全部走らせることもできる。指定したベンチだけ走らせる (`run_regression.py 3b 12`) 方法は [`regression_tests/README.md`](regression_tests/README.md)。

## 所要時間の目安

- script_tests 全部: 1 時間 45 分ほど。FITP の five の試験が約 10 分以上かかる。
- regression_tests: 全部で約 17 分、`QUICK=1` で約 7 分 (`run_regression.py` の各ベンチの目安の和。実測ではない)。

## 結果の読み方

- pytest: 「通過・skip・失敗」の件数が出る。失敗は 0 が期待値。skip は、環境変数が無い試験 (上記) と、POSIX の fork が要る試験などで、件数は環境で変わる (`-rs` で理由が出る)。参照ファイルの読み先を直し忘れると、`skipif(not exists)` 付きの試験が skip になって気づきにくいので、skip の件数が増えていないか確認する。
- 回帰: ベンチごとに PASS・`PASS*` (既知の不一致のみ)・FAIL・SKIP (`QUICK=1`) が一覧に出る。FAIL があれば終了コードが 1。
- 件数・時間の現在の値は、リポジトリ直下の [`README.md`](../../README.md) の「オリジナルの SAD (gs) との一致度」の節にある。
