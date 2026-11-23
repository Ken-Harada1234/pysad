# regression_tests/ (回帰テスト)

SAD (K64) との照合ベンチ 50 本 (Level 1〜14・16・22・24・26〜31) を、子プロセスで順に走らせて合否を一覧にする。
**基準はあらかじめ SAD (K64) で取って保存した出力 (`bench/bench_levelN.result`)** で、実行時に SAD (WSL) は要らない。pysad は互換モード k64 (`PYSAD_COMPAT` が未設定なら自動で k64) で走り、SAD の出力と数値を突き合わせる。
機能ごとの自動試験は姉妹フォルダ [`../script_tests/`](../script_tests/README.md)。上位の説明は [`../README.md`](../README.md)。

## 目次

1. [必要な環境](#必要な環境)
2. [実行手順](#実行手順)
3. [結果の読み方](#結果の読み方)
4. [50 ベンチの一覧](#50-ベンチの一覧)
5. [ファイルの構成と名前の規則](#ファイルの構成と名前の規則)
6. [基準の作り方](#基準の作り方)
7. [bench_extra/ (自動判定の対象外)](#bench_extra-自動判定の対象外)

## 必要な環境

- Python 3.12、pysad の依存 (numpy・numba など)。リポジトリのルートの pysad を使う (ベンチの pysad スクリプトは `sys.path` にルートを足す)。
- SAD (K64) は不要。ベンチの台本 `.sad` は、pysad の `run_sad.py` が実行する。
- 目安時間は `run_regression.py` の値 (実測の記録ではなく目安)。合計は全部で約 17 分、`QUICK=1` で約 7 分。

## 実行手順

リポジトリのルートから (どのフォルダからでも動く)。

```
PYTHONIOENCODING=utf-8 py -3.12 -u pysad/tests/regression_tests/run_regression.py            # 全部 (約 17 分)
QUICK=1 PYTHONIOENCODING=utf-8 py -3.12 -u pysad/tests/regression_tests/run_regression.py    # 遅い 4 本 (7・12t・13a・16a) を飛ばす (約 7 分)
py -3.12 -u pysad/tests/regression_tests/run_regression.py 3b 12 13a                          # 指定したものだけ (キーは下の表)
```

PowerShell では `$env:QUICK="1"; $env:PYTHONIOENCODING="utf-8"` を設定してから実行する。各ベンチの出力は `regression_tests/regression_out/<キー>.txt` に保存される (実行ごとの出力なので git には入れない)。終了コードは FAIL があれば 1。

## 結果の読み方

- **PASS**: 全項目が許容内。`PASS*`: 既知の不一致 (KNOWN) だけで、ほかは許容内。**FAIL**: 許容を超える項目がある (一覧の最後に項目と差・許容が出る)。**SKIP**: `QUICK=1` で省略。
- 許容 (`run_regression.py` の定数、`docs/27_benchmarks.md`): 追跡は座標の絶対差 1e-9。CALC (オンモメンタム、閉軌道 0) は 1e-9。CALC (閉軌道・結合・オフモメンタム経由) は相対 1e-7 (SAD の qcod の収束と pysad の数値微分の丸めのため)。Emittance[] (6D) は相対 1e-6、SAD の表示が 5 桁のときは相対 5e-5。マッチングは 1e-9。DA は地図の文字列と Score が一致 (損失周数は境界粒子 1 個の ±1 周まで既知)。Level 1・2 の Twiss は SAD の表示が 6 桁なので絶対 1e-6、RADINT は 7 桁表示なので相対 1e-6、色収差の 1 次係数は 1e-5。

## 50 ベンチの一覧

キーは実行時の指定名。「実行」の (A) は pysad スクリプト `bench/bench_levelN_pysad.py` を直接、(B) は `run_sad.py` に `bench/bench_levelN.sad` を渡す。(A) は 1〜12・12t・12f・27a・27b の 17 本、(B) は残り 33 本。

| キー | Level | 何を確かめるか | 実行 | 許容 | 目安時間 | QUICK |
|---|---|---|---|---|---|---|
| 1 | L1 | 線形光学・COD (QD に DX)・RADINT・Emittance[] (小リング) | (A) | Twiss 絶対 1e-6 (SAD 表示 6 桁)、COD 1e-10+相対 1e-7、RADINT 相対 1e-6、Emittance 相対 1e-6 | 3 s |  |
| 2 | L2 | オフモメンタム Twiss (dp=0, ±0.01, ±0.02) と K1 付き BEND (小リング) | (A) | Twiss 絶対 1e-6、RADINT 1e-6、Emittance 1e-6 | 6 s |  |
| 3a | L3A | 6 極と色収差 (DP0 スキャン 21 点、3 次当てはめ) | (A) | チューン 1e-7 (相対)、色収差 c1 は 1e-5 | 10 s |  |
| 3b | L3B | 6 極を通る COD (DELX/DELY 1 台) | (A) | 閉軌道 1e-12+相対 1e-7、光学 1e-8+相対 1e-7、Emittance 1e-6 | 6 s |  |
| 3c | L3C | x-y 結合 (DTHETA)・Emittance[] | (A) | 同上 (結合 R は 1e-7) | 4 s |  |
| 3d | L3D | 強さの誤差 (DK/DDK) | (A) | 同上 | 3 s |  |
| 4 | L4 | RF 空洞 (RFSW・RADCOD) の Emittance[] | (A) | Emittance 相対 1e-6、閉軌道 1e-9 | 4 s |  |
| 5 | L5 | 輸送ライン (INS) とフリンジ | (A) | Twiss 絶対 1e-6、軌道 1e-9、位相 1e-7 | 1 s |  |
| 6 | L6 | 追跡: RAD OFF・RAD ON NOFLUC (5 粒子、1/10/100 周) | (A) | 座標の絶対 1e-9 | 9 s |  |
| 7 | L7 | 追跡: 誤差入り要素・K1≠0 の BEND | (A) | 1e-9 | 2 min | 飛ばす |
| 8 | L8 | 追跡: 四極フリンジ・長い空洞・MULT | (A) | 1e-9 | 35 s |  |
| 9 | L9 | 追跡: ソレノイド (SOL BZ=2 T) | (A) | 1e-9 | 9 s |  |
| 10 | L10 | 追跡: RAD+FLUC (SAD の乱数列) | (A) | 1e-9 | 9 s |  |
| 11 | L11 | 追跡: MULT の K0/K1/K2 とずれ・回転 | (A) | 1e-9 | 17 s |  |
| 12 | L12 | ESRF 型の実ラティス (3 GeV、20 セル、3564 要素): 読み込み・CELL CALC・emit | (A) | 読み込み・CALC 1e-9、emit 相対 5e-5 (SAD 表示 5〜8 桁) | 10 s |  |
| 12t | L12 | 同: 追跡 (純 Python の track6) | (A) | 1e-9 | 1 min | 飛ばす |
| 12f | L12 | 同: 追跡 (numba の track6n) | (A) | 純 Python との一致を判定 | 20 s |  |
| 13a | L13A | マッチング FIT/FREE/GO (小リング)。13A-2 は収束しない場合の反復も比較 | (B) | 1e-9 (タグ: K1、残差、Twiss) | 3〜4 min | 飛ばす |
| 13b | L13B | マッチング (ESRF 型、7 ファミリー) | (B) | 1e-9 | 10 s〜 |  |
| 14a | L14A | SADScript: esrf/calc_chromaticity | (B) | (check_level14a) | 30 s |  |
| 16a | L16A | DynamicApertureSurvey (NORAD) | (B) | 地図の文字列と Score が一致。損失周は ±1 周まで既知 | 30 s〜 | 飛ばす |
| 27a | L27A | MAIN の TRACK: 規格化振幅→初期分布、smear、チューン (wtunem)、RAD/FLUC (SAD 側は SADScript 版) | (A) | (check_level27: 保存した SAD 版の出力と比較) | 10 s |  |
| 27b | L27B | MAIN の TRACK: 粒子が落ちる場合 (損失周、落ちた粒子の x,y,z) | (A) | 同上 | 10 s |  |
| 22e | L22E | INS/CELL: 四極ずれ・DP0≠0 (閉軌道≠0、ESRF 型) | (B) | 行比較 1e-10 | 10 s |  |
| 26a | L26A | EMIT 報告の全行 (ESRF 型。RFSW × RADCOD の 4 通り) | (B) | 報告行の数値 (雑音 1e-11、相対 1e-9) | 12 s |  |
| 26b | L26B | 四極 FRINGE=1〜3・ずれ・ロール付き四極・6 極 (CALC と Emittance[]) | (B) | 行比較 3e-9 | 3 s |  |
| 26c | L26C | Orbit dilation (FSHIFT): EMIT 報告・Emittance[]・CALC6D | (B) | 報告行+行比較 1e-9 | 8 s |  |
| 26d | L26D | RFSW OFF+RADCOD ON の Emittance[]・CALC6D (ESRF 型) | (B) | 行比較 1e-9 | 10 s |  |
| 28a | L28A | 不安定ライン (INTRES/HALFRES ON) の CELL CALC の続行・Unstable 表示・GO (nstab) (V・H・HV の 3 通り) | (B) | 1e-8 (一部のタグは 1e-3〜1e-6 に緩める) | 12 s |  |
| 28b | L28B | 不安定ライン: INTRES/HALFRES の 4 通りと GO | (B) | 1e-8 (一部緩める) | 12 s |  |
| 28c | L28C | 不安定な光学 (V: QD ×1.19) の EMIT 報告の全行 | (B) | 雑音 1e-10 | 10 s |  |
| 28d | L28D | 不安定な光学 (H: QF ×0.9) の EMIT 報告の全行 | (B) | 雑音 1e-4・相対 3e-5 (T の悪条件) | 10 s |  |
| 28f | L28F | 不安定な光学の Emittance[] (V・H): Stable=0、戻り値の全量 | (B) | 行比較 3e-5 | 15 s |  |
| 30a | L30A | INS+CALC6D (qcell6d の cell=.false. の枝): 先頭 MARK の 28 成分から 6D 伝搬 (5 通り) | (B) | 行比較 1e-9 | 20 s |  |
| 30b | L30B | RFSW ON+RADCOD OFF+FSHIFT≠0 の Emittance[] (tcod の 6D Newton) | (B) | 行比較 1e-10 | 60 s |  |
| 29a | L29A | RADINT コマンド (ESRF 型 6 通り: RFSW/EMIT 後/ずれ・ロール/DP0/COD) | (B) | 相対 1e-6 | 10 s |  |
| 29a2 | L29A | RADINT: 小リング (E1/E2 付き BEND、RF なし/あり、結合) | (B) | 相対 1e-4 | 3 s |  |
| 29b | L29B | DISP の位置範囲 (数値・名前・NAME.n±k・パターン・周回・ALL・R/O/E モード) | (B) | ブロック比較 1e-9 (D37 を除く) | 10 s |  |
| 29c | L29C | CALC6D;/CALC4D; 単独 (CALC なし) の挙動 | (B) | 行比較 1e-9 | 15 s |  |
| 29d | L29D | FFS 内の ON/OFF flag | (B) | ブロック比較 1e-10 (E17 を除く) | 5 s |  |
| 29e | L29E | FSHIFT の追跡: NORAD / RAD / RAD+FLUC (乱数も) | (B) | 行比較 5e-11 | 10 s |  |
| 31a | L31A | 誤差コマンド (DELX/DELY/DK/DTHETA/DDK、UNIFORM/COHERENT/ADD) と CALC (ESRF 型、128 タグ) | (B) | タグ比較 (1e-10, 1e-7) | 10 s |  |
| 31b | L31B | Twiss[]/LINE[] の位置パターンの返り値 (一致 1 個はスカラー、0 個は {}) (ESRF 型) | (B) | タグ比較 (1e-10, 1e-8) | 5 s |  |
| 24a | L24A | SADScript のリスト・文字列・行列関数 (559 タグ。ラティス不要) | (B) | タグ比較。既知の不一致 10 タグ (T118 T191 T195 T220 T221 T222 T229 T231 U143 W075) | 3 s |  |
| 24b | L24B | CHROMATICITY コマンド (ESRF 型、全 960 要素) | (B) | タグ比較 | 15 s |  |
| 24b2 | L24B | CHROMATICITY (MULT の K1/K2、DP0≠0) | (B) | タグ比較 | 5 s |  |
| 24c | L24C | FFS の REPEAT/UNTIL、IF/ELSEIF/ELSE/ENDIF | (B) | タグ比較 | 5 s |  |
| 24d | L24D | CalculateOptics (4 次元、cell/INS・区間・結合・分散) | (B) | タグ比較 (1e-8, 1e-7) | 5 s |  |
| 24d2 | L24D | CalculateOptics (ESRF 型、6 次元) | (B) | タグ比較 (1e-8, 1e-7) | 60 s |  |
| 24f | L24F | Element/SetElement/FLAG/BeamLine/WriteBeamLine/USE BeamLine | (B) | タグ比較。既知 1 タグ (G004) | 5 s |  |

QUICK=1 で飛ばすのは 7・12t・13a・16a の 4 本。ラティスの種別は、`run_regression.py` の説明とベンチ冒頭の記述による (26b・24 系の一部は未確認)。Level 1〜12 の小リングは Level 4・6 で共有する 2.5 GeV 前後のリング、ESRF 型の大きなリング (3 GeV、20 セル、3564 要素) は台本の冒頭に埋め込む。pysad スクリプトは 12・12t・12f・13b で `pysad/tests/kek-ls_v03_12.lat` を直接読む。

## ファイルの構成と名前の規則

```
regression_tests/
  run_regression.py         実行と判定 (BENCHES に 50 本の定義)
  cmp_tags.py cmp_lines.py cmp_blocks.py cmp_radint.py cmp_emit_report.py   判定の部品 (script_tests の一部の試験も import する)
  bench/                    50 本が使う台本・pysad スクリプト・基準結果 (+ lat_extra.sad・modules_general.n: 台本が相対パスで読む補助)
  bench_extra/              自動判定の対象外のもの (手動確認用)
  regression_out/           実行ごとの出力 (追跡しない)
```

N = Level 番号 + 枝 (3a、13b、26c など)。

| ファイル | 内容 |
|---|---|
| `bench_levelN.sad` | SAD (K64) 用の台本。`run_sad.py` で pysad にも走らせる |
| `bench_levelN.result` | SAD (K64) で取った基準の出力 (保存したもの。判定の基準) |
| `bench_levelN_pysad.py` | pysad 側のスクリプト。差を印字する。`.result` を読んで比べるものもある (1・2・7・12 ほか) |
| `bench_levelN_pysad.result` | pysad の出力の保存 (参考。判定には使わない。`bench_extra/`) |
| `bench_level12_track_pysad.py`・`bench_level12_track_fast.py` | 12t・12f 用 (純 Python と numba の追跡) |
| `bench_level27a.sad`・`bench_level27a_main.sad`・`bench_level27b_main.sad`・`bench_level27a_pysad.py` | 27a・27b 用。MAIN の TRACK は SAD (K64) では動かないため、SAD 側は SADScript で書いた版 (`bench_level27a.sad`)、pysad 側は MAIN 形式 (`_main.sad`) |
| `bench_levelN.pysad` | (24a〜24f の一部) pysad の出力の保存 |
| `cmp_*.py` | 判定の部品 (タグ・行・ブロック・RADINT・EMIT 報告の比較) |

実行時に生成される `bench_level27a_pysad.out`・`bench_level27b_pysad.out`・`datas.txt` は追跡しない (リポジトリ直下の `.gitignore`)。

## 基準の作り方

SAD (K64) がある環境 (WSL の `gs`) で行う。この手順は、このリポジトリでは実行していない (保存済みの基準を使うため)。

1. `bench/` へ移り、台本を SAD で走らせて画面出力を保存する (例: `gs bench_level26a.sad < /dev/null > bench_level26a.result`)。台本の相対パスの読み込み (`lat_extra.sad`・`modules_general.n`) は `bench/` が基準。
2. `run_regression.py` の `BENCHES` に、キー・実行するファイル・判定関数・許容を足す。判定関数は `check_tags` (タグ付き `Print`)、`check_lines`、`check_emit_report`、`check_radint`、`check_blocks` などから選ぶ。
3. `py -3.12 -u pysad/tests/regression_tests/run_regression.py <キー>` で PASS を確認する。SAD の値が小数 6〜7 桁表示のときは許容を表示桁に合わせる (上の許容)。

## bench_extra/ (自動判定の対象外)

`run_regression.py` の 50 本に入らないベンチの台本と基準 (245 ファイル)。自動判定は無く、必要なときに手で走らせる (`py -3.12 -u run_sad.py pysad/tests/regression_tests/bench_extra/bench_levelNN.sad`)。内容: Level 15・16b〜l・17〜21・23・25 の台本と基準 (`bench_level16*.sad` は DA の大規模な走査、`bench_level20_dasv_*` は DAsv[] の長い計算)、`bench_item6_*` (Emittance 表示・DISP の名前の台本。`script_tests/test_emit_jac.py` が読む)、`bench_level25a〜f` (`test_rotate.py` が読む)、調査用 (`bench_level14a_dbg*.sad`・`bench_match6d_check.py`)、実行の記録 (`.npy`・`.start`・`.end`・`.time`)。説明は `docs/27_benchmarks.md`。
