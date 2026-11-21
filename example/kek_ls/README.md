# KEK-LS の例集 (3 GeV、20 セル、ESRF 型 HMBA)

KEK-LS の設計・修正・パラメータ表示のための SAD スクリプト一式です。KEK-LS は ESRF 型の HMBA (Hybrid Multi-Bend Achromat) で、ESRF の Pantaleo Raimondi 氏と Simone Liuzzo 氏の協力で設計しました。ベースのラティスは設計版 `v03_68` (3 GeV、周長 570.72 m、標準セル 20 個、RF 2 MV・ハーモニック数 952、チューン (48.58, 17.62)、εx 約 134 pm rad)。元の設計用スクリプトの集まりから、ラティス・マッチング・パラメータ表示・誤差を入れた補正を取り出し、`v03_68` で動くように整え、説明を付けました。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点で、gs (SAD K64) と pysad の出力を並べて確認しています ([`example/README.md`](../README.md) と同じ作法)。

| フォルダ | 内容 |
|---|---|
| [`lattice/`](lattice/README.md) | ラティス `kek_ls_v03_68.sad`、構成図、パラメータ表、旧版 (v03_12) との違い、入射セル |
| [`design/`](design/README.md) | 設計・マッチング (セルのチューン、六極、入射セル、リングを閉じる、縮小、セル数、チューンサーベイ) |
| [`display/`](display/README.md) | パラメータ表示 (εx・σδ・ξ・αc・U0、オプティクス、チューンシフト) |
| [`errors/`](errors/README.md) | 誤差入りの計算と COD の補正 |
| `modules/` | 共通の関数 (`modules_general.n`・`modules_amp_dep_ts.n`・`modules_cell_match.n`・`modules_errors.n`) |
| [`ORIGINAL_MAP.md`](ORIGINAL_MAP.md) | 新しい名前と元のスクリプトの対応表 |

## 全体の流れ

```
ラティス (kek_ls_v03_68.sad)
   │  cell calc emit: 周長・チューン・εx          lattice/KEKLS_cell_calc_emit
   ▼
チューンのマッチング (標準セル 1 個の位相進み)       design/KEKLS_match_cell_tune
   ▼
六極のクロマティシティ調整                          design/KEKLS_match_chrom
   ▼
入射セルのマッチング                               design/KEKLS_match_injcell
   ▼
リングを閉じる (入射セルの長直線部)                  design/KEKLS_match_ring_geo
   ▼
パラメータ表示 (εx・ξ・オプティクス・チューンシフト)   display/*
   ▼
誤差・補正 (COD・分散)                             errors/*
```

設計の変形 (`design/KEKLS_scale_mini` 縮小、`KEKLS_scale_ncell` セル数、`KEKLS_tune_survey_cell` チューンサーベイ) は、この流れの途中のどこからでも、ラティスから独立に走ります。

## スクリプトの一覧

gs と pysad の一致 (k64 モード): ○ 一致 (数値の相対差 1e-5 以内)、△ 差あり、× 不一致。所要時間はこの PC での 1 回の計測 (gs は WSL、pysad は Python 3.12)。

| 例 | 元のスクリプト | 何をするか | 入力 | 出力 | 所要時間 (gs / pysad) | gs と pysad |
|---|---|---|---|---|---|---|
| [`lattice/KEKLS_cell_calc_emit`](lattice/KEKLS_cell_calc_emit.md) | - | RINGC・RING の周期解とエミッタンス | lattice/kek_ls_v03_68.sad | 周長・チューン・εx・σδ | 1 s / 3 s | ○ |
| [`lattice/KEKLS_cell_optics`](lattice/KEKLS_cell_optics.md) | disp_optics | 標準セルの配置・オプティクスの数表と図 | 同上 | 数表 160 行、PNG (`../tools/kek_ls_figs.py`) | 1 s / 2 s | ○ |
| [`design/KEKLS_match_cell_tune`](design/KEKLS_match_cell_tune.md) | matching2_01 前半 | セルのチューンのマッチング | 同上 | セルの位相進み・四極の K1 | 1 s / 2 s | △ |
| [`design/KEKLS_match_chrom`](design/KEKLS_match_chrom.md) | matching2_01 後半 | 六極のクロマティシティ調整 | 同上 | ξ・六極の K2 | 1 s / 9 s | △ |
| [`design/KEKLS_match_injcell`](design/KEKLS_match_injcell.md) | matching2_03_INJ | 入射セルのマッチング | 同上 | 入射点の Twiss・QS0〜QS5 の K1 | 1 s / 2 s | ○ |
| [`design/KEKLS_match_ring_geo`](design/KEKLS_match_ring_geo.md) | matching2_02_geo | 入射セルを含めてリングを閉じる (ジオメトリ) | 同上 | DRID1INJ・GX・GY | 1 s / 2 s | ○ |
| [`design/KEKLS_match_repeat`](design/KEKLS_match_repeat.md) | (新規) | マッチングの再現性の診断 | 同上 | K1・K2 (3 回の繰り返し) | 1 s / 20 s | △ |
| [`design/KEKLS_scale_mini`](design/KEKLS_scale_mini.md) | matching23_v03_35_mini | セルの一様な縮小 | 同上 | セル長・位相進み・K1 | 1 s / 4 s | ○ |
| [`design/KEKLS_scale_ncell`](design/KEKLS_scale_ncell.md) | matching21_v03_35 | セル数を変える | 同上 | j セルのリングの εx・チューン・周長 | 1 s / 3 s | ○ |
| [`design/KEKLS_tune_survey_cell`](design/KEKLS_tune_survey_cell.md) | matching05 + modules_tune_sv_v03_68.n | セルのチューンサーベイ | 同上 | ξ・εx・六極の強さの 25 点の地図 | 10 s / 160 s | △ |
| [`display/KEKLS_ring_params`](display/KEKLS_ring_params.md) | calc_emit.sad 前半、disp、chromaticity_plot | リングのパラメータ (εx・σδ・減衰・αc・U0・ξ) | 同上 | emit の表・ξ | 1 s / 4 s | ○ |
| [`display/KEKLS_optics_table`](display/KEKLS_optics_table.md) | disp_optics、disp、draw_cad の入力 | リングのオプティクス・配置の数表と図 | 同上 | 要素ごとの表 3226 行、PNG | 1 s / 3 s | ○ |
| [`display/KEKLS_calc_H`](display/KEKLS_calc_H.md) | calc_H.sad、calc_emit.sad 後半 | H 関数の積分によるエミッタンス | 同上 | εx (SAD と 2.3% 差) | 1 s / 3 s | ○ |
| [`display/KEKLS_chromaticity`](display/KEKLS_chromaticity.md) | chromaticity_plot.sad | チューンの運動量依存 | 同上 | ν(δ)・PNG | 1 s / 5 s | ○ |
| [`display/KEKLS_amp_dep_ts`](display/KEKLS_amp_dep_ts.md) | amp_dep_ts*.sad | 振幅依存チューンシフト | 同上 + modules/modules_amp_dep_ts.n | ν(振幅)・PNG | 1 s / 36 s | ○ |
| [`display/KEKLS_sextupole_tune_scan`](display/KEKLS_sextupole_tune_scan.md) | dp_sv_sx.sad | 六極を変えたチューンシフト | 同上 | Σ|Δν| の 5×5 の地図・PNG | 2 s / 18 s | △ |
| [`display/KEKLS_scan_target`](display/KEKLS_scan_target.md) | dp_sv.sad、dx_sv*.sad、dy_sv*.sad | マッチングの目標値の振り (δ 依存・ADTS) | 同上 + modules/modules_cell_match.n | ξ・ν(δ)・ADTS | 5 s / 119 s | △ |
| [`errors/KEKLS_errors_cod`](errors/KEKLS_errors_cod.md) | print_all_error.sad ほか print_* | 誤差の種類ごとの COD・分散・チューン | 同上 + modules/modules_errors.n | COD・分散・チューンのずれ | 1 s / 6 s | △ |
| [`errors/KEKLS_correct_cod`](errors/KEKLS_correct_cod.md) | correct_cod.sad | COD・分散の補正 | 同上 + modules/modules_general.n・modules_errors.n | 補正前後の COD・分散 | 1 s / 37 s | △ |
| [`errors/KEKLS_cod_response`](errors/KEKLS_cod_response.md) | write_cod_response.sad | COD 応答行列と解析式の照合 | 同上 | 応答行列の代表値・解析式との差 | 1 s / 9 s | ○ |

数値の判定は `tools/kek_ls_cmp.py` が、Print が出す行 (指定の接頭辞) の数値で自動で付けた (`kek_ls_cmp.py` の `SPEC`)。

## 使い方

```bash
# リポジトリ直下から。台本は自己完結 (先頭にラティスを貼ってある)。モジュールは台本のあるフォルダをカレントにして相対パスで読む
py -3.12 example/tools/run_kekls.py example/kek_ls/design/KEKLS_match_injcell.sad     # gs と pysad で走らせ、.gs.txt・.pysad.txt を作る
PYSAD_COMPAT=k64 PYTHONUTF8=1 py -3.12 -c "from pysad.ffscmd import run_sad_string; run_sad_string(open('KEKLS_match_injcell.sad', encoding='utf-8').read(), cwd='.')"   # pysad だけ (台本のフォルダで)
py -3.12 example/tools/mk_kekls.py                 # tools/kek_ls_src/*/*.body とラティスから .sad を作り直す
py -3.12 example/tools/kek_ls_figs.py              # 図 (PNG) を gs の出力から作る
py -3.12 example/tools/kek_ls_cmp.py               # gs と pysad の判定
```

- gs は WSL で `gs 台本.sad < /dev/null`。`gs` は `FFS;` の前にラインや要素の定義が必要で、`Get[]`・`READ` でラティスの台本は読めない (確認済み)。そのため各 `.sad` は先頭にラティスを貼っている (7 KB)。
- `Get["../modules/..."]` は台本のあるフォルダが基準。gs・pysad とも、そのフォルダをカレントにして走らせる。
- 並列数は `NPARA = 8` 以下。長い計算は例の設定 (格子・種・周回数を小さくした設定) で数分以内にしてあり、元の設定 (完全な設定) は各フォルダの README の「元の設定と例の設定」の表に、変数を変えるだけで走る形で書いた。
- pysad は `PYSAD_COMPAT=k64` (gs と同じ結果)。v2 (既定) モードは `Emittance[]` が数 % 違う (lattice/README.md の表)。
- 試験: `py -3.12 -m pytest pysad/tests/test_kek_ls_example.py -q` (ラティスの `cell calc emit` の値と、速い台本の数本)。

## 整理の方針

- **モジュール化**: 使う箇所が 1〜2 本のモジュールは台本に取り込んだ (`modules_tune_sv_v03_68.n` のうち `InitMag`・`MatchingARCA`・`ChromCor` は 2 本が使うので `modules/modules_cell_match.n`、`modules_general.n`・`modules_amp_dep_ts.n`・誤差の `modules_errors.n` も 2 本以上で共通なので `modules/`)。ほかの `modules_*.n` (ID 入りラティス・別のリング・PF・v03_35 用) は取り込まない。`Get[]` は相対パス。
- **コメント**: 有用な説明は日本語で台本の先頭に整理し、デバッグの残骸 (使っていない `!` の行、連続する `stop;`、コメントアウトされた旧い試行) は削除した。式・順序は元のまま。
- **元の挙動を変えない**: 整理の前後で、元の台本 (ラティスを貼って、元のモジュールと同じフォルダで、`NPARA` だけ 4 にして走らせたもの) の gs の出力と、整理後の gs の出力を比べた。比べられたもの: `KEKLS_match_cell_tune` (1 点目一致、残りは 5e-6 以内)、`KEKLS_match_injcell` (1e-5 以内)、`KEKLS_calc_H` (1e-10)、`KEKLS_amp_dep_ts` (5e-9)、`KEKLS_chromaticity` (δ ≤ 2.7% で 2e-6)。元が v03_35 用・動かない・別のラティス用のもの (`scale_*`・`tune_survey_cell`・`match_ring_geo`・`errors/*`) は直接は比べられていない (各 `.md` に理由)。
- **元の台本の問題**: 動かない・前提のラティスが無い・名前の誤りは直して、各フォルダの README と各 `.md` に書いた (`RINGI` → `RING`、`free BPI6 K1` が効かない、v03_35 の名前、`BPM.*` → `PM.*`、`DELX`・`DELY` が COD に効かない、補正磁石の K0 に `NOEXPAND` が要る、など)。

## 未確認・限界

- **pysad の FIT が gs と違う場面がある** (`KEKLS_match_repeat`・`tune_survey_cell`・`scan_target`。FitFunction を fit 条件と併用するマッチング。同じ入力の繰り返しで結果が変わる)。pysad のコードは変えていない。gs の結果が基準。
- 誤差入りの閉軌道は、位置ずれが大きい種 (種 17) で gs と pysad が別の解に落ちる。
- 解析式 (`KEKLS_calc_H`) は SAD の εx と 2.3% 違う (原因は未確認。解析式を SAD に合わせる補正はしていない)。
- ξ の値が作業メモ (RINGC の (+0.133, −0.658)) と違う (差分で +0.119, −0.684。メモの値の出し方は未確認)。
- 粒子の追跡 (DA を含む) は含めていない (チューンサーベイと誤差の元の台本は DA を計算するが、重いので外し、`DAsv[]` を `modules_general.n` に残した)。検証の基準は SAD の粒子追跡だが、この例集の計算は線形光学・周期解・マッチングで、追跡による検証はしていない。作動点が整数・半整数共鳴に近い点 (チューンサーベイの格子、`scale_ncell` の j = 14) は表に注記した。
- 磁場の上限・アパーチャは元のファイルに無く、確認していない。使った条件は 3 GeV、周長 570.72 m、RF 2 MV (誤差なしの設計ラティス)。
- 元のスクリプトの所要時間の記録は元のフォルダに無い (ファイルの更新日時は 2016 年〜2026 年)。完全な設定の所要時間は、例の設定からの比例の外挿 (DA は測っていない)。

## 参考文献

- KEK-LS Conceptual Design Report (CDR)、KEK、2017 年 5 月 22 日。<https://www2.kek.jp/imss/notice/assets/2017/05/22/KEKLS_CDR_170522.pdf>

CDR は設計の概要で、この例集の `v03_68` はその後に更新された版です。ここに書いた値は、すべて `v03_68` を gs・pysad で計算した実測値です。CDR の値との並べた表は、CDR の本文をこの作業では読んでいない (取得していない) ので作っていません (未確認)。値が CDR と違う場合も、違いの理由は推測しません。
