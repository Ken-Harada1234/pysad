# design: KEK-LS の設計・マッチング

v03_68 のセル・リング・入射セルのマッチングと、セルの縮小・セル数の変更・チューンサーベイ。元のスクリプトを `kek_ls_v03_68.sad` で動くように整えた。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点 (`NAME.time.txt` は所要時間の記録)。

## 流れ

1. **セルのチューン** `KEKLS_match_cell_tune`: 標準セル `ARCA` の位相進みを (48.58, 17.62)/20 に合わせる (四極の K1)。
2. **六極** `KEKLS_match_chrom`: リング `RINGC` の六極 `SDH`・`SFH` でクロマティシティを 0 にする。
3. **入射セル** `KEKLS_match_injcell`: `ARCINJH` の四極 `QS0`〜`QS5` で、入射点の位相・ツイスを合わせる。
4. **リングを閉じる** `KEKLS_match_ring_geo`: 入射セルの長直線部 `DRID1INJ` で、ジオメトリを閉じる (GX = GY = 0)。
5. 設計の変形: `KEKLS_scale_mini` (全体の一様な縮小)、`KEKLS_scale_ncell` (セル数の変更)。
6. セルのチューンサーベイ `KEKLS_tune_survey_cell`、再現性の診断 `KEKLS_match_repeat`。

## 一覧 (gs と pysad の一致: ○ 一致、△ 差あり、× 不一致または未対応)

| 例 | 内容 | 元のスクリプト | 出力 | 所要時間 | gs と pysad |
|---|---|---|---|---|---|
| [`KEKLS_match_cell_tune`](KEKLS_match_cell_tune.md) | セルのチューンのマッチング | matching2_01 前半 | セルの位相進み・四極の K1 | gs 1 s / pysad 2 s | △ |
| [`KEKLS_match_chrom`](KEKLS_match_chrom.md) | 六極のクロマティシティ調整 | matching2_01 後半 | ξ・六極の K2 | gs 1 s / pysad 9 s | △ |
| [`KEKLS_match_injcell`](KEKLS_match_injcell.md) | 入射セルのマッチング | matching2_03_INJ | 入射点の Twiss・QS0〜QS5 の K1 | gs 1 s / pysad 2 s | ○ |
| [`KEKLS_match_ring_geo`](KEKLS_match_ring_geo.md) | 入射セルを含めてリングを閉じる (ジオメトリ) | matching2_02_geo | DRID1INJ・GX・GY | gs 1 s / pysad 2 s | ○ |
| [`KEKLS_match_repeat`](KEKLS_match_repeat.md) | マッチングの再現性の診断 | (新規) | K1・K2 (3 回の繰り返し) | gs 1 s / pysad 20 s | △ |
| [`KEKLS_scale_mini`](KEKLS_scale_mini.md) | セルの一様な縮小 | matching23_v03_35_mini | セル長・位相進み・K1 | gs 1 s / pysad 4 s | ○ |
| [`KEKLS_scale_ncell`](KEKLS_scale_ncell.md) | セル数を変える | matching21_v03_35 | j セルのリングの εx・チューン・周長 | gs 1 s / pysad 3 s | ○ |
| [`KEKLS_tune_survey_cell`](KEKLS_tune_survey_cell.md) | セルのチューンサーベイ | matching05 + modules_tune_sv_v03_68.n | ξ・εx・六極の強さの 25 点の地図 | gs 10 s / pysad 160 s | △ |

マッチング後のラティスの保存: 各例は、結果を画面 (`type`・`Print`) に出す。元の台本の `SAVE ALL` はラティス全体を画面に書き出す (gs・pysad とも) ので、例では K1・K2 の値を `Print` する形にした。`SAVE ALL` の出力をファイルに取るには、台本の `SAVE ALL;` を FFS の最後に書いて標準出力を保存する (`FFS["SAVE ALL", ファイル]` の形はファイルが空になるので使えない)。

## 元の設定と例の設定

| 台本 | 元の設定 | 例の設定 | 完全な設定にするには | 見込みの所要時間 (完全) |
|---|---|---|---|---|
| `KEKLS_match_cell_tune` | `exmkc` 0.004 → 0 (0.001 刻み 5 点)、NPARA 32/16 | 0.004 → 0 (0.002 刻み 3 点)、NPARA 8 | `{exmkc,0.004,0.0,-0.001}` | gs 1 秒、pysad 数秒 |
| `KEKLS_match_injcell` | マッチングを 50 回繰り返す | 2 回 (1 回で収束) | `iterations = 50` | gs 数秒 |
| `KEKLS_scale_mini` | factor 1 → 0.5 を 0.01 刻み (51 点)、NPARA 8 | 0.1 刻み (6 点) | `{factor,1,0.5,-0.01}` | gs 約 10 秒、pysad 約 1 分 |
| `KEKLS_scale_ncell` | j = 20 → 2 (19 点)、NPARA 8 | j = 20, 18, 16, 12, 10 (5 点) | `names`・`ncs` を増やす | gs 数秒 |
| `KEKLS_tune_survey_cell` | tunex 46.05〜46.55 (6 点) × tuney 16.05〜19 (30 点) = 180 点、各点に DA (5000 周) | 48.30〜48.90 (5 点) × 17.30〜17.90 (5 点) = 25 点、DA なし | `tx0 = 46.05; dtx = 0.1; ntx = 6; ty0 = 16.05; dty = 0.1; nty = 30;` + `DAsv[]` | 25 点で gs 10 秒・pysad 約 160 秒から、DA なしで gs 約 1.2 分・pysad 約 19 分。DA は未測定 |

(見込みは、例の設定の所要時間からの比例の外挿。DA の時間は測っていない。元のスクリプトのファイルの更新日時は 2016 年で、当時の所要時間の記録は元のフォルダに無い。)

## 元の台本の問題と直したもの

- `matching2_02_geo`: `USE RINGI` (v03_68 に無い) → `RING`。
- `matching2_01`・`modules_tune_sv_v03_68.n`・`matching23`・`matching21`: `free BPI6 K1`・`free BPI7 K1` は v03_68 の名前 `BPI6H`・`BPI7H` に一致せず効かない。`matching23`・`matching21`・`tune_survey_cell` の例は `BPI6H`・`BPI7H` に直した (`freebend`)。`match_cell_tune` は元の挙動 (効かない) を既定に残した。
- `matching2_01`・`matching21`: FitFunction が `EMITX` を使うが、閉じていない ARCA では EMITX は NaN (`matching21` は `1^-9` の書き間違いも)。gs では何も効かず、pysad では fit が動かなくなる。この 2 本は FitFunction を除いた。
- `matching23`・`modules_tune_sv_v03_68.n`: v03_35 の `MKSFA2`・`SDA`・`SFA`・`QFA4A/B`・`ARCA1`・目標値が残っていた → v03_68 の名前・設計値に直した。
- `matching05`: `modules_tune_sv_v03_35.n` を読む。v03_68 用の `modules_tune_sv_v03_68.n` を読む台本は無い。
- 元の台本の `WriteGraphs[]` (`draw` コマンド)・`datas.txt` への書き出し・`SAVE ALL`・`stop; stop;` の後ろの旧い試行は、例では除いた (描画は PNG に別に、数表は画面に)。

## 要確認 (pysad)

- `FitFunction` を `fit` 条件と併用する `MatchingARCA` (`KEKLS_match_repeat`・`KEKLS_tune_survey_cell`・`KEKLS_scan_target`) で、pysad は gs と違う解に落ち、同じ入力の繰り返しでも結果が変わる (1 回目と 2 回目以降が違う)。NaN の FitFunction では fit が動かない。pysad のコードは変えていない (追跡に合わせる補正は足さない)。
- `GaussRandom[]` の `GCUT = 0` の扱い (gs は ≈ 0、pysad は打ち切りなし) は `../errors/` で確認 (台本は `GCUT = 1`)。
