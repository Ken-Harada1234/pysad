# errors: KEK-LS の誤差入りの計算・補正

誤差 (位置ずれ・強さ・回転) を入れたリングの COD・分散・チューンのずれと、COD の補正。元の `correct_cod.sad`・`print_*.sad`・`write_cod_response.sad` を、短く動く形に整えた。各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 点。

## 誤差の与え方と補正の流れ

- 乱数の種: `SeedRandom[seed]`。例は seeds = {17, 111, 3} (補正の例は {17, 111})。元の台本は `SeedRandom[187]` と `Random[1000]` から作った 115 個の種。`GCUT = 1` (σ 単位の打ち切り。元の台本と同じ)。
- 位置ずれ ae = 40 μm: 四極・六極の成分ごとの `LINE["DX"/"DY", 成分名] = ae GaussRandom[]` (`modules/modules_errors.n`)。元の `DELX`・`DELY` コマンドは、この K64 ビルドでは `cell calc` の閉軌道に効かない (DUMP には値が出るだけ。gs・pysad とも) ので使わない。
- 強さの誤差 fe = 4e-4 (`DK`)、回転の誤差 te = 100 μrad (`DTHETA`): FFS のコマンド (偏向磁石・四極・六極に)。
- 成分ごとの値 (DX・DY・K0) は、`USE RING NOEXPAND` と `cell calc noexpand` でないと同名の要素の値に揃えられてしまう (`docs/11_sad_quirks.md` A56)。HD・VD (補正磁石) は同じ名前の繰り返しなので K0 は NOEXPAND が必要。
- 補正の流れ: (1) COD 応答行列を解析式で作る (`KEKLS_cod_response` で SAD の COD と照合)、(2) 誤差を与える、(3) 応答行列の特異値分解の疑似逆行列で K0 を求め、0.5 倍ずつ繰り返しかける、(4) COD・分散の rms と最大値、チューンの変化を表示する。

## 一覧 (gs と pysad の一致: ○ 一致、△ 差あり、× 不一致または未対応)

| 例 | 内容 | 元のスクリプト | 出力 | 所要時間 | gs と pysad |
|---|---|---|---|---|---|
| [`KEKLS_errors_cod`](KEKLS_errors_cod.md) | 誤差の種類ごとの COD・分散・チューン | print_all_error.sad ほか print_* | COD・分散・チューンのずれ | gs 1 s / pysad 6 s | △ |
| [`KEKLS_correct_cod`](KEKLS_correct_cod.md) | COD・分散の補正 | correct_cod.sad | 補正前後の COD・分散 | gs 1 s / pysad 37 s | △ |
| [`KEKLS_cod_response`](KEKLS_cod_response.md) | COD 応答行列と解析式の照合 | write_cod_response.sad | 応答行列の代表値・解析式との差 | gs 1 s / pysad 9 s | ○ |

## 元の設定と例の設定

| 台本 | 元の設定 | 例の設定 | 完全な設定にするには | 見込みの所要時間 (完全) |
|---|---|---|---|---|
| `KEKLS_errors_cod` | 種 115 個 | 種 3 個 × 誤差 4 通り | `seeds` を 115 個に、`nseed = 115` | 3 つの種 (各 4 通り) で gs 1 秒、pysad 6 秒 → 115 個で gs 約 40 秒、pysad 約 4 分 |
| `KEKLS_correct_cod` | 種 115 個 (k = 1〜115)、補正 40 回、毎回 DA (5000 周) | 種 2 個、補正 6 回、DA なし | `seeds`・`niter = 40`・`doDA = 1; DAturns = 5000` | 種 2 個で gs 1 秒、pysad 37 秒 (応答行列の準備に約 17 秒) → 115 個で gs 約 1 分、pysad 約 20 分。DA は未測定 |
| `KEKLS_cod_response` | 全ステアリング分を fort.92 に書く | 代表値と照合 3 個 | `j0` を 1〜ICX | pysad 約 9 秒 |

## 元のスクリプトのうち取り込んでいないもの

- `correct_cod_best`・`_worst`・`_noSX`・`_GCUT1`・`_v03_70`: `correct_cod.sad` の設定違い (monitor 名 `BPM.*`、DA の範囲、誤差の大きさ ae = 5e-5・fe = 5e-4、六極を切る `S* 0`、チューンの窓 18〜18.5 など)。例の変数 (`ae`・`fe`・`te`・`doDA`、チューンの窓) を変えれば同じ形になる。
- `print_*.sad` 11 本: 同じ台本に表示だけ違う (`print_IEX_after`・`print_IEY_after`・`print_rep_after` は補正の固有値の数 IEX・IEY、`print_residual_cod`・`print_cod_kick_ratio` は補正後の COD・キック、`print_emittance`・`print_solution_ok`・`print_all_error_after`・`print_cod_form` は表示の形式)。誤差だけの `print_all_error`・`print_dk_effect`・`print_dt_error`・`print_cod_magnitude` の流れは `KEKLS_errors_cod`、補正は `KEKLS_correct_cod` に統合した。
- `modules_correctoptics.n`・`modules_write_correctdata.n`: `correct_cod.sad` は読み込んでいない (コメントアウトされた `Get` のみ)。光学補正 (ベータ関数の補正) の関数集で、v03_68 の補正には使っていない。`modules_correctoptics_pf.n` は PF 用。
