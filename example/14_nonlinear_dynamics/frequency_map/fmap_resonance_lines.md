# 周波数マップの共鳴線の予想 (3 次 + 4 次 + 非共鳴の補正) と判別力の評価

台本ごとの式・表・原因は `FMAP_summary.md` §13、導出は `docs/25` §12、試験は `pysad/tests/script_tests/test_v23_S3_fmappred.py` (12 件)、コードは `pysad/fmappred.py` (`ResonancePredictorS3`・評価の関数)と `pysad/fmapres.py` (`from_lattice(model='s3')`)。
ラティス: FODO (`example/00_lattice/FODO_ring.sad`、QF の K1 = 0.697・SF・SD の K2 = 4, −6、3 GeV、105.6 m)、ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m)、ASTRID2 MPW (`ASTRID2_v01_02_id.sad` の `RINGMPW`)、KEK-LS (`pysad/tests/kek-ls_v03_12.lat`、大きなリング A。ESRF 型の HMBA、3 GeV、603 m)。
条件: v2、周波数マップは保存済みの表 (1024 周、RFSW・RAD OFF、LOSSAMPL 0.1 m、41×31 の x–y)、`adts(fourth='r3', optics='calc')`。**当てはめの項は入れない** (摂動論の 3・4 次の h、非共鳴の補正、線形の ADTS、振り子の式だけ)。

| ファイル | 内容 |
|---|---|
| `fmap_predict.py` | `astrid2|kekls|fodo|mpw`: 3 つの予想 (1 次の RDT だけ・4 次の h も入れた・非共鳴の補正つき) の判別力 (cov・面積・対照・AUC・順位相関・上位 10%・振幅だけの AUC)、信頼度、最寄りの共鳴線の次数、調和ごと、チューンの固着、ν(J) の誤差 |
| `fmap_predict_<ring>.out.txt`・`.png` | 出力 (全数値) と図 (D + 1 次だけの予想・補正つきの予想の包絡) |

## 結果 (D > −5 を正例、解析の ν(J)、補正つき = 3 次 + 4 次 + 非共鳴の補正)

| リング | 1 次だけ cov/area | 4 次も cov/area | 補正つき cov/area | 対照 | 補正つきの AUC | 振幅だけの AUC | 偏相関 |
|---|---|---|---|---|---|---|---|
| FODO | 1.00/0.35 | 1.00/0.35 | 1.00/0.35 | 0.37±0.36 | 0.94 | 0.42 | +0.62 |
| ASTRID2 | 0.04/0.02 | 0.41/0.25 | 0.45/0.25 | 0.25±0.17 | 0.68 | 0.79 | −0.19 |
| KEK-LS | 0.32/0.17 | 0.69/0.56 | 0.66/0.49 | 0.65±0.10 | 0.60 | 0.77 | +0.09 |
| ASTRID2 MPW | 0.21/0.11 | 0.40/0.16 | 0.33/0.14 | (不適切) | 0.79 | 0.73 | +0.46 |

結論: 予想できるのは小振幅・孤立した低次の共鳴 (FODO の 3νx、MPW の生き残り領域)。ASTRID2・KEK-LS は改善しなかった (理由: 次数 5 以上・どの線にも近くない点が D の大きい点の約半分、島の中は規則的で D の大きい点は縁・重なり、KEK-LS は大振幅)。島の位置・幅は当たる (2νx−2νy の島の中の点のチューンは線に固着)。

## コマンド

```
set PYTHONUTF8=1
py -3.12 example/14_nonlinear_dynamics/frequency_map/fmap_predict.py astrid2   (各 数秒。保存済みの表を読む。出力を上書きする)
py -3.12 -m pytest pysad/tests/script_tests/test_v23_S3_fmappred.py -q
```

Python から: `ResonancePredictor.from_lattice(lat, loc, rows, model='s3')`、`pr.predict(x, y)`・`pr.warnings_for(m)`、`fmappred.evaluate_prediction(r, m)`。

## 使えない範囲・未確認

δ ≠ 0、次数 5 以上の駆動項、J² の項、PF・SP8 (大きなリング B)、誤差つき・RFSW ON・放射つき、AUC の誤差の見積もり (docs/25 §12.7)。
