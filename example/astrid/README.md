# ASTRID2 の例集 (0.58 GeV)

ASTRID2 のラティスを使った、SAD の基本的な使い方と計算の例です。KEK の SAD チュートリアルの話題を、ASTRID2 のラティスで作り直しました。各例は `.sad`、`.gs.txt`、`.pysad.txt`、`.md` の 4 点です。式と数値の照合は各 `.md` にあります。

## 索引

| フォルダ | 例 | 内容 | gs と pysad |
|---|---|---|---|
| [`lattice/`](lattice/README.md) | ラティス | ASTRID2 のラティス | |
| [`tune_survey/`](tune_survey/README.md) | [`tune_survey_ncell_ideal_ASTRID2`](tune_survey/tune_survey_ncell_ideal_ASTRID2.md) | ノーマルセルのチューンサーベイ (誤差なし、縮小版) | △ (チューン・εx 一致、DA は ±1〜2 σ。4 点・200 周に縮小) |
| [`tune_survey/`](tune_survey/README.md) | [`tune_survey_ncell_smallerror_ASTRID2`](tune_survey/tune_survey_ncell_smallerror_ASTRID2.md) | ノーマルセルのチューンサーベイ (小さい誤差あり、縮小版) | △ (乱数列が違う。チューン・DA の傾向は一致) |
| [`tune_radiation/`](tune_radiation/README.md) | [`amp_dep_ts_ASTRID2`](tune_radiation/amp_dep_ts_ASTRID2.md) | 振幅依存チューンシフト | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`calc_RadInt_ASTRID2`](tune_radiation/calc_RadInt_ASTRID2.md) | 放射積分 I1〜I5 と、エミッタンスの式による検算 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`mom_dep_ts_ASTRID2`](tune_radiation/mom_dep_ts_ASTRID2.md) | 運動量依存チューンシフト (チューンと運動量ずれの関係) | ○ |
| [`tracking/`](tracking/README.md) | [`phase_x_ASTRID2`](tracking/phase_x_ASTRID2.md) | 位相空間の軌跡 (水平 x-x') | ○ (失われた粒子の座標の符号だけ違う) |
| [`insertion_device/`](insertion_device/README.md) | [`MPW_v01_02_id_ASTRID2`](insertion_device/MPW_v01_02_id_ASTRID2.md) | 挿入光源 MPW の線形シミュレーション | ○ |
| [`lattice_building/`](lattice_building/README.md) | [`SetElement_BeamLine_Join_ASTRID2`](lattice_building/SetElement_BeamLine_Join_ASTRID2.md) | 周期性のあるラティスの作り方 (六極だけ全セル独立にする) | ○ |
