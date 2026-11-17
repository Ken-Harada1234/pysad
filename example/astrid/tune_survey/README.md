# ASTRID2 / tune_survey: チューンサーベイ

各例は、`NAME.sad` (台本。ラティスを先頭に含む)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明) の 4 点です。一致の凡例は ○ 一致、△ 差あり、× 未対応です。

| 例 | 内容 | gs と pysad |
|---|---|---|
| [`tune_survey_ncell_ideal_ASTRID2`](tune_survey_ncell_ideal_ASTRID2.md) | ノーマルセルのチューンサーベイ (誤差なし、縮小版) | △ (チューン・εx 一致、DA は ±1〜2 σ。4 点・200 周に縮小) |
| [`tune_survey_ncell_smallerror_ASTRID2`](tune_survey_ncell_smallerror_ASTRID2.md) | ノーマルセルのチューンサーベイ (小さい誤差あり、縮小版) | △ (乱数列が違う。チューン・DA の傾向は一致) |
