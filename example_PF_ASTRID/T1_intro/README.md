# T1_intro: SAD 入門マニュアルの話題を PF で
元: KEK SAD チュートリアル「SAD の入門マニュアル」(その 1〜8)。NSLS X-ring の例を、PF 2.5 GeV リング (`../lattice/v9_16_asym.sad`、周長 187.035 m、チューン 9.5975 / 5.2790) に置き換えた。説明は自分の言葉で書き直してある。

| ファイル (.sad/.gs.txt/.pysad.txt/.md) | 入門マニュアルの節 | 内容 | gs と pysad |
|---|---|---|---|
| `FFS_cell_calc_emit_PF` | その 1・2 | 実行、`cell calc emit`、`SAVE`、周長・チューン・エミッタンス | ○ |
| `FFS_disp_type_PF` | その 2 | `disp`・`type` | ○ |
| `fn_Twiss_LINE_PF` | その 2 | `Twiss[]`・`LINE[]` | ○ |
| `fn_PolynomialFit_PF` | その 3 | チューン対運動量、`DP0`、`PolynomialFit` (SAD の係数の癖あり) | ○ |
| `FFS_draw_OpticsPlot_PF` | その 4 | `draw`、`OpticsPlot` (図の数値と PNG) | ○ (draw は gs 無出力・pysad 未対応) |
| `FFS_FIT_chrom_PF` | その 5 | 色収差補正の六極フィット (設計値が戻る) | ○ |
| `FFS_FIT_modules_PF` | その 6 | `ChromCorrect[]`・`MatchingNcellTune[]` のモジュール化 | ○ |
| `fn_DynamicApertureSurvey_PF` | その 7 | DA サーベイ (1000 周) | ○ (表示記号に 1〜2 文字の差) |
| `fn_TrackParticles_PF` | その 8 | 四極 1 個の転送行列 | ○ |

- 走らせ方: `py -3.12 tools/runpf.py T1_intro/NAME.sad` (gs は WSL、pysad は `PYSAD_SAD_PACKAGES` つき)。`.sad` は `tools/build_T1.py` が作る。
- 追跡の例 (DA、TrackParticles) は `OFF RAD FLUC; ON COD RFSW RADCOD`。それ以外はチュートリアルと同じ `ON RAD FLUC COD RFSW RADCOD`。
- 入門マニュアルで扱っていない/省いたもの: ファイル読み書き (`OpenWrite`、`Write`、`Close`。出力は画面で見る)、`FitValue` の上下限つきマッチング。
