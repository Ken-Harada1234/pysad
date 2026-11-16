# T2_samples: サンプル集 (KEK SAD チュートリアル Tutorial2) の PF と ASTRID2 の例

各例は `.sad`・`.gs.txt`・`.pysad.txt`・`.md` の 4 つ組 (図 PNG は同じフォルダ)。`*_PF` は PF (2.5 GeV) の例、`*_ASTRID2` は ASTRID2 (0.58 GeV) の例。

## PF (2.5 GeV) の話題

元: KEK SAD チュートリアル「SAD サンプル集」(Tutorial2) の PF のスクリプト。PF 2.5 GeV リング (`../lattice/v9_16_asym.sad`、入射バンプだけ `v9_16_inj_noSF.sad`、周長 187.035 m、チューン 9.5975 / 5.2790)。ASTRID2 の話題 (チューンサーベイ・MPW など) は別の担当の `*_ASTRID2.*` (このフォルダの別のファイル)。

| ファイル (.sad/.gs.txt/.pysad.txt/.md) | サンプル集の節 | 内容 | gs と pysad |
|---|---|---|---|
| `calc_respx_PF` | calc_respx / calc_resp_analy | 軌道応答 (K0 → モニター) と解析式、符号と RFSW の確認 | ○ (応答は解析式と 0.08% で一致) |
| `calc_dispersion_PF` | calc_dispersion | 分散 + 差分と α_c の確認 | ○ |
| `calc_chromaticity_PF` | calc_chromaticity | `PolynomialFit` の 6 次、自然色収差、四極・六極の和の式 | ○ |
| `FFS_FIT_injbump_PF` | ローカルバンプ | 入射バンプ (`ins`、`fit INJ DX`、`free KK* K0`) と解析式 | ○ |
| `fn_DynamicApertureSurvey_DAsv_PF` | sample_DAsv | DA サーベイ (x・y 方向、運動量 9 点、2000 周) | △ (y の 1 行だけ gs 11 / pysad 13) |
| `fn_DynamicApertureSurvey_da_xy_PF` | da_xy | X-Y 面の DA (24 方向、1000 周) と図 | ○ |
| `FFS_mom_dep_ts_PF` | mom_dep_ts | 運動量依存チューンシフト (六極あり・なし) と図 | ○ |
| `FFS_amp_dep_ts_PF` | amp_dep_ts | 振幅依存チューンシフト (12 点平均) と図 | ○ |
| `calc_RadInt_PF` | calc_RadInt | 放射積分、εx・U0・σE・Jx の式と数値の照合 | ○ (εx は SAD と 0.24%) |
| `phase_PF` | phase | 位相空間 (60 周) と図 | ○ |
| `FFS_draw_topdrawer_PF` | topdrawer 形式の図 | topdrawer ファイルの手書き。`draw` は gs で動かず (出力なし/止まる)、pysad は未対応 | ○ (手書き部分) / `draw` は × |

- 走らせ方: `py -3.12 tools/runpf.py T2_samples/NAME.sad` (全部を回すときは `*_PF.sad` だけを指定すること。このフォルダは ASTRID2 の担当と共有)。`.sad` は `tools/build_T2.py`、図は `tools/make_figs.py`。
- 長い計算 (DA) は 5000 周 → 1000〜2000 周、方向 161 → 24 などに短縮した (各 `.md` に書いてある)。
- 未対応・書けなかったもの: `calc_respy` (垂直の応答)、`drawmagvw3d` (CAD 出力、CAD に読み込む `.vss` は作っていない)、`modules_*` の GUI 版。

## ASTRID2 (0.58 GeV) の話題

加速器: ASTRID2 (0.58 GeV、周長 45.7044 m、6 セル、チューン 5.17999 / 2.14018)。ラティスは `../lattice/ASTRID2_v01_01.sad` (MPW の例だけ `ASTRID2_v01_02_id.sad`)。各例は 4 ファイル (`NAME.sad`・`NAME.gs.txt`・`NAME.pysad.txt`・`NAME.md`)。走らせ方 (リポジトリの直下から): `py -3.12 example/tools/run_example.py example_PF_ASTRID/T2_samples/NAME.sad`。長い計算は点数・周回数を縮めた (md に書いた)。凡例: ○ 一致 / △ 差あり (md に理由) / × gs 失敗 / × pysad 未対応。

| 例 | チュートリアルの節 | 内容 | gs と pysad |
|---|---|---|---|
| [tune_survey_ncell_ideal_ASTRID2](tune_survey_ncell_ideal_ASTRID2.md) | tune_survey_ncell_ideal | セルのチューンを FIT、色収差補正、DA (リング 200 周)。4 点に縮小。ring tune = 6×(cell tune+0.0005) | △ (チューン・εx は一致、DA は ±1〜2 σ) |
| [tune_survey_ncell_smallerror_ASTRID2](tune_survey_ncell_smallerror_ASTRID2.md) | tune_survey_ncell_smallerror | 上に小さい誤差 (5 μm、5e-5、1e-5 rad) を入れて DA の減りを見る | △ (乱数列が違う。傾向は一致) |
| [MPW_v01_02_id_ASTRID2](MPW_v01_02_id_ASTRID2.md) | 挿入光源の線形シミュレーション | MPW (K=32.7、λ=0.175 m) の ΔNY=0.109、ΔU0=1634 eV と、薄肉レンズ・Cγ の式 | ○ |
| [calc_RadInt_ASTRID2](calc_RadInt_ASTRID2.md) | calc_RadInt.sad | 放射積分 I1..I5 と εx・U0・αc・σδ の式 (Emittance[] と 2e-5〜7e-5 で一致) | ○ |
| [mom_dep_ts_ASTRID2](mom_dep_ts_ASTRID2.md) | mom_dep_ts.sad | チューン対 δ (六極 ON/OFF)、PolynomialFit、六極の寄与の式 (1.3 % 以内) | ○ |
| [amp_dep_ts_ASTRID2](amp_dep_ts_ASTRID2.md) | amp_dep_ts.sad | 振幅依存チューンシフト (ins calc)。A² 則。トラッキングの約 2 倍 (Jacobian) | ○ |
| [phase_x_ASTRID2](phase_x_ASTRID2.md) | phase_x.sad (位相空間プロット) | 6 粒子の x-x' 位相空間 (300 周)。PNG。トラッキングのチューンのずれ | ○ (失われた粒子の座標の符号だけ違う) |

図: `mom_dep_ts_ASTRID2_tune_vs_dp.png`、`phase_x_ASTRID2_xpx.png` (gs の出力を matplotlib で描いた)。
作っていない ASTRID2 の話題: 軌道応答 (`calc_respx`)・分散 (`calc_dispersion`)・DA (`da_xy`・`DAsv`)・topdrawer 出力は、PF の例 (`*_PF`) で扱う。
