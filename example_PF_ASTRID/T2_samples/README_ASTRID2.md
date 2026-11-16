# T2_samples: ASTRID2 (0.58 GeV) の例 (KEK SAD チュートリアル Tutorial2、もう一人の担当の PF の例とは別の表)

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
