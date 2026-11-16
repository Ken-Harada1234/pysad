# example_PF_ASTRID (KEK SAD チュートリアルをもとにした PF と ASTRID2 の例集)

KEK の SAD チュートリアル 4 ページ (入門マニュアル、Tutorial2 サンプル集、Tutorial5 PF 入射軌道、Tutorial6 周期ラティスの作り方) の話題を、手元のラティスで書き直し、gs (SAD 本体) と pysad の両方で走らせた例集。本文を写さず、式と数値の照合を付けた。加速器は PF (2.5 GeV、周長 187.035 m、チューン 9.5975 / 5.2790) と ASTRID2 (0.58 GeV、周長 45.7044 m、チューン 5.17999 / 2.14018)。

## 構成
```
lattice/                 PF (v9_16_asym、v9_16_inj_noSF) と ASTRID2 (v01_01、v01_02_id の MPW 入り) のラティス (README_PF.md、README_ASTRID2.md)
T1_intro/                入門マニュアルの話題 (PF)
T2_samples/              サンプル集 (PF と ASTRID2)
T5_injection/            PF 入射軌道 (キッカー、関数化、GUI)
T6_lattice_building/     周期性のあるラティスの作り方 (ASTRID2)
tools/                   PF 側の .sad を作る・走らせる道具
```
各フォルダに README.md (そのフォルダの表。`lattice/` は `README_PF.md` と `README_ASTRID2.md`)。各例は `NAME.sad` (台本、ラティスを先頭に貼った自己完結のファイル)・`NAME.gs.txt` (gs の画面表示)・`NAME.pysad.txt` (pysad の画面表示)・`NAME.md` (説明、式と数値の照合表、限界)。

走らせ方 (リポジトリの直下 `Claude outputs` から): `py -3.12 example/tools/run_example.py example_PF_ASTRID/T2_samples/NAME.sad` (gs は WSL、pysad は `PYSAD_SAD_PACKAGES` を `Packages` に向けて走る)。凡例: ○ 一致 (表示の桁・書式の差は除く) / △ 差あり (各 `.md` に理由) / × 動かない・未対応。

## 全体の表

| 例 (フォルダ/ファイル名) | 元のチュートリアルの節 | 加速器 | gs と pysad の一致 |
|---|---|---|---|
| `lattice/v9_16_asym.sad`、`v9_16_inj_noSF.sad` | サンプル集「PF ラティス台本」、Tutorial5 の入射ラティス | PF 2.5 GeV | (ラティス。T1 の `FFS_cell_calc_emit_PF` で周長・チューン・εx 全桁一致) |
| `T1_intro/FFS_cell_calc_emit_PF` | 入門マニュアル その 1・2 (実行、cell・calc・emit) | PF 2.5 GeV | ○ |
| `T1_intro/FFS_disp_type_PF` | 入門マニュアル その 2 (disp・type) | PF 2.5 GeV | ○ |
| `T1_intro/fn_Twiss_LINE_PF` | 入門マニュアル その 2 (Twiss・LINE) | PF 2.5 GeV | ○ |
| `T1_intro/fn_PolynomialFit_PF` | 入門マニュアル その 3 (チューン対運動量、PolynomialFit) | PF 2.5 GeV | ○ (SAD の係数の癖を確認) |
| `T1_intro/FFS_draw_OpticsPlot_PF` | 入門マニュアル その 4 (draw、OpticsPlot) | PF 2.5 GeV | ○ (draw は gs 出力なし・pysad 未対応) |
| `T1_intro/FFS_FIT_chrom_PF` | 入門マニュアル その 5 (色収差補正 FIT) | PF 2.5 GeV | ○ |
| `T1_intro/FFS_FIT_modules_PF` | 入門マニュアル その 6 (マッチングと色収差補正のモジュール化) | PF 2.5 GeV (標準セル) | ○ |
| `T1_intro/fn_DynamicApertureSurvey_PF` | 入門マニュアル その 7 (DA サーベイ) | PF 2.5 GeV | ○ |
| `T1_intro/fn_TrackParticles_PF` | 入門マニュアル その 8 (トラッキング、転送行列) | PF 2.5 GeV | ○ |
| `T2_samples/calc_respx_PF` | サンプル集 calc_respx・calc_resp_analy (軌道応答と解析式) | PF 2.5 GeV | ○ |
| `T2_samples/calc_dispersion_PF` | サンプル集 calc_dispersion | PF 2.5 GeV | ○ |
| `T2_samples/calc_chromaticity_PF` | サンプル集 calc_chromaticity | PF 2.5 GeV | ○ |
| `T2_samples/FFS_FIT_injbump_PF` | サンプル集 ローカルバンプ (入射バンプ) | PF 2.5 GeV (`v9_16_inj_noSF`) | ○ |
| `T2_samples/fn_DynamicApertureSurvey_DAsv_PF` | サンプル集 sample_DAsv | PF 2.5 GeV | △ (y の 1 行だけ gs 11 / pysad 13) |
| `T2_samples/fn_DynamicApertureSurvey_da_xy_PF` | サンプル集 da_xy | PF 2.5 GeV | ○ |
| `T2_samples/FFS_mom_dep_ts_PF` | サンプル集 mom_dep_ts | PF 2.5 GeV | ○ |
| `T2_samples/FFS_amp_dep_ts_PF` | サンプル集 amp_dep_ts | PF 2.5 GeV | ○ |
| `T2_samples/calc_RadInt_PF` | サンプル集 calc_RadInt (放射積分) | PF 2.5 GeV | ○ |
| `T2_samples/phase_PF` | サンプル集 phase (位相空間) | PF 2.5 GeV | ○ |
| `T2_samples/FFS_draw_topdrawer_PF` | サンプル集 topdrawer 形式の図 | PF 2.5 GeV | ○ (手書き部分)、`draw` コマンドは × (gs は出力なし/止まる、pysad は未対応) |
| `lattice/ASTRID2_v01_01.sad`、`ASTRID2_v01_02_id.sad` | サンプル集「tune_survey_ncell_ideal」の ASTRID2 ラティス、「挿入光源の線形シミュレーション」の MPW 入りラティス | ASTRID2 0.58 GeV | (ラティス。`T2_samples/calc_RadInt_ASTRID2`・`MPW_v01_02_id_ASTRID2` で周長 45.7044 m・チューン・εx・U0 が gs と pysad で全桁一致) |
| `T2_samples/tune_survey_ncell_ideal_ASTRID2` | サンプル集 tune_survey_ncell_ideal (チューンサーベイ、誤差なし) | ASTRID2 0.58 GeV | △ (チューン・εx 一致、DA は ±1〜2 σ。4 点・200 周に縮小) |
| `T2_samples/tune_survey_ncell_smallerror_ASTRID2` | サンプル集 tune_survey_ncell_smallerror (小さい誤差あり) | ASTRID2 0.58 GeV | △ (乱数列が違う。チューン・DA の傾向は一致) |
| `T2_samples/MPW_v01_02_id_ASTRID2` | サンプル集 挿入光源の線形シミュレーション (MPW v01_02_id) | ASTRID2 0.58 GeV | ○ |
| `T2_samples/calc_RadInt_ASTRID2` | サンプル集 calc_RadInt (放射積分) | ASTRID2 0.58 GeV | ○ |
| `T2_samples/mom_dep_ts_ASTRID2` | サンプル集 mom_dep_ts (運動量依存チューンシフト) | ASTRID2 0.58 GeV | ○ |
| `T2_samples/amp_dep_ts_ASTRID2` | サンプル集 amp_dep_ts (振幅依存チューンシフト) | ASTRID2 0.58 GeV | ○ |
| `T2_samples/phase_x_ASTRID2` | サンプル集 位相空間プロット (phase_x) | ASTRID2 0.58 GeV | ○ (失われた粒子の座標の符号だけ違う) |
| `T5_injection/calc_traj_01_PF` | Tutorial5 計算部分の準備 (キッカー・入射軌道) | PF 2.5 GeV (`v9_16_inj_noSF`) | ○ |
| `T5_injection/calc_traj_02_PF` | Tutorial5 モジュール化 (関数化、CS 不変量) | PF 2.5 GeV | ○ |
| `T5_injection/calc_traj_03_PF` | Tutorial5 GUI 入出力化 (KBMainFrame・Canvas・OpticsPlot) | PF 2.5 GeV | ○ (計算も窓の項目数 114・115 も一致。2026-11-17 に修正) |
| `T5_injection/calc_traj_03b_PF` | Tutorial5 グラフ表示 (pysad でも動く ListPlot 版、pysad の窓の PNG) | PF 2.5 GeV | ○ (Canvas の項目数 382 → 384 が一致) |
| `T6_lattice_building/SetElement_BeamLine_Join_ASTRID2` | Tutorial6 周期性のあるラティスの作り方 (SetElement・BeamLine・Join) | ASTRID2 0.58 GeV (教材は別の 16 セルのリング) | ○ |

## 共通の注意
- 長い計算 (DA の周回数、チューンサーベイの点数、位相空間の周回数) は 1 本 60 秒以内に縮めた。縮めた量は各 `.md` に書いてある。
- ASTRID2 のチューンサーベイは、現在の SAD では 1 セル (RF なし) の emit が閉軌道を見つけられないため、DA の測定だけリング (6 セル + RF) で行った。元のスクリプトの NCELL での DA とは値が比べられない (`T2_samples/tune_survey_ncell_ideal_ASTRID2.md`)。
- `PolynomialFit` は点の中心が 0 に近いと SAD の癖で別の点まわりの係数を返す (`T2_samples/mom_dep_ts_ASTRID2.md`、`T1_intro/fn_PolynomialFit_PF.md`)。
- pysad が動かなかった部分: `OpticsPlot` に自前のデータ列を渡す使い方と、`Canvas$Widget = c1` で自前の枠に描く書き方 (`T5_injection/calc_traj_03_PF.md`)。`ListPlot` + `CanvasDrawer` は動く。gs では `calc_traj_03` の窓も 114 項目が描けた。
- 乱数を使う例 (`tune_survey_ncell_smallerror_ASTRID2`) は `SeedRandom` を固定したが、gs と pysad で乱数列が違うので値の一致は求めていない。
- KEK 公式サイトの参照は、今回の依頼で指定された 4 ページに限った。ラティスとスクリプトの出所は、ユーザーの手元のアーカイブのコピー (各 `lattice/README_*.md`)。
