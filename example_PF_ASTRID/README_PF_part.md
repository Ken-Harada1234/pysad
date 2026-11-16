# README の全体表: PF 側の行 (T1・T2・lattice)

(全体の `README.md` に結合する前の PF の分。ASTRID2・T5・T6 の行は別の担当。gs と pysad の一致: ○ 一致 / △ 差あり / × 動かない・未対応。)

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

各行のファイルは `.sad` / `.gs.txt` / `.pysad.txt` / `.md` の 4 つ組。図 (PNG) は同じフォルダ。計算を短縮した例 (DA の周回数・方向、amp_dep の点数、位相空間の周回数) は各 `.md` に書いてある。
