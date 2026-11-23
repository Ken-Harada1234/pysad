# example の旧フォルダ・旧ファイル名から新フォルダ・新ファイル名への対応

`example/` の版名フォルダ (`v2_3`・`v2_4`・`v2_B`・`v2_C`・`v2_next`) を、機能別のフォルダへ移した。この文書は、その前後の対応表である。パスはすべて `example/` からの相対。`docs/history/` と `CHANGELOG.md` に出てくる旧いパスは、この表で引く。

## 1. フォルダ・接頭辞ごとの対応

旧ファイル名の接頭辞 (`Q1_`・`R3_` など。作業の符号) は、新しい名前では機能名に変わった。接頭辞ごとに、移動先のフォルダと新しい名前の例を示す。

| 旧フォルダ | 旧接頭辞 | 件数 | 新フォルダ | 新ファイル名の例 |
|---|---|---|---|---|
| `v2_3/` | `ADTS_` | 2 | `14_nonlinear_dynamics/adts_rdt/` (2) | `ADTS_summary.md`、`ADTS_summary_rdt_check.py` |
| `v2_3/` | `FMAP_` | 1 | `14_nonlinear_dynamics/frequency_map/` (1) | `FMAP_summary.md` |
| `v2_3/` | `L2_` | 62 | `15_lma_touschek/lma/` (55)、`15_lma_touschek/touschek_elegant/` (7) | `lma_reliability.md`、`lma_cell.py` |
| `v2_3/` | `Q1_` | 38 | `14_nonlinear_dynamics/frequency_map/` (20)、`15_lma_touschek/lma/` (12)、`(削除・分割)/` (6) | `fmap_basic.md`、`fmap_astrid2.out.txt` |
| `v2_3/` | `Q2_` | 17 | `14_nonlinear_dynamics/adts_rdt/` (17) | `rdt_adts_basic.md`、`adts_gs.py` |
| `v2_3/` | `R1_` | 9 | `14_nonlinear_dynamics/adts_rdt/` (9) | `adts_gs_norfsw.md`、`adts_gs_astrid_norfsw.txt` |
| `v2_3/` | `R2_` | 32 | `14_nonlinear_dynamics/adts_rdt/` (32) | `adts_direct.md`、`adts_ablate_astrid.pysad.txt` |
| `v2_3/` | `R3_` | 23 | `14_nonlinear_dynamics/adts_rdt/` (23) | `adts_bend_cubic.md`、`adts_bend.png` |
| `v2_3/` | `R4_` | 36 | `14_nonlinear_dynamics/adts_rdt/` (36) | `rdt_tracking.md`、`rdt_chrom_astrid.pysad.json` |
| `v2_3/` | `R5_` | 37 | `14_nonlinear_dynamics/adts_rdt/` (30)、`14_nonlinear_dynamics/frequency_map/` (7) | `rdt4_normal_form.md`、`rdt4_adts_astrid.pysad.json` |
| `v2_3/` | `R6_` | 12 | `14_nonlinear_dynamics/adts_rdt/` (12) | `adts_exact_terms.md`、`adts_exact_adts.out.txt` |
| `v2_3/` | `README.md` | 1 | `(削除・分割)/` (1) |  |
| `v2_3/` | `S1_` | 51 | `14_nonlinear_dynamics/frequency_map/` (51) | `fmap_validation.md`、`fmap_delta_astrid2_xdp.png` |
| `v2_3/` | `S2_` | 38 | `14_nonlinear_dynamics/frequency_map/` (29)、`14_nonlinear_dynamics/frequency_map/fmap_check_data/` (9) | `fmap_reliability.md`、`fmap_check_coupling.py` |
| `v2_3/` | `S3_` | 10 | `14_nonlinear_dynamics/frequency_map/` (10) | `fmap_resonance_lines.md`、`fmap_predict.py` |
| `v2_3/` | `T1_` | 44 | `15_lma_touschek/touschek_elegant/` (44) | `elegant_comparison.md`、`astrid2_L2isr.lte` |
| `v2_3/` | `T2_` | 20 | `15_lma_touschek/touschek_elegant/` (20) | `touschek_formula.md`、`touschek_formula_astrid2_div.ele` |
| `v2_3/` | `U2_` | 7 | `15_lma_touschek/touschek_elegant/` (7) | `bend_k1_definition.md`、`bend_k1_def_astrid2_chrom.out.txt` |
| `v2_3/` | `primer_` | 19 | `14_nonlinear_dynamics/primer/` (19) | `primer_adts_fit.png`、`primer_adts_fit.py` |
| `v2_4/` | `A1_` | 3 | `14_nonlinear_dynamics/normal_form/` (3) | `ampnf_vs_tracking.py`、`ampnf_vs_tracking_astrid.out.txt` |
| `v2_4/` | `A2_` | 1 | `14_nonlinear_dynamics/normal_form/` (1) | `k0_bend_check.py` |
| `v2_4/` | `A3_` | 8 | `15_lma_touschek/touschek_elegant/` (8) | `touschek_bunch_mc.py`、`touschek_bunch_mc_astrid2.out.txt` |
| `v2_4/` | `A4_` | 5 | `15_lma_touschek/lma/` (5) | `lma_support_coupl.out.txt`、`lma_support_k64.out.txt` |
| `v2_4/` | `A5_` | 5 | `14_nonlinear_dynamics/normal_form/` (5) | `nu_of_J_fmap_grid_astrid2.out.txt`、`nu_of_J_fmap_grid_astrid2_D8.out.txt` |
| `v2_4/` | `A7_` | 2 | `14_nonlinear_dynamics/normal_form/` (2) | `ampnf_kekls.py`、`ampnf_kekls.out.txt` |
| `v2_4/` | `F1_` | 3 | `14_nonlinear_dynamics/normal_form/` (3) | `rotated_bend_astrid.out.txt`、`rotated_bend_fodo.out.txt` |
| `v2_4/` | `G_` | 6 | `14_nonlinear_dynamics/normal_form/` (6) | `k0_independent_flow.out.txt`、`k0_independent_flow.py` |
| `v2_4/` | `H1_` | 3 | `14_nonlinear_dynamics/normal_form/` (3) | `coupled_h4_astrid.out.txt`、`coupled_h4.py` |
| `v2_4/` | `H2_` | 44 | `_deletable/lma_corrected_kekls/` (40)、`15_lma_touschek/lma/` (4) | `lma_corrected_kekls_v03_12_corr_s111.out.txt`、`lma_corrected_kekls_v03_12_corr_s17.out.txt` |
| `v2_4/` | `H2b_` | 1 | `15_lma_touschek/lma/` (1) | `lma_corrected_kekls_v03_12.template.sad` |
| `v2_4/` | `I1_` | 4 | `14_nonlinear_dynamics/normal_form/` (4) | `nu_of_J_orders_astrid.out.txt`、`nu_of_J_orders_fodo.out.txt` |
| `v2_4/` | `I2_` | 2 | `14_nonlinear_dynamics/normal_form/` (2) | `island_vs_tracking_fodo.out.txt`、`island_vs_tracking.py` |
| `v2_4/` | `README.md` | 1 | `(削除・分割)/` (1) |  |
| `v2_B/` | `B1_` | 2 | `10_math_data/` (2) | `fn_SpecialFunctions_gs_errors.pysad.txt`、`fn_SpecialFunctions_gs_errors.sad` |
| `v2_B/` | `B2_` | 2 | `08_scripting/` (2) | `syntax_package_context.pysad.txt`、`syntax_package_context.sad` |
| `v2_B/` | `B3_` | 2 | `09_io_system/` (2) | `fn_SharedMemory_Parallelize_System_Fork.pysad.txt`、`fn_SharedMemory_Parallelize_System_Fork.sad` |
| `v2_B/` | `B4_` | 2 | `01_lattice_elements/` (2) | `FFS_small_options.pysad.txt`、`FFS_small_options.sad` |
| `v2_B/` | `README.md` | 1 | `(削除・分割)/` (1) |  |
| `v2_C/` | `C1_` | 2 | `01_lattice_elements/` (2) | `elem_UND.pysad.txt`、`elem_UND.sad` |
| `v2_C/` | `C2_` | 2 | `01_lattice_elements/` (2) | `elem_PHSROT.pysad.txt`、`elem_PHSROT.sad` |
| `v2_C/` | `C2a_` | 2 | `07_spin_solenoid_special/` (2) | `elem_SOL_tilt_offset.pysad.txt`、`elem_SOL_tilt_offset.sad` |
| `v2_C/` | `C2b_` | 2 | `07_spin_solenoid_special/` (2) | `elem_SOL_RF_TRPT.pysad.txt`、`elem_SOL_RF_TRPT.sad` |
| `v2_C/` | `C2c_` | 2 | `16_spin_polarization/` (2) | `spin_tracking_sol_und.pysad.txt`、`spin_tracking_sol_und.sad` |
| `v2_C/` | `C2d_` | 2 | `04_emittance_radiation/` (2) | `flag_RADTAPER_optics.pysad.txt`、`flag_RADTAPER_optics.sad` |
| `v2_C/` | `C2m_` | 4 | `02_optics_cod/` (4) | `fn_PUTMEASURE_GETMEASURE.pysad.txt`、`fn_PUTMEASURE_GETMEASURE.sad` |
| `v2_C/` | `C2s_` | 4 | `05_tracking/` (4) | `flag_SPAC_MULT_radius.pysad.txt`、`flag_SPAC_MULT_radius.sad` |
| `v2_C/` | `C2w_` | 6 | `02_optics_cod/` (2)、`01_lattice_elements/` (2)、`06_collective/` (2) | `flag_DELCOR_rms.pysad.txt`、`flag_DELCOR_rms.sad` |
| `v2_C/` | `C3_` | 2 | `07_spin_solenoid_special/` (2) | `elem_SOL_region_elements.pysad.txt`、`elem_SOL_region_elements.sad` |
| `v2_C/` | `C4_` | 2 | `05_tracking/` (2) | `flag_LOSSMAP_RADTAPER.pysad.txt`、`flag_LOSSMAP_RADTAPER.sad` |
| `v2_C/` | `C5_` | 2 | `07_spin_solenoid_special/` (2) | `flag_RADLIGHT_orbit_field.pysad.txt`、`flag_RADLIGHT_orbit_field.sad` |
| `v2_C/` | `C6_` | 2 | `02_optics_cod/` (2) | `FFS_misc_DELCOR_WAKE_MULT_skew.pysad.txt`、`FFS_misc_DELCOR_WAKE_MULT_skew.sad` |
| `v2_C/` | `README.md` | 1 | `(削除・分割)/` (1) |  |
| `v2_next/` | `A1_` | 29 | `01_lattice_elements/` (9)、`(削除・分割)/` (8)、`04_emittance_radiation/` (6)、`02_optics_cod/` (3)、`05_tracking/` (3) | `elem_BEND_F1_adjacent.gs.txt`、`elem_BEND_F1_adjacent.pysad.txt` |
| `v2_next/` | `A2_` | 16 | `16_spin_polarization/` (7)、`04_emittance_radiation/` (6)、`01_lattice_elements/` (3) | `elem_BEND_DX_offset.gs.txt`、`elem_BEND_DX_offset.pysad.txt` |
| `v2_next/` | `B1_` | 17 | `(削除・分割)/` (4)、`04_emittance_radiation/` (3)、`07_spin_solenoid_special/` (3)、`05_tracking/` (3)、`01_lattice_elements/` (2)、`16_spin_polarization/` (2) | `elem_BEND_F1_emittance.gs.txt`、`elem_BEND_F1_emittance.pysad.txt` |
| `v2_next/` | `B2_` | 29 | `08_scripting/` (9)、`(削除・分割)/` (7)、`09_io_system/` (4)、`03_matching/` (3)、`13_packages/` (3)、`10_math_data/` (3) | `syntax_context.gs.txt`、`syntax_context.pysad.txt` |
| `v2_next/` | `B3_` | 18 | `16_spin_polarization/` (6)、`(削除・分割)/` (4)、`03_matching/` (3)、`07_spin_solenoid_special/` (3)、`05_tracking/` (2) | `FFS_FIT_geometry_quad_dx.gs.txt`、`FFS_FIT_geometry_quad_dx.pysad.txt` |
| `v2_next/` | `C3_` | 25 | `04_emittance_radiation/` (9)、`(削除・分割)/` (7)、`05_tracking/` (3)、`07_spin_solenoid_special/` (3)、`02_optics_cod/` (3) | `flag_FLUC_nparallel.gs.txt`、`flag_FLUC_nparallel.pysad.txt` |
| `v2_next/` | `D1_` | 14 | `07_spin_solenoid_special/` (6)、`04_emittance_radiation/` (3)、`(削除・分割)/` (3)、`16_spin_polarization/` (2) | `flag_NORFSW_orbit.gs.txt`、`flag_NORFSW_orbit.pysad.txt` |
| `v2_next/` | `D2_` | 21 | `13_packages/` (12)、`(削除・分割)/` (6)、`01_lattice_elements/` (3) | `fn_CorrectLocal.gs.txt`、`fn_CorrectLocal.pysad.txt` |
| `v2_next/` | `E_` | 14 | `01_lattice_elements/` (8)、`16_spin_polarization/` (4)、`(削除・分割)/` (2) | `elem_BEND_steer_AE_poles.out.txt`、`elem_BEND_steer_AE_poles.py` |
| `v2_next/` | `F2_` | 2 | `07_spin_solenoid_special/` (2) | `flag_PHOTONS_coord_tcavi.gs.txt`、`flag_PHOTONS_coord_tcavi.sad` |
| `v2_next/` | `F4_` | 2 | `07_spin_solenoid_special/` (2) | `elem_SOL_bound0_track.gs.txt`、`elem_SOL_bound0_track.sad` |
| `v2_next/` | `P1_` | 20 | `16_spin_polarization/` (17)、`(削除・分割)/` (3) | `spin_bend_vertical_angle.out.txt`、`spin_bend_vertical_angle.py` |
| `v2_next/` | `P2_` | 16 | `04_emittance_radiation/` (12)、`(削除・分割)/` (2)、`01_lattice_elements/` (2) | `elem_BEND_K1_face.out.txt`、`elem_BEND_K1_face.py` |
| `v2_next/` | `P3_` | 26 | `01_lattice_elements/` (21)、`(削除・分割)/` (5) | `elem_AE_exact_geometry.out.txt`、`elem_AE_exact_geometry.py` |
| `v2_next/` | `P4_` | 29 | `05_tracking/` (14)、`(削除・分割)/` (7)、`03_matching/` (3)、`07_spin_solenoid_special/` (3)、`16_spin_polarization/` (2) | `fn_DynamicApertureSurvey_dapwidth_ring_a.out.txt`、`fn_DynamicApertureSurvey_dapwidth_window.gs.txt` |
| `v2_next/` | `P5_` | 82 | `02_optics_cod/` (24)、`(削除・分割)/` (20)、`13_packages/` (12)、`03_matching/` (9)、`04_emittance_radiation/` (9)、`12_gui_tkinter/` (5)、`10_math_data/` (3) | `flag_DETR_astrid2.gs.txt`、`flag_DETR_astrid2.pysad.txt` |
| `v2_next/` | `README.md` | 1 | `(削除・分割)/` (1) |  |

旧フォルダ全体の行き先: `v2_3` (459 件) → `14_nonlinear_dynamics/` (adts_rdt・frequency_map・primer) と `15_lma_touschek/` (lma・touschek_elegant)。`v2_4` (88 件) → `14_nonlinear_dynamics/normal_form/` と `15_lma_touschek/`。`v2_B` (9 件)・`v2_C` (35 件)・`v2_next` (361 件) → `01_lattice_elements/` 〜 `13_packages/` の機能別フォルダと `16_spin_polarization/`。旧い索引 (`v2_*/README.md`・`v2_next/A1_README.md` など) は分割し、移動先フォルダの `README.md` と各例の `NAME.md` に再編した。

## 2. ファイル単位の対応 (全件)

種別は元の表のとおり。新パスが無いものは、削除または分割の理由を右に書いた。

| 旧パス | 新パス | 種別 | 備考 |
|---|---|---|---|
| `v2_3/ADTS_summary.md` | `14_nonlinear_dynamics/adts_rdt/ADTS_summary.md` | md |  |
| `v2_3/ADTS_summary_rdt_check.py` | `14_nonlinear_dynamics/adts_rdt/ADTS_summary_rdt_check.py` | スクリプト |  |
| `v2_3/FMAP_summary.md` | `14_nonlinear_dynamics/frequency_map/FMAP_summary.md` | md |  |
| `v2_3/L2_README.md` | `15_lma_touschek/lma/lma_reliability.md` | md |  |
| `v2_3/L2_cell.py` | `15_lma_touschek/lma/lma_cell.py` | スクリプト |  |
| `v2_3/L2_cell_astrid2.lma.txt` | `15_lma_touschek/lma/lma_cell_astrid2.lma.txt` | 出力 |  |
| `v2_3/L2_cell_astrid2.out.txt` | `15_lma_touschek/lma/lma_cell_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_cell_astrid2_tol5e-06.lma.txt` | `15_lma_touschek/lma/lma_cell_astrid2_tol5e-06.lma.txt` | 出力 |  |
| `v2_3/L2_cell_astrid2_tol5e-6.out.txt` | `15_lma_touschek/lma/lma_cell_astrid2_tol5e-6.out.txt` | 出力 |  |
| `v2_3/L2_cell_kekls.lma.txt` | `15_lma_touschek/lma/lma_cell_kekls.lma.txt` | 出力 |  |
| `v2_3/L2_cell_kekls.out.txt` | `15_lma_touschek/lma/lma_cell_kekls.out.txt` | 出力 |  |
| `v2_3/L2_check_every.py` | `15_lma_touschek/lma/lma_check_every.py` | スクリプト |  |
| `v2_3/L2_check_every_astrid2.out.txt` | `15_lma_touschek/lma/lma_check_every_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_check_every_astrid2_ce0.lma.txt` | `15_lma_touschek/lma/lma_check_every_astrid2_ce0.lma.txt` | 出力 |  |
| `v2_3/L2_check_every_astrid2_ce1.lma.txt` | `15_lma_touschek/lma/lma_check_every_astrid2_ce1.lma.txt` | 出力 |  |
| `v2_3/L2_check_every_astrid2_ce5.lma.txt` | `15_lma_touschek/lma/lma_check_every_astrid2_ce5.lma.txt` | 出力 |  |
| `v2_3/L2_check_every_kekls.out.txt` | `15_lma_touschek/lma/lma_check_every_kekls.out.txt` | 出力 |  |
| `v2_3/L2_check_every_kekls_ce0.lma.txt` | `15_lma_touschek/lma/lma_check_every_kekls_ce0.lma.txt` | 出力 |  |
| `v2_3/L2_check_every_kekls_ce1.lma.txt` | `15_lma_touschek/lma/lma_check_every_kekls_ce1.lma.txt` | 出力 |  |
| `v2_3/L2_check_every_kekls_ce5.lma.txt` | `15_lma_touschek/lma/lma_check_every_kekls_ce5.lma.txt` | 出力 |  |
| `v2_3/L2_common.py` | `15_lma_touschek/lma/lma_common.py` | スクリプト |  |
| `v2_3/L2_elegant_isr.py` | `15_lma_touschek/touschek_elegant/elegant_isr.py` | スクリプト |  |
| `v2_3/L2_elegant_isr_astrid2.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_astrid2.lma.txt` | 出力 |  |
| `v2_3/L2_elegant_isr_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_elegant_isr_kekls.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_kekls.lma.txt` | 出力 |  |
| `v2_3/L2_elegant_isr_kekls.out.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_kekls.out.txt` | 出力 |  |
| `v2_3/L2_elegant_isr_kekls_n1500.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_kekls_n1500.lma.txt` | 出力 |  |
| `v2_3/L2_elegant_isr_kekls_n1500.out.txt` | `15_lma_touschek/touschek_elegant/elegant_isr_kekls_n1500.out.txt` | 出力 |  |
| `v2_3/L2_errors.py` | `15_lma_touschek/lma/lma_errors.py` | スクリプト |  |
| `v2_3/L2_errors_astrid2.out.txt` | `15_lma_touschek/lma/lma_errors_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_errors_mpw.lma.txt` | `15_lma_touschek/lma/lma_errors_mpw.lma.txt` | 出力 |  |
| `v2_3/L2_errors_mpw.out.txt` | `15_lma_touschek/lma/lma_errors_mpw.out.txt` | 出力 |  |
| `v2_3/L2_errors_seed1.lma.txt` | `15_lma_touschek/lma/lma_errors_seed1.lma.txt` | 出力 |  |
| `v2_3/L2_errors_seed2.lma.txt` | `15_lma_touschek/lma/lma_errors_seed2.lma.txt` | 出力 |  |
| `v2_3/L2_errors_seed3.lma.txt` | `15_lma_touschek/lma/lma_errors_seed3.lma.txt` | 出力 |  |
| `v2_3/L2_fluc.py` | `15_lma_touschek/lma/lma_fluc.py` | スクリプト |  |
| `v2_3/L2_fluc_astrid2.json` | `15_lma_touschek/lma/lma_fluc_astrid2.json` | 出力 |  |
| `v2_3/L2_fluc_astrid2.out.txt` | `15_lma_touschek/lma/lma_fluc_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_fluc_astrid2_long.out.txt` | `15_lma_touschek/lma/lma_fluc_astrid2_long.out.txt` | 出力 |  |
| `v2_3/L2_fluc_astrid2_tol5e-05_n13120.json` | `15_lma_touschek/lma/lma_fluc_astrid2_tol5e-05_n13120.json` | 出力 |  |
| `v2_3/L2_fluc_kekls.json` | `15_lma_touschek/lma/lma_fluc_kekls.json` | 出力 |  |
| `v2_3/L2_fluc_kekls.out.txt` | `15_lma_touschek/lma/lma_fluc_kekls.out.txt` | 出力 |  |
| `v2_3/L2_fluc_kekls_n1500.out.txt` | `15_lma_touschek/lma/lma_fluc_kekls_n1500.out.txt` | 出力 |  |
| `v2_3/L2_fluc_kekls_tol0.0001_n1500.json` | `15_lma_touschek/lma/lma_fluc_kekls_tol0.0001_n1500.json` | 出力 |  |
| `v2_3/L2_limit_types.py` | `15_lma_touschek/lma/lma_limit_types.py` | スクリプト |  |
| `v2_3/L2_limit_types_astrid2.lma.txt` | `15_lma_touschek/lma/lma_limit_types_astrid2.lma.txt` | 出力 |  |
| `v2_3/L2_limit_types_astrid2.out.txt` | `15_lma_touschek/lma/lma_limit_types_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_limit_types_astrid2.table.txt` | `15_lma_touschek/lma/lma_limit_types_astrid2.table.txt` | 出力 |  |
| `v2_3/L2_limit_types_kekls.lma.txt` | `15_lma_touschek/lma/lma_limit_types_kekls.lma.txt` | 出力 |  |
| `v2_3/L2_limit_types_kekls.out.txt` | `15_lma_touschek/lma/lma_limit_types_kekls.out.txt` | 出力 |  |
| `v2_3/L2_limit_types_kekls.table.txt` | `15_lma_touschek/lma/lma_limit_types_kekls.table.txt` | 出力 |  |
| `v2_3/L2_turns.py` | `15_lma_touschek/lma/lma_turns.py` | スクリプト |  |
| `v2_3/L2_turns_kekls.out.txt` | `15_lma_touschek/lma/lma_turns_kekls.out.txt` | 出力 |  |
| `v2_3/L2_turns_kekls_n1000.lma.txt` | `15_lma_touschek/lma/lma_turns_kekls_n1000.lma.txt` | 出力 |  |
| `v2_3/L2_turns_kekls_n1500.lma.txt` | `15_lma_touschek/lma/lma_turns_kekls_n1500.lma.txt` | 出力 |  |
| `v2_3/L2_turns_kekls_n1500_nobucket.lma.txt` | `15_lma_touschek/lma/lma_turns_kekls_n1500_nobucket.lma.txt` | 出力 |  |
| `v2_3/L2_turns_kekls_n2000.lma.txt` | `15_lma_touschek/lma/lma_turns_kekls_n2000.lma.txt` | 出力 |  |
| `v2_3/L2_turns_kekls_n491.lma.txt` | `15_lma_touschek/lma/lma_turns_kekls_n491.lma.txt` | 出力 |  |
| `v2_3/L2_turns_short.out.txt` | `15_lma_touschek/lma/lma_turns_short.out.txt` | 出力 |  |
| `v2_3/L2_turns_short.py` | `15_lma_touschek/lma/lma_turns_short.py` | スクリプト |  |
| `v2_3/L2_verify.py` | `15_lma_touschek/lma/lma_verify.py` | スクリプト |  |
| `v2_3/L2_verify_astrid2.lma.txt` | `15_lma_touschek/lma/lma_verify_astrid2.lma.txt` | 出力 |  |
| `v2_3/L2_verify_astrid2.out.txt` | `15_lma_touschek/lma/lma_verify_astrid2.out.txt` | 出力 |  |
| `v2_3/L2_verify_kekls.lma.txt` | `15_lma_touschek/lma/lma_verify_kekls.lma.txt` | 出力 |  |
| `v2_3/L2_verify_kekls.out.txt` | `15_lma_touschek/lma/lma_verify_kekls.out.txt` | 出力 |  |
| `v2_3/Q1_README.md` | `14_nonlinear_dynamics/frequency_map/fmap_basic.md` | md |  |
| `v2_3/Q1_fmap_astrid2.out.txt` | — | 出力 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_astrid2.out.txt`) |
| `v2_3/Q1_fmap_astrid2.png` | — | 図 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_astrid2.png`) |
| `v2_3/Q1_fmap_astrid2.table.txt` | — | 出力 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_astrid2.table.txt`) |
| `v2_3/Q1_fmap_astrid2_v241.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_astrid2.out.txt` | 出力 |  |
| `v2_3/Q1_fmap_astrid2_v241.png` | `14_nonlinear_dynamics/frequency_map/fmap_astrid2.png` | 図 |  |
| `v2_3/Q1_fmap_astrid2_v241.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_astrid2.table.txt` | 出力 |  |
| `v2_3/Q1_fmap_kekls.out.txt` | — | 出力 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_kekls.out.txt`) |
| `v2_3/Q1_fmap_kekls.png` | — | 図 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_kekls.png`) |
| `v2_3/Q1_fmap_kekls.table.txt` | — | 出力 | 削除 (現行は `14_nonlinear_dynamics/frequency_map/fmap_kekls.table.txt`) |
| `v2_3/Q1_fmap_kekls_v241.err.txt` | `14_nonlinear_dynamics/frequency_map/fmap_kekls.err.txt` | 出力 |  |
| `v2_3/Q1_fmap_kekls_v241.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_kekls.out.txt` | 出力 |  |
| `v2_3/Q1_fmap_kekls_v241.png` | `14_nonlinear_dynamics/frequency_map/fmap_kekls.png` | 図 |  |
| `v2_3/Q1_fmap_kekls_v241.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_kekls.table.txt` | 出力 |  |
| `v2_3/Q1_fmap_ring.py` | `14_nonlinear_dynamics/frequency_map/fmap_ring.py` | スクリプト |  |
| `v2_3/Q1_gs_fmap_astrid2.gs.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_astrid2.gs.txt` | 出力 |  |
| `v2_3/Q1_gs_fmap_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_astrid2.out.txt` | 出力 |  |
| `v2_3/Q1_gs_fmap_astrid2.pysad.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_astrid2.pysad.txt` | 出力 |  |
| `v2_3/Q1_gs_fmap_astrid2.sad` | `14_nonlinear_dynamics/frequency_map/gs_fmap_astrid2.sad` | スクリプト |  |
| `v2_3/Q1_gs_fmap_compare.py` | `14_nonlinear_dynamics/frequency_map/gs_fmap_compare_astrid2.py` | スクリプト |  |
| `v2_3/Q1_gs_loss_astrid2.gs.txt` | `15_lma_touschek/lma/gs_loss_astrid2.gs.txt` | 出力 |  |
| `v2_3/Q1_gs_loss_astrid2.out.txt` | `15_lma_touschek/lma/gs_loss_astrid2.out.txt` | 出力 |  |
| `v2_3/Q1_gs_loss_astrid2.pysad.txt` | `15_lma_touschek/lma/gs_loss_astrid2.pysad.txt` | 出力 |  |
| `v2_3/Q1_gs_loss_astrid2.sad` | `15_lma_touschek/lma/gs_loss_astrid2.sad` | スクリプト |  |
| `v2_3/Q1_gs_loss_compare.py` | `15_lma_touschek/lma/gs_loss_compare.py` | スクリプト |  |
| `v2_3/Q1_lma_astrid2.out.txt` | `15_lma_touschek/lma/lma_astrid2.out.txt` | 出力 |  |
| `v2_3/Q1_lma_astrid2.png` | `15_lma_touschek/lma/lma_astrid2.png` | 図 |  |
| `v2_3/Q1_lma_astrid2.table.txt` | `15_lma_touschek/lma/lma_astrid2.table.txt` | 出力 |  |
| `v2_3/Q1_lma_kekls.out.txt` | `15_lma_touschek/lma/lma_kekls.out.txt` | 出力 |  |
| `v2_3/Q1_lma_kekls.png` | `15_lma_touschek/lma/lma_kekls.png` | 図 |  |
| `v2_3/Q1_lma_kekls.table.txt` | `15_lma_touschek/lma/lma_kekls.table.txt` | 出力 |  |
| `v2_3/Q1_lma_ring.py` | `15_lma_touschek/lma/lma_ring.py` | スクリプト |  |
| `v2_3/Q1_script_commands.pysad.txt` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance.pysad.txt` | 出力 |  |
| `v2_3/Q1_script_commands.sad` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance.sad` | スクリプト |  |
| `v2_3/Q1_script_fmap.png` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance_fmap.png` | 図 |  |
| `v2_3/Q1_script_fmap.table.txt` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance_fmap.table.txt` | 出力 |  |
| `v2_3/Q1_script_lma.png` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance_lma.png` | 図 |  |
| `v2_3/Q1_script_lma.table.txt` | `14_nonlinear_dynamics/frequency_map/fn_FrequencyMap_MomentumAcceptance_lma.table.txt` | 出力 |  |
| `v2_3/Q2_README.md` | `14_nonlinear_dynamics/adts_rdt/rdt_adts_basic.md` | md |  |
| `v2_3/Q2_adts_gs.py` | `14_nonlinear_dynamics/adts_rdt/adts_gs.py` | スクリプト |  |
| `v2_3/Q2_adts_gs_astrid.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_astrid.txt` | 出力 |  |
| `v2_3/Q2_adts_gs_fodo.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_fodo.txt` | 出力 |  |
| `v2_3/Q2_adts_gs_kekls.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls.txt` | 出力 |  |
| `v2_3/Q2_adts_gs_kekls_oct0.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls_oct0.txt` | 出力 |  |
| `v2_3/Q2_adts_tracking.py` | `14_nonlinear_dynamics/adts_rdt/adts_tracking.py` | スクリプト |  |
| `v2_3/Q2_adts_tracking_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_tracking_astrid.pysad.txt` | 出力 |  |
| `v2_3/Q2_adts_tracking_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_tracking_fodo.pysad.txt` | 出力 |  |
| `v2_3/Q2_adts_tracking_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_tracking_kekls.pysad.txt` | 出力 |  |
| `v2_3/Q2_resdrive_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/fn_ResonanceDrivingTerms_fodo.pysad.txt` | 出力 |  |
| `v2_3/Q2_resdrive_fodo.sad` | `14_nonlinear_dynamics/adts_rdt/fn_ResonanceDrivingTerms_fodo.sad` | スクリプト |  |
| `v2_3/Q2_rings.py` | `14_nonlinear_dynamics/adts_rdt/rdt_rings.py` | スクリプト |  |
| `v2_3/Q2_rings.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_rings.pysad.txt` | 出力 |  |
| `v2_3/Q2_rings_elements.png` | `14_nonlinear_dynamics/adts_rdt/rdt_rings_elements.png` | 図 |  |
| `v2_3/Q2_thick_thin.py` | `14_nonlinear_dynamics/adts_rdt/rdt_thick_thin.py` | スクリプト |  |
| `v2_3/Q2_thick_thin.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_thick_thin.pysad.txt` | 出力 |  |
| `v2_3/R1_README.md` | `14_nonlinear_dynamics/adts_rdt/adts_gs_norfsw.md` | md |  |
| `v2_3/R1_adts_gs_astrid_norfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_astrid_norfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_astrid_rfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_astrid_rfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_astrid_rfsw_512.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_astrid_rfsw_512.txt` | 出力 |  |
| `v2_3/R1_adts_gs_fodo_norfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_fodo_norfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_kekls_norfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls_norfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_kekls_oct0_norfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls_oct0_norfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_kekls_oct0_rfsw.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls_oct0_rfsw.txt` | 出力 |  |
| `v2_3/R1_adts_gs_kekls_oct0_rfsw_512.txt` | `14_nonlinear_dynamics/adts_rdt/adts_gs_kekls_oct0_rfsw_512.txt` | 出力 |  |
| `v2_3/R2_README.md` | `14_nonlinear_dynamics/adts_rdt/adts_direct.md` | md |  |
| `v2_3/R2_ablate_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_ablate_astrid.pysad.txt` | 出力 |  |
| `v2_3/R2_ablate_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_ablate_fodo.pysad.txt` | 出力 |  |
| `v2_3/R2_ablate_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_ablate_kekls.pysad.txt` | 出力 |  |
| `v2_3/R2_adts_direct.py` | `14_nonlinear_dynamics/adts_rdt/adts_direct.py` | スクリプト |  |
| `v2_3/R2_adts_direct_astrid.gs.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid.gs.json` | 出力 |  |
| `v2_3/R2_adts_direct_astrid.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid.gs.txt` | 出力 |  |
| `v2_3/R2_adts_direct_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid.pysad.json` | 出力 |  |
| `v2_3/R2_adts_direct_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid.pysad.txt` | 出力 |  |
| `v2_3/R2_adts_direct_astrid_mpw.gs.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid_mpw.gs.json` | 出力 |  |
| `v2_3/R2_adts_direct_astrid_mpw.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid_mpw.gs.txt` | 出力 |  |
| `v2_3/R2_adts_direct_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R2_adts_direct_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R2_adts_direct_fodo.gs.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fodo.gs.json` | 出力 |  |
| `v2_3/R2_adts_direct_fodo.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fodo.gs.txt` | 出力 |  |
| `v2_3/R2_adts_direct_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fodo.pysad.json` | 出力 |  |
| `v2_3/R2_adts_direct_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fodo.pysad.txt` | 出力 |  |
| `v2_3/R2_adts_direct_kekls.gs.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_kekls.gs.json` | 出力 |  |
| `v2_3/R2_adts_direct_kekls.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_kekls.gs.txt` | 出力 |  |
| `v2_3/R2_adts_direct_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_direct_kekls.pysad.json` | 出力 |  |
| `v2_3/R2_adts_direct_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_direct_kekls.pysad.txt` | 出力 |  |
| `v2_3/R2_coupling.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_coupling.pysad.txt` | 出力 |  |
| `v2_3/R2_delta_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_delta_astrid.pysad.txt` | 出力 |  |
| `v2_3/R2_delta_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_delta_fodo.pysad.txt` | 出力 |  |
| `v2_3/R2_delta_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_delta_kekls.pysad.txt` | 出力 |  |
| `v2_3/R2_fit_astrid.png` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fit_astrid.png` | 図 |  |
| `v2_3/R2_fit_astrid_mpw.png` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fit_astrid_mpw.png` | 図 |  |
| `v2_3/R2_fit_fodo.png` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fit_fodo.png` | 図 |  |
| `v2_3/R2_fit_kekls.png` | `14_nonlinear_dynamics/adts_rdt/adts_direct_fit_kekls.png` | 図 |  |
| `v2_3/R2_octupole_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_octupole_kekls.pysad.txt` | 出力 |  |
| `v2_3/R2_resonance.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_resonance.pysad.txt` | 出力 |  |
| `v2_3/R2_unverified.py` | `14_nonlinear_dynamics/adts_rdt/adts_checks.py` | スクリプト |  |
| `v2_3/R3_README.md` | `14_nonlinear_dynamics/adts_rdt/adts_bend_cubic.md` | md |  |
| `v2_3/R3_adts_bend.png` | `14_nonlinear_dynamics/adts_rdt/adts_bend.png` | 図 |  |
| `v2_3/R3_adts_bend.py` | `14_nonlinear_dynamics/adts_rdt/adts_bend.py` | スクリプト |  |
| `v2_3/R3_adts_gs_astrid_lam0.5.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_gs_astrid_lam0.5.gs.txt` | 出力 |  |
| `v2_3/R3_adts_gs_fodo_lam2.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_gs_fodo_lam2.gs.txt` | 出力 |  |
| `v2_3/R3_adts_gs_kekls_lam2.gs.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_gs_kekls_lam2.gs.txt` | 出力 |  |
| `v2_3/R3_adts_scan.png` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan.png` | 図 |  |
| `v2_3/R3_conv_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_conv_astrid.pysad.txt` | 出力 |  |
| `v2_3/R3_conv_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_conv_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R3_conv_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_conv_fodo.pysad.txt` | 出力 |  |
| `v2_3/R3_conv_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_conv_kekls.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_astrid.pysad.json` | 出力 |  |
| `v2_3/R3_scan_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_astrid.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R3_scan_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo.pysad.json` | 出力 |  |
| `v2_3/R3_scan_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_fodo_asym.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo_asym.pysad.json` | 出力 |  |
| `v2_3/R3_scan_fodo_asym.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo_asym.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_fodo_rect.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo_rect.pysad.json` | 出力 |  |
| `v2_3/R3_scan_fodo_rect.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_fodo_rect.pysad.txt` | 出力 |  |
| `v2_3/R3_scan_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_kekls.pysad.json` | 出力 |  |
| `v2_3/R3_scan_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/adts_bend_scan_kekls.pysad.txt` | 出力 |  |
| `v2_3/R4_README.md` | `14_nonlinear_dynamics/adts_rdt/rdt_tracking.md` | md |  |
| `v2_3/R4_chrom_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_astrid.pysad.json` | 出力 |  |
| `v2_3/R4_chrom_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_astrid.pysad.txt` | 出力 |  |
| `v2_3/R4_chrom_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R4_chrom_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R4_chrom_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_fodo.pysad.json` | 出力 |  |
| `v2_3/R4_chrom_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_fodo.pysad.txt` | 出力 |  |
| `v2_3/R4_chrom_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_kekls.pysad.json` | 出力 |  |
| `v2_3/R4_chrom_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_chrom_kekls.pysad.txt` | 出力 |  |
| `v2_3/R4_gs_map.py` | `14_nonlinear_dynamics/adts_rdt/rdt_gs_map.py` | スクリプト |  |
| `v2_3/R4_gs_map_astrid.gs.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_gs_map_astrid.gs.txt` | 出力 |  |
| `v2_3/R4_gs_map_fodo.gs.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_gs_map_fodo.gs.txt` | 出力 |  |
| `v2_3/R4_gs_map_kekls.gs.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_gs_map_kekls.gs.txt` | 出力 |  |
| `v2_3/R4_rdt_table.png` | `14_nonlinear_dynamics/adts_rdt/rdt_table.png` | 図 |  |
| `v2_3/R4_rdt_tracking.py` | `14_nonlinear_dynamics/adts_rdt/rdt_tracking.py` | スクリプト |  |
| `v2_3/R4_res_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_res_astrid.pysad.json` | 出力 |  |
| `v2_3/R4_res_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_res_astrid.pysad.txt` | 出力 |  |
| `v2_3/R4_res_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_res_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R4_res_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_res_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R4_res_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_res_fodo.pysad.json` | 出力 |  |
| `v2_3/R4_res_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_res_fodo.pysad.txt` | 出力 |  |
| `v2_3/R4_res_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_res_kekls.pysad.json` | 出力 |  |
| `v2_3/R4_res_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_res_kekls.pysad.txt` | 出力 |  |
| `v2_3/R4_table_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_astrid.pysad.json` | 出力 |  |
| `v2_3/R4_table_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_astrid.pysad.txt` | 出力 |  |
| `v2_3/R4_table_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R4_table_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R4_table_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo.pysad.json` | 出力 |  |
| `v2_3/R4_table_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo.pysad.txt` | 出力 |  |
| `v2_3/R4_table_fodo_asym.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo_asym.pysad.json` | 出力 |  |
| `v2_3/R4_table_fodo_asym.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo_asym.pysad.txt` | 出力 |  |
| `v2_3/R4_table_fodo_rect.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo_rect.pysad.json` | 出力 |  |
| `v2_3/R4_table_fodo_rect.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_fodo_rect.pysad.txt` | 出力 |  |
| `v2_3/R4_table_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt_table_kekls.pysad.json` | 出力 |  |
| `v2_3/R4_table_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_table_kekls.pysad.txt` | 出力 |  |
| `v2_3/R4_thin_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt_thin_fodo.pysad.txt` | 出力 |  |
| `v2_3/R5_README.md` | `14_nonlinear_dynamics/adts_rdt/rdt4_normal_form.md` | md |  |
| `v2_3/R5_adts_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_astrid.pysad.json` | 出力 |  |
| `v2_3/R5_adts_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_astrid.pysad.txt` | 出力 |  |
| `v2_3/R5_adts_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R5_adts_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R5_adts_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_fodo.pysad.json` | 出力 |  |
| `v2_3/R5_adts_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_fodo.pysad.txt` | 出力 |  |
| `v2_3/R5_adts_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_kekls.pysad.json` | 出力 |  |
| `v2_3/R5_adts_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_adts_kekls.pysad.txt` | 出力 |  |
| `v2_3/R5_bend_generators.out.txt` | `14_nonlinear_dynamics/adts_rdt/bend_generators.out.txt` | 出力 |  |
| `v2_3/R5_bend_generators.py` | `14_nonlinear_dynamics/adts_rdt/bend_generators.py` | スクリプト |  |
| `v2_3/R5_fmap_rdt2.py` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4.py` | スクリプト |  |
| `v2_3/R5_fmap_rdt2_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_astrid2.out.txt` | 出力 |  |
| `v2_3/R5_fmap_rdt2_astrid2.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_astrid2.png` | 図 |  |
| `v2_3/R5_fmap_rdt2_fodo.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_fodo.out.txt` | 出力 |  |
| `v2_3/R5_fmap_rdt2_fodo.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_fodo.png` | 図 |  |
| `v2_3/R5_fmap_rdt2_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_kekls.out.txt` | 出力 |  |
| `v2_3/R5_fmap_rdt2_kekls.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt4_kekls.png` | 図 |  |
| `v2_3/R5_gs_map.py` | `14_nonlinear_dynamics/adts_rdt/rdt4_gs_map.py` | スクリプト |  |
| `v2_3/R5_gs_map_astrid.gs.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_gs_map_astrid.gs.txt` | 出力 |  |
| `v2_3/R5_gs_map_fodo.gs.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_gs_map_fodo.gs.txt` | 出力 |  |
| `v2_3/R5_nform.out.txt` | `14_nonlinear_dynamics/adts_rdt/normal_form_numeric.out.txt` | 出力 |  |
| `v2_3/R5_nform.py` | `14_nonlinear_dynamics/adts_rdt/normal_form_numeric.py` | スクリプト |  |
| `v2_3/R5_rdt4_table.png` | `14_nonlinear_dynamics/adts_rdt/rdt4_table.png` | 図 |  |
| `v2_3/R5_rdt4_tracking.py` | `14_nonlinear_dynamics/adts_rdt/rdt4_tracking.py` | スクリプト |  |
| `v2_3/R5_table_astrid.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_astrid.pysad.json` | 出力 |  |
| `v2_3/R5_table_astrid.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_astrid.pysad.txt` | 出力 |  |
| `v2_3/R5_table_astrid_mpw.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_astrid_mpw.pysad.json` | 出力 |  |
| `v2_3/R5_table_astrid_mpw.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_astrid_mpw.pysad.txt` | 出力 |  |
| `v2_3/R5_table_fodo.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo.pysad.json` | 出力 |  |
| `v2_3/R5_table_fodo.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo.pysad.txt` | 出力 |  |
| `v2_3/R5_table_fodo_asym.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo_asym.pysad.json` | 出力 |  |
| `v2_3/R5_table_fodo_asym.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo_asym.pysad.txt` | 出力 |  |
| `v2_3/R5_table_fodo_rect.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo_rect.pysad.json` | 出力 |  |
| `v2_3/R5_table_fodo_rect.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_fodo_rect.pysad.txt` | 出力 |  |
| `v2_3/R5_table_kekls.pysad.json` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_kekls.pysad.json` | 出力 |  |
| `v2_3/R5_table_kekls.pysad.txt` | `14_nonlinear_dynamics/adts_rdt/rdt4_table_kekls.pysad.txt` | 出力 |  |
| `v2_3/R6_README.md` | `14_nonlinear_dynamics/adts_rdt/adts_exact_terms.md` | md |  |
| `v2_3/R6_adts.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_adts.out.txt` | 出力 |  |
| `v2_3/R6_bend.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_bend.out.txt` | 出力 |  |
| `v2_3/R6_delta_astrid.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_astrid.out.txt` | 出力 |  |
| `v2_3/R6_delta_astrid.track.json` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_astrid.track.json` | 出力 |  |
| `v2_3/R6_delta_fodo.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_fodo.out.txt` | 出力 |  |
| `v2_3/R6_delta_fodo.track.json` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_fodo.track.json` | 出力 |  |
| `v2_3/R6_delta_kekls.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_kekls.out.txt` | 出力 |  |
| `v2_3/R6_delta_kekls.track.json` | `14_nonlinear_dynamics/adts_rdt/adts_exact_delta_kekls.track.json` | 出力 |  |
| `v2_3/R6_mpw.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_mpw.out.txt` | 出力 |  |
| `v2_3/R6_nf.out.txt` | `14_nonlinear_dynamics/adts_rdt/adts_exact_nf.out.txt` | 出力 |  |
| `v2_3/R6_resolve.py` | `14_nonlinear_dynamics/adts_rdt/adts_exact.py` | スクリプト |  |
| `v2_3/README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 14_nonlinear_dynamics/adts_rdt, 14_nonlinear_dynamics/frequency_map, 14_nonlinear_dynamics/frequency_map/fmap_check_data, 14_nonlinear_dynamics/primer, 15_lma_touschek/lma, 15_lma_touschek/touschek_elegant |
| `v2_3/S1_README.md` | `14_nonlinear_dynamics/frequency_map/fmap_validation.md` | md |  |
| `v2_3/S1_fmap_delta_astrid2_xdp.png` | `14_nonlinear_dynamics/frequency_map/fmap_delta_astrid2_xdp.png` | 図 |  |
| `v2_3/S1_fmap_delta_astrid2_xdp.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_delta_astrid2_xdp.table.txt` | 出力 |  |
| `v2_3/S1_fmap_delta_astrid2_ydp.png` | `14_nonlinear_dynamics/frequency_map/fmap_delta_astrid2_ydp.png` | 図 |  |
| `v2_3/S1_fmap_delta_astrid2_ydp.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_delta_astrid2_ydp.table.txt` | 出力 |  |
| `v2_3/S1_fmap_delta_kekls_xdp.png` | `14_nonlinear_dynamics/frequency_map/fmap_delta_kekls_xdp.png` | 図 |  |
| `v2_3/S1_fmap_delta_kekls_xdp.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_delta_kekls_xdp.table.txt` | 出力 |  |
| `v2_3/S1_fmap_errors_astrid2_L1.png` | `14_nonlinear_dynamics/frequency_map/fmap_errors_astrid2_L1.png` | 図 |  |
| `v2_3/S1_fmap_errors_astrid2_L1.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_errors_astrid2_L1.table.txt` | 出力 |  |
| `v2_3/S1_fmap_errors_astrid2_L3.png` | `14_nonlinear_dynamics/frequency_map/fmap_errors_astrid2_L3.png` | 図 |  |
| `v2_3/S1_fmap_errors_astrid2_L3.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_errors_astrid2_L3.table.txt` | 出力 |  |
| `v2_3/S1_fmap_grid_astrid2_81x61.png` | `14_nonlinear_dynamics/frequency_map/fmap_grid_astrid2_81x61.png` | 図 |  |
| `v2_3/S1_fmap_grid_astrid2_81x61.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_grid_astrid2_81x61.table.txt` | 出力 |  |
| `v2_3/S1_fmap_grid_summary.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_grid_summary.out.txt` | 出力 |  |
| `v2_3/S1_fmap_grid_summary.py` | `14_nonlinear_dynamics/frequency_map/fmap_grid_summary.py` | スクリプト |  |
| `v2_3/S1_fmap_k64_astrid2.png` | `14_nonlinear_dynamics/frequency_map/fmap_k64_astrid2.png` | 図 |  |
| `v2_3/S1_fmap_k64_astrid2.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_k64_astrid2.table.txt` | 出力 |  |
| `v2_3/S1_fmap_k64_kekls.png` | `14_nonlinear_dynamics/frequency_map/fmap_k64_kekls.png` | 図 |  |
| `v2_3/S1_fmap_k64_kekls.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_k64_kekls.table.txt` | 出力 |  |
| `v2_3/S1_fmap_lambda_astrid2_0.5.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_lambda_astrid2_0.5.table.txt` | 出力 |  |
| `v2_3/S1_fmap_lambda_astrid2_0.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_lambda_astrid2_0.table.txt` | 出力 |  |
| `v2_3/S1_fmap_lambda_astrid2_1.5.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_lambda_astrid2_1.5.table.txt` | 出力 |  |
| `v2_3/S1_fmap_lambda_astrid2_1.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_lambda_astrid2_1.table.txt` | 出力 |  |
| `v2_3/S1_fmap_mpw.png` | `14_nonlinear_dynamics/frequency_map/fmap_mpw.png` | 図 |  |
| `v2_3/S1_fmap_mpw.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_mpw.table.txt` | 出力 |  |
| `v2_3/S1_fmap_naff_astrid2.base.txt` | `14_nonlinear_dynamics/frequency_map/fmap_naff_astrid2.base.txt` | データ |  |
| `v2_3/S1_fmap_rdt.py` | `14_nonlinear_dynamics/frequency_map/fmap_rdt.py` | スクリプト |  |
| `v2_3/S1_fmap_rdt_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_astrid2.out.txt` | 出力 |  |
| `v2_3/S1_fmap_rdt_astrid2.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_astrid2.png` | 図 |  |
| `v2_3/S1_fmap_rdt_fodo.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_fodo.out.txt` | 出力 |  |
| `v2_3/S1_fmap_rdt_fodo.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_fodo.png` | 図 |  |
| `v2_3/S1_fmap_rdt_fodo.table.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_fodo.table.txt` | 出力 |  |
| `v2_3/S1_fmap_rdt_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_kekls.out.txt` | 出力 |  |
| `v2_3/S1_fmap_rdt_kekls.png` | `14_nonlinear_dynamics/frequency_map/fmap_rdt_kekls.png` | 図 |  |
| `v2_3/S1_fmap_variants.py` | `14_nonlinear_dynamics/frequency_map/fmap_variants.py` | スクリプト |  |
| `v2_3/S1_fmap_variants_delta_astrid2_xdp.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_delta_astrid2_xdp.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_delta_astrid2_ydp.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_delta_astrid2_ydp.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_delta_kekls_xdp.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_delta_kekls_xdp.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_errors_astrid2_1.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_errors_astrid2_1.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_errors_astrid2_3.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_errors_astrid2_3.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_grid_astrid2_81x61.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_grid_astrid2_81x61.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_k64_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_k64_astrid2.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_k64_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_k64_kekls.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_lambda_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_lambda_astrid2.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_mpw.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_mpw.out.txt` | 出力 |  |
| `v2_3/S1_fmap_variants_naff_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_variants_naff_astrid2.out.txt` | 出力 |  |
| `v2_3/S1_gs_fmap_compare.py` | `14_nonlinear_dynamics/frequency_map/gs_fmap_compare_kekls.py` | スクリプト |  |
| `v2_3/S1_gs_fmap_kekls.gs.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_kekls.gs.txt` | 出力 |  |
| `v2_3/S1_gs_fmap_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_kekls.out.txt` | 出力 |  |
| `v2_3/S1_gs_fmap_kekls.pysad.txt` | `14_nonlinear_dynamics/frequency_map/gs_fmap_kekls.pysad.txt` | 出力 |  |
| `v2_3/S1_gs_fmap_kekls.sad` | `14_nonlinear_dynamics/frequency_map/gs_fmap_kekls.sad` | スクリプト |  |
| `v2_3/S2_README.md` | `14_nonlinear_dynamics/frequency_map/fmap_reliability.md` | md |  |
| `v2_3/S2_coupling.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_coupling.py` | スクリプト |  |
| `v2_3/S2_data/S2_coupling.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/coupling.jsonl` | データ |  |
| `v2_3/S2_data/S2_errors_astrid2.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/errors_astrid2.jsonl` | データ |  |
| `v2_3/S2_data/S2_errors_kekls.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/errors_kekls.jsonl` | データ |  |
| `v2_3/S2_data/S2_k1opts.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/k1opts.jsonl` | データ |  |
| `v2_3/S2_data/S2_mpw.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/mpw.jsonl` | データ |  |
| `v2_3/S2_data/S2_nturn_astrid2.npz` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/nturn_astrid2.npz` | データ |  |
| `v2_3/S2_data/S2_nturn_astrid2_err_s1_k1.npz` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/nturn_astrid2_err_s1_k1.npz` | データ |  |
| `v2_3/S2_data/S2_nturn_kekls.npz` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/nturn_kekls.npz` | データ |  |
| `v2_3/S2_data/S2_rfsw.jsonl` | `14_nonlinear_dynamics/frequency_map/fmap_check_data/rfsw.jsonl` | データ |  |
| `v2_3/S2_errors_scan.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_errors_scan.py` | スクリプト |  |
| `v2_3/S2_errors_summary.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_errors_summary.py` | スクリプト |  |
| `v2_3/S2_errors_summary_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_errors_summary_astrid2.out.txt` | 出力 |  |
| `v2_3/S2_errors_summary_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_errors_summary_kekls.out.txt` | 出力 |  |
| `v2_3/S2_grid.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_grid.out.txt` | 出力 |  |
| `v2_3/S2_grid.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_grid.py` | スクリプト |  |
| `v2_3/S2_gs_errors.gs.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.gs.txt` | 出力 |  |
| `v2_3/S2_gs_errors.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.out.txt` | 出力 |  |
| `v2_3/S2_gs_errors.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.py` | スクリプト |  |
| `v2_3/S2_gs_errors.pysad.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.pysad.txt` | 出力 |  |
| `v2_3/S2_gs_errors.sad` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_errors.sad` | スクリプト |  |
| `v2_3/S2_gs_rfsw.gs.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.gs.txt` | 出力 |  |
| `v2_3/S2_gs_rfsw.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.out.txt` | 出力 |  |
| `v2_3/S2_gs_rfsw.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.py` | スクリプト |  |
| `v2_3/S2_gs_rfsw.pysad.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.pysad.txt` | 出力 |  |
| `v2_3/S2_gs_rfsw.sad` | `14_nonlinear_dynamics/frequency_map/fmap_check_gs_rfsw.sad` | スクリプト |  |
| `v2_3/S2_k1opts.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_k1opts.out.txt` | 出力 |  |
| `v2_3/S2_k1opts.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_k1opts.py` | スクリプト |  |
| `v2_3/S2_kekls_orbit.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_kekls_orbit.out.txt` | 出力 |  |
| `v2_3/S2_kekls_orbit.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_kekls_orbit.py` | スクリプト |  |
| `v2_3/S2_mpw.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_mpw.out.txt` | 出力 |  |
| `v2_3/S2_mpw.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_mpw.py` | スクリプト |  |
| `v2_3/S2_nturn.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_nturn.py` | スクリプト |  |
| `v2_3/S2_nturn_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_nturn_astrid2.out.txt` | 出力 |  |
| `v2_3/S2_nturn_astrid2_err_s1_k1.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_nturn_astrid2_err_s1_k1.out.txt` | 出力 |  |
| `v2_3/S2_nturn_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_check_nturn_kekls.out.txt` | 出力 |  |
| `v2_3/S2_rfsw.py` | `14_nonlinear_dynamics/frequency_map/fmap_check_rfsw.py` | スクリプト |  |
| `v2_3/S3_README.md` | `14_nonlinear_dynamics/frequency_map/fmap_resonance_lines.md` | md |  |
| `v2_3/S3_fmap_predict.py` | `14_nonlinear_dynamics/frequency_map/fmap_predict.py` | スクリプト |  |
| `v2_3/S3_fmap_predict_astrid2.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_predict_astrid2.out.txt` | 出力 |  |
| `v2_3/S3_fmap_predict_astrid2.png` | `14_nonlinear_dynamics/frequency_map/fmap_predict_astrid2.png` | 図 |  |
| `v2_3/S3_fmap_predict_fodo.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_predict_fodo.out.txt` | 出力 |  |
| `v2_3/S3_fmap_predict_fodo.png` | `14_nonlinear_dynamics/frequency_map/fmap_predict_fodo.png` | 図 |  |
| `v2_3/S3_fmap_predict_kekls.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_predict_kekls.out.txt` | 出力 |  |
| `v2_3/S3_fmap_predict_kekls.png` | `14_nonlinear_dynamics/frequency_map/fmap_predict_kekls.png` | 図 |  |
| `v2_3/S3_fmap_predict_mpw.out.txt` | `14_nonlinear_dynamics/frequency_map/fmap_predict_mpw.out.txt` | 出力 |  |
| `v2_3/S3_fmap_predict_mpw.png` | `14_nonlinear_dynamics/frequency_map/fmap_predict_mpw.png` | 図 |  |
| `v2_3/T1_README.md` | `15_lma_touschek/touschek_elegant/elegant_comparison.md` | md |  |
| `v2_3/T1_astrid2_L2isr.lte` | `15_lma_touschek/touschek_elegant/astrid2_L2isr.lte` | データ |  |
| `v2_3/T1_astrid2_L2isr.map.json` | `15_lma_touschek/touschek_elegant/astrid2_L2isr.map.json` | データ |  |
| `v2_3/T1_astrid2_all.lte` | `15_lma_touschek/touschek_elegant/astrid2_all.lte` | データ |  |
| `v2_3/T1_astrid2_all.twiss.ele` | `15_lma_touschek/touschek_elegant/astrid2_all.twiss.ele` | データ |  |
| `v2_3/T1_astrid2_ma.ele` | `15_lma_touschek/touschek_elegant/astrid2_ma.ele` | データ |  |
| `v2_3/T1_astrid2_ma.lte` | `15_lma_touschek/touschek_elegant/astrid2_ma.lte` | データ |  |
| `v2_3/T1_astrid2_ma_exactk1.lte` | `15_lma_touschek/touschek_elegant/astrid2_ma_exactk1.lte` | データ |  |
| `v2_3/T1_bend_k1_elegant.out.txt` | `15_lma_touschek/touschek_elegant/bend_k1_elegant.out.txt` | 出力 |  |
| `v2_3/T1_bend_k1_elegant.py` | `15_lma_touschek/touschek_elegant/bend_k1_elegant.py` | スクリプト |  |
| `v2_3/T1_elegant_astrid2.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_astrid2.lma.txt` | 出力 |  |
| `v2_3/T1_elegant_astrid2.touschek.txt` | `15_lma_touschek/touschek_elegant/elegant_astrid2.touschek.txt` | 出力 |  |
| `v2_3/T1_elegant_astrid2_exactk1.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_astrid2_exactk1.lma.txt` | 出力 |  |
| `v2_3/T1_elegant_kekls.lma.txt` | `15_lma_touschek/touschek_elegant/elegant_kekls.lma.txt` | 出力 |  |
| `v2_3/T1_elegant_kekls.touschek.txt` | `15_lma_touschek/touschek_elegant/elegant_kekls.touschek.txt` | 出力 |  |
| `v2_3/T1_elegant_lma.py` | `15_lma_touschek/touschek_elegant/elegant_lma.py` | スクリプト |  |
| `v2_3/T1_elegant_touschek.py` | `15_lma_touschek/touschek_elegant/elegant_touschek.py` | スクリプト |  |
| `v2_3/T1_gs_bend_k1.out.txt` | `15_lma_touschek/touschek_elegant/gs_bend_k1.out.txt` | 出力 |  |
| `v2_3/T1_gs_bend_k1.py` | `15_lma_touschek/touschek_elegant/gs_bend_k1.py` | スクリプト |  |
| `v2_3/T1_gs_bend_k1.sad` | `15_lma_touschek/touschek_elegant/gs_bend_k1.sad` | スクリプト |  |
| `v2_3/T1_kekls_L2isr.lte` | `15_lma_touschek/touschek_elegant/kekls_L2isr.lte` | データ |  |
| `v2_3/T1_kekls_L2isr.map.json` | `15_lma_touschek/touschek_elegant/kekls_L2isr.map.json` | データ |  |
| `v2_3/T1_kekls_all.lte` | `15_lma_touschek/touschek_elegant/kekls_all.lte` | データ |  |
| `v2_3/T1_kekls_all.twiss.ele` | `15_lma_touschek/touschek_elegant/kekls_all.twiss.ele` | データ |  |
| `v2_3/T1_kekls_ma.lte` | `15_lma_touschek/touschek_elegant/kekls_ma.lte` | データ |  |
| `v2_3/T1_kekls_ma_p000.ele` | `15_lma_touschek/touschek_elegant/kekls_ma_p000.ele` | データ |  |
| `v2_3/T1_lma_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/lma_compare_astrid2.out.txt` | 出力 |  |
| `v2_3/T1_lma_astrid2.png` | `15_lma_touschek/touschek_elegant/lma_compare_astrid2.png` | 図 |  |
| `v2_3/T1_lma_astrid2.table.txt` | `15_lma_touschek/touschek_elegant/lma_compare_astrid2.table.txt` | 出力 |  |
| `v2_3/T1_lma_compare.py` | `15_lma_touschek/touschek_elegant/lma_compare.py` | スクリプト |  |
| `v2_3/T1_lma_kekls.out.txt` | `15_lma_touschek/touschek_elegant/lma_compare_kekls.out.txt` | 出力 |  |
| `v2_3/T1_lma_kekls.png` | `15_lma_touschek/touschek_elegant/lma_compare_kekls.png` | 図 |  |
| `v2_3/T1_lma_kekls.table.txt` | `15_lma_touschek/touschek_elegant/lma_compare_kekls.table.txt` | 出力 |  |
| `v2_3/T1_pysad_astrid2.ref.json` | `15_lma_touschek/touschek_elegant/pysad_astrid2.ref.json` | 出力 |  |
| `v2_3/T1_pysad_astrid2.ref.out.txt` | `15_lma_touschek/touschek_elegant/pysad_astrid2.ref.out.txt` | 出力 |  |
| `v2_3/T1_pysad_kekls.lma.txt` | `15_lma_touschek/touschek_elegant/pysad_kekls.lma.txt` | 出力 |  |
| `v2_3/T1_pysad_kekls.ref.json` | `15_lma_touschek/touschek_elegant/pysad_kekls.ref.json` | 出力 |  |
| `v2_3/T1_pysad_kekls.ref.out.txt` | `15_lma_touschek/touschek_elegant/pysad_kekls.ref.out.txt` | 出力 |  |
| `v2_3/T1_pysad_ref.py` | `15_lma_touschek/touschek_elegant/pysad_ref.py` | スクリプト |  |
| `v2_3/T1_sad2ele.py` | `15_lma_touschek/touschek_elegant/sad2ele.py` | スクリプト |  |
| `v2_3/T1_summary.out.txt` | `15_lma_touschek/touschek_elegant/elegant_summary.out.txt` | 出力 |  |
| `v2_3/T1_summary.py` | `15_lma_touschek/touschek_elegant/elegant_summary.py` | スクリプト |  |
| `v2_3/T1_touschek_scan.py` | `15_lma_touschek/touschek_elegant/elegant_touschek_scan.py` | スクリプト |  |
| `v2_3/T1_touschek_scan_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/elegant_touschek_scan_astrid2.out.txt` | 出力 |  |
| `v2_3/T2_README.md` | `15_lma_touschek/touschek_elegant/touschek_formula.md` | md |  |
| `v2_3/T2_astrid2_div.ele` | `15_lma_touschek/touschek_elegant/touschek_formula_astrid2_div.ele` | データ |  |
| `v2_3/T2_cause.py` | `15_lma_touschek/touschek_elegant/touschek_formula_cause.py` | スクリプト |  |
| `v2_3/T2_cause_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_cause_astrid2.out.txt` | 出力 |  |
| `v2_3/T2_cause_kekls.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_cause_kekls.out.txt` | 出力 |  |
| `v2_3/T2_check3d.py` | `15_lma_touschek/touschek_elegant/touschek_formula_check3d.py` | スクリプト |  |
| `v2_3/T2_kekls_div.ele` | `15_lma_touschek/touschek_elegant/touschek_formula_kekls_div.ele` | データ |  |
| `v2_3/T2_limits.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_limits.out.txt` | 出力 |  |
| `v2_3/T2_limits.py` | `15_lma_touschek/touschek_elegant/touschek_formula_limits.py` | スクリプト |  |
| `v2_3/T2_lma_lifetime.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_lma_lifetime.out.txt` | 出力 |  |
| `v2_3/T2_lma_lifetime.py` | `15_lma_touschek/touschek_elegant/touschek_formula_lma_lifetime.py` | スクリプト |  |
| `v2_3/T2_mc.py` | `15_lma_touschek/touschek_elegant/touschek_formula_mc.py` | スクリプト |  |
| `v2_3/T2_scan.py` | `15_lma_touschek/touschek_elegant/touschek_formula_scan.py` | スクリプト |  |
| `v2_3/T2_scan_astrid2.json` | `15_lma_touschek/touschek_elegant/touschek_formula_scan_astrid2.json` | 出力 |  |
| `v2_3/T2_scan_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_scan_astrid2.out.txt` | 出力 |  |
| `v2_3/T2_scan_kekls.json` | `15_lma_touschek/touschek_elegant/touschek_formula_scan_kekls.json` | 出力 |  |
| `v2_3/T2_scan_kekls.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_scan_kekls.out.txt` | 出力 |  |
| `v2_3/T2_validate.py` | `15_lma_touschek/touschek_elegant/touschek_formula_validate.py` | スクリプト |  |
| `v2_3/T2_validate_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_validate_astrid2.out.txt` | 出力 |  |
| `v2_3/T2_validate_kekls.out.txt` | `15_lma_touschek/touschek_elegant/touschek_formula_validate_kekls.out.txt` | 出力 |  |
| `v2_3/U2_README.md` | `15_lma_touschek/touschek_elegant/bend_k1_definition.md` | md |  |
| `v2_3/U2_astrid2_chrom.out.txt` | `15_lma_touschek/touschek_elegant/bend_k1_def_astrid2_chrom.out.txt` | 出力 |  |
| `v2_3/U2_astrid2_chrom.py` | `15_lma_touschek/touschek_elegant/bend_k1_def_astrid2_chrom.py` | スクリプト |  |
| `v2_3/U2_lma_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/bend_k1_def_lma_astrid2.out.txt` | 出力 |  |
| `v2_3/U2_lma_astrid2.py` | `15_lma_touschek/touschek_elegant/bend_k1_def_lma_astrid2.py` | スクリプト |  |
| `v2_3/U2_verify.out.txt` | `15_lma_touschek/touschek_elegant/bend_k1_def_verify.out.txt` | 出力 |  |
| `v2_3/U2_verify.py` | `15_lma_touschek/touschek_elegant/bend_k1_def_verify.py` | スクリプト |  |
| `v2_3/primer_adts_fit.png` | `14_nonlinear_dynamics/primer/primer_adts_fit.png` | 図 |  |
| `v2_3/primer_adts_fit.py` | `14_nonlinear_dynamics/primer/primer_adts_fit.py` | スクリプト |  |
| `v2_3/primer_breakdown.png` | `14_nonlinear_dynamics/primer/primer_breakdown.png` | 図 |  |
| `v2_3/primer_breakdown.py` | `14_nonlinear_dynamics/primer/primer_breakdown.py` | スクリプト |  |
| `v2_3/primer_common.py` | `14_nonlinear_dynamics/primer/primer_common.py` | スクリプト |  |
| `v2_3/primer_fmap_astrid2.png` | `14_nonlinear_dynamics/primer/primer_fmap_astrid2.png` | 図 |  |
| `v2_3/primer_fmap_astrid2.table.txt` | `14_nonlinear_dynamics/primer/primer_fmap_astrid2.table.txt` | 出力 |  |
| `v2_3/primer_fmap_lma.py` | `14_nonlinear_dynamics/primer/primer_fmap_lma.py` | スクリプト |  |
| `v2_3/primer_lma_astrid2.png` | `14_nonlinear_dynamics/primer/primer_lma_astrid2.png` | 図 |  |
| `v2_3/primer_lma_astrid2.table.txt` | `14_nonlinear_dynamics/primer/primer_lma_astrid2.table.txt` | 出力 |  |
| `v2_3/primer_lma_scan.png` | `14_nonlinear_dynamics/primer/primer_lma_scan.png` | 図 |  |
| `v2_3/primer_naff.png` | `14_nonlinear_dynamics/primer/primer_naff.png` | 図 |  |
| `v2_3/primer_naff.py` | `14_nonlinear_dynamics/primer/primer_naff.py` | スクリプト |  |
| `v2_3/primer_octupole.py` | `14_nonlinear_dynamics/primer/primer_octupole.py` | スクリプト |  |
| `v2_3/primer_phase_space.png` | `14_nonlinear_dynamics/primer/primer_phase_space.png` | 図 |  |
| `v2_3/primer_phase_space.py` | `14_nonlinear_dynamics/primer/primer_phase_space.py` | スクリプト |  |
| `v2_3/primer_rdt.png` | `14_nonlinear_dynamics/primer/primer_rdt.png` | 図 |  |
| `v2_3/primer_rdt.py` | `14_nonlinear_dynamics/primer/primer_rdt.py` | スクリプト |  |
| `v2_3/primer_sext_even_odd.py` | `14_nonlinear_dynamics/primer/primer_sext_even_odd.py` | スクリプト |  |
| `v2_4/A1_ampnf_vs_tracking.py` | `14_nonlinear_dynamics/normal_form/ampnf_vs_tracking.py` | スクリプト |  |
| `v2_4/A1_astrid.out.txt` | `14_nonlinear_dynamics/normal_form/ampnf_vs_tracking_astrid.out.txt` | 出力 |  |
| `v2_4/A1_fodo.out.txt` | `14_nonlinear_dynamics/normal_form/ampnf_vs_tracking_fodo.out.txt` | 出力 |  |
| `v2_4/A2_k0_bend_check.py` | `14_nonlinear_dynamics/normal_form/k0_bend_check.py` | スクリプト |  |
| `v2_4/A3_mc_bunch.py` | `15_lma_touschek/touschek_elegant/touschek_bunch_mc.py` | スクリプト |  |
| `v2_4/A3_mc_bunch_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/touschek_bunch_mc_astrid2.out.txt` | 出力 |  |
| `v2_4/A3_mc_bunch_kekls.out.txt` | `15_lma_touschek/touschek_elegant/touschek_bunch_mc_kekls.out.txt` | 出力 |  |
| `v2_4/A3_qs_3d.py` | `15_lma_touschek/touschek_elegant/touschek_bunch_qs_3d.py` | スクリプト |  |
| `v2_4/A3_qs_3d_astrid2.out.txt` | `15_lma_touschek/touschek_elegant/touschek_bunch_qs_3d_astrid2.out.txt` | 出力 |  |
| `v2_4/A3_qs_3d_kekls.out.txt` | `15_lma_touschek/touschek_elegant/touschek_bunch_qs_3d_kekls.out.txt` | 出力 |  |
| `v2_4/A3_qs_ring.out.txt` | `15_lma_touschek/touschek_elegant/touschek_bunch_qs_ring.out.txt` | 出力 |  |
| `v2_4/A3_qs_ring.py` | `15_lma_touschek/touschek_elegant/touschek_bunch_qs_ring.py` | スクリプト |  |
| `v2_4/A4_coupl.out.txt` | `15_lma_touschek/lma/lma_support_coupl.out.txt` | 出力 |  |
| `v2_4/A4_k64.out.txt` | `15_lma_touschek/lma/lma_support_k64.out.txt` | 出力 |  |
| `v2_4/A4_kekls.out.txt` | `15_lma_touschek/lma/lma_support_kekls.out.txt` | 出力 |  |
| `v2_4/A4_kekls_level1.out.txt` | `15_lma_touschek/lma/lma_support_kekls_level1.out.txt` | 出力 |  |
| `v2_4/A4_lma_support.py` | `15_lma_touschek/lma/lma_support.py` | スクリプト |  |
| `v2_4/A5_astrid2.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_fmap_grid_astrid2.out.txt` | 出力 |  |
| `v2_4/A5_astrid2_D8.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_fmap_grid_astrid2_D8.out.txt` | 出力 |  |
| `v2_4/A5_astrid2_D8_skip.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_fmap_grid_astrid2_D8_skip.out.txt` | 出力 |  |
| `v2_4/A5_kekls.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_fmap_grid_kekls.out.txt` | 出力 |  |
| `v2_4/A5_nu_of_J_grid.py` | `14_nonlinear_dynamics/normal_form/nu_of_J_fmap_grid.py` | スクリプト |  |
| `v2_4/A7_ampnf_kekls.py` | `14_nonlinear_dynamics/normal_form/ampnf_kekls.py` | スクリプト |  |
| `v2_4/A7_kekls.out.txt` | `14_nonlinear_dynamics/normal_form/ampnf_kekls.out.txt` | 出力 |  |
| `v2_4/F1_astrid.out.txt` | `14_nonlinear_dynamics/normal_form/rotated_bend_astrid.out.txt` | 出力 |  |
| `v2_4/F1_fodo.out.txt` | `14_nonlinear_dynamics/normal_form/rotated_bend_fodo.out.txt` | 出力 |  |
| `v2_4/F1_rotated_bend.py` | `14_nonlinear_dynamics/normal_form/rotated_bend.py` | スクリプト |  |
| `v2_4/G_independent_flow.out.txt` | `14_nonlinear_dynamics/normal_form/k0_independent_flow.out.txt` | 出力 |  |
| `v2_4/G_independent_flow.py` | `14_nonlinear_dynamics/normal_form/k0_independent_flow.py` | スクリプト |  |
| `v2_4/G_k0_ampnf.out.txt` | `14_nonlinear_dynamics/normal_form/k0_ampnf.out.txt` | 出力 |  |
| `v2_4/G_k0_ampnf.py` | `14_nonlinear_dynamics/normal_form/k0_ampnf.py` | スクリプト |  |
| `v2_4/G_real_scale.out.txt` | `14_nonlinear_dynamics/normal_form/k0_real_scale.out.txt` | 出力 |  |
| `v2_4/G_real_scale.py` | `14_nonlinear_dynamics/normal_form/k0_real_scale.py` | スクリプト |  |
| `v2_4/H1_astrid.out.txt` | `14_nonlinear_dynamics/normal_form/coupled_h4_astrid.out.txt` | 出力 |  |
| `v2_4/H1_coupled_h4.py` | `14_nonlinear_dynamics/normal_form/coupled_h4.py` | スクリプト |  |
| `v2_4/H1_fodo.out.txt` | `14_nonlinear_dynamics/normal_form/coupled_h4_fodo.out.txt` | 出力 |  |
| `v2_4/H2_kekls12_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_corr_s17.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_corr_s17.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_ctrl_expand_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_ctrl_expand_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_corr_s17.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_corr_s17.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_corr_s3.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_corr_s3.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_ctrl_expand_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_ctrl_expand_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t1964_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t1964_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t1964_corr_s42.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t1964_corr_s42.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t1964_corr_s5.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t1964_corr_s5.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t1964_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t1964_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t3928_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t3928_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t3928_corr_s42.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t3928_corr_s42.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s17.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s17.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s3.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s3.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s42.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s42.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s5.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s5.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_corr_s99.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_corr_s99.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_t982_p9_uncorr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_t982_p9_uncorr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_rad_uncorr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_rad_uncorr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls12_t982_summary.out.txt` | `15_lma_touschek/lma/lma_corrected_kekls_v03_12_t982_summary.out.txt` | 出力 |  |
| `v2_4/H2_kekls12_uncorr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_12_uncorr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_corr_s111_turns856.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_corr_s111_turns856.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_corr_s17.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_corr_s17.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_corr_s3.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_corr_s3.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_ctrl_expand_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_ctrl_expand_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_lma_corrected.template.sad` | `15_lma_touschek/lma/lma_corrected_kekls_v03_68.template.sad` | スクリプト |  |
| `v2_4/H2_kekls_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_noerr_turns856.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_noerr_turns856.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_rad_corr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_rad_corr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_rad_ctrl_expand_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_rad_ctrl_expand_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_rad_noerr.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_rad_noerr.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_rad_uncorr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_rad_uncorr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_uncorr_s111.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_uncorr_s111.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_uncorr_s111_turns856.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_uncorr_s111_turns856.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_x1.5_s111_fail.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_x1.5_s111_fail.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_kekls_x2_s111_fail.out.txt` | `_deletable/lma_corrected_kekls/lma_corrected_kekls_v03_68_x2_s111_fail.out.txt` | 出力 | `_deletable/` に移した (削除してよい) |
| `v2_4/H2_run.py` | `15_lma_touschek/lma/lma_corrected_run.py` | スクリプト |  |
| `v2_4/H2_summary_t982.py` | `15_lma_touschek/lma/lma_corrected_summary.py` | スクリプト |  |
| `v2_4/H2b_kekls_v03_12_lma_corrected.template.sad` | `15_lma_touschek/lma/lma_corrected_kekls_v03_12.template.sad` | スクリプト |  |
| `v2_4/I1_astrid.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_orders_astrid.out.txt` | 出力 |  |
| `v2_4/I1_fodo.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_orders_fodo.out.txt` | 出力 |  |
| `v2_4/I1_kekls.out.txt` | `14_nonlinear_dynamics/normal_form/nu_of_J_orders_kekls.out.txt` | 出力 |  |
| `v2_4/I1_nu_of_J_orders.py` | `14_nonlinear_dynamics/normal_form/nu_of_J_orders.py` | スクリプト |  |
| `v2_4/I2_island_fodo.out.txt` | `14_nonlinear_dynamics/normal_form/island_vs_tracking_fodo.out.txt` | 出力 |  |
| `v2_4/I2_island_vs_tracking.py` | `14_nonlinear_dynamics/normal_form/island_vs_tracking.py` | スクリプト |  |
| `v2_4/README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 14_nonlinear_dynamics/normal_form, 15_lma_touschek/lma, 15_lma_touschek/touschek_elegant |
| `v2_B/B1_special_functions.pysad.txt` | `10_math_data/fn_SpecialFunctions_gs_errors.pysad.txt` | 出力 |  |
| `v2_B/B1_special_functions.sad` | `10_math_data/fn_SpecialFunctions_gs_errors.sad` | スクリプト |  |
| `v2_B/B2_package_context.pysad.txt` | `08_scripting/syntax_package_context.pysad.txt` | 出力 |  |
| `v2_B/B2_package_context.sad` | `08_scripting/syntax_package_context.sad` | スクリプト |  |
| `v2_B/B3_shared_system.pysad.txt` | `09_io_system/fn_SharedMemory_Parallelize_System_Fork.pysad.txt` | 出力 |  |
| `v2_B/B3_shared_system.sad` | `09_io_system/fn_SharedMemory_Parallelize_System_Fork.sad` | スクリプト |  |
| `v2_B/B4_ffs_v2.pysad.txt` | `01_lattice_elements/FFS_small_options.pysad.txt` | 出力 |  |
| `v2_B/B4_ffs_v2.sad` | `01_lattice_elements/FFS_small_options.sad` | スクリプト |  |
| `v2_B/README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 08_scripting, 09_io_system, 10_math_data |
| `v2_C/C1_undulator.pysad.txt` | `01_lattice_elements/elem_UND.pysad.txt` | 出力 |  |
| `v2_C/C1_undulator.sad` | `01_lattice_elements/elem_UND.sad` | スクリプト |  |
| `v2_C/C2_phsrot.pysad.txt` | `01_lattice_elements/elem_PHSROT.pysad.txt` | 出力 |  |
| `v2_C/C2_phsrot.sad` | `01_lattice_elements/elem_PHSROT.sad` | スクリプト |  |
| `v2_C/C2a_sol_tilt.pysad.txt` | `07_spin_solenoid_special/elem_SOL_tilt_offset.pysad.txt` | 出力 |  |
| `v2_C/C2a_sol_tilt.sad` | `07_spin_solenoid_special/elem_SOL_tilt_offset.sad` | スクリプト |  |
| `v2_C/C2b_sol_rf.pysad.txt` | `07_spin_solenoid_special/elem_SOL_RF_TRPT.pysad.txt` | 出力 |  |
| `v2_C/C2b_sol_rf.sad` | `07_spin_solenoid_special/elem_SOL_RF_TRPT.sad` | スクリプト |  |
| `v2_C/C2c_spin.pysad.txt` | `16_spin_polarization/spin_tracking_sol_und.pysad.txt` | 出力 |  |
| `v2_C/C2c_spin.sad` | `16_spin_polarization/spin_tracking_sol_und.sad` | スクリプト |  |
| `v2_C/C2d_taper_optics.pysad.txt` | `04_emittance_radiation/flag_RADTAPER_optics.pysad.txt` | 出力 |  |
| `v2_C/C2d_taper_optics.sad` | `04_emittance_radiation/flag_RADTAPER_optics.sad` | スクリプト |  |
| `v2_C/C2m_putget_measure.pysad.txt` | `02_optics_cod/fn_PUTMEASURE_GETMEASURE.pysad.txt` | 出力 |  |
| `v2_C/C2m_putget_measure.sad` | `02_optics_cod/fn_PUTMEASURE_GETMEASURE.sad` | スクリプト |  |
| `v2_C/C2m_trcod_reference.pysad.txt` | `02_optics_cod/flag_TRCOD_reference.pysad.txt` | 出力 |  |
| `v2_C/C2m_trcod_reference.sad` | `02_optics_cod/flag_TRCOD_reference.sad` | スクリプト |  |
| `v2_C/C2s_mult_radius.pysad.txt` | `05_tracking/flag_SPAC_MULT_radius.pysad.txt` | 出力 |  |
| `v2_C/C2s_mult_radius.sad` | `05_tracking/flag_SPAC_MULT_radius.sad` | スクリプト |  |
| `v2_C/C2s_wspac.pysad.txt` | `05_tracking/flag_WSPAC_tune_shift.pysad.txt` | 出力 |  |
| `v2_C/C2s_wspac.sad` | `05_tracking/flag_WSPAC_tune_shift.sad` | スクリプト |  |
| `v2_C/C2w_delcor_rms.pysad.txt` | `02_optics_cod/flag_DELCOR_rms.pysad.txt` | 出力 |  |
| `v2_C/C2w_delcor_rms.sad` | `02_optics_cod/flag_DELCOR_rms.sad` | スクリプト |  |
| `v2_C/C2w_main_script.pysad.txt` | `01_lattice_elements/syntax_MAIN_assign.pysad.txt` | 出力 |  |
| `v2_C/C2w_main_script.sad` | `01_lattice_elements/syntax_MAIN_assign.sad` | スクリプト |  |
| `v2_C/C2w_wake_table.pysad.txt` | `06_collective/flag_WAKE_command_table.pysad.txt` | 出力 |  |
| `v2_C/C2w_wake_table.sad` | `06_collective/flag_WAKE_command_table.sad` | スクリプト |  |
| `v2_C/C3_sol_region.pysad.txt` | `07_spin_solenoid_special/elem_SOL_region_elements.pysad.txt` | 出力 |  |
| `v2_C/C3_sol_region.sad` | `07_spin_solenoid_special/elem_SOL_region_elements.sad` | スクリプト |  |
| `v2_C/C4_loss_and_taper.pysad.txt` | `05_tracking/flag_LOSSMAP_RADTAPER.pysad.txt` | 出力 |  |
| `v2_C/C4_loss_and_taper.sad` | `05_tracking/flag_LOSSMAP_RADTAPER.sad` | スクリプト |  |
| `v2_C/C5_radlight.pysad.txt` | `07_spin_solenoid_special/flag_RADLIGHT_orbit_field.pysad.txt` | 出力 |  |
| `v2_C/C5_radlight.sad` | `07_spin_solenoid_special/flag_RADLIGHT_orbit_field.sad` | スクリプト |  |
| `v2_C/C6_misc.pysad.txt` | `02_optics_cod/FFS_misc_DELCOR_WAKE_MULT_skew.pysad.txt` | 出力 |  |
| `v2_C/C6_misc.sad` | `02_optics_cod/FFS_misc_DELCOR_WAKE_MULT_skew.sad` | スクリプト |  |
| `v2_C/README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 02_optics_cod, 04_emittance_radiation, 05_tracking, 06_collective, 07_spin_solenoid_special, 16_spin_polarization |
| `v2_next/A1_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 02_optics_cod, 04_emittance_radiation, 05_tracking |
| `v2_next/A1_bend_f1_adjacent.gs.txt` | `01_lattice_elements/elem_BEND_F1_adjacent.gs.txt` | 出力 |  |
| `v2_next/A1_bend_f1_adjacent.pysad.txt` | `01_lattice_elements/elem_BEND_F1_adjacent.pysad.txt` | 出力 |  |
| `v2_next/A1_bend_f1_adjacent.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_bend_f1_adjacent.sad` | `01_lattice_elements/elem_BEND_F1_adjacent.sad` | スクリプト |  |
| `v2_next/A1_cavity_offset.gs.txt` | `01_lattice_elements/elem_CAVI_offset.gs.txt` | 出力 |  |
| `v2_next/A1_cavity_offset.pysad.txt` | `01_lattice_elements/elem_CAVI_offset.pysad.txt` | 出力 |  |
| `v2_next/A1_cavity_offset.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_cavity_offset.sad` | `01_lattice_elements/elem_CAVI_offset.sad` | スクリプト |  |
| `v2_next/A1_n12_ins_detr.gs.txt` | `02_optics_cod/flag_INS_DETR.gs.txt` | 出力 |  |
| `v2_next/A1_n12_ins_detr.pysad.txt` | `02_optics_cod/flag_INS_DETR.pysad.txt` | 出力 |  |
| `v2_next/A1_n12_ins_detr.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_n12_ins_detr.sad` | `02_optics_cod/flag_INS_DETR.sad` | スクリプト |  |
| `v2_next/A1_proton_ring.gs.txt` | `04_emittance_radiation/var_MASS_CHARGE_proton_ring.gs.txt` | 出力 |  |
| `v2_next/A1_proton_ring.pysad.txt` | `04_emittance_radiation/var_MASS_CHARGE_proton_ring.pysad.txt` | 出力 |  |
| `v2_next/A1_proton_ring.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_proton_ring.sad` | `04_emittance_radiation/var_MASS_CHARGE_proton_ring.sad` | スクリプト |  |
| `v2_next/A1_pspac_resolution.gs.txt` | `05_tracking/flag_PSPAC_resolution.gs.txt` | 出力 |  |
| `v2_next/A1_pspac_resolution.pysad.txt` | `05_tracking/flag_PSPAC_resolution.pysad.txt` | 出力 |  |
| `v2_next/A1_pspac_resolution.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_pspac_resolution.sad` | `05_tracking/flag_PSPAC_resolution.sad` | スクリプト |  |
| `v2_next/A1_radcod_closure.gs.txt` | `04_emittance_radiation/flag_RADCOD_closure.gs.txt` | 出力 |  |
| `v2_next/A1_radcod_closure.pysad.txt` | `04_emittance_radiation/flag_RADCOD_closure.pysad.txt` | 出力 |  |
| `v2_next/A1_radcod_closure.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_radcod_closure.sad` | `04_emittance_radiation/flag_RADCOD_closure.sad` | スクリプト |  |
| `v2_next/A1_use_after_emit.gs.txt` | `01_lattice_elements/FFS_USE_after_EMIT.gs.txt` | 出力 |  |
| `v2_next/A1_use_after_emit.pysad.txt` | `01_lattice_elements/FFS_USE_after_EMIT.pysad.txt` | 出力 |  |
| `v2_next/A1_use_after_emit.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/A1_use_after_emit.sad` | `01_lattice_elements/FFS_USE_after_EMIT.sad` | スクリプト |  |
| `v2_next/A2_dx_bend.gs.txt` | `01_lattice_elements/elem_BEND_DX_offset.gs.txt` | 出力 |  |
| `v2_next/A2_dx_bend.pysad.txt` | `01_lattice_elements/elem_BEND_DX_offset.pysad.txt` | 出力 |  |
| `v2_next/A2_dx_bend.sad` | `01_lattice_elements/elem_BEND_DX_offset.sad` | スクリプト |  |
| `v2_next/A2_i4v_coupled.gs.txt` | `04_emittance_radiation/FFS_RADINT_I4V_coupled.gs.txt` | 出力 |  |
| `v2_next/A2_i4v_coupled.pysad.txt` | `04_emittance_radiation/FFS_RADINT_I4V_coupled.pysad.txt` | 出力 |  |
| `v2_next/A2_i4v_coupled.sad` | `04_emittance_radiation/FFS_RADINT_I4V_coupled.sad` | スクリプト |  |
| `v2_next/A2_pol_fluc_tracking.py` | `16_spin_polarization/pol_fluc_tracking.py` | スクリプト |  |
| `v2_next/A2_pol_mpw.gs.txt` | `16_spin_polarization/pol_mpw.gs.txt` | 出力 |  |
| `v2_next/A2_pol_mpw.pysad.txt` | `16_spin_polarization/pol_mpw.pysad.txt` | 出力 |  |
| `v2_next/A2_pol_mpw.sad` | `16_spin_polarization/pol_mpw.sad` | スクリプト |  |
| `v2_next/A2_radint_pf.gs.txt` | `04_emittance_radiation/FFS_RADINT_pf.gs.txt` | 出力 |  |
| `v2_next/A2_radint_pf.pysad.txt` | `04_emittance_radiation/FFS_RADINT_pf.pysad.txt` | 出力 |  |
| `v2_next/A2_radint_pf.sad` | `04_emittance_radiation/FFS_RADINT_pf.sad` | スクリプト |  |
| `v2_next/A2_spin_turn.gs.txt` | `16_spin_polarization/spin_turn.gs.txt` | 出力 |  |
| `v2_next/A2_spin_turn.pysad.txt` | `16_spin_polarization/spin_turn.pysad.txt` | 出力 |  |
| `v2_next/A2_spin_turn.sad` | `16_spin_polarization/spin_turn.sad` | スクリプト |  |
| `v2_next/B1_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 04_emittance_radiation, 05_tracking, 07_spin_solenoid_special, 16_spin_polarization |
| `v2_next/B1_bend_f1_emittance.gs.txt` | `04_emittance_radiation/elem_BEND_F1_emittance.gs.txt` | 出力 |  |
| `v2_next/B1_bend_f1_emittance.pysad.txt` | `04_emittance_radiation/elem_BEND_F1_emittance.pysad.txt` | 出力 |  |
| `v2_next/B1_bend_f1_emittance.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B1_bend_f1_emittance.sad` | `04_emittance_radiation/elem_BEND_F1_emittance.sad` | スクリプト |  |
| `v2_next/B1_k1_bend_edge.out.txt` | `01_lattice_elements/elem_BEND_K1_edge.out.txt` | 出力 |  |
| `v2_next/B1_k1_bend_edge.py` | `01_lattice_elements/elem_BEND_K1_edge.py` | スクリプト |  |
| `v2_next/B1_sol_calc6d_intra.gs.txt` | `07_spin_solenoid_special/elem_SOL_calc6d_intra.gs.txt` | 出力 |  |
| `v2_next/B1_sol_calc6d_intra.pysad.txt` | `07_spin_solenoid_special/elem_SOL_calc6d_intra.pysad.txt` | 出力 |  |
| `v2_next/B1_sol_calc6d_intra.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B1_sol_calc6d_intra.sad` | `07_spin_solenoid_special/elem_SOL_calc6d_intra.sad` | スクリプト |  |
| `v2_next/B1_sol_spin_bigbz.out.txt` | `16_spin_polarization/spin_sol_bigbz.out.txt` | 出力 |  |
| `v2_next/B1_sol_spin_bigbz.py` | `16_spin_polarization/spin_sol_bigbz.py` | スクリプト |  |
| `v2_next/B1_wspac_emittance.gs.txt` | `05_tracking/flag_WSPAC_emittance.gs.txt` | 出力 |  |
| `v2_next/B1_wspac_emittance.pysad.txt` | `05_tracking/flag_WSPAC_emittance.pysad.txt` | 出力 |  |
| `v2_next/B1_wspac_emittance.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B1_wspac_emittance.sad` | `05_tracking/flag_WSPAC_emittance.sad` | スクリプト |  |
| `v2_next/B2_context.gs.txt` | `08_scripting/syntax_context.gs.txt` | 出力 |  |
| `v2_next/B2_context.pysad.txt` | `08_scripting/syntax_context.pysad.txt` | 出力 |  |
| `v2_next/B2_context.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_context.sad` | `08_scripting/syntax_context.sad` | スクリプト |  |
| `v2_next/B2_free_component.gs.txt` | `03_matching/flag_FREE_component.gs.txt` | 出力 |  |
| `v2_next/B2_free_component.pysad.txt` | `03_matching/flag_FREE_component.pysad.txt` | 出力 |  |
| `v2_next/B2_free_component.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_free_component.sad` | `03_matching/flag_FREE_component.sad` | スクリプト |  |
| `v2_next/B2_packages_compare.gs.txt` | `13_packages/fn_packages_compare.gs.txt` | 出力 |  |
| `v2_next/B2_packages_compare.pysad.txt` | `13_packages/fn_packages_compare.pysad.txt` | 出力 |  |
| `v2_next/B2_packages_compare.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_packages_compare.sad` | `13_packages/fn_packages_compare.sad` | スクリプト |  |
| `v2_next/B2_parallelize.gs.txt` | `09_io_system/fn_Parallelize.gs.txt` | 出力 |  |
| `v2_next/B2_parallelize.pysad.txt` | `09_io_system/fn_Parallelize.pysad.txt` | 出力 |  |
| `v2_next/B2_parallelize.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_parallelize.pysad_wsl.txt` | `09_io_system/fn_Parallelize.pysad_wsl.txt` | 出力 |  |
| `v2_next/B2_parallelize.sad` | `09_io_system/fn_Parallelize.sad` | スクリプト |  |
| `v2_next/B2_special_functions_complex.gs.txt` | `10_math_data/fn_SpecialFunctions_complex.gs.txt` | 出力 |  |
| `v2_next/B2_special_functions_complex.pysad.txt` | `10_math_data/fn_SpecialFunctions_complex.pysad.txt` | 出力 |  |
| `v2_next/B2_special_functions_complex.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_special_functions_complex.sad` | `10_math_data/fn_SpecialFunctions_complex.sad` | スクリプト |  |
| `v2_next/B2_symbol_message.gs.txt` | `08_scripting/syntax_symbol_message.gs.txt` | 出力 |  |
| `v2_next/B2_symbol_message.pysad.txt` | `08_scripting/syntax_symbol_message.pysad.txt` | 出力 |  |
| `v2_next/B2_symbol_message.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_symbol_message.sad` | `08_scripting/syntax_symbol_message.sad` | スクリプト |  |
| `v2_next/B2_visit_bye.gs.txt` | `08_scripting/FFS_VISIT_BYE.gs.txt` | 出力 |  |
| `v2_next/B2_visit_bye.pysad.txt` | `08_scripting/FFS_VISIT_BYE.pysad.txt` | 出力 |  |
| `v2_next/B2_visit_bye.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B2_visit_bye.sad` | `08_scripting/FFS_VISIT_BYE.sad` | スクリプト |  |
| `v2_next/B3_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 03_matching, 05_tracking, 07_spin_solenoid_special, 16_spin_polarization |
| `v2_next/B3_geofit_quad_dx.gs.txt` | `03_matching/FFS_FIT_geometry_quad_dx.gs.txt` | 出力 |  |
| `v2_next/B3_geofit_quad_dx.pysad.txt` | `03_matching/FFS_FIT_geometry_quad_dx.pysad.txt` | 出力 |  |
| `v2_next/B3_geofit_quad_dx.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B3_geofit_quad_dx.sad` | `03_matching/FFS_FIT_geometry_quad_dx.sad` | スクリプト |  |
| `v2_next/B3_naff_tunes.out.txt` | `05_tracking/fn_NAFF_tunes.out.txt` | 出力 |  |
| `v2_next/B3_naff_tunes.py` | `05_tracking/fn_NAFF_tunes.py` | スクリプト |  |
| `v2_next/B3_photons_bend.gs.txt` | `07_spin_solenoid_special/flag_PHOTONS_bend.gs.txt` | 出力 |  |
| `v2_next/B3_photons_bend.pysad.txt` | `07_spin_solenoid_special/flag_PHOTONS_bend.pysad.txt` | 出力 |  |
| `v2_next/B3_photons_bend.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B3_photons_bend.sad` | `07_spin_solenoid_special/flag_PHOTONS_bend.sad` | スクリプト |  |
| `v2_next/B3_pol_sideband.gs.txt` | `16_spin_polarization/pol_sideband.gs.txt` | 出力 |  |
| `v2_next/B3_pol_sideband.pysad.txt` | `16_spin_polarization/pol_sideband.pysad.txt` | 出力 |  |
| `v2_next/B3_pol_sideband.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/B3_pol_sideband.sad` | `16_spin_polarization/pol_sideband.sad` | スクリプト |  |
| `v2_next/B3_radpol_roundoff.py` | `16_spin_polarization/radpol_roundoff.py` | スクリプト |  |
| `v2_next/B3_radpol_roundoff.pysad.txt` | `16_spin_polarization/radpol_roundoff.pysad.txt` | 出力 |  |
| `v2_next/B3_spin_mc.py` | `16_spin_polarization/pol_sideband_spin_mc.py` | スクリプト |  |
| `v2_next/C3_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 02_optics_cod, 04_emittance_radiation, 05_tracking, 07_spin_solenoid_special |
| `v2_next/C3_fluc_nparallel.gs.txt` | `05_tracking/flag_FLUC_nparallel.gs.txt` | 出力 |  |
| `v2_next/C3_fluc_nparallel.pysad.txt` | `05_tracking/flag_FLUC_nparallel.pysad.txt` | 出力 |  |
| `v2_next/C3_fluc_nparallel.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_fluc_nparallel.sad` | `05_tracking/flag_FLUC_nparallel.sad` | スクリプト |  |
| `v2_next/C3_radfield_bend.gs.txt` | `07_spin_solenoid_special/fn_RadiationField_bend.gs.txt` | 出力 |  |
| `v2_next/C3_radfield_bend.pysad.txt` | `07_spin_solenoid_special/fn_RadiationField_bend.pysad.txt` | 出力 |  |
| `v2_next/C3_radfield_bend.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_radfield_bend.sad` | `07_spin_solenoid_special/fn_RadiationField_bend.sad` | スクリプト |  |
| `v2_next/C3_radint_k1_bend_vdisp.gs.txt` | `04_emittance_radiation/FFS_RADINT_k1_bend_vdisp.gs.txt` | 出力 |  |
| `v2_next/C3_radint_k1_bend_vdisp.pysad.txt` | `04_emittance_radiation/FFS_RADINT_k1_bend_vdisp.pysad.txt` | 出力 |  |
| `v2_next/C3_radint_k1_bend_vdisp.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_radint_k1_bend_vdisp.sad` | `04_emittance_radiation/FFS_RADINT_k1_bend_vdisp.sad` | スクリプト |  |
| `v2_next/C3_radint_quad_offset.gs.txt` | `04_emittance_radiation/FFS_RADINT_quad_offset.gs.txt` | 出力 |  |
| `v2_next/C3_radint_quad_offset.pysad.txt` | `04_emittance_radiation/FFS_RADINT_quad_offset.pysad.txt` | 出力 |  |
| `v2_next/C3_radint_quad_offset.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_radint_quad_offset.sad` | `04_emittance_radiation/FFS_RADINT_quad_offset.sad` | スクリプト |  |
| `v2_next/C3_synchrobeta_coupled.gs.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_coupled.gs.txt` | 出力 |  |
| `v2_next/C3_synchrobeta_coupled.pysad.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_coupled.pysad.txt` | 出力 |  |
| `v2_next/C3_synchrobeta_coupled.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_synchrobeta_coupled.sad` | `04_emittance_radiation/fn_SynchroBetaEmittance_coupled.sad` | スクリプト |  |
| `v2_next/C3_trcod_orbit.gs.txt` | `02_optics_cod/flag_TRCOD_orbit.gs.txt` | 出力 |  |
| `v2_next/C3_trcod_orbit.pysad.txt` | `02_optics_cod/flag_TRCOD_orbit.pysad.txt` | 出力 |  |
| `v2_next/C3_trcod_orbit.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/C3_trcod_orbit.sad` | `02_optics_cod/flag_TRCOD_orbit.sad` | スクリプト |  |
| `v2_next/D1_norfsw_orbit.gs.txt` | `04_emittance_radiation/flag_NORFSW_orbit.gs.txt` | 出力 |  |
| `v2_next/D1_norfsw_orbit.pysad.txt` | `04_emittance_radiation/flag_NORFSW_orbit.pysad.txt` | 出力 |  |
| `v2_next/D1_norfsw_orbit.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D1_norfsw_orbit.sad` | `04_emittance_radiation/flag_NORFSW_orbit.sad` | スクリプト |  |
| `v2_next/D1_photons_bb_bend.gs.txt` | `07_spin_solenoid_special/flag_PHOTONS_beambeam_bend.gs.txt` | 出力 |  |
| `v2_next/D1_photons_bb_bend.pysad.txt` | `07_spin_solenoid_special/flag_PHOTONS_beambeam_bend.pysad.txt` | 出力 |  |
| `v2_next/D1_photons_bb_bend.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D1_photons_bb_bend.sad` | `07_spin_solenoid_special/flag_PHOTONS_beambeam_bend.sad` | スクリプト |  |
| `v2_next/D1_photons_mult.gs.txt` | `07_spin_solenoid_special/flag_PHOTONS_mult.gs.txt` | 出力 |  |
| `v2_next/D1_photons_mult.pysad.txt` | `07_spin_solenoid_special/flag_PHOTONS_mult.pysad.txt` | 出力 |  |
| `v2_next/D1_photons_mult.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D1_photons_mult.sad` | `07_spin_solenoid_special/flag_PHOTONS_mult.sad` | スクリプト |  |
| `v2_next/D1_spin_accel.out.txt` | `16_spin_polarization/spin_accel.out.txt` | 出力 |  |
| `v2_next/D1_spin_accel.py` | `16_spin_polarization/spin_accel.py` | スクリプト |  |
| `v2_next/D2_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 13_packages |
| `v2_next/D2_correctlocal.gs.txt` | `13_packages/fn_CorrectLocal.gs.txt` | 出力 |  |
| `v2_next/D2_correctlocal.pysad.txt` | `13_packages/fn_CorrectLocal.pysad.txt` | 出力 |  |
| `v2_next/D2_correctlocal.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D2_correctlocal.sad` | `13_packages/fn_CorrectLocal.sad` | スクリプト |  |
| `v2_next/D2_csroymatrix.gs.txt` | `13_packages/fn_CSROYMatrix.gs.txt` | 出力 |  |
| `v2_next/D2_csroymatrix.pysad.txt` | `13_packages/fn_CSROYMatrix.pysad.txt` | 出力 |  |
| `v2_next/D2_csroymatrix.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D2_csroymatrix.sad` | `13_packages/fn_CSROYMatrix.sad` | スクリプト |  |
| `v2_next/D2_main_expr.gs.txt` | `01_lattice_elements/syntax_MAIN_expression.gs.txt` | 出力 |  |
| `v2_next/D2_main_expr.pysad.txt` | `01_lattice_elements/syntax_MAIN_expression.pysad.txt` | 出力 |  |
| `v2_next/D2_main_expr.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D2_main_expr.sad` | `01_lattice_elements/syntax_MAIN_expression.sad` | スクリプト |  |
| `v2_next/D2_micado.gs.txt` | `13_packages/fn_CorrectOrbit_Micado.gs.txt` | 出力 |  |
| `v2_next/D2_micado.pysad.txt` | `13_packages/fn_CorrectOrbit_Micado.pysad.txt` | 出力 |  |
| `v2_next/D2_micado.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D2_micado.sad` | `13_packages/fn_CorrectOrbit_Micado.sad` | スクリプト |  |
| `v2_next/D2_periodicr.gs.txt` | `13_packages/fn_PeriodicR.gs.txt` | 出力 |  |
| `v2_next/D2_periodicr.pysad.txt` | `13_packages/fn_PeriodicR.pysad.txt` | 出力 |  |
| `v2_next/D2_periodicr.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/D2_periodicr.sad` | `13_packages/fn_PeriodicR.sad` | スクリプト |  |
| `v2_next/E_b108_poles.out.txt` | `01_lattice_elements/elem_BEND_steer_AE_poles.out.txt` | 出力 |  |
| `v2_next/E_b108_poles.py` | `01_lattice_elements/elem_BEND_steer_AE_poles.py` | スクリプト |  |
| `v2_next/E_mult_f1.gs.txt` | `01_lattice_elements/elem_MULT_F1.gs.txt` | 出力 |  |
| `v2_next/E_mult_f1.pysad.txt` | `01_lattice_elements/elem_MULT_F1.pysad.txt` | 出力 |  |
| `v2_next/E_mult_f1.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/E_mult_f1.sad` | `01_lattice_elements/elem_MULT_F1.sad` | スクリプト |  |
| `v2_next/E_pol_kick_response.out.txt` | `16_spin_polarization/pol_kick_response.out.txt` | 出力 |  |
| `v2_next/E_pol_kick_response.py` | `16_spin_polarization/pol_kick_response.py` | スクリプト |  |
| `v2_next/E_sol_tilt_spin.out.txt` | `16_spin_polarization/spin_sol_tilt.out.txt` | 出力 |  |
| `v2_next/E_sol_tilt_spin.py` | `16_spin_polarization/spin_sol_tilt.py` | スクリプト |  |
| `v2_next/E_steer_f1.gs.txt` | `01_lattice_elements/elem_BEND_steer_F1.gs.txt` | 出力 |  |
| `v2_next/E_steer_f1.pysad.txt` | `01_lattice_elements/elem_BEND_steer_F1.pysad.txt` | 出力 |  |
| `v2_next/E_steer_f1.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/E_steer_f1.sad` | `01_lattice_elements/elem_BEND_steer_F1.sad` | スクリプト |  |
| `v2_next/F2_photons_coord_tcavi.gs.txt` | `07_spin_solenoid_special/flag_PHOTONS_coord_tcavi.gs.txt` | 出力 |  |
| `v2_next/F2_photons_coord_tcavi.sad` | `07_spin_solenoid_special/flag_PHOTONS_coord_tcavi.sad` | スクリプト |  |
| `v2_next/F4_sol_bound0_track.gs.txt` | `07_spin_solenoid_special/elem_SOL_bound0_track.gs.txt` | 出力 |  |
| `v2_next/F4_sol_bound0_track.sad` | `07_spin_solenoid_special/elem_SOL_bound0_track.sad` | スクリプト |  |
| `v2_next/P1_bend_vertical_angle_spin.out.txt` | `16_spin_polarization/spin_bend_vertical_angle.out.txt` | 出力 |  |
| `v2_next/P1_bend_vertical_angle_spin.py` | `16_spin_polarization/spin_bend_vertical_angle.py` | スクリプト |  |
| `v2_next/P1_pol_chromatic.out.txt` | `16_spin_polarization/pol_chromatic.out.txt` | 出力 |  |
| `v2_next/P1_pol_chromatic.py` | `16_spin_polarization/pol_chromatic.py` | スクリプト |  |
| `v2_next/P1_pol_mpw_rotated.gs.txt` | `16_spin_polarization/pol_mpw_rotated.gs.txt` | 出力 |  |
| `v2_next/P1_pol_mpw_rotated.pysad.txt` | `16_spin_polarization/pol_mpw_rotated.pysad.txt` | 出力 |  |
| `v2_next/P1_pol_mpw_rotated.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P1_pol_mpw_rotated.sad` | `16_spin_polarization/pol_mpw_rotated.sad` | スクリプト |  |
| `v2_next/P1_pol_reversed_bend.gs.txt` | `16_spin_polarization/pol_reversed_bend.gs.txt` | 出力 |  |
| `v2_next/P1_pol_reversed_bend.pysad.txt` | `16_spin_polarization/pol_reversed_bend.pysad.txt` | 出力 |  |
| `v2_next/P1_pol_reversed_bend.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P1_pol_reversed_bend.sad` | `16_spin_polarization/pol_reversed_bend.sad` | スクリプト |  |
| `v2_next/P1_pol_vertical_bend.gs.txt` | `16_spin_polarization/pol_vertical_bend.gs.txt` | 出力 |  |
| `v2_next/P1_pol_vertical_bend.pysad.txt` | `16_spin_polarization/pol_vertical_bend.pysad.txt` | 出力 |  |
| `v2_next/P1_pol_vertical_bend.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P1_pol_vertical_bend.sad` | `16_spin_polarization/pol_vertical_bend.sad` | スクリプト |  |
| `v2_next/P1_response_samples_f1_M1M2chr.json` | `16_spin_polarization/pol_chromatic_response_samples_f1.json` | データ |  |
| `v2_next/P1_response_samples_f2_models.json` | `16_spin_polarization/pol_chromatic_response_samples_f2_models.json` | データ |  |
| `v2_next/P1_st_dk_check.out.txt` | `16_spin_polarization/pol_st_dk_check.out.txt` | 出力 |  |
| `v2_next/P1_st_dk_check.py` | `16_spin_polarization/pol_st_dk_check.py` | スクリプト |  |
| `v2_next/P2_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 04_emittance_radiation |
| `v2_next/P2_bend_k1_face.out.txt` | `01_lattice_elements/elem_BEND_K1_face.out.txt` | 出力 |  |
| `v2_next/P2_bend_k1_face.py` | `01_lattice_elements/elem_BEND_K1_face.py` | スクリプト |  |
| `v2_next/P2_emitdiv_rect.out.txt` | `04_emittance_radiation/flag_EMITDIV_rect.out.txt` | 出力 |  |
| `v2_next/P2_emitdiv_rect.py` | `04_emittance_radiation/flag_EMITDIV_rect.py` | スクリプト |  |
| `v2_next/P2_radcod_derivative.out.txt` | `04_emittance_radiation/flag_RADCOD_derivative.out.txt` | 出力 |  |
| `v2_next/P2_radcod_derivative.py` | `04_emittance_radiation/flag_RADCOD_derivative.py` | スクリプト |  |
| `v2_next/P2_sb_dispersion_start.gs.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_dispersion_start.gs.txt` | 出力 |  |
| `v2_next/P2_sb_dispersion_start.pysad.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_dispersion_start.pysad.txt` | 出力 |  |
| `v2_next/P2_sb_dispersion_start.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P2_sb_dispersion_start.sad` | `04_emittance_radiation/fn_SynchroBetaEmittance_dispersion_start.sad` | スクリプト |  |
| `v2_next/P2_sb_ring_a.gs.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_rings_ring_a.gs.txt` | 出力 |  |
| `v2_next/P2_sb_rings.out.txt` | `04_emittance_radiation/fn_SynchroBetaEmittance_rings.out.txt` | 出力 |  |
| `v2_next/P2_sb_rings.py` | `04_emittance_radiation/fn_SynchroBetaEmittance_rings.py` | スクリプト |  |
| `v2_next/P2_tradke_exact.out.txt` | `04_emittance_radiation/flag_RAD_tradke_exact.out.txt` | 出力 |  |
| `v2_next/P2_tradke_exact.py` | `04_emittance_radiation/flag_RAD_tradke_exact.py` | スクリプト |  |
| `v2_next/P3_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements |
| `v2_next/P3_b108_exact.out.txt` | `01_lattice_elements/elem_AE_exact_geometry.out.txt` | 出力 |  |
| `v2_next/P3_b108_exact.py` | `01_lattice_elements/elem_AE_exact_geometry.py` | スクリプト |  |
| `v2_next/P3_maxwell_ramp.out.txt` | `01_lattice_elements/elem_BEND_steer_F1_maxwell_ramp.out.txt` | 出力 |  |
| `v2_next/P3_maxwell_ramp.py` | `01_lattice_elements/elem_BEND_steer_F1_maxwell_ramp.py` | スクリプト |  |
| `v2_next/P3_mult_ae_emit.gs.txt` | `01_lattice_elements/elem_MULT_AE_emittance.gs.txt` | 出力 |  |
| `v2_next/P3_mult_ae_emit.pysad.txt` | `01_lattice_elements/elem_MULT_AE_emittance.pysad.txt` | 出力 |  |
| `v2_next/P3_mult_ae_emit.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P3_mult_ae_emit.sad` | `01_lattice_elements/elem_MULT_AE_emittance.sad` | スクリプト |  |
| `v2_next/P3_mult_f1_ae_track.gs.txt` | `01_lattice_elements/elem_MULT_F1_AE_track.gs.txt` | 出力 |  |
| `v2_next/P3_mult_f1_ae_track.pysad.txt` | `01_lattice_elements/elem_MULT_F1_AE_track.pysad.txt` | 出力 |  |
| `v2_next/P3_mult_f1_ae_track.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P3_mult_f1_ae_track.sad` | `01_lattice_elements/elem_MULT_F1_AE_track.sad` | スクリプト |  |
| `v2_next/P3_ring_a.out.txt` | `01_lattice_elements/elem_BEND_steer_ring_a.out.txt` | 出力 |  |
| `v2_next/P3_ring_a.py` | `01_lattice_elements/elem_BEND_steer_ring_a.py` | スクリプト |  |
| `v2_next/P3_ring_a_ra1.gs.txt` | `01_lattice_elements/elem_BEND_steer_ring_a_ra1.gs.txt` | 出力 |  |
| `v2_next/P3_ring_a_ra2.gs.txt` | `01_lattice_elements/elem_BEND_steer_ring_a_ra2.gs.txt` | 出力 |  |
| `v2_next/P3_ring_a_ra3.gs.txt` | `01_lattice_elements/elem_BEND_steer_ring_a_ra3.gs.txt` | 出力 |  |
| `v2_next/P3_steer_ae2_tune.gs.txt` | `01_lattice_elements/elem_BEND_steer_AE2_tune.gs.txt` | 出力 |  |
| `v2_next/P3_steer_ae2_tune.pysad.txt` | `01_lattice_elements/elem_BEND_steer_AE2_tune.pysad.txt` | 出力 |  |
| `v2_next/P3_steer_ae2_tune.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P3_steer_ae2_tune.sad` | `01_lattice_elements/elem_BEND_steer_AE2_tune.sad` | スクリプト |  |
| `v2_next/P3_steer_f1_track.gs.txt` | `01_lattice_elements/elem_BEND_steer_F1_track.gs.txt` | 出力 |  |
| `v2_next/P3_steer_f1_track.pysad.txt` | `01_lattice_elements/elem_BEND_steer_F1_track.pysad.txt` | 出力 |  |
| `v2_next/P3_steer_f1_track.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P3_steer_f1_track.sad` | `01_lattice_elements/elem_BEND_steer_F1_track.sad` | スクリプト |  |
| `v2_next/P4_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 03_matching, 05_tracking, 07_spin_solenoid_special, 16_spin_polarization |
| `v2_next/P4_dapwidth_ringA.out.txt` | `05_tracking/fn_DynamicApertureSurvey_dapwidth_ring_a.out.txt` | 出力 |  |
| `v2_next/P4_dapwidth_window.gs.txt` | `05_tracking/fn_DynamicApertureSurvey_dapwidth_window.gs.txt` | 出力 |  |
| `v2_next/P4_dapwidth_window.pysad.txt` | `05_tracking/fn_DynamicApertureSurvey_dapwidth_window.pysad.txt` | 出力 |  |
| `v2_next/P4_dapwidth_window.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_dapwidth_window.sad` | `05_tracking/fn_DynamicApertureSurvey_dapwidth_window.sad` | スクリプト |  |
| `v2_next/P4_free_noop_variable.gs.txt` | `03_matching/flag_FREE_noop_variable.gs.txt` | 出力 |  |
| `v2_next/P4_free_noop_variable.pysad.txt` | `03_matching/flag_FREE_noop_variable.pysad.txt` | 出力 |  |
| `v2_next/P4_free_noop_variable.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_free_noop_variable.sad` | `03_matching/flag_FREE_noop_variable.sad` | スクリプト |  |
| `v2_next/P4_pspac_astrid2.gs.txt` | `05_tracking/flag_PSPAC_astrid2.gs.txt` | 出力 |  |
| `v2_next/P4_pspac_astrid2.pysad.txt` | `05_tracking/flag_PSPAC_astrid2.pysad.txt` | 出力 |  |
| `v2_next/P4_pspac_astrid2.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_pspac_astrid2.sad` | `05_tracking/flag_PSPAC_astrid2.sad` | スクリプト |  |
| `v2_next/P4_pspac_k64_shape.gs.txt` | `05_tracking/flag_PSPAC_k64_shape.gs.txt` | 出力 |  |
| `v2_next/P4_pspac_k64_shape.pysad.txt` | `05_tracking/flag_PSPAC_k64_shape.pysad.txt` | 出力 |  |
| `v2_next/P4_pspac_k64_shape.pysad_k64.txt` | `05_tracking/flag_PSPAC_k64_shape.pysad_k64.txt` | 出力 |  |
| `v2_next/P4_pspac_k64_shape.sad` | `05_tracking/flag_PSPAC_k64_shape.sad` | スクリプト |  |
| `v2_next/P4_sol_photons.gs.txt` | `07_spin_solenoid_special/flag_PHOTONS_sol.gs.txt` | 出力 |  |
| `v2_next/P4_sol_photons.pysad.txt` | `07_spin_solenoid_special/flag_PHOTONS_sol.pysad.txt` | 出力 |  |
| `v2_next/P4_sol_photons.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_sol_photons.sad` | `07_spin_solenoid_special/flag_PHOTONS_sol.sad` | スクリプト |  |
| `v2_next/P4_sol_spin_quad_mult.out.txt` | `16_spin_polarization/spin_sol_quad_mult.out.txt` | 出力 |  |
| `v2_next/P4_sol_spin_quad_mult.out_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_sol_spin_quad_mult.py` | `16_spin_polarization/spin_sol_quad_mult.py` | スクリプト |  |
| `v2_next/P4_wspac_tracking_beam.gs.txt` | `05_tracking/flag_WSPAC_tracking_beam.gs.txt` | 出力 |  |
| `v2_next/P4_wspac_tracking_beam.pysad.txt` | `05_tracking/flag_WSPAC_tracking_beam.pysad.txt` | 出力 |  |
| `v2_next/P4_wspac_tracking_beam.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P4_wspac_tracking_beam.sad` | `05_tracking/flag_WSPAC_tracking_beam.sad` | スクリプト |  |
| `v2_next/P5_README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 02_optics_cod, 03_matching, 04_emittance_radiation, 10_math_data, 12_gui_tkinter, 13_packages |
| `v2_next/P5_detr_astrid2.gs.txt` | `02_optics_cod/flag_DETR_astrid2.gs.txt` | 出力 |  |
| `v2_next/P5_detr_astrid2.pysad.txt` | `02_optics_cod/flag_DETR_astrid2.pysad.txt` | 出力 |  |
| `v2_next/P5_detr_astrid2.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_detr_astrid2.sad` | `02_optics_cod/flag_DETR_astrid2.sad` | スクリプト |  |
| `v2_next/P5_detr_pf.gs.txt` | `02_optics_cod/flag_DETR_pf.gs.txt` | 出力 |  |
| `v2_next/P5_detr_pf.pysad.txt` | `02_optics_cod/flag_DETR_pf.pysad.txt` | 出力 |  |
| `v2_next/P5_detr_pf.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_detr_pf.sad` | `02_optics_cod/flag_DETR_pf.sad` | スクリプト |  |
| `v2_next/P5_detr_ring_a.gs.txt` | `02_optics_cod/flag_DETR_ring_a.gs.txt` | 出力 |  |
| `v2_next/P5_detr_ring_a.pysad.txt` | `02_optics_cod/flag_DETR_ring_a.pysad.txt` | 出力 |  |
| `v2_next/P5_detr_ring_a.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_detr_ring_a.sad` | `02_optics_cod/flag_DETR_ring_a.sad` | スクリプト |  |
| `v2_next/P5_fitp_beta.gs.txt` | `03_matching/flag_FITP_beta.gs.txt` | 出力 |  |
| `v2_next/P5_fitp_beta.pysad.txt` | `03_matching/flag_FITP_beta.pysad.txt` | 出力 |  |
| `v2_next/P5_fitp_beta.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_fitp_beta.sad` | `03_matching/flag_FITP_beta.sad` | スクリプト |  |
| `v2_next/P5_fitp_five.gs.txt` | `03_matching/flag_FITP_five.gs.txt` | 出力 |  |
| `v2_next/P5_fitp_five.pysad.txt` | `03_matching/flag_FITP_five.pysad.txt` | 出力 |  |
| `v2_next/P5_fitp_five.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_fitp_five.sad` | `03_matching/flag_FITP_five.sad` | スクリプト |  |
| `v2_next/P5_fitp_one.gs.txt` | `03_matching/flag_FITP_one.gs.txt` | 出力 |  |
| `v2_next/P5_fitp_one.pysad.txt` | `03_matching/flag_FITP_one.pysad.txt` | 出力 |  |
| `v2_next/P5_fitp_one.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_fitp_one.sad` | `03_matching/flag_FITP_one.sad` | スクリプト |  |
| `v2_next/P5_measure_astrid2.gs.txt` | `02_optics_cod/fn_MEASURE_astrid2.gs.txt` | 出力 |  |
| `v2_next/P5_measure_astrid2.pysad.txt` | `02_optics_cod/fn_MEASURE_astrid2.pysad.txt` | 出力 |  |
| `v2_next/P5_measure_astrid2.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_measure_astrid2.sad` | `02_optics_cod/fn_MEASURE_astrid2.sad` | スクリプト |  |
| `v2_next/P5_measure_default_np.gs.txt` | `02_optics_cod/fn_MEASURE_default_np.gs.txt` | 出力 |  |
| `v2_next/P5_measure_default_np.pysad.txt` | `02_optics_cod/fn_MEASURE_default_np.pysad.txt` | 出力 |  |
| `v2_next/P5_measure_default_np.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_measure_default_np.sad` | `02_optics_cod/fn_MEASURE_default_np.sad` | スクリプト |  |
| `v2_next/P5_measure_fodo.gs.txt` | `02_optics_cod/fn_MEASURE_fodo.gs.txt` | 出力 |  |
| `v2_next/P5_measure_fodo.pysad.txt` | `02_optics_cod/fn_MEASURE_fodo.pysad.txt` | 出力 |  |
| `v2_next/P5_measure_fodo.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_measure_fodo.sad` | `02_optics_cod/fn_MEASURE_fodo.sad` | スクリプト |  |
| `v2_next/P5_measure_pf.gs.txt` | `02_optics_cod/fn_MEASURE_pf.gs.txt` | 出力 |  |
| `v2_next/P5_measure_pf.pysad.txt` | `02_optics_cod/fn_MEASURE_pf.pysad.txt` | 出力 |  |
| `v2_next/P5_measure_pf.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_measure_pf.sad` | `02_optics_cod/fn_MEASURE_pf.sad` | スクリプト |  |
| `v2_next/P5_measure_ring_a.gs.txt` | `02_optics_cod/fn_MEASURE_ring_a.gs.txt` | 出力 |  |
| `v2_next/P5_measure_ring_a.pysad.txt` | `02_optics_cod/fn_MEASURE_ring_a.pysad.txt` | 出力 |  |
| `v2_next/P5_measure_ring_a.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_measure_ring_a.sad` | `02_optics_cod/fn_MEASURE_ring_a.sad` | スクリプト |  |
| `v2_next/P5_micado_bump.gs.txt` | `13_packages/fn_CorrectOrbit_Bump.gs.txt` | 出力 |  |
| `v2_next/P5_micado_bump.pysad.txt` | `13_packages/fn_CorrectOrbit_Bump.pysad.txt` | 出力 |  |
| `v2_next/P5_micado_bump.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_micado_bump.sad` | `13_packages/fn_CorrectOrbit_Bump.sad` | スクリプト |  |
| `v2_next/P5_micado_variants.gs.txt` | `13_packages/fn_CorrectOrbit_Micado_variants.gs.txt` | 出力 |  |
| `v2_next/P5_micado_variants.pysad.txt` | `13_packages/fn_CorrectOrbit_Micado_variants.pysad.txt` | 出力 |  |
| `v2_next/P5_micado_variants.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_micado_variants.sad` | `13_packages/fn_CorrectOrbit_Micado_variants.sad` | スクリプト |  |
| `v2_next/P5_oy_pipeline.gs.txt` | `13_packages/fn_OYStability_pipeline.gs.txt` | 出力 |  |
| `v2_next/P5_oy_pipeline.pysad.txt` | `13_packages/fn_OYStability_pipeline.pysad.txt` | 出力 |  |
| `v2_next/P5_oy_pipeline.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_oy_pipeline.sad` | `13_packages/fn_OYStability_pipeline.sad` | スクリプト |  |
| `v2_next/P5_oy_pipeline_unstable.gs.txt` | `13_packages/fn_OYStability_pipeline_unstable.gs.txt` | 出力 |  |
| `v2_next/P5_oy_pipeline_unstable.pysad.txt` | `13_packages/fn_OYStability_pipeline_unstable.pysad.txt` | 出力 |  |
| `v2_next/P5_oy_pipeline_unstable.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_oy_pipeline_unstable.sad` | `13_packages/fn_OYStability_pipeline_unstable.sad` | スクリプト |  |
| `v2_next/P5_proton_beta_p1p5.gs.txt` | `04_emittance_radiation/var_MASS_proton_beta_p1p5.gs.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p1p5.pysad.txt` | `04_emittance_radiation/var_MASS_proton_beta_p1p5.pysad.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p1p5.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_proton_beta_p1p5.sad` | `04_emittance_radiation/var_MASS_proton_beta_p1p5.sad` | スクリプト |  |
| `v2_next/P5_proton_beta_p3.gs.txt` | `04_emittance_radiation/var_MASS_proton_beta_p3.gs.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p3.pysad.txt` | `04_emittance_radiation/var_MASS_proton_beta_p3.pysad.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p3.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_proton_beta_p3.sad` | `04_emittance_radiation/var_MASS_proton_beta_p3.sad` | スクリプト |  |
| `v2_next/P5_proton_beta_p30.gs.txt` | `04_emittance_radiation/var_MASS_proton_beta_p30.gs.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p30.pysad.txt` | `04_emittance_radiation/var_MASS_proton_beta_p30.pysad.txt` | 出力 |  |
| `v2_next/P5_proton_beta_p30.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_proton_beta_p30.sad` | `04_emittance_radiation/var_MASS_proton_beta_p30.sad` | スクリプト |  |
| `v2_next/P5_specfn_cx_scan.gs.txt` | `10_math_data/fn_SpecialFunctions_complex_scan.gs.txt` | 出力 |  |
| `v2_next/P5_specfn_cx_scan.pysad.txt` | `10_math_data/fn_SpecialFunctions_complex_scan.pysad.txt` | 出力 |  |
| `v2_next/P5_specfn_cx_scan.pysad_k64.txt` | — | 出力 | 削除 (現在は存在しない) |
| `v2_next/P5_specfn_cx_scan.sad` | `10_math_data/fn_SpecialFunctions_complex_scan.sad` | スクリプト |  |
| `v2_next/P5_tk90_canvas.png` | `12_gui_tkinter/Tk90_canvas.png` | 図 |  |
| `v2_next/P5_tk90_entry_listbox_scale.png` | `12_gui_tkinter/Tk90_entry_listbox_scale.png` | 図 |  |
| `v2_next/P5_tk90_kbmainframe.png` | `12_gui_tkinter/Tk90_kbmainframe.png` | 図 |  |
| `v2_next/P5_tk90_shot.py` | `12_gui_tkinter/Tk90_shot.py` | スクリプト |  |
| `v2_next/P5_tk90_window_button.png` | `12_gui_tkinter/Tk90_window_button.png` | 図 |  |
| `v2_next/README.md` | — | md | 分割: 各フォルダの README・各例の NAME.md へ。索引・説明 md。新フォルダの README と各例の NAME.md に再編。内容の移動先: 01_lattice_elements, 02_optics_cod, 03_matching, 04_emittance_radiation, 05_tracking, 07_spin_solenoid_special, 08_scripting, 09_io_system, 10_math_data, 12_gui_tkinter, 13_packages, 16_spin_polarization |
