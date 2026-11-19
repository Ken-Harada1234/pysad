# example/v2_next: v2.1.0 の作業で作った例

v2.0.0 の後の作業 (v2.1.0) で、未対応・未確認の項目を追求したときの例です。各例は台本 (`.sad` または `.py`) と、出力 (`.gs.txt`: gs、`.pysad.txt`: pysad の v2 モード、`.k64.txt`: pysad の k64 モード、`.out.txt`: Python の出力) からなります。担当ごとに先頭の記号 (A1 など) が付いています。記録は `docs/history/v2next_<担当>.md`、試験は `pysad/tests/test_v2_<担当>_*.py` です。

gs と同じ結果が要るときは、環境変数 `PYSAD_COMPAT=k64` で走らせます (既定は v2)。

| 担当 | 内容 | 説明 | 例 (台本) |
|---|---|---|---|
| A1 | 不具合と gs との差 (DETR、MASS・CHARGE、F1・FRINGE、CAVI の DX・DY、USE と EMIT、RADTAPER、PSPAC) | [A1_README.md](A1_README.md) | [`bend_f1_adjacent`](A1_bend_f1_adjacent.sad)、[`cavity_offset`](A1_cavity_offset.sad)、[`n12_ins_detr`](A1_n12_ins_detr.sad)、[`proton_ring`](A1_proton_ring.sad)、[`pspac_resolution`](A1_pspac_resolution.sad)、[`radcod_closure`](A1_radcod_closure.sad)、[`use_after_emit`](A1_use_after_emit.sad) |
| A2 | 物理の検証 (スピン角、偏極度、RADINT、ずれた BEND) | (記録を参照) | [`dx_bend`](A2_dx_bend.sad)、[`i4v_coupled`](A2_i4v_coupled.sad)、[`pol_fluc_tracking`](A2_pol_fluc_tracking.py)、[`pol_mpw`](A2_pol_mpw.sad)、[`radint_pf`](A2_radint_pf.sad)、[`spin_turn`](A2_spin_turn.sad) |
| B1 | 光学・Emittance の実装 (F1・FRINGE、WSPAC、CALC6D と SOL、SOL 領域の POL) | [B1_README.md](B1_README.md) | [`bend_f1_emittance`](B1_bend_f1_emittance.sad)、[`k1_bend_edge`](B1_k1_bend_edge.py)、[`sol_calc6d_intra`](B1_sol_calc6d_intra.sad)、[`sol_spin_bigbz`](B1_sol_spin_bigbz.py)、[`wspac_emittance`](B1_wspac_emittance.sad) |
| B2 | スクリプト・関数 (成分ごとの FREE、複素特殊関数、Parallelize、文脈、Packages) | (記録を参照) | [`context`](B2_context.sad)、[`free_component`](B2_free_component.sad)、[`packages_compare`](B2_packages_compare.sad)、[`parallelize`](B2_parallelize.sad)、[`special_functions_complex`](B2_special_functions_complex.sad)、[`symbol_message`](B2_symbol_message.sad)、[`visit_bye`](B2_visit_bye.sad) |
| B3 | 偏極・光子・幾何 (側帯項、光子、tdgeo、NAFF) | [B3_README.md](B3_README.md) | [`geofit_quad_dx`](B3_geofit_quad_dx.sad)、[`naff_tunes`](B3_naff_tunes.py)、[`photons_bend`](B3_photons_bend.sad)、[`pol_sideband`](B3_pol_sideband.sad)、[`radpol_roundoff`](B3_radpol_roundoff.py)、[`spin_mc`](B3_spin_mc.py) |
| C3 | 独立な検証 (RADINT、RadiationField、MEASURE、FLUC の乱数列) | [C3_README.md](C3_README.md) | [`fluc_nparallel`](C3_fluc_nparallel.sad)、[`radfield_bend`](C3_radfield_bend.sad)、[`radint_k1_bend_vdisp`](C3_radint_k1_bend_vdisp.sad)、[`radint_quad_offset`](C3_radint_quad_offset.sad)、[`synchrobeta_coupled`](C3_synchrobeta_coupled.sad)、[`trcod_orbit`](C3_trcod_orbit.sad) |
| D1 | 物理の残り (NORFSW、WSPAC、SynchroBeta、POL+TRPT+加速、光子) | (記録を参照) | [`norfsw_orbit`](D1_norfsw_orbit.sad)、[`photons_bb_bend`](D1_photons_bb_bend.sad)、[`photons_mult`](D1_photons_mult.sad)、[`spin_accel`](D1_spin_accel.py) |
| D2 | パーサ・GUI・Packages (MAIN の式、Tk 9.0、Micado、CSROYMatrix) | [D2_README.md](D2_README.md) | [`correctlocal`](D2_correctlocal.sad)、[`csroymatrix`](D2_csroymatrix.sad)、[`main_expr`](D2_main_expr.sad)、[`micado`](D2_micado.sad)、[`periodicr`](D2_periodicr.sad) |
| E | 傾いたソレノイドの偏極、ANGLE=0 の BEND と MULT の F1、脱偏極率の限界 | (記録を参照) | [`b108_poles`](E_b108_poles.py)、[`mult_f1`](E_mult_f1.sad)、[`pol_kick_response`](E_pol_kick_response.py)、[`sol_tilt_spin`](E_sol_tilt_spin.py)、[`steer_f1`](E_steer_f1.sad) |

記録: [`docs/history/`](../../docs/history/) の `v2next_A1.md`・`v2next_A2.md`・`v2next_B1.md`・`v2next_B2.md`・`v2next_B3.md`・`v2next_C.md` (C3 の例)・`v2next_D1.md`・`v2next_D2.md`・`v2next_E.md`。
