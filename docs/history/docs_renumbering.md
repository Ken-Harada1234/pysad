# docs の旧ファイル名から新ファイル名への対応 (旧 51 本 → 新 36 本)

この文書は、`docs/` の整理 (番号の付け直し・統合・分割) の前後の対応表である。`docs/history/` と `CHANGELOG.md` は、履歴なので旧い番号・名前のまま残してある。そこに出てくる旧い名前 (`docs/38_wake_guide.md` など) は、この表で新しい文書を引く。

元の表は、整理の設計 (`docs_map.csv`) に基づく。節の対応は、`docs/history/*_merge_check.md` の照合記録からまとめた概略である (細かい対応は各 merge_check を見る)。

## 1. 旧番号から新番号への早見表

| 旧ファイル | 新ファイル | 処置 |
|---|---|---|
| `01_overview.md` | `01_overview.md` | そのまま |
| `02_lattice_main.md` | `02_lattice_main.md` | そのまま |
| `03_ffs_commands.md` | `03_ffs_commands.md` | そのまま |
| `04_sadscript.md` | `05_sadscript.md` | 改名 |
| `05_calc_optics.md` | `07_calc_optics.md` | 改名 |
| `06_matching.md` | `08_matching.md` | 改名 |
| `07_emittance.md` | `09_emittance.md` | 改名 |
| `08_tracking.md` | `10_tracking.md` | 改名 |
| `09_dynamic_aperture.md` | `11_dynamic_aperture.md` | 改名 |
| `10_benchmarks.md` | `27_benchmarks.md` | 改名 |
| `11_sad_quirks.md` | `29_sad_quirks.md` | 改名 |
| `12_pysad_api.md` | `12_pysad_api.md` | そのまま |
| `13_sad_help_index.md` | `34_sad_help_index.md` | 改名 |
| `14_pysad_vs_sadhelp.md` | `35_pysad_vs_sadhelp.md` | 改名 |
| `15_sadhelp_purpose_policy.md` | `36_sadhelp_purpose_policy.md` | 改名 |
| `16_fortran_status.md` | `33_fortran_status.md` | 統合 |
| `17_implementation_summary.md` | `01_overview.md` | 統合 |
| `18_physics_guide.md` | `23_physics_guide.md` | 改名 |
| `19_source_analysis.md` | `32_source_analysis.md` | 分割 |
| `20_translation_candidates.md` | `33_fortran_status.md` | 統合 |
| `21_sad_bug_report.md` | `30_sad_bug_report.md` | 改名 |
| `22_gui_tclbridge.md` | `20_tcltk_gui.md` | 統合 |
| `23_wake.md` | `13_wake.md` | 統合 |
| `24_beambeam.md` | `14_beambeam.md` | 統合 |
| `25_coord.md` | `19_other_elements.md` | 統合 |
| `26_tcav.md` | `19_other_elements.md` | 統合 |
| `27_und.md` | `19_other_elements.md` | 統合 |
| `28_mult_angle.md` | `19_other_elements.md` | 統合 |
| `29_sol_optics.md` | `17_sol.md` | 統合 |
| `30_map.md` | `19_other_elements.md` | 統合 |
| `31_delcor.md` | `04_ffs_small_commands.md` | 統合 |
| `32_gaussiancoulomb_bbbrem.md` | `06_sadscript_library.md` | 統合 |
| `33_specfn_fit_misc.md` | `06_sadscript_library.md` | 統合 |
| `34_sad_packages.md` | `06_sadscript_library.md` | 統合 |
| `35_ffs_small_commands.md` | `04_ffs_small_commands.md` | 改名 |
| `36_spin_csr.md` | `16_spin_pol.md` | 分割 |
| `37_plot_canvas.md` | `21_graphics.md` | 統合 |
| `38_wake_guide.md` | `13_wake.md` | 統合 |
| `39_beambeam_guide.md` | `14_beambeam.md` | 統合 |
| `40_csr_guide.md` | `15_csr.md` | 統合 |
| `41_tcltk_guide.md` | `20_tcltk_gui.md` | 統合 |
| `42_graphics_guide.md` | `21_graphics.md` | 統合 |
| `43_pol_guide.md` | `16_spin_pol.md` | 統合 |
| `44_sol_guide.md` | `17_sol.md` | 統合 |
| `45_radlight_guide.md` | `18_radlight.md` | 統合 |
| `46_unimplemented_and_unverified.md` | `31_unimplemented_and_unverified.md` | 改名 |
| `47_v2_vs_k64.md` | `28_v2_vs_k64.md` | 改名 |
| `48_frequency_map_lma.md` | `26_frequency_map_lma.md` | 改名 |
| `49_resonance_terms_adts.md` | `25_resonance_terms_adts.md` | 改名 |
| `50_adts_rdt_primer.md` | `24_adts_rdt_primer.md` | 改名 |
| `51_lattice_dxf.md` | `22_lattice_dxf.md` | 改名 |
| `sadhelp_summary.html` | `sadhelp_summary.html` | そのまま |

## 2. 新ファイルごとの統合元 (旧ファイル)

| 新ファイル | 旧ファイル (土台は先頭) |
|---|---|
| `01_overview.md` | `01_overview.md`、`17_implementation_summary.md` |
| `02_lattice_main.md` | `02_lattice_main.md` |
| `03_ffs_commands.md` | `03_ffs_commands.md` |
| `04_ffs_small_commands.md` | `31_delcor.md`、`35_ffs_small_commands.md` |
| `05_sadscript.md` | `04_sadscript.md` |
| `06_sadscript_library.md` | `32_gaussiancoulomb_bbbrem.md`、`33_specfn_fit_misc.md`、`34_sad_packages.md` |
| `07_calc_optics.md` | `05_calc_optics.md` |
| `08_matching.md` | `06_matching.md` |
| `09_emittance.md` | `07_emittance.md` |
| `10_tracking.md` | `08_tracking.md` |
| `11_dynamic_aperture.md` | `09_dynamic_aperture.md` |
| `12_pysad_api.md` | `12_pysad_api.md` |
| `13_wake.md` | `23_wake.md`、`38_wake_guide.md` |
| `14_beambeam.md` | `24_beambeam.md`、`39_beambeam_guide.md` |
| `15_csr.md` | `40_csr_guide.md` |
| `16_spin_pol.md` | `36_spin_csr.md`、`43_pol_guide.md` |
| `17_sol.md` | `29_sol_optics.md`、`44_sol_guide.md` |
| `18_radlight.md` | `45_radlight_guide.md` |
| `19_other_elements.md` | `25_coord.md`、`26_tcav.md`、`27_und.md`、`28_mult_angle.md`、`30_map.md` |
| `20_tcltk_gui.md` | `22_gui_tclbridge.md`、`41_tcltk_guide.md` |
| `21_graphics.md` | `37_plot_canvas.md`、`42_graphics_guide.md` |
| `22_lattice_dxf.md` | `51_lattice_dxf.md` |
| `23_physics_guide.md` | `18_physics_guide.md` |
| `24_adts_rdt_primer.md` | `50_adts_rdt_primer.md` |
| `25_resonance_terms_adts.md` | `49_resonance_terms_adts.md` |
| `26_frequency_map_lma.md` | `48_frequency_map_lma.md` |
| `27_benchmarks.md` | `10_benchmarks.md` |
| `28_v2_vs_k64.md` | `47_v2_vs_k64.md` |
| `29_sad_quirks.md` | `11_sad_quirks.md` |
| `30_sad_bug_report.md` | `21_sad_bug_report.md` |
| `31_unimplemented_and_unverified.md` | `46_unimplemented_and_unverified.md` |
| `32_source_analysis.md` | `19_source_analysis.md` |
| `33_fortran_status.md` | `16_fortran_status.md`、`20_translation_candidates.md` |
| `34_sad_help_index.md` | `13_sad_help_index.md` |
| `35_pysad_vs_sadhelp.md` | `14_pysad_vs_sadhelp.md` |
| `36_sadhelp_purpose_policy.md` | `15_sadhelp_purpose_policy.md` |
| `sadhelp_summary.html` | `sadhelp_summary.html` |

## 3. 節の対応 (概略)

統合・分割していない文書は章番号を保った。節番号が変わったものと、統合先の節を、旧ファイルごとに示す。旧い文書の `§` 付き参照は、この表で新しい文書の節に読み替える。

| 旧ファイル | 新ファイルと節の対応 |
|---|---|
| `01_overview.md` → `01_overview.md` | 章番号は同じ (旧 01 の §4 (古い目次) と §6 は 01 の新しい目次・31 への案内に置き換えた) |
| `02_lattice_main.md` → `02_lattice_main.md` | 章番号は同じ |
| `03_ffs_commands.md` → `03_ffs_commands.md` | 章番号は同じ |
| `04_sadscript.md` → `05_sadscript.md` | 章番号は同じ。§1.1 (動かし方) を追加し、旧 §7 の使い方の例をそこへ。§6 は 6.1〜6.3 に分けた |
| `05_calc_optics.md` → `07_calc_optics.md` | §1〜§6 は同じ。§0 (使い方と結果の読み方) を追加。旧 §5 (ベンチの要約) は §5 検証 (3 行) に |
| `06_matching.md` → `08_matching.md` | §1〜§7 は同じ。§0 を追加。旧 §7 (照合結果) は §7 検証 |
| `07_emittance.md` → `09_emittance.md` | §1〜§6 は同じ。§0 を追加。旧 §6 (照合の要約) は §6 検証 |
| `08_tracking.md` → `10_tracking.md` | 章番号は同じ (§1〜§9、9.x も同じ)。旧 §7 (照合結果) は §7 検証と限界 |
| `09_dynamic_aperture.md` → `11_dynamic_aperture.md` | 章番号は同じ。旧 §6 (照合結果) は §6 検証と見積り |
| `10_benchmarks.md` → `27_benchmarks.md` | 章番号は同じ。§0 を追加。§5.1・§6 は最新の値・見出しに |
| `11_sad_quirks.md` → `29_sad_quirks.md` | 章番号は同じ (A 符号も同じ) |
| `12_pysad_api.md` → `12_pysad_api.md` | 章番号は同じ |
| `13_sad_help_index.md` → `34_sad_help_index.md` | 章番号は同じ |
| `14_pysad_vs_sadhelp.md` → `35_pysad_vs_sadhelp.md` | 章番号は同じ |
| `15_sadhelp_purpose_policy.md` → `36_sadhelp_purpose_policy.md` | 章番号は同じ |
| `16_fortran_status.md` → `33_fortran_status.md` | 旧 §1〜§8 (ファイル別の表) は新 33 の §1〜§8。旧 §9 は 31 へ。旧 §10 は 33 §10 |
| `17_implementation_summary.md` → `01_overview.md` | 旧 §1・§2・§5 は 01 の §1・§3〜§7。旧 §3 は 27、旧 §4 は 29、旧 §6〜§8 は 31 |
| `18_physics_guide.md` → `23_physics_guide.md` | 旧 §1〜§23 は新 23 の §1〜§23 (同じ番号)。旧 §24〜§30 は機能別の文書 (13〜18) の物理の章へ。新 23 §24 は案内の表 |
| `19_source_analysis.md` → `32_source_analysis.md` | 旧 §1〜§4・§6・§8 は 32 の §1〜§4・§5・§6。旧 §5・§7 の表は 33 §1〜§8。旧 §9〜§11 は 31 |
| `20_translation_candidates.md` → `33_fortran_status.md` | 旧 §1〜§4 は 33 §9。旧 §5・§6 は 31 |
| `21_sad_bug_report.md` → `30_sad_bug_report.md` | 章番号は同じ (B 番号も同じ) |
| `22_gui_tclbridge.md` → `20_tcltk_gui.md` | 旧 §1 方針 → 20 §1、§2 → §2、§3 → §4.1、§4 → §4.3、§5 (Class) → §5、§6 → §7、§7 (試験の実行) → §7.1 |
| `23_wake.md` → `13_wake.md` | 旧 §1 → 13 §1・§2・§3.1、§2 → §3.1、§2a → §3.2、§2b → §3.3、§2c → §5.1・§5.2・§2.3、§3 → §5.3、§4 → 検証 (history)、§5 → §7、§6 → §6、§7 → §3.4・§4.2、§8 → §5.4 |
| `24_beambeam.md` → `14_beambeam.md` | 旧 §1 → 14 §2.1、§2 → §3.1、§2b → §3.2、§3 → §6、§4 → §5、§5 → §7 |
| `25_coord.md` → `19_other_elements.md` | 19 §1 |
| `26_tcav.md` → `19_other_elements.md` | 19 §2 |
| `27_und.md` → `19_other_elements.md` | 19 §3 |
| `28_mult_angle.md` → `19_other_elements.md` | 19 §4 |
| `29_sol_optics.md` → `17_sol.md` | 旧 作り方・Emittance[] の対応表 → 17 §8、旧 幾何 → §4.5、旧 限界 → §7、旧 v2 の要約 → §5.1 |
| `30_map.md` → `19_other_elements.md` | 19 §5 |
| `31_delcor.md` → `04_ffs_small_commands.md` | 旧 冒頭・書式 → 04 §2、k64 と v2 の違い → §2.1、作り方 → §2.2、v2 の根拠 → §2.3、限界 → §2.4 |
| `32_gaussiancoulomb_bbbrem.md` → `06_sadscript_library.md` | GaussianCoulomb 系 → 06 §4、BBBrem1 → 14 §4.4 |
| `33_specfn_fit_misc.md` → `06_sadscript_library.md` | 旧 §1 → 06 §1、§2 → §2、§3〜§5 → §3.1〜§3.3、§6 → §6 |
| `34_sad_packages.md` → `06_sadscript_library.md` | 旧 パッケージ表・UpValue・corr.py → 06 §5、§5.1〜§5.5、限界 → §7 |
| `35_ffs_small_commands.md` → `04_ffs_small_commands.md` | 旧 表 → 04 §1、動作 → §1.1、限界・K64 との差 → §3 |
| `36_spin_csr.md` → `16_spin_pol.md` | 旧 §1 (POL) → 16 §2・§3・§7、§2 (CSR) → 15 §4・§8、§3 (スピン追跡) → 16 §4、§4 (RADLIGHT) → 18 |
| `37_plot_canvas.md` → `21_graphics.md` | 旧 §1・§2 → 21 §2・§3・§7、§3 → §10.1、§4 → §10.2、§5 → §7・§9 |
| `38_wake_guide.md` → `13_wake.md` | 旧 §1 → 13 §1、§2.1 → §2.1、§2.2 → §2.2、§3.1 → §3.1、§3.2 → §3.4、§4 → §4.1、§5 → §2.3、§5a → §5.1、§6 → §5.3・§5.4、§7 → §8、§8 → §6、§9 → §7、§10 → §9 |
| `39_beambeam_guide.md` → `14_beambeam.md` | 旧 §0 → 14 §1、§1 → §2.1、§2 → §3.1、§3 → §3.2、§4 → §4.1、§5 → §4.2・§2.4、§6 → §4.3、§7 → §2.2・§2.3・§2.5、§8 → §4.4、§9 → §4.5、§10 → §6、§11 → §5、§12 → §7、§13 → §9 |
| `40_csr_guide.md` → `15_csr.md` | 旧 §1 → 15 §1、§2 → §2、§3 → §3、§4 → §5、§5 → §6、§6 → §7.1、§7 → §7.2、§8 → §9・§8 |
| `41_tcltk_guide.md` → `20_tcltk_gui.md` | 旧 §1〜§4 は 20 の §1〜§4 (同じ番号)、§5 → §7.1、§6 → §6、§7 → §7、§8 → §8、§9 → §9 |
| `42_graphics_guide.md` → `21_graphics.md` | 旧 §1〜§4 は 21 の §1〜§4、§5 → §7、§6 → §8、§7 → §9、§7b → §5、§7c → §6、§8 → §11 |
| `43_pol_guide.md` → `16_spin_pol.md` | 旧 §1〜§3 は 16 の §1〜§3、§4 → §4.1、§5.1〜§5.3 → §4.2〜§4.4、§5.4 → §5.1、§6 → §5.2、§7 → §6、§8 → §7、§8a → §5.3、§9 → §9 |
| `44_sol_guide.md` → `17_sol.md` | 旧 §1 → 17 §1、§2 → §2、§3 → §3、§4 → §4、§5 → §5.2、§6 → §6、§7 → §7、§7a → §5.3、§8 → §9 |
| `45_radlight_guide.md` → `18_radlight.md` | 旧 §1 → 18 §1、§2 → §2、§3 → §5.1、§3a → §5.2、§4 → §3、§5 → §6、§6 → §4、§7 → §7、§8 → §9 |
| `46_unimplemented_and_unverified.md` → `31_unimplemented_and_unverified.md` | 章番号は同じ (§1〜§9)。§6・§7 の見出しの語句のみ変更 |
| `47_v2_vs_k64.md` → `28_v2_vs_k64.md` | 章番号は同じ (§1〜§6) |
| `48_frequency_map_lma.md` → `26_frequency_map_lma.md` | 章番号は同じ (§1〜§10) |
| `49_resonance_terms_adts.md` → `25_resonance_terms_adts.md` | 章番号は同じ (§1〜§15) |
| `50_adts_rdt_primer.md` → `24_adts_rdt_primer.md` | 章番号は同じ |
| `51_lattice_dxf.md` → `22_lattice_dxf.md` | 旧 §2b → 22 §3、§2c → §4、§2d → §5、§3 → §6、§4 → §7、§5 → §8 (§1・§2 は同じ) |

## 4. 注意

- 新 23 (`23_physics_guide.md`) の §1〜§23 は旧 18 と同じ番号。旧 18 の §24〜§30 (POL・SOL・WAKE・BEAM-BEAM・CSR・RADLIGHT) は、機能別の文書 (16・17・13・14・15・18) の「物理と式」の章にある。
- 旧 46 (未対応・未確認の一覧) は新 31。旧 16 §9、旧 17 §6〜§8、旧 19 §9〜§11、旧 20 §5〜§6 の項目は、新 31 に 1 項目ずつ照合して寄せた。
- `docs/sadhelp_summary.html` は名前を変えず、中の参照だけ新 34・35・36 に合わせた。
- 旧 `docs/html/` (生成済みの旧い HTML) は、リポジトリの外 (追跡しない `local/old_docs_html/`) に移した。HTML は md と同じフォルダに同名で作る (`tools/build_html.py`)。
