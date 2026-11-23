# example/07_spin_solenoid_special の md から取り除いた履歴の記述
例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきでそのまま残す (行番号は整理前のもの)。取り除いたあとの現在の記述は、各 md を参照。
整理日: 2026-11-23。

## elem_MULT_ANGLE.md

元の 45 行目:
> - $K_1$ だけ。$K_0$ (双極の過不足)、$K_2$ 以上、スキュー成分 (`SK`)、`BZ`、`VOLT` (最初の版では pysad 未対応と書いた。v2 の SOL 領域内の MULT と領域の外の RF つき MULT の記述が `docs/47_v2_changes.md` にあるが、この例では確かめていない)、`F1`、`TAPER`、`L=0` は見ていない。

## elem_SOL.md

元の 47 行目:
> - 領域は DRIFT だけ (BZ=5 T、L=1 m、3 GeV)。この例の最初の版では、BEND・MULT・CAVI が領域内にあるリングと DX・DY・CHI のずれたソレノイドは pysad 未対応だった。v2 で対応した (`docs/47_v2_changes.md` の B69・B80 と「ずれて傾いたソレノイド」。追跡・CALC・Emittance[]・EMIT。bz = 0.3〜3 1/m で確認)。ただしこの例では使っておらず、gs との比較もしていない (gs は領域内の BEND の偏向をほぼ無視し、領域内の四極 + RAD も保存値と食い違う: B68・B69)。BOUND=0 で BZ≠0 の SOL は未確認。

元の 50 行目:
> - SOL を含む CALC6D と RADCOD つき CALC は pysad で拒否される (docs/44)。v2 での扱いの変化は未確認。

## flag_POL_RADPOL.md

元の 50 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行、要確認)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。スピンの追跡 500 周の `spin angle` が、gs と k64 モードの .8793923838 に対して v2 では 1.1515563959 になる (RADPOL・NORADPOL とも。`sy/|s|` は 1e-6 で同程度、`PolarizationTime`・`SpinTune` など `Emittance[]` の値は同じ)。原因は未特定: v2 は追跡の端の角の行路 (B75) や BEND の端の縦磁場の積分 (B58) を直したので、一周の軌道・スピンの歳差が変わった可能性があるが、この例では切り分けていない。`PYSAD_COMPAT=k64` では保存してあった旧い出力と同じで gs と一致。以下の「○ 一致」は k64 モードの結果。

## flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum.md

元の 41 行目:
> - 軌道の記録 (`RADLIGHT`)、`RadiationSpectrum`: ○ 一致 (軌道の点列、スペクトルの $k$、$c$ が同じ)。**ただし v2 (既定) では軌道の点列が gs と違う** (2026-11-17 の再実行): v2 は BEND の中も円弧で進めて記録する (B65。1/γ の角度あたり `RADLIGHTDIV` 点、既定 16) ので、軌道の点が 11 点から 13 点になり、`t(end)*c` が 5.27 m から 7.67 m、`RadiationField` の戻りの寸法が {10,11} から {10,13} になる。K64 は BEND の中を記録しない。`PYSAD_COMPAT=k64` では保存してあった旧い出力 (11 点) と同じで gs と一致。`.pysad.txt` は v2 の出力。以下の「11 点」は k64 の値。

## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### elem_SOL_calc6d_intra.sad (旧 v2_next/B1_sol_calc6d_intra.sad)
> ! v2next B1 項目 4: SOL 領域のあるリングの CALC6D、RADCOD つきの CALC、INTRA (IBS)。以前は SADError・NotImplementedError。
### fn_RadiationField_bend.sad (旧 v2_next/C3_radfield_bend.sad)
> ! v2next C3 項目 2: RadiationField を Liénard-Wiechert の式と比べる …
### flag_PHOTONS_sol.sad (旧 v2_next/P4_sol_photons.sad)
> ! v2next P4 項目 2: SOL 領域内の光子 (PhotonList)。…
### flag_RADLIGHT_orbit_field.sad (旧 v2_C/C5_radlight.sad)
> ! v2 (B65): RADLIGHT の軌道 {t, x, y, z} に横方向の座標が入り、BEND の中も記録する …
作業符号 (「v2next B1 項目 4」など) と、以前の状態 (SADError・NotImplementedError) の記述を取り除いた。
