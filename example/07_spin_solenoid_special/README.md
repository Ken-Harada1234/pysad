# 07_spin_solenoid_special (スピン・ソレノイド・特殊)

使ったラティスは `00_lattice/` の `FODO_ring_V4MV.sad` (POL)、新規の `sol_ring.sad`、`multa_ring.sad`、`coord_ring.sad`、`radlight_line.sad`。
各例は 4 ファイル。一致の記号: ○ 一致、△ 差あり、× gs 失敗または pysad 未対応。

| 名前 | 1 行説明 | gs と pysad |
|---|---|---|
| [flag_POL_RADPOL](flag_POL_RADPOL.md) | 偏極時間・スピンチューン・Sokolov-Ternov の立ち上がり | ○ (k64 モード。v2 の既定ではスピン追跡の `spin angle` が gs と違う: 原因は未確認) |
| [elem_SOL](elem_SOL.md) | ソレノイド領域の 4×4 伝達行列 (結合) | ○ |
| [elem_MULT_ANGLE](elem_MULT_ANGLE.md) | ANGLE 付き MULT (BEND と同じ・複合機能型) | ○ |
| [elem_COORD](elem_COORD.md) | 座標変換 (ずらし・回転)、挟んだリングと打ち消し | ○ |
| [flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum](flag_RADLIGHT_fn_RadiationField_fn_RadiationSpectrum.md) | 軌道の記録とスペクトル (ガウスパルスと照合) | △ (軌道・スペクトルは k64 モードで ○。v2 の既定は BEND の中も記録するので点数が 13 で gs の 11 と違う: B65。**RadiationField は gs が SIGSEGV で失敗**、pysad は動く) |
| [fn_SymplecticJ](fn_SymplecticJ.md) | シンプレクティック行列 J と M^T J M = J | ○ |
| elem_SOL_tilt_offset | ずれて傾いたソレノイド (GEO=1 の SOL の `DX DY DPX DPY`)、端の放射 `F1`、領域内の四極で GZ を合わせる幾何の fit。設計軌道の粒子は 1 周して 0 に戻り、F1 による U0 の増分 2.11 eV は式 (2/3) r_e γ³ (bz r/2)²/F1 の和と合う | v2 の機能 (gs の出力なし。k64 モードでは SOL のずれ・回転は未対応) |
| elem_SOL_RF_TRPT | SOL 領域 (BZ = 1 T) の中の CAVI・RF つき MULT の TRPT の加速、AUTOPHI、wake。CALC (TRPT) の軌道は追跡と同じで、AUTOPHI で設計軌道が公称の加速を受ける (DDP 7.4e-3 → 2e-10) | v2 の機能 (gs の出力なし) |
| elem_SOL_region_elements | SOL 領域の中の BEND・MULT・CAVI・APERT。CALC (INS) の軌道と追跡が全桁一致し、APERT で落ちた要素が LOSSMAP に出る。リングの `Emittance[]` も出す | v2 の機能 (gs の出力なし。K64 は領域内の偏向をほぼ無視する。v2 は一様磁場の中の運動方程式どおりに追う) |
| elem_SOL_calc6d_intra | SOL 領域のあるリング (3 GeV、周長 45.8 m の 12 セル FODO) の CALC6D、RADCOD つきの CALC、INTRA (IBS) | k64 モードでは gs と C4D・C6D・C6DRC が 1e-9 以内、INTRA の εx が 1e-7 で一致。v2 の値との差は `Emittance[]` の式 (B75) などによる (v2 の出力は `.pysad.txt`) |
| elem_SOL_bound0_track | `BOUND=0` の SOL (領域の入口だけ、終端は `SE`) を通る追跡 (1.7 GeV、BZ = 0.8 T、QUAD ×4、3 粒子 1 周、RAD OFF) | 6 座標 18 個が gs と 10 桁一致 (v2・k64 とも)。試験 `test_v2_F4_sol_bound0.py` |
| flag_PHOTONS_bend | BEND・QUAD の光子 (`PhotonList`)。4 粒子 × 2 セル、ずれ・回転のある BEND を含む | k64 モードで gs と光子 358 個が全列一致 (E は 3e-8)。v2 は偏光 xi1 が ±1 以内 (K64 は −3735 などの誤り)。試験 `test_v2_B3_photons.py` |
| flag_PHOTONS_mult | MULT (K0・K1・K2、ずれ・回転つき)・SEXT・OCT の光子 | k64 モードで gs と一致 (試験 `test_v2_D1_emit.py`) |
| flag_PHOTONS_beambeam_bend | BEAMBEAM と BEND を含むリングの光子 | k64 モードで gs と一致 (試験 `test_v2_D1_emit.py`) |
| flag_PHOTONS_coord_tcavi | COORD・TCAVI を同じラティスに置いたときの光子 (BEND ×4 + COORD の対 ×2 + TCAVI ×2、4 粒子 8 周、RAD・FLUC) | gs と光子 278 個が個数・全 12 列 (1e-10) で一致。K64 は COORD・TCAVI の中では光子を記録しない。試験 `test_v2_F2_photons_coord_tcavi.py`。`.pysad.txt` は無い |
| flag_PHOTONS_sol | SOL 領域内の光子。3 GeV、BZ = 1 T の領域の中の DRIFT・MULT (K0 = 0.01、L = 0.4)・QUAD (K1 = 0.8、L = 0.4) を 1500 粒子が通る。要素ごとの光子数とエネルギーの合計 (MULT は古典値 285 eV/粒子) | gs は MULT・DRIFT を記録し QUAD は記録しない (B68 の疑い)。領域内の MULT は偏向を無視する (B69)。k64 モードは領域内の MULT・BEND が未対応 |
| fn_RadiationField_bend | `RadiationField` の振幅を Liénard-Wiechert の式と比べる (10.2 MeV の電子、γ = 20、ρ = 1 m の BEND)。v2 の RADLIGHT の軌道 513 点でピークの比 1.004 | gs は SIGSEGV で落ちる (B49)。`.gs.txt` はその記録。独立な全比較は試験 `test_v2_C3_radfield.py` |
| flag_RADLIGHT_orbit_field | RADLIGHT の軌道 {t, x, y, z} に横方向の座標が入り、BEND の中も記録する (RADLIGHTDIV 点)。偏向磁石 (ρ = 1 m、γ = 20) の放射場のピークを、接線上 100 m の観測点で Liénard-Wiechert の式と比べる (0.3% 以内) | v2 の機能 (B65。gs の出力なし) |

図: `flag_RADLIGHT_spectrum.png`。

## このフォルダで分かったこと
- スピンチューン $a\gamma=6.808148$、Sokolov-Ternov の偏極時間 (143.3652 s、式 143.3651 s)、`sy` の立ち上がりは式と一致し、gs と pysad も同じ。`NORADPOL` で立ち上がりが止まる。
- ソレノイドの結合 (`SC`、`S²/K` の項) は式の 4×4 行列と 8 桁一致。SOL 内の四極と RAD の組み合わせは gs が信用できない (癖) ので避けた。
- COORD の逆変換 (`DIR=1`) は、`DX` を同じ値にし、回転角の符号を反転するだけで元に戻る (`tcoord.f` の `DIR≠0` の分岐では `DY`、`DZ` の符号が反転する)。
- `RadiationField` は gs (この K64 ビルド) が SIGSEGV で落ちるため、値は比べられない。
