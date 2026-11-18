# example/v2_next: 担当 C の例 (v2next、2026-11-18)

記録は `docs/history/v2next_C.md`。各台本の先頭に、条件と何を確かめる例かを書いた。実行は `py -3.12 run_sad.py C3_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)。
gs の出力は `example/tools/run_example.py` で取った (`.gs.txt`)。`.pysad_k64.txt` は k64 モードの出力。

| 例 | 項目 | pysad と gs の関係 |
|---|---|---|
| `C3_radint_k1_bend_vdisp.sad` | 5 (B107): K1≠0 の BEND + 結合の Y 方向の分散の RADINT | k64 は gs と全桁一致 (I4U + I4V = 0.420949)。v2 は独立な I4 (0.414392) と 1e-7 |
| `C3_radint_quad_offset.sad` | 5 (B107b): 据え付けのずれた四極の RADINT | k64 は gs と全桁一致 (ΔJx = +0.0028)。v2 は独立な積分・追跡・Emittance[] の ΔJx = +0.0078 |
| `C3_fluc_nparallel.sad` | 7: FLUC の乱数列 (NPARA=1 と 4) | k64 は gs と rms dp・最終の種まで全桁一致。v2 は 6e-6 (1 周)〜1e-3 (100 周) 違う (A1・A2 の追跡の修正) |
| `C3_radfield_bend.sad` | 2: RadiationField の振幅 (Liénard-Wiechert のピーク) | v2 の RADLIGHT の軌道 513 点で比 1.004。k64 は軌道が 1 点で NaN、gs は SIGSEGV (B49) |
| `C3_trcod_orbit.sad` | 1: TRCOD の閉軌道の測定 | CALC との差 7e-10、解析式との差 5e-6。gs は SIGSEGV (B42) |
| `C3_synchrobeta_coupled.sad` | 3: 結合のあるリングの SynchroBetaEmittance | gs も走り、k64 は 7 桁一致。v2 は 1.3% 違う |
