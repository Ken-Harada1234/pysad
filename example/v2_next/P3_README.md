# example/v2_next: 担当 P3 の例 (v2next、2026-11-19)

記録は `docs/history/v2next_P3.md`。実行は `py -3.12 run_sad.py example/v2_next/P3_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)、`.py` は
`py -3.12 example/v2_next/P3_xxx.py` (リポジトリ直下から)。gs の出力は WSL の gs (10-19 ビルド、2026-11-19) で取った (`.gs.txt`)。`.pysad.txt` は v2、`.pysad_k64.txt` は k64 モード。

| 例 | 項目 | pysad と gs の関係 |
|---|---|---|
| `P3_steer_f1_track.sad` | 1: ステアリング (ANGLE=0 の BEND) の F1・FB1・FB2・FRINGE・AE の追跡 (K64 の tsteer) | k64 は gs と 10 桁表示の全桁一致 (NORAD・RAD)。v2 は RAD で AE の行路の補正 (B108) だけ違う |
| `P3_mult_f1_ae_track.sad` | 1・2: MULT の F1・F2・F1K1F・FB1・FB2・FRINGE・AE・ROTATE・DX の追跡 (K64 の tmulti) | k64 は gs と全桁一致 |
| `P3_mult_ae_emit.sad` | 2: MULT の AE1・AE2 の CALC・Emittance[] (K64 の tmulte) | k64 は gs と全桁一致。v2 は出口の端の蹴りの位置と行路 (B108) が違う |
| `P3_steer_ae2_tune.sad` | 2・3: 追跡の tsteer は AE2 を FB2≠0 のときだけ入れる (K64 の不整合) | k64 は gs と全桁一致 (追跡のチューン 0.1340952)。v2 は AE2 を常に入れる (0.1338835) |
| `P3_maxwell_ramp.py` (`.out.txt`) | 1: F1 の追跡を x 方向に一様な Maxwell 場の端の厳密な運動と比較 | 鉛直の収束 (fb1+fb2)/(6ρ²)、fb1=fb2 の x・px、放射の積分 (gs は走らない) |
| `P3_b108_exact.py` (`.out.txt`) | 2・3: 端の角のある極 (MULT・K1≠0・BEND) の ΔJx を厳密な幾何 + ロビンソンの定理と比較 | v2 の比は 1.03〜1.04 (Emittance[]・追跡)、k64 は約 2 倍。残る 2.7% の切り分け |
| `P3_ring_a.py` (`.out.txt`、`P3_ring_a_ra*.gs.txt`) | 1・2: 大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m) に補正磁石 2 台 | k64 の追跡・CALC・Emittance[] は gs と一致。追跡のチューンと CALC のチューンの一致 (v2) |

注意: `.sad` の数値は、試験 (`pysad/tests/test_v2_P3_f1.py`・`test_v2_P3_ae.py`) が使う gs の値 (`pysad/tests/_p3_gsdata.py`) と同じ台本から取った。
`P3_ring_a_ra*.gs.txt` は大きなリング A のラティス (`kek-ls_v03_12.lat`) の全文の反映を含む。
