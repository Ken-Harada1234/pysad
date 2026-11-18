# example/v2_next: 担当 B1 の例 (v2next、2026-11-18)

記録は `docs/history/v2next_B1.md`。各台本の先頭に、条件と何を確かめる例かを書いた。実行は `py -3.12 run_sad.py B1_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)。
Python の例は `py -3.12 B1_xxx.py` (出力は `.out.txt`)。gs の出力は `local/wake/gsrun.sh` で取った (`.gs.txt`)。

| 例 | 項目 | pysad と gs の関係 |
|---|---|---|
| `B1_bend_f1_emittance.sad` | 1: Emittance[]・CALC6D の BEND の F1・FRINGE | k64 は gs と全桁一致 (`.pysad_k64.txt` と `.gs.txt`)。v2 は端の角の補正で U0 が 1e-3 違う |
| `B1_k1_bend_edge.py` | 2: K1≠0 の BEND の CALC と追跡 (どちらが正しいか) | 独立な厳密な写像 (ローレンツ力の積分) との差の表。誤りは K64 の追跡の薄い端 |
| `B1_wspac_emittance.sad` | 3: WSPAC の Emittance[]・CALC6D・CALC | Emittance[] のシフトは gs より 4% (x)・3% (y) 弱い (原因未特定、一次のラスレットの式とは 0.15%)。gs の CALC・CALC6D は NaN |
| `B1_sol_calc6d_intra.sad` | 4: SOL 領域の CALC6D・RADCOD つき CALC・INTRA | k64 は gs と 1e-9 (チューン・β)、INTRA の ε は 1e-7 で一致 |
| `B1_sol_spin_bigbz.py` | 6: SOL 領域の DRIFT のスピン (bz = 4.5・9 1/m)・POL の SpinTune・加速 + RAD | Thomas-BMT・解析式との差の表 |

(項目 5: 4 次元の CALC は RADTAPER を見ない。gs も同じ。台本は不要 (`test_calc_4d_ignores_radtaper_like_k64`)。)
