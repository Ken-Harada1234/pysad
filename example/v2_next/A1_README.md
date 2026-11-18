# example/v2_next: 担当 A1 の例 (v2next、2026-11-17)

記録は `docs/history/v2next_A1.md`。各台本の先頭に、条件 (エネルギー、周長、要素) と何を確かめる例かを書いた。実行は `py -3.12 run_sad.py A1_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)。
gs の出力は `local/wake/gsrun.sh` で取った (`.gs.txt`)。

| 台本 | 項目 | pysad v2 (`.pysad.txt`) と gs (`.gs.txt`) の関係 |
|---|---|---|
| `A1_n12_ins_detr.sad` | N12: INS の CALC の後の `Twiss["DETR"]` | 全桁一致 (DETR = 0.0021332) |
| `A1_proton_ring.sad` | TrackParticles・Emittance[] の `MASS`・`CHARGE` (30 GeV/c の陽子) | Qs、1 周・100 周の追跡が一致 (表示 6 桁)。**U0 は v2 が 6.89e-4 eV (古典式)、gs・k64 が 1.265 eV** (K64 の放射は電子の古典半径を使う) |
| `A1_bend_f1_adjacent.sad` | BEND の F1・FRINGE + RAD (隣り合う BEND) | v2 は gs と 1e-4 以内 (v2 は端の角の補正が入る。k64 は 1e-6 で一致: `.pysad_k64.txt`) |
| `A1_cavity_offset.sad` | CAVI の DX・DY・ROTATE・V1 V20 V11 V02 | v2 は gs と 1e-4 以内。k64 は全桁一致 |
| `A1_use_after_emit.sad` | EMIT のあとの USE (trf0 が残る) | gs と全桁一致 (A の値は v2 が 1e-4 違う: 端の角の補正) |
| `A1_radcod_closure.sad` | RADCOD・RADTAPER の閉軌道が追跡で閉じるか | **gs は閉じない (x 4e-7 m、RADTAPER あり 2e-6 m)、v2 は閉じる (1e-9 m)**。k64 は gs と同じ |
| `A1_pspac_resolution.sad` | PSPAC のメッシュの分解能 | gs の 2 行目以降は無効 (K64 の static 変数)。台本の先頭に別セッションで取った gs の値と解析式を書いた |

`.pysad_k64.txt` は `PYSAD_COMPAT=k64` の出力。PSPAC は k64 では追跡に入らない。
