# example/v2_next: 担当 P5 の例 (v2next、2026-11-19)

記録は `docs/history/v2next_P5.md`。実行は `py -3.12 run_sad.py example/v2_next/P5_xxx.sad` (既定は v2。`PYSAD_COMPAT=k64` で K64 の動作)。ラティスを連結する例 (`P5_measure_{ring_a,pf,astrid2}`、
`P5_detr_{ring_a,pf,astrid2}`) は台本の後半だけで、先頭に `pysad/tests/kek-ls_v03_12.lat` (大きなリング A、3 GeV、603 m)・`example/pf/lattice/v9_16_asym.sad` (PF、2.5 GeV、187 m)・
`example/astrid/lattice/ASTRID2_v01_01.sad` (ASTRID2、0.58 GeV) を連結して走らせる (`cat lattice tail > in.sad`)。gs の出力 (`.gs.txt`) は WSL の gs (10-19 ビルド、2026-11-19、
`timeout` 120〜900 s) で、走らせた後に `ps` で sad1.exe の残りが無いことを確かめた。`.pysad.txt` は v2、`.pysad_k64.txt` は k64。gs の画面にだけ台本のエコー行が出る。

| 例 | 項目 | pysad と gs の関係 |
|---|---|---|
| `P5_measure_fodo.sad` (他にラティス 3 つ) | 1: MEASURE・LOGMEASURE・TRCOD | k64 は gs の統計ブロック (4 桁) と全一致。gs は先頭 MARK の位置だけ出力する (ほかは無出力、TRCOD は SIGSEGV)。v2 は種を戻す・py の重心が 0 |
| `P5_detr_{astrid2,pf,ring_a}.sad` | 2: DETR >= 1 の CELL CALC・INS・Emittance[] | β・DETR・Emittance[] は k64・v2 とも gs と一致。一周の位相の整数部は k64 が ±1 ずれ (K64 の丸め誤差の符号)、v2 は gs と一致 |
| `P5_specfn_cx_scan.sad` | 3: 複素の特殊関数 732 式 | 追加後 641 式が gs と 1e-9。残りは gs の誤り・不正確 (scipy と比べて pysad が正しい) |
| `P5_oy_pipeline.sad`・`P5_oy_pipeline_unstable.sad` | 4: OYStability の後半 (CSROYMatrix・Eigensystem・成長率) | gs と 1e-12 (安定 −damp、不安定 0.0272796226)。通しの OYStability・OYScan は gs が動かず未確認 |
| `P5_micado_variants.sad` | 5: Micado の水平・2 平面・SVDmode->1 | k64 は gs と 5e-6 (B129 の K64 の誤りまで再現)。v2 は n で単調に下がる |
| `P5_micado_bump.sad` | 5: MakeBump と CorrectOrbit の Bump-> | MakeBump は gs と 1e-12。Bump-> は K64 の誤りで gs も pysad もエラー |
| `P5_fitp_{five,beta,one}.sad` | 6: FITP>1 (5 点・3 点 + β・FITP 1) | 最終の残差が gs と 2e-4 以内 (反復の経路は違う) |
| `P5_tk90_shot.py`・`P5_tk90_*.png` | 7: Tk 9.0 (Python 3.14) の見た目 | 窓の大きさは Tk 8.6 と同じ。画像は 3.14 + Tk 9.0 の PrintWindow |
| `P5_proton_beta_{p1p5,p3,p30}.sad` | 8: 陽子の放射の β の因子 (B112) | k64 は gs と全桁。v2 の U0 は古典式の 1/β³ 倍 (未修正の限界) |

項目 9 (文法の細部) は gs が最初の構文エラーで止まるため 1 式ずつ別の台本で比べた。試験 `test_v2_P5_syntax.py` の表が gs の実測 (`BEND B1=(L=1 ANGLE=式)` のあと `Element["ANGLE","B1"]`):
`E=0.5;` の `ANGLE=E` → .5、`Pi=3;` の `Pi/10` → .3、`Degree=0.1;` → .1、`A=2;` の `A*0.1` → .2、`1+` → 1、`2*3+` → 1、`45DEG`・`SIN[0.1]`・`pi/10`・`ZZZ*2` → 構文エラー、`3 DEG`・`PI/10` → 正しい値。

試験: `pysad/tests/test_v2_P5_{measure,detr,micado_fitp_oy,specfn_cx,proton_beta,syntax,tk90}.py`。gs の値は `.gs.txt` と `pysad/tests/_p5_gscx.py` (特殊関数) を読む。
