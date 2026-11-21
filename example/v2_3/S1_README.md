# example/v2_3: 担当 S1 の例 (周波数マップの検証、v2.3.0、2026-11-20)

台本ごとのまとめと全数値は [`FMAP_summary.md`](FMAP_summary.md)、解説は `docs/48_frequency_map_lma.md` (§3.5・§9)、記録は `docs/history/v23_S1.md`、試験は `pysad/tests/test_v23_S1_fmapres.py` (12 件)、コードは `pysad/fmapres.py`。
実行はリポジトリ直下から `py -3.12 example/v2_3/S1_xxx.py ...` (標準出力をファイルへ落とすときは `PYTHONUTF8=1`)。スクリプトは出力を `example/v2_3/` へ書く (保存済みの出力を残したいときは複製したフォルダで走らせる)。
使ったラティス: `example/astrid/lattice/ASTRID2_v01_01.sad`・`ASTRID2_v01_02_id.sad` (MPW)、`pysad/tests/kek-ls_v03_12.lat`、`example/00_lattice/FODO_ring.sad`。

| ファイル | 内容 |
|---|---|
| `S1_fmap_rdt.py` / `S1_fmap_rdt_{fodo,astrid2,kekls}.{out.txt,png}`、`S1_fmap_rdt_fodo.table.txt` | (a) RDT・ADTS から予想した共鳴の領域と、D の大きい領域の比較 |
| `S1_fmap_variants.py` / `S1_fmap_variants_<モード>_*.out.txt`、`S1_fmap_{delta,k64,errors,mpw,grid,lambda}_*.table.txt`、`.png` | (b) δ、(c) k64、(d) 誤差、(e) MPW、(f) NAFF・格子、(a の原因) 6 極の λ 倍。モード: `delta <ring> xdp\|ydp`、`k64 <ring>`、`errors <ring> <level>`、`mpw`、`naff <ring>`、`grid <ring> <nx> <ny>`、`lambda <ring>` |
| `S1_fmap_grid_summary.py` / `.out.txt` | (f) 81 × 61 の格子から部分格子の統計 |
| `S1_gs_fmap_compare.py` / `S1_gs_fmap_kekls.{sad,gs.txt,pysad.txt,out.txt}` | (g) KEK-LS の追跡の gs との比較 |
