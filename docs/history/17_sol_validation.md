# 17_sol の検証の詳細 (gs との比較)

docs/17_sol.md の「検証」の詳細。旧 docs/44_sol_guide.md §6 と CALC6D の確認の数値、旧 docs/29_sol_optics.md の検証を元の文脈つきで移した。

## 旧 44 §6 gs との一致の度合い


数値は [`docs/17_sol.md`](17_sol.md) と各試験から (いずれも gs = K64 2026-10-19 ビルドの出力を参照に取ったもの。特記なければ NORAD または RAD ON、領域内の要素は DRIFT・QUAD・SOL)。

| 項目 | 試験 | 条件 | 一致 |
| --- | --- | --- | --- |
| Twiss (チューン、β、分散) | [`test_sol_optics.py`](../pysad/tests/test_sol_optics.py) 4 本 | BZ=0.3 (1.7 GeV)、BZ=1.2 (3 GeV、回転四極・領域内 SOL)、INS・RAD ON・DP≠0 | チューン 1e-15、β 3e-14、領域内の β α 分散 結合 1e-11 |
| Emittance[] e1 | [`test_sol_emit.py`](../pysad/tests/test_sol_emit.py) | BZ=2、3 GeV、12 セルのリング、RFSW | ε_x 2e-12、ε_y 4e-10、ε_z 1e-7、チューン 1e-13、軌道 1e-11 (相対) |
| Emittance[] e2 | 同上 | BZ=5、K1=0.3 L=0.5、軌道 2.7 mm | ε_x 4e-9、ε_y 2e-8、ε_z 3e-6、チューン 1e-11 |
| 幾何 | [`test_sol_geo.py`](../pysad/tests/test_sol_geo.py) (sg2 sg3) | GEO が始端 (前向き) と終端 (後ろ向き)、DX DY DPX DPY あり、3 GeV | 18 要素の S GX GZ GCHI2 GCHI3 が相対 1e-9 (絶対 1e-13) |
| 追跡 NORAD | `bench_level9` (Level 9) | BZ=2、2.5 GeV、領域に D1 + QD | 5e-11 以下 |
| 追跡 RAD (NOFLUC) | 同上 | 同上 | 1.5e-10 以下 (**ただし保存済みの参照に対して。下の限界 1 を参照**) |
| スピン追跡 | [`test_spin_track.py`](../pysad/tests/test_spin_track.py) | SOL 領域内の DRIFT (BZ=0.05) | 4e-12 以下 (四極は照合できていない) |

(追跡の NORAD・RAD の 5e-11・1.5e-10 は [`docs/10_tracking.md`](10_tracking.md) の Level 9 の記録)


## 旧 44 §4 CALC6D の確認の数値 (gs、SOL 領域のあるリング、3 GeV、BZ = 2 T、回転四極、QFE の DX・DY)

CALC 4D 2.73475252978800・1.87531224259842 (pysad k64 2.73475252978794)、CALC6D 2.73495437251036・1.87554377922481・NZ −.01173531080122 (2.73495437251039 など。1e-9、NZ は絶対 1e-9)、RADCOD つき CALC6D 2.73497629185364 (1e-10)、RADCOD の CALC6D のあとの 4D の CALC 2.73465818826401 (pysad 2.73465818826743。CELL CALC6D のあと `DP0` を先頭の `DDP` にする)、β (1・30・70 番) も 1e-9。INTRA + SOL の Emittances は {2.706114556595e-7、2.245979910482e-8、2.778701736555e-5} で、gs と εx・εy 2e-8、εz 9e-8 (相対)。TouschekLifetime は gs・pysad とも NaN (この条件で確認できず)。

## 旧 29 検証 (gs との比較)


### Twiss (`pysad/tests/test_sol_optics.py` 4 本)

- 単純な領域 (BZ=0.3): チューン 1e-15、β 3e-14、分散 1e-14。
- 強い場 (BZ=1.2、領域内に回転四極 QS、領域内の SOL SM) の CELL: チューン 4e-16。領域内 8 要素の β、α、分散、結合 (R2 R3) が 1e-11 以内。
- INS・RAD ON・DP≠0: チューン 1e-15、領域内の Twiss が 1e-11 以内。

### Emittance[] (`pysad/tests/test_sol_emit.py` 3 本。12 セルのリングの 1 セルに領域、RFSW)

- e1 (BZ=2、回転四極、領域内の SOL、QFE の DX DY で軌道が領域を斜めに通る): ε_x 2e-12、ε_y 4e-10、ε_z 1e-7、チューン 1e-13、軌道 1e-11 (相対)。
- e2 (BZ=5、長い四極 K1=0.3 L=0.5、大きい軌道 2.7 mm): ε_x 4e-9、ε_y 2e-8、ε_z 3e-6、チューン 1e-11。
- e1z (BZ=0 の領域): 普通の要素になり gs と 1e-10。
- 要素ごとの切り分け (SA SB のみ、SA D3 SB、SA QS SB、回転なし、DISFRIN=1) もすべて一致した。


## 旧 29 幾何の検証


- 検証 ([`pysad/tests/test_sol_geo.py`](../pysad/tests/test_sol_geo.py)): GEO が始端 (前向き、回転四極 K1>0、領域内の SOL) と終端 (後ろ向き、K1<0) の 2 通りで、S・GX・GZ・GCHI2・GCHI3 の 18 要素が gs と 1e-9 (絶対 1e-13) で一致。
