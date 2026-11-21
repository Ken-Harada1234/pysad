# 02 線形光学と COD

CALC・チューン・クロマティシティ・転送行列・ビームサイズ・閉軌道・幾何。ラティスは `00_lattice/FODO_ring.sad`、閉軌道は新規の `00_lattice/errors_ring.sad` (ずれのある四極 2 台)。解析値の検算は `tools/ana_fodo_ring.py` (2x2 行列)。

各例は 4 ファイル: `NAME.sad` (スクリプト)、`NAME.gs.txt` (gs の画面表示)、`NAME.pysad.txt` (pysad の画面表示)、`NAME.md` (説明)。
走らせ方: `py -3.12 example/tools/run_example.py example/02_optics_cod/NAME.sad`。凡例: ○ 一致 (表示の桁・書式の差は除く)、△ 差あり (内容を `.md` に書いた)、× gs 失敗、× pysad 未対応。

| ファイル | 内容 | gs と pysad |
|---|---|---|
| [FFS_CALCULATE](FFS_CALCULATE.md) | CALC。チューン・Twiss と 2x2 行列の解析値 | ○ |
| [FFS_CHROMATICITY](FFS_CHROMATICITY.md) | CHRO と DP0 の有限差分 (ξ)、6 極の寄与 | ○ |
| [FFS_DISPLAY](FFS_DISPLAY.md) | DISP の各モード (標準 R O E G GA Z) | ○ |
| [FFS_DUMP](FFS_DUMP.md) | 誤差コマンド (DELX DELY DTHETA DL DK DDK) と DUMP | ○ (gs の DUMP は語を 20 文字で切る) |
| [FFS_MATRIX](FFS_MATRIX.md) | MATRIX。1 セルの行列、cos μ | ○ |
| [FFS_REFERENCE](FFS_REFERENCE.md) | REFERENCE。参照光学との差 | ○ |
| [FFS_COUPLE](FFS_COUPLE.md) | COUPLE と DECOUPLE | △ (DECOUPLE 後の gs の値が式と合わない。pysad は ATTRIBUTE に結合を出さない) |
| [FFS_BEAMSIZE](FFS_BEAMSIZE.md) | BEAMSIZE。σ² = εβ + (ησδ)² | ○ |
| [fn_Twiss](fn_Twiss.md) | Twiss 関数。補間・ALL・PEX・GMX | ○ |
| [fn_TransferMatrix_BeamMatrix](fn_TransferMatrix_BeamMatrix.md) | TransferMatrix と BeamMatrix | ○ |
| [flag_CELL_INS](flag_CELL_INS.md) | CELL (周期解) と INS (初期値からの伝搬) | ○ |
| [flag_COD](flag_COD.md) | 誤差入りリング (errors_ring) の閉軌道と 2x2 行列の式 | ○ |
| [fn_OrbitGeo_GeoBase](fn_OrbitGeo_GeoBase.md) | OrbitGeo と GeoBase | ○ |
| [flag_GEOCAL_GEOFIX](flag_GEOCAL_GEOFIX.md) | GEOCAL と GEOFIX。幾何の更新 | ○ |
| [flag_CALC4D_CALC6D](flag_CALC4D_CALC6D.md) | CALC4D と CALC6D。6x6 行列の M56 = -α_c C | ○ |

`flag_CALC4D_CALC6D` では `CALC6D` に RF の効果が入らず NZ に意味のある値が出なかった (gs・pysad とも)。シンクロトロンチューンは `01_lattice_elements/elem_CAVI` の `Emittance[]` で確認した。
