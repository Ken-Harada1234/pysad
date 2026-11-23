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
| [flag_INS_DETR](flag_INS_DETR.sad) | INS の CALC のあとの `Twiss["DETR"]` (結合行列の行列式)。45 度回した四極 1 個の結合リング | ○ (全桁一致、DETR = 0.0021332) |
| [flag_DETR_ring_a](flag_DETR_ring_a.sad) | DETR ≥ 1 (モードが入れ替わる) の CELL・CALC・INS・Emittance[] (大きなリング A、3 GeV、603 m)。台本の先頭にラティスを連結して走らせる (`pysad/tests/kek-ls_v03_12.lat`) | ○ (β・DETR・Emittance[] は一致。一周の位相の整数部は k64 が ±1 ずれ (丸め誤差の符号)、v2 は gs と一致) |
| [flag_DETR_pf](flag_DETR_pf.sad) | 同 (PF リング、2.5 GeV、187 m)。先頭に `pf/lattice/v9_16_asym.sad` を連結 | 同上 |
| [flag_DETR_astrid2](flag_DETR_astrid2.sad) | 同 (ASTRID2、0.58 GeV、45.7 m)。先頭に `astrid/lattice/ASTRID2_v01_01.sad` を連結 | 同上 |
| [fn_MEASURE_fodo](fn_MEASURE_fodo.sad) | MEASURE・LOGMEASURE・TRCOD (3 GeV、周長 48 m の 12 セル FODO + スキュー四極 1 台) | ○ (k64 モードは gs の統計ブロック (4 桁) と全一致。gs は先頭 MARK の位置だけ出力し、TRCOD は SIGSEGV。v2 は種を戻す・重心が 0) |
| [fn_MEASURE_default_np](fn_MEASURE_default_np.sad) | NP を指定しない MEASURE。gs は粒子数 1000 (NP の既定) | ○ (k64 も同じ。v2 は NP の指定を要求する) |
| [fn_MEASURE_ring_a](fn_MEASURE_ring_a.sad) | MEASURE・TRCOD (大きなリング A)。先頭に `pysad/tests/kek-ls_v03_12.lat` を連結 | 同 fn_MEASURE_fodo |
| [fn_MEASURE_pf](fn_MEASURE_pf.sad) | MEASURE・TRCOD (PF リング)。先頭に `pf/lattice/v9_16_asym.sad` を連結 | 同 fn_MEASURE_fodo |
| [fn_MEASURE_astrid2](fn_MEASURE_astrid2.sad) | MEASURE・TRCOD (ASTRID2)。先頭に `astrid/lattice/ASTRID2_v01_01.sad` を連結 | 同 fn_MEASURE_fodo |
| [flag_TRCOD_orbit](flag_TRCOD_orbit.sad) | TRCOD の閉軌道の測定。24 セルの平面 FODO (3 GeV)、リングの先頭に薄いキック (L=0 の BEND、K0=1e-5) | CALC との差 7e-10、解析式との差 5e-6。gs は SIGSEGV (B42) |
| [flag_TRCOD_reference](flag_TRCOD_reference.sad) | TRCOD FIXCOD / MVCOD / STNCOD (基準軌道の保存と差 `Twiss["DDX"]`)、TRCOD ALIGN、ORBITJITTER。1 GeV の輸送路 | v2 のみ (K64 では TRCOD が落ちる。gs の参照なし) |
| [fn_PUTMEASURE_GETMEASURE](fn_PUTMEASURE_GETMEASURE.sad) | PUTMEASURE / GETMEASURE (測定値をファイルで受け渡し、fit の目標値にする) と LWAKE ON の MEASURE。実行すると `fn_PUTMEASURE_GETMEASURE.measure.dat` ができる (消してよい) | v2 のみ (K64 では MEASURE が無出力。gs の参照なし) |
| [flag_DELCOR_rms](flag_DELCOR_rms.sad) | DELCOR の規格化。v2 は DX の rms = delx、DY の rms = dely (相関長 5 m と 5 km、150 回の統計) | v2 の規格化。k64 は 1/√2 (rms の比が 0.71 前後) |
| [FFS_misc_DELCOR_WAKE_MULT_skew](FFS_misc_DELCOR_WAKE_MULT_skew.sad) | DELCOR (CAVI・MARK にも DX が入る)、ResponseMatrix の TRPT (応答に sqrt(p(補正磁石)/p(モニタ)) が掛かる)、旧式の WAKE コマンド (`WAKE L 名前 s w … ENDWAKE`)、MULT のスキュー成分 (B40・B62・B66) | v2 の出力だけ |

上の表のうち、拡張子つきのリンクの行 (`.sad`・`.py`) は説明の `.md` を持たない。内容と結果の要点は表と、台本・スクリプトの先頭のコメントに書いた。`.py` の例は、リポジトリ直下から `py -3.12 example/02_optics_cod/NAME.py` で走らせる (出力は `NAME.out.txt`)。gs と同じ動作で走らせるときは `PYSAD_COMPAT=k64`。

`flag_CALC4D_CALC6D` では `CALC6D` に RF の効果が入らず NZ に意味のある値が出なかった (gs・pysad とも)。シンクロトロンチューンは `01_lattice_elements/elem_CAVI` の `Emittance[]` で確認した。
