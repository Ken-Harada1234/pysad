# example/02_optics_cod の md から取り除いた履歴の記述

例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきで残す (行番号は整理前のもの)。この節は、整理前の md (リポジトリの先頭の版) と現在の md の差分から再構成した。

## FFS_BEAMSIZE.md

元の 39 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

元の 40 行目:
> なお、保存してあった旧い `.pysad.txt` には `[pysad] DISP のモード B は未対応` と書いてあったが、これは「`DISP BEAM` は 2026-11-17 に対応」の前の出力で古かった。今回の再実行で `DISP BEAM` の表が出ている。

元の 43 行目:
> - `DISP BEAM` (`DISP B`、EMIT 後) は 2026-11-17 に対応 (gs と同じ見出しと桁。射影 Twiss の列 AXp BXp EMITXp… を含む。雑音 1e-16 の列は pysad が 0)。`LINE["SIGX",...]` は pysad でも gs と同じ値を返す (ただし pysad は EMIT を先に実行していないと `BEAMSIZE` が `EMIT を先に実行してください` と止まる)。

## FFS_DISPLAY.md

元の 26 行目:
> **○ 一致 (2026-11-17 に対応)。** 標準・`R`・`O`・`E`・`G` に加え `GA`・`Z`・`P`・`A`・`B` (`DISP BEAM`、EMIT 後)・`D` (28 光学関数の 1 行ダンプ)・`REFERENCE`・`DREFERENCE` の表を gs と同じ見出し・桁で出す (表の数値は全桁一致。雑音 1e-16 の列だけ pysad は 0 を出す)。`DISP B REFERENCE` は gs と同じ `???-FFS-Info-REF and DREF not implemented for DISP B.`。

## FFS_DUMP.md

元の 32 行目:
> **○ 一致 (2026-11-17 に対応)。** 差の原因は gs の `DUMP` が語を 20 文字 (`character*20`) で切ること: 長い `QF.1|QF.2|QF.3|QD.2|QD.3|D1.1|B.1|Q*.4` は先頭の `QF.1|QF.2|QF.3|QD.2|` だけが使われるので D1.1 (DL) と B.1 (DDK) は出なかった (短い `DUMP QF.1|D1.1|B.1;` は gs も出す)。pysad も 20 文字で切る。どの要素にも合わない語のときは全要素 (誤差のあるもの) を出す。

元の 33 行目:
> もう 1 点: `DDK` は K0 のキーワードのある要素 (BEND・MULT・TCAVI) だけに入る (gs で QUAD に入らないのを確認。pysad は QUAD にも入れていた)。`QF.1|QF.2` の 2 つ目に誤差が入らないのは、誤差コマンドが `*%{}` を含まない語では最初の一致だけに作用するため (gs・pysad とも)。

## FFS_REFERENCE.md

元の 23 行目:
> **○ 一致 (2026-11-17 に対応)。** `REFERENCE` (略 `REF`) は直前の CALC の光学を参照光学にする。`Twiss["RBX",…]` (参照)、`Twiss["DBX",…]` (現在 - 参照、BX BY BZ は比)、`Twiss["RALL"/"DALL",…]`、`DISP REFERENCE`・`DREFERENCE` (見出しの文字の置換も含む) が gs と全桁一致。参照光学の位置は整数位置だけ (小数位置はエラー)。REFERENCE 前の参照は 0。

## fn_OrbitGeo_GeoBase.md

元の 33 行目:
> - `OrbitGeo[$$$]` も一致 (2026-11-17 に対応: GX = -5.338e-05、GY = 1.2836e-4、GZ = -3.354e-4)。原因は `OrbitGeo` の式 (Optics.n) が一周の経路長のずれ `Twiss["DZ",$$$]` (= -5.338e-5) を進行方向に足すことで、pysad は最後の行 (`$$$`) の DZ を 0 にしていた。`$$$` の DZ は周回の出口の値にした (CELL の 4D の CALC)。

## fn_TransferMatrix_BeamMatrix.md

元の 38 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## fn_Twiss.md

元の 40 行目:
> 値はすべて一致 (全桁)。2026-11-17: `Twiss["LENGTH"]` (成分数。引数を付けると `General::narg`)、`Twiss["GAMMA",1]`・`"GAMMABETA"`・`"S"`・`"SIZE…"` (LINE と同じ値) に対応した。光学関数でない名前は `General::wrongval:  #1 (NAME) is must be to be name of optical function` (gs と同じ)。


# 台本 (`.sad`・`.py`) の冒頭コメントから、機能別フォルダへの移動時に取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `flag_INS_DETR.sad` 1 行目:
  - 元: `! v2next A1 項目 1 (N12): INS の CALC の後の Twiss["DETR"] (結合行列の行列式)。旧 pysad は 0 を返した (k64・v2 とも、pysad 自身の不具合)。`
  - 後: `! INS の CALC の後の Twiss["DETR"] (結合行列の行列式)。`
- `flag_TRCOD_orbit.sad` 1 行目:
  - 元: `! v2next C3 項目 1: TRCOD の閉軌道の測定。24 セルの平面 FODO (3 GeV)、リングの先頭に薄いキック (L=0 の BEND、K0=1e-5)。`
  - 後: `! TRCOD の閉軌道の測定。24 セルの平面 FODO (3 GeV)、リングの先頭に薄いキック (L=0 の BEND、K0=1e-5)。`
- `flag_DETR_astrid2.sad` 1 行目:
  - 元: `! v2next P5 項目 2: DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (ASTRID2 (0.58 GeV、周長 45.7 m))。このファイルは台本の後半だけ: 先頭に example/astrid/lattice/ASTRID2_v01_01.sad を連結して走らせる。`
  - 後: `! DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (ASTRID2 (0.58 GeV、周長 45.7 m))。このファイルは台本の後半だけ: 先頭に example/astrid/lattice/ASTRID2_v01_01.sad を連結して走らせる。`
- `flag_DETR_pf.sad` 1 行目:
  - 元: `! v2next P5 項目 2: DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (PF リング (2.5 GeV、187 m))。このファイルは台本の後半だけ: 先頭に example/pf/lattice/v9_16_asym.sad を連結して走らせる。`
  - 後: `! DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (PF リング (2.5 GeV、187 m))。このファイルは台本の後半だけ: 先頭に example/pf/lattice/v9_16_asym.sad を連結して走らせる。`
- `flag_DETR_ring_a.sad` 1 行目:
  - 元: `! v2next P5 項目 2: DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m))。このファイルは台本の後半だけ: 先頭に pysad/tests/kek-ls_v03_12.lat を連結して走らせる。`
  - 後: `! DETR >= 1 (モードが入れ替わる) の CALC・INS・Emittance[] の gs との比較 (大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m))。このファイルは台本の後半だけ: 先頭に pysad/tests/kek-ls_v03_12.lat を連結して走らせる。`
- `fn_MEASURE_astrid2.sad` 1 行目:
  - 元: `! v2next P5 項目 1: MEASURE / TRCOD の gs との比較 (ASTRID2 (0.58 GeV、周長 45.7 m))。このファイルは台本の後半だけ: 先頭に example/astrid/lattice/ASTRID2_v01_01.sad を連結して走らせる。`
  - 後: `! MEASURE / TRCOD の gs との比較 (ASTRID2 (0.58 GeV、周長 45.7 m))。このファイルは台本の後半だけ: 先頭に example/astrid/lattice/ASTRID2_v01_01.sad を連結して走らせる。`
- `fn_MEASURE_default_np.sad` 1 行目:
  - 元: `! v2next P5 項目 1: NP を指定しない MEASURE。gs は粒子数 1000 (NP の既定、tffsa.f の np0=1000)。k64 も同じ (v2 は NP の指定を要求する)。`
  - 後: `! NP を指定しない MEASURE。gs は粒子数 1000 (NP の既定、tffsa.f の np0=1000)。k64 も同じ (v2 は NP の指定を要求する)。`
- `fn_MEASURE_fodo.sad` 1 行目:
  - 元: `! v2next P5 項目 1: MEASURE / TRCOD (gs の実動作との比較)。条件: 3 GeV、周長 48 m の 12 セル FODO + スキュー四極 1 台 (結合あり)。`
  - 後: `! MEASURE / TRCOD (gs の実動作との比較)。条件: 3 GeV、周長 48 m の 12 セル FODO + スキュー四極 1 台 (結合あり)。`
- `fn_MEASURE_pf.sad` 1 行目:
  - 元: `! v2next P5 項目 1: MEASURE / TRCOD の gs との比較 (PF リング (2.5 GeV、187 m))。このファイルは台本の後半だけ: 先頭に example/pf/lattice/v9_16_asym.sad を連結して走らせる。`
  - 後: `! MEASURE / TRCOD の gs との比較 (PF リング (2.5 GeV、187 m))。このファイルは台本の後半だけ: 先頭に example/pf/lattice/v9_16_asym.sad を連結して走らせる。`
- `fn_MEASURE_ring_a.sad` 1 行目:
  - 元: `! v2next P5 項目 1: MEASURE / TRCOD の gs との比較 (大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m))。このファイルは台本の後半だけ: 先頭に pysad/tests/kek-ls_v03_12.lat を連結して走らせる。`
  - 後: `! MEASURE / TRCOD の gs との比較 (大きなリング A (kek-ls_v03_12.lat、3 GeV、603 m))。このファイルは台本の後半だけ: 先頭に pysad/tests/kek-ls_v03_12.lat を連結して走らせる。`
- `fn_PUTMEASURE_GETMEASURE.sad` 先頭付近:
  - 元: `! v2 (C2): PUTMEASURE`
  - 後: `! PUTMEASURE`
- `flag_TRCOD_reference.sad` 先頭付近:
  - 元: `! v2 (C2): TRCOD`
  - 後: `! TRCOD`
- `flag_DELCOR_rms.sad` 先頭付近:
  - 元: `! v2 (担当 C2): DELCOR`
  - 後: `! DELCOR`
- `FFS_misc_DELCOR_WAKE_MULT_skew.sad` 先頭付近:
  - 元: `! v2 の小さな修正: DELCOR`
  - 後: `! DELCOR`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- この移動で、旧 `v2_next` の `*.pysad_k64.txt` のうち、このフォルダに移る例の分を削除した (k64 モードの出力は `PYSAD_COMPAT=k64` で再現できる。`flag_PSPAC_k64_shape` だけは k64 オプションの例として残した)。
