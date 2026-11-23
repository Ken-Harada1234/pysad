# example/01_lattice_elements の md から取り除いた履歴の記述

例の md を「現在の動作と結果の説明」にするため、版・経過・以前の状態の記述を取り除いた。取り除いた元の行を、ファイルごと・元の行番号つきで残す (行番号は整理前のもの)。この節は、整理前の md (リポジトリの先頭の版) と現在の md の差分から再構成した。

## FFS_RENUMBER.md

元の 29 行目:
> **○ 一致 (2026-11-17 に対応)。** `RENUMBER comp` (略 `RENUM`) は comp を起点に周回順で同名の要素の番号を 1 から振り直す (gs の `QF.23 D2.45 B.23 QD.23 SF.23 DS.46` と同じ)。`LINE["MULT"]` も振り直した番号になる。起点が無い (`RENUMBER;`、存在しない名前) と ` ???-FFS-Error-Missing origin component for RENUM_BER ` を出す。

元の 32 行目:
> - `RENUMBER QF` (名前だけ) は最初の QF が起点なので変化しない (gs と一致)。要素の INDEX 属性による番号の初期値は、最初の版では未対応だった。v2 で対応した (`docs/47_v2_changes.md`。gs の `LINE["NAME"]`・`RENUMBER` 後の名前と一致と記録されている)。この例では INDEX を使っておらず、再実行していない。RENUMBER のあと `LINE["POSITION","QF.1"]` が gs は `LINE["NAME"]` の表示と食い違い、pysad は表示と一致する位置を返す (意図した差、docs/47)。

## FFS_USE_EXPAND.md

元の 34 行目:
> **○ 一致 (2026-11-17 に対応)。** `USE EXPAND/NOEXPAND RING`、`CALC NOEXPAND`・`CALC EXPAND`・`GO NOEXPAND` が gs と同じ動作になった: NOEXPAND で成分ごとのずれ (QF.5 の K1=0.6、NX=6.637402) が残り、`???-FFS-Info-Element values are not expanded.` と `NOCALEXP` の表示が出る。`NOEXPAND` は KEEPEXP なので以後の `CALC` にも残る (`CALC EXPAND` か `USE` で戻る)。`SAVE`・`EXPAND` は代表 (最初の 1 台) の値を全成分に揃える。TrackParticles は成分ごとの値をそのまま使う (CALC・Emittance[] で揃う)。

元の 38 行目:
> - `NOKEEPEXP` (CALC のたびに tfvcomp で EXPAND を決める): 最初の版では未対応だった。v2 で gs に合わせて対応した (`docs/47_v2_changes.md`。gs の実測と一致と記録されている)。ただし pysad の `tfvcomp` は常に偽 (成分ごとの FREE が無い)。この例では使っておらず、再実行していない。

## README.md

元の 16 行目:
> | [elem_APERT](elem_APERT.md) | APERT。追跡での開口 (flag) | ○ (v2 で対応。`TYPE` の種別表示だけ違う) |

## elem_APERT.md

元の 31 行目:
> **○ 一致 (2026-11-17、v2.0.0 の再実行)。** flag と落ちた粒子の x は 1 回目・2 回目とも gs と同じ (内側は x=1.715 mm、flag=1。外側は x=20 mm で flag=0)。この例の最初の版では pysad が APERT を落とさず、外側の粒子が x=37.8 mm まで進んでいた (v1 の未対応)。v2 で APERT の追跡を実装した (`docs/47_v2_changes.md`)。互換モードによらず同じ結果 (`PYSAD_COMPAT=k64` でも上の表のとおり)。

## elem_CAVI.md

元の 37 行目:
> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率の X と Z・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## elem_SEXT_OCT_DECA_DODECA.md

元の 40 行目:
> - `TYPE D6` の `K5 =10000` も一致 (2026-11-17 に修正: pysad は DODECA の K5 を長さあたりの内部値 (50000) のまま返していた。`LINE["K5"]`・`Element["K5"]`・TYPE・ATTRIBUTE は積分強度 10000)。

## fn_LINE.md

元の 33 行目:
> **○ 一致 (2026-11-17 に対応)。** `LINE["MULT"]` は gs と同じく同名の要素の通し番号 (QF.1..QF.3 が 1, 2, 3)。同名が 1 台だけの要素と `$$$` は 0。ほかの項目は全桁一致。


# 台本 (`.sad`・`.py`) の冒頭コメントから、機能別フォルダへの移動時に取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `elem_BEND_F1_adjacent.sad` 1 行目:
  - 元: `! v2next A1 項目 4: BEND の F1・FRINGE + RAD の追跡。K64 の tfgeo (inext・iprev) は、隣り合う同種の要素 (間が長さ 0 の要素だけ) を繋ぎ、`
  - 後: `! BEND の F1・FRINGE + RAD の追跡。K64 の tfgeo (inext・iprev) は、隣り合う同種の要素 (間が長さ 0 の要素だけ) を繋ぎ、`
- `elem_BEND_F1_adjacent.sad` 2 行目:
  - 元: `! tbend は繋がっていない端にだけ F1 の断片 (長さ F1/2 の 2 つの放射) を入れる。旧 pysad は常に両端に入れた (相対 3.6%、1 個だけの輸送路で 3%)。`
  - 後: `! tbend は繋がっていない端にだけ F1 の断片 (長さ F1/2 の 2 つの放射) を入れる。`
- `elem_CAVI_offset.sad` 1 行目:
  - 元: `! v2next A1 項目 5: CAVI の DX・DY・ROTATE と V1・V20・V11・V02 (tcav.f の TENT・TEXIT と v1a…v02a)。`
  - 後: `! CAVI の DX・DY・ROTATE と V1・V20・V11・V02 (tcav.f の TENT・TEXIT と v1a…v02a)。`
- `elem_CAVI_offset.sad` 2 行目:
  - 元: `! 旧 pysad は CAVI の DX・DY・ROTATE を追跡で無視し (x = y = 0 の粒子に入る縦の蹴り y 6e-9 が無く、y が 1e-3 違った)、V1 V20 V11 V02 はパーサが捨てた。`
  - 後: `!`
- `FFS_USE_after_EMIT.sad` 1 行目:
  - 元: `! v2next A1 項目 6: EMIT (Emittance[]・CALC6D) のあとの USE。K64 の trf0 (同期位置の時間) は大域変数で、temit と CALC6D だけが更新する。`
  - 後: `! EMIT (Emittance[]・CALC6D) のあとの USE。K64 の trf0 (同期位置の時間) は大域変数で、temit と CALC6D だけが更新する。`
- `FFS_USE_after_EMIT.sad` 2 行目:
  - 元: `! USE では 0 に戻らない。旧 pysad は USE で 0 に戻し、追跡の RF の位相が EMIT 前の値に戻った (δ が 87% 違った)。gs と全桁一致。`
  - 後: `! USE では 0 に戻らない。gs と全桁一致。`
- `elem_BEND_DX_offset.sad` 1 行目:
  - 元: `! v2next A2 項目 4: ずれた BEND の CALC (tbende の tchge) と追跡 (tbend の tbshift) の差。矩形磁石 (E1=E2=0.5) 1 台を DX=2 mm ずらす。`
  - 後: `! ずれた BEND の CALC (tbende の tchge) と追跡 (tbend の tbshift) の差。矩形磁石 (E1=E2=0.5) 1 台を DX=2 mm ずらす。`
- `elem_BEND_K1_edge.py` 2 行目:
  - 元: `"""v2next B1 項目 2: K1≠0 の BEND の CALC (tbende) と追跡 (tbendi) の写像の差。どちらが正しいかを、独立な「厳密な写像」(ローレンツ力の数値積分) との`
  - 後: `"""K1≠0 の BEND の CALC (tbende) と追跡 (tbendi) の写像の差。どちらが正しいかを、独立な「厳密な写像」(ローレンツ力の数値積分) との`
- `syntax_MAIN_expression.sad` 2 行目:
  - 元: `! v2next D2 項目 1: MAIN の要素定義の式 (PI・関数・単位・変数) を gs と同じ文法・同じ値で読む`
  - 後: `! MAIN の要素定義の式 (PI・関数・単位・変数) を gs と同じ文法・同じ値で読む`
- `elem_BEND_steer_AE_poles.py` 2 行目:
  - 元: `"""v2next E 項目 2 (B108): ANGLE=0 の BEND (補正磁石・挿入光源の極) の端の角 AE1・AE2 と放射減衰の分配 Jx。`
  - 後: `"""(B108) ANGLE=0 の BEND (補正磁石・挿入光源の極) の端の角 AE1・AE2 と放射減衰の分配 Jx。`
- `elem_MULT_F1.sad` 1 行目:
  - 元: `! v2next E 項目 2: MULT の F1・F2・F1K1F・F1K1B・FRINGE (四極の線形フリンジ tqlfre)・K0 つきの FB1・FB2 (tbfrme の線形フリンジ) の Emittance[] と CALC。`
  - 後: `! MULT の F1・F2・F1K1F・F1K1B・FRINGE (四極の線形フリンジ tqlfre)・K0 つきの FB1・FB2 (tbfrme の線形フリンジ) の Emittance[] と CALC。`
- `elem_MULT_F1.sad` 2 行目:
  - 元: `! 旧 pysad は F1・FRINGE があると NotImplementedError (FB1・FB2 は黙って無視) だった。tmulte の入口 (tqfrie, tbfrme(fb1), tqlfre(f1in, f2in)、そのあと tradke(f1in) か tsetr0)、`
  - 後: `!tmulte の入口 (tqfrie, tbfrme(fb1), tqlfre(f1in, f2in)、そのあと tradke(f1in) か tsetr0)、`
- `elem_BEND_steer_F1.sad` 1 行目:
  - 元: `! v2next E 項目 2: ANGLE=0 の BEND (補正磁石、K64 の tsteee) の F1・FB1・FB2・FRINGE・AE1・AE2 の Emittance[] と CALC。`
  - 後: `! ANGLE=0 の BEND (補正磁石、K64 の tsteee) の F1・FB1・FB2・FRINGE・AE1・AE2 の Emittance[] と CALC。`
- `elem_BEND_steer_F1.sad` 2 行目:
  - 元: `! 旧 pysad は F1 を黙って無視するか (Emittance[]) NotImplementedError だった。tsteee の tblfre (dxfr1 = -fb1^2/rhob/24 ...)、入口の断片 2 つ (n = -1, 0)、`
  - 後: `!tsteee の tblfre (dxfr1 = -fb1^2/rhob/24 ...)、入口の断片 2 つ (n = -1, 0)、`
- `elem_BEND_K1_face.py` 2 行目:
  - 元: `"""v2next P2 項目 5: K1≠0 の BEND の端の模型。K64 の tbende・tbendi は、双極磁場を端面から、四極 (K1) の磁場を半径方向の入口・出口の面から始める`
  - 後: `"""K1≠0 の BEND の端の模型。K64 の tbende・tbendi は、双極磁場を端面から、四極 (K1) の磁場を半径方向の入口・出口の面から始める`
- `elem_AE_exact_geometry.py` 2 行目:
  - 元: `"""v2next P3 項目 2・3 (B108 の残り): 端の角 AE1・AE2 のある極 (K0 を持つ ANGLE=0 の要素) の ΔJx を、厳密な幾何 (円弧・微分方程式) + ロビンソンの定理と比べる。`
  - 後: `"""(B108 の残り) 端の角 AE1・AE2 のある極 (K0 を持つ ANGLE=0 の要素) の ΔJx を、厳密な幾何 (円弧・微分方程式) + ロビンソンの定理と比べる。`
- `elem_BEND_steer_F1_maxwell_ramp.py` 2 行目:
  - 元: `"""v2next P3 項目 1: 補正磁石 (ANGLE=0 の BEND、K64 の tsteer) の F1 (端の磁場が直線で立ち上がる長さ) の追跡を、独立な基準 (x 方向に一様な Maxwell 場の端の厳密な運動) と比べる。`
  - 後: `"""補正磁石 (ANGLE=0 の BEND、K64 の tsteer) の F1 (端の磁場が直線で立ち上がる長さ) の追跡を、独立な基準 (x 方向に一様な Maxwell 場の端の厳密な運動) と比べる。`
- `elem_MULT_AE_emittance.sad` 1 行目:
  - 元: `! v2next P3 項目 2: MULT の端の角 AE1・AE2 の CALC・Emittance[] (K64 の tmulke の apsi)。旧 pysad は MULT の AE を読まず黙って無視していた (K64 の動作の追加)。`
  - 後: `! MULT の端の角 AE1・AE2 の CALC・Emittance[] (K64 の tmulke の apsi)。`
- `elem_MULT_F1_AE_track.sad` 1 行目:
  - 元: `! v2next P3 項目 1・2: TrackParticles の MULT の F1・F2・FB1・FB2・FRINGE (K64 の tmulti の tqlfre 相当・tblfri) と AE1・AE2 (tmulk)。旧 pysad は F1・FRINGE で`
  - 後: `! TrackParticles の MULT の F1・F2・FB1・FB2・FRINGE (K64 の tmulti の tqlfre 相当・tblfri) と AE1・AE2 (tmulk)。旧 pysad は F1・FRINGE で`
- `elem_BEND_steer_ring_a.py` 2 行目:
  - 元: `"""v2next P3: 大きなリング A (pysad/tests/kek-ls_v03_12.lat、3 GeV、周長 603 m、3563 要素) の先頭に補正磁石 2 台 (K0 = +-5e-4 rad、L = 0.2 m) を入れて確かめる。`
  - 後: `"""大きなリング A (pysad/tests/kek-ls_v03_12.lat、3 GeV、周長 603 m、3563 要素) の先頭に補正磁石 2 台 (K0 = +-5e-4 rad、L = 0.2 m) を入れて確かめる。`
- `elem_BEND_steer_AE2_tune.sad` 1 行目:
  - 元: `! v2next P3 項目 2・3: K64 の tsteer (追跡) は出口の端の角 AE2 を FB2 != 0 のときだけ入れる (CALC・Emittance[] の tsteee は常に入れる)。`
  - 後: `! K64 の tsteer (追跡) は出口の端の角 AE2 を FB2 != 0 のときだけ入れる (CALC・Emittance[] の tsteee は常に入れる)。`
- `elem_BEND_steer_F1_track.sad` 1 行目:
  - 元: `! v2next P3 項目 1: TrackParticles の F1 (ANGLE=0 の BEND = ステアリング、K64 の tsteer)。旧 pysad は NotImplementedError。`
  - 後: `! TrackParticles の F1 (ANGLE=0 の BEND = ステアリング、K64 の tsteer)。`
- `elem_UND.sad` 先頭付近:
  - 元: `! v2 (B37):`
  - 後: `! (B37)`
- `elem_PHSROT.sad` 先頭付近:
  - 元: `! v2 (B38):`
  - 後: `! (B38)`
- `syntax_MAIN_assign.sad` 先頭付近:
  - 元: `! (担当 C2) 台本の読み方`
  - 後: `! 台本の読み方`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- この移動で、旧 `v2_next` の `*.pysad_k64.txt` のうち、このフォルダに移る例の分を削除した (k64 モードの出力は `PYSAD_COMPAT=k64` で再現できる。`flag_PSPAC_k64_shape` だけは k64 オプションの例として残した)。
