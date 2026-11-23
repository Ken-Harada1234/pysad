# 34_sad_packages から取り除いた履歴の記述

`docs/34_sad_packages.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## # SAD の Packages/*.n を SADScript のまま使う (軌道補正・最適化・数学パッケージ)

取り除いた部分: 「最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。」

元の行:

> 最終更新: 2026-11-18 (v2.1.0 の後の現状のスナップショット。経緯は `history/` にある)

## # SAD の Packages/*.n を SADScript のまま使う (軌道補正・最適化・数学パッケージ)

取り除いた部分: 「「LDSolve・Optics.n」「v2.1.0 の確認」)」

元の行:

> | LDSolve.n、Optics.n | LDSolve、Optics の TransferMatrix・BeamMatrix・NormalCoordinates・OrbitGeo | gs と一致 (下記「LDSolve・Optics.n」「v2.1.0 の確認」)。PeriodicR も gs と 3e-16 |

## ## UpValue を定義するパッケージ (MakeBump・Rational・ElementValues)

取り除いた部分: 「は v2.1.0 で実装した (docs/04)。**未対応**」

元の行:

>   `VISIT` での退避と `BYE` での復元は v2.1.0 で実装した (docs/04)。**未対応**: フィットの変数との連動 (`CouplingMatrix`・`CoupledKeys` を使う Fortran 側 `tffssetupcouple` は呼ばない)。

## ## LDSolve・Optics.n (2026-11-16 夜に gs と比較)

取り除いた部分: 「pysad は `DiagonalMatrix` が複素数の列を受けず失敗していた (`float() argument must be … not 'complex'`) ので直した。あわせて行列の演算 (Dot・Inverse) の結果は、gs が要素ごとに虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻すようにした (以前は全要素が実のときだけ)。」

元の行:

> - LDSolve.n: pysad は `DiagonalMatrix` が複素数の列を受けず失敗していた (`float() argument must be … not 'complex'`) ので直した。あわせて行列の演算 (Dot・Inverse) の結果は、gs が要素ごとに虚部が 0 の複素数を実数に戻すので、pysad も要素ごとに戻すようにした (以前は全要素が実のときだけ)。

## ## LDSolve・Optics.n (2026-11-16 夜に gs と比較)

取り除いた部分: 「直したあと、LDSolve の終点」

元の行:

>   直したあと、LDSolve の終点 (cos 1, −sin 1) と 3 通りの呼び出しが gs と一致 (試験は [`test_evalfixes.py`](../pysad/tests/test_evalfixes.py) の 2 件)。

## ## LDSolve・Optics.n (2026-11-16 夜に gs と比較)

取り除いた部分: 「- `PeriodicR` は v2.1.0 で比較した (下の「v2.1.0 の確認」。」

元の行:

> - `PeriodicR` は v2.1.0 で比較した (下の「v2.1.0 の確認」。`PeriodicOptics` を先に呼ぶ必要がある)。

## ## v2.1.0 の確認 (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「## v2.1.0 の確認 (小さな台本を gs と行ごとに突き合わせた。k64 モード)」

元の行:

> ## v2.1.0 の確認 (小さな台本を gs と行ごとに突き合わせた。k64 モード)

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「| 一致の状態 | 直したもの・残った差 |」

元の行:

> | ファイル | 比べた関数 | 一致の状態 | 直したもの・残った差 |

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「**`Log2$`**: 組み込みが無く未評価だった → 追加 (gs は」

元の行:

> | `functions.n` | `Log2$`、`InverseErfc`、`Names`、`RotateRight`・`RotateLeft`、`Union` (`SameTest`)、`ListQ`、`MatrixQ`、`Diagonal`、`Count`、`MemberQ`・`FreeQ` (レベルつき)、`Nest`・`NestList`・`Fold`・`FoldList`、`Through`、`FractionalPart`、`Positive` 系、`IntegerQ`、`Subtract`・`Divide`、`Xor`、`Append`・`Prepend`、`HeldPart`、`Evaluate`、`StringQ`・`Characters`・`StringDrop`・`StringInsert`、`Timing`、`For` | `Log2$`・`Union` (`SameTest`)・`Names` 以外は一致 | **`Log2$`**: 組み込みが無く未評価だった → 追加 (gs は `Log[x]*(1/Log[b])`)。**`Union` の `SameTest`**: gs は並べ替えも比較もせず結合した先頭の要素だけを返す (K64 の誤り。B121) → k64 は gs と同じ、v2 は真値。**`Names`**: gs は大域に作った記号だけ → k64 は gs と同じ、v2 は組み込みも含める。残った差: `Through[{f,g}[1,2]]` はどちらもエラーだが文言が違う。`FixedPoint[Cos,1.,3]` は最後の桁 (libm) |

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「複合代入が `Element`・`LINE` のフックを通らず、利用者の規則を作って以後の代入が効かなかった pysad の不具合 → 直した。**`SetElement` の戻り値**は gs は MAIN の定義の値 → 合わせた。」

元の行:

> | `element.n` | `Element` の読み出し (LENGTH・NAME・POSITION・COMPONENT・DEFAULT・VALUE・TYPE・K1・L・ANGLE・GEO・DK・EXPAND)、代入 (`=`・`+=`・`-=`・`*=`・`/=`・リスト形・`VALUE`)、`PROFILE`、`SetElement` | 読み出しは一致 | **`Element[..]+=`・`LINE[..]+=`**: 複合代入が `Element`・`LINE` のフックを通らず、利用者の規則を作って以後の代入が効かなかった pysad の不具合 → 直した。**`SetElement` の戻り値**は gs は MAIN の定義の値 → 合わせた。残った差: `PROFILE` の代入は gs も `General::invset` で pysad は別の文言のエラー、MARK の Twiss のキー (AX など) は CALC のあとの `SetElement` で gs は MAIN の 0、pysad は CALC が書き戻した値、`Element["K1","QF"]=0.9` の読み戻しが 1 ulp (0.8999999999999999) |

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「`Fill$String`・`Remove$Paren` は組み込みが無かった → 追加 |」

元の行:

> | `beamline.n` | `BeamLine` の演算、`ExtractBeamLine`、`WriteBeamLine`・`PrintBeamLine` (FFS 形式・MAIN 形式・`Name` のリスト)、`Fill$String`、`Remove$Paren` | `BeamLine`・`ExtractBeamLine`・`WriteBeamLine` は全行一致 | **`WriteBeamLine` の `Name`**: gs は `Name->"L2"` (文字列) で `Insufficient number of Names` (K64 の誤り。`{"L2"}` と書く。B122) → k64 は gs と同じ、v2 は文字列も受ける。`Fill$String`・`Remove$Paren` は組み込みが無かった → 追加 |

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「pysad は無視していた → gs と同じく DampingRate 以降の規則を返さない。」

元の行:

> | `emit.n` | `Emittance[]` (既定、`OneTurnInformation`、`Matrix`、`Orbit`、`Emittance->False`、`Region`、`SaveEMIT`) の戻りの規則の名前・並び・値、`SymplecticJ`・`SymplecticInverse`、`TouschekF`、`TouschekLifetime` の未定義のメッセージ、`Beamstrahlung*` 4 つ | 値は 1e-12 以内 (3 GeV、周長 41.8 m 余りの FODO 8 セル。Emittances の y は gs が 2e-38、pysad は 0 の雑音) | **`Emittance->False`**: pysad は無視していた → gs と同じく DampingRate 以降の規則を返さない。`Touschek` のメッセージの形を `???Touschek::undef:  …` に合わせた。残った差: `SymplecticInverse` の `-0` の表示 |

## ## 小さな台本の gs との突き合わせ (小さな台本を gs と行ごとに突き合わせた。k64 モード)

取り除いた部分: 「gs は −1、pysad は 1 → −1 に。」

元の行:

> | `dapert.n` | `DynamicApertureSurvey` (`Output`・`Axes`・`ExpandElementValues`、不正な `Axes`)、`AccelerateParticles`、`TrackParticles` | `AccelerateParticles` は 1e-14、`TrackParticles` も一致、`DynamicApertureSurvey` の Score は一致 | **`TURNS` の既定**: gs は −1、pysad は 1 → −1 に。**残った差**: `DAPWIDTH` の既定 (gs は 7 で x の振幅を上から 7 個ずつ追跡し、続けて生き残れば残りは追跡せず `*` で表示。pysad は 51 で全点を追跡。`Maximum number of particles` が gs 42、pysad 306)、`Axes->"XZ"` で z の範囲が 3 点のとき gs はエラー、pysad は受ける、不正な `Axes` のメッセージ (`DALib::wrongmode`) は未比較 |

## ### `dapert.n` ほか gs との主な差の一覧

取り除いた部分: 「| 項目 | gs | pysad (v2.1.0) |」

元の行:

> | 項目 | gs | pysad (v2.1.0) |

## ### `dapert.n` ほか gs との主な差の一覧

取り除いた部分: 「| −1 | −1 (以前は 1) |」

元の行:

> | `DynamicApertureSurvey` の `TURNS` の既定 | −1 | −1 (以前は 1) |

## ### `dapert.n` ほか gs との主な差の一覧

取り除いた部分: 「| 同じ (以前は無視) |」

元の行:

> | `Emittance->False` | DampingRate 以降の規則を返さない | 同じ (以前は無視) |

## ## SAD 自身の Packages を読むときの注意

取り除いた部分: 「が通る (v2.1.0。Windows の Tk 9.0」

元の行:

> - Python 3.14 (Tcl/Tk 9.0) でも [`test_tclbridge.py`](../pysad/tests/test_tclbridge.py)・[`test_plot.py`](../pysad/tests/test_plot.py)・[`test_gui_lang.py`](../pysad/tests/test_gui_lang.py)・[`test_corr.py`](../pysad/tests/test_corr.py) が通る (v2.1.0。Windows の Tk 9.0 は `$DefaultFixedFont` を `Courier New` にする: docs/22)。numba のキャッシュは `py -3.12 -m pysad.nbcache` (docs/12)。

## ## 未対応・限界

取り除いた部分: 「gs と全桁: v2.0.0)。」

元の行:

> - `ResponseMatrix` の TRPT (加速、`gammab` の比) は実装済み (応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。gs と全桁: v2.0.0)。
