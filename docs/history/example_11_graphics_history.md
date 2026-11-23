# example/11_graphics の履歴 (取り除いた記述)

md から、版・日付つきの修正・経過の記述を取り除いた。取り除いた元の行 (整理前) を、ファイルごとにそのまま残す。現在の説明は各 md を見る。

## example/11_graphics/README.md

整理前の 11 行目:

> | [fn_ListContourPlot](fn_ListContourPlot.md) | ListContourPlot (sin cos のデータ) | ○ (2026-11-17 に lower の Tcl エラーを直した。項目 68 で同数) |

整理前の 12 行目:

> | [fn_ListContourPlot_SADHelp](fn_ListContourPlot_SADHelp.md) | ListContourPlot の SADHelp の例 | ○ (2026-11-17 に修正。項目 132 で同数) |

整理前の 15 行目:

> | [fn_ParametricPlot](fn_ParametricPlot.md) | ParametricPlot (範囲に Pi あり / なし) | ○ (2026-11-17 に修正。範囲に Pi を含む A も一致) |

整理前の 18 行目:

> | [fn_PointSymbol](fn_PointSymbol.md) | PointSymbol の一覧 (未確認だった) | ○ (2026-11-17 に修正。名前の文字も項目 30 も一致) |

整理前の 21 行目:

> | [fn_OpticsPlot_mpl](fn_OpticsPlot_mpl.md) | OpticsPlot・Plot・ListPlot を matplotlib に描く (`$DisplayFunction=MatplotlibDrawer;`、v2.4.1。pysad のみ) | gs には `MatplotlibDrawer` が無い。Graphics の要素数・座標と matplotlib 側のデータが一致 (docs/42 §7c) |

整理前の 22 行目:

> | [fn_DRAW](fn_DRAW.md) | FFS の DRAW (v2.4.1、matplotlib。FODO・ASTRID2) | gs は表示環境がなく描けない。描く前のデータは gs の Graphics と一致 (docs/42 §7b) |

## example/11_graphics/fn_BeamPlot.md

整理前の 27 行目:

> **v2 の注記 (2026-11-17、v2.0.0 の再実行)**: `.pysad.txt` は既定の v2 モードの出力に作り直した。v2 は `Emittance[]`・EMIT・追跡の Jx の端の項 (B75) と `MomentumCompaction` (B77) を直したので、この例の pysad の値のうち εx (約 2.4% 大きい)・εz・減衰率・σδ・バンチ長・α・バケット高さ・放射によるチューンシフトが gs と違う (K64 のバグを直したため。`docs/47_v2_changes.md`)。以下の表と本文の「pysad」の値と「一致」の判定は、v1 と同じ `PYSAD_COMPAT=k64` での結果 (旧い出力と同じことを 2026-11-17 に確認)。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

## example/11_graphics/fn_ListContourPlot.md

整理前の 27 行目:

> - **○ 一致 (2026-11-17 に修正)**: gs・pysad とも Tcl エラーなし、`Graphics` を返し、Canvas は 68 項目 (rectangle 20、polygon 8、line 21、text 19) で同数。

整理前の 28 行目:

> - 以前の pysad は、2 つ目 (`MeshRange` つき) の描画中に `tagOrId "1" doesn't match any items` / `.1r.2b lower ...` の Tcl エラー (`Tkinter::tclerror: 1501`) を 9 回出した。原因は `TkCanvasLastItem`: K64 では Tk の `nextId-1` (`delete all` のあとも id は続く) だが、pysad は `find all` の最大値を返したので、再描画で全項目を消した直後は 0 になり、Canvas.n の `StartID=CanvasLastID+1` が古い 1 のままで `lower 85 1` (BelowFrame の矩形を枠の下へ) が存在しない項目 1 を指した。Canvas ごとに見た最大の id を覚えるようにした (試験 `test_tkcanvaslastitem_keeps_counting_after_delete_all`)。

## example/11_graphics/fn_ListContourPlot_SADHelp.md

整理前の 27 行目:

> - **○ 一致 (2026-11-17 に修正)**: gs・pysad とも `Head=Graphics Length=2`、Canvas は 132 項目 (rectangle 20、polygon 28、line 61、text 23) で同数。

整理前の 28 行目:

> - 以前の pysad は `Index out of range in Second` で失敗した。等高線の入れ子を調べる `Position[ls[[k]],_?(CP$Inside[u1,Second[#]]&),1,1]` は、level 1 の検査で頭部 `List` にも `Second[#]` をかける。gs はこの検査の中のエラー (`Second[List]`) を不一致として黙って続けるが、pysad は文ごと落とした。パターンの検査 (`_?(…)`) の中のエラーは不一致として扱うようにした (試験 `test_pattern_test_error_is_mismatch_and_listcontourplot_sadhelp_data`)。fn_ListContourPlot のデータで動いていたのは、等高線が入れ子にならず検査が頭部に届かなかったため。

## example/11_graphics/fn_ParametricPlot.md

整理前の 26 行目:

> - **○ A も B も一致 (2026-11-17 に修正)**。A (範囲に `Pi` を含む曲線のリスト) は gs・pysad とも `List`、長さ 5、B は `Graphics`、長さ 2。Canvas は gs・pysad とも 38 項目 (rectangle 1、line 21、text 16)。

整理前の 27 行目:

> - 以前の pysad は範囲に `Pi` を含むと `ParametricPlot` が評価されず式のまま残った。`ParametricPlot` は HoldAll で、`{t,0,2*Pi}` の `2*Pi` が `xmax_Real` に照合しなかったため。gs は定数式 (Pi・E・Degree と数値関数だけの式) を値で照合するので、pysad もそうした。曲線のリストは `Evaluate[ReleaseHold[a]]` を保持引数の位置で評価する必要もあった (試験 `test_parametricplot_pi_in_range_and_curve_list_like_gs`)。

## example/11_graphics/fn_PointSymbol.md

整理前の 27 行目:

> - **○ 一致 (2026-11-17 に修正)**。記号の図形 (10 個) も名前の文字も gs と同じで、Canvas の項目は gs・pysad とも 30 (rectangle 2、line 2、text 17、oval 1、polygon 8)。

整理前の 28 行目:

> - 以前は名前の文字が `KL1`、`T6` のように崩れ、項目は pysad が 29 だった。台本の `(sym[[k]]//"")[2,-2]` の `//` は StringJoin で、gs はリストを連結せず `{"CI","CL","1O"}` の形の文字列にする (StringJoin[{"a","b"}] が `{"a","b"}`)。pysad は要素を連結 (`CICL1O`) していた。gs に合わせた (試験 `test_stringjoin_does_not_concatenate_lists_like_gs`)。

## example/11_graphics/fn_DRAW.md

整理前の 3 行目:

> 種別: FFS コマンド (描画、v2.4.1)・フォルダ 11_graphics。ラティスは FODO_ring (24 セル、3 GeV、周長 105.6 m) と ASTRID2 (0.58 GeV、周長 45.7 m)。`CALC` のみ。

整理前の 7 行目:

> K64 の `DRAW` は `Packages/CANVASDRAW.n` が `OpticsPlot` を呼び、`$DisplayFunction` が Tk のときだけ描く。pysad は v2.4.1 から matplotlib で描く (`pysad/draw_mpl.py`)。

## example/11_graphics/fn_OpticsPlot_mpl.md

整理前の 3 行目:

> 種別: 関数 (描画)・フォルダ 11_graphics・v2.4.1 の新機能 (pysad のみ)。ラティスは FODO_ring (24 セル、3 GeV、周長 105.6 m、`fn_OpticsPlot` と同じ)。


## 追加 (整理の 2 回目)

### example/11_graphics/fn_ListContourPlot_SADHelp.md

整理前の 3 行目:

> 種別: 関数 (描画)・フォルダ 11_graphics・ラティスなし (pysad は 2026-11-17 に修正して gs と一致)

整理前の 16 行目:

> gs・pysad とも `Head=Graphics Length=2`、窓に等高線が出る (Canvas は 132 項目)。修正前の pysad は `???General::index: Index out of range in Second` で `g` が代入されず (`Head=Symbol`)、窓も出なかった。

整理前の 31 行目:

> `Packages/ListPlot.n` の等高線の部分 (`Second` を呼ぶ箇所) を調べていない。図 (`fn_ListContourPlot_SADHelp_plot.png`) は修正前のもので、右 (pysad) は空のまま。

### example/11_graphics/fn_OpticsPlot.md

整理前の 27 行目:

> - **△ Canvas の項目が gs 308、pysad 313 (5 個差)**。2026-11-17 に再確認した。text が gs 69、pysad 74 で、rectangle 196・line 43 は同数。文字の大きさの違いによる名前の重なり消し (`RemoveOverlap`) の結果で、docs/42 §5 と同じ原因。

### example/11_graphics/fn_ParametricPlot.md

整理前の 15 行目:

> gs・pysad とも `A Head=List Length=5`、`B Head=Graphics Length=2`。修正前の pysad は A で `???Part: index 1 out of range (length 0)` が出て `A Head=Symbol` (描かれない) だった。

整理前の 28 行目:

> - 図 (`fn_ParametricPlot_plot.png`) は 2026-11-16 のもの (A が描けなかった時点)。

### example/11_graphics/fn_PointSymbol.md

整理前の 29 行目:

> - 図 (`fn_PointSymbol_plot.png`) は 2026-11-16 のもので、右 (pysad) の名前の文字は修正前の崩れたまま。

整理前の 3 行目:

> 種別: 描画のオプション・機能 (PointSymbol)・フォルダ 11_graphics・ラティスなし。pysad で未確認だったもの

