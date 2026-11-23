# グラフィック描画 の検証の詳細 (gs との比較)

本文 [../21_graphics.md](../21_graphics.md) から移した、K64 (gs) との細かい比較。元の文書: 42_graphics_guide (§5、§7b の gs との照合)、37_plot_canvas (§2、§5)。見出しの節番号は元の文書のもの。

---
元: 42_graphics_guide.md

## 5. gs との比較の方法と結果

方法: gs (WSLg の X11) と pysad で同じ台本を動かし、キャンバスの全項目について、種類・座標・タグを書き出して比べた。
`Graphics[…]` の式 (`Length` と中身) も比べた。両者の Tk のフォントが違うので、文字の外接矩形だけは一致しない。

| 対象 | 結果 |
|----|----|
| `ListPlot[{1,4,9,16}, PlotJoined->True]` | 43 項目 (rectangle 1、line 22、text 16、oval 4) が種類・座標・タグとも一致。式 `Graphics[…]` も同じ |
| `Plot[Sin[x],{x,0,6}]` | `Graphics` (長さ 2)。全項目が一致 |
| `Plot[{Sin[y],Cos[y]},{y,0,6}]`、`ParametricPlot[{Sin[t],Cos[2*t]},{t,0,6}]` | gs と項目の数が一致 (line 24・text 16、line 21・text 14) |
| `ListDensityPlot` | 全項目が一致 |
| `HistoPlot[Table[Sin[i*1.7]+Cos[i*0.3],{i,1,60}]]` | 全 38 項目が一致 |
| `ListContourPlot[Table[Sin[i/3.]*Cos[j/3.],{i,0,9},{j,0,9}]]` | 種類と数 (rectangle 20、polygon 10、line 25、text 17) が一致。等高線は同じ曲線だが始点・向きが違う (K64 は `?contl` の定義の表をハッシュ表の順に返し、pysad は書いた順の逆) |
| `FitPlot[…]` (§3 の例) | キャンバス 43 項目。当てはめ値 a=1.013299、b=1.996738、ChiSquare=0.04844 は gs (a=1.0132994、b=1.9967382、ChiSquare=.0484422) と 1e-8 の相対差 |
| `FitGaussian[d, Plot->True]` | 戻り値 5 要素、キャンバス 74 項目 (rectangle 1、line 40、text 13、oval 20) が gs と同じ。値は gs と 1e-10 |
| `OpticsPlot[{{"BX","BY"},{"EX"}}]` (12 セル、3 GeV の例) | 式は `Graphics` (長さ 2) で gs と同じ構造。項目は gs が 174 (rectangle 76、line 40、text 58)、pysad は 188 (rectangle 76、line 40、text 72)。枠・矩形・格子線の数と位置は同じ。text が 14 多いのは、文字の大きさが Tk のフォントで違い、`RemoveOverlap` と目盛りの間引きの結果が変わるため。`LatticeNames` の名前の表は全部同じ |
| `GeometryPlot[]` (同じ 12 セル) | line 74、polygon 37 で gs と同じ数・座標 (線の既定色だけ Tk の違い) |

- 試験: [`pysad/tests/test_plot.py`](../../pysad/tests/test_plot.py) (8 本)、[`test_evalfixes.py`](../../pysad/tests/test_evalfixes.py)。Tk を使う試験は `PYSAD_SAD_PACKAGES` と tkinter が要る (なければ skip)。
  全試験を回すと数回に 1 回、単独では通るのに落ちることがある (Tk の窓の更新の待ち)。落ちたら `py -3.12 -m pytest pysad/tests/test_plot.py -q` で単独に回して確かめる。
- 回帰 (`run_regression.py`) は描画を含まない (描画は pytest だけ)。
- 比べていないもの: 色・フォント・大きさ・見た目、クリック・拡大などの操作、実画面のスクリーンショットの比較。値の一致は「キャンバスの項目の座標とタグ」まで。

---
元: 42_graphics_guide.md (DRAW の gs との照合)

**gs との照合** (gs は画像を出せないので、描く前のデータを比べる): `OpticsPlot[…, DisplayFunction->Identity]` が返す Graphics (pysad は Canvas.n の経路、gs は同じ式) の `Line`・`Rectangle`・`Text` を取り出し、matplotlib に渡す点列・箱・名前と比べた。
- pysad の Graphics との一致: 点列は y を単位で割ったあとまで 1e-12 以内で全点一致、ラティスの箱・名前も一致 (FODO の 6 つの組合せ: 領域・要素名の領域・周回をまたぐ領域・名前のパターン、ASTRID2 全周)。試験 `test_v25_draw_mpl.py`。
- gs (WSL) の Graphics との一致: FODO (3 GeV、6 セル、領域 3〜20、名前 `B*`) と ASTRID2 全周 (0.58 GeV、周長 45.7 m、名前 `*`、箱 139・名前 199) で、点列の最大の相対差 1.7e-13 (x は 1.4e-14 m)、箱・名前の個数と位置が一致。参照は `pysad/tests/scripts/v25_draw_gs_*.gs.txt` (台本は同名の `.sad`)。
- 例: [`example/11_graphics/fn_DRAW.md`](../../example/11_graphics/fn_DRAW.md)。

---
元: 37_plot_canvas.md

## 2. gs との一致

`ListPlot[{1,4,9,16}, PlotJoined->True]` が作るキャンバスの項目を、gs (WSLg の X11) と pysad で全部書き出して比べた。

| | gs | pysad |
|----|----|----|
| 項目の数と種類 | 43 (rectangle 1、line 22、text 16、oval 4) | 同じ |
| 座標・タグ | — | 全 43 項目が同じ (文字の外接矩形だけ Tk のフォントで数ピクセル違う) |

`ListPlot` が返す `Graphics[…]` の式も gs と同じ。`Plot[Sin[x],{x,0,6}]` は `Graphics` (長さ 2) を返し、`PlotInterpolate` の点列を作る。

`OpticsPlot[{{"BX","BY"},{"EX"}}]` (12 セルのリング、`USE RING; CELL; CALC;` のあと) は、`Graphics` (長さ 2) を返し、gs と同じ構造 (3 枠: 光学関数 2 枚と
ラティス) の式を作る。キャンバスの項目は gs が 174 個 (rectangle 76、line 40、text 58)、pysad は 188 個 (rectangle 76、line 40、text 72) で、枠 (BG)・矩形・
格子線の数と位置は同じ。text が 14 個多いのは、文字の外接矩形が Tk のフォントで違い、`RemoveOverlap` (重なる名前を消す) と目盛の間引きの結果が変わるため
(`LatticeNames` が作る名前の表は gs と全部同じ)。試験は [`test_plot.py`](../../pysad/tests/test_plot.py) の `test_opticsplot_draws_the_same_geometry_as_gs`。

---
元: 37_plot_canvas.md §5 のうち比較の記述

## 5. 限界・未確認

- 描画の項目は ListPlot の折れ線、OpticsPlot (上記)、HistoPlot (38 項目が全部一致)、ListContourPlot (種類と数が一致、等高線の始点・向きが違う)、
  FitPlot (キャンバス 43 項目、当てはめの値は gs と 1e-8 の相対差) で比べた。`ListDensityPlot`・`Plot`・`ParametricPlot` は全項目が gs と一致。`GeometryPlot` も一致 (上記)。`BeamPlot` は gs が `SIGXX=0` で楕円を描かず、比べられない (docs/42 §6)。`UseBG->True` の文字の背景は描かない。
- オプション付きの `ListPlot`・`Plot` の戻りの Graphics 式 (2026-11-16 夜、10 例): PlotRange (指定値)、Scale->{Linear,Log}・{Log,Log}、PlotLabel、GridLines->None と {{1,2},{5,10}}、Frame->False、Legend、FontScale、Prolog・Epilog・Background、Plot->False、AspectRatio->1 を与えたとき、gs と pysad で**文字列として完全に一致**。窓に描いた見た目ではなく式の比較。
  差: `Options[g]` (g は Graphics 式) は gs が Graphics の全オプションを返し、pysad は `{}` を返す (未修正)。
