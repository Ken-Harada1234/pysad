# 11_graphics: 描画関数

Plot 系の関数の例。グラフは Tk の窓に出るので標準出力には出ない (`.gs.txt`・`.pysad.txt` には起動表示・台本・`Print`・エラーだけが残る)。図 (`NAME_plot.png`) は、Canvas の項目を Tcl で書き出し、matplotlib で再描画したもの (左 gs、右 pysad)。作り方と、gs の窓が撮れない理由は、各 `.md` の「図」と [`example/tools/canvas_render.py`](../tools/canvas_render.py) を見る。

一致の記号: ○ 一致 / △ 差あり (md に理由) / ×gs 失敗 / ×pysad 未対応・失敗

| 例 (NAME.sad・.gs.txt・.pysad.txt・.md) | 内容 | gs と pysad |
|---|---|---|
| [fn_Plot](fn_Plot.md) | Plot (関数 2 本、Identity) | ○ |
| [fn_ListPlot](fn_ListPlot.md) | ListPlot (折れ線・点・誤差棒・凡例) | ○ |
| [fn_ListContourPlot](fn_ListContourPlot.md) | ListContourPlot (sin cos のデータ) | ○ (2026-11-17 に lower の Tcl エラーを直した。項目 68 で同数) |
| [fn_ListContourPlot_SADHelp](fn_ListContourPlot_SADHelp.md) | ListContourPlot の SADHelp の例 | ○ (2026-11-17 に修正。項目 132 で同数) |
| [fn_ListDensityPlot](fn_ListDensityPlot.md) | ListDensityPlot | ○ |
| [fn_HistoPlot](fn_HistoPlot.md) | HistoPlot | ○ |
| [fn_ParametricPlot](fn_ParametricPlot.md) | ParametricPlot (範囲に Pi あり / なし) | ○ (2026-11-17 に修正。範囲に Pi を含む A も一致) |
| [fn_FitPlot](fn_FitPlot.md) | FitPlot | ○ |
| [fn_ColumnPlot](fn_ColumnPlot.md) | ColumnPlot (pysad で未確認だった) | ○ (動く) |
| [fn_PointSymbol](fn_PointSymbol.md) | PointSymbol の一覧 (未確認だった) | ○ (2026-11-17 に修正。名前の文字も項目 30 も一致) |
| [fn_LaTeX](fn_LaTeX.md) | LaTeX 文字列 (未確認だった) | ○ (両方 LaTeX 無効の状態で一致。LaTeX 有効時は未確認) |
| [fn_OpticsPlot](fn_OpticsPlot.md) | OpticsPlot (BX BY EX と配置) | △ (項目数 308 対 313、文字の重なり消し) |
| [fn_OpticsPlot_mpl](fn_OpticsPlot_mpl.md) | OpticsPlot・Plot・ListPlot を matplotlib に描く (`$DisplayFunction=MatplotlibDrawer;`、v2.4.1。pysad のみ) | gs には `MatplotlibDrawer` が無い。Graphics の要素数・座標と matplotlib 側のデータが一致 (docs/42 §7c) |
| [fn_DRAW](fn_DRAW.md) | FFS の DRAW (v2.4.1、matplotlib。FODO・ASTRID2) | gs は表示環境がなく描けない。描く前のデータは gs の Graphics と一致 (docs/42 §7b) |
| [fn_GeometryPlot](fn_GeometryPlot.md) | GeometryPlot | ○ |
| [fn_BeamPlot](fn_BeamPlot.md) | BeamPlot (未確認だった) | △ (gs は楕円なし、pysad は描く) |

gs の窓の画面そのものは撮れていない (WSL に撮る道具が無い)。gs が窓を開いたことは `winfo ismapped`=1 で確かめた。
