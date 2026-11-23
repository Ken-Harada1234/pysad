# 37_plot_canvas から取り除いた履歴の記述

`docs/37_plot_canvas.md` を現在のスナップショットだけの解説にするため取り除いた文章。元の文脈 (節の見出し) つきで、文章はそのまま。

## ## 1. 使い方

取り除いた部分: 「- v2.4.1: 環境変数 `PYSAD_DISPLAY=mpl`」

元の行:

> - v2.4.1: 環境変数 `PYSAD_DISPLAY=mpl` (または台本の `$DisplayFunction=MatplotlibDrawer;`) で、同じ `Graphics[…]` を Tk でなく matplotlib に描く (`pysad/mpl_display.py`。[docs/42 §7c](42_graphics_guide.md))。

## ## 4. 描画を動かすために gs に合わせた SADScript の評価

取り除いた部分: 「どれも gs で同じ式を実行して比べて直した (試験は」

元の行:

> どれも gs で同じ式を実行して比べて直した (試験は [`test_evalfixes.py`](../pysad/tests/test_evalfixes.py))。

## ## 4. 描画を動かすために gs に合わせた SADScript の評価

取り除いた部分: 「- OpticsPlot.n を動かすために直したもの (どれも」

元の行:

> - OpticsPlot.n を動かすために直したもの (どれも gs と比べた。試験は [`test_evalfixes.py`](../pysad/tests/test_evalfixes.py)):

## ## 4. 描画を動かすために gs に合わせた SADScript の評価

取り除いた部分: 「(Plots.n) を動かすために直したもの (gs と比べた。」

元の行:

> - ListContourPlot・HistoPlot・FitPlot (Plots.n) を動かすために直したもの (gs と比べた。試験は [`test_evalfixes.py`](../pysad/tests/test_evalfixes.py)・[`test_plot.py`](../pysad/tests/test_plot.py)):

## ## 4. 描画を動かすために gs に合わせた SADScript の評価

取り除いた部分: 「で動かすために直したもの」

元の行:

> - Plot・ParametricPlot を x 以外の変数名 (`Plot[Sin[t], {t,0,6}]`) や関数のリスト (`Plot[{Sin[y], Cos[y]}, {y,0,6}]`) で動かすために直したもの
  (`x` 以外の名前でも動く): HoldAll の関数の `x_Symbol` は `Held[y]` に束縛されるので、`Block[{x}]` はその記号 y を局所にする。

## ## 4. 描画を動かすために gs に合わせた SADScript の評価

取り除いた部分: 「- GeometryPlot (`GeometryPlot.n`) を動かすために直したもの:」

元の行:

> - GeometryPlot (`GeometryPlot.n`) を動かすために直したもの: `LINE["L$SUM", i]`・`LINE["ANGLE$SUM", i]`・`Element["ANGLE$SUM", name]` (tffs.f の tfkeyv:
