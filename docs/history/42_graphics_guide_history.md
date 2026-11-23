# 42_graphics_guide の履歴 (整理で本文から取り除いた記述)

本文 (`docs/42_graphics_guide.md`) を版のないスナップショットにするとき、取り除いた版・経過の記述を、元の文脈つきで残したもの。元の行をそのまま引用する (引用した行のうち、履歴にあたるのは版・経過を述べた部分)。

## 元 36 行 (節: 2. 仕組み)

>   環境変数 `PYSAD_DISPLAY=mpl` (v2.4.1) で `MatplotlibDrawer` (Graphics を matplotlib に描く。§7c)、`tk`・`none` で `CanvasDrawer`・`Identity` に固定。

## 元 170 行 (節: 7. 未対応・限界・未確認)

>   `OpticsPlot` に自前のデータ (`{{{line1},{line2},…, 選択肢}}`。x は要素番号) を渡せる (定義の順に合わせた規則の優先順位と `MapThread` の非リスト要素の繰り返しを直した)。

## 元 173 行 (節: 7. 未対応・限界・未確認)

> - `Options[g]` の差: gs は Graphics の全オプションを返し、pysad は `{}` を返す (未修正)。

## 元 179 行 (節: 7b. DRAW (matplotlib) — v2.4.1)

> ## 7b. DRAW (matplotlib) — v2.4.1

## 元 211 行 (節: 7c. MatplotlibDrawer: Graphics を matplotlib に描く `$DisplayFunction` — v2.4.1)

> ## 7c. MatplotlibDrawer: Graphics を matplotlib に描く `$DisplayFunction` — v2.4.1

## 元 215 行 (節: 7c. MatplotlibDrawer: Graphics を matplotlib に描く `$DisplayFunction` — v2.4.1)

> **仕組み** `ListPlot.n` の `Show1$` が `$DisplayFunction` を `df[opt, list, init, w, rect, last, ini]` の形で呼ぶ (K64 では `CanvasDrawer`)。pysad は同じ引数を受ける組み込み関数 `MatplotlibDrawer` を足した。`opt` の `PlotRange`・`DataRange`・`Scale`・`GridLines`・`FrameTicks` などは `Show1$` が決めた値。`Rectangle[min, max, Graphics[…]]` (OpticsPlot の枠) は、`Canvas.n` と同じく `Show1$[False, {min,max}, …]` を呼んで入れ子の枠を作る (`Position->Right[]`・`Left[]`・`Above[]`・`Below[]` も `RectPosition` どおり)。枠の位置は `Canvas.n` の `MakeWindow` の移植 (キャンバス 13 × 10 インチ、窓の原点 {2, 1.8}、幅 9.6、縦横比 `GoldenRatio`、`AspectRatio`・`PlotRegion`・`rect` に従う)、軸の範囲は `FrameRange` の移植 (`PlotRange->Automatic` は `DataRange` の外に 7% の余白、Log は log10 で 7%)。

## 元 221 行 (節: 7c. MatplotlibDrawer: Graphics を matplotlib に描く `$DisplayFunction` — v2.4.1)

> | 環境変数 `PYSAD_DISPLAY=mpl` (または `matplotlib`) | `$DisplayFunction` を `MatplotlibDrawer` にする (`tk` = `CanvasDrawer`、`none` = `Identity`。未設定は §2 の従来どおり) |

## 元 254 行 (節: 8. 関連文書)

> - [`docs/37_plot_canvas.md`](37_plot_canvas.md): 描画の gs との一致、Tk の組み込み関数の表、評価器の直した点の全部

