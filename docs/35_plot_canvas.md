# 35. 描画 (ListPlot・Plot と Tk のキャンバス)

最終更新: 2026-11-16 (現状のスナップショット。経緯は `history/` にある)

`ListPlot`、`Plot`、`ParametricPlot`、`OpticsPlot`、`Show`、`Graphics` など、SAD の描画関数を、SAD 自身の `Packages/ListPlot.n`・`Plots.n`・`Canvas.n` のまま
pysad のインタプリタで動かし、tkinter のキャンバスに描く。GUI の仕組み (docs/20) の上に載っている。
コードは `pysad/tkcanvas.py`、`pysad/vecfn.py`、`pkgload.py`、`tclbridge.py`、`sadclass.py`、`sadscript.py`。試験は `test_plot.py`・`test_evalfixes.py`。

## 1. 使い方

```
FFS;
g=ListPlot[{1,4,9,16},PlotJoined->True];     (* 窓が開き、キャンバスに描かれる。戻り値は Graphics[…] *)
h=Plot[Sin[x],{x,0,6}];
```

- `$DisplayFunction` は、Windows など画面のある環境と、Linux で `DISPLAY` があるときは `CanvasDrawer` (tkinter に描く)。画面がない環境や
  `PYSAD_NO_DISPLAY=1` では `Identity` で、`Graphics` の式を返すだけ。`DisplayFunction->Identity` を引数で与えても同じ。
- `PYSAD_SAD_PACKAGES` (SAD の Packages の場所) が要る。

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
(`LatticeNames` が作る名前の表は gs と全部同じ)。試験は `test_plot.py` の `test_opticsplot_draws_the_same_geometry_as_gs`。

## 3. Tk のキャンバスの組み込み関数 (tkcanvas.py)

K64 は C (`tfTkInter_.c`) と Fortran (`tfcanvasclip.f`) で Tk の内部を直接たたく。pysad は tkinter の Tcl に `create` を送る。

| 関数 | 内容 |
|----|----|
| `TkCanvasGetTypes`・`TkCanvasPointer`・`TkCanvasLastItem`・`TkCanvasFindEnclosed` | 項目の種類名の表、キャンバスの識別子 (パス名)、最後の項目、矩形に入る項目 |
| `TkCanvasCreateItem`・`TkCanvasCreateItemDirect` | 項目を作る。`Fill->"red"` などの規則は `TkOptionLabel` で option 名にする |
| `CanvasClipLine` | 枠で折れ線を切り取り、破線に分ける (`tfcanvasclip`) |
| `CanvasSymbolDirect` | 点の記号 (○ △ □ ◇ + ×、棒) を描く |
| `TkDefineBitmap` | Windows の Tcl は名前つきビットマップを登録できないので XBM ファイルに書き、`-stipple`・`-bitmap` の名前を `@ファイル` に置き換える |
| `XListFonts`・`XGetFontPath` | gs を X で動かして得た組み込みフォント 6 個 (`{"built-ins"}`) |
| `Vector`・`ReleaseVector` (vecfn.py) | ベクトル評価。`Plots.n` の `PlotInterpolate` が使う |

未対応: `CanvasSymbol` (Tcl 文字列を返す版)、`Canvas3DClipTriangle`・`Canvas3DLightTriangle`・`Canvas3DProjection` (Canvas3D.n)。

## 4. 描画を動かすために gs に合わせた SADScript の評価

どれも gs で同じ式を実行して比べて直した (試験は `test_evalfixes.py`)。

- 連鎖した比較 `a<x<b`、`a==b==c` (以前は左から 2 つずつ評価して範囲外でも真になった)。
- `Map`・`Table`・`MapThread`・`/.` の結果の `Null[]` (空の Sequence) は取り除き、`Null[a,b]` は a, b を並べる。
- `StringPosition[s, {"x","+"}]`・`"a"|"b"` (選択肢を与えた順に、重なる一致も数える)、`Thread[{a,{1,2}}]` (リストでない要素を繰り返す)、
  `ReplacePart[l, 0, {,2}]`・`Extract[l, {,2}]` (空の位置は「全部」)。
- `With[{f[x_] := …}, …]` は関数の定義 (クラスの中ではメソッド)。`With[{x}, …]` の x は局所変数ではなく「今の値を埋め込む」指定。
- 関数の引数 `z_` が記号 `x` に束縛されているとき、本体の `Block[{z}]`・`Module[{z}]`・`With[{z = …}]` は x を局所にする (SAD は引数を本体へ文字どおり代入する)。
  そのとき引数の式 `Sin[x]` の x にも効くよう、関数の中の `With` の内側では引数の式を評価し直す。
- `Literal[x]` は x そのもの (評価も改名もしない)。`Extract[メンバー, {1}, Unevaluated$]` は評価しないメンバー自身で、`CanvasW[Create$Line] = …` は
  `CanvasWidget[Create$Line] = …` と同じ。クラスの中の `Class@メンバー` で改名後の名前も引ける。
- リテラルだけの選択肢 (`""|" "`) を引数にもつ定義は、一般の定義より先に試す。
- 定義の本体へ引数の値を埋め込むとき、定義自身のパターン変数 (`f[x_] := …` の x) は外側の同名の値に置き換えない。別の関数の呼び出しの引数は見ない。
- 評価が `FromDate[]`・`MapAt[f, a->b, {2}]` に対応 (前回の変更で GUI が壊れていたのを直した。GUI の試験は環境変数がないと skip になるので気づかなかった)。
- OpticsPlot.n を動かすために直したもの (どれも gs と比べた。試験は `test_evalfixes.py`):
  - `` `X``・`` System`X``・`` Global`X`` は、`Begin[…]` で文脈を作らない限り大域の `X` (OpticsPlot.n は `` `OpticsPlot=System`OpticsPlot=OPlot$@OpticsPlot``)。
  - Module の局所変数への `h := Sqrt[Abs[#]]&` と `h[4]` (遅延した純関数の適用)、大域の `fun := …&` の `fun[4]`。
  - `SetAttributes[f, HoldNone]` は `HoldAll` などを外す。
  - `MapAt[f, v, {,,1}]` の空の位置、`Override[…]` の結果は左辺の標準順 (文字列、続いて記号を大小無視、同じなら小文字が先)。
  - `opt:((_Rule|_RuleDelayed)...)` は要素がすべて規則のときだけ合う (以前はなんでも合った)。同じ特殊さの定義は、クラスの中でも書いた順に試す
    (親子のクラスでは子が先)。`__`・`___` は Blank を並べた定義よりわずかに一般的 (`OpticsList[__]={}` は後の定義が先)。
  - `G[x_] = Round[{H[x], x}*3]` のように右辺を評価した値のリストが x を含むとき、呼び出しで x を置き換える (Canvas.n の `CoordLimited`)。
  - `PlotRange :> With[…]` の遅延規則の右辺にも、外側の `With[{s1, leng}, …]` の値を埋め込む。`With[{c := X}, f[c] = …]` (クラスの中) は左辺にも埋め込む。
    値がクラスのメンバー (`PlotColorDefault[$DisplayFunction]`) を使う `With` の宣言は、クラスの定義時に評価せず、インスタンスのメソッドの実行時に評価する。
  - `ListPlot.n` の読み込みの途中で `Tkinter.n` が読まれても、外側の「組み込みと同名の定義を捨てる」処理が `Options[Window]` などを消さない。
  - 既定の `$HaveLaTeX` は `""` (gs も LaTeX がなければ `""`)。`PkgConfig.n` は最初に使うときに読む。
  - K64 の Tk (改造版) の text 項目だけが持つ `-usebg` は、標準の Tk にないので落とす (文字の背景は塗らない)。
- ListContourPlot・HistoPlot・FitPlot (Plots.n) を動かすために直したもの (gs と比べた。試験は `test_evalfixes.py`・`test_plot.py`):
  - `Block[{mesh}, mesh[x_] := …]` は、入るとき局所にした記号の関数定義を退避して空にし、出るとき戻す (色目盛りの ListDensityPlot が同じ Block を入れ子で使う)。
  - `?名前` (Definition) は `Hold[{名前:=名前, 定義 …}]` を返す。`Extract[…, {…,0,…}]`・`Part[…,0]` は式の頭部。`contl[k][x] =.` は添字つきの定義だけを消す。
    K64 は定義をハッシュ表の順に返す (再現できない) ので、pysad は書いた順の逆 (具体的な定義、続いてパターンの定義 `contl[_][_]`)。
    等高線の線は gs と同じ曲線だが、始点・向きが違う。
  - `MinMax[ref, ll]` は全引数。`Apply[f, l, {1}]` で f が `Null[]` を返した要素は取り除く。
  - 定義の特殊さ: パターンの中のリテラル (`{Default, y_}`) は `{x_, y_}` より特殊、完全にリテラルのリスト (`{Automatic, Automatic}`) は一部だけリテラルの
    パターン (`{Automatic, r2_}`) より特殊。`Right[sep_:0]` のように Optional を引数にもつパターンは引数の数が変わる。
  - `Module[{vars}, {init}, body]` (本体が 2 つ以上) は順に評価して最後の値。HoldRest の `par0__` は使うときに呼び出し側の環境で 1 つずつ評価する。
    `Fit[data, fun, x, par0]` の `fun`・`x`・`par0` は、関数の引数が束縛している式・記号・列に置き換える。
- Plot・ParametricPlot を x 以外の変数名 (`Plot[Sin[t], {t,0,6}]`) や関数のリスト (`Plot[{Sin[y], Cos[y]}, {y,0,6}]`) で動かすために直したもの
  (以前は `Plot[Sin[x], {x,…}]` だけで確かめていた): HoldAll の関数の `x_Symbol` は `Held[y]` に束縛されるので、`Block[{x}]` はその記号 y を局所にする。
  `x_List` は保持された (評価前の) `{…}` にも合う。`Hold[Sin[y]] /. Hold -> Null` は `Null[Sin[y]]` = Sin[y]。
  関数の中の `With[{z = …}]` は `Hold[Sin[y]]` の中の記号にも値を埋め込む。gs と比べて ListDensityPlot (全項目一致)・Plot・ParametricPlot・Plot[リスト] も全項目一致。
- GeometryPlot (`GeometryPlot.n`) を動かすために直したもの: `LINE["L$SUM", i]`・`LINE["ANGLE$SUM", i]`・`Element["ANGLE$SUM", name]` (tffs.f の tfkeyv:
  `KEY$SUM` は KEY の値で、実数のリストならその和。pysad は SPLIT を持たないので値そのもの)、記号の直後の `(` は積 (`s2(_Real|_String)`)、
  `Fill->Null` などの値が Null の規則は Tk の項目に渡さない (gs は既定の色: 線は X11 で `#000000`、Windows の Tk で `SystemButtonText`)。
  12 セルのリングの `GeometryPlot[]` は gs と項目の種類・数 (line 74、polygon 37)・座標が同じ (線の既定色だけ Tk の違い)。
- `pkgload.py` の「組み込みと同名の定義を捨てる」処理が、GUI のパッケージ (`Window`・`Canvas` …) の定義まで捨てていた (FFS ドライバ経由だと GUI が動かなかった)。直した。

## 5. 限界・未確認

- 描画の項目は ListPlot の折れ線、OpticsPlot (上記)、HistoPlot (38 項目が全部一致)、ListContourPlot (種類と数が一致、等高線の始点・向きが違う)、
  FitPlot (キャンバス 43 項目、当てはめの値は gs と 1e-8 の相対差) で比べた。`ListDensityPlot`・`Plot`・`ParametricPlot` は全項目が gs と一致。`GeometryPlot` も一致 (上記)。`BeamPlot` は比べていない。`UseBG->True` の文字の背景は描かない。
- オプション付きの `ListPlot`・`Plot` の戻りの Graphics 式 (2026-11-16 夜、`local/wake/pl2.sad`・`pl3.sad` の 10 例): PlotRange (指定値)、Scale->{Linear,Log}・{Log,Log}、PlotLabel、GridLines->None と {{1,2},{5,10}}、Frame->False、Legend、FontScale、Prolog・Epilog・Background、Plot->False、AspectRatio->1 を与えたとき、gs と pysad で**文字列として完全に一致**。窓に描いた見た目ではなく式の比較。
  差: `Options[g]` (g は Graphics 式) は gs が Graphics の全オプションを返し、pysad は `{}` を返す (`pl1.sad`。未修正)。
- 見た目 (色・フォント・大きさ) と操作 (クリック・拡大) は比べていない。文字の大きさは Tk のフォントで違う。
- `Packages/CSR.n` は読み込め、`CSR[]` のメンバーの値は gs と一致する。`MakeZL` などの流れは `CSRInit` (2 次元 Poisson) と `CSRMatrix` mode 1000 (docs/34) で動くようになり、
  3 周波数の ZL の表が gs と相対 1e-10 以内で一致する (試験は `PYSAD_SLOW=1` の `test_csr_makezl_matches_gs`、約 7 分)。`ScanTrack`・`OYStability` などは未照合。
- GUI の試験 (test_plot・test_tclbridge) は、全試験を回すとき数回に 1 回、単独では通るのに落ちることがある (Tk の窓の更新の待ち)。落ちたら単独で回して確かめる。
- K64 の `Tkinter.n` の `FontSubstitute$X11` には `Returun[…]` の綴りの誤りがあり、「フォント表が見つからない」ときの早期復帰が効かない (結果は変わらない)。

解説: docs/40_graphics_guide.md (描画)、docs/39_tcltk_guide.md (Tcl/Tk)
