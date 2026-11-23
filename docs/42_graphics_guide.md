# 42. グラフィック描画の解説 (Plot・ListPlot・OpticsPlot ほかと Tk のキャンバス)

最終更新: 2026-11-22

SAD の描画関数 (`Plot`、`ListPlot`、`ListContourPlot`、`HistoPlot`、`FitPlot`、`OpticsPlot`、`LatticePlot`、`BeamPlot`、`GeometryPlot`、`FitGaussian` の `Plot->True`、
`Canvas` 系) を、pysad がどう動かし、どこまで gs (SAD K64) と一致するかを、初めて読む人向けにまとめる。
要点の一覧と SADScript の評価の直し方の細目は [`docs/37_plot_canvas.md`](37_plot_canvas.md)、GUI (Tk) の土台は [`docs/22_gui_tclbridge.md`](22_gui_tclbridge.md)。この文書は両者の案内と、使い方の例、限界をまとめたもの。

## 1. 何ができるか

- SAD 自身の `.n` (`Packages/ListPlot.n`、`Plots.n`、`Canvas.n`、`OpticsPlot.n`、`GeometryPlot.n`、`Fit.n`) を **書き直さず**、pysad のインタプリタで SADScript のまま実行する。
  描画の式 (`Graphics[…]`) は SAD と同じ構造で作られ、最後に `$DisplayFunction` が `CanvasDrawer` なら、`Canvas.n` が Tk の `create` を呼んで tkinter のキャンバスに描く。
- 描画の物理や座標の計算 (目盛りの選び方、格子線、`PlotInterpolate` の点列、等高線、`RemoveOverlap`) は SAD の `.n` がするので、pysad はそれが動く評価器と、
  K64 が C・Fortran で持つ Tk の組み込み関数 (`TkCanvas*`、`CanvasClipLine`、`CanvasSymbolDirect` など) だけを Python で書いた ([`pysad/tkcanvas.py`](../pysad/tkcanvas.py))。
- 対応する関数: `Plot`・`ParametricPlot`・`ListPlot`・`ListDensityPlot`・`ListContourPlot`・`HistoPlot`・`FitPlot`・`OpticsPlot`・`LatticePlot` (OpticsPlot の一部)・`GeometryPlot`・
  `BeamPlot` (動くが gs と比べられない: §6)・`FitGaussian[…, Plot->True]`・`Show`・`Graphics`。

## 2. 仕組み

```
SADScript (Plot[Sin[x],{x,0,6}] など)
   │  pysad のインタプリタ (sadscript.py ほか) が Packages/Plots.n・ListPlot.n・Canvas.n を読んで実行
   ▼
Graphics[{Line[…], Text[…], …}, オプション]      ← DisplayFunction->Identity ならここで終わり (戻り値)
   │  DisplayFunction->CanvasDrawer
   ▼
Canvas.n の描画 (座標変換、枠、目盛り、格子、文字)
   │  TkCanvasCreateItem[…] など (pysad/tkcanvas.py)
   ▼
tkinter の Tcl に  `.c create line …` を送る → キャンバス (tclbridge.py の Tk の窓)
```

- 読む場所: `PYSAD_SAD_PACKAGES` (SAD の `Packages` のフォルダ) が要る。`pkgload.py` が `init.n` と同じ自動読み込みの表を持ち、`Plot`・`ListPlot` を最初に使うときにそのファイルを読む。
- `$DisplayFunction`: Windows、または Linux で `DISPLAY` があれば `CanvasDrawer` (tkinter に描く)。画面がない環境と `PYSAD_NO_DISPLAY=1` では `Identity`
  (`Graphics[…]` の式を返すだけ)。引数の `DisplayFunction->Identity` でも同じ (`pkgload.py`)。
  環境変数 `PYSAD_DISPLAY=mpl` (v2.4.1) で `MatplotlibDrawer` (Graphics を matplotlib に描く。§7c)、`tk`・`none` で `CanvasDrawer`・`Identity` に固定。
- `PYSAD_TK_NOLOOP=1`: `TkEventLoop[]` が 1 回 `update` して戻る。自動試験用 (`tclbridge.py`)。普段は窓を閉じるまでイベントを回す。
- `tkcanvas.py` の組み込み関数 (K64 の `tfTkInter_.c`・`tfcanvasclip.f` に対応。表は `docs/37` §3): `TkCanvasGetTypes`・`TkCanvasPointer`・`TkCanvasCreateItem(Direct)`・
  `TkCanvasLastItem`・`TkCanvasFindEnclosed`・`CanvasClipLine` (枠で折れ線を切る・破線に分ける)・`CanvasSymbolDirect` (点の記号)・`TkDefineBitmap`・`XListFonts`・`XGetFontPath`。
  K64 は Tk の内部構造を直接たたくが、pysad は Tcl の `create` コマンドに文字列で渡す。Windows の Tcl は名前つきビットマップを登録できないので、
  XBM ファイルに書いて `-stipple`・`-bitmap` の値を `@ファイル` に置き換える。
- `Vector`・`ReleaseVector` (`vecfn.py`): `Plots.n` の `PlotInterpolate` が使うベクトル評価。`FitGaussian` の `Plot->True` が `Plot[InverseErf[x]…]` を使うので、
  `InverseErf` などの特殊関数も Vector に入れた (`docs/33`)。
- `fitroot.py` の Plot 関係: `FitGaussian[d, opt…]` は Fit.n と同じく、`Plot->True` のとき `{sig, 平均, {…}, ChiSquare, Show[Plot(当てはめの曲線), ListPlot(ビンのデータ), opt]}` の 5 要素のリストを返す
  (内部で `Plot` と `ListPlot[…, DisplayFunction->Identity]` を評価して `Show` に渡す)。`FitPlot` は `Fit[data, fun, x, par0]` を呼ぶので、`fun`・`x`・`par0` を
  関数の引数が束縛している式に置き換える処理が `fitroot.py` にある。
- 評価器の側で描画のために gs に合わせた細目 (連鎖した比較、`With`・`Block` の局所化、`Literal`、`Override` の順、定義の特殊さの順など) は `docs/37` §4 に全部ある。

## 3. 使い方の例 (動く SADScript)

画面のある環境では、窓が開いて描かれる。`PYSAD_SAD_PACKAGES` を設定して `FFS;` から。

```
FFS;
g=ListPlot[{1,4,9,16},PlotJoined->True];       (* 折れ線。戻り値は Graphics[…] (長さ 2) *)
h=Plot[Sin[x],{x,0,6}];                          (* Plot[f,{x,xmin,xmax}]。変数名は x でなくてよい *)
k=Plot[{Sin[y],Cos[y]},{y,0,6}];                 (* 関数のリスト *)
p=ParametricPlot[{Sin[t],Cos[2*t]},{t,0,6}];
Hs=HistoPlot[Table[Sin[i*1.7]+Cos[i*0.3],{i,1,60}]];
C=ListContourPlot[Table[Sin[i/3.]*Cos[j/3.],{i,0,9},{j,0,9}]];
F=FitPlot[Table[{i,2*i+1+0.1*Sin[i*3]},{i,1,10}],a+b*x,x,{a,1},{b,1}];   (* a,b と ChiSquare が g[[1]] に入る *)
d=Table[Sqrt[2]*InverseErf[2*(i-0.5)/400-1]*0.8+0.3+0.05*Sin[i*1.3],{i,400}];
r=FitGaussian[d,Plot->True];                     (* {sig, 平均, {…}, ChiSquare, Show[…]} *)
```

リングの光学関数・配置図 (`USE`、`CELL`、`CALC` のあと):

```
MOMENTUM = 3GEV;
MARK ST=();
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1);
FFS;
USE RING; CELL; CALC;
g=OpticsPlot[{{"BX","BY"},{"EX"}}];     (* 3 枠: BX・BY、EX、ラティス *)
q=GeometryPlot[];                        (* 平面の配置図 *)
```

画面なしで式だけを得る (自動処理・試験向け):

```
g=ListPlot[{1,4,9,16},PlotJoined->True,DisplayFunction->Identity];
Print[Head[g]," ",Length[g]];            (* Graphics 2 *)
```

BeamPlot (ビームの断面の楕円):

```
MOMENTUM = 3GEV;
(* …上と同じラティスに CAVI RFC = (VOLT=2e6 HARM=64); LINE RING = (ST 12*CELL1 RFC); … *)
FFS;
USE RING; CELL; RFSW; RAD; CODPLOT; CALC; EMIT;
g=BeamPlot["QF.5",{"X","Y"},DisplayFunction->Identity,SizeFunction->"SIG"];
```

試験の元の台本は [`pysad/tests/test_plot.py`](../pysad/tests/test_plot.py) (8 本)。
キャンバスの項目を書き出す方法 (gs と比べるのに使った) は、[`test_plot.py`](../pysad/tests/test_plot.py) の `_DUMP`: `TclEval` で全 Canvas の `find all`・`type`・`coords`・`gettags` を列挙してファイルに書く。

## 4. オプションの対応

オプションは `ListPlot.n` の定義 (`LP$def`、`Graphics$Options`) のまま pysad が評価する。pysad 側で特別に処理するものはなく、`Override[opt, 既定]` の規則どおりに効く。
次は `ListPlot.n` の既定値 (読み取った値) と、試験で動作を確認したか。

| オプション | 既定 (`ListPlot.n`) | 確認 |
|----|----|----|
| `PlotJoined` | `False` (折れ線にするとき `True`) | `ListPlot[…, PlotJoined->True]` で gs と全 43 項目一致 |
| `PlotRange` | `{Automatic,Automatic}` | `OpticsPlot` が `PlotRange :> With[…]` を使って動く。指定値を与えたときの戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較) |
| `GridLines` | `{Automatic,Automatic}` | 格子線の位置は gs と一致 (`G$X1`・`G$Y1` のタグ)。`GridLines->None` と `{{1,2},{5,10}}` を与えたときの戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較) |
| `AspectRatio` | `Automatic` (`BeamPlot` は `1`) | `AspectRatio->1` の戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較)。`BeamPlot` は式の生成だけ。窓の見た目は未確認 |
| `Scale` | `{Linear,Linear}` | `Scale->{Linear,Linear}`以外 (`{Linear,Log}`・`{Log,Log}`) の戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較)。窓の見た目は未確認 |
| `Frame`・`FrameLabel`・`FrameTicks`・`TickSize` | `True`・空・`{True,True,Ticks,Ticks}`・`Automatic` | `Frame->False` の戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較)。他は既定のみ |
| `PointSymbol`・`PointSize`・`PointColor`・`PlotColor`・`Dashing`・`Thickness` | 既定は表の色・記号・破線の巡回 | ListPlot の 4 個の `oval` (記号 `CI`)。他の記号は `CanvasSymbolDirect` の実装 (○ △ □ ◇ + ×、棒) があるが gs と比べていない |
| `PlotLabel`・`Legend`・`Prolog`・`Epilog`・`Background`・`FontScale` ほか | 既定 | `PlotLabel`・`Legend`・`Prolog`・`Epilog`・`Background`・`FontScale` を与えたときの戻りの Graphics 式は gs と文字列として完全に一致 (2026-11-16。窓の見た目ではなく式の比較)。窓の見た目は未確認 |
| `DisplayFunction` | `$DisplayFunction` (遅延) | `Identity` と `CanvasDrawer` の両方を試験 |
| `RemoveOverlap` | `""` | OpticsPlot が名前の重なりを消す。文字の大きさに依存 (§5) |
| `Plot->False` | (`Plot->True` が既定) | 戻りの Graphics 式が gs と一致 (同上) |
| `Orbit`・`SizeFunction` (BeamPlot)、`DataRange` | `True`・`"SIZE"`、`Default` | §6 |

`Fill->Null` のように値が `Null` の規則は Tk の項目に渡さない (gs は既定の色を使う。線は X11 で `#000000`、Windows の Tk で `SystemButtonText`)。

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

- 試験: [`pysad/tests/test_plot.py`](../pysad/tests/test_plot.py) (8 本)、[`test_evalfixes.py`](../pysad/tests/test_evalfixes.py)。Tk を使う試験は `PYSAD_SAD_PACKAGES` と tkinter が要る (なければ skip)。
  全試験を回すと数回に 1 回、単独では通るのに落ちることがある (Tk の窓の更新の待ち)。落ちたら `py -3.12 -m pytest pysad/tests/test_plot.py -q` で単独に回して確かめる。
- 回帰 (`run_regression.py`) は描画を含まない (描画は pytest だけ)。
- 比べていないもの: 色・フォント・大きさ・見た目、クリック・拡大などの操作、実画面のスクリーンショットの比較。値の一致は「キャンバスの項目の座標とタグ」まで。

## 6. K64 の癖

- **BeamPlot は gs の `SIG*`・`SIZE*` が 0 で、比較できない**。`CALC` のみ、`EMIT`、`Emittance[]`、`EMITX=…; EMITY=…; SIGE=…` を与えた場合の
  どれでも、gs の `LINE["SIZEXX","QF.5"]`・`LINE["SIGXX","QF.5"]` などが 0 を返す。このため `BeamPlot` の楕円の半径が 0 になり、
  gs が返す `Graphics` の線は `{座標, NaN}` の列 (`EMIT` や `Emittance[]` のあとは 33 点が `{-10.864…, NaN}`) で、描いても意味のある形にならない。
  pysad は、`EMIT` のあとは値を返す (`SIGXX` = 2.03e-6、`SIGYY` = 7.1e-9。gs は 0) ので、**結果が gs と食い違う**。どちらが正しいかは確認していない
  (K64 の `LINE["SIG…"]` の仕様か、この gs のビルド 2026-10-19 の不具合かは未確認)。
  また `EMIT` なしでは、pysad は `???Set: list assignment length mismatch` を出して `BeamPlot` が戻らない (gs は NaN の線を含む `Graphics` を返す)。この差は直していない。
- gs は X11 の `-usebg` を持つ改造版の Tk を使う。標準の Tk にない `-usebg` は pysad が落とす (文字の背景は塗らない)。
- `Tkinter.n` の `FontSubstitute$X11` に `Returun[…]` の綴りの誤りがあり、フォント表がないときの早期復帰が効かない (結果は変わらない)。
- 線の既定色は、gs の X11 では `#000000`、Windows の Tk では `SystemButtonText`。比較では「同じ」とみなした。
- 既定の `$HaveLaTeX` は `""` (LaTeX がなければ gs も同じ)。
- `?名前` (Definition) の定義の順は K64 がハッシュ表の順で、再現できない (等高線の始点・向きの差の原因)。

## 7. 未対応・限界・未確認

- 未対応: `CanvasSymbol` (Tcl 文字列を返す版。`CanvasSymbolDirect` は実装)、`Canvas3DClipTriangle`・`Canvas3DLightTriangle`・`Canvas3DProjection` (`Canvas3D.n`)、
  `UseBG->True` の文字の背景 (描かない)。
- 未確認: `BeamPlot` の描画 (上記)、`LatticePlot` を単独で呼ぶ使い方 (`OpticsPlot` の 3 枠目として動くことだけ確認)、上記の各オプションを窓に描いたときの見た目 (式は gs と一致するが、描画は未確認)、
  `CanvasCM.n`・`CANVASDRAW.n`、`TabFrame.n` の描画、色・フォント・大きさの見た目、マウス操作。
- gs の実動作に合わせた細目 (2026-11-17): `ParametricPlot[…,{t,0,2*Pi}]` のように範囲に Pi などの定数式を書いても `x_Real` に照合する (HoldAll の引数の定数式 (Pi・E・Degree と数値関数) は値で照合する)。
  曲線のリスト `ParametricPlot[{{…},{…}},…]` は List (曲線数+1 要素) を返す (保持引数位置の `Evaluate[…]` を評価)。`ListContourPlot` (MeshRange つきの SADHelp の例を含む) は Graphics を返し、窓への描画で Tcl エラーが出ない
  (`_?(…)` の検査の中のエラーは不一致として扱う)。`PointSymbol` の名前の文字は、`{"CI","CL"}//""` がリストを `{"CI","CL"}` の形の文字列にする (gs と同じ) ことで一致し、Canvas の項目は gs と同じ 30。
  `OpticsPlot` に自前のデータ (`{{{line1},{line2},…, 選択肢}}`。x は要素番号) を渡せる (定義の順に合わせた規則の優先順位と `MapThread` の非リスト要素の繰り返しを直した)。
  `Canvas$Widget = c1` を使う台本の窓は gs と同じ 114 項目 (再描画後 115)。
- `OpticsPlot` の text 項目の数の差 (FODO 24 セルの例で gs 69、pysad 74) は、文字の大きさ (Tk のフォントの違い) による `RemoveOverlap` の結果。直せないので未対応として残す。
- `Options[g]` の差: gs は Graphics の全オプションを返し、pysad は `{}` を返す (未修正)。
- 文字の数と外接矩形は Tk のフォントで決まるので、`OpticsPlot` の text 項目の数は gs と一致しない (pysad は +5〜14)。枠・矩形・線は一致する。
- 比較したのは 3 GeV の 12 セルの小さなリング (E1=E2=θ/2 の BEND を持つ FODO) と、数値の列だけ。実ラティスの `OpticsPlot` は比べていない。
- 描画の窓は tkinter。Windows のほかは試していない (gs の比較は WSLg の X11)。
- GUI の試験は全体実行で不安定なことがある (§5)。環境変数がないと skip になるので、評価器を変えたら `PYSAD_SAD_PACKAGES` を付けて試験を回すこと。

## 7b. DRAW (matplotlib) — v2.4.1

FFS の `DRAW` コマンドを、Tk のキャンバスを使わず matplotlib で描く ([`pysad/draw_mpl.py`](../pysad/draw_mpl.py)。`pip install matplotlib`。無くても DRAW 以外は動き、DRAW は警告して何もしない)。

**K64 の DRAW**: `src/tffsa.f` が `System`CANVASDRAW[]` を評価し、実体は `Packages/CANVASDRAW.n` が `OpticsPlot[list, Region->{r1,r2}, Lattice->…, Names->…, GridLines->{None,Automatic}, InfoLabel->True, PlotLabel->TITLE]` を呼ぶ。`$DisplayFunction` が `CanvasDrawer` でなければ何もしない。gs は表示環境がないと何も出さず、`FFS["draw ...",12]` は内部エラー (docs/46 §5)。

**構文** `DRAW [r1 [r2]] 関数 [関数 …] [& 関数 …] … [LAT|LATTICE|要素名];`
- r1・r2: 整数か要素名 (`LINE["POSITION"]`)。読めなければ r1=1 でその語は関数として読み、r2 の失敗は `LINE["LENGTH"]`。複数の要素に合う名前は最初の位置。
- 関数: `Twiss$Label` が空でない名前 (BX・BY・EX・EY・NX・NY・EPX・PEX・PEPX・AX・BZ・ZX・DETR・R1〜R4、R・D を付けた参照光学 (RBX・DBX …)、BMAGX、BZS・BZT・GX …、SIGX・SIZEXX など)。BX・BY・BZ は √β を描く。`&` か `&` を語に付けて (`BX&BY`) 図 (枠) を分け、縦に積んで s の軸を共有する。1 つの枠に単位の違う関数があるときは、先頭と同じ単位だけを描く (OpticsPlot::units)。
- `LAT`・`LATTICE`: ラティス図を一番下に付ける (枠の高さの 0.23 倍)。要素名 (パターン可、`b*`) を書くと、ラティス図 + その名前だけを表示 (`Names->w`)。LAT のときの名前は全要素 (`Names->"*"`) で、重なる文字は消す。それ以外の語は DRAW が読まずに次のコマンドへ返す (`DRAW BX CALC;` は CALC を実行する)。
- 単位は `Optics$Unit` と同じに振幅で選ぶ (EX は mm、位相 NX は 2π を単位に、正の関数 (BX など) は 0 から)。`TITLE` を図の題にする。
- ラティス図は OpticsPlot.n の `LatticePlot` の箱 (四極 0.8 limegreen、偏向 0.4 gold (ANGLE=0 は tomato)、6 極以上 0.6 skyblue、ドリフト 0.08 で塗りなし、MULT は TRPT で 0.5/0.8、RF 0.3 …)。これは **ラティス配置図 (長直線部を細線 0.55 pt、その他の基準軌道 1.8 pt、偏向磁石 3 pt で描く図) とは別物**で、OpticsPlot の枠の下の簡略図。
- InfoLabel: 図をクリックすると、その s の要素の入口の名前・位置・K・β・α・ν・η・η'・R を注釈に出す (OP$InfoLabel。右クリックで消す)。OpticsPlot.n の R3,4 の行は R4 を出すべきところに R2 を出している (誤記と見て R4)。Python API の `tunes=True` でチューンと周長を図に添える (SAD にない追加)。

**出力** 保存先があれば保存して表示しない: 環境変数 `PYSAD_DRAW_FILE` (拡張子 png・pdf・svg。同じ DRAW が 2 回目以降は `名前_2.png` のように番号を付け、`{n}` があれば番号に置き換える)、または SAD 側の `Draw$Option={File->"x.png"}`。保存先がなければ、対話の backend のときだけ `plt.show` (Agg など非対話では何も出さない)。`$DisplayFunction` が Identity (headless の既定) で保存先もなければ K64 と同じく何もしない。`Draw$Option` は CANVASDRAW.n と同じく最優先 (`Region`・`Lattice`・`Names`・`GridLines`・`InfoLabel`・`PlotLabel` と、pysad の拡張 `File`・`Show`・`Tunes`)。`$DisplayFunction` が CanvasDrawer か `MatplotlibDrawer` なら matplotlib で描く。対話表示した窓は、スクリプトの終了時に閉じるまで残す (`PYSAD_DRAW_BLOCK=0` で残さない)。

```python
from pysad.ffscmd import run_sad_string
from pysad import draw_mpl
ffs, I = run_sad_string(script)                     # USE・CALC まで
fig = draw_mpl.draw(ffs, "BX BY & EX LAT", region=(1, 40), save="optics.png")   # Figure を返す
```

**CALC していないとき**: pysad は Twiss[] と同じく先に光学を計算する。gs の `OpticsPlot` は計算していない光学の値 (0) をそのまま描く (`Twiss["BX",5]` が 0) ので、ここは gs と違う (計算してから描く方が意図に合うと判断)。USE していなければ ` ???-DRAW-Missing beamline in USE.`。

**gs との照合** (gs は画像を出せないので、描く前のデータを比べる): `OpticsPlot[…, DisplayFunction->Identity]` が返す Graphics (pysad は Canvas.n の経路、gs は同じ式) の `Line`・`Rectangle`・`Text` を取り出し、matplotlib に渡す点列・箱・名前と比べた。
- pysad の Graphics との一致: 点列は y を単位で割ったあとまで 1e-12 以内で全点一致、ラティスの箱・名前も一致 (FODO の 6 つの組合せ: 領域・要素名の領域・周回をまたぐ領域・名前のパターン、ASTRID2 全周)。試験 `test_v25_draw_mpl.py`。
- gs (WSL) の Graphics との一致: FODO (3 GeV、6 セル、領域 3〜20、名前 `B*`) と ASTRID2 全周 (0.58 GeV、周長 45.7 m、名前 `*`、箱 139・名前 199) で、点列の最大の相対差 1.7e-13 (x は 1.4e-14 m)、箱・名前の個数と位置が一致。参照は `pysad/tests/scripts/v25_draw_gs_*.gs.txt` (台本は同名の `.sad`)。
- 例: [`example/11_graphics/fn_DRAW.md`](../example/11_graphics/fn_DRAW.md)。

**未確認・限界**: 色・線の太さ・文字の大きさの見た目は Tk の窓と比べていない (比べたのはデータ)。文字の重なり消し (`RemoveOverlap`) は Tk のフォントの代わりに s 軸上の間隔で決めるので、gs の窓で残る名前とは違う。`PEX` などは pysad の `Twiss["PEX"]` (現状は EX と同じ値) をそのまま描く。SIG・SIZE の関数は EMIT 後の値。SADScript の `OpticsPlot[…]` の Graphics 自体を matplotlib に流すのは §7c の `MatplotlibDrawer` (DRAW は `draw_mpl` が自前でデータを組む)。対話表示の窓・クリックの InfoLabel は Agg では動かせず、手では確認していない (イベント処理の関数 `element_info` の文字だけ試験)。ラティス図の Element$Height は K64 の型の番号に従い、SOL (20) は 0.08 の細い箱になる。

## 7c. MatplotlibDrawer: Graphics を matplotlib に描く `$DisplayFunction` — v2.4.1

`Plot`・`ListPlot`・`HistoPlot`・`ListContourPlot`・`OpticsPlot`・`GeometryPlot` などが作る `Graphics[…]` を、Tk のキャンバスでなく matplotlib に描く出口 ([`pysad/mpl_display.py`](../pysad/mpl_display.py)。matplotlib は遅延 import で、無ければ警告するだけで評価は続く)。§7b の `DRAW` は OpticsPlot 相当のデータを自前で組むが、こちらは SAD 自身の `.n` (ListPlot.n・Plots.n・OpticsPlot.n・GeometryPlot.n) が作った `Graphics[…]` をそのまま描く。評価器・`.n` は変えていない。

**仕組み** `ListPlot.n` の `Show1$` が `$DisplayFunction` を `df[opt, list, init, w, rect, last, ini]` の形で呼ぶ (K64 では `CanvasDrawer`)。pysad は同じ引数を受ける組み込み関数 `MatplotlibDrawer` を足した。`opt` の `PlotRange`・`DataRange`・`Scale`・`GridLines`・`FrameTicks` などは `Show1$` が決めた値。`Rectangle[min, max, Graphics[…]]` (OpticsPlot の枠) は、`Canvas.n` と同じく `Show1$[False, {min,max}, …]` を呼んで入れ子の枠を作る (`Position->Right[]`・`Left[]`・`Above[]`・`Below[]` も `RectPosition` どおり)。枠の位置は `Canvas.n` の `MakeWindow` の移植 (キャンバス 13 × 10 インチ、窓の原点 {2, 1.8}、幅 9.6、縦横比 `GoldenRatio`、`AspectRatio`・`PlotRegion`・`rect` に従う)、軸の範囲は `FrameRange` の移植 (`PlotRange->Automatic` は `DataRange` の外に 7% の余白、Log は log10 で 7%)。

**切替**

| 方法 | 内容 |
|----|----|
| 環境変数 `PYSAD_DISPLAY=mpl` (または `matplotlib`) | `$DisplayFunction` を `MatplotlibDrawer` にする (`tk` = `CanvasDrawer`、`none` = `Identity`。未設定は §2 の従来どおり) |
| SADScript で `$DisplayFunction=MatplotlibDrawer;` | 環境変数なしでも切り替わる |
| `Plot[…, DisplayFunction->MatplotlibDrawer]`、`Show[g, DisplayFunction->MatplotlibDrawer]` | その 1 回だけ。`DisplayFunction->Identity` で作った `g` をあとから描くのにも使う |
| Python `mpl_display.render(I, g, save="x.png")` | Graphics の式 `g` を描いて Figure を返す (`mpl_display.last_figure()` は最後に描いた Figure) |

**出力** 環境変数 `PYSAD_DISPLAY_FILE` (拡張子 png・pdf・svg。同じ名前の 2 回目以降は `名前_2.png` のように番号を付け、`{n}` があれば番号に置き換える。DRAW の `PYSAD_DRAW_FILE` と同じ流儀で、番号の数え方も `draw_mpl` と共通)、または SAD 側の大域変数 `MatplotlibDrawer$Option={File->"x.png", Show->True, DPI->120}` (最優先。保存先のフォルダは作る。dpi の既定は環境変数 `PYSAD_DISPLAY_DPI` か 100)。保存先がなければ、対話の backend のときだけ `plt.show` (Agg など非対話では何も出さない)。保存は `bbox_inches="tight"` で、キャンバスの余白は切る。**`ListPlot.n` 自身が `Show` を 2 回呼ぶ** (ファイルの中でブロックが重複している) ので、直前と全く同じ内容の `Show` が続いたときは出力を 1 つにまとめる (同じ式を 2 回続けて描いても 1 枚。`mpl_display.render` は常に描く)。

**対応する要素とオプション**

| 要素・オプション | 扱い |
|----|----|
| `Line` | 折れ線 (`PlotColor`・`Thickness` (1 px = 0.8 pt)・`Dashing` (`{1}` = 実線、実数 d = {d, d/3}、リストはインチ)・`Clip`・`Plot->True` の点・`PlotJoined->True/False/Step` と `StepRatio`・`Legend`)。4 列 {x,y,dx,dy} は誤差棒 (`ErrorBarColor`・`ErrorBarTickSize`・`ErrorBarThickness`) |
| `Point` | 記号つきの点 (`PointSymbol` の `CI`・`LT`・`RT`・`DT`・`UT`・`SQ`・`DI`・`PL`・`TI` と `1O` 系、`PointSize`・`PointColor`・`PointBorderColor`) |
| `Rectangle` | 箱 (`FillColor` が文字列なら塗る、`PlotColor` が枠線、`""` は枠線なし)。中に `Graphics` を持つものは入れ子の枠 |
| `Polygon`・`Circle` | 多角形 (`FillColor`・`PlotColor`・`Thickness`・`Dashing`・`PlotJoined`)、楕円 (データ座標の半径) |
| `Text` | 位置は `Scaled` (データ座標)・`WindowScaled`・`RectScaled`・`CanvasScaled`・数 (キャンバスのインチ)。`Anchor`・`TextAlign` (Tk の anchor: n・s・e・w・c・ne・nw・se・sw)・`TextRotate`・`Rotation`・`TextSize`・`PlotColor`・`Justify`・`Clip`・`Tags`。SAD の文字の書式 (`` `f ``=ギリシャ、`` `d ``=下付き、`` `u ``=上付き、`` `n ``) は mathtext に |
| 色 | Tk の名前 (`dodger blue`・`tomato3` (X11 の 1〜4 = 明るさ 100・93.2・80.4・54.5%)・`gray40`・`#rrggbb`・`SystemButtonText`)。評価されていない `PlotColorDefault[…]` (OpticsPlot) は `Plot$DefaultColor` の表から枠の中の順に選ぶ |
| `PlotRange`・`DataRange`・`Scale`・`AspectRatio`・`PlotRegion` | `Canvas.n` の `OverrideOptions`・`FrameRange`・`MakeWindow` どおり (`PlotRange` は `Automatic`・`{y1,y2}`・`{{x1,x2},{y1,y2}}`・`{x範囲, Automatic}` など。`Scale` は Linear・Log) |
| `Frame`・`FrameTicks`・`TickSize`・`FrameColor`・`TickColor`・`Background` | `Frame->False` は枠・目盛り・格子なし。`FrameTicks` は `{下,左,上,右}` で `True` = 線と文字、`Ticks` = 線だけ、`False` = なし (下の辺の線は Canvas.n どおり常に引く)。目盛りは内向き。枠をまたいで重なる目盛りの文字は後の方を消す |
| `GridLines` | `Automatic`・`Major` = 主目盛りの点線 + 原点の線、`Minor`・`All`・`Both`・`Axis`・`Axis[x…]`・`None`。**拡張**: 位置のリスト `{{1,2},{5,10}}` はその位置に線を引く (`Canvas.n` の `DrawGrid` は何も描かない) |
| `FrameLabel`・`FrameUnitLabel`・`PlotLabel`・`Legend`・`LegendBG` | `FrameLabel` は `{下,左,上,右}` (文字の書式つき)。`FrameUnitLabel` は下 (x) の単位だけ。`PlotLabel` は一番上の枠の上。凡例は枠ごと (同じ名前が続けば 1 つ)、位置は matplotlib の `best` |
| `FontScale`・`TickFontScale`・`FrameFontScale`・`LegendFontScale`・`TitleFontScale` | 文字の大きさに掛ける |
| `RemoveOverlap` | `Canvas.n` の `Remove$Overlap` と同じ手順 (その `Tags` の文字を書いた順に、残した文字の外接矩形の和に重なるものを消す)。外接矩形は matplotlib のフォントなので、Tk の窓で残る名前とは違う (OpticsPlot の名前) |
| `Axes` (拡張) | `Canvas.n` にない。`Frame->False` のときだけ、原点を通る軸の線を引く |

**未対応 (警告して飛ばす。同じ種類の警告は 1 回)** `Hook[…]` (対話の束縛。InfoLabel など。黙って飛ばす)、`Arrow`・`ArrowShape` (線だけ描く)、`UseBG->True` の文字の背景、`TextFont` (フォントの名前・書体)、`MeshStyle` の stipple の模様 (塗りつぶしだけ)、`FrameTicks` の `{Ticks, 関数}` (目盛りの文字を関数で与える形。線だけ)、`FrameUnitLabel` の左 (y) 側、`Scale->Date`・`Date` の目盛り (Linear で描く)、`Text` の位置 `LeftTags`・`RightTags`・`AboveTags`・`BelowTags` (Tk の外接矩形が要る)、`Rectangle` の `Position->` の `Right`・`Left`・`Above`・`Below` 以外、3 列の誤差棒、`PointSymbol->Bar` の基線 `Offset` (0 からの垂線で近似)、`Initialize->False`・`Frame` (前の図に重ねる・再描画。新しい Figure を作る)、`Canvas3D`・`ListPlot3D` 系、`CanvasSymbol`、未知の要素 `Foo[…]` と文字列・数が要素の位置に来たもの、解釈できない色名 (黒にする)。マウス操作 (拡大・凡例の移動・InfoLabel のクリック) は無い。

**照合** Tk の窓と画像で比べたものはない。試験 [`test_v241_mpl_display.py`](../pysad/tests/test_v241_mpl_display.py) (32 件、Agg、Tk 不要) は、`Graphics[…]` の `Line`・`Rectangle`・`Text` の個数と座標が matplotlib 側のデータ (`Line2D`・`Patch`・`Text`) と一致することを数値で確かめる: `Plot[Sin[x],{x,0,6}]` (505 点が 1e-12 以内、軸の範囲が `DataRange` の 7% の余白と `PlotRange`)、`ListPlot` (折れ線・点・誤差棒・Step・Log・`PlotRange`・`GridLines`・`Frame`・`FrameTicks`・`FrameLabel`・凡例)、`HistoPlot` (箱の数)、`OpticsPlot` (3 本の線が全点 1e-12 以内、箱 192 個の座標、名前は Graphics の `Text` の順序つきの部分列で、残った名前は互いに重ならない、枠の積み方・x 範囲)、`GeometryPlot` (線 1 本 + 多角形 73 個)、手で組んだ `Graphics` (各要素の色・太さ・破線・文字の揃え・回転・クリップ)。ほかに切替、保存 (png・pdf・svg の先頭のバイト、2 回目以降の番号)、matplotlib が無いとき、未対応の警告。例は [`example/11_graphics/fn_OpticsPlot_mpl.md`](../example/11_graphics/fn_OpticsPlot_mpl.md)。

**未確認・限界** 色・線の太さ・文字の大きさの見た目は Tk の窓と比べていない (線の太さ 1 px = 0.8 pt、破線の長さ、文字の大きさ 8.5 pt × `TextSize` は目で決めた近似。キャンバス 13 × 10 インチを 0.65 倍して図にする)。`OpticsPlot` の実ラティス (大きなリング) と `BeamPlot`・`FitPlot`・`ListDensityPlot` は描いて見ていない (要素は上の表の範囲)。対話の backend の窓 (`plt.show`) は手で確認していない (Agg で保存と図のデータだけ)。PDF・SVG は保存できることだけ確認 (見た目は PNG のみ目で見た)。目盛りの位置と文字は matplotlib の自動で、`Canvas.n` の `AutoTick` とは違う (指数の表記も matplotlib 流)。凡例の位置は `Canvas.n` の `AdjustLegend` とは違う。`PYSAD_DISPLAY=mpl` のとき `DRAW` は §7b のまま (保存先は `PYSAD_DRAW_FILE`)。gs は表示環境がないので `Graphics[…]` の式の一致 (§5) までが gs との比較で、描いた図は比べられない。

## 8. 関連文書

- [`docs/37_plot_canvas.md`](37_plot_canvas.md): 描画の gs との一致、Tk の組み込み関数の表、評価器の直した点の全部
- [`docs/22_gui_tclbridge.md`](22_gui_tclbridge.md): GUI (Tkinter.n・KBFrame.n) を tkinter で動かす仕組み (描画の土台)
- [`docs/33_specfn_fit_misc.md`](33_specfn_fit_misc.md): `Fit`・特殊関数 (`FitPlot`・`FitGaussian` が使う)
- [`docs/04_sadscript.md`](04_sadscript.md): SADScript の評価器
- [`docs/11_sad_quirks.md`](11_sad_quirks.md): SAD の癖
- [`docs/history/work_unimpl.md`](history/work_unimpl.md): 未対応・未解決の最新の記録
