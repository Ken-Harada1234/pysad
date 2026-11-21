# 42. グラフィック描画の解説 (Plot・ListPlot・OpticsPlot ほかと Tk のキャンバス)

最終更新: 2026-11-17

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

## 8. 関連文書

- [`docs/37_plot_canvas.md`](37_plot_canvas.md): 描画の gs との一致、Tk の組み込み関数の表、評価器の直した点の全部
- [`docs/22_gui_tclbridge.md`](22_gui_tclbridge.md): GUI (Tkinter.n・KBFrame.n) を tkinter で動かす仕組み (描画の土台)
- [`docs/33_specfn_fit_misc.md`](33_specfn_fit_misc.md): `Fit`・特殊関数 (`FitPlot`・`FitGaussian` が使う)
- [`docs/04_sadscript.md`](04_sadscript.md): SADScript の評価器
- [`docs/11_sad_quirks.md`](11_sad_quirks.md): SAD の癖
- [`docs/history/work_unimpl.md`](history/work_unimpl.md): 未対応・未解決の最新の記録
