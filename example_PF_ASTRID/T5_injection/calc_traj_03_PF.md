# calc_traj_03_PF (PF 入射軌道 その 3: Tk の窓 (KBMainFrame・Canvas・OpticsPlot) と入力欄・ボタン)
元の話題: KEK SAD チュートリアル Tutorial5「GUI 入出力化」(`calc_traj_03`) と、「グラフ表示」の `OpticsPlot`。加速器・ラティス・線は `calc_traj_01_PF` と同じ (PF 2.5 GeV、`lattice/v9_16_inj_noSF.sad`、`INJL`)。

## 何を見る例か
軌道を `OpticsPlot` で描き (ラティスの絵つき)、キッカーの値を窓の入力欄 (`Entry`) とボタン (`Button`、Command で関数を呼ぶ) で変えて、軌道と Courant-Snyder 不変量を再計算し、グラフを描き直す。教材の書き方どおり、`KBMainFrame[...]` の枠の中に `Canvas` を置き、`$DisplayFunction = CanvasDrawer; Canvas$Widget = c1;` でそこに描く。

## スクリプト
```
w = KBMainFrame["plot", plotframe, Title->"Trajectory Plot"];
c1 = Canvas[plotframe, Width->800, Height->500];
$DisplayFunction = CanvasDrawer;  Canvas$Widget = c1;                         (* 描く先を c1 に *)
MakePlot[];      (* OpticsPlot[{{{line1},{line2},{line3}, PlotJoined->True, PlotColor->{...}, ...}}, Lattice->True, Names->"B*"] *)
k4b3 = Entry[f01, Width->9, TextVariable:>k0k4, Side->"left"];               (* 入力欄: 変数 k0k4 と結ぶ *)
k4b5 = Button[f01, Text->"Weaker", Command:>(k4deltaup[]; Recalc[])];         (* ボタンを押すと KK4 を 0.2 mrad 変えて再計算・再描画 *)
```
- 窓の入力欄 (`Entry`) の値は文字列になる。数にするには `ToExpression[ToString[x]*1.00]*1.00` を通す (入力が文字列でないと `ToExpression` がエラーになるので、いったん `ToString` する)。入力欄に使った変数は型を変えられないので、入力欄用 (`k0k4`) と計算用 (`k0k4v`) を別の名前にする。
- 描画の実際のボタン操作はできないので、Command の中身 (`k4deltaup[]; Recalc[]`) を台本から直接呼んで、押した代わりにした。
- 窓の項目数は、`TclEval` で `Canvas` の全項目の個数を読んで書いた (`canvas items = ...`)。

## 結果の読み方
`GUI-start` の行 (窓を作る前の計算: 蓄積ビームの終端の x と CS 不変量)、`canvas items` (窓の Canvas の項目数)、`after Weaker` の行 (KK4 を 0.2 mrad 弱くしたあとの値)。

| | gs | pysad |
|---|---|---|
| `GUI-start`: 終端の x [mm]、CS 不変量 [nm rad] | 0.518722、391.460 | 0.518722、391.460 |
| `canvas items` (窓を作った直後) | **114** | **114** |
| `after Weaker`: KK4 [mrad]、終端の x [mm]、CS 不変量 | -2.6、0.346732、74.531 | -2.6、0.346732、74.531 |
| `canvas items` (再描画後) | **115** | **115** |

KK4 を -2.8 → -2.6 mrad にすると、蓄積ビームの終端の x が 0.5187 → 0.3467 mm、CS 不変量が 391.5 → 74.5 nm rad に下がる (`calc_traj_02_PF` の SCAN と同じ値)。

## 物理・解析式
計算の中身は `calc_traj_02_PF` と同じ (式と数値の照合はそちら)。窓の描画項目の個数 (114) は、OpticsPlot がラティスの絵 (ベンド・四極・六極の矩形と名前、枠、目盛り、3 本の軌道) を Canvas の項目として描いた数。

## gs と pysad の比較
**○ 一致 (2026-11-17 に pysad を直した)**: 計算の値 (`GUI-start`・`after Weaker`) も、窓の Canvas の項目数 (114、再描画後 115) も gs と pysad で一致する。`.pysad.txt` にエラー表示は無い。以前は窓の部分が 2 つの理由で動かなかった:
1. `Canvas$Widget = c1` のあと Tcl が `Tkinter::tclerror: 1501 (invalid command name "plotframe")` で失敗した。原因は 2 つ: (a) `e=Entry[w]`・`c1=Canvas[plotframe,…]` の値が `Widget[…]` にならず (`Tkinter.n` の入れ子の規則で `With` が呼び出し元の引数名を見ていた)、(b) `Canvas$Widget = c1` の UpValue (Canvas.n の `(Canvas$Widget=Widget[c_,_Canvas])^:=…`) が、遅延した定義を持つ記号 `The$CanvasObject := …` を照合のために評価して、既定の Canvas の窓を作っていた。いまは gs と同じく、描く先が `c1` になる (docs/20)。
2. `OpticsPlot` に自前のデータ列を渡すと `General::equalleng: Equal lengths of elements in Thread are expected` で失敗した。原因は、OpticsPlot.n の `OpticsList1$[pr_,{f_}?(Length[Dimensions[#]]>=3&),opt___]` を、`OpticsList1$[{pr1_,pr2_},f_,opt___]` より先に試すこと (gs は定義の順。pysad は特殊さの点数で後にしていた) と、`Legend -> False` を `MapThread` に渡したとき、リストでない要素を全部に繰り返すこと (gs の `MapThread[List,{{1,2},3}]` は `{{1,3},{2,3}}`)。試験 `test_opticsplot_with_own_data_lists_like_gs` (自前データの x は要素番号で、この例は `TrackParticles` が返す位置番号)。

gs の窓の画像は撮れていない (gs の X11 の窓は撮る道具が無い)。pysad の窓の画像は `calc_traj_03b_PF` を参照。

## 未確認・限界
- 実際に窓のボタンを押す操作 (ユーザーの操作への反応) は確かめていない (Command の中身を呼んで代用)。
- `KBMainFrame` は 1 つのプロセスで 1 回しか呼べない (2 回目は `KBF::main: KBMainFrame can be called only once`)。窓のレイアウト (`f01`・`Entry` の並び) の見た目は、gs の画面で確認していない。
- チュートリアルの GUI には KK1〜KK4 それぞれの入力欄・Up/Down ボタン・入射位置の入力欄がある。ここは KK1 の入力欄と KK4 のボタンだけに絞った (残りは同じ部品の繰り返し)。
- `Canvas` のサイズ (800×500) は元 (800×600) から縮めた。
