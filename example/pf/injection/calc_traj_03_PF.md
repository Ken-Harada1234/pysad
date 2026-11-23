# calc_traj_03_PF (Tk の窓、入力欄・ボタン、OpticsPlot)
KEK の SAD チュートリアルの「GUI 入出力化」(`calc_traj_03`) と、「グラフ表示」の `OpticsPlot` をもとにした例。加速器・ラティス・線は `calc_traj_01_PF` と同じ (PF 2.5 GeV、`lattice/v9_16_inj_noSF.sad`、`INJL`)。

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
計算の中身は `calc_traj_02_PF` と同じ (式と数値の照合はそちら)。窓の描画項目の個数 (114) は、OpticsPlot がラティスの絵 (ベンド・四極・6 極の矩形と名前、枠、目盛り、3 本の軌道) を Canvas の項目として描いた数。

## gs と pysad の比較
**v2 モードと k64 モード**: `.pysad.txt` は既定の v2 モードの出力。v2 は K64 のバグを含まない計算なので (`docs/28_v2_vs_k64.md`。主に `Emittance[]`・EMIT・追跡の Jx の端の項 B75 と `MomentumCompaction` B77)、この例の出力のうち `Emittance[]` の表示 (εx が数 % 大きい、εz・減衰率・α・バケット高さなど) と、それに依る量が gs と違う。以下の「一致」の判定は k64 モードの結果。gs と同じ値が要るときは `PYSAD_COMPAT=k64` を付ける。

**○ 一致**: 計算の値 (`GUI-start`・`after Weaker`) も、窓の Canvas の項目数 (114、再描画後 115) も gs と pysad で一致する。`.pysad.txt` にエラー表示は無い。`Canvas$Widget = c1` の描く先は gs と同じく `c1` になる (docs/20)。`OpticsPlot` に自前のデータ列を渡す動作は試験 `test_opticsplot_with_own_data_lists_like_gs` が確かめる (自前データの x は要素番号で、この例は `TrackParticles` が返す位置番号)。

gs の窓の画像は撮れていない (gs の X11 の窓は撮る道具が無い)。pysad の窓の画像は `calc_traj_03b_PF` を参照。

## 未確認・限界
- 実際に窓のボタンを押す操作 (ユーザーの操作への反応) は確かめていない (Command の中身を呼んで代用)。
- `KBMainFrame` は 1 つのプロセスで 1 回しか呼べない (2 回目は `KBF::main: KBMainFrame can be called only once`)。窓のレイアウト (`f01`・`Entry` の並び) の見た目は、gs の画面で確認していない。
- チュートリアルの GUI には KK1〜KK4 それぞれの入力欄・Up/Down ボタン・入射位置の入力欄がある。ここは KK1 の入力欄と KK4 のボタンだけに絞った (残りは同じ部品の繰り返し)。
- `Canvas` のサイズ (800×500) は元 (800×600) から縮めた。
