# Tcl/Tk の GUI の検証の詳細 (gs との比較)

本文 [../20_tcltk_gui.md](../20_tcltk_gui.md) から移した、K64 (gs) との細かい比較。元の文書: 41_tcltk_guide (§7)、22_gui_tclbridge (§6)。見出しの節番号は元の文書のもの。

---
元: 22_gui_tclbridge.md

## 6. 検証と現状

基準は、WSL の SAD (gs、Tk 8.6.11 X11、WSLg の DISPLAY=:0) で同じ台本を走らせて取ったウィジェット木 (`winfo children` を辿った、クラス・`-text`・`-state`・`-title`)。

| 台本 | 結果 |
| --- | --- |
| Tkinter の `Window` / `Button` / `CheckButton` / `TextLabel` | オプションまで一致 ([`test_tclbridge.py`](../../pysad/tests/test_tclbridge.py)) |
| `KBMainFrame` | メニューバー (Menubutton 6 個: KEK、File、Edit、Window、Task、Help と、その Menu、仕切りの Frame、時計の Label)、状態行、進捗表示まで一致 |
| KBFrame の demo-LabelFrame (`KBFComponentFrame`、`KBFGroup[WidgetOption->{LabelWidget->cb[0]}]`、`KBFCheckButton`) | 木の形が一致 (36 行)。違いは時刻の文字列と、状態行の進捗表示 2 番目の Label が 1 個多いこと (SAD では作られていない。原因は未確認) |
| KBFrame の demo-KBFrame1 (`KBFFrameToggled` / `KBFFrameTabed` を含む入れ子の部品) | 途中まで動く (SAD の 57 行に対して 43 行)。入れ子の部品の配置 (`KBF$ArrangeComponents` の `w[iwsym]` を使う親の背景の参照) で、まだ失敗する |

サンプル一括実行 (2026-11-13): GUI のサンプル台本 161 本を、EPICS と `System[]` を無効にした複製で実行した。Tk を開く約 50 本は全て最後まで走り、Tcl・pysad ともエラー無し。
残り 81 本は Tk を使わない計算の台本。SAD との木の比較は上の表の 3 本だけで、他は「エラー無しで木ができた」までの確認。

試験: [`test_tclbridge.py`](../../pysad/tests/test_tclbridge.py) 7、[`test_sadclass.py`](../../pysad/tests/test_sadclass.py) 16、[`test_gui_lang.py`](../../pysad/tests/test_gui_lang.py) 10 の計 33 本 (2026-11-17 に収集数を確認)。GUI 系は `PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>` を付けたときだけ実行する。全試験は 575 通過・3 skip (578 件収集、2026-11-17)。全体を回すと `test_tclbridge`・`test_plot` が数回に 1 回落ちることがある (単独では通る。Tk の窓の更新の待ち)。

**Widget 値**: `w=Window[]`・`e=Entry[w]` の値は gs と同じ `Widget[…]` (`Head[w]=Widget`)。`e[Insert]={"end","XY"}`・`lb[Insert]`・`e[Get]` も効く。
- 仕組み: `DefineWidget` の `(a$_=w[opt$___])^:=(a$=w[Widget[Root$Window,Null],opt$])` が次の規則 `(a$$_=w[Widget[master$_,_],opt$___])^:=With[{a$=a$$,…},…]` を呼ぶ入れ子である。`With` が呼び出し元の規則の引数 `a$` (= ユーザーの `w`) を見つけ、`a$` でなく `w` を束縛しがちである (`sadscript._with` はこの探索を別の規則の呼び出しの境界 (`boundary` の frame) で止める)。
  探索を止めないと本体の `a$:=Widget[a$,…]` が別の名前に付いて、ユーザーの `w` に遅延した値が付かない。
- `Canvas$Widget = c1` (Canvas.n の UpValue) は、左辺の記号が最初の使用で読み込みの引き金になり、UpValue の照合で遅延した定義 (`The$CanvasObject := …`) を持つ記号を評価しない。教材の台本の `Canvas$Widget=c1` で tclerror が出ず、窓の Canvas に 114 項目 (再描画後 115) で gs と同じ。
- **残る差**: gs は最初の `Window[]` が Toplevel (`.5`、K64 起動時の `.` には「Return to SAD」ボタンがある) だが、pysad は `Window[]` が `.` 自身。木の段と `TextVariable` の名前の番号が違う。
- 補足: `Insert` の位置に "End" (大文字) を渡すと gs は Tcl エラー 1501 (Tcl の添字は "end")。

未確認・限界:
- ウィジェット木 (クラス、文字、状態) は比べたが、見た目 (色、大きさ、フォント) と、操作したときの動作 (ボタンを押す、メニューを開く) は SAD と比べていない。
- `Canvas.n` の描画は [`pysad/tkcanvas.py`](../../pysad/tkcanvas.py) で実装し、`ListPlot` の 43 項目を gs と比べて一致 (docs/37)。`TabFrame.n` の描画は未検証。
- EPICS 関係 (`CaGet` など) は動かさない。`System[]` も動かさない。
- Windows の Tk は X11 の Tk と違う (フォント名、`-scriptfont`、`wm` の挙動、X のウィンドウ木)。表示の違いはここから出る。

---
元: 41_tcltk_guide.md

## 7. 検証の度合い (どの部品を実際に確かめたか)

基準は WSL の SAD (gs、Tk 8.6.11 X11、WSLg の `DISPLAY=:0`) で同じ台本を走らせて取ったウィジェット木・キャンバス項目。

| 部品・台本 | 確かめた内容 | 結果 |
| --- | --- | --- |
| `Window`・`Button`・`CheckButton`・`TextLabel` | 木 (クラス、`-text`) | gs と一致 ([`test_tclbridge.py`](../../pysad/tests/test_tclbridge.py)) |
| `KBMainFrame` | メニューバー (Menubutton 6 個と Menu、Frame、時計の Label)、状態行、進捗表示 | 一致 |
| demo-LabelFrame (`KBFComponentFrame`、`KBFGroup`、`KBFCheckButton`) | 木 36 行 | 一致。差は時刻の文字列と、状態行の進捗表示の Label が 1 個多いこと (原因未確認) |
| demo-KBFrame1 (`KBFFrameToggled`・`KBFFrameTabed`) | 木 | **途中まで** (SAD 57 行に対し pysad 43 行)。入れ子の部品の配置 `KBF$ArrangeComponents` で失敗 |
| `Canvas` (ListPlot・Plot・ParametricPlot・ListDensityPlot・HistoPlot・GeometryPlot) | キャンバスの全項目 (種類・座標・タグ) | gs と一致 (文字の外接矩形だけ Tk のフォントで違う)。詳細は docs/37 §2、§5 |
| `OpticsPlot` | 枠・矩形・線 | 一致。text は gs 58 に対し pysad 72 (12 セルの例。FODO 24 セルでは 69 対 74。フォントで `RemoveOverlap` が変わる) |
| `ListContourPlot` | 種類と数 | 一致。等高線の始点・向きは違う |
| GUI のサンプル台本 161 本 | EPICS と `System[]` を無効にした複製で実行 | Tk を開く約 50 本は最後まで走りエラー無し。残り約 81 本は Tk を使わない計算の台本 |
| `Entry`・`ListBox`・`Scale`・`RadioButton`・`Canvas`・`TextEditor` (単体) | 作成し `Pack` して `update` (2026-11-16) | エラー無しで作成された (木のクラス名のみ確認) |
| `Button` の `Command` | Tcl から `invoke` | SAD の式が評価された (`n = 5`) |

試験: [`test_tclbridge.py`](../../pysad/tests/test_tclbridge.py) 7、[`test_sadclass.py`](../../pysad/tests/test_sadclass.py) 16、[`test_gui_lang.py`](../../pysad/tests/test_gui_lang.py) 10、[`test_plot.py`](../../pysad/tests/test_plot.py) 8。
