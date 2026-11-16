# 39. Tcl/Tk の GUI の解説 (SAD の Tkinter.n / KBFrame.n を pysad で動かす)

最終更新: 2026-11-16

SAD の GUI (窓、ボタン、入力欄、キャンバス、メニュー) を、SAD 自身の台本 (`Packages/Tkinter.n`、`KBFrame.n` など) を書き換えずに
pysad のインタプリタで解釈し、Python の tkinter (Tk 8.6) に表示する仕組みの解説。初めて読む人向け。
要素別の詳しい実装の記録は `docs/20_gui_tclbridge.md` (接点と言語機能)、`docs/35_plot_canvas.md` (描画)。経緯は `docs/history/work_gui_tclbridge.md`。

## 1. 何ができるか

- SAD の GUI 台本 (`w = Window[]`、`b = Button[w, Text->"Go", Command:>...]`、`KBMainFrame[...]`、`ListPlot[...]` など) が、そのまま動く。
- 方針は「GUI のコードを Python へ訳さない」こと。SAD が C / Fortran で持っているのは Tcl との接点 (`src/tfTkInter_.c`、`src/tfTclArg.f`、`src/tfXlib_.c`) だけで、
  `.n` のパッケージは SADScript で書かれている。pysad はその接点だけを Python で作り、`.n` は SAD のまま解釈する。
  (リポジトリ直下 `CLAUDE.md` の「例外: GUI は書き直さず .n を解釈する」に当たる。)
- 層の構成:

```
SAD 台本 (Window[]、KBMainFrame[...]、ListPlot[...])
   │   Packages/*.n を pysad の SADScript インタプリタが解釈 (Class は sadclass.py)
   ▼
Tkinter.n / TkinterCore.n / KBFrame.n / Canvas.n   →  TclEval["..."] などの組み込み関数を呼ぶ
   ▼
pysad/tclbridge.py (Tcl/Tk の接点)、pysad/tkcanvas.py (キャンバス項目の作成)
   ▼
Python の tkinter (Tk 8.6) が窓を表示
```

## 2. 使い方

SAD の `Packages` の場所を環境変数で渡す。**`PYSAD_SAD_PACKAGES` が無いと、GUI は何も有効にならず、GUI の自動試験は skip になる。**

```
set PYSAD_SAD_PACKAGES=D:/work/SAD-master/SAD-master/Packages
py -3.12 -c "from pysad.sadscript import Interp; I=Interp(); I.run(open('demo.sad').read(),'script')"
```

| 環境変数 | 意味 |
| --- | --- |
| `PYSAD_SAD_PACKAGES` | SAD の `Packages`。これがあるときだけ `tclbridge.install` が働く。`Tkinter.n` などの名前は最初に使われたときに読み込まれる (`init.n` の `AutoLoad` と同じ) |
| `PYSAD_TK_NOLOOP=1` | `TkEventLoop[]` が 1 回 `update` して戻る (自動試験・ウィジェット木の比較用)。無ければ本物の mainloop (窓を閉じるまで戻らない) |
| `PYSAD_NO_DISPLAY=1` | `$DisplayFunction` を `Identity` にする。`ListPlot` などは窓を開かず `Graphics[...]` の式を返すだけ (`pkgload.py`)。Linux で `DISPLAY` が無いときも同じ |
| `PYSAD_TRACE=1` | パッケージの読み込み中に失敗した文を標準エラーに出す |
| `PYSAD_TRACE_TCL=1` | Tcl に送った命令を標準エラーに出す |

FFS から使うときは `FFS;` の後に同じ台本を書く (`run_sad_string` 経由。`ListPlot` の例は docs/35)。

## 3. 使える部品

`Tkinter.n` の `Widgets` 一覧 (Tk の部品の対応) と、`KBFrame.n` / `TabFrame.n` の部品。

| SAD の部品 | Tk のクラス | 備考 |
| --- | --- | --- |
| `Frame`、`Window`、`LabelFrame` | Frame、Toplevel(`.`)、Labelframe | `Window[]` は最上位の窓 |
| `PanedWindow` | Panedwindow | |
| `TextLabel`、`Button`、`CheckButton`、`RadioButton` | Label、Button、Checkbutton、Radiobutton | `Command:>(...)` が押されたとき評価される |
| `MenuButton`、`Menu`、`OptionMenu` | Menubutton、Menu、optMenu.tcl | |
| `Entry`、`SpinBox` | Entry、Spinbox | `TextVariable:>変数` で変数と結ぶ |
| `ListBox` | Listbox | `ListVariable` |
| `Scale`、`ScrollBar` | Scale、Scrollbar | |
| `Canvas` | Canvas | 描画は `tkcanvas.py` (docs/35) |
| `TextEditor`、`TextMessage` | Text、Message | |
| `Pack`、`Grid`、`Place` | 配置 | オプションは `Side`、`Fill`、`Expand`、`Row`、`Column` など |
| `KBMainFrame`、`KBFrame`、`KBSubFrame`、`KBDialogFrame`、`KBFComponentFrame`、`KBFGroup`、`KBFCheckButton`、`KBFHook`、`KBFFrameToggled`、`KBFFrameTabed`、`KBFTabFrame` ほか | 上の部品の組み合わせ | KEK の標準の枠 (メニューバー、状態行、時計) |

### 3.1 動く SAD スクリプト例

下は試験 (`pysad/tests/test_tclbridge.py`) にある台本。`Interp()` に `run(src, 'script')` で渡す。

簡単な部品 (4 個。Window、Button、CheckButton、TextLabel):

```
w = Window[];
b = Button[w, Text->"Hi"];
c = CheckButton[w, Text->"Chk"];
l = TextLabel[w, Text->"Lab"];
TclEval["update"];
info = TclEval["set r {}; foreach c [winfo children .] {lappend r [winfo class $c] [$c cget -text]}; return $r", String];
```
`info` は `Button Hi Checkbutton Chk Label Lab` になる (SAD の gs と一致)。

メインフレーム (メニューバー付き):

```
mf = KBMainFrame["TkTest", f, Title->"Test"];
TclEval["update"];
```
`.` の下に Menubutton 6 個 (KEK、File、Edit、Window、Task、Help)、仕切り、時計の Label ができる。

LabelFrame のラベルにチェックボタンを埋め込む入れ子の部品:

```
w = KBMainFrame["LabelFrameDemo", f, Title->"LabelFrame Demo"];
cb[0] = CheckButton[Text->"Enable", Variable:>cbv[0],
  Command:>(Scan[(cb[#][State] = If[cbv[0], "normal", "disabled"])&, Range[3]])];
cf = KBFComponentFrame[f, Add->{
  KBFHook[Command:>(cbv[0] = True)],
  KBFGroup[WidgetOption->{LabelWidget->cb[0]}],
  KBFCheckButton[Text->"A", Variable:>cbv[1], WidgetVariable:>cb[1]],
  KBFCheckButton[Text->"B", Variable:>cbv[2], WidgetVariable:>cb[2]],
  KBFCheckButton[Text->"C", Variable:>cbv[3], WidgetVariable:>cb[3]]}];
TclEval["update"];
```
Checkbutton 4 個、Labelframe 2 個、Menubutton 6 個ができる (SAD と一致)。

描画 (`test_plot.py`。キャンバスに 43 項目、gs と全部一致):

```
FFS;
g=ListPlot[{1,4,9,16},PlotJoined->True];     (* 窓が開く。戻り値は Graphics[...] *)
h=Plot[Sin[x],{x,0,6}];
```

2026-11-16 に追加で動かして確かめた最小の例 (試験には入れていない。ウィジェットが作られ、`Command` が評価されるまで):

```
w = Window[];
ev = "abc";
e = Entry[w, TextVariable:>ev];
b = Button[w, Text->"Go", Command:>(n = 5)];
Scan[Pack, {e, b}];
TclEval["update"];
(* Tcl 側から押す: b を invoke すると n = 5 になった *)
```
同様に `Entry`、`ListBox`、`Scale`、`RadioButton`、`Canvas`、`TextEditor` を `Window` の下に作り `Pack` すると、`winfo children .` が
`Entry Listbox Scale Radiobutton Canvas Text Button` の順に返った (エラーなし)。
ウィジェットへの操作は `TclEval["... "]` で Tk の命令を直接送る形でも、`e[Insert]={"end","XY"}` (Entry・ListBox) の形でも書く (後者は `Tkinter.n` の定義)。Insert の位置は Tcl の添字なので **小文字の `"end"`** を使う。大文字の `"End"` は gs で Tcl エラー 1501 になる。
注意 (2026-11-16 の gs との比較): `e=Entry[w]` のあと、**gs では `e` が `Widget[e,Entry[w]]` を値に持ち**、`e[Insert]={0,"AB"}`・`e[Get]`・`e[Delete]`・`l[Insert]` (Listbox) が効く (Entry の内容 AB→CD、Listbox 2 項目、`v=e[Get]` が AB)。**pysad では `e` が裸のシンボルのまま** (`Head[w]=Symbol`、gs は `Widget`) で、`e[Insert]=…`・`e[Get]` は何も呼ばない (Tcl 呼び出し 0 回、エラーもなし)。ウィジェット自体は作られる (`winfo children` に Entry・Listbox が出る)。`w=Window[]; b=Button[w,...]` でも同じで、pysad は親が `w` でなく Root になる (`winfo` で `.` の直下に並ぶ。gs は Toplevel の下)。これは**既知の差 (未修正)**: `Tkinter.n` の `DefineWidget` 内の `a$:=Widget[a$,w[master$]]` が pysad で `w` の値にならない。root の Window は toplevel を作らず、`TkPathName[w]` が "" のまま。`f[x_]:=(x=5); f[aa]` のようなパターン変数への代入は gs と同じに働くので、原因は未特定 (玩具の再現は gs も pysad も再帰になり、実物が gs で再帰しない条件は未特定)。pysad はこの再帰で Python の `RecursionError` で落ちることがある (gs は `General::deep`)。

### 3.2 試験の台本の探し方

- `pysad/tests/test_tclbridge.py` (4 本): Tk の版と `$FontSystem`、簡単な部品、`KBMainFrame` のメニューバー、`KBFComponentFrame`。
- `pysad/tests/test_gui_lang.py` (10 本): GUI を使わない言語機能 (`Sym::tag`、`f@@[x,{1}]`、`Null[]`、`0x33`、UpValue、SubValue、Hold 属性、`Difference`、`SwitchCases`、`StringFill`、`BaseForm` など)。
  `Tkinter.n` を動かすために足したもので、`PYSAD_SAD_PACKAGES` が無くても走る。
- `pysad/tests/test_sadclass.py` (16 本): `Class` (docs/20 §5)。
- `pysad/tests/test_plot.py` (8 本): `ListPlot`、`OpticsPlot`、`HistoPlot`、`ListContourPlot`、`FitPlot`、`Plot`、`GeometryPlot`、`FitGaussian` のキャンバス。

## 4. 実装の構成

### 4.1 `pysad/tclbridge.py` (接点)

- `TclCreateInterp`、`Tcl` / `TclEval`、`TclSetResult`、`TclArg` / `TclArg1`、`TkEventLoop` / `TkQuitEventLoop`、`TkDefineBitmap`。
  Tcl に `sadeval` / `sadcommand` を `createcommand` で登録し、Tk のイベント (ボタン、メニュー) から SAD の式を呼ぶ。
- `TclArg` は `tfTclArg.f` の `tftclarggen` と同じ規則 (規則 `a->b` はラベルと値、リストは `{}`、`Widget[...]` はパス名、`String[...]` はそのまま)。
- 実数は整数なら整数で Tcl に渡す (`-width`・`after` は整数を要る)。`$FORM = "6.2"` の台本でも `1001.` にならないよう、`StandardForm` の中で評価する。
- Xlib の関数 (`XQueryTree`、`XInternAtom` など) は「他のアプリは見つからない」と答える空の実装。
- EPICS (`CaGet`、`CaPut`、`CaMonitor` …) と `System[]` は実行しない。
- 窓の更新は `TclEval["update"]`。台本の最後の `TkEventLoop[]` は、`PYSAD_TK_NOLOOP=1` なら 1 回の `update` で戻る。

### 4.2 `pysad/tkcanvas.py` (キャンバス)

`Canvas.n` が使う `TkCanvas*` (`TkCanvasGetTypes`・`TkCanvasPointer`・`TkCanvasCreateItem(Direct)`・`TkCanvasLastItem`・`TkCanvasFindEnclosed`)、
`CanvasClipLine`・`CanvasSymbolDirect` (`tfcanvasclip.f`) を実装。K64 は Tk の内部構造を C で直接たたくが、pysad は tkinter の Tcl に `create` を送る。
`TkDefineBitmap` は Windows の Tcl が名前つきビットマップを登録できないので XBM ファイルに書き、`-stipple` の名前を `@ファイル` に替える。詳しくは docs/35 §3。

### 4.3 言語機能 (`sadscript.py`、`sadlist.py`、`sadclass.py`)

`Tkinter.n`・`KBFrame.n` を動かすため、UpValue・SubValue・保持属性・パターン照合・`Class` などを K64 の実測に合わせて足した (docs/20 §4、§5)。
`Class` は「メンバーをインスタンスごとの記号 (`c$i12`) に改名する」方式。

## 5. ヘッドレスでの試験

- 通常: `set PYSAD_SAD_PACKAGES=...` のうえで
  `py -3.12 -m pytest pysad/tests/test_tclbridge.py pysad/tests/test_sadclass.py pysad/tests/test_gui_lang.py pysad/tests/test_plot.py -q`。
  全体は `py -3.12 -m pytest pysad/tests -q --ignore=pysad/tests/test_bend_feeddown.py --ignore=pysad/tests/test_coord.py --ignore=pysad/tests/test_sext_radiation.py` で 520 通過・3 skip (523 件収集、2026-11-16)。
- 試験の fixture は `PYSAD_TK_NOLOOP=1` を設定し、`tkinter.Tk()` が作れなければ skip する (画面が無い環境)。
- Windows の Tcl はまれに起動直後に "Can't find a usable tk.tcl" で失敗するので、`tkinter.Tk()` は最大 5 回やり直す (ブリッジ側は 4 回)。
- **窓を出さずに描画だけ確かめる**: `PYSAD_NO_DISPLAY=1` で `ListPlot` などを呼ぶと `Graphics[...]` の式 (`Head` = Graphics、`Length` = 2) が返る (`test_listplot_returns_graphics_and_plot_matches_gs`)。
  この方式では Tk を使わない。
- **ウィジェット木で比べる**: 台本を走らせたあと `winfo children` を再帰して、クラス・`-text`・`-state`・`-title` を書き出し、SAD (gs) の同じ台本の木と比べる。
  `local/gui_baseline/run_py_tree.py` がその道具 (`PYSAD_SAD_PACKAGES` と `PYTHONPATH` を設定して `py -3.12 run_py_tree.py demo.sad`)。
- **キャンバスの項目で比べる**: `$c find all` と `$c type`・`$c coords`・`$c gettags` を書き出す (`test_plot.py` の台本)。
- **サンプル一括**: `local/gui_samples/run_all.sh` が、サニタイズした台本の写し (EPICS と `System[]` を `$off` の名前に替えたもの) を 1 本ずつ 120 秒の制限で走らせる。
- 全体試験では GUI 系 (`test_plot`、`test_tclbridge`) が数回に 1 回、単独なら通るのに落ちることがある (Tk の窓の更新の待ち)。落ちたら単独で回し直す。

## 6. K64 の癖

pysad が合わせた K64 の挙動。

- SAD の Tk は改造版。`-scriptfont`・`-altfont`・`-altscriptfont` (ウィジェット全般) と、text 項目の `-usebg` は標準の Tk に無いので、pysad は送らない/落とす
  (文字の背景は塗らない)。
- 既定のフォント名は X11 のもの (`$DefaultSansSerifFont` = "helvetica neue" など、`$FontSystem` = "X11")。Windows の Tk と見た目が違う。
- `Tkinter.n` の `FontSubstitute$X11` に `Returun[...]` の綴りの誤りがあり、「フォント表が見つからない」ときの早期復帰が効かない (結果は変わらない)。
- SAD の `Sequence` は `Null`。`Null@@list` を引数列として展開しないと `-bg` などのオプションが消える。
- `KBMainFrame` の Menubutton は 6 個 (7 ではない)。`KBF$ConfDir` はソース配布に `share/KBFrame` が無いので `<SAD>/KBFrame` に落とす。
- ビットマップは XBM 形式の改行位置を守らないと `Image` が置換されず `DeleteImage` が再帰する。
- 値が `Null` の規則 (`Fill->Null`) は Tk の項目に渡さない。gs は既定の色 (X11 で `#000000`) になり、Windows の Tk では `SystemButtonText` になる (色の違い)。
- `DateString[]` は `mm/dd/yyyy hh:mm:ss`。`ToLowerCase` / `ToUpperCase` は文字列以外をそのまま返す。
- Tk の `create` が失敗しても SAD は `???Tkinter::tclerror: ...` のメッセージだけで `Null` を返して続ける (pysad も同じ)。

## 7. 検証の度合い (どの部品を実際に確かめたか)

基準は WSL の SAD (gs、Tk 8.6.11 X11、WSLg の `DISPLAY=:0`) で同じ台本を走らせて取ったウィジェット木・キャンバス項目。

| 部品・台本 | 確かめた内容 | 結果 |
| --- | --- | --- |
| `Window`・`Button`・`CheckButton`・`TextLabel` | 木 (クラス、`-text`) | gs と一致 (`test_tclbridge.py`) |
| `KBMainFrame` | メニューバー (Menubutton 6 個と Menu、Frame、時計の Label)、状態行、進捗表示 | 一致 |
| demo-LabelFrame (`KBFComponentFrame`、`KBFGroup`、`KBFCheckButton`) | 木 36 行 | 一致。差は時刻の文字列と、状態行の進捗表示の Label が 1 個多いこと (原因未確認) |
| demo-KBFrame1 (`KBFFrameToggled`・`KBFFrameTabed`) | 木 | **途中まで** (SAD 57 行に対し pysad 43 行)。入れ子の部品の配置 `KBF$ArrangeComponents` で失敗 |
| `Canvas` (ListPlot・Plot・ParametricPlot・ListDensityPlot・HistoPlot・GeometryPlot) | キャンバスの全項目 (種類・座標・タグ) | gs と一致 (文字の外接矩形だけ Tk のフォントで違う)。詳細は docs/35 §2、§5 |
| `OpticsPlot` | 枠・矩形・線 | 一致。text は gs 58 に対し pysad 72 (フォントで `RemoveOverlap` が変わる) |
| `ListContourPlot` | 種類と数 | 一致。等高線の始点・向きは違う |
| サンプル 161 本 (`D:\work\SAD\GUI`) | EPICS と `System[]` を無効にした複製で実行 | Tk を開く約 50 本は最後まで走りエラー無し。残り約 81 本は Tk を使わない計算の台本 |
| `Entry`・`ListBox`・`Scale`・`RadioButton`・`Canvas`・`TextEditor` (単体) | 作成し `Pack` して `update` (2026-11-16) | エラー無しで作成された (木のクラス名のみ確認) |
| `Button` の `Command` | Tcl から `invoke` | SAD の式が評価された (`n = 5`) |

試験: `test_tclbridge.py` 4、`test_sadclass.py` 16、`test_gui_lang.py` 10、`test_plot.py` 8。

## 8. 限界・未確認

- **見た目 (色、大きさ、フォント) と操作 (ボタンを押す、メニューを開く、拡大) は SAD と比べていない。** 比べたのはウィジェット木とキャンバスの項目の座標。
  Windows の Tk は X11 の Tk と違う (フォント名、`-scriptfont`、`wm` の挙動、X のウィンドウ木)。表示の違いはここから出る。
- 木を gs と比べたのは、上の表の台本 (Window 系、KBMainFrame、demo-LabelFrame、キャンバス系)。残りのサンプルは「エラー無しで木ができた」までで、gs とは未比較。
- `Entry`・`ListBox`・`Scale`・`RadioButton`・`TextEditor`・`SpinBox`・`PanedWindow`・`OptionMenu`・`TextMessage`・`ScrollBar` は、
  単体の台本で gs と木を比べていない。**ウィジェットのメソッド (`e[Insert]`・`e[Get]`・`e[Delete]`・`l[Insert]`) は gs では効くが pysad では効かない** (ウィジェットが `Widget` 値にならない既知の差、§3.1)。`TextVariable` で結んだ変数への値の反映は確認していない。
- `TabFrame.n` の描画 (`KBFTabFrame`、`KBFToggleFrame`) は未検証。demo-KBFrame1 は入れ子の部品の配置で失敗する (`cf[13]` が未定義。原因は未解決)。
- `Class` の未対応: `親` メンバーの書き方、上位値 (`h_[f_[arg___], b___] ^:= g`) の形のメソッド、`Members[]` / `Parents[]`。
- キャンバス: `CanvasSymbol` (Tcl 文字列を返す版)、`Canvas3DClipTriangle`・`Canvas3DLightTriangle`・`Canvas3DProjection` (Canvas3D.n) は未対応。`BeamPlot` は比べていない。`UseBG->True` の文字の背景は描かない。
- EPICS 関係 (`CaGet` など) と `System[]` は動かさない。
- 使った環境の制約: Windows、Python 3.12、Tk 8.6 (tkinter)。画面が無いと Tk の試験は skip になる。GUI の台本は、本物の mainloop で長く動かし続ける使い方 (ユーザーの操作に反応し続ける) は確認していない。

## 9. 関連文書

- `docs/20_gui_tclbridge.md` (接点・言語機能・Class・検証の詳細)
- `docs/35_plot_canvas.md` (描画と `tkcanvas.py`、gs との一致)
- `docs/history/work_gui_tclbridge.md` (実装の経緯、落とし穴、未解決)
- `docs/04_sadscript.md`、`docs/11_sad_quirks.md` (SADScript と K64 の癖)
- `docs/32_sad_packages.md` (SAD パッケージ)、`docs/34_spin_csr.md` (`CSR.n`)
- `docs/history/00_RESUME.md` (引き継ぎ、環境、試験の手順)
