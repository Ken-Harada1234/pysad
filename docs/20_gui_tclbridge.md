# 20. Tcl/Tk の GUI (Tkinter.n / KBFrame.n) を pysad で動かす仕組み

最終更新: 2026-11-13 (現状のスナップショット。経緯は `history/` にある)

対象: `pysad/tclbridge.py` (Tcl/Tk との接点)、`pysad/sadclass.py` (SAD の Class)、`pysad/sadscript.py` / `sadlist.py` の言語機能の追加分、
試験 `pysad/tests/test_tclbridge.py` / `test_sadclass.py` / `test_gui_lang.py`。

## 1. 方針

GUI のコードを Python へ訳すのではなく、SAD が自分で持っている SADScript のパッケージ
(`Packages/Tkinter.n`、`TkinterCore.n`、`KBFrame.n`、`KBMainFrame.n`、`TabFrame.n`、`Canvas.n` ほか) を pysad のインタプリタで解釈して動かす。
C / Fortran で書かれているのは、Tcl との接点 (`src/tfTkInter_.c`、`src/tfTclArg.f`、`src/tfXlib_.c`) だけなので、その分を Python の tkinter で補う。
表示は Python の tkinter (Tk 8.6) がする。SAD 側の台本 (`w = KBMainFrame[...]`、`cf = KBFComponentFrame[...]` など) は書き換えずに動くことを目標にしている。

## 2. 使い方

SAD のパッケージの場所を環境変数で渡す。渡さなければ従来どおり何も有効にならない (自動試験もこの場合は GUI の試験を飛ばす)。

| 環境変数 | 意味 |
| --- | --- |
| `PYSAD_SAD_PACKAGES` | SAD の `Packages` ディレクトリ (例: `<SAD のソース>/Packages`)。これがあるときだけ `tclbridge.install` が働く |
| `PYSAD_TK_NOLOOP=1` | `TkEventLoop[]` が 1 回 `update` して戻る (自動試験・ウィジェット木の比較用。無ければ本物の mainloop) |
| `PYSAD_TRACE=1` | パッケージの読み込み中に失敗した文を標準エラーに出す |
| `PYSAD_TRACE_TCL=1` | Tcl に送った命令を標準エラーに出す |

`Interp()` を作ると、`Tkinter.n` などの名前 (`Window`、`Button`、`KBMainFrame`、…) は `init.n` の `AutoLoad` と同じく、最初に使われたときに読み込まれる。

## 3. 実装した接点 (tclbridge.py)

- `TclCreateInterp`、`Tcl` / `TclEval`、`TclSetResult`、`TclArg` / `TclArg1`、`TkEventLoop` / `TkQuitEventLoop`、`TkDefineBitmap`。
  Tcl のコマンド `sadeval` / `sadcommand` を tkinter の `createcommand` で登録し、Tk のイベントから SAD の式を呼ぶ。
- `TclArg` は `tfTclArg.f` の `tftclarggen` と同じ規則 (規則 `a->b` はラベルと値、リストは `{}`、`Widget[...]` はパス名、`String[...]` はそのまま)。
  SAD が手を入れた Tk にしか無い `-scriptfont`、`-altfont`、`-altscriptfont` は標準の Tk に無いので送らない。
- Xlib の関数 (`XQueryTree`、`XInternAtom`、…) は「他のアプリは見つからない」と答える空の実装 (Windows の tkinter には X のウィンドウ木が無い)。
- SAD の既定値: X11 のフォント名 (`$DefaultSansSerifFont` = "helvetica neue"、serif = "utopia"、fixed = "courier"、…、WSL の gs で実測)、`$FontSystem` = "X11"、
  `ReturnToSAD` = 0、`KBF$ConfDir` = `<SAD>/KBFrame` (`share/KBFrame` が無いソース配布のとき)。
- 実数を Tcl に渡すときは `$FORM` に左右されず、整数なら整数 (Tk の `-width` や `after` は整数を要求する)、それ以外は最短の表記。
- Windows の Tcl はまれに起動直後の初期化 (tk.tcl の探索) に失敗するので、`tkinter.Tk()` は数回やり直す。
- EPICS (`CaGet`、`CaPut`、`CaMonitor` …) と `System[]` は実行しない。サンプルを動かすときは、複製の中でこれらの呼び出しを `$off` の名前に替えて無効にする。原本は書き換えない。

## 4. Tkinter.n / KBFrame.n を動かすために足した SADScript の言語機能

SAD (K64) の実測と SADHelp に合わせた。`test_gui_lang.py` と `test_sadclass.py` に試験がある。

- 式: `Sym::tag`、`f@@[x, {1}]`、`Null[]` と `Null[a, b]` と `Null@@list` (引数列として展開。SAD では Sequence = Null)、`0x33`、n 項の `<=>` と `===`、
  空の右辺の規則 `x->`、`a@b` (クラスのメンバー参照。右は記号だけ、`x@m[args]` は `(x@m)[args]`。インスタンスでなければ `f@g[x]` = `f[g[x]]` のまま)。
- 代入: UpValue (`(w_ = f[...]) ^:= …`、Set / SetDelayed の UpValue)、SubValue (`f[a][b]`、3 段の `f[a][b][c]` まで)、`a/:f[a, x] := …` (TagSet)、
  `f[x] = v` の値は v、`Unevaluated$[x]` を引数にすると変数 x そのもの (名前) を渡して関数の中の `i++` や `w[i] = …` が呼び出し側の x に届く、
  局所変数への代入で同名の大域変数の定義を消さない。
- 評価の規則: 保持属性 (HoldAll / HoldFirst / HoldRest) の引数に対する置き換え (関数の引数・With の変数を、その時点の値で定義へ埋め込む。内側の With / Module / Function /
  パターン変数が同じ名前を束縛し直したときは置き換えない)、`With[{a=1; b=2, c}, …]` (宣言の中の `;`、`{c}` は現在値)、`WithPriority` / `Priority` (先頭に定義を置く)、
  定義の並び (リテラル > 構造パターン > 型つき Blank・Alternatives > 素の Blank の順に試す)、`And` / `Or` の短絡評価、`Hold[x]` が x を評価しない、
  遅延した値を持つ変数が `f[x]` の `f` に来たとき (`fr := Widget[fr, …]` の `fr[Background]`) の適用。
- パターン: 構造パターン (`wid_[master_, opt___]`、`(Rule|RuleDelayed)[Options[w], _]`)、Rule オブジェクトとの照合 (`Rule[_String, _]`)、
  変数が持つパターン (`Switch[val, types, …]` の types)、型つき Blank、評価されていない `{x, y}` との照合、パターンを含む規則の `ReplaceAll` (`{x_->_, y_} :> …`)。
  `a:>b` の b は `ReplaceAll` でも評価しない。
- 関数: `Difference`、`SwitchCases` (リスト以外の式の引数も)、`StringFill`、`Symbol`、`Second` / `Third`、`BaseForm`、`Names`、`RGBColor` (`#RRGGBB`、tfloor.f と同じ丸め)、
  `StandardForm[expr]` (評価の間だけ `$FORM` を既定の書式にする。台本が `$FORM = "6.2"` としていても Tcl への命令の数が `1001.` にならない)、`Reverse[a->b]`、`Union` (同じ頭の式どうし)、`Extract` (規則の中)、`l[[, 2]]` (空の添字は全要素)、`DateString[]` は SAD と同じ `mm/dd/yyyy hh:mm:ss`、
  `ToLowerCase` / `ToUpperCase` は文字列以外をそのまま返す。
- Python の `TypeError` / `AttributeError` が組み込み関数から漏れたときは、SAD の式のエラー (`General::wrongtype`) にして止める。

## 5. Class (sadclass.py)

SADHelp の「Class」の節に沿った実装。`Class[親の列, クラス変数, インスタンス変数, メソッド…]`、`x = c[引数]` (規則 `a->1` はメンバーの初期値、他の規則は `Options`)、
`Constructor` / `Destructor` (`x=.` で呼ぶ)、`x@メンバー`、`x@メソッド[引数]`、`x@メンバー = 値`、`This`、親クラスの継承 (子の同じ形の定義が親のものを置き換える)、`x = c1` の複写 (同じインスタンス)。
「メンバーが出てきたら `This@メンバー` と読む」という規則は、メンバーの記号をインスタンスごとに一意な記号 (`c$i12`) に置き換えて実現している。メソッドの本体は普通の関数として評価され、
ウィジェットの生成 (`c = Canvas[...]`) や `coords[0] = …` のように記号に定義を付ける処理も、そのままインスタンスごとの記号に対して働く。
未対応: `親`メンバー の書き方、上位値 (`h_[f_[arg___], b___] ^:= g`) の形のメソッド、`Members[]` / `Parents[]` などの問い合わせ。

## 6. 検証と現状

基準は、WSL の SAD (gs、Tk 8.6.11 X11、WSLg の DISPLAY=:0) で同じ台本を走らせて取ったウィジェット木 (`winfo children` を辿った、クラス・`-text`・`-state`・`-title`)。

| 台本 | 結果 |
| --- | --- |
| Tkinter の `Window` / `Button` / `CheckButton` / `TextLabel` | オプションまで一致 (`test_tclbridge.py`) |
| `KBMainFrame` | メニューバー (Menubutton 6 個: KEK、File、Edit、Window、Task、Help と、その Menu、仕切りの Frame、時計の Label)、状態行、進捗表示まで一致 |
| KBFrame の demo-LabelFrame (`KBFComponentFrame`、`KBFGroup[WidgetOption->{LabelWidget->cb[0]}]`、`KBFCheckButton`) | 木の形が一致 (36 行)。違いは時刻の文字列と、状態行の進捗表示 2 番目の Label が 1 個多いこと (SAD では作られていない。原因は未確認) |
| KBFrame の demo-KBFrame1 (`KBFFrameToggled` / `KBFFrameTabed` を含む入れ子の部品) | 途中まで動く (SAD の 57 行に対して 43 行)。入れ子の部品の配置 (`KBF$ArrangeComponents` の `w[iwsym]` を使う親の背景の参照) で、まだ失敗する |

サンプル一括実行 (2026-11-13): `D:\work\SAD\GUI` の台本 161 本を、EPICS と `System[]` を無効にした複製で実行した。Tk を開く約 50 本は全て最後まで走り、Tcl・pysad ともエラー無し。
残り 81 本は Tk を使わない計算の台本。SAD との木の比較は上の表の 3 本だけで、他は「エラー無しで木ができた」までの確認。

試験: `test_tclbridge.py` 4、`test_sadclass.py` 16、`test_gui_lang.py` 10 の計 30 本が通る (2026-11-14 確認)。

未確認・限界:
- ウィジェット木 (クラス、文字、状態) は比べたが、見た目 (色、大きさ、フォント) と、操作したときの動作 (ボタンを押す、メニューを開く) は SAD と比べていない。
- `Canvas.n` と `TabFrame.n` の描画 (`c[Create$Polygon] = …` など) は、C の `TkCanvas*` が最小限の実装のため未検証。
- EPICS 関係 (`CaGet` など) は動かさない。`System[]` も動かさない。
- Windows の Tk は X11 の Tk と違う (フォント名、`-scriptfont`、`wm` の挙動、X のウィンドウ木)。表示の違いはここから出る。

## 7. 試験の実行

```
set PYSAD_SAD_PACKAGES=<SAD のソース>/Packages
py -3.12 -m pytest pysad/tests/test_tclbridge.py pysad/tests/test_sadclass.py pysad/tests/test_gui_lang.py -q
```

`PYSAD_SAD_PACKAGES` が無いか、tkinter が使えないときは GUI の試験は skip になる (GUI を使わない言語機能の試験は常に走る)。
