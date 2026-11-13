# GUI (Tkinter.n / KBFrame.n) の実装記録

作業日: 2026-11-12 〜 2026-11-14。現状の説明は `docs/20_gui_tclbridge.md`。ここは経緯と落とし穴。

## 依頼
「リスト処理、KBFrame 関数、Tcl/Tk の GUI 関数の実装へ進んで」。SAD の GUI 台本をそのまま動かすため、コードの翻訳ではなく、
SAD 自身の `Tkinter.n` / `KBFrame.n` / `KBMainFrame.n` / `TabFrame.n` / `Canvas.n` を pysad のインタプリタで解釈し、Python の tkinter で表示する方針にした。
サンプルは `D:\work\SAD\GUI` (EPICS 関連は実行しない指示)。マニュアルは `SAD-master\doc` と `SADTkinter.pdf` (後から提供。実装は .n のソースと gs の実測を基準にしたので PDF は未参照)。

## 進め方
1. `tclbridge.py` に C 側の接点 (`TclCreateInterp`、`Tcl`、`TclArg`、`TkEventLoop` ほか) を実装。
2. `sadclass.py` に `Class` を実装 (メンバーをインスタンスごとの記号へ改名する方式)。
3. KBMainFrame を読み込み、WSL の gs で同じ台本のウィジェット木を取って比べ、差が出た言語機能を一つずつ足した。
4. サンプル 161 本を、EPICS と `System[]` を `$off` の名前に替えた複製で一括実行 (`local/gui_samples/`)。

## 一括実行の結果 (2026-11-13 夜)
- 161 本中、Tk を開くのは約 50 本。すべて最後まで走り、Tcl のエラーも pysad のエラーも無し。
- 81 本は計算だけで Tk を使わない (木は空、実行確認にならない)。
- SAD の gs と木を比べたのは demo-LabelFrame と KBMainFrame だけ。残りは未比較。見た目と操作は未確認。

## 落とし穴 (再発防止)
- Windows の Tk は起動直後にまれに "Can't find a usable tk.tcl" で失敗する。`tkinter.Tk()` を最大 4 回やり直す (ブリッジと試験の fixture)。
- `$FORM = "6.2"` の台本では Tcl に送る数が `1001.` になる。`StandardForm` 内で評価して既定の書式にし、実数は整数なら整数で送る。
- SAD の `Sequence` は `Null`。`Null@@list` を引数列として展開しないと `-bg` などのオプションが消える。
- ビットマップは XBM 形式 (`static char x_bits[] = {...};`) の改行位置を守らないと `Image` が置換されず、`DeleteImage` が再帰する。
- `KBF$ConfDir` はソース配布に `share/KBFrame` が無いので `<SAD>/KBFrame` に落とす。
- 子クラスの定義が親の定義に消される不具合は、親を先に登録し、同じ具体性なら新しい定義を先頭に置いて解決。
- bash の heredoc で python を書くと文字が壊れる。スクリプトを `cat > file` で書いて `py -3.12 file` で実行する。Windows の python は `/tmp` を見えない。
- Menubutton の数は 7 ではなく 6 (SAD の実測)。

## 未解決
- demo-KBFrame1: `frame … -bg cf[13][Background]` で `cf[13]` が未定義。入れ子の部品 (KBFFrameToggled / KBFFrameTabed) の配置で、`w[iwsym]=wsym` (`wsym=Unevaluated$[tf1]`) の扱いか Class の相互作用を疑っている。
- demo-LabelFrame の状態行に、SAD に無い Label (ProgNum2) が 1 個多い。原因未確認。
- Canvas / TabFrame の描画 (`TkCanvas*`) は最小実装で未検証。

## 次の作業
利用者の指示: GUI が終わったら wake と beam-beam の実装に入る。実例の台本は無いので、SADHelp (CERN) と K64 のソースを根拠にする。
