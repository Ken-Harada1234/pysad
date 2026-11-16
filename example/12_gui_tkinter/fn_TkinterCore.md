# fn_TkinterCore (TkinterCore.n の Tcl インタプリタと環境)

種別: GUI の土台 (Packages/TkinterCore.n)・フォルダ 12_gui_tkinter・ラティスなし・**窓は開かない**

## 何を見る例か
`TkinterCore.n` が作る Tcl インタプリタ (`TclEval`)、Tk の版 (`Tk$FullVersion`)、フォント系 (`$FontSystem`)、既定のフォント、`Tk$Echo` (送った Tcl の文を表示)。GUI が無くても標準出力に出る部分。

## スクリプト
```
Print["Tk version=",Tk$FullVersion];  Print["FontSystem=",$FontSystem];  Print["DefaultFont=",$DefaultFont];
Print["expr 6*7=",TclEval["expr 6*7",String]];
Tk$Echo=1;  TclEval["set x 5"];  Tk$Echo=0;       (* "Tcl> set x 5" と表示される *)
```

## 結果の読み方
`Tk version=` は Tk の版、`FontSystem=X11`、`DefaultFont={"helvetica neue",12}`、`expr 6*7=42`、`Tcl> set x 5` (`Tk$Echo=1` のときの表示)、`y=10`。

## gs と pysad の比較
- 版以外は gs と pysad で一致 (○)。`Tk version`: gs は 8.6.11、pysad (Windows の Python 3.12 の tkinter) は 8.6.13。
- `$FontSystem` は pysad も `X11` を返す (gs に合わせた値。Windows の Tk の本当のフォント系ではない)。

## 未確認・限界
Tcl のエラーの形、`Tk$Echo=-1`・`-2` は試していない。
