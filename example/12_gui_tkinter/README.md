# 12_gui_tkinter: Tkinter・KBFrame の GUI

SAD の `Tkinter.n`・`KBFrame.n`・`TkinterCore.n` の最小の例。画面表示が主役。窓の画面 (`NAME_window.png`) は pysad の窓を `example/tools/gui_shot.py` で切り取ったもの。gs の窓は撮れていない。標準出力には部品の木 (`tree:`) と、ボタンを押した結果だけが出る。

一致の記号: ○ 一致 / △ 差あり (md に理由) / ×gs 失敗 / ×pysad 未対応・失敗

| 例 (NAME.sad・.gs.txt・.pysad.txt・.md) | 内容 | gs と pysad |
|---|---|---|
| [fn_Window_Button_TextLabel](fn_Window_Button_TextLabel.md) | Window・TextLabel・Button・CheckButton | ○ (窓の親と Head[w] の差は既知) |
| [fn_Entry_ListBox_Scale](fn_Entry_ListBox_Scale.md) | Entry・ListBox・Scale・RadioButton | △ (メソッドは 2026-11-17 に一致。窓の親 (Toplevel) と変数名の番号が違う) |
| [fn_Canvas](fn_Canvas.md) | Canvas と項目の描画 | ○ |
| [fn_Menu_MenuButton](fn_Menu_MenuButton.md) | MenuButton・Menu | ○ |
| [fn_KBMainFrame_KBFComponentFrame](fn_KBMainFrame_KBFComponentFrame.md) | KBMainFrame・KBFComponentFrame・KBFCheckButton | ○ |
| [fn_TkinterCore](fn_TkinterCore.md) | TkinterCore.n (Tcl インタプリタ、Tk の版)、窓なし | ○ (Tk の版は違う) |

pysad では `Window[]` が別の窓を作らず、ウィジェットが `Widget[...]` の値にならない (docs/39 の未解決点)。動かないものは md に書いた。
