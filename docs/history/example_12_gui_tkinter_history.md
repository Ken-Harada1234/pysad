# example/12_gui_tkinter の履歴 (取り除いた記述)

md から、版・日付つきの修正・経過の記述を取り除いた。取り除いた元の行 (整理前) を、ファイルごとにそのまま残す。現在の説明は各 md を見る。

## example/12_gui_tkinter/fn_Entry_ListBox_Scale.md

整理前の 3 行目:

> 種別: GUI (Packages/Tkinter.n)・フォルダ 12_gui_tkinter・ラティスなし。pysad でも動く (2026-11-17 に修正。親の違いが残る)

整理前の 28 行目:

> - **○ メソッドも一致 (2026-11-17 に修正)**: `e[Insert]={"end","XY"}` で Entry が `abcXY`、`lb[Insert]` 2 回で ListBox が 2 項目、`Head[e]` は `Widget` で gs と同じ。以前は `Head[e]` が `Symbol` で `e[...]` が何も呼ばなかった (`Tkinter.n` の入れ子の規則で `With` が呼び出し元の引数名 (ユーザーの `w`) を見て、`a$:=Widget[a$,…]` が別名に付いていた。docs/22、試験 `test_widget_values_and_methods_like_gs`)。

## example/12_gui_tkinter/fn_Menu_MenuButton.md

整理前の 26 行目:

> - **窓の親が違う (pysad の既知の差、docs/41 §3.1・§8)**: gs では `Window[]` が Toplevel を作り、部品はその下に並ぶ (ルート `.` の下には「Return to SAD」ボタンと Toplevel)。2026-11-17 の修正で Menu が Menubutton の下に入る段は gs と同じになった (`Menubutton File;   Menu`)。残る差は、pysad では `Window[]` が別の窓を作らず、部品がルート窓 `.` の直下に並ぶ。`tree:` の行のインデントの差はこのため。


## 追加 (整理の 2 回目)

### example/12_gui_tkinter/README.md

整理前の 10 行目:

> | [fn_Entry_ListBox_Scale](fn_Entry_ListBox_Scale.md) | Entry・ListBox・Scale・RadioButton | △ (メソッドは 2026-11-17 に一致。窓の親 (Toplevel) と変数名の番号が違う) |

### example/12_gui_tkinter/fn_Window_Button_TextLabel.md

整理前の 32 行目:

> - **`Head[w]`・`Head[b]`**: gs・pysad とも `Widget` (○。2026-11-17 に修正: `Tkinter.n` の入れ子の規則で `With` が呼び出し元の引数名を見ていた。docs/22)。


## 台本の冒頭コメントから取り除いた履歴 (v2_next・v2_B・v2_C から移したとき)

### Tk90_shot.py (旧 v2_next/P5_tk90_shot.py)
> """v2next P5 項目 7: Tk 9.0 (Python 3.14) の GUI の見た目を撮る。…
使い方の行の旧パス example/v2_next/P5_tk90_shot.py を現在のパスに直した。
