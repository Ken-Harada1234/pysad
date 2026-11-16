# fn_Switch
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Switch[x, pat1, v1, pat2, v2, ...]` はパターンに合う最初の値を返す。型 (`_String`, `_Real`, `_List`)、選択 (`2 | 3`)、既定 (`_`)。

## スクリプト
```
h[x_] := Switch[x, 1, "one", 2 | 3, "two-three", _String, "str", _, "other"];
Switch[9, 1, "one"];   (* どれにも合わない *)
Switch[f[3], f[x_], x*10, _, 0];
```

## 結果の読み方
{"one","two-three","str","other","other"}。`Switch[f[3], f[x_], x*10, _, 0]` は `(10 x)` (gs・pysad とも): パターンの名前 `x` は右辺に束縛されない (この SAD の Switch の癖)。型で分岐する書き方 `sw[x_] := Switch[x, _String, ..., _Real, ...]` は期待どおり動く。

## gs と pysad の比較
○ 一致する (2026-11-17 に修正: 合うものがないとき `Switch[9]` を返す。引数がすべて定数なら gs が読み込み時に畳んで `Switch[9]` になり、変数を含むと `Switch[5,1,2]` のように第 1 引数だけ評価した式を返す)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
