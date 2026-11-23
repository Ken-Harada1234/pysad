# fn_Apply_Thread (Apply, @@, Thread)
種別: 関数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Apply` (`@@`) は先頭 (ヘッド) を差し替え、`Thread` は引数のリストを要素ごとに展開する。

## スクリプト
```
Apply[Plus, {1,2,3}];  Times @@ {1,2,3,4};  f @@ {a,b};  Apply[f, {{1,2},{3,4}}, {1}];
Thread[f[{1,2,3}, x]];  Thread[f[{1,2},{a,b}]];  Thread[{x,y} -> {1,2}];  Thread[{1,2} == {1,3}];
```

## 結果の読み方
`Apply[Plus,{1,2,3}]`=6、`Times@@{1,2,3,4}`=24、`f@@{a,b}`=`f[a,b]`。`Thread[f[{1,2,3},x]]` は {f[1,x], f[2,x], f[3,x]}、規則は {x->1, y->2}。`Thread[{1,2}=={1,3}]` は {1,0}。

## gs と pysad の比較
○ 一致する (2026-11-17 確認: リストどうしの `==` は要素ごとに比べるので、`Thread[{1,2}=={1,3}]` は gs と同じ {1,0})。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
