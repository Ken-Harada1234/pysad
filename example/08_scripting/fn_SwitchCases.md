# fn_SwitchCases (SwitchCases, SelectCases)
種別: 関数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`SwitchCases[list, {case1,..}]` は 1 回の走査で、case ごとに合う要素のリストを返す (Cases と DeleteCases の同時実行)。`SelectCases` は検査関数版。最後に `_` や `True &` を置くと「どれにも合わない要素」が入る。

## スクリプト
```
l = {1, 2.5, "a", x, 3, "b", {1}, y};
SwitchCases[l, {_String, _Real}];  SwitchCases[l, {_String, _Real, _}];
SelectCases[Range[10], {EvenQ, # > 7 &}];  SelectCases[Range[10], {EvenQ, # > 7 &, True &}];
```

## 結果の読み方
`SwitchCases[l,{_String,_Real,_}]` = {{"a","b"}, {1,2.5,3}, {x,{1},y}}。`SelectCases[Range[10], {EvenQ, #>7&}]` = {{2,4,6,8,10},{9}}: 先に合った 8,10 は 2 つ目に入らない (最初に合った case に入る)。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
