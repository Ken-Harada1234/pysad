# fn_Cases_DeleteCases_Position
種別: 関数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
パターンに合う要素を取り出す (`Cases`)、除く (`DeleteCases`)、位置を返す (`Position`)、数える (`Count`)、含むか (`MemberQ`)。

## スクリプト
```
l = {1, 2, a, 3, b, 2, {4,5}};
Cases[l, _Real];  Cases[l, _Symbol];  Cases[{f[1],g[2],f[3]}, f[x_] -> x];
DeleteCases[l, _Real];  Position[l, 2];  Position[{{1,2},{3,2}}, 2];
```

## 結果の読み方
`Cases[l,_Real]`={1,2,3,2}、`_Symbol`={a,b}。`Cases[..., f[x_]->x]` は規則の右辺を返して {1,3}。`Position[l,2]`={{2},{6}} (添字のリスト)、入れ子は {{1,2},{2,2}}。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
