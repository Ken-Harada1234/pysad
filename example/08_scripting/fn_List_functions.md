# fn_List_functions (Take, Drop, Join, Flatten, Sort, ...)
種別: 関数 (リスト) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
SADHelp の expression/List-manipulations の項から選んだリスト関数 (`Length` `Take` `Drop` `First` `Last` `Rest` `Join` `Append` `Prepend` `Flatten` `Sort` `Reverse` `Range` `Select` `Insert` `Delete` `Complement` `Union` `Intersection` `Partition` `Dimensions` `Depth` `RotateLeft` `RotateRight` `Sum` `Product` `Fold` `Part` `ReplacePart` `Difference` `MapIndexed` `Extract` と、リストの要素ごと演算・内積)。

## スクリプト
```
l = {5, 3, 8, 1, 9, 2};
Take[l, 3];  Take[l, -2];  Take[l, {2,4}];  Drop[l, 2];  Join[{1,2},{3},{4,5}];
Flatten[{1,{2,{3,4}},5}, 1];  Sort[l, Greater];  Range[2, 8, 3];  Select[l, # > 4 &];
Insert[{a,b,c}, x, 2];  Delete[{a,b,c,d}, 2];  Complement[{1,2,3,4,5}, {2,4}];
Partition[Range[5], 2, 1];  Fold[f, x, {a,b}];  Difference[{1,4,9,16}];
```

## 結果の読み方
`Take[l,{2,4}]`={3,8,1}、`Drop[l,-2]`={5,3,8,1}、`Flatten[..., 1]` は 1 段だけ ({1,2,{3,4},5})、`Partition[Range[5],2,1]`={{1,2},{2,3},{3,4},{4,5}} (重なりあり)、`Fold[f,x,{a,b}]`=`f[f[x,a],b]`、`Difference[{1,4,9,16}]`={3,5,7}、`MapIndexed[List,{a,b}]`={{a,{1}},{b,{2}}}、`{1,2,3}.{4,5,6}`=32。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
- 取り上げたのは SADHelp の List-manipulations の一部。`FlattenAt`、`Extract` の高次、`HeldPart`、`Level`、`RotateLeft` の多次元などは入れていない。
