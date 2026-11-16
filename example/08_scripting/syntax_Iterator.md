# syntax_Iterator
種別: 文法 (反復子) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
反復子の書き方一覧: `{n}`、`{i,n}`、`{i,a,b}`、`{i,a,b,step}`、実数刻み、負の刻み、多重 (`{i,..},{j,..}`)、$i$ に依存する内側の範囲。`Table`・`Do`・`Sum`・`Product`・`Scan`・`Map` で使う。

## スクリプト
```
Table[x, {3}];  Table[i, {i, 4}];  Table[i, {i, 2, 5}];  Table[i, {i, 1, 10, 3}];
Table[i, {i, 0, 1, 0.25}];  Table[i, {i, 3, 0, -1}];  Table[{i,j}, {i,2}, {j,2}];
Table[i + j, {i, 3}, {j, 1, i}];  Sum[i*j, {i,3}, {j,3}];  Table[i, {i, 5, 1}];
```

## 結果の読み方
`{n}`→{x,x,x} (変数なし)、`{i,4}`→{1,2,3,4}、`{i,1,10,3}`→{1,4,7,10}、`{i,3,0,-1}`→{3,2,1,0}、`Sum[i*j]`=36、`Table[i,{i,5,1}]`={} (空)、`Sum` 空=0。リストを走る `{i, {a,b,c}}` は gs・pysad ともエラー (`Argument must be {symbol, n}, {symbol, start, stop}, {symbol, start, stop, step}`)。

## gs と pysad の比較
○ 一致 (リスト反復子のエラーは両方で出る。メッセージの体裁だけ違う)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
