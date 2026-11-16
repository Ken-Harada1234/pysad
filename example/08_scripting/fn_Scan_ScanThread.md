# fn_Scan_ScanThread
種別: 関数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Scan` は結果を返さずに各要素へ関数を適用 (副作用用)、`ScanThread` は `MapThread` の結果を返さない版。

## スクリプト
```
Scan[Print, {1,2,3}];  s = 0; Scan[(s += #) &, Range[10]];
ScanThread[Print["pair ", #1, " ", #2] &, {{1,2},{a,b}}];
ScanThread[Print, {{1,2,3},{10,20,30}}];
```

## 結果の読み方
合計 55。`ScanThread[Print, {{1,2,3},{10,20,30}}]` は `Print[1,10]` などを呼ぶので "110"、"220"、"330" (区切りなしで連結される)。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
