# fn_Do (Do, Break, Continue, Catch, Throw)
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Do` の反復子 `{n}`、`{i,n}`、`{i,a,b,step}`、多重反復、`Continue`・`Break`、`Catch`・`Throw` による大域脱出。

## スクリプト
```
Do[Print["Do i=", i], {i, 3}];  Do[..., {i, 2, 10, 4}];  Do[Print[i,",",j], {i, 2}, {j, i}];
s = 0; Do[s += i, {i, 100}];   (* 5050 *)
Do[If[i == 3, Break[]]; ..., {i, 10}];   Catch[Do[If[i == 3, Throw[i*100]], {i, 10}]];
```

## 結果の読み方
`{i,2,10,4}` は 2, 6, 10。`{j,i}` は $i$ に依存する内側の反復 (1,1 / 2,1 / 2,2)。`Catch[Do[... Throw[300]]]`=300。`Do` 自体の値は `Null`。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。
- **`{i, {5,7,11}}` (リストを走る反復子) は gs も pysad もエラー** (SAD の反復子は `{n}`、`{i,n}`、`{i,a,b}`、`{i,a,b,step}` のみ。リストを走るには `Scan`・`Map`)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
