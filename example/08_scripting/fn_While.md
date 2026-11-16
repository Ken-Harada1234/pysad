# fn_While
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`While[cond, body]`、`Break`、`Continue`、複合代入 (`+=`、`*=`、`/=`、`++`)。

## スクリプト
```
n = 0; While[n < 5, n++];
k = 1; While[True, k *= 2; If[k > 100, Break[]]];   (* k=128 *)
a = 100; b = 0; While[a > 1, a /= 2; b++];
m = 0; While[m < 6, m++; If[OddQ[m], Continue[]]; Print["even ", m]];
```

## 結果の読み方
n=5、k=128、100 を 1 以下になるまで半分にする回数 7 (a=0.78125)、偶数 2,4,6。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
