# fn_NumericalFunctions (Abs Ceiling Floor Round Max Min MinMax Mod Sign FractionalPart ほか)

種別: 関数 (SADHelp functions の Numerical-functions)・フォルダ 10_math_data・ラティスなし

## 何を見る例か
数値を整える・比べる関数の約束 (端数処理の向き、負数の `Mod`、`FractionalPart` の符号など) を、小さな入力で確かめる。

## スクリプト
```
Abs[-3.5] Abs[3+4I]   Ceiling[-2.3]   Floor[-2.3]   Round[2.5] Round[-2.6]   Max[{1,9},{3,4}]   MinMax[{3,-1,7,2}]
Mod[-7,3]   Mod[7.5,2]   Sign[-2.]   FractionalPart[-3.75]   Positive[2] Negative[2] NonNegative[0]
```

## 結果の読み方
各行が `名前={結果のリスト}`。真偽は 1 (True) と 0 (False)。

| 式 | 結果 (gs・pysad とも同じ) | 説明 |
|---|---|---|
| `Abs[-3.5]`, `Abs[2]`, `Abs[3+4I]` | 3.5, 2, 5 | 複素数は大きさ |
| `Ceiling[2.3]`, `Ceiling[-2.3]` | 3, -2 | 切り上げ |
| `Floor[2.7]`, `Floor[-2.3]` | 2, -3 | 切り捨て (負でも小さい側へ) |
| `Round[2.4]`, `Round[2.6]`, `Round[-2.6]`, `Round[2.5]`, `Round[3.5]` | 2, 3, -3, 3, 4 | 最近接。0.5 は絶対値が大きい側へ (Mathematica は偶数側へ丸めるので 2.5 は 2。違う) |
| `Max[1,5,3]`, `Max[{2,8,4}]`, `Max[{1,9},{3,4}]` | 5, 8, 9 | リストを平らにして最大 |
| `Min[1,5,3]`, `Min[{2,8,4}]` | 1, 2 | |
| `MinMax[{3,-1,7,2}]` | {-1, 7} | 最小と最大 (SADHelp の `*` つき、SAD 独自) |
| `Mod[7,3]`, `Mod[-7,3]`, `Mod[7.5,2]` | 1, 2, 1.5 | 結果は除数と同じ符号 (`Mod[-7,3]`=2) |
| `Sign[-2.]`, `Sign[0]`, `Sign[5]` | -1, 0, 1 | |
| `FractionalPart[3.75]`, `FractionalPart[-3.75]` | .75, .25 | $x-\lfloor x\rfloor$ (Mathematica の -.75 とは違い、常に 0 以上) |
| `Positive[2]`, `Negative[2]`, `NonNegative[0]` | 1, 0, 1 | |

## 物理・解析式
式で値が決まる性質: $\mathrm{Mod}(a,b)=a-b\lfloor a/b\rfloor$、$\mathrm{FractionalPart}(x)=x-\lfloor x\rfloor$、$|3+4i|=5$。表の結果はこれで確かめられる ($-7-3\lfloor-7/3\rfloor=-7+9=2$、$-3.75-(-4)=0.25$)。

## gs と pysad の比較
全ての行が一致 (○)。出力は全く同じ文字列。

## 未確認・限界
`Round` の「0.5 は絶対値が大きい側」は 2.5→3、3.5→4、-2.5 は未確認 (-2.6 だけ試した)。`Mod` の負の除数、`Max`・`Min` の引数なし、複素数の `Round` は試していない。
