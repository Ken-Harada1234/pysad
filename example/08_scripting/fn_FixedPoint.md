# fn_FixedPoint (FixedPoint, FixedPointList, Nest, NestList)
種別: 関数 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
関数を結果が変わらなくなるまで繰り返す。$\cos x=x$ の解と $\sqrt2$ のニュートン法。

## スクリプト
```
FixedPoint[Cos, 1.0];  FixedPoint[(# + 2/#)/2 &, 1.0];  FixedPoint[(# + 2/#)/2 &, 1.0, 3];
Nest[f, x, 3];  NestList[2*# &, 1, 5];  Length[FixedPointList[Cos, 1.0]];
```

## 結果の読み方
`FixedPoint[Cos,1.0]`=0.7390851332157367 (gs)。ニュートン法は 1.414213562373095 ($\sqrt2=1.4142135623730951$)、3 回で打ち切ると 1.4142156862745097。`FixedPointList[Cos,1.0]` は 68 点。

## 物理・解析式
$\cos x=x$ の解 $x=0.739085133215161$ (Dottie 数)。gs の値は 14 桁まで一致 (SADHelp の NearlySameQ の許容による打ち切り)。$\sqrt2$ の反復 $x_{n+1}=(x_n+2/x_n)/2$ は 2 次収束で 1→1.5→1.41667→1.414216 (3 回の値と一致)。

## gs と pysad の比較
○ 一致。`FixedPoint[Cos,1.0]` の最後の桁だけ違う (gs ...367、pysad ...368、1 ulp)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
- オプション `SameTest`、`Threshold` は使っていない。
