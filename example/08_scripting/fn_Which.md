# fn_Which
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Which[c1, v1, c2, v2, ...]` は最初に真になった条件の値を返す。

## スクリプト
```
g[x_] := Which[x < 0, "neg", x == 0, "zero", x < 10, "small", True, "large"];
Which[1 > 2, "a", 3 > 4, "b"];   (* どれも真でない *)
```

## 結果の読み方
{"neg","zero","small","large"}。最後の `True` が既定の枝。どれも真でないと、gs は評価されないまま `Which[1>2,"a",3>4,"b"]` を返す。

## gs と pysad の比較
○ 一致する (2026-11-17 確認: どれも真でないときは式をそのまま返す)。引数の個数が合わなければ `General::narg`。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
