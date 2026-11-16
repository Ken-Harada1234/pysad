# fn_Function_Pattern (defining-functions, pattern)
種別: 関数 (定義とパターン) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
関数の定義 `f[x_] := ...`、定数引数の定義 (`fib[0] = 0` が優先)、型つきパターン `x_Real`、`x_?test`、列パターン `x__` `x___`、既定値 `n_:2`、名前付き `x:pat`、`/.`・`//.`・`MatchQ`・`Alternatives`・`Repeated`、純関数、`TagSet`、メモ化、`Clear`。

## スクリプト
```
fib[0] = 0; fib[1] = 1; fib[n_] := fib[n-1] + fib[n-2];
kind[x_Real] := "real"; kind[x_String] := "string"; kind[x_List] := "list"; kind[x_] := "other";
pos[x_?Positive] := "positive";   cnt[x___] := Length[{x}];   pw[x_, n_:2] := x^n;
MatchQ[f[1,2,3], f[_Real..]];   mem[n_] := mem[n] = n^2;   tg /: tg[a_] := "tagged";
```

## 結果の読み方
`fib` は {0,1,1,2,3,5,8,13,21,34,55}、`kind` は {"real","string","list","other"}、`pw[3]`=9、`pw[3,3]`=27、`cnt[]`=0、`Replace[f[2], f[a:_] -> g[a,a]]`=`g[2,2]`、`MatchQ[f[], f[_Real...]]`=1。`Clear[sq]` の後の `sq[3]` は未定義でそのまま。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。
- **条件パターン `x_ /; test` (Condition) は SAD にない**: 試した結果 gs は黙って常に偽の枝に入り (`pos[2]` が "non-positive")、pysad は `unexpected token ';'` で失敗する。検査は `x_?test` (PatternTest) を使う。
- **`Function[x, body]` 形 (Mathematica の名前付き Function) も SAD にない** (gs は壊れた値、pysad は `5+n` のように未評価)。純関数は `# ... &`。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
- 左辺の `f[x_ : 既定値]`、`Optional`、`Pattern` の細部、`Orderless` は見ていない。
