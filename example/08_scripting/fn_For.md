# fn_For
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`For[ini, cond, incr, body]`。SAD では `Packages/functions.n` に `For[ini_,cond_,incr_,body_] := (ini; While[cond, body; incr];)` と定義されている。

## スクリプト
```
t = 0;
For[i = 0, i < 5, i++, t += i^2];   (* body は 4 番目 *)
For[i = 0, i < 3, i++, Print["body i=", i]];
```

## 結果の読み方
期待どおりなら t=30 (0+1+4+9+16)、i=5、`body i=0,1,2` を表示。**gs の実際: `t=0 i=5`、`Print` も出ない** (ループ回数は進む = `incr` は実行されるが、4 番目の引数の `body` が評価されない)。同じ gs で `For[i=0, i<3, Print["incr"]; i++]` (3 引数) も何も起きず、4 引数の 3 番目を `Print[..]; i++` にすると `incr` は表示される。

## gs と pysad の比較
△ **差あり (意図的に合わせない)**: pysad は `body` を実行して t=30、`body i=0,1,2` も表示する。gs (K64、2026-10-19 のビルド) は `body` を実行せず、3 番目の `incr` だけを繰り返す。`Packages/functions.n:34` の定義 `For[`ini_,`cond_,`incr_,`body_]:=(ini;While[cond,body;incr];)` は本体を実行する形で、`src/` に For の組み込みはなく、同じ式を別名で定義すると gs も本体を実行する (`Clear[For]` で再定義しても実行する)。したがって gs の動作は定義とも SAD の意図とも合わない K64 の癖 (読み込み時の名前解決と思われる。原因は未特定。docs/21 B71) と判断し、pysad は定義どおり本体を実行する。3 引数の `For` は gs と同じく未評価で返す。`For` は避けて `Do`・`While` を使うのが安全。

## 未確認・限界
- `For` の 3 引数形は gs で何も起きず (エラーも出ない)、pysad は実行する。
