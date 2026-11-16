# fn_If (If, 論理演算, FORTRAN 風 IF-ELSEIF-ELSE-ENDIF)
種別: 関数 (制御構文) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`If` の 2〜3 引数、条件は実数 (0 が偽)、`&&` `||` `~` の短絡評価、FFS の FORTRAN 風 `IF ... ELSEIF ... ELSE ... ENDIF`。

## スクリプト
```
f[x_] := If[x > 0, "pos", If[x < 0, "neg", "zero"]];
If[1 > 2, "yes"];       (* 偽で 2 引数: Null *)
IF n < 0 ... ELSEIF n == 0 ... ELSEIF n < 10 ... ELSE ... ENDIF;
```

## 結果の読み方
`If` は選んだ枝だけ評価する (使わない枝の `Print` は出ない)。`If[2,...]` は真、`If[0,...]` は偽。`(0 && Print[..])` は Print を評価しない。IF ブロックは `n=5` で "small positive" の枝だけ実行。

## gs と pysad の比較
○ 一致 (出力の全行が同じ)。

## 未確認・限界
- 出力は `.gs.txt` / `.pysad.txt` の通り。SADHelp の記述 (Expression の項) と gs の実動作を照合した。
- 乱数・時刻・ファイルの順序に依存する出力は入れていない。
