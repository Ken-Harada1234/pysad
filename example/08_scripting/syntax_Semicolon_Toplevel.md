# syntax_Semicolon_Toplevel (最上位の ;)
種別: 文法 / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
FFS の最上位では `;` は 1 行の中の文の区切りで、**末尾の `;` は結果の表示 (`Out[n]:= ...`) を抑える**。FFS コマンド (`USE RING; CELL; CALCULATE;`) と SADScript の式を 1 行に混ぜられる。ラティスは `FODO_ring`。

## スクリプト
```
1 + 2          (* gs: Out[1]:= 3 *)
1 + 2;         (* 表示なし *)
b = 1; c = 2; d = b + c     (* 行の最後の式の値 Out[3]:= 3 *)
z = 1 +
  2 +
  3;                (* 行が完結していなければ次行へ続く *)
USE RING; CELL; CALCULATE; Print["NX=", ...];
```

## 結果の読み方
gs は `;` のない行の式の値を `Out[n]:= 値` で表示する (1+2→`Out[1]:= 3`、`a = 5`→`Out[2]:= 5`、`x = 1` `y = 2` も)。`;` で終わる行は表示しない。1 行に複数文を書くと、表示されるのは最後の式 (`d = b + c`→3)。式が完結しない行 (`1 +`) は次の行に続く。`USE RING; CELL; CALCULATE;` のように FFS コマンドを 1 行に並べてよい。

## gs と pysad の比較
○ 一致する (2026-11-17 に修正: `;` のない最上位の式の値を `Out[n]:= 値` で表示する)。値の末尾の桁は表示の差 (NX=6.696328940901385 gs、...378 pysad)。pysad の CALCULATE は `Matched.` の行を表示する (gs も表示する。`.gs.txt` で確認)。

## 未確認・限界
- gs の `Out[n]` の番号は FFS 全体で数えられ、`In[n]:=` の表示は最後の行だけ。
- `Print` の出力は `Out` に出ない。
