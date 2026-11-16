# fn_Class (Class[], メンバー・メソッド・継承・コンストラクタ)
種別: 関数 (クラス化) / フォルダ: 08_scripting / FFS 内の SADScript (ラティスなし)。gs は SAD 1.4.4.2k64、pysad は `run_sad_string`。

## 何を見る例か
`Class[{スーパークラス}, {クラス変数}, {インスタンス変数}, メソッド...]` でクラスを定義する (SADHelp の Object-oriented-programing、`Packages/Class.n`)。コンストラクタ、`@` でメンバー・メソッドの呼び出し、クラス変数 (全インスタンスで共有)、継承とオーバーライド、`Parents[]`、破棄 `x=.`。

## スクリプト
```
Acc = Class[{}, {count = 0}, {v = 0, name = "none"},
  Constructor[a___] := (count++; v = Plus[a]);
  add[x_] := (v += x; v);  show[] := Print[name, ": v=", v, " count=", count];  sum2[] := v^2];
a1 = Acc[1, 2];  a1@add[10];  a1@show[];  a2 = Acc[100];  a2@name = "second";
Sub = Class[{Acc}, {}, {w = 7}, add[x_] := (v += 2*x; v); both[] := v + w];
s1 = Sub[1];   s1@add[3];   a1=.;
```

## 結果の読み方
`a1@v`=3 (コンストラクタで `Plus[1,2]`)、`a1@add[10]` 後の `show[]` は `none: v=13 count=1`。2 つ目のインスタンス `a2` を作ると `count`=2 でクラス変数 `Acc@count` も `a1@count` も 2。`Sub` は `Acc` を継承 (`s1@sum2[]`=49 は `Acc` のメソッド)、`add` を上書き (`s1@add[3]` は `v += 2*3` で 7)、`both[]`=`v + w`=14。`a1=.` で破棄しても `a2` は残る。

## gs と pysad の比較
○ 一致 (値)。違うのは `Sub@Parents[]` の表示だけ (gs `{Acc`}`、pysad `SadClass#2@4`)。

## 未確認・限界
- 多重継承、`Destructor` の再定義、`Short`、`This` を明示する書き方、`Members[]`・`AllMembers[]` は見ていない。
- pysad は SAD 自身の `Class.n` を解釈するのではなく、Python 側の実装 (`SadClass`) でクラスを扱う。
