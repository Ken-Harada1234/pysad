# fn_TableForm (tableform.n の TableForm)

種別: パッケージ (Packages/tableform.n)・フォルダ 13_packages・ラティスなし

## 何を見る例か
`TableForm[list]` はリスト (行列・ベクトル・スカラー) を、桁をそろえた表の文字列にする。1 次元リストは 1 行、スカラーは 1 セル。幅は `PageWidth` と列数で決まる (1 セルは 7〜17 文字)。

## スクリプト
```
m={{1.5,2.25,3},{4,5.125,6},{7,8,9.0625}};
Print[TableForm[m]];   Print[TableForm[{1,2,3,4}]];   Print[TableForm[3.5]];
```

## 結果の読み方
行列は 3 行 3 列で、各セルが 17 文字幅 (左詰め)。ベクトルは 4 列 1 行、スカラーは 1 セル。

## gs と pysad の比較
- 表の中身 (数字・桁・並び) は gs と pysad で一致 (○)。
- pysad は `Print[TableForm[...]]` のあとに空行が 1 行余分に出る (表示の差)。

## 未確認・限界
文字列の要素 (`{"x","y"}`) は gs・pysad とも表が出なかったので例から外した。`PageWidth` を変えたときの列幅は試していない。
