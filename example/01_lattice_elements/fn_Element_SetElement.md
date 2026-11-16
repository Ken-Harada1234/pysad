# fn_Element_SetElement (Element と SetElement)
種別: SADScript 関数。フォルダ: 01_lattice_elements。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
- `Element[キー, 要素名]`: **要素** (名前の単位。同名の成分をまとめた電源のようなもの) の値を返す。`"VALUE"` は既定キーワードの値、`"TYPENAME"`、`"POSITION"`、`"KEYWORDS"` (使えるキーワードの一覧) などの特別なキーがある。代入 `Element["K1","QF"] = 0.8;` は全成分に配られる。
- `SetElement[名前, 型, キー->値]`: MAIN 部の要素を作る・値を変える (型は変えられない)。戻り値は同じ要素を再作成できる形のリスト。
- 同名の成分 1 台だけを見るのは `LINE[...]` (fn_LINE)。

## スクリプト
```
Print["Element K1 QF = ", Element["K1","QF"], "  L = ", Element["L","QF"], "  VALUE = ", Element["VALUE","QF"]];
Print["Element TYPENAME QF B RF = ", Element["TYPENAME",{"QF","B","RF"}]];
Print["Element POSITION QF = ", Element["POSITION","QF"], "  NAME = ", Element["NAME", Element["POSITION","QF"]]];
Element["K1","QF"] = 0.8;   CALC;
r = SetElement["QF", "QUAD", "K1" -> 0.7];   CALC;
```

## 結果の読み方
- K1=0.7、L=0.3、TYPENAME は {"QUAD","BEND","CAVI"}、QF の位置は 5 (要素リストの中の番号)。`KEYWORDS B` は BEND のキーワード一覧。
- `Element["K1","QF"]=0.8` のあと NX = 8.046572 (全 24 台の QF が 0.8 になった)、`LINE["K1","QF.5"]` も 0.8。`SetElement` で 0.7 に戻すと NX = 6.696329。

## 物理・解析式
QF の K1 を 0.7→0.8 にしたときのチューンシフトの線形見積もり:
$$\Delta\nu_x\approx\frac{1}{4\pi}\sum_{QF}\beta_x\,\Delta K_1=\frac{24\times\bar\beta_{x}\times0.1}{4\pi}\ (\bar\beta_x\approx6.7\ \text{m})=1.28$$
実測は 1.35 で、見積もり (1.28) に 5 % の差で近い (K1 を 14 % 変えるので 2 次以上の項が残る)。

## gs と pysad の比較
全桁一致。

## 未確認・限界
- `Element[...] = v` の代入が EXPAND や `Saved->True` とどう関係するか (SADHelp の記述) は未確認。
- `SetElement` で新しい要素を作る使い方 (ラインに入れるには MAIN 部の定義が必要) は未確認。
