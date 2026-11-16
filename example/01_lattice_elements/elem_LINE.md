# elem_LINE (LINE の定義)
種別: elem (MAIN 部の定義)。フォルダ: 01_lattice_elements。ラティスは短い輸送路 (QUAD 2 種 + DRIFT。ファイル内)。エネルギー 3 GeV。

## 何を見る例か
`LINE 名前 = (要素 要素 …)` の書き方。**ラインの先頭は MARK** (FFS で使うため)。
- `n*要素` で n 回の繰り返し。`n*名前` はラインの繰り返し。
- `-名前` は逆向きの配置 (要素の並びが逆順になり、向き `DIR` が -1)。
- ラインの中にラインを入れ子にできる。`USE` したときは最下層まで展開され、同じ名前の要素は `Q.1`, `Q.2`, … と番号が付く。

## スクリプト
```
LINE HALF = (Q D QQ D);
LINE TWICE = (M 2*HALF);          ! M Q D QQ D Q D QQ D
LINE REVL  = (M -HALF);           ! 逆向き
LINE NEST  = (M HALF -HALF 2*Q);  ! 入れ子
FFS;
USE TWICE; INS; CALC; DISP ALL;   Print["nlat = ", LINE["LENGTH"]];
```
`DISP ALL;` は全要素を表示 (66 行ごとの停止を省く)。`LINE["LENGTH"]` は要素 (成分) の個数 ($$$ を含む)。全文は `.sad`。

## 結果の読み方
- TWICE: 10 行 (M、Q.1…D.4、$$$)、`s(m)` の最終 5.2 m = 2×(0.3+1+0.3+1)。
- REVL: 要素名の前に `-` が付き (`-D.1`、`-QQ` …)、並びが逆 (D QQ D Q の順)、`LINE["DIR"]` = {1,-1,-1,-1,-1,1}。
- NEST: 展開後の名前は `M Q.1 D.1 QQ.1 D.2 D.3 QQ.2 D.4 Q.2 Q.3 Q.4 $$$`、全長 5.8 m (HALF 2.6 + 逆向き HALF 2.6 + 2×0.3 = 5.8)。

## 物理・解析式
なし (ラインの展開規則そのもの。長さは足し算)。

## gs と pysad の比較
全桁一致 (表・`Print` の値とも)。反転配置の `-` の付き方 (`-QQ` のように ".1" が省かれる場合を含む) も同じ。

## 未確認・限界
- 要素ではなく LINE への繰り返し `n*(...)` の入れ子の括弧、`BeamLine[...]` オブジェクトによる動的なライン作成 (SADHelp の BeamLine の項) は別の例 (fn_Element_SetElement は要素だけ) でも確認していない。
- `INS` の初期値のため、先頭の MARK に BX=BY=1 を与えている (CELL では周期解が不安定な短いラインなので)。
