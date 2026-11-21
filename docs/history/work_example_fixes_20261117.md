# example 10〜13・PF T5 の gs と pysad の差の修正 (2026-11-17)

元の一覧は `local/wake/example_findings.md` の「10-13」と「ASTRID2/T5/T6」の項。基準は gs の実動作 (`local/wake/gsrun.sh`)。pysad の確認は `PYSAD_SAD_PACKAGES` を付けて行う
(付けないと GUI の Tkinter が読まれず、`TclEval[...]` が式のまま残って、別の不具合に見える)。試験は `pysad/tests/test_example_fixes.py` (14 件、`mkbump_ring.sad` を使う)、`test_tclbridge.py` に 3 件追加。

## 直したもの (gs と比べて一致を確認)
| 項目 | 原因 | 直し方 (ファイル) |
|---|---|---|
| (a) MakeBump の `list assignment` エラー | `b=Bump[…]` は bump.n の `(a_=Bump[…])^:=…` (Set の UpValue)。(1) `Bump` が自動読み込みの表に無かった (`Name[..]:=` の形でないため)。(2) `_protect_builtins` が、組み込みと同名の `Set` の定義 (= この UpValue) を読み込み後に捨てていた | `pkgload.py` (`_UPDEF`、Set 系の UpValue を残す)、`sadscript._set` (右辺・左辺の記号が読み込みの引き金のとき先に読む)。ybump=.0010000141、12 本で gs と一致 |
| (b) Couple の ElementValues | (1) `"K1"["QD"]` (文字列が頭部の式) が文字列の添字として評価され int('QD') で落ちた。(2) `With[{y},Hold[y]/.…]`・`With[{f},f&]`・`Hold[…]/.Hold->Null`・`EVCompiled[[1]]` (Hold の Part) の評価が gs と違った。(3) `Element$` が無かった。(4) CALC の前に SetCoupledElements を呼ばなかった。(5) USE で ElementValues が消えなかった | `sadscript.py` (文字列頭部の式は式のまま、With が Hold・純関数の本体に値を埋め込む、ReplaceAll の Hold・Func・Null、Part の評価し直し、`subst_slots` は Hold の中を置換しない)、`sadelem.py` (`Element$`)、`ffscmd.py` (`_apply_element_values`、USE で消す)。QD の K1 = -0.56、NX=6.788257 で gs と一致 |
| (c) Rational が和・積を評価しない | rational.n の `Plus[Rational[…],…]^:=…` などが、`_protect_builtins` で捨てられていた。`a-b`・`a/b` は Subtract・Divide で、規則に届かない | `pkgload.py` (Plus・Times・Power・Minus・Divide・Subtract の UpValue を残す)、`sadscript.py` (Rational を含む Subtract・Divide は Plus・Times へ回す)。13 式で gs と一致 (`1/3+1/6` は `Rational[3,6]`。gs も同じ) |
| (d) NIntegrate の複素数 | 被積分関数の値を float にしていた | `fitroot.py` (complex のまま。スプラインは線形) |
| (e) LogGamma1 | 未実装 | `specfn.py` (`LogGamma1[x]=LogGamma[x+1]`、gs と 1e-15) |
| (f) ParametricPlot の範囲に Pi | HoldAll の引数が `2*Pi` のままで `xmax_Real` に照合しなかった。曲線のリストは `Evaluate[ReleaseHold[a]]` が保持引数位置で評価されなかった | `sadscript._match_held` (定数式 = Pi・E・Degree と数値関数は値で照合)、保持引数の `Evaluate[x]` は評価 |
| (g) ListContourPlot の Tcl エラー、SADHelp の `Second` | (1) `Position[…,_?(Second[#]>3&),1,1]` の検査の中のエラー (Second[List]) で pysad は文が落ち、gs は不一致として続ける。(2) Tcl の `tagOrId "1" doesn't match any items` (`lower 85 1`): K64 の `TkCanvasLastItem` は Tk の nextId-1 で `delete all` のあとも id が続くが、pysad は `find all` の最大値 (全消去後は 0) を返し、Canvas.n の `StartID=CanvasLastID+1` が古い 1 になっていた | `sadscript._match_one`・`_match_held` (PatternTest の中の SADError は不一致)、`tkcanvas.last_item` (Canvas ごと (パス名と winfo id) に見た最大の id を覚える)。項目数は gs と同じ (68、132) |
| (h) OptimizeOptics Output->0 | `DownhillSimplex[p,f,Join[{opt},{VariableRange->limit}]]` は規則のリストを 1 引数で渡すのに、選択肢を平らに読まず Output 既定 6・MaxIteration 無視だった | `miscfn.py` (入れ子のリストを読む) |
| (i) PointSymbol の名前の文字 | `{"CI","CL","1O"}//""` を pysad は要素を連結 (`CICL1O`)、gs は `{"CI","CL","1O"}` の形の文字列 (StringJoin はリストを連結しない) | `sadscript.py` (StringJoin)。Canvas 項目 30 で gs と一致 |
| (j) Entry/ListBox/Scale のメソッド (Widget 値) | `Tkinter.n` の `(a$_=w[opt$___])^:=(a$=w[Widget[Root$Window,Null],opt$])` が次の規則 `(a$$_=w[Widget[master$_,_],opt$___])^:=With[{a$=a$$,…},…a$:=Widget[a$,…]]` を呼ぶ入れ子で、`With` が呼び出し元の規則の引数 `a$` (ユーザーの w) を見て、`a$` でなく `w` を束縛していた。本体の `a$:=Widget[…]` が別名に付き、w に値が付かなかった | `sadscript._with` (別の規則の呼び出し = boundary の frame で探索を止める)。`Head[e]=Widget`、`e[Insert]`・`lb[Insert]` が gs と同じ (abcXY、2 項目) |
| (k) calc_traj_03 の窓 | (1) `Canvas$Widget = c1` の UpValue は、左辺の記号が最初の使用で読み込まれ、かつ UpValue の照合が遅延した定義 (`The$CanvasObject := …`) を持つ記号を評価して既定の Canvas を作っていた (tclerror 1501 の原因)。(2) `OpticsPlot` の自前データ: `OpticsList1$[pr_,{f_}?(Length[Dimensions[#]]>=3&),opt___]` が `OpticsList1$[{pr1_,pr2_},f_,opt___]` より後に試されていた (pysad の特殊さの点数が、検査つきの構造パターンを低く見ていた)。`Legend->False` を MapThread に渡すと「List of Lists」 | `sadscript.py` (`_match_one` の遅延定義の記号、`_spec_score` の PatternTest)、`sadlist.py` (MapThread・ScanThread はリストでない要素を繰り返す)。canvas items は 114 と 115 で gs と同じ |

## 直さなかった・残る差
- `OpticsPlot` の text 項目の数 (FODO 24 セルの例で gs 69、pysad 74): `RemoveOverlap` が使う文字の外接矩形が、Tk のフォント (gs は X11、pysad は Windows) で決まる。フォントの測り方を合わせないと直らないので未対応。
- 最初の `Window[]` の親: gs は Toplevel (`.5`)、`.` に K64 起動時の「Return to SAD」ボタン。pysad は `Window[]` が `.` 自身 (木の段、`TextVariable` の名前の番号が違う)。起動時の窓の再現が要るので未対応。
- `LINE["POSITION",{{1,2},{3,4}}]` (入れ子のリスト) は gs が入れ子のまま返し、pysad は平らにする。`LINE["POSITION",0.5]`・`0.0` は gs が 244 (nlat+1)、pysad は 1。OpticsPlot の自前データには影響しない (x は整数の要素番号) ので未修正。
- `Hold[x]` の中の関数の仮引数 (HoldAll の `f[x_]:=Hold[x]` の x) は gs が値を埋め込む (`Hold[other,y]`)。pysad は `x` のまま (`With` の変数は埋め込む)。`Hold[1+2]` は gs が `Hold[3]` (算術を評価する) で pysad は `Hold[1+2]`。
- `OptimizeOptics` の先頭の `{{`,System`},System`}` (BeginPackage の戻り値) は文脈を持たない pysad で `{$ContextPath,"`"}`。
- BesselY_JZero・Bessel・HermiteH・Ei の gs の誤値は gs 側の問題なので直さない。BeamPlot (gs は SIGXX=0 で楕円なし) は未対応のまま。
- ElementValues の VISIT での退避・BYE での復元、フィット変数との連動 (`tffssetupcouple` の `CouplingMatrix`・`CoupledKeys`) は未対応。

## 確認の仕方
`PYSAD_SAD_PACKAGES=<SAD の Packages のフォルダ>` を付けて全試験 (旧式 3 本を除く) と回帰 (`QUICK=1`) を回した。全試験は 575 通過・3 skip (578 件収集、3 分 21 秒)。回帰の結果は下。
直した例の `.pysad.txt` は `example/tools/run_example.py` で再実行し、`.md` の「gs と pysad の比較」を更新した (10_math_data 2 本、11_graphics 5 本、12_gui_tkinter 6 本の出力と一部の説明、13_packages 4 本、PF の calc_traj_03)。図 (`*_plot.png`) は作り直していない (11/16 のまま)。

## 結果 (2026-11-17)
- 全試験 (旧式 3 本を除く、`PYSAD_SAD_PACKAGES` つき): 575 通過・3 skip (578 件収集)。
- 回帰 `QUICK=1 run_regression.py`: 50 ベンチ、FAIL 0 (所要 125 s)。
