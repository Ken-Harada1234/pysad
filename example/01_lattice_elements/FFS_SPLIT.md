# FFS_SPLIT (SPLIT と MARK の OFFSET)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティス: `FODO_ring.sad` と同じ光学のリングに、QF の中心を指す `MARK PQF = (OFFSET=-0.5)` を加えたもの (ファイル内、RF なし)。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
SADHelp の `SPLIT comp length` は要素を 2 つに割るコマンドだが、**この K64 ビルドの gs は `???-FFS-Error-SPLIT is obsolete. Use OFFSET of a marker.` と答える** (廃止)。
代わりに、MARK の `OFFSET` (相対位置。小数で要素の途中を指せる) を使うと、QF の中心の光学関数を MARK の位置から読める。

## スクリプト
```
MARK PQF = (OFFSET=-0.5);
LINE CELL = (SF DS QF PQF D2 B D2 QD DS SD D1);
...
SPLIT QF.1 0.1;                     ! gs: obsolete (エラー)、pysad: 未定義のコマンド
p = LINE["POSITION","QF.1"];
Print["BX at QF.1 entrance, QF.1 + 0.5, QF.1 exit = ", Twiss["BX",{p, p+0.5, p+1}]];
Print["BX at PQF.1 (OFFSET=-0.5)                    = ", Twiss["BX","PQF.1"]];
```
OFFSET=-0.5 は「この MARK の位置より 0.5 要素分手前」= QF の中心。`Twiss[関数, 位置]` の位置に小数を与えると要素の途中の値が補間される。

## 結果の読み方
`BX at QF.1 ...`: QF の入口 6.426926、中心 6.846232 (最大値付近)、`p+1` (= MARK PQF の位置) も 6.846232 となる (MARK が OFFSET の位置の値を持つため)。`BX at PQF.1` も同じ 6.846232。

## 物理・解析式
なし (補間そのもの)。参考に、QF の中心では $\alpha_x\approx0$ で β が極大 (6.846 m)。

## gs と pysad の比較
- `SPLIT`: gs はエラー (obsolete)、pysad は `?Undefined command or element: SPLIT`。**どちらも動かない** (メッセージが違うだけ)。
- OFFSET の MARK の Twiss は一致 (BX = 6.426926 / 6.846232 / 6.846232 / 6.846232)。
- **`LINE["S", p+0.5]` が違う**: gs は 0.375、pysad は 0.45。QF の入口 s=0.3 m、長さ 0.3 m なので中心は 0.45 m が自然 (pysad と、gs の `LINE["S","PQF.1"]`=0.45)。gs が 0.375 を返す理由は不明 (未調査)。

## 未確認・限界
- 要素を実際に 2 つに割る方法 (SADHelp の SPLIT の記述) は、この gs では動かないため確認できない。
- OFFSET が先頭・末尾の MARK にあるときの動作 (周期条件がずれる、SADHelp の OFFSET の項) は未確認。
