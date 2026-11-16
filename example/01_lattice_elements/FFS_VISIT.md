# FFS_VISIT (VISIT と BYE)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティスは `FODO_ring.sad` と同じ光学 + 1 セルだけのライン `ONECELL` (ファイル内、RF なし)。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`VISIT ライン` は、いまのラインの状態 (マッチング条件など) を保存したまま別のラインに移る。`BYE` で元のラインに戻る。`USE` は状態を捨てるのと対照的。
リングの途中で 1 セルだけの光学を見て戻る、という使い方。

## スクリプト
```
USE RING; CELL; CALC;      Print["RING     : NX=", ..., " nlat=", LINE["LENGTH"], " name=", BeamLineName[]];
VISIT ONECELL; CELL; CALC;  Print["ONECELL  : NX=", ...];
BYE;                        Print["after BYE: nlat=", ...];
CALC;                       Print["RING again: NX=", ...];
```
`BeamLineName[]` は今のラインの名前。

## 結果の読み方
RING: NX=6.696329 NY=4.906056、nlat=242。VISIT 中: ONECELL (nlat=12): NX=0.279014、NY=0.204419。`BYE` のあと nlat=242、name=RING、再計算で NX=6.696329。

## 物理・解析式
24 セルが同一なので 1 セルの位相進みはリング全体の 1/24:
$$\nu_{\rm cell}=\nu/24:\quad 0.2790137\ (x),\ 0.204419\ (y)$$
| 量 | 式の値 | gs | pysad |
|---|---|---|---|
| NX (1 セル) | 0.279014 | 0.279014 | 0.279014 |
| NY (1 セル) | 0.204419 | 0.204419 | 0.204419 |


## gs と pysad の比較
全桁一致。

## 未確認・限界
- VISIT 中に FIT 条件を設定して BYE で元に戻る、といった保存の動作は未確認 (docs/03 の `_visit_stack`)。
