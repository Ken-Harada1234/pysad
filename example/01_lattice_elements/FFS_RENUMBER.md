# FFS_RENUMBER (RENUMBER)
種別: FFS コマンド。フォルダ: 01_lattice_elements。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`RENUMBER comp` は、同名の要素につく番号 (`QF.1`, `QF.2`, …) を `comp` から数え直す (SADHelp: "RENUM comp renumbers the component number starting from a component comp")。
周期的なリングで、ある要素 (例: 入射点に近い QF) を `.1` にしたいときに使う。

## スクリプト
```
USE RING; CELL; CALC;
Print["before: ", LINE["NAME",{5,6,7,9,13,14}]];
RENUMBER QF.3;
Print["after : ", LINE["NAME",{5,6,7,9,13,14}]];
```

## 結果の読み方
| 位置 | 5 | 6 | 7 | 9 | 13 | 14 |
|---|---|---|---|---|---|---|
| 前 | QF.1 | D2.1 | B.1 | QD.1 | SF.2 | DS.3 |
| gs の後 | QF.23 | D2.45 | B.23 | QD.23 | SF.23 | DS.46 |
| pysad の後 | QF.23 | D2.45 | B.23 | QD.23 | SF.23 | DS.46 |

gs: `QF.3` を起点に番号が振り直され (QF.3 → QF.1)、それより前の要素は末尾側の番号 (QF.1 → QF.23 …) になる。`D2.45` のように、同名で 2 個ずつ現れる要素 (D2、DS) も数え直されている。

## 物理・解析式
なし。

## gs と pysad の比較
**○ 一致 (2026-11-17 に対応)。** `RENUMBER comp` (略 `RENUM`) は comp を起点に周回順で同名の要素の番号を 1 から振り直す (gs の `QF.23 D2.45 B.23 QD.23 SF.23 DS.46` と同じ)。`LINE["MULT"]` も振り直した番号になる。起点が無い (`RENUMBER;`、存在しない名前) と ` ???-FFS-Error-Missing origin component for RENUM_BER ` を出す。

## 未確認・限界
- `RENUMBER QF` (名前だけ) は最初の QF が起点なので変化しない (gs と一致)。要素の INDX 属性による番号の初期値は未対応。
- 番号の振り直し後に、`QF.1` を使う FIT の位置指定が新しい番号を参照するかは未確認 (`Twiss`・`LINE`・`DISP` の位置指定は新しい番号)。
