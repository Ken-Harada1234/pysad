# FFS_REFERENCE (REFERENCE: 参照光学)
種別: FFS コマンド。フォルダ: 02_optics_cod。ラティス: `00_lattice/FODO_ring.sad`。エネルギー 3 GeV。フラグ: CELL。

## 何を見る例か
`REFERENCE` (略 `REF`) は、現在の光学を**参照光学**として記憶する (`USE` 後の最初の `CALC`/`GO` でも自動で設定される)。以後の変更との差を、`Twiss["R関数", 位置]` (参照の値)、`Twiss["D関数", 位置]` (現在 − 参照)、`DISP REFERENCE` (参照の表)、`DISP DREFERENCE` ((β−βR)/βR などの表) で見られる。ここでは QF の K1 を 0.7→0.75 に変えて、β の変化を見る。

## スクリプト
```
USE RING; CELL; CALC;     Print["BX at 5,50 : ", Twiss["BX",{5,50}]];
REFERENCE;   QF K1 0.75;   CALC;
Print["RBX (reference) : ", Twiss["RBX",{5,50}]];   Print["DBX (BX-RBX)    : ", Twiss["DBX",{5,50}]];
DISP REFERENCE 1 8;      DISP DREFERENCE 1 8;
```

## 結果の読み方 (gs)
変更前の BX (位置 5 = QF.1 入口、位置 50): {6.42693, 1.30734}。変更後: {6.68626, 1.12741}。`RBX` は変更前と同じ値 {6.42693, 1.30734}、`DBX` = 変更後 − 変更前 = {0.040352, -0.13763}。
`DISP REFERENCE` の表は列名に `R` が付き (`AXR BXR NXR …`)、`DISP DREFERENCE` は `dAX dBX/BX dNX …` (β は相対差)。

## 物理・解析式
差 `DBX` は定義どおり BX − RBX: gs で 6.68626 − 6.42693 = 0.259337 (表の `DBX`: 0.040352)、1.12741 − 1.30734 = -0.179929 (表: -0.13763)。

## gs と pysad の比較
**○ 一致 (2026-11-17 に対応)。** `REFERENCE` (略 `REF`) は直前の CALC の光学を参照光学にする。`Twiss["RBX",…]` (参照)、`Twiss["DBX",…]` (現在 - 参照、BX BY BZ は比)、`Twiss["RALL"/"DALL",…]`、`DISP REFERENCE`・`DREFERENCE` (見出しの文字の置換も含む) が gs と全桁一致。参照光学の位置は整数位置だけ (小数位置はエラー)。REFERENCE 前の参照は 0。

## 未確認・限界
- `BMAGX` などの BMAG 係数 (`Twiss["BMAGX",...]`、参照光学とのミスマッチ) は gs でも出していない。
- `DRAW`、`OpticsPlot` での参照光学の使用は 11_graphics の担当。
