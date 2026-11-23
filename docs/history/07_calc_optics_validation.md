# 07_calc_optics: gs (K64) との比較・検証の詳細 (本文から移したもの)

本文 `docs/07_calc_optics.md` から、gs との比較の細部を元の文脈つきでそのまま移した。

## 元の節: INS (cell=.false.) と twmov (詳細: `history/work_ins_offset.md`)

ベンチ: `bench_level22a..g.sad` (SAD と 1e-9〜1e-10)。

## 元の節: INS (cell=.false.) と twmov (詳細: `history/work_ins_offset.md`)

`bench_level30a.sad` (実ラティス、5 通り: 4D の CELL → INS、6D の CELL → INS、MARK にオフセット、初期 Twiss・BZ・DDP を変更、RFSW RADCOD) で SAD と ≤1.6e-10。

## 元の節: 2. 要素ごとの CALC 用の写像 (`qtwiss.f:184-400` の分岐)

。CALC・INS・CELL・RAD ON・DP≠0 でチューン 1e-15、領域内の Twiss 1e-11 (gs と一致、[`test_sol_optics.py`](../pysad/tests/test_sol_optics.py))

## 元の節: 行列の作り方: 解析行列と数値 Jacobian

これで 14A の DP0 スキャンのチューンが 1e-10、3A のチューンが 1.5e-10、3B の EY が SAD と一致する。

## 元の節: 行列の作り方: 解析行列と数値 Jacobian

    SAD との差: 22e の CELL (DP0=0.004、四極ずれ、閉軌道≠0) **1.4e-11**、23a 1.7e-11、25d 4.5e-11。FRINGE=1〜3 の四極を含む `Emittance[]` (bench_level26b) は 1.3e-9 (残りは SAD の tcod 収束の差。数値微分だとチューンが 1.1e-5 違う)。

## 元の節: 4b. CALC6D と新しい要素の Twiss の範囲

gs と一致: 

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

gs と終端の BX AX BY AY EX EPX NX NY が全桁、

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

gs と同じ向きで、入口に分散のある加速区間の EX EPX は gs と全桁 (薄い空洞・厚い空洞。`test_dispersion_through_accelerating_cavity_matches_gs`)

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

 (pysad k64 2.73475252978794)

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

 (1e-10)、β も 1e-9

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

 (gs と 8e-8)

## 元の節: v2 モードで変わった CALC・光学の動作 (k64 モードは K64 のまま)

- **`Twiss["DETR"]`** (両モード): INS の CALC の後の値は `ffscmd._rows_from_tws` が `detr` を渡して得る。小さな結合リング (45 度回した四極 K1=0.8)、大きなリング A (QDA2 を 0.03 rad 回し空洞に DX=2e-4・DY=−3e-4: `Twiss["DETR", {1, 500, 1000, 2000}]` と R1 が gs と 10 桁)、PF リングで gs と全桁。DETR = R1R4 − R2R3 の恒等式が 1e-12。`MEASURE` の初期分布が読むのは先頭 MARK の生の値 (CELL CALC が書く) で、正しい。**DETR ≥ 1 (モードが入れ替わる) の INS と MEASURE は未確認** (DETR < 1 のリングは R1 = 0.3 まで入れ替わらなかった。入口を DETR > 1 の生の値で与える輸送路は gs と pysad で入口の規格化が違い比較できなかった)。

## 元の節: 5. ベンチマークの要約 (CALC 関連)

## 5. ベンチマークの要約 (CALC 関連)

| Level | 内容 | SAD との差 |
|---|---|---|
| 1 | 20 セル FODO: 閉軌道 (全 QD に DX=1e-4)、放射積分 I1–I5 | 閉軌道 入口 x −110.8245 µm 一致、I1–I5 表示桁一致 |
| 2 | オフモメンタム Twiss (dp=0, ±0.01, ±0.02)、複合機能 BEND | βx, αx, ηx ≤3e-8 / βy, Qx, Qy ≤6e-11 |
| 3A | 6 極と色収差 (DP0 スキャン 21 点 + 3 次当てはめ) | チューン ≤2.6e-9、ξ 相対 ≤3e-7 |
| 3B | 6 極を通る COD (DELX QD.1, DELY SF.1) | 閉軌道 9–10 桁、β/α/η/Q ≤2.4e-8、EY 8e-10、R1–R4 5–6 桁 |
| 3C | x–y 結合 (DTHETA QF.3, B.5) | 閉軌道 8–9 桁、光学 ≤1.6e-8 (EX/EPX は 3e-11)、R1–R4 7 桁 |
| 3D | 強さの誤差 (DK, DDK) | 閉軌道 9–10 桁、光学 ≤9e-8 |
| 5 | INS 輸送ライン、F1/FRINGE/DISFRIN のフリンジ | Twiss ≤2e-8、軌道 ≤1.3e-11、位相 ≤5e-11 |
| 12 | 実ラティス 20 セル/試験用のラティス (3563 要素、NLAT=3564) の CELL CALC | Twiss・チューン ≤4e-11 |
| 13/14 | マッチング後の Twiss、DP0 スキャン | K1 ≤1.1e-10、チューン ≤2e-8 |
| 28 | 不安定ライン (QD/QF を ×1.19, ×0.9, …): CELL CALC の続行、`Unstable.`、GO の nstab、INTRES/HALFRES | Twiss ≤2.5e-10 (β≲1e5)、GO の反復表・K1 一致 (docs/history/work_unstable.md) |

