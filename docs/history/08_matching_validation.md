# 08_matching: gs (K64) との比較・検証の詳細 (本文から移したもの)

本文 `docs/08_matching.md` から、gs との比較の細部を元の文脈つきでそのまま移した。

## 元の節: 幾何の条件 (GX GY GZ CHI1–3) と tdgeo

QF K1・QD K1 の 2 変数で NX・NY を 0.2 に fit した 3 通りが gs と全桁一致: 既定 K=−.0068645845/−.0057261776、UNIPOL K=.3289522375/−.3239879365、RELW K=−.0024047842/−.0107636499 ([`test_item4.py`](../pysad/tests/test_item4.py) 3 件)。この 2 変数の小さな例だけの確認で、

## 元の節: 幾何の条件 (GX GY GZ CHI1–3) と tdgeo

- **DELCOR で OCT を対象に含める**: チューンと QF の DX が gs と 12 桁で一致 (OC 自体の DX は gs も pysad も 0 を返すので、OCT の効果そのものを見る強い比較ではない)。

## 元の節: 7. 照合結果

## 7. 照合結果

| ベンチ | 内容 | 結果 |
|---|---|---|
| 13A-1 (小リング) | fit NX 7.55 NY 7.18, free QF QD | 2 反復 NEWTON、反復行 `2 0 1.6776E-06 (NEWTON) 1.000 2` が同一、Matched 5.0763e-16 vs SAD 5.0759e-16、K1 ≤5e-11 |
| 13A-2 | + ORG BX 17.5 BY 1.0, free QF QD QF2 (解なし) | SAD も収束せず " Residual = 2.5977E-03" で打ち切り。pysad も同じ経路 (qcod overflow → nstab) で残差・K1 (≤4.4e-11) が一致 |
| 13B (実ラティス 試験用のラティス, 7 ファミリー) | fit NX 47.6 NY 17.6, MKC BX/AX/BY/AY | 2 反復 Matched 3.3724e-10 (両者)、K1 28 台 ≤1.1e-10 |
| 不安定ライン | 不安定な試行点での直線探索・crate=0 | fmincube の γ と `_log10` の扱いで SAD の経路を再現 (docs/history/work_unstable.md) |

解が一意な問題 (13A/13B) では数値微分でも SAD と同じ解・同じ反復数になる。多解の問題で SAD と同じ経路を取るには解析微分 (既定の match6d、3 節) が要る。

