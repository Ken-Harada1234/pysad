# TCAVI (横偏向空洞・クラブ空洞)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `ttcav.f` (追跡)、`ttcave.f` (6 次元の写像と行列)。[`pysad/tcav.py`](../pysad/tcav.py)。

```
TCAVI TC = (L=0.4 K0=1e-4 FREQ=4.5e8 PHI=0.1 DX=2e-4 DY=-1e-4 ROTATE=0.05);
```

| キーワード | 意味 |
| --- | --- |
| K0 | 偏向の強さ (ak、kick = ak sin(位相)) |
| FREQ / HARM | RF の周波数。HARM≠0 は ω0·HARM/c。HARM の小数部は周回ごとに位相をずらす (追跡) |
| PHI | 位相 (電荷の符号を掛ける) |
| DX DY ROTATE | ずれと回転 (入口と出口の TENT/TEXIT) |
| RAD | 0 (既定) で放射する (追跡は RAD と FLUC の両方が要る) |
| RANKICK RANPHASE | 乱数の揺らぎ (tgauss)。**追跡は未対応** (NotImplementedError) |

RFSW OFF (CALC の既定) では長さ L のドリフト (tturn、ttcave とも)。

## 使われる所

| 計算 | 実装 |
| --- | --- |
| 追跡 | `track6n` の `K_TCAV` (`tcav.apply_tcav`)。TENT/TEXIT の回転は K64 の式の順 (更新した x を使う) のまま |
| Emittance[]、CALC6D、Touschek | `emit6d._ttcave` (tchge → [tsetr0] → 半ドリフト → キック (解析行列) → 半ドリフト [tradke] → tchge) |
| Twiss (CALC) | RFSW OFF はドリフト。**RFSW ON の CALC は未対応** (NotImplementedError。Twiss が momentum を持たないため) |

## 検証 (`pysad/tests/test_tcav.py` 8 本)

- 追跡 (gs と比較): RFSW ON、全周 3 周、200 粒子。FREQ、ずれ・回転付き、HARM の小数部 (周回で位相が動く)、L=0、すべて相対 0。RAD FLUC (20 粒子、乱数列まで) も相対 0。
- Emittance[]: L=0 の TCAVI で gs と一致 (エミッタンス 1e-11、チューン・運動量収縮因子・COD 1e-10 台)。ずれ・回転付きは gs の EMT が NaN なので、COD・チューン・運動量収縮因子だけ比べた。
- 追跡と 6 次元の写像が同じ写像であること (二つの直訳の整合): 1 粒子で 1e-19、キックの解析行列は数値微分と 4e-12 で一致。

## K64 の癖・不具合

- **L≠0 の TCAVI を含むリングの Emittance[] は gs が壊れる** (docs/21 B36)。K0=0 でも EnergyLossU0 が NaN、縦のチューンが 0、エミッタンスが桁違い (K0=1e-6 で 1.3e-13 など) になる。ドリフトに置き換えると正常。`ttcave` の 2 番目のドリフトが `tdrife0` に放射の旗 `krad` を渡し、長さ 0 の放射 (alr=0) で `tradke` を呼ぶことが疑わしいが、RAD=1 でも NaN が出るので、原因は特定できていない。pysad は gs の NaN を再現せず、物理として自然な値を出す (ドリフトと同じ結果から K0 の効果が加わる)。L=0 は gs と一致する。
- 追跡は RAD と FLUC の両方 (`rad .and. RAD_TCAV==0`) のときだけ放射する。Emittance[] は `RAD_TCAV==0` だけを見る (NORAD でも TCAVI は放射の旗が立つ)。pysad は Emittance[] の `enarad` に従う。

## 限界

- RANKICK・RANPHASE (追跡)。
- RFSW ON の Twiss (CALC、INS)。
- wake との併用、SOL 領域内の TCAVI (track6n の領域外の要素のみ対応) は確認していない。
