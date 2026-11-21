# UND (アンジュレータ)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `undulator.f` (K.Ohmi と E.Forest、KEK Report 92-14)。[`pysad/und.py`](../pysad/und.py)。

```
UND U1 = (L=1 KX=0.1 QY=0.2 POLE=10 SLICE=4 PHI=30);              ! 磁場なし (k64 と v2 の両方で書ける)
UND U2 = (L=1 FBX=1.8 POLE=20 SLICE=2000);                         ! 磁場つき (v2 だけ)
```

## モードによる違い

**K64 の UND は磁場を設定できない** (docs/21 B37): `initb1.f` は `kytbl(kwFBX,icUND)=ky_FBX_UND` と書くが、`kwFBX` は未定義の変数 (=0) で、`FBX` `FBY` `BX` `BY` はどれも "Undefined keyword" になる。そのため F0 = G0 = 0 で、undulator は磁場なしの (正準運動量の近似の) ドリフトを刻みごとに行う。
さらに `undulator.f` は陰的な刻みの右辺に dpds を二重に掛けている (docs/21 B78)。磁場があっても収束力が刻み幅に比例して消える (200 刻みと 2000 刻みで py が 10 倍違う)。

| | k64 モード (K64 と同じ) | v2 モード (既定) |
| --- | --- | --- |
| `FBX`・`FBY` | NotImplementedError | 受ける (磁場つきの UND) |
| 追跡 | 磁場なしのドリフト | 磁場つきのカーネル (`und.und_seg1`。dpds は 1 回だけ掛ける) |
| Twiss・Emittance[]・CALC6D | 恒等 (長さは幾何にだけ入る) | 追跡の写像の数値微分が入る。**磁場がなくても長さ L のドリフトとして光学に入る** (追跡と揃えた。`SLICE=0` は恒等) |
| 放射 (RAD・FLUC) | なし | 1 周期を 64 個以上に分けて `tradke`。RAD で放射する |

磁場なしの UND の追跡は両モードで同じ (gs と一致)。光学が違うのは「磁場なしで v2 はドリフト、k64 は恒等」の 1 点で、v2 のこの判断は要判断のまま (K64 の「光学では恒等」を仕様と見るなら、磁場つきのときだけにする選択もある)。

| キーワード | 意味 |
| --- | --- |
| L | 長さ |
| FBX | v2: 垂直磁場のピーク [T] (x 方向に蛇行)。軸上の By = (FBX·c/p) cos(kz s) |
| FBY | v2: 水平磁場のピーク [T]。Bx = (FBY·c/p) cos(kz s + PHI) |
| KX、QY | 横方向の波数 (kx = KX、ky = √(kz²+KX²)、qy = QY、qx = √(kz²+QY²)) |
| PHI | 位相 [度] |
| SLICE | 刻み数 (ds = L/SLICE) |
| POLE | 極の数 (kz = 2π POLE / L)。0 だと kz=0 で NaN になる (K64 も同じ) |

## 使われる所

| 計算 | k64 モード | v2 モード |
| --- | --- | --- |
| 追跡 (`TrackParticles`、MAIN の TRACK、DA。`track6n`) | `K_UND` (`und.k_und`)。入口・出口の TENT/TEXIT は dx=dy=θ=0 で何もしない | 磁場つきのカーネル |
| Twiss・Emittance[]・EMIT・CALC6D・Touschek | 恒等 (K64 の `qtwiss`・`tturne` は `case default` で何もしない) | 追跡のカーネルの数値微分。放射は 1 周期を 64 個以上に分けた `tradke` |
| 幾何 | `geo` が `UND` の長さだけ直進 | 同じ |
| スピン追跡 (track6 の経路。POL + RAD) | NotImplementedError | 刻みごとに [写像 → 放射とスピンの歳差]。縦磁場 Bs (軸外) と、周期の途中で終わる端のポテンシャルの飛び (薄い蹴り) の歳差も入る |

## 検証

- k64 モード (gs との比較、[`pysad/tests/test_und.py`](../pysad/tests/test_und.py) 4 本): 追跡 (20 粒子、3 周) 相対 0、INS の Twiss と幾何 2e-10 以下、Emittance[] 1e-12 台。
- v2 モード ([`pysad/tests/test_v2_C_und.py`](../pysad/tests/test_v2_C_und.py) 9 本。gs は磁場を設定できないので比べていない。平面アンジュレータの解析式で確かめた): 軸上の軌道 (出口で x = x' = 0)、行路長 −(F0/kz)² L/4 (1e-6)、垂直の収束 k = F0²/2 (行列要素 5e-4)、KX≠0 の kx・ky (5e-3)、シンプレクティック性 1e-7、放射損失 (2/3) r_e γ³ (F0²/2) L (1e-4)。リング: Δνy = ∫βk ds/4π (1%)、ΔU0 (2e-3)、減衰率の和 = 2U0/E (1e-3)、量子励起は同じ ⟨1/ρ²⟩ の BEND の列との比が ⟨|1/ρ|³⟩ の比 (2%)。
- v2 のスピン追跡 ([`test_v2_C2_spin.py`](../pysad/tests/test_v2_C2_spin.py)): B = rot A の中の Thomas-BMT の積分と 2e-6 (1 T、1 GeV、POLE = 0.25、5、5.25、中点面の外 2〜3 mm)。
- 例: [`example/v2_C/C1_undulator.sad`](../example/v2_C/C1_undulator.sad)。

## 使える範囲と条件 (v2)

- `SLICE` は `POLE` の 20 倍以上 (刻みは 1 次の陰的な蛙跳び。刻みが周期の整数倍だと sin(kz s) = 0 で磁場が見えない)。試験は 3 GeV (スピンは 1 GeV)、1 T・1.8 T、λ = 5 cm、20 周期。
- 実機の大きなリング・挿入光源のあるリングでは確かめていない。

## 限界

- k64 モード: 磁場がないので、アンジュレータの放射・軌道への効果は計算されない (K64 と同じ)。
- v2 モード: `BX`・`BY` という名前は K64 にも無い (ValueError)。端の磁極の整形は無い (軌道は x' = −(F0/kz) sin(kz s) で、中心が −F0/kz² ずれる)。PHI≠0 の入口では正準運動量をそのままつなぐので、ベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま)。
- track6 (低速追跡) は、k64 モードでは UND を扱わない。v2 では POL + RAD のスピン追跡の経路で扱う。
