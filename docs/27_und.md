# UND (アンジュレータ)

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

K64 の `undulator.f` (K.Ohmi と E.Forest、KEK Report 92-14)。`pysad/und.py`。追跡だけ。

```
UND U1 = (L=1 KX=0.1 QY=0.2 POLE=10 SLICE=4 PHI=30);
```

**K64 の UND は磁場を設定できない** (docs/21 B37): `initb1.f` は `kytbl(kwFBX,icUND)=ky_FBX_UND` と書くが、`kwFBX` は未定義の変数 (=0) で、`FBX` `FBY` `BX` `BY` はどれも "Undefined keyword" になる。そのため F0 = G0 = 0 で、undulator は磁場なしの (正準運動量の近似の) ドリフトを刻みごとに行う。直訳は磁場の項を残してあるが、利用者は設定できない。pysad も FBX 等を受けない (NotImplementedError)。

| キーワード | 意味 |
| --- | --- |
| L | 長さ |
| KX、QY | 横方向の波数 (kx = KX、ky = √(kz²+KX²)、qy = QY、qx = √(kz²+QY²)) |
| PHI | 位相 [度] |
| SLICE | 刻み数 (ds = L/SLICE) |
| POLE | 極の数 (kz = 2π POLE / L)。0 だと kz=0 で NaN になる (K64 も同じ) |

## 使われる所

| 計算 | 扱い |
| --- | --- |
| 追跡 | `track6n` の `K_UND` (`und.k_und`)。入口・出口の TENT/TEXIT は dx=dy=θ=0 で何もしない |
| Twiss・Emittance[]・CALC6D | **恒等** (K64 の `qtwiss`・`tturne` は `case default` で何もしない。長さは幾何 (s、LINE["GX"]) にだけ入る) |
| 幾何 | `geo` が `UND` の長さだけ直進 |

## 検証 (gs との比較、`pysad/tests/test_und.py` 4 本)

- 追跡 (20 粒子、3 周): 相対 0。
- INS の Twiss と幾何: 2e-10 以下 (恒等であることの確認を含む)。
- Emittance[]: 1e-12 台。

## 限界

- 磁場がないので、アンジュレータの放射・軌道への効果は計算されない (K64 と同じ)。
- track6 (低速追跡) は UND を扱わない。
