# その他の要素 (COORD・TCAVI・UND・MULT ANGLE・MAP) の検証の詳細 (gs との比較)

本文 [../19_other_elements.md](../19_other_elements.md) から移した、K64 (gs) との細かい比較。元の文書: 25_coord、26_tcav、27_und、28_mult_angle、30_map。見出しの節番号は元の文書のもの。

---
元: 25_coord.md

## 検証 (gs との比較、`pysad/tests/test_coord_elem.py` 12 本)

- 追跡: DX/DY/DZ のみ、CHI のみ、すべて、DIR=1、CHI3 のみ (200 粒子、3 周): 相対 0。
- Twiss と幾何 (INS、BEND を含む並びで 12 要素の BX BY NX R1 R3 EX DX DY と GX GY GZ GCHI1-3): 1e-11 以下。
- Emittance[] (リングの中に COORD を 2 個: 順変換と逆変換、または非対称): エミッタンス・チューン・COD・運動量収縮因子が 5e-10 以下。

---
元: 26_tcav.md

## 検証 (`pysad/tests/test_tcav.py` 8 本)

- 追跡 (gs と比較): RFSW ON、全周 3 周、200 粒子。FREQ、ずれ・回転付き、HARM の小数部 (周回で位相が動く)、L=0、すべて相対 0。RAD FLUC (20 粒子、乱数列まで) も相対 0。
- Emittance[]: L=0 の TCAVI で gs と一致 (エミッタンス 1e-11、チューン・運動量収縮因子・COD 1e-10 台)。ずれ・回転付きは gs の EMT が NaN なので、COD・チューン・運動量収縮因子だけ比べた。
- 追跡と 6 次元の写像が同じ写像であること (二つの直訳の整合): 1 粒子で 1e-19、キックの解析行列は数値微分と 4e-12 で一致。

---
元: 27_und.md

## 検証

- k64 モード (gs との比較、[`pysad/tests/test_und.py`](../../pysad/tests/test_und.py) 4 本): 追跡 (20 粒子、3 周) 相対 0、INS の Twiss と幾何 2e-10 以下、Emittance[] 1e-12 台。
- v2 モード ([`pysad/tests/test_v2_C_und.py`](../../pysad/tests/test_v2_C_und.py) 9 本。gs は磁場を設定できないので比べていない。平面アンジュレータの解析式で確かめた): 軸上の軌道 (出口で x = x' = 0)、行路長 −(F0/kz)² L/4 (1e-6)、垂直の収束 k = F0²/2 (行列要素 5e-4)、KX≠0 の kx・ky (5e-3)、シンプレクティック性 1e-7、放射損失 (2/3) r_e γ³ (F0²/2) L (1e-4)。リング: Δνy = ∫βk ds/4π (1%)、ΔU0 (2e-3)、減衰率の和 = 2U0/E (1e-3)、量子励起は同じ ⟨1/ρ²⟩ の BEND の列との比が ⟨|1/ρ|³⟩ の比 (2%)。
- v2 のスピン追跡 ([`test_v2_C2_spin.py`](../../pysad/tests/test_v2_C2_spin.py)): B = rot A の中の Thomas-BMT の積分と 2e-6 (1 T、1 GeV、POLE = 0.25、5、5.25、中点面の外 2〜3 mm)。
- 例: [`example/v2_C/C1_undulator.sad`](../../example/v2_C/C1_undulator.sad)。

---
元: 28_mult_angle.md

## 検証 (gs との比較、`pysad/tests/test_mult_angle.py` 21 本)

- 追跡 (11 通り、20 粒子、3 周): 相対 0 (K1≠0/=0、K0、スキュー SK2 SK3、E1 E2、フリンジ FB1 FB2 FRINGE=3、
  DX DY ROTATE DROTATE、CHI2、DISFRIN、RAD、RAD+FLUC (乱数も一致)、ANGLE=0.3 の大角度 (ndiv 大))。
- Twiss (CALC): チューン 1e-15、β・α 1e-13。
- Emittance[] (12 セルのリング、MULT が曲げ磁石): ε 1e-11、チューン 1e-13 (K1、K2、K3、K4、SK2、SK3、E1、E2、K0、DX DY ROTATE)。

---
元: 30_map.md

## 検証 (gs との比較、`pysad/tests/test_map.py` 5 本)

- CALC (OPTICS): チューン 1e-14、軌道 1e-11。
- Emittance[] (EMIT、`mm.cod+kick` の整合した写像、x-y 結合の回転 + せん断): ε_x ε_y ε_z、チューン、減衰率、軌道、σδ、σz が 1e-11 以内。4 つ組 (dtrans、dbeam つき) も同じ。
- TrackParticles (TRPT、20 粒子、4 周、MAP を 2 か所): 座標を書き換える写像、フラグ 0 で落とす写像、落ちた粒子の復活を含めて 10 桁一致。ExternalMap が受け取る nt、粒子数、座標の和も一致。
- 既定 (ExternalMap なし): 恒等、S は MAP の L を数えない (gs と同じ 4.1)。

