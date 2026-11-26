# v2.5.0 記録: POL のスピン行列の CHI2・DX・DZ の gs との食い違いの切り分けと BEND の CHI2 (翻訳の抜け)

日付: 2026-11-24。試験 `pysad/tests/script_tests/test_v250_bend_chi2.py` (18 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_bend_chi2__*.{sad,gs.txt}`。

## 1. 食い違いの切り分け

B56 の確認 (`test_v250_spin_tbrote_gs.py`) で、「CHI2・DX・DZ のずれのある偏向磁石の 1 周のスピン行列で、pysad の k64 が gs と 3.8e-3〜1.1e-2 違う」と記録されていた。gs で、3 GeV、24 セルの FODO リング (周長 105.6 m、矩形の偏向磁石、RF・RAD・POL。1 セルは QF D2 B D2 QD DS D1 DS MM、MM は L = 0.2 m の MULT (K1 = 0.05、K2 = 0.3)) の最後の 1 セルの 1 台に、ずれの 1 種類ずつを与えて、1 周のスピン行列 (`EMIT` の「Spin precession matrix」) を k64・v2 と比べた。

| 要素 | ずれ | gs の直交性の欠陥 | k64 − gs (最大) | 判定 |
|---|---|---|---|---|
| BEND | DROTATE = 0.05 | 9.7e-3 | 0 | 一致 (B56) |
| BEND | ROTATE = 0.05 | 5.4e-8 | 0 | 一致 |
| BEND | DX = 1e-3、DY = 1e-3 | 2.9e-8、7.7e-8 | 4e-20、0 | 一致 |
| BEND | DZ = 1e-3 | 1.1e-7 | 5e-21 | 一致 |
| QUAD | ROTATE = 0.05、DX・DY = 1e-3 | 1.3e-7、3.9e-8、1.1e-7 | 2e-16、2e-17、0 | 一致 |
| MULT (ANGLE = 0) | DX・DY・DZ = 1e-3、ROTATE = 0.05 | 6.3e-8、7.0e-8、1.3e-7、1.3e-7 | 1.6e-20、0、8e-31、1e-16 | 一致 (DZ は gs でも効かない) |
| MULT (ANGLE = 0) | CHI1 = 0.03 | 5.7e-8 | **3.8e-4** | pysad が CHI1 を使わない (未実装) |
| MULT (ANGLE = 0) | CHI2 = 0.03 | 1.3e-7 | **2.9e-5** | pysad が CHI2 を使わない (未実装) |
| BEND (ANGLE = 0.26) | CHI2 = ±0.03 | **1.9e-3** (符号によらない) | **2.8e-3** | pysad が BEND の CHI2 を読まない (翻訳の抜け) |
| MULT (ANGLE = 0.26) | CHI2 = 0.03 | 1.9e-3 | (pysad は POL で未対応) | 同上 |

結論: DX・DY・DZ・ROTATE・DROTATE のずれでは k64 は gs と一致する (翻訳の誤りではない)。3.8e-3〜1.1e-2 は、**CHI2 (BEND・MULT) と CHI1 (MULT) を pysad が使っていなかった**ことが原因。

## 2. 原因と修正

- BEND の CHI2: パーサ (`parser._apply_extra_params`) が MULT にだけ CHI2 を読み、BEND には読まなかった。追跡 (`track6n`) は BEND の `g(el, 'CHI2')` を既に引数に持つが 0 だった。CALC・Emittance[] の `tchge` は dchi2 = 0 固定。
- 修正 (k64・v2 とも、翻訳の抜け): パーサが BEND の CHI2 を読む。`lattice._bend_dchi2` (ANGLE ≠ 0 のときだけ。K64 は ANGLE = 0 の BEND の tsteee・tquade・tmulte に CHI2 を渡さない)。`_tchge_cod`・`_tchge_jac`・`spin.tchge_spin` に `dchi2` を足し、`_bend_element_cod_matrix`・`emit6d._tbende`・`_el_has_offset`・`_bend_cod_matrix` が渡す。数値 (numba) の BEND の経路 (`_bend_step_fast`・`elkern.k_bend`・`turnkern` のパック) は CHI2 ≠ 0 で Python の経路に回す。
- ANGLE ≠ 0 の MULT の CHI2: `multae.prep` の NotImplementedError をやめ、`_chg_in` が dchi2 を渡す。

## 3. gs との一致 (修正後、k64)

| 比較 | 結果 |
|---|---|
| 1 要素の追跡 (L = 1 m、ANGLE = 0.3、E1 = 0.5、E2 = 0.3、3 GeV/c、7 粒子) | BEND の CHI2 = 0.03 (DROTATE = 0.02・DX = 1e-3 つきも)、MULT (ANGLE ≠ 0) の CHI2: 1e-15 以内 |
| 1 周のスピン行列 (BEND の CHI2 = ±0.03) | gs と 7 桁で一致 (直交性の欠陥 1.9e-3 も同じ)。v2 は 7.5e-8 (B56 の修正) |
| Emittance[] (POL なし。41 個の数) | BEND の CHI2・MULT (ANGLE ≠ 0) の CHI2: 相対 5e-13 |
| Twiss の最後の ν | 6e-14 (修正前は BEND の CHI2 で 1.3e-4、DY の分散が 0) |

v2 の追跡は gs と 1e-9 以内で一致する (v2 の他の修正の分)。

## 4. 限界・未確認

- CHI2 = 0.03 rad・3 GeV (電子)・1 台だけ。
- **未実装のまま**: MULT (ANGLE = 0) の CHI1・CHI2・DZ (`tsolrot`、追跡・CALC・スピンとも指定しても結果が変わらない。gs では L = 1 m の MULT の CHI2 = 0.03 で追跡が 7.8e-4 動く)。B337 の符号の確認 (gs の MULT の CHI1・CHI2 の符号) を実装の基準にできる。
- ANGLE ≠ 0 の MULT のスピン (POL の Emittance[]) は未対応 (KeyError で止まる)。
- ampnf (TPSA の 1 周の写像) は BEND の CHI2 を入れない。
- ANGLE = 0 の BEND の CHI2 は K64 も使わないと読んだが (`tbende.f`)、gs では確かめていない。
