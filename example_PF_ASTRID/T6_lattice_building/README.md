# T6_lattice_building: 周期性のあるラティスの作り方 (KEK SAD チュートリアル Tutorial6)

`6*NCELL` のように周期的に作ったリングで、六極だけ全セル別の要素にしたいとき、LINE を全部書き下さずに `SetElement` (要素を新しく作る)・`BeamLine` (要素の並びをくくる)・`Join` (つなぐ) で組む例。教材は別の 16 セルのリングで、ここは ASTRID2 (0.58 GeV、周長 45.7044 m、6 セル) の要素で書き直した。ラティスは `../lattice/ASTRID2_v01_01.sad`。
各例は 4 ファイル。走らせ方 (リポジトリの直下から): `py -3.12 example/tools/run_example.py example_PF_ASTRID/T6_lattice_building/NAME.sad`。凡例: ○ 一致 / △ 差あり / × gs 失敗 / × pysad 未対応。

| 例 | チュートリアルの節 | 内容 | gs と pysad |
|---|---|---|---|
| [SetElement_BeamLine_Join_ASTRID2](SetElement_BeamLine_Join_ASTRID2.md) | Tutorial6 全体 | 24 台の六極を別要素 (SDA01…SDB06) にして 6 セルのリングを `Join` で組む。元の `RING` とチューン・εx が一致。セル 3 の SF だけ 0 にすると ξ が式 (±(1/4π)ΣK2ηβ) どおり変わる | ○ |

名前の連結は文字列のリストから (`aa = {"01",...}`、`"SFA"//aa[j]`)。`ToString[j]` や `"0"//id` を使うと名前が切れて全セルが同名になる (md に書いた)。
