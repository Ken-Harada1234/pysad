# lattice: ASTRID2 (0.58 GeV) のラティス

ASTRID2 (オーフス大学の放射光源リング) の台本 2 本。元の台本のコピーを改変せずに置いた (中身は同じ。拡張子 `.sad` を付け、名前に `ASTRID2_` を付けた)。PF のラティスは同じフォルダの `v9_16_*.sad` と `README_PF.md` (もう一人の担当)。

| ファイル | 内容 | 使う例 |
|---|---|---|
| `ASTRID2_v01_01.sad` | ASTRID2 の周回ラティス (誤差なし、設計値)。ライン `NCELL` (標準セル 1 個、7.6174 m)、`RING = (3*NCELL CAVI 3*NCELL)` (200 要素、周長 45.7044 m)。MARK `MKS`・`MKE` (セルの両端)、モニター `PM`、四極 `QF3`・`QF2`、六極 `SF`・`SD`、複合機能の偏向磁石 `BH`、空洞 `CAVI`、回転 0/90 度のステアリング `PH`・`PV` (`HV = (PH PV)`) | T2: tune_survey (ideal・smallerror)、mom_dep_ts、amp_dep_ts、calc_RadInt、phase_x。T6: SetElement_BeamLine_Join |
| `ASTRID2_v01_02_id.sad` | 上に MPW (挿入光源 11 周期、K=32.7、λ=0.175 m、B0=2.00 T) を加えたもの。MPW のライン `MPW`、MPW 入りセル `TCELL`、MPW 入りリング `RINGMPW = (2*NCELL CAVI TCELL 2*NCELL)` (223 要素)。MPW の両側の直線部 `D14837MPW` (0.5212 m) を短くして周長 45.7044 m を保つ | T2: MPW_v01_02_id |

## 出所
手元のラティスのアーカイブ (ASTRID2 のフォルダ) の `v01_01` と `v01_02_id`。チュートリアル (KEK SAD の Tutorial2) の `ASTRID2_lattice_12nmrad.sad` (チューンサーベイ用) と `v01_02_id.sad` (MPW 用) に当たる。`ASTRID2_lattice_12nmrad` という名前のファイルはアーカイブに無く、`v01_01` の εx が 12.68 nm なので、これを `ASTRID2_lattice_12nmrad` の代わりに使った (中身が同じかは確かめていない)。
アーカイブの ASTRID2 フォルダの `v9_16_asym` の中身は ASTRID2 でなく **PF のラティス** (見出しが "PF-RING Super-Low Emittance Lattice")。`v01_03id`・`v01_04_inj`・`v02_01` (見出しに別施設の名前があるものを含む) は使っていない。

## 確認した値 (gs と pysad で全桁一致、`T2_samples/MPW_v01_02_id_ASTRID2`・`calc_RadInt_ASTRID2`)
- エネルギー 0.58 GeV (`MOMENTUM= 0.58 GEV`)、周長 45.7044 m、偏向磁石 12 台 (BH: L=0.85 m、ANGLE=0.523599 rad、K1 の積分値 -1.6782×0.85 = -1.4265 (密度 -1.6782 /m²)、E1=E2=0 の扇形)、セル 6 個、直線部 6 か所 (各 2.9674 m、セルの両端の `D14837` = 1.4837 m の和)。
- チューン リング (5.17999, 2.14018)、セル (0.86333, 0.35670)。MKS での β = (3.8916, 1.6225) m、分散はゼロ (-8.3e-4 m)。
- εx = 12.680 nm、エネルギー拡がり 4.326e-4、運動量コンパクション率 0.010669、U0 = 6166.7 eV、減衰分配数 (Jx, Jy, Jz) = (1.375, 1.000, 1.625)、減衰時間 (x, y, z) = (20.9, 28.7, 17.6) ms。
- 空洞 `CAVI`: VOLT = 50 kV、FREQ = 104.95 MHz (ハーモニック数 16)、シンクロトロン周波数 10.0 kHz。
- 六極は台本の最後で上書きされている: `SF` K2=8.6814 (L=0.075 m)、`SD` K2=-11.7295 (L=0.15 m)。最初の定義 (`K2 = 125.62*0.075`、`-81.93*0.15`) は使われない (後の定義が勝つ)。このときのクロマティシティは (-2.1e-4, -2.8e-4) で、色収差補正された値。

## 使い方
- `FFS;` の前にラティスを置き、後ろに計算を続ける (各例の `.sad` は、このラティスを先頭に貼った自己完結のファイル)。
- チューンサーベイの例は、`USE NCELL` と `USE RING` を行き来する。**`USE` は保存値 (`SAVE ALL`) に戻す**ので、FIT などで磁石を変えたら `SAVE ALL` してから別のラインに `USE` する (`tune_survey_ncell_ideal_ASTRID2.md`)。

## 未確認・限界
- 元の台本は設計値で、誤差は入っていない (誤差入りは `tune_survey_ncell_smallerror_ASTRID2` が乱数で入れる)。
- `v01_01` と `ASTRID2_lattice_12nmrad` が同じ台本かは確かめていない。MPW のモデルの θ・長さはチュートリアルの値 (ASTRID2 の設計報告書そのものとの照合はしていない)。
