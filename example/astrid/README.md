# ASTRID2 の例集 (0.58 GeV)

ASTRID2 のラティスを使った、SAD の基本的な使い方と計算の例です。KEK の SAD チュートリアルの話題を、ASTRID2 のラティスで作り直しました。各例は `.sad`、`.gs.txt`、`.pysad.txt`、`.md` の 4 点です。式と数値の照合は各 `.md` にあります。

## 索引

| フォルダ | 例 | 内容 | gs と pysad |
|---|---|---|---|
| [`ASTRID2_in_notebook.ipynb`](ASTRID2_in_notebook.ipynb) | notebook 版 | ラティスを冒頭で読み込み、セルの計算・放射積分・振幅依存チューンシフト・pysad 独自の ADTS・RDT・運動量依存チューンを順に流す (出力つき) | ○ (`PYSAD_COMPAT=k64` で gs と一致。ADTS・RDT は gs に無い) |
| [`lattice/`](lattice/README.md) | ラティス | ASTRID2 のラティス | |
| [`tune_survey/`](tune_survey/README.md) | [`tune_survey_ncell_ideal_ASTRID2`](tune_survey/tune_survey_ncell_ideal_ASTRID2.md) | ノーマルセルのチューンサーベイ (誤差なし、縮小版) | △ (チューン・εx 一致、DA は ±1〜2 σ。4 点・200 周に縮小) |
| [`tune_survey/`](tune_survey/README.md) | [`tune_survey_ncell_smallerror_ASTRID2`](tune_survey/tune_survey_ncell_smallerror_ASTRID2.md) | ノーマルセルのチューンサーベイ (小さい誤差あり、縮小版) | △ (乱数列が違う。チューン・DA の傾向は一致) |
| [`tune_radiation/`](tune_radiation/README.md) | [`amp_dep_ts_ASTRID2`](tune_radiation/amp_dep_ts_ASTRID2.md) | 振幅依存チューンシフト | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`calc_RadInt_ASTRID2`](tune_radiation/calc_RadInt_ASTRID2.md) | 放射積分 I1〜I5 と、エミッタンスの式による検算 | ○ |
| [`tune_radiation/`](tune_radiation/README.md) | [`mom_dep_ts_ASTRID2`](tune_radiation/mom_dep_ts_ASTRID2.md) | 運動量依存チューンシフト (チューンと運動量ずれの関係) | ○ |
| [`tracking/`](tracking/README.md) | [`phase_x_ASTRID2`](tracking/phase_x_ASTRID2.md) | 位相空間の軌跡 (水平 x-x') | ○ (失われた粒子の座標の符号だけ違う) |
| [`insertion_device/`](insertion_device/README.md) | [`MPW_v01_02_id_ASTRID2`](insertion_device/MPW_v01_02_id_ASTRID2.md) | 挿入光源 MPW の線形シミュレーション | ○ |
| [`lattice_building/`](lattice_building/README.md) | [`SetElement_BeamLine_Join_ASTRID2`](lattice_building/SetElement_BeamLine_Join_ASTRID2.md) | 周期性のあるラティスの作り方 (6 極だけ全セル独立にする) | ○ |

## notebook 版

[`ASTRID2_in_notebook.ipynb`](ASTRID2_in_notebook.ipynb) は、このリングのラティスを冒頭で読み込み、このフォルダの例のコマンドを `%%sad` と `%%sad --continue` で順に流す notebook です (出力つきで保存してあります。`example/tools/make_ring_notebooks.py` が、元の `.sad` からコマンドを取り出して作ります。動かし方は [`../17_notebook/README.md`](../17_notebook/README.md))。

- 保存した出力は pysad の既定 (`v2`) の実行結果です。`PYSAD_COMPAT=k64` で実行すると、εx・σδ・ダンピング分配数・U0 は gs と全桁で一致しました (確認した量: εx 1.26804e-8 m、σδ 4.32609e-4、Jx 1.3751、U0 0.0061667 MV。出力の数値の行のうち、gs の出力に同じ値がある割合は 174 行中 144 行。残りは、書式の違いと 1e-12 以下の数値雑音と、gs にない出力 (pysad 独自の解析式など))。`v2` との差は、`v2` が K64 のバグを直した動作であることによります (docs/28)。
- 確認した範囲: Windows・Python 3.12・IPython 9.19.1・ipykernel 7.4.0・nbclient 0.11.0 で、本物のカーネルを起動して実行しました。ASTRID2 は 30° の複合機能偏向磁石 (K1≠0) を持つので、5 番目のセルの ADTS・RDT (pysad 独自) は、その本体と端の 3 次の項を含みます (追跡との比較は docs/25 §11.8)。

## 参考文献
- ASTRID2 のデザインレポート: *ASTRID2 Preliminary Design Report* (Aarhus University、ISA、2009 年 8 月)。<https://www.isa.au.dk/facilities/astrid2/docs/astrid2-preliminarydesign-report-august2009.pdf>
  ASTRID2 のリングの設計 (エネルギー 0.58 GeV、周長、セル構造、光学) の公開資料。このリポジトリの例のラティスは、設計値 (誤差なし) のラティスで、デザインレポートの値との照合は、一部の主要なパラメータ (エネルギー、周長、セル数) に限り、全部は確かめていない。
